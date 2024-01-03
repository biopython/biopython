# Copyright 2022 by Michiel de Hoon.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.Align support for alignment files in the bigBed format.

The bigBed format stores a series of pairwise alignments in a single indexed
binary file. Typically they are used for transcript to genome alignments. As
in the BED format, the alignment positions and alignment scores are stored,
but the aligned sequences are not.

See http://genome.ucsc.edu/goldenPath/help/bigBed.html for more information.

You are expected to use this module via the Bio.Align functions.
"""

# This parser was written based on the description of the bigBed file format in
# W. J. Kent, A. S. Zweig,* G. Barber, A. S. Hinrichs, and D. Karolchik:
# "BigWig and BigBed: enabling browsing of large distributed datasets."
# Bioinformatics 26(17): 2204–2207 (2010)
# in particular the tables in the supplemental materials listing the contents
# of a bigBed file byte-by-byte.

# The writer is based on files in the src/lib and src/utils/bedToBigBed
# directories in Jim Kent's source code from UCSC.

# These program files are provided by the University of California Santa Cruz
# under the following license:

# ------------------------------------------------------------------------------
# Copyright (C) 2001 UC Regents
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
# ------------------------------------------------------------------------------


import sys
import io
import copy
import array
import itertools
import struct
import zlib
from collections import namedtuple
from io import BytesIO
import numpy as np


from Bio.Align import Alignment
from Bio.Align import interfaces
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


Field = namedtuple("Field", ("as_type", "name", "comment"))


class AutoSQLTable(list):
    """AutoSQL table describing the columns of an (possibly extended) BED format."""

    default: "AutoSQLTable"

    def __init__(self, name, comment, fields):
        """Create an AutoSQL table describing the columns of an (extended) BED format."""
        self.name = name
        self.comment = comment
        self[:] = fields

    @classmethod
    def from_bytes(cls, data):
        """Return an AutoSQLTable initialized using the bytes object data."""
        assert data.endswith(b"\0")  # NULL-terminated string
        text = data[:-1].decode()
        word, text = text.split(None, 1)
        assert word == "table"
        name, text = text.split(None, 1)
        assert text.startswith('"')
        i = text.find('"', 1)
        comment = text[1:i]
        text = text[i + 1 :].strip()
        assert text.startswith("(")
        assert text.endswith(")")
        text = text[1:-1].strip()
        fields = []
        while text:
            i = text.index('"')
            j = text.index('"', i + 1)
            field_comment = text[i + 1 : j]
            definition = text[:i].strip()
            assert definition.endswith(";")
            field_type, field_name = definition[:-1].rsplit(None, 1)
            if field_type.endswith("]"):
                i = field_type.index("[")
                data_type = field_type[:i]
            else:
                data_type = field_type
            assert data_type in (
                "int",
                "uint",
                "short",
                "ushort",
                "byte",
                "ubyte",
                "float",
                "char",
                "string",
                "lstring",
            )
            field = Field(field_type, field_name, field_comment)
            fields.append(field)
            text = text[j + 1 :].strip()
        return AutoSQLTable(name, comment, fields)

    @classmethod
    def from_string(cls, data):
        """Return an AutoSQLTable initialized using the string object data."""
        return cls.from_bytes(data.encode() + b"\0")

    def __str__(self):
        type_width = max(len(str(field.as_type)) for field in self)
        name_width = max(len(field.name) for field in self) + 1
        lines = []
        lines.append("table %s\n" % self.name)
        lines.append('"%s"\n' % self.comment)
        lines.append("(\n")
        for field in self:
            name = field.name + ";"
            lines.append(
                f'   %-{type_width}s %-{name_width}s    "%s"\n'
                % (
                    field.as_type,
                    name,
                    field.comment,
                )
            )
        lines.append(")\n")
        return "".join(lines)

    def __bytes__(self):
        return str(self).encode() + b"\0"

    def __getitem__(self, i):
        if isinstance(i, slice):
            fields = super().__getitem__(i)
            return AutoSQLTable(self.name, self.comment, fields)
        else:
            return super().__getitem__(i)


AutoSQLTable.default = AutoSQLTable(
    "bed",
    "Browser Extensible Data",
    [
        Field(
            as_type="string",
            name="chrom",
            comment="Reference sequence chromosome or scaffold",
        ),
        Field(
            as_type="uint",
            name="chromStart",
            comment="Start position in chromosome",
        ),
        Field(
            as_type="uint",
            name="chromEnd",
            comment="End position in chromosome",
        ),
        Field(as_type="string", name="name", comment="Name of item."),
        Field(as_type="uint", name="score", comment="Score (0-1000)"),
        Field(as_type="char[1]", name="strand", comment="+ or - for strand"),
        Field(
            as_type="uint",
            name="thickStart",
            comment="Start of where display should be thick (start codon)",
        ),
        Field(
            as_type="uint",
            name="thickEnd",
            comment="End of where display should be thick (stop codon)",
        ),
        Field(
            as_type="uint",
            name="reserved",
            comment="Used as itemRgb as of 2004-11-22",
        ),
        Field(as_type="int", name="blockCount", comment="Number of blocks"),
        Field(
            as_type="int[blockCount]",
            name="blockSizes",
            comment="Comma separated list of block sizes",
        ),
        Field(
            as_type="int[blockCount]",
            name="chromStarts",
            comment="Start positions relative to chromStart",
        ),
    ],
)


class AlignmentWriter(interfaces.AlignmentWriter):
    """Alignment file writer for the bigBed file format."""

    fmt = "bigBed"
    mode = "wb"

    def __init__(
        self,
        target,
        bedN=12,
        declaration=None,
        targets=None,
        compress=True,
        extraIndex=(),
    ):
        """Create an AlignmentWriter object.

        Arguments:
         - target      - output stream or file name.
         - bedN        - number of columns in the BED file.
                         This must be between 3 and 12; default value is 12.
         - declaration - an AutoSQLTable object declaring the fields in the BED file.
                         Required only if the BED file contains extra (custom) fields.
                         Default value is None.
         - targets     - A list of SeqRecord objects with the chromosomes in the
                         order as they appear in the alignments. The sequence
                         contents in each SeqRecord may be undefined, but the
                         sequence length must be defined, as in this example:

                         SeqRecord(Seq(None, length=248956422), id="chr1")

                         If targets is None (the default value), the alignments
                         must have an attribute .targets providing the list of
                         SeqRecord objects.
         - compress    - If True (default), compress data using zlib.
                         If False, do not compress data.
         - extraIndex  - List of strings with the names of extra columns to be
                         indexed.
                         Default value is an empty list.
        """
        if bedN < 3 or bedN > 12:
            raise ValueError("bedN must be between 3 and 12")
        super().__init__(target)
        self.bedN = bedN
        self.declaration = declaration
        self.targets = targets
        self.compress = compress
        self.extraIndexNames = extraIndex
        self.itemsPerSlot = 512
        self.blockSize = 256

    def write_file(self, stream, alignments):
        """Write the alignments to the file strenm, and return the number of alignments.

        alignments - A list or iterator returning Alignment objects
        stream     - Output file stream.
        """
        if self.targets is None:
            targets = alignments.targets
        else:
            targets = self.targets
        header = _Header()
        header.definedFieldCount = self.bedN
        if self.declaration is None:
            try:
                self.declaration = alignments.declaration[: self.bedN]
            except AttributeError:
                self.declaration = AutoSQLTable.default[: self.bedN]
        declaration = self.declaration
        header.fieldCount = len(declaration)
        extra_indices = _ExtraIndices(self.extraIndexNames, declaration)
        chromUsageList, aveSize = self._get_chrom_usage(
            alignments, targets, extra_indices
        )
        stream.write(bytes(header.size))
        stream.write(bytes(_ZoomLevels.size))
        header.autoSqlOffset = stream.tell()
        stream.write(bytes(declaration))
        header.totalSummaryOffset = stream.tell()
        stream.write(bytes(_Summary.size))
        header.extraIndicesOffset = stream.tell()
        stream.write(bytes(extra_indices.size))
        header.chromosomeTreeOffset = stream.tell()
        _BPlusTreeFormatter().write(
            chromUsageList, min(self.blockSize, len(chromUsageList)), stream
        )
        header.fullDataOffset = stream.tell()
        reductions = _ZoomLevels.calculate_reductions(aveSize)
        stream.write(len(alignments).to_bytes(8, sys.byteorder))
        extra_indices.initialize(len(alignments))
        maxBlockSize, regions = self.write_alignments(
            alignments,
            stream,
            reductions,
            extra_indices,
        )
        if self.compress:
            header.uncompressBufSize = max(
                maxBlockSize, self.itemsPerSlot * _RegionSummary.size
            )
        else:
            header.uncompressBufSize = 0
        header.fullIndexOffset = stream.tell()
        _RTreeFormatter().write(
            regions, self.blockSize, 1, header.fullIndexOffset, stream
        )
        zoomList, totalSum = self._write_zoom_levels(
            alignments,
            stream,
            header.fullIndexOffset - header.fullDataOffset,
            chromUsageList,
            reductions,
        )
        header.zoomLevels = len(zoomList)
        for extra_index in extra_indices:
            extra_index.fileOffset = stream.tell()
            extra_index.chunks.sort()
            _BPlusTreeFormatter().write(extra_index.chunks, self.blockSize, stream)
        stream.seek(0)
        stream.write(bytes(header))
        stream.write(bytes(zoomList))
        stream.seek(header.totalSummaryOffset)
        stream.write(bytes(totalSum))
        assert header.extraIndicesOffset == stream.tell()
        extra_indices.tofile(stream)
        stream.seek(0, io.SEEK_END)
        data = header.signature.to_bytes(4, sys.byteorder)
        stream.write(data)

    def _get_chrom_usage(cls, alignments, targets, extra_indices):
        aveSize = 0
        chromId = 0
        totalBases = 0
        bedCount = 0
        name = ""
        chromUsageList = []
        keySize = 0
        chromSize = -1
        minDiff = sys.maxsize
        for alignment in alignments:
            chrom = alignment.sequences[0].id
            # Note that for MAF files, this may be a multiple alignment, so
            # alignment.target.id won't necessarily work.
            start = alignment.coordinates[0, 0]
            end = alignment.coordinates[0, -1]
            for extra_index in extra_indices:
                extra_index.updateMaxFieldSize(alignment)
            if start > end:
                raise ValueError(
                    f"end ({end}) before start ({start}) in alignment [{bedCount}]"
                )
            bedCount += 1
            totalBases += end - start
            if name != chrom:
                if name > chrom:
                    raise ValueError(
                        f"alignments are not sorted by target name at alignment [{bedCount}]"
                    )
                if name:
                    chromUsageList.append((name, chromId, chromSize))
                    chromId += 1
                for target in targets:
                    if target.id == chrom:
                        break
                else:
                    raise ValueError(
                        f"failed to find target '{chrom}' in target list at alignment [{bedCount}]"
                    )
                name = chrom
                keySize = max(keySize, len(chrom))
                chromSize = len(target)
                lastStart = -1
            if end > chromSize:
                raise ValueError(
                    f"end coordinate {end} bigger than {chrom} size of {chromSize} at alignment [{bedCount}]'"
                )
            if lastStart >= 0:
                diff = start - lastStart
                if diff < minDiff:
                    if diff < 0:
                        raise ValueError(
                            f"alignments are not sorted at alignment [{bedCount}]"
                        )
                    minDiff = diff
            lastStart = start
        if name:
            chromUsageList.append((name, chromId, chromSize))
        chromUsageList = np.array(
            chromUsageList,
            dtype=[("name", f"S{keySize}"), ("id", "=i4"), ("size", "=i4")],
        )
        if bedCount > 0:
            aveSize = totalBases / bedCount
        alignments._len = bedCount
        return chromUsageList, aveSize

    def _write_zoom_levels(
        self, alignments, output, dataSize, chromUsageList, reductions
    ):
        zoomList = _ZoomLevels()
        totalSum = _Summary()
        if len(alignments) == 0:
            totalSum.minVal = 0.0
            totalSum.maxVal = 0.0
        else:
            blockSize = self.blockSize
            doCompress = self.compress
            itemsPerSlot = self.itemsPerSlot
            maxReducedSize = dataSize / 2
            zoomList[0].dataOffset = output.tell()
            for initialReduction in reductions:
                reducedSize = initialReduction["size"] * _RegionSummary.size
                if doCompress:
                    reducedSize /= 2
                if reducedSize <= maxReducedSize:
                    break
            else:
                initialReduction = reductions[0]
            initialReduction["size"].tofile(output)
            size = itemsPerSlot * _RegionSummary.size
            if doCompress:
                buffer = _ZippedBufferedStream(output, size)
            else:
                buffer = _BufferedStream(output, size)
            regions = []
            rezoomedList = []
            trees = _RangeTree.generate(chromUsageList, alignments)
            scale = initialReduction["scale"]
            doubleReductionSize = scale * _ZoomLevels.bbiResIncrement
            for tree in trees:
                start = -sys.maxsize
                summaries = tree.generate_summaries(scale, totalSum)
                for summary in summaries:
                    buffer.write(summary)
                    regions.append(summary)
                    if start + doubleReductionSize < summary.end:
                        rezoomed = copy.copy(summary)
                        start = rezoomed.start
                        rezoomedList.append(rezoomed)
                    else:
                        rezoomed += summary
            buffer.flush()
            assert len(regions) == initialReduction["size"]
            zoomList[0].amount = initialReduction["scale"]
            indexOffset = output.tell()
            zoomList[0].indexOffset = indexOffset
            _RTreeFormatter().write(
                regions, blockSize, itemsPerSlot, indexOffset, output
            )
            if doCompress:
                buffer = _ZippedBufferedStream(output, size)
            else:
                buffer = _BufferedStream(output, _RegionSummary.size)
            zoomList.reduce(
                rezoomedList, initialReduction, buffer, blockSize, itemsPerSlot
            )
        return zoomList, totalSum

    def _extract_fields(self, alignment):
        bedN = self.bedN
        row = []
        chrom = alignment.target.id
        if len(chrom) >= 255:
            raise ValueError(
                f"alignment target name '{chrom}' is too long (must not exceed 254 characters)"
            )
        if len(chrom) < 1:
            raise ValueError("alignment target name cannot be blank or empty")
        chromStart = alignment.coordinates[0, 0]
        chromEnd = alignment.coordinates[0, -1]
        if chromEnd < chromStart:
            raise ValueError(f"chromStart after chromEnd ({chromEnd} > {chromStart})")
        if bedN > 3:
            name = alignment.query.id
            if name == "":
                name = "."
            elif len(name) > 255:
                raise ValueError(
                    f"alignment query name '{name}' is too long (must not exceed 255 characters"
                )
            row.append(name)
        if bedN > 4:
            try:
                score = alignment.score
            except AttributeError:
                score = "."
            else:
                if score < 0 or score > 1000:
                    raise ValueError(f"score ({score}) must be between 0 and 1000")
                score = str(score)
            row.append(score)
        if bedN > 5:
            if alignment.coordinates[1, 0] <= alignment.coordinates[1, -1]:
                strand = "+"
            else:
                strand = "-"
            row.append(strand)
        if bedN > 6:
            try:
                thickStart = alignment.thickStart
            except AttributeError:
                thickStart = chromStart
            row.append(str(thickStart))
        if bedN > 7:
            try:
                thickEnd = alignment.thickEnd
            except AttributeError:
                thickEnd = chromEnd
            else:
                if thickEnd < thickStart:
                    raise ValueError(
                        f"thickStart ({thickStart}) after thickEnd ({thickEnd})"
                    )
                if thickStart != 0 and (
                    thickStart < chromStart or thickStart > chromEnd
                ):
                    raise ValueError(
                        f"thickStart out of range for {name}:{chromStart}-{chromEnd}, thick:{thickStart}-{thickEnd}"
                    )
                if thickEnd != 0 and (thickEnd < chromStart or thickEnd > chromEnd):
                    raise ValueError(
                        f"thickEnd out of range for {name}:{chromStart}-{chromEnd}, thick:{thickStart}-{thickEnd}"
                    )
            row.append(str(thickEnd))
        if bedN > 8:
            try:
                itemRgb = alignment.itemRgb
            except AttributeError:
                itemRgb = "."
            else:
                colors = itemRgb.rstrip(",").split(",")
                if len(colors) == 3 and all(0 <= int(color) < 256 for color in colors):
                    pass
                elif 0 <= int(itemRgb) < (2 << 32):
                    pass
                else:
                    raise ValueError(
                        f"Expecting color to consist of r,g,b values from 0 to 255. Got '{itemRgb}'"
                    )
            row.append(itemRgb)
        if bedN > 9:
            steps = np.diff(alignment.coordinates)
            aligned = sum(steps != 0, 0) == 2
            blockSizes = steps.max(0)[aligned]
            blockCount = len(blockSizes)
            row.append(str(blockCount))
        if bedN > 10:
            row.append(",".join(str(blockSize) for blockSize in blockSizes) + ",")
        if bedN > 11:
            chromStarts = alignment.coordinates[0, :-1][aligned] - chromStart
            row.append(",".join(str(chromStart) for chromStart in chromStarts) + ",")
        if bedN > 12:
            if bedN != 15:
                raise ValueError(f"Unexpected value {bedN} for bedN in _extract_fields")
            expIds = alignment.annotations["expIds"]
            expScores = alignment.annotations["expScores"]
            expCount = len(expIds)
            assert expCount == len(expScores)
            row.append(str(expCount))
            row.append(",".join(expIds))
            row.append(",".join(str(expScore) for expScore in expScores))
        for field in self.declaration[bedN:]:
            value = alignment.annotations[field.name]
            if isinstance(value, str):
                row.append(value)
            elif isinstance(value, (int, float)):
                row.append(str(value))
            else:
                row.append(",".join(map(str, value)))
        rest = "\t".join(row).encode()
        return chrom, chromStart, chromEnd, rest

    def write_alignments(self, alignments, output, reductions, extra_indices):
        """Write alignments to the output file, and return the number of alignments.

        alignments - A list or iterator returning Alignment objects
        stream     - Output file stream.
        """
        itemsPerSlot = self.itemsPerSlot
        chromId = -1
        itemIx = 0
        sectionStartIx = 0
        sectionEndIx = 0
        currentChrom = None

        regions = []
        if self.compress is True:
            buffer = _ZippedStream()
        else:
            buffer = BytesIO()
        maxBlockSize = 0

        # Supplemental Table 12: Binary BED-data format
        # chromId     4 bytes, unsigned
        # chromStart  4 bytes, unsigned
        # chromEnd    4 bytes, unsigned
        # rest        zero-terminated string in tab-separated format
        formatter = struct.Struct("=III")

        done = False
        region = None
        alignments.rewind()
        while True:
            try:
                alignment = next(alignments)
            except StopIteration:
                itemIx = itemsPerSlot
                done = True
            else:
                chrom, start, end, rest = self._extract_fields(alignment)
                if chrom != currentChrom:
                    if currentChrom is not None:
                        itemIx = itemsPerSlot
                    currentChrom = chrom
                    chromId += 1
                    reductions["end"] = 0
            if itemIx == itemsPerSlot:
                blockStartOffset = output.tell()
                size = buffer.tell()
                if size > maxBlockSize:
                    maxBlockSize = size
                data = buffer.getvalue()
                output.write(data)
                buffer.seek(0)
                buffer.truncate(0)
                if extra_indices:
                    blockEndOffset = output.tell()
                    blockSize = blockEndOffset - blockStartOffset
                    for extra_index in extra_indices:
                        extra_index.addOffsetSize(
                            blockStartOffset,
                            blockSize,
                            sectionStartIx,
                            sectionEndIx,
                        )
                    sectionStartIx = sectionEndIx
                region.offset = blockStartOffset
                if done is True:
                    break
                itemIx = 0
            if itemIx == 0:
                region = _Region(chromId, start, end)
                regions.append(region)
            elif end > region.end:
                region.end = end
            itemIx += 1
            for row in reductions:
                if start >= row["end"]:
                    row["size"] += 1
                    row["end"] = start + row["scale"]
                while end > row["end"]:
                    row["size"] += 1
                    row["end"] += row["scale"]
            if extra_indices:
                for extra_index in extra_indices:
                    extra_index.addKeysFromRow(alignment, sectionEndIx)
                sectionEndIx += 1
            data = formatter.pack(chromId, start, end)
            buffer.write(data + rest + b"\0")
        return maxBlockSize, regions


class AlignmentIterator(interfaces.AlignmentIterator):
    """Alignment iterator for bigBed files.

    The pairwise alignments stored in the bigBed file are loaded and returned
    incrementally.  Additional alignment information is stored as attributes
    of each alignment.
    """

    fmt = "bigBed"
    mode = "b"

    def _read_header(self, stream):
        header = _Header.fromfile(stream)
        byteorder = header.byteorder
        autoSqlOffset = header.autoSqlOffset
        self.byteorder = byteorder
        fieldCount = header.fieldCount
        definedFieldCount = header.definedFieldCount
        fullDataOffset = header.fullDataOffset
        self.declaration = self._read_autosql(stream, header)
        stream.seek(fullDataOffset)
        (dataCount,) = struct.unpack(byteorder + "Q", stream.read(8))
        self._length = dataCount

        if header.uncompressBufSize > 0:
            self._compressed = True
        else:
            self._compressed = False

        stream.seek(header.chromosomeTreeOffset)
        self.targets = _BPlusTreeFormatter(byteorder).read(stream)

        stream.seek(header.fullIndexOffset)
        self.tree = _RTreeFormatter(byteorder).read(stream)

        self._data = self._iterate_index(stream)

    def _read_autosql(self, stream, header):
        autoSqlSize = header.totalSummaryOffset - header.autoSqlOffset
        fieldCount = header.fieldCount
        self.bedN = header.definedFieldCount
        stream.seek(header.autoSqlOffset)
        data = stream.read(autoSqlSize)
        declaration = AutoSQLTable.from_bytes(data)
        self._analyze_fields(declaration, fieldCount, self.bedN)
        return declaration

    def _analyze_fields(self, fields, fieldCount, definedFieldCount):
        names = (
            "chrom",
            "chromStart",
            "chromEnd",
            "name",
            "score",
            "strand",
            "thickStart",
            "thickEnd",
            "reserved",
            "blockCount",
            "blockSizes",
            "chromStarts",
        )
        for i in range(self.bedN):
            name = fields[i].name
            if name != names[i]:
                raise ValueError(
                    "Expected field name '%s'; found '%s'" % (names[i], name)
                )
        if fieldCount > definedFieldCount:
            self._custom_fields = []
        for i in range(definedFieldCount, fieldCount):
            field_name = fields[i].name
            field_type = fields[i].as_type
            if "[" in field_type and "]" in field_type:
                make_array = True
                field_type, _ = field_type.split("[")
                field_type = field_type.strip()
            else:
                make_array = False
            if field_type in ("int", "uint", "short", "ushort"):
                converter = int
            elif field_type in ("byte", "ubyte"):
                converter = bytes
            elif field_type == "float":
                converter = float
            elif field_type in ("float", "char", "string", "lstring"):
                converter = str
            else:
                raise Exception("Unknown field type %s" % field_type)
            if make_array:
                item_converter = converter

                def converter(data, item_converter=item_converter):
                    values = data.rstrip(",").split(",")
                    return [item_converter(value) for value in values]

            self._custom_fields.append([field_name, converter])

    def _iterate_index(self, stream):
        # Supplemental Table 12: Binary BED-data format
        # chromId     4 bytes, unsigned
        # chromStart  4 bytes, unsigned
        # chromEnd    4 bytes, unsigned
        # rest        zero-terminated string in tab-separated format
        formatter = struct.Struct(self.byteorder + "III")
        size = formatter.size
        node = self.tree
        while True:
            try:
                children = node.children
            except AttributeError:
                stream.seek(node.dataOffset)
                data = stream.read(node.dataSize)
                if self._compressed > 0:
                    data = zlib.decompress(data)
                while data:
                    chromId, chromStart, chromEnd = formatter.unpack(data[:size])
                    rest, data = data[size:].split(b"\00", 1)
                    yield (chromId, chromStart, chromEnd, rest)
                while True:
                    parent = node.parent
                    if parent is None:
                        return
                    for index, child in enumerate(parent.children):
                        if id(node) == id(child):
                            break
                    else:
                        raise RuntimeError("Failed to find child node")
                    try:
                        node = parent.children[index + 1]
                    except IndexError:
                        node = parent
                    else:
                        break
            else:
                node = children[0]

    def _search_index(self, stream, chromIx, start, end):
        # Supplemental Table 12: Binary BED-data format
        # chromId     4 bytes, unsigned
        # chromStart  4 bytes, unsigned
        # chromEnd    4 bytes, unsigned
        # rest        zero-terminated string in tab-separated format
        formatter = struct.Struct(self.byteorder + "III")
        size = formatter.size
        padded_start = start - 1
        padded_end = end + 1
        node = self.tree
        while True:
            try:
                children = node.children
            except AttributeError:
                stream.seek(node.dataOffset)
                data = stream.read(node.dataSize)
                if self._compressed > 0:
                    data = zlib.decompress(data)
                while data:
                    child_chromIx, child_chromStart, child_chromEnd = formatter.unpack(
                        data[:size]
                    )
                    rest, data = data[size:].split(b"\00", 1)
                    if child_chromIx != chromIx:
                        continue
                    if end <= child_chromStart or child_chromEnd <= start:
                        if child_chromStart != child_chromEnd:
                            continue
                        if child_chromStart != end and child_chromEnd != start:
                            continue
                    yield (child_chromIx, child_chromStart, child_chromEnd, rest)
            else:
                visit_child = False
                for child in children:
                    if (child.endChromIx, child.endBase) < (chromIx, padded_start):
                        continue
                    if (chromIx, padded_end) < (child.startChromIx, child.startBase):
                        continue
                    visit_child = True
                    break
                if visit_child:
                    node = child
                    continue
            while True:
                parent = node.parent
                if parent is None:
                    return
                for index, child in enumerate(parent.children):
                    if id(node) == id(child):
                        break
                else:
                    raise RuntimeError("Failed to find child node")
                try:
                    node = parent.children[index + 1]
                except IndexError:
                    node = parent
                else:
                    break

    def _read_next_alignment(self, stream):
        try:
            row = next(self._data)
        except StopIteration:
            return
        return self._create_alignment(row)

    def _create_alignment(self, row):
        chromId, chromStart, chromEnd, rest = row
        if rest:
            words = rest.decode().split("\t")
        else:
            words = []
        target_record = self.targets[chromId]
        if self.bedN > 3:
            name = words[0]
        else:
            name = None
        if self.bedN > 5:
            strand = words[2]
        else:
            strand = "+"
        if self.bedN > 9:
            blockCount = int(words[6])
            blockSizes = np.fromiter(words[7].rstrip(",").split(","), int)
            blockStarts = np.fromiter(words[8].rstrip(",").split(","), int)
            if len(blockSizes) != blockCount:
                raise ValueError(
                    "Inconsistent number of block sizes (%d found, expected %d)"
                    % (len(blockSizes), blockCount)
                )
            if len(blockStarts) != blockCount:
                raise ValueError(
                    "Inconsistent number of block start positions (%d found, expected %d)"
                    % (len(blockStarts), blockCount)
                )
            tPosition = 0
            qPosition = 0
            coordinates = [[tPosition, qPosition]]
            for blockSize, blockStart in zip(blockSizes, blockStarts):
                if blockStart != tPosition:
                    coordinates.append([blockStart, qPosition])
                    tPosition = blockStart
                tPosition += blockSize
                qPosition += blockSize
                coordinates.append([tPosition, qPosition])
            coordinates = np.array(coordinates).transpose()
            qSize = sum(blockSizes)
        else:
            blockSize = chromEnd - chromStart
            coordinates = np.array([[0, blockSize], [0, blockSize]])
            qSize = blockSize
        coordinates[0, :] += chromStart
        query_sequence = Seq(None, length=qSize)
        query_record = SeqRecord(query_sequence, id=name)
        records = [target_record, query_record]
        if strand == "-":
            coordinates[1, :] = qSize - coordinates[1, :]
        if chromStart != coordinates[0, 0]:
            raise ValueError(
                "Inconsistent chromStart found (%d, expected %d)"
                % (chromStart, coordinates[0, 0])
            )
        if chromEnd != coordinates[0, -1]:
            raise ValueError(
                "Inconsistent chromEnd found (%d, expected %d)"
                % (chromEnd, coordinates[0, -1])
            )
        alignment = Alignment(records, coordinates)
        if len(words) > self.bedN - 3:
            alignment.annotations = {}
            for word, custom_field in zip(words[self.bedN - 3 :], self._custom_fields):
                name, converter = custom_field
                alignment.annotations[name] = converter(word)
        if self.bedN <= 4:
            return alignment
        score = words[1]
        try:
            score = float(score)
        except ValueError:
            pass
        else:
            if score.is_integer():
                score = int(score)
        alignment.score = score
        if self.bedN <= 6:
            return alignment
        alignment.thickStart = int(words[3])
        if self.bedN <= 7:
            return alignment
        alignment.thickEnd = int(words[4])
        if self.bedN <= 8:
            return alignment
        alignment.itemRgb = words[5]
        return alignment

    def __len__(self):
        return self._length

    def search(self, chromosome=None, start=None, end=None):
        """Iterate over alignments overlapping the specified chromosome region..

        This method searches the index to find alignments to the specified
        chromosome that fully or partially overlap the chromosome region
        between start and end.

        Arguments:
         - chromosome - chromosome name. If None (default value), include all
           alignments.
         - start      - starting position on the chromosome. If None (default
           value), use 0 as the starting position.
         - end        - end position on the chromosome. If None (default value),
           use the length of the chromosome as the end position.

        """
        stream = self._stream
        if chromosome is None:
            if start is not None or end is not None:
                raise ValueError(
                    "start and end must both be None if chromosome is None"
                )
        else:
            for chromIx, target in enumerate(self.targets):
                if target.id == chromosome:
                    break
            else:
                raise ValueError("Failed to find %s in alignments" % chromosome)
            if start is None:
                if end is None:
                    start = 0
                    end = len(target)
                else:
                    raise ValueError("end must be None if start is None")
            elif end is None:
                end = start + 1
        data = self._search_index(stream, chromIx, start, end)
        for row in data:
            alignment = self._create_alignment(row)
            yield alignment


class _ZippedStream(io.BytesIO):
    def getvalue(self):
        data = super().getvalue()
        return zlib.compress(data)


class _BufferedStream:
    def __init__(self, output, size):
        self.buffer = BytesIO()
        self.output = output
        self.size = size

    def write(self, item):
        item.offset = self.output.tell()
        data = bytes(item)
        self.buffer.write(data)
        if self.buffer.tell() == self.size:
            self.output.write(self.buffer.getvalue())
            self.buffer.seek(0)
            self.buffer.truncate(0)

    def flush(self):
        self.output.write(self.buffer.getvalue())
        self.buffer.seek(0)
        self.buffer.truncate(0)


class _ZippedBufferedStream(_BufferedStream):
    def write(self, item):
        item.offset = self.output.tell()
        data = bytes(item)
        self.buffer.write(data)
        if self.buffer.tell() == self.size:
            self.output.write(zlib.compress(self.buffer.getvalue()))
            self.buffer.seek(0)
            self.buffer.truncate(0)

    def flush(self):
        self.output.write(zlib.compress(self.buffer.getvalue()))
        self.buffer.seek(0)
        self.buffer.truncate(0)


class _Header:
    __slots__ = (
        "byteorder",
        "zoomLevels",
        "chromosomeTreeOffset",
        "fullDataOffset",
        "fullIndexOffset",
        "fieldCount",
        "definedFieldCount",
        "autoSqlOffset",
        "totalSummaryOffset",
        "uncompressBufSize",
        "extraIndicesOffset",
    )

    # Supplemental Table 5: Common header
    # magic                  4 bytes, unsigned
    # version                2 bytes, unsigned
    # zoomLevels             2 bytes, unsigned
    # chromosomeTreeOffset   8 bytes, unsigned
    # fullDataOffset         8 bytes, unsigned
    # fullIndexOffset        8 bytes, unsigned
    # fieldCount             2 bytes, unsigned
    # definedFieldCount      2 bytes, unsigned
    # autoSqlOffset          8 bytes, unsigned
    # totalSummaryOffset     8 bytes, unsigned
    # uncompressBufSize      4 bytes, unsigned
    # reserved               8 bytes, unsigned

    formatter = struct.Struct("=IHHQQQHHQQIQ")
    size = formatter.size
    signature = 0x8789F2EB
    bbiCurrentVersion = 4

    @classmethod
    def fromfile(cls, stream):
        magic = stream.read(4)
        if int.from_bytes(magic, byteorder="little") == _Header.signature:
            byteorder = "<"
        elif int.from_bytes(magic, byteorder="big") == _Header.signature:
            byteorder = ">"
        else:
            raise ValueError("not a bigBed file")
        formatter = struct.Struct(byteorder + "HHQQQHHQQIQ")
        header = _Header()
        header.byteorder = byteorder
        size = formatter.size
        data = stream.read(size)
        (
            version,
            header.zoomLevels,
            header.chromosomeTreeOffset,
            header.fullDataOffset,
            header.fullIndexOffset,
            header.fieldCount,
            header.definedFieldCount,
            header.autoSqlOffset,
            header.totalSummaryOffset,
            header.uncompressBufSize,
            header.extraIndicesOffset,
        ) = formatter.unpack(data)
        assert version == _Header.bbiCurrentVersion
        definedFieldCount = header.definedFieldCount
        if definedFieldCount < 3 or definedFieldCount > 12:
            raise ValueError(
                "expected between 3 and 12 columns, found %d" % definedFieldCount
            )
        return header

    def __bytes__(self):
        return _Header.formatter.pack(
            _Header.signature,
            _Header.bbiCurrentVersion,
            self.zoomLevels,
            self.chromosomeTreeOffset,
            self.fullDataOffset,
            self.fullIndexOffset,
            self.fieldCount,
            self.definedFieldCount,
            self.autoSqlOffset,
            self.totalSummaryOffset,
            self.uncompressBufSize,
            self.extraIndicesOffset,
        )


class _ExtraIndex:
    __slots__ = ("indexField", "maxFieldSize", "fileOffset", "chunks", "get_value")

    formatter = struct.Struct("=xxHQxxxxHxx")

    def __init__(self, name, declaration):
        self.maxFieldSize = 0
        self.fileOffset = None
        for index, field in enumerate(declaration):
            if field.name == name:
                break
        else:
            raise ValueError(
                "extraIndex field %s not a standard bed field or found in 'as' file.",
                name,
            ) from None
        if field.as_type != "string":
            raise ValueError("Sorry for now can only index string fields.")
        self.indexField = index
        if name == "chrom":
            self.get_value = lambda alignment: alignment.target.id
        elif name == "name":
            self.get_value = lambda alignment: alignment.query.id
        else:
            self.get_value = lambda alignment: alignment.annotations[name]

    def updateMaxFieldSize(self, alignment):
        value = self.get_value(alignment)
        size = len(value)
        if size > self.maxFieldSize:
            self.maxFieldSize = size

    def addKeysFromRow(self, alignment, recordIx):
        value = self.get_value(alignment)
        self.chunks[recordIx]["name"] = value.encode()

    def addOffsetSize(self, offset, size, startIx, endIx):
        self.chunks[startIx:endIx]["offset"] = offset
        self.chunks[startIx:endIx]["size"] = size

    def __bytes__(self):
        indexFieldCount = 1
        return self.formatter.pack(indexFieldCount, self.fileOffset, self.indexField)


class _ExtraIndices(list):
    formatter = struct.Struct("=HHQ52x")

    def __init__(self, names, declaration):
        self[:] = [_ExtraIndex(name, declaration) for name in names]

    @property
    def size(self):
        return self.formatter.size + _ExtraIndex.formatter.size * len(self)

    def initialize(self, bedCount):
        if bedCount == 0:
            return
        for extra_index in self:
            keySize = extra_index.maxFieldSize
            # bbiFile.h: bbNamedFileChunk
            # name    keySize bytes
            # offset  8 bytes, unsigned
            # size    8 bytes, unsigned
            dtype = np.dtype(
                [("name", f"=S{keySize}"), ("offset", "=u8"), ("size", "=u8")]
            )
            extra_index.chunks = np.zeros(bedCount, dtype=dtype)

    def tofile(self, stream):
        size = self.formatter.size
        if len(self) > 0:
            offset = stream.tell() + size
            data = self.formatter.pack(size, len(self), offset)
            stream.write(data)
            for extra_index in self:
                stream.write(bytes(extra_index))
        else:
            data = self.formatter.pack(size, 0, 0)
            stream.write(data)


class _ZoomLevel:
    __slots__ = ["amount", "dataOffset", "indexOffset"]

    # Supplemental Table 6: The zoom header
    # reductionLevel   4 bytes, unsigned
    # reserved         4 bytes, unsigned
    # dataOffset       8 bytes, unsigned
    # indexOffset      8 bytes, unsigned

    formatter = struct.Struct("=IxxxxQQ")

    def __bytes__(self):
        return self.formatter.pack(self.amount, self.dataOffset, self.indexOffset)


class _ZoomLevels(list):
    bbiResIncrement = 4
    bbiMaxZoomLevels = 10
    size = _ZoomLevel.formatter.size * bbiMaxZoomLevels

    def __init__(self):
        self[:] = [_ZoomLevel() for i in range(_ZoomLevels.bbiMaxZoomLevels)]

    def __bytes__(self):
        data = b"".join(bytes(item) for item in self)
        data += bytes(_ZoomLevels.size - len(data))
        return data

    @classmethod
    def calculate_reductions(cls, aveSize):
        bbiMaxZoomLevels = _ZoomLevels.bbiMaxZoomLevels
        reductions = np.zeros(
            bbiMaxZoomLevels,
            dtype=[("scale", "=i4"), ("size", "=i4"), ("end", "=i4")],
        )
        minZoom = 10
        res = max(int(aveSize), minZoom)
        maxInt = np.iinfo(reductions.dtype["scale"]).max
        for resTry in range(bbiMaxZoomLevels):
            if res > maxInt:
                break
            reductions[resTry]["scale"] = res
            res *= _ZoomLevels.bbiResIncrement
        return reductions[:resTry]

    def reduce(self, summaries, initialReduction, buffer, blockSize, itemsPerSlot):
        zoomCount = initialReduction["size"]
        reduction = initialReduction["scale"] * _ZoomLevels.bbiResIncrement
        output = buffer.output
        formatter = _RTreeFormatter()
        for zoomLevels in range(1, _ZoomLevels.bbiMaxZoomLevels):
            rezoomCount = len(summaries)
            if rezoomCount >= zoomCount:
                break
            zoomCount = rezoomCount
            self[zoomLevels].dataOffset = output.tell()
            data = zoomCount.to_bytes(4, sys.byteorder)
            output.write(data)
            for summary in summaries:
                buffer.write(summary)
            buffer.flush()
            indexOffset = output.tell()
            formatter.write(summaries, blockSize, itemsPerSlot, indexOffset, output)
            self[zoomLevels].indexOffset = indexOffset
            self[zoomLevels].amount = reduction
            reduction *= _ZoomLevels.bbiResIncrement
            i = 0
            chromId = None
            for summary in summaries:
                if summary.chromId != chromId or summary.end > end:  # noqa: F821
                    end = summary.start + reduction
                    chromId = summary.chromId
                    currentSummary = summary
                    summaries[i] = currentSummary
                    i += 1
                else:
                    currentSummary += summary
            del summaries[i:]
        del self[zoomLevels:]


class _Summary:
    __slots__ = ("validCount", "minVal", "maxVal", "sumData", "sumSquares")

    formatter = struct.Struct("=Qdddd")
    size = formatter.size

    def __init__(self):
        self.validCount = 0
        self.minVal = sys.maxsize
        self.maxVal = -sys.maxsize
        self.sumData = 0.0
        self.sumSquares = 0.0

    def update(self, size, val):
        self.validCount += size
        if val < self.minVal:
            self.minVal = val
        if val > self.maxVal:
            self.maxVal = val
        self.sumData += val * size
        self.sumSquares += val * val * size

    def __bytes__(self):
        return self.formatter.pack(
            self.validCount,
            self.minVal,
            self.maxVal,
            self.sumData,
            self.sumSquares,
        )


class _Region:
    __slots__ = ("chromId", "start", "end", "offset")

    def __init__(self, chromId, start, end):
        self.chromId = chromId
        self.start = start
        self.end = end


class _RegionSummary(_Summary):
    __slots__ = _Region.__slots__ + _Summary.__slots__

    formatter = struct.Struct("=IIIIffff")
    size = formatter.size

    def __init__(self, chromId, start, end, value):
        self.chromId = chromId
        self.start = start
        self.end = end
        self.validCount = 0
        self.minVal = np.float32(value)
        self.maxVal = np.float32(value)
        self.sumData = np.float32(0.0)
        self.sumSquares = np.float32(0.0)
        self.offset = None

    def __iadd__(self, other):
        self.end = other.end
        self.validCount += other.validCount
        self.minVal = min(self.minVal, other.minVal)
        self.maxVal = max(self.maxVal, other.maxVal)
        self.sumData = np.float32(self.sumData + other.sumData)
        self.sumSquares = np.float32(self.sumSquares + other.sumSquares)
        return self

    def update(self, overlap, val):
        self.validCount += overlap
        if self.minVal > val:
            self.minVal = np.float32(val)
        if self.maxVal < val:
            self.maxVal = np.float32(val)
        self.sumData = np.float32(self.sumData + val * overlap)
        self.sumSquares = np.float32(self.sumSquares + val * val * overlap)

    def __bytes__(self):
        return self.formatter.pack(
            self.chromId,
            self.start,
            self.end,
            self.validCount,
            self.minVal,
            self.maxVal,
            self.sumData,
            self.sumSquares,
        )


class _RTreeNode:
    __slots__ = [
        "children",
        "parent",
        "startChromId",
        "startBase",
        "endChromId",
        "endBase",
        "startFileOffset",
        "endFileOffset",
    ]

    def __init__(self):
        self.parent = None
        self.children = []

    def calcLevelSizes(self, levelSizes, level):
        levelSizes[level] += 1
        level += 1
        if level == len(levelSizes):
            return
        for child in self.children:
            child.calcLevelSizes(levelSizes, level)


class _RangeTree:
    __slots__ = ("root", "n", "freeList", "stack", "chromId", "chromSize")

    def __init__(self, chromId, chromSize):
        self.root = None
        self.n = 0
        self.freeList = []
        self.chromId = chromId
        self.chromSize = chromSize

    @classmethod
    def generate(cls, chromUsageList, alignments):
        alignments.rewind()
        alignment = None
        for chromName, chromId, chromSize in chromUsageList:
            chromName = chromName.decode()
            tree = _RangeTree(chromId, chromSize)
            if alignment is not None:
                tree.addToCoverageDepth(alignment)
            for alignment in alignments:
                if alignment.target.id != chromName:
                    break
                tree.addToCoverageDepth(alignment)
            yield tree

    def generate_summaries(self, scale, totalSum):
        ranges = self.root.traverse()
        start, end, val = next(ranges)
        chromId = self.chromId
        chromSize = self.chromSize
        summary = _RegionSummary(
            chromId,
            start,
            min(start + scale, chromSize),
            val,
        )
        while True:
            size = max(end - start, 1)
            totalSum.update(size, val)
            if summary.end <= start:
                summary = _RegionSummary(
                    chromId,
                    start,
                    min(start + scale, chromSize),
                    val,
                )
            while end > summary.end:
                overlap = min(end, summary.end) - max(start, summary.start)
                assert overlap > 0
                summary.update(overlap, val)
                size -= overlap
                start = summary.end
                yield summary
                summary = _RegionSummary(
                    chromId,
                    start,
                    min(start + scale, chromSize),
                    val,
                )
            summary.update(size, val)
            try:
                start, end, val = next(ranges)
            except StopIteration:
                break
            if summary.end <= start:
                yield summary
        yield summary

    def find(self, start, end):
        p = self.root
        while p is not None:
            if end <= p.item.start:
                p = p.left
            elif p.item.end <= start:
                p = p.right
            else:
                return p.item

    def restructure(self, x, y, z):
        tos = len(self.stack)
        if y is x.left:
            if z is y.left:
                midNode = y
                y.left = z
                x.left = y.right
                y.right = x
            else:
                midNode = z
                y.right = z.left
                z.left = y
                x.left = z.right
                z.right = x
        else:
            if z is y.left:
                midNode = z
                x.right = z.left
                z.left = x
                y.left = z.right
                z.right = y
            else:
                midNode = y
                x.right = y.left
                y.left = x
                y.right = z
        if tos != 0:
            parent = self.stack[tos - 1]
            if x is parent.left:
                parent.left = midNode
            else:
                parent.right = midNode
        else:
            self.root = midNode
        return midNode

    def add(self, item):
        self.stack = []
        try:
            x = self.freeList.pop()
        except IndexError:
            x = _RedBlackTreeNode()
        else:
            self.freeList = x.right
        x.left = None
        x.right = None
        x.item = item
        x.color = None
        p = self.root
        if p is not None:
            while True:
                self.stack.append(p)
                if item.end <= p.item.start:
                    p = p.left
                    if p is None:
                        p = self.stack.pop()
                        p.left = x
                        break
                elif p.item.end <= item.start:
                    p = p.right
                    if p is None:
                        p = self.stack.pop()
                        p.right = x
                        break
                else:
                    return
            col = True
        else:
            self.root = x
            col = False
        x.color = col
        self.n += 1
        if len(self.stack) > 0:
            while p.color is True:
                m = self.stack.pop()
                if p == m.left:
                    q = m.right
                else:
                    q = m.left
                if q is None or q.color is False:
                    m = self.restructure(m, p, x)
                    m.color = False
                    m.left.color = True
                    m.right.color = True
                    break
                p.color = False
                q.color = False
                if len(self.stack) == 0:
                    break
                m.color = True
                x = m
                p = self.stack.pop()

    def addToCoverageDepth(self, alignment):
        start = alignment.coordinates[0, 0]
        end = alignment.coordinates[0, -1]
        if start > end:
            start, end = end, start
        existing = self.find(start, end)
        if existing is None:
            r = _Range(start, end, val=1)
            self.add(r)
        else:
            if existing.start <= start and existing.end >= end:
                if existing.start < start:
                    r = _Range(existing.start, start, existing.val)
                    existing.start = start
                    self.add(r)
                if existing.end > end:
                    r = _Range(end, existing.end, existing.val)
                    existing.end = end
                    self.add(r)
                existing.val += 1
            else:
                items = list(self.root.traverse_range(start, end))
                s = start
                e = end
                for item in items:
                    if s < item.start:
                        r = _Range(s, item.start, 1)
                        s = item.start
                        self.add(r)
                    elif s > item.start:
                        r = _Range(item.start, s, item.val)
                        item.start = s
                        self.add(r)
                    if item.start < end and item.end > end:
                        r = _Range(end, item.end, item.val)
                        item.end = end
                        self.add(r)
                    item.val += 1
                    s = item.end
                if s < e:
                    r = _Range(s, e, 1)
                    self.add(r)


class _Range:
    __slots__ = ("next", "start", "end", "val")

    def __init__(self, start, end, val):
        self.start = start
        self.end = end
        self.val = val

    def __iter__(self):
        return iter((self.start, self.end, self.val))


class _RedBlackTreeNode:
    __slots__ = ("left", "right", "color", "item")

    def traverse(self):
        if self.left is not None:
            yield from self.left.traverse()
        yield self.item
        if self.right is not None:
            yield from self.right.traverse()

    def traverse_range(self, start, end):
        if self.item.end <= start:
            if self.right is not None:
                yield from self.right.traverse_range(start, end)
        elif end <= self.item.start:
            if self.left is not None:
                yield from self.left.traverse_range(start, end)
        else:
            if self.left is not None:
                yield from self.left.traverse_range(start, end)
            yield self.item
            if self.right is not None:
                yield from self.right.traverse_range(start, end)


class _RTreeFormatter:
    signature = 0x2468ACE0

    def __init__(self, byteorder="="):
        # Supplemental Table 14: R tree index header
        # magic          4 bytes, unsigned
        # blockSize      4 bytes, unsigned
        # itemCount      8 bytes, unsigned
        # startChromIx   4 bytes, unsigned
        # startBase      4 bytes, unsigned
        # endChromIx     4 bytes, unsigned
        # endBase        4 bytes, unsigned
        # endFileOffset  8 bytes, unsigned
        # itemsPerSlot   4 bytes, unsigned
        # reserved       4 bytes, unsigned
        self.formatter_header = struct.Struct(byteorder + "IIQIIIIQIxxxx")

        # Supplemental Table 15: R tree node format
        # isLeaf    1 byte
        # reserved  1 byte
        # count     2 bytes, unsigned
        self.formatter_node = struct.Struct(byteorder + "?xH")

        # Supplemental Table 17: R tree non-leaf format
        # startChromIx   4 bytes, unsigned
        # startBase      4 bytes, unsigned
        # endChromIx     4 bytes, unsigned
        # endBase        4 bytes, unsigned
        # dataOffset     8 bytes, unsigned
        self.formatter_nonleaf = struct.Struct(byteorder + "IIIIQ")

        # Supplemental Table 16: R tree leaf format
        # startChromIx   4 bytes, unsigned
        # startBase      4 bytes, unsigned
        # endChromIx     4 bytes, unsigned
        # endBase        4 bytes, unsigned
        # dataOffset     8 bytes, unsigned
        # dataSize       8 bytes, unsigned
        self.formatter_leaf = struct.Struct(byteorder + "IIIIQQ")

    def read(self, stream):
        NonLeaf = namedtuple(
            "NonLeaf",
            [
                "parent",
                "children",
                "startChromIx",
                "startBase",
                "endChromIx",
                "endBase",
                "dataOffset",
            ],
        )

        Leaf = namedtuple(
            "Leaf",
            [
                "parent",
                "startChromIx",
                "startBase",
                "endChromIx",
                "endBase",
                "dataOffset",
                "dataSize",
            ],
        )

        data = stream.read(self.formatter_header.size)
        (
            magic,
            blockSize,
            itemCount,
            startChromIx,
            startBase,
            endChromIx,
            endBase,
            endFileOffset,
            itemsPerSlot,
        ) = self.formatter_header.unpack(data)
        assert magic == _RTreeFormatter.signature

        formatter_node = self.formatter_node
        formatter_nonleaf = self.formatter_nonleaf
        formatter_leaf = self.formatter_leaf

        root = NonLeaf(None, [], startChromIx, startBase, endChromIx, endBase, None)
        node = root
        itemsCounted = 0
        while True:
            data = stream.read(formatter_node.size)
            isLeaf, count = formatter_node.unpack(data)
            if isLeaf:
                children = node.children
                for i in range(count):
                    data = stream.read(formatter_leaf.size)
                    (
                        startChromIx,
                        startBase,
                        endChromIx,
                        endBase,
                        dataOffset,
                        dataSize,
                    ) = formatter_leaf.unpack(data)
                    child = Leaf(
                        node,
                        startChromIx,
                        startBase,
                        endChromIx,
                        endBase,
                        dataOffset,
                        dataSize,
                    )
                    children.append(child)
                itemsCounted += count
                while True:
                    parent = node.parent
                    if parent is None:
                        assert itemsCounted == itemCount
                        return node
                    for index, child in enumerate(parent.children):
                        if id(node) == id(child):
                            break
                    else:
                        raise RuntimeError("Failed to find child node")
                    try:
                        node = parent.children[index + 1]
                    except IndexError:
                        node = parent
                    else:
                        break
            else:
                children = node.children
                for i in range(count):
                    data = stream.read(formatter_nonleaf.size)
                    (
                        startChromIx,
                        startBase,
                        endChromIx,
                        endBase,
                        dataOffset,
                    ) = formatter_nonleaf.unpack(data)
                    child = NonLeaf(
                        node,
                        [],
                        startChromIx,
                        startBase,
                        endChromIx,
                        endBase,
                        dataOffset,
                    )
                    children.append(child)
                parent = node
                node = children[0]
            stream.seek(node.dataOffset)

    def rTreeFromChromRangeArray(self, blockSize, items, endFileOffset):
        itemCount = len(items)
        if itemCount == 0:
            return
        children = []
        nextOffset = items[0].offset
        oneSize = 0
        i = 0
        while i < itemCount:
            child = _RTreeNode()
            children.append(child)
            startItem = items[i]
            child.startChromId = child.endChromId = startItem.chromId
            child.startBase = startItem.start
            child.endBase = startItem.end
            child.startFileOffset = nextOffset
            oneSize = 1
            endItem = startItem
            for j in range(i + 1, itemCount):
                endItem = items[j]
                nextOffset = endItem.offset
                if nextOffset != child.startFileOffset:
                    break
                oneSize += 1
            else:
                nextOffset = endFileOffset
            child.endFileOffset = nextOffset
            for item in items[i + 1 : i + oneSize]:
                if item.chromId < child.startChromId:
                    child.startChromId = item.chromId
                    child.startBase = item.start
                elif (
                    item.chromId == child.startChromId and item.start < child.startBase
                ):
                    child.startBase = item.start
                if item.chromId > child.endChromId:
                    child.endChromId = item.chromId
                    child.endBase = item.end
                elif item.chromId == child.endChromId and item.end > child.endBase:
                    child.endBase = item.end
            i += oneSize
        levelCount = 1
        while True:
            parents = []
            slotsUsed = blockSize
            for child in children:
                if slotsUsed >= blockSize:
                    slotsUsed = 1
                    parent = _RTreeNode()
                    parent.parent = child.parent
                    parent.startChromId = child.startChromId
                    parent.startBase = child.startBase
                    parent.endChromId = child.endChromId
                    parent.endBase = child.endBase
                    parent.startFileOffset = child.startFileOffset
                    parent.endFileOffset = child.endFileOffset
                    parents.append(parent)
                else:
                    slotsUsed += 1
                    if child.startChromId < parent.startChromId:
                        parent.startChromId = child.startChromId
                        parent.startBase = child.startBase
                    elif (
                        child.startChromId == parent.startChromId
                        and child.startBase < parent.startBase
                    ):
                        parent.startBase = child.startBase
                    if child.endChromId > parent.endChromId:
                        parent.endChromId = child.endChromId
                        parent.endBase = child.endBase
                    elif (
                        child.endChromId == parent.endChromId
                        and child.endBase > parent.endBase
                    ):
                        parent.endBase = child.endBase
                parent.children.append(child)
                child.parent = parent
            levelCount += 1
            if len(parents) == 1:
                break
            children = parents
        return parent, levelCount

    def rWriteLeaves(self, itemsPerSlot, lNodeSize, tree, curLevel, leafLevel, output):
        formatter_leaf = self.formatter_leaf
        if curLevel == leafLevel:
            isLeaf = True
            data = self.formatter_node.pack(isLeaf, len(tree.children))
            output.write(data)
            for child in tree.children:
                data = formatter_leaf.pack(
                    child.startChromId,
                    child.startBase,
                    child.endChromId,
                    child.endBase,
                    child.startFileOffset,
                    child.endFileOffset - child.startFileOffset,
                )
                output.write(data)
            output.write(
                bytes((itemsPerSlot - len(tree.children)) * self.formatter_nonleaf.size)
            )
            # self.formatter_leaf.size seems more reasonable here, but this is
            # what bedToBigBed has here. See also
            # https://github.com/ucscGenomeBrowser/kent/issues/76.
        else:
            for child in tree.children:
                self.rWriteLeaves(
                    itemsPerSlot, lNodeSize, child, curLevel + 1, leafLevel, output
                )

    def rWriteIndexLevel(
        self, parent, blockSize, childNodeSize, curLevel, destLevel, offset, output
    ):
        previous_offset = offset
        formatter_nonleaf = self.formatter_nonleaf
        if curLevel == destLevel:
            isLeaf = False
            data = self.formatter_node.pack(isLeaf, len(parent.children))
            output.write(data)
            for child in parent.children:
                data = formatter_nonleaf.pack(
                    child.startChromId,
                    child.startBase,
                    child.endChromId,
                    child.endBase,
                    offset,
                )
                output.write(data)
                offset += childNodeSize
            output.write(
                bytes((blockSize - len(parent.children)) * self.formatter_nonleaf.size)
            )
        else:
            for child in parent.children:
                offset = self.rWriteIndexLevel(
                    child,
                    blockSize,
                    childNodeSize,
                    curLevel + 1,
                    destLevel,
                    offset,
                    output,
                )
        position = output.tell()
        if position != previous_offset:
            raise RuntimeError(
                f"Internal error: offset mismatch ({position} vs {previous_offset})"
            )
        return offset

    def write(self, items, blockSize, itemsPerSlot, endFileOffset, output):
        root, levelCount = self.rTreeFromChromRangeArray(
            blockSize, items, endFileOffset
        )

        data = self.formatter_header.pack(
            _RTreeFormatter.signature,
            blockSize,
            len(items),
            root.startChromId,
            root.startBase,
            root.endChromId,
            root.endBase,
            endFileOffset,
            itemsPerSlot,
        )
        output.write(data)

        if root is None:
            return

        levelSizes = np.zeros(levelCount, int)
        root.calcLevelSizes(levelSizes, level=0)
        size = self.formatter_node.size + self.formatter_nonleaf.size * blockSize
        levelOffset = output.tell()
        for i in range(levelCount - 2):
            levelOffset += levelSizes[i] * size
            if i == levelCount - 3:
                size = self.formatter_node.size + self.formatter_leaf.size * blockSize
            self.rWriteIndexLevel(root, blockSize, size, 0, i, levelOffset, output)
        leafLevel = levelCount - 2
        self.rWriteLeaves(blockSize, size, root, 0, leafLevel, output)


class _BPlusTreeFormatter:
    signature = 0x78CA8C91

    def __init__(self, byteorder="="):
        # Supplemental Table 8: Chromosome B+ tree header
        # magic     4 bytes, unsigned
        # blockSize 4 bytes, unsigned
        # keySize   4 bytes, unsigned
        # valSize   4 bytes, unsigned
        # itemCount 8 bytes, unsigned
        # reserved  8 bytes, unsigned
        self.formatter_header = struct.Struct(byteorder + "IIIIQxxxxxxxx")

        # Supplemental Table 9: Chromosome B+ tree node
        # isLeaf    1 byte
        # reserved  1 byte
        # count     2 bytes, unsigned
        self.formatter_node = struct.Struct(byteorder + "?xH")

        # Supplemental Table 11: Chromosome B+ tree non-leaf item
        # key          keySize bytes
        # childOffset  8 bytes, unsigned
        self.fmt_nonleaf = byteorder + "{keySize}sQ"

        self.byteorder = byteorder

    def read(self, stream):
        byteorder = self.byteorder

        formatter = self.formatter_header
        data = stream.read(formatter.size)
        magic, blockSize, keySize, valSize, itemCount = formatter.unpack(data)
        assert magic == _BPlusTreeFormatter.signature

        formatter_node = self.formatter_node
        formatter_nonleaf = struct.Struct(self.fmt_nonleaf.format(keySize=keySize))

        # Supplemental Table 10: Chromosome B+ tree leaf item format
        # key        keySize bytes
        # chromId    4 bytes, unsigned
        # chromSize  4 bytes, unsigned
        formatter_leaf = struct.Struct(f"{byteorder}{keySize}sII")
        assert keySize == formatter_leaf.size - valSize
        assert valSize == 8

        Node = namedtuple("Node", ["parent", "children"])

        targets = []
        node = None
        while True:
            data = stream.read(formatter_node.size)
            isLeaf, count = formatter_node.unpack(data)
            if isLeaf:
                for i in range(count):
                    data = stream.read(formatter_leaf.size)
                    key, chromId, chromSize = formatter_leaf.unpack(data)
                    name = key.rstrip(b"\x00").decode()
                    assert chromId == len(targets)
                    sequence = Seq(None, length=chromSize)
                    record = SeqRecord(sequence, id=name)
                    targets.append(record)
            else:
                children = []
                for i in range(count):
                    data = stream.read(formatter_nonleaf.size)
                    key, pos = formatter_nonleaf.unpack(data)
                    children.append(pos)
                parent = node
                node = Node(parent, children)
            while True:
                if node is None:
                    assert len(targets) == itemCount
                    return targets
                children = node.children
                try:
                    pos = children.pop(0)
                except IndexError:
                    node = node.parent
                else:
                    break
            stream.seek(pos)

    def write(self, items, blockSize, output):
        signature = _BPlusTreeFormatter.signature
        keySize = items.dtype["name"].itemsize
        valSize = items.itemsize - keySize
        itemCount = len(items)
        formatter = self.formatter_header
        data = formatter.pack(signature, blockSize, keySize, valSize, itemCount)
        output.write(data)

        formatter_node = self.formatter_node
        formatter_nonleaf = struct.Struct(self.fmt_nonleaf.format(keySize=keySize))

        levels = 1
        while itemCount > blockSize:
            itemCount = len(range(0, itemCount, blockSize))
            levels += 1
        itemCount = len(items)
        bytesInIndexBlock = formatter_node.size + blockSize * formatter_nonleaf.size
        bytesInLeafBlock = formatter_node.size + blockSize * items.itemsize
        isLeaf = False
        indexOffset = output.tell()
        for level in range(levels - 1, 0, -1):
            slotSizePer = blockSize**level
            nodeSizePer = slotSizePer * blockSize
            indices = range(0, itemCount, nodeSizePer)
            if level == 1:
                bytesInNextLevelBlock = bytesInLeafBlock
            else:
                bytesInNextLevelBlock = bytesInIndexBlock
            levelSize = len(indices) * bytesInIndexBlock
            endLevel = indexOffset + levelSize
            nextChild = endLevel
            for index in indices:
                block = items[index : index + nodeSizePer : slotSizePer]
                n = len(block)
                output.write(formatter_node.pack(isLeaf, n))
                for item in block:
                    data = formatter_nonleaf.pack(item["name"], nextChild)
                    output.write(data)
                    nextChild += bytesInNextLevelBlock
                data = bytes((blockSize - n) * formatter_nonleaf.size)
                output.write(data)
            indexOffset = endLevel
        isLeaf = True
        for index in itertools.count(0, blockSize):
            block = items[index : index + blockSize]
            n = len(block)
            if n == 0:
                break
            output.write(formatter_node.pack(isLeaf, n))
            block.tofile(output)
            data = bytes((blockSize - n) * items.itemsize)
            output.write(data)
