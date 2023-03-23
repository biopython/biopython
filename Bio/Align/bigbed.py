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

# The writer is based on the following files of the bedToBigBed program:
# bedToBigBed.c
# bbiWrite.c

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
import copy
import array
import struct
import zlib
from collections import namedtuple
from io import BytesIO
from operator import attrgetter
import numpy


from Bio.Align import Alignment
from Bio.Align import interfaces
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# flake8: noqa


class bbiSummaryElement:
    # See bbiSummaryElementWrite in bbiWrite.c
    __slots__ = ["validCount", "minVal", "maxVal", "sumData", "sumSquares"]

    formatter = struct.Struct("Qdddd")
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


class ExtHeader:
    # See bbFileCreate in bedToBigBed.c
    __slots__ = ("extraIndexCount", "extraIndexListOffset")

    formatter = struct.Struct("=HHQ52x")
    size = formatter.size

    def __bytes__(self):
        return self.formatter.pack(
            self.size, self.extraIndexCount, self.extraIndexListOffset
        )


class ExtraIndex:
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
        self.chunks[recordIx].name = value.encode()

    def addOffsetSize(self, offset, size, startIx, endIx):
        for chunk in self.chunks[startIx:endIx]:
            chunk.offset = offset
            chunk.size = size

    def __bytes__(self):
        indexFieldCount = 1
        return self.formatter.pack(indexFieldCount, self.fileOffset, self.indexField)


Field = namedtuple("Field", ("as_type", "name", "comment"))


class AutoSQLTable(list):
    """AutoSQL table describing the columns of an (possibly extended) BED format."""

    def __init__(self, name, comment, fields):
        """Create an AutoSQL table describing the columns of an (extended) BED format."""
        self.name = name
        self.comment = comment
        self[:] = fields

    @classmethod
    def from_bytes(cls, data):
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
                '   %s %s    "%s"\n'
                % (
                    field.as_type.ljust(type_width),
                    name.ljust(name_width),
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
        extraIndex=[],
    ):
        """Create an AlignmentWriter object.
        Arguments:
         - target      - output stream or file name
         - bedN        - number of columns in the BED file.
                         This must be between 3 and 12; default value is 12.
         - declaration - an AutoSQLTable object declaring the fields in the BED file.
                         Required only if the BED file contains extra (custom) fields.
                         Default value is None.
         - targets     - A list of SeqRecord objects with the chromosomes in the
                         order as they appear in the alignments. The sequence
                         contents in each SeqRecord may be undefined, but the
                         sequence length must be defined, as in this example:
                         SeqRecord(Seq(None, length=248956422), id="chr1").
                         If targets is None (the default value), the alignments
                         must have an attribute .targets providing the list of
                         SeqRecord objects.
         - compress    - If True (default), compress data using zlib.
                         If False, do not compress data.
         - extraIndex  - List of strings with the names of extra columns to be indexed.
                         Default value is an empty list.
        """
        if bedN < 3 or bedN > 12:
            raise ValueError("bedN must be between 3 and 12")
        super().__init__(target)
        self.bedN = bedN
        self.declaration = declaration
        self.targets = targets
        self.compress = compress
        self.extraIndex = extraIndex

    def write_file(self, stream, alignments):
        blockSize = 256
        itemsPerSlot = 512
        if self.targets is None:
            targets = alignments.targets
        else:
            targets = self.targets
        bedN = self.bedN
        if self.declaration is None:
            try:
                declaration = alignments.declaration[:bedN]
            except AttributeError:
                declaration = AutoSQLTable.default[:bedN]
        else:
            declaration = self.declaration
        extraIndex = self.extraIndex
        compress = self.compress
        # see bbFileCreate in bedToBigBed.c
        fieldCount = len(declaration)
        extHeader = ExtHeader()
        extHeader.extraIndexCount = len(extraIndex)
        eim = [ExtraIndex(name, declaration) for name in extraIndex]
        usageList, minDiff, aveSize, bedCount = bbiChromUsageFromBedFile(
            alignments, targets, eim
        )
        bbiCurrentVersion = 4
        stream.write(bytes(64))  # bbiWriteDummyHeader
        stream.write(bytes(bbiMaxZoomLevels * 24))  # bbiWriteDummyZooms
        asOffset = stream.tell()
        stream.write(bytes(declaration))  # asText
        totalSummaryOffset = stream.tell()
        stream.write(bytes(bbiSummaryElement.size))
        extHeaderOffset = stream.tell()
        stream.write(bytes(extHeader.size))
        if extraIndex:
            extHeader.extraIndexListOffset = stream.tell()
            extraIndexSize = 16 + 4 * 1
            stream.write(bytes(extraIndexSize * extHeader.extraIndexCount))
            extraIndexListEndOffset = stream.tell()
        else:
            extHeader.extraIndexListOffset = 0
            extraIndexListEndOffset = 0
        chromTreeOffset = stream.tell()
        bbiWriteChromInfo(usageList, blockSize, stream)
        dataOffset = stream.tell()
        resTryCount, resScales, resSizes = bbiCalcResScalesAndSizes(aveSize)
        stream.write(struct.pack("Q", bedCount))
        if bedCount > 0:
            for element in eim:
                element.chunks = [bbNamedFileChunk() for j in range(bedCount)]
            maxBlockSize, boundsArray = writeBlocks(
                alignments,
                declaration,
                itemsPerSlot,
                compress,
                stream,
                resTryCount,
                resScales,
                resSizes,
                eim,
                bedCount,
                fieldCount,
                bedN,
            )
        else:
            maxBlockSize = 0
            boundsArray = []
        indexOffset = stream.tell()
        cirTreeFileBulkIndexToOpenFile(boundsArray, blockSize, 1, indexOffset, stream)
        zoomAmounts = numpy.empty(bbiMaxZoomLevels, numpy.int32)
        zoomDataOffsets = numpy.empty(bbiMaxZoomLevels, numpy.int64)
        zoomIndexOffsets = numpy.empty(bbiMaxZoomLevels, numpy.int64)
        zoomLevels = 0
        if bedCount > 0:
            zoomLevels, totalSum = bbiWriteZoomLevels(
                alignments,
                stream,
                blockSize,
                itemsPerSlot,
                fieldCount,
                compress,
                indexOffset - dataOffset,
                usageList,
                resTryCount,
                resScales,
                resSizes,
                zoomAmounts,
                zoomDataOffsets,
                zoomIndexOffsets,
            )
        for element in eim:
            element.fileOffset = stream.tell()
            maxBedNameSize = element.maxFieldSize
            element.chunks.sort()
            bptFileBulkIndexToOpenFile(
                element.chunks,
                bedCount,
                blockSize,
                maxBedNameSize,
                stream,
            )
        if compress:
            uncompressBufSize = max(maxBlockSize, itemsPerSlot * bbiSummary.size)
        else:
            uncompressBufSize = 0
        stream.seek(0)
        signature = 0x8789F2EB
        data = struct.pack(
            "=IHHQQQHHQQIQ",
            signature,
            bbiCurrentVersion,
            zoomLevels,
            chromTreeOffset,
            dataOffset,
            indexOffset,
            fieldCount,
            bedN,
            asOffset,
            totalSummaryOffset,
            uncompressBufSize,
            extHeaderOffset,
        )
        stream.write(data)
        for i in range(zoomLevels):
            data = struct.pack(
                "IxxxxQQ", zoomAmounts[i], zoomDataOffsets[i], zoomIndexOffsets[i]
            )
            stream.write(data)
        stream.write(bytes(24 * (bbiMaxZoomLevels - zoomLevels)))
        stream.seek(totalSummaryOffset)
        stream.write(bytes(totalSum))
        stream.seek(extHeaderOffset)
        stream.write(bytes(extHeader))
        if extraIndex:
            stream.seek(extHeader.extraIndexListOffset)
            for element in eim:
                stream.write(bytes(element))
            assert stream.tell() == extraIndexListEndOffset
        stream.seek(0, 2)
        data = struct.pack("I", signature)
        stream.write(data)


class AlignmentIterator(interfaces.AlignmentIterator):
    """Alignment iterator for bigBed files.

    The pairwise alignments stored in the bigBed file are loaded and returned
    incrementally.  Additional alignment information is stored as attributes
    of each alignment.
    """

    fmt = "bigBed"
    mode = "b"

    def _read_header(self, stream):
        # Supplemental Table 5: Common header
        # magic                  4 bytes, unsigned
        # version                2 bytes, unsigned
        # zoomLevels             2 bytes, unsigned
        # chromosomeTreeOffset   8 bytes, unsigned
        # fullDataOffset         8 bytes, unsigned; points to dataCount
        # fullIndexOffset        8 bytes, unsigned
        # fieldCount             2 bytes, unsigned
        # definedFieldCount      2 bytes, unsigned
        # autoSqlOffset          8 bytes, unsigned
        # totalSummaryOffset     8 bytes, unsigned
        # uncompressBufSize      4 bytes, unsigned
        # reserved               8 bytes, unsigned
        signature = 0x8789F2EB
        magic = stream.read(4)
        for byteorder in ("little", "big"):
            if int.from_bytes(magic, byteorder=byteorder) == signature:
                break
        else:
            raise ValueError("not a bigBed file")
        self.byteorder = byteorder
        if byteorder == "little":
            byteorder_char = "<"
        elif byteorder == "big":
            byteorder_char = ">"
        else:
            raise ValueError("Unexpected byteorder '%s'" % byteorder)

        (
            version,
            zoomLevels,
            chromosomeTreeOffset,
            fullDataOffset,
            fullIndexOffset,
            fieldCount,
            definedFieldCount,
            autoSqlOffset,
            totalSummaryOffset,
            uncompressBufSize,
        ) = struct.unpack(byteorder_char + "HHQQQHHQQIxxxxxxxx", stream.read(60))

        autoSqlSize = totalSummaryOffset - autoSqlOffset
        self.declaration = self._read_autosql(
            stream, autoSqlOffset, autoSqlSize, fieldCount, definedFieldCount
        )

        stream.seek(fullDataOffset)
        dataCount = int.from_bytes(stream.read(8), byteorder=byteorder)
        self._length = dataCount

        if uncompressBufSize > 0:
            self._compressed = True
        else:
            self._compressed = False

        self.targets = self._read_chromosomes(stream, chromosomeTreeOffset)

        self.tree = self._read_index(stream, fullIndexOffset)

        self._data = self._iterate_index(stream)

    def _read_autosql(self, stream, pos, size, fieldCount, definedFieldCount):
        if definedFieldCount < 3 or definedFieldCount > 12:
            raise ValueError(
                "expected between 3 and 12 columns, found %d" % definedFieldCount
            )
        self.bedN = definedFieldCount
        stream.seek(pos)
        data = stream.read(size)
        declaration = AutoSQLTable.from_bytes(data)
        self._analyze_fields(declaration, fieldCount, definedFieldCount)
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

    def _read_chromosomes(self, stream, pos):
        byteorder = self.byteorder
        if byteorder == "little":
            byteorder_char = "<"
        elif byteorder == "big":
            byteorder_char = ">"
        else:
            raise ValueError("Unexpected byteorder '%s'" % byteorder)
        # Supplemental Table 8: Chromosome B+ tree header
        # magic     4 bytes, unsigned
        # blockSize 4 bytes, unsigned
        # keySize   4 bytes, unsigned
        # valSize   4 bytes, unsigned
        # itemCount 8 bytes, unsigned
        # reserved  8 bytes, unsigned
        stream.seek(pos)
        signature = 0x78CA8C91
        magic = int.from_bytes(stream.read(4), byteorder=byteorder)
        assert magic == signature
        blockSize, keySize, valSize, itemCount = struct.unpack(
            byteorder_char + "IIIQxxxxxxxx", stream.read(28)
        )
        assert valSize == 8

        Node = namedtuple("Node", ["parent", "children"])
        targets = []
        node = None
        while True:
            # Supplemental Table 9: Chromosome B+ tree node
            # isLeaf    1 byte
            # reserved  1 byte
            # count     2 bytes, unsigned
            isLeaf, count = struct.unpack(byteorder_char + "?xH", stream.read(4))
            if isLeaf:
                for i in range(count):
                    # Supplemental Table 10: Chromosome B+ tree leaf item format
                    # key        keySize bytes
                    # chromId    4 bytes, unsigned
                    # chromSize  4 bytes, unsigned
                    key = stream.read(keySize)
                    name = key.rstrip(b"\x00").decode()
                    chromId, chromSize = struct.unpack(
                        byteorder_char + "II", stream.read(valSize)
                    )
                    assert chromId == len(targets)
                    sequence = Seq(None, length=chromSize)
                    record = SeqRecord(sequence, id=name)
                    targets.append(record)
            else:
                children = []
                for i in range(count):
                    # Supplemental Table 11: Chromosome B+ tree non-leaf item
                    # key          keySize bytes
                    # childOffset  8 bytes, unsigned
                    key = stream.read(keySize)
                    pos = int.from_bytes(stream.read(8), byteorder)
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

    def _read_index(self, stream, pos):
        byteorder = self.byteorder
        if byteorder == "little":
            byteorder_char = "<"
        elif byteorder == "big":
            byteorder_char = ">"
        else:
            raise ValueError("Unexpected byteorder '%s'" % byteorder)
        Node = namedtuple(
            "Node",
            [
                "parent",
                "children",
                "startChromIx",
                "startBase",
                "endChromIx",
                "endBase",
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
        stream.seek(pos)
        signature = 0x2468ACE0
        magic = int.from_bytes(stream.read(4), byteorder=byteorder)
        assert magic == signature
        (
            blockSize,
            itemCount,
            startChromIx,
            startBase,
            endChromIx,
            endBase,
            endFileOffset,
            itemsPerSlot,
        ) = struct.unpack(byteorder_char + "IQIIIIQIxxxx", stream.read(44))
        root = Node(None, [], startChromIx, startBase, endChromIx, endBase)
        node = root
        itemsCounted = 0
        dataOffsets = {}
        while True:
            # Supplemental Table 15: R tree node format
            # isLeaf    1 byte
            # reserved  1 byte
            # count     2 bytes, unsigned
            isLeaf, count = struct.unpack(byteorder_char + "?xH", stream.read(4))
            if isLeaf:
                children = node.children
                for i in range(count):
                    # Supplemental Table 16: R tree leaf format
                    # startChromIx   4 bytes, unsigned
                    # startBase      4 bytes, unsigned
                    # endChromIx     4 bytes, unsigned
                    # endBase        4 bytes, unsigned
                    # dataOffset     8 bytes, unsigned
                    # dataSize       8 bytes, unsigned
                    (
                        startChromIx,
                        startBase,
                        endChromIx,
                        endBase,
                        dataOffset,
                        dataSize,
                    ) = struct.unpack(byteorder_char + "IIIIQQ", stream.read(32))
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
                        assert not dataOffsets
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
                    # Supplemental Table 17: R tree non-leaf format
                    # startChromIx   4 bytes, unsigned
                    # startBase      4 bytes, unsigned
                    # endChromIx     4 bytes, unsigned
                    # endBase        4 bytes, unsigned
                    # dataOffset     8 bytes, unsigned
                    (
                        startChromIx,
                        startBase,
                        endChromIx,
                        endBase,
                        dataOffset,
                    ) = struct.unpack(byteorder_char + "IIIIQ", stream.read(24))
                    child = Node(node, [], startChromIx, startBase, endChromIx, endBase)
                    dataOffsets[id(child)] = dataOffset
                    children.append(child)
                parent = node
                node = children[0]
            pos = dataOffsets.pop(id(node))
            stream.seek(pos)

    def _iterate_index(self, stream):
        byteorder = self.byteorder
        if byteorder == "little":
            byteorder_char = "<"
        elif byteorder == "big":
            byteorder_char = ">"
        else:
            raise ValueError("Unexpected byteorder '%s'" % byteorder)
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
                    # Supplemental Table 12: Binary BED-data format
                    # chromId     4 bytes, unsigned
                    # chromStart  4 bytes, unsigned
                    # chromEnd    4 bytes, unsigned
                    # rest        zero-terminated string in tab-separated format
                    chromId, chromStart, chromEnd = struct.unpack(
                        byteorder_char + "III", data[:12]
                    )
                    rest, data = data[12:].split(b"\00", 1)
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
        byteorder = self.byteorder
        if byteorder == "little":
            byteorder_char = "<"
        elif byteorder == "big":
            byteorder_char = ">"
        else:
            raise ValueError("Unexpected byteorder '%s'" % byteorder)
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
                    # Supplemental Table 12: Binary BED-data format
                    # chromId     4 bytes, unsigned
                    # chromStart  4 bytes, unsigned
                    # chromEnd    4 bytes, unsigned
                    # rest        zero-terminated string in tab-separated format
                    child_chromIx, child_chromStart, child_chromEnd = struct.unpack(
                        byteorder_char + "III", data[:12]
                    )
                    rest, data = data[12:].split(b"\00", 1)
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
            chunk = next(self._data)
        except StopIteration:
            return
        return self._create_alignment(chunk)

    def _create_alignment(self, chunk):
        chromId, chromStart, chromEnd, rest = chunk
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
            blockSizes = [
                int(blockSize) for blockSize in words[7].rstrip(",").split(",")
            ]
            blockStarts = [
                int(blockStart) for blockStart in words[8].rstrip(",").split(",")
            ]
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
            blockSizes = numpy.array(blockSizes)
            blockStarts = numpy.array(blockStarts)
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
            coordinates = numpy.array(coordinates).transpose()
            qSize = sum(blockSizes)
        else:
            blockSize = chromEnd - chromStart
            coordinates = numpy.array([[0, blockSize], [0, blockSize]])
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
        for chunk in data:
            alignment = self._create_alignment(chunk)
            yield alignment


bbiMaxZoomLevels = 10
bbiResIncrement = 4

nodeHeaderSize = 4
indexSlotSize = 24
leafSlotSize = 32
bptBlockHeaderSize = 4

bbiBoundsArray = namedtuple("bbiBoundsArray", ["offset", "chromId", "start", "end"])

bbiChromUsage = namedtuple("bbiChromUsage", ["name", "itemCount", "id", "size"])


class bbiSummary:
    formatter = struct.Struct("=IIIIffff")
    size = formatter.size

    def __init__(self, chromId, start, end, value):
        self.chromId = chromId
        self.start = start
        self.end = end
        self.validCount = 0
        self.values = numpy.array(
            (value, value, 0.0, 0.0),
            dtype=[
                ("minVal", "f4"),
                ("maxVal", "f4"),
                ("sumData", "f4"),
                ("sumSquares", "f4"),
            ],
        )
        self.offset = 0

    def update(self, overlap, val):
        self.validCount += overlap
        if self.values["minVal"] > val:
            self.values["minVal"] = val
        if self.values["maxVal"] < val:
            self.values["maxVal"] = val
        self.values["sumData"] += val * overlap
        self.values["sumSquares"] += val * val * overlap

    def __bytes__(self):
        return self.formatter.pack(
            self.chromId,
            self.start,
            self.end,
            self.validCount,
            self.values["minVal"],
            self.values["maxVal"],
            self.values["sumData"],
            self.values["sumSquares"],
        )


class bbiSumOutStream:
    __slots__ = ["buffer", "elCount", "allocCount", "output", "doCompress"]

    def __init__(self, allocCount, output, doCompress):
        self.buffer = BytesIO()
        self.elCount = 0
        self.allocCount = allocCount
        self.output = output
        self.doCompress = doCompress

    def write(self, summary):
        data = struct.pack(
            "IIIIffff",
            summary.chromId,
            summary.start,
            summary.end,
            summary.validCount,
            summary.values["minVal"],
            summary.values["maxVal"],
            summary.values["sumData"],
            summary.values["sumSquares"],
        )
        self.buffer.write(data)
        self.elCount += 1
        if self.elCount >= self.allocCount:
            self.flush()

    def flush(self):
        data = self.buffer.getvalue()
        if data:
            assert self.elCount > 0
            if self.doCompress:
                data = zlib.compress(data)
            self.output.write(data)
            self.buffer = BytesIO()
        else:
            assert self.elCount == 0
        self.elCount = 0


class rTree:
    __slots__ = [
        "children",  # rTree*
        "parent",  # rTree*
        "startChromId",  # bits32, unsigned
        "startBase",  # bits32, unsigned
        "endChromId",  # bits32, unsigned
        "endBase",  # bits32, unsigned
        "startFileOffset",  # bits64
        "endFileOffset",  # bits64
    ]

    formatter = struct.Struct("=IIII")
    size = formatter.size

    def __init__(self):
        self.parent = None
        self.children = None

    def __bytes__(self):
        return self.formatter.pack(
            self.startChromId, self.startBase, self.endChromId, self.endBase
        )

    def clone(self):
        tree = rTree()
        tree.children = self.children
        tree.parent = self.parent
        tree.startChromId = self.startChromId
        tree.startBase = self.startBase
        tree.endChromId = self.endChromId
        tree.endBase = self.endBase
        tree.startFileOffset = self.startFileOffset
        tree.endFileOffset = self.endFileOffset
        return tree

    def rWriteIndexLevel(
        self, blockSize, childNodeSize, curLevel, destLevel, offset, output
    ):
        # in cirTree.c
        if curLevel == destLevel:
            byteorder = sys.byteorder
            countOne = len(self.children)
            isLeaf = False
            output.write(struct.pack("?xH", isLeaf, countOne))
            for child in self.children:
                data = bytes(child) + struct.pack("Q", offset)  # FIXME
                assert len(data) == indexSlotSize
                output.write(data)
                offset += childNodeSize
            output.write(bytes((blockSize - countOne) * indexSlotSize))
        else:
            for child in self.children:
                offset = child.rWriteIndexLevel(
                    blockSize, childNodeSize, curLevel + 1, destLevel, offset, stream
                )
        return offset


class RangeTree:
    def __init__(self):
        self.root = None
        self.n = 0
        self.freeList = []


class Range:
    __slots__ = ("next", "start", "end", "val")

    def __init__(self, start, end, val):
        self.start = start
        self.end = end
        self.val = val


class rbTreeNode:
    __slots__ = ("left", "right", "color", "item")


class bbNamedFileChunk:
    __slots__ = ("name", "offset", "size")

    def __lt__(self, other):
        return self.name < other.name


def bbiChromUsageFromBedFile(alignments, targets, eim):
    aveSize = 0
    chromId = 0
    totalBases = 0
    bedCount = 0
    usageName = ""
    usage = []
    minDiff = sys.maxsize
    for alignment in alignments:
        chrom = alignment.target.id
        start = alignment.coordinates[0, 0]
        end = alignment.coordinates[0, -1]
        for element in eim:
            element.updateMaxFieldSize(alignment)
        if start > end:
            raise ValueError(
                f"end ({end}) before start ({start}) in alignment [{bedCount}]"
            )
        bedCount += 1
        totalBases += end - start
        if usageName != chrom:
            if usageName > chrom:
                raise ValueError(
                    f"alignments are not sorted by target name at alignment [{counter}]"
                )
            if usageName:
                usage.append(bbiChromUsage(usageName, itemCount, chromId, chromSize))
                chromId += 1
            for target in targets:
                if target.id == chrom:
                    break
            else:
                raise ValueError(
                    f"failed to find target '{target.name}' in target list at alignment [{counter}]"
                )
            usageName = chrom
            chromSize = len(target)
            itemCount = 1
            lastStart = -1
        else:
            itemCount += 1
        if end > chromSize:
            raise ValueError(
                f"end coordinate {end} bigger than {chrom} size of {chromSize} at alignment [{counter}]'"
            )
        if lastStart >= 0:
            diff = start - lastStart
            if diff < minDiff:
                if diff < 0:
                    raise ValueError(
                        f"alignments are not sorted at alignment [{counter}]"
                    )
                minDiff = diff
        lastStart = start
    if usageName:
        usage.append(bbiChromUsage(usageName, itemCount, chromId, chromSize))
    if bedCount > 0:
        aveSize = totalBases / bedCount
    return usage, minDiff, aveSize, bedCount


def bbiWriteChromInfo(usageList, blockSize, output):
    chromCount = len(usageList)
    chromInfoList = []
    maxChromNameSize = 0
    for chromId, usage in enumerate(usageList):
        name = usage.name
        maxChromNameSize = max(maxChromNameSize, len(name))
        itemCount = usage.itemCount
        assert usage.id == chromId
        chromSize = usage.size
        chromInfoList.append([name, chromId, chromSize])
    chromInfoList.sort()
    chromBlockSize = min(blockSize, chromCount)

    itemArray = chromInfoList
    itemCount = chromCount
    blockSize = chromBlockSize
    keySize = maxChromNameSize
    valSize = 8

    byteorder = sys.byteorder
    signature = 0x78CA8C91
    formatter = struct.Struct("=IIIIQxxxxxxxx")
    data = formatter.pack(signature, blockSize, keySize, valSize, itemCount)
    output.write(data)
    indexOffset = output.tell()
    levels = 1
    while itemCount > blockSize:
        itemCount = (itemCount + blockSize - 1) // blockSize
        levels += 1
    itemCount = chromCount
    level = levels - 1
    while level > 0:
        slotSizePer = blockSize**level
        nodeSizePer = slotSizePer * blockSize
        nodeCount = (itemCount + nodeSizePer - 1) // nodeSizePer
        bytesInIndexBlock = bptBlockHeaderSize + blockSize * (keySize + 8)
        bytesInLeafBlock = bptBlockHeaderSize + blockSize * (keySize + valSize)
        if level == 1:
            bytesInNextLevelBlock = bytesInLeafBlock
        else:
            bytesInNextLevelBlock = bytesInIndexBlock
        levelSize = nodeCount * bytesInIndexBlock
        endLevel = indexOffset + levelSize
        nextChild = endLevel
        isLeaf = False
        for i in range(0, itemCount, nodeSizePer):
            countOne = min((itemCount - i + slotSizePer - 1) // slotSizePer, blockSize)
            output.write(struct.pack("?xH", isLeaf, countOne))
            slotsUsed = 0
            endIx = min(i + nodeSizePer, itemCount)
            for item in itemArray[i:endIx:slotSizePer]:
                s = item[0].encode().ljust(keySize, b"0")
                data = struct.pack(f"={len(s)}sQ", s, nextChild)
                output.write(data)
                nextChild += bytesInNextLevelBlock
                slotsUsed += 1
            assert slotsUsed == countOne
            slotSize = keySize + 8
            for j in range(countOne, blockSize):
                output.write(bytes(slotSize))
        endLevelOffset = endLevel
        indexOffset = output.tell()
        assert endLevelOffset == indexOffset
        level -= 1
    isLeaf = True
    countLeft = itemCount
    i = 0
    while i < itemCount:
        if countLeft > blockSize:
            countOne = blockSize
        else:
            countOne = countLeft
        output.write(isLeaf.to_bytes(1, byteorder))
        output.write(bytes(1))
        output.write(countOne.to_bytes(2, byteorder))
        for j in range(countOne):
            assert i + j < itemCount
            item = itemArray[i + j]
            output.write(item[0].encode().ljust(keySize, b"\0"))
            output.write(item[1].to_bytes(4, byteorder))
            output.write(item[2].to_bytes(4, byteorder))
        slotSize = keySize + valSize
        for j in range(countOne, blockSize):
            output.write(bytes(slotSize))
        countLeft -= countOne
        i += countOne


def bbiCalcResScalesAndSizes(aveSize):
    resScales = numpy.zeros(bbiMaxZoomLevels, int)
    resSizes = array.array("i", [0] * bbiMaxZoomLevels)
    resTryCount = bbiMaxZoomLevels
    resIncrement = bbiResIncrement
    minZoom = 10
    res = max(int(aveSize), minZoom)
    for resTry in range(resTryCount):
        resScales[resTry] = res
        if res > sys.maxsize // bbiResIncrement:
            resTryCount = resTry + 1
            break
        res *= resIncrement
    return resTryCount, resScales, resSizes


def loadAndValidateBed(alignment, bedFieldCount, fieldCount, declaration):
    row = []
    BB_MAX_CHROM_STRING = 255
    chrom = alignment.target.id
    if len(chrom) >= BB_MAX_CHROM_STRING:
        raise ValueError(
            f"alignment target name '{chrom}' is too long (must not exceed {BB_MAX_CHROM_STRING-1} characters)"
        )
    if len(chrom) < 1:
        raise ValueError("alignment target name cannot be blank or empty")
    chromStart = alignment.coordinates[0, 0]
    chromEnd = alignment.coordinates[0, -1]
    if chromEnd < chromStart:
        raise ValueError(f"chromStart after chromEnd ({chromEnd} > {chromStart})")
    bed = (chrom, chromStart, chromEnd, row)
    if bedFieldCount > 3:
        name = alignment.query.id
        if name == "":
            name = "."
        elif len(name) > 255:
            raise ValueError(
                f"alignment query name '{name}' is too long (must not exceed 255 characters"
            )
        row.append(name)
    if bedFieldCount > 4:
        try:
            score = alignment.score
        except AttributeError:
            score = "."
        else:
            if score < 0 or score > 1000:
                raise ValueError(f"score ({score}) must be between 0 and 1000")
            score = str(score)
        row.append(score)
    if bedFieldCount > 5:
        if alignment.coordinates[1, 0] <= alignment.coordinates[1, -1]:
            strand = "+"
        else:
            strand = "-"
        row.append(strand)
    if bedFieldCount > 6:
        try:
            thickStart = alignment.thickStart
        except AttributeError:
            thickStart = chromStart
        row.append(str(thickStart))
    if bedFieldCount > 7:
        try:
            thickEnd = alignment.thickEnd
        except AttributeError:
            thickEnd = chromEnd
        else:
            if thickEnd < thickStart:
                raise ValueError(
                    f"thickStart ({thickStart}) after thickEnd ({thickEnd})"
                )
            if thickStart != 0 and (thickStart < chromStart or thickStart > chromEnd):
                raise ValueError(
                    f"thickStart out of range for {name}:{chromStart}-{chromEnd}, thick:{thickStart}-{thickEnd}"
                )
            if thickEnd != 0 and (thickEnd < chromStart or thickEnd > chromEnd):
                raise ValueError(
                    f"thickEnd out of range for {name}:{chromStart}-{chromEnd}, thick:{thickStart}-{thickEnd}"
                )
        row.append(str(thickEnd))
    if bedFieldCount > 8:
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
    if bedFieldCount > 9:
        steps = numpy.diff(alignment.coordinates)
        aligned = sum(steps != 0, 0) == 2
        blockSizes = steps.max(0)[aligned]
        blockCount = len(blockSizes)
        row.append(str(blockCount))
    if bedFieldCount > 10:
        row.append(",".join(str(blockSize) for blockSize in blockSizes) + ",")
    if bedFieldCount > 11:
        chromStarts = alignment.coordinates[0, :-1][aligned] - chromStart
        row.append(",".join(str(chromStart) for chromStart in chromStarts) + ",")
    if bedFieldCount > 12:
        if bedFieldCount != 15:
            raise ValueError(
                f"Unexpected value {bedFieldCount} for bedFieldCount in loadAndValidateBed"
            )
        expIds = alignment.annotations["expIds"]
        expScores = alignment.annotations["expScores"]
        expCount = len(expIds)
        assert expCount == len(expScores)
        row.append(str(expCount))
        row.append(",".join(expIds))
        row.append(",".join(str(expScore) for expScore in expScores))
    # parse the declaration
    for i in range(bedFieldCount, fieldCount):
        name = declaration[i].name
        value = alignment.annotations[name]
        if isinstance(value, str):
            row.append(value)
        elif isinstance(value, (int, float)):
            row.append(str(value))
        else:
            row.append(",".join(map(str, value)))
    return bed


def writeBlocks(
    alignments,
    declaration,
    itemsPerSlot,
    doCompress,
    output,
    resTryCount,
    resScales,
    resSizes,
    eim,
    bedCount,
    fieldCount,
    bedN,
):
    maxBlockSize = 0
    lastField = fieldCount - 1
    itemIx = 0
    blockStartOffset = 0
    startPos = 0
    endPos = 0
    chromId = -1
    resEnds = numpy.zeros(resTryCount, int)
    atEnd = False
    start = 0
    end = 0
    sectionStartIx = 0
    sectionEndIx = 0
    bounds = []
    currentChrom = None
    stream = BytesIO()

    def write_data(sectionStartIx):
        data = stream.getvalue()
        size = len(data)
        if doCompress:
            data = zlib.compress(data)
        output.write(data)
        if eim:
            blockEndOffset = output.tell()
            blockSize = blockEndOffset - blockStartOffset
            for element in eim:
                element.addOffsetSize(
                    blockStartOffset,
                    blockSize,
                    sectionStartIx,
                    sectionEndIx,
                )
            sectionStartIx = sectionEndIx
        bounds.append(bbiBoundsArray(blockStartOffset, chromId, startPos, endPos))
        return sectionStartIx, size

    alignments.rewind()
    for alignment in alignments:
        bed = loadAndValidateBed(alignment, bedN, fieldCount, declaration)
        chrom, start, end, row = bed
        if currentChrom is not None:
            if chrom != currentChrom or itemIx >= itemsPerSlot:
                sectionStartIx, size = write_data(sectionStartIx)
                if size > maxBlockSize:
                    maxBlockSize = size
                itemIx = 0
                stream = BytesIO()
        if chrom != currentChrom:
            currentChrom = chrom
            resEnds[:] = 0
            chromId += 1

        if itemIx == 0:
            blockStartOffset = output.tell()
            startPos = start
            endPos = end
        else:
            endPos = max(endPos, end)

        if eim:
            for element in eim:
                element.addKeysFromRow(alignment, sectionEndIx)
            sectionEndIx += 1

        data = struct.pack("III", chromId, start, end)
        stream.write(data)
        if fieldCount > 3:
            rest = "\t".join(row)
            stream.write(rest.encode())
        stream.write(b"\0")

        itemIx += 1
        for resTry in range(resTryCount):
            resEnd = resEnds[resTry]
            if start >= resEnd:
                resSizes[resTry] += 1
                resEnd = start + resScales[resTry]
            while end > resEnd:
                resSizes[resTry] += 1
                resEnd += resScales[resTry]
            resEnds[resTry] = resEnd

    sectionStartIx, size = write_data(sectionStartIx)
    if size > maxBlockSize:
        maxBlockSize = size

    return maxBlockSize, bounds


def rTreeFromChromRangeArray(blockSize, items, endFileOffset):
    itemCount = len(items)
    if itemCount == 0:
        return
    children = []
    nextOffset = items[0].offset
    oneSize = 0
    i = 0
    while i < itemCount:
        child = rTree()
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
        for j in range(1, oneSize):
            item = items[i + j]
            if item.chromId < child.startChromId:
                child.startChromId = item.chromId
                child.startBase = item.start
            elif item.chromId == child.startChromId and item.start < child.startBase:
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
                parent = child.clone()
                parent.children = [child]
                child.parent = parent
                parents.append(parent)
            else:
                slotsUsed += 1
                parent.children.append(child)
                child.parent = parent
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
        levelCount += 1
        if len(parents) == 1:
            break
        children = parents

    return parent, levelCount


def calcLevelSizes(tree, levelSizes, level, maxLevel):
    levelSizes[level] += 1
    if level < maxLevel:
        level += 1
        for el in tree.children:
            calcLevelSizes(el, levelSizes, level, maxLevel)


def rWriteLeaves(itemsPerSlot, lNodeSize, tree, curLevel, leafLevel, output):
    byteorder = sys.byteorder
    if curLevel == leafLevel:
        reserved = 0
        isLeaf = True
        countOne = len(tree.children)
        output.write(struct.pack("?xH", isLeaf, countOne))
        for el in tree.children:
            output.write(bytes(el))  # FIXME
            data = struct.pack(
                "QQ", el.startFileOffset, el.endFileOffset - el.startFileOffset
            )
            output.write(data)
        output.write(bytes((itemsPerSlot - countOne) * indexSlotSize))
    else:
        for el in tree.children:
            rWriteLeaves(itemsPerSlot, lNodeSize, el, curLevel + 1, leafLevel, output)


def writeTreeToOpenFile(tree, blockSize, levelCount, output):
    level = 0
    levelSizes = numpy.zeros(levelCount, int)
    calcLevelSizes(tree, levelSizes, level, levelCount - 1)
    iNodeSize = nodeHeaderSize + indexSlotSize * blockSize
    lNodeSize = nodeHeaderSize + leafSlotSize * blockSize
    levelOffsets = numpy.zeros(levelCount, numpy.int64) + output.tell()
    levelOffsets[1:] += numpy.cumsum(levelSizes[:-1]) * iNodeSize
    finalLevel = levelCount - 3
    for i in range(finalLevel + 1):
        if i == finalLevel:
            childNodeSize = lNodeSize
        else:
            childNodeSize = iNodeSize
        tree.rWriteIndexLevel(
            blockSize, childNodeSize, 0, i, levelOffsets[i + 1], output
        )
        if output.tell() != levelOffsets[i + 1]:
            raise RuntimeError(
                "Internal error: offset mismatch (%d vs %d)"
                % (output.tell(), levelOffsets[i + 1])
            )

    leafLevel = levelCount - 2
    rWriteLeaves(blockSize, lNodeSize, tree, 0, leafLevel, output)


def cirTreeFileBulkIndexToOpenFile(
    itemArray, blockSize, itemsPerSlot, endFileOffset, output
):
    tree, levelCount = rTreeFromChromRangeArray(blockSize, itemArray, endFileOffset)
    signature = 0x2468ACE0
    data = struct.pack("=IIQ", signature, blockSize, len(itemArray))
    output.write(data)
    output.write(bytes(tree))  # FIXME
    data = struct.pack("QIxxxx", endFileOffset, itemsPerSlot)
    output.write(data)
    if tree is not None:
        writeTreeToOpenFile(tree, blockSize, levelCount, output)


def rbTreeFind(tree, start, end):
    p = tree.root
    while p is not None:
        if end <= p.item.start:
            p = p.left
        elif p.item.end <= start:
            p = p.right
        else:
            return p.item


def rangeCmp(a, b):
    if a.end <= b.start:
        return -1
    if b.end <= a.start:
        return +1
    return 0


def rTreeTraverseRangeWithContext(n, minIt, maxIt, rangeList):
    # n: struct rbTreeNode*
    if n is not None:
        minCmp = rangeCmp(n.item, minIt)
        maxCmp = rangeCmp(n.item, maxIt)
        if minCmp >= 0:
            rTreeTraverseRangeWithContext(n.left, minIt, maxIt, rangeList)
        if minCmp >= 0 and maxCmp <= 0:
            rangeList.append(n.item)
        if maxCmp <= 0:
            rTreeTraverseRangeWithContext(n.right, minIt, maxIt, rangeList)


def rangeTreeAllOverlapping(tree, start, end):
    rangeList = []
    tempR = Range(start, end, None)
    rTreeTraverseRangeWithContext(tree.root, tempR, tempR, rangeList)
    return rangeList


def rTreeTraverseWithContext(n, rangeList):
    if n is not None:
        rTreeTraverseWithContext(n.left, rangeList)
        rangeList.append(n.item)
        rTreeTraverseWithContext(n.right, rangeList)


def rangeTreeList(tree):
    rangeList = []
    rTreeTraverseWithContext(tree.root, rangeList)
    return rangeList


def restructure(t, tos, x, y, z):
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
        parent = t.stack[tos - 1]
        if x is parent.left:
            parent.left = midNode
        else:
            parent.right = midNode
    else:
        t.root = midNode
    return midNode


def rbTreeAdd(tree, item):
    tree.stack = []
    try:
        x = tree.freeList.pop()
    except IndexError:
        x = rbTreeNode()
    else:
        tree.freeList = x.right
    x.left = None
    x.right = None
    x.item = item
    x.color = None
    p = tree.root
    if p is not None:
        while True:
            tree.stack.append(p)
            if item.end <= p.item.start:
                p = p.left
                if p is None:
                    p = tree.stack.pop()
                    p.left = x
                    break
            elif p.item.end <= item.start:
                p = p.right
                if p is None:
                    p = tree.stack.pop()
                    p.right = x
                    break
            else:
                return p.item
        col = "red"
    else:
        tree.root = x
        col = "black"

    x.color = col

    tree.n += 1
    if len(tree.stack) > 0:
        while p.color == "red":
            m = tree.stack.pop()
            if p == m.left:
                q = m.right
            else:
                q = m.left
            if q is None or q.color == "black":
                m = restructure(tree, len(tree.stack), m, p, x)
                m.color = "black"
                m.left.color = m.right.color = "red"
                break
            p.color = "black"
            q.color = "black"
            if len(tree.stack) == 0:
                break
            m.color = "red"
            x = m
            p = tree.stack.pop()


def rangeTreeAddToCoverageDepth(tree, start, end):
    existing = rbTreeFind(tree, start, end)
    if existing is None:
        r = Range(start, end, val=1)
        rbTreeAdd(tree, r)
    else:
        if existing.start <= start and existing.end >= end:
            if existing.start < start:
                r = Range(existing.start, start, existing.val)
                existing.start = start
                rbTreeAdd(tree, r)
            if existing.end > end:
                r = Range(end, existing.end, existing.val)
                existing.end = end
                rbTreeAdd(tree, r)
            existing.val += 1
        else:
            existingList = rangeTreeAllOverlapping(tree, start, end)
            s = start
            for existing in existingList:
                if s < existing.start:
                    r = Range(s, existing.start, 1)
                    s = existing.start
                    rbTreeAdd(tree, r)
                elif s > existing.start:
                    r = Range(existing.start, s, existing.val)
                    existing.start = s
                    rbTreeAdd(tree, r)
                if existing.start < end and existing.end > end:
                    r = Range(end, existing.end, existing.val)
                    existing.end = end
                    rbTreeAdd(tree, r)
                existing.val += 1
                s = existing.end
            if s < end:
                r = Range(s, end, 1)
                rbTreeAdd(tree, r)


def rangeTreeGenerator(alignments):
    name = None
    alignments.rewind()
    for alignment in alignments:
        chrom = alignment.target.id
        if name != chrom:
            if name is not None:
                yield name, tree
            name = chrom
            tree = RangeTree()
        start = alignment.coordinates[0, 0]
        end = alignment.coordinates[0, -1]
        if start > end:
            start, end = end, start
        rangeTreeAddToCoverageDepth(tree, start, end)
    if name is not None:
        yield name, tree


def bbiOutputOneSummaryFurtherReduce(
    summary, twiceReducedList, doubleReductionSize, boundsArray, stream
):
    offset = stream.output.tell()
    bounds = bbiBoundsArray(offset, summary.chromId, summary.start, summary.end)
    boundsArray.append(bounds)

    stream.write(summary)
    try:
        twiceReduced = twiceReducedList[-1]
    except IndexError:
        pass
    else:
        if (
            twiceReduced.chromId == summary.chromId
            and twiceReduced.start + doubleReductionSize >= summary.end
        ):
            twiceReduced.end = summary.end
            twiceReduced.validCount += summary.validCount
            twiceReduced.values["minVal"] = min(
                twiceReduced.values["minVal"], summary.values["minVal"]
            )
            twiceReduced.values["maxVal"] = max(
                twiceReduced.values["maxVal"], summary.values["maxVal"]
            )
            twiceReduced.values["sumData"] += summary.values["sumData"]
            twiceReduced.values["sumSquares"] += summary.values["sumSquares"]
            return
    twiceReduced = copy.copy(summary)
    twiceReduced.values = summary.values.copy()
    twiceReducedList.append(twiceReduced)


def bedWriteReducedOnceReturnReducedTwice(
    usageList,
    fieldCount,
    alignments,
    initialReduction,
    initialReductionCount,
    zoomIncrement,
    blockSize,
    itemsPerSlot,
    doCompress,
    output,
):
    byteorder = sys.byteorder
    twiceReducedList = []
    doubleReductionSize = initialReduction * zoomIncrement
    usage = usageList
    boundsArray = []

    dataStart = output.tell()
    output.write(initialReductionCount.to_bytes(4, byteorder))

    stream = bbiSumOutStream(itemsPerSlot, output, doCompress)
    totalSum = bbiSummaryElement()
    rangeTrees = rangeTreeGenerator(alignments)
    for usage in usageList:
        summary = None
        name, rangeTree = next(rangeTrees)
        assert name == usage.name
        rangeList = rangeTreeList(rangeTree)
        for range in rangeList:
            val = range.val
            start = range.start
            end = range.end
            size = end - start
            if size == 0:
                size = 1

            totalSum.update(size, val)

            if (
                summary is not None
                and summary.end <= start
                and summary.end < usage.size
            ):
                bbiOutputOneSummaryFurtherReduce(
                    summary, twiceReducedList, doubleReductionSize, boundsArray, stream
                )
                summary = None
            if summary is None:
                summary = bbiSummary(
                    usage.id, start, min(start + initialReduction, usage.size), val
                )
            while end > summary.end:
                overlap = min(end, summary.end) - max(start, summary.start)
                assert overlap > 0
                summary.update(overlap, val)
                bbiOutputOneSummaryFurtherReduce(
                    summary, twiceReducedList, doubleReductionSize, boundsArray, stream
                )
                size -= overlap
                summary.start = start = summary.end
                summary.end = min(start + initialReduction, usage.size)
                summary.values["minVal"] = summary.values["maxVal"] = val
                summary.values["sumData"] = summary.values["sumSquares"] = 0.0
                summary.validCount = 0

            summary.update(size, val)
        if summary is not None:
            bbiOutputOneSummaryFurtherReduce(
                summary, twiceReducedList, doubleReductionSize, boundsArray, stream
            )
    stream.flush()

    assert len(boundsArray) == initialReductionCount
    indexOffset = output.tell()
    cirTreeFileBulkIndexToOpenFile(
        boundsArray, blockSize, itemsPerSlot, indexOffset, output
    )

    return twiceReducedList, dataStart, indexOffset, totalSum


def bbiSummarySimpleReduce(summaries, reduction):
    newSummaries = []
    newSummary = None
    for summary in summaries:
        if (
            newSummary is None
            or newSummary.chromId != summary.chromId
            or summary.end > newSummary.start + reduction
        ):
            newSummary = copy.copy(summary)
            newSummary.values = summary.values.copy()
            newSummaries.append(newSummary)
        else:
            assert newSummary.end < summary.end
            newSummary.end = summary.end
            newSummary.validCount += summary.validCount
            if newSummary.values["minVal"] > summary.values["minVal"]:
                newSummary.values["minVal"] = summary.values["minVal"]
            if newSummary.values["maxVal"] < summary.values["maxVal"]:
                newSummary.values["maxVal"] = summary.values["maxVal"]
            newSummary.values["sumData"] += summary.values["sumData"]
            newSummary.values["sumSquares"] += summary.values["sumSquares"]
    return newSummaries


def bbiWriteSummary(summaryList, itemsPerSlot, doCompress, output):
    # See bbiWriteSummaryAndIndexUnc, bbiWriteSummaryAndIndexComp in bbiWrite.c
    count = len(summaryList)
    data = struct.pack("I", count)
    output.write(data)
    if doCompress:
        for start in range(0, count, itemsPerSlot):
            buffer = BytesIO()
            filePos = output.tell()
            for summary in summaryList[start : start + itemsPerSlot]:
                buffer.write(bytes(summary))
                summary.offset = filePos
            data = zlib.compress(buffer.getvalue())
            output.write(data)
    else:
        for summary in summaryList:
            summary.offset = output.tell()
            output.write(bytes(summary))


def bbiWriteZoomLevels(
    alignments,  # struct lineFile *lf,    /* Input file. */
    output,  # FILE *f,                /* Output. */
    blockSize,  # int blockSize,          /* Size of index block */
    itemsPerSlot,  # int itemsPerSlot,       /* Number of data points bundled at lowest level. */
    fieldCount,  # int fieldCount,         /* Number of fields in bed (4 for bedGraph) */
    doCompress,  # boolean doCompress,     /* Do we compress.  Answer really should be yes! */
    dataSize,  # bits64 dataSize,        /* Size of data on disk (after compression if any). */
    usageList,  # struct bbiChromUsage *usageList, /* Result from bbiChromUsageFromBedFile */
    resTryCount,
    resScales,
    resSizes,  # int resTryCount, int resScales[], int resSizes[],   /* How much to zoom at each level */
    zoomAmounts,  # bits32 zoomAmounts[bbiMaxZoomLevels],      /* Fills in amount zoomed at each level. */
    zoomDataOffsets,  # bits64 zoomDataOffsets[bbiMaxZoomLevels],  /* Fills in where data starts for each zoom level. */
    zoomIndexOffsets,  # bits64 zoomIndexOffsets[bbiMaxZoomLevels], /* Fills in where index starts for each level. */
):
    assert resTryCount > 0
    maxReducedSize = dataSize / 2
    initialReduction = 0
    initialReducedCount = 0

    for resTry in range(resTryCount):
        reducedSize = resSizes[resTry] * bbiSummary.size
        if doCompress:
            reducedSize /= 2
        if reducedSize <= maxReducedSize:
            initialReduction = resScales[resTry]
            initialReducedCount = resSizes[resTry]
            break
    else:
        initialReduction = resScales[0]
        initialReducedCount = resSizes[0]
    zoomIncrement = bbiResIncrement
    (
        rezoomedList,
        zoomDataOffsets[0],
        zoomIndexOffsets[0],
        totalSum,
    ) = bedWriteReducedOnceReturnReducedTwice(
        usageList,
        fieldCount,
        alignments,
        initialReduction,
        initialReducedCount,
        zoomIncrement,
        blockSize,
        itemsPerSlot,
        doCompress,
        output,
    )
    zoomAmounts[0] = initialReduction
    zoomLevels = 1
    zoomCount = initialReducedCount
    reduction = initialReduction * zoomIncrement
    while zoomLevels < bbiMaxZoomLevels:
        rezoomCount = len(rezoomedList)
        if rezoomCount >= zoomCount:
            break
        zoomCount = rezoomCount
        zoomDataOffsets[zoomLevels] = output.tell()
        bbiWriteSummary(rezoomedList, itemsPerSlot, doCompress, output)
        indexOffset = output.tell()
        cirTreeFileBulkIndexToOpenFile(
            rezoomedList, blockSize, itemsPerSlot, indexOffset, output
        )
        zoomIndexOffsets[zoomLevels] = indexOffset
        zoomAmounts[zoomLevels] = reduction
        zoomLevels += 1
        reduction *= zoomIncrement
        rezoomedList = bbiSummarySimpleReduce(rezoomedList, reduction)
    return zoomLevels, totalSum


def bptCountLevels(maxBlockSize, itemCount):
    levels = 1
    while itemCount > maxBlockSize:
        itemCount = (itemCount + maxBlockSize - 1) / maxBlockSize
        levels += 1
    return levels


def writeIndexLevel(
    blockSize, chunkArray, itemCount, indexOffset, level, keySize, valSize, output
):
    # in bPlusTree.c
    byteorder = sys.byteorder
    slotSizePer = pow(blockSize, level)
    nodeSizePer = slotSizePer * blockSize
    nodeCount = (itemCount + nodeSizePer - 1) // nodeSizePer
    bytesInIndexBlock = bptBlockHeaderSize + blockSize * (keySize + 8)
    bytesInLeafBlock = bptBlockHeaderSize + blockSize * (keySize + valSize)
    if level == 1:
        bytesInNextLevelBlock = bytesInLeafBlock
    else:
        bytesInNextLevelBlock = bytesInIndexBlock
    levelSize = nodeCount * bytesInIndexBlock
    endLevel = indexOffset + levelSize
    nextChild = endLevel
    for i in range(0, itemCount, nodeSizePer):
        countOne = min((itemCount - i + slotSizePer - 1) // slotSizePer, blockSize)
        output.write(b"\0")  # isLeaf = False
        output.write(b"\0")  # reserved
        output.write(countOne.to_bytes(2, byteorder))
        slotsUsed = 0
        endIx = min(i + nodeSizePer, itemCount)
        for chunk in chunkArray[i:endIx:slotSizePer]:  # bbNamedFileChunk
            output.write(chunk.name.ljust(keySize, b"\0"))
            output.write(nextChild.to_bytes(8, byteorder))
            nextChild += bytesInNextLevelBlock
            slotsUsed += 1
        assert slotsUsed == countOne
        output.write(bytes((keySize + 8) * (blockSize - countOne)))
    return endLevel


def writeLeafLevel(blockSize, chunkArray, itemCount, keySize, valSize, output):
    byteorder = sys.byteorder
    countLeft = itemCount
    i = 0
    while i < itemCount:
        countOne = min(blockSize, countLeft)
        output.write(b"\1")  # isLeaf = True
        output.write(b"\0")  # reserved
        output.write(countOne.to_bytes(2, byteorder))
        for j in range(countOne):
            assert i + j < itemCount
            chunk = chunkArray[i + j]
            output.write(chunk.name.ljust(keySize, b"\0"))
            output.write(chunk.offset.to_bytes(8, byteorder))
            output.write(chunk.size.to_bytes(8, byteorder))
        output.write(bytes((keySize + valSize) * (blockSize - countOne)))
        countLeft -= countOne
        i += countOne


def bptFileBulkIndexToOpenFile(chunkArray, itemCount, blockSize, keySize, output):
    signature = 0x78CA8C91
    valSize = 16
    formatter = struct.Struct("=IIIIQxxxxxxxx")
    data = formatter.pack(signature, blockSize, keySize, valSize, itemCount)
    output.write(data)
    indexOffset = output.tell()
    levels = bptCountLevels(blockSize, itemCount)
    for i in range(levels - 1, 0, -1):
        endLevelOffset = writeIndexLevel(
            blockSize, chunkArray, itemCount, indexOffset, i, keySize, valSize, output
        )  # in bPlusTree.c
        indexOffset = output.tell()
        assert endLevelOffset == indexOffset
    writeLeafLevel(blockSize, chunkArray, itemCount, keySize, valSize, output)
