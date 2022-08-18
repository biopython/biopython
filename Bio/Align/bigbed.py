# Copyright 2022 by Michiel de Hoon.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.Align support for BED (Browser Extensible Data) files.

The Browser Extensible Data (BED) format, stores a series of pairwise
alignments in a single file. Typically they are used for transcript to genome
alignments. BED files store the alignment positions and alignment scores, but
not the aligned sequences.

See http://genome.ucsc.edu/FAQ/FAQformat.html#format1

You are expected to use this module via the Bio.Align functions.

Coordinates in the BED format are defined in terms of zero-based start
positions (like Python) and aligning region sizes.

A minimal aligned region of length one and starting at first position in the
source sequence would have ``start == 0`` and ``size == 1``.

As we can see in this example, ``start + size`` will give one more than the
zero-based end position. We can therefore manipulate ``start`` and
``start + size`` as python list slice boundaries.
"""
import numpy
import io
import sys
import zlib
import struct
from collections import namedtuple


from Bio.Align import Alignment
from Bio.Align import interfaces
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import BiopythonExperimentalWarning

import warnings

warnings.warn(
    "Bio.Align.bed is an experimental module which may undergo "
    "significant changes prior to its future official release.",
    BiopythonExperimentalWarning,
)


class AlignmentWriter(interfaces.AlignmentWriter):
    """Alignment file writer for the Browser Extensible Data (BED) file format."""

    def __init__(self, target, bedN=12):
        """Create an AlignmentWriter object.

        Arguments:
         - target    - output stream or file name
         - bedN      - number of columns in the BED file.
                       This must be between 3 and 12; default value is 12.

        """
        if bedN < 3 or bedN > 12:
            raise ValueError("bedN must be between 3 and 12")
        super().__init__(target, mode="w")
        self.bedN = bedN

    def format_alignment(self, alignment):
        """Return a string with one alignment formatted as a BED line."""
        if not isinstance(alignment, Alignment):
            raise TypeError("Expected an Alignment object")
        coordinates = alignment.coordinates
        if not coordinates.size:  # alignment consists of gaps only
            return ""
        bedN = self.bedN
        target, query = alignment.sequences
        try:
            chrom = target.id
        except AttributeError:
            chrom = "target"
        assert coordinates[0, 0] < coordinates[0, -1]
        if coordinates[1, 0] > coordinates[1, -1]:
            # DNA/RNA mapped to reverse strand of DNA/RNA
            strand = "-"
        else:
            # mapped to forward strand
            strand = "+"
        # variable names follow those in the BED file format specification
        blockSizes = []
        blockStarts = []
        tStart, qStart = coordinates[:, 0]
        for tEnd, qEnd in coordinates[:, 1:].transpose():
            if tStart == tEnd:
                qStart = qEnd
            elif qStart == qEnd:
                tStart = tEnd
            else:
                blockSize = tEnd - tStart
                blockStarts.append(tStart)
                blockSizes.append(blockSize)
                tStart = tEnd
                qStart = qEnd
        chromStart = blockStarts[0]  # start of alignment in target
        chromEnd = blockStarts[-1] + blockSize  # end of alignment in target
        fields = [chrom, str(chromStart), str(chromEnd)]
        if bedN == 3:
            return "\t".join(fields) + "\n"
        try:
            name = query.id
        except AttributeError:
            name = "query"
        fields.append(name)
        if bedN == 4:
            return "\t".join(fields) + "\n"
        try:
            score = alignment.score
        except AttributeError:
            score = 0
        fields.append(str(score))
        if bedN == 5:
            return "\t".join(fields) + "\n"
        fields.append(strand)
        if bedN == 6:
            return "\t".join(fields) + "\n"
        try:
            thickStart = alignment.thickStart
        except AttributeError:
            thickStart = chromStart
        fields.append(str(thickStart))
        if bedN == 7:
            return "\t".join(fields) + "\n"
        try:
            thickEnd = alignment.thickEnd
        except AttributeError:
            thickEnd = chromEnd
        fields.append(str(thickEnd))
        if bedN == 8:
            return "\t".join(fields) + "\n"
        try:
            itemRgb = alignment.itemRgb
        except AttributeError:
            itemRgb = "0"
        fields.append(str(itemRgb))
        if bedN == 9:
            return "\t".join(fields) + "\n"
        blockCount = len(blockSizes)
        fields.append(str(blockCount))
        if bedN == 10:
            return "\t".join(fields) + "\n"
        fields.append(",".join(map(str, blockSizes)) + ",")
        if bedN == 11:
            return "\t".join(fields) + "\n"
        blockStarts -= chromStart
        fields.append(",".join(map(str, blockStarts)) + ",")
        return "\t".join(fields) + "\n"


class AlignmentIterator(interfaces.AlignmentIterator):
    """Alignment iterator bigBed files.

    The pairwise alignments stored in the bigBed file are loaded and returned
    incrementally.  Additional alignment information is stored as attributes
    of each alignment.
    """

    def __init__(self, source):
        """Create an AlignmentIterator object.

        Arguments:
         - source   - input data or file name

        """
        super().__init__(source, mode="b", fmt="bigBed")

    def _read_header(self, stream):

        # Supplemental Table 5: Common header
        # magic                  4 bytes
        # version                2 bytes
        # zoomLevels             2 bytes
        # chromosomeTreeOffset   8 bytes
        # fullDataOffset         8 bytes; points to dataCount
        # fullIndexOffset        8 bytes
        # fieldCount             2 bytes
        # definedFieldCount      2 bytes
        # autoSqlOffset          8 bytes
        # totalSummaryOffset     8 bytes
        # uncompressBufSize      4 bytes
        # reserved               8 bytes
        signature = 0x8789F2EB
        magic = stream.read(4)
        for byteorder in ("little", "big"):
            if int.from_bytes(magic, byteorder=byteorder) == signature:
                break
        else:
            raise ValueError("not a bigBed file")

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
        ) = struct.unpack("<hhqqqhhqqixxxxxxxx", stream.read(60))

        stream.seek(fullDataOffset)
        dataCount = int.from_bytes(stream.read(8), byteorder=byteorder)

        if uncompressBufSize > 0:
            self._compressed = True
        else:
            self._compressed = False

        # Supplemental Table 8: Chromosome B+ tree header
        # magic     4 bytes
        # blockSize 4 bytes
        # keySize   4 bytes
        # valSize   4 bytes
        # itemCOunt 8 bytes
        # reserved  8 bytes
        signature = 0x78CA8C91
        stream.seek(chromosomeTreeOffset)
        magic = int.from_bytes(stream.read(4), byteorder=byteorder)
        assert magic == signature
        blockSize, keySize, valSize, itemCount = struct.unpack(
            "<iiiqxxxxxxxx", stream.read(28)
        )
        assert valSize == 8

        while True:
            # Supplemental Table 9: Chromosome B+ tree node
            # isLeaf    1 byte
            # reserved  1 byte
            # count     2 bytes
            isLeaf, count = struct.unpack("<?xh", stream.read(4))
            if isLeaf:
                break
            assert count > 0
            # Supplemental Table 11: Chromosome B+ tree non-leaf item
            # key          keySize bytes
            # childOffset  8 bytes
            key = stream.read(keySize)
            childOffset = int.from_bytes(stream.read(8), byteorder)
            stream.seek(childOffset)

        if byteorder == "little":
            fmt = "<II"
        elif byteorder == "big":
            fmt = ">II"

        targets = {}

        while True:
            for i in range(count):
                # Supplemental Table 10: Chromosome B+ tree leaf item format
                # key        keySize bytes
                # chromId    4 bytes
                # chromSize  4 bytes
                key = stream.read(keySize)
                chromName = key.split(b"\x00", 1).pop(0)
                chromId, chromSize = struct.unpack(fmt, stream.read(8))
                chromName = chromName.decode()
                sequence = Seq(None, length=chromSize)
                target = SeqRecord(sequence, id=chromName)
                targets[chromName] = target
            if chromId + 1 == itemCount:
                break
            # Supplemental Table 9: Chromosome B+ tree node
            # isLeaf    1 byte
            # reserved  1 byte
            # count     2 bytes
            isLeaf, count = struct.unpack("<?xh", stream.read(4))
            assert isLeaf

        # Supplemental Table 14: R tree index header
        # magic            4 bytes
        # blockSize        4 bytes
        # itemCount        8 bytes
        # startChromIx     4 bytes
        # startBase        4 bytes
        # endChromIx       4 bytes
        # endBase          4 bytes
        # endFileOffset    8 bytes
        # itemsPerSlot     4 bytes
        # reserved         4 bytes
        signature = 0x2468ACE0
        stream.seek(fullIndexOffset)
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
        ) = struct.unpack("<iqiiiiqixxxx", stream.read(44))

        self.targets = targets
        self._names = list(targets.keys())
        self._index = 0
        self._length = dataCount
        self._cache = ("", 0)

    def _read_next_alignment(self, stream):
        if self._index == self._length:
            return
        data, count = self._cache
        if not data:
            if not count:
                while True:
                    # Supplemental Table 15: R tree node format
                    # isLeaf     1 byte
                    # reserved   1 byte
                    # count      2 bytes
                    isLeaf, count = struct.unpack("<?xh", stream.read(4))
                    if isLeaf:
                        assert count > 0
                        break
                    else:
                        stream.seek(24 * count, 1)
            # Supplemental Table 16: R tree leaf item format
            # startChromIx   4 bytes
            # startBase      4 bytes
            # endChromIx     4 bytes
            # endBase        4 bytes
            # dataOffset     8 bytes
            # dataSize       8 bytes
            (
                startChromIx,
                startBase,
                endChromIx,
                endBase,
                dataOffset,
                dataSize,
            ) = struct.unpack("<iiiiqq", stream.read(32))
            assert startChromIx == endChromIx
            filepos = stream.tell()
            stream.seek(dataOffset)
            data = stream.read(dataSize)
            stream.seek(filepos)
            if self._compressed > 0:
                data = zlib.decompress(data)
            count -= 1

        # Supplemental Table 12: Binary BED-data format
        # chromId     4 bytes
        # chromStart  4 bytes
        # chromEnd    4 bytes
        # rest        zero-terminated string in tab-separated format
        chromId, chromStart, chromEnd = struct.unpack("<III", data[:12])
        rest, data = data[12:].split(b"\00", 1)
        words = rest.decode().split("\t")
        bedN = len(words) + 3
        if bedN < 3 or bedN > 12:
            raise ValueError("expected between 3 and 12 columns, found %d" % bedN)
        chrom = self._names[chromId]
        if bedN > 3:
            name = words[0]
        else:
            name = None
        if bedN > 5:
            strand = words[2]
        else:
            strand = "+"
        if bedN > 9:
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
        target_record = SeqRecord(None, id=chrom)
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
        self._index += 1
        self._cache = (data, count)
        if bedN <= 4:
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
        if bedN <= 6:
            return alignment
        alignment.thickStart = int(words[3])
        if bedN <= 7:
            return alignment
        alignment.thickEnd = int(words[4])
        if bedN <= 8:
            return alignment
        alignment.itemRgb = words[5]
        return alignment
