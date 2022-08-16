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


BptFile = namedtuple(
    "BptFile",
    (
        "blockSize",  # 4 bytes
        "keySize",  # 4 bytes
        "valSize",  # 4 bytes
        "itemCount",  # 8 bytes
        "rootOffset",  # 8 bytes
    ),
)

BbiFile = namedtuple(
    "BbiFile",
    (
        "chromBpt",  # BptFile
        "version",  # 2 bytes
        "zoomLevels",  # 2 bytes
        "chromTreeOffset",  # 8 bytes
        "unzoomedDataOffset",  # 8 bytes
        "unzoomedIndexOffset",  # 8 bytes
        "fieldCount",  # 2 bytes
        "definedFieldCount",  # 2 bytes
        "asOffset",  # 8 bytes
        "totalSummaryOffset",  # 8 bytes
        "uncompressBufSize",  # 4 bytes
        "extensionOffset",  # 8 bytes
        "unzoomedCir",  # CirTreeFile
        "extensionSize",  # 2 bytes
        "extraIndexCount",  # 2 bytes
        "extraIndexListOffset",  # 8 bytes
    ),
)

CirTreeFile = namedtuple(
    "CirTreeFile",
    (
        "rootOffset",  # 8 bytes
        "blockSize",  # 4 bytes
        "itemCount",  # 8 bytes
        "startChromIx",  # 4 bytes
        "startBase",  # 4 bytes
        "endChromIx",  # 4 bytes
        "endBase",  # 4 bytes
        "fileSize",  # 8 bytes
        "itemsPerSlot",  # 4 bytes
    ),
)

BbiChromInfo = namedtuple(
    "BbiChromInfo",
    (
        "name",  # str
        "id",  # 4 bytes
        "size",  # 4 bytes
    ),
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
        signature = 0x8789F2EB
        magic = stream.read(4)
        for byteorder, byteorder_character in zip(("little", "big"), "<>"):
            if int.from_bytes(magic, byteorder=byteorder) == signature:
                self._byteorder = byteorder
                self._byteorder_character = byteorder_character
                break
        else:
            raise ValueError("not a bigBed file")
        (
            version,
            zoomLevels,
            chromTreeOffset,
            unzoomedDataOffset,
            unzoomedIndexOffset,
            fieldCount,
            definedFieldCount,
            asOffset,
            totalSummaryOffset,
            uncompressBufSize,
            extensionOffset,
        ) = struct.unpack(byteorder_character + "hhqqqhhqqiq", stream.read(60))

        if extensionOffset != 0:
            stream.seek(extensionOffset)
            extensionSize, extraIndexCount, extraIndexListOffset = struct.unpack(
                byteorder_character + "iii", stream.read(12)
            )

        signature = 0x78CA8C91
        stream.seek(chromTreeOffset)
        magic = int.from_bytes(stream.read(4), byteorder=byteorder)
        assert magic == signature
        blockSize, keySize, valSize, itemCount = struct.unpack(
            byteorder_character + "iiiqxxxxxxxx", stream.read(28)
        )
        rootOffset = stream.tell()
        chromBpt = BptFile(blockSize, keySize, valSize, itemCount, rootOffset)

        signature = 0x2468ACE0
        stream.seek(unzoomedIndexOffset)
        magic = int.from_bytes(stream.read(4), byteorder=byteorder)
        assert magic == signature
        (
            blockSize,
            itemCount,
            startChromIx,
            startBase,
            endChromIx,
            endBase,
            fileSize,
            itemsPerSlot,
        ) = struct.unpack(byteorder_character + "iqiiiiqixxxx", stream.read(44))

        rootOffset = stream.tell()
        unzoomedCir = CirTreeFile(
            rootOffset,
            blockSize,
            itemCount,
            startChromIx,
            startBase,
            endChromIx,
            endBase,
            fileSize,
            itemsPerSlot,
        )

        bbi = BbiFile(
            chromBpt,
            version,
            zoomLevels,
            chromTreeOffset,
            unzoomedDataOffset,
            unzoomedIndexOffset,
            fieldCount,
            definedFieldCount,
            asOffset,
            totalSummaryOffset,
            uncompressBufSize,
            extensionOffset,
            unzoomedCir,
            extensionSize,
            extraIndexCount,
            extraIndexListOffset,
        )

        chromList = []
        bpt = bbi.chromBpt
        self._rTraverse(stream, bpt, bpt.rootOffset, chromList)
        targets = {}

        for item in sorted(chromList, key=lambda item: item.id):
            name = item.name.decode()
            length = item.size
            sequence = Seq(None, length=length)
            target = SeqRecord(sequence, id=name)
            targets[name] = target

        self.targets = targets

    def _read_next_alignment(self, stream):
        try:
            line = next(stream)
        except StopIteration:
            return None
        words = line.split()
        bedN = len(words)
        if bedN < 3 or bedN > 12:
            raise ValueError("expected between 3 and 12 columns, found %d" % bedN)
        chrom = words[0]
        chromStart = int(words[1])
        chromEnd = int(words[2])
        if bedN > 3:
            name = words[3]
        else:
            name = None
        if bedN > 5:
            strand = words[5]
        else:
            strand = "+"
        if bedN > 9:
            blockCount = int(words[9])
            blockSizes = [
                int(blockSize) for blockSize in words[10].rstrip(",").split(",")
            ]
            blockStarts = [
                int(blockStart) for blockStart in words[11].rstrip(",").split(",")
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
        if bedN <= 4:
            return alignment
        score = words[4]
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
        alignment.thickStart = int(words[6])
        if bedN <= 7:
            return alignment
        alignment.thickEnd = int(words[7])
        if bedN <= 8:
            return alignment
        alignment.itemRgb = words[8]
        return alignment

    def _rTraverse(self, stream, bpt, blockStart, chromList):
        byteorder = self._byteorder
        byteorder_character = self._byteorder_character
        stream.seek(blockStart)
        isLeaf, childCount = struct.unpack(byteorder_character + "?xh", stream.read(4))
        keySize = bpt.keySize
        valSize = bpt.valSize
        fmt = byteorder_character + "II"
        if isLeaf:
            for i in range(childCount):
                key = stream.read(keySize)
                val = stream.read(valSize)
                name = key.split(b"\x00", 1).pop(0)
                chromId, chromSize = struct.unpack(fmt, val)
                info = BbiChromInfo(name, chromId, chromSize)
                chromList.append(info)
        else:
            fileOffsets = numpy.empty(childCount, "i8")
            for i in range(childCount):
                key = stream.read(keySize)
                fileOffsets[i] = int.from_bytes(stream.read(8), byteorder)
            for i in range(childCount):
                self._rTraverse(stream, bpt, fileOffsets[i], chromList)


code = """\
clChrom = sys.argv[1]
clStart = int(sys.argv[2])
clEnd = int(sys.argv[3])



Block = namedtuple("Block", ("offset", "size"))



def rFind(bpt, blockStart, name):
    stream.seek(blockStart)
    isLeaf, childCount = struct.unpack("<?xh", stream.read(4))
    if isLeaf:
        for i in range(childCount):
            key = stream.read(bpt.keySize)
            val = stream.read(bpt.valSize)
            if name == key.split(b"\x00", 1).pop(0):
                return val
    else:
        stream.read(bpt.keySize)
        fileOffset = int.from_bytes(stream.read(8), byteorder=byteorder)
        for i in range(1, childCount):
            key = stream.read(bpt.keySize)
            if name < key.split(b"\x00", 1).pop(0):
                break
            fileOffset = int.from_bytes(stream.read(8), byteorder=byteorder)
        return rFind(bpt, fileOffset, name)


def cirTreeOverlaps(qChrom, qStart, qEnd, rStartChrom, rStartBase, rEndChrom, rEndBase):
    return (qChrom, qStart) < (rEndChrom, rEndBase) and (rStartChrom, rStartBase) < (
        qChrom,
        qEnd,
    )


def rFindOverlappingBlocks(crt, level, indexFileOffset, chromIx, start, end, blockList):
    stream.seek(indexFileOffset)
    isLeaf, childCount = struct.unpack("<?xh", stream.read(4))
    if isLeaf:
        for i in range(childCount):
            startChromIx, startBase, endChromIx, endBase, offset, size = struct.unpack(
                "<iiiiqq", stream.read(32)
            )
            if cirTreeOverlaps(
                chromIx, start, end, startChromIx, startBase, endChromIx, endBase
            ):
                block = Block(offset, size)
                blockList.append(block)
    else:
        startChromIx = numpy.empty(childCount, "i4")
        startBase = numpy.empty(childCount, "i4")
        endChromIx = numpy.empty(childCount, "i4")
        endBase = numpy.empty(childCount, "i4")
        offset = numpy.empty(childCount, "i8")
        for i in range(childCount):
            (
                startChromIx[i],
                startBase[i],
                endChromIx[i],
                endBase[i],
                offset[i],
            ) = struct.unpack("<iiiiq", stream.read(24))
        for i in range(childCount):
            if cirTreeOverlaps(
                chromIx,
                start,
                end,
                startChromIx[i],
                startBase[i],
                endChromIx[i],
                endBase[i],
            ):
                rFindOverlappingBlocks(
                    crt, level + 1, offset[i], chromIx, start, end, blockList
                )

def fileOffsetSizeFindGap(blocks):
    beforeGap = blocks[0]
    for afterGap in blocks[1:]:
        if afterGap.offset != beforeGap.offset + beforeGap.size:
            break
        beforeGap = afterGap
    else:
        afterGap = None
    return beforeGap, afterGap

def bigBedIntervalQuery(bbi, chromName, start, end):
    intervals = []
    paddedStart = max(start - 1, 0)
    paddedEnd = end + 1

    bpt = bbi.chromBpt
    if len(chromName) > bpt.keySize:
        return
    if bpt.valSize != 8:
        raise ValueError(
            "Value size mismatch in bptFileFind (valSize=%; expected 8)" % bpt.valSize
        )
    value = rFind(bpt, bpt.rootOffset, chromName)
    if value is None:
        return intervals
    if byteorder == "little":
        fmt = "<II"
    elif byteorder == "big":
        fmt = ">II"
    chromId, chromSize = struct.unpack(fmt, value)
    crt = bbi.unzoomedCir
    blocks = []
    level = 0
    rFindOverlappingBlocks(
        crt, level, crt.rootOffset, chromId, paddedStart, paddedEnd, blocks
    )

    while blocks:
        beforeGap, afterGap = fileOffsetSizeFindGap(blocks)
        mergedOffset = blocks[0].offset
        mergedSize = beforeGap.offset + beforeGap.size - mergedOffset
        stream.seek(mergedOffset)
        blockBuf = io.BytesIO()
        blockBuf.write(stream.read(mergedSize))
        blockBuf.seek(0)
        while True:
            try:
                block = blocks.pop(0)
            except IndexError:
                assert afterGap is None
                break
            else:
                if block == afterGap:
                    break
            if bbi.uncompressBufSize > 0:
                data = blockBuf.read(block.size)
                data = zlib.decompress(data)
            else:
                data = blockBuf.read()
            while data:
                ci, s, e = struct.unpack("<III", data[:12])
                rest, data = data[12:].split(b'\00', 1)
                if ci == chromId and (
                    (s < end and e > start) or (s == e and (s == end or e == start))
                ):
                    interval = (chromId, s, e, rest)
                    intervals.append(interval)
    return intervals

if clStart > 0:
    start = clStart
if clEnd > 0:
    end = clEnd

intervals = bigBedIntervalQuery(bbi, chromName, start, end)

names = list(targets.keys())
for interval in intervals:
    chromId, start, end, extra = interval
    name = names[chromId]
    print("%s\t%d\t%d\t%s" % (name, start, end, extra.decode()))
"""
