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
# "BigWig and BigBed: enabling browsing of large distributed datasets:
# Bioinformatics 26(17): 2204–2207 (2010)
# in particular the tables in the supplemental materials listing the contents
# of a bigBed file byte-by-byte.


import numpy
import zlib
import struct


from Bio.Align import Alignment
from Bio.Align import interfaces
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import BiopythonExperimentalWarning

import warnings

warnings.warn(
    "Bio.Align.bigbed is an experimental module which may undergo "
    "significant changes prior to its future official release.",
    BiopythonExperimentalWarning,
)


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
        for byteorder, byteorder_char in (("little", "<"), ("big", ">")):
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
        ) = struct.unpack(byteorder_char + "hhqqqhhqqixxxxxxxx", stream.read(60))
        if definedFieldCount < 3 or definedFieldCount > 12:
            raise ValueError(
                "expected between 3 and 12 columns, found %d" % definedFieldCount
            )
        self.bedN = definedFieldCount

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
            byteorder_char + "iiiqxxxxxxxx", stream.read(28)
        )
        assert valSize == 8

        while True:
            # Supplemental Table 9: Chromosome B+ tree node
            # isLeaf    1 byte
            # reserved  1 byte
            # count     2 bytes
            isLeaf, count = struct.unpack(byteorder_char + "?xh", stream.read(4))
            if isLeaf:
                break
            assert count > 0
            # Supplemental Table 11: Chromosome B+ tree non-leaf item
            # key          keySize bytes
            # childOffset  8 bytes
            key = stream.read(keySize)
            childOffset = int.from_bytes(stream.read(8), byteorder)
            stream.seek(childOffset)

        targets = []

        while True:
            for i in range(count):
                # Supplemental Table 10: Chromosome B+ tree leaf item format
                # key        keySize bytes
                # chromId    4 bytes
                # chromSize  4 bytes
                key = stream.read(keySize)
                chromName = key.split(b"\x00", 1).pop(0)
                chromId, chromSize = struct.unpack(
                    byteorder_char + "II", stream.read(8)
                )
                chromName = chromName.decode()
                sequence = Seq(None, length=chromSize)
                target = SeqRecord(sequence, id=chromName)
                targets.append(target)
            if chromId + 1 == itemCount:
                break
            # Supplemental Table 9: Chromosome B+ tree node
            # isLeaf    1 byte
            # reserved  1 byte
            # count     2 bytes
            isLeaf, count = struct.unpack(byteorder_char + "?xh", stream.read(4))
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
        ) = struct.unpack(byteorder_char + "iqiiiiqixxxx", stream.read(44))

        self.targets = targets
        self._index = 0
        self._length = dataCount
        self._cache = ("", 0)
        self._byteorder_char = byteorder_char

    def _read_next_alignment(self, stream):
        if self._index == self._length:
            return
        byteorder_char = self._byteorder_char
        data, count = self._cache
        if not data:
            if not count:
                while True:
                    # Supplemental Table 15: R tree node format
                    # isLeaf     1 byte
                    # reserved   1 byte
                    # count      2 bytes
                    isLeaf, count = struct.unpack(
                        byteorder_char + "?xh", stream.read(4)
                    )
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
            ) = struct.unpack(byteorder_char + "iiiiqq", stream.read(32))
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
        chromId, chromStart, chromEnd = struct.unpack(byteorder_char + "III", data[:12])
        rest, data = data[12:].split(b"\00", 1)
        words = rest.decode().split("\t")
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
        self._index += 1
        self._cache = (data, count)
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
