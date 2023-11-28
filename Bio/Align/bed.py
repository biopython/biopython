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
import sys
import numpy as np


from Bio.Align import Alignment
from Bio.Align import interfaces
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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
        super().__init__(target)
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
        else:
            if chrom is None:
                chrom = "target"
        if coordinates[0, 0] > coordinates[0, -1]:
            # read the alignment right-to-left to be consistent with BED
            coordinates = coordinates[:, ::-1]
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
        try:
            chromStart = blockStarts[0]  # start of alignment in target
            chromEnd = blockStarts[-1] + blockSize  # end of alignment in target
        except IndexError:  # no aligned blocks
            chromStart = chromEnd = tStart
        fields = [chrom, str(chromStart), str(chromEnd)]
        if bedN == 3:
            return "\t".join(fields) + "\n"
        try:
            name = query.id
        except AttributeError:
            name = "query"
        else:
            if name is None:
                name = "query"
        fields.append(name)
        if bedN == 4:
            return "\t".join(fields) + "\n"
        try:
            score = alignment.score
        except AttributeError:
            score = 0
        fields.append(format(score, "g"))
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
    """Alignment iterator for Browser Extensible Data (BED) files.

    Each line in the file contains one pairwise alignment, which are loaded
    and returned incrementally.  Additional alignment information is stored as
    attributes of each alignment.
    """

    fmt = "BED"

    def _read_next_alignment(self, stream):
        for line in stream:
            # note that we cannot extract one line by calling next, as stream
            # may be iterable but not an iterator (for example, TemporaryFile
            # or NamedTemporaryFile objects in tempfile).
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
                blockSizes = np.array(blockSizes)
                blockStarts = np.array(blockStarts)
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
            query_record = SeqRecord(query_sequence, id=name, description="")
            target_sequence = Seq(None, length=sys.maxsize)
            target_record = SeqRecord(target_sequence, id=chrom, description="")
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
