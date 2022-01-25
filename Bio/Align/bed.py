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
from itertools import chain
import numpy


from Bio.Align import Alignment
from Bio.Align import interfaces
from Bio.SeqRecord import SeqRecord


class AlignmentWriter(interfaces.AlignmentWriter):
    """Alignment file writer for the Browser Extensible Data (BED) file format."""

    def write_alignment(self, alignment):
        """Write a complete alignment as one BED line."""
        if not isinstance(alignment, Alignment):
            raise TypeError("Expected an Alignment object")
        coordinates = alignment.coordinates
        if not coordinates.size:  # alignment consists of gaps only
            return ""
        target, query = alignment.sequences
        try:
            name = query.id
        except AttributeError:
            name = "query"
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
        try:
            score = alignment.score
        except AttributeError:
            score = 0
        chromStart = blockStarts[0]  # start of alignment in target
        chromEnd = blockStarts[-1] + blockSize  # end of alignment in target
        blockStarts -= chromStart
        blockCount = len(blockSizes)
        blockSizes = ",".join(map(str, blockSizes)) + ","
        blockStarts = ",".join(map(str, blockStarts)) + ","
        try:
            thickStart = alignment.thickStart
        except AttributeError:
            thickStart = chromStart
        try:
            thickEnd = alignment.thickEnd
        except AttributeError:
            thickEnd = chromEnd
        try:
            itemRgb = alignment.itemRgb
        except AttributeError:
            itemRgb = "0"
        words = [
            chrom,
            str(chromStart),
            str(chromEnd),
            str(name),
            str(score),
            strand,
            str(thickStart),
            str(thickEnd),
            itemRgb,
            str(blockCount),
            blockSizes,
            blockStarts,
        ]
        line = "\t".join(words) + "\n"
        self.stream.write(line)


class AlignmentIterator(interfaces.AlignmentIterator):
    """Alignment iterator for Browser Extensible Data (BED) files.

    Each line in the file contains one pairwise alignment, which are loaded
    and returned incrementally.  Additional alignment information is stored as
    attributes of each alignment.
    """

    def __init__(self, source):
        """Create an AlignmentIterator object.

        Arguments:
         - source   - input data or file name

        """
        super().__init__(source, mode="t", fmt="BED")
        stream = self.stream
        line = next(stream)
        if line.startswith("psLayout "):
            words = line.split()
            if words[1] != "version":
                raise ValueError("Unexpected word '%s' in header line" % words[1])
            self.metadata = {"version": words[2]}
            line = next(stream)
            line = next(stream)
            line = next(stream)
            line = next(stream)
            if line.lstrip("-").strip() != "":
                raise ValueError("End of header not found")
            self.line = None
        else:
            self.line = line

    def parse(self, stream):
        """Parse the next alignment from the stream."""
        if stream is None:
            raise StopIteration

        line = self.line
        self.line = None
        if line is not None:
            lines = chain([line], stream)
        else:
            lines = stream
        for line in lines:
            words = line.split()
            if len(words) != 12:
                raise ValueError("line has %d columns; expected 21" % len(words))
            chrom = words[0]
            chromStart = int(words[1])
            chromEnd = int(words[2])
            name = words[3]
            score = words[4]
            try:
                score = float(score)
            except ValueError:
                pass
            else:
                if score.is_integer():
                    score = int(score)
            strand = words[5]
            blockCount = int(words[9])
            blockSizes = [
                int(blockSize) for blockSize in words[10].rstrip(",").split(",")
            ]
            blockStarts = [int(blockStart) for blockStart in words[11].rstrip(",").split(",")]
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
            target_record = SeqRecord(None, id=chrom)
            query_record = SeqRecord(None, id=name)
            records = [target_record, query_record]
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
            coordinates[0, :] += chromStart
            if strand == "-":
                coordinates[1, :] = sum(blockSizes) - coordinates[1, :]
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
            alignment.score = score
            yield alignment
