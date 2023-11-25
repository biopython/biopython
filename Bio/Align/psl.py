# Copyright 2022 by Michiel de Hoon.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.Align support for the "psl" pairwise alignment format.

The Pattern Space Layout (PSL) format, described by UCSC, stores a series
of pairwise alignments in a single file. Typically they are used for
transcript to genome alignments. PSL files store the alignment positions
and alignment scores, but do not store the aligned sequences.

See http://genome.ucsc.edu/FAQ/FAQformat.html#format2

You are expected to use this module via the Bio.Align functions.

Coordinates in the PSL format are defined in terms of zero-based start
positions (like Python) and aligning region sizes.

A minimal aligned region of length one and starting at first position in the
source sequence would have ``start == 0`` and ``size == 1``.

As we can see in this example, ``start + size`` will give one more than the
zero-based end position. We can therefore manipulate ``start`` and
``start + size`` as python list slice boundaries.
"""
from itertools import chain
import numpy as np


from Bio.Align import Alignment
from Bio.Align import interfaces
from Bio.Seq import Seq, reverse_complement, UndefinedSequenceError
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, ExactPosition, SimpleLocation, CompoundLocation


class AlignmentWriter(interfaces.AlignmentWriter):
    """Alignment file writer for the Pattern Space Layout (PSL) file format."""

    fmt = "PSL"

    def __init__(self, target, header=True, mask=None, wildcard="N"):
        """Create an AlignmentWriter object.

        Arguments:
         - target    - output stream or file name
         - header    - If True (default), write the PSL header consisting of
                       five lines containing the PSL format version and a
                       header for each column.
                       If False, suppress the PSL header, resulting in a simple
                       tab-delimited file.
         - mask      - Specify if repeat regions in the target sequence are
                       masked and should be reported in the `repMatches` field
                       of the PSL file instead of in the `matches` field.
                       Acceptable values are
                       None   : no masking (default);
                       "lower": masking by lower-case characters;
                       "upper": masking by upper-case characters.
         - wildcard  - Report alignments to the wildcard character in the
                       target or query sequence in the `nCount` field of the
                       PSL file instead of in the `matches`, `misMatches`, or
                       `repMatches` fields.
                       Default value is 'N'.

        """
        super().__init__(target)
        self.header = header
        if wildcard is not None:
            if mask == "upper":
                wildcard = ord(wildcard.lower())
            else:
                wildcard = ord(wildcard.upper())
        self.wildcard = wildcard
        self.mask = mask

    def write_header(self, stream, alignments):
        """Write the PSL header."""
        if not self.header:
            return
        try:
            metadata = alignments.metadata
        except AttributeError:
            version = "3"
        else:
            version = metadata.get("psLayout version", "3")
        # fmt: off
        stream.write(
            f"""\
psLayout version {version}

match	mis- 	rep. 	N's	Q gap	Q gap	T gap	T gap	strand	Q        	Q   	Q    	Q  	T        	T   	T    	T  	block	blockSizes 	qStarts	 tStarts
     	match	match	   	count	bases	count	bases	      	name     	size	start	end	name     	size	start	end	count
---------------------------------------------------------------------------------------------------------------------------------------------------------------
"""  # noqa: W191, E101
        )
        # fmt: on

    def format_alignment(self, alignment):
        """Return a string with a single alignment formatted as one PSL line."""
        if not isinstance(alignment, Alignment):
            raise TypeError("Expected an Alignment object")
        coordinates = alignment.coordinates
        if not coordinates.size:  # alignment consists of gaps only
            return ""
        target, query = alignment.sequences
        try:
            qName = query.id
        except AttributeError:
            qName = "query"
        try:
            query = query.seq
        except AttributeError:
            pass
        try:
            tName = target.id
        except AttributeError:
            tName = "target"
        try:
            target = target.seq
        except AttributeError:
            pass
        tSize = len(target)
        qSize = len(query)
        steps = abs(np.diff(coordinates, 1))
        aligned = steps.min(0) > 0
        tCount, qCount = steps[:, aligned].sum(1)
        if tCount == qCount:
            dnax = False
            if coordinates[0, 0] > coordinates[0, -1]:
                coordinates = coordinates[:, ::-1]
            if coordinates[1, 0] > coordinates[1, -1]:
                # DNA/RNA mapped to reverse strand of DNA/RNA
                strand = "-"
                query = reverse_complement(query)
                coordinates = coordinates.copy()
                coordinates[1, :] = qSize - coordinates[1, :]
            else:
                # mapped to forward strand
                strand = "+"
        elif tCount == 3 * qCount:
            dnax = True
            if coordinates[0, 0] > coordinates[0, -1]:
                # protein mapped to reverse strand of DNA
                strand = "-"
                target = reverse_complement(target)
                coordinates = coordinates.copy()
                coordinates[0, :] = tSize - coordinates[0, :]
            else:
                # mapped to forward strand
                strand = "+"
        else:
            raise ValueError("inconsistent step sizes %d, %d" % (tCount, qCount))
        # fmt: on
        wildcard = self.wildcard
        mask = self.mask
        # variable names follow those in the PSL file format specification
        matches = 0
        misMatches = 0
        repMatches = 0
        nCount = 0
        qNumInsert = 0
        qBaseInsert = 0
        tNumInsert = 0
        tBaseInsert = 0
        blockSizes = []
        qStarts = []
        tStarts = []
        tStart, qStart = coordinates[:, 0]
        for tEnd, qEnd in coordinates[:, 1:].transpose():
            if tStart == tEnd:
                if qStart > 0 and qEnd < qSize:
                    qNumInsert += 1
                    qBaseInsert += qEnd - qStart
                qStart = qEnd
            elif qStart == qEnd:
                if tStart > 0 and tEnd < tSize:
                    tNumInsert += 1
                    tBaseInsert += tEnd - tStart
                tStart = tEnd
            else:
                tCount = tEnd - tStart
                qCount = qEnd - qStart
                tStarts.append(tStart)
                qStarts.append(qStart)
                blockSizes.append(qCount)
                if tCount == qCount:
                    assert dnax is False
                else:
                    # translated DNA aligned to protein, typically generated by
                    # blat -t=dnax -q=prot
                    assert tCount == 3 * qCount
                    assert dnax is True
                tSeq = target[tStart:tEnd]
                qSeq = query[qStart:qEnd]
                try:
                    tSeq = bytes(tSeq)
                except TypeError:  # string
                    tSeq = bytes(tSeq, "ASCII")
                except UndefinedSequenceError:  # sequence contents is unknown
                    tSeq = None
                try:
                    qSeq = bytes(qSeq)
                except TypeError:  # string
                    qSeq = bytes(qSeq, "ASCII")
                except UndefinedSequenceError:  # sequence contents is unknown
                    qSeq = None
                if tSeq is None or qSeq is None:
                    # contents of at least one sequence is unknown;
                    # count all aligned letters as matches:
                    matches += qCount
                else:
                    if mask == "lower":
                        for u1, u2, c1 in zip(tSeq.upper(), qSeq.upper(), tSeq):
                            if u1 == wildcard or u2 == wildcard:
                                nCount += 1
                            elif u1 == u2:
                                if u1 == c1:
                                    matches += 1
                                else:
                                    repMatches += 1
                            else:
                                misMatches += 1
                    elif mask == "upper":
                        for u1, u2, c1 in zip(tSeq.lower(), qSeq.lower(), tSeq):
                            if u1 == wildcard or u2 == wildcard:
                                nCount += 1
                            elif u1 == u2:
                                if u1 == c1:
                                    matches += 1
                                else:
                                    repMatches += 1
                            else:
                                misMatches += 1
                    else:
                        for u1, u2 in zip(tSeq.upper(), qSeq.upper()):
                            if u1 == wildcard or u2 == wildcard:
                                nCount += 1
                            elif u1 == u2:
                                matches += 1
                            else:
                                misMatches += 1
                tStart = tEnd
                qStart = qEnd
        try:
            matches = alignment.matches
        except AttributeError:
            pass
        try:
            misMatches = alignment.misMatches
        except AttributeError:
            pass
        try:
            repMatches = alignment.repMatches
        except AttributeError:
            pass
        try:
            nCount = alignment.nCount
        except AttributeError:
            pass
        tStart = tStarts[0]  # start of alignment in target
        qStart = qStarts[0]  # start of alignment in query
        tEnd = tStarts[-1] + tCount  # end of alignment in target
        qEnd = qStarts[-1] + qCount  # end of alignment in query
        if strand == "-":
            if dnax is True:
                tStart, tEnd = tSize - tEnd, tSize - tStart
            else:
                qStart, qEnd = qSize - qEnd, qSize - qStart
        blockCount = len(blockSizes)
        blockSizes = ",".join(map(str, blockSizes)) + ","
        qStarts = ",".join(map(str, qStarts)) + ","
        tStarts = ",".join(map(str, tStarts)) + ","
        if dnax:
            strand = "+" + strand
        words = [
            str(matches),
            str(misMatches),
            str(repMatches),
            str(nCount),
            str(qNumInsert),
            str(qBaseInsert),
            str(tNumInsert),
            str(tBaseInsert),
            strand,
            qName,
            str(qSize),
            str(qStart),
            str(qEnd),
            tName,
            str(tSize),
            str(tStart),
            str(tEnd),
            str(blockCount),
            blockSizes,
            qStarts,
            tStarts,
        ]
        line = "\t".join(words) + "\n"
        return line


class AlignmentIterator(interfaces.AlignmentIterator):
    """Alignment iterator for Pattern Space Layout (PSL) files.

    Each line in the file contains one pairwise alignment, which are loaded
    and returned incrementally.  Alignment score information such as the number
    of matches and mismatches are stored as attributes of each alignment.
    """

    fmt = "PSL"

    def _read_header(self, stream):
        line = next(stream)
        if line.startswith("psLayout "):
            words = line.split()
            if words[1] != "version":
                raise ValueError("Unexpected word '%s' in header line" % words[1])
            self.metadata = {"psLayout version": words[2]}
            line = next(stream)
            line = next(stream)
            line = next(stream)
            line = next(stream)
            if line.lstrip("-").strip() != "":
                raise ValueError("End of header not found")
        else:
            self._line = line

    def _read_next_alignment(self, stream):
        try:
            line = self._line
        except AttributeError:
            lines = stream
        else:
            del self._line
            lines = chain([line], stream)
        for line in lines:
            words = line.split()
            if len(words) == 23:
                pslx = True
            elif len(words) == 21:
                pslx = False
            else:
                raise ValueError("line has %d columns; expected 21 or 23" % len(words))
            strand = words[8]
            qName = words[9]
            qSize = int(words[10])
            tName = words[13]
            tSize = int(words[14])
            blockCount = int(words[17])
            blockSizes = [
                int(blockSize) for blockSize in words[18].rstrip(",").split(",")
            ]
            qStarts = [int(start) for start in words[19].rstrip(",").split(",")]
            tStarts = [int(start) for start in words[20].rstrip(",").split(",")]
            if len(blockSizes) != blockCount:
                raise ValueError(
                    "Inconsistent number of blocks (%d found, expected %d)"
                    % (len(blockSizes), blockCount)
                )
            if len(qStarts) != blockCount:
                raise ValueError(
                    "Inconsistent number of query start positions (%d found, expected %d)"
                    % (len(qStarts), blockCount)
                )
            if len(tStarts) != blockCount:
                raise ValueError(
                    "Inconsistent number of target start positions (%d found, expected %d)"
                    % (len(tStarts), blockCount)
                )
            qStarts = np.array(qStarts)
            tStarts = np.array(tStarts)
            qBlockSizes = np.array(blockSizes)
            if strand in ("++", "+-"):
                # protein sequence aligned against translated DNA sequence
                tBlockSizes = 3 * qBlockSizes
            else:
                tBlockSizes = qBlockSizes
            qPosition = qStarts[0]
            tPosition = tStarts[0]
            coordinates = [[tPosition, qPosition]]
            for tBlockSize, qBlockSize, tStart, qStart in zip(
                tBlockSizes, qBlockSizes, tStarts, qStarts
            ):
                if tStart != tPosition:
                    coordinates.append([tStart, qPosition])
                    tPosition = tStart
                if qStart != qPosition:
                    coordinates.append([tPosition, qStart])
                    qPosition = qStart
                tPosition += tBlockSize
                qPosition += qBlockSize
                coordinates.append([tPosition, qPosition])
            coordinates = np.array(coordinates).transpose()
            qNumInsert = 0
            qBaseInsert = 0
            tNumInsert = 0
            tBaseInsert = 0
            tStart, qStart = coordinates[:, 0]
            for tEnd, qEnd in coordinates[:, 1:].transpose():
                tCount = tEnd - tStart
                qCount = qEnd - qStart
                if tCount == 0:
                    if qStart > 0 and qEnd < qSize:
                        qNumInsert += 1
                        qBaseInsert += qCount
                    qStart = qEnd
                elif qCount == 0:
                    if tStart > 0 and tEnd < tSize:
                        tNumInsert += 1
                        tBaseInsert += tCount
                    tStart = tEnd
                else:
                    tStart = tEnd
                    qStart = qEnd
            if qNumInsert != int(words[4]):
                raise ValueError(
                    "Inconsistent qNumInsert found (%s, expected %d)"
                    % (words[4], qNumInsert)
                )
            if qBaseInsert != int(words[5]):
                raise ValueError(
                    "Inconsistent qBaseInsert found (%s, expected %d)"
                    % (words[5], qBaseInsert)
                )
            if tNumInsert != int(words[6]):
                raise ValueError(
                    "Inconsistent tNumInsert found (%s, expected %d)"
                    % (words[6], tNumInsert)
                )
            if tBaseInsert != int(words[7]):
                raise ValueError(
                    "Inconsistent tBaseInsert found (%s, expected %d)"
                    % (words[7], tBaseInsert)
                )
            qStart = int(words[11])
            qEnd = int(words[12])
            tStart = int(words[15])
            tEnd = int(words[16])
            if strand == "-":
                qStart, qEnd = qEnd, qStart
                coordinates[1, :] = qSize - coordinates[1, :]
            elif strand == "+-":
                tStart, tEnd = tEnd, tStart
                coordinates[0, :] = tSize - coordinates[0, :]
            if tStart != coordinates[0, 0]:
                raise ValueError(
                    "Inconsistent tStart found (%d, expected %d)"
                    % (tStart, coordinates[0, 0])
                )
            if tEnd != coordinates[0, -1]:
                raise ValueError(
                    "Inconsistent tEnd found (%d, expected %d)"
                    % (tEnd, coordinates[0, -1])
                )
            if qStart != coordinates[1, 0]:
                raise ValueError(
                    "Inconsistent qStart found (%d, expected %d)"
                    % (qStart, coordinates[1, 0])
                )
            if qEnd != coordinates[1, -1]:
                raise ValueError(
                    "Inconsistent qEnd found (%d, expected %d)"
                    % (qEnd, coordinates[1, -1])
                )
            feature = None
            if pslx is True:
                qSeqs = words[21].rstrip(",").split(",")
                tSeqs = words[22].rstrip(",").split(",")
                qSeq = dict(zip(qStarts, qSeqs))
                if strand in ("++", "+-"):
                    # protein sequence aligned against translated DNA sequence
                    target_sequence = Seq(None, length=tSize)
                    query_sequence = Seq(qSeq, length=qSize)
                    if strand == "++":
                        tStart, qStart = coordinates[:, 0]
                        locations = []
                        for tEnd, qEnd in coordinates[:, 1:].transpose():
                            if qStart < qEnd and tStart < tEnd:
                                location = SimpleLocation(
                                    ExactPosition(tStart),
                                    ExactPosition(tEnd),
                                    strand=+1,
                                )
                                locations.append(location)
                            qStart = qEnd
                            tStart = tEnd
                        if len(locations) > 1:
                            location = CompoundLocation(locations, "join")
                        tSeq = "".join(tSeqs)
                        qualifiers = {"translation": [tSeq]}
                        feature = SeqFeature(
                            location, type="CDS", qualifiers=qualifiers
                        )
                    elif strand == "+-":
                        tEnd, qStart = coordinates[:, 0]
                        locations = []
                        for tStart, qEnd in coordinates[:, 1:].transpose():
                            if qStart < qEnd and tStart < tEnd:
                                location = SimpleLocation(
                                    ExactPosition(tStart),
                                    ExactPosition(tEnd),
                                    strand=-1,
                                )
                                locations.append(location)
                            tEnd = tStart
                            qStart = qEnd
                        if len(locations) > 1:
                            location = CompoundLocation(locations, "join")
                        tSeq = "".join(tSeqs)
                        qualifiers = {"translation": [tSeq]}
                        feature = SeqFeature(
                            location, type="CDS", qualifiers=qualifiers
                        )
                else:
                    tSeq = dict(zip(tStarts, tSeqs))
                    target_sequence = Seq(tSeq, length=tSize)
                    query_sequence = Seq(qSeq, length=qSize)
                    if strand == "-":
                        query_sequence = query_sequence.reverse_complement()
            else:
                target_sequence = Seq(None, length=tSize)
                query_sequence = Seq(None, length=qSize)
            target_record = SeqRecord(target_sequence, id=tName, description="")
            query_record = SeqRecord(query_sequence, id=qName, description="")
            if feature is not None:
                target_record.features.append(feature)
            records = [target_record, query_record]
            alignment = Alignment(records, coordinates)
            alignment.matches = int(words[0])
            alignment.misMatches = int(words[1])
            alignment.repMatches = int(words[2])
            alignment.nCount = int(words[3])
            return alignment
