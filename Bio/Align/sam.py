# Copyright 2022 by Michiel de Hoon.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.Align support for the "sam" pairwise alignment format.

The Sequence Alignment/Map (SAM) format, created by Heng Li and Richard Durbin
at the Wellcome Trust Sanger Institute, stores a series of alignments to the
genome in a single file. Typically they are used for next-generation sequencing
data. SAM files store the alignment positions for mapped sequences, and may
also store the aligned sequences and other information associated with the
sequence.

See http://www.htslib.org/ for more information.

You are expected to use this module via the Bio.Align functions.

Coordinates in the SAM format are defined in terms of one-based start
positions; the parser converts these to zero-based coordinates to be consistent
with Python and other alignment formats.
"""
from itertools import chain
import numpy


from Bio.Align import Alignment
from Bio.Align import interfaces
from Bio.Seq import Seq, reverse_complement, UndefinedSequenceError
from Bio.SeqRecord import SeqRecord


class AlignmentWriter(interfaces.AlignmentWriter):
    """Alignment file writer for the Sequence Alignment/Map (SAM) ile format."""

    def write_header(self, alignments):
        """Write the SAM header."""
        if not self.header:
            return
        try:
            metadata = alignments.metadata
        except AttributeError:
            version = "3"
        else:
            version = metadata.get("version", "3")
        # fmt: off
        self.stream.write(
            f"""\
psLayout version {version}

match	mis- 	rep. 	N's	Q gap	Q gap	T gap	T gap	strand	Q        	Q   	Q    	Q  	T        	T   	T    	T  	block	blockSizes 	qStarts	 tStarts
     	match	match	   	count	bases	count	bases	      	name     	size	start	end	name     	size	start	end	count
---------------------------------------------------------------------------------------------------------------------------------------------------------------
"""  # noqa: W191, E101
        )
        # fmt: on

    def format_alignment(self, alignment):
        """Return a string with a single alignment formatted as one SAM line."""
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
        # fmt: off
        dnax = None  # set to True for translated DNA aligned to protein,
                     # and to False for DNA/RNA aligned to DNA/RNA  # noqa: E114, E116
        if coordinates[1, 0] > coordinates[1, -1]:
            # DNA/RNA mapped to reverse strand of DNA/RNA
            strand = "-"
            query = reverse_complement(query, inplace=False)
            coordinates = coordinates.copy()
            coordinates[1, :] = qSize - coordinates[1, :]
        elif coordinates[0, 0] > coordinates[0, -1]:
            # protein mapped to reverse strand of DNA
            strand = "-"
            target = reverse_complement(target, inplace=False)
            coordinates = coordinates.copy()
            coordinates[0, :] = tSize - coordinates[0, :]
            dnax = True
        else:
            # mapped to forward strand
            strand = "+"
        # fmt: on
        wildcard = self.wildcard
        mask = self.mask
        # variable names follow those in the SAM file format specification
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
                    assert dnax is not True
                    dnax = False
                else:
                    # translated DNA aligned to protein, typically generated by
                    # blat -t=dnax -q=prot
                    assert tCount == 3 * qCount
                    assert dnax is not False
                    dnax = True
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
    """Alignment iterator for Sequence Alignment/Map (SAM) files.

    Each line in the file contains one genomic alignment, which are loaded
    and returned incrementally.  The following columns are stored as attributes
    of the alignment:

      - flag: The FLAG combination of bitwise flags;
      - mapq: Mapping Quality (only stored if available)
      - rnext: Reference sequence name of the primary alignment of the next read
               in the alignment (only stored if available)
      - pnext: Zero-based position of the primary alignment of the next read in
               the template (only stored if available)
      - tlen: signed observed template length (only stored if available)
      - qual: ASCII of base quality plus 33 (only stored if available)

    Other information associated with the alignment by its tags are stored in
    the annotations attribute of each alignment.
    """

    def __init__(self, source):
        """Create an AlignmentIterator object.

        Arguments:
         - source   - input data or file name

        """
        super().__init__(source, mode="t", fmt="SAM")
        stream = self.stream
        self.metadata = {}
        self.targets = {}
        line = next(stream)
        if not line.startswith("@"):
            raise ValueError("Header is missing")
        self._parse_header_line(line)
        for line in stream:
            if not line.startswith("@"):
                self.line = line
                break
            self._parse_header_line(line)

    def _parse_header_line(self, line):
        words = line[1:].strip().split("\t")
        tag = words[0]
        values = {}
        if tag == "SQ":
            annotations = {}
            description = None
            for word in words[1:]:
                key, value = word.split(":", 1)
                assert len(key) == 2
                if key == "SN":
                    rname = value
                elif key == "LN":
                    length = int(value)
                elif key == "AH":
                    annotations["alternate_locus"] = value
                elif key == "AH":
                    annotations["names"] = value.split(",")
                elif key == "AS":
                    annotations["assembly"] = value
                elif key == "DS":
                    description = value
                elif key == "M5":
                    annotations["MD5"] = value
                elif key == "SP":
                    annotations["species"] = value
                elif key == "TP":
                    assert value in ("linear", "circular")
                    annotations["topology"] = value
                elif key == "URI":
                    annotations["URI"] = value
                else:
                    annotations[key] = value
            assert rname not in self.targets
            sequence = Seq(None, length=length)
            record = SeqRecord(sequence, id=rname, annotations=annotations)
            if description is not None:
                record.description = description
            self.targets[rname] = record
        else:
            for word in words[1:]:
                key, value = word.split(":", 1)
                assert len(key) == 2
                values[key] = value
            if tag == "HD":
                self.metadata[tag] = values
            else:
                if tag not in self.metadata:
                    self.metadata[tag] = []
                self.metadata[tag].append(values)

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
            if len(words) < 11:
                raise ValueError("line has %d columns; expected at least 11" % len(words))
            qname = words[0]
            flag = int(words[1])
            rname = words[2]
            pos = int(words[3]) - 1
            mapq = int(words[4])
            cigar = words[5]
            rnext = words[6]
            pnext = int(words[7]) - 1
            tlen = int(words[8])
            seq = words[9]
            qual = words[10]
            number = ""
            query_pos = 0
            coordinates = [[pos, query_pos]]
            for letter in cigar:
                if letter in "M=X":
                    number = int(number)
                    pos += number
                    query_pos += number
                    coordinates.append([pos, query_pos])
                    number = ""
                elif letter in "IS":
                    number = int(number)
                    query_pos += number
                    coordinates.append([pos, query_pos])
                    number = ""
                elif letter in "DN":
                    number = int(number)
                    pos += number
                    coordinates.append([pos, query_pos])
                    number = ""
                elif letter == "H":
                    # hard clipping (clipped sequences not present in sequence)
                    number = int(number)
                    number = ""
                elif letter == "P":
                    number = int(number)
                    raise NotImplementedError("padding operator is not yet implemented")
                    number = ""
                else:
                    number += letter
            coordinates = numpy.array(coordinates).transpose()
            target = self.targets[rname]
            if seq == "*":
                sequence = Seq(None, length=query_pos)
            else:
                sequence = Seq(seq)
                if flag & 0x10:
                    sequence = sequence.reverse_complement()
            if flag & 0x10:
                coordinates[1, :] = query_pos - coordinates[1, :]
            query = SeqRecord(sequence, id=qname)
            records = [target, query]
            alignment = Alignment(records, coordinates)
            alignment.flag = flag
            if mapq != 255:
                alignment.mapq = mapq
            if rnext == "=":
                alignment.rnext = rname
            elif rnext != "*":
                alignment.rnext = rnext
            if pnext >= 0:
                alignment.pnext = pnext
            if tlen > 0:
                alignment.tlen = tlen
            if qual != "*":
                alignment.qual = qual
            yield alignment
