# Copyright 2022 by Michiel de Hoon.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.Align support for Exonerate output format.

This module provides support for Exonerate outputs. Exonerate is a generic
tool for pairwise sequence comparison that allows you to align sequences using
several different models.

Bio.Align.exonerate was tested on the following Exonerate versions and models:

    - version: 2.2
    - models:
      - affine:local                - cdna2genome
      - coding2coding               - est2genome
      - genome2genome               - ner
      - protein2dna                 - protein2genome
      - ungapped                    - ungapped:translated

Although model testing were not exhaustive, ExonerateIO should be able to cope
with all Exonerate models. Please file a bug report if you stumble upon an
unparsable file.

More information on Exonerate is available on its home page at
www.ebi.ac.uk/~guy/exonerate/

You are expected to use this module via the Bio.Align functions.
"""
from itertools import chain
import numpy


from Bio.Align import Alignment
from Bio.Align import interfaces
from Bio.Seq import Seq, reverse_complement, UndefinedSequenceError
from Bio.SeqRecord import SeqRecord
from Bio import BiopythonExperimentalWarning

import warnings

warnings.warn(
    "Bio.Align.exonerate is an experimental module which may undergo "
    "significant changes prior to its future official release.",
    BiopythonExperimentalWarning,
)


class AlignmentWriter(interfaces.AlignmentWriter):
    """Alignment file writer for the Exonerate cigar and vulgar file format."""

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
        super().__init__(target, mode="w")
        self.header = header
        if wildcard is not None:
            if mask == "upper":
                wildcard = ord(wildcard.lower())
            else:
                wildcard = ord(wildcard.upper())
        self.wildcard = wildcard
        self.mask = mask

    def write_header(self, alignments):
        """Write the PSL header."""
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
    """Alignment iterator for the Exonerate text, cigar, and vulgar formats.

    Each line in the file contains one pairwise alignment, which are loaded
    and returned incrementally.  Alignment score information such as the number
    of matches and mismatches are stored as attributes of each alignment.
    """

    def __init__(self, source):
        """Create an AlignmentIterator object.

        Arguments:
         - source   - input data or file name

        """
        super().__init__(source, mode="t", fmt="PSL")
        stream = self.stream
        self.program = "exonerate"
        line = next(stream)
        prefix = "Command line: "
        assert line.startswith(prefix)
        commandline = line[len(prefix):].strip()
        assert commandline.startswith("[")
        assert commandline.endswith("]")
        self.commandline = commandline[1:-1]
        line = next(stream)
        prefix = "Hostname: "
        assert line.startswith(prefix)
        hostname = line[len(prefix):].strip()
        assert hostname.startswith("[")
        assert hostname.endswith("]")
        self.hostname = hostname[1:-1]

    @staticmethod
    def _parse_cigar(words):
        query_id = words[0]
        query_start = int(words[1])
        query_end = int(words[2])
        query_strand = words[3]
        target_id = words[4]
        target_start = int(words[5])
        target_end = int(words[6])
        target_strand = words[7]
        score = int(words[8])
        target_seq = Seq(None, length=target_end)
        query_seq = Seq(None, length=query_end)
        target = SeqRecord(target_seq, id=target_id)
        query = SeqRecord(query_seq, id=query_id)
        qs = 0
        ts = 0
        n = (len(words) - 8) // 2
        coordinates =  numpy.empty((2, n+1), int)
        coordinates[0, 0] = ts
        coordinates[1, 0] = qs
        for i, (operation, step) in enumerate(zip(words[9::2], words[10::2])):
            step = int(step)
            if operation == "M":  # match
                ts += step
                qs += step
            elif operation == "I":  # insertion
                if query_strand == "." and target_strand != ".":
                    qs += step * 3
                else:
                    qs += step
            elif operation == "D":  # deletion
                if target_strand == "." and query_strand != ".":
                    ts += step * 3
                else:
                    ts += step
            else:
                raise ValueError("Unknown operation %s in cigar string" % operation)
            coordinates[0, i+1] = ts
            coordinates[1, i+1] = qs
        if target_strand == "+":
            coordinates[0, :] += target_start
        elif target_strand == "-":
            coordinates[0, :] = target_start - coordinates[0, :]
        elif target_strand == ".":  # protein
            if query_strand != ".":
                # dna to protein alignment; integer division, but round up:
                coordinates[0, :] = (coordinates[0, :] + 2 ) // 3
            coordinates[0, :] += target_start
        if query_strand == "+":
            coordinates[1, :] += query_start
        elif query_strand == "-":
            coordinates[1, :] = query_start - coordinates[1, :]
        elif query_strand == ".":  # protein
            if target_strand != ".":
                # protein to dna alignment; integer division, but round up:
                coordinates[1, :] = (coordinates[1, :] + 2 ) // 3
            coordinates[1, :] += query_start
        alignment = Alignment([target, query], coordinates)
        alignment.score = score
        return alignment

    @staticmethod
    def _parse_vulgar(words):
        query_id = words[0]
        query_start = int(words[1])
        query_end = int(words[2])
        query_strand = words[3]
        target_id = words[4]
        target_start = int(words[5])
        target_end = int(words[6])
        target_strand = words[7]
        score = int(words[8])
        target_seq = Seq(None, length=target_end)
        query_seq = Seq(None, length=query_end)
        target = SeqRecord(target_seq, id=target_id)
        query = SeqRecord(query_seq, id=query_id)
        qs = 0
        ts = 0
        n = (len(words) - 8) // 3
        coordinates =  numpy.empty((2, n+1), int)
        coordinates[0, 0] = ts
        coordinates[1, 0] = qs
        operations = bytearray(n)
        for i, (operation, query_step, target_step) in enumerate(zip(words[9::3], words[10::3], words[11::3])):
            query_step = int(query_step)
            target_step = int(target_step)
            ts += target_step
            qs += query_step
            if operation == "M":  # Match
                pass
            elif operation == "5":  # 5' splice site
                pass
            elif operation == "I":  # Intron
                operation = "N"  # use SAM/BAM definitions of operations
            elif operation == "3":  # 3' splice site
                pass
            elif operation == "C":  # Codon
                pass
            elif operation == "G":  # Gap
                pass
            elif operation == "N":  # Non-equivalenced region
                operation = "R"  # 'N' is alread used for introns in SAM/BAM
            elif operation == "S":  # Split codon
                pass
            elif operation == "F":  # Frame shift
                pass
            else:
                raise ValueError("Unknown operation %s in vulgar string" % operation)
            coordinates[0, i+1] = ts
            coordinates[1, i+1] = qs
            operations[i] = ord(operation)
        if target_strand == "+":
            coordinates[0, :] += target_start
        elif target_strand == "-":
            coordinates[0, :] = target_start - coordinates[0, :]
        elif target_strand == ".":  # protein
            coordinates[0, :] += target_start
        if query_strand == "+":
            coordinates[1, :] += query_start
        elif query_strand == "-":
            coordinates[1, :] = query_start - coordinates[1, :]
        elif query_strand == ".":  # protein
            coordinates[1, :] += query_start
        alignment = Alignment([target, query], coordinates)
        alignment.operations = operations
        alignment.score = score
        return alignment

    def parse(self, stream):
        """Parse the next alignment from the stream."""
        if stream is None:
            raise StopIteration

        for line in stream:
            line = line.strip()
            if line == "-- completed exonerate analysis":
                try:
                    next(stream)
                except StopIteration:
                    return
                raise ValueError("Found additional data after 'completed exonerate analysis'; corrupt file?")
            elif line.startswith("vulgar: "):
                words = line[8:].split()
                alignment = self._parse_vulgar(words)
                yield alignment
            elif line.startswith("cigar: "):
                words = line[7:].split()
                alignment = self._parse_cigar(words)
                yield alignment
        raise ValueError("Failed to find 'completed exonerate analysis'; truncated file?")
