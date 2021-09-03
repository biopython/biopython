# Copyright 2008-2016 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.Align support for output from William Pearson's FASTA alignment tools.

This module contains a parser for output from the FASTA programs generated with
the '-m 8CB' or '-m 8CC' output formats.
"""
import re
import enum
import numpy
from Bio.Align import Alignment
from Bio.Align import interfaces
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class State(enum.Enum):
    MATCH = enum.auto()
    QUERY_GAP = enum.auto()
    TARGET_GAP = enum.auto()
    NONE = enum.auto()


class AlignmentIterator(interfaces.AlignmentIterator):
    """FASTA output alignment iterator.

    For reading the (pairwise) alignments from the FASTA alignment programs
    using the '-m 8CB' or '-m 8CC' output formats.
    """

    def __init__(self, source):
        """Create an Iterator object.

        Arguments:
         - source   - input data or file name

        """
        super().__init__(source, mode="t", fmt="EMBOSS")
        stream = self.stream
        try:
            line = next(stream)
        except StopIteration:
            raise ValueError("Empty file.") from None

        assert line.startswith("# ")
        self.commandline = line[2:].strip()

    def parse(self, stream):
        """Parse the next alignment from the stream."""
        if stream is None:
            raise StopIteration

        for line in stream:
            if line.startswith("# "):
                line = line.strip()
                if line.startswith("# FASTA processed ") and line.endswith(" queries"):
                    return
                self._program = line[2:]
                line = next(stream)
                prefix = "# Query: "
                assert line.startswith(prefix)
                query_line, query_size = line[len(prefix) :].strip().rsplit(" - ", 1)
                query_size, unit = query_size.split()
                self._query_size = int(query_size)
                assert unit in ("nt", "aa")
                try:
                    self._query_id, self._query_description = query_line.split(None, 1)
                except ValueError:
                    self._query_id = query_line.strip()
                    self._query_description = None
                line = next(stream)
                prefix = "# Database: "
                assert line.startswith(prefix)
                self._database = line[len(prefix) :].strip()
                line = next(stream)
                prefix = "# Fields: "
                assert line.startswith(prefix)
                fields = line[len(prefix) :].strip().split(", ")
                assert fields[0] == "query id"
                assert fields[1] == "subject id"
                assert fields[2] == "% identity"
                assert fields[3] == "alignment length"
                assert fields[4] == "mismatches"
                assert fields[5] == "gap opens"
                assert fields[6] == "q. start"
                assert fields[7] == "q. end"
                assert fields[8] == "s. start"
                assert fields[9] == "s. end"
                assert fields[10] == "evalue"
                assert fields[11] == "bit score"
                if fields[12] == "BTOP":
                    self._alignment_representation = "BTOP"
                elif fields[12] == "aln_code":
                    self._alignment_representation = "CIGAR"
                else:
                    raise ValueError("Unexpected field '%s'" % fields[12])
                line = next(stream)
                line = line.strip()
                assert line.startswith("# ")
                suffix = " hits found"
                assert line.endswith(suffix)
                hits = int(line[2 : -len(suffix)])
            else:
                yield self.create_alignment(line)

    def create_alignment(self, line):
        columns = line.split()
        assert len(columns) == 13
        annotations = {}
        annotations["program"] = self._program
        annotations["database"] = self._database
        if self._query_id is not None:
            assert columns[0] == self._query_id
        query_id = columns[0]
        target_id = columns[1]
        percentage_identity = float(columns[2])
        alignment_length = int(columns[3])
        mismatches = int(columns[4])
        matches = alignment_length - mismatches
        difference = abs(100 * matches / alignment_length - percentage_identity)
        assert difference < 0.015
        gap_opens = int(columns[5])
        query_start = int(columns[6]) - 1
        query_end = int(columns[7])
        target_start = int(columns[8]) - 1
        target_end = int(columns[9])
        annotations["mismatches"] = mismatches
        annotations["evalue"] = float(columns[10])
        annotations["bit_score"] = float(columns[11])
        if self._alignment_representation == "BTOP":
            coordinates = self.parse_btop(columns[12])
        elif self._alignment_representation == "CIGAR":
            coordinates = self.parse_cigar(columns[12])
        coordinates[0, :] += target_start
        if query_start < query_end:
            coordinates[1, :] += query_start
        else:
            # mapped to reverse strand
            coordinates[1, :] = query_start - coordinates[1, :] + 1
        if self._query_size is None:
            query_size = query_end
        else:
            query_size = self._query_size
        query_sequence = Seq(None, length=query_size)
        query = SeqRecord(query_sequence, id=query_id)
        if self._query_description is not None:
            query.description = self._query_description
        target_sequence = Seq(None, length=target_end)
        target = SeqRecord(target_sequence, id=target_id)
        records = [target, query]
        alignment = Alignment(records, coordinates)
        alignment.annotations = annotations
        return alignment

    def parse_btop(self, btop):
        target_coordinates = []
        query_coordinates = []
        target_coordinates.append(0)
        query_coordinates.append(0)
        state = State.NONE
        tokens = re.findall("([A-Z-]{2}|\d+)", btop)
        # each token is now
        # - an integer
        # - a pair of characters, which may include dashes
        for token in tokens:
            if token.startswith("-"):
                if state != State.QUERY_GAP:
                    target_coordinates.append(target_coordinates[-1])
                    query_coordinates.append(query_coordinates[-1])
                    state = State.QUERY_GAP
                target_coordinates[-1] += 1
            elif token.endswith("-"):
                if state != State.TARGET_GAP:
                    target_coordinates.append(target_coordinates[-1])
                    query_coordinates.append(query_coordinates[-1])
                    state = State.TARGET_GAP
                query_coordinates[-1] += 1
            else:
                try:
                    length = int(token)
                except ValueError:
                    # pair of mismatched letters
                    length = 1
                if state == State.MATCH:
                    target_coordinates[-1] += length
                    query_coordinates[-1] += length
                else:
                    target_coordinates.append(target_coordinates[-1] + length)
                    query_coordinates.append(query_coordinates[-1] + length)
                    state = State.MATCH
        coordinates = numpy.array([target_coordinates, query_coordinates])
        return coordinates

    def parse_cigar(self, cigar):
        target_coordinates = []
        query_coordinates = []
        target_coordinate = 0
        query_coordinate = 0
        target_coordinates.append(target_coordinate)
        query_coordinates.append(query_coordinate)
        state = State.NONE
        tokens = re.findall("(M|D|I|\d+)", cigar)
        # each token is now
        # - the length of the operation
        # - the operation
        for length, operation in zip(tokens[::2], tokens[1::2]):
            length = int(length)
            if operation == "M":
                target_coordinate += length
                query_coordinate += length
            elif operation == "I":
                target_coordinate += length
            elif operation == "D":
                query_coordinate += length
            else:
                raise ValueError("Unexpected operation '%s'" % operation)
            target_coordinates.append(target_coordinate)
            query_coordinates.append(query_coordinate)
        coordinates = numpy.array([target_coordinates, query_coordinates])
        return coordinates
