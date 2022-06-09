# Copyright 2021 by Michiel de Hoon.  All rights reserved.
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
from Bio import BiopythonExperimentalWarning


import warnings

warnings.warn(
    "Bio.Align.tabular is an experimental module which may undergo "
    "significant changes prior to its future official release.",
    BiopythonExperimentalWarning,
)


class State(enum.Enum):
    """Enumerate alignment states needed when parsing a BTOP string."""

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
        """Create an AlignmentIterator object.

        Arguments:
         - source   - input data or file name

        """
        super().__init__(source, mode="t", fmt="FASTA")
        stream = self.stream
        try:
            line = next(stream)
        except StopIteration:
            raise ValueError("Empty file.") from None
        assert line.startswith("# ")
        line = line.rstrip()
        self._parse_header(stream, line)

    def _parse_header(self, stream, line):
        if line[2:].startswith("TBLASTN "):
            self.program, self.version = line[2:].split(None, 1)
            self.final_prefix = "# BLAST processed "
        else:
            self.commandline = line[2:]
            line = next(stream)
            assert line.startswith("# ")
            self.program, self.version = line[2:].rstrip().split(None, 1)
            self.final_prefix = "# FASTA processed "
        line = next(stream)
        line = line.strip()
        prefix = "# Query: "
        assert line.startswith(prefix)
        if self.program == 'FASTA':
            query_line, query_size = line[len(prefix) :].strip().rsplit(" - ", 1)
            query_size, unit = query_size.split()
            self._query_size = int(query_size)
            assert unit in ("nt", "aa")
        else:
            query_line = line[len(prefix) :].strip()
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
        if line.startswith(prefix):
            self.fields = line[len(prefix) :].strip().split(", ")
            line = next(stream)
        line = line.strip()
        assert line.startswith("# ")
        suffix = " hits found"
        assert line.endswith(suffix)
        hits = int(line[2 : -len(suffix)])

    def parse(self, stream):
        """Parse the next alignment from the stream."""
        if stream is None:
            raise StopIteration

        for line in stream:
            line = line.rstrip()
            if line.startswith("# "):
                if line.startswith(self.final_prefix) and line.endswith(" queries"):
                    return
                self._parse_header(stream, line)
            else:
                yield self.create_alignment(line)

    def create_alignment(self, line):
        """Parse one line of output and return an Alignment object."""
        percentage_identity = None
        alignment_length = None
        identical = None
        gaps = None
        mismatches = None
        btop = None
        cigar = None
        score = None
        query_sequence = None
        target_sequence = None
        target_length = None
        try:
            query_size = self._query_size
        except AttributeError:
            query_size = None
        columns = line.split()
        assert len(columns) == len(self.fields)
        annotations = {}
        annotations["database"] = self._database
        for column, field in zip(columns, self.fields):
            if field == "query id":
                query_id = columns[0]
                if self._query_id is not None:
                    assert query_id == self._query_id
            elif field == "subject id":
                target_id = columns[1]
            elif field == "% identity":
                percentage_identity = float(column)
            elif field == "alignment length":
                alignment_length = int(column)
            elif field == "mismatches":
                mismatches = int(column)
                annotations[field] = mismatches
            elif field == "gap opens":
                gap_opens = int(column)
            elif field == "q. start":
                query_start = int(column) - 1
            elif field == "q. end":
                query_end = int(column)
            elif field == "s. start":
                target_start = int(column) - 1
            elif field == "s. end":
                target_end = int(column)
            elif field == "evalue":
                annotations["evalue"] = float(column)
            elif field == "bit score":
                annotations["bit_score"] = float(column)
            elif field == "BTOP":
                coordinates = self.parse_btop(column)
            elif field == "aln_code":
                coordinates = self.parse_cigar(column)
            elif field == "query gi":
                annotations[field] = column
            elif field == "query acc.":
                annotations[field] = column
            elif field == "query acc.ver":
                annotations[field] = column
            elif field == "query length":
                if query_size is None:
                    query_size = int(column)
                else:
                    assert query_size == int(column)
            elif field == "subject ids":
                annotations[field] = column
            elif field == "subject gi":
                annotations[field] = column
            elif field == "subject gis":
                annotations[field] = column
            elif field == "subject acc.":
                annotations[field] = column
            elif field == "subject acc.ver":
                annotations[field] = column
            elif field == "subject length":
                target_length = int(column)
            elif field == "query seq":
                query_sequence = column
            elif field == "subject seq":
                target_sequence = column
            elif field == "score":
                score = int(column)
            elif field == "identical":
                identical = int(column)
                annotations[field] = identical
            elif field == "positives":
                annotations[field] = int(column)
            elif field == "gaps":
                gaps = int(column)
                annotations[field] = gaps
            elif field == "% positives":
                annotations[field] = float(column)
            elif field == "query/sbjct frames":
                annotations[field] = column
            elif field == "query frame":
                annotations[field] = column
            elif field == "sbjct frame":
                annotations[field] = column
            else:
                raise ValueError("Unexpected field '%s'" % field)
        if alignment_length is not None and percentage_identity is not None:
            if identical is None:
                if mismatches is not None and gaps is not None:
                    identical = alignment_length - mismatches - gaps
            if identical is not None:
                difference = abs(100 * identical / alignment_length - percentage_identity)
                assert difference < 0.015
        coordinates[0, :] += target_start
        if query_start < query_end:
            coordinates[1, :] += query_start
        else:
            # mapped to reverse strand
            coordinates[1, :] = query_start - coordinates[1, :] + 1
        if query_sequence is None:
            query_seq = Seq(None, length=query_size)
        else:
            query_sequence = query_sequence.replace("-", "")
            assert len(query_sequence) == query_end - query_start
            query_seq = Seq({query_start: query_sequence}, length=query_size)
        query = SeqRecord(query_seq, id=query_id)
        if self._query_description is not None:
            query.description = self._query_description
        if target_sequence is None:
            if target_length is None:
                target_seq = Seq(None, length=target_end)
            else:
                target_seq = Seq(None, length=target_length)
        else:
            target_sequence = target_sequence.replace("-", "")
            if target_start is not None and target_end is not None:
                if self.program == "TBLASTN":
                    target_length = (target_end - target_start) // 3
                else:
                    target_length = target_end - target_start
                assert len(target_sequence) == target_length
            target_seq = Seq({target_start: target_sequence}, length=target_end)
        target = SeqRecord(target_seq, id=target_id)
        records = [target, query]
        alignment = Alignment(records, coordinates)
        alignment.annotations = annotations
        if score is not None:
            alignment.score = score
        return alignment

    def parse_btop(self, btop):
        """Parse a BTOP string and return alignment coordinates.

        A BTOP (Blast trace-back operations) string is used by BLAST to
        describe a sequence alignment.
        """
        target_coordinates = []
        query_coordinates = []
        target_coordinates.append(0)
        query_coordinates.append(0)
        state = State.NONE
        tokens = re.findall("([A-Z-]{2}|\\d+)", btop)
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
        """Parse a CIGAR string and return alignment coordinates.

        A CIGAR string, as defined by the SAM Sequence Alignment/Map format,
        describes a sequence alignment as a series of lengths and operation
        (alignment/insertion/deletion) codes.
        """
        target_coordinates = []
        query_coordinates = []
        target_coordinate = 0
        query_coordinate = 0
        target_coordinates.append(target_coordinate)
        query_coordinates.append(query_coordinate)
        state = State.NONE
        tokens = re.findall("(M|D|I|\\d+)", cigar)
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
            target_coordinates.append(target_coordinate)
            query_coordinates.append(query_coordinate)
        coordinates = numpy.array([target_coordinates, query_coordinates])
        return coordinates
