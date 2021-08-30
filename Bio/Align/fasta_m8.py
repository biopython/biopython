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
import numpy
from Bio.Align import Alignment
from Bio.Align import interfaces
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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
                query_line, query_size = line[len(prefix):].strip().rsplit(" - ", 1)
                query_size, unit = query_size.split()
                self._query_size = int(query_size)
                assert unit in ("nt", "aa")
                self._query_id, self._query_description = query_line.split(None, 1)
                line = next(stream)
                prefix = "# Database: "
                assert line.startswith(prefix)
                database = line[len(prefix):].strip()
                line = next(stream)
                prefix = "# Fields: "
                assert line.startswith(prefix)
                fields = line[len(prefix):].strip().split(", ")
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
                hits = int(line[2:-len(suffix)])
            else:
                yield self.create_alignment(line)

    def create_alignment(self, line):
        columns = line.split()
        assert len(columns) == 13
        annotations = {}
        if self._program is not None:
            annotations['program'] = self._program
        if self._query_id is not None:
            assert columns[0] == self._query_id
        query_id = columns[0]
        target_id = columns[1]
        percentage_identity = float(columns[2])
        alignment_length = int(columns[3])
        mismatches = int(columns[4])
        gap_opens = int(columns[5])
        query_start = int(columns[6]) - 1
        query_end = int(columns[7])
        target_start = int(columns[8]) - 1
        target_end = int(columns[9])
        annotations['evalue'] = float(columns[10])
        annotations['bit_score'] = float(columns[11])
        if self._alignment_representation == "BTOP":
            coordinates = self.parse_btop(columns[12])
        elif self._alignment_representation == "CIGAR":
            coordinates = self.parse_cigar(columns[12])
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
        target_position = 0
        query_position = 0
        target_coordinates.append(target_position)
        query_coordinates.append(query_position)
        digits = ""
        aligned = 0
        letters = iter(btop)
        state = None
        for letter in letters:
            if letter in '0123456789':
                state = "match"
                digits += letter
            else:
                query_letter = letter
                target_letter = next(letters)
                if query_letter == "-":
                    if state == "match":
                        number = aligned
                        if digits:
                            number += int(digits)
                        target_position += number
                        query_position += number
                        target_coordinates.append(target_position)
                        query_coordinates.append(query_position)
                        digits = ""
                        aligned = 0
                    elif state == "target_gap":
                        target_coordinates.append(target_position)
                        query_coordinates.append(query_position)
                    state = "query_gap"
                    target_position += 1
                elif target_letter == "-":
                    if state == "match":
                        number = aligned
                        if digits:
                            number += int(digits)
                        target_position += number
                        query_position += number
                        target_coordinates.append(target_position)
                        query_coordinates.append(query_position)
                        digits = ""
                        aligned = 0
                    elif state == "query_gap":
                        target_coordinates.append(target_position)
                        query_coordinates.append(query_position)
                    state = "target_gap"
                    query_position += 1
                else:
                    # mismatch; stay in match state
                    if state == "match":
                        if digits:
                            aligned += int(digits)
                            digits = ""
                    else:
                        query_coordinates.append(query_position)
                        target_coordinates.append(target_position)
                        state = "match"
                        aligned = 0
                    aligned += 1
        target_coordinates.append(target_position)
        query_coordinates.append(query_position)
        return numpy.array([target_coordinates, query_coordinates])

    def parse_cigar(self, cigar):
        return None
