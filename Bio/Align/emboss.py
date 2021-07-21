# Copyright 2008-2016 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.AlignIO support for "emboss" alignment output from EMBOSS tools.

You are expected to use this module via the Bio.Align functions (or the
Bio.SeqIO functions if you are interested in the sequences only).

This module contains a parser for the EMBOSS pairs/simple file format, for
example from the alignret, water and needle tools.
"""
import Bio
from Bio.Align import Alignment
from Bio.Align import interfaces
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class AlignmentIterator(interfaces.AlignmentIterator):
    """Emboss alignment iterator.

    For reading the (pairwise) alignments from EMBOSS tools in what they
    call the "pairs" and "simple" formats.
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
        if line.rstrip() != "########################################":
            raise ValueError("Unexpected line: %s") % line

        commandline = None
        for line in stream:
            if line.rstrip() == "########################################":
                break
            if not line.startswith("# "):
                raise ValueError("Unexpected line: %s") % line
            if commandline is not None:
                if line.startswith("#    "):
                    commandline += " " + line[1:].strip()
                    continue
                self.commandline = commandline
                commandline = None
            key, value = line[2:].split(":", 1)
            if key == "Program":
                self.program = value.strip()
            elif key == "Rundate":
                self.rundate = value.strip()
            elif key == "Report_file":
                self.report_file = value.strip()
            elif key == "Align_format":
                self.align_format = value.strip()
            elif key == "Commandline":
                commandline = value.strip()

    def parse(self, stream):
        """Parse the next alignment from the stream."""
        if stream is None:
            raise StopIteration

        identifiers = None
        number_of_sequences = None
        for line in stream:
            line = line.rstrip("\r\n")
            if identifiers is None:
                # searching for alignment metadata start
                if not line:
                    continue
                elif line.startswith("#---------------------------------------"):
                    # may appear between alignments
                    continue
                elif line.startswith("#======================================="):
                    # found the alignment metadata start
                    identifiers = []
                    ncols = None
                    sequences = None
                    matrix = None
                    gap_penalty = None
                    extend_penalty = None
                    identity = None
                    similarity = None
                    gaps = None
                    score = None
                else:
                    raise ValueError("Unexpected line: %s" % line)
            elif sequences is None:
                # parsing the alignment metadata
                if line == "#=======================================":
                    # reached the end of alignment metadata
                    if len(identifiers) == 0:
                        raise ValueError("Number of sequences missing!")
                    if ncols is None:
                        raise ValueError("Length of alignment missing!")
                    sequences = [""] * number_of_sequences
                    aligned_sequences = [""] * number_of_sequences
                    consensus = ""
                    starts = [0] * number_of_sequences
                    ends = [0] * number_of_sequences
                    column = 0
                    index = 0
                    continue
                if line.strip() == "#":
                    continue
                if not line.startswith("# "):
                    raise ValueError("Unexpected line: %s") % line
                key, value = line[2:].split(":", 1)
                if key == "Aligned_sequences":
                    number_of_sequences = int(value.strip())
                    assert len(identifiers) == 0
                    # Should now expect the record identifiers...
                    for i, line in enumerate(stream):
                        if not line.startswith("# "):
                            raise ValueError("Unexpected line: %s") % line
                        number, identifier = line[2:].split(":")
                        assert i + 1 == int(number)
                        identifiers.append(identifier.strip())
                        if len(identifiers) == number_of_sequences:
                            break
                elif key == "Matrix":
                    matrix = value.strip()
                elif key == "Gap_penalty":
                    gap_penalty = float(value.strip())
                elif key == "Extend_penalty":
                    extend_penalty = float(value.strip())
                elif key == "Length":
                    ncols = int(value.strip())
                elif key == "Identity":
                    identity = int(value.strip().split("/")[0])
                elif key == "Similarity":
                    similarity = int(value.strip().split("/")[0])
                elif key == "Gaps":
                    gaps = int(value.strip().split("/")[0])
                elif key == "Score":
                    score = float(value.strip())
            else:
                # parse the sequences
                if not line:
                    # empty line
                    if index == number_of_sequences:
                        # reached the end of an alignment block
                        index = 0
                        if column == ncols:
                            # reached the end of the sequences
                            coordinates = Alignment.infer_coordinates(aligned_sequences)
                            for i, start in enumerate(starts):
                                start -= 1  # Python counting
                                coordinates[i, :] += start
                            sequences = [Seq(sequence) for sequence in sequences]
                            records = [
                                SeqRecord(sequence, id=identifier)
                                for sequence, identifier in zip(sequences, identifiers)
                            ]
                            alignment = Alignment(records, coordinates)
                            if matrix is not None:
                                alignment.matrix = matrix
                            if gap_penalty is not None:
                                alignment.gap_penalty = gap_penalty
                            if extend_penalty is not None:
                                alignment.extend_penalty = extend_penalty
                            if identity is not None:
                                alignment.identity = identity
                            if similarity is not None:
                                alignment.similarity = similarity
                            if gaps is not None:
                                alignment.gaps = gaps
                            if score is not None:
                                alignment.score = score
                            if consensus:
                                alignment.column_annotations = {
                                    "emboss_consensus": consensus
                                }
                            yield alignment
                            identifiers = None
                    continue
                prefix = line[:21].strip()
                if prefix == "":
                    # match line
                    consensus += line[21:71]
                else:
                    identifier, start = prefix.split(None, 1)
                    aligned_sequence, end = line[21:].split(None, 1)
                    start = int(start)
                    end = int(end)
                    sequence = aligned_sequence.replace("-", "")
                    if len(sequences[index]) > 0:
                        length = len(sequence)
                        if length == 0:
                            assert start == ends[index]
                            assert end == ends[index]
                        else:
                            assert start == ends[index] + 1
                            assert end == ends[index] + length
                    assert identifiers[index].startswith(identifier)
                    if starts[index] == 0:
                        # Record the start and end
                        starts[index] = start
                    ends[index] = end
                    sequences[index] += sequence
                    aligned_sequences[index] += aligned_sequence
                    if index == 0:
                        column += len(aligned_sequence)
                    else:
                        assert column == len(aligned_sequences[index])
                    index += 1
