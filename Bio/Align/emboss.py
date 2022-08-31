# Copyright 2008-2016 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.Align support for "emboss" alignment output from EMBOSS tools.

This module contains a parser for the EMBOSS srspair/pair/simple file format,
for example from the needle, water, and stretcher tools.
"""
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

    def _read_header(self, stream):
        try:
            line = next(stream)
        except StopIteration:
            raise ValueError("Empty file.") from None
        if line.rstrip() != "########################################":
            raise ValueError("Unexpected line: %s") % line

        # assume srspair format (default) if not specified explicitly in
        # the output file
        self.metadata = {}
        self.metadata["Align_format"] = "srspair"
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
                self.metadata["Command line"] = commandline
                commandline = None
            key, value = line[2:].split(":", 1)
            if key == "Program":
                self.metadata["Program"] = value.strip()
            elif key == "Rundate":
                self.metadata["Rundate"] = value.strip()
            elif key == "Report_file":
                self.metadata["Report_file"] = value.strip()
            elif key == "Align_format":
                self.metadata["Align_format"] = value.strip()
            elif key == "Commandline":
                commandline = value.strip()

    def _read_next_alignment(self, stream):
        identifiers = None
        number_of_sequences = None
        annotations = {}
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
                    column = 0
                    index = 0
                    continue
                if line.strip() == "#":
                    continue
                if not line.startswith("# "):
                    raise ValueError("Unexpected line: %s") % line
                try:
                    key, value = line[2:].split(":", 1)
                except ValueError:
                    # An equal sign is used for Longest_Identity,
                    # Longest_Similarity, Shortest_Identity, and
                    # Shortest_Similarity, which are included if command line
                    # argument -nobrief was used.
                    key, value = line[2:].split(" = ", 1)
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
                    annotations[key] = value.strip()
                elif key == "Gap_penalty":
                    annotations[key] = float(value.strip())
                elif key == "Extend_penalty":
                    annotations[key] = float(value.strip())
                elif key == "Length":
                    ncols = int(value.strip())
                elif key == "Identity":
                    annotations[key] = int(value.strip().split("/")[0])
                elif key == "Similarity":
                    annotations[key] = int(value.strip().split("/")[0])
                elif key == "Gaps":
                    annotations[key] = int(value.strip().split("/")[0])
                elif key == "Score":
                    annotations[key] = float(value.strip())
                # TODO:
                # The following are generated if the -nobrief command line
                # argument used. We could simply calculate them from the
                # alignment, but then we have to define what we mean by
                # "similar". For now, simply store them as an annotation.
                elif key == "Longest_Identity":
                    annotations[key] = value.strip()
                elif key == "Longest_Similarity":
                    annotations[key] = value.strip()
                elif key == "Shortest_Identity":
                    annotations[key] = value.strip()
                elif key == "Shortest_Similarity":
                    annotations[key] = value.strip()
                else:
                    raise ValueError("Failed to parse line '%s'" % line)
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
                            records = []
                            n = len(sequences)
                            for i in range(n):
                                start = starts[i]
                                if start == 0:
                                    sequence = Seq(sequences[i])
                                else:
                                    coordinates[i, :] += start
                                    # create a partially defined sequence
                                    length = start + len(sequences[i])
                                    data = {start: sequences[i]}
                                    sequence = Seq(data, length=length)
                                record = SeqRecord(sequence, identifiers[i])
                                records.append(record)
                            alignment = Alignment(records, coordinates)
                            if annotations:
                                alignment.annotations = annotations
                            if consensus:
                                alignment.column_annotations = {
                                    "emboss_consensus": consensus
                                }
                            return alignment
                    continue
                prefix = line[:21].strip()
                if prefix == "":
                    # match line
                    consensus += line[21:71]
                else:
                    identifier, start = prefix.split(None, 1)
                    assert identifiers[index].startswith(identifier)
                    aligned_sequence, end = line[21:].split(None, 1)
                    start = int(start) - 1  # Python counting
                    end = int(end)
                    length = len(sequences[index])
                    sequence = aligned_sequence.replace("-", "")
                    if length == 0 and len(sequence) > 0:
                        # Record the start
                        starts[index] = start
                    else:
                        if (
                            self.metadata["Align_format"] == "srspair"
                            and len(sequence) == 0
                        ):
                            start += 1
                        assert start == starts[index] + length
                    assert end == start + len(sequence)
                    sequences[index] += sequence
                    aligned_sequences[index] += aligned_sequence
                    if index == 0:
                        column += len(aligned_sequence)
                    else:
                        assert column == len(aligned_sequences[index])
                    index += 1
