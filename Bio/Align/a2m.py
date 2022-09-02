# Copyright 2022 by Michiel de Hoon.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.Align support for A2M files.

A2M files are alignment files created by align2model or hmmscore in the SAM
Sequence Alignment and Modeling Software System.
"""
from Bio.Align import Alignment
from Bio.Align import interfaces
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class AlignmentWriter(interfaces.AlignmentWriter):
    """Alignment file writer for the A2M file format."""

    def __init__(self, target):
        """Create an AlignmentWriter object.

        Arguments:
         - target    - output stream or file name

        """
        super().__init__(target, mode="w")

    def format_alignment(self, alignment):
        """Return a string with the alignment in the A2M file format."""
        if not isinstance(alignment, Alignment):
            raise TypeError("Expected an Alignment object")
        lines = []
        state = alignment.column_annotations["state"]
        for sequence, line in zip(alignment.sequences, alignment):
            lines.append(f">{sequence.id} {sequence.description}")
            s = ""
            for c, m in zip(line, state):
                if m == "D":
                    s += c.upper()
                elif m == "I":
                    if c == "-":
                        s += "."
                    else:
                        s += c.lower()
            lines.append(s)
        return "\n".join(lines) + "\n"


class AlignmentIterator(interfaces.AlignmentIterator):
    """Alignment iterator for files in the A2M file format.

    An A2M file contains one multiple alignment. Matches are represented by
    upper case letters and deletions by dashes in alignment columns containing
    matches or deletions only. Insertions are represented by lower case letters,
    with gaps aligned to the insertion shown as periods.  Header lines start
    with '>' followed by the name of the sequence, and optionally a description.
    """

    def __init__(self, source):
        """Create an AlignmentIterator object.

        Arguments:
         - source   - input data or file name

        """
        super().__init__(source, mode="t", fmt="A2M")

    def _read_next_alignment(self, stream):
        names = []
        descriptions = []
        lines = []
        for line in stream:
            if line.startswith(">"):
                parts = line[1:].rstrip().split(None, 1)
                try:
                    name, description = parts
                except ValueError:
                    name = parts[0]
                    description = None
                names.append(name)
                descriptions.append(description)
                lines.append("")
            else:
                lines[-1] += line.strip()
        if not lines:
            raise ValueError("Empty file.")
        state = ""
        for c in lines[0]:
            if c == "-" or c.isupper():
                state += "D"  # Match/deletion state
            elif c == "." or c.islower():
                state += "I"  # Insertion state
            else:
                raise Exception("Unexpected letter '%s' in alignment" % c)
        for line in lines[1:]:
            for c, m in zip(line, state):
                if m == "D":  # Match/deletion state
                    assert c == "-" or c.isupper()
                elif m == "I":  # Insertion state
                    assert c == "." or c.islower()
                else:
                    raise Exception("Unexpected letter '%s' in alignment" % c)
        for i, line in enumerate(lines):
            lines[i] = line.upper().replace(".", "-")
        coordinates = Alignment.infer_coordinates(lines)
        records = []
        for name, description, line in zip(names, descriptions, lines):
            line = line.replace("-", "")
            sequence = Seq(line)
            if description is None:
                record = SeqRecord(sequence, name)
            else:
                record = SeqRecord(sequence, name, description=description)
            records.append(record)
        alignment = Alignment(records, coordinates)
        alignment.column_annotations = {}
        alignment.column_annotations["state"] = state
        self._close()  # a2m files contain only one alignment
        return alignment
