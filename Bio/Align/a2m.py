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

    fmt = "A2M"

    def format_alignment(self, alignment):
        """Return a string with the alignment in the A2M file format."""
        if not isinstance(alignment, Alignment):
            raise TypeError("Expected an Alignment object")
        lines = []
        state = alignment.column_annotations["state"]
        for sequence, line in zip(alignment.sequences, alignment):
            try:
                name = sequence.id
            except AttributeError:
                name = ""
            try:
                description = sequence.description
            except AttributeError:
                description = ""
            if description:
                lines.append(f">{name} {description}")
            else:
                lines.append(f">{name}")
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

    write_alignments = interfaces.AlignmentWriter.write_single_alignment


class AlignmentIterator(interfaces.AlignmentIterator):
    """Alignment iterator for files in the A2M file format.

    An A2M file contains one multiple alignment. Matches are represented by
    upper case letters and deletions by dashes in alignment columns containing
    matches or deletions only. Insertions are represented by lower case letters,
    with gaps aligned to the insertion shown as periods.  Header lines start
    with '>' followed by the name of the sequence, and optionally a description.
    """

    fmt = "A2M"

    def _read_next_alignment(self, stream):
        names = []
        descriptions = []
        lines = []
        for line in stream:
            if line.startswith(">"):
                parts = line[1:].rstrip().split(None, 1)
                try:
                    name = parts[0]
                except IndexError:
                    name = ""
                try:
                    description = parts[1]
                except IndexError:
                    description = ""
                names.append(name)
                descriptions.append(description)
                lines.append("")
            else:
                lines[-1] += line.strip()
        if not lines:
            if self._stream.tell() == 0:
                raise ValueError("Empty file.")
            return
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
            record = SeqRecord(sequence, name, description=description)
            records.append(record)
        alignment = Alignment(records, coordinates)
        alignment.column_annotations = {}
        alignment.column_annotations["state"] = state
        return alignment
