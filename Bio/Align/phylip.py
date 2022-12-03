# Copyright 2006-2016 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.Align support for the alignment format for input files for PHYLIP tools.

You are expected to use this module via the Bio.Align functions.
"""
from Bio.Align import Alignment
from Bio.Align import interfaces
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


_PHYLIP_ID_WIDTH = 10


class AlignmentWriter(interfaces.AlignmentWriter):
    """Clustalw alignment writer."""

    fmt = "PHYLIP"

    def format_alignment(self, alignment):
        """Return a string with a single alignment in the Phylip format."""
        names = []
        for record in alignment.sequences:
            try:
                name = record.id
            except AttributeError:
                name = ""
            else:
                name = name.strip()
                for char in "[](),":
                    name = name.replace(char, "")
                for char in ":;":
                    name = name.replace(char, "|")
                name = name[:_PHYLIP_ID_WIDTH]
            names.append(name)

        lines = []
        nseqs, length = alignment.shape
        if nseqs == 0:
            raise ValueError("Must have at least one sequence")
        if length == 0:
            raise ValueError("Non-empty sequences are required")
        line = "%d %d\n" % (nseqs, length)
        lines.append(line)

        # From experimentation, the use of tabs is not understood by the
        # EMBOSS suite.  The nature of the expected white space is not
        # defined in the PHYLIP documentation, simply "These are in free
        # format, separated by blanks".  We'll use spaces to keep EMBOSS
        # happy.
        for name, sequence in zip(names, alignment):
            # Write the entire sequence to one line
            line = name[:_PHYLIP_ID_WIDTH].ljust(_PHYLIP_ID_WIDTH) + sequence + "\n"
            lines.append(line)
        return "".join(lines)


class AlignmentIterator(interfaces.AlignmentIterator):
    """Reads a Phylip alignment file and returns an Alignment iterator.

    Record names are limited to at most 10 characters.

    The parser determines from the file contents if the file format is
    sequential or interleaved, and parses the file accordingly.

    For more information on the file format, please see:
    http://evolution.genetics.washington.edu/phylip/doc/sequence.html
    http://evolution.genetics.washington.edu/phylip/doc/main.html#inputfiles
    """

    fmt = "PHYLIP"

    def _read_header(self, stream):
        try:
            line = next(stream)
        except StopIteration:
            raise ValueError("Empty file.") from None

        words = line.split()
        if len(words) == 2:
            try:
                self._number_of_seqs = int(words[0])
                self._length_of_seqs = int(words[1])
                return
            except ValueError:
                pass
        raise ValueError(
            "Expected two integers in the first line, received '%s'" % line
        )

    def _parse_interleaved_first_block(self, lines, seqs, names):
        for line in lines:
            line = line.rstrip()
            name = line[:_PHYLIP_ID_WIDTH].strip()
            seq = line[_PHYLIP_ID_WIDTH:].strip().replace(" ", "")
            names.append(name)
            seqs.append([seq])

    def _parse_interleaved_other_blocks(self, stream, seqs):
        i = 0
        for line in stream:
            line = line.rstrip()
            if not line:
                assert i == self._number_of_seqs
                i = 0
            else:
                seq = line.replace(" ", "")
                seqs[i].append(seq)
                i += 1
        if i != 0 and i != self._number_of_seqs:
            raise ValueError("Unexpected file format")

    def _parse_sequential(self, lines, seqs, names, length):
        for line in lines:
            if length == 0:
                line = line.rstrip()
                name = line[:_PHYLIP_ID_WIDTH].strip()
                seq = line[_PHYLIP_ID_WIDTH:].strip()
                names.append(name)
                seqs.append([])
            else:
                seq = line.strip()
            seq = seq.replace(" ", "")
            seqs[-1].append(seq)
            length += len(seq)
            if length == self._length_of_seqs:
                length = 0
        return length

    def _read_file(self, stream):
        names = []
        seqs = []
        lines = [next(stream) for i in range(self._number_of_seqs)]
        try:
            line = next(stream)
        except StopIteration:
            pass
        else:
            if line.rstrip():
                # sequential file format
                lines.append(line)
                length = self._parse_sequential(lines, seqs, names, 0)
                self._parse_sequential(stream, seqs, names, length)
                return names, seqs
        # interleaved file format
        self._parse_interleaved_first_block(lines, seqs, names)
        self._parse_interleaved_other_blocks(stream, seqs)
        return names, seqs

    def _read_next_alignment(self, stream):
        names, seqs = self._read_file(stream)

        seqs = ["".join(seq) for seq in seqs]
        if len(seqs) != self._number_of_seqs:
            raise ValueError(
                "Found %i records in this alignment, told to expect %i"
                % (len(seqs), self._number_of_seqs)
            )
        for seq in seqs:
            if len(seq) != self._length_of_seqs:
                raise ValueError(
                    "Expected all sequences to have length %d; found %d"
                    % (self._length_of_seqs, len(seq))
                )
            if "." in seq:
                raise ValueError("PHYLIP format no longer allows dots in sequence")

        coordinates = Alignment.infer_coordinates(seqs)
        seqs = [seq.replace("-", "") for seq in seqs]
        records = [
            SeqRecord(Seq(seq), id=name, description="")
            for (name, seq) in zip(names, seqs)
        ]
        alignment = Alignment(records, coordinates)
        del self._number_of_seqs
        del self._length_of_seqs
        self._close()
        return alignment
