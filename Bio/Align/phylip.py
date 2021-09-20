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

    def write_alignment(self, alignment):
        """Use this to write (another) single alignment to an open file."""
        names = []
        for record in alignment.sequences:
            name = record.id.strip()
            for char in "[](),":
                name = name.replace(char, "")
            for char in ":;":
                name = name.replace(char, "|")
            name = name[:_PHYLIP_ID_WIDTH]
            names.append(name)
            if "." in record.seq:
                raise ValueError("PHYLIP format no longer allows dots in sequence")

        stream = self.stream

        nseqs, length = alignment.shape
        if nseqs == 0:
            raise ValueError("Must have at least one sequence")
        if length == 0:
            raise ValueError("Non-empty sequences are required")
        line = "%d %d\n" % (nseqs, length)
        stream.write(line)

        # From experimentation, the use of tabs is not understood by the
        # EMBOSS suite.  The nature of the expected white space is not
        # defined in the PHYLIP documentation, simply "These are in free
        # format, separated by blanks".  We'll use spaces to keep EMBOSS
        # happy.
        for name, sequence in zip(names, alignment):
            stream.write(name[:_PHYLIP_ID_WIDTH].ljust(_PHYLIP_ID_WIDTH))
            # Write the entire sequence to one line (see sequential format
            # notes in the AlignmentIterator docstring)
            stream.write(sequence)
            stream.write("\n")


class AlignmentIterator(interfaces.AlignmentIterator):
    """Reads a Phylip alignment file and returns an Alignment iterator.

    Record names are limited to at most 10 characters.

    It only copes with interleaved phylip files!  Sequential files won't work
    where the sequences are split over multiple lines.

    For more information on the file format, please see:
    http://evolution.genetics.washington.edu/phylip/doc/sequence.html
    http://evolution.genetics.washington.edu/phylip/doc/main.html#inputfiles
    """

    def __init__(self, source):
        """Create an AlignmentIterator object.

        Arguments:
         - source   - input data or file name

        """
        super().__init__(source, mode="t", fmt="PHYLIP")
        stream = self.stream
        try:
            line = next(stream)
        except StopIteration:
            raise ValueError("Empty file.") from None

        words = line.split()
        if len(words) == 2:
            try:
                self.number_of_seqs = int(words[0])
                self.length_of_seqs = int(words[1])
                return
            except ValueError:
                pass
        raise ValueError("Expected two integers in the first line, received '%s'" % line)

    def parse(self, stream):
        """Parse the next alignment from the stream."""

        names = []
        seqs = []
        i = 0
        for line in stream:
            line = line.rstrip()
            if not line:
                break
            name = line[: _PHYLIP_ID_WIDTH].strip()
            seq = line[_PHYLIP_ID_WIDTH :].strip().replace(" ", "")
            names.append(name)
            if "." in seq:
                raise ValueError("PHYLIP format no longer allows dots in sequence")
            seqs.append([seq])
            i += 1

        if i != self.number_of_seqs:
            raise ValueError("Unexpected file format")
        i = 0
        for line in stream:
            line = line.rstrip()
            if not line:
                assert i == self.number_of_seqs
                i = 0
            else:
                seq = line.replace(" ", "")
                seqs[i].append(seq)
                i += 1
        if i != 0 and i != self.number_of_seqs:
            raise ValueError("Unexpected file format")

        seqs = ["".join(seq) for seq in seqs]
        if len(seqs) != self.number_of_seqs:
            raise ValueError(
                "Found %i records in this alignment, told to expect %i"
                % (len(seqs), number_of_seqs)
            )
        for seq in seqs:
            if len(seq) != self.length_of_seqs:
                raise ValueError("Expected all sequences to have length %d; found %d" % (length_of_seqs, len(seq)))
            if "." in seq:
                raise ValueError("PHYLIP format no longer allows dots in sequence")

        coordinates = Alignment.infer_coordinates(seqs)
        seqs = [seq.replace("-", "") for seq in seqs]
        records = [SeqRecord(Seq(seq), id=name) for (name, seq) in zip(names, seqs)]
        alignment = Alignment(records, coordinates)
        yield alignment
