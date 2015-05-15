# Copyright 2015-2015 by Eric Rasche.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.
"""Bio.AlignIO support for "xmfa" output from Mauve/ProgressiveMauve.

You are expected to use this module via the Bio.AlignIO functions (or the
Bio.SeqIO functions if you want to work directly with the gapped sequences).

For example, consider a progressiveMauve alignment file containing the following::

    #FormatVersion Mauve1
    #Sequence1File	three.fa
    #Sequence1Entry	1
    #Sequence1Format	FastA
    #Sequence2File	three.fa
    #Sequence2Entry	2
    #Sequence2Format	FastA
    #Sequence3File	three.fa
    #Sequence3Entry	3
    #Sequence3Format	FastA
    #BackboneFile	three.xmfa.bbcols
    > 1:0-0 + three.fa
    --------------------------------------------------------------------------------
    --------------------------------------------------------------------------------
    --------------------------------------------------------------------------------
    > 2:5417-5968 + three.fa
    TTTAAACATCCCTCGGCCCGTCGCCCTTTTATAATAGCAGTACGTGAGAGGAGCGCCCTAAGCTTTGGGAAATTCAAGC-
    --------------------------------------------------------------------------------
    CTGGAACGTACTTGCTGGTTTCGCTACTATTTCAAACAAGTTAGAGGCCGTTACCTCGGGCGAACGTATAAACCATTCTG
    > 3:9476-10076 - three.fa
    TTTAAACACCTTTTTGGATG--GCCCAGTTCGTTCAGTTGTG-GGGAGGAGATCGCCCCAAACGTATGGTGAGTCGGGCG
    TTTCCTATAGCTATAGGACCAATCCACTTACCATACGCCCGGCGTCGCCCAGTCCGGTTCGGTACCCTCCATGACCCACG
    ---------------------------------------------------------AAATGAGGGCCCAGGGTATGCTT
    =
    > 2:5969-6015 + three.fa
    -----------------------
    GGGCGAACGTATAAACCATTCTG
    > 3:9429-9476 - three.fa
    TTCGGTACCCTCCATGACCCACG
    AAATGAGGGCCCAGGGTATGCTT

This is a multiple sequence alignment with multiple aligned sections, so you
would probably load this using the Bio.AlignIO.parse() function:

    >>> from Bio import AlignIO
    >>> align = AlignIO.parse("Mauve/simple.xmfa", "mauve")
    >>> alignments = list(align)
    >>> for aln in alignments:
    ...     print(align)
    SingleLetterAlphabet() alignment with 3 rows and 240 columns
    --------------------------------------------...--- 1
    TTTAAACATCCCTCGGCCCGTCGCCCTTTTATAATAGCAGTACG...CTG 2
    TTTAAACACCTTTTTGGATG--GCCCAGTTCGTTCAGTTGTG-G...CTT 3
    SingleLetterAlphabet() alignment with 3 rows and 46 columns
    ---------------------------------------------- 1
    -----------------------GGGCGAACGTATAAACCATTCTG 2
    TTCGGTACCCTCCATGACCCACGAAATGAGGGCCCAGGGTATGCTT 3

Additional information is extracted from the XMFA file and available through
the annotation attribute of each record::

    >>> for record in alignments[0]:
    ...   print record.id, len(record), record.annotations
    1 240 {'start': 0, 'end': 0, 'strand': 1}
    2 240 {'start': 5417, 'end': 5968, 'strand': 1}
    3 240 {'start': 9476, 'end': 10076, 'strand': -1}

"""

from __future__ import print_function

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from .Interfaces import AlignmentIterator
from .Interfaces import SequentialAlignmentWriter

__docformat__ = "restructuredtext en"


class MauveWriter(SequentialAlignmentWriter):
    """Mauve/XMFA alignment writer."""
    _wrote_header = False
    _wrote_first = False


    def write_alignment(self, alignment):
        """Use this to write (another) single alignment to an open file.

        Note that sequences and their annotation are recorded
        together (rather than having a block of annotation followed
        by a block of aligned sequences).
        """
        count = len(alignment)

        self._length_of_sequences = alignment.get_alignment_length()

        # NOTE - For now, the alignment object does not hold any per column
        # or per alignment annotation - only per sequence.

        if count == 0:
            raise ValueError("Must have at least one sequence")
        if self._length_of_sequences == 0:
            raise ValueError("Non-empty sequences are required")

        if not self._wrote_header:
            self._wrote_header = True
            self.handle.write("#FormatVersion Mauve1\n")

        for record in alignment:
            self._write_record(record)
        self.handle.write('=\n')

    def _write_record(self, record):
        """Write a single SeqRecord to the file"""
        if self._length_of_sequences != len(record.seq):
            raise ValueError("Sequences must all be the same length")

        seq_name = record.id

        if "start" in record.annotations \
                and "end" in record.annotations \
                and "strand" in record.annotations:
            id_line = "> %s:%s-%s %s unknown.fa\n" % (
                seq_name,
                record.annotations["start"] + 1,
                record.annotations["end"],
                "+" if record.annotations["strand"] == 1 else "-"
            )
        else:
            id_line = "> %s:0-0 + unknown.fa\n" % seq_name

        self.handle.write(id_line)
        for i in range(0, len(record.seq), 80):
            self.handle.write("%s\n" % str(record.seq[i:i + 80]))
        else:
            if not self._wrote_first:
                self._wrote_first = True
                # The first LCB we write out is special, and must list ALL
                # sequences, for the Mauve GUI
                # http://darlinglab.org/mauve/user-guide/files.html#non-standard-xmfa-formatting-used-by-the-mauve-gui
                id_line = "> %s:0-0 + \n\n" % seq_name
                self.handle.write(id_line)
            # Alignments lacking a start/stop/strand were generated by
            # BioPython on load, and shouldn't exist according to XMFA


class MauveIterator(AlignmentIterator):
    """Mauve xmfa alignment iterator."""

    _ids = []  # for caching IDs between __next__ calls

    def __next__(self):
        handle = self.handle
        line = handle.readline()

        if not line:
            raise StopIteration

        # Strip out header comments
        while line and line.strip().startswith('#'):
            line = handle.readline()

        seqs = {}
        seq_regions = {}
        passed_end_alignment = False

        latest_id = None
        while True:
            if not line:
                break  # end of file
            line = line.strip()

            if line.startswith('='):
                # There may be more data, but we've reached the end of this
                # alignment
                break
            elif line.startswith('>'):
                parts = line.split()
                id, start, end = self._identifier_split(parts[1])

                if id not in self._ids:
                    self._ids.append(id)

                seq_regions[id] = (start, end, parts[2])

                seqs.setdefault(id, '')
                latest_id = id
            else:
                assert not passed_end_alignment
                if latest_id is None:
                    raise ValueError("Saw sequence before definition line")
                seqs[latest_id] += line
            line = handle.readline()

        assert len(seqs) <= len(self._ids)

        self.ids = self._ids
        self.sequences = seqs

        if self._ids and seqs:
            alignment_length = len(list(seqs.values())[0])
            records = []
            for id in self._ids:
                if id not in seqs or len(seqs[id]) == 0:
                    seq = '-' * alignment_length
                else:
                    seq = seqs[id]
                if alignment_length != len(seq):
                    raise ValueError("Sequences have different lengths, or repeated identifier")
                record = SeqRecord(Seq(seq, self.alphabet), id=id, name=id,
                                   description=id)

                if id in seq_regions:
                    record.annotations["start"] = seq_regions[id][0]
                    record.annotations["end"] = seq_regions[id][1]
                    record.annotations["strand"] = 1 if seq_regions[id][2] == '+' else -1

                records.append(record)
            alignment = MultipleSeqAlignment(records, self.alphabet)
            return alignment
        else:
            raise StopIteration

    def _identifier_split(self, identifier):
        """Returns (name, start, end) string tuple from an identifier"""
        id, loc = identifier.split(':')
        start, end = map(int, loc.split('-'))
        start -= 1
        return id, start, end
