# Copyright 2006-2014 by Peter Cock.  All rights reserved.
# Revisions copyright 2011 Brandon Invergo. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.Align.mauve module."""
import os
import unittest
from io import StringIO

from Bio import Align
from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.Seq import Seq

try:
    import numpy as np
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align.mauve."
    ) from None


class TestCombinedFile(unittest.TestCase):
    # Generate the output file combined.xmfa by running
    # progressiveMauve combined.fa --output=combined.xmfa

    filename = "combined.fa"
    path = os.path.join("Mauve", filename)
    records = SeqIO.parse(path, "fasta")
    sequences = {str(index): record.seq for index, record in enumerate(records)}
    del filename
    del path
    del records

    def test_parse(self):
        path = os.path.join("Mauve", "combined.xmfa")
        with open(path) as stream:
            alignments = Align.parse(stream, "mauve")
            self.check_alignments(alignments)
            alignments = iter(alignments)
            self.check_alignments(alignments)
        with Align.parse(path, "mauve") as alignments:
            self.check_alignments(alignments)
        with self.assertRaises(AttributeError):
            alignments._stream
        with Align.parse(path, "mauve") as alignments:
            pass
        with self.assertRaises(AttributeError):
            alignments._stream

    def check_alignments(self, alignments):
        saved_alignments = []
        metadata = alignments.metadata
        self.assertEqual(len(metadata), 3)
        self.assertEqual(metadata["FormatVersion"], "Mauve1")
        self.assertEqual(metadata["File"], "combined.fa")
        self.assertEqual(metadata["BackboneFile"], "combined.xmfa.bbcols")
        identifiers = alignments.identifiers
        alignment = next(alignments)
        saved_alignments.append(alignment)
        self.assertEqual(len(alignment), 3)
        self.assertEqual(len(alignment.sequences), 3)
        self.assertEqual(alignment.sequences[0].id, "0")
        self.assertEqual(
            repr(alignment.sequences[0].seq),
            "Seq({1: 'AAAAGGAAAGTACGGCCCGGCCACTCCGGGTGTGTGCTAGGAGGGCTT'}, length=49)",
        )
        sequence = self.sequences[alignment.sequences[0].id]
        start = len(sequence) - alignment.coordinates[0, 0]
        end = len(sequence) - alignment.coordinates[0, -1]
        self.assertEqual(start, 1)
        self.assertEqual(end, 49)
        self.assertEqual(alignment.sequences[0].seq[start:end], sequence[start:end])
        self.assertEqual(alignment.sequences[1].id, "1")
        self.assertEqual(alignment.sequences[1].seq, "")
        start = alignment.coordinates[1, 0]
        end = alignment.coordinates[1, -1]
        self.assertEqual(start, 0)
        self.assertEqual(end, 0)
        self.assertEqual(alignment.sequences[2].id, "2")
        self.assertEqual(
            repr(alignment.sequences[2].seq),
            "Seq({1: 'AAGCCCTGCGCGCTCAGCCGGAGTGTCCCGGGCCCTGCTTTCCTTTT'}, length=48)",
        )
        start = alignment.coordinates[2, 0]
        end = alignment.coordinates[2, -1]
        self.assertEqual(start, 1)
        self.assertEqual(end, 48)
        sequence = self.sequences[alignment.sequences[2].id][start:end]
        self.assertEqual(alignment.sequences[2].seq[start:end], sequence)
        self.assertEqual(
            alignment[0], "AAGCCCTCCTAGCACACACCCGGAGTGG-CCGGGCCGTACTTTCCTTTT"
        )
        self.assertEqual(
            alignment[1], "-------------------------------------------------"
        )
        self.assertEqual(
            alignment[2], "AAGCCCTGC--GCGCTCAGCCGGAGTGTCCCGGGCCCTGCTTTCCTTTT"
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                np.array(
                    [
                        [49, 40, 38, 21, 21, 1],
                        [0, 0, 0, 0, 0, 0],
                        [1, 10, 10, 27, 28, 48],
                    ]
                ),
            )
        )
        self.assertEqual(
            str(alignment),
            """\
0                49 AAGCCCTCCTAGCACACACCCGGAGTGG-CCGGGCCGTACTTTCCTTTT  1
1                 0 -------------------------------------------------  0
2                 1 AAGCCCTGC--GCGCTCAGCCGGAGTGTCCCGGGCCCTGCTTTCCTTTT 48
""",
        )
        self.assertEqual(
            alignment.format("mauve", metadata, identifiers),
            """\
> 1:2-49 - combined.fa
AAGCCCTCCTAGCACACACCCGGAGTGG-CCGGGCCGTACTTTCCTTTT
> 2:0-0 + combined.fa
-------------------------------------------------
> 3:2-48 + combined.fa
AAGCCCTGC--GCGCTCAGCCGGAGTGTCCCGGGCCCTGCTTTCCTTTT
=
""",
        )
        self.assertTrue(
            np.array_equal(
                np.array(alignment, "U"),
                # fmt: off
np.array([['A', 'A', 'G', 'C', 'C', 'C', 'T', 'C', 'C', 'T', 'A', 'G', 'C',
           'A', 'C', 'A', 'C', 'A', 'C', 'C', 'C', 'G', 'G', 'A', 'G', 'T',
           'G', 'G', '-', 'C', 'C', 'G', 'G', 'G', 'C', 'C', 'G', 'T', 'A',
           'C', 'T', 'T', 'T', 'C', 'C', 'T', 'T', 'T', 'T'],
          ['-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
           '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
           '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
           '-', '-', '-', '-', '-', '-', '-', '-', '-', '-'],
          ['A', 'A', 'G', 'C', 'C', 'C', 'T', 'G', 'C', '-', '-', 'G', 'C',
           'G', 'C', 'T', 'C', 'A', 'G', 'C', 'C', 'G', 'G', 'A', 'G', 'T',
           'G', 'T', 'C', 'C', 'C', 'G', 'G', 'G', 'C', 'C', 'C', 'T', 'G',
           'C', 'T', 'T', 'T', 'C', 'C', 'T', 'T', 'T', 'T']], dtype='U')
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(counts.left_insertions, 47)
        self.assertEqual(counts.left_deletions, 48)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 1)
        self.assertEqual(counts.internal_deletions, 2)
        self.assertEqual(counts.left_gaps, 95)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 3)
        self.assertEqual(counts.insertions, 48)
        self.assertEqual(counts.deletions, 50)
        self.assertEqual(counts.gaps, 98)
        self.assertEqual(counts.aligned, 46)
        self.assertEqual(counts.identities, 39)
        self.assertEqual(counts.mismatches, 7)
        alignment = next(alignments)
        saved_alignments.append(alignment)
        self.assertEqual(len(alignment), 1)
        self.assertEqual(len(alignment.sequences), 1)
        self.assertEqual(alignment.sequences[0].id, "0")
        self.assertEqual(alignment.sequences[0].seq, "G")
        sequence = self.sequences[alignment.sequences[0].id]
        start = alignment.coordinates[0, 0]
        end = alignment.coordinates[0, -1]
        self.assertEqual(alignment.sequences[0].seq[start:end], sequence[start:end])
        self.assertEqual(alignment[0], "G")
        self.assertTrue(np.array_equal(alignment.coordinates, np.array([[0, 1]])))
        self.assertEqual(
            str(alignment),
            """\
0                 0 G 1
""",
        )
        self.assertEqual(
            alignment.format("mauve", metadata, identifiers),
            """\
> 1:1-1 + combined.fa
G
=
""",
        )
        self.assertTrue(
            np.array_equal(np.array(alignment, "U"), np.array([["G"]], dtype="U"))
        )
        counts = alignment.counts()
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 0)
        self.assertEqual(counts.identities, 0)
        self.assertEqual(counts.mismatches, 0)
        alignment = next(alignments)
        saved_alignments.append(alignment)
        self.assertEqual(len(alignment), 1)
        self.assertEqual(len(alignment.sequences), 1)
        self.assertEqual(alignment.sequences[0].id, "0")
        self.assertEqual(repr(alignment.sequences[0].seq), "Seq({49: 'A'}, length=50)")
        sequence = self.sequences[alignment.sequences[0].id]
        start = alignment.coordinates[0, 0]
        end = alignment.coordinates[0, -1]
        self.assertEqual(alignment.sequences[0].seq[start:end], sequence[start:end])
        self.assertEqual(alignment[0], "A")
        self.assertTrue(np.array_equal(alignment.coordinates, np.array([[49, 50]])))
        self.assertEqual(
            str(alignment),
            """\
0                49 A 50
""",
        )
        self.assertEqual(
            alignment.format("mauve", metadata, identifiers),
            """\
> 1:50-50 + combined.fa
A
=
""",
        )
        self.assertTrue(
            np.array_equal(np.array(alignment, "U"), np.array([["A"]], dtype="U"))
        )
        counts = alignment.counts()
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 0)
        self.assertEqual(counts.identities, 0)
        self.assertEqual(counts.mismatches, 0)
        alignment = next(alignments)
        saved_alignments.append(alignment)
        self.assertEqual(len(alignment), 1)
        self.assertEqual(len(alignment.sequences), 1)
        self.assertEqual(alignment.sequences[0].id, "1")
        self.assertEqual(
            alignment.sequences[0].seq, "GAAGAGGAAAAGTAGATCCCTGGCGTCCGGAGCTGGGACGT"
        )
        sequence = self.sequences[alignment.sequences[0].id]
        start = alignment.coordinates[0, 0]
        end = alignment.coordinates[0, -1]
        self.assertEqual(alignment.sequences[0].seq[start:end], sequence[start:end])
        self.assertEqual(alignment[0], "GAAGAGGAAAAGTAGATCCCTGGCGTCCGGAGCTGGGACGT")
        self.assertTrue(np.array_equal(alignment.coordinates, np.array([[0, 41]])))
        self.assertEqual(
            str(alignment),
            """\
1                 0 GAAGAGGAAAAGTAGATCCCTGGCGTCCGGAGCTGGGACGT 41
""",
        )
        self.assertEqual(
            alignment.format("mauve", metadata, identifiers),
            """\
> 2:1-41 + combined.fa
GAAGAGGAAAAGTAGATCCCTGGCGTCCGGAGCTGGGACGT
=
""",
        )
        self.assertTrue(
            np.array_equal(
                np.array(alignment, "U"),
                # fmt: off
np.array([['G', 'A', 'A', 'G', 'A', 'G', 'G', 'A', 'A', 'A', 'A', 'G', 'T',
           'A', 'G', 'A', 'T', 'C', 'C', 'C', 'T', 'G', 'G', 'C', 'G', 'T',
           'C', 'C', 'G', 'G', 'A', 'G', 'C', 'T', 'G', 'G', 'G', 'A', 'C',
           'G', 'T']], dtype='U')
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 0)
        self.assertEqual(counts.identities, 0)
        self.assertEqual(counts.mismatches, 0)
        alignment = next(alignments)
        saved_alignments.append(alignment)
        self.assertEqual(len(alignment), 1)
        self.assertEqual(len(alignment.sequences), 1)
        self.assertEqual(alignment.sequences[0].id, "2")
        self.assertEqual(alignment.sequences[0].seq, "C")
        sequence = self.sequences[alignment.sequences[0].id]
        start = alignment.coordinates[0, 0]
        end = alignment.coordinates[0, -1]
        self.assertEqual(alignment.sequences[0].seq[start:end], sequence[start:end])
        self.assertEqual(alignment[0], "C")
        self.assertTrue(np.array_equal(alignment.coordinates, np.array([[0, 1]])))
        self.assertEqual(
            str(alignment),
            """\
2                 0 C 1
""",
        )
        self.assertEqual(
            alignment.format("mauve", metadata, identifiers),
            """\
> 3:1-1 + combined.fa
C
=
""",
        )
        self.assertTrue(
            np.array_equal(np.array(alignment, "U"), np.array([["C"]], dtype="U"))
        )
        counts = alignment.counts()
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 0)
        self.assertEqual(counts.identities, 0)
        self.assertEqual(counts.mismatches, 0)
        alignment = next(alignments)
        saved_alignments.append(alignment)
        self.assertEqual(len(alignment), 1)
        self.assertEqual(len(alignment.sequences), 1)
        self.assertEqual(alignment.sequences[0].id, "2")
        self.assertEqual(repr(alignment.sequences[0].seq), "Seq({48: 'C'}, length=49)")
        sequence = self.sequences[alignment.sequences[0].id]
        start = alignment.coordinates[0, 0]
        end = alignment.coordinates[0, -1]
        self.assertEqual(alignment.sequences[0].seq[start:end], sequence[start:end])
        self.assertEqual(alignment[0], "C")
        self.assertTrue(np.array_equal(alignment.coordinates, np.array([[48, 49]])))
        self.assertEqual(
            str(alignment),
            """\
2                48 C 49
""",
        )
        self.assertEqual(
            alignment.format("mauve", metadata, identifiers),
            """\
> 3:49-49 + combined.fa
C
=
""",
        )
        self.assertTrue(
            np.array_equal(np.array(alignment, "U"), np.array([["C"]], dtype="U"))
        )
        counts = alignment.counts()
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 0)
        self.assertEqual(counts.identities, 0)
        self.assertEqual(counts.mismatches, 0)
        self.assertRaises(StopIteration, next, alignments)
        # As each nucleotide in each sequence is stored exactly once in an XMFA
        # file, we can reconstitute the full sequences:
        self.assertEqual(len(saved_alignments), 6)
        maxindex = -1
        for alignment in saved_alignments:
            for record in alignment.sequences:
                index = int(record.id)
                if index > maxindex:
                    maxindex = index
        n = maxindex + 1
        self.assertEqual(n, 3)
        lengths = [0] * n
        for alignment in saved_alignments:
            for record in alignment.sequences:
                index = int(record.id)
                length = len(record.seq)
                if length > lengths[index]:
                    lengths[index] = length
        self.assertEqual(lengths[0], 50)
        self.assertEqual(lengths[1], 41)
        self.assertEqual(lengths[2], 49)
        sequences = [None] * 3
        for index, length in enumerate(lengths):
            sequences[index] = MutableSeq("N" * length)
        # Now fill up the sequences:
        for alignment in saved_alignments:
            for row, record in zip(alignment.coordinates, alignment.sequences):
                index = int(record.id)
                start = row[0]
                end = row[-1]
                if start > end:
                    start, end = end, start
                sequences[index][start:end] = record.seq[start:end]
        # Confirm that the fully defined sequences agree with the Fasta file:
        for index, sequence in enumerate(sequences):
            sequences[index] = Seq(sequence)
            key = str(index)
            self.assertEqual(sequences[index], self.sequences[key])
        # Make sure we can replace the partially defined sequences by these
        # fully defined sequences, and get the same alignment:
        alignment = saved_alignments[0]
        for record in alignment.sequences:
            index = int(record.id)
            record.seq = sequences[index]
            self.assertEqual(
                alignment[0], "AAGCCCTCCTAGCACACACCCGGAGTGG-CCGGGCCGTACTTTCCTTTT"
            )
            self.assertEqual(
                alignment[1], "-------------------------------------------------"
            )
            self.assertEqual(
                alignment[2], "AAGCCCTGC--GCGCTCAGCCGGAGTGTCCCGGGCCCTGCTTTCCTTTT"
            )
        alignment = saved_alignments[1]
        for record in alignment.sequences:
            index = int(record.id)
            record.seq = sequences[index]
            self.assertEqual(alignment[0], "G")
        alignment = saved_alignments[2]
        for record in alignment.sequences:
            index = int(record.id)
            record.seq = sequences[index]
            self.assertEqual(alignment[0], "A")
        alignment = saved_alignments[3]
        for record in alignment.sequences:
            index = int(record.id)
            record.seq = sequences[index]
            self.assertEqual(alignment[0], "GAAGAGGAAAAGTAGATCCCTGGCGTCCGGAGCTGGGACGT")
        alignment = saved_alignments[4]
        for record in alignment.sequences:
            index = int(record.id)
            record.seq = sequences[index]
            self.assertEqual(alignment[0], "C")
        alignment = saved_alignments[5]
        for record in alignment.sequences:
            index = int(record.id)
            record.seq = sequences[index]
            self.assertEqual(alignment[0], "C")

    def test_write_read(self):
        path = os.path.join("Mauve", "combined.xmfa")
        with open(path) as stream:
            data = stream.read()

        stream = StringIO()
        stream.write(data)
        stream.seek(0)
        alignments = Align.parse(stream, "mauve")
        output = StringIO()
        n = Align.write(alignments, output, "mauve")
        self.assertEqual(n, 6)
        output.seek(0)
        self.assertEqual(output.read(), data)


class TestSeparateFiles(unittest.TestCase):
    # Generate the output file separate.xmfa by running
    # progressiveMauve --solid-seeds equCab1.fa canFam2.fa mm9.fa --output=separate.xmfa

    sequences = {}
    for species in ("equCab1", "canFam2", "mm9"):
        filename = f"{species}.fa"
        path = os.path.join("Mauve", filename)
        record = SeqIO.read(path, "fasta")
        sequences[filename] = record.seq
        del filename
        del path
        del record

    def test_parse(self):
        path = os.path.join("Mauve", "separate.xmfa")
        saved_alignments = []
        with open(path) as stream:
            alignments = Align.parse(stream, "mauve")
            metadata = alignments.metadata
            self.assertEqual(len(metadata), 2)
            self.assertEqual(metadata["FormatVersion"], "Mauve1")
            self.assertEqual(metadata["BackboneFile"], "separate.xmfa.bbcols")
            identifiers = alignments.identifiers
            alignment = next(alignments)
            saved_alignments.append(alignment)
            self.assertEqual(len(alignment), 3)
            self.assertEqual(len(alignment.sequences), 3)
            self.assertEqual(alignment.sequences[0].id, "equCab1.fa")
            self.assertEqual(
                alignment.sequences[0].seq,
                Seq("GAAAAGGAAAGTACGGCCCGGCCACTCCGGGTGTGTGCTAGGAGGGCTTA"),
            )
            start = alignment.coordinates[0, 0]
            end = alignment.coordinates[0, -1]
            self.assertEqual(start, 50)
            self.assertEqual(end, 0)
            self.assertEqual(alignment.sequences[1].id, "canFam2.fa")
            self.assertEqual(
                alignment.sequences[1].seq,
                Seq("CAAGCCCTGCGCGCTCAGCCGGAGTGTCCCGGGCCCTGCTTTCCTTTTC"),
            )
            start = alignment.coordinates[1, 0]
            end = alignment.coordinates[1, -1]
            sequence = self.sequences[alignment.sequences[1].id]
            self.assertEqual(start, 0)
            self.assertEqual(end, 49)
            self.assertEqual(alignment.sequences[1].seq[start:end], sequence[start:end])
            self.assertEqual(alignment.sequences[2].id, "mm9.fa")

            sequence = alignment.sequences[2].seq
            start = len(sequence) - alignment.coordinates[2, 0]
            end = len(sequence) - alignment.coordinates[2, -1]
            self.assertEqual(start, 0)
            self.assertEqual(end, 19)
            sequence = self.sequences[alignment.sequences[2].id][start:end]
            self.assertEqual(alignment.sequences[2].seq[start:end], sequence)
            self.assertEqual(
                alignment[0], "TAAGCCCTCCTAGCACACACCCGGAGTGGCC-GGGCCGTAC-TTTCCTTTTC"
            )
            self.assertEqual(
                alignment[1], "CAAGCCCTGC--GCGCTCAGCCGGAGTGTCCCGGGCCCTGC-TTTCCTTTTC"
            )
            self.assertEqual(
                alignment[2], "---------------------------------GGATCTACTTTTCCTCTTC"
            )
            self.assertEqual(
                str(alignment),
                """\
equCab1.f        50 TAAGCCCTCCTAGCACACACCCGGAGTGGCC-GGGCCGTAC-TTTCCTTTTC  0
canFam2.f         0 CAAGCCCTGC--GCGCTCAGCCGGAGTGTCCCGGGCCCTGC-TTTCCTTTTC 49
mm9.fa           19 ---------------------------------GGATCTACTTTTCCTCTTC  0
""",
            )
            self.assertTrue(
                np.array_equal(
                    alignment.coordinates,
                    # fmt: off
                    np.array([[50, 40, 38, 19, 19, 18, 10, 10,  0],
                              [ 0, 10, 10, 29, 30, 31, 39, 39, 49],
                              [19, 19, 19, 19, 19, 19, 11, 10,  0]]),
                    # fmt: on
                )
            )
            self.assertEqual(
                alignment.format("mauve", metadata, identifiers),
                """\
> 1:1-50 - equCab1.fa
TAAGCCCTCCTAGCACACACCCGGAGTGGCC-GGGCCGTAC-TTTCCTTTTC
> 2:1-49 + canFam2.fa
CAAGCCCTGC--GCGCTCAGCCGGAGTGTCCCGGGCCCTGC-TTTCCTTTTC
> 3:1-19 - mm9.fa
---------------------------------GGATCTACTTTTCCTCTTC
=
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['T', 'A', 'A', 'G', 'C', 'C', 'C', 'T', 'C', 'C', 'T', 'A', 'G',
           'C', 'A', 'C', 'A', 'C', 'A', 'C', 'C', 'C', 'G', 'G', 'A', 'G',
           'T', 'G', 'G', 'C', 'C', '-', 'G', 'G', 'G', 'C', 'C', 'G', 'T',
           'A', 'C', '-', 'T', 'T', 'T', 'C', 'C', 'T', 'T', 'T', 'T', 'C'],
          ['C', 'A', 'A', 'G', 'C', 'C', 'C', 'T', 'G', 'C', '-', '-', 'G',
           'C', 'G', 'C', 'T', 'C', 'A', 'G', 'C', 'C', 'G', 'G', 'A', 'G',
           'T', 'G', 'T', 'C', 'C', 'C', 'G', 'G', 'G', 'C', 'C', 'C', 'T',
           'G', 'C', '-', 'T', 'T', 'T', 'C', 'C', 'T', 'T', 'T', 'T', 'C'],
          ['-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
           '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
           '-', '-', '-', '-', '-', '-', '-', 'G', 'G', 'A', 'T', 'C', 'T',
           'A', 'C', 'T', 'T', 'T', 'T', 'C', 'C', 'T', 'C', 'T', 'T', 'C']],
         dtype='U')
                    # fmt: on
                )
            )
            counts = alignment.counts()
            self.assertEqual(counts.left_insertions, 0)
            self.assertEqual(counts.left_deletions, 63)
            self.assertEqual(counts.right_insertions, 0)
            self.assertEqual(counts.right_deletions, 0)
            self.assertEqual(counts.internal_insertions, 3)
            self.assertEqual(counts.internal_deletions, 2)
            self.assertEqual(counts.left_gaps, 63)
            self.assertEqual(counts.right_gaps, 0)
            self.assertEqual(counts.internal_gaps, 5)
            self.assertEqual(counts.insertions, 3)
            self.assertEqual(counts.deletions, 65)
            self.assertEqual(counts.gaps, 68)
            self.assertEqual(counts.aligned, 84)
            self.assertEqual(counts.identities, 68)
            self.assertEqual(counts.mismatches, 16)
            alignment = next(alignments)
            saved_alignments.append(alignment)
            self.assertEqual(len(alignment), 1)
            self.assertEqual(len(alignment.sequences), 1)
            self.assertEqual(alignment.sequences[0].id, "mm9.fa")
            self.assertEqual(
                repr(alignment.sequences[0].seq),
                "Seq({19: 'CTGGCGTCCGGAGCTGGGACGT'}, length=41)",
            )
            sequence = self.sequences[alignment.sequences[0].id]
            start = alignment.coordinates[0, 0]
            end = alignment.coordinates[0, -1]
            self.assertEqual(alignment.sequences[0].seq[start:end], sequence[start:end])
            self.assertEqual(alignment[0], "CTGGCGTCCGGAGCTGGGACGT")
            self.assertTrue(np.array_equal(alignment.coordinates, np.array([[19, 41]])))
            self.assertEqual(
                str(alignment),
                """\
mm9.fa           19 CTGGCGTCCGGAGCTGGGACGT 41
""",
            )
            self.assertEqual(
                alignment.format("mauve", metadata, identifiers),
                """\
> 3:20-41 + mm9.fa
CTGGCGTCCGGAGCTGGGACGT
=
""",
            )
            self.assertTrue(
                np.array_equal(
                    np.array(alignment, "U"),
                    # fmt: off
np.array([['C', 'T', 'G', 'G', 'C', 'G', 'T', 'C', 'C', 'G', 'G',
           'A', 'G', 'C', 'T', 'G', 'G', 'G', 'A', 'C', 'G', 'T']], dtype='U')
                    # fmt: on
                )
            )
            counts = alignment.counts()
            self.assertEqual(counts.left_insertions, 0)
            self.assertEqual(counts.left_deletions, 0)
            self.assertEqual(counts.right_insertions, 0)
            self.assertEqual(counts.right_deletions, 0)
            self.assertEqual(counts.internal_insertions, 0)
            self.assertEqual(counts.internal_deletions, 0)
            self.assertEqual(counts.left_gaps, 0)
            self.assertEqual(counts.right_gaps, 0)
            self.assertEqual(counts.internal_gaps, 0)
            self.assertEqual(counts.insertions, 0)
            self.assertEqual(counts.deletions, 0)
            self.assertEqual(counts.gaps, 0)
            self.assertEqual(counts.aligned, 0)
            self.assertEqual(counts.identities, 0)
            self.assertEqual(counts.mismatches, 0)
            self.assertRaises(StopIteration, next, alignments)
        # As each nucleotide in each sequence is stored exactly once in an XMFA
        # file, we can reconstitute the full sequences:
        self.assertEqual(len(saved_alignments), 2)
        filenames = []
        for alignment in saved_alignments:
            for record in alignment.sequences:
                filename = record.id
                filenames.append(filename)
        filenames = set(filenames)
        n = len(filenames)
        self.assertEqual(n, 3)
        lengths = dict.fromkeys(filenames, 0)
        for alignment in saved_alignments:
            for record in alignment.sequences:
                filename = record.id
                length = len(record.seq)
                if length > lengths[filename]:
                    lengths[filename] = length
        self.assertEqual(lengths["equCab1.fa"], 50)
        self.assertEqual(lengths["canFam2.fa"], 49)
        self.assertEqual(lengths["mm9.fa"], 41)
        sequences = {}
        for filename, length in lengths.items():
            sequences[filename] = MutableSeq("N" * length)
        # Now fill up the sequences:
        for alignment in saved_alignments:
            for row, record in zip(alignment.coordinates, alignment.sequences):
                filename = record.id
                start = row[0]
                end = row[-1]
                if start > end:
                    start, end = end, start
                sequences[filename][start:end] = record.seq[start:end]
        # Confirm that the fully defined sequences agree with the Fasta file:
        for filename, sequence in sequences.items():
            sequences[filename] = Seq(sequence)
            self.assertEqual(sequences[filename], self.sequences[filename])
        # Make sure we can replace the partially defined sequences by these
        # fully defined sequences, and get the same alignment:
        alignment = saved_alignments[0]
        for record in alignment.sequences:
            filename = record.id
            record.seq = sequences[filename]
            self.assertEqual(
                alignment[0], "TAAGCCCTCCTAGCACACACCCGGAGTGGCC-GGGCCGTAC-TTTCCTTTTC"
            )
            self.assertEqual(
                alignment[1], "CAAGCCCTGC--GCGCTCAGCCGGAGTGTCCCGGGCCCTGC-TTTCCTTTTC"
            )
            self.assertEqual(
                alignment[2], "---------------------------------GGATCTACTTTTCCTCTTC"
            )
        alignment = saved_alignments[1]
        for record in alignment.sequences:
            filename = record.id
            record.seq = sequences[filename]
            self.assertEqual(alignment[0], "CTGGCGTCCGGAGCTGGGACGT")

    def test_write_read(self):
        path = os.path.join("Mauve", "separate.xmfa")
        with open(path) as stream:
            data = stream.read()

        stream = StringIO()
        stream.write(data)
        stream.seek(0)
        alignments = Align.parse(stream, "mauve")
        output = StringIO()
        n = Align.write(alignments, output, "mauve")
        self.assertEqual(n, 2)
        output.seek(0)
        self.assertEqual(output.read(), data)


class TestMauveBasic(unittest.TestCase):
    def test_empty(self):
        stream = StringIO()
        with self.assertRaisesRegex(ValueError, "Empty file."):
            Align.parse(stream, "mauve")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
