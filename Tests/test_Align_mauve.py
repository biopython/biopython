# Copyright 2006-2014 by Peter Cock.  All rights reserved.
# Revisions copyright 2011 Brandon Invergo. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.Align.mauve module."""
import os
import unittest

from io import StringIO

from Bio import SeqIO
from Bio.Align.mauve import AlignmentIterator
# from Bio.AlignIO.MauveIO import MauveWriter


try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align.mauve."
    ) from None


class TestCombinedFile(unittest.TestCase):

    def setUp(self):
        filename = "combined.fa"
        path = os.path.join("Mauve", filename)
        records = SeqIO.parse(path, "fasta")
        self.sequences = {f"{filename}:{index}": record.seq for index, record in enumerate(records)}

    def test_parse(self):
        ids = []
        path = os.path.join("Mauve", "combined.xmfa")
        with open(path) as handle:
            alignments = AlignmentIterator(handle)
            self.assertEqual(len(alignments.metadata), 2)
            self.assertEqual(alignments.metadata["FormatVersion"], "Mauve1")
            self.assertEqual(alignments.metadata["BackboneFile"], "combined.xmfa.bbcols")
            alignment = next(alignments)
            self.assertEqual(len(alignment), 2)
            self.assertEqual(len(alignment.sequences), 2)
            self.assertEqual(alignment.sequences[0].id, "combined.fa:0")
            self.assertEqual(
                repr(alignment.sequences[0].seq),
                "Seq({1: 'AAAAGGAAAGTACGGCCCGGCCACTCCGGGTGTGTGCTAGGAGGGCTT'}, length=49)"
            )
            sequence = self.sequences[alignment.sequences[0].id]
            start = alignment.coordinates[0, 0]
            end = alignment.coordinates[0, -1]
            self.assertEqual(alignment.sequences[0].seq[start:end], sequence[start:end])
            self.assertEqual(alignment.sequences[1].id, "combined.fa:2")
            self.assertEqual(
                repr(alignment.sequences[1].seq),
                "Seq({1: 'AAGCCCTGCGCGCTCAGCCGGAGTGTCCCGGGCCCTGCTTTCCTTTT'}, length=48)"
            )
            start = alignment.coordinates[1, 0]
            end = alignment.coordinates[1, -1]
            if start < end:
                sequence = self.sequences[alignment.sequences[1].id][start:end]
                print(repr(alignment.sequences[1].seq[start:end]), start, end)
                self.assertEqual(alignment.sequences[1].seq[start:end], sequence)
            else:
                sequence = self.sequences[alignment.sequences[1].id][end:start]
                self.assertEqual(alignment.sequences[1].seq[end:start], sequence)
            print(alignment.coordinates)
            self.assertEqual(alignment[0], "AAGCCCTCCTAGCACACACCCGGAGTGG-CCGGGCCGTACTTTCCTTTT")
            self.assertEqual(alignment[1], "AAGCCCTGC--GCGCTCAGCCGGAGTGTCCCGGGCCCTGCTTTCCTTTT")
            self.assertTrue(
                numpy.array_equal(
                    alignment.coordinates,
                    numpy.array([[0, 15, 24, 74, 84, 131], [0, 15, 15, 65, 65, 112]]),
                )
            )
            alignment = next(alignments)
            self.assertEqual(len(alignment), 1)
            self.assertEqual(len(alignment.sequences), 1)
            self.assertEqual(alignment.sequences[0].id, "combined.fa:0")
            self.assertEqual(alignment.sequences[0].seq, "G")
            sequence = self.sequences[alignment.sequences[0].id]
            start = alignment.coordinates[0, 0]
            end = alignment.coordinates[0, -1]
            self.assertEqual(alignment.sequences[0].seq[start:end], sequence[start:end])
            self.assertEqual(alignment[0], "G")
            print(alignment.coordinates)
            self.assertTrue(
                numpy.array_equal(
                    alignment.coordinates,
                    numpy.array([[0, 15, 24, 74, 84, 131], [0, 15, 15, 65, 65, 112]]),
                )
            )
            alignment = next(alignments)
            self.assertEqual(len(alignment), 1)
            self.assertEqual(len(alignment.sequences), 1)
            self.assertEqual(alignment.sequences[0].id, "combined.fa:0")
            self.assertEqual(
                repr(alignment.sequences[0].seq),
                "Seq({49: 'A'}, length=50)"
            )
            sequence = self.sequences[alignment.sequences[0].id]
            start = alignment.coordinates[0, 0]
            end = alignment.coordinates[0, -1]
            self.assertEqual(alignment.sequences[0].seq[start:end], sequence[start:end])
            self.assertEqual(alignment[0], "A")
            print(alignment.coordinates)
            self.assertTrue(
                numpy.array_equal(
                    alignment.coordinates,
                    numpy.array([[0, 15, 24, 74, 84, 131], [0, 15, 15, 65, 65, 112]]),
                )
            )
            alignment = next(alignments)
            self.assertEqual(len(alignment), 1)
            self.assertEqual(len(alignment.sequences), 1)
            self.assertEqual(alignment.sequences[0].id, "combined.fa:1")
            self.assertEqual(alignment.sequences[0].seq, "GAAGAGGAAAAGTAGATCCCTGGCGTCCGGAGCTGGGACGT")
            sequence = self.sequences[alignment.sequences[0].id]
            start = alignment.coordinates[0, 0]
            end = alignment.coordinates[0, -1]
            self.assertEqual(alignment.sequences[0].seq[start:end], sequence[start:end])
            self.assertEqual(alignment[0], "C")
            print(alignment.coordinates)
            self.assertTrue(
                numpy.array_equal(
                    alignment.coordinates,
                    numpy.array([[0, 15, 24, 74, 84, 131], [0, 15, 15, 65, 65, 112]]),
                )
            )
            alignment = next(alignments)
            self.assertEqual(len(alignment), 1)
            self.assertEqual(len(alignment.sequences), 1)
            self.assertEqual(alignment.sequences[0].id, "combined.fa:2")
            self.assertEqual(alignment.sequences[0].seq, "C")
            sequence = self.sequences[alignment.sequences[0].id]
            start = alignment.coordinates[0, 0]
            end = alignment.coordinates[0, -1]
            self.assertEqual(alignment.sequences[0].seq[start:end], sequence[start:end])
            self.assertEqual(alignment[0], "C")
            print(alignment.coordinates)
            self.assertTrue(
                numpy.array_equal(
                    alignment.coordinates,
                    numpy.array([[0, 15, 24, 74, 84, 131], [0, 15, 15, 65, 65, 112]]),
                )
            )
            alignment = next(alignments)
            self.assertEqual(len(alignment), 1)
            self.assertEqual(len(alignment.sequences), 1)
            self.assertEqual(alignment.sequences[0].id, "combined.fa:2")
            self.assertEqual(
                repr(alignment.sequences[0].seq),
                "Seq({48: 'C'}, length=49)"
            )
            sequence = self.sequences[alignment.sequences[0].id]
            start = alignment.coordinates[0, 0]
            end = alignment.coordinates[0, -1]
            self.assertEqual(alignment.sequences[0].seq[start:end], sequence[start:end])
            self.assertEqual(alignment[0], "AAGCCCTCCTAGCACACACCCGGAGTGG-CCGGGCCGTACTTTCCTTTT")
            self.assertEqual(alignment[1], "AAGCCCTGC--GCGCTCAGCCGGAGTGTCCCGGGCCCTGCTTTCCTTTT")
            print(alignment.coordinates)
            self.assertTrue(
                numpy.array_equal(
                    alignment.coordinates,
                    numpy.array([[0, 15, 24, 74, 84, 131], [0, 15, 15, 65, 65, 112]]),
                )
            )
            self.assertRaises(StopIteration, next, alignments)


    def test_write_read(self):
        return
        with open(self.SIMPLE_XMFA) as handle:
            aln_list = list(MauveIterator(handle))

        handle = StringIO()
        MauveWriter(handle).write_file(aln_list)
        handle.seek(0)
        aln_list_out = list(MauveIterator(handle))

        for a1, a2 in zip(aln_list, aln_list_out):
            self.assertEqual(len(a1), len(a2))
            for r1, r2 in zip(a1, a2):
                self.assertEqual(r1.id, r2.id)
                self.assertEqual(r1.seq, r2.seq)


class TestSeparateFiles(unittest.TestCase):

    def test_parse(self):
        ids = []
        path = os.path.join("Mauve", "separate.xmfa")
        with open(path) as handle:
            alignments = AlignmentIterator(handle)
            self.assertEqual(len(alignments.metadata), 2)
            self.assertEqual(alignments.metadata["FormatVersion"], "Mauve1")
            self.assertEqual(alignments.metadata["BackboneFile"], "separate.xmfa.bbcols")
            for alignment in alignments:
                for sequence in alignment.sequences:
                    ids.append(sequence.id)

        self.assertEqual(
            ids, ["canFam2.fa", "mm9.fa", "equCab1.fa", "canFam2.fa", "mm9.fa"],
        )
        return

        expected = expected.replace(" ", "").replace("\n", "")
        self.assertEqual(
            repr(sequence),
            "Seq({11410: '%s'})" % expected,
        )

    def test_sequence_positions(self):
        seqs = {}
        for species in ("canFam2", "equCab1", "mm9"):
            filename = "%s.fa" % species
            path = os.path.join("Mauve", filename)
            record = SeqIO.read(path, "fasta")
            seqs[filename] = record.seq

        path = "Mauve/separate.xmfa"
        with open(path) as stream:
            aln_list = list(AlignmentIterator(stream))

        for aln in aln_list:
            for index, record in enumerate(aln.sequences):
                filename = record.id
                fasta_seq = seqs[filename]
                start = aln.coordinates[index, 0]
                end = aln.coordinates[index, -1]
                self.assertEqual(record.seq[start:end], fasta_seq[start:end])

    def test_write_read(self):
        return
        with open(self.SIMPLE_XMFA) as handle:
            aln_list = list(MauveIterator(handle))

        handle = StringIO()
        MauveWriter(handle).write_file(aln_list)
        handle.seek(0)
        aln_list_out = list(MauveIterator(handle))

        for a1, a2 in zip(aln_list, aln_list_out):
            self.assertEqual(len(a1), len(a2))
            for r1, r2 in zip(a1, a2):
                self.assertEqual(r1.id, r2.id)
                self.assertEqual(r1.seq, r2.seq)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
