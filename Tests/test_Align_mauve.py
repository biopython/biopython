# Copyright 2006-2014 by Peter Cock.  All rights reserved.
# Revisions copyright 2011 Brandon Invergo. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.Align.mauve module."""
import os
import unittest

from io import StringIO

from Bio.Seq import Seq, MutableSeq
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
        path = os.path.join("Mauve", "combined.xmfa")
        saved_alignments = []
        with open(path) as stream:
            alignments = AlignmentIterator(stream)
            self.assertEqual(len(alignments.metadata), 2)
            self.assertEqual(alignments.metadata["FormatVersion"], "Mauve1")
            self.assertEqual(alignments.metadata["BackboneFile"], "combined.xmfa.bbcols")
            alignment = next(alignments)
            saved_alignments.append(alignment)
            self.assertEqual(len(alignment), 2)
            self.assertEqual(len(alignment.sequences), 2)
            self.assertEqual(alignment.sequences[0].id, "combined.fa:0")
            self.assertEqual(
                repr(alignment.sequences[0].seq),
                "Seq({1: 'AAAAGGAAAGTACGGCCCGGCCACTCCGGGTGTGTGCTAGGAGGGCTT'}, length=49)"
            )
            sequence = self.sequences[alignment.sequences[0].id]
            start = len(sequence) - alignment.coordinates[0, 0]
            end = len(sequence) - alignment.coordinates[0, -1]
            self.assertEqual(start, 1)
            self.assertEqual(end, 49)
            self.assertEqual(alignment.sequences[0].seq[start:end], sequence[start:end])
            self.assertEqual(alignment.sequences[1].id, "combined.fa:2")
            self.assertEqual(
                repr(alignment.sequences[1].seq),
                "Seq({1: 'AAGCCCTGCGCGCTCAGCCGGAGTGTCCCGGGCCCTGCTTTCCTTTT'}, length=48)"
            )
            start = alignment.coordinates[1, 0]
            end = alignment.coordinates[1, -1]
            self.assertEqual(start, 1)
            self.assertEqual(end, 48)
            sequence = self.sequences[alignment.sequences[1].id][start:end]
            self.assertEqual(alignment.sequences[1].seq[start:end], sequence)
            self.assertEqual(alignment[0], "AAGCCCTCCTAGCACACACCCGGAGTGG-CCGGGCCGTACTTTCCTTTT")
            self.assertEqual(alignment[1], "AAGCCCTGC--GCGCTCAGCCGGAGTGTCCCGGGCCCTGCTTTCCTTTT")
            self.assertTrue(
                numpy.array_equal(
                    alignment.coordinates,
                    numpy.array([[49, 40, 38, 21, 21, 1], [1, 10, 10, 27, 28, 48]]),
                )
            )
            alignment = next(alignments)
            saved_alignments.append(alignment)
            self.assertEqual(len(alignment), 1)
            self.assertEqual(len(alignment.sequences), 1)
            self.assertEqual(alignment.sequences[0].id, "combined.fa:0")
            self.assertEqual(alignment.sequences[0].seq, "G")
            sequence = self.sequences[alignment.sequences[0].id]
            start = alignment.coordinates[0, 0]
            end = alignment.coordinates[0, -1]
            self.assertEqual(alignment.sequences[0].seq[start:end], sequence[start:end])
            self.assertEqual(alignment[0], "G")
            self.assertTrue(
                numpy.array_equal(alignment.coordinates, numpy.array([[0, 1]]))
            )
            alignment = next(alignments)
            saved_alignments.append(alignment)
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
            self.assertTrue(
                numpy.array_equal(alignment.coordinates, numpy.array([[49, 50]]))
            )
            alignment = next(alignments)
            saved_alignments.append(alignment)
            self.assertEqual(len(alignment), 1)
            self.assertEqual(len(alignment.sequences), 1)
            self.assertEqual(alignment.sequences[0].id, "combined.fa:1")
            self.assertEqual(alignment.sequences[0].seq, "GAAGAGGAAAAGTAGATCCCTGGCGTCCGGAGCTGGGACGT")
            sequence = self.sequences[alignment.sequences[0].id]
            start = alignment.coordinates[0, 0]
            end = alignment.coordinates[0, -1]
            self.assertEqual(alignment.sequences[0].seq[start:end], sequence[start:end])
            self.assertEqual(alignment[0], "GAAGAGGAAAAGTAGATCCCTGGCGTCCGGAGCTGGGACGT")
            self.assertTrue(
                numpy.array_equal(alignment.coordinates, numpy.array([[0, 41]]))
            )
            alignment = next(alignments)
            saved_alignments.append(alignment)
            self.assertEqual(len(alignment), 1)
            self.assertEqual(len(alignment.sequences), 1)
            self.assertEqual(alignment.sequences[0].id, "combined.fa:2")
            self.assertEqual(alignment.sequences[0].seq, "C")
            sequence = self.sequences[alignment.sequences[0].id]
            start = alignment.coordinates[0, 0]
            end = alignment.coordinates[0, -1]
            self.assertEqual(alignment.sequences[0].seq[start:end], sequence[start:end])
            self.assertEqual(alignment[0], "C")
            self.assertTrue(
                numpy.array_equal(alignment.coordinates, numpy.array([[0, 1]]))
            )
            alignment = next(alignments)
            saved_alignments.append(alignment)
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
            self.assertEqual(alignment[0], "C")
            self.assertTrue(
                numpy.array_equal(alignment.coordinates, numpy.array([[48, 49]]))
            )
            self.assertRaises(StopIteration, next, alignments)
        # As each nucleotide in each sequence is stored exactly once in an XMFA
        # file, we can reconstitute the full sequences:
        self.assertEqual(len(saved_alignments), 6)
        maxindex = -1
        for alignment in saved_alignments:
            for record in alignment.sequences:
                filename, index = record.id.split(":")
                self.assertEqual(filename, "combined.fa")
                index = int(index)
                if index > maxindex:
                    maxindex = index
        n = maxindex + 1
        self.assertEqual(n, 3)
        lengths = [0] * n
        for alignment in saved_alignments:
            for record in alignment.sequences:
                filename, index = record.id.split(":")
                index = int(index)
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
                filename, index = record.id.split(":")
                index = int(index)
                start = row[0]
                end = row[-1]
                if start > end:
                    start, end = end, start
                sequences[index][start:end] = record.seq[start:end]
        # Confirm that the fully defined sequences agree with the Fasta file:
        filename = "combined.fa"
        for index, sequence in enumerate(sequences):
            sequences[index] = Seq(sequence)
            key = f"{filename}:{index}"
            self.assertEqual(sequences[index], self.sequences[key])
        # Make sure we can replace the partially defined sequences by these
        # fully defined sequences, and get the same alignment:
        alignment = saved_alignments[0]
        for record in alignment.sequences:
            filename, index = record.id.split(":")
            index = int(index)
            record.seq = sequences[index]
            self.assertEqual(alignment[0], "AAGCCCTCCTAGCACACACCCGGAGTGG-CCGGGCCGTACTTTCCTTTT")
            self.assertEqual(alignment[1], "AAGCCCTGC--GCGCTCAGCCGGAGTGTCCCGGGCCCTGCTTTCCTTTT")
        alignment = saved_alignments[1]
        for record in alignment.sequences:
            filename, index = record.id.split(":")
            index = int(index)
            record.seq = sequences[index]
            self.assertEqual(alignment[0], "G")
        alignment = saved_alignments[2]
        for record in alignment.sequences:
            filename, index = record.id.split(":")
            index = int(index)
            record.seq = sequences[index]
            self.assertEqual(alignment[0], "A")
        alignment = saved_alignments[3]
        for record in alignment.sequences:
            filename, index = record.id.split(":")
            index = int(index)
            record.seq = sequences[index]
            self.assertEqual(alignment[0], "GAAGAGGAAAAGTAGATCCCTGGCGTCCGGAGCTGGGACGT")
        alignment = saved_alignments[4]
        for record in alignment.sequences:
            filename, index = record.id.split(":")
            index = int(index)
            record.seq = sequences[index]
            self.assertEqual(alignment[0], "C")
        alignment = saved_alignments[5]
        for record in alignment.sequences:
            filename, index = record.id.split(":")
            index = int(index)
            record.seq = sequences[index]
            self.assertEqual(alignment[0], "C")

    def test_write_read(self):
        return
        with open(self.SIMPLE_XMFA) as stream:
            aln_list = list(MauveIterator(stream))

        stream = StringIO()
        MauveWriter(stream).write_file(aln_list)
        stream.seek(0)
        aln_list_out = list(MauveIterator(stream))

        for a1, a2 in zip(aln_list, aln_list_out):
            self.assertEqual(len(a1), len(a2))
            for r1, r2 in zip(a1, a2):
                self.assertEqual(r1.id, r2.id)
                self.assertEqual(r1.seq, r2.seq)


class TestSeparateFiles(unittest.TestCase):

    def setUp(self):
        self.sequences = {}
        for species in ("equCab1", "canFam2", "mm9"):
            filename = f"{species}.fa"
            path = os.path.join("Mauve", filename)
            record = SeqIO.read(path, "fasta")
            self.sequences[filename] = record.seq

    def test_parse(self):
        path = os.path.join("Mauve", "separate.xmfa")
        saved_alignments = []
        with open(path) as stream:
            alignments = AlignmentIterator(stream)
            self.assertEqual(len(alignments.metadata), 2)
            self.assertEqual(alignments.metadata["FormatVersion"], "Mauve1")
            self.assertEqual(alignments.metadata["BackboneFile"], "separate.xmfa.bbcols")
            alignment = next(alignments)
            saved_alignments.append(alignment)
            self.assertEqual(len(alignment), 2)
            self.assertEqual(len(alignment.sequences), 2)
            self.assertEqual(alignment.sequences[0].id, "canFam2.fa")
            self.assertEqual(
                repr(alignment.sequences[0].seq),
                "Seq({25: 'GTCCCGGGCCCTGCTTTCCTTTTC'}, length=49)"
            )
            start = alignment.coordinates[0, 0]
            end = alignment.coordinates[0, -1]
            sequence = self.sequences[alignment.sequences[0].id]
            self.assertEqual(start, 25)
            self.assertEqual(end, 49)
            self.assertEqual(alignment.sequences[0].seq[start:end], sequence[start:end])
            self.assertEqual(alignment.sequences[1].id, "mm9.fa")

            sequence = alignment.sequences[1].seq
            start = len(sequence) - alignment.coordinates[1, 0]
            end = len(sequence) - alignment.coordinates[1, -1]
            self.assertEqual(start, 0)
            self.assertEqual(end, 24)
            sequence = self.sequences[alignment.sequences[1].id][start:end]
            self.assertEqual(alignment.sequences[1].seq[start:end], sequence)
            self.assertEqual(alignment[0], "GTCCCGGGCCCTGCTTTCCTTTTC")
            self.assertEqual(alignment[1], "GCCAGGGATCTACTTTTCCTCTTC")
            self.assertTrue(
                numpy.array_equal(
                    alignment.coordinates, numpy.array([[25, 49], [24, 0]]),
                )
            )
            alignment = next(alignments)
            saved_alignments.append(alignment)
            self.assertEqual(len(alignment), 1)
            self.assertEqual(len(alignment.sequences), 1)
            self.assertEqual(alignment.sequences[0].id, "equCab1.fa")
            self.assertEqual(alignment.sequences[0].seq, "GAAAAGGAAAGTACGGCCCGGCCACTCCGGGTGTGTGCTAGGAGGGCTTA")
            sequence = self.sequences[alignment.sequences[0].id]
            start = alignment.coordinates[0, 0]
            end = alignment.coordinates[0, -1]
            self.assertEqual(alignment.sequences[0].seq[start:end], sequence[start:end])
            self.assertEqual(alignment[0], "GAAAAGGAAAGTACGGCCCGGCCACTCCGGGTGTGTGCTAGGAGGGCTTA")
            self.assertTrue(
                numpy.array_equal(alignment.coordinates, numpy.array([[0, 50]]))
            )
            alignment = next(alignments)
            saved_alignments.append(alignment)
            self.assertEqual(len(alignment), 1)
            self.assertEqual(len(alignment.sequences), 1)
            self.assertEqual(alignment.sequences[0].id, "canFam2.fa")
            self.assertEqual(alignment.sequences[0].seq, "CAAGCCCTGCGCGCTCAGCCGGAGT")
            sequence = self.sequences[alignment.sequences[0].id]
            start = alignment.coordinates[0, 0]
            end = alignment.coordinates[0, -1]
            self.assertEqual(alignment.sequences[0].seq[start:end], sequence[start:end])
            self.assertEqual(alignment[0], "CAAGCCCTGCGCGCTCAGCCGGAGT")
            self.assertTrue(
                numpy.array_equal(alignment.coordinates, numpy.array([[0, 25]]))
            )
            alignment = next(alignments)
            saved_alignments.append(alignment)
            self.assertEqual(len(alignment), 1)
            self.assertEqual(len(alignment.sequences), 1)
            self.assertEqual(alignment.sequences[0].id, "mm9.fa")
            sequence = self.sequences[alignment.sequences[0].id]
            start = alignment.coordinates[0, 0]
            end = alignment.coordinates[0, -1]
            self.assertEqual(start, 24)
            self.assertEqual(end, 41)
            self.assertEqual(alignment.sequences[0].seq[start:end], "GTCCGGAGCTGGGACGT")
            sequence = self.sequences[alignment.sequences[0].id]
            start = alignment.coordinates[0, 0]
            end = alignment.coordinates[0, -1]
            self.assertEqual(alignment.sequences[0].seq[start:end], sequence[start:end])
            self.assertEqual(alignment[0], "GTCCGGAGCTGGGACGT")
            self.assertTrue(
                numpy.array_equal(alignment.coordinates, numpy.array([[24, 41]]))
            )
            self.assertRaises(StopIteration, next, alignments)
        # As each nucleotide in each sequence is stored exactly once in an XMFA
        # file, we can reconstitute the full sequences:
        self.assertEqual(len(saved_alignments), 4)
        filenames = []
        for alignment in saved_alignments:
            for record in alignment.sequences:
                filename = record.id
                filenames.append(filename)
        filenames = set(filenames)
        n = len(filenames)
        self.assertEqual(n, 3)
        lengths = {filename: 0 for filename in filenames} 
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
            self.assertEqual(alignment[0], "GTCCCGGGCCCTGCTTTCCTTTTC")
            self.assertEqual(alignment[1], "GCCAGGGATCTACTTTTCCTCTTC")
        alignment = saved_alignments[1]
        for record in alignment.sequences:
            filename = record.id
            record.seq = sequences[filename]
            self.assertEqual(alignment[0], "GAAAAGGAAAGTACGGCCCGGCCACTCCGGGTGTGTGCTAGGAGGGCTTA")
        alignment = saved_alignments[2]
        for record in alignment.sequences:
            filename = record.id
            record.seq = sequences[filename]
            self.assertEqual(alignment[0], "CAAGCCCTGCGCGCTCAGCCGGAGT")
        alignment = saved_alignments[3]
        for record in alignment.sequences:
            filename = record.id
            record.seq = sequences[filename]
            self.assertEqual(alignment[0], "GTCCGGAGCTGGGACGT")

    def test_write_read(self):
        return
        with open(self.SIMPLE_XMFA) as stream:
            aln_list = list(MauveIterator(stream))

        stream = StringIO()
        MauveWriter(stream).write_file(aln_list)
        stream.seek(0)
        aln_list_out = list(MauveIterator(stream))

        for a1, a2 in zip(aln_list, aln_list_out):
            self.assertEqual(len(a1), len(a2))
            for r1, r2 in zip(a1, a2):
                self.assertEqual(r1.id, r2.id)
                self.assertEqual(r1.seq, r2.seq)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
