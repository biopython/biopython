# Copyright 2022 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Align.psl module."""
import unittest


from Bio.Align import Alignment, psl
from Bio.Seq import Seq
from Bio import SeqIO


try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align.psl."
    ) from None


class TestAlign_reading(unittest.TestCase):
    def setUp(self):
        records = SeqIO.parse("dna.fa", "fasta")
        self.dna = {record.id: record.seq for record in records}
        records = SeqIO.parse("rna.fa", "fasta")
        self.rna = {record.id: record.seq for record in records}

    def test_reading_dna_rna(self):
        """Test parsing dna_rna.psl."""
        path = "dna_rna.psl"
        alignments = psl.AlignmentIterator(path)
        self.assertEqual(alignments.metadata["version"], "3")
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 207)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 3707))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chrY:22316600-22320400")
        self.assertEqual(alignment.query.id, "NR_104151.1")
        self.assertEqual(len(alignment.target.seq), 3800)
        self.assertEqual(len(alignment.query.seq), 207)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[  77,  147, 2182, 2297, 3762, 3784],
                             [   0,   70,   70,  185,  185,  207]]),
            )
        )
        alignment.target.seq = self.dna[alignment.target.id]
        alignment.query.seq = self.rna[alignment.query.id]
        self.assertTrue(numpy.array_equal(
            alignment.substitutions,
            numpy.array([[64.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0., 44.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0., 52.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  0., 47.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                        ])
                )
        )
        self.assertEqual(alignment.substitutions.alphabet, "ACGTacgt")
        matches = sum(alignment.substitutions[c,c] for c in alignment.substitutions.alphabet)
        repMatches = sum(alignment.substitutions[c,c.swapcase()] for c in alignment.substitutions.alphabet)
        self.assertEqual(matches, alignment.matches)
        self.assertEqual(repMatches, alignment.repMatches)
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 175)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 6)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 1711))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr3:42530800-42532700")
        self.assertEqual(alignment.query.id, "NR_046654.1")
        self.assertEqual(len(alignment.target.seq), 1900)
        self.assertEqual(len(alignment.query.seq), 181)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[  95,  158, 1220, 1295, 1763, 1806],
                             [ 181,  118,  118,   43,   43,    0]]),
            )
        )
        alignment.target.seq = self.dna[alignment.target.id]
        alignment.query.seq = self.rna[alignment.query.id]
        self.assertTrue(numpy.array_equal(
            alignment.substitutions,
            numpy.array([[36.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0., 40.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0., 57.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  0., 42.,  0.,  0.,  0.,  0.],
                         [ 2.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  1.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  3.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                        ])
            )
        )
        self.assertEqual(alignment.substitutions.alphabet, "ACGTacgt")
        matches = sum(alignment.substitutions[c,c] for c in alignment.substitutions.alphabet)
        repMatches = sum(alignment.substitutions[c,c.swapcase()] for c in alignment.substitutions.alphabet)
        self.assertEqual(matches, alignment.matches)
        self.assertEqual(repMatches, alignment.repMatches)
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 207)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 3707))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chrY:22316600-22320400")
        self.assertEqual(alignment.query.id, "NR_104151.1_extended")
        self.assertEqual(len(alignment.target.seq), 3800)
        self.assertEqual(len(alignment.query.seq), 215)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[  77,  147, 2182, 2297, 3762, 3784],
                             [   3,   73,   73,  188,  188,  210]]),
            )
        )
        alignment.target.seq = self.dna[alignment.target.id]
        alignment.query.seq = self.rna[alignment.query.id]
        self.assertTrue(numpy.array_equal(
            alignment.substitutions,
            numpy.array([[64.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0., 44.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0., 52.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  0., 47.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                        ]),
            )
        )
        self.assertEqual(alignment.substitutions.alphabet, "ACGTXacgt")
        matches = sum(alignment.substitutions[c,c] for c in alignment.substitutions.alphabet)
        repMatches = sum(alignment.substitutions[c,c.swapcase()] for c in alignment.substitutions.alphabet if c != "X")
        self.assertEqual(matches, alignment.matches)
        self.assertEqual(repMatches, alignment.repMatches)
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 175)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 6)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 1711))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr3:42530800-42532700")
        self.assertEqual(alignment.query.id, "NR_046654.1_extended")
        self.assertEqual(len(alignment.target.seq), 1900)
        self.assertEqual(len(alignment.query.seq), 194)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[  95,  158, 1220, 1295, 1763, 1806],
                             [ 184,  121,  121,   46,   46,    3]]),
            )
        )
        alignment.target.seq = self.dna[alignment.target.id]
        alignment.query.seq = self.rna[alignment.query.id]
        self.assertTrue(numpy.array_equal(
            alignment.substitutions,
            numpy.array([[36.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0., 40.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0., 57.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  0., 42.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 2.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  3.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                        ]),
                )
        )
        self.assertEqual(alignment.substitutions.alphabet, "ACGTXacgt")
        matches = sum(alignment.substitutions[c,c] for c in alignment.substitutions.alphabet)
        repMatches = sum(alignment.substitutions[c,c.swapcase()] for c in alignment.substitutions.alphabet if c != "X")
        self.assertEqual(matches, alignment.matches)
        self.assertEqual(repMatches, alignment.repMatches)
        self.assertRaises(StopIteration, next, alignments)

    def test_reading_psl_34_001(self):
        """Test parsing psl_34_001.psl."""
        path = "psl_34_001.psl"
        alignments = psl.AlignmentIterator(path)
        self.assertEqual(alignments.metadata["version"], "3")
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 16)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 16))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 191154276)
        self.assertEqual(len(alignment.query.seq), 33)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[61646095, 61646111],
                             [      11,       27]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 33)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 33))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 33)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[10271783, 10271816],
                             [       0,       33]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 17)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 17))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 33)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[53575980, 53575997],
                             [      25,        8]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 38)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr9")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 141213431)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[85737865, 85737906],
                             [       9,       50]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 41)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr8")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 146364022)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[95160479, 95160520],
                             [       8,       49]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 33)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 51304566)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[42144400, 42144436],
                             [      11,       47]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 43)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 48))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[183925984, 183925990, 183925990, 183926028],
                             [        1,         7,        11,        49]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 34)
        self.assertEqual(alignment.misMatches, 2)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 170))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[35483340, 35483365, 35483499, 35483510],
                             [      10,       35,       35,       46]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 39)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 78077248)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[23891310, 23891349],
                             [      10,       49]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 27)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 28))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 78077248)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[43252217, 43252245],
                             [      21,       49]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 44)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 54))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 115169878)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[52759147, 52759154, 52759160, 52759160, 52759198],
                             [       1,        8,        8,       11,       49]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 50)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 50))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[1207056, 1207106],
                             [      0,      50]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 31)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[61700837, 61700871],
                             [       1,       35]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 28)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 44))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 191154276)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[37558157, 37558167, 37558173, 37558173, 37558191],
                             [      49,       39,       39,       29,       11]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 35)
        self.assertEqual(alignment.misMatches, 2)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 37))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 51304566)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[48997405, 48997442],
                             [      49,       12]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 35)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[120641740, 120641776],
                             [       49,        13]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 39)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[54017130, 54017169],
                             [      49,       10]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 36)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[553742, 553781],
                             [    49,     10]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 33)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 135534747)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[99388555, 99388591],
                             [      49,       13]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 24)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 25))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 135534747)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[112178171, 112178196],
                             [       35,        10]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 35)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[39368490, 39368526],
                             [      49,       13]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 33)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[220325687, 220325721],
                             [       47,        13]]),
            )
        )
        self.assertRaises(StopIteration, next, alignments)

    def test_reading_psl_34_002(self):
        """Test parsing psl_34_002.psl."""
        path = "psl_34_002.psl"
        alignments = psl.AlignmentIterator(path)
        self.assertEqual(alignments.metadata["version"], "3")
        self.assertRaises(StopIteration, next, alignments)

    def test_reading_psl_34_003(self):
        """Test parsing psl_34_003.psl."""
        path = "psl_34_003.psl"
        alignments = psl.AlignmentIterator(path)
        self.assertEqual(alignments.metadata["version"], "3")
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 16)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 16))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 191154276)
        self.assertEqual(len(alignment.query.seq), 33)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[61646095, 61646111],
                             [      11,       27]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 33)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 33))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 33)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[10271783, 10271816],
                             [       0,       33]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 17)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 17))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 33)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[53575980, 53575997],
                             [      25,        8]]),
            )
        )
        self.assertRaises(StopIteration, next, alignments)

    def test_reading_psl_34_004(self):
        """Test parsing psl_34_004.psl."""
        path = "psl_34_004.psl"
        alignments = psl.AlignmentIterator(path)
        self.assertEqual(alignments.metadata["version"], "3")
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 38)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr9")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 141213431)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[85737865, 85737906],
                             [       9,       50]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 41)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr8")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 146364022)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[95160479, 95160520],
                             [       8,       49]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 33)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 51304566)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[42144400, 42144436],
                             [      11,       47]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 43)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 48))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[183925984, 183925990, 183925990, 183926028],
                             [        1,         7,        11,        49]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 34)
        self.assertEqual(alignment.misMatches, 2)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 170))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[35483340, 35483365, 35483499, 35483510],
                             [      10,       35,       35,       46]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 39)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 78077248)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[23891310, 23891349],
                             [      10,       49]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 27)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 28))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 78077248)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[43252217, 43252245],
                             [      21,       49]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 44)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 54))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 115169878)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[52759147, 52759154, 52759160, 52759160, 52759198],
                             [       1,        8,        8,       11,       49]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 50)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 50))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[1207056, 1207106],
                             [      0,      50]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 31)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[61700837, 61700871],
                             [       1,       35]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 28)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 44))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 191154276)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[37558157, 37558167, 37558173, 37558173, 37558191],
                             [      49,       39,       39,       29,       11]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 35)
        self.assertEqual(alignment.misMatches, 2)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 37))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 51304566)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[48997405, 48997442],
                             [      49,       12]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 35)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[120641740, 120641776],
                             [       49,        13]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 39)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[54017130, 54017169],
                             [      49,       10]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 36)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[553742, 553781],
                             [    49,     10]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 33)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 135534747)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[99388555, 99388591],
                             [      49,       13]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 24)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 25))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 135534747)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[112178171, 112178196],
                             [       35,        10]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 35)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[39368490, 39368526],
                             [      49,       13]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 33)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[220325687, 220325721],
                             [       47,        13]]),
            )
        )
        self.assertRaises(StopIteration, next, alignments)

    def test_reading_psl_34_005(self):
        """Test parsing psl_34_005.psl."""
        path = "psl_34_005.psl"
        alignments = psl.AlignmentIterator(path)
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 16)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 16))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 191154276)
        self.assertEqual(len(alignment.query.seq), 33)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[61646095, 61646111],
                             [      11,       27]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 33)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 33))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 33)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[10271783, 10271816],
                             [       0,       33]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 17)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 17))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 33)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[53575980, 53575997],
                             [      25,        8]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 38)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr9")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 141213431)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[85737865, 85737906],
                             [       9,       50]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 41)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr8")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 146364022)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[95160479, 95160520],
                             [       8,       49]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 33)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 51304566)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[42144400, 42144436],
                             [      11,       47]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 43)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 48))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[183925984, 183925990, 183925990, 183926028],
                             [        1,         7,        11,        49]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 34)
        self.assertEqual(alignment.misMatches, 2)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 170))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[35483340, 35483365, 35483499, 35483510],
                             [      10,       35,       35,       46]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 39)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 78077248)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[23891310, 23891349],
                             [      10,       49]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 27)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 28))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 78077248)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[43252217, 43252245],
                             [      21,       49]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 44)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 54))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 115169878)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[52759147, 52759154, 52759160, 52759160, 52759198],
                             [       1,        8,        8,       11,       49]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 50)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 50))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[1207056, 1207106],
                             [      0,      50]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 31)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[61700837, 61700871],
                             [       1,       35]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 28)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 44))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 191154276)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[37558157, 37558167, 37558173, 37558173, 37558191],
                             [      49,       39,       39,       29,       11]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 35)
        self.assertEqual(alignment.misMatches, 2)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 37))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 51304566)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[48997405, 48997442],
                             [      49,       12]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 35)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[120641740, 120641776],
                             [       49,        13]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 39)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[54017130, 54017169],
                             [      49,       10]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 36)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[553742, 553781],
                             [    49,     10]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 33)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 135534747)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[99388555, 99388591],
                             [      49,       13]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 24)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 25))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 135534747)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[112178171, 112178196],
                             [       35,        10]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 35)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[39368490, 39368526],
                             [      49,       13]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 33)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[220325687, 220325721],
                             [       47,        13]]),
            )
        )
        self.assertRaises(StopIteration, next, alignments)

    def test_reading_psl_35_001(self):
        """Test parsing psl_35_001.psl."""
        path = "psl_35_001.psl"
        alignments = psl.AlignmentIterator(path)
        self.assertEqual(alignments.metadata["version"], "3")
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 52)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertEqual(len(alignment.target.seq), 114364328)
        self.assertEqual(len(alignment.query.seq), 230)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[75566694, 75566850],
                             [      61,      113]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 44)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertEqual(len(alignment.target.seq), 114364328)
        self.assertEqual(len(alignment.query.seq), 230)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[75560749, 75560881],
                             [      17,       61]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 44)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertEqual(len(alignment.target.seq), 114364328)
        self.assertEqual(len(alignment.query.seq), 230)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[75549820, 75549865, 75567225, 75567225, 75567312],
                             [       0,       15,       15,      113,      142]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 47)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertEqual(len(alignment.target.seq), 114364328)
        self.assertEqual(len(alignment.query.seq), 230)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[75604767, 75604827, 75605728, 75605809],
                             [     183,      203,      203,      230]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 25)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertEqual(len(alignment.target.seq), 114364328)
        self.assertEqual(len(alignment.query.seq), 230)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[75594914, 75594989],
                             [     158,      183]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 16)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertEqual(len(alignment.target.seq), 114364328)
        self.assertEqual(len(alignment.query.seq), 230)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[75569459, 75569507],
                             [     142,      158]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 26)
        self.assertEqual(alignment.misMatches, 8)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertEqual(len(alignment.target.seq), 190214555)
        self.assertEqual(len(alignment.query.seq), 230)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[41260685, 41260787],
                             [      76,      110]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 37)
        self.assertEqual(alignment.misMatches, 26)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertEqual(len(alignment.target.seq), 190214555)
        self.assertEqual(len(alignment.query.seq), 230)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[41257605, 41257731, 41263227, 41263227, 41263290],
                             [      17,       59,       59,      162,      183]]),
            )
        )
        self.assertRaises(StopIteration, next, alignments)

    def test_reading_psl_35_002(self):
        """Test parsing psl_35_002.psl."""
        path = "psl_35_002.psl"
        alignments = psl.AlignmentIterator(path)
        self.assertEqual(alignments.metadata["version"], "3")
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 210)
        self.assertEqual(alignment.misMatches, 3)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "KI537979")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertEqual(len(alignment.target.seq), 14052872)
        self.assertEqual(len(alignment.query.seq), 230)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[9712654, 9712786, 9715941, 9716097, 9716445, 9716532, 9718374,
                              9718422, 9739264, 9739339, 9743706, 9743766, 9744511, 9744592],
                             [     17,      61,      61,     113,     113,     142,     142,
                                  158,     158,     183,     183,     203,     203,     230]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 207)
        self.assertEqual(alignment.misMatches, 22)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "KI538594")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertEqual(len(alignment.target.seq), 7819582)
        self.assertEqual(len(alignment.query.seq), 230)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[2103463, 2103523, 2103522, 2103522, 2104149],
                             [      0,      20,      20,      21,     230]]),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 204)
        self.assertEqual(alignment.misMatches, 6)
        self.assertEqual(alignment.repMatches, 0)
        self.assertEqual(alignment.nCount, 0)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "KI537194")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertEqual(len(alignment.target.seq), 37111980)
        self.assertEqual(len(alignment.query.seq), 230)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[20873021, 20872472, 20872471, 20872471, 20872390],
                             [       0,      183,      183,      203,      230]]),
            )
        )
        self.assertRaises(StopIteration, next, alignments)
        # These alignments were generated using
        # blat -t=dnax -q=prot balAcu1.2bit CAG33136.1.fasta psl_35_002.psl
        # As a bonus, for the last alignment, let's extract the translated DNA
        # to protein alignment.
        length = len(alignment.sequences[0].seq)
        self.assertEqual(length, 37111980)
        # To avoid having to include a 37Mb Fasta file in Biopython, we used
        # twoBitToFa -noMask balAcu1.2bit:KI537194:20872390-20873021 KI537194_part.fa
        # to create a Fasta file containing only the relevant DNA sequence.
        filename = "KI537194_part.fa"
        dna = SeqIO.read(filename, "fasta")
        name, start_end = dna.id.split(":")
        self.assertEqual(alignment.sequences[0].id, name)
        start, end = start_end.split("-")
        start = int(start)
        end = int(end)
        self.assertEqual(start, 20872390)
        self.assertEqual(end, 20873021)
        sequence = str(dna.seq)
        dna = Seq({start: sequence}, length=length)
        # Now we have a partially defined sequence with the appropriate length,
        # and with the sequence contents starting at the correct position in the
        # genome.
        alignment.sequences[0] = dna
        # Load the protein sequence:
        filename = "CAG33136.1.fasta"
        protein = SeqIO.read(filename, "fasta")
        self.assertEqual(alignment.sequences[1].id, protein.id)
        alignment.sequences[1] = protein.seq
        # The alignment is on the reverse strand of the DNA sequence:
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        # so we take the reverse complement:
        alignment.coordinates[0, :] = len(dna) - alignment.coordinates[0, :]
        alignment.sequences[0] = dna.reverse_complement()
        # Now the alignment is on the forward strand:
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        # The protein alignment was already in the forward orientation:
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        # Now extract the aligned sequences:
        aligned_dna = ""
        aligned_protein = ""
        for start, end in alignment.aligned[0]:
            aligned_dna += alignment.sequences[0][start:end]
        for start, end in alignment.aligned[1]:
            aligned_protein += alignment.sequences[1][start:end]
        self.assertEqual(len(aligned_dna), 630)
        self.assertEqual(len(aligned_protein), 210)
        # Create a new alignment including the aligned sequences only:
        sequences = [aligned_dna, aligned_protein]
        coordinates = numpy.array([[0, len(aligned_dna)],
                                   [0, len(aligned_protein)]])
        alignment = Alignment(sequences, coordinates)
        # and translate the aligned DNA sequence:
        alignment.sequences[0] = alignment.sequences[0].translate()
        alignment.coordinates[0, :] //= 3
        # Now we can print the alignment:
        self.assertEqual(str(alignment), """\
MESQRWLPLEANPEVTNQFLKQLGLHPNWQCVDVYGMDPELLSMVPRPVCAVLLLFPITEKYEIFRTEEEEKTKSQGQDVTSSVYFMKQTISNACGTIGLIHAIANNKDKMHFESGSTLKKFLEESASMSPEERARYLENYDAIRVTHETSAHEGQTEAPNIDEKVDLHFIALVHVDGHLYELDAIEVCKKFMERDPDELRFNAIALSAA
||.|||||||||||||||||||||||||||.||||||||||||||||||||||||||||||||.||||||||.|||||||||||||||||||||||||||||||||||||||||||||||||||||.|||||||||||||||||||||||||||||||||.|||||||||||||||||||||||||||||||||||||||||||||||||
MEGQRWLPLEANPEVTNQFLKQLGLHPNWQFVDVYGMDPELLSMVPRPVCAVLLLFPITEKYEVFRTEEEEKIKSQGQDVTSSVYFMKQTISNACGTIGLIHAIANNKDKMHFESGSTLKKFLEESVSMSPEERARYLENYDAIRVTHETSAHEGQTEAPSIDEKVDLHFIALVHVDGHLYELDAIEVCKKFMERDPDELRFNAIALSAA
""")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
