# Copyright 2022 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Align.bigbed module."""
import unittest
import os
import warnings
from io import StringIO


from Bio.Align import Alignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import BiopythonExperimentalWarning

with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonExperimentalWarning)
    from Bio.Align import bigbed


try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align.bigbed."
    ) from None


class TestAlign_dna_rna(unittest.TestCase):

    # The bigBed file dna_rna.bb was generated using the commands
    # sort -k1,1 -k2,2n dna_rna.bed > dna_rna.sorted.bed
    # twoBitInfo hg38.2bit hg38.chrom.sizes
    # bedToBigBed dna_rna.sorted.bed hg38.chrom.sizes dna_rna.bb

    def setUp(self):
        data = {}
        records = SeqIO.parse("Blat/dna.fa", "fasta")
        for record in records:
            name, start_end = record.id.split(":")
            assert name == "chr3"
            start, end = start_end.split("-")
            start = int(start)
            end = int(end)
            sequence = str(record.seq)
            assert len(sequence) == end - start
            data[start] = sequence
        self.dna = Seq(data, length=198295559)  # hg38 chr3
        records = SeqIO.parse("Blat/rna.fa", "fasta")
        self.rna = {record.id: record.seq for record in records}

    def test_reading(self):
        """Test parsing dna_rna.bb."""
        path = "Blat/dna_rna.bb"
        alignments = bigbed.AlignmentIterator(path)
        self.assertEqual(len(alignments.targets), 1)
        self.assertEqual(alignments.targets[0].id, "chr3")
        self.assertEqual(len(alignments.targets[0]), 198295559)
        self.assertEqual(len(alignments), 4)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 1711))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr3")
        self.assertEqual(alignment.query.id, "NR_046654.1")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[42530895, 42530958, 42532020,
                              42532095, 42532563, 42532606],
                             [     181,      118,      118,
                                    43,       43,        0]])
                # fmt: on
            )
        )
        alignment.target.seq = self.dna
        alignment.query.seq = self.rna[alignment.query.id]
        self.assertTrue(
            numpy.array_equal(
                alignment.substitutions,
                # fmt: off
# flake8: noqa
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
        # The modified RNAs have gaps in their sequence. As this information is
        # not stored in a BED file, we cannot calculate the substitution matrix.
        alignment = next(alignments)
        self.assertEqual(alignment.score, 978)
        self.assertEqual(alignment.shape, (2, 1711))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr3")
        self.assertEqual(alignment.query.id, "NR_046654.1_modified")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[42530895, 42530922, 42530958, 42532020, 42532037,
                              42532039, 42532095, 42532563, 42532606],
                             [     179,      152,      116,      116,       99,
                                    99,       43,       43,        0]])
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 5407))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr3")
        self.assertEqual(alignment.query.id, "NR_111921.1")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[48663767, 48663813, 48665640,
                              48665722, 48669098, 48669174],
                             [       0,       46,       46,
                                   128,      128,      204]])
                # fmt: on
            )
        )
        alignment.target.seq = self.dna
        alignment.query.seq = self.rna[alignment.query.id]
        self.assertTrue(
            numpy.array_equal(
                alignment.substitutions,
                # fmt: off
# flake8: noqa
            numpy.array([[53.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0., 35.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0., 50.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  0., 27.,  0.,  0.,  0.,  0.],
                         [ 9.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  7.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0., 16.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  0.,  7.,  0.,  0.,  0.,  0.],
                        ])
            )
        )
        self.assertEqual(alignment.substitutions.alphabet, "ACGTacgt")
        alignment = next(alignments)
        self.assertEqual(alignment.score, 972)
        self.assertEqual(alignment.shape, (2, 5407))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr3")
        self.assertEqual(alignment.query.id, "NR_111921.1_modified")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[48663767, 48663795, 48663796, 48663813, 48665640,
                              48665716, 48665722, 48669098, 48669174],
                             [       0,       28,       28,       45,       45,
                                   121,      127,      127,      203]])
                # fmt: on
            )
        )
        self.assertRaises(StopIteration, next, alignments)


class TestAlign_dna(unittest.TestCase):
    def test_reading_psl_34_001(self):
        """Test parsing psl_34_001.bb."""

        # The bigBed file psl_34_001.bb was generated using the commands
        # sort -k1,1 -k2,2n psl_34_001.bed > psl_34_001.sorted.bed
        # twoBitInfo hg19.2bit hg19.chrom.sizes
        # bedToBigBed psl_34_001.sorted.bed hg19.chrom.sizes psl_34_001.bb

        path = "Blat/psl_34_001.bb"
        alignments = bigbed.AlignmentIterator(path)
        self.assertEqual(len(alignments.targets), 10)
        self.assertEqual(alignments.targets[0].id, "chr1")
        self.assertEqual(len(alignments.targets[0]), 249250621)
        self.assertEqual(alignments.targets[1].id, "chr10")
        self.assertEqual(len(alignments.targets[1]), 135534747)
        self.assertEqual(alignments.targets[2].id, "chr13")
        self.assertEqual(len(alignments.targets[2]), 115169878)
        self.assertEqual(alignments.targets[3].id, "chr18")
        self.assertEqual(len(alignments.targets[3]), 78077248)
        self.assertEqual(alignments.targets[4].id, "chr19")
        self.assertEqual(len(alignments.targets[4]), 59128983)
        self.assertEqual(alignments.targets[5].id, "chr2")
        self.assertEqual(len(alignments.targets[5]), 243199373)
        self.assertEqual(alignments.targets[6].id, "chr22")
        self.assertEqual(len(alignments.targets[6]), 51304566)
        self.assertEqual(alignments.targets[7].id, "chr4")
        self.assertEqual(len(alignments.targets[7]), 191154276)
        self.assertEqual(alignments.targets[8].id, "chr8")
        self.assertEqual(len(alignments.targets[8]), 146364022)
        self.assertEqual(alignments.targets[9].id, "chr9")
        self.assertEqual(len(alignments.targets[9]), 141213431)
        self.assertEqual(len(alignments), 22)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 50))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[1207056, 1207106],
                             [      0,      50]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 33))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[10271783, 10271816],
                             [       0,       33]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 946)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[39368490, 39368526],
                             [      36,        0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 824)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[61700837, 61700871],
                             [       0,       34]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 942)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[220325687, 220325721],
                             [       34,         0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 834)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[99388555, 99388591],
                             [      36,        0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 920)
        self.assertEqual(alignment.shape, (2, 25))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[112178171, 112178196],
                             [       25,         0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 912)
        self.assertEqual(alignment.shape, (2, 51))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[52759147, 52759154, 52759160, 52759198],
                             [       0,        7,        7,       45]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[23891310, 23891349],
                             [       0,       39]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 930)
        self.assertEqual(alignment.shape, (2, 28))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[43252217, 43252245],
                             [       0,       28]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 848)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[553742, 553781],
                             [    39,      0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 890)
        self.assertEqual(alignment.shape, (2, 170))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[35483340, 35483365, 35483499, 35483510],
                             [       0,       25,       25,       36]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[54017130, 54017169],
                             [      39,        0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 17))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[53575980, 53575997],
                             [      17,        0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 946)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[120641740, 120641776],
                             [       36,         0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 682)
        self.assertEqual(alignment.shape, (2, 44))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[183925984, 183925990, 183926028],
                             [        0,         6,        44]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 834)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[42144400, 42144436],
                             [       0,       36]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 892)
        self.assertEqual(alignment.shape, (2, 37))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[48997405, 48997442],
                             [      37,        0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 572)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[37558157, 37558167, 37558173, 37558191],
                             [      28,       18,       18,        0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 16))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[61646095, 61646111],
                             [       0,       16]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr8")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[95160479, 95160520],
                             [       0,       41]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 854)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr9")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[85737865, 85737906],
                             [       0,       41]]),
                # fmt: on
            )
        )
        self.assertRaises(StopIteration, next, alignments)

    def test_reading_psl_34_003(self):
        """Test parsing psl_34_003.bb."""

        # The bigBed file psl_34_003.bb was generated using the commands
        # sort -k1,1 -k2,2n psl_34_003.bed > psl_34_003.sorted.bed
        # twoBitInfo hg19.2bit hg19.chrom.sizes
        # bedToBigBed psl_34_003.sorted.bed hg19.chrom.sizes psl_34_003.bb

        path = "Blat/psl_34_003.bb"
        alignments = bigbed.AlignmentIterator(path)
        self.assertEqual(len(alignments.targets), 3)
        self.assertEqual(alignments.targets[0].id, "chr1")
        self.assertEqual(len(alignments.targets[0]), 249250621)
        self.assertEqual(alignments.targets[1].id, "chr2")
        self.assertEqual(len(alignments.targets[1]), 243199373)
        self.assertEqual(alignments.targets[2].id, "chr4")
        self.assertEqual(len(alignments.targets[2]), 191154276)
        self.assertEqual(len(alignments), 3)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 33))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[10271783, 10271816],
                             [       0,       33]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 17))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[53575980, 53575997],
                             [      17,        0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 16))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[61646095, 61646111],
                             [       0,       16]]),
                # fmt: on
            )
        )
        self.assertRaises(StopIteration, next, alignments)

    def test_reading_psl_34_004(self):
        """Test parsing psl_34_004.bb."""

        # The bigBed file psl_34_004.bb was generated using the commands
        # sort -k1,1 -k2,2n psl_34_004.bed > psl_34_004.sorted.bed
        # twoBitInfo hg19.2bit hg19.chrom.sizes
        # bedToBigBed psl_34_004.sorted.bed hg19.chrom.sizes psl_34_004.bb

        path = "Blat/psl_34_004.bb"
        alignments = bigbed.AlignmentIterator(path)
        self.assertEqual(len(alignments.targets), 10)
        self.assertEqual(alignments.targets[0].id, "chr1")
        self.assertEqual(len(alignments.targets[0]), 249250621)
        self.assertEqual(alignments.targets[1].id, "chr10")
        self.assertEqual(len(alignments.targets[1]), 135534747)
        self.assertEqual(alignments.targets[2].id, "chr13")
        self.assertEqual(len(alignments.targets[2]), 115169878)
        self.assertEqual(alignments.targets[3].id, "chr18")
        self.assertEqual(len(alignments.targets[3]), 78077248)
        self.assertEqual(alignments.targets[4].id, "chr19")
        self.assertEqual(len(alignments.targets[4]), 59128983)
        self.assertEqual(alignments.targets[5].id, "chr2")
        self.assertEqual(len(alignments.targets[5]), 243199373)
        self.assertEqual(alignments.targets[6].id, "chr22")
        self.assertEqual(len(alignments.targets[6]), 51304566)
        self.assertEqual(alignments.targets[7].id, "chr4")
        self.assertEqual(len(alignments.targets[7]), 191154276)
        self.assertEqual(alignments.targets[8].id, "chr8")
        self.assertEqual(len(alignments.targets[8]), 146364022)
        self.assertEqual(alignments.targets[9].id, "chr9")
        self.assertEqual(len(alignments.targets[9]), 141213431)
        self.assertEqual(len(alignments), 19)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 50))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[1207056, 1207106],
                             [      0,      50]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 946)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[39368490, 39368526],
                             [      36,        0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 824)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[61700837, 61700871],
                             [       0,       34]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 942)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[220325687, 220325721],
                             [       34,         0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 834)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[99388555, 99388591],
                             [      36,        0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 920)
        self.assertEqual(alignment.shape, (2, 25))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[112178171, 112178196],
                             [       25,         0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 912)
        self.assertEqual(alignment.shape, (2, 51))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[52759147, 52759154, 52759160, 52759198],
                             [       0,        7,        7,       45]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[23891310, 23891349],
                             [       0,       39]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 930)
        self.assertEqual(alignment.shape, (2, 28))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[43252217, 43252245],
                             [       0,       28]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 848)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[553742, 553781],
                             [    39,      0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 890)
        self.assertEqual(alignment.shape, (2, 170))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[35483340, 35483365, 35483499, 35483510],
                             [       0,       25,       25,       36]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[54017130, 54017169],
                             [      39,        0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 946)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[120641740, 120641776],
                             [       36,         0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 682)
        self.assertEqual(alignment.shape, (2, 44))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[183925984, 183925990, 183926028],
                             [        0,         6,        44]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 834)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[42144400, 42144436],
                             [       0,       36]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 892)
        self.assertEqual(alignment.shape, (2, 37))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[48997405, 48997442],
                             [      37,        0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 572)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[37558157, 37558167, 37558173, 37558191],
                             [      28,       18,       18,        0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr8")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[95160479, 95160520],
                             [       0,       41]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 854)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr9")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[85737865, 85737906],
                             [       0,       41]]),
                # fmt: on
            )
        )
        self.assertRaises(StopIteration, next, alignments)

    def test_reading_psl_34_005(self):
        """Test parsing psl_34_005.bb."""
        path = "Blat/psl_34_005.bb"
        alignments = bigbed.AlignmentIterator(path)
        self.assertEqual(len(alignments.targets), 10)
        self.assertEqual(alignments.targets[0].id, "chr1")
        self.assertEqual(len(alignments.targets[0]), 249250621)
        self.assertEqual(alignments.targets[1].id, "chr10")
        self.assertEqual(len(alignments.targets[1]), 135534747)
        self.assertEqual(alignments.targets[2].id, "chr13")
        self.assertEqual(len(alignments.targets[2]), 115169878)
        self.assertEqual(alignments.targets[3].id, "chr18")
        self.assertEqual(len(alignments.targets[3]), 78077248)
        self.assertEqual(alignments.targets[4].id, "chr19")
        self.assertEqual(len(alignments.targets[4]), 59128983)
        self.assertEqual(alignments.targets[5].id, "chr2")
        self.assertEqual(len(alignments.targets[5]), 243199373)
        self.assertEqual(alignments.targets[6].id, "chr22")
        self.assertEqual(len(alignments.targets[6]), 51304566)
        self.assertEqual(alignments.targets[7].id, "chr4")
        self.assertEqual(len(alignments.targets[7]), 191154276)
        self.assertEqual(alignments.targets[8].id, "chr8")
        self.assertEqual(len(alignments.targets[8]), 146364022)
        self.assertEqual(alignments.targets[9].id, "chr9")
        self.assertEqual(len(alignments.targets[9]), 141213431)
        self.assertEqual(len(alignments), 22)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 50))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[1207056, 1207106],
                             [      0,      50]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 33))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[10271783, 10271816],
                             [       0,       33]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 946)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
                numpy.array([[39368490, 39368526],
                             [      36,        0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 824)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[61700837, 61700871],
                             [       0,       34]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 942)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[220325687, 220325721],
                             [       34,         0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 834)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[99388555, 99388591],
                             [      36,        0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 920)
        self.assertEqual(alignment.shape, (2, 25))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[112178171, 112178196],
                             [       25,         0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 912)
        self.assertEqual(alignment.shape, (2, 51))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[52759147, 52759154, 52759160, 52759198],
                             [       0,        7,        7,       45]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[23891310, 23891349],
                             [       0,       39]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 930)
        self.assertEqual(alignment.shape, (2, 28))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[43252217, 43252245],
                             [       0,       28]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 848)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[553742, 553781],
                             [    39,      0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 890)
        self.assertEqual(alignment.shape, (2, 170))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[35483340, 35483365, 35483499, 35483510],
                             [       0,       25,       25,       36]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[54017130, 54017169],
                             [      39,        0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 17))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[53575980, 53575997],
                             [      17,        0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 946)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[120641740, 120641776],
                             [       36,         0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 682)
        self.assertEqual(alignment.shape, (2, 44))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[183925984, 183925990, 183926028],
                             [        0,         6,        44]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 834)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[42144400, 42144436],
                             [       0,       36]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 892)
        self.assertEqual(alignment.shape, (2, 37))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[48997405, 48997442],
                             [      37,        0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 572)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[37558157, 37558167, 37558173, 37558191],
                             [      28,       18,       18,        0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 16))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[61646095, 61646111],
                             [       0,       16]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr8")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[95160479, 95160520],
                             [       0,       41]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 854)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr9")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[85737865, 85737906],
                             [       0,       41]]),
                # fmt: on
            )
        )
        self.assertRaises(StopIteration, next, alignments)


class TestAlign_dnax_prot(unittest.TestCase):
    def test_reading_psl_35_001(self):
        """Test parsing psl_35_001.bb."""

        # The bigBed file psl_35_001.bb was generated using the commands
        # sort -k1,1 -k2,2n psl_35_001.bed > psl_35_001.sorted.bed
        # twoBitInfo hg38.2bit hg38.chrom.sizes
        # bedToBigBed psl_35_001.sorted.bed hg38.chrom.sizes psl_35_001.bb

        path = "Blat/psl_35_001.bb"
        alignments = bigbed.AlignmentIterator(path)
        self.assertEqual(len(alignments.targets), 2)
        self.assertEqual(alignments.targets[0].id, "chr13")
        self.assertEqual(len(alignments.targets[0]), 114364328)
        self.assertEqual(alignments.targets[1].id, "chr4")
        self.assertEqual(len(alignments.targets[1]), 190214555)
        self.assertEqual(len(alignments), 8)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 986)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[75549820, 75549865, 75567225, 75567312],
                             [       0,       45,       45,      132]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[75560749, 75560881],
                             [       0,      132]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[75566694, 75566850],
                             [       0,      156]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[75569459, 75569507],
                             [       0,       48]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[75594914, 75594989],
                             [       0,       75]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[75604767, 75604827, 75605728, 75605809],
                             [       0,       60,       60,      141]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 166)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[41257605, 41257731, 41263227, 41263290],
                             [       0,      126,      126,      189]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 530)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[41260685, 41260787],
                             [       0,      102]]),
                # fmt: on
            )
        )
        self.assertRaises(StopIteration, next, alignments)


class TestAlign_bed12(unittest.TestCase):

    # The bigBed files were generated using the commands
    # twoBitInfo hg19.2bit hg19.chrom.sizes
    # bedToBigBed bed4.bed hg19.chrom.sizes bed4.bb
    # bedToBigBed bed5.bed hg19.chrom.sizes bed5.bb
    # bedToBigBed bed6.bed hg19.chrom.sizes bed6.bb
    # bedToBigBed bed7.bed hg19.chrom.sizes bed7.bb
    # bedToBigBed bed8.bed hg19.chrom.sizes bed8.bb
    # bedToBigBed bed9.bed hg19.chrom.sizes bed9.bb
    # bedToBigBed bed12.bed hg19.chrom.sizes bed12.bb

    def test_reading(self):
        """Test parsing alignments in file formats BED3 through BED12."""
        for bedN in (3, 4, 5, 6, 7, 8, 9, 12):
            filename = "bed%d.bb" % bedN
            path = os.path.join("Blat", filename)
            alignments = bigbed.AlignmentIterator(path)
            self.assertEqual(len(alignments), 2)
            alignment = next(alignments)
            if bedN >= 5:
                self.assertEqual(alignment.score, 960, msg=filename)
            self.assertEqual(alignment.shape, (2, 4000), msg=filename)
            self.assertLess(
                alignment.coordinates[0, 0], alignment.coordinates[0, -1], msg=filename
            )
            self.assertLess(
                alignment.coordinates[1, 0], alignment.coordinates[1, -1], msg=filename
            )
            self.assertEqual(len(alignment), 2, msg=filename)
            self.assertIs(alignment.sequences[0], alignment.target, msg=filename)
            self.assertIs(alignment.sequences[1], alignment.query, msg=filename)
            self.assertEqual(alignment.target.id, "chr22", msg=filename)
            if bedN >= 4:
                self.assertEqual(alignment.query.id, "mRNA1", msg=filename)
            else:
                self.assertIsNone(alignment.query.id, msg=filename)
            if bedN == 12:
                self.assertTrue(
                    numpy.array_equal(
                        alignment.coordinates,
                        # fmt: off
# flake8: noqa
                        numpy.array([[1000, 1567, 4512, 5000],
                                     [   0,  567,  567, 1055]]),
                        # fmt: on
                    ),
                    msg=filename,
                )
            else:
                self.assertTrue(
                    numpy.array_equal(
                        alignment.coordinates,
                        numpy.array([[1000, 5000], [0, 4000]]),
                    ),
                    msg=filename,
                )
            if bedN >= 7:
                self.assertEqual(alignment.thickStart, 1200, msg=filename)
            if bedN >= 8:
                self.assertEqual(alignment.thickEnd, 4900, msg=filename)
            if bedN >= 9:
                self.assertEqual(alignment.itemRgb, "255,0,0", msg=filename)
            alignment = next(alignments)
            if bedN >= 5:
                self.assertEqual(alignment.score, 900, msg=filename)
            self.assertEqual(alignment.shape, (2, 4000), msg=filename)
            self.assertLess(
                alignment.coordinates[0, 0], alignment.coordinates[0, -1], msg=filename
            )
            if bedN >= 6:
                self.assertGreater(
                    alignment.coordinates[1, 0],
                    alignment.coordinates[1, -1],
                    msg=filename,
                )
            else:
                self.assertLess(
                    alignment.coordinates[1, 0],
                    alignment.coordinates[1, -1],
                    msg=filename,
                )
            self.assertEqual(len(alignment), 2, msg=filename)
            self.assertIs(alignment.sequences[0], alignment.target, msg=filename)
            self.assertIs(alignment.sequences[1], alignment.query, msg=filename)
            self.assertEqual(alignment.target.id, "chr22", msg=filename)
            if bedN >= 4:
                self.assertEqual(alignment.query.id, "mRNA2", msg=filename)
            else:
                self.assertIsNone(alignment.query.id, msg=filename)
            if bedN == 12:
                self.assertTrue(
                    numpy.array_equal(
                        alignment.coordinates,
                        # fmt: off
# flake8: noqa
                        numpy.array([[2000, 2433, 5601, 6000],
                                     [ 832,  399,  399,    0]])
                        # fmt: on
                    ),
                    msg=filename,
                )
            elif bedN >= 6:
                self.assertTrue(
                    numpy.array_equal(
                        alignment.coordinates,
                        numpy.array([[2000, 6000], [4000, 0]]),
                    ),
                    msg=filename,
                )
            else:
                self.assertTrue(
                    numpy.array_equal(
                        alignment.coordinates,
                        numpy.array([[2000, 6000], [0, 4000]]),
                    ),
                    msg=filename,
                )
            if bedN >= 7:
                self.assertEqual(alignment.thickStart, 2300, msg=filename)
            if bedN >= 8:
                self.assertEqual(alignment.thickEnd, 5960, msg=filename)
            if bedN >= 9:
                self.assertEqual(alignment.itemRgb, "0,255,0", msg=filename)
            with self.assertRaises(StopIteration) as cm:
                next(alignments)
                self.fail(f"More than two alignments reported in {filename}")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
