# Copyright 2022 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Align.bigbed module."""
import unittest
import os
import warnings
from io import StringIO


from Bio import Align
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
        alignments = Align.parse(path, "bigbed")
        self.assertEqual(
            str(alignments.declaration),
            """\
table bed
"Browser Extensible Data"
(
   string          chrom;          "Reference sequence chromosome or scaffold"
   uint            chromStart;     "Start position in chromosome"
   uint            chromEnd;       "End position in chromosome"
   string          name;           "Name of item."
   uint            score;          "Score (0-1000)"
   char[1]         strand;         "+ or - for strand"
   uint            thickStart;     "Start of where display should be thick (start codon)"
   uint            thickEnd;       "End of where display should be thick (stop codon)"
   uint            reserved;       "Used as itemRgb as of 2004-11-22"
   int             blockCount;     "Number of blocks"
   int[blockCount] blockSizes;     "Comma separated list of block sizes"
   int[blockCount] chromStarts;    "Start positions relative to chromStart"
)
""",
        )
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
        alignments = Align.parse(path, "bigbed")
        self.assertEqual(
            str(alignments.declaration),
            """\
table bed
"Browser Extensible Data"
(
   string          chrom;          "Reference sequence chromosome or scaffold"
   uint            chromStart;     "Start position in chromosome"
   uint            chromEnd;       "End position in chromosome"
   string          name;           "Name of item."
   uint            score;          "Score (0-1000)"
   char[1]         strand;         "+ or - for strand"
   uint            thickStart;     "Start of where display should be thick (start codon)"
   uint            thickEnd;       "End of where display should be thick (stop codon)"
   uint            reserved;       "Used as itemRgb as of 2004-11-22"
   int             blockCount;     "Number of blocks"
   int[blockCount] blockSizes;     "Comma separated list of block sizes"
   int[blockCount] chromStarts;    "Start positions relative to chromStart"
)
""",
        )
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
        alignments = Align.parse(path, "bigbed")
        self.assertEqual(
            str(alignments.declaration),
            """\
table bed
"Browser Extensible Data"
(
   string          chrom;          "Reference sequence chromosome or scaffold"
   uint            chromStart;     "Start position in chromosome"
   uint            chromEnd;       "End position in chromosome"
   string          name;           "Name of item."
   uint            score;          "Score (0-1000)"
   char[1]         strand;         "+ or - for strand"
   uint            thickStart;     "Start of where display should be thick (start codon)"
   uint            thickEnd;       "End of where display should be thick (stop codon)"
   uint            reserved;       "Used as itemRgb as of 2004-11-22"
   int             blockCount;     "Number of blocks"
   int[blockCount] blockSizes;     "Comma separated list of block sizes"
   int[blockCount] chromStarts;    "Start positions relative to chromStart"
)
""",
        )
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
        alignments = Align.parse(path, "bigbed")
        self.assertEqual(
            str(alignments.declaration),
            """\
table bed
"Browser Extensible Data"
(
   string          chrom;          "Reference sequence chromosome or scaffold"
   uint            chromStart;     "Start position in chromosome"
   uint            chromEnd;       "End position in chromosome"
   string          name;           "Name of item."
   uint            score;          "Score (0-1000)"
   char[1]         strand;         "+ or - for strand"
   uint            thickStart;     "Start of where display should be thick (start codon)"
   uint            thickEnd;       "End of where display should be thick (stop codon)"
   uint            reserved;       "Used as itemRgb as of 2004-11-22"
   int             blockCount;     "Number of blocks"
   int[blockCount] blockSizes;     "Comma separated list of block sizes"
   int[blockCount] chromStarts;    "Start positions relative to chromStart"
)
""",
        )
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
        alignments = Align.parse(path, "bigbed")
        self.assertEqual(
            str(alignments.declaration),
            """\
table bed
"Browser Extensible Data"
(
   string          chrom;          "Reference sequence chromosome or scaffold"
   uint            chromStart;     "Start position in chromosome"
   uint            chromEnd;       "End position in chromosome"
   string          name;           "Name of item."
   uint            score;          "Score (0-1000)"
   char[1]         strand;         "+ or - for strand"
   uint            thickStart;     "Start of where display should be thick (start codon)"
   uint            thickEnd;       "End of where display should be thick (stop codon)"
   uint            reserved;       "Used as itemRgb as of 2004-11-22"
   int             blockCount;     "Number of blocks"
   int[blockCount] blockSizes;     "Comma separated list of block sizes"
   int[blockCount] chromStarts;    "Start positions relative to chromStart"
)
""",
        )
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
        alignments = Align.parse(path, "bigbed")
        self.assertEqual(
            str(alignments.declaration),
            """\
table bed
"Browser Extensible Data"
(
   string          chrom;          "Reference sequence chromosome or scaffold"
   uint            chromStart;     "Start position in chromosome"
   uint            chromEnd;       "End position in chromosome"
   string          name;           "Name of item."
   uint            score;          "Score (0-1000)"
   char[1]         strand;         "+ or - for strand"
   uint            thickStart;     "Start of where display should be thick (start codon)"
   uint            thickEnd;       "End of where display should be thick (stop codon)"
   uint            reserved;       "Used as itemRgb as of 2004-11-22"
   int             blockCount;     "Number of blocks"
   int[blockCount] blockSizes;     "Comma separated list of block sizes"
   int[blockCount] chromStarts;    "Start positions relative to chromStart"
)
""",
        )
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

    def check_autosql(self, declaration, bedN):
        if bedN == 3:
            self.assertEqual(
                str(declaration),
                """\
table bed
"Browser Extensible Data"
(
   string chrom;         "Reference sequence chromosome or scaffold"
   uint   chromStart;    "Start position in chromosome"
   uint   chromEnd;      "End position in chromosome"
)
""",
            )
        elif bedN == 4:
            self.assertEqual(
                str(declaration),
                """\
table bed
"Browser Extensible Data"
(
   string chrom;         "Reference sequence chromosome or scaffold"
   uint   chromStart;    "Start position in chromosome"
   uint   chromEnd;      "End position in chromosome"
   string name;          "Name of item."
)
""",
            )
        elif bedN == 5:
            self.assertEqual(
                str(declaration),
                """\
table bed
"Browser Extensible Data"
(
   string chrom;         "Reference sequence chromosome or scaffold"
   uint   chromStart;    "Start position in chromosome"
   uint   chromEnd;      "End position in chromosome"
   string name;          "Name of item."
   uint   score;         "Score (0-1000)"
)
""",
            )
        elif bedN == 6:
            self.assertEqual(
                str(declaration),
                """\
table bed
"Browser Extensible Data"
(
   string  chrom;         "Reference sequence chromosome or scaffold"
   uint    chromStart;    "Start position in chromosome"
   uint    chromEnd;      "End position in chromosome"
   string  name;          "Name of item."
   uint    score;         "Score (0-1000)"
   char[1] strand;        "+ or - for strand"
)
""",
            )
        elif bedN == 7:
            self.assertEqual(
                str(declaration),
                """\
table bed
"Browser Extensible Data"
(
   string  chrom;         "Reference sequence chromosome or scaffold"
   uint    chromStart;    "Start position in chromosome"
   uint    chromEnd;      "End position in chromosome"
   string  name;          "Name of item."
   uint    score;         "Score (0-1000)"
   char[1] strand;        "+ or - for strand"
   uint    thickStart;    "Start of where display should be thick (start codon)"
)
""",
            )
        elif bedN == 8:
            self.assertEqual(
                str(declaration),
                """\
table bed
"Browser Extensible Data"
(
   string  chrom;         "Reference sequence chromosome or scaffold"
   uint    chromStart;    "Start position in chromosome"
   uint    chromEnd;      "End position in chromosome"
   string  name;          "Name of item."
   uint    score;         "Score (0-1000)"
   char[1] strand;        "+ or - for strand"
   uint    thickStart;    "Start of where display should be thick (start codon)"
   uint    thickEnd;      "End of where display should be thick (stop codon)"
)
""",
            )
        elif bedN == 9:
            self.assertEqual(
                str(declaration),
                """\
table bed
"Browser Extensible Data"
(
   string  chrom;         "Reference sequence chromosome or scaffold"
   uint    chromStart;    "Start position in chromosome"
   uint    chromEnd;      "End position in chromosome"
   string  name;          "Name of item."
   uint    score;         "Score (0-1000)"
   char[1] strand;        "+ or - for strand"
   uint    thickStart;    "Start of where display should be thick (start codon)"
   uint    thickEnd;      "End of where display should be thick (stop codon)"
   uint    reserved;      "Used as itemRgb as of 2004-11-22"
)
""",
            )
        elif bedN == 12:
            self.assertEqual(
                str(declaration),
                """\
table bed
"Browser Extensible Data"
(
   string          chrom;          "Reference sequence chromosome or scaffold"
   uint            chromStart;     "Start position in chromosome"
   uint            chromEnd;       "End position in chromosome"
   string          name;           "Name of item."
   uint            score;          "Score (0-1000)"
   char[1]         strand;         "+ or - for strand"
   uint            thickStart;     "Start of where display should be thick (start codon)"
   uint            thickEnd;       "End of where display should be thick (stop codon)"
   uint            reserved;       "Used as itemRgb as of 2004-11-22"
   int             blockCount;     "Number of blocks"
   int[blockCount] blockSizes;     "Comma separated list of block sizes"
   int[blockCount] chromStarts;    "Start positions relative to chromStart"
)
""",
            )

    def test_reading(self):
        """Test parsing alignments in file formats BED3 through BED12."""
        for bedN in (3, 4, 5, 6, 7, 8, 9, 12):
            filename = "bed%d.bb" % bedN
            path = os.path.join("Blat", filename)
            alignments = Align.parse(path, "bigbed")
            self.check_autosql(alignments.declaration, bedN)
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


class TestAlign_extended_bed(unittest.TestCase):

    # The bigBed file bigbed_extended.bb is a BED9+2 file, with nine predefined
    # BED fields and 2 extra (custom) fields. It was created by running
    #
    # bedToBigBed -as=bedExample2.as -type=bed9+2 -extraIndex=name,geneSymbol bedExample2.bed hg18.chrom.sizes bigbed_extended.bb
    #
    # where bedExample2.bed contains 10 lines selected from the example BED9+2
    # file bedExample2.bed, and bedExample2.as the associated autoSql file from
    # UCSC declaring the nine predefined BED fields and the two extra fields.

    def test_reading(self):
        """Test parsing bigbed_extended.bb."""
        path = "Blat/bigbed_extended.bb"
        alignments = Align.parse(path, "bigbed")
        self.maxDiff = None
        self.assertEqual(
            str(alignments.declaration),
            """\
table hg18KGchr7
"UCSC Genes for chr7 with color plus GeneSymbol and SwissProtID"
(
   string  chrom;         "Reference sequence chromosome or scaffold"
   uint    chromStart;    "Start position of feature on chromosome"
   uint    chromEnd;      "End position of feature on chromosome"
   string  name;          "Name of gene"
   uint    score;         "Score"
   char[1] strand;        "+ or - for strand"
   uint    thickStart;    "Coding region start"
   uint    thickEnd;      "Coding region end"
   uint    reserved;      "Green on + strand, Red on - strand"
   string  geneSymbol;    "Gene Symbol"
   string  spID;          "SWISS-PROT protein Accession number"
)
""",
        )
        self.assertEqual(len(alignments.targets), 1)
        self.assertEqual(alignments.targets[0].id, "chr7")
        self.assertEqual(len(alignments.targets[0]), 158821424)
        self.assertEqual(len(alignments), 10)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.thickStart, 60328)
        self.assertEqual(alignment.thickEnd, 60328)
        self.assertEqual(alignment.itemRgb, "255,0,0")
        self.assertEqual(alignment.shape, (2, 1241))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr7")
        self.assertEqual(alignment.query.id, "uc010krx.1")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[60328, 61569], [1241, 0]])
            )
        )
        self.assertEqual(alignment.annotations["geneSymbol"], ".")
        self.assertEqual(alignment.annotations["spID"], "PDGFA")
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.thickStart, 506606)
        self.assertEqual(alignment.thickEnd, 525164)
        self.assertEqual(alignment.itemRgb, "255,0,0")
        self.assertEqual(alignment.shape, (2, 22585))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr7")
        self.assertEqual(alignment.query.id, "uc003sir.1")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[503422, 526007], [22585, 0]])
            )
        )
        self.assertEqual(alignment.annotations["geneSymbol"], "PDGFA")
        self.assertEqual(alignment.annotations["spID"], "P04085")
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.thickStart, 504726)
        self.assertEqual(alignment.thickEnd, 525164)
        self.assertEqual(alignment.itemRgb, "255,0,0")
        self.assertEqual(alignment.shape, (2, 22585))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr7")
        self.assertEqual(alignment.query.id, "uc003sis.1")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[503422, 526007], [22585, 0]])
            )
        )
        self.assertEqual(alignment.annotations["geneSymbol"], "PDGFA")
        self.assertEqual(alignment.annotations["spID"], "P04085-2")
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.thickStart, 507195)
        self.assertEqual(alignment.thickEnd, 518820)
        self.assertEqual(alignment.itemRgb, "255,0,0")
        self.assertEqual(alignment.shape, (2, 12690))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr7")
        self.assertEqual(alignment.query.id, "uc003sit.1")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[506940, 519630], [12690, 0]])
            )
        )
        self.assertEqual(alignment.annotations["geneSymbol"], "PDGFA")
        self.assertEqual(alignment.annotations["spID"], "Q32M96")
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.thickStart, 556592)
        self.assertEqual(alignment.thickEnd, 717668)
        self.assertEqual(alignment.itemRgb, "255,0,0")
        self.assertEqual(alignment.shape, (2, 162747))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr7")
        self.assertEqual(alignment.query.id, "uc003siu.1")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[555912, 718659], [162747, 0]])
            )
        )
        self.assertEqual(alignment.annotations["geneSymbol"], "PRKAR1B")
        self.assertEqual(alignment.annotations["spID"], "Q8N422")
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.thickStart, 556592)
        self.assertEqual(alignment.thickEnd, 717668)
        self.assertEqual(alignment.itemRgb, "255,0,0")
        self.assertEqual(alignment.shape, (2, 163357))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr7")
        self.assertEqual(alignment.query.id, "uc003siv.1")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[555912, 719269], [163357, 0]])
            )
        )
        self.assertEqual(alignment.annotations["geneSymbol"], "PRKAR1B")
        self.assertEqual(alignment.annotations["spID"], "Q8N422")
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.thickStart, 556592)
        self.assertEqual(alignment.thickEnd, 717668)
        self.assertEqual(alignment.itemRgb, "255,0,0")
        self.assertEqual(alignment.shape, (2, 177901))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr7")
        self.assertEqual(alignment.query.id, "uc003siw.1")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[555912, 733813], [177901, 0]])
            )
        )
        self.assertEqual(alignment.annotations["geneSymbol"], "PRKAR1B")
        self.assertEqual(alignment.annotations["spID"], "Q8N422")
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.thickStart, 585418)
        self.assertEqual(alignment.thickEnd, 585418)
        self.assertEqual(alignment.itemRgb, "255,0,0")
        self.assertEqual(alignment.shape, (2, 22329))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr7")
        self.assertEqual(alignment.query.id, "uc003six.1")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[585418, 607747], [22329, 0]])
            )
        )
        self.assertEqual(alignment.annotations["geneSymbol"], ".")
        self.assertEqual(alignment.annotations["spID"], "PRKAR1B")
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.thickStart, 733217)
        self.assertEqual(alignment.thickEnd, 791816)
        self.assertEqual(alignment.itemRgb, "0,255,0")
        self.assertEqual(alignment.shape, (2, 59779))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr7")
        self.assertEqual(alignment.query.id, "uc003siz.2")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[732863, 792642], [0, 59779]])
            )
        )
        self.assertEqual(alignment.annotations["geneSymbol"], ".")
        self.assertEqual(alignment.annotations["spID"], "DKFZp762F1415")
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.thickStart, 732883)
        self.assertEqual(alignment.thickEnd, 791816)
        self.assertEqual(alignment.itemRgb, "0,255,0")
        self.assertEqual(alignment.shape, (2, 59779))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr7")
        self.assertEqual(alignment.query.id, "uc010krz.1")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[732863, 792642], [0, 59779]])
            )
        )
        self.assertEqual(alignment.annotations["geneSymbol"], "HEATR2")
        self.assertEqual(alignment.annotations["spID"], "Q86Y56")


class TestAlign_searching(unittest.TestCase):

    # The bigBed file bigbedtest.bb contains the following data:
    # chr1     10     100     name1   1       +
    # chr1     29      39     name2   2       -
    # chr1    200     300     name3   3       +
    # chr2     50      50     name4   6       +
    # chr2    100     110     name5   4       +
    # chr2    200     210     name6   5       +
    # chr2    220     220     name7   6       +
    # chr3      0       0     name8   7       -

    def test_search_chromosome(self):
        path = "Blat/bigbedtest.bb"
        alignments = Align.parse(path, "bigbed")
        self.assertEqual(
            str(alignments.declaration),
            """\
table bed
"Browser Extensible Data"
(
   string  chrom;         "Reference sequence chromosome or scaffold"
   uint    chromStart;    "Start position in chromosome"
   uint    chromEnd;      "End position in chromosome"
   string  name;          "Name of item."
   uint    score;         "Score (0-1000)"
   char[1] strand;        "+ or - for strand"
)
""",
        )
        selected_alignments = alignments.search("chr2")
        names = [alignment.query.id for alignment in selected_alignments]
        self.assertEqual(names, ["name4", "name5", "name6", "name7"])

    def test_search_region(self):
        path = "Blat/bigbedtest.bb"
        alignments = Align.parse(path, "bigbed")
        selected_alignments = alignments.search("chr2", 105, 1000)
        names = [alignment.query.id for alignment in selected_alignments]
        self.assertEqual(names, ["name5", "name6", "name7"])
        selected_alignments = alignments.search("chr2", 110, 1000)
        names = [alignment.query.id for alignment in selected_alignments]
        self.assertEqual(names, ["name6", "name7"])
        selected_alignments = alignments.search("chr2", 40, 50)
        names = [alignment.query.id for alignment in selected_alignments]
        self.assertEqual(names, ["name4"])
        selected_alignments = alignments.search("chr2", 50, 50)
        names = [alignment.query.id for alignment in selected_alignments]
        self.assertEqual(names, ["name4"])
        selected_alignments = alignments.search("chr2", 50, 200)
        names = [alignment.query.id for alignment in selected_alignments]
        self.assertEqual(names, ["name4", "name5"])
        selected_alignments = alignments.search("chr2", 200, 220)
        names = [alignment.query.id for alignment in selected_alignments]
        self.assertEqual(names, ["name6", "name7"])
        selected_alignments = alignments.search("chr2", 220, 220)
        names = [alignment.query.id for alignment in selected_alignments]
        self.assertEqual(names, ["name7"])

    def test_search_position(self):
        path = "Blat/bigbedtest.bb"
        alignments = Align.parse(path, "bigbed")
        selected_alignments = alignments.search("chr1", 250)
        names = [alignment.query.id for alignment in selected_alignments]
        self.assertEqual(names, ["name3"])

    def test_three_iterators(self):
        """Create three iterators and use them concurrently."""
        path = "Blat/bigbedtest.bb"
        alignments1 = Align.parse(path, "bigbed")
        alignments2 = alignments1.search("chr2")
        alignments3 = alignments1.search("chr2", 110, 1000)
        alignment1 = next(alignments1)
        self.assertEqual(alignment1.query.id, "name1")
        alignment1 = next(alignments1)
        self.assertEqual(alignment1.query.id, "name2")
        alignment2 = next(alignments2)
        self.assertEqual(alignment2.query.id, "name4")
        alignment2 = next(alignments2)
        self.assertEqual(alignment2.query.id, "name5")
        alignment2 = next(alignments2)
        self.assertEqual(alignment2.query.id, "name6")
        alignment3 = next(alignments3)
        self.assertEqual(alignment3.query.id, "name6")
        alignment3 = next(alignments3)
        self.assertEqual(alignment3.query.id, "name7")
        alignment1 = next(alignments1)
        self.assertEqual(alignment1.query.id, "name3")
        alignment1 = next(alignments1)
        self.assertEqual(alignment1.query.id, "name4")
        alignment1 = next(alignments1)
        self.assertEqual(alignment1.query.id, "name5")
        alignment2 = next(alignments2)
        self.assertEqual(alignment2.query.id, "name7")
        self.assertRaises(StopIteration, next, alignments2)
        alignment1 = next(alignments1)
        self.assertEqual(alignment1.query.id, "name6")
        alignment1 = next(alignments1)
        self.assertEqual(alignment1.query.id, "name7")
        self.assertRaises(StopIteration, next, alignments3)
        alignment1 = next(alignments1)
        self.assertEqual(alignment1.query.id, "name8")
        self.assertRaises(StopIteration, next, alignments1)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
