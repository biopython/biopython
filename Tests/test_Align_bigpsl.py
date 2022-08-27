# Copyright 2022 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Align.bigpsl module."""
import unittest
import warnings
from io import StringIO


from Bio.Align import Alignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, ExactPosition, CompoundLocation
from Bio import SeqIO
from Bio import BiopythonExperimentalWarning

with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonExperimentalWarning)
    from Bio.Align import bigpsl


try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align.bigpsl."
    ) from None


class TestAlign_dna_rna(unittest.TestCase):

    # The bigPsl file dna_rna.psl.bb was generated using these commands:
    # pslToBigPsl dna_rna.psl stdout | sort -k1,1 -k2,2n > dna_rna.bigPslInput
    # bedToBigBed -type=bed12+13 -tab -as=bigPsl.as dna_rna.bigPslInput hg38.chrom.sizes dna_rna.psl.bb

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
        self.dna = data
        records = SeqIO.parse("Blat/rna.fa", "fasta")
        self.rna = {record.id: record.seq for record in records}

    def test_reading(self):
        """Test parsing dna_rna.psl.bb."""
        path = "Blat/dna_rna.psl.bb"
        alignments = bigpsl.AlignmentIterator(path)
        self.assertEqual(
            str(alignments.declaration),
            """\
table bigPsl
"bigPsl pairwise alignment"
(
   string          chrom;           "Reference sequence chromosome or scaffold"
   uint            chromStart;      "Start position in chromosome"
   uint            chromEnd;        "End position in chromosome"
   string          name;            "Name or ID of item, ideally both human readable and unique"
   uint            score;           "Score (0-1000)"
   char[1]         strand;          "+ or - indicates whether the query aligns to the + or - strand on the reference"
   uint            thickStart;      "Start of where display should be thick (start codon)"
   uint            thickEnd;        "End of where display should be thick (stop codon)"
   uint            reserved;        "RGB value (use R,G,B string in input file)"
   int             blockCount;      "Number of blocks"
   int[blockCount] blockSizes;      "Comma separated list of block sizes"
   int[blockCount] chromStarts;     "Start positions relative to chromStart"
   uint            oChromStart;     "Start position in other chromosome"
   uint            oChromEnd;       "End position in other chromosome"
   char[1]         oStrand;         "+ or -, - means that psl was reversed into BED-compatible coordinates"
   uint            oChromSize;      "Size of other chromosome."
   int[blockCount] oChromStarts;    "Start positions relative to oChromStart or from oChromStart+oChromSize depending on strand"
   lstring         oSequence;       "Sequence on other chrom (or edit list, or empty)"
   string          oCDS;            "CDS in NCBI format"
   uint            chromSize;       "Size of target chromosome"
   uint            match;           "Number of bases matched."
   uint            misMatch;        " Number of bases that don't match "
   uint            repMatch;        " Number of bases that match but are part of repeats "
   uint            nCount;          " Number of 'N' bases "
   uint            seqType;         "0=empty, 1=nucleotide, 2=amino_acid"
)
""",
        )
        self.assertEqual(len(alignments.targets), 1)
        self.assertEqual(alignments.targets[0].id, "chr3")
        self.assertEqual(len(alignments.targets[0]), 198295559)
        self.assertEqual(len(alignments), 4)
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 175)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 6)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.thickStart, 42530895)
        self.assertEqual(alignment.thickEnd, 42532606)
        self.assertEqual(alignment.itemRgb, "0")
        self.assertEqual(alignment.shape, (2, 1711))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr3")
        self.assertEqual(alignment.query.id, "NR_046654.1")
        self.assertEqual(len(alignment.target.seq), 198295559)
        self.assertEqual(len(alignment.query.seq), 181)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[42530895, 42530958, 42532020, 42532095, 42532563, 42532606],
                             [     181,      118,      118,       43,       43,        0]])
                # fmt: on
            )
        )
        dna = Seq(self.dna, length=len(alignment.target))
        alignment.target.seq = dna
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
        matches = sum(
            alignment.substitutions[c, c] for c in alignment.substitutions.alphabet
        )
        repMatches = sum(
            alignment.substitutions[c, c.swapcase()]
            for c in alignment.substitutions.alphabet
        )
        self.assertEqual(matches, alignment.matches)
        self.assertEqual(repMatches, alignment.repMatches)
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 172)
        self.assertEqual(alignment.misMatches, 1)
        self.assertEqual(alignment.repMatches, 6)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.thickStart, 42530895)
        self.assertEqual(alignment.thickEnd, 42532606)
        self.assertEqual(alignment.itemRgb, "0")
        self.assertEqual(alignment.shape, (2, 1714))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr3")
        self.assertEqual(alignment.query.id, "NR_046654.1_modified")
        self.assertEqual(len(alignment.target.seq), 198295559)
        self.assertEqual(len(alignment.query.seq), 190)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[42530895, 42530922, 42530922, 42530958, 42532020,
                              42532037, 42532039, 42532095, 42532563, 42532606],
                             [     185,      158,      155,      119,      119,
                                   102,      102,       46,       46,        3],
                            ])
                # fmt: on
            )
        )
        dna = Seq(self.dna, length=len(alignment.target))
        alignment.target.seq = dna
        alignment.query.seq = self.rna[alignment.query.id]
        self.assertTrue(
            numpy.array_equal(
                alignment.substitutions,
                # fmt: off
# flake8: noqa
            numpy.array([[34.,  0.,  0.,  1.,  0.,  0.,  0.,  0.],
                         [ 0., 40.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0., 57.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  0., 41.,  0.,  0.,  0.,  0.],
                         [ 2.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  1.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  3.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                        ]),
            )
        )
        self.assertEqual(alignment.substitutions.alphabet, "ACGTacgt")
        matches = sum(
            alignment.substitutions[c, c] for c in alignment.substitutions.alphabet
        )
        repMatches = sum(
            alignment.substitutions[c, c.swapcase()]
            for c in alignment.substitutions.alphabet
            if c != "X"
        )
        self.assertEqual(matches, alignment.matches)
        self.assertEqual(repMatches, alignment.repMatches)
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 165)
        self.assertEqual(alignment.misMatches, 0)
        self.assertEqual(alignment.repMatches, 39)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.thickStart, 48663767)
        self.assertEqual(alignment.thickEnd, 48669174)
        self.assertEqual(alignment.itemRgb, "0")
        self.assertEqual(alignment.shape, (2, 5407))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr3")
        self.assertEqual(alignment.query.id, "NR_111921.1")
        self.assertEqual(len(alignment.target.seq), 198295559)
        self.assertEqual(len(alignment.query.seq), 216)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array( [[48663767, 48663813, 48665640, 48665722, 48669098, 48669174],
                              [       0,        46,      46,      128,      128,      204]]),
                # fmt: on
            )
        )
        dna = Seq(self.dna, length=len(alignment.target.seq))
        alignment.target.seq = dna
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
        matches = sum(
            alignment.substitutions[c, c] for c in alignment.substitutions.alphabet
        )
        repMatches = sum(
            alignment.substitutions[c, c.swapcase()]
            for c in alignment.substitutions.alphabet
        )
        self.assertEqual(matches, alignment.matches)
        self.assertEqual(repMatches, alignment.repMatches)
        alignment = next(alignments)
        self.assertEqual(alignment.matches, 162)
        self.assertEqual(alignment.misMatches, 2)
        self.assertEqual(alignment.repMatches, 39)
        self.assertEqual(alignment.nCount, 0)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.thickStart, 48663767)
        self.assertEqual(alignment.thickEnd, 48669174)
        self.assertEqual(alignment.itemRgb, "0")
        self.assertEqual(alignment.shape, (2, 5409))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr3")
        self.assertEqual(alignment.query.id, "NR_111921.1_modified")
        self.assertEqual(len(alignment.target.seq), 198295559)
        self.assertEqual(len(alignment.query.seq), 220)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[48663767, 48663795, 48663796, 48663813, 48665640,
                              48665716, 48665716, 48665722, 48669098, 48669174],
                             [       3,       31,       31,       48,       48,
                                   124,      126,      132,      132,      208]
                            ])
                # fmt: on
            )
        )
        dna = Seq(self.dna, length=len(alignment.target))
        alignment.target.seq = dna
        alignment.query.seq = self.rna[alignment.query.id]
        self.assertTrue(
            numpy.array_equal(
                alignment.substitutions,
                # fmt: off
# flake8: noqa
            numpy.array([[53.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0., 34.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  2., 48.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  0., 27.,  0.,  0.,  0.,  0.],
                         [ 9.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  7.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0., 16.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  0.,  7.,  0.,  0.,  0.,  0.],
                        ]),
            )
        )
        self.assertEqual(alignment.substitutions.alphabet, "ACGTacgt")
        matches = sum(
            alignment.substitutions[c, c] for c in alignment.substitutions.alphabet
        )
        repMatches = sum(
            alignment.substitutions[c, c.swapcase()]
            for c in alignment.substitutions.alphabet
            if c != "X"
        )
        self.assertEqual(matches, alignment.matches)
        self.assertEqual(repMatches, alignment.repMatches)
        self.assertRaises(StopIteration, next, alignments)


class TestAlign_dna(unittest.TestCase):

    queries = {
        record.id: record.seq for record in SeqIO.parse("Blat/fasta_34.fa", "fasta")
    }

    def test_reading_psl_34_001(self):
        """Test parsing psl_34_001.psl.bb."""

        # The bigPsl file psl_34_001.psl.bb was generated using these commands:
        # pslToBigPsl -fa=fasta_34.fa psl_34_001.psl stdout | sort -k1,1 -k2,2n > psl_34_001.bigPslInput
        # bedToBigBed -type=bed12+13 -tab -as=bigPsl.as psl_34_001.bigPslInput hg19.chrom.sizes psl_34_001.psl.bb

        path = "Blat/psl_34_001.psl.bb"
        alignments = bigpsl.AlignmentIterator(path)
        self.assertEqual(
            str(alignments.declaration),
            """\
table bigPsl
"bigPsl pairwise alignment"
(
   string          chrom;           "Reference sequence chromosome or scaffold"
   uint            chromStart;      "Start position in chromosome"
   uint            chromEnd;        "End position in chromosome"
   string          name;            "Name or ID of item, ideally both human readable and unique"
   uint            score;           "Score (0-1000)"
   char[1]         strand;          "+ or - indicates whether the query aligns to the + or - strand on the reference"
   uint            thickStart;      "Start of where display should be thick (start codon)"
   uint            thickEnd;        "End of where display should be thick (stop codon)"
   uint            reserved;        "RGB value (use R,G,B string in input file)"
   int             blockCount;      "Number of blocks"
   int[blockCount] blockSizes;      "Comma separated list of block sizes"
   int[blockCount] chromStarts;     "Start positions relative to chromStart"
   uint            oChromStart;     "Start position in other chromosome"
   uint            oChromEnd;       "End position in other chromosome"
   char[1]         oStrand;         "+ or -, - means that psl was reversed into BED-compatible coordinates"
   uint            oChromSize;      "Size of other chromosome."
   int[blockCount] oChromStarts;    "Start positions relative to oChromStart or from oChromStart+oChromSize depending on strand"
   lstring         oSequence;       "Sequence on other chrom (or edit list, or empty)"
   string          oCDS;            "CDS in NCBI format"
   uint            chromSize;       "Size of target chromosome"
   uint            match;           "Number of bases matched."
   uint            misMatch;        " Number of bases that don't match "
   uint            repMatch;        " Number of bases that match but are part of repeats "
   uint            nCount;          " Number of 'N' bases "
   uint            seqType;         "0=empty, 1=nucleotide, 2=amino_acid"
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[39368490, 39368526],
                             [      49,       13]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[61700837, 61700871],
                             [       1,       35]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[220325687, 220325721],
                             [       47,        13]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[99388555, 99388591],
                             [      49,       13]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[112178171, 112178196],
                             [       35,        10]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[52759147, 52759154, 52759160, 52759160, 52759198],
                             [       1,        8,        8,       11,       49]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[23891310, 23891349],
                             [      10,       49]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[43252217, 43252245],
                             [      21,       49]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[553742, 553781],
                             [    49,     10]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[35483340, 35483365, 35483499, 35483510],
                             [      10,       35,       35,       46]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[54017130, 54017169],
                             [      49,       10]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[53575980, 53575997],
                             [      25,        8]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[120641740, 120641776],
                             [       49,        13]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[183925984, 183925990, 183925990, 183926028],
                             [        1,         7,        11,        49]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[42144400, 42144436],
                             [      11,       47]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[48997405, 48997442],
                             [      49,       12]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[37558157, 37558167, 37558173, 37558173, 37558191],
                             [      49,       39,       39,       29,       11]]),
                # fmt: on
            )
        )
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[61646095, 61646111],
                             [      11,       27]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[95160479, 95160520],
                             [       8,       49]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[85737865, 85737906],
                             [       9,       50]]),
                # fmt: on
            )
        )
        self.assertRaises(StopIteration, next, alignments)

    def test_reading_psl_34_003(self):
        """Test parsing psl_34_003.psl.bb."""

        # The bigPsl file psl_34_003.psl.bb was generated using these commands:
        # pslToBigPsl -fa=fasta_34.fa psl_34_003.psl stdout | sort -k1,1 -k2,2n > psl_34_003.bigPslInput
        # bedToBigBed -type=bed12+13 -tab -as=bigPsl.as psl_34_003.bigPslInput hg19.chrom.sizes psl_34_003.psl.bb

        path = "Blat/psl_34_003.psl.bb"
        alignments = bigpsl.AlignmentIterator(path)
        self.assertEqual(
            str(alignments.declaration),
            """\
table bigPsl
"bigPsl pairwise alignment"
(
   string          chrom;           "Reference sequence chromosome or scaffold"
   uint            chromStart;      "Start position in chromosome"
   uint            chromEnd;        "End position in chromosome"
   string          name;            "Name or ID of item, ideally both human readable and unique"
   uint            score;           "Score (0-1000)"
   char[1]         strand;          "+ or - indicates whether the query aligns to the + or - strand on the reference"
   uint            thickStart;      "Start of where display should be thick (start codon)"
   uint            thickEnd;        "End of where display should be thick (stop codon)"
   uint            reserved;        "RGB value (use R,G,B string in input file)"
   int             blockCount;      "Number of blocks"
   int[blockCount] blockSizes;      "Comma separated list of block sizes"
   int[blockCount] chromStarts;     "Start positions relative to chromStart"
   uint            oChromStart;     "Start position in other chromosome"
   uint            oChromEnd;       "End position in other chromosome"
   char[1]         oStrand;         "+ or -, - means that psl was reversed into BED-compatible coordinates"
   uint            oChromSize;      "Size of other chromosome."
   int[blockCount] oChromStarts;    "Start positions relative to oChromStart or from oChromStart+oChromSize depending on strand"
   lstring         oSequence;       "Sequence on other chrom (or edit list, or empty)"
   string          oCDS;            "CDS in NCBI format"
   uint            chromSize;       "Size of target chromosome"
   uint            match;           "Number of bases matched."
   uint            misMatch;        " Number of bases that don't match "
   uint            repMatch;        " Number of bases that match but are part of repeats "
   uint            nCount;          " Number of 'N' bases "
   uint            seqType;         "0=empty, 1=nucleotide, 2=amino_acid"
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[53575980, 53575997],
                             [      25,        8]]),
                # fmt: on
            )
        )
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[61646095, 61646111],
                             [      11,       27]]),
                # fmt: on
            )
        )
        self.assertRaises(StopIteration, next, alignments)

    def test_reading_psl_34_004(self):
        """Test parsing psl_34_004.psl.bb."""

        # The bigPsl file psl_34_004.psl.bb was generated using these commands:
        # pslToBigPsl -fa=fasta_34.fa psl_34_004.psl stdout | sort -k1,1 -k2,2n > psl_34_004.bigPslInput
        # bedToBigBed -type=bed12+13 -tab -as=bigPsl.as psl_34_004.bigPslInput hg19.chrom.sizes psl_34_004.psl.bb

        path = "Blat/psl_34_004.psl.bb"
        alignments = bigpsl.AlignmentIterator(path)
        self.assertEqual(
            str(alignments.declaration),
            """\
table bigPsl
"bigPsl pairwise alignment"
(
   string          chrom;           "Reference sequence chromosome or scaffold"
   uint            chromStart;      "Start position in chromosome"
   uint            chromEnd;        "End position in chromosome"
   string          name;            "Name or ID of item, ideally both human readable and unique"
   uint            score;           "Score (0-1000)"
   char[1]         strand;          "+ or - indicates whether the query aligns to the + or - strand on the reference"
   uint            thickStart;      "Start of where display should be thick (start codon)"
   uint            thickEnd;        "End of where display should be thick (stop codon)"
   uint            reserved;        "RGB value (use R,G,B string in input file)"
   int             blockCount;      "Number of blocks"
   int[blockCount] blockSizes;      "Comma separated list of block sizes"
   int[blockCount] chromStarts;     "Start positions relative to chromStart"
   uint            oChromStart;     "Start position in other chromosome"
   uint            oChromEnd;       "End position in other chromosome"
   char[1]         oStrand;         "+ or -, - means that psl was reversed into BED-compatible coordinates"
   uint            oChromSize;      "Size of other chromosome."
   int[blockCount] oChromStarts;    "Start positions relative to oChromStart or from oChromStart+oChromSize depending on strand"
   lstring         oSequence;       "Sequence on other chrom (or edit list, or empty)"
   string          oCDS;            "CDS in NCBI format"
   uint            chromSize;       "Size of target chromosome"
   uint            match;           "Number of bases matched."
   uint            misMatch;        " Number of bases that don't match "
   uint            repMatch;        " Number of bases that match but are part of repeats "
   uint            nCount;          " Number of 'N' bases "
   uint            seqType;         "0=empty, 1=nucleotide, 2=amino_acid"
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[39368490, 39368526],
                             [      49,       13]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[61700837, 61700871],
                             [       1,       35]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[220325687, 220325721],
                             [       47,        13]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[99388555, 99388591],
                             [      49,       13]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[112178171, 112178196],
                             [       35,        10]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[52759147, 52759154, 52759160, 52759160, 52759198],
                             [       1,        8,        8,       11,       49]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[23891310, 23891349],
                             [      10,       49]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[43252217, 43252245],
                             [      21,       49]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[553742, 553781],
                             [    49,     10]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[35483340, 35483365, 35483499, 35483510],
                             [      10,       35,       35,       46]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[54017130, 54017169],
                             [      49,       10]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[120641740, 120641776],
                             [       49,        13]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[183925984, 183925990, 183925990, 183926028],
                             [        1,         7,        11,        49]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[42144400, 42144436],
                             [      11,       47]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[48997405, 48997442],
                             [      49,       12]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[37558157, 37558167, 37558173, 37558173, 37558191],
                             [      49,       39,       39,       29,       11]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[95160479, 95160520],
                             [       8,       49]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[85737865, 85737906],
                             [       9,       50]]),
                # fmt: on
            )
        )
        self.assertRaises(StopIteration, next, alignments)

    def test_reading_psl_34_005(self):
        """Test parsing psl_34_005.psl.bb."""

        # The bigPsl file psl_34_005.psl.bb was generated using these commands:
        # pslToBigPsl -fa=fasta_34.fa psl_34_005.psl stdout | sort -k1,1 -k2,2n > psl_34_005.bigPslInput
        # bedToBigBed -type=bed12+13 -tab -as=bigPsl.as psl_34_005.bigPslInput hg19.chrom.sizes psl_34_005.psl.bb

        path = "Blat/psl_34_005.psl.bb"
        alignments = bigpsl.AlignmentIterator(path)
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
                numpy.array([[39368490, 39368526],
                             [      49,       13]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[61700837, 61700871],
                             [       1,       35]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[220325687, 220325721],
                             [       47,        13]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[99388555, 99388591],
                             [      49,       13]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[112178171, 112178196],
                             [       35,        10]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[52759147, 52759154, 52759160, 52759160, 52759198],
                             [       1,        8,        8,       11,       49]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[23891310, 23891349],
                             [      10,       49]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[43252217, 43252245],
                             [      21,       49]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[553742, 553781],
                             [    49,     10]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[35483340, 35483365, 35483499, 35483510],
                             [      10,       35,       35,       46]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[54017130, 54017169],
                             [      49,       10]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[53575980, 53575997],
                             [      25,        8]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[120641740, 120641776],
                             [       49,        13]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[183925984, 183925990, 183925990, 183926028],
                             [        1,         7,        11,        49]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[42144400, 42144436],
                             [      11,       47]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[48997405, 48997442],
                             [      49,       12]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[37558157, 37558167, 37558173, 37558173, 37558191],
                             [      49,       39,       39,       29,       11]]),
                # fmt: on
            )
        )
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[61646095, 61646111],
                             [      11,       27]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[95160479, 95160520],
                             [       8,       49]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[85737865, 85737906],
                             [       9,       50]]),
                # fmt: on
            )
        )
        self.assertRaises(StopIteration, next, alignments)


class TestAlign_dnax_prot(unittest.TestCase):

    queries = {
        record.id: record.seq
        for record in SeqIO.parse("Blat/CAG33136.1.fasta", "fasta")
    }

    def test_reading_psl_35_001(self):
        """Test parsing psl_35_001.psl.bb."""

        # The bigPsl file psl_35_001.psl.bb was generated using these commands:
        # pslToBigPsl -fa=CAG33136.1.fasta psl_35_001.psl stdout | sort -k1,1 -k2,2n > psl_35_001.bigPslInput
        # bedToBigBed -type=bed12+13 -tab -as=bigPsl.as psl_35_001.bigPslInput hg38.chrom.sizes psl_35_001.psl.bb

        path = "Blat/psl_35_001.psl.bb"
        alignments = bigpsl.AlignmentIterator(path)
        self.assertEqual(
            str(alignments.declaration),
            """\
table bigPsl
"bigPsl pairwise alignment"
(
   string          chrom;           "Reference sequence chromosome or scaffold"
   uint            chromStart;      "Start position in chromosome"
   uint            chromEnd;        "End position in chromosome"
   string          name;            "Name or ID of item, ideally both human readable and unique"
   uint            score;           "Score (0-1000)"
   char[1]         strand;          "+ or - indicates whether the query aligns to the + or - strand on the reference"
   uint            thickStart;      "Start of where display should be thick (start codon)"
   uint            thickEnd;        "End of where display should be thick (stop codon)"
   uint            reserved;        "RGB value (use R,G,B string in input file)"
   int             blockCount;      "Number of blocks"
   int[blockCount] blockSizes;      "Comma separated list of block sizes"
   int[blockCount] chromStarts;     "Start positions relative to chromStart"
   uint            oChromStart;     "Start position in other chromosome"
   uint            oChromEnd;       "End position in other chromosome"
   char[1]         oStrand;         "+ or -, - means that psl was reversed into BED-compatible coordinates"
   uint            oChromSize;      "Size of other chromosome."
   int[blockCount] oChromStarts;    "Start positions relative to oChromStart or from oChromStart+oChromSize depending on strand"
   lstring         oSequence;       "Sequence on other chrom (or edit list, or empty)"
   string          oCDS;            "CDS in NCBI format"
   uint            chromSize;       "Size of target chromosome"
   uint            match;           "Number of bases matched."
   uint            misMatch;        " Number of bases that don't match "
   uint            repMatch;        " Number of bases that match but are part of repeats "
   uint            nCount;          " Number of 'N' bases "
   uint            seqType;         "0=empty, 1=nucleotide, 2=amino_acid"
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[75549820, 75549865, 75567225, 75567225, 75567312],
                             [       0,       15,       15,      113,      142]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[75560749, 75560881],
                             [      17,       61]]),
                # fmt: on
            )
        )
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[75566694, 75566850],
                             [      61,      113]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[75569459, 75569507],
                             [     142,      158]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[75594914, 75594989],
                             [     158,      183]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[75604767, 75604827, 75605728, 75605809],
                             [     183,      203,      203,      230]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[41257605, 41257731, 41263227, 41263227, 41263290],
                             [      17,       59,       59,      162,      183]]),
                # fmt: on
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
        self.assertEqual(alignment.query.seq, self.queries[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[41260685, 41260787],
                             [      76,      110]]),
                # fmt: on
            )
        )
        self.assertRaises(StopIteration, next, alignments)

    def test_reading_psl_35_002(self):
        """Test parsing psl_35_002.psl.bb."""

        # The bigPsl file psl_35_002.psl.bb was generated using these commands:
        # pslToBigPsl -fa=CAG33136.1.fasta psl_35_002.psl stdout | grep -v KI538594 | sort -k1,1 -k2,2n > psl_35_002.bigPslInput
        # (where we excluded KI538594 because its alignment has a negative gap)
        # bedToBigBed -type=bed12+13 -tab -as=bigPsl.as psl_35_002.bigPslInput balAcu1.chrom.sizes psl_35_002.psl.bb

        # See below for a description of the file balAcu1.fa.
        # We use this file here so we can check the SeqFeatures.
        records = SeqIO.parse("Blat/balAcu1.fa", "fasta")
        self.dna = {}
        for record in records:
            name, start_end = record.id.split(":")
            start, end = start_end.split("-")
            start = int(start)
            end = int(end)
            sequence = str(record.seq)
            self.dna[name] = Seq({start: sequence}, length=end)
        path = "Blat/psl_35_002.psl.bb"
        alignments = bigpsl.AlignmentIterator(path)
        self.assertEqual(
            str(alignments.declaration),
            """\
table bigPsl
"bigPsl pairwise alignment"
(
   string          chrom;           "Reference sequence chromosome or scaffold"
   uint            chromStart;      "Start position in chromosome"
   uint            chromEnd;        "End position in chromosome"
   string          name;            "Name or ID of item, ideally both human readable and unique"
   uint            score;           "Score (0-1000)"
   char[1]         strand;          "+ or - indicates whether the query aligns to the + or - strand on the reference"
   uint            thickStart;      "Start of where display should be thick (start codon)"
   uint            thickEnd;        "End of where display should be thick (stop codon)"
   uint            reserved;        "RGB value (use R,G,B string in input file)"
   int             blockCount;      "Number of blocks"
   int[blockCount] blockSizes;      "Comma separated list of block sizes"
   int[blockCount] chromStarts;     "Start positions relative to chromStart"
   uint            oChromStart;     "Start position in other chromosome"
   uint            oChromEnd;       "End position in other chromosome"
   char[1]         oStrand;         "+ or -, - means that psl was reversed into BED-compatible coordinates"
   uint            oChromSize;      "Size of other chromosome."
   int[blockCount] oChromStarts;    "Start positions relative to oChromStart or from oChromStart+oChromSize depending on strand"
   lstring         oSequence;       "Sequence on other chrom (or edit list, or empty)"
   string          oCDS;            "CDS in NCBI format"
   uint            chromSize;       "Size of target chromosome"
   uint            match;           "Number of bases matched."
   uint            misMatch;        " Number of bases that don't match "
   uint            repMatch;        " Number of bases that match but are part of repeats "
   uint            nCount;          " Number of 'N' bases "
   uint            seqType;         "0=empty, 1=nucleotide, 2=amino_acid"
)
""",
        )
        self.assertEqual(len(alignments.targets), 2)
        self.assertEqual(alignments.targets[0].id, "KI537194")
        self.assertEqual(len(alignments.targets[0]), 37111980)
        self.assertEqual(alignments.targets[1].id, "KI537979")
        self.assertEqual(len(alignments.targets[1]), 14052872)
        self.assertEqual(len(alignments), 2)
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
                # fmt: off
# flake8: noqa
                numpy.array([[20873021, 20872472, 20872471, 20872471, 20872390],
                             [       0,      183,      183,      203,      230]]),
                # fmt: on
            )
        )
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
                # fmt: off
# flake8: noqa
                numpy.array([[9712654, 9712786, 9715941, 9716097, 9716445, 9716532, 9718374,
                              9718422, 9739264, 9739339, 9743706, 9743766, 9744511, 9744592],
                             [     17,      61,      61,     113,     113,     142,     142,
                                  158,     158,     183,     183,     203,     203,     230]]),
                # fmt: on
            )
        )
        self.assertRaises(StopIteration, next, alignments)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
