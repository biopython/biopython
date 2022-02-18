# Copyright 2022 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Align.sam module."""
import unittest
from io import StringIO


from Bio.Align import Alignment, sam
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align.sam."
    ) from None


class TestAlign_dna_rna(unittest.TestCase):

    # The SAM file dna_rna.sam was generated using these commands:
    # twoBitToFa hg38.2bit stdout | samtools dict -a hg38 -s "Homo sapiens" | grep -v chrUn | grep -v alt | grep -v random > dna_rna.sam
    # psl2sam.pl dna_rna.psl >> dna_rna.sam
    # The CIGAR string was then edited to replace D by N for introns and H by S
    # where appropriate.

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
        self.rna["NR_111921.1"] = self.rna["NR_111921.1"][:-12]
        self.rna["NR_111921.1_modified"] = self.rna["NR_111921.1_modified"][:-12]
        # Last 12 nucleotides were clipped by Blat as the poly(A) tail

    def test_reading(self):
        """Test parsing dna_rna.sam."""
        path = "Blat/dna_rna.sam"
        alignments = sam.AlignmentIterator(path)
        self.assertEqual(list(alignments.metadata), ["HD", "SQ"])
        self.assertEqual(alignments.metadata["HD"], {"VN": "1.0", "SO": "unsorted"})
        self.assertEqual(len(alignments.metadata["SQ"]), 25)
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 5407))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr3")
        self.assertEqual(alignment.query.id, "NR_111921.1")
        self.assertEqual(len(alignment.target.seq), 198295559)
        self.assertEqual(len(alignment.query.seq), len(self.rna[alignment.query.id]))
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
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 1711))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr3")
        self.assertEqual(alignment.query.id, "NR_046654.1")
        self.assertEqual(len(alignment.target.seq), 198295559)
        self.assertEqual(len(alignment.query.seq), len(self.rna[alignment.query.id]))
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
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 5412))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr3")
        self.assertEqual(alignment.query.id, "NR_111921.1_modified")
        self.assertEqual(len(alignment.target.seq), 198295559)
        self.assertEqual(len(alignment.query.seq), len(self.rna[alignment.query.id]))
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[48663767, 48663767, 48663795, 48663796, 48663813,
                              48665640, 48665716, 48665716, 48665722, 48669098,
                              48669174],
                             [       0,        3,       31,       31,       48,
                                    48,      124,      126,      132,      132,
                                   208],
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
                         [ 0.,  0., 50.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  0., 27.,  0.,  0.,  0.,  0.],
                         [ 9.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  7.,  0.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0., 16.,  0.,  0.,  0.,  0.,  0.],
                         [ 0.,  0.,  0.,  7.,  0.,  0.,  0.,  0.],
                        ]),
            )
        )
        self.assertEqual(alignment.substitutions.alphabet, "ACGTacgt")
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 1722))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr3")
        self.assertEqual(alignment.query.id, "NR_046654.1_modified")
        self.assertEqual(len(alignment.target.seq), 198295559)
        self.assertEqual(len(alignment.query.seq), len(self.rna[alignment.query.id]))
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[42530895, 42530895, 42530922, 42530922, 42530958,
                              42532020, 42532037, 42532039, 42532095, 42532563,
                              42532606, 42532606],
                             [     190,      185,      158,      155,      119,
                                   119,      102,      102,       46,       46,
                                     3,        0],
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
            numpy.array([[35.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
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
        self.assertRaises(StopIteration, next, alignments)

    def test_writing(self):
        """Test writing the alignments in dna_rna.sam."""
        path = "Blat/dna_rna.sam"
        with open(path) as stream:
            original_data = stream.read()
        alignments = sam.AlignmentIterator(path)
        stream = StringIO()
        writer = sam.AlignmentWriter(stream)
        n = writer.write_file(alignments, mincount=4, maxcount=4)
        self.assertEqual(n, 4)
        stream.seek(0)
        written_data = stream.read()
        stream.close()
        self.assertEqual(original_data, written_data)
        # Try this again. This time, we first strip the matches, misMatches,
        # repMatches, and nCount attributes from each alignment, and insert the
        # appropriate sequence data in each alignment. The writer will then
        # recalculate the matches, misMatches, repMatches, and nCount values
        # from the sequence data and the alignment, and store those values in
        # the PSL file.
        alignments = []
        for alignment in sam.AlignmentIterator(path):
            del alignment.matches
            del alignment.misMatches
            del alignment.repMatches
            del alignment.nCount
            dna = Seq(self.dna, length=len(alignment.target))
            alignment.target.seq = dna
            alignment.query.seq = self.rna[alignment.sequences[1].id]
            alignments.append(alignment)
        stream = StringIO()
        writer = sam.AlignmentWriter(stream, mask="lower")
        n = writer.write_file(alignments, mincount=4, maxcount=4)
        self.assertEqual(n, 4)
        stream.seek(0)
        written_data = stream.read()
        stream.close()
        self.assertEqual(original_data, written_data)


class TestAlign_dna(unittest.TestCase):

    # The SAM files were generated using these commands:
    # twoBitInfo hg18.2bit stdout | grep -v _random | grep -v _hap | sort -n -k 2 -r > hg18.chrom.sizes
    # twoBitInfo hg19.2bit stdout | grep -v chrUn | grep -v _random | grep -v _hap |  sort -n -k 2 -r > hg19.chrom.sizes

    # psl2sam.pl psl_34_001.psl | samtools view -h -t hg18.chrom.sizes - > psl_34_001.sam
    # psl2sam.pl psl_34_003.psl | samtools view -h -t hg18.chrom.sizes - > psl_34_003.sam
    # psl2sam.pl psl_34_004.psl | samtools view -h -t hg19.chrom.sizes - > psl_34_004.sam
    # psl2sam.pl psl_34_005.psl | samtools view -h -t hg19.chrom.sizes - > psl_34_005.sam

    # Note that psl_34_002 was not included as the SAM format no longer allows
    # an empty SAM file.

    def test_reading_psl_34_001(self):
        """Test parsing psl_34_001.sam."""
        path = "Blat/psl_34_001.sam"
        alignments = sam.AlignmentIterator(path)
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
        self.assertEqual(len(alignment.query.seq), len(self.rna[alignment.query.id]))
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
                # fmt: off
# flake8: noqa
                numpy.array([[85737865, 85737906],
                             [       9,       50]]),
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
                # fmt: off
# flake8: noqa
                numpy.array([[42144400, 42144436],
                             [      11,       47]]),
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
                # fmt: off
# flake8: noqa
                numpy.array([[52759147, 52759154, 52759160, 52759160, 52759198],
                             [       1,        8,        8,       11,       49]]),
                # fmt: on
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
                # fmt: off
# flake8: noqa
                numpy.array([[1207056, 1207106],
                             [      0,      50]]),
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
                # fmt: off
# flake8: noqa
                numpy.array([[37558157, 37558167, 37558173, 37558173, 37558191],
                             [      49,       39,       39,       29,       11]]),
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
                # fmt: off
# flake8: noqa
                numpy.array([[120641740, 120641776],
                             [       49,        13]]),
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
                # fmt: off
# flake8: noqa
                numpy.array([[553742, 553781],
                             [    49,     10]]),
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
                # fmt: off
# flake8: noqa
                numpy.array([[39368490, 39368526],
                             [      49,       13]]),
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
        self.assertRaises(StopIteration, next, alignments)

    def test_writing_psl_34_001(self):
        """Test writing the alignments in psl_34_001.sam."""
        path = "Blat/psl_34_001.sam"
        with open(path) as stream:
            original_data = stream.read()
        alignments = sam.AlignmentIterator(path)
        stream = StringIO()
        writer = sam.AlignmentWriter(stream)
        n = writer.write_file(alignments, mincount=22, maxcount=22)
        self.assertEqual(n, 22)
        stream.seek(0)
        written_data = stream.read()
        stream.close()
        self.assertEqual(original_data, written_data)

    def test_reading_psl_34_002(self):
        """Test parsing psl_34_002.sam."""
        path = "Blat/psl_34_002.sam"
        alignments = sam.AlignmentIterator(path)
        self.assertEqual(alignments.metadata["version"], "3")
        self.assertRaises(StopIteration, next, alignments)

    def test_writing_psl_34_002(self):
        """Test writing the alignments in psl_34_002.sam."""
        path = "Blat/psl_34_002.sam"
        with open(path) as stream:
            original_data = stream.read()
        alignments = sam.AlignmentIterator(path)
        stream = StringIO()
        writer = sam.AlignmentWriter(stream)
        n = writer.write_file(alignments, mincount=0, maxcount=0)
        self.assertEqual(n, 0)
        stream.seek(0)
        written_data = stream.read()
        stream.close()
        self.assertEqual(original_data, written_data)

    def test_reading_psl_34_003(self):
        """Test parsing psl_34_003.sam."""
        path = "Blat/psl_34_003.sam"
        alignments = sam.AlignmentIterator(path)
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
                # fmt: off
# flake8: noqa
                numpy.array([[61646095, 61646111],
                             [      11,       27]]),
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
        self.assertRaises(StopIteration, next, alignments)

    def test_writing_psl_34_003(self):
        """Test writing the alignments in psl_34_003.sam."""
        path = "Blat/psl_34_003.sam"
        with open(path) as stream:
            original_data = stream.read()
        alignments = sam.AlignmentIterator(path)
        stream = StringIO()
        writer = sam.AlignmentWriter(stream)
        n = writer.write_file(alignments, mincount=3, maxcount=3)
        self.assertEqual(n, 3)
        stream.seek(0)
        written_data = stream.read()
        stream.close()
        self.assertEqual(original_data, written_data)

    def test_reading_psl_34_004(self):
        """Test parsing psl_34_004.sam."""
        path = "Blat/psl_34_004.sam"
        alignments = sam.AlignmentIterator(path)
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
                # fmt: off
# flake8: noqa
                numpy.array([[85737865, 85737906],
                             [       9,       50]]),
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
                # fmt: off
# flake8: noqa
                numpy.array([[42144400, 42144436],
                             [      11,       47]]),
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
                # fmt: off
# flake8: noqa
                numpy.array([[52759147, 52759154, 52759160, 52759160, 52759198],
                             [       1,        8,        8,       11,       49]]),
                # fmt: on
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
                # fmt: off
# flake8: noqa
                numpy.array([[1207056, 1207106],
                             [      0,      50]]),
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
                # fmt: off
# flake8: noqa
                numpy.array([[37558157, 37558167, 37558173, 37558173, 37558191],
                             [      49,       39,       39,       29,       11]]),
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
                # fmt: off
# flake8: noqa
                numpy.array([[120641740, 120641776],
                             [       49,        13]]),
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
                # fmt: off
# flake8: noqa
                numpy.array([[553742, 553781],
                             [    49,     10]]),
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
                # fmt: off
# flake8: noqa
                numpy.array([[39368490, 39368526],
                             [      49,       13]]),
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
        self.assertRaises(StopIteration, next, alignments)

    def test_writing_psl_34_004(self):
        """Test writing the alignments in psl_34_004.sam."""
        path = "Blat/psl_34_004.sam"
        with open(path) as stream:
            original_data = stream.read()
        alignments = sam.AlignmentIterator(path)
        stream = StringIO()
        writer = sam.AlignmentWriter(stream)
        n = writer.write_file(alignments, mincount=19, maxcount=19)
        self.assertEqual(n, 19)
        stream.seek(0)
        written_data = stream.read()
        stream.close()
        self.assertEqual(original_data, written_data)

    def test_reading_psl_34_005(self):
        """Test parsing psl_34_005.sam."""
        path = "Blat/psl_34_005.sam"
        alignments = sam.AlignmentIterator(path)
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
                # fmt: off
# flake8: noqa
                numpy.array([[61646095, 61646111],
                             [      11,       27]]),
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
                # fmt: off
# flake8: noqa
                numpy.array([[85737865, 85737906],
                             [       9,       50]]),
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
                # fmt: off
# flake8: noqa
                numpy.array([[42144400, 42144436],
                             [      11,       47]]),
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
                # fmt: off
# flake8: noqa
                numpy.array([[52759147, 52759154, 52759160, 52759160, 52759198],
                             [       1,        8,        8,       11,       49]]),
                # fmt: on
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
                # fmt: off
# flake8: noqa
                numpy.array([[1207056, 1207106],
                             [      0,      50]]),
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
                # fmt: off
# flake8: noqa
                numpy.array([[37558157, 37558167, 37558173, 37558173, 37558191],
                             [      49,       39,       39,       29,       11]]),
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
                # fmt: off
# flake8: noqa
                numpy.array([[120641740, 120641776],
                             [       49,        13]]),
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
                # fmt: off
# flake8: noqa
                numpy.array([[553742, 553781],
                             [    49,     10]]),
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
                # fmt: off
                numpy.array([[39368490, 39368526],
                             [      49,       13]]),
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
        self.assertRaises(StopIteration, next, alignments)

    def test_writing_psl_34_005(self):
        """Test writing the alignments in psl_34_005.sam."""
        path = "Blat/psl_34_005.sam"
        with open(path) as stream:
            original_data = stream.read()
        alignments = sam.AlignmentIterator(path)
        stream = StringIO()
        writer = sam.AlignmentWriter(stream, header=False)
        n = writer.write_file(alignments, mincount=22, maxcount=22)
        self.assertEqual(n, 22)
        stream.seek(0)
        written_data = stream.read()
        stream.close()
        self.assertEqual(original_data, written_data)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
