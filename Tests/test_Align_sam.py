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
        self.assertEqual(alignments.metadata["SQ"][0], {"SN": "chr1", "LN": 248956422, "M5": "2648ae1bacce4ec4b6cf337dcae37816", "AS": "hg38", "SP": "Homo sapiens"})
        self.assertEqual(alignments.metadata["SQ"][1], {"SN": "chr10", "LN": 133797422, "M5": "907112d17fcb73bcab1ed1c72b97ce68", "AS": "hg38", "SP": "Homo sapiens"})
        self.assertEqual(alignments.metadata["SQ"][2], {"SN": "chr11", "LN": 135086622, "M5": "1511375dc2dd1b633af8cf439ae90cec", "AS": "hg38", "SP": "Homo sapiens"})
        self.assertEqual(alignments.metadata["SQ"][3], {"SN": "chr12", "LN": 133275309, "M5": "e81e16d3f44337034695a29b97708fce", "AS": "hg38", "SP": "Homo sapiens"})
        self.assertEqual(alignments.metadata["SQ"][4], {"SN": "chr13", "LN": 114364328, "M5": "17dab79b963ccd8e7377cef59a54fe1c", "AS": "hg38", "SP": "Homo sapiens"})
        self.assertEqual(alignments.metadata["SQ"][5], {"SN": "chr14", "LN": 107043718, "M5": "acbd9552c059d9b403e75ed26c1ce5bc", "AS": "hg38", "SP": "Homo sapiens"})
        self.assertEqual(alignments.metadata["SQ"][6], {"SN": "chr15", "LN": 101991189, "M5": "f036bd11158407596ca6bf3581454706", "AS": "hg38", "SP": "Homo sapiens"})
        self.assertEqual(alignments.metadata["SQ"][7], {"SN": "chr16", "LN": 90338345, "M5": "24e7cabfba3548a2bb4dff582b9ee870", "AS": "hg38", "SP": "Homo sapiens"})
        self.assertEqual(alignments.metadata["SQ"][8], {"SN": "chr17", "LN": 83257441, "M5": "a8499ca51d6fb77332c2d242923994eb", "AS": "hg38", "SP": "Homo sapiens"})
        self.assertEqual(alignments.metadata["SQ"][9], {"SN": "chr18", "LN": 80373285, "M5": "11eeaa801f6b0e2e36a1138616b8ee9a", "AS": "hg38", "SP": "Homo sapiens"})
        self.assertEqual(alignments.metadata["SQ"][10], {"SN": "chr19", "LN": 58617616, "M5": "b0eba2c7bb5c953d1e06a508b5e487de", "AS": "hg38", "SP": "Homo sapiens"})
        self.assertEqual(alignments.metadata["SQ"][11], {"SN": "chr2", "LN": 242193529, "M5": "4bb4f82880a14111eb7327169ffb729b", "AS": "hg38", "SP": "Homo sapiens"})
        self.assertEqual(alignments.metadata["SQ"][12], {"SN": "chr20", "LN": 64444167, "M5": "b18e6c531b0bd70e949a7fc20859cb01", "AS": "hg38", "SP": "Homo sapiens"})
        self.assertEqual(alignments.metadata["SQ"][13], {"SN": "chr21", "LN": 46709983, "M5": "2f45a3455007b7e271509161e52954a9", "AS": "hg38", "SP": "Homo sapiens"})
        self.assertEqual(alignments.metadata["SQ"][14], {"SN": "chr22", "LN": 50818468, "M5": "221733a2a15e2de66d33e73d126c5109", "AS": "hg38", "SP": "Homo sapiens"})
        self.assertEqual(alignments.metadata["SQ"][15], {"SN": "chr3", "LN": 198295559, "M5": "a48af509898d3736ba95dc0912c0b461", "AS": "hg38", "SP": "Homo sapiens"})
        self.assertEqual(alignments.metadata["SQ"][16], {"SN": "chr4", "LN": 190214555, "M5": "3210fecf1eb92d5489da4346b3fddc6e", "AS": "hg38", "SP": "Homo sapiens"})
        self.assertEqual(alignments.metadata["SQ"][17], {"SN": "chr5", "LN": 181538259, "M5": "f7f05fb7ceea78cbc32ce652c540ff2d", "AS": "hg38", "SP": "Homo sapiens"})
        self.assertEqual(alignments.metadata["SQ"][18], {"SN": "chr6", "LN": 170805979, "M5": "6a48dfa97e854e3c6f186c8ff973f7dd", "AS": "hg38", "SP": "Homo sapiens"})
        self.assertEqual(alignments.metadata["SQ"][19], {"SN": "chr7", "LN": 159345973, "M5": "94eef2b96fd5a7c8db162c8c74378039", "AS": "hg38", "SP": "Homo sapiens"})
        self.assertEqual(alignments.metadata["SQ"][20], {"SN": "chr8", "LN": 145138636, "M5": "c67955b5f7815a9a1edfaa15893d3616", "AS": "hg38", "SP": "Homo sapiens"})
        self.assertEqual(alignments.metadata["SQ"][21], {"SN": "chr9", "LN": 138394717, "M5": "addd2795560986b7491c40b1faa3978a", "AS": "hg38", "SP": "Homo sapiens"})
        self.assertEqual(alignments.metadata["SQ"][22], {"SN": "chrM", "LN": 16569, "M5": "c68f52674c9fb33aef52dcf399755519", "AS": "hg38", "SP": "Homo sapiens"})
        self.assertEqual(alignments.metadata["SQ"][23], {"SN": "chrX", "LN": 156040895, "M5": "49527016a48497d9d1cbd8e4a9049bd3", "AS": "hg38", "SP": "Homo sapiens"})
        self.assertEqual(alignments.metadata["SQ"][24], {"SN": "chrY", "LN": 57227415, "M5": "b2b7e6369564d89059e763cd6e736837", "AS": "hg38", "SP": "Homo sapiens"})
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
    # twoBitInfo hg19.2bit stdout | grep -v chrUn | grep -v _random | grep -v _hap |  sort -n -k 2 -r > hg19.chrom.sizes

    # psl2sam.pl psl_34_001.psl | samtools view -h -t hg19.chrom.sizes - > psl_34_001.sam
    # psl2sam.pl psl_34_003.psl | samtools view -h -t hg19.chrom.sizes - > psl_34_003.sam
    # psl2sam.pl psl_34_004.psl | samtools view -h -t hg19.chrom.sizes - > psl_34_004.sam
    # psl2sam.pl psl_34_005.psl | samtools view -h -t hg19.chrom.sizes - > psl_34_005.sam

    # Note that psl_34_002 was not included as the SAM format no longer allows
    # an empty SAM file.

    # The hard clipping symbols H were replaced by soft clipping symbols S in
    # the file psl_34_005.sam.

    def test_reading_psl_34_001(self):
        """Test parsing psl_34_001.sam."""
        path = "Blat/psl_34_001.sam"
        alignments = sam.AlignmentIterator(path)
        self.assertEqual(list(alignments.metadata), ["SQ", "PG"])
        self.assertEqual(len(alignments.metadata["SQ"]), 25)
        self.assertEqual(alignments.metadata["SQ"][0], {"SN": "chr1", "LN": 249250621})
        self.assertEqual(alignments.metadata["SQ"][1], {"SN": "chr2", "LN": 243199373})
        self.assertEqual(alignments.metadata["SQ"][2], {"SN": "chr3", "LN": 198022430})
        self.assertEqual(alignments.metadata["SQ"][3], {"SN": "chr4", "LN": 191154276})
        self.assertEqual(alignments.metadata["SQ"][4], {"SN": "chr5", "LN": 180915260})
        self.assertEqual(alignments.metadata["SQ"][5], {"SN": "chr6", "LN": 171115067})
        self.assertEqual(alignments.metadata["SQ"][6], {"SN": "chr7", "LN": 159138663})
        self.assertEqual(alignments.metadata["SQ"][7], {"SN": "chrX", "LN": 155270560})
        self.assertEqual(alignments.metadata["SQ"][8], {"SN": "chr8", "LN": 146364022})
        self.assertEqual(alignments.metadata["SQ"][9], {"SN": "chr9", "LN": 141213431})
        self.assertEqual(alignments.metadata["SQ"][10], {"SN": "chr10", "LN": 135534747})
        self.assertEqual(alignments.metadata["SQ"][11], {"SN": "chr11", "LN": 135006516})
        self.assertEqual(alignments.metadata["SQ"][12], {"SN": "chr12", "LN": 133851895})
        self.assertEqual(alignments.metadata["SQ"][13], {"SN": "chr13", "LN": 115169878})
        self.assertEqual(alignments.metadata["SQ"][14], {"SN": "chr14", "LN": 107349540})
        self.assertEqual(alignments.metadata["SQ"][15], {"SN": "chr15", "LN": 102531392})
        self.assertEqual(alignments.metadata["SQ"][16], {"SN": "chr16", "LN": 90354753})
        self.assertEqual(alignments.metadata["SQ"][17], {"SN": "chr17", "LN": 81195210})
        self.assertEqual(alignments.metadata["SQ"][18], {"SN": "chr18", "LN": 78077248})
        self.assertEqual(alignments.metadata["SQ"][19], {"SN": "chr20", "LN": 63025520})
        self.assertEqual(alignments.metadata["SQ"][20], {"SN": "chrY", "LN": 59373566})
        self.assertEqual(alignments.metadata["SQ"][21], {"SN": "chr19", "LN": 59128983})
        self.assertEqual(alignments.metadata["SQ"][22], {"SN": "chr22", "LN": 51304566})
        self.assertEqual(alignments.metadata["SQ"][23], {"SN": "chr21", "LN": 48129895})
        self.assertEqual(alignments.metadata["SQ"][24], {"SN": "chrM", "LN": 16571})
        self.assertEqual(len(alignments.metadata["PG"]), 1)
        self.assertEqual(alignments.metadata["PG"][0], {"ID": "samtools", "PN": "samtools", "VN": "1.14", "CL": "samtools view -h -t hg19.chrom.sizes -"})
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 16))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 191154276)
        self.assertEqual(len(alignment.query.seq), 16)
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
        self.assertEqual(alignment.shape, (2, 17))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 17)
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
        self.assertEqual(alignment.shape, (2, 41))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr9")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 141213431)
        self.assertEqual(len(alignment.query.seq), 41)
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
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr8")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 146364022)
        self.assertEqual(len(alignment.query.seq), 41)
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
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 51304566)
        self.assertEqual(len(alignment.query.seq), 36)
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
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 48)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[183925984, 183925990, 183925990, 183926028],
                             [        0,         6,        10,        48]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 170))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 36)
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
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 78077248)
        self.assertEqual(len(alignment.query.seq), 39)
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
        self.assertEqual(alignment.shape, (2, 28))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 78077248)
        self.assertEqual(len(alignment.query.seq), 28)
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
        self.assertEqual(alignment.shape, (2, 51))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 115169878)
        self.assertEqual(len(alignment.query.seq), 48)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[52759147, 52759157, 52759160, 52759198],
                             [       0,       10,       10,       48]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
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
        self.assertEqual(alignment.shape, (2, 34))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 34)
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
        self.assertEqual(alignment.shape, (2, 38))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 191154276)
        self.assertEqual(len(alignment.query.seq), 38)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[37558157, 37558173, 37558173, 37558191],
                             [      38,       22,       18,        0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 51304566)
        self.assertEqual(len(alignment.query.seq), 37)
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
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 36)
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
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 39)
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
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 39)
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
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 135534747)
        self.assertEqual(len(alignment.query.seq), 36)
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
        self.assertEqual(alignment.shape, (2, 25))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 135534747)
        self.assertEqual(len(alignment.query.seq), 25)
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
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 36)
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
        self.assertEqual(alignment.shape, (2, 34))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 34)
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

    def test_reading_psl_34_003(self):
        """Test parsing psl_34_003.sam."""
        path = "Blat/psl_34_003.sam"
        alignments = sam.AlignmentIterator(path)
        self.assertEqual(list(alignments.metadata), ["SQ", "PG"])
        self.assertEqual(len(alignments.metadata["SQ"]), 25)
        self.assertEqual(alignments.metadata["SQ"][0], {"SN": "chr1", "LN": 249250621})
        self.assertEqual(alignments.metadata["SQ"][1], {"SN": "chr2", "LN": 243199373})
        self.assertEqual(alignments.metadata["SQ"][2], {"SN": "chr3", "LN": 198022430})
        self.assertEqual(alignments.metadata["SQ"][3], {"SN": "chr4", "LN": 191154276})
        self.assertEqual(alignments.metadata["SQ"][4], {"SN": "chr5", "LN": 180915260})
        self.assertEqual(alignments.metadata["SQ"][5], {"SN": "chr6", "LN": 171115067})
        self.assertEqual(alignments.metadata["SQ"][6], {"SN": "chr7", "LN": 159138663})
        self.assertEqual(alignments.metadata["SQ"][7], {"SN": "chrX", "LN": 155270560})
        self.assertEqual(alignments.metadata["SQ"][8], {"SN": "chr8", "LN": 146364022})
        self.assertEqual(alignments.metadata["SQ"][9], {"SN": "chr9", "LN": 141213431})
        self.assertEqual(alignments.metadata["SQ"][10], {"SN": "chr10", "LN": 135534747})
        self.assertEqual(alignments.metadata["SQ"][11], {"SN": "chr11", "LN": 135006516})
        self.assertEqual(alignments.metadata["SQ"][12], {"SN": "chr12", "LN": 133851895})
        self.assertEqual(alignments.metadata["SQ"][13], {"SN": "chr13", "LN": 115169878})
        self.assertEqual(alignments.metadata["SQ"][14], {"SN": "chr14", "LN": 107349540})
        self.assertEqual(alignments.metadata["SQ"][15], {"SN": "chr15", "LN": 102531392})
        self.assertEqual(alignments.metadata["SQ"][16], {"SN": "chr16", "LN": 90354753})
        self.assertEqual(alignments.metadata["SQ"][17], {"SN": "chr17", "LN": 81195210})
        self.assertEqual(alignments.metadata["SQ"][18], {"SN": "chr18", "LN": 78077248})
        self.assertEqual(alignments.metadata["SQ"][19], {"SN": "chr20", "LN": 63025520})
        self.assertEqual(alignments.metadata["SQ"][20], {"SN": "chrY", "LN": 59373566})
        self.assertEqual(alignments.metadata["SQ"][21], {"SN": "chr19", "LN": 59128983})
        self.assertEqual(alignments.metadata["SQ"][22], {"SN": "chr22", "LN": 51304566})
        self.assertEqual(alignments.metadata["SQ"][23], {"SN": "chr21", "LN": 48129895})
        self.assertEqual(alignments.metadata["SQ"][24], {"SN": "chrM", "LN": 16571})
        self.assertEqual(len(alignments.metadata["PG"]), 1)
        self.assertEqual(alignments.metadata["PG"][0], {"ID": "samtools", "PN": "samtools", "VN": "1.14", "CL": "samtools view -h -t hg19.chrom.sizes -"})
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 16))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 191154276)
        self.assertEqual(len(alignment.query.seq), 16)
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
        self.assertEqual(alignment.shape, (2, 17))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 17)
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
        self.assertEqual(list(alignments.metadata), ["SQ", "PG"])
        self.assertEqual(len(alignments.metadata["SQ"]), 25)
        self.assertEqual(alignments.metadata["SQ"][0], {"SN": "chr1", "LN": 249250621})
        self.assertEqual(alignments.metadata["SQ"][1], {"SN": "chr2", "LN": 243199373})
        self.assertEqual(alignments.metadata["SQ"][2], {"SN": "chr3", "LN": 198022430})
        self.assertEqual(alignments.metadata["SQ"][3], {"SN": "chr4", "LN": 191154276})
        self.assertEqual(alignments.metadata["SQ"][4], {"SN": "chr5", "LN": 180915260})
        self.assertEqual(alignments.metadata["SQ"][5], {"SN": "chr6", "LN": 171115067})
        self.assertEqual(alignments.metadata["SQ"][6], {"SN": "chr7", "LN": 159138663})
        self.assertEqual(alignments.metadata["SQ"][7], {"SN": "chrX", "LN": 155270560})
        self.assertEqual(alignments.metadata["SQ"][8], {"SN": "chr8", "LN": 146364022})
        self.assertEqual(alignments.metadata["SQ"][9], {"SN": "chr9", "LN": 141213431})
        self.assertEqual(alignments.metadata["SQ"][10], {"SN": "chr10", "LN": 135534747})
        self.assertEqual(alignments.metadata["SQ"][11], {"SN": "chr11", "LN": 135006516})
        self.assertEqual(alignments.metadata["SQ"][12], {"SN": "chr12", "LN": 133851895})
        self.assertEqual(alignments.metadata["SQ"][13], {"SN": "chr13", "LN": 115169878})
        self.assertEqual(alignments.metadata["SQ"][14], {"SN": "chr14", "LN": 107349540})
        self.assertEqual(alignments.metadata["SQ"][15], {"SN": "chr15", "LN": 102531392})
        self.assertEqual(alignments.metadata["SQ"][16], {"SN": "chr16", "LN": 90354753})
        self.assertEqual(alignments.metadata["SQ"][17], {"SN": "chr17", "LN": 81195210})
        self.assertEqual(alignments.metadata["SQ"][18], {"SN": "chr18", "LN": 78077248})
        self.assertEqual(alignments.metadata["SQ"][19], {"SN": "chr20", "LN": 63025520})
        self.assertEqual(alignments.metadata["SQ"][20], {"SN": "chrY", "LN": 59373566})
        self.assertEqual(alignments.metadata["SQ"][21], {"SN": "chr19", "LN": 59128983})
        self.assertEqual(alignments.metadata["SQ"][22], {"SN": "chr22", "LN": 51304566})
        self.assertEqual(alignments.metadata["SQ"][23], {"SN": "chr21", "LN": 48129895})
        self.assertEqual(alignments.metadata["SQ"][24], {"SN": "chrM", "LN": 16571})
        self.assertEqual(len(alignments.metadata["PG"]), 1)
        self.assertEqual(alignments.metadata["PG"][0], {"ID": "samtools", "PN": "samtools", "VN": "1.14", "CL": "samtools view -h -t hg19.chrom.sizes -"})
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr9")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 141213431)
        self.assertEqual(len(alignment.query.seq), 41)
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
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr8")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 146364022)
        self.assertEqual(len(alignment.query.seq), 41)
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
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 51304566)
        self.assertEqual(len(alignment.query.seq), 36)
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
        self.assertEqual(alignment.shape, (2, 48))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 48)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[183925984, 183925990, 183925990, 183926028],
                             [        0,         6,        10,        48]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 170))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 36)
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
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 78077248)
        self.assertEqual(len(alignment.query.seq), 39)
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
        self.assertEqual(alignment.shape, (2, 28))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 78077248)
        self.assertEqual(len(alignment.query.seq), 28)
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
        self.assertEqual(alignment.shape, (2, 51))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 115169878)
        self.assertEqual(len(alignment.query.seq), 48)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[52759147, 52759157, 52759160, 52759198],
                             [       0,       10,       10,       48]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
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
        self.assertEqual(alignment.shape, (2, 34))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 34)
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
        self.assertEqual(alignment.shape, (2, 38))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 191154276)
        self.assertEqual(len(alignment.query.seq), 38)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[37558157, 37558173, 37558173, 37558191],
                             [      38,       22,       18,        0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 37))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 51304566)
        self.assertEqual(len(alignment.query.seq), 37)
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
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 243199373)
        self.assertEqual(len(alignment.query.seq), 36)
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
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 39)
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
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 59128983)
        self.assertEqual(len(alignment.query.seq), 39)
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
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 135534747)
        self.assertEqual(len(alignment.query.seq), 36)
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
        self.assertEqual(alignment.shape, (2, 25))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 135534747)
        self.assertEqual(len(alignment.query.seq), 25)
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
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 36)
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
        self.assertEqual(alignment.shape, (2, 34))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertEqual(len(alignment.target.seq), 249250621)
        self.assertEqual(len(alignment.query.seq), 34)
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
        self.assertEqual(list(alignments.metadata), ["SQ", "PG"])
        self.assertEqual(len(alignments.metadata["SQ"]), 25)
        self.assertEqual(alignments.metadata["SQ"][0], {"SN": "chr1", "LN": 249250621})
        self.assertEqual(alignments.metadata["SQ"][1], {"SN": "chr2", "LN": 243199373})
        self.assertEqual(alignments.metadata["SQ"][2], {"SN": "chr3", "LN": 198022430})
        self.assertEqual(alignments.metadata["SQ"][3], {"SN": "chr4", "LN": 191154276})
        self.assertEqual(alignments.metadata["SQ"][4], {"SN": "chr5", "LN": 180915260})
        self.assertEqual(alignments.metadata["SQ"][5], {"SN": "chr6", "LN": 171115067})
        self.assertEqual(alignments.metadata["SQ"][6], {"SN": "chr7", "LN": 159138663})
        self.assertEqual(alignments.metadata["SQ"][7], {"SN": "chrX", "LN": 155270560})
        self.assertEqual(alignments.metadata["SQ"][8], {"SN": "chr8", "LN": 146364022})
        self.assertEqual(alignments.metadata["SQ"][9], {"SN": "chr9", "LN": 141213431})
        self.assertEqual(alignments.metadata["SQ"][10], {"SN": "chr10", "LN": 135534747})
        self.assertEqual(alignments.metadata["SQ"][11], {"SN": "chr11", "LN": 135006516})
        self.assertEqual(alignments.metadata["SQ"][12], {"SN": "chr12", "LN": 133851895})
        self.assertEqual(alignments.metadata["SQ"][13], {"SN": "chr13", "LN": 115169878})
        self.assertEqual(alignments.metadata["SQ"][14], {"SN": "chr14", "LN": 107349540})
        self.assertEqual(alignments.metadata["SQ"][15], {"SN": "chr15", "LN": 102531392})
        self.assertEqual(alignments.metadata["SQ"][16], {"SN": "chr16", "LN": 90354753})
        self.assertEqual(alignments.metadata["SQ"][17], {"SN": "chr17", "LN": 81195210})
        self.assertEqual(alignments.metadata["SQ"][18], {"SN": "chr18", "LN": 78077248})
        self.assertEqual(alignments.metadata["SQ"][19], {"SN": "chr20", "LN": 63025520})
        self.assertEqual(alignments.metadata["SQ"][20], {"SN": "chrY", "LN": 59373566})
        self.assertEqual(alignments.metadata["SQ"][21], {"SN": "chr19", "LN": 59128983})
        self.assertEqual(alignments.metadata["SQ"][22], {"SN": "chr22", "LN": 51304566})
        self.assertEqual(alignments.metadata["SQ"][23], {"SN": "chr21", "LN": 48129895})
        self.assertEqual(alignments.metadata["SQ"][24], {"SN": "chrM", "LN": 16571})
        self.assertEqual(len(alignments.metadata["PG"]), 1)
        self.assertEqual(alignments.metadata["PG"][0], {"ID": "samtools", "PN": "samtools", "VN": "1.14", "CL": "samtools view -h -t hg19.chrom.sizes -"})
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 33))
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
                numpy.array([[61646095, 61646095, 61646111, 61646111],
                             [       0,       11,       27,       33]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
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
        self.assertEqual(alignment.shape, (2, 33))
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
                numpy.array([[53575980, 53575980, 53575997, 53575997],
                             [      33,       25,        8,        0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 50))
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
                numpy.array([[85737865, 85737865, 85737906],
                             [       0,        9,       50]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 50))
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
                numpy.array([[95160479, 95160479, 95160520, 95160520],
                             [       0,        8,       49,       50]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 50))
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
                numpy.array([[42144400, 42144400, 42144436, 42144436],
                             [       0,       11,       47,       50]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 50))
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
                numpy.array([[183925984, 183925984, 183925990, 183925990, 183926028, 183926028],
                             [        0,         1,         7,        11,        49,        50]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 184))
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
                numpy.array([[35483340, 35483340, 35483365, 35483499, 35483510, 35483510],
                             [       0,       10,       35,       35,       46,       50]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 50))
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
                numpy.array([[23891310, 23891310, 23891349, 23891349],
                             [       0,       10,       49,       50]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 50))
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
                numpy.array([[43252217, 43252217, 43252245, 43252245],
                             [       0,       21,       49,       50]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 53))
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
                numpy.array([[52759147, 52759147, 52759157, 52759160, 52759198, 52759198],
                             [       0,        1,       11,       11,       49,       50]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
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
                numpy.array([[61700837, 61700837, 61700871, 61700871],
                             [       0,        1,       35,       50]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 50))
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
                numpy.array([[37558157, 37558157, 37558173, 37558173, 37558191, 37558191],
                             [      50,       49,       33,       29,       11,        0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 50))
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
                numpy.array([[48997405, 48997405, 48997442, 48997442],
                             [      50,       49,       12,        0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 50))
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
                numpy.array([[120641740, 120641740, 120641776, 120641776],
                             [       50,        49,        13,         0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 50))
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
                numpy.array([[54017130, 54017130, 54017169, 54017169],
                             [      50,       49,       10,        0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 50))
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
                numpy.array([[553742, 553742, 553781, 553781],
                             [    50,     49,     10,      0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 50))
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
                numpy.array([[99388555, 99388555, 99388591, 99388591],
                             [      50,       49,       13,        0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 50))
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
                numpy.array([[112178171, 112178171, 112178196, 112178196],
                             [       50,        35,        10,         0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 50))
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
                numpy.array([[39368490, 39368490, 39368526, 39368526],
                             [      50,       49,       13,        0]]),
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 50))
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
                numpy.array([[220325687, 220325687, 220325721, 220325721],
                             [       50,        47,        13,         0]]),
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
