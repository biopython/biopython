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
    # The alignment scores (AS tag) were copied from the BED file dna_rna.bed.

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

    def check_alignments(self, alignments):
        """Check the alignments."""
        self.assertEqual(list(alignments.metadata), ["HD"])
        self.assertEqual(alignments.metadata["HD"], {"VN": "1.0", "SO": "unsorted"})
        self.assertEqual(len(alignments.targets), 25)
        self.assertEqual(alignments.targets["chr1"].id, "chr1")
        self.assertEqual(len(alignments.targets["chr1"]), 248956422)
        self.assertEqual(alignments.targets["chr1"].annotations, {"MD5": "2648ae1bacce4ec4b6cf337dcae37816", "assembly": "hg38", "species": "Homo sapiens"})
        self.assertEqual(alignments.targets["chr10"].id, "chr10")
        self.assertEqual(len(alignments.targets["chr10"]), 133797422)
        self.assertEqual(alignments.targets["chr10"].annotations, {"MD5": "907112d17fcb73bcab1ed1c72b97ce68", "assembly": "hg38", "species": "Homo sapiens"})
        self.assertEqual(alignments.targets["chr11"].id, "chr11")
        self.assertEqual(len(alignments.targets["chr11"]), 135086622)
        self.assertEqual(alignments.targets["chr11"].annotations, {"MD5": "1511375dc2dd1b633af8cf439ae90cec", "assembly": "hg38", "species": "Homo sapiens"})
        self.assertEqual(alignments.targets["chr12"].id, "chr12")
        self.assertEqual(len(alignments.targets["chr12"]), 133275309)
        self.assertEqual(alignments.targets["chr12"].annotations, {"MD5": "e81e16d3f44337034695a29b97708fce", "assembly": "hg38", "species": "Homo sapiens"})
        self.assertEqual(alignments.targets["chr13"].id, "chr13")
        self.assertEqual(len(alignments.targets["chr13"]), 114364328)
        self.assertEqual(alignments.targets["chr13"].annotations, {"MD5": "17dab79b963ccd8e7377cef59a54fe1c", "assembly": "hg38", "species": "Homo sapiens"})
        self.assertEqual(alignments.targets["chr14"].id, "chr14")
        self.assertEqual(len(alignments.targets["chr14"]), 107043718)
        self.assertEqual(alignments.targets["chr14"].annotations, {"MD5": "acbd9552c059d9b403e75ed26c1ce5bc", "assembly": "hg38", "species": "Homo sapiens"})
        self.assertEqual(alignments.targets["chr15"].id, "chr15")
        self.assertEqual(len(alignments.targets["chr15"]), 101991189)
        self.assertEqual(alignments.targets["chr15"].annotations, {"MD5": "f036bd11158407596ca6bf3581454706", "assembly": "hg38", "species": "Homo sapiens"})
        self.assertEqual(alignments.targets["chr16"].id, "chr16")
        self.assertEqual(len(alignments.targets["chr16"]), 90338345)
        self.assertEqual(alignments.targets["chr16"].annotations, {"MD5": "24e7cabfba3548a2bb4dff582b9ee870", "assembly": "hg38", "species": "Homo sapiens"})
        self.assertEqual(alignments.targets["chr17"].id, "chr17")
        self.assertEqual(len(alignments.targets["chr17"]), 83257441)
        self.assertEqual(alignments.targets["chr17"].annotations, {"MD5": "a8499ca51d6fb77332c2d242923994eb", "assembly": "hg38", "species": "Homo sapiens"})
        self.assertEqual(alignments.targets["chr18"].id, "chr18")
        self.assertEqual(len(alignments.targets["chr18"]), 80373285)
        self.assertEqual(alignments.targets["chr18"].annotations, {"MD5": "11eeaa801f6b0e2e36a1138616b8ee9a", "assembly": "hg38", "species": "Homo sapiens"})
        self.assertEqual(alignments.targets["chr19"].id, "chr19")
        self.assertEqual(len(alignments.targets["chr19"]), 58617616)
        self.assertEqual(alignments.targets["chr19"].annotations, {"MD5": "b0eba2c7bb5c953d1e06a508b5e487de", "assembly": "hg38", "species": "Homo sapiens"})
        self.assertEqual(alignments.targets["chr2"].id, "chr2")
        self.assertEqual(len(alignments.targets["chr2"]), 242193529)
        self.assertEqual(alignments.targets["chr2"].annotations, {"MD5": "4bb4f82880a14111eb7327169ffb729b", "assembly": "hg38", "species": "Homo sapiens"})
        self.assertEqual(alignments.targets["chr20"].id, "chr20")
        self.assertEqual(len(alignments.targets["chr20"]), 64444167)
        self.assertEqual(alignments.targets["chr20"].annotations, {"MD5": "b18e6c531b0bd70e949a7fc20859cb01", "assembly": "hg38", "species": "Homo sapiens"})
        self.assertEqual(alignments.targets["chr21"].id, "chr21")
        self.assertEqual(len(alignments.targets["chr21"]), 46709983)
        self.assertEqual(alignments.targets["chr21"].annotations, {"MD5": "2f45a3455007b7e271509161e52954a9", "assembly": "hg38", "species": "Homo sapiens"})
        self.assertEqual(alignments.targets["chr22"].id, "chr22")
        self.assertEqual(len(alignments.targets["chr22"]), 50818468)
        self.assertEqual(alignments.targets["chr22"].annotations, {"MD5": "221733a2a15e2de66d33e73d126c5109", "assembly": "hg38", "species": "Homo sapiens"})
        self.assertEqual(alignments.targets["chr3"].id, "chr3")
        self.assertEqual(len(alignments.targets["chr3"]), 198295559)
        self.assertEqual(alignments.targets["chr3"].annotations, {"MD5": "a48af509898d3736ba95dc0912c0b461", "assembly": "hg38", "species": "Homo sapiens"})
        self.assertEqual(alignments.targets["chr4"].id, "chr4")
        self.assertEqual(len(alignments.targets["chr4"]), 190214555)
        self.assertEqual(alignments.targets["chr4"].annotations, {"MD5": "3210fecf1eb92d5489da4346b3fddc6e", "assembly": "hg38", "species": "Homo sapiens"})
        self.assertEqual(alignments.targets["chr5"].id, "chr5")
        self.assertEqual(len(alignments.targets["chr5"]), 181538259)
        self.assertEqual(alignments.targets["chr5"].annotations, {"MD5": "f7f05fb7ceea78cbc32ce652c540ff2d", "assembly": "hg38", "species": "Homo sapiens"})
        self.assertEqual(alignments.targets["chr6"].id, "chr6")
        self.assertEqual(len(alignments.targets["chr6"]), 170805979)
        self.assertEqual(alignments.targets["chr6"].annotations, {"MD5": "6a48dfa97e854e3c6f186c8ff973f7dd", "assembly": "hg38", "species": "Homo sapiens"})
        self.assertEqual(alignments.targets["chr7"].id, "chr7")
        self.assertEqual(len(alignments.targets["chr7"]), 159345973)
        self.assertEqual(alignments.targets["chr7"].annotations, {"MD5": "94eef2b96fd5a7c8db162c8c74378039", "assembly": "hg38", "species": "Homo sapiens"})
        self.assertEqual(alignments.targets["chr8"].id, "chr8")
        self.assertEqual(len(alignments.targets["chr8"]), 145138636)
        self.assertEqual(alignments.targets["chr8"].annotations, {"MD5": "c67955b5f7815a9a1edfaa15893d3616", "assembly": "hg38", "species": "Homo sapiens"})
        self.assertEqual(alignments.targets["chr9"].id, "chr9")
        self.assertEqual(len(alignments.targets["chr9"]), 138394717)
        self.assertEqual(alignments.targets["chr9"].annotations, {"MD5": "addd2795560986b7491c40b1faa3978a", "assembly": "hg38", "species": "Homo sapiens"})
        self.assertEqual(alignments.targets["chrM"].id, "chrM")
        self.assertEqual(len(alignments.targets["chrM"]), 16569)
        self.assertEqual(alignments.targets["chrM"].annotations, {"MD5": "c68f52674c9fb33aef52dcf399755519", "assembly": "hg38", "species": "Homo sapiens"})
        self.assertEqual(alignments.targets["chrX"].id, "chrX")
        self.assertEqual(len(alignments.targets["chrX"]), 156040895)
        self.assertEqual(alignments.targets["chrX"].annotations, {"MD5": "49527016a48497d9d1cbd8e4a9049bd3", "assembly": "hg38", "species": "Homo sapiens"})
        self.assertEqual(alignments.targets["chrY"].id, "chrY")
        self.assertEqual(len(alignments.targets["chrY"]), 57227415)
        self.assertEqual(alignments.targets["chrY"].annotations, {"MD5": "b2b7e6369564d89059e763cd6e736837", "assembly": "hg38", "species": "Homo sapiens"})
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
        self.assertEqual(alignment.query.seq, self.rna[alignment.query.id])
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
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.annotations["NM"], 0)
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
        self.assertEqual(alignment.query.seq, self.rna[alignment.query.id])
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
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.annotations["NM"], 0)
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
        self.assertEqual(alignment.query.seq, self.rna[alignment.query.id])
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
        self.assertEqual(alignment.score, 972)
        self.assertEqual(alignment.annotations["NM"], 5)
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
        self.assertEqual(alignment.query.seq, self.rna[alignment.query.id])
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
        self.assertEqual(alignment.score, 978)
        self.assertEqual(alignment.annotations["NM"], 6)
        self.assertRaises(StopIteration, next, alignments)

    def test_reading(self):
        """Test parsing dna_rna.sam."""
        path = "Blat/dna_rna.sam"
        alignments = sam.AlignmentIterator(path)
        self.check_alignments(alignments)

    def test_writing(self):
        """Test writing the alignments in dna_rna.sam."""
        path = "Blat/dna_rna.sam"
        alignments = sam.AlignmentIterator(path)
        stream = StringIO()
        writer = sam.AlignmentWriter(stream)
        n = writer.write_file(alignments, mincount=4, maxcount=4)
        self.assertEqual(n, 4)
        stream.seek(0)
        alignments = sam.AlignmentIterator(stream)
        self.check_alignments(alignments)
        stream.close()


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

    def check_alignments_psl_34_001(self, alignments):
        """Check the alignments for psl_34_001/sam."""
        self.assertEqual(list(alignments.metadata), ["PG"])
        self.assertEqual(len(alignments.targets), 25)
        self.assertEqual(alignments.targets["chr1"].id, "chr1")
        self.assertEqual(len(alignments.targets["chr1"]), 249250621)
        self.assertEqual(alignments.targets["chr2"].id, "chr2")
        self.assertEqual(len(alignments.targets["chr2"]), 243199373)
        self.assertEqual(alignments.targets["chr3"].id, "chr3")
        self.assertEqual(len(alignments.targets["chr3"]), 198022430)
        self.assertEqual(alignments.targets["chr4"].id, "chr4")
        self.assertEqual(len(alignments.targets["chr4"]), 191154276)
        self.assertEqual(alignments.targets["chr5"].id, "chr5")
        self.assertEqual(len(alignments.targets["chr5"]), 180915260)
        self.assertEqual(alignments.targets["chr6"].id, "chr6")
        self.assertEqual(len(alignments.targets["chr6"]), 171115067)
        self.assertEqual(alignments.targets["chr7"].id, "chr7")
        self.assertEqual(len(alignments.targets["chr7"]), 159138663)
        self.assertEqual(alignments.targets["chrX"].id, "chrX")
        self.assertEqual(len(alignments.targets["chrX"]), 155270560)
        self.assertEqual(alignments.targets["chr8"].id, "chr8")
        self.assertEqual(len(alignments.targets["chr8"]), 146364022)
        self.assertEqual(alignments.targets["chr9"].id, "chr9")
        self.assertEqual(len(alignments.targets["chr9"]), 141213431)
        self.assertEqual(alignments.targets["chr10"].id, "chr10")
        self.assertEqual(len(alignments.targets["chr10"]), 135534747)
        self.assertEqual(alignments.targets["chr11"].id, "chr11")
        self.assertEqual(len(alignments.targets["chr11"]), 135006516)
        self.assertEqual(alignments.targets["chr12"].id, "chr12")
        self.assertEqual(len(alignments.targets["chr12"]), 133851895)
        self.assertEqual(alignments.targets["chr13"].id, "chr13")
        self.assertEqual(len(alignments.targets["chr13"]), 115169878)
        self.assertEqual(alignments.targets["chr14"].id, "chr14")
        self.assertEqual(len(alignments.targets["chr14"]), 107349540)
        self.assertEqual(alignments.targets["chr15"].id, "chr15")
        self.assertEqual(len(alignments.targets["chr15"]), 102531392)
        self.assertEqual(alignments.targets["chr16"].id, "chr16")
        self.assertEqual(len(alignments.targets["chr16"]), 90354753)
        self.assertEqual(alignments.targets["chr17"].id, "chr17")
        self.assertEqual(len(alignments.targets["chr17"]), 81195210)
        self.assertEqual(alignments.targets["chr18"].id, "chr18")
        self.assertEqual(len(alignments.targets["chr18"]), 78077248)
        self.assertEqual(alignments.targets["chr20"].id, "chr20")
        self.assertEqual(len(alignments.targets["chr20"]), 63025520)
        self.assertEqual(alignments.targets["chrY"].id, "chrY")
        self.assertEqual(len(alignments.targets["chrY"]), 59373566)
        self.assertEqual(alignments.targets["chr19"].id, "chr19")
        self.assertEqual(len(alignments.targets["chr19"]), 59128983)
        self.assertEqual(alignments.targets["chr22"].id, "chr22")
        self.assertEqual(len(alignments.targets["chr22"]), 51304566)
        self.assertEqual(alignments.targets["chr21"].id, "chr21")
        self.assertEqual(len(alignments.targets["chr21"]), 48129895)
        self.assertEqual(alignments.targets["chrM"].id, "chrM")
        self.assertEqual(len(alignments.targets["chrM"]), 16571)
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

    def test_reading_psl_34_001(self):
        """Test parsing psl_34_001.sam."""
        path = "Blat/psl_34_001.sam"
        alignments = sam.AlignmentIterator(path)
        self.check_alignments_psl_34_001(alignments)

    def test_writing_psl_34_001(self):
        """Test writing the alignments in psl_34_001.sam."""
        path = "Blat/psl_34_001.sam"
        alignments = sam.AlignmentIterator(path)
        stream = StringIO()
        writer = sam.AlignmentWriter(stream)
        n = writer.write_file(alignments, mincount=22, maxcount=22)
        self.assertEqual(n, 22)
        stream.seek(0)
        alignments = sam.AlignmentIterator(stream)
        self.check_alignments_psl_34_001(alignments)
        stream.close()

    def check_alignments_psl_34_003(self, alignments):
        """Check the alignments for psl_34_003/sam."""
        self.assertEqual(list(alignments.metadata), ["PG"])
        self.assertEqual(len(alignments.targets), 25)
        self.assertEqual(alignments.targets["chr1"].id, "chr1")
        self.assertEqual(len(alignments.targets["chr1"]), 249250621)
        self.assertEqual(alignments.targets["chr2"].id, "chr2")
        self.assertEqual(len(alignments.targets["chr2"]), 243199373)
        self.assertEqual(alignments.targets["chr3"].id, "chr3")
        self.assertEqual(len(alignments.targets["chr3"]), 198022430)
        self.assertEqual(alignments.targets["chr4"].id, "chr4")
        self.assertEqual(len(alignments.targets["chr4"]), 191154276)
        self.assertEqual(alignments.targets["chr5"].id, "chr5")
        self.assertEqual(len(alignments.targets["chr5"]), 180915260)
        self.assertEqual(alignments.targets["chr6"].id, "chr6")
        self.assertEqual(len(alignments.targets["chr6"]), 171115067)
        self.assertEqual(alignments.targets["chr7"].id, "chr7")
        self.assertEqual(len(alignments.targets["chr7"]), 159138663)
        self.assertEqual(alignments.targets["chrX"].id, "chrX")
        self.assertEqual(len(alignments.targets["chrX"]), 155270560)
        self.assertEqual(alignments.targets["chr8"].id, "chr8")
        self.assertEqual(len(alignments.targets["chr8"]), 146364022)
        self.assertEqual(alignments.targets["chr9"].id, "chr9")
        self.assertEqual(len(alignments.targets["chr9"]), 141213431)
        self.assertEqual(alignments.targets["chr10"].id, "chr10")
        self.assertEqual(len(alignments.targets["chr10"]), 135534747)
        self.assertEqual(alignments.targets["chr11"].id, "chr11")
        self.assertEqual(len(alignments.targets["chr11"]), 135006516)
        self.assertEqual(alignments.targets["chr12"].id, "chr12")
        self.assertEqual(len(alignments.targets["chr12"]), 133851895)
        self.assertEqual(alignments.targets["chr13"].id, "chr13")
        self.assertEqual(len(alignments.targets["chr13"]), 115169878)
        self.assertEqual(alignments.targets["chr14"].id, "chr14")
        self.assertEqual(len(alignments.targets["chr14"]), 107349540)
        self.assertEqual(alignments.targets["chr15"].id, "chr15")
        self.assertEqual(len(alignments.targets["chr15"]), 102531392)
        self.assertEqual(alignments.targets["chr16"].id, "chr16")
        self.assertEqual(len(alignments.targets["chr16"]), 90354753)
        self.assertEqual(alignments.targets["chr17"].id, "chr17")
        self.assertEqual(len(alignments.targets["chr17"]), 81195210)
        self.assertEqual(alignments.targets["chr18"].id, "chr18")
        self.assertEqual(len(alignments.targets["chr18"]), 78077248)
        self.assertEqual(alignments.targets["chr20"].id, "chr20")
        self.assertEqual(len(alignments.targets["chr20"]), 63025520)
        self.assertEqual(alignments.targets["chrY"].id, "chrY")
        self.assertEqual(len(alignments.targets["chrY"]), 59373566)
        self.assertEqual(alignments.targets["chr19"].id, "chr19")
        self.assertEqual(len(alignments.targets["chr19"]), 59128983)
        self.assertEqual(alignments.targets["chr22"].id, "chr22")
        self.assertEqual(len(alignments.targets["chr22"]), 51304566)
        self.assertEqual(alignments.targets["chr21"].id, "chr21")
        self.assertEqual(len(alignments.targets["chr21"]), 48129895)
        self.assertEqual(alignments.targets["chrM"].id, "chrM")
        self.assertEqual(len(alignments.targets["chrM"]), 16571)
        self.assertEqual(len(alignments.targets), 25)
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

    def test_reading_psl_34_003(self):
        """Test parsing psl_34_003.sam."""
        path = "Blat/psl_34_003.sam"
        alignments = sam.AlignmentIterator(path)
        self.check_alignments_psl_34_003(alignments)

    def test_writing_psl_34_003(self):
        """Test writing the alignments in psl_34_003.sam."""
        path = "Blat/psl_34_003.sam"
        alignments = sam.AlignmentIterator(path)
        stream = StringIO()
        writer = sam.AlignmentWriter(stream)
        n = writer.write_file(alignments, mincount=3, maxcount=3)
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = sam.AlignmentIterator(stream)
        self.check_alignments_psl_34_003(alignments)
        stream.close()

    def check_alignments_psl_34_004(self, alignments):
        """Check the alignments for psl_34_004/sam."""
        self.assertEqual(list(alignments.metadata), ["PG"])
        self.assertEqual(len(alignments.targets), 25)
        self.assertEqual(alignments.targets["chr1"].id, "chr1")
        self.assertEqual(len(alignments.targets["chr1"]), 249250621)
        self.assertEqual(alignments.targets["chr2"].id, "chr2")
        self.assertEqual(len(alignments.targets["chr2"]), 243199373)
        self.assertEqual(alignments.targets["chr3"].id, "chr3")
        self.assertEqual(len(alignments.targets["chr3"]), 198022430)
        self.assertEqual(alignments.targets["chr4"].id, "chr4")
        self.assertEqual(len(alignments.targets["chr4"]), 191154276)
        self.assertEqual(alignments.targets["chr5"].id, "chr5")
        self.assertEqual(len(alignments.targets["chr5"]), 180915260)
        self.assertEqual(alignments.targets["chr6"].id, "chr6")
        self.assertEqual(len(alignments.targets["chr6"]), 171115067)
        self.assertEqual(alignments.targets["chr7"].id, "chr7")
        self.assertEqual(len(alignments.targets["chr7"]), 159138663)
        self.assertEqual(alignments.targets["chrX"].id, "chrX")
        self.assertEqual(len(alignments.targets["chrX"]), 155270560)
        self.assertEqual(alignments.targets["chr8"].id, "chr8")
        self.assertEqual(len(alignments.targets["chr8"]), 146364022)
        self.assertEqual(alignments.targets["chr9"].id, "chr9")
        self.assertEqual(len(alignments.targets["chr9"]), 141213431)
        self.assertEqual(alignments.targets["chr10"].id, "chr10")
        self.assertEqual(len(alignments.targets["chr10"]), 135534747)
        self.assertEqual(alignments.targets["chr11"].id, "chr11")
        self.assertEqual(len(alignments.targets["chr11"]), 135006516)
        self.assertEqual(alignments.targets["chr12"].id, "chr12")
        self.assertEqual(len(alignments.targets["chr12"]), 133851895)
        self.assertEqual(alignments.targets["chr13"].id, "chr13")
        self.assertEqual(len(alignments.targets["chr13"]), 115169878)
        self.assertEqual(alignments.targets["chr14"].id, "chr14")
        self.assertEqual(len(alignments.targets["chr14"]), 107349540)
        self.assertEqual(alignments.targets["chr15"].id, "chr15")
        self.assertEqual(len(alignments.targets["chr15"]), 102531392)
        self.assertEqual(alignments.targets["chr16"].id, "chr16")
        self.assertEqual(len(alignments.targets["chr16"]), 90354753)
        self.assertEqual(alignments.targets["chr17"].id, "chr17")
        self.assertEqual(len(alignments.targets["chr17"]), 81195210)
        self.assertEqual(alignments.targets["chr18"].id, "chr18")
        self.assertEqual(len(alignments.targets["chr18"]), 78077248)
        self.assertEqual(alignments.targets["chr20"].id, "chr20")
        self.assertEqual(len(alignments.targets["chr20"]), 63025520)
        self.assertEqual(alignments.targets["chrY"].id, "chrY")
        self.assertEqual(len(alignments.targets["chrY"]), 59373566)
        self.assertEqual(alignments.targets["chr19"].id, "chr19")
        self.assertEqual(len(alignments.targets["chr19"]), 59128983)
        self.assertEqual(alignments.targets["chr22"].id, "chr22")
        self.assertEqual(len(alignments.targets["chr22"]), 51304566)
        self.assertEqual(alignments.targets["chr21"].id, "chr21")
        self.assertEqual(len(alignments.targets["chr21"]), 48129895)
        self.assertEqual(alignments.targets["chrM"].id, "chrM")
        self.assertEqual(len(alignments.targets["chrM"]), 16571)
        self.assertEqual(len(alignments.targets), 25)
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

    def test_reading_psl_34_004(self):
        """Test parsing psl_34_004.sam."""
        path = "Blat/psl_34_004.sam"
        alignments = sam.AlignmentIterator(path)
        self.check_alignments_psl_34_004(alignments)

    def test_writing_psl_34_004(self):
        """Test writing the alignments in psl_34_004.sam."""
        path = "Blat/psl_34_004.sam"
        alignments = sam.AlignmentIterator(path)
        stream = StringIO()
        writer = sam.AlignmentWriter(stream)
        n = writer.write_file(alignments, mincount=19, maxcount=19)
        self.assertEqual(n, 19)
        stream.seek(0)
        alignments = sam.AlignmentIterator(stream)
        self.check_alignments_psl_34_004(alignments)
        stream.close()

    def check_alignments_psl_34_005(self, alignments):
        """Check the alignments for psl_34_005/sam."""
        self.assertEqual(list(alignments.metadata), ["PG"])
        self.assertEqual(len(alignments.targets), 25)
        self.assertEqual(alignments.targets["chr1"].id, "chr1")
        self.assertEqual(len(alignments.targets["chr1"]), 249250621)
        self.assertEqual(alignments.targets["chr2"].id, "chr2")
        self.assertEqual(len(alignments.targets["chr2"]), 243199373)
        self.assertEqual(alignments.targets["chr3"].id, "chr3")
        self.assertEqual(len(alignments.targets["chr3"]), 198022430)
        self.assertEqual(alignments.targets["chr4"].id, "chr4")
        self.assertEqual(len(alignments.targets["chr4"]), 191154276)
        self.assertEqual(alignments.targets["chr5"].id, "chr5")
        self.assertEqual(len(alignments.targets["chr5"]), 180915260)
        self.assertEqual(alignments.targets["chr6"].id, "chr6")
        self.assertEqual(len(alignments.targets["chr6"]), 171115067)
        self.assertEqual(alignments.targets["chr7"].id, "chr7")
        self.assertEqual(len(alignments.targets["chr7"]), 159138663)
        self.assertEqual(alignments.targets["chrX"].id, "chrX")
        self.assertEqual(len(alignments.targets["chrX"]), 155270560)
        self.assertEqual(alignments.targets["chr8"].id, "chr8")
        self.assertEqual(len(alignments.targets["chr8"]), 146364022)
        self.assertEqual(alignments.targets["chr9"].id, "chr9")
        self.assertEqual(len(alignments.targets["chr9"]), 141213431)
        self.assertEqual(alignments.targets["chr10"].id, "chr10")
        self.assertEqual(len(alignments.targets["chr10"]), 135534747)
        self.assertEqual(alignments.targets["chr11"].id, "chr11")
        self.assertEqual(len(alignments.targets["chr11"]), 135006516)
        self.assertEqual(alignments.targets["chr12"].id, "chr12")
        self.assertEqual(len(alignments.targets["chr12"]), 133851895)
        self.assertEqual(alignments.targets["chr13"].id, "chr13")
        self.assertEqual(len(alignments.targets["chr13"]), 115169878)
        self.assertEqual(alignments.targets["chr14"].id, "chr14")
        self.assertEqual(len(alignments.targets["chr14"]), 107349540)
        self.assertEqual(alignments.targets["chr15"].id, "chr15")
        self.assertEqual(len(alignments.targets["chr15"]), 102531392)
        self.assertEqual(alignments.targets["chr16"].id, "chr16")
        self.assertEqual(len(alignments.targets["chr16"]), 90354753)
        self.assertEqual(alignments.targets["chr17"].id, "chr17")
        self.assertEqual(len(alignments.targets["chr17"]), 81195210)
        self.assertEqual(alignments.targets["chr18"].id, "chr18")
        self.assertEqual(len(alignments.targets["chr18"]), 78077248)
        self.assertEqual(alignments.targets["chr20"].id, "chr20")
        self.assertEqual(len(alignments.targets["chr20"]), 63025520)
        self.assertEqual(alignments.targets["chrY"].id, "chrY")
        self.assertEqual(len(alignments.targets["chrY"]), 59373566)
        self.assertEqual(alignments.targets["chr19"].id, "chr19")
        self.assertEqual(len(alignments.targets["chr19"]), 59128983)
        self.assertEqual(alignments.targets["chr22"].id, "chr22")
        self.assertEqual(len(alignments.targets["chr22"]), 51304566)
        self.assertEqual(alignments.targets["chr21"].id, "chr21")
        self.assertEqual(len(alignments.targets["chr21"]), 48129895)
        self.assertEqual(alignments.targets["chrM"].id, "chrM")
        self.assertEqual(len(alignments.targets["chrM"]), 16571)
        self.assertEqual(len(alignments.targets), 25)
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

    def test_reading_psl_34_005(self):
        """Test parsing psl_34_005.sam."""
        path = "Blat/psl_34_005.sam"
        alignments = sam.AlignmentIterator(path)
        self.check_alignments_psl_34_005(alignments)

    def test_writing_psl_34_005(self):
        """Test writing the alignments in psl_34_005.sam."""
        path = "Blat/psl_34_005.sam"
        alignments = sam.AlignmentIterator(path)
        stream = StringIO()
        writer = sam.AlignmentWriter(stream)
        n = writer.write_file(alignments, mincount=22, maxcount=22)
        self.assertEqual(n, 22)
        stream.seek(0)
        alignments = sam.AlignmentIterator(stream)
        self.check_alignments_psl_34_005(alignments)
        stream.close()


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
