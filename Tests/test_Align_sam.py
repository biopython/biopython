# Copyright 2022 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Align.sam module."""
import unittest
import warnings
from io import StringIO


from Bio.Align import Alignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import Align


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
            sequence = str(record.seq).upper()
            assert len(sequence) == end - start
            data[start] = sequence
        self.dna = Seq(data, length=198295559)
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
        self.assertEqual(alignments.targets[0].id, "chr1")
        self.assertEqual(len(alignments.targets[0]), 248956422)
        self.assertEqual(
            alignments.targets[0].annotations,
            {
                "MD5": "2648ae1bacce4ec4b6cf337dcae37816",
                "assembly": "hg38",
                "species": "Homo sapiens",
            },
        )
        self.assertEqual(alignments.targets[1].id, "chr10")
        self.assertEqual(len(alignments.targets[1]), 133797422)
        self.assertEqual(
            alignments.targets[1].annotations,
            {
                "MD5": "907112d17fcb73bcab1ed1c72b97ce68",
                "assembly": "hg38",
                "species": "Homo sapiens",
            },
        )
        self.assertEqual(alignments.targets[2].id, "chr11")
        self.assertEqual(len(alignments.targets[2]), 135086622)
        self.assertEqual(
            alignments.targets[2].annotations,
            {
                "MD5": "1511375dc2dd1b633af8cf439ae90cec",
                "assembly": "hg38",
                "species": "Homo sapiens",
            },
        )
        self.assertEqual(alignments.targets[3].id, "chr12")
        self.assertEqual(len(alignments.targets[3]), 133275309)
        self.assertEqual(
            alignments.targets[3].annotations,
            {
                "MD5": "e81e16d3f44337034695a29b97708fce",
                "assembly": "hg38",
                "species": "Homo sapiens",
            },
        )
        self.assertEqual(alignments.targets[4].id, "chr13")
        self.assertEqual(len(alignments.targets[4]), 114364328)
        self.assertEqual(
            alignments.targets[4].annotations,
            {
                "MD5": "17dab79b963ccd8e7377cef59a54fe1c",
                "assembly": "hg38",
                "species": "Homo sapiens",
            },
        )
        self.assertEqual(alignments.targets[5].id, "chr14")
        self.assertEqual(len(alignments.targets[5]), 107043718)
        self.assertEqual(
            alignments.targets[5].annotations,
            {
                "MD5": "acbd9552c059d9b403e75ed26c1ce5bc",
                "assembly": "hg38",
                "species": "Homo sapiens",
            },
        )
        self.assertEqual(alignments.targets[6].id, "chr15")
        self.assertEqual(len(alignments.targets[6]), 101991189)
        self.assertEqual(
            alignments.targets[6].annotations,
            {
                "MD5": "f036bd11158407596ca6bf3581454706",
                "assembly": "hg38",
                "species": "Homo sapiens",
            },
        )
        self.assertEqual(alignments.targets[7].id, "chr16")
        self.assertEqual(len(alignments.targets[7]), 90338345)
        self.assertEqual(
            alignments.targets[7].annotations,
            {
                "MD5": "24e7cabfba3548a2bb4dff582b9ee870",
                "assembly": "hg38",
                "species": "Homo sapiens",
            },
        )
        self.assertEqual(alignments.targets[8].id, "chr17")
        self.assertEqual(len(alignments.targets[8]), 83257441)
        self.assertEqual(
            alignments.targets[8].annotations,
            {
                "MD5": "a8499ca51d6fb77332c2d242923994eb",
                "assembly": "hg38",
                "species": "Homo sapiens",
            },
        )
        self.assertEqual(alignments.targets[9].id, "chr18")
        self.assertEqual(len(alignments.targets[9]), 80373285)
        self.assertEqual(
            alignments.targets[9].annotations,
            {
                "MD5": "11eeaa801f6b0e2e36a1138616b8ee9a",
                "assembly": "hg38",
                "species": "Homo sapiens",
            },
        )
        self.assertEqual(alignments.targets[10].id, "chr19")
        self.assertEqual(len(alignments.targets[10]), 58617616)
        self.assertEqual(
            alignments.targets[10].annotations,
            {
                "MD5": "b0eba2c7bb5c953d1e06a508b5e487de",
                "assembly": "hg38",
                "species": "Homo sapiens",
            },
        )
        self.assertEqual(alignments.targets[11].id, "chr2")
        self.assertEqual(len(alignments.targets[11]), 242193529)
        self.assertEqual(
            alignments.targets[11].annotations,
            {
                "MD5": "4bb4f82880a14111eb7327169ffb729b",
                "assembly": "hg38",
                "species": "Homo sapiens",
            },
        )
        self.assertEqual(alignments.targets[12].id, "chr20")
        self.assertEqual(len(alignments.targets[12]), 64444167)
        self.assertEqual(
            alignments.targets[12].annotations,
            {
                "MD5": "b18e6c531b0bd70e949a7fc20859cb01",
                "assembly": "hg38",
                "species": "Homo sapiens",
            },
        )
        self.assertEqual(alignments.targets[13].id, "chr21")
        self.assertEqual(len(alignments.targets[13]), 46709983)
        self.assertEqual(
            alignments.targets[13].annotations,
            {
                "MD5": "2f45a3455007b7e271509161e52954a9",
                "assembly": "hg38",
                "species": "Homo sapiens",
            },
        )
        self.assertEqual(alignments.targets[14].id, "chr22")
        self.assertEqual(len(alignments.targets[14]), 50818468)
        self.assertEqual(
            alignments.targets[14].annotations,
            {
                "MD5": "221733a2a15e2de66d33e73d126c5109",
                "assembly": "hg38",
                "species": "Homo sapiens",
            },
        )
        self.assertEqual(alignments.targets[15].id, "chr3")
        self.assertEqual(len(alignments.targets[15]), 198295559)
        self.assertEqual(
            alignments.targets[15].annotations,
            {
                "MD5": "a48af509898d3736ba95dc0912c0b461",
                "assembly": "hg38",
                "species": "Homo sapiens",
            },
        )
        self.assertEqual(alignments.targets[16].id, "chr4")
        self.assertEqual(len(alignments.targets[16]), 190214555)
        self.assertEqual(
            alignments.targets[16].annotations,
            {
                "MD5": "3210fecf1eb92d5489da4346b3fddc6e",
                "assembly": "hg38",
                "species": "Homo sapiens",
            },
        )
        self.assertEqual(alignments.targets[17].id, "chr5")
        self.assertEqual(len(alignments.targets[17]), 181538259)
        self.assertEqual(
            alignments.targets[17].annotations,
            {
                "MD5": "f7f05fb7ceea78cbc32ce652c540ff2d",
                "assembly": "hg38",
                "species": "Homo sapiens",
            },
        )
        self.assertEqual(alignments.targets[18].id, "chr6")
        self.assertEqual(len(alignments.targets[18]), 170805979)
        self.assertEqual(
            alignments.targets[18].annotations,
            {
                "MD5": "6a48dfa97e854e3c6f186c8ff973f7dd",
                "assembly": "hg38",
                "species": "Homo sapiens",
            },
        )
        self.assertEqual(alignments.targets[19].id, "chr7")
        self.assertEqual(len(alignments.targets[19]), 159345973)
        self.assertEqual(
            alignments.targets[19].annotations,
            {
                "MD5": "94eef2b96fd5a7c8db162c8c74378039",
                "assembly": "hg38",
                "species": "Homo sapiens",
            },
        )
        self.assertEqual(alignments.targets[20].id, "chr8")
        self.assertEqual(len(alignments.targets[20]), 145138636)
        self.assertEqual(
            alignments.targets[20].annotations,
            {
                "MD5": "c67955b5f7815a9a1edfaa15893d3616",
                "assembly": "hg38",
                "species": "Homo sapiens",
            },
        )
        self.assertEqual(alignments.targets[21].id, "chr9")
        self.assertEqual(len(alignments.targets[21]), 138394717)
        self.assertEqual(
            alignments.targets[21].annotations,
            {
                "MD5": "addd2795560986b7491c40b1faa3978a",
                "assembly": "hg38",
                "species": "Homo sapiens",
            },
        )
        self.assertEqual(alignments.targets[22].id, "chrM")
        self.assertEqual(len(alignments.targets[22]), 16569)
        self.assertEqual(
            alignments.targets[22].annotations,
            {
                "MD5": "c68f52674c9fb33aef52dcf399755519",
                "assembly": "hg38",
                "species": "Homo sapiens",
            },
        )
        self.assertEqual(alignments.targets[23].id, "chrX")
        self.assertEqual(len(alignments.targets[23]), 156040895)
        self.assertEqual(
            alignments.targets[23].annotations,
            {
                "MD5": "49527016a48497d9d1cbd8e4a9049bd3",
                "assembly": "hg38",
                "species": "Homo sapiens",
            },
        )
        self.assertEqual(alignments.targets[24].id, "chrY")
        self.assertEqual(len(alignments.targets[24]), 57227415)
        self.assertEqual(
            alignments.targets[24].annotations,
            {
                "MD5": "b2b7e6369564d89059e763cd6e736837",
                "assembly": "hg38",
                "species": "Homo sapiens",
            },
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 5407))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr3")
        self.assertEqual(alignment.query.id, "NR_111921.1")
        self.assertEqual(len(alignment.target.seq), len(self.dna))
        self.assertEqual(
            alignment.target.seq.defined_ranges,
            ((48663767, 48663813), (48665640, 48665722), (48669098, 48669174)),
        )
        for start, end in alignment.target.seq.defined_ranges:
            self.assertEqual(alignment.target.seq[start:end], self.dna[start:end])
        self.assertEqual(alignment.query.seq, self.rna[alignment.query.id])
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
        self.assertTrue(
            numpy.array_equal(
                alignment.substitutions,
                # fmt: off
                numpy.array([[62.,  0.,  0.,  0.],
                             [ 0., 42.,  0.,  0.],
                             [ 0.,  0., 66.,  0.],
                             [ 0.,  0.,  0., 34.],
                            ])
                # fmt: on
            )
        )
        self.assertEqual(alignment.substitutions.alphabet, "ACGT")
        self.assertEqual(alignment.mapq, 0)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.annotations["NM"], 0)
        self.assertNotIn("hard_clip_left", alignment.query.annotations)
        self.assertEqual(alignment.query.annotations["hard_clip_right"], 12)
        self.assertEqual(alignment.operations, bytearray(b"MNMNM"))
        self.assertEqual(
            format(alignment, "sam"),
            """\
NR_111921.1	0	chr3	48663768	0	46M1827N82M3376N76M12H	*	0	0	CACGAGAGGAGCGGAGGCGAGGGGTGAACGCGGAGCACTCCAATCGCTCCCAACTAGAGGTCCACCCAGGACCCAGAGACCTGGATTTGAGGCTGCTGGGCGGCAGATGGAGCGATCAGAAGACCAGGAGACGGGAGCTGGAGTGCAGTGGCTGTTCACAAGCGTGAAAGCAAAGATTAAAAAATTTGTTTTTATATTAAAAAA	*	AS:i:1000	NM:i:0
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 1711))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr3")
        self.assertEqual(alignment.query.id, "NR_046654.1")
        self.assertEqual(len(alignment.target.seq), len(self.dna))
        self.assertEqual(
            alignment.target.seq.defined_ranges,
            ((42530895, 42530958), (42532020, 42532095), (42532563, 42532606)),
        )
        for start, end in alignment.target.seq.defined_ranges:
            self.assertEqual(alignment.target.seq[start:end], self.dna[start:end])
        self.assertEqual(alignment.query.seq, self.rna[alignment.query.id])
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
        self.assertTrue(
            numpy.array_equal(
                alignment.substitutions,
                # fmt: off
# flake8: noqa
                numpy.array([[38.,  0.,  0.,  0.],
                             [ 0., 41.,  0.,  0.],
                             [ 0.,  0., 60.,  0.],
                             [ 0.,  0.,  0., 42.],
                            ])
                # fmt: on
            )
        )
        self.assertEqual(alignment.substitutions.alphabet, "ACGT")
        self.assertEqual(alignment.mapq, 0)
        matches = sum(
            alignment.substitutions[c, c] for c in alignment.substitutions.alphabet
        )
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.annotations["NM"], 0)
        self.assertNotIn("hard_clip_left", alignment.query.annotations)
        self.assertNotIn("hard_clip_right", alignment.query.annotations)
        self.assertEqual(alignment.operations, bytearray(b"MNMNM"))
        self.assertEqual(
            format(alignment, "sam"),
            """\
NR_046654.1	16	chr3	42530896	0	63M1062N75M468N43M	*	0	0	CGGAAGTACTTCTGGGGGTACATACTCATCGGCTGGGGTATGGTACCAGGGAGGGCTTCCAGGCAGTTCTTCCTTGAGCGTAAGCGGATTGGGAGCACAGTCCTTAGGGATTTGAAGGAGGTAGAGTTCCCGGATGACCTAGCATCCTTCCCAGGTATGCATCTGCTGCCAAGCCAGGGAG	*	AS:i:1000	NM:i:0
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 5409))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr3")
        self.assertEqual(alignment.query.id, "NR_111921.1_modified")
        self.assertEqual(len(alignment.target.seq), len(self.dna))
        self.assertEqual(
            alignment.target.seq.defined_ranges,
            ((48663767, 48663813), (48665640, 48665722), (48669098, 48669174)),
        )
        for start, end in alignment.target.seq.defined_ranges:
            self.assertEqual(alignment.target.seq[start:end], self.dna[start:end])
        self.assertEqual(alignment.query.seq, self.rna[alignment.query.id])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[48663767, 48663795, 48663796, 48663813, 48665640,
                              48665716, 48665716, 48665722, 48669098, 48669174],
                             [       3,       31,       31,       48,       48,
                                   124,      126,      132,      132,      208],
                            ])
                # fmt: on
            )
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.substitutions,
                # fmt: off
                numpy.array([[62.,  0.,  0.,  0.],
                             [ 0., 41.,  0.,  0.],
                             [ 0.,  2., 64.,  0.],
                             [ 0.,  0.,  0., 34.],
                            ]),
                # fmt: on
            )
        )
        self.assertEqual(alignment.substitutions.alphabet, "ACGT")
        self.assertEqual(alignment.mapq, 0)
        self.assertEqual(alignment.score, 972)
        self.assertEqual(alignment.annotations["NM"], 5)
        self.assertNotIn("hard_clip_left", alignment.query.annotations)
        self.assertEqual(alignment.query.annotations["hard_clip_right"], 12)
        self.assertEqual(alignment.operations, bytearray(b"MDMNMIMNM"))
        self.assertEqual(
            format(alignment, "sam"),
            """\
NR_111921.1_modified	0	chr3	48663768	0	3S28M1D17M1827N76M2I6M3376N76M12H	*	0	0	AAACACGAGAGGAGCGGAGGCGAGGGGTGAAGCGGAGCACTCCAATCGCTCCCAACTAGAGGTCCACCCAGGACCCAGAGACCTGGATTTGAGGCTGCTGCCCGGCAGATGGAGCGATCAGAAGCCACCAGGAGACGGGAGCTGGAGTGCAGTGGCTGTTCACAAGCGTGAAAGCAAAGATTAAAAAATTTGTTTTTATATTAAAAAA	*	AS:i:972	NM:i:5
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (2, 1714))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr3")
        self.assertEqual(alignment.query.id, "NR_046654.1_modified")
        self.assertEqual(len(alignment.target.seq), len(self.dna))
        self.assertEqual(
            alignment.target.seq.defined_ranges,
            ((42530895, 42530958), (42532020, 42532095), (42532563, 42532606)),
        )
        for start, end in alignment.target.seq.defined_ranges:
            self.assertEqual(alignment.target.seq[start:end], self.dna[start:end])
        self.assertEqual(alignment.query.seq, self.rna[alignment.query.id])
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
        self.assertTrue(
            numpy.array_equal(
                alignment.substitutions,
                # fmt: off
                numpy.array([[36.,  0.,  0.,  1.],
                             [ 0., 41.,  0.,  0.],
                             [ 0.,  0., 60.,  0.],
                             [ 0.,  0.,  0., 41.],
                            ]),
                # fmt: on
            )
        )
        self.assertEqual(alignment.substitutions.alphabet, "ACGT")
        self.assertEqual(alignment.mapq, 0)
        self.assertEqual(alignment.score, 978)
        self.assertEqual(alignment.annotations["NM"], 6)
        self.assertNotIn("hard_clip_left", alignment.query.annotations)
        self.assertNotIn("hard_clip_right", alignment.query.annotations)
        self.assertEqual(alignment.operations, bytearray(b"MIMNMDMNM"))
        self.assertEqual(
            format(alignment, "sam"),
            """\
NR_046654.1_modified	16	chr3	42530896	0	5S27M3I36M1062N17M2D56M468N43M3S	*	0	0	AAAAACGGAAGTACTTCTGGGGGTACATACTCCCCATCGGCTGGGGTATGGTACCAGGGAGGGCTTCCAGGCAGTTCTTCCTTGAGCGAGCGGATTGGGTGCACAGTCCTTAGGGATTTGAAGGAGGTAGAGTTCCCGGATGACCTAGCATCCTTCCCAGGTATGCATCTGCTGCCAAGCCAGGGAGAAA	*	AS:i:978	NM:i:6
""",
        )
        self.assertRaises(StopIteration, next, alignments)

    def test_reading(self):
        """Test parsing dna_rna.sam."""
        path = "Blat/dna_rna.sam"
        alignments = Align.parse(path, "sam")
        self.check_alignments(alignments)

    def test_reading_psl_comparison(self):
        """Test parsing dna_rna.sam and comparing to dna_rna.psl."""
        path = "Blat/dna_rna.sam"
        sam_alignments = Align.parse(path, "sam")
        path = "Blat/dna_rna.psl"
        psl_alignments = Align.parse(path, "psl")
        for sam_alignment, psl_alignment in zip(sam_alignments, psl_alignments):
            self.assertEqual(sam_alignment.target.id, psl_alignment.target.id)
            self.assertEqual(sam_alignment.query.id, psl_alignment.query.id)
            self.assertTrue(
                numpy.array_equal(sam_alignment.coordinates, psl_alignment.coordinates)
            )

    def test_writing(self):
        """Test writing the alignments in dna_rna.sam."""
        path = "Blat/dna_rna.sam"
        alignments = Align.parse(path, "sam")
        stream = StringIO()
        n = Align.write(alignments, stream, "sam", md=True)
        self.assertEqual(n, 4)
        stream.seek(0)
        alignments = Align.parse(stream, "sam")
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
        self.assertEqual(alignments.targets[0].id, "chr1")
        self.assertEqual(len(alignments.targets[0]), 249250621)
        self.assertEqual(alignments.targets[1].id, "chr2")
        self.assertEqual(len(alignments.targets[1]), 243199373)
        self.assertEqual(alignments.targets[2].id, "chr3")
        self.assertEqual(len(alignments.targets[2]), 198022430)
        self.assertEqual(alignments.targets[3].id, "chr4")
        self.assertEqual(len(alignments.targets[3]), 191154276)
        self.assertEqual(alignments.targets[4].id, "chr5")
        self.assertEqual(len(alignments.targets[4]), 180915260)
        self.assertEqual(alignments.targets[5].id, "chr6")
        self.assertEqual(len(alignments.targets[5]), 171115067)
        self.assertEqual(alignments.targets[6].id, "chr7")
        self.assertEqual(len(alignments.targets[6]), 159138663)
        self.assertEqual(alignments.targets[7].id, "chrX")
        self.assertEqual(len(alignments.targets[7]), 155270560)
        self.assertEqual(alignments.targets[8].id, "chr8")
        self.assertEqual(len(alignments.targets[8]), 146364022)
        self.assertEqual(alignments.targets[9].id, "chr9")
        self.assertEqual(len(alignments.targets[9]), 141213431)
        self.assertEqual(alignments.targets[10].id, "chr10")
        self.assertEqual(len(alignments.targets[10]), 135534747)
        self.assertEqual(alignments.targets[11].id, "chr11")
        self.assertEqual(len(alignments.targets[11]), 135006516)
        self.assertEqual(alignments.targets[12].id, "chr12")
        self.assertEqual(len(alignments.targets[12]), 133851895)
        self.assertEqual(alignments.targets[13].id, "chr13")
        self.assertEqual(len(alignments.targets[13]), 115169878)
        self.assertEqual(alignments.targets[14].id, "chr14")
        self.assertEqual(len(alignments.targets[14]), 107349540)
        self.assertEqual(alignments.targets[15].id, "chr15")
        self.assertEqual(len(alignments.targets[15]), 102531392)
        self.assertEqual(alignments.targets[16].id, "chr16")
        self.assertEqual(len(alignments.targets[16]), 90354753)
        self.assertEqual(alignments.targets[17].id, "chr17")
        self.assertEqual(len(alignments.targets[17]), 81195210)
        self.assertEqual(alignments.targets[18].id, "chr18")
        self.assertEqual(len(alignments.targets[18]), 78077248)
        self.assertEqual(alignments.targets[19].id, "chr20")
        self.assertEqual(len(alignments.targets[19]), 63025520)
        self.assertEqual(alignments.targets[20].id, "chrY")
        self.assertEqual(len(alignments.targets[20]), 59373566)
        self.assertEqual(alignments.targets[21].id, "chr19")
        self.assertEqual(len(alignments.targets[21]), 59128983)
        self.assertEqual(alignments.targets[22].id, "chr22")
        self.assertEqual(len(alignments.targets[22]), 51304566)
        self.assertEqual(alignments.targets[23].id, "chr21")
        self.assertEqual(len(alignments.targets[23]), 48129895)
        self.assertEqual(alignments.targets[24].id, "chrM")
        self.assertEqual(len(alignments.targets[24]), 16571)
        self.assertEqual(len(alignments.metadata["PG"]), 1)
        self.assertEqual(
            alignments.metadata["PG"][0],
            {
                "ID": "samtools",
                "PN": "samtools",
                "VN": "1.14",
                "CL": "samtools view -h -t hg19.chrom.sizes -",
            },
        )
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg18_dna	0	chr4	61646096	0	11H16M6H	*	0	0	*	*	AS:i:16
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg18_dna	0	chr1	10271784	0	33M	*	0	0	*	*	AS:i:33
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg18_dna	16	chr2	53575981	0	8H17M8H	*	0	0	*	*	AS:i:17
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	0	chr9	85737866	0	9H41M	*	0	0	*	*	AS:i:29
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	0	chr8	95160480	0	8H41M1H	*	0	0	*	*	AS:i:41
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	0	chr22	42144401	0	11H36M3H	*	0	0	*	*	AS:i:24
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	0	chr2	183925985	0	1H6M4I38M1H	*	0	0	*	*	AS:i:27
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	0	chr19	35483341	0	10H25M134D11M4H	*	0	0	*	*	AS:i:0
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	0	chr18	23891311	0	10H39M1H	*	0	0	*	*	AS:i:39
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	0	chr18	43252218	0	21H28M1H	*	0	0	*	*	AS:i:24
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	0	chr13	52759148	0	1H10M3D38M1H	*	0	0	*	*	AS:i:30
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	0	chr1	1207057	0	50M	*	0	0	*	*	AS:i:50
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	0	chr1	61700838	0	1H34M15H	*	0	0	*	*	AS:i:22
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	16	chr4	37558158	0	1H16M4I18M11H	*	0	0	*	*	AS:i:15
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	16	chr22	48997406	0	1H37M12H	*	0	0	*	*	AS:i:29
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	16	chr2	120641741	0	1H36M13H	*	0	0	*	*	AS:i:32
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	16	chr19	54017131	0	1H39M10H	*	0	0	*	*	AS:i:39
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	16	chr19	553743	0	1H39M10H	*	0	0	*	*	AS:i:27
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	16	chr10	99388556	0	1H36M13H	*	0	0	*	*	AS:i:24
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	16	chr10	112178172	0	15H25M10H	*	0	0	*	*	AS:i:21
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	16	chr1	39368491	0	1H36M13H	*	0	0	*	*	AS:i:32
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	16	chr1	220325688	0	3H34M13H	*	0	0	*	*	AS:i:30
""",
        )
        self.assertRaises(StopIteration, next, alignments)

    def test_reading_psl_34_001(self):
        """Test parsing psl_34_001.sam."""
        path = "Blat/psl_34_001.sam"
        alignments = Align.parse(path, "sam")
        self.check_alignments_psl_34_001(alignments)

    def test_writing_psl_34_001(self):
        """Test writing the alignments in psl_34_001.sam."""
        path = "Blat/psl_34_001.sam"
        alignments = Align.parse(path, "sam")
        stream = StringIO()
        n = Align.write(alignments, stream, "sam")
        self.assertEqual(n, 22)
        stream.seek(0)
        alignments = Align.parse(stream, "sam")
        self.check_alignments_psl_34_001(alignments)
        stream.close()

    def check_alignments_psl_34_003(self, alignments):
        """Check the alignments for psl_34_003/sam."""
        self.assertEqual(list(alignments.metadata), ["PG"])
        self.assertEqual(len(alignments.targets), 25)
        self.assertEqual(alignments.targets[0].id, "chr1")
        self.assertEqual(len(alignments.targets[0]), 249250621)
        self.assertEqual(alignments.targets[1].id, "chr2")
        self.assertEqual(len(alignments.targets[1]), 243199373)
        self.assertEqual(alignments.targets[2].id, "chr3")
        self.assertEqual(len(alignments.targets[2]), 198022430)
        self.assertEqual(alignments.targets[3].id, "chr4")
        self.assertEqual(len(alignments.targets[3]), 191154276)
        self.assertEqual(alignments.targets[4].id, "chr5")
        self.assertEqual(len(alignments.targets[4]), 180915260)
        self.assertEqual(alignments.targets[5].id, "chr6")
        self.assertEqual(len(alignments.targets[5]), 171115067)
        self.assertEqual(alignments.targets[6].id, "chr7")
        self.assertEqual(len(alignments.targets[6]), 159138663)
        self.assertEqual(alignments.targets[7].id, "chrX")
        self.assertEqual(len(alignments.targets[7]), 155270560)
        self.assertEqual(alignments.targets[8].id, "chr8")
        self.assertEqual(len(alignments.targets[8]), 146364022)
        self.assertEqual(alignments.targets[9].id, "chr9")
        self.assertEqual(len(alignments.targets[9]), 141213431)
        self.assertEqual(alignments.targets[10].id, "chr10")
        self.assertEqual(len(alignments.targets[10]), 135534747)
        self.assertEqual(alignments.targets[11].id, "chr11")
        self.assertEqual(len(alignments.targets[11]), 135006516)
        self.assertEqual(alignments.targets[12].id, "chr12")
        self.assertEqual(len(alignments.targets[12]), 133851895)
        self.assertEqual(alignments.targets[13].id, "chr13")
        self.assertEqual(len(alignments.targets[13]), 115169878)
        self.assertEqual(alignments.targets[14].id, "chr14")
        self.assertEqual(len(alignments.targets[14]), 107349540)
        self.assertEqual(alignments.targets[15].id, "chr15")
        self.assertEqual(len(alignments.targets[15]), 102531392)
        self.assertEqual(alignments.targets[16].id, "chr16")
        self.assertEqual(len(alignments.targets[16]), 90354753)
        self.assertEqual(alignments.targets[17].id, "chr17")
        self.assertEqual(len(alignments.targets[17]), 81195210)
        self.assertEqual(alignments.targets[18].id, "chr18")
        self.assertEqual(len(alignments.targets[18]), 78077248)
        self.assertEqual(alignments.targets[19].id, "chr20")
        self.assertEqual(len(alignments.targets[19]), 63025520)
        self.assertEqual(alignments.targets[20].id, "chrY")
        self.assertEqual(len(alignments.targets[20]), 59373566)
        self.assertEqual(alignments.targets[21].id, "chr19")
        self.assertEqual(len(alignments.targets[21]), 59128983)
        self.assertEqual(alignments.targets[22].id, "chr22")
        self.assertEqual(len(alignments.targets[22]), 51304566)
        self.assertEqual(alignments.targets[23].id, "chr21")
        self.assertEqual(len(alignments.targets[23]), 48129895)
        self.assertEqual(alignments.targets[24].id, "chrM")
        self.assertEqual(len(alignments.targets[24]), 16571)
        self.assertEqual(len(alignments.targets), 25)
        self.assertEqual(len(alignments.metadata["PG"]), 1)
        self.assertEqual(
            alignments.metadata["PG"][0],
            {
                "ID": "samtools",
                "PN": "samtools",
                "VN": "1.14",
                "CL": "samtools view -h -t hg19.chrom.sizes -",
            },
        )
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg18_dna	0	chr4	61646096	0	11H16M6H	*	0	0	*	*	AS:i:16
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg18_dna	0	chr1	10271784	0	33M	*	0	0	*	*	AS:i:33
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg18_dna	16	chr2	53575981	0	8H17M8H	*	0	0	*	*	AS:i:17
""",
        )
        self.assertRaises(StopIteration, next, alignments)

    def test_reading_psl_34_003(self):
        """Test parsing psl_34_003.sam."""
        path = "Blat/psl_34_003.sam"
        alignments = Align.parse(path, "sam")
        self.check_alignments_psl_34_003(alignments)

    def test_writing_psl_34_003(self):
        """Test writing the alignments in psl_34_003.sam."""
        path = "Blat/psl_34_003.sam"
        alignments = Align.parse(path, "sam")
        stream = StringIO()
        n = Align.write(alignments, stream, "sam")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "sam")
        self.check_alignments_psl_34_003(alignments)
        stream.close()

    def check_alignments_psl_34_004(self, alignments):
        """Check the alignments for psl_34_004/sam."""
        self.assertEqual(list(alignments.metadata), ["PG"])
        self.assertEqual(len(alignments.targets), 25)
        self.assertEqual(alignments.targets[0].id, "chr1")
        self.assertEqual(len(alignments.targets[0]), 249250621)
        self.assertEqual(alignments.targets[1].id, "chr2")
        self.assertEqual(len(alignments.targets[1]), 243199373)
        self.assertEqual(alignments.targets[2].id, "chr3")
        self.assertEqual(len(alignments.targets[2]), 198022430)
        self.assertEqual(alignments.targets[3].id, "chr4")
        self.assertEqual(len(alignments.targets[3]), 191154276)
        self.assertEqual(alignments.targets[4].id, "chr5")
        self.assertEqual(len(alignments.targets[4]), 180915260)
        self.assertEqual(alignments.targets[5].id, "chr6")
        self.assertEqual(len(alignments.targets[5]), 171115067)
        self.assertEqual(alignments.targets[6].id, "chr7")
        self.assertEqual(len(alignments.targets[6]), 159138663)
        self.assertEqual(alignments.targets[7].id, "chrX")
        self.assertEqual(len(alignments.targets[7]), 155270560)
        self.assertEqual(alignments.targets[8].id, "chr8")
        self.assertEqual(len(alignments.targets[8]), 146364022)
        self.assertEqual(alignments.targets[9].id, "chr9")
        self.assertEqual(len(alignments.targets[9]), 141213431)
        self.assertEqual(alignments.targets[10].id, "chr10")
        self.assertEqual(len(alignments.targets[10]), 135534747)
        self.assertEqual(alignments.targets[11].id, "chr11")
        self.assertEqual(len(alignments.targets[11]), 135006516)
        self.assertEqual(alignments.targets[12].id, "chr12")
        self.assertEqual(len(alignments.targets[12]), 133851895)
        self.assertEqual(alignments.targets[13].id, "chr13")
        self.assertEqual(len(alignments.targets[13]), 115169878)
        self.assertEqual(alignments.targets[14].id, "chr14")
        self.assertEqual(len(alignments.targets[14]), 107349540)
        self.assertEqual(alignments.targets[15].id, "chr15")
        self.assertEqual(len(alignments.targets[15]), 102531392)
        self.assertEqual(alignments.targets[16].id, "chr16")
        self.assertEqual(len(alignments.targets[16]), 90354753)
        self.assertEqual(alignments.targets[17].id, "chr17")
        self.assertEqual(len(alignments.targets[17]), 81195210)
        self.assertEqual(alignments.targets[18].id, "chr18")
        self.assertEqual(len(alignments.targets[18]), 78077248)
        self.assertEqual(alignments.targets[19].id, "chr20")
        self.assertEqual(len(alignments.targets[19]), 63025520)
        self.assertEqual(alignments.targets[20].id, "chrY")
        self.assertEqual(len(alignments.targets[20]), 59373566)
        self.assertEqual(alignments.targets[21].id, "chr19")
        self.assertEqual(len(alignments.targets[21]), 59128983)
        self.assertEqual(alignments.targets[22].id, "chr22")
        self.assertEqual(len(alignments.targets[22]), 51304566)
        self.assertEqual(alignments.targets[23].id, "chr21")
        self.assertEqual(len(alignments.targets[23]), 48129895)
        self.assertEqual(alignments.targets[24].id, "chrM")
        self.assertEqual(len(alignments.targets[24]), 16571)
        self.assertEqual(len(alignments.targets), 25)
        self.assertEqual(len(alignments.metadata["PG"]), 1)
        self.assertEqual(
            alignments.metadata["PG"][0],
            {
                "ID": "samtools",
                "PN": "samtools",
                "VN": "1.14",
                "CL": "samtools view -h -t hg19.chrom.sizes -",
            },
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	0	chr9	85737866	0	9H41M	*	0	0	*	*	AS:i:29
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	0	chr8	95160480	0	8H41M1H	*	0	0	*	*	AS:i:41
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	0	chr22	42144401	0	11H36M3H	*	0	0	*	*	AS:i:24
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	0	chr2	183925985	0	1H6M4I38M1H	*	0	0	*	*	AS:i:27
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	0	chr19	35483341	0	10H25M134D11M4H	*	0	0	*	*	AS:i:0
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	0	chr18	23891311	0	10H39M1H	*	0	0	*	*	AS:i:39
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	0	chr18	43252218	0	21H28M1H	*	0	0	*	*	AS:i:24
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	0	chr13	52759148	0	1H10M3D38M1H	*	0	0	*	*	AS:i:30
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	0	chr1	1207057	0	50M	*	0	0	*	*	AS:i:50
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	0	chr1	61700838	0	1H34M15H	*	0	0	*	*	AS:i:22
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	16	chr4	37558158	0	1H16M4I18M11H	*	0	0	*	*	AS:i:15
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	16	chr22	48997406	0	1H37M12H	*	0	0	*	*	AS:i:29
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	16	chr2	120641741	0	1H36M13H	*	0	0	*	*	AS:i:32
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	16	chr19	54017131	0	1H39M10H	*	0	0	*	*	AS:i:39
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	16	chr19	553743	0	1H39M10H	*	0	0	*	*	AS:i:27
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	16	chr10	99388556	0	1H36M13H	*	0	0	*	*	AS:i:24
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	16	chr10	112178172	0	15H25M10H	*	0	0	*	*	AS:i:21
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	16	chr1	39368491	0	1H36M13H	*	0	0	*	*	AS:i:32
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	16	chr1	220325688	0	3H34M13H	*	0	0	*	*	AS:i:30
""",
        )
        self.assertRaises(StopIteration, next, alignments)

    def test_reading_psl_34_004(self):
        """Test parsing psl_34_004.sam."""
        path = "Blat/psl_34_004.sam"
        alignments = Align.parse(path, "sam")
        self.check_alignments_psl_34_004(alignments)

    def test_writing_psl_34_004(self):
        """Test writing the alignments in psl_34_004.sam."""
        path = "Blat/psl_34_004.sam"
        alignments = Align.parse(path, "sam")
        stream = StringIO()
        n = Align.write(alignments, stream, "sam")
        self.assertEqual(n, 19)
        stream.seek(0)
        alignments = Align.parse(stream, "sam")
        self.check_alignments_psl_34_004(alignments)
        stream.close()

    def check_alignments_psl_34_005(self, alignments):
        """Check the alignments for psl_34_005.sam."""
        self.assertEqual(list(alignments.metadata), ["PG"])
        self.assertEqual(len(alignments.targets), 25)
        self.assertEqual(alignments.targets[0].id, "chr1")
        self.assertEqual(len(alignments.targets[0]), 249250621)
        self.assertEqual(alignments.targets[1].id, "chr2")
        self.assertEqual(len(alignments.targets[1]), 243199373)
        self.assertEqual(alignments.targets[2].id, "chr3")
        self.assertEqual(len(alignments.targets[2]), 198022430)
        self.assertEqual(alignments.targets[3].id, "chr4")
        self.assertEqual(len(alignments.targets[3]), 191154276)
        self.assertEqual(alignments.targets[4].id, "chr5")
        self.assertEqual(len(alignments.targets[4]), 180915260)
        self.assertEqual(alignments.targets[5].id, "chr6")
        self.assertEqual(len(alignments.targets[5]), 171115067)
        self.assertEqual(alignments.targets[6].id, "chr7")
        self.assertEqual(len(alignments.targets[6]), 159138663)
        self.assertEqual(alignments.targets[7].id, "chrX")
        self.assertEqual(len(alignments.targets[7]), 155270560)
        self.assertEqual(alignments.targets[8].id, "chr8")
        self.assertEqual(len(alignments.targets[8]), 146364022)
        self.assertEqual(alignments.targets[9].id, "chr9")
        self.assertEqual(len(alignments.targets[9]), 141213431)
        self.assertEqual(alignments.targets[10].id, "chr10")
        self.assertEqual(len(alignments.targets[10]), 135534747)
        self.assertEqual(alignments.targets[11].id, "chr11")
        self.assertEqual(len(alignments.targets[11]), 135006516)
        self.assertEqual(alignments.targets[12].id, "chr12")
        self.assertEqual(len(alignments.targets[12]), 133851895)
        self.assertEqual(alignments.targets[13].id, "chr13")
        self.assertEqual(len(alignments.targets[13]), 115169878)
        self.assertEqual(alignments.targets[14].id, "chr14")
        self.assertEqual(len(alignments.targets[14]), 107349540)
        self.assertEqual(alignments.targets[15].id, "chr15")
        self.assertEqual(len(alignments.targets[15]), 102531392)
        self.assertEqual(alignments.targets[16].id, "chr16")
        self.assertEqual(len(alignments.targets[16]), 90354753)
        self.assertEqual(alignments.targets[17].id, "chr17")
        self.assertEqual(len(alignments.targets[17]), 81195210)
        self.assertEqual(alignments.targets[18].id, "chr18")
        self.assertEqual(len(alignments.targets[18]), 78077248)
        self.assertEqual(alignments.targets[19].id, "chr20")
        self.assertEqual(len(alignments.targets[19]), 63025520)
        self.assertEqual(alignments.targets[20].id, "chrY")
        self.assertEqual(len(alignments.targets[20]), 59373566)
        self.assertEqual(alignments.targets[21].id, "chr19")
        self.assertEqual(len(alignments.targets[21]), 59128983)
        self.assertEqual(alignments.targets[22].id, "chr22")
        self.assertEqual(len(alignments.targets[22]), 51304566)
        self.assertEqual(alignments.targets[23].id, "chr21")
        self.assertEqual(len(alignments.targets[23]), 48129895)
        self.assertEqual(alignments.targets[24].id, "chrM")
        self.assertEqual(len(alignments.targets[24]), 16571)
        self.assertEqual(len(alignments.targets), 25)
        self.assertEqual(len(alignments.metadata["PG"]), 1)
        self.assertEqual(
            alignments.metadata["PG"][0],
            {
                "ID": "samtools",
                "PN": "samtools",
                "VN": "1.14",
                "CL": "samtools view -h -t hg19.chrom.sizes -",
            },
        )
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg18_dna	0	chr4	61646096	0	11S16M6S	*	0	0	*	*	AS:i:16
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg18_dna	0	chr1	10271784	0	33M	*	0	0	*	*	AS:i:33
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg18_dna	16	chr2	53575981	0	8S17M8S	*	0	0	*	*	AS:i:17
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	0	chr9	85737866	0	9S41M	*	0	0	*	*	AS:i:29
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	0	chr8	95160480	0	8S41M1S	*	0	0	*	*	AS:i:41
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	0	chr22	42144401	0	11S36M3S	*	0	0	*	*	AS:i:24
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	0	chr2	183925985	0	1S6M4I38M1S	*	0	0	*	*	AS:i:27
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	0	chr19	35483341	0	10S25M134D11M4S	*	0	0	*	*	AS:i:0
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	0	chr18	23891311	0	10S39M1S	*	0	0	*	*	AS:i:39
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	0	chr18	43252218	0	21S28M1S	*	0	0	*	*	AS:i:24
""",
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
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[52759147, 52759157, 52759160, 52759198],
                             [       1,       11,       11,       49]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	0	chr13	52759148	0	1S10M3D38M1S	*	0	0	*	*	AS:i:30
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	0	chr1	1207057	0	50M	*	0	0	*	*	AS:i:50
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	0	chr1	61700838	0	1S34M15S	*	0	0	*	*	AS:i:22
""",
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
        self.assertEqual(len(alignment.query.seq), 50)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[37558157, 37558173, 37558173, 37558191],
                             [      49,       33,       29,       11]]),
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	16	chr4	37558158	0	1S16M4I18M11S	*	0	0	*	*	AS:i:15
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	16	chr22	48997406	0	1S37M12S	*	0	0	*	*	AS:i:29
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	16	chr2	120641741	0	1S36M13S	*	0	0	*	*	AS:i:32
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	16	chr19	54017131	0	1S39M10S	*	0	0	*	*	AS:i:39
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	16	chr19	553743	0	1S39M10S	*	0	0	*	*	AS:i:27
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	16	chr10	99388556	0	1S36M13S	*	0	0	*	*	AS:i:24
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	16	chr10	112178172	0	15S25M10S	*	0	0	*	*	AS:i:21
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	16	chr1	39368491	0	1S36M13S	*	0	0	*	*	AS:i:32
""",
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
        self.assertEqual(
            format(alignment, "sam"),
            """\
hg19_dna	16	chr1	220325688	0	3S34M13S	*	0	0	*	*	AS:i:30
""",
        )
        self.assertRaises(StopIteration, next, alignments)

    def test_reading_psl_34_005(self):
        """Test parsing psl_34_005.sam."""
        path = "Blat/psl_34_005.sam"
        alignments = Align.parse(path, "sam")
        self.check_alignments_psl_34_005(alignments)

    def test_writing_psl_34_005(self):
        """Test writing the alignments in psl_34_005.sam."""
        path = "Blat/psl_34_005.sam"
        alignments = Align.parse(path, "sam")
        stream = StringIO()
        n = Align.write(alignments, stream, "sam")
        self.assertEqual(n, 22)
        stream.seek(0)
        alignments = Align.parse(stream, "sam")
        self.check_alignments_psl_34_005(alignments)
        stream.close()


class TestAlign_sambam(unittest.TestCase):
    def test_ex1(self):
        alignments = Align.parse("SamBam/ex1.sam", "sam")
        n = 0
        for alignment in alignments:
            n += 1
        self.assertEqual(n, 3270)
        self.assertEqual(alignment.sequences[0].id, "chr2")
        self.assertEqual(alignment.sequences[1].id, "EAS114_26:7:37:79:581")
        self.assertEqual(
            alignment.sequences[1].seq, "TTTTCTGGCATGAAAAAAAAAAAAAAAAAAAAAAA"
        )
        self.assertEqual(alignment.flag, 83)
        self.assertEqual(alignment.mapq, 68)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1532, 1567], [35, 0]])
            )
        )
        self.assertEqual(alignment.rnext, "chr2")
        self.assertEqual(alignment.pnext, 1348)
        self.assertEqual(alignment.tlen, -219)
        self.assertEqual(
            alignment.sequences[1].letter_annotations["phred_quality"],
            "3,,,===6===<===<;=====-============",
        )
        self.assertEqual(len(alignment.annotations), 6)
        self.assertEqual(alignment.annotations["MF"], 18)
        self.assertEqual(alignment.annotations["Aq"], 27)
        self.assertEqual(alignment.annotations["NM"], 2)
        self.assertEqual(alignment.annotations["UQ"], 23)
        self.assertEqual(alignment.annotations["H0"], 0)
        self.assertEqual(alignment.annotations["H1"], 1)

    def test_ex1_header(self):
        alignments = Align.parse("SamBam/ex1_header.sam", "sam")
        self.assertEqual(alignments.metadata["HD"], {"VN": "1.3", "SO": "coordinate"})
        self.assertEqual(len(alignments.targets), 2)
        self.assertEqual(alignments.targets[0].id, "chr1")
        self.assertEqual(len(alignments.targets[0].seq), 1575)
        self.assertEqual(alignments.targets[1].id, "chr2")
        self.assertEqual(len(alignments.targets[1].seq), 1584)
        n = 0
        for alignment in alignments:
            n += 1
        self.assertEqual(n, 3270)
        self.assertEqual(alignment.sequences[0].id, "chr2")
        self.assertEqual(len(alignment.sequences[0].seq), 1584)
        self.assertEqual(alignment.sequences[1].id, "EAS114_26:7:37:79:581")
        self.assertEqual(
            alignment.sequences[1].seq, "TTTTCTGGCATGAAAAAAAAAAAAAAAAAAAAAAA"
        )
        self.assertEqual(alignment.flag, 83)
        self.assertEqual(alignment.mapq, 68)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1532, 1567], [35, 0]])
            )
        )
        self.assertEqual(alignment.rnext, "chr2")
        self.assertEqual(alignment.pnext, 1348)
        self.assertEqual(alignment.tlen, -219)
        self.assertEqual(
            alignment.sequences[1].letter_annotations["phred_quality"],
            "3,,,===6===<===<;=====-============",
        )
        self.assertEqual(len(alignment.annotations), 6)
        self.assertEqual(alignment.annotations["MF"], 18)
        self.assertEqual(alignment.annotations["Aq"], 27)
        self.assertEqual(alignment.annotations["NM"], 2)
        self.assertEqual(alignment.annotations["UQ"], 23)
        self.assertEqual(alignment.annotations["H0"], 0)
        self.assertEqual(alignment.annotations["H1"], 1)

    def test_sam1(self):
        alignments = Align.parse("SamBam/sam1.sam", "sam")
        self.assertEqual(len(alignments.targets), 1)
        self.assertEqual(alignments.targets[0].id, "1")
        self.assertEqual(len(alignments.targets[0].seq), 239940)
        self.assertEqual(
            alignments.metadata["PG"][0],
            {
                "ID": "bwa",
                "PN": "bwa",
                "VN": "0.6.2-r126",
            },
        )
        n = 0
        for alignment in alignments:
            n += 1
        self.assertEqual(n, 200)
        self.assertIsNone(alignment.sequences[0])
        self.assertEqual(
            alignment.sequences[1].id, "HWI-1KL120:88:D0LRBACXX:1:1101:5516:2195"
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "GGCCCAACCGTCCTATATGAGATGTAGCATGGTACAGAACAAACTGCTTACACAGGTCTCACTAGTTAGAAACCTGTGGGCCATGGAGGTCAGACATCCAT",
        )
        self.assertEqual(alignment.flag, 141)
        self.assertEqual(alignment.mapq, 0)
        self.assertIsNone(alignment.coordinates)
        self.assertEqual(
            alignment.sequences[1].letter_annotations["phred_quality"],
            "B?1ADDDDAFFFDEGFEGHEED?D?EB<EGB;F>FHI>GEBHEF@@<BF>D?F<FB=C>F;C@FC7@=;=E=7=?@;;;856?@;;;;559(,,5?3>5>@",
        )

    def test_sam2(self):
        alignments = Align.parse("SamBam/sam2.sam", "sam")
        self.assertEqual(len(alignments.targets), 1)
        self.assertEqual(alignments.targets[0].id, "1")
        self.assertEqual(len(alignments.targets[0].seq), 239940)
        self.assertEqual(
            alignments.metadata["PG"][0],
            {
                "ID": "bwa",
                "PN": "bwa",
                "VN": "0.6.2-r126",
            },
        )
        n = 0
        for alignment in alignments:
            if n == 8:
                self.assertEqual(alignment.sequences[0].id, "1")
                self.assertEqual(len(alignment.sequences[0].seq), 239940)
                self.assertEqual(
                    alignment.sequences[0].seq.defined_ranges, ((132615, 132716),)
                )
                self.assertEqual(
                    alignment.sequences[0].seq[132615:132716],
                    "GGTCACACCCTGTCCTCCTCCTACACATACTCGGATGCTTCCTCCTCAACCTTGGCACCCACCTCCTTCTTACTGGGCCCAGGAGCCTTCAAAGCCCAGGA",
                )
                self.assertEqual(
                    alignment.sequences[1].id,
                    "HWI-1KL120:88:D0LRBACXX:1:1101:2205:2204",
                )
                self.assertEqual(
                    alignment.sequences[1].seq,
                    "TCCTGGGCATTGAAGGCTCCTGGGCCCAGTAAGAAGGAGGTGGGTGCCAAGGTTGAGGAGGAAGCATCCGAGTATGTGTAGGAGGAGGACAAGGTGGGACC",
                )
                self.assertEqual(alignment.flag, 83)
                self.assertEqual(alignment.mapq, 60)
                self.assertTrue(
                    numpy.array_equal(
                        alignment.coordinates, numpy.array([[132615, 132716], [101, 0]])
                    )
                )
                self.assertEqual(alignment.rnext, "1")
                self.assertEqual(alignment.pnext, 132490)
                self.assertEqual(alignment.tlen, -226)
                self.assertEqual(
                    alignment.sequences[1].letter_annotations["phred_quality"],
                    "BBB@?C>???CBBDDDDDDDCC>C>>C???=DEEEDFEDBGGHIIEED=HFAGIIHDHGD?GIJJJIHIGFDHFIJJJIJJJJJJJJJHHHHFFFFFFCCC",
                )
                self.assertEqual(len(alignment.annotations), 9)
                self.assertEqual(alignment.annotations["XT"], "U")
                self.assertEqual(alignment.annotations["NM"], 3)
                self.assertEqual(alignment.annotations["SM"], 37)
                self.assertEqual(alignment.annotations["AM"], 37)
                self.assertEqual(alignment.annotations["X0"], 1)
                self.assertEqual(alignment.annotations["X1"], 0)
                self.assertEqual(alignment.annotations["XM"], 3)
                self.assertEqual(alignment.annotations["XO"], 0)
                self.assertEqual(alignment.annotations["XG"], 0)
            elif n == 9:
                self.assertEqual(alignment.sequences[0].id, "1")
                self.assertEqual(len(alignment.sequences[0].seq), 239940)
                self.assertEqual(
                    alignment.sequences[0].seq.defined_ranges, ((132490, 132591),)
                )
                self.assertEqual(
                    alignment.sequences[0].seq[132490:132591],
                    "GCAACAAGGGCTTTGGTGGGAAGGTATTTGCACCTGTCATTCCTTCCTCCTTTACTCCTGCCGCCCCTTGCTGGATCCTGAGCCCCCAGGGTCCCCCGATC",
                )
                self.assertEqual(
                    alignment.sequences[1].id,
                    "HWI-1KL120:88:D0LRBACXX:1:1101:2205:2204",
                )
                self.assertEqual(
                    alignment.sequences[1].seq,
                    "GCAACAAGGGCTTTGGTGGGAAGGTATCTGCACCTGTCATTCCTTCCTCCTTTACTCCTGCCGCCCCTTGCTGGATCCTGAGCCCCCAGGGTCCCCCGATC",
                )
                self.assertEqual(alignment.flag, 163)
                self.assertEqual(alignment.mapq, 60)
                self.assertTrue(
                    numpy.array_equal(
                        alignment.coordinates, numpy.array([[132490, 132591], [0, 101]])
                    )
                )
                self.assertEqual(alignment.rnext, "1")
                self.assertEqual(alignment.pnext, 132615)
                self.assertEqual(alignment.tlen, 226)
                self.assertEqual(
                    alignment.sequences[1].letter_annotations["phred_quality"],
                    "CCCDFFFFHHHHHJJHEHIJIIIJ?EHIJIIJJJGHFGCHGJJIIIJJJJJJHIIIIIJJJJIJHFFBECEEDDBDDDDDDDDDDDD>@<59ABDDBB###",
                )
                self.assertEqual(len(alignment.annotations), 9)
                self.assertEqual(alignment.annotations["XT"], "U")
                self.assertEqual(alignment.annotations["NM"], 1)
                self.assertEqual(alignment.annotations["SM"], 37)
                self.assertEqual(alignment.annotations["AM"], 37)
                self.assertEqual(alignment.annotations["X0"], 1)
                self.assertEqual(alignment.annotations["X1"], 0)
                self.assertEqual(alignment.annotations["XM"], 1)
                self.assertEqual(alignment.annotations["XO"], 0)
                self.assertEqual(alignment.annotations["XG"], 0)
            elif n == 100:
                self.assertEqual(alignment.sequences[0].id, "1")
                self.assertEqual(len(alignment.sequences[0].seq), 239940)
                self.assertEqual(
                    alignment.sequences[0].seq.defined_ranges, ((137538, 137639),)
                )
                self.assertEqual(
                    alignment.sequences[0].seq[137538:137639],
                    "AAAGTTCGGGGCCTACAAAGGCGGTTGGGAGCTGGGCAGGAGTTGAGCCAAAAGAGCTTGCTTACTTGCTGGGAGGCAGGGCCGGGAGAGCCCGACTTCAG",
                )
                self.assertEqual(
                    alignment.sequences[1].id,
                    "HWI-1KL120:88:D0LRBACXX:1:1101:4673:2125",
                )
                self.assertEqual(
                    alignment.sequences[1].seq,
                    "AAAGTTCGGGGCCTACAAAGGCGGTTGGGAGCTGGGCAGGAGTTGAGCCAAAAGAGCTTGCTTACTTGCTGGGAGGCAGGACCGGGAGAGGCCGACTTCAG",
                )
                self.assertEqual(alignment.flag, 97)
                self.assertEqual(alignment.mapq, 37)
                self.assertTrue(
                    numpy.array_equal(
                        alignment.coordinates, numpy.array([[137538, 137639], [0, 101]])
                    )
                )
                self.assertEqual(alignment.rnext, "1")
                self.assertEqual(alignment.pnext, 135649)
                self.assertEqual(alignment.tlen, -1788)
                self.assertEqual(
                    alignment.sequences[1].letter_annotations["phred_quality"],
                    "CCCFFFFFHHHHHJJJJJJJJJJJGHGIJIIIJJIJIHHHFFCCCEECEEDDBDDDCDDDCDCDDDDDDDDDD?B@BDDDDDDDDDDDBDDDB>@DB@CCD",
                )
                self.assertEqual(len(alignment.annotations), 9)
                self.assertEqual(alignment.annotations["XT"], "U")
                self.assertEqual(alignment.annotations["NM"], 2)
                self.assertEqual(alignment.annotations["SM"], 37)
                self.assertEqual(alignment.annotations["AM"], 37)
                self.assertEqual(alignment.annotations["X0"], 1)
                self.assertEqual(alignment.annotations["X1"], 0)
                self.assertEqual(alignment.annotations["XM"], 2)
                self.assertEqual(alignment.annotations["XO"], 0)
                self.assertEqual(alignment.annotations["XG"], 0)
            elif n == 101:
                self.assertEqual(alignment.sequences[0].id, "1")
                self.assertEqual(len(alignment.sequences[0].seq), 239940)
                self.assertEqual(
                    alignment.sequences[0].seq.defined_ranges, ((135649, 135750),)
                )
                self.assertEqual(
                    alignment.sequences[0].seq[135649:135750],
                    "TGGAGAGGCCACCGCGAGGCCTGAGCTGGGCCTGGGGAGCTTGGCTTAGGGAAGTTGTGGGCCTACCAGGGCCGCTGGGAGCTGGGCAGGAGCTGAGTCCA",
                )
                self.assertEqual(
                    alignment.sequences[1].id,
                    "HWI-1KL120:88:D0LRBACXX:1:1101:4673:2125",
                )
                self.assertEqual(
                    alignment.sequences[1].seq,
                    "TGGACTCAGCTCCTGCCCAGCTCCCAGCGGCCCTGGTAGGCCCACAACTTCCCGAAGCCAAGCTCCCCAGGCCCAGCTCAGGCCTCACGGTGGCCTCTCCA",
                )
                self.assertEqual(alignment.flag, 145)
                self.assertEqual(alignment.mapq, 37)
                self.assertTrue(
                    numpy.array_equal(
                        alignment.coordinates, numpy.array([[135649, 135750], [101, 0]])
                    )
                )
                self.assertEqual(alignment.rnext, "1")
                self.assertEqual(alignment.pnext, 137538)
                self.assertEqual(alignment.tlen, 1788)
                self.assertEqual(
                    alignment.sequences[1].letter_annotations["phred_quality"],
                    "CCCABCAABB@BBDDDDBDCDCDDBDDDDDB?DDDDDCBDECEDFFHFHIIJJIJJJIJJIIHHGJIJIJIJIGJJJJJJIJIIJJJJHHHHHFFFFFCCC",
                )
                self.assertEqual(len(alignment.annotations), 9)
                self.assertEqual(alignment.annotations["XT"], "U")
                self.assertEqual(alignment.annotations["NM"], 2)
                self.assertEqual(alignment.annotations["SM"], 37)
                self.assertEqual(alignment.annotations["AM"], 37)
                self.assertEqual(alignment.annotations["X0"], 1)
                self.assertEqual(alignment.annotations["X1"], 0)
                self.assertEqual(alignment.annotations["XM"], 2)
                self.assertEqual(alignment.annotations["XO"], 0)
                self.assertEqual(alignment.annotations["XG"], 0)
            else:
                self.assertIsNone(alignment.sequences[0])
                self.assertEqual(alignment.mapq, 0)
                self.assertIsNone(alignment.coordinates)
            n += 1
        self.assertEqual(n, 200)


class TestAlign_clippping(unittest.TestCase):
    def test_6M(self):
        """Test alignment starting at non-zero position."""
        target_seq = Seq("AAAAAAAACCCCCC")
        query_seq = Seq("CCCCCC")
        target = SeqRecord(target_seq, id="target")
        query = SeqRecord(query_seq, id="query")
        sequences = [target, query]
        coordinates = numpy.array([[8, 14], [0, 6]])
        alignment = Alignment(sequences, coordinates)
        self.assertEqual(
            str(alignment),
            """\
target            8 CCCCCC 14
                  0 ||||||  6
query             0 CCCCCC  6
""",
        )
        line = alignment.format("sam")
        self.assertEqual(line, "query\t0\ttarget\t9\t255\t6M\t*\t0\t0\tCCCCCC\t*\n")
        fields = line.split()
        pos = int(fields[3]) - 1
        self.assertEqual(pos, 8)
        cigar = fields[5]
        self.assertEqual(cigar, "6M")
        stream = StringIO(line)
        alignments = Align.parse(stream, "sam")
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['C', 'C', 'C', 'C', 'C', 'C'],
             ['C', 'C', 'C', 'C', 'C', 'C']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        stream.close()
        self.assertTrue(numpy.array_equal(alignment.coordinates, coordinates))

    def test_8D6M_ex1(self):
        """Test alignment starting with deletion."""
        target_seq = Seq("AAAAAAAACCCCCC")
        query_seq = Seq("CCCCCC")
        target = SeqRecord(target_seq, id="target")
        query = SeqRecord(query_seq, id="query")
        sequences = [target, query]
        coordinates = numpy.array([[0, 8, 14], [0, 0, 6]])
        alignment = Alignment(sequences, coordinates)
        self.assertEqual(
            str(alignment),
            """\
target            0 AAAAAAAACCCCCC 14
                  0 --------|||||| 14
query             0 --------CCCCCC  6
""",
        )
        line = alignment.format("sam")
        self.assertEqual(line, "query\t0\ttarget\t1\t255\t8D6M\t*\t0\t0\tCCCCCC\t*\n")
        fields = line.split()
        pos = int(fields[3]) - 1
        self.assertEqual(pos, 0)
        cigar = fields[5]
        self.assertEqual(cigar, "8D6M")
        stream = StringIO(line)
        alignments = Align.parse(stream, "sam")
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'C', 'C', 'C', 'C', 'C',
              'C'],
             ['-', '-', '-', '-', '-', '-', '-', '-', 'C', 'C', 'C', 'C', 'C',
              'C']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        stream.close()
        self.assertTrue(numpy.array_equal(alignment.coordinates, coordinates))

    def test_8D6M_ex2(self):
        """Test alignment starting with deletion at non-zero position."""
        target_seq = Seq("GGGGAAAAAAAACCCCCC")
        query_seq = Seq("CCCCCC")
        target = SeqRecord(target_seq, id="target")
        query = SeqRecord(query_seq, id="query")
        sequences = [target, query]
        coordinates = numpy.array([[4, 12, 18], [0, 0, 6]])
        alignment = Alignment(sequences, coordinates)
        self.assertEqual(
            str(alignment),
            """\
target            4 AAAAAAAACCCCCC 18
                  0 --------|||||| 14
query             0 --------CCCCCC  6
""",
        )
        line = alignment.format("sam")
        self.assertEqual(line, "query\t0\ttarget\t5\t255\t8D6M\t*\t0\t0\tCCCCCC\t*\n")
        fields = line.split()
        pos = int(fields[3]) - 1
        self.assertEqual(pos, 4)
        cigar = fields[5]
        self.assertEqual(cigar, "8D6M")
        stream = StringIO(line)
        alignments = Align.parse(stream, "sam")
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'C', 'C', 'C', 'C', 'C',
              'C'],
             ['-', '-', '-', '-', '-', '-', '-', '-', 'C', 'C', 'C', 'C', 'C',
              'C']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        stream.close()
        self.assertTrue(numpy.array_equal(alignment.coordinates, coordinates))

    def test_8I6M_ex1(self):
        """Test alignment starting with insertion."""
        target_seq = Seq("CCCCCC")
        query_seq = Seq("AAAAAAAACCCCCC")
        target = SeqRecord(target_seq, id="target")
        query = SeqRecord(query_seq, id="query")
        sequences = [target, query]
        coordinates = numpy.array([[0, 0, 6], [0, 8, 14]])
        alignment = Alignment(sequences, coordinates)
        self.assertEqual(
            str(alignment),
            """\
target            0 --------CCCCCC  6
                  0 --------|||||| 14
query             0 AAAAAAAACCCCCC 14
""",
        )
        line = alignment.format("sam")
        self.assertEqual(
            line, "query\t0\ttarget\t1\t255\t8I6M\t*\t0\t0\tAAAAAAAACCCCCC\t*\n"
        )
        fields = line.split()
        pos = int(fields[3]) - 1
        self.assertEqual(pos, 0)
        cigar = fields[5]
        self.assertEqual(cigar, "8I6M")
        stream = StringIO(line)
        alignments = Align.parse(stream, "sam")
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['-', '-', '-', '-', '-', '-', '-', '-', 'C', 'C', 'C', 'C', 'C',
              'C'],
             ['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'C', 'C', 'C', 'C', 'C',
              'C']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        stream.close()
        self.assertTrue(numpy.array_equal(alignment.coordinates, coordinates))

    def test_8I6M_ex2(self):
        """Test alignment starting with insertion at non-zero position."""
        target_seq = Seq("GGGGCCCCCC")
        query_seq = Seq("AAAAAAAACCCCCC")
        target = SeqRecord(target_seq, id="target")
        query = SeqRecord(query_seq, id="query")
        sequences = [target, query]
        coordinates = numpy.array([[4, 4, 10], [0, 8, 14]])
        alignment = Alignment(sequences, coordinates)
        self.assertEqual(
            str(alignment),
            """\
target            4 --------CCCCCC 10
                  0 --------|||||| 14
query             0 AAAAAAAACCCCCC 14
""",
        )
        line = alignment.format("sam")
        self.assertEqual(
            line, "query\t0\ttarget\t5\t255\t8I6M\t*\t0\t0\tAAAAAAAACCCCCC\t*\n"
        )
        fields = line.split()
        pos = int(fields[3]) - 1
        self.assertEqual(pos, 4)
        cigar = fields[5]
        self.assertEqual(cigar, "8I6M")
        stream = StringIO(line)
        alignments = Align.parse(stream, "sam")
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['-', '-', '-', '-', '-', '-', '-', '-', 'C', 'C', 'C', 'C', 'C',
              'C'],
             ['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'C', 'C', 'C', 'C', 'C',
              'C']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        stream.close()
        self.assertTrue(numpy.array_equal(alignment.coordinates, coordinates))

    def test_8S6M(self):
        """Test alignment starting with soft clip."""
        target_seq = Seq("CCCCCC")
        query_seq = Seq("AAAAAAAACCCCCC")
        target = SeqRecord(target_seq, id="target")
        query = SeqRecord(query_seq, id="query")
        sequences = [target, query]
        coordinates = numpy.array([[0, 6], [8, 14]])
        alignment = Alignment(sequences, coordinates)
        self.assertEqual(
            str(alignment),
            """\
target            0 CCCCCC  6
                  0 ||||||  6
query             8 CCCCCC 14
""",
        )
        line = alignment.format("sam")
        self.assertEqual(
            line, "query\t0\ttarget\t1\t255\t8S6M\t*\t0\t0\tAAAAAAAACCCCCC\t*\n"
        )
        fields = line.split()
        pos = int(fields[3]) - 1
        self.assertEqual(pos, 0)
        cigar = fields[5]
        self.assertEqual(cigar, "8S6M")
        stream = StringIO(line)
        alignments = Align.parse(stream, "sam")
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['C', 'C', 'C', 'C', 'C', 'C'],
             ['C', 'C', 'C', 'C', 'C', 'C']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        stream.close()
        self.assertTrue(numpy.array_equal(alignment.coordinates, coordinates))

    def test_4S8D6M(self):
        """Test alignment starting with soft clip followed by deletion."""
        target_seq = Seq("AAAAAAAACCCCCC")
        query_seq = Seq("GGGGCCCCCC")
        target = SeqRecord(target_seq, id="target")
        query = SeqRecord(query_seq, id="query")
        sequences = [target, query]
        coordinates = numpy.array([[0, 8, 14], [4, 4, 10]])
        alignment = Alignment(sequences, coordinates)
        self.assertEqual(
            str(alignment),
            """\
target            0 AAAAAAAACCCCCC 14
                  0 --------|||||| 14
query             4 --------CCCCCC 10
""",
        )
        line = alignment.format("sam")
        self.assertEqual(
            line, "query\t0\ttarget\t1\t255\t4S8D6M\t*\t0\t0\tGGGGCCCCCC\t*\n"
        )
        fields = line.split()
        pos = int(fields[3]) - 1
        self.assertEqual(pos, 0)
        cigar = fields[5]
        self.assertEqual(cigar, "4S8D6M")
        stream = StringIO(line)
        alignments = Align.parse(stream, "sam")
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'C', 'C', 'C', 'C', 'C',
              'C'],
             ['-', '-', '-', '-', '-', '-', '-', '-', 'C', 'C', 'C', 'C', 'C',
              'C']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        stream.close()
        self.assertTrue(numpy.array_equal(alignment.coordinates, coordinates))

    def test_4I8D6M(self):
        """Test alignment starting with insertion followed by deletion."""
        target_seq = Seq("AAAAAAAACCCCCC")
        query_seq = Seq("GGGGCCCCCC")
        target = SeqRecord(target_seq, id="target")
        query = SeqRecord(query_seq, id="query")
        sequences = [target, query]
        coordinates = numpy.array([[0, 0, 8, 14], [0, 4, 4, 10]])
        alignment = Alignment(sequences, coordinates)
        self.assertEqual(
            str(alignment),
            """\
target            0 ----AAAAAAAACCCCCC 14
                  0 ------------|||||| 18
query             0 GGGG--------CCCCCC 10
""",
        )
        line = alignment.format("sam")
        self.assertEqual(
            line, "query\t0\ttarget\t1\t255\t4I8D6M\t*\t0\t0\tGGGGCCCCCC\t*\n"
        )
        fields = line.split()
        pos = int(fields[3]) - 1
        self.assertEqual(pos, 0)
        cigar = fields[5]
        self.assertEqual(cigar, "4I8D6M")
        stream = StringIO(line)
        alignments = Align.parse(stream, "sam")
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['-', '-', '-', '-', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'C',
              'C', 'C', 'C', 'C', 'C'],
             ['G', 'G', 'G', 'G', '-', '-', '-', '-', '-', '-', '-', '-', 'C',
              'C', 'C', 'C', 'C', 'C']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        stream.close()
        self.assertTrue(numpy.array_equal(alignment.coordinates, coordinates))

    def test_4S6M(self):
        """Test alignment starting with soft clip at non-zero position."""
        target_seq = Seq("AAAAAAAACCCCCC")
        query_seq = Seq("GGGGCCCCCC")
        target = SeqRecord(target_seq, id="target")
        query = SeqRecord(query_seq, id="query")
        sequences = [target, query]
        coordinates = numpy.array([[8, 14], [4, 10]])
        alignment = Alignment(sequences, coordinates)
        self.assertEqual(
            str(alignment),
            """\
target            8 CCCCCC 14
                  0 ||||||  6
query             4 CCCCCC 10
""",
        )
        line = alignment.format("sam")
        self.assertEqual(
            line, "query\t0\ttarget\t9\t255\t4S6M\t*\t0\t0\tGGGGCCCCCC\t*\n"
        )
        fields = line.split()
        pos = int(fields[3]) - 1
        self.assertEqual(pos, 8)
        cigar = fields[5]
        self.assertEqual(cigar, "4S6M")
        stream = StringIO(line)
        alignments = Align.parse(stream, "sam")
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['C', 'C', 'C', 'C', 'C', 'C'],
             ['C', 'C', 'C', 'C', 'C', 'C']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        stream.close()
        self.assertTrue(numpy.array_equal(alignment.coordinates, coordinates))

    def test_4D8I6M(self):
        """Test alignment starting with deletion followed by insertion."""
        target_seq = Seq("GGGGCCCCCC")
        query_seq = Seq("AAAAAAAACCCCCC")
        target = SeqRecord(target_seq, id="target")
        query = SeqRecord(query_seq, id="query")
        sequences = [target, query]
        coordinates = numpy.array([[0, 4, 4, 10], [0, 0, 8, 14]])
        alignment = Alignment(sequences, coordinates)
        self.assertEqual(
            str(alignment),
            """\
target            0 GGGG--------CCCCCC 10
                  0 ------------|||||| 18
query             0 ----AAAAAAAACCCCCC 14
""",
        )
        line = alignment.format("sam")
        self.assertEqual(
            line, "query\t0\ttarget\t1\t255\t4D8I6M\t*\t0\t0\tAAAAAAAACCCCCC\t*\n"
        )
        fields = line.split()
        pos = int(fields[3]) - 1
        self.assertEqual(pos, 0)
        cigar = fields[5]
        self.assertEqual(cigar, "4D8I6M")
        stream = StringIO(line)
        alignments = Align.parse(stream, "sam")
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['G', 'G', 'G', 'G', '-', '-', '-', '-', '-', '-', '-', '-', 'C',
              'C', 'C', 'C', 'C', 'C'],
             ['-', '-', '-', '-', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'C',
              'C', 'C', 'C', 'C', 'C']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        stream.close()
        self.assertTrue(numpy.array_equal(alignment.coordinates, coordinates))

    def test_4S8I6M(self):
        """Test alignment starting with soft clip followed by insertion."""
        target_seq = Seq("CCCCCC")
        query_seq = Seq("GGGGAAAAAAAACCCCCC")
        target = SeqRecord(target_seq, id="target")
        query = SeqRecord(query_seq, id="query")
        sequences = [target, query]
        coordinates = numpy.array([[0, 0, 6], [4, 12, 18]])
        alignment = Alignment(sequences, coordinates)
        self.assertEqual(
            str(alignment),
            """\
target            0 --------CCCCCC  6
                  0 --------|||||| 14
query             4 AAAAAAAACCCCCC 18
""",
        )
        line = alignment.format("sam")
        self.assertEqual(
            line, "query\t0\ttarget\t1\t255\t4S8I6M\t*\t0\t0\tGGGGAAAAAAAACCCCCC\t*\n"
        )
        fields = line.split()
        pos = int(fields[3]) - 1
        self.assertEqual(pos, 0)
        cigar = fields[5]
        self.assertEqual(cigar, "4S8I6M")
        stream = StringIO(line)
        alignments = Align.parse(stream, "sam")
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['-', '-', '-', '-', '-', '-', '-', '-', 'C', 'C', 'C', 'C', 'C',
              'C'],
             ['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'C', 'C', 'C', 'C', 'C',
              'C']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        stream.close()
        self.assertTrue(numpy.array_equal(alignment.coordinates, coordinates))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
