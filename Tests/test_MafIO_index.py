# Copyright 2012 by Andrew Sczesnak.  All rights reserved.
# Revisions Copyright 2017 by Blaise Li.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Unit tests for Bio.AlignIO.MafIO.MafIndex()."""

try:
    import sqlite3
except ImportError:
    # skip most tests if sqlite is not available
    sqlite3 = None

import os
import unittest
import tempfile
import shutil
import sys

from Bio.AlignIO.MafIO import MafIndex
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from seq_tests_common import SeqRecordTestBaseClass


class StaticMethodTest(unittest.TestCase):
    """Test static UCSC binning-related functions."""

    def test_region2bin(self):
        data = [
            (25079603, 25079787, {0, 1, 11, 96, 776}),
            (25128173, 25128248, {0, 1, 11, 96, 776}),
            (50312474, 50312703, {0, 1, 968, 14, 120}),
            (41905591, 41906101, {0, 1, 904, 13, 112}),
            (16670899, 16673060, {0, 1, 10, 712, 88}),
            (75495356, 75495494, {0, 1, 2, 1160, 144, 17}),
            (92259501, 92261053, {0, 1, 2, 1288, 160, 19}),
            (83834063, 83838132, {0, 1, 2, 1224, 18, 152}),
            (7309597, 7310411, {0, 1, 640, 79, 9}),
            (6190410, 6190999, {0, 1, 632, 78, 9}),
        ]

        for x, y, z in data:
            self.assertEqual(MafIndex._region2bin(x, y), z)

        for x, y, z in data:
            self.assertRaises(TypeError, MafIndex._region2bin, str(x), str(y))

    def test_ucscbin(self):
        data = [
            (25079603, 25079787, 776),
            (25128173, 25128248, 776),
            (50312474, 50312703, 968),
            (41905591, 41906101, 904),
            (16670899, 16673060, 712),
            (75495356, 75495494, 1160),
            (92259501, 92261053, 1288),
            (83834063, 83838132, 1224),
            (7309597, 7310411, 640),
            (6190410, 6190999, 632),
        ]

        for x, y, z in data:
            self.assertEqual(MafIndex._ucscbin(x, y), z)

        for x, y, z in data:
            self.assertRaises(TypeError, MafIndex._ucscbin, str(x), str(y))


if sqlite3:

    class PreBuiltIndexTest(unittest.TestCase):
        """Test loading of prebuilt indices."""

        def test_old(self):
            idx = MafIndex(
                "MAF/ucsc_mm9_chr10.mafindex", "MAF/ucsc_mm9_chr10.maf", "mm9.chr10"
            )
            self.assertEqual(len(idx), 48)

        def test_old_wrong_target_seqname(self):
            self.assertRaises(
                ValueError,
                MafIndex,
                "MAF/ucsc_mm9_chr10.mafindex",
                "MAF/ucsc_mm9_chr10.maf",
                "mm9.chr11",
            )

        def test_old_wrong_filename(self):
            self.assertRaises(
                ValueError,
                MafIndex,
                "MAF/ucsc_mm9_chr10.mafindex",
                "MAF/humor.maf",
                "mm9.chr10",
            )

        def test_old_file_not_found(self):
            self.assertRaises(
                FileNotFoundError,
                MafIndex,
                "MAF/ucsc_mm9_chr11.mafindex",
                "MAF/ucsc_mm9_chr11.maf",
                "mm9.chr11",
            )

        def test_old_wrong_version(self):
            self.assertRaises(
                ValueError,
                MafIndex,
                "MAF/wrong_version.idx",
                "MAF/ucsc_mm9_chr10.maf",
                "mm9.chr10",
            )

        def test_old_unfinished_index(self):
            self.assertRaises(
                ValueError,
                MafIndex,
                "MAF/unfinished.idx",
                "MAF/ucsc_mm9_chr10.maf",
                "mm9.chr10",
            )

        def test_old_corrupt_index(self):
            self.assertRaises(
                ValueError,
                MafIndex,
                "MAF/corrupt.idx",
                "MAF/ucsc_mm9_chr10.maf",
                "mm9.chr10",
            )

        def test_old_invalid_sqlite(self):
            self.assertRaises(
                ValueError,
                MafIndex,
                "MAF/invalid.idx",
                "MAF/ucsc_mm9_chr10.maf",
                "mm9.chr10",
            )

    class NewIndexTest(unittest.TestCase):
        """Test creation of new indices."""

        def setUp(self):
            self.tmpdir = tempfile.mkdtemp()
            self.tmpfile = self.tmpdir + "/database.sqlite3"

        def tearDown(self):
            if os.path.isdir(self.tmpdir):
                shutil.rmtree(self.tmpdir)

        def test_good_small(self):
            idx = MafIndex(self.tmpfile, "MAF/ucsc_mm9_chr10.maf", "mm9.chr10")
            self.assertEqual(len(idx), 48)

        def test_good_big(self):
            idx = MafIndex(self.tmpfile, "MAF/ucsc_mm9_chr10_big.maf", "mm9.chr10")
            self.assertEqual(len(idx), 983)

        def test_bundle_without_target(self):
            self.assertRaises(
                ValueError,
                MafIndex,
                self.tmpfile,
                "MAF/bundle_without_target.maf",
                "mm9.chr10",
            )

        def test_length_coords_mismatch(self):
            self.assertRaises(
                ValueError,
                MafIndex,
                self.tmpfile,
                "MAF/length_coords_mismatch.maf",
                "mm9.chr10",
            )

    class TestGetRecord(SeqRecordTestBaseClass):
        """Make sure we can seek and fetch records properly."""

        def setUp(self):
            self.idx = MafIndex(
                "MAF/ucsc_mm9_chr10.mafindex", "MAF/ucsc_mm9_chr10.maf", "mm9.chr10"
            )
            self.assertEqual(len(self.idx), 48)

        def test_records_begin(self):
            rec1 = SeqRecord(
                Seq(
                    "TCATAGGTATTTATTTTTAAATATGGTTTGCTTTATGGCTAGAA"
                    "CACACCGATTACTTAAAATAGGATTAACC--CCCATACACTTTA"
                    "AAAATGATTAAACAACATTTCTGCTGCTCGCTCACATTCTTCAT"
                    "AGAAGATGACATAATGTATTTTCCTTTTGGTT"
                ),
                id="mm9.chr10",
                name="mm9.chr10",
                description="",
                annotations={
                    "start": 3009319,
                    "srcSize": 129993255,
                    "strand": 1,
                    "size": 162,
                },
            )

            rec2 = SeqRecord(
                Seq(
                    "TCACAGATATTTACTATTAAATATGGTTTGTTATATGGTTACGG"
                    "TTCATAGGTTACTTGGAATTGGATTAACCTTCTTATTCATTGCA"
                    "GAATTGGTTACACTGTGTTCTTGACCTTTGCTTGTTTTCTCCAT"
                    "GGAAACTGATGTCAAATACTTTCCCTTTGGTT"
                ),
                id="oryCun1.scaffold_133159",
                name="oryCun1.scaffold_133159",
                description="",
                annotations={
                    "start": 11087,
                    "srcSize": 13221,
                    "strand": 1,
                    "size": 164,
                },
            )

            recs = [rec1, rec2]

            fetched_recs = self.idx._get_record(34)

            self.compare_records(recs, fetched_recs)

        def test_records_end(self):
            rec1 = SeqRecord(
                Seq("TGTTTAGTACC----ATGCTTAGGAATGATAAACTCACTTAGTGtt"),
                id="mm9.chr10",
                name="mm9.chr10",
                description="",
                annotations={
                    "start": 3021494,
                    "srcSize": 129993255,
                    "strand": 1,
                    "size": 42,
                },
            )

            rec2 = SeqRecord(
                Seq("TGTTGCATGTCCTTTATTCTTTGGCGTGATAGGCTCACCCAATCTT"),
                id="ponAbe2.chr6",
                name="ponAbe2.chr6",
                description="",
                annotations={
                    "start": 16173516,
                    "srcSize": 174210431,
                    "strand": -1,
                    "size": 46,
                },
            )

            rec3 = SeqRecord(
                Seq("TGTTGCATATCCTTTATTCTTTGGCGTGATAGGCTCACCCAATCTT"),
                id="panTro2.chr6",
                name="panTro2.chr6",
                description="",
                annotations={
                    "start": 16393864,
                    "srcSize": 173908612,
                    "strand": -1,
                    "size": 46,
                },
            )

            rec4 = SeqRecord(
                Seq("TGTTGCATGTCGTTTATTCTTTGGCGTGATAGGCTCACCCAATCTT"),
                id="hg18.chr6",
                name="hg18.chr6",
                description="",
                annotations={
                    "start": 15875298,
                    "srcSize": 170899992,
                    "strand": -1,
                    "size": 46,
                },
            )

            rec5 = SeqRecord(
                Seq("TGTTAAGTCTCACTTGCTGTTCAAAGTGATAGCTTCACTCCATCAT"),
                id="canFam2.chr1",
                name="canFam2.chr1",
                description="",
                annotations={
                    "start": 78072287,
                    "srcSize": 125616256,
                    "strand": -1,
                    "size": 46,
                },
            )

            rec6 = SeqRecord(
                Seq("TGTTTAAAATG----ATTGCTAGAACTTCTA--CTCACTGGA----"),
                id="ornAna1.chr2",
                name="ornAna1.chr2",
                description="",
                annotations={
                    "start": 14757144,
                    "srcSize": 54797317,
                    "strand": -1,
                    "size": 36,
                },
            )

            recs = [rec1, rec2, rec3, rec4, rec5, rec6]

            fetched_recs = self.idx._get_record(99228)

            self.compare_records(recs, fetched_recs)

    class TestSearchGoodMAF(unittest.TestCase):
        """Test index searching on a properly-formatted MAF."""

        def setUp(self):
            self.idx = MafIndex(
                "MAF/ucsc_mm9_chr10.mafindex", "MAF/ucsc_mm9_chr10.maf", "mm9.chr10"
            )
            self.assertEqual(len(self.idx), 48)

        def test_invalid_type_1(self):
            search = self.idx.search((500, 1000), ("string", 1500))
            self.assertRaises(TypeError, next, search)

        def test_invalid_type_2(self):
            search = self.idx.search((500, 1000), (750, 1500.25))
            self.assertRaises(TypeError, next, search)

        def test_invalid_exon_count(self):
            search = self.idx.search((0, 1000, 2000), (500, 1500))
            self.assertRaises(ValueError, next, search)

        def test_invalid_exon_schema(self):
            search = self.idx.search((0, 1000, 2000), (250, 500, 2500))
            self.assertRaises(ValueError, next, search)

        def test_correct_retrieval_1(self):
            """Correct retrieval of Cnksr3 in mouse."""
            search = self.idx.search((3014742, 3018161), (3015028, 3018644))
            results = list(search)

            self.assertEqual(len(results), 4 + 4)

            self.assertEqual({len(x) for x in results}, {4, 1, 9, 10, 4, 3, 5, 1})

            # Code formatting note:
            # Expected start coordinates are grouped by alignment blocks
            # Turn black code style off
            # fmt: off
            self.assertEqual(
                {x.annotations["start"] for y in results for x in y},
                {
                    3014742, 6283, 184202, 1257,
                    3014778,
                    3014795, 184257, 6365, 15871286, 16389854, 16169492, 171521, 7816, 1309,
                    3014842, 1371, 7842, 171548, 16169512, 16389874, 15871306, 6404, 184317, 14750994,
                    3018161, 16390178, 15871611, 16169818,
                    3018230, 15871676, 16390243,
                    3018359, 16390338, 15871771, 184712, 16169976, 3018482
                }
            )
            # Turn black code style on
            # fmt: on

        def test_correct_retrieval_2(self):
            search = self.idx.search((3009319, 3021421), (3012566, 3021536))
            results = list(search)

            self.assertEqual(len(results), 6)

            self.assertEqual({len(x) for x in results}, {2, 4, 5, 14, 7, 6})

            # Code formatting note:
            # Expected start coordinates are grouped by alignment blocks
            # Turn black code style off
            # fmt: off
            self.assertEqual(
                {x.annotations["start"] for y in results for x in y},
                {
                    3009319, 11087,
                    3012076, 16160203, 16379004, 15860456,
                    3012441, 15860899, 16379447, 16160646, 180525,
                    3021421, 9910, 996, 16173434, 16393782, 15875216, 11047, 175213, 3552, 677, 78072203, 3590, 95587, 14757054,
                    3021465, 9957, 16173483, 16393831, 15875265, 78072243, 14757099,
                    3021494, 16173516, 16393864, 15875298, 78072287, 14757144
                }
            )
            # Turn black code style on
            # fmt: on

        def test_correct_retrieval_3(self):
            """Following issue 1083.

            https://github.com/biopython/biopython/issues/1083
            """
            search = self.idx.search(
                (3012076, 3012076 + 300), (3012076 + 100, 3012076 + 400)
            )
            results = list(search)

            self.assertEqual(len(results), 2)

            self.assertEqual({len(x) for x in results}, {4, 5})

            # Code formatting note:
            # Expected start coordinates are grouped by alignment blocks
            # Turn black code style off
            # fmt: off
            self.assertEqual(
                {x.annotations["start"] for y in results for x in y},
                {
                    3012076, 16160203, 16379004, 15860456,
                    3012441, 15860899, 16379447, 16160646, 180525
                }
            )
            # Turn black code style on
            # fmt: on

        def test_correct_block_boundary(self):
            """Following issues 504 and 1086.

            https://github.com/biopython/biopython/pull/504
            https://github.com/biopython/biopython/pull/1086#issuecomment-285080702

            We test what happens at the boundary between these two MAF blocks:

            a score=19159.000000
            s mm9.chr10                         3014644 45 + 129993255 CCTGTACC---CTTTGGTGAGAATTTTTGTTTCAGTGTTAAAAGTTTG
            s hg18.chr6                        15870786 46 - 170899992 CCTATACCTTTCTTTTATGAGAA-TTTTGTTTTAATCCTAAAC-TTTT
            i hg18.chr6                        I 9085 C 0
            s panTro2.chr6                     16389355 46 - 173908612 CCTATACCTTTCTTTTATGAGAA-TTTTGTTTTAATCCTAAAC-TTTT
            q panTro2.chr6                                             99999999999999999999999-9999999999999999999-9999
            i panTro2.chr6                     I 9106 C 0
            s calJac1.Contig6394                   6182 46 +    133105 CCTATACCTTTCTTTCATGAGAA-TTTTGTTTGAATCCTAAAC-TTTT
            i calJac1.Contig6394               N 0 C 0
            s loxAfr1.scaffold_75566               1167 34 -     10574 ------------TTTGGTTAGAA-TTATGCTTTAATTCAAAAC-TTCC
            q loxAfr1.scaffold_75566                                   ------------99999699899-9999999999999869998-9997
            i loxAfr1.scaffold_75566           N 0 C 0
            e tupBel1.scaffold_114895.1-498454   167376 4145 -    498454 I
            e echTel1.scaffold_288249             87661 7564 +    100002 I
            e otoGar1.scaffold_334.1-359464      181217 2931 -    359464 I
            e ponAbe2.chr6                     16161448 8044 - 174210431 I

            a score=40840.000000
            s mm9.chr10                         3014689 53 + 129993255 GGGAGCATAAAACTCTAAATCTGCTAAATGTCTTGTCCCT-TTGGAAAGAGTTG
            s hg18.chr6                        15870832 53 - 170899992 GGGATCATAAACCATTTAATCTGTGAAATATCTAATCTTT-TGGGAAATAGTGG
            i hg18.chr6                        C 0 I 401
            s panTro2.chr6                     16389401 53 - 173908612 GGGATCATAAACCATTTAATCTGTGAAATATCTAATCTTT-TGGGAAATAGTGG
            q panTro2.chr6                                             9999999999999999999999999999999999999999-9999999999999
            i panTro2.chr6                     C 0 I 400
            s calJac1.Contig6394                   6228 53 +    133105 GGGATCATAAGCCATTTAATCTGTGAAATGTGAAATCTTT-TGGGAAACAGTGG
            i calJac1.Contig6394               C 0 I 2
            s otoGar1.scaffold_334.1-359464      184148 52 -    359464 GGAAGCATAAACT-TTTAATCTATGAAATATCAAATCACT-TGGGCAATAGCTG
            q otoGar1.scaffold_334.1-359464                            7455455669566-99665699769895555689997599-9984787795599
            i otoGar1.scaffold_334.1-359464    I 2931 I 2
            s loxAfr1.scaffold_75566               1201 54 -     10574 GGGAGTATAAACCATTTAGTCTGCGAAATGCCAAATCTTCAGGGGAAAAAGCTG
            q loxAfr1.scaffold_75566                                   899989799999979999999999999999797999999999999999999999
            i loxAfr1.scaffold_75566           C 0 I 2
            e tupBel1.scaffold_114895.1-498454   167376 4145 -    498454 I
            e echTel1.scaffold_288249             87661 7564 +    100002 I
            e ponAbe2.chr6                     16161448 8044 - 174210431 I
            """
            # Segments ending at the end of the first block
            search = self.idx.search([3014687], [3014689])
            self.assertEqual(len(list(search)), 1)
            search = self.idx.search([3014688], [3014689])
            self.assertEqual(len(list(search)), 1)

            # Segments starting at the beginning of the second block
            search = self.idx.search([3014689], [3014690])
            self.assertEqual(len(list(search)), 1)
            search = self.idx.search([3014689], [3014691])
            self.assertEqual(len(list(search)), 1)

            # Segments overlapping the 2 blocks
            search = self.idx.search([3014688], [3014690])
            self.assertEqual(len(list(search)), 2)
            search = self.idx.search([3014687], [3014690])
            self.assertEqual(len(list(search)), 2)
            search = self.idx.search([3014687], [3014691])
            self.assertEqual(len(list(search)), 2)

        def test_correct_block_length(self):
            """Following issues 504 and 1086.

            https://github.com/biopython/biopython/pull/504
            https://github.com/biopython/biopython/pull/1086#issuecomment-285080702

            We get the alignment corresponding to the following whole MAF block
            and check that the lengths of its sequences are correct:

            a score=40840.000000
            s mm9.chr10                         3014689 53 + 129993255 GGGAGCATAAAACTCTAAATCTGCTAAATGTCTTGTCCCT-TTGGAAAGAGTTG
            s hg18.chr6                        15870832 53 - 170899992 GGGATCATAAACCATTTAATCTGTGAAATATCTAATCTTT-TGGGAAATAGTGG
            i hg18.chr6                        C 0 I 401
            s panTro2.chr6                     16389401 53 - 173908612 GGGATCATAAACCATTTAATCTGTGAAATATCTAATCTTT-TGGGAAATAGTGG
            q panTro2.chr6                                             9999999999999999999999999999999999999999-9999999999999
            i panTro2.chr6                     C 0 I 400
            s calJac1.Contig6394                   6228 53 +    133105 GGGATCATAAGCCATTTAATCTGTGAAATGTGAAATCTTT-TGGGAAACAGTGG
            i calJac1.Contig6394               C 0 I 2
            s otoGar1.scaffold_334.1-359464      184148 52 -    359464 GGAAGCATAAACT-TTTAATCTATGAAATATCAAATCACT-TGGGCAATAGCTG
            q otoGar1.scaffold_334.1-359464                            7455455669566-99665699769895555689997599-9984787795599
            i otoGar1.scaffold_334.1-359464    I 2931 I 2
            s loxAfr1.scaffold_75566               1201 54 -     10574 GGGAGTATAAACCATTTAGTCTGCGAAATGCCAAATCTTCAGGGGAAAAAGCTG
            q loxAfr1.scaffold_75566                                   899989799999979999999999999999797999999999999999999999
            i loxAfr1.scaffold_75566           C 0 I 2
            e tupBel1.scaffold_114895.1-498454   167376 4145 -    498454 I
            e echTel1.scaffold_288249             87661 7564 +    100002 I
            e ponAbe2.chr6                     16161448 8044 - 174210431 I
            """
            ali = self.idx.get_spliced([3014689], [3014689 + 53])
            seq_dict = {seqrec.id: seqrec.seq for seqrec in ali}
            correct_lengths = {
                "mm9.chr10": 53,
                "hg18.chr6": 53,
                "panTro2.chr6": 53,
                "calJac1.Contig6394": 53,
                "otoGar1.scaffold_334.1-359464": 52,
                "loxAfr1.scaffold_75566": 54,
            }
            for seq_id, length in correct_lengths.items():
                self.assertEqual(len(seq_dict[seq_id].replace("-", "")), length)

        def test_correct_spliced_sequences_1(self):
            """Checking that spliced sequences are correct.

            We get the alignment corresponding to the following whole MAF block
            and check that the sequences are correct:

            a score=40840.000000
            s mm9.chr10                         3014689 53 + 129993255 GGGAGCATAAAACTCTAAATCTGCTAAATGTCTTGTCCCT-TTGGAAAGAGTTG
            s hg18.chr6                        15870832 53 - 170899992 GGGATCATAAACCATTTAATCTGTGAAATATCTAATCTTT-TGGGAAATAGTGG
            i hg18.chr6                        C 0 I 401
            s panTro2.chr6                     16389401 53 - 173908612 GGGATCATAAACCATTTAATCTGTGAAATATCTAATCTTT-TGGGAAATAGTGG
            q panTro2.chr6                                             9999999999999999999999999999999999999999-9999999999999
            i panTro2.chr6                     C 0 I 400
            s calJac1.Contig6394                   6228 53 +    133105 GGGATCATAAGCCATTTAATCTGTGAAATGTGAAATCTTT-TGGGAAACAGTGG
            i calJac1.Contig6394               C 0 I 2
            s otoGar1.scaffold_334.1-359464      184148 52 -    359464 GGAAGCATAAACT-TTTAATCTATGAAATATCAAATCACT-TGGGCAATAGCTG
            q otoGar1.scaffold_334.1-359464                            7455455669566-99665699769895555689997599-9984787795599
            i otoGar1.scaffold_334.1-359464    I 2931 I 2
            s loxAfr1.scaffold_75566               1201 54 -     10574 GGGAGTATAAACCATTTAGTCTGCGAAATGCCAAATCTTCAGGGGAAAAAGCTG
            q loxAfr1.scaffold_75566                                   899989799999979999999999999999797999999999999999999999
            i loxAfr1.scaffold_75566           C 0 I 2
            e tupBel1.scaffold_114895.1-498454   167376 4145 -    498454 I
            e echTel1.scaffold_288249             87661 7564 +    100002 I
            e ponAbe2.chr6                     16161448 8044 - 174210431 I
            """
            ali = self.idx.get_spliced([3014689], [3014689 + 53])
            seq_dict = {seqrec.id: seqrec.seq for seqrec in ali}
            correct_sequences = {
                "mm9.chr10": "GGGAGCATAAAACTCTAAATCTGCTAAATGTCTTGTCCCTTTGGAAAGAGTTG",
                "hg18.chr6": "GGGATCATAAACCATTTAATCTGTGAAATATCTAATCTTTTGGGAAATAGTGG",
                "panTro2.chr6": "GGGATCATAAACCATTTAATCTGTGAAATATCTAATCTTTTGGGAAATAGTGG",
                "calJac1.Contig6394": "GGGATCATAAGCCATTTAATCTGTGAAATGTGAAATCTTTTGGGAAACAGTGG",
                "otoGar1.scaffold_334.1-359464": "GGAAGCATAAACTTTTAATCTATGAAATATCAAATCACTTGGGCAATAGCTG",
                "loxAfr1.scaffold_75566": "GGGAGTATAAACCATTTAGTCTGCGAAATGCCAAATCTTCAGGGGAAAAAGCTG",
            }
            for seq_id, sequence in correct_sequences.items():
                self.assertEqual(seq_dict[seq_id].replace("-", ""), sequence)

        def test_correct_spliced_sequences_2(self):
            """Checking that spliced sequences are correct.

            We get spliced alignements from following MAF blocks
            and check that the sequences are correct:

            a score=19159.000000
            s mm9.chr10                         3014644 45 + 129993255 CCTGTACC---CTTTGGTGAGAATTTTTGTTTCAGTGTTAAAAGTTTG
            s hg18.chr6                        15870786 46 - 170899992 CCTATACCTTTCTTTTATGAGAA-TTTTGTTTTAATCCTAAAC-TTTT
            i hg18.chr6                        I 9085 C 0
            s panTro2.chr6                     16389355 46 - 173908612 CCTATACCTTTCTTTTATGAGAA-TTTTGTTTTAATCCTAAAC-TTTT
            q panTro2.chr6                                             99999999999999999999999-9999999999999999999-9999
            i panTro2.chr6                     I 9106 C 0
            s calJac1.Contig6394                   6182 46 +    133105 CCTATACCTTTCTTTCATGAGAA-TTTTGTTTGAATCCTAAAC-TTTT
            i calJac1.Contig6394               N 0 C 0
            s loxAfr1.scaffold_75566               1167 34 -     10574 ------------TTTGGTTAGAA-TTATGCTTTAATTCAAAAC-TTCC
            q loxAfr1.scaffold_75566                                   ------------99999699899-9999999999999869998-9997
            i loxAfr1.scaffold_75566           N 0 C 0
            e tupBel1.scaffold_114895.1-498454   167376 4145 -    498454 I
            e echTel1.scaffold_288249             87661 7564 +    100002 I
            e otoGar1.scaffold_334.1-359464      181217 2931 -    359464 I
            e ponAbe2.chr6                     16161448 8044 - 174210431 I

            a score=40840.000000
            s mm9.chr10                         3014689 53 + 129993255 GGGAGCATAAAACTCTAAATCTGCTAAATGTCTTGTCCCT-TTGGAAAGAGTTG
            s hg18.chr6                        15870832 53 - 170899992 GGGATCATAAACCATTTAATCTGTGAAATATCTAATCTTT-TGGGAAATAGTGG
            i hg18.chr6                        C 0 I 401
            s panTro2.chr6                     16389401 53 - 173908612 GGGATCATAAACCATTTAATCTGTGAAATATCTAATCTTT-TGGGAAATAGTGG
            q panTro2.chr6                                             9999999999999999999999999999999999999999-9999999999999
            i panTro2.chr6                     C 0 I 400
            s calJac1.Contig6394                   6228 53 +    133105 GGGATCATAAGCCATTTAATCTGTGAAATGTGAAATCTTT-TGGGAAACAGTGG
            i calJac1.Contig6394               C 0 I 2
            s otoGar1.scaffold_334.1-359464      184148 52 -    359464 GGAAGCATAAACT-TTTAATCTATGAAATATCAAATCACT-TGGGCAATAGCTG
            q otoGar1.scaffold_334.1-359464                            7455455669566-99665699769895555689997599-9984787795599
            i otoGar1.scaffold_334.1-359464    I 2931 I 2
            s loxAfr1.scaffold_75566               1201 54 -     10574 GGGAGTATAAACCATTTAGTCTGCGAAATGCCAAATCTTCAGGGGAAAAAGCTG
            q loxAfr1.scaffold_75566                                   899989799999979999999999999999797999999999999999999999
            i loxAfr1.scaffold_75566           C 0 I 2
            e tupBel1.scaffold_114895.1-498454   167376 4145 -    498454 I
            e echTel1.scaffold_288249             87661 7564 +    100002 I
            e ponAbe2.chr6                     16161448 8044 - 174210431 I
            """
            ali = self.idx.get_spliced([3014644, 3014689], [3014644 + 45, 3014689 + 53])
            seq_dict = {seqrec.id: seqrec.seq for seqrec in ali}
            correct_sequences = {
                "mm9.chr10": "CCTGTACCCTTTGGTGAGAATTTTTGTTTCAGTGTTAAAAGTTTGGGGAGCATAAAACTCTAAATCTGCTAAATGTCTTGTCCCTTTGGAAAGAGTTG",
                "hg18.chr6": "CCTATACCTTTCTTTTATGAGAATTTTGTTTTAATCCTAAACTTTTGGGATCATAAACCATTTAATCTGTGAAATATCTAATCTTTTGGGAAATAGTGG",
                "panTro2.chr6": "CCTATACCTTTCTTTTATGAGAATTTTGTTTTAATCCTAAACTTTTGGGATCATAAACCATTTAATCTGTGAAATATCTAATCTTTTGGGAAATAGTGG",
                "calJac1.Contig6394": "CCTATACCTTTCTTTCATGAGAATTTTGTTTGAATCCTAAACTTTTGGGATCATAAGCCATTTAATCTGTGAAATGTGAAATCTTTTGGGAAACAGTGG",
                "otoGar1.scaffold_334.1-359464": "GGAAGCATAAACTTTTAATCTATGAAATATCAAATCACTTGGGCAATAGCTG",
                "loxAfr1.scaffold_75566": "TTTGGTTAGAATTATGCTTTAATTCAAAACTTCCGGGAGTATAAACCATTTAGTCTGCGAAATGCCAAATCTTCAGGGGAAAAAGCTG",
            }
            for seq_id, sequence in correct_sequences.items():
                self.assertEqual(seq_dict[seq_id].replace("-", ""), sequence)

    class TestSearchBadMAF(unittest.TestCase):
        """Test index searching on an incorrectly-formatted MAF."""

        def setUp(self):
            self.idx = MafIndex(
                "MAF/ucsc_mm9_chr10_bad.mafindex",
                "MAF/ucsc_mm9_chr10_bad.maf",
                "mm9.chr10",
            )
            self.assertEqual(len(self.idx), 48)

        def test_incorrect_bundle_coords(self):
            search = self.idx.search((3013219,), (3013319,))
            self.assertRaises(ValueError, next, search)

    class TestSpliceGoodMAF(unittest.TestCase):
        """Test in silico splicing on a correctly-formatted MAF."""

        def setUp(self):
            self.idx = MafIndex(
                "MAF/ucsc_mm9_chr10_big.mafindex",
                "MAF/ucsc_mm9_chr10_big.maf",
                "mm9.chr10",
            )
            self.assertEqual(len(self.idx), 983)

        def test_invalid_strand(self):
            self.assertRaises(
                ValueError, self.idx.get_spliced, (0, 1000), (500, 1500), "."
            )

        def test_no_alignment(self):
            result = self.idx.get_spliced((0, 1000), (500, 1500), 1)

            self.assertEqual(len(result), 1)
            self.assertEqual(result[0].seq, "N" * 1000)

        def test_correct_retrieval_1(self):
            """Correct retrieval of Cnksr3 in mouse.

            This is the real thing. We're pulling the spliced alignment of
            an actual gene (Cnksr3) in mouse. It should perfectly match the
            spliced transcript pulled independently from UCSC.
            """
            if sys.platform == "win32":
                # TODO - fix this hack to get other tests to pass
                # See https://github.com/biopython/biopython/issues/3640
                # and https://github.com/biopython/biopython/issues/1149
                return
            result = self.idx.get_spliced(
                (
                    3134303,
                    3185733,
                    3192055,
                    3193589,
                    3203538,
                    3206102,
                    3208126,
                    3211424,
                    3211872,
                    3217393,
                    3219697,
                    3220356,
                    3225954,
                ),
                (
                    3134909,
                    3185897,
                    3192258,
                    3193677,
                    3203580,
                    3206222,
                    3208186,
                    3211493,
                    3212019,
                    3217518,
                    3219906,
                    3220446,
                    3227479,
                ),
                1,
            )

            cnksr3 = SeqIO.read("MAF/cnksr3.fa", "fasta").seq.upper()
            mm9_seq = "".join(
                [str(x.seq) for x in result if x.id.startswith("mm9")]
            ).replace("-", "")

            self.assertEqual(mm9_seq, cnksr3)

    class TestSpliceBadMAF(unittest.TestCase):
        """Test in silico splicing on an incorrectly-formatted MAF."""

        def setUp(self):
            self.idx = MafIndex(
                "MAF/ucsc_mm9_chr10_bad.mafindex",
                "MAF/ucsc_mm9_chr10_bad.maf",
                "mm9.chr10",
            )
            self.assertEqual(len(self.idx), 48)

        def test_inconsistent_strand(self):
            self.assertRaises(
                ValueError, self.idx.get_spliced, (0, 3021421), (1000, 3022000), 1
            )

        def test_bundle_without_target(self):
            self.assertRaises(
                ValueError, self.idx.get_spliced, (3009319,), (3009900,), 1
            )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
