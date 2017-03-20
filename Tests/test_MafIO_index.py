# Copyright 2012 by Andrew Sczesnak.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for Bio.AlignIO.MafIO.MafIndex()"""

try:
    import sqlite3
except ImportError:
    # skip most tests if sqlite is not available
    sqlite3 = None

import os
import unittest
import tempfile
import shutil

from Bio.AlignIO.MafIO import MafIndex
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from seq_tests_common import compare_record


class StaticMethodTest(unittest.TestCase):
    """Test static UCSC binning-related functions"""

    def test_region2bin(self):
        data = [(25079603, 25079787, set([0, 1, 11, 96, 776])),
                (25128173, 25128248, set([0, 1, 11, 96, 776])),
                (50312474, 50312703, set([0, 1, 968, 14, 120])),
                (41905591, 41906101, set([0, 1, 904, 13, 112])),
                (16670899, 16673060, set([0, 1, 10, 712, 88])),
                (75495356, 75495494, set([0, 1, 2, 1160, 144, 17])),
                (92259501, 92261053, set([0, 1, 2, 1288, 160, 19])),
                (83834063, 83838132, set([0, 1, 2, 1224, 18, 152])),
                (7309597, 7310411, set([0, 1, 640, 79, 9])),
                (6190410, 6190999, set([0, 1, 632, 78, 9]))]

        for x, y, z in data:
            self.assertEqual(MafIndex._region2bin(x, y), z)

        for x, y, z in data:
            self.assertRaises(TypeError, MafIndex._region2bin, str(x), str(y))

    def test_ucscbin(self):
        data = [(25079603, 25079787, 776),
                (25128173, 25128248, 776),
                (50312474, 50312703, 968),
                (41905591, 41906101, 904),
                (16670899, 16673060, 712),
                (75495356, 75495494, 1160),
                (92259501, 92261053, 1288),
                (83834063, 83838132, 1224),
                (7309597, 7310411, 640),
                (6190410, 6190999, 632)]

        for x, y, z in data:
            self.assertEqual(MafIndex._ucscbin(x, y), z)

        for x, y, z in data:
            self.assertRaises(TypeError, MafIndex._ucscbin, str(x), str(y))


if sqlite3:
    class PreBuiltIndexTest(unittest.TestCase):
        """Test loading of prebuilt indices"""

        def test_old(self):
            idx = MafIndex("MAF/ucsc_mm9_chr10.mafindex",
                           "MAF/ucsc_mm9_chr10.maf", "mm9.chr10")
            self.assertEqual(len(idx), 48)

        def test_old_wrong_target_seqname(self):
            self.assertRaises(ValueError,
                              MafIndex,
                              "MAF/ucsc_mm9_chr10.mafindex",
                              "MAF/ucsc_mm9_chr10.maf",
                              "mm9.chr11")

        def test_old_wrong_filename(self):
            self.assertRaises(ValueError,
                              MafIndex,
                              "MAF/ucsc_mm9_chr10.mafindex",
                              "MAF/humor.maf",
                              "mm9.chr10")

        def test_old_file_not_found(self):
            # TODO: Switch to FileNotFoundError once we drop Python 2 support
            # Under Python 2, we expect IOError.
            # Under Python 3, we expect FileNotFoundError which is a subclass
            # of OSError, and that has IOError as an alias so this works.
            self.assertRaises(IOError,
                              MafIndex,
                              "MAF/ucsc_mm9_chr11.mafindex",
                              "MAF/ucsc_mm9_chr11.maf",
                              "mm9.chr11")

        def test_old_wrong_version(self):
            self.assertRaises(ValueError,
                              MafIndex,
                              "MAF/wrong_version.idx",
                              "MAF/ucsc_mm9_chr10.maf",
                              "mm9.chr10")

        def test_old_unfinished_index(self):
            self.assertRaises(ValueError,
                              MafIndex,
                              "MAF/unfinished.idx",
                              "MAF/ucsc_mm9_chr10.maf",
                              "mm9.chr10")

        def test_old_corrupt_index(self):
            self.assertRaises(ValueError,
                              MafIndex,
                              "MAF/corrupt.idx",
                              "MAF/ucsc_mm9_chr10.maf",
                              "mm9.chr10")

        def test_old_invalid_sqlite(self):
            self.assertRaises(ValueError,
                              MafIndex,
                              "MAF/invalid.idx",
                              "MAF/ucsc_mm9_chr10.maf",
                              "mm9.chr10")

    class NewIndexTest(unittest.TestCase):
        """Test creation of new indices"""

        def setUp(self):
            self.tmpdir = tempfile.mkdtemp()
            self.tmpfile = self.tmpdir + "/database.sqlite3"

        def tearDown(self):
            if os.path.isdir(self.tmpdir):
                shutil.rmtree(self.tmpdir)

        def test_good_small(self):
            idx = MafIndex(self.tmpfile, "MAF/ucsc_mm9_chr10.maf", "mm9.chr10")
            self.assertEquals(len(idx), 48)

        def test_good_big(self):
            idx = MafIndex(self.tmpfile, "MAF/ucsc_mm9_chr10_big.maf", "mm9.chr10")
            self.assertEquals(len(idx), 983)

        def test_bundle_without_target(self):
            self.assertRaises(ValueError,
                              MafIndex,
                              self.tmpfile,
                              "MAF/bundle_without_target.maf",
                              "mm9.chr10")

        def test_length_coords_mismatch(self):
            self.assertRaises(ValueError,
                              MafIndex,
                              self.tmpfile,
                              "MAF/length_coords_mismatch.maf",
                              "mm9.chr10")

    class TestGetRecord(unittest.TestCase):
        """Make sure we can seek and fetch records properly"""

        def setUp(self):
            self.idx = MafIndex("MAF/ucsc_mm9_chr10.mafindex",
                                "MAF/ucsc_mm9_chr10.maf", "mm9.chr10")
            self.assertEqual(len(self.idx), 48)

        def test_records_begin(self):
            recs = {}

            recs[0] = SeqRecord(Seq("TCATAGGTATTTATTTTTAAATATGGTTTGCTTTATGGCTAGAA"
                                    "CACACCGATTACTTAAAATAGGATTAACC--CCCATACACTTTA"
                                    "AAAATGATTAAACAACATTTCTGCTGCTCGCTCACATTCTTCAT"
                                    "AGAAGATGACATAATGTATTTTCCTTTTGGTT"),
                                id="mm9.chr10",
                                name="mm9.chr10",
                                description="",
                                annotations={"start": 3009319,
                                             "srcSize": 129993255,
                                             "strand": 1,
                                             "size": 162})

            recs[1] = SeqRecord(Seq("TCACAGATATTTACTATTAAATATGGTTTGTTATATGGTTACGG"
                                    "TTCATAGGTTACTTGGAATTGGATTAACCTTCTTATTCATTGCA"
                                    "GAATTGGTTACACTGTGTTCTTGACCTTTGCTTGTTTTCTCCAT"
                                    "GGAAACTGATGTCAAATACTTTCCCTTTGGTT"),
                                id="oryCun1.scaffold_133159",
                                name="oryCun1.scaffold_133159",
                                description="",
                                annotations={"start": 11087,
                                             "srcSize": 13221,
                                             "strand": 1,
                                             "size": 164})

            fetched_recs = self.idx._get_record(34)

            for i in range(2):
                self.assertTrue(compare_record(recs[i], fetched_recs[i]))

        def test_records_end(self):
            recs = {}

            recs[0] = SeqRecord(Seq("TGTTTAGTACC----ATGCTTAGGAATGATAAACTCACTTAGTGtt"),
                                id="mm9.chr10",
                                name="mm9.chr10",
                                description="",
                                annotations={"start": 3021494,
                                             "srcSize": 129993255,
                                             "strand": 1,
                                             "size": 42})

            recs[1] = SeqRecord(Seq("TGTTGCATGTCCTTTATTCTTTGGCGTGATAGGCTCACCCAATCTT"),
                                id="ponAbe2.chr6",
                                name="ponAbe2.chr6",
                                description="",
                                annotations={"start": 16173516,
                                             "srcSize": 174210431,
                                             "strand": -1,
                                             "size": 46})

            recs[2] = SeqRecord(Seq("TGTTGCATATCCTTTATTCTTTGGCGTGATAGGCTCACCCAATCTT"),
                                id="panTro2.chr6",
                                name="panTro2.chr6",
                                description="",
                                annotations={"start": 16393864,
                                             "srcSize": 173908612,
                                             "strand": -1,
                                             "size": 46})

            recs[3] = SeqRecord(Seq("TGTTGCATGTCGTTTATTCTTTGGCGTGATAGGCTCACCCAATCTT"),
                                id="hg18.chr6",
                                name="hg18.chr6",
                                description="",
                                annotations={"start": 15875298,
                                             "srcSize": 170899992,
                                             "strand": -1,
                                             "size": 46})

            recs[4] = SeqRecord(Seq("TGTTAAGTCTCACTTGCTGTTCAAAGTGATAGCTTCACTCCATCAT"),
                                id="canFam2.chr1",
                                name="canFam2.chr1",
                                description="",
                                annotations={"start": 78072287,
                                             "srcSize": 125616256,
                                             "strand": -1,
                                             "size": 46})

            recs[5] = SeqRecord(Seq("TGTTTAAAATG----ATTGCTAGAACTTCTA--CTCACTGGA----"),
                                id="ornAna1.chr2",
                                name="ornAna1.chr2",
                                description="",
                                annotations={"start": 14757144,
                                             "srcSize": 54797317,
                                             "strand": -1,
                                             "size": 36})

            fetched_recs = self.idx._get_record(99228)

            for i in range(6):
                self.assertTrue(compare_record(recs[i], fetched_recs[i]))

    class TestSearchGoodMAF(unittest.TestCase):
        """Test index searching on a properly-formatted MAF"""

        def setUp(self):
            self.idx = MafIndex("MAF/ucsc_mm9_chr10.mafindex",
                                "MAF/ucsc_mm9_chr10.maf", "mm9.chr10")
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
            search = self.idx.search((3014742, 3018161), (3015028, 3018644))
            results = [x for x in search]

            self.assertEqual(len(results), 12)

            self.assertEqual(set([len(x) for x in results]),
                             set([5, 10, 7, 6, 3, 1, 1, 1, 2, 4, 4, 9]))

            self.assertEqual(set([x.annotations["start"] for y in results
                                  for x in y]),
                             set([3018359, 16390338, 15871771, 184712,
                                  16169512, 16169976, 3014842, 1371, 7842,
                                  171548, 16389874, 15871306, 6404, 184317,
                                  14750994, 3015028, 1616, 8040, 171763,
                                  16169731, 6627, 184539, 3014689, 15870832,
                                  16389401, 6228, 184148, 1201, 3018230,
                                  15871676, 16390243, 3014778, 3018482, 3017743,
                                  3018644, 78070420, 3014742, 6283, 184202,
                                  1257, 3018161, 16390178, 15871611, 16169818,
                                  3014795, 184257, 6365, 15871286, 16389854,
                                  16169492, 171521, 7816, 1309]))

        def test_correct_retrieval_2(self):
            search = self.idx.search((3009319, 3021421), (3012566, 3021536))
            results = [x for x in search]

            self.assertEqual(len(results), 8)

            self.assertEqual(set([len(x) for x in results]),
                             set([14, 5, 2, 6, 7, 15, 6, 4]))

            self.assertEqual(set([x.annotations["start"] for y in results
                                  for x in y]),
                             set([3021421, 9910, 996, 16173434, 16393782,
                                  15875216, 11047, 175213, 3552, 677, 78072203,
                                  3590, 95587, 14757054, 3012441, 15860899,
                                  16379447, 16160646, 180525, 3009319, 11087,
                                  3012566, 15861013, 16379561, 16160760, 180626,
                                  310, 3021465, 9957, 16173483, 16393831,
                                  15875265, 78072243, 14757099, 3021275, 9741,
                                  838, 16173265, 16393613, 15875047, 10878,
                                  175057, 3382, 521, 78072035, 73556, 3422,
                                  95418, 14756885, 3021494, 16173516, 16393864,
                                  15875298, 78072287, 14757144, 3012076,
                                  16160203, 16379004, 15860456]))

        def test_correct_retrieval_3(self):
            search = self.idx.search((3012076, 3012076 + 300), (3012076 + 100, 3012076 + 400))
            results = [x for x in search]

            self.assertEqual(len(results), 2)

            self.assertEqual(set([len(x) for x in results]),
                             set([4, 5]))

            # Code formatting note:
            # Expected start coordinates are grouped by alignment blocks
            self.assertEqual(
                set([x.annotations["start"] for y in results for x in y]),
                set([
                    3012076, 16160203, 16379004, 15860456,
                    3012441, 15860899, 16379447, 16160646, 180525]))

    class TestSearchBadMAF(unittest.TestCase):
        """Test index searching on an incorrectly-formatted MAF"""

        def setUp(self):
            self.idx = MafIndex("MAF/ucsc_mm9_chr10_bad.mafindex",
                                "MAF/ucsc_mm9_chr10_bad.maf", "mm9.chr10")
            self.assertEqual(len(self.idx), 48)

        def test_incorrect_bundle_coords(self):
            search = self.idx.search((3013219,), (3013319,))
            self.assertRaises(ValueError, next, search)

    class TestSpliceGoodMAF(unittest.TestCase):
        """Test in silico splicing on a correctly-formatted MAF"""

        def setUp(self):
            self.idx = MafIndex("MAF/ucsc_mm9_chr10_big.mafindex",
                                "MAF/ucsc_mm9_chr10_big.maf", "mm9.chr10")
            self.assertEqual(len(self.idx), 983)

        def test_invalid_strand(self):
            self.assertRaises(ValueError,
                              self.idx.get_spliced,
                              (0, 1000), (500, 1500), ".")

        def test_no_alignment(self):
            result = self.idx.get_spliced((0, 1000), (500, 1500), 1)

            self.assertEqual(len(result), 1)
            self.assertEqual(len(result[0].seq), 1000)
            self.assertEqual(str(result[0].seq), "N" * 1000)

        def test_correct_retrieval_1(self):
            """
            This is the real thing. We're pulling the spliced alignment of
            an actual gene (Cnksr3) in mouse. It should perfectly match the
            spliced transcript pulled independently from UCSC.
            """

            result = self.idx.get_spliced((3134303, 3185733, 3192055, 3193589,
                                           3203538, 3206102, 3208126, 3211424,
                                           3211872, 3217393, 3219697, 3220356,
                                           3225954),
                                          (3134909, 3185897, 3192258, 3193677,
                                           3203580, 3206222, 3208186, 3211493,
                                           3212019, 3217518, 3219906, 3220446,
                                           3227479), 1)

            cnksr3 = str(SeqIO.read("MAF/cnksr3.fa", "fasta").seq).upper()
            mm9_seq = "".join([str(x.seq) for x in result
                               if x.id.startswith("mm9")]).replace("-", "")

            self.assertEqual(mm9_seq, cnksr3)

    class TestSpliceBadMAF(unittest.TestCase):
        """Test in silico splicing on an incorrectly-formatted MAF"""

        def setUp(self):
            self.idx = MafIndex("MAF/ucsc_mm9_chr10_bad.mafindex",
                                "MAF/ucsc_mm9_chr10_bad.maf", "mm9.chr10")
            self.assertEqual(len(self.idx), 48)

        def test_inconsistent_strand(self):
            self.assertRaises(ValueError,
                              self.idx.get_spliced,
                              (0, 3021421), (1000, 3022000), 1)

        def test_bundle_without_target(self):
            self.assertRaises(ValueError,
                              self.idx.get_spliced,
                              (3009319,), (3009900,), 1)

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
