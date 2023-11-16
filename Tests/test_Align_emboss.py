# Copyright 2008-2014 by Peter Cock.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.Align.emboss module."""
import unittest

from io import StringIO

from Bio import Align


try:
    import numpy as np
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align.emboss."
    ) from None


class TestEmboss(unittest.TestCase):
    def test_pair_example(self):
        # Alignment file obtained from EMBOSS:
        # http://emboss.sourceforge.net/docs/themes/alnformats/align.pair
        path = "Emboss/water.txt"
        alignments = Align.parse(path, "emboss")
        self.assertEqual(alignments.metadata["Program"], "water")
        self.assertEqual(alignments.metadata["Rundate"], "Wed Jan 16 17:23:19 2002")
        self.assertEqual(alignments.metadata["Report_file"], "stdout")
        alignment = next(alignments)
        self.assertEqual(alignment.annotations["Matrix"], "EBLOSUM62")
        self.assertAlmostEqual(alignment.annotations["Gap_penalty"], 10.0)
        self.assertAlmostEqual(alignment.annotations["Extend_penalty"], 0.5)
        self.assertEqual(alignment.annotations["Identity"], 112)
        self.assertEqual(alignment.annotations["Similarity"], 112)
        self.assertEqual(alignment.annotations["Gaps"], 19)
        self.assertAlmostEqual(alignment.annotations["Score"], 591.5)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 131))
        self.assertEqual(alignment.sequences[0].id, "IXI_234")
        self.assertEqual(alignment.sequences[1].id, "IXI_235")
        self.assertEqual(
            alignment.sequences[0].seq,
            "TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "TSPASIRPPAGPSSRRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                np.array([[0, 15, 24, 74, 84, 131], [0, 15, 15, 65, 65, 112]]),
            )
        )
        self.assertEqual(
            alignment[0],
            "TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertEqual(
            alignment[1],
            "TSPASIRPPAGPSSR---------RPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGW----------RASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            "|||||||||||||||         ||||||||||||||||||||||||||||||||||||||||||||||||||          |||||||||||||||||||||||||||||||||||||||||||||||",
        )
        self.assertEqual(
            str(alignment),
            """\
IXI_234           0 TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTC
                  0 |||||||||||||||---------||||||||||||||||||||||||||||||||||||
IXI_235           0 TSPASIRPPAGPSSR---------RPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTC

IXI_234          60 TTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTG
                 60 ||||||||||||||----------||||||||||||||||||||||||||||||||||||
IXI_235          51 TTSTSTRHRGRSGW----------RASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTG

IXI_234         120 PPAWAGDRSHE 131
                120 ||||||||||| 131
IXI_235         101 PPAWAGDRSHE 112
""",
        )
        self.assertTrue(
            np.array_equal(
                np.array(alignment, "U"),
                # fmt: off
np.array([['T', 'S', 'P', 'A', 'S', 'I', 'R', 'P', 'P', 'A', 'G', 'P', 'S', 'S',
           'R', 'P', 'A', 'M', 'V', 'S', 'S', 'R', 'R', 'T', 'R', 'P', 'S', 'P',
           'P', 'G', 'P', 'R', 'R', 'P', 'T', 'G', 'R', 'P', 'C', 'C', 'S', 'A',
           'A', 'P', 'R', 'R', 'P', 'Q', 'A', 'T', 'G', 'G', 'W', 'K', 'T', 'C',
           'S', 'G', 'T', 'C', 'T', 'T', 'S', 'T', 'S', 'T', 'R', 'H', 'R', 'G',
           'R', 'S', 'G', 'W', 'S', 'A', 'R', 'T', 'T', 'T', 'A', 'A', 'C', 'L',
           'R', 'A', 'S', 'R', 'K', 'S', 'M', 'R', 'A', 'A', 'C', 'S', 'R', 'S',
           'A', 'G', 'S', 'R', 'P', 'N', 'R', 'F', 'A', 'P', 'T', 'L', 'M', 'S',
           'S', 'C', 'I', 'T', 'S', 'T', 'T', 'G', 'P', 'P', 'A', 'W', 'A', 'G',
           'D', 'R', 'S', 'H', 'E'],
          ['T', 'S', 'P', 'A', 'S', 'I', 'R', 'P', 'P', 'A', 'G', 'P', 'S', 'S',
           'R', '-', '-', '-', '-', '-', '-', '-', '-', '-', 'R', 'P', 'S', 'P',
           'P', 'G', 'P', 'R', 'R', 'P', 'T', 'G', 'R', 'P', 'C', 'C', 'S', 'A',
           'A', 'P', 'R', 'R', 'P', 'Q', 'A', 'T', 'G', 'G', 'W', 'K', 'T', 'C',
           'S', 'G', 'T', 'C', 'T', 'T', 'S', 'T', 'S', 'T', 'R', 'H', 'R', 'G',
           'R', 'S', 'G', 'W', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
           'R', 'A', 'S', 'R', 'K', 'S', 'M', 'R', 'A', 'A', 'C', 'S', 'R', 'S',
           'A', 'G', 'S', 'R', 'P', 'N', 'R', 'F', 'A', 'P', 'T', 'L', 'M', 'S',
           'S', 'C', 'I', 'T', 'S', 'T', 'T', 'G', 'P', 'P', 'A', 'W', 'A', 'G',
           'D', 'R', 'S', 'H', 'E']], dtype='U')
                # fmt: on
            )
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_local_water2(self):
        """Test parsing a local alignment."""
        path = "Emboss/water2.txt"
        alignments = Align.parse(path, "emboss")
        self.assertEqual(alignments.metadata["Program"], "water")
        self.assertEqual(alignments.metadata["Rundate"], "Sat Apr 04 2009 22:08:44")
        self.assertEqual(
            alignments.metadata["Command line"],
            "water -asequence asis:ACACACTCACACACACTTGGTCAGAGATGCTGTGCTTCTTGGAAGCAAGGNCTCAAAGGCAAGGTGCACGCAGAGGGACGTTTGAGTCTGGGATGAAGCATGTNCGTATTATTTATATGATGGAATTTCACGTTTTTATG -bsequence asis:CGTTTGAGTACTGGGATG -gapopen 10 -gapextend 0.5 -filter",
        )
        self.assertEqual(alignments.metadata["Align_format"], "srspair")
        self.assertEqual(alignments.metadata["Report_file"], "stdout")
        alignment = next(alignments)
        self.assertEqual(alignment.annotations["Matrix"], "EDNAFULL")
        self.assertAlmostEqual(alignment.annotations["Gap_penalty"], 10.0)
        self.assertAlmostEqual(alignment.annotations["Extend_penalty"], 0.5)
        self.assertEqual(alignment.annotations["Identity"], 17)
        self.assertEqual(alignment.annotations["Similarity"], 17)
        self.assertEqual(alignment.annotations["Gaps"], 1)
        self.assertAlmostEqual(alignment.annotations["Score"], 75.0)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 18))
        self.assertEqual(alignment.sequences[0].id, "asis")
        self.assertEqual(alignment.sequences[1].id, "asis")
        self.assertEqual(
            repr(alignment.sequences[0].seq),
            "Seq({78: 'CGTTTGAGTCTGGGATG'}, length=95)",
        )
        self.assertEqual(alignment.sequences[0].seq[78:95], "CGTTTGAGTCTGGGATG")
        self.assertEqual(alignment.sequences[1].seq[0:18], "CGTTTGAGTACTGGGATG")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates, np.array([[78, 87, 87, 95], [0, 9, 10, 18]])
            )
        )
        self.assertEqual(alignment[0], "CGTTTGAGT-CTGGGATG")
        self.assertEqual(alignment[1], "CGTTTGAGTACTGGGATG")
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"], "||||||||| ||||||||"
        )
        self.assertEqual(
            str(alignment),
            """\
asis             78 CGTTTGAGT-CTGGGATG 95
                  0 |||||||||-|||||||| 18
asis              0 CGTTTGAGTACTGGGATG 18
""",
        )
        self.assertTrue(
            np.array_equal(
                np.array(alignment, "U"),
                # fmt: off
np.array([['C', 'G', 'T', 'T', 'T', 'G', 'A', 'G', 'T', '-', 'C', 'T', 'G',
           'G', 'G', 'A', 'T', 'G'],
          ['C', 'G', 'T', 'T', 'T', 'G', 'A', 'G', 'T', 'A', 'C', 'T', 'G',
           'G', 'G', 'A', 'T', 'G']], dtype='U')
                # fmt: on
            )
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_matcher_simple(self):
        path = "Emboss/matcher_simple.txt"
        alignments = Align.parse(path, "emboss")
        self.assertEqual(alignments.metadata["Program"], "matcher")
        self.assertEqual(alignments.metadata["Rundate"], "Tue  8 Dec 2009 11:48:35")
        self.assertEqual(
            alignments.metadata["Command line"],
            "matcher [-asequence] rose.pro [-bsequence] rosemary.pro [-outfile] matcher_simple.txt -auto -sprotein -aformat simple",
        )
        self.assertEqual(alignments.metadata["Align_format"], "simple")
        self.assertEqual(alignments.metadata["Report_file"], "matcher_simple.txt")
        alignment = next(alignments)
        self.assertEqual(alignment.annotations["Matrix"], "EBLOSUM62")
        self.assertAlmostEqual(alignment.annotations["Gap_penalty"], 14)
        self.assertAlmostEqual(alignment.annotations["Extend_penalty"], 4)
        self.assertEqual(alignment.annotations["Identity"], 7)
        self.assertEqual(alignment.annotations["Similarity"], 8)
        self.assertEqual(alignment.annotations["Gaps"], 0)
        self.assertAlmostEqual(alignment.annotations["Score"], 29)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 16))
        self.assertEqual(alignment.sequences[0].id, "AF069992_1")
        self.assertEqual(alignment.sequences[1].id, "CAA85685.1")
        self.assertEqual(
            repr(alignment.sequences[0].seq),
            "Seq({72: 'GPPPQSPDENRAGESS'}, length=88)",
        )
        self.assertEqual(
            repr(alignment.sequences[1].seq),
            "Seq({46: 'GVPPEEAGAAVAAESS'}, length=62)",
        )
        self.assertEqual(alignment.sequences[0].seq[72:88], "GPPPQSPDENRAGESS")
        self.assertEqual(alignment.sequences[1].seq[46:62], "GVPPEEAGAAVAAESS")
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[72, 88], [46, 62]]))
        )
        self.assertEqual(alignment[0], "GPPPQSPDENRAGESS")
        self.assertEqual(alignment[1], "GVPPEEAGAAVAAESS")
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"], "|.||:......|.|||"
        )
        self.assertEqual(
            str(alignment),
            """\
AF069992_        72 GPPPQSPDENRAGESS 88
                  0 |.||.......|.||| 16
CAA85685.        46 GVPPEEAGAAVAAESS 62
""",
        )
        self.assertTrue(
            np.array_equal(
                np.array(alignment, "U"),
                # fmt: off
np.array([['G', 'P', 'P', 'P', 'Q', 'S', 'P', 'D', 'E', 'N', 'R', 'A', 'G',
           'E', 'S', 'S'],
          ['G', 'V', 'P', 'P', 'E', 'E', 'A', 'G', 'A', 'A', 'V', 'A', 'A',
           'E', 'S', 'S']], dtype='U')
                # fmt: on
            )
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_matcher_pair(self):
        path = "Emboss/matcher_pair.txt"
        alignments = Align.parse(path, "emboss")
        self.check_matcher_pair(alignments)
        alignments.rewind()
        self.check_matcher_pair(alignments)
        with Align.parse(path, "emboss") as alignments:
            self.check_matcher_pair(alignments)
        with self.assertRaises(AttributeError):
            alignments._stream
        with Align.parse(path, "emboss") as alignments:
            pass
        with self.assertRaises(AttributeError):
            alignments._stream

    def check_matcher_pair(self, alignments):
        self.assertEqual(alignments.metadata["Program"], "matcher")
        self.assertEqual(alignments.metadata["Rundate"], "Tue  8 Dec 2009 12:01:34")
        self.assertEqual(
            alignments.metadata["Command line"],
            "matcher [-asequence] hba_human.fasta [-bsequence] hbb_human.fasta [-outfile] matcher_pair.txt -alternatives 5 -aformat pair -sprotein",
        )
        self.assertEqual(alignments.metadata["Align_format"], "pair")
        self.assertEqual(alignments.metadata["Report_file"], "matcher_pair.txt")
        alignment = next(alignments)
        self.assertEqual(alignment.annotations["Matrix"], "EBLOSUM62")
        self.assertAlmostEqual(alignment.annotations["Gap_penalty"], 14)
        self.assertAlmostEqual(alignment.annotations["Extend_penalty"], 4)
        self.assertEqual(alignment.annotations["Identity"], 63)
        self.assertEqual(alignment.annotations["Similarity"], 88)
        self.assertEqual(alignment.annotations["Gaps"], 8)
        self.assertAlmostEqual(alignment.annotations["Score"], 264)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 145))
        self.assertEqual(alignment.sequences[0].id, "HBA_HUMAN")
        self.assertEqual(alignment.sequences[1].id, "HBB_HUMAN")
        self.assertEqual(
            repr(alignment.sequences[0].seq),
            "Seq({2: 'LSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQV...SKY'}, length=141)",
        )
        self.assertEqual(
            repr(alignment.sequences[1].seq),
            "Seq({3: 'LTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMG...HKY'}, length=146)",
        )
        self.assertEqual(
            alignment.sequences[0].seq[2:141],
            "LSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKY",
        )
        self.assertEqual(
            alignment.sequences[1].seq[3:146],
            "LTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKY",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                np.array(
                    [[2, 18, 20, 47, 47, 51, 51, 141], [3, 19, 19, 46, 47, 51, 56, 146]]
                ),
            )
        )
        self.assertEqual(
            alignment[0],
            "LSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF-DLSH-----GSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKY",
        )
        self.assertEqual(
            alignment[1],
            "LTPEEKSAVTALWGKV--NVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKY",
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            "|:|.:|:.|.|.||||  :..|.|.|||.|:.:.:|.|:.:|..| |||.     |:.:||.|||||..|.::.:||:|::....:.||:||..||.|||.||:||.:.|:..||.|...||||.|.|:..|.:|.|:..|..||",
        )
        self.assertEqual(
            str(alignment),
            """\
HBA_HUMAN         2 LSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF-DLSH-----GSAQV
                  0 |.|..|..|.|.||||--...|.|.|||.|.....|.|...|..|-|||.-----|...|
HBB_HUMAN         3 LTPEEKSAVTALWGKV--NVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKV

HBA_HUMAN        56 KGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPA
                 60 |.|||||..|.....||.|........||.||..||.|||.||.||...|...||.|...
HBB_HUMAN        61 KAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGK

HBA_HUMAN       116 EFTPAVHASLDKFLASVSTVLTSKY 141
                120 ||||.|.|...|..|.|...|..|| 145
HBB_HUMAN       121 EFTPPVQAAYQKVVAGVANALAHKY 146
""",
        )
        self.assertTrue(
            np.array_equal(
                np.array(alignment, "U"),
                # fmt: off
np.array([['L', 'S', 'P', 'A', 'D', 'K', 'T', 'N', 'V', 'K', 'A', 'A', 'W',
           'G', 'K', 'V', 'G', 'A', 'H', 'A', 'G', 'E', 'Y', 'G', 'A', 'E',
           'A', 'L', 'E', 'R', 'M', 'F', 'L', 'S', 'F', 'P', 'T', 'T', 'K',
           'T', 'Y', 'F', 'P', 'H', 'F', '-', 'D', 'L', 'S', 'H', '-', '-',
           '-', '-', '-', 'G', 'S', 'A', 'Q', 'V', 'K', 'G', 'H', 'G', 'K',
           'K', 'V', 'A', 'D', 'A', 'L', 'T', 'N', 'A', 'V', 'A', 'H', 'V',
           'D', 'D', 'M', 'P', 'N', 'A', 'L', 'S', 'A', 'L', 'S', 'D', 'L',
           'H', 'A', 'H', 'K', 'L', 'R', 'V', 'D', 'P', 'V', 'N', 'F', 'K',
           'L', 'L', 'S', 'H', 'C', 'L', 'L', 'V', 'T', 'L', 'A', 'A', 'H',
           'L', 'P', 'A', 'E', 'F', 'T', 'P', 'A', 'V', 'H', 'A', 'S', 'L',
           'D', 'K', 'F', 'L', 'A', 'S', 'V', 'S', 'T', 'V', 'L', 'T', 'S',
           'K', 'Y'],
          ['L', 'T', 'P', 'E', 'E', 'K', 'S', 'A', 'V', 'T', 'A', 'L', 'W',
           'G', 'K', 'V', '-', '-', 'N', 'V', 'D', 'E', 'V', 'G', 'G', 'E',
           'A', 'L', 'G', 'R', 'L', 'L', 'V', 'V', 'Y', 'P', 'W', 'T', 'Q',
           'R', 'F', 'F', 'E', 'S', 'F', 'G', 'D', 'L', 'S', 'T', 'P', 'D',
           'A', 'V', 'M', 'G', 'N', 'P', 'K', 'V', 'K', 'A', 'H', 'G', 'K',
           'K', 'V', 'L', 'G', 'A', 'F', 'S', 'D', 'G', 'L', 'A', 'H', 'L',
           'D', 'N', 'L', 'K', 'G', 'T', 'F', 'A', 'T', 'L', 'S', 'E', 'L',
           'H', 'C', 'D', 'K', 'L', 'H', 'V', 'D', 'P', 'E', 'N', 'F', 'R',
           'L', 'L', 'G', 'N', 'V', 'L', 'V', 'C', 'V', 'L', 'A', 'H', 'H',
           'F', 'G', 'K', 'E', 'F', 'T', 'P', 'P', 'V', 'Q', 'A', 'A', 'Y',
           'Q', 'K', 'V', 'V', 'A', 'G', 'V', 'A', 'N', 'A', 'L', 'A', 'H',
           'K', 'Y']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.annotations["Matrix"], "EBLOSUM62")
        self.assertAlmostEqual(alignment.annotations["Gap_penalty"], 14)
        self.assertAlmostEqual(alignment.annotations["Extend_penalty"], 4)
        self.assertEqual(alignment.annotations["Identity"], 6)
        self.assertEqual(alignment.annotations["Similarity"], 9)
        self.assertEqual(alignment.annotations["Gaps"], 0)
        self.assertAlmostEqual(alignment.annotations["Score"], 32)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 13))
        self.assertEqual(alignment.sequences[0].id, "HBA_HUMAN")
        self.assertEqual(alignment.sequences[1].id, "HBB_HUMAN")
        self.assertEqual(
            repr(alignment.sequences[0].seq),
            "Seq({60: 'KKVADALTNAVAH'}, length=73)",
        )
        self.assertEqual(
            repr(alignment.sequences[1].seq),
            "Seq({131: 'QKVVAGVANALAH'}, length=144)",
        )
        self.assertEqual(alignment.sequences[0].seq[60:73], "KKVADALTNAVAH")
        self.assertEqual(alignment.sequences[1].seq[131:144], "QKVVAGVANALAH")
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[60, 73], [131, 144]]))
        )
        self.assertEqual(alignment[0], "KKVADALTNAVAH")
        self.assertEqual(alignment[1], "QKVVAGVANALAH")
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"], ":||...:.||:||"
        )
        self.assertEqual(
            str(alignment),
            """\
HBA_HUMAN        60 KKVADALTNAVAH  73
                  0 .||.....||.||  13
HBB_HUMAN       131 QKVVAGVANALAH 144
""",
        )
        self.assertTrue(
            np.array_equal(
                np.array(alignment, "U"),
                # fmt: off
np.array([['K', 'K', 'V', 'A', 'D', 'A', 'L', 'T', 'N', 'A', 'V', 'A', 'H'],
          ['Q', 'K', 'V', 'V', 'A', 'G', 'V', 'A', 'N', 'A', 'L', 'A', 'H']],
         dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.annotations["Matrix"], "EBLOSUM62")
        self.assertAlmostEqual(alignment.annotations["Gap_penalty"], 14)
        self.assertAlmostEqual(alignment.annotations["Extend_penalty"], 4)
        self.assertEqual(alignment.annotations["Identity"], 7)
        self.assertEqual(alignment.annotations["Similarity"], 10)
        self.assertEqual(alignment.annotations["Gaps"], 0)
        self.assertAlmostEqual(alignment.annotations["Score"], 28)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 18))
        self.assertEqual(alignment.sequences[0].id, "HBA_HUMAN")
        self.assertEqual(alignment.sequences[1].id, "HBB_HUMAN")
        self.assertEqual(
            repr(alignment.sequences[0].seq),
            "Seq({90: 'KLRVDPVNFKLLSHCLLV'}, length=108)",
        )
        self.assertEqual(
            repr(alignment.sequences[1].seq),
            "Seq({17: 'KVNVDEVGGEALGRLLVV'}, length=35)",
        )
        self.assertEqual(alignment.sequences[0].seq[90:108], "KLRVDPVNFKLLSHCLLV")
        self.assertEqual(alignment.sequences[1].seq[17:35], "KVNVDEVGGEALGRLLVV")
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[90, 108], [17, 35]]))
        )
        self.assertEqual(alignment[0], "KLRVDPVNFKLLSHCLLV")
        self.assertEqual(alignment[1], "KVNVDEVGGEALGRLLVV")
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"], "|:.||.|..:.|...|:|"
        )
        self.assertEqual(
            str(alignment),
            """\
HBA_HUMAN        90 KLRVDPVNFKLLSHCLLV 108
                  0 |..||.|....|...|.|  18
HBB_HUMAN        17 KVNVDEVGGEALGRLLVV  35
""",
        )
        self.assertTrue(
            np.array_equal(
                np.array(alignment, "U"),
                # fmt: off
np.array([['K', 'L', 'R', 'V', 'D', 'P', 'V', 'N', 'F', 'K', 'L', 'L', 'S',
           'H', 'C', 'L', 'L', 'V'],
          ['K', 'V', 'N', 'V', 'D', 'E', 'V', 'G', 'G', 'E', 'A', 'L', 'G',
           'R', 'L', 'L', 'V', 'V']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.annotations["Matrix"], "EBLOSUM62")
        self.assertAlmostEqual(alignment.annotations["Gap_penalty"], 14)
        self.assertAlmostEqual(alignment.annotations["Extend_penalty"], 4)
        self.assertEqual(alignment.annotations["Identity"], 6)
        self.assertEqual(alignment.annotations["Similarity"], 6)
        self.assertEqual(alignment.annotations["Gaps"], 0)
        self.assertAlmostEqual(alignment.annotations["Score"], 23)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 10))
        self.assertEqual(alignment.sequences[0].id, "HBA_HUMAN")
        self.assertEqual(alignment.sequences[1].id, "HBB_HUMAN")
        self.assertEqual(
            repr(alignment.sequences[0].seq),
            "Seq({80: 'LSALSDLHAH'}, length=90)",
        )
        self.assertEqual(
            repr(alignment.sequences[1].seq),
            "Seq({68: 'LGAFSDGLAH'}, length=78)",
        )
        self.assertEqual(alignment.sequences[0].seq[80:90], "LSALSDLHAH")
        self.assertEqual(alignment.sequences[1].seq[68:78], "LGAFSDGLAH")
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[80, 90], [68, 78]]))
        )
        self.assertEqual(alignment[0], "LSALSDLHAH")
        self.assertEqual(alignment[1], "LGAFSDGLAH")
        self.assertEqual(alignment.column_annotations["emboss_consensus"], "|.|.||..||")
        self.assertEqual(
            str(alignment),
            """\
HBA_HUMAN        80 LSALSDLHAH 90
                  0 |.|.||..|| 10
HBB_HUMAN        68 LGAFSDGLAH 78
""",
        )
        self.assertTrue(
            np.array_equal(
                np.array(alignment, "U"),
                # fmt: off
np.array([['L', 'S', 'A', 'L', 'S', 'D', 'L', 'H', 'A', 'H'],
          ['L', 'G', 'A', 'F', 'S', 'D', 'G', 'L', 'A', 'H']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.annotations["Matrix"], "EBLOSUM62")
        self.assertAlmostEqual(alignment.annotations["Gap_penalty"], 14)
        self.assertAlmostEqual(alignment.annotations["Extend_penalty"], 4)
        self.assertEqual(alignment.annotations["Identity"], 6)
        self.assertEqual(alignment.annotations["Similarity"], 8)
        self.assertEqual(alignment.annotations["Gaps"], 0)
        self.assertAlmostEqual(alignment.annotations["Score"], 23)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 10))
        self.assertEqual(alignment.sequences[0].id, "HBA_HUMAN")
        self.assertEqual(alignment.sequences[1].id, "HBB_HUMAN")
        self.assertEqual(
            repr(alignment.sequences[0].seq),
            "Seq({10: 'VKAAWGKVGA'}, length=20)",
        )
        self.assertEqual(
            repr(alignment.sequences[1].seq),
            "Seq({126: 'VQAAYQKVVA'}, length=136)",
        )
        self.assertEqual(alignment.sequences[0].seq[10:20], "VKAAWGKVGA")
        self.assertEqual(alignment.sequences[1].seq[126:136], "VQAAYQKVVA")
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[10, 20], [126, 136]]))
        )
        self.assertEqual(alignment[0], "VKAAWGKVGA")
        self.assertEqual(alignment[1], "VQAAYQKVVA")
        self.assertEqual(alignment.column_annotations["emboss_consensus"], "|:||:.||.|")
        self.assertEqual(
            str(alignment),
            """\
HBA_HUMAN        10 VKAAWGKVGA  20
                  0 |.||..||.|  10
HBB_HUMAN       126 VQAAYQKVVA 136
""",
        )
        self.assertTrue(
            np.array_equal(
                np.array(alignment, "U"),
                # fmt: off
np.array([['V', 'K', 'A', 'A', 'W', 'G', 'K', 'V', 'G', 'A'],
          ['V', 'Q', 'A', 'A', 'Y', 'Q', 'K', 'V', 'V', 'A']], dtype='U')
                # fmt: on
            )
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_pair_example_nobrief(self):
        # Variation on the alignment file obtained from EMBOSS
        # (http://emboss.sourceforge.net/docs/themes/alnformats/align.pair)
        # if we include 3 sequences to align against, and we use the -nobrief
        # command line option.
        path = "Emboss/needle_nobrief_multiple.pair"
        alignments = Align.parse(path, "emboss")
        self.assertEqual(alignments.metadata["Program"], "needle")
        self.assertEqual(alignments.metadata["Rundate"], "Fri 23 Jul 2021 22:45:41")
        self.assertEqual(
            alignments.metadata["Command line"],
            "needle -asequence seqa.fa -bsequence seqb.fa -datafile EBLOSUM62 -gapopen 10 -gapextend 0.5 -nobrief -outfile stdout",
        )
        self.assertEqual(alignments.metadata["Align_format"], "srspair")
        self.assertEqual(alignments.metadata["Report_file"], "stdout")
        alignment = next(alignments)
        self.assertEqual(alignment.annotations["Matrix"], "EBLOSUM62")
        self.assertAlmostEqual(alignment.annotations["Gap_penalty"], 10.0)
        self.assertAlmostEqual(alignment.annotations["Extend_penalty"], 0.5)
        self.assertEqual(alignment.annotations["Identity"], 112)
        self.assertEqual(alignment.annotations["Similarity"], 112)
        self.assertEqual(alignment.annotations["Gaps"], 19)
        self.assertAlmostEqual(alignment.annotations["Score"], 591.5)
        self.assertEqual(alignment.annotations["Longest_Identity"], "100.00%")
        self.assertEqual(alignment.annotations["Longest_Similarity"], "100.00%")
        self.assertEqual(alignment.annotations["Shortest_Identity"], "85.50%")
        self.assertEqual(alignment.annotations["Shortest_Similarity"], "85.50%")
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 131))
        self.assertEqual(alignment.sequences[0].id, "IXI_234")
        self.assertEqual(alignment.sequences[1].id, "IXI_235")
        self.assertEqual(
            alignment.sequences[0].seq,
            "TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "TSPASIRPPAGPSSRRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                np.array([[0, 15, 24, 74, 84, 131], [0, 15, 15, 65, 65, 112]]),
            )
        )
        self.assertEqual(
            alignment[0],
            "TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertEqual(
            alignment[1],
            "TSPASIRPPAGPSSR---------RPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGW----------RASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            "|||||||||||||||         ||||||||||||||||||||||||||||||||||||||||||||||||||          |||||||||||||||||||||||||||||||||||||||||||||||",
        )
        self.assertEqual(
            str(alignment),
            """\
IXI_234           0 TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTC
                  0 |||||||||||||||---------||||||||||||||||||||||||||||||||||||
IXI_235           0 TSPASIRPPAGPSSR---------RPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTC

IXI_234          60 TTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTG
                 60 ||||||||||||||----------||||||||||||||||||||||||||||||||||||
IXI_235          51 TTSTSTRHRGRSGW----------RASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTG

IXI_234         120 PPAWAGDRSHE 131
                120 ||||||||||| 131
IXI_235         101 PPAWAGDRSHE 112
""",
        )
        self.assertTrue(
            np.array_equal(
                np.array(alignment, "U"),
                # fmt: off
np.array([['T', 'S', 'P', 'A', 'S', 'I', 'R', 'P', 'P', 'A', 'G', 'P', 'S',
           'S', 'R', 'P', 'A', 'M', 'V', 'S', 'S', 'R', 'R', 'T', 'R', 'P',
           'S', 'P', 'P', 'G', 'P', 'R', 'R', 'P', 'T', 'G', 'R', 'P', 'C',
           'C', 'S', 'A', 'A', 'P', 'R', 'R', 'P', 'Q', 'A', 'T', 'G', 'G',
           'W', 'K', 'T', 'C', 'S', 'G', 'T', 'C', 'T', 'T', 'S', 'T', 'S',
           'T', 'R', 'H', 'R', 'G', 'R', 'S', 'G', 'W', 'S', 'A', 'R', 'T',
           'T', 'T', 'A', 'A', 'C', 'L', 'R', 'A', 'S', 'R', 'K', 'S', 'M',
           'R', 'A', 'A', 'C', 'S', 'R', 'S', 'A', 'G', 'S', 'R', 'P', 'N',
           'R', 'F', 'A', 'P', 'T', 'L', 'M', 'S', 'S', 'C', 'I', 'T', 'S',
           'T', 'T', 'G', 'P', 'P', 'A', 'W', 'A', 'G', 'D', 'R', 'S', 'H',
           'E'],
          ['T', 'S', 'P', 'A', 'S', 'I', 'R', 'P', 'P', 'A', 'G', 'P', 'S',
           'S', 'R', '-', '-', '-', '-', '-', '-', '-', '-', '-', 'R', 'P',
           'S', 'P', 'P', 'G', 'P', 'R', 'R', 'P', 'T', 'G', 'R', 'P', 'C',
           'C', 'S', 'A', 'A', 'P', 'R', 'R', 'P', 'Q', 'A', 'T', 'G', 'G',
           'W', 'K', 'T', 'C', 'S', 'G', 'T', 'C', 'T', 'T', 'S', 'T', 'S',
           'T', 'R', 'H', 'R', 'G', 'R', 'S', 'G', 'W', '-', '-', '-', '-',
           '-', '-', '-', '-', '-', '-', 'R', 'A', 'S', 'R', 'K', 'S', 'M',
           'R', 'A', 'A', 'C', 'S', 'R', 'S', 'A', 'G', 'S', 'R', 'P', 'N',
           'R', 'F', 'A', 'P', 'T', 'L', 'M', 'S', 'S', 'C', 'I', 'T', 'S',
           'T', 'T', 'G', 'P', 'P', 'A', 'W', 'A', 'G', 'D', 'R', 'S', 'H',
           'E']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.annotations["Matrix"], "EBLOSUM62")
        self.assertAlmostEqual(alignment.annotations["Gap_penalty"], 10.0)
        self.assertAlmostEqual(alignment.annotations["Extend_penalty"], 0.5)
        self.assertEqual(alignment.annotations["Identity"], 120)
        self.assertEqual(alignment.annotations["Similarity"], 120)
        self.assertEqual(alignment.annotations["Gaps"], 4)
        self.assertAlmostEqual(alignment.annotations["Score"], 618.0)
        self.assertEqual(alignment.annotations["Longest_Identity"], "94.49%")
        self.assertEqual(alignment.annotations["Longest_Similarity"], "94.49%")
        self.assertEqual(alignment.annotations["Shortest_Identity"], "91.60%")
        self.assertEqual(alignment.annotations["Shortest_Similarity"], "91.60%")
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 131))
        self.assertEqual(alignment.sequences[0].id, "IXI_234")
        self.assertEqual(alignment.sequences[1].id, "IXI_236")
        self.assertEqual(
            alignment.sequences[0].seq,
            "TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "TSPASIRPPAGPSSRPAMVSSRRPSPPPPRRPPGRPCCSAAPPRPQATGGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRGSRPPRFAPPLMSSCITSTTGPPPPAGDRSHE",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                np.array([[0, 22, 24, 97, 99, 131], [0, 22, 22, 95, 95, 127]]),
            )
        )
        self.assertEqual(
            alignment[0],
            "TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertEqual(
            alignment[1],
            "TSPASIRPPAGPSSRPAMVSSR--RPSPPPPRRPPGRPCCSAAPPRPQATGGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSR--GSRPPRFAPPLMSSCITSTTGPPPPAGDRSHE",
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            "||||||||||||||||||||||  |||||.||||.|||||||||.||||||||||||||||||||||||||||||||||||||||||||||||||||  ||||.||||.|||||||||||||..|||||||",
        )
        self.assertEqual(
            str(alignment),
            """\
IXI_234           0 TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTC
                  0 ||||||||||||||||||||||--|||||.||||.|||||||||.|||||||||||||||
IXI_236           0 TSPASIRPPAGPSSRPAMVSSR--RPSPPPPRRPPGRPCCSAAPPRPQATGGWKTCSGTC

IXI_234          60 TTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTG
                 60 |||||||||||||||||||||||||||||||||||||--||||.||||.|||||||||||
IXI_236          58 TTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSR--GSRPPRFAPPLMSSCITSTTG

IXI_234         120 PPAWAGDRSHE 131
                120 ||..||||||| 131
IXI_236         116 PPPPAGDRSHE 127
""",
        )
        self.assertTrue(
            np.array_equal(
                np.array(alignment, "U"),
                # fmt: off
np.array([['T', 'S', 'P', 'A', 'S', 'I', 'R', 'P', 'P', 'A', 'G', 'P', 'S',
           'S', 'R', 'P', 'A', 'M', 'V', 'S', 'S', 'R', 'R', 'T', 'R', 'P',
           'S', 'P', 'P', 'G', 'P', 'R', 'R', 'P', 'T', 'G', 'R', 'P', 'C',
           'C', 'S', 'A', 'A', 'P', 'R', 'R', 'P', 'Q', 'A', 'T', 'G', 'G',
           'W', 'K', 'T', 'C', 'S', 'G', 'T', 'C', 'T', 'T', 'S', 'T', 'S',
           'T', 'R', 'H', 'R', 'G', 'R', 'S', 'G', 'W', 'S', 'A', 'R', 'T',
           'T', 'T', 'A', 'A', 'C', 'L', 'R', 'A', 'S', 'R', 'K', 'S', 'M',
           'R', 'A', 'A', 'C', 'S', 'R', 'S', 'A', 'G', 'S', 'R', 'P', 'N',
           'R', 'F', 'A', 'P', 'T', 'L', 'M', 'S', 'S', 'C', 'I', 'T', 'S',
           'T', 'T', 'G', 'P', 'P', 'A', 'W', 'A', 'G', 'D', 'R', 'S', 'H',
           'E'],
          ['T', 'S', 'P', 'A', 'S', 'I', 'R', 'P', 'P', 'A', 'G', 'P', 'S',
           'S', 'R', 'P', 'A', 'M', 'V', 'S', 'S', 'R', '-', '-', 'R', 'P',
           'S', 'P', 'P', 'P', 'P', 'R', 'R', 'P', 'P', 'G', 'R', 'P', 'C',
           'C', 'S', 'A', 'A', 'P', 'P', 'R', 'P', 'Q', 'A', 'T', 'G', 'G',
           'W', 'K', 'T', 'C', 'S', 'G', 'T', 'C', 'T', 'T', 'S', 'T', 'S',
           'T', 'R', 'H', 'R', 'G', 'R', 'S', 'G', 'W', 'S', 'A', 'R', 'T',
           'T', 'T', 'A', 'A', 'C', 'L', 'R', 'A', 'S', 'R', 'K', 'S', 'M',
           'R', 'A', 'A', 'C', 'S', 'R', '-', '-', 'G', 'S', 'R', 'P', 'P',
           'R', 'F', 'A', 'P', 'P', 'L', 'M', 'S', 'S', 'C', 'I', 'T', 'S',
           'T', 'T', 'G', 'P', 'P', 'P', 'P', 'A', 'G', 'D', 'R', 'S', 'H',
           'E']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.annotations["Matrix"], "EBLOSUM62")
        self.assertAlmostEqual(alignment.annotations["Gap_penalty"], 10.0)
        self.assertAlmostEqual(alignment.annotations["Extend_penalty"], 0.5)
        self.assertEqual(alignment.annotations["Identity"], 119)
        self.assertEqual(alignment.annotations["Similarity"], 124)
        self.assertEqual(alignment.annotations["Gaps"], 7)
        self.assertAlmostEqual(alignment.annotations["Score"], 609.0)
        self.assertEqual(alignment.annotations["Longest_Identity"], "95.97%")
        self.assertEqual(alignment.annotations["Longest_Similarity"], "100.00%")
        self.assertEqual(alignment.annotations["Shortest_Identity"], "90.84%")
        self.assertEqual(alignment.annotations["Shortest_Similarity"], "94.66%")
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 131))
        self.assertEqual(alignment.sequences[0].id, "IXI_234")
        self.assertEqual(alignment.sequences[1].id, "IXI_237")
        self.assertEqual(
            alignment.sequences[0].seq,
            "TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "TSPASLRPPAGPSSRPAMVSSRRRPSPPGPRRPTCSAAPRRPQATGGYKTCSGTCTTSTSTRHRGRSGYSARTTTAACLRASRKSMRAACSRGSRPNRFAPTLMSSCLTSTTGPPAYAGDRSHE",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                np.array(
                    [[0, 23, 24, 35, 39, 97, 99, 131], [0, 23, 23, 34, 34, 92, 92, 124]]
                ),
            )
        )
        self.assertEqual(
            alignment[0],
            "TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertEqual(
            alignment[1],
            "TSPASLRPPAGPSSRPAMVSSRR-RPSPPGPRRPT----CSAAPRRPQATGGYKTCSGTCTTSTSTRHRGRSGYSARTTTAACLRASRKSMRAACSR--GSRPNRFAPTLMSSCLTSTTGPPAYAGDRSHE",
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            "|||||:||||||||||||||||| |||||||||||    |||||||||||||:||||||||||||||||||||:|||||||||||||||||||||||  |||||||||||||||:||||||||:|||||||",
        )
        self.assertEqual(
            str(alignment),
            """\
IXI_234           0 TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTC
                  0 |||||.|||||||||||||||||-|||||||||||----|||||||||||||.|||||||
IXI_237           0 TSPASLRPPAGPSSRPAMVSSRR-RPSPPGPRRPT----CSAAPRRPQATGGYKTCSGTC

IXI_234          60 TTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTG
                 60 |||||||||||||.|||||||||||||||||||||||--|||||||||||||||.|||||
IXI_237          55 TTSTSTRHRGRSGYSARTTTAACLRASRKSMRAACSR--GSRPNRFAPTLMSSCLTSTTG

IXI_234         120 PPAWAGDRSHE 131
                120 |||.||||||| 131
IXI_237         113 PPAYAGDRSHE 124
""",
        )
        self.assertTrue(
            np.array_equal(
                np.array(alignment, "U"),
                # fmt: off
np.array([['T', 'S', 'P', 'A', 'S', 'I', 'R', 'P', 'P', 'A', 'G', 'P', 'S',
           'S', 'R', 'P', 'A', 'M', 'V', 'S', 'S', 'R', 'R', 'T', 'R', 'P',
           'S', 'P', 'P', 'G', 'P', 'R', 'R', 'P', 'T', 'G', 'R', 'P', 'C',
           'C', 'S', 'A', 'A', 'P', 'R', 'R', 'P', 'Q', 'A', 'T', 'G', 'G',
           'W', 'K', 'T', 'C', 'S', 'G', 'T', 'C', 'T', 'T', 'S', 'T', 'S',
           'T', 'R', 'H', 'R', 'G', 'R', 'S', 'G', 'W', 'S', 'A', 'R', 'T',
           'T', 'T', 'A', 'A', 'C', 'L', 'R', 'A', 'S', 'R', 'K', 'S', 'M',
           'R', 'A', 'A', 'C', 'S', 'R', 'S', 'A', 'G', 'S', 'R', 'P', 'N',
           'R', 'F', 'A', 'P', 'T', 'L', 'M', 'S', 'S', 'C', 'I', 'T', 'S',
           'T', 'T', 'G', 'P', 'P', 'A', 'W', 'A', 'G', 'D', 'R', 'S', 'H',
           'E'],
          ['T', 'S', 'P', 'A', 'S', 'L', 'R', 'P', 'P', 'A', 'G', 'P', 'S',
           'S', 'R', 'P', 'A', 'M', 'V', 'S', 'S', 'R', 'R', '-', 'R', 'P',
           'S', 'P', 'P', 'G', 'P', 'R', 'R', 'P', 'T', '-', '-', '-', '-',
           'C', 'S', 'A', 'A', 'P', 'R', 'R', 'P', 'Q', 'A', 'T', 'G', 'G',
           'Y', 'K', 'T', 'C', 'S', 'G', 'T', 'C', 'T', 'T', 'S', 'T', 'S',
           'T', 'R', 'H', 'R', 'G', 'R', 'S', 'G', 'Y', 'S', 'A', 'R', 'T',
           'T', 'T', 'A', 'A', 'C', 'L', 'R', 'A', 'S', 'R', 'K', 'S', 'M',
           'R', 'A', 'A', 'C', 'S', 'R', '-', '-', 'G', 'S', 'R', 'P', 'N',
           'R', 'F', 'A', 'P', 'T', 'L', 'M', 'S', 'S', 'C', 'L', 'T', 'S',
           'T', 'T', 'G', 'P', 'P', 'A', 'Y', 'A', 'G', 'D', 'R', 'S', 'H',
           'E']], dtype='U')
                # fmt: on
            )
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_pair_example2(self):
        path = "Emboss/needle.txt"
        alignments = Align.parse(path, "emboss")
        self.assertEqual(alignments.metadata["Program"], "needle")
        self.assertEqual(alignments.metadata["Rundate"], "Sun 27 Apr 2007 17:20:35")
        self.assertEqual(
            alignments.metadata["Command line"],
            "needle [-asequence] Spo0F.faa [-bsequence] paired_r.faa -sformat2 pearson",
        )
        self.assertEqual(alignments.metadata["Align_format"], "srspair")
        self.assertEqual(alignments.metadata["Report_file"], "ref_rec .needle")
        alignment = next(alignments)
        self.assertEqual(alignment.annotations["Matrix"], "EBLOSUM62")
        self.assertAlmostEqual(alignment.annotations["Gap_penalty"], 10.0)
        self.assertAlmostEqual(alignment.annotations["Extend_penalty"], 0.5)
        self.assertEqual(alignment.annotations["Identity"], 32)
        self.assertEqual(alignment.annotations["Similarity"], 64)
        self.assertEqual(alignment.annotations["Gaps"], 17)
        self.assertAlmostEqual(alignment.annotations["Score"], 112.0)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 124))
        self.assertEqual(alignment.sequences[0].id, "ref_rec")
        self.assertEqual(alignment.sequences[1].id, "gi|94968718|receiver")
        self.assertEqual(
            alignment.sequences[0].seq,
            "KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIVTKERPDLVLLDMKIPGMDGIEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHFAKPFDIDEIRDAV",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "VLLADDHALVRRGFRLMLEDDPEIEIVAEAGDGAQAVKLAGELHPRVVVMDCAMPGMSGMDATKQIRTQWPDIAVLMLTMHSEDTWVRLALEAGANGYILKSAIDLDLIQAVRRVANGET",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                np.array(
                    [
                        [0, 1, 7, 7, 17, 19, 100, 100, 108, 109, 111, 111],
                        [0, 0, 6, 10, 20, 20, 101, 102, 110, 110, 112, 120],
                    ]
                ),
            )
        )
        self.assertEqual(
            alignment[0],
            "KILIVDD----QYGIRILLNEVFNKEGYQTFQAANGLQALDIVTKERPDLVLLDMKIPGMDGIEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHFAK-PFDIDEIRDAV--------",
        )
        self.assertEqual(
            alignment[1],
            "-VLLADDHALVRRGFRLMLED--DPEIEIVAEAGDGAQAVKLAGELHPRVVVMDCAMPGMSGMDATKQIRTQWPDIAVLMLTMHSEDTWVRLALEAGANGYILKSAIDLDLIQ-AVRRVANGET",
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            " :|:.||    :.|.|::|.:  :.|.....:|.:|.||:.:..:..|.:|::|..:|||.|::..|:::....:|.|:::|.:.|...::.:.|.||..:..| ..|:|.|: ||        ",
        )
        self.assertEqual(
            str(alignment),
            """\
ref_rec           0 KILIVDD----QYGIRILLNEVFNKEGYQTFQAANGLQALDIVTKERPDLVLLDMKIPGM
                  0 -.|..||----..|.|..|..--..|......|..|.||........|..|..|...|||
gi|949687         0 -VLLADDHALVRRGFRLMLED--DPEIEIVAEAGDGAQAVKLAGELHPRVVVMDCAMPGM

ref_rec          56 DGIEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHFAK-PFDIDEIRDAV----
                 60 .|....|........|.|...|...|........|.||.....|-..|.|.|.-||----
gi|949687        57 SGMDATKQIRTQWPDIAVLMLTMHSEDTWVRLALEAGANGYILKSAIDLDLIQ-AVRRVA

ref_rec         111 ---- 111
                120 ---- 124
gi|949687       116 NGET 120
""",
        )
        self.assertTrue(
            np.array_equal(
                np.array(alignment, "U"),
                # fmt: off
np.array([['K', 'I', 'L', 'I', 'V', 'D', 'D', '-', '-', '-', '-', 'Q', 'Y',
           'G', 'I', 'R', 'I', 'L', 'L', 'N', 'E', 'V', 'F', 'N', 'K', 'E',
           'G', 'Y', 'Q', 'T', 'F', 'Q', 'A', 'A', 'N', 'G', 'L', 'Q', 'A',
           'L', 'D', 'I', 'V', 'T', 'K', 'E', 'R', 'P', 'D', 'L', 'V', 'L',
           'L', 'D', 'M', 'K', 'I', 'P', 'G', 'M', 'D', 'G', 'I', 'E', 'I',
           'L', 'K', 'R', 'M', 'K', 'V', 'I', 'D', 'E', 'N', 'I', 'R', 'V',
           'I', 'I', 'M', 'T', 'A', 'Y', 'G', 'E', 'L', 'D', 'M', 'I', 'Q',
           'E', 'S', 'K', 'E', 'L', 'G', 'A', 'L', 'T', 'H', 'F', 'A', 'K',
           '-', 'P', 'F', 'D', 'I', 'D', 'E', 'I', 'R', 'D', 'A', 'V', '-',
           '-', '-', '-', '-', '-', '-', '-'],
          ['-', 'V', 'L', 'L', 'A', 'D', 'D', 'H', 'A', 'L', 'V', 'R', 'R',
           'G', 'F', 'R', 'L', 'M', 'L', 'E', 'D', '-', '-', 'D', 'P', 'E',
           'I', 'E', 'I', 'V', 'A', 'E', 'A', 'G', 'D', 'G', 'A', 'Q', 'A',
           'V', 'K', 'L', 'A', 'G', 'E', 'L', 'H', 'P', 'R', 'V', 'V', 'V',
           'M', 'D', 'C', 'A', 'M', 'P', 'G', 'M', 'S', 'G', 'M', 'D', 'A',
           'T', 'K', 'Q', 'I', 'R', 'T', 'Q', 'W', 'P', 'D', 'I', 'A', 'V',
           'L', 'M', 'L', 'T', 'M', 'H', 'S', 'E', 'D', 'T', 'W', 'V', 'R',
           'L', 'A', 'L', 'E', 'A', 'G', 'A', 'N', 'G', 'Y', 'I', 'L', 'K',
           'S', 'A', 'I', 'D', 'L', 'D', 'L', 'I', 'Q', '-', 'A', 'V', 'R',
           'R', 'V', 'A', 'N', 'G', 'E', 'T']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.annotations["Matrix"], "EBLOSUM62")
        self.assertAlmostEqual(alignment.annotations["Gap_penalty"], 10.0)
        self.assertAlmostEqual(alignment.annotations["Extend_penalty"], 0.5)
        self.assertEqual(alignment.annotations["Identity"], 34)
        self.assertEqual(alignment.annotations["Similarity"], 58)
        self.assertEqual(alignment.annotations["Gaps"], 9)
        self.assertAlmostEqual(alignment.annotations["Score"], 154.0)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 119))
        self.assertEqual(alignment.sequences[0].id, "ref_rec")
        self.assertEqual(alignment.sequences[1].id, "gi|94968761|receiver")
        self.assertEqual(
            alignment.sequences[0].seq,
            "KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIVTKERPDLVLLDMKIPGMDGIEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHFAKPFDIDEIRDAV",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "ILIVDDEANTLASLSRAFRLAGHEATVCDNAVRALEIAKSKPFDLILSDVVMPGRDGLTLLEDLKTAGVQAPVVMMSGQAHIEMAVKATRLGALDFLEKPLSTDKLLLTVENALKLKR",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates, np.array([[0, 1, 111, 111], [0, 0, 110, 118]])
            )
        )
        self.assertEqual(
            alignment[0],
            "KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIVTKERPDLVLLDMKIPGMDGIEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHFAKPFDIDEIRDAV--------",
        )
        self.assertEqual(
            alignment[1],
            "-ILIVDDEANTLASLSRAFRLAGHEATVCDNAVRALEIAKSKPFDLILSDVVMPGRDGLTLLEDLKTAGVQAPVVMMSGQAHIEMAVKATRLGALDFLEKPLSTDKLLLTVENALKLKR",
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            " ||||||:......|:..|...|::.....|.::||:|...:..||:|.|:.:||.||:.:|:.:|.......|::|:....::|..::..||||....||...|::...|        ",
        )
        self.assertEqual(
            str(alignment),
            """\
ref_rec           0 KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIVTKERPDLVLLDMKIPGMDGIE
                  0 -||||||.......|...|...|.......|...||.|......||.|.|...||.||..
gi|949687         0 -ILIVDDEANTLASLSRAFRLAGHEATVCDNAVRALEIAKSKPFDLILSDVVMPGRDGLT

ref_rec          60 ILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHFAKPFDIDEIRDAV--------
                 60 .|...|.......|..|.......|......||||....||...|.....|--------
gi|949687        59 LLEDLKTAGVQAPVVMMSGQAHIEMAVKATRLGALDFLEKPLSTDKLLLTVENALKLKR

ref_rec         111
                119
gi|949687       118
""",
        )
        self.assertTrue(
            np.array_equal(
                np.array(alignment, "U"),
                # fmt: off
np.array([['K', 'I', 'L', 'I', 'V', 'D', 'D', 'Q', 'Y', 'G', 'I', 'R', 'I',
           'L', 'L', 'N', 'E', 'V', 'F', 'N', 'K', 'E', 'G', 'Y', 'Q', 'T',
           'F', 'Q', 'A', 'A', 'N', 'G', 'L', 'Q', 'A', 'L', 'D', 'I', 'V',
           'T', 'K', 'E', 'R', 'P', 'D', 'L', 'V', 'L', 'L', 'D', 'M', 'K',
           'I', 'P', 'G', 'M', 'D', 'G', 'I', 'E', 'I', 'L', 'K', 'R', 'M',
           'K', 'V', 'I', 'D', 'E', 'N', 'I', 'R', 'V', 'I', 'I', 'M', 'T',
           'A', 'Y', 'G', 'E', 'L', 'D', 'M', 'I', 'Q', 'E', 'S', 'K', 'E',
           'L', 'G', 'A', 'L', 'T', 'H', 'F', 'A', 'K', 'P', 'F', 'D', 'I',
           'D', 'E', 'I', 'R', 'D', 'A', 'V', '-', '-', '-', '-', '-', '-',
           '-', '-'],
          ['-', 'I', 'L', 'I', 'V', 'D', 'D', 'E', 'A', 'N', 'T', 'L', 'A',
           'S', 'L', 'S', 'R', 'A', 'F', 'R', 'L', 'A', 'G', 'H', 'E', 'A',
           'T', 'V', 'C', 'D', 'N', 'A', 'V', 'R', 'A', 'L', 'E', 'I', 'A',
           'K', 'S', 'K', 'P', 'F', 'D', 'L', 'I', 'L', 'S', 'D', 'V', 'V',
           'M', 'P', 'G', 'R', 'D', 'G', 'L', 'T', 'L', 'L', 'E', 'D', 'L',
           'K', 'T', 'A', 'G', 'V', 'Q', 'A', 'P', 'V', 'V', 'M', 'M', 'S',
           'G', 'Q', 'A', 'H', 'I', 'E', 'M', 'A', 'V', 'K', 'A', 'T', 'R',
           'L', 'G', 'A', 'L', 'D', 'F', 'L', 'E', 'K', 'P', 'L', 'S', 'T',
           'D', 'K', 'L', 'L', 'L', 'T', 'V', 'E', 'N', 'A', 'L', 'K', 'L',
           'K', 'R']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.annotations["Matrix"], "EBLOSUM62")
        self.assertAlmostEqual(alignment.annotations["Gap_penalty"], 10.0)
        self.assertAlmostEqual(alignment.annotations["Extend_penalty"], 0.5)
        self.assertEqual(alignment.annotations["Identity"], 29)
        self.assertEqual(alignment.annotations["Similarity"], 53)
        self.assertEqual(alignment.annotations["Gaps"], 9)
        self.assertAlmostEqual(alignment.annotations["Score"], 121.0)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 120))
        self.assertEqual(alignment.sequences[0].id, "ref_rec")
        self.assertEqual(alignment.sequences[1].id, "gi|94967506|receiver")
        self.assertEqual(
            alignment.sequences[0].seq,
            "KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIVTKERPDLVLLDMKIPGMDGIEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHFAKPFDIDEIRDAV",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "LHIVVVDDDPGTCVYIESVFAELGHTCKSFVRPEAAEEYILTHPVDLAIVDVYLGSTTGVEVLRRCRVHRPKLYAVIITGQISLEMAARSIAEGAVDYIQKPIDIDALLNIAERALEHKE",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates, np.array([[0, 0, 111, 111], [0, 1, 112, 120]])
            )
        )
        self.assertEqual(
            alignment[0],
            "-KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIVTKERPDLVLLDMKIPGMDGIEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHFAKPFDIDEIRDAV--------",
        )
        self.assertEqual(
            alignment[1],
            "LHIVVVDDDPGTCVYIESVFAELGHTCKSFVRPEAAEEYILTHPVDLAIVDVYLGSTTGVEVLRRCRVHRPKLYAVIITGQISLEMAARSIAEGAVDYIQKPIDIDALLNIAERALEHKE",
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            " .|::|||..|..:.:..||.:.|:..........|.:.:.....||.::|:.:....|:|:|:|.:|....:..:|:|....|:|...|...||:.:..||.|||.:.:..        ",
        )
        self.assertEqual(
            str(alignment),
            """\
ref_rec           0 -KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIVTKERPDLVLLDMKIPGMDGI
                  0 -.|..|||..|.......||...|...........|.........||...|.......|.
gi|949675         0 LHIVVVDDDPGTCVYIESVFAELGHTCKSFVRPEAAEEYILTHPVDLAIVDVYLGSTTGV

ref_rec          59 EILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHFAKPFDIDEIRDAV--------
                 60 |.|.|..|........|.|....|.|...|...||.....||.|||......--------
gi|949675        60 EVLRRCRVHRPKLYAVIITGQISLEMAARSIAEGAVDYIQKPIDIDALLNIAERALEHKE

ref_rec         111 
                120 
gi|949675       120 
""",
        )
        self.assertTrue(
            np.array_equal(
                np.array(alignment, "U"),
                # fmt: off
np.array([['-', 'K', 'I', 'L', 'I', 'V', 'D', 'D', 'Q', 'Y', 'G', 'I', 'R',
           'I', 'L', 'L', 'N', 'E', 'V', 'F', 'N', 'K', 'E', 'G', 'Y', 'Q',
           'T', 'F', 'Q', 'A', 'A', 'N', 'G', 'L', 'Q', 'A', 'L', 'D', 'I',
           'V', 'T', 'K', 'E', 'R', 'P', 'D', 'L', 'V', 'L', 'L', 'D', 'M',
           'K', 'I', 'P', 'G', 'M', 'D', 'G', 'I', 'E', 'I', 'L', 'K', 'R',
           'M', 'K', 'V', 'I', 'D', 'E', 'N', 'I', 'R', 'V', 'I', 'I', 'M',
           'T', 'A', 'Y', 'G', 'E', 'L', 'D', 'M', 'I', 'Q', 'E', 'S', 'K',
           'E', 'L', 'G', 'A', 'L', 'T', 'H', 'F', 'A', 'K', 'P', 'F', 'D',
           'I', 'D', 'E', 'I', 'R', 'D', 'A', 'V', '-', '-', '-', '-', '-',
           '-', '-', '-'],
          ['L', 'H', 'I', 'V', 'V', 'V', 'D', 'D', 'D', 'P', 'G', 'T', 'C',
           'V', 'Y', 'I', 'E', 'S', 'V', 'F', 'A', 'E', 'L', 'G', 'H', 'T',
           'C', 'K', 'S', 'F', 'V', 'R', 'P', 'E', 'A', 'A', 'E', 'E', 'Y',
           'I', 'L', 'T', 'H', 'P', 'V', 'D', 'L', 'A', 'I', 'V', 'D', 'V',
           'Y', 'L', 'G', 'S', 'T', 'T', 'G', 'V', 'E', 'V', 'L', 'R', 'R',
           'C', 'R', 'V', 'H', 'R', 'P', 'K', 'L', 'Y', 'A', 'V', 'I', 'I',
           'T', 'G', 'Q', 'I', 'S', 'L', 'E', 'M', 'A', 'A', 'R', 'S', 'I',
           'A', 'E', 'G', 'A', 'V', 'D', 'Y', 'I', 'Q', 'K', 'P', 'I', 'D',
           'I', 'D', 'A', 'L', 'L', 'N', 'I', 'A', 'E', 'R', 'A', 'L', 'E',
           'H', 'K', 'E']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.annotations["Matrix"], "EBLOSUM62")
        self.assertAlmostEqual(alignment.annotations["Gap_penalty"], 10.0)
        self.assertAlmostEqual(alignment.annotations["Extend_penalty"], 0.5)
        self.assertEqual(alignment.annotations["Identity"], 30)
        self.assertEqual(alignment.annotations["Similarity"], 64)
        self.assertEqual(alignment.annotations["Gaps"], 9)
        self.assertAlmostEqual(alignment.annotations["Score"], 126.0)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 118))
        self.assertEqual(alignment.sequences[0].id, "ref_rec")
        self.assertEqual(alignment.sequences[1].id, "gi|94970045|receiver")
        self.assertEqual(
            alignment.sequences[0].seq,
            "KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIVTKERPDLVLLDMKIPGMDGIEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHFAKPFDIDEIRDAV",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "VLLVEDEEALRAAAGDFLETRGYKIMTARDGTEALSMASKFAERIDVLITDLVMPGISGRVLAQELVKIHPETKVMYMSGYDDETVMVNGEIDSSSAFLRKPFRMDALSAKIREVL",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                np.array(
                    [
                        [0, 1, 41, 41, 82, 83, 98, 98, 105, 105, 111],
                        [0, 0, 40, 42, 83, 83, 98, 99, 106, 110, 116],
                    ]
                ),
            )
        )
        self.assertEqual(
            alignment[0],
            "KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIVTK--ERPDLVLLDMKIPGMDGIEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHF-AKPFDID----EIRDAV",
        )
        self.assertEqual(
            alignment[1],
            "-VLLVEDEEALRAAAGDFLETRGYKIMTARDGTEALSMASKFAERIDVLITDLVMPGISGRVLAQELVKIHPETKVMYMSGYDD-ETVMVNGEIDSSSAFLRKPFRMDALSAKIREVL",
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            " :|:|:|:..:|....:.....||:...|.:|.:||.:.:|  ||.|:::.|:.:||:.|..:.:.:..|....:|:.|:.|.: :.:..:.|:.:.:.| .|||.:|    :||:.:",
        )
        self.assertEqual(
            str(alignment),
            """\
ref_rec           0 KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIVTK--ERPDLVLLDMKIPGMDG
                  0 -.|.|.|....|..........||....|..|..||....|--||.|....|...||..|
gi|949700         0 -VLLVEDEEALRAAAGDFLETRGYKIMTARDGTEALSMASKFAERIDVLITDLVMPGISG

ref_rec          58 IEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHF-AKPFDID----EIRDAV
                 60 .........|.....|..|..|..-.......|......|-.|||..|----.||...
gi|949700        59 RVLAQELVKIHPETKVMYMSGYDD-ETVMVNGEIDSSSAFLRKPFRMDALSAKIREVL

ref_rec         111
                118
gi|949700       116
""",
        )
        self.assertTrue(
            np.array_equal(
                np.array(alignment, "U"),
                # fmt: off
np.array([['K', 'I', 'L', 'I', 'V', 'D', 'D', 'Q', 'Y', 'G', 'I', 'R', 'I',
           'L', 'L', 'N', 'E', 'V', 'F', 'N', 'K', 'E', 'G', 'Y', 'Q', 'T',
           'F', 'Q', 'A', 'A', 'N', 'G', 'L', 'Q', 'A', 'L', 'D', 'I', 'V',
           'T', 'K', '-', '-', 'E', 'R', 'P', 'D', 'L', 'V', 'L', 'L', 'D',
           'M', 'K', 'I', 'P', 'G', 'M', 'D', 'G', 'I', 'E', 'I', 'L', 'K',
           'R', 'M', 'K', 'V', 'I', 'D', 'E', 'N', 'I', 'R', 'V', 'I', 'I',
           'M', 'T', 'A', 'Y', 'G', 'E', 'L', 'D', 'M', 'I', 'Q', 'E', 'S',
           'K', 'E', 'L', 'G', 'A', 'L', 'T', 'H', 'F', '-', 'A', 'K', 'P',
           'F', 'D', 'I', 'D', '-', '-', '-', '-', 'E', 'I', 'R', 'D', 'A',
           'V'],
          ['-', 'V', 'L', 'L', 'V', 'E', 'D', 'E', 'E', 'A', 'L', 'R', 'A',
           'A', 'A', 'G', 'D', 'F', 'L', 'E', 'T', 'R', 'G', 'Y', 'K', 'I',
           'M', 'T', 'A', 'R', 'D', 'G', 'T', 'E', 'A', 'L', 'S', 'M', 'A',
           'S', 'K', 'F', 'A', 'E', 'R', 'I', 'D', 'V', 'L', 'I', 'T', 'D',
           'L', 'V', 'M', 'P', 'G', 'I', 'S', 'G', 'R', 'V', 'L', 'A', 'Q',
           'E', 'L', 'V', 'K', 'I', 'H', 'P', 'E', 'T', 'K', 'V', 'M', 'Y',
           'M', 'S', 'G', 'Y', 'D', 'D', '-', 'E', 'T', 'V', 'M', 'V', 'N',
           'G', 'E', 'I', 'D', 'S', 'S', 'S', 'A', 'F', 'L', 'R', 'K', 'P',
           'F', 'R', 'M', 'D', 'A', 'L', 'S', 'A', 'K', 'I', 'R', 'E', 'V',
           'L']], dtype='U')
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.annotations["Matrix"], "EBLOSUM62")
        self.assertAlmostEqual(alignment.annotations["Gap_penalty"], 10.0)
        self.assertAlmostEqual(alignment.annotations["Extend_penalty"], 0.5)
        self.assertEqual(alignment.annotations["Identity"], 35)
        self.assertEqual(alignment.annotations["Similarity"], 70)
        self.assertEqual(alignment.annotations["Gaps"], 18)
        self.assertAlmostEqual(alignment.annotations["Score"], 156.5)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 125))
        self.assertEqual(alignment.sequences[0].id, "ref_rec")
        self.assertEqual(alignment.sequences[1].id, "gi|94970041|receiver")
        self.assertEqual(
            alignment.sequences[0].seq,
            "KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIVTKERPDLVLLDMKIPGMDGIEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHFAKPFDIDEIRDAV",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "TVLLVEDEEGVRKLVRGILSRQGYHVLEATSGEEALEIVRESTQKIDMLLSDVVLVGMSGRELSERLRIQMPSLKVIYMSGYTDDAIVRHGVLTESAEFLQKPFTSDSLLRKVRAVLQKRQ",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                np.array(
                    [
                        [0, 39, 39, 88, 92, 99, 99, 111, 111],
                        [0, 39, 41, 90, 90, 97, 101, 113, 121],
                    ]
                ),
            )
        )
        self.assertEqual(
            alignment[0],
            "KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIV--TKERPDLVLLDMKIPGMDGIEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHFA----KPFDIDEIRDAV--------",
        )
        self.assertEqual(
            alignment[1],
            "TVLLVEDEEGVRKLVRGILSRQGYHVLEATSGEEALEIVRESTQKIDMLLSDVVLVGMSGRELSERLRIQMPSLKVIYMSGYTDDAIVRH----GVLTESAEFLQKPFTSDSLLRKVRAVLQKRQ",
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            ".:|:|:|:.|:|.|:..:.:::||...:|.:|.:||:||  :.::.|::|.|:.:.||.|.|:.:|:::...:::||.|:.|.:..:::.    |.||..|    |||..|.:...|        ",
        )
        self.assertEqual(
            str(alignment),
            """\
ref_rec           0 KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIV--TKERPDLVLLDMKIPGMDG
                  0 ..|.|.|..|.|.|........||....|..|..||.||--.....|..|.|....||.|
gi|949700         0 TVLLVEDEEGVRKLVRGILSRQGYHVLEATSGEEALEIVRESTQKIDMLLSDVVLVGMSG

ref_rec          58 IEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHFA----KPFDIDEIRDAV---
                 60 .|...|.........||.|..|........----|.||..|----|||..|.....|---
gi|949700        60 RELSERLRIQMPSLKVIYMSGYTDDAIVRH----GVLTESAEFLQKPFTSDSLLRKVRAV

ref_rec         111 ----- 111
                120 ----- 125
gi|949700       116 LQKRQ 121
""",
        )
        self.assertTrue(
            np.array_equal(
                np.array(alignment, "U"),
                # fmt: off
np.array([['K', 'I', 'L', 'I', 'V', 'D', 'D', 'Q', 'Y', 'G', 'I', 'R', 'I',
           'L', 'L', 'N', 'E', 'V', 'F', 'N', 'K', 'E', 'G', 'Y', 'Q', 'T',
           'F', 'Q', 'A', 'A', 'N', 'G', 'L', 'Q', 'A', 'L', 'D', 'I', 'V',
           '-', '-', 'T', 'K', 'E', 'R', 'P', 'D', 'L', 'V', 'L', 'L', 'D',
           'M', 'K', 'I', 'P', 'G', 'M', 'D', 'G', 'I', 'E', 'I', 'L', 'K',
           'R', 'M', 'K', 'V', 'I', 'D', 'E', 'N', 'I', 'R', 'V', 'I', 'I',
           'M', 'T', 'A', 'Y', 'G', 'E', 'L', 'D', 'M', 'I', 'Q', 'E', 'S',
           'K', 'E', 'L', 'G', 'A', 'L', 'T', 'H', 'F', 'A', '-', '-', '-',
           '-', 'K', 'P', 'F', 'D', 'I', 'D', 'E', 'I', 'R', 'D', 'A', 'V',
           '-', '-', '-', '-', '-', '-', '-', '-'],
          ['T', 'V', 'L', 'L', 'V', 'E', 'D', 'E', 'E', 'G', 'V', 'R', 'K',
           'L', 'V', 'R', 'G', 'I', 'L', 'S', 'R', 'Q', 'G', 'Y', 'H', 'V',
           'L', 'E', 'A', 'T', 'S', 'G', 'E', 'E', 'A', 'L', 'E', 'I', 'V',
           'R', 'E', 'S', 'T', 'Q', 'K', 'I', 'D', 'M', 'L', 'L', 'S', 'D',
           'V', 'V', 'L', 'V', 'G', 'M', 'S', 'G', 'R', 'E', 'L', 'S', 'E',
           'R', 'L', 'R', 'I', 'Q', 'M', 'P', 'S', 'L', 'K', 'V', 'I', 'Y',
           'M', 'S', 'G', 'Y', 'T', 'D', 'D', 'A', 'I', 'V', 'R', 'H', '-',
           '-', '-', '-', 'G', 'V', 'L', 'T', 'E', 'S', 'A', 'E', 'F', 'L',
           'Q', 'K', 'P', 'F', 'T', 'S', 'D', 'S', 'L', 'L', 'R', 'K', 'V',
           'R', 'A', 'V', 'L', 'Q', 'K', 'R', 'Q']], dtype='U')
                # fmt: on
            )
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_pair_example3(self):
        path = "Emboss/needle_overhang.txt"
        alignments = Align.parse(path, "emboss")
        self.assertEqual(alignments.metadata["Program"], "needle")
        self.assertEqual(alignments.metadata["Rundate"], "Mon 14 Jul 2008 11:45:42")
        self.assertEqual(
            alignments.metadata["Command line"],
            "needle [-asequence] asis:TGTGGTTAGGTTTGGTTTTATTGGGGGCTTGGTTTGGGCCCACCCCAAATAGGGAGTGGGGGTATGACCTCAGATAGACGAGCTTATTTTAGGGCGGCGACTATAATTATTTCGTTTCCTACAAGGATTAAAGTTTTTTCTTTTACTGTGGGAGGGGGTTTGGTATTAAGAAACGCTAGTCCGGATGTGGCTCTCCATGATACTTATTGTGTAGTAGCTCATTTTCATTATGTTCTTCGAATGGGAGCAGTCATTGGTATTTTTTTGGTTTTTTTTTGAAATTTTTAGGTTATTTAGACCATTTTTTTTTGTTTCGCTAATTAGAATTTTATTAGCCTTTGGTTTTTTTTTATTTTTTGGGGTTAAGACAAGGTGTCGTTGAATTAGTTTAGCAAAATACTGCTTAAGGTAGGCTATAGGATCTACCTTTTATCTTTCTAATCTTTTGTTTTAGTATAATTGGTCTTCGATTCAACAATTTTTAGTCTTCAGTCTTTTTTTTTATTTTGAAAAGGTTTTAACACTCTTGGTTTTGGAGGCTTTGGCTTTCTTCTTACTCTTAGGAGGATGGGCGCTAGAAAGAGTTTTAAGAGGGTGTGAAAGGGGGTTAATAGC [-bsequence] asis:TTATTAATCTTATGGTTTTGCCGTAAAATTTCTTTCTTTATTTTTTATTGTTAGGATTTTGTTGATTTTATTTTTCTCAAGAATTTTTAGGTCAATTAGACCGGCTTATTTTTTTGTCAGTGTTTAAAGTTTTATTAATTTTTGGGGGGGGGGGGAGACGGGGTGTTATCTGAATTAGTTTTTGGGAGTCTCTAGACATCTCATGGGTTGGCCGGGGGCCTGCCGTCTATAGTTCTTATTCCTTTTAAGGGAGTAAGAATTTCGATTCAGCAACTTTAGTTCACAGTCTTTTTTTTTATTAAGAAAGGTTT -filter",
        )
        self.assertEqual(alignments.metadata["Align_format"], "srspair")
        self.assertEqual(alignments.metadata["Report_file"], "stdout")
        alignment = next(alignments)
        self.assertEqual(alignment.annotations["Matrix"], "EDNAFULL")
        self.assertAlmostEqual(alignment.annotations["Gap_penalty"], 10.0)
        self.assertAlmostEqual(alignment.annotations["Extend_penalty"], 0.5)
        self.assertEqual(alignment.annotations["Identity"], 210)
        self.assertEqual(alignment.annotations["Similarity"], 210)
        self.assertEqual(alignment.annotations["Gaps"], 408)
        self.assertAlmostEqual(alignment.annotations["Score"], 561.0)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 667))
        self.assertEqual(alignment.sequences[0].id, "asis")
        self.assertEqual(alignment.sequences[1].id, "asis")
        self.assertEqual(
            alignment.sequences[0].seq,
            "TGTGGTTAGGTTTGGTTTTATTGGGGGCTTGGTTTGGGCCCACCCCAAATAGGGAGTGGGGGTATGACCTCAGATAGACGAGCTTATTTTAGGGCGGCGACTATAATTATTTCGTTTCCTACAAGGATTAAAGTTTTTTCTTTTACTGTGGGAGGGGGTTTGGTATTAAGAAACGCTAGTCCGGATGTGGCTCTCCATGATACTTATTGTGTAGTAGCTCATTTTCATTATGTTCTTCGAATGGGAGCAGTCATTGGTATTTTTTTGGTTTTTTTTTGAAATTTTTAGGTTATTTAGACCATTTTTTTTTGTTTCGCTAATTAGAATTTTATTAGCCTTTGGTTTTTTTTTATTTTTTGGGGTTAAGACAAGGTGTCGTTGAATTAGTTTAGCAAAATACTGCTTAAGGTAGGCTATAGGATCTACCTTTTATCTTTCTAATCTTTTGTTTTAGTATAATTGGTCTTCGATTCAACAATTTTTAGTCTTCAGTCTTTTTTTTTATTTTGAAAAGGTTTTAACACTCTTGGTTTTGGAGGCTTTGGCTTTCTTCTTACTCTTAGGAGGATGGGCGCTAGAAAGAGTTTTAAGAGGGTGTGAAAGGGGGTTAATAGC",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "TTATTAATCTTATGGTTTTGCCGTAAAATTTCTTTCTTTATTTTTTATTGTTAGGATTTTGTTGATTTTATTTTTCTCAAGAATTTTTAGGTCAATTAGACCGGCTTATTTTTTTGTCAGTGTTTAAAGTTTTATTAATTTTTGGGGGGGGGGGGAGACGGGGTGTTATCTGAATTAGTTTTTGGGAGTCTCTAGACATCTCATGGGTTGGCCGGGGGCCTGCCGTCTATAGTTCTTATTCCTTTTAAGGGAGTAAGAATTTCGATTCAGCAACTTTAGTTCACAGTCTTTTTTTTTATTAAGAAAGGTTT",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[  0, 162, 169, 201, 210, 210, 220, 222, 236, 240,
                           244, 253, 277, 277, 278, 279, 300, 300, 310, 310,
                           314, 320, 334, 351, 357, 357, 379, 379, 390, 403,
                           405, 407, 418, 418, 442, 442, 447, 447, 448, 452,
                           455, 455, 460, 465, 478, 479, 509, 510, 518, 615],
                          [  0,   0,   7,   7,  16,  22,  32,  32,  46,  46,
                            50,  50,  74,  80,  81,  81, 102, 107, 117, 119,
                           123, 123, 137, 137, 143, 147, 169, 170, 181, 181,
                           183, 183, 194, 215, 239, 241, 246, 250, 251, 251,
                           254, 255, 260, 260, 273, 273, 303, 303, 311, 311]])
                # fmt: on
            )
        )
        self.assertEqual(
            alignment[0],
            "TGTGGTTAGGTTTGGTTTTATTGGGGGCTTGGTTTGGGCCCACCCCAAATAGGGAGTGGGGGTATGACCTCAGATAGACGAGCTTATTTTAGGGCGGCGACTATAATTATTTCGTTTCCTACAAGGATTAAAGTTTTTTCTTTTACTGTGGGAGGGGGTTTGGTATTAAGAAACGCTAGTCCGGATGTGGCTCTCCATGATACTTATTGT------GTAGTAGCTCATTTTCATTATGTTCTTCGAATGGGAGCAGTCATTGGTATTTTTTTGGTTTTTTTTT------GAAATTTTTAGGTTATTTAGACC-----ATTTTTTTTT--GTTTCGCTAATTAGAATTTTATTAGCCTTTGGTTTTTTTTTATTTTT----TGGGGTTAAGACAAGGTGTCGT-TGAATTAGTTTAGCAAAATACTGCTTAAGGTAGGCTATA---------------------GGATCTACCTTTTATCTTTCTAAT--CTTTT----GTTTTAGT-ATAATTGGTCTTCGATTCAACAATTTTTAGTCTTCAGTCTTTTTTTTTATTTTGAAAAGGTTTTAACACTCTTGGTTTTGGAGGCTTTGGCTTTCTTCTTACTCTTAGGAGGATGGGCGCTAGAAAGAGTTTTAAGAGGGTGTGAAAGGGGGTTAATAGC",
        )
        self.assertEqual(
            alignment[1],
            "------------------------------------------------------------------------------------------------------------------------------------------------------------------TTATTAA--------------------------------TCTTATGGTTTTGCCGTAAAATTTC--TTTCTTTATTTTTT----ATTG---------TTAGGATTTTGTTGATTTTATTTTTCTCAAG-AATTTTTAGGTCAATTAGACCGGCTTATTTTTTTGTCAGTGT------TTAAAGTTTTATTA-----------------ATTTTTGGGGGGGGGGGGAGACGGGGTGTTATCTGAATTAGTTT-------------TT--GGGAGTCTCTAGACATCTCATGGGTTGGCCGGGGGCCTGCCGTCTATAGTTCTTATTCCTTTTAAGGG----AGTAAGAAT-----TTCGATTCAGCAA-CTTTAGTTCACAGTCTTTTTTTTTATTAAG-AAAGGTTT-------------------------------------------------------------------------------------------------",
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            "                                                                                                                                                                  .||||||                                .|||||.||      |||..|..||  ||||.||||.||.|    ||.|         ||.|.|||||.|||.||||.||||      | |||||||||||.|.|||||||     ||||||||.|  ||.|      |||.|.||||||||                 ||||||    .||||...||||..|||||..| |||||||||||             ||  ||.||.||.||                     ||..||.||.|.|||..||||.||  |||||    |    ||| |.|||     |||||||||.||| .||||||...|||||||||||||||||..| ||||||||                                                                                                 ",
        )
        self.assertEqual(
            str(alignment),
            """\
asis              0 TGTGGTTAGGTTTGGTTTTATTGGGGGCTTGGTTTGGGCCCACCCCAAATAGGGAGTGGG
                  0 ------------------------------------------------------------
asis              0 ------------------------------------------------------------

asis             60 GGTATGACCTCAGATAGACGAGCTTATTTTAGGGCGGCGACTATAATTATTTCGTTTCCT
                 60 ------------------------------------------------------------
asis              0 ------------------------------------------------------------

asis            120 ACAAGGATTAAAGTTTTTTCTTTTACTGTGGGAGGGGGTTTGGTATTAAGAAACGCTAGT
                120 ------------------------------------------.||||||-----------
asis              0 ------------------------------------------TTATTAA-----------

asis            180 CCGGATGTGGCTCTCCATGATACTTATTGT------GTAGTAGCTCATTTTCATTATGTT
                180 ---------------------.|||||.||------|||..|..||--||||.||||.||
asis              7 ---------------------TCTTATGGTTTTGCCGTAAAATTTC--TTTCTTTATTTT

asis            234 CTTCGAATGGGAGCAGTCATTGGTATTTTTTTGGTTTTTTTTT------GAAATTTTTAG
                240 .|----||.|---------||.|.|||||.|||.||||.||||------|-|||||||||
asis             44 TT----ATTG---------TTAGGATTTTGTTGATTTTATTTTTCTCAAG-AATTTTTAG

asis            288 GTTATTTAGACC-----ATTTTTTTTT--GTTTCGCTAATTAGAATTTTATTAGCCTTTG
                300 ||.|.|||||||-----||||||||.|--||.|------|||.|.||||||||-------
asis             90 GTCAATTAGACCGGCTTATTTTTTTGTCAGTGT------TTAAAGTTTTATTA-------

asis            341 GTTTTTTTTTATTTTT----TGGGGTTAAGACAAGGTGTCGT-TGAATTAGTTTAGCAAA
                360 ----------||||||----.||||...||||..|||||..|-|||||||||||------
asis            137 ----------ATTTTTGGGGGGGGGGGGAGACGGGGTGTTATCTGAATTAGTTT------

asis            396 ATACTGCTTAAGGTAGGCTATA---------------------GGATCTACCTTTTATCT
                420 -------||--||.||.||.||---------------------||..||.||.|.|||..
asis            181 -------TT--GGGAGTCTCTAGACATCTCATGGGTTGGCCGGGGGCCTGCCGTCTATAG

asis            435 TTCTAAT--CTTTT----GTTTTAGT-ATAATTGGTCTTCGATTCAACAATTTTTAGTCT
                480 ||||.||--|||||----|----|||-|.|||-----|||||||||.|||-.||||||..
asis            232 TTCTTATTCCTTTTAAGGG----AGTAAGAAT-----TTCGATTCAGCAA-CTTTAGTTC

asis            488 TCAGTCTTTTTTTTTATTTTGAAAAGGTTTTAACACTCTTGGTTTTGGAGGCTTTGGCTT
                540 .|||||||||||||||||..|-||||||||------------------------------
asis            282 ACAGTCTTTTTTTTTATTAAG-AAAGGTTT------------------------------

asis            548 TCTTCTTACTCTTAGGAGGATGGGCGCTAGAAAGAGTTTTAAGAGGGTGTGAAAGGGGGT
                600 ------------------------------------------------------------
asis            311 ------------------------------------------------------------

asis            608 TAATAGC 615
                660 ------- 667
asis            311 ------- 311
""",
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_needle_asis(self):
        path = "Emboss/needle_asis.txt"
        alignments = Align.parse(path, "emboss")
        self.assertEqual(alignments.metadata["Program"], "needle")
        self.assertEqual(alignments.metadata["Rundate"], "Mon 14 Jul 2008 11:37:15")
        self.assertEqual(
            alignments.metadata["Command line"],
            "needle [-asequence] asis:TATTTTTTGGATTTTTTTCTAGATTTTCTAGGTTATTTAAACCGTTTTTTTTTAATTTAGTGTTTGAGTTTTGACAGGTCTCCACTTTGGGGGCTCCATCGCAAGGAAATTAGAATTCTTATACTTGGTTCTCTTTCCCAGGGACTCCAAGGATCTTTTCATTAGTTTGGATTTTGGTGTTTTCTTTAATTTTGTTAAGAAACAAATCCTTTCTAGAGTTTTTTCTAGCATTATGTTTTTTTTTCTCCTTATCTAAGGGGGTTTGTCGAGGTTTCTTAAATCTTTTTTTCTCTGGGTTTTAAAATTGTTTAAATTTTTTTGACCGAGGGGTTGGGGTGGTTTTCTCATGATAACAGGGGCTGGTGCTTTAGATCCTACCTCTACTGACCCGGGGTCTGCTACTGTGGCTTCTGATGAAGATCCACAGTATGCGCCTACGGAARCTCGGCAGTTTGGTGTTCGAAATCCAGCCCCTCGAATTAATACTCTTGTGCAGGTGGTTGACGAGCGCGGTATCGAATTGCAAAATTTGGGGCGGGACCCCGCTGTTCCGCCTGTTGCTCCGGGGGGGGCAGGTTAATCCTCCAGTCGTCTCCTTTTGGGGGCGTCTTTGACGGGGGTTTAAATCTTTCTTTGGTTGTGGATAGGATTTTTTTTCTAATATCGATCCTACCTGTTTTGGCGGGGCTATTACTTTGTTACTTTTGACCGAAATTTTAATGGAAATTTCTTTGATTCAAATGAATCCCTTAGTTTTCCAACACTTTTTTTTGGTTTTTTTAGGGATAGTCTACGCTGTGGTTAGGTTTGGTTTTATTGGGGGCTTGGTTTGGGCCCACCCCAAATAGGGAGTGGGGGTATGACCTCAGATAGACGAGCTTATTTTAGGGCGGCGACTATAATTATTTCGTTTCCTACAAGGATTAAAGTTTTTTCTTTTACTGTGGGAGGGGGTTTGGTATTAAGAAACGCTAGTCCGGATGTGGCTCTCCATGATACTTATTGTGTAGTAGCTCATTTTCATTATGTTCTTCGAATGGGAGCAGTCATTGGTATTTTTTTGGTTTTTTTTTGAAATTTTTAGGTTATTTAGACCATTTTTTTTTGTTTCGCTAATTAGAATTTTATTAGCCTTTGGTTTTTTTTTATTTTTTGGGGTTAAGACAAGGTGTCGTTGAATTAGTTTAGCAAAATACTGCTTAAGGTAGGCTATAGGATCTACCTTTTATCTTTCTAATCTTTTGTTTTAGTATAATTGGTCTTCGATTCAACAATTTTTAGTCTTCAGTCTTTTTTTTTATTTTGAAAAGGTTTTAACACTCTTGGTTTTGGAGGCTTTGGCTTTCTTCTTACTCTTAGGAGGATGGGCGCTAGAAAGAGTTTTAAGAGGGTGTGAAAGGGGGTTAATAGCAGGATTTGCTTTTTTAACTTATACTGGTTCGTAACGCATTAGCTCAACTCTCTCTTGTAGTTCTAGCAGCCGCCTTTTCTTTGTTGGGGGAGGGTTTAGGAGGAGTCTTTTTTTTCCTAACCCAAGGTGTTTCTTTCTTTTTTTCTTTAAAGTTCTTGACTGTTGGCACTTGTCTCCATAAATTTTCTTTCTTGTAAAGGGCTCCTAAGGCTTCTTGTTTCTGAATTCCTCTTTTCTTTTATTCTGTTTTGAGCTTATTTTTCTTGTTAGCTATTACGTAGGCATAGGGCAAATAATTTTTTTTTCTGCTCTCATTATTCCTTCTCCCTGCTTGTTTCACCCTGTGGGCTCTTTGAGCCCCACTAAGTGAGCGGGGCTCCTGCTTCCGCTCAATTAAATTTTGGTGGGTATTGAGTCTCAGAGGGACTATGATATAGGTTCAGATTGATGGACCTAATCAATCAATTGTATCGCTATACAATCTAGTACCCCTACCAGGGTACCAAGAGAGAGATAACTAGGGTGAATACTACGACTTAGATGTAGTGTTTAAGTTTCTACGGGCTACAGAGAAGCTACCCGCAGGGTAYATATTTGTTCATTACATATTTGTTGACTTTTCTATCTCTGCTTTTACTTTTTTATTTATTTTTAAATCTTTTTAACTTCAGCTGTTTTTCCTTATCTATTTGACGTAGGCATAGGAAAGTTAACGAATTTTGTAATATTTTTAATTATTTTGTATAGTATACAGGGTAGTGGTATGTAATAGGTAAATTCCATAAGTTCATTATAGTTTATCAGTTGAGAGGAATTTAGTATAAGAAGGCCCATTGGGGCTCTTGTCTTATCCAAGAACTGGTAAGATTTAATTCTACCGGGACGGTAGAAATCGGGCAGAGCATGATCTATTTCTTCGGGTATGGCTATAGGGACTAGGTGCTAGGGAGGTATTAGGGCACCGCTCTTTATACAATCTCCATAGATACAACCAGGTCAACTAGGACAACGGAGGACGTTGACAGAGCATAAATAGCGATAGCGTACAAGATAWAATAGGGGCAGTGGTAGCGAAGCGTAGAAGAAAAAATAAGAGTATTGTTTGTAAATAATTCTTTTTTTAGTTTTTAAATATTCTTTTTTTAGGTGGTGTGTGGTTAGGTATGGGGTTAAGGGTGTGGCAAAGAGAAATGTTTATTAAACATTCTTATGGCCGTAGATAGCATATCGATTATACGAGACCTTCGTAAGATCAATCCCCACTAGCATTGCTCATACAGGTTAACTCAATAGGAGGAGCTGGGGTAGAACGTTTCTAGTTCGGGGGTAACCGCAGTTCAATGAAAGTGACGACGTCGGATGGAACAAACTTAATACCACCAGTTGTGCTAACGATTGTTATCTCAATCTATCCCAACAGGCCCCCAGGTAGTGATGAGTGGTGGAATGGTACAGGGTACCAGTGGGTGAAGAGCGTCACGAACCAGGGAATACGGAGTACAGAGTTGAGCGCCCGGGGCTCCGCCCCCGGCTTTTATAGCGCGAGACGTGGTCAGTCGATTCAGCGTTAGGTTTAAACTCCTTTGGCAAAGATTGACTCTAGCGATCCAGAGACCCTGCCTGGCATAAAAGTCTTTATWAACACCAGTAGGTTCAATAAGGTAGTAATCCAATAGAATGGAAAACTCAAGATCTAATCTCTCGAYTTCCTAGTGTCATGGAAATCAGCCAGGTTCTCTTCATCTGCAACAGTAGAAGAAGAAGAGAGACTAGCGAGAGAGTCTTATGGCGGAGACGCTAAGGCTTAAATGTAATGTAGATAACCCCTTACGGAACACTAGAGTGCGACGTAGACTACATAATCCCTCAGGGATATTAGCTCTGCTCGATTAACAATAGCATACTTTGTTACACGGAGTGTATCTAGGGGGAATAATACTAACTTACTTAGCACTATCGCGATGCTACGCATTCGCTCTTTCGCTAAATAAGATACGACGATGAGTGGTTGGTGGAGAGAATAACCGATTCTAACTTGATAATTCGCATGAAATAATTTTTTTATTTTGTTTTTTTTTTGCTCTTAATTTTAGWGGGRGTGTTTATTTTTATTCTAATAAAAAGGATCCGTTGAA [-bsequence] asis:TTATTAATCTTATGGTTTTGCCGTAAAATTTCTTTCTTTATTTTTTATTGTTAGGATTTTGTTGATTTTATTTTTCTCAAGAATTTTTAGGTCAATTAGACCGGCTTATTTTTTTGTCAGTGTTTAAAGTTTTATTAATTTTTGGGGGGGGGGGGAGACGGGGTGTTATCTGAATTAGTTTTTGGGAGTCTCTAGACATCTCATGGGTTGGCCGGGGGCCTGCCGTCTATAGTTCTTATTCCTTTTAAGGGAGTAAGAATTTCGATTCAGCAACTTTAGTTCACAGTCTTTTTTTTTATTAAGAAAGGTTTTAATATTCTTGTGGTTTTGAACCTTTAGGTTTCTTTCTTTACCTTCGAGGGATTGGGCACTAGAATGAGTTTTAAGAGTGTGTGAAAGGGGGCTTGATAGCAGGGGAATGCTTTTTTAACTTATACTGGCTCGTAACGCATCAGTTCAACTCTCTCTTGCAGTTCTAGCAGCCGCCTTTTTTTTGTTGGGGGGGGGTTAAGAGAGTGTTTTTTTTCTAATCCAAGGGTCTTACTTTCTTTCTTTCTTTAAAAATTCTTTGGCTGTCGACACCTTTCTCTCCCGTCAGTCTCATGGTTTCTGGCTCTCTTGGGCTTTTTTTGTTTGTGAATGCCTCTTTTTTTTATTCTGTTTTGAGCTTATTTTTCTTGTTTACTATTACGTAGGTATAGGGCAAATAATTTTTTTTTCGCGTCTCTTGGCATGCCCATTACTCTAGTTTTATTCCCGGGCTTCTTCTCTCACCCTAGAGGGCTCTTTGAGCCCACACTCAAGTGAGCGGGGCTCCCGCTTCCGCTCAATTAAATTTGGTGGGTATTGAGTCTCAGAGGGACTATGATATAGGTTCAGATTGATGGACCTAGTCAATCAATTGTATCGCTATACAATCTAGTACCCCTACCAGGGTACCAGGAGAGAGATAACTAGGGTGAATACTACGACTTAGATGTACTGTTTAAGTTTCTACGGGCTACAGAGAAGCTACCCGCAGGGTATATATTTGCTCATTACATATTTGTTGATTTTTCTATGTCCGCTTTACTTTTTATATTTTTTTAACTTCAGCTGTTTTTCCTTATCTATTTGACGTAGGCATAGGAAAGTTAACGAATTTTGTAATATTTTTAATTATTTTGTATAGTATACAGGGTAGTGGTATGTAATAGGTAAATTCCATAAGTTCATTATAGTCTATCAGTTGAGAGGAATTTAGTATAAGAAAGCCTGTCAGGGCTCTTGCCTTATCCAAGAACTGGTAAGGATTTCTTGACAGAGGGACTCTGTCAAATCGGGCAGAGCATGATCTATTTCTTCGGGTATGGTTATAAGGCTTAGGTGCTTGGAGGGTATTAGGGCACCGCTCTTAATACAGTCTCCATAGGTGTAACCAGGTCAACTAGGACAACGGAGGACGTTGACAAAGCATGGATAGCGATAGCGTAGAAGATAAAATGGGGCAGTGGTAGCGAAGCGTAGAAGAAAAAATAAGAGTATTGTTTGTAAATAATTCTTTTTTTAGTTTTTAAATATTCTTTTTTTAGGTGGTGTGTGGTTAGGTATGGGGTTAGGGGAGTGGCAAAGAGAAGTGTTTATTAAACATTCTTATGGCCGTAGATAGCATATCGATTATACGAGACCTTCGTAAGATCAATCCCCACTAGCATTGCTCATACAGGTTAACTCAATAGGAGGAGCTGGGGTAGAACGTATCTAGTTCGGGGGTAACCGCAGTTCAATGAAAGTGACGACGTCGGATGGAACAAACTTAATACCACCAGTTGTGCTAACGATTGTTATCTCAATCTATCCCAACAGGCCCCCAGGTAGTGATGAGTGGTGGAATGGTACAGGGTACCAGTGGGTGAAGAGCGTCACGAACCAGGGAATACGGAGTACAGAGTTGAGCGCCCGGGGCTCCGCCCCCGGCTTTTATAGCGCGAGACGTGGTCAGTCGATTCAGCGTTAGGTTTTAAACTCCTTTGGCAAAGATTGATTCTAGCGATCCAGAGACCCTGCCTGGCATAAAAGTCTTTATTAGCACCAGTAGGTTCAATAAGGTAGTAGTCCAATAGAATGGAAAACTCGAGATCTAATCTCTCGATTTCCTAGTGTCATGGAAATCAGCCAGGTTCTCTTCATCTGCAACAGTAGAAGAAGAAGAGAGGCTAGCGAGAGAGTCTTATGGCGGAGACGCTAAGGCTTAAATGTAATGTAGATAACCCCTTACGGAACACTTGAGTGCGACGTAGACTACATAATCCCTCAGGGATATTAGCTCTGCTCGATTAACAATAGCATACTTTGTTACACGGAGTGTATCTGGGGGGAATAATACTAACTTACTTAGCACTATCGCGATGCTACGCATTCGCTCTTTCGCTAAATAAGATACGACGATGAGTGGTTGGTGGAGAGAATAACCGATTCTAACTTGATAATTCGCATGAAATAATTTTTTATTTGTTTTTTTTTTTGCTCTTAATTTTAGAGGATGTTTATTTTTATTCTAATAAAAAGGATCCGTTGAA -filter",
        )
        self.assertEqual(alignments.metadata["Align_format"], "srspair")
        self.assertEqual(alignments.metadata["Report_file"], "stdout")
        alignment = next(alignments)
        self.assertEqual(alignment.annotations["Matrix"], "EDNAFULL")
        self.assertAlmostEqual(alignment.annotations["Gap_penalty"], 10.0)
        self.assertAlmostEqual(alignment.annotations["Extend_penalty"], 0.5)
        self.assertEqual(alignment.annotations["Identity"], 2296)
        self.assertEqual(alignment.annotations["Similarity"], 2301)
        self.assertEqual(alignment.annotations["Gaps"], 1202)
        self.assertAlmostEqual(alignment.annotations["Score"], 10155.0)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 3653))
        self.assertEqual(alignment.sequences[0].id, "asis")
        self.assertEqual(alignment.sequences[1].id, "asis")
        self.assertEqual(
            alignment.sequences[0].seq,
            "TATTTTTTGGATTTTTTTCTAGATTTTCTAGGTTATTTAAACCGTTTTTTTTTAATTTAGTGTTTGAGTTTTGACAGGTCTCCACTTTGGGGGCTCCATCGCAAGGAAATTAGAATTCTTATACTTGGTTCTCTTTCCCAGGGACTCCAAGGATCTTTTCATTAGTTTGGATTTTGGTGTTTTCTTTAATTTTGTTAAGAAACAAATCCTTTCTAGAGTTTTTTCTAGCATTATGTTTTTTTTTCTCCTTATCTAAGGGGGTTTGTCGAGGTTTCTTAAATCTTTTTTTCTCTGGGTTTTAAAATTGTTTAAATTTTTTTGACCGAGGGGTTGGGGTGGTTTTCTCATGATAACAGGGGCTGGTGCTTTAGATCCTACCTCTACTGACCCGGGGTCTGCTACTGTGGCTTCTGATGAAGATCCACAGTATGCGCCTACGGAARCTCGGCAGTTTGGTGTTCGAAATCCAGCCCCTCGAATTAATACTCTTGTGCAGGTGGTTGACGAGCGCGGTATCGAATTGCAAAATTTGGGGCGGGACCCCGCTGTTCCGCCTGTTGCTCCGGGGGGGGCAGGTTAATCCTCCAGTCGTCTCCTTTTGGGGGCGTCTTTGACGGGGGTTTAAATCTTTCTTTGGTTGTGGATAGGATTTTTTTTCTAATATCGATCCTACCTGTTTTGGCGGGGCTATTACTTTGTTACTTTTGACCGAAATTTTAATGGAAATTTCTTTGATTCAAATGAATCCCTTAGTTTTCCAACACTTTTTTTTGGTTTTTTTAGGGATAGTCTACGCTGTGGTTAGGTTTGGTTTTATTGGGGGCTTGGTTTGGGCCCACCCCAAATAGGGAGTGGGGGTATGACCTCAGATAGACGAGCTTATTTTAGGGCGGCGACTATAATTATTTCGTTTCCTACAAGGATTAAAGTTTTTTCTTTTACTGTGGGAGGGGGTTTGGTATTAAGAAACGCTAGTCCGGATGTGGCTCTCCATGATACTTATTGTGTAGTAGCTCATTTTCATTATGTTCTTCGAATGGGAGCAGTCATTGGTATTTTTTTGGTTTTTTTTTGAAATTTTTAGGTTATTTAGACCATTTTTTTTTGTTTCGCTAATTAGAATTTTATTAGCCTTTGGTTTTTTTTTATTTTTTGGGGTTAAGACAAGGTGTCGTTGAATTAGTTTAGCAAAATACTGCTTAAGGTAGGCTATAGGATCTACCTTTTATCTTTCTAATCTTTTGTTTTAGTATAATTGGTCTTCGATTCAACAATTTTTAGTCTTCAGTCTTTTTTTTTATTTTGAAAAGGTTTTAACACTCTTGGTTTTGGAGGCTTTGGCTTTCTTCTTACTCTTAGGAGGATGGGCGCTAGAAAGAGTTTTAAGAGGGTGTGAAAGGGGGTTAATAGCAGGATTTGCTTTTTTAACTTATACTGGTTCGTAACGCATTAGCTCAACTCTCTCTTGTAGTTCTAGCAGCCGCCTTTTCTTTGTTGGGGGAGGGTTTAGGAGGAGTCTTTTTTTTCCTAACCCAAGGTGTTTCTTTCTTTTTTTCTTTAAAGTTCTTGACTGTTGGCACTTGTCTCCATAAATTTTCTTTCTTGTAAAGGGCTCCTAAGGCTTCTTGTTTCTGAATTCCTCTTTTCTTTTATTCTGTTTTGAGCTTATTTTTCTTGTTAGCTATTACGTAGGCATAGGGCAAATAATTTTTTTTTCTGCTCTCATTATTCCTTCTCCCTGCTTGTTTCACCCTGTGGGCTCTTTGAGCCCCACTAAGTGAGCGGGGCTCCTGCTTCCGCTCAATTAAATTTTGGTGGGTATTGAGTCTCAGAGGGACTATGATATAGGTTCAGATTGATGGACCTAATCAATCAATTGTATCGCTATACAATCTAGTACCCCTACCAGGGTACCAAGAGAGAGATAACTAGGGTGAATACTACGACTTAGATGTAGTGTTTAAGTTTCTACGGGCTACAGAGAAGCTACCCGCAGGGTAYATATTTGTTCATTACATATTTGTTGACTTTTCTATCTCTGCTTTTACTTTTTTATTTATTTTTAAATCTTTTTAACTTCAGCTGTTTTTCCTTATCTATTTGACGTAGGCATAGGAAAGTTAACGAATTTTGTAATATTTTTAATTATTTTGTATAGTATACAGGGTAGTGGTATGTAATAGGTAAATTCCATAAGTTCATTATAGTTTATCAGTTGAGAGGAATTTAGTATAAGAAGGCCCATTGGGGCTCTTGTCTTATCCAAGAACTGGTAAGATTTAATTCTACCGGGACGGTAGAAATCGGGCAGAGCATGATCTATTTCTTCGGGTATGGCTATAGGGACTAGGTGCTAGGGAGGTATTAGGGCACCGCTCTTTATACAATCTCCATAGATACAACCAGGTCAACTAGGACAACGGAGGACGTTGACAGAGCATAAATAGCGATAGCGTACAAGATAWAATAGGGGCAGTGGTAGCGAAGCGTAGAAGAAAAAATAAGAGTATTGTTTGTAAATAATTCTTTTTTTAGTTTTTAAATATTCTTTTTTTAGGTGGTGTGTGGTTAGGTATGGGGTTAAGGGTGTGGCAAAGAGAAATGTTTATTAAACATTCTTATGGCCGTAGATAGCATATCGATTATACGAGACCTTCGTAAGATCAATCCCCACTAGCATTGCTCATACAGGTTAACTCAATAGGAGGAGCTGGGGTAGAACGTTTCTAGTTCGGGGGTAACCGCAGTTCAATGAAAGTGACGACGTCGGATGGAACAAACTTAATACCACCAGTTGTGCTAACGATTGTTATCTCAATCTATCCCAACAGGCCCCCAGGTAGTGATGAGTGGTGGAATGGTACAGGGTACCAGTGGGTGAAGAGCGTCACGAACCAGGGAATACGGAGTACAGAGTTGAGCGCCCGGGGCTCCGCCCCCGGCTTTTATAGCGCGAGACGTGGTCAGTCGATTCAGCGTTAGGTTTAAACTCCTTTGGCAAAGATTGACTCTAGCGATCCAGAGACCCTGCCTGGCATAAAAGTCTTTATWAACACCAGTAGGTTCAATAAGGTAGTAATCCAATAGAATGGAAAACTCAAGATCTAATCTCTCGAYTTCCTAGTGTCATGGAAATCAGCCAGGTTCTCTTCATCTGCAACAGTAGAAGAAGAAGAGAGACTAGCGAGAGAGTCTTATGGCGGAGACGCTAAGGCTTAAATGTAATGTAGATAACCCCTTACGGAACACTAGAGTGCGACGTAGACTACATAATCCCTCAGGGATATTAGCTCTGCTCGATTAACAATAGCATACTTTGTTACACGGAGTGTATCTAGGGGGAATAATACTAACTTACTTAGCACTATCGCGATGCTACGCATTCGCTCTTTCGCTAAATAAGATACGACGATGAGTGGTTGGTGGAGAGAATAACCGATTCTAACTTGATAATTCGCATGAAATAATTTTTTTATTTTGTTTTTTTTTTGCTCTTAATTTTAGWGGGRGTGTTTATTTTTATTCTAATAAAAAGGATCCGTTGAA",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "TTATTAATCTTATGGTTTTGCCGTAAAATTTCTTTCTTTATTTTTTATTGTTAGGATTTTGTTGATTTTATTTTTCTCAAGAATTTTTAGGTCAATTAGACCGGCTTATTTTTTTGTCAGTGTTTAAAGTTTTATTAATTTTTGGGGGGGGGGGGAGACGGGGTGTTATCTGAATTAGTTTTTGGGAGTCTCTAGACATCTCATGGGTTGGCCGGGGGCCTGCCGTCTATAGTTCTTATTCCTTTTAAGGGAGTAAGAATTTCGATTCAGCAACTTTAGTTCACAGTCTTTTTTTTTATTAAGAAAGGTTTTAATATTCTTGTGGTTTTGAACCTTTAGGTTTCTTTCTTTACCTTCGAGGGATTGGGCACTAGAATGAGTTTTAAGAGTGTGTGAAAGGGGGCTTGATAGCAGGGGAATGCTTTTTTAACTTATACTGGCTCGTAACGCATCAGTTCAACTCTCTCTTGCAGTTCTAGCAGCCGCCTTTTTTTTGTTGGGGGGGGGTTAAGAGAGTGTTTTTTTTCTAATCCAAGGGTCTTACTTTCTTTCTTTCTTTAAAAATTCTTTGGCTGTCGACACCTTTCTCTCCCGTCAGTCTCATGGTTTCTGGCTCTCTTGGGCTTTTTTTGTTTGTGAATGCCTCTTTTTTTTATTCTGTTTTGAGCTTATTTTTCTTGTTTACTATTACGTAGGTATAGGGCAAATAATTTTTTTTTCGCGTCTCTTGGCATGCCCATTACTCTAGTTTTATTCCCGGGCTTCTTCTCTCACCCTAGAGGGCTCTTTGAGCCCACACTCAAGTGAGCGGGGCTCCCGCTTCCGCTCAATTAAATTTGGTGGGTATTGAGTCTCAGAGGGACTATGATATAGGTTCAGATTGATGGACCTAGTCAATCAATTGTATCGCTATACAATCTAGTACCCCTACCAGGGTACCAGGAGAGAGATAACTAGGGTGAATACTACGACTTAGATGTACTGTTTAAGTTTCTACGGGCTACAGAGAAGCTACCCGCAGGGTATATATTTGCTCATTACATATTTGTTGATTTTTCTATGTCCGCTTTACTTTTTATATTTTTTTAACTTCAGCTGTTTTTCCTTATCTATTTGACGTAGGCATAGGAAAGTTAACGAATTTTGTAATATTTTTAATTATTTTGTATAGTATACAGGGTAGTGGTATGTAATAGGTAAATTCCATAAGTTCATTATAGTCTATCAGTTGAGAGGAATTTAGTATAAGAAAGCCTGTCAGGGCTCTTGCCTTATCCAAGAACTGGTAAGGATTTCTTGACAGAGGGACTCTGTCAAATCGGGCAGAGCATGATCTATTTCTTCGGGTATGGTTATAAGGCTTAGGTGCTTGGAGGGTATTAGGGCACCGCTCTTAATACAGTCTCCATAGGTGTAACCAGGTCAACTAGGACAACGGAGGACGTTGACAAAGCATGGATAGCGATAGCGTAGAAGATAAAATGGGGCAGTGGTAGCGAAGCGTAGAAGAAAAAATAAGAGTATTGTTTGTAAATAATTCTTTTTTTAGTTTTTAAATATTCTTTTTTTAGGTGGTGTGTGGTTAGGTATGGGGTTAGGGGAGTGGCAAAGAGAAGTGTTTATTAAACATTCTTATGGCCGTAGATAGCATATCGATTATACGAGACCTTCGTAAGATCAATCCCCACTAGCATTGCTCATACAGGTTAACTCAATAGGAGGAGCTGGGGTAGAACGTATCTAGTTCGGGGGTAACCGCAGTTCAATGAAAGTGACGACGTCGGATGGAACAAACTTAATACCACCAGTTGTGCTAACGATTGTTATCTCAATCTATCCCAACAGGCCCCCAGGTAGTGATGAGTGGTGGAATGGTACAGGGTACCAGTGGGTGAAGAGCGTCACGAACCAGGGAATACGGAGTACAGAGTTGAGCGCCCGGGGCTCCGCCCCCGGCTTTTATAGCGCGAGACGTGGTCAGTCGATTCAGCGTTAGGTTTTAAACTCCTTTGGCAAAGATTGATTCTAGCGATCCAGAGACCCTGCCTGGCATAAAAGTCTTTATTAGCACCAGTAGGTTCAATAAGGTAGTAGTCCAATAGAATGGAAAACTCGAGATCTAATCTCTCGATTTCCTAGTGTCATGGAAATCAGCCAGGTTCTCTTCATCTGCAACAGTAGAAGAAGAAGAGAGGCTAGCGAGAGAGTCTTATGGCGGAGACGCTAAGGCTTAAATGTAATGTAGATAACCCCTTACGGAACACTTGAGTGCGACGTAGACTACATAATCCCTCAGGGATATTAGCTCTGCTCGATTAACAATAGCATACTTTGTTACACGGAGTGTATCTGGGGGGAATAATACTAACTTACTTAGCACTATCGCGATGCTACGCATTCGCTCTTTCGCTAAATAAGATACGACGATGAGTGGTTGGTGGAGAGAATAACCGATTCTAACTTGATAATTCGCATGAAATAATTTTTTATTTGTTTTTTTTTTTGCTCTTAATTTTAGAGGATGTTTATTTTTATTCTAATAAAAAGGATCCGTTGAA",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[   0,  958,  965,  997, 1006, 1006, 1016, 1018,
                           1032, 1036, 1040, 1049, 1073, 1073, 1074, 1075,
                           1096, 1096, 1106, 1106, 1110, 1116, 1130, 1147,
                           1153, 1153, 1175, 1175, 1186, 1199, 1201, 1203,
                           1214, 1214, 1238, 1238, 1243, 1243, 1244, 1248,
                           1251, 1251, 1256, 1261, 1274, 1275, 1305, 1306,
                           1323, 1323, 1330, 1331, 1348, 1348, 1353, 1354,
                           1364, 1364, 1403, 1403, 1414, 1414, 1415, 1417,
                           1510, 1512, 1526, 1527, 1536, 1536, 1559, 1559,
                           1566, 1566, 1580, 1580, 1598, 1598, 1603, 1610,
                           1615, 1615, 1622, 1622, 1646, 1647, 1717, 1718,
                           1719, 1719, 1723, 1723, 1728, 1728, 1735, 1736,
                           1738, 1738, 1754, 1754, 1771, 1771, 1775, 1775,
                           1809, 1810, 2042, 2043, 2053, 2053, 2057, 2069,
                           2276, 2276, 2279, 2283, 2287, 2287, 2294, 2294,
                           2298, 2299, 2468, 2469, 2983, 2983, 3467, 3468,
                           3507, 3509, 3546],
                          [   0,    0,    7,    7,   16,   22,   32,   32,
                             46,   46,   50,   50,   74,   80,   81,   81,
                            102,  107,  117,  119,  123,  123,  137,  137,
                            143,  147,  169,  170,  181,  181,  183,  183,
                            194,  215,  239,  241,  246,  250,  251,  251,
                            254,  255,  260,  260,  273,  273,  303,  303,
                            320,  322,  329,  329,  346,  348,  353,  353,
                            363,  364,  403,  404,  415,  418,  419,  419,
                            512,  512,  526,  526,  535,  536,  559,  560,
                            567,  568,  582,  584,  602,  606,  611,  611,
                            616,  617,  624,  626,  650,  650,  720,  720,
                            721,  724,  728,  737,  742,  749,  756,  756,
                            758,  761,  777,  778,  795,  796,  800,  801,
                            835,  835, 1067, 1067, 1077, 1078, 1082, 1082,
                           1289, 1290, 1293, 1293, 1297, 1301, 1308, 1310,
                           1314, 1314, 1483, 1483, 1997, 1998, 2482, 2482,
                           2521, 2521, 2558]])
                # fmt: on
            )
        )
        self.assertEqual(
            alignment[0],
            "TATTTTTTGGATTTTTTTCTAGATTTTCTAGGTTATTTAAACCGTTTTTTTTTAATTTAGTGTTTGAGTTTTGACAGGTCTCCACTTTGGGGGCTCCATCGCAAGGAAATTAGAATTCTTATACTTGGTTCTCTTTCCCAGGGACTCCAAGGATCTTTTCATTAGTTTGGATTTTGGTGTTTTCTTTAATTTTGTTAAGAAACAAATCCTTTCTAGAGTTTTTTCTAGCATTATGTTTTTTTTTCTCCTTATCTAAGGGGGTTTGTCGAGGTTTCTTAAATCTTTTTTTCTCTGGGTTTTAAAATTGTTTAAATTTTTTTGACCGAGGGGTTGGGGTGGTTTTCTCATGATAACAGGGGCTGGTGCTTTAGATCCTACCTCTACTGACCCGGGGTCTGCTACTGTGGCTTCTGATGAAGATCCACAGTATGCGCCTACGGAARCTCGGCAGTTTGGTGTTCGAAATCCAGCCCCTCGAATTAATACTCTTGTGCAGGTGGTTGACGAGCGCGGTATCGAATTGCAAAATTTGGGGCGGGACCCCGCTGTTCCGCCTGTTGCTCCGGGGGGGGCAGGTTAATCCTCCAGTCGTCTCCTTTTGGGGGCGTCTTTGACGGGGGTTTAAATCTTTCTTTGGTTGTGGATAGGATTTTTTTTCTAATATCGATCCTACCTGTTTTGGCGGGGCTATTACTTTGTTACTTTTGACCGAAATTTTAATGGAAATTTCTTTGATTCAAATGAATCCCTTAGTTTTCCAACACTTTTTTTTGGTTTTTTTAGGGATAGTCTACGCTGTGGTTAGGTTTGGTTTTATTGGGGGCTTGGTTTGGGCCCACCCCAAATAGGGAGTGGGGGTATGACCTCAGATAGACGAGCTTATTTTAGGGCGGCGACTATAATTATTTCGTTTCCTACAAGGATTAAAGTTTTTTCTTTTACTGTGGGAGGGGGTTTGGTATTAAGAAACGCTAGTCCGGATGTGGCTCTCCATGATACTTATTGT------GTAGTAGCTCATTTTCATTATGTTCTTCGAATGGGAGCAGTCATTGGTATTTTTTTGGTTTTTTTTT------GAAATTTTTAGGTTATTTAGACC-----ATTTTTTTTT--GTTTCGCTAATTAGAATTTTATTAGCCTTTGGTTTTTTTTTATTTTT----TGGGGTTAAGACAAGGTGTCGT-TGAATTAGTTTAGCAAAATACTGCTTAAGGTAGGCTATA---------------------GGATCTACCTTTTATCTTTCTAAT--CTTTT----GTTTTAGT-ATAATTGGTCTTCGATTCAACAATTTTTAGTCTTCAGTCTTTTTTTTTATTTTGAAAAGGTTTTAACACTCT--TGGTTTTGGAGGCTTTGGCTTTCTT--CTTACTCTTAGGAGGA-TGGGCGCTAGAAAGAGTTTTAAGAGGGTGTGAAAGGGGG-TTAATAGCAGG---ATTTGCTTTTTTAACTTATACTGGTTCGTAACGCATTAGCTCAACTCTCTCTTGTAGTTCTAGCAGCCGCCTTTTCTTTGTTGGGGGAGGGTTTAGGAGGAGTCTTTTTTTTCCTAACCCAA-GGTGTTTCTTTCTTTTTTTCTTT-AAAGTTC-TTGACTGTTGGCAC--TTGTCTCCATAAATTTTC----TTTCTTGTAAAGGGCTC-CTAAGGC--TTCTTGTTTCTGAATTCCTCTTTTCTTTTATTCTGTTTTGAGCTTATTTTTCTTGTTAGCTATTACGTAGGCATAGGGCAAATAATTTTTTTTTCTG---CTCT---------CATTA-------TTCCTTCTCC---CTGCTTGTTTCACCCT-GTGGGCTCTTTGAGCCC-CACT-AAGTGAGCGGGGCTCCTGCTTCCGCTCAATTAAATTTTGGTGGGTATTGAGTCTCAGAGGGACTATGATATAGGTTCAGATTGATGGACCTAATCAATCAATTGTATCGCTATACAATCTAGTACCCCTACCAGGGTACCAAGAGAGAGATAACTAGGGTGAATACTACGACTTAGATGTAGTGTTTAAGTTTCTACGGGCTACAGAGAAGCTACCCGCAGGGTAYATATTTGTTCATTACATATTTGTTGACTTTTCTATCTCTGCTTTTACTTTTT-TATTTATTTTTAAATCTTTTTAACTTCAGCTGTTTTTCCTTATCTATTTGACGTAGGCATAGGAAAGTTAACGAATTTTGTAATATTTTTAATTATTTTGTATAGTATACAGGGTAGTGGTATGTAATAGGTAAATTCCATAAGTTCATTATAGTTTATCAGTTGAGAGGAATTTAGTATAAGAAGGCCCATTGGGGCTCTTGTCTTATCCAAGAACTGGTAA-GATTTAATTCT----ACCGGGA--CGGTAGAAATCGGGCAGAGCATGATCTATTTCTTCGGGTATGGCTATAGGGACTAGGTGCTAGGGAGGTATTAGGGCACCGCTCTTTATACAATCTCCATAGATACAACCAGGTCAACTAGGACAACGGAGGACGTTGACAGAGCATAAATAGCGATAGCGTACAAGATAWAATAGGGGCAGTGGTAGCGAAGCGTAGAAGAAAAAATAAGAGTATTGTTTGTAAATAATTCTTTTTTTAGTTTTTAAATATTCTTTTTTTAGGTGGTGTGTGGTTAGGTATGGGGTTAAGGGTGTGGCAAAGAGAAATGTTTATTAAACATTCTTATGGCCGTAGATAGCATATCGATTATACGAGACCTTCGTAAGATCAATCCCCACTAGCATTGCTCATACAGGTTAACTCAATAGGAGGAGCTGGGGTAGAACGTTTCTAGTTCGGGGGTAACCGCAGTTCAATGAAAGTGACGACGTCGGATGGAACAAACTTAATACCACCAGTTGTGCTAACGATTGTTATCTCAATCTATCCCAACAGGCCCCCAGGTAGTGATGAGTGGTGGAATGGTACAGGGTACCAGTGGGTGAAGAGCGTCACGAACCAGGGAATACGGAGTACAGAGTTGAGCGCCCGGGGCTCCGCCCCCGGCTTTTATAGCGCGAGACGTGGTCAGTCGATTCAGCGTTAGG-TTTAAACTCCTTTGGCAAAGATTGACTCTAGCGATCCAGAGACCCTGCCTGGCATAAAAGTCTTTATWAACACCAGTAGGTTCAATAAGGTAGTAATCCAATAGAATGGAAAACTCAAGATCTAATCTCTCGAYTTCCTAGTGTCATGGAAATCAGCCAGGTTCTCTTCATCTGCAACAGTAGAAGAAGAAGAGAGACTAGCGAGAGAGTCTTATGGCGGAGACGCTAAGGCTTAAATGTAATGTAGATAACCCCTTACGGAACACTAGAGTGCGACGTAGACTACATAATCCCTCAGGGATATTAGCTCTGCTCGATTAACAATAGCATACTTTGTTACACGGAGTGTATCTAGGGGGAATAATACTAACTTACTTAGCACTATCGCGATGCTACGCATTCGCTCTTTCGCTAAATAAGATACGACGATGAGTGGTTGGTGGAGAGAATAACCGATTCTAACTTGATAATTCGCATGAAATAATTTTTTTATTTTGTTTTTTTTTTGCTCTTAATTTTAGWGGGRGTGTTTATTTTTATTCTAATAAAAAGGATCCGTTGAA",
        )
        self.assertEqual(
            alignment[1],
            "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TTATTAA--------------------------------TCTTATGGTTTTGCCGTAAAATTTC--TTTCTTTATTTTTT----ATTG---------TTAGGATTTTGTTGATTTTATTTTTCTCAAG-AATTTTTAGGTCAATTAGACCGGCTTATTTTTTTGTCAGTGT------TTAAAGTTTTATTA-----------------ATTTTTGGGGGGGGGGGGAGACGGGGTGTTATCTGAATTAGTTT-------------TT--GGGAGTCTCTAGACATCTCATGGGTTGGCCGGGGGCCTGCCGTCTATAGTTCTTATTCCTTTTAAGGG----AGTAAGAAT-----TTCGATTCAGCAA-CTTTAGTTCACAGTCTTTTTTTTTATTAAG-AAAGGTTTTAATATTCTTGTGGTTTT-GAACCTTTAGGTTTCTTTCTTTAC-CTTCGAGGGATTGGGCACTAGAATGAGTTTTAAGAGTGTGTGAAAGGGGGCTTGATAGCAGGGGAA--TGCTTTTTTAACTTATACTGGCTCGTAACGCATCAGTTCAACTCTCTCTTGCAGTTCTAGCAGCCGCCTTTTTTTTGTTGGGGGGGGGTTAAG--AGAGTGTTTTTTTT-CTAATCCAAGGGTCTTACTTTCTTTCTTTCTTTAAAAATTCTTTGGCTGTCGACACCTTTCTCTCCCGTCAGTCTCATGGTTTCT-------GGCTCTCTTGGGCTTTTTTTGTTTGTGAATGCCTCTTTT-TTTTATTCTGTTTTGAGCTTATTTTTCTTGTTTACTATTACGTAGGTATAGGGCAAATAATTTTTTTTTC-GCGTCTCTTGGCATGCCCATTACTCTAGTTTTATTC-CCGGGCTTCTTCTCTCACCCTAGAGGGCTCTTTGAGCCCACACTCAAGTGAGCGGGGCTCCCGCTTCCGCTCAATTAAA-TTTGGTGGGTATTGAGTCTCAGAGGGACTATGATATAGGTTCAGATTGATGGACCTAGTCAATCAATTGTATCGCTATACAATCTAGTACCCCTACCAGGGTACCAGGAGAGAGATAACTAGGGTGAATACTACGACTTAGATGTACTGTTTAAGTTTCTACGGGCTACAGAGAAGCTACCCGCAGGGTATATATTTGCTCATTACATATTTGTTGATTTTTCTATGTCCGC-TTTACTTTTTATATT------------TTTTTAACTTCAGCTGTTTTTCCTTATCTATTTGACGTAGGCATAGGAAAGTTAACGAATTTTGTAATATTTTTAATTATTTTGTATAGTATACAGGGTAGTGGTATGTAATAGGTAAATTCCATAAGTTCATTATAGTCTATCAGTTGAGAGGAATTTAGTATAAGAAAGCCTGTCAGGGCTCTTGCCTTATCCAAGAACTGGTAAGGAT----TTCTTGACAGAGGGACTCTGT-CAAATCGGGCAGAGCATGATCTATTTCTTCGGGTATGGTTATAAGGCTTAGGTGCTTGGAGGGTATTAGGGCACCGCTCTTAATACAGTCTCCATAGGTGTAACCAGGTCAACTAGGACAACGGAGGACGTTGACAAAGCATGGATAGCGATAGCGTAGAAGATAAAAT-GGGGCAGTGGTAGCGAAGCGTAGAAGAAAAAATAAGAGTATTGTTTGTAAATAATTCTTTTTTTAGTTTTTAAATATTCTTTTTTTAGGTGGTGTGTGGTTAGGTATGGGGTTAGGGGAGTGGCAAAGAGAAGTGTTTATTAAACATTCTTATGGCCGTAGATAGCATATCGATTATACGAGACCTTCGTAAGATCAATCCCCACTAGCATTGCTCATACAGGTTAACTCAATAGGAGGAGCTGGGGTAGAACGTATCTAGTTCGGGGGTAACCGCAGTTCAATGAAAGTGACGACGTCGGATGGAACAAACTTAATACCACCAGTTGTGCTAACGATTGTTATCTCAATCTATCCCAACAGGCCCCCAGGTAGTGATGAGTGGTGGAATGGTACAGGGTACCAGTGGGTGAAGAGCGTCACGAACCAGGGAATACGGAGTACAGAGTTGAGCGCCCGGGGCTCCGCCCCCGGCTTTTATAGCGCGAGACGTGGTCAGTCGATTCAGCGTTAGGTTTTAAACTCCTTTGGCAAAGATTGATTCTAGCGATCCAGAGACCCTGCCTGGCATAAAAGTCTTTATTAGCACCAGTAGGTTCAATAAGGTAGTAGTCCAATAGAATGGAAAACTCGAGATCTAATCTCTCGATTTCCTAGTGTCATGGAAATCAGCCAGGTTCTCTTCATCTGCAACAGTAGAAGAAGAAGAGAGGCTAGCGAGAGAGTCTTATGGCGGAGACGCTAAGGCTTAAATGTAATGTAGATAACCCCTTACGGAACACTTGAGTGCGACGTAGACTACATAATCCCTCAGGGATATTAGCTCTGCTCGATTAACAATAGCATACTTTGTTACACGGAGTGTATCTGGGGGGAATAATACTAACTTACTTAGCACTATCGCGATGCTACGCATTCGCTCTTTCGCTAAATAAGATACGACGATGAGTGGTTGGTGGAGAGAATAACCGATTCTAACTTGATAATTCGCATGAAATAA-TTTTTTATTTGTTTTTTTTTTTGCTCTTAATTTTAGAGG--ATGTTTATTTTTATTCTAATAAAAAGGATCCGTTGAA",
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            "                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              .||||||                                .|||||.||      |||..|..||  ||||.||||.||.|    ||.|         ||.|.|||||.|||.||||.||||      | |||||||||||.|.|||||||     ||||||||.|  ||.|      |||.|.||||||||                 ||||||    .||||...||||..|||||..| |||||||||||             ||  ||.||.||.||                     ||..||.||.|.|||..||||.||  |||||    |    ||| |.|||     |||||||||.||| .||||||...|||||||||||||||||..| |||||||||||.|.|||  ||||||| ||..||||.|.||||||  .|||| |||.|..||| |||||.||||||.||||||||||||.||||||||||||| ||.||||||||   |  |||||||||||||||||||||.|||||||||||.||.||||||||||||||.||||||||||||||||||||.|||||||||||.|||||.||  .||||.|||||||| ||||.|||| |||.||.||||||||.||||||| |||.||| |||.||||.|.|||  ||.|||||....|.|.||    |||||       ||||| ||..|||  ||.||||||.|||||.|||||||| ||||||||||||||||||||||||||||||||..||||||||||||.||||||||||||||||||||||| |   ||||         |||||       ||..||| ||   ||.|||.|.||||||| |.||||||||||||||| |||| ||||||||||||||||.||||||||||||||||| |||||||||||||||||||||||||||||||||||||||||||||||||||||||||.||||||||||||||||||||||||||||||||||||||||||||||||.|||||||||||||||||||||||||||||||||||||||.|||||||||||||||||||||||||||||||||||||||||||.|||||||.||||||||||||||||||.||||||||.||.|| |||||||||| ||||            |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||.|||||||||||||||||||||||||||||.|||..|..|||||||||.||||||||||||||||||| |||    ||||    |..||||  |.|| .|||||||||||||||||||||||||||||||||||||.||||.||..||||||||.||..||||||||||||||||||||.|||||.|||||||||.|..|||||||||||||||||||||||||||||||||||.|||||..||||||||||||||.||||||.||| ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||.|||.|||||||||||||.||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||.|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| |||||||||||||||||||||||||.|||||||||||||||||||||||||||||||||||||||||.|.|||||||||||||||||||||||||.||||||||||||||||||||.||||||||||||||||.||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||.||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||.|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||.|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ||||||||||..||||||||||||||||||||||||.||  .||||||||||||||||||||||||||||||||||||",
        )
        self.assertEqual(
            str(alignment),
            """\
asis              0 TATTTTTTGGATTTTTTTCTAGATTTTCTAGGTTATTTAAACCGTTTTTTTTTAATTTAG
                  0 ------------------------------------------------------------
asis              0 ------------------------------------------------------------

asis             60 TGTTTGAGTTTTGACAGGTCTCCACTTTGGGGGCTCCATCGCAAGGAAATTAGAATTCTT
                 60 ------------------------------------------------------------
asis              0 ------------------------------------------------------------

asis            120 ATACTTGGTTCTCTTTCCCAGGGACTCCAAGGATCTTTTCATTAGTTTGGATTTTGGTGT
                120 ------------------------------------------------------------
asis              0 ------------------------------------------------------------

asis            180 TTTCTTTAATTTTGTTAAGAAACAAATCCTTTCTAGAGTTTTTTCTAGCATTATGTTTTT
                180 ------------------------------------------------------------
asis              0 ------------------------------------------------------------

asis            240 TTTTCTCCTTATCTAAGGGGGTTTGTCGAGGTTTCTTAAATCTTTTTTTCTCTGGGTTTT
                240 ------------------------------------------------------------
asis              0 ------------------------------------------------------------

asis            300 AAAATTGTTTAAATTTTTTTGACCGAGGGGTTGGGGTGGTTTTCTCATGATAACAGGGGC
                300 ------------------------------------------------------------
asis              0 ------------------------------------------------------------

asis            360 TGGTGCTTTAGATCCTACCTCTACTGACCCGGGGTCTGCTACTGTGGCTTCTGATGAAGA
                360 ------------------------------------------------------------
asis              0 ------------------------------------------------------------

asis            420 TCCACAGTATGCGCCTACGGAARCTCGGCAGTTTGGTGTTCGAAATCCAGCCCCTCGAAT
                420 ------------------------------------------------------------
asis              0 ------------------------------------------------------------

asis            480 TAATACTCTTGTGCAGGTGGTTGACGAGCGCGGTATCGAATTGCAAAATTTGGGGCGGGA
                480 ------------------------------------------------------------
asis              0 ------------------------------------------------------------

asis            540 CCCCGCTGTTCCGCCTGTTGCTCCGGGGGGGGCAGGTTAATCCTCCAGTCGTCTCCTTTT
                540 ------------------------------------------------------------
asis              0 ------------------------------------------------------------

asis            600 GGGGGCGTCTTTGACGGGGGTTTAAATCTTTCTTTGGTTGTGGATAGGATTTTTTTTCTA
                600 ------------------------------------------------------------
asis              0 ------------------------------------------------------------

asis            660 ATATCGATCCTACCTGTTTTGGCGGGGCTATTACTTTGTTACTTTTGACCGAAATTTTAA
                660 ------------------------------------------------------------
asis              0 ------------------------------------------------------------

asis            720 TGGAAATTTCTTTGATTCAAATGAATCCCTTAGTTTTCCAACACTTTTTTTTGGTTTTTT
                720 ------------------------------------------------------------
asis              0 ------------------------------------------------------------

asis            780 TAGGGATAGTCTACGCTGTGGTTAGGTTTGGTTTTATTGGGGGCTTGGTTTGGGCCCACC
                780 ------------------------------------------------------------
asis              0 ------------------------------------------------------------

asis            840 CCAAATAGGGAGTGGGGGTATGACCTCAGATAGACGAGCTTATTTTAGGGCGGCGACTAT
                840 ------------------------------------------------------------
asis              0 ------------------------------------------------------------

asis            900 AATTATTTCGTTTCCTACAAGGATTAAAGTTTTTTCTTTTACTGTGGGAGGGGGTTTGGT
                900 ----------------------------------------------------------.|
asis              0 ----------------------------------------------------------TT

asis            960 ATTAAGAAACGCTAGTCCGGATGTGGCTCTCCATGATACTTATTGT------GTAGTAGC
                960 |||||--------------------------------.|||||.||------|||..|..
asis              2 ATTAA--------------------------------TCTTATGGTTTTGCCGTAAAATT

asis           1014 TCATTTTCATTATGTTCTTCGAATGGGAGCAGTCATTGGTATTTTTTTGGTTTTTTTTT-
               1020 ||--||||.||||.||.|----||.|---------||.|.|||||.|||.||||.||||-
asis             30 TC--TTTCTTTATTTTTT----ATTG---------TTAGGATTTTGTTGATTTTATTTTT

asis           1073 -----GAAATTTTTAGGTTATTTAGACC-----ATTTTTTTTT--GTTTCGCTAATTAGA
               1080 -----|-|||||||||||.|.|||||||-----||||||||.|--||.|------|||.|
asis             75 CTCAAG-AATTTTTAGGTCAATTAGACCGGCTTATTTTTTTGTCAGTGT------TTAAA

asis           1121 ATTTTATTAGCCTTTGGTTTTTTTTTATTTTT----TGGGGTTAAGACAAGGTGTCGT-T
               1140 .||||||||-----------------||||||----.||||...||||..|||||..|-|
asis            128 GTTTTATTA-----------------ATTTTTGGGGGGGGGGGGAGACGGGGTGTTATCT

asis           1176 GAATTAGTTTAGCAAAATACTGCTTAAGGTAGGCTATA---------------------G
               1200 ||||||||||-------------||--||.||.||.||---------------------|
asis            171 GAATTAGTTT-------------TT--GGGAGTCTCTAGACATCTCATGGGTTGGCCGGG

asis           1215 GATCTACCTTTTATCTTTCTAAT--CTTTT----GTTTTAGT-ATAATTGGTCTTCGATT
               1260 |..||.||.|.|||..||||.||--|||||----|----|||-|.|||-----|||||||
asis            216 GGCCTGCCGTCTATAGTTCTTATTCCTTTTAAGGG----AGTAAGAAT-----TTCGATT

asis           1268 CAACAATTTTTAGTCTTCAGTCTTTTTTTTTATTTTGAAAAGGTTTTAACACTCT--TGG
               1320 ||.|||-.||||||...|||||||||||||||||..|-|||||||||||.|.|||--|||
asis            267 CAGCAA-CTTTAGTTCACAGTCTTTTTTTTTATTAAG-AAAGGTTTTAATATTCTTGTGG

asis           1326 TTTTGGAGGCTTTGGCTTTCTT--CTTACTCTTAGGAGGA-TGGGCGCTAGAAAGAGTTT
               1380 ||||-||..||||.|.||||||--.||||-|||.|..|||-|||||.||||||.||||||
asis            325 TTTT-GAACCTTTAGGTTTCTTTCTTTAC-CTTCGAGGGATTGGGCACTAGAATGAGTTT

asis           1383 TAAGAGGGTGTGAAAGGGGG-TTAATAGCAGG---ATTTGCTTTTTTAACTTATACTGGT
               1440 ||||||.|||||||||||||-||.||||||||---|--|||||||||||||||||||||.
asis            383 TAAGAGTGTGTGAAAGGGGGCTTGATAGCAGGGGAA--TGCTTTTTTAACTTATACTGGC

asis           1439 TCGTAACGCATTAGCTCAACTCTCTCTTGTAGTTCTAGCAGCCGCCTTTTCTTTGTTGGG
               1500 |||||||||||.||.||||||||||||||.||||||||||||||||||||.|||||||||
asis            441 TCGTAACGCATCAGTTCAACTCTCTCTTGCAGTTCTAGCAGCCGCCTTTTTTTTGTTGGG

asis           1499 GGAGGGTTTAGGAGGAGTCTTTTTTTTCCTAACCCAA-GGTGTTTCTTTCTTTTTTTCTT
               1560 ||.|||||.||--.||||.||||||||-||||.||||-|||.||.||||||||.||||||
asis            501 GGGGGGTTAAG--AGAGTGTTTTTTTT-CTAATCCAAGGGTCTTACTTTCTTTCTTTCTT

asis           1558 T-AAAGTTC-TTGACTGTTGGCAC--TTGTCTCCATAAATTTTC----TTTCTTGTAAAG
               1620 |-|||.|||-|||.||||.|.|||--||.|||||....|.|.||----|||||-------
asis            558 TAAAAATTCTTTGGCTGTCGACACCTTTCTCTCCCGTCAGTCTCATGGTTTCT-------

asis           1610 GGCTC-CTAAGGC--TTCTTGTTTCTGAATTCCTCTTTTCTTTTATTCTGTTTTGAGCTT
               1680 |||||-||..|||--||.||||||.|||||.||||||||-||||||||||||||||||||
asis            611 GGCTCTCTTGGGCTTTTTTTGTTTGTGAATGCCTCTTTT-TTTTATTCTGTTTTGAGCTT

asis           1667 ATTTTTCTTGTTAGCTATTACGTAGGCATAGGGCAAATAATTTTTTTTTCTG---CTCT-
               1740 ||||||||||||..||||||||||||.|||||||||||||||||||||||-|---||||-
asis            670 ATTTTTCTTGTTTACTATTACGTAGGTATAGGGCAAATAATTTTTTTTTC-GCGTCTCTT

asis           1723 --------CATTA-------TTCCTTCTCC---CTGCTTGTTTCACCCT-GTGGGCTCTT
               1800 --------|||||-------||..|||-||---||.|||.|.|||||||-|.||||||||
asis            729 GGCATGCCCATTACTCTAGTTTTATTC-CCGGGCTTCTTCTCTCACCCTAGAGGGCTCTT

asis           1764 TGAGCCC-CACT-AAGTGAGCGGGGCTCCTGCTTCCGCTCAATTAAATTTTGGTGGGTAT
               1860 |||||||-||||-||||||||||||||||.|||||||||||||||||-||||||||||||
asis            788 TGAGCCCACACTCAAGTGAGCGGGGCTCCCGCTTCCGCTCAATTAAA-TTTGGTGGGTAT

asis           1822 TGAGTCTCAGAGGGACTATGATATAGGTTCAGATTGATGGACCTAATCAATCAATTGTAT
               1920 |||||||||||||||||||||||||||||||||||||||||||||.||||||||||||||
asis            847 TGAGTCTCAGAGGGACTATGATATAGGTTCAGATTGATGGACCTAGTCAATCAATTGTAT

asis           1882 CGCTATACAATCTAGTACCCCTACCAGGGTACCAAGAGAGAGATAACTAGGGTGAATACT
               1980 ||||||||||||||||||||||||||||||||||.|||||||||||||||||||||||||
asis            907 CGCTATACAATCTAGTACCCCTACCAGGGTACCAGGAGAGAGATAACTAGGGTGAATACT

asis           1942 ACGACTTAGATGTAGTGTTTAAGTTTCTACGGGCTACAGAGAAGCTACCCGCAGGGTAYA
               2040 ||||||||||||||.|||||||||||||||||||||||||||||||||||||||||||.|
asis            967 ACGACTTAGATGTACTGTTTAAGTTTCTACGGGCTACAGAGAAGCTACCCGCAGGGTATA

asis           2002 TATTTGTTCATTACATATTTGTTGACTTTTCTATCTCTGCTTTTACTTTTT-TATTTATT
               2100 ||||||.||||||||||||||||||.||||||||.||.||-||||||||||-||||----
asis           1027 TATTTGCTCATTACATATTTGTTGATTTTTCTATGTCCGC-TTTACTTTTTATATT----

asis           2061 TTTAAATCTTTTTAACTTCAGCTGTTTTTCCTTATCTATTTGACGTAGGCATAGGAAAGT
               2160 --------||||||||||||||||||||||||||||||||||||||||||||||||||||
asis           1082 --------TTTTTAACTTCAGCTGTTTTTCCTTATCTATTTGACGTAGGCATAGGAAAGT

asis           2121 TAACGAATTTTGTAATATTTTTAATTATTTTGTATAGTATACAGGGTAGTGGTATGTAAT
               2220 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
asis           1134 TAACGAATTTTGTAATATTTTTAATTATTTTGTATAGTATACAGGGTAGTGGTATGTAAT

asis           2181 AGGTAAATTCCATAAGTTCATTATAGTTTATCAGTTGAGAGGAATTTAGTATAAGAAGGC
               2280 |||||||||||||||||||||||||||.|||||||||||||||||||||||||||||.||
asis           1194 AGGTAAATTCCATAAGTTCATTATAGTCTATCAGTTGAGAGGAATTTAGTATAAGAAAGC

asis           2241 CCATTGGGGCTCTTGTCTTATCCAAGAACTGGTAA-GATTTAATTCT----ACCGGGA--
               2340 |..|..|||||||||.|||||||||||||||||||-|||----||||----|..||||--
asis           1254 CTGTCAGGGCTCTTGCCTTATCCAAGAACTGGTAAGGAT----TTCTTGACAGAGGGACT

asis           2294 CGGTAGAAATCGGGCAGAGCATGATCTATTTCTTCGGGTATGGCTATAGGGACTAGGTGC
               2400 |.||-.|||||||||||||||||||||||||||||||||||||.||||.||..|||||||
asis           1310 CTGT-CAAATCGGGCAGAGCATGATCTATTTCTTCGGGTATGGTTATAAGGCTTAGGTGC

asis           2354 TAGGGAGGTATTAGGGCACCGCTCTTTATACAATCTCCATAGATACAACCAGGTCAACTA
               2460 |.||..||||||||||||||||||||.|||||.|||||||||.|..||||||||||||||
asis           1369 TTGGAGGGTATTAGGGCACCGCTCTTAATACAGTCTCCATAGGTGTAACCAGGTCAACTA

asis           2414 GGACAACGGAGGACGTTGACAGAGCATAAATAGCGATAGCGTACAAGATAWAATAGGGGC
               2520 |||||||||||||||||||||.|||||..||||||||||||||.||||||.|||-|||||
asis           1429 GGACAACGGAGGACGTTGACAAAGCATGGATAGCGATAGCGTAGAAGATAAAAT-GGGGC

asis           2474 AGTGGTAGCGAAGCGTAGAAGAAAAAATAAGAGTATTGTTTGTAAATAATTCTTTTTTTA
               2580 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
asis           1488 AGTGGTAGCGAAGCGTAGAAGAAAAAATAAGAGTATTGTTTGTAAATAATTCTTTTTTTA

asis           2534 GTTTTTAAATATTCTTTTTTTAGGTGGTGTGTGGTTAGGTATGGGGTTAAGGGTGTGGCA
               2640 |||||||||||||||||||||||||||||||||||||||||||||||||.|||.||||||
asis           1548 GTTTTTAAATATTCTTTTTTTAGGTGGTGTGTGGTTAGGTATGGGGTTAGGGGAGTGGCA

asis           2594 AAGAGAAATGTTTATTAAACATTCTTATGGCCGTAGATAGCATATCGATTATACGAGACC
               2700 |||||||.||||||||||||||||||||||||||||||||||||||||||||||||||||
asis           1608 AAGAGAAGTGTTTATTAAACATTCTTATGGCCGTAGATAGCATATCGATTATACGAGACC

asis           2654 TTCGTAAGATCAATCCCCACTAGCATTGCTCATACAGGTTAACTCAATAGGAGGAGCTGG
               2760 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
asis           1668 TTCGTAAGATCAATCCCCACTAGCATTGCTCATACAGGTTAACTCAATAGGAGGAGCTGG

asis           2714 GGTAGAACGTTTCTAGTTCGGGGGTAACCGCAGTTCAATGAAAGTGACGACGTCGGATGG
               2820 ||||||||||.|||||||||||||||||||||||||||||||||||||||||||||||||
asis           1728 GGTAGAACGTATCTAGTTCGGGGGTAACCGCAGTTCAATGAAAGTGACGACGTCGGATGG

asis           2774 AACAAACTTAATACCACCAGTTGTGCTAACGATTGTTATCTCAATCTATCCCAACAGGCC
               2880 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
asis           1788 AACAAACTTAATACCACCAGTTGTGCTAACGATTGTTATCTCAATCTATCCCAACAGGCC

asis           2834 CCCAGGTAGTGATGAGTGGTGGAATGGTACAGGGTACCAGTGGGTGAAGAGCGTCACGAA
               2940 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
asis           1848 CCCAGGTAGTGATGAGTGGTGGAATGGTACAGGGTACCAGTGGGTGAAGAGCGTCACGAA

asis           2894 CCAGGGAATACGGAGTACAGAGTTGAGCGCCCGGGGCTCCGCCCCCGGCTTTTATAGCGC
               3000 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
asis           1908 CCAGGGAATACGGAGTACAGAGTTGAGCGCCCGGGGCTCCGCCCCCGGCTTTTATAGCGC

asis           2954 GAGACGTGGTCAGTCGATTCAGCGTTAGG-TTTAAACTCCTTTGGCAAAGATTGACTCTA
               3060 |||||||||||||||||||||||||||||-|||||||||||||||||||||||||.||||
asis           1968 GAGACGTGGTCAGTCGATTCAGCGTTAGGTTTTAAACTCCTTTGGCAAAGATTGATTCTA

asis           3013 GCGATCCAGAGACCCTGCCTGGCATAAAAGTCTTTATWAACACCAGTAGGTTCAATAAGG
               3120 |||||||||||||||||||||||||||||||||||||.|.||||||||||||||||||||
asis           2028 GCGATCCAGAGACCCTGCCTGGCATAAAAGTCTTTATTAGCACCAGTAGGTTCAATAAGG

asis           3073 TAGTAATCCAATAGAATGGAAAACTCAAGATCTAATCTCTCGAYTTCCTAGTGTCATGGA
               3180 |||||.||||||||||||||||||||.||||||||||||||||.||||||||||||||||
asis           2088 TAGTAGTCCAATAGAATGGAAAACTCGAGATCTAATCTCTCGATTTCCTAGTGTCATGGA

asis           3133 AATCAGCCAGGTTCTCTTCATCTGCAACAGTAGAAGAAGAAGAGAGACTAGCGAGAGAGT
               3240 ||||||||||||||||||||||||||||||||||||||||||||||.|||||||||||||
asis           2148 AATCAGCCAGGTTCTCTTCATCTGCAACAGTAGAAGAAGAAGAGAGGCTAGCGAGAGAGT

asis           3193 CTTATGGCGGAGACGCTAAGGCTTAAATGTAATGTAGATAACCCCTTACGGAACACTAGA
               3300 |||||||||||||||||||||||||||||||||||||||||||||||||||||||||.||
asis           2208 CTTATGGCGGAGACGCTAAGGCTTAAATGTAATGTAGATAACCCCTTACGGAACACTTGA

asis           3253 GTGCGACGTAGACTACATAATCCCTCAGGGATATTAGCTCTGCTCGATTAACAATAGCAT
               3360 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
asis           2268 GTGCGACGTAGACTACATAATCCCTCAGGGATATTAGCTCTGCTCGATTAACAATAGCAT

asis           3313 ACTTTGTTACACGGAGTGTATCTAGGGGGAATAATACTAACTTACTTAGCACTATCGCGA
               3420 |||||||||||||||||||||||.||||||||||||||||||||||||||||||||||||
asis           2328 ACTTTGTTACACGGAGTGTATCTGGGGGGAATAATACTAACTTACTTAGCACTATCGCGA

asis           3373 TGCTACGCATTCGCTCTTTCGCTAAATAAGATACGACGATGAGTGGTTGGTGGAGAGAAT
               3480 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
asis           2388 TGCTACGCATTCGCTCTTTCGCTAAATAAGATACGACGATGAGTGGTTGGTGGAGAGAAT

asis           3433 AACCGATTCTAACTTGATAATTCGCATGAAATAATTTTTTTATTTTGTTTTTTTTTTGCT
               3540 ||||||||||||||||||||||||||||||||||-||||||||||..|||||||||||||
asis           2448 AACCGATTCTAACTTGATAATTCGCATGAAATAA-TTTTTTATTTGTTTTTTTTTTTGCT

asis           3493 CTTAATTTTAGWGGGRGTGTTTATTTTTATTCTAATAAAAAGGATCCGTTGAA 3546
               3600 |||||||||||.||--.|||||||||||||||||||||||||||||||||||| 3653
asis           2507 CTTAATTTTAGAGG--ATGTTTATTTTTATTCTAATAAAAAGGATCCGTTGAA 2558
""",
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_water_reverse1(self):
        # water -asequence seqA.fa -bsequence seqB.fa -gapopen 10 -gapextend 0.5 -sreverse1 -outfile water_reverse1.txt
        path = "Emboss/water_reverse1.txt"
        alignments = Align.parse(path, "emboss")
        self.assertEqual(alignments.metadata["Program"], "water")
        self.assertEqual(alignments.metadata["Rundate"], "Sat 22 Oct 2022 23:47:41")
        self.assertEqual(
            alignments.metadata["Command line"],
            "water -asequence seqA.fa -bsequence seqB.fa -gapopen 0.001 -gapextend 0.001 -sreverse1 -outfile water_reverse1.txt",
        )
        self.assertEqual(alignments.metadata["Align_format"], "srspair")
        self.assertEqual(alignments.metadata["Report_file"], "water_reverse1.txt")

        alignment = next(alignments)
        self.assertEqual(alignment.annotations["Matrix"], "EDNAFULL")
        self.assertAlmostEqual(alignment.annotations["Gap_penalty"], 0.001)
        self.assertAlmostEqual(alignment.annotations["Extend_penalty"], 0.001)
        self.assertEqual(alignment.annotations["Identity"], 32)
        self.assertEqual(alignment.annotations["Similarity"], 32)
        self.assertEqual(alignment.annotations["Gaps"], 89)
        self.assertAlmostEqual(alignment.annotations["Score"], 159.911)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 121))
        self.assertEqual(alignment.sequences[0].id, "seqA")
        self.assertEqual(alignment.sequences[1].id, "seqB")
        self.assertEqual(
            alignment.sequences[0].seq,
            "GGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCC",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "GGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCC",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[121, 102, 13,  0],
                          [  0,  19, 19, 32]])
                # fmt: on
            )
        )
        self.assertEqual(
            alignment[0],
            "GGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCC",
        )
        self.assertEqual(
            alignment[1],
            "GGGGGGGGGGGGGGGGGGG-----------------------------------------------------------------------------------------CCCCCCCCCCCCC",
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            "|||||||||||||||||||                                                                                         |||||||||||||",
        )
        self.assertEqual(
            str(alignment),
            """\
seqA            121 GGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
                  0 |||||||||||||||||||-----------------------------------------
seqB              0 GGGGGGGGGGGGGGGGGGG-----------------------------------------

seqA             61 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCCCCCCCCCCCC
                 60 ------------------------------------------------||||||||||||
seqB             19 ------------------------------------------------CCCCCCCCCCCC

seqA              1 C   0
                120 | 121
seqB             31 C  32
""",
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_water_reverse2(self):
        # water -asequence seqA.fa -bsequence seqB.fa -gapopen 10 -gapextend 0.5 -sreverse2 -outfile water_reverse2.txt
        path = "Emboss/water_reverse2.txt"
        alignments = Align.parse(path, "emboss")
        self.assertEqual(alignments.metadata["Program"], "water")
        self.assertEqual(alignments.metadata["Rundate"], "Sun 23 Oct 2022 00:06:18")
        self.assertEqual(
            alignments.metadata["Command line"],
            "water -asequence seqA.fa -bsequence seqB.fa -gapopen 0.001 -gapextend 0.001 -sreverse2 -outfile water_reverse2.txt",
        )
        self.assertEqual(alignments.metadata["Align_format"], "srspair")
        self.assertEqual(alignments.metadata["Report_file"], "water_reverse2.txt")
        alignment = next(alignments)
        self.assertEqual(alignment.annotations["Matrix"], "EDNAFULL")
        self.assertAlmostEqual(alignment.annotations["Gap_penalty"], 0.001)
        self.assertAlmostEqual(alignment.annotations["Extend_penalty"], 0.001)
        self.assertEqual(alignment.annotations["Identity"], 32)
        self.assertEqual(alignment.annotations["Similarity"], 32)
        self.assertEqual(alignment.annotations["Gaps"], 89)
        self.assertAlmostEqual(alignment.annotations["Score"], 159.911)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 121))
        self.assertEqual(alignment.sequences[0].id, "seqA")
        self.assertEqual(alignment.sequences[1].id, "seqB")
        self.assertEqual(
            alignment.sequences[0].seq,
            "GGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCC",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "GGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCC",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[ 0, 13, 102, 121],
                          [32, 19,  19,   0]])
                # fmt: on
            )
        )
        self.assertEqual(
            alignment[0],
            "GGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCC",
        )
        self.assertEqual(
            alignment[1],
            "GGGGGGGGGGGGG-----------------------------------------------------------------------------------------CCCCCCCCCCCCCCCCCCC",
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            "|||||||||||||                                                                                         |||||||||||||||||||",
        )
        self.assertEqual(
            str(alignment),
            """\
seqA              0 GGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                  0 |||||||||||||-----------------------------------------------
seqB             32 GGGGGGGGGGGGG-----------------------------------------------

seqA             60 AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCC
                 60 ------------------------------------------||||||||||||||||||
seqB             19 ------------------------------------------CCCCCCCCCCCCCCCCCC

seqA            120 C 121
                120 | 121
seqB              1 C   0
""",
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_water_reverse3(self):
        # water -asequence seqA.fa -bsequence seqB.fa -gapopen 10 -gapextend 0.5 -sreverse1 -outfile water_reverse3.txt
        path = "Emboss/water_reverse3.txt"
        alignments = Align.parse(path, "emboss")
        self.assertEqual(alignments.metadata["Program"], "water")
        self.assertEqual(alignments.metadata["Rundate"], "Sat 22 Oct 2022 22:56:03")
        self.assertEqual(
            alignments.metadata["Command line"],
            "water -asequence seqA.fa -bsequence seqB.fa -gapopen 1 -gapextend 0.5 -sreverse1 -outfile water_reverse3.txt",
        )
        self.assertEqual(alignments.metadata["Align_format"], "srspair")
        self.assertEqual(alignments.metadata["Report_file"], "water_reverse3.txt")

        alignment = next(alignments)
        self.assertEqual(alignment.annotations["Matrix"], "EDNAFULL")
        self.assertAlmostEqual(alignment.annotations["Gap_penalty"], 1.0)
        self.assertAlmostEqual(alignment.annotations["Extend_penalty"], 0.5)
        self.assertEqual(alignment.annotations["Identity"], 16)
        self.assertEqual(alignment.annotations["Similarity"], 16)
        self.assertEqual(alignment.annotations["Gaps"], 3)
        self.assertAlmostEqual(alignment.annotations["Score"], 77.5)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 19))
        self.assertEqual(alignment.sequences[0].id, "seqA")
        self.assertEqual(alignment.sequences[1].id, "seqB")
        self.assertEqual(
            repr(alignment.sequences[0].seq),
            "Seq({2: 'GGGCCCGGTTTAAAAAAA'}, length=20)",
        )
        self.assertEqual(
            repr(alignment.sequences[1].seq),
            "Seq({2: 'TTTTTTTACCCGGGCCC'}, length=19)",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[20, 13, 11, 10, 10,  2],
                          [ 2,  9,  9, 10, 11, 19]])
                # fmt: on
            )
        )
        self.assertEqual(alignment[0], "TTTTTTTAAA-CCGGGCCC")
        self.assertEqual(alignment[1], "TTTTTTT--ACCCGGGCCC")
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            "|||||||  | ||||||||",
        )
        self.assertEqual(
            str(alignment),
            """\
seqA             20 TTTTTTTAAA-CCGGGCCC  2
                  0 |||||||--|-|||||||| 19
seqB              2 TTTTTTT--ACCCGGGCCC 19
""",
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_water_reverse4(self):
        # water -asequence seqA.fa -bsequence seqB.fa -gapopen 10 -gapextend 0.5 -sreverse2 -outfile water_reverse4.txt
        path = "Emboss/water_reverse4.txt"
        alignments = Align.parse(path, "emboss")
        self.assertEqual(alignments.metadata["Program"], "water")
        self.assertEqual(alignments.metadata["Rundate"], "Sat 22 Oct 2022 22:56:15")
        self.assertEqual(
            alignments.metadata["Command line"],
            "water -asequence seqA.fa -bsequence seqB.fa -gapopen 1 -gapextend 0.5 -sreverse2 -outfile water_reverse4.txt",
        )
        self.assertEqual(alignments.metadata["Align_format"], "srspair")
        self.assertEqual(alignments.metadata["Report_file"], "water_reverse4.txt")
        alignment = next(alignments)
        self.assertEqual(alignment.annotations["Matrix"], "EDNAFULL")
        self.assertAlmostEqual(alignment.annotations["Gap_penalty"], 1.0)
        self.assertAlmostEqual(alignment.annotations["Extend_penalty"], 0.5)
        self.assertEqual(alignment.annotations["Identity"], 16)
        self.assertEqual(alignment.annotations["Similarity"], 16)
        self.assertEqual(alignment.annotations["Gaps"], 3)
        self.assertAlmostEqual(alignment.annotations["Score"], 77.5)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 19))
        self.assertEqual(alignment.sequences[0].id, "seqA")
        self.assertEqual(alignment.sequences[1].id, "seqB")
        self.assertEqual(
            repr(alignment.sequences[0].seq),
            "Seq({2: 'GGGCCCGGTTTAAAAAAA'}, length=20)",
        )
        self.assertEqual(
            repr(alignment.sequences[1].seq),
            "Seq({2: 'TTTTTTTACCCGGGCCC'}, length=19)",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[ 2, 10, 12, 12, 20],
                          [19, 11, 11, 10,  2]])
                # fmt: on
            )
        )
        self.assertEqual(alignment[0], "GGGCCCGGTT-TAAAAAAA")
        self.assertEqual(alignment[1], "GGGCCCGG--GTAAAAAAA")
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            "||||||||   ||||||||",
        )
        self.assertEqual(
            str(alignment),
            """\
seqA              2 GGGCCCGGTT-TAAAAAAA 20
                  0 ||||||||---|||||||| 19
seqB             19 GGGCCCGG--GTAAAAAAA  2
""",
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_pair_aln_full_blank_line(self):
        path = "Emboss/emboss_pair_aln_full_blank_line.txt"
        alignments = Align.parse(path, "emboss")
        self.assertEqual(alignments.metadata["Program"], "stretcher")
        self.assertEqual(alignments.metadata["Rundate"], "Tue 15 May 2018 17:01:31")
        self.assertEqual(
            alignments.metadata["Command line"],
            "stretcher -auto -stdout -asequence emboss_stretcher-I20180515-170128-0371-22292969-p1m.aupfile -bsequence emboss_stretcher-I20180515-170128-0371-22292969-p1m.bupfile -datafile EDNAFULL -gapopen 16 -gapextend 4 -aformat3 pair -snucleotide1 -snucleotide2",
        )
        self.assertEqual(alignments.metadata["Align_format"], "pair")
        self.assertEqual(alignments.metadata["Report_file"], "stdout")
        alignment = next(alignments)
        self.assertEqual(alignment.annotations["Matrix"], "EDNAFULL")
        self.assertAlmostEqual(alignment.annotations["Gap_penalty"], 16)
        self.assertAlmostEqual(alignment.annotations["Extend_penalty"], 4)
        self.assertEqual(alignment.annotations["Identity"], 441)
        self.assertEqual(alignment.annotations["Similarity"], 441)
        self.assertEqual(alignment.annotations["Gaps"], 847)
        self.assertAlmostEqual(alignment.annotations["Score"], -2623)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 1450))
        self.assertEqual(
            alignment.sequences[0].id, "hg38_chrX_131691529_131830643_47210_48660"
        )
        self.assertEqual(
            alignment.sequences[1].id, "mm10_chrX_50555743_50635321_27140_27743"
        )
        self.assertEqual(
            alignment.sequences[0].seq,
            "GGCAGGTGCATAGCTTGAGCCTAGGAGTTCAAGTCCAGCCCTGACAATGTAGAGAGACCCCGTCTCTTCAAAAAATACAAAAAATAGCCAGGCATGGTGACCTACAATGGAAGCCCTAGCTACGTAGGAGGCGGAAATGGGAGGATCACCTCAGCCCAGGGAGGCTGATGTTGCAGTGAGCCATGATCATGCCTCTACACTCCACCCTGGGCAACAGAGTAAGATGCTGTCTAAAATATATATATATGCATATCTGTGTGTATATATATATATATATATGTGTGTGTGTGTGTGTGTATATACATATGTGTGTGTATATACATATATGTGTATATATATATGTGTGTATATATACATATACATATTCAGCATCACCTTATATTCTTTGAATATATCTACATCAATACATACTTTTGAGTGCTTGAAATTTTTTATATTTTACTCTAGAAGAACTGTAAGAAATTATAAAGTAGAAAACTTGTGGTAGGTCAAACATAGTAAGAAGAAATAATCACTTTTTAAAGGTCTGTGCTAGGTACTATGATCTGTTCCCTATATATACATAATATGGACTTTTATAAACTAATGTTCAAATTCCCCTGTAGTATAACTTCTTGTTGTTGTTTATTTTTTTTTTTTTGTATTTTTCATTTTAGATATGGGGTTTCACTCTGTTGACCAGGCTGATCTCGAACCACTGGTCTCAAGCGATCCTCCCATCTTGGACTCCCAAAGTGCTAGGATTACAGGCACGAGGCACCTTGACTGGCCACCATGTACTATAGCTGTTAAAACAAGTTTGTTTCACTGATAACTGGAGTACTTTTCAAATATAATTAATAATTCATGGAAATAATGATAGCTTTAAAAGTATTGGCACTTTTAAAAACTGAGTTTGTAAACTTCATATAACATAAAATTAACCATTAAAATGTATTAATTTCAATGGCATTTAGGACACTCACAATGCAGTGCAAGCATTACCACTATGTAGTGGCAAATCATTTTCACTACCACAAAAGAAAATCCTGGACCCATTAGTTAGTCATTCCCCATTCCACTCTCTGCCCAGCCCCTGGCAAACACTCATCTGATTTCCCTCACTACTGATCATCACAACAAGTGGCCTTGTTCATCTTGTTGTGGGAACCAGGAGACCAGAGAGACCAATGGGTGGAACAGGAGGATTTTACTAGGTGGTCACCGACTCAGCAGATTAACATCCAAAGGCTGAGCCCCAAACCAAGAGAGGGCTTGACTTTTATACATATATCTGAAAAGGGCCCAAAACCTGTAAGGCCGGTAAGCAAGCTTACAGCAGAACAAAGGCAGTTTATCAAACAGTGACAGGTTTTACAGTTCAGGCATGTCTTGTGACCTTTGCCATAACTGCACAGCTGGAAAACAGGAACTTACAAAATCCTTACAAGCTTGCAGAAACAGTTACAAA",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "GTTCAAGGCCATCCGGGATTAAAGGTGTGGTAGAACTCTTCTGATGGAGACAATATAAGGACATTGGAAGAAGGGAGTCTTGCCCTTGCTCCTTCGCCTACTTGCTGTGTAAGACTGAGTAACTCCTAGACCCTTGGACTTCCATTTCAGCCACTACTGAACCATTGTTGGGAATTGGGCTGCAGACTGTAAGTCATCAATAAATTCCTTTACTATATAGAGACTATCCATAAATTCTGTGACTCTAGAGAACCCTGACAATACAACTGGGAAGCACGGACATCCTCTTTGAGATATAATTATCAACTGGCAAGTGTTTGTTTATTGATATTTTACTTAAGACAAAGTTAAACCTACTCCTGTCCTCTGGGCATGGTAGCATGGACTTATTCTGGAACTACCAGAGGAAAAGACAGAAGCCTACTGGAAAGGCCCAGGCCATCCTGCCTCTTGTAGTTCACTAGGACCAGGGCTCAGCATAGTCCTTGGCTTCTAAATCTGCTACCATATCTTTATCATGTAAAACTGACACAAAATTAAACATATCAAAATTTTATGAAAACCATTAAGTATCTGGAAAAGAAAAAAATCAACAGTTATAAA",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[   0,    1,   27,   45,   49,   54,   61,   68,
                             79,   91,   93,  101,  105,  113,  133,  145,
                            157,  163,  168,  174,  180,  199,  224,  236,
                            253,  259,  370,  377,  378,  391,  408,  420,
                            440,  466,  469,  479,  483,  497,  509,  518,
                            537,  543,  576,  585,  595,  598,  606,  616,
                            647,  659,  668,  673,  676,  682,  693, 701,
                            753,  760,  767,  778,  825,  839,  844,  847,
                            854,  858,  861,  863,  868,  877,  881,  885,
                            895,  898,  899,  927,  939,  944,  946,  957,
                            962,  964,  972,  983,  985,  988,  989,  997,
                           1008, 1019, 1021, 1036, 1042, 1044, 1064, 1080,
                           1087, 1091, 1106, 1109, 1133, 1144, 1189, 1197,
                           1210, 1215, 1228, 1236, 1247, 1263, 1269, 1279,
                           1285, 1297, 1299, 1305, 1308, 1327, 1329, 1331,
                           1334, 1341, 1357, 1368, 1376, 1381, 1386, 1391,
                           1395, 1407, 1414, 1422, 1438, 1450],
                          [   0,    1,    1,   19,   19,   24,   24,   31,
                             31,   43,   43,   51,   51,   59,   59,   71,
                             71,   77,   77,   83,   83,  102,  102,  114,
                            114,  120,  120,  127,  127,  140,  140,  152,
                            152,  178,  178,  188,  188,  202,  202,  211,
                            211,  217,  217,  226,  226,  229,  229,  239,
                            239,  251,  251,  256,  256,  262,  262,  270,
                            270,  277,  277,  288,  288,  302,  302,  305,
                            305,  309,  309,  311,  311,  320,  320,  324,
                            324,  327,  327,  355,  355,  360,  360,  371,
                            371,  373,  373,  384,  384,  387,  387,  395,
                            395,  406,  406,  421,  421,  423,  423,  439,
                            439,  443,  443,  446,  446,  457,  457,  465,
                            465,  470,  470,  478,  478,  494,  494,  504,
                            504,  516,  516,  522,  522,  541,  541,  543,
                            543,  550,  550,  561,  561,  566,  566,  571,
                            571,  583,  583,  591,  591,  603]])
                # fmt: on
            )
        )
        self.assertEqual(
            alignment[0],
            "GGCAGGTGCATAGCTTGAGCCTAGGAGTTCAAGTCCAGCCCTGACAATGTAGAGAGACCCCGTCTCTTCAAAAAATACAAAAAATAGCCAGGCATGGTGACCTACAATGGAAGCCCTAGCTACGTAGGAGGCGGAAATGGGAGGATCACCTCAGCCCAGGGAGGCTGATGTTGCAGTGAGCCATGATCATGCCTCTACACTCCACCCTGGGCAACAGAGTAAGATGCTGTCTAAAATATATATATATGCATATCTGTGTGTATATATATATATATATATGTGTGTGTGTGTGTGTGTATATACATATGTGTGTGTATATACATATATGTGTATATATATATGTGTGTATATATACATATACATATTCAGCATCACCTTATATTCTTTGAATATATCTACATCAATACATACTTTTGAGTGCTTGAAATTTTTTATATTTTACTCTAGAAGAACTGTAAGAAATTATAAAGTAGAAAACTTGTGGTAGGTCAAACATAGTAAGAAGAAATAATCACTTTTTAAAGGTCTGTGCTAGGTACTATGATCTGTTCCCTATATATACATAATATGGACTTTTATAAACTAATGTTCAAATTCCCCTGTAGTATAACTTCTTGTTGTTGTTTATTTTTTTTTTTTTGTATTTTTCATTTTAGATATGGGGTTTCACTCTGTTGACCAGGCTGATCTCGAACCACTGGTCTCAAGCGATCCTCCCATCTTGGACTCCCAAAGTGCTAGGATTACAGGCACGAGGCACCTTGACTGGCCACCATGTACTATAGCTGTTAAAACAAGTTTGTTTCACTGATAACTGGAGTACTTTTCAAATATAATTAATAATTCATGGAAATAATGATAGCTTTAAAAGTATTGGCACTTTTAAAAACTGAGTTTGTAAACTTCATATAACATAAAATTAACCATTAAAATGTATTAATTTCAATGGCATTTAGGACACTCACAATGCAGTGCAAGCATTACCACTATGTAGTGGCAAATCATTTTCACTACCACAAAAGAAAATCCTGGACCCATTAGTTAGTCATTCCCCATTCCACTCTCTGCCCAGCCCCTGGCAAACACTCATCTGATTTCCCTCACTACTGATCATCACAACAAGTGGCCTTGTTCATCTTGTTGTGGGAACCAGGAGACCAGAGAGACCAATGGGTGGAACAGGAGGATTTTACTAGGTGGTCACCGACTCAGCAGATTAACATCCAAAGGCTGAGCCCCAAACCAAGAGAGGGCTTGACTTTTATACATATATCTGAAAAGGGCCCAAAACCTGTAAGGCCGGTAAGCAAGCTTACAGCAGAACAAAGGCAGTTTATCAAACAGTGACAGGTTTTACAGTTCAGGCATGTCTTGTGACCTTTGCCATAACTGCACAGCTGGAAAACAGGAACTTACAAAATCCTTACAAGCTTGCAGAAACAGTTACAAA",
        )
        self.assertEqual(
            alignment[1],
            "G--------------------------TTCAAGGCCATCCGGGAT----TAAAG-------GTGTGGT-----------AGAACTCTTCTG--ATGGAGAC----AATATAAG--------------------GACATTGGAAGA------------AGGGAG-----TCTTGC------CCTTGCTCCTTCGCCTACT-------------------------TGCTGTGTAAGA-----------------CTGAGT---------------------------------------------------------------------------------------------------------------AACTCCT-AGACCCTTGGACT-----------------TCCATTTCAGCC--------------------ACTACTGAACCATTGTTGGGAATTGG---GCTGCAGACT----GTAAGTCATCAATA------------AATTCCTTT-------------------ACTATA---------------------------------TAGAGACTA----------TCC--------ATAAATTCTG-------------------------------TGACTCTAGAGA---------ACCCT---GACAAT-----------ACAACTGG----------------------------------------------------GAAGCAC-------GGACATCCTCT-----------------------------------------------TTGAGATATAATTA-----TCA-------ACTG---GC-----AAGTGTTTG----TTTA----------TTG-ATATTTTACTTAAGACAAAGTTAAACCT------------ACTCC--TGTCCTCTGGG-----CA--------TGGTAGCATGG--ACT-TATTCTGG-----------AACTACCAGAG--GAAAAGACAGAAGCC------TA--------------------CTGGAAAGGCCCAGGC-------CATC---------------CTG------------------------CCTCTTGTAGT---------------------------------------------TCACTAGG-------------ACCAG-------------GGCTCAGC-----------ATAGTCCTTGGCTTCT------AAATCTGCTA------CCATATCTTTAT--CATGTA---AAACTGACACAAAATTAAA--CA---TATCAAA----------------ATTTTATGAAA--------ACCAT-----TAAGT----ATCTGGAAAAGA-------AAAAAATC----------------AACAGTTATAAA",
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            "|                          ||||||.|||.||..||.    ||.||       ||.|..|           |.||.|...|.|  ||||.|||    |||..|||                    ||.||.|||.||            ||||||     |.||||      ||.||.||.|.|..||||.                         ||||||.|||.|                 |||.||                                                                                                               |.|.||| |.|..|||.||.|                 |.|.|||.||..                    |||...|||..|.|||..|.||||..   |..|.|.|||    |||.||||...|||            |||..||||                   |||||.                                 ||.|.||||          |||        ||||.||||.                               |.|.|.||||.|         ||.||   |||.|.           ||.|||||                                                    ||.||||       ||.||.|.|.|                                               ||.|.|||||||||     |||       |.||   ||     ||||.||.|    ||||          ||| |.|.||.|..|||.|.|||.||||.|.|            |.|.|  ||.|.|.|.||     ||        ||..|||||..  ||| |.|..|||           .|||||||.|.  |||||..|.|.|.||      ||                    |||...||.|||.|||       ||||               |||                        |.||||||.||                                             |.||||||             |.|||             ||||.|||           |.||..||||.|||.|      |.|||||..|      |.|.|.||.||.  |..|||   ||.||.|||..|.|..|||  ||   |||||||                |.||.|.|.|.        |||.|     |||.|    |.||||||||.|       |.||||||                ||||||||.|||",
        )
        self.assertEqual(
            str(alignment),
            """\
hg38_chrX         0 GGCAGGTGCATAGCTTGAGCCTAGGAGTTCAAGTCCAGCCCTGACAATGTAGAGAGACCC
                  0 |--------------------------||||||.|||.||..||.----||.||------
mm10_chrX         0 G--------------------------TTCAAGGCCATCCGGGAT----TAAAG------

hg38_chrX        60 CGTCTCTTCAAAAAATACAAAAAATAGCCAGGCATGGTGACCTACAATGGAAGCCCTAGC
                 60 -||.|..|-----------|.||.|...|.|--||||.|||----|||..|||-------
mm10_chrX        24 -GTGTGGT-----------AGAACTCTTCTG--ATGGAGAC----AATATAAG-------

hg38_chrX       120 TACGTAGGAGGCGGAAATGGGAGGATCACCTCAGCCCAGGGAGGCTGATGTTGCAGTGAG
                120 -------------||.||.|||.||------------||||||-----|.||||------
mm10_chrX        59 -------------GACATTGGAAGA------------AGGGAG-----TCTTGC------

hg38_chrX       180 CCATGATCATGCCTCTACACTCCACCCTGGGCAACAGAGTAAGATGCTGTCTAAAATATA
                180 ||.||.||.|.|..||||.-------------------------||||||.|||.|----
mm10_chrX        83 CCTTGCTCCTTCGCCTACT-------------------------TGCTGTGTAAGA----

hg38_chrX       240 TATATATGCATATCTGTGTGTATATATATATATATATATGTGTGTGTGTGTGTGTGTATA
                240 -------------|||.||-----------------------------------------
mm10_chrX       114 -------------CTGAGT-----------------------------------------

hg38_chrX       300 TACATATGTGTGTGTATATACATATATGTGTATATATATATGTGTGTATATATACATATA
                300 ------------------------------------------------------------
mm10_chrX       120 ------------------------------------------------------------

hg38_chrX       360 CATATTCAGCATCACCTTATATTCTTTGAATATATCTACATCAATACATACTTTTGAGTG
                360 ----------|.|.|||-|.|..|||.||.|-----------------|.|.|||.||..
mm10_chrX       120 ----------AACTCCT-AGACCCTTGGACT-----------------TCCATTTCAGCC

hg38_chrX       420 CTTGAAATTTTTTATATTTTACTCTAGAAGAACTGTAAGAAATTATAAAGTAGAAAACTT
                420 --------------------|||...|||..|.|||..|.||||..---|..|.|.|||-
mm10_chrX       152 --------------------ACTACTGAACCATTGTTGGGAATTGG---GCTGCAGACT-

hg38_chrX       480 GTGGTAGGTCAAACATAGTAAGAAGAAATAATCACTTTTTAAAGGTCTGTGCTAGGTACT
                480 ---|||.||||...|||------------|||..||||-------------------|||
mm10_chrX       188 ---GTAAGTCATCAATA------------AATTCCTTT-------------------ACT

hg38_chrX       540 ATGATCTGTTCCCTATATATACATAATATGGACTTTTATAAACTAATGTTCAAATTCCCC
                540 ||.---------------------------------||.|.||||----------|||--
mm10_chrX       214 ATA---------------------------------TAGAGACTA----------TCC--

hg38_chrX       600 TGTAGTATAACTTCTTGTTGTTGTTTATTTTTTTTTTTTTGTATTTTTCATTTTAGATAT
                600 ------||||.||||.-------------------------------|.|.|.||||.|-
mm10_chrX       229 ------ATAAATTCTG-------------------------------TGACTCTAGAGA-

hg38_chrX       660 GGGGTTTCACTCTGTTGACCAGGCTGATCTCGAACCACTGGTCTCAAGCGATCCTCCCAT
                660 --------||.||---|||.|.-----------||.|||||-------------------
mm10_chrX       251 --------ACCCT---GACAAT-----------ACAACTGG-------------------

hg38_chrX       720 CTTGGACTCCCAAAGTGCTAGGATTACAGGCACGAGGCACCTTGACTGGCCACCATGTAC
                720 ---------------------------------||.||||-------||.||.|.|.|--
mm10_chrX       270 ---------------------------------GAAGCAC-------GGACATCCTCT--

hg38_chrX       780 TATAGCTGTTAAAACAAGTTTGTTTCACTGATAACTGGAGTACTTTTCAAATATAATTAA
                780 ---------------------------------------------||.|.|||||||||-
mm10_chrX       288 ---------------------------------------------TTGAGATATAATTA-

hg38_chrX       840 TAATTCATGGAAATAATGATAGCTTTAAAAGTATTGGCACTTTTAAAAACTGAGTTTGTA
                840 ----|||-------|.||---||-----||||.||.|----||||----------|||-|
mm10_chrX       302 ----TCA-------ACTG---GC-----AAGTGTTTG----TTTA----------TTG-A

hg38_chrX       900 AACTTCATATAACATAAAATTAACCATTAAAATGTATTAATTTCAATGGCATTTAGGACA
                900 .|.||.|..|||.|.|||.||||.|.|------------|.|.|--||.|.|.|.||---
mm10_chrX       328 TATTTTACTTAAGACAAAGTTAAACCT------------ACTCC--TGTCCTCTGGG---

hg38_chrX       960 CTCACAATGCAGTGCAAGCATTACCACTATGTAGTGGCAAATCATTTTCACTACCACAAA
                960 --||--------||..|||||..--|||-|.|..|||-----------.|||||||.|.-
mm10_chrX       371 --CA--------TGGTAGCATGG--ACT-TATTCTGG-----------AACTACCAGAG-

hg38_chrX      1020 AGAAAATCCTGGACCCATTAGTTAGTCATTCCCCATTCCACTCTCTGCCCAGCCCCTGGC
               1020 -|||||..|.|.|.||------||--------------------|||...||.|||.|||
mm10_chrX       406 -GAAAAGACAGAAGCC------TA--------------------CTGGAAAGGCCCAGGC

hg38_chrX      1080 AAACACTCATCTGATTTCCCTCACTACTGATCATCACAACAAGTGGCCTTGTTCATCTTG
               1080 -------||||---------------|||------------------------|.|||||
mm10_chrX       439 -------CATC---------------CTG------------------------CCTCTTG

hg38_chrX      1140 TTGTGGGAACCAGGAGACCAGAGAGACCAATGGGTGGAACAGGAGGATTTTACTAGGTGG
               1140 |.||---------------------------------------------|.||||||---
mm10_chrX       453 TAGT---------------------------------------------TCACTAGG---

hg38_chrX      1200 TCACCGACTCAGCAGATTAACATCCAAAGGCTGAGCCCCAAACCAAGAGAGGGCTTGACT
               1200 ----------|.|||-------------||||.|||-----------|.||..||||.||
mm10_chrX       465 ----------ACCAG-------------GGCTCAGC-----------ATAGTCCTTGGCT

hg38_chrX      1260 TTTATACATATATCTGAAAAGGGCCCAAAACCTGTAAGGCCGGTAAGCAAGCTTACAGCA
               1260 |.|------|.|||||..|------|.|.|.||.||.--|..|||---||.||.|||..|
mm10_chrX       491 TCT------AAATCTGCTA------CCATATCTTTAT--CATGTA---AAACTGACACAA

hg38_chrX      1320 GAACAAAGGCAGTTTATCAAACAGTGACAGGTTTTACAGTTCAGGCATGTCTTGTGACCT
               1320 .|..|||--||---|||||||----------------|.||.|.|.|.--------|||.
mm10_chrX       534 AATTAAA--CA---TATCAAA----------------ATTTTATGAAA--------ACCA

hg38_chrX      1380 TTGCCATAACTGCACAGCTGGAAAACAGGAACTTACAAAATCCTTACAAGCTTGCAGAAA
               1380 |-----|||.|----|.||||||||.|-------|.||||||----------------||
mm10_chrX       565 T-----TAAGT----ATCTGGAAAAGA-------AAAAAATC----------------AA

hg38_chrX      1440 CAGTTACAAA 1450
               1440 ||||||.||| 1450
mm10_chrX       593 CAGTTATAAA  603
""",
        )
        with self.assertRaises(StopIteration):
            next(alignments)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
