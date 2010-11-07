#!/usr/bin/env python

# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import unittest

from Bio import pairwise2


class TestPairwiseGlobal(unittest.TestCase):

    def test_globalxx_simple(self):
        aligns = pairwise2.align.globalxx("GAACT", "GAT")
        self.assertEqual(len(aligns), 2)
        aligns.sort()
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(alignment, """\
GAACT
|||||
G-A-T
  Score=3
""")
        seq1, seq2, score, begin, end = aligns[1]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(alignment, """\
GAACT
|||||
GA--T
  Score=3
""")


class TestPairwiseLocal(unittest.TestCase):

    def test_localxs(self):
        aligns = pairwise2.align.localxs("AxBx", "zABz", -0.1, 0)
        aligns.sort()
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(alignment, """\
-AxBx
 |||
zA-Bz
  Score=1.9
""")
        seq1, seq2, score, begin, end = aligns[1]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(alignment, """\
-AxBx
 ||||
zA-Bz
  Score=1.9
""")


class TestPairwiseOpenPenalty(unittest.TestCase):

    def test_match_score_open_penalty1(self):
        aligns = pairwise2.align.globalms("AA", "A", 2.0, -1, -0.1, 0)
        self.assertEqual(len(aligns), 2)
        aligns.sort()
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(alignment, """\
AA
||
-A
  Score=1.9
""")
        seq1, seq2, score, begin, end = aligns[1]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(alignment, """\
AA
||
A-
  Score=1.9
""")

    def test_match_score_open_penalty2(self):
        aligns = pairwise2.align.globalms("GAA", "GA", 1.5, 0, -0.1, 0)
        self.assertEqual(len(aligns), 2)
        aligns.sort()
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(alignment, """\
GAA
|||
G-A
  Score=2.9
""")
        seq1, seq2, score, begin, end = aligns[1]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(alignment, """\
GAA
|||
GA-
  Score=2.9
""")

    def test_match_score_open_penalty3(self):
        aligns = pairwise2.align.globalxs("GAACT", "GAT", -0.1, 0)
        self.assertEqual(len(aligns), 1)
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(alignment, """\
GAACT
|||||
GA--T
  Score=2.9
""")

    def test_match_score_open_penalty4(self):
        aligns = pairwise2.align.globalms("GCT", "GATA", 1, -2, -0.1, 0)
        self.assertEqual(len(aligns), 1)
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(alignment, """\
GCT-
||||
GATA
  Score=-0.1
""")


class TestPairwiseExtendPenalty(unittest.TestCase):

    def test_extend_penalty1(self):
        aligns = pairwise2.align.globalxs("GACT", "GT", -0.2, -0.5)
        self.assertEqual(len(aligns), 1)
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(alignment, """\
GACT
||||
G--T
  Score=1.3
""")

    def test_extend_penalty2(self):
        aligns = pairwise2.align.globalxs("GACT", "GT", -0.2, -1.5)
        self.assertEqual(len(aligns), 2)
        aligns.sort()
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(alignment, """\
GACT
||||
-G-T
  Score=0.6
""")
        seq1, seq2, score, begin, end = aligns[1]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(alignment, """\
GACT
||||
G-T-
  Score=0.6
""")


class TestPairwisePenalizeExtendWhenOpening(unittest.TestCase):

    def test_penalize_extend_when_opening(self):
        aligns = pairwise2.align.globalxs("GACT", "GT", -0.2, -1.5, penalize_extend_when_opening=1)
        self.assertEqual(len(aligns), 1)
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(alignment, """\
GACT
||||
G--T
  Score=-1.2
""")


class TestPairwisePenalizeEndgaps(unittest.TestCase):

    def test_penalize_end_gaps(self):
        aligns = pairwise2.align.globalxs("GACT", "GT", -0.2, -0.8, penalize_end_gaps=0)
        self.assertEqual(len(aligns), 3)
        aligns.sort()
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(alignment, """\
GACT
||||
--GT
  Score=1
""")
        seq1, seq2, score, begin, end = aligns[1]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(alignment, """\
GACT
||||
G--T
  Score=1
""")
        seq1, seq2, score, begin, end = aligns[2]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(alignment, """\
GACT
||||
GT--
  Score=1
""")


class TestPairwiseSeparateGapPenalties(unittest.TestCase):

    def test_separate_gap_penalties1(self):
        aligns = pairwise2.align.localxd("GAT", "GTCT", -0.3, 0, -0.8, 0)
        self.assertEqual(len(aligns), 2)
        aligns.sort()
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(alignment, """\
G-AT
||||
GTCT
  Score=1.7
""")
        seq1, seq2, score, begin, end = aligns[1]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(alignment, """\
GA-T
||||
GTCT
  Score=1.7
""")

    def test_separate_gap_penalties2(self):
        aligns = pairwise2.align.localxd("GAT", "GTCT", -0.5, 0, -0.2, 0)
        self.assertEqual(len(aligns), 1)
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(alignment, """\
GAT--
|||
G-TCT
  Score=1.8
""")



class TestPairwiseSeparateGapPenaltiesWithExtension(unittest.TestCase):

    def test_separate_gap_penalties_with_extension(self):
        aligns = pairwise2.align.localxd(list("GAAT"), list("GTCCT"), -0.1, 0, -0.1, -0.1, gap_char=["-"])
        self.assertEqual(len(aligns), 3)
        aligns.sort()
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(alignment, """\
['G', '-', 'A', 'A', 'T']
|||||
['G', 'T', 'C', 'C', 'T']
  Score=1.9
""")
        seq1, seq2, score, begin, end = aligns[1]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(alignment, """\
['G', 'A', '-', 'A', 'T']
|||||
['G', 'T', 'C', 'C', 'T']
  Score=1.9
""")
        seq1, seq2, score, begin, end = aligns[2]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(alignment, """\
['G', 'A', 'A', '-', 'T']
|||||
['G', 'T', 'C', 'C', 'T']
  Score=1.9
""")


class TestPairwiseMatchDictionary(unittest.TestCase):

    match_dict = {
        ("A", "A") : 1.5,
        ("A", "T") : 0.5,
        ("T", "T") : 1.0
        }

    def test_match_dictionary1(self):
        aligns = pairwise2.align.localds("ATAT", "ATT", self.match_dict, -.5, 0)
        self.assertEqual(len(aligns), 2)
        aligns.sort()
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(alignment, """\
ATAT
||||
AT-T
  Score=3
""")
        seq1, seq2, score, begin, end = aligns[1]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(alignment, """\
ATAT
|||
ATT-
  Score=3
""")

    def test_match_dictionary2(self):
        aligns = pairwise2.align.localds("ATAT", "ATT", self.match_dict, -1, 0)
        self.assertEqual(len(aligns), 1)
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(alignment, """\
ATAT
|||
ATT-
  Score=3
""")

    def test_match_dictionary3(self):
        aligns = pairwise2.align.localds("ATT", "ATAT", self.match_dict, -1, 0)
        self.assertEqual(len(aligns), 1)
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(alignment, """\
ATT-
|||
ATAT
  Score=3
""")


class TestPairwiseOneCharacter(unittest.TestCase):

    def test_align_one_char1(self):
        aligns = pairwise2.align.localxs("abcde", "c", -0.3, -0.1)
        self.assertEqual(len(aligns), 1)
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(alignment, """\
abcde
  |
--c--
  Score=1
""")

    def test_align_one_char2(self):
        aligns = pairwise2.align.localxs("abcce", "c", -0.3, -0.1)
        self.assertEqual(len(aligns), 2)
        aligns.sort()
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(alignment, """\
abcce
   |
---c-
  Score=1
""")
        seq1, seq2, score, begin, end = aligns[1]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(alignment, """\
abcce
  |
--c--
  Score=1
""")

    def test_align_one_char3(self):
        aligns = pairwise2.align.globalxs("abcde", "c", -0.3, -0.1)
        self.assertEqual(len(aligns), 1)
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(alignment, """\
abcde
|||||
--c--
  Score=0.2
""")


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
