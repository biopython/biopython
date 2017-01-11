#!/usr/bin/env python

# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
import unittest

from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62


class TestPairwiseErrorConditions(unittest.TestCase):
    """Test several error conditions"""

    def test_function_name(self):
        """Test for wrong function names"""
        # Function name pattern must be globalXX or localXX
        self.assertRaises(AttributeError, lambda: pairwise2.align.globalxxx)
        self.assertRaises(AttributeError, lambda: pairwise2.align.localxxx)
        self.assertRaises(AttributeError, lambda: pairwise2.align.glocalxx)
        # First X must be from (x, m, d, c), second from (x, s, d, c)
        self.assertRaises(AttributeError, lambda: pairwise2.align.globalax)
        self.assertRaises(AttributeError, lambda: pairwise2.align.globalxa)

    def test_function_parameters(self):
        """Test for number of parameteres"""
        # globalxx takes two parameters
        self.assertRaises(TypeError, pairwise2.align.globalxx, 'A')
        # matrix_only is no keyword argument
        self.assertRaises(TypeError, pairwise2.align.globalxx, 'A', 'C',
                          {'matrix_only': True})
        # Both sequences must be either strings or lists
        self.assertRaises(TypeError, pairwise2.align.globalxx, 'A', ['C'])
        # If both sequences are lists, gap_char must also be set as list
        self.assertRaises(TypeError, pairwise2.align.globalxx, ['A'], ['C'])

        # If one or both sequences are empty, there is no alignment
        alignment = pairwise2.align.globalxx('A', '')
        self.assertEqual(alignment, [])

        # Gap scores must be negativ
        self.assertRaises(ValueError, pairwise2.align.globalxs, 'A', 'C',
                          5, -1)
        self.assertRaises(ValueError, pairwise2.align.globalxs, 'A', 'C',
                          -5, 1)
        # Gap open penalty must be higher than gap extension penalty
        self.assertRaises(ValueError, pairwise2.align.globalxs, 'A', 'C',
                          -1, -5)


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

    def test_globalxx_simple2(self):
        """Do the same test with sequence order reversed"""
        aligns = pairwise2.align.globalxx("GAT", "GAACT")
        self.assertEqual(len(aligns), 2)
        aligns.sort()
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(alignment, """\
G-A-T
|||||
GAACT
  Score=3
""")
        seq1, seq2, score, begin, end = aligns[1]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(alignment, """\
GA--T
|||||
GAACT
  Score=3
""")


class TestPairwiseLocal(unittest.TestCase):

    def test_localxs(self):
        aligns = sorted(pairwise2.align.localxs("AxBx", "zABz", -0.1, 0))
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

    def test_localms(self):
        """Two different local alignments"""
        aligns = sorted(pairwise2.align.localms("xxxABCDxxx", "zzzABzzCDz", 1,
                                                -0.5, -3, -1))
        alignment = pairwise2.format_alignment(*aligns[0])
        self.assertEqual(alignment, """\
--xxxABCDxxx
       ||
zzzABzzCDz--
  Score=2
""")
        alignment = pairwise2.format_alignment(*aligns[1])
        self.assertEqual(alignment, """\
xxxABCDxxx
   ||
zzzABzzCDz
  Score=2
""")

    def test_blosum62(self):
        """Test localds with blosum62."""
        self.assertEqual(1, blosum62[('K', 'Q')])
        self.assertEqual(4, blosum62[('A', 'A')])
        self.assertEqual(8, blosum62[('H', 'H')])
        alignments = pairwise2.align.localds('VKAHGKKV', 'FQAHCAGV',
                                             blosum62, -4, -4)
        for a in alignments:
            self.assertEqual(pairwise2.format_alignment(*a),
                             "VKAHGKKV\n |||\nFQAHCAGV\n  Score=13\n")


class TestScoreOnly(unittest.TestCase):
    """Test paramater ``score_only``"""

    def test_score_only_global(self):
        """Test ``score_only`` in a global alignment"""
        aligns1 = pairwise2.align.globalxx("GAACT", "GAT")
        aligns2 = pairwise2.align.globalxx("GAACT", "GAT", score_only=True)
        self.assertEqual(aligns1[0][2], aligns2)

    def test_score_only_local(self):
        """Test ``score_only`` in a local alignment"""
        aligns1 = pairwise2.align.localms("xxxABCDxxx", "zzzABzzCDz", 1, -0.5,
                                          -3, -1)
        aligns2 = pairwise2.align.localms("xxxABCDxxx", "zzzABzzCDz", 1, -0.5,
                                          -3, -1, score_only=True)
        self.assertEqual(aligns1[0][2], aligns2)


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
GC-T-
|||||
G-ATA
  Score=1.7
""")


class TestPairwiseExtendPenalty(unittest.TestCase):

    def test_extend_penalty1(self):
        aligns = pairwise2.align.globalxs("GACT", "GT", -0.5, -0.2)
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
        aligns = pairwise2.align.globalxs("GACT", "GT", -1.5, -0.2)
        self.assertEqual(len(aligns), 1)
        aligns.sort()
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(alignment, """\
GACT
||||
G--T
  Score=0.3
""")


class TestPairwisePenalizeExtendWhenOpening(unittest.TestCase):

    def test_penalize_extend_when_opening(self):
        aligns = pairwise2.align.globalxs("GACT", "GT", -0.2, -1.5,
                                          penalize_extend_when_opening=1)
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
        aligns = pairwise2.align.globalxs("GACT", "GT", -0.8, -0.2,
                                          penalize_end_gaps=0)
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

    def test_penalize_end_gaps2(self):
        """Do the same, but use the generic method (with the same resutlt)"""
        aligns = pairwise2.align.globalxs("GACT", "GT", -0.8, -0.2,
                                          penalize_end_gaps=0,
                                          force_generic=True)
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
        aligns = pairwise2.align.localxd(list("GAAT"), list("GTCCT"),
                                         -0.1, 0, -0.1, -0.1, gap_char=["-"])
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
        ("A", "A"): 1.5,
        ("A", "T"): 0.5,
        ("T", "T"): 1.0
        }

    def test_match_dictionary1(self):
        aligns = pairwise2.align.localds("ATAT", "ATT", self.match_dict,
                                         -.5, 0)
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


class TestPersiteGapPenalties(unittest.TestCase):
    """Check gap penalty callbacks use correct gap opening position.

    This tests that the gap penalty callbacks are really being used
    with the correct gap opening position.
    """

    def test_gap_here_only_1(self):
        seq1 = "AAAABBBAAAACCCCCCCCCCCCCCAAAABBBAAAA"
        seq2 = "AABBBAAAACCCCAAAABBBAA"

        def no_gaps(x, y):
            """Very expensive to open a gap in seq1."""

            x = 0  # fool QuantifiedCode, x is not used here
            return -2000 - y

        def specific_gaps(x, y):
            """Very expensive to open a gap in seq2

            ...unless it is in one of the allowed positions:
            """
            breaks = [0, 11, len(seq2)]
            return (-2 - y) if x in breaks else (-2000 - y)

        alignments = pairwise2.align.globalmc(seq1, seq2, 1, -1, no_gaps,
                                              specific_gaps)
        self.assertEqual(len(alignments), 1)
        formatted = pairwise2.format_alignment(*alignments[0])
        self.assertEqual(formatted, """\
AAAABBBAAAACCCCCCCCCCCCCCAAAABBBAAAA
||||||||||||||||||||||||||||||||||||
--AABBBAAAACC----------CCAAAABBBAA--
  Score=2
""")

    def test_gap_here_only_2(self):
        """Force a bad alignment.

        Forces a bad alignment by having a very expensive gap penalty
        where one would normally expect a gap, and a cheap gap penalty
        in another place.
        """
        seq1 = "AAAABBBAAAACCCCCCCCCCCCCCAAAABBBAAAA"
        seq2 = "AABBBAAAACCCCAAAABBBAA"

        def no_gaps(x, y):
            """Very expensive to open a gap in seq1."""
            x = 0  # fool QuantifiedCode, x is not used here
            return -2000 - y

        def specific_gaps(x, y):
            """Very expensive to open a gap in seq2

            ...unless it is in one of the allowed positions:
            """
            breaks = [0, 3, len(seq2)]
            return (-2 - y) if x in breaks else (-2000 - y)

        alignments = pairwise2.align.globalmc(seq1, seq2, 1, -1, no_gaps,
                                              specific_gaps)
        self.assertEqual(len(alignments), 1)
        formatted = pairwise2.format_alignment(*alignments[0])
        self.assertEqual(formatted, """\
AAAABBBAAAACCCCCCCCCCCCCCAAAABBBAAAA
||||||||||||||||||||||||||||||||||||
--AAB----------BBAAAACCCCAAAABBBAA--
  Score=-10
""")


class TestOtherFunctions(unittest.TestCase):
    """Test remaining non-tested private methods."""

    def test_clean_alignments(self):
        """``_clean_alignments`` removes redundant alignments."""
        alns = [
            ('ACCGT', 'AC-G-', 3.0, 0, 4),
            ('ACCGT', 'A-CG-', 3.0, 0, 4),
            ('ACCGT', 'AC-G-', 3.0, 0, 4),
            ('ACCGT', 'A-CG-', 3.0, 0, 4),
        ]
        expected = [
            ('ACCGT', 'AC-G-', 3.0, 0, 4),
            ('ACCGT', 'A-CG-', 3.0, 0, 4),
        ]
        result = pairwise2._clean_alignments(alns)
        self.assertEqual(expected, result)

    def test_print_matrix(self):
        """``print_matrix`` prints nested lists as nice matrices."""
        import sys

        try:  # Python 2
            from StringIO import StringIO
        except ImportError:  # Python 3
            from io import StringIO
        out = StringIO()
        sys.stdout = out
        pairwise2.print_matrix([[0.0, -1.0, -1.5, -2.0], [-1.0, 4.0, 3.0, 2.5],
                                [-1.5, 3.0, 8.0, 7.0], [-2.0, 2.5, 7.0, 6.0],
                                [-2.5, 2.0, 6.5, 11.0],
                                [-3.0, 1.5, 6.0, 10.0]])
        self.assertEqual(out.getvalue(),
                         ' 0.0  -1.0  -1.5  -2.0 \n'
                         '-1.0   4.0   3.0   2.5 \n'
                         '-1.5   3.0   8.0   7.0 \n'
                         '-2.0   2.5   7.0   6.0 \n'
                         '-2.5   2.0   6.5  11.0 \n'
                         '-3.0   1.5   6.0  10.0 \n')
        sys.stdout = sys.__stdout__


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
