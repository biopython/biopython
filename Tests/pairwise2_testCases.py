#!/usr/bin/env python

# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for pairwise2 module.

Put new test case here, the classes here will be imported and run
as TestCases in ``test_pairwise2.py`` and ``test_pairwise2_no_C.py``
with or without complementing C extensions.

"""
import pickle
import unittest
import warnings

from Bio import BiopythonWarning
from Bio import pairwise2
from Bio.Align import substitution_matrices


class TestPairwiseErrorConditions(unittest.TestCase):
    """Test several error conditions."""

    def test_function_name(self):
        """Test for wrong function names."""
        # Function name pattern must be globalXX or localXX
        self.assertRaises(AttributeError, lambda: pairwise2.align.globalxxx)
        self.assertRaises(AttributeError, lambda: pairwise2.align.localxxx)
        self.assertRaises(AttributeError, lambda: pairwise2.align.glocalxx)
        # First X must be from (x, m, d, c), second from (x, s, d, c)
        self.assertRaises(AttributeError, lambda: pairwise2.align.globalax)
        self.assertRaises(AttributeError, lambda: pairwise2.align.globalxa)

    def test_function_parameters(self):
        """Test for number of parameters."""
        # globalxx takes two parameters
        self.assertRaises(TypeError, pairwise2.align.globalxx, "A")
        # matrix_only is no keyword argument
        self.assertRaises(
            TypeError, pairwise2.align.globalxx, "A", "C", {"matrix_only": True}
        )
        # Both sequences must be either strings or lists
        self.assertRaises(TypeError, pairwise2.align.globalxx, "A", ["C"])
        # If both sequences are lists, gap_char must also be set as list
        self.assertRaises(TypeError, pairwise2.align.globalxx, ["A"], ["C"])

        # If one or both sequences are empty, there is no alignment
        alignment = pairwise2.align.globalxx("A", "")
        self.assertEqual(alignment, [])

        # Gap scores must be negative
        self.assertRaises(ValueError, pairwise2.align.globalxs, "A", "C", 5, -1)
        self.assertRaises(ValueError, pairwise2.align.globalxs, "A", "C", -5, 1)
        # Gap open penalty must be higher than gap extension penalty
        self.assertRaises(ValueError, pairwise2.align.globalxs, "A", "C", -1, -5)

    def test_param_names(self):
        """Test for unknown parameter in parameter names."""
        a = pairwise2.align.alignment_function("globalxx")
        a.param_names = ["Hello"]
        self.assertRaises(ValueError, a.decode, "Bye")

    def test_warnings(self):
        """Test for warnings."""
        with warnings.catch_warnings(record=True) as w:
            # Cause all warnings to always be triggered.
            warnings.simplefilter("always")
            # Trigger a warning.
            pairwise2.align.localxx("GA", "CGA", penalize_end_gaps=True)
            # Verify some things
            self.assertEqual(len(w), 1)
            self.assertEqual(w[-1].category, BiopythonWarning)
            self.assertIn("should not", str(w[-1].message))


class TestPairwiseKeywordUsage(unittest.TestCase):
    """Tests for keyword usage."""

    def test_keywords(self):
        """Test equality of calls with and without keywords."""
        aligns = pairwise2.align.globalxx("GAACT", "GAT")
        aligns_kw = pairwise2.align.globalxx(sequenceA="GAACT", sequenceB="GAT")
        self.assertEqual(aligns, aligns_kw)

        aligns = pairwise2.align.globalmx("GAACT", "GAT", 5, -4)
        aligns_kw = pairwise2.align.globalmx(
            sequenceA="GAACT", sequenceB="GAT", match=5, mismatch=-4
        )
        self.assertEqual(aligns, aligns_kw)


class TestPairwiseGlobal(unittest.TestCase):
    """Test some usual global alignments."""

    def test_globalxx_simple(self):
        """Test globalxx."""
        aligns = pairwise2.align.globalxx("GAACT", "GAT")
        self.assertEqual(len(aligns), 2)
        aligns.sort()
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(
            alignment,
            """\
GAACT
| | |
G-A-T
  Score=3
""",
        )
        seq1, seq2, score, begin, end = aligns[1]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(
            alignment,
            """\
GAACT
||  |
GA--T
  Score=3
""",
        )

    def test_globalxx_simple2(self):
        """Do the same test with sequence order reversed."""
        aligns = pairwise2.align.globalxx("GAT", "GAACT")
        self.assertEqual(len(aligns), 2)
        aligns.sort()
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(
            alignment,
            """\
G-A-T
| | |
GAACT
  Score=3
""",
        )
        seq1, seq2, score, begin, end = aligns[1]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(
            alignment,
            """\
GA--T
||  |
GAACT
  Score=3
""",
        )

    def test_one_alignment_only(self):
        """Test one_alignment_only parameter."""
        aligns = pairwise2.align.globalxx("ACCGT", "ACG")
        self.assertEqual(len(aligns), 2)
        aligns = pairwise2.align.globalxx("ACCGT", "ACG", one_alignment_only=True)
        self.assertEqual(len(aligns), 1)

    def test_list_input(self):
        """Do a global alignment with sequences supplied as lists."""
        aligns = pairwise2.align.globalxx(
            ["Gly", "Ala", "Thr"], ["Gly", "Ala", "Ala", "Cys", "Thr"], gap_char=["---"]
        )
        aligns.sort()
        seq1, seq2, score, begin, end = aligns[0]
        self.assertEqual(score, 3)
        self.assertEqual(seq1, ["Gly", "---", "Ala", "---", "Thr"])
        self.assertEqual(seq2, ["Gly", "Ala", "Ala", "Cys", "Thr"])


class TestPairwiseLocal(unittest.TestCase):
    """Test some simple local alignments."""

    def setUp(self):
        self.blosum62 = substitution_matrices.load("BLOSUM62")

    def test_localxs_1(self):
        """Test localxx."""
        aligns = sorted(pairwise2.align.localxs("AxBx", "zABz", -0.1, 0))
        # From Biopython 1.74 on this should only give one alignment, since
        # we disallow leading and trailing 'zero-extensions'
        self.assertEqual(len(aligns), 1)
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(
            alignment,
            """\
1 AxB
  | |
2 A-B
  Score=1.9
""",
        )

    def test_localxs_2(self):
        """Test localxx with ``full_sequences=True``."""
        aligns = sorted(pairwise2.align.localxs("AxBx", "zABz", -0.1, 0))
        # From Biopython 1.74 on this should only give one alignment, since
        # we disallow leading and trailing 'zero-extensions'
        self.assertEqual(len(aligns), 1)
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(
            seq1, seq2, score, begin, end, full_sequences=True
        )
        self.assertEqual(
            alignment,
            """\
-AxBx
 | | 
zA-Bz
  Score=1.9
""",  # noqa: W291
        )

    def test_localds_zero_score_segments_symmetric(self):
        """Test if alignment is independent on direction of sequence."""
        aligns1 = pairwise2.align.localds(
            "CWHISLKM", "CWHGISGLKM", self.blosum62, -11, -1
        )
        aligns2 = pairwise2.align.localds(
            "MKLSIHWC", "MKLGSIGHWC", self.blosum62, -11, -1
        )
        self.assertEqual(len(aligns1), len(aligns2))

    def test_localxs_generic(self):
        """Test the generic method with local alignments."""
        aligns = sorted(
            pairwise2.align.localxs("AxBx", "zABz", -0.1, 0, force_generic=True)
        )
        # From Biopython 1.74 on this should only give one alignment, since
        # we disallow leading and trailing 'zero-extensions'
        self.assertEqual(len(aligns), 1)
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(
            alignment,
            """\
1 AxB
  | |
2 A-B
  Score=1.9
""",
        )

    def test_localms(self):
        """Two different local alignments."""
        aligns = sorted(
            pairwise2.align.localms("xxxABCDxxx", "zzzABzzCDz", 1, -0.5, -3, -1)
        )
        alignment = pairwise2.format_alignment(*aligns[0])
        self.assertEqual(
            alignment,
            """\
6 CD
  ||
8 CD
  Score=2
""",
        )
        alignment = pairwise2.format_alignment(*aligns[1])
        self.assertEqual(
            alignment,
            """\
4 AB
  ||
4 AB
  Score=2
""",
        )

    def test_blosum62(self):
        """Test localds with blosum62."""
        self.assertEqual(1, self.blosum62[("K", "Q")])
        self.assertEqual(4, self.blosum62[("A", "A")])
        self.assertEqual(8, self.blosum62[("H", "H")])
        alignments = pairwise2.align.localds(
            "VKAHGKKV", "FQAHCAGV", self.blosum62, -4, -4
        )
        for a in alignments:
            self.assertEqual(
                pairwise2.format_alignment(*a), "2 KAH\n  .||\n2 QAH\n  Score=13\n"
            )

    def test_empty_result(self):
        """Return no alignment."""
        self.assertEqual(pairwise2.align.localxx("AT", "GC"), [])


class TestScoreOnly(unittest.TestCase):
    """Test parameter ``score_only``."""

    def test_score_only_global(self):
        """Test ``score_only`` in a global alignment."""
        aligns1 = pairwise2.align.globalxx("GAACT", "GAT")
        aligns2 = pairwise2.align.globalxx("GAACT", "GAT", score_only=True)
        self.assertEqual(aligns1[0][2], aligns2)

    def test_score_only_local(self):
        """Test ``score_only`` in a local alignment."""
        aligns1 = pairwise2.align.localms("xxxABCDxxx", "zzzABzzCDz", 1, -0.5, -3, -1)
        aligns2 = pairwise2.align.localms(
            "xxxABCDxxx", "zzzABzzCDz", 1, -0.5, -3, -1, score_only=True
        )
        self.assertEqual(aligns1[0][2], aligns2)


class TestPairwiseOpenPenalty(unittest.TestCase):
    """Alignments with gap-open penalty."""

    def test_match_score_open_penalty1(self):
        """Test 1."""
        aligns = pairwise2.align.globalms("AA", "A", 2.0, -1, -0.1, 0)
        self.assertEqual(len(aligns), 2)
        aligns.sort()
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(
            alignment,
            """\
AA
 |
-A
  Score=1.9
""",
        )
        seq1, seq2, score, begin, end = aligns[1]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(
            alignment,
            """\
AA
| 
A-
  Score=1.9
""",  # noqa: W291
        )

    def test_match_score_open_penalty2(self):
        """Test 2."""
        aligns = pairwise2.align.globalms("GAA", "GA", 1.5, 0, -0.1, 0)
        self.assertEqual(len(aligns), 2)
        aligns.sort()
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(
            alignment,
            """\
GAA
| |
G-A
  Score=2.9
""",
        )
        seq1, seq2, score, begin, end = aligns[1]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(
            alignment,
            """\
GAA
|| 
GA-
  Score=2.9
""",  # noqa: W291
        )

    def test_match_score_open_penalty3(self):
        """Test 3."""
        aligns = pairwise2.align.globalxs("GAACT", "GAT", -0.1, 0)
        self.assertEqual(len(aligns), 1)
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(
            alignment,
            """\
GAACT
||  |
GA--T
  Score=2.9
""",
        )

    def test_match_score_open_penalty4(self):
        """Test 4."""
        aligns = pairwise2.align.globalms("GCT", "GATA", 1, -2, -0.1, 0)
        self.assertEqual(len(aligns), 1)
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(
            alignment,
            """\
GC-T-
|  | 
G-ATA
  Score=1.7
""",  # noqa: W291
        )


class TestPairwiseExtendPenalty(unittest.TestCase):
    """Alignments with gap-extend penalties."""

    def test_extend_penalty1(self):
        """Test 1."""
        aligns = pairwise2.align.globalxs("GACT", "GT", -0.5, -0.2)
        self.assertEqual(len(aligns), 1)
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(
            alignment,
            """\
GACT
|  |
G--T
  Score=1.3
""",
        )

    def test_extend_penalty2(self):
        """Test 2."""
        aligns = pairwise2.align.globalxs("GACT", "GT", -1.5, -0.2)
        self.assertEqual(len(aligns), 1)
        aligns.sort()
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(
            alignment,
            """\
GACT
|  |
G--T
  Score=0.3
""",
        )


class TestPairwisePenalizeExtendWhenOpening(unittest.TestCase):
    """Alignment with ``penalize_extend_when_opening``."""

    def test_penalize_extend_when_opening(self):
        """Add gap-extend penalty to gap-opening penalty."""
        aligns = pairwise2.align.globalxs(
            "GACT", "GT", -0.2, -1.5, penalize_extend_when_opening=1
        )
        self.assertEqual(len(aligns), 1)
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(
            alignment,
            """\
GACT
|  |
G--T
  Score=-1.2
""",
        )


class TestPairwisePenalizeEndgaps(unittest.TestCase):
    """Alignments with end-gaps penalized or not."""

    def test_penalize_end_gaps(self):
        """Turn off end-gap penalties."""
        aligns = pairwise2.align.globalxs("GACT", "GT", -0.8, -0.2, penalize_end_gaps=0)
        self.assertEqual(len(aligns), 3)
        aligns.sort()
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(
            alignment,
            """\
GACT
  .|
--GT
  Score=1
""",
        )
        seq1, seq2, score, begin, end = aligns[1]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(
            alignment,
            """\
GACT
|  |
G--T
  Score=1
""",
        )
        seq1, seq2, score, begin, end = aligns[2]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(
            alignment,
            """\
GACT
|.  
GT--
  Score=1
""",  # noqa: W291
        )

    def test_penalize_end_gaps2(self):
        """Do the same, but use the generic method (with the same result)."""
        aligns = pairwise2.align.globalxs(
            "GACT", "GT", -0.8, -0.2, penalize_end_gaps=0, force_generic=True
        )
        self.assertEqual(len(aligns), 3)
        aligns.sort()
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(
            alignment,
            """\
GACT
  .|
--GT
  Score=1
""",
        )
        seq1, seq2, score, begin, end = aligns[1]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(
            alignment,
            """\
GACT
|  |
G--T
  Score=1
""",
        )
        seq1, seq2, score, begin, end = aligns[2]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(
            alignment,
            """\
GACT
|.  
GT--
  Score=1
""",  # noqa: W291
        )

    def test_separate_penalize_end_gaps(self):
        """Test alignment where end-gaps are differently penalized."""
        align = pairwise2.align.globalms(
            "AT", "AGG", 1.0, -0.5, -1.75, -0.25, penalize_end_gaps=(True, False)
        )
        self.assertEqual(align[0], ("A--T", "AGG-", -1.0, 0, 4))


class TestPairwiseSeparateGapPenalties(unittest.TestCase):
    """Alignments with separate gap-open penalties for both sequences."""

    def test_separate_gap_penalties1(self):
        """Test 1."""
        aligns = pairwise2.align.localxd("GAT", "GTCT", -0.3, 0, -0.8, 0)
        self.assertEqual(len(aligns), 2)
        aligns.sort()
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(
            alignment,
            """\
G-AT
| .|
GTCT
  Score=1.7
""",
        )
        seq1, seq2, score, begin, end = aligns[1]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(
            alignment,
            """\
GA-T
|. |
GTCT
  Score=1.7
""",
        )

    def test_separate_gap_penalties2(self):
        """Test 2."""
        aligns = pairwise2.align.localxd("GAT", "GTCT", -0.5, 0, -0.2, 0)
        self.assertEqual(len(aligns), 1)
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(
            alignment,
            """\
1 GAT
  | |
1 G-T
  Score=1.8
""",
        )


class TestPairwiseSeparateGapPenaltiesWithExtension(unittest.TestCase):
    """Alignments with separate gap-extension penalties for both sequences."""

    def test_separate_gap_penalties_with_extension(self):
        """Test separate gap-extension penalties and list input."""
        aligns = pairwise2.align.localxd(
            list("GAAT"), list("GTCCT"), -0.1, 0, -0.1, -0.1, gap_char=["-"]
        )
        self.assertEqual(len(aligns), 3)
        aligns.sort()
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(
            alignment,
            """\
G - A A T 
|   . . | 
G T C C T 
  Score=1.9
""",  # noqa: W291
        )
        seq1, seq2, score, begin, end = aligns[1]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(
            alignment,
            """\
G A - A T 
| .   . | 
G T C C T 
  Score=1.9
""",  # noqa: W291
        )
        seq1, seq2, score, begin, end = aligns[2]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(
            alignment,
            """\
G A A - T 
| . .   | 
G T C C T 
  Score=1.9
""",  # noqa: W291
        )


class TestPairwiseMatchDictionary(unittest.TestCase):
    """Alignments with match dictionaries."""

    match_dict = {("A", "A"): 1.5, ("A", "T"): 0.5, ("T", "T"): 1.0}

    def test_match_dictionary1(self):
        """Test 1."""
        aligns = pairwise2.align.localds("ATAT", "ATT", self.match_dict, -0.5, 0)
        self.assertEqual(len(aligns), 2)
        aligns.sort()
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(
            alignment,
            """\
ATAT
|| |
AT-T
  Score=3
""",
        )
        seq1, seq2, score, begin, end = aligns[1]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(
            alignment,
            """\
1 ATA
  ||.
1 ATT
  Score=3
""",
        )

    def test_match_dictionary2(self):
        """Test 2."""
        aligns = pairwise2.align.localds("ATAT", "ATT", self.match_dict, -1, 0)
        self.assertEqual(len(aligns), 1)
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(
            alignment,
            """\
1 ATA
  ||.
1 ATT
  Score=3
""",
        )

    def test_match_dictionary3(self):
        """Test 3."""
        aligns = pairwise2.align.localds("ATT", "ATAT", self.match_dict, -1, 0)
        self.assertEqual(len(aligns), 1)
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(
            alignment,
            """\
1 ATT
  ||.
1 ATA
  Score=3
""",
        )


class TestPairwiseOneCharacter(unittest.TestCase):
    """Alignments where one sequence has length 1."""

    def test_align_one_char1(self):
        """Test sequence with only one match."""
        aligns = pairwise2.align.localxs("abcde", "c", -0.3, -0.1)
        self.assertEqual(len(aligns), 1)
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(
            alignment,
            """\
3 c
  |
1 c
  Score=1
""",
        )

    def test_align_one_char2(self):
        """Test sequences with two possible match positions."""
        aligns = pairwise2.align.localxs("abcce", "c", -0.3, -0.1)
        self.assertEqual(len(aligns), 2)
        aligns.sort()
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(
            alignment,
            """\
4 c
  |
1 c
  Score=1
""",
        )
        seq1, seq2, score, begin, end = aligns[1]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(
            alignment,
            """\
3 c
  |
1 c
  Score=1
""",
        )

    def test_align_one_char3(self):
        """Like test 1, but global alignment."""
        aligns = pairwise2.align.globalxs("abcde", "c", -0.3, -0.1)
        self.assertEqual(len(aligns), 1)
        seq1, seq2, score, begin, end = aligns[0]
        alignment = pairwise2.format_alignment(seq1, seq2, score, begin, end)
        self.assertEqual(
            alignment,
            """\
abcde
  |  
--c--
  Score=0.2
""",  # noqa: W291
        )


class TestPersiteGapPenalties(unittest.TestCase):
    """Check gap penalty callbacks use correct gap opening position.

    This tests that the gap penalty callbacks are really being used
    with the correct gap opening position.
    """

    def test_gap_here_only_1(self):
        """Open a gap in second sequence only."""
        seq1 = "AAAABBBAAAACCCCCCCCCCCCCCAAAABBBAAAA"
        seq2 = "AABBBAAAACCCCAAAABBBAA"

        def no_gaps(x, y):
            """Very expensive to open a gap in seq1."""
            return -2000 - y

        def specific_gaps(x, y):
            """Very expensive to open a gap in seq2.

            ...unless it is in one of the allowed positions:
            """
            breaks = [0, 11, len(seq2)]
            return (-2 - y) if x in breaks else (-2000 - y)

        alignments = pairwise2.align.globalmc(seq1, seq2, 1, -1, no_gaps, specific_gaps)
        self.assertEqual(len(alignments), 1)
        formatted = pairwise2.format_alignment(*alignments[0])
        self.assertEqual(
            formatted,
            """\
AAAABBBAAAACCCCCCCCCCCCCCAAAABBBAAAA
  |||||||||||          |||||||||||  
--AABBBAAAACC----------CCAAAABBBAA--
  Score=2
""",  # noqa: W291
        )

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
            return -2000 - y

        def specific_gaps(x, y):
            """Very expensive to open a gap in seq2.

            ...unless it is in one of the allowed positions:
            """
            breaks = [0, 3, len(seq2)]
            return (-2 - y) if x in breaks else (-2000 - y)

        alignments = pairwise2.align.globalmc(seq1, seq2, 1, -1, no_gaps, specific_gaps)
        self.assertEqual(len(alignments), 1)
        formatted = pairwise2.format_alignment(*alignments[0])
        self.assertEqual(
            formatted,
            """\
AAAABBBAAAACCCCCCCCCCCCCCAAAABBBAAAA
  |||          ......|||||||||||||  
--AAB----------BBAAAACCCCAAAABBBAA--
  Score=-10
""",  # noqa: W291
        )


class TestOtherFunctions(unittest.TestCase):
    """Test remaining non-tested private methods."""

    def test_clean_alignments(self):
        """``_clean_alignments`` removes redundant and wrong alignments."""
        alns = [
            ("ACCGT", "AC-G-", 3.0, 0, 4),
            ("ACCGT", "AC-G-", 3.0, 1, 1),
            ("ACCGT", "A-CG-", 3.0, 0, 4),
            ("ACCGT", "AC-G-", 3.0, 0, 4),
            ("ACCGT", "A-CG-", 3.0, 0, 4),
        ]
        expected = [("ACCGT", "AC-G-", 3.0, 0, 4), ("ACCGT", "A-CG-", 3.0, 0, 4)]
        result = pairwise2._clean_alignments(alns)
        self.assertEqual(expected, result)

    def test_alignments_can_be_pickled(self):
        alns = [("ACCGT", "AC-G-", 3.0, 0, 4)]
        expected = [("ACCGT", "AC-G-", 3.0, 0, 4)]
        result = pickle.loads(pickle.dumps(pairwise2._clean_alignments(alns)))
        self.assertEqual(expected, result)

    def test_print_matrix(self):
        """``print_matrix`` prints nested lists as nice matrices."""
        import sys
        from io import StringIO

        out = StringIO()
        sys.stdout = out
        pairwise2.print_matrix(
            [
                [0.0, -1.0, -1.5, -2.0],
                [-1.0, 4.0, 3.0, 2.5],
                [-1.5, 3.0, 8.0, 7.0],
                [-2.0, 2.5, 7.0, 6.0],
                [-2.5, 2.0, 6.5, 11.0],
                [-3.0, 1.5, 6.0, 10.0],
            ]
        )
        self.assertEqual(
            out.getvalue(),
            " 0.0  -1.0  -1.5  -2.0 \n"
            "-1.0   4.0   3.0   2.5 \n"
            "-1.5   3.0   8.0   7.0 \n"
            "-2.0   2.5   7.0   6.0 \n"
            "-2.5   2.0   6.5  11.0 \n"
            "-3.0   1.5   6.0  10.0 \n",
        )
        sys.stdout = sys.__stdout__

    def test_recover_alignments(self):
        """One possible start position in local alignment is not a match."""
        self.assertEqual(len(pairwise2.align.localxx("AC", "GA")), 1)


if __name__ == "__main__":
    if pairwise2.rint != pairwise2._python_rint:
        # This uses the default C extensions, if import didn't fail.
        runner = unittest.TextTestRunner(verbosity=2)
        unittest.main(testRunner=runner, exit=False)
    else:
        print(
            "Import of C functions failed. Only testing pure Python "
            "fallback functions."
        )

    # Now, we switch explicitly to the fallback Python functions:
    pairwise2._make_score_matrix_fast = pairwise2._python_make_score_matrix_fast
    pairwise2.rint = pairwise2._python_rint

    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
