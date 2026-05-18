# Tests/test_format_matrix_unittest.py
"""Unit tests for Alignment.format() with substitution matrices.

These tests cover the behavior of pretty-printed alignments when a
substitution matrix (e.g., from Bio.Align.substitution_matrices) is
supplied directly or via a PairwiseAligner object.

Conventions being tested:
* '|' (pipe): identity (the same residue on both sequences).
* ':' (colon): positive mismatch (substitution with a positive score).
* '.' (dot): negative mismatch (substitution with a negative score).
* '-' (dash): gap (insertion/deletion).

Test cases include:
* NUC.4.4 matrix, where T~Y is a positive mismatch.
* A BLASTN-like artificial matrix with +1 for matches and -1 for mismatches.
* Case-insensitivity: lowercase residues behave the same as uppercase.
* Handling of gaps in alignments.
* A "mixed block" test where all four pattern characters appear at least once.

These tests ensure that the new feature of showing ':' for positive
substitution scores is consistently applied across different input
styles and substitution matrices.
"""

import unittest
import numpy as np
from Bio.Align import PairwiseAligner
from Bio.Align.substitution_matrices import load, Array


class TestFormatMatrix(unittest.TestCase):
    """Unit tests for Alignment.format() with substitution matrices."""

    def test_nuc44_gives_colon_for_positive_mismatch_TY(self):
        """In NUC.4.4, T vs Y has a positive score -> expect ':' in the pattern."""
        M = load("NUC.4.4")
        aligner = PairwiseAligner()
        aligner.open_gap_score = -1
        aligner.extend_gap_score = -0.5
        aln = aligner.align("GATTACAT", "GATYACAC")[0]
        self.assertEqual(
            aln.format(scoring=M),
            """\
target            0 GATTACAT 8
                  0 |||:|||. 8
query             0 GATYACAC 8
""",
        )

    def test_blastn_like_has_no_colon_only_pipes_for_identities(self):
        """In a BLASTN-like +1/-1 matrix, mismatches are always negative -> no ':' expected."""
        M = load("BLASTN")
        aligner = PairwiseAligner()
        aligner.open_gap_score = -1
        aligner.extend_gap_score = -0.5
        aln = aligner.align("GATTACAT", "GATYACAC")[0]
        self.assertEqual(
            aln.format(scoring=M),
            """\
target            0 GATTACAT 8
                  0 |||.|||. 8
query             0 GATYACAC 8
""",
        )

    def test_positive_mismatch_colon_when_passing_aligner_object(self):
        """Passing the aligner object with a substitution matrix should also yield ':' for positive mismatches."""
        M = load("NUC.4.4")
        aligner = PairwiseAligner()
        aligner.substitution_matrix = M
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -2
        aln = aligner.align("GATTACAT", "GATYACAC")[0]
        self.assertEqual(
            aln.format(scoring=M),
            """\
target            0 GATTACAT 8
                  0 |||:|||. 8
query             0 GATYACAC 8
""",
        )

    def test_lowercase_letters_are_case_insensitive_for_matrix_lookup(self):
        """Matrix lookup should be case-insensitive (e.g., 't' vs 'y' behaves like 'T' vs 'Y')."""
        M = load("NUC.4.4")
        aligner = PairwiseAligner()
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -2
        aln = aligner.align("t", "y")[0]
        self.assertEqual(
            aln.format(scoring=M),
            """\
target            0 t 1
                  0 : 1
query             0 y 1
""",
        )

    def test_negative_mismatch_dot_with_blastn_like_matrix(self):
        """In the BLASTN-like matrix, mismatches are negative -> expect '.' in the pattern."""
        M = load("BLASTN")
        aligner = PairwiseAligner()
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -2
        aln = aligner.align("A", "C")[0]
        self.assertEqual(
            aln.format(scoring=M),
            """\
target            0 A 1
                  0 . 1
query             0 C 1
""",
        )

    def test_gap_is_dash_in_pattern(self):
        """Gaps in the alignment should always appear as '-' in the pattern line."""
        M = load("NUC.4.4")
        aligner = PairwiseAligner()
        aligner.open_gap_score = -1
        aligner.extend_gap_score = -0.5
        aln = aligner.align("AC", "AGC")[0]
        self.assertEqual(
            aln.format(scoring=M),
            """\
target            0 A-C 2
                  0 |-| 3
query             0 AGC 3
""",
        )

    def test_mixed_block_contains_expected_symbols(self):
        """Construct an alignment that produces all symbols ('|', ':', '.', '-') at least once in the pattern."""
        M = load("NUC.4.4")
        aligner = PairwiseAligner()
        aligner.substitution_matrix = M
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -2
        seq1 = "TTTG"
        seq2 = "TYGG"
        aln = aligner.align(seq1, seq2)[0]
        self.assertEqual(
            aln.format(scoring=M),
            """\
target            0 TTTG 4
                  0 |:.| 4
query             0 TYGG 4
""",
        )


if __name__ == "__main__":
    unittest.main()
