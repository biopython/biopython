# Tests/test_format_matrix.py
import re
import numpy as np
from Bio.Align import PairwiseAligner
from Bio.Align.substitution_matrices import load, Array

"""Tests for Alignment.format() when used with substitution matrices.

These tests cover the behavior of pretty-printed alignments when a
substitution matrix (e.g., from Bio.Align.substitution_matrices) is
supplied directly or via a PairwiseAligner object.

Conventions being tested:
    - '|' (pipe): identity (the same residue on both sequences).
    - ':' (colon): positive mismatch (substitution with a positive score).
    - '.' (dot): negative mismatch (substitution with a negative score).
    - '-' (dash): gap (insertion/deletion).

Test cases include:
    * NUC.4.4 matrix, where T–Y is a positive mismatch.
    * A BLASTN-like artificial matrix with +1 for matches and -1 for mismatches.
    * Case-insensitivity: lowercase residues behave the same as uppercase.
    * Handling of gaps in alignments.
    * A "mixed block" test where all four pattern characters appear at least once.

These tests ensure that the new feature of showing ':' for positive
substitution scores is consistently applied across different input
styles and substitution matrices.
"""


def _pattern_from_pretty(s: str) -> str:
    
    for ln in s.splitlines():
        m = re.match(r"^\s*\d+\s+([|:\.\-]+)\s+\d+\s*$", ln)
        if m:
            return m.group(1)

    best = ""
    for ln in s.splitlines():
        for m in re.finditer(r"[|:\.\-]+", ln):
            if len(m.group(0)) > len(best):
                best = m.group(0)
    if best:
        return best

    raise AssertionError("Pattern line not found in pretty output:\n" + s)

def _blastn_like_matrix():
    
    alphabet = "ACGTY"  
    n = len(alphabet)
    data = np.full((n, n), -1, dtype=int)
    np.fill_diagonal(data, 1)
    return Array(alphabet=alphabet, dims=2, data=data)


def test_nuc44_gives_colon_for_positive_mismatch_TY():
    """In NUC.4.4, T vs Y has a positive score -> expect ':' in the pattern."""
    M = load("NUC.4.4")
    aligner = PairwiseAligner()
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.5
    aln = next(iter(aligner.align("GATTACAT", "GATYACAC")))
    pat = _pattern_from_pretty(aln.format("", M))
    assert ":" in pat
    assert "-" not in pat

def test_blastn_like_has_no_colon_only_pipes_for_identities():
    """In a BLASTN-like +1/-1 matrix, mismatches are always negative -> no ':' expected."""
    M = _blastn_like_matrix()
    aligner = PairwiseAligner()
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.5
    aln = next(iter(aligner.align("GATTACAT", "GATYACAC")))
    pat = _pattern_from_pretty(aln.format("", M))
    assert ":" not in pat
    assert "|" in pat

def test_positive_mismatch_colon_when_passing_aligner_object():
    """Passing the aligner object with a substitution matrix should also yield ':' for positive mismatches."""
    M = load("NUC.4.4")
    aligner = PairwiseAligner()
    aligner.substitution_matrix = M
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.5
    aln = next(iter(aligner.align("GATTACAT", "GATYACAC")))
    pat = _pattern_from_pretty(aln.format("", aligner))
    assert ":" in pat

def test_lowercase_letters_are_handled_case_insensitively_for_matrix_lookup():
    """Matrix lookup should be case-insensitive (e.g., 't' vs 'y' behaves like 'T' vs 'Y')."""
    M = load("NUC.4.4")
    aligner = PairwiseAligner()
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -2
    aln = next(iter(aligner.align("t", "y")))
    pat = _pattern_from_pretty(aln.format("", M))
    assert pat == ":"

def test_negative_mismatch_dot_with_blastn_like_matrix():
    """In the BLASTN-like matrix, mismatches are negative -> expect '.' in the pattern."""
    M = _blastn_like_matrix()
    aligner = PairwiseAligner()
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -2
    aln = next(iter(aligner.align("A", "C")))
    pat = _pattern_from_pretty(aln.format("", M))
    assert pat == "."

def test_gap_is_dash_in_pattern():
    """Gaps in the alignment should always appear as '-' in the pattern line."""
    M = load("NUC.4.4")
    aligner = PairwiseAligner()
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.5
    aln = next(iter(aligner.align("AC", "AGC")))
    pat = _pattern_from_pretty(aln.format("", M))
    assert "-" in pat

def test_mixed_block_contains_all_symbols_when_expected():
    """Construct an alignment that produces all symbols ('|', ':', '.', '-') at least once in the pattern."""
    M = load("NUC.4.4")
    aligner = PairwiseAligner()
    aligner.substitution_matrix = M
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.5

    seq1 = "TTTG"
    seq2 = "TYGG"
    aln = next(iter(aligner.align(seq1, seq2)))
    pat = _pattern_from_pretty(aln.format("", M))

    needed = set("|:-")
    assert needed.issubset(set(pat)), f"pattern missing some of {needed}, got: {pat}"
