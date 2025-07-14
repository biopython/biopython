# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


"""Tests for the PairwiseAligner in Bio.Align."""

import array
import os
import sys
import unittest

try:
    import numpy as np
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align."
    ) from None

from Bio import BiopythonDeprecationWarning
from Bio import BiopythonWarning
from Bio import Align
from Bio.Align.substitution_matrices import Array
from Bio import SeqIO
from Bio.Seq import reverse_complement
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class TestAlignerProperties(unittest.TestCase):
    def test_aligner_property_epsilon(self):
        aligner = Align.PairwiseAligner()
        self.assertAlmostEqual(aligner.epsilon, 1.0e-6)
        aligner.epsilon = 1.0e-4
        self.assertAlmostEqual(aligner.epsilon, 1.0e-4)
        aligner.epsilon = 1.0e-8
        self.assertAlmostEqual(aligner.epsilon, 1.0e-8, places=8)
        with self.assertRaises(TypeError):
            aligner.epsilon = "not a number"
        with self.assertRaises(TypeError):
            aligner.epsilon = None

    def test_aligner_property_mode(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        self.assertEqual(aligner.mode, "global")
        aligner.mode = "local"
        self.assertEqual(aligner.mode, "local")
        with self.assertRaises(ValueError):
            aligner.mode = "wrong"

    def test_aligner_property_match_mismatch(self):
        aligner = Align.PairwiseAligner()
        aligner.match_score = 3.0
        self.assertAlmostEqual(aligner.match_score, 3.0)
        aligner.mismatch_score = -2.0
        self.assertAlmostEqual(aligner.mismatch_score, -2.0)
        with self.assertRaises(ValueError):
            aligner.match_score = "not a number"
        with self.assertRaises(ValueError):
            aligner.mismatch_score = "not a number"
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  wildcard: None
  match_score: 3.000000
  mismatch_score: -2.000000
  open_internal_insertion_score: 0.000000
  extend_internal_insertion_score: 0.000000
  open_left_insertion_score: 0.000000
  extend_left_insertion_score: 0.000000
  open_right_insertion_score: 0.000000
  extend_right_insertion_score: 0.000000
  open_internal_deletion_score: 0.000000
  extend_internal_deletion_score: 0.000000
  open_left_deletion_score: 0.000000
  extend_left_deletion_score: 0.000000
  open_right_deletion_score: 0.000000
  extend_right_deletion_score: 0.000000
  mode: global
""",
        )

    def test_aligner_property_gapscores(self):
        aligner = Align.PairwiseAligner()
        open_score, extend_score = (-5, -1)
        aligner.open_insertion_score = open_score
        aligner.extend_insertion_score = extend_score
        self.assertAlmostEqual(aligner.open_insertion_score, open_score)
        self.assertAlmostEqual(aligner.extend_insertion_score, extend_score)
        open_score, extend_score = (-6, -7)
        aligner.open_deletion_score = open_score
        aligner.extend_deletion_score = extend_score
        self.assertAlmostEqual(aligner.open_deletion_score, open_score)
        self.assertAlmostEqual(aligner.extend_deletion_score, extend_score)
        open_score, extend_score = (-3, -9)
        aligner.open_end_insertion_score = open_score
        aligner.extend_end_insertion_score = extend_score
        self.assertAlmostEqual(aligner.open_end_insertion_score, open_score)
        self.assertAlmostEqual(aligner.extend_end_insertion_score, extend_score)
        open_score, extend_score = (-1, -2)
        aligner.open_end_deletion_score = open_score
        aligner.extend_end_deletion_score = extend_score
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  wildcard: None
  match_score: 1.000000
  mismatch_score: 0.000000
  open_internal_insertion_score: -5.000000
  extend_internal_insertion_score: -1.000000
  open_left_insertion_score: -3.000000
  extend_left_insertion_score: -9.000000
  open_right_insertion_score: -3.000000
  extend_right_insertion_score: -9.000000
  open_internal_deletion_score: -6.000000
  extend_internal_deletion_score: -7.000000
  open_left_deletion_score: -1.000000
  extend_left_deletion_score: -2.000000
  open_right_deletion_score: -1.000000
  extend_right_deletion_score: -2.000000
  mode: global
""",
        )
        self.assertAlmostEqual(aligner.open_end_deletion_score, open_score)
        self.assertAlmostEqual(aligner.extend_end_deletion_score, extend_score)
        score = -3
        aligner.insertion_score = score
        self.assertAlmostEqual(aligner.insertion_score, score)
        self.assertAlmostEqual(aligner.open_insertion_score, score)
        self.assertAlmostEqual(aligner.extend_insertion_score, score)
        score = -2
        aligner.deletion_score = score
        self.assertAlmostEqual(aligner.deletion_score, score)
        self.assertAlmostEqual(aligner.open_deletion_score, score)
        self.assertAlmostEqual(aligner.extend_deletion_score, score)
        score = -4
        aligner.end_insertion_score = score
        self.assertAlmostEqual(aligner.end_insertion_score, score)
        self.assertAlmostEqual(aligner.open_end_insertion_score, score)
        self.assertAlmostEqual(aligner.extend_end_insertion_score, score)
        self.assertAlmostEqual(aligner.left_insertion_score, score)
        self.assertAlmostEqual(aligner.open_left_insertion_score, score)
        self.assertAlmostEqual(aligner.extend_left_insertion_score, score)
        self.assertAlmostEqual(aligner.right_insertion_score, score)
        self.assertAlmostEqual(aligner.open_right_insertion_score, score)
        self.assertAlmostEqual(aligner.extend_right_insertion_score, score)
        score = -5
        aligner.end_deletion_score = score
        self.assertAlmostEqual(aligner.end_deletion_score, score)
        self.assertAlmostEqual(aligner.open_end_deletion_score, score)
        self.assertAlmostEqual(aligner.extend_end_deletion_score, score)
        self.assertAlmostEqual(aligner.left_deletion_score, score)
        self.assertAlmostEqual(aligner.open_left_deletion_score, score)
        self.assertAlmostEqual(aligner.extend_left_deletion_score, score)
        self.assertAlmostEqual(aligner.right_deletion_score, score)
        self.assertAlmostEqual(aligner.open_right_deletion_score, score)
        self.assertAlmostEqual(aligner.extend_right_deletion_score, score)
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  wildcard: None
  match_score: 1.000000
  mismatch_score: 0.000000
  open_internal_insertion_score: -3.000000
  extend_internal_insertion_score: -3.000000
  open_left_insertion_score: -4.000000
  extend_left_insertion_score: -4.000000
  open_right_insertion_score: -4.000000
  extend_right_insertion_score: -4.000000
  open_internal_deletion_score: -2.000000
  extend_internal_deletion_score: -2.000000
  open_left_deletion_score: -5.000000
  extend_left_deletion_score: -5.000000
  open_right_deletion_score: -5.000000
  extend_right_deletion_score: -5.000000
  mode: global
""",
        )
        with self.assertRaises(ValueError):
            aligner.insertion_score = "wrong"
        with self.assertRaises(ValueError):
            aligner.deletion_score = "wrong"
        with self.assertRaises(TypeError):
            aligner.end_insertion_score = "wrong"
        with self.assertRaises(TypeError):
            aligner.end_deletion_score = "wrong"

    def test_aligner_property_gapscores_deprecated(self):
        aligner = Align.PairwiseAligner()
        value = 1
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.target_left_open_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.target_left_open_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 2
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.query_left_open_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.query_left_open_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 3
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.left_open_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.left_open_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 4
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.target_internal_open_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.target_internal_open_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 5
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.query_internal_open_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.query_internal_open_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 6
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.internal_open_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.internal_open_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 7
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.target_right_open_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.target_right_open_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 8
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.query_right_open_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.query_right_open_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 9
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.right_open_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.right_open_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 10
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.target_end_open_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.target_end_open_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 11
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.query_end_open_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.query_end_open_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 12
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.end_open_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.end_open_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 13
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.target_open_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.target_open_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 14
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.query_open_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.query_open_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 15
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.target_left_extend_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.target_left_extend_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 16
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.query_left_extend_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.query_left_extend_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 17
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.left_extend_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.left_extend_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 18
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.target_internal_extend_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.target_internal_extend_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 19
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.query_internal_extend_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.query_internal_extend_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 20
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.internal_extend_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.internal_extend_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 21
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.target_right_extend_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.target_right_extend_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 22
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.query_right_extend_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.query_right_extend_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 23
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.right_extend_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.right_extend_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 24
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.target_end_extend_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.target_end_extend_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 25
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.query_end_extend_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.query_end_extend_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 26
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.end_extend_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.end_extend_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 27
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.target_extend_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.target_extend_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 28
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.query_extend_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.query_extend_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 29
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.target_left_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.target_left_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 30
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.query_left_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.query_left_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 31
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.target_internal_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.target_internal_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 32
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.query_internal_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.query_internal_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 33
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.target_right_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.target_right_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 34
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.query_right_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.query_right_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 35
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.target_end_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.target_end_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 36
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.query_end_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.query_end_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 37
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.target_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.target_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 38
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.query_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.query_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 39
        aligner.open_left_insertion_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.target_left_open_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 40
        aligner.open_left_deletion_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.query_left_open_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 41
        aligner.open_left_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.left_open_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 42
        aligner.open_internal_insertion_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.target_internal_open_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 43
        aligner.open_internal_deletion_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.query_internal_open_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 44
        aligner.open_internal_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.internal_open_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 45
        aligner.open_right_insertion_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.target_right_open_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 46
        aligner.open_right_deletion_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.query_right_open_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 47
        aligner.open_right_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.right_open_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 48
        aligner.open_end_insertion_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.target_end_open_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 49
        aligner.open_end_deletion_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.query_end_open_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 50
        aligner.open_end_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.end_open_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 51
        aligner.open_insertion_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.target_open_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 52
        aligner.open_deletion_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.query_open_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 53
        aligner.extend_left_insertion_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.target_left_extend_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 54
        aligner.extend_left_deletion_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.query_left_extend_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 55
        aligner.extend_left_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.left_extend_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 56
        aligner.extend_internal_insertion_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.target_internal_extend_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 57
        aligner.extend_internal_deletion_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.query_internal_extend_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 58
        aligner.extend_internal_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.internal_extend_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 59
        aligner.extend_right_insertion_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.target_right_extend_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 60
        aligner.extend_right_deletion_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.query_right_extend_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 61
        aligner.extend_right_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.right_extend_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 62
        aligner.extend_end_insertion_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.target_end_extend_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 63
        aligner.extend_end_deletion_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.query_end_extend_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 64
        aligner.extend_end_gap_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.end_extend_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 65
        aligner.extend_insertion_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.target_extend_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 66
        aligner.extend_deletion_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.query_extend_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 67
        aligner.left_insertion_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.target_left_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 68
        aligner.left_deletion_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.query_left_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 69
        aligner.internal_insertion_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.target_internal_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 70
        aligner.internal_deletion_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.query_internal_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 71
        aligner.right_insertion_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.target_right_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 72
        aligner.right_deletion_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.query_right_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 73
        aligner.end_insertion_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.target_end_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 74
        aligner.end_deletion_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.query_end_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 75
        aligner.insertion_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.target_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 76
        aligner.deletion_score = value
        with self.assertWarns(BiopythonDeprecationWarning):
            stored_value = aligner.query_gap_score
        self.assertAlmostEqual(stored_value, value)
        value = 77
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.target_left_open_gap_score = value
        self.assertAlmostEqual(aligner.open_left_insertion_score, value)
        value = 78
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.query_left_open_gap_score = value
        self.assertAlmostEqual(aligner.open_left_deletion_score, value)
        value = 79
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.left_open_gap_score = value
        self.assertAlmostEqual(aligner.open_left_gap_score, value)
        value = 80
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.target_internal_open_gap_score = value
        self.assertAlmostEqual(aligner.open_internal_insertion_score, value)
        value = 81
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.query_internal_open_gap_score = value
        self.assertAlmostEqual(aligner.open_internal_deletion_score, value)
        value = 82
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.internal_open_gap_score = value
        self.assertAlmostEqual(aligner.open_internal_gap_score, value)
        value = 83
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.target_right_open_gap_score = value
        self.assertAlmostEqual(aligner.open_right_insertion_score, value)
        value = 84
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.query_right_open_gap_score = value
        self.assertAlmostEqual(aligner.open_right_deletion_score, value)
        value = 85
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.right_open_gap_score = value
        self.assertAlmostEqual(aligner.open_right_gap_score, value)
        value = 86
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.target_end_open_gap_score = value
        self.assertAlmostEqual(aligner.open_end_insertion_score, value)
        value = 87
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.query_end_open_gap_score = value
        self.assertAlmostEqual(aligner.open_end_deletion_score, value)
        value = 88
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.end_open_gap_score = value
        self.assertAlmostEqual(aligner.open_end_gap_score, value)
        value = 89
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.target_open_gap_score = value
        self.assertAlmostEqual(aligner.open_insertion_score, value)
        value = 90
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.query_open_gap_score = value
        self.assertAlmostEqual(aligner.open_deletion_score, value)
        value = 91
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.target_left_extend_gap_score = value
        self.assertAlmostEqual(aligner.extend_left_insertion_score, value)
        value = 92
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.query_left_extend_gap_score = value
        self.assertAlmostEqual(aligner.extend_left_deletion_score, value)
        value = 93
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.left_extend_gap_score = value
        self.assertAlmostEqual(aligner.extend_left_gap_score, value)
        value = 94
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.target_internal_extend_gap_score = value
        self.assertAlmostEqual(aligner.extend_internal_insertion_score, value)
        value = 95
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.query_internal_extend_gap_score = value
        self.assertAlmostEqual(aligner.extend_internal_deletion_score, value)
        value = 96
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.internal_extend_gap_score = value
        self.assertAlmostEqual(aligner.extend_internal_gap_score, value)
        value = 97
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.target_right_extend_gap_score = value
        self.assertAlmostEqual(aligner.extend_right_insertion_score, value)
        value = 98
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.query_right_extend_gap_score = value
        self.assertAlmostEqual(aligner.extend_right_deletion_score, value)
        value = 99
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.right_extend_gap_score = value
        self.assertAlmostEqual(aligner.extend_right_gap_score, value)
        value = 100
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.target_end_extend_gap_score = value
        self.assertAlmostEqual(aligner.extend_end_insertion_score, value)
        value = 101
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.query_end_extend_gap_score = value
        self.assertAlmostEqual(aligner.extend_end_deletion_score, value)
        value = 102
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.end_extend_gap_score = value
        self.assertAlmostEqual(aligner.extend_end_gap_score, value)
        value = 103
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.target_extend_gap_score = value
        self.assertAlmostEqual(aligner.extend_insertion_score, value)
        value = 104
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.query_extend_gap_score = value
        self.assertAlmostEqual(aligner.extend_deletion_score, value)
        value = 105
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.target_left_gap_score = value
        self.assertAlmostEqual(aligner.left_insertion_score, value)
        value = 106
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.query_left_gap_score = value
        self.assertAlmostEqual(aligner.left_deletion_score, value)
        value = 107
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.target_internal_gap_score = value
        self.assertAlmostEqual(aligner.internal_insertion_score, value)
        value = 108
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.query_internal_gap_score = value
        self.assertAlmostEqual(aligner.internal_deletion_score, value)
        value = 109
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.target_right_gap_score = value
        self.assertAlmostEqual(aligner.right_insertion_score, value)
        value = 110
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.query_right_gap_score = value
        self.assertAlmostEqual(aligner.right_deletion_score, value)
        value = 111
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.target_end_gap_score = value
        self.assertAlmostEqual(aligner.end_insertion_score, value)
        value = 112
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.query_end_gap_score = value
        self.assertAlmostEqual(aligner.end_deletion_score, value)
        value = 113
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.target_gap_score = value
        self.assertAlmostEqual(aligner.insertion_score, value)
        value = 114
        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.query_gap_score = value
        self.assertAlmostEqual(aligner.deletion_score, value)

        def gap_function1(x, y):
            return x + y

        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.query_gap_score = gap_function1
        gap_function = aligner.deletion_score
        self.assertEqual(gap_function, gap_function1)

        def gap_function2(x, y):
            return x * y

        aligner.deletion_score = gap_function2
        self.assertEqual(aligner.deletion_score, gap_function2)

        def gap_function3(x, y):
            return x / y

        aligner.deletion_score = gap_function3
        with self.assertWarns(BiopythonDeprecationWarning):
            gap_function = aligner.query_gap_score
        self.assertEqual(gap_function, gap_function3)

        def gap_function4(x, y):
            return x / y - 9

        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.query_gap_score = gap_function4
        with self.assertWarns(BiopythonDeprecationWarning):
            gap_function = aligner.query_gap_score
        self.assertEqual(gap_function, gap_function4)

        def gap_function5(x, y):
            return x + 2 * y

        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.target_gap_score = gap_function5
        gap_function = aligner.insertion_score
        self.assertEqual(gap_function, gap_function5)

        def gap_function6(x, y):
            return x * y + 2

        aligner.insertion_score = gap_function6
        self.assertEqual(aligner.insertion_score, gap_function6)

        def gap_function7(x, y):
            return x / y - 2

        aligner.insertion_score = gap_function7
        with self.assertWarns(BiopythonDeprecationWarning):
            gap_function = aligner.target_gap_score
        self.assertEqual(gap_function, gap_function7)

        def gap_function8(x, y):
            return x / y * 2

        with self.assertWarns(BiopythonDeprecationWarning):
            aligner.target_gap_score = gap_function8
        with self.assertWarns(BiopythonDeprecationWarning):
            gap_function = aligner.target_gap_score
        self.assertEqual(gap_function, gap_function8)

        def gap_function9(x, y):
            return x + 9 * y

        aligner.gap_score = gap_function9
        gap_function = aligner.gap_score
        self.assertEqual(gap_function, gap_function9)

    def test_aligner_nonexisting_property(self):
        aligner = Align.PairwiseAligner()
        with self.assertRaises(AttributeError) as cm:
            aligner.no_such_property
        self.assertEqual(
            str(cm.exception),
            "'PairwiseAligner' object has no attribute 'no_such_property'",
        )
        with self.assertRaises(AttributeError) as cm:
            aligner.no_such_property = 1
        self.assertEqual(
            str(cm.exception),
            "'PairwiseAligner' object has no attribute 'no_such_property'",
        )


class TestPairwiseGlobal(unittest.TestCase):
    def test_needlemanwunsch_simple1(self):
        seq1 = "GAACT"
        seq2 = "GAT"
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  wildcard: None
  match_score: 1.000000
  mismatch_score: 0.000000
  open_internal_insertion_score: 0.000000
  extend_internal_insertion_score: 0.000000
  open_left_insertion_score: 0.000000
  extend_left_insertion_score: 0.000000
  open_right_insertion_score: 0.000000
  extend_right_insertion_score: 0.000000
  open_internal_deletion_score: 0.000000
  extend_internal_deletion_score: 0.000000
  open_left_deletion_score: 0.000000
  extend_left_deletion_score: 0.000000
  open_right_deletion_score: 0.000000
  extend_right_deletion_score: 0.000000
  mode: global
""",
        )
        self.assertEqual(aligner.algorithm, "Needleman-Wunsch")
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 3.0)
        score = aligner.score(seq1, reverse_complement(seq2), "-")
        self.assertAlmostEqual(score, 3.0)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (2 alignments; score=3) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 GAACT 5
                  0 ||--| 5
query             0 GA--T 3
""",
        )
        self.assertEqual(alignment.shape, (2, 5))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 2], [4, 5]], [[0, 2], [2, 3]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (3 aligned letters; 3 identities; 0 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 3:
        identities = 3,
        mismatches = 0.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 GAACT 5
                  0 |-|-| 5
query             0 G-A-T 3
""",
        )
        self.assertEqual(alignment.shape, (2, 5))
        self.assertTrue(
            np.array_equal(
                alignment.aligned,
                np.array([[[0, 1], [2, 3], [4, 5]], [[0, 1], [1, 2], [2, 3]]]),
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (3 aligned letters; 3 identities; 0 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 3:
        identities = 3,
        mismatches = 0.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 2,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (2 alignments; score=3) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 GAACT 5
                  0 ||--| 5
query             3 GA--T 0
""",
        )
        self.assertEqual(alignment.shape, (2, 5))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 2], [4, 5]], [[3, 1], [1, 0]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (3 aligned letters; 3 identities; 0 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 3:
        identities = 3,
        mismatches = 0.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 GAACT 5
                  0 |-|-| 5
query             3 G-A-T 0
""",
        )
        self.assertEqual(alignment.shape, (2, 5))
        self.assertTrue(
            np.array_equal(
                alignment.aligned,
                np.array([[[0, 1], [2, 3], [4, 5]], [[3, 2], [2, 1], [1, 0]]]),
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (3 aligned letters; 3 identities; 0 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 3:
        identities = 3,
        mismatches = 0.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 2,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)

    def test_align_affine1_score(self):
        seq1 = "CC"
        seq2 = "ACCT"
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score = 0
        aligner.mismatch_score = -1
        aligner.open_gap_score = -5
        aligner.extend_gap_score = -1
        self.assertEqual(aligner.algorithm, "Gotoh global alignment algorithm")
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  wildcard: None
  match_score: 0.000000
  mismatch_score: -1.000000
  open_internal_insertion_score: -5.000000
  extend_internal_insertion_score: -1.000000
  open_left_insertion_score: -5.000000
  extend_left_insertion_score: -1.000000
  open_right_insertion_score: -5.000000
  extend_right_insertion_score: -1.000000
  open_internal_deletion_score: -5.000000
  extend_internal_deletion_score: -1.000000
  open_left_deletion_score: -5.000000
  extend_left_deletion_score: -1.000000
  open_right_deletion_score: -5.000000
  extend_right_deletion_score: -1.000000
  mode: global
""",
        )
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, -7.0)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, -7.0)

    def test_fogsaa_simple1(self):
        seq1 = "GAACT"
        seq2 = "GAT"
        aligner = Align.PairwiseAligner(mode="fogsaa")
        self.assertEqual(
            aligner.algorithm, "Fast Optimal Global Sequence Alignment Algorithm"
        )
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 3.0)
        score = aligner.score(seq1, reverse_complement(seq2), "-")
        self.assertAlmostEqual(score, 3.0)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=3) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 GAACT 5
                  0 ||--| 5
query             0 GA--T 3
""",
        )
        self.assertEqual(alignment.shape, (2, 5))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 2], [4, 5]], [[0, 2], [2, 3]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (3 aligned letters; 3 identities; 0 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 3:
        identities = 3,
        mismatches = 0.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)

        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=3) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 GAACT 5
                  0 ||--| 5
query             3 GA--T 0
""",
        )
        self.assertEqual(alignment.shape, (2, 5))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 2], [4, 5]], [[3, 1], [1, 0]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (3 aligned letters; 3 identities; 0 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 3:
        identities = 3,
        mismatches = 0.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)

    def test_fogsaa_affine1(self):
        seq1 = "CC"
        seq2 = "ACCT"
        aligner = Align.PairwiseAligner(mode="fogsaa")
        aligner.match_score = 0
        aligner.mismatch_score = -1
        aligner.open_gap_score = -5
        aligner.extend_gap_score = -1
        self.assertEqual(
            aligner.algorithm, "Fast Optimal Global Sequence Alignment Algorithm"
        )
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, -7.0)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, -7.0)

    def test_fogsaa_confirms_needleman_wunsch(self):
        seq1 = "CCCCC"
        seq2 = "ACCCCCT"

        aligner_fogsaa = Align.PairwiseAligner(mode="fogsaa")
        aligner_fogsaa.match_score = 1.1
        aligner_fogsaa.mismatch_score = -1.83
        aligner_fogsaa.gap_score = -2
        self.assertEqual(
            aligner_fogsaa.algorithm, "Fast Optimal Global Sequence Alignment Algorithm"
        )

        aligner_nw = Align.PairwiseAligner(mode="global")
        aligner_nw.match_score = 1.1
        aligner_nw.mismatch_score = -1.83
        aligner_nw.gap_score = -2
        self.assertEqual(aligner_nw.algorithm, "Needleman-Wunsch")

        score_fogsaa = aligner_fogsaa.score(seq1, seq2)
        score_nw = aligner_nw.score(seq1, seq2)
        self.assertAlmostEqual(score_fogsaa, score_nw)

    def test_fogsaa_matrix_scoring(self):
        seq1 = "AAAAAAAAAAA"
        seq2 = "AAAAAAATAAA"
        aligner = Align.PairwiseAligner(mode="fogsaa", scoring="blastn")
        self.assertEqual(
            aligner.algorithm, "Fast Optimal Global Sequence Alignment Algorithm"
        )
        with self.assertWarns(BiopythonWarning):
            score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 17.0)


class TestPairwiseLocal(unittest.TestCase):
    def test_smithwaterman(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = "local"
        aligner.gap_score = -0.1
        self.assertEqual(aligner.algorithm, "Smith-Waterman")
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  wildcard: None
  match_score: 1.000000
  mismatch_score: 0.000000
  open_internal_insertion_score: -0.100000
  extend_internal_insertion_score: -0.100000
  open_left_insertion_score: -0.100000
  extend_left_insertion_score: -0.100000
  open_right_insertion_score: -0.100000
  extend_right_insertion_score: -0.100000
  open_internal_deletion_score: -0.100000
  extend_internal_deletion_score: -0.100000
  open_left_deletion_score: -0.100000
  extend_left_deletion_score: -0.100000
  open_right_deletion_score: -0.100000
  extend_right_deletion_score: -0.100000
  mode: local
""",
        )
        score = aligner.score("AwBw", "zABz")
        self.assertAlmostEqual(score, 1.9)
        alignments = aligner.align("AwBw", "zABz")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=1.9) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1.9)
        self.assertEqual(
            str(alignment),
            """\
target            0 AwB 3
                  0 |-| 3
query             1 A-B 3
""",
        )
        self.assertEqual(alignment.shape, (2, 3))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 1], [2, 3]], [[1, 2], [2, 3]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 2 identities; 0 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)

    def test_gotoh_local(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = "local"
        aligner.open_gap_score = -0.1
        aligner.extend_gap_score = 0.0
        self.assertEqual(aligner.algorithm, "Gotoh local alignment algorithm")
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  wildcard: None
  match_score: 1.000000
  mismatch_score: 0.000000
  open_internal_insertion_score: -0.100000
  extend_internal_insertion_score: 0.000000
  open_left_insertion_score: -0.100000
  extend_left_insertion_score: 0.000000
  open_right_insertion_score: -0.100000
  extend_right_insertion_score: 0.000000
  open_internal_deletion_score: -0.100000
  extend_internal_deletion_score: 0.000000
  open_left_deletion_score: -0.100000
  extend_left_deletion_score: 0.000000
  open_right_deletion_score: -0.100000
  extend_right_deletion_score: 0.000000
  mode: local
""",
        )
        score = aligner.score("AwBw", "zABz")
        self.assertAlmostEqual(score, 1.9)
        alignments = aligner.align("AwBw", "zABz")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=1.9) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1.9)
        self.assertEqual(
            str(alignment),
            """\
target            0 AwB 3
                  0 |-| 3
query             1 A-B 3
""",
        )
        self.assertEqual(alignment.shape, (2, 3))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 1], [2, 3]], [[1, 2], [2, 3]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 2 identities; 0 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)


class TestUnknownCharacter(unittest.TestCase):
    def test_needlemanwunsch_simple1(self):
        seq1 = "GACT"
        seq2 = "GA?T"
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        aligner.gap_score = -1.0
        aligner.mismatch_score = -1.0
        aligner.wildcard = "?"
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 3.0)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, 3.0)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=3) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 GACT 4
                  0 ||.| 4
query             0 GA?T 4
""",
        )
        self.assertEqual(alignment.shape, (2, 4))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[0, 4]], [[0, 4]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (4 aligned letters; 3 identities; 1 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 4:
        identities = 3,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts("?")
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (4 aligned letters; 3 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 4:
        identities = 3,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 3.0; substitution score = 3.0; gap score = 0.0; 4 aligned letters; 3 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 3.0:
        substitution_score = 3.0,
        gap_score = 0.0.
    aligned = 4:
        identities = 3,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 3.0)
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=3) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 GACT 4
                  0 ||.| 4
query             4 GA?T 0
""",
        )
        self.assertEqual(alignment.shape, (2, 4))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[0, 4]], [[4, 0]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (4 aligned letters; 3 identities; 1 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 4:
        identities = 3,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts("?")
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (4 aligned letters; 3 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 4:
        identities = 3,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 3.0; substitution score = 3.0; gap score = 0.0; 4 aligned letters; 3 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 3.0:
        substitution_score = 3.0,
        gap_score = 0.0.
    aligned = 4:
        identities = 3,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 3.0)
        seq2 = "GAXT"
        aligner.wildcard = "X"
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 3.0)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, 3.0)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=3) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 GACT 4
                  0 ||.| 4
query             0 GAXT 4
""",
        )
        self.assertEqual(alignment.shape, (2, 4))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[0, 4]], [[0, 4]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (4 aligned letters; 3 identities; 1 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 4:
        identities = 3,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts("X")
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (4 aligned letters; 3 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 4:
        identities = 3,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 3.0; substitution score = 3.0; gap score = 0.0; 4 aligned letters; 3 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 3.0:
        substitution_score = 3.0,
        gap_score = 0.0.
    aligned = 4:
        identities = 3,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 3.0)
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=3) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 GACT 4
                  0 ||.| 4
query             4 GAXT 0
""",
        )
        self.assertEqual(alignment.shape, (2, 4))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[0, 4]], [[4, 0]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (4 aligned letters; 3 identities; 1 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 4:
        identities = 3,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts("X")
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (4 aligned letters; 3 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 4:
        identities = 3,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 3.0; substitution score = 3.0; gap score = 0.0; 4 aligned letters; 3 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 3.0:
        substitution_score = 3.0,
        gap_score = 0.0.
    aligned = 4:
        identities = 3,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 3.0)
        aligner.wildcard = None
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 2.0)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, 2.0)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=2) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 2.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 GACT 4
                  0 ||.| 4
query             0 GAXT 4
""",
        )
        self.assertEqual(alignment.shape, (2, 4))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[0, 4]], [[0, 4]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (4 aligned letters; 3 identities; 1 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 4:
        identities = 3,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts("X")
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (4 aligned letters; 3 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 4:
        identities = 3,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 2.0; substitution score = 2.0; gap score = 0.0; 4 aligned letters; 3 identities; 1 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 2.0:
        substitution_score = 2.0,
        gap_score = 0.0.
    aligned = 4:
        identities = 3,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 2.0)
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=2) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 2.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 GACT 4
                  0 ||.| 4
query             4 GAXT 0
""",
        )
        self.assertEqual(alignment.shape, (2, 4))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[0, 4]], [[4, 0]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (4 aligned letters; 3 identities; 1 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 4:
        identities = 3,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts("X")
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (4 aligned letters; 3 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 4:
        identities = 3,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 2.0; substitution score = 2.0; gap score = 0.0; 4 aligned letters; 3 identities; 1 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 2.0:
        substitution_score = 2.0,
        gap_score = 0.0.
    aligned = 4:
        identities = 3,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 2.0)

    def test_needlemanwunsch_simple2(self):
        seq1 = "GA?AT"
        seq2 = "GAA?T"
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        aligner.wildcard = "?"
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 3.0)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, 3.0)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=3) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 GA?AT 5
                  0 ||..| 5
query             0 GAA?T 5
""",
        )
        self.assertEqual(alignment.shape, (2, 5))
        self.assertTrue(
            np.array_equal(
                alignment.aligned,
                np.array([[[0, 5]], [[0, 5]]]),
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (5 aligned letters; 3 identities; 2 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 5:
        identities = 3,
        mismatches = 2.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 2)
        self.assertIsNone(counts.score)
        counts = alignment.counts("?")
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (5 aligned letters; 3 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 5:
        identities = 3,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 3.0; substitution score = 3.0; gap score = 0.0; 5 aligned letters; 3 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 3.0:
        substitution_score = 3.0,
        gap_score = 0.0.
    aligned = 5:
        identities = 3,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 3.0)
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=3) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 GA?AT 5
                  0 ||..| 5
query             5 GAA?T 0
""",
        )
        self.assertEqual(alignment.shape, (2, 5))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[0, 5]], [[5, 0]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (5 aligned letters; 3 identities; 2 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 5:
        identities = 3,
        mismatches = 2.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 2)
        self.assertIsNone(counts.score)
        counts = alignment.counts("?")
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (5 aligned letters; 3 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 5:
        identities = 3,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 3.0; substitution score = 3.0; gap score = 0.0; 5 aligned letters; 3 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 3.0:
        substitution_score = 3.0,
        gap_score = 0.0.
    aligned = 5:
        identities = 3,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 3.0)
        seq1 = "GAXAT"
        seq2 = "GAAXT"
        aligner.wildcard = "X"
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 3.0)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, 3.0)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=3) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 GAXAT 5
                  0 ||..| 5
query             0 GAAXT 5
""",
        )
        self.assertEqual(alignment.shape, (2, 5))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[0, 5]], [[0, 5]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (5 aligned letters; 3 identities; 2 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 5:
        identities = 3,
        mismatches = 2.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 2)
        self.assertIsNone(counts.score)
        counts = alignment.counts("?")
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (5 aligned letters; 3 identities; 2 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 5:
        identities = 3,
        mismatches = 2.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 2)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 3.0; substitution score = 3.0; gap score = 0.0; 5 aligned letters; 3 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 3.0:
        substitution_score = 3.0,
        gap_score = 0.0.
    aligned = 5:
        identities = 3,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 3.0)
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=3) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 GAXAT 5
                  0 ||..| 5
query             5 GAAXT 0
""",
        )
        self.assertEqual(alignment.shape, (2, 5))
        self.assertTrue(
            np.array_equal(
                alignment.aligned,
                np.array([[[0, 5]], [[5, 0]]]),
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (5 aligned letters; 3 identities; 2 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 5:
        identities = 3,
        mismatches = 2.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 2)
        self.assertIsNone(counts.score)
        counts = alignment.counts("?")
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (5 aligned letters; 3 identities; 2 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 5:
        identities = 3,
        mismatches = 2.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 2)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 3.0; substitution score = 3.0; gap score = 0.0; 5 aligned letters; 3 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 3.0:
        substitution_score = 3.0,
        gap_score = 0.0.
    aligned = 5:
        identities = 3,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 3.0)

    def test_fogsaa_simple2(self):
        seq1 = "GA?AT"
        seq2 = "GAA?T"
        aligner = Align.PairwiseAligner()
        aligner.mode = "fogsaa"
        aligner.wildcard = "?"
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 3.0)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, 3.0)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=3) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 GA?AT 5
                  0 ||..| 5
query             0 GAA?T 5
""",
        )
        self.assertEqual(alignment.shape, (2, 5))
        self.assertTrue(
            np.array_equal(
                alignment.aligned,
                np.array([[[0, 5]], [[0, 5]]]),
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (5 aligned letters; 3 identities; 2 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 5:
        identities = 3,
        mismatches = 2.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 2)
        self.assertIsNone(counts.score)
        counts = alignment.counts("?")
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (5 aligned letters; 3 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 5:
        identities = 3,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 3.0; substitution score = 3.0; gap score = 0.0; 5 aligned letters; 3 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 3.0:
        substitution_score = 3.0,
        gap_score = 0.0.
    aligned = 5:
        identities = 3,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 3.0)
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=3) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 GA?AT 5
                  0 ||..| 5
query             5 GAA?T 0
""",
        )
        self.assertEqual(alignment.shape, (2, 5))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[0, 5]], [[5, 0]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (5 aligned letters; 3 identities; 2 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 5:
        identities = 3,
        mismatches = 2.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 2)
        self.assertIsNone(counts.score)
        counts = alignment.counts("?")
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (5 aligned letters; 3 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 5:
        identities = 3,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 3.0; substitution score = 3.0; gap score = 0.0; 5 aligned letters; 3 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 3.0:
        substitution_score = 3.0,
        gap_score = 0.0.
    aligned = 5:
        identities = 3,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 3.0)
        seq1 = "GAXAT"
        seq2 = "GAAXT"
        aligner.wildcard = "X"
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 3.0)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, 3.0)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=3) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 GAXAT 5
                  0 ||..| 5
query             0 GAAXT 5
""",
        )
        self.assertEqual(alignment.shape, (2, 5))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[0, 5]], [[0, 5]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (5 aligned letters; 3 identities; 2 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 5:
        identities = 3,
        mismatches = 2.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 2)
        self.assertIsNone(counts.score)
        counts = alignment.counts("X")
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (5 aligned letters; 3 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 5:
        identities = 3,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 3.0; substitution score = 3.0; gap score = 0.0; 5 aligned letters; 3 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 3.0:
        substitution_score = 3.0,
        gap_score = 0.0.
    aligned = 5:
        identities = 3,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 3.0)
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=3) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 GAXAT 5
                  0 ||..| 5
query             5 GAAXT 0
""",
        )
        self.assertEqual(alignment.shape, (2, 5))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[0, 5]], [[5, 0]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (5 aligned letters; 3 identities; 2 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 5:
        identities = 3,
        mismatches = 2.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 2)
        self.assertIsNone(counts.score)
        counts = alignment.counts("X")
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (5 aligned letters; 3 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 5:
        identities = 3,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 3.0; substitution score = 3.0; gap score = 0.0; 5 aligned letters; 3 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 3.0:
        substitution_score = 3.0,
        gap_score = 0.0.
    aligned = 5:
        identities = 3,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 3.0)


class TestPairwiseOpenPenalty(unittest.TestCase):
    def test_match_score_open_penalty1(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -0.1
        aligner.extend_gap_score = 0.0
        self.assertEqual(aligner.algorithm, "Gotoh global alignment algorithm")
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  wildcard: None
  match_score: 2.000000
  mismatch_score: -1.000000
  open_internal_insertion_score: -0.100000
  extend_internal_insertion_score: 0.000000
  open_left_insertion_score: -0.100000
  extend_left_insertion_score: 0.000000
  open_right_insertion_score: -0.100000
  extend_right_insertion_score: 0.000000
  open_internal_deletion_score: -0.100000
  extend_internal_deletion_score: 0.000000
  open_left_deletion_score: -0.100000
  extend_left_deletion_score: 0.000000
  open_right_deletion_score: -0.100000
  extend_right_deletion_score: 0.000000
  mode: global
""",
        )
        seq1 = "AA"
        seq2 = "A"
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 1.9)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, 1.9)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (2 alignments; score=1.9) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1.9)
        self.assertEqual(
            str(alignment),
            """\
target            0 AA 2
                  0 -| 2
query             0 -A 1
""",
        )
        self.assertEqual(alignment.shape, (2, 2))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[1, 2]], [[0, 1]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (1 aligned letters; 1 identities; 0 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 1:
        identities = 1,
        mismatches = 0.
    gaps = 1:
        left_gaps = 1:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 1:
                open_left_deletions = 1,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 1)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.9; substitution score = 2.0; gap score = -0.1; 1 aligned letters; 1 identities; 0 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.9:
        substitution_score = 2.0,
        gap_score = -0.1.
    aligned = 1:
        identities = 1,
        mismatches = 0.
    gaps = 1:
        left_gaps = 1:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 1:
                open_left_deletions = 1,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 1)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 1.9)
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 1.9)
        self.assertEqual(
            str(alignment),
            """\
target            0 AA 2
                  0 |- 2
query             0 A- 1
""",
        )
        self.assertEqual(alignment.shape, (2, 2))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[0, 1]], [[0, 1]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (1 aligned letters; 1 identities; 0 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 1:
        identities = 1,
        mismatches = 0.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 1:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 1:
                open_right_deletions = 1,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 1)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.9; substitution score = 2.0; gap score = -0.1; 1 aligned letters; 1 identities; 0 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.9:
        substitution_score = 2.0,
        gap_score = -0.1.
    aligned = 1:
        identities = 1,
        mismatches = 0.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 1:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 1:
                open_right_deletions = 1,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 1)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 1.9)
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (2 alignments; score=1.9) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1.9)
        self.assertEqual(
            str(alignment),
            """\
target            0 AA 2
                  0 -| 2
query             1 -A 0
""",
        )
        self.assertEqual(alignment.shape, (2, 2))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[1, 2]], [[1, 0]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (1 aligned letters; 1 identities; 0 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 1:
        identities = 1,
        mismatches = 0.
    gaps = 1:
        left_gaps = 1:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 1:
                open_left_deletions = 1,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 1)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.9; substitution score = 2.0; gap score = -0.1; 1 aligned letters; 1 identities; 0 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.9:
        substitution_score = 2.0,
        gap_score = -0.1.
    aligned = 1:
        identities = 1,
        mismatches = 0.
    gaps = 1:
        left_gaps = 1:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 1:
                open_left_deletions = 1,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 1)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 1.9)
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 1.9)
        self.assertEqual(
            str(alignment),
            """\
target            0 AA 2
                  0 |- 2
query             1 A- 0
""",
        )
        self.assertEqual(alignment.shape, (2, 2))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[0, 1]], [[1, 0]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (1 aligned letters; 1 identities; 0 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 1:
        identities = 1,
        mismatches = 0.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 1:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 1:
                open_right_deletions = 1,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 1)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.9; substitution score = 2.0; gap score = -0.1; 1 aligned letters; 1 identities; 0 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.9:
        substitution_score = 2.0,
        gap_score = -0.1.
    aligned = 1:
        identities = 1,
        mismatches = 0.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 1:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 1:
                open_right_deletions = 1,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 1)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 1.9)

    def test_match_score_open_penalty2(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score = 1.5
        aligner.mismatch_score = 0.0
        aligner.open_gap_score = -0.1
        aligner.extend_gap_score = 0.0
        self.assertEqual(aligner.algorithm, "Gotoh global alignment algorithm")
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  wildcard: None
  match_score: 1.500000
  mismatch_score: 0.000000
  open_internal_insertion_score: -0.100000
  extend_internal_insertion_score: 0.000000
  open_left_insertion_score: -0.100000
  extend_left_insertion_score: 0.000000
  open_right_insertion_score: -0.100000
  extend_right_insertion_score: 0.000000
  open_internal_deletion_score: -0.100000
  extend_internal_deletion_score: 0.000000
  open_left_deletion_score: -0.100000
  extend_left_deletion_score: 0.000000
  open_right_deletion_score: -0.100000
  extend_right_deletion_score: 0.000000
  mode: global
""",
        )
        seq1 = "GAA"
        seq2 = "GA"
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 2.9)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, 2.9)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (2 alignments; score=2.9) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 2.9)
        self.assertEqual(
            str(alignment),
            """\
target            0 GAA 3
                  0 |-| 3
query             0 G-A 2
""",
        )
        self.assertEqual(alignment.shape, (2, 3))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 1], [2, 3]], [[0, 1], [1, 2]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 2 identities; 0 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 2.9; substitution score = 3.0; gap score = -0.1; 2 aligned letters; 2 identities; 0 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 2.9:
        substitution_score = 3.0,
        gap_score = -0.1.
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 2.9)
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 2.9)
        self.assertEqual(
            str(alignment),
            """\
target            0 GAA 3
                  0 ||- 3
query             0 GA- 2
""",
        )
        self.assertEqual(alignment.shape, (2, 3))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[0, 2]], [[0, 2]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 2 identities; 0 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 1:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 1:
                open_right_deletions = 1,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 2.9; substitution score = 3.0; gap score = -0.1; 2 aligned letters; 2 identities; 0 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 2.9:
        substitution_score = 3.0,
        gap_score = -0.1.
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 1:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 1:
                open_right_deletions = 1,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 2.9)
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (2 alignments; score=2.9) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 2.9)
        self.assertEqual(
            str(alignment),
            """\
target            0 GAA 3
                  0 |-| 3
query             2 G-A 0
""",
        )
        self.assertEqual(alignment.shape, (2, 3))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 1], [2, 3]], [[2, 1], [1, 0]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 2 identities; 0 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 2.9; substitution score = 3.0; gap score = -0.1; 2 aligned letters; 2 identities; 0 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 2.9:
        substitution_score = 3.0,
        gap_score = -0.1.
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 2.9)
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 2.9)
        self.assertEqual(
            str(alignment),
            """\
target            0 GAA 3
                  0 ||- 3
query             2 GA- 0
""",
        )
        self.assertEqual(alignment.shape, (2, 3))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[0, 2]], [[2, 0]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 2 identities; 0 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 1:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 1:
                open_right_deletions = 1,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 2.9; substitution score = 3.0; gap score = -0.1; 2 aligned letters; 2 identities; 0 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 2.9:
        substitution_score = 3.0,
        gap_score = -0.1.
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 1:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 1:
                open_right_deletions = 1,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 2.9)

    def test_match_score_open_penalty3(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        aligner.open_deletion_score = -0.1
        aligner.extend_deletion_score = 0.0
        self.assertEqual(aligner.algorithm, "Gotoh global alignment algorithm")
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  wildcard: None
  match_score: 1.000000
  mismatch_score: 0.000000
  open_internal_insertion_score: 0.000000
  extend_internal_insertion_score: 0.000000
  open_left_insertion_score: 0.000000
  extend_left_insertion_score: 0.000000
  open_right_insertion_score: 0.000000
  extend_right_insertion_score: 0.000000
  open_internal_deletion_score: -0.100000
  extend_internal_deletion_score: 0.000000
  open_left_deletion_score: -0.100000
  extend_left_deletion_score: 0.000000
  open_right_deletion_score: -0.100000
  extend_right_deletion_score: 0.000000
  mode: global
""",
        )
        seq1 = "GAACT"
        seq2 = "GAT"
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 2.9)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, 2.9)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=2.9) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 2.9)
        self.assertEqual(
            str(alignment),
            """\
target            0 GAACT 5
                  0 ||--| 5
query             0 GA--T 3
""",
        )
        self.assertEqual(alignment.shape, (2, 5))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 2], [4, 5]], [[0, 2], [2, 3]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (3 aligned letters; 3 identities; 0 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 3:
        identities = 3,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 2.9; substitution score = 3.0; gap score = -0.1; 3 aligned letters; 3 identities; 0 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 2.9:
        substitution_score = 3.0,
        gap_score = -0.1.
    aligned = 3:
        identities = 3,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 2.9)
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=2.9) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 2.9)
        self.assertEqual(
            str(alignment),
            """\
target            0 GAACT 5
                  0 ||--| 5
query             3 GA--T 0
""",
        )
        self.assertEqual(alignment.shape, (2, 5))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 2], [4, 5]], [[3, 1], [1, 0]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (3 aligned letters; 3 identities; 0 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 3:
        identities = 3,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 2.9; substitution score = 3.0; gap score = -0.1; 3 aligned letters; 3 identities; 0 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 2.9:
        substitution_score = 3.0,
        gap_score = -0.1.
    aligned = 3:
        identities = 3,
        mismatches = 0.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 2.9)

    def test_match_score_open_penalty3_fogsaa(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = "fogsaa"
        aligner.open_deletion_score = -0.1
        aligner.extend_deletion_score = 0.0
        self.assertEqual(
            aligner.algorithm, "Fast Optimal Global Sequence Alignment Algorithm"
        )
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  wildcard: None
  match_score: 1.000000
  mismatch_score: 0.000000
  open_internal_insertion_score: 0.000000
  extend_internal_insertion_score: 0.000000
  open_left_insertion_score: 0.000000
  extend_left_insertion_score: 0.000000
  open_right_insertion_score: 0.000000
  extend_right_insertion_score: 0.000000
  open_internal_deletion_score: -0.100000
  extend_internal_deletion_score: 0.000000
  open_left_deletion_score: -0.100000
  extend_left_deletion_score: 0.000000
  open_right_deletion_score: -0.100000
  extend_right_deletion_score: 0.000000
  mode: fogsaa
""",
        )
        seq1 = "GAACT"
        seq2 = "GAT"
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 2.9)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, 2.9)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=2.9) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 2.9)
        self.assertEqual(
            str(alignment),
            """\
target            0 GAACT 5
                  0 ||--| 5
query             0 GA--T 3
""",
        )
        self.assertEqual(alignment.shape, (2, 5))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 2], [4, 5]], [[0, 2], [2, 3]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (3 aligned letters; 3 identities; 0 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 3:
        identities = 3,
        mismatches = 0.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 2.9; substitution score = 3.0; gap score = -0.1; 3 aligned letters; 3 identities; 0 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 2.9:
        substitution_score = 3.0,
        gap_score = -0.1.
    aligned = 3:
        identities = 3,
        mismatches = 0.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 2.9)
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=2.9) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 2.9)
        self.assertEqual(
            str(alignment),
            """\
target            0 GAACT 5
                  0 ||--| 5
query             3 GA--T 0
""",
        )
        self.assertEqual(alignment.shape, (2, 5))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 2], [4, 5]], [[3, 1], [1, 0]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (3 aligned letters; 3 identities; 0 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 3:
        identities = 3,
        mismatches = 0.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 2.9; substitution score = 3.0; gap score = -0.1; 3 aligned letters; 3 identities; 0 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 2.9:
        substitution_score = 3.0,
        gap_score = -0.1.
    aligned = 3:
        identities = 3,
        mismatches = 0.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 2.9)

    def test_match_score_open_penalty4(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        aligner.mismatch_score = -2.0
        aligner.open_gap_score = -0.1
        aligner.extend_gap_score = 0.0
        self.assertEqual(aligner.algorithm, "Gotoh global alignment algorithm")
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  wildcard: None
  match_score: 1.000000
  mismatch_score: -2.000000
  open_internal_insertion_score: -0.100000
  extend_internal_insertion_score: 0.000000
  open_left_insertion_score: -0.100000
  extend_left_insertion_score: 0.000000
  open_right_insertion_score: -0.100000
  extend_right_insertion_score: 0.000000
  open_internal_deletion_score: -0.100000
  extend_internal_deletion_score: 0.000000
  open_left_deletion_score: -0.100000
  extend_left_deletion_score: 0.000000
  open_right_deletion_score: -0.100000
  extend_right_deletion_score: 0.000000
  mode: global
""",
        )
        seq1 = "GCT"
        seq2 = "GATA"
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 1.7)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, 1.7)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (2 alignments; score=1.7) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1.7)
        self.assertEqual(
            str(alignment),
            """\
target            0 G-CT- 3
                  0 |--|- 5
query             0 GA-TA 4
""",
        )
        self.assertEqual(alignment.shape, (2, 5))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 1], [2, 3]], [[0, 1], [2, 3]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 2 identities; 0 mismatches; 3 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 3:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 1:
            right_insertions = 1:
                open_right_insertions = 1,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.7; substitution score = 2.0; gap score = -0.30000000000000004; 2 aligned letters; 2 identities; 0 mismatches; 3 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.7:
        substitution_score = 2.0,
        gap_score = -0.30000000000000004.
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 3:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 1:
            right_insertions = 1:
                open_right_insertions = 1,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 1.7)
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 1.7)
        self.assertEqual(
            str(alignment),
            """\
target            0 GC-T- 3
                  0 |--|- 5
query             0 G-ATA 4
""",
        )
        self.assertEqual(alignment.shape, (2, 5))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 1], [2, 3]], [[0, 1], [2, 3]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 2 identities; 0 mismatches; 3 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 3:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 1:
            right_insertions = 1:
                open_right_insertions = 1,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.7; substitution score = 2.0; gap score = -0.30000000000000004; 2 aligned letters; 2 identities; 0 mismatches; 3 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.7:
        substitution_score = 2.0,
        gap_score = -0.30000000000000004.
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 3:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 1:
            right_insertions = 1:
                open_right_insertions = 1,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 1.7)
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (2 alignments; score=1.7) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1.7)
        self.assertEqual(
            str(alignment),
            """\
target            0 G-CT- 3
                  0 |--|- 5
query             4 GA-TA 0
""",
        )
        self.assertEqual(alignment.shape, (2, 5))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 1], [2, 3]], [[4, 3], [2, 1]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 2 identities; 0 mismatches; 3 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 3:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 1:
            right_insertions = 1:
                open_right_insertions = 1,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.7; substitution score = 2.0; gap score = -0.30000000000000004; 2 aligned letters; 2 identities; 0 mismatches; 3 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.7:
        substitution_score = 2.0,
        gap_score = -0.30000000000000004.
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 3:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 1:
            right_insertions = 1:
                open_right_insertions = 1,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 1.7)
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 1.7)
        self.assertEqual(
            str(alignment),
            """\
target            0 GC-T- 3
                  0 |--|- 5
query             4 G-ATA 0
""",
        )
        self.assertEqual(alignment.shape, (2, 5))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 1], [2, 3]], [[4, 3], [2, 1]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 2 identities; 0 mismatches; 3 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 3:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 1:
            right_insertions = 1:
                open_right_insertions = 1,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.7; substitution score = 2.0; gap score = -0.30000000000000004; 2 aligned letters; 2 identities; 0 mismatches; 3 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.7:
        substitution_score = 2.0,
        gap_score = -0.30000000000000004.
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 3:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 1:
            right_insertions = 1:
                open_right_insertions = 1,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 1.7)


class TestPairwiseExtendPenalty(unittest.TestCase):
    def test_extend_penalty1(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        aligner.open_gap_score = -0.2
        aligner.extend_gap_score = -0.5
        self.assertEqual(aligner.algorithm, "Gotoh global alignment algorithm")
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  wildcard: None
  match_score: 1.000000
  mismatch_score: 0.000000
  open_internal_insertion_score: -0.200000
  extend_internal_insertion_score: -0.500000
  open_left_insertion_score: -0.200000
  extend_left_insertion_score: -0.500000
  open_right_insertion_score: -0.200000
  extend_right_insertion_score: -0.500000
  open_internal_deletion_score: -0.200000
  extend_internal_deletion_score: -0.500000
  open_left_deletion_score: -0.200000
  extend_left_deletion_score: -0.500000
  open_right_deletion_score: -0.200000
  extend_right_deletion_score: -0.500000
  mode: global
""",
        )
        seq1 = "GACT"
        seq2 = "GT"
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 1.3)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, 1.3)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=1.3) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1.3)
        self.assertEqual(
            str(alignment),
            """\
target            0 GACT 4
                  0 |--| 4
query             0 G--T 2
""",
        )
        self.assertEqual(alignment.shape, (2, 4))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 1], [3, 4]], [[0, 1], [1, 2]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 2 identities; 0 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.3; substitution score = 2.0; gap score = -0.7; 2 aligned letters; 2 identities; 0 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.3:
        substitution_score = 2.0,
        gap_score = -0.7.
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 1.3)
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=1.3) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1.3)
        self.assertEqual(
            str(alignment),
            """\
target            0 GACT 4
                  0 |--| 4
query             2 G--T 0
""",
        )
        self.assertEqual(alignment.shape, (2, 4))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 1], [3, 4]], [[2, 1], [1, 0]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 2 identities; 0 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.3; substitution score = 2.0; gap score = -0.7; 2 aligned letters; 2 identities; 0 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.3:
        substitution_score = 2.0,
        gap_score = -0.7.
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 1.3)

    def test_extend_penalty2(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        aligner.open_gap_score = -0.2
        aligner.extend_gap_score = -1.5
        self.assertEqual(aligner.algorithm, "Gotoh global alignment algorithm")
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  wildcard: None
  match_score: 1.000000
  mismatch_score: 0.000000
  open_internal_insertion_score: -0.200000
  extend_internal_insertion_score: -1.500000
  open_left_insertion_score: -0.200000
  extend_left_insertion_score: -1.500000
  open_right_insertion_score: -0.200000
  extend_right_insertion_score: -1.500000
  open_internal_deletion_score: -0.200000
  extend_internal_deletion_score: -1.500000
  open_left_deletion_score: -0.200000
  extend_left_deletion_score: -1.500000
  open_right_deletion_score: -0.200000
  extend_right_deletion_score: -1.500000
  mode: global
""",
        )
        seq1 = "GACT"
        seq2 = "GT"
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 0.6)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, 0.6)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (2 alignments; score=0.6) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 0.6)
        self.assertEqual(
            str(alignment),
            """\
target            0 GACT 4
                  0 -.-| 4
query             0 -G-T 2
""",
        )
        self.assertEqual(alignment.shape, (2, 4))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[1, 2], [3, 4]], [[0, 1], [1, 2]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 1 identities; 1 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 1,
        mismatches = 1.
    gaps = 2:
        left_gaps = 1:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 1:
                open_left_deletions = 1,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 0.6; substitution score = 1.0; gap score = -0.4; 2 aligned letters; 1 identities; 1 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 0.6:
        substitution_score = 1.0,
        gap_score = -0.4.
    aligned = 2:
        identities = 1,
        mismatches = 1.
    gaps = 2:
        left_gaps = 1:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 1:
                open_left_deletions = 1,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 0.6)
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 0.6)
        self.assertEqual(
            str(alignment),
            """\
target            0 GACT 4
                  0 |-.- 4
query             0 G-T- 2
""",
        )
        self.assertEqual(alignment.shape, (2, 4))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 1], [2, 3]], [[0, 1], [1, 2]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 1 identities; 1 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 1,
        mismatches = 1.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 1:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 1:
                open_right_deletions = 1,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 0.6; substitution score = 1.0; gap score = -0.4; 2 aligned letters; 1 identities; 1 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 0.6:
        substitution_score = 1.0,
        gap_score = -0.4.
    aligned = 2:
        identities = 1,
        mismatches = 1.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 1:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 1:
                open_right_deletions = 1,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 0.6)
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (2 alignments; score=0.6) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 0.6)
        self.assertEqual(
            str(alignment),
            """\
target            0 GACT 4
                  0 -.-| 4
query             2 -G-T 0
""",
        )
        self.assertEqual(alignment.shape, (2, 4))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[1, 2], [3, 4]], [[2, 1], [1, 0]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 1 identities; 1 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 1,
        mismatches = 1.
    gaps = 2:
        left_gaps = 1:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 1:
                open_left_deletions = 1,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 0.6; substitution score = 1.0; gap score = -0.4; 2 aligned letters; 1 identities; 1 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 0.6:
        substitution_score = 1.0,
        gap_score = -0.4.
    aligned = 2:
        identities = 1,
        mismatches = 1.
    gaps = 2:
        left_gaps = 1:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 1:
                open_left_deletions = 1,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 0.6)
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 0.6)
        self.assertEqual(
            str(alignment),
            """\
target            0 GACT 4
                  0 |-.- 4
query             2 G-T- 0
""",
        )
        self.assertEqual(alignment.shape, (2, 4))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 1], [2, 3]], [[2, 1], [1, 0]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 1 identities; 1 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 1,
        mismatches = 1.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 1:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 1:
                open_right_deletions = 1,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 0.6; substitution score = 1.0; gap score = -0.4; 2 aligned letters; 1 identities; 1 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 0.6:
        substitution_score = 1.0,
        gap_score = -0.4.
    aligned = 2:
        identities = 1,
        mismatches = 1.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 1:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 1:
                open_right_deletions = 1,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 0.6)

    def test_extend_penalty2_fogsaa(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = "fogsaa"
        aligner.open_gap_score = -0.2
        aligner.extend_gap_score = -1.5
        self.assertEqual(
            aligner.algorithm, "Fast Optimal Global Sequence Alignment Algorithm"
        )
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  wildcard: None
  match_score: 1.000000
  mismatch_score: 0.000000
  open_internal_insertion_score: -0.200000
  extend_internal_insertion_score: -1.500000
  open_left_insertion_score: -0.200000
  extend_left_insertion_score: -1.500000
  open_right_insertion_score: -0.200000
  extend_right_insertion_score: -1.500000
  open_internal_deletion_score: -0.200000
  extend_internal_deletion_score: -1.500000
  open_left_deletion_score: -0.200000
  extend_left_deletion_score: -1.500000
  open_right_deletion_score: -0.200000
  extend_right_deletion_score: -1.500000
  mode: fogsaa
""",
        )
        seq1 = "GACT"
        seq2 = "GT"
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 0.6)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, 0.6)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=0.6) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 0.6)
        self.assertEqual(
            str(alignment),
            """\
target            0 GACT 4
                  0 -.-| 4
query             0 -G-T 2
""",
        )
        self.assertEqual(alignment.shape, (2, 4))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[1, 2], [3, 4]], [[0, 1], [1, 2]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 1 identities; 1 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 1,
        mismatches = 1.
    gaps = 2:
        left_gaps = 1:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 1:
                open_left_deletions = 1,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 0.6; substitution score = 1.0; gap score = -0.4; 2 aligned letters; 1 identities; 1 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 0.6:
        substitution_score = 1.0,
        gap_score = -0.4.
    aligned = 2:
        identities = 1,
        mismatches = 1.
    gaps = 2:
        left_gaps = 1:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 1:
                open_left_deletions = 1,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 0.6)
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=0.6) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 0.6)
        self.assertEqual(
            str(alignment),
            """\
target            0 GACT 4
                  0 -.-| 4
query             2 -G-T 0
""",
        )
        self.assertEqual(alignment.shape, (2, 4))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[1, 2], [3, 4]], [[2, 1], [1, 0]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 1 identities; 1 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 1,
        mismatches = 1.
    gaps = 2:
        left_gaps = 1:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 1:
                open_left_deletions = 1,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 0.6; substitution score = 1.0; gap score = -0.4; 2 aligned letters; 1 identities; 1 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 0.6:
        substitution_score = 1.0,
        gap_score = -0.4.
    aligned = 2:
        identities = 1,
        mismatches = 1.
    gaps = 2:
        left_gaps = 1:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 1:
                open_left_deletions = 1,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 0.6)


class TestPairwisePenalizeExtendWhenOpening(unittest.TestCase):
    def test_penalize_extend_when_opening(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        aligner.open_gap_score = -1.7
        aligner.extend_gap_score = -1.5
        self.assertEqual(aligner.algorithm, "Gotoh global alignment algorithm")
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  wildcard: None
  match_score: 1.000000
  mismatch_score: 0.000000
  open_internal_insertion_score: -1.700000
  extend_internal_insertion_score: -1.500000
  open_left_insertion_score: -1.700000
  extend_left_insertion_score: -1.500000
  open_right_insertion_score: -1.700000
  extend_right_insertion_score: -1.500000
  open_internal_deletion_score: -1.700000
  extend_internal_deletion_score: -1.500000
  open_left_deletion_score: -1.700000
  extend_left_deletion_score: -1.500000
  open_right_deletion_score: -1.700000
  extend_right_deletion_score: -1.500000
  mode: global
""",
        )
        seq1 = "GACT"
        seq2 = "GT"
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, -1.2)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, -1.2)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=-1.2) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, -1.2)
        self.assertEqual(
            str(alignment),
            """\
target            0 GACT 4
                  0 |--| 4
query             0 G--T 2
""",
        )
        self.assertEqual(alignment.shape, (2, 4))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 1], [3, 4]], [[0, 1], [1, 2]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 2 identities; 0 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = -1.2000000000000002; substitution score = 2.0; gap score = -3.2; 2 aligned letters; 2 identities; 0 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = -1.2000000000000002:
        substitution_score = 2.0,
        gap_score = -3.2.
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, -1.2)
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=-1.2) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, -1.2)
        self.assertEqual(
            str(alignment),
            """\
target            0 GACT 4
                  0 |--| 4
query             2 G--T 0
""",
        )
        self.assertEqual(alignment.shape, (2, 4))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 1], [3, 4]], [[2, 1], [1, 0]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 2 identities; 0 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = -1.2000000000000002; substitution score = 2.0; gap score = -3.2; 2 aligned letters; 2 identities; 0 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = -1.2000000000000002:
        substitution_score = 2.0,
        gap_score = -3.2.
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, -1.2)

    def test_penalize_extend_when_opening_fogsaa(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = "fogsaa"
        aligner.open_gap_score = -1.7
        aligner.extend_gap_score = -1.5
        self.assertEqual(
            aligner.algorithm, "Fast Optimal Global Sequence Alignment Algorithm"
        )
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  wildcard: None
  match_score: 1.000000
  mismatch_score: 0.000000
  open_internal_insertion_score: -1.700000
  extend_internal_insertion_score: -1.500000
  open_left_insertion_score: -1.700000
  extend_left_insertion_score: -1.500000
  open_right_insertion_score: -1.700000
  extend_right_insertion_score: -1.500000
  open_internal_deletion_score: -1.700000
  extend_internal_deletion_score: -1.500000
  open_left_deletion_score: -1.700000
  extend_left_deletion_score: -1.500000
  open_right_deletion_score: -1.700000
  extend_right_deletion_score: -1.500000
  mode: fogsaa
""",
        )
        seq1 = "GACT"
        seq2 = "GT"
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, -1.2)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, -1.2)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=-1.2) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, -1.2)
        self.assertEqual(
            str(alignment),
            """\
target            0 GACT 4
                  0 |--| 4
query             0 G--T 2
""",
        )
        self.assertEqual(alignment.shape, (2, 4))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 1], [3, 4]], [[0, 1], [1, 2]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 2 identities; 0 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = -1.2000000000000002; substitution score = 2.0; gap score = -3.2; 2 aligned letters; 2 identities; 0 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = -1.2000000000000002:
        substitution_score = 2.0,
        gap_score = -3.2.
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, -1.2)
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=-1.2) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, -1.2)
        self.assertEqual(
            str(alignment),
            """\
target            0 GACT 4
                  0 |--| 4
query             2 G--T 0
""",
        )
        self.assertEqual(alignment.shape, (2, 4))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 1], [3, 4]], [[2, 1], [1, 0]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 2 identities; 0 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = -1.2000000000000002; substitution score = 2.0; gap score = -3.2; 2 aligned letters; 2 identities; 0 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = -1.2000000000000002:
        substitution_score = 2.0,
        gap_score = -3.2.
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, -1.2)


class TestPairwisePenalizeEndgaps(unittest.TestCase):
    def test_penalize_end_gaps(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        aligner.open_gap_score = -0.2
        aligner.extend_gap_score = -0.8
        aligner.end_gap_score = 0.0
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  wildcard: None
  match_score: 1.000000
  mismatch_score: 0.000000
  open_internal_insertion_score: -0.200000
  extend_internal_insertion_score: -0.800000
  open_left_insertion_score: 0.000000
  extend_left_insertion_score: 0.000000
  open_right_insertion_score: 0.000000
  extend_right_insertion_score: 0.000000
  open_internal_deletion_score: -0.200000
  extend_internal_deletion_score: -0.800000
  open_left_deletion_score: 0.000000
  extend_left_deletion_score: 0.000000
  open_right_deletion_score: 0.000000
  extend_right_deletion_score: 0.000000
  mode: global
""",
        )
        self.assertEqual(aligner.algorithm, "Gotoh global alignment algorithm")
        seq1 = "GACT"
        seq2 = "GT"
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 1.0)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, 1.0)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (3 alignments; score=1) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 3)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 GACT 4
                  0 --.| 4
query             0 --GT 2
""",
        )
        self.assertEqual(alignment.shape, (2, 4))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[2, 4]], [[0, 2]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 1 identities; 1 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 1,
        mismatches = 1.
    gaps = 2:
        left_gaps = 2:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 2:
                open_left_deletions = 1,
                extend_left_deletions = 1;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.0; substitution score = 1.0; gap score = 0.0; 2 aligned letters; 1 identities; 1 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.0:
        substitution_score = 1.0,
        gap_score = 0.0.
    aligned = 2:
        identities = 1,
        mismatches = 1.
    gaps = 2:
        left_gaps = 2:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 2:
                open_left_deletions = 1,
                extend_left_deletions = 1;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 1.0)
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 1.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 GACT 4
                  0 |--| 4
query             0 G--T 2
""",
        )
        self.assertEqual(alignment.shape, (2, 4))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 1], [3, 4]], [[0, 1], [1, 2]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 2 identities; 0 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.0; substitution score = 2.0; gap score = -1.0; 2 aligned letters; 2 identities; 0 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.0:
        substitution_score = 2.0,
        gap_score = -1.0.
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 1.0)
        alignment = alignments[2]
        self.assertAlmostEqual(alignment.score, 1.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 GACT 4
                  0 |.-- 4
query             0 GT-- 2
""",
        )
        self.assertEqual(alignment.shape, (2, 4))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[0, 2]], [[0, 2]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 1 identities; 1 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 1,
        mismatches = 1.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 2:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 2:
                open_right_deletions = 1,
                extend_right_deletions = 1.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.0; substitution score = 1.0; gap score = 0.0; 2 aligned letters; 1 identities; 1 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.0:
        substitution_score = 1.0,
        gap_score = 0.0.
    aligned = 2:
        identities = 1,
        mismatches = 1.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 2:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 2:
                open_right_deletions = 1,
                extend_right_deletions = 1.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 1.0)
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (3 alignments; score=1) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 3)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 GACT 4
                  0 --.| 4
query             2 --GT 0
""",
        )
        self.assertEqual(alignment.shape, (2, 4))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[2, 4]], [[2, 0]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 1 identities; 1 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 1,
        mismatches = 1.
    gaps = 2:
        left_gaps = 2:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 2:
                open_left_deletions = 1,
                extend_left_deletions = 1;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.0; substitution score = 1.0; gap score = 0.0; 2 aligned letters; 1 identities; 1 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.0:
        substitution_score = 1.0,
        gap_score = 0.0.
    aligned = 2:
        identities = 1,
        mismatches = 1.
    gaps = 2:
        left_gaps = 2:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 2:
                open_left_deletions = 1,
                extend_left_deletions = 1;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 1.0)
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 1.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 GACT 4
                  0 |--| 4
query             2 G--T 0
""",
        )
        self.assertEqual(alignment.shape, (2, 4))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 1], [3, 4]], [[2, 1], [1, 0]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 2 identities; 0 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.0; substitution score = 2.0; gap score = -1.0; 2 aligned letters; 2 identities; 0 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.0:
        substitution_score = 2.0,
        gap_score = -1.0.
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 1.0)
        alignment = alignments[2]
        self.assertAlmostEqual(alignment.score, 1.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 GACT 4
                  0 |.-- 4
query             2 GT-- 0
""",
        )
        self.assertEqual(alignment.shape, (2, 4))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[0, 2]], [[2, 0]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 1 identities; 1 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 1,
        mismatches = 1.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 2:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 2:
                open_right_deletions = 1,
                extend_right_deletions = 1.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.0; substitution score = 1.0; gap score = 0.0; 2 aligned letters; 1 identities; 1 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.0:
        substitution_score = 1.0,
        gap_score = 0.0.
    aligned = 2:
        identities = 1,
        mismatches = 1.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 2:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 2:
                open_right_deletions = 1,
                extend_right_deletions = 1.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 1.0)

    def test_penalize_end_gaps_fogsaa(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = "fogsaa"
        aligner.open_gap_score = -0.2
        aligner.extend_gap_score = -0.8
        aligner.end_gap_score = 0.0
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  wildcard: None
  match_score: 1.000000
  mismatch_score: 0.000000
  open_internal_insertion_score: -0.200000
  extend_internal_insertion_score: -0.800000
  open_left_insertion_score: 0.000000
  extend_left_insertion_score: 0.000000
  open_right_insertion_score: 0.000000
  extend_right_insertion_score: 0.000000
  open_internal_deletion_score: -0.200000
  extend_internal_deletion_score: -0.800000
  open_left_deletion_score: 0.000000
  extend_left_deletion_score: 0.000000
  open_right_deletion_score: 0.000000
  extend_right_deletion_score: 0.000000
  mode: fogsaa
""",
        )
        self.assertEqual(
            aligner.algorithm, "Fast Optimal Global Sequence Alignment Algorithm"
        )
        seq1 = "GACT"
        seq2 = "GT"
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 1.0)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, 1.0)

        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=1) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 GACT 4
                  0 |.-- 4
query             0 GT-- 2
""",
        )
        self.assertEqual(alignment.shape, (2, 4))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[0, 2]], [[0, 2]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 1 identities; 1 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 1,
        mismatches = 1.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 2:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 2:
                open_right_deletions = 1,
                extend_right_deletions = 1.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.0; substitution score = 1.0; gap score = 0.0; 2 aligned letters; 1 identities; 1 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.0:
        substitution_score = 1.0,
        gap_score = 0.0.
    aligned = 2:
        identities = 1,
        mismatches = 1.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 2:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 2:
                open_right_deletions = 1,
                extend_right_deletions = 1.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 1.0)

        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=1) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 GACT 4
                  0 |.-- 4
query             2 GT-- 0
""",
        )
        self.assertEqual(alignment.shape, (2, 4))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[0, 2]], [[2, 0]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 1 identities; 1 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 1,
        mismatches = 1.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 2:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 2:
                open_right_deletions = 1,
                extend_right_deletions = 1.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.0; substitution score = 1.0; gap score = 0.0; 2 aligned letters; 1 identities; 1 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.0:
        substitution_score = 1.0,
        gap_score = 0.0.
    aligned = 2:
        identities = 1,
        mismatches = 1.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 2:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 2:
                open_right_deletions = 1,
                extend_right_deletions = 1.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 1.0)


class TestPairwiseSeparateGapPenalties(unittest.TestCase):
    def test_separate_gap_penalties1(self):
        seq1 = "GAT"
        seq2 = "GTCT"
        aligner = Align.PairwiseAligner()
        aligner.mode = "local"
        open_score, extend_score = (-0.3, 0)
        aligner.open_insertion_score = open_score
        aligner.extend_insertion_score = extend_score
        open_score, extend_score = (-0.8, 0)
        aligner.open_deletion_score = open_score
        aligner.extend_deletion_score = extend_score
        self.assertEqual(aligner.algorithm, "Gotoh local alignment algorithm")
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  wildcard: None
  match_score: 1.000000
  mismatch_score: 0.000000
  open_internal_insertion_score: -0.300000
  extend_internal_insertion_score: 0.000000
  open_left_insertion_score: -0.300000
  extend_left_insertion_score: 0.000000
  open_right_insertion_score: -0.300000
  extend_right_insertion_score: 0.000000
  open_internal_deletion_score: -0.800000
  extend_internal_deletion_score: 0.000000
  open_left_deletion_score: -0.800000
  extend_left_deletion_score: 0.000000
  open_right_deletion_score: -0.800000
  extend_right_deletion_score: 0.000000
  mode: local
""",
        )
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 1.7)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, 1.7)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (2 alignments; score=1.7) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1.7)
        self.assertEqual(
            str(alignment),
            """\
target            0 G-AT 3
                  0 |-.| 4
query             0 GTCT 4
""",
        )
        self.assertEqual(alignment.shape, (2, 4))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 1], [1, 3]], [[0, 1], [2, 4]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (3 aligned letters; 2 identities; 1 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 3:
        identities = 2,
        mismatches = 1.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.7; substitution score = 2.0; gap score = -0.3; 3 aligned letters; 2 identities; 1 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.7:
        substitution_score = 2.0,
        gap_score = -0.3.
    aligned = 3:
        identities = 2,
        mismatches = 1.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 1.7)
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 1.7)
        self.assertEqual(
            str(alignment),
            """\
target            0 GA-T 3
                  0 |.-| 4
query             0 GTCT 4
""",
        )
        self.assertEqual(alignment.shape, (2, 4))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 2], [2, 3]], [[0, 2], [3, 4]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (3 aligned letters; 2 identities; 1 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 3:
        identities = 2,
        mismatches = 1.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.7; substitution score = 2.0; gap score = -0.3; 3 aligned letters; 2 identities; 1 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.7:
        substitution_score = 2.0,
        gap_score = -0.3.
    aligned = 3:
        identities = 2,
        mismatches = 1.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 1.7)
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (2 alignments; score=1.7) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1.7)
        self.assertEqual(
            str(alignment),
            """\
target            0 G-AT 3
                  0 |-.| 4
query             4 GTCT 0
""",
        )
        self.assertEqual(alignment.shape, (2, 4))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 1], [1, 3]], [[4, 3], [2, 0]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (3 aligned letters; 2 identities; 1 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 3:
        identities = 2,
        mismatches = 1.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.7; substitution score = 2.0; gap score = -0.3; 3 aligned letters; 2 identities; 1 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.7:
        substitution_score = 2.0,
        gap_score = -0.3.
    aligned = 3:
        identities = 2,
        mismatches = 1.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 1.7)
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 1.7)
        self.assertEqual(
            str(alignment),
            """\
target            0 GA-T 3
                  0 |.-| 4
query             4 GTCT 0
""",
        )
        self.assertEqual(alignment.shape, (2, 4))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 2], [2, 3]], [[4, 2], [1, 0]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (3 aligned letters; 2 identities; 1 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 3:
        identities = 2,
        mismatches = 1.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.7; substitution score = 2.0; gap score = -0.3; 3 aligned letters; 2 identities; 1 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.7:
        substitution_score = 2.0,
        gap_score = -0.3.
    aligned = 3:
        identities = 2,
        mismatches = 1.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 1.7)

    def test_separate_gap_penalties2(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = "local"
        aligner.open_insertion_score = -0.3
        aligner.extend_insertion_score = 0.0
        aligner.open_deletion_score = -0.2
        aligner.extend_deletion_score = 0.0
        self.assertEqual(aligner.algorithm, "Gotoh local alignment algorithm")
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  wildcard: None
  match_score: 1.000000
  mismatch_score: 0.000000
  open_internal_insertion_score: -0.300000
  extend_internal_insertion_score: 0.000000
  open_left_insertion_score: -0.300000
  extend_left_insertion_score: 0.000000
  open_right_insertion_score: -0.300000
  extend_right_insertion_score: 0.000000
  open_internal_deletion_score: -0.200000
  extend_internal_deletion_score: 0.000000
  open_left_deletion_score: -0.200000
  extend_left_deletion_score: 0.000000
  open_right_deletion_score: -0.200000
  extend_right_deletion_score: 0.000000
  mode: local
""",
        )
        seq1 = "GAT"
        seq2 = "GTCT"
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 1.8)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, 1.8)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=1.8) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1.8)
        self.assertEqual(
            str(alignment),
            """\
target            0 GAT 3
                  0 |-| 3
query             0 G-T 2
""",
        )
        self.assertEqual(alignment.shape, (2, 3))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 1], [2, 3]], [[0, 1], [1, 2]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 2 identities; 0 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.8; substitution score = 2.0; gap score = -0.2; 2 aligned letters; 2 identities; 0 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.8:
        substitution_score = 2.0,
        gap_score = -0.2.
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 1.8)
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=1.8) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1.8)
        self.assertEqual(
            str(alignment),
            """\
target            0 GAT 3
                  0 |-| 3
query             4 G-T 2
""",
        )
        self.assertEqual(alignment.shape, (2, 3))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 1], [2, 3]], [[4, 3], [3, 2]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 2 identities; 0 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.8; substitution score = 2.0; gap score = -0.2; 2 aligned letters; 2 identities; 0 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.8:
        substitution_score = 2.0,
        gap_score = -0.2.
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 1.8)


class TestPairwiseSeparateGapPenaltiesWithExtension(unittest.TestCase):
    def test_separate_gap_penalties_with_extension(self):
        seq1 = "GAAT"
        seq2 = "GTCCT"
        aligner = Align.PairwiseAligner()
        aligner.mode = "local"
        open_score, extend_score = (-0.1, 0)
        aligner.open_insertion_score = open_score
        aligner.extend_insertion_score = extend_score
        score = -0.1
        aligner.deletion_score = score
        self.assertEqual(aligner.algorithm, "Gotoh local alignment algorithm")
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  wildcard: None
  match_score: 1.000000
  mismatch_score: 0.000000
  open_internal_insertion_score: -0.100000
  extend_internal_insertion_score: 0.000000
  open_left_insertion_score: -0.100000
  extend_left_insertion_score: 0.000000
  open_right_insertion_score: -0.100000
  extend_right_insertion_score: 0.000000
  open_internal_deletion_score: -0.100000
  extend_internal_deletion_score: -0.100000
  open_left_deletion_score: -0.100000
  extend_left_deletion_score: -0.100000
  open_right_deletion_score: -0.100000
  extend_right_deletion_score: -0.100000
  mode: local
""",
        )
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 1.9)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, 1.9)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (3 alignments; score=1.9) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 3)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1.9)
        self.assertEqual(
            str(alignment),
            """\
target            0 G-AAT 4
                  0 |-..| 5
query             0 GTCCT 5
""",
        )
        self.assertEqual(alignment.shape, (2, 5))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 1], [1, 4]], [[0, 1], [2, 5]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (5 aligned letters; 2 identities; 2 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 4:
        identities = 2,
        mismatches = 2.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 2)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.9; substitution score = 2.0; gap score = -0.1; 4 aligned letters; 2 identities; 2 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.9:
        substitution_score = 2.0,
        gap_score = -0.1.
    aligned = 4:
        identities = 2,
        mismatches = 2.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 2)
        self.assertAlmostEqual(counts.score, 1.9)
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 1.9)
        self.assertEqual(
            str(alignment),
            """\
target            0 GA-AT 4
                  0 |.-.| 5
query             0 GTCCT 5
""",
        )
        self.assertEqual(alignment.shape, (2, 5))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 2], [2, 4]], [[0, 2], [3, 5]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (4 aligned letters; 2 identities; 2 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 4:
        identities = 2,
        mismatches = 2.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 2)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.9; substitution score = 2.0; gap score = -0.1; 4 aligned letters; 2 identities; 2 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.9:
        substitution_score = 2.0,
        gap_score = -0.1.
    aligned = 4:
        identities = 2,
        mismatches = 2.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 2)
        self.assertAlmostEqual(counts.score, 1.9)
        alignment = alignments[2]
        self.assertAlmostEqual(alignment.score, 1.9)
        self.assertEqual(
            str(alignment),
            """\
target            0 GAA-T 4
                  0 |..-| 5
query             0 GTCCT 5
""",
        )
        self.assertEqual(alignment.shape, (2, 5))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 3], [3, 4]], [[0, 3], [4, 5]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (4 aligned letters; 2 identities; 2 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 4:
        identities = 2,
        mismatches = 2.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 2)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.9; substitution score = 2.0; gap score = -0.1; 4 aligned letters; 2 identities; 2 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.9:
        substitution_score = 2.0,
        gap_score = -0.1.
    aligned = 4:
        identities = 2,
        mismatches = 2.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 2)
        self.assertAlmostEqual(counts.score, 1.9)
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (3 alignments; score=1.9) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 3)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1.9)
        self.assertEqual(
            str(alignment),
            """\
target            0 G-AAT 4
                  0 |-..| 5
query             5 GTCCT 0
""",
        )
        self.assertEqual(alignment.shape, (2, 5))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 1], [1, 4]], [[5, 4], [3, 0]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (4 aligned letters; 2 identities; 2 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 4:
        identities = 2,
        mismatches = 2.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 2)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.9; substitution score = 2.0; gap score = -0.1; 4 aligned letters; 2 identities; 2 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.9:
        substitution_score = 2.0,
        gap_score = -0.1.
    aligned = 4:
        identities = 2,
        mismatches = 2.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 2)
        self.assertAlmostEqual(counts.score, 1.9)
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 1.9)
        self.assertEqual(
            str(alignment),
            """\
target            0 GA-AT 4
                  0 |.-.| 5
query             5 GTCCT 0
""",
        )
        self.assertEqual(alignment.shape, (2, 5))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 2], [2, 4]], [[5, 3], [2, 0]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (4 aligned letters; 2 identities; 2 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 4:
        identities = 2,
        mismatches = 2.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 2)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.9; substitution score = 2.0; gap score = -0.1; 4 aligned letters; 2 identities; 2 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.9:
        substitution_score = 2.0,
        gap_score = -0.1.
    aligned = 4:
        identities = 2,
        mismatches = 2.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 2)
        self.assertAlmostEqual(counts.score, 1.9)
        alignment = alignments[2]
        self.assertAlmostEqual(alignment.score, 1.9)
        self.assertEqual(
            str(alignment),
            """\
target            0 GAA-T 4
                  0 |..-| 5
query             5 GTCCT 0
""",
        )
        self.assertEqual(alignment.shape, (2, 5))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 3], [3, 4]], [[5, 2], [1, 0]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (4 aligned letters; 2 identities; 2 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 4:
        identities = 2,
        mismatches = 2.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 2)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.9; substitution score = 2.0; gap score = -0.1; 4 aligned letters; 2 identities; 2 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.9:
        substitution_score = 2.0,
        gap_score = -0.1.
    aligned = 4:
        identities = 2,
        mismatches = 2.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 2)
        self.assertAlmostEqual(counts.score, 1.9)


class TestPairwiseMatchDictionary(unittest.TestCase):
    match_dict = {("A", "A"): 1.5, ("A", "T"): 0.5, ("T", "A"): 0.5, ("T", "T"): 1.0}

    def test_match_dictionary1(self):
        try:
            from Bio.Align import substitution_matrices
        except ImportError:
            return
        substitution_matrix = substitution_matrices.Array(data=self.match_dict)
        seq1 = "ATAT"
        seq2 = "ATT"
        aligner = Align.PairwiseAligner()
        aligner.mode = "local"
        aligner.substitution_matrix = substitution_matrix
        aligner.open_gap_score = -0.5
        aligner.extend_gap_score = 0.0
        self.assertEqual(aligner.algorithm, "Gotoh local alignment algorithm")
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  substitution_matrix: <Array object at 0x%x>
  open_internal_insertion_score: -0.500000
  extend_internal_insertion_score: 0.000000
  open_left_insertion_score: -0.500000
  extend_left_insertion_score: 0.000000
  open_right_insertion_score: -0.500000
  extend_right_insertion_score: 0.000000
  open_internal_deletion_score: -0.500000
  extend_internal_deletion_score: 0.000000
  open_left_deletion_score: -0.500000
  extend_left_deletion_score: 0.000000
  open_right_deletion_score: -0.500000
  extend_right_deletion_score: 0.000000
  mode: local
"""
            % id(substitution_matrix),
        )
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 3.0)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, 3.0)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (2 alignments; score=3) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 ATA 3
                  0 ||. 3
query             0 ATT 3
""",
        )
        self.assertEqual(alignment.shape, (2, 3))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[0, 3]], [[0, 3]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (3 aligned letters; 2 identities; 1 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 3:
        identities = 2,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        counts = alignment.counts(substitution_matrix)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (substitution score = 3.0; 3 aligned letters; 2 identities; 1 mismatches; 3 positives; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    substitution_score = 3.0,
    aligned = 3:
        identities = 2,
        positives = 3,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 3.0; substitution score = 3.0; gap score = 0.0; 3 aligned letters; 2 identities; 1 mismatches; 3 positives; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 3.0:
        substitution_score = 3.0,
        gap_score = 0.0.
    aligned = 3:
        identities = 2,
        positives = 3,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 3.0)
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 ATAT 4
                  0 ||-| 4
query             0 AT-T 3
""",
        )
        self.assertEqual(alignment.shape, (2, 4))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 2], [3, 4]], [[0, 2], [2, 3]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (3 aligned letters; 3 identities; 0 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 3:
        identities = 3,
        mismatches = 0.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(substitution_matrix)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (substitution score = 3.5; 3 aligned letters; 3 identities; 0 mismatches; 3 positives; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    substitution_score = 3.5,
    aligned = 3:
        identities = 3,
        positives = 3,
        mismatches = 0.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 3.0; substitution score = 3.5; gap score = -0.5; 3 aligned letters; 3 identities; 0 mismatches; 3 positives; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 3.0:
        substitution_score = 3.5,
        gap_score = -0.5.
    aligned = 3:
        identities = 3,
        positives = 3,
        mismatches = 0.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 3.0)
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (2 alignments; score=3) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 ATA 3
                  0 ||. 3
query             3 ATT 0
""",
        )
        self.assertEqual(alignment.shape, (2, 3))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[0, 3]], [[3, 0]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (3 aligned letters; 2 identities; 1 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 3:
        identities = 2,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts(substitution_matrix)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (substitution score = 3.0; 3 aligned letters; 2 identities; 1 mismatches; 3 positives; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    substitution_score = 3.0,
    aligned = 3:
        identities = 2,
        positives = 3,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 3.0; substitution score = 3.0; gap score = 0.0; 3 aligned letters; 2 identities; 1 mismatches; 3 positives; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 3.0:
        substitution_score = 3.0,
        gap_score = 0.0.
    aligned = 3:
        identities = 2,
        positives = 3,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 3.0)
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 ATAT 4
                  0 ||-| 4
query             3 AT-T 0
""",
        )
        self.assertEqual(alignment.shape, (2, 4))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 2], [3, 4]], [[3, 1], [1, 0]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (3 aligned letters; 3 identities; 0 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 3:
        identities = 3,
        mismatches = 0.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(substitution_matrix)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (substitution score = 3.5; 3 aligned letters; 3 identities; 0 mismatches; 3 positives; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    substitution_score = 3.5,
    aligned = 3:
        identities = 3,
        positives = 3,
        mismatches = 0.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 3.0; substitution score = 3.5; gap score = -0.5; 3 aligned letters; 3 identities; 0 mismatches; 3 positives; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 3.0:
        substitution_score = 3.5,
        gap_score = -0.5.
    aligned = 3:
        identities = 3,
        positives = 3,
        mismatches = 0.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 3.0)
        self.assertEqual(alignment.sequences[0], "ATAT")
        alignment.sequences[0] = "ATCG"
        with self.assertRaises(ValueError) as cm:
            alignment.counts(aligner)
        self.assertEqual(
            str(cm.exception), "sequence contains letters not in the alphabet"
        )

    def test_match_dictionary2(self):
        try:
            from Bio.Align import substitution_matrices
        except ImportError:
            return
        substitution_matrix = substitution_matrices.Array(data=self.match_dict)
        seq1 = "ATAT"
        seq2 = "ATT"
        aligner = Align.PairwiseAligner()
        aligner.mode = "local"
        aligner.substitution_matrix = substitution_matrix
        aligner.open_gap_score = -1.0
        aligner.extend_gap_score = 0.0
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  substitution_matrix: <Array object at 0x%x>
  open_internal_insertion_score: -1.000000
  extend_internal_insertion_score: 0.000000
  open_left_insertion_score: -1.000000
  extend_left_insertion_score: 0.000000
  open_right_insertion_score: -1.000000
  extend_right_insertion_score: 0.000000
  open_internal_deletion_score: -1.000000
  extend_internal_deletion_score: 0.000000
  open_left_deletion_score: -1.000000
  extend_left_deletion_score: 0.000000
  open_right_deletion_score: -1.000000
  extend_right_deletion_score: 0.000000
  mode: local
"""
            % id(substitution_matrix),
        )
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 3.0)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, 3.0)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=3) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 ATA 3
                  0 ||. 3
query             0 ATT 3
""",
        )
        self.assertEqual(alignment.shape, (2, 3))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[0, 3]], [[0, 3]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (3 aligned letters; 2 identities; 1 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 3:
        identities = 2,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts(substitution_matrix)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (substitution score = 3.0; 3 aligned letters; 2 identities; 1 mismatches; 3 positives; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    substitution_score = 3.0,
    aligned = 3:
        identities = 2,
        positives = 3,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 3.0; substitution score = 3.0; gap score = 0.0; 3 aligned letters; 2 identities; 1 mismatches; 3 positives; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 3.0:
        substitution_score = 3.0,
        gap_score = 0.0.
    aligned = 3:
        identities = 2,
        positives = 3,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 3.0)
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=3) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 ATA 3
                  0 ||. 3
query             3 ATT 0
""",
        )
        self.assertEqual(alignment.shape, (2, 3))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[0, 3]], [[3, 0]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (3 aligned letters; 2 identities; 1 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 3:
        identities = 2,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts(substitution_matrix)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (substitution score = 3.0; 3 aligned letters; 2 identities; 1 mismatches; 3 positives; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    substitution_score = 3.0,
    aligned = 3:
        identities = 2,
        positives = 3,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 3.0; substitution score = 3.0; gap score = 0.0; 3 aligned letters; 2 identities; 1 mismatches; 3 positives; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 3.0:
        substitution_score = 3.0,
        gap_score = 0.0.
    aligned = 3:
        identities = 2,
        positives = 3,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 3.0)

    def test_match_dictionary3(self):
        try:
            from Bio.Align import substitution_matrices
        except ImportError:
            return
        substitution_matrix = substitution_matrices.Array(data=self.match_dict)
        seq1 = "ATT"
        seq2 = "ATAT"
        aligner = Align.PairwiseAligner()
        aligner.mode = "local"
        aligner.substitution_matrix = substitution_matrix
        aligner.open_gap_score = -1.0
        aligner.extend_gap_score = 0.0
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  substitution_matrix: <Array object at 0x%x>
  open_internal_insertion_score: -1.000000
  extend_internal_insertion_score: 0.000000
  open_left_insertion_score: -1.000000
  extend_left_insertion_score: 0.000000
  open_right_insertion_score: -1.000000
  extend_right_insertion_score: 0.000000
  open_internal_deletion_score: -1.000000
  extend_internal_deletion_score: 0.000000
  open_left_deletion_score: -1.000000
  extend_left_deletion_score: 0.000000
  open_right_deletion_score: -1.000000
  extend_right_deletion_score: 0.000000
  mode: local
"""
            % id(substitution_matrix),
        )
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 3.0)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, 3.0)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=3) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 ATT 3
                  0 ||. 3
query             0 ATA 3
""",
        )
        self.assertEqual(alignment.shape, (2, 3))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[0, 3]], [[0, 3]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (3 aligned letters; 2 identities; 1 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 3:
        identities = 2,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts(substitution_matrix)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (substitution score = 3.0; 3 aligned letters; 2 identities; 1 mismatches; 3 positives; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    substitution_score = 3.0,
    aligned = 3:
        identities = 2,
        positives = 3,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 3.0; substitution score = 3.0; gap score = 0.0; 3 aligned letters; 2 identities; 1 mismatches; 3 positives; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 3.0:
        substitution_score = 3.0,
        gap_score = 0.0.
    aligned = 3:
        identities = 2,
        positives = 3,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 3.0)
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=3) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 ATT 3
                  0 ||. 3
query             4 ATA 1
""",
        )
        self.assertEqual(alignment.shape, (2, 3))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[0, 3]], [[4, 1]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (3 aligned letters; 2 identities; 1 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 3:
        identities = 2,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts(substitution_matrix)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (substitution score = 3.0; 3 aligned letters; 2 identities; 1 mismatches; 3 positives; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    substitution_score = 3.0,
    aligned = 3:
        identities = 2,
        positives = 3,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 3.0; substitution score = 3.0; gap score = 0.0; 3 aligned letters; 2 identities; 1 mismatches; 3 positives; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 3.0:
        substitution_score = 3.0,
        gap_score = 0.0.
    aligned = 3:
        identities = 2,
        positives = 3,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 3.0)

    def test_match_dictionary4(self):
        try:
            from Bio.Align import substitution_matrices
        except ImportError:
            return
        substitution_matrix = substitution_matrices.Array(alphabet="AT", dims=2)
        self.assertEqual(substitution_matrix.shape, (2, 2))
        substitution_matrix.update(self.match_dict)
        seq1 = "ATAT"
        seq2 = "ATT"
        aligner = Align.PairwiseAligner()
        aligner.mode = "local"
        aligner.substitution_matrix = substitution_matrix
        aligner.open_gap_score = -0.5
        aligner.extend_gap_score = 0.0
        self.assertEqual(aligner.algorithm, "Gotoh local alignment algorithm")
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  substitution_matrix: <Array object at 0x%x>
  open_internal_insertion_score: -0.500000
  extend_internal_insertion_score: 0.000000
  open_left_insertion_score: -0.500000
  extend_left_insertion_score: 0.000000
  open_right_insertion_score: -0.500000
  extend_right_insertion_score: 0.000000
  open_internal_deletion_score: -0.500000
  extend_internal_deletion_score: 0.000000
  open_left_deletion_score: -0.500000
  extend_left_deletion_score: 0.000000
  open_right_deletion_score: -0.500000
  extend_right_deletion_score: 0.000000
  mode: local
"""
            % id(substitution_matrix),
        )
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 3.0)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, 3.0)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (2 alignments; score=3) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 ATA 3
                  0 ||. 3
query             0 ATT 3
""",
        )
        self.assertEqual(alignment.shape, (2, 3))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[0, 3]], [[0, 3]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (3 aligned letters; 2 identities; 1 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 3:
        identities = 2,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        counts = alignment.counts(substitution_matrix)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (substitution score = 3.0; 3 aligned letters; 2 identities; 1 mismatches; 3 positives; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    substitution_score = 3.0,
    aligned = 3:
        identities = 2,
        positives = 3,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 3.0; substitution score = 3.0; gap score = 0.0; 3 aligned letters; 2 identities; 1 mismatches; 3 positives; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 3.0:
        substitution_score = 3.0,
        gap_score = 0.0.
    aligned = 3:
        identities = 2,
        positives = 3,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 3.0)
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 ATAT 4
                  0 ||-| 4
query             0 AT-T 3
""",
        )
        self.assertEqual(alignment.shape, (2, 4))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 2], [3, 4]], [[0, 2], [2, 3]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (3 aligned letters; 3 identities; 0 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 3:
        identities = 3,
        mismatches = 0.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(substitution_matrix)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (substitution score = 3.5; 3 aligned letters; 3 identities; 0 mismatches; 3 positives; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    substitution_score = 3.5,
    aligned = 3:
        identities = 3,
        positives = 3,
        mismatches = 0.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 3.0; substitution score = 3.5; gap score = -0.5; 3 aligned letters; 3 identities; 0 mismatches; 3 positives; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 3.0:
        substitution_score = 3.5,
        gap_score = -0.5.
    aligned = 3:
        identities = 3,
        positives = 3,
        mismatches = 0.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 3.0)
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (2 alignments; score=3) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 ATA 3
                  0 ||. 3
query             3 ATT 0
""",
        )
        self.assertEqual(alignment.shape, (2, 3))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[0, 3]], [[3, 0]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (3 aligned letters; 2 identities; 1 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 3:
        identities = 2,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts(substitution_matrix)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (substitution score = 3.0; 3 aligned letters; 2 identities; 1 mismatches; 3 positives; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    substitution_score = 3.0,
    aligned = 3:
        identities = 2,
        positives = 3,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 3.0; substitution score = 3.0; gap score = 0.0; 3 aligned letters; 2 identities; 1 mismatches; 3 positives; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 3.0:
        substitution_score = 3.0,
        gap_score = 0.0.
    aligned = 3:
        identities = 2,
        positives = 3,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 3.0)
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 ATAT 4
                  0 ||-| 4
query             3 AT-T 0
""",
        )
        self.assertEqual(alignment.shape, (2, 4))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 2], [3, 4]], [[3, 1], [1, 0]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (3 aligned letters; 3 identities; 0 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 3:
        identities = 3,
        mismatches = 0.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(substitution_matrix)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (substitution score = 3.5; 3 aligned letters; 3 identities; 0 mismatches; 3 positives; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    substitution_score = 3.5,
    aligned = 3:
        identities = 3,
        positives = 3,
        mismatches = 0.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 3.0; substitution score = 3.5; gap score = -0.5; 3 aligned letters; 3 identities; 0 mismatches; 3 positives; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 3.0:
        substitution_score = 3.5,
        gap_score = -0.5.
    aligned = 3:
        identities = 3,
        positives = 3,
        mismatches = 0.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 3.0)

    def test_match_dictionary5(self):
        try:
            from Bio.Align import substitution_matrices
        except ImportError:
            return
        substitution_matrix = substitution_matrices.Array(alphabet="AT", dims=2)
        self.assertEqual(substitution_matrix.shape, (2, 2))
        substitution_matrix.update(self.match_dict)
        seq1 = "ATAT"
        seq2 = "ATT"
        aligner = Align.PairwiseAligner()
        aligner.mode = "local"
        aligner.substitution_matrix = substitution_matrix
        aligner.open_gap_score = -1.0
        aligner.extend_gap_score = 0.0
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  substitution_matrix: <Array object at 0x%x>
  open_internal_insertion_score: -1.000000
  extend_internal_insertion_score: 0.000000
  open_left_insertion_score: -1.000000
  extend_left_insertion_score: 0.000000
  open_right_insertion_score: -1.000000
  extend_right_insertion_score: 0.000000
  open_internal_deletion_score: -1.000000
  extend_internal_deletion_score: 0.000000
  open_left_deletion_score: -1.000000
  extend_left_deletion_score: 0.000000
  open_right_deletion_score: -1.000000
  extend_right_deletion_score: 0.000000
  mode: local
"""
            % id(substitution_matrix),
        )
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 3.0)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, 3.0)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=3) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 ATA 3
                  0 ||. 3
query             0 ATT 3
""",
        )
        self.assertEqual(alignment.shape, (2, 3))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[0, 3]], [[0, 3]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (3 aligned letters; 2 identities; 1 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 3:
        identities = 2,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts(substitution_matrix)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (substitution score = 3.0; 3 aligned letters; 2 identities; 1 mismatches; 3 positives; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    substitution_score = 3.0,
    aligned = 3:
        identities = 2,
        positives = 3,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 3.0; substitution score = 3.0; gap score = 0.0; 3 aligned letters; 2 identities; 1 mismatches; 3 positives; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 3.0:
        substitution_score = 3.0,
        gap_score = 0.0.
    aligned = 3:
        identities = 2,
        positives = 3,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 3.0)
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=3) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 ATA 3
                  0 ||. 3
query             3 ATT 0
""",
        )
        self.assertEqual(alignment.shape, (2, 3))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[0, 3]], [[3, 0]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (3 aligned letters; 2 identities; 1 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 3:
        identities = 2,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts(substitution_matrix)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (substitution score = 3.0; 3 aligned letters; 2 identities; 1 mismatches; 3 positives; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    substitution_score = 3.0,
    aligned = 3:
        identities = 2,
        positives = 3,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 3.0; substitution score = 3.0; gap score = 0.0; 3 aligned letters; 2 identities; 1 mismatches; 3 positives; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 3.0:
        substitution_score = 3.0,
        gap_score = 0.0.
    aligned = 3:
        identities = 2,
        positives = 3,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 3.0)

    def test_match_dictionary6(self):
        try:
            from Bio.Align import substitution_matrices
        except ImportError:
            return
        substitution_matrix = substitution_matrices.Array(alphabet="AT", dims=2)
        self.assertEqual(substitution_matrix.shape, (2, 2))
        substitution_matrix.update(self.match_dict)
        seq1 = "ATT"
        seq2 = "ATAT"
        aligner = Align.PairwiseAligner()
        aligner.mode = "local"
        aligner.substitution_matrix = substitution_matrix
        aligner.open_gap_score = -1.0
        aligner.extend_gap_score = 0.0
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  substitution_matrix: <Array object at 0x%x>
  open_internal_insertion_score: -1.000000
  extend_internal_insertion_score: 0.000000
  open_left_insertion_score: -1.000000
  extend_left_insertion_score: 0.000000
  open_right_insertion_score: -1.000000
  extend_right_insertion_score: 0.000000
  open_internal_deletion_score: -1.000000
  extend_internal_deletion_score: 0.000000
  open_left_deletion_score: -1.000000
  extend_left_deletion_score: 0.000000
  open_right_deletion_score: -1.000000
  extend_right_deletion_score: 0.000000
  mode: local
"""
            % id(substitution_matrix),
        )
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 3.0)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, 3.0)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=3) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 ATT 3
                  0 ||. 3
query             0 ATA 3
""",
        )
        self.assertEqual(alignment.shape, (2, 3))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[0, 3]], [[0, 3]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (3 aligned letters; 2 identities; 1 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 3:
        identities = 2,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts(substitution_matrix)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (substitution score = 3.0; 3 aligned letters; 2 identities; 1 mismatches; 3 positives; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    substitution_score = 3.0,
    aligned = 3:
        identities = 2,
        positives = 3,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 3.0; substitution score = 3.0; gap score = 0.0; 3 aligned letters; 2 identities; 1 mismatches; 3 positives; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 3.0:
        substitution_score = 3.0,
        gap_score = 0.0.
    aligned = 3:
        identities = 2,
        positives = 3,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 3.0)
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=3) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 ATT 3
                  0 ||. 3
query             4 ATA 1
""",
        )
        self.assertEqual(alignment.shape, (2, 3))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[0, 3]], [[4, 1]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (3 aligned letters; 2 identities; 1 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 3:
        identities = 2,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts(substitution_matrix)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (substitution score = 3.0; 3 aligned letters; 2 identities; 1 mismatches; 3 positives; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    substitution_score = 3.0,
    aligned = 3:
        identities = 2,
        positives = 3,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 3.0; substitution score = 3.0; gap score = 0.0; 3 aligned letters; 2 identities; 1 mismatches; 3 positives; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 3.0:
        substitution_score = 3.0,
        gap_score = 0.0.
    aligned = 3:
        identities = 2,
        positives = 3,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 3.0)


class TestPairwiseOneCharacter(unittest.TestCase):
    def test_align_one_char1(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = "local"
        aligner.open_gap_score = -0.3
        aligner.extend_gap_score = -0.1
        self.assertEqual(aligner.algorithm, "Gotoh local alignment algorithm")
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  wildcard: None
  match_score: 1.000000
  mismatch_score: 0.000000
  open_internal_insertion_score: -0.300000
  extend_internal_insertion_score: -0.100000
  open_left_insertion_score: -0.300000
  extend_left_insertion_score: -0.100000
  open_right_insertion_score: -0.300000
  extend_right_insertion_score: -0.100000
  open_internal_deletion_score: -0.300000
  extend_internal_deletion_score: -0.100000
  open_left_deletion_score: -0.300000
  extend_left_deletion_score: -0.100000
  open_right_deletion_score: -0.300000
  extend_right_deletion_score: -0.100000
  mode: local
""",
        )
        score = aligner.score("abcde", "c")
        self.assertAlmostEqual(score, 1)
        alignments = aligner.align("abcde", "c")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=1) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1.0)
        self.assertEqual(
            str(alignment),
            """\
target            2 c 3
                  0 | 1
query             0 c 1
""",
        )
        self.assertEqual(alignment.shape, (2, 1))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[2, 3]], [[0, 1]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (1 aligned letters; 1 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 1:
        identities = 1,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 1)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.0; substitution score = 1.0; gap score = 0.0; 1 aligned letters; 1 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.0:
        substitution_score = 1.0,
        gap_score = 0.0.
    aligned = 1:
        identities = 1,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 1)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 1.0)

    def test_align_one_char2(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = "local"
        aligner.open_gap_score = -0.3
        aligner.extend_gap_score = -0.1
        self.assertEqual(aligner.algorithm, "Gotoh local alignment algorithm")
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  wildcard: None
  match_score: 1.000000
  mismatch_score: 0.000000
  open_internal_insertion_score: -0.300000
  extend_internal_insertion_score: -0.100000
  open_left_insertion_score: -0.300000
  extend_left_insertion_score: -0.100000
  open_right_insertion_score: -0.300000
  extend_right_insertion_score: -0.100000
  open_internal_deletion_score: -0.300000
  extend_internal_deletion_score: -0.100000
  open_left_deletion_score: -0.300000
  extend_left_deletion_score: -0.100000
  open_right_deletion_score: -0.300000
  extend_right_deletion_score: -0.100000
  mode: local
""",
        )
        score = aligner.score("abcce", "c")
        self.assertAlmostEqual(score, 1)
        alignments = aligner.align("abcce", "c")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (2 alignments; score=1) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1.0)
        self.assertEqual(
            str(alignment),
            """\
target            2 c 3
                  0 | 1
query             0 c 1
""",
        )
        self.assertEqual(alignment.shape, (2, 1))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[2, 3]], [[0, 1]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (1 aligned letters; 1 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 1:
        identities = 1,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 1)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.0; substitution score = 1.0; gap score = 0.0; 1 aligned letters; 1 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.0:
        substitution_score = 1.0,
        gap_score = 0.0.
    aligned = 1:
        identities = 1,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 1)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 1.0)
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 1.0)
        self.assertEqual(
            str(alignment),
            """\
target            3 c 4
                  0 | 1
query             0 c 1
""",
        )
        self.assertEqual(alignment.shape, (2, 1))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[3, 4]], [[0, 1]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (1 aligned letters; 1 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 1:
        identities = 1,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 1)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.0; substitution score = 1.0; gap score = 0.0; 1 aligned letters; 1 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.0:
        substitution_score = 1.0,
        gap_score = 0.0.
    aligned = 1:
        identities = 1,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 1)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 1.0)

    def test_align_one_char3(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        aligner.open_gap_score = -0.3
        aligner.extend_gap_score = -0.1
        self.assertEqual(aligner.algorithm, "Gotoh global alignment algorithm")
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  wildcard: None
  match_score: 1.000000
  mismatch_score: 0.000000
  open_internal_insertion_score: -0.300000
  extend_internal_insertion_score: -0.100000
  open_left_insertion_score: -0.300000
  extend_left_insertion_score: -0.100000
  open_right_insertion_score: -0.300000
  extend_right_insertion_score: -0.100000
  open_internal_deletion_score: -0.300000
  extend_internal_deletion_score: -0.100000
  open_left_deletion_score: -0.300000
  extend_left_deletion_score: -0.100000
  open_right_deletion_score: -0.300000
  extend_right_deletion_score: -0.100000
  mode: global
""",
        )
        seq1 = "abcde"
        seq2 = "c"
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 0.2)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=0.2) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 0.2)
        self.assertEqual(
            str(alignment),
            """\
target            0 abcde 5
                  0 --|-- 5
query             0 --c-- 1
""",
        )
        self.assertEqual(alignment.shape, (2, 5))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[2, 3]], [[0, 1]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (1 aligned letters; 1 identities; 0 mismatches; 4 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 1:
        identities = 1,
        mismatches = 0.
    gaps = 4:
        left_gaps = 2:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 2:
                open_left_deletions = 1,
                extend_left_deletions = 1;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 2:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 2:
                open_right_deletions = 1,
                extend_right_deletions = 1.
""",
        )
        self.assertEqual(counts.aligned, 1)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 0.20000000000000007; substitution score = 1.0; gap score = -0.7999999999999999; 1 aligned letters; 1 identities; 0 mismatches; 4 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 0.20000000000000007:
        substitution_score = 1.0,
        gap_score = -0.7999999999999999.
    aligned = 1:
        identities = 1,
        mismatches = 0.
    gaps = 4:
        left_gaps = 2:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 2:
                open_left_deletions = 1,
                extend_left_deletions = 1;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 2:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 2:
                open_right_deletions = 1,
                extend_right_deletions = 1.
""",
        )
        self.assertEqual(counts.aligned, 1)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 0.2)

    def test_align_one_char_score3(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        aligner.open_gap_score = -0.3
        aligner.extend_gap_score = -0.1
        self.assertEqual(aligner.algorithm, "Gotoh global alignment algorithm")
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  wildcard: None
  match_score: 1.000000
  mismatch_score: 0.000000
  open_internal_insertion_score: -0.300000
  extend_internal_insertion_score: -0.100000
  open_left_insertion_score: -0.300000
  extend_left_insertion_score: -0.100000
  open_right_insertion_score: -0.300000
  extend_right_insertion_score: -0.100000
  open_internal_deletion_score: -0.300000
  extend_internal_deletion_score: -0.100000
  open_left_deletion_score: -0.300000
  extend_left_deletion_score: -0.100000
  open_right_deletion_score: -0.300000
  extend_right_deletion_score: -0.100000
  mode: global
""",
        )
        score = aligner.score("abcde", "c")
        self.assertAlmostEqual(score, 0.2)

    def test_align_one_char_score3_fogsaa(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = "fogsaa"
        aligner.open_gap_score = -0.3
        aligner.extend_gap_score = -0.1
        self.assertEqual(
            aligner.algorithm, "Fast Optimal Global Sequence Alignment Algorithm"
        )
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  wildcard: None
  match_score: 1.000000
  mismatch_score: 0.000000
  open_internal_insertion_score: -0.300000
  extend_internal_insertion_score: -0.100000
  open_left_insertion_score: -0.300000
  extend_left_insertion_score: -0.100000
  open_right_insertion_score: -0.300000
  extend_right_insertion_score: -0.100000
  open_internal_deletion_score: -0.300000
  extend_internal_deletion_score: -0.100000
  open_left_deletion_score: -0.300000
  extend_left_deletion_score: -0.100000
  open_right_deletion_score: -0.300000
  extend_right_deletion_score: -0.100000
  mode: fogsaa
""",
        )
        score = aligner.score("abcde", "c")
        self.assertAlmostEqual(score, 0.2)


class TestPerSiteGapPenalties(unittest.TestCase):
    """Check gap penalty callbacks use correct gap opening position.

    This tests that the gap penalty callbacks are really being used
    with the correct gap opening position.
    """

    def test_gap_here_only_1(self):
        seq1 = "AAAABBBAAAACCCCCCCCCCCCCCAAAABBBAAAA"
        seq2 = "AABBBAAAACCCCAAAABBBAA"
        breaks = [0, 11, len(seq2)]

        # Very expensive to open a gap in seq1:
        def nogaps(x, y):
            return -2000 - y

        # Very expensive to open a gap in seq2 unless it is in one of the allowed positions:
        def specificgaps(x, y):
            if x in breaks:
                return -2 - y
            else:
                return -2000 - y

        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score = 1
        aligner.mismatch_score = -1
        aligner.insertion_score = nogaps
        aligner.deletion_score = specificgaps
        self.assertEqual(
            str(aligner),
            f"""Pairwise sequence aligner with parameters
  wildcard: None
  match_score: 1.000000
  mismatch_score: -1.000000
  insertion_score_function: {nogaps}
  deletion_score_function: {specificgaps}
  mode: global
""",
        )
        self.assertEqual(
            aligner.algorithm, "Waterman-Smith-Beyer global alignment algorithm"
        )
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 2)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, 2)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=2) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 2.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 AAAABBBAAAACCCCCCCCCCCCCCAAAABBBAAAA 36
                  0 --|||||||||||----------|||||||||||-- 36
query             0 --AABBBAAAACC----------CCAAAABBBAA-- 22
""",
        )
        self.assertEqual(alignment.shape, (2, 36))
        self.assertTrue(
            np.array_equal(
                alignment.aligned,
                np.array([[[2, 13], [23, 34]], [[0, 11], [11, 22]]]),
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (22 aligned letters; 22 identities; 0 mismatches; 14 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 22:
        identities = 22,
        mismatches = 0.
    gaps = 14:
        left_gaps = 2:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 2:
                open_left_deletions = 1,
                extend_left_deletions = 1;
        internal_gaps = 10:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 10:
                open_internal_deletions = 1,
                extend_internal_deletions = 9;
        right_gaps = 2:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 2:
                open_right_deletions = 1,
                extend_right_deletions = 1.
""",
        )
        self.assertEqual(counts.aligned, 22)
        self.assertEqual(counts.identities, 22)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 2.0; substitution score = 22.0; gap score = -20.0; 22 aligned letters; 22 identities; 0 mismatches; 14 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 2.0:
        substitution_score = 22.0,
        gap_score = -20.0.
    aligned = 22:
        identities = 22,
        mismatches = 0.
    gaps = 14:
        left_gaps = 2:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 2:
                open_left_deletions = 1,
                extend_left_deletions = 1;
        internal_gaps = 10:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 10:
                open_internal_deletions = 1,
                extend_internal_deletions = 9;
        right_gaps = 2:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 2:
                open_right_deletions = 1,
                extend_right_deletions = 1.
""",
        )
        self.assertEqual(counts.aligned, 22)
        self.assertEqual(counts.identities, 22)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 2.0)
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=2) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 2.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 AAAABBBAAAACCCCCCCCCCCCCCAAAABBBAAAA 36
                  0 --|||||||||||----------|||||||||||-- 36
query            22 --AABBBAAAACC----------CCAAAABBBAA--  0
""",
        )
        self.assertEqual(alignment.shape, (2, 36))
        self.assertTrue(
            np.array_equal(
                alignment.aligned,
                np.array([[[2, 13], [23, 34]], [[22, 11], [11, 0]]]),
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (22 aligned letters; 22 identities; 0 mismatches; 14 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 22:
        identities = 22,
        mismatches = 0.
    gaps = 14:
        left_gaps = 2:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 2:
                open_left_deletions = 1,
                extend_left_deletions = 1;
        internal_gaps = 10:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 10:
                open_internal_deletions = 1,
                extend_internal_deletions = 9;
        right_gaps = 2:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 2:
                open_right_deletions = 1,
                extend_right_deletions = 1.
""",
        )
        self.assertEqual(counts.aligned, 22)
        self.assertEqual(counts.identities, 22)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 2.0; substitution score = 22.0; gap score = -20.0; 22 aligned letters; 22 identities; 0 mismatches; 14 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 2.0:
        substitution_score = 22.0,
        gap_score = -20.0.
    aligned = 22:
        identities = 22,
        mismatches = 0.
    gaps = 14:
        left_gaps = 2:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 2:
                open_left_deletions = 1,
                extend_left_deletions = 1;
        internal_gaps = 10:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 10:
                open_internal_deletions = 1,
                extend_internal_deletions = 9;
        right_gaps = 2:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 2:
                open_right_deletions = 1,
                extend_right_deletions = 1.
""",
        )
        self.assertEqual(counts.aligned, 22)
        self.assertEqual(counts.identities, 22)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 2.0)

    def test_gap_here_only_2(self):
        # Force a bad alignment.
        #
        # Forces a bad alignment by having a very expensive gap penalty
        # where one would normally expect a gap, and a cheap gap penalty
        # in another place.
        seq1 = "AAAABBBAAAACCCCCCCCCCCCCCAAAABBBAAAA"
        seq2 = "AABBBAAAACCCCAAAABBBAA"
        breaks = [0, 3, len(seq2)]

        # Very expensive to open a gap in seq1:
        def nogaps(x, y):
            return -2000 - y

        # Very expensive to open a gap in seq2 unless it is in one of the allowed positions:
        def specificgaps(x, y):
            if x in breaks:
                return -2 - y
            else:
                return -2000 - y

        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score = 1
        aligner.mismatch_score = -1
        aligner.insertion_score = nogaps
        aligner.deletion_score = specificgaps
        self.assertEqual(
            str(aligner),
            f"""Pairwise sequence aligner with parameters
  wildcard: None
  match_score: 1.000000
  mismatch_score: -1.000000
  insertion_score_function: {nogaps}
  deletion_score_function: {specificgaps}
  mode: global
""",
        )
        self.assertEqual(
            aligner.algorithm, "Waterman-Smith-Beyer global alignment algorithm"
        )
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, -10)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, -10)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (2 alignments; score=-10) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, -10.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 AAAABBBAAAACCCCCCCCCCCCCCAAAABBBAAAA 36
                  0 --|||----------......|||||||||||||-- 36
query             0 --AAB----------BBAAAACCCCAAAABBBAA-- 22
""",
        )
        self.assertEqual(alignment.shape, (2, 36))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[2, 5], [15, 34]], [[0, 3], [3, 22]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (22 aligned letters; 16 identities; 6 mismatches; 14 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 22:
        identities = 16,
        mismatches = 6.
    gaps = 14:
        left_gaps = 2:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 2:
                open_left_deletions = 1,
                extend_left_deletions = 1;
        internal_gaps = 10:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 10:
                open_internal_deletions = 1,
                extend_internal_deletions = 9;
        right_gaps = 2:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 2:
                open_right_deletions = 1,
                extend_right_deletions = 1.
""",
        )
        self.assertEqual(counts.aligned, 22)
        self.assertEqual(counts.identities, 16)
        self.assertEqual(counts.mismatches, 6)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = -10.0; substitution score = 10.0; gap score = -20.0; 22 aligned letters; 16 identities; 6 mismatches; 14 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = -10.0:
        substitution_score = 10.0,
        gap_score = -20.0.
    aligned = 22:
        identities = 16,
        mismatches = 6.
    gaps = 14:
        left_gaps = 2:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 2:
                open_left_deletions = 1,
                extend_left_deletions = 1;
        internal_gaps = 10:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 10:
                open_internal_deletions = 1,
                extend_internal_deletions = 9;
        right_gaps = 2:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 2:
                open_right_deletions = 1,
                extend_right_deletions = 1.
""",
        )
        self.assertEqual(counts.aligned, 22)
        self.assertEqual(counts.identities, 16)
        self.assertEqual(counts.mismatches, 6)
        self.assertAlmostEqual(counts.score, -10.0)
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, -10.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 AAAABBBAAAACCCCCCCCCCCCCCAAAABBBAAAA 36
                  0 ||.------------......|||||||||||||-- 36
query             0 AAB------------BBAAAACCCCAAAABBBAA-- 22
""",
        )
        self.assertEqual(alignment.shape, (2, 36))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 3], [15, 34]], [[0, 3], [3, 22]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (22 aligned letters; 15 identities; 7 mismatches; 14 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 22:
        identities = 15,
        mismatches = 7.
    gaps = 14:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 12:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 12:
                open_internal_deletions = 1,
                extend_internal_deletions = 11;
        right_gaps = 2:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 2:
                open_right_deletions = 1,
                extend_right_deletions = 1.
""",
        )
        self.assertEqual(counts.aligned, 22)
        self.assertEqual(counts.identities, 15)
        self.assertEqual(counts.mismatches, 7)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = -10.0; substitution score = 8.0; gap score = -18.0; 22 aligned letters; 15 identities; 7 mismatches; 14 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = -10.0:
        substitution_score = 8.0,
        gap_score = -18.0.
    aligned = 22:
        identities = 15,
        mismatches = 7.
    gaps = 14:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 12:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 12:
                open_internal_deletions = 1,
                extend_internal_deletions = 11;
        right_gaps = 2:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 2:
                open_right_deletions = 1,
                extend_right_deletions = 1.
""",
        )
        self.assertEqual(counts.aligned, 22)
        self.assertEqual(counts.identities, 15)
        self.assertEqual(counts.mismatches, 7)
        self.assertAlmostEqual(counts.score, -10.0)
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (2 alignments; score=-10) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, -10.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 AAAABBBAAAACCCCCCCCCCCCCCAAAABBBAAAA 36
                  0 --|||||||||||||......------------.|| 36
query            22 --AABBBAAAACCCCAAAABB------------BAA  0
""",
        )
        self.assertEqual(alignment.shape, (2, 36))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[2, 21], [33, 36]], [[22, 3], [3, 0]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (22 aligned letters; 15 identities; 7 mismatches; 14 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 22:
        identities = 15,
        mismatches = 7.
    gaps = 14:
        left_gaps = 2:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 2:
                open_left_deletions = 1,
                extend_left_deletions = 1;
        internal_gaps = 12:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 12:
                open_internal_deletions = 1,
                extend_internal_deletions = 11;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 22)
        self.assertEqual(counts.identities, 15)
        self.assertEqual(counts.mismatches, 7)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = -10.0; substitution score = 8.0; gap score = -18.0; 22 aligned letters; 15 identities; 7 mismatches; 14 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = -10.0:
        substitution_score = 8.0,
        gap_score = -18.0.
    aligned = 22:
        identities = 15,
        mismatches = 7.
    gaps = 14:
        left_gaps = 2:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 2:
                open_left_deletions = 1,
                extend_left_deletions = 1;
        internal_gaps = 12:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 12:
                open_internal_deletions = 1,
                extend_internal_deletions = 11;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 22)
        self.assertEqual(counts.identities, 15)
        self.assertEqual(counts.mismatches, 7)
        self.assertAlmostEqual(counts.score, -10.0)
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, -10.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 AAAABBBAAAACCCCCCCCCCCCCCAAAABBBAAAA 36
                  0 --|||||||||||||......----------|||-- 36
query            22 --AABBBAAAACCCCAAAABB----------BAA--  0
""",
        )
        self.assertEqual(alignment.shape, (2, 36))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[2, 21], [31, 34]], [[22, 3], [3, 0]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (22 aligned letters; 16 identities; 6 mismatches; 14 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 22:
        identities = 16,
        mismatches = 6.
    gaps = 14:
        left_gaps = 2:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 2:
                open_left_deletions = 1,
                extend_left_deletions = 1;
        internal_gaps = 10:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 10:
                open_internal_deletions = 1,
                extend_internal_deletions = 9;
        right_gaps = 2:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 2:
                open_right_deletions = 1,
                extend_right_deletions = 1.
""",
        )
        self.assertEqual(counts.aligned, 22)
        self.assertEqual(counts.identities, 16)
        self.assertEqual(counts.mismatches, 6)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = -10.0; substitution score = 10.0; gap score = -20.0; 22 aligned letters; 16 identities; 6 mismatches; 14 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = -10.0:
        substitution_score = 10.0,
        gap_score = -20.0.
    aligned = 22:
        identities = 16,
        mismatches = 6.
    gaps = 14:
        left_gaps = 2:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 2:
                open_left_deletions = 1,
                extend_left_deletions = 1;
        internal_gaps = 10:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 10:
                open_internal_deletions = 1,
                extend_internal_deletions = 9;
        right_gaps = 2:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 2:
                open_right_deletions = 1,
                extend_right_deletions = 1.
""",
        )
        self.assertEqual(counts.aligned, 22)
        self.assertEqual(counts.identities, 16)
        self.assertEqual(counts.mismatches, 6)
        self.assertAlmostEqual(counts.score, -10.0)

    def test_gap_here_only_3(self):
        # Check if gap open and gap extend penalties are handled correctly.
        seq1 = "TTCCAA"
        seq2 = "TTGGAA"

        def gap_score(i, n):
            if i == 3:
                return -10
            if n == 1:
                return -1
            return -10

        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score = 1
        aligner.mismatch_score = -10
        aligner.insertion_score = gap_score
        self.assertEqual(
            aligner.algorithm, "Waterman-Smith-Beyer global alignment algorithm"
        )
        self.assertEqual(
            str(aligner),
            f"""Pairwise sequence aligner with parameters
  wildcard: None
  match_score: 1.000000
  mismatch_score: -10.000000
  insertion_score_function: {gap_score}
  open_internal_deletion_score: 0.000000
  extend_internal_deletion_score: 0.000000
  open_left_deletion_score: 0.000000
  extend_left_deletion_score: 0.000000
  open_right_deletion_score: 0.000000
  extend_right_deletion_score: 0.000000
  mode: global
""",
        )
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 2.0)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, 2.0)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=2) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 2.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 TT-CC-AA 6
                  0 ||----|| 8
query             0 TTG--GAA 6
""",
        )
        self.assertEqual(alignment.shape, (2, 8))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 2], [4, 6]], [[0, 2], [4, 6]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (4 aligned letters; 4 identities; 0 mismatches; 4 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 4:
        identities = 4,
        mismatches = 0.
    gaps = 4:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 4:
            internal_insertions = 2:
                open_internal_insertions = 2,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 4)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 2.0; substitution score = 4.0; gap score = -2.0; 4 aligned letters; 4 identities; 0 mismatches; 4 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 2.0:
        substitution_score = 4.0,
        gap_score = -2.0.
    aligned = 4:
        identities = 4,
        mismatches = 0.
    gaps = 4:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 4:
            internal_insertions = 2:
                open_internal_insertions = 2,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 4)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 2.0)
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=2) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 2.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 TT-CC-AA 6
                  0 ||----|| 8
query             6 TTG--GAA 0
""",
        )
        self.assertEqual(alignment.shape, (2, 8))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 2], [4, 6]], [[6, 4], [2, 0]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (4 aligned letters; 4 identities; 0 mismatches; 4 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 4:
        identities = 4,
        mismatches = 0.
    gaps = 4:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 4:
            internal_insertions = 2:
                open_internal_insertions = 2,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 4)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 2.0; substitution score = 4.0; gap score = -2.0; 4 aligned letters; 4 identities; 0 mismatches; 4 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 2.0:
        substitution_score = 4.0,
        gap_score = -2.0.
    aligned = 4:
        identities = 4,
        mismatches = 0.
    gaps = 4:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 4:
            internal_insertions = 2:
                open_internal_insertions = 2,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 4)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 2.0)
        aligner.deletion_score = gap_score
        self.assertEqual(
            str(aligner),
            f"""Pairwise sequence aligner with parameters
  wildcard: None
  match_score: 1.000000
  mismatch_score: -10.000000
  insertion_score_function: {gap_score}
  deletion_score_function: {gap_score}
  mode: global
""",
        )
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, -8.0)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, -8.0)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (4 alignments; score=-8) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 4)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, -8.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 TT-CCAA 6
                  0 ||-.-|| 7
query             0 TTGG-AA 6
""",
        )
        self.assertEqual(alignment.shape, (2, 7))
        self.assertTrue(
            np.array_equal(
                alignment.aligned,
                np.array([[[0, 2], [2, 3], [4, 6]], [[0, 2], [3, 4], [4, 6]]]),
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (5 aligned letters; 4 identities; 1 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 5:
        identities = 4,
        mismatches = 1.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 4)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = -8.0; substitution score = -6.0; gap score = -2.0; 5 aligned letters; 4 identities; 1 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = -8.0:
        substitution_score = -6.0,
        gap_score = -2.0.
    aligned = 5:
        identities = 4,
        mismatches = 1.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 4)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, -8.0)
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, -8.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 TTC--CAA 6
                  0 ||----|| 8
query             0 TT-GG-AA 6
""",
        )
        self.assertEqual(alignment.shape, (2, 8))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 2], [4, 6]], [[0, 2], [4, 6]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (4 aligned letters; 4 identities; 0 mismatches; 4 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 4:
        identities = 4,
        mismatches = 0.
    gaps = 4:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 4:
            internal_insertions = 2:
                open_internal_insertions = 1,
                extend_internal_insertions = 1;
            internal_deletions = 2:
                open_internal_deletions = 2,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 4)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = -8.0; substitution score = 4.0; gap score = -12.0; 4 aligned letters; 4 identities; 0 mismatches; 4 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = -8.0:
        substitution_score = 4.0,
        gap_score = -12.0.
    aligned = 4:
        identities = 4,
        mismatches = 0.
    gaps = 4:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 4:
            internal_insertions = 2:
                open_internal_insertions = 1,
                extend_internal_insertions = 1;
            internal_deletions = 2:
                open_internal_deletions = 2,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 4)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, -8.0)
        alignment = alignments[2]
        self.assertAlmostEqual(alignment.score, -8.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 TTCC-AA 6
                  0 ||-.-|| 7
query             0 TT-GGAA 6
""",
        )
        self.assertEqual(alignment.shape, (2, 7))
        self.assertTrue(
            np.array_equal(
                alignment.aligned,
                np.array([[[0, 2], [3, 4], [4, 6]], [[0, 2], [2, 3], [4, 6]]]),
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (5 aligned letters; 4 identities; 1 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 5:
        identities = 4,
        mismatches = 1.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 4)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = -8.0; substitution score = -6.0; gap score = -2.0; 5 aligned letters; 4 identities; 1 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = -8.0:
        substitution_score = -6.0,
        gap_score = -2.0.
    aligned = 5:
        identities = 4,
        mismatches = 1.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 4)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, -8.0)
        alignment = alignments[3]
        self.assertAlmostEqual(alignment.score, -8.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 TT-CC-AA 6
                  0 ||----|| 8
query             0 TTG--GAA 6
""",
        )
        self.assertEqual(alignment.shape, (2, 8))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 2], [4, 6]], [[0, 2], [4, 6]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (4 aligned letters; 4 identities; 0 mismatches; 4 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 4:
        identities = 4,
        mismatches = 0.
    gaps = 4:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 4:
            internal_insertions = 2:
                open_internal_insertions = 2,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 4)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = -8.0; substitution score = 4.0; gap score = -12.0; 4 aligned letters; 4 identities; 0 mismatches; 4 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = -8.0:
        substitution_score = 4.0,
        gap_score = -12.0.
    aligned = 4:
        identities = 4,
        mismatches = 0.
    gaps = 4:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 4:
            internal_insertions = 2:
                open_internal_insertions = 2,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 4)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, -8.0)
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (4 alignments; score=-8) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 4)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, -8.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 TT-CCAA 6
                  0 ||-.-|| 7
query             6 TTGG-AA 0
""",
        )
        self.assertEqual(alignment.shape, (2, 7))
        self.assertTrue(
            np.array_equal(
                alignment.aligned,
                np.array([[[0, 2], [2, 3], [4, 6]], [[6, 4], [3, 2], [2, 0]]]),
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (5 aligned letters; 4 identities; 1 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 5:
        identities = 4,
        mismatches = 1.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 4)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = -8.0; substitution score = -6.0; gap score = -2.0; 5 aligned letters; 4 identities; 1 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = -8.0:
        substitution_score = -6.0,
        gap_score = -2.0.
    aligned = 5:
        identities = 4,
        mismatches = 1.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 4)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, -8.0)
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, -8.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 TTC--CAA 6
                  0 ||----|| 8
query             6 TT-GG-AA 0
""",
        )
        self.assertEqual(alignment.shape, (2, 8))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 2], [4, 6]], [[6, 4], [2, 0]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (4 aligned letters; 4 identities; 0 mismatches; 4 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 4:
        identities = 4,
        mismatches = 0.
    gaps = 4:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 4:
            internal_insertions = 2:
                open_internal_insertions = 1,
                extend_internal_insertions = 1;
            internal_deletions = 2:
                open_internal_deletions = 2,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 4)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = -8.0; substitution score = 4.0; gap score = -12.0; 4 aligned letters; 4 identities; 0 mismatches; 4 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = -8.0:
        substitution_score = 4.0,
        gap_score = -12.0.
    aligned = 4:
        identities = 4,
        mismatches = 0.
    gaps = 4:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 4:
            internal_insertions = 2:
                open_internal_insertions = 1,
                extend_internal_insertions = 1;
            internal_deletions = 2:
                open_internal_deletions = 2,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 4)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, -8.0)
        alignment = alignments[2]
        self.assertAlmostEqual(alignment.score, -8.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 TTCC-AA 6
                  0 ||-.-|| 7
query             6 TT-GGAA 0
""",
        )
        self.assertEqual(alignment.shape, (2, 7))
        self.assertTrue(
            np.array_equal(
                alignment.aligned,
                np.array([[[0, 2], [3, 4], [4, 6]], [[6, 4], [4, 3], [2, 0]]]),
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (5 aligned letters; 4 identities; 1 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 5:
        identities = 4,
        mismatches = 1.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 4)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = -8.0; substitution score = -6.0; gap score = -2.0; 5 aligned letters; 4 identities; 1 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = -8.0:
        substitution_score = -6.0,
        gap_score = -2.0.
    aligned = 5:
        identities = 4,
        mismatches = 1.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 4)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, -8.0)
        alignment = alignments[3]
        self.assertAlmostEqual(alignment.score, -8.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 TT-CC-AA 6
                  0 ||----|| 8
query             6 TTG--GAA 0
""",
        )
        self.assertEqual(alignment.shape, (2, 8))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 2], [4, 6]], [[6, 4], [2, 0]]])
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (4 aligned letters; 4 identities; 0 mismatches; 4 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 4:
        identities = 4,
        mismatches = 0.
    gaps = 4:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 4:
            internal_insertions = 2:
                open_internal_insertions = 2,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 4)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = -8.0; substitution score = 4.0; gap score = -12.0; 4 aligned letters; 4 identities; 0 mismatches; 4 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = -8.0:
        substitution_score = 4.0,
        gap_score = -12.0.
    aligned = 4:
        identities = 4,
        mismatches = 0.
    gaps = 4:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 4:
            internal_insertions = 2:
                open_internal_insertions = 2,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 4)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, -8.0)

    def test_gap_here_only_local_1(self):
        seq1 = "AAAABBBAAAACCCCCCCCCCCCCCAAAABBBAAAA"
        seq2 = "AABBBAAAACCCCAAAABBBAA"
        breaks = [0, 11, len(seq2)]

        # Very expensive to open a gap in seq1:
        def nogaps(x, y):
            return -2000 - y

        # Very expensive to open a gap in seq2 unless it is in one of the allowed positions:
        def specificgaps(x, y):
            if x in breaks:
                return -2 - y
            else:
                return -2000 - y

        aligner = Align.PairwiseAligner()
        aligner.mode = "local"
        aligner.match_score = 1
        aligner.mismatch_score = -1
        aligner.insertion_score = nogaps
        aligner.deletion_score = specificgaps
        self.assertEqual(
            aligner.algorithm, "Waterman-Smith-Beyer local alignment algorithm"
        )
        self.assertEqual(
            str(aligner),
            f"""Pairwise sequence aligner with parameters
  wildcard: None
  match_score: 1.000000
  mismatch_score: -1.000000
  insertion_score_function: {nogaps}
  deletion_score_function: {specificgaps}
  mode: local
""",
        )
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 13)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, 13)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (2 alignments; score=13) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 13)
        self.assertEqual(
            str(alignment),
            """\
target            2 AABBBAAAACCCC 15
                  0 ||||||||||||| 13
query             0 AABBBAAAACCCC 13
""",
        )
        self.assertEqual(alignment.shape, (2, 13))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[2, 15]], [[0, 13]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (13 aligned letters; 13 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 13:
        identities = 13,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 13)
        self.assertEqual(counts.identities, 13)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 13.0; substitution score = 13.0; gap score = 0.0; 13 aligned letters; 13 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 13.0:
        substitution_score = 13.0,
        gap_score = 0.0.
    aligned = 13:
        identities = 13,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 13)
        self.assertEqual(counts.identities, 13)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 13.0)
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 13.0)
        self.assertEqual(
            str(alignment),
            """\
target           21 CCCCAAAABBBAA 34
                  0 ||||||||||||| 13
query             9 CCCCAAAABBBAA 22
""",
        )
        self.assertEqual(alignment.shape, (2, 13))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[21, 34]], [[9, 22]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (13 aligned letters; 13 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 13:
        identities = 13,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 13)
        self.assertEqual(counts.identities, 13)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 13.0; substitution score = 13.0; gap score = 0.0; 13 aligned letters; 13 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 13.0:
        substitution_score = 13.0,
        gap_score = 0.0.
    aligned = 13:
        identities = 13,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 13)
        self.assertEqual(counts.identities, 13)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 13.0)
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (2 alignments; score=13) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 13)
        self.assertEqual(
            str(alignment),
            """\
target            2 AABBBAAAACCCC 15
                  0 ||||||||||||| 13
query            22 AABBBAAAACCCC  9
""",
        )
        self.assertEqual(alignment.shape, (2, 13))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[2, 15]], [[22, 9]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (13 aligned letters; 13 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 13:
        identities = 13,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 13)
        self.assertEqual(counts.identities, 13)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 13.0; substitution score = 13.0; gap score = 0.0; 13 aligned letters; 13 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 13.0:
        substitution_score = 13.0,
        gap_score = 0.0.
    aligned = 13:
        identities = 13,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 13)
        self.assertEqual(counts.identities, 13)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 13.0)
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 13)
        self.assertEqual(
            str(alignment),
            """\
target           21 CCCCAAAABBBAA 34
                  0 ||||||||||||| 13
query            13 CCCCAAAABBBAA  0
""",
        )
        self.assertEqual(alignment.shape, (2, 13))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[21, 34]], [[13, 0]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (13 aligned letters; 13 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 13:
        identities = 13,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 13)
        self.assertEqual(counts.identities, 13)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 13.0; substitution score = 13.0; gap score = 0.0; 13 aligned letters; 13 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 13.0:
        substitution_score = 13.0,
        gap_score = 0.0.
    aligned = 13:
        identities = 13,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 13)
        self.assertEqual(counts.identities, 13)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 13.0)

    def test_gap_here_only_local_2(self):
        # Force a bad alignment.
        #
        # Forces a bad alignment by having a very expensive gap penalty
        # where one would normally expect a gap, and a cheap gap penalty
        # in another place.
        seq1 = "AAAABBBAAAACCCCCCCCCCCCCCAAAABBBAAAA"
        seq2 = "AABBBAAAACCCCAAAABBBAA"
        breaks = [0, 3, len(seq2)]

        # Very expensive to open a gap in seq1:
        def nogaps(x, y):
            return -2000 - y

        # Very expensive to open a gap in seq2 unless it is in one of the allowed positions:
        def specificgaps(x, y):
            if x in breaks:
                return -2 - y
            else:
                return -2000 - y

        aligner = Align.PairwiseAligner()
        aligner.mode = "local"
        aligner.match_score = 1
        aligner.mismatch_score = -1
        aligner.insertion_score = nogaps
        aligner.deletion_score = specificgaps
        self.assertEqual(
            str(aligner),
            f"""Pairwise sequence aligner with parameters
  wildcard: None
  match_score: 1.000000
  mismatch_score: -1.000000
  insertion_score_function: {nogaps}
  deletion_score_function: {specificgaps}
  mode: local
""",
        )
        self.assertEqual(
            aligner.algorithm, "Waterman-Smith-Beyer local alignment algorithm"
        )
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 13)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, 13)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (2 alignments; score=13) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 13)
        self.assertEqual(
            str(alignment),
            """\
target            2 AABBBAAAACCCC 15
                  0 ||||||||||||| 13
query             0 AABBBAAAACCCC 13
""",
        )
        self.assertEqual(alignment.shape, (2, 13))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[2, 15]], [[0, 13]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (13 aligned letters; 13 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 13:
        identities = 13,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 13)
        self.assertEqual(counts.identities, 13)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 13.0; substitution score = 13.0; gap score = 0.0; 13 aligned letters; 13 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 13.0:
        substitution_score = 13.0,
        gap_score = 0.0.
    aligned = 13:
        identities = 13,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 13)
        self.assertEqual(counts.identities, 13)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 13.0)
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 13.0)
        self.assertEqual(
            str(alignment),
            """\
target           21 CCCCAAAABBBAA 34
                  0 ||||||||||||| 13
query             9 CCCCAAAABBBAA 22
""",
        )
        self.assertEqual(alignment.shape, (2, 13))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[21, 34]], [[9, 22]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (13 aligned letters; 13 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 13:
        identities = 13,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 13)
        self.assertEqual(counts.identities, 13)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 13.0; substitution score = 13.0; gap score = 0.0; 13 aligned letters; 13 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 13.0:
        substitution_score = 13.0,
        gap_score = 0.0.
    aligned = 13:
        identities = 13,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 13)
        self.assertEqual(counts.identities, 13)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 13.0)
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (2 alignments; score=13) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 13)
        self.assertEqual(
            str(alignment),
            """\
target            2 AABBBAAAACCCC 15
                  0 ||||||||||||| 13
query            22 AABBBAAAACCCC  9
""",
        )
        self.assertEqual(alignment.shape, (2, 13))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[2, 15]], [[22, 9]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (13 aligned letters; 13 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 13:
        identities = 13,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 13)
        self.assertEqual(counts.identities, 13)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 13.0; substitution score = 13.0; gap score = 0.0; 13 aligned letters; 13 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 13.0:
        substitution_score = 13.0,
        gap_score = 0.0.
    aligned = 13:
        identities = 13,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 13)
        self.assertEqual(counts.identities, 13)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 13.0)
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 13)
        self.assertEqual(
            str(alignment),
            """\
target           21 CCCCAAAABBBAA 34
                  0 ||||||||||||| 13
query            13 CCCCAAAABBBAA  0
""",
        )
        self.assertEqual(alignment.shape, (2, 13))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[21, 34]], [[13, 0]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (13 aligned letters; 13 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 13:
        identities = 13,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 13)
        self.assertEqual(counts.identities, 13)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 13.0; substitution score = 13.0; gap score = 0.0; 13 aligned letters; 13 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 13.0:
        substitution_score = 13.0,
        gap_score = 0.0.
    aligned = 13:
        identities = 13,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 13)
        self.assertEqual(counts.identities, 13)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 13.0)

    def test_gap_here_only_local_3(self):
        # Check if gap open and gap extend penalties are handled correctly.
        seq1 = "TTCCAA"
        seq2 = "TTGGAA"

        def gap_score(i, n):
            if i == 3:
                return -10
            if n == 1:
                return -1
            return -10

        aligner = Align.PairwiseAligner()
        aligner.mode = "local"
        aligner.match_score = 1
        aligner.mismatch_score = -10
        aligner.insertion_score = gap_score
        self.assertEqual(
            aligner.algorithm, "Waterman-Smith-Beyer local alignment algorithm"
        )
        self.assertEqual(
            str(aligner),
            f"""Pairwise sequence aligner with parameters
  wildcard: None
  match_score: 1.000000
  mismatch_score: -10.000000
  insertion_score_function: {gap_score}
  open_internal_deletion_score: 0.000000
  extend_internal_deletion_score: 0.000000
  open_left_deletion_score: 0.000000
  extend_left_deletion_score: 0.000000
  open_right_deletion_score: 0.000000
  extend_right_deletion_score: 0.000000
  mode: local
""",
        )
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 2.0)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, 2.0)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (2 alignments; score=2) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 2.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 TT 2
                  0 || 2
query             0 TT 2
""",
        )
        self.assertEqual(alignment.shape, (2, 2))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[0, 2]], [[0, 2]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 2 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 2.0; substitution score = 2.0; gap score = 0.0; 2 aligned letters; 2 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 2.0:
        substitution_score = 2.0,
        gap_score = 0.0.
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 2.0)
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 2.0)
        self.assertEqual(
            str(alignment),
            """\
target            4 AA 6
                  0 || 2
query             4 AA 6
""",
        )
        self.assertEqual(alignment.shape, (2, 2))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[4, 6]], [[4, 6]]]))
        )
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (2 alignments; score=2) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 2.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 TT 2
                  0 || 2
query             6 TT 4
""",
        )
        self.assertEqual(alignment.shape, (2, 2))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[0, 2]], [[6, 4]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 2 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 2.0; substitution score = 2.0; gap score = 0.0; 2 aligned letters; 2 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 2.0:
        substitution_score = 2.0,
        gap_score = 0.0.
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 2.0)
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 2.0)
        self.assertEqual(
            str(alignment),
            """\
target            4 AA 6
                  0 || 2
query             2 AA 0
""",
        )
        self.assertEqual(alignment.shape, (2, 2))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[4, 6]], [[2, 0]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 2 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 2.0; substitution score = 2.0; gap score = 0.0; 2 aligned letters; 2 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 2.0:
        substitution_score = 2.0,
        gap_score = 0.0.
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 2.0)
        aligner.deletion_score = gap_score
        self.assertEqual(
            str(aligner),
            f"""Pairwise sequence aligner with parameters
  wildcard: None
  match_score: 1.000000
  mismatch_score: -10.000000
  insertion_score_function: {gap_score}
  deletion_score_function: {gap_score}
  mode: local
""",
        )
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 2.0)
        score = aligner.score(seq1, reverse_complement(seq2), strand="-")
        self.assertAlmostEqual(score, 2.0)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (2 alignments; score=2) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 2.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 TT 2
                  0 || 2
query             0 TT 2
""",
        )
        self.assertEqual(alignment.shape, (2, 2))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[0, 2]], [[0, 2]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 2 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 2.0; substitution score = 2.0; gap score = 0.0; 2 aligned letters; 2 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 2.0:
        substitution_score = 2.0,
        gap_score = 0.0.
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 2.0)
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 2.0)
        self.assertEqual(
            str(alignment),
            """\
target            4 AA 6
                  0 || 2
query             4 AA 6
""",
        )
        self.assertEqual(alignment.shape, (2, 2))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[4, 6]], [[4, 6]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 2 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 2.0; substitution score = 2.0; gap score = 0.0; 2 aligned letters; 2 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 2.0:
        substitution_score = 2.0,
        gap_score = 0.0.
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 2.0)
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (2 alignments; score=2) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 2.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 TT 2
                  0 || 2
query             6 TT 4
""",
        )
        self.assertEqual(alignment.shape, (2, 2))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[0, 2]], [[6, 4]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 2 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 2.0; substitution score = 2.0; gap score = 0.0; 2 aligned letters; 2 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 2.0:
        substitution_score = 2.0,
        gap_score = 0.0.
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 2.0)
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 2.0)
        self.assertEqual(
            str(alignment),
            """\
target            4 AA 6
                  0 || 2
query             2 AA 0
""",
        )
        self.assertEqual(alignment.shape, (2, 2))
        self.assertTrue(
            np.array_equal(alignment.aligned, np.array([[[4, 6]], [[2, 0]]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (2 aligned letters; 2 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 2.0; substitution score = 2.0; gap score = 0.0; 2 aligned letters; 2 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 2.0:
        substitution_score = 2.0,
        gap_score = 0.0.
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 2.0)

    def test_broken_gap_function(self):
        # Check if an Exception is propagated if the gap function raises one
        seq1 = "TTCCAA"
        seq2 = "TTGGAA"

        def gap_score(i, n):
            raise RuntimeError("broken gap function")

        aligner = Align.PairwiseAligner()
        aligner.insertion_score = gap_score
        aligner.deletion_score = -1
        aligner.mode = "global"
        with self.assertRaises(RuntimeError):
            aligner.score(seq1, seq2)
        with self.assertRaises(RuntimeError):
            aligner.score(seq1, reverse_complement(seq2), strand="-")
        with self.assertRaises(RuntimeError):
            alignments = aligner.align(seq1, seq2)
            alignments = list(alignments)
        with self.assertRaises(RuntimeError):
            alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
            alignments = list(alignments)
        aligner.mode = "local"
        with self.assertRaises(RuntimeError):
            aligner.score(seq1, seq2)
        with self.assertRaises(RuntimeError):
            aligner.score(seq1, reverse_complement(seq2), strand="-")
        with self.assertRaises(RuntimeError):
            alignments = aligner.align(seq1, seq2)
            alignments = list(alignments)
        with self.assertRaises(RuntimeError):
            alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
            alignments = list(alignments)
        aligner.insertion_score = -1
        aligner.deletion_score = gap_score
        aligner.mode = "global"
        with self.assertRaises(RuntimeError):
            aligner.score(seq1, seq2)
        with self.assertRaises(RuntimeError):
            aligner.score(seq1, reverse_complement(seq2), strand="-")
        with self.assertRaises(RuntimeError):
            alignments = aligner.align(seq1, seq2)
            alignments = list(alignments)
        with self.assertRaises(RuntimeError):
            alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
            alignments = list(alignments)
        aligner.mode = "local"
        with self.assertRaises(RuntimeError):
            aligner.score(seq1, seq2)
        with self.assertRaises(RuntimeError):
            aligner.score(seq1, reverse_complement(seq2), strand="-")
        with self.assertRaises(RuntimeError):
            alignments = aligner.align(seq1, seq2)
            alignments = list(alignments)
        with self.assertRaises(RuntimeError):
            alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
            alignments = list(alignments)

    def test_mixing_affine_and_function(self):
        # Check if we can use a gap function for one sequence, and an affine
        # gap score for the other sequence.

        def gap_score_function(position, length):
            if position == 4:
                return 0
            else:
                return -100

        aligner = Align.PairwiseAligner()
        aligner.mismatch_score = -1
        aligner.deletion_score = gap_score_function
        aligner.internal_insertion_score = -10
        seqA = "AAAAAAAAAAA"
        seqB = "TTAAAAA"
        alignments = aligner.align(seqA, seqB)
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 5.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 --AAAAAAAAAAA 11
                  0 --||------||| 13
query             0 TTAA------AAA  7
""",
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 5.0; substitution score = 5.0; gap score = 0.0; 5 aligned letters; 5 identities; 0 mismatches; 8 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 5.0:
        substitution_score = 5.0,
        gap_score = 0.0.
    aligned = 5:
        identities = 5,
        mismatches = 0.
    gaps = 8:
        left_gaps = 2:
            left_insertions = 2:
                open_left_insertions = 1,
                extend_left_insertions = 1;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 6:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 6:
                open_internal_deletions = 1,
                extend_internal_deletions = 5;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 5)
        self.assertEqual(counts.mismatches, 0)
        self.assertEqual(counts.insertions, 2)
        self.assertEqual(counts.deletions, 6)
        self.assertAlmostEqual(counts.score, 5.0)
        seqA = "AAAAAAAAAAA"
        seqB = "AAAAATT"
        alignments = aligner.align(seqA, seqB)
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 5.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 AAAAAAAAAAA-- 11
                  0 ||||------|-- 13
query             0 AAAA------ATT  7
""",
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 5.0; substitution score = 5.0; gap score = 0.0; 5 aligned letters; 5 identities; 0 mismatches; 8 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 5.0:
        substitution_score = 5.0,
        gap_score = 0.0.
    aligned = 5:
        identities = 5,
        mismatches = 0.
    gaps = 8:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 6:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 6:
                open_internal_deletions = 1,
                extend_internal_deletions = 5;
        right_gaps = 2:
            right_insertions = 2:
                open_right_insertions = 1,
                extend_right_insertions = 1;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 5)
        self.assertEqual(counts.mismatches, 0)
        self.assertEqual(counts.insertions, 2)
        self.assertEqual(counts.deletions, 6)
        self.assertAlmostEqual(counts.score, 5.0)

        # Same thing, switching A and B:
        aligner = Align.PairwiseAligner()
        aligner.mismatch_score = -1
        aligner.insertion_score = gap_score_function
        aligner.internal_deletion_score = -10
        seqA = "TTAAAAA"
        seqB = "AAAAAAAAAAA"
        alignments = aligner.align(seqA, seqB)
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 5.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 TTAA------AAA  7
                  0 --||------||| 13
query             0 --AAAAAAAAAAA 11
""",
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 5.0; substitution score = 5.0; gap score = 0.0; 5 aligned letters; 5 identities; 0 mismatches; 8 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 5.0:
        substitution_score = 5.0,
        gap_score = 0.0.
    aligned = 5:
        identities = 5,
        mismatches = 0.
    gaps = 8:
        left_gaps = 2:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 2:
                open_left_deletions = 1,
                extend_left_deletions = 1;
        internal_gaps = 6:
            internal_insertions = 6:
                open_internal_insertions = 1,
                extend_internal_insertions = 5;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 5)
        self.assertEqual(counts.mismatches, 0)
        self.assertEqual(counts.insertions, 6)
        self.assertEqual(counts.deletions, 2)
        self.assertAlmostEqual(counts.score, 5.0)
        seqA = "AAAAATT"
        seqB = "AAAAAAAAAAA"
        alignments = aligner.align(seqA, seqB)
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 5.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 AAAA------ATT  7
                  0 ||||------|-- 13
query             0 AAAAAAAAAAA-- 11
""",
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 5.0; substitution score = 5.0; gap score = 0.0; 5 aligned letters; 5 identities; 0 mismatches; 8 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 5.0:
        substitution_score = 5.0,
        gap_score = 0.0.
    aligned = 5:
        identities = 5,
        mismatches = 0.
    gaps = 8:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 6:
            internal_insertions = 6:
                open_internal_insertions = 1,
                extend_internal_insertions = 5;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 2:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 2:
                open_right_deletions = 1,
                extend_right_deletions = 1.
""",
        )
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 5)
        self.assertEqual(counts.mismatches, 0)
        self.assertEqual(counts.insertions, 6)
        self.assertEqual(counts.deletions, 2)
        self.assertAlmostEqual(counts.score, 5.0)


class TestAlignerInput(unittest.TestCase):
    """Check aligning sequences provided as lists, str, Seq, or SeqRecord objects."""

    def test_three_letter_amino_acids_global(self):
        """Test aligning sequences provided as lists of three-letter amino acids."""
        seq1 = ["Gly", "Ala", "Thr"]
        seq2 = ["Gly", "Ala", "Ala", "Cys", "Thr"]
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 3.0)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (2 alignments; score=3) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 2)
        self.assertEqual(
            str(alignments[0]),
            """\
Gly Ala --- --- Thr
||| ||| --- --- |||
Gly Ala Ala Cys Thr
""",
        )
        self.assertEqual(
            str(alignments[1]),
            """\
Gly --- Ala --- Thr
||| --- ||| --- |||
Gly Ala Ala Cys Thr
""",
        )
        self.assertAlmostEqual(alignments[0].score, 3.0)
        self.assertAlmostEqual(alignments[1].score, 3.0)
        counts = alignments[0].counts()
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignments[0].counts(aligner)
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 3.0)
        counts = alignments[1].counts(aligner)
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertAlmostEqual(counts.score, 3.0)

        seq1 = ["Pro", "Pro", "Gly", "Ala", "Thr"]
        seq2 = ["Gly", "Ala", "Ala", "Cys", "Thr", "Asn", "Asn"]
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 3.0)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (2 alignments; score=3) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
Pro Pro Gly Ala --- --- Thr --- ---
--- --- ||| ||| --- --- ||| --- ---
--- --- Gly Ala Ala Cys Thr Asn Asn
""",
        )
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(
            alignment[0], ["Pro", "Pro", "Gly", "Ala", None, None, "Thr", None, None]
        )
        self.assertEqual(
            alignment[0, :], ["Pro", "Pro", "Gly", "Ala", None, None, "Thr", None, None]
        )
        self.assertEqual(
            alignment[0, 1:], ["Pro", "Gly", "Ala", None, None, "Thr", None, None]
        )
        self.assertEqual(alignment[0, ::2], ["Pro", "Gly", None, "Thr", None])
        self.assertEqual(
            alignment[1], [None, None, "Gly", "Ala", "Ala", "Cys", "Thr", "Asn", "Asn"]
        )
        self.assertEqual(
            alignment[1, :],
            [None, None, "Gly", "Ala", "Ala", "Cys", "Thr", "Asn", "Asn"],
        )
        self.assertEqual(
            alignment[1, 1:], [None, "Gly", "Ala", "Ala", "Cys", "Thr", "Asn", "Asn"]
        )
        self.assertEqual(alignment[1, ::2], [None, "Gly", "Ala", "Thr", "Asn"])
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (3 aligned letters; 3 identities; 0 mismatches; 6 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 3:
        identities = 3,
        mismatches = 0.
    gaps = 6:
        left_gaps = 2:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 2:
                open_left_deletions = 1,
                extend_left_deletions = 1;
        internal_gaps = 2:
            internal_insertions = 2:
                open_internal_insertions = 1,
                extend_internal_insertions = 1;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 2:
            right_insertions = 2:
                open_right_insertions = 1,
                extend_right_insertions = 1;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.score)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 3.0; substitution score = 3.0; gap score = 0.0; 3 aligned letters; 3 identities; 0 mismatches; 6 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 3.0:
        substitution_score = 3.0,
        gap_score = 0.0.
    aligned = 3:
        identities = 3,
        mismatches = 0.
    gaps = 6:
        left_gaps = 2:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 2:
                open_left_deletions = 1,
                extend_left_deletions = 1;
        internal_gaps = 2:
            internal_insertions = 2:
                open_internal_insertions = 1,
                extend_internal_insertions = 1;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 2:
            right_insertions = 2:
                open_right_insertions = 1,
                extend_right_insertions = 1;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertEqual(counts.deletions, 2)
        self.assertEqual(counts.insertions, 4)
        self.assertAlmostEqual(counts.score, 3.0)
        alignment = alignments[1]
        self.assertEqual(
            str(alignment),
            """\
Pro Pro Gly --- Ala --- Thr --- ---
--- --- ||| --- ||| --- ||| --- ---
--- --- Gly Ala Ala Cys Thr Asn Asn
""",
        )
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(
            alignment[0], ["Pro", "Pro", "Gly", None, "Ala", None, "Thr", None, None]
        )
        self.assertEqual(
            alignment[0, :], ["Pro", "Pro", "Gly", None, "Ala", None, "Thr", None, None]
        )
        self.assertEqual(
            alignment[0, 1:-1], ["Pro", "Gly", None, "Ala", None, "Thr", None]
        )
        self.assertEqual(alignment[0, 1::2], ["Pro", None, None, None])
        self.assertEqual(
            alignment[1], [None, None, "Gly", "Ala", "Ala", "Cys", "Thr", "Asn", "Asn"]
        )
        self.assertEqual(
            alignment[1, :],
            [None, None, "Gly", "Ala", "Ala", "Cys", "Thr", "Asn", "Asn"],
        )
        self.assertEqual(
            alignment[1, 1:-1], [None, "Gly", "Ala", "Ala", "Cys", "Thr", "Asn"]
        )
        self.assertEqual(alignment[1, 1::2], [None, "Ala", "Cys", "Asn"])
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 3.0; substitution score = 3.0; gap score = 0.0; 3 aligned letters; 3 identities; 0 mismatches; 6 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 3.0:
        substitution_score = 3.0,
        gap_score = 0.0.
    aligned = 3:
        identities = 3,
        mismatches = 0.
    gaps = 6:
        left_gaps = 2:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 2:
                open_left_deletions = 1,
                extend_left_deletions = 1;
        internal_gaps = 2:
            internal_insertions = 2:
                open_internal_insertions = 2,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 2:
            right_insertions = 2:
                open_right_insertions = 1,
                extend_right_insertions = 1;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertEqual(counts.deletions, 2)
        self.assertEqual(counts.insertions, 4)
        self.assertAlmostEqual(counts.score, 3.0)

    def test_three_letter_amino_acids_local(self):
        seq1 = ["Asn", "Asn", "Gly", "Ala", "Thr", "Glu", "Glu"]
        seq2 = ["Pro", "Pro", "Gly", "Ala", "Ala", "Cys", "Thr", "Leu"]
        aligner = Align.PairwiseAligner()
        aligner.mode = "local"
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 3.0)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (2 alignments; score=3) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
Gly Ala --- --- Thr
||| ||| --- --- |||
Gly Ala Ala Cys Thr
""",
        )
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(alignment[0], ["Gly", "Ala", None, None, "Thr"])
        self.assertEqual(alignment[0, :], ["Gly", "Ala", None, None, "Thr"])
        self.assertEqual(alignment[0, 1:], ["Ala", None, None, "Thr"])
        self.assertEqual(alignment[0, :-1], ["Gly", "Ala", None, None])
        self.assertEqual(alignment[0, ::2], ["Gly", None, "Thr"])
        self.assertEqual(alignment[1], ["Gly", "Ala", "Ala", "Cys", "Thr"])
        self.assertEqual(alignment[1, :], ["Gly", "Ala", "Ala", "Cys", "Thr"])
        self.assertEqual(alignment[1, 1:], ["Ala", "Ala", "Cys", "Thr"])
        self.assertEqual(alignment[1, :-1], ["Gly", "Ala", "Ala", "Cys"])
        self.assertEqual(alignment[1, ::2], ["Gly", "Ala", "Thr"])
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 3.0; substitution score = 3.0; gap score = 0.0; 3 aligned letters; 3 identities; 0 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 3.0:
        substitution_score = 3.0,
        gap_score = 0.0.
    aligned = 3:
        identities = 3,
        mismatches = 0.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 2:
                open_internal_insertions = 1,
                extend_internal_insertions = 1;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.insertions, 2)
        self.assertAlmostEqual(counts.score, 3.0)
        alignment = alignments[1]
        self.assertEqual(
            str(alignment),
            """\
Gly --- Ala --- Thr
||| --- ||| --- |||
Gly Ala Ala Cys Thr
""",
        )
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(alignment[0], ["Gly", None, "Ala", None, "Thr"])
        self.assertEqual(alignment[0, :], ["Gly", None, "Ala", None, "Thr"])
        self.assertEqual(alignment[0, 1:], [None, "Ala", None, "Thr"])
        self.assertEqual(alignment[0, :-1], ["Gly", None, "Ala", None])
        self.assertEqual(alignment[0, ::2], ["Gly", "Ala", "Thr"])
        self.assertEqual(alignment[1], ["Gly", "Ala", "Ala", "Cys", "Thr"])
        self.assertEqual(alignment[1, :], ["Gly", "Ala", "Ala", "Cys", "Thr"])
        self.assertEqual(alignment[1, 1:], ["Ala", "Ala", "Cys", "Thr"])
        self.assertEqual(alignment[1, :-1], ["Gly", "Ala", "Ala", "Cys"])
        self.assertEqual(alignment[1, ::2], ["Gly", "Ala", "Thr"])
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 3.0; substitution score = 3.0; gap score = 0.0; 3 aligned letters; 3 identities; 0 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 3.0:
        substitution_score = 3.0,
        gap_score = 0.0.
    aligned = 3:
        identities = 3,
        mismatches = 0.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 2:
                open_internal_insertions = 2,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.insertions, 2)
        self.assertAlmostEqual(counts.score, 3.0)

    def test_str_seq_seqrecord(self):
        """Test aligning sequences provided as str, Seq, or SeqRecord objects."""
        aligner = Align.PairwiseAligner("blastn")
        t1 = "ACGT"
        t2 = "CGTT"
        s1 = Seq(t1)
        s2 = Seq(t2)
        r1 = SeqRecord(s1, id="first", description="1st sequence")
        r2 = SeqRecord(s2, id="second", description="2nd sequence")
        alignments = aligner.align(t1, t2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=-7) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
target            0 ACGT 4
                  0 ...| 4
query             0 CGTT 4
""",
        )
        self.assertEqual(
            format(alignment, "fasta"),
            """\
>
ACGT
>
CGTT
""",
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = -7.0; substitution score = -7.0; gap score = 0.0; 4 aligned letters; 1 identities; 3 mismatches; 1 positives; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = -7.0:
        substitution_score = -7.0,
        gap_score = 0.0.
    aligned = 4:
        identities = 1,
        positives = 1,
        mismatches = 3.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 3)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.positives, 1)
        self.assertAlmostEqual(counts.score, -7.0)
        alignments = aligner.align(s1, s2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=-7) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
target            0 ACGT 4
                  0 ...| 4
query             0 CGTT 4
""",
        )
        self.assertEqual(
            format(alignment, "fasta"),
            """\
>
ACGT
>
CGTT
""",
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = -7.0; substitution score = -7.0; gap score = 0.0; 4 aligned letters; 1 identities; 3 mismatches; 1 positives; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = -7.0:
        substitution_score = -7.0,
        gap_score = 0.0.
    aligned = 4:
        identities = 1,
        positives = 1,
        mismatches = 3.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 3)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.positives, 1)
        self.assertAlmostEqual(counts.score, -7.0)
        alignments = aligner.align(r1, r2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=-7) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
first             0 ACGT 4
                  0 ...| 4
second            0 CGTT 4
""",
        )
        self.assertEqual(
            format(alignment, "fasta"),
            """\
>first 1st sequence
ACGT
>second 2nd sequence
CGTT
""",
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = -7.0; substitution score = -7.0; gap score = 0.0; 4 aligned letters; 1 identities; 3 mismatches; 1 positives; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = -7.0:
        substitution_score = -7.0,
        gap_score = 0.0.
    aligned = 4:
        identities = 1,
        positives = 1,
        mismatches = 3.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 4)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 3)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertAlmostEqual(counts.score, -7.0)


class TestArgumentErrors(unittest.TestCase):
    def test_aligner_string_errors(self):
        aligner = Align.PairwiseAligner()
        message = "^'int' object is not iterable$"
        with self.assertRaisesRegex(TypeError, message):
            aligner.score("AAA", 3)
        message = "^sequence has zero length$"
        with self.assertRaisesRegex(ValueError, message):
            aligner.score("AAA", "")
        with self.assertRaisesRegex(ValueError, message):
            aligner.score("AAA", "", strand="-")

    def test_aligner_array_errors(self):
        aligner = Align.PairwiseAligner()
        s1 = "GGG"
        s2 = array.array("i", [ord("G"), ord("A"), ord("G")])
        score = aligner.score(s1, s2)
        self.assertAlmostEqual(score, 2.0)
        s2 = array.array("f", [1.0, 0.0, 1.0])
        message = "^sequence has incorrect data type 'f'$"
        with self.assertRaisesRegex(ValueError, message):
            aligner.score(s1, s2)
        aligner.wildcard = chr(99)
        s1 = array.array("i", [1, 5, 6])
        s2 = array.array("i", [1, 8, 6])
        s2a = array.array("i", [1, 8, 99])
        s2b = array.array("i", [1, 28, 6])
        aligner.match = 3.0
        aligner.mismatch = -2.0
        aligner.gap_score = -10.0
        score = aligner.score(s1, s2)
        self.assertAlmostEqual(score, 4.0)
        # the following two are valid as we are using match/mismatch scores
        # instead of a substitution matrix:
        score = aligner.score(s1, s2a)
        # since we set the wildcard character to chr(99), the number 99
        # is interpreted as an unknown character, and gets a zero score:
        self.assertAlmostEqual(score, 1.0)
        score = aligner.score(s1, s2b)
        self.assertAlmostEqual(score, 4.0)
        try:
            import numpy as np
        except ImportError:
            return
        aligner = Align.PairwiseAligner()
        aligner.wildcard = chr(99)
        s1 = "GGG"
        s2 = np.array([ord("G"), ord("A"), ord("G")], np.int32)
        score = aligner.score(s1, s2)
        self.assertAlmostEqual(score, 2.0)
        alignments = aligner.align(s1, s2)
        self.assertEqual(len(alignments), 5)
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
G  -- G  G
.  -- .  -
71 65 71 -
""",
        )
        self.assertAlmostEqual(alignment.score, 2.0)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 2.0; substitution score = 2.0; gap score = 0.0; 2 aligned letters; 2 identities; 0 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 2.0:
        substitution_score = 2.0,
        gap_score = 0.0.
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 1:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 1:
                open_right_deletions = 1,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertEqual(counts.deletions, 1)
        self.assertEqual(counts.insertions, 1)
        self.assertAlmostEqual(counts.score, 2.0)

        s2 = np.array([1.0, 0.0, 1.0])
        message = "^sequence has incorrect data type 'd'$"
        with self.assertRaisesRegex(ValueError, message):
            aligner.score(s1, s2)
        with self.assertRaisesRegex(ValueError, message):
            aligner.align(s1, s2)
        s2 = np.zeros((3, 2), np.int32)
        message = "^sequence has incorrect rank \\(2 expected 1\\)$"
        with self.assertRaisesRegex(ValueError, message):
            aligner.score(s1, s2)
        with self.assertRaisesRegex(ValueError, message):
            aligner.align(s1, s2)
        s1 = np.array([1, 5, 6], np.int32)
        s2 = np.array([1, 8, 6], np.int32)
        s2a = np.array([1, 8, 99], np.int32)
        s2b = np.array([1, 28, 6], np.int32)
        s2c = np.array([1, 8, -6], np.int32)
        aligner.match = 3.0
        aligner.mismatch = -2.0
        aligner.gap_score = -10.0
        score = aligner.score(s1, s2)
        self.assertAlmostEqual(score, 4.0)
        alignments = aligner.align(s1, s2)
        self.assertAlmostEqual(alignments.score, 4.0)
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 4.0)
        self.assertEqual(
            str(alignment),
            """\
1 5 6
| . |
1 8 6
""",
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 4.0; substitution score = 4.0; gap score = 0.0; 3 aligned letters; 2 identities; 1 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 4.0:
        substitution_score = 4.0,
        gap_score = 0.0.
    aligned = 3:
        identities = 2,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 1)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertAlmostEqual(counts.score, 4.0)
        # alignments are valid as we are using match/mismatch scores
        # instead of a substitution matrix:
        score = aligner.score(s1, s2a)
        self.assertAlmostEqual(score, 1.0)
        alignments = aligner.align(s1, s2a)
        self.assertEqual(len(alignments), 1)
        self.assertAlmostEqual(alignments.score, 1.0)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1.0)
        self.assertEqual(
            str(alignment),
            """\
1 5 6 
| . . 
1 8 99
""",
        )
        score = aligner.score(s1, s2b)
        self.assertAlmostEqual(score, 4.0)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.0; substitution score = 1.0; gap score = 0.0; 3 aligned letters; 1 identities; 1 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.0:
        substitution_score = 1.0,
        gap_score = 0.0.
    aligned = 3:
        identities = 1,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 1)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertAlmostEqual(counts.score, 1.0)
        alignments = aligner.align(s1, s2b)
        self.assertEqual(len(alignments), 1)
        self.assertAlmostEqual(alignments.score, 4.0)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 4.0)
        self.assertEqual(
            str(alignment),
            """\
1 5  6
| .  |
1 28 6
""",
        )
        self.assertAlmostEqual(score, 4.0)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 4.0; substitution score = 4.0; gap score = 0.0; 3 aligned letters; 2 identities; 1 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 4.0:
        substitution_score = 4.0,
        gap_score = 0.0.
    aligned = 3:
        identities = 2,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 1)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertAlmostEqual(counts.score, 4.0)
        # when using a substitution matrix, all indices should be between 0
        # and the size of the substitution matrix:
        m = 5 * np.eye(10)
        aligner.substitution_matrix = m
        score = aligner.score(s1, s2)  # no ValueError
        self.assertAlmostEqual(score, 10.0)
        alignments = aligner.align(s1, s2)
        self.assertAlmostEqual(alignments.score, 10.0)
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 10.0)
        self.assertEqual(
            str(alignment),
            """\
1 5 6
| . |
1 8 6
""",
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 10.0; substitution score = 10.0; gap score = 0.0; 3 aligned letters; 2 identities; 1 mismatches; 2 positives; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 10.0:
        substitution_score = 10.0,
        gap_score = 0.0.
    aligned = 3:
        identities = 2,
        positives = 2,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 1)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertAlmostEqual(counts.score, 10.0)
        counts = alignment.counts(m)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (substitution score = 10.0; 3 aligned letters; 2 identities; 1 mismatches; 2 positives; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    substitution_score = 10.0,
    aligned = 3:
        identities = 2,
        positives = 2,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        alignment.sequences[1] = s2b
        with self.assertRaises(ValueError) as cm:
            alignment.counts(m)
        self.assertEqual(
            str(cm.exception), "sequence[1][1] is out of bound (28, should be < 10)"
        )
        with self.assertRaises(ValueError) as cm:
            alignment.counts(aligner)
        self.assertEqual(
            str(cm.exception), "sequence[1][1] is out of bound (28, should be < 10)"
        )
        alignment.sequences[1] = s2c
        with self.assertRaises(ValueError) as cm:
            alignment.counts(m)
        self.assertEqual(str(cm.exception), "sequences[1][2] is negative (-6)")
        with self.assertRaises(ValueError) as cm:
            alignment.counts(aligner)
        self.assertEqual(str(cm.exception), "sequences[1][2] is negative (-6)")
        message = "^sequence item 2 is negative \\(-6\\)$"
        with self.assertRaisesRegex(ValueError, message):
            aligner.score(s1, s2c)
        with self.assertRaisesRegex(ValueError, message):
            aligner.align(s1, s2c)
        message = "^sequence item 1 is out of bound \\(28, should be < 10\\)$"
        with self.assertRaisesRegex(ValueError, message):
            aligner.score(s1, s2b)
        with self.assertRaisesRegex(ValueError, message):
            aligner.align(s1, s2b)
        # note that the wildcard character is ignored when using a substitution
        # matrix, so 99 is interpreted as an index here:
        message = "^sequence item 2 is out of bound \\(99, should be < 10\\)$"
        with self.assertRaisesRegex(ValueError, message):
            aligner.score(s1, s2a)
        with self.assertRaisesRegex(ValueError, message):
            aligner.align(s1, s2a)


class TestOverflowError(unittest.TestCase):
    def test_align_overflow_error(self):
        aligner = Align.PairwiseAligner()
        path = os.path.join("Align", "bsubtilis.fa")
        record = SeqIO.read(path, "fasta")
        seq1 = record.seq
        path = os.path.join("Align", "ecoli.fa")
        record = SeqIO.read(path, "fasta")
        seq2 = record.seq
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (>{sys.maxsize} alignments; score=1286) at {hex(id(alignments))}>""",
        )
        self.assertAlmostEqual(alignments.score, 1286.0)
        message = "^number of optimal alignments is larger than (%d|%d)$" % (
            2147483647,  # on 32-bit systems
            9223372036854775807,
        )  # on 64-bit systems
        with self.assertRaisesRegex(OverflowError, message):
            n = len(alignments)
        # confirm that we can still pull out individual alignments
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
target            0 ATTTA-TC-GGA-GAGTTTGATCC-TGGCTCAGGAC--GAACGCTGGCGGC-GTGCCTAA
                  0 |---|-|--|-|-||||||||||--||||||||-|---|||||||||||||-|-||||||
query             0 A---AAT-TG-AAGAGTTTGATC-ATGGCTCAG-A-TTGAACGCTGGCGGCAG-GCCTAA

target           53 T-ACATGCAAGTCGAG-CGG-A-CAG-AT-GGGA-GCTTGCT-C----CCTGAT-GTTAG
                 60 --|||||||||||||--|||-|-|||-|--|--|-|||||||-|----|-|||--|--||
query            51 -CACATGCAAGTCGA-ACGGTAACAGGA-AG--AAGCTTGCTTCTTTGC-TGA-CG--AG

target          100 C-GGCGGACGGGTGAGTAACAC-GT--GGGTAA-CCTGCCTGTAA-G-ACTGGG--ATAA
                120 --|||||||||||||||||----||--|||-||-|-||||||-|--|-|--|||--||||
query           102 -TGGCGGACGGGTGAGTAA---TGTCTGGG-AAAC-TGCCTG-A-TGGA--GGGGGATAA

target          151 CT-CC-GGGAAACCGG--GGCTAATACCGG-ATGGTTGTTTGAACCGCAT-GGTTCAA-A
                180 ||-|--||-||||-||--|-|||||||||--||---------|||-|--|-|---|||-|
query           152 CTAC-TGG-AAAC-GGTAG-CTAATACCG-CAT---------AAC-G--TCG---CAAGA

target          204 C-ATAA-AAGGTGG--C-TTCGG-C-TACCACTTA-C-A--G-ATG-GACCC-GC--GGC
                240 |-|-||-|-||-||--|-|||||-|-|-|---||--|-|--|-|||-|-|||-|---||-
query           192 CCA-AAGA-GG-GGGACCTTCGGGCCT-C---TT-GCCATCGGATGTG-CCCAG-ATGG-

target          248 GCATTAGCTAGTT-GGTGAGG-TAACGGCTCACC-AAGGCGACGATGCG--TAGCC-GA-
                300 |-||||||||||--||||-||-||||||||||||-|-|||||||||-|---||||--|--
query           241 G-ATTAGCTAGT-AGGTG-GGGTAACGGCTCACCTA-GGCGACGAT-C-CCTAGC-TG-G

target          301 -CCTGAGAGGG-TGATC--GGCCACACTGGGA-CTGAGACACGG-CCCAGACTCCTACGG
                360 -|-|||||||--|||-|--|-|||||||||-|-|||||||||||-||-||||||||||||
query           293 TC-TGAGAGG-ATGA-CCAG-CCACACTGG-AACTGAGACACGGTCC-AGACTCCTACGG

target          355 GAGGCAGCAGTAGGG-AATC-TTCCGCA-A-TGGA-CG-AAAGTC-TGAC-GG-AGCAAC
                420 |||||||||||-|||-|||--||--|||-|-|||--||-||-|-|-|||--|--|||--|
query           347 GAGGCAGCAGT-GGGGAAT-ATT--GCACAATGG-GCGCAA-G-CCTGA-TG-CAGC--C

target          406 --GCCGCGTG-AGTGAT-GAAGG--TTTTCGGA-TC-GTAAAGCT-CTGTTGTT-AG-GG
                480 --||||||||-|-|||--|||||--||--|||--|--||||||-|-||-||----||-||
query           396 ATGCCGCGTGTA-TGA-AGAAGGCCTT--CGG-GT-TGTAAAG-TACT-TT---CAGCGG

target          455 --A--A-G--A--ACAAGTGCCGTTCGAATAGGGC----GG-TACC-TTGACGGT-ACCT
                540 --|--|-|--|--|-||||----|---||||---|----|--|-|--||||||-|-|||-
query           445 GGAGGAAGGGAGTA-AAGT----T---AATA---CCTTTG-CT-C-ATTGACG-TTACC-

target          499 AAC-CAGAA-A-GCCAC-GGCTAACTAC-GTGCCAGCAGCCGCGGTAATACGT-AGG-TG
                600 --|-|||||-|-||-||-||||||||-|-|||||||||||||||||||||||--|||-||
query           489 --CGCAGAAGAAGC-ACCGGCTAACT-CCGTGCCAGCAGCCGCGGTAATACG-GAGGGTG

target          552 GCAAGCGTTG--TCCGGAATTA-TTGGGCGTAAAG-GGCT-CGCAGGCGGTTTC-TTAAG
                660 -||||||||---||-|||||||-|-||||||||||-|-|--||||||||||||--|||||
query           544 -CAAGCGTT-AATC-GGAATTACT-GGGCGTAAAGCG-C-ACGCAGGCGGTTT-GTTAAG

target          606 TCT-GATGTGAAAG-CCCCCGG-CTCAACC-GGGGAGGG--T-CAT-TGGA-AACTGGGG
                720 ||--|||||||||--||||-||-|||||||-|||-|-----|-|||-||-|-|-||||--
query           597 TC-AGATGTGAAA-TCCCC-GGGCTCAACCTGGG-A---ACTGCATCTG-ATA-CTGG--

target          657 -AA-CTTGAGTGCA--G-AAGAGGAGAGTGG-A-A-TTCCACG-TGTAGCGGTGAAATGC
                780 -||-|||||||-|---|-|-||||-|-|-||-|-|-|||||-|-||||||||||||||||
query           646 CAAGCTTGAGT-C-TCGTA-GAGG-G-G-GGTAGAATTCCA-GGTGTAGCGGTGAAATGC

target          708 GTAGAGATG-TGGAGGAAC-ACCAG-TGGCGAAGGCGA-CTCTC--TGGT-CTGTAA--C
                840 ||||||||--||||||||--|||-|-|||||||||||--|-|-|--|||--|-|-||--|
query           699 GTAGAGAT-CTGGAGGAA-TACC-GGTGGCGAAGGCG-GC-C-CCCTGG-AC-G-AAGAC

target          759 TGACGCTG-AGGA-GCGAAAGCGTGGGGAGCGAA-CAGGATTAGATACCCTGGTAGTCCA
                900 |||||||--|||--|||||||||||||||||-||-|||||||||||||||||||||||||
query           750 TGACGCT-CAGG-TGCGAAAGCGTGGGGAGC-AAACAGGATTAGATACCCTGGTAGTCCA

target          816 CGCCGTAAACGATGAGT-G-CTAAGTGTT-AGGGGGTT-TCCGCCCCTT-AGTGC-TG-C
                960 ||||||||||||||--|-|-||---||---|||---||-|--||||-||-||-||-||-|
query           807 CGCCGTAAACGATG--TCGACT---TG--GAGG---TTGT--GCCC-TTGAG-GCGTGGC

target          869 ------AGCTAACGCA-TTAAG-C-ACTCCGCCTGGGGAGTACGGTC-GCAAGACTG--A
               1020 ------|||||||||--|||||-|-||-|-|||||||||||||||-|-|||||---|--|
query           853 TTCCGGAGCTAACGC-GTTAAGTCGAC-C-GCCTGGGGAGTACGG-CCGCAAG---GTTA

target          917 AA-CTCAAA-GGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAA
               1080 ||-||||||-|-|||||||||||||||||||||||||||||||||||||||||||||||-
query           906 AAACTCAAATG-AATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGA-

target          975 -GCAACGCGAAGAACCTTACCA-GGTCTTGACATCCTCTGACA-A--T--CCTAGAGATA
               1140 -||||||||||||||||||||--|||||||||||||----|||-|--|--||-||||||-
query           964 TGCAACGCGAAGAACCTTACC-TGGTCTTGACATCC----ACAGAACTTTCC-AGAGAT-

target         1028 GGAC--G-T-CCCCTTCGGGGGCAGA--GTGA--CAGGTGG-TGCATGG-TTGTCGTCAG
               1200 |||---|-|-||--||||||---|-|--||||--||||||--|||||||-|-||||||||
query          1017 GGA-TTGGTGCC--TTCGGG---A-ACTGTGAGACAGGTG-CTGCATGGCT-GTCGTCAG

target         1078 CTCGTGTC-GTGAGA-TGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTGATCTTA--G
               1260 |||||||--||||-|-||||||||||||||||||||||||||||||||||-|||||---|
query          1068 CTCGTGT-TGTGA-AATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTT-ATCTT-TTG

target         1134 TTGCCAGCA--TTCA-GTTG--GGC-A-CTCTAA-GGT-GACTGCC-GGTGAC-AAACC-
               1320 ||||||||---|-|--|--|--||--|-|||-||-||--|||||||-|-|||--||||--
query          1124 TTGCCAGC-GGT-C-CG--GCCGG-GAACTC-AAAGG-AGACTGCCAG-TGA-TAAAC-T

target         1182 GGAGGAAGGTGGGGATGACGTCAAA-TCATCATG-CCCCTTAT-GACCT-GGGCTACACA
               1380 ||||||||||||||||||||||||--||||||||-|||-|||--||||--||||||||||
query          1173 GGAGGAAGGTGGGGATGACGTCAA-GTCATCATGGCCC-TTA-CGACC-AGGGCTACACA

target         1238 CGTGCTACAATGGACAG-A-ACAAAG-GGCA-GCGAAACC--GCGAG-GTT-AAGCC--A
               1440 |||||||||||||-|-|-|-||||||-|--|-||||--||--|||||-|---||||---|
query          1229 CGTGCTACAATGG-C-GCATACAAAGAG--AAGCGA--CCTCGCGAGAG--CAAGC-GGA

target         1288 ATCC-CAC-AAA-T-CTGTTC-TCAGTTC-GGATC-GC-AGTCTGCAACTCGACTGCG--
               1500 --||-||--|||-|-|-||-|-|-|||-|-||||--|--||||||||||||||||-|---
query          1280 --CCTCA-TAAAGTGC-GT-CGT-AGT-CCGGAT-TG-GAGTCTGCAACTCGACT-C-CA

target         1338 TGAAGCT-GGAATCGCTAGTAATCGC-GGATCAGCA-TGCCG-CGGTGAATACGTTCCCG
               1560 |||||-|-|||||||||||||||||--|||||||-|-||||--|||||||||||||||||
query          1329 TGAAG-TCGGAATCGCTAGTAATCG-TGGATCAG-AATGCC-ACGGTGAATACGTTCCCG

target         1394 GGCCTTGTACACACCGCCCGTCACACCAC-GAG-AGT---TTGT-AACACCC-GAAGTC-
               1620 ||||||||||||||||||||||||||||--|-|-|||---|||--||-|----|||||--
query          1385 GGCCTTGTACACACCGCCCGTCACACCA-TG-GGAGTGGGTTG-CAA-A---AGAAGT-A

target         1446 GGTGAGG-T-AACCTTTTA-GG-AG--C-C--AGCCG-CC---GAAGGTGGGA--CAGAT
               1680 |||-||--|-||||||----||-||--|-|--|-||--|----|----||--|--||--|
query          1437 GGT-AG-CTTAACCTT---CGGGAGGGCGCTTA-CC-AC-TTTG----TG--ATTCA--T

target         1491 GA-TTGGGGTGAAGTCGTAACAAGGTAG-CCGTATCGGAAGG----TGCGGCT-GGATCA
               1740 ||-|-||||||||||||||||||||||--|||||--||--||----|||||-|-||||||
query          1481 GACT-GGGGTGAAGTCGTAACAAGGTA-ACCGTA--GG--GGAACCTGCGG-TTGGATCA

target         1544 CCTCCTTTCTA 1555
               1800 |||||||---| 1811
query          1534 CCTCCTT---A 1542
""",
        )
        self.assertEqual(alignment.shape, (2, 1811))
        self.assertAlmostEqual(alignment.score, 1286.0)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1286.0; substitution score = 1286.0; gap score = 0.0; 1286 aligned letters; 1286 identities; 0 mismatches; 525 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1286.0:
        substitution_score = 1286.0,
        gap_score = 0.0.
    aligned = 1286:
        identities = 1286,
        mismatches = 0.
    gaps = 525:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 525:
            internal_insertions = 256:
                open_internal_insertions = 202,
                extend_internal_insertions = 54;
            internal_deletions = 269:
                open_internal_deletions = 199,
                extend_internal_deletions = 70;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 1286)
        self.assertEqual(counts.identities, 1286)
        self.assertEqual(counts.mismatches, 0)
        self.assertEqual(counts.deletions, 269)
        self.assertEqual(counts.insertions, 256)
        self.assertAlmostEqual(counts.score, 1286.0)
        alignments = aligner.align(seq1, reverse_complement(seq2), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (>{sys.maxsize} alignments; score=1286) at {hex(id(alignments))}>""",
        )
        self.assertAlmostEqual(alignments.score, 1286.0)
        message = "^number of optimal alignments is larger than (%d|%d)$" % (
            2147483647,  # on 32-bit systems
            9223372036854775807,
        )  # on 64-bit systems
        with self.assertRaisesRegex(OverflowError, message):
            n = len(alignments)
        # confirm that we can still pull out individual alignments
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
target            0 ATTTA-TC-GGA-GAGTTTGATCC-TGGCTCAGGAC--GAACGCTGGCGGC-GTGCCTAA
                  0 |---|-|--|-|-||||||||||--||||||||-|---|||||||||||||-|-||||||
query          1542 A---AAT-TG-AAGAGTTTGATC-ATGGCTCAG-A-TTGAACGCTGGCGGCAG-GCCTAA

target           53 T-ACATGCAAGTCGAG-CGG-A-CAG-AT-GGGA-GCTTGCT-C----CCTGAT-GTTAG
                 60 --|||||||||||||--|||-|-|||-|--|--|-|||||||-|----|-|||--|--||
query          1491 -CACATGCAAGTCGA-ACGGTAACAGGA-AG--AAGCTTGCTTCTTTGC-TGA-CG--AG

target          100 C-GGCGGACGGGTGAGTAACAC-GT--GGGTAA-CCTGCCTGTAA-G-ACTGGG--ATAA
                120 --|||||||||||||||||----||--|||-||-|-||||||-|--|-|--|||--||||
query          1440 -TGGCGGACGGGTGAGTAA---TGTCTGGG-AAAC-TGCCTG-A-TGGA--GGGGGATAA

target          151 CT-CC-GGGAAACCGG--GGCTAATACCGG-ATGGTTGTTTGAACCGCAT-GGTTCAA-A
                180 ||-|--||-||||-||--|-|||||||||--||---------|||-|--|-|---|||-|
query          1390 CTAC-TGG-AAAC-GGTAG-CTAATACCG-CAT---------AAC-G--TCG---CAAGA

target          204 C-ATAA-AAGGTGG--C-TTCGG-C-TACCACTTA-C-A--G-ATG-GACCC-GC--GGC
                240 |-|-||-|-||-||--|-|||||-|-|-|---||--|-|--|-|||-|-|||-|---||-
query          1350 CCA-AAGA-GG-GGGACCTTCGGGCCT-C---TT-GCCATCGGATGTG-CCCAG-ATGG-

target          248 GCATTAGCTAGTT-GGTGAGG-TAACGGCTCACC-AAGGCGACGATGCG--TAGCC-GA-
                300 |-||||||||||--||||-||-||||||||||||-|-|||||||||-|---||||--|--
query          1301 G-ATTAGCTAGT-AGGTG-GGGTAACGGCTCACCTA-GGCGACGAT-C-CCTAGC-TG-G

target          301 -CCTGAGAGGG-TGATC--GGCCACACTGGGA-CTGAGACACGG-CCCAGACTCCTACGG
                360 -|-|||||||--|||-|--|-|||||||||-|-|||||||||||-||-||||||||||||
query          1249 TC-TGAGAGG-ATGA-CCAG-CCACACTGG-AACTGAGACACGGTCC-AGACTCCTACGG

target          355 GAGGCAGCAGTAGGG-AATC-TTCCGCA-A-TGGA-CG-AAAGTC-TGAC-GG-AGCAAC
                420 |||||||||||-|||-|||--||--|||-|-|||--||-||-|-|-|||--|--|||--|
query          1195 GAGGCAGCAGT-GGGGAAT-ATT--GCACAATGG-GCGCAA-G-CCTGA-TG-CAGC--C

target          406 --GCCGCGTG-AGTGAT-GAAGG--TTTTCGGA-TC-GTAAAGCT-CTGTTGTT-AG-GG
                480 --||||||||-|-|||--|||||--||--|||--|--||||||-|-||-||----||-||
query          1146 ATGCCGCGTGTA-TGA-AGAAGGCCTT--CGG-GT-TGTAAAG-TACT-TT---CAGCGG

target          455 --A--A-G--A--ACAAGTGCCGTTCGAATAGGGC----GG-TACC-TTGACGGT-ACCT
                540 --|--|-|--|--|-||||----|---||||---|----|--|-|--||||||-|-|||-
query          1097 GGAGGAAGGGAGTA-AAGT----T---AATA---CCTTTG-CT-C-ATTGACG-TTACC-

target          499 AAC-CAGAA-A-GCCAC-GGCTAACTAC-GTGCCAGCAGCCGCGGTAATACGT-AGG-TG
                600 --|-|||||-|-||-||-||||||||-|-|||||||||||||||||||||||--|||-||
query          1053 --CGCAGAAGAAGC-ACCGGCTAACT-CCGTGCCAGCAGCCGCGGTAATACG-GAGGGTG

target          552 GCAAGCGTTG--TCCGGAATTA-TTGGGCGTAAAG-GGCT-CGCAGGCGGTTTC-TTAAG
                660 -||||||||---||-|||||||-|-||||||||||-|-|--||||||||||||--|||||
query           998 -CAAGCGTT-AATC-GGAATTACT-GGGCGTAAAGCG-C-ACGCAGGCGGTTT-GTTAAG

target          606 TCT-GATGTGAAAG-CCCCCGG-CTCAACC-GGGGAGGG--T-CAT-TGGA-AACTGGGG
                720 ||--|||||||||--||||-||-|||||||-|||-|-----|-|||-||-|-|-||||--
query           945 TC-AGATGTGAAA-TCCCC-GGGCTCAACCTGGG-A---ACTGCATCTG-ATA-CTGG--

target          657 -AA-CTTGAGTGCA--G-AAGAGGAGAGTGG-A-A-TTCCACG-TGTAGCGGTGAAATGC
                780 -||-|||||||-|---|-|-||||-|-|-||-|-|-|||||-|-||||||||||||||||
query           896 CAAGCTTGAGT-C-TCGTA-GAGG-G-G-GGTAGAATTCCA-GGTGTAGCGGTGAAATGC

target          708 GTAGAGATG-TGGAGGAAC-ACCAG-TGGCGAAGGCGA-CTCTC--TGGT-CTGTAA--C
                840 ||||||||--||||||||--|||-|-|||||||||||--|-|-|--|||--|-|-||--|
query           843 GTAGAGAT-CTGGAGGAA-TACC-GGTGGCGAAGGCG-GC-C-CCCTGG-AC-G-AAGAC

target          759 TGACGCTG-AGGA-GCGAAAGCGTGGGGAGCGAA-CAGGATTAGATACCCTGGTAGTCCA
                900 |||||||--|||--|||||||||||||||||-||-|||||||||||||||||||||||||
query           792 TGACGCT-CAGG-TGCGAAAGCGTGGGGAGC-AAACAGGATTAGATACCCTGGTAGTCCA

target          816 CGCCGTAAACGATGAGT-G-CTAAGTGTT-AGGGGGTT-TCCGCCCCTT-AGTGC-TG-C
                960 ||||||||||||||--|-|-||---||---|||---||-|--||||-||-||-||-||-|
query           735 CGCCGTAAACGATG--TCGACT---TG--GAGG---TTGT--GCCC-TTGAG-GCGTGGC

target          869 ------AGCTAACGCA-TTAAG-C-ACTCCGCCTGGGGAGTACGGTC-GCAAGACTG--A
               1020 ------|||||||||--|||||-|-||-|-|||||||||||||||-|-|||||---|--|
query           689 TTCCGGAGCTAACGC-GTTAAGTCGAC-C-GCCTGGGGAGTACGG-CCGCAAG---GTTA

target          917 AA-CTCAAA-GGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAA
               1080 ||-||||||-|-|||||||||||||||||||||||||||||||||||||||||||||||-
query           636 AAACTCAAATG-AATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGA-

target          975 -GCAACGCGAAGAACCTTACCA-GGTCTTGACATCCTCTGACA-A--T--CCTAGAGATA
               1140 -||||||||||||||||||||--|||||||||||||----|||-|--|--||-||||||-
query           578 TGCAACGCGAAGAACCTTACC-TGGTCTTGACATCC----ACAGAACTTTCC-AGAGAT-

target         1028 GGAC--G-T-CCCCTTCGGGGGCAGA--GTGA--CAGGTGG-TGCATGG-TTGTCGTCAG
               1200 |||---|-|-||--||||||---|-|--||||--||||||--|||||||-|-||||||||
query           525 GGA-TTGGTGCC--TTCGGG---A-ACTGTGAGACAGGTG-CTGCATGGCT-GTCGTCAG

target         1078 CTCGTGTC-GTGAGA-TGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTGATCTTA--G
               1260 |||||||--||||-|-||||||||||||||||||||||||||||||||||-|||||---|
query           474 CTCGTGT-TGTGA-AATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTT-ATCTT-TTG

target         1134 TTGCCAGCA--TTCA-GTTG--GGC-A-CTCTAA-GGT-GACTGCC-GGTGAC-AAACC-
               1320 ||||||||---|-|--|--|--||--|-|||-||-||--|||||||-|-|||--||||--
query           418 TTGCCAGC-GGT-C-CG--GCCGG-GAACTC-AAAGG-AGACTGCCAG-TGA-TAAAC-T

target         1182 GGAGGAAGGTGGGGATGACGTCAAA-TCATCATG-CCCCTTAT-GACCT-GGGCTACACA
               1380 ||||||||||||||||||||||||--||||||||-|||-|||--||||--||||||||||
query           369 GGAGGAAGGTGGGGATGACGTCAA-GTCATCATGGCCC-TTA-CGACC-AGGGCTACACA

target         1238 CGTGCTACAATGGACAG-A-ACAAAG-GGCA-GCGAAACC--GCGAG-GTT-AAGCC--A
               1440 |||||||||||||-|-|-|-||||||-|--|-||||--||--|||||-|---||||---|
query           313 CGTGCTACAATGG-C-GCATACAAAGAG--AAGCGA--CCTCGCGAGAG--CAAGC-GGA

target         1288 ATCC-CAC-AAA-T-CTGTTC-TCAGTTC-GGATC-GC-AGTCTGCAACTCGACTGCG--
               1500 --||-||--|||-|-|-||-|-|-|||-|-||||--|--||||||||||||||||-|---
query           262 --CCTCA-TAAAGTGC-GT-CGT-AGT-CCGGAT-TG-GAGTCTGCAACTCGACT-C-CA

target         1338 TGAAGCT-GGAATCGCTAGTAATCGC-GGATCAGCA-TGCCG-CGGTGAATACGTTCCCG
               1560 |||||-|-|||||||||||||||||--|||||||-|-||||--|||||||||||||||||
query           213 TGAAG-TCGGAATCGCTAGTAATCG-TGGATCAG-AATGCC-ACGGTGAATACGTTCCCG

target         1394 GGCCTTGTACACACCGCCCGTCACACCAC-GAG-AGT---TTGT-AACACCC-GAAGTC-
               1620 ||||||||||||||||||||||||||||--|-|-|||---|||--||-|----|||||--
query           157 GGCCTTGTACACACCGCCCGTCACACCA-TG-GGAGTGGGTTG-CAA-A---AGAAGT-A

target         1446 GGTGAGG-T-AACCTTTTA-GG-AG--C-C--AGCCG-CC---GAAGGTGGGA--CAGAT
               1680 |||-||--|-||||||----||-||--|-|--|-||--|----|----||--|--||--|
query           105 GGT-AG-CTTAACCTT---CGGGAGGGCGCTTA-CC-AC-TTTG----TG--ATTCA--T

target         1491 GA-TTGGGGTGAAGTCGTAACAAGGTAG-CCGTATCGGAAGG----TGCGGCT-GGATCA
               1740 ||-|-||||||||||||||||||||||--|||||--||--||----|||||-|-||||||
query            61 GACT-GGGGTGAAGTCGTAACAAGGTA-ACCGTA--GG--GGAACCTGCGG-TTGGATCA

target         1544 CCTCCTTTCTA 1555
               1800 |||||||---| 1811
query             8 CCTCCTT---A    0
""",
        )
        self.assertAlmostEqual(alignment.score, 1286.0)
        self.assertEqual(alignment.shape, (2, 1811))
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1286.0; substitution score = 1286.0; gap score = 0.0; 1286 aligned letters; 1286 identities; 0 mismatches; 525 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1286.0:
        substitution_score = 1286.0,
        gap_score = 0.0.
    aligned = 1286:
        identities = 1286,
        mismatches = 0.
    gaps = 525:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 525:
            internal_insertions = 256:
                open_internal_insertions = 202,
                extend_internal_insertions = 54;
            internal_deletions = 269:
                open_internal_deletions = 199,
                extend_internal_deletions = 70;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 1286)
        self.assertEqual(counts.identities, 1286)
        self.assertEqual(counts.mismatches, 0)
        self.assertEqual(counts.deletions, 269)
        self.assertEqual(counts.insertions, 256)
        self.assertAlmostEqual(counts.score, 1286.0)


class TestKeywordArgumentsConstructor(unittest.TestCase):
    def test_confusing_arguments(self):
        aligner = Align.PairwiseAligner(
            mode="local",
            open_gap_score=-0.3,
            extend_gap_score=-0.1,
            open_insertion_score=-0.2,
        )
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  wildcard: None
  match_score: 1.000000
  mismatch_score: 0.000000
  open_internal_insertion_score: -0.200000
  extend_internal_insertion_score: -0.100000
  open_left_insertion_score: -0.200000
  extend_left_insertion_score: -0.100000
  open_right_insertion_score: -0.200000
  extend_right_insertion_score: -0.100000
  open_internal_deletion_score: -0.300000
  extend_internal_deletion_score: -0.100000
  open_left_deletion_score: -0.300000
  extend_left_deletion_score: -0.100000
  open_right_deletion_score: -0.300000
  extend_right_deletion_score: -0.100000
  mode: local
""",
        )


class TestPredefinedScoringSchemes(unittest.TestCase):
    def test_blastn(self):
        aligner = Align.PairwiseAligner(scoring="blastn")
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  substitution_matrix: <Array object at 0x%x>
  open_internal_insertion_score: -7.000000
  extend_internal_insertion_score: -2.000000
  open_left_insertion_score: -7.000000
  extend_left_insertion_score: -2.000000
  open_right_insertion_score: -7.000000
  extend_right_insertion_score: -2.000000
  open_internal_deletion_score: -7.000000
  extend_internal_deletion_score: -2.000000
  open_left_deletion_score: -7.000000
  extend_left_deletion_score: -2.000000
  open_right_deletion_score: -7.000000
  extend_right_deletion_score: -2.000000
  mode: global
"""
            % id(aligner.substitution_matrix),
        )
        self.assertEqual(
            str(aligner.substitution_matrix[:, :]),
            """\
     A    T    G    C    S    W    R    Y    K    M    B    V    H    D    N
A  2.0 -3.0 -3.0 -3.0 -3.0 -1.0 -1.0 -3.0 -3.0 -1.0 -3.0 -1.0 -1.0 -1.0 -2.0
T -3.0  2.0 -3.0 -3.0 -3.0 -1.0 -3.0 -1.0 -1.0 -3.0 -1.0 -3.0 -1.0 -1.0 -2.0
G -3.0 -3.0  2.0 -3.0 -1.0 -3.0 -1.0 -3.0 -1.0 -3.0 -1.0 -1.0 -3.0 -1.0 -2.0
C -3.0 -3.0 -3.0  2.0 -1.0 -3.0 -3.0 -1.0 -3.0 -1.0 -1.0 -1.0 -1.0 -3.0 -2.0
S -3.0 -3.0 -1.0 -1.0 -1.0 -3.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -2.0
W -1.0 -1.0 -3.0 -3.0 -3.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -2.0
R -1.0 -3.0 -1.0 -3.0 -1.0 -1.0 -1.0 -3.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -2.0
Y -3.0 -1.0 -3.0 -1.0 -1.0 -1.0 -3.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -2.0
K -3.0 -1.0 -1.0 -3.0 -1.0 -1.0 -1.0 -1.0 -1.0 -3.0 -1.0 -1.0 -1.0 -1.0 -2.0
M -1.0 -3.0 -3.0 -1.0 -1.0 -1.0 -1.0 -1.0 -3.0 -1.0 -1.0 -1.0 -1.0 -1.0 -2.0
B -3.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -2.0
V -1.0 -3.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -2.0
H -1.0 -1.0 -3.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -2.0
D -1.0 -1.0 -1.0 -3.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -2.0
N -2.0 -2.0 -2.0 -2.0 -2.0 -2.0 -2.0 -2.0 -2.0 -2.0 -2.0 -2.0 -2.0 -2.0 -2.0
""",
        )

    def test_megablast(self):
        aligner = Align.PairwiseAligner(scoring="megablast")
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  substitution_matrix: <Array object at 0x%x>
  open_internal_insertion_score: -2.500000
  extend_internal_insertion_score: -2.500000
  open_left_insertion_score: -2.500000
  extend_left_insertion_score: -2.500000
  open_right_insertion_score: -2.500000
  extend_right_insertion_score: -2.500000
  open_internal_deletion_score: -2.500000
  extend_internal_deletion_score: -2.500000
  open_left_deletion_score: -2.500000
  extend_left_deletion_score: -2.500000
  open_right_deletion_score: -2.500000
  extend_right_deletion_score: -2.500000
  mode: global
"""
            % id(aligner.substitution_matrix),
        )
        self.assertEqual(
            str(aligner.substitution_matrix[:, :]),
            """\
     A    T    G    C    S    W    R    Y    K    M    B    V    H    D    N
A  1.0 -2.0 -2.0 -2.0 -2.0 -1.0 -1.0 -2.0 -2.0 -1.0 -2.0 -1.0 -1.0 -1.0 -1.0
T -2.0  1.0 -2.0 -2.0 -2.0 -1.0 -2.0 -1.0 -1.0 -2.0 -1.0 -2.0 -1.0 -1.0 -1.0
G -2.0 -2.0  1.0 -2.0 -1.0 -2.0 -1.0 -2.0 -1.0 -2.0 -1.0 -1.0 -2.0 -1.0 -1.0
C -2.0 -2.0 -2.0  1.0 -1.0 -2.0 -2.0 -1.0 -2.0 -1.0 -1.0 -1.0 -1.0 -2.0 -1.0
S -2.0 -2.0 -1.0 -1.0 -1.0 -2.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0
W -1.0 -1.0 -2.0 -2.0 -2.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0
R -1.0 -2.0 -1.0 -2.0 -1.0 -1.0 -1.0 -2.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0
Y -2.0 -1.0 -2.0 -1.0 -1.0 -1.0 -2.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0
K -2.0 -1.0 -1.0 -2.0 -1.0 -1.0 -1.0 -1.0 -1.0 -2.0 -1.0 -1.0 -1.0 -1.0 -1.0
M -1.0 -2.0 -2.0 -1.0 -1.0 -1.0 -1.0 -1.0 -2.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0
B -2.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0
V -1.0 -2.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0
H -1.0 -1.0 -2.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0
D -1.0 -1.0 -1.0 -2.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0
N -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0
""",
        )

    def test_blastp(self):
        aligner = Align.PairwiseAligner(scoring="blastp")
        self.assertEqual(
            str(aligner),
            """\
Pairwise sequence aligner with parameters
  substitution_matrix: <Array object at 0x%x>
  open_internal_insertion_score: -12.000000
  extend_internal_insertion_score: -1.000000
  open_left_insertion_score: -12.000000
  extend_left_insertion_score: -1.000000
  open_right_insertion_score: -12.000000
  extend_right_insertion_score: -1.000000
  open_internal_deletion_score: -12.000000
  extend_internal_deletion_score: -1.000000
  open_left_deletion_score: -12.000000
  extend_left_deletion_score: -1.000000
  open_right_deletion_score: -12.000000
  extend_right_deletion_score: -1.000000
  mode: global
"""
            % id(aligner.substitution_matrix),
        )
        self.assertEqual(
            str(aligner.substitution_matrix[:, :]),
            """\
     A    B    C    D    E    F    G    H    I    J    K    L    M    N    O    P    Q    R    S    T    U    V    W    X    Y    Z    *
A  4.0 -2.0  0.0 -2.0 -1.0 -2.0  0.0 -2.0 -1.0 -1.0 -1.0 -1.0 -1.0 -2.0 -1.0 -1.0 -1.0 -1.0  1.0  0.0  0.0  0.0 -3.0 -1.0 -2.0 -1.0 -4.0
B -2.0  4.0 -3.0  4.0  1.0 -3.0 -1.0  0.0 -3.0 -3.0  0.0 -4.0 -3.0  4.0 -1.0 -2.0  0.0 -1.0  0.0 -1.0 -3.0 -3.0 -4.0 -1.0 -3.0  0.0 -4.0
C  0.0 -3.0  9.0 -3.0 -4.0 -2.0 -3.0 -3.0 -1.0 -1.0 -3.0 -1.0 -1.0 -3.0 -1.0 -3.0 -3.0 -3.0 -1.0 -1.0  9.0 -1.0 -2.0 -1.0 -2.0 -3.0 -4.0
D -2.0  4.0 -3.0  6.0  2.0 -3.0 -1.0 -1.0 -3.0 -3.0 -1.0 -4.0 -3.0  1.0 -1.0 -1.0  0.0 -2.0  0.0 -1.0 -3.0 -3.0 -4.0 -1.0 -3.0  1.0 -4.0
E -1.0  1.0 -4.0  2.0  5.0 -3.0 -2.0  0.0 -3.0 -3.0  1.0 -3.0 -2.0  0.0 -1.0 -1.0  2.0  0.0  0.0 -1.0 -4.0 -2.0 -3.0 -1.0 -2.0  4.0 -4.0
F -2.0 -3.0 -2.0 -3.0 -3.0  6.0 -3.0 -1.0  0.0  0.0 -3.0  0.0  0.0 -3.0 -1.0 -4.0 -3.0 -3.0 -2.0 -2.0 -2.0 -1.0  1.0 -1.0  3.0 -3.0 -4.0
G  0.0 -1.0 -3.0 -1.0 -2.0 -3.0  6.0 -2.0 -4.0 -4.0 -2.0 -4.0 -3.0  0.0 -1.0 -2.0 -2.0 -2.0  0.0 -2.0 -3.0 -3.0 -2.0 -1.0 -3.0 -2.0 -4.0
H -2.0  0.0 -3.0 -1.0  0.0 -1.0 -2.0  8.0 -3.0 -3.0 -1.0 -3.0 -2.0  1.0 -1.0 -2.0  0.0  0.0 -1.0 -2.0 -3.0 -3.0 -2.0 -1.0  2.0  0.0 -4.0
I -1.0 -3.0 -1.0 -3.0 -3.0  0.0 -4.0 -3.0  4.0  3.0 -3.0  2.0  1.0 -3.0 -1.0 -3.0 -3.0 -3.0 -2.0 -1.0 -1.0  3.0 -3.0 -1.0 -1.0 -3.0 -4.0
J -1.0 -3.0 -1.0 -3.0 -3.0  0.0 -4.0 -3.0  3.0  3.0 -3.0  3.0  2.0 -3.0 -1.0 -3.0 -2.0 -2.0 -2.0 -1.0 -1.0  2.0 -2.0 -1.0 -1.0 -3.0 -4.0
K -1.0  0.0 -3.0 -1.0  1.0 -3.0 -2.0 -1.0 -3.0 -3.0  5.0 -2.0 -1.0  0.0 -1.0 -1.0  1.0  2.0  0.0 -1.0 -3.0 -2.0 -3.0 -1.0 -2.0  1.0 -4.0
L -1.0 -4.0 -1.0 -4.0 -3.0  0.0 -4.0 -3.0  2.0  3.0 -2.0  4.0  2.0 -3.0 -1.0 -3.0 -2.0 -2.0 -2.0 -1.0 -1.0  1.0 -2.0 -1.0 -1.0 -3.0 -4.0
M -1.0 -3.0 -1.0 -3.0 -2.0  0.0 -3.0 -2.0  1.0  2.0 -1.0  2.0  5.0 -2.0 -1.0 -2.0  0.0 -1.0 -1.0 -1.0 -1.0  1.0 -1.0 -1.0 -1.0 -1.0 -4.0
N -2.0  4.0 -3.0  1.0  0.0 -3.0  0.0  1.0 -3.0 -3.0  0.0 -3.0 -2.0  6.0 -1.0 -2.0  0.0  0.0  1.0  0.0 -3.0 -3.0 -4.0 -1.0 -2.0  0.0 -4.0
O -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -4.0
P -1.0 -2.0 -3.0 -1.0 -1.0 -4.0 -2.0 -2.0 -3.0 -3.0 -1.0 -3.0 -2.0 -2.0 -1.0  7.0 -1.0 -2.0 -1.0 -1.0 -3.0 -2.0 -4.0 -1.0 -3.0 -1.0 -4.0
Q -1.0  0.0 -3.0  0.0  2.0 -3.0 -2.0  0.0 -3.0 -2.0  1.0 -2.0  0.0  0.0 -1.0 -1.0  5.0  1.0  0.0 -1.0 -3.0 -2.0 -2.0 -1.0 -1.0  4.0 -4.0
R -1.0 -1.0 -3.0 -2.0  0.0 -3.0 -2.0  0.0 -3.0 -2.0  2.0 -2.0 -1.0  0.0 -1.0 -2.0  1.0  5.0 -1.0 -1.0 -3.0 -3.0 -3.0 -1.0 -2.0  0.0 -4.0
S  1.0  0.0 -1.0  0.0  0.0 -2.0  0.0 -1.0 -2.0 -2.0  0.0 -2.0 -1.0  1.0 -1.0 -1.0  0.0 -1.0  4.0  1.0 -1.0 -2.0 -3.0 -1.0 -2.0  0.0 -4.0
T  0.0 -1.0 -1.0 -1.0 -1.0 -2.0 -2.0 -2.0 -1.0 -1.0 -1.0 -1.0 -1.0  0.0 -1.0 -1.0 -1.0 -1.0  1.0  5.0 -1.0  0.0 -2.0 -1.0 -2.0 -1.0 -4.0
U  0.0 -3.0  9.0 -3.0 -4.0 -2.0 -3.0 -3.0 -1.0 -1.0 -3.0 -1.0 -1.0 -3.0 -1.0 -3.0 -3.0 -3.0 -1.0 -1.0  9.0 -1.0 -2.0 -1.0 -2.0 -3.0 -4.0
V  0.0 -3.0 -1.0 -3.0 -2.0 -1.0 -3.0 -3.0  3.0  2.0 -2.0  1.0  1.0 -3.0 -1.0 -2.0 -2.0 -3.0 -2.0  0.0 -1.0  4.0 -3.0 -1.0 -1.0 -2.0 -4.0
W -3.0 -4.0 -2.0 -4.0 -3.0  1.0 -2.0 -2.0 -3.0 -2.0 -3.0 -2.0 -1.0 -4.0 -1.0 -4.0 -2.0 -3.0 -3.0 -2.0 -2.0 -3.0 11.0 -1.0  2.0 -2.0 -4.0
X -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -4.0
Y -2.0 -3.0 -2.0 -3.0 -2.0  3.0 -3.0  2.0 -1.0 -1.0 -2.0 -1.0 -1.0 -2.0 -1.0 -3.0 -1.0 -2.0 -2.0 -2.0 -2.0 -1.0  2.0 -1.0  7.0 -2.0 -4.0
Z -1.0  0.0 -3.0  1.0  4.0 -3.0 -2.0  0.0 -3.0 -3.0  1.0 -3.0 -1.0  0.0 -1.0 -1.0  4.0  0.0  0.0 -1.0 -3.0 -2.0 -2.0 -1.0 -2.0  4.0 -4.0
* -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0 -4.0  1.0
""",
        )


class TestUnicodeStrings(unittest.TestCase):
    def test_needlemanwunsch_simple1(self):
        seq1 = "ĞĀĀČŦ"
        seq2 = "ĞĀŦ"
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        self.assertEqual(aligner.algorithm, "Needleman-Wunsch")
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 1.0)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (2 alignments; score=1) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1.0)
        self.assertEqual(
            str(alignment),
            """\
ĞĀĀČŦ
||--|
ĞĀ--Ŧ
""",
        )
        self.assertEqual(alignment.shape, (2, 5))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 2], [4, 5]], [[0, 2], [2, 3]]])
            )
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.0; substitution score = 3.0; gap score = -2.0; 3 aligned letters; 3 identities; 0 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.0:
        substitution_score = 3.0,
        gap_score = -2.0.
    aligned = 3:
        identities = 3,
        mismatches = 0.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertEqual(counts.deletions, 2)
        self.assertEqual(counts.insertions, 0)
        self.assertAlmostEqual(counts.score, 1.0)
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 1.0)
        self.assertEqual(
            str(alignment),
            """\
ĞĀĀČŦ
|-|-|
Ğ-Ā-Ŧ
""",
        )
        self.assertEqual(alignment.shape, (2, 5))
        self.assertTrue(
            np.array_equal(
                alignment.aligned,
                np.array([[[0, 1], [2, 3], [4, 5]], [[0, 1], [1, 2], [2, 3]]]),
            )
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.0; substitution score = 3.0; gap score = -2.0; 3 aligned letters; 3 identities; 0 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.0:
        substitution_score = 3.0,
        gap_score = -2.0.
    aligned = 3:
        identities = 3,
        mismatches = 0.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 2,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertEqual(counts.deletions, 2)
        self.assertEqual(counts.insertions, 0)
        self.assertAlmostEqual(counts.score, 1.0)

    def test_needlemanwunsch_simple1_fogsaa(self):
        seq1 = "ĞĀĀČŦ"
        seq2 = "ĞĀŦ"
        aligner = Align.PairwiseAligner()
        aligner.mode = "fogsaa"
        self.assertEqual(
            aligner.algorithm, "Fast Optimal Global Sequence Alignment Algorithm"
        )
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 1.0)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=1) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1.0)
        self.assertEqual(
            str(alignment),
            """\
ĞĀĀČŦ
||--|
ĞĀ--Ŧ
""",
        )
        self.assertEqual(alignment.shape, (2, 5))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 2], [4, 5]], [[0, 2], [2, 3]]])
            )
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.0; substitution score = 3.0; gap score = -2.0; 3 aligned letters; 3 identities; 0 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.0:
        substitution_score = 3.0,
        gap_score = -2.0.
    aligned = 3:
        identities = 3,
        mismatches = 0.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 2:
                open_internal_deletions = 1,
                extend_internal_deletions = 1;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertEqual(counts.deletions, 2)
        self.assertEqual(counts.insertions, 0)
        self.assertAlmostEqual(counts.score, 1.0)

    def test_align_affine1_score(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score = 0
        aligner.mismatch_score = -1
        aligner.open_gap_score = -5
        aligner.extend_gap_score = -1
        self.assertEqual(aligner.algorithm, "Gotoh global alignment algorithm")
        score = aligner.score("いい", "あいいう")
        self.assertAlmostEqual(score, -7.0)
        alignments = aligner.align("いい", "あいいう")
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
--いい
--|.
あいいう
""",
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = -7.0; substitution score = -1.0; gap score = -6.0; 2 aligned letters; 1 identities; 1 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = -7.0:
        substitution_score = -1.0,
        gap_score = -6.0.
    aligned = 2:
        identities = 1,
        mismatches = 1.
    gaps = 2:
        left_gaps = 2:
            left_insertions = 2:
                open_left_insertions = 1,
                extend_left_insertions = 1;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 1)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.insertions, 2)
        self.assertAlmostEqual(counts.score, -7.0)
        alignment = alignments[1]
        self.assertEqual(
            str(alignment),
            """\
いい--
.|--
あいいう
""",
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = -7.0; substitution score = -1.0; gap score = -6.0; 2 aligned letters; 1 identities; 1 mismatches; 2 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = -7.0:
        substitution_score = -1.0,
        gap_score = -6.0.
    aligned = 2:
        identities = 1,
        mismatches = 1.
    gaps = 2:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 2:
            right_insertions = 2:
                open_right_insertions = 1,
                extend_right_insertions = 1;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 1)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.insertions, 2)
        self.assertAlmostEqual(counts.score, -7.0)

    def test_align_affine1_score_fogsaa(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = "fogsaa"
        aligner.match_score = 0
        aligner.mismatch_score = -1
        aligner.open_gap_score = -5
        aligner.extend_gap_score = -1
        self.assertEqual(
            aligner.algorithm, "Fast Optimal Global Sequence Alignment Algorithm"
        )
        score = aligner.score("いい", "あいいう")
        self.assertAlmostEqual(score, -7.0)

    def test_smithwaterman(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = "local"
        aligner.gap_score = -0.1
        self.assertEqual(aligner.algorithm, "Smith-Waterman")
        score = aligner.score("ℵℷℶℷ", "ℸℵℶℸ")
        self.assertAlmostEqual(score, 1.9)
        alignments = aligner.align("ℵℷℶℷ", "ℸℵℶℸ")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=1.9) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1.9)
        self.assertEqual(
            str(alignment),
            """\
ℵℷℶ
|-|
ℵ-ℶ
""",
        )
        self.assertEqual(alignment.shape, (2, 3))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 1], [2, 3]], [[1, 2], [2, 3]]])
            )
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.9; substitution score = 2.0; gap score = -0.1; 2 aligned letters; 2 identities; 0 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.9:
        substitution_score = 2.0,
        gap_score = -0.1.
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertEqual(counts.deletions, 1)
        self.assertEqual(counts.insertions, 0)
        self.assertAlmostEqual(counts.score, 1.9)

    def test_gotoh_local(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = "local"
        aligner.open_gap_score = -0.1
        aligner.extend_gap_score = 0.0
        self.assertEqual(aligner.algorithm, "Gotoh local alignment algorithm")
        score = aligner.score("生物科物", "学生科学")
        self.assertAlmostEqual(score, 1.9)
        alignments = aligner.align("生物科物", "学生科学")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=1.9) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1.9)
        self.assertEqual(
            str(alignment),
            """\
生物科
|-|
生-科
""",
        )
        self.assertEqual(alignment.shape, (2, 3))
        self.assertTrue(
            np.array_equal(
                alignment.aligned, np.array([[[0, 1], [2, 3]], [[1, 2], [2, 3]]])
            )
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1.9; substitution score = 2.0; gap score = -0.1; 2 aligned letters; 2 identities; 0 mismatches; 1 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1.9:
        substitution_score = 2.0,
        gap_score = -0.1.
    aligned = 2:
        identities = 2,
        mismatches = 0.
    gaps = 1:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 2)
        self.assertEqual(counts.identities, 2)
        self.assertEqual(counts.mismatches, 0)
        self.assertEqual(counts.deletions, 1)
        self.assertEqual(counts.insertions, 0)
        self.assertAlmostEqual(counts.score, 1.9)


class TestAlignerPickling(unittest.TestCase):
    def test_pickle_aligner_match_mismatch(self):
        import pickle

        aligner = Align.PairwiseAligner()
        aligner.wildcard = "X"
        aligner.match_score = 3
        aligner.mismatch_score = -2
        aligner.open_internal_insertion_score = -2.5
        aligner.extend_internal_insertion_score = -3.5
        aligner.open_left_insertion_score = -2.5
        aligner.extend_left_insertion_score = -3.5
        aligner.open_right_insertion_score = -4
        aligner.extend_right_insertion_score = -4
        aligner.open_internal_deletion_score = -0.1
        aligner.extend_internal_deletion_score = -2
        aligner.open_left_deletion_score = -9
        aligner.extend_left_deletion_score = +1
        aligner.open_right_deletion_score = -1
        aligner.extend_right_deletion_score = -2
        aligner.mode = "local"
        state = pickle.dumps(aligner)
        pickled_aligner = pickle.loads(state)
        self.assertEqual(aligner.wildcard, pickled_aligner.wildcard)
        self.assertAlmostEqual(aligner.match_score, pickled_aligner.match_score)
        self.assertAlmostEqual(aligner.mismatch_score, pickled_aligner.mismatch_score)
        self.assertIsNone(pickled_aligner.substitution_matrix)
        self.assertAlmostEqual(
            aligner.open_internal_insertion_score,
            pickled_aligner.open_internal_insertion_score,
        )
        self.assertAlmostEqual(
            aligner.extend_internal_insertion_score,
            pickled_aligner.extend_internal_insertion_score,
        )
        self.assertAlmostEqual(
            aligner.open_left_insertion_score,
            pickled_aligner.open_left_insertion_score,
        )
        self.assertAlmostEqual(
            aligner.extend_left_insertion_score,
            pickled_aligner.extend_left_insertion_score,
        )
        self.assertAlmostEqual(
            aligner.open_right_insertion_score,
            pickled_aligner.open_right_insertion_score,
        )
        self.assertAlmostEqual(
            aligner.extend_right_insertion_score,
            pickled_aligner.extend_right_insertion_score,
        )
        self.assertAlmostEqual(
            aligner.open_internal_deletion_score,
            pickled_aligner.open_internal_deletion_score,
        )
        self.assertAlmostEqual(
            aligner.extend_internal_deletion_score,
            pickled_aligner.extend_internal_deletion_score,
        )
        self.assertAlmostEqual(
            aligner.open_left_deletion_score, pickled_aligner.open_left_deletion_score
        )
        self.assertAlmostEqual(
            aligner.extend_left_deletion_score,
            pickled_aligner.extend_left_deletion_score,
        )
        self.assertAlmostEqual(
            aligner.open_right_deletion_score,
            pickled_aligner.open_right_deletion_score,
        )
        self.assertAlmostEqual(
            aligner.extend_right_deletion_score,
            pickled_aligner.extend_right_deletion_score,
        )
        self.assertEqual(aligner.mode, pickled_aligner.mode)

    def test_pickle_aligner_substitution_matrix(self):
        try:
            from Bio.Align import substitution_matrices
        except ImportError:
            return
        import pickle

        aligner = Align.PairwiseAligner()
        aligner.wildcard = "N"
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM80")
        aligner.open_internal_insertion_score = -5
        aligner.extend_internal_insertion_score = -3
        aligner.open_left_insertion_score = -2
        aligner.extend_left_insertion_score = -3
        aligner.open_right_insertion_score = -4.5
        aligner.extend_right_insertion_score = -4.3
        aligner.open_internal_deletion_score = -2
        aligner.extend_internal_deletion_score = -2.5
        aligner.open_left_deletion_score = -9.1
        aligner.extend_left_deletion_score = +1.7
        aligner.open_right_deletion_score = -1.9
        aligner.extend_right_deletion_score = -2.0
        aligner.mode = "global"
        state = pickle.dumps(aligner)
        pickled_aligner = pickle.loads(state)
        self.assertEqual(aligner.wildcard, pickled_aligner.wildcard)
        self.assertIsNone(pickled_aligner.match_score)
        self.assertIsNone(pickled_aligner.mismatch_score)
        self.assertTrue(
            (aligner.substitution_matrix == pickled_aligner.substitution_matrix).all()
        )
        self.assertEqual(
            aligner.substitution_matrix.alphabet,
            pickled_aligner.substitution_matrix.alphabet,
        )
        self.assertAlmostEqual(
            aligner.open_internal_insertion_score,
            pickled_aligner.open_internal_insertion_score,
        )
        self.assertAlmostEqual(
            aligner.extend_internal_insertion_score,
            pickled_aligner.extend_internal_insertion_score,
        )
        self.assertAlmostEqual(
            aligner.open_left_insertion_score,
            pickled_aligner.open_left_insertion_score,
        )
        self.assertAlmostEqual(
            aligner.extend_left_insertion_score,
            pickled_aligner.extend_left_insertion_score,
        )
        self.assertAlmostEqual(
            aligner.open_right_insertion_score,
            pickled_aligner.open_right_insertion_score,
        )
        self.assertAlmostEqual(
            aligner.extend_right_insertion_score,
            pickled_aligner.extend_right_insertion_score,
        )
        self.assertAlmostEqual(
            aligner.open_internal_deletion_score,
            pickled_aligner.open_internal_deletion_score,
        )
        self.assertAlmostEqual(
            aligner.extend_internal_deletion_score,
            pickled_aligner.extend_internal_deletion_score,
        )
        self.assertAlmostEqual(
            aligner.open_left_deletion_score, pickled_aligner.open_left_deletion_score
        )
        self.assertAlmostEqual(
            aligner.extend_left_deletion_score,
            pickled_aligner.extend_left_deletion_score,
        )
        self.assertAlmostEqual(
            aligner.open_right_deletion_score,
            pickled_aligner.open_right_deletion_score,
        )
        self.assertAlmostEqual(
            aligner.extend_right_deletion_score,
            pickled_aligner.extend_right_deletion_score,
        )
        self.assertEqual(aligner.mode, pickled_aligner.mode)


class TestAlignmentFormat(unittest.TestCase):
    def test_alignment_simple(self):
        chromosome = "ACGATCAGCGAGCATNGAGCACTACGACAGCGAGTGACCACTATTCGCGATCAGGAGCAGATACTTTACGAGCATCGGC"
        transcript = "AGCATCGAGCGACTTGAGTACTATTCATACTTTCGAGC"
        aligner = Align.PairwiseAligner()
        aligner.extend_deletion_score = 0
        aligner.open_deletion_score = -3
        aligner.insertion_score = -3
        aligner.end_gap_score = 0
        aligner.mismatch = -1
        alignments = aligner.align(chromosome, transcript)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=19) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 19.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 ACGATCAGCGAGCATNGAGC-ACTACGACAGCGAGTGACCACTATTCGCGATCAGGAGCA
                  0 ----------|||||.||||-|||-----------|||..|||||||-------------
query             0 ----------AGCATCGAGCGACT-----------TGAGTACTATTC-------------

target           59 GATACTTTACGAGCATCGGC 79
                 60 -|||||||-|||||------ 80
query            26 -ATACTTT-CGAGC------ 38
""",
        )
        self.assertEqual(alignment.shape, (2, 80))
        self.assertEqual(
            alignment.format("psl"),
            """\
34	2	0	1	1	1	3	26	+	query	38	0	38	target	79	10	73	5	10,3,12,7,5,	0,11,14,26,33,	10,20,34,60,68,
""",
        )
        self.assertEqual(
            alignment.format("bed"),
            """\
target	10	73	query	19	+	10	73	0	5	10,3,12,7,5,	0,10,24,50,58,
""",
        )
        self.assertEqual(
            alignment.format("sam"),
            """\
query	0	target	1	255	10D10M1I3M11D12M14D7M1D5M6D	*	0	0	AGCATCGAGCGACTTGAGTACTATTCATACTTTCGAGC	*	AS:i:19
""",
        )
        self.assertAlmostEqual(alignment.counts(aligner).score, 19.0)
        alignments = aligner.align(
            chromosome, reverse_complement(transcript), strand="-"
        )
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=19) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 19.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 ACGATCAGCGAGCATNGAGC-ACTACGACAGCGAGTGACCACTATTCGCGATCAGGAGCA
                  0 ----------|||||.||||-|||-----------|||..|||||||-------------
query            38 ----------AGCATCGAGCGACT-----------TGAGTACTATTC-------------

target           59 GATACTTTACGAGCATCGGC 79
                 60 -|||||||-|||||------ 80
query            12 -ATACTTT-CGAGC------  0
""",
        )
        self.assertEqual(alignment.shape, (2, 80))
        self.assertEqual(
            alignment.format("psl"),
            """\
34	2	0	1	1	1	3	26	-	query	38	0	38	target	79	10	73	5	10,3,12,7,5,	0,11,14,26,33,	10,20,34,60,68,
""",
        )
        self.assertEqual(
            alignment.format("bed"),
            """\
target	10	73	query	19	-	10	73	0	5	10,3,12,7,5,	0,10,24,50,58,
""",
        )
        self.assertEqual(
            alignment.format("sam"),
            """\
query	16	target	1	255	10D10M1I3M11D12M14D7M1D5M6D	*	0	0	AGCATCGAGCGACTTGAGTACTATTCATACTTTCGAGC	*	AS:i:19
""",
        )
        self.assertAlmostEqual(alignment.counts(aligner).score, 19.0)

    def test_alignment_end_gap(self):
        aligner = Align.PairwiseAligner()
        aligner.gap_score = -1
        aligner.end_gap_score = 0
        aligner.mismatch = -10
        alignments = aligner.align("ACGTAGCATCAGC", "CCCCACGTAGCATCAGC")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=13) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        self.assertAlmostEqual(alignments.score, 13.0)
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
target            0 ----ACGTAGCATCAGC 13
                  0 ----||||||||||||| 17
query             0 CCCCACGTAGCATCAGC 17
""",
        )
        self.assertEqual(alignment.shape, (2, 17))
        self.assertEqual(
            alignment.format("psl"),
            """\
13	0	0	0	0	0	0	0	+	query	17	4	17	target	13	0	13	1	13,	4,	0,
""",
        )
        self.assertEqual(
            alignment.format("bed"),
            """\
target	0	13	query	13	+	0	13	0	1	13,	0,
""",
        )
        self.assertEqual(
            alignment.format("sam"),
            """\
query	0	target	1	255	4I13M	*	0	0	CCCCACGTAGCATCAGC	*	AS:i:13
""",
        )
        self.assertAlmostEqual(alignment.counts(aligner).score, 13.0)
        alignments = aligner.align(
            "ACGTAGCATCAGC", reverse_complement("CCCCACGTAGCATCAGC"), strand="-"
        )
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=13) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        self.assertAlmostEqual(alignments.score, 13.0)
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
target            0 ----ACGTAGCATCAGC 13
                  0 ----||||||||||||| 17
query            17 CCCCACGTAGCATCAGC  0
""",
        )
        self.assertEqual(alignment.shape, (2, 17))
        self.assertEqual(
            alignment.format("psl"),
            """\
13	0	0	0	0	0	0	0	-	query	17	0	13	target	13	0	13	1	13,	4,	0,
""",
        )
        self.assertEqual(
            alignment.format("bed"),
            """\
target	0	13	query	13	-	0	13	0	1	13,	0,
""",
        )
        self.assertEqual(
            alignment.format("sam"),
            """\
query	16	target	1	255	4I13M	*	0	0	CCCCACGTAGCATCAGC	*	AS:i:13
""",
        )
        self.assertAlmostEqual(alignment.counts(aligner).score, 13.0)
        alignments = aligner.align("CCCCACGTAGCATCAGC", "ACGTAGCATCAGC")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=13) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 13.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 CCCCACGTAGCATCAGC 17
                  0 ----||||||||||||| 17
query             0 ----ACGTAGCATCAGC 13
""",
        )
        self.assertEqual(alignment.shape, (2, 17))
        self.assertEqual(
            alignment.format("psl"),
            """\
13	0	0	0	0	0	0	0	+	query	13	0	13	target	17	4	17	1	13,	0,	4,
""",
        )
        self.assertEqual(
            alignment.format("bed"),
            """\
target	4	17	query	13	+	4	17	0	1	13,	0,
""",
        )
        self.assertEqual(
            alignment.format("sam"),
            """\
query	0	target	1	255	4D13M	*	0	0	ACGTAGCATCAGC	*	AS:i:13
""",
        )
        self.assertAlmostEqual(alignment.counts(aligner).score, 13.0)
        alignments = aligner.align(
            "CCCCACGTAGCATCAGC", reverse_complement("ACGTAGCATCAGC"), strand="-"
        )
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=13) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 13.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 CCCCACGTAGCATCAGC 17
                  0 ----||||||||||||| 17
query            13 ----ACGTAGCATCAGC  0
""",
        )
        self.assertEqual(alignment.shape, (2, 17))
        self.assertEqual(
            alignment.format("psl"),
            """\
13	0	0	0	0	0	0	0	-	query	13	0	13	target	17	4	17	1	13,	0,	4,
""",
        )
        self.assertEqual(
            alignment.format("bed"),
            """\
target	4	17	query	13	-	4	17	0	1	13,	0,
""",
        )
        self.assertEqual(
            alignment.format("sam"),
            """\
query	16	target	1	255	4D13M	*	0	0	ACGTAGCATCAGC	*	AS:i:13
""",
        )
        self.assertAlmostEqual(alignment.counts(aligner).score, 13.0)
        alignments = aligner.align("ACGTAGCATCAGC", "ACGTAGCATCAGCGGGG")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=13) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
target            0 ACGTAGCATCAGC---- 13
                  0 |||||||||||||---- 17
query             0 ACGTAGCATCAGCGGGG 17
""",
        )
        self.assertEqual(alignment.shape, (2, 17))
        self.assertEqual(
            alignment.format("psl"),
            """\
13	0	0	0	0	0	0	0	+	query	17	0	13	target	13	0	13	1	13,	0,	0,
""",
        )
        self.assertEqual(
            alignment.format("bed"),
            """\
target	0	13	query	13	+	0	13	0	1	13,	0,
""",
        )
        self.assertEqual(
            alignment.format("sam"),
            """\
query	0	target	1	255	13M4I	*	0	0	ACGTAGCATCAGCGGGG	*	AS:i:13
""",
        )
        self.assertAlmostEqual(alignment.counts(aligner).score, 13.0)
        alignments = aligner.align(
            "ACGTAGCATCAGC", reverse_complement("ACGTAGCATCAGCGGGG"), strand="-"
        )
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=13) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
target            0 ACGTAGCATCAGC---- 13
                  0 |||||||||||||---- 17
query            17 ACGTAGCATCAGCGGGG  0
""",
        )
        self.assertEqual(alignment.shape, (2, 17))
        self.assertEqual(
            alignment.format("psl"),
            """\
13	0	0	0	0	0	0	0	-	query	17	4	17	target	13	0	13	1	13,	0,	0,
""",
        )
        self.assertEqual(
            alignment.format("bed"),
            """\
target	0	13	query	13	-	0	13	0	1	13,	0,
""",
        )
        self.assertEqual(
            alignment.format("sam"),
            """\
query	16	target	1	255	13M4I	*	0	0	ACGTAGCATCAGCGGGG	*	AS:i:13
""",
        )
        self.assertAlmostEqual(alignment.counts(aligner).score, 13.0)
        alignments = aligner.align("ACGTAGCATCAGCGGGG", "ACGTAGCATCAGC")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=13) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 13.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 ACGTAGCATCAGCGGGG 17
                  0 |||||||||||||---- 17
query             0 ACGTAGCATCAGC---- 13
""",
        )
        self.assertEqual(alignment.shape, (2, 17))
        self.assertEqual(
            alignment.format("psl"),
            """\
13	0	0	0	0	0	0	0	+	query	13	0	13	target	17	0	13	1	13,	0,	0,
""",
        )
        self.assertEqual(
            alignment.format("bed"),
            """\
target	0	13	query	13	+	0	13	0	1	13,	0,
""",
        )
        self.assertEqual(
            alignment.format("sam"),
            """\
query	0	target	1	255	13M4D	*	0	0	ACGTAGCATCAGC	*	AS:i:13
""",
        )
        self.assertAlmostEqual(alignment.counts(aligner).score, 13.0)
        alignments = aligner.align(
            "ACGTAGCATCAGCGGGG", reverse_complement("ACGTAGCATCAGC"), strand="-"
        )
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=13) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 13.0)
        self.assertEqual(
            str(alignment),
            """\
target            0 ACGTAGCATCAGCGGGG 17
                  0 |||||||||||||---- 17
query            13 ACGTAGCATCAGC----  0
""",
        )
        self.assertEqual(alignment.shape, (2, 17))
        self.assertEqual(
            alignment.format("psl"),
            """\
13	0	0	0	0	0	0	0	-	query	13	0	13	target	17	0	13	1	13,	0,	0,
""",
        )
        self.assertEqual(
            alignment.format("bed"),
            """\
target	0	13	query	13	-	0	13	0	1	13,	0,
""",
        )
        self.assertEqual(
            alignment.format("sam"),
            """\
query	16	target	1	255	13M4D	*	0	0	ACGTAGCATCAGC	*	AS:i:13
""",
        )
        self.assertAlmostEqual(alignment.counts(aligner).score, 13.0)

    def test_alignment_wildcard(self):
        aligner = Align.PairwiseAligner()
        aligner.gap_score = -10
        aligner.mismatch = -2
        aligner.wildcard = "N"
        target = "TTTTTNACGCTCGAGCAGCTACG"
        query = "ACGATCGAGCNGCTACGCCCNC"
        # local alignment
        aligner.mode = "local"
        # use strings for target and query
        alignments = aligner.align(target, query)
        self.assertAlmostEqual(alignments.score, 13.0)
        self.assertEqual(len(alignments), 1)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=13) at {hex(id(alignments))}>""",
        )
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
target            6 ACGCTCGAGCAGCTACG 23
                  0 |||.||||||.|||||| 17
query             0 ACGATCGAGCNGCTACG 17
""",
        )
        self.assertEqual(alignment.shape, (2, 17))
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (17 aligned letters; 15 identities; 2 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 17:
        identities = 15,
        mismatches = 2.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 17)
        self.assertEqual(counts.identities, 15)
        self.assertEqual(counts.mismatches, 2)
        counts = alignment.counts("N")
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (17 aligned letters; 15 identities; 1 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 17:
        identities = 15,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 17)
        self.assertEqual(counts.identities, 15)
        self.assertEqual(counts.mismatches, 1)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 13.0; substitution score = 13.0; gap score = 0.0; 17 aligned letters; 15 identities; 1 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 13.0:
        substitution_score = 13.0,
        gap_score = 0.0.
    aligned = 17:
        identities = 15,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 17)
        self.assertEqual(counts.identities, 15)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 13.0)
        self.assertEqual(
            alignment.format("psl"),
            """\
15	1	0	1	0	0	0	0	+	query	22	0	17	target	23	6	23	1	17,	0,	6,
""",
        )
        self.assertEqual(
            alignment.format("bed"),
            """\
target	6	23	query	13	+	6	23	0	1	17,	0,
""",
        )
        self.assertEqual(
            alignment.format("sam"),
            """\
query	0	target	7	255	17M5S	*	0	0	ACGATCGAGCNGCTACGCCCNC	*	AS:i:13
""",
        )
        alignments = aligner.align(target, reverse_complement(query), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=13) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
target            6 ACGCTCGAGCAGCTACG 23
                  0 |||.||||||.|||||| 17
query            22 ACGATCGAGCNGCTACG  5
""",
        )
        self.assertEqual(alignment.shape, (2, 17))
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (17 aligned letters; 15 identities; 2 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 17:
        identities = 15,
        mismatches = 2.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 17)
        self.assertEqual(counts.identities, 15)
        self.assertEqual(counts.mismatches, 2)
        counts = alignment.counts("N")
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (17 aligned letters; 15 identities; 1 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 17:
        identities = 15,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 17)
        self.assertEqual(counts.identities, 15)
        self.assertEqual(counts.mismatches, 1)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 13.0; substitution score = 13.0; gap score = 0.0; 17 aligned letters; 15 identities; 1 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 13.0:
        substitution_score = 13.0,
        gap_score = 0.0.
    aligned = 17:
        identities = 15,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 17)
        self.assertEqual(counts.identities, 15)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 13.0)
        self.assertEqual(
            alignment.format("psl"),
            """\
15	1	0	1	0	0	0	0	-	query	22	5	22	target	23	6	23	1	17,	0,	6,
""",
        )
        self.assertEqual(
            alignment.format("bed"),
            """\
target	6	23	query	13	-	6	23	0	1	17,	0,
""",
        )
        self.assertEqual(
            alignment.format("sam"),
            """\
query	16	target	7	255	17M5S	*	0	0	ACGATCGAGCNGCTACGCCCNC	*	AS:i:13
""",
        )
        # use Seq objects for target and query
        alignments = aligner.align(Seq(target), Seq(query))
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=13) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
target            6 ACGCTCGAGCAGCTACG 23
                  0 |||.||||||.|||||| 17
query             0 ACGATCGAGCNGCTACG 17
""",
        )
        self.assertEqual(alignment.shape, (2, 17))
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (17 aligned letters; 15 identities; 2 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 17:
        identities = 15,
        mismatches = 2.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 17)
        self.assertEqual(counts.identities, 15)
        self.assertEqual(counts.mismatches, 2)
        counts = alignment.counts("N")
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (17 aligned letters; 15 identities; 1 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 17:
        identities = 15,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 17)
        self.assertEqual(counts.identities, 15)
        self.assertEqual(counts.mismatches, 1)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 13.0; substitution score = 13.0; gap score = 0.0; 17 aligned letters; 15 identities; 1 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 13.0:
        substitution_score = 13.0,
        gap_score = 0.0.
    aligned = 17:
        identities = 15,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 17)
        self.assertEqual(counts.identities, 15)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 13.0)
        self.assertEqual(
            alignment.format("psl"),
            """\
15	1	0	1	0	0	0	0	+	query	22	0	17	target	23	6	23	1	17,	0,	6,
""",
        )
        self.assertEqual(
            alignment.format("bed"),
            """\
target	6	23	query	13	+	6	23	0	1	17,	0,
""",
        )
        self.assertEqual(
            alignment.format("sam"),
            """\
query	0	target	7	255	17M5S	*	0	0	ACGATCGAGCNGCTACGCCCNC	*	AS:i:13
""",
        )
        alignments = aligner.align(
            Seq(target), Seq(query).reverse_complement(), strand="-"
        )
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=13) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
target            6 ACGCTCGAGCAGCTACG 23
                  0 |||.||||||.|||||| 17
query            22 ACGATCGAGCNGCTACG  5
""",
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (17 aligned letters; 15 identities; 2 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 17:
        identities = 15,
        mismatches = 2.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 17)
        self.assertEqual(counts.identities, 15)
        self.assertEqual(counts.mismatches, 2)
        counts = alignment.counts("N")
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (17 aligned letters; 15 identities; 1 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 17:
        identities = 15,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 17)
        self.assertEqual(counts.identities, 15)
        self.assertEqual(counts.mismatches, 1)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 13.0; substitution score = 13.0; gap score = 0.0; 17 aligned letters; 15 identities; 1 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 13.0:
        substitution_score = 13.0,
        gap_score = 0.0.
    aligned = 17:
        identities = 15,
        mismatches = 1.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 17)
        self.assertEqual(counts.identities, 15)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 13.0)
        self.assertEqual(alignment.shape, (2, 17))
        self.assertEqual(
            alignment.format("psl"),
            """\
15	1	0	1	0	0	0	0	-	query	22	5	22	target	23	6	23	1	17,	0,	6,
""",
        )
        self.assertEqual(
            alignment.format("bed"),
            """\
target	6	23	query	13	-	6	23	0	1	17,	0,
""",
        )
        self.assertEqual(
            alignment.format("sam"),
            """\
query	16	target	7	255	17M5S	*	0	0	ACGATCGAGCNGCTACGCCCNC	*	AS:i:13
""",
        )
        # global alignment
        aligner.mode = "global"
        aligner.end_gap_score = 0
        # use strings for target and query
        alignments = aligner.align(target, query)
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=13) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
target            0 TTTTTNACGCTCGAGCAGCTACG----- 23
                  0 ------|||.||||||.||||||----- 28
query             0 ------ACGATCGAGCNGCTACGCCCNC 22
""",
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (17 aligned letters; 15 identities; 2 mismatches; 11 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 17:
        identities = 15,
        mismatches = 2.
    gaps = 11:
        left_gaps = 6:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 6:
                open_left_deletions = 1,
                extend_left_deletions = 5;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 5:
            right_insertions = 5:
                open_right_insertions = 1,
                extend_right_insertions = 4;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 17)
        self.assertEqual(counts.identities, 15)
        self.assertEqual(counts.mismatches, 2)
        counts = alignment.counts("N")
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (17 aligned letters; 15 identities; 1 mismatches; 11 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 17:
        identities = 15,
        mismatches = 1.
    gaps = 11:
        left_gaps = 6:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 6:
                open_left_deletions = 1,
                extend_left_deletions = 5;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 5:
            right_insertions = 5:
                open_right_insertions = 1,
                extend_right_insertions = 4;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 17)
        self.assertEqual(counts.identities, 15)
        self.assertEqual(counts.mismatches, 1)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 13.0; substitution score = 13.0; gap score = 0.0; 17 aligned letters; 15 identities; 1 mismatches; 11 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 13.0:
        substitution_score = 13.0,
        gap_score = 0.0.
    aligned = 17:
        identities = 15,
        mismatches = 1.
    gaps = 11:
        left_gaps = 6:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 6:
                open_left_deletions = 1,
                extend_left_deletions = 5;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 5:
            right_insertions = 5:
                open_right_insertions = 1,
                extend_right_insertions = 4;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 17)
        self.assertEqual(counts.identities, 15)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 13.0)
        self.assertEqual(alignment.shape, (2, 28))
        self.assertEqual(
            alignment.format("psl"),
            """\
15	1	0	1	0	0	0	0	+	query	22	0	17	target	23	6	23	1	17,	0,	6,
""",
        )
        self.assertEqual(
            alignment.format("bed"),
            """\
target	6	23	query	13	+	6	23	0	1	17,	0,
""",
        )
        self.assertEqual(
            alignment.format("sam"),
            """\
query	0	target	1	255	6D17M5I	*	0	0	ACGATCGAGCNGCTACGCCCNC	*	AS:i:13
""",
        )
        alignments = aligner.align(target, reverse_complement(query), strand="-")
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=13) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
target            0 TTTTTNACGCTCGAGCAGCTACG----- 23
                  0 ------|||.||||||.||||||----- 28
query            22 ------ACGATCGAGCNGCTACGCCCNC  0
""",
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (17 aligned letters; 15 identities; 2 mismatches; 11 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 17:
        identities = 15,
        mismatches = 2.
    gaps = 11:
        left_gaps = 6:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 6:
                open_left_deletions = 1,
                extend_left_deletions = 5;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 5:
            right_insertions = 5:
                open_right_insertions = 1,
                extend_right_insertions = 4;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 17)
        self.assertEqual(counts.identities, 15)
        self.assertEqual(counts.mismatches, 2)
        counts = alignment.counts("N")
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (17 aligned letters; 15 identities; 1 mismatches; 11 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 17:
        identities = 15,
        mismatches = 1.
    gaps = 11:
        left_gaps = 6:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 6:
                open_left_deletions = 1,
                extend_left_deletions = 5;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 5:
            right_insertions = 5:
                open_right_insertions = 1,
                extend_right_insertions = 4;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 17)
        self.assertEqual(counts.identities, 15)
        self.assertEqual(counts.mismatches, 1)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 13.0; substitution score = 13.0; gap score = 0.0; 17 aligned letters; 15 identities; 1 mismatches; 11 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 13.0:
        substitution_score = 13.0,
        gap_score = 0.0.
    aligned = 17:
        identities = 15,
        mismatches = 1.
    gaps = 11:
        left_gaps = 6:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 6:
                open_left_deletions = 1,
                extend_left_deletions = 5;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 5:
            right_insertions = 5:
                open_right_insertions = 1,
                extend_right_insertions = 4;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 17)
        self.assertEqual(counts.identities, 15)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 13.0)
        self.assertEqual(alignment.shape, (2, 28))
        self.assertEqual(
            alignment.format("psl"),
            """\
15	1	0	1	0	0	0	0	-	query	22	5	22	target	23	6	23	1	17,	0,	6,
""",
        )
        self.assertEqual(
            alignment.format("bed"),
            """\
target	6	23	query	13	-	6	23	0	1	17,	0,
""",
        )
        self.assertEqual(
            alignment.format("sam"),
            """\
query	16	target	1	255	6D17M5I	*	0	0	ACGATCGAGCNGCTACGCCCNC	*	AS:i:13
""",
        )
        # use Seq objects for target and query
        alignments = aligner.align(Seq(target), Seq(query))
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=13) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
target            0 TTTTTNACGCTCGAGCAGCTACG----- 23
                  0 ------|||.||||||.||||||----- 28
query             0 ------ACGATCGAGCNGCTACGCCCNC 22
""",
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (17 aligned letters; 15 identities; 2 mismatches; 11 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 17:
        identities = 15,
        mismatches = 2.
    gaps = 11:
        left_gaps = 6:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 6:
                open_left_deletions = 1,
                extend_left_deletions = 5;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 5:
            right_insertions = 5:
                open_right_insertions = 1,
                extend_right_insertions = 4;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 17)
        self.assertEqual(counts.identities, 15)
        self.assertEqual(counts.mismatches, 2)
        counts = alignment.counts("N")
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (17 aligned letters; 15 identities; 1 mismatches; 11 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 17:
        identities = 15,
        mismatches = 1.
    gaps = 11:
        left_gaps = 6:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 6:
                open_left_deletions = 1,
                extend_left_deletions = 5;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 5:
            right_insertions = 5:
                open_right_insertions = 1,
                extend_right_insertions = 4;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 17)
        self.assertEqual(counts.identities, 15)
        self.assertEqual(counts.mismatches, 1)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 13.0; substitution score = 13.0; gap score = 0.0; 17 aligned letters; 15 identities; 1 mismatches; 11 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 13.0:
        substitution_score = 13.0,
        gap_score = 0.0.
    aligned = 17:
        identities = 15,
        mismatches = 1.
    gaps = 11:
        left_gaps = 6:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 6:
                open_left_deletions = 1,
                extend_left_deletions = 5;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 5:
            right_insertions = 5:
                open_right_insertions = 1,
                extend_right_insertions = 4;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 17)
        self.assertEqual(counts.identities, 15)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 13.0)
        self.assertEqual(alignment.shape, (2, 28))
        self.assertEqual(
            alignment.format("psl"),
            """\
15	1	0	1	0	0	0	0	+	query	22	0	17	target	23	6	23	1	17,	0,	6,
""",
        )
        self.assertEqual(
            alignment.format("bed"),
            """\
target	6	23	query	13	+	6	23	0	1	17,	0,
""",
        )
        self.assertEqual(
            alignment.format("sam"),
            """\
query	0	target	1	255	6D17M5I	*	0	0	ACGATCGAGCNGCTACGCCCNC	*	AS:i:13
""",
        )
        alignments = aligner.align(
            Seq(target), Seq(query).reverse_complement(), strand="-"
        )
        self.assertEqual(
            repr(alignments),
            f"""\
<PairwiseAlignments object (1 alignment; score=13) at {hex(id(alignments))}>""",
        )
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
target            0 TTTTTNACGCTCGAGCAGCTACG----- 23
                  0 ------|||.||||||.||||||----- 28
query            22 ------ACGATCGAGCNGCTACGCCCNC  0
""",
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (17 aligned letters; 15 identities; 2 mismatches; 11 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 17:
        identities = 15,
        mismatches = 2.
    gaps = 11:
        left_gaps = 6:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 6:
                open_left_deletions = 1,
                extend_left_deletions = 5;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 5:
            right_insertions = 5:
                open_right_insertions = 1,
                extend_right_insertions = 4;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 17)
        self.assertEqual(counts.identities, 15)
        self.assertEqual(counts.mismatches, 2)
        counts = alignment.counts("N")
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (17 aligned letters; 15 identities; 1 mismatches; 11 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 17:
        identities = 15,
        mismatches = 1.
    gaps = 11:
        left_gaps = 6:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 6:
                open_left_deletions = 1,
                extend_left_deletions = 5;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 5:
            right_insertions = 5:
                open_right_insertions = 1,
                extend_right_insertions = 4;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 17)
        self.assertEqual(counts.identities, 15)
        self.assertEqual(counts.mismatches, 1)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 13.0; substitution score = 13.0; gap score = 0.0; 17 aligned letters; 15 identities; 1 mismatches; 11 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 13.0:
        substitution_score = 13.0,
        gap_score = 0.0.
    aligned = 17:
        identities = 15,
        mismatches = 1.
    gaps = 11:
        left_gaps = 6:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 6:
                open_left_deletions = 1,
                extend_left_deletions = 5;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 5:
            right_insertions = 5:
                open_right_insertions = 1,
                extend_right_insertions = 4;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.aligned, 17)
        self.assertEqual(counts.identities, 15)
        self.assertEqual(counts.mismatches, 1)
        self.assertAlmostEqual(counts.score, 13.0)
        self.assertEqual(alignment.shape, (2, 28))
        self.assertEqual(
            alignment.format("psl"),
            """\
15	1	0	1	0	0	0	0	-	query	22	5	22	target	23	6	23	1	17,	0,	6,
""",
        )
        self.assertEqual(
            alignment.format("bed"),
            """\
target	6	23	query	13	-	6	23	0	1	17,	0,
""",
        )
        self.assertEqual(
            alignment.format("sam"),
            """\
query	16	target	1	255	6D17M5I	*	0	0	ACGATCGAGCNGCTACGCCCNC	*	AS:i:13
""",
        )


class TestAlgorithmRestrictions(unittest.TestCase):
    def test_fogsaa_restrictions(self):
        aligner = Align.PairwiseAligner(mode="fogsaa")
        aligner.match_score = -1
        with self.assertWarns(BiopythonWarning):
            aligner.score("AAAAAAAAAAAA", "AAAAATAAAAAA")
        aligner.mismatch_score = 1
        with self.assertWarns(BiopythonWarning):
            aligner.score("AAAAAAAAAAAA", "AAAAATAAAAAA")
        aligner.gap_score = 1
        with self.assertWarns(BiopythonWarning):
            aligner.score("AAAAAAAAAAAA", "AAAAATAAAAAA")


class TestCounts(unittest.TestCase):

    def check_counts(self, counts):
        self.assertEqual(counts.left_insertions, 2)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 1)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 7)
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 4)
        self.assertEqual(counts.mismatches, 1)

    def check_counts_and_score(self, counts):
        self.check_counts(counts)
        self.assertAlmostEqual(counts.score, -7.0)

    def test_string_bytes(self):
        aligner = Align.PairwiseAligner()
        aligner.mismatch_score = -1
        aligner.gap_score = -1
        alignments = aligner.align("TTACGTCCCCCCC", "ACTTTGT")
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
target            0 --TTACGTCCCCCCC 13
                  0 --||.-||------- 15
query             0 ACTTT-GT-------  7
""",
        )
        self.assertAlmostEqual(alignment.score, -7.0)
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (5 aligned letters; 4 identities; 1 mismatches; 10 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 5:
        identities = 4,
        mismatches = 1.
    gaps = 10:
        left_gaps = 2:
            left_insertions = 2:
                open_left_insertions = 1,
                extend_left_insertions = 1;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 7:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 7:
                open_right_deletions = 1,
                extend_right_deletions = 6.
""",
        )
        self.check_counts(counts)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = -7.0; substitution score = 3.0; gap score = -10.0; 5 aligned letters; 4 identities; 1 mismatches; 10 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = -7.0:
        substitution_score = 3.0,
        gap_score = -10.0.
    aligned = 5:
        identities = 4,
        mismatches = 1.
    gaps = 10:
        left_gaps = 2:
            left_insertions = 2:
                open_left_insertions = 1,
                extend_left_insertions = 1;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 7:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 7:
                open_right_deletions = 1,
                extend_right_deletions = 6.
""",
        )
        self.check_counts_and_score(counts)
        alignments = aligner.align("TTACGTCCCCCCC", b"ACTTTGT")
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
-- -- T  T  A  C G  T  C C C C C C C
-- -- .  .  .  - .  .  - - - - - - -
65 67 84 84 84 - 71 84 - - - - - - -
""",
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (5 aligned letters; 4 identities; 1 mismatches; 10 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 5:
        identities = 4,
        mismatches = 1.
    gaps = 10:
        left_gaps = 2:
            left_insertions = 2:
                open_left_insertions = 1,
                extend_left_insertions = 1;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 7:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 7:
                open_right_deletions = 1,
                extend_right_deletions = 6.
""",
        )
        self.check_counts(counts)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = -7.0; substitution score = 3.0; gap score = -10.0; 5 aligned letters; 4 identities; 1 mismatches; 10 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = -7.0:
        substitution_score = 3.0,
        gap_score = -10.0.
    aligned = 5:
        identities = 4,
        mismatches = 1.
    gaps = 10:
        left_gaps = 2:
            left_insertions = 2:
                open_left_insertions = 1,
                extend_left_insertions = 1;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 7:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 7:
                open_right_deletions = 1,
                extend_right_deletions = 6.
""",
        )
        self.check_counts_and_score(counts)
        alignments = aligner.align(b"TTACGTCCCCCCC", "ACTTTGT")
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
- - 84 84 65 67 71 84 67 67 67 67 67 67 67
- - .  .  .  -- .  .  -- -- -- -- -- -- --
A C T  T  T  -- G  T  -- -- -- -- -- -- --
""",
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (5 aligned letters; 4 identities; 1 mismatches; 10 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 5:
        identities = 4,
        mismatches = 1.
    gaps = 10:
        left_gaps = 2:
            left_insertions = 2:
                open_left_insertions = 1,
                extend_left_insertions = 1;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 7:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 7:
                open_right_deletions = 1,
                extend_right_deletions = 6.
""",
        )
        self.check_counts(counts)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = -7.0; substitution score = 3.0; gap score = -10.0; 5 aligned letters; 4 identities; 1 mismatches; 10 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = -7.0:
        substitution_score = 3.0,
        gap_score = -10.0.
    aligned = 5:
        identities = 4,
        mismatches = 1.
    gaps = 10:
        left_gaps = 2:
            left_insertions = 2:
                open_left_insertions = 1,
                extend_left_insertions = 1;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 7:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 7:
                open_right_deletions = 1,
                extend_right_deletions = 6.
""",
        )
        self.check_counts_and_score(counts)
        alignments = aligner.align(b"TTACGTCCCCCCC", b"ACTTTGT")
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
-- -- 84 84 65 67 71 84 67 67 67 67 67 67 67
-- -- || || .. -- || || -- -- -- -- -- -- --
65 67 84 84 84 -- 71 84 -- -- -- -- -- -- --
""",
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (5 aligned letters; 4 identities; 1 mismatches; 10 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 5:
        identities = 4,
        mismatches = 1.
    gaps = 10:
        left_gaps = 2:
            left_insertions = 2:
                open_left_insertions = 1,
                extend_left_insertions = 1;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 7:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 7:
                open_right_deletions = 1,
                extend_right_deletions = 6.
""",
        )
        self.check_counts(counts)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = -7.0; substitution score = 3.0; gap score = -10.0; 5 aligned letters; 4 identities; 1 mismatches; 10 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = -7.0:
        substitution_score = 3.0,
        gap_score = -10.0.
    aligned = 5:
        identities = 4,
        mismatches = 1.
    gaps = 10:
        left_gaps = 2:
            left_insertions = 2:
                open_left_insertions = 1,
                extend_left_insertions = 1;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 7:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 7:
                open_right_deletions = 1,
                extend_right_deletions = 6.
""",
        )
        self.check_counts_and_score(counts)
        alignments = aligner.align(Seq("TTACGTCCCCCCC"), Seq("ACTTTGT"))
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
target            0 --TTACGTCCCCCCC 13
                  0 --||.-||------- 15
query             0 ACTTT-GT-------  7
""",
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (5 aligned letters; 4 identities; 1 mismatches; 10 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 5:
        identities = 4,
        mismatches = 1.
    gaps = 10:
        left_gaps = 2:
            left_insertions = 2:
                open_left_insertions = 1,
                extend_left_insertions = 1;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 7:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 7:
                open_right_deletions = 1,
                extend_right_deletions = 6.
""",
        )
        self.check_counts(counts)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = -7.0; substitution score = 3.0; gap score = -10.0; 5 aligned letters; 4 identities; 1 mismatches; 10 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = -7.0:
        substitution_score = 3.0,
        gap_score = -10.0.
    aligned = 5:
        identities = 4,
        mismatches = 1.
    gaps = 10:
        left_gaps = 2:
            left_insertions = 2:
                open_left_insertions = 1,
                extend_left_insertions = 1;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 7:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 7:
                open_right_deletions = 1,
                extend_right_deletions = 6.
""",
        )
        self.check_counts_and_score(counts)
        alignments = aligner.align(Seq("TTACGTCCCCCCC"), "ACTTTGT")
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
target            0 --TTACGTCCCCCCC 13
                  0 --||.-||------- 15
query             0 ACTTT-GT-------  7
""",
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (5 aligned letters; 4 identities; 1 mismatches; 10 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 5:
        identities = 4,
        mismatches = 1.
    gaps = 10:
        left_gaps = 2:
            left_insertions = 2:
                open_left_insertions = 1,
                extend_left_insertions = 1;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 7:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 7:
                open_right_deletions = 1,
                extend_right_deletions = 6.
""",
        )
        self.check_counts(counts)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = -7.0; substitution score = 3.0; gap score = -10.0; 5 aligned letters; 4 identities; 1 mismatches; 10 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = -7.0:
        substitution_score = 3.0,
        gap_score = -10.0.
    aligned = 5:
        identities = 4,
        mismatches = 1.
    gaps = 10:
        left_gaps = 2:
            left_insertions = 2:
                open_left_insertions = 1,
                extend_left_insertions = 1;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 7:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 7:
                open_right_deletions = 1,
                extend_right_deletions = 6.
""",
        )
        self.check_counts_and_score(counts)
        alignments = aligner.align("TTACGTCCCCCCC", Seq("ACTTTGT"))
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
target            0 --TTACGTCCCCCCC 13
                  0 --||.-||------- 15
query             0 ACTTT-GT-------  7
""",
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (5 aligned letters; 4 identities; 1 mismatches; 10 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 5:
        identities = 4,
        mismatches = 1.
    gaps = 10:
        left_gaps = 2:
            left_insertions = 2:
                open_left_insertions = 1,
                extend_left_insertions = 1;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 7:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 7:
                open_right_deletions = 1,
                extend_right_deletions = 6.
""",
        )
        self.check_counts(counts)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = -7.0; substitution score = 3.0; gap score = -10.0; 5 aligned letters; 4 identities; 1 mismatches; 10 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = -7.0:
        substitution_score = 3.0,
        gap_score = -10.0.
    aligned = 5:
        identities = 4,
        mismatches = 1.
    gaps = 10:
        left_gaps = 2:
            left_insertions = 2:
                open_left_insertions = 1,
                extend_left_insertions = 1;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 7:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 7:
                open_right_deletions = 1,
                extend_right_deletions = 6.
""",
        )
        self.check_counts_and_score(counts)
        alignments = aligner.align(Seq("TTACGTCCCCCCC"), b"ACTTTGT")
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
-- -- T  T  A  C G  T  C C C C C C C
-- -- .  .  .  - .  .  - - - - - - -
65 67 84 84 84 - 71 84 - - - - - - -
""",
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (5 aligned letters; 4 identities; 1 mismatches; 10 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 5:
        identities = 4,
        mismatches = 1.
    gaps = 10:
        left_gaps = 2:
            left_insertions = 2:
                open_left_insertions = 1,
                extend_left_insertions = 1;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 7:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 7:
                open_right_deletions = 1,
                extend_right_deletions = 6.
""",
        )
        self.check_counts(counts)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = -7.0; substitution score = 3.0; gap score = -10.0; 5 aligned letters; 4 identities; 1 mismatches; 10 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = -7.0:
        substitution_score = 3.0,
        gap_score = -10.0.
    aligned = 5:
        identities = 4,
        mismatches = 1.
    gaps = 10:
        left_gaps = 2:
            left_insertions = 2:
                open_left_insertions = 1,
                extend_left_insertions = 1;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 7:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 7:
                open_right_deletions = 1,
                extend_right_deletions = 6.
""",
        )
        self.check_counts_and_score(counts)
        alignments = aligner.align(b"TTACGTCCCCCCC", Seq("ACTTTGT"))
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
- - 84 84 65 67 71 84 67 67 67 67 67 67 67
- - .  .  .  -- .  .  -- -- -- -- -- -- --
A C T  T  T  -- G  T  -- -- -- -- -- -- --
""",
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (5 aligned letters; 4 identities; 1 mismatches; 10 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 5:
        identities = 4,
        mismatches = 1.
    gaps = 10:
        left_gaps = 2:
            left_insertions = 2:
                open_left_insertions = 1,
                extend_left_insertions = 1;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 7:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 7:
                open_right_deletions = 1,
                extend_right_deletions = 6.
""",
        )
        self.check_counts(counts)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = -7.0; substitution score = 3.0; gap score = -10.0; 5 aligned letters; 4 identities; 1 mismatches; 10 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = -7.0:
        substitution_score = 3.0,
        gap_score = -10.0.
    aligned = 5:
        identities = 4,
        mismatches = 1.
    gaps = 10:
        left_gaps = 2:
            left_insertions = 2:
                open_left_insertions = 1,
                extend_left_insertions = 1;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 7:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 7:
                open_right_deletions = 1,
                extend_right_deletions = 6.
""",
        )
        self.check_counts_and_score(counts)

    def check_incomplete_nucleotide_sequence(self, counts):
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 4)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.aligned, 13)
        self.assertEqual(counts.identities, 12)
        self.assertEqual(counts.mismatches, 1)

    def check_incomplete_nucleotide_sequence_switched(self, counts):
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.internal_insertions, 4)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.aligned, 13)
        self.assertEqual(counts.identities, 12)
        self.assertEqual(counts.mismatches, 1)

    def test_incomplete_nucleotide_sequence(self):
        aligner = Align.PairwiseAligner()
        aligner_blastn = Align.PairwiseAligner("blastn")
        seqA = Seq({10: "TTACGT", 20: "CCCCCCC"}, length=50)
        seqB = Seq("TTACGTCCCTCCC")
        coordinates = np.array([[10, 16, 20, 27], [0, 6, 6, 13]])
        alignment = Align.Alignment([seqA, seqB], coordinates)
        self.assertEqual(
            str(alignment),
            """\
target           10 TTACGT????CCCCCCC 27
                  0 ||||||----|||.||| 17
query             0 TTACGT----CCCTCCC 13
""",
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 12.0; substitution score = 12.0; gap score = 0.0; 13 aligned letters; 12 identities; 1 mismatches; 4 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 12.0:
        substitution_score = 12.0,
        gap_score = 0.0.
    aligned = 13:
        identities = 12,
        mismatches = 1.
    gaps = 4:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 4:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 4:
                open_internal_deletions = 1,
                extend_internal_deletions = 3;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.check_incomplete_nucleotide_sequence(counts)
        self.assertIsNone(counts.positives)
        self.assertAlmostEqual(counts.score, 12.0)
        counts = alignment.counts(aligner_blastn)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 8.0; substitution score = 21.0; gap score = -13.0; 13 aligned letters; 12 identities; 1 mismatches; 12 positives; 4 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 8.0:
        substitution_score = 21.0,
        gap_score = -13.0.
    aligned = 13:
        identities = 12,
        positives = 12,
        mismatches = 1.
    gaps = 4:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 4:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 4:
                open_internal_deletions = 1,
                extend_internal_deletions = 3;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.check_incomplete_nucleotide_sequence(counts)
        self.assertEqual(counts.positives, 12)
        self.assertAlmostEqual(counts.score, 8.0)
        alignment = Align.Alignment([seqB, seqA], coordinates[::-1])
        self.assertEqual(
            str(alignment),
            """\
target            0 TTACGT----CCCTCCC 13
                  0 ||||||----|||.||| 17
query            10 TTACGT????CCCCCCC 27
""",
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 12.0; substitution score = 12.0; gap score = 0.0; 13 aligned letters; 12 identities; 1 mismatches; 4 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 12.0:
        substitution_score = 12.0,
        gap_score = 0.0.
    aligned = 13:
        identities = 12,
        mismatches = 1.
    gaps = 4:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 4:
            internal_insertions = 4:
                open_internal_insertions = 1,
                extend_internal_insertions = 3;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.check_incomplete_nucleotide_sequence_switched(counts)
        self.assertIsNone(counts.positives)
        counts = alignment.counts(aligner_blastn)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 8.0; substitution score = 21.0; gap score = -13.0; 13 aligned letters; 12 identities; 1 mismatches; 12 positives; 4 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 8.0:
        substitution_score = 21.0,
        gap_score = -13.0.
    aligned = 13:
        identities = 12,
        positives = 12,
        mismatches = 1.
    gaps = 4:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 4:
            internal_insertions = 4:
                open_internal_insertions = 1,
                extend_internal_insertions = 3;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.check_incomplete_nucleotide_sequence_switched(counts)
        self.assertEqual(counts.positives, 12)
        seqA = Seq({10: "TTACGT", 20: "CCCTCCC"}, length=50)
        seqB = Seq({100: "TTACGTCCCCCCC"}, length=200)
        coordinates = np.array([[10, 16, 20, 27], [100, 106, 106, 113]])
        alignment = Align.Alignment([seqA, seqB], coordinates)
        self.assertEqual(
            str(alignment),
            """\
target           10 TTACGT????CCCTCCC  27
                  0 ||||||----|||.|||  17
query           100 TTACGT----CCCCCCC 113
""",
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 12.0; substitution score = 12.0; gap score = 0.0; 13 aligned letters; 12 identities; 1 mismatches; 4 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 12.0:
        substitution_score = 12.0,
        gap_score = 0.0.
    aligned = 13:
        identities = 12,
        mismatches = 1.
    gaps = 4:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 4:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 4:
                open_internal_deletions = 1,
                extend_internal_deletions = 3;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.check_incomplete_nucleotide_sequence(counts)
        self.assertIsNone(counts.positives)
        self.assertAlmostEqual(counts.score, 12.0)
        counts = alignment.counts(aligner_blastn)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 8.0; substitution score = 21.0; gap score = -13.0; 13 aligned letters; 12 identities; 1 mismatches; 12 positives; 4 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 8.0:
        substitution_score = 21.0,
        gap_score = -13.0.
    aligned = 13:
        identities = 12,
        positives = 12,
        mismatches = 1.
    gaps = 4:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 4:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 4:
                open_internal_deletions = 1,
                extend_internal_deletions = 3;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.check_incomplete_nucleotide_sequence(counts)
        self.assertEqual(counts.positives, 12)
        self.assertAlmostEqual(counts.score, 8.0)
        seqA_rc = seqA.reverse_complement()
        seqB_rc = seqB.reverse_complement()
        coordinates_rc = np.array(coordinates)
        coordinates_rc[0, :] = len(seqA) - coordinates[0, :]
        alignment = Align.Alignment([seqA_rc, seqB], coordinates_rc)
        self.assertEqual(
            str(alignment),
            """\
target           40 TTACGT????CCCTCCC  23
                  0 ||||||----|||.|||  17
query           100 TTACGT----CCCCCCC 113
""",
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 12.0; substitution score = 12.0; gap score = 0.0; 13 aligned letters; 12 identities; 1 mismatches; 4 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 12.0:
        substitution_score = 12.0,
        gap_score = 0.0.
    aligned = 13:
        identities = 12,
        mismatches = 1.
    gaps = 4:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 4:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 4:
                open_internal_deletions = 1,
                extend_internal_deletions = 3;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.check_incomplete_nucleotide_sequence(counts)
        self.assertIsNone(counts.positives)
        self.assertAlmostEqual(counts.score, 12.0)
        counts = alignment.counts(aligner_blastn)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 8.0; substitution score = 21.0; gap score = -13.0; 13 aligned letters; 12 identities; 1 mismatches; 12 positives; 4 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 8.0:
        substitution_score = 21.0,
        gap_score = -13.0.
    aligned = 13:
        identities = 12,
        positives = 12,
        mismatches = 1.
    gaps = 4:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 4:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 4:
                open_internal_deletions = 1,
                extend_internal_deletions = 3;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.check_incomplete_nucleotide_sequence(counts)
        self.assertEqual(counts.positives, 12)
        self.assertAlmostEqual(counts.score, 8.0)
        coordinates_rc = np.array(coordinates)
        coordinates_rc[1, :] = len(seqB) - coordinates[1, :]
        alignment = Align.Alignment([seqA, seqB_rc], coordinates_rc)
        self.assertEqual(
            str(alignment),
            """\
target           10 TTACGT????CCCTCCC 27
                  0 ||||||----|||.||| 17
query           100 TTACGT----CCCCCCC 87
""",
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 12.0; substitution score = 12.0; gap score = 0.0; 13 aligned letters; 12 identities; 1 mismatches; 4 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 12.0:
        substitution_score = 12.0,
        gap_score = 0.0.
    aligned = 13:
        identities = 12,
        mismatches = 1.
    gaps = 4:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 4:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 4:
                open_internal_deletions = 1,
                extend_internal_deletions = 3;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.check_incomplete_nucleotide_sequence(counts)
        self.assertIsNone(counts.positives)
        self.assertAlmostEqual(counts.score, 12.0)
        counts = alignment.counts(aligner_blastn)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 8.0; substitution score = 21.0; gap score = -13.0; 13 aligned letters; 12 identities; 1 mismatches; 12 positives; 4 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 8.0:
        substitution_score = 21.0,
        gap_score = -13.0.
    aligned = 13:
        identities = 12,
        positives = 12,
        mismatches = 1.
    gaps = 4:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 4:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 4:
                open_internal_deletions = 1,
                extend_internal_deletions = 3;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.check_incomplete_nucleotide_sequence(counts)
        self.assertEqual(counts.positives, 12)
        self.assertAlmostEqual(counts.score, 8.0)
        coordinates_rc = np.array(coordinates)
        coordinates_rc[0, :] = len(seqA) - coordinates[0, :]
        coordinates_rc[1, :] = len(seqB) - coordinates[1, :]
        alignment = Align.Alignment([seqA_rc, seqB_rc], coordinates_rc)
        self.assertEqual(
            str(alignment),
            """\
target           40 TTACGT????CCCTCCC 23
                  0 ||||||----|||.||| 17
query           100 TTACGT----CCCCCCC 87
""",
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 12.0; substitution score = 12.0; gap score = 0.0; 13 aligned letters; 12 identities; 1 mismatches; 4 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 12.0:
        substitution_score = 12.0,
        gap_score = 0.0.
    aligned = 13:
        identities = 12,
        mismatches = 1.
    gaps = 4:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 4:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 4:
                open_internal_deletions = 1,
                extend_internal_deletions = 3;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.check_incomplete_nucleotide_sequence(counts)
        self.assertIsNone(counts.positives)
        self.assertAlmostEqual(counts.score, 12.0)
        counts = alignment.counts(aligner_blastn)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 8.0; substitution score = 21.0; gap score = -13.0; 13 aligned letters; 12 identities; 1 mismatches; 12 positives; 4 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 8.0:
        substitution_score = 21.0,
        gap_score = -13.0.
    aligned = 13:
        identities = 12,
        positives = 12,
        mismatches = 1.
    gaps = 4:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 4:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 4:
                open_internal_deletions = 1,
                extend_internal_deletions = 3;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.check_incomplete_nucleotide_sequence(counts)
        self.assertEqual(counts.positives, 12)
        self.assertAlmostEqual(counts.score, 8.0)

    def check_blastp(self, counts):
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.internal_insertions, 2)
        self.assertEqual(counts.internal_deletions, 34)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.aligned, 444)
        self.assertEqual(counts.identities, 237)
        self.assertEqual(counts.mismatches, 207)
        self.assertEqual(counts.gaps, 36)

    def test_blastp(self):
        aligner = Align.PairwiseAligner()
        aligner_blastp = Align.PairwiseAligner("blastp")
        seqA = "MRPRPRSAPGKPRRRSRARLRSSRTPSGGASGGGGSSSSSSTATGGSGSSTGSPGGAASAPAPAPAGMYRSGERLLGSHALPAEQRDFLPLETTNNNNNHHQPGAWARRAGSSASSPPSASSSPHPSAAVPAADPADSASGSSNKRKRDNKASTYGLNYSLLQPSGGRAAGGGRADGGGVVYSGTPWKRRNYNQGVVGLHEEISDFYEYMSPRPEEEKMRMEVVNRIESVIKELWPSADVQIFGSFKTGLYLPTSDIDLVVFGKWENLPLWTLEEALRKHKVADEDSVKVLDKATVPIIKLTDSFTEVKVDISFNVQNGVRAADLIKDFTKKYPVLPYLVLVLKQFLLQRDLNEVFTGGIGSYSLFLMAVSFLQLHPREDACIPNTNYGVLLIEFFELYGRHFNYLKTGIRIKDGGSYVAKDEVQKNMLDGYRPSMLYIEDPLQPGNDVGRSSYGAMQVKQAFDYAYVVLSHAVSPIAKYYPNNETESILGRIIRVTDEVATYRDWISKQWGLKNRPEPSCNGNGVTLIVDTQQLDKCNNNLSEENEALGKCRSKTSESLSKHSSNSSSGPVSSSSATQSSSSDVDSDATPCKTPKQLLCRPSTGNRVGSQDVSLESSQAVGKMQSTQTTNTSNSTNKSQHGSARLFRSSSKGFQGTTQTSHGSLMTNKQHQGKSNNQYYHGKKRKHKRDAPLSDLCR"
        seqB = "MDPFVGWFQKEQEGPSLCTWLKIWETNAQEMGALLNQQLQQATTINGSTSSSSSSSNSGNNNNNNNNNIINNTITNTTNNTGNNSSAKPYLSRPYSSLNRVLNFRADSLEILQQQQQQQQLNGTTQRNSTNINTTSGGSTSSSADSTTNRDNNSPANSSSTNGPGAGTGTSTGAGGTGTNSPATTASSTAATTTGPATSMSDTSNNPPQSTTTPASRTNSIYYNPSRKKRPENKAGGAHYYMNNHMEMIAKYKGEPWRKPDYPYGEGVIGLHEEIEHFYQYVLPTPCEHAIRNEVVKRIEAVVHSIWPQAVVEIFGSFRTGLFLPTSDIDLVVLGLWEKLPLRTLEFELVSRGIAEACTVRVLDKASVPIIKLTDRETQVKVDISFNMQSGVQSAELIKKFKRDYPVLEKLVLVLKQFLLLRDLNEVFTGGISSYSLILMCISFLQMHPRGIYHDTANLGVLLLEFFELYGRRFNYMKIGISIKNGGRYMPKDELQRDMVDGHRPSLLCIEDPLTPGNDIGRSSYGVFQVQQAFKCAYRVLALAVSPLNLLGIDPRVNSILGRIIHITDDVIDYREWIRENFEHLVVVDRISPLPTAAPTAYATANGAPKYVIMPSGAVVQQLYHHPVQVPTAHGHSHAHSHSHGHAHPGAHLCQPYVTGTTVSAVTTTTTMAVVTVGVSAGGVQQQQQQQNATAHTHSQQQTQNQSQSRHRRGSTSSGDDSEDSKDGDVVETTSSAQEVVDIALSTPNGLANMSMPMPVHAVGMPASNSWSGNGNGNGNSSSSTGSSPEIAHIAAQEMDPELEDQQQQQQHQETSGGNGFIRPGDVGTGSNRGGGDGSGGRNYNQRNNHNSSGYYHQQYYVPPPMQQQLSKSNSSSNYHQQHHHSHSHGNHSHRQQHHHQQQHHHQQRPQHLRVGGGNRYQKSLGGSPIISAGNASNSSSNCSNSSSSSGSNNSRLPPLRGTLVNSSSAISIISISSESSIASSSSSSSRSGQDQQRDER"
        sA = {33: seqA[33:511]}
        sB = {136: seqB[136:582]}
        # fmt: off
        coordinates = np.array([[ 33,  44,  46,  80,  89,  97, 106, 153, 155,
                                 170, 181, 192, 192, 378, 379, 511],
                                [136, 147, 147, 181, 181, 189, 189, 236, 236,
                                 251, 251, 262, 264, 450, 450, 582]])
        # fmt: on
        sequences = [Seq(seqA), Seq(seqB)]
        alignment = Align.Alignment(sequences, coordinates)
        self.assertEqual(
            str(alignment),
            """\
target           33 GGSSSSSSTATGGSGSSTGSPGGAASAPAPAPAGMYRSGERLLGSHALPAEQRDFLPLET
                  0 |||.|||...|--......||....|...|........|....|...---------|..|
query           136 GGSTSSSADST--TNRDNNSPANSSSTNGPGAGTGTSTGAGGTGTNS---------PATT

target           93 TNNNNNHHQPGAWARRAGSSASSPPSASSSPHPSAAVPAADPADSASGSSNKRKRDNKAS
                 60 ....---------|......|.|....|..|..|...||..........|.|....|||.
query           185 ASST---------AATTTGPATSMSDTSNNPPQSTTTPASRTNSIYYNPSRKKRPENKAG

target          153 TYGLNYSLLQPSGGRAAGGGRADGGGVVYSGTPWKRRNY--NQGVVGLHEEISDFYEYMS
                120 --|..|.........|.-----------|.|.||....|--..||.||||||..||.|..
query           236 --GAHYYMNNHMEMIAK-----------YKGEPWRKPDYPYGEGVIGLHEEIEHFYQYVL

target          211 PRPEEEKMRMEVVNRIESVIKELWPSADVQIFGSFKTGLYLPTSDIDLVVFGKWENLPLW
                180 |.|.|...|.|||.|||.|....||.|.|.|||||.|||.||||||||||.|.||.|||.
query           283 PTPCEHAIRNEVVKRIEAVVHSIWPQAVVEIFGSFRTGLFLPTSDIDLVVLGLWEKLPLR

target          271 TLEEALRKHKVADEDSVKVLDKATVPIIKLTDSFTEVKVDISFNVQNGVRAADLIKDFTK
                240 |||..|.....|....|.|||||.||||||||..|.||||||||.|.||..|.|||.|..
query           343 TLEFELVSRGIAEACTVRVLDKASVPIIKLTDRETQVKVDISFNMQSGVQSAELIKKFKR

target          331 KYPVLPYLVLVLKQFLLQRDLNEVFTGGIGSYSLFLMAVSFLQLHPREDACIPNTNYGVL
                300 .||||..||||||||||.|||||||||||.||||.||..||||.|||-.......|.|||
query           403 DYPVLEKLVLVLKQFLLLRDLNEVFTGGISSYSLILMCISFLQMHPR-GIYHDTANLGVL

target          391 LIEFFELYGRHFNYLKTGIRIKDGGSYVAKDEVQKNMLDGYRPSMLYIEDPLQPGNDVGR
                360 |.||||||||.|||.|.||.||.||.|..|||.|..|.||.|||.|.|||||.||||.||
query           462 LLEFFELYGRRFNYMKIGISIKNGGRYMPKDELQRDMVDGHRPSLLCIEDPLTPGNDIGR

target          451 SSYGAMQVKQAFDYAYVVLSHAVSPIAKYYPNNETESILGRIIRVTDEVATYRDWISKQW
                420 ||||..||.|||..||.||..||||...........|||||||..||.|..||.||....
query           522 SSYGVFQVQQAFKCAYRVLALAVSPLNLLGIDPRVNSILGRIIHITDDVIDYREWIRENF

target          511 
                480 
query           582 
""",
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 237.0; substitution score = 237.0; gap score = 0.0; 444 aligned letters; 237 identities; 207 mismatches; 36 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 237.0:
        substitution_score = 237.0,
        gap_score = 0.0.
    aligned = 444:
        identities = 237,
        mismatches = 207.
    gaps = 36:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 36:
            internal_insertions = 2:
                open_internal_insertions = 1,
                extend_internal_insertions = 1;
            internal_deletions = 34:
                open_internal_deletions = 6,
                extend_internal_deletions = 28;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.check_blastp(counts)
        self.assertIsNone(counts.positives)
        self.assertAlmostEqual(counts.score, 237.0)
        counts = alignment.counts(aligner_blastp)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1071.0; substitution score = 1184.0; gap score = -113.0; 444 aligned letters; 237 identities; 207 mismatches; 306 positives; 36 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1071.0:
        substitution_score = 1184.0,
        gap_score = -113.0.
    aligned = 444:
        identities = 237,
        positives = 306,
        mismatches = 207.
    gaps = 36:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 36:
            internal_insertions = 2:
                open_internal_insertions = 1,
                extend_internal_insertions = 1;
            internal_deletions = 34:
                open_internal_deletions = 6,
                extend_internal_deletions = 28;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.check_blastp(counts)
        self.assertEqual(counts.positives, 306)
        self.assertAlmostEqual(counts.score, 1071.0)
        sequences = [Seq(seqA), Seq(seqB)]
        alignment = Align.Alignment(sequences, coordinates)
        counts = alignment.counts(aligner_blastp)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1071.0; substitution score = 1184.0; gap score = -113.0; 444 aligned letters; 237 identities; 207 mismatches; 306 positives; 36 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1071.0:
        substitution_score = 1184.0,
        gap_score = -113.0.
    aligned = 444:
        identities = 237,
        positives = 306,
        mismatches = 207.
    gaps = 36:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 36:
            internal_insertions = 2:
                open_internal_insertions = 1,
                extend_internal_insertions = 1;
            internal_deletions = 34:
                open_internal_deletions = 6,
                extend_internal_deletions = 28;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.check_blastp(counts)
        sequences = [Seq(seqA), Seq(sB, length=len(seqB))]
        alignment = Align.Alignment(sequences, coordinates)
        counts = alignment.counts(aligner_blastp)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1071.0; substitution score = 1184.0; gap score = -113.0; 444 aligned letters; 237 identities; 207 mismatches; 306 positives; 36 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1071.0:
        substitution_score = 1184.0,
        gap_score = -113.0.
    aligned = 444:
        identities = 237,
        positives = 306,
        mismatches = 207.
    gaps = 36:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 36:
            internal_insertions = 2:
                open_internal_insertions = 1,
                extend_internal_insertions = 1;
            internal_deletions = 34:
                open_internal_deletions = 6,
                extend_internal_deletions = 28;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.check_blastp(counts)
        sequences = [Seq(sA, length=len(seqA)), Seq(seqB)]
        alignment = Align.Alignment(sequences, coordinates)
        counts = alignment.counts(aligner_blastp)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1071.0; substitution score = 1184.0; gap score = -113.0; 444 aligned letters; 237 identities; 207 mismatches; 306 positives; 36 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1071.0:
        substitution_score = 1184.0,
        gap_score = -113.0.
    aligned = 444:
        identities = 237,
        positives = 306,
        mismatches = 207.
    gaps = 36:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 36:
            internal_insertions = 2:
                open_internal_insertions = 1,
                extend_internal_insertions = 1;
            internal_deletions = 34:
                open_internal_deletions = 6,
                extend_internal_deletions = 28;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.check_blastp(counts)
        sequences = [Seq(sA, length=len(seqA)), Seq(sB, length=len(seqB))]
        alignment = Align.Alignment(sequences, coordinates)
        counts = alignment.counts(aligner_blastp)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 1071.0; substitution score = 1184.0; gap score = -113.0; 444 aligned letters; 237 identities; 207 mismatches; 306 positives; 36 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 1071.0:
        substitution_score = 1184.0,
        gap_score = -113.0.
    aligned = 444:
        identities = 237,
        positives = 306,
        mismatches = 207.
    gaps = 36:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 36:
            internal_insertions = 2:
                open_internal_insertions = 1,
                extend_internal_insertions = 1;
            internal_deletions = 34:
                open_internal_deletions = 6,
                extend_internal_deletions = 28;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.check_blastp(counts)

    def test_string(self):
        aligner = Align.PairwiseAligner()
        aligner.gap_score = -2
        seqA = "あいうえおかきくけこ"
        seqB = "いうきあけ"
        alignments = aligner.align(seqA, seqB)
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertEqual(alignment[0], "あいうえおかきくけこ")
        self.assertEqual(alignment[1], "-いう---きあけ-")
        self.assertAlmostEqual(alignment.score, -6.0)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = -6.0; substitution score = 4.0; gap score = -10.0; 5 aligned letters; 4 identities; 1 mismatches; 5 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = -6.0:
        substitution_score = 4.0,
        gap_score = -10.0.
    aligned = 5:
        identities = 4,
        mismatches = 1.
    gaps = 5:
        left_gaps = 1:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 1:
                open_left_deletions = 1,
                extend_left_deletions = 0;
        internal_gaps = 3:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 3:
                open_internal_deletions = 1,
                extend_internal_deletions = 2;
        right_gaps = 1:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 1:
                open_right_deletions = 1,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 1)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 3)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 1)
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 4)
        self.assertEqual(counts.mismatches, 1)
        self.assertIsNone(counts.positives)
        self.assertAlmostEqual(counts.score, -6.0)

    def test_greek(self):
        aligner = Align.PairwiseAligner()
        seqA = "αβγγβα"
        seqB = "ABCBAαβγ"
        alignments = aligner.align(seqA, seqB)
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
-----αβγγβα
-----|||---
ABCBAαβγ---
""",
        )
        self.assertAlmostEqual(alignment.score, 3.0)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 3.0; substitution score = 3.0; gap score = 0.0; 3 aligned letters; 3 identities; 0 mismatches; 8 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 3.0:
        substitution_score = 3.0,
        gap_score = 0.0.
    aligned = 3:
        identities = 3,
        mismatches = 0.
    gaps = 8:
        left_gaps = 5:
            left_insertions = 5:
                open_left_insertions = 1,
                extend_left_insertions = 4;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 3:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 3:
                open_right_deletions = 1,
                extend_right_deletions = 2.
""",
        )
        self.assertEqual(counts.left_insertions, 5)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 3)
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.positives)
        self.assertAlmostEqual(counts.score, 3.0)
        alignments = aligner.align(seqA, seqB)
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
-----αβγγβα
-----|||---
ABCBAαβγ---
""",
        )
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 3.0; substitution score = 3.0; gap score = 0.0; 3 aligned letters; 3 identities; 0 mismatches; 8 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 3.0:
        substitution_score = 3.0,
        gap_score = 0.0.
    aligned = 3:
        identities = 3,
        mismatches = 0.
    gaps = 8:
        left_gaps = 5:
            left_insertions = 5:
                open_left_insertions = 1,
                extend_left_insertions = 4;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 3:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 3:
                open_right_deletions = 1,
                extend_right_deletions = 2.
""",
        )
        self.assertEqual(counts.left_insertions, 5)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 3)
        self.assertEqual(counts.aligned, 3)
        self.assertEqual(counts.identities, 3)
        self.assertEqual(counts.mismatches, 0)
        self.assertIsNone(counts.positives)
        self.assertAlmostEqual(counts.score, 3.0)
        substitution_matrix = Array("ABCαβγ", dims=2)
        alphabet = substitution_matrix.alphabet
        for c in alphabet:
            substitution_matrix[c, c] = 2.0
        for c1, c2 in zip(alphabet[:3], alphabet[3:]):
            substitution_matrix[c1, c2] = 1.0
            substitution_matrix[c2, c1] = 1.0
        substitution_matrix["B", "β"] = 5.0
        substitution_matrix["β", "B"] = 5.0
        substitution_matrix["C", "γ"] = -0.1
        substitution_matrix["γ", "C"] = -0.1
        self.assertEqual(
            str(substitution_matrix),
            """\
    A   B    C   α   β    γ
A 2.0 0.0  0.0 1.0 0.0  0.0
B 0.0 2.0  0.0 0.0 5.0  0.0
C 0.0 0.0  2.0 0.0 0.0 -0.1
α 1.0 0.0  0.0 2.0 0.0  0.0
β 0.0 5.0  0.0 0.0 2.0  0.0
γ 0.0 0.0 -0.1 0.0 0.0  2.0
""",
        )
        aligner.substitution_matrix = substitution_matrix
        aligner.gap_score = -1.0
        alignments = aligner.align(seqA, seqB)
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
αβγγβ-α--
...-.-|--
ABC-BAαβγ
""",
        )
        self.assertAlmostEqual(alignment.score, 8.9)
        counts = alignment.counts(aligner)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (score = 8.9; substitution score = 12.9; gap score = -4.0; 5 aligned letters; 1 identities; 4 mismatches; 4 positives; 4 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    score = 8.9:
        substitution_score = 12.9,
        gap_score = -4.0.
    aligned = 5:
        identities = 1,
        positives = 4,
        mismatches = 4.
    gaps = 4:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 2:
            right_insertions = 2:
                open_right_insertions = 1,
                extend_right_insertions = 1;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.internal_insertions, 1)
        self.assertEqual(counts.internal_deletions, 1)
        self.assertEqual(counts.right_insertions, 2)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.aligned, 5)
        self.assertEqual(counts.identities, 1)
        self.assertEqual(counts.mismatches, 4)
        self.assertEqual(counts.positives, 4)
        self.assertAlmostEqual(counts.score, 8.9)
        counts = alignment.counts(substitution_matrix)
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (substitution score = 12.9; 5 aligned letters; 1 identities; 4 mismatches; 4 positives; 4 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    substitution_score = 12.9,
    aligned = 5:
        identities = 1,
        positives = 4,
        mismatches = 4.
    gaps = 4:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 2:
            internal_insertions = 1:
                open_internal_insertions = 1,
                extend_internal_insertions = 0;
            internal_deletions = 1:
                open_internal_deletions = 1,
                extend_internal_deletions = 0;
        right_gaps = 2:
            right_insertions = 2:
                open_right_insertions = 1,
                extend_right_insertions = 1;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
