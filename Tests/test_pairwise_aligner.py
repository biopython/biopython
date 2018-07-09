#!/usr/bin/env python

# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


import unittest

from Bio import Align


class TestAlignerProperties(unittest.TestCase):

    def test_aligner_property_epsilon(self):
        aligner = Align.PairwiseAligner()
        self.assertAlmostEqual(aligner.epsilon, 1.e-6)
        aligner.epsilon = 1.e-4
        self.assertAlmostEqual(aligner.epsilon, 1.e-4)
        aligner.epsilon = 1.e-8
        self.assertAlmostEqual(aligner.epsilon, 1.e-8)
        with self.assertRaises(TypeError):
            aligner.epsilon = 'not a number'
        with self.assertRaises(TypeError):
            aligner.epsilon = None

    def test_aligner_property_mode(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        self.assertEqual(aligner.mode, 'global')
        aligner.mode = 'local'
        self.assertEqual(aligner.mode, 'local')
        with self.assertRaises(ValueError):
            aligner.mode = 'wrong'

    def test_aligner_property_match_mismatch(self):
        aligner = Align.PairwiseAligner()
        aligner.match = 3.0
        self.assertAlmostEqual(aligner.match, 3.0)
        aligner.mismatch = -2.0
        self.assertAlmostEqual(aligner.mismatch, -2.0)
        with self.assertRaises(ValueError):
            aligner.match = 'not a number'
        with self.assertRaises(ValueError):
            aligner.mismatch = 'not a number'
        self.assertEqual(str(aligner), """\
Pairwise sequence aligner with parameters
  match score: 3.000000
  mismatch score: -2.000000
  target open gap score: 0.000000
  target extend gap score: 0.000000
  target left open gap score: 0.000000
  target left extend gap score: 0.000000
  target right open gap score: 0.000000
  target right extend gap score: 0.000000
  query open gap score: 0.000000
  query extend gap score: 0.000000
  query left open gap score: 0.000000
  query left extend gap score: 0.000000
  query right open gap score: 0.000000
  query right extend gap score: 0.000000
  mode: global
""")

    def test_aligner_property_gapscores(self):
        aligner = Align.PairwiseAligner()
        open_score, extend_score = (-5, -1)
        aligner.target_gap_open = open_score
        aligner.target_gap_extend = extend_score
        self.assertAlmostEqual(aligner.target_gap_open, open_score)
        self.assertAlmostEqual(aligner.target_gap_extend, extend_score)
        open_score, extend_score = (-6, -7)
        aligner.query_open_gap_score = open_score
        aligner.query_extend_gap_score = extend_score
        self.assertAlmostEqual(aligner.query_open_gap_score, open_score)
        self.assertAlmostEqual(aligner.query_extend_gap_score, extend_score)
        open_score, extend_score = (-3, -9)
        aligner.target_end_open = open_score
        aligner.target_end_extend = extend_score
        self.assertAlmostEqual(aligner.target_end_open, open_score)
        self.assertAlmostEqual(aligner.target_end_extend, extend_score)
        open_score, extend_score = (-1, -2)
        aligner.query_end_open_gap_score = open_score
        aligner.query_end_extend_gap_score = extend_score
        self.assertEqual(str(aligner), """\
Pairwise sequence aligner with parameters
  match score: 1.000000
  mismatch score: 0.000000
  target open gap score: 0.000000
  target extend gap score: 0.000000
  target left open gap score: 0.000000
  target left extend gap score: 0.000000
  target right open gap score: 0.000000
  target right extend gap score: 0.000000
  query open gap score: -6.000000
  query extend gap score: -7.000000
  query left open gap score: -1.000000
  query left extend gap score: -2.000000
  query right open gap score: -1.000000
  query right extend gap score: -2.000000
  mode: global
""")
        self.assertAlmostEqual(aligner.query_end_open_gap_score, open_score)
        self.assertAlmostEqual(aligner.query_end_extend_gap_score, extend_score)
        score = -3
        aligner.target_gap_score = score
        self.assertAlmostEqual(aligner.target_gap_score, score)
        self.assertAlmostEqual(aligner.target_open_gap_score, score)
        self.assertAlmostEqual(aligner.target_extend_gap_score, score)
        score = -2
        aligner.query_gap_score = score
        self.assertAlmostEqual(aligner.query_gap_score, score)
        self.assertAlmostEqual(aligner.query_open_gap_score, score)
        self.assertAlmostEqual(aligner.query_extend_gap_score, score)
        score = -4
        aligner.target_end_gap_score = score
        self.assertAlmostEqual(aligner.target_end_gap_score, score)
        self.assertAlmostEqual(aligner.target_end_open_gap_score, score)
        self.assertAlmostEqual(aligner.target_end_extend_gap_score, score)
        self.assertAlmostEqual(aligner.target_left_gap_score, score)
        self.assertAlmostEqual(aligner.target_left_open_gap_score, score)
        self.assertAlmostEqual(aligner.target_left_extend_gap_score, score)
        self.assertAlmostEqual(aligner.target_right_gap_score, score)
        self.assertAlmostEqual(aligner.target_right_open_gap_score, score)
        self.assertAlmostEqual(aligner.target_right_extend_gap_score, score)
        score = -5
        aligner.query_end_gap_score = score
        self.assertAlmostEqual(aligner.query_end_gap_score, score)
        self.assertAlmostEqual(aligner.query_end_open_gap_score, score)
        self.assertAlmostEqual(aligner.query_end_extend_gap_score, score)
        self.assertAlmostEqual(aligner.query_left_gap_score, score)
        self.assertAlmostEqual(aligner.query_left_open_gap_score, score)
        self.assertAlmostEqual(aligner.query_left_extend_gap_score, score)
        self.assertAlmostEqual(aligner.query_right_gap_score, score)
        self.assertAlmostEqual(aligner.query_right_open_gap_score, score)
        self.assertAlmostEqual(aligner.query_right_extend_gap_score, score)
        self.assertEqual(str(aligner), """\
Pairwise sequence aligner with parameters
  match score: 1.000000
  mismatch score: 0.000000
  target open gap score: -3.000000
  target extend gap score: -3.000000
  target left open gap score: -4.000000
  target left extend gap score: -4.000000
  target right open gap score: -4.000000
  target right extend gap score: -4.000000
  query open gap score: -2.000000
  query extend gap score: -2.000000
  query left open gap score: -5.000000
  query left extend gap score: -5.000000
  query right open gap score: -5.000000
  query right extend gap score: -5.000000
  mode: global
""")
        with self.assertRaises(ValueError):
            aligner.target_gap_score = 'wrong'
        with self.assertRaises(ValueError):
            aligner.query_gap_score = 'wrong'
        with self.assertRaises(TypeError):
            aligner.target_end_gap_score = 'wrong'
        with self.assertRaises(TypeError):
            aligner.query_end_gap_score = 'wrong'


class TestPairwiseGlobal(unittest.TestCase):

    def test_needlemanwunsch_simple1(self):
        seq1 = "GAACT"
        seq2 = "GAT"
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        self.assertEqual(str(aligner), """\
Pairwise sequence aligner with parameters
  match score: 1.000000
  mismatch score: 0.000000
  target open gap score: 0.000000
  target extend gap score: 0.000000
  target left open gap score: 0.000000
  target left extend gap score: 0.000000
  target right open gap score: 0.000000
  target right extend gap score: 0.000000
  query open gap score: 0.000000
  query extend gap score: 0.000000
  query left open gap score: 0.000000
  query left extend gap score: 0.000000
  query right open gap score: 0.000000
  query right extend gap score: 0.000000
  mode: global
""")
        self.assertEqual(aligner.algorithm, "Needleman-Wunsch")
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 3.0)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(str(alignment), """\
GAACT
||--|
GA--T
""")
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 3.0)
        self.assertEqual(str(alignment), """\
GAACT
|-|-|
G-A-T
""")

    def test_align_affine1_score(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aligner.match = 0
        aligner.mismatch = -1
        aligner.open_gap_score = -5
        aligner.extend_gap_score = -1
        self.assertEqual(aligner.algorithm, 'Gotoh global alignment algorithm')
        self.assertEqual(str(aligner), """\
Pairwise sequence aligner with parameters
  match score: 0.000000
  mismatch score: -1.000000
  target open gap score: -5.000000
  target extend gap score: -1.000000
  target left open gap score: -5.000000
  target left extend gap score: -1.000000
  target right open gap score: -5.000000
  target right extend gap score: -1.000000
  query open gap score: -5.000000
  query extend gap score: -1.000000
  query left open gap score: -5.000000
  query left extend gap score: -1.000000
  query right open gap score: -5.000000
  query right extend gap score: -1.000000
  mode: global
""")
        score = aligner.score("CC", "ACCT")
        self.assertAlmostEqual(score, -7.0)


class TestPairwiseLocal(unittest.TestCase):

    def test_smithwaterman(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'
        aligner.gap_score = -0.1
        self.assertEqual(aligner.algorithm, 'Smith-Waterman')
        self.assertEqual(str(aligner), """\
Pairwise sequence aligner with parameters
  match score: 1.000000
  mismatch score: 0.000000
  target open gap score: -0.100000
  target extend gap score: -0.100000
  target left open gap score: -0.100000
  target left extend gap score: -0.100000
  target right open gap score: -0.100000
  target right extend gap score: -0.100000
  query open gap score: -0.100000
  query extend gap score: -0.100000
  query left open gap score: -0.100000
  query left extend gap score: -0.100000
  query right open gap score: -0.100000
  query right extend gap score: -0.100000
  mode: local
""")
        score = aligner.score("AxBx", "zABz")
        self.assertAlmostEqual(score, 1.9)
        alignments = aligner.align("AxBx", "zABz")
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1.9)
        self.assertEqual(str(alignment), """\
.AxBx
.|-|.
zA-Bz
""")
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 1.9)
        self.assertEqual(str(alignment), """\
.AxBx
.|-|X
zA-Bz
""")

    def test_gotoh_local(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'
        aligner.open_gap_score = -0.1
        aligner.extend_gap_score = 0.0
        self.assertEqual(aligner.algorithm, 'Gotoh local alignment algorithm')
        self.assertEqual(str(aligner), """\
Pairwise sequence aligner with parameters
  match score: 1.000000
  mismatch score: 0.000000
  target open gap score: -0.100000
  target extend gap score: 0.000000
  target left open gap score: -0.100000
  target left extend gap score: 0.000000
  target right open gap score: -0.100000
  target right extend gap score: 0.000000
  query open gap score: -0.100000
  query extend gap score: 0.000000
  query left open gap score: -0.100000
  query left extend gap score: 0.000000
  query right open gap score: -0.100000
  query right extend gap score: 0.000000
  mode: local
""")
        score = aligner.score("AxBx", "zABz")
        self.assertAlmostEqual(score, 1.9)
        alignments = aligner.align("AxBx", "zABz")
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1.9)
        self.assertEqual(str(alignment), """\
.AxBx
.|-|.
zA-Bz
""")
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 1.9)
        self.assertEqual(str(alignment), """\
.AxBx
.|-|X
zA-Bz
""")


class TestPairwiseOpenPenalty(unittest.TestCase):

    def test_match_score_open_penalty1(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aligner.match = 2
        aligner.mismatch = -1
        aligner.open_gap_score = -0.1
        aligner.extend_gap_score = 0.0
        self.assertEqual(aligner.algorithm, 'Gotoh global alignment algorithm')
        self.assertEqual(str(aligner), """\
Pairwise sequence aligner with parameters
  match score: 2.000000
  mismatch score: -1.000000
  target open gap score: -0.100000
  target extend gap score: 0.000000
  target left open gap score: -0.100000
  target left extend gap score: 0.000000
  target right open gap score: -0.100000
  target right extend gap score: 0.000000
  query open gap score: -0.100000
  query extend gap score: 0.000000
  query left open gap score: -0.100000
  query left extend gap score: 0.000000
  query right open gap score: -0.100000
  query right extend gap score: 0.000000
  mode: global
""")
        seq1 = "AA"
        seq2 = "A"
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 1.9)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1.9)
        self.assertEqual(str(alignment), """\
AA
-|
-A
""")
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 1.9)
        self.assertEqual(str(alignment), """\
AA
|-
A-
""")

    def test_match_score_open_penalty2(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aligner.match = 1.5
        aligner.mismatch = 0.0
        aligner.open_gap_score = -0.1
        aligner.extend_gap_score = 0.0
        self.assertEqual(aligner.algorithm, 'Gotoh global alignment algorithm')
        self.assertEqual(str(aligner), """\
Pairwise sequence aligner with parameters
  match score: 1.500000
  mismatch score: 0.000000
  target open gap score: -0.100000
  target extend gap score: 0.000000
  target left open gap score: -0.100000
  target left extend gap score: 0.000000
  target right open gap score: -0.100000
  target right extend gap score: 0.000000
  query open gap score: -0.100000
  query extend gap score: 0.000000
  query left open gap score: -0.100000
  query left extend gap score: 0.000000
  query right open gap score: -0.100000
  query right extend gap score: 0.000000
  mode: global
""")
        seq1 = "GAA"
        seq2 = "GA"
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 2.9)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 2.9)
        self.assertEqual(str(alignment), """\
GAA
|-|
G-A
""")
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 2.9)
        self.assertEqual(str(alignment), """\
GAA
||-
GA-
""")

    def test_match_score_open_penalty3(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aligner.query_open_gap_score = -0.1
        aligner.query_extend_gap_score = 0.0
        self.assertEqual(aligner.algorithm, 'Gotoh global alignment algorithm')
        self.assertEqual(str(aligner), """\
Pairwise sequence aligner with parameters
  match score: 1.000000
  mismatch score: 0.000000
  target open gap score: 0.000000
  target extend gap score: 0.000000
  target left open gap score: 0.000000
  target left extend gap score: 0.000000
  target right open gap score: 0.000000
  target right extend gap score: 0.000000
  query open gap score: -0.100000
  query extend gap score: 0.000000
  query left open gap score: -0.100000
  query left extend gap score: 0.000000
  query right open gap score: -0.100000
  query right extend gap score: 0.000000
  mode: global
""")
        seq1 = "GAACT"
        seq2 = "GAT"
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 2.9)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 2.9)
        self.assertEqual(str(alignment), """\
GAACT
||--|
GA--T
""")

    def test_match_score_open_penalty4(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aligner.mismatch = -2.0
        aligner.open_gap_score = -0.1
        aligner.extend_gap_score = 0.0
        self.assertEqual(aligner.algorithm, 'Gotoh global alignment algorithm')
        self.assertEqual(str(aligner), """\
Pairwise sequence aligner with parameters
  match score: 1.000000
  mismatch score: -2.000000
  target open gap score: -0.100000
  target extend gap score: 0.000000
  target left open gap score: -0.100000
  target left extend gap score: 0.000000
  target right open gap score: -0.100000
  target right extend gap score: 0.000000
  query open gap score: -0.100000
  query extend gap score: 0.000000
  query left open gap score: -0.100000
  query left extend gap score: 0.000000
  query right open gap score: -0.100000
  query right extend gap score: 0.000000
  mode: global
""")
        seq1 = "GCT"
        seq2 = "GATA"
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 1.7)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1.7)
        self.assertEqual(str(alignment), """\
G-CT-
|--|-
GA-TA
""")
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 1.7)
        self.assertEqual(str(alignment), """\
GC-T-
|--|-
G-ATA
""")


class TestPairwiseExtendPenalty(unittest.TestCase):

    def test_extend_penalty1(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aligner.open_gap_score = -0.2
        aligner.extend_gap_score = -0.5
        self.assertEqual(aligner.algorithm, 'Gotoh global alignment algorithm')
        self.assertEqual(str(aligner), """\
Pairwise sequence aligner with parameters
  match score: 1.000000
  mismatch score: 0.000000
  target open gap score: -0.200000
  target extend gap score: -0.500000
  target left open gap score: -0.200000
  target left extend gap score: -0.500000
  target right open gap score: -0.200000
  target right extend gap score: -0.500000
  query open gap score: -0.200000
  query extend gap score: -0.500000
  query left open gap score: -0.200000
  query left extend gap score: -0.500000
  query right open gap score: -0.200000
  query right extend gap score: -0.500000
  mode: global
""")
        seq1 = "GACT"
        seq2 = "GT"
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 1.3)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1.3)
        self.assertEqual(str(alignment), """\
GACT
|--|
G--T
""")

    def test_extend_penalty2(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aligner.open_gap_score = -0.2
        aligner.extend_gap_score = -1.5
        self.assertEqual(aligner.algorithm, 'Gotoh global alignment algorithm')
        self.assertEqual(str(aligner), """\
Pairwise sequence aligner with parameters
  match score: 1.000000
  mismatch score: 0.000000
  target open gap score: -0.200000
  target extend gap score: -1.500000
  target left open gap score: -0.200000
  target left extend gap score: -1.500000
  target right open gap score: -0.200000
  target right extend gap score: -1.500000
  query open gap score: -0.200000
  query extend gap score: -1.500000
  query left open gap score: -0.200000
  query left extend gap score: -1.500000
  query right open gap score: -0.200000
  query right extend gap score: -1.500000
  mode: global
""")
        seq1 = "GACT"
        seq2 = "GT"
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 0.6)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 0.6)
        self.assertEqual(str(alignment), """\
GACT
-X-|
-G-T
""")
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 0.6)
        self.assertEqual(str(alignment), """\
GACT
|-X-
G-T-
""")


class TestPairwisePenalizeExtendWhenOpening(unittest.TestCase):

    def test_penalize_extend_when_opening(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aligner.open_gap_score = -1.7
        aligner.extend_gap_score = -1.5
        self.assertEqual(aligner.algorithm, 'Gotoh global alignment algorithm')
        self.assertEqual(str(aligner), """\
Pairwise sequence aligner with parameters
  match score: 1.000000
  mismatch score: 0.000000
  target open gap score: -1.700000
  target extend gap score: -1.500000
  target left open gap score: -1.700000
  target left extend gap score: -1.500000
  target right open gap score: -1.700000
  target right extend gap score: -1.500000
  query open gap score: -1.700000
  query extend gap score: -1.500000
  query left open gap score: -1.700000
  query left extend gap score: -1.500000
  query right open gap score: -1.700000
  query right extend gap score: -1.500000
  mode: global
""")
        seq1 = "GACT"
        seq2 = "GT"
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, -1.2)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, -1.2)
        self.assertEqual(str(alignment), """\
GACT
|--|
G--T
""")


class TestPairwisePenalizeEndgaps(unittest.TestCase):

    def test_penalize_end_gaps(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aligner.open_gap_score = -0.2
        aligner.extend_gap_score = -0.8
        end_score = 0.0
        aligner.target_end_gap_score = end_score
        aligner.query_end_gap_score = end_score
        self.assertEqual(str(aligner), """\
Pairwise sequence aligner with parameters
  match score: 1.000000
  mismatch score: 0.000000
  target open gap score: -0.200000
  target extend gap score: -0.800000
  target left open gap score: 0.000000
  target left extend gap score: 0.000000
  target right open gap score: 0.000000
  target right extend gap score: 0.000000
  query open gap score: -0.200000
  query extend gap score: -0.800000
  query left open gap score: 0.000000
  query left extend gap score: 0.000000
  query right open gap score: 0.000000
  query right extend gap score: 0.000000
  mode: global
""")
        self.assertEqual(aligner.algorithm, 'Gotoh global alignment algorithm')
        seq1 = "GACT"
        seq2 = "GT"
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 1.0)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(len(alignments), 3)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1.0)
        self.assertEqual(str(alignment), """\
GACT
--X|
--GT
""")
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 1.0)
        self.assertEqual(str(alignment), """\
GACT
|--|
G--T
""")
        alignment = alignments[2]
        self.assertAlmostEqual(alignment.score, 1.0)
        self.assertEqual(str(alignment), """\
GACT
|X--
GT--
""")


class TestPairwiseSeparateGapPenalties(unittest.TestCase):

    def test_separate_gap_penalties1(self):
        seq1 = "GAT"
        seq2 = "GTCT"
        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'
        open_score, extend_score = (-0.3, 0)
        aligner.target_open_gap_score = open_score
        aligner.target_extend_gap_score = extend_score
        aligner.target_end_open_gap_score = open_score
        aligner.target_end_extend_gap_score = extend_score
        open_score, extend_score = (-0.8, 0)
        aligner.query_open_gap_score = open_score
        aligner.query_extend_gap_score = extend_score
        aligner.query_end_open_gap_score = open_score
        aligner.query_end_extend_gap_score = extend_score
        self.assertEqual(aligner.algorithm, 'Gotoh local alignment algorithm')
        self.assertEqual(str(aligner), """\
Pairwise sequence aligner with parameters
  match score: 1.000000
  mismatch score: 0.000000
  target open gap score: -0.300000
  target extend gap score: 0.000000
  target left open gap score: -0.300000
  target left extend gap score: 0.000000
  target right open gap score: -0.300000
  target right extend gap score: 0.000000
  query open gap score: -0.800000
  query extend gap score: 0.000000
  query left open gap score: -0.800000
  query left extend gap score: 0.000000
  query right open gap score: -0.800000
  query right extend gap score: 0.000000
  mode: local
""")
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 1.7)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1.7)
        self.assertEqual(str(alignment), """\
G-AT
|-X|
GTCT
""")
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 1.7)
        self.assertEqual(str(alignment), """\
GA-T
|X-|
GTCT
""")

    def test_separate_gap_penalties2(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'
        aligner.target_open_gap_score = -0.3
        aligner.target_extend_gap_score = 0.0
        aligner.query_open_gap_score = -0.2
        aligner.query_extend_gap_score = 0.0
        self.assertEqual(aligner.algorithm, "Gotoh local alignment algorithm")
        self.assertEqual(str(aligner), """\
Pairwise sequence aligner with parameters
  match score: 1.000000
  mismatch score: 0.000000
  target open gap score: -0.300000
  target extend gap score: 0.000000
  target left open gap score: -0.300000
  target left extend gap score: 0.000000
  target right open gap score: -0.300000
  target right extend gap score: 0.000000
  query open gap score: -0.200000
  query extend gap score: 0.000000
  query left open gap score: -0.200000
  query left extend gap score: 0.000000
  query right open gap score: -0.200000
  query right extend gap score: 0.000000
  mode: local
""")
        seq1 = "GAT"
        seq2 = "GTCT"
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 1.8)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1.8)
        self.assertEqual(str(alignment), """\
GAT..
|-|..
G-TCT
""")


class TestPairwiseSeparateGapPenaltiesWithExtension(unittest.TestCase):

    def test_separate_gap_penalties_with_extension(self):
        seq1 = "GAAT"
        seq2 = "GTCCT"
        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'
        open_score, extend_score = (-0.1, 0)
        aligner.target_open_gap_score = open_score
        aligner.target_extend_gap_score = extend_score
        aligner.target_end_open_gap_score = open_score
        aligner.target_end_extend_gap_score = extend_score
        score = -0.1
        aligner.query_gap_score = score
        aligner.query_end_gap_score = score
        self.assertEqual(aligner.algorithm, "Gotoh local alignment algorithm")
        self.assertEqual(str(aligner), """\
Pairwise sequence aligner with parameters
  match score: 1.000000
  mismatch score: 0.000000
  target open gap score: -0.100000
  target extend gap score: 0.000000
  target left open gap score: -0.100000
  target left extend gap score: 0.000000
  target right open gap score: -0.100000
  target right extend gap score: 0.000000
  query open gap score: -0.100000
  query extend gap score: -0.100000
  query left open gap score: -0.100000
  query left extend gap score: -0.100000
  query right open gap score: -0.100000
  query right extend gap score: -0.100000
  mode: local
""")
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 1.9)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(len(alignments), 3)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1.9)
        self.assertEqual(str(alignment), """\
G-AAT
|-XX|
GTCCT
""")
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 1.9)
        self.assertEqual(str(alignment), """\
GA-AT
|X-X|
GTCCT
""")
        alignment = alignments[2]
        self.assertAlmostEqual(alignment.score, 1.9)
        self.assertEqual(str(alignment), """\
GAA-T
|XX-|
GTCCT
""")
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(len(alignments), 3)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1.9)
        self.assertEqual(str(alignment), """\
G-AAT
|-XX|
GTCCT
""")
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 1.9)
        self.assertEqual(str(alignment), """\
GA-AT
|X-X|
GTCCT
""")
        alignment = alignments[2]
        self.assertAlmostEqual(alignment.score, 1.9)
        self.assertEqual(str(alignment), """\
GAA-T
|XX-|
GTCCT
""")


class TestPairwiseMatchDictionary(unittest.TestCase):

    match_dict = {
        ("A", "A"): 1.5,
        ("A", "T"): 0.5,
        ("T", "A"): 0.5,
        ("T", "T"): 1.0
        }

    def test_match_dictionary1(self):
        seq1 = "ATAT"
        seq2 = "ATT"
        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'
        aligner.substitution_matrix = self.match_dict
        aligner.open_gap_score = -0.5
        aligner.extend_gap_score = 0.0
        self.assertEqual(aligner.algorithm, "Gotoh local alignment algorithm")
        self.assertEqual(str(aligner), """\
Pairwise sequence aligner with parameters
  match/mismatch score: <substitution matrix>
  target open gap score: -0.500000
  target extend gap score: 0.000000
  target left open gap score: -0.500000
  target left extend gap score: 0.000000
  target right open gap score: -0.500000
  target right extend gap score: 0.000000
  query open gap score: -0.500000
  query extend gap score: 0.000000
  query left open gap score: -0.500000
  query left extend gap score: 0.000000
  query right open gap score: -0.500000
  query right extend gap score: 0.000000
  mode: local
""")
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 3.0)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(len(alignments), 2)
        self.assertAlmostEqual(alignments[0].score, 3.0)
        self.assertEqual(str(alignments[0]), """\
ATAT
||X.
ATT.
""")
        self.assertAlmostEqual(alignments[1].score, 3.0)
        self.assertEqual(str(alignments[1]), """\
ATAT
||-|
AT-T
""")

    def test_match_dictionary2(self):
        seq1 = "ATAT"
        seq2 = "ATT"
        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'
        aligner.substitution_matrix = self.match_dict
        aligner.open_gap_score = -1.0
        aligner.extend_gap_score = 0.0
        self.assertEqual(str(aligner), """\
Pairwise sequence aligner with parameters
  match/mismatch score: <substitution matrix>
  target open gap score: -1.000000
  target extend gap score: 0.000000
  target left open gap score: -1.000000
  target left extend gap score: 0.000000
  target right open gap score: -1.000000
  target right extend gap score: 0.000000
  query open gap score: -1.000000
  query extend gap score: 0.000000
  query left open gap score: -1.000000
  query left extend gap score: 0.000000
  query right open gap score: -1.000000
  query right extend gap score: 0.000000
  mode: local
""")
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 3.0)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(len(alignments), 1)
        self.assertAlmostEqual(alignments[0].score, 3.0)
        self.assertEqual(str(alignments[0]), """\
ATAT
||X.
ATT.
""")

    def test_match_dictionary3(self):
        seq1 = "ATT"
        seq2 = "ATAT"
        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'
        aligner.substitution_matrix = self.match_dict
        aligner.open_gap_score = -1.0
        aligner.extend_gap_score = 0.0
        self.assertEqual(str(aligner), """\
Pairwise sequence aligner with parameters
  match/mismatch score: <substitution matrix>
  target open gap score: -1.000000
  target extend gap score: 0.000000
  target left open gap score: -1.000000
  target left extend gap score: 0.000000
  target right open gap score: -1.000000
  target right extend gap score: 0.000000
  query open gap score: -1.000000
  query extend gap score: 0.000000
  query left open gap score: -1.000000
  query left extend gap score: 0.000000
  query right open gap score: -1.000000
  query right extend gap score: 0.000000
  mode: local
""")
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 3.0)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(len(alignments), 1)
        self.assertAlmostEqual(alignments[0].score, 3.0)
        self.assertEqual(str(alignments[0]), """\
ATT.
||X.
ATAT
""")


class TestPairwiseOneCharacter(unittest.TestCase):

    def test_align_one_char1(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'
        aligner.open_gap_score = -0.3
        aligner.extend_gap_score = -0.1
        self.assertEqual(aligner.algorithm, 'Gotoh local alignment algorithm')
        self.assertEqual(str(aligner), """\
Pairwise sequence aligner with parameters
  match score: 1.000000
  mismatch score: 0.000000
  target open gap score: -0.300000
  target extend gap score: -0.100000
  target left open gap score: -0.300000
  target left extend gap score: -0.100000
  target right open gap score: -0.300000
  target right extend gap score: -0.100000
  query open gap score: -0.300000
  query extend gap score: -0.100000
  query left open gap score: -0.300000
  query left extend gap score: -0.100000
  query right open gap score: -0.300000
  query right extend gap score: -0.100000
  mode: local
""")
        score = aligner.score("abcde", "c")
        self.assertAlmostEqual(score, 1)
        alignments = aligner.align("abcde", "c")
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1)
        self.assertEqual(str(alignment), """\
abcde
..|..
..c..
""")

    def test_align_one_char2(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'
        aligner.open_gap_score = -0.3
        aligner.extend_gap_score = -0.1
        self.assertEqual(aligner.algorithm, 'Gotoh local alignment algorithm')
        self.assertEqual(str(aligner), """\
Pairwise sequence aligner with parameters
  match score: 1.000000
  mismatch score: 0.000000
  target open gap score: -0.300000
  target extend gap score: -0.100000
  target left open gap score: -0.300000
  target left extend gap score: -0.100000
  target right open gap score: -0.300000
  target right extend gap score: -0.100000
  query open gap score: -0.300000
  query extend gap score: -0.100000
  query left open gap score: -0.300000
  query left extend gap score: -0.100000
  query right open gap score: -0.300000
  query right extend gap score: -0.100000
  mode: local
""")
        score = aligner.score("abcce", "c")
        self.assertAlmostEqual(score, 1)
        alignments = aligner.align("abcce", "c")
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 1)
        self.assertEqual(str(alignment), """\
abcce
..|..
..c..
""")
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 1)
        self.assertEqual(str(alignment), """\
abcce
...|.
...c.
""")

    def test_align_one_char3(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aligner.open_gap_score = -0.3
        aligner.extend_gap_score = -0.1
        self.assertEqual(aligner.algorithm, 'Gotoh global alignment algorithm')
        self.assertEqual(str(aligner), """\
Pairwise sequence aligner with parameters
  match score: 1.000000
  mismatch score: 0.000000
  target open gap score: -0.300000
  target extend gap score: -0.100000
  target left open gap score: -0.300000
  target left extend gap score: -0.100000
  target right open gap score: -0.300000
  target right extend gap score: -0.100000
  query open gap score: -0.300000
  query extend gap score: -0.100000
  query left open gap score: -0.300000
  query left extend gap score: -0.100000
  query right open gap score: -0.300000
  query right extend gap score: -0.100000
  mode: global
""")
        seq1 = "abcde"
        seq2 = "c"
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 0.2)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 0.2)
        self.assertEqual(str(alignment), """\
abcde
--|--
--c--
""")

    def test_align_one_char_score3(self):
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aligner.open_gap_score = -0.3
        aligner.extend_gap_score = -0.1
        self.assertEqual(aligner.algorithm, 'Gotoh global alignment algorithm')
        self.assertEqual(str(aligner), """\
Pairwise sequence aligner with parameters
  match score: 1.000000
  mismatch score: 0.000000
  target open gap score: -0.300000
  target extend gap score: -0.100000
  target left open gap score: -0.300000
  target left extend gap score: -0.100000
  target right open gap score: -0.300000
  target right extend gap score: -0.100000
  query open gap score: -0.300000
  query extend gap score: -0.100000
  query left open gap score: -0.300000
  query left extend gap score: -0.100000
  query right open gap score: -0.300000
  query right extend gap score: -0.100000
  mode: global
""")
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
        nogaps = lambda x, y: -2000 - y
        # Very expensive to open a gap in seq2 unless it is in one of the allowed positions
        specificgaps = lambda x, y: (-2 - y) if x in breaks else (-2000 - y)
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aligner.match = 1
        aligner.mismatch = -1
        aligner.target_gap_score = nogaps
        aligner.query_gap_score = specificgaps
        self.assertEqual(str(aligner), """\
Pairwise sequence aligner with parameters
  match score: 1.000000
  mismatch score: -1.000000
  target gap function: %s
  query gap function: %s
  mode: global
""" % (nogaps, specificgaps))
        self.assertEqual(aligner.algorithm,
                         "Waterman-Smith-Beyer global alignment algorithm")
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 2)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 2)
        self.assertEqual(str(alignment), """\
AAAABBBAAAACCCCCCCCCCCCCCAAAABBBAAAA
--|||||||||||----------|||||||||||--
--AABBBAAAACC----------CCAAAABBBAA--
""")

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
        nogaps = lambda x, y: -2000 - y
        # Very expensive to open a gap in seq2 unless it is in one of the allowed positions:
        specificgaps = lambda x, y: (-2 - y) if x in breaks else (-2000 - y)
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aligner.match = 1
        aligner.mismatch = -1
        aligner.target_gap_score = nogaps
        aligner.query_gap_score = specificgaps
        self.assertEqual(str(aligner), """\
Pairwise sequence aligner with parameters
  match score: 1.000000
  mismatch score: -1.000000
  target gap function: %s
  query gap function: %s
  mode: global
""" % (nogaps, specificgaps))
        self.assertEqual(aligner.algorithm,
                         "Waterman-Smith-Beyer global alignment algorithm")
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, -10)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, -10)
        self.assertEqual(str(alignment), """\
AAAABBBAAAACCCCCCCCCCCCCCAAAABBBAAAA
--|||----------XXXXXX|||||||||||||--
--AAB----------BBAAAACCCCAAAABBBAA--
""")
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, -10)
        self.assertEqual(str(alignment), """\
AAAABBBAAAACCCCCCCCCCCCCCAAAABBBAAAA
||X------------XXXXXX|||||||||||||--
AAB------------BBAAAACCCCAAAABBBAA--
""")

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
        aligner.mode = 'global'
        aligner.match = 1
        aligner.mismatch = -10
        aligner.target_gap_score = gap_score
        self.assertEqual(aligner.algorithm,
                         "Waterman-Smith-Beyer global alignment algorithm")
        self.assertEqual(str(aligner), """\
Pairwise sequence aligner with parameters
  match score: 1.000000
  mismatch score: -10.000000
  target gap function: %s
  query open gap score: 0.000000
  query extend gap score: 0.000000
  query left open gap score: 0.000000
  query left extend gap score: 0.000000
  query right open gap score: 0.000000
  query right extend gap score: 0.000000
  mode: global
""" % gap_score)
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 2.0)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 2.0)
        self.assertEqual(str(alignment), """\
TT-CC-AA
||----||
TTG--GAA
""")
        aligner.query_gap_score = gap_score
        self.assertEqual(str(aligner), """\
Pairwise sequence aligner with parameters
  match score: 1.000000
  mismatch score: -10.000000
  target gap function: %s
  query gap function: %s
  mode: global
""" % (gap_score, gap_score))
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, -8.0)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(len(alignments), 4)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, -8.0)
        self.assertEqual(str(alignment), """\
TT-CCAA
||-X-||
TTGG-AA
""")
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, -8.0)
        self.assertEqual(str(alignment), """\
TTC--CAA
||----||
TT-GG-AA
""")
        alignment = alignments[2]
        self.assertAlmostEqual(alignment.score, -8.0)
        self.assertEqual(str(alignment), """\
TTCC-AA
||-X-||
TT-GGAA
""")
        alignment = alignments[3]
        self.assertAlmostEqual(alignment.score, -8.0)
        self.assertEqual(str(alignment), """\
TT-CC-AA
||----||
TTG--GAA
""")

    def test_gap_here_only_local_1(self):
        seq1 = "AAAABBBAAAACCCCCCCCCCCCCCAAAABBBAAAA"
        seq2 = "AABBBAAAACCCCAAAABBBAA"
        breaks = [0, 11, len(seq2)]
        # Very expensive to open a gap in seq1:
        nogaps = lambda x, y: -2000 - y
        # Very expensive to open a gap in seq2 unless it is in one of the allowed positions
        specificgaps = lambda x, y: (-2 - y) if x in breaks else (-2000 - y)
        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'
        aligner.match = 1
        aligner.mismatch = -1
        aligner.target_gap_score = nogaps
        aligner.query_gap_score = specificgaps
        self.assertEqual(aligner.algorithm,
                         "Waterman-Smith-Beyer local alignment algorithm")
        self.assertEqual(str(aligner), """\
Pairwise sequence aligner with parameters
  match score: 1.000000
  mismatch score: -1.000000
  target gap function: %s
  query gap function: %s
  mode: local
""" % (nogaps, specificgaps))
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 13)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 13)
        self.assertEqual(str(alignment), """\
AAAABBBAAAACCCCCCCCCCCCCCAAAABBBAAAA
..|||||||||||||.....................
..AABBBAAAACCCCAAAABBBAA............
""")
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 13)
        self.assertEqual(str(alignment), """\
AAAABBBAAAACCCCCCCCCCCCCCAAAABBBAAAA
.....................|||||||||||||..
............AABBBAAAACCCCAAAABBBAA..
""")

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
        nogaps = lambda x, y: -2000 - y
        # Very expensive to open a gap in seq2 unless it is in one of the allowed positions:
        specificgaps = lambda x, y: (-2 - y) if x in breaks else (-2000 - y)
        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'
        aligner.match = 1
        aligner.mismatch = -1
        aligner.target_gap_score = nogaps
        aligner.query_gap_score = specificgaps
        self.assertEqual(str(aligner), """\
Pairwise sequence aligner with parameters
  match score: 1.000000
  mismatch score: -1.000000
  target gap function: %s
  query gap function: %s
  mode: local
""" % (nogaps, specificgaps))
        self.assertEqual(aligner.algorithm,
                         "Waterman-Smith-Beyer local alignment algorithm")
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 13)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 13)
        self.assertEqual(str(alignment), """\
AAAABBBAAAACCCCCCCCCCCCCCAAAABBBAAAA
..|||||||||||||.....................
..AABBBAAAACCCCAAAABBBAA............
""")
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 13)
        self.assertEqual(str(alignment), """\
AAAABBBAAAACCCCCCCCCCCCCCAAAABBBAAAA
.....................|||||||||||||..
............AABBBAAAACCCCAAAABBBAA..
""")

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
        aligner.mode = 'local'
        aligner.match = 1
        aligner.mismatch = -10
        aligner.target_gap_score = gap_score
        self.assertEqual(aligner.algorithm,
                         "Waterman-Smith-Beyer local alignment algorithm")
        self.assertEqual(str(aligner), """\
Pairwise sequence aligner with parameters
  match score: 1.000000
  mismatch score: -10.000000
  target gap function: %s
  query open gap score: 0.000000
  query extend gap score: 0.000000
  query left open gap score: 0.000000
  query left extend gap score: 0.000000
  query right open gap score: 0.000000
  query right extend gap score: 0.000000
  mode: local
""" % gap_score)
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 2.0)
        alignments = aligner.align(seq1, seq2)
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 2.0)
        self.assertEqual(str(alignment), """\
TTCCAA
||....
TTGGAA
""")
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 2.0)
        self.assertEqual(str(alignment), """\
TTCCAA
....||
TTGGAA
""")
        aligner.query_gap = gap_score
        self.assertEqual(str(aligner), """\
Pairwise sequence aligner with parameters
  match score: 1.000000
  mismatch score: -10.000000
  target gap function: %s
  query open gap score: 0.000000
  query extend gap score: 0.000000
  query left open gap score: 0.000000
  query left extend gap score: 0.000000
  query right open gap score: 0.000000
  query right extend gap score: 0.000000
  mode: local
""" % gap_score)
        alignments = aligner.align(seq1, seq2)
        score = aligner.score(seq1, seq2)
        self.assertAlmostEqual(score, 2.0)
        self.assertEqual(len(alignments), 2)
        alignment = alignments[0]
        self.assertAlmostEqual(alignment.score, 2.0)
        self.assertEqual(str(alignment), """\
TTCCAA
||....
TTGGAA
""")
        alignment = alignments[1]
        self.assertAlmostEqual(alignment.score, 2.0)
        self.assertEqual(str(alignment), """\
TTCCAA
....||
TTGGAA
""")

    def test_broken_gap_function(self):
        # Check if an Exception is propagated if the gap function raises one
        seq1 = "TTCCAA"
        seq2 = "TTGGAA"

        def gap_score(i, n):
            raise RuntimeError('broken gap function')
        aligner = Align.PairwiseAligner()
        aligner.target_gap_score = gap_score
        aligner.query_gap_score = -1
        aligner.mode = 'global'
        with self.assertRaises(RuntimeError):
            aligner.score(seq1, seq2)
        with self.assertRaises(RuntimeError):
            alignments = aligner.align(seq1, seq2)
            alignments = list(alignments)
        aligner.mode = 'local'
        with self.assertRaises(RuntimeError):
            aligner.score(seq1, seq2)
        with self.assertRaises(RuntimeError):
            alignments = aligner.align(seq1, seq2)
            alignments = list(alignments)
        aligner.target_gap_score = -1
        aligner.query_gap_score = gap_score
        aligner.mode = 'global'
        with self.assertRaises(RuntimeError):
            aligner.score(seq1, seq2)
        with self.assertRaises(RuntimeError):
            alignments = aligner.align(seq1, seq2)
            alignments = list(alignments)
        aligner.mode = 'local'
        with self.assertRaises(RuntimeError):
            aligner.score(seq1, seq2)
        with self.assertRaises(RuntimeError):
            alignments = aligner.align(seq1, seq2)
            alignments = list(alignments)


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
