# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SubsMat module."""

try:
    import numpy
    del numpy
except ImportError:
    from Bio import MissingExternalDependencyError
    raise MissingExternalDependencyError(
        "Install NumPy if you want to use Bio.SubsMat.")

try:
    import cPickle as pickle  # Only available on Python 2
except ImportError:
    import pickle

import os
import unittest


from Bio._py3k import StringIO

from Bio import SubsMat
from Bio.SubsMat import FreqTable, MatrixInfo


class TestGeo(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        ftab_file = os.path.join('SubsMat', 'protein_count.txt')
        with open(ftab_file) as handle:
            cls.ftab_prot = FreqTable.read_count(handle)
        pickle_file = os.path.join('SubsMat', 'acc_rep_mat.pik')
        # Don't want to use text mode on Python 3,
        with open(pickle_file, 'rb') as handle:
            dictionary = pickle.load(handle)
        cls.acc_rep_mat = SubsMat.AcceptedReplacementsMatrix(dictionary)
        cls.obs_freq_mat = SubsMat._build_obs_freq_mat(cls.acc_rep_mat)

    def checkMatrix(self, mat, expected):
        handle = StringIO()
        mat.print_mat(f=handle)
        text = handle.getvalue()
        self.assertEqual(text, expected)

    def test_protein_freq(self):
        letters = "ACDEFGHIKLMNPQRSTVWY"
        ctab_file = os.path.join('SubsMat', 'protein_freq.txt')
        with open(ctab_file) as handle:
            ctab_prot = FreqTable.read_freq(handle)
        self.assertEqual(self.ftab_prot.alphabet.letters, letters)
        for letter in letters:
            difference = self.ftab_prot[letter] - ctab_prot[letter]
            self.assertAlmostEqual(abs(difference), 0, places=3)

    def test_freqtable(self):
        ftab_prot2 = SubsMat._exp_freq_table_from_obs_freq(self.obs_freq_mat)
        self.assertEqual(len(ftab_prot2), 20)
        for letter in 'ACDEFGHIKLMNPQRSTVWY':
            difference = self.ftab_prot[letter] - ftab_prot2[letter]
            self.assertAlmostEqual(abs(difference), 0, places=1)

    def test_obs_freq_mat(self):
        self.assertEqual(len(self.obs_freq_mat), 210)
        self.assertAlmostEqual(self.obs_freq_mat[('A', 'A')], 0.010561, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('A', 'C')], 0.002452, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('C', 'C')], 0.002102, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('A', 'D')], 0.008559, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('C', 'D')], 0.000751, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('D', 'D')], 0.007558, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('A', 'E')], 0.010360, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('C', 'E')], 0.001301, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('D', 'E')], 0.008158, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('E', 'E')], 0.005956, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('A', 'F')], 0.006106, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('C', 'F')], 0.001602, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('D', 'F')], 0.003153, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('E', 'F')], 0.002503, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('F', 'F')], 0.003954, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('A', 'G')], 0.011461, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('C', 'G')], 0.001952, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('D', 'G')], 0.007157, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('E', 'G')], 0.007307, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('F', 'G')], 0.005105, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('G', 'G')], 0.016066, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('A', 'H')], 0.003303, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('C', 'H')], 0.000501, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('D', 'H')], 0.002603, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('E', 'H')], 0.001902, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('F', 'H')], 0.001952, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('G', 'H')], 0.002653, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('H', 'H')], 0.002803, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('A', 'I')], 0.010110, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('C', 'I')], 0.002102, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('D', 'I')], 0.003303, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('E', 'I')], 0.003153, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('F', 'I')], 0.007558, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('G', 'I')], 0.004605, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('H', 'I')], 0.001552, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('I', 'I')], 0.007257, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('A', 'K')], 0.009560, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('C', 'K')], 0.001201, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('D', 'K')], 0.008458, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('E', 'K')], 0.008208, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('F', 'K')], 0.003353, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('G', 'K')], 0.007407, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('H', 'K')], 0.003153, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('I', 'K')], 0.004004, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('K', 'K')], 0.006056, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('A', 'L')], 0.015365, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('C', 'L')], 0.003103, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('D', 'L')], 0.006006, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('E', 'L')], 0.006657, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('F', 'L')], 0.011812, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('G', 'L')], 0.007708, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('H', 'L')], 0.003504, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('I', 'L')], 0.017417, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('K', 'L')], 0.005606, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('L', 'L')], 0.016016, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('A', 'M')], 0.003954, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('C', 'M')], 0.000701, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('D', 'M')], 0.001001, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('E', 'M')], 0.001652, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('F', 'M')], 0.003353, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('G', 'M')], 0.002452, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('H', 'M')], 0.000801, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('I', 'M')], 0.004404, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('K', 'M')], 0.001602, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('L', 'M')], 0.007658, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('M', 'M')], 0.001101, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('A', 'N')], 0.006957, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('C', 'N')], 0.001201, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('D', 'N')], 0.008158, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('E', 'N')], 0.005756, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('F', 'N')], 0.002753, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('G', 'N')], 0.007357, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('H', 'N')], 0.002853, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('I', 'N')], 0.003253, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('K', 'N')], 0.006957, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('L', 'N')], 0.005055, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('M', 'N')], 0.001702, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('N', 'N')], 0.003303, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('A', 'P')], 0.005956, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('C', 'P')], 0.001051, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('D', 'P')], 0.005255, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('E', 'P')], 0.004555, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('F', 'P')], 0.002503, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('G', 'P')], 0.005405, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('H', 'P')], 0.001151, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('I', 'P')], 0.003353, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('K', 'P')], 0.004705, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('L', 'P')], 0.003954, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('M', 'P')], 0.001001, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('N', 'P')], 0.004705, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('P', 'P')], 0.006006, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('A', 'Q')], 0.006456, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('C', 'Q')], 0.000801, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('D', 'Q')], 0.005806, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('E', 'Q')], 0.005255, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('F', 'Q')], 0.001952, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('G', 'Q')], 0.004454, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('H', 'Q')], 0.002202, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('I', 'Q')], 0.003103, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('K', 'Q')], 0.005405, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('L', 'Q')], 0.003954, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('M', 'Q')], 0.001351, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('N', 'Q')], 0.003854, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('P', 'Q')], 0.003403, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('Q', 'Q')], 0.002102, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('A', 'R')], 0.007608, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('C', 'R')], 0.000751, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('D', 'R')], 0.005756, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('E', 'R')], 0.005906, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('F', 'R')], 0.002853, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('G', 'R')], 0.004054, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('H', 'R')], 0.002202, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('I', 'R')], 0.002903, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('K', 'R')], 0.007558, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('L', 'R')], 0.005455, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('M', 'R')], 0.002002, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('N', 'R')], 0.004505, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('P', 'R')], 0.003604, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('Q', 'R')], 0.003654, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('R', 'R')], 0.004805, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('A', 'S')], 0.009960, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('C', 'S')], 0.001151, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('D', 'S')], 0.008358, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('E', 'S')], 0.007157, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('F', 'S')], 0.004605, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('G', 'S')], 0.008058, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('H', 'S')], 0.002953, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('I', 'S')], 0.004204, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('K', 'S')], 0.007207, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('L', 'S')], 0.007958, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('M', 'S')], 0.002052, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('N', 'S')], 0.005756, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('P', 'S')], 0.005155, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('Q', 'S')], 0.005055, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('R', 'S')], 0.004705, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('S', 'S')], 0.005806, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('A', 'T')], 0.009760, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('C', 'T')], 0.001502, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('D', 'T')], 0.006256, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('E', 'T')], 0.007508, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('F', 'T')], 0.003453, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('G', 'T')], 0.007658, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('H', 'T')], 0.003453, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('I', 'T')], 0.004404, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('K', 'T')], 0.007007, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('L', 'T')], 0.009159, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('M', 'T')], 0.001752, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('N', 'T')], 0.005255, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('P', 'T')], 0.004605, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('Q', 'T')], 0.004855, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('R', 'T')], 0.005706, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('S', 'T')], 0.008609, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('T', 'T')], 0.004204, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('A', 'V')], 0.012613, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('C', 'V')], 0.002953, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('D', 'V')], 0.004755, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('E', 'V')], 0.005506, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('F', 'V')], 0.008008, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('G', 'V')], 0.006557, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('H', 'V')], 0.002703, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('I', 'V')], 0.015516, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('K', 'V')], 0.006557, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('L', 'V')], 0.022372, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('M', 'V')], 0.004755, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('N', 'V')], 0.004655, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('P', 'V')], 0.005455, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('Q', 'V')], 0.003453, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('R', 'V')], 0.004404, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('S', 'V')], 0.005706, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('T', 'V')], 0.007257, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('V', 'V')], 0.011962, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('A', 'W')], 0.002352, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('C', 'W')], 0.000400, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('D', 'W')], 0.001101, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('E', 'W')], 0.001051, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('F', 'W')], 0.002152, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('G', 'W')], 0.001752, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('H', 'W')], 0.000651, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('I', 'W')], 0.001902, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('K', 'W')], 0.001051, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('L', 'W')], 0.003704, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('M', 'W')], 0.001251, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('N', 'W')], 0.001451, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('P', 'W')], 0.000901, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('Q', 'W')], 0.000751, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('R', 'W')], 0.001151, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('S', 'W')], 0.001101, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('T', 'W')], 0.001752, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('V', 'W')], 0.002202, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('W', 'W')], 0.001602, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('A', 'Y')], 0.006206, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('C', 'Y')], 0.000951, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('D', 'Y')], 0.003203, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('E', 'Y')], 0.003003, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('F', 'Y')], 0.007107, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('G', 'Y')], 0.005055, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('H', 'Y')], 0.002152, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('I', 'Y')], 0.004705, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('K', 'Y')], 0.003303, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('L', 'Y')], 0.008659, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('M', 'Y')], 0.002052, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('N', 'Y')], 0.003754, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('P', 'Y')], 0.003003, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('Q', 'Y')], 0.002052, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('R', 'Y')], 0.003103, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('S', 'Y')], 0.003654, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('T', 'Y')], 0.004555, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('V', 'Y')], 0.006707, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('W', 'Y')], 0.002553, places=6)
        self.assertAlmostEqual(self.obs_freq_mat[('Y', 'Y')], 0.004104, places=6)

    def test_obs_freq_mat_sum(self):
        counts = self.obs_freq_mat.sum()
        self.assertEqual(len(counts), 20)
        self.assertAlmostEqual(counts['A'], 0.0851, places=3)
        self.assertAlmostEqual(counts['C'], 0.0153, places=3)
        self.assertAlmostEqual(counts['D'], 0.0565, places=3)
        self.assertAlmostEqual(counts['E'], 0.0544, places=3)
        self.assertAlmostEqual(counts['F'], 0.0449, places=3)
        self.assertAlmostEqual(counts['G'], 0.0701, places=3)
        self.assertAlmostEqual(counts['H'], 0.0239, places=3)
        self.assertAlmostEqual(counts['I'], 0.0580, places=3)
        self.assertAlmostEqual(counts['K'], 0.0572, places=3)
        self.assertAlmostEqual(counts['L'], 0.0936, places=3)
        self.assertAlmostEqual(counts['M'], 0.0238, places=3)
        self.assertAlmostEqual(counts['N'], 0.0463, places=3)
        self.assertAlmostEqual(counts['P'], 0.0409, places=3)
        self.assertAlmostEqual(counts['Q'], 0.0360, places=3)
        self.assertAlmostEqual(counts['R'], 0.0437, places=3)
        self.assertAlmostEqual(counts['S'], 0.0575, places=3)
        self.assertAlmostEqual(counts['T'], 0.0565, places=3)
        self.assertAlmostEqual(counts['V'], 0.0780, places=3)
        self.assertAlmostEqual(counts['W'], 0.0162, places=3)
        self.assertAlmostEqual(counts['Y'], 0.0420, places=3)

    def test_log_odds_matris(self):
        lo_mat_prot = SubsMat.make_log_odds_matrix(acc_rep_mat=self.acc_rep_mat,
                                                   round_digit=1)
        handle = StringIO()
        lo_mat_prot.print_mat(f=handle, format=" %d",
                              alphabet='AVILMCFWYHSTNQKRDEGP')
        text = handle.getvalue()
        self.assertEqual(text, """\
A 0
V 0 1
I 0 0 1
L 0 0 0 0
M 0 0 0 0 1
C 0 0 0 0 0 3
F 0 0 0 0 0 0 1
W 0 0 0 0 0 0 0 2
Y 0 0 0 0 0 0 0 0 1
H 0 0 0 0 0 0 0 0 0 2
S 0 0 0 0 0 0 0 0 0 0 0
T 0 0 0 0 0 0 0 0 0 0 0 0
N 0 0 0 0 0 0 0 0 0 0 0 0 0
Q 0 0 0 0 0 0 0 0 0 0 0 0 0 0
K 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
R 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1
D 0 0 -1 0 -1 -1 0 0 0 0 0 0 0 0 0 0 1
E 0 0 -1 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 1
G 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1
P 0 0 0 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1
   A   V   I   L   M   C   F   W   Y   H   S   T   N   Q   K   R   D   E   G   P
""")
        handle = StringIO()
        lo_mat_prot.print_full_mat(f=handle, format=" %d",
                                   alphabet='AVILMCFWYHSTNQKRDEGP')
        text = handle.getvalue()
        self.assertEqual(text, """\
   A   V   I   L   M   C   F   W   Y   H   S   T   N   Q   K   R   D   E   G   P
A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
V 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
I 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 -1 0 0
L 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1
M 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 -1
C 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 -1 0 0 0
F 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 -1 0 0
W 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0
Y 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0
H 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0
S 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
T 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
N 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
Q 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
K 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
R 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0
D 0 0 -1 0 -1 -1 0 0 0 0 0 0 0 0 0 0 1 0 0 0
E 0 0 -1 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 1 0 0
G 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0
P 0 0 0 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1
""")
        relative_entropy = lo_mat_prot.calculate_relative_entropy(self.obs_freq_mat)
        self.assertAlmostEqual(relative_entropy, 0.162, places=3)

    def test_matrices(self):
        self.assertEqual(len(MatrixInfo.available_matrices), 40)
        mat = SubsMat.SeqMat(MatrixInfo.benner6)
        self.assertEqual(len(mat), 210)
        self.checkMatrix(mat, """\
A   2
C  -1  12
D   0  -3   5
E   0  -4   4   5
F  -3   0  -5  -6   8
G   0  -1   0   0  -5   5
H  -2  -1   0   0   0  -2   6
I   0  -3  -4  -4   0  -3  -3   4
K  -1  -2   0   0  -6  -1   0  -3   5
L  -1  -3  -5  -5   2  -4  -2   2  -4   4
M   0  -3  -4  -4   0  -3  -3   4  -2  -2   4
N   0  -1   2   1  -3   0   1  -2   1  -3  -2   3
P   1  -2  -2  -2  -3  -1   0  -2  -2   0  -1  -1   6
Q  -1  -3   0   2  -4  -1   3  -3   2  -2  -3   0   0   5
R  -1   0  -1   0  -4   0   1  -3   4  -3  -3   0  -1   2   5
S   1   0   0  -1  -1   0   0  -1  -1  -1  -1   1   1  -1   0   2
T   1  -1  -1  -1  -2   0  -1   0  -1   0   0   0   0  -1  -1   1   2
V   0  -3  -3  -3   0  -2  -3   3  -3   1   3  -2  -1  -3  -3   0   0   4
W  -4   1  -6  -5  -1  -1  -2  -5  -1  -3  -4  -4  -4  -2   2  -2  -2  -4  14
Y  -4   2  -2  -4   5  -4   4  -3  -4  -1  -3   0  -3  -1  -2  -1  -3  -3   0   9
   A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
""")
        mat = SubsMat.SeqMat(MatrixInfo.benner22)
        self.assertEqual(len(mat), 210)
        self.checkMatrix(mat, """\
A   2
C  -1  12
D   0  -3   4
E   0  -4   3   4
F  -3   0  -5  -5   7
G   0  -1   0   0  -5   6
H  -1  -1   0   0   0  -2   6
I   0  -2  -4  -3   0  -3  -3   4
K  -1  -3   0   1  -5  -1   0  -3   4
L  -1  -2  -4  -4   2  -4  -2   2  -3   4
M   0  -2  -3  -3   0  -3  -2   3  -2   3   4
N   0  -1   2   1  -3   0   1  -2   1  -3  -2   3
P   0  -3  -1  -1  -3  -1   0  -2  -1  -1  -2  -1   7
Q   0  -3   0   1  -3  -1   2  -2   2  -2  -1   0   0   4
R  -1  -1  -1   0  -4   0   1  -3   3  -2  -2   0  -1   2   5
S   1   0   0   0  -2   0   0  -1   0  -2  -1   1   1   0   0   2
T   1  -1   0   0  -2   0  -1   0   0  -1   0   0   0   0   0   1   2
V   0  -1  -3  -2   0  -2  -3   3  -2   2   2  -2  -1  -2  -2   0   0   3
W  -5   0  -6  -6   0  -4  -2  -4  -3  -1  -2  -5  -5  -3  -1  -3  -4  -4  15
Y  -3   0  -3  -4   5  -4   3  -2  -3   0  -1  -1  -3  -1  -2  -1  -3  -2   1   9
   A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
""")
        mat = SubsMat.SeqMat(MatrixInfo.benner74)
        self.assertEqual(len(mat), 210)
        self.checkMatrix(mat, """\
A   2
C   0  11
D   0  -3   4
E   0  -3   2   3
F  -2   0  -4  -4   7
G   0  -2   0   0  -5   6
H  -1  -1   0   0   0  -1   6
I   0  -1  -3  -2   0  -4  -2   4
K   0  -2   0   1  -3  -1   0  -2   3
L  -1  -1  -4  -3   2  -4  -1   2  -2   4
M   0  -1  -3  -2   1  -3  -1   2  -1   2   4
N   0  -1   2   1  -3   0   1  -2   0  -3  -2   3
P   0  -3  -1   0  -3  -1  -1  -2   0  -2  -2  -1   7
Q   0  -2   0   1  -2  -1   1  -2   1  -1  -1   0   0   3
R   0  -2   0   0  -3  -1   1  -2   2  -2  -1   0   0   1   4
S   1   0   0   0  -2   0   0  -1   0  -2  -1   0   0   0   0   2
T   0   0   0   0  -2  -1   0   0   0  -1   0   0   0   0   0   1   2
V   0   0  -2  -2   0  -3  -2   3  -1   1   1  -2  -1  -1  -2  -1   0   3
W  -4   0  -5  -4   3  -4  -1  -2  -3   0  -1  -4  -5  -2  -1  -3  -3  -2  14
Y  -2   0  -2  -3   5  -4   2  -1  -2   0   0  -1  -3  -1  -2  -1  -2  -1   3   8
   A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
""")
        mat = SubsMat.SeqMat(MatrixInfo.blosum100)
        self.assertEqual(len(mat), 276)
        self.checkMatrix(mat, """\
A   5
B  -3   4
C  -1  -5   9
D  -3   4  -5   7
E  -2   0  -6   1   6
F  -4  -5  -3  -5  -5   7
G  -1  -2  -5  -3  -4  -5   6
H  -3  -1  -5  -2  -1  -2  -4   9
I  -3  -5  -2  -6  -5  -1  -6  -5   5
K  -2  -1  -5  -2   0  -4  -3  -2  -4   6
L  -3  -5  -3  -6  -5   0  -5  -4   1  -4   5
M  -2  -4  -3  -5  -4  -1  -5  -3   1  -2   2   8
N  -2   4  -4   1  -1  -5  -2   0  -5  -1  -5  -4   7
P  -1  -3  -5  -3  -3  -5  -4  -3  -4  -2  -4  -4  -4   8
Q  -1  -1  -5  -2   1  -4  -3   0  -4   1  -3  -1  -1  -2   7
R  -2  -2  -5  -3  -2  -4  -4  -1  -4   2  -4  -2  -1  -3   0   7
S   1  -1  -2  -1  -1  -3  -1  -2  -4  -1  -4  -3   0  -2  -1  -2   6
T  -1  -2  -2  -2  -2  -3  -3  -3  -2  -2  -3  -2  -1  -3  -2  -2   1   6
V  -1  -5  -2  -5  -3  -2  -5  -5   2  -4   0   0  -4  -4  -3  -4  -3  -1   5
W  -4  -6  -5  -7  -5   0  -5  -3  -4  -5  -4  -3  -6  -6  -3  -4  -4  -5  -4  11
X  -1  -2  -3  -3  -2  -3  -3  -2  -2  -2  -2  -2  -2  -3  -2  -2  -1  -1  -2  -4  -2
Y  -4  -4  -4  -5  -4   3  -6   1  -3  -4  -3  -3  -3  -5  -3  -3  -3  -3  -3   1  -3   8
Z  -2   1  -6   0   5  -5  -4  -1  -4   0  -4  -3  -1  -3   3  -1  -1  -2  -3  -4  -2  -4   4
   A   B   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   X   Y   Z
""")
        mat = SubsMat.SeqMat(MatrixInfo.blosum30)
        self.assertEqual(len(mat), 276)
        self.checkMatrix(mat, """\
A   4
B   0   5
C  -3  -2  17
D   0   5  -3   9
E   0   0   1   1   6
F  -2  -3  -3  -5  -4  10
G   0   0  -4  -1  -2  -3   8
H  -2  -2  -5  -2   0  -3  -3  14
I   0  -2  -2  -4  -3   0  -1  -2   6
K   0   0  -3   0   2  -1  -1  -2  -2   4
L  -1  -1   0  -1  -1   2  -2  -1   2  -2   4
M   1  -2  -2  -3  -1  -2  -2   2   1   2   2   6
N   0   4  -1   1  -1  -1   0  -1   0   0  -2   0   8
P  -1  -2  -3  -1   1  -4  -1   1  -3   1  -3  -4  -3  11
Q   1  -1  -2  -1   2  -3  -2   0  -2   0  -2  -1  -1   0   8
R  -1  -2  -2  -1  -1  -1  -2  -1  -3   1  -2   0  -2  -1   3   8
S   1   0  -2   0   0  -1   0  -1  -1   0  -2  -2   0  -1  -1  -1   4
T   1   0  -2  -1  -2  -2  -2  -2   0  -1   0   0   1   0   0  -3   2   5
V   1  -2  -2  -2  -3   1  -3  -3   4  -2   1   0  -2  -4  -3  -1  -1   1   5
W  -5  -5  -2  -4  -1   1   1  -5  -3  -2  -2  -3  -7  -3  -1   0  -3  -5  -3  20
X   0  -1  -2  -1  -1  -1  -1  -1   0   0   0   0   0  -1   0  -1   0   0   0  -2  -1
Y  -4  -3  -6  -1  -2   3  -3   0  -1  -1   3  -1  -4  -2  -1   0  -2  -1   1   5  -1   9
Z   0   0   0   0   5  -4  -2   0  -3   1  -1  -1  -1   0   4   0  -1  -1  -3  -1   0  -2   4
   A   B   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   X   Y   Z
""")
        mat = SubsMat.SeqMat(MatrixInfo.blosum35)
        self.assertEqual(len(mat), 276)
        self.checkMatrix(mat, """\
A   5
B  -1   5
C  -2  -2  15
D  -1   5  -3   8
E  -1   0  -1   2   6
F  -2  -2  -4  -3  -3   8
G   0   0  -3  -2  -2  -3   7
H  -2   0  -4   0  -1  -3  -2  12
I  -1  -2  -4  -3  -3   1  -3  -3   5
K   0   0  -2  -1   1  -1  -1  -2  -2   5
L  -2  -2  -2  -2  -1   2  -3  -2   2  -2   5
M   0  -2  -4  -3  -2   0  -1   1   1   0   3   6
N  -1   4  -1   1  -1  -1   1   1  -1   0  -2  -1   7
P  -2  -1  -4  -1   0  -4  -2  -1  -1   0  -3  -3  -2  10
Q   0   0  -3  -1   2  -4  -2  -1  -2   0  -2  -1   1   0   7
R  -1  -1  -3  -1  -1  -1  -2  -1  -3   2  -2   0  -1  -2   2   8
S   1   0  -3  -1   0  -1   1  -1  -2   0  -2  -1   0  -2   0  -1   4
T   0  -1  -1  -1  -1  -1  -2  -2  -1   0   0   0   0   0   0  -2   2   5
V   0  -2  -2  -2  -2   1  -3  -4   4  -2   2   1  -2  -3  -3  -1  -1   1   5
W  -2  -3  -5  -3  -1   1  -1  -4  -1   0   0   1  -2  -4  -1   0  -2  -2  -2  16
X   0  -1  -2  -1  -1  -1  -1  -1   0   0   0   0   0  -1  -1  -1   0   0   0  -1  -1
Y  -1  -2  -5  -2  -1   3  -2   0   0  -1   0   0  -2  -3   0   0  -1  -2   0   3  -1   8
Z  -1   0  -2   1   5  -3  -2  -1  -3   1  -2  -2   0   0   4   0   0  -1  -2  -1   0  -1   4
   A   B   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   X   Y   Z
""")
        mat = SubsMat.SeqMat(MatrixInfo.blosum40)
        self.assertEqual(len(mat), 276)
        self.checkMatrix(mat, """\
A   5
B  -1   5
C  -2  -2  16
D  -1   6  -2   9
E  -1   1  -2   2   7
F  -3  -3  -2  -4  -3   9
G   1  -1  -3  -2  -3  -3   8
H  -2   0  -4   0   0  -2  -2  13
I  -1  -3  -4  -4  -4   1  -4  -3   6
K  -1   0  -3   0   1  -3  -2  -1  -3   6
L  -2  -3  -2  -3  -2   2  -4  -2   2  -2   6
M  -1  -3  -3  -3  -2   0  -2   1   1  -1   3   7
N  -1   4  -2   2  -1  -3   0   1  -2   0  -3  -2   8
P  -2  -2  -5  -2   0  -4  -1  -2  -2  -1  -4  -2  -2  11
Q   0   0  -4  -1   2  -4  -2   0  -3   1  -2  -1   1  -2   8
R  -2  -1  -3  -1  -1  -2  -3   0  -3   3  -2  -1   0  -3   2   9
S   1   0  -1   0   0  -2   0  -1  -2   0  -3  -2   1  -1   1  -1   5
T   0   0  -1  -1  -1  -1  -2  -2  -1   0  -1  -1   0   0  -1  -2   2   6
V   0  -3  -2  -3  -3   0  -4  -4   4  -2   2   1  -3  -3  -3  -2  -1   1   5
W  -3  -4  -6  -5  -2   1  -2  -5  -3  -2  -1  -2  -4  -4  -1  -2  -5  -4  -3  19
X   0  -1  -2  -1  -1  -1  -1  -1  -1  -1  -1   0  -1  -2  -1  -1   0   0  -1  -2  -1
Y  -2  -3  -4  -3  -2   4  -3   2   0  -1   0   1  -2  -3  -1  -1  -2  -1  -1   3  -1   9
Z  -1   2  -3   1   5  -4  -2   0  -4   1  -2  -2   0  -1   4   0   0  -1  -3  -2  -1  -2   5
   A   B   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   X   Y   Z
""")
        mat = SubsMat.SeqMat(MatrixInfo.blosum45)
        self.assertEqual(len(mat), 276)
        self.checkMatrix(mat, """\
A   5
B  -1   4
C  -1  -2  12
D  -2   5  -3   7
E  -1   1  -3   2   6
F  -2  -3  -2  -4  -3   8
G   0  -1  -3  -1  -2  -3   7
H  -2   0  -3   0   0  -2  -2  10
I  -1  -3  -3  -4  -3   0  -4  -3   5
K  -1   0  -3   0   1  -3  -2  -1  -3   5
L  -1  -3  -2  -3  -2   1  -3  -2   2  -3   5
M  -1  -2  -2  -3  -2   0  -2   0   2  -1   2   6
N  -1   4  -2   2   0  -2   0   1  -2   0  -3  -2   6
P  -1  -2  -4  -1   0  -3  -2  -2  -2  -1  -3  -2  -2   9
Q  -1   0  -3   0   2  -4  -2   1  -2   1  -2   0   0  -1   6
R  -2  -1  -3  -1   0  -2  -2   0  -3   3  -2  -1   0  -2   1   7
S   1   0  -1   0   0  -2   0  -1  -2  -1  -3  -2   1  -1   0  -1   4
T   0   0  -1  -1  -1  -1  -2  -2  -1  -1  -1  -1   0  -1  -1  -1   2   5
V   0  -3  -1  -3  -3   0  -3  -3   3  -2   1   1  -3  -3  -3  -2  -1   0   5
W  -2  -4  -5  -4  -3   1  -2  -3  -2  -2  -2  -2  -4  -3  -2  -2  -4  -3  -3  15
X   0  -1  -2  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1   0   0  -1  -2  -1
Y  -2  -2  -3  -2  -2   3  -3   2   0  -1   0   0  -2  -3  -1  -1  -2  -1  -1   3  -1   8
Z  -1   2  -3   1   4  -3  -2   0  -3   1  -2  -1   0  -1   4   0   0  -1  -3  -2  -1  -2   4
   A   B   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   X   Y   Z
""")
        mat = SubsMat.SeqMat(MatrixInfo.blosum50)
        self.assertEqual(len(mat), 276)
        self.checkMatrix(mat, """\
A   5
B  -2   5
C  -1  -3  13
D  -2   5  -4   8
E  -1   1  -3   2   6
F  -3  -4  -2  -5  -3   8
G   0  -1  -3  -1  -3  -4   8
H  -2   0  -3  -1   0  -1  -2  10
I  -1  -4  -2  -4  -4   0  -4  -4   5
K  -1   0  -3  -1   1  -4  -2   0  -3   6
L  -2  -4  -2  -4  -3   1  -4  -3   2  -3   5
M  -1  -3  -2  -4  -2   0  -3  -1   2  -2   3   7
N  -1   4  -2   2   0  -4   0   1  -3   0  -4  -2   7
P  -1  -2  -4  -1  -1  -4  -2  -2  -3  -1  -4  -3  -2  10
Q  -1   0  -3   0   2  -4  -2   1  -3   2  -2   0   0  -1   7
R  -2  -1  -4  -2   0  -3  -3   0  -4   3  -3  -2  -1  -3   1   7
S   1   0  -1   0  -1  -3   0  -1  -3   0  -3  -2   1  -1   0  -1   5
T   0   0  -1  -1  -1  -2  -2  -2  -1  -1  -1  -1   0  -1  -1  -1   2   5
V   0  -4  -1  -4  -3  -1  -4  -4   4  -3   1   1  -3  -3  -3  -3  -2   0   5
W  -3  -5  -5  -5  -3   1  -3  -3  -3  -3  -2  -1  -4  -4  -1  -3  -4  -3  -3  15
X  -1  -1  -2  -1  -1  -2  -2  -1  -1  -1  -1  -1  -1  -2  -1  -1  -1   0  -1  -3  -1
Y  -2  -3  -3  -3  -2   4  -3   2  -1  -2  -1   0  -2  -3  -1  -1  -2  -2  -1   2  -1   8
Z  -1   2  -3   1   5  -4  -2   0  -3   1  -3  -1   0  -1   4   0   0  -1  -3  -2  -1  -2   5
   A   B   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   X   Y   Z
""")
        mat = SubsMat.SeqMat(MatrixInfo.blosum55)
        self.assertEqual(len(mat), 276)
        self.checkMatrix(mat, """\
A   5
B  -2   5
C  -1  -3  13
D  -2   5  -4   8
E  -1   1  -3   2   6
F  -3  -4  -2  -5  -3   8
G   0  -1  -3  -1  -3  -4   8
H  -2   0  -3  -1   0  -1  -2  10
I  -1  -4  -2  -4  -4   0  -4  -4   5
K  -1   0  -3  -1   1  -4  -2   0  -3   6
L  -2  -4  -2  -4  -3   1  -4  -3   2  -3   5
M  -1  -3  -2  -4  -2   0  -3  -1   2  -2   3   7
N  -1   4  -2   2   0  -4   0   1  -3   0  -4  -2   7
P  -1  -2  -4  -1  -1  -4  -2  -2  -3  -1  -4  -3  -2  10
Q  -1   0  -3   0   2  -4  -2   1  -3   2  -2   0   0  -1   7
R  -2  -1  -4  -2   0  -3  -3   0  -4   3  -3  -2  -1  -3   1   7
S   1   0  -1   0  -1  -3   0  -1  -3   0  -3  -2   1  -1   0  -1   5
T   0   0  -1  -1  -1  -2  -2  -2  -1  -1  -1  -1   0  -1  -1  -1   2   5
V   0  -4  -1  -4  -3  -1  -4  -4   4  -3   1   1  -3  -3  -3  -3  -2   0   5
W  -3  -5  -5  -5  -3   1  -3  -3  -3  -3  -2  -1  -4  -4  -1  -3  -4  -3  -3  15
X  -1  -1  -2  -1  -1  -2  -2  -1  -1  -1  -1  -1  -1  -2  -1  -1  -1   0  -1  -3  -1
Y  -2  -3  -3  -3  -2   4  -3   2  -1  -2  -1   0  -2  -3  -1  -1  -2  -2  -1   2  -1   8
Z  -1   2  -3   1   5  -4  -2   0  -3   1  -3  -1   0  -1   4   0   0  -1  -3  -2  -1  -2   5
   A   B   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   X   Y   Z
""")
        mat = SubsMat.SeqMat(MatrixInfo.blosum60)
        self.assertEqual(len(mat), 276)
        self.checkMatrix(mat, """\
A   4
B  -2   4
C   0  -3   9
D  -2   4  -3   6
E  -1   1  -3   2   5
F  -2  -3  -2  -3  -3   6
G   0  -1  -2  -1  -2  -3   6
H  -2   0  -3  -1   0  -1  -2   7
I  -1  -3  -1  -3  -3   0  -3  -3   4
K  -1   0  -3  -1   1  -3  -1  -1  -3   4
L  -1  -3  -1  -3  -3   0  -4  -3   2  -2   4
M  -1  -3  -1  -3  -2   0  -2  -1   1  -1   2   5
N  -1   3  -2   1   0  -3   0   1  -3   0  -3  -2   6
P  -1  -2  -3  -1  -1  -4  -2  -2  -3  -1  -3  -2  -2   7
Q  -1   0  -3   0   2  -3  -2   1  -3   1  -2   0   0  -1   5
R  -1  -1  -3  -1   0  -3  -2   0  -3   2  -2  -1   0  -2   1   5
S   1   0  -1   0   0  -2   0  -1  -2   0  -2  -1   1  -1   0  -1   4
T   0   0  -1  -1  -1  -2  -2  -2  -1  -1  -1  -1   0  -1  -1  -1   1   4
V   0  -3  -1  -3  -2  -1  -3  -3   3  -2   1   1  -3  -2  -2  -2  -2   0   4
W  -3  -4  -2  -4  -3   1  -2  -2  -2  -3  -2  -1  -4  -4  -2  -3  -3  -2  -3  10
X   0  -1  -2  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -2  -1  -1   0   0  -1  -2  -1
Y  -2  -2  -2  -3  -2   3  -3   2  -1  -2  -1  -1  -2  -3  -1  -2  -2  -2  -1   2  -1   6
Z  -1   1  -3   1   4  -3  -2   0  -3   1  -2  -1   0  -1   3   0   0  -1  -2  -2  -1  -2   3
   A   B   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   X   Y   Z
""")
        mat = SubsMat.SeqMat(MatrixInfo.blosum62)
        self.assertEqual(len(mat), 276)
        self.checkMatrix(mat, """\
A   4
B  -2   4
C   0  -3   9
D  -2   4  -3   6
E  -1   1  -4   2   5
F  -2  -3  -2  -3  -3   6
G   0  -1  -3  -1  -2  -3   6
H  -2   0  -3  -1   0  -1  -2   8
I  -1  -3  -1  -3  -3   0  -4  -3   4
K  -1   0  -3  -1   1  -3  -2  -1  -3   5
L  -1  -4  -1  -4  -3   0  -4  -3   2  -2   4
M  -1  -3  -1  -3  -2   0  -3  -2   1  -1   2   5
N  -2   3  -3   1   0  -3   0   1  -3   0  -3  -2   6
P  -1  -2  -3  -1  -1  -4  -2  -2  -3  -1  -3  -2  -2   7
Q  -1   0  -3   0   2  -3  -2   0  -3   1  -2   0   0  -1   5
R  -1  -1  -3  -2   0  -3  -2   0  -3   2  -2  -1   0  -2   1   5
S   1   0  -1   0   0  -2   0  -1  -2   0  -2  -1   1  -1   0  -1   4
T   0  -1  -1  -1  -1  -2  -2  -2  -1  -1  -1  -1   0  -1  -1  -1   1   5
V   0  -3  -1  -3  -2  -1  -3  -3   3  -2   1   1  -3  -2  -2  -3  -2   0   4
W  -3  -4  -2  -4  -3   1  -2  -2  -3  -3  -2  -1  -4  -4  -2  -3  -3  -2  -3  11
X   0  -1  -2  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -2  -1  -1   0   0  -1  -2  -1
Y  -2  -3  -2  -3  -2   3  -3   2  -1  -2  -1  -1  -2  -3  -1  -2  -2  -2  -1   2  -1   7
Z  -1   1  -3   1   4  -3  -2   0  -3   1  -3  -1   0  -1   3   0   0  -1  -2  -3  -1  -2   4
   A   B   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   X   Y   Z
""")
        mat = SubsMat.SeqMat(MatrixInfo.blosum65)
        self.assertEqual(len(mat), 276)
        self.checkMatrix(mat, """\
A   4
B  -2   4
C   0  -3   9
D  -2   4  -4   6
E  -1   1  -4   2   5
F  -2  -3  -2  -4  -3   6
G   0  -1  -3  -1  -2  -3   6
H  -2   0  -3  -1   0  -1  -2   8
I  -1  -3  -1  -3  -3   0  -4  -3   4
K  -1   0  -3  -1   1  -3  -2  -1  -3   5
L  -2  -4  -1  -4  -3   0  -4  -3   2  -3   4
M  -1  -3  -2  -3  -2   0  -3  -2   1  -2   2   6
N  -2   3  -3   1   0  -3  -1   1  -3   0  -4  -2   6
P  -1  -2  -3  -2  -1  -4  -2  -2  -3  -1  -3  -3  -2   8
Q  -1   0  -3   0   2  -3  -2   1  -3   1  -2   0   0  -1   6
R  -1  -1  -4  -2   0  -3  -2   0  -3   2  -2  -2   0  -2   1   6
S   1   0  -1   0   0  -2   0  -1  -2   0  -3  -2   1  -1   0  -1   4
T   0  -1  -1  -1  -1  -2  -2  -2  -1  -1  -1  -1   0  -1  -1  -1   1   5
V   0  -3  -1  -3  -3  -1  -3  -3   3  -2   1   1  -3  -2  -2  -3  -2   0   4
W  -3  -4  -2  -5  -3   1  -3  -2  -2  -3  -2  -2  -4  -4  -2  -3  -3  -3  -3  10
X  -1  -1  -2  -1  -1  -2  -2  -1  -1  -1  -1  -1  -1  -2  -1  -1  -1  -1  -1  -2  -1
Y  -2  -3  -2  -3  -2   3  -3   2  -1  -2  -1  -1  -2  -3  -2  -2  -2  -2  -1   2  -1   7
Z  -1   1  -4   1   4  -3  -2   0  -3   1  -3  -2   0  -1   3   0   0  -1  -2  -3  -1  -2   4
   A   B   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   X   Y   Z
""")
        mat = SubsMat.SeqMat(MatrixInfo.blosum70)
        self.assertEqual(len(mat), 276)
        self.checkMatrix(mat, """\
A   4
B  -2   4
C  -1  -4   9
D  -2   4  -4   6
E  -1   1  -4   1   5
F  -2  -4  -2  -4  -4   6
G   0  -1  -3  -2  -2  -4   6
H  -2  -1  -4  -1   0  -1  -2   8
I  -2  -4  -1  -4  -4   0  -4  -4   4
K  -1  -1  -4  -1   1  -3  -2  -1  -3   5
L  -2  -4  -2  -4  -3   0  -4  -3   2  -3   4
M  -1  -3  -2  -3  -2   0  -3  -2   1  -2   2   6
N  -2   3  -3   1   0  -3  -1   0  -4   0  -4  -2   6
P  -1  -2  -3  -2  -1  -4  -3  -2  -3  -1  -3  -3  -2   8
Q  -1   0  -3  -1   2  -3  -2   1  -3   1  -2   0   0  -2   6
R  -2  -1  -4  -2   0  -3  -3   0  -3   2  -3  -2  -1  -2   1   6
S   1   0  -1   0   0  -3  -1  -1  -3   0  -3  -2   0  -1   0  -1   4
T   0  -1  -1  -1  -1  -2  -2  -2  -1  -1  -2  -1   0  -1  -1  -1   1   5
V   0  -3  -1  -4  -3  -1  -4  -3   3  -3   1   1  -3  -3  -2  -3  -2   0   4
W  -3  -4  -3  -5  -4   1  -3  -2  -3  -3  -2  -2  -4  -4  -2  -3  -3  -3  -3  11
X  -1  -1  -2  -2  -1  -2  -2  -1  -1  -1  -1  -1  -1  -2  -1  -1  -1  -1  -1  -3  -1
Y  -2  -3  -3  -4  -3   3  -4   2  -1  -2  -1  -1  -2  -3  -2  -2  -2  -2  -2   2  -2   7
Z  -1   0  -4   1   4  -4  -2   0  -3   1  -3  -2   0  -1   3   0   0  -1  -3  -3  -1  -2   4
   A   B   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   X   Y   Z
""")
        mat = SubsMat.SeqMat(MatrixInfo.blosum75)
        self.assertEqual(len(mat), 276)
        self.checkMatrix(mat, """\
A   4
B  -2   4
C  -1  -4   9
D  -2   4  -4   6
E  -1   1  -5   1   5
F  -3  -4  -2  -4  -4   6
G   0  -1  -3  -2  -3  -4   6
H  -2  -1  -4  -1   0  -2  -2   8
I  -2  -4  -1  -4  -4   0  -5  -4   4
K  -1  -1  -4  -1   1  -4  -2  -1  -3   5
L  -2  -4  -2  -4  -4   0  -4  -3   1  -3   4
M  -1  -3  -2  -4  -2   0  -3  -2   1  -2   2   6
N  -2   3  -3   1  -1  -4  -1   0  -4   0  -4  -3   6
P  -1  -2  -4  -2  -1  -4  -3  -2  -3  -1  -3  -3  -3   8
Q  -1   0  -3  -1   2  -4  -2   1  -3   1  -3   0   0  -2   6
R  -2  -1  -4  -2   0  -3  -3   0  -3   2  -3  -2  -1  -2   1   6
S   1   0  -1  -1   0  -3  -1  -1  -3   0  -3  -2   0  -1   0  -1   5
T   0  -1  -1  -1  -1  -2  -2  -2  -1  -1  -2  -1   0  -1  -1  -1   1   5
V   0  -4  -1  -4  -3  -1  -4  -4   3  -3   1   1  -3  -3  -2  -3  -2   0   4
W  -3  -5  -3  -5  -4   1  -3  -2  -3  -4  -2  -2  -4  -5  -2  -3  -3  -3  -3  11
X  -1  -2  -2  -2  -1  -2  -2  -1  -2  -1  -1  -1  -1  -2  -1  -1  -1  -1  -1  -3  -1
Y  -2  -3  -3  -4  -3   3  -4   2  -2  -2  -1  -2  -3  -4  -2  -2  -2  -2  -2   2  -2   7
Z  -1   0  -4   1   4  -4  -2   0  -4   1  -3  -2   0  -2   3   0   0  -1  -3  -3  -1  -3   4
   A   B   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   X   Y   Z
""")
        mat = SubsMat.SeqMat(MatrixInfo.blosum80)
        self.assertEqual(len(mat), 276)
        self.checkMatrix(mat, """\
A   5
B  -2   4
C  -1  -4   9
D  -2   4  -4   6
E  -1   1  -5   1   6
F  -3  -4  -3  -4  -4   6
G   0  -1  -4  -2  -3  -4   6
H  -2  -1  -4  -2   0  -2  -3   8
I  -2  -4  -2  -4  -4  -1  -5  -4   5
K  -1  -1  -4  -1   1  -4  -2  -1  -3   5
L  -2  -4  -2  -5  -4   0  -4  -3   1  -3   4
M  -1  -3  -2  -4  -2   0  -4  -2   1  -2   2   6
N  -2   4  -3   1  -1  -4  -1   0  -4   0  -4  -3   6
P  -1  -2  -4  -2  -2  -4  -3  -3  -4  -1  -3  -3  -3   8
Q  -1   0  -4  -1   2  -4  -2   1  -3   1  -3   0   0  -2   6
R  -2  -2  -4  -2  -1  -4  -3   0  -3   2  -3  -2  -1  -2   1   6
S   1   0  -2  -1   0  -3  -1  -1  -3  -1  -3  -2   0  -1   0  -1   5
T   0  -1  -1  -1  -1  -2  -2  -2  -1  -1  -2  -1   0  -2  -1  -1   1   5
V   0  -4  -1  -4  -3  -1  -4  -4   3  -3   1   1  -4  -3  -3  -3  -2   0   4
W  -3  -5  -3  -6  -4   0  -4  -3  -3  -4  -2  -2  -4  -5  -3  -4  -4  -4  -3  11
X  -1  -2  -3  -2  -1  -2  -2  -2  -2  -1  -2  -1  -1  -2  -1  -1  -1  -1  -1  -3  -1
Y  -2  -3  -3  -4  -3   3  -4   2  -2  -3  -2  -2  -3  -4  -2  -3  -2  -2  -2   2  -2   7
Z  -1   0  -4   1   4  -4  -3   0  -4   1  -3  -2   0  -2   3   0   0  -1  -3  -4  -1  -3   4
   A   B   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   X   Y   Z
""")
        mat = SubsMat.SeqMat(MatrixInfo.blosum85)
        self.assertEqual(len(mat), 276)
        self.checkMatrix(mat, """\
A   5
B  -2   4
C  -1  -4   9
D  -2   4  -5   7
E  -1   0  -5   1   6
F  -3  -4  -3  -4  -4   7
G   0  -1  -4  -2  -3  -4   6
H  -2  -1  -5  -2  -1  -2  -3   8
I  -2  -5  -2  -5  -4  -1  -5  -4   5
K  -1  -1  -4  -1   0  -4  -2  -1  -3   6
L  -2  -5  -2  -5  -4   0  -5  -3   1  -3   4
M  -2  -4  -2  -4  -3  -1  -4  -3   1  -2   2   7
N  -2   4  -4   1  -1  -4  -1   0  -4   0  -4  -3   7
P  -1  -3  -4  -2  -2  -4  -3  -3  -4  -2  -4  -3  -3   8
Q  -1  -1  -4  -1   2  -4  -3   1  -4   1  -3   0   0  -2   6
R  -2  -2  -4  -2  -1  -4  -3   0  -4   2  -3  -2  -1  -2   1   6
S   1   0  -2  -1  -1  -3  -1  -1  -3  -1  -3  -2   0  -1  -1  -1   5
T   0  -1  -2  -2  -1  -3  -2  -2  -1  -1  -2  -1   0  -2  -1  -2   1   5
V  -1  -4  -1  -4  -3  -1  -4  -4   3  -3   0   0  -4  -3  -3  -3  -2   0   5
W  -3  -5  -4  -6  -4   0  -4  -3  -3  -5  -3  -2  -5  -5  -3  -4  -4  -4  -3  11
X  -1  -2  -3  -2  -1  -2  -2  -2  -2  -1  -2  -1  -2  -2  -1  -2  -1  -1  -1  -3  -2
Y  -3  -4  -3  -4  -4   3  -5   2  -2  -3  -2  -2  -3  -4  -2  -3  -2  -2  -2   2  -2   7
Z  -1   0  -5   1   4  -4  -3   0  -4   1  -4  -2  -1  -2   4   0  -1  -1  -3  -4  -1  -3   4
   A   B   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   X   Y   Z
""")
        mat = SubsMat.SeqMat(MatrixInfo.blosum90)
        self.assertEqual(len(mat), 276)
        self.checkMatrix(mat, """\
A   5
B  -2   4
C  -1  -4   9
D  -3   4  -5   7
E  -1   0  -6   1   6
F  -3  -4  -3  -5  -5   7
G   0  -2  -4  -2  -3  -5   6
H  -2  -1  -5  -2  -1  -2  -3   8
I  -2  -5  -2  -5  -4  -1  -5  -4   5
K  -1  -1  -4  -1   0  -4  -2  -1  -4   6
L  -2  -5  -2  -5  -4   0  -5  -4   1  -3   5
M  -2  -4  -2  -4  -3  -1  -4  -3   1  -2   2   7
N  -2   4  -4   1  -1  -4  -1   0  -4   0  -4  -3   7
P  -1  -3  -4  -3  -2  -4  -3  -3  -4  -2  -4  -3  -3   8
Q  -1  -1  -4  -1   2  -4  -3   1  -4   1  -3   0   0  -2   7
R  -2  -2  -5  -3  -1  -4  -3   0  -4   2  -3  -2  -1  -3   1   6
S   1   0  -2  -1  -1  -3  -1  -2  -3  -1  -3  -2   0  -2  -1  -1   5
T   0  -1  -2  -2  -1  -3  -3  -2  -1  -1  -2  -1   0  -2  -1  -2   1   6
V  -1  -4  -2  -5  -3  -2  -5  -4   3  -3   0   0  -4  -3  -3  -3  -2  -1   5
W  -4  -6  -4  -6  -5   0  -4  -3  -4  -5  -3  -2  -5  -5  -3  -4  -4  -4  -3  11
X  -1  -2  -3  -2  -2  -2  -2  -2  -2  -1  -2  -1  -2  -2  -1  -2  -1  -1  -2  -3  -2
Y  -3  -4  -4  -4  -4   3  -5   1  -2  -3  -2  -2  -3  -4  -3  -3  -3  -2  -3   2  -2   8
Z  -1   0  -5   0   4  -4  -3   0  -4   1  -4  -2  -1  -2   4   0  -1  -1  -3  -4  -1  -3   4
   A   B   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   X   Y   Z
""")
        mat = SubsMat.SeqMat(MatrixInfo.blosum95)
        self.assertEqual(len(mat), 276)
        self.checkMatrix(mat, """\
A   5
B  -3   4
C  -1  -4   9
D  -3   4  -5   7
E  -1   0  -6   1   6
F  -3  -5  -3  -5  -5   7
G  -1  -2  -5  -2  -3  -5   6
H  -3  -1  -5  -2  -1  -2  -3   9
I  -2  -5  -2  -5  -4  -1  -6  -4   5
K  -1  -1  -5  -2   0  -4  -3  -1  -4   6
L  -2  -5  -3  -5  -4   0  -5  -4   1  -3   5
M  -2  -4  -3  -5  -3  -1  -4  -3   1  -2   2   7
N  -2   4  -4   1  -1  -4  -1   0  -4   0  -5  -3   7
P  -1  -3  -5  -3  -2  -5  -4  -3  -4  -2  -4  -3  -3   8
Q  -1  -1  -4  -1   2  -4  -3   1  -4   1  -3  -1   0  -2   7
R  -2  -2  -5  -3  -1  -4  -4  -1  -4   2  -3  -2  -1  -3   0   7
S   1  -1  -2  -1  -1  -3  -1  -2  -3  -1  -3  -3   0  -2  -1  -2   5
T   0  -1  -2  -2  -2  -3  -3  -2  -2  -1  -2  -1  -1  -2  -1  -2   1   6
V  -1  -5  -2  -5  -3  -2  -5  -4   3  -3   0   0  -4  -4  -3  -4  -3  -1   5
W  -4  -6  -4  -6  -5   0  -5  -3  -4  -5  -3  -2  -5  -5  -3  -4  -4  -4  -3  11
X  -1  -2  -3  -2  -2  -2  -3  -2  -2  -1  -2  -2  -2  -3  -1  -2  -1  -1  -2  -4  -2
Y  -3  -4  -4  -5  -4   3  -5   1  -2  -3  -2  -3  -3  -5  -3  -3  -3  -3  -3   2  -2   8
Z  -1   0  -5   0   4  -4  -3   0  -4   0  -4  -2  -1  -2   4  -1  -1  -2  -3  -4  -1  -4   4
   A   B   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   X   Y   Z
""")
        mat = SubsMat.SeqMat(MatrixInfo.feng)
        self.assertEqual(len(mat), 210)
        self.checkMatrix(mat, """\
A   6
C   2   6
D   4   1   6
E   4   0   5   6
F   2   3   1   0   6
G   5   3   4   4   1   6
H   2   2   3   2   2   1   6
I   2   2   1   1   4   2   1   6
K   3   0   3   4   0   2   3   2   6
L   2   2   1   1   4   2   3   5   2   6
M   2   2   0   1   2   1   1   4   2   5   6
N   3   2   5   3   1   3   4   2   4   1   1   6
P   5   2   2   3   2   3   3   2   2   3   2   2   6
Q   3   1   4   4   1   2   4   1   4   2   2   3   3   6
R   2   2   2   2   1   3   4   2   5   2   2   2   3   3   6
S   5   4   3   3   3   5   3   2   3   2   1   5   4   3   3   6
T   5   2   2   3   1   2   2   3   4   2   3   4   4   3   3   5   6
V   5   2   3   4   4   4   1   5   3   5   4   2   3   2   2   2   3   6
W   2   3   0   1   3   3   1   2   1   4   3   0   2   1   2   2   1   3   6
Y   2   3   2   1   5   2   3   3   1   3   2   3   2   2   1   3   2   3   3   6
   A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
""")
        mat = SubsMat.SeqMat(MatrixInfo.fitch)
        self.assertEqual(len(mat), 171)
        self.checkMatrix(mat, """\
A   3
C   1   3
E   1   1   3
F   1   2   1   3
H   2   1   1   1   3
I   1   0   2   0   1   3
L   2   1   2   1   1   1   3
M   0   0   2   1   0   2   1   3
N   2   1   2   1   2   2   1   1   3
O   2   2   1   2   2   1   1   0   2   3
Q   1   1   1   1   2   2   1   1   1   2   3
R   1   2   2   1   2   2   1   2   1   1   2   3
S   1   2   2   2   1   1   2   1   2   2   2   2   3
T   0   2   1   1   0   1   1   1   0   1   2   2   2   3
U   1   1   1   2   2   1   1   2   1   1   2   2   2   2   3
V   2   2   1   2   1   1   2   2   1   1   1   1   1   1   2   3
W   1   1   2   2   1   1   1   2   2   1   0   1   2   0   2   2   3
Y   2   2   1   1   1   1   2   1   1   1   1   2   2   2   2   2   1   3
   A   C   E   F   H   I   L   M   N   O   Q   R   S   T   U   V   W   Y
""")
        mat = SubsMat.SeqMat(MatrixInfo.genetic)
        self.assertEqual(len(mat), 210)
        self.checkMatrix(mat, """\
A   4
C  -1   5
D   1  -1   4
E   1  -3   3   5
F  -2   1  -1  -2   4
G   1   1   1   1  -1   4
H  -2  -1   1   0  -1  -2   4
I  -1  -1  -2  -2   1  -2  -1   4
K  -1  -3   0   2  -2  -2   0   0   5
L  -2  -1  -2  -2   2  -2   0   1  -2   3
M  -2  -2  -2  -1   0  -2  -1   3   1   1   5
N  -1  -1   1   0  -1  -2   1   0   3  -2   0   4
P   0  -1  -2  -2  -1  -1   0  -1  -1   0  -1  -1   3
Q  -2  -3   0   2  -2  -2   3  -1   2   0  -1   0   1   5
R  -1   0  -2  -2  -1   0   3  -1   0   0   0  -1   0   0   2
S   0   1  -2  -2   0   0  -1   0  -1  -1  -1   0   0  -2   0   2
T   0  -1  -2  -2  -2  -2  -1   0   1  -1   0   0   1  -1   0   1   4
V   1  -2   1   1   1   1  -2   1  -2   1   1  -2  -2  -2  -2  -2  -2   4
W  -2   4  -2  -3   0   1  -2  -2  -3   0  -2  -3  -1  -2   1   0  -2  -2   7
Y  -2   2   2   0   2  -1   2  -1   0  -1  -2   2  -2   0  -1   0  -2  -2   0   6
   A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
""")
        mat = SubsMat.SeqMat(MatrixInfo.gonnet)
        self.assertEqual(len(mat), 210)
        self.checkMatrix(mat, """\
A   2
C   0  11
D   0  -3   4
E   0  -3   2   3
F  -2   0  -4  -3   7
G   0  -2   0   0  -5   6
H   0  -1   0   0   0  -1   6
I   0  -1  -3  -2   1  -4  -2   4
K   0  -2   0   1  -3  -1   0  -2   3
L  -1  -1  -4  -2   2  -4  -1   2  -2   4
M   0   0  -3  -2   1  -3  -1   2  -1   2   4
N   0  -1   2   0  -3   0   1  -2   0  -3  -2   3
P   0  -3   0   0  -3  -1  -1  -2   0  -2  -2   0   7
Q   0  -2   0   1  -2  -1   1  -1   1  -1  -1   0   0   2
R   0  -2   0   0  -3  -1   0  -2   2  -2  -1   0   0   1   4
S   1   0   0   0  -2   0   0  -1   0  -2  -1   0   0   0   0   2
T   0   0   0   0  -2  -1   0   0   0  -1   0   0   0   0   0   1   2
V   0   0  -2  -1   0  -3  -2   3  -1   1   1  -2  -1  -1  -2  -1   0   3
W  -3  -1  -5  -4   3  -4   0  -1  -3   0  -1  -3  -5  -2  -1  -3  -3  -2  14
Y  -2   0  -2  -2   5  -4   2   0  -2   0   0  -1  -3  -1  -1  -1  -1  -1   4   7
   A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
""")
        mat = SubsMat.SeqMat(MatrixInfo.grant)
        self.assertEqual(len(mat), 210)
        self.checkMatrix(mat, """\
A 215
C  20 215
D  89  61 215
E 108  45 170 215
F 102  10  38  75 215
G 155  56 121 117  62 215
H 129  41 134 175 115 117 215
I 121  17  47  81 194  80 121 215
K 109  13 114 159 113  88 183 113 215
L 119  17  43  77 193  77 116 210 108 215
M 131  19  55  89 187  88 128 205 120 200 215
N 104  76 192 173  57 135 147  66 121  62  73 215
P 188  46 107 122 101 173 138 120 112 117 128 124 215
Q 124  61 154 186  99 128 191 106 162 102 114 169 139 215
R 103  35 119 161 118  90 186 118 189 113 124 129 112 172 215
S 116 103 150 135  60 159 126  73  94  70  80 169 141 147 105 215
T 157  66 130 150 112 156 168 126 137 123 134 150 177 173 144 157 215
V 151  23  63  94 165 106 131 186 118 183 194  82 147 119 119  91 146 215
W  67   0  34  63 175  31 100 154 105 154 148  41  68  85 114  38  87 127 215
Y 103  21  55  93 193  68 132 182 130 179 179  72 105 116 138  71 123 160 178 215
   A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
""")
        mat = SubsMat.SeqMat(MatrixInfo.ident)
        self.checkMatrix(mat, """\
A   6
C  -1   6
D  -1  -1   6
E  -1  -1  -1   6
F  -1  -1  -1  -1   6
G  -1  -1  -1  -1  -1   6
H  -1  -1  -1  -1  -1  -1   6
I  -1  -1  -1  -1  -1  -1  -1   6
K  -1  -1  -1  -1  -1  -1  -1  -1   6
L  -1  -1  -1  -1  -1  -1  -1  -1  -1   6
M  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1   6
N  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1   6
P  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1   6
Q  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1   6
R  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1   6
S  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1   6
T  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1   6
V  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1   6
W  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1   6
Y  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1   6
   A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
""")
        mat = SubsMat.SeqMat(MatrixInfo.johnson)
        self.assertEqual(len(mat), 210)
        self.checkMatrix(mat, """\
A   6
C  -3  16
D  -1  -9   8
E   0  -6   2   8
F  -3  -4  -7  -6  10
G   0  -8  -2  -2  -8   8
H  -3  -8   0  -2  -1  -3  12
I  -2  -7  -4  -4   0  -5  -5   8
K   0  -8  -1   1  -5  -3   0  -4   7
L  -3  -8  -8  -5   1  -7  -4   2  -3   7
M  -1  -4  -5  -2   0  -5  -2   2  -1   4  11
N  -1  -7   2   0  -3  -1   1  -4   0  -4  -3   8
P  -1  -8  -1  -1  -5  -2  -4  -5   0  -2  -9  -2  10
Q   0  -6  -1   2  -6  -2   1  -7   1  -4   0   0  -3   9
R  -1  -5  -3   0  -6  -2   0  -5   3  -3  -4  -1  -3   2  10
S   0  -7   0  -2  -4  -1  -2  -4  -1  -5  -4   1  -1  -1   0   5
T   0  -6  -1   0  -5  -3  -3  -3   0  -4  -3   0  -2   0  -1   2   6
V   0  -4  -5  -4  -1  -5  -3   3  -3   1   0  -5  -5  -3  -4  -4  -1   7
W  -5  -9  -6  -7   3  -6  -4  -3  -5  -1   0  -6  -7  -8  -3  -6  -9  -4  15
Y  -4  -7  -3  -3   3  -5   0  -2  -3  -2  -1  -1  -7  -5  -2  -3  -2  -1   2  10
   A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
""")
        mat = SubsMat.SeqMat(MatrixInfo.levin)
        self.assertEqual(len(mat), 210)
        self.checkMatrix(mat, """\
A   2
C   0   2
D   0   0   2
E   1   0   1   2
F  -1  -1  -1  -1   2
G   0   0   0   0  -1   2
H   0   0   0   0  -1   0   2
I   0   0  -1  -1   1  -1  -1   2
K   0   0   0   0  -1   0   0  -1   2
L   0   0  -1  -1   0  -1  -1   0  -1   2
M   0   0  -1  -1   0  -1  -1   0  -1   2   2
N   0   0   1   0  -1   0   0  -1   1  -1  -1   3
P  -1   0   0  -1  -1   0   0  -1   0  -1  -1   0   3
Q   0   0   0   1  -1   0   0  -1   0  -1  -1   1   0   2
R   0   0   0   0  -1   0   0  -1   1  -1  -1   0   0   0   2
S   1   0   0   0  -1   0   0  -1   0  -1  -1   0   0   0   0   2
T   0   0   0   0  -1   0   0   0   0   0   0   0   0   0   0   0   2
V   0   0  -1  -1   0  -1  -1   1  -1   1   0  -1  -1  -1  -1  -1   0   2
W  -1  -1  -1  -1   0  -1  -1   0  -1   0   0  -1  -1  -1   0  -1  -1   0   2
Y  -1  -1  -1  -1   1  -1   0   0  -1   0   0  -1  -1  -1  -1  -1  -1   0   0   2
   A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
""")
        mat = SubsMat.SeqMat(MatrixInfo.mclach)
        self.assertEqual(len(mat), 210)
        self.checkMatrix(mat, """\
A   8
C   1   9
D   3   1   8
E   4   0   5   8
F   1   0   1   0   9
G   3   1   3   3   0   8
H   3   3   4   2   4   2   8
I   2   1   1   1   3   1   2   8
K   3   0   3   4   0   3   4   1   8
L   2   0   1   1   5   1   2   5   2   8
M   3   3   2   1   5   1   3   5   1   6   8
N   3   1   5   4   0   3   4   1   4   1   2   8
P   4   0   3   4   1   3   3   1   3   1   1   1   8
Q   3   0   4   5   0   2   4   0   4   3   3   4   3   8
R   2   1   1   3   1   3   5   1   5   2   1   3   3   5   8
S   4   2   3   4   2   3   3   2   3   2   2   5   3   4   4   8
T   3   2   3   4   1   2   4   3   3   3   3   3   3   3   3   5   8
V   3   1   1   2   3   2   2   5   2   5   4   1   2   2   2   2   3   8
W   1   2   0   1   6   1   3   3   1   3   1   0   0   2   3   3   2   2   9
Y   1   1   1   2   6   0   4   3   1   3   2   2   0   1   2   3   1   3   6   9
   A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
""")
        mat = SubsMat.SeqMat(MatrixInfo.miyata)
        self.assertEqual(len(mat), 210)
        self.checkMatrix(mat, """\
A   1
C   0   1
D  -1  -2   1
E  -1  -2   0   1
F  -1   0  -3  -2   1
G   0   0  -1  -1  -2   1
H   0  -1   0   0  -1  -1   1
I  -1   0  -2  -2   0  -2  -1   1
K  -1  -2   0   0  -1  -2   0  -1   1
L  -1   0  -2  -2   0  -2  -1   1  -1   1
M  -1   0  -2  -1   0  -2   0   0  -1   0   1
N   0  -1   0   0  -2   0   0  -2   0  -2  -1   1
P   1   0  -1  -1  -1   0   0  -1  -1  -1  -1   0   1
Q   0  -1   0   0  -1  -1   0  -1   0  -1  -1   0   0   1
R  -1  -1  -1   0  -1  -2   0  -1   0  -1  -1   0  -1   0   1
S   0   0   0   0  -2   0   0  -1  -1  -1  -1   0   0   0  -1   1
T   0   0   0   0  -1   0   0   0   0  -1   0   0   0   0   0   0   1
V   0   0  -2  -1   0  -1   0   0  -1   0   0  -1   0   0  -1   0   0   1
W  -2  -2  -3  -2   0  -3  -1   0  -1   0   0  -3  -2  -2  -1  -3  -2  -1   1
Y  -1  -1  -2  -1   0  -2  -1   0  -1   0   0  -2  -1  -1   0  -2  -1   0   0   1
   A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
""")
        mat = SubsMat.SeqMat(MatrixInfo.nwsgappep)
        self.assertEqual(len(mat), 253)
        self.checkMatrix(mat, """\
A   1
B   0   1
C   0   0   1
D   0   1   0   1
E   0   0   0   1   1
F   0   0   0  -1   0   1
G   0   0   0   0   0   0   1
H   0   0   0   0   0   0   0   1
I   0   0   0   0   0   0   0   0   1
K   0   0   0   0   0   0   0   0   0   1
L   0   0   0   0   0   1   0   0   0   0   1
M   0   0   0   0   0   0   0   0   0   0   1   1
N   0   1   0   0   0   0   0   0   0   0   0   0   1
P   0   0   0   0   0   0   0   0   0   0   0   0   0   1
Q   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1
R   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1
S   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1
T   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1
V   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   1
W   0   0  -1  -1  -1   1  -1   0   0   0   0   0   0   0   0   1   0   0   0   1
Y   0   0   1   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1
Z   0   0   0   0   1   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   1
   A   B   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y   Z
""")
        mat = SubsMat.SeqMat(MatrixInfo.pam120)
        self.assertEqual(len(mat), 276)
        self.checkMatrix(mat, """\
A   3
B   0   4
C  -3  -6   9
D   0   4  -7   5
E   0   3  -7   3   5
F  -4  -5  -6  -7  -7   8
G   1   0  -4   0  -1  -5   5
H  -3   1  -4   0  -1  -3  -4   7
I  -1  -3  -3  -3  -3   0  -4  -4   6
K  -2   0  -7  -1  -1  -7  -3  -2  -3   5
L  -3  -4  -7  -5  -4   0  -5  -3   1  -4   5
M  -2  -4  -6  -4  -3  -1  -4  -4   1   0   3   8
N  -1   3  -5   2   1  -4   0   2  -2   1  -4  -3   4
P   1  -2  -4  -3  -2  -5  -2  -1  -3  -2  -3  -3  -2   6
Q  -1   0  -7   1   2  -6  -3   3  -3   0  -2  -1   0   0   6
R  -3  -2  -4  -3  -3  -5  -4   1  -2   2  -4  -1  -1  -1   1   6
S   1   0   0   0  -1  -3   1  -2  -2  -1  -4  -2   1   1  -2  -1   3
T   1   0  -3  -1  -2  -4  -1  -3   0  -1  -3  -1   0  -1  -2  -2   2   4
V   0  -3  -3  -3  -3  -3  -2  -3   3  -4   1   1  -3  -2  -3  -3  -2   0   5
W  -7  -6  -8  -8  -8  -1  -8  -3  -6  -5  -3  -6  -4  -7  -6   1  -2  -6  -8  12
X  -1  -1  -4  -2  -1  -3  -2  -2  -1  -2  -2  -2  -1  -2  -1  -2  -1  -1  -1  -5  -2
Y  -4  -3  -1  -5  -5   4  -6  -1  -2  -5  -2  -4  -2  -6  -5  -5  -3  -3  -3  -2  -3   8
Z  -1   2  -7   3   4  -6  -2   1  -3  -1  -3  -2   0  -1   4  -1  -1  -2  -3  -7  -1  -5   4
   A   B   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   X   Y   Z
""")
        mat = SubsMat.SeqMat(MatrixInfo.pam180)
        self.assertEqual(len(mat), 276)
        self.checkMatrix(mat, """\
A   3
B   0   4
C  -3  -6  13
D   0   4  -7   5
E   0   3  -7   4   5
F  -5  -6  -6  -8  -7  10
G   1   0  -5   0   0  -6   6
H  -2   1  -4   0   0  -3  -3   8
I  -1  -3  -3  -3  -3   1  -4  -4   6
K  -2   0  -7   0  -1  -7  -3  -1  -3   6
L  -3  -5  -8  -6  -5   1  -6  -3   2  -4   7
M  -2  -3  -7  -4  -3   0  -4  -3   2   1   4   9
N   0   3  -5   3   2  -5   0   2  -3   1  -4  -3   4
P   1  -2  -4  -2  -1  -6  -1  -1  -3  -2  -4  -3  -1   8
Q  -1   1  -7   2   3  -6  -2   4  -3   0  -2  -1   0   0   6
R  -3  -2  -5  -3  -2  -6  -4   2  -3   4  -4  -1  -1  -1   1   8
S   1   1   0   0  -1  -4   1  -2  -2  -1  -4  -2   1   1  -1  -1   3
T   2   0  -3  -1  -1  -4  -1  -2   0   0  -3  -1   0   0  -2  -2   2   4
V   0  -3  -3  -3  -3  -2  -2  -3   5  -4   2   2  -3  -2  -3  -4  -2   0   6
W  -8  -7 -10  -9  -9   0  -9  -4  -7  -5  -3  -6  -5  -7  -6   2  -3  -7  -8  18
X  -1  -1  -4  -1  -1  -3  -2  -1  -1  -1  -2  -1  -1  -1  -1  -2   0  -1  -1  -6  -1
Y  -5  -4   0  -6  -6   7  -7   0  -2  -6  -2  -4  -2  -7  -6  -6  -4  -4  -4  -1  -3  11
Z   0   3  -7   3   5  -7  -1   2  -3   0  -3  -2   1  -1   5   0  -1  -1  -3  -8  -1  -6   5
   A   B   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   X   Y   Z
""")
        mat = SubsMat.SeqMat(MatrixInfo.pam250)
        self.assertEqual(len(mat), 276)
        self.checkMatrix(mat, """\
A   2
B   0   3
C  -2  -4  12
D   0   3  -5   4
E   0   3  -5   3   4
F  -3  -4  -4  -6  -5   9
G   1   0  -3   1   0  -5   5
H  -1   1  -3   1   1  -2  -2   6
I  -1  -2  -2  -2  -2   1  -3  -2   5
K  -1   1  -5   0   0  -5  -2   0  -2   5
L  -2  -3  -6  -4  -3   2  -4  -2   2  -3   6
M  -1  -2  -5  -3  -2   0  -3  -2   2   0   4   6
N   0   2  -4   2   1  -3   0   2  -2   1  -3  -2   2
P   1  -1  -3  -1  -1  -5   0   0  -2  -1  -3  -2   0   6
Q   0   1  -5   2   2  -5  -1   3  -2   1  -2  -1   1   0   4
R  -2  -1  -4  -1  -1  -4  -3   2  -2   3  -3   0   0   0   1   6
S   1   0   0   0   0  -3   1  -1  -1   0  -3  -2   1   1  -1   0   2
T   1   0  -2   0   0  -3   0  -1   0   0  -2  -1   0   0  -1  -1   1   3
V   0  -2  -2  -2  -2  -1  -1  -2   4  -2   2   2  -2  -1  -2  -2  -1   0   4
W  -6  -5  -8  -7  -7   0  -7  -3  -5  -3  -2  -4  -4  -6  -5   2  -2  -5  -6  17
X   0  -1  -3  -1  -1  -2  -1  -1  -1  -1  -1  -1   0  -1  -1  -1   0   0  -1  -4  -1
Y  -3  -3   0  -4  -4   7  -5   0  -1  -4  -1  -2  -2  -5  -4  -4  -3  -3  -2   0  -2  10
Z   0   2  -5   3   3  -5   0   2  -2   0  -3  -2   1   0   3   0   0  -1  -2  -6  -1  -4   3
   A   B   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   X   Y   Z
""")
        mat = SubsMat.SeqMat(MatrixInfo.pam30)
        self.assertEqual(len(mat), 276)
        self.checkMatrix(mat, """\
A   6
B  -3   6
C  -6 -12  10
D  -3   6 -14   8
E  -2   1 -14   2   8
F  -8 -10 -13 -15 -14   9
G  -2  -3  -9  -3  -4  -9   6
H  -7  -1  -7  -4  -5  -6  -9   9
I  -5  -6  -6  -7  -5  -2 -11  -9   8
K  -7  -2 -14  -4  -4 -14  -7  -6  -6   7
L  -6  -9 -15 -12  -9  -3 -10  -6  -1  -8   7
M  -5 -10 -13 -11  -7  -4  -8 -10  -1  -2   1  11
N  -4   6 -11   2  -2  -9  -3   0  -5  -1  -7  -9   8
P  -2  -7  -8  -8  -5 -10  -6  -4  -8  -6  -7  -8  -6   8
Q  -4  -3 -14  -2   1 -13  -7   1  -8  -3  -5  -4  -3  -3   8
R  -7  -7  -8 -10  -9  -9  -9  -2  -5   0  -8  -4  -6  -4  -2   8
S   0  -1  -3  -4  -4  -6  -2  -6  -7  -4  -8  -5   0  -2  -5  -3   6
T  -1  -3  -8  -5  -6  -9  -6  -7  -2  -3  -7  -4  -2  -4  -5  -6   0   7
V  -2  -8  -6  -8  -6  -8  -5  -6   2  -9  -2  -1  -8  -6  -7  -8  -6  -3   7
W -13 -10 -15 -15 -17  -4 -15  -7 -14 -12  -6 -13  -8 -14 -13  -2  -5 -13 -15  13
X  -3  -5  -9  -5  -5  -8  -5  -5  -5  -5  -6  -5  -3  -5  -5  -6  -3  -4  -5 -11  -5
Y  -8  -6  -4 -11  -8   2 -14  -3  -6  -9  -7 -11  -4 -13 -12 -10  -7  -6  -7  -5  -7  10
Z  -3   0 -14   1   6 -13  -5  -1  -6  -4  -7  -5  -3  -4   6  -4  -5  -6  -6 -14  -5  -9   6
   A   B   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   X   Y   Z
""")
        mat = SubsMat.SeqMat(MatrixInfo.pam300)
        self.assertEqual(len(mat), 276)
        self.checkMatrix(mat, """\
A   2
B   0   3
C  -2  -5  15
D   0   3  -6   4
E   0   3  -6   4   4
F  -4  -5  -5  -6  -6  11
G   2   1  -4   1   0  -5   5
H  -1   1  -4   1   1  -2  -2   7
I   0  -2  -3  -2  -2   1  -3  -2   5
K  -1   1  -6   0   0  -6  -2   0  -2   5
L  -2  -4  -7  -4  -4   3  -4  -2   3  -3   7
M  -1  -2  -6  -3  -2   1  -3  -2   3   0   4   6
N   0   2  -4   2   2  -4   1   2  -2   1  -3  -2   2
P   1   0  -3  -1   0  -5   0   0  -2  -1  -3  -2   0   6
Q   0   2  -6   2   3  -5  -1   3  -2   1  -2  -1   1   0   4
R  -1   0  -4  -1  -1  -5  -2   2  -2   4  -3   0   0   0   2   7
S   1   1   0   0   0  -4   1  -1  -1   0  -3  -2   1   1   0   0   1
T   1   0  -2   0   0  -3   0  -1   0   0  -2  -1   0   1  -1  -1   1   2
V   0  -2  -2  -2  -2  -1  -1  -2   4  -2   2   2  -2  -1  -2  -3  -1   0   5
W  -6  -6  -9  -7  -8   1  -8  -3  -6  -4  -2  -5  -5  -6  -5   3  -3  -6  -7  22
X   0   0  -3  -1  -1  -2  -1   0  -1  -1  -1  -1   0  -1   0  -1   0   0   0  -4  -1
Y  -4  -4   1  -5  -5   9  -6   0  -1  -5   0  -2  -2  -5  -4  -5  -3  -3  -3   0  -2  12
Z   0   2  -6   3   3  -5   0   2  -2   1  -3  -2   1   0   3   0   0   0  -2  -6  -1  -5   3
   A   B   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   X   Y   Z
""")
        mat = SubsMat.SeqMat(MatrixInfo.pam60)
        self.assertEqual(len(mat), 276)
        self.checkMatrix(mat, """\
A   5
B  -2   5
C  -5  -9   9
D  -2   5 -10   7
E  -1   2 -10   3   7
F  -6  -8  -9 -11 -10   8
G   0  -2  -7  -2  -2  -7   6
H  -5   0  -6  -2  -3  -4  -6   8
I  -3  -4  -4  -5  -4  -1  -7  -6   7
K  -5  -1 -10  -2  -3 -10  -5  -4  -4   6
L  -4  -7 -11  -9  -7  -1  -8  -4   0  -6   6
M  -3  -6 -10  -7  -5  -2  -6  -7   1   0   2  10
N  -2   5  -7   2   0  -6  -1   1  -4   0  -5  -6   6
P   0  -4  -6  -5  -3  -7  -4  -2  -6  -4  -5  -6  -4   7
Q  -3  -1 -10  -1   2  -9  -5   2  -5  -1  -3  -2  -2  -1   7
R  -5  -5  -6  -6  -6  -7  -7   0  -4   2  -6  -2  -3  -2   0   8
S   1   0  -1  -2  -2  -5   0  -4  -4  -2  -6  -4   1   0  -3  -2   5
T   1  -2  -5  -3  -4  -6  -3  -5  -1  -2  -5  -2  -1  -2  -4  -4   1   6
V  -1  -5  -4  -6  -4  -5  -4  -5   3  -6  -1   0  -5  -4  -5  -5  -4  -1   6
W -10  -8 -12 -11 -12  -3 -11  -5 -10  -8  -4  -9  -6 -10  -9   0  -4  -9 -11  13
X  -2  -3  -6  -3  -3  -5  -3  -3  -3  -3  -4  -3  -2  -3  -3  -4  -2  -2  -3  -8  -3
Y  -6  -5  -2  -8  -7   3 -10  -2  -4  -7  -5  -7  -3 -10  -8  -8  -5  -5  -5  -3  -5   9
Z  -2   1 -10   2   5 -10  -3   0  -4  -2  -5  -4  -1  -2   6  -2  -3  -4  -5 -11  -3  -7   5
   A   B   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   X   Y   Z
""")
        mat = SubsMat.SeqMat(MatrixInfo.pam90)
        self.assertEqual(len(mat), 253)
        self.checkMatrix(mat, """\
A   4
B  -1   4
C  -3  -7   9
D  -1   5  -8   6
E   0   2  -8   4   6
F  -5  -6  -7  -8  -8   8
G   0  -1  -5  -1  -1  -6   5
H  -4   1  -5  -1  -1  -3  -5   8
I  -2  -3  -3  -4  -3   0  -5  -5   6
K  -3   0  -8  -2  -2  -8  -4  -2  -3   5
L  -3  -5  -9  -7  -5   0  -6  -3   1  -5   6
M  -2  -5  -8  -5  -4  -1  -5  -5   1   0   2   9
N  -1   4  -6   3   0  -5  -1   2  -3   1  -4  -4   5
P   0  -3  -5  -4  -2  -6  -3  -2  -4  -3  -4  -4  -2   7
Q  -2   0  -8   0   2  -7  -3   2  -4  -1  -3  -2  -1  -1   6
R  -4  -3  -5  -5  -4  -6  -5   1  -3   2  -5  -2  -2  -1   0   7
S   1   0  -1  -1  -2  -4   0  -3  -3  -1  -5  -3   1   0  -2  -1   4
T   1  -1  -4  -2  -2  -5  -2  -3   0  -1  -3  -2   0  -1  -3  -3   2   5
V   0  -4  -3  -4  -3  -4  -3  -4   3  -5   0   1  -4  -3  -4  -4  -3  -1   6
W  -8  -7 -10  -9 -10  -2  -9  -4  -8  -6  -3  -7  -5  -8  -7   0  -3  -7  -9  13
Y  -5  -4  -1  -6  -6   4  -8  -1  -3  -6  -3  -6  -2  -8  -6  -6  -4  -4  -4  -2   9
Z  -1   2  -8   3   5  -8  -2   1  -3  -1  -4  -3   0  -2   5  -1  -2  -2  -3  -8  -6   5
   A   B   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y   Z
""")
        mat = SubsMat.SeqMat(MatrixInfo.rao)
        self.assertEqual(len(mat), 210)
        self.checkMatrix(mat, """\
A  16
C  11  16
D   9   8  16
E  10   9  11  16
F  10  10   4   6  16
G   8   8   9   6   7  16
H  11  10   9  11   9   7  16
I   9   8   3   4  12   6   8  16
K  10   9  11  11   6   7  11   4  16
L  11  11   6   7  11   6  10  10   7  16
M  11  10   5   8  10   4  10   9   8  11  16
N   9   9  11  10   6  10  10   5  11   7   6  16
P   6   7   8   5   4  11   5   3   6   4   2   9  16
Q  11  10  11  11   7   8  11   6  12   9   9  11   7  16
R   8   8  10   9   5   7  10   4  11   6   6  10   6  10  16
S  10  10  10   9   8  11  10   8  10   8   7  11  10  10   9  16
T  10  10   9   8  10  10  10  10   9   9   8  10   8  10   9  11  16
V   9   8   3   4  11   6   9  12   5  10   9   5   3   6   5   8  10  16
W  11  11   6   7  11   8  10  11   7  11  10   8   6   9   7  10  11  11  16
Y   9  10   7   6  10  10   9  10   7   9   8   8   8   8   7  11  11  10  11  16
   A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
""")
        mat = SubsMat.SeqMat(MatrixInfo.risler)
        self.assertEqual(len(mat), 210)
        self.checkMatrix(mat, """\
A   2
C  -1   2
D   0  -1   2
E   1  -1   1   2
F   0  -1   0   0   2
G   0  -1   0   0   0   2
H   0  -1  -1   0  -1  -1   2
I   1  -1   0   1   1   0   0   2
K   1  -1   0   1   0   0  -1   1   2
L   1  -1   0   0   1   0   0   2   0   2
M   1  -1   0   0   0   0  -1   0   0   1   2
N   1  -1   0   1   0   0   0   0   1   0   0   2
P   0  -1  -1   0  -1  -1  -1   0   0   0  -1  -1   2
Q   1  -1   0   2   0   0   0   1   1   1   1   1   0   2
R   1  -1   0   1   0   0   0   1   2   1   1   1   0   2   2
S   2  -1   0   1   0   0   0   1   1   1   0   1   0   1   2   2
T   1  -1   0   1   0   0   0   1   1   1   0   1   0   1   1   2   2
V   2  -1   0   1   0   0   0   2   1   2   0   1   0   1   1   1   1   2
W   0  -1  -1  -1   0  -1  -1   0  -1   0  -1  -1  -1  -1   0   0  -1   0   2
Y   0  -1   0   0   2   0   0   0   0   0   0   0  -1   0   0   0   0   0   0   2
   A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
""")
        mat = SubsMat.SeqMat(MatrixInfo.structure)
        self.assertEqual(len(mat), 210)
        self.checkMatrix(mat, """\
A   4
C  -2  11
D  -1  -7   6
E   0  -3   2   5
F  -3  -2  -5  -4   7
G   0  -6  -1  -2  -6   5
H  -2  -6   0  -2  -2  -3   8
I  -2  -4  -3  -3   1  -5  -5   6
K  -1  -4  -1   1  -3  -3   0  -3   5
L  -2  -6  -6  -4   2  -5  -3   2  -2   5
M   0  -5  -4  -2   0  -4  -2   1  -1   3   8
N  -1  -6   2   0  -3  -1   2  -3   0  -3  -2   5
P  -1  -8  -1  -1  -5  -2  -3  -4  -1  -3  -6  -2   7
Q   0  -3   0   2  -4  -2   0  -5   1  -3   1   0  -2   6
R  -1  -2  -2   0  -4  -2   0  -3   2  -3  -4  -1  -2   1   7
S   0  -4   0  -1  -3  -1  -2  -3  -1  -4  -4   0  -1  -1   0   4
T  -1  -5  -1   0  -3  -3  -2  -2   0  -3  -2   0  -1   0  -1   1   5
V   0  -4  -4  -2  -1  -4  -2   2  -3   1   0  -4  -4  -2  -3  -3  -1   5
W  -3  -6  -6  -6   2  -4  -3  -2  -3  -1  -2  -5  -4  -5  -2  -5  -5  -4  10
Y  -3  -6  -3  -2   3  -3   0  -1  -2  -2  -1  -1  -6  -3  -1  -2  -2  -1   2   7
   A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
""")

    def test_matrix_correlations(self):
        blosum90 = SubsMat.SeqMat(MatrixInfo.blosum90)
        blosum30 = SubsMat.SeqMat(MatrixInfo.blosum30)
        correlation = SubsMat.two_mat_correlation(blosum30, blosum90)
        self.assertAlmostEqual(correlation, 0.878, places=3)
        correlation = SubsMat.two_mat_correlation(blosum90, blosum30)
        self.assertAlmostEqual(correlation, 0.878, places=3)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
