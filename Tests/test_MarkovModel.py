# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

try:
    from numpy import array
    from numpy import random #missing in PyPy's micronumpy
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(\
        "Install NumPy if you want to use Bio.MarkovModel.")

import unittest

import warnings
#Silence this warning:
#For optimal speed, please update to Numpy version 1.3 or later
warnings.filterwarnings("ignore", category=UserWarning)
from Bio import MarkovModel
warnings.filters.pop()


class TestMarkovModel(unittest.TestCase):

    def test_train_visible(self):
        states = ["0", "1", "2", "3"]
        alphabet = ["A", "C", "G", "T"]
        training_data = [
            ("AACCCGGGTTTTTTT", "001112223333333"),
            ("ACCGTTTTTTT", "01123333333"),
            ("ACGGGTTTTTT", "01222333333"),
            ("ACCGTTTTTTTT", "011233333333"),
            ]
        markov_model = MarkovModel.train_visible(states, alphabet, training_data)
        states = MarkovModel.find_states(markov_model, "AACGTT")
        self.assertEqual(len(states), 1)
        state_list, state_float = states[0]
        self.assertEqual(state_list, ['0', '0', '1', '2', '3', '3'])
        self.assertAlmostEqual(state_float, 0.0082128906)
        self.assertEqual(markov_model.states, ['0', '1', '2', '3'])
        self.assertEqual(markov_model.alphabet, ['A', 'C', 'G', 'T'])
        self.assertEqual(len(markov_model.p_initial), 4)
        self.assertAlmostEqual(markov_model.p_initial[0], 1.0)
        self.assertAlmostEqual(markov_model.p_initial[1], 0.0)
        self.assertAlmostEqual(markov_model.p_initial[2], 0.0)
        self.assertAlmostEqual(markov_model.p_initial[3], 0.0)
        self.assertEqual(len(markov_model.p_transition), 4)
        self.assertEqual(len(markov_model.p_transition[0]), 4)
        self.assertEqual(len(markov_model.p_transition[1]), 4)
        self.assertEqual(len(markov_model.p_transition[2]), 4)
        self.assertEqual(len(markov_model.p_transition[3]), 4)
        self.assertAlmostEqual(markov_model.p_transition[0][0], 0.2)
        self.assertAlmostEqual(markov_model.p_transition[0][1], 0.8)
        self.assertAlmostEqual(markov_model.p_transition[0][2], 0.0)
        self.assertAlmostEqual(markov_model.p_transition[0][3], 0.0)
        self.assertAlmostEqual(markov_model.p_transition[1][0], 0.0)
        self.assertAlmostEqual(markov_model.p_transition[1][1], 0.5)
        self.assertAlmostEqual(markov_model.p_transition[1][2], 0.5)
        self.assertAlmostEqual(markov_model.p_transition[1][3], 0.0)
        self.assertAlmostEqual(markov_model.p_transition[2][0], 0.0)
        self.assertAlmostEqual(markov_model.p_transition[2][1], 0.0)
        self.assertAlmostEqual(markov_model.p_transition[2][2], 0.5)
        self.assertAlmostEqual(markov_model.p_transition[2][3], 0.5)
        self.assertAlmostEqual(markov_model.p_transition[3][0], 0.0)
        self.assertAlmostEqual(markov_model.p_transition[3][1], 0.0)
        self.assertAlmostEqual(markov_model.p_transition[3][2], 0.0)
        self.assertAlmostEqual(markov_model.p_transition[3][3], 1.0)
        self.assertEqual(len(markov_model.p_emission), 4)
        self.assertEqual(len(markov_model.p_emission[0]), 4)
        self.assertEqual(len(markov_model.p_emission[1]), 4)
        self.assertEqual(len(markov_model.p_emission[2]), 4)
        self.assertEqual(len(markov_model.p_emission[3]), 4)
        self.assertAlmostEqual(markov_model.p_emission[0][0], 0.666667,
                               places=4)
        self.assertAlmostEqual(markov_model.p_emission[0][1], 0.111111,
                               places=4)
        self.assertAlmostEqual(markov_model.p_emission[0][2], 0.111111,
                               places=4)
        self.assertAlmostEqual(markov_model.p_emission[0][3], 0.111111,
                               places=4)
        self.assertAlmostEqual(markov_model.p_emission[1][0], 0.083333,
                               places=4)
        self.assertAlmostEqual(markov_model.p_emission[1][1], 0.750000,
                               places=4)
        self.assertAlmostEqual(markov_model.p_emission[1][2], 0.083333,
                               places=4)
        self.assertAlmostEqual(markov_model.p_emission[1][3], 0.083333,
                               places=4)
        self.assertAlmostEqual(markov_model.p_emission[2][0], 0.083333,
                               places=4)
        self.assertAlmostEqual(markov_model.p_emission[2][1], 0.083333,
                               places=4)
        self.assertAlmostEqual(markov_model.p_emission[2][2], 0.750000,
                               places=4)
        self.assertAlmostEqual(markov_model.p_emission[2][3], 0.083333,
                               places=4)
        self.assertAlmostEqual(markov_model.p_emission[3][0], 0.031250,
                               places=4)
        self.assertAlmostEqual(markov_model.p_emission[3][1], 0.031250,
                               places=4)
        self.assertAlmostEqual(markov_model.p_emission[3][2], 0.031250,
                               places=4)
        self.assertAlmostEqual(markov_model.p_emission[3][3], 0.906250,
                               places=4)

    def test_baum_welch(self):
        states = ["CP", "IP"]
        alphabet = ["cola", "ice_t", "lem"]
        outputs = [
            (2, 1, 0)
            ]
        p_initial = [1.0, 0.0000001]
        p_transition = [[0.7, 0.3],
                        [0.5, 0.5]]
        p_emission = [[0.6, 0.1, 0.3],
                      [0.1, 0.7, 0.2]]
        N, M = len(states), len(alphabet)
        x = MarkovModel._baum_welch(N, M, outputs,
                                    p_initial=p_initial,
                                    p_transition=p_transition,
                                    p_emission=p_emission
                                    )
        p_initial, p_transition, p_emission = x
        markov_model = MarkovModel.MarkovModel(states, alphabet,
                                     p_initial, p_transition, p_emission)
        self.assertEqual(markov_model.states, ['CP', 'IP'])
        self.assertEqual(markov_model.alphabet, ['cola', 'ice_t', 'lem'])
        self.assertEqual(len(markov_model.p_initial), 2)
        self.assertAlmostEqual(markov_model.p_initial[0], 1.0,
                               places=4)
        self.assertAlmostEqual(markov_model.p_initial[1], 0.0,
                               places=4)
        self.assertEqual(len(markov_model.p_transition), 2)
        self.assertEqual(len(markov_model.p_transition[0]), 2)
        self.assertEqual(len(markov_model.p_transition[1]), 2)
        self.assertAlmostEqual(markov_model.p_transition[0][0], 0.02460365,
                               places=4)
        self.assertAlmostEqual(markov_model.p_transition[0][1], 0.97539634,
                               places=4)
        self.assertAlmostEqual(markov_model.p_transition[1][0], 1.0,
                               places=4)
        self.assertAlmostEqual(markov_model.p_transition[1][1], 0.0,
                               places=4)
        self.assertEqual(len(markov_model.p_emission), 2)
        self.assertEqual(len(markov_model.p_emission[0]), 3)
        self.assertEqual(len(markov_model.p_emission[1]), 3)
        self.assertAlmostEqual(markov_model.p_emission[0][0], 0.5)
        self.assertAlmostEqual(markov_model.p_emission[0][1], 0.0)
        self.assertAlmostEqual(markov_model.p_emission[0][2], 0.5)
        self.assertAlmostEqual(markov_model.p_emission[1][0], 0.0)
        self.assertAlmostEqual(markov_model.p_emission[1][1], 1.0)
        self.assertAlmostEqual(markov_model.p_emission[1][2], 0.0)

    # Do some tests from the topcoder competition.

    def test_topcoder1(self):
        # NNNN
        states = "NR"
        alphabet = "AGTC"
        p_initial = array([1.0, 0.0])
        p_transition = array([[0.90, 0.10],
                              [0.20, 0.80]])
        p_emission = array([[0.30, 0.20, 0.30, 0.20],
                            [0.10, 0.40, 0.10, 0.40]])
        markov_model = MarkovModel.MarkovModel(
            states, alphabet, p_initial, p_transition, p_emission)
        states = MarkovModel.find_states(markov_model, "TGCC")
        self.assertEqual(len(states), 1)
        state_list, state_float = states[0]
        self.assertEqual(state_list, ['N', 'N', 'N', 'N'])

    def test_topcoder2(self):
        # NNNRRRNNRRNRRN
        states = "NR"
        alphabet = "AGTC"
        p_initial = array([1.0, 0.0])
        p_transition = array([[0.56, 0.44],
                              [0.25, 0.75]])
        p_emission = array([[0.04, 0.14, 0.62, 0.20],
                            [0.39, 0.15, 0.04, 0.42]])
        markov_model = MarkovModel.MarkovModel(
            states, alphabet, p_initial, p_transition, p_emission)
        states = MarkovModel.find_states(markov_model, "CCTGAGTTAGTCGT")
        self.assertEqual(len(states), 1)
        state_list, state_float = states[0]
        self.assertEqual(state_list, ['N', 'N', 'N', 'R', 'R', 'R', 'N', 'N', 'R', 'R', 'N', 'R', 'R', 'N'])

    def test_topcoder3(self):
        # NRRRRRRRRRRRNNNNRRRRRRRRR
        states = "NR"
        alphabet = "AGTC"
        p_initial = array([1.0, 0.0])
        p_transition = array([[0.75, 0.25],
                              [0.25, 0.75]])
        p_emission = array([[0.45, 0.36, 0.06, 0.13],
                            [0.24, 0.18, 0.12, 0.46]])
        markov_model = MarkovModel.MarkovModel(
            states, alphabet, p_initial, p_transition, p_emission)
        states = MarkovModel.find_states(markov_model, "CCGTACTTACCCAGGACCGCAGTCC")
        self.assertEqual(len(states), 1)
        state_list, state_float = states[0]
        self.assertEqual(state_list, ['N', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'N', 'N', 'N', 'N', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R'])

    def test_topcoder4(self):
        # NRRRRRRRRRR
        states = "NR"
        alphabet = "AGTC"
        p_initial = array([1.0, 0.0])
        p_transition = array([[0.55, 0.45],
                              [0.15, 0.85]])
        p_emission = array([[0.75, 0.03, 0.01, 0.21],
                            [0.34, 0.11, 0.39, 0.16]])
        markov_model = MarkovModel.MarkovModel(
            states, alphabet, p_initial, p_transition, p_emission)
        states = MarkovModel.find_states(markov_model, "TTAGCAGTGCG")
        self.assertEqual(len(states), 1)
        state_list, state_float = states[0]
        self.assertEqual(state_list, ['N','R','R','R','R','R','R','R','R','R','R'])

    def test_topcoder5(self):
        # N
        states = "NR"
        alphabet = "AGTC"
        p_initial = array([1.0, 0.0])
        p_transition = array([[0.84, 0.16],
                              [0.25, 0.75]])
        p_emission = array([[0.26, 0.37, 0.08, 0.29],
                            [0.31, 0.13, 0.33, 0.23]])
        markov_model = MarkovModel.MarkovModel(
            states, alphabet, p_initial, p_transition, p_emission)
        states = MarkovModel.find_states(markov_model, "T")
        self.assertEqual(len(states), 1)
        state_list, state_float = states[0]
        self.assertEqual(state_list, ["N"])


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
