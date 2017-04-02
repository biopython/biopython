# Copyright 2002 by Jeffrey Chang.  All rights reserved.
# Revisions copyright 2008 by Brad Chapman. All rights reserved.
# Revisions copyright 2008 by Michiel de Hoon. All rights reserved.
# Revisions copyright 2008-2010,2013-2014 by Peter Cock. All rights reserved.
# Revisions copyright 2012 by Christian Brueffer. All rights reserved.
# Revisions copyright 2017 by Maximilian Greil. All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

import warnings
import unittest

from Bio._py3k import StringIO

try:
    from numpy import array
    from numpy import random  # missing in PyPy's micronumpy
    from numpy import array_equal
    from numpy import around
    from numpy import log
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.MarkovModel.")

with warnings.catch_warnings():
    # Silence this warning:
    # For optimal speed, please update to Numpy version 1.3 or later
    warnings.simplefilter("ignore", UserWarning)
    from Bio import MarkovModel


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
        self.assertEqual(state_list, ['N', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R'])

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

    def test_readline_and_check_start(self):
        states = "NR"
        alphabet = "AGTC"
        markov_model = MarkovModel.MarkovModel(states, alphabet)

        line = "This is a \n string with two lines \n"
        handle = StringIO(line)
        start = "This is a \n"
        self.assertEqual(
            start, MarkovModel._readline_and_check_start(handle, start))

    def test_save_and_load(self):
        states = "NR"
        alphabet = "AGTC"
        p_initial = array([1.0, 0.0])
        p_transition = array([[0.75, 0.25], [0.25, 0.75]])
        p_emission = array(
            [[0.45, 0.36, 0.06, 0.13], [0.24, 0.18, 0.12, 0.46]])
        markov_model_save = MarkovModel.MarkovModel(
            states,
            alphabet,
            p_initial,
            p_transition,
            p_emission)

        handle = StringIO()
        MarkovModel.save(markov_model_save, handle)
        handle.seek(0)
        markov_model_load = MarkovModel.load(handle)

        self.assertEqual(''.join(markov_model_load.states), states)
        self.assertEqual(''.join(markov_model_load.alphabet), alphabet)
        self.assertTrue(array_equal(markov_model_load.p_initial, p_initial))
        self.assertTrue(array_equal
                        (markov_model_load.p_transition, p_transition))
        self.assertTrue(array_equal(markov_model_load.p_emission, p_emission))

    def test_train_bw(self):
        random.seed(0)
        states = ["0", "1", "2", "3"]
        alphabet = ["A", "C", "G", "T"]
        training_data = ["AACCCGGGTTTTTTT", "ACCGTTTTTTT",
                         "ACGGGTTTTTT", "ACCGTTTTTTTT"]

        output_p_initial = array([0.2275677, 0.29655611,
                                  0.24993822, 0.22593797])
        output_p_transition = array(
            [[5.16919807e-001, 3.65825814e-033, 4.83080193e-001, 9.23220689e-042],
             [3.65130247e-001,
              1.00000000e-300,
              6.34869753e-001,
              1.00000000e-300],
             [8.68776164e-001,
              1.02254304e-034,
              1.31223836e-001,
              6.21835051e-047],
             [3.33333333e-301, 3.33333333e-001, 3.33333333e-301, 6.66666667e-001]])
        output_p_emission = array(
            [[2.02593570e-301, 2.02593570e-301, 2.02593570e-301, 1.00000000e+000],
             [1.00000000e-300,
              1.00000000e-300,
              1.00000000e+000,
              1.09629016e-259],
             [3.26369779e-301,
              3.26369779e-301,
              3.26369779e-301,
              1.00000000e+000],
             [3.33333333e-001, 6.66666667e-001, 3.33333333e-301, 3.33333333e-301]])

        markov_model = MarkovModel.train_bw(states, alphabet, training_data)
        self.assertEqual(''.join(markov_model.states), ''.join(states))
        self.assertEqual(''.join(markov_model.alphabet), ''.join(alphabet))
        self.assertTrue(array_equal(
            around(markov_model.p_initial, decimals=3),
            around(output_p_initial, decimals=3)))
        self.assertTrue(array_equal(around(
            markov_model.p_transition, decimals=3),
            around(output_p_transition, decimals=3)))
        self.assertTrue(array_equal(around(
            markov_model.p_emission, decimals=3),
            around(output_p_emission, decimals=3)))

    def test_forward(self):
        states = ["CP", "IP"]
        outputs = [2, 1, 0]
        lp_initial = log([1.0, 0.0000001])
        lp_transition = log([[0.7, 0.3], [0.5, 0.5]])
        lp_emission = log([[0.6, 0.1, 0.3], [0.1, 0.7, 0.2]])

        matrix = array([[0., -1.5606477, -3.07477539, -3.84932984],
                        [-16.11809565, -2.4079455, -3.27544608, -4.5847794]])
        self.assertTrue(
            array_equal(around(MarkovModel._forward(len(states), len(outputs),
                                                    lp_initial,
                                                    lp_transition,
                                                    lp_emission,
                                                    outputs), decimals=3),
                        around(matrix, decimals=3)))

    def test_backward(self):
        states = ["CP", "IP"]
        outputs = [2, 1, 0]
        lp_transition = log([[0.7, 0.3], [0.5, 0.5]])
        lp_emission = log([[0.6, 0.1, 0.3], [0.1, 0.7, 0.2]])

        matrix = array([[-3.45776773, -3.10109279, -0.51082562, 0.],
                        [-3.54045945, -1.40649707, -2.30258509, 0.]])
        self.assertTrue(
            array_equal(around(MarkovModel._backward(
                len(states), len(outputs), lp_transition, lp_emission, outputs), decimals=3),
                around(matrix, decimals=3)))

    def test_mle(self):
        states = ["0", "1", "2", "3"]
        alphabet = ["A", "C", "G", "T"]
        training_data = [("AACCCGGGTTTTTTT", "001112223333333"),
                         ("ACCGTTTTTTT", "01123333333"),
                         ("ACGGGTTTTTT", "01222333333"),
                         ("ACCGTTTTTTTT", "011233333333"), ]
        training_outputs = array([[0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3], [
                                 0, 1, 1, 2, 3, 3, 3, 3, 3, 3, 3], [0, 1, 2, 2, 2, 3, 3, 3, 3, 3, 3], [0, 1, 1, 2, 3, 3, 3, 3, 3, 3, 3, 3]])
        training_states = array([[0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3], [
                                0, 1, 1, 2, 3, 3, 3, 3, 3, 3, 3], [0, 1, 2, 2, 2, 3, 3, 3, 3, 3, 3], [0, 1, 1, 2, 3, 3, 3, 3, 3, 3, 3, 3]])

        p_initial = array([1., 0., 0., 0.])
        p_transition = array([[0.2, 0.8, 0., 0.],
                              [0., 0.5, 0.5, 0.],
                              [0., 0., 0.5, 0.5],
                              [0., 0., 0., 1.]])
        p_emission = array(
            [[0.66666667, 0.11111111, 0.11111111, 0.11111111],
             [0.08333333, 0.75, 0.08333333, 0.08333333],
             [0.08333333, 0.08333333, 0.75, 0.08333333],
             [0.03125, 0.03125, 0.03125, 0.90625]])
        p_initial_out, p_transition_out, p_emission_out = MarkovModel._mle(
            len(states), len(alphabet), training_outputs, training_states, None, None, None)
        self.assertTrue(
            array_equal(around(p_initial_out, decimals=3), around(p_initial, decimals=3)))
        self.assertTrue(
            array_equal(around(p_transition_out, decimals=3), around(p_transition, decimals=3)))
        self.assertTrue(
            array_equal(around(p_emission_out, decimals=3), around(p_emission, decimals=3)))

    def test_argmaxes(self):
        matrix = array([[4, 5, 6], [9, 7, 8], [1, 2, 3]])
        output = [3]
        self.assertEqual(len(MarkovModel._argmaxes(matrix)), len(output))
        self.assertEqual(MarkovModel._argmaxes(matrix)[0], output[0])

    def test_viterbi(self):
        states = ["CP", "IP"]
        outputs = [2, 1, 0]
        lp_initial = log([1.0, 0.0000001])
        lp_transition = log([[0.7, 0.3], [0.5, 0.5]])
        lp_emission = log([[0.6, 0.1, 0.3], [0.1, 0.7, 0.2]])

        output1 = [0, 1, 0]
        output2 = -3.968593356916541

        viterbi_output = MarkovModel._viterbi(
            len(states), lp_initial, lp_transition,
            lp_emission, outputs)
        self.assertEqual(len(viterbi_output[0][0]), 3)
        self.assertEqual(viterbi_output[0][0][0], output1[0])
        self.assertEqual(viterbi_output[0][0][1], output1[1])
        self.assertEqual(viterbi_output[0][0][2], output1[2])
        self.assertEqual(
            float('%.3f' % viterbi_output[0][1]),
            float('%.3f' % output2))

    def test_normalize_and_copy_and_check(self):
        matrix_in1 = array(
            [[1.1, 2.2, 3.3], [4.4, 5.5, 6.6], [7.7, 8.8, 9.9]])
        matrix_in2 = array([1, 2, 3])

        matrix_out1 = array(
            [[0.16666667, 0.33333333, 0.5], [0.26666667, 0.33333333, 0.4], [0.29166667, 0.33333333, 0.375]])
        matrix_out2 = array([0.16666667, 0.33333333, 0.5])
        self.assertTrue(
            array_equal(around(MarkovModel._normalize(matrix_in1), decimals=3), around(matrix_out1, decimals=3)))
        self.assertTrue(
            array_equal(around(MarkovModel._normalize(matrix_in2), decimals=3), around(matrix_out2, decimals=3)))

        shape1 = (3, 3)
        shape2 = (3,)
        self.assertTrue(
            array_equal(around(MarkovModel._copy_and_check(matrix_out1, shape1), decimals=3), around(matrix_out1, decimals=3)))
        self.assertTrue(
            array_equal(around(MarkovModel._copy_and_check(matrix_out2, shape2), decimals=3), around(matrix_out2, decimals=3)))

    def test_uniform_norm(self):
        shape = (4, 3)
        matrix = array([[0.33333333, 0.33333333, 0.33333333],
                        [0.33333333, 0.33333333, 0.33333333],
                        [0.33333333, 0.33333333, 0.33333333],
                        [0.33333333, 0.33333333, 0.33333333]])
        self.assertTrue(
            array_equal(around(MarkovModel._uniform_norm(shape), decimals=3), around(matrix, decimals=3)))

    def test_random_norm(self):
        random.seed(0)
        shape = (4, 3)
        matrix = array([[0.29399155, 0.38311672, 0.32289173],
                        [0.33750765, 0.26241723, 0.40007512],
                        [0.1908342, 0.38890714, 0.42025866],
                        [0.22501625, 0.46461061, 0.31037314]])
        self.assertTrue(
            array_equal(around(MarkovModel._random_norm(shape), decimals=3), around(matrix, decimals=3)))

    def test_logsum_and_exp_logsum(self):
        matrix = array(
            [[1.1, 2.2, 3.3], [4.4, 5.5, 6.6], [7.7, 8.8, 9.9]])
        matrix1 = array([1, 2, 3])

        output = 10.304721798
        output1 = 3.40760596444
        self.assertEqual(
            float('%.3f' % MarkovModel._logsum(matrix)),
            float('%.3f' % output))
        self.assertEqual(
            float('%.3f' % MarkovModel._logsum(matrix1)),
            float('%.3f' % output1))

        output2 = 29873.342245
        output3 = 30.1928748506
        self.assertEqual(
            float('%.3f' % MarkovModel._exp_logsum(matrix)),
            float('%.3f' % output2))
        self.assertEqual(
            float('%.3f' % MarkovModel._exp_logsum(matrix1)),
            float('%.3f' % output3))

    def test_logvecadd(self):
        vec1 = log(array([1, 2, 3, 4]))
        vec2 = log(array([5, 6, 7, 8]))

        sumvec = array([1.79175947, 2.07944154, 2.30258509, 2.48490665])
        self.assertTrue(
            array_equal(around(MarkovModel._logvecadd(vec1, vec2), decimals=3), around(sumvec, decimals=3)))

if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
