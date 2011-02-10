#!/usr/bin/env python
"""Test the HMM.MarkovModel and HMM.DynamicProgramming modules.

Also tests Training methods.
"""
# standard modules
import unittest
import math

# biopython
from Bio import Alphabet
from Bio.Seq import Seq


# stuff we are testing
from Bio.HMM import MarkovModel
from Bio.HMM import DynamicProgramming
from Bio.HMM import Trainer

# create some simple alphabets
class NumberAlphabet(Alphabet.Alphabet):
    """Numbers as the states of the model.
    """
    letters = ['1', '2']

class LetterAlphabet(Alphabet.Alphabet):
    """Letters as the emissions of the model.
    """
    letters = ['A', 'B']

# -- helper functions
def test_assertion(name, result, expected):
    """Helper function to test an assertion and print out a reasonable error.
    """
    assert result == expected, "Expected %s, got %s for %s" \
           % (expected, result, name)
    
class MarkovModelBuilderTest(unittest.TestCase):
    def setUp(self):
        self.mm_builder = MarkovModel.MarkovModelBuilder(NumberAlphabet(),
                                                         LetterAlphabet())

    def test_test_initialize(self):
        """Making sure MarkovModelBuilder is initialized correctly.
        """
        expected_transition_prob = {}
        expected_transition_pseudo = {}

        expected_emission_prob = {('2', 'A'): 0, ('1', 'A'): 0,
                                  ('1', 'B'): 0, ('2', 'B'): 0}
        expected_emission_pseudo = {('2', 'A'): 1, ('1', 'A'): 1,
                                    ('1', 'B'): 1, ('2', 'B'): 1}

        assertions = []
        test_assertion("Transition prob", self.mm_builder.transition_prob,
                          expected_transition_prob)
        test_assertion("Transition pseudo",
                          self.mm_builder.transition_pseudo,
                          expected_transition_pseudo)
        test_assertion("Emission prob", self.mm_builder.emission_prob,
                           expected_emission_prob)
        test_assertion("Emission pseudo", self.mm_builder.emission_pseudo,
                           expected_emission_pseudo)


    def test_allow_all_transitions(self):
        """Testing allow_all_transitions.
        """
        self.mm_builder.allow_all_transitions()

        expected_prob = {('2', '1'): 0, ('1', '1'): 0,
                         ('1', '2'): 0, ('2', '2'): 0}

        expected_pseudo = {('2', '1'): 1, ('1', '1'): 1,
                           ('1', '2'): 1, ('2', '2'): 1}

        test_assertion("Probabilities", self.mm_builder.transition_prob,
                       expected_prob)
        
        test_assertion("Pseudo counts",  self.mm_builder.transition_pseudo,
                       expected_pseudo)

class HiddenMarkovModelTest(unittest.TestCase):
    def setUp(self):
        self.mm_builder = MarkovModel.MarkovModelBuilder(NumberAlphabet(),
                                                    LetterAlphabet())

    def test_transitions_from(self):
        """Testing the calculation of transitions_from
        """
        self.mm_builder.allow_transition('1', '2', 1.0)
        self.mm_builder.allow_transition('2', '1', 0.5)
        self.mm_builder.allow_transition('2', '2', 0.5)
        self.mm = self.mm_builder.get_markov_model()

        state_1 = self.mm.transitions_from("1")
        expected_state_1 = ["2"]
        state_1.sort()
        expected_state_1.sort()
        test_assertion("States reached by transitions from state 1",
                       state_1, expected_state_1)

        state_2 = self.mm.transitions_from("2")
        expected_state_2 = ["1", "2"]
        state_2.sort()
        expected_state_2.sort()
        test_assertion("States reached by transitions from state 2",
                       state_2, expected_state_2)

        fake_state = self.mm.transitions_from("Fake")
        expected_fake_state = []
        test_assertion("States reached by transitions from a fake transition",
                       fake_state, expected_fake_state)

    def test_transitions_to(self):
        """Testing the calculation of transitions_to
        """
        self.mm_builder.allow_transition('1', '1', 0.5)
        self.mm_builder.allow_transition('1', '2', 0.5)
        self.mm_builder.allow_transition('2', '1', 1.0)
        self.mm = self.mm_builder.get_markov_model()

        state_1 = self.mm.transitions_to("1")
        expected_state_1 = ["1", "2"]
        state_1.sort()
        expected_state_1.sort()
        test_assertion("States with transitions to state 1",
                       state_1, expected_state_1)

        state_2 = self.mm.transitions_to("2")
        expected_state_2 = ["1"]
        state_2.sort()
        expected_state_2.sort()
        test_assertion("States with transitions to state 2",
                       state_2, expected_state_2)

        fake_state = self.mm.transitions_to("Fake")
        expected_fake_state = []
        test_assertion("States with transitions to a fake transition",
                       fake_state, expected_fake_state)

    def test_allow_transition(self):
        """Testing allow_transition
        """
        self.mm_builder.allow_transition('1', '2', 1.0)
        self.mm = self.mm_builder.get_markov_model()

        state_1 = self.mm.transitions_from("1")
        expected_state_1 = ["2"]
        state_1.sort()
        expected_state_1.sort()
        test_assertion("States reached by transitions from state 1",
                       state_1, expected_state_1)

        state_2 = self.mm.transitions_from("2")
        expected_state_2 = []
        state_2.sort()
        expected_state_2.sort()
        test_assertion("States reached by transitions from state 2",
                       state_2, expected_state_2)

        state_1 = self.mm.transitions_to("1")
        expected_state_1 = []
        state_1.sort()
        expected_state_1.sort()
        test_assertion("States with transitions to state 1",
                       state_1, expected_state_1)

        state_2 = self.mm.transitions_to("2")
        expected_state_2 = ["1"]
        state_2.sort()
        expected_state_2.sort()
        test_assertion("States with transitions to state 2",
                       state_2, expected_state_2)


    def test_non_ergodic(self):
        """Test a non-ergodic model (meaning that some transitions are not
        allowed).
        """

        # probabilities of transitioning from state 1 to 1, and 1 to 2
        prob_1_to_1 = 0.5
        prob_1_to_2 = 0.5

        # set up allowed transitions
        self.mm_builder.allow_transition('1', '1', prob_1_to_1)
        self.mm_builder.allow_transition('1', '2', prob_1_to_2)

        # Emission probabilities
        # In state 1 the most likely emission is A, in state 2 the most
        # likely emission is B. (Would be simpler just to use 1.0 and 0.0
        # emission probabilities here, but the algorithm blows up on zero
        # probabilities because of the conversion to log space.)
        prob_emit_A_in_state_1 = 0.95
        prob_emit_B_in_state_1 = 0.05
        prob_emit_A_in_state_2 = 0.05
        prob_emit_B_in_state_2 = 0.95

        # set emission probabilities
        self.mm_builder.set_emission_score('1', 'A', prob_emit_A_in_state_1)
        self.mm_builder.set_emission_score('1', 'B', prob_emit_B_in_state_1)
        self.mm_builder.set_emission_score('2', 'A', prob_emit_A_in_state_2)
        self.mm_builder.set_emission_score('2', 'B', prob_emit_B_in_state_2)

        # run the Viterbi algorithm to find the most probable state path
        model = self.mm_builder.get_markov_model()
        observed_emissions = ['A', 'B']
        viterbi=model.viterbi(observed_emissions, NumberAlphabet)
        seq=viterbi[0]
        prob=viterbi[1]

        # the most probable path must be from state 1 to state 2
        test_assertion("most probable path", str(seq), '12')

        # The probability of that path is the probability of transitioning
        # from state 1 to state 1 (1 -> 1), then emitting an A, then
        # transitioning 1 -> 2, then emitting a B. That first hidden
        # transition (not part of the published state sequence) from
        # 1 -> 1 is not a real transition, rather it's an implementation
        # detail that determines the initial state of the sequence.
        # Note that probabilities are converted into log space.
        expected_prob = math.log(prob_1_to_1)\
        + math.log(prob_emit_A_in_state_1)\
        + math.log(prob_1_to_2)\
        + math.log(prob_emit_B_in_state_2)
        test_assertion("log probability of most probable path",
                       prob, expected_prob)

class ScaledDPAlgorithmsTest(unittest.TestCase):
    def setUp(self):
        # set up our Markov Model
        mm_builder = MarkovModel.MarkovModelBuilder(NumberAlphabet(),
                                                    LetterAlphabet())
        mm_builder.allow_all_transitions()
        mm_builder.set_equal_probabilities()

        mm = mm_builder.get_markov_model()

        # now set up a test sequence
        emission_seq = Seq("ABB", LetterAlphabet())
        state_seq = Seq("", NumberAlphabet())
        training_seq = Trainer.TrainingSequence(emission_seq, state_seq)

        # finally set up the DP
        self.dp = DynamicProgramming.ScaledDPAlgorithms(mm, training_seq)
        
    def test_calculate_s_value(self):
        """Testing the calculation of s values.
        """
        previous_vars = {('1', 0) : .5,
                         ('2', 0) : .7}
        s_value = self.dp._calculate_s_value(1, previous_vars)

        # print s_value

class AbstractTrainerTest(unittest.TestCase):
    def setUp(self):
        # set up a bogus HMM and our trainer
        hmm = MarkovModel.HiddenMarkovModel({}, {}, {}, {})
        self.test_trainer = Trainer.AbstractTrainer(hmm)
    
    def test_ml_estimator(self):
        """Test the maximum likelihood estimator for simple cases.
        """
        # set up a simple dictionary
        counts = {('A', 'A') : 10,
                  ('A', 'B') : 20,
                  ('A', 'C') : 15,
                  ('B', 'B') : 5,
                  ('C', 'A') : 15,
                  ('C', 'C') : 10}

        results = self.test_trainer.ml_estimator(counts)

        # now make sure we are getting back the right thing
        result_tests = []
        result_tests.append([('A', 'A'), float(10) / float(45)])
        result_tests.append([('A', 'B'), float(20) / float(45)])
        result_tests.append([('A', 'C'), float(15) / float(45)])
        result_tests.append([('B', 'B'), float(5) / float(5)])
        result_tests.append([('C', 'A'), float(15) / float(25)])
        result_tests.append([('C', 'C'), float(10) / float(25)])

        for test_result in result_tests:
            assert results[test_result[0]] == test_result[1], \
                   "Got %f, expected %f for %s" % (results[test_result[0]],
                                                   test_result[1],
                                                   test_result[0])

    def test_log_likelihood(self):
        """Calculate log likelihood.
        """
        probs = [.25, .13, .12, .17]

        log_prob = self.test_trainer.log_likelihood(probs)
        expected_log_prob = -7.31873556778
        assert abs(expected_log_prob - log_prob) < 0.1, \
          "Bad probability calculated: %s" % log_prob

# run the tests
if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
