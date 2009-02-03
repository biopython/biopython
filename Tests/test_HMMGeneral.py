#!/usr/bin/env python
"""Test the HMM.MarkovModel and HMM.DynamicProgramming modules.

Also tests Training methods.
"""
# standard modules
import unittest

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
        """Testing allow_all_transtions.
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
        mm_builder = MarkovModel.MarkovModelBuilder(NumberAlphabet(),
                                                    LetterAlphabet())
        mm_builder.allow_all_transitions()

        self.mm = mm_builder.get_markov_model()

    def test_transitions_from(self):
        """Testing the calculation of transitions_from
        """
        state_1 = self.mm.transitions_from("1")
        expected_state_1 = ["1", "2"]
        state_1.sort()
        expected_state_1.sort()
        test_assertion("State 1 transitions", state_1, expected_state_1)

        state_2 = self.mm.transitions_from("2")
        expected_state_2 = ["1", "2"]
        state_2.sort()
        expected_state_2.sort()
        test_assertion("State 2 transitions", state_2, expected_state_2)

        fake_state = self.mm.transitions_from("Fake")
        expected_fake_state = []
        test_assertion("Fake transition", fake_state, expected_fake_state)

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
