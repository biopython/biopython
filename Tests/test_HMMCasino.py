# Copyright 2001 Brad Chapman.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Test out HMMs using the Occasionally Dishonest Casino.

This uses the occasionally dishonest casino example from Biological
Sequence Analysis by Durbin et al.

In this example, we are dealing with a casino that has two types of
dice, a fair dice that has 1/6 probability of rolling any number and
a loaded dice that has 1/2 probability to roll a 6, and 1/10 probability
to roll any other number. The probability of switching from the fair to
loaded dice is .05 and the probability of switching from loaded to fair is
.1.
"""


# standard modules
import os
import random
import unittest

# biopython
from Bio import Alphabet
from Bio.Seq import MutableSeq
from Bio.Seq import Seq

# HMM stuff we are testing
from Bio.HMM import MarkovModel
from Bio.HMM import Trainer
from Bio.HMM import Utilities


if os.name == "java":
    from Bio import MissingExternalDependencyError

    # This is a slight miss-use of MissingExternalDependencyError,
    # but it will do in the short term to skip this unit test on Jython
    raise MissingExternalDependencyError(
        "This test can cause a fatal error on Jython with some versions of Java"
    )


# whether we should print everything out. Set this to zero for
# regression testing
VERBOSE = 0


# -- set up our alphabets
class DiceRollAlphabet(Alphabet.Alphabet):
    letters = ["1", "2", "3", "4", "5", "6"]


class DiceTypeAlphabet(Alphabet.Alphabet):
    letters = ["F", "L"]


# -- useful functions
def _loaded_dice_roll(chance_num, cur_state):
    """Generate a loaded dice roll based on the state and a random number."""
    if cur_state == "F":
        if chance_num <= (float(1) / float(6)):
            return "1"
        elif chance_num <= (float(2) / float(6)):
            return "2"
        elif chance_num <= (float(3) / float(6)):
            return "3"
        elif chance_num <= (float(4) / float(6)):
            return "4"
        elif chance_num <= (float(5) / float(6)):
            return "5"
        else:
            return "6"
    elif cur_state == "L":
        if chance_num <= (float(1) / float(10)):
            return "1"
        elif chance_num <= (float(2) / float(10)):
            return "2"
        elif chance_num <= (float(3) / float(10)):
            return "3"
        elif chance_num <= (float(4) / float(10)):
            return "4"
        elif chance_num <= (float(5) / float(10)):
            return "5"
        else:
            return "6"
    else:
        raise ValueError("Unexpected cur_state %s" % cur_state)


def generate_rolls(num_rolls):
    """Generate a bunch of rolls corresponding to the casino probabilities.

    Returns:
    - The generate roll sequence
    - The state sequence that generated the roll.

    """
    # start off in the fair state
    cur_state = "F"
    roll_seq = MutableSeq("", DiceRollAlphabet())
    state_seq = MutableSeq("", DiceTypeAlphabet())
    # generate the sequence
    for roll in range(num_rolls):
        state_seq.append(cur_state)
        # generate a random number
        chance_num = random.random()
        # add on a new roll to the sequence
        new_roll = _loaded_dice_roll(chance_num, cur_state)
        roll_seq.append(new_roll)
        # now give us a chance to switch to a new state
        chance_num = random.random()
        if cur_state == "F":
            if chance_num <= 0.05:
                cur_state = "L"
        elif cur_state == "L":
            if chance_num <= 0.1:
                cur_state = "F"
    return roll_seq.toseq(), state_seq.toseq()


class TestHMMCasino(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mm_builder = MarkovModel.MarkovModelBuilder(
            DiceTypeAlphabet(), DiceRollAlphabet()
        )
        cls.mm_builder.allow_all_transitions()
        cls.mm_builder.set_random_probabilities()
        # get a sequence of rolls to train the markov model with
        cls.rolls, cls.states = generate_rolls(3000)

    def test_baum_welch_training_standard(self):
        """Standard Training with known states."""
        known_training_seq = Trainer.TrainingSequence(self.rolls, self.states)
        standard_mm = self.mm_builder.get_markov_model()
        trainer = Trainer.KnownStateTrainer(standard_mm)
        trained_mm = trainer.train([known_training_seq])
        if VERBOSE:
            print(trained_mm.transition_prob)
            print(trained_mm.emission_prob)
        test_rolls, test_states = generate_rolls(300)
        predicted_states, prob = trained_mm.viterbi(test_rolls, DiceTypeAlphabet())
        if VERBOSE:
            print("Prediction probability: %f" % prob)
            Utilities.pretty_print_prediction(test_rolls, test_states, predicted_states)

    def test_baum_welch_training_without(self):
        """Baum-Welch training without known state sequences."""
        training_seq = Trainer.TrainingSequence(self.rolls, Seq("", DiceTypeAlphabet()))

        def stop_training(log_likelihood_change, num_iterations):
            """Tell the training model when to stop."""
            if VERBOSE:
                print("ll change: %f" % log_likelihood_change)
            if log_likelihood_change < 0.01:
                return 1
            elif num_iterations >= 10:
                return 1
            else:
                return 0

        baum_welch_mm = self.mm_builder.get_markov_model()
        trainer = Trainer.BaumWelchTrainer(baum_welch_mm)
        trained_mm = trainer.train([training_seq], stop_training)
        if VERBOSE:
            print(trained_mm.transition_prob)
            print(trained_mm.emission_prob)
        test_rolls, test_states = generate_rolls(300)
        predicted_states, prob = trained_mm.viterbi(test_rolls, DiceTypeAlphabet())
        if VERBOSE:
            print("Prediction probability: %f" % prob)
            Utilities.pretty_print_prediction(
                self.test_rolls, test_states, predicted_states
            )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
