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

# HMM stuff we are testing
from Bio.HMM import MarkovModel
from Bio.HMM import Trainer
from Bio.HMM import Utilities


# whether we should print everything out. Set this to zero for
# regression testing
VERBOSE = 0


# -- set up our alphabets
dice_roll_alphabet = ("1", "2", "3", "4", "5", "6")
dice_type_alphabet = ("F", "L")


def generate_rolls(num_rolls):
    """Generate a bunch of rolls corresponding to the casino probabilities.

    Returns:
    - The generate roll sequence
    - The state sequence that generated the roll.

    """
    # start off in the fair state
    cur_state = "F"
    roll_seq = []
    state_seq = []
    loaded_weights = [0.1, 0.1, 0.1, 0.1, 0.1, 0.5]
    # generate the sequence
    for roll in range(num_rolls):
        state_seq.append(cur_state)
        # add on a new roll to the sequence
        if cur_state == "F":
            new_rolls = random.choices(dice_roll_alphabet)
        elif cur_state == "L":
            new_rolls = random.choices(dice_roll_alphabet, weights=loaded_weights)
        new_roll = new_rolls[0]

        roll_seq.append(new_roll)
        # now give us a chance to switch to a new state
        chance_num = random.random()
        if cur_state == "F":
            if chance_num <= 0.05:
                cur_state = "L"
        elif cur_state == "L":
            if chance_num <= 0.1:
                cur_state = "F"
    return roll_seq, state_seq


class TestHMMCasino(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mm_builder = MarkovModel.MarkovModelBuilder(
            dice_type_alphabet, dice_roll_alphabet
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
        predicted_states, prob = trained_mm.viterbi(test_rolls, dice_type_alphabet)
        if VERBOSE:
            print(f"Prediction probability: {prob:f}")
            Utilities.pretty_print_prediction(test_rolls, test_states, predicted_states)

    def test_baum_welch_training_without(self):
        """Baum-Welch training without known state sequences."""
        training_seq = Trainer.TrainingSequence(self.rolls, ())

        def stop_training(log_likelihood_change, num_iterations):
            """Tell the training model when to stop."""
            if VERBOSE:
                print(f"ll change: {log_likelihood_change:f}")
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
        predicted_states, prob = trained_mm.viterbi(test_rolls, dice_type_alphabet)
        if VERBOSE:
            print(f"Prediction probability: {prob:f}")
            Utilities.pretty_print_prediction(
                self.test_rolls, test_states, predicted_states
            )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
