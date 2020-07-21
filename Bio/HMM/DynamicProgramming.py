# Copyright 2001 Brad Chapman.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Dynamic Programming algorithms for general usage.

This module contains classes which implement Dynamic Programming
algorithms that can be used generally.
"""


class AbstractDPAlgorithms:
    """An abstract class to calculate forward and backward probabilities.

    This class should not be instantiated directly, but should be used
    through a derived class which implements proper scaling of variables.

    This class is just meant to encapsulate the basic forward and backward
    algorithms, and allow derived classes to deal with the problems of
    multiplying probabilities.

    Derived class of this must implement:

    - _forward_recursion -- Calculate the forward values in the recursion
      using some kind of technique for preventing underflow errors.
    - _backward_recursion -- Calculate the backward values in the recursion
      step using some technique to prevent underflow errors.

    """

    def __init__(self, markov_model, sequence):
        """Initialize to calculate forward and backward probabilities.

        Arguments:
         - markov_model -- The current Markov model we are working with.
         - sequence -- A training sequence containing a set of emissions.

        """
        self._mm = markov_model
        self._seq = sequence

    def _forward_recursion(self, cur_state, sequence_pos, forward_vars):
        """Calculate the forward recursion value (PRIVATE)."""
        raise NotImplementedError("Subclasses must implement")

    def forward_algorithm(self):
        """Calculate sequence probability using the forward algorithm.

        This implements the forward algorithm, as described on p57-58 of
        Durbin et al.

        Returns:
         - A dictionary containing the forward variables. This has keys of the
           form (state letter, position in the training sequence), and values
           containing the calculated forward variable.
         - The calculated probability of the sequence.

        """
        # all of the different letters that the state path can be in
        state_letters = self._mm.state_alphabet

        # -- initialize the algorithm
        #
        # NOTE: My index numbers are one less than what is given in Durbin
        # et al, since we are indexing the sequence going from 0 to
        # (Length - 1) not 1 to Length, like in Durbin et al.
        #
        forward_var = {}
        # f_{0}(0) = 1
        forward_var[(state_letters[0], -1)] = 1
        # f_{k}(0) = 0, for k > 0
        for k in range(1, len(state_letters)):
            forward_var[(state_letters[k], -1)] = 0

        # -- now do the recursion step
        # loop over the training sequence
        # Recursion step: (i = 1 .. L)
        for i in range(len(self._seq.emissions)):
            # now loop over the letters in the state path
            for main_state in state_letters:
                # calculate the forward value using the appropriate
                # method to prevent underflow errors
                forward_value = self._forward_recursion(main_state, i, forward_var)

                if forward_value is not None:
                    forward_var[(main_state, i)] = forward_value

        # -- termination step - calculate the probability of the sequence
        first_state = state_letters[0]
        seq_prob = 0

        for state_item in state_letters:
            # f_{k}(L)
            forward_value = forward_var[(state_item, len(self._seq.emissions) - 1)]
            # a_{k0}
            transition_value = self._mm.transition_prob[(state_item, first_state)]

            seq_prob += forward_value * transition_value

        return forward_var, seq_prob

    def _backward_recursion(self, cur_state, sequence_pos, forward_vars):
        """Calculate the backward recursion value (PRIVATE)."""
        raise NotImplementedError("Subclasses must implement")

    def backward_algorithm(self):
        """Calculate sequence probability using the backward algorithm.

        This implements the backward algorithm, as described on p58-59 of
        Durbin et al.

        Returns:
         - A dictionary containing the backwards variables. This has keys
           of the form (state letter, position in the training sequence),
           and values containing the calculated backward variable.

        """
        # all of the different letters that the state path can be in
        state_letters = self._mm.state_alphabet

        # -- initialize the algorithm
        #
        # NOTE: My index numbers are one less than what is given in Durbin
        # et al, since we are indexing the sequence going from 0 to
        # (Length - 1) not 1 to Length, like in Durbin et al.
        #
        backward_var = {}

        first_letter = state_letters[0]
        # b_{k}(L) = a_{k0} for all k
        for state in state_letters:
            backward_var[
                (state, len(self._seq.emissions) - 1)
            ] = self._mm.transition_prob[(state, state_letters[0])]

        # -- recursion
        # first loop over the training sequence backwards
        # Recursion step: (i = L - 1 ... 1)
        all_indexes = list(range(len(self._seq.emissions) - 1))
        all_indexes.reverse()
        for i in all_indexes:
            # now loop over the letters in the state path
            for main_state in state_letters:
                # calculate the backward value using the appropriate
                # method to prevent underflow errors
                backward_value = self._backward_recursion(main_state, i, backward_var)

                if backward_value is not None:
                    backward_var[(main_state, i)] = backward_value

        # skip the termination step to avoid recalculations -- you should
        # get sequence probabilities using the forward algorithm

        return backward_var


class ScaledDPAlgorithms(AbstractDPAlgorithms):
    """Implement forward and backward algorithms using a rescaling approach.

    This scales the f and b variables, so that they remain within a
    manageable numerical interval during calculations. This approach is
    described in Durbin et al. on p 78.

    This approach is a little more straightforward then log transformation
    but may still give underflow errors for some types of models. In these
    cases, the LogDPAlgorithms class should be used.
    """

    def __init__(self, markov_model, sequence):
        """Initialize the scaled approach to calculating probabilities.

        Arguments:
         - markov_model -- The current Markov model we are working with.
         - sequence -- A TrainingSequence object that must have a
           set of emissions to work with.

        """
        AbstractDPAlgorithms.__init__(self, markov_model, sequence)

        self._s_values = {}

    def _calculate_s_value(self, seq_pos, previous_vars):
        """Calculate the next scaling variable for a sequence position (PRIVATE).

        This utilizes the approach of choosing s values such that the
        sum of all of the scaled f values is equal to 1.

        Arguments:
         - seq_pos -- The current position we are at in the sequence.
         - previous_vars -- All of the forward or backward variables
           calculated so far.

        Returns:
         - The calculated scaling variable for the sequence item.

        """
        # all of the different letters the state can have
        state_letters = self._mm.state_alphabet

        # loop over all of the possible states
        s_value = 0
        for main_state in state_letters:
            emission = self._mm.emission_prob[
                (main_state, self._seq.emissions[seq_pos])
            ]

            # now sum over all of the previous vars and transitions
            trans_and_var_sum = 0
            for second_state in self._mm.transitions_from(main_state):
                # the value of the previous f or b value
                var_value = previous_vars[(second_state, seq_pos - 1)]

                # the transition probability
                trans_value = self._mm.transition_prob[(second_state, main_state)]

                trans_and_var_sum += var_value * trans_value

            s_value += emission * trans_and_var_sum

        return s_value

    def _forward_recursion(self, cur_state, sequence_pos, forward_vars):
        """Calculate the value of the forward recursion (PRIVATE).

        Arguments:
         - cur_state -- The letter of the state we are calculating the
           forward variable for.
         - sequence_pos -- The position we are at in the training seq.
         - forward_vars -- The current set of forward variables

        """
        # calculate the s value, if we haven't done so already (ie. during
        # a previous forward or backward recursion)
        if sequence_pos not in self._s_values:
            self._s_values[sequence_pos] = self._calculate_s_value(
                sequence_pos, forward_vars
            )

        # e_{l}(x_{i})
        seq_letter = self._seq.emissions[sequence_pos]
        cur_emission_prob = self._mm.emission_prob[(cur_state, seq_letter)]
        # divide by the scaling value
        scale_emission_prob = float(cur_emission_prob) / float(
            self._s_values[sequence_pos]
        )

        # loop over all of the possible states at the position
        state_pos_sum = 0
        have_transition = 0
        for second_state in self._mm.transitions_from(cur_state):
            have_transition = 1

            # get the previous forward_var values
            # f_{k}(i - 1)
            prev_forward = forward_vars[(second_state, sequence_pos - 1)]

            # a_{kl}
            cur_trans_prob = self._mm.transition_prob[(second_state, cur_state)]
            state_pos_sum += prev_forward * cur_trans_prob

        # if we have the possibility of having a transition
        # return the recursion value
        if have_transition:
            return scale_emission_prob * state_pos_sum
        else:
            return None

    def _backward_recursion(self, cur_state, sequence_pos, backward_vars):
        """Calculate the value of the backward recursion (PRIVATE).

        Arguments:
         - cur_state -- The letter of the state we are calculating the
           forward variable for.
         - sequence_pos -- The position we are at in the training seq.
         - backward_vars -- The current set of backward variables

        """
        # calculate the s value, if we haven't done so already (ie. during
        # a previous forward or backward recursion)
        if sequence_pos not in self._s_values:
            self._s_values[sequence_pos] = self._calculate_s_value(
                sequence_pos, backward_vars
            )

        # loop over all of the possible states at the position
        state_pos_sum = 0
        have_transition = 0
        for second_state in self._mm.transitions_from(cur_state):
            have_transition = 1
            # e_{l}(x_{i + 1})
            seq_letter = self._seq.emissions[sequence_pos + 1]
            cur_emission_prob = self._mm.emission_prob[(cur_state, seq_letter)]

            # get the previous backward_var value
            # b_{l}(i + 1)
            prev_backward = backward_vars[(second_state, sequence_pos + 1)]

            # the transition probability -- a_{kl}
            cur_transition_prob = self._mm.transition_prob[(cur_state, second_state)]

            state_pos_sum += cur_emission_prob * prev_backward * cur_transition_prob

        # if we have a probability for a transition, return it
        if have_transition:
            return state_pos_sum / float(self._s_values[sequence_pos])
        # otherwise we have no probability (ie. we can't do this transition)
        # and return None
        else:
            return None


class LogDPAlgorithms(AbstractDPAlgorithms):
    """Implement forward and backward algorithms using a log approach.

    This uses the approach of calculating the sum of log probabilities
    using a lookup table for common values.

    XXX This is not implemented yet!
    """

    def __init__(self, markov_model, sequence):
        """Initialize."""
        raise NotImplementedError("Haven't coded this yet...")
