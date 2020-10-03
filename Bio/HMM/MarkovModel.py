# Copyright 2001 Brad Chapman.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Deal with representations of Markov Models."""
# standard modules
import copy
import math
import random
from collections import defaultdict

from Bio.Seq import Seq


def _gen_random_array(n):
    """Return an array of n random numbers summing to 1.0 (PRIVATE)."""
    randArray = [random.random() for _ in range(n)]
    total = sum(randArray)
    return [x / total for x in randArray]


def _calculate_emissions(emission_probs):
    """Calculate which symbols can be emitted in each state (PRIVATE)."""
    # loop over all of the state-symbol duples, mapping states to
    # lists of emitted symbols
    emissions = defaultdict(list)
    for state, symbol in emission_probs:
        emissions[state].append(symbol)

    return emissions


def _calculate_from_transitions(trans_probs):
    """Calculate which 'from transitions' are allowed for each state (PRIVATE).

    This looks through all of the trans_probs, and uses this dictionary
    to determine allowed transitions. It converts this information into
    a dictionary, whose keys are source states and whose values are
    lists of destination states reachable from the source state via a
    transition.
    """
    transitions = defaultdict(list)
    for from_state, to_state in trans_probs:
        transitions[from_state].append(to_state)

    return transitions


def _calculate_to_transitions(trans_probs):
    """Calculate which 'to transitions' are allowed for each state (PRIVATE).

    This looks through all of the trans_probs, and uses this dictionary
    to determine allowed transitions. It converts this information into
    a dictionary, whose keys are destination states and whose values are
    lists of source states from which the destination is reachable via a
    transition.
    """
    transitions = defaultdict(list)
    for from_state, to_state in trans_probs:
        transitions[to_state].append(from_state)

    return transitions


class MarkovModelBuilder:
    """Interface to build up a Markov Model.

    This class is designed to try to separate the task of specifying the
    Markov Model from the actual model itself. This is in hopes of making
    the actual Markov Model classes smaller.

    So, this builder class should be used to create Markov models instead
    of trying to initiate a Markov Model directly.
    """

    # the default pseudo counts to use
    DEFAULT_PSEUDO = 1

    def __init__(self, state_alphabet, emission_alphabet):
        """Initialize a builder to create Markov Models.

        Arguments:
         - state_alphabet -- An iterable (e.g., tuple or list) containing
           all of the letters that can appear in the states
         - emission_alphabet -- An iterable (e.g., tuple or list) containing
           all of the letters for states that can be emitted by the HMM.

        """
        self._state_alphabet = tuple(state_alphabet)
        self._emission_alphabet = tuple(emission_alphabet)

        # probabilities for the initial state, initialized by calling
        # set_initial_probabilities (required)
        self.initial_prob = {}

        # the probabilities for transitions and emissions
        # by default we have no transitions and all possible emissions
        self.transition_prob = {}
        self.emission_prob = self._all_blank(state_alphabet, emission_alphabet)

        # the default pseudocounts for transition and emission counting
        self.transition_pseudo = {}
        self.emission_pseudo = self._all_pseudo(state_alphabet, emission_alphabet)

    def _all_blank(self, first_alphabet, second_alphabet):
        """Return a dictionary with all counts set to zero (PRIVATE).

        This uses the letters in the first and second alphabet to create
        a dictionary with keys of two tuples organized as
        (letter of first alphabet, letter of second alphabet). The values
        are all set to 0.
        """
        all_blank = {}
        for first_state in first_alphabet:
            for second_state in second_alphabet:
                all_blank[(first_state, second_state)] = 0

        return all_blank

    def _all_pseudo(self, first_alphabet, second_alphabet):
        """Return a dictionary with all counts set to a default value (PRIVATE).

        This takes the letters in first alphabet and second alphabet and
        creates a dictionary with keys of two tuples organized as:
        (letter of first alphabet, letter of second alphabet). The values
        are all set to the value of the class attribute DEFAULT_PSEUDO.
        """
        all_counts = {}
        for first_state in first_alphabet:
            for second_state in second_alphabet:
                all_counts[(first_state, second_state)] = self.DEFAULT_PSEUDO

        return all_counts

    def get_markov_model(self):
        """Return the markov model corresponding with the current parameters.

        Each markov model returned by a call to this function is unique
        (ie. they don't influence each other).
        """
        # user must set initial probabilities
        if not self.initial_prob:
            raise Exception(
                "set_initial_probabilities must be called to "
                "fully initialize the Markov model"
            )

        initial_prob = copy.deepcopy(self.initial_prob)
        transition_prob = copy.deepcopy(self.transition_prob)
        emission_prob = copy.deepcopy(self.emission_prob)
        transition_pseudo = copy.deepcopy(self.transition_pseudo)
        emission_pseudo = copy.deepcopy(self.emission_pseudo)

        return HiddenMarkovModel(
            self._state_alphabet,
            self._emission_alphabet,
            initial_prob,
            transition_prob,
            emission_prob,
            transition_pseudo,
            emission_pseudo,
        )

    def set_initial_probabilities(self, initial_prob):
        """Set initial state probabilities.

        initial_prob is a dictionary mapping states to probabilities.
        Suppose, for example, that the state alphabet is ('A', 'B'). Call
        set_initial_prob({'A': 1}) to guarantee that the initial
        state will be 'A'. Call set_initial_prob({'A': 0.5, 'B': 0.5})
        to make each initial state equally probable.

        This method must now be called in order to use the Markov model
        because the calculation of initial probabilities has changed
        incompatibly; the previous calculation was incorrect.

        If initial probabilities are set for all states, then they should add up
        to 1. Otherwise the sum should be <= 1. The residual probability is
        divided up evenly between all the states for which the initial
        probability has not been set. For example, calling
        set_initial_prob({}) results in P('A') = 0.5 and P('B') = 0.5,
        for the above example.
        """
        self.initial_prob = copy.copy(initial_prob)

        # ensure that all referenced states are valid
        for state in initial_prob:
            if state not in self._state_alphabet:
                raise ValueError(
                    "State %s was not found in the sequence alphabet" % state
                )

        # distribute the residual probability, if any
        num_states_not_set = len(self._state_alphabet) - len(self.initial_prob)
        if num_states_not_set < 0:
            raise Exception("Initial probabilities can't exceed # of states")
        prob_sum = sum(self.initial_prob.values())
        if prob_sum > 1.0:
            raise Exception("Total initial probability cannot exceed 1.0")
        if num_states_not_set > 0:
            prob = (1.0 - prob_sum) / num_states_not_set
            for state in self._state_alphabet:
                if state not in self.initial_prob:
                    self.initial_prob[state] = prob

    def set_equal_probabilities(self):
        """Reset all probabilities to be an average value.

        Resets the values of all initial probabilities and all allowed
        transitions and all allowed emissions to be equal to 1 divided by the
        number of possible elements.

        This is useful if you just want to initialize a Markov Model to
        starting values (ie. if you have no prior notions of what the
        probabilities should be -- or if you are just feeling too lazy
        to calculate them :-).

        Warning 1 -- this will reset all currently set probabilities.

        Warning 2 -- This just sets all probabilities for transitions and
        emissions to total up to 1, so it doesn't ensure that the sum of
        each set of transitions adds up to 1.
        """
        # set initial state probabilities
        new_initial_prob = float(1) / float(len(self.transition_prob))
        for state in self._state_alphabet:
            self.initial_prob[state] = new_initial_prob

        # set the transitions
        new_trans_prob = float(1) / float(len(self.transition_prob))
        for key in self.transition_prob:
            self.transition_prob[key] = new_trans_prob

        # set the emissions
        new_emission_prob = float(1) / float(len(self.emission_prob))
        for key in self.emission_prob:
            self.emission_prob[key] = new_emission_prob

    def set_random_initial_probabilities(self):
        """Set all initial state probabilities to a randomly generated distribution.

        Returns the dictionary containing the initial probabilities.
        """
        initial_freqs = _gen_random_array(len(self._state_alphabet))
        for state in self._state_alphabet:
            self.initial_prob[state] = initial_freqs.pop()

        return self.initial_prob

    def set_random_transition_probabilities(self):
        """Set all allowed transition probabilities to a randomly generated distribution.

        Returns the dictionary containing the transition probabilities.
        """
        if not self.transition_prob:
            raise Exception(
                "No transitions have been allowed yet. "
                "Allow some or all transitions by calling "
                "allow_transition or allow_all_transitions first."
            )

        transitions_from = _calculate_from_transitions(self.transition_prob)
        for from_state in transitions_from:
            freqs = _gen_random_array(len(transitions_from[from_state]))
            for to_state in transitions_from[from_state]:
                self.transition_prob[(from_state, to_state)] = freqs.pop()

        return self.transition_prob

    def set_random_emission_probabilities(self):
        """Set all allowed emission probabilities to a randomly generated distribution.

        Returns the dictionary containing the emission probabilities.
        """
        if not self.emission_prob:
            raise Exception(
                "No emissions have been allowed yet. Allow some or all emissions."
            )

        emissions = _calculate_emissions(self.emission_prob)
        for state in emissions:
            freqs = _gen_random_array(len(emissions[state]))
            for symbol in emissions[state]:
                self.emission_prob[(state, symbol)] = freqs.pop()

        return self.emission_prob

    def set_random_probabilities(self):
        """Set all probabilities to randomly generated numbers.

        Resets probabilities of all initial states, transitions, and
        emissions to random values.
        """
        self.set_random_initial_probabilities()
        self.set_random_transition_probabilities()
        self.set_random_emission_probabilities()

    # --- functions to deal with the transitions in the sequence

    def allow_all_transitions(self):
        """Create transitions between all states.

        By default all transitions within the alphabet are disallowed;
        this is a convenience function to change this to allow all
        possible transitions.
        """
        # first get all probabilities and pseudo counts set
        # to the default values
        all_probs = self._all_blank(self._state_alphabet, self._state_alphabet)

        all_pseudo = self._all_pseudo(self._state_alphabet, self._state_alphabet)

        # now set any probabilities and pseudo counts that
        # were previously set
        for set_key in self.transition_prob:
            all_probs[set_key] = self.transition_prob[set_key]

        for set_key in self.transition_pseudo:
            all_pseudo[set_key] = self.transition_pseudo[set_key]

        # finally reinitialize the transition probs and pseudo counts
        self.transition_prob = all_probs
        self.transition_pseudo = all_pseudo

    def allow_transition(
        self, from_state, to_state, probability=None, pseudocount=None
    ):
        """Set a transition as being possible between the two states.

        probability and pseudocount are optional arguments
        specifying the probabilities and pseudo counts for the transition.
        If these are not supplied, then the values are set to the
        default values.

        Raises:
        KeyError -- if the two states already have an allowed transition.

        """
        # check the sanity of adding these states
        for state in [from_state, to_state]:
            if state not in self._state_alphabet:
                raise ValueError(
                    "State %s was not found in the sequence alphabet" % state
                )

        # ensure that the states are not already set
        if (from_state, to_state) not in self.transition_prob and (
            from_state,
            to_state,
        ) not in self.transition_pseudo:
            # set the initial probability
            if probability is None:
                probability = 0
            self.transition_prob[(from_state, to_state)] = probability

            # set the initial pseudocounts
            if pseudocount is None:
                pseudocount = self.DEFAULT_PSEUDO
            self.transition_pseudo[(from_state, to_state)] = pseudocount
        else:
            raise KeyError(
                "Transition from %s to %s is already allowed." % (from_state, to_state)
            )

    def destroy_transition(self, from_state, to_state):
        """Restrict transitions between the two states.

        Raises:
        KeyError if the transition is not currently allowed.

        """
        try:
            del self.transition_prob[(from_state, to_state)]
            del self.transition_pseudo[(from_state, to_state)]
        except KeyError:
            raise KeyError(
                "Transition from %s to %s is already disallowed."
                % (from_state, to_state)
            )

    def set_transition_score(self, from_state, to_state, probability):
        """Set the probability of a transition between two states.

        Raises:
        KeyError if the transition is not allowed.

        """
        if (from_state, to_state) in self.transition_prob:
            self.transition_prob[(from_state, to_state)] = probability
        else:
            raise KeyError(
                "Transition from %s to %s is not allowed." % (from_state, to_state)
            )

    def set_transition_pseudocount(self, from_state, to_state, count):
        """Set the default pseudocount for a transition.

        To avoid computational problems, it is helpful to be able to
        set a 'default' pseudocount to start with for estimating
        transition and emission probabilities (see p62 in Durbin et al
        for more discussion on this. By default, all transitions have
        a pseudocount of 1.

        Raises:
        KeyError if the transition is not allowed.

        """
        if (from_state, to_state) in self.transition_pseudo:
            self.transition_pseudo[(from_state, to_state)] = count
        else:
            raise KeyError(
                "Transition from %s to %s is not allowed." % (from_state, to_state)
            )

    # --- functions to deal with emissions from the sequence

    def set_emission_score(self, seq_state, emission_state, probability):
        """Set the probability of a emission from a particular state.

        Raises:
        KeyError if the emission from the given state is not allowed.

        """
        if (seq_state, emission_state) in self.emission_prob:
            self.emission_prob[(seq_state, emission_state)] = probability
        else:
            raise KeyError(
                "Emission of %s from %s is not allowed." % (emission_state, seq_state)
            )

    def set_emission_pseudocount(self, seq_state, emission_state, count):
        """Set the default pseudocount for an emission.

        To avoid computational problems, it is helpful to be able to
        set a 'default' pseudocount to start with for estimating
        transition and emission probabilities (see p62 in Durbin et al
        for more discussion on this. By default, all emissions have
        a pseudocount of 1.

        Raises:
        KeyError if the emission from the given state is not allowed.

        """
        if (seq_state, emission_state) in self.emission_pseudo:
            self.emission_pseudo[(seq_state, emission_state)] = count
        else:
            raise KeyError(
                "Emission of %s from %s is not allowed." % (emission_state, seq_state)
            )


class HiddenMarkovModel:
    """Represent a hidden markov model that can be used for state estimation."""

    def __init__(
        self,
        state_alphabet,
        emission_alphabet,
        initial_prob,
        transition_prob,
        emission_prob,
        transition_pseudo,
        emission_pseudo,
    ):
        """Initialize a Markov Model.

        Note: You should use the MarkovModelBuilder class instead of
        initiating this class directly.

        Arguments:
         - state_alphabet -- A tuple containing all of the letters that can
           appear in the states.
         - emission_alphabet -- A tuple containing all of the letters for
           states that can be emitted by the HMM.
         - initial_prob - A dictionary of initial probabilities for all states.
         - transition_prob -- A dictionary of transition probabilities for all
           possible transitions in the sequence.
         - emission_prob -- A dictionary of emission probabilities for all
           possible emissions from the sequence states.
         - transition_pseudo -- Pseudo-counts to be used for the transitions,
           when counting for purposes of estimating transition probabilities.
         - emission_pseudo -- Pseudo-counts to be used for the emissions,
           when counting for purposes of estimating emission probabilities.

        """
        self.state_alphabet = state_alphabet
        self.emission_alphabet = emission_alphabet

        self.initial_prob = initial_prob

        self._transition_pseudo = transition_pseudo
        self._emission_pseudo = emission_pseudo

        self.transition_prob = transition_prob
        self.emission_prob = emission_prob

        # a dictionary of the possible transitions from each state
        # each key is a source state, mapped to a list of the destination states
        # that are reachable from the source state via a transition
        self._transitions_from = _calculate_from_transitions(self.transition_prob)

        # a dictionary of the possible transitions to each state
        # each key is a destination state, mapped to a list of source states
        # from which the destination is reachable via a transition
        self._transitions_to = _calculate_to_transitions(self.transition_prob)

    def get_blank_transitions(self):
        """Get the default transitions for the model.

        Returns a dictionary of all of the default transitions between any
        two letters in the sequence alphabet. The dictionary is structured
        with keys as (letter1, letter2) and values as the starting number
        of transitions.
        """
        return self._transition_pseudo

    def get_blank_emissions(self):
        """Get the starting default emmissions for each sequence.

        This returns a dictionary of the default emmissions for each
        letter. The dictionary is structured with keys as
        (seq_letter, emmission_letter) and values as the starting number
        of emmissions.
        """
        return self._emission_pseudo

    def transitions_from(self, state_letter):
        """Get all destination states which can transition from source state_letter.

        This returns all letters which the given state_letter can transition
        to, i.e. all the destination states reachable from state_letter.

        An empty list is returned if state_letter has no outgoing transitions.
        """
        if state_letter in self._transitions_from:
            return self._transitions_from[state_letter]
        else:
            return []

    def transitions_to(self, state_letter):
        """Get all source states which can transition to destination state_letter.

        This returns all letters which the given state_letter is reachable
        from, i.e. all the source states which can reach state_later

        An empty list is returned if state_letter is unreachable.
        """
        if state_letter in self._transitions_to:
            return self._transitions_to[state_letter]
        else:
            return []

    def viterbi(self, sequence, state_alphabet):
        """Calculate the most probable state path using the Viterbi algorithm.

        This implements the Viterbi algorithm (see pgs 55-57 in Durbin et
        al for a full explanation -- this is where I took my implementation
        ideas from), to allow decoding of the state path, given a sequence
        of emissions.

        Arguments:
         - sequence -- A Seq object with the emission sequence that we
           want to decode.
         - state_alphabet -- An iterable (e.g., tuple or list) containing
           all of the letters that can appear in the states

        """
        # calculate logarithms of the initial, transition, and emission probs
        log_initial = self._log_transform(self.initial_prob)
        log_trans = self._log_transform(self.transition_prob)
        log_emission = self._log_transform(self.emission_prob)

        viterbi_probs = {}
        pred_state_seq = {}

        # --- recursion
        # loop over the training squence (i = 1 .. L)
        # NOTE: My index numbers are one less than what is given in Durbin
        # et al, since we are indexing the sequence going from 0 to
        # (Length - 1) not 1 to Length, like in Durbin et al.
        for i in range(0, len(sequence)):
            # loop over all of the possible i-th states in the state path
            for cur_state in state_alphabet:
                # e_{l}(x_{i})
                emission_part = log_emission[(cur_state, sequence[i])]

                max_prob = 0
                if i == 0:
                    # for the first state, use the initial probability rather
                    # than looking back to previous states
                    max_prob = log_initial[cur_state]
                else:
                    # loop over all possible (i-1)-th previous states
                    possible_state_probs = {}
                    for prev_state in self.transitions_to(cur_state):
                        # a_{kl}
                        trans_part = log_trans[(prev_state, cur_state)]

                        # v_{k}(i - 1)
                        viterbi_part = viterbi_probs[(prev_state, i - 1)]
                        cur_prob = viterbi_part + trans_part

                        possible_state_probs[prev_state] = cur_prob

                    # calculate the viterbi probability using the max
                    max_prob = max(possible_state_probs.values())

                # v_{k}(i)
                viterbi_probs[(cur_state, i)] = emission_part + max_prob

                if i > 0:
                    # get the most likely prev_state leading to cur_state
                    for state in possible_state_probs:
                        if possible_state_probs[state] == max_prob:
                            pred_state_seq[(i - 1, cur_state)] = state
                            break

        # --- termination
        # calculate the probability of the state path
        # loop over all states
        all_probs = {}
        for state in state_alphabet:
            # v_{k}(L)
            all_probs[state] = viterbi_probs[(state, len(sequence) - 1)]

        state_path_prob = max(all_probs.values())

        # find the last pointer we need to trace back from
        last_state = ""
        for state in all_probs:
            if all_probs[state] == state_path_prob:
                last_state = state

        assert last_state != "", "Didn't find the last state to trace from!"

        # --- traceback
        traceback_seq = []

        loop_seq = list(range(1, len(sequence)))
        loop_seq.reverse()

        # last_state is the last state in the most probable state sequence.
        # Compute that sequence by walking backwards in time. From the i-th
        # state in the sequence, find the (i-1)-th state as the most
        # probable state preceding the i-th state.
        state = last_state
        traceback_seq.append(state)
        for i in loop_seq:
            state = pred_state_seq[(i - 1, state)]
            traceback_seq.append(state)

        # put the traceback sequence in the proper orientation
        traceback_seq.reverse()
        traceback_seq = "".join(traceback_seq)

        return Seq(traceback_seq), state_path_prob

    def _log_transform(self, probability):
        """Return log transform of the given probability dictionary (PRIVATE).

        When calculating the Viterbi equation, add logs of probabilities rather
        than multiplying probabilities, to avoid underflow errors. This method
        returns a new dictionary with the same keys as the given dictionary
        and log-transformed values.
        """
        log_prob = copy.copy(probability)
        for key in log_prob:
            prob = log_prob[key]
            if prob > 0:
                log_prob[key] = math.log(log_prob[key])
            else:
                log_prob[key] = -math.inf

        return log_prob
