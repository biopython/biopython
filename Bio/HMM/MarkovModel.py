"""Deal with representations of Markov Models.
"""
# standard modules
import copy
import math
import random

# biopython
from Bio.Seq import MutableSeq

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

        o state_alphabet -- An alphabet containing all of the letters that
        can appear in the states
       
        o emission_alphabet -- An alphabet containing all of the letters for
        states that can be emitted by the HMM.
        """
        self._state_alphabet = state_alphabet
        self._emission_alphabet = emission_alphabet
        
        # the probabilities for transitions and emissions
        # by default we have no transitions and all possible emissions
        self.transition_prob = {}
        self.emission_prob = self._all_blank(state_alphabet,
                                             emission_alphabet)

        # the default pseudocounts for transition and emission counting
        self.transition_pseudo = {}
        self.emission_pseudo = self._all_pseudo(state_alphabet,
                                                emission_alphabet)

    def _all_blank(self, first_alphabet, second_alphabet):
        """Return a dictionary with all counts set to zero.

        This uses the letters in the first and second alphabet to create
        a dictionary with keys of two tuples organized as
        (letter of first alphabet, letter of second alphabet). The values
        are all set to 0.
        """
        all_blank = {}
        for first_state in first_alphabet.letters:
            for second_state in second_alphabet.letters:
                all_blank[(first_state, second_state)] = 0

        return all_blank

    def _all_pseudo(self, first_alphabet, second_alphabet):
        """Return a dictionary with all counts set to a default value.

        This takes the letters in first alphabet and second alphabet and
        creates a dictionary with keys of two tuples organized as:
        (letter of first alphabet, letter of second alphabet). The values
        are all set to the value of the class attribute DEFAULT_PSEUDO.
        """
        all_counts = {}
        for first_state in first_alphabet.letters:
            for second_state in second_alphabet.letters:
                all_counts[(first_state, second_state)] = self.DEFAULT_PSEUDO

        return all_counts
                
    def get_markov_model(self):
        """Return the markov model corresponding with the current parameters.

        Each markov model returned by a call to this function is unique
        (ie. they don't influence each other).
        """
        transition_prob = copy.deepcopy(self.transition_prob)
        emission_prob = copy.deepcopy(self.emission_prob)
        transition_pseudo = copy.deepcopy(self.transition_pseudo)
        emission_pseudo = copy.deepcopy(self.emission_pseudo)
        
        return HiddenMarkovModel(transition_prob, emission_prob,
                                 transition_pseudo, emission_pseudo)

    def set_equal_probabilities(self):
        """Reset all probabilities to be an average value.

        This resets the values of all allowed transitions and all allowed
        emissions to be equal to 1 divided by the number of possible elements.

        This is useful if you just want to initialize a Markov Model to
        starting values (ie. if you have no prior notions of what the
        probabilities should be -- or if you are just feeling too lazy
        to calculate them :-).

        Warning 1 -- this will reset all currently set probabilities.

        Warning 2 -- This just sets all probabilities for transitions and
        emissions to total up to 1, so it doesn't ensure that the sum of
        each set of transitions adds up to 1.
        """
        # first set the transitions
        new_trans_prob = float(1) / float(len(self.transition_prob))
        for key in self.transition_prob:
            self.transition_prob[key] = new_trans_prob

        # now set the emissions
        new_emission_prob = float(1) / float(len(self.emission_prob))
        for key in self.emission_prob:
            self.emission_prob[key] = new_emission_prob
            

    def set_random_probabilities(self):
        """Set all probabilities to randomly generated numbers.

        This will reset the value of all allowed transitions and emissions
        to random values.

        Warning 1 -- This will reset any currently set probabibilities.

        Warning 2 -- This does not check to ensure that the sum of
        all of the probabilities is less then 1. It just randomly assigns
        a probability to each
        """
        for key in self.transition_prob:
            self.transition_prob[key] = random.random()

        for key in self.emission_prob:
            self.emission_prob[key] = random.random()

    # --- functions to deal with the transitions in the sequence

    def allow_all_transitions(self):
        """A convenience function to create transitions between all states.

        By default all transitions within the alphabet are disallowed; this
        is a way to change this to allow all possible transitions.
        """
        # first get all probabilities and pseudo counts set
        # to the default values
        all_probs = self._all_blank(self._state_alphabet,
                                    self._state_alphabet)

        all_pseudo = self._all_pseudo(self._state_alphabet,
                                      self._state_alphabet)

        # now set any probabilities and pseudo counts that
        # were previously set
        for set_key in self.transition_prob:
            all_probs[set_key] = self.transition_prob[set_key]

        for set_key in self.transition_pseudo:
            all_pseudo[set_key] = self.transition_pseudo[set_key]

        # finally reinitialize the transition probs and pseudo counts
        self.transition_prob = all_probs
        self.transition_pseudo = all_pseudo

    def allow_transition(self, from_state, to_state, probability = None,
                         pseudocount = None):
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
            assert state in self._state_alphabet.letters, \
                   "State %s was not found in the sequence alphabet" % state

        # ensure that the states are not already set
        if ((from_state, to_state) not in self.transition_prob and 
            (from_state, to_state) not in self.transition_pseudo):
            # set the initial probability
            if probability is None:
                probability = 0
            self.transition_prob[(from_state, to_state)] = probability

            # set the initial pseudocounts
            if pseudocount is None:
                pseudcount = self.DEFAULT_PSEUDO
            self.transition_pseudo[(from_state, to_state)] = pseudocount 
        else:
            raise KeyError("Transtion from %s to %s is already allowed."
                           % (from_state, to_state))

    def destroy_transition(self, from_state, to_state):
        """Restrict transitions between the two states.

        Raises:
        KeyError if the transition is not currently allowed.
        """
        try:
            del self.transition_prob[(from_state, to_state)]
            del self.transition_pseudo[(from_state, to_state)]
        except KeyError:
            raise KeyError("Transition from %s to %s is already disallowed."
                           % (from_state, to_state))

    def set_transition_score(self, from_state, to_state, probability):
        """Set the probability of a transition between two states.

        Raises:
        KeyError if the transition is not allowed.
        """
        if (from_state, to_state) in self.transition_prob:
            self.transition_prob[(from_state, to_state)] = probability
        else:
            raise KeyError("Transition from %s to %s is not allowed."
                           % (from_state, to_state))

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
            raise KeyError("Transition from %s to %s is not allowed."
                           % (from_state, to_state))

    # --- functions to deal with emissions from the sequence

    def set_emission_score(self, seq_state, emission_state, probability):
        """Set the probability of a emission from a particular state.

        Raises:
        KeyError if the emission from the given state is not allowed.
        """
        if (seq_state, emission_state) in self.emission_prob:
            self.emission_prob[(seq_state, emission_state)] = probability
        else:
            raise KeyError("Emission of %s from %s is not allowed."
                           % (emission_state, seq_state))

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
            raise KeyError("Emission of %s from %s is not allowed."
                           % (emission_state, seq_state))

class HiddenMarkovModel:
    """Represent a hidden markov model that can be used for state estimation.
    """
    def __init__(self, transition_prob, emission_prob, transition_pseudo,
                 emission_pseudo):
        """Initialize a Markov Model.

        Note: You should use the MarkovModelBuilder class instead of
        initiating this class directly.

        Arguments:

        o transition_prob -- A dictionary of transition probabilities for all
        possible transitions in the sequence.

        o emission_prob -- A dictionary of emission probabilities for all
        possible emissions from the sequence states.

        o transition_pseudo -- Pseudo-counts to be used for the transitions,
        when counting for purposes of estimating transition probabilities.

        o emission_pseudo -- Pseudo-counts to be used for the emissions,
        when counting for purposes of estimating emission probabilities.
        """
        self._transition_pseudo = transition_pseudo
        self._emission_pseudo = emission_pseudo
        
        self.transition_prob = transition_prob
        self.emission_prob = emission_prob

        # a dictionary of the possible transitions from each state
        # each key is a source state, mapped to a list of the destination states
        # that are reachable from the source state via a transition
        self._transitions_from = \
           self._calculate_from_transitions(self.transition_prob)

        # a dictionary of the possible transitions to each state
        # each key is a destination state, mapped to a list of source states
        # from which the destination is reachable via a transition
        self._transitions_to = \
           self._calculate_to_transitions(self.transition_prob)

    def _calculate_from_transitions(self, trans_probs):
        """Calculate which 'from transitions' are allowed for each state

        This looks through all of the trans_probs, and uses this dictionary
        to determine allowed transitions. It converts this information into
        a dictionary, whose keys are source states and whose values are
        lists of destination states reachable from the source state via a
        transition.
        """
        from_transitions = {}

        # loop over all of the transitions, mapping source states to lists
        # of destination states
        for trans_key in trans_probs:
            from_state = trans_key[0]
            to_state = trans_key[1]
            # Add to_state to the list of destination states that are reachable
            # from from_state. Map from_state if we haven't seen it before.
            if from_state in from_transitions:
                from_transitions[from_state].append(to_state)
            # otherwise create the list and add the state
            else:
                from_transitions[from_state] = [to_state]

        return from_transitions

    def _calculate_to_transitions(self, trans_probs):
        """Calculate which 'to transitions' are allowed for each state

        This looks through all of the trans_probs, and uses this dictionary
        to determine allowed transitions. It converts this information into
        a dictionary, whose keys are destination states and whose values are
        lists of source states from which the destination is reachable via a
        transition.
        """
        to_transitions = {}

        # loop over all of the transitions, mapping destination states to lists
        # of source states
        for trans_key in trans_probs:
            from_state = trans_key[0]
            to_state = trans_key[1]
            # Add from_state to the list of source states from which to_state
            # is reachable. Map to_state if we haven't seen it before.
            if to_state in to_transitions:
                to_transitions[to_state].append(from_state)
            # otherwise create the list and add the state
            else:
                to_transitions[to_state] = [from_state]

        return to_transitions

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
        """Get all destination states to which there are transitions from the
        state_letter source state.

        This returns all letters which the given state_letter can transition
        to. An empty list is returned if state_letter has no outgoing
        transitions.
        """
        if state_letter in self._transitions_from:
            return self._transitions_from[state_letter]
        else:
            return []

    def transitions_to(self, state_letter):
        """Get all source states from which there are transitions to the
        state_letter destination state.

        This returns all letters which the given state_letter is reachable
        from. An empty list is returned if state_letter is unreachable.
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

        o sequence -- A Seq object with the emission sequence that we
        want to decode.

        o state_alphabet -- The alphabet of the possible state sequences
        that can be generated.
        """
        # calculate logarithms of the transition and emission probs
        log_trans = self._log_transform(self.transition_prob)
        log_emission = self._log_transform(self.emission_prob)

        viterbi_probs = {}
        pred_state_seq = {}
        state_letters = state_alphabet.letters
        # --- initialization
        #
        # NOTE: My index numbers are one less than what is given in Durbin
        # et al, since we are indexing the sequence going from 0 to
        # (Length - 1) not 1 to Length, like in Durbin et al.
        #
        # v_{0}(0) = 0
        viterbi_probs[(state_letters[0], -1)] = 0
        # v_{k}(0) = 0 for k > 0
        for state_letter in state_letters[1:]:
            viterbi_probs[(state_letter, -1)] = 0

        # --- recursion
        # loop over the training squence (i = 1 .. L)
        for i in range(0, len(sequence)):
            # loop over all of the possible i-th states in the state path
            for cur_state in state_letters:
                # e_{l}(x_{i})
                emission_part = log_emission[(cur_state, sequence[i])]

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
                viterbi_probs[(cur_state, i)] = (emission_part + max_prob)

                # get the most likely prev_state leading to cur_state
                for state in possible_state_probs:
                    if possible_state_probs[state] == max_prob:
                        pred_state_seq[(i - 1, cur_state)] = state
                        break
                    
        # --- termination
        # calculate the probability of the state path
        # loop over all states
        all_probs = {}
        for state in state_letters:
            # v_{k}(L)
            all_probs[state] = viterbi_probs[(state, len(sequence) - 1)]

        state_path_prob = max(all_probs.values())

        # find the last pointer we need to trace back from
        last_state = ''
        for state in all_probs:
            if all_probs[state] == state_path_prob:
                last_state = state

        assert last_state != '', "Didn't find the last state to trace from!"
                
        # --- traceback
        traceback_seq = MutableSeq('', state_alphabet)
        
        loop_seq = range(0, len(sequence))
        loop_seq.reverse()

        # last_state is the last state in the most probable state sequence.
        # Compute that sequence by walking backwards in time. From the i-th
        # state in the sequence, find the (i-1)-th state as the most
        # probable state preceding the i-th state.
        state = last_state
        for i in loop_seq:
            traceback_seq.append(state)            
            state = pred_state_seq[(i - 1, state)]

        # put the traceback sequence in the proper orientation
        traceback_seq.reverse()

        return traceback_seq.toseq(), state_path_prob

    def _log_transform(self, probability):
        """Return log transform of the given probability dictionary.

        When calculating the Viterbi equation, we need to deal with things
        as sums of logs instead of products of probabilities, so that we
        don't get underflow errors.. This copies the given probability
        dictionary and returns the same dictionary with everything
        transformed with a log.
        """
        log_prob = copy.copy(probability)

        for key in log_prob:
            log_prob[key] = math.log(log_prob[key])

        return log_prob
    
