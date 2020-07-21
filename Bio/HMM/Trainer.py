# Copyright 2001 Brad Chapman.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Provide trainers which estimate parameters based on training sequences.

These should be used to 'train' a Markov Model prior to actually using
it to decode state paths. When supplied training sequences and a model
to work from, these classes will estimate parameters of the model.

This aims to estimate two parameters:

- a_{kl} -- the number of times there is a transition from k to l in the
  training data.
- e_{k}(b) -- the number of emissions of the state b from the letter k
  in the training data.

"""
# standard modules
import math

# local stuff
from .DynamicProgramming import ScaledDPAlgorithms


class TrainingSequence:
    """Hold a training sequence with emissions and optionally, a state path."""

    def __init__(self, emissions, state_path):
        """Initialize a training sequence.

        Arguments:
         - emissions - An iterable (e.g., a tuple, list, or Seq object)
           containing the sequence of emissions in the training sequence.
         - state_path - An iterable (e.g., a tuple or list) containing the
           sequence of states. If there is no known state path, then the
           sequence of states should be an empty iterable.

        """
        if len(state_path) > 0 and len(emissions) != len(state_path):
            raise ValueError("State path does not match associated emissions.")
        self.emissions = emissions
        self.states = state_path


class AbstractTrainer:
    """Provide generic functionality needed in all trainers."""

    def __init__(self, markov_model):
        """Initialize."""
        self._markov_model = markov_model

    def log_likelihood(self, probabilities):
        """Calculate the log likelihood of the training seqs.

        Arguments:
         - probabilities -- A list of the probabilities of each training
           sequence under the current parameters, calculated using the
           forward algorithm.

        """
        total_likelihood = 0
        for probability in probabilities:
            total_likelihood += math.log(probability)

        return total_likelihood

    def estimate_params(self, transition_counts, emission_counts):
        """Get a maximum likelihood estimation of transition and emmission.

        Arguments:
         - transition_counts -- A dictionary with the total number of counts
           of transitions between two states.
         - emissions_counts -- A dictionary with the total number of counts
           of emmissions of a particular emission letter by a state letter.

        This then returns the maximum likelihood estimators for the
        transitions and emissions, estimated by formulas 3.18 in
        Durbin et al::

            a_{kl} = A_{kl} / sum(A_{kl'})
            e_{k}(b) = E_{k}(b) / sum(E_{k}(b'))

        Returns:
        Transition and emission dictionaries containing the maximum
        likelihood estimators.

        """
        # now calculate the information
        ml_transitions = self.ml_estimator(transition_counts)
        ml_emissions = self.ml_estimator(emission_counts)

        return ml_transitions, ml_emissions

    def ml_estimator(self, counts):
        """Calculate the maximum likelihood estimator.

        This can calculate maximum likelihoods for both transitions
        and emissions.

        Arguments:
         - counts -- A dictionary of the counts for each item.

        See estimate_params for a description of the formula used for
        calculation.

        """
        # get an ordered list of all items
        all_ordered = sorted(counts)

        ml_estimation = {}

        # the total counts for the current letter we are on
        cur_letter = None
        cur_letter_counts = 0

        for cur_item in all_ordered:
            # if we are on a new letter (ie. the first letter of the tuple)
            if cur_item[0] != cur_letter:
                # set the new letter we are working with
                cur_letter = cur_item[0]

                # count up the total counts for this letter
                cur_letter_counts = counts[cur_item]

                # add counts for all other items with the same first letter
                cur_position = all_ordered.index(cur_item) + 1

                # keep adding while we have the same first letter or until
                # we get to the end of the ordered list
                while (
                    cur_position < len(all_ordered)
                    and all_ordered[cur_position][0] == cur_item[0]
                ):
                    cur_letter_counts += counts[all_ordered[cur_position]]
                    cur_position += 1
            # otherwise we've already got the total counts for this letter
            else:
                pass

            # now calculate the ml and add it to the estimation
            cur_ml = float(counts[cur_item]) / float(cur_letter_counts)
            ml_estimation[cur_item] = cur_ml

        return ml_estimation


class BaumWelchTrainer(AbstractTrainer):
    """Trainer that uses the Baum-Welch algorithm to estimate parameters.

    These should be used when a training sequence for an HMM has unknown
    paths for the actual states, and you need to make an estimation of the
    model parameters from the observed emissions.

    This uses the Baum-Welch algorithm, first described in
    Baum, L.E. 1972. Inequalities. 3:1-8
    This is based on the description in 'Biological Sequence Analysis' by
    Durbin et al. in section 3.3

    This algorithm is guaranteed to converge to a local maximum, but not
    necessarily to the global maxima, so use with care!
    """

    def __init__(self, markov_model):
        """Initialize the trainer.

        Arguments:
         - markov_model - The model we are going to estimate parameters for.
           This should have the parameters with some initial estimates, that
           we can build from.

        """
        AbstractTrainer.__init__(self, markov_model)

    def train(self, training_seqs, stopping_criteria, dp_method=ScaledDPAlgorithms):
        """Estimate the parameters using training sequences.

        The algorithm for this is taken from Durbin et al. p64, so this
        is a good place to go for a reference on what is going on.

        Arguments:
         - training_seqs -- A list of TrainingSequence objects to be used
           for estimating the parameters.
         - stopping_criteria -- A function, that when passed the change
           in log likelihood and threshold, will indicate if we should stop
           the estimation iterations.
         - dp_method -- A class instance specifying the dynamic programming
           implementation we should use to calculate the forward and
           backward variables. By default, we use the scaling method.

        """
        prev_log_likelihood = None
        num_iterations = 1

        while True:
            transition_count = self._markov_model.get_blank_transitions()
            emission_count = self._markov_model.get_blank_emissions()

            # remember all of the sequence probabilities
            all_probabilities = []

            for training_seq in training_seqs:
                # calculate the forward and backward variables
                DP = dp_method(self._markov_model, training_seq)
                forward_var, seq_prob = DP.forward_algorithm()
                backward_var = DP.backward_algorithm()

                all_probabilities.append(seq_prob)

                # update the counts for transitions and emissions
                transition_count = self.update_transitions(
                    transition_count, training_seq, forward_var, backward_var, seq_prob
                )
                emission_count = self.update_emissions(
                    emission_count, training_seq, forward_var, backward_var, seq_prob
                )

            # update the markov model with the new probabilities
            ml_transitions, ml_emissions = self.estimate_params(
                transition_count, emission_count
            )
            self._markov_model.transition_prob = ml_transitions
            self._markov_model.emission_prob = ml_emissions

            cur_log_likelihood = self.log_likelihood(all_probabilities)

            # if we have previously calculated the log likelihood (ie.
            # not the first round), see if we can finish
            if prev_log_likelihood is not None:
                # XXX log likelihoods are negatives -- am I calculating
                # the change properly, or should I use the negatives...
                # I'm not sure at all if this is right.
                log_likelihood_change = abs(
                    abs(cur_log_likelihood) - abs(prev_log_likelihood)
                )

                # check whether we have completed enough iterations to have
                # a good estimation
                if stopping_criteria(log_likelihood_change, num_iterations):
                    break

            # set up for another round of iterations
            prev_log_likelihood = cur_log_likelihood
            num_iterations += 1

        return self._markov_model

    def update_transitions(
        self,
        transition_counts,
        training_seq,
        forward_vars,
        backward_vars,
        training_seq_prob,
    ):
        """Add the contribution of a new training sequence to the transitions.

        Arguments:
         - transition_counts -- A dictionary of the current counts for the
           transitions
         - training_seq -- The training sequence we are working with
         - forward_vars -- Probabilities calculated using the forward
           algorithm.
         - backward_vars -- Probabilities calculated using the backwards
           algorithm.
         - training_seq_prob - The probability of the current sequence.

        This calculates A_{kl} (the estimated transition counts from state
        k to state l) using formula 3.20 in Durbin et al.

        """
        # set up the transition and emission probabilities we are using
        transitions = self._markov_model.transition_prob
        emissions = self._markov_model.emission_prob

        # loop over the possible combinations of state path letters
        for k in self._markov_model.state_alphabet:
            for l in self._markov_model.transitions_from(k):
                estimated_counts = 0
                # now loop over the entire training sequence
                for i in range(len(training_seq.emissions) - 1):
                    # the forward value of k at the current position
                    forward_value = forward_vars[(k, i)]

                    # the backward value of l in the next position
                    backward_value = backward_vars[(l, i + 1)]

                    # the probability of a transition from k to l
                    trans_value = transitions[(k, l)]

                    # the probability of getting the emission at the next pos
                    emm_value = emissions[(l, training_seq.emissions[i + 1])]

                    estimated_counts += (
                        forward_value * trans_value * emm_value * backward_value
                    )

                # update the transition approximation
                transition_counts[(k, l)] += float(estimated_counts) / training_seq_prob

        return transition_counts

    def update_emissions(
        self,
        emission_counts,
        training_seq,
        forward_vars,
        backward_vars,
        training_seq_prob,
    ):
        """Add the contribution of a new training sequence to the emissions.

        Arguments:
         - emission_counts -- A dictionary of the current counts for the
           emissions
         - training_seq -- The training sequence we are working with
         - forward_vars -- Probabilities calculated using the forward
           algorithm.
         - backward_vars -- Probabilities calculated using the backwards
           algorithm.
         - training_seq_prob - The probability of the current sequence.

        This calculates E_{k}(b) (the estimated emission probability for
        emission letter b from state k) using formula 3.21 in Durbin et al.

        """
        # loop over the possible combinations of state path letters
        for k in self._markov_model.state_alphabet:
            # now loop over all of the possible emissions
            for b in self._markov_model.emission_alphabet:
                expected_times = 0
                # finally loop over the entire training sequence
                for i in range(len(training_seq.emissions)):
                    # only count the forward and backward probability if the
                    # emission at the position is the same as b
                    if training_seq.emissions[i] == b:
                        # f_{k}(i) b_{k}(i)
                        expected_times += forward_vars[(k, i)] * backward_vars[(k, i)]

                # add to E_{k}(b)
                emission_counts[(k, b)] += float(expected_times) / training_seq_prob

        return emission_counts


class KnownStateTrainer(AbstractTrainer):
    """Estimate probabilities with known state sequences.

    This should be used for direct estimation of emission and transition
    probabilities when both the state path and emission sequence are
    known for the training examples.
    """

    def __init__(self, markov_model):
        """Initialize."""
        AbstractTrainer.__init__(self, markov_model)

    def train(self, training_seqs):
        """Estimate the Markov Model parameters with known state paths.

        This trainer requires that both the state and the emissions are
        known for all of the training sequences in the list of
        TrainingSequence objects.
        This training will then count all of the transitions and emissions,
        and use this to estimate the parameters of the model.
        """
        # count up all of the transitions and emissions
        transition_counts = self._markov_model.get_blank_transitions()
        emission_counts = self._markov_model.get_blank_emissions()

        for training_seq in training_seqs:
            emission_counts = self._count_emissions(training_seq, emission_counts)
            transition_counts = self._count_transitions(
                training_seq.states, transition_counts
            )

        # update the markov model from the counts
        ml_transitions, ml_emissions = self.estimate_params(
            transition_counts, emission_counts
        )
        self._markov_model.transition_prob = ml_transitions
        self._markov_model.emission_prob = ml_emissions

        return self._markov_model

    def _count_emissions(self, training_seq, emission_counts):
        """Add emissions from the training sequence to the current counts (PRIVATE).

        Arguments:
         - training_seq -- A TrainingSequence with states and emissions
           to get the counts from
         - emission_counts -- The current emission counts to add to.

        """
        for index in range(len(training_seq.emissions)):
            cur_state = training_seq.states[index]
            cur_emission = training_seq.emissions[index]

            try:
                emission_counts[(cur_state, cur_emission)] += 1
            except KeyError:
                raise KeyError(
                    "Unexpected emission (%s, %s)" % (cur_state, cur_emission)
                )
        return emission_counts

    def _count_transitions(self, state_seq, transition_counts):
        """Add transitions from the training sequence to the current counts (PRIVATE).

        Arguments:
         - state_seq -- A Seq object with the states of the current training
           sequence.
         - transition_counts -- The current transition counts to add to.

        """
        for cur_pos in range(len(state_seq) - 1):
            cur_state = state_seq[cur_pos]
            next_state = state_seq[cur_pos + 1]

            try:
                transition_counts[(cur_state, next_state)] += 1
            except KeyError:
                raise KeyError(
                    "Unexpected transition (%s, %s)" % (cur_state, next_state)
                )

        return transition_counts
