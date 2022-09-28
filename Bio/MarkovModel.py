# Copyright 2002 by Jeffrey Chang.
# All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""A state-emitting MarkovModel.

Note terminology similar to Manning and Schutze is used.


Functions:
train_bw        Train a markov model using the Baum-Welch algorithm.
train_visible   Train a visible markov model using MLE.
find_states     Find the a state sequence that explains some observations.

load            Load a MarkovModel.
save            Save a MarkovModel.

Classes:
MarkovModel     Holds the description of a markov model
"""

import numpy


try:
    logaddexp = numpy.logaddexp
except AttributeError:
    # Numpy versions older than 1.3 do not contain logaddexp.
    # Once we require Numpy version 1.3 or later, we should revisit this
    # module to see if we can simplify some of the other functions in
    # this module.
    import warnings

    warnings.warn(
        "For optimal speed, please update to Numpy version 1.3 or later (current version is %s)"
        % numpy.__version__
    )

    def logaddexp(logx, logy):
        """Implement logaddexp method if Numpy version is older than 1.3."""
        if logy - logx > 100:
            return logy
        elif logx - logy > 100:
            return logx
        minxy = min(logx, logy)
        return minxy + numpy.log(numpy.exp(logx - minxy) + numpy.exp(logy - minxy))


def itemindex(values):
    """Return a dictionary of values with their sequence offset as keys."""
    d = {}
    entries = enumerate(values[::-1])
    n = len(values) - 1
    for index, key in entries:
        d[key] = n - index
    return d


numpy.random.seed()

VERY_SMALL_NUMBER = 1e-300
LOG0 = numpy.log(VERY_SMALL_NUMBER)


class MarkovModel:
    """Create a state-emitting MarkovModel object."""

    def __init__(
        self, states, alphabet, p_initial=None, p_transition=None, p_emission=None
    ):
        """Initialize the class."""
        self.states = states
        self.alphabet = alphabet
        self.p_initial = p_initial
        self.p_transition = p_transition
        self.p_emission = p_emission

    def __str__(self):
        """Create a string representation of the MarkovModel object."""
        from io import StringIO

        handle = StringIO()
        save(self, handle)
        handle.seek(0)
        return handle.read()


def _readline_and_check_start(handle, start):
    """Read the first line and evaluate that begisn with the correct start (PRIVATE)."""
    line = handle.readline()
    if not line.startswith(start):
        raise ValueError(f"I expected {start!r} but got {line!r}")
    return line


def load(handle):
    """Parse a file handle into a MarkovModel object."""
    # Load the states.
    line = _readline_and_check_start(handle, "STATES:")
    states = line.split()[1:]

    # Load the alphabet.
    line = _readline_and_check_start(handle, "ALPHABET:")
    alphabet = line.split()[1:]

    mm = MarkovModel(states, alphabet)
    N, M = len(states), len(alphabet)

    # Load the initial probabilities.
    mm.p_initial = numpy.zeros(N)
    line = _readline_and_check_start(handle, "INITIAL:")
    for i in range(len(states)):
        line = _readline_and_check_start(handle, f"  {states[i]}:")
        mm.p_initial[i] = float(line.split()[-1])

    # Load the transition.
    mm.p_transition = numpy.zeros((N, N))
    line = _readline_and_check_start(handle, "TRANSITION:")
    for i in range(len(states)):
        line = _readline_and_check_start(handle, f"  {states[i]}:")
        mm.p_transition[i, :] = [float(v) for v in line.split()[1:]]

    # Load the emission.
    mm.p_emission = numpy.zeros((N, M))
    line = _readline_and_check_start(handle, "EMISSION:")
    for i in range(len(states)):
        line = _readline_and_check_start(handle, f"  {states[i]}:")
        mm.p_emission[i, :] = [float(v) for v in line.split()[1:]]

    return mm


def save(mm, handle):
    """Save MarkovModel object into handle."""
    # This will fail if there are spaces in the states or alphabet.
    w = handle.write
    w(f"STATES: {' '.join(mm.states)}\n")
    w(f"ALPHABET: {' '.join(mm.alphabet)}\n")
    w("INITIAL:\n")
    for i in range(len(mm.p_initial)):
        w(f"  {mm.states[i]}: {mm.p_initial[i]:g}\n")
    w("TRANSITION:\n")
    for i in range(len(mm.p_transition)):
        w(f"  {mm.states[i]}: {' '.join(str(x) for x in mm.p_transition[i])}\n")
    w("EMISSION:\n")
    for i in range(len(mm.p_emission)):
        w(f"  {mm.states[i]}: {' '.join(str(x) for x in mm.p_emission[i])}\n")


# XXX allow them to specify starting points
def train_bw(
    states,
    alphabet,
    training_data,
    pseudo_initial=None,
    pseudo_transition=None,
    pseudo_emission=None,
    update_fn=None,
):
    """Train a MarkovModel using the Baum-Welch algorithm.

    Train a MarkovModel using the Baum-Welch algorithm.  states is a list
    of strings that describe the names of each state.  alphabet is a
    list of objects that indicate the allowed outputs.  training_data
    is a list of observations.  Each observation is a list of objects
    from the alphabet.

    pseudo_initial, pseudo_transition, and pseudo_emission are
    optional parameters that you can use to assign pseudo-counts to
    different matrices.  They should be matrices of the appropriate
    size that contain numbers to add to each parameter matrix, before
    normalization.

    update_fn is an optional callback that takes parameters
    (iteration, log_likelihood).  It is called once per iteration.
    """
    N, M = len(states), len(alphabet)
    if not training_data:
        raise ValueError("No training data given.")
    if pseudo_initial is not None:
        pseudo_initial = numpy.asarray(pseudo_initial)
        if pseudo_initial.shape != (N,):
            raise ValueError("pseudo_initial not shape len(states)")
    if pseudo_transition is not None:
        pseudo_transition = numpy.asarray(pseudo_transition)
        if pseudo_transition.shape != (N, N):
            raise ValueError("pseudo_transition not shape len(states) X len(states)")
    if pseudo_emission is not None:
        pseudo_emission = numpy.asarray(pseudo_emission)
        if pseudo_emission.shape != (N, M):
            raise ValueError("pseudo_emission not shape len(states) X len(alphabet)")

    # Training data is given as a list of members of the alphabet.
    # Replace those with indexes into the alphabet list for easier
    # computation.
    training_outputs = []
    indexes = itemindex(alphabet)
    for outputs in training_data:
        training_outputs.append([indexes[x] for x in outputs])

    # Do some sanity checking on the outputs.
    lengths = [len(x) for x in training_outputs]
    if min(lengths) == 0:
        raise ValueError("I got training data with outputs of length 0")

    # Do the training with baum welch.
    x = _baum_welch(
        N,
        M,
        training_outputs,
        pseudo_initial=pseudo_initial,
        pseudo_transition=pseudo_transition,
        pseudo_emission=pseudo_emission,
        update_fn=update_fn,
    )
    p_initial, p_transition, p_emission = x
    return MarkovModel(states, alphabet, p_initial, p_transition, p_emission)


MAX_ITERATIONS = 1000


def _baum_welch(
    N,
    M,
    training_outputs,
    p_initial=None,
    p_transition=None,
    p_emission=None,
    pseudo_initial=None,
    pseudo_transition=None,
    pseudo_emission=None,
    update_fn=None,
):
    """Implement the Baum-Welch algorithm to evaluate unknown parameters in the MarkovModel object (PRIVATE)."""
    if p_initial is None:
        p_initial = _random_norm(N)
    else:
        p_initial = _copy_and_check(p_initial, (N,))

    if p_transition is None:
        p_transition = _random_norm((N, N))
    else:
        p_transition = _copy_and_check(p_transition, (N, N))
    if p_emission is None:
        p_emission = _random_norm((N, M))
    else:
        p_emission = _copy_and_check(p_emission, (N, M))

    # Do all the calculations in log space to avoid underflows.
    lp_initial = numpy.log(p_initial)
    lp_transition = numpy.log(p_transition)
    lp_emission = numpy.log(p_emission)
    if pseudo_initial is not None:
        lpseudo_initial = numpy.log(pseudo_initial)
    else:
        lpseudo_initial = None
    if pseudo_transition is not None:
        lpseudo_transition = numpy.log(pseudo_transition)
    else:
        lpseudo_transition = None
    if pseudo_emission is not None:
        lpseudo_emission = numpy.log(pseudo_emission)
    else:
        lpseudo_emission = None

    # Iterate through each sequence of output, updating the parameters
    # to the HMM.  Stop when the log likelihoods of the sequences
    # stops varying.
    prev_llik = None
    for i in range(MAX_ITERATIONS):
        llik = LOG0
        for outputs in training_outputs:
            llik += _baum_welch_one(
                N,
                M,
                outputs,
                lp_initial,
                lp_transition,
                lp_emission,
                lpseudo_initial,
                lpseudo_transition,
                lpseudo_emission,
            )
        if update_fn is not None:
            update_fn(i, llik)
        if prev_llik is not None and numpy.fabs(prev_llik - llik) < 0.1:
            break
        prev_llik = llik
    else:
        raise RuntimeError("HMM did not converge in %d iterations" % MAX_ITERATIONS)

    # Return everything back in normal space.
    return [numpy.exp(_) for _ in (lp_initial, lp_transition, lp_emission)]


def _baum_welch_one(
    N,
    M,
    outputs,
    lp_initial,
    lp_transition,
    lp_emission,
    lpseudo_initial,
    lpseudo_transition,
    lpseudo_emission,
):
    """Execute one step for Baum-Welch algorithm (PRIVATE).

    Do one iteration of Baum-Welch based on a sequence of output.
    Changes the value for lp_initial, lp_transition and lp_emission in place.
    """
    T = len(outputs)
    fmat = _forward(N, T, lp_initial, lp_transition, lp_emission, outputs)
    bmat = _backward(N, T, lp_transition, lp_emission, outputs)

    # Calculate the probability of traversing each arc for any given
    # transition.
    lp_arc = numpy.zeros((N, N, T))
    for t in range(T):
        k = outputs[t]
        lp_traverse = numpy.zeros((N, N))  # P going over one arc.
        for i in range(N):
            for j in range(N):
                # P(getting to this arc)
                # P(making this transition)
                # P(emitting this character)
                # P(going to the end)
                lp = (
                    fmat[i][t]
                    + lp_transition[i][j]
                    + lp_emission[i][k]
                    + bmat[j][t + 1]
                )
                lp_traverse[i][j] = lp
        # Normalize the probability for this time step.
        lp_arc[:, :, t] = lp_traverse - _logsum(lp_traverse)

    # Sum of all the transitions out of state i at time t.
    lp_arcout_t = numpy.zeros((N, T))
    for t in range(T):
        for i in range(N):
            lp_arcout_t[i][t] = _logsum(lp_arc[i, :, t])

    # Sum of all the transitions out of state i.
    lp_arcout = numpy.zeros(N)
    for i in range(N):
        lp_arcout[i] = _logsum(lp_arcout_t[i, :])

    # UPDATE P_INITIAL.
    lp_initial = lp_arcout_t[:, 0]
    if lpseudo_initial is not None:
        lp_initial = _logvecadd(lp_initial, lpseudo_initial)
        lp_initial = lp_initial - _logsum(lp_initial)

    # UPDATE P_TRANSITION.  p_transition[i][j] is the sum of all the
    # transitions from i to j, normalized by the sum of the
    # transitions out of i.
    for i in range(N):
        for j in range(N):
            lp_transition[i][j] = _logsum(lp_arc[i, j, :]) - lp_arcout[i]
        if lpseudo_transition is not None:
            lp_transition[i] = _logvecadd(lp_transition[i], lpseudo_transition)
            lp_transition[i] = lp_transition[i] - _logsum(lp_transition[i])

    # UPDATE P_EMISSION.  lp_emission[i][k] is the sum of all the
    # transitions out of i when k is observed, divided by the sum of
    # the transitions out of i.
    for i in range(N):
        ksum = numpy.zeros(M) + LOG0  # ksum[k] is the sum of all i with k.
        for t in range(T):
            k = outputs[t]
            for j in range(N):
                ksum[k] = logaddexp(ksum[k], lp_arc[i, j, t])
        ksum = ksum - _logsum(ksum)  # Normalize
        if lpseudo_emission is not None:
            ksum = _logvecadd(ksum, lpseudo_emission[i])
            ksum = ksum - _logsum(ksum)  # Renormalize
        lp_emission[i, :] = ksum

    # Calculate the log likelihood of the output based on the forward
    # matrix.  Since the parameters of the HMM has changed, the log
    # likelihoods are going to be a step behind, and we might be doing
    # one extra iteration of training.  The alternative is to rerun
    # the _forward algorithm and calculate from the clean one, but
    # that may be more expensive than overshooting the training by one
    # step.
    return _logsum(fmat[:, T])


def _forward(N, T, lp_initial, lp_transition, lp_emission, outputs):
    """Implement forward algorithm (PRIVATE).

    Calculate a Nx(T+1) matrix, where the last column is the total
    probability of the output.
    """
    matrix = numpy.zeros((N, T + 1))

    # Initialize the first column to be the initial values.
    matrix[:, 0] = lp_initial
    for t in range(1, T + 1):
        k = outputs[t - 1]
        for j in range(N):
            # The probability of the state is the sum of the
            # transitions from all the states from time t-1.
            lprob = LOG0
            for i in range(N):
                lp = matrix[i][t - 1] + lp_transition[i][j] + lp_emission[i][k]
                lprob = logaddexp(lprob, lp)
            matrix[j][t] = lprob
    return matrix


def _backward(N, T, lp_transition, lp_emission, outputs):
    """Implement backward algorithm (PRIVATE)."""
    matrix = numpy.zeros((N, T + 1))
    for t in range(T - 1, -1, -1):
        k = outputs[t]
        for i in range(N):
            # The probability of the state is the sum of the
            # transitions from all the states from time t+1.
            lprob = LOG0
            for j in range(N):
                lp = matrix[j][t + 1] + lp_transition[i][j] + lp_emission[i][k]
                lprob = logaddexp(lprob, lp)
            matrix[i][t] = lprob
    return matrix


def train_visible(
    states,
    alphabet,
    training_data,
    pseudo_initial=None,
    pseudo_transition=None,
    pseudo_emission=None,
):
    """Train a visible MarkovModel using maximum likelihoood estimates for each of the parameters.

    Train a visible MarkovModel using maximum likelihoood estimates
    for each of the parameters.  states is a list of strings that
    describe the names of each state.  alphabet is a list of objects
    that indicate the allowed outputs.  training_data is a list of
    (outputs, observed states) where outputs is a list of the emission
    from the alphabet, and observed states is a list of states from
    states.

    pseudo_initial, pseudo_transition, and pseudo_emission are
    optional parameters that you can use to assign pseudo-counts to
    different matrices.  They should be matrices of the appropriate
    size that contain numbers to add to each parameter matrix.
    """
    N, M = len(states), len(alphabet)
    if pseudo_initial is not None:
        pseudo_initial = numpy.asarray(pseudo_initial)
        if pseudo_initial.shape != (N,):
            raise ValueError("pseudo_initial not shape len(states)")
    if pseudo_transition is not None:
        pseudo_transition = numpy.asarray(pseudo_transition)
        if pseudo_transition.shape != (N, N):
            raise ValueError("pseudo_transition not shape len(states) X len(states)")
    if pseudo_emission is not None:
        pseudo_emission = numpy.asarray(pseudo_emission)
        if pseudo_emission.shape != (N, M):
            raise ValueError("pseudo_emission not shape len(states) X len(alphabet)")

    # Training data is given as a list of members of the alphabet.
    # Replace those with indexes into the alphabet list for easier
    # computation.
    training_states, training_outputs = [], []
    states_indexes = itemindex(states)
    outputs_indexes = itemindex(alphabet)
    for toutputs, tstates in training_data:
        if len(tstates) != len(toutputs):
            raise ValueError("states and outputs not aligned")
        training_states.append([states_indexes[x] for x in tstates])
        training_outputs.append([outputs_indexes[x] for x in toutputs])

    x = _mle(
        N,
        M,
        training_outputs,
        training_states,
        pseudo_initial,
        pseudo_transition,
        pseudo_emission,
    )
    p_initial, p_transition, p_emission = x

    return MarkovModel(states, alphabet, p_initial, p_transition, p_emission)


def _mle(
    N,
    M,
    training_outputs,
    training_states,
    pseudo_initial,
    pseudo_transition,
    pseudo_emission,
):
    """Implement Maximum likelihood estimation algorithm (PRIVATE)."""
    # p_initial is the probability that a sequence of states starts
    # off with a particular one.
    p_initial = numpy.zeros(N)
    if pseudo_initial:
        p_initial = p_initial + pseudo_initial
    for states in training_states:
        p_initial[states[0]] += 1
    p_initial = _normalize(p_initial)

    # p_transition is the probability that a state leads to the next
    # one.  C(i,j)/C(i) where i and j are states.
    p_transition = numpy.zeros((N, N))
    if pseudo_transition:
        p_transition = p_transition + pseudo_transition
    for states in training_states:
        for n in range(len(states) - 1):
            i, j = states[n], states[n + 1]
            p_transition[i, j] += 1
    for i in range(len(p_transition)):
        p_transition[i, :] = p_transition[i, :] / sum(p_transition[i, :])

    # p_emission is the probability of an output given a state.
    # C(s,o)|C(s) where o is an output and s is a state.
    p_emission = numpy.zeros((N, M))
    if pseudo_emission:
        p_emission = p_emission + pseudo_emission
    p_emission = numpy.ones((N, M))
    for outputs, states in zip(training_outputs, training_states):
        for o, s in zip(outputs, states):
            p_emission[s, o] += 1
    for i in range(len(p_emission)):
        p_emission[i, :] = p_emission[i, :] / sum(p_emission[i, :])

    return p_initial, p_transition, p_emission


def _argmaxes(vector, allowance=None):
    """Return indices of the maximum values aong the vector (PRIVATE)."""
    return [numpy.argmax(vector)]


def find_states(markov_model, output):
    """Find states in the given Markov model output.

    Returns a list of (states, score) tuples.
    """
    mm = markov_model
    N = len(mm.states)

    # _viterbi does calculations in log space.  Add a tiny bit to the
    # matrices so that the logs will not break.
    lp_initial = numpy.log(mm.p_initial + VERY_SMALL_NUMBER)
    lp_transition = numpy.log(mm.p_transition + VERY_SMALL_NUMBER)
    lp_emission = numpy.log(mm.p_emission + VERY_SMALL_NUMBER)
    # Change output into a list of indexes into the alphabet.
    indexes = itemindex(mm.alphabet)
    output = [indexes[x] for x in output]

    # Run the viterbi algorithm.
    results = _viterbi(N, lp_initial, lp_transition, lp_emission, output)

    for i in range(len(results)):
        states, score = results[i]
        results[i] = [mm.states[x] for x in states], numpy.exp(score)
    return results


def _viterbi(N, lp_initial, lp_transition, lp_emission, output):
    """Implement Viterbi algorithm to find most likely states for a given input (PRIVATE)."""
    T = len(output)
    # Store the backtrace in a NxT matrix.
    backtrace = []  # list of indexes of states in previous timestep.
    for i in range(N):
        backtrace.append([None] * T)

    # Store the best scores.
    scores = numpy.zeros((N, T))
    scores[:, 0] = lp_initial + lp_emission[:, output[0]]
    for t in range(1, T):
        k = output[t]
        for j in range(N):
            # Find the most likely place it came from.
            i_scores = scores[:, t - 1] + lp_transition[:, j] + lp_emission[j, k]
            indexes = _argmaxes(i_scores)
            scores[j, t] = i_scores[indexes[0]]
            backtrace[j][t] = indexes

    # Do the backtrace.  First, find a good place to start.  Then,
    # we'll follow the backtrace matrix to find the list of states.
    # In the event of ties, there may be multiple paths back through
    # the matrix, which implies a recursive solution.  We'll simulate
    # it by keeping our own stack.
    in_process = []  # list of (t, states, score)
    results = []  # return values.  list of (states, score)
    indexes = _argmaxes(scores[:, T - 1])  # pick the first place
    for i in indexes:
        in_process.append((T - 1, [i], scores[i][T - 1]))
    while in_process:
        t, states, score = in_process.pop()
        if t == 0:
            results.append((states, score))
        else:
            indexes = backtrace[states[0]][t]
            for i in indexes:
                in_process.append((t - 1, [i] + states, score))
    return results


def _normalize(matrix):
    """Normalize matrix object (PRIVATE)."""
    if len(matrix.shape) == 1:
        matrix = matrix / sum(matrix)
    elif len(matrix.shape) == 2:
        # Normalize by rows.
        for i in range(len(matrix)):
            matrix[i, :] = matrix[i, :] / sum(matrix[i, :])
    else:
        raise ValueError("I cannot handle matrixes of that shape")
    return matrix


def _uniform_norm(shape):
    """Normalize a uniform matrix (PRIVATE)."""
    matrix = numpy.ones(shape)
    return _normalize(matrix)


def _random_norm(shape):
    """Normalize a random matrix (PRIVATE)."""
    matrix = numpy.random.random(shape)
    return _normalize(matrix)


def _copy_and_check(matrix, desired_shape):
    """Copy a matrix and check its dimension. Normalize at the end (PRIVATE)."""
    # Copy the matrix.
    matrix = numpy.array(matrix, copy=1)
    # Check the dimensions.
    if matrix.shape != desired_shape:
        raise ValueError("Incorrect dimension")
    # Make sure it's normalized.
    if len(matrix.shape) == 1:
        if numpy.fabs(sum(matrix) - 1.0) > 0.01:
            raise ValueError("matrix not normalized to 1.0")
    elif len(matrix.shape) == 2:
        for i in range(len(matrix)):
            if numpy.fabs(sum(matrix[i]) - 1.0) > 0.01:
                raise ValueError("matrix %d not normalized to 1.0" % i)
    else:
        raise ValueError("I don't handle matrices > 2 dimensions")
    return matrix


def _logsum(matrix):
    """Implement logsum for a matrix object (PRIVATE)."""
    if len(matrix.shape) > 1:
        vec = numpy.reshape(matrix, (numpy.product(matrix.shape),))
    else:
        vec = matrix
    sum = LOG0
    for num in vec:
        sum = logaddexp(sum, num)
    return sum


def _logvecadd(logvec1, logvec2):
    """Implement a log sum for two vector objects (PRIVATE)."""
    assert len(logvec1) == len(logvec2), "vectors aren't the same length"
    sumvec = numpy.zeros(len(logvec1))
    for i in range(len(logvec1)):
        sumvec[i] = logaddexp(logvec1[i], logvec2[i])
    return sumvec


def _exp_logsum(numbers):
    """Return the exponential of a logsum (PRIVATE)."""
    sum = _logsum(numbers)
    return numpy.exp(sum)
