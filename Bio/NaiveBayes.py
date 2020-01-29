# Copyright 2000 by Jeffrey Chang.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""General Naive Bayes learner.

Naive Bayes is a supervised classification algorithm that uses Bayes
rule to compute the fit between a new observation and some previously
observed data.  The observations are discrete feature vectors, with
the Bayes assumption that the features are independent.  Although this
is hardly ever true, the classifier works well enough in practice.

Glossary:
 - observation - A feature vector of discrete data.
 - class       - A possible classification for an observation.

Classes:
 - NaiveBayes - Holds information for a naive Bayes classifier.

Functions:
 - train     - Train a new naive Bayes classifier.
 - calculate - Calculate the probabilities of each class,
   given an observation.
 - classify  - Classify an observation into a class.

"""


try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.MaxEntropy."
    )


def _contents(items):
    """Return a dictionary where the key is the item and the value is the probablity associated (PRIVATE)."""
    term = 1.0 / len(items)
    counts = {}
    for item in items:
        counts[item] = counts.get(item, 0) + term
    return counts


class NaiveBayes:
    """Hold information for a NaiveBayes classifier.

    Attributes:
     - classes        - List of the possible classes of data.
     - p_conditional  - CLASS x DIM array of dicts of value -> ``P(value|class,dim)``
     - p_prior        - List of the prior probabilities for every class.
     - dimensionality - Dimensionality of the data.

    """

    def __init__(self):
        """Initialize the class."""
        self.classes = []
        self.p_conditional = None
        self.p_prior = []
        self.dimensionality = None


def calculate(nb, observation, scale=False):
    """Calculate the logarithmic conditional probability for each class.

    Arguments:
     - nb          - A NaiveBayes classifier that has been trained.
     - observation - A list representing the observed data.
     - scale       - Boolean to indicate whether the probability should be
       scaled by ``P(observation)``.  By default, no scaling is done.

    A dictionary is returned where the key is the class and the value is
    the log probability of the class.
    """
    # P(class|observation) = P(observation|class)*P(class)/P(observation)
    # Taking the log:
    # lP(class|observation) = lP(observation|class)+lP(class)-lP(observation)

    # Make sure the observation has the right dimensionality.
    if len(observation) != nb.dimensionality:
        raise ValueError(
            f"observation in {len(observation)} dimension,"
            f" but classifier in {nb.dimensionality}"
        )

    # Calculate log P(observation|class) for every class.
    n = len(nb.classes)
    lp_observation_class = numpy.zeros(n)  # array of log P(observation|class)
    for i in range(n):
        # log P(observation|class) = SUM_i log P(observation_i|class)
        probs = [None] * len(observation)
        for j in range(len(observation)):
            probs[j] = nb.p_conditional[i][j].get(observation[j], 0)
        lprobs = numpy.log(numpy.clip(probs, 1.0e-300, 1.0e300))
        lp_observation_class[i] = sum(lprobs)

    # Calculate log P(class).
    lp_prior = numpy.log(nb.p_prior)

    # Calculate log P(observation).
    lp_observation = 0.0  # P(observation)
    if scale:  # Only calculate this if requested.
        # log P(observation) = log SUM_i P(observation|class_i)P(class_i)
        obs = numpy.exp(numpy.clip(lp_prior + lp_observation_class, -700, +700))
        lp_observation = numpy.log(sum(obs))

    # Calculate log P(class|observation).
    lp_class_observation = {}  # Dict of class : log P(class|observation)
    for i in range(len(nb.classes)):
        lp_class_observation[nb.classes[i]] = (
            lp_observation_class[i] + lp_prior[i] - lp_observation
        )

    return lp_class_observation


def classify(nb, observation):
    """Classify an observation into a class."""
    # The class is the one with the highest probability.
    probs = calculate(nb, observation, scale=False)
    max_prob = max_class = None
    for klass in nb.classes:
        if max_prob is None or probs[klass] > max_prob:
            max_prob, max_class = probs[klass], klass
    return max_class


def train(training_set, results, priors=None, typecode=None):
    """Train a NaiveBayes classifier on a training set.

    Arguments:
     - training_set - List of observations.
     - results      - List of the class assignments for each observation.
       Thus, training_set and results must be the same length.
     - priors       - Optional dictionary specifying the prior probabilities
       for each type of result. If not specified, the priors will
       be estimated from the training results.

    """
    if not len(training_set):
        raise ValueError("No data in the training set.")
    if len(training_set) != len(results):
        raise ValueError("training_set and results should be parallel lists.")

    # If no typecode is specified, try to pick a reasonable one.  If
    # training_set is a Numeric array, then use that typecode.
    # Otherwise, choose a reasonable default.
    # XXX NOT IMPLEMENTED

    # Check to make sure each vector in the training set has the same
    # dimensionality.
    dimensions = [len(x) for x in training_set]
    if min(dimensions) != max(dimensions):
        raise ValueError("observations have different dimensionality")

    nb = NaiveBayes()
    nb.dimensionality = dimensions[0]

    # Get a list of all the classes, and
    # estimate the prior probabilities for the classes.
    if priors is not None:
        percs = priors
        nb.classes = list(set(results))
    else:
        class_freq = _contents(results)
        nb.classes = list(class_freq.keys())
        percs = class_freq
    nb.classes.sort()  # keep it tidy

    nb.p_prior = numpy.zeros(len(nb.classes))
    for i in range(len(nb.classes)):
        nb.p_prior[i] = percs[nb.classes[i]]

    # Collect all the observations in class.  For each class, make a
    # matrix of training instances versus dimensions.  I might be able
    # to optimize this with Numeric, if the training_set parameter
    # were guaranteed to be a matrix.  However, this may not be the
    # case, because the client may be hacking up a sparse matrix or
    # something.
    c2i = {}  # class to index of class
    for index, key in enumerate(nb.classes):
        c2i[key] = index
    observations = [[] for c in nb.classes]  # separate observations by class
    for i in range(len(results)):
        klass, obs = results[i], training_set[i]
        observations[c2i[klass]].append(obs)
    # Now make the observations Numeric matrix.
    for i in range(len(observations)):
        # XXX typecode must be specified!
        observations[i] = numpy.asarray(observations[i], typecode)

    # Calculate P(value|class,dim) for every class.
    # This is a good loop to optimize.
    nb.p_conditional = []
    for i in range(len(nb.classes)):
        class_observations = observations[i]  # observations for this class
        nb.p_conditional.append([None] * nb.dimensionality)
        for j in range(nb.dimensionality):
            # Collect all the values in this dimension.
            values = class_observations[:, j]

            # Add pseudocounts here.  This needs to be parameterized.
            # values = list(values) + range(len(nb.classes))  # XXX add 1

            # Estimate P(value|class,dim)
            nb.p_conditional[i][j] = _contents(values)
    return nb
