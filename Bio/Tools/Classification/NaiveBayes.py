# Copyright 2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""This provides code for a general Naive Bayes learner.

Naive Bayes is a supervised classification algorithm that is simple to
implement and performs reasonably well in practice.  It is used to
classify multi-dimensional observations into 1 of several classes, with
the assumptions that the classes are mutually exclusive and comprehensive.

In this file, an observation is represented as a list of discrete data
points.  Discrete data is data that can be represented as one of a finite,
usually small, number of values.  Examples include text words or bins in a
bucket.  A class is a possible classification of an observation.  Since
'class' is a reserved word in Python, variables are instead named
'class_'.


Classes:
NaiveBayes   Holds information for a naive bayes classifier.

Functions:
train        Train a new NaiveBayes classifier.
calculate    Calculate the probabilities of each class, given an observation.
classify     Classify an observation into a class.

_estimate_probs  Estimate the probabilities for each value in a data set.
_count           Count the number of occurrences of each value in a data set.
_safe_log        Calculate a natural log, with some safety checks.
_uniq            Return the unique members of a list

"""
# To Do:
# add code to help discretize data

import math

class NaiveBayes:
    """Holds information for a NaiveBayes classifier.

    Members:
    classes         List of the possible classes of data.
    p_conditional   Dict with keys (value, dim, class) -> P(value,dim|class).
    p_prior         Dict with the class as the key and P(class) as value.
    dimensionality  Dimensionality of the data.

    """
    def __init__(self):
        self.classes = {}
        self.p_conditional = {}
        self.p_prior = {}
        self.dimensionality = None

def calculate(nb, observation, scale=0):
    """calculate(nb, observation, scale=0) -> dict of probabilities

    Calculate log P(class|observation) for each class.  nb is a NaiveBayes
    classifier that has been trained.  observation is a list representing
    the observed data.  scale is whether the probability should be
    scaled by P(observation).  Returns a dictionary where the keys is the
    class and the value is the log probability of the class.

    """
    # P(observation|class) = P(class|observation)*P(observation)/P(class)
    # Taking the log:
    # lP(observation|class) = lP(class|observation)+P(observation)-P(class)

    # Make sure the observation has the right dimensionality.
    if len(observation) != nb.dimensionality:
        raise ValueError, "observation in %d dimension, but classifier in %d" \
              % (len(observation), nb.dimensionality)

    # Calculate log P(observation|class) for every class.
    lp_observation_class = {}     # Dict of class : log P(observation|class)
    for class_ in nb.classes:
        # log P(observation|class) = SUM_i log P(observation_i|class)
        lprob = 0.0
        for dim in range(len(observation)):
            value = observation[dim]
            p = nb.p_conditional.get((value, dim, class_), 0.0)
            lprob = lprob + _safe_log(p)
        lp_observation_class[class_] = lprob

    # Calculate log P(class).
    lp_prior = {}                 # Dict of class : log P(class)
    for class_ in nb.classes:
        lp_prior[class_] = math.log(nb.p_prior[class_])

    # Calculate log P(observation).
    lp_observation = 0.0          # P(observation)
    if scale:   # Only calculate this if requested.
        # log P(observation) = log SUM_i P(observation|class_i)P(class_i)
        sum = 0.0
        for class_ in nb.classes:
            p = math.exp(lp_prior[class_] + lp_observation_class[class_])
            sum = sum + p
        lp_observation = math.log(sum)

    # Calculate log P(class|observation).
    lp_class_observation = {}      # Dict of class : log P(class|observation)
    for class_ in nb.classes:
        lp_class_observation[class_] = \
            lp_observation_class[class_] + lp_prior[class_] - lp_observation

    return lp_class_observation

def classify(nb, observation):
    """classify(nb, observation) -> class

    Classify an observation into a class.

    """
    # The class is the one with the highest probability.
    probs = calculate(nb, observation, scale=0)
    max_prob = None
    max_class = None
    for class_ in nb.classes:
        if max_prob is None or probs[class_] > max_prob:
            max_prob = probs[class_]
            max_class = class_
    return max_class

def train(training_set, results):
    """train(training_set, results) -> NaiveBayes

    Train a naive bayes classifier on a training set.  training_set is a
    list of observations.  results is a list of the class assignments
    for each observation.  Thus, training_set and results must be the same
    length.

    """
    if not len(training_set):
        raise ValueError, "No data in the training set."
    if len(training_set) != len(results):
        raise ValueError, "training_set and results should be parallel lists."
    dim = None
    for obs in training_set:
        if dim is None:
            dim = len(obs)
        elif dim != len(obs):
            raise ValueError, "observations have different dimensionality"
    
    nb = NaiveBayes()
    nb.dimensionality = len(training_set[0])
    
    # Get a list of all the classes.
    nb.classes = _uniq(results)
    
    # Estimate the prior probabilities for the classes.
    nb.p_prior = _estimate_probs(results)

    # Calculate P(value,dim|class) for every class.
    for class_ in nb.classes:
        # Collect all the observations in class_.
        observations = []
        for i in range(len(results)):
            if results[i] == class_:
                observations.append(training_set[i])

        for dim in range(nb.dimensionality):
            # Collect all the values in this dimension.
            values = map(lambda x,dim=dim: x[dim], observations)
            # Estimate P(value,dim|class)
            p_values = _estimate_probs(values)
            for value in p_values.keys():
                nb.p_conditional[(value, dim, class_)] = p_values[value]
    return nb

def _count(values):
    """_count(values) -> dict of counts 

    Count the number of times each unique value appears in a list of data.

    """
    counts = {}
    for value in values:
        counts[value] = counts.get(value, 0) + 1
    return counts

def _estimate_probs(values):
    """_estimate_probs(values) -> dict of probabilities

    Calculate the major likelihood estimate of probability for each
    unique value in a list of data.

    """
    # Divide the number of time an observation occurs by the total
    # number of observations.
    total = float(len(values))
    counts = _count(values)
    probs = {}
    for value in counts.keys():
        probs[value] = counts[value] / total
    return probs

def _safe_log(n):
    """_safe_log(n) -> log of n

    Calculate the natural log of n.  If n is 0, then calculates the log
    of a really small number.

    """
    try:
        return math.log(n)
    except OverflowError:
        # Return the log of an arbitrarily low number.
        return math.log(1E-100)
    assert 0, "How did I get here?"

def _uniq(l):
    """uniq(l) -> list of unique members of l"""
    dict = {}
    for i in l:
        dict[i] = 1
    return dict.keys()

