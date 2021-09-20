# Copyright 2002 by Jeffrey Chang.
# All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Code for doing k-nearest-neighbors classification.

k Nearest Neighbors is a supervised learning algorithm that classifies
a new observation based the classes in its surrounding neighborhood.

Glossary:
 - distance   The distance between two points in the feature space.
 - weight     The importance given to each point for classification.

Classes:
 - kNN           Holds information for a nearest neighbors classifier.


Functions:
 - train        Train a new kNN classifier.
 - calculate    Calculate the probabilities of each class, given an observation.
 - classify     Classify an observation into a class.

Weighting Functions:
 - equal_weight    Every example is given a weight of 1.

"""

import numpy


class kNN:
    """Holds information necessary to do nearest neighbors classification.

    Attributes:
     - classes  Set of the possible classes.
     - xs       List of the neighbors.
     - ys       List of the classes that the neighbors belong to.
     - k        Number of neighbors to look at.
    """

    def __init__(self):
        """Initialize the class."""
        self.classes = set()
        self.xs = []
        self.ys = []
        self.k = None


def equal_weight(x, y):
    """Return integer one (dummy method for equally weighting)."""
    # everything gets 1 vote
    return 1


def train(xs, ys, k, typecode=None):
    """Train a k nearest neighbors classifier on a training set.

    xs is a list of observations and ys is a list of the class assignments.
    Thus, xs and ys should contain the same number of elements.  k is
    the number of neighbors that should be examined when doing the
    classification.
    """
    knn = kNN()
    knn.classes = set(ys)
    knn.xs = numpy.asarray(xs, typecode)
    knn.ys = ys
    knn.k = k
    return knn


def calculate(knn, x, weight_fn=None, distance_fn=None):
    """Calculate the probability for each class.

    Arguments:
     - x is the observed data.
     - weight_fn is an optional function that takes x and a training
       example, and returns a weight.
     - distance_fn is an optional function that takes two points and
       returns the distance between them.  If distance_fn is None (the
       default), the Euclidean distance is used.

    Returns a dictionary of the class to the weight given to the class.
    """
    if weight_fn is None:
        weight_fn = equal_weight

    x = numpy.asarray(x)

    order = []  # list of (distance, index)
    if distance_fn:
        for i in range(len(knn.xs)):
            dist = distance_fn(x, knn.xs[i])
            order.append((dist, i))
    else:
        # Default: Use a fast implementation of the Euclidean distance
        temp = numpy.zeros(len(x))
        # Predefining temp allows reuse of this array, making this
        # function about twice as fast.
        for i in range(len(knn.xs)):
            temp[:] = x - knn.xs[i]
            dist = numpy.sqrt(numpy.dot(temp, temp))
            order.append((dist, i))
    order.sort()

    # first 'k' are the ones I want.
    weights = {}  # class -> number of votes
    for k in knn.classes:
        weights[k] = 0.0
    for dist, i in order[: knn.k]:
        klass = knn.ys[i]
        weights[klass] = weights[klass] + weight_fn(x, knn.xs[i])

    return weights


def classify(knn, x, weight_fn=None, distance_fn=None):
    """Classify an observation into a class.

    If not specified, weight_fn will give all neighbors equal weight.
    distance_fn is an optional function that takes two points and returns
    the distance between them.  If distance_fn is None (the default),
    the Euclidean distance is used.
    """
    if weight_fn is None:
        weight_fn = equal_weight

    weights = calculate(knn, x, weight_fn=weight_fn, distance_fn=distance_fn)

    most_class = None
    most_weight = None
    for klass, weight in weights.items():
        if most_class is None or weight > most_weight:
            most_class = klass
            most_weight = weight
    return most_class
