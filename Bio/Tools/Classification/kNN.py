#!/usr/bin/env python

"""kNN.py

This module provides code for doing k-nearest-neighbors classification.

k Nearest Neighbors is a supervised learning algorithm that classifies
a new observation based the classes in its surrounding neighborhood.

Glossary:
distance   The distance between two points in the feature space.
weight     The importance given to each point for classification. 


Classes:
kNN           Holds information for a nearest neighbors classifier.


Functions:
train        Train a new kNN classifier.
calculate    Calculate the probabilities of each class, given an observation.
classify     Classify an observation into a class.

    Distance Functions:
euclidean_dist  The euclidean distance between two points.

    Weighting Functions:
equal_weight    Every example is given a weight of 1.

"""
import math

try:
    from Numeric import *
except ImportError, x:
    raise ImportError, "This module requires NumPy"

from Bio.Tools import listfns

class kNN:
    """Holds information necessary to do nearest neighbors classification.

    Members:
    classes  List of the possible classes.
    xs       List of the neighbors.
    ys       List of the classes that the neighbors belong to.
    k        Number of neighbors to look at.

    """
    def __init__(self):
        """kNN()"""
        self.classes = []
        self.xs = []
        self.ys = []
        self.k = None

def euclidean_dist(x, y):
    """euclidean_dist(x, y) -> euclidean distance between x and y"""
    if len(x) != len(y):
        raise ValueError, "vectors must be same length"
    return math.sqrt(sum((x-y)**2))

def euclidean_dist_py(x, y):
    """euclidean_dist_py(x, y) -> euclidean distance between x and y"""
    # lightly modified from implementation by Thomas Sicheritz-Ponten.
    # This works faster than the Numeric implementation on shorter
    # vectors.
    if len(x) != len(y):
        raise ValueError, "vectors must be same length"
    sum = 0
    for i in range(len(x)):
        sum += (x[i]-y[i])**2
    return math.sqrt(sum)

def equal_weight(x, y):
    """equal_weight(x, y) -> 1"""
    # everything gets 1 vote
    return 1

def train(xs, ys, k, typecode=None):
    """train(xs, ys, k) -> kNN
    
    Train a k nearest neighbors classifier on a training set.  xs is a
    list of observations and ys is a list of the class assignments.
    Thus, xs and ys should contain the same number of elements.  k is
    the number of neighbors that should be examined when doing the
    classification.
    
    """
    knn = kNN()
    knn.classes = listfns.items(ys)
    knn.xs = asarray(xs, typecode)
    knn.ys = ys
    knn.k = k
    return knn

def calculate(knn, x, weight_fn=equal_weight, distance_fn=euclidean_dist):
    """calculate(knn, x[, weight_fn][, distance_fn]) -> weight dict

    Calculate the probability for each class.  knn is a kNN object.  x
    is the observed data.  weight_fn is an optional function that
    takes x and a training example, and returns a weight.  distance_fn
    is an optional function that takes two points and returns the
    distance between them.  Returns a dictionary of the class to the
    weight given to the class.
    
    """
    x = asarray(x)

    order = []  # list of (distance, index)
    for i in range(len(knn.xs)):
        dist = distance_fn(x, knn.xs[i])
        order.append((dist, i))
    order.sort()

    # first 'k' are the ones I want.
    weights = {}  # class -> number of votes
    for k in knn.classes:
        weights[k] = 0.0
    for dist, i in order[:knn.k]:
        klass = knn.ys[i]
        weights[klass] = weights[klass] + weight_fn(x, knn.xs[i])

    return weights

def classify(knn, x, weight_fn=equal_weight, distance_fn=euclidean_dist):
    """classify(knn, x[, weight_fn][, distance_fn]) -> class

    Classify an observation into a class.  If not specified, weight_fn will
    give all neighbors equal weight and distance_fn will be the euclidean
    distance.

    """
    weights = calculate(
        knn, x, weight_fn=weight_fn, distance_fn=distance_fn)

    most_class = None
    most_weight = None
    for klass, weight in weights.items():
        if most_class is None or weight > most_weight:
            most_class = klass
            most_weight = weight
    return most_class
