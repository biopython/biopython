#!/usr/bin/env python

"""
This module provides code for doing logistic regressions.


Classes:
LogisticRegression    Holds information for a LogisticRegression classifier.


Functions:
train        Train a new classifier.
calculate    Calculate the probabilities of each class, given an observation.
classify     Classify an observation into a class.
"""
try:
    from Numeric import *
    from LinearAlgebra import *  # inverse
except ImportError, x:
    raise ImportError, "This module requires NumPy with the LinearAlgebra lib"

from Bio import listfns

class LogisticRegression:
    """Holds information necessary to do logistic regression
    classification.

    Members:
    beta    List of the weights for each dimension.

    """
    def __init__(self):
        """LogisticRegression()"""
        beta = []

def train(xs, ys, update_fn=None, typecode=None):
    """train(xs, ys[, update_fn]) -> LogisticRegression
    
    Train a logistic regression classifier on a training set.  xs is a
    list of observations and ys is a list of the class assignments,
    which should be 0 or 1.  xs and ys should contain the same number
    of elements.  update_fn is an optional callback function that
    takes as parameters that iteration number and log likelihood.
    
    """
    if len(xs) != len(ys):
        raise ValueError, "xs and ys should be the same length."
    if not xs or not xs[0]:
        raise ValueError, "No observations or observation of 0 dimension."
    classes = listfns.items(ys)
    classes.sort()
    if classes != [0, 1]:
        raise ValueError, "Classes should be 0's and 1's"
    if typecode is None:
        typecode = Float

    # Dimensionality of the data is the dimensionality of the
    # observations plus a constant dimension.
    N, ndims = len(xs), len(xs[0]) + 1

    # Make an X array, with a constant first dimension.
    X = ones((N, ndims), typecode)
    X[:, 1:] = xs
    Xt = transpose(X)
    y = asarray(ys, typecode)

    # Initialize the beta parameter to 0.
    beta = zeros(ndims, typecode)

    MAX_ITERATIONS = 500
    CONVERGE_THRESHOLD = 0.01
    stepsize = 1.0
    # Now iterate using Newton-Raphson until the log-likelihoods
    # converge.
    iter = 0
    old_beta = old_llik = None
    while iter < MAX_ITERATIONS:
        # Calculate the probabilities.  p = e^(beta X) / (1+e^(beta X))
        ebetaX = exp(dot(beta, Xt))
        p = ebetaX / (1+ebetaX)
        
        # Find the log likelihood score and see if I've converged.
        logp = y*log(p) + (1-y)*log(1-p)
        llik = sum(logp)
        if update_fn is not None:
            update_fn(iter, llik)
        # Check to see if the likelihood decreased.  If it did, then
        # restore the old beta parameters and half the step size.
        if llik < old_llik:
            stepsize = stepsize / 2.0
            beta = old_beta
        # If I've converged, then stop.
        if old_llik is not None and fabs(llik-old_llik) <= CONVERGE_THRESHOLD:
            break
        old_llik, old_beta = llik, beta
        iter += 1

        W = identity(N) * p
        Xtyp = dot(Xt, y-p)         # Calculate the first derivative.
        XtWX = dot(dot(Xt, W), X)   # Calculate the second derivative.
        #u, s, vt = singular_value_decomposition(XtWX)
        #print "U", u
        #print "S", s
        delta = dot(inverse(XtWX), Xtyp)
        if fabs(stepsize-1.0) > 0.001:
            delta = delta * stepsize
        beta = beta + delta                 # Update beta.
    else:
        raise AssertionError, "Didn't converge."

    lr = LogisticRegression()
    lr.beta = map(float, beta)   # Convert back to regular array.
    return lr

def calculate(lr, x):
    """calculate(lr, x) -> list of probabilities

    Calculate the probability for each class.  lr is a
    LogisticRegression object.  x is the observed data.  Returns a
    list of the probability that it fits each class.

    """
    # Insert a constant term for x.
    x = asarray([1.0] + x)
    # Calculate the probability.  p = e^(beta X) / (1+e^(beta X))
    ebetaX = exp(dot(lr.beta, x))
    p = ebetaX / (1+ebetaX)
    return [1-p, p]

def classify(lr, x):
    """classify(lr, x) -> 1 or 0

    Classify an observation into a class.

    """
    probs = calculate(lr, x)
    if probs[0] > probs[1]:
        return 0
    return 1
