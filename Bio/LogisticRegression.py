# Copyright 2002 by Jeffrey Chang.
# All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Code for doing logistic regressions (DEPRECATED).

Classes:
 - LogisticRegression    Holds information for a LogisticRegression classifier.

Functions:
 - train        Train a new classifier.
 - calculate    Calculate the probabilities of each class, given an observation.
 - classify     Classify an observation into a class.

This module has been deprecated, please consider an alternative like scikit-learn
insead.
"""

import warnings
from Bio import BiopythonDeprecationWarning

warnings.warn(
    "The 'Bio.LogisticRegression' module is deprecated and will be removed in a future "
    "release of Biopython. Consider using scikit-learn instead.",
    BiopythonDeprecationWarning,
)

try:
    import numpy as np
    import numpy.linalg
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Please install NumPy if you want to use Bio.LogisticRegression. "
        "See http://www.numpy.org/"
    ) from None


class LogisticRegression:
    """Holds information necessary to do logistic regression classification.

    Attributes:
     - beta - List of the weights for each dimension.

    """

    def __init__(self):
        """Initialize the class."""
        self.beta = []


def train(xs, ys, update_fn=None, typecode=None):
    """Train a logistic regression classifier on a training set.

    Argument xs is a list of observations and ys is a list of the class
    assignments, which should be 0 or 1.  xs and ys should contain the
    same number of elements.  update_fn is an optional callback function
    that takes as parameters that iteration number and log likelihood.
    """
    if len(xs) != len(ys):
        raise ValueError("xs and ys should be the same length.")
    classes = set(ys)
    if classes != {0, 1}:
        raise ValueError("Classes should be 0's and 1's")
    if typecode is None:
        typecode = "d"

    # Dimensionality of the data is the dimensionality of the
    # observations plus a constant dimension.
    N, ndims = len(xs), len(xs[0]) + 1
    if N == 0 or ndims == 1:
        raise ValueError("No observations or observation of 0 dimension.")

    # Make an X array, with a constant first dimension.
    X = np.ones((N, ndims), typecode)
    X[:, 1:] = xs
    Xt = np.transpose(X)
    y = np.asarray(ys, typecode)

    # Initialize the beta parameter to 0.
    beta = np.zeros(ndims, typecode)

    MAX_ITERATIONS = 500
    CONVERGE_THRESHOLD = 0.01
    stepsize = 1.0
    # Now iterate using Newton-Raphson until the log-likelihoods
    # converge.
    i = 0
    old_beta = old_llik = None
    while i < MAX_ITERATIONS:
        # Calculate the probabilities.  p = e^(beta X) / (1+e^(beta X))
        ebetaX = np.exp(np.dot(beta, Xt))
        p = ebetaX / (1 + ebetaX)

        # Find the log likelihood score and see if I've converged.
        logp = y * np.log(p) + (1 - y) * np.log(1 - p)
        llik = sum(logp)
        if update_fn is not None:
            update_fn(iter, llik)
        if old_llik is not None:
            # Check to see if the likelihood decreased.  If it did, then
            # restore the old beta parameters and half the step size.
            if llik < old_llik:
                stepsize /= 2.0
                beta = old_beta
            # If I've converged, then stop.
            if np.fabs(llik - old_llik) <= CONVERGE_THRESHOLD:
                break
        old_llik, old_beta = llik, beta
        i += 1

        W = np.identity(N) * p
        Xtyp = np.dot(Xt, y - p)  # Calculate the first derivative.
        XtWX = np.dot(np.dot(Xt, W), X)  # Calculate the second derivative.
        delta = numpy.linalg.solve(XtWX, Xtyp)
        if np.fabs(stepsize - 1.0) > 0.001:
            delta *= stepsize
        beta += delta  # Update beta.
    else:
        raise RuntimeError("Didn't converge.")

    lr = LogisticRegression()
    lr.beta = list(beta)
    return lr


def calculate(lr, x):
    """Calculate the probability for each class.

    Arguments:
     - lr is a LogisticRegression object.
     - x is the observed data.

    Returns a list of the probability that it fits each class.
    """
    # Insert a constant term for x.
    x = np.asarray([1.0] + x)
    # Calculate the probability.  p = e^(beta X) / (1+e^(beta X))
    ebetaX = np.exp(np.dot(lr.beta, x))
    p = ebetaX / (1 + ebetaX)
    return [1 - p, p]


def classify(lr, x):
    """Classify an observation into a class."""
    probs = calculate(lr, x)
    if probs[0] > probs[1]:
        return 0
    return 1
