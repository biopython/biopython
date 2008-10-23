"""
This module provides code for various distance measures.

Functions:
euclidean       Euclidean distance between two points
euclidean_py    Pure Python implementation of euclidean.

"""
# XXX cosine distance

import warnings
warnings.warn("Bio.distance is deprecated. If you use this module, please notify the Biopython developers at biopython-dev@biopython.org", DeprecationWarning)

from numpy import *

def euclidean(x, y):
    """euclidean(x, y) -> euclidean distance between x and y"""
    if len(x) != len(y):
        raise ValueError("vectors must be same length")
    #return sqrt(sum((x-y)**2))
    # Optimization by John Corradi (JCorradi@msn.com)
    d = x-y
    return sqrt(dot(d, d))

def euclidean_py(x, y):
    """euclidean_py(x, y) -> euclidean distance between x and y"""
    # lightly modified from implementation by Thomas Sicheritz-Ponten.
    # This works faster than the Numeric implementation on shorter
    # vectors.
    if len(x) != len(y):
        raise ValueError("vectors must be same length")
    sum = 0
    for i in range(len(x)):
        sum += (x[i]-y[i])**2
    return sqrt(sum)
