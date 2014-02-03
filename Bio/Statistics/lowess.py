# Copyright 2004-2008 by M de Hoon.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Implements the Lowess function for nonparametric regression.

Functions:
lowess        Fit a smooth nonparametric regression curve to a scatterplot.

For more information, see

William S. Cleveland: "Robust locally weighted regression and smoothing
scatterplots", Journal of the American Statistical Association, December 1979,
volume 74, number 368, pp. 829-836.

William S. Cleveland and Susan J. Devlin: "Locally weighted regression: An
approach to regression analysis by local fitting", Journal of the American
Statistical Association, September 1988, volume 83, number 403, pp. 596-610.
"""

from __future__ import print_function

from Bio._py3k import range

import numpy

try:
    from Bio.Cluster import median
    # The function median in Bio.Cluster is faster than the function median
    # in NumPy, as it does not require a full sort.
except ImportError as x:
    # Use the median function in NumPy if Bio.Cluster is not available
    from numpy import median


def lowess(x, y, f=2. / 3., iter=3):
    """lowess(x, y, f=2./3., iter=3) -> yest

    Lowess smoother: Robust locally weighted regression.
    The lowess function fits a nonparametric regression curve to a scatterplot.
    The arrays x and y contain an equal number of elements; each pair
    (x[i], y[i]) defines a data point in the scatterplot. The function returns
    the estimated (smooth) values of y.

    The smoothing span is given by f. A larger value for f will result in a
    smoother curve. The number of robustifying iterations is given by iter. The
    function will run faster with a smaller number of iterations.

    x and y should be numpy float arrays of equal length.  The return value is
    also a numpy float array of that length.

    e.g.
    >>> import numpy
    >>> x = numpy.array([4,  4,  7,  7,  8,  9, 10, 10, 10, 11, 11, 12, 12, 12,
    ...                 12, 13, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15, 16, 16,
    ...                 17, 17, 17, 18, 18, 18, 18, 19, 19, 19, 20, 20, 20, 20,
    ...                 20, 22, 23, 24, 24, 24, 24, 25], numpy.float)
    >>> y = numpy.array([2, 10,  4, 22, 16, 10, 18, 26, 34, 17, 28, 14, 20, 24,
    ...                 28, 26, 34, 34, 46, 26, 36, 60, 80, 20, 26, 54, 32, 40,
    ...                 32, 40, 50, 42, 56, 76, 84, 36, 46, 68, 32, 48, 52, 56,
    ...                 64, 66, 54, 70, 92, 93, 120, 85], numpy.float)
    >>> result = lowess(x, y)
    >>> len(result)
    50
    >>> print("[%0.2f, ..., %0.2f]" % (result[0], result[-1]))
    [4.85, ..., 84.98]
    """
    n = len(x)
    r = int(numpy.ceil(f * n))
    h = [numpy.sort(abs(x - x[i]))[r] for i in range(n)]
    w = numpy.clip(abs(([x] - numpy.transpose([x])) / h), 0.0, 1.0)
    w = 1 - w * w * w
    w = w * w * w
    yest = numpy.zeros(n)
    delta = numpy.ones(n)
    for iteration in range(iter):
        for i in range(n):
            weights = delta * w[:, i]
            weights_mul_x = weights * x
            b1 = numpy.dot(weights, y)
            b2 = numpy.dot(weights_mul_x, y)
            A11 = sum(weights)
            A12 = sum(weights_mul_x)
            A21 = A12
            A22 = numpy.dot(weights_mul_x, x)
            determinant = A11 * A22 - A12 * A21
            beta1 = (A22 * b1 - A12 * b2) / determinant
            beta2 = (A11 * b2 - A21 * b1) / determinant
            yest[i] = beta1 + beta2 * x[i]
        residuals = y - yest
        s = median(abs(residuals))
        delta[:] = numpy.clip(residuals / (6 * s), -1, 1)
        delta[:] = 1 - delta * delta
        delta[:] = delta * delta
    return yest


def _test():
    """Run the Bio.Statistics.lowess module's doctests."""
    print("Running doctests...")
    import doctest
    doctest.testmod()
    print("Done")

if __name__ == "__main__":
    _test()
