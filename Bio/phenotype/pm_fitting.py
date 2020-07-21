# Copyright 2014-2016 by Marco Galardini.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Growth curves fitting and parameters extraction for phenotype data.

This module provides functions to perform sigmoid functions fitting to
Phenotype Microarray data. This module depends on scipy curve_fit function.
If not available, a warning is raised.

Functions:
logistic           Logistic growth model.
gompertz           Gompertz growth model.
richards           Richards growth model.
guess_plateau      Guess the plateau point to improve sigmoid fitting.
guess_lag          Guess the lag point to improve sigmoid fitting.
fit                Sigmoid functions fit.
get_area           Calculate the area under the PM curve.
"""

import numpy as np

try:
    from scipy.optimize.minpack import curve_fit
    from scipy.integrate import trapz
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError("Install scipy to extract curve parameters.")


def logistic(x, A, u, d, v, y0):
    """Logistic growth model.

    Proposed in Zwietering et al., 1990 (PMID: 16348228)
    """
    y = (A / (1 + np.exp((((4 * u) / A) * (d - x)) + 2))) + y0
    return y


def gompertz(x, A, u, d, v, y0):
    """Gompertz growth model.

    Proposed in Zwietering et al., 1990 (PMID: 16348228)
    """
    y = (A * np.exp(-np.exp((((u * np.e) / A) * (d - x)) + 1))) + y0
    return y


def richards(x, A, u, d, v, y0):
    """Richards growth model (equivalent to Stannard).

    Proposed in Zwietering et al., 1990 (PMID: 16348228)
    """
    y = (
        A
        * pow(
            1
            + (
                v
                + (np.exp(1 + v) * np.exp((u / A) * (1 + v) * (1 + (1 / v)) * (d - x)))
            ),
            -(1 / v),
        )
    ) + y0
    return y


def guess_lag(x, y):
    """Given two axes returns a guess of the lag point.

    The lag point is defined as the x point where the difference in y
    with the next point is higher then the mean differences between
    the points plus one standard deviation. If such point is not found
    or x and y have different lengths the function returns zero.
    """
    if len(x) != len(y):
        return 0

    diffs = []
    indexes = range(len(x))

    for i in indexes:
        if i + 1 not in indexes:
            continue
        diffs.append(y[i + 1] - y[i])
    diffs = np.array(diffs)

    flex = x[-1]
    for i in indexes:
        if i + 1 not in indexes:
            continue
        if (y[i + 1] - y[i]) > (diffs.mean() + (diffs.std())):
            flex = x[i]
            break

    return flex


def guess_plateau(x, y):
    """Given two axes returns a guess of the plateau point.

    The plateau point is defined as the x point where the y point
    is near one standard deviation of the differences between the y points to
    the maximum y value. If such point is not found or x and y have
    different lengths the function returns zero.
    """
    if len(x) != len(y):
        return 0

    diffs = []
    indexes = range(len(y))

    for i in indexes:
        if i + 1 not in indexes:
            continue
        diffs.append(y[i + 1] - y[i])
    diffs = np.array(diffs)

    ymax = y[-1]
    for i in indexes:
        if y[i] > (ymax - diffs.std()) and y[i] < (ymax + diffs.std()):
            ymax = y[i]
            break

    return ymax


def fit(function, x, y):
    """Fit the provided function to the x and y values.

    The function parameters and the parameters covariance.
    """
    # Compute guesses for the parameters
    # This is necessary to get significant fits
    p0 = [guess_plateau(x, y), 4.0, guess_lag(x, y), 0.1, min(y)]

    params, pcov = curve_fit(function, x, y, p0=p0)
    return params, pcov


def get_area(y, x):
    """Get the area under the curve."""
    return trapz(y=y, x=x)
