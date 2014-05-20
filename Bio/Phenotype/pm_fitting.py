# Copyright 2014 by Marco Galardini.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""This module provides functions to perform sigmoid functions fitting to
Phenotype Microarray data. This module depends on scipy curve_fit function.
If not available, a warning is raised.

Functions:
logistic           Logistic growth model.
gompertz           Gompertz growth model.
richards           Richards growth model.
fit                Sigmoid functions fit.
get_area           Calculate the area under the PM curve."""

import warnings
import numpy as np

from scipy.optimize.minpack import curve_fit
from scipy.integrate import trapz

def logistic(x, A, u, d, v, y0):
    """Logistic growth model
    
    Proposed in Zwietering et al., 1990 (PMID: 16348228)
    """
    y = (A / (1 + np.exp( ( ((4 * u)/A) * (d - x) ) + 2 ))) + y0
    return y
    
def gompertz(x, A, u, d, v, y0):
    """Gompertz growth model
    
    Proposed in Zwietering et al., 1990 (PMID: 16348228)
    """
    y = (A * np.exp( -np.exp( (((u * np.e)/A) * (d - x)) + 1 ) ) ) + y0
    return y
    
def richards(x, A, u, d, v, y0):
    """Gompertz growth model (equivalent to Stannard)
    
    Proposed in Zwietering et al., 1990 (PMID: 16348228)
    """
    y = (A * pow(1 + (v + (np.exp(1 + v) * np.exp( (u/A) * (1 + v) * (1 + (1/v)) * (d - x) ) ) ),-(1/v))) + y0
    return y
    
def fit(function, x, y):
    """Fit the provided functrion to the x and y values.
    
    The function parameters and the parameters covariance."""
    params, pcov = curve_fit(function, x, y)
    return params, pcov
    
def get_area(y, x):
    """Get the area under the curve"""
    return trapz(y = y, x = x)
