# Copyright 2024 by the Biopython team. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Rate equations for chemical kinetics.

This module provides functions for calculating concentrations over
time for common reaction orders, and the Arrhenius equation for
temperature-dependent rate constants.

Functions
---------
- firstOrder  : Concentration at time t for a first-order reaction.
- secondOrder : Concentration at time t for a second-order reaction.
- zeroOrder   : Concentration at time t for a zero-order reaction.
- Arrenheius  : Rate constant at temperature T via the Arrhenius equation.
"""

import math

from scipy.constants import gas_constant


def firstOrder(k, A0, t):
    """Concentration at time t for a first-order reaction.

    For a first-order reaction A -> products, the concentration of
    reactant A at time t is given by:

        A(t) = A0 * exp(-k * t)

    Parameters
    ----------
    k : float
        Rate constant (units of 1/time).
    A0 : float
        Initial concentration of reactant.
    t : float
        Time at which to evaluate the concentration.

    Returns
    -------
    float
        Concentration of reactant at time t.

    Examples
    --------
    >>> from Bio.Kinetics.Arrhenius import firstOrder
    >>> firstOrder(0.1, 100, 10)
    36.787944117144235
    """
    return A0 * (math.exp(-1 * t * k))


def secondOrder(k, A0, t):
    """Concentration at time t for a second-order reaction.

    For a second-order reaction 2A -> products (or A + B -> products
    with equal initial concentrations), the concentration of reactant
    at time t is given by:

        A(t) = 1 / (1/A0 + k * t)

    Parameters
    ----------
    k : float
        Rate constant (units of 1/(concentration * time)).
    A0 : float
        Initial concentration of reactant. Must be non-zero.
    t : float
        Time at which to evaluate the concentration.

    Returns
    -------
    float
        Concentration of reactant at time t.

    Raises
    ------
    ZeroDivisionError
        If A0 is zero.

    Examples
    --------
    >>> from Bio.Kinetics.Arrhenius import secondOrder
    >>> secondOrder(0.01, 50, 5)
    14.285714285714285
    """
    return 1 / ((1 / A0) + k * t)


def zeroOrder(k, A0, t):
    """Concentration at time t for a zero-order reaction.

    For a zero-order reaction A -> products, the concentration of
    reactant at time t decreases linearly:

        A(t) = A0 - k * t

    Parameters
    ----------
    k : float
        Rate constant (units of concentration/time).
    A0 : float
        Initial concentration of reactant.
    t : float
        Time at which to evaluate the concentration.

    Returns
    -------
    float
        Concentration of reactant at time t. May be negative
        if t exceeds the depletion time (A0 / k).

    Examples
    --------
    >>> from Bio.Kinetics.Arrhenius import zeroOrder
    >>> zeroOrder(0.5, 100, 30)
    85.0
    """
    return (-1 * k * t) + A0


def Arrenheius(A, Ea, T):
    """Rate constant at temperature T via the Arrhenius equation.

    The Arrhenius equation describes the temperature dependence of
    reaction rates:

        k = A * exp(-Ea / (R * T))

    where R is the universal gas constant (8.314462618 J/(mol*K)).

    Parameters
    ----------
    A : float
        Pre-exponential factor (frequency factor) in same units as
        the desired rate constant.
    Ea : float
        Activation energy in J/mol.
    T : float
        Absolute temperature in Kelvin.

    Returns
    -------
    float
        The rate constant at temperature T.

    Examples
    --------
    >>> from Bio.Kinetics.Arrhenius import Arrenheius
    >>> k = Arrenheius(1e13, 75000, 298.15)
    >>> print(f"{k:.3f}")
    0.725
    """
    return A * math.exp(-1 * Ea / (gas_constant * T))
