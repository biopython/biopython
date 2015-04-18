"""Python implementation of chisqprob, to avoid SciPy dependency.

Adapted from SciPy: scipy/special/cephes/{chdtr,igam}.
"""

import math
import sys

__docformat__ = "restructuredtext en"

# Cephes Math Library Release 2.0:  April, 1987
# Copyright 1985, 1987 by Stephen L. Moshier
# Direct inquiries to 30 Frost Street, Cambridge, MA 02140
MACHEP = 0.0000001     # the machine roundoff error / tolerance
BIG = 4.503599627370496e15
BIGINV = 2.22044604925031308085e-16


def chisqprob(x, df):
    """Probability value (1-tail) for the Chi^2 probability distribution.

    Broadcasting rules apply.

    Parameters
    ----------
    x : array_like or float > 0

    df : array_like or float, probably int >= 1

    Returns
    -------
    chisqprob : ndarray
        The area from `chisq` to infinity under the Chi^2 probability
        distribution with degrees of freedom `df`.

    """
    if x <= 0:
        return 1.0
    if x == 0:
        return 0.0
    if df <= 0:
        raise ValueError("Domain error.")
    if x < 1.0 or x < df:
        return 1.0 - _igam(0.5 * df, 0.5 * x)
    return _igamc(0.5 * df, 0.5 * x)


def _igamc(a, x):
    """Complemented incomplete Gamma integral.

    SYNOPSIS:

    double a, x, y, igamc();

    y = igamc( a, x );

    DESCRIPTION:

    The function is defined by::

        igamc(a,x)   =   1 - igam(a,x)

                                inf.
                                   -
                          1       | |  -t  a-1
                    =   -----     |   e   t   dt.
                         -      | |
                        | (a)    -
                                    x

    In this implementation both arguments must be positive.
    The integral is evaluated by either a power series or
    continued fraction expansion, depending on the relative
    values of a and x.
    """
    # Compute  x**a * exp(-x) / Gamma(a)
    ax = math.exp(a * math.log(x) - x - math.lgamma(a))

    # Continued fraction
    y = 1.0 - a
    z = x + y + 1.0
    c = 0.0
    pkm2 = 1.0
    qkm2 = x
    pkm1 = x + 1.0
    qkm1 = z * x
    ans = pkm1 / qkm1
    while True:
        c += 1.0
        y += 1.0
        z += 2.0
        yc = y * c
        pk = pkm1 * z - pkm2 * yc
        qk = qkm1 * z - qkm2 * yc
        if qk != 0:
            r = pk / qk
            t = abs((ans - r) / r)
            ans = r
        else:
            t = 1.0
        pkm2 = pkm1
        pkm1 = pk
        qkm2 = qkm1
        qkm1 = qk
        if abs(pk) > BIG:
                pkm2 *= BIGINV
                pkm1 *= BIGINV
                qkm2 *= BIGINV
                qkm1 *= BIGINV
        if t <= MACHEP:
            return ans * ax


def _igam(a, x):
    """Left tail of incomplete Gamma function.

    Computes this formula::

                 inf.      k
          a  -x   -       x
         x  e     >   ----------
                  -     -
                k=0   | (a+k+1)
    """

    # Compute  x**a * exp(-x) / Gamma(a)
    ax = math.exp(a * math.log(x) - x - math.lgamma(a))

    # Power series
    r = a
    c = 1.0
    ans = 1.0

    while True:
        r += 1.0
        c *= x / r
        ans += c
        if c / ans <= MACHEP:
            return ans * ax / a


# --- Speed ---

# try:
#    from scipy.stats import chisqprob
# except ImportError:
#    pass
