# Copyright 2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""This provides useful general math tools (OBSOLETE).

This module and its C code equivalent are considered to be obsolete, and
are likely to be deprecated in a future release of Biopython, before being
removed.  Please get in touch via the mailing list if this will affect you.

Functions:
fcmp       Compare two floating point numbers, up to a specified precision.
intd       Represent a floating point number as an integer.
safe_log   log, but returns an arbitrarily small number for log(0).
safe_exp   exp, but returns a large or small number instead of overflows.

"""
import math

def fcmp(x, y, precision):
    """fcmp(x, y, precision) -> -1, 0, or 1"""
    if math.fabs(x-y) < precision:
        return 0
    elif x < y:
        return -1
    return 1

def intd(x, digits_after_decimal=0):
    """intd(x[, digits_after_decimal]) -> int x, rounded

    Represent a floating point number with some digits after the
    decimal point as an integer.  This is useful when floating point
    comparisons are failing due to precision problems.  e.g.
    intd(5.35, 1) -> 54.

    """
    precision = 10.**digits_after_decimal
    if x >= 0:
        x = int(x * precision + 0.5)
    else:
        x = int(x * precision - 0.5)
    return x

def safe_log(n, zero=None, neg=None):
    """safe_log(n, zero=None, neg=None) -> log(n)

    Calculate the log of n.  If n is 0, returns the value of zero.  If n is
    negative, returns the value of neg.

    """
    if n < 0:
        return neg
    elif n < 1E-100:
        return zero
    return math.log(n)

LOG2 = math.log(2)
def safe_log2(n, zero=None, neg=None):
    """safe_log2(n, zero=None, neg=None) -> log(n)

    Calculate the log base 2 of n.  If n is 0, returns the value of
    zero.  If n is negative, returns the value of neg.

    """
    l = safe_log(n, zero=zero, neg=neg)
    if l is None:
        return l
    return l/LOG2

def safe_exp(n, under=None, over=None):
    """safe_exp(n, under=None, over=None) -> e**n

    Guaranteed not to overflow.  Instead of overflowing, it returns
    the values of 'under' for underflows or 'over' for overflows.

    """
    try:
        return math.exp(n)
    except OverflowError:
        if n < 0:
            return under
        return over
    raise "How did I get here?"

# Try and load C implementations of functions.  If I can't,
# then just ignore and use the pure python implementations.
try:
    from cmathfns import *
except ImportError:
    pass
