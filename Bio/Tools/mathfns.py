# Copyright 2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""This provides useful general math tools.

Functions:
safe_log   log, but returns an arbitrarily small number for log(0).
safe_exp   exp, but returns a large or small number instead of overflows.

"""

import math
from cmathfns import *

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
