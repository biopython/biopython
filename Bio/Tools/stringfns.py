# Copyright 2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""This provides useful general functions for working with strings.

Functions:
split        Split a string using many delimiters.
starts_with  Check whether a string starts with another string.

"""

def split(s, sep=" \011\012\013\014\015", maxsplit=None, negate=0):
    """split(str [,sep [,maxsplit [,negate]]]) -> list of strings

    Split a string.  Similar to string.split, except that this considers
    any one of the characters in sep to be a delimiter.  If negate is
    true, then everything but sep will be a separator.

    """
    strlist = []
    prev = 0
    for i in range(len(s)):
        if maxsplit is not None and len(strlist) >= maxsplit:
            break
        if (s[i] in sep) == (not negate):
            strlist.append(s[prev:i])
            prev = i+1
    strlist.append(s[prev:])
    return strlist

def starts_with(s, start):
    """starts_with(s, start) -> 1/0

    Return whether s begins with start.

    """
    return s[:len(start)] == start

# Try and load C implementations of functions.  If I can't,
# then just ignore and use the pure python implementations.
try:
    from cstringfns import *
except ImportError:
    pass
