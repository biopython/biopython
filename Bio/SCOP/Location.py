# Copyright 2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Location.py

This module provides code to work with protein locations described in
SCOP files.

parse    Parse a location string into a list of (chain, start, end).
str      Represent a location as a string.

"""
import string
import re

_location_re = re.compile(r"(\w:)?(\w*)-?(\w*)")
def parse(str):
    """parse(str) -> list of tuples (chain, start, end).

    str is the string that describes the location of the domain.

    """
    # The format of this is not described, but here are some examples:
    # f:
    # 87-116
    # a:103-131
    # c:,d:
    # 1:1-120,1:255-349
    # 4a-99a
    if str == '-':  # no location, whole sequence
        return([])
    locations = []
    for l in string.split(str, ","):
        m = _location_re.match(l)
        if m is None:
            raise SyntaxError, "I don't understand the format of %s" % l
        chain, start, end = m.groups()

        if chain:
            if chain[-1] != ':':
                raise SyntaxError, "I don't understand the chain in %s" % l
            chain = chain[:-1]   # chop off the ':'
        try:
            start = int(start)
            end = int(end)
        except ValueError:  # can't convert to int
            pass
        locations.append((chain, start, end))
    return locations

def str(locations):
    """str(locations) -> string

    Represent as a string.  locations is a list of tuples (chain, start, end).

    """
    if not locations:
        return '-'
    strs = []
    for chain, start, end in locations:
        s = []
        if chain:
            s.append("%s:" % chain)
        if start:
            s.append("%s-%s" % (start, end))
        strs.append(string.join(s, ''))
    return string.join(strs, ',')
