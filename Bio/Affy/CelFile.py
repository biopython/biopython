# Copyright 2004 by Harry Zuzan.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
Classes for accessing the information in Affymetrix cel files.

Functions:
read      Read a cel file and store its contents in a Record

Classes:
Record    Contains the information from a cel file
"""

import numpy


class Record(object):
    """
    Stores the information in a cel file

    >>> with open('Tests/Affy/affy_v3_ex.CEL.bz2', 'rb') as handle:
    ...     cel_data = bz2.open('Tests/Affy/affy_v3_ex.CEL.bz2', 'rt', encoding='ascii')
    ...
    >>> c = CelFile.read(cel_data)
    >>> print(c.ncols, c.nrows)
    478 478
    >>> print(c.intensities)
    [[   234.  22428.    249. ...,  19179.    229.  19592.]
     [ 22730.    302.  22274. ...,    266.  18839.    241.]
     [   284.  22154.    287. ...,  18696.    332.  20289.]
     ...,
     [ 20743.    364.  19723. ...,    242.  19834.    220.]
     [   329.  20431.    368. ...,  19328.    266.  19990.]
     [ 20522.    325.  20713. ...,    224.  19963.    207.]]
    >>> print(c.stdevs)
    [[   24.   2378.1    27.7 ...,  2723.4    62.5  2771.6]
     [ 2075.3    46.7  2633.9 ...,    48.3  2208.4    49.5]
     [   73.3  2327.6    35.  ...,  2856.6    64.8  2280.5]
     ...,
     [ 3147.3    63.2  2707.7 ...,    52.9  2557.9    89.9]
     [   49.3  2170.9    73.7 ...,  3128.9    53.   3018.2]
     [ 2913.3    38.7  2341.4 ...,    50.6  3532.9    45.8]]
    >>> print(c.npix)
    [[25 25 25 ..., 25 25 25]
     [25 25 25 ..., 25 25 25]
     [25 25 25 ..., 25 25 25]
     ...,
     [25 25 25 ..., 25 25 25]
     [25 25 25 ..., 25 25 25]
     [25 25 25 ..., 25 25 25]]

    """
    def __init__(self):
        self.intensities = None
        self.stdevs = None
        self.npix = None
        self.nrows = None
        self.ncols = None


def read(handle):
    """
    Read the information in a cel file, and store it in a Record.
    """
    # Needs error handling.
    # Needs to know the chip design.
    record = Record()
    section = ""
    for line in handle:
        if not line.strip():
            continue
        if line[:8] == "[HEADER]":
            section = "HEADER"
        elif line[:11] == "[INTENSITY]":
            section = "INTENSITY"
            record.intensities = numpy.zeros((record.nrows, record.ncols))
            record.stdevs = numpy.zeros((record.nrows, record.ncols))
            record.npix = numpy.zeros((record.nrows, record.ncols), int)
        elif line[0] == "[":
            section = ""
        elif section == "HEADER":
            keyword, value = line.split("=", 1)
            if keyword == "Cols":
                record.ncols = int(value)
            elif keyword == "Rows":
                record.nrows = int(value)
        elif section == "INTENSITY":
            if "=" in line:
                continue
            words = line.split()
            y = int(words[0])
            x = int(words[1])
            record.intensities[x, y] = float(words[2])
            record.stdevs[x, y] = float(words[3])
            record.npix[x, y] = int(words[4])
    return record
