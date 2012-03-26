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
        if line[:8]=="[HEADER]":
            section = "HEADER"
        elif line[:11]=="[INTENSITY]":
            section = "INTENSITY"
            record.intensities  = numpy.zeros((record.nrows, record.ncols))
            record.stdevs = numpy.zeros((record.nrows, record.ncols))
            record.npix = numpy.zeros((record.nrows, record.ncols), int)
        elif line[0]=="[":
            section = ""
        elif section=="HEADER":
            keyword, value = line.split("=", 1)
            if keyword=="Cols":
                record.ncols = int(value)
            elif keyword=="Rows":
                record.nrows = int(value)
        elif section=="INTENSITY":
            if "=" in line:
                continue
            words = line.split()
            y, x = map(int, words[:2])
            record.intensities[x,y]  = float(words[2])
            record.stdevs[x,y] = float(words[3])
            record.npix[x,y]  = int(words[4])
    return record
