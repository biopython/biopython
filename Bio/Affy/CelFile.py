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

# We use print in the doctests
from __future__ import print_function

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.Affy.CelFile")

__docformat__ = "restructuredtext en"


class Record(object):
    """Stores the information in a cel file

    Example usage:

    >>> from Bio.Affy import CelFile
    >>> with open('Affy/affy_v3_example.CEL') as handle:
    ...     c = CelFile.read(handle)
    ...
    >>> print(c.ncols, c.nrows)
    5 5
    >>> print(c.intensities)
    [[   234.    170.  22177.    164.  22104.]
     [   188.    188.  21871.    168.  21883.]
     [   188.    193.  21455.    198.  21300.]
     [   188.    182.  21438.    188.  20945.]
     [   193.  20370.    174.  20605.    168.]]
    >>> print(c.stdevs)
    [[   24.     34.5  2669.     19.7  3661.2]
     [   29.8    29.8  2795.9    67.9  2792.4]
     [   29.8    88.7  2976.5    62.   2914.5]
     [   29.8    76.2  2759.5    49.2  2762. ]
     [   38.8  2611.8    26.6  2810.7    24.1]]
    >>> print(c.npix)
    [[25 25 25 25 25]
     [25 25 25 25 25]
     [25 25 25 25 25]
     [25 25 25 25 25]
     [25 25 25 25 25]]

    """
    def __init__(self):
        self.version = None
        self.GridCornerUL = None
        self.GridCornerUR = None
        self.GridCornerLR = None
        self.GridCornerLL = None
        self.DatHeader = None
        self.Algorithm = None
        self.AlgorithmParameters = None
        self.NumberCells = None
        self.intensities = None
        self.stdevs = None
        self.npix = None
        self.nrows = None
        self.ncols = None
        self.nmask = None
        self.mask = None
        self.noutliers = None
        self.outliers = None
        self.modified = None


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
        # Set current section
        if line[:5] == "[CEL]":
            section = "CEL"
        elif line[:8] == "[HEADER]":
            section = "HEADER"
        elif line[:11] == "[INTENSITY]":
            section = "INTENSITY"
            record.intensities = numpy.zeros((record.nrows, record.ncols))
            record.stdevs = numpy.zeros((record.nrows, record.ncols))
            record.npix = numpy.zeros((record.nrows, record.ncols), int)
        elif line[:7] == "[MASKS]":
            section = "MASKS"
            record.mask = numpy.zeros((record.nrows, record.ncols))
        elif line[:10] == "[OUTLIERS]":
            section = "OUTLIERS"
            record.outliers = numpy.zeros((record.nrows, record.ncols))
        elif line[:10] == "[MODIFIED]":
            section = "MODIFIED"
            record.modified = numpy.zeros((record.nrows, record.ncols))
        elif line[0] == "[":
            # This would be an unknown section
            section = ""
        elif section == "CEL":
            keyword, value = line.split("=", 1)
            if keyword == 'Version':
                record.version = int(value)
        elif section == "HEADER":
            # Set record.ncols and record.nrows, remaining data goes into
            # record.header dict
            keyword, value = line.split("=", 1)
            if keyword == "Cols":
                record.ncols = int(value)
            elif keyword == "Rows":
                record.nrows = int(value)
            elif keyword == 'GridCornerUL':
                x, y = value.split()
                record.GridCornerUL = (int(x), int(y))
            elif keyword == 'GridCornerUR':
                x, y = value.split()
                record.GridCornerUR = (int(x), int(y))
            elif keyword == 'GridCornerLR':
                x, y = value.split()
                record.GridCornerLR = (int(x), int(y))
            elif keyword == 'GridCornerLL':
                x, y = value.split()
                record.GridCornerLL = (int(x), int(y))
            elif keyword == 'DatHeader':
                record.DatHeader = value.strip('\n\r')
            elif keyword == 'Algorithm':
                record.Algorithm = value.strip('\n\r')
            elif keyword == 'AlgorithmParameters':
                record.AlgorithmParameters = value.strip('\n\r')
        elif section == "INTENSITY":
            if "NumberCells" in line:
                record.NumberCells = int(line.split("=", 1)[1])
            elif "CellHeader" in line:
                pass
            else:
                words = line.split()
                y = int(words[0])
                x = int(words[1])
                record.intensities[x, y] = float(words[2])
                record.stdevs[x, y] = float(words[3])
                record.npix[x, y] = int(words[4])
        elif section == "MASKS":
            if "NumberCells" in line:
                record.nmask = int(line.split("=", 1)[1])
            elif "CellHeader" in line:
                pass
            else:
                words = line.split()
                y = int(words[0])
                x = int(words[1])
                record.mask[x, y] = int(1)
        elif section == "OUTLIERS":
            if "NumberCells" in line:
                record.noutliers = int(line.split("=", 1)[1])
            elif "CellHeader" in line:
                pass
            else:
                words = line.split()
                y = int(words[0])
                x = int(words[1])
                record.outliers[x, y] = int(1)
        elif section == "MODIFIED":
            if "NumberCells" in line:
                record.nmodified = int(line.split("=", 1)[1])
            elif "CellHeader" in line:
                pass
            else:
                words = line.split()
                y = int(words[0])
                x = int(words[1])
                record.modified[x, y] = float(words[2])
        else:
            continue
    return record

if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()
