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


The following classes are obsolete:

class CelParser: parses cel files
class CelRecord: stores the information from a cel file

"""

import numpy

class Record:
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


# Everything below is considered obsolete

from Bio.ParserSupport import AbstractConsumer
from numpy import *

class CelScanner:
    """Scannner for Affymetrix CEL files.

    Methods:
    feed     Feed data into the scanner.

    The scanner generates (and calls the consumer) the following
    types of events:

    Rows - the number of rows on the microarray
    Cols - the number of columns on the microarray
    StartIntensity - generated when the section [INTENSITY] is found
    ReadIntensity - one line in the section [INTENSITY]

    """
    def feed(self, handle, consumer):
        """scanner.feed(handle, consumer)

        Feed in a handle to a Cel file for scanning.  handle is a file-like
        object that contains the Cel file.  consumer is a Consumer
        object that will receive events as the report is scanned.
        """
        section = ""
        for line in handle:
            if line.strip()=="": continue
            if line[0]=="[":
                section = ""
                if line[:8]=="[HEADER]":
                    section = "HEADER"
                elif line[:11]=="[INTENSITY]":
                    section = "INTENSITY"
                    consumer.StartIntensity()
                continue
            if section=="HEADER":
                keyword, value = line.split("=", 1)
                if keyword=="Cols": consumer.Cols(value)
                if keyword=="Rows": consumer.Rows(value)
                continue
            elif section=="INTENSITY":
                if "=" in line: continue
                consumer.ReadIntensity(line)


class CelConsumer(AbstractConsumer):

    def __init__(self):
        self._mean  = None
        self._stdev = None
        self._npix  = None

    def Cols(self, value):
        self._cols = int(value)

    def Rows(self, value):
        self._rows = int(value)

    def StartIntensity(self):
        self._mean  = zeros((self._rows, self._cols))
        self._stdev = zeros((self._rows, self._cols))
        self._npix  = zeros((self._rows, self._cols), int)

    def ReadIntensity(self, line):
        y, x, mean, stdev, npix = map(float, line.split())
        x = int(x)
        y = int(y)
        self._mean[x,y]  = mean
        self._stdev[x,y] = stdev
        self._npix[x,y]  = int(npix)

class CelRecord:
    """
    Stores the information in a cel file

    Needs error handling.
    Needs to know the chip design.
    """


    def __init__(self, data_dict):
        """
        Pass the data attributes as a dictionary.
        """
        from copy import deepcopy as dcopy

        self._intensities = dcopy(data_dict['intensities'])
        self._stdevs      = dcopy(data_dict['stdevs'])
        self._npix        = dcopy(data_dict['npix'])

        self._nrows, self._ncols = self._intensities.shape


    def intensities(self):
        """
        Return a two dimensional array of probe cell intensities.
        Dimension 1 -> rows
        Dimension 2 -> columns
        """
        return self._intensities


    def stdevs(self):
        """
        Return a two dimensional array of probe cell standard deviations.
        Dimension 1 -> rows
        Dimension 2 -> columns
        """
        return self._stdevs


    def npix(self):
        """
        Return a two dimensional array of the number of pixels in a probe cell.
        Dimension 1 -> rows
        Dimension 2 -> columns
        """
        return self._npix


    def nrows(self):
        """
        The number of rows of probe cells in an array.
        """
        return self._nrows

    def ncols(self):
        """
        The number of columns of probe cells in an array.
        """
        return self._ncols

    def size(self):
        """
        The size of the probe cell array as a tuple (nrows,ncols).
        """
        return self._nrows, self._ncols



class CelParser:
    """
    Takes a handle to an Affymetrix cel file, parses the file and
    returns an instance of a CelRecord

    This class needs error handling.
    """

    def __init__(self, handle=None):
        """
        Usually load the class with the cel file (not file name) as
        an argument.
        """
        
        self._intensities = None
        self._stdevs      = None
        self._npix        = None

        if handle is not None: self.parse(handle)


    def parse(self, handle):
        """
        Takes a handle to a cel file, parses it
        and stores it in the three arrays.

        There is more information in the cel file that could be retrieved
        and stored in CelRecord.  The chip type should be a priority.
        """

        # (self._intensities, self._stdevs, self._npix) = _cel.parse(data)
        scanner = CelScanner()
        consumer = CelConsumer()
        scanner.feed(handle, consumer)
        self._intensities = consumer._mean
        self._stdevs = consumer._stdev
        self._npix = consumer._npix
        self._nrows = self._intensities.shape[0]
        self._ncols = self._intensities.shape[1]


    def __call__(self):
        """
        Returns the parsed data as a CelRecord.
        """

        record_dict = {}
        record_dict['intensities'] = self._intensities
        record_dict['stdevs'] = self._stdevs
        record_dict['npix'] = self._npix

        return CelRecord(record_dict)

