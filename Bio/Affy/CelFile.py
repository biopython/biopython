# Copyright 2004 by Harry Zuzan.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
No version number yet.

Classes for accessing the information in Affymetrix cel files.

class CelParser: parses cel files
class CelRecord: stores the information from a cel file

"""

import _cel

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
    Parses an Affymetrix cel file passed in as a string and returns
    an instance of a CelRecord

    This class needs error handling.
    """

    def __init__(self, data=None):
        """
        Usually load the class with the cel file (not file name) as
        an argument.
        """
        
        self._intensities = None
        self._stdevs      = None
        self._npix        = None

        if data is not None: self.parse(data)


    def parse(self, data):
        """
        Takes the contents of a cel file passed as a string, parses it
        and stores it in the three arrays.

        There is more information in the cel file that could be retrieved
        and stored in CelRecord.  The chip type should be a priority.
        """

        (self._intensities, self._stdevs, self._npix) = _cel.parse(data)
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

