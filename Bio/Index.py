# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Index.py

This module provides a way to create indexes to text files.

Classes:
Index     Dictionary-like class used to store index information.

"""

import os
import shelve
import string
import UserDict

class _ShelveIndex(UserDict.UserDict):
    """An index file wrapped around shelve.

    """
    # Without a good dbm module installed, this is pretty slow and
    # generates large files.  When generating an index on a FASTA-
    # formatted file with 82000 sequences (37Mb), the 
    # index 'dat' file is 42Mb and 'dir' file is 8Mb.
    __version = 2
    __version_key = '__version'

    def __init__(self, indexname, truncate=None):
        UserDict.UserDict.__init__(self)
        try:
            if truncate:
                # In python 1.52 and before, dumbdbm (under shelve)
                # doesn't clear the old database.
                files = [indexname + '.dir',
                         indexname + '.dat',
                         indexname + '.bak'
                         ]
                for file in files:
                    if os.path.exists(file):
                        os.unlink(file)
                raise "open a new shelf"
            self.data = shelve.open(indexname, flag='r')
        except:
            # No database exists.
            self.data = shelve.open(indexname, flag='n')
            self.data[self.__version_key] = self.__version
        else:
            # Check to make sure the database is the correct version.
            version = self.data.get(self.__version_key, None)
            if version is None:
                raise IOError, "Unrecognized index format"
            elif version != self.__version:
                raise IOError, "Version %s doesn't match my version %s" % \
                      (version, self.__version)
            
    def __del__(self):
        if self.__dict__.has_key('data'):
            self.data.close()

class _CustomIndex(UserDict.UserDict):
    """This creates an in-memory index file.

    # Format:
    version
    key value

    """
    __version = 3
    __version_key = '__version'

    def __init__(self, indexname, truncate=None):
        raise NotImplementedError
        UserDict.UserDict.__init__(self)
        if truncate and os.path.exists(indexname):
            os.unlink(indexname)

        if os.path.exists(indexname):
            handle = open(indexname)
            version = string.rstrip(handle.readline())
            if version != self.__version:
                raise IOError, "Version %s doesn't match my version %s" % \
                      (version, self.__version)
            lines = handle.readlines()
            lines = map(string.split, lines)
            for key, value in lines:
                pass
            
    def __del__(self):
        # XXX save file
        if self.__dict__.has_key('data'):
            self.data.close()

Index = _ShelveIndex
