# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Index.py

This module provides a way to create indexes to text files.

"""

# Using anydbm or shelve, pickled or not pickled
# fasta file with 82000 sequences, 37Mb
# index 'dat' file is 42Mb and 'dir' file is 8Mb

import os
import shelve
import UserDict

class Index(UserDict.UserDict):
    """Create 

    """
    
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
            self.data[Index.__version_key] = Index.__version
        else:
            # Check to make sure the database is the correct version.
            version = self.data.get(Index.__version_key, None)
            if version is None:
                raise IOError, "Unrecognized index format"
            elif version != Index.__version:
                raise IOError, "Version %d doesn't match my version %s" % \
                      (version, Index.__version)
            
    def __del__(self):
        if self.__dict__.has_key('data'):
            self.data.close()
