# Copyright 2001 by Gavin E. Crooks.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This functionality may be of general use, in which case this module should
# be moved out of the SCOP package.

import warnings
warnings.warn("Bio.SCOP.FileIndex was deprecated, as it does not seem to have any users. If you do use this module, please contact the Biopython developers at biopython-dev@biopython.org to avoid permanent removal of this module", DeprecationWarning)



class defaultdict(dict):

    def __init__(self, default=None):
        dict.__init__(self)
        self.default = default

    def __getitem__(self, key):
        try:
            return dict.__getitem__(self, key)
        except KeyError:
            return self.default

class FileIndex(dict):
    """ An in memory index that allows rapid random access into a file.

    The class can be used to turn a file into a read-only
    database.
    """
    def __init__(self, filename, iterator_gen, key_gen ):
        """
        Arguments:
        
          filename  -- The file to index

          iterator_gen --   A function that eats a file handle, and returns
            a file iterator. The iterator has a method next()
            that returns the next item to be indexed from the file.

          key_gen -- A function that generates an index key from the items
                     created by the iterator. 
        """
        dict.__init__(self)
        
        self.filename = filename
        self.iterator_gen = iterator_gen
        
        f = open(self.filename)
        try:
            loc = 0
            i = self.iterator_gen(f)
            while 1:
                next_thing = i.next()
                if next_thing is None : break
                key = key_gen(next_thing)
                if key != None:
                    self[key]=loc
                loc = f.tell()
        finally:
            f.close()

    def __getitem__(self, key):
        """ Return an item from the indexed file. """
        loc = dict.__getitem__(self,key)

        f = open(self.filename)
        try:
            f.seek(loc)
            thing = self.iterator_gen(f).next()
        finally:
            f.close()
        return thing
