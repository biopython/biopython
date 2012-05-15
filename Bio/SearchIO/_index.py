# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Custom indexing for Bio.SearchIO objects (PRIVATE).

"""

from sqlite3 import dbapi2 as sqlite
from sqlite3 import IntegrityError, OperationalError

try:
    from collections import UserDict as _dict_base
except ImportError:
    from UserDict import DictMixin as _dict_base


class IndexedSearch(_dict_base):

    """Dictionary-like object for implementing Search indexing.

    """
    def __init__(self, handle, format, indexer, key_function):
        """Initializes IndexedSearch instance.

        handle -- The source filename as string.
        format -- Lower case string denoting one of the supported formats.
        indexer -- Format-specific Indexer class that.
        key_function -- Optional callbak function which when given a Result
                        should return a unique key for the dictionary.

        """
        self._source = handle
        self._format = format
        self._key_function = key_function

        indexed_obj = indexer(handle)
        self._indexer = indexed_obj

        # default key function is lambda rec: rec.id
        offset_iter = ((key_function(key), offset, length) for \
                (key, offset, length) in indexed_obj)

        index = {}
        for key, offset, length in offset_iter:
            if key in index:
                self._indexer._handle.close()
                raise ValueError("Duplicate key '%s'" % key)
            else:
                index[key] = offset

        self._index = index

    def __repr__(self):
        return "IndexedSearch('%r', '%r', key_function=%r)" % \
                (self._source, self._format, self._key_function)

    def __str__(self):
        if self:
            return "{%s: Result(...), ...}" % repr(self.keys()[0])
        else:
            return "{}"

    def __contains__(self, key):
        return key in self._index

    def __len__(self):
        return len(self._index)

    # dynamic method assignments to deal with different dict behavior in
    # python 2 and python 3
    # case 1: python2
    if hasattr(dict, 'iteritems'):
        def values(self):
            raise NotImplementedError("Due to memory concerns, when indexing "
                    "a search output file you cannot access all the records "
                    "at once.")

        def items(self):
            raise NotImplementedError("Due to memory concerns, when indexing "
                    "a search output file you cannot access all the records "
                    "at once.")

        def keys(self):
            return self._index.keys()

        def itervalues(self):
            for key in self.__iter__():
                yield self.__getitem__(key)

        def iteritems(self):
            for key in self.__iter__():
                yield key, self.__getitem__(key)

        def iterkeys(self):
            return self.__iter__()

    # case 2: python3, where dict.{items,keys,values} returns a generator-like
    # object
    else:
        def values(self):
            for key in self.__iter__():
                yield self.__getitem__(key)

        def items(self):
            for key in self.__iter__():
                yield key, self.__getitem__(key)

        def keys(self):
            return self.__iter__()

    def __iter__(self):
        return iter(self._index)

    def __getitem__(self, key):
        result = self._indexer.get(self._index[key])
        key2 = self._key_function(result.id)

        if key != key2:
            raise ValueError("Key did not match (%s vs %s)" % (key, key2))

        return result

    def get(self, key, default=None):
        try:
            return self.__getitem__(key)
        except KeyError:
            return default

    def get_raw(self, key):
        """Similar to get, but returns the raw string of the Result object.

        Note that any text or information in the search output file located
        prior to the first result is ignored.

        """
        return self._indexer.get_raw(self._index[key])

    def __delitem__(self, key, value):
        raise NotImplementedError("An indexed search output file is read-only.")

    def __setitem__(self, key, value):
        raise NotImplementedError("An indexed search output file is read-only.")

    def clear(self, *args, **kwargs):
        raise NotImplementedError("An indexed search output file is read-only.")

    def pop(self, key, default=None):
        raise NotImplementedError("An indexed search output file is read-only.")

    def popitem(self):
        raise NotImplementedError("An indexed search output file is read-only.")

    def update(self, *args, **kwargs):
        raise NotImplementedError("An indexed search output file is read-only.")

    def copy(self, *args, **kwargs):
        raise NotImplementedError("An indexed search output file does not "
                "support this.")

    def fromkeys(self, *args, **kwargs):
        raise NotImplementedError("An indexed search output file does not "
                "support this.")
