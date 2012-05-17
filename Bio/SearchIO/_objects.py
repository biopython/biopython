# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO objects to model homology search program outputs (PRIVATE).

"""

try:
    from collections import OrderedDict
except ImportError:
    from Bio.SearchIO._aux import OrderedDict

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment


class Result(object):

    """Class representing search results from a single query.

    The Result object represents a search result from a single query. It
    is a subclass of Python's built in OrderedDict class, so it
    behaves quite similar to OrderedDict, with the following caveats:

    * You can get the corresponding Hit object by using the Hit ID as key
    * Hit ID keys must be a string. This limitation was put to enable
      indexing and/or slicing of the Result object.
    * Iteration over Result returns the Hit objects contained within.

    """

    def __init__(self, program, query, target, header={}, *args, **kwargs):
        """Initializes a Result object.

        program -- String of search program name.
        query -- String of query sequence ID.
        target -- String of database name to search against.
        header -- Dictionary of additional information about the search.
                  This is the information stored in the header of the
                  search output file (anything prior to the first result,
                  if it exists) and varies depending on the search program
                  used.

        """
        # header must be a dict
        if not isinstance(header, dict):
            raise TypeError("Header argument must be a dictionary object.")

        self.program = program
        self.query = query
        self.target = target
        self.header = header

        # We're implementing Result as a wrapper for OrderedDict;
        # it could be implemented as a subclass of OrderedDict, but iterating
        # over Result would then return the dict keys instead of the Hit objects
        # contained within. I haven't found an elegant way to make it return
        # the Hit objects without breaking the other methods (e.g. keys())
        # The tradeoff of implementing it as and OrderedDict wrapper is that
        # there's more methods to define, so the code is more verbose. However
        # this also makes it more malleable, which could be useful in the future.
        self._hits = OrderedDict(*args, **kwargs)

        # check if there's any non-Hit objects
        # TODO: is there a way to do this before self._hits creation?
        if not all([isinstance(x, Hit) for x in self.hits]):
            raise TypeError("Result object values must be a Hit object.")

    def __repr__(self):
        return "Result(program='%s', query='%s', target='%s', %i hits)" % \
                (self.program, self.query, self.target, len(self._hits))

    # handle Python 2 OrderedDict behavior
    if hasattr(OrderedDict, 'iteritems'):

        @property
        def hits(self):
            """Returns a list of Hit objects contained by this object."""
            return self._hits.values()

        @property
        def hit_ids(self):
            """Returns a list of Hit IDs contained by this object."""
            return self._hits.keys()

        @property
        def items(self):
            """Returns a tuple of Hit ID and Hit object contained by this object."""
            return self._hits.items()

        def iterhits(self):
            """Returns an iterator over the Hit objects."""
            for hit in self._hits.itervalues():
                yield hit

        def iterhit_ids(self):
            """Returns an iterator over the ID of the Hit objects."""
            for hit_id in self._hits.iterkeys():
                yield hit_id

        def iteritems(self):
            """Returns an iterator of a tuple of Hit ID and Hit objects."""
            for item in self._hits.iteritems():
                yield item

        def __iter__(self):
            return iter(self.iterhits)

    else:

        @property
        def hits(self):
            """Returns an iterator over the Hit objects contained by this object."""
            for hit in self._hits.values():
                yield hit

        @property
        def hit_ids(self):
            """Returns an iterator over the Hit IDs contained by this object."""
            for hit_id in  self._hits.keys():
                yield hit_id

        @property
        def items(self):
            """Returns an iterator over the Hit ID and Hit object contained by this object."""
            for item in self._hits.items():
                yield item

        def __iter__(self):
            return iter(self.hits)

    def rank(self, key):
        """Returns the rank of a given Hit ID, 0-based.

        Also accepts a Hit object as the argument, which returns the rank of
        the Hit object ID. If the given key is not found, returns -1 instead.

        """
        try:
            if isinstance(key, Hit):
                return list(self.hit_ids).index(key.id)
            return list(self.hit_ids).index(key)
        except ValueError:
            return -1

    def sort(self, cmp=None, key=lambda hit: hit.evalue, reverse=False):
        """Sorts the Hit objects.

        The sort creates a new Hit container object, but appears to be
        in-place since the new Hit container replaces the old one.

        By default, sorting is based on the expect values of the Hit objects,
        from the smallest to the largest. If the Hit objects does not have any
        expect values (e.g. BLAT Hit objects), then sorting is based on the
        Hit IDs.

        """
        # create the new sorted OrderedDict
        try:
            sorted_hits = OrderedDict(sorted(self.items, cmp, key, reverse))
        except AttributeError:
            key = lambda hit: hit.id
            sorted_hits = OrderedDict(sorted(self.items, cmp, key, reverse))

        # and replace the old one
        self._hits = sorted_hits

    # marker for default self.pop() return value
    # this method is adapted from Python's built in OrderedDict.pop
    # implementation
    __marker = object()

    def pop(self, key, default=__marker):
        """Removes the specified key and return its value.

        Similar to __getitem__, pop also allows Hit retrieval (removal in this
        case) by its index in addition to its ID.

        If the key does not exist and default is given, return the default
        value. Otherwise a KeyError will be raised.

        """
        # if key is an integer (index)
        # get the ID for the Hit object at that index
        if isinstance(key, int):
            key = list(self.hit_ids)[key]

        try:
            return self._hits.pop(key)
        except KeyError:
            # if key doesn't exist and no default is set, raise a KeyError
            if default is self.__marker:
                raise KeyError(key)
        # if key doesn't exist but a default is set, return the default value
        return default

    def __contains__(self, key):
        """Checks whether a Hit object or a Hit object with the given ID exists."""
        if isinstance(key, Hit):
            return key.id in self._hits
        return key in self._hits

    def __len__(self):
        return len(self._hits)

    def __nonzero__(self):
        return bool(self._hits)

    def __reversed__(self):
        items = reversed(list(self._hits.items()))
        return self.__class__(self.program, self.query, self.target, \
                self.header, items)

    def __setitem__(self, key, value):
        """Custom Search __setitem__.

        Key must be a string and value must be a Hit object.

        """
        if not isinstance(key, basestring):
            raise TypeError("Result object keys must be a string.")
        if not isinstance(value, Hit):
            raise TypeError("Result object values must be a Hit object.")

        self._hits[key] = value

    def __delitem__(self, key):
        """Custom Search __delitem__

        If key is a string, then the method will delete the Hit object whose
        ID matches the string. If key is an integer or a slice object, then
        the Hit objects within that range will be deleted.

        Also accepts a list of strings, which will delete all the Hit whose ID
        matches the strings in the list.

        """
        # if it's an integer or slice, get the corresponding key first
        # and put it into a list
        if isinstance(key, int):
            keys = [list(self.hit_ids)[key]]
        # the same, if it's a slice
        elif isinstance(key, slice):
            keys = list(self.hit_ids)[key]
        # if it's already a list or tuple, we just reassign it to a new name
        elif isinstance(key, (list, tuple)):
            keys = key
        # otherwise put it in a list
        else:
            keys = [key]

        for k in keys:
            del self._hits[k]

    def __getitem__(self, key):
        """Custom Search __getitem__.

        Allows value retrieval by its key, location index, or a slice of
        location index. Also allows direct HSP retrieval if key is a tuple
        or list of hit ID and HSP index.

        """
        # allow for direct HSP retrieval if key is list / tuple
        # of two items (hit index / id, hsp rank)
        # retrieve Hit object first, then get the HSP according
        # to the given index
        if isinstance(key, (tuple, list)):
            # ensure the tuple / list has at least 2 items
            try:
                hit_id = key[0]
                # note that this allows the second item to be a slice as well
                hsp_rank = key[1]
            except IndexError:
                raise ValueError("Expected tuple or list with 2 items, found: "
                        "%i item(s)." % len(key))

            # allow for hit retrieval by string ID or integer index
            if isinstance(hit_id, basestring):
                # if key is string, we can retrieve directly from the dict
                hit = self._hits[hit_id]
            elif isinstance(hit_id, int):
                # if not, retrieve from the list
                hit = list(self.hits)[hit_id]
            else:
                raise TypeError("Hit index must be a string or an integer if "
                        "double-index slicing is performed.")

            return hit[hsp_rank]

        # retrieval using slice objects returns another Result object
        elif isinstance(key, slice):
            # should we return just a list of Hits instead of a full blown
            # Result object if it's a slice?
            items = list(self._hits.items())[key]
            return self.__class__(self.program, self.query, self.target, \
                    self.header, items)

        # if key is an int, then retrieve the Hit at the int index
        elif isinstance(key, int):
            return list(self.hits)[key]

        # if key is a string, then do a regular dictionary retrieval
        return self._hits[key]


class Hit(object):

    """Class representing the entire database entry of a sequence match.

    """


class HSP(object):

    """Class representing high-scoring alignment regions of the query and hit.

    """


class SearchIndexer(object):

    """Iterator that returns file positions of results in a Search output file.

    """
