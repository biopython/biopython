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

    # TODO: Improve Docstrings
    # TODO: Check for self.filter()? Or implement this in SearchIO.parse?

    """Class representing search results from a single query.

    The Result object represents a search result from a single query, with
    the following behaviors:

    * You can get the corresponding Hit object by using the Hit ID as key
    * Hit ID keys must be a string. This limitation was put to enable
      indexing and/or slicing of the Result object.
    * Iteration over Result returns the Hit objects within.

    """

    def __init__(self, program, target, query_id, meta={}, *args, **kwargs):
        """Initializes a Result object.

        program -- String of search program name.
        query_id -- String of query sequence ID.
        target -- String of database name to search against.
        meta -- Dictionary of additional information about the search.
                  This is the information stored in the header of the
                  search output file (anything prior to the first result,
                  if it exists) and varies depending on the search program
                  used.

        """
        # meta must be a dict
        if not isinstance(meta, dict):
            raise TypeError("Meta argument must be a dictionary object.")

        self.program = program
        self.id = query_id
        self.target = target
        self.meta = meta

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
        # TODO: is there a way to do this before self._hits creation that's
        # compatible with OrderedDict creation?
        for hit in self.hits:
            self._validate_hit(hit)

    def __repr__(self):
        return "Result(program='%s', target='%s', id='%s', %i hits)" % \
                (self.program, self.target, self.id, len(self._hits))

    # handle Python 2 OrderedDict behavior
    if hasattr(OrderedDict, 'iteritems'):

        def __iter__(self):
            return iter(self.iterhits())

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

    else:

        def __iter__(self):
            return iter(self.hits)

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

    def __contains__(self, hit_key):
        """Checks whether a Hit object or a Hit object with the given ID exists."""
        if isinstance(hit_key, Hit):
            return hit_key.id in self._hits
        return hit_key in self._hits

    def __len__(self):
        return len(self._hits)

    def __nonzero__(self):
        return bool(self._hits)

    def __reversed__(self):
        items = reversed(list(self._hits.items()))
        return self.__class__(self.program, self.target, self.id, \
                self.meta, items)

    def __setitem__(self, hit_key, hit):
        """Custom Search __setitem__.

        Hit key must be a string and hit must be a Hit object.

        """
        if not isinstance(hit_key, basestring):
            raise TypeError("Result object keys must be a string.")
        self._validate_hit(hit)

        self._hits[hit_key] = hit

    def __getitem__(self, hit_key):
        """Custom Search __getitem__.

        Allows value retrieval by its key, location index, or a slice of
        location index. Also allows direct HSP retrieval if hit key is a tuple
        or list of hit ID and HSP index.

        """
        # allow for direct HSP retrieval if key is list / tuple
        # of two items (hit index / id, hsp rank)
        # retrieve Hit object first, then get the HSP according
        # to the given index
        if isinstance(hit_key, (tuple, list)):
            # ensure the tuple / list has at least 2 items
            try:
                hit_id = hit_key[0]
                # note that this allows the second item to be a slice as well
                hsp_rank = hit_key[1]
            except IndexError:
                raise ValueError("Expected tuple or list with 2 items, found: "
                        "%i item(s)." % len(hit_key))

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
        elif isinstance(hit_key, slice):
            # should we return just a list of Hits instead of a full blown
            # Result object if it's a slice?
            items = list(self._hits.items())[hit_key]
            return self.__class__(self.program, self.target, self.id, \
                    self.meta, items)

        # if key is an int, then retrieve the Hit at the int index
        elif isinstance(hit_key, int):
            return list(self.hits)[hit_key]

        # if key is a string, then do a regular dictionary retrieval
        return self._hits[hit_key]

    def __delitem__(self, hit_key):
        """Custom Search __delitem__

        If key is a string, then the method will delete the Hit object whose
        ID matches the string. If key is an integer or a slice object, then
        the Hit objects within that range will be deleted.

        Also accepts a list of strings, which will delete all the Hit whose ID
        matches the strings in the list.

        """
        # if it's an integer or slice, get the corresponding key first
        # and put it into a list
        if isinstance(hit_key, int):
            hit_keys = [list(self.hit_ids)[hit_key]]
        # the same, if it's a slice
        elif isinstance(hit_key, slice):
            hit_keys = list(self.hit_ids)[hit_key]
        # if it's already a list or tuple, we just reassign it to a new name
        elif isinstance(hit_key, (list, tuple)):
            hit_keys = hit_key
        # otherwise put it in a list
        else:
            hit_keys = [hit_key]

        for key in hit_keys:
            del self._hits[key]

    def _validate_hit(self, hit):
        """Checks whether the Hit object is of the correct type and has the right query ID."""
        if not isinstance(hit, Hit):
            raise TypeError("Result objects can only contain Hit objects.")
        if hit.query_id != self.id:
            raise ValueError("Expected Hit with query ID '%s', found '%s' "
                    "instead." % (self.id, hit.query_id))

    # marker for default self.pop() return value
    # this method is adapted from Python's built in OrderedDict.pop
    # implementation
    __marker = object()

    def pop(self, hit_key=-1, default=__marker):
        """Removes the specified key and return its value.

        Similar to __getitem__, pop also allows Hit retrieval (removal in this
        case) by its index in addition to its ID.

        If the key does not exist and default is given, return the default
        value. Otherwise a KeyError will be raised.

        """
        # if key is an integer (index)
        # get the ID for the Hit object at that index
        if isinstance(hit_key, int):
            # raise the appropriate error if there is no hit
            if not self:
                raise IndexError("pop from empty list")
            hit_key = list(self.hit_ids)[hit_key]

        try:
            return self._hits.pop(hit_key)
        except KeyError:
            # if key doesn't exist and no default is set, raise a KeyError
            if default is self.__marker:
                raise KeyError(hit_key)
        # if key doesn't exist but a default is set, return the default value
        return default

    def rank(self, hit_key):
        """Returns the rank of a given Hit ID, 0-based.

        Also accepts a Hit object as the argument, which returns the rank of
        the Hit object ID. If the given key is not found, returns -1 instead.

        """
        try:
            if isinstance(hit_key, Hit):
                return list(self.hit_ids).index(hit_key.id)
            return list(self.hit_ids).index(hit_key)
        except ValueError:
            return -1

    def sort(self, key=None, reverse=False):
        # no cmp argument to make sort more Python 3-like
        """Sorts the Hit objects.

        key -- Function used to sort the Hit objects.
        reverse -- Boolean, whether to reverse the sorting or not.

        By default, sorting is based on the expect values of the Hit objects,
        from the smallest to the largest. If the Hit objects do not have any
        expect values (e.g. BLAT Hit objects), then no sorting is performed.

        The sort creates a new Hit container object, but appears to be
        in-place since the new Hit container replaces the old one.

        """
        # if no key is specified, attempt to sort by Hit evalue
        if key is None:
            try:
                sorted_hits = OrderedDict(sorted(self.items, \
                        key=lambda hit: hit.evalue, reverse=reverse))
            # handle cases where the Hit objects doesn't have evalue
            except AttributeError:
                # if reverse is set to True, reverse the ordering
                # we don't use sorted() since no __eq__ etc. magic methods
                # are defined for Hit objects
                if reverse:
                    sorted_hits = OrderedDict(self.items[::-1])
                # otherwise the object is the same as the old one
                else:
                    sorted_hits = self._hits
        # otherwise try to sort using the given parameters
        # and let any exceptions rise to the top
        else:
            sorted_hits = OrderedDict(sorted(self.items, key=key, \
                    reverse=reverse))

        self._hits = sorted_hits


class Hit(object):

    """Class representing the entire database entry of a sequence match.

    """

    def __init__(self, hit_id, query_id, hsps=[]):
        """Initializes a Hit object.

        query_id -- String of the query name used to obtain this hit.
        hit_id -- String of unique identifier for this hit.
        hsps -- List containing HSP objects.

        """
        self.query_id= query_id
        self.id = hit_id

        if not hsps:
            raise ValueError("Hit object must contain at least one HSP object.")

        self._hsps = []
        for hsp in hsps:
            self.append(hsp)

    def __repr__(self):
        return "Hit(id='%s', %i alignments)" % (self.id, len(self))

    @property
    def hsps(self):
        return self._hsps

    def __iter__(self):
        return iter(self._hsps)

    def __len__(self):
        return len(self._hsps)

    def __nonzero__(self):
        return bool(self._hsps)

    def __reversed__(self):
        return self.__class__(self.id, self.query_id, reversed(self._hsps))

    def __setitem__(self, idx, hsps):
        self._validate_hsps(hsps)
        self._hsps[idx] = hsps

    def __getitem__(self, idx):
        # if key is slice, return a new Hit instance
        if isinstance(idx, slice):
            return self.__class__(self.id, self.query_id, self._hsps[idx])
        return self._hsps[idx]

    def __delitem__(self, idx):
        del self._hsps[idx]

    def _validate_hsp(self, hsp):
        if not isinstance(hsp, HSP):
            raise TypeError("Hit objects can only contain HSP objects.")
        if hsp.hit_id != self.id:
            raise ValueError("Only HSP objects from the same Hit can be added.")

    def append(self, hsp):
        self._validate_hsp(hsp)
        self._hsps.append(hsp)

    def pop(self, index=-1):
        return self._hsps.pop(index)

    def reverse(self):
        self._hsps.reverse()

    def sort(self, key=None, reverse=False):
        self._hsps.sort(key=key, reverse=reverse)


class HSP(object):

    """Class representing high-scoring alignment regions of the query and hit.

    """


class SearchIndexer(object):

    """Iterator that returns file positions of results in a Search output file.

    """
