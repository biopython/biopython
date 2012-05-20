# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO objects to model homology search program outputs (PRIVATE).

"""

from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import single_letter_alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio._py3k import OrderedDict


class _StickyObject(object):

    """Abstract class that defines a method to transfer instance attributes."""

    _NON_STICKY_ATTRS = ()

    def _transfer_attrs(self, obj):
        """Transfer instance attributes to the given object.

        This method is used to transfer attributes set externally (for example
        using setattr()) to a new object created from this one (for example
        from slicing).

        The reason this method is necessary is because different parsers will
        set different attributes for each Result, Hit, or HSP object they use,
        depending on the attributes they found in the search output file.
        Ideally, we want these attributes to 'stick' with any new instance
        object created from the original one.

        """
        # list of attribute names we don't want to transfer
        for attr in self.__dict__.keys():
            if attr not in self._NON_STICKY_ATTRS:
                setattr(obj, attr, self.__dict__[attr])


class Result(_StickyObject):

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

    # attributes we don't want to transfer when creating a new Result class
    # from this one
    _NON_STICKY_ATTRS = ('_hits',)

    def __init__(self, query_id, hits=[], meta={}, program='<unknown>', \
            target='<unknown>', hit_key_function=lambda hit: hit.id):
        """Initializes a Result object.

        query_id -- String of query sequence ID.
        hits -- List of Hit objects.
        meta -- Dictionary of additional information about the search. This is
                the information stored in the header of the search output file
                (anything prior to the first result, if it exists) and varies
                depending on the search program used.
        program -- String of search program name.
        target -- String of database name to search against.
        hit_key_function -- Function to define hit keys.

        """
        # meta must be a dict
        if not isinstance(meta, dict):
            raise TypeError("Meta argument must be a dictionary object.")

        self.id = query_id
        self.meta = meta
        self.program = program
        self.target = target
        self._hit_key_function = hit_key_function
        self._hits = OrderedDict()

        # validate Hit objects and fill up self._hits
        for hit in hits:
            # validation is handled by __setitem__
            self.append(hit)

    def __repr__(self):
        hit = 'hit'
        if len(self) != 1:
            hit += 's'
        return "Result(program='%s', target='%s', id='%s', %i %s)" % \
                (self.program, self.target, self.id, len(self), hit)

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
            return self._hit_key_function(hit_key) in self._hits
        return hit_key in self._hits

    def __len__(self):
        return len(self._hits)

    def __nonzero__(self):
        return bool(self._hits)

    def __reversed__(self):
        hits = reversed(list(self.hits))
        obj =  self.__class__(self.id, hits, self.meta, self.program, \
                self.target, self._hit_key_function)
        self._transfer_attrs(obj)
        return obj

    def __setitem__(self, hit_key, hit):
        """Custom Search __setitem__.

        Hit key must be a string and hit must be a Hit object.

        """
        # only accept string keys
        if not isinstance(hit_key, basestring):
            raise TypeError("Result object keys must be a string.")
        # hit must be a Hit object
        if not isinstance(hit, Hit):
            raise TypeError("Result objects can only contain Hit objects.")
        # and it must have the same query ID as this object's ID
        if hit.query_id != self.id:
            raise ValueError("Expected Hit with query ID '%s', found '%s' "
                    "instead." % (self.id, hit.query_id))

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

            # return a list of hsps if hsp_rank is a slice object, instead
            # of simply slicing the hit object (which would return another hit
            # object)
            return hit.hsps[hsp_rank]

        # retrieval using slice objects returns another Result object
        elif isinstance(hit_key, slice):
            # should we return just a list of Hits instead of a full blown
            # Result object if it's a slice?
            hits = list(self.hits)[hit_key]
            obj =  self.__class__(self.id, hits, self.meta, self.program, \
                    self.target, self._hit_key_function)
            self._transfer_attrs(obj)
            return obj

        # if key is an int, then retrieve the Hit at the int index
        elif isinstance(hit_key, int):
            return list(self.hits)[hit_key]

        # if key is a string, then do a regular dictionary retrieval
        return self._hits[hit_key]

    def __delitem__(self, hit_key):
        """Custom Search __delitem__

        If hit_key is a string, then the method will delete the Hit object whose
        ID matches the string. If key is an integer or a slice object, then
        the Hit objects within that range will be deleted.

        The method can also delete HSP objects within a given Hit object if
        hit_key is given as a tuple or list of two items: the Hit index and
        the HSP index.

        """
        if isinstance(hit_key, (list, tuple)):
            # if hit_key is a list or tuple, deletion is done at the hsp level
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

            # delete the HSPs from the hit and assign it to replace the old hit
            # TODO: this allows cases where a Hit object can have 0 HSPs ~ is
            # it ok to allow this? should we delete the hit object instead?
            del hit[hsp_rank]
            self._hits[self._hit_key_function(hit)] = hit
            return

        else:
            # if hit_key an integer or slice, get the corresponding key first
            # and put it into a list
            if isinstance(hit_key, int):
                hit_keys = [list(self.hit_ids)[hit_key]]
            # the same, if it's a slice
            elif isinstance(hit_key, slice):
                hit_keys = list(self.hit_ids)[hit_key]
            # otherwise put it in a list
            else:
                hit_keys = [hit_key]

            for key in hit_keys:
                del self._hits[key]
            return

    def append(self, hit):
        """Adds a Hit object to the end of Result.

        """
        # if a custom hit_key_function is supplied, use it to define th hit key
        if self._hit_key_function is not None:
            hit_key = self._hit_key_function(hit)
        else:
            hit_key = hit.id

        if hit_key not in self:
            self[hit_key] = hit
        else:
            raise ValueError("Hit '%s' already present in this Result." % \
                    hit_key)

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


class Hit(_StickyObject):

    """Class representing the entire database entry of a sequence match.

    """

    # attributes we don't want to transfer when creating a new Hit class
    # from this one
    _NON_STICKY_ATTRS = ('_hsps',)

    def __init__(self, hit_id, query_id, hsps=[]):
        """Initializes a Hit object.

        query_id -- String of the query name used to obtain this hit.
        hit_id -- String of unique identifier for this hit.
        hsps -- Iterable returning HSP objects.

        """
        self.query_id= query_id
        self.id = hit_id

        # Hit must contain a minimum of one HSP object
        if not hsps:
            raise ValueError("Hit object must contain at least one HSP object.")

        self._hsps = []
        for hsp in hsps:
            # validate each HSP
            self._validate_hsp(hsp)
            # and store it them as an instance attribute
            self._hsps.append(hsp)

    def __repr__(self):
        al = 'alignment'
        if len(self) != 1:
            al += 's'
        return "Hit(id='%s', query_id='%s', %i %s)" % (self.id, \
                self.query_id, len(self), al)

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
        obj = self.__class__(self.id, self.query_id, reversed(self._hsps))
        self._transfer_attrs(obj)
        return obj

    def __setitem__(self, idx, hsps):
        # handle case if hsps is a list of hsp
        if isinstance(hsps, (list, tuple)):
            for hsp in hsps:
                self._validate_hsp(hsp)
        else:
            self._validate_hsp(hsps)

        self._hsps[idx] = hsps

    def __getitem__(self, idx):
        # if key is slice, return a new Hit instance
        if isinstance(idx, slice):
            obj = self.__class__(self.id, self.query_id, self._hsps[idx])
            self._transfer_attrs(obj)
            return obj
        return self._hsps[idx]

    def __delitem__(self, idx):
        del self._hsps[idx]

    def _validate_hsp(self, hsp):
        if not isinstance(hsp, HSP):
            raise TypeError("Hit objects can only contain HSP objects.")
        if hsp.hit_id != self.id:
            raise ValueError("Expected HSP with hit ID '%s', found '%s' "
                    "instead." % (self.id, hsp.hit_id))
        if hsp.query_id != self.query_id:
            raise ValueError("Expected HSP with query ID '%s', found '%s' "
                    "instead." % (self.query_id, hsp.query_id))

    def append(self, hsp):
        self._validate_hsp(hsp)
        self._hsps.append(hsp)

    def pop(self, index=-1):
        return self._hsps.pop(index)

    def reverse(self):
        self._hsps.reverse()

    def sort(self, key=None, reverse=False):
        self._hsps.sort(key=key, reverse=reverse)


class HSP(_StickyObject):

    """Class representing high-scoring alignment regions of the query and hit.

    """

    # attributes we don't want to transfer when creating a new HSP class
    # from this one
    _NON_STICKY_ATTRS = ('hit', 'query', 'alignment',)

    def __init__(self, hit_id, query_id, hit_seq='', query_seq='', \
            alphabet=single_letter_alphabet):
        """Initializes an HSP object.

        hit_id -- String, Hit ID of the HSP object.
        query_id -- String of the search query ID.
        hit_seq -- String or SeqRecord object of the aligned Hit sequence.
        query_seq -- String or SeqRecord object of the aligned query sequence.

        """
        self.hit_id = hit_id
        self.query_id = query_id
        self._alphabet = alphabet

        # only accept hit_seq and query_seq as string or SeqRecord objects
        if not isinstance(hit_seq, (SeqRecord, basestring)) or \
                not isinstance(query_seq, (SeqRecord, basestring)):
            raise TypeError("HSP sequence must be a string or a "
                    "SeqRecord object.")

        # only initialize MultipleSeqAlignment if hit_seq and query_seq is given
        if hit_seq and query_seq:

            # if hit_seq is a string, create a new SeqRecord object
            if isinstance(hit_seq, basestring):
                hit = SeqRecord(Seq(hit_seq, alphabet), id=hit_id, \
                        name='hit', description='aligned hit sequence')
            # otherwise hit is the hit_seq
            else:
                hit = hit_seq

            # same thing for query
            if isinstance(query_seq, basestring):
                query = SeqRecord(Seq(query_seq, alphabet), id=query_id, \
                        name='query', description='aligned query sequence')
            else:
                query = query_seq

            self.query = query
            self.hit = hit
            self.alignment = MultipleSeqAlignment([query, hit], alphabet)

        else:
            self.query, self.hit, self.alignment = None, None, None

    def __repr__(self):
        info = "hit_id='%s', query_id='%s'" % (self.hit_id, self.query_id)

        try:
            info += ", evalue=%s" % str(self.evalue)
        except AttributeError:
            pass

        try:
            info += ", %i-column alignment" % len(self)
        except TypeError:
            pass

        return "HSP(%s)" % (info)

    def __len__(self):
        if hasattr(self, 'alignment'):
            return len(self.query)
        else:
            try:
                return self.length
            except AttributeError:
                raise TypeError("HSP objects without alignment does not have any length.")

    def __getitem__(self, idx):
        if hasattr(self, 'alignment'):
            obj = self.__class__(self.hit_id, self.query_id, self.hit[idx], \
                    self.query[idx], self._alphabet)
            self._transfer_attrs(obj)
            return obj
        else:
            raise TypeError("Slicing for HSP objects without alignment is not supported.")

    def __delitem__(self, idx):
        raise TypeError("HSP objects are read-only.")

    def __setitem__(self, idx, value):
        raise TypeError("HSP objects are read-only.")

    def __iter__(self):
        raise TypeError("HSP objects do not support iteration.")


class SearchIndexer(object):

    """Iterator that returns file positions of results in a Search output file.

    """
