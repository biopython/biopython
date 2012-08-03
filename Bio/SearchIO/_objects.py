# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO objects to model homology search program outputs (PRIVATE)."""

import warnings
from itertools import chain

from Bio import BiopythonWarning
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import single_letter_alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio._py3k import OrderedDict


class BaseSearchObject(object):

    """Abstract class for SearchIO objects."""

    _NON_STICKY_ATTRS = ()

    def _transfer_attrs(self, obj):
        """Transfer instance attributes to the given object.

        This method is used to transfer attributes set externally (for example
        using setattr()) to a new object created from this one (for example
        from slicing).

        The reason this method is necessary is because different parsers will
        set different attributes for each QueryResult, Hit, or HSP object they use,
        depending on the attributes they found in the search output file.
        Ideally, we want these attributes to 'stick' with any new instance
        object created from the original one.

        """
        # list of attribute names we don't want to transfer
        for attr in self.__dict__.keys():
            if attr not in self._NON_STICKY_ATTRS:
                setattr(obj, attr, self.__dict__[attr])

    def _concat_display(string, max_len, concat_char):
        """Concatenates the given string for display."""
        if len(string) > max_len:
            return string[:max_len - len(concat_char)] + concat_char
        return string

    _concat_display = staticmethod(_concat_display)

    def _attr_display(obj, attr, fmt=None, fallback='?'):
        """Returns a string of the given object's attribute."""
        if hasattr(obj, attr):
            if fmt is not None:
                return fmt % getattr(obj, attr)
            return str(getattr(obj, attr))
        return fallback

    _attr_display = staticmethod(_attr_display)


class QueryResult(BaseSearchObject):

    """Class representing search results from a single query.

    The QueryResult object is a container for storing all search hits from a single
    search query. It is the top-level object returned by SearchIO's two main
    functions, SearchIO.read:

    >>> from Bio import SearchIO
    >>> qresult = SearchIO.read('tblastx_human_wnts.xml', 'blast-xml')
    >>> qresult
    QueryResult(program='TBLASTX', target='refseq_mrna', id='gi|195230749:301-1383', 5 hits)

    and SearchIO.parse:

    >>> qresults = SearchIO.parse('tblastx_human_wnts.xml', 'blast-xml')
    >>> qresult = qresults.next()
    >>> qresult
    QueryResult(program='TBLASTX', target='refseq_mrna', id='gi|195230749:301-1383', 5 hits)

    QueryResult is basically a container of the hits (see Hit objects) from one
    query sequence. Its length is how many hits it has and iteration over a
    QueryResult object returns Hit objects.

    >>> len(qresult)
    5
    >>> for hit in qresult:
    ...     print hit.id
    ...
    gi|195230749|ref|NM_003391.2|
    gi|281183280|ref|NM_001168718.1|
    gi|281182577|ref|NM_001168597.1|
    gi|274325896|ref|NM_001168687.1|
    gi|209529663|ref|NM_001135848.1|

    QueryResult objects behaves like a hybrid of Python's built-in list and
    dictionary, enabling retrieval of search hits using its index (integer) or
    its key (string, defaults to ID).

    Indexing using integers works exactly the same as Python lists:

    >>> first_hit = qresult[0]
    >>> first_hit
    Hit(id='gi|195230749|ref|NM_003391.2|', query_id='gi|195230749:301-1383', 10 alignments)

    >>> last_hit = qresult[-1]
    >>> last_hit
    Hit(id='gi|209529663|ref|NM_001135848.1|', query_id='gi|195230749:301-1383', 10 alignments)

    Indexing using hit IDs works just like Python dictionaries. This is useful
    if you know what you are expecting from the search beforehand.

    >>> qresult['gi|195230749|ref|NM_003391.2|']
    Hit(id='gi|195230749|ref|NM_003391.2|', query_id='gi|195230749:301-1383', 10 alignments)

    To get a list of all the hits contained in a QueryResult object, you can use
    the hits attribute. To obtain all the hit keys, the hit_keys attribute
    is used.

    >>> qresult.hits
    [Hit(id='gi|195230749|ref|NM_003391.2|', query_id='gi|195230749:301-1383', 10 alignments), Hit(id='gi|281183280|ref|NM_001168718.1|', query_id='gi|195230749:301-1383', 10 alignments), Hit(id='gi|281182577|ref|NM_001168597.1|', query_id='gi|195230749:301-1383', 10 alignments), Hit(id='gi|274325896|ref|NM_001168687.1|', query_id='gi|195230749:301-1383', 10 alignments), Hit(id='gi|209529663|ref|NM_001135848.1|', query_id='gi|195230749:301-1383', 10 alignments)]

    >>> qresult.hit_keys
    [u'gi|195230749|ref|NM_003391.2|', u'gi|281183280|ref|NM_001168718.1|', u'gi|281182577|ref|NM_001168597.1|', u'gi|274325896|ref|NM_001168687.1|', u'gi|209529663|ref|NM_001135848.1|']

    Similar to Python lists, you can also slice QueryResult objects. However,
    instead of returning a list, slicing will return a new QueryResult object.
    The new QueryResult object will have its hits sliced accordingly and attributes
    present in the unsliced QueryResult object are retained.

    >>> qresult
    QueryResult(program='TBLASTX', target='refseq_mrna', id='gi|195230749:301-1383', 5 hits)
    >>> sliced_qresult = qresult[:3]
    >>> sliced_qresult
    QueryResult(program='TBLASTX', target='refseq_mrna', id='gi|195230749:301-1383', 3 hits)
    >>> len(qresult)
    5
    >>> len(sliced_qresult)
    3
    >>> qresult.program
    u'TBLASTX'
    >>> qresult.program == sliced_qresult.program
    True
    >>> qresult[0] in sliced_qresult
    True
    >>> qresult[4] in sliced_qresult
    False

    You can check whether a hit is present in a QueryResult object using its key
    (defaults to hit ID) or the Hit object itself.

    >>> hit = qresult[0]
    >>> hit in qresult
    True
    >>> hit.id in qresult
    True

    Finally, QueryResult objects has other methods normally used in Python lists:
    append(), index(), pop(), and sort(). Consult their documentation for more
    information.

    """

    # attributes we don't want to transfer when creating a new QueryResult class
    # from this one
    _NON_STICKY_ATTRS = ('_hits',)

    def __init__(self, id=None, hits=[], hit_key_function=lambda hit: hit.id):
        """Initializes a QueryResult object.

        Arguments:
        query_id -- String of query sequence ID.
        hits     -- Iterator returning Hit objects.
        hit_key_function -- Function to define hit keys, defaults to a function
                            that return Hit object IDs.

        """
        if id is None:
            raise ValueError("Query ID string is required for QueryResult "
                    "creation")

        self._id = id
        self._hit_key_function = hit_key_function
        self._hits = OrderedDict()
        # default program, target, version, and description
        self.program = '<unknown>'
        self.target = '<unknown>'
        self.version = '<unknown>'
        self._description = ''

        # validate Hit objects and fill up self._hits
        for hit in hits:
            # validation is handled by __setitem__
            self.append(hit)

    # handle Python 2 OrderedDict behavior
    if hasattr(OrderedDict, 'iteritems'):

        def __iter__(self):
            return iter(self.iterhits())

        def _hits_get(self):
            return self._hits.values()

        hits = property(fget=_hits_get, \
                doc="""Returns a list of Hit objects contained by this object.""")

        def _hit_keys_get(self):
            return self._hits.keys()

        hit_keys = property(fget=_hit_keys_get, \
                doc="""Returns a list of Hit IDs contained by this object.""")

        def _items_get(self):
            return self._hits.items()

        items = property(fget=_items_get, \
            doc="""Returns a list of tuples of Hit ID and Hit object contained by this object.""")

        def iterhits(self):
            """Returns an iterator over the Hit objects."""
            for hit in self._hits.itervalues():
                yield hit

        def iterhit_keys(self):
            """Returns an iterator over the ID of the Hit objects."""
            for hit_id in self._hits.iterkeys():
                yield hit_id

        def iteritems(self):
            """Returns an iterator of tuples of Hit ID and Hit objects."""
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
        def hit_keys(self):
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

    def __repr__(self):
        return "QueryResult(id=%r, %r hits)" % (self.id, len(self))

    def __str__(self):
        lines = []

        # set program and version line
        lines.append('Program: %s (%s)' % (self.program, self.version))

        # set query id line
        qid_line = '  Query: %s' % self.id
        if hasattr(self, 'seq_len'):
            qid_line += ' (%i)' % self.seq_len
        if self.description:
            qid_line += QueryResult._concat_display('\n         %s' % \
                    self.description, 80, '...')
        lines.append(qid_line)

        # set target line
        lines.append(' Target: %s' % self.target)

        # set hit lines
        if not self.hits:
            lines.append('   Hits: 0')
        else:
            lines.append('   Hits: %s  %s  %s' % ('-'*4, '-'*5, '-'*58))
            pattern = '%13s  %5s  %56s'
            lines.append(pattern % ('#', '# HSP', 'ID + description'.ljust(58)))
            lines.append(pattern % ('-'*4, '-'*5, '-'*58))
            for idx, hit in enumerate(self.hits):
                if idx < 30:
                    hid_line = '%s  %s' % (hit.id, hit.description)
                    if len(hid_line) > 58:
                        hid_line = hid_line[:55] + '...'
                    lines.append(pattern % (idx, str(len(hit)), hid_line.ljust(58)))
                elif idx > len(self.hits) - 4:
                    hid_line = '%s  %s' % (hit.id, hit.description)
                    if len(hid_line) > 58:
                        hid_line = hid_line[:55] + '...'
                    lines.append(pattern % (idx, str(len(hit)), hid_line.ljust(58)))
                elif idx == 30:
                    lines.append('%14s' % '~~~')

        return '\n'.join(lines)

    def __reversed__(self):
        hits = reversed(list(self.hits))
        obj =  self.__class__(self.id, hits, self._hit_key_function)
        self._transfer_attrs(obj)
        return obj

    def __setitem__(self, hit_key, hit):
        """Custom Search object item assignment.

        Hit key must be a string and hit must be a Hit object.

        """
        # only accept string keys
        if not isinstance(hit_key, basestring):
            raise TypeError("QueryResult object keys must be a string.")
        # hit must be a Hit object
        if not isinstance(hit, Hit):
            raise TypeError("QueryResult objects can only contain Hit objects.")
        # and it must have the same query ID as this object's ID
        if hit.query_id != self.id:
            raise ValueError("Expected Hit with query ID '%s', found '%s' "
                    "instead." % (self.id, hit.query_id))

        self._hits[hit_key] = hit

    def __getitem__(self, hit_key):
        """Custom Search object item retrieval.

        Allows value retrieval by its key, location index, or a slice of
        location index.

        """
        # retrieval using slice objects returns another QueryResult object
        if isinstance(hit_key, slice):
            # should we return just a list of Hits instead of a full blown
            # QueryResult object if it's a slice?
            hits = list(self.hits)[hit_key]
            obj =  self.__class__(self.id, hits, self._hit_key_function)
            self._transfer_attrs(obj)
            return obj

        # if key is an int, then retrieve the Hit at the int index
        elif isinstance(hit_key, int):
            return list(self.hits)[hit_key]

        # if key is a string, then do a regular dictionary retrieval
        return self._hits[hit_key]

    def __delitem__(self, hit_key):
        """Custom Search object item deletion.

        If hit_key is a string, then the method will delete the Hit object whose
        ID matches the string. If key is an integer or a slice object, then
        the Hit objects within that range will be deleted.

        """
        # if hit_key an integer or slice, get the corresponding key first
        # and put it into a list
        if isinstance(hit_key, int):
            hit_keys = [list(self.hit_keys)[hit_key]]
        # the same, if it's a slice
        elif isinstance(hit_key, slice):
            hit_keys = list(self.hit_keys)[hit_key]
        # otherwise put it in a list
        else:
            hit_keys = [hit_key]

        for key in hit_keys:
            del self._hits[key]
        return

    def _description_get(self):
        return self._description

    def _description_set(self, value):
        self._description = value
        # try to set descriptions of hsp.query.seq within
        for hit in self.hits:
            hit.query_description = value

    description = property(fget=_description_get, fset=_description_set)

    def _id_get(self):
        return self._id

    def _id_set(self, value):
        self._id = value
        # set all Hit IDs contained to have the new Query ID
        for hit in self.hits:
            hit.query_id = value

    id = property(fget=_id_get, fset=_id_set)

    def append(self, hit):
        """Adds a Hit object to the end of QueryResult.

        Argument:
        hit -- Hit object.

        The hit key used for the appended Hit object depends on the
        hit_key_function used to initialize the QueryResult object. Any Hit object
        appended must have the same query_id attribute as the QueryResult object's
        id attribute. If the hit key already exists, a ValueError will be
        raised.

        """
        # if a custom hit_key_function is supplied, use it to define th hit key
        if self._hit_key_function is not None:
            hit_key = self._hit_key_function(hit)
        else:
            hit_key = hit.id

        if hit_key not in self:
            self[hit_key] = hit
        else:
            raise ValueError("Hit '%s' already present in this QueryResult." % \
                    hit_key)

    def hit_filter(self, func=None):
        """Creates a new QueryResult object whose Hit objects pass the filter function."""
        hits = filter(func, self.hits)
        obj =  self.__class__(self.id, hits, self._hit_key_function)
        self._transfer_attrs(obj)
        return obj

    def hit_map(self, func=None):
        """Creates a new QueryResult object, mapping the given function to its Hits."""
        hits = (hit[:] for hit in self.hits)
        hits = map(func, hits)
        obj =  self.__class__(self.id, hits, self._hit_key_function)
        self._transfer_attrs(obj)
        return obj

    def hsp_filter(self, func=None):
        """Creates a new QueryResult object whose HSP objects pass the filter function."""
        hits = filter(None, (hit.filter(func) for hit in self.hits))
        obj =  self.__class__(self.id, hits, self._hit_key_function)
        self._transfer_attrs(obj)
        return obj

    def hsp_map(self, func=None):
        """Creates a new QueryResult object, mapping the given function to its HSPs."""
        hits = filter(None, (hit.map(func) for hit in self.hits[:]))
        obj =  self.__class__(self.id, hits, self._hit_key_function)
        self._transfer_attrs(obj)
        return obj

    # marker for default self.pop() return value
    # this method is adapted from Python's built in OrderedDict.pop
    # implementation
    __marker = object()

    def pop(self, hit_key=-1, default=__marker):
        """Removes the specified hit key and return the Hit object.

        Arguments:
        hit_key -- Integer index or string of hit key that points to a Hit
                   object.
        default -- Value that will be returned if the Hit object with the
                   specified index or hit key is not found.

        By default, pop will remove and return the last Hit object in the
        QueryResult object. To remove specific Hit objects, you can use its integer
        index or its hit key.

        >>> from Bio import SearchIO
        >>> qresult = SearchIO.read('tblastx_human_wnts.xml', 'blast-xml')
        >>> len(qresult)
        5
        >>> for hit in qresult:
        ...     print hit.id
        ...
        gi|195230749|ref|NM_003391.2|
        gi|281183280|ref|NM_001168718.1|
        gi|281182577|ref|NM_001168597.1|
        gi|274325896|ref|NM_001168687.1|
        gi|209529663|ref|NM_001135848.1|

        >>> qresult.pop()
        Hit(id='gi|209529663|ref|NM_001135848.1|', query_id='gi|195230749:301-1383', 10 alignments)

        >>> qresult.pop(0)
        Hit(id='gi|195230749|ref|NM_003391.2|', query_id='gi|195230749:301-1383', 10 alignments)

        >>> qresult.pop('gi|281182577|ref|NM_001168597.1|')
        Hit(id='gi|281182577|ref|NM_001168597.1|', query_id='gi|195230749:301-1383', 10 alignments)

        >>> len(qresult)
        2

        """
        # if key is an integer (index)
        # get the ID for the Hit object at that index
        if isinstance(hit_key, int):
            # raise the appropriate error if there is no hit
            if not self:
                raise IndexError("pop from empty list")
            hit_key = list(self.hit_keys)[hit_key]

        try:
            return self._hits.pop(hit_key)
        except KeyError:
            # if key doesn't exist and no default is set, raise a KeyError
            if default is self.__marker:
                raise KeyError(hit_key)
        # if key doesn't exist but a default is set, return the default value
        return default

    def index(self, hit_key):
        """Returns the index of a given hit key, zero-based.

        Argument:
        hit_key -- String of hit key or Hit object.

        >>> from Bio import SearchIO
        >>> qresult = SearchIO.read('tblastx_human_wnts.xml', 'blast-xml')
        >>> qresult.index('gi|209529663|ref|NM_001135848.1|')
        4
        >>> hit = qresult['gi|209529663|ref|NM_001135848.1|']
        >>> qresult.index(hit)
        4
        >>> qresult.index('my_key')
        -1

        This method is useful for finding out the integer index (usually
        correlated with search rank) of a given hit key.

        """
        try:
            if isinstance(hit_key, Hit):
                return list(self.hit_keys).index(hit_key.id)
            return list(self.hit_keys).index(hit_key)
        except ValueError:
            return -1

    def sort(self, key=None, reverse=False, in_place=True):
        # no cmp argument to make sort more Python 3-like
        """Sorts the Hit objects.

        Arguments:
        key -- Function used to sort the Hit objects.
        reverse -- Boolean, whether to reverse the sorting or not.

        >>> from Bio import SearchIO
        >>> qresult = SearchIO.read('tblastx_human_wnts.xml', 'blast-xml')
        >>> for hit in qresult:
        ...     print hit.id
        ...
        gi|195230749|ref|NM_003391.2|
        gi|281183280|ref|NM_001168718.1|
        gi|281182577|ref|NM_001168597.1|
        gi|274325896|ref|NM_001168687.1|
        gi|209529663|ref|NM_001135848.1|

        >>> qresult.sort(reverse=True)
        >>> for hit in qresult:
        ...     print hit.id
        ...
        gi|209529663|ref|NM_001135848.1|
        gi|274325896|ref|NM_001168687.1|
        gi|281182577|ref|NM_001168597.1|
        gi|281183280|ref|NM_001168718.1|
        gi|195230749|ref|NM_003391.2|

        >>> qresult.sort(key=lambda hit: hit.id)
        >>> for hit in qresult:
        ...     print hit.id
        ...
        gi|195230749|ref|NM_003391.2|
        gi|209529663|ref|NM_001135848.1|
        gi|274325896|ref|NM_001168687.1|
        gi|281182577|ref|NM_001168597.1|
        gi|281183280|ref|NM_001168718.1|

        By default, sorting is based on the expect values of the Hit objects,
        from the smallest to the largest. If the Hit objects do not have any
        expect values (e.g. BLAT Hit objects), then no sorting is performed.

        The sort creates a new Hit container object, but appears to be
        in-place since the new Hit container replaces the old one.

        """
        if key is None:
            # if reverse is True, reverse the hits
            if reverse:
                sorted_hits = self.hits[::-1]
            # otherwise (default options) make a copy of the hits
            else:
                sorted_hits = self.hits[:]
        else:
            sorted_hits = sorted(self.hits, key=key, reverse=reverse)

        # if sorting is in-place, don't create a new QueryResult object
        if in_place:
            new_hits = OrderedDict()
            for hit in sorted_hits:
                new_hits[self._hit_key_function(hit)] = hit
            self._hits = new_hits
        # otherwise, return a new sorted QueryResult object
        else:
            obj =  self.__class__(self.id, sorted_hits, self._hit_key_function)
            self._transfer_attrs(obj)
            return obj



class Hit(BaseSearchObject):

    """Class representing a single database hit of a search result.

    Hit objects are the second-level container in the SearchIO module. They
    are the objects contained within a QueryResult object (see QueryResult). Each Hit
    object is uniquely identified by its ID and the query ID that results in
    its creation.

    >>> from Bio import SearchIO
    >>> qresult = SearchIO.read('tblastx_human_wnts.xml', 'blast-xml')
    >>> hit = qresult[0]
    >>> hit
    Hit(id='gi|195230749|ref|NM_003391.2|', query_id='gi|195230749:301-1383', 10 alignments)

    Hit objects themselves are container for the basic SearchIO unit: the
    HSP object (see HSP). Since these HSPs usually do not have any unique
    IDs, Hit objects behave very similar to a Python built-in list.

    The length of a Hit object is how many HSPs it contains and iteration
    over Hit objects return these HSPs:

    >>> len(hit)
    10
    >>> for hsp in hit:
    ...     print hsp.hit_id, hsp.evalue, len(hsp)
    gi|195230749|ref|NM_003391.2| 0.0 340
    gi|195230749|ref|NM_003391.2| 0.0 253
    gi|195230749|ref|NM_003391.2| 0.0 69
    gi|195230749|ref|NM_003391.2| 0.0 361
    gi|195230749|ref|NM_003391.2| 0.0 178
    gi|195230749|ref|NM_003391.2| 0.0 161
    gi|195230749|ref|NM_003391.2| 0.0 237
    gi|195230749|ref|NM_003391.2| 0.0 106
    gi|195230749|ref|NM_003391.2| 0.0 288
    gi|195230749|ref|NM_003391.2| 0.0 28

    Like built-in Python lists, you can index Hit objects with integers or
    slice them. Slicing a Hit object will return a new Hit object with its
    HSPs properly sliced but any other attributes retained.

    >>> hit[0]
    HSP(hit_id='gi|195230749|ref|NM_003391.2|', query_id='gi|195230749:301-1383', evalue=0.0, 340-column alignment)
    >>> hit[-1]
    HSP(hit_id='gi|195230749|ref|NM_003391.2|', query_id='gi|195230749:301-1383', evalue=0.0, 28-column alignment)

    >>> hit
    Hit(id='gi|195230749|ref|NM_003391.2|', query_id='gi|195230749:301-1383', 10 alignments)
    >>> sliced_hit = hit[:3]
    >>> sliced_hit
    Hit(id='gi|195230749|ref|NM_003391.2|', query_id='gi|195230749:301-1383', 3 alignments)
    >>> len(hit)
    10
    >>> len(sliced_hit)
    3
    >>> hit.description
    u'Homo sapiens wingless-type MMTV integration site family member 2 (WNT2), transcript variant 1, mRNA'
    >>> hit.description == sliced_hit.description
    True
    >>> hit[0] in sliced_hit
    True
    >>> hit[6] in sliced_hit
    False

    You can check whether an hsp is present in a QueryResult object using the HSP
    object itself.

    >>> hsp = hit[0]
    >>> hsp in hit
    True

    Finally, similar to Python built-in list, the Hit object also has the
    append(), pop(), reverse(), and sort() method, which behave similar to
    their list method counterparts.

    """

    # attributes we don't want to transfer when creating a new Hit class
    # from this one
    _NON_STICKY_ATTRS = ('_items',)

    def __init__(self, id=None, query_id=None, hsps=[]):
        """Initializes a Hit object.

        Arguments:
        query_id -- String of the query name used to obtain this hit.
        hit_id   -- String of unique identifier for this hit.
        hsps     -- Iterable returning HSP objects.

        """
        if id is None:
            raise ValueError("Hit ID string is required for Hit creation")
        if query_id is None:
            raise ValueError("Query ID string is required for Hit creation")

        self._id = id
        self._query_id= query_id
        self._description = ''
        self._query_description = ''

        self._items = []
        for hsp in hsps:
            # validate each HSP
            self._validate_hsp(hsp)
            # and store it them as an instance attribute
            self.append(hsp)

    def __repr__(self):
        return "Hit(id=%r, query_id=%r, %r hsps)" % (self.id, self.query_id, \
                len(self))

    def __iter__(self):
        return iter(self.hsps)

    def __len__(self):
        return len(self.hsps)

    def __nonzero__(self):
        return bool(self.hsps)

    def __contains__(self, hsp):
        return hsp in self._items

    def __str__(self):
        lines = []

        # set query id line
        lines.append('Query: %s' % self.query_id)

        # set hit id line
        hid_line = '  Hit: %s' % self.id
        if hasattr(self, 'seq_len'):
            hid_line += ' (%i)' % self.seq_len
        if self.description:
            hid_line += Hit._concat_display('\n       %s' % self.description, \
                    80, '...')
        lines.append(hid_line)

        # set hsp line and table
        if not self.hsps:
            lines.append(' HSPs: ?')
        else:
            lines.append(' HSPs: %s  %s  %s  %s  %s  %s' % \
                    ('-'*4, '-'*8, '-'*9, '-'*6, '-'*18, '-'*18))
            pattern = '%11s  %8s  %9s  %6s  %18s  %18s'
            lines.append(pattern % ('#', 'E-value', 'Bit score', 'Span', \
                    'Query range', 'Hit range'))
            lines.append(pattern % ('-'*4, '-'*8, '-'*9, '-'*6, '-'*18, '-'*18))
            for idx, hsp in enumerate(self.hsps):
                # evalue
                evalue = Hit._attr_display(hsp, 'evalue', fmt='%.2g')
                # bitscore
                bitscore = Hit._attr_display(hsp, 'bitscore', fmt='%.2f')
                # alignment length
                aln_span = Hit._attr_display(hsp, 'aln_len')
                # query region
                query_start = Hit._attr_display(hsp, 'query_start')
                query_end = Hit._attr_display(hsp, 'query_end')
                query_range = '%s:%s' % (query_start, query_end)
                # max column length is 18
                query_range = Hit._concat_display(query_range, 18, '~')
                # hit region
                hit_start = Hit._attr_display(hsp, 'hit_start')
                hit_end = Hit._attr_display(hsp, 'hit_end')
                hit_range = '%s:%s' % (hit_start, hit_end)
                hit_range = Hit._concat_display(hit_range, 18, '~')
                # append the hsp row
                lines.append(pattern % (str(idx), evalue, bitscore, aln_span, \
                        query_range, hit_range))

        return '\n'.join(lines)

    def __reversed__(self):
        obj = self.__class__(self.id, self.query_id, reversed(self._items))
        self._transfer_attrs(obj)
        return obj

    def __getitem__(self, idx):
        # if key is slice, return a new Hit instance
        if isinstance(idx, slice):
            print idx
            obj = self.__class__(self.id, self.query_id, self.hsps[idx])
            self._transfer_attrs(obj)
            return obj
        return self._items[idx]

    def __setitem__(self, idx, hsps):
        # handle case if hsps is a list of hsp
        if isinstance(hsps, (list, tuple)):
            for hsp in hsps:
                self._validate_hsp(hsp)
        else:
            self._validate_hsp(hsps)

        self._items[idx] = hsps

    def __delitem__(self, idx):
        del self._items[idx]

    ## hsp properties ##
    def _validate_hsp(self, hsp):
        """Validates an HSP object.

        Valid HSP objects have the same hit_id as the Hit object ID and the
        same query_id as the Hit object's query_id.

        """
        if not isinstance(hsp, HSP):
            raise TypeError("Hit objects can only contain HSP objects.")
        if hsp.hit_id != self.id:
            raise ValueError("Expected HSP with hit ID %r, " \
                    "found %r instead." % (self.id, hsp.hit_id))
        if hsp.query_id != self.query_id:
            raise ValueError("Expected HSP with query ID %r, " \
                    "found %r instead." % (self.query_id, hsp.query_id))

    def _hsps_get(self):
        return self._items

    hsps = property(fget=_hsps_get)

    def _fragments_get(self):
        return list(chain.from_iterable(self._items))

    fragments = property(fget=_fragments_get)

    ## id and description properties ##
    def _id_get(self):
        return self._id

    def _id_set(self, value):
        self._id = value
        # set all HSP IDs contained to have the new Hit ID
        for hsp in self._items:
            hsp.hit_id = value

    id = property(fget=_id_get, fset=_id_set)

    def _query_id_get(self):
        return self._query_id

    def _query_id_set(self, value):
        self._query_id = value
        # set all HSP query IDs contained to have the new query ID
        for hsp in self._items:
            hsp.query_id = value

    query_id = property(fget=_query_id_get, fset=_query_id_set)

    def _description_get(self):
        return self._description

    def _description_set(self, value):
        self._description = value
        # cascade through contained HSP hit descriptions
        for hsp in self._items:
            hsp.hit_description = value

    description = property(fget=_description_get, fset=_description_set)

    def _query_description_get(self):
        return self._query_description

    def _query_description_set(self, value):
        self._query_description = value
        for hsp in self._items:
            hsp.query_description = value

    query_description = property(fget=_query_description_get, \
            fset=_query_description_set)

    ## public methods ##
    def append(self, hsp):
        self._validate_hsp(hsp)
        self._items.append(hsp)

    def filter(self, func=None):
        """Creates a new Hit object whose HSP objects pass the filter function."""
        hsps = filter(func, self.hsps)
        if hsps:
            obj = self.__class__(self.id, self.query_id, hsps)
            self._transfer_attrs(obj)
            return obj

    def index(self, hsp):
        return self._items.index(hsp)

    def map(self, func=None):
        """Creates a new Hit object, mapping the given function to its HSPs."""
        hsps = map(func, self.hsps[:])  # this creates a shallow copy
        if hsps:
            obj = self.__class__(self.id, self.query_id, hsps)
            self._transfer_attrs(obj)
            return obj

    def pop(self, index=-1):
        return self._items.pop(index)

    def sort(self, key=None, reverse=False, in_place=True):
        if in_place:
            self._items.sort(key=key, reverse=reverse)
        else:
            hsps = self.hsps[:]
            hsps.sort(key=key, reverse=reverse)
            obj = self.__class__(self.id, self.query_id, hsps)
            self._transfer_attrs(obj)
            return obj


class BaseHSP(BaseSearchObject):

    """Abstract base class for HSP objects."""

    def _display_hsp_header(self):
        """Prints the alignment header info."""
        lines = []
        # set query id line
        qid_line = self._concat_display('      Query: %s %s' % \
                (self.query_id, self.query_description), 80, '...')
        # set hit id line
        hid_line = self._concat_display('        Hit: %s %s' % \
                (self.hit_id, self.hit_description), 80, '...')
        lines.append(qid_line)
        lines.append(hid_line)

        # coordinates
        query_start = BaseHSP._attr_display(self, 'query_start')
        query_end = BaseHSP._attr_display(self, 'query_end')
        hit_start = BaseHSP._attr_display(self, 'hit_start')
        hit_end = BaseHSP._attr_display(self, 'hit_end')

        lines.append('Query range: %s:%s (%r)' % (query_start, query_end, \
                self.query_strand))
        lines.append('  Hit range: %s:%s (%r)' % (hit_start, hit_end, \
                self.hit_strand))

        return '\n'.join(lines)


class HSP(BaseHSP):

    """Class representing high-scoring region between query and hit."""

    # attributes we don't want to transfer when creating a new Hit class
    # from this one
    _NON_STICKY_ATTRS = ('_fragments', '_aln_len')

    def __init__(self, fragments):
        """Initializes an HSP object.

        Arguments:
        fragments -- List of HSPFragment objects.

        """
        if not fragments:
            raise ValueError("HSP objects must have at least one HSPFragment " \
                    "object.")
        # check that all fragments contain the same IDs, descriptions, and alphabet
        for attr in ('query_id', 'query_description', 'hit_id', \
                'hit_description', 'alphabet'):
            if len(set([getattr(fragment, attr) for fragment in fragments])) != 1:
                raise ValueError("HSP object can not contain fragments with " \
                        "more than one %s." % attr)

        self._items = []
        for fragment in fragments:
            self._validate_fragment(fragment)
            self._items.append(fragment)

    def __repr__(self):
        return "%s(%r)" % (self.__class__.__name__, self._items)

    def __iter__(self):
        return iter(self._items)

    def __len__(self):
        return len(self._items)

    def __nonzero__(self):
        return bool(self._items)

    def __str__(self):

        lines = []
        # set hsp info line
        statline = []
        # evalue
        evalue = HSP._attr_display(self, 'evalue', fmt='%.2g')
        statline.append('evalue ' +  evalue)
        # bitscore
        bitscore = HSP._attr_display(self, 'bitscore', fmt='%.2f')
        statline.append('bitscore ' +  bitscore)
        lines.append('Quick stats: ' + '; '.join(statline))

        if len(self.fragments) == 1:
            return '\n'.join([self._display_hsp_header(), '\n'.join(lines), \
                    self.fragments[0]._display_aln()])
        else:
            lines.append('  Fragments: %s  %s  %s  %s' % \
                    ('-'*3, '-'*18, '-'*18, '-'*18))
            pattern = '%16s  %18s  %18s  %18s'
            lines.append(pattern % ('#', 'Length', 'Query range', 'Hit range'))
            lines.append(pattern % ('-'*3, '-'*18, '-'*18, '-'*18))
            for idx, block in enumerate(self.fragments):
                # set hsp line and table
                # alignment span
                aln_span = HSP._attr_display(block, 'aln_len')
                # query region
                query_start = HSP._attr_display(block, 'query_start')
                query_end = HSP._attr_display(block, 'query_end')
                query_range = '%s:%s' % (query_start, query_end)
                # max column length is 18
                query_range = HSP._concat_display(query_range, 18, '~')
                # hit region
                hit_start = HSP._attr_display(block, 'hit_start')
                hit_end = HSP._attr_display(block, 'hit_end')
                hit_range = '%s:%s' % (hit_start, hit_end)
                hit_range = HSP._concat_display(hit_range, 18, '~')
                # append the hsp row
                lines.append(pattern % (str(idx), aln_span, query_range, hit_range))

            return self._display_hsp_header() + '\n' + '\n'.join(lines)

    def __getitem__(self, idx):
        # if key is slice, return a new HSP instance
        if isinstance(idx, slice):
            obj = self.__class__(self._items[idx])
            self._transfer_attrs(obj)
            return obj
        return self._items[idx]

    def __setitem__(self, idx, fragments):
        # handle case if hsps is a list of hsp
        if isinstance(fragments, (list, tuple)):
            for fragment in fragments:
                self._validate_fragment(fragment)
        else:
            self._validate_fragment(fragments)

        self._items[idx] = fragments

    def __delitem__(self, idx):
        # note that this may result in an empty HSP object, which should be
        # invalid
        del self._items[idx]

    def __contains__(self, fragment):
        return fragment in self._items

    def _validate_fragment(self, fragment):
        if not isinstance(fragment, HSPFragment):
            raise TypeError("HSP objects can only contain HSPFragment " \
                    "objects.")

    ## sequence / fragment properties ##
    def _fragments_get(self):
        return self._items

    fragments = property(fget=_fragments_get)

    def _fragment_get(self):
        return self._items[0]

    fragment = property(fget=_fragment_get)

    def _is_fragmented_get(self):
        return len(self) > 1

    is_fragmented = property(fget=_is_fragmented_get)

    def _hits_get(self):
        return [fragment.hit for fragment in self.fragments]

    hits = property(fget=_hits_get)

    def _hit_get(self):
        return self._items[0].hit

    hit = property(fget=_hit_get)

    def _queries_get(self):
        return [fragment.query for fragment in self.fragments]

    queries = property(fget=_queries_get)

    def _query_get(self):
        return self._items[0].query

    query = property(fget=_query_get)

    def _alignments_get(self):
        return [fragment.alignment for fragment in self.fragments]

    alignments = property(fget=_alignments_get)

    def _alignment_get(self):
        return self._items[0].alignment

    alignment = property(fget=_alignment_get)

    def _aln_len_get(self):
        # length of all alignments
        # alignment span can be its own attribute, or computed from
        # query / hit length
        if len(self) == 1:
            self._aln_len = self._items[0].aln_len
        else:
            if not hasattr(self, '_aln_len'):
                if all([query.seq for query in self.queries]):
                    self._aln_len = sum([len(query.seq) for query in self.queries])
                elif all([hit.seq for hit in self.hits]):
                    self._aln_len = sum([len(hit.seq) for hit in self.hits])
                else:
                    self._aln_len = None

        return self._aln_len

    def _aln_len_set(self, value):
        self._aln_len = value

    aln_len = property(fget=_aln_len_get, fset=_aln_len_set)

    ## id and description properties ##
    def _set_id_or_desc(self, value, seq_type, attr):
        assert seq_type in ('query', 'hit')
        assert attr in ('id', 'description')
        attr_name = '%s_%s' % (seq_type, attr)
        # set attr for fragments
        for fragment in self.fragments:
            setattr(fragment, attr_name, value)

    def _hit_description_get(self):
        return self._items[0].hit_description

    def _hit_description_set(self, value):
        self._set_id_or_desc(value, 'hit', 'description')

    hit_description = property(fget=_hit_description_get, \
            fset=_hit_description_set)

    def _query_description_get(self):
        return self._items[0].query_description

    def _query_description_set(self, value):
        self._set_id_or_desc(value, 'query', 'description')

    query_description = property(fget=_query_description_get, \
            fset=_query_description_set)

    def _hit_id_get(self):
        return self._items[0].hit_id

    def _hit_id_set(self, value):
        self._set_id_or_desc(value, 'hit', 'id')

    hit_id = property(fget=_hit_id_get, \
            fset=_hit_id_set)

    def _query_id_get(self):
        return self._items[0].query_id

    def _query_id_set(self, value):
        self._set_id_or_desc(value, 'query', 'id')

    query_id = property(fget=_query_id_get, \
            fset=_query_id_set)

    ## strand properties ##
    def _hit_strands_get(self):
        return [fragment.hit_strand for fragment in self.fragments]

    hit_strands = property(fget=_hit_strands_get)

    def _hit_strand_get(self):
        return self._items[0].hit_strand

    hit_strand = property(fget=_hit_strand_get)

    def _query_strands_get(self):
        return [fragment.query_strand for fragment in self.fragments]

    query_strands = property(fget=_query_strands_get)

    def _query_strand_get(self):
        return self._items[0].query_strand

    query_strand = property(fget=_query_strand_get)

    ## frame properties ##
    def _hit_frames_get(self):
        return [fragment.hit_frame for fragment in self.fragments]

    hit_frames = property(fget=_hit_frames_get)

    def _hit_frame_get(self):
        return self._items[0].hit_frame

    hit_frame = property(fget=_hit_frame_get)

    def _query_frames_get(self):
        return [fragment.query_frame for fragment in self.fragments]

    query_frames = property(fget=_query_frames_get)

    def _query_frame_get(self):
        return self._items[0].query_frame

    query_frame = property(fget=_query_frame_get)

    ## coordinate properties ##
    def _get_coords(self, seq_type, coord_type):
        assert seq_type in ('hit', 'query')
        assert coord_type in ('start', 'end')
        coord_name = '%s_%s' % (seq_type, coord_type)
        coords = [getattr(fragment, coord_name) for fragment in self.fragments]
        if None in coords:
            warnings.warn("'None' exist in %s coordinates; ignored" % \
                    (coord_name), BiopythonWarning)
        return coords

    def _hit_start_get(self):
        return min(self._get_coords('hit', 'start'))

    hit_start = property(fget=_hit_start_get)

    def _query_start_get(self):
        return min(self._get_coords('query', 'start'))

    query_start = property(fget=_query_start_get)

    def _hit_end_get(self):
        return max(self._get_coords('hit', 'end'))

    hit_end = property(fget=_hit_end_get)

    def _query_end_get(self):
        return max(self._get_coords('query', 'end'))

    query_end = property(fget=_query_end_get)

    ## coordinate-dependent properties ##
    def _hit_span_get(self):
        try:
            return self.hit_end - self.hit_start
        except TypeError:  # triggered if any of the coordinates are None
            return None

    hit_span = property(fget=_hit_span_get)

    def _query_span_get(self):
        try:
            return self.query_end - self.query_start
        except TypeError:  # triggered if any of the coordinates are None
            return None

    query_span = property(fget=_query_span_get)

    def _hit_ranges_get(self):
        return [fragment.hit_range for fragment in self.fragments]

    hit_ranges = property(fget=_hit_ranges_get)

    def _hit_range_get(self):
        return (self.hit_start, self.hit_end)

    hit_range = property(fget=_hit_range_get)

    def _query_ranges_get(self):
        return [fragment.query_range for fragment in self.fragments]

    query_ranges = property(fget=_query_ranges_get)

    def _query_range_get(self):
        return (self.query_start, self.query_end)

    query_range = property(fget=_query_range_get)


class HSPFragment(BaseHSP):

    """Class representing a fragment of matching hit-query sequence."""

    def __init__(self, hit_id=None, query_id=None, hit='', query='', \
            alphabet=single_letter_alphabet, aln_annotation={}):

        if hit_id is None:
            raise ValueError("Hit ID string is required for HSPFragment creation.")
        if query_id is None:
            raise ValueError("Query ID string is required for HSPFragment creation.")

        self.alphabet = alphabet
        for seq_type in ('query', 'hit'):
            # self.query or self.hit
            if eval(seq_type):
                setattr(self, seq_type, eval(seq_type))
            else:
                setattr(self, seq_type, None)
            # query or hit attributes
            for attr in ('strand', 'frame', 'start', 'end'):
                setattr(self, '%s_%s' % (seq_type, attr), None)
            attr_name = '%s_id' % seq_type
            setattr(self, attr_name, eval(attr_name))
        self.alignment_annotation = aln_annotation
        self.hit_description = ''
        self.query_description = ''

    def __repr__(self):
        info = "hit_id=%r, query_id=%r" % (self.hit_id, self.query_id)
        try:
            info += ", %i columns" % len(self)
        except TypeError:
            pass
        return "%s(%s)" % (self.__class__.__name__, info)

    def __len__(self):
        return self.aln_len

    def __str__(self):
        return self._display_hsp_header() + '\n' + self._display_aln()

    def __getitem__(self, idx):
        if self.alignment is not None:
            obj = self.__class__(
                    hit_id=self.hit_id, query_id=self.query_id, \
                    alphabet=self.alphabet)
            # transfer query and hit attributes
            if self.query is not None:
                obj.query = self.query[idx]
            if self.hit is not None:
                obj.hit = self.hit[idx]
            if self.query_description:
                obj.query_description = self.query_description
            if self.hit_description:
                obj.hit_description = self.hit_description
            # alignment annotation should be transferred, since we can compute
            # the resulting annotation
            if hasattr(self, 'alignment_annotation'):
                obj.alignment_annotation = {}
                for key, value in self.alignment_annotation.items():
                    assert len(value[idx]) == len(obj)
                    obj.alignment_annotation[key] = value[idx]
            return obj
        else:
            raise TypeError("Slicing for HSP objects without "
                    "alignment is not supported.")

    def _display_aln(self):
        lines = []
        # alignment length
        aln_len = HSPFragment._attr_display(self, 'aln_len')
        lines.append('   Fragment: %s columns' % aln_len)
        # sequences
        if hasattr(self, 'query') and hasattr(self, 'hit'):
            try:
                qseq = str(self.query.seq)
            except AttributeError:  # query is None
                qseq = '?'
            try:
                hseq = str(self.hit.seq)
            except AttributeError:  # hit is None
                hseq = '?'

            # homology line
            homol = ''
            if 'homology' in self.alignment_annotation:
                homol = self.alignment_annotation['homology']

            if self.aln_len <= 67:
                lines.append("%10s - %s" % ('Query', qseq))
                if homol:
                    lines.append("             %s" % homol)
                lines.append("%10s - %s" % ('Hit', hseq))
            else:
                # adjust continuation character length, so we don't display
                # the same residues twice
                if self.aln_len - 66 > 3:
                    cont = '~' * 3
                else:
                    cont = '~' * (self.aln_span - 66)
                lines.append("%10s - %s%s%s" % ('Query', \
                                qseq[:59], cont, qseq[-5:]))
                if homol:
                    lines.append("             %s%s%s" % \
                            (homol[:59], cont, homol[-5:]))
                lines.append("%10s - %s%s%s" % ('Hit', \
                                hseq[:59], cont, hseq[-5:]))

        return '\n'.join(lines)

    ## sequence properties ##
    def _prep_seq(self, seq, seq_type):
        """Transforms a sequence into a SeqRecord object).

        Argument:
        seq -- String of sequence.
        seq_type -- String of sequence type, must be 'hit' or 'query'

        """
        assert seq_type in ('hit', 'query')
        if seq is None: return seq # return immediately if seq is None
        # check length if the opposite sequence is not None
        opp_type = 'query' if seq_type == 'query' else 'query'
        opp_seq = getattr(self, opp_type, None)
        if opp_seq is not None:
            if len(seq) != len(opp_seq):
                raise ValueError("Sequence lengths do not match. Expected: " \
                        "%r (%s); found: %r (%s)." % (len(opp_seq), opp_type, \
                        len(seq), seq_type))

        seq_name = 'aligned %s sequence' % seq_type
        if isinstance(seq, SeqRecord):
            seq.name = seq_name
            return seq
        elif isinstance(seq, basestring):
            return SeqRecord(Seq(seq, self.alphabet), name=seq_name)
        else:
            raise TypeError("%s sequence must be a string or a "
                    "SeqRecord object." % seq_type.capitalize())

    def _hit_get(self):
        return self._hit

    def _hit_set(self, value):
        self._hit = self._prep_seq(value, 'hit')

    hit = property(fget=_hit_get, fset=_hit_set)

    def _query_get(self):
        return self._query

    def _query_set(self, value):
        self._query = self._prep_seq(value, 'query')

    query = property(fget=_query_get, fset=_query_set)

    def _alignment_get(self):
        if self.query is None and self.hit is None:
            return None
        elif self.hit is None:
            return MultipleSeqAlignment([self.query], self.alphabet)
        elif self.query is None:
            return MultipleSeqAlignment([self.hit], self.alphabet)
        else:
            return MultipleSeqAlignment([self.query, self.hit], self.alphabet)

    alignment = property(fget=_alignment_get)

    def _aln_len_get(self):
        # length of alignment (gaps included)
        # alignment span can be its own attribute, or computed from
        # query / hit length
        if not hasattr(self, '_aln_len'):
            if self.query is not None:
                self._aln_len = len(self.query)
            elif self.hit is not None:
                self._aln_len = len(self.hit)
            else:
                self._aln_len = None

        return self._aln_len

    def _aln_len_set(self, value):
        self._aln_len = value

    aln_len = property(fget=_aln_len_get, fset=_aln_len_set)

    ## id and description properties ##
    def _set_id_or_desc(self, value, seq_type, attr):
        assert seq_type in ('query', 'hit')
        assert attr in ('id', 'description')
        # set own attribute
        setattr(self, '_%s_%s' % (seq_type, attr), value)
        # and the seqrecord object's, if it exists
        if getattr(self, seq_type) is not None:
            seq = getattr(self, seq_type)
            setattr(seq, attr, value)

    def _hit_description_get(self):
        return self._hit_description

    def _hit_description_set(self, value):
        self._set_id_or_desc(value, 'hit', 'description')

    hit_description = property(fget=_hit_description_get, \
            fset=_hit_description_set)

    def _query_description_get(self):
        return self._query_description

    def _query_description_set(self, value):
        self._set_id_or_desc(value, 'query', 'description')

    query_description = property(fget=_query_description_get, \
            fset=_query_description_set)

    def _hit_id_get(self):
        return self._hit_id

    def _hit_id_set(self, value):
        self._set_id_or_desc(value, 'hit', 'id')

    hit_id = property(fget=_hit_id_get, \
            fset=_hit_id_set)

    def _query_id_get(self):
        return self._query_id

    def _query_id_set(self, value):
        self._set_id_or_desc(value, 'query', 'id')

    query_id = property(fget=_query_id_get, \
            fset=_query_id_set)

    ## strand properties ##
    def _prep_strand(self, strand):
        # follow SeqFeature's convention
        if not strand in (-1, 0, 1, None):
            raise ValueError("Strand should be -1, 0, 1, or None; not %r" % \
                    strand)
        return strand

    def _get_strand(self, seq_type):
        assert seq_type in ('hit', 'query')
        strand = getattr(self, '_%s_strand' % seq_type)

        if strand is None:
            # try to compute strand from frame
            frame = getattr(self, '%s_frame' % seq_type)
            if frame is not None:
                try:
                    strand = frame / abs(frame)
                except ZeroDivisionError:
                    strand = 0
                setattr(self, '%s_strand' % seq_type, strand)

        return strand

    def _hit_strand_get(self):
        return self._get_strand('hit')

    def _hit_strand_set(self, value):
        self._hit_strand = self._prep_strand(value)

    hit_strand = property(fget=_hit_strand_get, fset=_hit_strand_set)

    def _query_strand_get(self):
        return self._get_strand('query')

    def _query_strand_set(self, value):
        self._query_strand = self._prep_strand(value)

    query_strand = property(fget=_query_strand_get, fset=_query_strand_set)

    ## frame properties ##
    def _prep_frame(self, frame):
        if not frame in (-3, -2, -1, 0, 1, 2, 3, None):
            raise ValueError("Strand should be an integer between -3 and 3, "
                    "or None; not %r" % frame)
        return frame

    def _hit_frame_get(self):
        return self._hit_frame

    def _hit_frame_set(self, value):
        self._hit_frame = self._prep_frame(value)

    hit_frame = property(fget=_hit_frame_get, fset=_hit_frame_set)

    def _query_frame_get(self):
        return self._query_frame

    def _query_frame_set(self, value):
        self._query_frame = self._prep_frame(value)

    query_frame = property(fget=_query_frame_get, fset=_query_frame_set)

    ## coordinate properties ##
    def _prep_coord(self, coord):
        if not isinstance(coord, int):
            if coord is not None:
                raise ValueError("Coordinate must be an integer or None.")
        return coord

    def _hit_start_get(self):
        return self._hit_start

    def _hit_start_set(self, value):
        self._hit_start = self._prep_coord(value)

    hit_start = property(fget=_hit_start_get, fset=_hit_start_set)

    def _query_start_get(self):
        return self._query_start

    def _query_start_set(self, value):
        self._query_start = self._prep_coord(value)

    query_start = property(fget=_query_start_get, fset=_query_start_set)

    def _hit_end_get(self):
        return self._hit_end

    def _hit_end_set(self, value):
        self._hit_end = self._prep_coord(value)

    hit_end = property(fget=_hit_end_get, fset=_hit_end_set)

    def _query_end_get(self):
        return self._query_end

    def _query_end_set(self, value):
        self._query_end = self._prep_coord(value)

    query_end = property(fget=_query_end_get, fset=_query_end_set)

    ## coordinate-dependent properties ##
    def _hit_span_get(self):
        try:
            return self.hit_end - self.hit_start
        except TypeError:  # triggered if any of the coordinates are None
            return None

    hit_span = property(fget=_hit_span_get)

    def _query_span_get(self):
        try:
            return self.query_end - self.query_start
        except TypeError:  # triggered if any of the coordinates are None
            return None

    query_span = property(fget=_query_span_get)

    def _hit_range_get(self):
        return (self.hit_start, self.hit_end)

    hit_range = property(fget=_hit_range_get)

    def _query_range_get(self):
        return (self.query_start, self.query_end)

    query_range = property(fget=_query_range_get)


def _test():
    """Run the Bio.SearchIO._object module's doctests.

    This will try and locate the unit tests directory, and run the doctests
    from there in order that the relative paths used in the examples work.
    """
    import doctest
    import os

    test_dir = os.path.join('Tests', 'Blast')

    if os.path.isdir(os.path.join('..', '..', test_dir)):
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join('..', '..', test_dir))
        doctest.testmod()
        os.chdir(cur_dir)
        print "Done"


# if not used as a module, run the doctest
if __name__ == "__main__":
    _test()
