# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO objects to model homology search program outputs.

The SearchIO object model is made up of a hierarchy of four nested objects:

    * QueryResult, to represent a search query.

      This is the top-level object returned by the main SearchIO `parse` and
      `read` functions. QueryResult objects may contain zero or more Hit
      objects, each accessible by its ID string (like in Python dictionaries)
      or integer index (like in Python lists).

    * Hit, to represent a database entry containing a full or partial sequence
      match with the query sequence.

      Hit objects contain one or more HSP objects, each accessible by its integer
      index. They behave very similar to a Python list.

    * HSP, to represent a region of significant alignment(s) between the query
      and hit sequences.

      HSP objects contain one or more HSPFragment objects, each accessible by
      its integer index. In most cases, the HSP objects are where the bulk of
      search result statistics (e.g. e-value, bitscore) are stored. Like Hit
      objects, HSPs also behave very similar to a Python list.

    * HSPFragment, to represent a single contiguous alignment between the query
      and hit sequences.

      HSPFragment objects may store hit and query sequences resulting from the
      sequence search. If present, these sequences are stored as SeqRecord
      objects (see SeqRecord). If both of them are present, HSPFragment will
      create a MultipleSeqAlignment object from both sequences.

      Most search programs only have HSPs with one HSPFragment in them, making
      these two objects inseparable. However, there are programs (e.g. BLAT and
      Exonerate) which may have more than one HSPFragment objects in any given
      HSP. If you are not using these programs, you can safely consider HSP and
      HSPFragment as a single union.

"""

import warnings
from copy import deepcopy
from itertools import chain

from Bio import BiopythonWarning
from Bio._py3k import OrderedDict
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import single_letter_alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


__all__ = ['QueryResult', 'Hit', 'HSP', 'HSPFragment']


# helper property functions
def _singleitem(attr=None, doc=''):
    """Returns a property that fetches the given attribute from
    the first item in a SearchIO container object."""
    def getter(self):
        if len(self._items) > 1:
            raise ValueError("More than one HSPFragment objects "
                    "found in HSP")
        if attr is None:
            return self._items[0]
        return getattr(self._items[0], attr)
    return property(fget=getter, doc=doc)


def _allitems(attr=None, doc=''):
    """Returns a property that fetches the given attributes from
    all items in a SearchIO container object."""
    def getter(self):
        if attr is None:
            return self._items
        return [getattr(frag, attr) for frag in self._items]
    return property(fget=getter, doc=doc)


def _partialcascade(cont_attr, item_attr, doc=''):
    """Returns a getter property with a cascading setter.

    This is used for the `id` and `description` properties of the container
    objects. These items have their own private attributes that stores query
    and/or hit ID and description. To keep the container items' query and/or
    hit ID and description in-sync, the setter cascades any new value given
    to the items' values as well.

    """
    def getter(self):
        return getattr(self, cont_attr)

    def setter(self, value):
        setattr(self, cont_attr, value)
        for item in self:
            setattr(item, item_attr, value)

    return property(fget=getter, fset=setter, doc=doc)


def _fullcascade(attr, doc=''):
    """Returns a getter property with a cascading setter.

    This is similar to `_partialcascade`, but for SearchIO containers that have
    at least one item (Hit and HSP). The getter always retrieves the attribute
    value from the first item. If the items have more than one attribute values,
    an error will be raised. The setter behaves like `_partialcascade`.

    """
    def getter(self):
        attrset = set([getattr(item, attr) for item in self._items])
        if len(attrset) > 1:
            raise ValueError("More than one value present in the contained "
                    "items: %r" % list(attrset))
        return getattr(self._items[0], attr)

    def setter(self, value):
        for item in self:
            setattr(item, attr, value)

    return property(fget=getter, fset=setter, doc=doc)


class _BaseSearchObject(object):

    """Abstract class for SearchIO objects."""

    _NON_STICKY_ATTRS = ()

    def _transfer_attrs(self, obj):
        """Transfer instance attributes to the given object.

        This method is used to transfer attributes set externally (for example
        using `setattr`) to a new object created from this one (for example
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

    def _trunc_display(string, max_len, concat_char):
        """Truncates the given string for display."""
        if len(string) > max_len:
            return string[:max_len - len(concat_char)] + concat_char
        return string

    _trunc_display = staticmethod(_trunc_display)

    def _attr_display(obj, attr, fmt=None, fallback='?'):
        """Returns a string of the given object's attribute."""
        if hasattr(obj, attr):
            if fmt is not None:
                return fmt % getattr(obj, attr)
            return str(getattr(obj, attr))
        return fallback

    _attr_display = staticmethod(_attr_display)


class QueryResult(_BaseSearchObject):

    """Class representing search results from a single query.

    QueryResult is the container object that stores all search hits from a
    single search query. It is the top-level object returned by SearchIO's two
    main functions, `read` and `parse`. Depending on the search results and
    search output format, a QueryResult object contains zero or more Hit
    objects (see Hit).

    You can take a quick look at a QueryResult's contents and attributes by
    invoking `print` on it:

    >>> from Bio import SearchIO
    >>> qresult = SearchIO.parse('Blast/mirna.xml', 'blast-xml').next()
    >>> print qresult[:12]      # show the first 11 hits
    Program: blastn (2.2.27+)
      Query: 33211 (61)
             mir_1
     Target: refseq_rna
       Hits: ----  -----  ----------------------------------------------------------
                #  # HSP  ID + description                                          
             ----  -----  ----------------------------------------------------------
                0      1  gi|262205317|ref|NR_030195.1|  Homo sapiens microRNA 52...
                1      1  gi|301171311|ref|NR_035856.1|  Pan troglodytes microRNA...
                2      1  gi|270133242|ref|NR_032573.1|  Macaca mulatta microRNA ...
                3      2  gi|301171322|ref|NR_035857.1|  Pan troglodytes microRNA...
                4      1  gi|301171267|ref|NR_035851.1|  Pan troglodytes microRNA...
                5      2  gi|262205330|ref|NR_030198.1|  Homo sapiens microRNA 52...
                6      1  gi|262205302|ref|NR_030191.1|  Homo sapiens microRNA 51...
                7      1  gi|301171259|ref|NR_035850.1|  Pan troglodytes microRNA...
                8      1  gi|262205451|ref|NR_030222.1|  Homo sapiens microRNA 51...
                9      2  gi|301171447|ref|NR_035871.1|  Pan troglodytes microRNA...
               10      1  gi|301171276|ref|NR_035852.1|  Pan troglodytes microRNA...
               11      1  gi|262205290|ref|NR_030188.1|  Homo sapiens microRNA 51...

    If you just want to know how many hits a QueryResult has, you can invoke
    `len` on it. Alternatively, you can simply type its name in the interpreter:

    >>> len(qresult)
    100
    >>> qresult
    QueryResult(id='33211', 100 hits)

    QueryResult behaves like a hybrid of Python's built-in list and dictionary.
    You can retrieve its items (Hit objects) using the integer index of the
    item, just like regular Python lists:

    >>> first_hit = qresult[0]
    >>> first_hit
    Hit(id='gi|262205317|ref|NR_030195.1|', query_id='33211', 1 hsps)

    You can slice QueryResult objects as well. However, instead of returning a
    list, slicing will return a new QueryResult object containing only the
    sliced hits:

    >>> sliced_qresult = qresult[:3]    # slice the first three hits
    >>> len(qresult)
    100
    >>> len(sliced_qresult)
    3
    >>> print sliced_qresult
    Program: blastn (2.2.27+)
      Query: 33211 (61)
             mir_1
     Target: refseq_rna
       Hits: ----  -----  ----------------------------------------------------------
                #  # HSP  ID + description                                          
             ----  -----  ----------------------------------------------------------
                0      1  gi|262205317|ref|NR_030195.1|  Homo sapiens microRNA 52...
                1      1  gi|301171311|ref|NR_035856.1|  Pan troglodytes microRNA...
                2      1  gi|270133242|ref|NR_032573.1|  Macaca mulatta microRNA ...

    It is also possible to retrieve hits using the hit's ID. This is useful for
    retrieving hits that you know should exist in a given search:

    >>> hit = qresult['gi|262205317|ref|NR_030195.1|']
    >>> hit
    Hit(id='gi|262205317|ref|NR_030195.1|', query_id='33211', 1 hsps)

    You can also replace a Hit in QueryResult with another Hit using either the
    integer index or hit key string. Note that the replacing object must be a
    Hit that has the same `query_id` property as the QueryResult object.

    If you're not sure whether a QueryResult contains a particular hit, you can
    use the hit ID to check for membership first:

    >>> 'gi|262205317|ref|NR_030195.1|' in qresult
    True
    >>> 'gi|262380031|ref|NR_023426.1|' in qresult
    False

    Or, if you just want to know the rank / position of a given hit, you can
    use the hit ID as an argument for the `index` method. Note that the values
    returned will be zero-based. So zero (0) means the hit is the first in the
    QueryResult, three (3) means the hit is the fourth item, and so on. If the
    hit does not exist in the QueryResult, a `ValueError` will be raised.

    >>> qresult.index('gi|262205317|ref|NR_030195.1|')
    0
    >>> qresult.index('gi|262205330|ref|NR_030198.1|')
    5
    >>> qresult.index('gi|262380031|ref|NR_023426.1|')
    Traceback (most recent call last):
    ...
    ValueError: ...

    To ease working with a large number of hits, QueryResult has several filter
    and map methods, analogous to Python's built-in functions with the same
    names. There are filter and map methods available for operations over Hit
    objects or HSP objects. As an example, here we are using the `hit_map`
    method to rename all hit IDs within a QueryResult:

    >>> def renamer(hit):
    ...     hit.id = hit.id.split('|')[3]
    ...     return hit
    >>> mapped_qresult = qresult.hit_map(renamer)
    >>> print mapped_qresult[:3]
    Program: blastn (2.2.27+)
      Query: 33211 (61)
             mir_1
     Target: refseq_rna
       Hits: ----  -----  ----------------------------------------------------------
                #  # HSP  ID + description                                          
             ----  -----  ----------------------------------------------------------
                0      1  NR_030195.1  Homo sapiens microRNA 520b (MIR520B), micr...
                1      1  NR_035856.1  Pan troglodytes microRNA mir-520b (MIR520B...
                2      1  NR_032573.1  Macaca mulatta microRNA mir-519a (MIR519A)...

    The principle for other map and filter methods are similar: they take a
    function, applies it, and returns a new QueryResult object.

    There are also other methods useful for working with list-like objects:
    `append`, `pop`, and `sort`. More details and examples are available in
    their respective documentations.

    Finally, just like Python lists and dictionaries, QueryResult objects are
    iterable. Iteration over QueryResults will yield Hit objects:

    >>> for hit in qresult[:4]:     # iterate over the first four items
    ...     hit
    ...
    Hit(id='gi|262205317|ref|NR_030195.1|', query_id='33211', 1 hsps)
    Hit(id='gi|301171311|ref|NR_035856.1|', query_id='33211', 1 hsps)
    Hit(id='gi|270133242|ref|NR_032573.1|', query_id='33211', 1 hsps)
    Hit(id='gi|301171322|ref|NR_035857.1|', query_id='33211', 2 hsps)

    If you need access to all the hits in a QueryResult object, you can access
    them in a list using the `hits` property. Similary, access to all hit IDs is
    available through the `hit_keys` property.

    >>> qresult.hits
    [Hit(id='gi|262205317|ref|NR_030195.1|', query_id='33211', 1 hsps), ...]
    >>> qresult.hit_keys
    ['gi|262205317|ref|NR_030195.1|', 'gi|301171311|ref|NR_035856.1|', ...]

    """

    # attributes we don't want to transfer when creating a new QueryResult class
    # from this one
    _NON_STICKY_ATTRS = ('_items',)

    def __init__(self, id='<unknown id>', hits=[], \
            hit_key_function=lambda hit: hit.id):
        """Initializes a QueryResult object.

        Arguments:
        id -- String of query sequence ID.
        hits -- Iterator returning Hit objects.
        hit_key_function -- Function to define hit keys, defaults to a function
                            that return Hit object IDs.

        """
        if id is None:
            raise ValueError("Query ID string is required for QueryResult "
                    "creation")

        self._id = id
        self._hit_key_function = hit_key_function
        self._items = OrderedDict()
        self._description = '<unknown description>'
        self.program = '<unknown program>'
        self.target = '<unknown target>'
        self.version = '<unknown version>'

        # validate Hit objects and fill up self._items
        for hit in hits:
            # validation is handled by __setitem__
            self.append(hit)

    # handle Python 2 OrderedDict behavior
    if hasattr(OrderedDict, 'iteritems'):

        def __iter__(self):
            return iter(self.iterhits())

        @property
        def hits(self):
            """Hit objects contained in the QueryResult."""
            return self._items.values()

        @property
        def hit_keys(self):
            """Hit IDs of the Hit objects contained in the QueryResult."""
            return self._items.keys()

        @property
        def items(self):
            """List of tuples of Hit IDs and Hit objects."""
            return self._items.items()

        def iterhits(self):
            """Returns an iterator over the Hit objects."""
            for hit in self._items.itervalues():
                yield hit

        def iterhit_keys(self):
            """Returns an iterator over the ID of the Hit objects."""
            for hit_id in self._items.iterkeys():
                yield hit_id

        def iteritems(self):
            """Returns an iterator yielding tuples of Hit ID and Hit objects."""
            for item in self._items.iteritems():
                yield item

    else:

        def __iter__(self):
            return iter(self.hits)

        @property
        def hits(self):
            """Hit objects contained in the QueryResult."""
            return list(self._items.values())

        @property
        def hit_keys(self):
            """Hit IDs of the Hit objects contained in the QueryResult."""
            return list(self._items.keys())

        @property
        def items(self):
            """List of tuples of Hit IDs and Hit objects."""
            return list(self._items.items())

        def iterhits(self):
            """Returns an iterator over the Hit objects."""
            for hit in self._items.values():
                yield hit

        def iterhit_keys(self):
            """Returns an iterator over the ID of the Hit objects."""
            for hit_id in self._items.keys():
                yield hit_id

        def iteritems(self):
            """Returns an iterator yielding tuples of Hit ID and Hit objects."""
            for item in self._items.items():
                yield item

    def __contains__(self, hit_key):
        if isinstance(hit_key, Hit):
            return self._hit_key_function(hit_key) in self._items
        return hit_key in self._items

    def __len__(self):
        return len(self._items)

    def __nonzero__(self):
        return bool(self._items)

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
            qid_line += QueryResult._trunc_display('\n         %s' %
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
            lines.append(pattern % ('#', '# HSP',
                'ID + description'.ljust(58)))
            lines.append(pattern % ('-'*4, '-'*5, '-'*58))
            for idx, hit in enumerate(self.hits):
                if idx < 30:
                    hid_line = '%s  %s' % (hit.id, hit.description)
                    if len(hid_line) > 58:
                        hid_line = hid_line[:55] + '...'
                    lines.append(pattern % (idx, str(len(hit)),
                        hid_line.ljust(58)))
                elif idx > len(self.hits) - 4:
                    hid_line = '%s  %s' % (hit.id, hit.description)
                    if len(hid_line) > 58:
                        hid_line = hid_line[:55] + '...'
                    lines.append(pattern % (idx, str(len(hit)),
                        hid_line.ljust(58)))
                elif idx == 30:
                    lines.append('%14s' % '~~~')

        return '\n'.join(lines)

    def __getitem__(self, hit_key):
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
            return list(self._items.values())[hit_key]

        # if key is a string, then do a regular dictionary retrieval
        return self._items[hit_key]

    def __setitem__(self, hit_key, hit):
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

        self._items[hit_key] = hit

    def __delitem__(self, hit_key):
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
            del self._items[key]
        return

    ## properties ##
    id = _partialcascade('_id', 'query_id', """QueryResult ID string""")
    description = _partialcascade('_description', 'query_description',
            """QueryResult description""")

    ## public methods ##
    def absorb(self, hit):
        """Adds a Hit object to the end of QueryResult. If the QueryResult
        already has a Hit with the same ID, append the new Hit's HSPs into
        the existing Hit.

        Arguments:
        hit -- Hit object to absorb.

        This method is used for file formats that may output the same Hit in
        separate places, such as BLAT or Exonerate. In both formats, Hit
        with different strands are put in different places. However, SearchIO
        considers them to be the same as a Hit object should be all database
        entries with the same ID, regardless of strand orientation.

        """
        try:
            self.append(hit)
        except ValueError:
            assert hit.id in self
            for hsp in hit:
                self[hit.id].append(hsp)

    def append(self, hit):
        """Adds a Hit object to the end of QueryResult.

        Parameters
        hit -- Hit object to append.

        Any Hit object appended must have the same `query_id` property as the
        QueryResult's `id` property. If the hit key already exists, a
        `ValueError` will be raised.

        """
        # if a custom hit_key_function is supplied, use it to define th hit key
        if self._hit_key_function is not None:
            hit_key = self._hit_key_function(hit)
        else:
            hit_key = hit.id

        if hit_key not in self:
            self[hit_key] = hit
        else:
            raise ValueError("Hit '%s' already present in this QueryResult." %
                    hit_key)

    def hit_filter(self, func=None):
        """Creates a new QueryResult object whose Hit objects pass the filter
        function.

        Arguments:
        func -- Callback function that accepts a Hit object as its parameter,
                does a boolean check, and returns True or False

        Here is an example of using `hit_filter` to select Hits whose
        description begins with the string 'Homo sapiens', case sensitive:

        >>> from Bio import SearchIO
        >>> qresult = SearchIO.parse('Blast/mirna.xml', 'blast-xml').next()
        >>> def desc_filter(hit):
        ...     return hit.description.startswith('Homo sapiens')
        ...
        >>> len(qresult)
        100
        >>> filtered = qresult.hit_filter(desc_filter)
        >>> len(filtered)
        39
        >>> print filtered[:4]
        Program: blastn (2.2.27+)
          Query: 33211 (61)
                 mir_1
         Target: refseq_rna
           Hits: ----  -----  ----------------------------------------------------------
                    #  # HSP  ID + description                                          
                 ----  -----  ----------------------------------------------------------
                    0      1  gi|262205317|ref|NR_030195.1|  Homo sapiens microRNA 52...
                    1      2  gi|262205330|ref|NR_030198.1|  Homo sapiens microRNA 52...
                    2      1  gi|262205302|ref|NR_030191.1|  Homo sapiens microRNA 51...
                    3      1  gi|262205451|ref|NR_030222.1|  Homo sapiens microRNA 51...

        Note that instance attributes (other than the hits) from the unfiltered
        QueryResult are retained in the filtered object.

            >>> qresult.program == filtered.program
            True
            >>> qresult.target == filtered.target
            True

        """
        hits = filter(func, self.hits)
        obj =  self.__class__(self.id, hits, self._hit_key_function)
        self._transfer_attrs(obj)
        return obj

    def hit_map(self, func=None):
        """Creates a new QueryResult object, mapping the given function to its
        Hits.

        Arguments:
        func -- Callback function that accepts a Hit object as its parameter and
                also returns a Hit object.

        Here is an example of using `hit_map` with a function that discards all
        HSPs in a Hit except for the first one:

        >>> from Bio import SearchIO
        >>> qresult = SearchIO.parse('Blast/mirna.xml', 'blast-xml').next()
        >>> print qresult[:8]
        Program: blastn (2.2.27+)
          Query: 33211 (61)
                 mir_1
         Target: refseq_rna
           Hits: ----  -----  ----------------------------------------------------------
                    #  # HSP  ID + description                                          
                 ----  -----  ----------------------------------------------------------
                    0      1  gi|262205317|ref|NR_030195.1|  Homo sapiens microRNA 52...
                    1      1  gi|301171311|ref|NR_035856.1|  Pan troglodytes microRNA...
                    2      1  gi|270133242|ref|NR_032573.1|  Macaca mulatta microRNA ...
                    3      2  gi|301171322|ref|NR_035857.1|  Pan troglodytes microRNA...
                    4      1  gi|301171267|ref|NR_035851.1|  Pan troglodytes microRNA...
                    5      2  gi|262205330|ref|NR_030198.1|  Homo sapiens microRNA 52...
                    6      1  gi|262205302|ref|NR_030191.1|  Homo sapiens microRNA 51...
                    7      1  gi|301171259|ref|NR_035850.1|  Pan troglodytes microRNA...

        >>> top_hsp = lambda hit: hit[:1]
        >>> mapped_qresult = qresult.hit_map(top_hsp)
        >>> print mapped_qresult[:8]
        Program: blastn (2.2.27+)
          Query: 33211 (61)
                 mir_1
         Target: refseq_rna
           Hits: ----  -----  ----------------------------------------------------------
                    #  # HSP  ID + description                                          
                 ----  -----  ----------------------------------------------------------
                    0      1  gi|262205317|ref|NR_030195.1|  Homo sapiens microRNA 52...
                    1      1  gi|301171311|ref|NR_035856.1|  Pan troglodytes microRNA...
                    2      1  gi|270133242|ref|NR_032573.1|  Macaca mulatta microRNA ...
                    3      1  gi|301171322|ref|NR_035857.1|  Pan troglodytes microRNA...
                    4      1  gi|301171267|ref|NR_035851.1|  Pan troglodytes microRNA...
                    5      1  gi|262205330|ref|NR_030198.1|  Homo sapiens microRNA 52...
                    6      1  gi|262205302|ref|NR_030191.1|  Homo sapiens microRNA 51...
                    7      1  gi|301171259|ref|NR_035850.1|  Pan troglodytes microRNA...

        """
        hits = [deepcopy(hit) for hit in self.hits]
        if func is not None:
            hits = map(func, hits)
        obj =  self.__class__(self.id, hits, self._hit_key_function)
        self._transfer_attrs(obj)
        return obj

    def hsp_filter(self, func=None):
        """Creates a new QueryResult object whose HSP objects pass the filter
        function.

        `hsp_filter` is the same as `hit_filter`, except that it filters
        directly on each HSP object in every Hit. If a the filtering removes
        all HSP object in a given Hit, the entire Hit will be discarded. This
        will result in the QueryResult having less Hit after filtering.

        """
        hits = filter(None, (hit.filter(func) for hit in self.hits))
        obj =  self.__class__(self.id, hits, self._hit_key_function)
        self._transfer_attrs(obj)
        return obj

    def hsp_map(self, func=None):
        """Creates a new QueryResult object, mapping the given function to its
        HSPs.

        `hsp_map` is the same as `hit_map`, except that it applies the given
        function to all HSP objects in every Hit, instead of the Hit objects.

        """
        hits = filter(None, (hit.map(func) for hit in list(self.hits)[:]))
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

        By default, `pop` will remove and return the last Hit object in the
        QueryResult object. To remove specific Hit objects, you can use its
        integer index or hit key.

        >>> from Bio import SearchIO
        >>> qresult = SearchIO.parse('Blast/mirna.xml', 'blast-xml').next()
        >>> len(qresult)
        100
        >>> for hit in qresult[:5]:
        ...     print hit.id
        ...
        gi|262205317|ref|NR_030195.1|
        gi|301171311|ref|NR_035856.1|
        gi|270133242|ref|NR_032573.1|
        gi|301171322|ref|NR_035857.1|
        gi|301171267|ref|NR_035851.1|

        # remove the last hit
        >>> qresult.pop()
        Hit(id='gi|397513516|ref|XM_003827011.1|', query_id='33211', 1 hsps)

        # remove the first hit
        >>> qresult.pop(0)
        Hit(id='gi|262205317|ref|NR_030195.1|', query_id='33211', 1 hsps)

        # remove hit with the given ID
        >>> qresult.pop('gi|301171322|ref|NR_035857.1|')
        Hit(id='gi|301171322|ref|NR_035857.1|', query_id='33211', 2 hsps)

        """
        # if key is an integer (index)
        # get the ID for the Hit object at that index
        if isinstance(hit_key, int):
            # raise the appropriate error if there is no hit
            if not self:
                raise IndexError("pop from empty list")
            hit_key = list(self.hit_keys)[hit_key]

        try:
            return self._items.pop(hit_key)
        except KeyError:
            # if key doesn't exist and no default is set, raise a KeyError
            if default is self.__marker:
                raise KeyError(hit_key)
        # if key doesn't exist but a default is set, return the default value
        return default

    def index(self, hit_key):
        """Returns the index of a given hit key, zero-based.

        Arguments:
        hit_key -- Hit ID string to look up.

        This method is useful for finding out the integer index (usually
        correlated with search rank) of a given hit key.

        >>> from Bio import SearchIO
        >>> qresult = SearchIO.parse('Blast/mirna.xml', 'blast-xml').next()
        >>> qresult.index('gi|301171259|ref|NR_035850.1|')
        7

        """
        if isinstance(hit_key, Hit):
            return list(self.hit_keys).index(hit_key.id)
        return list(self.hit_keys).index(hit_key)

    def sort(self, key=None, reverse=False, in_place=True):
        # no cmp argument to make sort more Python 3-like
        """Sorts the Hit objects.

        Arguments:
        key -- Function used to sort the Hit objects.
        reverse -- Boolean, whether to reverse the sorting or not.
        in_place -- Boolean, whether to perform sorting in place (in the same
                    object) or not (creating a new object).

        `sort` defaults to sorting in-place, to mimick Python's `list.sort`
        method. If you set the `in_place` argument to False, it will treat
        return a new, sorted QueryResult object and keep the initial one
        unsorted.

        """
        if key is None:
            # if reverse is True, reverse the hits
            if reverse:
                sorted_hits = list(self.hits)[::-1]
            # otherwise (default options) make a copy of the hits
            else:
                sorted_hits = list(self.hits)[:]
        else:
            sorted_hits = sorted(self.hits, key=key, reverse=reverse)

        # if sorting is in-place, don't create a new QueryResult object
        if in_place:
            new_hits = OrderedDict()
            for hit in sorted_hits:
                new_hits[self._hit_key_function(hit)] = hit
            self._items = new_hits
        # otherwise, return a new sorted QueryResult object
        else:
            obj =  self.__class__(self.id, sorted_hits, self._hit_key_function)
            self._transfer_attrs(obj)
            return obj


class Hit(_BaseSearchObject):

    """Class representing a single database hit of a search result.

    Hit objects are the second-level container in the SearchIO module. They
    are the objects contained within a QueryResult (see QueryResult). They
    themselves are container for HSP objects and will contain at least one
    HSP.

    To have a quick look at a Hit and its contents, invoke `print` on it:

    >>> from Bio import SearchIO
    >>> qresult = SearchIO.parse('Blast/mirna.xml', 'blast-xml').next()
    >>> hit = qresult[3]
    >>> print hit
    Query: 33211
           mir_1
      Hit: gi|301171322|ref|NR_035857.1| (86)
           Pan troglodytes microRNA mir-520c (MIR520C), microRNA
     HSPs: ----  --------  ---------  ------  ---------------  ---------------------
              #   E-value  Bit score    Span      Query range              Hit range
           ----  --------  ---------  ------  ---------------  ---------------------
              0   8.9e-20     100.47      60           [1:61]                [13:73]
              1   3.3e-06      55.39      60           [0:60]                [13:73]

    You can invoke `len` on a Hit object to see how many HSP objects it contains:

    >>> len(hit)
    2

    Hit objects behave very similar to Python lists. You can retrieve the HSP
    object inside a Hit using the HSP's integer index. Hit objects can also be
    sliced, which will return a new Hit objects containing only the sliced HSPs:

    # HSP items inside the Hit can be retrieved using its integer index
    >>> hit[0]
    HSP(hit_id='gi|301171322|ref|NR_035857.1|', query_id='33211', 1 fragments)

    # slicing returns a new Hit
    >>> hit
    Hit(id='gi|301171322|ref|NR_035857.1|', query_id='33211', 2 hsps)
    >>> hit[:1]
    Hit(id='gi|301171322|ref|NR_035857.1|', query_id='33211', 1 hsps)
    >>> print hit[1:]
    Query: 33211
           mir_1
      Hit: gi|301171322|ref|NR_035857.1| (86)
           Pan troglodytes microRNA mir-520c (MIR520C), microRNA
     HSPs: ----  --------  ---------  ------  ---------------  ---------------------
              #   E-value  Bit score    Span      Query range              Hit range
           ----  --------  ---------  ------  ---------------  ---------------------
              0   3.3e-06      55.39      60           [0:60]                [13:73]

    Hit objects provide `filter` and `map` methods, which are analogous to
    Python's built-in `filter` and `map` except that they return a new Hit
    object instead of a list.

    Here is an example of using `filter` to select for HSPs whose e-value is
    less than 1e-10:

    >>> evalue_filter = lambda hsp: hsp.evalue < 1e-10
    >>> filtered_hit = hit.filter(evalue_filter)
    >>> len(hit)
    2
    >>> len(filtered_hit)
    1
    >>> print filtered_hit
    Query: 33211
           mir_1
      Hit: gi|301171322|ref|NR_035857.1| (86)
           Pan troglodytes microRNA mir-520c (MIR520C), microRNA
     HSPs: ----  --------  ---------  ------  ---------------  ---------------------
              #   E-value  Bit score    Span      Query range              Hit range
           ----  --------  ---------  ------  ---------------  ---------------------
              0   8.9e-20     100.47      60           [1:61]                [13:73]

    There are also other methods which are counterparts of Python lists' methods
    with the same names: `append`, `index`, `pop`, and `sort`. Consult their
    respective documentations for more details and examples of their usage.

    """

    # attributes we don't want to transfer when creating a new Hit class
    # from this one
    _NON_STICKY_ATTRS = ('_items', )

    def __init__(self, hsps=[]):
        """Initializes a Hit object.

        Arguments:
        hsps -- List containing HSP objects.

        Hit objects must be initialized with a list containing at least one HSP
        object. If multiple HSP objects are used for initialization, they must
        all have the same `query_id`, `query_description`, `hit_id`, and
        `hit_description` properties.

        """
        if not hsps:
            raise ValueError("Hit objects must have at least one HSP object.")
        # check that all fragments contain the same IDs, descriptions
        for attr in ('query_id', 'query_description', 'hit_id', \
                'hit_description'):
            if len(set([getattr(hsp, attr) for hsp in hsps])) != 1:
                raise ValueError("Hit object can not contain HSPs with "
                        "more than one %s." % attr)

        self._items = []
        for hsp in hsps:
            # validate each HSP
            self._validate_hsp(hsp)
            # and store it them as an instance attribute
            self.append(hsp)

    def __repr__(self):
        return "Hit(id=%r, query_id=%r, %r hsps)" % (self.id, self.query_id,
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
        qid_line = 'Query: %s' % self.query_id
        if self.query_description:
            qid_line += Hit._trunc_display('\n       %s' %
                    self.query_description, 80, '...')
        lines.append(qid_line)

        # set hit id line
        hid_line = '  Hit: %s' % self.id
        if hasattr(self, 'seq_len'):
            hid_line += ' (%i)' % self.seq_len
        if self.description:
            hid_line += Hit._trunc_display('\n       %s' % self.description,
                    80, '...')
        lines.append(hid_line)

        # set hsp line and table
        if not self.hsps:
            lines.append(' HSPs: ?')
        else:
            lines.append(' HSPs: %s  %s  %s  %s  %s  %s' % \
                    ('-'*4, '-'*8, '-'*9, '-'*6, '-'*15, '-'*21))
            pattern = '%11s  %8s  %9s  %6s  %15s  %21s'
            lines.append(pattern % ('#', 'E-value', 'Bit score', 'Span',
                    'Query range', 'Hit range'))
            lines.append(pattern % ('-'*4, '-'*8, '-'*9, '-'*6, '-'*15, '-'*21))
            for idx, hsp in enumerate(self.hsps):
                # evalue
                evalue = Hit._attr_display(hsp, 'evalue', fmt='%.2g')
                # bitscore
                bitscore = Hit._attr_display(hsp, 'bitscore', fmt='%.2f')
                # alignment length
                aln_span = Hit._attr_display(hsp, 'aln_span')
                # query region
                query_start = Hit._attr_display(hsp, 'query_start')
                query_end = Hit._attr_display(hsp, 'query_end')
                query_range = '[%s:%s]' % (query_start, query_end)
                # max column length is 18
                query_range = Hit._trunc_display(query_range, 15, '~]')
                # hit region
                hit_start = Hit._attr_display(hsp, 'hit_start')
                hit_end = Hit._attr_display(hsp, 'hit_end')
                hit_range = '[%s:%s]' % (hit_start, hit_end)
                hit_range = Hit._trunc_display(hit_range, 21, '~]')
                # append the hsp row
                lines.append(pattern % (str(idx), evalue, bitscore, aln_span,
                        query_range, hit_range))

        return '\n'.join(lines)

    def __getitem__(self, idx):
        # if key is slice, return a new Hit instance
        if isinstance(idx, slice):
            obj = self.__class__(self.hsps[idx])
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
        # HACK: to make validation during __init__ work
        if self._items:
            if hsp.hit_id != self.id:
                raise ValueError("Expected HSP with hit ID %r, " \
                        "found %r instead." % (self.id, hsp.hit_id))
            if hsp.query_id != self.query_id:
                raise ValueError("Expected HSP with query ID %r, " \
                        "found %r instead." % (self.query_id, hsp.query_id))

    ## properties ##
    description = _fullcascade('hit_description', """Hit description""")
    query_description = _fullcascade('query_description',
            """Description of the query that produced the hit""")
    id = _fullcascade('hit_id', """Hit ID string.""")
    query_id = _fullcascade('query_id',
            """ID string of the query that produced the hit""")
    # returns all hsps
    hsps = _allitems(doc="""HSP objects contained in the Hit""")
    # returns all fragments
    fragments = property(lambda self: list(chain(*self._items)), \
            doc="""HSPFragment objects contained in the Hit""")

    ## public methods ##
    def append(self, hsp):
        """Adds a HSP object to the end of Hit.

        Parameters
        hsp -- HSP object to append.

        Any HSP object appended must have the same `hit_id` property as the
        Hit object's `id` property and the same `query_id` property as the
        Hit object's `query_id` property.

        """
        self._validate_hsp(hsp)
        self._items.append(hsp)

    def filter(self, func=None):
        """Creates a new Hit object whose HSP objects pass the filter
        function.

        Arguments:
        func -- Callback function that accepts a HSP object as its parameter,
                does a boolean check, and returns True or False.

        `filter` is analogous to Python's built-in `filter` function, except
        that instead of returning a list it returns a `Hit` object. Here is an
        example of using `filter` to select for HSPs having bitscores bigger
        than 60:

        >>> from Bio import SearchIO
        >>> qresult = SearchIO.parse('Blast/mirna.xml', 'blast-xml').next()
        >>> hit = qresult[3]
        >>> evalue_filter = lambda hsp: hsp.bitscore > 60
        >>> filtered_hit = hit.filter(evalue_filter)
        >>> len(hit)
        2
        >>> len(filtered_hit)
        1
        >>> print filtered_hit
        Query: 33211
               mir_1
          Hit: gi|301171322|ref|NR_035857.1| (86)
               Pan troglodytes microRNA mir-520c (MIR520C), microRNA
         HSPs: ----  --------  ---------  ------  ---------------  ---------------------
                  #   E-value  Bit score    Span      Query range              Hit range
               ----  --------  ---------  ------  ---------------  ---------------------
                  0   8.9e-20     100.47      60           [1:61]                [13:73]

        """
        hsps = filter(func, self.hsps)
        if hsps:
            obj = self.__class__(hsps)
            self._transfer_attrs(obj)
            return obj

    def index(self, hsp):
        """Returns the index of a given HSP object, zero-based.

        Arguments:
        hsp -- HSP object to be looked up.

        """
        return self._items.index(hsp)

    def map(self, func=None):
        """Creates a new Hit object, mapping the given function to its HSPs.

        Arguments:
        func -- Callback function that accepts a HSP object as its parameter and
                also returns a HSP object.

        `map` is analogous to Python's built-in `map` function. It is applied to
        all HSPs contained in the Hit object and returns a new Hit object.

        """
        if func is not None:
            hsps = map(func, self.hsps[:])  # this creates a shallow copy
        else:
            hsps = self.hsps[:]
        if hsps:
            obj = self.__class__(hsps)
            self._transfer_attrs(obj)
            return obj

    def pop(self, index=-1):
        """Removes and returns the HSP object at the specified index.

        Arguments:
        index -- Integer denoting the index of the HSP object to remove.

        """
        return self._items.pop(index)

    def sort(self, key=None, reverse=False, in_place=True):
        """Sorts the HSP objects.

        Arguments:
        key -- Function used to sort the HSP objects.
        reverse -- Boolean, whether to reverse the sorting or not.
        in_place -- Boolean, whether to perform sorting in place (in the same
                    object) or not (creating a new object).

        `sort` defaults to sorting in-place, to mimick Python's `list.sort`
        method. If you set the `in_place` argument to False, it will treat
        return a new, sorted Hit object and keep the initial one unsorted

        """
        if in_place:
            self._items.sort(key=key, reverse=reverse)
        else:
            hsps = self.hsps[:]
            hsps.sort(key=key, reverse=reverse)
            obj = self.__class__(hsps)
            self._transfer_attrs(obj)
            return obj


class _BaseHSP(_BaseSearchObject):

    """Abstract base class for HSP objects."""

    def _str_hsp_header(self):
        """Prints the alignment header info."""
        lines = []
        # set query id line
        qid_line = self._trunc_display('      Query: %s %s' %
                (self.query_id, self.query_description), 80, '...')
        # set hit id line
        hid_line = self._trunc_display('        Hit: %s %s' %
                (self.hit_id, self.hit_description), 80, '...')
        lines.append(qid_line)
        lines.append(hid_line)

        # coordinates
        query_start = _BaseHSP._attr_display(self, 'query_start')
        query_end = _BaseHSP._attr_display(self, 'query_end')
        hit_start = _BaseHSP._attr_display(self, 'hit_start')
        hit_end = _BaseHSP._attr_display(self, 'hit_end')

        # strands
        try:
            qstrand = self.query_strand
            hstrand = self.hit_strand
        except ValueError:
            qstrand = self.query_strands[0]
            hstrand = self.hit_strands[0]
        lines.append('Query range: [%s:%s] (%r)' % (query_start, query_end,
                qstrand))
        lines.append('  Hit range: [%s:%s] (%r)' % (hit_start, hit_end,
                hstrand))

        return '\n'.join(lines)


class HSP(_BaseHSP):

    """Class representing high-scoring region(s) between query and hit.

    HSP (high-scoring pair) objects are contained by Hit objects (see Hit).
    In most cases, HSP objects store the bulk of the statistics and results
    (e.g. e-value, bitscores, query sequence, etc.) produced by a search
    program.

    Depending on the search output file format, a given HSP will contain one
    or more HSPFragment object(s). Examples of search programs that produce HSP
    with one HSPFragments are BLAST, HMMER, and FASTA. Other programs such as
    BLAT or Exonerate may produce HSPs containing more than one HSPFragment.
    However, their native terminologies may differ: in BLAT these fragments
    are called 'blocks' while in in Exonerate they are called exons or NER.

    Here are examples from each type of HSP. The first one comes from a BLAST
    search:

    >>> from Bio import SearchIO
    >>> blast_qresult = SearchIO.parse('Blast/mirna.xml', 'blast-xml').next()
    >>> blast_hsp = blast_qresult[1][0]     # the first HSP from the second hit
    >>> blast_hsp
    HSP(hit_id='gi|301171311|ref|NR_035856.1|', query_id='33211', 1 fragments)
    >>> print blast_hsp
          Query: 33211 mir_1
            Hit: gi|301171311|ref|NR_035856.1| Pan troglodytes microRNA mir-520b ...
    Query range: [1:61] (1)
      Hit range: [0:60] (1)
    Quick stats: evalue 1.7e-22; bitscore 109.49
      Fragments: 1 (60 columns)
         Query - CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG
                 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
           Hit - CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG

    For HSPs with a single HSPFragment, you can invoke `print` on it and see the
    underlying sequence alignment, if it exists. This is not the case for HSPs
    with more than one HSPFragment. Below is an example, using an HSP from a
    BLAT search. Invoking `print` on these HSPs will instead show a table of the
    HSPFragment objects it contains:

    >>> blat_qresult = SearchIO.read('Blat/mirna.pslx', 'blat-psl', pslx=True)
    >>> blat_hsp = blat_qresult[1][0]       # the first HSP from the second hit
    >>> blat_hsp
    HSP(hit_id='chr11', query_id='blat_1', 2 fragments)
    >>> print blat_hsp
          Query: blat_1 <unknown description>
            Hit: chr11 <unknown description>
    Query range: [42:67] (-1)
      Hit range: [59018929:59018955] (1)
    Quick stats: evalue ?; bitscore ?
      Fragments: ---  --------------  ----------------------  ----------------------
                   #            Span             Query range               Hit range
                 ---  --------------  ----------------------  ----------------------
                   0               6                 [61:67]     [59018929:59018935]
                   1              16                 [42:58]     [59018939:59018955]

    Notice that in HSPs with more than one HSPFragments, the HSP's `query_range`
    `hit_range` properties encompasses all fragments it contains.

    You can check whether an HSP has more than one HSPFragments or not using the
    `is_fragmented` property:

    >>> blast_hsp.is_fragmented
    False
    >>> blat_hsp.is_fragmented
    True

    Since HSP objects are also containers similar to Python lists, you can
    access a single fragment in an HSP using its integer index:

    >>> blat_fragment = blat_hsp[0]
    >>> print blat_fragment
          Query: blat_1 <unknown description>
            Hit: chr11 <unknown description>
    Query range: [61:67] (-1)
      Hit range: [59018929:59018935] (1)
      Fragments: 1 (6 columns)
         Query - tatagt
           Hit - tatagt

    This applies to HSPs objects with a single fragment as well:

    >>> blast_fragment = blast_hsp[0]
    >>> print blast_fragment
          Query: 33211 mir_1
            Hit: gi|301171311|ref|NR_035856.1| Pan troglodytes microRNA mir-520b ...
    Query range: [1:61] (1)
      Hit range: [0:60] (1)
      Fragments: 1 (60 columns)
         Query - CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG
                 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
           Hit - CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG

    Regardless of the search output file format, HSP objects provide the
    properties listed below. These properties always return values in a list,
    due to the HSP object itself being a list-like container. However, for
    HSP objects with a single HSPFragment, shortcut properties that fetches
    the item from the list are also provided.

    ---------------------  --------------------  -------------------------------
    Property               Shortcut              Value
    ---------------------  --------------------  -------------------------------
    alignments             alignment             HSP alignments as
                                                 MultipleSeqAlignment object

    alignment_annotations  alignment_annotation  Dictionary of annotation(s) of
                                                 all fragments' alignments

    fragments              fragment              HSPFragment objects

    hits                   hit                   Hit sequence as SeqRecord
                                                 objects

    hit_starts             hit_start*            Start coordinates of the hit
                                                 fragments

    hit_ends               hit_end*              End coordinates of the hit
                                                 fragments

    hit_spans              hit_span*             Sizes of each hit fragments

    hit_strands            hit_strand            Strand orientations of the hit
                                                 fragments

    hit_frames             hit_frame             Reading frames of the hit
                                                 fragments

    hit_ranges             hit_range             Tuples of start and end
                                                 coordinates of each hit
                                                 fragment

    queries                query                 Query sequence as SeqRecord
                                                 object

    query_starts           query_start*          Start coordinates of the query
                                                 fragments

    query_ends             query_end*            End coordinates of the query
                                                 fragments

    query_spans            query_span*           Sizes of each query fragments

    query_strands          query_strand          Strand orientations of the
                                                 query fragments

    query_frames           query_frame           Reading frames of the query
                                                 fragments

    query_ranges           query_range           Tuples of start and end
                                                 coordinates of each query
                                                 fragment
    ----------------------------------------------------------------------------
    * may be used in HSPs with multiple fragments

    For all types of HSP objects, the property will return the values in a list.
    Shorcuts are only applicable for HSPs with one fragment. Except the ones
    noted, if they are used on an HSP with more than one fragments, an exception
    will be raised.

    For properties that may be used in HSPs with multiple or single fragments
    (`*_start`, `*_end`, and `*_span` properties), their interpretation depends
    on how many fragment the HSP has:

    -----------  ------------------------------------------------
    Property     Value
    -----------  ------------------------------------------------
    hit_start    Smallest coordinate value of all hit fragments

    hit_end      Largest coordinate value of all hit fragments

    hit_span     Difference between `hit_start` and `hit_end`

    query_start  Smallest coordinate value of all query fragments

    query_end    Largest coordinate value of all query fragments

    query_span   Difference between `query_start` and `query_end`
    --------------------------------------------------------------

    In addition to the objects listed above, HSP objects also provide the
    following properties:

    ------------------  --------------------------------------------------------
    Property            Value
    ------------------  --------------------------------------------------------
    aln_span            Total number of residues in all HSPFragment objects

    alphabet            Alphabet used in hit and query SeqRecord objects

    is_fragmented       Boolean, whether the HSP has multiple fragments or not

    hit_id              ID of the hit sequence

    hit_description     Description of the hit sequence

    hit_inter_ranges    List of hit sequence coordinates of the regions between
                        fragments

    query_id            ID of the query sequence

    query_description   Description of the query sequence

    query_inter_ranges  List of query sequence coordinates of the regions
                        between fragments
    ------------------  --------------------------------------------------------

    """
    # attributes we don't want to transfer when creating a new Hit class
    # from this one
    _NON_STICKY_ATTRS = ('_items', '_aln_span')

    def __init__(self, fragments=[]):
        """Initializes an HSP object.

        Arguments:
        fragments -- List of HSPFragment objects.

        HSP objects must be initialized with a list containing at least one
        HSPFragment object. If multiple HSPFragment objects are used for
        initialization, they must all have the same `query_id`,
        `query_description`, `hit_id`, `hit_description`, and alphabet
        properties.

        """
        if not fragments:
            raise ValueError("HSP objects must have at least one HSPFragment "
                    "object.")
        # check that all fragments contain the same IDs, descriptions, alphabet
        for attr in ('query_id', 'query_description', 'hit_id',
                'hit_description', 'alphabet'):
            if len(set([getattr(frag, attr) for frag in fragments])) != 1:
                raise ValueError("HSP object can not contain fragments with "
                        "more than one %s." % attr)

        self._items = []
        for fragment in fragments:
            self._validate_fragment(fragment)
            self._items.append(fragment)

    def __repr__(self):
        return "%s(hit_id=%r, query_id=%r, %r fragments)" % \
                (self.__class__.__name__, self.hit_id, self.query_id, len(self))

    def __iter__(self):
        return iter(self._items)

    def __contains__(self, fragment):
        return fragment in self._items

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
            return '\n'.join([self._str_hsp_header(), '\n'.join(lines),
                    self.fragments[0]._str_aln()])
        else:
            lines.append('  Fragments: %s  %s  %s  %s' %
                    ('-'*3, '-'*14, '-'*22, '-'*22))
            pattern = '%16s  %14s  %22s  %22s'
            lines.append(pattern % ('#', 'Span', 'Query range', 'Hit range'))
            lines.append(pattern % ('-'*3, '-'*14, '-'*22, '-'*22))
            for idx, block in enumerate(self.fragments):
                # set hsp line and table
                # alignment span
                aln_span = HSP._attr_display(block, 'aln_span')
                # query region
                query_start = HSP._attr_display(block, 'query_start')
                query_end = HSP._attr_display(block, 'query_end')
                query_range = '[%s:%s]' % (query_start, query_end)
                # max column length is 20
                query_range = HSP._trunc_display(query_range, 22, '~]')
                # hit region
                hit_start = HSP._attr_display(block, 'hit_start')
                hit_end = HSP._attr_display(block, 'hit_end')
                hit_range = '[%s:%s]' % (hit_start, hit_end)
                hit_range = HSP._trunc_display(hit_range, 22, '~]')
                # append the hsp row
                lines.append(pattern % (str(idx), aln_span, query_range, hit_range))

            return self._str_hsp_header() + '\n' + '\n'.join(lines)

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

    def _validate_fragment(self, fragment):
        if not isinstance(fragment, HSPFragment):
            raise TypeError("HSP objects can only contain HSPFragment "
                    "objects.")

    def _aln_span_get(self):
        # length of all alignments
        # alignment span can be its own attribute, or computed from
        # query / hit length
        if len(self) == 1:
            self._aln_span = self._items[0].aln_span
        else:
            if not hasattr(self, '_aln_span'):
                self._aln_span = sum([frg.aln_span for frg in self.fragments])

        return self._aln_span

    def _aln_span_set(self, value):
        self._aln_span = value

    aln_span = property(fget=_aln_span_get, fset=_aln_span_set, \
            doc="""Total number of columns in all HSPFragment objects.""")

    ## coordinate properties ##
    def _get_coords(self, seq_type, coord_type):
        assert seq_type in ('hit', 'query')
        assert coord_type in ('start', 'end')
        coord_name = '%s_%s' % (seq_type, coord_type)
        coords = [getattr(frag, coord_name) for frag in self.fragments]
        if None in coords:
            warnings.warn("'None' exist in %s coordinates; ignored" %
                    (coord_name), BiopythonWarning)
        return coords

    def _hit_start_get(self):
        return min(self._get_coords('hit', 'start'))

    hit_start = property(fget=_hit_start_get, \
            doc="""Smallest coordinate value of all hit fragments""")

    def _query_start_get(self):
        return min(self._get_coords('query', 'start'))

    query_start = property(fget=_query_start_get, \
            doc="""Smallest coordinate value of all query fragments""")

    def _hit_end_get(self):
        return max(self._get_coords('hit', 'end'))

    hit_end = property(fget=_hit_end_get, \
            doc="""Largest coordinate value of all hit fragments""")

    def _query_end_get(self):
        return max(self._get_coords('query', 'end'))

    query_end = property(fget=_query_end_get, \
            doc="""Largest coordinate value of all hit fragments""")

    ## coordinate-dependent properties ##
    def _hit_span_get(self):
        try:
            return self.hit_end - self.hit_start
        except TypeError:  # triggered if any of the coordinates are None
            return None

    hit_span = property(fget=_hit_span_get, \
            doc="""The number of hit residues covered by the HSP.""")

    def _query_span_get(self):
        try:
            return self.query_end - self.query_start
        except TypeError:  # triggered if any of the coordinates are None
            return None

    query_span = property(fget=_query_span_get, \
            doc="""The number of query residues covered by the HSP.""")

    def _hit_range_get(self):
        return (self.hit_start, self.hit_end)

    hit_range = property(fget=_hit_range_get, \
            doc="""Tuple of HSP hit start and end coordinates.""")

    def _query_range_get(self):
        return (self.query_start, self.query_end)

    query_range = property(fget=_query_range_get, \
            doc="""Tuple of HSP query start and end coordinates.""")

    def _inter_ranges_get(self, seq_type):
        # this property assumes that there are no mixed strands in a hit/query
        assert seq_type in ('query', 'hit')
        strand = getattr(self, '%s_strands' % seq_type)[0]
        coords = getattr(self, '%s_ranges' % seq_type)
        # determine function used to set inter range
        # start and end coordinates, given two pairs
        # of fragment start and end coordinates
        if strand == -1:
            startfunc, endfunc = min, max
        else:
            startfunc, endfunc = max, min
        inter_coords = []
        for idx, coord in enumerate(coords[:-1]):
            start = startfunc(coords[idx])
            end = endfunc(coords[idx+1])
            inter_coords.append((min(start, end), max(start, end)))

        return inter_coords

    def _hit_inter_ranges_get(self):
        return self._inter_ranges_get('hit')

    hit_inter_ranges = property(fget=_hit_inter_ranges_get,
        doc="""Hit sequence coordinates of the regions between fragments""")

    def _query_inter_ranges_get(self):
        return self._inter_ranges_get('query')

    query_inter_ranges = property(fget=_query_inter_ranges_get,
        doc="""Query sequence coordinates of the regions between fragments""")

    ## shortcuts for fragments' properties ##

    # bool check if there's more than one fragments
    is_fragmented = property(lambda self: len(self) > 1,
            doc="""Whether the HSP has more than one HSPFragment objects""")

    # first item properties with setters
    hit_description = _fullcascade('hit_description',
            doc="""Description of the hit sequence""")

    query_description = _fullcascade('query_description',
            doc="""Description of the query sequence""")

    hit_id = _fullcascade('hit_id',
            doc="""ID of the hit sequence""")

    query_id = _fullcascade('query_id',
            doc="""ID of the query sequence""")

    alphabet = _fullcascade('alphabet',
            doc="""Alphabet used in hit and query SeqRecord objects""")

    # properties for single-fragment HSPs
    fragment = _singleitem(\
            doc="""HSPFragment object, first fragment""")

    hit = _singleitem('hit',
            doc="""Hit sequence as a SeqRecord object, first fragment""")

    query = _singleitem('query',
            doc="""Query sequence as a SeqRecord object, first fragment""")

    alignment = _singleitem('alignment',
            doc="""Alignment of the first fragment as a MultipleSeqAlignment object""")

    alignment_annotation = _singleitem('alignment_annotation',
            doc="""Dictionary of annotation(s) of the first fragment's alignment""")

    hit_strand = _singleitem('hit_strand',
            doc="""Hit strand orientation, first fragment""")

    query_strand = _singleitem('query_strand',
            doc="""Query strand orientation, first fragment""")

    hit_frame = _singleitem('hit_frame',
            doc="""Hit sequence reading frame, first fragment""")

    query_frame = _singleitem('query_frame',
            doc="""Query sequence reading frame, first fragment""")

    # properties for multi-fragment HSPs
    fragments = _allitems(doc="""List of all HSPFragment objects""")

    hits = _allitems('hit',
            doc="""List of all fragments' hit sequences as SeqRecord objects""")

    queries = _allitems('query',
            doc="""List of all fragments' query sequences as SeqRecord objects""")

    alignments = _allitems('alignment',
            doc="""List of all fragments' alignments as MultipleSeqAlignment objects""")

    alignment_annotations = _allitems('alignment_annotation',
            doc="""Dictionary of annotation(s) of all fragments' alignments""")

    hit_strands = _allitems('hit_strand',
            doc="""List of all fragments' hit sequence strands""")

    query_strands = _allitems('query_strand',
            doc="""List of all fragments' query sequence strands""")

    hit_frames = _allitems('hit_frame',
            doc="""List of all fragments' hit sequence reading frames""")

    query_frames = _allitems('query_frame',
            doc="""List of all fragments' query sequence reading frames""")

    hit_starts = _allitems('hit_start',
            doc="""List of all fragments' hit start coordinates""")

    query_starts = _allitems('query_starts',
            doc="""List of all fragments' query start coordinates""")

    hit_ends = _allitems('hit_ends',
            doc="""List of all fragments' hit end coordinates""")

    query_ends = _allitems('query_ends',
            doc="""List of all fragments' query end coordinates""")

    hit_spans = _allitems('hit_span',
            doc="""List of all fragments' hit sequence size""")

    query_spans = _allitems('query_span',
            doc="""List of all fragments' query sequence size""")

    hit_ranges = _allitems('hit_range',
            doc="""List of all fragments' hit start and end coordinates""")

    query_ranges = _allitems('query_range',
            doc="""List of all fragments' query start and end coordinates""")


class HSPFragment(_BaseHSP):

    """Class representing a contiguous alignment of hit-query sequence.

    HSPFragment forms the core of any parsed search output file. Depending on
    the search output file format, it may contain the actual query and/or hit
    sequences that produces the search hits. These sequences are stored as
    SeqRecord objects (see SeqRecord):

    >>> from Bio import SearchIO
    >>> qresult = SearchIO.parse('Blast/mirna.xml', 'blast-xml').next()
    >>> fragment = qresult[0][0][0]   # first hit, first hsp, first fragment
    >>> print fragment
          Query: 33211 mir_1
            Hit: gi|262205317|ref|NR_030195.1| Homo sapiens microRNA 520b (MIR520...
    Query range: [0:61] (1)
      Hit range: [0:61] (1)
      Fragments: 1 (61 columns)
         Query - CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG
                 |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
           Hit - CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG

    # the query sequence is a SeqRecord object
    >>> fragment.query.__class__
    <class 'Bio.SeqRecord.SeqRecord'>
    >>> print fragment.query
    ID: 33211
    Name: aligned query sequence
    Description: mir_1
    Number of features: 0
    Seq('CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTT...GGG', DNAAlphabet())

    # the hit sequence is a SeqRecord object as well
    >>> fragment.hit.__class__
    <class 'Bio.SeqRecord.SeqRecord'>
    >>> print fragment.hit
    ID: gi|262205317|ref|NR_030195.1|
    Name: aligned hit sequence
    Description: Homo sapiens microRNA 520b (MIR520B), microRNA
    Number of features: 0
    Seq('CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTT...GGG', DNAAlphabet())

    # when both query and hit are present, we get a MultipleSeqAlignment object
    >>> fragment.alignment.__class__
    <class 'Bio.Align.MultipleSeqAlignment'>
    >>> print fragment.alignment
    DNAAlphabet() alignment with 2 rows and 61 columns
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAG...GGG 33211
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAG...GGG gi|262205317|ref|NR_030195.1|

    """

    def __init__(self, hit_id='<unknown id>', query_id='<unknown id>',
            hit='', query='',
            hit_description='<unknown description>',
            query_description='<unknown description>',
            aln_annotation=None, alphabet=single_letter_alphabet):

        # no callables in default args!
        if aln_annotation is None:
            self.alignment_annotation = {}
        else:
            self.alignment_annotation = aln_annotation

        self._hit_id = hit_id
        self._query_id = query_id
        self._hit_description = hit_description
        self._query_description = query_description
        self._alphabet = alphabet

        for seq_type in ('query', 'hit'):
            # self.query or self.hit
            if eval(seq_type):
                setattr(self, seq_type, eval(seq_type))
            else:
                setattr(self, seq_type, None)
            # query or hit attributes
            for attr in ('strand', 'frame', 'start', 'end'):
                setattr(self, '%s_%s' % (seq_type, attr), None)

    def __repr__(self):
        info = "hit_id=%r, query_id=%r" % (self.hit_id, self.query_id)
        try:
            info += ", %i columns" % len(self)
        except AttributeError:
            pass
        return "%s(%s)" % (self.__class__.__name__, info)

    def __len__(self):
        return self.aln_span

    def __str__(self):
        return self._str_hsp_header() + '\n' + self._str_aln()

    def __getitem__(self, idx):
        if self.alignment is not None:
            obj = self.__class__(
                    hit_id=self.hit_id, query_id=self.query_id,
                    hit_description=self.hit_description,
                    query_description=self.query_description,
                    alphabet=self.alphabet)
            # transfer query and hit attributes
            if self.query is not None:
                obj.query = self.query[idx]
            if self.hit is not None:
                obj.hit = self.hit[idx]
            # and strand
            obj.hit_strand = self.hit_strand
            obj.query_strand = self.query_strand
            # alignment annotation should be transferred, since we can compute
            # the resulting annotation
            obj.alignment_annotation = {}
            for key, value in self.alignment_annotation.items():
                assert len(value[idx]) == len(obj)
                obj.alignment_annotation[key] = value[idx]
            return obj
        else:
            raise TypeError("Slicing for HSP objects without "
                    "alignment is not supported.")

    def _str_aln(self):
        lines = []
        # alignment length
        aln_span = HSPFragment._attr_display(self, 'aln_span')
        lines.append('  Fragments: 1 (%s columns)' % aln_span)
        # sequences
        if self.query is not None and self.hit is not None:
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

            if self.aln_span <= 67:
                lines.append("%10s - %s" % ('Query', qseq))
                if homol:
                    lines.append("             %s" % homol)
                lines.append("%10s - %s" % ('Hit', hseq))
            else:
                # adjust continuation character length, so we don't display
                # the same residues twice
                if self.aln_span - 66 > 3:
                    cont = '~' * 3
                else:
                    cont = '~' * (self.aln_span - 66)
                lines.append("%10s - %s%s%s" % ('Query',
                                qseq[:59], cont, qseq[-5:]))
                if homol:
                    lines.append("             %s%s%s" %
                            (homol[:59], cont, homol[-5:]))
                lines.append("%10s - %s%s%s" % ('Hit',
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
        opp_type = 'hit' if seq_type == 'query' else 'query'
        opp_seq = getattr(self, opp_type, None)
        if opp_seq is not None:
            if len(seq) != len(opp_seq):
                raise ValueError("Sequence lengths do not match. Expected: "
                        "%r (%s); found: %r (%s)." % (len(opp_seq), opp_type,
                        len(seq), seq_type))

        seq_id = getattr(self, '%s_id' % seq_type)
        seq_desc = getattr(self, '%s_description' % seq_type)
        seq_name = 'aligned %s sequence' % seq_type

        if isinstance(seq, SeqRecord):
            seq.id = seq_id
            seq.description = seq_desc
            seq.name = seq_name
            seq.seq.alphabet = self.alphabet
            return seq
        elif isinstance(seq, basestring):
            return SeqRecord(Seq(seq, self.alphabet), id=seq_id, name=seq_name,
                    description=seq_desc)
        else:
            raise TypeError("%s sequence must be a string or a "
                    "SeqRecord object." % seq_type.capitalize())

    def _hit_get(self):
        return self._hit

    def _hit_set(self, value):
        self._hit = self._prep_seq(value, 'hit')

    hit = property(fget=_hit_get, fset=_hit_set,
            doc="""Hit sequence as a SeqRecord object, defaults to None""")

    def _query_get(self):
        return self._query

    def _query_set(self, value):
        self._query = self._prep_seq(value, 'query')

    query = property(fget=_query_get, fset=_query_set,
            doc="""Query sequence as a SeqRecord object, defaults to None""")

    def _alignment_get(self):
        if self.query is None and self.hit is None:
            return None
        elif self.hit is None:
            return MultipleSeqAlignment([self.query], self.alphabet)
        elif self.query is None:
            return MultipleSeqAlignment([self.hit], self.alphabet)
        else:
            return MultipleSeqAlignment([self.query, self.hit], self.alphabet)

    alignment = property(fget=_alignment_get,
            doc="""Query-hit alignment as a MultipleSeqAlignment object,
            defaults to None""")

    def _alphabet_get(self):
        return self._alphabet

    def _alphabet_set(self, value):
        self._alphabet = value
        # try to set SeqRecords' alphabets if they exist
        if self.query is not None:
            self.query.seq.alphabet = value
        if self.hit is not None:
            self.hit.seq.alphabet = value

    alphabet = property(fget=_alphabet_get, fset=_alphabet_set,
            doc="""Alphabet object used in the fragment's sequences and alignment,
            defaults to single_letter_alphabet""")

    def _aln_span_get(self):
        # length of alignment (gaps included)
        # alignment span can be its own attribute, or computed from
        # query / hit length
        if not hasattr(self, '_aln_span'):
            if self.query is not None:
                self._aln_span = len(self.query)
            elif self.hit is not None:
                self._aln_span = len(self.hit)

        return self._aln_span

    def _aln_span_set(self, value):
        self._aln_span = value

    aln_span = property(fget=_aln_span_get, fset=_aln_span_set,
            doc="""The number of alignment columns covered by the fragment""")

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

    hit_description = property(fget=_hit_description_get,
            fset=_hit_description_set,
            doc="""Hit sequence description""")

    def _query_description_get(self):
        return self._query_description

    def _query_description_set(self, value):
        self._set_id_or_desc(value, 'query', 'description')

    query_description = property(fget=_query_description_get,
            fset=_query_description_set,
            doc="""Query sequence description""")

    def _hit_id_get(self):
        return self._hit_id

    def _hit_id_set(self, value):
        self._set_id_or_desc(value, 'hit', 'id')

    hit_id = property(fget=_hit_id_get, fset=_hit_id_set,
            doc="""Hit ID string""")

    def _query_id_get(self):
        return self._query_id

    def _query_id_set(self, value):
        self._set_id_or_desc(value, 'query', 'id')

    query_id = property(fget=_query_id_get, fset=_query_id_set,
            doc="""Query ID string""")

    ## strand properties ##
    def _prep_strand(self, strand):
        # follow SeqFeature's convention
        if not strand in (-1, 0, 1, None):
            raise ValueError("Strand should be -1, 0, 1, or None; not %r" %
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
                    strand = frame // abs(frame)
                except ZeroDivisionError:
                    strand = 0
                setattr(self, '%s_strand' % seq_type, strand)

        return strand

    def _hit_strand_get(self):
        return self._get_strand('hit')

    def _hit_strand_set(self, value):
        self._hit_strand = self._prep_strand(value)

    hit_strand = property(fget=_hit_strand_get, fset=_hit_strand_set,
            doc="""Hit sequence strand, defaults to None""")

    def _query_strand_get(self):
        return self._get_strand('query')

    def _query_strand_set(self, value):
        self._query_strand = self._prep_strand(value)

    query_strand = property(fget=_query_strand_get, fset=_query_strand_set,
            doc="""Query sequence strand, defaults to None""")

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

    hit_frame = property(fget=_hit_frame_get, fset=_hit_frame_set,
            doc="""Hit sequence reading frame, defaults to None""")

    def _query_frame_get(self):
        return self._query_frame

    def _query_frame_set(self, value):
        self._query_frame = self._prep_frame(value)

    query_frame = property(fget=_query_frame_get, fset=_query_frame_set,
            doc="""Query sequence reading frame, defaults to None""")

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

    hit_start = property(fget=_hit_start_get, fset=_hit_start_set,
            doc="""Hit sequence start coordinate, defaults to None""")

    def _query_start_get(self):
        return self._query_start

    def _query_start_set(self, value):
        self._query_start = self._prep_coord(value)

    query_start = property(fget=_query_start_get, fset=_query_start_set,
            doc="""Query sequence start coordinate, defaults to None""")

    def _hit_end_get(self):
        return self._hit_end

    def _hit_end_set(self, value):
        self._hit_end = self._prep_coord(value)

    hit_end = property(fget=_hit_end_get, fset=_hit_end_set,
            doc="""Hit sequence start coordinate, defaults to None""")

    def _query_end_get(self):
        return self._query_end

    def _query_end_set(self, value):
        self._query_end = self._prep_coord(value)

    query_end = property(fget=_query_end_get, fset=_query_end_set,
            doc="""Query sequence end coordinate, defaults to None""")

    ## coordinate-dependent properties ##
    def _hit_span_get(self):
        try:
            return self.hit_end - self.hit_start
        except TypeError:  # triggered if any of the coordinates are None
            return None

    hit_span = property(fget=_hit_span_get,
            doc="""The number of residues covered by the hit sequence""")

    def _query_span_get(self):
        try:
            return self.query_end - self.query_start
        except TypeError:  # triggered if any of the coordinates are None
            return None

    query_span = property(fget=_query_span_get,
            doc="""The number of residues covered by the query sequence""")

    def _hit_range_get(self):
        return (self.hit_start, self.hit_end)

    hit_range = property(fget=_hit_range_get,
            doc="""Tuple of hit start and end coordinates""")

    def _query_range_get(self):
        return (self.query_start, self.query_end)

    query_range = property(fget=_query_range_get,
            doc="""Tuple of query start and end coordinates""")


def _test():
    """Run the Bio.SearchIO._object module's doctests.

    This will try and locate the unit tests directory, and run the doctests
    from there in order that the relative paths used in the examples work.
    """
    import doctest
    import os

    test_dir = 'Tests'

    if os.path.isdir(os.path.join('..', '..', test_dir)):
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join('..', '..', test_dir))
        doctest.testmod(optionflags=doctest.ELLIPSIS)
        os.chdir(cur_dir)
        print "Done"


# if not used as a module, run the doctest
if __name__ == "__main__":
    _test()
