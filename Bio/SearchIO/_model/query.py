# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.SearchIO object to model search results from a single query."""


from copy import deepcopy
from itertools import chain

from Bio.SearchIO._utils import optionalcascade

from ._base import _BaseSearchObject
from .hit import Hit


class QueryResult(_BaseSearchObject):
    """Class representing search results from a single query.

    QueryResult is the container object that stores all search hits from a
    single search query. It is the top-level object returned by SearchIO's two
    main functions, ``read`` and ``parse``. Depending on the search results and
    search output format, a QueryResult object will contain zero or more Hit
    objects (see Hit).

    You can take a quick look at a QueryResult's contents and attributes by
    invoking ``print`` on it::

        >>> from Bio import SearchIO
        >>> qresult = next(SearchIO.parse('Blast/mirna.xml', 'blast-xml'))
        >>> print(qresult)
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
        ...

    If you just want to know how many hits a QueryResult has, you can invoke
    ``len`` on it. Alternatively, you can simply type its name in the interpreter::

        >>> len(qresult)
        100
        >>> qresult
        QueryResult(id='33211', 100 hits)

    QueryResult behaves like a hybrid of Python's built-in list and dictionary.
    You can retrieve its items (Hit objects) using the integer index of the
    item, just like regular Python lists::

        >>> first_hit = qresult[0]
        >>> first_hit
        Hit(id='gi|262205317|ref|NR_030195.1|', query_id='33211', 1 hsps)

    You can slice QueryResult objects as well. Slicing will return a new
    QueryResult object containing only the sliced hits::

        >>> sliced_qresult = qresult[:3]    # slice the first three hits
        >>> len(qresult)
        100
        >>> len(sliced_qresult)
        3
        >>> print(sliced_qresult)
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

    Like Python dictionaries, you can also retrieve hits using the hit's ID.
    This is useful for retrieving hits that you know should exist in a given
    search::

        >>> hit = qresult['gi|262205317|ref|NR_030195.1|']
        >>> hit
        Hit(id='gi|262205317|ref|NR_030195.1|', query_id='33211', 1 hsps)

    You can also replace a Hit in QueryResult with another Hit using either the
    integer index or hit key string. Note that the replacing object must be a
    Hit that has the same ``query_id`` property as the QueryResult object.

    If you're not sure whether a QueryResult contains a particular hit, you can
    use the hit ID to check for membership first::

        >>> 'gi|262205317|ref|NR_030195.1|' in qresult
        True
        >>> 'gi|262380031|ref|NR_023426.1|' in qresult
        False

    Or, if you just want to know the rank / position of a given hit, you can
    use the hit ID as an argument for the ``index`` method. Note that the values
    returned will be zero-based. So zero (0) means the hit is the first in the
    QueryResult, three (3) means the hit is the fourth item, and so on. If the
    hit does not exist in the QueryResult, a ``ValueError`` will be raised.

        >>> qresult.index('gi|262205317|ref|NR_030195.1|')
        0
        >>> qresult.index('gi|262205330|ref|NR_030198.1|')
        5
        >>> qresult.index('gi|262380031|ref|NR_023426.1|')
        Traceback (most recent call last):
        ...
        ValueError: ...

    To ease working with a large number of hits, QueryResult has several
    ``filter`` and ``map`` methods, analogous to Python's built-in functions with
    the same names. There are ``filter`` and ``map`` methods available for
    operations over both Hit objects or HSP objects. As an example, here we are
    using the ``hit_map`` method to rename all hit IDs within a QueryResult::

        >>> def renamer(hit):
        ...     hit.id = hit.id.split('|')[3]
        ...     return hit
        >>> mapped_qresult = qresult.hit_map(renamer)
        >>> print(mapped_qresult)
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
        ...

    The principle for other ``map`` and ``filter`` methods are similar: they accept
    a function, applies it, and returns a new QueryResult object.

    There are also other methods useful for working with list-like objects:
    ``append``, ``pop``, and ``sort``. More details and examples are available in
    their respective documentations.

    Finally, just like Python lists and dictionaries, QueryResult objects are
    iterable. Iteration over QueryResults will yield Hit objects::

        >>> for hit in qresult[:4]:     # iterate over the first four items
        ...     hit
        ...
        Hit(id='gi|262205317|ref|NR_030195.1|', query_id='33211', 1 hsps)
        Hit(id='gi|301171311|ref|NR_035856.1|', query_id='33211', 1 hsps)
        Hit(id='gi|270133242|ref|NR_032573.1|', query_id='33211', 1 hsps)
        Hit(id='gi|301171322|ref|NR_035857.1|', query_id='33211', 2 hsps)

    If you need access to all the hits in a QueryResult object, you can get
    them in a list using the ``hits`` property. Similarly, access to all hit IDs is
    available through the ``hit_keys`` property.

        >>> qresult.hits
        [Hit(id='gi|262205317|ref|NR_030195.1|', query_id='33211', 1 hsps), ...]
        >>> qresult.hit_keys
        ['gi|262205317|ref|NR_030195.1|', 'gi|301171311|ref|NR_035856.1|', ...]

    """

    # attributes we don't want to transfer when creating a new QueryResult class
    # from this one
    _NON_STICKY_ATTRS = ("_items", "__alt_hit_ids")

    def __init__(self, hits=(), id=None, hit_key_function=None):
        """Initialize a QueryResult object.

        :param id: query sequence ID
        :type id: string
        :param hits: iterator yielding Hit objects
        :type hits: iterable
        :param hit_key_function: function to define hit keys
        :type hit_key_function: callable, accepts Hit objects, returns string

        """
        # default values
        self._id = id
        self._hit_key_function = hit_key_function or _hit_key_func
        self._items = {}
        self._description = None
        self.__alt_hit_ids = {}
        self.program = "<unknown program>"
        self.target = "<unknown target>"
        self.version = "<unknown version>"

        # validate Hit objects and fill up self._items
        for hit in hits:
            # validation is handled by __setitem__
            self.append(hit)

    def __iter__(self):
        """Iterate over hits."""
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
        """Return an iterator over the Hit objects."""
        yield from self._items.values()

    def iterhit_keys(self):
        """Return an iterator over the ID of the Hit objects."""
        yield from self._items

    def iteritems(self):
        """Return an iterator yielding tuples of Hit ID and Hit objects."""
        yield from self._items.items()

    def __contains__(self, hit_key):
        """Return True if hit key in items or alternative hit identifiers."""
        if isinstance(hit_key, Hit):
            return self._hit_key_function(hit_key) in self._items
        return hit_key in self._items or hit_key in self.__alt_hit_ids

    def __len__(self):
        """Return the number of items."""
        return len(self._items)

    def __bool__(self):
        """Return True if there are items."""
        return bool(self._items)

    def __repr__(self):
        """Return string representation of the QueryResult object."""
        return "QueryResult(id=%r, %r hits)" % (self.id, len(self))

    def __str__(self):
        """Return a human readable summary of the QueryResult object."""
        lines = []

        # set program and version line
        lines.append("Program: %s (%s)" % (self.program, self.version))

        # set query id line
        qid_line = "  Query: %s" % self.id
        try:
            seq_len = self.seq_len
        except AttributeError:
            pass
        else:
            qid_line += " (%i)" % seq_len
        lines.append(qid_line)
        if self.description:
            line = "         %s" % self.description
            line = line[:77] + "..." if len(line) > 80 else line
            lines.append(line)

        # set target line
        lines.append(" Target: %s" % self.target)

        # set hit lines
        if not self.hits:
            lines.append("   Hits: 0")
        else:
            lines.append("   Hits: %s  %s  %s" % ("-" * 4, "-" * 5, "-" * 58))
            pattern = "%13s  %5s  %s"
            lines.append(pattern % ("#", "# HSP", "ID + description"))
            lines.append(pattern % ("-" * 4, "-" * 5, "-" * 58))
            for idx, hit in enumerate(self.hits):
                if idx < 30:
                    hid_line = "%s  %s" % (hit.id, hit.description)
                    if len(hid_line) > 58:
                        hid_line = hid_line[:55] + "..."
                    lines.append(pattern % (idx, len(hit), hid_line))
                elif idx > len(self.hits) - 4:
                    hid_line = "%s  %s" % (hit.id, hit.description)
                    if len(hid_line) > 58:
                        hid_line = hid_line[:55] + "..."
                    lines.append(pattern % (idx, len(hit), hid_line))
                elif idx == 30:
                    lines.append("%14s" % "~~~")

        return "\n".join(lines)

    def __getitem__(self, hit_key):
        """Return a QueryResult object that matches the hit_key."""
        # retrieval using slice objects returns another QueryResult object
        if isinstance(hit_key, slice):
            # should we return just a list of Hits instead of a full blown
            # QueryResult object if it's a slice?
            hits = list(self.hits)[hit_key]
            obj = self.__class__(hits, self.id, self._hit_key_function)
            self._transfer_attrs(obj)
            return obj

        # if key is an int, then retrieve the Hit at the int index
        elif isinstance(hit_key, int):
            length = len(self)
            if 0 <= hit_key < length:
                for idx, item in enumerate(self.iterhits()):
                    if idx == hit_key:
                        return item
            elif -1 * length <= hit_key < 0:
                for idx, item in enumerate(self.iterhits()):
                    if length + hit_key == idx:
                        return item
            raise IndexError("list index out of range")

        # if key is a string, then do a regular dictionary retrieval
        # falling back on alternative hit IDs
        try:
            return self._items[hit_key]
        except KeyError:
            return self._items[self.__alt_hit_ids[hit_key]]

    def __setitem__(self, hit_key, hit):
        """Add an item of key hit_key and value hit."""
        # only accept string keys
        if not isinstance(hit_key, str):
            raise TypeError("QueryResult object keys must be a string.")
        # hit must be a Hit object
        if not isinstance(hit, Hit):
            raise TypeError("QueryResult objects can only contain Hit objects.")
        qid = self.id
        hqid = hit.query_id
        # and it must have the same query ID as this object's ID
        # unless it's the query ID is None (default for empty objects), in which
        # case we want to use the hit's query ID as the query ID
        if qid is not None:
            if hqid != qid:
                raise ValueError(
                    "Expected Hit with query ID %r, found %r instead." % (qid, hqid)
                )
        else:
            self.id = hqid
        # same thing with descriptions
        qdesc = self.description
        hqdesc = hit.query_description
        if qdesc is not None:
            if hqdesc != qdesc:
                raise ValueError(
                    "Expected Hit with query description %r, found %r instead."
                    % (qdesc, hqdesc)
                )
        else:
            self.description = hqdesc

        # remove existing alt_id references, if hit_key already exists
        if hit_key in self._items:
            for alt_key in self._items[hit_key].id_all[1:]:
                del self.__alt_hit_ids[alt_key]

        # if hit_key is already present as an alternative ID
        # delete it from the alternative ID dict
        if hit_key in self.__alt_hit_ids:
            del self.__alt_hit_ids[hit_key]

        self._items[hit_key] = hit
        for alt_id in hit.id_all[1:]:
            self.__alt_hit_ids[alt_id] = hit_key

    def __delitem__(self, hit_key):
        """Delete item of key hit_key."""
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
            deleted = False
            if key in self._items:
                del self._items[key]
                deleted = True
            if key in self.__alt_hit_ids:
                del self._items[self.__alt_hit_ids[key]]
                del self.__alt_hit_ids[key]
                deleted = True
            if not deleted:
                raise KeyError(repr(key))

    # properties #
    id = optionalcascade("_id", "query_id", """QueryResult ID string""")
    description = optionalcascade(
        "_description", "query_description", """QueryResult description"""
    )

    @property
    def hsps(self):
        """Access the HSP objects contained in the QueryResult."""
        return sorted(
            (hsp for hsp in chain(*self.hits)), key=lambda hsp: hsp.output_index
        )

    @property
    def fragments(self):
        """Access the HSPFragment objects contained in the QueryResult."""
        return list(chain(*self.hsps))

    # public methods #
    def absorb(self, hit):
        """Add a Hit object to the end of QueryResult.

        If the QueryResult already has a Hit with the same ID, append the new
        Hit's HSPs into the existing Hit.

        :param hit: object to absorb
        :type hit: Hit

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
        """Add a Hit object to the end of QueryResult.

        :param hit: object to append
        :type hit: Hit

        Any Hit object appended must have the same ``query_id`` property as the
        QueryResult's ``id`` property. If the hit key already exists, a
        ``ValueError`` will be raised.

        """
        # if a custom hit_key_function is supplied, use it to define th hit key
        if self._hit_key_function is not None:
            hit_key = self._hit_key_function(hit)
        else:
            hit_key = hit.id

        if hit_key not in self and all(pid not in self for pid in hit.id_all[1:]):
            self[hit_key] = hit
        else:
            raise ValueError(
                "The ID or alternative IDs of Hit %r exists in this QueryResult."
                % hit_key
            )

    def hit_filter(self, func=None):
        """Create new QueryResult object whose Hit objects pass the filter function.

        :param func: filter function
        :type func: callable, accepts Hit, returns bool

        Here is an example of using ``hit_filter`` to select Hits whose
        description begins with the string 'Homo sapiens', case sensitive::

            >>> from Bio import SearchIO
            >>> qresult = next(SearchIO.parse('Blast/mirna.xml', 'blast-xml'))
            >>> def desc_filter(hit):
            ...     return hit.description.startswith('Homo sapiens')
            ...
            >>> len(qresult)
            100
            >>> filtered = qresult.hit_filter(desc_filter)
            >>> len(filtered)
            39
            >>> print(filtered[:4])
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
        hits = list(filter(func, self.hits))
        obj = self.__class__(hits, self.id, self._hit_key_function)
        self._transfer_attrs(obj)
        return obj

    def hit_map(self, func=None):
        """Create new QueryResult object, mapping the given function to its Hits.

        :param func: map function
        :type func: callable, accepts Hit, returns Hit

        Here is an example of using ``hit_map`` with a function that discards all
        HSPs in a Hit except for the first one::

            >>> from Bio import SearchIO
            >>> qresult = next(SearchIO.parse('Blast/mirna.xml', 'blast-xml'))
            >>> print(qresult[:8])
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
            >>> print(mapped_qresult[:8])
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
            hits = [func(x) for x in hits]
        obj = self.__class__(hits, self.id, self._hit_key_function)
        self._transfer_attrs(obj)
        return obj

    def hsp_filter(self, func=None):
        """Create new QueryResult object whose HSP objects pass the filter function.

        ``hsp_filter`` is the same as ``hit_filter``, except that it filters
        directly on each HSP object in every Hit. If the filtering removes
        all HSP objects in a given Hit, the entire Hit will be discarded. This
        will result in the QueryResult having less Hit after filtering.
        """
        hits = [x for x in (hit.filter(func) for hit in self.hits) if x]
        obj = self.__class__(hits, self.id, self._hit_key_function)
        self._transfer_attrs(obj)
        return obj

    def hsp_map(self, func=None):
        """Create new QueryResult object, mapping the given function to its HSPs.

        ``hsp_map`` is the same as ``hit_map``, except that it applies the given
        function to all HSP objects in every Hit, instead of the Hit objects.
        """
        hits = [x for x in (hit.map(func) for hit in list(self.hits)[:]) if x]
        obj = self.__class__(hits, self.id, self._hit_key_function)
        self._transfer_attrs(obj)
        return obj

    # marker for default self.pop() return value
    # this method is adapted from Python's built in OrderedDict.pop
    # implementation
    __marker = object()

    def pop(self, hit_key=-1, default=__marker):
        """Remove the specified hit key and return the Hit object.

        :param hit_key: key of the Hit object to return
        :type hit_key: int or string
        :param default: return value if no Hit exists with the given key
        :type default: object

        By default, ``pop`` will remove and return the last Hit object in the
        QueryResult object. To remove specific Hit objects, you can use its
        integer index or hit key.

            >>> from Bio import SearchIO
            >>> qresult = next(SearchIO.parse('Blast/mirna.xml', 'blast-xml'))
            >>> len(qresult)
            100
            >>> for hit in qresult[:5]:
            ...     print(hit.id)
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
            hit = self._items.pop(hit_key)
            # remove all alternative IDs of the popped hit
            for alt_id in hit.id_all[1:]:
                self.__alt_hit_ids.pop(alt_id, None)
        except KeyError:
            try:
                hit = self.pop(self.__alt_hit_ids[hit_key])
            except KeyError:
                # hit_key is not a valid id
                # use the default if it has been set
                if default is not self.__marker:
                    hit = default
                else:
                    raise KeyError(hit_key) from None
        return hit

    def index(self, hit_key):
        """Return the index of a given hit key, zero-based.

        :param hit_key: hit ID
        :type hit_key: string

        This method is useful for finding out the integer index (usually
        correlated with search rank) of a given hit key.

            >>> from Bio import SearchIO
            >>> qresult = next(SearchIO.parse('Blast/mirna.xml', 'blast-xml'))
            >>> qresult.index('gi|301171259|ref|NR_035850.1|')
            7

        """
        if isinstance(hit_key, Hit):
            return list(self.hit_keys).index(hit_key.id)
        try:
            return list(self.hit_keys).index(hit_key)
        except ValueError:
            if hit_key in self.__alt_hit_ids:
                return self.index(self.__alt_hit_ids[hit_key])
            raise

    def sort(self, key=None, reverse=False, in_place=True):
        """Sort the Hit objects.

        :param key: sorting function
        :type key: callable, accepts Hit, returns key for sorting
        :param reverse: whether to reverse sorting results or no
        :type reverse: bool
        :param in_place: whether to do in-place sorting or no
        :type in_place: bool

        ``sort`` defaults to sorting in-place, to mimic Python's ``list.sort``
        method. If you set the ``in_place`` argument to False, it will treat
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
            self._items = {self._hit_key_function(hit): hit for hit in sorted_hits}
        # otherwise, return a new sorted QueryResult object
        else:
            obj = self.__class__(sorted_hits, self.id, self._hit_key_function)
            self._transfer_attrs(obj)
            return obj


def _hit_key_func(hit):
    """Map hit to its identifier (PRIVATE).

    Default hit key function for QueryResult.__init__ use.
    """
    return hit.id


# if not used as a module, run the doctest
if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
