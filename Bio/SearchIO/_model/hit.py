# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO object to model a single database hit."""

from __future__ import print_function

from itertools import chain

from Bio._py3k import filter

from Bio._utils import getattr_str, trim_str
from Bio.SearchIO._utils import allitems, optionalcascade

from ._base import _BaseSearchObject
from .hsp import HSP


__docformat__ = "restructuredtext en"


class Hit(_BaseSearchObject):

    """Class representing a single database hit of a search result.

    Hit objects are the second-level container in the SearchIO module. They
    are the objects contained within a QueryResult (see QueryResult). They
    themselves are container for HSP objects and will contain at least one
    HSP.

    To have a quick look at a Hit and its contents, invoke ``print`` on it::

        >>> from Bio import SearchIO
        >>> qresult = next(SearchIO.parse('Blast/mirna.xml', 'blast-xml'))
        >>> hit = qresult[3]
        >>> print(hit)
        Query: 33211
               mir_1
          Hit: gi|301171322|ref|NR_035857.1| (86)
               Pan troglodytes microRNA mir-520c (MIR520C), microRNA
         HSPs: ----  --------  ---------  ------  ---------------  ---------------------
                  #   E-value  Bit score    Span      Query range              Hit range
               ----  --------  ---------  ------  ---------------  ---------------------
                  0   8.9e-20     100.47      60           [1:61]                [13:73]
                  1   3.3e-06      55.39      60           [0:60]                [13:73]

    You can invoke ``len`` on a Hit object to see how many HSP objects it contains::

        >>> len(hit)
        2

    Hit objects behave very similar to Python lists. You can retrieve the HSP
    object inside a Hit using the HSP's integer index. Hit objects can also be
    sliced, which will return a new Hit objects containing only the sliced HSPs::

        # HSP items inside the Hit can be retrieved using its integer index
        >>> hit[0]
        HSP(hit_id='gi|301171322|ref|NR_035857.1|', query_id='33211', 1 fragments)

        # slicing returns a new Hit
        >>> hit
        Hit(id='gi|301171322|ref|NR_035857.1|', query_id='33211', 2 hsps)
        >>> hit[:1]
        Hit(id='gi|301171322|ref|NR_035857.1|', query_id='33211', 1 hsps)
        >>> print(hit[1:])
        Query: 33211
               mir_1
          Hit: gi|301171322|ref|NR_035857.1| (86)
               Pan troglodytes microRNA mir-520c (MIR520C), microRNA
         HSPs: ----  --------  ---------  ------  ---------------  ---------------------
                  #   E-value  Bit score    Span      Query range              Hit range
               ----  --------  ---------  ------  ---------------  ---------------------
                  0   3.3e-06      55.39      60           [0:60]                [13:73]

    Hit objects provide ``filter`` and ``map`` methods, which are analogous to
    Python's built-in ``filter`` and ``map`` except that they return a new Hit
    object instead of a list.

    Here is an example of using ``filter`` to select for HSPs whose e-value is
    less than 1e-10::

        >>> evalue_filter = lambda hsp: hsp.evalue < 1e-10
        >>> filtered_hit = hit.filter(evalue_filter)
        >>> len(hit)
        2
        >>> len(filtered_hit)
        1
        >>> print(filtered_hit)
        Query: 33211
               mir_1
          Hit: gi|301171322|ref|NR_035857.1| (86)
               Pan troglodytes microRNA mir-520c (MIR520C), microRNA
         HSPs: ----  --------  ---------  ------  ---------------  ---------------------
                  #   E-value  Bit score    Span      Query range              Hit range
               ----  --------  ---------  ------  ---------------  ---------------------
                  0   8.9e-20     100.47      60           [1:61]                [13:73]

    There are also other methods which are counterparts of Python lists' methods
    with the same names: ``append``, ``index``, ``pop``, and ``sort``. Consult their
    respective documentations for more details and examples of their usage.

    """

    # attributes we don't want to transfer when creating a new Hit class
    # from this one
    _NON_STICKY_ATTRS = ('_items', )

    def __init__(self, hsps=[], id=None, query_id=None):
        """Initializes a Hit object.

        :param hsps: HSP objects contained in the Hit object
        :type hsps: iterable yielding HSP
        :param id: hit ID
        :type id: string
        :param query_id: query ID
        :type query_id: string

        If multiple HSP objects are used for initialization, they must all
        have the same ``query_id``, ``query_description``, ``hit_id``, and
        ``hit_description`` properties.
        """
        # default attribute values
        self._id = id
        self._id_alt = []
        self._query_id = query_id
        self._description = None
        self._description_alt = []
        self._query_description = None

        for attr in ('query_id', 'query_description', 'hit_id',
                'hit_description'):
            # HACK: setting the if clause to '> 1' allows for empty hit objects.
            # This makes it easier to work with file formats with unpredictable
            # hit-hsp ordering. The empty hit object itself is nonfunctional,
            # however, since all its cascading properties are empty.
            if len(set(getattr(hsp, attr) for hsp in hsps)) > 1:
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

    # Python 3:
    def __bool__(self):
        return bool(self.hsps)

    # Python 2:
    __nonzero__= __bool__

    def __contains__(self, hsp):
        return hsp in self._items

    def __str__(self):
        lines = []

        # set query id line
        qid_line = 'Query: %s' % self.query_id
        if self.query_description:
            qid_line += trim_str('\n       %s' %
                    self.query_description, 80, '...')
        lines.append(qid_line)

        # set hit id line
        hid_line = '  Hit: %s' % self.id
        if hasattr(self, 'seq_len'):
            hid_line += ' (%i)' % self.seq_len
        if self.description:
            hid_line += trim_str('\n       %s' % self.description,
                    80, '...')
        lines.append(hid_line)

        # set hsp line and table
        if not self.hsps:
            lines.append(' HSPs: ?')
        else:
            lines.append(' HSPs: %s  %s  %s  %s  %s  %s' %
                    ('-'*4, '-'*8, '-'*9, '-'*6, '-'*15, '-'*21))
            pattern = '%11s  %8s  %9s  %6s  %15s  %21s'
            lines.append(pattern % ('#', 'E-value', 'Bit score', 'Span',
                    'Query range', 'Hit range'))
            lines.append(pattern % ('-'*4, '-'*8, '-'*9, '-'*6, '-'*15, '-'*21))
            for idx, hsp in enumerate(self.hsps):
                # evalue
                evalue = getattr_str(hsp, 'evalue', fmt='%.2g')
                # bitscore
                bitscore = getattr_str(hsp, 'bitscore', fmt='%.2f')
                # alignment length
                aln_span = getattr_str(hsp, 'aln_span')
                # query region
                query_start = getattr_str(hsp, 'query_start')
                query_end = getattr_str(hsp, 'query_end')
                query_range = '[%s:%s]' % (query_start, query_end)
                # max column length is 18
                query_range = trim_str(query_range, 15, '~]')
                # hit region
                hit_start = getattr_str(hsp, 'hit_start')
                hit_end = getattr_str(hsp, 'hit_end')
                hit_range = '[%s:%s]' % (hit_start, hit_end)
                hit_range = trim_str(hit_range, 21, '~]')
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

    # hsp properties #
    def _validate_hsp(self, hsp):
        """Validates an HSP object.

        Valid HSP objects have the same hit_id as the Hit object ID and the
        same query_id as the Hit object's query_id.

        """
        if not isinstance(hsp, HSP):
            raise TypeError("Hit objects can only contain HSP objects.")
        # HACK: to make validation during __init__ work
        if self._items:
            if self.id is not None:
                if hsp.hit_id != self.id:
                    raise ValueError("Expected HSP with hit ID %r, "
                            "found %r instead." % (self.id, hsp.hit_id))
            else:
                self.id = hsp.hit_id

            if self.description is not None:
                if hsp.hit_description != self.description:
                    raise ValueError("Expected HSP with hit description %r, "
                            "found %r instead." % (self.description,
                        hsp.hit_description))
            else:
                self.description = hsp.hit_description

            if self.query_id is not None:
                if hsp.query_id != self.query_id:
                    raise ValueError("Expected HSP with query ID %r, "
                            "found %r instead." % (self.query_id, hsp.query_id))
            else:
                self.query_id = hsp.query_id

            if self.query_description is not None:
                if hsp.query_description != self.query_description:
                    raise ValueError("Expected HSP with query description %r, "
                            "found %r instead." % (self.query_description,
                            hsp.query_description))
            else:
                self.query_description = hsp.query_description

    # properties #
    description = optionalcascade('_description', 'hit_description',
            """Hit description""")
    query_description = optionalcascade('_query_description',
            'query_description',
            """Description of the query that produced the hit""")
    id = optionalcascade('_id', 'hit_id', """Hit ID string.""")
    query_id = optionalcascade('_query_id', 'query_id',
            """ID string of the query that produced the hit""")
    # returns all hsps
    hsps = allitems(doc="""HSP objects contained in the Hit""")

    @property
    def id_all(self):
        """Alternative ID(s) of the Hit"""
        return [self.id] + self._id_alt

    @property
    def description_all(self):
        """Alternative descriptions of the Hit"""
        return [self.description] + self._description_alt

    @property
    def fragments(self):
        """HSPFragment objects contained in the Hit"""
        return [frag for frag in chain(*self._items)]

    # public methods #
    def append(self, hsp):
        """Adds a HSP object to the end of Hit.

        Parameters
        hsp -- HSP object to append.

        Any HSP object appended must have the same ``hit_id`` property as the
        Hit object's ``id`` property and the same ``query_id`` property as the
        Hit object's ``query_id`` property.

        """
        self._validate_hsp(hsp)
        self._items.append(hsp)

    def filter(self, func=None):
        """Creates a new Hit object whose HSP objects pass the filter
        function.

        :param func: function for filtering
        :type func: callable, accepts HSP, returns bool

        ``filter`` is analogous to Python's built-in ``filter`` function, except
        that instead of returning a list it returns a ``Hit`` object. Here is an
        example of using ``filter`` to select for HSPs having bitscores bigger
        than 60::

            >>> from Bio import SearchIO
            >>> qresult = next(SearchIO.parse('Blast/mirna.xml', 'blast-xml'))
            >>> hit = qresult[3]
            >>> evalue_filter = lambda hsp: hsp.bitscore > 60
            >>> filtered_hit = hit.filter(evalue_filter)
            >>> len(hit)
            2
            >>> len(filtered_hit)
            1
            >>> print(filtered_hit)
            Query: 33211
                   mir_1
              Hit: gi|301171322|ref|NR_035857.1| (86)
                   Pan troglodytes microRNA mir-520c (MIR520C), microRNA
             HSPs: ----  --------  ---------  ------  ---------------  ---------------------
                      #   E-value  Bit score    Span      Query range              Hit range
                   ----  --------  ---------  ------  ---------------  ---------------------
                      0   8.9e-20     100.47      60           [1:61]                [13:73]

        """
        hsps = list(filter(func, self.hsps))
        if hsps:
            obj = self.__class__(hsps)
            self._transfer_attrs(obj)
            return obj

    def index(self, hsp):
        """Returns the index of a given HSP object, zero-based.

        :param hsp: object to look up
        :type hsp: HSP

        """
        return self._items.index(hsp)

    def map(self, func=None):
        """Creates a new Hit object, mapping the given function to its HSPs.

        :param func: function for mapping
        :type func: callable, accepts HSP, returns HSP

        ``map`` is analogous to Python's built-in ``map`` function. It is applied to
        all HSPs contained in the Hit object and returns a new Hit object.

        """
        if func is not None:
            hsps = [func(x) for x in self.hsps[:]] # this creates a shallow copy
        else:
            hsps = self.hsps[:]
        if hsps:
            obj = self.__class__(hsps)
            self._transfer_attrs(obj)
            return obj

    def pop(self, index=-1):
        """Removes and returns the HSP object at the specified index.

        :param index: index of HSP object to pop
        :type index: int

        """
        return self._items.pop(index)

    def sort(self, key=None, reverse=False, in_place=True):
        """Sorts the HSP objects.

        :param key: sorting function
        :type key: callable, accepts HSP, returns key for sorting
        :param reverse: whether to reverse sorting results or no
        :type reverse: bool
        :param in_place: whether to do in-place sorting or no
        :type in_place: bool

        ``sort`` defaults to sorting in-place, to mimick Python's ``list.sort``
        method. If you set the ``in_place`` argument to False, it will treat
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


# if not used as a module, run the doctest
if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()
