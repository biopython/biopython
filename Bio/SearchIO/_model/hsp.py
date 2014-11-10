# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO objects to model high scoring regions between query and hit."""

from __future__ import print_function
from Bio._py3k import basestring

import warnings
from operator import ge, le

from Bio import BiopythonWarning
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import single_letter_alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from Bio._utils import getattr_str, trim_str
from Bio.SearchIO._utils import singleitem, allitems, fullcascade, \
        fragcascade

from ._base import _BaseHSP


__docformat__ = "restructuredtext en"


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
    are called 'blocks' while in Exonerate they are called exons or NER.

    Here are examples from each type of HSP. The first one comes from a BLAST
    search::

        >>> from Bio import SearchIO
        >>> blast_qresult = next(SearchIO.parse('Blast/mirna.xml', 'blast-xml'))
        >>> blast_hsp = blast_qresult[1][0]     # the first HSP from the second hit
        >>> blast_hsp
        HSP(hit_id='gi|301171311|ref|NR_035856.1|', query_id='33211', 1 fragments)
        >>> print(blast_hsp)
              Query: 33211 mir_1
                Hit: gi|301171311|ref|NR_035856.1| Pan troglodytes microRNA mir-520b ...
        Query range: [1:61] (1)
          Hit range: [0:60] (1)
        Quick stats: evalue 1.7e-22; bitscore 109.49
          Fragments: 1 (60 columns)
             Query - CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG
                     ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
               Hit - CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG

    For HSPs with a single HSPFragment, you can invoke ``print`` on it and see the
    underlying sequence alignment, if it exists. This is not the case for HSPs
    with more than one HSPFragment. Below is an example, using an HSP from a
    BLAT search. Invoking ``print`` on these HSPs will instead show a table of the
    HSPFragment objects it contains::

        >>> blat_qresult = SearchIO.read('Blat/mirna.pslx', 'blat-psl', pslx=True)
        >>> blat_hsp = blat_qresult[1][0]       # the first HSP from the second hit
        >>> blat_hsp
        HSP(hit_id='chr11', query_id='blat_1', 2 fragments)
        >>> print(blat_hsp)
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

    Notice that in HSPs with more than one HSPFragments, the HSP's ``query_range``
    ``hit_range`` properties encompasses all fragments it contains.

    You can check whether an HSP has more than one HSPFragments or not using the
    ``is_fragmented`` property::

        >>> blast_hsp.is_fragmented
        False
        >>> blat_hsp.is_fragmented
        True

    Since HSP objects are also containers similar to Python lists, you can
    access a single fragment in an HSP using its integer index::

        >>> blat_fragment = blat_hsp[0]
        >>> print(blat_fragment)
              Query: blat_1 <unknown description>
                Hit: chr11 <unknown description>
        Query range: [61:67] (-1)
          Hit range: [59018929:59018935] (1)
          Fragments: 1 (6 columns)
             Query - tatagt
               Hit - tatagt

    This applies to HSPs objects with a single fragment as well::

        >>> blast_fragment = blast_hsp[0]
        >>> print(blast_fragment)
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

    +----------------------+---------------------+-----------------------------+
    | Property             | Shortcut            | Value                       |
    +======================+=====================+=============================+
    | aln_all              | aln                 | HSP alignments as           |
    |                      |                     | MultipleSeqAlignment object |
    +----------------------+---------------------+-----------------------------+
    | aln_annotation_all   | aln_annotation      | dictionary of annotation(s) |
    |                      |                     | of all fragments' alignments|
    +----------------------+---------------------+-----------------------------+
    | fragments            | fragment            | HSPFragment objects         |
    +----------------------+---------------------+-----------------------------+
    | hit_all              | hit                 | hit sequence as SeqRecord   |
    |                      |                     | objects                     |
    +----------------------+---------------------+-----------------------------+
    | hit_features_all     | hit_features        | SeqFeatures of all hit      |
    |                      |                     | fragments                   |
    +----------------------+---------------------+-----------------------------+
    | hit_start_all        | hit_start*          | start coordinates of the    |
    |                      |                     | hit fragments               |
    +----------------------+---------------------+-----------------------------+
    | hit_end_all          | hit_end*            | end coordinates of the hit  |
    |                      |                     | fragments                   |
    +----------------------+---------------------+-----------------------------+
    | hit_span_all         | hit_span*           | sizes of each hit fragments |
    +----------------------+---------------------+-----------------------------+
    | hit_strand_all       | hit_strand          | strand orientations of the  |
    |                      |                     | hit fragments               |
    +----------------------+---------------------+-----------------------------+
    | hit_frame_all        | hit_frame           | reading frames of the hit   |
    |                      |                     | fragments                   |
    +----------------------+---------------------+-----------------------------+
    | hit_range_all        | hit_range           | tuples of start and end     |
    |                      |                     | coordinates of each hit     |
    |                      |                     | fragment                    |
    +----------------------+---------------------+-----------------------------+
    | query_all            | query               | query sequence as SeqRecord |
    |                      |                     | object                      |
    +----------------------+---------------------+-----------------------------+
    | query_features_all   | query_features      | SeqFeatures of all query    |
    |                      |                     | fragments                   |
    +----------------------+---------------------+-----------------------------+
    | query_start_all      | query_start*        | start coordinates of the    |
    |                      |                     | fragments                   |
    +----------------------+---------------------+-----------------------------+
    | query_end_all        | query_end*          | end coordinates of the      |
    |                      |                     | query fragments             |
    +----------------------+---------------------+-----------------------------+
    | query_span_all       | query_span*         | sizes of each query         |
    |                      |                     | fragments                   |
    +----------------------+---------------------+-----------------------------+
    | query_strand_all     | query_strand        | strand orientations of the  |
    |                      |                     | query fragments             |
    +----------------------+---------------------+-----------------------------+
    | query_frame_all      | query_frame         | reading frames of the query |
    |                      |                     | fragments                   |
    +----------------------+---------------------+-----------------------------+
    | query_range_all      | query_range         | tuples of start and end     |
    |                      |                     | coordinates of each query   |
    |                      |                     | fragment                    |
    +----------------------+---------------------+-----------------------------+

    For all types of HSP objects, the property will return the values in a list.
    Shorcuts are only applicable for HSPs with one fragment. Except the ones
    noted, if they are used on an HSP with more than one fragments, an exception
    will be raised.

    For properties that may be used in HSPs with multiple or single fragments
    (``*_start``, ``*_end``, and ``*_span`` properties), their interpretation depends
    on how many fragment the HSP has:

    +------------+---------------------------------------------------+
    | Property   | Value                                             |
    +============+===================================================+
    | hit_start  | smallest coordinate value of all hit fragments    |
    +------------+---------------------------------------------------+
    | hit_end    | largest coordinate value of all hit fragments     |
    +------------+---------------------------------------------------+
    | hit_span   | difference between ``hit_start`` and ``hit_end``  |
    +------------+---------------------------------------------------+
    | query_start| smallest coordinate value of all query fragments  |
    +------------+---------------------------------------------------+
    | query_end  | largest coordinate value of all query fragments   |
    +------------+---------------------------------------------------+
    | query_span | difference between ``query_start`` and            |
    |            | ``query_end``                                     |
    +------------+---------------------------------------------------+

    In addition to the objects listed above, HSP objects also provide the
    following properties:

    +--------------------+------------------------------------------------------+
    | Property           | Value                                                |
    +====================+======================================================+
    | aln_span           | total number of residues in all HSPFragment objects  |
    +--------------------+------------------------------------------------------+
    | alphabet           | alphabet used in hit and query SeqRecord objects     |
    +--------------------+------------------------------------------------------+
    | is_fragmented      | boolean, whether there are multiple fragments or not |
    +--------------------+------------------------------------------------------+
    | hit_id             | ID of the hit sequence                               |
    +--------------------+------------------------------------------------------+
    | hit_description    | description of the hit sequence                      |
    +--------------------+------------------------------------------------------+
    | hit_inter_ranges   | list of hit sequence coordinates of the regions      |
    |                    | between fragments                                    |
    +--------------------+------------------------------------------------------+
    | hit_inter_spans    | list of lengths of the regions between hit fragments |
    +--------------------+------------------------------------------------------+
    | query_id           | ID of the query sequence                             |
    +--------------------+------------------------------------------------------+
    | query_description  | description of the query sequence                    |
    +--------------------+------------------------------------------------------+
    | query_inter_ranges | list of query sequence coordinates of the regions    |
    |                    | between fragments                                    |
    +--------------------+------------------------------------------------------+
    | query_inter_spans  | list of lengths of the regions between query         |
    |                    | fragments                                            |
    +--------------------+------------------------------------------------------+

    .. [1] may be used in HSPs with multiple fragments

    """
    # attributes we don't want to transfer when creating a new Hit class
    # from this one
    _NON_STICKY_ATTRS = ('_items', )

    def __init__(self, fragments=[]):
        """Initializes an HSP object.

        :param fragments: fragments contained in the HSP object
        :type fragments: iterable yielding HSPFragment

        HSP objects must be initialized with a list containing at least one
        HSPFragment object. If multiple HSPFragment objects are used for
        initialization, they must all have the same ``query_id``,
        ``query_description``, ``hit_id``, ``hit_description``, and alphabet
        properties.

        """
        if not fragments:
            raise ValueError("HSP objects must have at least one HSPFragment "
                    "object.")
        # check that all fragments contain the same IDs, descriptions, alphabet
        for attr in ('query_id', 'query_description', 'hit_id',
                'hit_description', 'alphabet'):
            if len(set(getattr(frag, attr) for frag in fragments)) != 1:
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

    # Python 3:
    def __bool__(self):
        return bool(self._items)

    # Python 2:
    __nonzero__= __bool__

    def __str__(self):

        lines = []
        # set hsp info line
        statline = []
        # evalue
        evalue = getattr_str(self, 'evalue', fmt='%.2g')
        statline.append('evalue ' + evalue)
        # bitscore
        bitscore = getattr_str(self, 'bitscore', fmt='%.2f')
        statline.append('bitscore ' + bitscore)
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
                aln_span = getattr_str(block, 'aln_span')
                # query region
                query_start = getattr_str(block, 'query_start')
                query_end = getattr_str(block, 'query_end')
                query_range = '[%s:%s]' % (query_start, query_end)
                # max column length is 20
                query_range = trim_str(query_range, 22, '~]')
                # hit region
                hit_start = getattr_str(block, 'hit_start')
                hit_end = getattr_str(block, 'hit_end')
                hit_range = '[%s:%s]' % (hit_start, hit_end)
                hit_range = trim_str(hit_range, 22, '~]')
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
        # HACK: to make validation during __init__ work
        if self._items:
            if fragment.hit_id != self.hit_id:
                raise ValueError("Expected HSPFragment with hit ID %r, "
                        "found %r instead." % (self.id, fragment.hit_id))

            if fragment.hit_description != self.hit_description:
                raise ValueError("Expected HSPFragment with hit "
                        "description %r, found %r instead." %
                        (self.description, fragment.hit_description))

            if fragment.query_id != self.query_id:
                raise ValueError("Expected HSPFragment with query ID %r, "
                        "found %r instead." % (self.query_id,
                        fragment.query_id))

            if fragment.query_description != self.query_description:
                raise ValueError("Expected HSP with query description %r, "
                        "found %r instead." % (self.query_description,
                        fragment.query_description))

    def _aln_span_get(self):
        # length of all alignments
        # alignment span can be its own attribute, or computed from
        # query / hit length
        return sum(frg.aln_span for frg in self.fragments)

    aln_span = property(fget=_aln_span_get,
            doc="""Total number of columns in all HSPFragment objects.""")

    # coordinate properties #
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

    hit_start = property(fget=_hit_start_get,
            doc="""Smallest coordinate value of all hit fragments""")

    def _query_start_get(self):
        return min(self._get_coords('query', 'start'))

    query_start = property(fget=_query_start_get,
            doc="""Smallest coordinate value of all query fragments""")

    def _hit_end_get(self):
        return max(self._get_coords('hit', 'end'))

    hit_end = property(fget=_hit_end_get,
            doc="""Largest coordinate value of all hit fragments""")

    def _query_end_get(self):
        return max(self._get_coords('query', 'end'))

    query_end = property(fget=_query_end_get,
            doc="""Largest coordinate value of all hit fragments""")

    # coordinate-dependent properties #
    def _hit_span_get(self):
        try:
            return self.hit_end - self.hit_start
        except TypeError:  # triggered if any of the coordinates are None
            return None

    hit_span = property(fget=_hit_span_get,
            doc="""The number of hit residues covered by the HSP.""")

    def _query_span_get(self):
        try:
            return self.query_end - self.query_start
        except TypeError:  # triggered if any of the coordinates are None
            return None

    query_span = property(fget=_query_span_get,
            doc="""The number of query residues covered by the HSP.""")

    def _hit_range_get(self):
        return (self.hit_start, self.hit_end)

    hit_range = property(fget=_hit_range_get,
            doc="""Tuple of HSP hit start and end coordinates.""")

    def _query_range_get(self):
        return (self.query_start, self.query_end)

    query_range = property(fget=_query_range_get,
            doc="""Tuple of HSP query start and end coordinates.""")

    def _inter_ranges_get(self, seq_type):
        # this property assumes that there are no mixed strands in a hit/query
        assert seq_type in ('query', 'hit')
        strand = getattr(self, '%s_strand_all' % seq_type)[0]
        coords = getattr(self, '%s_range_all' % seq_type)
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

    def _inter_spans_get(self, seq_type):
        assert seq_type in ('query', 'hit')
        attr_name = '%s_inter_ranges' % seq_type
        return [coord[1] - coord[0] for coord in getattr(self, attr_name)]

    def _hit_inter_spans_get(self):
        return self._inter_spans_get('hit')

    hit_inter_spans = property(fget=_hit_inter_spans_get,
        doc="""Lengths of regions between hit fragments""")

    def _query_inter_spans_get(self):
        return self._inter_spans_get('query')

    query_inter_spans = property(fget=_query_inter_spans_get,
        doc="""Lengths of regions between query fragments""")

    # shortcuts for fragments' properties #

    # bool check if there's more than one fragments
    is_fragmented = property(lambda self: len(self) > 1,
            doc="""Whether the HSP has more than one HSPFragment objects""")

    # first item properties with setters
    hit_description = fullcascade('hit_description',
            doc="""Description of the hit sequence""")

    query_description = fullcascade('query_description',
            doc="""Description of the query sequence""")

    hit_id = fullcascade('hit_id',
            doc="""ID of the hit sequence""")

    query_id = fullcascade('query_id',
            doc="""ID of the query sequence""")

    alphabet = fullcascade('alphabet',
            doc="""Alphabet used in hit and query SeqRecord objects""")

    # properties for single-fragment HSPs
    fragment = singleitem(
            doc="""HSPFragment object, first fragment""")

    hit = singleitem('hit',
            doc="""Hit sequence as a SeqRecord object, first fragment""")

    query = singleitem('query',
            doc="""Query sequence as a SeqRecord object, first fragment""")

    aln = singleitem('aln',
            doc="""Alignment of the first fragment as a MultipleSeqAlignment object""")

    aln_annotation = singleitem('aln_annotation',
            doc="""Dictionary of annotation(s) of the first fragment's alignment""")

    hit_features = singleitem('hit_features',
            doc="""Hit sequence features, first fragment""")

    query_features = singleitem('query_features',
            doc="""Query sequence features, first fragment""")

    hit_strand = singleitem('hit_strand',
            doc="""Hit strand orientation, first fragment""")

    query_strand = singleitem('query_strand',
            doc="""Query strand orientation, first fragment""")

    hit_frame = singleitem('hit_frame',
            doc="""Hit sequence reading frame, first fragment""")

    query_frame = singleitem('query_frame',
            doc="""Query sequence reading frame, first fragment""")

    # properties for multi-fragment HSPs
    fragments = allitems(doc="""List of all HSPFragment objects""")

    hit_all = allitems('hit',
            doc="""List of all fragments' hit sequences as SeqRecord objects""")

    query_all = allitems('query',
            doc="""List of all fragments' query sequences as SeqRecord objects""")

    aln_all = allitems('aln',
            doc="""List of all fragments' alignments as MultipleSeqAlignment objects""")

    aln_annotation_all = allitems('aln_annotation',
            doc="""Dictionary of annotation(s) of all fragments' alignments""")

    hit_features_all = allitems('hit_features',
            doc="""List of all hit sequence features""")

    query_features_all = allitems('query_features',
            doc="""List of all query sequence features""")

    hit_strand_all = allitems('hit_strand',
            doc="""List of all fragments' hit sequence strands""")

    query_strand_all = allitems('query_strand',
            doc="""List of all fragments' query sequence strands""")

    hit_frame_all = allitems('hit_frame',
            doc="""List of all fragments' hit sequence reading frames""")

    query_frame_all = allitems('query_frame',
            doc="""List of all fragments' query sequence reading frames""")

    hit_start_all = allitems('hit_start',
            doc="""List of all fragments' hit start coordinates""")

    query_start_all = allitems('query_start',
            doc="""List of all fragments' query start coordinates""")

    hit_end_all = allitems('hit_end',
            doc="""List of all fragments' hit end coordinates""")

    query_end_all = allitems('query_end',
            doc="""List of all fragments' query end coordinates""")

    hit_span_all = allitems('hit_span',
            doc="""List of all fragments' hit sequence size""")

    query_span_all = allitems('query_span',
            doc="""List of all fragments' query sequence size""")

    hit_range_all = allitems('hit_range',
            doc="""List of all fragments' hit start and end coordinates""")

    query_range_all = allitems('query_range',
            doc="""List of all fragments' query start and end coordinates""")


class HSPFragment(_BaseHSP):

    """Class representing a contiguous alignment of hit-query sequence.

    HSPFragment forms the core of any parsed search output file. Depending on
    the search output file format, it may contain the actual query and/or hit
    sequences that produces the search hits. These sequences are stored as
    SeqRecord objects (see SeqRecord):

    >>> from Bio import SearchIO
    >>> qresult = next(SearchIO.parse('Blast/mirna.xml', 'blast-xml'))
    >>> fragment = qresult[0][0][0]   # first hit, first hsp, first fragment
    >>> print(fragment)
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
    >>> print(fragment.query)
    ID: 33211
    Name: aligned query sequence
    Description: mir_1
    Number of features: 0
    Seq('CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTT...GGG', DNAAlphabet())

    # the hit sequence is a SeqRecord object as well
    >>> fragment.hit.__class__
    <class 'Bio.SeqRecord.SeqRecord'>
    >>> print(fragment.hit)
    ID: gi|262205317|ref|NR_030195.1|
    Name: aligned hit sequence
    Description: Homo sapiens microRNA 520b (MIR520B), microRNA
    Number of features: 0
    Seq('CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTT...GGG', DNAAlphabet())

    # when both query and hit are present, we get a MultipleSeqAlignment object
    >>> fragment.aln.__class__
    <class 'Bio.Align.MultipleSeqAlignment'>
    >>> print(fragment.aln)
    DNAAlphabet() alignment with 2 rows and 61 columns
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAG...GGG 33211
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAG...GGG gi|262205317|ref|NR_030195.1|

    """

    def __init__(self, hit_id='<unknown id>', query_id='<unknown id>',
            hit=None, query=None, alphabet=single_letter_alphabet):

        self._alphabet = alphabet
        self.aln_annotation = {}

        self._hit_id = hit_id
        self._query_id = query_id

        for seq_type in ('query', 'hit'):
            # query or hit attributes default attributes
            setattr(self, '_%s_description' % seq_type, '<unknown description>')
            setattr(self, '_%s_features' % seq_type, [])
            # query or hit attributes whose default attribute is None
            for attr in ('strand', 'frame', 'start', 'end'):
                setattr(self, '%s_%s' % (seq_type, attr), None)
            # self.query or self.hit
            if eval(seq_type):
                setattr(self, seq_type, eval(seq_type))
            else:
                setattr(self, seq_type, None)

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
        if self.aln is not None:
            obj = self.__class__(
                    hit_id=self.hit_id, query_id=self.query_id,
                    alphabet=self.alphabet)
            # transfer query and hit attributes
            # let SeqRecord handle feature slicing, then retrieve the sliced
            # features into the sliced HSPFragment
            if self.query is not None:
                obj.query = self.query[idx]
                obj.query_features = obj.query.features
            if self.hit is not None:
                obj.hit = self.hit[idx]
                obj.hit_features = obj.hit.features
            # description, strand, frame
            for attr in ('description', 'strand', 'frame'):
                for seq_type in ('hit', 'query'):
                    attr_name = '%s_%s' % (seq_type, attr)
                    self_val = getattr(self, attr_name)
                    setattr(obj, attr_name, self_val)
            # alignment annotation should be transferred, since we can compute
            # the resulting annotation
            obj.aln_annotation = {}
            for key, value in self.aln_annotation.items():
                assert len(value[idx]) == len(obj)
                obj.aln_annotation[key] = value[idx]
            return obj
        else:
            raise TypeError("Slicing for HSP objects without "
                    "alignment is not supported.")

    def _str_aln(self):
        lines = []
        # alignment length
        aln_span = getattr_str(self, 'aln_span')
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

            # similarity line
            simil = ''
            if 'similarity' in self.aln_annotation and \
                    isinstance(self.aln_annotation.get('similarity'), basestring):
                simil = self.aln_annotation['similarity']

            if self.aln_span <= 67:
                lines.append("%10s - %s" % ('Query', qseq))
                if simil:
                    lines.append("             %s" % simil)
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
                if simil:
                    lines.append("             %s%s%s" %
                            (simil[:59], cont, simil[-5:]))
                lines.append("%10s - %s%s%s" % ('Hit',
                                hseq[:59], cont, hseq[-5:]))

        return '\n'.join(lines)

    # sequence properties #
    def _set_seq(self, seq, seq_type):
        """Checks the given sequence for attribute setting

        :param seq: sequence to check
        :type seq: string or SeqRecord
        :param seq_type: sequence type
        :type seq_type: string, choice of 'hit' or 'query'

        """
        assert seq_type in ('hit', 'query')
        if seq is None:
            return seq # return immediately if seq is None
        else:
            if not isinstance(seq, (basestring, SeqRecord)):
                raise TypeError("%s sequence must be a string or a SeqRecord"
                        " object." % seq_type)
        # check length if the opposite sequence is not None
        opp_type = 'hit' if seq_type == 'query' else 'query'
        opp_seq = getattr(self, '_%s' % opp_type, None)
        if opp_seq is not None:
            if len(seq) != len(opp_seq):
                raise ValueError("Sequence lengths do not match. Expected: "
                        "%r (%s); found: %r (%s)." % (len(opp_seq), opp_type,
                        len(seq), seq_type))

        seq_id = getattr(self, '%s_id' % seq_type)
        seq_desc = getattr(self, '%s_description' % seq_type)
        seq_feats = getattr(self, '%s_features' % seq_type)
        seq_name = 'aligned %s sequence' % seq_type

        if isinstance(seq, SeqRecord):
            seq.id = seq_id
            seq.description = seq_desc
            seq.name = seq_name
            seq.features = seq_feats
            seq.seq.alphabet = self.alphabet
        elif isinstance(seq, basestring):
            seq = SeqRecord(Seq(seq, self.alphabet), id=seq_id, name=seq_name,
                    description=seq_desc, features=seq_feats)

        return seq

    def _hit_get(self):
        return self._hit

    def _hit_set(self, value):
        self._hit = self._set_seq(value, 'hit')

    hit = property(fget=_hit_get, fset=_hit_set,
            doc="""Hit sequence as a SeqRecord object, defaults to None""")

    def _query_get(self):
        return self._query

    def _query_set(self, value):
        self._query = self._set_seq(value, 'query')

    query = property(fget=_query_get, fset=_query_set,
            doc="""Query sequence as a SeqRecord object, defaults to None""")

    def _aln_get(self):
        if self.query is None and self.hit is None:
            return None
        elif self.hit is None:
            return MultipleSeqAlignment([self.query], self.alphabet)
        elif self.query is None:
            return MultipleSeqAlignment([self.hit], self.alphabet)
        else:
            return MultipleSeqAlignment([self.query, self.hit], self.alphabet)

    aln = property(fget=_aln_get,
            doc="""Query-hit alignment as a MultipleSeqAlignment object,
            defaults to None""")

    def _alphabet_get(self):
        return self._alphabet

    def _alphabet_set(self, value):
        self._alphabet = value
        try:
            self.query.seq.alphabet = value
        except AttributeError:
            pass
        try:
            self.hit.seq.alphabet = value
        except AttributeError:
            pass

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

    # id, description, and features properties #
    hit_description = fragcascade('description', 'hit',
            doc="""Hit sequence description""")

    query_description = fragcascade('description', 'query',
            doc="""Query sequence description""")

    hit_id = fragcascade('id', 'hit',
            doc="""Hit sequence ID""")

    query_id = fragcascade('id', 'query',
            doc="""Query sequence ID""")

    hit_features = fragcascade('features', 'hit',
            doc="""Hit sequence features""")

    query_features = fragcascade('features', 'query',
            doc="""Query sequence features""")

    # strand properties #
    def _prep_strand(self, strand):
        # follow SeqFeature's convention
        if strand not in (-1, 0, 1, None):
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

    # frame properties #
    def _prep_frame(self, frame):
        if frame not in (-3, -2, -1, 0, 1, 2, 3, None):
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

    # coordinate properties #
    def _prep_coord(self, coord, opp_coord_name, op):
        # coord must either be None or int
        if coord is None:
            return coord
        assert isinstance(coord, int)
        # try to get opposite coordinate, if it's not present, return
        try:
            opp_coord = getattr(self, opp_coord_name)
        except AttributeError:
            return coord
        # if opposite coordinate is None, return
        if opp_coord is None:
            return coord
        # otherwise compare it to coord ('>=' or '<=')
        else:
            assert op(coord, opp_coord)
        return coord

    def _hit_start_get(self):
        return self._hit_start

    def _hit_start_set(self, value):
        self._hit_start = self._prep_coord(value, 'hit_end', le)

    hit_start = property(fget=_hit_start_get, fset=_hit_start_set,
            doc="""Hit sequence start coordinate, defaults to None""")

    def _query_start_get(self):
        return self._query_start

    def _query_start_set(self, value):
        self._query_start = self._prep_coord(value, 'query_end', le)

    query_start = property(fget=_query_start_get, fset=_query_start_set,
            doc="""Query sequence start coordinate, defaults to None""")

    def _hit_end_get(self):
        return self._hit_end

    def _hit_end_set(self, value):
        self._hit_end = self._prep_coord(value, 'hit_start', ge)

    hit_end = property(fget=_hit_end_get, fset=_hit_end_set,
            doc="""Hit sequence start coordinate, defaults to None""")

    def _query_end_get(self):
        return self._query_end

    def _query_end_set(self, value):
        self._query_end = self._prep_coord(value, 'query_start', ge)

    query_end = property(fget=_query_end_get, fset=_query_end_set,
            doc="""Query sequence end coordinate, defaults to None""")

    # coordinate-dependent properties #
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


# if not used as a module, run the doctest
if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()
