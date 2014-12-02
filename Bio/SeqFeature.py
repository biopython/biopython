# Copyright 2000-2003 Jeff Chang.
# Copyright 2001-2008 Brad Chapman.
# Copyright 2005-2012 by Peter Cock.
# Copyright 2006-2009 Michiel de Hoon.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Represent a Sequence Feature holding info about a part of a sequence.

This is heavily modeled after the Biocorba SeqFeature objects, and
may be pretty biased towards GenBank stuff since I'm writing it
for the GenBank parser output...

What's here:

Base class to hold a Feature
----------------------------

classes:

    - SeqFeature

Hold information about a Reference
----------------------------------

This is an attempt to create a General class to hold Reference type
information.

classes:

    - Reference

Specify locations of a feature on a Sequence
--------------------------------------------

This aims to handle, in Ewan's words, 'the dreaded fuzziness issue' in
much the same way as Biocorba. This has the advantages of allowing us
to handle fuzzy stuff in case anyone needs it, and also be compatible
with Biocorba.

classes:

    - FeatureLocation - Specify the start and end location of a feature.
    - CompoundLocation - Collection of FeatureLocation objects (for joins etc).

    - ExactPosition - Specify the position as being exact.
    - WithinPosition - Specify a position occuring within some range.
    - BetweenPosition - Specify a position occuring between a range (OBSOLETE?).
    - BeforePosition - Specify the position as being found before some base.
    - AfterPosition - Specify the position as being found after some base.
    - OneOfPosition - Specify a position where the location can be multiple positions.
    - UnknownPosition - Represents missing information like '?' in UniProt.
"""

from __future__ import print_function

from Bio.Seq import MutableSeq, reverse_complement

__docformat__ = "restructuredtext en"

class SeqFeature(object):
    """Represent a Sequence Feature on an object.

    Attributes:

        - location - the location of the feature on the sequence (FeatureLocation)
        - type - the specified type of the feature (ie. CDS, exon, repeat...)
        - location_operator - a string specifying how this SeqFeature may
          be related to others. For example, in the example GenBank feature
          shown below, the location_operator would be "join". This is a proxy
          for feature.location.operator and only applies to compound locations.
        - strand - A value specifying on which strand (of a DNA sequence, for
          instance) the feature deals with. 1 indicates the plus strand, -1
          indicates the minus strand, 0 indicates stranded but unknown (? in GFF3),
          while the default of None indicates that strand doesn't apply (dot in GFF3,
          e.g. features on proteins). Note this is a shortcut for accessing the
          strand property of the feature's location.
        - id - A string identifier for the feature.
        - ref - A reference to another sequence. This could be an accession
          number for some different sequence. Note this is a shortcut for the
          reference property of the feature's location.
        - ref_db - A different database for the reference accession number.
          Note this is a shortcut for the reference property of the location
        - qualifiers - A dictionary of qualifiers on the feature. These are
          analogous to the qualifiers from a GenBank feature table. The keys of
          the dictionary are qualifier names, the values are the qualifier
          values.
        - sub_features - Obsolete list of additional SeqFeatures which was
          used for holding compound locations (e.g. joins in GenBank/EMBL).
          This is now superceded by a CompoundFeatureLocation as the location,
          and should not be used (DEPRECATED).
    """

    def __init__(self, location=None, type='', location_operator='',
                 strand=None, id="<unknown id>",
                 qualifiers=None, sub_features=None,
                 ref=None, ref_db=None):
        """Initialize a SeqFeature on a Sequence.

        location can either be a FeatureLocation (with strand argument also
        given if required), or None.

        e.g. With no strand, on the forward strand, and on the reverse strand:

        >>> from Bio.SeqFeature import SeqFeature, FeatureLocation
        >>> f1 = SeqFeature(FeatureLocation(5, 10), type="domain")
        >>> f1.strand == f1.location.strand == None
        True
        >>> f2 = SeqFeature(FeatureLocation(7, 110, strand=1), type="CDS")
        >>> f2.strand == f2.location.strand == +1
        True
        >>> f3 = SeqFeature(FeatureLocation(9, 108, strand=-1), type="CDS")
        >>> f3.strand == f3.location.strand == -1
        True

        An invalid strand will trigger an exception:

        >>> f4 = SeqFeature(FeatureLocation(50, 60), strand=2)
        Traceback (most recent call last):
           ...
        ValueError: Strand should be +1, -1, 0 or None, not 2

        Similarly if set via the FeatureLocation directly:

        >>> loc4 = FeatureLocation(50, 60, strand=2)
        Traceback (most recent call last):
           ...
        ValueError: Strand should be +1, -1, 0 or None, not 2

        For exact start/end positions, an integer can be used (as shown above)
        as shorthand for the ExactPosition object. For non-exact locations, the
        FeatureLocation must be specified via the appropriate position objects.

        Note that the strand, ref and ref_db arguments to the SeqFeature are
        now obsolete and will be deprecated in a future release (which will
        give warning messages) and later removed. Set them via the location
        object instead.

        Note that location_operator and sub_features arguments can no longer
        be used, instead do this via the CompoundLocation object.
        """
        if location is not None and not isinstance(location, FeatureLocation) \
                and not isinstance(location, CompoundLocation):
            raise TypeError(
                "FeatureLocation, CompoundLocation (or None) required for the location")
        self.location = location
        self.type = type
        if location_operator:
            # TODO - Deprecation warning
            self.location_operator = location_operator
        if strand is not None:
            # TODO - Deprecation warning
            self.strand = strand
        self.id = id
        if qualifiers is None:
            qualifiers = {}
        self.qualifiers = qualifiers
        if sub_features is None:
            sub_features = []
        else:
            import warnings
            from Bio import BiopythonDeprecationWarning
            warnings.warn("Rather than sub_features, use a CompoundFeatureLocation",
                          BiopythonDeprecationWarning)
        self._sub_features = sub_features
        if ref is not None:
            # TODO - Deprecation warning
            self.ref = ref
        if ref_db is not None:
            # TODO - Deprecation warning
            self.ref_db = ref_db

    def _get_sub_features(self):
        if self._sub_features:
            import warnings
            from Bio import BiopythonDeprecationWarning
            warnings.warn("Rather using f.sub_features, f.location should be a CompoundFeatureLocation",
                          BiopythonDeprecationWarning)
        return self._sub_features

    def _set_sub_features(self, value):
        if value:
            import warnings
            from Bio import BiopythonDeprecationWarning
            warnings.warn("Rather than f.sub_features, use a CompoundFeatureLocation for f.location",
                          BiopythonDeprecationWarning)
        self._sub_features = value
    sub_features = property(fget=_get_sub_features, fset=_set_sub_features,
                            doc="Obsolete representation of compound locations (DEPRECATED).")

    def _get_strand(self):
        return self.location.strand

    def _set_strand(self, value):
        try:
            self.location.strand = value
        except AttributeError:
            if self.location is None:
                if value is not None:
                    raise ValueError("Can't set strand without a location.")
            else:
                raise

    strand = property(fget=_get_strand, fset=_set_strand,
                      doc="""Feature's strand

                          This is a shortcut for feature.location.strand
                          """)

    def _get_ref(self):
        try:
            return self.location.ref
        except AttributeError:
            return None

    def _set_ref(self, value):
        try:
            self.location.ref = value
        except AttributeError:
            if self.location is None:
                if value is not None:
                    raise ValueError("Can't set ref without a location.")
            else:
                raise
    ref = property(fget=_get_ref, fset=_set_ref,
                   doc="""Feature location reference (e.g. accession).

                       This is a shortcut for feature.location.ref
                       """)

    def _get_ref_db(self):
        try:
            return self.location.ref_db
        except AttributeError:
            return None

    def _set_ref_db(self, value):
        self.location.ref_db = value
    ref_db = property(fget=_get_ref_db, fset=_set_ref_db,
                      doc="""Feature location reference's database.

                          This is a shortcut for feature.location.ref_db
                          """)

    def _get_location_operator(self):
        try:
            return self.location.operator
        except AttributeError:
            return None

    def _set_location_operator(self, value):
        if value:
            if isinstance(self.location, CompoundLocation):
                self.location.operator = value
            elif self.location is None:
                raise ValueError(
                    "Location is None so can't set its operator (to %r)" % value)
            else:
                raise ValueError(
                    "Only CompoundLocation gets an operator (%r)" % value)
    location_operator = property(fget=_get_location_operator, fset=_set_location_operator,
                                 doc="Location operator for compound locations (e.g. join).")

    def __repr__(self):
        """A string representation of the record for debugging."""
        answer = "%s(%s" % (self.__class__.__name__, repr(self.location))
        if self.type:
            answer += ", type=%s" % repr(self.type)
        if self.location_operator:
            answer += ", location_operator=%s" % repr(self.location_operator)
        if self.id and self.id != "<unknown id>":
            answer += ", id=%s" % repr(self.id)
        if self.ref:
            answer += ", ref=%s" % repr(self.ref)
        if self.ref_db:
            answer += ", ref_db=%s" % repr(self.ref_db)
        answer += ")"
        return answer

    def __str__(self):
        """A readable summary of the feature intended to be printed to screen.
        """
        out = "type: %s\n" % self.type
        out += "location: %s\n" % self.location
        if self.id and self.id != "<unknown id>":
            out += "id: %s\n" % self.id
        out += "qualifiers:\n"
        for qual_key in sorted(self.qualifiers):
            out += "    Key: %s, Value: %s\n" % (qual_key,
                                                 self.qualifiers[qual_key])
        # TODO - Remove this from __str__ since deprecated
        if len(self._sub_features) != 0:
            out += "Sub-Features\n"
            for sub_feature in self._sub_features:
                out += "%s\n" % sub_feature
        return out

    def _shift(self, offset):
        """Returns a copy of the feature with its location shifted (PRIVATE).

        The annotation qaulifiers are copied."""
        answer = SeqFeature(location=self.location._shift(offset),
                            type=self.type,
                            location_operator=self.location_operator,
                            id=self.id,
                            qualifiers=dict(self.qualifiers.items()))
        # This is to avoid the deprecation warning:
        answer._sub_features = [f._shift(offset) for f in self._sub_features]
        return answer

    def _flip(self, length):
        """Returns a copy of the feature with its location flipped (PRIVATE).

        The argument length gives the length of the parent sequence. For
        example a location 0..20 (+1 strand) with parent length 30 becomes
        after flipping 10..30 (-1 strand). Strandless (None) or unknown
        strand (0) remain like that - just their end points are changed.

        The annotation qaulifiers are copied.
        """
        answer = SeqFeature(location=self.location._flip(length),
                            type=self.type,
                            location_operator=self.location_operator,
                            id=self.id,
                            qualifiers=dict(self.qualifiers.items()))
        # This is to avoid the deprecation warning:
        answer._sub_features = [f._flip(length)
                                for f in self._sub_features[::-1]]
        return answer

    def extract(self, parent_sequence):
        """Extract feature sequence from the supplied parent sequence.

        The parent_sequence can be a Seq like object or a string, and will
        generally return an object of the same type. The exception to this is
        a MutableSeq as the parent sequence will return a Seq object.

        This should cope with complex locations including complements, joins
        and fuzzy positions. Even mixed strand features should work! This
        also covers features on protein sequences (e.g. domains), although
        here reverse strand features are not permitted.

        >>> from Bio.Seq import Seq
        >>> from Bio.Alphabet import generic_protein
        >>> from Bio.SeqFeature import SeqFeature, FeatureLocation
        >>> seq = Seq("MKQHKAMIVALIVICITAVVAAL", generic_protein)
        >>> f = SeqFeature(FeatureLocation(8, 15), type="domain")
        >>> f.extract(seq)
        Seq('VALIVIC', ProteinAlphabet())

        Note - currently only sub-features of type "join" are supported.
        """
        return self.location.extract(parent_sequence)

    # Python 3:
    def __bool__(self):
        """Boolean value of an instance of this class (True).

        This behaviour is for backwards compatibility, since until the
        __len__ method was added, a SeqFeature always evaluated as True.

        Note that in comparison, Seq objects, strings, lists, etc, will all
        evaluate to False if they have length zero.

        WARNING: The SeqFeature may in future evaluate to False when its
        length is zero (in order to better match normal python behaviour)!
        """
        return True

    # Python 2:
    __nonzero__ = __bool__

    def __len__(self):
        """Returns the length of the region described by a feature.

        >>> from Bio.Seq import Seq
        >>> from Bio.Alphabet import generic_protein
        >>> from Bio.SeqFeature import SeqFeature, FeatureLocation
        >>> seq = Seq("MKQHKAMIVALIVICITAVVAAL", generic_protein)
        >>> f = SeqFeature(FeatureLocation(8, 15), type="domain")
        >>> len(f)
        7
        >>> f.extract(seq)
        Seq('VALIVIC', ProteinAlphabet())
        >>> len(f.extract(seq))
        7

        This is a proxy for taking the length of the feature's location:

        >>> len(f.location)
        7

        For simple features this is the same as the region spanned (end
        position minus start position using Pythonic counting). However, for
        a compound location (e.g. a CDS as the join of several exons) the
        gaps are not counted (e.g. introns). This ensures that len(f) matches
        len(f.extract(parent_seq)), and also makes sure things work properly
        with features wrapping the origin etc.
        """
        return len(self.location)

    def __iter__(self):
        """Iterate over the parent positions within the feature.

        The iteration order is strand aware, and can be thought of as moving
        along the feature using the parent sequence coordinates:

        >>> from Bio.SeqFeature import SeqFeature, FeatureLocation
        >>> f = SeqFeature(FeatureLocation(5, 10), type="domain", strand=-1)
        >>> len(f)
        5
        >>> for i in f: print(i)
        9
        8
        7
        6
        5
        >>> list(f)
        [9, 8, 7, 6, 5]

        This is a proxy for iterating over the location,

        >>> list(f.location)
        [9, 8, 7, 6, 5]
        """
        return iter(self.location)

    def __contains__(self, value):
        """Check if an integer position is within the feature.

        >>> from Bio.SeqFeature import SeqFeature, FeatureLocation
        >>> f = SeqFeature(FeatureLocation(5, 10), type="domain", strand=-1)
        >>> len(f)
        5
        >>> [i for i in range(15) if i in f]
        [5, 6, 7, 8, 9]

        For example, to see which features include a SNP position, you could
        use this:

        >>> from Bio import SeqIO
        >>> record = SeqIO.read("GenBank/NC_000932.gb", "gb")
        >>> for f in record.features:
        ...     if 1750 in f:
        ...         print("%s %s" % (f.type, f.location))
        source [0:154478](+)
        gene [1716:4347](-)
        tRNA join{[4310:4347](-), [1716:1751](-)}

        Note that for a feature defined as a join of several subfeatures (e.g.
        the union of several exons) the gaps are not checked (e.g. introns).
        In this example, the tRNA location is defined in the GenBank file as
        complement(join(1717..1751,4311..4347)), so that position 1760 falls
        in the gap:

        >>> for f in record.features:
        ...     if 1760 in f:
        ...         print("%s %s" % (f.type, f.location))
        source [0:154478](+)
        gene [1716:4347](-)

        Note that additional care may be required with fuzzy locations, for
        example just before a BeforePosition:

        >>> from Bio.SeqFeature import SeqFeature, FeatureLocation
        >>> from Bio.SeqFeature import BeforePosition
        >>> f = SeqFeature(FeatureLocation(BeforePosition(3), 8), type="domain")
        >>> len(f)
        5
        >>> [i for i in range(10) if i in f]
        [3, 4, 5, 6, 7]

        Note that is is a proxy for testing membership on the location.

        >>> [i for i in range(10) if i in f.location]
        [3, 4, 5, 6, 7]
        """
        return value in self.location


# --- References


# TODO -- Will this hold PubMed and Medline information decently?
class Reference(object):
    """Represent a Generic Reference object.

    Attributes:
    o location - A list of Location objects specifying regions of
    the sequence that the references correspond to. If no locations are
    specified, the entire sequence is assumed.
    o authors - A big old string, or a list split by author, of authors
    for the reference.
    o title - The title of the reference.
    o journal - Journal the reference was published in.
    o medline_id - A medline reference for the article.
    o pubmed_id - A pubmed reference for the article.
    o comment - A place to stick any comments about the reference.
    """

    def __init__(self):
        self.location = []
        self.authors = ''
        self.consrtm = ''
        self.title = ''
        self.journal = ''
        self.medline_id = ''
        self.pubmed_id = ''
        self.comment = ''

    def __str__(self):
        """Output an informative string for debugging.
        """
        out = ""
        for single_location in self.location:
            out += "location: %s\n" % single_location
        out += "authors: %s\n" % self.authors
        if self.consrtm:
            out += "consrtm: %s\n" % self.consrtm
        out += "title: %s\n" % self.title
        out += "journal: %s\n" % self.journal
        out += "medline id: %s\n" % self.medline_id
        out += "pubmed id: %s\n" % self.pubmed_id
        out += "comment: %s\n" % self.comment
        return out

    def __repr__(self):
        # TODO - Update this is __init__ later accpets values
        return "%s(title=%s, ...)" % (self.__class__.__name__,
                                      repr(self.title))


# --- Handling feature locations

class FeatureLocation(object):
    """Specify the location of a feature along a sequence.

    The FeatureLocation is used for simple continous features, which can
    be described as running from a start position to and end position
    (optionally with a strand and reference information).  More complex
    locations made up from several non-continuous parts (e.g. a coding
    sequence made up of several exons) are currently described using a
    SeqFeature with sub-features.

    Note that the start and end location numbering follow Python's scheme,
    thus a GenBank entry of 123..150 (one based counting) becomes a location
    of [122:150] (zero based counting).

    >>> from Bio.SeqFeature import FeatureLocation
    >>> f = FeatureLocation(122, 150)
    >>> print(f)
    [122:150]
    >>> print(f.start)
    122
    >>> print(f.end)
    150
    >>> print(f.strand)
    None

    Note the strand defaults to None. If you are working with nucleotide
    sequences you'd want to be explicit if it is the forward strand:

    >>> from Bio.SeqFeature import FeatureLocation
    >>> f = FeatureLocation(122, 150, strand=+1)
    >>> print(f)
    [122:150](+)
    >>> print(f.strand)
    1

    Note that for a parent sequence of length n, the FeatureLocation
    start and end must satisfy the inequality 0 <= start <= end <= n.
    This means even for features on the reverse strand of a nucleotide
    sequence, we expect the 'start' coordinate to be less than the
    'end'.

    >>> from Bio.SeqFeature import FeatureLocation
    >>> r = FeatureLocation(122, 150, strand=-1)
    >>> print(r)
    [122:150](-)
    >>> print(r.start)
    122
    >>> print(r.end)
    150
    >>> print(r.strand)
    -1

    i.e. Rather than thinking of the 'start' and 'end' biologically in a
    strand aware manor, think of them as the 'left most' or 'minimum'
    boundary, and the 'right most' or 'maximum' boundary of the region
    being described. This is particularly important with compound
    locations describing non-continuous regions.

    In the example above we have used standard exact positions, but there
    are also specialised position objects used to represent fuzzy positions
    as well, for example a GenBank location like complement(<123..150)
    would use a BeforePosition object for the start.
    """

    def __init__(self, start, end, strand=None, ref=None, ref_db=None):
        """Specify the start, end, strand etc of a sequence feature.

        start and end arguments specify the values where the feature begins
        and ends. These can either by any of the ``*Position`` objects that
        inherit from AbstractPosition, or can just be integers specifying the
        position. In the case of integers, the values are assumed to be
        exact and are converted in ExactPosition arguments. This is meant
        to make it easy to deal with non-fuzzy ends.

        i.e. Short form:

        >>> from Bio.SeqFeature import FeatureLocation
        >>> loc = FeatureLocation(5, 10, strand=-1)
        >>> print(loc)
        [5:10](-)

        Explicit form:

        >>> from Bio.SeqFeature import FeatureLocation, ExactPosition
        >>> loc = FeatureLocation(ExactPosition(5), ExactPosition(10), strand=-1)
        >>> print(loc)
        [5:10](-)

        Other fuzzy positions are used similarly,

        >>> from Bio.SeqFeature import FeatureLocation
        >>> from Bio.SeqFeature import BeforePosition, AfterPosition
        >>> loc2 = FeatureLocation(BeforePosition(5), AfterPosition(10), strand=-1)
        >>> print(loc2)
        [<5:>10](-)

        For nucleotide features you will also want to specify the strand,
        use 1 for the forward (plus) strand, -1 for the reverse (negative)
        strand, 0 for stranded but strand unknown (? in GFF3), or None for
        when the strand does not apply (dot in GFF3), e.g. features on
        proteins.

        >>> loc = FeatureLocation(5, 10, strand=+1)
        >>> print(loc)
        [5:10](+)
        >>> print(loc.strand)
        1

        Normally feature locations are given relative to the parent
        sequence you are working with, but an explicit accession can
        be given with the optional ref and db_ref strings:

        >>> loc = FeatureLocation(105172, 108462, ref="AL391218.9", strand=1)
        >>> print(loc)
        AL391218.9[105172:108462](+)
        >>> print(loc.ref)
        AL391218.9

        """
        # TODO - Check 0 <= start <= end (<= length of reference)
        if isinstance(start, AbstractPosition):
            self._start = start
        elif isinstance(start, int) or isinstance(start, long):
            self._start = ExactPosition(start)
        else:
            raise TypeError("start=%r %s" % (start, type(start)))
        if isinstance(end, AbstractPosition):
            self._end = end
        elif isinstance(end, int) or isinstance(end, long):
            self._end = ExactPosition(end)
        else:
            raise TypeError("end=%r %s" % (end, type(end)))
        self.strand = strand
        self.ref = ref
        self.ref_db = ref_db

    def _get_strand(self):
        return self._strand

    def _set_strand(self, value):
        if value not in [+1, -1, 0, None]:
            raise ValueError("Strand should be +1, -1, 0 or None, not %r"
                             % value)
        self._strand = value

    strand = property(fget=_get_strand, fset=_set_strand,
                      doc="Strand of the location (+1, -1, 0 or None).")

    def __str__(self):
        """Returns a representation of the location (with python counting).

        For the simple case this uses the python splicing syntax, [122:150]
        (zero based counting) which GenBank would call 123..150 (one based
        counting).
        """
        answer = "[%s:%s]" % (self._start, self._end)
        if self.ref and self.ref_db:
            answer = "%s:%s%s" % (self.ref_db, self.ref, answer)
        elif self.ref:
            answer = self.ref + answer
        # Is ref_db without ref meaningful?
        if self.strand is None:
            return answer
        elif self.strand == +1:
            return answer + "(+)"
        elif self.strand == -1:
            return answer + "(-)"
        else:
            # strand = 0, stranded but strand unknown, ? in GFF3
            return answer + "(?)"

    def __repr__(self):
        """A string representation of the location for debugging."""
        optional = ""
        if self.strand is not None:
            optional += ", strand=%r" % self.strand
        if self.ref is not None:
            optional += ", ref=%r" % self.ref
        if self.ref_db is not None:
            optional += ", ref_db=%r" % self.ref_db
        return "%s(%r, %r%s)" \
            % (self.__class__.__name__, self.start, self.end, optional)

    def __add__(self, other):
        """Combine location with another feature location, or shift it.

        You can add two feature locations to make a join CompoundLocation:

        >>> from Bio.SeqFeature import FeatureLocation
        >>> f1 = FeatureLocation(5, 10)
        >>> f2 = FeatureLocation(20, 30)
        >>> combined = f1 + f2
        >>> print(combined)
        join{[5:10], [20:30]}

        This is thus equivalent to:

        >>> from Bio.SeqFeature import CompoundLocation
        >>> join = CompoundLocation([f1, f2])
        >>> print(join)
        join{[5:10], [20:30]}

        You can also use sum(...) in this way:

        >>> join = sum([f1, f2])
        >>> print(join)
        join{[5:10], [20:30]}

        Furthermore, you can combine a FeatureLocation with a CompoundLocation
        in this way.

        Separately, adding an integer will give a new FeatureLocation with
        its start and end offset by that amount. For example:

        >>> print(f1)
        [5:10]
        >>> print(f1 + 100)
        [105:110]
        >>> print(200 + f1)
        [205:210]

        This can be useful when editing annotation.
        """
        if isinstance(other, FeatureLocation):
            return CompoundLocation([self, other])
        elif isinstance(other, int):
            return self._shift(other)
        else:
            # This will allow CompoundLocation's __radd__ to be called:
            return NotImplemented

    def __radd__(self, other):
        if isinstance(other, int):
            return self._shift(other)
        else:
            return NotImplemented

    def __nonzero__(self):
        """Returns True regardless of the length of the feature.

        This behaviour is for backwards compatibility, since until the
        __len__ method was added, a FeatureLocation always evaluated as True.

        Note that in comparison, Seq objects, strings, lists, etc, will all
        evaluate to False if they have length zero.

        WARNING: The FeatureLocation may in future evaluate to False when its
        length is zero (in order to better match normal python behaviour)!
        """
        return True

    def __len__(self):
        """Returns the length of the region described by the FeatureLocation.

        Note that extra care may be needed for fuzzy locations, e.g.

        >>> from Bio.SeqFeature import FeatureLocation
        >>> from Bio.SeqFeature import BeforePosition, AfterPosition
        >>> loc = FeatureLocation(BeforePosition(5), AfterPosition(10))
        >>> len(loc)
        5
        """
        return int(self._end) - int(self._start)

    def __contains__(self, value):
        """Check if an integer position is within the FeatureLocation.

        Note that extra care may be needed for fuzzy locations, e.g.

        >>> from Bio.SeqFeature import FeatureLocation
        >>> from Bio.SeqFeature import BeforePosition, AfterPosition
        >>> loc = FeatureLocation(BeforePosition(5), AfterPosition(10))
        >>> len(loc)
        5
        >>> [i for i in range(15) if i in loc]
        [5, 6, 7, 8, 9]
        """
        if not isinstance(value, int):
            raise ValueError("Currently we only support checking for integer "
                             "positions being within a FeatureLocation.")
        if value < self._start or value >= self._end:
            return False
        else:
            return True

    def __iter__(self):
        """Iterate over the parent positions within the FeatureLocation.

        >>> from Bio.SeqFeature import FeatureLocation
        >>> from Bio.SeqFeature import BeforePosition, AfterPosition
        >>> loc = FeatureLocation(BeforePosition(5), AfterPosition(10))
        >>> len(loc)
        5
        >>> for i in loc: print(i)
        5
        6
        7
        8
        9
        >>> list(loc)
        [5, 6, 7, 8, 9]
        >>> [i for i in range(15) if i in loc]
        [5, 6, 7, 8, 9]

        Note this is strand aware:

        >>> loc = FeatureLocation(BeforePosition(5), AfterPosition(10), strand = -1)
        >>> list(loc)
        [9, 8, 7, 6, 5]
        """
        if self.strand == -1:
            for i in range(self._end - 1, self._start - 1, -1):
                yield i
        else:
            for i in range(self._start, self._end):
                yield i

    def _shift(self, offset):
        """Returns a copy of the location shifted by the offset (PRIVATE)."""
        # TODO - What if offset is a fuzzy position?
        if self.ref or self.ref_db:
            # TODO - Return self?
            raise ValueError("Feature references another sequence.")
        return FeatureLocation(start=self._start._shift(offset),
                               end=self._end._shift(offset),
                               strand=self.strand)

    def _flip(self, length):
        """Returns a copy of the location after the parent is reversed (PRIVATE)."""
        if self.ref or self.ref_db:
            # TODO - Return self?
            raise ValueError("Feature references another sequence.")
        # Note this will flip the start and end too!
        if self.strand == +1:
            flip_strand = -1
        elif self.strand == -1:
            flip_strand = +1
        else:
            # 0 or None
            flip_strand = self.strand
        return FeatureLocation(start=self._end._flip(length),
                               end=self._start._flip(length),
                               strand=flip_strand)

    @property
    def parts(self):
        """Read only list of parts (always one, the Feature Location).

        This is a convience property allowing you to write code handling
        both simple FeatureLocation objects (with one part) and more complex
        CompoundLocation objects (with multiple parts) interchangably.
        """
        return [self]

    @property
    def start(self):
        """Start location (integer like, possibly a fuzzy position, read only)."""
        return self._start

    @property
    def end(self):
        """End location (integer like, possibly a fuzzy position, read only)."""
        return self._end

    @property
    def nofuzzy_start(self):
        """Start position (integer, approximated if fuzzy, read only) (OBSOLETE).

        This is now an alias for int(feature.start), which should be
        used in preference -- unless you are trying to support old
        versions of Biopython.
        """
        try:
            return int(self._start)
        except TypeError:
            if isinstance(self._start, UnknownPosition):
                return None
            raise

    @property
    def nofuzzy_end(self):
        """End position (integer, approximated if fuzzy, read only) (OBSOLETE).

        This is now an alias for int(feature.end), which should be
        used in preference -- unless you are trying to support old
        versions of Biopython.
        """
        try:
            return int(self._end)
        except TypeError:
            if isinstance(self._end, UnknownPosition):
                return None
            raise

    def extract(self, parent_sequence):
        """Extract feature sequence from the supplied parent sequence."""
        if self.ref or self.ref_db:
            # TODO - Take a dictionary as an optional argument?
            raise ValueError("Feature references another sequence.")
        if isinstance(parent_sequence, MutableSeq):
            # This avoids complications with reverse complements
            # (the MutableSeq reverse complement acts in situ)
            parent_sequence = parent_sequence.toseq()
        f_seq = parent_sequence[self.nofuzzy_start:self.nofuzzy_end]
        if self.strand == -1:
            try:
                f_seq = f_seq.reverse_complement()
            except AttributeError:
                assert isinstance(f_seq, str)
                f_seq = reverse_complement(f_seq)
        return f_seq


class CompoundLocation(object):
    """For handling joins etc where a feature location has several parts."""

    def __init__(self, parts, operator="join"):
        """Create a compound location with several parts.

        >>> from Bio.SeqFeature import FeatureLocation, CompoundLocation
        >>> f1 = FeatureLocation(10, 40, strand=+1)
        >>> f2 = FeatureLocation(50, 59, strand=+1)
        >>> f = CompoundLocation([f1, f2])
        >>> len(f) == len(f1) + len(f2) == 39 == len(list(f))
        True
        >>> print(f.operator)
        join
        >>> 5 in f
        False
        >>> 15 in f
        True
        >>> f.strand
        1

        Notice that the strand of the compound location is computed
        automatically - in the case of mixed strands on the sub-locations
        the overall strand is set to None.

        >>> f = CompoundLocation([FeatureLocation(3, 6, strand=+1),
        ...                       FeatureLocation(10, 13, strand=-1)])
        >>> print(f.strand)
        None
        >>> len(f)
        6
        >>> list(f)
        [3, 4, 5, 12, 11, 10]

        The example above doing list(f) iterates over the coordinates within the
        feature. This allows you to use max and min on the location, to find the
        range covered:

        >>> min(f)
        3
        >>> max(f)
        12

        More generally, you can use the compound location's start and end which
        give the full range covered, 0 <= start <= end <= full sequence length.

        >>> f.start == min(f)
        True
        >>> f.end == max(f) + 1
        True

        This is consistent with the behaviour of the simple FeatureLocation for
        a single region, where again the 'start' and 'end' do not necessarily
        give the biological start and end, but rather the 'minimal' and 'maximal'
        coordinate boundaries.

        Note that adding locations provides a more intuitive method of
        construction:

        >>> f = FeatureLocation(3, 6, strand=+1) + FeatureLocation(10, 13, strand=-1)
        >>> len(f)
        6
        >>> list(f)
        [3, 4, 5, 12, 11, 10]
        """
        self.operator = operator
        self.parts = list(parts)
        for loc in self.parts:
            if not isinstance(loc, FeatureLocation):
                raise ValueError("CompoundLocation should be given a list of "
                                 "FeatureLocation objects, not %s" % loc.__class__)
        if len(parts) < 2:
            raise ValueError(
                "CompoundLocation should have at least 2 parts, not %r" % parts)

    def __str__(self):
        """Returns a representation of the location (with python counting)."""
        return "%s{%s}" % (self.operator, ", ".join(str(loc) for loc in self.parts))

    def __repr__(self):
        """String representation of the location for debugging."""
        return "%s(%r, %r)" % (self.__class__.__name__,
                               self.parts, self.operator)

    def _get_strand(self):
        # Historically a join on the reverse strand has been represented
        # in Biopython with both the parent SeqFeature and its children
        # (the exons for a CDS) all given a strand of -1.  Likewise, for
        # a join feature on the forward strand they all have strand +1.
        # However, we must also consider evil mixed strand examples like
        # this, join(complement(69611..69724),139856..140087,140625..140650)
        if len(set(loc.strand for loc in self.parts)) == 1:
            return self.parts[0].strand
        else:
            return None  # i.e. mixed strands

    def _set_strand(self, value):
        # Should this be allowed/encouraged?
        for loc in self.parts:
            loc.strand = value
    strand = property(fget=_get_strand, fset=_set_strand,
                      doc="""Overall strand of the compound location.

        If all the parts have the same strand, that is returned. Otherwise
        for mixed strands, this returns None.

        >>> from Bio.SeqFeature import FeatureLocation, CompoundLocation
        >>> f1 = FeatureLocation(15, 17, strand=1)
        >>> f2 = FeatureLocation(20, 30, strand=-1)
        >>> f = f1 + f2
        >>> f1.strand
        1
        >>> f2.strand
        -1
        >>> f.strand
        >>> f.strand is None
        True

        If you set the strand of a CompoundLocation, this is applied to
        all the parts - use with caution:

        >>> f.strand = 1
        >>> f1.strand
        1
        >>> f2.strand
        1
        >>> f.strand
        1

        """)

    def __add__(self, other):
        """Combine locations, or shift the location by an integer offset.

        >>> from Bio.SeqFeature import FeatureLocation, CompoundLocation
        >>> f1 = FeatureLocation(15, 17) + FeatureLocation(20, 30)
        >>> print(f1)
        join{[15:17], [20:30]}

        You can add another FeatureLocation:

        >>> print(f1 + FeatureLocation(40, 50))
        join{[15:17], [20:30], [40:50]}
        >>> print(FeatureLocation(5, 10) + f1)
        join{[5:10], [15:17], [20:30]}

        You can also add another CompoundLocation:

        >>> f2 = FeatureLocation(40, 50) + FeatureLocation(60, 70)
        >>> print(f2)
        join{[40:50], [60:70]}
        >>> print(f1 + f2)
        join{[15:17], [20:30], [40:50], [60:70]}

        Also, as with the FeatureLocation, adding an integer shifts the
        location's co-ordinates by that offset:

        >>> print(f1 + 100)
        join{[115:117], [120:130]}
        >>> print(200 + f1)
        join{[215:217], [220:230]}
        >>> print(f1 + (-5))
        join{[10:12], [15:25]}
        """
        if isinstance(other, FeatureLocation):
            return CompoundLocation(self.parts + [other], self.operator)
        elif isinstance(other, CompoundLocation):
            if self.operator != other.operator:
                # Handle join+order -> order as a special case?
                raise ValueError("Mixed operators %s and %s"
                                 % (self.operator, other.operator))
            return CompoundLocation(self.parts + other.parts, self.operator)
        elif isinstance(other, int):
            return self._shift(other)
        else:
            raise NotImplementedError

    def __radd__(self, other):
        """Combine locations."""
        if isinstance(other, FeatureLocation):
            return CompoundLocation([other] + self.parts, self.operator)
        elif isinstance(other, int):
            return self._shift(other)
        else:
            raise NotImplementedError

    def __contains__(self, value):
        """Check if an integer position is within the location."""
        for loc in self.parts:
            if value in loc:
                return True
        return False

    def __nonzero__(self):
        """Returns True regardless of the length of the feature.

        This behaviour is for backwards compatibility, since until the
        __len__ method was added, a FeatureLocation always evaluated as True.

        Note that in comparison, Seq objects, strings, lists, etc, will all
        evaluate to False if they have length zero.

        WARNING: The FeatureLocation may in future evaluate to False when its
        length is zero (in order to better match normal python behaviour)!
        """
        return True

    def __len__(self):
        return sum(len(loc) for loc in self.parts)

    def __iter__(self):
        for loc in self.parts:
            for pos in loc:
                yield pos

    def _shift(self, offset):
        """Returns a copy of the location shifted by the offset (PRIVATE)."""
        return CompoundLocation([loc._shift(offset) for loc in self.parts],
                                self.operator)

    def _flip(self, length):
        """Returns a copy of the location after the parent is reversed (PRIVATE).

        Note that the order of the parts is reversed too.
        """
        return CompoundLocation([loc._flip(length) for loc in self.parts[::-1]],
                                self.operator)

    @property
    def start(self):
        """Start location (integer like, possibly a fuzzy position, read only)."""
        return min(loc.start for loc in self.parts)

    @property
    def end(self):
        """End location (integer like, possibly a fuzzy position, read only)."""
        return max(loc.end for loc in self.parts)

    @property
    def nofuzzy_start(self):
        """Start position (integer, approximated if fuzzy, read only) (OBSOLETE).

        This is an alias for int(feature.start), which should be used in
        preference -- unless you are trying to support old versions of
        Biopython.
        """
        try:
            return int(self.start)
        except TypeError:
            if isinstance(self.start, UnknownPosition):
                return None
            raise

    @property
    def nofuzzy_end(self):
        """End position (integer, approximated if fuzzy, read only) (OBSOLETE).

        This is an alias for int(feature.end), which should be used in
        preference -- unless you are trying to support old versions of
        Biopython.
        """
        try:
            return int(self.end)
        except TypeError:
            if isinstance(self.end, UnknownPosition):
                return None
            raise

    @property
    def ref(self):
        """CompoundLocation's don't have a ref (dummy method for API compatibility)."""
        return None

    @property
    def ref_db(self):
        """CompoundLocation's don't have a ref_db (dummy method for API compatibility)."""
        return None

    def extract(self, parent_sequence):
        """Extract feature sequence from the supplied parent sequence."""
        # This copes with mixed strand features & all on reverse:
        parts = [loc.extract(parent_sequence) for loc in self.parts]
        # We use addition rather than a join to avoid alphabet issues:
        f_seq = parts[0]
        for part in parts[1:]:
            f_seq += part
        return f_seq


class AbstractPosition(object):
    """Abstract base class representing a position."""

    def __repr__(self):
        """String representation of the location for debugging."""
        return "%s(...)" % (self.__class__.__name__)


class ExactPosition(int, AbstractPosition):
    """Specify the specific position of a boundary.

    o position - The position of the boundary.
    o extension - An optional argument which must be zero since we don't
    have an extension. The argument is provided so that the same number of
    arguments can be passed to all position types.

    In this case, there is no fuzziness associated with the position.

    >>> p = ExactPosition(5)
    >>> p
    ExactPosition(5)
    >>> print(p)
    5

    >>> isinstance(p, AbstractPosition)
    True
    >>> isinstance(p, int)
    True

    Integer comparisons and operations should work as expected:

    >>> p == 5
    True
    >>> p < 6
    True
    >>> p <= 5
    True
    >>> p + 10
    15

    """
    def __new__(cls, position, extension=0):
        if extension != 0:
            raise AttributeError("Non-zero extension %s for exact position."
                                 % extension)
        return int.__new__(cls, position)

    def __repr__(self):
        """String representation of the ExactPosition location for debugging."""
        return "%s(%i)" % (self.__class__.__name__, int(self))

    @property
    def position(self):
        """Legacy attribute to get position as integer (OBSOLETE)."""
        return int(self)

    @property
    def extension(self):
        """Legacy attribute to get extension (zero) as integer (OBSOLETE)."""
        return 0

    def _shift(self, offset):
        # By default preserve any subclass
        return self.__class__(int(self) + offset)

    def _flip(self, length):
        # By default perserve any subclass
        return self.__class__(length - int(self))


class UncertainPosition(ExactPosition):
    """Specify a specific position which is uncertain.

    This is used in UniProt, e.g. ?222 for uncertain position 222, or in the
    XML format explicitly marked as uncertain. Does not apply to GenBank/EMBL.
    """
    pass


class UnknownPosition(AbstractPosition):
    """Specify a specific position which is unknown (has no position).

    This is used in UniProt, e.g. ? or in the XML as unknown.
    """

    def __repr__(self):
        """String representation of the UnknownPosition location for debugging."""
        return "%s()" % self.__class__.__name__

    def __hash__(self):
        return hash(None)

    @property
    def position(self):
        """Legacy attribute to get position (None) (OBSOLETE)."""
        return None

    @property
    def extension(self):
        """Legacy attribute to get extension (zero) as integer (OBSOLETE)."""
        return 0

    def _shift(self, offset):
        return self

    def _flip(self, length):
        return self


class WithinPosition(int, AbstractPosition):
    """Specify the position of a boundary within some coordinates.

    Arguments:
    o position - The default integer position
    o left - The start (left) position of the boundary
    o right - The end (right) position of the boundary

    This allows dealing with a position like ((1.4)..100). This
    indicates that the start of the sequence is somewhere between 1
    and 4. Since this is a start coordinate, it should acts like
    it is at position 1 (or in Python counting, 0).

    >>> p = WithinPosition(10, 10, 13)
    >>> p
    WithinPosition(10, left=10, right=13)
    >>> print(p)
    (10.13)
    >>> int(p)
    10

    Basic integer comparisons and operations should work as though
    this were a plain integer:

    >>> p == 10
    True
    >>> p in [9, 10, 11]
    True
    >>> p < 11
    True
    >>> p + 10
    20

    >>> isinstance(p, WithinPosition)
    True
    >>> isinstance(p, AbstractPosition)
    True
    >>> isinstance(p, int)
    True

    Note this also applies for comparison to other position objects,
    where again the integer behaviour is used:

    >>> p == 10
    True
    >>> p == ExactPosition(10)
    True
    >>> p == BeforePosition(10)
    True
    >>> p == AfterPosition(10)
    True

    If this were an end point, you would want the position to be 13:

    >>> p2 = WithinPosition(13, 10, 13)
    >>> p2
    WithinPosition(13, left=10, right=13)
    >>> print(p2)
    (10.13)
    >>> int(p2)
    13
    >>> p2 == 13
    True
    >>> p2 == ExactPosition(13)
    True

    The old legacy properties of position and extension give the
    starting/lower/left position as an integer, and the distance
    to the ending/higher/right position as an integer. Note that
    the position object will act like either the left or the right
    end-point depending on how it was created:

    >>> p.position == p2.position == 10
    True
    >>> p.extension == p2.extension == 3
    True
    >>> int(p) == int(p2)
    False
    >>> p == 10
    True
    >>> p2 == 13
    True

    """
    def __new__(cls, position, left, right):
        assert position == left or position == right, \
            "WithinPosition: %r should match left %r or right %r" \
            (position, left, right)
        obj = int.__new__(cls, position)
        obj._left = left
        obj._right = right
        return obj

    def __repr__(self):
        """String representation of the WithinPosition location for debugging."""
        return "%s(%i, left=%i, right=%i)" \
               % (self.__class__.__name__, int(self),
                  self._left, self._right)

    def __str__(self):
        return "(%s.%s)" % (self._left, self._right)

    @property
    def position(self):
        """Legacy attribute to get (left) position as integer (OBSOLETE)."""
        return self._left

    @property
    def extension(self):
        """Legacy attribute to get extension (from left to right) as an integer (OBSOLETE)."""
        return self._right - self._left

    def _shift(self, offset):
        return self.__class__(int(self) + offset,
                              self._left + offset,
                              self._right + offset)

    def _flip(self, length):
        return self.__class__(length - int(self),
                              length - self._right,
                              length - self._left)


class BetweenPosition(int, AbstractPosition):
    """Specify the position of a boundary between two coordinates (OBSOLETE?).

    Arguments:
    o position - The default integer position
    o left - The start (left) position of the boundary
    o right - The end (right) position of the boundary

    This allows dealing with a position like 123^456. This
    indicates that the start of the sequence is somewhere between
    123 and 456. It is up to the parser to set the position argument
    to either boundary point (depending on if this is being used as
    a start or end of the feature). For example as a feature end:

    >>> p = BetweenPosition(456, 123, 456)
    >>> p
    BetweenPosition(456, left=123, right=456)
    >>> print(p)
    (123^456)
    >>> int(p)
    456

    Integer equality and comparison use the given position,

    >>> p == 456
    True
    >>> p in [455, 456, 457]
    True
    >>> p > 300
    True

    The old legacy properties of position and extension give the
    starting/lower/left position as an integer, and the distance
    to the ending/higher/right position as an integer. Note that
    the position object will act like either the left or the right
    end-point depending on how it was created:

    >>> p2 = BetweenPosition(123, left=123, right=456)
    >>> p.position == p2.position == 123
    True
    >>> p.extension
    333
    >>> p2.extension
    333
    >>> p.extension == p2.extension == 333
    True
    >>> int(p) == int(p2)
    False
    >>> p == 456
    True
    >>> p2 == 123
    True

    Note this potentially surprising behaviour:

    >>> BetweenPosition(123, left=123, right=456) == ExactPosition(123)
    True
    >>> BetweenPosition(123, left=123, right=456) == BeforePosition(123)
    True
    >>> BetweenPosition(123, left=123, right=456) == AfterPosition(123)
    True

    i.e. For equality (and sorting) the position objects behave like
    integers.
    """
    def __new__(cls, position, left, right):
        assert position == left or position == right
        obj = int.__new__(cls, position)
        obj._left = left
        obj._right = right
        return obj

    def __repr__(self):
        """String representation of the WithinPosition location for debugging."""
        return "%s(%i, left=%i, right=%i)" \
               % (self.__class__.__name__, int(self),
                  self._left, self._right)

    def __str__(self):
        return "(%s^%s)" % (self._left, self._right)

    @property
    def position(self):
        """Legacy attribute to get (left) position as integer (OBSOLETE)."""
        return self._left

    @property
    def extension(self):
        """Legacy attribute to get extension (from left to right) as an integer (OBSOLETE)."""
        return self._right - self._left

    def _shift(self, offset):
        return self.__class__(int(self) + offset,
                              self._left + offset,
                              self._right + offset)

    def _flip(self, length):
        return self.__class__(length - int(self),
                              length - self._right,
                              length - self._left)


class BeforePosition(int, AbstractPosition):
    """Specify a position where the actual location occurs before it.

    Arguments:
    o position - The upper boundary of where the location can occur.
    o extension - An optional argument which must be zero since we don't
    have an extension. The argument is provided so that the same number of
    arguments can be passed to all position types.

    This is used to specify positions like (<10..100) where the location
    occurs somewhere before position 10.

    >>> p = BeforePosition(5)
    >>> p
    BeforePosition(5)
    >>> print(p)
    <5
    >>> int(p)
    5
    >>> p + 10
    15

    Note this potentially surprising behaviour:

    >>> p == ExactPosition(5)
    True
    >>> p == AfterPosition(5)
    True

    Just remember that for equality and sorting the position objects act
    like integers.
    """
    # Subclasses int so can't use __init__
    def __new__(cls, position, extension=0):
        if extension != 0:
            raise AttributeError("Non-zero extension %s for exact position."
                                 % extension)
        return int.__new__(cls, position)

    @property
    def position(self):
        """Legacy attribute to get position as integer (OBSOLETE)."""
        return int(self)

    @property
    def extension(self):
        """Legacy attribute to get extension (zero) as integer (OBSOLETE)."""
        return 0

    def __repr__(self):
        """A string representation of the location for debugging."""
        return "%s(%i)" % (self.__class__.__name__, int(self))

    def __str__(self):
        return "<%s" % self.position

    def _shift(self, offset):
        return self.__class__(int(self) + offset)

    def _flip(self, length):
        return AfterPosition(length - int(self))


class AfterPosition(int, AbstractPosition):
    """Specify a position where the actual location is found after it.

    Arguments:
    o position - The lower boundary of where the location can occur.
    o extension - An optional argument which must be zero since we don't
    have an extension. The argument is provided so that the same number of
    arguments can be passed to all position types.

    This is used to specify positions like (>10..100) where the location
    occurs somewhere after position 10.

    >>> p = AfterPosition(7)
    >>> p
    AfterPosition(7)
    >>> print(p)
    >7
    >>> int(p)
    7
    >>> p + 10
    17

    >>> isinstance(p, AfterPosition)
    True
    >>> isinstance(p, AbstractPosition)
    True
    >>> isinstance(p, int)
    True

    Note this potentially surprising behaviour:

    >>> p == ExactPosition(7)
    True
    >>> p == BeforePosition(7)
    True

    Just remember that for equality and sorting the position objects act
    like integers.
    """
    # Subclasses int so can't use __init__
    def __new__(cls, position, extension=0):
        if extension != 0:
            raise AttributeError("Non-zero extension %s for exact position."
                                 % extension)
        return int.__new__(cls, position)

    @property
    def position(self):
        """Legacy attribute to get position as integer (OBSOLETE)."""
        return int(self)

    @property
    def extension(self):
        """Legacy attribute to get extension (zero) as integer (OBSOLETE)."""
        return 0

    def __repr__(self):
        """A string representation of the location for debugging."""
        return "%s(%i)" % (self.__class__.__name__, int(self))

    def __str__(self):
        return ">%s" % self.position

    def _shift(self, offset):
        return self.__class__(int(self) + offset)

    def _flip(self, length):
        return BeforePosition(length - int(self))


class OneOfPosition(int, AbstractPosition):
    """Specify a position where the location can be multiple positions.

    This models the GenBank 'one-of(1888,1901)' function, and tries
    to make this fit within the Biopython Position models. If this was
    a start position it should act like 1888, but as an end position 1901.

    >>> p = OneOfPosition(1888, [ExactPosition(1888), ExactPosition(1901)])
    >>> p
    OneOfPosition(1888, choices=[ExactPosition(1888), ExactPosition(1901)])
    >>> int(p)
    1888

    Interget comparisons and operators act like using int(p),

    >>> p == 1888
    True
    >>> p <= 1888
    True
    >>> p > 1888
    False
    >>> p + 100
    1988

    >>> isinstance(p, OneOfPosition)
    True
    >>> isinstance(p, AbstractPosition)
    True
    >>> isinstance(p, int)
    True

    The old legacy properties of position and extension give the
    starting/lowest/left-most position as an integer, and the
    distance to the ending/highest/right-most position as an integer.
    Note that the position object will act like one of the list of
    possible locations depending on how it was created:

    >>> p2 = OneOfPosition(1901, [ExactPosition(1888), ExactPosition(1901)])
    >>> p.position == p2.position == 1888
    True
    >>> p.extension == p2.extension == 13
    True
    >>> int(p) == int(p2)
    False
    >>> p == 1888
    True
    >>> p2 == 1901
    True

    """
    def __new__(cls, position, choices):
        """Initialize with a set of posssible positions.

        position_list is a list of AbstractPosition derived objects,
        specifying possible locations.

        position is an integer specifying the default behaviour.
        """
        assert position in choices, \
            "OneOfPosition: %r should match one of %r" % (position, choices)
        obj = int.__new__(cls, position)
        obj.position_choices = choices
        return obj

    @property
    def position(self):
        """Legacy attribute to get (left) position as integer (OBSOLETE)."""
        return min(int(pos) for pos in self.position_choices)

    @property
    def extension(self):
        """Legacy attribute to get extension as integer (OBSOLETE)."""
        positions = [int(pos) for pos in self.position_choices]
        return max(positions) - min(positions)

    def __repr__(self):
        """String representation of the OneOfPosition location for debugging."""
        return "%s(%i, choices=%r)" % (self.__class__.__name__,
                                       int(self), self.position_choices)

    def __str__(self):
        out = "one-of("
        for position in self.position_choices:
            out += "%s," % position
        # replace the last comma with the closing parenthesis
        out = out[:-1] + ")"
        return out

    def _shift(self, offset):
        return self.__class__(int(self) + offset,
                              [p._shift(offset) for p in self.position_choices])

    def _flip(self, length):
        return self.__class__(length - int(self),
                              [p._flip(length) for p in self.position_choices[::-1]])


class PositionGap(object):
    """Simple class to hold information about a gap between positions."""

    def __init__(self, gap_size):
        """Intialize with a position object containing the gap information.
        """
        self.gap_size = gap_size

    def __repr__(self):
        """A string representation of the position gap for debugging."""
        return "%s(%s)" % (self.__class__.__name__, repr(self.gap_size))

    def __str__(self):
        out = "gap(%s)" % self.gap_size
        return out


if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()
