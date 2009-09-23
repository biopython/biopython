"""Represent a Sequence Feature holding info about a part of a sequence.

This is heavily modeled after the Biocorba SeqFeature objects, and
may be pretty biased towards GenBank stuff since I'm writing it
for the GenBank parser output...

What's here:

Base class to hold a Feature.
----------------------------
classes:
o SeqFeature

Hold information about a Reference.
----------------------------------

This is an attempt to create a General class to hold Reference type
information.

classes:
o Reference

Specify locations of a feature on a Sequence.
---------------------------------------------

This aims to handle, in Ewan's words, 'the dreaded fuzziness issue' in
much the same way as Biocorba. This has the advantages of allowing us
to handle fuzzy stuff in case anyone needs it, and also be compatible
with Biocorba.

classes:
o FeatureLocation - Specify the start and end location of a feature.

o ExactPosition - Specify the position as being exact.
o WithinPosition - Specify a position occuring within some range.
o BetweenPosition - Specify a position occuring between a range.
o BeforePosition - Specify the position as being found before some base.
o AfterPosition - Specify the position as being found after some base.
"""

class SeqFeature(object):
    """Represent a Sequence Feature on an object.

    Attributes:
    o location - the location of the feature on the sequence
    o type - the specified type of the feature (ie. CDS, exon, repeat...)
    o location_operator - a string specifying how this SeqFeature may
    be related to others. For example, in the example GenBank feature
    shown below, the location_operator would be "join"
    o strand - A value specifying on which strand (of a DNA sequence, for
    instance) the feature deals with. 1 indicates the plus strand, -1 
    indicates the minus strand, 0 indicates both strands, and None indicates
    that strand doesn't apply (ie. for proteins) or is not known.
    o id - A string identifier for the feature.
    o ref - A reference to another sequence. This could be an accession
    number for some different sequence.
    o ref_db - A different database for the reference accession number.
    o qualifiers - A dictionary of qualifiers on the feature. These are
    analagous to the qualifiers from a GenBank feature table. The keys of
    the dictionary are qualifier names, the values are the qualifier
    values.
    o sub_features - Additional SeqFeatures which fall under this 'parent'
    feature. For instance, if we having something like:

    CDS    join(1..10,30..40,50..60)

    The the top level feature would be a CDS from 1 to 60, and the sub
    features would be of 'CDS_join' type and would be from 1 to 10, 30 to
    40 and 50 to 60, respectively.
    """
    def __init__(self, location = None, type = '', location_operator = '',
                 strand = None, id = "<unknown id>", 
                 qualifiers = None, sub_features = None,
                 ref = None, ref_db = None):
        """Initialize a SeqFeature on a Sequence.
        """
        self.location = location

        self.type = type
        self.location_operator = location_operator
        self.strand = strand
        self.id = id
        if qualifiers is None:
            qualifiers = {}
        self.qualifiers = qualifiers
        if sub_features is None:
            sub_features = []
        self.sub_features = sub_features
        self.ref = ref 
        self.ref_db = ref_db

    def __repr__(self):
        """A string representation of the record for debugging."""
        answer = "%s(%s" % (self.__class__.__name__, repr(self.location))
        if self.type :
            answer += ", type=%s" % repr(self.type)
        if self.location_operator :
            answer += ", location_operator=%s" % repr(self.location_operator)
        if self.strand :
            answer += ", strand=%s" % repr(self.strand)
        if self.id and self.id != "<unknown id>" :
            answer += ", id=%s" % repr(self.id)
        if self.ref :
            answer += ", ref=%s" % repr(self.ref)
        if self.ref_db :
            answer += ", ref_db=%s" % repr(self.ref_db)
        answer += ")"
        return answer

    def __str__(self):
        """A readable summary of the feature intended to be printed to screen.
        """
        out = "type: %s\n" % self.type
        out += "location: %s\n" % self.location
        out += "ref: %s:%s\n" % (self.ref, self.ref_db)
        out += "strand: %s\n" % self.strand
        out += "qualifiers: \n"
        qualifier_keys = self.qualifiers.keys()
        qualifier_keys.sort()
        for qual_key in qualifier_keys:
            out += "    Key: %s, Value: %s\n" % (qual_key,
                                               self.qualifiers[qual_key])
        if len(self.sub_features) != 0:
            out += "Sub-Features\n"
            for sub_feature in self.sub_features:
                out +="%s\n" % sub_feature

        return out

    def _shift(self, offset) :
        """Returns a copy of the feature with its location shifted (PRIVATE).

        The annotation qaulifiers are copied."""
        answer = SeqFeature(location = self.location._shift(offset),
                   type = self.type,
                   location_operator = self.location_operator,
                   strand = self.strand,
                   id = self.id,
                   #qualifiers = dict(self.qualifiers.iteritems()),
                   #sub_features = [f._shift(offset) for f in self.sub_features],
                   ref = self.ref,
                   ref_db = self.ref_db)
        #TODO - Sort out the use of sub_feature and qualifiers in __init___
        answer.sub_features = [f._shift(offset) for f in self.sub_features]
        answer.qualifiers = dict(self.qualifiers.iteritems())
        return answer

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
        #TODO - Update this is __init__ later accpets values
        return "%s(title=%s, ...)" % (self.__class__.__name__,
                                      repr(self.title))

# --- Handling feature locations

class FeatureLocation(object):
    """Specify the location of a feature along a sequence.

    This attempts to deal with fuzziness of position ends, but also
    make it easy to get the start and end in the 'normal' case (no
    fuzziness).

    You should access the start and end attributes with
    your_location.start and your_location.end. If the start and
    end are exact, this will return the positions, if not, we'll return
    the approriate Fuzzy class with info about the position and fuzziness.

    Note that the start and end location numbering follow Python's scheme,
    thus a GenBank entry of 123..150 (one based counting) becomes a location
    of [122:150] (zero based counting).
    """
    def __init__(self, start, end):
        """Specify the start and end of a sequence feature.

        start and end arguments specify the values where the feature begins
        and ends. These can either by any of the *Position objects that
        inherit from AbstractPosition, or can just be integers specifying the
        position. In the case of integers, the values are assumed to be
        exact and are converted in ExactPosition arguments. This is meant
        to make it easy to deal with non-fuzzy ends.
        """
        if isinstance(start, AbstractPosition):
            self._start = start
        else:
            self._start = ExactPosition(start)

        if isinstance(end, AbstractPosition):
            self._end = end
        else:
            self._end = ExactPosition(end)

    def __str__(self):
        """Returns a representation of the location (with python counting).

        For the simple case this uses the python splicing syntax, [122:150]
        (zero based counting) which GenBank would call 123..150 (one based
        counting).
        """
        return "[%s:%s]" % (self._start, self._end)

    def __repr__(self):
        """A string representation of the location for debugging."""
        return "%s(%s,%s)" \
               % (self.__class__.__name__, repr(self.start), repr(self.end))

    def _shift(self, offset) :
        """Returns a copy of the location shifted by the offset (PRIVATE)."""
        return FeatureLocation(start = self._start._shift(offset),
                               end = self._end._shift(offset))

    start = property(fget= lambda self : self._start,
                 doc="Start location (possibly a fuzzy position, read only).")

    end = property(fget= lambda self : self._end,
                   doc="End location (possibly a fuzzy position, read only).")

    def _get_nofuzzy_start(self) :
        #TODO - Do we still use the BetweenPosition class?
        if ((self._start == self._end) and isinstance(self._start,
             BetweenPosition)):
            return self._start.position
        else:
            return min(self._start.position,
                       self._start.position + self._start.extension)
    nofuzzy_start = property(fget=_get_nofuzzy_start,
        doc="""Start position (integer, approximated if fuzzy, read only).

        To get non-fuzzy attributes (ie. the position only) ask for
        'location.nofuzzy_start', 'location.nofuzzy_end'. These should return
        the largest range of the fuzzy position. So something like:
        (10.20)..(30.40) should return 10 for start, and 40 for end.
        """)

    def _get_nofuzzy_end(self) :
        #TODO - Do we still use the BetweenPosition class?
        if ((self._start == self._end) and isinstance(self._start,
             BetweenPosition)):
            return self._end.position
        else:
            return max(self._end.position,
                       self._end.position + self._end.extension)
    nofuzzy_end = property(fget=_get_nofuzzy_end,
        doc="""End position (integer, approximated if fuzzy, read only).

        To get non-fuzzy attributes (ie. the position only) ask for
        'location.nofuzzy_start', 'location.nofuzzy_end'. These should return
        the largest range of the fuzzy position. So something like:
        (10.20)..(30.40) should return 10 for start, and 40 for end.
        """)

class AbstractPosition(object):
    """Abstract base class representing a position.
    """
    def __init__(self, position, extension):
        self.position = position
        self.extension = extension

    def __repr__(self) :
        """String representation of the location for debugging."""
        return "%s(%s,%s)" % (self.__class__.__name__, \
                              repr(self.position), repr(self.extension))

    def __cmp__(self, other):
        """A simple comparison function for positions.

        This is very simple-minded and just compares the position attribute
        of the features; extensions are not considered at all. This could
        potentially be expanded to try to take advantage of extensions.
        """
        assert isinstance(other, AbstractPosition), \
          "We can only do comparisons between Biopython Position objects."

        return cmp(self.position, other.position)

    def _shift(self, offset) :
        #We want this to maintain the subclass when called from a subclass
        return self.__class__(self.position + offset, self.extension)
            
class ExactPosition(AbstractPosition):
    """Specify the specific position of a boundary.

    o position - The position of the boundary.
    o extension - An optional argument which must be zero since we don't
    have an extension. The argument is provided so that the same number of
    arguments can be passed to all position types.

    In this case, there is no fuzziness associated with the position.
    """
    def __init__(self, position, extension = 0):
        if extension != 0:
            raise AttributeError("Non-zero extension %s for exact position."
                                 % extension)
        AbstractPosition.__init__(self, position, 0)

    def __repr__(self) :
        """String representation of the ExactPosition location for debugging."""
        assert self.extension == 0
        return "%s(%s)" % (self.__class__.__name__, repr(self.position))

    def __str__(self):
        return str(self.position)

class WithinPosition(AbstractPosition):
    """Specify the position of a boundary within some coordinates.

    Arguments:
    o position - The start position of the boundary
    o extension - The range to which the boundary can extend.

    This allows dealing with a position like ((1.4)..100). This
    indicates that the start of the sequence is somewhere between 1
    and 4. To represent that with this class we would set position as
    1 and extension as 3.
    """
    def __init__(self, position, extension = 0):
        AbstractPosition.__init__(self, position, extension)

    def __str__(self):
        return "(%s.%s)" % (self.position, self.position + self.extension)

class BetweenPosition(AbstractPosition):
    """Specify the position of a boundary between two coordinates.

    Arguments:
    o position - The start position of the boundary.
    o extension - The range to the other position of a boundary.

    This specifies a coordinate which is found between the two positions.
    So this allows us to deal with a position like ((1^2)..100). To
    represent that with this class we set position as 1 and the
    extension as 1.
    """
    def __init__(self, position, extension = 0):
        AbstractPosition.__init__(self, position, extension)

    def __str__(self):
        return "(%s^%s)" % (self.position, self.position + self.extension)

class BeforePosition(AbstractPosition):
    """Specify a position where the actual location occurs before it.

    Arguments:
    o position - The upper boundary of where the location can occur.
    o extension - An optional argument which must be zero since we don't
    have an extension. The argument is provided so that the same number of
    arguments can be passed to all position types.

    This is used to specify positions like (<10..100) where the location
    occurs somewhere before position 10.
    """
    def __init__(self, position, extension = 0):
        if extension != 0:
            raise AttributeError("Non-zero extension %s for exact position."
                                 % extension)
        AbstractPosition.__init__(self, position, 0)

    def __repr__(self) :
        """A string representation of the location for debugging."""
        assert self.extension == 0
        return "%s(%s)" % (self.__class__.__name__, repr(self.position))

    def __str__(self):
        return "<%s" % self.position

class AfterPosition(AbstractPosition):
    """Specify a position where the actual location is found after it.

    Arguments:
    o position - The lower boundary of where the location can occur.
    o extension - An optional argument which must be zero since we don't
    have an extension. The argument is provided so that the same number of
    arguments can be passed to all position types.

    This is used to specify positions like (>10..100) where the location
    occurs somewhere after position 10.
    """
    def __init__(self, position, extension = 0):
        if extension != 0:
            raise AttributeError("Non-zero extension %s for exact position."
                                 % extension)
        AbstractPosition.__init__(self, position, 0)

    def __repr__(self) :
        """A string representation of the location for debugging."""
        assert self.extension == 0
        return "%s(%s)" % (self.__class__.__name__, repr(self.position))

    def __str__(self):
        return ">%s" % self.position

class OneOfPosition(AbstractPosition):
    """Specify a position where the location can be multiple positions.

    This models the GenBank 'one-of(1888,1901)' function, and tries
    to make this fit within the Biopython Position models. In our case
    the position of the "one-of" is set as the lowest choice, and the
    extension is the range to the highest choice.
    """
    def __init__(self, position_list):
        """Initialize with a set of posssible positions.

        position_list is a list of AbstractPosition derived objects,
        specifying possible locations.
        """
        # unique attribute for this type of positions
        self.position_choices = position_list
        # find the smallest and largest position in the choices
        smallest = None
        largest = None
        for position_choice in self.position_choices:
            assert isinstance(position_choice, AbstractPosition), \
              "Expected position objects, got %r" % position_choice
            if smallest is None and largest is None:
                smallest = position_choice.position
                largest = position_choice.position
            elif position_choice.position > largest:
                largest = position_choice.position
            elif position_choice.position < smallest:
                smallest = position_choice.position
        # initialize with our definition of position and extension
        AbstractPosition.__init__(self, smallest, largest - smallest)

    def __repr__(self) :
        """String representation of the OneOfPosition location for debugging."""
        return "%s(%s)" % (self.__class__.__name__, \
                           repr(self.position_choices))

    def __str__(self):
        out = "one-of("
        for position in self.position_choices:
            out += "%s," % position
        # replace the last comma with the closing parenthesis
        out = out[:-1] + ")"
        return out

class PositionGap(object):
    """Simple class to hold information about a gap between positions.
    """
    def __init__(self, gap_size):
        """Intialize with a position object containing the gap information.
        """
        self.gap_size = gap_size

    def __repr__(self) :
        """A string representation of the position gap for debugging."""
        return "%s(%s)" % (self.__class__.__name__, repr(self.gap_size))
    
    def __str__(self):
        out = "gap(%s)" % self.gap_size
        return out
