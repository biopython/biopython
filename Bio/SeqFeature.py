# Copyright 2000-2003 Jeff Chang.
# Copyright 2001-2008 Brad Chapman.
# Copyright 2005-2024 by Peter Cock.
# Copyright 2006-2009 Michiel de Hoon.
# All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Represent a Sequence Feature holding info about a part of a sequence.

This is heavily modeled after the Biocorba SeqFeature objects, and
may be pretty biased towards GenBank stuff since I'm writing it
for the GenBank parser output...

What's here:

Base class to hold a Feature
----------------------------

Classes:
 - SeqFeature

Hold information about a Reference
----------------------------------

This is an attempt to create a General class to hold Reference type
information.

Classes:
 - Reference

Specify locations of a feature on a Sequence
--------------------------------------------

This aims to handle, in Ewan Birney's words, 'the dreaded fuzziness issue'.
This has the advantages of allowing us to handle fuzzy stuff in case anyone
needs it, and also be compatible with BioPerl etc and BioSQL.

Classes:
 - Location - abstract base class of SimpleLocation and CompoundLocation.
 - SimpleLocation - Specify the start and end location of a feature.
 - CompoundLocation - Collection of SimpleLocation objects (for joins etc).
 - Position - abstract base class of ExactPosition, WithinPosition,
   BetweenPosition, AfterPosition, OneOfPosition, UncertainPosition, and
   UnknownPosition.
 - ExactPosition - Specify the position as being exact.
 - WithinPosition - Specify a position occurring within some range.
 - BetweenPosition - Specify a position occurring between a range (OBSOLETE?).
 - BeforePosition - Specify the position as being found before some base.
 - AfterPosition - Specify the position as being found after some base.
 - OneOfPosition - Specify a position consisting of multiple alternative positions.
 - UncertainPosition - Specify a specific position which is uncertain.
 - UnknownPosition - Represents missing information like '?' in UniProt.


Exceptions:
 - LocationParserError - Exception indicating a failure to parse a location
   string.

"""
import functools
import re
import warnings
from abc import ABC, abstractmethod

from Bio import BiopythonDeprecationWarning
from Bio import BiopythonParserWarning
from Bio.Seq import MutableSeq
from Bio.Seq import reverse_complement
from Bio.Seq import Seq


# Regular expressions for location parsing

_reference = r"(?:[a-zA-Z][a-zA-Z0-9_\.\|]*[a-zA-Z0-9]?\:)"
_oneof_position = r"one\-of\(\d+[,\d+]+\)"

_oneof_location = rf"[<>]?(?:\d+|{_oneof_position})\.\.[<>]?(?:\d+|{_oneof_position})"

_any_location = rf"({_reference}?{_oneof_location}|complement\({_oneof_location}\)|[^,]+|complement\([^,]+\))"

_split = re.compile(_any_location).split

assert _split("123..145")[1::2] == ["123..145"]
assert _split("123..145,200..209")[1::2] == ["123..145", "200..209"]
assert _split("one-of(200,203)..300")[1::2] == ["one-of(200,203)..300"]
assert _split("complement(123..145),200..209")[1::2] == [
    "complement(123..145)",
    "200..209",
]
assert _split("123..145,one-of(200,203)..209")[1::2] == [
    "123..145",
    "one-of(200,203)..209",
]
assert _split("123..145,one-of(200,203)..one-of(209,211),300")[1::2] == [
    "123..145",
    "one-of(200,203)..one-of(209,211)",
    "300",
]
assert _split("123..145,complement(one-of(200,203)..one-of(209,211)),300")[1::2] == [
    "123..145",
    "complement(one-of(200,203)..one-of(209,211))",
    "300",
]
assert _split("123..145,200..one-of(209,211),300")[1::2] == [
    "123..145",
    "200..one-of(209,211)",
    "300",
]
assert _split("123..145,200..one-of(209,211)")[1::2] == [
    "123..145",
    "200..one-of(209,211)",
]
assert _split(
    "complement(149815..150200),complement(293787..295573),NC_016402.1:6618..6676,181647..181905"
)[1::2] == [
    "complement(149815..150200)",
    "complement(293787..295573)",
    "NC_016402.1:6618..6676",
    "181647..181905",
]


_pair_location = r"[<>]?-?\d+\.\.[<>]?-?\d+"

_between_location = r"\d+\^\d+"

_within_position = r"\(\d+\.\d+\)"
_within_location = r"([<>]?\d+|%s)\.\.([<>]?\d+|%s)" % (
    _within_position,
    _within_position,
)
_within_position = r"\((\d+)\.(\d+)\)"
_re_within_position = re.compile(_within_position)
assert _re_within_position.match("(3.9)")

_oneof_location = r"([<>]?\d+|%s)\.\.([<>]?\d+|%s)" % (_oneof_position, _oneof_position)
_oneof_position = r"one\-of\((\d+[,\d+]+)\)"
_re_oneof_position = re.compile(_oneof_position)
assert _re_oneof_position.match("one-of(6,9)")
assert not _re_oneof_position.match("one-of(3)")
assert _re_oneof_position.match("one-of(3,6)")
assert _re_oneof_position.match("one-of(3,6,9)")

_solo_location = r"[<>]?\d+"
_solo_bond = r"bond\(%s\)" % _solo_location

_re_location_category = re.compile(
    r"^(?P<pair>%s)|(?P<between>%s)|(?P<within>%s)|(?P<oneof>%s)|(?P<bond>%s)|(?P<solo>%s)$"
    % (
        _pair_location,
        _between_location,
        _within_location,
        _oneof_location,
        _solo_bond,
        _solo_location,
    )
)


class LocationParserError(ValueError):
    """Could not parse a feature location string."""


class SeqFeature:
    """Represent a Sequence Feature on an object.

    Attributes:
     - location - the location of the feature on the sequence (SimpleLocation)
     - type - the specified type of the feature (ie. CDS, exon, repeat...)
     - id - A string identifier for the feature.
     - qualifiers - A dictionary of qualifiers on the feature. These are
       analogous to the qualifiers from a GenBank feature table. The keys of
       the dictionary are qualifier names, the values are the qualifier
       values.

    """

    def __init__(
        self,
        location=None,
        type="",
        id="<unknown id>",
        qualifiers=None,
        sub_features=None,
    ):
        """Initialize a SeqFeature on a sequence.

        location can either be a SimpleLocation (with strand argument also
        given if required), or None.

        e.g. With no strand, on the forward strand, and on the reverse strand:

        >>> from Bio.SeqFeature import SeqFeature, SimpleLocation
        >>> f1 = SeqFeature(SimpleLocation(5, 10), type="domain")
        >>> f1.location.strand == None
        True
        >>> f2 = SeqFeature(SimpleLocation(7, 110, strand=1), type="CDS")
        >>> f2.location.strand == +1
        True
        >>> f3 = SeqFeature(SimpleLocation(9, 108, strand=-1), type="CDS")
        >>> f3.location.strand == -1
        True

        For exact start/end positions, an integer can be used (as shown above)
        as shorthand for the ExactPosition object. For non-exact locations, the
        SimpleLocation must be specified via the appropriate position objects.
        """
        if (
            location is not None
            and not isinstance(location, SimpleLocation)
            and not isinstance(location, CompoundLocation)
        ):
            raise TypeError(
                "SimpleLocation, CompoundLocation (or None) required for the location"
            )
        self.location = location
        self.type = type
        self.id = id
        self.qualifiers = {}
        if qualifiers is not None:
            self.qualifiers.update(qualifiers)
        if sub_features is not None:
            raise TypeError("Rather than sub_features, use a CompoundLocation")

    def _get_strand(self):
        """Get function for the strand property (PRIVATE)."""
        warnings.warn(
            "Please use .location.strand rather than .strand",
            BiopythonDeprecationWarning,
        )
        return self.location.strand

    def _set_strand(self, value):
        """Set function for the strand property (PRIVATE)."""
        warnings.warn(
            "Please use .location.strand rather than .strand",
            BiopythonDeprecationWarning,
        )
        try:
            self.location.strand = value
        except AttributeError:
            if self.location is None:
                if value is not None:
                    raise ValueError("Can't set strand without a location.") from None
            else:
                raise

    strand = property(
        fget=_get_strand,
        fset=_set_strand,
        doc="Alias for the location's strand (DEPRECATED).",
    )

    def _get_ref(self):
        """Get function for the reference property (PRIVATE)."""
        warnings.warn(
            "Please use .location.ref rather than .ref",
            BiopythonDeprecationWarning,
        )
        try:
            return self.location.ref
        except AttributeError:
            return None

    def _set_ref(self, value):
        """Set function for the reference property (PRIVATE)."""
        warnings.warn(
            "Please use .location.ref rather than .ref",
            BiopythonDeprecationWarning,
        )
        try:
            self.location.ref = value
        except AttributeError:
            if self.location is None:
                if value is not None:
                    raise ValueError("Can't set ref without a location.") from None
            else:
                raise

    ref = property(
        fget=_get_ref,
        fset=_set_ref,
        doc="Alias for the location's ref (DEPRECATED).",
    )

    def _get_ref_db(self):
        """Get function for the database reference property (PRIVATE)."""
        warnings.warn(
            "Please use .location.ref_db rather than .ref_db",
            BiopythonDeprecationWarning,
        )
        try:
            return self.location.ref_db
        except AttributeError:
            return None

    def _set_ref_db(self, value):
        """Set function for the database reference property (PRIVATE)."""
        warnings.warn(
            "Please use .location.ref_db rather than .ref_db",
            BiopythonDeprecationWarning,
        )
        self.location.ref_db = value

    ref_db = property(
        fget=_get_ref_db,
        fset=_set_ref_db,
        doc="Alias for the location's ref_db (DEPRECATED).",
    )

    def __eq__(self, other):
        """Check if two SeqFeature objects should be considered equal."""
        return (
            isinstance(other, SeqFeature)
            and self.id == other.id
            and self.type == other.type
            and self.location == other.location
            and self.qualifiers == other.qualifiers
        )

    def __repr__(self):
        """Represent the feature as a string for debugging."""
        answer = f"{self.__class__.__name__}({self.location!r}"
        if self.type:
            answer += f", type={self.type!r}"
        if self.id and self.id != "<unknown id>":
            answer += f", id={self.id!r}"
        if self.qualifiers:
            answer += ", qualifiers=..."
        answer += ")"
        return answer

    def __str__(self):
        """Return the full feature as a python string."""
        out = f"type: {self.type}\n"
        out += f"location: {self.location}\n"
        if self.id and self.id != "<unknown id>":
            out += f"id: {self.id}\n"
        out += "qualifiers:\n"
        for qual_key in sorted(self.qualifiers):
            out += f"    Key: {qual_key}, Value: {self.qualifiers[qual_key]}\n"
        return out

    def _shift(self, offset):
        """Return a copy of the feature with its location shifted (PRIVATE).

        The annotation qualifiers are copied.
        """
        return SeqFeature(
            location=self.location._shift(offset),
            type=self.type,
            id=self.id,
            qualifiers=self.qualifiers.copy(),
        )

    def _flip(self, length):
        """Return a copy of the feature with its location flipped (PRIVATE).

        The argument length gives the length of the parent sequence. For
        example a location 0..20 (+1 strand) with parent length 30 becomes
        after flipping 10..30 (-1 strand). Strandless (None) or unknown
        strand (0) remain like that - just their end points are changed.

        The annotation qualifiers are copied.
        """
        return SeqFeature(
            location=self.location._flip(length),
            type=self.type,
            id=self.id,
            qualifiers=self.qualifiers.copy(),
        )

    def extract(self, parent_sequence, references=None):
        """Extract the feature's sequence from supplied parent sequence.

        The parent_sequence can be a Seq like object or a string, and will
        generally return an object of the same type. The exception to this is
        a MutableSeq as the parent sequence will return a Seq object.

        This should cope with complex locations including complements, joins
        and fuzzy positions. Even mixed strand features should work! This
        also covers features on protein sequences (e.g. domains), although
        here reverse strand features are not permitted. If the
        location refers to other records, they must be supplied in the
        optional dictionary references.

        >>> from Bio.Seq import Seq
        >>> from Bio.SeqFeature import SeqFeature, SimpleLocation
        >>> seq = Seq("MKQHKAMIVALIVICITAVVAAL")
        >>> f = SeqFeature(SimpleLocation(8, 15), type="domain")
        >>> f.extract(seq)
        Seq('VALIVIC')

        If the SimpleLocation is None, e.g. when parsing invalid locus
        locations in the GenBank parser, extract() will raise a ValueError.

        >>> from Bio.Seq import Seq
        >>> from Bio.SeqFeature import SeqFeature
        >>> seq = Seq("MKQHKAMIVALIVICITAVVAAL")
        >>> f = SeqFeature(None, type="domain")
        >>> f.extract(seq)
        Traceback (most recent call last):
           ...
        ValueError: The feature's .location is None. Check the sequence file for a valid location.

        Note - currently only compound features of type "join" are supported.
        """
        if self.location is None:
            raise ValueError(
                "The feature's .location is None. Check the "
                "sequence file for a valid location."
            )
        return self.location.extract(parent_sequence, references=references)

    def translate(
        self,
        parent_sequence,
        table="Standard",
        start_offset=None,
        stop_symbol="*",
        to_stop=False,
        cds=None,
        gap=None,
    ):
        """Get a translation of the feature's sequence.

        This method is intended for CDS or other features that code proteins
        and is a shortcut that will both extract the feature and
        translate it, taking into account the codon_start and transl_table
        qualifiers, if they are present. If they are not present the
        value of the arguments "table" and "start_offset" are used.

        The "cds" parameter is set to "True" if the feature is of type
        "CDS" but can be overridden by giving an explicit argument.

        The arguments stop_symbol, to_stop and gap have the same meaning
        as Seq.translate, refer to that documentation for further information.

        Arguments:
         - parent_sequence - A DNA or RNA sequence.
         - table - Which codon table to use if there is no transl_table
           qualifier for this feature. This can be either a name
           (string), an NCBI identifier (integer), or a CodonTable
           object (useful for non-standard genetic codes).  This
           defaults to the "Standard" table.
         - start_offset - offset at which the first complete codon of a
           coding feature can be found, relative to the first base of
           that feature. Has a valid value of 0, 1 or 2. NOTE: this
           uses python's 0-based numbering whereas the codon_start
           qualifier in files from NCBI use 1-based numbering.
           Will override a codon_start qualifier

        >>> from Bio.Seq import Seq
        >>> from Bio.SeqFeature import SeqFeature, SimpleLocation
        >>> seq = Seq("GGTTACACTTACCGATAATGTCTCTGATGA")
        >>> f = SeqFeature(SimpleLocation(0, 30), type="CDS")
        >>> f.qualifiers['transl_table'] = [11]

        Note that features of type CDS are subject to the usual
        checks at translation. But you can override this behavior
        by giving explicit arguments:

        >>> f.translate(seq, cds=False)
        Seq('GYTYR*CL**')

        Now use the start_offset argument to change the frame. Note
        this uses python 0-based numbering.

        >>> f.translate(seq, start_offset=1, cds=False)
        Seq('VTLTDNVSD')

        Alternatively use the codon_start qualifier to do the same
        thing. Note: this uses 1-based numbering, which is found
        in files from NCBI.

        >>> f.qualifiers['codon_start'] = [2]
        >>> f.translate(seq, cds=False)
        Seq('VTLTDNVSD')
        """
        # see if this feature should be translated in a different
        # frame using the "codon_start" qualifier
        if start_offset is None:
            try:
                start_offset = int(self.qualifiers["codon_start"][0]) - 1
            except KeyError:
                start_offset = 0

        if start_offset not in [0, 1, 2]:
            raise ValueError(
                "The start_offset must be 0, 1, or 2. "
                f"The supplied value is {start_offset}. "
                "Check the value of either the codon_start qualifier "
                "or the start_offset argument"
            )

        feat_seq = self.extract(parent_sequence)[start_offset:]
        codon_table = self.qualifiers.get("transl_table", [table])[0]

        if cds is None:
            cds = self.type == "CDS"

        return feat_seq.translate(
            table=codon_table,
            stop_symbol=stop_symbol,
            to_stop=to_stop,
            cds=cds,
            gap=gap,
        )

    def __bool__(self):
        """Boolean value of an instance of this class (True).

        This behavior is for backwards compatibility, since until the
        __len__ method was added, a SeqFeature always evaluated as True.

        Note that in comparison, Seq objects, strings, lists, etc, will all
        evaluate to False if they have length zero.

        WARNING: The SeqFeature may in future evaluate to False when its
        length is zero (in order to better match normal python behavior)!
        """
        return True

    def __len__(self):
        """Return the length of the region where the feature is located.

        >>> from Bio.Seq import Seq
        >>> from Bio.SeqFeature import SeqFeature, SimpleLocation
        >>> seq = Seq("MKQHKAMIVALIVICITAVVAAL")
        >>> f = SeqFeature(SimpleLocation(8, 15), type="domain")
        >>> len(f)
        7
        >>> f.extract(seq)
        Seq('VALIVIC')
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

        >>> from Bio.SeqFeature import SeqFeature, SimpleLocation
        >>> f = SeqFeature(SimpleLocation(5, 10, strand=-1), type="domain")
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

        >>> from Bio.SeqFeature import SeqFeature, SimpleLocation
        >>> f = SeqFeature(SimpleLocation(5, 10, strand=-1), type="domain")
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

        >>> from Bio.SeqFeature import SeqFeature, SimpleLocation
        >>> from Bio.SeqFeature import BeforePosition
        >>> f = SeqFeature(SimpleLocation(BeforePosition(3), 8), type="domain")
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
class Reference:
    """Represent a Generic Reference object.

    Attributes:
     - location - A list of Location objects specifying regions of
       the sequence that the references correspond to. If no locations are
       specified, the entire sequence is assumed.
     - authors - A big old string, or a list split by author, of authors
       for the reference.
     - title - The title of the reference.
     - journal - Journal the reference was published in.
     - medline_id - A medline reference for the article.
     - pubmed_id - A pubmed reference for the article.
     - comment - A place to stick any comments about the reference.

    """

    def __init__(self):
        """Initialize the class."""
        self.location = []
        self.authors = ""
        self.consrtm = ""
        self.title = ""
        self.journal = ""
        self.medline_id = ""
        self.pubmed_id = ""
        self.comment = ""

    def __str__(self):
        """Return the full Reference object as a python string."""
        out = ""
        for single_location in self.location:
            out += f"location: {single_location}\n"
        out += f"authors: {self.authors}\n"
        if self.consrtm:
            out += f"consrtm: {self.consrtm}\n"
        out += f"title: {self.title}\n"
        out += f"journal: {self.journal}\n"
        out += f"medline id: {self.medline_id}\n"
        out += f"pubmed id: {self.pubmed_id}\n"
        out += f"comment: {self.comment}\n"
        return out

    def __repr__(self):
        """Represent the Reference object as a string for debugging."""
        # TODO - Update this is __init__ later accepts values
        return f"{self.__class__.__name__}(title={self.title!r}, ...)"

    def __eq__(self, other):
        """Check if two Reference objects should be considered equal.

        Note prior to Biopython 1.70 the location was not compared, as
        until then __eq__ for the SimpleLocation class was not defined.
        """
        return (
            self.authors == other.authors
            and self.consrtm == other.consrtm
            and self.title == other.title
            and self.journal == other.journal
            and self.medline_id == other.medline_id
            and self.pubmed_id == other.pubmed_id
            and self.comment == other.comment
            and self.location == other.location
        )


# --- Handling feature locations


class Location(ABC):
    """Abstract base class representing a location."""

    @abstractmethod
    def __repr__(self):
        """Represent the Location object as a string for debugging."""
        return f"{self.__class__.__name__}(...)"

    def fromstring(text, length=None, circular=False, stranded=True):
        """Create a Location object from a string.

        This should accept any valid location string in the INSDC Feature Table
        format (https://www.insdc.org/submitting-standards/feature-table/) as
        used in GenBank, DDBJ and EMBL files.

        Simple examples:

        >>> Location.fromstring("123..456", 1000)
        SimpleLocation(ExactPosition(122), ExactPosition(456), strand=1)
        >>> Location.fromstring("complement(<123..>456)", 1000)
        SimpleLocation(BeforePosition(122), AfterPosition(456), strand=-1)

        A more complex location using within positions,

        >>> Location.fromstring("(9.10)..(20.25)", 1000)
        SimpleLocation(WithinPosition(8, left=8, right=9), WithinPosition(25, left=20, right=25), strand=1)

        Notice how that will act as though it has overall start 8 and end 25.

        Zero length between feature,

        >>> Location.fromstring("123^124", 1000)
        SimpleLocation(ExactPosition(123), ExactPosition(123), strand=1)

        The expected sequence length is needed for a special case, a between
        position at the start/end of a circular genome:

        >>> Location.fromstring("1000^1", 1000)
        SimpleLocation(ExactPosition(1000), ExactPosition(1000), strand=1)

        Apart from this special case, between positions P^Q must have P+1==Q,

        >>> Location.fromstring("123^456", 1000)
        Traceback (most recent call last):
           ...
        Bio.SeqFeature.LocationParserError: invalid feature location '123^456'

        You can optionally provide a reference name:

        >>> Location.fromstring("AL391218.9:105173..108462", 2000000)
        SimpleLocation(ExactPosition(105172), ExactPosition(108462), strand=1, ref='AL391218.9')

        >>> Location.fromstring("<2644..159", 2868, "circular")
        CompoundLocation([SimpleLocation(BeforePosition(2643), ExactPosition(2868), strand=1), SimpleLocation(ExactPosition(0), ExactPosition(159), strand=1)], 'join')
        """
        if text.startswith("complement("):
            if text[-1] != ")":
                raise ValueError(f"closing bracket missing in '{text}'")
            text = text[11:-1]
            strand = -1
        elif stranded:
            strand = 1
        else:
            strand = None

        # Determine if we have a simple location or a compound location
        if text.startswith("join("):
            operator = "join"
            parts = _split(text[5:-1])[1::2]
            # assert parts[0] == "" and parts[-1] == ""
        elif text.startswith("order("):
            operator = "order"
            parts = _split(text[6:-1])[1::2]
            # assert parts[0] == "" and parts[-1] == ""
        elif text.startswith("bond("):
            operator = "bond"
            parts = _split(text[5:-1])[1::2]
            # assert parts[0] == "" and parts[-1] == ""
        else:
            loc = SimpleLocation.fromstring(text, length, circular)
            loc.strand = strand
            if strand == -1:
                loc.parts.reverse()
            return loc
        locs = []
        for part in parts:
            loc = SimpleLocation.fromstring(part, length, circular)
            if loc is None:
                break
            if loc.strand == -1:
                if strand == -1:
                    raise LocationParserError("double complement in '{text}'?")
            else:
                loc.strand = strand
            locs.extend(loc.parts)
        else:
            if len(locs) == 1:
                return loc
            # Historically a join on the reverse strand has been represented
            # in Biopython with both the parent SeqFeature and its children
            # (the exons for a CDS) all given a strand of -1.  Likewise, for
            # a join feature on the forward strand they all have strand +1.
            # However, we must also consider evil mixed strand examples like
            # this, join(complement(69611..69724),139856..140087,140625..140650)
            if strand == -1:
                # Whole thing was wrapped in complement(...)
                for loc in locs:
                    assert loc.strand == -1
                # Reverse the backwards order used in GenBank files
                # with complement(join(...))
                locs = locs[::-1]
            return CompoundLocation(locs, operator=operator)
        # Not recognized
        if "order" in text and "join" in text:
            # See Bug 3197
            raise LocationParserError(
                f"failed to parse feature location '{text}' containing a combination of 'join' and 'order' (nested operators) are illegal"
            )

        # See issue #937. Note that NCBI has already fixed this record.
        if ",)" in text:
            warnings.warn(
                "Dropping trailing comma in malformed feature location",
                BiopythonParserWarning,
            )
            text = text.replace(",)", ")")
            return Location.fromstring(text)

        raise LocationParserError(f"failed to parse feature location '{text}'")


class SimpleLocation(Location):
    """Specify the location of a feature along a sequence.

    The SimpleLocation is used for simple continuous features, which can
    be described as running from a start position to and end position
    (optionally with a strand and reference information).  More complex
    locations made up from several non-continuous parts (e.g. a coding
    sequence made up of several exons) are described using a SeqFeature
    with a CompoundLocation.

    Note that the start and end location numbering follow Python's scheme,
    thus a GenBank entry of 123..150 (one based counting) becomes a location
    of [122:150] (zero based counting).

    >>> from Bio.SeqFeature import SimpleLocation
    >>> f = SimpleLocation(122, 150)
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

    >>> from Bio.SeqFeature import SimpleLocation
    >>> f = SimpleLocation(122, 150, strand=+1)
    >>> print(f)
    [122:150](+)
    >>> print(f.strand)
    1

    Note that for a parent sequence of length n, the SimpleLocation
    start and end must satisfy the inequality 0 <= start <= end <= n.
    This means even for features on the reverse strand of a nucleotide
    sequence, we expect the 'start' coordinate to be less than the
    'end'.

    >>> from Bio.SeqFeature import SimpleLocation
    >>> r = SimpleLocation(122, 150, strand=-1)
    >>> print(r)
    [122:150](-)
    >>> print(r.start)
    122
    >>> print(r.end)
    150
    >>> print(r.strand)
    -1

    i.e. Rather than thinking of the 'start' and 'end' biologically in a
    strand aware manner, think of them as the 'left most' or 'minimum'
    boundary, and the 'right most' or 'maximum' boundary of the region
    being described. This is particularly important with compound
    locations describing non-continuous regions.

    In the example above we have used standard exact positions, but there
    are also specialised position objects used to represent fuzzy positions
    as well, for example a GenBank location like complement(<123..150)
    would use a BeforePosition object for the start.
    """

    def __init__(self, start, end, strand=None, ref=None, ref_db=None):
        """Initialize the class.

        start and end arguments specify the values where the feature begins
        and ends. These can either by any of the ``*Position`` objects that
        inherit from Position, or can just be integers specifying the position.
        In the case of integers, the values are assumed to be exact and are
        converted in ExactPosition arguments. This is meant to make it easy
        to deal with non-fuzzy ends.

        i.e. Short form:

        >>> from Bio.SeqFeature import SimpleLocation
        >>> loc = SimpleLocation(5, 10, strand=-1)
        >>> print(loc)
        [5:10](-)

        Explicit form:

        >>> from Bio.SeqFeature import SimpleLocation, ExactPosition
        >>> loc = SimpleLocation(ExactPosition(5), ExactPosition(10), strand=-1)
        >>> print(loc)
        [5:10](-)

        Other fuzzy positions are used similarly,

        >>> from Bio.SeqFeature import SimpleLocation
        >>> from Bio.SeqFeature import BeforePosition, AfterPosition
        >>> loc2 = SimpleLocation(BeforePosition(5), AfterPosition(10), strand=-1)
        >>> print(loc2)
        [<5:>10](-)

        For nucleotide features you will also want to specify the strand,
        use 1 for the forward (plus) strand, -1 for the reverse (negative)
        strand, 0 for stranded but strand unknown (? in GFF3), or None for
        when the strand does not apply (dot in GFF3), e.g. features on
        proteins.

        >>> loc = SimpleLocation(5, 10, strand=+1)
        >>> print(loc)
        [5:10](+)
        >>> print(loc.strand)
        1

        Normally feature locations are given relative to the parent
        sequence you are working with, but an explicit accession can
        be given with the optional ref and db_ref strings:

        >>> loc = SimpleLocation(105172, 108462, ref="AL391218.9", strand=1)
        >>> print(loc)
        AL391218.9[105172:108462](+)
        >>> print(loc.ref)
        AL391218.9

        """
        # TODO - Check 0 <= start <= end (<= length of reference)
        if isinstance(start, Position):
            self._start = start
        elif isinstance(start, int):
            self._start = ExactPosition(start)
        else:
            raise TypeError(f"start={start!r} {type(start)}")
        if isinstance(end, Position):
            self._end = end
        elif isinstance(end, int):
            self._end = ExactPosition(end)
        else:
            raise TypeError(f"end={end!r} {type(end)}")
        if (
            isinstance(self.start, int)
            and isinstance(self.end, int)
            and self.start > self.end
        ):
            raise ValueError(
                f"End location ({self.end}) must be greater than "
                f"or equal to start location ({self.start})"
            )
        self.strand = strand
        self.ref = ref
        self.ref_db = ref_db

    @staticmethod
    def fromstring(text, length=None, circular=False):
        """Create a SimpleLocation object from a string."""
        if text.startswith("complement("):
            text = text[11:-1]
            strand = -1
        else:
            strand = None
        # Try simple cases first for speed
        try:
            s, e = text.split("..")
            s = int(s) - 1
            e = int(e)
        except ValueError:
            pass
        else:
            if 0 <= s < e:
                return SimpleLocation(s, e, strand)
        # Try general case
        try:
            ref, text = text.split(":")
        except ValueError:
            ref = None
        m = _re_location_category.match(text)
        if m is None:
            raise LocationParserError(f"Could not parse feature location '{text}'")
        for key, value in m.groupdict().items():
            if value is not None:
                break
        assert value == text
        if key == "bond":
            # e.g. bond(196)
            warnings.warn(
                "Dropping bond qualifier in feature location",
                BiopythonParserWarning,
            )
            text = text[5:-1]
            s_pos = Position.fromstring(text, -1)
            e_pos = Position.fromstring(text)
        elif key == "solo":
            # e.g. "123"
            s_pos = Position.fromstring(text, -1)
            e_pos = Position.fromstring(text)
        elif key in ("pair", "within", "oneof"):
            s, e = text.split("..")
            # Attempt to fix features that span the origin
            s_pos = Position.fromstring(s, -1)
            e_pos = Position.fromstring(e)
            if s_pos >= e_pos:
                # There is likely a problem with origin wrapping.
                # Create a CompoundLocation of the wrapped feature,
                # consisting of two SimpleLocation objects to extend to
                # the list of feature locations.
                if not circular:
                    raise LocationParserError(
                        f"it appears that '{text}' is a feature that spans the origin, but the sequence topology is undefined"
                    )
                warnings.warn(
                    "Attempting to fix invalid location %r as "
                    "it looks like incorrect origin wrapping. "
                    "Please fix input file, this could have "
                    "unintended behavior." % text,
                    BiopythonParserWarning,
                )

                f1 = SimpleLocation(s_pos, length, strand)
                f2 = SimpleLocation(0, e_pos, strand)

                if strand == -1:
                    # For complementary features spanning the origin
                    return f2 + f1
                else:
                    return f1 + f2
        elif key == "between":
            # A between location like "67^68" (one based counting) is a
            # special case (note it has zero length). In python slice
            # notation this is 67:67, a zero length slice.  See Bug 2622
            # Further more, on a circular genome of length N you can have
            # a location N^1 meaning the junction at the origin. See Bug 3098.
            # NOTE - We can imagine between locations like "2^4", but this
            # is just "3".  Similarly, "2^5" is just "3..4"
            s, e = text.split("^")
            s = int(s)
            e = int(e)
            if s + 1 == e or (s == length and e == 1):
                s_pos = ExactPosition(s)
                e_pos = s_pos
            else:
                raise LocationParserError(f"invalid feature location '{text}'")
        if s_pos < 0:
            raise LocationParserError(
                f"negative starting position in feature location '{text}'"
            )
        return SimpleLocation(s_pos, e_pos, strand, ref=ref)

    def _get_strand(self):
        """Get function for the strand property (PRIVATE)."""
        return self._strand

    def _set_strand(self, value):
        """Set function for the strand property (PRIVATE)."""
        if value not in [+1, -1, 0, None]:
            raise ValueError(f"Strand should be +1, -1, 0 or None, not {value!r}")
        self._strand = value

    strand = property(
        fget=_get_strand,
        fset=_set_strand,
        doc="Strand of the location (+1, -1, 0 or None).",
    )

    def __str__(self):
        """Return a representation of the SimpleLocation object (with python counting).

        For the simple case this uses the python splicing syntax, [122:150]
        (zero based counting) which GenBank would call 123..150 (one based
        counting).
        """
        answer = f"[{self._start}:{self._end}]"
        if self.ref and self.ref_db:
            answer = f"{self.ref_db}:{self.ref}{answer}"
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
        """Represent the SimpleLocation object as a string for debugging."""
        optional = ""
        if self.strand is not None:
            optional += f", strand={self.strand!r}"
        if self.ref is not None:
            optional += f", ref={self.ref!r}"
        if self.ref_db is not None:
            optional += f", ref_db={self.ref_db!r}"
        return f"{self.__class__.__name__}({self.start!r}, {self.end!r}{optional})"

    def __add__(self, other):
        """Combine location with another SimpleLocation object, or shift it.

        You can add two feature locations to make a join CompoundLocation:

        >>> from Bio.SeqFeature import SimpleLocation
        >>> f1 = SimpleLocation(5, 10)
        >>> f2 = SimpleLocation(20, 30)
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

        Furthermore, you can combine a SimpleLocation with a CompoundLocation
        in this way.

        Separately, adding an integer will give a new SimpleLocation with
        its start and end offset by that amount. For example:

        >>> print(f1)
        [5:10]
        >>> print(f1 + 100)
        [105:110]
        >>> print(200 + f1)
        [205:210]

        This can be useful when editing annotation.
        """
        if isinstance(other, SimpleLocation):
            return CompoundLocation([self, other])
        elif isinstance(other, int):
            return self._shift(other)
        else:
            # This will allow CompoundLocation's __radd__ to be called:
            return NotImplemented

    def __radd__(self, other):
        """Return a SimpleLocation object by shifting the location by an integer amount."""
        if isinstance(other, int):
            return self._shift(other)
        else:
            return NotImplemented

    def __sub__(self, other):
        """Subtracting an integer will shift the start and end by that amount.

        >>> from Bio.SeqFeature import SimpleLocation
        >>> f1 = SimpleLocation(105, 150)
        >>> print(f1)
        [105:150]
        >>> print(f1 - 100)
        [5:50]

        This can be useful when editing annotation. You can also add an integer
        to a feature location (which shifts in the opposite direction).
        """
        if isinstance(other, int):
            return self._shift(-other)
        else:
            return NotImplemented

    def __nonzero__(self):
        """Return True regardless of the length of the feature.

        This behavior is for backwards compatibility, since until the
        __len__ method was added, a SimpleLocation always evaluated as True.

        Note that in comparison, Seq objects, strings, lists, etc, will all
        evaluate to False if they have length zero.

        WARNING: The SimpleLocation may in future evaluate to False when its
        length is zero (in order to better match normal python behavior)!
        """
        return True

    def __len__(self):
        """Return the length of the region described by the SimpleLocation object.

        Note that extra care may be needed for fuzzy locations, e.g.

        >>> from Bio.SeqFeature import SimpleLocation
        >>> from Bio.SeqFeature import BeforePosition, AfterPosition
        >>> loc = SimpleLocation(BeforePosition(5), AfterPosition(10))
        >>> len(loc)
        5
        """
        return int(self._end) - int(self._start)

    def __contains__(self, value):
        """Check if an integer position is within the SimpleLocation object.

        Note that extra care may be needed for fuzzy locations, e.g.

        >>> from Bio.SeqFeature import SimpleLocation
        >>> from Bio.SeqFeature import BeforePosition, AfterPosition
        >>> loc = SimpleLocation(BeforePosition(5), AfterPosition(10))
        >>> len(loc)
        5
        >>> [i for i in range(15) if i in loc]
        [5, 6, 7, 8, 9]
        """
        if not isinstance(value, int):
            raise ValueError(
                "Currently we only support checking for integer "
                "positions being within a SimpleLocation."
            )
        if value < self._start or value >= self._end:
            return False
        else:
            return True

    def __iter__(self):
        """Iterate over the parent positions within the SimpleLocation object.

        >>> from Bio.SeqFeature import SimpleLocation
        >>> from Bio.SeqFeature import BeforePosition, AfterPosition
        >>> loc = SimpleLocation(BeforePosition(5), AfterPosition(10))
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

        >>> loc = SimpleLocation(BeforePosition(5), AfterPosition(10), strand = -1)
        >>> list(loc)
        [9, 8, 7, 6, 5]
        """
        if self.strand == -1:
            yield from range(self._end - 1, self._start - 1, -1)
        else:
            yield from range(self._start, self._end)

    def __eq__(self, other):
        """Implement equality by comparing all the location attributes."""
        if not isinstance(other, SimpleLocation):
            return False
        return (
            self._start == other.start
            and self._end == other.end
            and self._strand == other.strand
            and self.ref == other.ref
            and self.ref_db == other.ref_db
        )

    def _shift(self, offset):
        """Return a copy of the SimpleLocation shifted by an offset (PRIVATE).

        Returns self when location is relative to an external reference.
        """
        # TODO - What if offset is a fuzzy position?
        if self.ref or self.ref_db:
            return self
        return SimpleLocation(
            start=self._start + offset,
            end=self._end + offset,
            strand=self.strand,
        )

    def _flip(self, length):
        """Return a copy of the location after the parent is reversed (PRIVATE).

        Returns self when location is relative to an external reference.
        """
        if self.ref or self.ref_db:
            return self
        # Note this will flip the start and end too!
        if self.strand == +1:
            flip_strand = -1
        elif self.strand == -1:
            flip_strand = +1
        else:
            # 0 or None
            flip_strand = self.strand
        return SimpleLocation(
            start=self._end._flip(length),
            end=self._start._flip(length),
            strand=flip_strand,
        )

    @property
    def parts(self):
        """Read only list of sections (always one, the SimpleLocation object).

        This is a convenience property allowing you to write code handling
        both SimpleLocation objects (with one part) and more complex
        CompoundLocation objects (with multiple parts) interchangeably.
        """
        return [self]

    @property
    def start(self):
        """Start location - left most (minimum) value, regardless of strand.

        Read only, returns an integer like position object, possibly a fuzzy
        position.
        """
        return self._start

    @property
    def end(self):
        """End location - right most (maximum) value, regardless of strand.

        Read only, returns an integer like position object, possibly a fuzzy
        position.
        """
        return self._end

    def extract(self, parent_sequence, references=None):
        """Extract the sequence from supplied parent sequence using the SimpleLocation object.

        The parent_sequence can be a Seq like object or a string, and will
        generally return an object of the same type. The exception to this is
        a MutableSeq as the parent sequence will return a Seq object.
        If the location refers to other records, they must be supplied
        in the optional dictionary references.

        >>> from Bio.Seq import Seq
        >>> from Bio.SeqFeature import SimpleLocation
        >>> seq = Seq("MKQHKAMIVALIVICITAVVAAL")
        >>> feature_loc = SimpleLocation(8, 15)
        >>> feature_loc.extract(seq)
        Seq('VALIVIC')

        """
        if self.ref or self.ref_db:
            if not references:
                raise ValueError(
                    f"Feature references another sequence ({self.ref}),"
                    " references mandatory"
                )
            elif self.ref not in references:
                # KeyError?
                raise ValueError(
                    f"Feature references another sequence ({self.ref}),"
                    " not found in references"
                )
            parent_sequence = references[self.ref]
        f_seq = parent_sequence[int(self.start) : int(self.end)]
        if isinstance(f_seq, MutableSeq):
            f_seq = Seq(f_seq)
        if self.strand == -1:
            f_seq = reverse_complement(f_seq)
        return f_seq


FeatureLocation = SimpleLocation  # OBSOLETE; for backward compatability only.


class CompoundLocation(Location):
    """For handling joins etc where a feature location has several parts."""

    def __init__(self, parts, operator="join"):
        """Initialize the class.

        >>> from Bio.SeqFeature import SimpleLocation, CompoundLocation
        >>> f1 = SimpleLocation(10, 40, strand=+1)
        >>> f2 = SimpleLocation(50, 59, strand=+1)
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

        >>> f = CompoundLocation([SimpleLocation(3, 6, strand=+1),
        ...                       SimpleLocation(10, 13, strand=-1)])
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
        give the full span covered, 0 <= start <= end <= full sequence length.

        >>> f.start == min(f)
        True
        >>> f.end == max(f) + 1
        True

        This is consistent with the behavior of the SimpleLocation for a single
        region, where again the 'start' and 'end' do not necessarily give the
        biological start and end, but rather the 'minimal' and 'maximal'
        coordinate boundaries.

        Note that adding locations provides a more intuitive method of
        construction:

        >>> f = SimpleLocation(3, 6, strand=+1) + SimpleLocation(10, 13, strand=-1)
        >>> len(f)
        6
        >>> list(f)
        [3, 4, 5, 12, 11, 10]
        """
        self.operator = operator
        self.parts = list(parts)
        for loc in self.parts:
            if not isinstance(loc, SimpleLocation):
                raise ValueError(
                    "CompoundLocation should be given a list of "
                    "SimpleLocation objects, not %s" % loc.__class__
                )
        if len(parts) < 2:
            raise ValueError(
                f"CompoundLocation should have at least 2 parts, not {parts!r}"
            )

    def __str__(self):
        """Return a representation of the CompoundLocation object (with python counting)."""
        return "%s{%s}" % (self.operator, ", ".join(str(loc) for loc in self.parts))

    def __repr__(self):
        """Represent the CompoundLocation object as string for debugging."""
        return f"{self.__class__.__name__}({self.parts!r}, {self.operator!r})"

    def _get_strand(self):
        """Get function for the strand property (PRIVATE)."""
        # Historically a join on the reverse strand has been represented
        # in Biopython with both the parent SeqFeature and its children
        # (the exons for a CDS) all given a strand of -1.  Likewise, for
        # a join feature on the forward strand they all have strand +1.
        # However, we must also consider evil mixed strand examples like
        # this, join(complement(69611..69724),139856..140087,140625..140650)
        if len({loc.strand for loc in self.parts}) == 1:
            return self.parts[0].strand
        else:
            return None  # i.e. mixed strands

    def _set_strand(self, value):
        """Set function for the strand property (PRIVATE)."""
        # Should this be allowed/encouraged?
        for loc in self.parts:
            loc.strand = value

    strand = property(
        fget=_get_strand,
        fset=_set_strand,
        doc="""Overall strand of the compound location.

        If all the parts have the same strand, that is returned. Otherwise
        for mixed strands, this returns None.

        >>> from Bio.SeqFeature import SimpleLocation, CompoundLocation
        >>> f1 = SimpleLocation(15, 17, strand=1)
        >>> f2 = SimpleLocation(20, 30, strand=-1)
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

        """,
    )

    def __add__(self, other):
        """Combine locations, or shift the location by an integer offset.

        >>> from Bio.SeqFeature import SimpleLocation
        >>> f1 = SimpleLocation(15, 17) + SimpleLocation(20, 30)
        >>> print(f1)
        join{[15:17], [20:30]}

        You can add another SimpleLocation:

        >>> print(f1 + SimpleLocation(40, 50))
        join{[15:17], [20:30], [40:50]}
        >>> print(SimpleLocation(5, 10) + f1)
        join{[5:10], [15:17], [20:30]}

        You can also add another CompoundLocation:

        >>> f2 = SimpleLocation(40, 50) + SimpleLocation(60, 70)
        >>> print(f2)
        join{[40:50], [60:70]}
        >>> print(f1 + f2)
        join{[15:17], [20:30], [40:50], [60:70]}

        Also, as with the SimpleLocation, adding an integer shifts the
        location's coordinates by that offset:

        >>> print(f1 + 100)
        join{[115:117], [120:130]}
        >>> print(200 + f1)
        join{[215:217], [220:230]}
        >>> print(f1 + (-5))
        join{[10:12], [15:25]}
        """
        if isinstance(other, SimpleLocation):
            return CompoundLocation(self.parts + [other], self.operator)
        elif isinstance(other, CompoundLocation):
            if self.operator != other.operator:
                # Handle join+order -> order as a special case?
                raise ValueError(
                    f"Mixed operators {self.operator} and {other.operator}"
                )
            return CompoundLocation(self.parts + other.parts, self.operator)
        elif isinstance(other, int):
            return self._shift(other)
        else:
            raise NotImplementedError

    def __radd__(self, other):
        """Add a feature to the left."""
        if isinstance(other, SimpleLocation):
            return CompoundLocation([other] + self.parts, self.operator)
        elif isinstance(other, int):
            return self._shift(other)
        else:
            raise NotImplementedError

    def __contains__(self, value):
        """Check if an integer position is within the CompoundLocation object."""
        for loc in self.parts:
            if value in loc:
                return True
        return False

    def __nonzero__(self):
        """Return True regardless of the length of the feature.

        This behavior is for backwards compatibility, since until the
        __len__ method was added, a SimpleLocation always evaluated as True.

        Note that in comparison, Seq objects, strings, lists, etc, will all
        evaluate to False if they have length zero.

        WARNING: The SimpleLocation may in future evaluate to False when its
        length is zero (in order to better match normal python behavior)!
        """
        return True

    def __len__(self):
        """Return the length of the CompoundLocation object."""
        return sum(len(loc) for loc in self.parts)

    def __iter__(self):
        """Iterate over the parent positions within the CompoundLocation object."""
        for loc in self.parts:
            yield from loc

    def __eq__(self, other):
        """Check if all parts of CompoundLocation are equal to all parts of other CompoundLocation."""
        if not isinstance(other, CompoundLocation):
            return False
        if len(self.parts) != len(other.parts):
            return False
        if self.operator != other.operator:
            return False
        for self_part, other_part in zip(self.parts, other.parts):
            if self_part != other_part:
                return False
        return True

    def _shift(self, offset):
        """Return a copy of the CompoundLocation shifted by an offset (PRIVATE)."""
        return CompoundLocation(
            [loc._shift(offset) for loc in self.parts], self.operator
        )

    def _flip(self, length):
        """Return a copy of the locations after the parent is reversed (PRIVATE).

        Note that the order of the parts is NOT reversed too. Consider a CDS
        on the forward strand with exons small, medium and large (in length).
        Once we change the frame of reference to the reverse complement strand,
        the start codon is still part of the small exon, and the stop codon
        still part of the large exon - so the part order remains the same!

        Here is an artificial example, were the features map to the two upper
        case regions and the lower case runs of n are not used:

        >>> from Bio.Seq import Seq
        >>> from Bio.SeqFeature import SimpleLocation
        >>> dna = Seq("nnnnnAGCATCCTGCTGTACnnnnnnnnGAGAMTGCCATGCCCCTGGAGTGAnnnnn")
        >>> small = SimpleLocation(5, 20, strand=1)
        >>> large = SimpleLocation(28, 52, strand=1)
        >>> location = small + large
        >>> print(small)
        [5:20](+)
        >>> print(large)
        [28:52](+)
        >>> print(location)
        join{[5:20](+), [28:52](+)}
        >>> for part in location.parts:
        ...     print(len(part))
        ...
        15
        24

        As you can see, this is a silly example where each "exon" is a word:

        >>> print(small.extract(dna).translate())
        SILLY
        >>> print(large.extract(dna).translate())
        EXAMPLE*
        >>> print(location.extract(dna).translate())
        SILLYEXAMPLE*
        >>> for part in location.parts:
        ...     print(part.extract(dna).translate())
        ...
        SILLY
        EXAMPLE*

        Now, let's look at this from the reverse strand frame of reference:

        >>> flipped_dna = dna.reverse_complement()
        >>> flipped_location = location._flip(len(dna))
        >>> print(flipped_location.extract(flipped_dna).translate())
        SILLYEXAMPLE*
        >>> for part in flipped_location.parts:
        ...     print(part.extract(flipped_dna).translate())
        ...
        SILLY
        EXAMPLE*

        The key point here is the first part of the CompoundFeature is still the
        small exon, while the second part is still the large exon:

        >>> for part in flipped_location.parts:
        ...     print(len(part))
        ...
        15
        24
        >>> print(flipped_location)
        join{[37:52](-), [5:29](-)}

        Notice the parts are not reversed. However, there was a bug here in older
        versions of Biopython which would have given join{[5:29](-), [37:52](-)}
        and the translation would have wrongly been "EXAMPLE*SILLY" instead.

        """
        return CompoundLocation(
            [loc._flip(length) for loc in self.parts], self.operator
        )

    @property
    def start(self):
        """Start location - left most (minimum) value, regardless of strand.

        Read only, returns an integer like position object, possibly a fuzzy
        position.

        For the special case of a CompoundLocation wrapping the origin of a
        circular genome, this will return zero.
        """
        return min(loc.start for loc in self.parts)

    @property
    def end(self):
        """End location - right most (maximum) value, regardless of strand.

        Read only, returns an integer like position object, possibly a fuzzy
        position.

        For the special case of a CompoundLocation wrapping the origin of
        a circular genome this will match the genome length.
        """
        return max(loc.end for loc in self.parts)

    @property
    def ref(self):
        """Not present in CompoundLocation, dummy method for API compatibility."""
        return None

    @property
    def ref_db(self):
        """Not present in CompoundLocation, dummy method for API compatibility."""
        return None

    def extract(self, parent_sequence, references=None):
        """Extract the sequence from supplied parent sequence using the CompoundLocation object.

        The parent_sequence can be a Seq like object or a string, and will
        generally return an object of the same type. The exception to this is
        a MutableSeq as the parent sequence will return a Seq object.
        If the location refers to other records, they must be supplied
        in the optional dictionary references.

        >>> from Bio.Seq import Seq
        >>> from Bio.SeqFeature import SimpleLocation, CompoundLocation
        >>> seq = Seq("MKQHKAMIVALIVICITAVVAAL")
        >>> fl1 = SimpleLocation(2, 8)
        >>> fl2 = SimpleLocation(10, 15)
        >>> fl3 = CompoundLocation([fl1,fl2])
        >>> fl3.extract(seq)
        Seq('QHKAMILIVIC')

        """
        # This copes with mixed strand features & all on reverse:
        parts = [
            loc.extract(parent_sequence, references=references) for loc in self.parts
        ]
        f_seq = functools.reduce(lambda x, y: x + y, parts)
        return f_seq


class Position(ABC):
    """Abstract base class representing a position."""

    @abstractmethod
    def __repr__(self):
        """Represent the Position object as a string for debugging."""
        return f"{self.__class__.__name__}(...)"

    @staticmethod
    def fromstring(text, offset=0):
        """Build a Position object from the text string.

        For an end position, leave offset as zero (default):

        >>> Position.fromstring("5")
        ExactPosition(5)

        For a start position, set offset to minus one (for Python counting):

        >>> Position.fromstring("5", -1)
        ExactPosition(4)

        This also covers fuzzy positions:

        >>> p = Position.fromstring("<5")
        >>> p
        BeforePosition(5)
        >>> print(p)
        <5
        >>> int(p)
        5

        >>> Position.fromstring(">5")
        AfterPosition(5)

        By default assumes an end position, so note the integer behavior:

        >>> p = Position.fromstring("one-of(5,8,11)")
        >>> p
        OneOfPosition(11, choices=[ExactPosition(5), ExactPosition(8), ExactPosition(11)])
        >>> print(p)
        one-of(5,8,11)
        >>> int(p)
        11

        >>> Position.fromstring("(8.10)")
        WithinPosition(10, left=8, right=10)

        Fuzzy start positions:

        >>> p = Position.fromstring("<5", -1)
        >>> p
        BeforePosition(4)
        >>> print(p)
        <4
        >>> int(p)
        4

        Notice how the integer behavior changes too!

        >>> p = Position.fromstring("one-of(5,8,11)", -1)
        >>> p
        OneOfPosition(4, choices=[ExactPosition(4), ExactPosition(7), ExactPosition(10)])
        >>> print(p)
        one-of(4,7,10)
        >>> int(p)
        4

        """
        if offset != 0 and offset != -1:
            raise ValueError(
                "To convert one-based indices to zero-based indices, offset must be either 0 (for end positions) or -1 (for start positions)."
            )
        if text == "?":
            return UnknownPosition()
        if text.startswith("?"):
            return UncertainPosition(int(text[1:]) + offset)
        if text.startswith("<"):
            return BeforePosition(int(text[1:]) + offset)
        if text.startswith(">"):
            return AfterPosition(int(text[1:]) + offset)
        m = _re_within_position.match(text)
        if m is not None:
            s, e = m.groups()
            s = int(s) + offset
            e = int(e) + offset
            if offset == -1:
                default = s
            else:
                default = e
            return WithinPosition(default, left=s, right=e)
        m = _re_oneof_position.match(text)
        if m is not None:
            positions = m.groups()[0]
            parts = [ExactPosition(int(pos) + offset) for pos in positions.split(",")]
            if offset == -1:
                default = min(int(pos) for pos in parts)
            else:
                default = max(int(pos) for pos in parts)
            return OneOfPosition(default, choices=parts)
        return ExactPosition(int(text) + offset)


class ExactPosition(int, Position):
    """Specify the specific position of a boundary.

    Arguments:
     - position - The position of the boundary.
     - extension - An optional argument which must be zero since we don't
       have an extension. The argument is provided so that the same number
       of arguments can be passed to all position types.

    In this case, there is no fuzziness associated with the position.

    >>> p = ExactPosition(5)
    >>> p
    ExactPosition(5)
    >>> print(p)
    5

    >>> isinstance(p, Position)
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
    ExactPosition(15)

    """

    def __new__(cls, position, extension=0):
        """Create an ExactPosition object."""
        if extension != 0:
            raise AttributeError(f"Non-zero extension {extension} for exact position.")
        return int.__new__(cls, position)

    # Must define this on Python 3.8 onwards because we redefine __repr__
    def __str__(self):
        """Return a representation of the ExactPosition object (with python counting)."""
        return str(int(self))

    def __repr__(self):
        """Represent the ExactPosition object as a string for debugging."""
        return "%s(%i)" % (self.__class__.__name__, int(self))

    def __add__(self, offset):
        """Return a copy of the position object with its location shifted (PRIVATE)."""
        # By default preserve any subclass
        return self.__class__(int(self) + offset)

    def _flip(self, length):
        """Return a copy of the location after the parent is reversed (PRIVATE)."""
        # By default preserve any subclass
        return self.__class__(length - int(self))


class UncertainPosition(ExactPosition):
    """Specify a specific position which is uncertain.

    This is used in UniProt, e.g. ?222 for uncertain position 222, or in the
    XML format explicitly marked as uncertain. Does not apply to GenBank/EMBL.
    """


class UnknownPosition(Position):
    """Specify a specific position which is unknown (has no position).

    This is used in UniProt, e.g. ? or in the XML as unknown.
    """

    def __repr__(self):
        """Represent the UnknownPosition object as a string for debugging."""
        return f"{self.__class__.__name__}()"

    def __hash__(self):
        """Return the hash value of the UnknownPosition object."""
        return hash(None)

    def __add__(self, offset):
        """Return a copy of the position object with its location shifted (PRIVATE)."""
        return self

    def _flip(self, length):
        """Return a copy of the location after the parent is reversed (PRIVATE)."""
        return self


class WithinPosition(int, Position):
    """Specify the position of a boundary within some coordinates.

    Arguments:
    - position - The default integer position
    - left - The start (left) position of the boundary
    - right - The end (right) position of the boundary

    This allows dealing with a location like ((11.14)..100). This
    indicates that the start of the sequence is somewhere between 11
    and 14. Since this is a start coordinate, it should act like
    it is at position 11 (or in Python counting, 10).

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
    WithinPosition(20, left=20, right=23)

    >>> isinstance(p, WithinPosition)
    True
    >>> isinstance(p, Position)
    True
    >>> isinstance(p, int)
    True

    Note this also applies for comparison to other position objects,
    where again the integer behavior is used:

    >>> p == 10
    True
    >>> p == ExactPosition(10)
    True
    >>> p == BeforePosition(10)
    True
    >>> p == AfterPosition(10)
    True

    If this were an end point, you would want the position to be 13
    (the right/larger value, not the left/smaller value as above):

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

    """

    def __new__(cls, position, left, right):
        """Create a WithinPosition object."""
        if not (position == left or position == right):
            raise RuntimeError(
                "WithinPosition: %r should match left %r or "
                "right %r" % (position, left, right)
            )
        obj = int.__new__(cls, position)
        obj._left = left
        obj._right = right
        return obj

    def __getnewargs__(self):
        """Return the arguments accepted by __new__.

        Necessary to allow pickling and unpickling of class instances.
        """
        return (int(self), self._left, self._right)

    def __repr__(self):
        """Represent the WithinPosition object as a string for debugging."""
        return "%s(%i, left=%i, right=%i)" % (
            self.__class__.__name__,
            int(self),
            self._left,
            self._right,
        )

    def __str__(self):
        """Return a representation of the WithinPosition object (with python counting)."""
        return f"({self._left}.{self._right})"

    def __add__(self, offset):
        """Return a copy of the position object with its location shifted."""
        return self.__class__(
            int(self) + offset, self._left + offset, self._right + offset
        )

    def _flip(self, length):
        """Return a copy of the location after the parent is reversed (PRIVATE)."""
        return self.__class__(
            length - int(self), length - self._right, length - self._left
        )


class BetweenPosition(int, Position):
    """Specify the position of a boundary between two coordinates (OBSOLETE?).

    Arguments:
     - position - The default integer position
     - left - The start (left) position of the boundary
     - right - The end (right) position of the boundary

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
    >>> int(p) == int(p2)
    False
    >>> p == 456
    True
    >>> p2 == 123
    True

    Note this potentially surprising behavior:

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
        """Create a new instance in BetweenPosition object."""
        assert position == left or position == right
        # TODO - public API for getting left/right, especially the unknown one
        obj = int.__new__(cls, position)
        obj._left = left
        obj._right = right
        return obj

    def __getnewargs__(self):
        """Return the arguments accepted by __new__.

        Necessary to allow pickling and unpickling of class instances.
        """
        return (int(self), self._left, self._right)

    def __repr__(self):
        """Represent the BetweenPosition object as a string for debugging."""
        return "%s(%i, left=%i, right=%i)" % (
            self.__class__.__name__,
            int(self),
            self._left,
            self._right,
        )

    def __str__(self):
        """Return a representation of the BetweenPosition object (with python counting)."""
        return f"({self._left}^{self._right})"

    def __add__(self, offset):
        """Return a copy of the position object with its location shifted (PRIVATE)."""
        return self.__class__(
            int(self) + offset, self._left + offset, self._right + offset
        )

    def _flip(self, length):
        """Return a copy of the location after the parent is reversed (PRIVATE)."""
        return self.__class__(
            length - int(self), length - self._right, length - self._left
        )


class BeforePosition(int, Position):
    """Specify a position where the actual location occurs before it.

    Arguments:
     - position - The upper boundary of where the location can occur.
     - extension - An optional argument which must be zero since we don't
       have an extension. The argument is provided so that the same number
       of arguments can be passed to all position types.

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
    BeforePosition(15)

    Note this potentially surprising behavior:

    >>> p == ExactPosition(5)
    True
    >>> p == AfterPosition(5)
    True

    Just remember that for equality and sorting the position objects act
    like integers.
    """

    # Subclasses int so can't use __init__
    def __new__(cls, position, extension=0):
        """Create a new instance in BeforePosition object."""
        if extension != 0:
            raise AttributeError(f"Non-zero extension {extension} for exact position.")
        return int.__new__(cls, position)

    def __repr__(self):
        """Represent the location as a string for debugging."""
        return "%s(%i)" % (self.__class__.__name__, int(self))

    def __str__(self):
        """Return a representation of the BeforePosition object (with python counting)."""
        return f"<{int(self)}"

    def __add__(self, offset):
        """Return a copy of the position object with its location shifted (PRIVATE)."""
        return self.__class__(int(self) + offset)

    def _flip(self, length):
        """Return a copy of the location after the parent is reversed (PRIVATE)."""
        return AfterPosition(length - int(self))


class AfterPosition(int, Position):
    """Specify a position where the actual location is found after it.

    Arguments:
     - position - The lower boundary of where the location can occur.
     - extension - An optional argument which must be zero since we don't
       have an extension. The argument is provided so that the same number
       of arguments can be passed to all position types.

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
    AfterPosition(17)

    >>> isinstance(p, AfterPosition)
    True
    >>> isinstance(p, Position)
    True
    >>> isinstance(p, int)
    True

    Note this potentially surprising behavior:

    >>> p == ExactPosition(7)
    True
    >>> p == BeforePosition(7)
    True

    Just remember that for equality and sorting the position objects act
    like integers.
    """

    # Subclasses int so can't use __init__
    def __new__(cls, position, extension=0):
        """Create a new instance of the AfterPosition object."""
        if extension != 0:
            raise AttributeError(f"Non-zero extension {extension} for exact position.")
        return int.__new__(cls, position)

    def __repr__(self):
        """Represent the location as a string for debugging."""
        return "%s(%i)" % (self.__class__.__name__, int(self))

    def __str__(self):
        """Return a representation of the AfterPosition object (with python counting)."""
        return f">{int(self)}"

    def __add__(self, offset):
        """Return a copy of the position object with its location shifted (PRIVATE)."""
        return self.__class__(int(self) + offset)

    def _flip(self, length):
        """Return a copy of the location after the parent is reversed (PRIVATE)."""
        return BeforePosition(length - int(self))


class OneOfPosition(int, Position):
    """Specify a position where the location can be multiple positions.

    This models the GenBank 'one-of(1888,1901)' function, and tries
    to make this fit within the Biopython Position models. If this was
    a start position it should act like 1888, but as an end position 1901.

    >>> p = OneOfPosition(1888, [ExactPosition(1888), ExactPosition(1901)])
    >>> p
    OneOfPosition(1888, choices=[ExactPosition(1888), ExactPosition(1901)])
    >>> int(p)
    1888

    Integer comparisons and operators act like using int(p),

    >>> p == 1888
    True
    >>> p <= 1888
    True
    >>> p > 1888
    False
    >>> p + 100
    OneOfPosition(1988, choices=[ExactPosition(1988), ExactPosition(2001)])

    >>> isinstance(p, OneOfPosition)
    True
    >>> isinstance(p, Position)
    True
    >>> isinstance(p, int)
    True

    """

    def __new__(cls, position, choices):
        """Initialize with a set of possible positions.

        choices is a list of Position derived objects, specifying possible
        locations.

        position is an integer specifying the default behavior.
        """
        if position not in choices:
            raise ValueError(
                f"OneOfPosition: {position!r} should match one of {choices!r}"
            )
        obj = int.__new__(cls, position)
        obj.position_choices = choices
        return obj

    def __getnewargs__(self):
        """Return the arguments accepted by __new__.

        Necessary to allow pickling and unpickling of class instances.
        """
        return (int(self), self.position_choices)

    def __repr__(self):
        """Represent the OneOfPosition object as a string for debugging."""
        return "%s(%i, choices=%r)" % (
            self.__class__.__name__,
            int(self),
            self.position_choices,
        )

    def __str__(self):
        """Return a representation of the OneOfPosition object (with python counting)."""
        out = "one-of("
        for position in self.position_choices:
            out += f"{position},"
        # replace the last comma with the closing parenthesis
        return out[:-1] + ")"

    def __add__(self, offset):
        """Return a copy of the position object with its location shifted (PRIVATE)."""
        return self.__class__(
            int(self) + offset, [p + offset for p in self.position_choices]
        )

    def _flip(self, length):
        """Return a copy of the location after the parent is reversed (PRIVATE)."""
        return self.__class__(
            length - int(self), [p._flip(length) for p in self.position_choices[::-1]]
        )


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
