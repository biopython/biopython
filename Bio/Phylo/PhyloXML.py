# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Classes corresponding to phyloXML elements.

See Also
--------
Official specification:
   http://phyloxml.org/
Journal article:
    Han and Zmasek (2009), https://doi.org/10.1186/1471-2105-10-356

"""

import re
import warnings

from Bio.Align import Alignment, MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from Bio import BiopythonWarning

from Bio.Phylo import BaseTree


class PhyloXMLWarning(BiopythonWarning):
    """Warning for non-compliance with the phyloXML specification."""

    pass


def _check_str(text, testfunc):
    """Check a string using testfunc, and warn if there's no match (PRIVATE)."""
    if text is not None and not testfunc(text):
        warnings.warn(
            f"String {text} doesn't match the given regexp",
            PhyloXMLWarning,
            stacklevel=2,
        )


# Core elements


class PhyloElement(BaseTree.TreeElement):
    """Base class for all PhyloXML objects."""


class Phyloxml(PhyloElement):
    """Root node of the PhyloXML document.

    Contains an arbitrary number of Phylogeny elements, possibly followed by
    elements from other namespaces.

    :Parameters:
        attributes : dict
            (XML namespace definitions)
        phylogenies : list
            The phylogenetic trees
        other : list
            Arbitrary non-phyloXML elements, if any

    """

    def __init__(self, attributes, phylogenies=None, other=None):
        """Initialize parameters for PhyloXML object."""
        self.attributes = {
            # standard
            "xmlns:xsi": "http://www.w3.org/2001/XMLSchema-instance",
            "xmlns": "http://www.phyloxml.org",
            "xsi:schemaLocation": "http://www.phyloxml.org http://www.phyloxml.org/1.10/phyloxml.xsd",
        }
        if attributes:
            self.attributes.update(attributes)
        self.phylogenies = phylogenies or []
        self.other = other or []

    def __getitem__(self, index):
        """Get a phylogeny by index or name."""
        if isinstance(index, (int, slice)):
            return self.phylogenies[index]
        if not isinstance(index, str):
            raise KeyError(f"can't use {type(index)} as an index")
        for tree in self.phylogenies:
            if tree.name == index:
                return tree
        else:
            raise KeyError(f"no phylogeny found with name {index!r}")

    def __iter__(self):
        """Iterate through the phylogenetic trees in this object."""
        return iter(self.phylogenies)

    def __len__(self):
        """Return the number of phylogenetic trees in this object."""
        return len(self.phylogenies)

    def __str__(self):
        """Return name of phylogenies in the object."""
        return "%s([%s])" % (
            self.__class__.__name__,
            ",\n".join(map(str, self.phylogenies)),
        )


class Other(PhyloElement):
    """Container for non-phyloXML elements in the tree.

    Usually, an Other object will have either a 'value' or a non-empty list
    of 'children', but not both. This is not enforced here, though.

    :Parameters:
        tag : string
            local tag for the XML node
        namespace : string
            XML namespace for the node -- should not be the default phyloXML
            namespace.
        attributes : dict of strings
            attributes on the XML node
        value : string
            text contained directly within this XML node
        children : list
            child nodes, if any (also ``Other`` instances)

    """

    def __init__(self, tag, namespace=None, attributes=None, value=None, children=None):
        """Initialize values for non-phyloXML elements."""
        self.tag = tag
        self.namespace = namespace
        self.attributes = attributes or {}
        self.value = value
        self.children = children or []

    def __iter__(self):
        """Iterate through the children of this object (if any)."""
        return iter(self.children)


class Phylogeny(PhyloElement, BaseTree.Tree):
    """A phylogenetic tree.

    :Parameters:
        root : Clade
            the root node/clade of this tree
        rooted : bool
            True if this tree is rooted
        rerootable : bool
            True if this tree is rerootable
        branch_length_unit : string
            unit for branch_length values on clades
        name : string
            identifier for this tree, not required to be unique
        id : Id
            unique identifier for this tree
        description : string
            plain-text description
        date : Date
            date for the root node of this tree
        confidences : list
            Confidence objects for this tree
        clade_relations : list
            CladeRelation objects
        sequence_relations : list
            SequenceRelation objects
        properties : list
            Property objects
        other : list
            non-phyloXML elements (type ``Other``)

    """

    def __init__(
        self,
        root=None,
        rooted=True,
        rerootable=None,
        branch_length_unit=None,
        type=None,
        # Child nodes
        name=None,
        id=None,
        description=None,
        date=None,
        # Collections
        confidences=None,
        clade_relations=None,
        sequence_relations=None,
        properties=None,
        other=None,
    ):
        """Initialize values for phylogenetic tree object."""
        assert isinstance(rooted, bool)
        self.root = root
        self.rooted = rooted
        self.rerootable = rerootable
        self.branch_length_unit = branch_length_unit
        self.type = type
        self.name = name
        self.id = id
        self.description = description
        self.date = date
        self.confidences = confidences or []
        self.clade_relations = clade_relations or []
        self.sequence_relations = sequence_relations or []
        self.properties = properties or []
        self.other = other or []

    @classmethod
    def from_tree(cls, tree, **kwargs):
        """Create a new Phylogeny given a Tree (from Newick/Nexus or BaseTree).

        Keyword arguments are the usual ``Phylogeny`` constructor parameters.
        """
        phy = cls(
            root=Clade.from_clade(tree.root),
            rooted=tree.rooted,
            name=tree.name,
            id=(tree.id is not None) and Id(str(tree.id)) or None,
        )
        phy.__dict__.update(kwargs)
        return phy

    @classmethod
    def from_clade(cls, clade, **kwargs):
        """Create a new Phylogeny given a Newick or BaseTree Clade object.

        Keyword arguments are the usual ``PhyloXML.Clade`` constructor parameters.
        """
        return Clade.from_clade(clade).to_phylogeny(**kwargs)

    def as_phyloxml(self):
        """Return this tree, a PhyloXML-compatible Phylogeny object.

        Overrides the ``BaseTree`` method.
        """
        return self

    def to_phyloxml_container(self, **kwargs):
        """Create a new Phyloxml object containing just this phylogeny."""
        return Phyloxml(kwargs, phylogenies=[self])

    def to_alignment(self):
        """Construct a MultipleSeqAlignment from the aligned sequences in this tree."""

        def is_aligned_seq(elem):
            if isinstance(elem, Sequence) and elem.mol_seq.is_aligned:
                return True
            return False

        seqs = self._filter_search(is_aligned_seq, "preorder", True)
        records = (seq.to_seqrecord() for seq in seqs)
        return MultipleSeqAlignment(records)

    @property
    def alignment(self):
        """Construct an Alignment object from the aligned sequences in this tree."""

        def is_aligned_seq(elem):
            if isinstance(elem, Sequence) and elem.mol_seq.is_aligned:
                return True
            return False

        seqs = self._filter_search(is_aligned_seq, "preorder", True)
        records = []
        lines = []
        for seq in seqs:
            record = seq.to_seqrecord()
            lines.append(str(record.seq))
            record.seq = record.seq.replace("-", "")
            records.append(record)
        if lines:
            coordinates = Alignment.infer_coordinates(lines)
        else:
            coordinates = None
        return Alignment(records, coordinates)

    # Singular property for plural attribute
    def _get_confidence(self):
        """Equivalent to self.confidences[0] if there is only 1 value (PRIVATE).

        See Also: ``Clade.confidence``, ``Clade.taxonomy``

        """
        if len(self.confidences) == 0:
            return None
        if len(self.confidences) > 1:
            raise AttributeError(
                "more than 1 confidence value available; use Phylogeny.confidences"
            )
        return self.confidences[0]

    def _set_confidence(self, value):
        if value is None:
            # Special case: mirror the behavior of _get_confidence
            self.confidences = []
            return
        if isinstance(value, (float, int)):
            value = Confidence(value)
        elif not isinstance(value, Confidence):
            raise ValueError("value must be a number or Confidence instance")
        if len(self.confidences) == 0:
            self.confidences.append(value)
        elif len(self.confidences) == 1:
            self.confidences[0] = value
        else:
            raise ValueError(
                "multiple confidence values already exist; "
                "use Phylogeny.confidences instead"
            )

    def _del_confidence(self):
        self.confidences = []

    confidence = property(_get_confidence, _set_confidence, _del_confidence)


class Clade(PhyloElement, BaseTree.Clade):
    """Describes a branch of the current phylogenetic tree.

    Used recursively, describes the topology of a phylogenetic tree.

    Both ``color`` and ``width`` elements should be interpreted by client code
    as applying to the whole clade, including all descendents, unless
    overwritten in-sub clades. This module doesn't automatically assign these
    attributes to sub-clades to achieve this cascade -- and neither should you.

    :Parameters:
        branch_length
            parent branch length of this clade
        id_source
            link other elements to a clade (on the xml-level)
        name : string
            short label for this clade
        confidences : list of Confidence objects
            used to indicate the support for a clade/parent branch.
        width : float
            branch width for this clade (including branch from parent)
        color : BranchColor
            color used for graphical display of this clade
        node_id
            unique identifier for the root node of this clade
        taxonomies : list
            Taxonomy objects
        sequences : list
            Sequence objects
        events : Events
            describe such events as gene-duplications at the root node/parent
            branch of this clade
        binary_characters : BinaryCharacters
            binary characters
        distributions : list of Distribution objects
            distribution(s) of this clade
        date : Date
            a date for the root node of this clade
        references : list
            Reference objects
        properties : list
            Property objects
        clades : list Clade objects
            Sub-clades
        other : list of Other objects
            non-phyloXML objects

    """

    def __init__(
        self,
        # Attributes
        branch_length=None,
        id_source=None,
        # Child nodes
        name=None,
        width=None,
        color=None,
        node_id=None,
        events=None,
        binary_characters=None,
        date=None,
        # Collections
        confidences=None,
        taxonomies=None,
        sequences=None,
        distributions=None,
        references=None,
        properties=None,
        clades=None,
        other=None,
    ):
        """Initialize value for the Clade object."""
        self.branch_length = branch_length
        self.id_source = id_source
        self.name = name
        self.width = width
        self.color = color
        self.node_id = node_id
        self.events = events
        self.binary_characters = binary_characters
        self.date = date
        self.confidences = confidences or []
        self.taxonomies = taxonomies or []
        self.sequences = sequences or []
        self.distributions = distributions or []
        self.references = references or []
        self.properties = properties or []
        self.clades = clades or []
        self.other = other or []

    @classmethod
    def from_clade(cls, clade, **kwargs):
        """Create a new PhyloXML Clade from a Newick or BaseTree Clade object.

        Keyword arguments are the usual PhyloXML Clade constructor parameters.
        """
        new_clade = cls(branch_length=clade.branch_length, name=clade.name)
        new_clade.clades = [cls.from_clade(c) for c in clade]
        new_clade.confidence = clade.confidence
        new_clade.width = clade.width
        new_clade.color = (
            BranchColor(clade.color.red, clade.color.green, clade.color.blue)
            if clade.color
            else None
        )
        new_clade.__dict__.update(kwargs)
        return new_clade

    def to_phylogeny(self, **kwargs):
        """Create a new phylogeny containing just this clade."""
        phy = Phylogeny(root=self, date=self.date)
        phy.__dict__.update(kwargs)
        return phy

    # Shortcuts for list attributes that are usually only 1 item
    # NB: Duplicated from Phylogeny class
    def _get_confidence(self):
        """Return confidence values (PRIVATE)."""
        if len(self.confidences) == 0:
            return None
        if len(self.confidences) > 1:
            raise AttributeError(
                "more than 1 confidence value available; use Clade.confidences"
            )
        return self.confidences[0]

    def _set_confidence(self, value):
        """Set the confidence value (PRIVATE)."""
        if value is None:
            # Special case: mirror the behavior of _get_confidence
            self.confidences = []
            return
        if isinstance(value, (float, int)):
            value = Confidence(value)
        elif not isinstance(value, Confidence):
            raise ValueError("value must be a number or Confidence instance")
        if len(self.confidences) == 0:
            self.confidences.append(value)
        elif len(self.confidences) == 1:
            self.confidences[0] = value
        else:
            raise ValueError(
                "multiple confidence values already exist; "
                "use Phylogeny.confidences instead"
            )

    def _del_confidence(self):
        """Delete confidences values (PRIVATE)."""
        self.confidences = []

    confidence = property(_get_confidence, _set_confidence, _del_confidence)

    def _get_taxonomy(self):
        """Get taxonomy list for the clade (PRIVATE)."""
        if len(self.taxonomies) == 0:
            return None
        if len(self.taxonomies) > 1:
            raise AttributeError(
                "more than 1 taxonomy value available; use Clade.taxonomies"
            )
        return self.taxonomies[0]

    def _set_taxonomy(self, value):
        """Set a taxonomy for the clade (PRIVATE)."""
        if not isinstance(value, Taxonomy):
            raise ValueError("assigned value must be a Taxonomy instance")
        if len(self.taxonomies) == 0:
            self.taxonomies.append(value)
        elif len(self.taxonomies) == 1:
            self.taxonomies[0] = value
        else:
            raise ValueError(
                "multiple taxonomy values already exist; "
                "use Phylogeny.taxonomies instead"
            )

    taxonomy = property(_get_taxonomy, _set_taxonomy)


# PhyloXML wrapper for a special BaseTree attribute


class BranchColor(PhyloElement, BaseTree.BranchColor):
    """Manage Tree branch's color."""

    def __init__(self, *args, **kwargs):
        """Initialize parameters for the BranchColor object."""
        BaseTree.BranchColor.__init__(self, *args, **kwargs)


# PhyloXML-specific complex types


class Accession(PhyloElement):
    """Captures the local part in a sequence identifier.

    Example: In ``UniProtKB:P17304``, the Accession instance attribute ``value``
    is 'P17304' and the ``source`` attribute is 'UniProtKB'.
    """

    def __init__(self, value, source):
        """Initialize value for Accession object."""
        self.value = value
        self.source = source

    def __str__(self):
        """Show the class name and an identifying attribute."""
        return f"{self.source}:{self.value}"


class Annotation(PhyloElement):
    """The annotation of a molecular sequence.

    It is recommended to annotate by using the optional 'ref' attribute.

    :Parameters:
        ref : string
            reference string, e.g. 'GO:0008270',
            'KEGG:Tetrachloroethene degradation', 'EC:1.1.1.1'
        source : string
            plain-text source for this annotation
        evidence : str
            describe evidence as free text (e.g. 'experimental')
        desc : string
            free text description
        confidence : Confidence
            state the type and value of support (type Confidence)
        properties : list
            typed and referenced annotations from external resources
        uri : Uri
            link

    """

    re_ref = re.compile(r"[a-zA-Z0-9_]+:[a-zA-Z0-9_\.\-\s]+")

    def __init__(
        self,
        # Attributes
        ref=None,
        source=None,
        evidence=None,
        type=None,
        # Child nodes
        desc=None,
        confidence=None,
        uri=None,
        # Collection
        properties=None,
    ):
        """Initialize value for the Annotation object."""
        _check_str(ref, self.re_ref.match)
        self.ref = ref
        self.source = source
        self.evidence = evidence
        self.type = type
        self.desc = desc
        self.confidence = confidence
        self.uri = uri
        self.properties = properties or []


class BinaryCharacters(PhyloElement):
    """Binary characters at the root of a clade.

    The names and/or counts of binary characters present, gained, and lost
    at the root of a clade.
    """

    def __init__(
        self,
        # Attributes
        type=None,
        gained_count=None,
        lost_count=None,
        present_count=None,
        absent_count=None,
        # Child nodes (flattened into collections)
        gained=None,
        lost=None,
        present=None,
        absent=None,
    ):
        """Initialize values for the BinaryCharacters object."""
        self.type = type
        self.gained_count = gained_count
        self.lost_count = lost_count
        self.present_count = present_count
        self.absent_count = absent_count
        self.gained = gained or []
        self.lost = lost or []
        self.present = present or []
        self.absent = absent or []


class CladeRelation(PhyloElement):
    """Expresses a typed relationship between two clades.

    For example, this could be used to describe multiple parents of a clade.

    :type id_ref_0: str
    :type id_ref_1: str
    :type distance: str
    :type type: str

    :type confidence: Confidence
    """

    def __init__(self, type, id_ref_0, id_ref_1, distance=None, confidence=None):
        """Initialize values for the CladeRelation object."""
        self.distance = distance
        self.type = type
        self.id_ref_0 = id_ref_0
        self.id_ref_1 = id_ref_1
        self.confidence = confidence


class Confidence(float, PhyloElement):
    """A general purpose confidence element.

    For example, this can be used to express the bootstrap support value of a
    clade (in which case the ``type`` attribute is 'bootstrap').

    :Parameters:
        value : float
            confidence value
        type : string
            label for the type of confidence, e.g. 'bootstrap'

    """

    def __new__(cls, value, type="unknown"):
        """Create and return a Confidence object with the specified value and type."""
        obj = super(Confidence, cls).__new__(cls, value)
        obj.type = type
        return obj

    @property
    def value(self):
        """Return the float value of the Confidence object."""
        return float(self)


class Date(PhyloElement):
    """A date associated with a clade/node.

    Its value can be numerical by using the 'value' element and/or free text
    with the 'desc' element' (e.g. 'Silurian'). If a numerical value is used, it
    is recommended to employ the 'unit' attribute.

    :Parameters:
        unit : string
            type of numerical value (e.g. 'mya' for 'million years ago')
        value : float
            the date value
        desc : string
            plain-text description of the date
        minimum : float
            lower bound on the date value
        maximum : float
            upper bound on the date value

    """

    def __init__(self, value=None, unit=None, desc=None, minimum=None, maximum=None):
        """Initialize values of the Date object."""
        self.value = value
        self.unit = unit
        self.desc = desc
        self.minimum = minimum
        self.maximum = maximum

    def __str__(self):
        """Show the class name and the human-readable date."""
        if self.unit and self.value is not None:
            return f"{self.value} {self.unit}"
        if self.desc is not None:
            return self.desc
        return self.__class__.__name__


class Distribution(PhyloElement):
    """Geographic distribution of the items of a clade (species, sequences).

    Intended for phylogeographic applications.

    :Parameters:
        desc : string
            free-text description of the location
        points : list of ``Point`` objects
            coordinates (similar to the 'Point' element in Google's KML format)
        polygons : list of ``Polygon`` objects
            coordinate sets defining geographic regions

    """

    def __init__(self, desc=None, points=None, polygons=None):
        """Initialize values of Distribution object."""
        self.desc = desc
        self.points = points or []
        self.polygons = polygons or []


class DomainArchitecture(PhyloElement):
    """Domain architecture of a protein.

    :Parameters:
        length : int
            total length of the protein sequence
        domains : list ProteinDomain objects
            the domains within this protein

    """

    def __init__(self, length=None, domains=None):
        """Initialize values of the DomainArchitecture object."""
        self.length = length
        self.domains = domains


class Events(PhyloElement):
    """Events at the root node of a clade (e.g. one gene duplication).

    All attributes are set to None by default, but this object can also be
    treated as a dictionary, in which case None values are treated as missing
    keys and deleting a key resets that attribute's value back to None.
    """

    ok_type = {
        "transfer",
        "fusion",
        "speciation_or_duplication",
        "other",
        "mixed",
        "unassigned",
    }

    def __init__(
        self,
        type=None,
        duplications=None,
        speciations=None,
        losses=None,
        confidence=None,
    ):
        """Initialize values of the Events object."""
        _check_str(type, self.ok_type.__contains__)
        self.type = type
        self.duplications = duplications
        self.speciations = speciations
        self.losses = losses
        self.confidence = confidence

    def items(self):
        """Return Event's items."""
        return [(k, v) for k, v in self.__dict__.items() if v is not None]

    def keys(self):
        """Return Event's keys."""
        return [k for k, v in self.__dict__.items() if v is not None]

    def values(self):
        """Return values from a key-value pair in an Events dict."""
        return [v for v in self.__dict__.values() if v is not None]

    def __len__(self):
        """Return number of Events."""
        # TODO - Better way to do this?
        return len(self.values())

    def __getitem__(self, key):
        """Get value of Event with the given key."""
        try:
            val = getattr(self, key)
        except AttributeError:
            raise KeyError(key) from None
        if val is None:
            raise KeyError(f"{key!r} has not been set in this object")
        return val

    def __setitem__(self, key, val):
        """Add item to Event dict."""
        setattr(self, key, val)

    def __delitem__(self, key):
        """Delete Event with given key."""
        setattr(self, key, None)

    def __iter__(self):
        """Iterate over the keys present in a Events dict."""
        return iter(self.keys())

    def __contains__(self, key):
        """Return True if Event dict contains key."""
        try:
            return getattr(self, key) is not None
        except AttributeError:
            return False


class Id(PhyloElement):
    """A general-purpose identifier element.

    Allows to indicate the provider (or authority) of an identifier, e.g. NCBI,
    along with the value itself.
    """

    def __init__(self, value, provider=None):
        """Initialize values for the identifier object."""
        self.value = value
        self.provider = provider

    def __str__(self):
        """Return identifier as a string."""
        if self.provider is not None:
            return f"{self.provider}:{self.value}"
        return self.value


class MolSeq(PhyloElement):
    """Store a molecular sequence.

    :Parameters:
        value : string
            the sequence itself
        is_aligned : bool
            True if this sequence is aligned with the others (usually meaning
            all aligned seqs are the same length and gaps may be present)

    """

    re_value = re.compile(r"[a-zA-Z\.\-\?\*_]+")

    def __init__(self, value, is_aligned=None):
        """Initialize parameters for the MolSeq object."""
        _check_str(value, self.re_value.match)
        self.value = value
        self.is_aligned = is_aligned

    def __str__(self):
        """Return the value of the Molecular Sequence object."""
        return self.value


class Point(PhyloElement):
    """Geographic coordinates of a point, with an optional altitude.

    Used by element 'Distribution'.

    :Parameters:
        geodetic_datum : string, required
            the geodetic datum (also called 'map datum'). For example, Google's
            KML uses 'WGS84'.
        lat : numeric
            latitude
        long : numeric
            longitude
        alt : numeric
            altitude
        alt_unit : string
            unit for the altitude (e.g. 'meter')

    """

    def __init__(self, geodetic_datum, lat, long, alt=None, alt_unit=None):
        """Initialize value for the Point object."""
        self.geodetic_datum = geodetic_datum
        self.lat = lat
        self.long = long
        self.alt = alt
        self.alt_unit = alt_unit


class Polygon(PhyloElement):
    """A polygon defined by a list of 'Points' (used by element 'Distribution').

    :param points: list of 3 or more points representing vertices.

    """

    def __init__(self, points=None):
        """Initialize value for the Polygon object."""
        self.points = points or []

    def __str__(self):
        """Return list of points as a string."""
        return "%s([%s])" % (self.__class__.__name__, ",\n".join(map(str, self.points)))


class Property(PhyloElement):
    """A typed and referenced property from an external resources.

    Can be attached to ``Phylogeny``, ``Clade``, and ``Annotation`` objects.

    :Parameters:
        value : string
            the value of the property
        ref : string
            reference to an external resource, e.g. "NOAA:depth"
        applies_to : string
            indicates the item to which a property applies to (e.g.  'node' for
            the parent node of a clade, 'parent_branch' for the parent branch of
            a clade, or just 'clade').
        datatype : string
            the type of a property; limited to xsd-datatypes
            (e.g. 'xsd:string', 'xsd:boolean', 'xsd:integer', 'xsd:decimal',
            'xsd:float', 'xsd:double', 'xsd:date', 'xsd:anyURI').
        unit : string (optional)
            the unit of the property, e.g. "METRIC:m"
        id_ref : Id (optional)
            allows to attached a property specifically to one element (on the
            xml-level)

    """

    re_ref = re.compile(r"[a-zA-Z0-9_]+:[a-zA-Z0-9_\.\-\s]+")
    ok_applies_to = {
        "phylogeny",
        "clade",
        "node",
        "annotation",
        "parent_branch",
        "other",
    }
    ok_datatype = {
        "xsd:string",
        "xsd:boolean",
        "xsd:decimal",
        "xsd:float",
        "xsd:double",
        "xsd:duration",
        "xsd:dateTime",
        "xsd:time",
        "xsd:date",
        "xsd:gYearMonth",
        "xsd:gYear",
        "xsd:gMonthDay",
        "xsd:gDay",
        "xsd:gMonth",
        "xsd:hexBinary",
        "xsd:base64Binary",
        "xsd:anyURI",
        "xsd:normalizedString",
        "xsd:token",
        "xsd:integer",
        "xsd:nonPositiveInteger",
        "xsd:negativeInteger",
        "xsd:long",
        "xsd:int",
        "xsd:short",
        "xsd:byte",
        "xsd:nonNegativeInteger",
        "xsd:unsignedLong",
        "xsd:unsignedInt",
        "xsd:unsignedShort",
        "xsd:unsignedByte",
        "xsd:positiveInteger",
    }

    def __init__(self, value, ref, applies_to, datatype, unit=None, id_ref=None):
        """Initialize value for the Property object."""
        _check_str(ref, self.re_ref.match)
        _check_str(applies_to, self.ok_applies_to.__contains__)
        _check_str(datatype, self.ok_datatype.__contains__)
        _check_str(unit, self.re_ref.match)
        self.unit = unit
        self.id_ref = id_ref
        self.value = value
        self.ref = ref
        self.applies_to = applies_to
        self.datatype = datatype


class ProteinDomain(PhyloElement):
    """Represents an individual domain in a domain architecture.

    The locations use 0-based indexing, as most Python objects including
    SeqFeature do, rather than the usual biological convention starting at 1.
    This means the start and end attributes can be used directly as slice
    indexes on Seq objects.

    :Parameters:
        start : non-negative integer
            start of the domain on the sequence, using 0-based indexing
        end : non-negative integer
            end of the domain on the sequence
        confidence : float
            can be used to store e.g. E-values
        id : string
            unique identifier/name

    """

    def __init__(self, value, start, end, confidence=None, id=None):
        """Initialize value for a ProteinDomain object."""
        self.value = value
        self.start = start
        self.end = end
        self.confidence = confidence
        self.id = id

    @classmethod
    def from_seqfeature(cls, feat):
        """Create ProteinDomain object from SeqFeature."""
        return ProteinDomain(
            feat.id,
            feat.location.start,
            feat.location.end,
            confidence=feat.qualifiers.get("confidence"),
        )

    def to_seqfeature(self):
        """Create a SeqFeature from the ProteinDomain Object."""
        feat = SeqFeature(location=SimpleLocation(self.start, self.end), id=self.value)
        try:
            confidence = self.confidence
        except AttributeError:
            pass
        else:
            feat.qualifiers["confidence"] = confidence
        return feat


class Reference(PhyloElement):
    """Literature reference for a clade.

    NB: Whenever possible, use the ``doi`` attribute instead of the free-text
    ``desc`` element.
    """

    re_doi = re.compile(r"[a-zA-Z0-9_\.]+/[a-zA-Z0-9_\.]+")

    def __init__(self, doi=None, desc=None):
        """Initialize elements of the Reference class object."""
        _check_str(doi, self.re_doi.match)
        self.doi = doi
        self.desc = desc


class Sequence(PhyloElement):
    """A molecular sequence (Protein, DNA, RNA) associated with a node.

    One intended use for ``id_ref`` is to link a sequence to a taxonomy (via the
    taxonomy's ``id_source``) in case of multiple sequences and taxonomies per
    node.

    :Parameters:
        type : {'dna', 'rna', 'protein'}
            type of molecule this sequence represents
        id_ref : string
            reference to another resource
        id_source : string
            source for the reference
        symbol : string
            short symbol of the sequence, e.g. 'ACTM' (max. 10 chars)
        accession : Accession
            accession code for this sequence.
        name : string
            full name of the sequence, e.g. 'muscle Actin'
        location
            location of a sequence on a genome/chromosome.
        mol_seq : MolSeq
            the molecular sequence itself
        uri : Uri
            link
        annotations : list of Annotation objects
            annotations on this sequence
        domain_architecture : DomainArchitecture
            protein domains on this sequence
        other : list of Other objects
            non-phyloXML elements

    """

    types = {"dna", "rna", "protein"}
    re_symbol = re.compile(r"\S{1,10}")

    def __init__(
        self,
        # Attributes
        type=None,
        id_ref=None,
        id_source=None,
        # Child nodes
        symbol=None,
        accession=None,
        name=None,
        location=None,
        mol_seq=None,
        uri=None,
        domain_architecture=None,
        # Collections
        annotations=None,
        other=None,
    ):
        """Initialize value for a Sequence object."""
        _check_str(type, self.types.__contains__)
        _check_str(symbol, self.re_symbol.match)
        self.type = type
        self.id_ref = id_ref
        self.id_source = id_source
        self.symbol = symbol
        self.accession = accession
        self.name = name
        self.location = location
        self.mol_seq = mol_seq
        self.uri = uri
        self.domain_architecture = domain_architecture
        self.annotations = annotations or []
        self.other = other or []

    @classmethod
    def from_seqrecord(cls, record, is_aligned=None):
        """Create a new PhyloXML Sequence from a SeqRecord object."""
        if is_aligned is None:
            is_aligned = "-" in record.seq
        params = {
            "accession": Accession(record.id, ""),
            "symbol": record.name,
            "name": record.description,
            "mol_seq": MolSeq(str(record.seq), is_aligned),
        }
        molecule_type = record.annotations.get("molecule_type")
        if molecule_type is not None:
            if "DNA" in molecule_type:
                params["type"] = "dna"
            elif "RNA" in molecule_type:
                params["type"] = "rna"
            elif "protein" in molecule_type:
                params["type"] = "protein"

        # Unpack record.annotations
        for key in ("id_ref", "id_source", "location"):
            if key in record.annotations:
                params[key] = record.annotations[key]
        if isinstance(record.annotations.get("uri"), dict):
            params["uri"] = Uri(**record.annotations["uri"])
        # Build a Sequence.annotation object
        if record.annotations.get("annotations"):
            params["annotations"] = []
            for annot in record.annotations["annotations"]:
                ann_args = {}
                for key in ("ref", "source", "evidence", "type", "desc"):
                    if key in annot:
                        ann_args[key] = annot[key]
                if isinstance(annot.get("confidence"), list):
                    ann_args["confidence"] = Confidence(*annot["confidence"])
                if isinstance(annot.get("properties"), list):
                    ann_args["properties"] = [
                        Property(**prop)
                        for prop in annot["properties"]
                        if isinstance(prop, dict)
                    ]
                params["annotations"].append(Annotation(**ann_args))

        # Unpack record.features
        if record.features:
            params["domain_architecture"] = DomainArchitecture(
                length=len(record.seq),
                domains=[
                    ProteinDomain.from_seqfeature(feat) for feat in record.features
                ],
            )

        return Sequence(**params)

    def to_seqrecord(self):
        """Create a SeqRecord object from this Sequence instance.

        The seqrecord.annotations dictionary is packed like so::

            { # Sequence attributes with no SeqRecord equivalent:
              'id_ref': self.id_ref,
              'id_source': self.id_source,
              'location': self.location,
              'uri': { 'value': self.uri.value,
                              'desc': self.uri.desc,
                              'type': self.uri.type },
              # Sequence.annotations attribute (list of Annotations)
              'annotations': [{'ref': ann.ref,
                               'source': ann.source,
                               'evidence': ann.evidence,
                               'type': ann.type,
                               'confidence': [ann.confidence.value,
                                              ann.confidence.type],
                               'properties': [{'value': prop.value,
                                                'ref': prop.ref,
                                                'applies_to': prop.applies_to,
                                                'datatype': prop.datatype,
                                                'unit': prop.unit,
                                                'id_ref': prop.id_ref}
                                               for prop in ann.properties],
                              } for ann in self.annotations],
            }

        """

        def clean_dict(dct):
            """Remove None-valued items from a dictionary."""
            return {key: val for key, val in dct.items() if val is not None}

        seqrec = SeqRecord(
            Seq(self.mol_seq.value),
            **clean_dict(
                {
                    "id": str(self.accession),
                    "name": self.symbol,
                    "description": self.name,
                    # 'dbxrefs': None,
                }
            ),
        )
        if self.domain_architecture:
            seqrec.features = [
                dom.to_seqfeature() for dom in self.domain_architecture.domains
            ]
        # Sequence attributes with no SeqRecord equivalent
        if self.type == "dna":
            molecule_type = "DNA"
        elif self.type == "rna":
            molecule_type = "RNA"
        elif self.type == "protein":
            molecule_type = "protein"
        else:
            molecule_type = None
        seqrec.annotations = clean_dict(
            {
                "id_ref": self.id_ref,
                "id_source": self.id_source,
                "location": self.location,
                "uri": self.uri
                and clean_dict(
                    {
                        "value": self.uri.value,
                        "desc": self.uri.desc,
                        "type": self.uri.type,
                    }
                ),
                "molecule_type": molecule_type,
                "annotations": self.annotations
                and [
                    clean_dict(
                        {
                            "ref": ann.ref,
                            "source": ann.source,
                            "evidence": ann.evidence,
                            "type": ann.type,
                            "confidence": ann.confidence
                            and [ann.confidence.value, ann.confidence.type],
                            "properties": [
                                clean_dict(
                                    {
                                        "value": prop.value,
                                        "ref": prop.ref,
                                        "applies_to": prop.applies_to,
                                        "datatype": prop.datatype,
                                        "unit": prop.unit,
                                        "id_ref": prop.id_ref,
                                    }
                                )
                                for prop in ann.properties
                            ],
                        }
                    )
                    for ann in self.annotations
                ],
            }
        )
        return seqrec


class SequenceRelation(PhyloElement):
    """Express a typed relationship between two sequences.

    For example, this could be used to describe an orthology (in which case
    attribute 'type' is 'orthology').

    :Parameters:
        id_ref_0 : Id
            first sequence reference identifier
        id_ref_1 : Id
            second sequence reference identifier
        distance : float
            distance between the two sequences
        type : restricted string
            describe the type of relationship
        confidence : Confidence
            confidence value for this relation

    """

    ok_type = {
        "orthology",
        "one_to_one_orthology",
        "super_orthology",
        "paralogy",
        "ultra_paralogy",
        "xenology",
        "unknown",
        "other",
    }

    def __init__(self, type, id_ref_0, id_ref_1, distance=None, confidence=None):
        """Initialize the class."""
        _check_str(type, self.ok_type.__contains__)
        self.distance = distance
        self.type = type
        self.id_ref_0 = id_ref_0
        self.id_ref_1 = id_ref_1
        self.confidence = confidence


class Taxonomy(PhyloElement):
    """Describe taxonomic information for a clade.

    :Parameters:
        id_source : Id
            link other elements to a taxonomy (on the XML level)
        id : Id
            unique identifier of a taxon, e.g. Id('6500',
            provider='ncbi_taxonomy') for the California sea hare
        code : restricted string
            store UniProt/Swiss-Prot style organism codes, e.g. 'APLCA' for the
            California sea hare 'Aplysia californica'
        scientific_name : string
            the standard scientific name for this organism, e.g. 'Aplysia
            californica' for the California sea hare
        authority : string
            keep the authority, such as 'J. G. Cooper, 1863', associated with
            the 'scientific_name'
        common_names : list of strings
            common names for this organism
        synonyms : list of strings
            synonyms for this taxon?
        rank : restricted string
            taxonomic rank
        uri : Uri
            link
        other : list of Other objects
            non-phyloXML elements

    """

    re_code = re.compile(r"[a-zA-Z0-9_]{2,10}")
    ok_rank = {
        "domain",
        "kingdom",
        "subkingdom",
        "branch",
        "infrakingdom",
        "superphylum",
        "phylum",
        "subphylum",
        "infraphylum",
        "microphylum",
        "superdivision",
        "division",
        "subdivision",
        "infradivision",
        "superclass",
        "class",
        "subclass",
        "infraclass",
        "superlegion",
        "legion",
        "sublegion",
        "infralegion",
        "supercohort",
        "cohort",
        "subcohort",
        "infracohort",
        "superorder",
        "order",
        "suborder",
        "superfamily",
        "family",
        "subfamily",
        "supertribe",
        "tribe",
        "subtribe",
        "infratribe",
        "genus",
        "subgenus",
        "superspecies",
        "species",
        "subspecies",
        "variety",
        "subvariety",
        "form",
        "subform",
        "cultivar",
        "unknown",
        "other",
    }

    def __init__(
        self,
        # Attributes
        id_source=None,
        # Child nodes
        id=None,
        code=None,
        scientific_name=None,
        authority=None,
        rank=None,
        uri=None,
        # Collections
        common_names=None,
        synonyms=None,
        other=None,
    ):
        """Initialize the class."""
        _check_str(code, self.re_code.match)
        _check_str(rank, self.ok_rank.__contains__)
        self.id_source = id_source
        self.id = id
        self.code = code
        self.scientific_name = scientific_name
        self.authority = authority
        self.rank = rank
        self.uri = uri
        self.common_names = common_names or []
        self.synonyms = synonyms or []
        self.other = other or []

    def __str__(self):
        """Show the class name and an identifying attribute."""
        if self.code is not None:
            return self.code
        if self.scientific_name is not None:
            return self.scientific_name
        if self.rank is not None:
            return self.rank
        if self.id is not None:
            return str(self.id)
        return self.__class__.__name__


class Uri(PhyloElement):
    """A uniform resource identifier.

    In general, this is expected to be an URL (for example, to link to an image
    on a website, in which case the ``type`` attribute might be 'image' and
    ``desc`` might be 'image of a California sea hare').
    """

    def __init__(self, value, desc=None, type=None):
        """Initialize the class."""
        self.value = value
        self.desc = desc
        self.type = type

    def __str__(self):
        """Return string representation of Uri."""
        if self.value:
            return self.value
        return repr(self)
