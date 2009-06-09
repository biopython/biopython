# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Objects instantiated from elements in a parsed PhyloXML file.

"""

import warnings

try:
    from xml.etree import cElementTree as ElementTree
except ImportError:
    try:
        from xml.etree import ElementTree as ElementTree
    except ImportError:
        # Python 2.4 -- check for 3rd-party implementations
        try:
            from lxml.etree import ElementTree
        except ImportError:
            try:
                import cElementTree as ElementTree
            except ImportError:
                try:
                    from elementtree import ElementTree
                except ImportError:
                    from Bio import MissingExternalDependencyError
                    raise MissingExternalDependencyError(
                            "No ElementTree module was found. " \
                            "Use Python 2.5+, lxml or elementtree if you " \
                            "want to use Bio.PhyloXML.")


# Lookup table used to instantiate elements by XML tag
tags_to_classes = {
        ## special cases for parsing
#         'phylogeny':    Phylogeny,
#         'phyloxml':     Phyloxml,
#         'clade':        Clade,
        ## no special handling
#         'absent':       BinaryCharacterList,
#         'accession':    Accession,
#         'alt':          float,  # decimal
#         'annotation':   Annotation,
#         'bc':           str,
#         'binary_characters': BinaryCharacters,
#         'blue':         int, # unsignedByte
#         'branch_length': float,  # double
#         'clade_relation': CladeRelation,
#         'code':         TaxonomyCode,
#         'color':        BranchColor,
#         'common_name':  str,
#         'confidence':   Confidence,
#         'date':         Date,
#         'desc':         str,
#         'description':  str,
#         'distribution': Distribution,
#         'domain':       ProteinDomain,
#         'domain_architecture': DomainArchitecture,
#         'duplications': int, # nonNegativeInteger
#         'events':       Events,
#         'gained':       BinaryCharacterList,
#         'green':        int, # unsignedByte
#         'id':           Id,
#         'lat':          float,  # decimal
#         'location':     str,
#         'long':         float,  # decimal
#         'losses':       int, # nonNegativeInteger
#         'lost':         BinaryCharacterList,
#         'mol_seq':      MolSeq,
#         'name':         str,
#         'node_id':      Id,
#         'point':        Point,
#         'polygon':      Polygon,
#         'present':      BinaryCharacterList,
#         'property':     Property,
#         'rank':         Rank,
#         'red':          int, # unsignedByte
#         'reference':    Reference,
#         'scientific_name': str,
#         'sequence':     Sequence,
#         'sequence_relation': SequenceRelation,
#         'speciations':  int, # nonNegativeInteger
#         'symbol':       SequenceSymbol,
#         'taxonomy':     Taxonomy,
#         'type':         EventType,
#         'uri':          Uri,
#         'value':        float, # decimal
#         'width':        float, # double
        }

# Functions I wish ElementTree had

def local(tag):
    if tag[0] == '{':
        return tag.rsplit('}', 1)[1]
    return tag


def get_elem_text(elem, tag, default=None):
    val = elem.find(tag)
    if val is None:
        return default
    return elem.text


# Debugging helper

def _dump_tags(handle):
    """Extract tags from an XML document and print them to standard output.
    
    This function is meant for testing and debugging only.
    """
    events = ('start', 'end')
    for event, elem in ElementTree.iterparse(handle, events=events):
        if event == 'start':
            print local(elem.tag)
        else:
            elem.clear()


def read(handle):
    """Parse a phyloXML file or stream and build a tree of Biopython objects.

    The children of the root node are phylogenies and possibly other arbitrary
    (non-phyloXML) objects.

    To minimize memory use, the tree of ElementTree parsing events is cleared
    after completing each phylogeny, clade, and top-level 'other' element.
    Elements below the clade level are kept in memory until parsing of the
    current clade is finished -- this shouldn't be a problem because clade is
    the main recursive element, and non-clade nodes below this level are of
    bounded size.
    """
    # get an iterable context for XML parsing events
    context = iter(ElementTree.iterparse(handle, events=('start', 'end')))
    event, root = context.next()
    phyloxml = Phyloxml(root.attrib)
    for event, elem in context:
        if event == 'start' and local(elem.tag) == 'phylogeny':
            phylogeny = _parse_phylogeny(elem, context)
            # warnings.warn('Built a phylogeny %s, contents:' % repr(phylogeny))
            # warnings.warn(repr(phylogeny.__dict__))
            phyloxml.phylogenies.append(phylogeny)
            continue
        if event == 'end' and local(elem.tag) != 'phyloxml':
            # Deal with items not specified by phyloXML
            otr = Other.from_element(elem)
            phyloxml.other.append(otr)
            root.clear()
    return phyloxml


def _parse_phylogeny(parent, context):
    """Parse a single phylogeny within the phyloXML tree.

    Recursively builds a phylogenetic tree with help from parse_clade, then
    clears the XML event history for the phylogeny element and returns control
    to the top-level parsing function.
    """
    phylogeny = Phylogeny(parent.attrib)
    for event, elem in context:
        tag = local(elem.tag)
        if event == 'start' and tag == 'clade':
            assert phylogeny.clade is None, \
                    "Phylogeny object should only have 1 clade"
            phylogeny.clade = _parse_clade(elem, context)
            continue
        if event == 'end':
            if tag == 'phylogeny':
                parent.clear()
                break
            # Handle the other non-recursive children
            elif tag == 'id':
                phylogeny.id = Id(text=elem.text)
            elif tag == 'date':
                phylogeny.date = Date(text=elem.text)
            elif tag == 'clade_relation':
                phylogeny.clade_relation = CladeRelation(text=elem.text)
            elif tag == 'sequence_relation':
                phylogeny.sequence_relation = SequenceRelation(text=elem.text)
            # Simple types
            if tag == 'name': 
                phylogeny.name = str(elem.text)
            elif tag == 'description':
                phylogeny.description = str(elem.text)
            # Collections
            elif tag == 'confidence':
                phylogeny.confidences.append(Confidence.from_element(elem))
            elif tag == 'property':
                phylogeny.properties = [Property(text=elem.text)]
            # Unknown tags
            else:
                phylogeny.other.append(Other.from_element(elem))
            elem.clear()
    return phylogeny


def _parse_clade(parent, context):
    clade = Clade(parent.attrib)
    for event, elem in context:
        tag = local(elem.tag)
        if event == 'start' and tag == 'clade':
            subclade = _parse_clade(elem, context)
            clade.clades.append(subclade)
            continue
        if event == 'end':
            if tag == 'clade':
                parent.clear()
                break
            # Handle the other non-recursive children
            elif tag == 'branch_length':
                # NB: possible collision with the attribute
                if hasattr(clade, 'branch_length'):
                    warnings.warn('Attribute branch_length was already set for '
                                  'this Clade; overwriting the previous value.')
                clade.branch_length = elem.text
            elif tag == 'color':
                clade.color == BranchColor.from_element(elem)
            elif tag == 'node_id':
                clade.node_id == Id(text=elem.text)
            elif tag == 'events':
                clade.events == Events(text=elem.text)
            elif tag == 'binary_characters':
                clade.binary_characters == BinaryCharacters(text=elem.text)
            elif tag == 'date':
                clade.date == Date(text=elem.text)
            # Simple types
            if tag == 'name': 
                clade.name = str(elem.text)
            elif tag == 'width':
                clade.width == float(elem.text)
            # Collections
            elif tag == 'confidence':
                clade.confidences.append(Confidence(elem))
            elif tag == 'taxonomy':
                clade.taxonomies.append(Taxonomy.from_element(elem))
            elif tag == 'sequence':
                clade.sequences.append(Sequence(elem))
            elif tag == 'distributions':
                clade.distributions.append(Distribution(elem))
            elif tag == 'reference':
                clade.references.append(Reference(elem))
            elif tag == 'property':
                clade.properties.append(Property(elem))
            # Unknown tags
            else:
                clade.other.append(Other.from_element(elem))
            elem.clear()
    return clade


# ---------------------------------------------------------------------
# Classes instantiated from phyloXML nodes

# XXX maybe not needed
class PhyloElement(object):
    """Base class for all PhyloXML objects."""
    def __init__(self, attrib=None, text=None, **kwargs):
        if attrib is not None:
            # self._attrib = attrib
            self.__dict__.update(attrib)
        # if text is not None:
        #     self._text = text
        self.__dict__.update(kwargs)

    @classmethod
    def from_element(cls, elem):
        raise NotImplementedError("This method should be implemented by " \
                                  "the derived class %s." % cls)


class Other(PhyloElement):
    """Container for non-phyloXML elements in the tree."""
    # ENH: assert that the tag namespace is not phyloxml's
    def __init__(self, tag, attributes=None, text=None, children=[]):
        self.tag = tag
        self.attributes = attributes
        self.text = text
        self.children = children

    def __str__(self):
        return '<Other %s at %s>' % (self.tag, hex(id(self)))

    @classmethod
    def from_element(cls, elem):
        obj = Other(elem.tag, elem.attrib, elem.text)
        for child in elem:
            obj.children.append(Other.from_element(child))
        return obj


# Core elements

class Phyloxml(PhyloElement):
    """Root node of the PhyloXML document.

    Contains an arbitrary number of Phylogeny elements, possibly followed by
    elements from other namespaces.

    Attributes:
        (namespace definitions)

    Children:
        phylogenies []
        other []
    """
    def __init__(self, attributes, phylogenies=[]):
        PhyloElement.__init__(self, attributes, phylogenies=phylogenies)
        self.other = []

    def __iter__(self):
        """Iterate through the phylogenetic trees in this object."""
        return iter(self.phylogenies)

    def __len__(self):
        """Number of phylogenetic trees in this object."""
        return len(self.phylogenies)


class Phylogeny(PhyloElement):
    """A phylogenetic tree.

    Attributes:
        rooted
        rerootable
        branch_length_unit
        type

    Children:
        name
        id
        description
        date
        confidences []
        clade
        clade_relations []
        sequence_relations []
        properties []
        other []
    """
    def __init__(self, attributes, **kwargs):
        PhyloElement.__init__(self, attributes, **kwargs)
        # Single values
        for attr in (
                # Node attributes
                'rooted', 'rerootable', 'branch_length_unit', 'type',
                # Child nodes
                'name', 'id', 'description', 'date', 'clade',
                ):
            if not hasattr(self, attr):
                setattr(self, attr, None)
        # Lists
        for attr in (
                'confidences', 'clade_relations', 'sequence_relations',
                'properties', 'other',
                ):
            if not hasattr(self, attr):
                setattr(self, attr, [])

    # def __iter__(self):
    #     """Iterate through the clades (branches) within this phylogeny."""
    #     return iter(self.clade)

    # def __len__(self):
    #     """Number of clades directly under this element."""
        # XXX should only have 1 clade
        # ENH: count all branches within this tree?

    # From Bioperl's Bio::Tree::TreeI

    def get_leaf_nodes(self):
        """Request the taxa (leaves of the tree)."""
        raise NotImplementedError

    def get_root_node(self):
        """Get the root node of this tree."""
        return self

    def total_branch_length(self):
        """Get the total length of this tree (sum of all branch lengths)."""
        raise NotImplementedError

    # From Bioperl's Bio::Tree::TreeFunctionsI

    # find_node -- by name or other standard field
    # remove_node
    # get_lca (lowest common ancestor)
    # distance (between 2 nodes, specified however)
    # is_monophyletic
    # is_paraphyletic
    # reroot


class Clade(PhyloElement):
    """Describes a branch of the current phylogenetic tree.

    Used recursively, describes the topology of a phylogenetic tree.

    The parent branch length of a clade can be described either with the
    'branch_length' element or the 'branch_length' attribute (it is not
    recommended to use both at the same time, though). Usage of the
    'branch_length' attribute allows for a less verbose description.

    Element 'confidence' is used to indicate the support for a clade/parent
    branch.

    Element 'events' is used to describe such events as gene-duplications at the
    root node/parent branch of a clade.

    Element 'width' is the branch width for this clade (including parent
    branch). Both 'color' and 'width' elements apply for the whole clade unless
    overwritten in-sub clades.

    Attribute 'id_source' is used to link other elements to a clade (on the
    xml-level).

    Attributes:
        branch_length
        id_source

    Children:
        name
        branch_length -- equivalent to the attribute
        confidences []
        width
        color
        node_id
        taxonomies []
        sequences []
        events
        binary_characters
        distributions []
        date
        references []
        properties []
        clades [] -- recursive
        other []
    """
    def __init__(self, attributes, **kwargs):
        PhyloElement.__init__(self, attributes, **kwargs)
        # Single values
        for attr in (
                # Attributes
                'branch_length', 'id_source',
                # Child nodes
                'name', 'width', 'color', 'node_id', 
                'events', 'binary_characters', 'date',
                ):
            if not hasattr(self, attr):
                setattr(self, attr, None)
        # Collections
        for attr in (
                'confidences', 'taxonomies', 'sequences', 'distributions',
                'references', 'properties', 'clades', 'other',
                ):
            if not hasattr(self, attr):
                setattr(self, attr, [])

    def __iter__(self):
        """Iterate through the clades (sub-nodes) within this clade."""
        return iter(self.clades)

    def __len__(self):
        """Number of clades directy under this element."""
        return len(self.clades)


# Complex types

class Accession(PhyloElement):
    """
    """

class Annotation(PhyloElement):
    """
    """

class BinaryCharacterList(PhyloElement):
    """
    """

class BinaryCharacters(PhyloElement):
    """
    """

class BranchColor(PhyloElement):
    """Indicates the color of a clade when rendered graphically.

    The color applies to the whole clade unless overwritten by the color(s) of
    sub-clades.
    """
    def __init__(self, red, green, blue):
        assert isinstance(red, int)
        assert isinstance(green, int)
        assert isinstance(blue, int)
        self.red = red
        self.green = green
        self.blue = blue

    @classmethod
    def from_element(cls, elem):
        red, green, blue = (int(elem.find(color).text) for color in
                            ('red', 'green', 'blue'))
        return BranchColor(red, green, blue)


class CladeRelation(PhyloElement):
    """
    """
    def __init__(self, text):
        PhyloElement.__init__(self, text=text)

class Confidence(PhyloElement):
    """
    """
    def __init__(self, text):
        PhyloElement.__init__(self, text=text)

class Date(PhyloElement):
    """
    """
    def __init__(self, text):
        PhyloElement.__init__(self, text=text)

class Distribution(PhyloElement):
    """
    """

class DomainArchitecture(PhyloElement):
    """
    """

class Events(PhyloElement):
    """
    """

class Id(PhyloElement):
    """
    """
    def __init__(self, text):
        PhyloElement.__init__(self, text=text)

class Point(PhyloElement):
    """
    """

class Polygon(PhyloElement):
    """
    """

class Property(PhyloElement):
    """
    """
    def __init__(self, text):
        PhyloElement.__init__(self, text=text)

class ProteinDomain(PhyloElement):
    """
    """

class Reference(PhyloElement):
    """
    """

class Sequence(PhyloElement):
    """
    """

class SequenceRelation(PhyloElement):
    """
    """
    def __init__(self, text):
        PhyloElement.__init__(self, text=text)

class Taxonomy(PhyloElement):
    """Describe taxonomic information for a clade.

    Element 'code' is intended to store UniProt/Swiss-Prot style organism codes
    (e.g. 'APLCA' for the California sea hare 'Aplysia californica').

    Element 'id' is used for a unique identifier of a taxon (for example '6500'
    with 'ncbi_taxonomy' as 'type' for the California sea hare).

    Attributes:
        type
        id_source -- link other elements to a taxonomy (on the XML level)

    Children:
        id
        code
        scientific_name
        common_names []
        rank
        uri
        other []
    """
    def __init__(self, attributes, **kwargs):
        PhyloElement.__init__(self, attributes, **kwargs)
        # Single values
        for attr in (
                # Attributes
                'type', 'id_source',
                # Child nodes
                'id', 'code', 'scientific_name', 'rank', 'uri',
                ):
            if not hasattr(self, attr):
                setattr(self, attr, None)
        # Collections
        for attr in (
                'common_names', 'other',
                ):
            if not hasattr(self, attr):
                setattr(self, attr, [])

    @classmethod
    def from_element(cls, elem):
        return Taxonomy(elem.attrib, 
                id=Id(get_elem_text(elem, 'id')),
                code=TaxonomyCode(get_elem_text(elem, 'code')),
                scientific_name=get_elem_text(elem, ('scientific_name')),
                common_names=[e.text for e in elem.findall('common_name')],
                rank=Rank(get_elem_text(elem, 'rank')),
                uri=Uri(get_elem_text(elem, 'uri')),
                )


class Uri(PhyloElement):
    """
    """


# Simple types

class AppliesTo(PhyloElement):
    pass

class Doi(PhyloElement):
    pass

class EventType(PhyloElement):
    pass

class MolSeq(PhyloElement):
    pass

class PropertyDataType(PhyloElement):
    pass

class Rank(PhyloElement):
    pass

class SequenceRelationType(PhyloElement):
    pass

class SequenceSymbol(PhyloElement):
    pass

class SequenceType(PhyloElement):
    pass

class TaxonomyCode(PhyloElement):
    pass

class id_ref(PhyloElement):
    pass

class id_source(PhyloElement):
    pass

class ref(PhyloElement):
    pass

