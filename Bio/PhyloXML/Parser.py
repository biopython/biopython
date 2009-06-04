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
#         'bc':           Token,
#         'binary_characters': BinaryCharacters,
#         'blue':         int, # unsignedByte
#         'branch_length': float,  # double
#         'clade_relation': CladeRelation,
#         'code':         TaxonomyCode,
#         'color':        BranchColor,
#         'common_name':  Token,
#         'confidence':   Confidence,
#         'date':         Date,
#         'desc':         Token,
#         'description':  Token,
#         'distribution': Distribution,
#         'domain':       ProteinDomain,
#         'domain_architecture': DomainArchitecture,
#         'duplications': int, # nonNegativeInteger
#         'events':       Events,
#         'gained':       BinaryCharacterList,
#         'green':        int, # unsignedByte
#         'id':           Id,
#         'lat':          float,  # decimal
#         'location':     Token,
#         'long':         float,  # decimal
#         'losses':       int, # nonNegativeInteger
#         'lost':         BinaryCharacterList,
#         'mol_seq':      MolSeq,
#         'name':         Token,
#         'node_id':      Id,
#         'point':        Point,
#         'polygon':      Polygon,
#         'present':      BinaryCharacterList,
#         'property':     Property,
#         'rank':         Rank,
#         'red':          int, # unsignedByte
#         'reference':    Reference,
#         'scientific_name': Token,
#         'seq':          Sequence,
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

def local(tag):
    if tag[0] == '{':
        return tag.rsplit('}', 1)[1]
    return tag


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
    phyloxml = Phyloxml(attrib=root.attrib)
    for event, elem in context:
        if event == 'start' and local(elem.tag) == 'phylogeny':
            phylogeny = parse_phylogeny(elem, context)
            # warnings.warn('Built a phylogeny %s, contents:' % repr(phylogeny))
            # warnings.warn(repr(phylogeny.__dict__))
            phyloxml.phylogenies.append(phylogeny)
            continue
        if event == 'end':
            # Deal with Other items
            root.clear()
            pass
    return phyloxml


def parse_phylogeny(parent, context):
    """Parse a single phylogeny within the phyloXML tree.

    Recursively builds a phylogenetic tree with help from parse_clade, then
    clears the XML event history for the phylogeny element and returns control
    to the top-level parsing function.
    """
    phylogeny = Phylogeny(attrib=parent.attrib)
    for event, elem in context:
        if event == 'start' and local(elem.tag) == 'clade':
            clade = parse_clade(elem, context)
            phylogeny.clades.append(clade)
            continue
        if event == 'end':
            if local(elem.tag) == 'phylogeny':
                parent.clear()
                break
            # Handle the other non-recursive children
            if local(elem.tag) == 'name': 
                phylogeny.name = Token('name', elem.text)
                elem.clear()
            elif local(elem.tag) == 'id':
                phylogeny.id = Id(text=elem.text)
                elem.clear()
            elif local(elem.tag) == 'description':
                phylogeny.description = Token('description', elem.text)
                elem.clear()
            elif local(elem.tag) == 'date':
                phylogeny.date = Date(text=elem.text)
                elem.clear()
            elif local(elem.tag) == 'confidence':
                phylogeny.confidence = Confidence(text=elem.text)
                elem.clear()
            elif local(elem.tag) == 'clade_relation':
                phylogeny.clade_relation = CladeRelation(text=elem.text)
                elem.clear()
            elif local(elem.tag) == 'sequence_relation':
                phylogeny.sequence_relation = SequenceRelation(text=elem.text)
                elem.clear()
            elif local(elem.tag) == 'property':
                phylogeny.properties = [Property(text=elem.text)]
                elem.clear()
            else:
                # Unknown tag
                if hasattr(phylogeny, 'other'):
                    phylogeny.other.append(Other(elem.attrib, elem.text))
                else:
                    phylogeny.other = [Other(elem.attrib, elem.text)]
                elem.clear()
    return phylogeny


def parse_clade(parent, context):
    clade = Clade(attrib=parent.attrib)
    # clade_depth = 0
    # node = tags_to_classes.get(local(parent.tag), Other)(parent)
    for event, elem in context:
        if event == 'start' and local(elem.tag) == 'clade':
            # clade_depth += 1
            subclade = parse_clade(elem, context)
            clade.clades.append(subclade)
            continue
        if event == 'end' and local(elem.tag) == 'clade':
            # clade_depth -= 1
            parent.clear()
            break
    return clade


def parse_other(parent, context):
    node = Other(attrib=parent.attrib, text=parent.text)
    tag_depth = 0
    for event, elem in context:
        if event == 'start' and local(elem.tag) == local(parent.tag):
            tag_depth += 1
            continue
        if local(elem.tag) == local(parent.tag):
            tag_depth -= 1
            if tag_depth == 0:
                parent.clear()
                break
        else:
            constructor = tags_to_classes.get(local(elem.tag), Other)
            if constructor in (int, float, Token):
                obj = constructor(local(elem.tag), elem.text)
            else:
                obj = constructor(elem.attrib, elem.text)
            node.add_child(obj)
            elem.clear()
    return node

# ---------------------------------------------------------------------
# Classes instantiated from phyloXML nodes

class PhyloElement(object):
    """Base class for all PhyloXML objects."""
    def __init__(self, attrib=None, text=None):
        if attrib is not None:
            self._attrib = attrib
        if text is not None:
            self._text = text
        # self._children = []
        # Munge each of these into properties
        self._expose()

    def _expose(self):
        """Produce a useful interface for protected data.

        Overridden by classes that inherit from this.
        """
        pass


class Other(PhyloElement):
    """Container for non-phyloXML elements in the tree."""
    def add_child(self, tag, obj):
        if tag in self.__dict__:
            if isinstance(getattr(self, tag), list):
                getattr(self, tag).append(obj)
            else:
                setattr(self, tag, [getattr(self, tag), obj])
        else:
            setattr(self, tag, obj)


class Token(object):
    """Simple class for nodes containing only a string."""
    def __init__(self, tag, text):
        self.tag = tag
        self.text = text


# Core elements

class Phyloxml(PhyloElement):
    """Root node of the PhyloXML document.
    
    Contains an arbitrary number of Phylogeny elements, possibly followed by
    elements from other namespaces.
    """
    def _expose(self):
        self.phylogenies = []
        self.other = []


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
        confidence
        clade
        clade_relation
        sequence_relation
        property
        [other]
    """
    def _expose(self):
        self.__dict__.update(self._attrib)
        self.clades = []

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
        branch_length   (equivalent to the attribute)
        confidence
        width
        color
        node_id
        taxonomy
        sequence
        events
        binary_characters
        distribution
        date
        reference
        property
        clade   (recursive)
    """
    def _expose(self):
        self.__dict__.update(self._attrib)
        self.clades = []


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
    """
    """

class CladeRelation(PhyloElement):
    """
    """

class Confidence(PhyloElement):
    """
    """

class Date(PhyloElement):
    """
    """

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

class Point(PhyloElement):
    """
    """

class Polygon(PhyloElement):
    """
    """

class Property(PhyloElement):
    """
    """

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

class Taxonomy(PhyloElement):
    """
    """

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

