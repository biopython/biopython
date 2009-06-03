# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Objects instantiated from elements in a parsed PhyloXML file.

"""

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
        'absent':       BinaryCharacterList,
        'accession':    Accession,
        'alt':          float,  # decimal
        'annotation':   Annotation,
        'bc':           Token,
        'binary_characters': BinaryCharacters,
        'blue':         int, # unsignedByte
        'branch_length': float,  # double
        'clade':        Clade,
        'clade_relation': CladeRelation,
        'code':         TaxonomyCode,
        'color':        BranchColor,
        'common_name':  Token,
        'confidence':   Confidence,
        'date':         Date,
        'desc':         Token,
        'description':  Token,
        'distribution': Distribution,
        'domain':       ProteinDomain,
        'domain_architecture': DomainArchitecture,
        'duplications': int, # nonNegativeInteger
        'events':       Events,
        'gained':       BinaryCharacterList,
        'green':        int, # unsignedByte
        'id':           Id,
        'lat':          float,  # decimal
        'location':     Token,
        'long':         float,  # decimal
        'losses':       int, # nonNegativeInteger
        'lost':         BinaryCharacterList,
        'mol_seq':      MolSeq,
        'name':         Token,
        'node_id':      Id,
        'phylogeny':    Phylogeny,
        'phyloxml':     Phyloxml,
        'point':        Point,
        'polygon':      Polygon,
        'present':      BinaryCharacterList,
        'property':     Property,
        'rank':         Rank,
        'red':          int, # unsignedByte
        'reference':    Reference,
        'scientific_name': Token,
        'sequence':     Sequence,
        'sequence_relation': SequenceRelation,
        'speciations':  int, # nonNegativeInteger
        'symbol':       SequenceSymbol,
        'taxonomy':     Taxonomy,
        'type':         EventType,
        'uri':          Uri,
        'value':        float, # decimal
        'width':        float, # double
        }



def _dump_tags(source):
    """Extract tags from an XML document and print them to standard output.
    
    This function is meant for testing and debugging only.
    """
    events = ('start', 'end')
    for event, elem in ElementTree.iterparse(source, events=events):
        if event == 'start':
            print elem.tag
        else:
            elem.clear()


def read(source):
    """Parse a phyloXML file or stream and build a tree of Biopython objects.

    The children of the root node are phylogenies and possibly other arbitrary
    (non-phyloXML) objects.
    """
    events = ('start', 'end')
    for event, elem in ElementTree.iterparse(source, events=events):
        if event == 'start':
            pass
        else:
            # dispatch by elem.tag
            obj = get_phylo_obj(elem.tag, elem.attrib, elem.text)
            # instantiate a node


def get_phylo_obj(tag, attrib, text):
    # look up the class by tag
    constructor = tags_to_classes.get(tag, Other)
    return constructor(attrib, text, None)


# ---------------------------------------------------------------------
# Classes instantiated from phyloXML nodes

class PhyloElement(object):
    """Base class for all PhyloXML objects."""
    def __init__(self, attrib, text, children):
        self._attrib = attrib
        self._text = text
        self._children = children
        # Munge each of these into properties
        self._expose

    def _expose(self):
        """Produce a useful interface for protected data.

        Overridden by classes that inherit from this.
        """
        pass


class Other(PhyloElement):
    """Container for non-phyloXML elements in the tree."""


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

