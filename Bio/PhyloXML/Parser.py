# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""PhyloXML reader/parser and associated functions.

Instantiates Tree elements from a parsed PhyloXML file.
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

import Tree
from Exceptions import PhyloXMLError, PhyloXMLWarning

NAMESPACES = {
        'phy':  'http://www.phyloxml.org',
        'xml':  'http://www.w3.org/XML/1998/namespace',
        'xs':   'http://www.w3.org/2001/XMLSchema',
        }

# ---------------------------------------------------------
# Functions I wish ElementTree had

def local(tag):
    if tag[0] is '{':
        return tag[tag.index('}')+1:]
    return tag

def namespace(tag):
    try:
        if tag[0] is '{':
            return tag[1:tag.index('}')]
    finally:
        return ''

def split_namespace(tag):
    parts = tag.split('}')
    if len(parts) == 1:
        return ('', parts[0])
    return (parts[0][1:], parts[1])

def get_child_as(parent, tag, construct):
    child = parent.find("{%s}%s" % (NAMESPACES['phy'], tag))
    if child is not None:
        return construct(child)

def get_child_text(parent, tag, construct=unicode):
    child = parent.find("{%s}%s" % (NAMESPACES['phy'], tag))
    if child is not None:
        return child.text and construct(child.text.strip()) or None

def get_children_as(parent, tag, construct):
    return [construct(child) for child in 
            parent.findall("{%s}%s" % (NAMESPACES['phy'], tag))]

def get_children_text(parent, tag, construct=unicode):
    return [construct(child.text.strip()) for child in 
            parent.findall("{%s}%s" % (NAMESPACES['phy'], tag))
            if child.text]

# ---------------------------------------------------------
# Utilities

import re

def check_str(text, regexp):
    """Compare a string to a regexp, and warn if there's no match."""
    if text is not None and re.match(regexp, text) is None:
        warnings.warn("String %s doesn't match regexp %s"
                        % (text, regexp), PhyloXMLWarning)
    return text

def str2bool(text):
    if text == 'true':
        return True
    if text == 'false':
        return False
    raise ValueError('String could not be converted to boolean: ' + text)

def dict_str2bool(dct, keys):
    out = dct.copy()
    for key in keys:
        if key in out:
            out[key] = str2bool(out[key])
    return out


# ---------------------------------------------------------

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
    phyloxml = Tree.Phyloxml(dict((local(key), val)
                             for key, val in root.items()))
    other_depth = 0
    for event, elem in context:
        # print elem.tag, event
        xmlns, localtag = split_namespace(elem.tag)
        if event == 'start':
            if localtag == 'phylogeny' and xmlns == NAMESPACES['phy']:
                phylogeny = Parser._parse_phylogeny(elem, context)
                phyloxml.phylogenies.append(phylogeny)
                continue
            elif xmlns != NAMESPACES['phy']:
                other_depth += 1
        if event == 'end' and xmlns != NAMESPACES['phy']:
            # Deal with items not specified by phyloXML
            other_depth -= 1
            if other_depth == 0:
                # We're directly under the root node -- evaluate
                otr = Parser.to_other(elem)
                phyloxml.other.append(otr)
                root.clear()
    return phyloxml


def parse(handle):
    """Iterate over the phylogenetic trees in a phyloXML file.

    This ignores any additional data stored at the top level, but may be more
    memory-efficient than the read() function.
    """
    context = iter(ElementTree.iterparse(handle, events=('start', 'end')))
    event, root = context.next()
    for event, elem in context:
        if event == 'start' and local(elem.tag) == 'phylogeny':
            yield Parser._parse_phylogeny(elem, context)


class Parser(object):

    ## Index of phyloxml tags and corresponding classes

    ## special cases for parsing
    # tags_to_special = {
    #         'phylogeny':    Phylogeny,
    #         'phyloxml':     Phyloxml,
    #         'clade':        Clade,
    #         }

    ## no special handling
    # tags_to_classes = {
    #         'absent':       BinaryCharacterList,
    #         'accession':    Accession,
    #         'annotation':   Annotation,
    #         'binary_characters': BinaryCharacters,
    #         'clade_relation': CladeRelation,
    #         'code':         TaxonomyCode,
    #         'color':        BranchColor,
    #         'confidence':   Confidence,
    #         'date':         Date,
    #         'distribution': Distribution,
    #         'domain':       ProteinDomain,
    #         'domain_architecture': DomainArchitecture,
    #         'events':       Events,
    #         'gained':       BinaryCharacterList,
    #         'id':           Id,
    #         'lost':         BinaryCharacterList,
    #         'mol_seq':      MolSeq,
    #         'node_id':      Id,
    #         'point':        Point,
    #         'polygon':      Polygon,
    #         'present':      BinaryCharacterList,
    #         'property':     Property,
    #         'rank':         Rank,
    #         'reference':    Reference,
    #         'sequence':     Sequence,
    #         'sequence_relation': SequenceRelation,
    #         'symbol':       SequenceSymbol,
    #         'taxonomy':     Taxonomy,
    #         'type':         EventType,
    #         'uri':          Uri,
    #         }

    ## primitive types
    # tags_to_float = {
    #         'alt':          float,  # decimal
    #         'branch_length': float, # double
    #         'lat':          float,  # decimal
    #         'long':         float,  # decimal
    #         'value':        float,  # decimal
    #         'width':        float,  # double
    #         }

    # tags_to_int = {
    #         'blue':         int, # unsignedByte
    #         'duplications': int, # nonNegativeInteger
    #         'green':        int, # unsignedByte
    #         'losses':       int, # nonNegativeInteger
    #         'red':          int, # unsignedByte
    #         'speciations':  int, # nonNegativeInteger
    #         }

    # tags_to_str = {
    #         'bc':           str,
    #         'common_name':  str,
    #         'desc':         str,
    #         'description':  str,
    #         'location':     str,
    #         'name':         str,
    #         'scientific_name': str,
    #         }

    @classmethod
    def _parse_phylogeny(cls, parent, context):
        """Parse a single phylogeny within the phyloXML tree.

        Recursively builds a phylogenetic tree with help from parse_clade, then
        clears the XML event history for the phylogeny element and returns
        control to the top-level parsing function.
        """
        phylogeny = Tree.Phylogeny(**dict_str2bool(parent.attrib,
                                                   ['rooted', 'rerootable']))
        complex_types = ['date', 'id']
        list_types = {
                # XML tag, plural attribute
                'confidence':   'confidences',
                'property':     'properties',
                'clade_relation': 'clade_relations',
                'sequence_relation': 'sequence_relations',
                }
        for event, elem in context:
            namespace, tag = split_namespace(elem.tag)
            if event == 'start' and tag == 'clade':
                assert phylogeny.clade is None, \
                        "Phylogeny object should only have 1 clade"
                phylogeny.clade = Parser._parse_clade(elem, context)
                continue
            if event == 'end':
                if tag == 'phylogeny':
                    parent.clear()
                    break
                # Handle the other non-recursive children
                if tag in complex_types:
                    setattr(phylogeny, tag, getattr(cls, 'to_'+tag)(elem))
                elif tag in list_types:
                    getattr(phylogeny, list_types[tag]).append(
                            getattr(cls, 'to_'+tag)(elem))
                # Simple types
                elif tag == 'name': 
                    phylogeny.name = elem.text and elem.text.strip()
                elif tag == 'description':
                    phylogeny.description = elem.text and elem.text.strip()
                # Unknown tags
                elif namespace != NAMESPACES['phy']:
                    phylogeny.other.append(cls.to_other(elem))
                    parent.clear()
                else:
                    # NB: This shouldn't happen in valid files
                    raise PhyloXMLError('Misidentified tag: ' + tag)
        return phylogeny

    @classmethod
    def _parse_clade(cls, parent, context):
        if 'branch_length' in parent.keys():
            parent.set('branch_length', float(parent.get('branch_length')))
        clade = Tree.Clade(**parent.attrib)
        simple_types = ['branch_length', 'name', 'node_id', 'width']
        complex_types = ['color', 'events', 'binary_characters', 'date']
        list_types = {
                # XML tag, plural attribute
                'confidence':   'confidences',
                'taxonomy':     'taxonomies',
                'sequence':     'sequences',
                'distribution': 'distributions',
                'reference':    'references',
                'property':     'properties',
                }
        # NB: Only evaluate nodes at the current level
        tracked_tags = set(simple_types + complex_types + list_types.keys())
        last_tag = None
        for event, elem in context:
            namespace, tag = split_namespace(elem.tag)
            if event == 'start':
                if tag == 'clade':
                    subclade = Parser._parse_clade(elem, context)
                    clade.clades.append(subclade)
                    continue
                if tag in tracked_tags:
                    last_tag = tag
            if event == 'end':
                if tag == 'clade':
                    elem.clear()
                    break
                if tag != last_tag or tag not in tracked_tags:
                    continue
                # Handle the other non-recursive children
                if tag in complex_types:
                    setattr(clade, tag, getattr(cls, 'to_'+tag)(elem))
                elif tag in list_types:
                    getattr(clade, list_types[tag]).append(
                            getattr(cls, 'to_'+tag)(elem))
                # Simple types
                elif tag in simple_types:
                    if tag == 'branch_length':
                        # NB: possible collision with the attribute
                        if hasattr(clade, 'branch_length') \
                                and clade.branch_length is not None:
                            warnings.warn(
                                    'Attribute branch_length was already set for '
                                    'this Clade; overwriting the previous value.',
                                    PhyloXMLWarning)
                        clade.branch_length = float(elem.text.strip())
                    elif tag == 'width':
                        clade.width = float(elem.text)
                    else: # name or node_id
                        setattr(clade, tag, elem.text and elem.text.strip())
                # Unknown tags
                elif namespace != NAMESPACES['phy']:
                    clade.other.append(cls.to_other(elem))
                    elem.clear()
                else:
                    # NB: This shouldn't happen in valid files
                    raise PhyloXMLError('Misidentified tag: ' + tag)
        return clade

    @classmethod
    def to_other(cls, elem):
        namespace, localtag = split_namespace(elem.tag)
        return Tree.Other(localtag, namespace, elem.attrib,
                  value=(elem.text and elem.text.strip() or None),
                  children=[cls.to_other(child) for child in elem])

    # Complex types

    @classmethod
    def to_accession(cls, elem):
        return Tree.Accession(elem.text.strip(), elem.get('source'))

    @classmethod
    def to_annotation(cls, elem):
        return Tree.Annotation(
                desc=get_child_text(elem, 'desc'),
                confidence=get_child_as(elem, 'confidence', cls.to_confidence),
                properties=get_children_as(elem, 'property', cls.to_property),
                uri=get_child_as(elem, 'uri', cls.to_uri),
                **elem.attrib)

    @classmethod
    def to_binary_characters(cls, elem):
        raise NotImplementedError
        return Tree.BinaryCharacters(
                # TODO
                )

    @classmethod
    def to_clade_relation(cls, elem):
        return Tree.CladeRelation(
                elem.get('type'), elem.get('id_ref_0'), elem.get('id_ref_1'),
                distance=elem.get('distance'),
                confidence=get_child_as(elem, 'confidence', cls.to_confidence))

    @classmethod
    def to_color(cls, elem):
        red, green, blue = (get_child_text(elem, color, int) for color in
                            ('red', 'green', 'blue'))
        return Tree.BranchColor(red, green, blue)

    @classmethod
    def to_confidence(cls, elem):
        return Tree.Confidence(float(elem.text), elem.get('type'))

    @classmethod
    def to_date(cls, elem):
        return Tree.Date(
                value=get_child_text(elem, 'value', float),
                desc=get_child_text(elem, 'desc'),
                unit=elem.get('unit'),
                range=(('range' in elem.keys())
                       and float(elem.get('range')) or None))

    @classmethod
    def to_distribution(cls, elem):
        return Tree.Distribution(
                desc=get_child_text(elem, 'desc'),
                points=get_children_as(elem, 'point', cls.to_point),
                polygons=get_children_as(elem, 'polygon', cls.to_polygon))

    @classmethod
    def to_domain(cls, elem):
        return Tree.ProteinDomain(elem.text.strip(),
                int(elem.get('from')), int(elem.get('to')),
                confidence=(('confidence' in elem.keys())
                            and float(elem.get('confidence')) or None),
                id=elem.get('id'))

    @classmethod
    def to_domain_architecture(cls, elem):
        return Tree.DomainArchitecture(
                length=int(elem.get('length')),
                domains=get_children_as(elem, 'domain', cls.to_domain))

    @classmethod
    def to_events(cls, elem):
        return Tree.Events(
                type=check_str(get_child_text(elem, 'type'),
                               r'(%s)' % '|'.join((
                                   'transfer', 'fusion',
                                   'speciation_or_duplication', 'other',
                                   'mixed', 'unassigned',
                                   ))),
                duplications=get_child_text(elem, 'duplications', int),
                speciations=get_child_text(elem, 'speciations', int),
                losses=get_child_text(elem, 'losses', int),
                confidence=get_child_as(elem, 'confidence', cls.to_confidence))

    @classmethod
    def to_id(cls, elem):
        return Tree.Id(elem.text.strip(), elem.get('type'))

    @classmethod
    def to_point(cls, elem):
        return Tree.Point(
                elem.get('geodetic_datum'),
                get_child_text(elem, 'lat', float),
                get_child_text(elem, 'long', float),
                alt=get_child_text(elem, 'alt', float),
                )

    @classmethod
    def to_polygon(cls, elem):
        raise NotImplementedError
        return Tree.Polygon(
                # TODO
                )

    @classmethod
    def to_property(cls, elem):
        return Tree.Property(elem.text.strip(),
                elem.get('ref'), elem.get('applies_to'), elem.get('datatype'),
                unit=elem.get('unit'),
                id_ref=elem.get('id_ref'))

    @classmethod
    def to_reference(cls, elem):
        raise NotImplementedError
        return Tree.Reference(
                # TODO
                )

    @classmethod
    def to_sequence(cls, elem):
        return Tree.Sequence(
                symbol=check_str(get_child_text(elem, 'symbol'), r'\S{1,10}'),
                accession=get_child_as(elem, 'accession', cls.to_accession),
                name=get_child_text(elem, 'name'),
                location=get_child_text(elem, 'location'),
                mol_seq=check_str(get_child_text(elem, 'mol_seq'),
                                  r'[a-zA-Z\.\-\?\*_]+'),
                uri=get_child_as(elem, 'uri', cls.to_uri),
                domain_architecture=get_child_as(elem, 'domain_architecture',
                                                 cls.to_domain_architecture),
                annotations=get_children_as(elem, 'annotation',
                    cls.to_annotation),
                # TODO: handle "other"
                other=[],
                **elem.attrib)

    @classmethod
    def to_sequence_relation(cls, elem):
        return Tree.SequenceRelation(
                elem.get('type'), elem.get('id_ref_0'), elem.get('id_ref_1'),
                distance=elem.get('distance'),
                confidence=get_child_as(elem, 'confidence', cls.to_confidence))

    @classmethod
    def to_taxonomy(cls, elem):
        return Tree.Taxonomy(
                id=get_child_as(elem, 'id', cls.to_id),
                code=check_str(get_child_text(elem, 'code'),
                               r'[a-zA-Z0-9_]{2,10}'),
                scientific_name=get_child_text(elem, 'scientific_name'),
                common_names=get_children_text(elem, 'common_name'),
                rank=check_str(get_child_text(elem, 'rank'),
                               r'(%s)' % '|'.join((
                                   'domain', 'kingdom', 'subkingdom', 'branch',
                                   'infrakingdom', 'superphylum', 'phylum',
                                   'subphylum', 'infraphylum', 'microphylum',
                                   'superdivision', 'division', 'subdivision',
                                   'infradivision', 'superclass', 'class',
                                   'subclass', 'infraclass', 'superlegion',
                                   'legion', 'sublegion', 'infralegion',
                                   'supercohort', 'cohort', 'subcohort',
                                   'infracohort', 'superorder', 'order',
                                   'suborder', 'superfamily', 'family',
                                   'subfamily', 'supertribe', 'tribe',
                                   'subtribe', 'infratribe', 'genus',
                                   'subgenus', 'superspecies', 'species',
                                   'subspecies', 'variety', 'subvariety',
                                   'form', 'subform', 'cultivar', 'unknown',
                                   'other'))),
                uri=get_child_as(elem, 'uri', cls.to_uri),
                **elem.attrib)

    @classmethod
    def to_uri(cls, elem):
        return Tree.Uri(elem.text.strip(),
                desc=elem.get('desc'), type=elem.get('type'))

