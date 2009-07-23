# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""PhyloXML writer and associated functions.

Constructs an XML file from a Tree.PhyloXML object.
"""

import warnings

from Bio.Tree import PhyloXMLTree as Tree

from Parser import ElementTree, NAMESPACES

# Keep the standard namespace prefixes when writing
# See http://effbot.org/zone/element-namespaces.htm
try:
    register_namespace = ElementTree.register_namespace
except AttributeError:
    if not hasattr(ElementTree, '_namespace_map'):
        # cElementTree needs the pure-Python xml.etree.ElementTree
        # Py2.4 support: the exception handler can go away when Py2.4 does
        try:
            from xml.etree import ElementTree as ET_py
            ElementTree._namespace_map = ET_py._namespace_map
        except ImportError:
            warnings.warn("Couldn't import xml.etree.ElementTree; "
                    "phyloXML namespaces may have unexpected abbreviations "
                    "in the output.", RuntimeWarning)
            ElementTree._namespace_map = {}

    def register_namespace(prefix, uri):
        ElementTree._namespace_map[uri] = prefix

for prefix, uri in NAMESPACES.iteritems():
    register_namespace(prefix, uri)


def write(phyloxml, file, encoding=None):
    """Write a phyloXML file.

    The file argument can be either an open handle or a file name.
    """
    etree = Writer(phyloxml).get_etree()
    if encoding is not None:
        etree.write(file, encoding)
    else:
        etree.write(file)


# Helpers

def _ns(tag, namespace=NAMESPACES['phy']):
    """Format an XML tag with the given namespace."""
    return '{%s}%s' % (namespace, tag)


def serialize(value):
    if isinstance(value, float):
        return unicode(value).upper()
    elif isinstance(value, bool):
        return unicode(value).lower()
    return unicode(value)


def _clean_attrib(obj, attrs):
    """Create a dictionary from an object's specified, non-None attributes."""
    out = {}
    for key in attrs:
        val = getattr(obj, key)
        if val is not None:
            out[key] = serialize(val)
    return out


def _handle_complex(tag, attribs, subnodes, has_text=False):
    def wrapped(self, obj):
        elem = ElementTree.Element(tag, _clean_attrib(obj, attribs))
        for subn in subnodes:
            if isinstance(subn, basestring):
                # singular object: method and attribute names are the same
                if getattr(obj, subn) is not None:
                    elem.append(getattr(self, subn)(getattr(obj, subn)))
            else:
                # list: singular method, pluralized attribute name
                method, plural = subn
                for item in getattr(obj, plural):
                    elem.append(getattr(self, method)(item))
        if has_text:
            elem.text = serialize(obj.value)
        return elem
    wrapped.__doc__ = "Serialize a %s and its subnodes, in order." % tag
    return wrapped


def _handle_simple(tag):
    def wrapped(self, obj):
        elem = ElementTree.Element(tag)
        elem.text = serialize(obj)
        return elem
    wrapped.__doc__ = "Serialize a simple %s node." % tag
    return wrapped


class Writer(object):
    """
    """
    def __init__(self, phyloxml):
        """Build an ElementTree from a phyloXML object."""
        assert isinstance(phyloxml, Tree.Phyloxml)
        self._tree = ElementTree.ElementTree(self.phyloxml(phyloxml))

    def get_etree(self):
        return self._tree

    # Convert classes to ETree elements

    def phyloxml(self, obj):
        elem = ElementTree.Element(_ns('phyloxml'),
                # XXX not sure about this
                # {_ns('schemaLocation', NAMESPACES['xs']):
                #     obj.attributes['schemaLocation'],
                #     }
                )
        for tree in obj.phylogenies:
            elem.append(self.phylogeny(tree))
        for otr in obj.other:
            elem.append(self.other(otr))
        return elem

    def other(self, obj):
        elem = ElementTree.Element(_ns(obj.tag, obj.namespace), obj.attributes)
        elem.text = obj.value
        for child in obj.children:
            elem.append(self.other(child))
        return elem

    phylogeny = _handle_complex(_ns('phylogeny'),
            ('rooted', 'rerootable', 'branch_length_unit', 'type'),
            ( 'name',
              'id',
              'description',
              'date',
              ('confidence',        'confidences'),
              'clade',
              ('clade_relation',    'clade_relations'),
              ('sequence_relation', 'sequence_relations'),
              ('property',          'properties'),
              ('other',             'other'),
              ))

    clade = _handle_complex(_ns('clade'), ('id_source',),
            ( 'name',
              'branch_length',
              ('confidence',    'confidences'),
              'width',
              'color',
              'node_id',
              ('taxonomy',      'taxonomies'),
              ('sequence',      'sequences'),
              'events',
              'binary_characters',
              ('distribution',  'distributions'),
              'date',
              ('reference',     'references'),
              ('property',      'properties'),
              ('clade',         'clades'),
              ('other',         'other'),
              ))

    accession = _handle_complex(_ns('accession'), ('source',),
            (), has_text=True)

    annotation = _handle_complex(_ns('annotation'),
            ('ref', 'source', 'evidence', 'type'),
            ( 'desc',
              'confidence',
              ('property',   'properties'),
              'uri',
              ))

    def binary_characters(self, obj):
        """Serialize a binary_characters node and its subnodes."""
        elem = ElementTree.Element(_ns('binary_characters'),
                _clean_attrib(obj,
                    ('type', 'gained_count', 'lost_count',
                        'present_count', 'absent_count')))
        for subn in ('gained', 'lost', 'present', 'absent'):
            subelem = ElementTree.Element(_ns(subn))
            for token in getattr(obj, subn):
                subelem.append(self.bc(token))
            elem.append(subelem)
        return elem

    clade_relation = _handle_complex(_ns('clade_relation'),
            ('id_ref_0', 'id_ref_1', 'distance', 'type'),
            ('confidence',))

    color = _handle_complex(_ns('color'), (), ('red', 'green', 'blue'))

    confidence = _handle_complex(_ns('confidence'), ('type',),
            (), has_text=True)

    date = _handle_complex(_ns('date'), ('unit', 'range'),
            ('desc', 'value'))

    distribution = _handle_complex(_ns('distribution'), (),
            ( 'desc',
              ('point',     'points'),
              ('polygon',   'polygons'),
              ))

    def domain(self, obj):
        """Serialize a domain node."""
        elem = ElementTree.Element(_ns('domain'),
                {'from': str(obj.start + 1), 'to': str(obj.end)})
        if obj.confidence is not None:
            elem.set('confidence', serialize(obj.confidence))
        if obj.id is not None:
            elem.set('id', obj.id)
        elem.text = serialize(obj.value)
        return elem

    domain_architecture = _handle_complex(_ns('domain_architecture'),
            ('length',),
            (('domain', 'domains'),))

    events = _handle_complex(_ns('events'), (),
            ( 'type',
              'duplications',
              'speciations',
              'losses',
              'confidence',
              ))

    id = _handle_complex(_ns('id'), ('type',), (), has_text=True)

    node_id = _handle_complex(_ns('node_id'), ('type',), (), has_text=True)

    point = _handle_complex(_ns('point'), ('geodetic_datum',),
            ('lat', 'long', 'alt'))

    polygon = _handle_complex(_ns('polygon'), (), (('point', 'points'),))

    property = _handle_complex(_ns('property'),
            ('ref', 'unit', 'datatype', 'applies_to', 'id_ref'),
            (), has_text=True)

    reference = _handle_complex(_ns('reference'), ('doi',), ('desc',))

    sequence = _handle_complex(_ns('sequence'),
            ('type', 'id_ref', 'id_source'),
            ( 'symbol',
              'accession',
              'name',
              'location',
              'mol_seq',
              'uri',
              ('annotation', 'annotations'),
              'domain_architecture',
              ('other', 'other'),
              ))

    sequence_relation = _handle_complex(_ns('sequence_relation'),
            ('id_ref_0', 'id_ref_1', 'distance', 'type'),
            ('confidence',))

    taxonomy = _handle_complex(_ns('taxonomy'),
            ('id_source',),
            ( 'id',
              'code',
              'scientific_name',
              ('common_name',   'common_names'),
              'rank',
              'uri',
              ('other',         'other'),
              ))

    uri = _handle_complex(_ns('uri'), ('desc', 'type'), (), has_text=True)

    # Primitive types

    # Floating point
    alt = _handle_simple(_ns('alt'))
    branch_length = _handle_simple(_ns('branch_length'))
    lat = _handle_simple(_ns('lat'))
    long = _handle_simple(_ns('long'))
    value = _handle_simple(_ns('value'))
    width = _handle_simple(_ns('width'))

    # Integers
    blue = _handle_simple(_ns('blue'))
    duplications = _handle_simple(_ns('duplications'))
    green = _handle_simple(_ns('green'))
    losses = _handle_simple(_ns('losses'))
    red = _handle_simple(_ns('red'))
    speciations = _handle_simple(_ns('speciations'))

    # Strings
    bc = _handle_simple(_ns('bc'))
    code = _handle_simple(_ns('code'))      # TaxonomyCode
    common_name = _handle_simple(_ns('common_name'))
    desc = _handle_simple(_ns('desc'))
    description = _handle_simple(_ns('description'))
    location = _handle_simple(_ns('location'))
    mol_seq = _handle_simple(_ns('mol_seq'))
    name = _handle_simple(_ns('name'))
    rank = _handle_simple(_ns('rank')) # Rank
    scientific_name = _handle_simple(_ns('scientific_name'))
    symbol = _handle_simple(_ns('symbol'))
    type = _handle_simple(_ns('type')) # EventType

