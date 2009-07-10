# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""PhyloXML writer and associated functions.

Constructs an XML file from a Tree.PhyloXML object.
"""

import Tree
from Parser import ElementTree, NAMESPACES


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


def _handle_complex(tag, attribs, complex_types, list_types, has_text=False):
    def wrapped(self, obj):
        elem = ElementTree.Element(tag, _clean_attrib(obj, attribs))
        for ct in complex_types:
            if getattr(obj, ct) is not None:
                elem.append(getattr(self, ct)(getattr(obj, ct)))
        for lt, plural in list_types:
            for item in getattr(obj, plural):
                elem.append(getattr(self, lt)(item))
        if has_text:
            elem.text = serialize(obj.value)
        return elem
    return wrapped


def _handle_simple(tag):
    def wrapped(self, obj):
        elem = ElementTree.Element(tag)
        elem.text = serialize(obj)
        return elem
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
                {_ns('schemaLocation', NAMESPACES['xs']):
                    obj.attributes['schemaLocation'],
                    })
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
            ('name', 'id', 'description', 'date', 'clade'),
            ( ('confidence',        'confidences'),
              ('clade_relation',    'clade_relations'),
              ('sequence_relation', 'sequence_relations'),
              ('property',          'properties'),
              ('other',             'other'),
              ))

    clade = _handle_complex(_ns('clade'),
            ('id_source',),
            ('name', 'branch_length', 'width', 'color', 'node_id', 'events',
                'binary_characters', 'date'),
            ( ('confidence',    'confidences'),
              ('taxonomy',      'taxonomies'),
              ('sequence',      'sequences'),
              ('distribution',  'distributions'),
              ('reference',     'references'),
              ('property',      'properties'),
              ('clade',         'clades'),
              ('other',         'other'),
              ))

    # absent = _handle_complex(_ns('absent')) # BinaryCharacterList

    accession = _handle_complex(_ns('accession'),
            ('source',), (), (), has_text=True)

    annotation = _handle_complex(_ns('annotation'),
            ('ref', 'source', 'evidence', 'type'),
            ('desc', 'confidence', 'uri'),
            ( ('property',   'properties'),
              ))

    # binary_characters = _handle_complex(_ns('binary_characters'))

    clade_relation = _handle_complex(_ns('clade_relation'),
            ('id_ref_0', 'id_ref_1', 'distance', 'type'),
            ('confidence',), ())

    # color = _handle_complex(_ns('color'))   # BranchColor,

    confidence = _handle_complex(_ns('confidence'),
            ('type',), (), (), has_text=True)

    date = _handle_complex(_ns('date'),
            ('unit', 'range'), ('desc', 'value'), (), has_text=True)

    distribution = _handle_complex(_ns('distribution'), (), ('desc',),
            ( ('point',     'points'),
              ('polygon',   'polygons'),
              ))

    def domain(self, obj):
        elem = ElementTree.Element(_ns('domain'),
                {'from': str(obj.start + 1), 'to': str(obj.end)})
        if obj.confidence is not None:
            elem.set('confidence', serialize(obj.confidence))
        if obj.id is not None:
            elem.set('id', obj.id)
        elem.text = serialize(obj.value)
        return elem

    domain_architecture = _handle_complex(_ns('domain_architecture'),
            ('length',), (),
            ( ('domain', 'domains'),
              ))

    events = _handle_complex(_ns('events'), (),
            ('type', 'duplications', 'speciations', 'losses', 'confidence'),
            ())

    # gained = _handle_complex(_ns('gained')) # BinaryCharacterList,

    id = _handle_complex(_ns('id'), ('type',), (), (), has_text=True)

    node_id = _handle_complex(_ns('node_id'), ('type',), (), (), has_text=True)

    # lost = _handle_complex(_ns('lost')) # BinaryCharacterList,

    point = _handle_complex(_ns('point'), ('geodetic_datum',),
            ('lat', 'long', 'alt'), ())

    # polygon = _handle_complex(_ns('polygon'))

    # present = _handle_complex(_ns('present')) # BinaryCharacterList,

    property = _handle_complex(_ns('property'),
            ('ref', 'unit', 'datatype', 'applies_to', 'id_ref'),
            (), (), has_text=True)

    # reference = _handle_complex(_ns('reference'))

    sequence = _handle_complex(_ns('sequence'),
            ('type', 'id_ref', 'id_source'),
            ('symbol', 'accession', 'name', 'location', 'mol_seq', 'uri',
                'domain_architecture'),
            ( ('annotation', 'annotations'),
              ('other', 'other'),
              ))

    sequence_relation = _handle_complex(_ns('sequence_relation'),
            ('id_ref_0', 'id_ref_1', 'distance', 'type'),
            ('confidence',), ())

    taxonomy = _handle_complex(_ns('taxonomy'),
            ('id_source',),
            ('id', 'code', 'scientific_name', 'rank', 'uri'),
            ( ('common_name',   'common_names'),
              ('other',         'other'),
              ))

    uri = _handle_complex(_ns('uri'), ('desc', 'type'), (), (), has_text=True)

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

