# Copyright (C) 2013 by Ben Morris (ben@bendmorris.com)
# Based on Bio.Nexus, copyright 2005-2008 by Frank Kauff & Cymon J. Cox
# and Bio.Phylo.Newick, copyright 2009 by Eric Talevich.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""I/O function wrappers for the NeXML file format.

See: http://www.nexml.org
"""

__docformat__ = "restructuredtext en"

from Bio._py3k import StringIO

from Bio.Phylo import NeXML
from xml.dom import minidom
import sys
from ._cdao_owl import cdao_elements, cdao_namespaces, resolve_uri


# For speed try to use cElementTree rather than ElementTree
try:
    if (3, 0) <= sys.version_info[:2] <= (3, 1):
        # Workaround for bug in python 3.0 and 3.1,
        # see http://bugs.python.org/issue9257
        from xml.etree import ElementTree as ElementTree
    else:
        from xml.etree import cElementTree as ElementTree
except ImportError:
    from xml.etree import ElementTree as ElementTree

NAMESPACES = {
    'xsi': 'http://www.w3.org/2001/XMLSchema-instance',
    'xml': 'http://www.w3.org/XML/1998/namespace',
    'nex': 'http://www.nexml.org/2009',
    'xsd': 'http://www.w3.org/2001/XMLSchema#',
}
NAMESPACES.update(cdao_namespaces)
DEFAULT_NAMESPACE = NAMESPACES['nex']
VERSION = '0.9'
SCHEMA = 'http://www.nexml.org/2009/nexml/xsd/nexml.xsd'


try:
    register_namespace = ElementTree.register_namespace
except AttributeError:
    if not hasattr(ElementTree, '_namespace_map'):
        # cElementTree needs the pure-Python xml.etree.ElementTree
        from xml.etree import ElementTree as ET_py
        ElementTree._namespace_map = ET_py._namespace_map

    def register_namespace(prefix, uri):
        ElementTree._namespace_map[uri] = prefix

for prefix, uri in NAMESPACES.items():
    register_namespace(prefix, uri)


def qUri(s):
    """Given a prefixed URI, return the full URI."""
    return resolve_uri(s, namespaces=NAMESPACES, xml_style=True)


def cdao_to_obo(s):
    """Optionally converts a CDAO-prefixed URI into an OBO-prefixed URI."""
    return 'obo:%s' % cdao_elements[s[len('cdao:'):]]


def matches(s):
    """Check for matches in both CDAO and OBO namespaces."""
    if s.startswith('cdao:'):
        return (s, cdao_to_obo(s))
    else:
        return (s,)


class NeXMLError(Exception):
    """Exception raised when NeXML object construction cannot continue."""
    pass


# ---------------------------------------------------------
# Public API

def parse(handle, **kwargs):
    """Iterate over the trees in a NeXML file handle.

    :returns: generator of Bio.Phylo.NeXML.Tree objects.
    """
    return Parser(handle).parse(**kwargs)


def write(trees, handle, plain=False, **kwargs):
    """Write a trees in NeXML format to the given file handle.

    :returns: number of trees written.
    """
    return Writer(trees).write(handle, plain=plain, **kwargs)


# ---------------------------------------------------------
# Input

class Parser(object):
    """Parse a NeXML tree given a file handle.

    Based on the parser in `Bio.Nexus.Trees`.
    """

    def __init__(self, handle):
        self.handle = handle

    @classmethod
    def from_string(cls, treetext):
        handle = StringIO(treetext)
        return cls(handle)

    def add_annotation(self, node_dict, meta_node):
        if 'property' in meta_node.attrib:
            prop = meta_node.attrib['property']
        else:
            prop = 'meta'

        if prop in matches('cdao:has_Support_Value'):
            node_dict['confidence'] = float(meta_node.text)
        else:
            node_dict[prop] = meta_node.text

    def parse(self, values_are_confidence=False, rooted=False):
        """Parse the text stream this object was initialized with."""

        nexml_doc = ElementTree.iterparse(self.handle, events=('end',))

        for event, node in nexml_doc:
            if node.tag == qUri('nex:tree'):
                node_dict = {}
                node_children = {}
                root = None

                child_tags = node.getchildren()
                nodes = []
                edges = []
                for child in child_tags:
                    if child.tag == qUri('nex:node'):
                        nodes.append(child)
                    if child.tag == qUri('nex:edge'):
                        edges.append(child)

                for node in nodes:
                    node_id = node.attrib['id']
                    this_node = node_dict[node_id] = {}
                    if 'otu' in node.attrib and node.attrib['otu']:
                        this_node['name'] = node.attrib['otu']
                    if 'root' in node.attrib and node.attrib['root'] == 'true':
                        root = node_id

                    for child in node.getchildren():
                        if child.tag == qUri('nex:meta'):
                            self.add_annotation(node_dict[node_id], child)

                srcs = set()
                tars = set()
                for edge in edges:
                    src, tar = edge.attrib['source'], edge.attrib['target']
                    srcs.add(src)
                    tars.add(tar)
                    if src not in node_children:
                        node_children[src] = set()

                    node_children[src].add(tar)
                    if 'length' in edge.attrib:
                        node_dict[tar]['branch_length'] = float(edge.attrib['length'])
                    if 'property' in edge.attrib and edge.attrib['property'] in matches('cdao:has_Support_Value'):
                        node_dict[tar]['confidence'] = float(edge.attrib['content'])

                    for child in edge.getchildren():
                        if child.tag == qUri('nex:meta'):
                            self.add_annotation(node_dict[tar], child)

                if root is None:
                    # if no root specified, start the recursive tree creation function
                    # with the first node that's not a child of any other nodes
                    rooted = False
                    possible_roots = (node.attrib['id'] for node in nodes
                                      if node.attrib['id'] in srcs
                                      and not node.attrib['id'] in tars)
                    root = next(possible_roots)
                else:
                    rooted = True

                yield NeXML.Tree(root=self._make_tree(root, node_dict, node_children), rooted=rooted)

    @classmethod
    def _make_tree(cls, node, node_dict, children):
        """Traverse the tree creating a nested clade structure.

        Return a NeXML.Clade, and calls itself recursively for each child,
        traversing the  entire tree and creating a nested structure of NeXML.Clade
        objects.
        """

        this_node = node_dict[node]
        clade = NeXML.Clade(**this_node)

        if node in children:
            clade.clades = [cls._make_tree(child, node_dict, children)
                            for child in children[node]]

        return clade

# ---------------------------------------------------------
# Output


class Writer(object):
    """Based on the writer in Bio.Nexus.Trees (str, to_string)."""

    def __init__(self, trees):
        self.trees = trees

        self.node_counter = 0
        self.edge_counter = 0
        self.tree_counter = 0

    def new_label(self, obj_type):
        counter = '%s_counter' % obj_type
        setattr(self, counter, getattr(self, counter) + 1)
        return '%s%s' % (obj_type, getattr(self, counter))

    def write(self, handle, cdao_to_obo=True, **kwargs):
        """Write this instance's trees to a file handle."""

        self.cdao_to_obo = cdao_to_obo

        # set XML namespaces
        root_node = ElementTree.Element('nex:nexml')
        root_node.set('version', VERSION)
        root_node.set('xmlns', DEFAULT_NAMESPACE)
        root_node.set('xsi:schemaLocation', SCHEMA)

        for prefix, uri in NAMESPACES.items():
            root_node.set('xmlns:%s' % prefix, uri)

        otus = ElementTree.SubElement(root_node, 'otus',
                                      **{'id': 'tax', 'label': 'RootTaxaBlock'})

        # create trees
        trees = ElementTree.SubElement(root_node, 'trees',
                                       **{'id': 'Trees', 'label': 'TreesBlockFromXML', 'otus': 'tax'})
        count = 0
        tus = set()
        for tree in self.trees:
            this_tree = ElementTree.SubElement(trees, 'tree',
                                               **{'id': self.new_label('tree')})

            first_clade = tree.clade
            tus.update(self._write_tree(first_clade, this_tree, rooted=tree.rooted))

            count += 1

        # create OTUs
        for tu in tus:
            otu = ElementTree.SubElement(otus, 'otu', **{'id': tu})

        # write XML document to file handle
        # xml_doc = ElementTree.ElementTree(root_node)
        # xml_doc.write(handle,
        #              xml_declaration=True, encoding='utf-8',
        #              method='xml')

        # use xml.dom.minodom for pretty printing
        rough_string = ElementTree.tostring(root_node, 'utf-8')
        reparsed = minidom.parseString(rough_string)
        try:
            handle.write(reparsed.toprettyxml(indent="  "))
        except TypeError:
            # for compatibility with Python 3
            handle.write(bytes(reparsed.toprettyxml(indent="  "), 'utf8'))

        return count

    def _write_tree(self, clade, tree, parent=None, rooted=False):
        """Recursively process tree, adding nodes and edges to Tree object.

        Returns a set of all OTUs encountered.
        """
        tus = set()

        convert_uri = cdao_to_obo if self.cdao_to_obo else (lambda s: s)

        node_id = self.new_label('node')
        clade.node_id = node_id
        attrib = {'id': node_id, 'label': node_id}
        root = rooted and parent is None
        if root:
            attrib['root'] = 'true'
        if clade.name:
            tus.add(clade.name)
            attrib['otu'] = clade.name
        node = ElementTree.SubElement(tree, 'node', **attrib)

        if parent is not None:
            edge_id = self.new_label('edge')
            attrib = {
                'id': edge_id, 'source': parent.node_id, 'target': node_id,
                'length': str(clade.branch_length),
                'typeof': convert_uri('cdao:Edge'),
            }
            if hasattr(clade, 'confidence') and clade.confidence is not None:
                attrib.update({
                    'property': convert_uri('cdao:has_Support_Value'),
                    'datatype': 'xsd:float',
                    'content': '%1.2f' % clade.confidence,
                })
            node = ElementTree.SubElement(tree, 'edge', **attrib)

        if not clade.is_terminal():
            for new_clade in clade.clades:
                tus.update(self._write_tree(new_clade, tree, parent=clade))

        del clade.node_id

        return tus
