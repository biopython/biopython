# Copyright (C) 2012 by Ben Morris (ben@bendmorris.com)
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

from cStringIO import StringIO

from Bio.Phylo import NeXML
import xml.etree.ElementTree as ET
from xml.dom import minidom

NAMESPACES = {
              'xsi': 'http://www.w3.org/2001/XMLSchema-instance',
              #'xmi': 'http://www.w3.org/XML/1998/namespace',
              'nex': 'http://www.nexml.org/2009',
              'cdao': 'http://purl.obolibrary.org/obo/cdao.owl#',
              'xsd': 'http://www.w3.org/2001/XMLSchema#',
              }
DEFAULT_NAMESPACE = NAMESPACES['nex']

for prefix, uri in NAMESPACES.items():
    ET.register_namespace(prefix, uri)
    
def qUri(s):
    for prefix, uri in NAMESPACES.items():
        s = s.replace('%s:' % prefix, '{%s}' % uri)
    return s


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

    def parse(self, values_are_confidence=False, rooted=False):
        """Parse the text stream this object was initialized with."""

        nexml_doc = ET.iterparse(self.handle, events=('end',))
        
        for event, node in nexml_doc:
            if node.tag == qUri('nex:tree'):
                node_dict = {}
                node_children = {}
                root = None
                
                child_tags = node._children
                nodes = []                
                edges = []
                for child in child_tags:
                    if child.tag == qUri('nex:node'): nodes.append(child)
                    if child.tag == qUri('nex:edge'): edges.append(child)
                    
                for node in nodes:
                    node.id = node.attrib['id']
                    this_node = node_dict[node.id] = {}
                    if 'otu' in node.attrib and node.attrib['otu']: this_node['name'] = node.attrib['otu']
                    if 'root' in node.attrib and node.attrib['root'] == 'true': root = node.id
                    
                srcs = set()
                tars = set()
                for edge in edges:
                    src, tar = edge.attrib['source'], edge.attrib['target']
                    srcs.add(src)
                    tars.add(tar)
                    if not src in node_children: node_children[src] = set()
                    
                    node_children[src].add(tar)
                    if 'length' in edge.attrib: node_dict[tar]['branch_length'] = edge.attrib['length']
                    
                if root is None:
                    # if no root specified, start the recursive tree creation function
                    # with the first node that's not a child of any other nodes
                    rooted = False
                    possible_roots = (node.id for node in nodes if node.id in srcs and not node.id in tars)
                    root = possible_roots.next()
                else:
                    rooted = True
                    
                yield NeXML.Tree(root=self._make_tree(root, node_dict, node_children), rooted=rooted)
                
            
    @classmethod
    def _make_tree(cls, node, node_dict, children):
        '''Return a NeXML.Clade, and calls itself recursively for each child, 
        traversing the  entire tree and creating a nested structure of NeXML.Clade 
        objects.'''
        
        this_node = node_dict[node]
        clade = NeXML.Clade(**this_node)
        
        if node in children:
            clade.clades = [cls._make_tree(child, node_dict, children) for child in children[node]]
        
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

    def write(self, handle, **kwargs):
        """Write this instance's trees to a file handle."""
        
        # TODO: this is not handling XML namespaces and the root nex:nexml node correctly
        
        # set XML namespaces
        root_node = ET.Element('nex:nexml')
        root_node.set('version', '0.9')
        root_node.set('xmlns', DEFAULT_NAMESPACE)
        root_node.set('xsi:schemaLocation', 'http://www.nexml.org/2009/nexml/xsd/nexml.xsd')

        for prefix, uri in NAMESPACES.items():
            root_node.set('xmlns:%s' % prefix, uri)

        otus = ET.SubElement(root_node, 'otus', attrib={'id': 'tax', 'label': 'RootTaxaBlock'})
        
        # create trees
        trees = ET.SubElement(root_node, 'trees', attrib={'id':'Trees', 'label':'TreesBlockFromXML', 'otus': 'tax'})
        count = 0
        tus = set()
        for tree in self.trees:
            this_tree = ET.SubElement(trees, 'tree', attrib={'id':self.new_label('tree')})
            
            first_clade = tree.clade
            tus.update(self._write_tree(first_clade, this_tree, rooted=tree.rooted))

            count += 1
        
        # create OTUs
        for tu in tus:
            otu = ET.SubElement(otus, 'otu', attrib={'id':tu})
        
        # write XML document to file handle
        xml_doc = ET.ElementTree(root_node)
        #xml_doc.write(handle,
        #              xml_declaration=True, encoding='utf-8',
        #              method='xml')

        # use xml.dom.minodom for pretty printing
        rough_string = ET.tostring(root_node, 'utf-8')
        reparsed = minidom.parseString(rough_string)
        handle.write(reparsed.toprettyxml(indent="  "))
        
        return count
    
    def _write_tree(self, clade, tree, parent=None, rooted=False):
        '''Recursively process tree, adding nodes and edges to Tree object. 
        Returns a set of all OTUs encountered.'''
        tus = set()
        
        node_id = self.new_label('node')
        clade.node_id = node_id
        attrib={'id':node_id, 'label':node_id}
        root = rooted and parent is None
        if root: attrib['root'] = 'true'
        if clade.name:
            tus.add(clade.name)
            attrib['otu'] = clade.name
        node = ET.SubElement(tree, 'node', attrib=attrib)
        
        if not parent is None:
            edge_id = self.new_label('edge')
            node = ET.SubElement(tree, 'edge', attrib={'id':edge_id, 'source':parent.node_id, 'target':node_id,
                                                       'length':str(clade.branch_length)})
    
        if not clade.is_terminal():
            for new_clade in clade.clades:
                tus.update(self._write_tree(new_clade, tree, parent=clade))
                
        return tus
