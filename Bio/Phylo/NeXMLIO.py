# Copyright (C) 2012 by Ben Morris (ben@bendmorris.com)
# Based on Bio.Nexus, copyright 2005-2008 by Frank Kauff & Cymon J. Cox
# and Bio.Newick, copyright 2009 by Eric Talevich.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""I/O function wrappers for the NeXML file format.

See: http://www.nexml.org
"""
__docformat__ = "restructuredtext en"

from cStringIO import StringIO

from Bio.Phylo import Newick, _nexml_gds as gds


class NeXMLError(Exception):
    """Exception raised when NeXML object construction cannot continue."""
    pass


# ---------------------------------------------------------
# Public API

def parse(handle, **kwargs):
    """Iterate over the trees in a Newick file handle.

    :returns: generator of Bio.Phylo.Newick.Tree objects.
    """
    return Parser(handle).parse(**kwargs)


def write(trees, handle, plain=False, **kwargs):
    """Write a trees in Newick format to the given file handle.

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

        nexml_doc = gds.parseString(self.handle.read())
        trees = nexml_doc.get_trees()[0].get_tree()
        for tree in trees:
            print tree
            node_dict = {}
            children = {}
            
            # create dictionary of all nodes in this tree
            nodes = tree.get_node()
            root = None
            for node in nodes:
                this_node = node_dict[node.id] = {}
                if hasattr(node, 'otu') and node.otu: this_node['name'] = node.otu
                if node.root: root = node.id
            
            # create dictionary linking each node to all of its children
            edges = tree.get_edge()
            srcs = set()
            tars = set()
            for edge in edges:
                src, tar = edge.source, edge.target
                srcs.add(src)
                tars.add(tar)
                if not src in children: children[src] = set()
                
                children[src].add(tar)
                node_dict[tar]['branch_length'] = edge.length
                
            if root is None:
                # if no root specified, start the recursive tree creation function
                # with the first node that's not a child of any other nodes
                rooted = False
                possible_roots = (node.id for node in nodes if node.id in srcs and not node.id in tars)
                root = possible_roots.next()
            else:
                rooted = True
                
            yield Newick.Tree(root=self._make_tree(root, node_dict, children), rooted=rooted)
            
    @classmethod
    def _make_tree(cls, node, node_dict, children):
        '''Return a Newick.Clade, and calls itself recursively for each child, 
        traversing the  entire tree and creating a nested structure of Newick.Clade 
        objects.'''
        
        this_node = node_dict[node]
        clade = Newick.Clade(**this_node)
        
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
        self.tu_counter = 0
        
    def new_label(self, obj_type):
        counter = '%s_counter' % obj_type
        setattr(self, counter, getattr(self, counter) + 1)
        return '%s%s' % (obj_type, getattr(self, counter))

    def write(self, handle, **kwargs):
        """Write this instance's trees to a file handle."""

        xml_doc = gds.Nexml()
        trees = gds.Trees()
        count = 0
        tus = set()
        for tree in self.trees:
            this_tree = gds.FloatTree(id=self.new_label('tree'))
            
            first_clade = tree.clade
            tus.update(self._write_tree(first_clade, this_tree))
            
            trees.add_tree(this_tree)
            count += 1
            
        for tu in tus:
            # TODO: add OTUs to trees object
            pass
            
        xml_doc.set_trees([trees])
        xml_doc.export(outfile=handle, level=0)

        return count
    
    def _write_tree(self, clade, tree, parent=None):
        '''Recursively process tree, adding nodes and edges to Tree object. 
        Returns a set of all OTUs encountered.'''
        tus = set()
        if clade.name:
            tus.add(clade.name)
        
        node_id = self.new_label('node')
        # TODO: create new node, add to tree with tree.add_node
        
        if not parent is None:
            edge_id = self.new_label('edge')
            # TODO: create new edge, add to tree with tree.add_edge
    
        if not clade.is_terminal():
            for new_clade in clade.clades:
                tus.update(self._write_tree(new_clade, tree, parent=clade))
                
        return tus
