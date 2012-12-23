# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# Based on Bio.Nexus, copyright 2005-2008 by Frank Kauff & Cymon J. Cox.
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
            node_dict = {}
            children = {}
            
            nodes = tree.get_node()
            root = None
            for node in nodes:
                this_node = node_dict[node.id] = {}
                if hasattr(node, 'label') and node.label: this_node['name'] = node.label
                if node.root: root = node.id

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
                rooted = False
                possible_roots = (node.id for node in nodes if node.id in srcs and not node.id in tars)
                root = possible_roots.next()
            else:
                rooted = True
                
            yield Newick.Tree(root=self._make_tree(root, node_dict, children), rooted=rooted)
            
    @classmethod
    def _make_tree(cls, node, node_dict, children):
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

    def write(self, handle, **kwargs):
        """Write this instance's trees to a file handle."""

        # TODO: write trees to handle

        return count
