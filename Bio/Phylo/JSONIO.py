# Copyright (C) 2015 by Fabio Zanini (fabio.zanini __AT__ fastmail.fm)
# Based on Bio.Phylo.NewickIO, Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)

from __future__ import print_function
from Bio.Phylo import Tree


# ---------------------------------------------------------
# Public API
def parse(handle, **kwargs):
    """Parse a tree in JSON format

    :returns: a generator to a Bio.Phylo.Tree object.
    """
    return Parser(handle).parse(**kwargs)


def write(tree, handle, indent=2, **kwargs):
    """Write a tree in JSON format to the given file handle."""
    return Writer(tree).write(handle, indent, **kwargs)


# ---------------------------------------------------------
# Input
class Parser(object):
    """Parse a JSON tree given a file handle.

    Based on the Newick parser.
    """

    def __init__(self, handle):
        self.handle = handle

    @classmethod
    def from_string(cls, treetext):
        handle = StringIO(treetext)
        return cls(handle)


    def parse(self, **kwargs):
        yield self.read(**kwargs)
        

    def read(self,
             children_attrname="children",
             metadata_nodes={'name': str, 'default': str},
             metadata_tree={'name': str, 'default': str}):
        import json

        def node_from_json(node_dict, node):
            for attr in node_dict:
                if (children_attrname != 'clades') and (attr == 'clades'):
                    print('"Clades" metadata conflicts with Biopython, skipping')
                    continue

                if attr != children_attrname:
                    if attr in metadata_nodes:
                        fmt = metadata_nodes[attr]
                        setattr(node, attr, fmt(node_dict[attr]))
                        continue

                    fmt = metadata_nodes['default']
                    try:
                        setattr(node, attr, fmt(node_dict[attr]))
                    except:
                        setattr(node, attr, node_dict[attr])

                else:
                    for child_dict in node_dict[attr]:
                        child = Phylo.BaseTree.Clade()
                        node.clades.append(child)
                        node_from_json(child_dict, child)

        tree_dict = json.load(handle)
        tree = Phylo.BaseTree.Tree()

        if 'tree' in tree_dict:
            root_dict = tree_dict['tree']
            for attr in tree_dict:
                if attr == 'tree':
                    continue
                if attr in metadata_tree:
                    fmt = metadata_tree[fmt]
                    setattr(tree, attr, fmt(tree_dict[attr]))
                    continue
                fmt = metadata_tree['default']
                try:
                    setattr(tree, attr, fmt(tree_dict[attr]))
                except:
                    setattr(tree, attr, tree_dict[attr])

        else:
            root_dict = tree_dict

        node_from_json(root_dict, tree.root)
        return tree


# ---------------------------------------------------------
# Output
class Writer(object):

    def __init__(self, tree):
        self.tree = tree


    def write(self, handle, indent, **kwargs):
        """Write this instance's tree to a file handle."""
        import json
        json.dump(self.to_dict(self.tree, **kwargs), handle, indent=indent)


    def to_dict(self,
                children_attrname="children",
                metadata_nodes=[],
                metadata_tree=[],
                tree_layer=True):

        not_metadata = ['clades', children_attrname]
        metadata_nodes=[m for m in metadata_nodes if m not in not_metadata]

        def convert_to_dict(node):
            d = {}
            for field in metadata_nodes:
                if hasattr(node, field):
                    d[field] = getattr(node, field)
            d[children_attrname] = [convert_to_dict(c) for c in node.clades]
            return d

        tree_dict = convert_to_dict(tree.root)
        if tree_layer:
            not_metadata = ['tree']
            metadata_tree = [m for m in metadata_tree if m not in not_metadata]
            tree_dict = {'tree': tree_dict}
            for field in metadata_tree:
                if hasattr(tree, field):
                    tree_dict[field] = getattr(tree, field)

        return tree_dict
