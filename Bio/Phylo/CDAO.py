# Copyright (C) 2013 by Ben Morris (ben@bendmorris.com)
# based on code by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Classes corresponding to CDAO trees.

See classes in `Bio.Nexus`: Trees.Tree, Trees.NodeData, and Nodes.Chain.
"""
__docformat__ = "restructuredtext en"

from Bio.Phylo import BaseTree


class Tree(BaseTree.Tree):
    """CDAO Tree object."""

    def __init__(self, root=None, rooted=False, id=None, name=None, weight=1.0):
        BaseTree.Tree.__init__(self, root=root or Clade(),
                               rooted=rooted, id=id, name=name)
        self.weight = weight
        # a list of (predicate, object) pairs, containing additional triples
        # using this tree as subject
        self.attributes = []


class Clade(BaseTree.Clade):
    """CDAO Clade (sub-tree) object."""

    def __init__(self, branch_length=1.0, name=None, clades=None,
                 confidence=None, comment=None):
        BaseTree.Clade.__init__(self, branch_length=branch_length,
                                name=name, clades=clades, confidence=confidence)
        self.comment = comment
        # a list of (predicate, object) pairs, containing additional triples
        # using this clade as subject
        self.attributes = []
        self.tu_attributes = []
        self.edge_attributes = []
