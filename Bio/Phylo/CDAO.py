# Copyright (C) 2013 by Ben Morris (ben@bendmorris.com)
# based on code by Eric Talevich (eric.talevich@gmail.com)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Classes corresponding to CDAO trees.

See classes in ``Bio.Nexus``: Trees.Tree, Trees.NodeData, and Nodes.Chain.
"""

from Bio.Phylo import BaseTree


class Tree(BaseTree.Tree):
    """CDAO Tree object."""

    def __init__(self, root=None, rooted=False, id=None, name=None, weight=1.0):
        """Initialize value of for the CDAO tree object."""
        BaseTree.Tree.__init__(
            self, root=root or Clade(), rooted=rooted, id=id, name=name
        )
        self.weight = weight
        # a list of (predicate, object) pairs, containing additional triples
        # using this tree as subject
        self.attributes = []


class Clade(BaseTree.Clade):
    """CDAO Clade (sub-tree) object."""

    def __init__(
        self, branch_length=1.0, name=None, clades=None, confidence=None, comment=None
    ):
        """Initialize values for the CDAO Clade object."""
        BaseTree.Clade.__init__(
            self,
            branch_length=branch_length,
            name=name,
            clades=clades,
            confidence=confidence,
        )
        self.comment = comment
        # a list of (predicate, object) pairs, containing additional triples
        # using this clade as subject
        self.attributes = []
        self.tu_attributes = []
        self.edge_attributes = []
