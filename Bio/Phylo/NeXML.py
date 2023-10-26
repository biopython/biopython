# Copyright (C) 2013 Ben Morris (ben@bendmorris.com)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Classes corresponding to NeXML trees.

See classes in ``Bio.Nexus``: Trees.Tree, Trees.NodeData, and Nodes.Chain.
"""

from Bio.Phylo import BaseTree


class Tree(BaseTree.Tree):
    """NeXML Tree object."""

    def __init__(self, root=None, rooted=False, id=None, name=None, weight=1.0):
        """Instantiate a NeXML tree object with the given parameters."""
        BaseTree.Tree.__init__(
            self, root=root or Clade(), rooted=rooted, id=id, name=name
        )
        self.weight = weight


class Clade(BaseTree.Clade):
    """NeXML Clade (sub-tree) object."""

    def __init__(
        self,
        branch_length=1.0,
        name=None,
        clades=None,
        confidence=None,
        comment=None,
        **kwargs,
    ):
        """Initialize parameters for NeXML Clade object."""
        BaseTree.Clade.__init__(
            self,
            branch_length=branch_length,
            name=name,
            clades=clades,
            confidence=confidence,
        )
        self.comment = comment

        for key, value in kwargs.items():
            setattr(self, key, value)
