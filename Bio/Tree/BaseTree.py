# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Base classes for Bio.Tree objects.
"""
# More candidates (see BioSQL's PhyloDB):
#     TreeRoot (Term) -- is_alternate, significance
#     TreeQualifierValue (Term) -- value, rank
#     TreeDbxref (Dbxref)
#     Edge
#     EdgeQualifierValue (Term) -- value, rank
#     NodeQualifierValue (Term) -- value, rank
#     NodePath -- distance
#     NodeTaxon (Taxon) -- rank
#     NodeBioentry (Bioentry) -- rank
#     NodeDbxref (Dbxref)


class TreeElement(object):
    """Base class for all Bio.Tree classes."""
    pass


class Tree(TreeElement):
    # name, identifier, is_rooted

    def total_branch_length(self):
        """Get the total length of this tree (sum of all branch lengths)."""
        raise NotImplementedError


class Node(TreeElement):
    # label, left_idx, right_idx

    # From Bioperl's Bio::Tree::TreeI

    def get_leaf_nodes(self):
        """Request the taxa (leaves of the tree)."""
        raise NotImplementedError

    def get_root_node(self):
        """Get the root node of this tree."""
        raise NotImplementedError

    # From Bioperl's Bio::Tree::TreeFunctionsI

    # remove_node
    # get_lca (lowest common ancestor)
    # distance (between 2 nodes, specified however)
    # is_monophyletic
    # is_paraphyletic
    # reroot
