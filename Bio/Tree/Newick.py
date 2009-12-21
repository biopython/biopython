# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Classes corresponding to Newick trees.

See classes in Bio.Nexus: Trees.Tree, Trees.NodeData, and Nodes.Chain.
"""
__docformat__ = "epytext en"

import warnings

import BaseTree

# Py2.4 compatibility
try:
    from functools import wraps
except ImportError:
    # From Python 2.5+ functools module
    def wraps(wrapped):
        def update_wrapper(wrapper):
            for attr in ('__module__', '__name__', '__doc__'):
                setattr(wrapper, attr, getattr(wrapped, attr))
            wrapper.__dict__.update(wrapped.__dict__)
            return wrapper
        return update_wrapper


def deprecated(hint):
    """Decorator for deprecated Nexus functions.

    'hint' is the recommended replacement function or method.
    The warning is triggered when the deprecated function (or property) is
    called, *not* when it is defined.
    """
    message = "use %s instead" % hint
    def deprecate(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            warnings.warn(message, DeprecationWarning, stacklevel=3)
            return func(*args, **kwargs)
        return wrapper
    return deprecate


# XXX from Bio.Nexus.Trees
# move to Utils?
def consensus(trees, threshold=0.5,outgroup=None):
    """Compute a majority rule consensus tree of all clades with relative
    frequency>=threshold from a list of trees.
    """

class _Shim(object):
    """Shim for compatibility with Bio.Nexus.Trees.
    """
    # Methods with deprecated arguments -- duplicated in Bio.Tree.BaseTree

    def is_terminal(self, node=None):
        """Returns True if all direct descendents are terminal."""
        # Deprecated Newick behavior
        if node is not None:
            warnings.warn("use node.is_terminal() method instead",
                          DeprecationWarning, stacklevel=2)
            return node.is_terminal()
        return (not self.clades)

    def is_preterminal(self, node=None):
        """Returns True if all direct descendents are terminal."""
        # Deprecated Newick behavior
        if node is not None:
            warnings.warn("use node.tree.is_preterminal() method instead",
                          DeprecationWarning, stacklevel=2)
            return node.tree.is_preterminal()
        if self.is_terminal():
            return False
        for clade in self.clades:
            if not clade.is_terminal():
                return False
        return True

    def count_terminals(self, node=None):
        if node is not None:
            warnings.warn("use node.tree.count_terminals() directly",
                          DeprecationWarning, stacklevel=2)
            return node.tree.count_terminals()
        counter = 0
        for i, leaf in enumerate(self.get_terminals()):
            counter = i
        return counter + 1

    # Deprecated methods from Bio.Nexus.Trees.Tree

    @deprecated("\"not node.is_terminal()\"")
    def is_internal(self, node):
        """Returns True if node is an internal node."""
        return not node.is_terminal()

    @deprecated("root.tree.distance(node)")
    def sum_branchlength(self, root, node):
        """Adds up the branchlengths from root (default self.root) to node."""
        return root.branch_length_to(node)

    @deprecated("node.tree.find_clades()")
    def get_taxa(self, node_id=None):
        """Return a list of all OTUs downwards from a node (self, node_id)."""
        if node_id is None:
            node_id = self
        return list(node_id.find_clades())

    @deprecated("node.tree.find(name=taxon)")
    def search_taxon(self, taxon):
        """Returns the first matching taxon in self.data.taxon.

        Not restricted to terminal nodes.

        node_id = search_taxon(self,taxon)
        """
        return self.find(name=taxon)



class Tree(BaseTree.Tree, _Shim):
    """Newick Tree object.
    """
    def __init__(self, root=None, rooted=False, id=None, name='', weight=1.0):
        BaseTree.Tree.__init__(self, root=root or Clade(),
                rooted=rooted, id=id, name=name)
        self.weight = weight

    # Ported from Bio.Nexus.Trees.Tree

    # TODO - port the rest of these methods to Tree or BaseTree.Tree
    # See unit tests

    # XXX from Nexus.Nodes.Chain

    def is_parent_of(self, parent, grandchild):
        """Check if grandchild is a subnode of parent."""
        # XXX direct descendent? or "parent.get_path(grandchild) is not None"?

    # XXX from Nexus.Trees.Tree
    # """Get information about trees (monphyly of taxon sets, congruence between
    # trees, common ancestors,...) and to manipulate trees (reroot trees, split
    # terminal nodes)."""

    def collapse_genera(self,space_equals_underscore=True):
        """Collapses all subtrees which belong to the same genus.

        (i.e share the same first word in their taxon name.
        """
        # XXX whoa! this sounds error-prone


    def split(self, parent_id=None, n=2, branchlength=1.0):
        """Speciation: generates n (default two) descendants of a node.

        [new ids] = split(self,parent_id=None,n=2,branchlength=1.0):
        """ 

    def prune(self, taxon):
        """Prunes a terminal taxon from the tree.

        If taxon is from a bifurcation, the connectiong node will be collapsed
        and its branchlength added to remaining terminal node. This might be no
        longer a meaningful value'

        @return previous node
        """

    def set_subtree(self, node):
        """Return subtree as a set of nested sets."""

    def is_identical(self,tree2):
        """Compare tree and tree2 for identity."""
        return self.set_subtree(self.root)==tree2.set_subtree(tree2.root)

    def is_compatible(self, tree2, threshold, strict=True):
        """Compares branches with support>threshold for compatibility."""

    def is_monophyletic(self, taxon_list):
        """Return node_id of common ancestor if taxon_list is monophyletic, -1 otherwise."""

    def is_bifurcating(self, node=None):
        """Return True if tree downstream of node is strictly bifurcating."""

    def branchlength2support(self):
        """Move values stored in data.branchlength to data.support, and set
        branchlength to 0.0

        This is necessary when support has been stored as branchlength (e.g.
        paup), and has thus been read in as branchlength. 
        """

    def convert_absolute_support(self, nrep):
        """Convert absolute support (clade-count) to rel. frequencies.

        Some software (e.g. PHYLIP consense) just calculate how often clades
        appear, instead of calculating relative frequencies.
        """

    def has_support(self, node=None):
        """Returns True if any of the nodes has data.support != None."""

    def randomize(self, ntax=None, taxon_list=None, branchlength=1.0, branchlength_sd=None, bifurcate=True):
        """Generates a random tree with ntax taxa and/or taxa from taxlabels.

        Trees are bifurcating by default. (Polytomies not yet supported).

        @return new tree
        """

    def display(self):
        """Quick and dirty lists of all nodes."""

    def unroot(self):
        """Define a unrooted Tree structure, using data of a rooted Tree."""

    def root_with_outgroup(self, outgroup=None):
        """???

        Hint:
            Hook subtree starting with node child to parent.
        """

    def merge_with_support(self, bstrees=None, constree=None, threshold=0.5, outgroup=None):
        """Merge clade support with phylogeny.

        From consensus or list of bootstrap-trees.

        tree=merge_bootstrap(phylo,bs_tree=<list_of_trees>)
        or
        tree=merge_bootstrap(phylo,consree=consensus_tree with clade support)
        """


class Clade(BaseTree.Subtree, _Shim):
    """Newick Clade (subtree) object.
    """
    def __init__(self, branch_length=1.0, name=None, clades=None,
            support=None, comment=None):
        BaseTree.Subtree.__init__(self, branch_length=branch_length,
                name=name, clades=clades)
        self.support = support
        self.comment = comment

    # Deprecated attributes from Bio.Nexus.Trees

    @property
    @deprecated('Clade.name')
    def taxon(self):
        return self.name

    @property
    @deprecated('Clade.name')
    def id(self):
        return self.name

    @property
    @deprecated('Clade.clades')
    def nodes(self):
        return self.clades

    @property
    @deprecated("the Clade object's attributes")
    def data(self):
        return _NodeData(taxon=self.name,
                        branchlength=self.branch_length,
                        support=self.support,
                        comment=self.comment)

class _NodeData:
    """Stores tree-relevant data associated with nodes (e.g. branches or OTUs).

    This exists only for backward compatibility with Bio.Nexus, and is
    deprecated.
    """
    def __init__(self, taxon, branchlength, support, comment):
        self.taxon = taxon
        self.branchlength = branchlength
        self.support = support
        self.comment = comment

