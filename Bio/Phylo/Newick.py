# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Classes corresponding to Newick trees, also used for Nexus trees.

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


class _TreeShim(object):
    """Shim for compatibility with the Bio.Nexus.Trees.Tree class.

    This class and its use in Tree (below) can eventually be deleted.
    """
    # * Methods with deprecated arguments -- duplicated in Bio.Phylo.BaseTree *
    # Each of these method checks for usage of deprecated arguments, issues a
    # warning if so, then performs the usual BaseTree procedure, adjusting
    # arguments as needed.

    # from Bio.Nexus.Nodes.Chain
    def is_parent_of(self, target, child=None):
        """Check if grandchild is a subnode of parent."""
        if child is not None:
            warnings.warn("use parent.is_parent_of(child) directly",
                          DeprecationWarning, stacklevel=2)
            return (target.get_path(child) is not None)
        return (self.get_path(target) is not None)

    # from Bio.Nexus.Trees.Tree
    def count_terminals(self, node=None):
        if node is not None:
            warnings.warn("use node.count_terminals() directly",
                          DeprecationWarning, stacklevel=2)
            return node.tree.count_terminals()
        counter = 0
        for i, leaf in enumerate(self.get_terminals()):
            counter = i
        return counter + 1

    def is_bifurcating(self, node=None):
        """True if tree downstream of node is strictly bifurcating."""
        if node is not None:
            warnings.warn("use node.is_bifurcating() directly instead",
                          DeprecationWarning, stacklevel=2)
            return node.is_bifurcating()
        # Root can be trifurcating, because it has no ancestor
        if isinstance(self, BaseTree.Tree) and len(self.root) == 3:
            return (self.clade[0].is_bifurcating()
                    and self.clade[1].is_bifurcating()
                    and self.clade[2].is_bifurcating())
        if len(self.root) == 2:
            return (self.clade[0].is_bifurcating()
                    and self.clade[1].is_bifurcating())
        if len(self.root) == 0:
            return True
        return False

    def is_monophyletic(self, taxon_list):
        """Return node_id of common ancestor if taxon_list is monophyletic, -1 otherwise."""
        # Validation
        if isinstance(taxon_list, basestring):
            warnings.warn("argument should be a list, not a string",
                          DeprecationWarning, stacklevel=2)
            target_set = set([taxon_list])
        else:
            target_set = set(taxon_list)
        # Try narrower subclades until a complete match or mismatch is found
        current = self.root
        while True:
            if set(current.get_terminals()) == target_set:
                return current
            for subclade in current.clades:
                if set(subclade.get_terminals()).issuperset(target_set):
                    current = subclade
                    break
                else:
                    return False

    def is_preterminal(self, node=None):
        """Returns True if all direct descendents are terminal."""
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

    def is_terminal(self, node=None):
        """Returns True if all direct descendents are terminal."""
        # Deprecated Newick behavior
        if node is not None:
            warnings.warn("use node.is_terminal() method instead",
                          DeprecationWarning, stacklevel=2)
            return node.is_terminal()
        return (not self.clades)

    def split(self,
            parent_id=None, # deprecated
            n=2,
            branchlength=None, # deprecated
            branch_length=1.0):
        """Speciation: generates n (default two) descendants from parent."""
        # Warn deprecated arguments
        if branchlength is not None:
            warnings.warn("use branch_length argument instead of branchlength",
                          DeprecationWarning, stacklevel=2)
            branch_length = branchlength
        if parent_id is not None:
            warnings.warn("use parent_id's split() method directly",
                          DeprecationWarning, stacklevel=2)
            parent = parent_id
        else:
            parent = self.root
        for i in range(n):
            node = Clade(name=parent.name+str(i),
                         branch_length=branch_length)
            parent.clades.append(node)
        return parent.clades[-n:]

    # * Deprecated methods from Bio.Nexus.Trees.Tree *
    # It is not deemed necessary to implement these in BaseTree.

    @deprecated("\"not node.is_terminal()\"")
    def is_internal(self, node):
        """Returns True if node is an internal node."""
        return not node.is_terminal()

    @deprecated("root.distance(node)")
    def sum_branchlength(self, root, node):
        """Adds up the branchlengths from root (default self.root) to node."""
        return root.distance(node)

    @deprecated("node.find_clades()")
    def get_taxa(self, node_id=None):
        """Return a list of all OTUs downwards from a node (self, node_id)."""
        if node_id is None:
            node_id = self.root
        return list(node_id.find_clades())

    @deprecated('self.randomized()')
    def randomize(self, ntax=None, taxon_list=None, branchlength=1.0, branchlength_sd=None, bifurcate=True):
        """Generates a random tree with ntax taxa and/or taxa from taxlabels.

        Trees are bifurcating regardless of the value of the bifurcate argument
        (polytomies not yet supported).
        """
        if not ntax and taxon_list:
            ntax=len(taxon_list)
        elif not taxon_list and ntax:
            taxon_list=['taxon'+str(i+1) for i in range(ntax)]
        elif not ntax and not taxon_list:
            raise ValueError('Either numer of taxa or list of taxa must be specified.')
        elif ntax != len(taxon_list):
            raise ValueError('Length of taxon list must correspond to ntax.')
        self = self.randomized(taxa=taxon_list, branch_length=branchlength,
                               branch_stdev=branchlength_sd)
        return

    @deprecated("node.find(name=taxon)")
    def search_taxon(self, taxon):
        """Returns the first matching taxon in this tree.

        Not restricted to terminal nodes.
        """
        return self.find(name=taxon)

    # --- TODO: port these --- 
    # """Get information about trees (monphyly of taxon sets, congruence between
    # trees, common ancestors,...) and to manipulate trees (reroot trees, split
    # terminal nodes)."""

    def collapse_genera(self,space_equals_underscore=True):
        """Collapses all subtrees which belong to the same genus.

        (i.e share the same first word in their taxon name.
        """
        # XXX whoa! this sounds error-prone

    def set_subtree(self, node):
        """Return subtree as a set of nested sets."""

    def is_identical(self,tree2):
        """Compare tree and tree2 for identity."""
        return self.set_subtree(self.root)==tree2.set_subtree(tree2.root)

    def is_compatible(self, tree2, threshold, strict=True):
        """Compares branches with support>threshold for compatibility."""

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



class _NodeShim(object):
    """Shim for compatibility with the Bio.Nexus.Nodes.Chain class.
    """

    @property
    @deprecated('Clade.name')
    def id(self):
        return self.name

    @property
    @deprecated("the Clade object's attributes")
    def data(self):
        return _NodeData(taxon=self.name,
                        branchlength=self.branch_length,
                        support=self.support,
                        comment=self.comment)

    # TODO?
    # prev
    # succ

class _NodeData(object):
    """Stores tree-relevant data associated with nodes (e.g. branches or OTUs).

    This exists only for backward compatibility with Bio.Nexus, and is
    deprecated.
    """
    def __init__(self, taxon, branchlength, support, comment):
        self.taxon = taxon
        self.branchlength = branchlength
        self.support = support
        self.comment = comment

# /end of shims


class Tree(BaseTree.Tree, _TreeShim):
    """Newick Tree object."""

    def __init__(self, root=None, rooted=False, id=None, name='', weight=1.0):
        BaseTree.Tree.__init__(self, root=root or Clade(),
                rooted=rooted, id=id, name=name)
        self.weight = weight


class Clade(BaseTree.Subtree, _NodeShim):
    """Newick Clade (subtree) object."""

    def __init__(self, branch_length=1.0, name=None, clades=None,
            support=None, comment=None):
        BaseTree.Subtree.__init__(self, branch_length=branch_length,
                name=name, clades=clades)
        self.support = support
        self.comment = comment

