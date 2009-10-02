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


def deprecated(hint):
    """Decorator for deprecated Nexus functions.

    'hint' is the recommended replacement function or method.
    """
    message = "use %s instead" % hint
    def wrapper(func):
        warnings.warn(message, DeprecationWarning, stacklevel=3)
        return func
    return wrapper


# XXX from Bio.Nexus.Trees
# move to Utils?
def consensus(trees, threshold=0.5,outgroup=None):
    """Compute a majority rule consensus tree of all clades with relative frequency>=threshold from a list of trees."""


class NHTree(BaseTree.Tree):
    """Newick Tree object.
    """
    def __init__(self, root=None, clades=None, rooted=False, id=None, name='',
            weight=1.0,
            # values_are_support=False, max_support=1.0
            ):
        BaseTree.Tree.__init__(self,
                root=root or NHNode(), # originally the NodeData class
                clades=clades,            # list of NHNodes
                rooted=rooted,
                id=id,
                name=name)
        self.weight = weight
        # self.__values_are_support = values_are_support
        # self.max_support = max_support

    # Ported from Bio.Nexus.Trees.Tree

    def is_preterminal(self,node):
        """Returns True if all direct descendents are terminal."""
        # Py2.5+:
        # return (not self.is_terminal()) and all(t.is_terminal() for t in self)
        if self.is_terminal():
            return False
        for node in self.nodes:
            if not node.is_terminal():
                return False
        return True

    def count_terminals(self,node=None):
        """Counts the number of terminal nodes below this tree."""
        counter = 0
        for i, leaf in enumerate(self.get_leaves()):
            counter = i
        return counter + 1

    def collapse_genera(self,space_equals_underscore=True):
        """Collapses all subtrees which belong to the same genus.

        (i.e share the same first word in their taxon name.
        """

    # Deprecated methods from Bio.Nexus.Trees.Tree

    # XXX is_terminal in Nexus takes the node as an arg
    #   -- will have to monkeypatch there

    @deprecated("get_leaves")
    def get_terminals(self):
        """Return an iterable of all terminal nodes."""
        return self.get_leaves()

    @deprecated("\"not is_terminal\"")
    def is_internal(self,node):
        """Returns True if node is an internal node."""
        return not self.is_terminal()

    def sum_branchlength(self,root=None,node=None):
        """Adds up the branchlengths from root (default self.root) to node.
        
        sum = sum_branchlength(self,root=None,node=None)
        """

    # TODO - port the rest of these methods to NHTree or BaseTree
    # See unit tests

    # XXX from Nexus.Trees.Tree
    # """Get information about trees (monphyly of taxon sets, congruence between
    # trees, common ancestors,...) and to manipulate trees (reroot trees, split
    # terminal nodes)."""

    def split(self, parent_id=None, n=2, branchlength=1.0):
        """Speciation: generates n (default two) descendants of a node.

        [new ids] = split(self,parent_id=None,n=2,branchlength=1.0):
        """ 

    def search_taxon(self, taxon):
        """Returns the first matching taxon in self.data.taxon.

        Not restricted to terminal nodes.

        node_id = search_taxon(self,taxon)
        """

    def prune(self,taxon):
        """Prunes a terminal taxon from the tree.

        id_of_previous_node = prune(self,taxon)
        If taxon is from a bifurcation, the connectiong node will be collapsed
        and its branchlength added to remaining terminal node. This might be no
        longer a meaningful value'
        """

    def get_taxa(self,node_id=None):
        """Return a list of all OTUs downwards from a node (self, node_id).

        nodes = get_taxa(self,node_id=None)
        """


    def set_subtree(self,node):
        """Return subtree as a set of nested sets.

        sets = set_subtree(self,node)
        """

    def is_identical(self,tree2):
        """Compare tree and tree2 for identity.

        result = is_identical(self,tree2)
        """
        return self.set_subtree(self.root)==tree2.set_subtree(tree2.root)

    def is_compatible(self,tree2,threshold,strict=True):
        """Compares branches with support>threshold for compatibility.
        
        result = is_compatible(self,tree2,threshold)
        """

    def common_ancestor(self,node1,node2):
        """Return the common ancestor that connects two nodes.
        
        node_id = common_ancestor(self,node1,node2)
        """

    def distance(self,node1,node2):
        """Add and return the sum of the branchlengths between two nodes.
        dist = distance(self,node1,node2)
        """

    def is_monophyletic(self,taxon_list):
        """Return node_id of common ancestor if taxon_list is monophyletic, -1 otherwise.
        
        result = is_monophyletic(self,taxon_list)
        """

    def is_bifurcating(self,node=None):
        """Return True if tree downstream of node is strictly bifurcating."""

    def branchlength2support(self):
        """Move values stored in data.branchlength to data.support, and set branchlength to 0.0

        This is necessary when support has been stored as branchlength (e.g. paup), and has thus
        been read in as branchlength. 
        """

    def convert_absolute_support(self,nrep):
        """Convert absolute support (clade-count) to rel. frequencies.
        
        Some software (e.g. PHYLIP consense) just calculate how often clades appear, instead of
        calculating relative frequencies."""

    def has_support(self,node=None):
        """Returns True if any of the nodes has data.support != None."""

    def randomize(self,ntax=None,taxon_list=None,branchlength=1.0,branchlength_sd=None,bifurcate=True):
        """Generates a random tree with ntax taxa and/or taxa from taxlabels.
    
        new_tree = randomize(self,ntax=None,taxon_list=None,branchlength=1.0,branchlength_sd=None,bifurcate=True)
        Trees are bifurcating by default. (Polytomies not yet supported).
        """

    def display(self):
        """Quick and dirty lists of all nodes."""

    def to_string(self,support_as_branchlengths=False,branchlengths_only=False,plain=True,plain_newick=False,ladderize=None):
        """Return a paup compatible tree line.
       
        to_string(self,support_as_branchlengths=False,branchlengths_only=False,plain=True)
        """

    def unroot(self):
        """Defines a unrooted Tree structure, using data of a rooted Tree."""

    def root_with_outgroup(self,outgroup=None):
        """???

        Hint:
            Hook subtree starting with node child to parent.
        """

    def merge_with_support(self,bstrees=None,constree=None,threshold=0.5,outgroup=None):
        """Merges clade support (from consensus or list of bootstrap-trees) with phylogeny.

        tree=merge_bootstrap(phylo,bs_tree=<list_of_trees>)
        or
        tree=merge_bootstrap(phylo,consree=consensus_tree with clade support)
        """

    # XXX from Nexus.Nodes.Chain

    def collapse(self,id):
        """Deletes node from chain and relinks successors to predecessor: collapse(self, id)."""

    def is_parent_of(self,parent,grandchild):
        """Check if grandchild is a subnode of parent: is_parent_of(self,parent,grandchild)."""

    def trace(self,start,finish):
        """Returns a list of all node_ids between two nodes (excluding start, including end): trace(start,end)."""


class NHNode(BaseTree.Node):
    """Newick Node object.
    """
    def __init__(self, tree=None, label=None, branch_length=1.0,
            support=None, comment=None):
        BaseTree.Node.__init__(self,
                tree=tree or NHTree(), # a.k.a. taxon; self.tree.root == self
                label=label,
                branch_length=branch_length)
        self.support = support
        self.comment = comment

    # Deprecated attributes from Bio.Nexus.Trees

    @property
    @deprecated('Node.label')
    def id(self):
        return self.label

    @property
    @deprecated("the Node object's attributes")
    def data(self):
        return _NodeData(taxon=self.tree,
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

