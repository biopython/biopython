# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Utilities for handling and checking PhyloXML object trees.

Third-party libraries are loaded when the corresponding function is called.
"""
__docformat__ = "epytext en"

import sys

import BaseTree


def pretty_print(treeobj, file=sys.stdout, show_all=False, indent=0):
    """Print a summary of the structure of a PhyloXML file.

    With the show_all option, also prints the primitive (native Python instead
    of PhyloXML) objects in the object tree.
    """
    assert isinstance(treeobj, BaseTree.TreeElement)
    if show_all:
        show = repr
    else:
        def show(obj):
            return '%s: %s' % (obj.__class__.__name__, obj)

    # Closing over file
    def print_indented(text, indent):
        """Write an indented string of text to file."""
        file.write('\t'*indent + text + '\n')

    def print_phylo(obj, indent):
        """Recursively print a PhyloElement object tree."""
        print_indented(show(obj), indent)
        indent += 1
        for attr in obj.__dict__:
            child = getattr(obj, attr)
            if isinstance(child, BaseTree.TreeElement):
                print_phylo(child, indent)
            elif isinstance(child, list):
                for elem in child:
                    if isinstance(elem, BaseTree.TreeElement):
                        print_phylo(elem, indent)

    print_phylo(treeobj, indent)


def to_networkx(tree):
    """Convert a Tree object to a networkx graph.

    The result is useful for plotting with pylab, matplotlib or pygraphviz,
    though the result is not quite a proper dendrogram for representing a
    phylogeny.

    Example:

        >>> import networkx, pylab
        >>> tree = TreeIO.read('example.xml', 'phyloxml')
        >>> G = to_networkx(tree)
        >>> networkx.draw_graphviz(G)
    """
    try:
        import networkx
    except ImportError:
        from Bio import MissingExternalDependencyError
        raise MissingExternalDependencyError(
                "The networkx library is not installed.")

    if tree.rooted:
        G = networkx.DiGraph()
    else:
        G = networkx.Graph()

    def get_label(node):
        label = str(node)
        if label == node.__class__.__name__:
            return '%s <%d>' % (label, id(node))
        return label

    # Walk down the Tree, building graphs, edges and nodes.
    def build_subgraph(top):
        rlabel = get_label(top)
        G.add_node(rlabel)
        for node in top.nodes:
            nlabel = get_label(node)
            G.add_node(nlabel)
            G.add_edge(rlabel, nlabel)
            build_subgraph(node)

    build_subgraph(tree.root)
    return G

