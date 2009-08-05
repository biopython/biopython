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


def to_networkx(tree, graphviz=False):
    """Convert a Tree object to a networkx graph.

    The result is useful for plotting with pylab, matplotlib or pygraphviz,
    though the result is not quite a proper dendrogram for representing a
    phylogeny.

    Requires the networkx package.
    """
    # TODO: solve the graphviz labeling issue
    try:
        import networkx
    except ImportError:
        from Bio import MissingExternalDependencyError
        raise MissingExternalDependencyError(
                "The networkx library is not installed.")

    def get_label(node):
        """Get a unique hashable node value.

        Quirk: draw_graphviz doesn't honor the labels set in LabeledGraph, so
        str(node) must be unique for that engine. The regular nx.draw() doesn't
        have that problem.
        """
        if not graphviz:
            return node
        label = str(node)
        if label == node.__class__.__name__:
            return '<%d>' % id(node)
        return label

    def add_node(graph, node):
        if graphviz:
            graph.add_node(get_label(node))
        else:
            graph.add_node(node, get_label(node))

    def add_edge(graph, node1, node2):
        n1, n2 = map(get_label, [node1, node2])
        if node2.branch_length is not None:
            graph.add_edge(n1, n2, node2.branch_length)
        else:
            graph.add_edge(n1, n2)

    def build_subgraph(graph, top):
        """Walk down the Tree, building graphs, edges and nodes."""
        for node in top.nodes:
            add_node(graph, node)
            add_edge(graph, top, node)
            build_subgraph(graph, node)

    if tree.rooted:
        G = networkx.LabeledDiGraph()
    else:
        G = networkx.LabeledGraph()
    add_node(G, tree.root)
    build_subgraph(G, tree.root)
    return G


def draw_graphviz(tree):
    """Display a Tree object as a networkx graph, using the graphviz engine.

    Example:

        >>> import pylab
        >>> from Bio import Tree, TreeIO
        >>> tree = TreeIO.read('example.xml', 'phyloxml')
        >>> Tree.draw_graphviz(tree)
        >>> pylab.show()

    Requires networkx, matplotlib and pygraphviz.
    """
    try:
        import networkx
    except ImportError:
        from Bio import MissingExternalDependencyError
        raise MissingExternalDependencyError(
                "The networkx library is not installed.")

    G = to_networkx(tree, graphviz=True)
    networkx.draw_graphviz(G, prog='twopi')

