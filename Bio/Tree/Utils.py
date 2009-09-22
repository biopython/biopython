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

    The result is useful for graph-oriented analysis, and also interactive
    plotting with pylab, matplotlib or pygraphviz, though the result is not
    quite a proper dendrogram for representing a phylogeny.

    Requires the networkx-1.0 package.
    """
    try:
        import networkx
    except ImportError:
        from Bio import MissingExternalDependencyError
        raise MissingExternalDependencyError(
                "The networkx library is not installed.")

    def add_edge(graph, n1, n2):
        graph.add_edge(n1, n2, (n2.branch_length or 1.0))
        # ENH: add branch colors
        # if hasattr(n2, 'color') and n2.color is not None:
        #     graph[n1][n2]['color'] = n2.color

    def build_subgraph(graph, top):
        """Walk down the Tree, building graphs, edges and nodes."""
        for node in top.nodes:
            graph.add_node(node)
            add_edge(graph, top, node)
            build_subgraph(graph, node)

    if tree.rooted:
        G = networkx.DiGraph()
    else:
        G = networkx.Graph()
    G.add_node(tree.root)
    build_subgraph(G, tree.root)
    return G


def draw_graphviz(tree, prog='twopi', args='', node_color='#c0deff', **kwargs):
    """Display a Tree object as a networkx graph, using the graphviz engine.

    Requires NetworkX, matplotlib, Graphviz and either PyGraphviz or pydot.

    Example:

        >>> import pylab
        >>> from Bio import Tree, TreeIO
        >>> tree = TreeIO.read('example.xml', 'phyloxml')
        >>> Tree.draw_graphviz(tree)
        >>> pylab.show()
        >>> pylab.savefig('example.png')

    The second and third parameters apply to Graphviz, and the remaining
    arbitrary keyword arguments are passed directly to networkx.draw().  which
    in turn mostly wraps matplotlib/pylab.  See the documentation for Graphviz
    and networkx for detailed explanations.

    Graphviz parameters:

    @param prog: The graphviz program to use when rendering the graph (to file
        or screen). 'twopi' behaves the best for large graphs, reliably avoiding
        crossing edges, but for smaller graphs 'neato' can also look nice. 
        For small directed graphs, 'dot' may produce the most normal-looking
        phylogram, but is liable to cross edges in larger graphs. ('circo' and
        'fdp' are valid, but not recommended.)

    @param args: String of options passed to the external graphviz program.
        Normally not needed, but offered here for completeness.

    The NetworkX/matplotlib parameters are described in the docstrings for
    networkx.draw() and pylab.scatter(), but the most reasonable options to try
    are: 
        alpha, node_color, node_size, node_shape, edge_color, style,
        font_size, font_color, font_weight, font_family
    """
    try:
        import networkx
    except ImportError:
        from Bio import MissingExternalDependencyError
        raise MissingExternalDependencyError(
                "The networkx library is not installed.")

    G = to_networkx(tree)
    Gi = networkx.convert_node_labels_to_integers(G, discard_old_labels=False)

    try:
        posi = networkx.pygraphviz_layout(Gi, prog, args=args)
    except ImportError:
        try:
            posi = networkx.pydot_layout(Gi, prog)
        except ImportError:
            raise MissingExternalDependencyError(
                    "Neither PyGraphviz nor Pydot is installed.")

    posn = dict((node, posi[Gi.node_labels[node]]) for node in G)
    labels = dict((n, str(n)) for n in G.nodes()
                  if str(n) != n.__class__.__name__)
    networkx.draw(G, posn, nodelist=labels.keys(), labels=labels,
                  node_color=node_color, **kwargs)

