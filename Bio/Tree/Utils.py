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
    assert isinstance(treeobj, BaseTree.TreeElement), \
            "%s is not a valid TreeElement" % repr(treeobj)
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

    Requires NetworkX version 0.99 or 1.0.
    """
    try:
        import networkx
    except ImportError:
        from Bio import MissingExternalDependencyError
        raise MissingExternalDependencyError(
                "The networkx library is not installed.")

    def add_edge(graph, n1, n2):
        # NB (9/2009): the networkx API is hella unstable
        # Ubuntu Karmic uses 0.99, newest is 1.0rc1, let's support both
        if networkx.__version__ >= '1.0':
            graph.add_edge(n1, n2, weight=str(n2.branch_length or 1.0))
            if hasattr(n2, 'color') and n2.color is not None:
                graph[n1][n2]['color'] = n2.color.to_hex()
        elif networkx.__version__ >= '0.99':
            graph.add_edge(n1, n2, (n2.branch_length or 1.0))
        else:
            graph.add_edge(n1, n2)

    def build_subgraph(graph, top):
        """Walk down the Tree, building graphs, edges and nodes."""
        for clade in top:
            graph.add_node(clade.root)
            add_edge(graph, top.root, clade.root)
            build_subgraph(graph, clade)

    if tree.rooted:
        G = networkx.DiGraph()
    else:
        G = networkx.Graph()
    G.add_node(tree.root)
    build_subgraph(G, tree)
    return G


def draw_graphviz(tree, label_func=str, prog='neato', args='',
        node_color='#c0deff', **kwargs):
    """Display a Tree object as a dendrogram, using the graphviz engine.

    Requires NetworkX, matplotlib, Graphviz and either PyGraphviz or pydot.

    Example:

        >>> import pylab
        >>> from Bio import Tree, TreeIO
        >>> tree = TreeIO.read('example.xml', 'phyloxml')
        >>> Tree.draw_graphviz(tree)
        >>> pylab.show()
        >>> pylab.savefig('example.png')

    @param label_func: A function to extract a label from a node. By default
        this is str(), but you can use a different function to select another
        string associated with each node. If this function returns None for a
        node, no label will be shown for that node.

        The label will also be silently skipped if the throws an exception
        related to ordinary attribute access (LookupError, AttributeError,
        ValueError); all other exception types will still be raised. This
        means you can use a lambda expression that simply attempts to look up
        the desired value without checking if the intermediate attributes are
        available:

        >>> Tree.draw_graphviz(tree, lambda n: n.taxonomies[0].code)

    The third and fourth parameters apply to Graphviz, and the remaining
    arbitrary keyword arguments are passed directly to networkx.draw(), which
    in turn mostly wraps matplotlib/pylab.  See the documentation for Graphviz
    and networkx for detailed explanations.

    Graphviz parameters:

    @param prog: The Graphviz program to use when rendering the graph. 'twopi'
        behaves the best for large graphs, reliably avoiding crossing edges, but
        for moderate graphs 'neato' looks a bit nicer.  For small directed
        graphs, 'dot' may produce the most normal-looking phylogram, but will
        cross and distort edges in larger graphs. (The programs 'circo' and
        'fdp' are not recommended.)

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
    posn = dict((n, posi[Gi.node_labels[n]]) for n in G)

    def get_label_mapping(G, selection):
        for node in G.nodes():
            if (selection is None) or (node in selection):
                try:
                    label = label_func(node)
                    if label not in (None, node.__class__.__name__):
                        yield (node, label)
                except (LookupError, AttributeError, ValueError):
                    pass

    if 'nodelist' in kwargs:
        labels = dict(get_label_mapping(G, set(kwargs['nodelist'])))
    else:
        labels = dict(get_label_mapping(G, None))
    kwargs['nodelist'] = labels.keys()
    if 'edge_color' not in kwargs:
        kwargs['edge_color'] = [isinstance(e[2], dict)
                                and e[2].get('color', 'k') or 'k'
                                for e in G.edges(data=True)]
    networkx.draw(G, posn, labels=labels, node_color=node_color, **kwargs)
