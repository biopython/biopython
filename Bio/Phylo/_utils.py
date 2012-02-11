# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Utilities for handling, displaying and exporting Phylo trees.

Third-party libraries are loaded when the corresponding function is called.
"""
__docformat__ = "restructuredtext en"

import math
import sys


def to_networkx(tree):
    """Convert a Tree object to a networkx graph.

    The result is useful for graph-oriented analysis, and also interactive
    plotting with pylab, matplotlib or pygraphviz, though the resulting diagram
    is usually not ideal for displaying a phylogeny.

    Requires NetworkX version 0.99 or 1.0.
    """
    try:
        import networkx
    except ImportError:
        from Bio import MissingPythonDependencyError
        raise MissingPythonDependencyError(
                "Install NetworkX if you want to use to_networkx.")

    def add_edge(graph, n1, n2):
        # NB (1/2010): the networkx API congealed recently
        # Ubuntu Lucid uses v0.99, newest is v1.0.1, let's support both
        if networkx.__version__ >= '1.0':
            graph.add_edge(n1, n2, weight=str(n2.branch_length or 1.0))
            # Copy branch color value as hex, if available
            if hasattr(n2, 'color') and n2.color is not None:
                graph[n1][n2]['color'] = n2.color.to_hex()
            elif hasattr(n1, 'color') and n1.color is not None:
                # Cascading color attributes
                graph[n1][n2]['color'] = n1.color.to_hex()
                n2.color = n1.color
            # Copy branch weight value (float) if available
            if hasattr(n2, 'width') and n2.width is not None:
                graph[n1][n2]['width'] = n2.width
            elif hasattr(n1, 'width') and n1.width is not None:
                # Cascading width attributes
                graph[n1][n2]['width'] = n1.width
                n2.width = n1.width
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
    build_subgraph(G, tree.root)
    return G


def draw_graphviz(tree, label_func=str, prog='twopi', args='',
        node_color='#c0deff', **kwargs):
    """Display a tree or clade as a graph, using the graphviz engine.

    Requires NetworkX, matplotlib, Graphviz and either PyGraphviz or pydot.

    The third and fourth parameters apply to Graphviz, and the remaining
    arbitrary keyword arguments are passed directly to networkx.draw(), which
    in turn mostly wraps matplotlib/pylab.  See the documentation for Graphviz
    and networkx for detailed explanations.

    The NetworkX/matplotlib parameters are described in the docstrings for
    networkx.draw() and pylab.scatter(), but the most reasonable options to try
    are: *alpha, node_color, node_size, node_shape, edge_color, style,
    font_size, font_color, font_weight, font_family*

    :Parameters:

        label_func : callable
            A function to extract a label from a node. By default this is str(),
            but you can use a different function to select another string
            associated with each node. If this function returns None for a node,
            no label will be shown for that node.

            The label will also be silently skipped if the throws an exception
            related to ordinary attribute access (LookupError, AttributeError,
            ValueError); all other exception types will still be raised. This
            means you can use a lambda expression that simply attempts to look
            up the desired value without checking if the intermediate attributes
            are available:

                >>> Phylo.draw_graphviz(tree, lambda n: n.taxonomies[0].code)

        prog : string
            The Graphviz program to use when rendering the graph. 'twopi'
            behaves the best for large graphs, reliably avoiding crossing edges,
            but for moderate graphs 'neato' looks a bit nicer.  For small
            directed graphs, 'dot' may produce the most normal-looking
            phylogram, but will cross and distort edges in larger graphs. (The
            programs 'circo' and 'fdp' are not recommended.)
        args : string
            Options passed to the external graphviz program.  Normally not
            needed, but offered here for completeness.

    Example
    -------

    >>> import pylab
    >>> from Bio import Phylo
    >>> tree = Phylo.read('ex/apaf.xml', 'phyloxml')
    >>> Phylo.draw_graphviz(tree)
    >>> pylab.show()
    >>> pylab.savefig('apaf.png')
    """
    try:
        import networkx
    except ImportError:
        from Bio import MissingPythonDependencyError
        raise MissingPythonDependencyError(
                "Install NetworkX if you want to use to_networkx.")

    G = to_networkx(tree)
    Gi = networkx.convert_node_labels_to_integers(G, discard_old_labels=False)
    try:
        posi = networkx.pygraphviz_layout(Gi, prog, args=args)
    except ImportError:
        try:
            posi = networkx.pydot_layout(Gi, prog)
        except ImportError:
            raise MissingPythonDependencyError(
                    "Install PyGraphviz or Pydot if you want to use "
                    "draw_graphviz.")
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
        kwargs['edge_color'] = [isinstance(e[2], dict) and
                                e[2].get('color', 'k') or 'k'
                                for e in G.edges(data=True)]
    if 'width' not in kwargs:
        kwargs['width'] = [isinstance(e[2], dict) and
                           e[2].get('width', 1.0) or 1.0
                           for e in G.edges(data=True)]
    networkx.draw(G, posn, labels=labels, node_color=node_color, **kwargs)


def draw_ascii(tree, file=sys.stdout, column_width=80):
    """Draw an ascii-art phylogram of the given tree.

    The printed result looks like::

                                        _________ Orange
                         ______________|
                        |              |______________ Tangerine
          ______________|
         |              |          _________________________ Grapefruit
        _|              |_________|
         |                        |______________ Pummelo
         |
         |__________________________________ Apple


    :Parameters:
        file : file-like object
            File handle opened for writing the output drawing.
        column_width : int
            Total number of text columns used by the drawing.
    """
    taxa = tree.get_terminals()
    # Some constants for the drawing calculations
    max_label_width = max(len(str(taxon)) for taxon in taxa)
    drawing_width = column_width - max_label_width - 1
    drawing_height = 2 * len(taxa) - 1

    def get_col_positions(tree):
        """Create a mapping of each clade to its column position."""
        depths = tree.depths()
        # If there are no branch lengths, assume unit branch lengths
        if not max(depths.itervalues()):
            depths = tree.depths(unit_branch_lengths=True)
        # Potential drawing overflow due to rounding -- 1 char per tree layer
        fudge_margin = int(math.ceil(math.log(len(taxa), 2)))
        cols_per_branch_unit = ((drawing_width - fudge_margin)
                                / float(max(depths.itervalues())))
        return dict((clade, int(round(blen*cols_per_branch_unit + 0.5)))
                    for clade, blen in depths.iteritems())

    def get_row_positions(tree):
        positions = dict((taxon, 2*idx) for idx, taxon in enumerate(taxa))
        def calc_row(clade):
            for subclade in clade:
                if subclade not in positions:
                    calc_row(subclade)
            positions[clade] = (positions[clade.clades[0]] +
                                positions[clade.clades[-1]]) / 2
        calc_row(tree.root)
        return positions

    col_positions = get_col_positions(tree)
    row_positions = get_row_positions(tree)
    char_matrix = [[' ' for x in range(drawing_width)]
                    for y in range(drawing_height)]

    def draw_clade(clade, startcol):
        thiscol = col_positions[clade]
        thisrow = row_positions[clade]
        # Draw a horizontal line
        for col in range(startcol, thiscol):
            char_matrix[thisrow][col] = '_'
        if clade.clades:
            # Draw a vertical line
            toprow = row_positions[clade.clades[0]]
            botrow = row_positions[clade.clades[-1]]
            for row in range(toprow+1, botrow+1):
                char_matrix[row][thiscol] = '|'
            # NB: Short terminal branches need something to stop rstrip()
            if (col_positions[clade.clades[0]] - thiscol) < 2:
                char_matrix[toprow][thiscol] = ','
            # Draw descendents
            for child in clade:
                draw_clade(child, thiscol+1)

    draw_clade(tree.root, 0)
    # Print the complete drawing
    for idx, row in enumerate(char_matrix):
        line = ''.join(row).rstrip()
        # Add labels for terminal taxa in the right margin
        if idx % 2 == 0:
            line += ' ' + str(taxa[idx/2])
        file.write(line + '\n')
    file.write('\n')


def draw(tree, label_func=str, do_show=True, show_confidence=True, axes=None):
    """Plot the given tree using matplotlib (or pylab).

    The graphic is a rooted tree, drawn with roughly the same algorithm as
    draw_ascii.

    :Parameters:
        label_func : callable
            A function to extract a label from a node. By default this is str(),
            but you can use a different function to select another string
            associated with each node. If this function returns None for a node,
            no label will be shown for that node.
        do_show : bool
            Whether to show() the plot automatically.
        show_confidence : bool
            Whether to display confidence values, if present on the tree.
        axes : matplotlib/pylab axes
            If a valid matplotlib.axes.Axes instance, the phylogram is plotted
            in that Axes. By default (None), a new figure is created.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        try:
            import pylab as plt
        except ImportError:
            from Bio import MissingPythonDependencyError
            raise MissingPythonDependencyError(
                    "Install matplotlib or pylab if you want to use draw.")

    def get_x_positions(tree):
        """Create a mapping of each clade to its horizontal position.

        Dict of {clade: x-coord}
        """
        depths = tree.depths()
        # If there are no branch lengths, assume unit branch lengths
        if not max(depths.itervalues()):
            depths = tree.depths(unit_branch_lengths=True)
        return depths

    def get_y_positions(tree):
        """Create a mapping of each clade to its vertical position.

        Dict of {clade: y-coord}.
        Coordinates are negative, and integers for tips.
        """
        maxheight = tree.count_terminals()
        # Rows are defined by the tips
        heights = dict((tip, maxheight - i)
                for i, tip in enumerate(reversed(tree.get_terminals())))
        # Internal nodes: place at midpoint of children
        def calc_row(clade):
            for subclade in clade:
                if subclade not in heights:
                    calc_row(subclade)
            # Closure over heights
            heights[clade] = (heights[clade.clades[0]] +
                                heights[clade.clades[-1]]) / 2.0

        if tree.root.clades:
            calc_row(tree.root)
        return heights

    x_posns = get_x_positions(tree)
    y_posns = get_y_positions(tree)
    # The function draw_clade closes over the axes object
    if axes is None:
        fig = plt.figure()
        axes = fig.add_subplot(1,1,1)
    elif not isinstance(axes,plt.matplotlib.axes.Axes):
        raise ValueError("Invalid argument for axes: %s" % axes)

    def draw_clade(clade, x_start, color='k', lw=1):
        """Recursively draw a tree, down from the given clade."""
        x_here = x_posns[clade]
        y_here = y_posns[clade]
        # phyloXML-only graphics annotations
        if hasattr(clade, 'color') and clade.color is not None:
            color = clade.color.to_hex()
        if hasattr(clade, 'width') and clade.width is not None:
            lw = clade.width
        # Draw a horizontal line from start to here
        axes.hlines(y_here, x_start, x_here, color=color, lw=lw)
        # Add node/taxon labels
        label = label_func(clade)
        if label not in (None, clade.__class__.__name__):
            axes.text(x_here, y_here, ' ' + label,
                    fontsize=10, verticalalignment='center')
        # Add confidence
        if hasattr(clade, 'confidences'):
            # phyloXML supports multiple confidences
            conf_label = '/'.join(map(str, map(float, clade.confidences)))
        elif clade.confidence is not None:
            conf_label = str(clade.confidence)
        else:
            conf_label = None
        if conf_label and show_confidence:
            plt.text(x_start, y_here, conf_label, fontsize=9)
        if clade.clades:
            # Draw a vertical line connecting all children
            y_top = y_posns[clade.clades[0]]
            y_bot = y_posns[clade.clades[-1]]
            # Only apply widths to horizontal lines, like Archaeopteryx
            axes.vlines(x_here, y_bot, y_top, color=color)
            # Draw descendents
            for child in clade:
                draw_clade(child, x_here, color, lw)

    draw_clade(tree.root, 0)
    if hasattr(tree, 'name') and tree.name:
        axes.set_title(tree.name)
    axes.set_xlabel('branch length')
    axes.set_ylabel('taxa')
    # Add margins around the tree to prevent overlapping the axes
    xmin, xmax = axes.get_xlim()
    pad = 0.05 * xmax
    axes.set_xlim(-pad, xmax + pad)
    # Also invert the y-axis (origin at the top)
    # Add a small vertical margin, but avoid including 0 and N+1 on the y axis
    axes.set_ylim(max(y_posns.itervalues()) + 0.8, 0.2)
    if do_show:
        plt.show()

