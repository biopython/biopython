# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Utilities for handling, displaying and exporting Phylo trees.

Third-party libraries are loaded when the corresponding function is called.
"""

import math
import sys

from Bio import MissingPythonDependencyError


def to_networkx(tree):
    """Convert a Tree object to a networkx graph.

    The result is useful for graph-oriented analysis, and also interactive
    plotting with pylab, matplotlib or pygraphviz, though the resulting diagram
    is usually not ideal for displaying a phylogeny.

    Requires NetworkX version 0.99 or later.
    """
    try:
        import networkx
    except ImportError:
        raise MissingPythonDependencyError(
            "Install NetworkX if you want to use to_networkx."
        ) from None

    # NB (1/2010): the networkx API stabilized at v.1.0
    # 1.0+: edges accept arbitrary data as kwargs, weights are floats
    # 0.99: edges accept weight as a string, nothing else
    # pre-0.99: edges accept no additional data
    # Ubuntu Lucid LTS uses v0.99, let's support everything
    if networkx.__version__ >= "1.0":

        def add_edge(graph, n1, n2):
            graph.add_edge(n1, n2, weight=n2.branch_length or 1.0)
            # Copy branch color value as hex, if available
            if hasattr(n2, "color") and n2.color is not None:
                graph[n1][n2]["color"] = n2.color.to_hex()
            elif hasattr(n1, "color") and n1.color is not None:
                # Cascading color attributes
                graph[n1][n2]["color"] = n1.color.to_hex()
                n2.color = n1.color
            # Copy branch weight value (float) if available
            if hasattr(n2, "width") and n2.width is not None:
                graph[n1][n2]["width"] = n2.width
            elif hasattr(n1, "width") and n1.width is not None:
                # Cascading width attributes
                graph[n1][n2]["width"] = n1.width
                n2.width = n1.width

    elif networkx.__version__ >= "0.99":

        def add_edge(graph, n1, n2):
            graph.add_edge(n1, n2, (n2.branch_length or 1.0))

    else:

        def add_edge(graph, n1, n2):
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


def draw_graphviz(
    tree, label_func=str, prog="twopi", args="", node_color="#c0deff", **kwargs
):
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
            are available::

                from Bio import Phylo, AlignIO
                from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
                constructor = DistanceTreeConstructor()
                aln = AlignIO.read(open('TreeConstruction/msa.phy'), 'phylip')
                calculator = DistanceCalculator('identity')
                dm = calculator.get_distance(aln)
                tree = constructor.upgma(dm)
                Phylo.draw_graphviz(tree, lambda n: n.taxonomies[0].code)

        prog : string
            The Graphviz program to use when rendering the graph. 'twopi'
            behaves the best for large graphs, reliably avoiding crossing edges,
            but for moderate graphs 'neato' looks a bit nicer.  For small
            directed graphs, 'dot' may produce a normal-looking cladogram, but
            will cross and distort edges in larger graphs. (The programs 'circo'
            and 'fdp' are not recommended.)
        args : string
            Options passed to the external graphviz program.  Normally not
            needed, but offered here for completeness.

    Examples
    --------
    Load a PhyloXML format tree, and draw a PNG using GraphViz::

        import pylab
        from Bio import Phylo
        tree = Phylo.read('PhyloXML/apaf.xml', 'phyloxml')
        Phylo.draw_graphviz(tree)
        pylab.show()
        pylab.savefig('apaf.png')

    """
    # Deprecated in Biopython 1.70 (#1247)
    import warnings
    from Bio import BiopythonDeprecationWarning

    warnings.warn(
        "draw_graphviz is deprecated; use Bio.Phylo.draw instead",
        BiopythonDeprecationWarning,
    )

    try:
        import networkx
    except ImportError:
        raise MissingPythonDependencyError(
            "Install NetworkX if you want to use to_networkx."
        ) from None

    G = to_networkx(tree)
    try:
        # NetworkX version 1.8 or later (2013-01-20)
        Gi = networkx.convert_node_labels_to_integers(G, label_attribute="label")
        int_labels = {}
        for integer, nodeattrs in Gi.node.items():
            int_labels[nodeattrs["label"]] = integer
    except TypeError:
        # Older NetworkX versions (before 1.8)
        Gi = networkx.convert_node_labels_to_integers(G, discard_old_labels=False)
        int_labels = Gi.node_labels

    try:
        try:
            # networkx versions before 1.11 (#1247)
            graphviz_layout = networkx.graphviz_layout
        except AttributeError:
            # networkx version 1.11
            graphviz_layout = networkx.drawing.nx_agraph.graphviz_layout
        posi = graphviz_layout(Gi, prog, args=args)
    except ImportError:
        raise MissingPythonDependencyError(
            "Install PyGraphviz or pydot if you want to use draw_graphviz."
        ) from None

    def get_label_mapping(G, selection):
        """Apply the user-specified node relabeling."""
        for node in G.nodes():
            if (selection is None) or (node in selection):
                try:
                    label = label_func(node)
                    if label not in (None, node.__class__.__name__):
                        yield (node, label)
                except (LookupError, AttributeError, ValueError):
                    pass

    if "nodelist" in kwargs:
        labels = dict(get_label_mapping(G, set(kwargs["nodelist"])))
    else:
        labels = dict(get_label_mapping(G, None))
    kwargs["nodelist"] = list(labels.keys())
    if "edge_color" not in kwargs:
        kwargs["edge_color"] = [
            isinstance(e[2], dict) and e[2].get("color", "k") or "k"
            for e in G.edges(data=True)
        ]
    if "width" not in kwargs:
        kwargs["width"] = [
            isinstance(e[2], dict) and e[2].get("width", 1.0) or 1.0
            for e in G.edges(data=True)
        ]

    posn = {n: posi[int_labels[n]] for n in G}
    networkx.draw(
        G, posn, labels=labels, with_labels=True, node_color=node_color, **kwargs
    )


def draw_ascii(tree, file=None, column_width=80):
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
            File handle opened for writing the output drawing. (Default:
            standard output)
        column_width : int
            Total number of text columns used by the drawing.

    """
    if file is None:
        file = sys.stdout

    taxa = tree.get_terminals()
    # Some constants for the drawing calculations
    max_label_width = max(len(str(taxon)) for taxon in taxa)
    drawing_width = column_width - max_label_width - 1
    drawing_height = 2 * len(taxa) - 1

    def get_col_positions(tree):
        """Create a mapping of each clade to its column position."""
        depths = tree.depths()
        # If there are no branch lengths, assume unit branch lengths
        if max(depths.values()) == 0:
            depths = tree.depths(unit_branch_lengths=True)
        # Potential drawing overflow due to rounding -- 1 char per tree layer
        fudge_margin = int(math.ceil(math.log(len(taxa), 2)))
        cols_per_branch_unit = (drawing_width - fudge_margin) / float(
            max(depths.values())
        )
        return {
            clade: int(blen * cols_per_branch_unit + 1.0)
            for clade, blen in depths.items()
        }

    def get_row_positions(tree):
        positions = {taxon: 2 * idx for idx, taxon in enumerate(taxa)}

        def calc_row(clade):
            for subclade in clade:
                if subclade not in positions:
                    calc_row(subclade)
            positions[clade] = (
                positions[clade.clades[0]] + positions[clade.clades[-1]]
            ) // 2

        calc_row(tree.root)
        return positions

    col_positions = get_col_positions(tree)
    row_positions = get_row_positions(tree)
    char_matrix = [[" " for x in range(drawing_width)] for y in range(drawing_height)]

    def draw_clade(clade, startcol):
        thiscol = col_positions[clade]
        thisrow = row_positions[clade]
        # Draw a horizontal line
        for col in range(startcol, thiscol):
            char_matrix[thisrow][col] = "_"
        if clade.clades:
            # Draw a vertical line
            toprow = row_positions[clade.clades[0]]
            botrow = row_positions[clade.clades[-1]]
            for row in range(toprow + 1, botrow + 1):
                char_matrix[row][thiscol] = "|"
            # NB: Short terminal branches need something to stop rstrip()
            if (col_positions[clade.clades[0]] - thiscol) < 2:
                char_matrix[toprow][thiscol] = ","
            # Draw descendents
            for child in clade:
                draw_clade(child, thiscol + 1)

    draw_clade(tree.root, 0)
    # Print the complete drawing
    for idx, row in enumerate(char_matrix):
        line = "".join(row).rstrip()
        # Add labels for terminal taxa in the right margin
        if idx % 2 == 0:
            line += " " + str(taxa[idx // 2])
        file.write(line + "\n")
    file.write("\n")


def build_ete3tree(tree, ete_root):
    """Convert a Tree object to a ete3's Tree."""
    for node in tree.root:
        ete_node = ete_root.add_child(name=node.name, dist=node.branch_length)
        ete_node.add_features(branch_length=node.branch_length)

        if hasattr(node, "confidences"):
            if len(node.confidences) == 1:
                confidence = node.confidences[0].value
                ete_node.add_features(confidence=confidence)
            else:
                ete_node.add_features(confidences=node.confidences)

        if hasattr(node, "taxonomy"):
            ete_node.add_features(taxonomy=node.taxonomy)

        build_ete3tree(node, ete_node)


def draw(
    tree,
    label_func=str,
    do_show=True,
    show_confidence=True,
    branch_labels=None,
    label_colors=None,
    **kwargs
):
    """Plot the given tree using ete3 or matplotlib (or pylab).

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
        branch_labels : dict or callable
            A mapping of each node to the label that will be shown along the
            branch leading to it. By default this is the confidence value(s) of
            the node, taken from the ``confidence`` attribute, and can be
            easily toggled off with this function's ``show_confidence`` option.
            But if you would like to alter the formatting of confidence values,
            or label the branches with something other than confidence, then use
            this option.
        label_colors : dict or callable
            A function or a dictionary specifying the color of the tip label.
            If the tip label can't be found in the dict or label_colors is
            None, the label will be shown in black.

    """
    try:
        import ete3
    except ImportError:
        raise MissingPythonDependencyError(
            "Install ete3 if you want to use the draw method."
        ) from None

    from ete3 import Tree as EteTree, TextFace

    # Create the ete3's tree

    ete_tree = EteTree(dist=0, support=0)
    ete_tree.add_features(taxonomy="", branch_length=0)
    build_ete3tree(tree, ete_tree)

    # Options for displaying branch labels / confidence

    if not branch_labels:
        if show_confidence:

            def format_branch_label(node):
                if hasattr(node, "confidences") and len(node.confidences) > 1:
                    return "/".join(str(cnf.value) for cnf in node.confidences)

                if hasattr(node, "confidence") and node.confidence is not None:
                    return str(node.confidence)

                return None

        else:

            def format_branch_label(node):
                return None

    elif isinstance(branch_labels, dict):

        def format_branch_label(node):
            return branch_labels.get(node)

    else:
        if not callable(branch_labels):
            raise TypeError(
                "branch_labels must be either a dict or a callable (function)"
            )
        format_branch_label = branch_labels

    # Options for displaying label colors.

    if label_colors:
        if callable(label_colors):

            def get_label_color(label):
                return label_colors(label)

        else:
            # label_colors is presumed to be a dict
            def get_label_color(label):
                return label_colors.get(label, "black")

    else:

        def get_label_color(label):
            # if label_colors is not specified, use black
            return "black"

    # Create the labels to add

    # In ete3 a TextFace must be created to color the leaf label
    # The show_leaf_name property must be set to False to avoid
    # duplicate labels (the TextFace one and the default label)

    leaf_faces = {}
    for node in ete_tree.traverse():
        # Branch labels
        branch_face = TextFace(format_branch_label(node) or "", fsize=8)
        node.add_face(branch_face, column=0, position="branch-bottom")

        # Ignore the default value of the __str__ method
        if label_func is str:
            label = node.name
        else:
            label = label_func(node)

        # Node labels
        node_face = TextFace(label, fgcolor=get_label_color(node), fsize=10)

        if node.is_leaf():
            leaf_faces[node] = node_face
        else:
            node.add_face(node_face, column=0)

    if do_show:
        try:
            from ete3 import TreeStyle

            ts = TreeStyle()
            ts.show_leaf_name = False

            # Parse and process key word arguments as treeStyle options
            show_leaves = True
            for key, value in kwargs.items():
                if not hasattr(ts, str(key)):
                    raise AttributeError(
                        "TreeStyle instance has no attribute {}".format(str(key))
                    )
                setattr(ts, str(key), value)

                if str(key) == "show_leaf_name":
                    show_leaves = ts.show_leaf_name
                    ts.show_leaf_name = False

            if show_leaves:
                for key, value in leaf_faces.items():
                    key.add_face(value, column=0)

            ete_tree.show(tree_style=ts)
        except ImportError:
            try:
                import matplotlib.pyplot as plt
            except ImportError:
                try:
                    import pylab as plt
                except ImportError:
                    raise MissingPythonDependencyError(
                        "Install qt or matplotlib or pylab if you want to use the draw method."
                    ) from None
            else:
                # Show using matplotlib / pylab
                import tempfile

                with tempfile.NamedTemporaryFile(mode="w", suffix=".png") as tmp:
                    ete_tree.render(tmp.name)
                    im = plt.imread(tmp.name)
                    plt.imshow(im)
                    plt.show()
                    tmp.flush()
