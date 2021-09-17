# Copyright (C) 2013 by Ben Morris (ben@bendmorris.com)
# Based on Bio.Nexus, copyright 2005-2008 by Frank Kauff & Cymon J. Cox
# and Bio.Phylo.Newick, copyright 2009 by Eric Talevich.
# All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""I/O function wrappers for the RDF/CDAO file format.

This is an RDF format that conforms to the Comparative Data Analysis Ontology (CDAO).
See: http://evolutionaryontology.org/cdao

This module requires the librdf Python bindings (http://www.librdf.org)

The CDAOIO.Parser, in addition to parsing text files, can also parse directly
from a triple store that implements the Redland storage interface; similarly,
the CDAOIO.Writer can store triples in a triple store instead of serializing
them to a file.
"""

import os

from io import StringIO

from Bio import MissingPythonDependencyError

from Bio.Phylo import CDAO
from ._cdao_owl import cdao_namespaces, resolve_uri

# import of cdao_elements from ._cdao_owl removed in Biopython 1.74


try:
    import rdflib

    rdfver = rdflib.__version__
    if rdfver[0] in ["1", "2"] or (rdfver in ["3.0.0", "3.1.0", "3.2.0"]):
        raise MissingPythonDependencyError(
            "Support for CDAO tree format requires RDFlib v3.2.1 or later."
        )
except ImportError:
    raise MissingPythonDependencyError(
        "Support for CDAO tree format requires RDFlib."
    ) from None

RDF_NAMESPACES = {
    "owl": "http://www.w3.org/2002/07/owl#",
    "rdf": "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
    "rdfs": "http://www.w3.org/2000/01/rdf-schema#",
}
RDF_NAMESPACES.update(cdao_namespaces)
# pad node ids with zeroes until they're at least this length
ZEROES = 8


def qUri(x):
    """Resolve URI for librdf."""
    return resolve_uri(x, namespaces=RDF_NAMESPACES)


def format_label(x):
    """Format label for librdf."""
    return x.replace("_", " ")


# ---------------------------------------------------------
# Public API


def parse(handle, **kwargs):
    """Iterate over the trees in a CDAO file handle.

    :returns: generator of Bio.Phylo.CDAO.Tree objects.

    """
    return Parser(handle).parse(**kwargs)


def write(trees, handle, plain=False, **kwargs):
    """Write a trees in CDAO format to the given file handle.

    :returns: number of trees written.

    """
    return Writer(trees).write(handle, plain=plain, **kwargs)


# ---------------------------------------------------------
# Input


class Parser:
    """Parse a CDAO tree given a file handle."""

    def __init__(self, handle=None):
        """Initialize CDAO tree parser."""
        self.handle = handle
        self.graph = None
        self.node_info = None
        self.children = {}
        self.rooted = False

    @classmethod
    def from_string(cls, treetext):
        """Instantiate the class from the given string."""
        handle = StringIO(treetext)
        return cls(handle)

    def parse(self, **kwargs):
        """Parse the text stream this object was initialized with."""
        self.parse_handle_to_graph(**kwargs)
        return self.parse_graph()

    def parse_handle_to_graph(
        self, rooted=False, parse_format="turtle", context=None, **kwargs
    ):
        """Parse self.handle into RDF model self.model."""
        if self.graph is None:
            self.graph = rdflib.Graph()
        graph = self.graph

        for k, v in RDF_NAMESPACES.items():
            graph.bind(k, v)

        self.rooted = rooted

        if "base_uri" in kwargs:
            base_uri = kwargs["base_uri"]
        else:
            # Windows style slashes cannot be used in an RDF URI
            base_uri = "file://" + os.path.abspath(self.handle.name).replace("\\", "/")

        graph.parse(file=self.handle, publicID=base_uri, format=parse_format)

        return self.parse_graph(graph, context=context)

    def parse_graph(self, graph=None, context=None):
        """Iterate over RDF model yielding CDAO.Tree instances."""
        if graph is None:
            graph = self.graph

        # look up branch lengths/TUs for all nodes
        self.get_node_info(graph, context=context)

        for root_node in self.tree_roots:
            clade = self.parse_children(root_node)

            yield CDAO.Tree(root=clade, rooted=self.rooted)

    def new_clade(self, node):
        """Return a CDAO.Clade object for a given named node."""
        result = self.node_info[node]

        kwargs = {}
        if "branch_length" in result:
            kwargs["branch_length"] = result["branch_length"]
        if "label" in result:
            kwargs["name"] = result["label"].replace("_", " ")
        if "confidence" in result:
            kwargs["confidence"] = result["confidence"]

        clade = CDAO.Clade(**kwargs)

        return clade

    def get_node_info(self, graph, context=None):
        """Create a dictionary containing information about all nodes in the tree."""
        self.node_info = {}
        self.obj_info = {}
        self.children = {}
        self.nodes = set()
        self.tree_roots = set()

        assignments = {
            qUri("cdao:has_Parent"): "parent",
            qUri("cdao:belongs_to_Edge_as_Child"): "edge",
            qUri("cdao:has_Annotation"): "annotation",
            qUri("cdao:has_Value"): "value",
            qUri("cdao:represents_TU"): "tu",
            qUri("rdfs:label"): "label",
            qUri("cdao:has_Support_Value"): "confidence",
        }

        for s, v, o in graph:
            # process each RDF triple in the graph sequentially

            s, v, o = str(s), str(v), str(o)

            if s not in self.obj_info:
                self.obj_info[s] = {}
            this = self.obj_info[s]

            try:
                # if the predicate is one we care about, store information for
                # later
                this[assignments[v]] = o
            except KeyError:
                pass

            if v == qUri("rdf:type"):
                if o in (qUri("cdao:AncestralNode"), qUri("cdao:TerminalNode")):
                    # this is a tree node; store it in set of all nodes
                    self.nodes.add(s)
            if v == qUri("cdao:has_Root"):
                # this is a tree; store its root in set of all tree roots
                self.tree_roots.add(o)

        for node in self.nodes:
            # for each node, look up all information needed to create a
            # CDAO.Clade
            self.node_info[node] = {}
            node_info = self.node_info[node]

            obj = self.obj_info[node]
            if "edge" in obj:
                # if this object points to an edge, we need a branch length from
                # the annotation on that edge
                edge = self.obj_info[obj["edge"]]
                if "annotation" in edge:
                    annotation = self.obj_info[edge["annotation"]]
                    if "value" in annotation:
                        node_info["branch_length"] = float(annotation["value"])

            if "tu" in obj:
                # if this object points to a TU, we need the label of that TU
                tu = self.obj_info[obj["tu"]]
                if "label" in tu:
                    node_info["label"] = tu["label"]

            if "parent" in obj:
                # store this node as a child of its parent, if it has one,
                # so that the tree can be traversed from parent to children
                parent = obj["parent"]
                if parent not in self.children:
                    self.children[parent] = []
                self.children[parent].append(node)

    def parse_children(self, node):
        """Traverse the tree to create a nested clade structure.

        Return a CDAO.Clade, and calls itself recursively for each child,
        traversing the entire tree and creating a nested structure of CDAO.Clade
        objects.
        """
        clade = self.new_clade(node)

        children = self.children[node] if node in self.children else []
        clade.clades = [self.parse_children(child_node) for child_node in children]

        return clade


# ---------------------------------------------------------
# Output


class Writer:
    """Based on the writer in Bio.Nexus.Trees (str, to_string)."""

    prefixes = RDF_NAMESPACES

    def __init__(self, trees):
        """Initialize parameters for writing a CDAO tree."""
        self.trees = trees

        self.node_counter = 0
        self.edge_counter = 0
        self.tu_counter = 0
        self.tree_counter = 0

    def write(
        self,
        handle,
        tree_uri="",
        record_complete_ancestry=False,
        rooted=False,
        **kwargs,
    ):
        """Write this instance's trees to a file handle."""
        self.rooted = rooted
        self.record_complete_ancestry = record_complete_ancestry

        if tree_uri and not tree_uri.endswith("/"):
            tree_uri += "/"

        trees = self.trees

        if tree_uri:
            handle.write(f"@base <{tree_uri}>\n")
        for k, v in self.prefixes.items():
            handle.write(f"@prefix {k}: <{v}> .\n")

        handle.write(f"<{self.prefixes['cdao']}> a owl:Ontology .\n")

        for tree in trees:
            self.tree_counter += 1
            self.tree_uri = "tree%s"

            first_clade = tree.clade
            statements = self.process_clade(first_clade, root=tree)
            for stmt in statements:
                self.add_stmt_to_handle(handle, stmt)

    def add_stmt_to_handle(self, handle, stmt):
        """Add URI prefix to handle."""
        # apply URI prefixes
        stmt_strings = []
        for n, part in enumerate(stmt):
            if isinstance(part, rdflib.URIRef):
                node_uri = str(part)
                changed = False
                for prefix, uri in self.prefixes.items():
                    if node_uri.startswith(uri):
                        node_uri = node_uri.replace(uri, f"{prefix}:", 1)
                        if node_uri == "rdf:type":
                            node_uri = "a"
                        changed = True
                if changed or ":" in node_uri:
                    stmt_strings.append(node_uri)
                else:
                    stmt_strings.append(f"<{node_uri}>")

            elif isinstance(part, rdflib.Literal):
                stmt_strings.append(part.n3())

            else:
                stmt_strings.append(str(part))

        handle.write(f"{' '.join(stmt_strings)} .\n")

    def process_clade(self, clade, parent=None, root=False):
        """Recursively generate triples describing a tree of clades."""
        self.node_counter += 1
        clade.uri = f"node{str(self.node_counter).zfill(ZEROES)}"
        if parent:
            clade.ancestors = parent.ancestors + [parent.uri]
        else:
            clade.ancestors = []

        def nUri(s):
            # nUri = lambda s: rdflib.URIRef(s)
            return rdflib.URIRef(s)

        def pUri(s):
            # pUri = lambda s: rdflib.URIRef(qUri(s))
            return rdflib.URIRef(qUri(s))

        tree_id = nUri("")

        statements = []

        if root is not False:
            # create a cdao:RootedTree with reference to the tree root
            tree_type = (
                pUri("cdao:RootedTree") if self.rooted else pUri("cdao:UnrootedTree")
            )

            statements += [
                (tree_id, pUri("rdf:type"), tree_type),
                (tree_id, pUri("cdao:has_Root"), nUri(clade.uri)),
            ]

            try:
                tree_attributes = root.attributes
            except AttributeError:
                tree_attributes = []

            for predicate, obj in tree_attributes:
                statements.append((tree_id, predicate, obj))

        if clade.name:
            # create TU
            self.tu_counter += 1
            tu_uri = f"tu{str(self.tu_counter).zfill(ZEROES)}"

            statements += [
                (nUri(tu_uri), pUri("rdf:type"), pUri("cdao:TU")),
                (nUri(clade.uri), pUri("cdao:represents_TU"), nUri(tu_uri)),
                (
                    nUri(tu_uri),
                    pUri("rdfs:label"),
                    rdflib.Literal(format_label(clade.name)),
                ),
            ]

            try:
                tu_attributes = clade.tu_attributes
            except AttributeError:
                tu_attributes = []

            for predicate, obj in tu_attributes:
                yield (nUri(tu_uri), predicate, obj)

        # create this node
        node_type = "cdao:TerminalNode" if clade.is_terminal() else "cdao:AncestralNode"
        statements += [
            (nUri(clade.uri), pUri("rdf:type"), pUri(node_type)),
            (nUri(clade.uri), pUri("cdao:belongs_to_Tree"), tree_id),
        ]

        if parent is not None:
            # create edge from the parent node to this node
            self.edge_counter += 1
            edge_uri = f"edge{str(self.edge_counter).zfill(ZEROES)}"

            statements += [
                (nUri(edge_uri), pUri("rdf:type"), pUri("cdao:DirectedEdge")),
                (nUri(edge_uri), pUri("cdao:belongs_to_Tree"), tree_id),
                (nUri(edge_uri), pUri("cdao:has_Parent_Node"), nUri(parent.uri)),
                (nUri(edge_uri), pUri("cdao:has_Child_Node"), nUri(clade.uri)),
                (
                    nUri(clade.uri),
                    pUri("cdao:belongs_to_Edge_as_Child"),
                    nUri(edge_uri),
                ),
                (nUri(clade.uri), pUri("cdao:has_Parent"), nUri(parent.uri)),
                (
                    nUri(parent.uri),
                    pUri("cdao:belongs_to_Edge_as_Parent"),
                    nUri(edge_uri),
                ),
            ]

            try:
                confidence = clade.confidence
            except AttributeError:
                pass
            else:
                if confidence is not None:
                    confidence = rdflib.Literal(
                        confidence, datatype="http://www.w3.org/2001/XMLSchema#decimal"
                    )

                    statements += [
                        (nUri(clade.uri), pUri("cdao:has_Support_Value"), confidence)
                    ]

            if self.record_complete_ancestry and len(clade.ancestors) > 0:
                statements += [
                    (nUri(clade.uri), pUri("cdao:has_Ancestor"), nUri(ancestor))
                    for ancestor in clade.ancestors
                ]

            if clade.branch_length is not None:
                # add branch length
                edge_ann_uri = f"edge_annotation{str(self.edge_counter).zfill(ZEROES)}"

                branch_length = rdflib.Literal(
                    clade.branch_length,
                    datatype=rdflib.URIRef("http://www.w3.org/2001/XMLSchema#decimal"),
                )
                statements += [
                    (nUri(edge_ann_uri), pUri("rdf:type"), pUri("cdao:EdgeLength")),
                    (nUri(edge_uri), pUri("cdao:has_Annotation"), nUri(edge_ann_uri)),
                    (nUri(edge_ann_uri), pUri("cdao:has_Value"), branch_length),
                ]

            try:
                edge_attributes = clade.edge_attributes
            except AttributeError:
                edge_attributes = []

            for predicate, obj in edge_attributes:
                yield (nUri(edge_uri), predicate, obj)

        yield from statements

        try:
            clade_attributes = clade.attributes
        except AttributeError:
            clade_attributes = []

        for predicate, obj in clade_attributes:
            yield (nUri(clade.uri), predicate, obj)

        if not clade.is_terminal():
            for new_clade in clade.clades:
                yield from self.process_clade(new_clade, parent=clade, root=False)
