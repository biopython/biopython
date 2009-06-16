# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Utilities for handling and checking PhyloXML object trees.

"""

import sys
import warnings

from Bio.PhyloXML import Parser
from Bio.PhyloXML.Parser import ElementTree


def dump_tags(handle, output=sys.stdout):
    """Extract tags from an XML document, writing them to stdout by default.

    This utility is meant for testing and debugging.
    """
    for event, elem in ElementTree.iterparse(handle, events=('start', 'end')):
        if event == 'start':
            output.write(elem.tag + '\n')
        else:
            elem.clear()


def pretty_print(source, indent=0, show_all=False, output=sys.stdout):
    """Print a summary of the structure of a PhyloXML file.

    With the show_all option, also prints the primitive (native Python instead
    of PhyloXML) objects in the object tree.
    """
    # Closing over output
    def print_indented(text, indent):
        """Write an indented string of text to output."""
        output.write("%s%s\n" %('\t'*indent, text))

    # Closing over show_all
    def print_phylo(obj, indent):
        """Recursively print a PhyloElement object tree."""
        print_indented(obj.__class__.__name__, indent)
        indent += 1
        if show_all:
            simple_objs = []
        for attr in obj.__dict__:
            child = getattr(obj, attr)
            if isinstance(child, Parser.PhyloElement):
                print_indented(child.__class__.__name__, indent)
                print_phylo(child, indent)
            elif show_all:
                simple_objs.append(attr)
        if show_all and simple_objs:
            print_indented(', '.join(simple_objs), indent)

    def print_other(obj, indent):
        """Recursively print a tree of Other objects."""
        if show_all:
            print_indented("%s -- %s" % (obj, obj.attributes.keys()), indent)
        else:
            print_indented(obj, indent)
        indent += 1
        for child in obj.children:
            if isinstance(child, Parser.Other):
                print_other(child, indent)
            else:
                print_indented(child, indent)

    if isinstance(source, Parser.Phylogeny):
        print_phylo(source)
        return

    if isinstance(source, Parser.Phyloxml):
        phyloxml = source
    else:
        phyloxml = Parser.read(source)
    print_indented(phyloxml.__class__.__name__, indent)
    indent += 1
    for tree in phyloxml.phylogenies:
        print_phylo(tree, indent)
    for otr in phyloxml.other:
        print_other(otr, indent)

