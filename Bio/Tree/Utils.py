# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Utilities for handling and checking PhyloXML object trees.
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
    show = show_all and repr or str

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

