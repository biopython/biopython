# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""I/O function wrappers for phylogenetic tree formats.
"""
__docformat__ = "epytext en"

import NewickIO
import NexusIO
import PhyloXMLIO

supported_formats = {
        'newick':   NewickIO,
        'nexus':    NexusIO,
        'phyloxml': PhyloXMLIO,
        }


def parse(file, format):
    """Iteratively parse a file and return each of the trees it contains.

    This is only supported for formats that can represent multiple phylogenetic
    trees in a single file.
    """
    return getattr(supported_formats[format], 'parse')(file)


def read(file, format):
    """Parse a file in the given format and return a single tree.

    Raises a ValueError if there are zero or multiple trees -- if this occurs,
    use parse() instead to get the complete sequence of trees.
    """
    try:
        tree_gen = parse(file, format)
        tree = tree_gen.next()
    except StopIteration:
        raise ValueError("There are no trees in this file.")
    try:
        tree_gen.next()
    except StopIteration:
        return tree
    else:
        raise ValueError(
                "There are multiple trees in this file; use parse() instead.")

    return getattr(supported_formats[format], 'read')(file)


def write(trees, file, format, **kwargs):
    """Write a sequence of trees to file in the given format."""
    if not hasattr(trees, '__iter__'):
        # Probably passed a single tree instead of a sequence -- that's OK
        trees = [trees]
    return getattr(supported_formats[format], 'write')(trees, file, **kwargs)


def convert(in_file, in_format, out_file, out_format, **kwargs):
    """Convert between two tree file formats."""
    trees = parse(in_file, in_format)
    write(trees, out_file, out_format, **kwargs)
