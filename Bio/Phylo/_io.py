# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""I/O function wrappers for phylogenetic tree formats.

This API follows the same semantics as Biopython's `SeqIO` and `AlignIO`.
"""
__docformat__ = "restructuredtext en"

from Bio.Phylo import BaseTree, NewickIO, NexusIO

# Python 2.4 doesn't have ElementTree, which PhyloXMLIO needs
try:
    from Bio.Phylo import PhyloXMLIO
except ImportError:
    # TODO: should we issue a warning? the installer will have already whined
    # raise MissingPythonDependencyError(
    #         "Install an ElementTree implementation if you want to use "
    #         "Bio.Phylo to parse phyloXML files.")
    supported_formats = {
            'newick':   NewickIO,
            'nexus':    NexusIO,
            }
else:
    supported_formats = {
            'newick':   NewickIO,
            'nexus':    NexusIO,
            'phyloxml': PhyloXMLIO,
            }


def parse(file, format):
    """Iteratively parse a file and return each of the trees it contains.

    If a file only contains one tree, this still returns an iterable object that
    contains one element.

    Example
    -------

    >>> trees = parse('../../Tests/PhyloXML/apaf.xml', 'phyloxml')
    >>> for tree in trees:
    ...     print tree.rooted
    True
    """
    do_close = False
    if isinstance(file, basestring):
        file = open(file, 'r')
        do_close = True
    # Py2.4 compatibility: this should be in a try/finally block
    # try:
    for tree in getattr(supported_formats[format], 'parse')(file):
        yield tree
    # finally:
    if do_close:
        file.close()


def read(file, format):
    """Parse a file in the given format and return a single tree.

    Raises a `ValueError` if there are zero or multiple trees -- if this occurs,
    use `parse` instead to get the complete sequence of trees.
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


def write(trees, file, format, **kwargs):
    """Write a sequence of trees to file in the given format."""
    if isinstance(trees, BaseTree.Tree):
        # Passed a single tree instead of an iterable -- that's OK
        trees = [trees]
    do_close = False
    if isinstance(file, basestring):
        file = open(file, 'w+')
        do_close = True
    try:
        n = getattr(supported_formats[format], 'write')(trees, file, **kwargs)
    finally:
        if do_close:
            file.close()
    return n


def convert(in_file, in_format, out_file, out_format, **kwargs):
    """Convert between two tree file formats."""
    trees = parse(in_file, in_format)
    return write(trees, out_file, out_format, **kwargs)
