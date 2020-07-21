# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""I/O function wrappers for phylogenetic tree formats.

This API follows the same semantics as Biopython's ``SeqIO`` and
``AlignIO``.
"""


from Bio import File
from Bio.Phylo import BaseTree, NewickIO, NexusIO, PhyloXMLIO, NeXMLIO

supported_formats = {
    "newick": NewickIO,
    "nexus": NexusIO,
    "phyloxml": PhyloXMLIO,
    "nexml": NeXMLIO,
}

try:
    from Bio.Phylo import CDAOIO

    supported_formats["cdao"] = CDAOIO
except ImportError:
    pass


def parse(file, format, **kwargs):
    """Parse a file iteratively, and yield each of the trees it contains.

    If a file only contains one tree, this still returns an iterable object that
    contains one element.

    Examples
    --------
    >>> import Bio.Phylo
    >>> trees = Bio.Phylo.parse('PhyloXML/apaf.xml', 'phyloxml')
    >>> for tree in trees:
    ...     print(tree.rooted)
    True

    """
    with File.as_handle(file) as fp:
        yield from getattr(supported_formats[format], "parse")(fp, **kwargs)


def read(file, format, **kwargs):
    """Parse a file in the given format and return a single tree.

    Raises a ``ValueError`` if there are zero or multiple trees -- if this
    occurs, use ``parse`` instead to get the complete sequence of trees.
    """
    try:
        tree_gen = parse(file, format, **kwargs)
        tree = next(tree_gen)
    except StopIteration:
        raise ValueError("There are no trees in this file.") from None
    try:
        next(tree_gen)
    except StopIteration:
        return tree
    else:
        raise ValueError("There are multiple trees in this file; use parse() instead.")


def write(trees, file, format, **kwargs):
    """Write a sequence of trees to file in the given format."""
    if isinstance(trees, (BaseTree.Tree, BaseTree.Clade)):
        # Passed a single tree instead of an iterable -- that's OK
        trees = [trees]
    with File.as_handle(file, "w+") as fp:
        n = getattr(supported_formats[format], "write")(trees, fp, **kwargs)
    return n


def convert(in_file, in_format, out_file, out_format, parse_args=None, **kwargs):
    """Convert between two tree file formats."""
    if parse_args is None:
        parse_args = {}
    trees = parse(in_file, in_format, **parse_args)
    return write(trees, out_file, out_format, **kwargs)
