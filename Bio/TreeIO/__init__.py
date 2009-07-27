# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""I/O function wrappers for phylogenetic tree formats.
"""
__docformat__ = "epytext en"

import NexusIO
import PhyloXMLIO

supported_formats = {
        'phyloxml': PhyloXMLIO,
        'nexus':    NexusIO,
        }

def read(file, format):
    """Parse a file in the given format and return a single object."""
    return getattr(supported_formats[format], 'read')(file)

def parse(file, format):
    """Iteratively parse a file and return each of the trees it contains.

    This is only supported for formats that can represent multiple phylogenetic
    trees in a single file.
    """
    return getattr(supported_formats[format], 'parse')(file)

def write(obj, file, format, **kwargs):
    """Serialize a Tree object into the given format and write to file."""
    return getattr(supported_formats[format], 'write')(file, **kwargs)
