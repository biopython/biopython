# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""I/O function wrappers for phylogenetic tree formats.
"""

import NexusIO
import PhyloXMLIO

supported_formats = {
        'phyloxml': PhyloXMLIO,
        'nexus':    NexusIO,
        }

def read(file, format):
    return getattr(supported_formats[format], 'read')(file)

def parse(file, format):
    return getattr(supported_formats[format], 'parse')(file)

def write(obj, file, format, **kwargs):
    return getattr(supported_formats[format], 'write')(file, **kwargs)
