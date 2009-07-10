# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""I/O function wrappers for phylogenetic tree formats.
"""

from Bio import PhyloXML

def read(file, format):
    if format == 'phyloxml':
        return PhyloXML.read(file)

def write(obj, file, format, encoding=None):
    if format == 'phyloxml':
        return PhyloXML.write(obj, file, encoding)
