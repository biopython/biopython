# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""PhyloXML writer and associated functions.

Constructs an XML file from a Tree.PhyloXML object.
"""

import warnings

import Tree
from Exceptions import PhyloXMLError, PhyloXMLWarning
from Parser import ElementTree, NAMESPACES


def write(phyloxml, file, encoding='utf-8'):
    """Write a phyloXML file.

    The file argument can be either an open handle or a file name.
    """
    etree = Writer(phyloxml)
    etree.write(file, encoding)


class Writer(object):
    """
    """
    def __init__(self, phyloxml):
        """Build an ElementTree from a phyloXML object."""
        pass

