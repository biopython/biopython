# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Modules for parsing and writing files in the PhyloXML format.

See www.phyloxml.org for more info about the format.
"""

import zipfile

from Parser import parse


def read(fname):
    if zipfile.is_zipfile(fname):
        z = zipfile.ZipFile(fname)
        tree = parse(z.open(z.filelist[0].filename))
        z.close()
    else:
        tree = parse(fname)
    return tree

