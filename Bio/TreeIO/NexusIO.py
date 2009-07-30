# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""I/O function wrappers for Bio.Nexus trees.
"""
__docformat__ = "epytext en"

from Bio.Nexus import Nexus


def parse(file):
    nex = Nexus.Nexus(file)
    return iter(nex.trees)


def write(obj, file, **kwargs):
    do_close = False
    if not hasattr(file, 'write'):
        file = open(file, 'w+')
        do_close = True
    try:
        ohandle = Nexus.Nexus.write_nexus_data(obj, file, **kwargs)
    finally:
        if do_close:
            file.close()
