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
    raise NotImplementedError("This function doesn't work yet.")
    nex = Nexus.Nexus()
    if isinstance(obj, list):
        nex.trees = obj 
    else:
        nex.trees = list(obj)
    Nexus.Nexus.write_nexus_data(nex, file, **kwargs)
    return len(nex.trees)
