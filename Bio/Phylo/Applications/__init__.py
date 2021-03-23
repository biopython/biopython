# Copyright 2011 by Eric Talevich.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Phylogenetics command line tool wrappers (OBSOLETE).

We have decided to remove this module in future, and instead recommend
building your command and invoking it via the subprocess module directly.
"""

from ._Phyml import PhymlCommandline
from ._Raxml import RaxmlCommandline
from ._Fasttree import FastTreeCommandline

# Make this explicit, then they show up in the API docs
__all__ = ("PhymlCommandline", "RaxmlCommandline", "FastTreeCommandline")
