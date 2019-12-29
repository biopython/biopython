# Copyright 2009 by Peter Cock & Cymon J. Cox.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Alignment command line tool wrappers."""

from ._Muscle import MuscleCommandline
from ._Clustalw import ClustalwCommandline
from ._ClustalOmega import ClustalOmegaCommandline
from ._Prank import PrankCommandline
from ._Mafft import MafftCommandline
from ._Dialign import DialignCommandline
from ._Probcons import ProbconsCommandline
from ._TCoffee import TCoffeeCommandline
from ._MSAProbs import MSAProbsCommandline

# Make this explicit, then they show up in the API docs
__all__ = (
    "MuscleCommandline",
    "ClustalwCommandline",
    "ClustalOmegaCommandline",
    "PrankCommandline",
    "MafftCommandline",
    "DialignCommandline",
    "ProbconsCommandline",
    "TCoffeeCommandline",
    "MSAProbsCommandline",
)
