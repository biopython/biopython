# Copyright 2009 by Peter Cock & Cymon J. Cox.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Alignment command line tool wrappers."""

__docformat__ = "restructuredtext en"  # Don't just use plain text in epydoc API pages!

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
__all__ = ["MuscleCommandline",
           "ClustalwCommandline",
           "ClustalOmegaCommandline",
           "PrankCommandline",
           "MafftCommandline",
           "DialignCommandline",
           "ProbconsCommandline",
           "TCoffeeCommandline",
           "MSAProbsCommandline",
           ]
