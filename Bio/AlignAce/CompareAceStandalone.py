# Copyright 2003 by Bartek Wilczynski.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""

This module provides code to work with the standalone version of CompareAce, 
for motif comparison

CompareACE homepage:

http://atlas.med.harvard.edu/

functions:
CompareAce - runs the AlignACE standalone prgram and returns the ApplicationResult object
"""

import os
import re

from Bio import File
from Applications import CompareAceCommandline

import Scanner
import Parser


def CompareAce( cmd="CompareACE", **keywds):
    """Runs CompareACE and returns data.

    motif1, motif2 == files containing AlignACE motifs
    """

    if not os.path.exists(cmd):
        raise IOError("Executable does not exist at %s" % cmd)
    
    CompareCmd = CompareAceCommandline(cmd)

    for (par,val) in keywds.iteritems():
        CompareCmd.set_parameter(par,val)

    return CompareCmd.run()
