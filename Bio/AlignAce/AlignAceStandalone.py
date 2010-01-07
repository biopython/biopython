# Copyright 2003 by Bartek Wilczynski.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""

This module provides code to work with the standalone version of AlignACE, 
for motif search in DNA sequences.

AlignACE homepage:

http://atlas.med.harvard.edu/

AlignACE Citations:

Computational identification of cis-regulatory elements associated with 
groups of functionally related genes in Saccharomyces cerevisiae, 
Hughes, JD, Estep, PW, Tavazoie S, & GM Church, Journal of Molecular 
Biology 2000 Mar 10;296(5):1205-14.

Finding DNA Regulatory Motifs within Unaligned Non-Coding Sequences 
Clustered by Whole-Genome mRNA Quantitation, 
Roth, FR, Hughes, JD, Estep, PE & GM Church, Nature Biotechnology 
1998 Oct;16(10):939-45. 

functions:
AlignAce - runs the AlignACE standalone prgram and returns the 
ApplicationResult object
"""

import os

from Applications import AlignAceCommandline



def AlignAce(infile, cmd="AlignACE", **keywds):
    """Runs AlignACE and returns data.

    cmd == AlignACE executable
    infile == sequence file to process
    
    You may pass more parameters to **keywds to change the behavior of
    the search.  Otherwise, optional values will be chosen by blastall.

    numcols      number of columns to align (10)
    expect       number of sites expected in model (10)
    gcback       background fractional GC content of input sequence (0.38)
    minpass      minimum number of non-improved passes in phase 1 (200)
    seed         set seed for random number generator (time)
    undersample  possible sites / (expect * numcols * seedings) (1)
    oversample   1/undersample (1)
    """

    if not os.path.exists(cmd):
        raise IOError("Executable does not exist at %s" % cmd)

    if not os.path.exists(infile):
        raise IOError("Input file does not exist at %s" % infile)
    
    AlignCmd = AlignAceCommandline(cmd)

    AlignCmd.set_parameter("input",infile)
    
    for (par,val) in keywds.iteritems():
        AlignCmd.set_parameter(par,val)

    return AlignCmd.run()
