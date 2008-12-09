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
CompareAce - runs the AlignACE standalone prgram and returns the ApplicationResult object
"""

import os

from Bio import File
from Bio.ParserSupport import *
from Bio import Application
from Bio.Application import _Option,_Argument
import AlignAceParser


class AlignAceCommandline(Application.AbstractCommandline):
    """Create a commandline for the AlignAce program.

    XXX This could use more checking for valid paramters to the program.
    """
    def __init__(self, cmd = "AlignACE"):

        Application.AbstractCommandline.__init__(self)
        self.program_name = cmd

        self.parameters = \
          [
            _Option(["-i","input","Sequence File"],["input"],lambda x : x.__class__== str,1,
                    "Input Sequence file in FASTA format."),
            
            _Option(["-numcols","numcols","number of columns to align"],["input"],lambda x : x.__class__== int,0,
                    "Number of columns to align"),

            _Option(["-expect","expect","number of sites expected in model "],["input"],lambda x : x.__class__== int,0,
                    "number of sites expected in model "),
            
            _Option(["-gcback","gcback","background fractional GC content of input sequence"],["input"],lambda x : x.__class__== float,0,
                    "background fractional GC content of input sequence"),
            
            _Option(["-minpass","minpass","minimum number of non-improved passes in phase 1"],["input"],lambda x : x.__class__== int,0,
                    "minimum number of non-improved passes in phase 1"),
            
            _Option(["-seed","seed","set seed for random number generator (time)"],["input"],lambda x : x.__class__== int,0,
                    "set seed for random number generator (time)"),
            
            _Option(["-undersample","undersample","possible sites / (expect * numcols * seedings)"],["input"],lambda x : x.__class__== int,0,
                    "possible sites / (expect * numcols * seedings)"),

            _Option(["-oversample","oversample","1/undersample"],["input"],lambda x : x.__class__== int,0,
                    "1/undersample"),
          ]

    def run(self):
        return Application.generic_run(self)



class CompareAceCommandline(Application.AbstractCommandline):
    """Create a commandline for the CompareAce program.

    XXX This could use more checking for valid paramters to the program.
    """
    def __init__(self, cmd = "CompareACE"):

        import os.path
        Application.AbstractCommandline.__init__(self)
        self.program_name = cmd

        self.parameters = \
          [
            _Argument(["motif1"],["input","file"], os.path.exists,1,"name of file containing motif 1"),
            _Argument(["motif2"],["input","file"], os.path.exists,1,"name of file containing motif 2"),
          ]

    def run(self):
        return Application.generic_run(self)


def AlignAce(infile, cmd="AlignACE", **keywds):
    """Runs AlignACE and returns data.

    cmd == AlignACE executable
    infile == sequence file to process
    
    You may pass more parameters to **keywds to change the behavior of
    the search.  Otherwise, optional values will be chosen by blastall.

    numcols    	number of columns to align (10)
    expect     	number of sites expected in model (10)
    gcback     	background fractional GC content of input sequence (0.38)
    minpass    	minimum number of non-improved passes in phase 1 (200)
    seed       	set seed for random number generator (time)
    undersample	possible sites / (expect * numcols * seedings) (1)
    oversample	        1/undersample (1)
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

