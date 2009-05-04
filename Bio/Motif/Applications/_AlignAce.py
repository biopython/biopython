# Copyright 2003-2009 by Bartek Wilczynski.  All rights reserved.
# Revisions copyright 2009 by Peter Cock.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""This module provides code to work with the standalone version of AlignACE, 
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

"""
from Bio.Application import AbstractCommandline, _Option, _Argument


class AlignAceCommandline(AbstractCommandline):
    """Create a commandline for the AlignAce program.

    XXX This could use more checking for valid paramters to the program.
    """
    def __init__(self, cmd="AlignACE", **kwargs):
        self.parameters = \
          [
            _Option(["-i","input"],["input"],lambda x : x.__class__== str,1,
                    "Input Sequence file in FASTA format."),
            
            _Option(["-numcols","numcols"],["input"],lambda x : x.__class__== int,0,
                    "Number of columns to align"),

            _Option(["-expect","expect"],["input"],lambda x : x.__class__== int,0,
                    "number of sites expected in model "),
            
            _Option(["-gcback","gcback"],["input"],lambda x : x.__class__== float,0,
                    "background fractional GC content of input sequence"),
            
            _Option(["-minpass","minpass"],["input"],lambda x : x.__class__== int,0,
                    "minimum number of non-improved passes in phase 1"),
            
            _Option(["-seed","seed"],["input"],lambda x : x.__class__== int,0,
                    "set seed for random number generator (time)"),
            
            _Option(["-undersample","undersample"],["input"],lambda x : x.__class__== int,0,
                    "possible sites / (expect * numcols * seedings)"),

            _Option(["-oversample","oversample"],["input"],lambda x : x.__class__== int,0,
                    "1/undersample"),
          ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


class CompareAceCommandline(AbstractCommandline):
    """Create a commandline for the CompareAce program.

    XXX This could use more checking for valid paramters to the program.
    """
    def __init__(self, cmd="CompareACE", **kwargs):
        import os.path
        self.parameters = \
          [
            _Argument(["motif1"],["input","file"], os.path.exists,1,"name of file containing motif 1"),
            _Argument(["motif2"],["input","file"], os.path.exists,1,"name of file containing motif 2"),
          ]
        AbstractCommandline.__init__(self, cmd, **kwargs)
