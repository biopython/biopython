# Copyright 2009 by Cymon J. Cox and Brad Chapman. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Command line wrapper for the multiple alignment program TCOFFEE.

http://www.tcoffee.org/Projects_home_page/t_coffee_home_page.html

T-Coffee: A novel method for multiple sequence alignments.
Notredame, Higgins, Heringa, JMB,302(205-217) 2000

Last checked against: Version_6.92
"""

import types
from Bio.Application import _Option, _Switch, _Argument, AbstractCommandline

class TCoffeeCommandline(AbstractCommandline):
    """Commandline object for the TCoffee alignment program.
    
    Implements a VERY limited number of options.
    """
    SEQ_TYPES = ["dna","protein","dna_protein"]

    def __init__(self, cmd="t_coffee", **kwargs):
        self.parameters = \
          [_Option(["-output", "output"], ["output"],
                    None,
                    0,
                    "Specify the output type. "
                    "One (or more separated by a comma) of: "
                    "'clustalw_aln', 'clustalw', 'gcg', 'msf_aln', "
                    "'pir_aln', 'fasta_aln', 'phylip', 'pir_seq', 'fasta_seq'"
                    "Note that biopython will only read clustalw, pir, and fasta",
                    0),
           _Option(["-infile", "infile"], ["input"],
                    None,
                    1,
                    "Specify the input file.",
                    0,),
           #Indicates the name of the alignment output by t_coffee. If the
           #default is used, the alignment is named <your sequences>.aln
           _Option(["-outfile", "outfile"], ["output"],
                    None,
                    0,
                    "Specify the output file. Default: <your sequences>.aln",
                    0),
           _Switch(["-convert", "convert"], ["input"],
                    "Specify you want to perform a file conversion"),
           _Option(["-type"], ["input"],
                    lambda x: x in SEQ_TYPES,
                    0,
                    "Specify the type of sequence being aligned",
                    0),
           _Option(["-outorder", "outorder"], ["input"],
                    None,
                    0,
                    "Specify the order of sequence to output"
                    "Either 'input', 'aligned' or <filename> of "
                    "Fasta file with sequence order",
                    0),
           _Option(["-matrix", "matrix"], ["input"],
                    None,
                    0,
                    "Specify the filename of the substitution matrix to use."
                    "Default: blosum62mt",
                    0),
           _Option(["-gapopen", "gapopen"], ["input"],
                    lambda x: isinstance(x, types.IntType),
                    0,
                    "Indicates the penalty applied for opening a gap "
                    "(negative integer)",
                    0),
           _Option(["-gapext", "gapext"], ["input"],
                    lambda x: isinstance(x, types.IntType),
                    0,
                    "Indicates the penalty applied for extending a "
                    "gap. (negative integer)",
                    0),
           _Switch(["-quiet", "quiet"], ["input"],
                    "Turn off log output"),
           _Option(["-mode", "mode"], ["input"],
                    None,
                    0,
                    "Specifies a special mode: genome, quickaln, dali, 3dcoffee",
                    0)
           ]
        AbstractCommandline.__init__(self, cmd, **kwargs)           
