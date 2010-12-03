# Copyright 2009 by Cymon J. Cox and Brad Chapman. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Command line wrapper for the multiple alignment program TCOFFEE.
"""

from Bio.Application import _Option, _Switch, _Argument, AbstractCommandline

class TCoffeeCommandline(AbstractCommandline):
    """Commandline object for the TCoffee alignment program.

    http://www.tcoffee.org/Projects_home_page/t_coffee_home_page.html
    
    The T-Coffee command line tool has a lot of switches and options.
    This wrapper implements a VERY limited number of options - if you
    would like to help improve it please get in touch.

    Example:

    To align a FASTA file (unaligned.fasta) with the output in ClustalW
    format (file aligned.aln), and otherwise default settings, use:

    >>> from Bio.Align.Applications import TCoffeeCommandline
    >>> tcoffee_cline = TCoffeeCommandline(infile="unaligned.fasta",
    ...                                    output="clustalw",
    ...                                    outfile="aligned.aln")
    >>> print tcoffee_cline
    t_coffee -output clustalw -infile unaligned.fasta -outfile aligned.aln

    You would typically run the command line with tcoffee_cline() or via
    the Python subprocess module, as described in the Biopython tutorial.
    
    Citation:

    T-Coffee: A novel method for multiple sequence alignments.
    Notredame, Higgins, Heringa, JMB,302(205-217) 2000

    Last checked against: Version_6.92
    """
    SEQ_TYPES = ["dna","protein","dna_protein"]

    def __init__(self, cmd="t_coffee", **kwargs):
        self.parameters = \
          [_Option(["-output", "output"], [],
                    None,
                    0,
                    """Specify the output type.

                    One (or more separated by a comma) of:
                    'clustalw_aln', 'clustalw', 'gcg', 'msf_aln',
                    'pir_aln', 'fasta_aln', 'phylip', 'pir_seq', 'fasta_seq'

                    Note that of these Biopython's AlignIO module will only
                    read clustalw, pir, and fasta.
                    """,
                    0),
           _Option(["-infile", "infile"], ["file"],
                    None,
                    1,
                    "Specify the input file.",
                    0,),
           #Indicates the name of the alignment output by t_coffee. If the
           #default is used, the alignment is named <your sequences>.aln
           _Option(["-outfile", "outfile"], ["file"],
                    None,
                    0,
                    "Specify the output file. Default: <your sequences>.aln",
                    0),
           _Switch(["-convert", "convert"],
                    "Specify you want to perform a file conversion"),
           _Option(["-type", "type"], [],
                    lambda x: x in self.SEQ_TYPES,
                    0,
                    "Specify the type of sequence being aligned",
                    0),
           _Option(["-outorder", "outorder"], [],
                    None,
                    0,
                    "Specify the order of sequence to output"
                    "Either 'input', 'aligned' or <filename> of "
                    "Fasta file with sequence order",
                    0),
           _Option(["-matrix", "matrix"], [],
                    None,
                    0,
                    "Specify the filename of the substitution matrix to use."
                    "Default: blosum62mt",
                    0),
           _Option(["-gapopen", "gapopen"], [],
                    lambda x: isinstance(x, int),
                    0,
                    "Indicates the penalty applied for opening a gap "
                    "(negative integer)",
                    0),
           _Option(["-gapext", "gapext"], [],
                    lambda x: isinstance(x, int),
                    0,
                    "Indicates the penalty applied for extending a "
                    "gap. (negative integer)",
                    0),
           _Switch(["-quiet", "quiet"],
                    "Turn off log output"),
           _Option(["-mode", "mode"], [],
                    None,
                    0,
                    "Specifies a special mode: genome, quickaln, dali, 3dcoffee",
                    0)
           ]
        AbstractCommandline.__init__(self, cmd, **kwargs)           
