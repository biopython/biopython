# Copyright 2009 by Cymon J. Cox and Brad Chapman. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Command line wrapper for the multiple alignment program TCOFFEE.
"""

__docformat__ = "epytext en" #Don't just use plain text in epydoc API pages!

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
        self.parameters = [
           _Option(["-output", "output"],
                   """Specify the output type.
                   One (or more separated by a comma) of:
                   'clustalw_aln', 'clustalw', 'gcg', 'msf_aln',
                   'pir_aln', 'fasta_aln', 'phylip', 'pir_seq', 'fasta_seq'

                   Note that of these Biopython's AlignIO module will only
                   read clustalw, pir, and fasta.
                   """, #TODO - Can we read the PHYLIP output?
                   equate=False),
           _Option(["-infile", "infile"],
                   "Specify the input file.",
                   filename=True,
                   is_required=True,
                   equate=False),
           #Indicates the name of the alignment output by t_coffee. If the
           #default is used, the alignment is named <your sequences>.aln
           _Option(["-outfile", "outfile"],
                   "Specify the output file. Default: <your sequences>.aln",
                   filename=True,
                   equate=False),
           _Switch(["-convert", "convert"],
                   "Specify you want to perform a file conversion"),
           _Option(["-type", "type"],
                   "Specify the type of sequence being aligned",
                   checker_function=lambda x: x in self.SEQ_TYPES,
                   equate=False),
           _Option(["-outorder", "outorder"],
                   "Specify the order of sequence to output"
                   "Either 'input', 'aligned' or <filename> of "
                   "Fasta file with sequence order",
                   equate=False),
           _Option(["-matrix", "matrix"],
                   "Specify the filename of the substitution matrix to use."
                   "Default: blosum62mt",
                   equate=False),
           _Option(["-gapopen", "gapopen"],
                   "Indicates the penalty applied for opening a gap "
                   "(negative integer)",
                   checker_function=lambda x: isinstance(x, int),
                   equate=False),
           _Option(["-gapext", "gapext"],
                   "Indicates the penalty applied for extending a "
                   "gap. (negative integer)",
                   checker_function=lambda x: isinstance(x, int),
                   equate=False),
           _Switch(["-quiet", "quiet"],
                   "Turn off log output"),
           _Option(["-mode", "mode"],
                   "Specifies a special mode: genome, quickaln, dali, 3dcoffee",
                   equate=False),
           ]
        AbstractCommandline.__init__(self, cmd, **kwargs)           

def _test():
    """Run the module's doctests (PRIVATE)."""
    print "Runing modules doctests..."
    import doctest
    doctest.testmod()
    print "Done"

if __name__ == "__main__":
    _test()
