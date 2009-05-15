"""Definitions for interacting with AlignAce.
"""
from Bio.Application import _Option,_Argument, AbstractCommandline
from Bio.Application import generic_run

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

    def run(self):
        return generic_run(self)


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

    def run(self):
        return generic_run(self)
