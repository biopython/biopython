"""Definitions for interacting with AlignAce.
"""
from Bio import Application
from Bio.Application import _Option,_Argument

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
