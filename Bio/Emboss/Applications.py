"""Code to interact with and run various EMBOSS programs.

These classes follow the AbstractCommandline interfaces for running
programs.
"""

from Bio import Application
from Bio.Application import _Option

class Primer3Commandline(Application.AbstractCommandline):
    """Commandline object for the Primer3 interface from EMBOSS.
    """
    def __init__(self, cmd = "primer3"):
        Application.AbstractCommandline.__init__(self)
        self.program_name = cmd

        self.parameters = \
          [_Option(["-sequence"], ["input"], None, 1, 
                   "Sequence to choose primers from"),
           _Option(["-outfile"], ["output", "file"], None, 1,
                   "Output file name"),
           _Option(["-task"], ["input"], None, 0),
           _Option(["-numreturn"], ["input"], None, 0),
           _Option(["-includedregion"], ["input"], None, 0),
           _Option(["-target"], ["input"], None, 0),
           _Option(["-excludedregion"], ["input"], None, 0),
           _Option(["-forwardinput"], ["input"], None, 0),
           _Option(["-reverseinput"], ["input"], None, 0),
           _Option(["-gcclamp"], ["input"], None, 0),
           _Option(["-osize"], ["input"], None, 0),
           _Option(["-minsize"], ["input"], None, 0),
           _Option(["-maxsize"], ["input"], None, 0),
           _Option(["-otm"], ["input"], None, 0),
           _Option(["-mintm"], ["input"], None, 0),
           _Option(["-maxtm"], ["input"], None, 0),
           _Option(["-maxdifftm"], ["input"], None, 0),
           _Option(["-ogcpercent"], ["input"], None, 0),
           _Option(["-mingc"], ["input"], None, 0),
           _Option(["-maxgc"], ["input"], None, 0),
           _Option(["-saltconc"], ["input"], None, 0),
           _Option(["-dnaconc"], ["input"], None, 0),
           _Option(["-maxployx"], ["input"], None, 0),
           _Option(["-productosize"], ["input"], None, 0),
           _Option(["-productsizerange"], ["input"], None, 0),
           _Option(["-productotm"], ["input"], None, 0),
           _Option(["-productmintm"], ["input"], None, 0),
           _Option(["-productmaxtm"], ["input"], None, 0),
           _Option(["-oligoexcluderegion"], ["input"], None, 0),
           _Option(["-oligoinput"], ["input"], None, 0),
           _Option(["-oligosize"], ["input"], None, 0),
           _Option(["-oligominsize"], ["input"], None, 0),
           _Option(["-oligomaxsize"], ["input"], None, 0),
           _Option(["-oligotm"], ["input"], None, 0),
           _Option(["-oligomintm"], ["input"], None, 0),
           _Option(["-oligomaxtm"], ["input"], None, 0),
           _Option(["-oligoogcpercent"], ["input"], None, 0),
           _Option(["-oligomingc"], ["input"], None, 0),
           _Option(["-oligomaxgc"], ["input"], None, 0),
           _Option(["-oligosaltconc"], ["input"], None, 0),
           _Option(["-oligodnaconc"], ["input"], None, 0),
           _Option(["-oligoselfany"], ["input"], None, 0),
           _Option(["-oligoselfend"], ["input"], None, 0),
           _Option(["-oligomaxpolyx"], ["input"], None, 0)]

class PrimerSearchCommandline(Application.AbstractCommandline):
    """Commandline object for the primersearch program from EMBOSS.
    """
    def __init__(self, cmd = "primersearch"):
        Application.AbstractCommandline.__init__(self)
        self.program_name = cmd

        self.parameters = \
         [_Option(["-sequences"], ["input"], None, 1,
                  "Sequence to look for the primer pairs in."),
          _Option(["-primers"], ["input", "file"], None, 1,
                  "File containing the primer pairs to search for."),
          _Option(["-out"], ["output", "file"], None, 1,
                  "Name of the output file."),
          _Option(["-mismatchpercent"], ["input"], None, 1,
                  "Allowed percentage mismatch.")]


