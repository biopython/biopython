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

class EProtDistCommandline(Application.AbstractCommandline):
    """Commandline object for the eprotdist program from EMBOSS.

    This is an EMBOSS wrapper around protdist from PHYLIP.
    """
    def __init__(self, cmd = "eprotdist"):
        Application.AbstractCommandline.__init__(self)
        self.program_name = cmd

        self.parameters = \
         [_Option(["-msf"], ["input"], None, 1,
                  "File containing sequences"),
          _Option(["-outfile"], ["output"], None, 1,
                  "Output file name"),
          _Option(["-method"], ["input"], None, 1,
                  "Choose the method to use"),
          _Option(["-categ"], ["input"], None, 0,
                  "Choose the categorie to use"),
          _Option(["-gencode"], ["input"], None, 0,
                  "Which genetic code"),
          _Option(["-prob"], ["input"], None, 0,
                  "Prob change category (1.0=easy)"),
          _Option(["-tranrate"], ["input"], None, 0,
                  "Transition/transversion ratio"),
          _Option(["-freqa"], ["input"], None, 0,
                  "Frequency for A"),
          _Option(["-freqc"], ["input"], None, 0,
                  "Frequency for C"),
          _Option(["-freqg"], ["input"], None, 0,
                  "Frequency for G"),
          _Option(["-freqt"], ["input"], None, 0,
                  "Frequency for T"),
          _Option(["-printdata"], ["input"], None, 0,
                  "Print out the data at start of run"),
          _Option(["-progress"], ["input"], None, 0,
                  "Print indications of progress of run"),
          _Option(["-basefrequency"], ["input"], None, 0,
                  "Use empirical base frequencies")]

class ENeighborCommandline(Application.AbstractCommandline):
    """Commandline object for the eneighbor program from EMBOSS.

    This is an EMBOSS wrapper around neighbor from PHYLIP.
    """
    def __init__(self, cmd = "eneighbor"):
        Application.AbstractCommandline.__init__(self)
        self.program_name = cmd

        self.parameters = \
         [_Option(["-infile"], ["input"], None, 1,
                  "infile value"),
          _Option(["-outfile"], ["output"], None, 1,
                  "Output file name"),
          _Option(["-trout"], ["input"], None, 1,
                  "Create a tree file"),
          _Option(["-treefile"], ["input"], None, 1,
                  "Tree file name"),
          _Option(["-nj"], ["input"], None, 1,
                  "Neighbor-joining"),
          _Option(["-noog"], ["input"], None, 1,
                  "Outgroup root"),
          _Option(["-outgnum"], ["input"], None, 0,
                  "number of the outgroup"),
          _Option(["-randseed"], ["input"], None, 0,
                  "Random number seed (must be odd)"),
          _Option(["-datasets"], ["input"], None, 0,
                  "How many data sets"),
          _Option(["-drawtree"], ["input"], None, 0,
                  "Draw tree"),
          _Option(["-lt"], ["input"], None, 0,
                  "Lower-triangular data matrix"),
          _Option(["-ut"], ["input"], None, 0,
                  "Upper-triangular data matrix"),
          _Option(["-sr"], ["input"], None, 0,
                  "Subreplicates"),
          _Option(["-random"], ["input"], None, 0,
                  "Randomize input order of species"),
          _Option(["-multsets"], ["input"], None, 0,
                  "Analyze multiple data sets"),
          _Option(["-printdata"], ["input"], None, 0,
                  "Print out the data at start of run"),
          _Option(["-progress"], ["input"], None, 0,
                  "Print indications of progress of run")]

class EProtParsCommandline(Application.AbstractCommandline):
    """Commandline object for the eprotpars program from EMBOSS.

    This is an EMBOSS wrapper around protpars from PHYLIP.
    """
    def __init__(self, cmd = "eprotpars"):
        Application.AbstractCommandline.__init__(self)
        self.program_name = cmd

        self.parameters = \
         [_Option(["-msf"], ["input", "file"], None, 1,
                  "Sequences file to be read in"),
          _Option(["-outfile"], ["output", "file"], None, 1,
                  "Output file"),
          _Option(["-besttree"], ["input"], None, 0,
                  "Search for the best tree"),
          _Option(["-random"], ["input"], None, 0,
                  "Randomize input order of species"),
          _Option(["-norandom"], ["input"], None, 0,
                  "Do not randomize input order of species"),
          _Option(["-randseed"], ["input"], None, 0,
                  "Random number seed (must be odd)"),
          _Option(["-randtimes"], ["input"], None, 0,
                  "How many times to randomize"),
          _Option(["-og"], ["input"], None, 0,
                  "Use an outgroup root"),
          _Option(["-noog"], ["input"], None, 0,
                  "Do not use an outgroup root"),
          _Option(["-outgnum"], ["input"], None, 0,
                  "Number of the outgroup"),
          _Option(["-thresh"], ["input"], None, 0,
                  "Use Threshold parsimony"),
          _Option(["-valthresh"], ["input"], None, 0,
                  "threshold value"),
          _Option(["-printdata"], ["input"], None, 0,
                  "Print out the data at start of run"),
          _Option(["-progress"], ["input"], None, 0,
                  "Print indications of progress of run"),
          _Option(["-steps"], ["input"], None, 0,
                  "Print out steps in each site"),
          _Option(["-seqatnodes"], ["input"], None, 0,
                  "Print sequences at all nodes of tree"),
          _Option(["-drawtree"], ["input"], None, 0,
                  "Draw tree"),
          _Option(["-trout"], ["input"], None, 0,
                  "Create a tree file"),
          _Option(["-notrout"], ["input"], None, 0,
                  "Do not create a tree file"),
          _Option(["-treefile"], ["output", "file"], None, 0,
                  "Output treefile name")]

class EConsenseCommandline(Application.AbstractCommandline):
    """Commandline object for the econsense program from EMBOSS.

    This is an EMBOSS wrapper around consense from PHYLIP.
    """
    def __init__(self, cmd = "econsense"):
        Application.AbstractCommandline.__init__(self)
        self.program_name = cmd

        self.parameters = \
         [_Option(["-infile"], ["input", "file"], None, 1,
                  "file to read in (New Hampshire standard form)"),
          _Option(["-outfile"], ["output", "file"], None, 1,
                  "Output file name"),
          _Option(["-notrout"], ["input"], None, 0,
                  "Do not create a tree file"),
          _Option(["-trout"], ["input"], None, 0,
                  "Create a tree file"),
          _Option(["-treefile"], ["output", "file"], None, 0,
                  "tree file name"),
          _Option(["-noog"], ["input"], None, 0,
                  "Do not use an outgroup"),
          _Option(["-og"], ["input"], None, 0,
                  "Use an outgroup"),
          _Option(["-outgnum"], ["input"], None, 0,
                  "number of the outgroup"),
          _Option(["-nodrawtree"], ["input"], None, 0,
                  "Do not draw a tree"),
          _Option(["-drawtree"], ["input"], None, 0,
                  "Draw tree"),
          _Option(["-root"], ["input"], None, 0,
                  "Trees to be treated as Rooted"),
          _Option(["-progress"], ["input"], None, 0,
                  "Print indications of the progress of run"),
          _Option(["-noprintsets"], ["input"], None, 0,
                  "Do not print out the sets of species"),
          _Option(["-printsets"], ["input"], None, 0,
                  "Print out the sets of species")]
                  
