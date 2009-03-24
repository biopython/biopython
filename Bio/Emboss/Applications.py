"""Code to interact with and run various EMBOSS programs.

These classes follow the AbstractCommandline interfaces for running
programs.
"""

from Bio import Application
from Bio.Application import _Option

class Primer3Commandline(Application.AbstractCommandline):
    """Commandline object for the Primer3 interface from EMBOSS.
    """
    def __init__(self, cmd = "eprimer3"):
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
           _Option(["-oligomaxpolyx"], ["input"], None, 0),
           _Option(["-mispriminglibraryfile"], ["input"], None, 0),
           _Option(["-maxmispriming"], ["input"], None, 0),
           _Option(["-oligomishyblibraryfile"], ["input"], None, 0),
           _Option(["-oligomaxmishyb"], ["input"], None, 0),
           _Option(["-explainflag"], ["input"], None, 0),
           ]

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

class ESeqBootCommandline(Application.AbstractCommandline):
    """Commandline object for the eseqboot program from EMBOSS.

    This is an EMBOSS wrapper around seqboot from PHYLIP.
    """
    def __init__(self, cmd = "eseqboot"):
        Application.AbstractCommandline.__init__(self)
        self.program_name = cmd

        self.parameters = \
         [_Option(["-datafile"], ["input", "file"], None, 1,
                  "Input file"),
          _Option(["-outfile"], ["output", "file"], None, 1,
                  "Output file name"),
          _Option(["-randseed"], ["input"], None, 1,
                  "Random number seed (must be odd)"),
          _Option(["-method"], ["input"], None, 1,
                  "Choose the method"),
          _Option(["-test"], ["input"], None, 1,
                  "Choose test"),
          _Option(["-reps"], ["input"], None, 1,
                  "How many replicates"),
          _Option(["-inter"], ["input"], None, 0,
                  "Interleaved input"),
          _Option(["-enzymes"], ["input"], None, 0,
                  "Present in input file"),
          _Option(["-all"], ["input"], None, 0,
                  "All alleles present at each locus"),
          _Option(["-printdata"], ["input"], None, 0,
                  "Print out the data at start of run"),
          _Option(["-progress"], ["input"], None, 0,
                  "Print indications of progress of run")]

class WaterCommandline(Application.AbstractCommandline):
    """Commandline object for the water program from EMBOSS.
    """
    def __init__(self, cmd = "water"):
        Application.AbstractCommandline.__init__(self)
        self.program_name = cmd

        self.parameters = \
         [_Option(["-asequence"], ["input", "file"], None, 1,
                  "First sequence to align"),
         _Option(["-bsequence"], ["input", "file"], None, 1,
                  "Second sequence to align"),
         _Option(["-gapopen"], ["input"], None, 1,
                 "Gap open penalty"),
         _Option(["-gapextend"], ["input"], None, 1,
                 "Gap extension penalty"),
         _Option(["-outfile"], ["output", "file"], None, 1,
                 "Output file for the alignment"),
         _Option(["-datafile"], ["input", "file"], None, 0,
                 "Matrix file"),
         _Option(["-similarity"], ["input"], None, 0,
                 "Display percent identity and similarity"),
         _Option(["-nosimilarity"], ["input"], None, 0,
                 "Do not display percent identity and similarity"),
         _Option(["-aformat"], ["input"], None, 0,
                 "Display output in a different specified output format")]

class FuzznucCommandline(Application.AbstractCommandline):
    """Commandline object for the fuzznuc program from EMBOSS.
    """
    def __init__(self, cmd = "fuzznuc"):
        Application.AbstractCommandline.__init__(self)
        self.program_name = cmd

        self.parameters = [
         _Option(["-sequence"], ["input"], None, 1,
                 "Sequence database USA"),
         _Option(["-pattern"], ["input"], None, 1,
                 "Search pattern, using standard IUPAC one-letter codes"),
         _Option(["-mismatch"], ["input"], None, 1,
                 "Number of mismatches"),
         _Option(["-outfile"], ["output", "file"], None, 1,
                 "Output report file name"),
         _Option(["-complement"], ["input"], None, 0,
                 "Search complementary strand"),
         _Option(["-rformat"], ["input"], None, 0,
                 "Specify the report format to output in.")]

class Est2GenomeCommandline(Application.AbstractCommandline):
    """Commandline object for the est2genome program from EMBOSS.
    """
    def __init__(self, cmd = "est2genome"):
        Application.AbstractCommandline.__init__(self)
        self.program_name = cmd

        self.parameters = [
         _Option(["-est"], ["input"], None, 1,
                 "EST sequence(s)"),
         _Option(["-genome"], ["input"], None, 1,
                 "Genomic sequence"),
         _Option(["-outfile"], ["output", "file"], None, 1,
                 "Output file name"),
         _Option(["-match"], ["input"], None, 0, 
                 "Score for matching two bases"),
         _Option(["-mismatch"], ["input"], None, 0,
                 "Cost for mismatching two bases"),
         _Option(["-gappenalty"], ["input"], None, 0,
                 "Cost for deleting a single base in either sequence, " + \
                 "excluding introns"),
         _Option(["-intronpenalty"], ["input"], None, 0,
                 "Cost for an intron, independent of length."),
         _Option(["-splicepenalty"], ["input"], None, 0,
                 "Cost for an intron, independent of length " + \
                 "and starting/ending on donor-acceptor sites"),
         _Option(["-minscore"], ["input"], None, 0,
                 "Exclude alignments with scores below this threshold score."),
         _Option(["-reverse"], ["input"], None, 0,
                 "Reverse the orientation of the EST sequence"),
         _Option(["-splice"], ["input"], None, 0,
                 "Use donor and acceptor splice sites."),
         _Option(["-mode"], ["input"], None, 0,
                 "This determines the comparion mode. 'both', 'forward' " + \
                 "'reverse'"),
         _Option(["-best"], ["input"], None, 0,
                 "You can print out all comparisons instead of just the best"),
         _Option(["-space"], ["input"], None, 0,
                 "for linear-space recursion."),
         _Option(["-shuffle"], ["input"], None, 0,
                 "Shuffle"),
         _Option(["-seed"], ["input"], None, 0,
                 "Random number seed"),
         _Option(["-align"], ["input"], None, 0,
                 "Show the alignment."),
         _Option(["-width"], ["input"], None, 0,
                 "Alignment width")
        ]

class ETandemCommandline(Application.AbstractCommandline):
    """Commandline object for the etandem program from EMBOSS.
    """
    def __init__(self, cmd = "etandem"):
        Application.AbstractCommandline.__init__(self)
        self.program_name = cmd

        self.parameters = [
         _Option(["-sequence"], ["input", "file"], None, 1,
                 "Sequence"),
         _Option(["-minrepeat"], ["input"], None, 1,
                 "Minimum repeat size"),
         _Option(["-maxrepeat"], ["input"], None, 1,
                 "Maximum repeat size"),
         _Option(["-outfile"], ["output", "file"] , None, 1,
                 "Output report file name"),
         _Option(["-threshold"], ["input"], None, 0,
                 "Threshold score"),
         _Option(["-mismatch"], ["input"], None, 0,
                   "Allow N as a mismatch"),
         _Option(["-uniform"], ["input"], None, 0,
                   "Allow uniform consensus"),
         _Option(["-rformat"], ["output"], None, 0,
                 "Output report format")]

class EInvertedCommandline(Application.AbstractCommandline):
    """Commandline object for the einverted program from EMBOSS.
    """
    def __init__(self, cmd = "einverted"):
        Application.AbstractCommandline.__init__(self)
        self.program_name = cmd

        self.parameters = [
         _Option(["-sequence"], ["input", "file"], None, 1,
                 "Sequence"),
         _Option(["-gap"], ["input", "file"], None, 1,
                 "Gap penalty"),
         _Option(["-threshold"], ["input"], None, 1,
                 "Minimum score threshold"),
         _Option(["-match"], ["input"], None, 1,
                 "Match score"),
         _Option(["-mismatch"], ["input"], None, 1,
                   "Mismatch score"),
         _Option(["-outfile"], ["output", "file"] , None, 1,
                 "Output report file name"),
         _Option(["-maxrepeat"], ["input"], None, 0,
                 "Maximum separation between the start and end of repeat"),
         ]

class PalindromeCommandline(Application.AbstractCommandline):
    """Commandline object for the palindrome program from EMBOSS.
    """
    def __init__(self, cmd = "palindrome"):
        Application.AbstractCommandline.__init__(self)
        self.program_name = cmd

        self.parameters = [
         _Option(["-sequence"], ["input", "file"], None, 1,
                 "Sequence"),
         _Option(["-minpallen"], ["input"], None, 1,
                 "Minimum palindrome length"),
         _Option(["-maxpallen"], ["input"], None, 1,
                 "Maximum palindrome length"),
         _Option(["-gaplimit"], ["input"], None, 1,
                 "Maximum gap between repeats"),
         _Option(["-nummismatches"], ["input"], None, 1,
                 "Number of mismatches allowed"),
         _Option(["-overlap"], ["input"], None, 1,
                 "Report overlapping matches"),
         _Option(["-outfile"], ["output", "file"] , None, 1,
                 "Output report file name"),
         ]

class TranalignCommandline(Application.AbstractCommandline):
    """Commandline object for the tranalign program from EMBOSS.
    """
    def __init__(self, cmd = "tranalign"):
        Application.AbstractCommandline.__init__(self)
        self.program_name = cmd

        self.parameters = [
         _Option(["-asequence"], ["input", "file"], None, 1,
                 "Nucleotide sequences to be aligned."),
         _Option(["-bsequence"], ["input", "file"], None, 1,
                 "Protein sequence alignment"),
         _Option(["-outseq"], ["output", "file"], None, 1,
                 "Output sequence file."),
         _Option(["-table"], ["input"], None, 0,
                 "Code to use")]

class DiffseqCommandline(Application.AbstractCommandline):
    """Commandline object for the diffseq program from EMBOSS.
    """
    def __init__(self, cmd = "diffseq"):
        Application.AbstractCommandline.__init__(self)
        self.program_name = cmd

        self.parameters = [
         _Option(["-asequence"], ["input", "file"], None, 1,
                 "First sequence to compare"),
         _Option(["-bsequence"], ["input", "file"], None, 1,
                 "Second sequence to compare"),
         _Option(["-wordsize"], ["input"], None, 1,
                 "Word size to use for comparisons (10 default)"),
         _Option(["-outfile"], ["output", "file"], None, 1,
                "Output report file name"),
         _Option(["-aoutfeat"], ["output", "file"], None, 1,
                "File for output of first sequence's features"),
         _Option(["-boutfeat"], ["output", "file"], None, 1,
                "File for output of second sequence's features"),
         _Option(["-rformat"], ["output"], None, 0,
                 "Output report file format")
         ]

class IepCommandline(Application.AbstractCommandline):
    """Commandline for EMBOSS iep: calculated isoelectric point and charge.
    """
    def __init__(self, cmd = "iep"):
        Application.AbstractCommandline.__init__(self)
        self.program_name = cmd

        self.parameters = [
         _Option(["-sequence"], ["input", "file"], None, 1,
                "Protein sequence(s) filename"),
         _Option(["-outfile"], ["output", "file"], None, 1,
                "Output report file name"),
         _Option(["-amino"], ["input"], None, 0),
         _Option(["-lysinemodified"], ["input"], None, 0),
         _Option(["-disulphides"], ["input"], None, 0),
         _Option(["-notermini"], ["input"], None, 0),
         ]
