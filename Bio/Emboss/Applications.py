# Copyright 2001-2009 Brad Chapman.
# Revisions copyright 2009-2016 by Peter Cock.
# Revisions copyright 2009 by David Winter.
# Revisions copyright 2009-2010 by Leighton Pritchard.
# All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Code to interact with and run various EMBOSS programs (OBSOLETE).

These classes follow the AbstractCommandline interfaces for running
programs.

We have decided to remove this module in future, and instead recommend
building your command and invoking it via the subprocess module directly.
"""


from Bio.Application import _Option, _Switch, AbstractCommandline


class _EmbossMinimalCommandLine(AbstractCommandline):
    """Base Commandline object for EMBOSS wrappers (PRIVATE).

    This is provided for subclassing, it deals with shared options
    common to all the EMBOSS tools:

    Attributes:
     - auto               Turn off prompts
     - stdout             Write standard output
     - filter             Read standard input, write standard output
     - options            Prompt for standard and additional values
     - debug              Write debug output to program.dbg
     - verbose            Report some/full command line options
     - help               Report command line options. More
                          information on associated and general
                          qualifiers can be found with -help -verbose
     - warning            Report warnings
     - error              Report errors
     - fatal              Report fatal errors
     - die                Report dying program messages

    """

    def __init__(self, cmd=None, **kwargs):
        assert cmd is not None
        extra_parameters = [
            _Switch(
                ["-auto", "auto"],
                "Turn off prompts.\n\n"
                "Automatic mode disables prompting, so we recommend you set this "
                "argument all the time when calling an EMBOSS tool from Biopython.",
            ),
            _Switch(["-stdout", "stdout"], "Write standard output."),
            _Switch(
                ["-filter", "filter"], "Read standard input, write standard output."
            ),
            _Switch(
                ["-options", "options"],
                "Prompt for standard and additional values.\n\n"
                "If you are calling an EMBOSS tool from within Biopython, "
                "we DO NOT recommend using this option.",
            ),
            _Switch(["-debug", "debug"], "Write debug output to program.dbg."),
            _Switch(["-verbose", "verbose"], "Report some/full command line options"),
            _Switch(
                ["-help", "help"],
                "Report command line options.\n\n"
                "More information on associated and general qualifiers "
                "can be found with -help -verbose",
            ),
            _Switch(["-warning", "warning"], "Report warnings."),
            _Switch(["-error", "error"], "Report errors."),
            _Switch(["-die", "die"], "Report dying program messages."),
        ]
        try:
            # Insert extra parameters - at the start just in case there
            # are any arguments which must come last:
            self.parameters = extra_parameters + self.parameters
        except AttributeError:
            # Should we raise an error?  The subclass should have set this up!
            self.parameters = extra_parameters
        AbstractCommandline.__init__(self, cmd, **kwargs)


class _EmbossCommandLine(_EmbossMinimalCommandLine):
    """Base Commandline object for EMBOSS wrappers (PRIVATE).

    This is provided for subclassing, it deals with shared options
    common to all the EMBOSS tools plus:

     - outfile            Output filename

    """

    def __init__(self, cmd=None, **kwargs):
        assert cmd is not None
        extra_parameters = [
            _Option(["-outfile", "outfile"], "Output filename", filename=True)
        ]
        try:
            # Insert extra parameters - at the start just in case there
            # are any arguments which must come last:
            self.parameters = extra_parameters + self.parameters
        except AttributeError:
            # Should we raise an error?  The subclass should have set this up!
            self.parameters = extra_parameters
        _EmbossMinimalCommandLine.__init__(self, cmd, **kwargs)

    def _validate(self):
        # Check the outfile, filter, or stdout option has been set.
        # We can't simply do this via the required flag for the outfile
        # output - this seems the simplest solution.
        if not (self.outfile or self.filter or self.stdout):
            raise ValueError(
                "You must either set outfile (output filename), "
                "or enable filter or stdout (output to stdout)."
            )
        return _EmbossMinimalCommandLine._validate(self)


class Primer3Commandline(_EmbossCommandLine):
    """Commandline object for the Primer3 interface from EMBOSS.

    The precise set of supported arguments depends on your version of EMBOSS.
    This version accepts arguments current at EMBOSS 6.1.0:

    >>> cline = Primer3Commandline(sequence="mysequence.fas", auto=True, hybridprobe=True)
    >>> cline.explainflag = True
    >>> cline.osizeopt=20
    >>> cline.psizeopt=200
    >>> cline.outfile = "myresults.out"
    >>> cline.bogusparameter = 1967  # Invalid parameter
    Traceback (most recent call last):
        ...
    ValueError: Option name bogusparameter was not found.
    >>> print(cline)
    eprimer3 -auto -outfile=myresults.out -sequence=mysequence.fas -hybridprobe=True -psizeopt=200 -osizeopt=20 -explainflag=True

    """

    def __init__(self, cmd="eprimer3", **kwargs):
        """Initialize the class."""
        self.parameters = [
            _Option(
                ["-sequence", "sequence"],
                "Sequence to choose primers from.",
                is_required=True,
            ),
            _Option(["-task", "task"], "Tell eprimer3 what task to perform."),
            _Option(
                ["-hybridprobe", "hybridprobe"],
                "Find an internal oligo to use as a hyb probe.",
            ),
            _Option(
                ["-numreturn", "numreturn"], "Maximum number of primer pairs to return."
            ),
            _Option(
                ["-includedregion", "includedregion"],
                "Subregion of the sequence in which to pick primers.",
            ),
            _Option(["-target", "target"], "Sequence to target for flanking primers."),
            _Option(
                ["-excludedregion", "excludedregion"],
                "Regions to exclude from primer picking.",
            ),
            _Option(
                ["-forwardinput", "forwardinput"],
                "Sequence of a forward primer to check.",
            ),
            _Option(
                ["-reverseinput", "reverseinput"],
                "Sequence of a reverse primer to check.",
            ),
            _Option(
                ["-gcclamp", "gcclamp"],
                "The required number of Gs and Cs at the 3' of each primer.",
            ),
            _Option(["-osize", "osize"], "Optimum length of a primer oligo."),
            _Option(["-minsize", "minsize"], "Minimum length of a primer oligo."),
            _Option(["-maxsize", "maxsize"], "Maximum length of a primer oligo."),
            _Option(
                ["-otm", "otm"],
                "Melting temperature for primer oligo (OBSOLETE).\n\n"
                "Option replaced in EMBOSS 6.6.0 by -opttm",
            ),
            _Option(
                ["-opttm", "opttm"],
                "Optimum melting temperature for a primer oligo.\n\n"
                "Option added in EMBOSS 6.6.0, replacing -otm",
            ),
            _Option(
                ["-mintm", "mintm"], "Minimum melting temperature for a primer oligo."
            ),
            _Option(
                ["-maxtm", "maxtm"], "Maximum melting temperature for a primer oligo."
            ),
            _Option(
                ["-maxdifftm", "maxdifftm"],
                "Maximum difference in melting temperatures between "
                "forward and reverse primers.",
            ),
            _Option(["-ogcpercent", "ogcpercent"], "Optimum GC% for a primer."),
            _Option(["-mingc", "mingc"], "Minimum GC% for a primer."),
            _Option(["-maxgc", "maxgc"], "Maximum GC% for a primer."),
            _Option(
                ["-saltconc", "saltconc"], "Millimolar salt concentration in the PCR."
            ),
            _Option(
                ["-dnaconc", "dnaconc"],
                "Nanomolar concentration of annealing oligos in the PCR.",
            ),
            _Option(
                ["-maxpolyx", "maxpolyx"],
                "Maximum allowable mononucleotide repeat length in a primer.",
            ),
            # Primer length:
            _Option(["-psizeopt", "psizeopt"], "Optimum size for the PCR product."),
            _Option(
                ["-prange", "prange"], "Acceptable range of length for the PCR product."
            ),
            # Primer temperature:
            _Option(
                ["-ptmopt", "ptmopt"],
                "Optimum melting temperature for the PCR product.",
            ),
            _Option(
                ["-ptmmin", "ptmmin"],
                "Minimum allowed melting temperature for the amplicon.",
            ),
            _Option(
                ["-ptmmax", "ptmmax"],
                "Maximum allowed melting temperature for the amplicon.",
            ),
            # Note to self, should be -oexcludedregion not -oexcluderegion
            _Option(
                ["-oexcludedregion", "oexcludedregion"],
                "Do not pick internal oligos in this region.",
            ),
            _Option(["-oligoinput", "oligoinput"], "Sequence of the internal oligo."),
            # Oligo length:
            _Option(["-osizeopt", "osizeopt"], "Optimum length of internal oligo."),
            _Option(["-ominsize", "ominsize"], "Minimum length of internal oligo."),
            _Option(["-omaxsize", "omaxsize"], "Maximum length of internal oligo."),
            # Oligo GC temperature:
            _Option(
                ["-otmopt", "otmopt"], "Optimum melting temperature of internal oligo."
            ),
            _Option(
                ["-otmmin", "otmmin"], "Minimum melting temperature of internal oligo."
            ),
            _Option(
                ["-otmmax", "otmmax"], "Maximum melting temperature of internal oligo."
            ),
            # Oligo GC percent:
            _Option(["-ogcopt", "ogcopt"], "Optimum GC% for internal oligo."),
            _Option(["-ogcmin", "ogcmin"], "Minimum GC% for internal oligo."),
            _Option(["-ogcmax", "ogcmax"], "Maximum GC% for internal oligo."),
            # Oligo salt concentration:
            _Option(
                ["-osaltconc", "osaltconc"],
                "Millimolar concentration of salt in the hybridisation.",
            ),
            _Option(
                ["-odnaconc", "odnaconc"],
                "Nanomolar concentration of internal oligo in the hybridisation.",
            ),
            # Oligo self complementarity
            _Option(
                ["-oanyself", "oanyself"],
                "Maximum allowable alignment score for self-complementarity.",
            ),
            _Option(
                ["-oendself", "oendself"],
                "Max 3'-anchored self-complementarity global alignment score.",
            ),
            _Option(
                ["-opolyxmax", "opolyxmax"],
                "Maximum length of mononucleotide repeat in internal oligo.",
            ),
            _Option(
                ["-mispriminglibraryfile", "mispriminglibraryfile"],
                "File containing library of sequences to avoid amplifying",
            ),
            _Option(
                ["-maxmispriming", "maxmispriming"],
                "Maximum allowed similarity of primers to sequences in "
                "library specified by -mispriminglibrary",
            ),
            _Option(
                ["-omishybmax", "omishybmax"],
                "Maximum alignment score for hybridisation of internal oligo to "
                "library specified by -mishyblibraryfile.",
            ),
            _Option(
                ["-mishyblibraryfile", "mishyblibraryfile"],
                "Library file of seqs to avoid internal oligo hybridisation.",
            ),
            _Option(
                ["-explainflag", "explainflag"],
                "Produce output tags with eprimer3 statistics",
            ),
        ]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class PrimerSearchCommandline(_EmbossCommandLine):
    """Commandline object for the primersearch program from EMBOSS."""

    def __init__(self, cmd="primersearch", **kwargs):
        """Initialize the class."""
        self.parameters = [
            _Option(
                ["-seqall", "-sequences", "sequences", "seqall"],
                "Sequence to look for the primer pairs in.",
                is_required=True,
            ),
            # When this wrapper was written primersearch used -sequences
            # as the argument name. Since at least EMBOSS 5.0 (and
            # perhaps earlier) this has been -seqall instead.
            _Option(
                ["-infile", "-primers", "primers", "infile"],
                "File containing the primer pairs to search for.",
                filename=True,
                is_required=True,
            ),
            # When this wrapper was written primersearch used -primers
            # as the argument name. Since at least EMBOSS 5.0 (and
            # perhaps earlier) this has been -infile instead.
            _Option(
                ["-mismatchpercent", "mismatchpercent"],
                "Allowed percentage mismatch (any integer value, default 0).",
                is_required=True,
            ),
            _Option(
                ["-snucleotide", "snucleotide"], "Sequences are nucleotide (boolean)"
            ),
            _Option(["-sprotein", "sprotein"], "Sequences are protein (boolean)"),
        ]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class FDNADistCommandline(_EmbossCommandLine):
    """Commandline object for the fdnadist program from EMBOSS.

    fdnadist is an EMBOSS wrapper for the PHYLIP program dnadist for
    calulating distance matrices from DNA sequence files.
    """

    def __init__(self, cmd="fdnadist", **kwargs):
        """Initialize the class."""
        self.parameters = [
            _Option(
                ["-sequence", "sequence"],
                "seq file to use (phylip)",
                filename=True,
                is_required=True,
            ),
            _Option(["-method", "method"], "sub. model [f,k,j,l,s]", is_required=True),
            _Option(["-gamma", "gamma"], "gamma [g, i,n]"),
            _Option(
                ["-ncategories", "ncategories"], "number of rate catergories (1-9)"
            ),
            _Option(["-rate", "rate"], "rate for each category"),
            _Option(
                ["-categories", "categories"], "File of substitution rate categories"
            ),
            _Option(["-weights", "weights"], "weights file"),
            _Option(
                ["-gammacoefficient", "gammacoefficient"], "value for gamma (> 0.001)"
            ),
            _Option(["-invarfrac", "invarfrac"], "proportoin of invariant sites"),
            _Option(["-ttratio", "ttratio"], "ts/tv ratio"),
            _Option(["-freqsfrom", "freqsfrom"], "use emprical base freqs"),
            _Option(["-basefreq", "basefreq"], "specify basefreqs"),
            _Option(["-lower", "lower"], "lower triangle matrix (y/N)"),
        ]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class FTreeDistCommandline(_EmbossCommandLine):
    """Commandline object for the ftreedist program from EMBOSS.

    ftreedist is an EMBOSS wrapper for the PHYLIP program treedist used for
    calulating distance measures between phylogentic trees.
    """

    def __init__(self, cmd="ftreedist", **kwargs):
        """Initialize the class."""
        self.parameters = [
            _Option(
                ["-intreefile", "intreefile"],
                "tree file to score (phylip)",
                filename=True,
                is_required=True,
            ),
            _Option(["-dtype", "dtype"], "distance type ([S]ymetric, [b]ranch score)"),
            _Option(
                ["-pairing", "pairing"],
                "tree pairing method ([A]djacent pairs, all [p]ossible pairs)",
            ),
            _Option(["-style", "style"], "output style - [V]erbose, [f]ill, [s]parse"),
            _Option(["-noroot", "noroot"], "treat trees as rooted [N/y]"),
            _Option(
                ["-outgrno", "outgrno"],
                "which taxon to root the trees with (starts from 0)",
            ),
        ]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class FNeighborCommandline(_EmbossCommandLine):
    """Commandline object for the fneighbor program from EMBOSS.

    fneighbor is an EMBOSS wrapper for the PHYLIP program neighbor used for
    calulating neighbor-joining or UPGMA trees from distance matrices.
    """

    def __init__(self, cmd="fneighbor", **kwargs):
        """Initialize the class."""
        self.parameters = [
            _Option(
                ["-datafile", "datafile"],
                "dist file to use (phylip)",
                filename=True,
                is_required=True,
            ),
            _Option(
                ["-matrixtype", "matrixtype"],
                "is martrix [S]quare pr [u]pper or [l]ower",
            ),
            _Option(["-treetype", "treetype"], "nj or UPGMA tree (n/u)"),
            _Option(["-outgrno", "outgrno"], "taxon to use as OG"),
            _Option(["-jumble", "jumble"], "randommise input order (Y/n)"),
            _Option(["-seed", "seed"], "provide a random seed"),
            _Option(["-trout", "trout"], "write tree (Y/n)"),
            _Option(["-outtreefile", "outtreefile"], "filename for output tree"),
            _Option(["-progress", "progress"], "print progress (Y/n)"),
            _Option(["-treeprint", "treeprint"], "print tree (Y/n)"),
        ]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class FSeqBootCommandline(_EmbossCommandLine):
    """Commandline object for the fseqboot program from EMBOSS.

    fseqboot is an EMBOSS wrapper for the PHYLIP program seqboot used to
    pseudo-sample alignment files.
    """

    def __init__(self, cmd="fseqboot", **kwargs):
        """Initialize the class."""
        self.parameters = [
            _Option(
                ["-sequence", "sequence"],
                "seq file to sample (phylip)",
                filename=True,
                is_required=True,
            ),
            _Option(["-categories", "catergories"], "file of input categories"),
            _Option(["-weights", "weights"], " weights file"),
            _Option(["-test", "test"], "specify operation, default is bootstrap"),
            _Option(["-regular", "regular"], "absolute number to resample"),
            _Option(["-fracsample", "fracsample"], "fraction to resample"),
            _Option(
                ["-rewriteformat", "rewriteformat"],
                "output format ([P]hyilp, [n]exus, [x]ml",
            ),
            _Option(["-seqtype", "seqtype"], "output format ([D]na, [p]rotein, [r]na"),
            _Option(["-blocksize", "blocksize"], "print progress (Y/n)"),
            _Option(["-reps", "reps"], "how many replicates, defaults to 100)"),
            _Option(
                ["-justweights", "jusweights"],
                "what to write out [D]atasets of just [w]eights",
            ),
            _Option(["-seed", "seed"], "specify random seed"),
            _Option(["-dotdiff", "dotdiff"], "Use dot-differencing? [Y/n]"),
        ]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class FDNAParsCommandline(_EmbossCommandLine):
    """Commandline object for the fdnapars program from EMBOSS.

    fdnapars is an EMBOSS version of the PHYLIP program dnapars, for
    estimating trees from DNA sequences using parsiomny. Calling this command
    without providing a value for the option "-intreefile" will invoke
    "interactive mode" (and as a result fail if called with subprocess) if
    "-auto" is not set to true.
    """

    def __init__(self, cmd="fdnapars", **kwargs):
        """Initialize the class."""
        self.parameters = [
            _Option(
                ["-sequence", "sequence"],
                "seq file to use (phylip)",
                filename=True,
                is_required=True,
            ),
            _Option(["-intreefile", "intreefile"], "Phylip tree file"),
            _Option(["-weights", "weights"], "weights file"),
            _Option(["-maxtrees", "maxtrees"], "max trees to save during run"),
            _Option(["-thorough", "thorough"], "more thorough search (Y/n)"),
            _Option(["-rearrange", "rearrange"], "Rearrange on just 1 best tree (Y/n)"),
            _Option(
                ["-transversion", "transversion"], "Use tranversion parsimony (y/N)"
            ),
            _Option(
                ["-njumble", "njumble"],
                "number of times to randomise input order (default is 0)",
            ),
            _Option(["-seed", "seed"], "provide random seed"),
            _Option(["-outgrno", "outgrno"], "Specify outgroup"),
            _Option(["-thresh", "thresh"], "Use threshold parsimony (y/N)"),
            _Option(["-threshold", "threshold"], "Threshold value"),
            _Option(["-trout", "trout"], "Write trees to file (Y/n)"),
            _Option(["-outtreefile", "outtreefile"], "filename for output tree"),
            _Option(["-dotdiff", "dotdiff"], "Use dot-differencing? [Y/n]"),
        ]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class FProtParsCommandline(_EmbossCommandLine):
    """Commandline object for the fdnapars program from EMBOSS.

    fprotpars is an EMBOSS version of the PHYLIP program protpars, for
    estimating trees from protein  sequences using parsiomny. Calling this
    command without providing a value for the option "-intreefile" will invoke
    "interactive mode" (and as a result fail if called with subprocess) if
    "-auto" is not set to true.
    """

    def __init__(self, cmd="fprotpars", **kwargs):
        """Initialize the class."""
        self.parameters = [
            _Option(
                ["-sequence", "sequence"],
                "seq file to use (phylip)",
                filename=True,
                is_required=True,
            ),
            _Option(["-intreefile", "intreefile"], "Phylip tree file to score"),
            _Option(
                ["-outtreefile", "outtreefile"],
                "phylip tree output file",
                filename=True,
                is_required=True,
            ),
            _Option(["-weights", "weights"], "weights file"),
            _Option(["-whichcode", "whichcode"], "which genetic code, [U,M,V,F,Y]]"),
            _Option(
                ["-njumble", "njumble"],
                "number of times to randomise input order (default is 0)",
            ),
            _Option(["-seed", "seed"], "provide random seed"),
            _Option(["-outgrno", "outgrno"], "Specify outgroup"),
            _Option(["-thresh", "thresh"], "Use threshold parsimony (y/N)"),
            _Option(["-threshold", "threshold"], "Threshold value"),
            _Option(["-trout", "trout"], "Write trees to file (Y/n)"),
            _Option(["-dotdiff", "dotdiff"], "Use dot-differencing? [Y/n]"),
        ]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class FProtDistCommandline(_EmbossCommandLine):
    """Commandline object for the fprotdist program from EMBOSS.

    fprotdist is an EMBOSS wrapper for the PHYLIP program protdist used to
    estimate trees from protein sequences using parsimony
    """

    def __init__(self, cmd="fprotdist", **kwargs):
        """Initialize the class."""
        self.parameters = [
            _Option(
                ["-sequence", "sequence"],
                "seq file to use (phylip)",
                filename=True,
                is_required=True,
            ),
            _Option(
                ["-ncategories", "ncategories"], "number of rate catergories (1-9)"
            ),
            _Option(["-rate", "rate"], "rate for each category"),
            _Option(["-catergories", "catergories"], "file of rates"),
            _Option(["-weights", "weights"], "weights file"),
            _Option(["-method", "method"], "sub. model [j,h,d,k,s,c]"),
            _Option(["-gamma", "gamma"], "gamma [g, i,c]"),
            _Option(
                ["-gammacoefficient", "gammacoefficient"], "value for gamma (> 0.001)"
            ),
            _Option(
                ["-invarcoefficient", "invarcoefficient"],
                "float for variation of substitution rate among sites",
            ),
            _Option(["-aacateg", "aacateg"], "Choose the category to use [G,C,H]"),
            _Option(["-whichcode", "whichcode"], "genetic code [c,m,v,f,y]"),
            _Option(["-ease", "ease"], "Pob change catergory (float between -0 and 1)"),
            _Option(["-ttratio", "ttratio"], "Transition/transversion ratio (0-1)"),
            _Option(
                ["-basefreq", "basefreq"], "DNA base frequencies (space separated list)"
            ),
        ]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class FConsenseCommandline(_EmbossCommandLine):
    """Commandline object for the fconsense program from EMBOSS.

    fconsense is an EMBOSS wrapper for the PHYLIP program consense used to
    calculate consensus trees.
    """

    def __init__(self, cmd="fconsense", **kwargs):
        """Initialize the class."""
        self.parameters = [
            _Option(
                ["-intreefile", "intreefile"],
                "file with phylip trees to make consensus from",
                filename=True,
                is_required=True,
            ),
            _Option(["-method", "method"], "consensus method [s, mr, MRE, ml]"),
            _Option(
                ["-mlfrac", "mlfrac"],
                "cut-off freq for branch to appear in consensus (0.5-1.0)",
            ),
            _Option(["-root", "root"], "treat trees as rooted (YES, no)"),
            _Option(["-outgrno", "outgrno"], "OTU to use as outgroup (starts from 0)"),
            _Option(["-trout", "trout"], "treat trees as rooted (YES, no)"),
            _Option(
                ["-outtreefile", "outtreefile"], "Phylip tree output file (optional)"
            ),
        ]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class WaterCommandline(_EmbossCommandLine):
    """Commandline object for the water program from EMBOSS."""

    def __init__(self, cmd="water", **kwargs):
        """Initialize the class."""
        self.parameters = [
            _Option(
                ["-asequence", "asequence"],
                "First sequence to align",
                filename=True,
                is_required=True,
            ),
            _Option(
                ["-bsequence", "bsequence"],
                "Second sequence to align",
                filename=True,
                is_required=True,
            ),
            _Option(["-gapopen", "gapopen"], "Gap open penalty", is_required=True),
            _Option(
                ["-gapextend", "gapextend"], "Gap extension penalty", is_required=True
            ),
            _Option(["-datafile", "datafile"], "Matrix file", filename=True),
            _Switch(
                ["-nobrief", "nobrief"], "Display extended identity and similarity"
            ),
            _Switch(["-brief", "brief"], "Display brief identity and similarity"),
            _Option(
                ["-similarity", "similarity"], "Display percent identity and similarity"
            ),
            _Option(
                ["-snucleotide", "snucleotide"], "Sequences are nucleotide (boolean)"
            ),
            _Option(["-sprotein", "sprotein"], "Sequences are protein (boolean)"),
            _Option(
                ["-aformat", "aformat"],
                "Display output in a different specified output format",
            ),
        ]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class NeedleCommandline(_EmbossCommandLine):
    """Commandline object for the needle program from EMBOSS."""

    def __init__(self, cmd="needle", **kwargs):
        """Initialize the class."""
        self.parameters = [
            _Option(
                ["-asequence", "asequence"],
                "First sequence to align",
                filename=True,
                is_required=True,
            ),
            _Option(
                ["-bsequence", "bsequence"],
                "Second sequence to align",
                filename=True,
                is_required=True,
            ),
            _Option(["-gapopen", "gapopen"], "Gap open penalty", is_required=True),
            _Option(
                ["-gapextend", "gapextend"], "Gap extension penalty", is_required=True
            ),
            _Option(["-datafile", "datafile"], "Matrix file", filename=True),
            _Option(["-endweight", "endweight"], "Apply And gap penalties"),
            _Option(
                ["-endopen", "endopen"],
                "The score taken away when an end gap is created.",
            ),
            _Option(
                ["-endextend", "endextend"],
                "The score added to the end gap penality for each base or "
                "residue in the end gap.",
            ),
            _Switch(
                ["-nobrief", "nobrief"], "Display extended identity and similarity"
            ),
            _Switch(["-brief", "brief"], "Display brief identity and similarity"),
            _Option(
                ["-similarity", "similarity"], "Display percent identity and similarity"
            ),
            _Option(
                ["-snucleotide", "snucleotide"], "Sequences are nucleotide (boolean)"
            ),
            _Option(["-sprotein", "sprotein"], "Sequences are protein (boolean)"),
            _Option(
                ["-aformat", "aformat"],
                "Display output in a different specified output format",
            ),
        ]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class NeedleallCommandline(_EmbossCommandLine):
    """Commandline object for the needleall program from EMBOSS."""

    def __init__(self, cmd="needleall", **kwargs):
        """Initialize the class."""
        self.parameters = [
            _Option(
                ["-asequence", "asequence"],
                "First sequence to align",
                filename=True,
                is_required=True,
            ),
            _Option(
                ["-bsequence", "bsequence"],
                "Second sequence to align",
                filename=True,
                is_required=True,
            ),
            _Option(["-gapopen", "gapopen"], "Gap open penalty", is_required=True),
            _Option(
                ["-gapextend", "gapextend"], "Gap extension penalty", is_required=True
            ),
            _Option(["-datafile", "datafile"], "Matrix file", filename=True),
            _Option(
                ["-minscore", "minscore"],
                "Exclude alignments with scores below this threshold score.",
            ),
            _Option(["-errorfile", "errorfile"], "Error file to be written to."),
            _Option(["-endweight", "endweight"], "Apply And gap penalties"),
            _Option(
                ["-endopen", "endopen"],
                "The score taken away when an end gap is created.",
            ),
            _Option(
                ["-endextend", "endextend"],
                "The score added to the end gap penality for each base or "
                "residue in the end gap.",
            ),
            _Switch(
                ["-nobrief", "nobrief"], "Display extended identity and similarity"
            ),
            _Switch(["-brief", "brief"], "Display brief identity and similarity"),
            _Option(
                ["-similarity", "similarity"], "Display percent identity and similarity"
            ),
            _Option(
                ["-snucleotide", "snucleotide"], "Sequences are nucleotide (boolean)"
            ),
            _Option(["-sprotein", "sprotein"], "Sequences are protein (boolean)"),
            _Option(
                ["-aformat", "aformat"],
                "Display output in a different specified output format",
            ),
        ]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class StretcherCommandline(_EmbossCommandLine):
    """Commandline object for the stretcher program from EMBOSS."""

    def __init__(self, cmd="stretcher", **kwargs):
        """Initialize the class."""
        self.parameters = [
            _Option(
                ["-asequence", "asequence"],
                "First sequence to align",
                filename=True,
                is_required=True,
            ),
            _Option(
                ["-bsequence", "bsequence"],
                "Second sequence to align",
                filename=True,
                is_required=True,
            ),
            _Option(
                ["-gapopen", "gapopen"],
                "Gap open penalty",
                is_required=True,
                checker_function=lambda value: isinstance(value, int),
            ),
            _Option(
                ["-gapextend", "gapextend"],
                "Gap extension penalty",
                is_required=True,
                checker_function=lambda value: isinstance(value, int),
            ),
            _Option(["-datafile", "datafile"], "Matrix file", filename=True),
            _Option(
                ["-snucleotide", "snucleotide"], "Sequences are nucleotide (boolean)"
            ),
            _Option(["-sprotein", "sprotein"], "Sequences are protein (boolean)"),
            _Option(
                ["-aformat", "aformat"],
                "Display output in a different specified output format",
            ),
        ]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class FuzznucCommandline(_EmbossCommandLine):
    """Commandline object for the fuzznuc program from EMBOSS."""

    def __init__(self, cmd="fuzznuc", **kwargs):
        """Initialize the class."""
        self.parameters = [
            _Option(
                ["-sequence", "sequence"], "Sequence database USA", is_required=True
            ),
            _Option(
                ["-pattern", "pattern"],
                "Search pattern, using standard IUPAC one-letter codes",
                is_required=True,
            ),
            _Option(["-pmismatch", "pmismatch"], "Number of mismatches"),
            _Option(["-complement", "complement"], "Search complementary strand"),
            _Option(["-rformat", "rformat"], "Specify the report format to output in."),
        ]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class FuzzproCommandline(_EmbossCommandLine):
    """Commandline object for the fuzzpro program from EMBOSS."""

    def __init__(self, cmd="fuzzpro", **kwargs):
        """Initialize the class."""
        self.parameters = [
            _Option(
                ["-sequence", "sequence"], "Sequence database USA", is_required=True
            ),
            _Option(
                ["-pattern", "pattern"],
                "Search pattern, using standard IUPAC one-letter codes",
                is_required=True,
            ),
            _Option(["-pmismatch", "pmismatch"], "Number of mismatches"),
            _Option(["-rformat", "rformat"], "Specify the report format to output in."),
        ]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class Est2GenomeCommandline(_EmbossCommandLine):
    """Commandline object for the est2genome program from EMBOSS."""

    def __init__(self, cmd="est2genome", **kwargs):
        """Initialize the class."""
        self.parameters = [
            _Option(["-est", "est"], "EST sequence(s)", is_required=True),
            _Option(["-genome", "genome"], "Genomic sequence", is_required=True),
            _Option(["-match", "match"], "Score for matching two bases"),
            _Option(["-mismatch", "mismatch"], "Cost for mismatching two bases"),
            _Option(
                ["-gappenalty", "gappenalty"],
                "Cost for deleting a single base in either sequence, "
                "excluding introns",
            ),
            _Option(
                ["-intronpenalty", "intronpenalty"],
                "Cost for an intron, independent of length.",
            ),
            _Option(
                ["-splicepenalty", "splicepenalty"],
                "Cost for an intron, independent of length "
                "and starting/ending on donor-acceptor sites",
            ),
            _Option(
                ["-minscore", "minscore"],
                "Exclude alignments with scores below this threshold score.",
            ),
            _Option(
                ["-reverse", "reverse"], "Reverse the orientation of the EST sequence"
            ),
            _Option(["-splice", "splice"], "Use donor and acceptor splice sites."),
            _Option(
                ["-mode", "mode"],
                "This determines the comparion mode. 'both', 'forward', or 'reverse'",
            ),
            _Option(
                ["-best", "best"],
                "You can print out all comparisons instead of just the best",
            ),
            _Option(["-space", "space"], "for linear-space recursion."),
            _Option(["-shuffle", "shuffle"], "Shuffle"),
            _Option(["-seed", "seed"], "Random number seed"),
            _Option(["-align", "align"], "Show the alignment."),
            _Option(["-width", "width"], "Alignment width"),
        ]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class ETandemCommandline(_EmbossCommandLine):
    """Commandline object for the etandem program from EMBOSS."""

    def __init__(self, cmd="etandem", **kwargs):
        """Initialize the class."""
        self.parameters = [
            _Option(
                ["-sequence", "sequence"], "Sequence", filename=True, is_required=True
            ),
            _Option(
                ["-minrepeat", "minrepeat"], "Minimum repeat size", is_required=True
            ),
            _Option(
                ["-maxrepeat", "maxrepeat"], "Maximum repeat size", is_required=True
            ),
            _Option(["-threshold", "threshold"], "Threshold score"),
            _Option(["-mismatch", "mismatch"], "Allow N as a mismatch"),
            _Option(["-uniform", "uniform"], "Allow uniform consensus"),
            _Option(["-rformat", "rformat"], "Output report format"),
        ]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class EInvertedCommandline(_EmbossCommandLine):
    """Commandline object for the einverted program from EMBOSS."""

    def __init__(self, cmd="einverted", **kwargs):
        """Initialize the class."""
        self.parameters = [
            _Option(
                ["-sequence", "sequence"], "Sequence", filename=True, is_required=True
            ),
            _Option(["-gap", "gap"], "Gap penalty", filename=True, is_required=True),
            _Option(
                ["-threshold", "threshold"], "Minimum score threshold", is_required=True
            ),
            _Option(["-match", "match"], "Match score", is_required=True),
            _Option(["-mismatch", "mismatch"], "Mismatch score", is_required=True),
            _Option(
                ["-maxrepeat", "maxrepeat"],
                "Maximum separation between the start and end of repeat",
            ),
        ]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class PalindromeCommandline(_EmbossCommandLine):
    """Commandline object for the palindrome program from EMBOSS."""

    def __init__(self, cmd="palindrome", **kwargs):
        """Initialize the class."""
        self.parameters = [
            _Option(
                ["-sequence", "sequence"], "Sequence", filename=True, is_required=True
            ),
            _Option(
                ["-minpallen", "minpallen"],
                "Minimum palindrome length",
                is_required=True,
            ),
            _Option(
                ["-maxpallen", "maxpallen"],
                "Maximum palindrome length",
                is_required=True,
            ),
            _Option(
                ["-gaplimit", "gaplimit"],
                "Maximum gap between repeats",
                is_required=True,
            ),
            _Option(
                ["-nummismatches", "nummismatches"],
                "Number of mismatches allowed",
                is_required=True,
            ),
            _Option(
                ["-overlap", "overlap"], "Report overlapping matches", is_required=True
            ),
        ]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class TranalignCommandline(_EmbossCommandLine):
    """Commandline object for the tranalign program from EMBOSS."""

    def __init__(self, cmd="tranalign", **kwargs):
        """Initialize the class."""
        self.parameters = [
            _Option(
                ["-asequence", "asequence"],
                "Nucleotide sequences to be aligned.",
                filename=True,
                is_required=True,
            ),
            _Option(
                ["-bsequence", "bsequence"],
                "Protein sequence alignment",
                filename=True,
                is_required=True,
            ),
            _Option(
                ["-outseq", "outseq"],
                "Output sequence file.",
                filename=True,
                is_required=True,
            ),
            _Option(["-table", "table"], "Code to use"),
        ]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class DiffseqCommandline(_EmbossCommandLine):
    """Commandline object for the diffseq program from EMBOSS."""

    def __init__(self, cmd="diffseq", **kwargs):
        """Initialize the class."""
        self.parameters = [
            _Option(
                ["-asequence", "asequence"],
                "First sequence to compare",
                filename=True,
                is_required=True,
            ),
            _Option(
                ["-bsequence", "bsequence"],
                "Second sequence to compare",
                filename=True,
                is_required=True,
            ),
            _Option(
                ["-wordsize", "wordsize"],
                "Word size to use for comparisons (10 default)",
                is_required=True,
            ),
            _Option(
                ["-aoutfeat", "aoutfeat"],
                "File for output of first sequence's features",
                filename=True,
                is_required=True,
            ),
            _Option(
                ["-boutfeat", "boutfeat"],
                "File for output of second sequence's features",
                filename=True,
                is_required=True,
            ),
            _Option(["-rformat", "rformat"], "Output report file format"),
        ]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class IepCommandline(_EmbossCommandLine):
    """Commandline for EMBOSS iep: calculated isoelectric point and charge.

    Examples
    --------
    >>> from Bio.Emboss.Applications import IepCommandline
    >>> iep_cline = IepCommandline(sequence="proteins.faa",
    ...                            outfile="proteins.txt")
    >>> print(iep_cline)
    iep -outfile=proteins.txt -sequence=proteins.faa

    You would typically run the command line with iep_cline() or via the
    Python subprocess module, as described in the Biopython tutorial.

    """

    def __init__(self, cmd="iep", **kwargs):
        """Initialize the class."""
        self.parameters = [
            _Option(
                ["-sequence", "sequence"],
                "Protein sequence(s) filename",
                filename=True,
                is_required=True,
            ),
            _Option(
                ["-amino", "amino"],
                """Number of N-termini

                    Integer 0 (default) or more.
                    """,
            ),
            _Option(
                ["-carboxyl", "carboxyl"],
                """Number of C-termini

                    Integer 0 (default) or more.
                    """,
            ),
            _Option(
                ["-lysinemodified", "lysinemodified"],
                """Number of modified lysines

                    Integer 0 (default) or more.
                    """,
            ),
            _Option(
                ["-disulphides", "disulphides"],
                """Number of disulphide bridges

                    Integer 0 (default) or more.
                    """,
            ),
            # Should we implement the -termini switch as well?
            _Option(
                ["-notermini", "notermini"],
                "Exclude (True) or include (False) charge at N and C terminus.",
            ),
        ]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


# seqret uses -outseq, not -outfile, so use the base class:
class SeqretCommandline(_EmbossMinimalCommandLine):
    """Commandline object for the seqret program from EMBOSS.

    This tool allows you to interconvert between different sequence file
    formats (e.g. GenBank to FASTA). Combining Biopython's Bio.SeqIO module
    with seqret using a suitable intermediate file format can allow you to
    read/write to an even wider range of file formats.

    This wrapper currently only supports the core functionality, things like
    feature tables (in EMBOSS 6.1.0 onwards) are not yet included.
    """

    def __init__(self, cmd="seqret", **kwargs):
        """Initialize the class."""
        self.parameters = [
            _Option(
                ["-sequence", "sequence"], "Input sequence(s) filename", filename=True
            ),
            _Option(["-outseq", "outseq"], "Output sequence file.", filename=True),
            _Option(
                ["-sformat", "sformat"],
                "Input sequence(s) format (e.g. fasta, genbank)",
            ),
            _Option(
                ["-osformat", "osformat"],
                "Output sequence(s) format (e.g. fasta, genbank)",
            ),
        ]
        _EmbossMinimalCommandLine.__init__(self, cmd, **kwargs)

    def _validate(self):
        # Check the outfile, filter, or stdout option has been set.
        # We can't simply do this via the required flag for the outfile
        # output - this seems the simplest solution.
        if not (self.outseq or self.filter or self.stdout):
            raise ValueError(
                "You must either set outfile (output filename), "
                "or enable filter or stdout (output to stdout)."
            )
        if not (self.sequence or self.filter or self.stdint):
            raise ValueError(
                "You must either set sequence (input filename), "
                "or enable filter or stdin (input from stdin)."
            )
        return _EmbossMinimalCommandLine._validate(self)


class SeqmatchallCommandline(_EmbossCommandLine):
    """Commandline object for the seqmatchall program from EMBOSS.

    e.g.
    >>> cline = SeqmatchallCommandline(sequence="opuntia.fasta", outfile="opuntia.txt")
    >>> cline.auto = True
    >>> cline.wordsize = 18
    >>> cline.aformat = "pair"
    >>> print(cline)
    seqmatchall -auto -outfile=opuntia.txt -sequence=opuntia.fasta -wordsize=18 -aformat=pair

    """

    def __init__(self, cmd="seqmatchall", **kwargs):
        """Initialize the class."""
        self.parameters = [
            _Option(
                ["-sequence", "sequence"],
                "Readable set of sequences",
                filename=True,
                is_required=True,
            ),
            _Option(
                ["-wordsize", "wordsize"], "Word size (Integer 2 or more, default 4)"
            ),
            _Option(
                ["-aformat", "aformat"],
                "Display output in a different specified output format",
            ),
        ]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
