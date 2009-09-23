# Copyright 2001-2009 Brad Chapman.
# Revisions copyright 2009 by Peter Cock.
# Revisions copyright 2009 by David Winter.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Code to interact with and run various EMBOSS programs.

These classes follow the AbstractCommandline interfaces for running
programs.
"""

from Bio.Application import _Option, _Switch, AbstractCommandline

class _EmbossMinimalCommandLine(AbstractCommandline) :
    """Base Commandline object for EMBOSS wrappers (PRIVATE).

    This is provided for subclassing, it deals with shared options
    common to all the EMBOSS tools:

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
        extra_parameters = [\
           _Switch(["-auto","auto"], [],
                   """Turn off prompts.
           
                   Automatic mode disables prompting, so we recommend you set
                   this argument all the time when calling an EMBOSS tool from
                   Biopython.
                   """),
           _Switch(["-stdout","stdout"], [],
                   "Write standard output."),
           _Switch(["-filter","filter"], [],
                   "Read standard input, write standard output."),
           _Switch(["-options","options"], [],
                   """Prompt for standard and additional values.

                   If you are calling an EMBOSS tool from within Biopython,
                   we DO NOT recommend using this option.
                   """),
           _Switch(["-debug","debug"], [],
                   "Write debug output to program.dbg."),
           _Switch(["-verbose","verbose"], [],
                   "Report some/full command line options"),
           _Switch(["-help","help"], [],
                   """Report command line options.

                   More information on associated and general qualifiers can
                   be found with -help -verbose
                   """),
           _Switch(["-warning","warning"], [],
                   "Report warnings."),
           _Switch(["-error","error"], [],
                   "Report errors."),
           _Switch(["-die","die"], [],
                   "Report dying program messages."),
            ]
        try :
            #Insert extra parameters - at the start just in case there
            #are any arguments which must come last:
            self.parameters = extra_parameters + self.parameters
        except AttributeError:
            #Should we raise an error?  The subclass should have set this up!
            self.parameters = extra_parameters
        AbstractCommandline.__init__(self, cmd, **kwargs)

class _EmbossCommandLine(_EmbossMinimalCommandLine) :
    """Base Commandline object for EMBOSS wrappers (PRIVATE).

    This is provided for subclassing, it deals with shared options
    common to all the EMBOSS tools plus:
    
     - outfile            Output filename

    """
    def __init__(self, cmd=None, **kwargs):
        assert cmd is not None
        extra_parameters = [\
           _Option(["-outfile","outfile"], ["output", "file"], None, 0,
                   "Output filename"),
            ]
        try :
            #Insert extra parameters - at the start just in case there
            #are any arguments which must come last:
            self.parameters = extra_parameters + self.parameters
        except AttributeError:
            #Should we raise an error?  The subclass should have set this up!
            self.parameters = extra_parameters
        _EmbossMinimalCommandLine.__init__(self, cmd, **kwargs)

    def _validate(self) :
        #Check the outfile, filter, or stdout option has been set.
        #We can't simply do this via the required flag for the outfile
        #output - this seems the simplest solution.
        if not (self.outfile or self.filter or self.stdout) :
            raise ValueError("You must either set outfile (output filename), "
                             "or enable filter or stdout (output to stdout).")
        return _EmbossMinimalCommandLine._validate(self)

        
class Primer3Commandline(_EmbossCommandLine):
    """Commandline object for the Primer3 interface from EMBOSS.
    """
    def __init__(self, cmd="eprimer3", **kwargs):
        self.parameters = \
          [_Option(["-sequence","sequence"], ["input"], None, 1,
                   "Sequence to choose primers from"),
           _Option(["-task","task"], ["input"], None, 0),
           _Option(["-numreturn","numreturn"], ["input"], None, 0),
           _Option(["-includedregion","includedregion"], ["input"], None, 0),
           _Option(["-target","target"], ["input"], None, 0),
           _Option(["-excludedregion","excludedregion"], ["input"], None, 0),
           _Option(["-forwardinput","forwardinput"], ["input"], None, 0),
           _Option(["-reverseinput","reverseinput"], ["input"], None, 0),
           _Option(["-gcclamp","gcclamp"], ["input"], None, 0),
           _Option(["-osize","osize"], ["input"], None, 0),
           _Option(["-minsize","minsize"], ["input"], None, 0),
           _Option(["-maxsize","maxsize"], ["input"], None, 0),
           _Option(["-otm","otm"], ["input"], None, 0),
           _Option(["-mintm","mintm"], ["input"], None, 0),
           _Option(["-maxtm","maxtm"], ["input"], None, 0),
           _Option(["-maxdifftm","maxdifftm"], ["input"], None, 0),
           _Option(["-ogcpercent","ogcpercent"], ["input"], None, 0),
           _Option(["-mingc","mingc"], ["input"], None, 0),
           _Option(["-maxgc","maxgc"], ["input"], None, 0),
           _Option(["-saltconc","saltconc"], ["input"], None, 0),
           _Option(["-dnaconc","dnaconc"], ["input"], None, 0),
           _Option(["-maxployx","maxployx"], ["input"], None, 0),
           _Option(["-productosize","productosize"], ["input"], None, 0),
           _Option(["-productsizerange","productsizerange"], ["input"], None, 0),
           _Option(["-productotm","productotm"], ["input"], None, 0),
           _Option(["-productmintm","productmintm"], ["input"], None, 0),
           _Option(["-productmaxtm","productmaxtm"], ["input"], None, 0),
           _Option(["-oligoexcluderegion","oligoexcluderegion"], ["input"], None, 0),
           _Option(["-oligoinput","oligoinput"], ["input"], None, 0),
           _Option(["-oligosize","oligosize"], ["input"], None, 0),
           _Option(["-oligominsize","oligominsize"], ["input"], None, 0),
           _Option(["-oligomaxsize","oligomaxsize"], ["input"], None, 0),
           _Option(["-oligotm","oligotm"], ["input"], None, 0),
           _Option(["-oligomintm","oligomintm"], ["input"], None, 0),
           _Option(["-oligomaxtm","oligomaxtm"], ["input"], None, 0),
           _Option(["-oligoogcpercent","oligoogcpercent"], ["input"], None, 0),
           _Option(["-oligomingc","oligomingc"], ["input"], None, 0),
           _Option(["-oligomaxgc","oligomaxgc"], ["input"], None, 0),
           _Option(["-oligosaltconc","oligosaltconc"], ["input"], None, 0),
           _Option(["-oligodnaconc","oligodnaconc"], ["input"], None, 0),
           _Option(["-oligoselfany","oligoselfany"], ["input"], None, 0),
           _Option(["-oligoselfend","oligoselfend"], ["input"], None, 0),
           _Option(["-oligomaxpolyx","oligomaxpolyx"], ["input"], None, 0),
           _Option(["-mispriminglibraryfile","mispriminglibraryfile"], ["input"], None, 0),
           _Option(["-maxmispriming","maxmispriming"], ["input"], None, 0),
           _Option(["-oligomishyblibraryfile","oligomishyblibraryfile"], ["input"], None, 0),
           _Option(["-oligomaxmishyb","oligomaxmishyb"], ["input"], None, 0),
           _Option(["-explainflag","explainflag"], ["input"], None, 0),
           ]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class PrimerSearchCommandline(_EmbossCommandLine):
    """Commandline object for the primersearch program from EMBOSS.
    """
    def __init__(self, cmd="primersearch", **kwargs):
        self.parameters = \
         [_Option(["-sequences","sequences"], ["input"], None, 1,
                  "Sequence to look for the primer pairs in."),
          _Option(["-primers","primers"], ["input", "file"], None, 1,
                  "File containing the primer pairs to search for."),
          #Including -out and out for backwards compatibility only!
          #_Option(["-outfile","-out","out","outfile"], ["output", "file"], None, 0,
          #        "Name of the output file."),
          _Option(["-mismatchpercent","mismatchpercent"], ["input"], None, 1,
                  "Allowed percentage mismatch.")]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)

    def set_parameter(self, name, value=None) :
        #Due to a historical inconsistency, the PrimerSearchCommandline
        #wrapper used -out and out, rather than -output and output like all
        #the other EMBOSS wrappers.  I want to implement this paramter via
        #the common _EmbossCommandLine base class, hence this hack for
        #backwards compatibility:
        if name in ["out", "-out"] :
            import warnings
            warnings.warn('Aliases "-out" and "out" are deprecated, please use '
                          'either "-outfile" or "outfile" with set_parameter '
                          'instead, or use the outfile property.',
                          DeprecationWarning)
            name = "outfile"
        _EmbossCommandLine.set_parameter(self, name, value)


class EProtDistCommandline(_EmbossCommandLine):
    """Commandline object for the eprotdist program from EMBOSS (OBSOLETE).

    This is an EMBOSS wrapper around protdist from PHYLIP.

    It has been replaced by "fprotdist", see FProtDistCommandline.
    """
    def __init__(self, cmd="eprotdist", **kwargs):
        self.parameters = \
         [_Option(["-msf","msf"], ["input"], None, 1,
                  "File containing sequences"),
          _Option(["-method","method"], ["input"], None, 1,
                  "Choose the method to use"),
          _Option(["-categ","categ"], ["input"], None, 0,
                  "Choose the category to use"),
          _Option(["-gencode","gencode"], ["input"], None, 0,
                  "Which genetic code"),
          _Option(["-prob","prob"], ["input"], None, 0,
                  "Prob change category (1.0=easy)"),
          _Option(["-tranrate","tranrate"], ["input"], None, 0,
                  "Transition/transversion ratio"),
          _Option(["-freqa","freqa"], ["input"], None, 0,
                  "Frequency for A"),
          _Option(["-freqc","freqc"], ["input"], None, 0,
                  "Frequency for C"),
          _Option(["-freqg","freqg"], ["input"], None, 0,
                  "Frequency for G"),
          _Option(["-freqt","freqt"], ["input"], None, 0,
                  "Frequency for T"),
          _Option(["-printdata","printdata"], ["input"], None, 0,
                  "Print out the data at start of run"),
          _Option(["-progress","progress"], ["input"], None, 0,
                  "Print indications of progress of run"),
          _Option(["-basefrequency","basefrequency"], ["input"], None, 0,
                  "Use empirical base frequencies")]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class ENeighborCommandline(_EmbossCommandLine):
    """Commandline object for the eneighbor program from EMBOSS (OBSOLETE).

    This is an EMBOSS wrapper around neighbor from PHYLIP.

    It has been replaced by "fneighbor", see FNeighborCommandline.
    """
    def __init__(self, cmd="eneighbor", **kwargs):
        self.parameters = \
         [_Option(["-infile","infile"], ["input"], None, 1,
                  "infile value"),
          _Option(["-trout","trout"], ["input"], None, 1,
                  "Create a tree file"),
          _Option(["-treefile","treefile"], ["input"], None, 1,
                  "Tree file name"),
          _Option(["-nj","nj"], ["input"], None, 1,
                  "Neighbor-joining"),
          _Option(["-noog","noog"], ["input"], None, 1,
                  "Outgroup root"),
          _Option(["-outgnum","outgnum"], ["input"], None, 0,
                  "number of the outgroup"),
          _Option(["-randseed","randseed"], ["input"], None, 0,
                  "Random number seed (must be odd)"),
          _Option(["-datasets","datasets"], ["input"], None, 0,
                  "How many data sets"),
          _Option(["-drawtree","drawtree"], ["input"], None, 0,
                  "Draw tree"),
          _Option(["-lt","lt"], ["input"], None, 0,
                  "Lower-triangular data matrix"),
          _Option(["-ut","ut"], ["input"], None, 0,
                  "Upper-triangular data matrix"),
          _Option(["-sr","sr"], ["input"], None, 0,
                  "Subreplicates"),
          _Option(["-random","random"], ["input"], None, 0,
                  "Randomize input order of species"),
          _Option(["-multsets","multsets"], ["input"], None, 0,
                  "Analyze multiple data sets"),
          _Option(["-printdata","printdata"], ["input"], None, 0,
                  "Print out the data at start of run"),
          _Option(["-progress","progress"], ["input"], None, 0,
                  "Print indications of progress of run")]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class EProtParsCommandline(_EmbossCommandLine):
    """Commandline object for the eprotpars program from EMBOSS (OBSOLETE).

    This is an EMBOSS wrapper around protpars from PHYLIP.

    It has been replaced by "fprotpars", see FProtParsCommandline.
    """
    def __init__(self, cmd="eprotpars", **kwargs):
        self.parameters = \
         [_Option(["-msf","msf"], ["input", "file"], None, 1,
                  "Sequences file to be read in"),
          _Option(["-besttree","besttree"], ["input"], None, 0,
                  "Search for the best tree"),
          _Option(["-random","random"], ["input"], None, 0,
                  "Randomize input order of species"),
          _Option(["-norandom","norandom"], ["input"], None, 0,
                  "Do not randomize input order of species"),
          _Option(["-randseed","randseed"], ["input"], None, 0,
                  "Random number seed (must be odd)"),
          _Option(["-randtimes","randtimes"], ["input"], None, 0,
                  "How many times to randomize"),
          _Option(["-og","og"], ["input"], None, 0,
                  "Use an outgroup root"),
          _Option(["-noog","noog"], ["input"], None, 0,
                  "Do not use an outgroup root"),
          _Option(["-outgnum","outgnum"], ["input"], None, 0,
                  "Number of the outgroup"),
          _Option(["-thresh","thresh"], ["input"], None, 0,
                  "Use Threshold parsimony"),
          _Option(["-valthresh","valthresh"], ["input"], None, 0,
                  "threshold value"),
          _Option(["-printdata","printdata"], ["input"], None, 0,
                  "Print out the data at start of run"),
          _Option(["-progress","progress"], ["input"], None, 0,
                  "Print indications of progress of run"),
          _Option(["-steps","steps"], ["input"], None, 0,
                  "Print out steps in each site"),
          _Option(["-seqatnodes","seqatnodes"], ["input"], None, 0,
                  "Print sequences at all nodes of tree"),
          _Option(["-drawtree","drawtree"], ["input"], None, 0,
                  "Draw tree"),
          _Option(["-trout","trout"], ["input"], None, 0,
                  "Create a tree file"),
          _Option(["-notrout","notrout"], ["input"], None, 0,
                  "Do not create a tree file"),
          _Option(["-treefile","treefile"], ["output", "file"], None, 0,
                  "Output treefile name")]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class EConsenseCommandline(_EmbossCommandLine):
    """Commandline object for the econsense program from EMBOSS (OBSOLETE).

    This is an EMBOSS wrapper around consense from PHYLIP.

    It has been replaced by "fconsense", see FConsenseCommandline.
    """
    def __init__(self, cmd="econsense", **kwargs):
        self.parameters = \
         [_Option(["-infile","infile"], ["input", "file"], None, 1,
                  "file to read in (New Hampshire standard form)"),
          _Option(["-notrout","notrout"], ["input"], None, 0,
                  "Do not create a tree file"),
          _Option(["-trout","trout"], ["input"], None, 0,
                  "Create a tree file"),
          _Option(["-treefile","treefile"], ["output", "file"], None, 0,
                  "tree file name"),
          _Option(["-noog","noog"], ["input"], None, 0,
                  "Do not use an outgroup"),
          _Option(["-og","og"], ["input"], None, 0,
                  "Use an outgroup"),
          _Option(["-outgnum","outgnum"], ["input"], None, 0,
                  "number of the outgroup"),
          _Option(["-nodrawtree","nodrawtree"], ["input"], None, 0,
                  "Do not draw a tree"),
          _Option(["-drawtree","drawtree"], ["input"], None, 0,
                  "Draw tree"),
          _Option(["-root","root"], ["input"], None, 0,
                  "Trees to be treated as Rooted"),
          _Option(["-progress","progress"], ["input"], None, 0,
                  "Print indications of the progress of run"),
          _Option(["-noprintsets","noprintsets"], ["input"], None, 0,
                  "Do not print out the sets of species"),
          _Option(["-printsets","printsets"], ["input"], None, 0,
                  "Print out the sets of species")]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class ESeqBootCommandline(_EmbossCommandLine):
    """Commandline object for the eseqboot program from EMBOSS (OBSOLETE).

    This is an EMBOSS wrapper around seqboot from PHYLIP.

    It has been replaced by "fseqboot", see FSeqBootCommandline.
    """
    def __init__(self, cmd="eseqboot", **kwargs):
        self.parameters = \
         [_Option(["-datafile","datafile"], ["input", "file"], None, 1,
                  "Input file"),
          _Option(["-randseed","randseed"], ["input"], None, 1,
                  "Random number seed (must be odd)"),
          _Option(["-method","method"], ["input"], None, 1,
                  "Choose the method"),
          _Option(["-test","test"], ["input"], None, 1,
                  "Choose test"),
          _Option(["-reps","reps"], ["input"], None, 1,
                  "How many replicates"),
          _Option(["-inter","inter"], ["input"], None, 0,
                  "Interleaved input"),
          _Option(["-enzymes","enzymes"], ["input"], None, 0,
                  "Present in input file"),
          _Option(["-all","all"], ["input"], None, 0,
                  "All alleles present at each locus"),
          _Option(["-printdata","printdata"], ["input"], None, 0,
                  "Print out the data at start of run"),
          _Option(["-progress","progress"], ["input"], None, 0,
                  "Print indications of progress of run")]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class FDNADistCommandline(_EmbossCommandLine):
    """Commandline object for the fdnadist program from EMBOSS.

    fdnadist is an EMBOSS wrapper for the PHYLIP program dnadist for
    calulating distance matrices from DNA sequence files
    """
    def __init__(self, cmd = "fdnadist", **kwargs):
        self.parameters = \
        [_Option(["-sequence", "sequence"], ["input"], None, 1,
                  "seq file to use (phylip)"),
        _Option(["-method", "method"], ["input"], None, 1,
                 "sub. model [f,k,j,l,s]"),
        _Option(["-gamma", "gamma"], ["input"], None, 0,
                 "gamma [g, i,n]"),
        _Option(["-ncategories", "ncategories"], ["input"], None, 0,
                 "number of rate catergories (1-9)"),
        _Option(["-rate", "rate"], ["input"], None, 0,
                 "rate for each category"),
        _Option(["-categories","categories"], ["input"], None, 0,
                 "File of substitution rate categories"),
        _Option(["-weights", "weights"], ["input"], None, 0,
                 "weights file"),
        _Option(["-gammacoefficient", "gammacoefficient"], ["input"], None, 0,
                 "value for gamma (> 0.001)"),
        _Option(["-invarfrac", "invarfrac"], ["input"], None, 0,
                 "proportoin of invariant sites"),
        _Option(["-ttratio", "ttratio"], ["input"], None, 0,
                 "ts/tv ratio"),
        _Option(["-freqsfrom", "freqsfrom"], ["input"], None, 0,
                 "use emprical base freqs"),
        _Option(["-basefreq", "basefreq"], ["input"], None, 0,
                 "specify basefreqs"),
        _Option(["-lower", "lower"], ["input"], None, 0,
                 "lower triangle matrix (y/N)")]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class FTreeDistCommandline(_EmbossCommandLine):
    """Commandline object for the ftreedist program from EMBOSS.

    ftreedist is an EMBOSS wrapper for the PHYLIP program treedist used for
    calulating distance measures between phylogentic trees
    """
    def __init__(self, cmd = "ftreedist", **kwargs):
        self.parameters = \
        [_Option(["-intreefile", "intreefile"], ["input"], None, 1,
                  "tree file to score (phylip)"),
        _Option(["-dtype", "dtype"], ["input"], None, 0,
                 "distance type ([S]ymetric, [b]ranch score)"),
        _Option(["-pairing", "pairing"], ["input"], None, 0,
                 "tree pairing method ([A]djacent pairs, all [p]ossible pairs)"),
        _Option(["-style", "style"], ["input"], None, 0,
                 "output style - [V]erbose, [f]ill, [s]parse"),
        _Option(["-noroot", "noroot"], ["input"], None, 0,
                 "treat trees as rooted [N/y]"),
        _Option(["-outgrno", "outgrno"], ["input"], None, 0,
                 "which taxon to root the trees with (starts from 0)")]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class FNeighborCommandline(_EmbossCommandLine):
    """Commandline object for the fneighbor program from EMBOSS.

    fneighbor is an EMBOSS wrapper for the PHYLIP program neighbor used for
    calulating neighbor-joining or UPGMA trees from distance matrices 
    """
    def __init__(self, cmd = "fneighbor", **kwargs):
        self.parameters = \
        [_Option(["-datafile", "datafile"], ["input"], None, 1,
                  "dist file to use (phylip)"),
        _Option(["-matrixtype", "matrixtype"], ["input"], None, 0,
                 "is martrix [S]quare pr [u]pper or [l]ower"),
        _Option(["-treetype", "treetype"], ["input"], None, 0,
                 "nj or UPGMA tree (n/u)"),
        _Option(["-outgrno","outgrno" ], ["input"], None, 0,
                 "taxon to use as OG"),
        _Option(["-jumble", "jumble"], ["input"], None, 0,
                 "randommise input order (Y/n)"),
        _Option(["-seed", "seed"], ["input"], None, 0,
                 "provide a random seed"),
        _Option(["-trout", "trout"], ["input"], None, 0,
                 "write tree (Y/n)"),
        _Option(["-outtreefile", "outtreefile"], ["input"], None, 0,
                 "filename for output tree"),
        _Option(["-progress", "progress"], ["input"], None, 0,
                 "print progress (Y/n)"),
        _Option(["-treeprint", "treeprint"], ["input"], None, 0,
                 "print tree (Y/n)")]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class FSeqBootCommandline(_EmbossCommandLine):
    """Commandline object for the fseqboot program from EMBOSS.

    fseqboot is an EMBOSS wrapper for the PHYLIP program seqboot used to
    pseudo-sample alignment files
    """
    def __init__(self, cmd = "fseqboot", **kwargs):
        self.parameters = \
        [_Option(["-sequence", "sequence"], ["input"], None, 1,
                  "seq file to sample (phylip)"),
        _Option(["-categories", "catergories"], ["input"], None, 0,
                 "file of input categories"),
        _Option(["-weights", "weights"], ["input"], None, 0,
                 " weights file"),
        _Option(["-test", "test"], ["input"], None, 0,
                 "specify operation, default is bootstrap"),
        _Option(["-regular", "regular"], ["input"], None, 0,
                 "absolute number to resample"),
        _Option(["-fracsample", "fracsample"], ["input"], None, 0,
                 "fraction to resample"),
        _Option(["-rewriteformat", "rewriteformat"], ["input"], None, 0,
                 "output format ([P]hyilp, [n]exus, [x]ml"),
        _Option(["-seqtype", "seqtype"], ["input"], None, 0,
                 "output format ([D]na, [p]rotein, [r]na"),
        _Option(["-blocksize", "blocksize"], ["input"], None, 0,
                 "print progress (Y/n)"),
        _Option(["-reps", "reps"], ["input"], None, 0,
                 "how many replicates, defaults to 100)"),
        _Option(["-justweights", "jusweights"], ["input"], None, 0,
                 "what to write out [D]atasets of just [w]eights"),
        _Option(["-seed", "seed"], ["input"], None, 0,
                 "specify random seed"),
        _Option(["-dotdiff", "dotdiff"], ["input"], None, 0,
                 "Use dot-differencing? [Y/n]"),]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class FDNAParsCommandline(_EmbossCommandLine):
    """Commandline object for the fdnapars program from EMBOSS.

    fdnapars is an EMBOSS version of the PHYLIP program dnapars, for
    estimating trees from DNA sequences using parsiomny. Calling this command
    without providing a value for the option "-intreefile" will invoke
    "interactive mode" (and as a result fail if called with subprocess) if
    "-auto" is not set to true
    """
    def __init__(self, cmd = "fdnapars", **kwargs):
        self.parameters = \
        [_Option(["-sequence", "sequence"], ["input"], None, 1,
                  "seq file to use (phylip)"),
        _Option(["-intreefile", "intreefile"], ["input"], None, 0,
                 "Phylip tree file"),
        _Option(["-weights", "weights"], ["input"], None, 0,
                 "weights file"),
        _Option(["-maxtrees", "maxtrees"], ["input"], None, 0,
                 "max trees to save during run"),
        _Option(["-thorough", "thorough"], ["input"], None, 0,
                 "more thorough search (Y/n)"),
        _Option(["-rearrange", "rearrange"], ["input"], None, 0,
                 "Rearrange on jsut 1 best tree (Y/n)"),
        _Option(["-transversion", "transversion"], ["input"], None, 0,
                 "Use tranversion parsimony (y/N)"),
        _Option(["-njumble", "njumble"], ["input"], None, 0,
                 "number of times to randomise input order (default is 0)"),
        _Option(["-seed", "seed"], ["input"], None, 0,
                 "provde random seed"),
        _Option(["-outgrno", "outgrno"], ["input"], None, 0,
                 "Specify outgroup"),
        _Option(["-thresh", "thresh"], ["input"], None, 0,
                 "Use threshold parsimony (y/N)"),
        _Option(["-threshold", "threshold"], ["input"], None, 0,
                 "Threshold value"),
        _Option(["-trout", "trout"], ["input"], None, 0,
                 "Write trees to file (Y/n)"),
        _Option(["-outtreefile", "outtreefile"], ["input"], None, 0,
                 "filename for output tree"),
        _Option(["-dotdiff", "dotdiff"], ["input"], None, 0,
                 "Use dot-differencing? [Y/n]")]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class FProtParsCommandline(_EmbossCommandLine):
    """Commandline object for the fdnapars program from EMBOSS.

    fprotpars is an EMBOSS version of the PHYLIP program protpars, for
    estimating trees from protein  sequences using parsiomny. Calling this
    command without providing a value for the option "-intreefile" will invoke
    "interactive mode" (and as a result fail if called with subprocess) if
    "-auto" is not set to true
    """
    def __init__(self, cmd = "fprotpars", **kwargs):
        self.parameters = \
        [_Option(["-sequence", "sequence"], ["input"], None, 1,
                  "seq file to use (phylip)"),
        _Option(["-intreefile", "intreefile"], ["input"], None, 0,
                 "Phylip tree file to score"),
        _Option(["-outtreefile", "outtreefile"], ["input"], None, 1,
                 "phylip tree output file"),
        _Option(["-weights", "weights"], ["input"], None, 0,
                 "weights file"),
        _Option(["-whichcode", "whichcode"], ["input"], None, 0,
                 "which genetic code, [U,M,V,F,Y]]"),
        _Option(["-njumble", "njumble"], ["input"], None, 0,
                 "number of times to randomise input order (default is 0)"),
        _Option(["-seed", "seed"], ["input"], None, 0,
                 "provde random seed"),
        _Option(["-outgrno", "outgrno"], ["input"], None, 0,
                 "Specify outgroup"),
        _Option(["-thresh", "thresh"], ["input"], None, 0,
                 "Use threshold parsimony (y/N)"),
        _Option(["-threshold", "threshold"], ["input"], None, 0,
                 "Threshold value"),
        _Option(["-trout", "trout"], ["input"], None, 0,
                 "Write trees to file (Y/n)"),
        _Option(["-dotdiff", "dotdiff"], ["input"], None, 0,
                 "Use dot-differencing? [Y/n]")]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class FProtDistCommandline(_EmbossCommandLine):
    """Commandline object for the fprotdist program from EMBOSS.

    fprotdist is an EMBOSS wrapper for the PHYLIP program protdist used to
    estimate trees from protein sequences using parsimony
    """
    def __init__(self, cmd = "fprotdist", **kwargs):
        self.parameters = \
        [_Option(["-sequence", "sequence"], ["input"], None, 1,
                  "seq file to use (phylip)"),
        _Option(["-ncategories", "ncategories"], ["input"], None, 0,
                 "number of rate catergories (1-9)"),
        _Option(["-rate", "rate"], ["input"], None, 0,
                 "rate for each category"),
        _Option(["-catergories","catergories"], ["input"], None, 0,
                 "file of rates"),
        _Option(["-weights", "weights"], ["input"], None, 0,
                 "weights file"),
        _Option(["-method", "method"], ["input"], None, 0,
                 "sub. model [j,h,d,k,s,c]"),
        _Option(["-gamma", "gamma"], ["input"], None, 0,
                 "gamma [g, i,c]"),
        _Option(["-gammacoefficient", "gammacoefficient"], ["input"], None, 0,
                 "value for gamma (> 0.001)"),
        _Option(["-invarcoefficient", "invarcoefficient"], ["input"], None, 0,
                 "float for variation of substitution rate among sites"),
        _Option(["-aacateg", "aacateg"], ["input"], None, 0,
                 "Choose the category to use [G,C,H]"),
        _Option(["-whichcode", "whichcode"], ["input"], None, 0,
                 "genetic code [c,m,v,f,y]"),
        _Option(["-ease", "ease"], ["input"], None, 0,
                 "Pob change catergory (float between -0 and 1)"),
        _Option(["-ttratio", "ttratio"], ["input"], None, 0,
                 "Transition/transversion ratio (0-1)"),
        _Option(["-basefreq", "basefreq"], ["input"], None, 0,
                 "DNA base frequencies (space seperated list)")]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class FConsenseCommandline(_EmbossCommandLine):
    """Commandline object for the fconsense program from EMBOSS.

    fconsense is an EMBOSS wrapper for the PHYLIP program consense used to
    calculate consensus trees.
    """
    def __init__(self, cmd = "fconsense", **kwargs):
        self.parameters = \
        [_Option(["-intreefile", "intreefile"], ["input"], None, 1,
                  "file with phylip trees to make consensus from"),
        _Option(["-method", "method"], ["input"], None, 0,
                 "consensus method [s, mr, MRE, ml]"),
        _Option(["-mlfrac", "mlfrac"], ["input"], None, 0,
                 "cut-off freq for a branch to appear in consensus (0.5-1.0)"),
        _Option(["-root", "root"], ["input"], None, 0,
                 "treat trees as rooted (YES, no)"),
        _Option(["-outgrno", "outgrno"], ["input"], None, 0,
                 "OTU to use as outgroup (starts from 0)"),
        _Option(["-trout", "trout"], ["input"], None, 0,
                 "treat trees as rooted (YES, no)"),
        _Option(["-outtreefile", "outtreefile"], ["input"], None, 0,
                 "Phylip tree output file (optional)")]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class WaterCommandline(_EmbossCommandLine):
    """Commandline object for the water program from EMBOSS.
    """
    def __init__(self, cmd="water", **kwargs):
        self.parameters = \
         [_Option(["-asequence","asequence"], ["input", "file"], None, 1,
                  "First sequence to align"),
         _Option(["-bsequence","bsequence"], ["input", "file"], None, 1,
                  "Second sequence to align"),
         _Option(["-gapopen","gapopen"], ["input"], None, 1,
                 "Gap open penalty"),
         _Option(["-gapextend","gapextend"], ["input"], None, 1,
                 "Gap extension penalty"),
         _Option(["-datafile","datafile"], ["input", "file"], None, 0,
                 "Matrix file"),
         _Option(["-similarity","similarity"], ["input"], None, 0,
                 "Display percent identity and similarity"),
         _Option(["-snucleotide","snucleotide"], ["input"], None, 0,
                 "Sequences are nucleotide (boolean)"),
         _Option(["-sprotein","sprotein"], ["input"], None, 0,
                 "Sequences are protein (boolean)"),
         _Option(["-aformat","aformat"], ["input"], None, 0,
                 "Display output in a different specified output format")]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class NeedleCommandline(_EmbossCommandLine):
    """Commandline object for the needle program from EMBOSS.
    """
    def __init__(self, cmd="needle", **kwargs):
        self.parameters = \
         [_Option(["-asequence","asequence"], ["input", "file"], None, 1,
                  "First sequence to align"),
         _Option(["-bsequence","bsequence"], ["input", "file"], None, 1,
                  "Second sequence to align"),
         _Option(["-gapopen","gapopen"], ["input"], None, 1,
                 "Gap open penalty"),
         _Option(["-gapextend","gapextend"], ["input"], None, 1,
                 "Gap extension penalty"),
         _Option(["-datafile","datafile"], ["input", "file"], None, 0,
                 "Matrix file"),
         _Option(["-similarity","similarity"], ["input"], None, 0,
                 "Display percent identity and similarity"),
         _Option(["-snucleotide","snucleotide"], ["input"], None, 0,
                 "Sequences are nucleotide (boolean)"),
         _Option(["-sprotein","sprotein"], ["input"], None, 0,
                 "Sequences are protein (boolean)"),
         _Option(["-aformat","aformat"], ["input"], None, 0,
                 "Display output in a different specified output format")]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class FuzznucCommandline(_EmbossCommandLine):
    """Commandline object for the fuzznuc program from EMBOSS.
    """
    def __init__(self, cmd="fuzznuc", **kwargs):
        self.parameters = [
         _Option(["-sequence","sequence"], ["input"], None, 1,
                 "Sequence database USA"),
         _Option(["-pattern","pattern"], ["input"], None, 1,
                 "Search pattern, using standard IUPAC one-letter codes"),
         _Option(["-mismatch","mismatch"], ["input"], None, 1,
                 "Number of mismatches"),
         _Option(["-complement","complement"], ["input"], None, 0,
                 "Search complementary strand"),
         _Option(["-rformat","rformat"], ["input"], None, 0,
                 "Specify the report format to output in.")]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class Est2GenomeCommandline(_EmbossCommandLine):
    """Commandline object for the est2genome program from EMBOSS.
    """
    def __init__(self, cmd="est2genome", **kwargs):
        self.parameters = [
         _Option(["-est","est"], ["input"], None, 1,
                 "EST sequence(s)"),
         _Option(["-genome","genome"], ["input"], None, 1,
                 "Genomic sequence"),
         _Option(["-match","match"], ["input"], None, 0, 
                 "Score for matching two bases"),
         _Option(["-mismatch","mismatch"], ["input"], None, 0,
                 "Cost for mismatching two bases"),
         _Option(["-gappenalty","gappenalty"], ["input"], None, 0,
                 "Cost for deleting a single base in either sequence, " + \
                 "excluding introns"),
         _Option(["-intronpenalty","intronpenalty"], ["input"], None, 0,
                 "Cost for an intron, independent of length."),
         _Option(["-splicepenalty","splicepenalty"], ["input"], None, 0,
                 "Cost for an intron, independent of length " + \
                 "and starting/ending on donor-acceptor sites"),
         _Option(["-minscore","minscore"], ["input"], None, 0,
                 "Exclude alignments with scores below this threshold score."),
         _Option(["-reverse","reverse"], ["input"], None, 0,
                 "Reverse the orientation of the EST sequence"),
         _Option(["-splice","splice"], ["input"], None, 0,
                 "Use donor and acceptor splice sites."),
         _Option(["-mode","mode"], ["input"], None, 0,
                 "This determines the comparion mode. 'both', 'forward' " + \
                 "'reverse'"),
         _Option(["-best","best"], ["input"], None, 0,
                 "You can print out all comparisons instead of just the best"),
         _Option(["-space","space"], ["input"], None, 0,
                 "for linear-space recursion."),
         _Option(["-shuffle","shuffle"], ["input"], None, 0,
                 "Shuffle"),
         _Option(["-seed","seed"], ["input"], None, 0,
                 "Random number seed"),
         _Option(["-align","align"], ["input"], None, 0,
                 "Show the alignment."),
         _Option(["-width","width"], ["input"], None, 0,
                 "Alignment width")
        ]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class ETandemCommandline(_EmbossCommandLine):
    """Commandline object for the etandem program from EMBOSS.
    """
    def __init__(self, cmd="etandem", **kwargs):
        self.parameters = [
         _Option(["-sequence","sequence"], ["input", "file"], None, 1,
                 "Sequence"),
         _Option(["-minrepeat","minrepeat"], ["input"], None, 1,
                 "Minimum repeat size"),
         _Option(["-maxrepeat","maxrepeat"], ["input"], None, 1,
                 "Maximum repeat size"),
         _Option(["-threshold","threshold"], ["input"], None, 0,
                 "Threshold score"),
         _Option(["-mismatch","mismatch"], ["input"], None, 0,
                   "Allow N as a mismatch"),
         _Option(["-uniform","uniform"], ["input"], None, 0,
                   "Allow uniform consensus"),
         _Option(["-rformat","rformat"], ["output"], None, 0,
                 "Output report format")]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class EInvertedCommandline(_EmbossCommandLine):
    """Commandline object for the einverted program from EMBOSS.
    """
    def __init__(self, cmd="einverted", **kwargs):
        self.parameters = [
         _Option(["-sequence","sequence"], ["input", "file"], None, 1,
                 "Sequence"),
         _Option(["-gap","gap"], ["input", "file"], None, 1,
                 "Gap penalty"),
         _Option(["-threshold","threshold"], ["input"], None, 1,
                 "Minimum score threshold"),
         _Option(["-match","match"], ["input"], None, 1,
                 "Match score"),
         _Option(["-mismatch","mismatch"], ["input"], None, 1,
                   "Mismatch score"),
         _Option(["-maxrepeat","maxrepeat"], ["input"], None, 0,
                 "Maximum separation between the start and end of repeat"),
         ]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class PalindromeCommandline(_EmbossCommandLine):
    """Commandline object for the palindrome program from EMBOSS.
    """
    def __init__(self, cmd="palindrome", **kwargs):
        self.parameters = [
         _Option(["-sequence","sequence"], ["input", "file"], None, 1,
                 "Sequence"),
         _Option(["-minpallen","minpallen"], ["input"], None, 1,
                 "Minimum palindrome length"),
         _Option(["-maxpallen","maxpallen"], ["input"], None, 1,
                 "Maximum palindrome length"),
         _Option(["-gaplimit","gaplimit"], ["input"], None, 1,
                 "Maximum gap between repeats"),
         _Option(["-nummismatches","nummismatches"], ["input"], None, 1,
                 "Number of mismatches allowed"),
         _Option(["-overlap","overlap"], ["input"], None, 1,
                 "Report overlapping matches"),
         ]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class TranalignCommandline(_EmbossCommandLine):
    """Commandline object for the tranalign program from EMBOSS.
    """
    def __init__(self, cmd="tranalign", **kwargs):
        self.parameters = [
         _Option(["-asequence","asequence"], ["input", "file"], None, 1,
                 "Nucleotide sequences to be aligned."),
         _Option(["-bsequence","bsequence"], ["input", "file"], None, 1,
                 "Protein sequence alignment"),
         _Option(["-outseq","outseq"], ["output", "file"], None, 1,
                 "Output sequence file."),
         _Option(["-table","table"], ["input"], None, 0,
                 "Code to use")]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class DiffseqCommandline(_EmbossCommandLine):
    """Commandline object for the diffseq program from EMBOSS.
    """
    def __init__(self, cmd="diffseq", **kwargs):
        self.parameters = [
         _Option(["-asequence","asequence"], ["input", "file"], None, 1,
                 "First sequence to compare"),
         _Option(["-bsequence","bsequence"], ["input", "file"], None, 1,
                 "Second sequence to compare"),
         _Option(["-wordsize","wordsize"], ["input"], None, 1,
                 "Word size to use for comparisons (10 default)"),
         _Option(["-aoutfeat","aoutfeat"], ["output", "file"], None, 1,
                "File for output of first sequence's features"),
         _Option(["-boutfeat","boutfeat"], ["output", "file"], None, 1,
                "File for output of second sequence's features"),
         _Option(["-rformat","rformat"], ["output"], None, 0,
                 "Output report file format")
         ]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class IepCommandline(_EmbossCommandLine):
    """Commandline for EMBOSS iep: calculated isoelectric point and charge.
    """
    def __init__(self, cmd="iep", **kwargs):
        self.parameters = [
         _Option(["-sequence","sequence"], ["input", "file"], None, 1,
                "Protein sequence(s) filename"),
         _Option(["-amino","amino"], ["input"], None, 0),
         _Option(["-lysinemodified","lysinemodified"], ["input"], None, 0),
         _Option(["-disulphides","disulphides"], ["input"], None, 0),
         _Option(["-notermini","notermini"], ["input"], None, 0),
         ]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


#seqret uses -outseq, not -outfile, so use the base class:
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
        self.parameters = [
         _Option(["-sequence","sequence"], ["input", "file"], None, 0,
                 "Input sequence(s) filename"),
         _Option(["-outseq","outseq"], ["output", "file"], None, 0,
                 "Output sequence file."),
         _Option(["-sformat","sformat"], ["input"], None, 0,
                 "Input sequence(s) format (e.g. fasta, genbank)"),
         _Option(["-osformat","osformat"], ["input"], None, 0,
                 "Output sequence(s) format (e.g. fasta, genbank)"),
         ]
        _EmbossMinimalCommandLine.__init__(self, cmd, **kwargs)

    def _validate(self) :
        #Check the outfile, filter, or stdout option has been set.
        #We can't simply do this via the required flag for the outfile
        #output - this seems the simplest solution.
        if not (self.outseq or self.filter or self.stdout) :
            raise ValueError("You must either set outfile (output filename), "
                             "or enable filter or stdout (output to stdout).")
        if not (self.sequence or self.filter or self.stdint) :
            raise ValueError("You must either set sequence (input filename), "
                             "or enable filter or stdin (input from stdin).")
        return _EmbossMinimalCommandLine._validate(self)
