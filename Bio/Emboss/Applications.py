# Copyright 2001-2009 Brad Chapman.
# Revisions copyright 2009-2010 by Peter Cock.
# Revisions copyright 2009 by David Winter.
# Revisions copyright 2009-2010 by Leighton Pritchard.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Code to interact with and run various EMBOSS programs.

These classes follow the AbstractCommandline interfaces for running
programs.
"""

from Bio.Application import _Option, _Switch, AbstractCommandline

class _EmbossMinimalCommandLine(AbstractCommandline):
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
           _Switch(["-auto","auto"], 
                   """Turn off prompts.
           
                   Automatic mode disables prompting, so we recommend you set
                   this argument all the time when calling an EMBOSS tool from
                   Biopython.
                   """),
           _Switch(["-stdout","stdout"],
                   "Write standard output."),
           _Switch(["-filter","filter"],
                   "Read standard input, write standard output."),
           _Switch(["-options","options"],
                   """Prompt for standard and additional values.

                   If you are calling an EMBOSS tool from within Biopython,
                   we DO NOT recommend using this option.
                   """),
           _Switch(["-debug","debug"],
                   "Write debug output to program.dbg."),
           _Switch(["-verbose","verbose"],
                   "Report some/full command line options"),
           _Switch(["-help","help"],
                   """Report command line options.

                   More information on associated and general qualifiers can
                   be found with -help -verbose
                   """),
           _Switch(["-warning","warning"],
                   "Report warnings."),
           _Switch(["-error","error"],
                   "Report errors."),
           _Switch(["-die","die"],
                   "Report dying program messages."),
            ]
        try:
            #Insert extra parameters - at the start just in case there
            #are any arguments which must come last:
            self.parameters = extra_parameters + self.parameters
        except AttributeError:
            #Should we raise an error?  The subclass should have set this up!
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
        extra_parameters = [\
           _Option(["-outfile","outfile"], ["file"], None, 0,
                   "Output filename"),
            ]
        try:
            #Insert extra parameters - at the start just in case there
            #are any arguments which must come last:
            self.parameters = extra_parameters + self.parameters
        except AttributeError:
            #Should we raise an error?  The subclass should have set this up!
            self.parameters = extra_parameters
        _EmbossMinimalCommandLine.__init__(self, cmd, **kwargs)

    def _validate(self):
        #Check the outfile, filter, or stdout option has been set.
        #We can't simply do this via the required flag for the outfile
        #output - this seems the simplest solution.
        if not (self.outfile or self.filter or self.stdout):
            raise ValueError("You must either set outfile (output filename), "
                             "or enable filter or stdout (output to stdout).")
        return _EmbossMinimalCommandLine._validate(self)

class Primer3Commandline(_EmbossCommandLine):
    """Commandline object for the Primer3 interface from EMBOSS.

    The precise set of supported arguments depends on your version of EMBOSS.
    This version accepts arguments current at EMBOSS 6.1.0, but in order to
    remain backwards compatible also support the old argument names as well.

    e.g. Using EMBOSS 6.1.0 or later,
    
    >>> cline = Primer3Commandline(sequence="mysequence.fas", auto=True, hybridprobe=True)
    >>> cline.explainflag = True
    >>> cline.osizeopt=20
    >>> cline.psizeopt=200
    >>> cline.outfile = "myresults.out"
    >>> cline.bogusparameter = 1967     # Invalid parameter
    Traceback (most recent call last):
        ...
    ValueError: Option name bogusparameter was not found.
    >>> print cline
    eprimer3 -auto -outfile=myresults.out -sequence=mysequence.fas -hybridprobe=True -psizeopt=200 -osizeopt=20 -explainflag=True

    The equivalent for anyone still using an older version of EMBOSS would be:

    >>> cline = Primer3Commandline(sequence="mysequence.fas", auto=True, hybridprobe=True)
    >>> cline.explainflag = True
    >>> cline.oligosize=20              # Old EMBOSS, instead of osizeopt
    >>> cline.productosize=200          # Old EMBOSS, instead of psizeopt
    >>> cline.outfile = "myresults.out"
    >>> print cline
    eprimer3 -auto -outfile=myresults.out -sequence=mysequence.fas -hybridprobe=True -productosize=200 -oligosize=20 -explainflag=True
    
    """
    def __init__(self, cmd="eprimer3", **kwargs):
        self.parameters = \
          [_Option(["-sequence","sequence"], [], None, 1,
                   "Sequence to choose primers from."),
           _Option(["-task","task"], [], None, 0,
                   "Tell eprimer3 what task to perform."),
           _Option(["-hybridprobe","hybridprobe"], [], None, 0,
                   "Find an internal oligo to use as a hyb probe."),
           _Option(["-numreturn","numreturn"], [], None, 0,
                   "Maximum number of primer pairs to return."),
           _Option(["-includedregion","includedregion"], [], None, 0,
                   "Subregion of the sequence in which to pick primers."),
           _Option(["-target","target"], [], None, 0,
                   "Sequence to target for flanking primers."),
           _Option(["-excludedregion","excludedregion"], [], None, 0,
                   "Regions to exclude from primer picking."),
           _Option(["-forwardinput","forwardinput"], [], None, 0,
                   "Sequence of a forward primer to check."),
           _Option(["-reverseinput","reverseinput"], [], None, 0,
                   "Sequence of a reverse primer to check."),
           _Option(["-gcclamp","gcclamp"], [], None, 0,
                   "The required number of Gs and Cs at the 3' of each primer."),
           _Option(["-osize","osize"], [], None, 0,
                   "Optimum length of a primer oligo."),
           _Option(["-minsize","minsize"], [], None, 0,
                   "Minimum length of a primer oligo."),
           _Option(["-maxsize","maxsize"], [], None, 0,
                   "Maximum length of a primer oligo."),
           _Option(["-otm","otm"], [], None, 0,
                   "Optimum melting temperature for a primer oligo."),
           _Option(["-mintm","mintm"], [], None, 0,
                   "Minimum melting temperature for a primer oligo."),
           _Option(["-maxtm","maxtm"], [], None, 0,
                   "Maximum melting temperature for a primer oligo."),
           _Option(["-maxdifftm","maxdifftm"], [], None, 0,
                   "Maximum difference in melting temperatures between forward and " +\
                   "reverse primers."),
           _Option(["-ogcpercent","ogcpercent"], [], None, 0,
                   "Optimum GC% for a primer."),
           _Option(["-mingc","mingc"], [], None, 0,
                   "Minimum GC% for a primer."),
           _Option(["-maxgc","maxgc"], [], None, 0,
                   "Maximum GC% for a primer."),
           _Option(["-saltconc","saltconc"], [], None, 0,
                   "Millimolar salt concentration in the PCR."),
           _Option(["-dnaconc","dnaconc"], [], None, 0,
                   "Nanomolar concentration of annealing oligos in the PCR."),
           _Option(["-maxployx","maxployx"], [], None, 0,
                   "Maximum allowable mononucleotide repeat length in a primer."),
           #Primer length:
           _Option(["-productosize","productosize"], [], None, 0,
                   """Optimum size for the PCR product (OBSOLETE).

                   Option replaced in EMBOSS 6.1.0 by -psizeopt
                   """),
           _Option(["-psizeopt", "psizeopt"], [], None, 0, 
                   """Optimum size for the PCR product.

                   Option added in EMBOSS 6.1.0, replacing -productosize
                   """),
           _Option(["-productsizerange","productsizerange"], [], None, 0,
                   """Acceptable range of length for the PCR product (OBSOLETE).

                   Option replaced in EMBOSS 6.1.0 by -prange
                   """),
           _Option(["-prange", "prange"], [], None, 0, 
                   """Acceptable range of length for the PCR product.

                   Option added in EMBOSS 6.1.0, replacing -productsizerange
                   """),
           #Primer temperature:
           _Option(["-productotm","productotm"], [], None, 0,
                   """Optimum melting temperature for the PCR product (OBSOLETE).

                   Option replaced in EMBOSS 6.1.0 by -ptmopt
                   """),
           _Option(["-ptmopt", "ptmopt"], [], None, 0, 
                   """Optimum melting temperature for the PCR product.

                   Option added in EMBOSS 6.1.0, replacing -productotm
                   """),
           _Option(["-productmintm","productmintm"], [], None, 0,
                   """Minimum allowed melting temperature for the amplicon (OBSOLETE)

                   Option replaced in EMBOSS 6.1.0 by -ptmmin
                   """),
           _Option(["-ptmmin", "ptmmin"], [], None, 0, 
                   """Minimum allowed melting temperature for the amplicon."),

                   Option added in EMBOSS 6.1.0, replacing -productmintm
                   """),
           _Option(["-productmaxtm","productmaxtm"], [], None, 0,
                   """Maximum allowed melting temperature for the amplicon (OBSOLETE).

                   Option replaced in EMBOSS 6.1.0 by -ptmmax
                   """),
           _Option(["-ptmmax", "ptmmax"], [], None, 0, 
                   """Maximum allowed melting temperature for the amplicon."),

                   Option added in EMBOSS 6.1.0, replacing -productmaxtm
                   """),
           #Note to self, should be -oexcludedregion not -oexcluderegion
           _Option(["-oexcludedregion", "oexcludedregion"], [], None, 0, 
                   """Do not pick internal oligos in this region."),

                   Option added in EMBOSS 6.1.0, replacing -oligoexcludedregion.
                   """),
           _Option(["-oligoexcludedregion", "oligoexcludedregion"], [],
                   None, 0,
                   """Do not pick internal oligos in this region (OBSOLETE)."),

                   Option replaced in EMBOSS 6.1.0 by -oexcluderegion.
                   """),
           _Option(["-oligoinput","oligoinput"], [], None, 0,
                   "Sequence of the internal oligo."),
           #Oligo length:
           _Option(["-oligosize","oligosize"], [], None, 0, 
                   """Optimum length of internal oligo (OBSOLETE).
 
                   Option replaced in EMBOSS 6.1.0 by -osizeopt.
                   """),
           _Option(["-osizeopt", "osizeopt"], [], None, 0, 
                   """Optimum length of internal oligo.

                   Option added in EMBOSS 6.1.0, replaces -oligosize
                   """),
           _Option(["-oligominsize","oligominsize"], [], None, 0, 
                   """Minimum length of internal oligo (OBSOLETE)."),
 
                   Option replaced in EMBOSS 6.1.0 by -ominsize.
                   """),
           _Option(["-ominsize", "ominsize"], [], None, 0, 
                   """Minimum length of internal oligo."

                   Option added in EMBOSS 6.1.0, replaces -oligominsize
                   """),
           _Option(["-oligomaxsize","oligomaxsize"], [], None, 0, 
                   """Maximum length of internal oligo (OBSOLETE).
 
                   Option replaced in EMBOSS 6.1.0 by -omaxsize.
                   """),
           _Option(["-omaxsize", "omaxsize"], [], None, 0, 
                   """Maximum length of internal oligo.

                   Option added in EMBOSS 6.1.0, replaces -oligomaxsize
                   """),
           #Oligo GC temperature:
           _Option(["-oligotm","oligotm"], [], None, 0, 
                   """Optimum melting temperature of internal oligo (OBSOLETE).
 
                   Option replaced in EMBOSS 6.1.0 by -otmopt.
                   """),
           _Option(["-otmopt", "otmopt"], [], None, 0, 
                   """Optimum melting temperature of internal oligo.

                   Option added in EMBOSS 6.1.0.
                   """),
           _Option(["-oligomintm","oligomintm"], [], None, 0, 
                   """Minimum melting temperature of internal oligo (OBSOLETE).
 
                   Option replaced in EMBOSS 6.1.0 by -otmmin.
                   """),
           _Option(["-otmmin", "otmmin"], [], None, 0, 
                   """Minimum melting temperature of internal oligo.

                   Option added in EMBOSS 6.1.0, replacing -oligomintm
                   """),
           _Option(["-oligomaxtm","oligomaxtm"], [], None, 0, 
                   """Maximum melting temperature of internal oligo (OBSOLETE).
 
                   Option replaced in EMBOSS 6.1.0 by -otmmax.
                   """),
           _Option(["-otmmax", "otmmax"], [], None, 0, 
                   """Maximum melting temperature of internal oligo.

                   Option added in EMBOSS 6.1.0, replacing -oligomaxtm
                   """),
           #Oligo GC percent:
           _Option(["-oligoogcpercent","oligoogcpercent"], [], None, 0, 
                   """Optimum GC% for internal oligo (OBSOLETE).
 
                   Option replaced in EMBOSS 6.1.0 by -ogcopt.
                   """),
           _Option(["-ogcopt", "ogcopt"], [], None, 0, 
                   """Optimum GC% for internal oligo."

                   Option added in EMBOSS 6.1.0, replacing -oligoogcpercent
                   """),
           _Option(["-oligomingc","oligomingc"], [], None, 0, 
                   """Minimum GC% for internal oligo (OBSOLETE).
 
                   Option replaced in EMBOSS 6.1.0 by -ogcmin.
                   """),
           _Option(["-ogcmin", "ogcmin"], [], None, 0, 
                   """Minimum GC% for internal oligo.

                   Option added in EMBOSS 6.1.0, replacing -oligomingc
                   """),
           _Option(["-oligomaxgc","oligomaxgc"], [], None, 0, 
                   """Maximum GC% for internal oligo.
 
                   Option replaced in EMBOSS 6.1.0 by -ogcmax
                   """),
           _Option(["-ogcmax", "ogcmax"], [], None, 0, 
                   """Maximum GC% for internal oligo."),

                   Option added in EMBOSS 6.1.0, replacing -oligomaxgc
                   """),
           #Oligo salt concentration:
           _Option(["-oligosaltconc","oligosaltconc"], [], None, 0, 
                   """Millimolar concentration of salt in the hybridisation."),
 
                   Option replaced in EMBOSS 6.1.0 by -osaltconc
                   """),
           _Option(["-osaltconc", "osaltconc"], [], None, 0, 
                   """Millimolar concentration of salt in the hybridisation."),

                   Option added in EMBOSS 6.1.0, replacing -oligosaltconc
                   """),
           _Option(["-oligodnaconc","oligodnaconc"], [], None, 0, 
                   """Nanomolar concentration of internal oligo in the hybridisation.
 
                   Option replaced in EMBOSS 6.1.0 by -odnaconc
                   """),
           _Option(["-odnaconc", "odnaconc"], [], None, 0, 
                   """Nanomolar concentration of internal oligo in the hybridisation.

                   Option added in EMBOSS 6.1.0, replacing -oligodnaconc
                   """),
           #Oligo self complementarity
           _Option(["-oligoselfany","oligoselfany"], [], None, 0, 
                   """Maximum allowable alignment score for self-complementarity (OBSOLETE).
 
                   Option replaced in EMBOSS 6.1.0 by -oanyself
                   """),
           _Option(["-oanyself", "oanyself"], [], None, 0, 
                   """Maximum allowable alignment score for self-complementarity."),

                   Option added in EMBOSS 6.1.0, replacing -oligoselfany
                   """),
           _Option(["-oligoselfend","oligoselfend"], [], None, 0, 
                   """Maximum allowable 3`-anchored global alignment score for " +\
                   self-complementarity (OBSOLETE).
 
                   Option replaced in EMBOSS 6.1.0 by -oendself
                   """),
           _Option(["-oendself", "oendself"], [], None, 0, 
                   """Max 3`-anchored self-complementarity global alignment score.

                   Option added in EMBOSS 6.1.0, replacing -oligoselfend
                   """),
           _Option(["-oligomaxpolyx","oligomaxpolyx"], [], None, 0, 
                   """Maximum length of mononucleotide repeat in internal oligo (OBSOLETE).

                   Option replaced in EMBOSS 6.1.0 by -opolyxmax
                   """),
           _Option(["-opolyxmax", "opolyxmax"], [], None, 0, 
                   """Maximum length of mononucleotide repeat in internal oligo."),

                   Option added in EMBOSS 6.1.0, replacing -oligomaxpolyx
                   """),
           _Option(["-mispriminglibraryfile","mispriminglibraryfile"], [], None, 0,
                    "File containing library of sequences to avoid amplifying"),
           _Option(["-maxmispriming","maxmispriming"], [], None, 0,
                   "Maximum allowed similarity of primers to sequences in " +\
                   "library specified by -mispriminglibrary"),
           _Option(["-oligomaxmishyb","oligomaxmishyb"], [], None, 0, 
                   """Maximum alignment score for hybridisation of internal oligo to
                   library specified by -oligomishyblibraryfile (OBSOLETE).

                   Option replaced in EMBOSS 6.1.0 by -omishybmax
                   """),
           _Option(["-omishybmax", "omishybmax"], [], None, 0, 
                   """Maximum alignment score for hybridisation of internal oligo to
                   library specified by -mishyblibraryfile.

                   Option added in EMBOSS 6.1.0, replacing -oligomaxmishyb
                   """),
           _Option(["-oligomishyblibraryfile", "oligomishyblibraryfile"],
                   [], None, 0,
                    """Library file of seqs to avoid internal oligo hybridisation (OBSOLETE).

                   Option replaced in EMBOSS 6.1.0 by -mishyblibraryfile
                   """),
           _Option(["-mishyblibraryfile", "mishyblibraryfile"], [], None, 0,
                    """Library file of seqs to avoid internal oligo hybridisation.

                   Option added in EMBOSS 6.1.0, replacing -oligomishyblibraryfile
                   """),
           _Option(["-explainflag","explainflag"], [], None, 0,
                   "Produce output tags with eprimer3 statistics"),
           ]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)
        

class PrimerSearchCommandline(_EmbossCommandLine):
    """Commandline object for the primersearch program from EMBOSS.
    """
    def __init__(self, cmd="primersearch", **kwargs):
        self.parameters = \
         [_Option(["-seqall","-sequences","sequences","seqall"], [],
                  None, 1, "Sequence to look for the primer pairs in."),
                  #When this wrapper was written primersearch used -sequences
                  #as the argument name. Since at least EMBOSS 5.0 (and
                  #perhaps earlier) this has been -seqall instead.
          _Option(["-infile","-primers","primers","infile"], ["file"],
                  None, 1, "File containing the primer pairs to search for."),
                  #When this wrapper was written primersearch used -primers
                  #as the argument name. Since at least EMBOSS 5.0 (and
                  #perhaps earlier) this has been -infile instead.
          _Option(["-mismatchpercent","mismatchpercent"], [], None, 1,
                  "Allowed percentage mismatch (any integer value, default 0)."),
          _Option(["-snucleotide","snucleotide"], [], None, 0,
                  "Sequences are nucleotide (boolean)"),
          _Option(["-sprotein","sprotein"], [], None, 0,
                  "Sequences are protein (boolean)"),
          ]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class EProtDistCommandline(_EmbossCommandLine):
    """Commandline object for the eprotdist program from EMBOSS (DEPRECATED).

    This is an EMBOSS wrapper around protdist from PHYLIP.

    It has been replaced by "fprotdist", see FProtDistCommandline.
    """
    def __init__(self, cmd="eprotdist", **kwargs):
        import warnings
        import Bio
        warnings.warn("Bio.Emboss.Application.EProtDistCommandline has been deprecated; please use 'fprotdist' instead (see FProtDistCommandline).", Bio.BiopythonDeprecationWarning)
        self.parameters = \
         [_Option(["-msf","msf"], [], None, 1,
                  "File containing sequences"),
          _Option(["-method","method"], [], None, 1,
                  "Choose the method to use"),
          _Option(["-categ","categ"], [], None, 0,
                  "Choose the category to use"),
          _Option(["-gencode","gencode"], [], None, 0,
                  "Which genetic code"),
          _Option(["-prob","prob"], [], None, 0,
                  "Prob change category (1.0=easy)"),
          _Option(["-tranrate","tranrate"], [], None, 0,
                  "Transition/transversion ratio"),
          _Option(["-freqa","freqa"], [], None, 0,
                  "Frequency for A"),
          _Option(["-freqc","freqc"], [], None, 0,
                  "Frequency for C"),
          _Option(["-freqg","freqg"], [], None, 0,
                  "Frequency for G"),
          _Option(["-freqt","freqt"], [], None, 0,
                  "Frequency for T"),
          _Option(["-printdata","printdata"], [], None, 0,
                  "Print out the data at start of run"),
          _Option(["-progress","progress"], [], None, 0,
                  "Print indications of progress of run"),
          _Option(["-basefrequency","basefrequency"], [], None, 0,
                  "Use empirical base frequencies")]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class ENeighborCommandline(_EmbossCommandLine):
    """Commandline object for the eneighbor program from EMBOSS (DEPRECATED).

    This is an EMBOSS wrapper around neighbor from PHYLIP.

    It has been replaced by "fneighbor", see FNeighborCommandline.
    """
    def __init__(self, cmd="eneighbor", **kwargs):
        import warnings
        import Bio
        warnings.warn("Bio.Emboss.Application.ENeighborCommandline has been deprecated; please use 'fneighbor' instead (see FNeighborCommandline).", Bio.BiopythonDeprecationWarning)
        self.parameters = \
         [_Option(["-infile","infile"], [], None, 1,
                  "infile value"),
          _Option(["-trout","trout"], [], None, 1,
                  "Create a tree file"),
          _Option(["-treefile","treefile"], [], None, 1,
                  "Tree file name"),
          _Option(["-nj","nj"], [], None, 1,
                  "Neighbor-joining"),
          _Option(["-noog","noog"], [], None, 1,
                  "Outgroup root"),
          _Option(["-outgnum","outgnum"], [], None, 0,
                  "number of the outgroup"),
          _Option(["-randseed","randseed"], [], None, 0,
                  "Random number seed (must be odd)"),
          _Option(["-datasets","datasets"], [], None, 0,
                  "How many data sets"),
          _Option(["-drawtree","drawtree"], [], None, 0,
                  "Draw tree"),
          _Option(["-lt","lt"], [], None, 0,
                  "Lower-triangular data matrix"),
          _Option(["-ut","ut"], [], None, 0,
                  "Upper-triangular data matrix"),
          _Option(["-sr","sr"], [], None, 0,
                  "Subreplicates"),
          _Option(["-random","random"], [], None, 0,
                  "Randomize input order of species"),
          _Option(["-multsets","multsets"], [], None, 0,
                  "Analyze multiple data sets"),
          _Option(["-printdata","printdata"], [], None, 0,
                  "Print out the data at start of run"),
          _Option(["-progress","progress"], [], None, 0,
                  "Print indications of progress of run")]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class EProtParsCommandline(_EmbossCommandLine):
    """Commandline object for the eprotpars program from EMBOSS (DEPRECATED).

    This is an EMBOSS wrapper around protpars from PHYLIP.

    It has been replaced by "fprotpars", see FProtParsCommandline.
    """
    def __init__(self, cmd="eprotpars", **kwargs):
        import warnings
        import Bio
        warnings.warn("Bio.Emboss.Application.EProtParsCommandline has been deprecated; please use 'fprotpars' instead (see FProtParsCommandline).", Bio.BiopythonDeprecationWarning)
        self.parameters = \
         [_Option(["-msf","msf"], ["file"], None, 1,
                  "Sequences file to be read in"),
          _Option(["-besttree","besttree"], [], None, 0,
                  "Search for the best tree"),
          _Option(["-random","random"], [], None, 0,
                  "Randomize input order of species"),
          _Option(["-norandom","norandom"], [], None, 0,
                  "Do not randomize input order of species"),
          _Option(["-randseed","randseed"], [], None, 0,
                  "Random number seed (must be odd)"),
          _Option(["-randtimes","randtimes"], [], None, 0,
                  "How many times to randomize"),
          _Option(["-og","og"], [], None, 0,
                  "Use an outgroup root"),
          _Option(["-noog","noog"], [], None, 0,
                  "Do not use an outgroup root"),
          _Option(["-outgnum","outgnum"], [], None, 0,
                  "Number of the outgroup"),
          _Option(["-thresh","thresh"], [], None, 0,
                  "Use Threshold parsimony"),
          _Option(["-valthresh","valthresh"], [], None, 0,
                  "threshold value"),
          _Option(["-printdata","printdata"], [], None, 0,
                  "Print out the data at start of run"),
          _Option(["-progress","progress"], [], None, 0,
                  "Print indications of progress of run"),
          _Option(["-steps","steps"], [], None, 0,
                  "Print out steps in each site"),
          _Option(["-seqatnodes","seqatnodes"], [], None, 0,
                  "Print sequences at all nodes of tree"),
          _Option(["-drawtree","drawtree"], [], None, 0,
                  "Draw tree"),
          _Option(["-trout","trout"], [], None, 0,
                  "Create a tree file"),
          _Option(["-notrout","notrout"], [], None, 0,
                  "Do not create a tree file"),
          _Option(["-treefile","treefile"], ["file"], None, 0,
                  "Output treefile name")]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class EConsenseCommandline(_EmbossCommandLine):
    """Commandline object for the econsense program from EMBOSS (DEPRECATED).

    This is an EMBOSS wrapper around consense from PHYLIP.

    It has been replaced by "fconsense", see FConsenseCommandline.
    """
    def __init__(self, cmd="econsense", **kwargs):
        import warnings
        import Bio
        warnings.warn("Bio.Emboss.Application.EConsenseCommandline has been deprecated; please use 'fconsense' instead (see FConsenseCommandline).", Bio.BiopythonDeprecationWarning)
        self.parameters = \
         [_Option(["-infile","infile"], ["file"], None, 1,
                  "file to read in (New Hampshire standard form)"),
          _Option(["-notrout","notrout"], [], None, 0,
                  "Do not create a tree file"),
          _Option(["-trout","trout"], [], None, 0,
                  "Create a tree file"),
          _Option(["-treefile","treefile"], ["file"], None, 0,
                  "tree file name"),
          _Option(["-noog","noog"], [], None, 0,
                  "Do not use an outgroup"),
          _Option(["-og","og"], [], None, 0,
                  "Use an outgroup"),
          _Option(["-outgnum","outgnum"], [], None, 0,
                  "number of the outgroup"),
          _Option(["-nodrawtree","nodrawtree"], [], None, 0,
                  "Do not draw a tree"),
          _Option(["-drawtree","drawtree"], [], None, 0,
                  "Draw tree"),
          _Option(["-root","root"], [], None, 0,
                  "Trees to be treated as Rooted"),
          _Option(["-progress","progress"], [], None, 0,
                  "Print indications of the progress of run"),
          _Option(["-noprintsets","noprintsets"], [], None, 0,
                  "Do not print out the sets of species"),
          _Option(["-printsets","printsets"], [], None, 0,
                  "Print out the sets of species")]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class ESeqBootCommandline(_EmbossCommandLine):
    """Commandline object for the eseqboot program from EMBOSS (DEPRECATED).

    This is an EMBOSS wrapper around seqboot from PHYLIP.

    It has been replaced by "fseqboot", see FSeqBootCommandline.
    """
    def __init__(self, cmd="eseqboot", **kwargs):
        import warnings
        import Bio
        warnings.warn("Bio.Emboss.Application.ESeqBootCommandline has been deprecated; please use 'fseqboot' instead (see FSeqBootCommandline).", Bio.BiopythonDeprecationWarning)
        self.parameters = \
         [_Option(["-datafile","datafile"], ["file"], None, 1,
                  "Input file"),
          _Option(["-randseed","randseed"], [], None, 1,
                  "Random number seed (must be odd)"),
          _Option(["-method","method"], [], None, 1,
                  "Choose the method"),
          _Option(["-test","test"], [], None, 1,
                  "Choose test"),
          _Option(["-reps","reps"], [], None, 1,
                  "How many replicates"),
          _Option(["-inter","inter"], [], None, 0,
                  "Interleaved input"),
          _Option(["-enzymes","enzymes"], [], None, 0,
                  "Present in input file"),
          _Option(["-all","all"], [], None, 0,
                  "All alleles present at each locus"),
          _Option(["-printdata","printdata"], [], None, 0,
                  "Print out the data at start of run"),
          _Option(["-progress","progress"], [], None, 0,
                  "Print indications of progress of run")]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class FDNADistCommandline(_EmbossCommandLine):
    """Commandline object for the fdnadist program from EMBOSS.

    fdnadist is an EMBOSS wrapper for the PHYLIP program dnadist for
    calulating distance matrices from DNA sequence files.
    """
    def __init__(self, cmd = "fdnadist", **kwargs):
        self.parameters = \
        [_Option(["-sequence", "sequence"], [], None, 1,
                  "seq file to use (phylip)"),
        _Option(["-method", "method"], [], None, 1,
                 "sub. model [f,k,j,l,s]"),
        _Option(["-gamma", "gamma"], [], None, 0,
                 "gamma [g, i,n]"),
        _Option(["-ncategories", "ncategories"], [], None, 0,
                 "number of rate catergories (1-9)"),
        _Option(["-rate", "rate"], [], None, 0,
                 "rate for each category"),
        _Option(["-categories","categories"], [], None, 0,
                 "File of substitution rate categories"),
        _Option(["-weights", "weights"], [], None, 0,
                 "weights file"),
        _Option(["-gammacoefficient", "gammacoefficient"], [], None, 0,
                 "value for gamma (> 0.001)"),
        _Option(["-invarfrac", "invarfrac"], [], None, 0,
                 "proportoin of invariant sites"),
        _Option(["-ttratio", "ttratio"], [], None, 0,
                 "ts/tv ratio"),
        _Option(["-freqsfrom", "freqsfrom"], [], None, 0,
                 "use emprical base freqs"),
        _Option(["-basefreq", "basefreq"], [], None, 0,
                 "specify basefreqs"),
        _Option(["-lower", "lower"], [], None, 0,
                 "lower triangle matrix (y/N)")]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class FTreeDistCommandline(_EmbossCommandLine):
    """Commandline object for the ftreedist program from EMBOSS.

    ftreedist is an EMBOSS wrapper for the PHYLIP program treedist used for
    calulating distance measures between phylogentic trees.
    """
    def __init__(self, cmd = "ftreedist", **kwargs):
        self.parameters = \
        [_Option(["-intreefile", "intreefile"], [], None, 1,
                  "tree file to score (phylip)"),
        _Option(["-dtype", "dtype"], [], None, 0,
                 "distance type ([S]ymetric, [b]ranch score)"),
        _Option(["-pairing", "pairing"], [], None, 0,
                 "tree pairing method ([A]djacent pairs, all [p]ossible pairs)"),
        _Option(["-style", "style"], [], None, 0,
                 "output style - [V]erbose, [f]ill, [s]parse"),
        _Option(["-noroot", "noroot"], [], None, 0,
                 "treat trees as rooted [N/y]"),
        _Option(["-outgrno", "outgrno"], [], None, 0,
                 "which taxon to root the trees with (starts from 0)")]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class FNeighborCommandline(_EmbossCommandLine):
    """Commandline object for the fneighbor program from EMBOSS.

    fneighbor is an EMBOSS wrapper for the PHYLIP program neighbor used for
    calulating neighbor-joining or UPGMA trees from distance matrices.
    """
    def __init__(self, cmd = "fneighbor", **kwargs):
        self.parameters = \
        [_Option(["-datafile", "datafile"], [], None, 1,
                  "dist file to use (phylip)"),
        _Option(["-matrixtype", "matrixtype"], [], None, 0,
                 "is martrix [S]quare pr [u]pper or [l]ower"),
        _Option(["-treetype", "treetype"], [], None, 0,
                 "nj or UPGMA tree (n/u)"),
        _Option(["-outgrno","outgrno" ], [], None, 0,
                 "taxon to use as OG"),
        _Option(["-jumble", "jumble"], [], None, 0,
                 "randommise input order (Y/n)"),
        _Option(["-seed", "seed"], [], None, 0,
                 "provide a random seed"),
        _Option(["-trout", "trout"], [], None, 0,
                 "write tree (Y/n)"),
        _Option(["-outtreefile", "outtreefile"], [], None, 0,
                 "filename for output tree"),
        _Option(["-progress", "progress"], [], None, 0,
                 "print progress (Y/n)"),
        _Option(["-treeprint", "treeprint"], [], None, 0,
                 "print tree (Y/n)")]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class FSeqBootCommandline(_EmbossCommandLine):
    """Commandline object for the fseqboot program from EMBOSS.

    fseqboot is an EMBOSS wrapper for the PHYLIP program seqboot used to
    pseudo-sample alignment files.
    """
    def __init__(self, cmd = "fseqboot", **kwargs):
        self.parameters = \
        [_Option(["-sequence", "sequence"], [], None, 1,
                  "seq file to sample (phylip)"),
        _Option(["-categories", "catergories"], [], None, 0,
                 "file of input categories"),
        _Option(["-weights", "weights"], [], None, 0,
                 " weights file"),
        _Option(["-test", "test"], [], None, 0,
                 "specify operation, default is bootstrap"),
        _Option(["-regular", "regular"], [], None, 0,
                 "absolute number to resample"),
        _Option(["-fracsample", "fracsample"], [], None, 0,
                 "fraction to resample"),
        _Option(["-rewriteformat", "rewriteformat"], [], None, 0,
                 "output format ([P]hyilp, [n]exus, [x]ml"),
        _Option(["-seqtype", "seqtype"], [], None, 0,
                 "output format ([D]na, [p]rotein, [r]na"),
        _Option(["-blocksize", "blocksize"], [], None, 0,
                 "print progress (Y/n)"),
        _Option(["-reps", "reps"], [], None, 0,
                 "how many replicates, defaults to 100)"),
        _Option(["-justweights", "jusweights"], [], None, 0,
                 "what to write out [D]atasets of just [w]eights"),
        _Option(["-seed", "seed"], [], None, 0,
                 "specify random seed"),
        _Option(["-dotdiff", "dotdiff"], [], None, 0,
                 "Use dot-differencing? [Y/n]"),]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class FDNAParsCommandline(_EmbossCommandLine):
    """Commandline object for the fdnapars program from EMBOSS.

    fdnapars is an EMBOSS version of the PHYLIP program dnapars, for
    estimating trees from DNA sequences using parsiomny. Calling this command
    without providing a value for the option "-intreefile" will invoke
    "interactive mode" (and as a result fail if called with subprocess) if
    "-auto" is not set to true.
    """
    def __init__(self, cmd = "fdnapars", **kwargs):
        self.parameters = \
        [_Option(["-sequence", "sequence"], [], None, 1,
                  "seq file to use (phylip)"),
        _Option(["-intreefile", "intreefile"], [], None, 0,
                 "Phylip tree file"),
        _Option(["-weights", "weights"], [], None, 0,
                 "weights file"),
        _Option(["-maxtrees", "maxtrees"], [], None, 0,
                 "max trees to save during run"),
        _Option(["-thorough", "thorough"], [], None, 0,
                 "more thorough search (Y/n)"),
        _Option(["-rearrange", "rearrange"], [], None, 0,
                 "Rearrange on jsut 1 best tree (Y/n)"),
        _Option(["-transversion", "transversion"], [], None, 0,
                 "Use tranversion parsimony (y/N)"),
        _Option(["-njumble", "njumble"], [], None, 0,
                 "number of times to randomise input order (default is 0)"),
        _Option(["-seed", "seed"], [], None, 0,
                 "provde random seed"),
        _Option(["-outgrno", "outgrno"], [], None, 0,
                 "Specify outgroup"),
        _Option(["-thresh", "thresh"], [], None, 0,
                 "Use threshold parsimony (y/N)"),
        _Option(["-threshold", "threshold"], [], None, 0,
                 "Threshold value"),
        _Option(["-trout", "trout"], [], None, 0,
                 "Write trees to file (Y/n)"),
        _Option(["-outtreefile", "outtreefile"], [], None, 0,
                 "filename for output tree"),
        _Option(["-dotdiff", "dotdiff"], [], None, 0,
                 "Use dot-differencing? [Y/n]")]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class FProtParsCommandline(_EmbossCommandLine):
    """Commandline object for the fdnapars program from EMBOSS.

    fprotpars is an EMBOSS version of the PHYLIP program protpars, for
    estimating trees from protein  sequences using parsiomny. Calling this
    command without providing a value for the option "-intreefile" will invoke
    "interactive mode" (and as a result fail if called with subprocess) if
    "-auto" is not set to true.
    """
    def __init__(self, cmd = "fprotpars", **kwargs):
        self.parameters = \
        [_Option(["-sequence", "sequence"], [], None, 1,
                  "seq file to use (phylip)"),
        _Option(["-intreefile", "intreefile"], [], None, 0,
                 "Phylip tree file to score"),
        _Option(["-outtreefile", "outtreefile"], [], None, 1,
                 "phylip tree output file"),
        _Option(["-weights", "weights"], [], None, 0,
                 "weights file"),
        _Option(["-whichcode", "whichcode"], [], None, 0,
                 "which genetic code, [U,M,V,F,Y]]"),
        _Option(["-njumble", "njumble"], [], None, 0,
                 "number of times to randomise input order (default is 0)"),
        _Option(["-seed", "seed"], [], None, 0,
                 "provde random seed"),
        _Option(["-outgrno", "outgrno"], [], None, 0,
                 "Specify outgroup"),
        _Option(["-thresh", "thresh"], [], None, 0,
                 "Use threshold parsimony (y/N)"),
        _Option(["-threshold", "threshold"], [], None, 0,
                 "Threshold value"),
        _Option(["-trout", "trout"], [], None, 0,
                 "Write trees to file (Y/n)"),
        _Option(["-dotdiff", "dotdiff"], [], None, 0,
                 "Use dot-differencing? [Y/n]")]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class FProtDistCommandline(_EmbossCommandLine):
    """Commandline object for the fprotdist program from EMBOSS.

    fprotdist is an EMBOSS wrapper for the PHYLIP program protdist used to
    estimate trees from protein sequences using parsimony
    """
    def __init__(self, cmd = "fprotdist", **kwargs):
        self.parameters = \
        [_Option(["-sequence", "sequence"], [], None, 1,
                  "seq file to use (phylip)"),
        _Option(["-ncategories", "ncategories"], [], None, 0,
                 "number of rate catergories (1-9)"),
        _Option(["-rate", "rate"], [], None, 0,
                 "rate for each category"),
        _Option(["-catergories","catergories"], [], None, 0,
                 "file of rates"),
        _Option(["-weights", "weights"], [], None, 0,
                 "weights file"),
        _Option(["-method", "method"], [], None, 0,
                 "sub. model [j,h,d,k,s,c]"),
        _Option(["-gamma", "gamma"], [], None, 0,
                 "gamma [g, i,c]"),
        _Option(["-gammacoefficient", "gammacoefficient"], [], None, 0,
                 "value for gamma (> 0.001)"),
        _Option(["-invarcoefficient", "invarcoefficient"], [], None, 0,
                 "float for variation of substitution rate among sites"),
        _Option(["-aacateg", "aacateg"], [], None, 0,
                 "Choose the category to use [G,C,H]"),
        _Option(["-whichcode", "whichcode"], [], None, 0,
                 "genetic code [c,m,v,f,y]"),
        _Option(["-ease", "ease"], [], None, 0,
                 "Pob change catergory (float between -0 and 1)"),
        _Option(["-ttratio", "ttratio"], [], None, 0,
                 "Transition/transversion ratio (0-1)"),
        _Option(["-basefreq", "basefreq"], [], None, 0,
                 "DNA base frequencies (space seperated list)")]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class FConsenseCommandline(_EmbossCommandLine):
    """Commandline object for the fconsense program from EMBOSS.

    fconsense is an EMBOSS wrapper for the PHYLIP program consense used to
    calculate consensus trees.
    """
    def __init__(self, cmd = "fconsense", **kwargs):
        self.parameters = \
        [_Option(["-intreefile", "intreefile"], [], None, 1,
                  "file with phylip trees to make consensus from"),
        _Option(["-method", "method"], [], None, 0,
                 "consensus method [s, mr, MRE, ml]"),
        _Option(["-mlfrac", "mlfrac"], [], None, 0,
                 "cut-off freq for a branch to appear in consensus (0.5-1.0)"),
        _Option(["-root", "root"], [], None, 0,
                 "treat trees as rooted (YES, no)"),
        _Option(["-outgrno", "outgrno"], [], None, 0,
                 "OTU to use as outgroup (starts from 0)"),
        _Option(["-trout", "trout"], [], None, 0,
                 "treat trees as rooted (YES, no)"),
        _Option(["-outtreefile", "outtreefile"], [], None, 0,
                 "Phylip tree output file (optional)")]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class WaterCommandline(_EmbossCommandLine):
    """Commandline object for the water program from EMBOSS.
    """
    def __init__(self, cmd="water", **kwargs):
        self.parameters = \
         [_Option(["-asequence","asequence"], ["file"], None, 1,
                  "First sequence to align"),
         _Option(["-bsequence","bsequence"], ["file"], None, 1,
                  "Second sequence to align"),
         _Option(["-gapopen","gapopen"], [], None, 1,
                 "Gap open penalty"),
         _Option(["-gapextend","gapextend"], [], None, 1,
                 "Gap extension penalty"),
         _Option(["-datafile","datafile"], ["file"], None, 0,
                 "Matrix file"),
         _Option(["-similarity","similarity"], [], None, 0,
                 "Display percent identity and similarity"),
         _Option(["-snucleotide","snucleotide"], [], None, 0,
                 "Sequences are nucleotide (boolean)"),
         _Option(["-sprotein","sprotein"], [], None, 0,
                 "Sequences are protein (boolean)"),
         _Option(["-aformat","aformat"], [], None, 0,
                 "Display output in a different specified output format")]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class NeedleCommandline(_EmbossCommandLine):
    """Commandline object for the needle program from EMBOSS.
    """
    def __init__(self, cmd="needle", **kwargs):
        self.parameters = \
         [_Option(["-asequence","asequence"], ["file"], None, 1,
                  "First sequence to align"),
         _Option(["-bsequence","bsequence"], ["file"], None, 1,
                  "Second sequence to align"),
         _Option(["-gapopen","gapopen"], [], None, 1,
                 "Gap open penalty"),
         _Option(["-gapextend","gapextend"], [], None, 1,
                 "Gap extension penalty"),
         _Option(["-datafile","datafile"], ["file"], None, 0,
                 "Matrix file"),
         _Option(["-similarity","similarity"], [], None, 0,
                 "Display percent identity and similarity"),
         _Option(["-snucleotide","snucleotide"], [], None, 0,
                 "Sequences are nucleotide (boolean)"),
         _Option(["-sprotein","sprotein"], [], None, 0,
                 "Sequences are protein (boolean)"),
         _Option(["-aformat","aformat"], [], None, 0,
                 "Display output in a different specified output format")]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class FuzznucCommandline(_EmbossCommandLine):
    """Commandline object for the fuzznuc program from EMBOSS.
    """
    def __init__(self, cmd="fuzznuc", **kwargs):
        self.parameters = [
         _Option(["-sequence","sequence"], [], None, 1,
                 "Sequence database USA"),
         _Option(["-pattern","pattern"], [], None, 1,
                 "Search pattern, using standard IUPAC one-letter codes"),
         _Option(["-mismatch","mismatch"], [], None, 1,
                 "Number of mismatches"),
         _Option(["-complement","complement"], [], None, 0,
                 "Search complementary strand"),
         _Option(["-rformat","rformat"], [], None, 0,
                 "Specify the report format to output in.")]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class Est2GenomeCommandline(_EmbossCommandLine):
    """Commandline object for the est2genome program from EMBOSS.
    """
    def __init__(self, cmd="est2genome", **kwargs):
        self.parameters = [
         _Option(["-est","est"], [], None, 1,
                 "EST sequence(s)"),
         _Option(["-genome","genome"], [], None, 1,
                 "Genomic sequence"),
         _Option(["-match","match"], [], None, 0, 
                 "Score for matching two bases"),
         _Option(["-mismatch","mismatch"], [], None, 0,
                 "Cost for mismatching two bases"),
         _Option(["-gappenalty","gappenalty"], [], None, 0,
                 "Cost for deleting a single base in either sequence, " + \
                 "excluding introns"),
         _Option(["-intronpenalty","intronpenalty"], [], None, 0,
                 "Cost for an intron, independent of length."),
         _Option(["-splicepenalty","splicepenalty"], [], None, 0,
                 "Cost for an intron, independent of length " + \
                 "and starting/ending on donor-acceptor sites"),
         _Option(["-minscore","minscore"], [], None, 0,
                 "Exclude alignments with scores below this threshold score."),
         _Option(["-reverse","reverse"], [], None, 0,
                 "Reverse the orientation of the EST sequence"),
         _Option(["-splice","splice"], [], None, 0,
                 "Use donor and acceptor splice sites."),
         _Option(["-mode","mode"], [], None, 0,
                 "This determines the comparion mode. 'both', 'forward' " + \
                 "'reverse'"),
         _Option(["-best","best"], [], None, 0,
                 "You can print out all comparisons instead of just the best"),
         _Option(["-space","space"], [], None, 0,
                 "for linear-space recursion."),
         _Option(["-shuffle","shuffle"], [], None, 0,
                 "Shuffle"),
         _Option(["-seed","seed"], [], None, 0,
                 "Random number seed"),
         _Option(["-align","align"], [], None, 0,
                 "Show the alignment."),
         _Option(["-width","width"], [], None, 0,
                 "Alignment width")
        ]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class ETandemCommandline(_EmbossCommandLine):
    """Commandline object for the etandem program from EMBOSS.
    """
    def __init__(self, cmd="etandem", **kwargs):
        self.parameters = [
         _Option(["-sequence","sequence"], ["file"], None, 1,
                 "Sequence"),
         _Option(["-minrepeat","minrepeat"], [], None, 1,
                 "Minimum repeat size"),
         _Option(["-maxrepeat","maxrepeat"], [], None, 1,
                 "Maximum repeat size"),
         _Option(["-threshold","threshold"], [], None, 0,
                 "Threshold score"),
         _Option(["-mismatch","mismatch"], [], None, 0,
                   "Allow N as a mismatch"),
         _Option(["-uniform","uniform"], [], None, 0,
                   "Allow uniform consensus"),
         _Option(["-rformat","rformat"], ["output"], None, 0,
                 "Output report format")]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class EInvertedCommandline(_EmbossCommandLine):
    """Commandline object for the einverted program from EMBOSS.
    """
    def __init__(self, cmd="einverted", **kwargs):
        self.parameters = [
         _Option(["-sequence","sequence"], ["file"], None, 1,
                 "Sequence"),
         _Option(["-gap","gap"], ["file"], None, 1,
                 "Gap penalty"),
         _Option(["-threshold","threshold"], [], None, 1,
                 "Minimum score threshold"),
         _Option(["-match","match"], [], None, 1,
                 "Match score"),
         _Option(["-mismatch","mismatch"], [], None, 1,
                   "Mismatch score"),
         _Option(["-maxrepeat","maxrepeat"], [], None, 0,
                 "Maximum separation between the start and end of repeat"),
         ]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class PalindromeCommandline(_EmbossCommandLine):
    """Commandline object for the palindrome program from EMBOSS.
    """
    def __init__(self, cmd="palindrome", **kwargs):
        self.parameters = [
         _Option(["-sequence","sequence"], ["file"], None, 1,
                 "Sequence"),
         _Option(["-minpallen","minpallen"], [], None, 1,
                 "Minimum palindrome length"),
         _Option(["-maxpallen","maxpallen"], [], None, 1,
                 "Maximum palindrome length"),
         _Option(["-gaplimit","gaplimit"], [], None, 1,
                 "Maximum gap between repeats"),
         _Option(["-nummismatches","nummismatches"], [], None, 1,
                 "Number of mismatches allowed"),
         _Option(["-overlap","overlap"], [], None, 1,
                 "Report overlapping matches"),
         ]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class TranalignCommandline(_EmbossCommandLine):
    """Commandline object for the tranalign program from EMBOSS.
    """
    def __init__(self, cmd="tranalign", **kwargs):
        self.parameters = [
         _Option(["-asequence","asequence"], ["file"], None, 1,
                 "Nucleotide sequences to be aligned."),
         _Option(["-bsequence","bsequence"], ["file"], None, 1,
                 "Protein sequence alignment"),
         _Option(["-outseq","outseq"], ["file"], None, 1,
                 "Output sequence file."),
         _Option(["-table","table"], [], None, 0,
                 "Code to use")]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)


class DiffseqCommandline(_EmbossCommandLine):
    """Commandline object for the diffseq program from EMBOSS.
    """
    def __init__(self, cmd="diffseq", **kwargs):
        self.parameters = [
         _Option(["-asequence","asequence"], ["file"], None, 1,
                 "First sequence to compare"),
         _Option(["-bsequence","bsequence"], ["file"], None, 1,
                 "Second sequence to compare"),
         _Option(["-wordsize","wordsize"], [], None, 1,
                 "Word size to use for comparisons (10 default)"),
         _Option(["-aoutfeat","aoutfeat"], ["file"], None, 1,
                "File for output of first sequence's features"),
         _Option(["-boutfeat","boutfeat"], ["file"], None, 1,
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
         _Option(["-sequence","sequence"], ["file"], None, 1,
                "Protein sequence(s) filename"),
         _Option(["-amino","amino"], [], None, 0),
         _Option(["-lysinemodified","lysinemodified"], [], None, 0),
         _Option(["-disulphides","disulphides"], [], None, 0),
         _Option(["-notermini","notermini"], [], None, 0),
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
         _Option(["-sequence","sequence"], ["file"], None, 0,
                 "Input sequence(s) filename"),
         _Option(["-outseq","outseq"], ["file"], None, 0,
                 "Output sequence file."),
         _Option(["-sformat","sformat"], [], None, 0,
                 "Input sequence(s) format (e.g. fasta, genbank)"),
         _Option(["-osformat","osformat"], [], None, 0,
                 "Output sequence(s) format (e.g. fasta, genbank)"),
         ]
        _EmbossMinimalCommandLine.__init__(self, cmd, **kwargs)

    def _validate(self):
        #Check the outfile, filter, or stdout option has been set.
        #We can't simply do this via the required flag for the outfile
        #output - this seems the simplest solution.
        if not (self.outseq or self.filter or self.stdout):
            raise ValueError("You must either set outfile (output filename), "
                             "or enable filter or stdout (output to stdout).")
        if not (self.sequence or self.filter or self.stdint):
            raise ValueError("You must either set sequence (input filename), "
                             "or enable filter or stdin (input from stdin).")
        return _EmbossMinimalCommandLine._validate(self)

class SeqmatchallCommandline(_EmbossCommandLine):
    """ Commandline object for the seqmatchall program from EMBOSS

    e.g.
    >>> cline = SeqmatchallCommandline(sequence="opuntia.fasta", outfile="opuntia.txt")
    >>> cline.auto = True
    >>> cline.wordsize = 18
    >>> cline.aformat = "pair"
    >>> print cline
    seqmatchall -auto -outfile=opuntia.txt -sequence=opuntia.fasta -wordsize=18 -aformat=pair

    """
    def __init__(self, cmd="seqmatchall", **kwargs):
        self.parameters = [
          _Option(["-sequence", "sequence"], ["file"],
                  None, 1, "Readable set of sequences"),
          _Option(["-wordsize", "wordsize"], [],
                  None, 0, "Word size (Integer 2 or more, default 4)"),
          _Option(["-aformat","aformat"], [], None, 0,
                  "Display output in a different specified output format"),
        ]
        _EmbossCommandLine.__init__(self, cmd, **kwargs)

def _test():
    """Run the Bio.Emboss.Applications module doctests."""
    import doctest
    doctest.testmod(verbose=1)

if __name__ == "__main__":
    #Run the doctests
    _test()
