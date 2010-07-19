# Copyright 2009 by Cymon J. Cox.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Command line wrapper for the multiple alignment program MUSCLE.

http://www.drive5.com/muscle/

Citations:

Edgar, Robert C. (2004), MUSCLE: multiple sequence alignment with high accuracy
and high throughput, Nucleic Acids Research 32(5), 1792-97.

Edgar, R.C. (2004) MUSCLE: a multiple sequence alignment method with reduced
time and space complexity. BMC Bioinformatics 5(1): 113.

Last checked against version: 3.7, briefly against 3.8
"""

from Bio.Application import _Option, _Switch, AbstractCommandline

class MuscleCommandline(AbstractCommandline):
    """Command line wrapper for the multiple alignment program MUSCLE."""
    def __init__(self, cmd="muscle", **kwargs):
        CLUSTERING_ALGORITHMS   = ["upgma", "upgmb", "neighborjoining"]
        DISTANCE_MEASURES_ITER1 = ["kmer6_6", "kmer20_3", "kmer20_4", "kbit20_3",
                                   "kmer4_6"]
        DISTANCE_MEASURES_ITER2 = DISTANCE_MEASURES_ITER1 + \
                                  ["pctid_kimura", "pctid_log"]
        OBJECTIVE_SCORES        = ["sp", "ps", "dp", "xp", "spf", "spm"]
        TREE_ROOT_METHODS       = ["pseudo", "midlongestspan", "minavgleafdist"]
        SEQUENCE_TYPES          = ["protein", "nucleo", "auto"]
        WEIGHTING_SCHEMES       = ["none", "clustalw", "henikoff", "henikoffpb",
                                   "gsc", "threeway"]
        self.parameters = \
           [
            #Can't use "in" as the final alias as this is a reserved word in python:
            _Option(["-in", "in", "input"], ["input", "file"],
                    None, 0, "Input filename",
                    0), #No equate
            _Option(["-out", "out"], ["output", "file"],
                    None, 0, "Output filename",
                    0), #No equate
            _Switch(["-diags", "diags"], ["input"],
                    "Find diagonals (faster for similar sequences)"),
            _Switch(["-profile", "profile"], ["input"],
                    "Perform a profile alignment"),
            _Option(["-in1", "in1"], ["input", "file"],
                    None, 0,
                    "First input filename for profile alignment",
                    0),
            _Option(["-in2", "in2"], ["input", "file"],
                    None, 0,
                    "Second input filename for a profile alignment",
                    0),
            #anchorspacing   Integer              32                 Minimum spacing between
            _Option(["-anchorspacing", "anchorspacing"], ["input"],
                    lambda x: isinstance(x, int),
                    0,
                    "Minimum spacing between anchor columns",
                    0),
            #center          Floating point       [1]                Center parameter.
            #                                                        Should be negative.
            _Option(["-center", "center"], ["input"],
                    lambda x: isinstance(x, float),
                    0,
                    "Center parameter - should be negative",
                    0),
            #cluster1        upgma                upgmb              Clustering method.
            _Option(["-cluster1", "cluster1"], ["input"],
                    lambda x: x in CLUSTERING_ALGORITHMS, 0,
                    "Clustering method used in iteration 1",
                    0),
            #cluster2        upgmb                                   cluster1 is used in
            #                neighborjoining                         iteration 1 and 2,
            #                                                        cluster2 in later
            #                                                        iterations.
            _Option(["-cluster2", "cluster2"], ["input"],
                    lambda x: x in CLUSTERING_ALGORITHMS, 0,
                    "Clustering method used in iteration 2",
                    0),
            #diaglength      Integer              24                 Minimum length of
            #                                                        diagonal.
            _Option(["-diaglength", "diaglength"], ["input"],
                    lambda x: isinstance(x, int),
                    0,
                    "Minimum length of diagonal",
                    0),
            #diagmargin      Integer              5                  Discard this many
            #                                                        positions at ends of
            #                                                        diagonal.
            _Option(["-diagmargin", "diagmargin"], ["input"],
                    lambda x: isinstance(x, int),
                    0,
                    "Discard this many positions at ends of diagonal",
                    0),
            #distance1       kmer6_6              Kmer6_6 (amino) or Distance measure for
            #                kmer20_3             Kmer4_6 (nucleo)   iteration 1.
            #                kmer20_4
            #                kbit20_3
            #                kmer4_6
            _Option(["-distance1", "distance1"], ["input"],
                    lambda x: x in DISTANCE_MEASURES_ITER1,
                    0,
                    "Distance measure for iteration 1",
                    0),
            #distance2       kmer6_6              pctid_kimura       Distance measure for
            #                kmer20_3                                iterations 2, 3 ...
            #                kmer20_4
            #                kbit20_3
            #                pctid_kimura
            #                pctid_log
            _Option(["-distance2", "distance2"], ["input"],
                    lambda x: x in DISTANCE_MEASURES_ITER2,
                    0,
                    "Distance measure for iteration 2",
                    0),
            #gapopen         Floating point       [1]                The gap open score.
            #                                                        Must be negative.
            _Option(["-gapopen", "gapopen"], ["input"],
                    lambda x: isinstance(x, float),
                    0,
                    "Gap open score - negative number",
                    0),
            #hydro           Integer              5                  Window size for
            #                                                        determining whether a
            #                                                        region is hydrophobic.
            _Option(["-hydro", "hydro"], ["input"],
                    lambda x: isinstance(x, int),
                    0,
                    "Window size for hydrophobic region",
                    0),
            #hydrofactor     Floating point       1.2                Multiplier for gap
            #                                                        open/close penalties in
            #                                                        hydrophobic regions.
            _Option(["-hydrofactor", "hydrofactor"], ["input"],
                    lambda x: isinstance(x, float),
                    0,
                    "Multiplier for gap penalties in hydrophobic regions",
                    0),
            #log             File name            None.              Log file name (delete
            #                                                        existing file).
            _Option(["-log", "log"], ["output", "file"],
                    None, 0,
                    "Log file name",
                    0),
            #loga            File name            None.              Log file name (append
            #                                                        to existing file).
            _Option(["-loga", "loga"], ["output", "file"],
                    None, 0,
                    "Log file name (append to existing file)",
                    0),
            #maxdiagbreak    Integer              1                  Maximum distance
            #                                                        between two diagonals
            #                                                        that allows them to
            #                                                        merge into one
            #                                                        diagonal.
            _Option(["-maxdiagbreak", "maxdiagbreak"], ["input"],
                    lambda x: isinstance(x, int),
                    0,
                    "Maximum distance between two diagonals that allows " + \
                    "them to merge into one diagonal",
                    0),
            #maxhours        Floating point       None.              Maximum time to run in
            #                                                        hours. The actual time
            #                                                        may exceed the
            #                                                        requested limit by a
            #                                                        few minutes. Decimals
            #                                                        are allowed, so 1.5
            #                                                        means one hour and 30
            #                                                        minutes.
            _Option(["-maxhours", "maxhours"], ["input"],
                    lambda x: isinstance(x, float),
                    0,
                    "Maximum time to run in hours",
                    0),
            #maxiters        Integer 1, 2 ...     16                 Maximum number of
            #                                                        iterations.
            _Option(["-maxiters", "maxiters"], ["input"],
                    lambda x: isinstance(x, int),
                    0,
                    "Maximum number of iterations",
                    0),
            #maxtrees        Integer              1                  Maximum number of new
            #                                                        trees to build in
            #                                                        iteration 2.
            _Option(["-maxtrees", "maxtrees"], ["input"],
                    lambda x: isinstance(x, int),
                    0,
                    "Maximum number of trees to build in iteration 2",
                    0),
            #minbestcolscore Floating point       [1]                Minimum score a column
            #                                                        must have to be an
            #                                                        anchor.
            _Option(["-minbestcolscore", "minbestcolscore"], ["input"],
                    lambda x: isinstance(x, float),
                    0,
                    "Minimum score a column must have to be an anchor",
                    0),
            #minsmoothscore  Floating point       [1]                Minimum smoothed score
            #                                                        a column must have to
            #                                                        be an anchor.
            _Option(["-minsmoothscore", "minsmoothscore"], ["input"],
                    lambda x: isinstance(x, float),
                    0,
                    "Minimum smoothed score a column must have to " + \
                    "be an anchor",
                    0),
            #objscore        sp                   spm                Objective score used by
            #                ps                                      tree dependent
            #                dp                                      refinement.
            #                xp                                      sp=sum-of-pairs score.
            #                spf                                     spf=sum-of-pairs score
            #                spm                                     (dimer approximation)
            #                                                        spm=sp for < 100 seqs,
            #                                                        otherwise spf
            #                                                        dp=dynamic programming
            #                                                        score.
            #                                                        ps=average profile-
            #                                                        sequence score.
            #                                                        xp=cross profile score.
            _Option(["-objscore", "objscore"], ["input"],
                    lambda x: x in OBJECTIVE_SCORES,
                    0,
                    "Objective score used by tree dependent refinement",
                    0),
            #root1           pseudo               psuedo             Method used to root
            _Option(["-root1", "root1"], ["input"],
                    lambda x: x in TREE_ROOT_METHODS,
                    0,
                    "Method used to root tree in iteration 1",
                    0),
            #root2           midlongestspan                          tree; root1 is used in
            #                minavgleafdist                          iteration 1 and 2,
            #                                                        root2 in later
            #                                                        iterations.
            _Option(["-root2", "root2"], ["input"],
                    lambda x: x in TREE_ROOT_METHODS,
                    0,
                    "Method used to root tree in iteration 2",
                    0),
            #seqtype         protein              auto               Sequence type.
            #                nucleo
            #                auto
            _Option(["-seqtype", "seqtype"], ["input"],
                    lambda x: x in SEQUENCE_TYPES,
                    0,
                    "Sequence type",
                    0),
            #smoothscoreceil Floating point       [1]                Maximum value of column
            #                                                        score for smoothing
            #                                                        purposes.
            _Option(["-smoothscoreceil", "smoothscoreceil"], ["input"],
                    lambda x: isinstance(x, float),
                    0,
                    "Maximum value of column score for smoothing",
                    0),
            #smoothwindow    Integer              7                  Window used for anchor
            #                                                        column smoothing.
            _Option(["-smoothwindow", "smoothwindow"], ["input"],
                    lambda x: isinstance(x, int),
                    0,
                    "Window used for anchor column smoothing",
                    0),
            #SUEFF           Floating point value 0.1                Constant used in UPGMB
            #                between 0 and 1.                        clustering. Determines
            #                                                        the relative fraction
            #                                                        of average linkage
            #                                                        (SUEFF) vs. nearest-
            #                                                        neighbor linkage (1
            #                                                        SUEFF).
            _Option(["-sueff", "sueff"], ["input"],
                    lambda x: isinstance(x, float),
                    0,
                    "Constant used in UPGMB clustering",
                    0),
            #tree1           File name            None               Save tree produced in
            _Option(["-tree1", "tree1"], ["input"],
                    None, 0,
                    "Save Newick tree from iteration 1",
                    0),
            #tree2                                                   first or second
            #                                                        iteration to given file
            #                                                        in Newick (Phylip-
            #                                                        compatible) format.
            _Option(["-tree2", "tree2"], ["input"],
                    None, 0,
                    "Save Newick tree from iteration 2",
                    0),
            #weight1         none                 clustalw           Sequence weighting
            _Option(["-weight1", "weight1"], ["input"],
                    lambda x: x in WEIGHTING_SCHEMES,
                    0,
                    "Weighting scheme used in iteration 1",
                    0),
            #weight2         henikoff                                scheme.
            #                henikoffpb                              weight1 is used in
            #                gsc                                     iterations 1 and 2.
            #                clustalw                                weight2 is used for
            #                threeway                                tree-dependent
            #                                                        refinement.
            #                                                        none=all sequences have
            #                                                        equal weight.
            #                                                        henikoff=Henikoff &
            #                                                        Henikoff weighting
            #                                                        scheme.
            #                                                        henikoffpb=Modified
            #                                                        Henikoff scheme as used
            #                                                        in PSI-BLAST.
            #                                                        clustalw=CLUSTALW
            #                                                        method.
            #                                                        threeway=Gotoh three-
            #                                                        way method.
            _Option(["-weight2", "weight2"], ["input"],
                    lambda x: x in WEIGHTING_SCHEMES,
                    0,
                    "Weighting scheme used in iteration 2",
                    0),
            #################### FORMATS #######################################
            # Multiple formats can be specified on the command line
            # If -msf appears it will be used regardless of other formats
            # specified. If -clw appears (and not -msf), clustalw format will be
            # used regardless of other formats specified. If both -clw and
            # -clwstrict are specified -clwstrict will be used regardless of
            # other formats specified. If -fasta is specified and not -msf,
            # -clw, or clwstrict, fasta will be used. If -fasta and -html are
            # specified -fasta will be used. Only if -html is specified alone
            # will html be used. I kid ye not.
            #clw                no              Write output in CLUSTALW format (default is
            #                                   FASTA).
            _Switch(["-clw", "clw"], ["input"],
                    "Write output in CLUSTALW format (with a MUSCLE header)"),
            #clwstrict          no              Write output in CLUSTALW format with the
            #                                   "CLUSTAL W (1.81)" header rather than the
            #                                   MUSCLE version. This is useful when a post-
            #                                   processing step is picky about the file
            #                                   header.
            _Switch(["-clwstrict", "clwstrict"], ["input"],
                    "Write output in CLUSTALW format with version 1.81 header"),
            #fasta              yes             Write output in FASTA format. Alternatives
            #                                   include clw,
            #                                   clwstrict, msf and html.
            _Switch(["-fasta", "fasta"], ["input"],
                    "Write output in FASTA format"),
            #html               no              Write output in HTML format (default is
            #                                   FASTA).
            _Switch(["-html", "html"], ["input"],
                    "Write output in HTML format"),
            #msf                no              Write output in MSF format (default is
            #                                   FASTA).
            _Switch(["-msf", "msf"], ["input"],
                    "Write output in MSF format"),
            #Phylip interleaved - undocumented as of 3.7
            _Switch(["-phyi", "phyi"], ["input"],
                    "Write output in PHYLIP interleaved format"),
            #Phylip sequential - undocumented as of 3.7
            _Switch(["-phys", "phys"], ["input"],
                    "Write output in PHYLIP sequential format"),
            ################## Additional specified output files #########
            _Option(["-phyiout", "phyiout"], ["output", "file"],
                    None, 0,
                    "Write PHYLIP interleaved output to specified filename",
                    0), #No equate
            _Option(["-physout", "physout"], ["output", "file"],
                    None, 0,
                    "Write PHYLIP sequential format to specified filename",
                    0), #No equate
            _Option(["-htmlout", "htmlout"], ["output", "file"],
                    None, 0,
                    "Write HTML output to specified filename",
                    0), #No equate
            _Option(["-clwout", "clwout"], ["output", "file"],
                    None, 0,
                    "Write CLUSTALW output (with MUSCLE header) to specified "
                    "filename",
                    0), #No equate
            _Option(["-clwstrictout", "clwstrictout"], ["output", "file"],
                    None, 0,
                    "Write CLUSTALW output (with version 1.81 header) to "
                    "specified filename",
                    0), #No equate
            _Option(["-msfout", "msfout"], ["output", "file"],
                    None, 0,
                    "Write MSF format output to specified filename",
                    0), #No equate
            _Option(["-fastaout", "fastaout"], ["output", "file"],
                    None, 0,
                    "Write FASTA format output to specified filename",
                    0), #No equate
            ############## END FORMATS ###################################
            #anchors            yes             Use anchor optimization in tree dependent
            #                                   refinement iterations.
            _Switch(["-anchors", "anchors"], ["input"],
                    "Use anchor optimisation in tree dependent " + \
                    "refinement iterations"),
            #noanchors          no              Disable anchor optimization. Default is
            #                                   anchors.
            _Switch(["-noanchors", "noanchors"], ["input"],
                    "Do not use anchor optimisation in tree dependent " + \
                    "refinement iterations"),
            #group              yes             Group similar sequences together in the
            #                                   output. This is the default. See also
            #                                   stable.
            _Switch(["-group", "group"], ["input"],
                    "Group similar sequences in output"),
            #stable             no              Preserve input order of sequences in output
            #                                   file. Default is to group sequences by
            #                                   similarity (group).
            _Switch(["-stable", "stable"], ["input"],
                    "Do not group similar sequences in output (not supported in v3.8)"),
            ############## log-expectation profile score ######################
            # One of either -le, -sp, or -sv
            #
            # According to the doc, spn is default and the only option for
            # nucleotides: this doesnt appear to be true. -le, -sp, and -sv can
            # be used and produce numerically different logs (what is going on?)
            #
            #spn fails on proteins
            #le                 maybe           Use log-expectation profile score (VTML240).
            #                                    Alternatives are to use sp or sv. This is
            #                                    the default for amino acid sequences.
            _Switch(["-le", "le"], ["input"],
                    "Use log-expectation profile score (VTML240)"),
            #sv                 no              Use sum-of-pairs profile score (VTML240).
            #                                   Default is le.
            _Switch(["-sv", "sv"], ["input"],
                    "Use sum-of-pairs profile score (VTML240)"),
            #sp                 no              Use sum-of-pairs protein profile score
            #                                   (PAM200). Default is le.
            _Switch(["-sp", "sp"], ["input"],
                    "Use sum-of-pairs protein profile score (PAM200)"),
            #spn                maybe           Use sum-of-pairs nucleotide profile score
            #                                   (BLASTZ parameters). This is the only option
            #                                   for nucleotides, and is therefore the
            #                                   default.
            _Switch(["-spn", "spn"], ["input"],
                    "Use sum-of-pairs protein nucleotide profile score"),
            ############## END log-expectation profile score ######################
            #quiet              no              Do not display progress messages.
            _Switch(["-quiet", "quiet"], ["input"],
                    "Use sum-of-pairs protein nucleotide profile score"),
            #refine             no              Input file is already aligned, skip first
            #                                   two iterations and begin tree dependent
            #                                   refinement.
            _Switch(["-refine", "refine"], ["input"],
                    "Only do tree dependent refinement"),
            #core               yes in muscle,  Do not catch exceptions.
            #                   no in muscled.
            _Switch(["-core", "core"], ["input"],
                    "Catch exceptions"),
            #nocore             no in muscle,   Catch exceptions and give an error message
            #                   yes in muscled. if possible.
            _Switch(["-nocore", "nocore"], ["input"],
                    "Do not catch exceptions"),
            #termgapsfull       no              Terminal gaps penalized with full penalty.
            #                                   [1] Not fully supported in this version.
            #
            #termgapshalf       yes             Terminal gaps penalized with half penalty.
            #                                   [1] Not fully supported in this version.
            #
            #termgapshalflonger no              Terminal gaps penalized with half penalty if
            #                                   gap relative to
            #                                   longer sequence, otherwise with full
            #                                   penalty.
            #                                   [1] Not fully supported in this version.
            #verbose            no              Write parameter settings and progress
            #                                   messages to log file.
            _Switch(["-verbose", "verbose"], ["input"],
                    "Write parameter settings and progress"),
            #version            no              Write version string to stdout and exit.
            _Switch(["-version", "version"], ["input"],
                    "Write version string to stdout and exit"),
           ]
        AbstractCommandline.__init__(self, cmd, **kwargs)
