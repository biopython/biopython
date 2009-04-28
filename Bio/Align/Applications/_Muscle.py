# Copyright 2009 by Cymon J. Cox.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
Bio.Application command line for the multiple alignment programme MUSCLE

http://www.drive5.com/muscle/

Citations:

Edgar, Robert C. (2004), MUSCLE: multiple sequence alignment with high accuracy
and high throughput, Nucleic Acids Research 32(5), 1792-97. 

Edgar, R.C. (2004) MUSCLE: a multiple sequence alignment method with reduced
time and space complexity. BMC Bioinformatics 5(1): 113.

Last checked against version: 3.7
"""

import types
from Bio import Application
from Bio.Application import _Option
from Bio.Application import _Argument

class MuscleCommandline(Application.AbstractCommandline):

    def __init__(self, cmd = "muscle"):

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

        Application.AbstractCommandline.__init__(self)
        self.program_name = cmd
        self.parameters = \
           [
            _Option(["-in", "in"], ["input", "file"],
                    None, 0, "Input filename",
                    0), #No equate

            _Option(["-out", "out"], ["output", "file"],
                    None, 0, "Output filename",
                    0), #No equate

            _Option(["-diags", "diags"], ["input"],
                    lambda x: 0, #Does not take a value
                    0, "Find diagonals (faster for similar sequences)",
                    0), #No equate
            
            _Option(["-profile", "profile"], ["input"],
                    lambda x: 0, #Does not take a value
                    0, "Perform a profile alignment",
                    0), #No equate

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
                    lambda x: isinstance(x, types.IntType),
                    0,
                    "Minimum spacing between anchor columns",
                    0),

            #center          Floating point       [1]                Center parameter.
            #                                                        Should be negative.
            _Option(["-center", "center"], ["input"],
                    lambda x: isinstance(x, types.FloatType),
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
                    lambda x: isinstance(x, types.IntType),
                    0,
                    "Minimum length of diagonal",
                    0),

            #diagmargin      Integer              5                  Discard this many
            #                                                        positions at ends of
            #                                                        diagonal.
            _Option(["-diagmargin", "diagmargin"], ["input"],
                    lambda x: isinstance(x, types.IntType),
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
                    lambda x: isinstance(x, types.FloatType),
                    0,
                    "Gap open score - negative number",
                    0),

            #hydro           Integer              5                  Window size for
            #                                                        determining whether a
            #                                                        region is hydrophobic.
            _Option(["-hydro", "hydro"], ["input"],
                    lambda x: isinstance(x, types.IntType),
                    0,
                    "Window size for hydrophobic region",
                    0),

            #hydrofactor     Floating point       1.2                Multiplier for gap
            #                                                        open/close penalties in
            #                                                        hydrophobic regions.
            _Option(["-hydrofactor", "hydrofactor"], ["input"],
                    lambda x: isinstance(x, types.FloatType),
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
                    lambda x: isinstance(x, types.IntType),
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
                    lambda x: isinstance(x, types.FloatType),
                    0,
                    "Maximum time to run in hours",
                    0),

            #maxiters        Integer 1, 2 ...     16                 Maximum number of
            #                                                        iterations.
            _Option(["-maxiters", "maxiters"], ["input"],
                    lambda x: isinstance(x, types.IntType),
                    0,
                    "Maximum number of iterations",
                    0),

            #maxtrees        Integer              1                  Maximum number of new
            #                                                        trees to build in
            #                                                        iteration 2.
            _Option(["-maxtrees", "maxtrees"], ["input"],
                    lambda x: isinstance(x, types.IntType),
                    0,
                    "Maximum number of trees to build in iteration 2",
                    0),

            #minbestcolscore Floating point       [1]                Minimum score a column
            #                                                        must have to be an
            #                                                        anchor.
            _Option(["-minbestcolscore", "minbestcolscore"], ["input"],
                    lambda x: isinstance(x, types.FloatType),
                    0,
                    "Minimum score a column must have to be an anchor",
                    0),

            #minsmoothscore  Floating point       [1]                Minimum smoothed score
            #                                                        a column must have to
            #                                                        be an anchor.
            _Option(["-minsmoothscore", "minsmoothscore"], ["input"],
                    lambda x: isinstance(x, types.FloatType),
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
                    lambda x: isinstance(x, types.FloatType),
                    0,
                    "Maximum value of column score for smoothing",
                    0),

            #smoothwindow    Integer              7                  Window used for anchor
            #                                                        column smoothing.
            _Option(["-smoothwindow", "smoothwindow"], ["input"],
                    lambda x: isinstance(x, types.IntType),
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
                    lambda x: isinstance(x, types.FloatType),
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
            _Option(["-clw", "clw"], ["input"],
                    lambda x: 0, #Does not take a value,
                    0,
                    "Write output in CLUSTALW format",
                    0),

            #clwstrict          no              Write output in CLUSTALW format with the
            #                                   "CLUSTAL W (1.81)" header rather than the
            #                                   MUSCLE version. This is useful when a post-
            #                                   processing step is picky about the file
            #                                   header.
            _Option(["-clwstrict", "clwstrict"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Write output in CLUSTALW format with vers. 1.81 header",
                    0),

            #fasta              yes             Write output in FASTA format. Alternatives
            #                                   include clw,
            #                                   clwstrict, msf and html.
            _Option(["-fasta", "fasta"], ["input"],
                    lambda x: 0, #Does not take a value 
                    0,
                    "Write output in FASTA format",
                    0),

            #html               no              Write output in HTML format (default is
            #                                   FASTA).
            _Option(["-html", "html"], ["input"], 
                    lambda x: 0, #Does not take a value 
                    0,
                    "Write output in HTML format",
                    0),

            #msf                no              Write output in MSF format (default is
            #                                   FASTA).
            _Option(["-msf", "msf"], ["input"], 
                    lambda x: 0, #Does not take a value 
                    0,
                    "Write output in MSF format",
                    0),
            ############## END FORMATS ###################################

            #anchors            yes             Use anchor optimization in tree dependent
            #                                   refinement iterations.
            _Option(["-anchors", "anchors"], ["input"], 
                    lambda x: 0, #Does not take a value 
                    0,
                    "Use anchor optimisation in tree dependent " + \
                    "refinement iterations",
                    0),

            #noanchors          no              Disable anchor optimization. Default is
            #                                   anchors.
            _Option(["-noanchors", "noanchors"], ["input"], 
                    lambda x: 0, #Does not take a value 
                    0,
                    "Do not use anchor optimisation in tree dependent " + \
                    "refinement iterations",
                    0),

            #group              yes             Group similar sequences together in the
            #                                   output. This is the default. See also
            #                                   stable.
            _Option(["-group", "group"], ["input"], 
                    lambda x: 0, #Does not take a value 
                    0,
                    "Group similar sequences in output",
                    0),

            #stable             no              Preserve input order of sequences in output
            #                                   file. Default is to group sequences by
            #                                   similarity (group).
            _Option(["-stable", "stable"], ["input"], 
                    lambda x: 0, #Does not take a value 
                    0,
                    "Do not group similar sequences in output",
                    0),

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
            _Option(["-le", "le"], ["input"], 
                    lambda x: 0, #Does not take a value 
                    0,
                    "Use log-expectation profile score (VTML240)",
                    0),

            #sv                 no              Use sum-of-pairs profile score (VTML240).
            #                                   Default is le.
            _Option(["-sv", "sv"], ["input"], 
                    lambda x: 0, #Does not take a value 
                    0,
                    "Use sum-of-pairs profile score (VTML240)",
                    0),

            #sp                 no              Use sum-of-pairs protein profile score
            #                                   (PAM200). Default is le.
            _Option(["-sp", "sp"], ["input"], 
                    lambda x: 0, #Does not take a value 
                    0,
                    "Use sum-of-pairs protein profile score (PAM200)",
                    0),

            #spn                maybe           Use sum-of-pairs nucleotide profile score
            #                                   (BLASTZ parameters). This is the only option
            #                                   for nucleotides, and is therefore the
            #                                   default.
            _Option(["-spn", "spn"], ["input"], 
                    lambda x: 0, #Does not take a value 
                    0,
                    "Use sum-of-pairs protein nucleotide profile score",
                    0),
            ############## END log-expectation profile score ######################

            #quiet              no              Do not display progress messages.
            _Option(["-quiet", "quiet"], ["input"], 
                    lambda x: 0, #Does not take a value 
                    0,
                    "Use sum-of-pairs protein nucleotide profile score",
                    0),

            #refine             no              Input file is already aligned, skip first
            #                                   two iterations and begin tree dependent
            #                                   refinement.
            _Option(["-refine", "refine"], ["input"],
                    lambda x: 0, #Does not take a value 
                    0,
                    "Only do tree dependent refinement",
                    0),

            #core               yes in muscle,  Do not catch exceptions.
            #                   no in muscled.
            _Option(["-core", "core"], ["input"], 
                    lambda x: 0, #Does not take a value 
                    0,
                    "Catch exceptions",
                    0),

            #nocore             no in muscle,   Catch exceptions and give an error message
            #                   yes in muscled. if possible.
            _Option(["-nocore", "nocore"], ["input"],
                    lambda x: 0, #Does not take a value 
                    0,
                    "Do not catch exceptions",
                    0),

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

            _Option(["-verbose", "verbose"], ["input"],
                    lambda x: 0, #Does not take a value
                    0,
                    "Write parameter settings and progress",
                    0),

            #version            no              Write version string to stdout and exit.
            _Option(["-version", "version"], ["input"],
                    lambda x: 0, #Does not take a value,
                    0,
                    "Write version string to stdout and exit",
                    0)
           ]

