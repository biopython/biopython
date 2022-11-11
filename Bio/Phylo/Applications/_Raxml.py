# Copyright 2012 by Eric Talevich.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Command-line wrapper for the tree inference program RAxML.

Derived from the help page for RAxML version 7.3 by Alexandros Stamatakis, but
should work for any version 7.X (and probably earlier for most options).
"""

from Bio.Application import _Option, _Switch, AbstractCommandline


class RaxmlCommandline(AbstractCommandline):
    """Command-line wrapper for the tree inference program RAxML.

    The required parameters are 'sequences' (-s), 'model' (-m) and 'name' (-n).
    The parameter 'parsimony_seed' (-p) must also be set for RAxML, but if you
    do not specify it, this wrapper will set the seed to 10000 for you.

    References
    ----------
    Stamatakis A.
    RAxML-VI-HPC: Maximum Likelihood-based Phylogenetic Analyses with
    Thousands of Taxa and Mixed Models.
    Bioinformatics 2006, 22(21):2688-2690.

    Homepage: http://sco.h-its.org/exelixis/software.html

    Examples
    --------
    >>> from Bio.Phylo.Applications import RaxmlCommandline
    >>> raxml_cline = RaxmlCommandline(sequences="Tests/Phylip/interlaced2.phy",
    ...                                model="PROTCATWAG", name="interlaced2")
    >>> print(raxml_cline)
    raxmlHPC -m PROTCATWAG -n interlaced2 -p 10000 -s Tests/Phylip/interlaced2.phy

    You would typically run the command line with raxml_cline() or via
    the Python subprocess module, as described in the Biopython tutorial.

    """

    def __init__(self, cmd="raxmlHPC", **kwargs):
        """Initialize the class."""
        self.parameters = [
            _Option(
                ["-a", "weight_filename"],
                "Name of a column weight file to assign individual weights "
                "to each column of the alignment. Those weights must be "
                "integers separated by any type and number of whitespaces "
                "within a separate file.",
                filename=True,
                equate=False,
            ),
            _Option(
                ["-b", "bootstrap_seed"], "Random seed for bootstrapping.", equate=False
            ),
            _Option(
                ["-c", "num_categories"],
                "Number of distinct rate categories for RAxML when "
                "evolution model is set to GTRCAT or GTRMIX."
                "Individual per-site rates are categorized into this "
                "many rate categories to accelerate computations. "
                "Default: 25.",
                equate=False,
            ),
            _Switch(
                ["-d", "random_starting_tree"],
                "Start ML optimization from random starting tree.",
            ),
            _Option(
                ["-e", "epsilon"],
                "Set model optimization precision in log likelihood units "
                "for final optimization of tree topology under MIX/MIXI "
                "or GAMMA/GAMMAI."
                "Default: 0.1 for models not using proportion of "
                "invariant sites estimate; 0.001 for models using "
                "proportion of invariant sites estimate.",
                equate=False,
            ),
            _Option(
                ["-E", "exclude_filename"],
                "An exclude file name, containing a specification of "
                "alignment positions you wish to exclude.  Format is "
                "similar to Nexus, the file shall contain entries like "
                "'100-200 300-400'; to exclude a single column write, "
                "e.g., '100-100'. If you use a mixed model, an "
                "appropriately adapted model file will be written.",
                filename=True,
                equate=False,
            ),
            _Option(
                ["-f", "algorithm"],
                r"""
                        Select algorithm:

                        a: Rapid Bootstrap analysis and search for best-scoring ML
                        tree in one program run.

                        b: Draw bipartition information on a tree provided with '-t'
                        based on multiple trees (e.g. form a bootstrap) in a file
                        specified by '-z'.

                        c: Check if the alignment can be properly read by RAxML.

                        d: New rapid hill-climbing (DEFAULT).

                        e: Optimize model+branch lengths for given input tree under
                        GAMMA/GAMMAI only.

                        g: Compute per site log Likelihoods for one ore more trees
                        passed via '-z' and write them to a file that can be read
                        by CONSEL.

                        h: Compute log likelihood test (SH-test) between best tree
                        passed via '-t' and a bunch of other trees passed via '-z'.

                        i: Perform a really thorough bootstrap, refinement of final
                        bootstrap tree under GAMMA and a more exhaustive algorithm.

                        j: Generate a bunch of bootstrapped alignment files from an
                        original alignment file.

                        m: Compare bipartitions between two bunches of trees passed
                        via '-t' and '-z' respectively. This will return the
                        Pearson correlation between all bipartitions found in the
                        two tree files. A file called
                        RAxML_bipartitionFrequencies.outputFileName will be
                        printed that contains the pair-wise bipartition
                        frequencies of the two sets.

                        n: Compute the log likelihood score of all trees contained
                        in a tree file provided by '-z' under GAMMA or
                        GAMMA+P-Invar.

                        o: Old and slower rapid hill-climbing.

                        p: Perform pure stepwise MP addition of new sequences to an
                        incomplete starting tree.

                        s: Split up a multi-gene partitioned alignment into the
                        respective subalignments.

                        t: Do randomized tree searches on one fixed starting tree.

                        w: Compute ELW test on a bunch of trees passed via '-z'.

                        x: Compute pair-wise ML distances, ML model parameters will
                        be estimated on an MP starting tree or a user-defined
                        tree passed via '-t', only allowed for GAMMA-based models
                        of rate heterogeneity.
                        """,
                checker_function=(lambda x: isinstance(x, str) and len(x) == 1),
                equate=False,
            ),
            _Option(
                ["-g", "grouping_constraint"],
                "File name of a multifurcating constraint tree. "
                "this tree does not need to be comprehensive, i.e. "
                "contain all taxa.",
                filename=True,
                equate=False,
            ),
            _Option(
                ["-i", "rearrangements"],
                "Initial rearrangement setting for the subsequent "
                "application of topological changes phase.",
                equate=False,
            ),
            _Switch(
                ["-j", "checkpoints"],
                "Write checkpoints (intermediate tree topologies).",
            ),
            _Switch(
                ["-k", "bootstrap_branch_lengths"],
                "Print bootstrapped trees with branch lengths. "
                "The bootstraps will run a bit longer, because model "
                "parameters will be optimized at the end of each run. "
                "Use with CATMIX/PROTMIX or GAMMA/GAMMAI.",
            ),
            _Option(
                ["-l", "cluster_threshold"],
                "Threshold for sequence similarity clustering. "
                "RAxML will then print out an alignment to a file "
                "called sequenceFileName.reducedBy.threshold that "
                "only contains sequences <= the specified threshold "
                "that must be between 0.0 and 1.0. RAxML uses the "
                "QT-clustering algorithm to perform this task. "
                "In addition, a file called "
                "RAxML_reducedList.outputFileName will be written "
                "that contains clustering information.",
                equate=False,
            ),
            _Option(
                ["-L", "cluster_threshold_fast"],
                "Same functionality as '-l', but uses a less "
                "exhaustive and thus faster clustering algorithm. "
                "This is intended for very large datasets with more "
                "than 20,000-30,000 sequences.",
                equate=False,
            ),
            _Option(
                ["-m", "model"],
                r"""Model of Nucleotide or Amino Acid Substitution:

                        NUCLEOTIDES:

                        GTRCAT         : GTR + Optimization of substitution rates + Optimization of site-specific
                        evolutionary rates which are categorized into numberOfCategories distinct
                        rate categories for greater computational efficiency
                        if you do a multiple analysis with  '-#' or '-N' but without bootstrapping the program
                        will use GTRMIX instead

                        GTRGAMMA       : GTR + Optimization of substitution rates + GAMMA model of rate
                        heterogeneity (alpha parameter will be estimated)

                        GTRMIX         : Inference of the tree under GTRCAT
                        and thereafter evaluation of the final tree topology under GTRGAMMA

                        GTRCAT_GAMMA   : Inference of the tree with site-specific evolutionary rates.
                        However, here rates are categorized using the 4 discrete GAMMA rates.
                        Evaluation of the final tree topology under GTRGAMMA

                        GTRGAMMAI      : Same as GTRGAMMA, but with estimate of proportion of invariable sites

                        GTRMIXI        : Same as GTRMIX, but with estimate of proportion of invariable sites

                        GTRCAT_GAMMAI  : Same as GTRCAT_GAMMA, but with estimate of proportion of invariable sites

                        AMINO ACIDS:

                        PROTCATmatrixName[F]    : specified AA matrix + Optimization of substitution rates + Optimization of site-specific
                        evolutionary rates which are categorized into numberOfCategories distinct
                        rate categories for greater computational efficiency
                        if you do a multiple analysis with  '-#' or '-N' but without bootstrapping the program
                        will use PROTMIX... instead

                        PROTGAMMAmatrixName[F]  : specified AA matrix + Optimization of substitution rates + GAMMA model of rate
                        heterogeneity (alpha parameter will be estimated)

                        PROTMIXmatrixName[F]    : Inference of the tree under specified AA matrix + CAT
                        and thereafter evaluation of the final tree topology under specified AA matrix + GAMMA

                        PROTCAT_GAMMAmatrixName[F] : Inference of the tree under specified AA matrix and site-specific evolutionary rates.
                        However, here rates are categorized using the 4 discrete GAMMA rates.
                        Evaluation of the final tree topology under specified AA matrix + GAMMA

                        PROTGAMMAImatrixName[F] : Same as PROTGAMMAmatrixName[F], but with estimate of proportion of invariable sites

                        PROTMIXImatrixName[F]   : Same as PROTMIXmatrixName[F], but with estimate of proportion of invariable sites

                        PROTCAT_GAMMAImatrixName[F] : Same as PROTCAT_GAMMAmatrixName[F], but with estimate of proportion of invariable sites

                        Available AA substitution models: DAYHOFF, DCMUT, JTT, MTREV, WAG, RTREV, CPREV, VT, BLOSUM62, MTMAM, GTR
                        With the optional 'F' appendix you can specify if you want to use empirical base frequencies
                        Please not that for mixed models you can in addition specify the per-gene AA model in
                        the mixed model file (see manual for details)
                        """,
                equate=False,
            ),
            _Switch(
                ["-M", "partition_branch_lengths"],
                "Switch on estimation of individual per-partition "
                "branch lengths. Only has effect when used in "
                "combination with 'partition_filename' ('-q'). "
                "Branch lengths for individual partitions will be "
                "printed to separate files.  A weighted average of the "
                "branch lengths is computed by using the respective "
                "partition lengths. ",
            ),
            _Option(
                ["-n", "name"],
                "Name used in the output files.",
                filename=True,
                equate=False,
            ),
            _Option(
                ["-o", "outgroup"],
                "Name of a single outgroup or a comma-separated list "
                "of outgroups, eg '-o Rat' or '-o Rat,Mouse'. In case "
                "that multiple outgroups are not monophyletic the "
                "first name in the list will be selected as outgroup. "
                "Don't leave spaces between taxon names!",
                checker_function=lambda x: len(x.split()) == 1,
                equate=False,
            ),
            _Option(
                ["-q", "partition_filename"],
                "File name containing the assignment of models to "
                "alignment partitions for multiple models of "
                "substitution. For the syntax of this file please "
                "consult the RAxML manual.",
                filename=True,
                equate=False,
            ),
            _Option(
                ["-p", "parsimony_seed"],
                "Random number seed for the parsimony inferences. "
                "This allows you to reproduce your results and will "
                "help developers debug the program. This option HAS "
                "NO EFFECT in the parallel MPI version.",
                equate=False,
            ),
            _Option(
                ["-P", "protein_model"],
                "File name of a user-defined AA (Protein) substitution "
                "model. This file must contain 420 entries, the first "
                "400 being the AA substitution rates (this must be a "
                "symmetric matrix) and the last 20 are the empirical "
                "base frequencies.",
                filename=True,
                equate=False,
            ),
            _Option(
                ["-r", "binary_constraint"],
                "File name of a binary constraint tree. "
                "This tree does not need to be comprehensive, i.e. "
                "contain all taxa.",
                filename=True,
                equate=False,
            ),
            _Option(
                ["-s", "sequences"],
                "Name of the alignment data file, in PHYLIP format.",
                filename=True,
                equate=False,
            ),
            _Option(
                ["-t", "starting_tree"],
                "File name of a user starting tree, in Newick format.",
                filename=True,
                equate=False,
            ),
            _Option(
                ["-T", "threads"],
                "Number of threads to run. "
                "PTHREADS VERSION ONLY! "
                "Make sure to set this at most the number of CPUs "
                "you have on your machine, otherwise, there will be "
                "a huge performance decrease!",
                equate=False,
            ),
            _Option(
                ["-u", "num_bootstrap_searches"],
                "Number of multiple bootstrap searches per replicate. "
                "Use this to obtain better ML trees for each "
                "replicate. Default: 1 ML search per bootstrap "
                "replicate.",
                equate=False,
            ),
            _Switch(["-v", "version"], "Display version information."),
            _Option(
                ["-w", "working_dir"],
                "Name of the working directory where RAxML will "
                "write its output files. Default: current directory.",
                filename=True,
                equate=False,
            ),
            _Option(
                ["-x", "rapid_bootstrap_seed"],
                "Random seed for rapid bootstrapping.",
                equate=False,
            ),
            _Switch(
                ["-y", "parsimony"],
                "Only compute a parsimony starting tree, then exit.",
            ),
            _Option(
                ["-z", "bipartition_filename"],
                "Name of a file containing multiple trees, e.g. from "
                "a bootstrap run, that shall be used to draw "
                "bipartition values onto a tree provided with '-t'. "
                "It can also be used to compute per-site log "
                "likelihoods in combination with '-f g', and to read "
                "a bunch of trees for a couple of other options "
                "('-f h', '-f m', '-f n').",
                filename=True,
                equate=False,
            ),
            _Option(
                ["-N", "-#", "num_replicates"],
                "Number of alternative runs on distinct starting trees. "
                "In combination with the '-b' option, this will invoke a "
                "multiple bootstrap analysis. "
                "DEFAULT: 1 single analysis."
                "Note that '-N' has been added as an alternative since "
                "'-#' sometimes caused problems with certain MPI job "
                "submission systems, since '-#' is often used to start "
                "comments. ",
                equate=False,
            ),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)
        # ENH: enforce -s, -n and -m
        if not self.parsimony_seed:
            self.parsimony_seed = 10000


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
