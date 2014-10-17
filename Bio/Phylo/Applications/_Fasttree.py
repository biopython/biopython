# Copyright 2013 by Nate Sutton.  Based on code in _Phyml.py by Eric Talevich.  All rights reserved.
# This code is part of the Biopython distribution and governed by its license.
# Please see the LICENSE file that should have been included as part of this
# package.
"""Command-line wrapper for tree inference program Fasttree."""
from __future__ import print_function

__docformat__ = "restructuredtext en"

from Bio.Application import _Option, _Switch, _Argument, AbstractCommandline


def _is_int(x):
    """Test whether the argument can be serialized as an integer."""
    return isinstance(x, int) or str(x).isdigit()


def _is_numeric(x):
    """Test whether the argument can be serialized as a number."""
    try:
        float(str(x))
        return True
    except ValueError:
        return False


class FastTreeCommandline(AbstractCommandline):
    """Command-line wrapper for FastTree.

    Homepage: http://www.microbesonline.org/fasttree/

    Citations:

    Price, M.N., Dehal, P.S., and Arkin, A.P. (2010) FastTree 2 -- Approximately
    Maximum-Likelihood Trees for Large Alignments. PLoS ONE, 5(3):e9490.
    doi:10.1371/journal.pone.0009490.

    Example usage:

    >>> import _Fasttree
    >>> fasttree_exe = r"C:\FasttreeWin32\fasttree.exe"
    >>> cmd = _Fasttree.FastTreeCommandline(fasttree_exe, input=r'C:\Input\ExampleAlignment.fsa', out='C:\Output\ExampleTree.tree')
    >>> print(cmd)
    >>> out, err = cmd()
    >>> print(out)
    >>> print(err)

    Usage advice:
    the only parameters needed are (fasttree_exe, input='<InputFile>' out='<OutputFile>')

    parameters that use values are added this way: (fasttree_exe, parameter=value, input='<InputFile>' out='<OutputFile>')
    parameters that don't use values are added this way: (fasttree_exe, parameter=True, input='<InputFile>' out='<OutputFile>')

    from the command line use 'fasttree.exe -help' or 'fasttree.exe -expert' for more explanation of usage options
    """

    def __init__(self, cmd='fasttree', **kwargs):
        self.parameters = [
            _Switch(['-nt', 'nt'],
                """By default FastTree expects protein alignments, use -nt for nucleotides""",
                ),
            _Option(['-n', 'n'],
                """-n -- read N multiple alignments in.

                    This only works with phylip interleaved format. For example, you can
                    use it with the output from phylip's seqboot. If you use -n, FastTree
                    will write 1 tree per line to standard output.""",
                checker_function=_is_int,
                equate=False,
                ),
            _Switch(['-quote', 'quote'],
                """-quote -- add quotes to sequence names in output.

                    Quote sequence names in the output and allow spaces, commas,
                    parentheses, and colons in them but not ' characters (fasta files only).""",
                ),
            _Option(['-pseudo', 'pseudo'],
                """-pseudo [weight] -- Pseudocounts are used with sequence distance estimation.

                    Use pseudocounts to estimate distances between sequences with little or no
                    overlap. (Off by default.) Recommended if analyzing the alignment has
                    sequences with little or no overlap.
                    If the weight is not specified, it is 1.0 """,
                checker_function=_is_numeric,
                equate=False,
                ),
            _Option(['-boot', 'boot'],
                """Specify the number of resamples for support values.

                    Support value options:
                    By default, FastTree computes local support values by resampling the site
                    likelihoods 1,000 times and the Shimodaira Hasegawa test. If you specify -nome,
                    it will compute minimum-evolution bootstrap supports instead
                    In either case, the support values are proportions ranging from 0 to 1

                    Use -nosupport to turn off support values or -boot 100 to use just 100 resamples.""",
                checker_function=_is_int,
                equate=False,
                ),
            _Switch(['-nosupport', 'nosupport'],
                """Turn off support values.

                    Support value options:
                    By default, FastTree computes local support values by resampling the site
                    likelihoods 1,000 times and the Shimodaira Hasegawa test. If you specify -nome,
                    it will compute minimum-evolution bootstrap supports instead
                    In either case, the support values are proportions ranging from 0 to 1

                    Use -nosupport to turn off support values or -boot 100 to use just 100 resamples.""",
                ),
            _Option(['-intree', 'intree'],
                """-intree newickfile -- read the starting tree in from newickfile.

                    Any branch lengths in the starting trees are ignored.
                    -intree with -n will read a separate starting tree for each alignment.""",
                filename=True,
                equate=False,
                ),
            _Option(['-intree1', 'intree1'],
                """-intree1 newickfile -- read the same starting tree for each alignment.""",
                filename=True,
                equate=False,
                ),
            _Switch(['-quiet', 'quiet'],
                """-quiet -- do not write to standard error during normal operation

                    (no progress indicator, no options summary, no likelihood values, etc.)""",
                ),
            _Switch(['-nopr', 'nopr'],
                """-nopr -- do not write the progress indicator to stderr.""",
                ),
            _Option(['-nni', 'nni'],
                """Set the rounds of minimum-evolution nearest-neighbor interchanges

                    Topology refinement:
                    By default, FastTree tries to improve the tree with up to 4*log2(N)
                    rounds of minimum-evolution nearest-neighbor interchanges (NNI),
                    where N is the number of unique sequences, 2 rounds of
                    subtree-prune-regraft (SPR) moves (also min. evo.), and
                    up to 2*log(N) rounds of maximum-likelihood NNIs.
                    Use -nni to set the number of rounds of min. evo. NNIs.""",
                checker_function=_is_int,
                equate=False,
                ),
            _Option(['-spr', 'spr'],
                """Set the rounds of subtree-prune-regraft moves

                    Topology refinement:
                    By default, FastTree tries to improve the tree with up to 4*log2(N)
                    rounds of minimum-evolution nearest-neighbor interchanges (NNI),
                    where N is the number of unique sequences, 2 rounds of
                    subtree-prune-regraft (SPR) moves (also min. evo.), and
                    up to 2*log(N) rounds of maximum-likelihood NNIs.
                    Use -nni to set the number of rounds of min. evo. NNIs,
                    and -spr to set the rounds of SPRs.""",
                checker_function=_is_int,
                equate=False,
                ),
            _Switch(['-noml', 'noml'],
                """Deactivate min-evo NNIs and SPRs.

                    Topology refinement:
                    By default, FastTree tries to improve the tree with up to 4*log2(N)
                    rounds of minimum-evolution nearest-neighbor interchanges (NNI),
                    where N is the number of unique sequences, 2 rounds of
                    subtree-prune-regraft (SPR) moves (also min. evo.), and
                    up to 2*log(N) rounds of maximum-likelihood NNIs.
                    Use -nni to set the number of rounds of min. evo. NNIs,
                    and -spr to set the rounds of SPRs.
                    Use -noml to turn off both min-evo NNIs and SPRs (useful if refining
                    an approximately maximum-likelihood tree with further NNIs) """,
                ),
            _Switch(['-mllen', 'mllen'],
                """Optimize branch lengths on a fixed topology.

                    Topology refinement:
                    By default, FastTree tries to improve the tree with up to 4*log2(N)
                    rounds of minimum-evolution nearest-neighbor interchanges (NNI),
                    where N is the number of unique sequences, 2 rounds of
                    subtree-prune-regraft (SPR) moves (also min. evo.), and
                    up to 2*log(N) rounds of maximum-likelihood NNIs.
                    Use -nni to set the number of rounds of min. evo. NNIs,
                    and -spr to set the rounds of SPRs.
                    Use -mllen to optimize branch lengths without ML NNIs
                    Use -mllen -nome with -intree to optimize branch lengths on a fixed topology.""",
                ),
            _Switch(['-nome', 'nome'],
                """Changes support values calculation to a minimum-evolution bootstrap method.

                    Topology refinement:
                    By default, FastTree tries to improve the tree with up to 4*log2(N)
                    rounds of minimum-evolution nearest-neighbor interchanges (NNI),
                    where N is the number of unique sequences, 2 rounds of
                    subtree-prune-regraft (SPR) moves (also min. evo.), and
                    up to 2*log(N) rounds of maximum-likelihood NNIs.
                    Use -nni to set the number of rounds of min. evo. NNIs,
                    and -spr to set the rounds of SPRs.
                    Use -mllen to optimize branch lengths without ML NNIs
                    Use -mllen -nome with -intree to optimize branch lengths on a fixed topology

                    Support value options:
                    By default, FastTree computes local support values by resampling the site
                    likelihoods 1,000 times and the Shimodaira Hasegawa test. If you specify -nome,
                    it will compute minimum-evolution bootstrap supports instead
                    In either case, the support values are proportions ranging from 0 to 1.""",
                ),
            _Option(['-mlnni', 'mlnni'],
                """Set the number of rounds of maximum-likelihood NNIs.

                    Topology refinement:
                    By default, FastTree tries to improve the tree with up to 4*log2(N)
                    rounds of minimum-evolution nearest-neighbor interchanges (NNI),
                    where N is the number of unique sequences, 2 rounds of
                    subtree-prune-regraft (SPR) moves (also min. evo.), and
                    up to 2*log(N) rounds of maximum-likelihood NNIs.
                    Use -nni to set the number of rounds of min. evo. NNIs,
                    and -spr to set the rounds of SPRs.
                    Use -mlnni to set the number of rounds of maximum-likelihood NNIs.""",
                checker_function=_is_int,
                equate=False,
                ),
            _Option(['-mlacc', 'mlacc'],
                """Option for optimization of branches at each NNI.

                    Topology refinement:
                    By default, FastTree tries to improve the tree with up to 4*log2(N)
                    rounds of minimum-evolution nearest-neighbor interchanges (NNI),
                    where N is the number of unique sequences, 2 rounds of
                    subtree-prune-regraft (SPR) moves (also min. evo.), and
                    up to 2*log(N) rounds of maximum-likelihood NNIs.
                    Use -nni to set the number of rounds of min. evo. NNIs,
                    and -spr to set the rounds of SPRs.
                    Use -mlacc 2 or -mlacc 3 to always optimize all 5 branches at each NNI,
                    and to optimize all 5 branches in 2 or 3 rounds.""",
                checker_function=_is_int,
                equate=False,
                ),
            _Switch(['-slownni', 'slownni'],
                """Turn off heuristics to avoid constant subtrees with NNIs.

                    Topology refinement:
                    By default, FastTree tries to improve the tree with up to 4*log2(N)
                    rounds of minimum-evolution nearest-neighbor interchanges (NNI),
                    where N is the number of unique sequences, 2 rounds of
                    subtree-prune-regraft (SPR) moves (also min. evo.), and
                    up to 2*log(N) rounds of maximum-likelihood NNIs.
                    Use -nni to set the number of rounds of min. evo. NNIs,
                    and -spr to set the rounds of SPRs.
                    Use -slownni to turn off heuristics to avoid constant subtrees (affects both
                    ML and ME NNIs).""",
                ),
            _Switch(['-wag', 'wag'],
                """Maximum likelihood model options: Whelan-And-Goldman 2001 model instead of (default) Jones-Taylor-Thorton 1992 model (a.a. only)""",
                ),
            _Switch(['-gtr', 'gtr'],
                """Maximum likelihood model options: Use generalized time-reversible instead of (default) Jukes-Cantor (nt only)""",
                ),
            _Option(['-cat', 'cat'],
                """Maximum likelihood model options: Specify the number of rate categories of sites (default 20).""",
                checker_function=_is_int,
                equate=False,
                ),
            _Switch(['-nocat', 'nocat'],
                """Maximum likelihood model options: No CAT model (just 1 category)""",
                ),
            _Switch(['-gamma', 'gamma'],
                """Report the likelihood under the discrete gamma model.

                    Maximum likelihood model options:
                    -gamma -- after the final round of optimizing branch lengths with the CAT model,
                    report the likelihood under the discrete gamma model with the same
                    number of categories. FastTree uses the same branch lengths but
                    optimizes the gamma shape parameter and the scale of the lengths.
                    The final tree will have rescaled lengths. Used with -log, this
                    also generates per-site likelihoods for use with CONSEL, see
                    GammaLogToPaup.pl and documentation on the FastTree web site.""",
                ),
            _Switch(['-slow', 'slow'],
                """Use an exhaustive search.

                    Searching for the best join:
                    By default, FastTree combines the 'visible set' of fast neighbor-joining with
                    local hill-climbing as in relaxed neighbor-joining
                    -slow -- exhaustive search (like NJ or BIONJ, but different gap handling)
                    -slow takes half an hour instead of 8 seconds for 1,250 proteins""",
                ),
            _Switch(['-fastest', 'fastest'],
                """Search the visible set (the top hit for each node) only.

                    Searching for the best join:
                    By default, FastTree combines the 'visible set' of fast neighbor-joining with
                    local hill-climbing as in relaxed neighbor-joining
                    -fastest -- search the visible set (the top hit for each node) only
                    Unlike the original fast neighbor-joining, -fastest updates visible(C)
                    after joining A and B if join(AB,C) is better than join(C,visible(C))
                    -fastest also updates out-distances in a very lazy way,
                    -fastest sets -2nd on as well, use -fastest -no2nd to avoid this""",
                ),
            _Switch(['-2nd', 'second'],
                """Turn 2nd-level top hits heuristic on.

                    Top-hit heuristics:
                    By default, FastTree uses a top-hit list to speed up search
                    Use -notop (or -slow) to turn this feature off
                    and compare all leaves to each other,
                    and all new joined nodes to each other

                    -2nd or -no2nd to turn 2nd-level top hits heuristic on or off
                    This reduces memory usage and running time but may lead to
                    marginal reductions in tree quality.
                    (By default, -fastest turns on -2nd.)""",
                ),
            _Switch(['-no2nd', 'no2nd'],
                """Turn 2nd-level top hits heuristic off.

                    Top-hit heuristics:
                    By default, FastTree uses a top-hit list to speed up search
                    Use -notop (or -slow) to turn this feature off
                    and compare all leaves to each other,
                    and all new joined nodes to each other

                    -2nd or -no2nd to turn 2nd-level top hits heuristic on or off
                    This reduces memory usage and running time but may lead to
                    marginal reductions in tree quality.
                    (By default, -fastest turns on -2nd.)""",
                ),
            _Option(['-seed', 'seed'],
                """Use -seed to initialize the random number generator.

                    Support value options:
                    By default, FastTree computes local support values by resampling the site
                    likelihoods 1,000 times and the Shimodaira Hasegawa test. If you specify -nome,
                    it will compute minimum-evolution bootstrap supports instead
                    In either case, the support values are proportions ranging from 0 to 1""",
                checker_function=_is_int,
                equate=False,
                ),
            _Switch(['-top', 'top'],
                """Top-hit list to speed up search

                    Top-hit heuristics:
                    By default, FastTree uses a top-hit list to speed up search
                    Use -notop (or -slow) to turn this feature off
                    and compare all leaves to each other,
                    and all new joined nodes to each other""",
                ),
            _Switch(['-notop', 'notop'],
                """Turn off top-hit list to speed up search

                    Top-hit heuristics:
                    By default, FastTree uses a top-hit list to speed up search
                    Use -notop (or -slow) to turn this feature off
                    and compare all leaves to each other,
                    and all new joined nodes to each other""",
                ),
            _Option(['-topm', 'topm'],
                """Change the top hits calculation method

                    Top-hit heuristics:
                    By default, FastTree uses a top-hit list to speed up search
                    -topm 1.0 -- set the top-hit list size to parameter*sqrt(N)
                    FastTree estimates the top m hits of a leaf from the
                    top 2*m hits of a 'close' neighbor, where close is
                    defined as d(seed,close) < 0.75 * d(seed, hit of rank 2*m),
                    and updates the top-hits as joins proceed""",
                checker_function=_is_numeric,
                equate=False,
                ),
            _Option(['-close', 'close'],
                """Modify the close heuristic for the top-hit list

                    Top-hit heuristics:
                    By default, FastTree uses a top-hit list to speed up search
                    -close 0.75 -- modify the close heuristic, lower is more conservative""",
                checker_function=_is_numeric,
                equate=False,
                ),
            _Option(['-refresh', 'refresh'],
                """Parameter for conditions that joined nodes are compared to other nodes

                    Top-hit heuristics:
                    By default, FastTree uses a top-hit list to speed up search
                    -refresh 0.8 -- compare a joined node to all other nodes if its
                    top-hit list is less than 80% of the desired length,
                    or if the age of the top-hit list is log2(m) or greater""",
                checker_function=_is_numeric,
                equate=False,
                ),
            _Option(['-matrix', 'matrix'],
                """Specify a matrix for nucleotide or amino acid distances

                    Distances:
                    Default: For protein sequences, log-corrected distances and an
                    amino acid dissimilarity matrix derived from BLOSUM45
                    or for nucleotide sequences, Jukes-Cantor distances
                    To specify a different matrix, use -matrix FilePrefix or -nomatrix""",
                filename=True,
                equate=False,
                ),
            _Switch(['-nomatrix', 'nomatrix'],
                """Specify that no matrix should be used for nucleotide or amino acid distances

                    Distances:
                    Default: For protein sequences, log-corrected distances and an
                    amino acid dissimilarity matrix derived from BLOSUM45
                    or for nucleotide sequences, Jukes-Cantor distances
                    To specify a different matrix, use -matrix FilePrefix or -nomatrix""",
                ),
            _Switch(['-nj', 'nj'],
                """Join options: regular (unweighted) neighbor-joining (default)""",
                ),
            _Switch(['-bionj', 'bionj'],
                """Join options: weighted joins as in BIONJ.  FastTree will also weight joins during NNIs""",
                ),
            _Option(['-gtrrates', 'gtrrates'],
                """-gtrrates ac ag at cg ct gt""",
                equate=False,
                ),
            _Option(['-gtrfreq', 'gtrfreq'],
                """-gtrfreq A C G T""",
                equate=False,
                ),
            _Option(['-constraints', 'constraints'],
                """Specifies an alignment file for use with constrained topology searching

                    Constrained topology search options:
                    -constraints alignmentfile -- an alignment with values of 0, 1, and -
                    Not all sequences need be present. A column of 0s and 1s defines a
                    constrained split. Some constraints may be violated
                    (see 'violating constraints:' in standard error).""",
                filename=True,
                equate=False,
                ),
            _Option(['-constraintWeight', 'constraintWeight'],
                """Weight strength of contraints in topology searching

                    Constrained topology search options:
                    -constraintWeight -- how strongly to weight the constraints. A value of 1
                    means a penalty of 1 in tree length for violating a constraint
                    Default: 100.0""",
                checker_function=_is_numeric,
                equate=False,
                ),
            _Option(['-log', 'log'],
                """Create log files of data such as intermediate trees and per-site rates

                    -log logfile -- save intermediate trees so you can extract
                    the trees and restart long-running jobs if they crash
                    -log also reports the per-site rates (1 means slowest category).""",
                filename=True,
                equate=False,
                ),
            _Option(['-makematrix', 'makematrix'],
                """-makematrix [alignment]""",
                filename=True,
                equate=False,
                ),
            _Switch(['-rawdist', 'rawdist'],
                """Use -rawdist to turn the log-correction off or to use %different instead of Jukes-Cantor in AA or NT distances

                    Distances:
                    Default: For protein sequences, log-corrected distances and an
                    amino acid dissimilarity matrix derived from BLOSUM45
                    or for nucleotide sequences, Jukes-Cantor distances
                    To specify a different matrix, use -matrix FilePrefix or -nomatrix""",
                ),
            _Option(['-sprlength', 'sprlength'],
                """Use -sprlength set the maximum length of a SPR move (default 10) in topology refinement

                    Topology refinement:
                    By default, FastTree tries to improve the tree with up to 4*log2(N)
                    rounds of minimum-evolution nearest-neighbor interchanges (NNI),
                    where N is the number of unique sequences, 2 rounds of
                    subtree-prune-regraft (SPR) moves (also min. evo.), and
                    up to 2*log(N) rounds of maximum-likelihood NNIs.
                    Use -nni to set the number of rounds of min. evo. NNIs,
                    and -spr to set the rounds of SPRs.""",
                checker_function=_is_int,
                equate=False,
                ),
             _Switch(['-help', 'help'],
                """Show the help"""
                ),
             _Switch(['-expert', 'expert'],
                """Show the expert level help"""
                ),
             _Option(['-out', 'out'],
                """Enter <output file>

                    The path to a Newick Tree output file needs to be specified.""",
                filename=True,
                equate=False,
                ),
             _Argument(['input'],
                """Enter <input file>

                    An input file of sequence alignments in fasta or phylip format is needed.  By default FastTree expects protein
                    alignments, use -nt for nucleotides""",
                filename=True,
                is_required=True,
                ),
                ]

        AbstractCommandline.__init__(self, cmd, **kwargs)
