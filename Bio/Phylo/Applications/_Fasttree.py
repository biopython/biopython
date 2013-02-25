# Copyright 2013 by Nate Sutton.  Based on code in _Phyml.py by Eric Talevich.  All rights reserved.
# This code is part of the Biopython distribution and governed by its license.
# Please see the LICENSE file that should have been included as part of this
# package.
"""Command-line wrapper for tree inference program Fasttree."""
__docformat__ = "restructuredtext en"

from Bio.Application import _Option, _Switch, AbstractCommandline


class FastTreeCommandline(AbstractCommandline):
    """Command-line wrapper for FastTree.

    Homepage: http://www.microbesonline.org/fasttree/

    Citations:

    Price, M.N., Dehal, P.S., and Arkin, A.P. (2010) FastTree 2 -- Approximately 
    Maximum-Likelihood Trees for Large Alignments. PLoS ONE, 5(3):e9490. 
    doi:10.1371/journal.pone.0009490.
    
    Example usage:
    ######
    import _Fasttree

    fasttree_exe = "C:\\FasttreeWin32\\fasttree.exe"
    cmd = _Fasttree.FastTreeCommandline(fasttree_exe, output_input='C:\\Output\\ExampleTree.tree C:\\Input\\ExampleAlignment.fsa')
    print cmd
    out, err = cmd()
    print out
    print err
    ######
    
    Usage advice:
    the only parameters needed are (fasttree_exe, output_input='<OutputFile> <InputFile>')
    
    parameters that use values are added this way: (fasttree_exe, parameter=value, output_input='<OutputFile> <InputFile>')
    parameters that don't use values are added this way: (fasttree_exe, parameter=True, output_input='<OutputFile> <InputFile>')
    
    from the command line use 'fasttree.exe -help' or 'fasttree.exe -expert' for more explanation of usage options
    """

    def __init__(self, cmd='fasttree', **kwargs):
        self.parameters = [                                                     
            _Switch(['-nt', 'nt'],
                """By default FastTree expects protein alignments,  use -nt for nucleotides""",
                ),
            _Option(['-n', 'n'],
                """-n -- read in multiple alignments in. This only
                    works with phylip interleaved format. For example, you can
                    use it with the output from phylip's seqboot. If you use -n, FastTree
                    will write 1 tree per line to standard output.""",
                equate=False,
                ),                               
            _Switch(['-quote', 'quote'],
                """-quote -- quote sequence names in the output and allow spaces, commas,
                    parentheses, and colons in them but not ' characters (fasta files only) """,           
                ),
            _Option(['-pseudo', 'pseudo'],
                """-pseudo [weight] -- Use pseudocounts to estimate distances between
                    sequences with little or no overlap. (Off by default.) Recommended
                    if analyzing the alignment has sequences with little or no overlap.
                    If the weight is not specified, it is 1.0 """,                 
                equate=False,                          
                ),
            _Option(['-boot', 'boot'],
                """Support value options:
                    By default, FastTree computes local support values by resampling the site
                    likelihoods 1,000 times and the Shimodaira Hasegawa test. If you specify -nome
                    ,
                    it will compute minimum-evolution bootstrap supports instead
                    In either case, the support values are proportions ranging from 0 to 1

                    Use -nosupport to turn off support values or -boot 100 to use just 100 resamples """,
                equate=False,                    
                ),
            _Switch(['-nosupport', 'nosupport'],
                """Support value options:
                    By default, FastTree computes local support values by resampling the site
                    likelihoods 1,000 times and the Shimodaira Hasegawa test. If you specify -nome
                    ,
                    it will compute minimum-evolution bootstrap supports instead
                    In either case, the support values are proportions ranging from 0 to 1

                    Use -nosupport to turn off support values or -boot 100 to use just 100 resamples """,           
                ),         
            _Option(['-intree', 'intree'],
                """-intree newickfile -- read the starting tree in from newickfile.
                    Any branch lengths in the starting trees are ignored.
                    -intree with -n will read a separate starting tree for each alignment. """,
                equate=False,                    
                ),                  
            _Option(['-intree1', 'intree1'],
                """-intree1 newickfile -- read the same starting tree for each alignment """,
                equate=False,                
                ),  
            _Switch(['-quiet', 'quiet'],
                """-quiet -- do not write to standard error during normal operation (no progress
                    indicator, no options summary, no likelihood values, etc.) """,      
                ),  
            _Switch(['-nopr', 'nopr'],
                """-nopr -- do not write the progress indicator to stderr """, 
                ),  
            _Option(['-nni', 'nni'],
                """Topology refinement:
                    By default, FastTree tries to improve the tree with up to 4*log2(N)
                    rounds of minimum-evolution nearest-neighbor interchanges (NNI),
                    where N is the number of unique sequences, 2 rounds of
                    subtree-prune-regraft (SPR) moves (also min. evo.), and
                    up to 2*log(N) rounds of maximum-likelihood NNIs.
                    Use -nni to set the number of rounds of min. evo. NNIs """,
                equate=False,                    
                ),               
            _Option(['-spr', 'spr'],
                """Topology refinement:
                    By default, FastTree tries to improve the tree with up to 4*log2(N)
                    rounds of minimum-evolution nearest-neighbor interchanges (NNI),
                    where N is the number of unique sequences, 2 rounds of
                    subtree-prune-regraft (SPR) moves (also min. evo.), and
                    up to 2*log(N) rounds of maximum-likelihood NNIs.
                    Use -nni to set the number of rounds of min. evo. NNIs,
                    and -spr to set the rounds of SPRs. """,
                equate=False,                    
                ),   
            _Switch(['-noml', 'noml'],
                """Topology refinement:
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
                """Topology refinement:
                    By default, FastTree tries to improve the tree with up to 4*log2(N)
                    rounds of minimum-evolution nearest-neighbor interchanges (NNI),
                    where N is the number of unique sequences, 2 rounds of
                    subtree-prune-regraft (SPR) moves (also min. evo.), and
                    up to 2*log(N) rounds of maximum-likelihood NNIs.
                    Use -nni to set the number of rounds of min. evo. NNIs,
                    and -spr to set the rounds of SPRs.
                    Use -mllen to optimize branch lengths without ML NNIs
                    Use -mllen -nome with -intree to optimize branch lengths on a fixed topology """,     
                ),    
            _Switch(['-nome', 'nome'],
                ###########################################################
                # Not sure if this should be a _Switch or _Option instead #  
                ###########################################################
                """Topology refinement:
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
                    In either case, the support values are proportions ranging from 0 to 1""",        
                ),
            _Option(['-mlnni', 'mlnni'],
                """Topology refinement:
                    By default, FastTree tries to improve the tree with up to 4*log2(N)
                    rounds of minimum-evolution nearest-neighbor interchanges (NNI),
                    where N is the number of unique sequences, 2 rounds of
                    subtree-prune-regraft (SPR) moves (also min. evo.), and
                    up to 2*log(N) rounds of maximum-likelihood NNIs.
                    Use -nni to set the number of rounds of min. evo. NNIs,
                    and -spr to set the rounds of SPRs.
                    Use -mlnni to set the number of rounds of maximum-likelihood NNIs """,
                equate=False,                    
                ),          
            _Option(['-mlacc', 'mlacc'],
                """Topology refinement:
                    By default, FastTree tries to improve the tree with up to 4*log2(N)
                    rounds of minimum-evolution nearest-neighbor interchanges (NNI),
                    where N is the number of unique sequences, 2 rounds of
                    subtree-prune-regraft (SPR) moves (also min. evo.), and
                    up to 2*log(N) rounds of maximum-likelihood NNIs.
                    Use -nni to set the number of rounds of min. evo. NNIs,
                    and -spr to set the rounds of SPRs.
                    Use -mlacc 2 or -mlacc 3 to always optimize all 5 branches at each NNI,
                    and to optimize all 5 branches in 2 or 3 rounds """,
                equate=False,                    
                ),  
            _Switch(['-slownni', 'slownni'],
                """Topology refinement:
                    By default, FastTree tries to improve the tree with up to 4*log2(N)
                    rounds of minimum-evolution nearest-neighbor interchanges (NNI),
                    where N is the number of unique sequences, 2 rounds of
                    subtree-prune-regraft (SPR) moves (also min. evo.), and
                    up to 2*log(N) rounds of maximum-likelihood NNIs.
                    Use -nni to set the number of rounds of min. evo. NNIs,
                    and -spr to set the rounds of SPRs.
                    Use -slownni to turn off heuristics to avoid constant subtrees (affects both
                    ML and ME NNIs) """,  
                ),            
            _Switch(['-wag', 'wag'],
                """Maximum likelihood model options:
                    -wag -- Whelan-And-Goldman 2001 model instead of (default) Jones-Taylor-Thorton
                    1992 model (a.a. only)""",
                ),
            _Switch(['-gtr', 'gtr'],
                """Maximum likelihood model options:
                    -gtr -- generalized time-reversible instead of (default) Jukes-Cantor (nt only)""",
                ),                      
            _Option(['-cat', 'cat'],
                """Maximum likelihood model options:
                    -cat # -- specify the number of rate categories of sites (default 20) """,
                equate=False,                    
                ),      
            _Switch(['-nocat', 'nocat'],
                """Maximum likelihood model options:
                    -nocat -- no CAT model (just 1 category)""",   
                ),
            _Switch(['-gamma', 'gamma'],
                """Maximum likelihood model options:
                    -gamma -- after the final round of optimizing branch lengths with the CAT model,
                    report the likelihood under the discrete gamma model with the same
                    number of categories. FastTree uses the same branch lengths but
                    optimizes the gamma shape parameter and the scale of the lengths.
                    The final tree will have rescaled lengths. Used with -log, this
                    also generates per-site likelihoods for use with CONSEL, see
                    GammaLogToPaup.pl and documentation on the FastTree web site.""", 
                ),
            _Switch(['-slow', 'slow'],
                """Searching for the best join:
                    By default, FastTree combines the 'visible set' of fast neighbor-joining with
                    local hill-climbing as in relaxed neighbor-joining
                    -slow -- exhaustive search (like NJ or BIONJ, but different gap handling)
                    -slow takes half an hour instead of 8 seconds for 1,250 proteins""",      
                ),
            _Switch(['-fastest', 'fastest'],
                """Searching for the best join:
                    By default, FastTree combines the 'visible set' of fast neighbor-joining with
                    local hill-climbing as in relaxed neighbor-joining
                    -fastest -- search the visible set (the top hit for each node) only
                    Unlike the original fast neighbor-joining, -fastest updates visible(C)
                    after joining A and B if join(AB,C) is better than join(C,visible(C))
                    -fastest also updates out-distances in a very lazy way,
                    -fastest sets -2nd on as well, use -fastest -no2nd to avoid this""",   
                ),     
            _Switch(['-2nd', 'Second'],
                """Top-hit heuristics:
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
                """Top-hit heuristics:
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
                """Support value options:
                    By default, FastTree computes local support values by resampling the site
                    likelihoods 1,000 times and the Shimodaira Hasegawa test. If you specify -nome,
                    it will compute minimum-evolution bootstrap supports instead
                    In either case, the support values are proportions ranging from 0 to 1

                    Use -seed to initialize the random number generator""",
                equate=False,                    
                ),      
            _Switch(['-top', 'top'],
                """Top-hit heuristics:
                    By default, FastTree uses a top-hit list to speed up search
                    Use -notop (or -slow) to turn this feature off
                    and compare all leaves to each other,
                    and all new joined nodes to each other""",           
                ),
            _Switch(['-notop', 'notop'],
                """Top-hit heuristics:
                    By default, FastTree uses a top-hit list to speed up search
                    Use -notop (or -slow) to turn this feature off
                    and compare all leaves to each other,
                    and all new joined nodes to each other""",           
                ),  
            _Option(['-topm', 'topm'],
                """Top-hit heuristics:
                    By default, FastTree uses a top-hit list to speed up search
                    -topm 1.0 -- set the top-hit list size to parameter*sqrt(N)
                    FastTree estimates the top m hits of a leaf from the
                    top 2*m hits of a 'close' neighbor, where close is
                    defined as d(seed,close) < 0.75 * d(seed, hit of rank 2*m),
                    and updates the top-hits as joins proceed""",
                equate=False,                    
                ),
            _Option(['-close', 'close'],
                """Top-hit heuristics:
                    By default, FastTree uses a top-hit list to speed up search
                    -close 0.75 -- modify the close heuristic, lower is more conservative""",
                equate=False,                    
                ),
            _Option(['-refresh', 'refresh'],
                """Top-hit heuristics:
                    By default, FastTree uses a top-hit list to speed up search
                    -refresh 0.8 -- compare a joined node to all other nodes if its
                    top-hit list is less than 80% of the desired length,
                    or if the age of the top-hit list is log2(m) or greater""",
                equate=False,                    
                ),
            _Option(['-matrix', 'matrix'],
                """Distances:
                    Default: For protein sequences, log-corrected distances and an
                    amino acid dissimilarity matrix derived from BLOSUM45
                    or for nucleotide sequences, Jukes-Cantor distances
                    To specify a different matrix, use -matrix FilePrefix or -nomatrix""",
                equate=False,                    
                ),
            _Switch(['-nomatrix', 'nomatrix'],
                """Distances:
                    Default: For protein sequences, log-corrected distances and an
                    amino acid dissimilarity matrix derived from BLOSUM45
                    or for nucleotide sequences, Jukes-Cantor distances
                    To specify a different matrix, use -matrix FilePrefix or -nomatrix""",            
                ),
            _Switch(['-nj', 'nj'],
                """Join options:
                    -nj: regular (unweighted) neighbor-joining (default)""",              
                ),
            _Switch(['-bionj', 'bionj'],
                """Join options:
                    -bionj: weighted joins as in BIONJ
                    FastTree will also weight joins during NNIs""",           
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
                """Constrained topology search options:
                    -constraints alignmentfile -- an alignment with values of 0, 1, and -
                    Not all sequences need be present. A column of 0s and 1s defines a
                    constrained split. Some constraints may be violated
                    (see 'violating constraints:' in standard error).""",
                equate=False,                    
                ),
            _Option(['-constraintWeight', 'constraintWeight'],
                """Constrained topology search options:
                    -constraintWeight -- how strongly to weight the constraints. A value of 1
                    means a penalty of 1 in tree length for violating a constraint
                    Default: 100.0""",
                equate=False,                    
                ),          
            _Option(['-log', 'log'],
                """-log logfile -- save intermediate trees so you can extract
                    the trees and restart long-running jobs if they crash
                    -log also reports the per-site rates (1 means slowest category) """,
                equate=False,                    
                ),      
            _Option(['-makematrix', 'makematrix'],
                """-makematrix [alignment]""",
                equate=False,                
                ), 
            _Switch(['-rawdist', 'rawdist'],
                """Distances:
                    Default: For protein sequences, log-corrected distances and an
                    amino acid dissimilarity matrix derived from BLOSUM45
                    or for nucleotide sequences, Jukes-Cantor distances
                    To specify a different matrix, use -matrix FilePrefix or -nomatrix
                    Use -rawdist to turn the log-correction off
                    or to use %different instead of Jukes-Cantor""",    
                ),
            _Option(['-sprlength', 'sprlength'],
                """Topology refinement:
                    By default, FastTree tries to improve the tree with up to 4*log2(N)
                    rounds of minimum-evolution nearest-neighbor interchanges (NNI),
                    where N is the number of unique sequences, 2 rounds of
                    subtree-prune-regraft (SPR) moves (also min. evo.), and
                    up to 2*log(N) rounds of maximum-likelihood NNIs.
                    Use -nni to set the number of rounds of min. evo. NNIs,
                    and -spr to set the rounds of SPRs.
                    Use -sprlength set the maximum length of a SPR move (default 10)""",
                equate=False,                    
                ),
            _Switch(['-help', 'help'],  
                """FastTree 2.1.7 No SSE3:
  FastTree protein_alignment > tree
  FastTree < protein_alignment > tree
  FastTree -out tree protein_alignment
  FastTree -nt nucleotide_alignment > tree
  FastTree -nt -gtr < nucleotide_alignment > tree
  FastTree < nucleotide_alignment > tree
FastTree accepts alignments in fasta or phylip interleaved formats

Common options (must be before the alignment file):
  -quiet to suppress reporting information
  -nopr to suppress progress indicator
  -log logfile -- save intermediate trees, settings, and model details
  -fastest -- speed up the neighbor joining phase & reduce memory usage
        (recommended for >50,000 sequences)
  -n <number> to analyze multiple alignments (phylip format only)
        (use for global bootstrap, with seqboot and CompareToBootstrap.pl)
  -nosupport to not compute support values
  -intree newick_file to set the starting tree(s)
  -intree1 newick_file to use this starting tree for all the alignments
        (for faster global bootstrap on huge alignments)
  -pseudo to use pseudocounts (recommended for highly gapped sequences)
  -gtr -- generalized time-reversible model (nucleotide alignments only)
  -wag -- Whelan-And-Goldman 2001 model (amino acid alignments only)
  -quote -- allow spaces and other restricted characters (but not ' ) in
           sequence names and quote names in the output tree (fasta input only;
           FastTree will not be able to read these trees back in)
  -noml to turn off maximum-likelihood
  -nome to turn off minimum-evolution NNIs and SPRs
        (recommended if running additional ML NNIs with -intree)
  -nome -mllen with -intree to optimize branch lengths for a fixed topology
  -cat # to specify the number of rate categories of sites (default 20)
  -gamma -- after optimizing the tree under the CAT approximation,
      rescale the lengths to optimize the Gamma20 likelihood
  -constraints constraintAlignment to constrain the topology search
       constraintAlignment should have 1s or 0s to indicates splits
  -expert -- see more options
For more information, see http://www.microbesonline.org/fasttree/ """,
                ),
            _Switch(['-expert', 'expert'],
                """Detailed usage for FastTree 2.1.7 No SSE3:
FastTree [-nt] [-n 100] [-quote] [-pseudo | -pseudo 1.0]
           [-boot 1000 | -nosupport]
           [-intree starting_trees_file | -intree1 starting_tree_file]
           [-quiet | -nopr]
           [-nni 10] [-spr 2] [-noml | -mllen | -mlnni 10]
           [-mlacc 2] [-cat 20 | -nocat] [-gamma]
           [-slow | -fastest] [-2nd | -no2nd] [-slownni] [-seed 1253]
           [-top | -notop] [-topm 1.0 [-close 0.75] [-refresh 0.8]]
           [-matrix Matrix | -nomatrix] [-nj | -bionj]
           [-wag] [-nt] [-gtr] [-gtrrates ac ag at cg ct gt] [-gtrfreq A C G T]
           [ -constraints constraintAlignment [ -constraintWeight 100.0 ] ]
           [-log logfile]
         [ alignment_file ]
        [ -out output_newick_file | > newick_tree]

or

FastTree [-nt] [-matrix Matrix | -nomatrix] [-rawdist] -makematrix [alignment]
    [-n 100] > phylip_distance_matrix

  FastTree supports fasta or phylip interleaved alignments
  By default FastTree expects protein alignments,  use -nt for nucleotides
  FastTree reads standard input if no alignment file is given

Input/output options:
  -n -- read in multiple alignments in. This only
    works with phylip interleaved format. For example, you can
    use it with the output from phylip's seqboot. If you use -n, FastTree
    will write 1 tree per line to standard output.
  -intree newickfile -- read the starting tree in from newickfile.
     Any branch lengths in the starting trees are ignored.
    -intree with -n will read a separate starting tree for each alignment.
  -intree1 newickfile -- read the same starting tree for each alignment
  -quiet -- do not write to standard error during normal operation (no progress
     indicator, no options summary, no likelihood values, etc.)
  -nopr -- do not write the progress indicator to stderr
  -log logfile -- save intermediate trees so you can extract
    the trees and restart long-running jobs if they crash
    -log also reports the per-site rates (1 means slowest category)
  -quote -- quote sequence names in the output and allow spaces, commas,
    parentheses, and colons in them but not ' characters (fasta files only)

Distances:
  Default: For protein sequences, log-corrected distances and an
     amino acid dissimilarity matrix derived from BLOSUM45
  or for nucleotide sequences, Jukes-Cantor distances
  To specify a different matrix, use -matrix FilePrefix or -nomatrix
  Use -rawdist to turn the log-correction off
  or to use %different instead of Jukes-Cantor

  -pseudo [weight] -- Use pseudocounts to estimate distances between
      sequences with little or no overlap. (Off by default.) Recommended
      if analyzing the alignment has sequences with little or no overlap.
      If the weight is not specified, it is 1.0

Topology refinement:
  By default, FastTree tries to improve the tree with up to 4*log2(N)
  rounds of minimum-evolution nearest-neighbor interchanges (NNI),
  where N is the number of unique sequences, 2 rounds of
  subtree-prune-regraft (SPR) moves (also min. evo.), and
  up to 2*log(N) rounds of maximum-likelihood NNIs.
  Use -nni to set the number of rounds of min. evo. NNIs,
  and -spr to set the rounds of SPRs.
  Use -noml to turn off both min-evo NNIs and SPRs (useful if refining
       an approximately maximum-likelihood tree with further NNIs)
  Use -sprlength set the maximum length of a SPR move (default 10)
  Use -mlnni to set the number of rounds of maximum-likelihood NNIs
  Use -mlacc 2 or -mlacc 3 to always optimize all 5 branches at each NNI,
      and to optimize all 5 branches in 2 or 3 rounds
  Use -mllen to optimize branch lengths without ML NNIs
  Use -mllen -nome with -intree to optimize branch lengths on a fixed topology
  Use -slownni to turn off heuristics to avoid constant subtrees (affects both
       ML and ME NNIs)

Maximum likelihood model options:
  -wag -- Whelan-And-Goldman 2001 model instead of (default) Jones-Taylor-Thorto
n 1992 model (a.a. only)
  -gtr -- generalized time-reversible instead of (default) Jukes-Cantor (nt only
)
  -cat # -- specify the number of rate categories of sites (default 20)
  -nocat -- no CAT model (just 1 category)
  -gamma -- after the final round of optimizing branch lengths with the CAT mode
l,
            report the likelihood under the discrete gamma model with the same
            number of categories. FastTree uses the same branch lengths but
            optimizes the gamma shape parameter and the scale of the lengths.
            The final tree will have rescaled lengths. Used with -log, this
            also generates per-site likelihoods for use with CONSEL, see
            GammaLogToPaup.pl and documentation on the FastTree web site.

Support value options:
  By default, FastTree computes local support values by resampling the site
  likelihoods 1,000 times and the Shimodaira Hasegawa test. If you specify -nome
,
  it will compute minimum-evolution bootstrap supports instead
  In either case, the support values are proportions ranging from 0 to 1

  Use -nosupport to turn off support values or -boot 100 to use just 100 resampl
es
  Use -seed to initialize the random number generator

Searching for the best join:
  By default, FastTree combines the 'visible set' of fast neighbor-joining with
      local hill-climbing as in relaxed neighbor-joining
  -slow -- exhaustive search (like NJ or BIONJ, but different gap handling)
      -slow takes half an hour instead of 8 seconds for 1,250 proteins
  -fastest -- search the visible set (the top hit for each node) only
      Unlike the original fast neighbor-joining, -fastest updates visible(C)
      after joining A and B if join(AB,C) is better than join(C,visible(C))
      -fastest also updates out-distances in a very lazy way,
      -fastest sets -2nd on as well, use -fastest -no2nd to avoid this

Top-hit heuristics:
  By default, FastTree uses a top-hit list to speed up search
  Use -notop (or -slow) to turn this feature off
         and compare all leaves to each other,
         and all new joined nodes to each other
  -topm 1.0 -- set the top-hit list size to parameter*sqrt(N)
         FastTree estimates the top m hits of a leaf from the
         top 2*m hits of a 'close' neighbor, where close is
         defined as d(seed,close) < 0.75 * d(seed, hit of rank 2*m),
         and updates the top-hits as joins proceed
  -close 0.75 -- modify the close heuristic, lower is more conservative
  -refresh 0.8 -- compare a joined node to all other nodes if its
         top-hit list is less than 80% of the desired length,
         or if the age of the top-hit list is log2(m) or greater
   -2nd or -no2nd to turn 2nd-level top hits heuristic on or off
      This reduces memory usage and running time but may lead to
      marginal reductions in tree quality.
      (By default, -fastest turns on -2nd.)

Join options:
  -nj: regular (unweighted) neighbor-joining (default)
  -bionj: weighted joins as in BIONJ
          FastTree will also weight joins during NNIs

Constrained topology search options:
  -constraints alignmentfile -- an alignment with values of 0, 1, and -
       Not all sequences need be present. A column of 0s and 1s defines a
       constrained split. Some constraints may be violated
       (see 'violating constraints:' in standard error).
  -constraintWeight -- how strongly to weight the constraints. A value of 1
       means a penalty of 1 in tree length for violating a constraint
       Default: 100.0

For more information, see http://www.microbesonline.org/fasttree/
   or the comments in the source code""",
                ),                               
             _Option(['-out', 'output_input'],
                """The path to a Newick Tree output file needs to be specified first.  After that, using a separation of a space, 
                    an input file of sequence alignments in fasta or phylip format is needed.  By default FastTree expects protein 
                    alignments,  use -nt for nucleotides""",
                is_required=True,
                equate=False,
                ),                             
                ]        

        AbstractCommandline.__init__(self, cmd, **kwargs)
