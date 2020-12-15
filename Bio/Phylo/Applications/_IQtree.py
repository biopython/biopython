
# Based on code in _Phyml.py by Eric Talevich.
# All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Command-line wrapper for the phylogenomic inference IQ-Tree"""

from Bio.Application import _Option, _Switch, _Argument, AbstractCommandline


def _is_int(x):
	"""Checker function for Integer required inputs"""
	return isinstance(x, int) or x.isdigit()

def _is_number(x):
	"""Checker function for numeric inputs that include float"""
	try:
		float(x)
		return True
	except ValueError:
		return False


class IQTreeCommandline(AbstractCommandline):
	r"""Command-line wrapper for IQTree.

	Specify command is mandatory, used to specify input alignment file.

	From the terminal command line use ``iqtree -h`` or ``iqtree -?``
	for more explanation of usage options and commands

	Homepage http://www.iqtree.org/

	References:

	Usage:

		iqtree -s <alignment> [OPTIONS]

	Example on Windows:

		import _IQtree

		iqtree_exe = r"C:\iqtree-Windows\bin\iqtree.exe"
		input = r"C:\Documents\example.phy"

		cmd = _IQtree.IQTreeCommandline(iqtree_exe, Specify = input)
		print(cmd)
		cmd()

	"""

	def __init__(self, cmd="iqtree", **kwargs):
		"""initialize the class"""

		self.parameters = [
			_Option(						#INPUT 
				["-s", "s"],   
				"""Specify input alignment file in PHYLIP, FASTA, NEXUS
				   CLUSTAL or MSF format
				   """,
				is_required = True,
				equate = False,
				filename=True
				),  
			_Option(
				["-st", "st"],
				"""Specify sequence type as either DNA, AA, BIN, MORPH, 
			       CODON or NT2AA for DNA, ammino-acid, binary, morphological
				   codon or DNA-to-AA translate sequences
				   
				   Necessary only if IQ-Tree did not detect the sequence correctly
				   
				   Note: This is always necessary when using codon models otherwise
				         IQ-tree applies DNA models
				   """,
				equate = False,
				checker_function= lambda x: x in ("DNA", "AA", "BIN", "MORPH", "CODON", "NT2AA", "DNA-to-AA")
				),
			_Option(
				["-t", "t"],
				"""Specify a file containing starting tree for tree search.
					BIONJ starts a tree search from BIONJ tree
					RANDOM starts tree search from completely random tree
					
					Default: 100 parsimony trees + BIONJ tree
					""",
				equate = False,
				filename = (lambda x: not(x in ("BIONJ", "RANDOM","PARS", "PLLPARS"))),                            #Not sure if this works
				checker_function= lambda x: True if filename else (x in ("BIONJ", "RANDOM","PARS", "PLLPARS"))     #Need tree file to test
				),
			_Option(
				["-te", "te"],
				"""Like -t but fixing user tree, no tree search is performed
				   and the program computes the log-likelihood of the fixed user tree

				   Specify a user-defined tree to determine ancestral sequences along this tree.
				   If the nodes do not have names IQTree will assign node names as Node1, Node2, etc.
				   """,
				equate = False,
				filename = (lambda x: not(x in ("BIONJ", "RANDOM","PARS", "PLLPARS"))),                            #Not sure if this works
				checker_function= lambda x: True if filename else (x in ("BIONJ", "RANDOM","PARS", "PLLPARS"))     #Need tree file to test
				),
			_Option(
				["-o", "o"],
				"""Specify an outgroup taxon name to root the tree
				   Output tree will be rooted accordingly

				   Default: first taxon in alignment
				   """,
				equate = False
				),
			_Option(
				["-pre", "pre"],                        #This outputs the files in the project folder for some reason
				"""Specify a prefix for all output files

				   Default: either alignment file name (-s) or partition file name
				   (-q, -spp or -sp)
				   """,
				equate = False
					),
			_Option(
				["-nt", "nt"],
				""" Specify the number of CPU cores for the multicore version.
					A special option -nt AUTO will tell IQTree to automatically determine
					the best number of cores given the current data and computer
					""",
				checker_function = (lambda x: isinstance(x, int) or x.isdigit() or x in ("AUTO")),
				equate = False
				),
			_Option(
				["-ntmax", "ntmax"],
				"""Specify the maximal number of CPU cores.
					Default: # CPU cores on the current machine
					""",
				checker_function = (lambda x: isinstance(x, int) or x.isdigit() or x in ("AUTO")),
				equate= False
				),
			_Switch(
				["-v", "verbose"],
				"""Turn on verbose mode for printing more messages to screen, used for debugging
					Default: OFF
					"""
				),
			_Switch(
				["-quiet", "quiet"],
				"""Silent mode, suppress printing to the screen, note that .log file is still 
					written
					"""
				),
			_Switch(
				["-keep-ident", "keepident"],
				"""Keep identical sequences in the alignment.
				   Default: IQTree will remove identical sequences during the analysis and add
				   them at the end
				   """
				),
			_Switch(
				["-safe","safe"],
				"""Turn on safe numerical mode to avoid numerical underflow for large data sets
				   with many sequences.
				   This mode is automatically turned on when having more than 2000 sequences
				   """
				),
			_Option(
				["-mem", "mem"],
				"""Specify maximal RAM usage 
				   Example:
						-mem 64G to use at most 64 GB of RAM
						-mem 200M to use at most 200 MB of RAM
						-and so on

				   Default: IQTree will not exceed the computer RAM size
				   """,
				equate= False
				),
			_Switch(
				["-redo", "redo"],
				"""Redo the entire analysis no matter if it was stopped or successful

				   Note: This option will overwrite all existing output files
				   """
				),
			_Option(
				["-cptime", "cptime"],
				"""Specify the minimum checkpoint time interval in seconds

				   Default 20 seconds
				   """,
				equate = False,
				checker_function = _is_int
				),

			### Likelihood mapping analysis comands ###

			_Option(
				["-lmap", "lmap"],
				"""Specify the number of quartets to be randomly drawn.
				   If you specify -lmap ALL, all unique quartets will be drawn, instead.
				   """,
				equate = False,
				checker_function = lambda x: isinstance(x, int) or x.isdigit() or x in ("ALL")
				),
			_Option(
				["-lmclust", "lmclust"],
				"""Specify a NEXUS file containing taxon clusters for quartet mapping analysis
				""",
				filename = True,
				equate = False
				),
			_Switch(
				["-wql", "wql"],
				"""Write quartet log-likelihoods into '.lmap.quartelh' file
				"""
				),
			_Option(
				["-n", "n"],
				"""Fix number of iterations to stop (this option overrides -nstop criterion)
				   
				   0 - To skip subsequent tree search, useful when you want to assess the
				       phylogenetic information of the alignment

				   Default: auto
				   """,

				equate = False,
				checker_function = _is_int
				),

			### Automatic model selection ###

			_Option(
				["-m", "m"],
				"""The default model may not fit well to the data, therefore IQTree allows to
				   automatically determine the best-fit model via a series of -m TEST options:
				   - TESTONLY        (Standard model selection)
				   - TEST            (Like TESTONLY but immediately followed by tree reconstruction using
				                      the best-fit model)
				   - TESTNEWONLY     (or MF, Perform an extended model selection that additionally includes
				                      FreeRate model compared with -m TESTONLY)
				   - TESTNEW         (or MFP, like -m MF but immediately followed by tree reconstruction 
					                  using the best-fit model found)

				  IQTree version 1.6 or later allows to additionally test Lie Markov DNA models:
				   - LM              (Additionally consider all Lie Markov models
				   - LMRY            (Additionally consider all Lie Markov models with RY symmetry)
				   - LMWS            (Additionally consider all Lie Markov models with WS symmetry)
				   - LMMK            (Additionally consider all Lie Markov models with MK symmetry)
				   - LMSS            (Additionally consider all strand-symmetric Lie Markov models)
				   
				  When a partition file is specified (-q, -spp, -sp) you can append MERGE keyword into -m
				  option to find the best-fit partitioning scheme like PartitionFinder:
				  - TESTMERGEONLY    (Select best-fit partitioning scheme by possibly merging partitions to
									  reduce over-parameterization and increase model fit)
				  - TESTMERGE	     (Like TESTMERGEONLY but immediately followed by tree reconstruction using
									  the best partitioning scheme found)
				  - TESTNEWMERGEONLY (Like TESTMERGEONLY but additionally includes FreeRate model)
				  - TESTNEWMERGE	 (Like MF+MERGE but immediately followed by tree reconstruction using the
									  best partitioning scheme found)
									  
				  -m is a powerful tool option to specify substitution models, state frequency and rate
				  heterogeneity type, the general syntax is:

				     -m MODEL+FreqType+RateType

				  where MODEL is a model name, +FreqType (Optional) is the frequency type and +RateType (Optional)
				  is the rate heterogeneity type
				  (For more information on supported Model names or Frequency Typings available consult IQTree
				   documentation at http://www.iqtree.org/doc/Command-Reference)
				   """,
				equate = False,
				),   
			_Option(
				["-rcluster", "rcluster"],
				"""Specify the percentage for the relaxed clustering algorithm to speed up the computation
				   instead of the default slow greedy algorithm.

				   Example: '-rcluster 10' only the top 10% partition schemes are considered to save 
							computations
							""",
				equate = False,
				checker_function = _is_int
				),
			_Option(
				["-rclusterf", "rclusterf"],
				"""Similar to -rcluster but using the fast relaxed clustering algorithm of PartitionFinder2
				""",
				equate = False,
				checker_function= _is_int
				),
			_Option(
				["-rcluster-max", "rclustermax"],
				"""Specify the absolute maximum number of partition pairs in the partition merging phase.

				   Default: the larger of 1000 and 10 times the number of partitions
				   """,
				equate = False,
				checker_function  =_is_int
				),
			_Option(
				["-mset", "mset"],
				"""Specify the name of a program to restrict to only those models supported by the specified
				   program.
				   Alternatively, one can specify a comma-separated list of base models

				   Example: -mset WAG,LG,JTT will restrict model selection to WAG, LG and JTT instead of
				            all 18 AA models to save computations
				   """,
				equate = False
				),
			_Option(
				["-msub", "msub"],
				"""Specify either nuclear, mitochondrial, chloroplast or viral to restrict those AA models
				   designed for specified source
				   """,
				equate = False,
				checker_function = lambda x: x in ("nuclear", "mitochondrial", "chloroplast", "viral")
				),
			_Option(
				["-mfreq", "mfreq"],
				"""Specify a comma-separated list of frequency types for model selection.
			   
				   Default: -mfreq FU,F for protein models
				            -mfreq ,F1x4,F3x4,F for codon models
				   """,
				equate = False
				),
			_Option(
				["-mrate", "mrate"],
				"""Specify a comma-separated list of rate heterogeneity types for model selection

				   Default: -mrate E,I,G,I+G for standard procedure
				            -mrate E,I,G,I+G,R for new selection procedure
				   """,
				equate = False,
				),
			_Option(
				["-cmin", "cmin"],
				"""Specify minimum number of categories of FreeRate model

				   Default: 2
				   """,
				equate = False,
				checker_function = _is_int,               #this should be int
				),
			_Option(
				["-cmax", "cmax"],
				"""Specify maximum number of categories for FreeRate model.
				   It is recommended to increase if alignment is long enough

				   Default: 10
				   """,
				equate = False,
				checker_function = _is_int,              #this should be int
				),
			_Option(
				["-merit", "merit"],
				"""Specify either AIC, AICc or BIC for the optimality criterion to apply for new
				   procedure.
				   
				   Default: all three criteria are considered
				   """,
				equate = False,
				checker_function = lambda x: x in ("AIC", "AICc", "BIC")
				),
			_Switch(
				["-mtree", "mtree"],
				"""Turn on full tree search for  each model considered, to obtain a more accurate result.
				   Only recommended if enouhg computational resources are available.
				   
				   Default: Fixed starting tree
				   """
				),
			_Switch(
				["-mredo","mredo"],
				"""Ignore model checkpoint file computed earlier.
				
				   Default: model checkpoint file (if exists) is loaded to reuse previous computations
				   """
				),
			_Option(
				["-madd", "madd"],
				"""Specify a comma-separated list of mixture models to additionally consider for model
				   selection.
				   
				   Example: '-madd LG4M,LG4X' to additionally include these two protein mixture models
				   """,
				equate = False,
				),
			_Option(
				["-mdef","mdef"],
				"""Specify a NEXUS model file to define new models
				""",
				equate = False,
				filename = True
				),
			_Switch(
				["-mwopt", "mwopt"],
				"""Turn on optimizing weights for mixture models
				""",
				),
			_Option(
				["-a", "a"],    #idk if this is always a number or what
				"""Specify the Gamma shape parameter
				   Default: estimate""",
				equate = False,
				),
			_Option(
				["-amin", "amin"],
				"""Minimum Gamma shape parameter for site rates
				   Default: 0.02""",
				equate = False,
				checker_function = _is_number,
				),
			_Switch(
				["-gmedian", "gmedian"],
				"""Perform the median approximation for Gamma rate heterogeneity
				   Default: mean approximation""",
				),
			_Option(
				["-i", "i"],
				"""Specify the proportion of invariable sites
				   Default: estimate""",
				equate = False,
				),
			_Switch(
				["--opt-gamma-inv", "optgammainv"],
				"""Perform a more thorough estimation for +I+G model parameters
				""",
				),
			_Switch(
				["-wsr", "wsr"],
				"""Write per -site rates to .rate file
				""",
				),
			_Switch(
				["-mh", "mh"],                    
				"""Computing site-specific rates to .mhrate file using Meyer & von Haeseler method""",
				),

			### Partition model options ###

			_Option(
				["-q", "q"],
				"""Specify partition file for edge-equal partition model.
				   That means, all partitions share the same set of branch lengths
				   Edge-linked partition model (file in NEXUS/RAxML format)""",
				equate = False,
				filename = True,
				),
			_Option(
				["-spp", "spp"],
				"""Like -q but allowing partitions to have different evolutionary speeds
				   (Edge proportional partition model)
				   """,
				equate = False,
				filename = True,
				),
			_Option(
				["-ft", "ft"],
				"""Specify a guide tree (in Newick format) to infer site frequency profiles
				""",
				equate = False,
				filename = True,
				),
			_Option(
				["-fs", "fs"],
				"""Specify a site frequency file
				   Example: the .sitefreq file obtained from -ft run. This will save memory used
				   for the first phase of the analysis
				   """,
				equate = False,
				filename = True,
				),
			_Switch(
				["-fmax", "fmax"],
				"""Switch to posterior maximum mode for obtaining site-specific profiles
				   Default: Posterior Mean.
				   """,
				),

			### Tree search parameters ###

			_Switch(
				["-allnni", "allnni"],
				"""Turn on more thorough and slower NNI search.
				   It means that IQTree will consider all possible NNIs instead of only those in the
				   vicinity of previously applied NNIs.

				   Default: OFF.
				   """,
				),
			_Switch(
				["-djc", "djc"],
				"""Avoid computing ML pairwise distances and BIONJ tree""",
				),
			_Switch(
				["-fast", "fast"],
				"""Turn on the fast tree search mode, where IQTree will just construct two startign trees:
				   maximum parsimony and BIONJ, which are then optimized by nearest neighbor interchange
				   (NNI).
				   """,
				   ),
			_Option(
				["-g", "g"],
				"""Specify a topological constraint tree file in NEWICK format.
				   The constraint tree can be a multifurcating tree and need to not includ all taxa.
				   """,
				equate = False,
				filename = True,
				),
			_Option(
				["-ninit", "ninit"],
				"""Specify number of initial parsimony trees
				   The PLL Library is used, which implements the randomized stepwise addition and parsimony
				   subtree pruning and regafting (SPR)

				   Default: 100
				   """,
				equate = False,
				checker_function = _is_int,      #im assuming you can't input half a tree so it should be int
				),
			_Option(
				["-ntop", "ntop"],
				"""Specify number of top initial parsimony trees to optimize with ML nearest neighbor interchange
				   (NNI) search to initialize the candidate set.

				   Default: 20
				   """,
				equate = False,
				checker_function = _is_int,                   #assuming this has to be an int
				),
			_Option(
				["-nbest", "nbest"],
				"""Specify number of trees in the candidate set to maintain during ML tree search.

				   Default: 5
				   """,
				equate = False,
				checker_function = _is_int,
				),
			_Option(
				["-nstop", "nstop"],
				"""Specify number of unsuccessful iterations to stop.

				   Default: 100
				   """,
				equate = False,
				checker_function = _is_int,
				),
			_Option(
				["-pers", "pers"],
				"""Specify perturbation strength (between 0 and 1) for randomized NNI.

				   Default: 0.5
				   """,
				equate = False,
				checker_function = lambda x: True if (x in range(0, 1)) else False
				),
			_Option(
				["-sprrad", "sprrad"],
				"""Specify SPR radius for the initial parsimony tree serach

				   Default: 6
				   """,
				equate = False,
				checker_function= _is_number,
				),

			### Ultrafast booststrap parameters ###


			_Option(
				["-bb", "bb"],
				"""Specify a number of bootstrap replicates (>= 1000)
				   """,
				checker_function = _is_number   #not sure if it must be integer or not
				),
			_Option(
				["-bcor", "bcor"],
				"""Specify minimum correlation coefficient for UFBoot convergence criterion

				   Default: 0.99
				   """,
				equate = False,
				checker_function = _is_number
				),
			_Option(
				["-beps", "beps"],
				"""Specify a small epsilon to break tie in RELL evaluation for bootstrap trees.

				   Default: 0.5
				   """,
				equate = False,
				checker_function = _is_number
				),
			_Switch(
				["-bnni", "bnni"],
				"""Perform an additional step to further optimize UFBoot trees by nearest neighbor interchange
				   (NNI) based directly on bootstrap alignments.
				   This option is recommended in the presence of severe model violations.
				   It increases computing times by 2-fold but reduces the risk of overestimating branch supports
				   due to severe model violations
				   """,
				),
			_Option(
				["-bsam", "bsam"],
				"""Specify the resampling strategies for partitioned analysis.
				   By default IQTree resamples alignment sites within partitions.
				     -bsam GENE           IQTree will resample partitions
					 -bsam GENESITE       IQTree will resample partitions and then resample sites within 
										  resampled partitons
				   """,
				equate = False,
				checker_function = lambda x: x in ("GENE", "GENESITE"),
				),
			_Option(
				["-nm", "nm"],
				"""Specify maximum number of iterations to stop.

				   Default: 1000
				   """,
				equate = False,
				checker_function = _is_int
				),
			_Option(
				["-nstep", "nstep"],
				"""Specify iteration interval checking for UFBoot convergence.

				   Default: every 1000 iterations.
				   """,
				equate = False,
				checker_function = _is_int
				),
			_Switch(
				["-wbt", "wbt"],
				"""Turn on writing bootstrap tree to .ufboot file

				   Default: OFF
				   """,
				),
			_Switch(
				["-wbtl", "wbtl"],
				"""Like -wbt but booststrap trees written with branch lengths.""",
				),
			_Option(                                 
				["-j", "j"],
				"""Proportion of sites for jackknife

				   Default: NONE
				   """,
				equate = False,
				),

			### Nonpoarametric bootstrap ###

			_Option(
				["-b", "b"],
				"""Specify number of bootstrap replicates (recommended >= 100).
				   This will perfom both bootstrap and analysis on original alignment and
				   provide a consensus tree.
				   """,
				equate = False,
				checker_function = _is_int      #im assuming this requires an integer
				),
			_Option(
				["-bs", "bs"],
				"""Like -b but omit analysis on original alignment.""",
				equate = False,
				checker_function = _is_int     #same as above
				),
			_Option(
				["-bo", "bo"],
				"""Like -b but only perform bootstrap analysis (no analysis on original alignment
				   and no consensus tree)
				   """,
				equate = False,
				checker_function = _is_int,
				),
			_Option(
				["-bc", "bc"],
				"""Like -b but with no ML tree, only consensus tree""",
				equate = False,
				checker_function = _is_int,
				),

			### Single branch tests ###

			_Option(
				["-alrt", "alrt"],
				"""Specify number of replicates (>=1000) to perform SH-like approximate likelihood
				   ratio test (SH-aLRT).
				   If number of replicates is set to 0 (-alrt 0), then the parametric aLRT test is
				   performed, instead of SH-aLRT.
				   """,
				equate = False,
				checker_function = _is_int,
				),
			_Switch(
				["-abayes", "abayes"],
				"""Perform approximate Bayes test""",
				),
			_Option(
				["-lbp", "lbp"],
				"""Specify number of replicates (>=1000) to perform fast local bootstrap probability
				   method""",
				equate = False,
				checker_function = _is_int,
				),

			### Ancestral sequence reconstruction ###

			_Switch(
				["-asr", "asr"],
				"""Write ancestral sequences (by empirical Bayesian method) for all nodes of the
				   tree to .state file
				"""
				),
			_Option(
				["-asr-min","asrmin"],
				"""Specify the minimum threshold of posterior probability to determine the best ancestral
				   state.

				   Default: observed state frequency from the alignment""",
				equate = False,
				checker_function = _is_number,       #assuming probability is a float
				),

			### Tree topology tests ###

			_Option(
				["-z", "z"],
				"""Specify a file containing a set of trees.
				   IQTree will compute the log-likelihoods of all trees
				   """,
				equate = False,
				filename = True,
				),
			_Option(
				["-zb", "zb"],
				"""Specify the number of RELL replicates (>=1000) to perform several tere topology tests for
				   all trees passed via -z.
				   The tests include bootstrap proportion (BP), KH test, SH test and expected likelihood weights
				   (ELW).
				   """,
				equate = False,
				checker_function = _is_int,    #Assuming this is an int
				),
			_Switch(
				["-zw", "zw"],
				"""Used together with -zb to additionally perform the weighted-KH and weighted-SH tests.
				   """,
				),
			_Switch(
				["-au", "au"],
				"""Used together with -zb to additionally perform the approximately unbiased (AU) test.
				   Note that you have to specify the number of replicates for the AU test via -zb
				   """,
				),

			### Constructing consensus tree ###

			_Switch(
				["-con", "con"],
				"""Compute consensus tree of the trees passed via -t.
				   Resulting consensus tree is written to .contree file.
				   """
				),
			_Switch(
			    ["-net", "net"],
				"""Compute consensus network of the trees passed via -t.
				   Resulting consensus network is written to .nex file.
				   """
				),
			_Option(
				["-minsup", "minsup"],
				"""Specify a minimum threshold (between 0 and 1) to keep branches in the consensus tree.
				   -minsup 0.5 means to compute majority-rule consensus tree.

				   Default: 0 to compute extended majority-rule consensus.
				   """,
				equate = False,
				checker_function = lambda x: True if (x in range(0, 1)) else False,
				),
			_Option(
				["-bi", "bi"],
				"""Specify a burn-in, which is the number of beginning trees passed via -t to discard
				   before consensus construction.
				   This is useful e.g. when summarizing trees from MrBayes analysis.
				   """,
				equate = False,
				checker_function = _is_int,
				),
			_Option(
				["-sup", "sup"],
				"""Specify an input "target" tree file.
				   That means, support values are first extracted from the trees passed via -t,
				   and then mapped onto the target tree.

				   Resulting tree with assigned support values is written to .suptree file.

				   This option is useful to map and compare support values from different approaches onto
				   a single tree
				   """,
				equate = False,
				filename = True,
				),
			_Option(
				["-suptag", "suptag"],
				"""Specify name of a node in -sup target tree.
				   The corresponding node of .suptree will then be assigned with IDs of trees where this node
				   appears.

				   Special option -suptag ALL will assign such IDs for all nodes of the target tree.
				   """,
				equate = False,
				),
			_Option(
				["-scale", "scale"],
				"""Set the scaling factor of support values for -sup option

				   Default: 100, i.e. percentages
				   """,
				equate = False,
				checker_function = _is_int,
				),
			_Option(
				["-prec", "prec"],
				"""Set the precision of support values for -sup option.

				   Default: 0
				   """,
				equate = False,
				checker_function = _is_int
				),
			_Option(
				["-rf", "rf"],
				"""Specify a second set of trees.
				   IQTree computes all pairwise RF distances between two tree sets passed via -t and -rf
				   """,
				equate = False,
				filename = True,
				),
			_Switch(
				["-rf_all", "rf_all"],
				"""Compute all-to-all RF distances between all tress passed via -t""",
				),
			_Switch(
				["-rf_adj", "rf_adj"],
				"""Compute RF distances between adjacent trees passed via -t""",
				),

			### Generating random trees ###

			_Option(
				["-r", "r"],
				"""Specify number of taxa.
				   IQTree will create a random tree under Yule-Harding model with specified number of taxa
				   """,
				equate = False,
				checker_function = _is_number  #not sure if it's supposed to be an int
				),
			_Option(
				["-ru", "ru"],
				"""Like -r but a random tree is created under uniform model""",
				equate = False,
				checker_function = _is_number  #not sure if it's supposed to be an int
				),
			_Option(
				["-rcat", "rcat"],
				"""Like -r, but a random caterpillar tree is created""",
				equate = False,
				checker_function = _is_number
				),
			_Option(
				["-rbal", "rbal"],
				"""Like -r, but a random balanced tree is created""",
				equate = False,
				checker_function = _is_number
				),
			_Option(
				["-rcsg", "rcsg"],
				"""Like -r, but a random circular split network is created""",
				equate = False,
				checker_function = _is_number
				),
			_Option(
				["-rlen", "rlen"],
				"""Specify three numbers: minimum, mean and maximum branch lengths of the random tree

				   Defaul: -rlen 0.001 0.1 0.999
				   """,
				equate = False,        
				),

			### Miscellaneous options ###

			_Switch(
				["-alninfo", "alninfo"],
				"""Print alignment site statistics to .alninfo file""",
				),
			_Switch(
				["-blfix", "blfix"],
				"""Fix branch lengths of tree passed via -t or -te.
				   This is useful to evaluate the log-likelihood of an input tree with fixed topology
				   and branch lengths.

				   Default: OFF
				   """,
				),
			_Option(
				["-blmin", "blmin"],
				"""Specify minimum branch length.

				   Default: the smaller of 0.000001 and 0.1/alignment_length
				   """,
				equate = False,
				checker_function = _is_number,
				),
			_Option(
				["-blmax", "blmax"],
				"""Specify the maximum branch length.

				   Default: 10
				   """,
				equate = False,
				checker_function = _is_number,
				),
			_Switch(
				["-czb", "czb"],
				"""Collapse near zero branches, so that the final tree may be multifurcating.
				   This is useful for bootstrapping in the presence of polytomy to reduce bootstrap
				   supports of short branches.
				   """,
				),
			_Option(
				["-me", "me"],
				"""Specify the log-likelihood epsilon for final model parameter estimation.
				   With -fast option the epsilon is raised to 0.05

				   Default: 0.01
				   """,
				equate = False,
				checker_function = _is_number,
				),
			_Switch(
				["-wpl", "wpl"],
				"""Write partition log-likelihoods to .partlh file.
				   Only effective with partition model.
				   """,
				),
			_Switch(
				["-wspr", "wspr"],
				"""Write site posterior probabilities per rate category to .siteprob file""",
				),
			_Switch(
				["-wspm", "wspm"],
				"""Write site posterior probabilities per mixture class to .siteprob file.""",
				),
			_Switch(
				["-wspmr", "wspmr"],
				"""Write site posterior probabilities per mixture class and rate category to .siteprob file""",
				),
			_Switch(
				["-wsl", "wsl"],
				"""Write size log-likelihoods to .sitelh file in TREE-PUZZLE format.
				   Such file can then be passed on to CONSEL for further tree tests.
				   """,
				),
			_Switch(
				["-wslr", "wslr"],
				"""Write site log likelihoods per rate category to .sitelh file""",
				),
			_Switch(
				["-wslm", "wslm"],
				"""Write site log likelihoods per mixture class to .sitelh file.""",
				),
			_Switch(
				["-wslmr", "wslmr"],
				"""Write site log likelihoods per mixture class and rate category to .sitelh file""",
				),
			_Switch(
				["-wt", "wt"],
				"""Turn on writing all locally optimal trees into .treels file.""",
				),
			_Option(
				["-fconst", "fconst"],
				"""Specify a list of comma-separated integer numbers.
				   The number of entries should be equal to the number of states in the model
				   (e.g. 4 for DNA and 20 for protein).
				   IQTree will then add a number of constant sites accordingly.

				   For example, -fconst 10,20,15,40 will add 10 constant sites of all A, 20 constant sites of all C,
				   15 constant sites of all G and 40 constant sites of all T into the alignment.
				   """,
				equate = False,
				checker_function = lambda x: True if (all(y.isdigit() for y in x.split(","))) else False,
				),

			#These weren't on the documentation but are present on -h list of commands

			_Switch(
				["--show-lh", "showlh"],
				"""Compute tree likelihood without optimisation""",
				),
			_Switch(
				["--eigenlib", "eigenlib"],
				"""Use Eigen3 library""",
				),
			_Switch(
				["--no-outfiles", "nooutfiles"],
				"""Suppress printing output files""",
				),
			_Option(
				["-nni-eval", "nnieval"],
				"""Specify n times to loop for NNI evaluation 
				   
				   Default: 1 
				   """,
				equate = False,
				checker_function = _is_int,
				),
			_Option(
				["--runs", "runs"],
				"""Number of independent runs

				   Default: 1
				   """,
				equate = False,
				checker_function = _is_int,
				),
		]

		AbstractCommandline.__init__(self, cmd, **kwargs)
