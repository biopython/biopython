
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


class IQTreeCommandline(AbstractCommandline):
	r"""Command-line wrapper for IQTree.

	Specify command is mandatory, used to specify input alignment file.

	From the terminal command line use ``iqtree -h`` or ``iqtree -?``
	for more explanation of usage options and commands

	Homepage http://www.iqtree.org/

	References:


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
				["-s", "specify"],   
				"""Specify input alignment file in PHYLIP, FASTA, NEXUS
				   CLUSTAL or MSF format
				   """,
				is_required = True,
				equate = False,
				filename=True
				),  
			_Option(
				["-st", "sequencetype"],
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
				["-t", "tree"],
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
				["-te", "usertree"],
				"""Like -t but fixing user tree, no tree search is performed
				   and the program computes the log-likelihood of the fixed user tree
				   """,
				equate = False,
				filename = (lambda x: not(x in ("BIONJ", "RANDOM","PARS", "PLLPARS"))),                            #Not sure if this works
				checker_function= lambda x: True if filename else (x in ("BIONJ", "RANDOM","PARS", "PLLPARS"))     #Need tree file to test
				),
			_Option(
				["-o", "outgroup"],
				"""Specify an outgroup taxon name to root the tree
				   Output tree will be rooted accordingly

				   Default: first taxon in alignment
				   """,
				equate = False
				),
			_Option(
				["-pre", "prefix"],                        #This outputs the files in the project folder for some reason
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
				["-keep-ident", "keep"],
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
				["-mem", "memory"],
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
				["-cptime", "checkpoint_interval"],
				"""Specify the minimum checkpoint time interval in seconds

				   Default 20 seconds
				   """,
				equate = False,
				checker_function = _is_int
				),

			#Likelihood mapping analysis comands

			_Option(
				["-lmap", "lmap"],
				"""Specify the number of quartets to be randomly drawn.
				   If you specify -lmap ALL, all unique quartets will be drawn, instead.
				   """,
				equate = False,
				checker_function = lambda x: isinstance(x, int) or x.isdigit() or x in ("ALL")
				),
			_Option(
				["-lmclust", "taxon_clusters"],
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

			#Automatic model selection

			_Option(
				["-m", "model"],
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
				["-rcluster-max", "rcluster-max"],
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

				   Example: -mset WAG, LG, JTT will restrict model selection to WAG, LG and JTT instead of
				            all 18 AA models to save computations
				   """,

				equate = False
				#might need a checker function
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
				),    #idk about a checker function for this one
			_Option(
				["-mrate", "mrate"],
				"""Specify a comma-separated list of rate heterogeneity types for model selection

				   Default: -mrate E,I,G,I+G for standard procedure
				            -mrate E,I,G,I+G,R for new selection procedure
				   """,
				equate = False
				),
			_Option(
				["-cmin", "cmin"],
				"""Specify minimum number of categories of FreeRate model

				   Default: 2
				   """,
				equate = False
				),
			_Option(
				["-cmax", "cmax"],
				"""Specify maximum number of categories for FreeRate model.
				   It is recommended to increase if alignment is long enough

				   Default: 10
				   """,
				equate = False
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
				["-a", "a"],    #Maybe change a to gamma also idk if this is always a number or what
				"""Specify the Gamma shape parameter
				   Default: estimate""",
				equate = False,
				),
			_Option(
				["-amin", "amin"],
				"""Minimum Gamma shape parameter for site rates
				   Default: 0.02""",
				equate = False,
				checker_function = lambda x: x.isnumeric(),
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
				["-mh", "mh"],                     #not in documentation, found in -h list of commands
				"""Computing site-specific rates to .mhrate file using Meyer & von Haeseler method""",
				),

			### Partition model options 

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
				["-spp", "spp"],
				"""Specify partition file for edge-unlinked partition model.
				   That means, each partition has its own set of branch lengths.
				   This is the mot parameter-rich partition model to accomodate heterotachy
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

			#Tree search parameters

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
		]

		AbstractCommandline.__init__(self, cmd, **kwargs)
