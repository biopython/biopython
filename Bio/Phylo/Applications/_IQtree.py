# Based on code in _Phyml.py by Eric Talevich.
# All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Command-line wrapper for the phylogenomic inference IQ-Tree"""

from Bio.Application import _Option, _Switch, _Argument, AbstractCommandline

class IQTreeCommandline(AbstractCommandline):
	r"""Command-line wrapper for IQTree.

	Input and Output parameters are mandatory.

	From the terminal command line use ``iqtree -h`` or ``iqtree -?``
	for more explanation of usage options and commands

	Homepage http://www.iqtree.org/

	References:


	"""

	def __init__(self, cmd="iqtree", **kwargs):
		"""initialize the class"""

		self.parameters = [
			_Option(						#INPUT
				["-s", "Specify"],   
				"""Specify input alignment file in PHYLIP, FASTA, NEXUS
				   CLUSTAL or MSF format""",
				   is_required = True,
				   equate = False,
				   filename=True,
				),  
			_Option(
				["-st", "SequenceType"],
				"""Specify sequence type as either DNA, AA, BIN, MORPH, 
			       CODON or NT2AA for DNA, ammino-acid, binary, morphological
				   codon or DNA-to-AA translate sequences
				   
				   Necessary only if IQ-Tree did not detect the sequence correctly
				   
				   Note: This is always necessary when using codon models otherwise
				         IQ-tree applies DNA models""",
				equate = False,
				checker_function= lambda x: x in ("DNA", "AA", "BIN", "MORPH", "CODON", "NT2AA", "DNA-to-AA"),
				),
			_Option(
				["-t", "Tree"],
				"""Specify a file containing starting tree for tree search.
					BIONJ starts a tree search from BIONJ tree
					RANDOM starts tree search from completely random tree
					
					Default: 100 parsimony trees + BIONJ tree""",
					#TODO: checker function with filename if in it
				filename=True,
				),
			_Option(
				["-te", "UserTree"],
				"""Like -t but fixing user tree, no tree search is performed
				   and the program computes the log-likelihood of the fixed user tree""",
				   filename = True,
				   #TODO: checker function with filename if in it like -t
				),
			_Option(
				["-o", "outgroup"],
				"""Specify an outgroup taxon name to root the tree
				   Output tree will be rooted accordingly

				   Default: first taxon in alignment""",
				equate = False,
				),
			_Option(
				["-pre", "prefix"],                       #Not really sure how this one works
				"""Specify a prefix for all output files

				   Default: either alignment file name (-s) or partition file name
				   (-q, -spp or -sp)""",
				equate = False,
					),
			_Option(
				["-nt", "nt"],
				""" Specify the number of CPU cores for the multicore version.
					A special option -nt AUTO will tell IQTree to automatically determine
					the best number of cores given the current data and computer""",
				checker_function = (lambda x: isinstance(x, int) or x.isdigit()),
				equate = False,
				),
			_Option(
				["-ntmax", "ntmax"],
				"""Specify the maximal number of CPU cores.
					Default: # CPU cores on the current machine""",
				checker_function = (lambda x: isinstance(x, int) or x.isdigit()),
				equate= False,
				),
			_Switch(
				["-v", "verbose"],
				"""Turn on verbose mode for printing more messages to screen, used for debugging
					Default: OFF""",
				),
			_Switch(
				["-quiet", "quiet"],
				"""Silent mode, suppress printing to the screen, note that .log file is still 
					written""",
				),
			_Switch(
				["-keep-ident", "keep"],
				"""Keep identical sequences in the alignment.
				   Default: IQTree will remove identical sequences during the analysis and add
				   them at the end""",
				),
			_Switch(
				["-safe","safe"],
				"""Turn on safe numerical mode to avoid numerical underflow for large data sets
				   with many sequences.
				   This mode is automatically turned on when having more than 2000 sequences""",
				),
			_Option(
				["-mem", "memory"],
				"""Specify maximal RAM usage 
				   Example:
						-mem 64G to use at most 64 GB of RAM
						-mem 200M to use at most 200 MB of RAM
						-and so on

				   Default: IQTree will not exceed the computer RAM size""",
				equate= False,
				),

		]

		AbstractCommandline.__init__(self, cmd, **kwargs)
