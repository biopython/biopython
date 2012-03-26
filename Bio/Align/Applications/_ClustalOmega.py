# Copyright 2011 by Andreas Wilm. All rights reserved.
# Based on ClustalW wrapper copyright 2009 by Cymon J. Cox.
#
# Wrapper for Clustal Omega by Andreas Wilm (2011). Used _Clustalw.py
# as template.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Command line wrapper for the multiple alignment program Clustal Omega.
"""

from Bio.Application import _Option, _Switch, AbstractCommandline

class ClustalOmegaCommandline(AbstractCommandline):
    """Command line wrapper for clustal omega

    http://www.clustal.org/omega

    Example:

    >>> from Bio.Align.Applications import ClustalOmegaCommandline
    >>> in_file = "unaligned.fasta"
    >>> out_file = "aligned.fasta"
    >>> clustalomega_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, verbose=True, auto=True)
    >>> print clustalomega_cline
    clustalo -i unaligned.fasta -o aligned.fasta --auto -v


    You would typically run the command line with clustalomega_cline() or via
    the Python subprocess module, as described in the Biopython tutorial.

    Citation:

    Sievers F, Wilm A, Dineen D, Gibson TJ, Karplus K, Li W, Lopez R,
    McWilliam H, Remmert R, Soding J, Thompson JD Higgins DG
    Fast, scalable generation of high-quality protein multiple
    sequence alignments using Clustal Omega.
    Molecular Systems Biology 2011; accepted. 

    Last checked against versions: 1.0.3
    """
    def __init__(self, cmd="clustalo", **kwargs):
        # order parameters in the same order as clustalo --help
        self.parameters = \
            [
            # Sequence Input
            _Option(["-i", "--in", "--infile", "infile"],
                    "Multiple sequence input",
                    filename=True,
                    equate=False),
            _Switch(["--dealign", "dealign"],
                    "Dealign input sequences"),
            _Option(["--hmm-in", "HMM input", "hmm_input"],
                    "HMM input files",
                    filename=True,
                    equate=False),
            _Option(["--profile1", "--p1", "profile1"],
                    "Pre-aligned multiple sequence file (aligned columns will be kept fix).",
                    filename=True,
                    equate=False),
            _Option(["--profile2", "--p2", "profile2"],
                    "Pre-aligned multiple sequence file (aligned columns will be kept fix).",
                    filename=True,
                    equate=False),

            # Clustering
            _Option(["--distmat-in", "distmat_in"],
                    "Pairwise distance matrix input file (skips distance computation).",
                    filename=True,
                    equate=False),
            _Option(["--distmat-out", "distmat_out"],
                    "Pairwise distance matrix output file.",
                    filename=True,
                    equate=False),
            _Option(["--guidetree-in", "guidetree_in"],
                    "Guide tree input file (skips distance computation and guide-tree clustering step).",
                    filename=True,
                    equate=False),
            _Option(["--guidetree-out", "guidetree_out"],
                    "Guide tree output file.",
                    filename=True,
                    equate=False),
            _Switch(["--full", "distmat_full"],
                    "Use full distance matrix for guide-tree calculation (might be slow; mBed is default)"),
            _Switch(["--full-iter", "distmat_full_iter"],
                    "Use full distance matrix for guide-tree calculation during iteration (might be slowish; mBed is default)"),

            # Alignment Output
            _Option(["-o", "--out", "--outfile", "outfile"],
                    "Multiple sequence alignment output file (default: stdout).",
                    filename=True,
                    equate=False),
            _Option(["--outfmt", "outfmt"],
                    "MSA output file format:"
                    " a2m=fa[sta],clu[stal],msf,phy[lip],selex,st[ockholm],vie[nna]"
                    " (default: fasta).",
                    equate=False,
                    checker_function=lambda x: x in ["a2m", "fa", "fasta", 
                                                     "clu", "clustal", 
                                                     "msf",
                                                     "phy", "phylip", 
                                                     "selex",
                                                     "st", "stockholm",
                                                     "vie", "vienna"]),
           # Iteration
            _Option(["--iterations", "--iter", "iterations"],
                    "Number of (combined guide-tree/HMM) iterations",
                    equate=False,
                    checker_function=lambda x: isinstance(x, int)),
            _Option(["--max-guidetree-iterations", "max_guidetree_iterations"],
                    "Maximum number of guidetree iterations",
                    equate=False,
                    checker_function=lambda x: isinstance(x, int)),
            _Option(["--max-hmm-iterations", "max_hmm_iterations"],
                    "Maximum number of HMM iterations",
                    equate=False,
                    checker_function=lambda x: isinstance(x, int)),

            # Limits (will exit early, if exceeded):
            _Option(["--maxnumseq", "maxnumseq"],
                    "Maximum allowed number of sequences",
                    equate=False,
                    checker_function=lambda x: isinstance(x, int)),
            _Option(["--maxseqlen", "maxseqlen"],
                    "Maximum allowed sequence length",
                    equate=False,
                    checker_function=lambda x: isinstance(x, int)),

            # Miscellaneous:

            _Switch(["--auto", "auto"],
                    "Set options automatically (might overwrite some of your options)"),
            _Option(["--threads", "threads"],
                    "Number of processors to use",
                    equate=False,
                    checker_function=lambda x: isinstance(x, int)),
            _Option(["-l", "--log", "log"],
                    "Log all non-essential output to this file.",
                    filename=True,
                    equate=False),
            _Switch(["-h", "--help", "help"],
                    "Outline the command line params."),
            _Switch(["-v", "--verbose", "verbose"],
                    "Verbose output"),
            _Switch(["--version", "version"],
                    "Print version information and exit"),
            _Switch(["--long-version", "long_version"],
                    "Print long version information and exit"),
            _Switch(["--force", "force"],
                    "Force file overwriting."),

            ]
        AbstractCommandline.__init__(self, cmd, **kwargs)

def _test():
    """Run the module's doctests (PRIVATE)."""
    print "Runing ClustalW doctests..."
    import doctest
    doctest.testmod()
    print "Done"

if __name__ == "__main__":
    _test()
