# Copyright 2011 by Andreas Wilm. All rights reserved.
# Based on ClustalW wrapper copyright 2009 by Cymon J. Cox.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Command line wrapper for the multiple alignment program Clustal Omega."""


from Bio.Application import _Option, _Switch, AbstractCommandline


class ClustalOmegaCommandline(AbstractCommandline):
    """Command line wrapper for clustal omega.

    http://www.clustal.org/omega

    Notes
    -----
    Last checked against version: 1.2.0

    References
    ----------
    Sievers F, Wilm A, Dineen DG, Gibson TJ, Karplus K, Li W, Lopez R,
    McWilliam H, Remmert M, Söding J, Thompson JD, Higgins DG (2011).
    Fast, scalable generation of high-quality protein multiple
    sequence alignments using Clustal Omega.
    Molecular Systems Biology 7:539 https://doi.org/10.1038/msb.2011.75

    Examples
    --------
    >>> from Bio.Align.Applications import ClustalOmegaCommandline
    >>> in_file = "unaligned.fasta"
    >>> out_file = "aligned.fasta"
    >>> clustalomega_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, verbose=True, auto=True)
    >>> print(clustalomega_cline)
    clustalo -i unaligned.fasta -o aligned.fasta --auto -v

    You would typically run the command line with clustalomega_cline() or via
    the Python subprocess module, as described in the Biopython tutorial.

    """

    def __init__(self, cmd="clustalo", **kwargs):
        """Initialize the class."""
        # order parameters in the same order as clustalo --help
        self.parameters = [
            # Sequence Input
            _Option(
                ["-i", "--in", "--infile", "infile"],
                "Multiple sequence input file",
                filename=True,
                equate=False,
            ),
            _Option(
                ["--hmm-in", "HMM input", "hmm_input"],
                "HMM input files",
                filename=True,
                equate=False,
            ),
            _Switch(["--dealign", "dealign"], "Dealign input sequences"),
            _Option(
                ["--profile1", "--p1", "profile1"],
                "Pre-aligned multiple sequence file (aligned columns will be kept fix).",
                filename=True,
                equate=False,
            ),
            _Option(
                ["--profile2", "--p2", "profile2"],
                "Pre-aligned multiple sequence file (aligned columns will be kept fix).",
                filename=True,
                equate=False,
            ),
            _Option(
                ["-t", "--seqtype", "seqtype"],
                "{Protein, RNA, DNA} Force a sequence type (default: auto).",
                equate=False,
                checker_function=lambda x: x
                in ["protein", "rna", "dna", "Protein", "RNA", "DNA", "PROTEIN"],
            ),
            _Switch(
                ["--is-profile", "isprofile"],
                "disable check if profile, force profile (default no)",
            ),
            _Option(
                ["--infmt", "infmt"],
                """Forced sequence input file format (default: auto)

                    Allowed values: a2m, fa[sta], clu[stal], msf, phy[lip], selex, st[ockholm], vie[nna]
                    """,
                equate=False,
                checker_function=lambda x: x
                in [
                    "a2m",
                    "fa",
                    "fasta",
                    "clu",
                    "clustal",
                    "msf",
                    "phy",
                    "phylip",
                    "selex",
                    "st",
                    "stockholm",
                    "vie",
                    "vienna",
                ],
            ),
            # Clustering
            _Option(
                ["--distmat-in", "distmat_in"],
                "Pairwise distance matrix input file (skips distance computation).",
                filename=True,
                equate=False,
            ),
            _Option(
                ["--distmat-out", "distmat_out"],
                "Pairwise distance matrix output file.",
                filename=True,
                equate=False,
            ),
            _Option(
                ["--guidetree-in", "guidetree_in"],
                "Guide tree input file (skips distance computation and guide-tree clustering step).",
                filename=True,
                equate=False,
            ),
            _Option(
                ["--guidetree-out", "guidetree_out"],
                "Guide tree output file.",
                filename=True,
                equate=False,
            ),
            _Switch(
                ["--full", "distmat_full"],
                "Use full distance matrix for guide-tree calculation (slow; mBed is default)",
            ),
            _Switch(
                ["--full-iter", "distmat_full_iter"],
                "Use full distance matrix for guide-tree calculation during iteration (mBed is default)",
            ),
            _Option(
                ["--cluster-size", "clustersize"],
                "soft maximum of sequences in sub-clusters",
                checker_function=lambda x: isinstance(x, int),
            ),
            _Option(
                ["--clustering-out", "clusteringout"],
                "Clustering output file",
                filename=True,
            ),
            _Switch(
                ["--use-kimura", "usekimura"],
                "use Kimura distance correction for aligned sequences (default no)",
            ),
            _Switch(
                ["--percent-id", "percentid"],
                "convert distances into percent identities (default no)",
            ),
            # Alignment Output
            _Option(
                ["-o", "--out", "--outfile", "outfile"],
                "Multiple sequence alignment output file (default: stdout).",
                filename=True,
                equate=False,
            ),
            _Option(
                ["--outfmt", "outfmt"],
                "MSA output file format:"
                " a2m=fa[sta],clu[stal],msf,phy[lip],selex,st[ockholm],vie[nna]"
                " (default: fasta).",
                equate=False,
                checker_function=lambda x: x
                in [
                    "a2m",
                    "fa",
                    "fasta",
                    "clu",
                    "clustal",
                    "msf",
                    "phy",
                    "phylip",
                    "selex",
                    "st",
                    "stockholm",
                    "vie",
                    "vienna",
                ],
            ),
            _Switch(
                ["--residuenumber", "--resno", "residuenumber"],
                "in Clustal format print residue numbers (default no)",
            ),
            _Option(
                ["--wrap", "wrap"],
                "number of residues before line-wrap in output",
                checker_function=lambda x: isinstance(x, int),
            ),
            _Option(
                ["--output-order", "outputorder"],
                "MSA output order like in input/guide-tree",
                checker_function=lambda x: x in ["input-order", "tree-order"],
            ),
            # Iteration
            _Option(
                ["--iterations", "--iter", "iterations"],
                "Number of (combined guide-tree/HMM) iterations",
                equate=False,
                checker_function=lambda x: isinstance(x, int),
            ),
            _Option(
                ["--max-guidetree-iterations", "max_guidetree_iterations"],
                "Maximum number of guidetree iterations",
                equate=False,
                checker_function=lambda x: isinstance(x, int),
            ),
            _Option(
                ["--max-hmm-iterations", "max_hmm_iterations"],
                "Maximum number of HMM iterations",
                equate=False,
                checker_function=lambda x: isinstance(x, int),
            ),
            # Limits (will exit early, if exceeded):
            _Option(
                ["--maxnumseq", "maxnumseq"],
                "Maximum allowed number of sequences",
                equate=False,
                checker_function=lambda x: isinstance(x, int),
            ),
            _Option(
                ["--maxseqlen", "maxseqlen"],
                "Maximum allowed sequence length",
                equate=False,
                checker_function=lambda x: isinstance(x, int),
            ),
            # Miscellaneous:
            _Switch(
                ["--auto", "auto"],
                "Set options automatically (might overwrite some of your options)",
            ),
            _Option(
                ["--threads", "threads"],
                "Number of processors to use",
                equate=False,
                checker_function=lambda x: isinstance(x, int),
            ),
            _Option(
                ["-l", "--log", "log"],
                "Log all non-essential output to this file.",
                filename=True,
                equate=False,
            ),
            _Switch(["-h", "--help", "help"], "Print help and exit."),
            _Switch(["-v", "--verbose", "verbose"], "Verbose output"),
            _Switch(["--version", "version"], "Print version information and exit"),
            _Switch(
                ["--long-version", "long_version"],
                "Print long version information and exit",
            ),
            _Switch(["--force", "force"], "Force file overwriting."),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
