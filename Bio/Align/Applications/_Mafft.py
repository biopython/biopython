# Copyright 2009 by Cymon J. Cox.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Command line wrapper for the multiple alignment programme MAFFT."""


from Bio.Application import _Option, _Switch, _Argument, AbstractCommandline


class MafftCommandline(AbstractCommandline):
    """Command line wrapper for the multiple alignment program MAFFT.

    http://align.bmr.kyushu-u.ac.jp/mafft/software/

    Notes
    -----
    Last checked against version: MAFFT v6.717b (2009/12/03)

    References
    ----------
    Katoh, Toh (BMC Bioinformatics 9:212, 2008) Improved accuracy of
    multiple ncRNA alignment by incorporating structural information into
    a MAFFT-based framework (describes RNA structural alignment methods)

    Katoh, Toh (Briefings in Bioinformatics 9:286-298, 2008) Recent
    developments in the MAFFT multiple sequence alignment program
    (outlines version 6)

    Katoh, Toh (Bioinformatics 23:372-374, 2007)  Errata PartTree: an
    algorithm to build an approximate tree from a large number of
    unaligned sequences (describes the PartTree algorithm)

    Katoh, Kuma, Toh, Miyata (Nucleic Acids Res. 33:511-518, 2005) MAFFT
    version 5: improvement in accuracy of multiple sequence alignment
    (describes [ancestral versions of] the G-INS-i, L-INS-i and E-INS-i
    strategies)

    Katoh, Misawa, Kuma, Miyata (Nucleic Acids Res. 30:3059-3066, 2002)

    Examples
    --------
    >>> from Bio.Align.Applications import MafftCommandline
    >>> mafft_exe = "/opt/local/mafft"
    >>> in_file = "../Doc/examples/opuntia.fasta"
    >>> mafft_cline = MafftCommandline(mafft_exe, input=in_file)
    >>> print(mafft_cline)
    /opt/local/mafft ../Doc/examples/opuntia.fasta

    If the mafft binary is on the path (typically the case on a Unix style
    operating system) then you don't need to supply the executable location:

    >>> from Bio.Align.Applications import MafftCommandline
    >>> in_file = "../Doc/examples/opuntia.fasta"
    >>> mafft_cline = MafftCommandline(input=in_file)
    >>> print(mafft_cline)
    mafft ../Doc/examples/opuntia.fasta

    You would typically run the command line with mafft_cline() or via
    the Python subprocess module, as described in the Biopython tutorial.

    Note that MAFFT will write the alignment to stdout, which you may
    want to save to a file and then parse, e.g.::

        stdout, stderr = mafft_cline()
        with open("aligned.fasta", "w") as handle:
            handle.write(stdout)
        from Bio import AlignIO
        align = AlignIO.read("aligned.fasta", "fasta")

    Alternatively, to parse the output with AlignIO directly you can
    use StringIO to turn the string into a handle::

        stdout, stderr = mafft_cline()
        from io import StringIO
        from Bio import AlignIO
        align = AlignIO.read(StringIO(stdout), "fasta")

    """

    def __init__(self, cmd="mafft", **kwargs):
        """Initialize the class."""
        BLOSUM_MATRICES = ["30", "45", "62", "80"]
        self.parameters = [
            # **** Algorithm ****
            # Automatically selects an appropriate strategy from L-INS-i, FFT-NS-
            # i and FFT-NS-2, according to data size. Default: off (always FFT-NS-2)
            _Switch(["--auto", "auto"], "Automatically select strategy. Default off."),
            # Distance is calculated based on the number of shared 6mers. Default: on
            _Switch(
                ["--6merpair", "6merpair", "sixmerpair"],
                "Distance is calculated based on the number of shared "
                "6mers. Default: on",
            ),
            # All pairwise alignments are computed with the Needleman-Wunsch
            # algorithm. More accurate but slower than --6merpair. Suitable for a
            # set of globally alignable sequences. Applicable to up to ~200
            # sequences. A combination with --maxiterate 1000 is recommended (G-
            # INS-i). Default: off (6mer distance is used)
            _Switch(
                ["--globalpair", "globalpair"],
                "All pairwise alignments are computed with the "
                "Needleman-Wunsch algorithm. Default: off",
            ),
            # All pairwise alignments are computed with the Smith-Waterman
            # algorithm. More accurate but slower than --6merpair. Suitable for a
            # set of locally alignable sequences. Applicable to up to ~200
            # sequences. A combination with --maxiterate 1000 is recommended (L-
            # INS-i). Default: off (6mer distance is used)
            _Switch(
                ["--localpair", "localpair"],
                "All pairwise alignments are computed with the "
                "Smith-Waterman algorithm. Default: off",
            ),
            # All pairwise alignments are computed with a local algorithm with
            # the generalized affine gap cost (Altschul 1998). More accurate but
            # slower than --6merpair. Suitable when large internal gaps are
            # expected. Applicable to up to ~200 sequences. A combination with --
            # maxiterate 1000 is recommended (E-INS-i). Default: off (6mer
            # distance is used)
            _Switch(
                ["--genafpair", "genafpair"],
                "All pairwise alignments are computed with a local "
                "algorithm with the generalized affine gap cost "
                "(Altschul 1998). Default: off",
            ),
            # All pairwise alignments are computed with FASTA (Pearson and Lipman
            # 1988). FASTA is required. Default: off (6mer distance is used)
            _Switch(
                ["--fastapair", "fastapair"],
                "All pairwise alignments are computed with FASTA "
                "(Pearson and Lipman 1988). Default: off",
            ),
            # Weighting factor for the consistency term calculated from pairwise
            # alignments. Valid when either of --blobalpair, --localpair, --
            # genafpair, --fastapair or --blastpair is selected. Default: 2.7
            _Option(
                ["--weighti", "weighti"],
                "Weighting factor for the consistency term calculated "
                "from pairwise alignments. Default: 2.7",
                checker_function=lambda x: isinstance(x, float),
                equate=False,
            ),
            # Guide tree is built number times in the progressive stage. Valid
            # with 6mer distance. Default: 2
            _Option(
                ["--retree", "retree"],
                "Guide tree is built number times in the progressive "
                "stage. Valid with 6mer distance. Default: 2",
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            # Number cycles of iterative refinement are performed. Default: 0
            _Option(
                ["--maxiterate", "maxiterate"],
                "Number cycles of iterative refinement are performed. Default: 0",
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            # Number of threads to use. Default: 1
            _Option(
                ["--thread", "thread"],
                "Number of threads to use. Default: 1",
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            # Use FFT approximation in group-to-group alignment. Default: on
            _Switch(
                ["--fft", "fft"],
                "Use FFT approximation in group-to-group alignment. Default: on",
            ),
            # Do not use FFT approximation in group-to-group alignment. Default:
            # off
            _Switch(
                ["--nofft", "nofft"],
                "Do not use FFT approximation in group-to-group "
                "alignment. Default: off",
            ),
            # Alignment score is not checked in the iterative refinement stage.
            # Default: off (score is checked)
            _Switch(
                ["--noscore", "noscore"],
                "Alignment score is not checked in the iterative "
                "refinement stage. Default: off (score is checked)",
            ),
            # Use the Myers-Miller (1988) algorithm. Default: automatically
            # turned on when the alignment length exceeds 10,000 (aa/nt).
            _Switch(
                ["--memsave", "memsave"],
                "Use the Myers-Miller (1988) algorithm. Default: "
                "automatically turned on when the alignment length "
                "exceeds 10,000 (aa/nt).",
            ),
            # Use a fast tree-building method (PartTree, Katoh and Toh 2007) with
            # the 6mer distance. Recommended for a large number (> ~10,000) of
            # sequences are input. Default: off
            _Switch(
                ["--parttree", "parttree"],
                "Use a fast tree-building method with the 6mer "
                "distance. Default: off",
            ),
            # The PartTree algorithm is used with distances based on DP. Slightly
            # more accurate and slower than --parttree. Recommended for a large
            # number (> ~10,000) of sequences are input. Default: off
            _Switch(
                ["--dpparttree", "dpparttree"],
                "The PartTree algorithm is used with distances "
                "based on DP. Default: off",
            ),
            # The PartTree algorithm is used with distances based on FASTA.
            # Slightly more accurate and slower than --parttree. Recommended for
            # a large number (> ~10,000) of sequences are input. FASTA is
            # required. Default: off
            _Switch(
                ["--fastaparttree", "fastaparttree"],
                "The PartTree algorithm is used with distances based "
                "on FASTA. Default: off",
            ),
            # The number of partitions in the PartTree algorithm. Default: 50
            _Option(
                ["--partsize", "partsize"],
                "The number of partitions in the PartTree algorithm. Default: 50",
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            # Do not make alignment larger than number sequences. Valid only with
            # the --*parttree options. Default: the number of input sequences
            _Switch(
                ["--groupsize", "groupsize"],
                "Do not make alignment larger than number sequences. "
                "Default: the number of input sequences",
            ),
            # Adjust direction according to the first sequence
            # Mafft V6 beta function
            _Switch(
                ["--adjustdirection", "adjustdirection"],
                "Adjust direction according to the first sequence. Default off.",
            ),
            # Adjust direction according to the first sequence
            # for highly diverged data; very slow
            # Mafft V6 beta function
            _Switch(
                ["--adjustdirectionaccurately", "adjustdirectionaccurately"],
                "Adjust direction according to the first sequence,"
                "for highly diverged data; very slow"
                "Default off.",
            ),
            # **** Parameter ****
            # Gap opening penalty at group-to-group alignment. Default: 1.53
            _Option(
                ["--op", "op"],
                "Gap opening penalty at group-to-group alignment. Default: 1.53",
                checker_function=lambda x: isinstance(x, float),
                equate=False,
            ),
            # Offset value, which works like gap extension penalty, for group-to-
            # group alignment. Deafult: 0.123
            _Option(
                ["--ep", "ep"],
                "Offset value, which works like gap extension penalty, "
                "for group-to- group alignment. Default: 0.123",
                checker_function=lambda x: isinstance(x, float),
                equate=False,
            ),
            # Gap opening penalty at local pairwise alignment. Valid when the --
            # localpair or --genafpair option is selected. Default: -2.00
            _Option(
                ["--lop", "lop"],
                "Gap opening penalty at local pairwise alignment. Default: 0.123",
                checker_function=lambda x: isinstance(x, float),
                equate=False,
            ),
            # Offset value at local pairwise alignment. Valid when the --
            # localpair or --genafpair option is selected. Default: 0.1
            _Option(
                ["--lep", "lep"],
                "Offset value at local pairwise alignment. Default: 0.1",
                checker_function=lambda x: isinstance(x, float),
                equate=False,
            ),
            # Gap extension penalty at local pairwise alignment. Valid when the -
            # -localpair or --genafpair option is selected. Default: -0.1
            _Option(
                ["--lexp", "lexp"],
                "Gap extension penalty at local pairwise alignment. Default: -0.1",
                checker_function=lambda x: isinstance(x, float),
                equate=False,
            ),
            # Gap opening penalty to skip the alignment. Valid when the --
            # genafpair option is selected. Default: -6.00
            _Option(
                ["--LOP", "LOP"],
                "Gap opening penalty to skip the alignment. Default: -6.00",
                checker_function=lambda x: isinstance(x, float),
                equate=False,
            ),
            # Gap extension penalty to skip the alignment. Valid when the --
            # genafpair option is selected. Default: 0.00
            _Option(
                ["--LEXP", "LEXP"],
                "Gap extension penalty to skip the alignment. Default: 0.00",
                checker_function=lambda x: isinstance(x, float),
                equate=False,
            ),
            # BLOSUM number matrix (Henikoff and Henikoff 1992) is used.
            # number=30, 45, 62 or 80. Default: 62
            _Option(
                ["--bl", "bl"],
                "BLOSUM number matrix is used. Default: 62",
                checker_function=lambda x: x in BLOSUM_MATRICES,
                equate=False,
            ),
            # JTT PAM number (Jones et al. 1992) matrix is used. number>0.
            # Default: BLOSUM62
            _Option(
                ["--jtt", "jtt"],
                "JTT PAM number (Jones et al. 1992) matrix is used. "
                "number>0. Default: BLOSUM62",
                equate=False,
            ),
            # Transmembrane PAM number (Jones et al. 1994) matrix is used.
            # number>0. Default: BLOSUM62
            _Option(
                ["--tm", "tm"],
                "Transmembrane PAM number (Jones et al. 1994) "
                "matrix is used. number>0. Default: BLOSUM62",
                filename=True,  # to ensure spaced inputs are quoted
                equate=False,
            ),
            # Use a user-defined AA scoring matrix. The format of matrixfile is
            # the same to that of BLAST. Ignored when nucleotide sequences are
            # input. Default: BLOSUM62
            _Option(
                ["--aamatrix", "aamatrix"],
                "Use a user-defined AA scoring matrix. Default: BLOSUM62",
                filename=True,  # to ensure spaced inputs are quoted
                equate=False,
            ),
            # Incorporate the AA/nuc composition information into the scoring
            # matrix. Default: off
            _Switch(
                ["--fmodel", "fmodel"],
                "Incorporate the AA/nuc composition information into "
                "the scoring matrix (True) or not (False, default)",
            ),
            # **** Output ****
            # Name length for CLUSTAL and PHYLIP format output
            _Option(
                ["--namelength", "namelength"],
                """Name length in CLUSTAL and PHYLIP output.

                    MAFFT v6.847 (2011) added --namelength for use with
                    the --clustalout option for CLUSTAL output.

                    MAFFT v7.024 (2013) added support for this with the
                    --phylipout option for PHYLIP output (default 10).
                    """,
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            # Output format: clustal format. Default: off (fasta format)
            _Switch(
                ["--clustalout", "clustalout"],
                "Output format: clustal (True) or fasta (False, default)",
            ),
            # Output format: phylip format.
            # Added in beta with v6.847, fixed in v6.850 (2011)
            _Switch(
                ["--phylipout", "phylipout"],
                "Output format: phylip (True), or fasta (False, default)",
            ),
            # Output order: same as input. Default: on
            _Switch(
                ["--inputorder", "inputorder"],
                "Output order: same as input (True, default) or alignment "
                "based (False)",
            ),
            # Output order: aligned. Default: off (inputorder)
            _Switch(
                ["--reorder", "reorder"],
                "Output order: aligned (True) or in input order (False, default)",
            ),
            # Guide tree is output to the input.tree file. Default: off
            _Switch(
                ["--treeout", "treeout"],
                "Guide tree is output to the input.tree file (True) or "
                "not (False, default)",
            ),
            # Do not report progress. Default: off
            _Switch(
                ["--quiet", "quiet"],
                "Do not report progress (True) or not (False, default).",
            ),
            # **** Input ****
            # Assume the sequences are nucleotide. Deafult: auto
            _Switch(
                ["--nuc", "nuc"],
                "Assume the sequences are nucleotide (True/False). Default: auto",
            ),
            # Assume the sequences are amino acid. Deafult: auto
            _Switch(
                ["--amino", "amino"],
                "Assume the sequences are amino acid (True/False). Default: auto",
            ),
            # MAFFT has multiple --seed commands where the unaligned input is
            # aligned to the seed alignment. There can be multiple seeds in the
            # form: "mafft --seed align1 --seed align2 [etc] input"
            # Effectively for n number of seed alignments.
            # TODO - Can we use class _ArgumentList here?
            _Option(
                ["--seed", "seed"],
                "Seed alignments given in alignment_n (fasta format) "
                "are aligned with sequences in input.",
                filename=True,
                equate=False,
            ),
            # The input (must be FASTA format)
            _Argument(["input"], "Input file name", filename=True, is_required=True),
            # mafft-profile takes a second alignment input as an argument:
            # mafft-profile align1 align2
            _Argument(
                ["input1"],
                "Second input file name for the mafft-profile command",
                filename=True,
            ),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
