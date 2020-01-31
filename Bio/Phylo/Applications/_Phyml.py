# Copyright 2011 by Eric Talevich.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Command-line wrapper for the tree inference program PhyML."""

from Bio.Application import _Option, _Switch, AbstractCommandline


class PhymlCommandline(AbstractCommandline):
    """Command-line wrapper for the tree inference program PhyML.

    Homepage: http://www.atgc-montpellier.fr/phyml

    References
    ----------
    Guindon S, Gascuel O.
    A simple, fast, and accurate algorithm to estimate large phylogenies by maximum
    likelihood.
    Systematic Biology, 2003 Oct;52(5):696-704.
    PubMed PMID: 14530136.

    Guindon S, Dufayard JF, Lefort V, Anisimova M, Hordijk W, Gascuel O.
    New Algorithms and Methods to Estimate Maximum-Likelihood Phylogenies: Assessing
    the Performance of PhyML 3.0.
    Systematic Biology, 2010 59(3):307-21.

    """

    def __init__(self, cmd="phyml", **kwargs):
        """Initialize the class."""
        self.parameters = [
            _Option(
                ["-i", "--input", "input"],
                "PHYLIP format input nucleotide or amino-acid sequence filenam.",
                filename=True,
                is_required=True,
                equate=False,
            ),
            _Option(
                ["-d", "--datatype", "datatype"],
                "Datatype 'nt' for nucleotide (default) or 'aa' for amino-acids.",
                checker_function=lambda x: x in ("nt", "aa"),
                equate=False,
            ),
            _Switch(
                ["-q", "--sequential", "sequential"],
                "Changes interleaved format (default) to sequential format.",
            ),
            _Option(
                ["-n", "--multiple", "multiple"],
                "Number of data sets to analyse (integer).",
                checker_function=(lambda x: isinstance(x, int) or x.isdigit()),
                equate=False,
            ),
            _Switch(
                ["-p", "--pars", "pars"],
                """Use a minimum parsimony starting tree.

                    This option is taken into account when the '-u' option is absent
                    and when tree topology modifications are to be done.
                    """,
            ),
            _Option(
                ["-b", "--bootstrap", "bootstrap"],
                r"""Number of bootstrap replicates, if value is > 0.

                    Otherwise:

                    0: neither approximate likelihood ratio test nor bootstrap
                    values are computed.

                    -1: approximate likelihood ratio test returning aLRT statistics.

                    -2: approximate likelihood ratio test returning Chi2-based
                    parametric branch supports.

                    -4: SH-like branch supports alone.
                    """,
                equate=False,
            ),
            _Option(
                ["-m", "--model", "model"],
                """Substitution model name.

                    Nucleotide-based models:

                    HKY85 (default) | JC69 | K80 | F81 | F84 | TN93 | GTR | custom

                    For the custom option, a string of six digits identifies the
                    model. For instance, 000000 corresponds to F81 (or JC69,
                    provided the distribution of nucleotide frequencies is uniform).
                    012345 corresponds to GTR. This option can be used for encoding
                    any model that is a nested within GTR.

                    Amino-acid based models:

                    LG (default) | WAG | JTT | MtREV | Dayhoff | DCMut | RtREV |
                    CpREV | VT | Blosum62 | MtMam | MtArt | HIVw | HIVb | custom
                    """,
                checker_function=(
                    lambda x: x
                    in (
                        # Nucleotide models:
                        "HKY85",
                        "JC69",
                        "K80",
                        "F81",
                        "F84",
                        "TN93",
                        "GTR",
                        # Amino acid models:
                        "LG",
                        "WAG",
                        "JTT",
                        "MtREV",
                        "Dayhoff",
                        "DCMut",
                        "RtREV",
                        "CpREV",
                        "VT",
                        "Blosum62",
                        "MtMam",
                        "MtArt",
                        "HIVw",
                        "HIVb",
                    )
                    or isinstance(x, int)
                ),
                equate=False,
            ),
            _Option(
                ["-f", "frequencies"],
                """Character frequencies.

                    -f e, m, or "fA fC fG fT"

                    e : Empirical frequencies, determined as follows :

                        - Nucleotide sequences: (Empirical) the equilibrium base
                          frequencies are estimated by counting the occurrence
                          of the different bases in the alignment.
                        - Amino-acid sequences: (Empirical) the equilibrium
                          amino-acid frequencies are estimated by counting the
                          occurrence of the different amino-acids in the alignment.

                    m : ML/model-based frequencies, determined as follows :

                        - Nucleotide sequences: (ML) the equilibrium base
                          frequencies are estimated using maximum likelihood
                        - Amino-acid sequences: (Model) the equilibrium amino-acid
                          frequencies are estimated using the frequencies defined by
                          the substitution model.

                    "fA fC fG fT" : only valid for nucleotide-based models.
                    fA, fC, fG and fT are floating-point numbers that correspond
                    to the frequencies of A, C, G and T, respectively.
                    """,
                filename=True,  # ensure ".25 .25 .25 .25" stays quoted
                equate=False,
            ),
            _Option(
                ["-t", "--ts/tv", "ts_tv_ratio"],
                """Transition/transversion ratio. (DNA sequences only.)

                    Can be a fixed positive value (ex:4.0) or e to get the
                    maximum-likelihood estimate.
                    """,
                equate=False,
            ),
            _Option(
                ["-v", "--pinv", "prop_invar"],
                """Proportion of invariable sites.

                    Can be a fixed value in the range [0,1], or 'e' to get the
                    maximum-likelihood estimate.
                    """,
                equate=False,
            ),
            _Option(
                ["-c", "--nclasses", "nclasses"],
                """Number of relative substitution rate categories.

                    Default 1. Must be a positive integer.
                    """,
                equate=False,
            ),
            _Option(
                ["-a", "--alpha", "alpha"],
                """Distribution of the gamma distribution shape parameter.

                    Can be a fixed positive value, or 'e' to get the
                    maximum-likelihood estimate.
                    """,
                equate=False,
            ),
            _Option(
                ["-s", "--search", "search"],
                """Tree topology search operation option.

                    Can be one of:

                        NNI : default, fast

                        SPR : a bit slower than NNI

                        BEST : best of NNI and SPR search
                    """,
                checker_function=lambda x: x in ("NNI", "SPR", "BEST"),
                equate=False,
            ),
            # alt name: user_tree_file
            _Option(
                ["-u", "--inputtree", "input_tree"],
                "Starting tree filename. The tree must be in Newick format.",
                filename=True,
                equate=False,
            ),
            _Option(
                ["-o", "optimize"],
                r"""Specific parameter optimisation.

                    tlr : tree topology (t), branch length (l) and
                    rate parameters (r) are optimised.

                    tl  : tree topology and branch length are optimised.

                    lr  : branch length and rate parameters are optimised.

                    l   : branch length are optimised.

                    r   : rate parameters are optimised.

                    n   : no parameter is optimised.
                    """,
                equate=False,
            ),
            _Switch(
                ["--rand_start", "rand_start"],
                """Sets the initial tree to random.

                    Only valid if SPR searches are to be performed.
                    """,
            ),
            _Option(
                ["--n_rand_starts", "n_rand_starts"],
                """Number of initial random trees to be used.

                    Only valid if SPR searches are to be performed.
                    """,
                equate=False,
            ),
            _Option(
                ["--r_seed", "r_seed"],
                """Seed used to initiate the random number generator.

                    Must be an integer.
                    """,
                equate=False,
            ),
            _Switch(
                ["--print_site_lnl", "print_site_lnl"],
                r"Print the likelihood for each site in file \*_phyml_lk.txt.",
            ),
            _Switch(
                ["--print_trace", "print_trace"],
                r"""
                    Print each phylogeny explored during the tree search process
                    in file \*_phyml_trace.txt.""",
            ),
            _Option(
                ["--run_id", "run_id"],
                """Append the given string at the end of each PhyML output file.

                    This option may be useful when running simulations involving
                    PhyML.
                    """,
                checker_function=lambda x: isinstance(x, str),
                equate=False,
            ),
            # XXX should this always be set to True?
            _Switch(
                ["--quiet", "quiet"],
                "No interactive questions (for running in batch mode).",
            ),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)
