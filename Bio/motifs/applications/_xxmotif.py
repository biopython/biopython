# Copyright 2012 by Christian Brueffer.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Command line wrapper for the motif finding program XXmotif."""


import os
from Bio.Application import AbstractCommandline, _Option, _Switch, _Argument


class XXmotifCommandline(AbstractCommandline):
    """Command line wrapper for XXmotif.

    http://xxmotif.genzentrum.lmu.de/

    Notes
    -----
    Last checked against version: 1.3

    References
    ----------
    Luehr S, Hartmann H, and Söding J. The XXmotif web server for eXhaustive,
    weight matriX-based motif discovery in nucleotide sequences,
    Nucleic Acids Res. 40: W104-W109 (2012).

    Hartmann H, Guthoehrlein EW, Siebert M., Luehr S, and Söding J. P-value
    based regulatory motif discovery using positional weight matrices,
    Genome Res. 23: 181–194 (2013)

    Examples
    --------
    >>> from Bio.motifs.applications import XXmotifCommandline
    >>> out_dir = "results"
    >>> in_file = "sequences.fasta"
    >>> xxmotif_cline = XXmotifCommandline(outdir=out_dir, seqfile=in_file, revcomp=True)
    >>> print(xxmotif_cline)
    XXmotif results sequences.fasta --revcomp

    You would typically run the command line with xxmotif_cline() or via
    the Python subprocess module, as described in the Biopython tutorial.

    """

    def __init__(self, cmd="XXmotif", **kwargs):
        """Initialize the class."""
        # order of parameters is the same as in XXmotif --help
        _valid_alphabet = set("ACGTNX")

        self.parameters = [
            _Argument(
                ["outdir", "OUTDIR"],
                "output directory for all results",
                filename=True,
                is_required=True,
                # XXmotif currently does not accept spaces in the outdir name
                checker_function=lambda x: " " not in x,
            ),
            _Argument(
                ["seqfile", "SEQFILE"],
                "file name with sequences from positive set in FASTA format",
                filename=True,
                is_required=True,
                # XXmotif currently only accepts a pure filename
                checker_function=lambda x: os.path.split(x)[0] == "",
            ),
            # Options
            _Option(
                ["--negSet", "negSet", "NEGSET", "negset"],
                "sequence set which has to be used as a reference set",
                filename=True,
                equate=False,
            ),
            _Switch(
                ["--zoops", "ZOOPS", "zoops"],
                "use zero-or-one occurrence per sequence model (DEFAULT)",
            ),
            _Switch(
                ["--mops", "MOPS", "mops"], "use multiple occurrence per sequence model"
            ),
            _Switch(
                ["--oops", "OOPS", "oops"], "use one occurrence per sequence model"
            ),
            _Switch(
                ["--revcomp", "REVCOMP", "revcomp"],
                "search in reverse complement of sequences as well (DEFAULT: NO)",
            ),
            _Option(
                [
                    "--background-model-order",
                    "background-model-order",
                    "BACKGROUND-MODEL-ORDER",
                    "background_model_order",
                ],
                "order of background distribution (DEFAULT: 2, 8(--negset) )",
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            _Option(
                ["--pseudo", "PSEUDO", "pseudo"],
                "percentage of pseudocounts used (DEFAULT: 10)",
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            _Option(
                ["-g", "--gaps", "GAPS", "gaps"],
                "maximum number of gaps used for start seeds [0-3] (DEFAULT: 0)",
                checker_function=lambda x: x in [0 - 3],
                equate=False,
            ),
            _Option(
                ["--type", "TYPE", "type"],
                "defines what kind of start seeds are used (DEFAULT: ALL)"
                "possible types: ALL, FIVEMERS, PALINDROME, TANDEM, NOPALINDROME, NOTANDEM",
                checker_function=lambda x: x
                in [
                    "ALL",
                    "all",
                    "FIVEMERS",
                    "fivemers",
                    "PALINDROME",
                    "palindrome",
                    "TANDEM",
                    "tandem",
                    "NOPALINDROME",
                    "nopalindrome",
                    "NOTANDEM",
                    "notandem",
                ],
                equate=False,
            ),
            _Option(
                [
                    "--merge-motif-threshold",
                    "merge-motif-threshold",
                    "MERGE-MOTIF-THRESHOLD",
                    "merge_motif_threshold",
                ],
                "defines the similarity threshold for merging motifs (DEFAULT: HIGH)"
                "possible modes: LOW, MEDIUM, HIGH",
                checker_function=lambda x: x
                in ["LOW", "low", "MEDIUM", "medium", "HIGH", "high"],
                equate=False,
            ),
            _Switch(
                [
                    "--no-pwm-length-optimization",
                    "no-pwm-length-optimization",
                    "NO-PWM-LENGTH-OPTIMIZATION",
                    "no_pwm_length_optimization",
                ],
                "do not optimize length during iterations (runtime advantages)",
            ),
            _Option(
                [
                    "--max-match-positions",
                    "max-match-positions",
                    "MAX-MATCH-POSITIONS",
                    "max_match_positions",
                ],
                "max number of positions per motif (DEFAULT: 17, higher values will lead to very long runtimes)",
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            _Switch(
                ["--batch", "BATCH", "batch"],
                "suppress progress bars (reduce output size for batch jobs)",
            ),
            _Option(
                ["--maxPosSetSize", "maxPosSetSize", "MAXPOSSETSIZE", "maxpossetsize"],
                "maximum number of sequences from the positive set used [DEFAULT: all]",
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            # does not make sense in biopython
            # _Switch(["--help", "help", "HELP"],
            #         "print this help page"),
            _Option(
                ["--trackedMotif", "trackedMotif", "TRACKEDMOTIF", "trackedmotif"],
                "inspect extensions and refinement of a given seed (DEFAULT: not used)",
                checker_function=lambda x: any((c in _valid_alphabet) for c in x),
                equate=False,
            ),
            # Using conservation information
            _Option(
                ["--format", "FORMAT", "format"],
                "defines what kind of format the input sequences have (DEFAULT: FASTA)",
                checker_function=lambda x: x in ["FASTA", "fasta", "MFASTA", "mfasta"],
                equate=False,
            ),
            _Option(
                [
                    "--maxMultipleSequences",
                    "maxMultipleSequences",
                    "MAXMULTIPLESEQUENCES",
                    "maxmultiplesequences",
                ],
                "maximum number of sequences used in an alignment [DEFAULT: all]",
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            # Using localization information
            _Switch(
                ["--localization", "LOCALIZATION", "localization"],
                "use localization information to calculate combined P-values"
                "(sequences should have all the same length)",
            ),
            _Option(
                ["--downstream", "DOWNSTREAM", "downstream"],
                "number of residues in positive set downstream of anchor point (DEFAULT: 0)",
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            # Start with self defined motif
            _Option(
                ["-m", "--startMotif", "startMotif", "STARTMOTIF", "startmotif"],
                "Start motif (IUPAC characters)",
                checker_function=lambda x: any((c in _valid_alphabet) for c in x),
                equate=False,
            ),
            _Option(
                ["-p", "--profileFile", "profileFile", "PROFILEFILE", "profilefile"],
                "profile file",
                filename=True,
                equate=False,
            ),
            _Option(
                ["--startRegion", "startRegion", "STARTREGION", "startregion"],
                "expected start position for motif occurrences relative to anchor point (--localization)",
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            _Option(
                ["--endRegion", "endRegion", "ENDREGION", "endregion"],
                "expected end position for motif occurrences relative to anchor point (--localization)",
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            # XXmotif wrapper options
            _Switch(
                ["--XXmasker", "masker"],
                "mask the input sequences for homology, repeats and low complexity regions",
            ),
            _Switch(
                ["--XXmasker-pos", "maskerpos"],
                "mask only the positive set for homology, repeats and low complexity regions",
            ),
            _Switch(
                ["--no-graphics", "nographics"], "run XXmotif without graphical output"
            ),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
