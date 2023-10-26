# Copyright 2009 by Cymon J. Cox.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Command line wrapper for the multiple alignment program DIALIGN2-2."""

from Bio.Application import _Option, _Argument, _Switch, AbstractCommandline


class DialignCommandline(AbstractCommandline):
    """Command line wrapper for the multiple alignment program DIALIGN2-2.

    http://bibiserv.techfak.uni-bielefeld.de/dialign/welcome.html

    Notes
    -----
    Last checked against version: 2.2

    References
    ----------
    B. Morgenstern (2004). DIALIGN: Multiple DNA and Protein Sequence
    Alignment at BiBiServ. Nucleic Acids Research 32, W33-W36.

    Examples
    --------
    To align a FASTA file (unaligned.fasta) with the output files names
    aligned.* including a FASTA output file (aligned.fa), use:

    >>> from Bio.Align.Applications import DialignCommandline
    >>> dialign_cline = DialignCommandline(input="unaligned.fasta",
    ...                                    fn="aligned", fa=True)
    >>> print(dialign_cline)
    dialign2-2 -fa -fn aligned unaligned.fasta

    You would typically run the command line with dialign_cline() or via
    the Python subprocess module, as described in the Biopython tutorial.

    """

    def __init__(self, cmd="dialign2-2", **kwargs):
        """Initialize the class."""
        self.program_name = cmd
        self.parameters = [
            _Switch(
                ["-afc", "afc"],
                r"Creates additional output file '\*.afc' "
                "containing data of all fragments considered "
                "for alignment WARNING: this file can be HUGE !",
            ),
            _Switch(
                ["-afc_v", "afc_v"],
                "Like '-afc' but verbose: fragments are explicitly "
                "printed. WARNING: this file can be EVEN BIGGER !",
            ),
            _Switch(
                ["-anc", "anc"],
                "Anchored alignment. Requires a file <seq_file>.anc "
                "containing anchor points.",
            ),
            _Switch(
                ["-cs", "cs"],
                "If segments are translated, not only the 'Watson "
                "strand' but also the 'Crick strand' is looked at.",
            ),
            _Switch(["-cw", "cw"], "Additional output file in CLUSTAL W format."),
            _Switch(
                ["-ds", "ds"],
                "'dna alignment speed up' - non-translated nucleic acid "
                "fragments are taken into account only if they start "
                "with at least two matches. Speeds up DNA alignment at "
                "the expense of sensitivity.",
            ),
            _Switch(["-fa", "fa"], "Additional output file in FASTA format."),
            _Switch(
                ["-ff", "ff"],
                r"Creates file \*.frg containing information about all "
                "fragments that are part of the respective optimal "
                "pairwise alignmnets plus information about "
                "consistency in the multiple alignment",
            ),
            _Option(
                ["-fn", "fn"],
                "Output files are named <out_file>.<extension>.",
                equate=False,
            ),
            _Switch(
                ["-fop", "fop"],
                r"Creates file \*.fop containing coordinates of all "
                "fragments that are part of the respective pairwise alignments.",
            ),
            _Switch(
                ["-fsm", "fsm"],
                r"Creates file \*.fsm containing coordinates of all "
                "fragments that are part of the final alignment",
            ),
            _Switch(
                ["-iw", "iw"],
                "Overlap weights switched off (by default, overlap "
                "weights are used if up to 35 sequences are aligned). "
                "This option speeds up the alignment but may lead "
                "to reduced alignment quality.",
            ),
            _Switch(
                ["-lgs", "lgs"],
                "'long genomic sequences' - combines the following "
                "options: -ma, -thr 2, -lmax 30, -smin 8, -nta, -ff, "
                "-fop, -ff, -cs, -ds, -pst ",
            ),
            _Switch(
                ["-lgs_t", "lgs_t"],
                "Like '-lgs' but with all segment pairs assessed "
                "at the peptide level (rather than 'mixed alignments' "
                "as with the '-lgs' option). Therefore faster than "
                "-lgs but not very sensitive for non-coding regions.",
            ),
            _Option(
                ["-lmax", "lmax"],
                "Maximum fragment length = x  (default: x = 40 or "
                "x = 120 for 'translated' fragments). Shorter x "
                "speeds up the program but may affect alignment quality.",
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            _Switch(
                ["-lo", "lo"],
                r"(Long Output) Additional file \*.log with information "
                "about fragments selected for pairwise alignment and "
                "about consistency in multi-alignment procedure.",
            ),
            _Switch(
                ["-ma", "ma"],
                "'mixed alignments' consisting of P-fragments and "
                "N-fragments if nucleic acid sequences are aligned.",
            ),
            _Switch(
                ["-mask", "mask"],
                "Residues not belonging to selected fragments are "
                r"replaced by '\*' characters in output alignment "
                "(rather than being printed in lower-case characters)",
            ),
            _Switch(
                ["-mat", "mat"],
                r"Creates file \*mat with substitution counts derived "
                "from the fragments that have been selected for alignment.",
            ),
            _Switch(
                ["-mat_thr", "mat_thr"],
                "Like '-mat' but only fragments with weight score "
                "> t are considered",
            ),
            _Switch(
                ["-max_link", "max_link"],
                "'maximum linkage' clustering used to construct "
                "sequence tree (instead of UPGMA).",
            ),
            _Switch(["-min_link", "min_link"], "'minimum linkage' clustering used."),
            _Option(["-mot", "mot"], "'motif' option.", equate=False),
            _Switch(["-msf", "msf"], "Separate output file in MSF format."),
            _Switch(
                ["-n", "n"],
                "Input sequences are nucleic acid sequences. "
                "No translation of fragments.",
            ),
            _Switch(
                ["-nt", "nt"],
                "Input sequences are nucleic acid sequences and "
                "'nucleic acid segments' are translated to 'peptide "
                "segments'.",
            ),
            _Switch(
                ["-nta", "nta"],
                "'no textual alignment' - textual alignment suppressed. "
                "This option makes sense if other output files are of "
                "interest -- e.g. the fragment files created with -ff, "
                "-fop, -fsm or -lo.",
            ),
            _Switch(
                ["-o", "o"],
                "Fast version, resulting alignments may be slightly different.",
            ),
            _Switch(
                ["-ow", "ow"],
                "Overlap weights enforced (By default, overlap weights "
                "are used only if up to 35 sequences are aligned since "
                "calculating overlap weights is time consuming).",
            ),
            _Switch(
                ["-pst", "pst"],
                r"'print status'. Creates and updates a file \*.sta with "
                "information about the current status of the program "
                "run.  This option is recommended if large data sets "
                "are aligned since it allows the user to estimate the "
                "remaining running time.",
            ),
            _Switch(
                ["-smin", "smin"],
                "Minimum similarity value for first residue pair "
                "(or codon pair) in fragments. Speeds up protein "
                "alignment or alignment of translated DNA fragments "
                "at the expense of sensitivity.",
            ),
            _Option(
                ["-stars", "stars"],
                r"Maximum number of '\*' characters indicating degree "
                "of local similarity among sequences. By default, no "
                "stars are used but numbers between 0 and 9, instead.",
                checker_function=lambda x: x in range(0, 10),
                equate=False,
            ),
            _Switch(["-stdo", "stdo"], "Results written to standard output."),
            _Switch(
                ["-ta", "ta"],
                "Standard textual alignment printed (overrides "
                "suppression of textual alignments in special "
                "options, e.g. -lgs)",
            ),
            _Option(
                ["-thr", "thr"],
                "Threshold T = x.",
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            _Switch(
                ["-xfr", "xfr"],
                "'exclude fragments' - list of fragments can be "
                "specified that are NOT considered for pairwise alignment",
            ),
            _Argument(
                ["input"],
                "Input file name. Must be FASTA format",
                filename=True,
                is_required=True,
            ),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
