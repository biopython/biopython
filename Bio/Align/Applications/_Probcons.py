# Copyright 2009 by Cymon J. Cox.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Command line wrapper for the multiple alignment program PROBCONS."""

from Bio.Application import _Option, _Switch, _Argument, AbstractCommandline


class ProbconsCommandline(AbstractCommandline):
    """Command line wrapper for the multiple alignment program PROBCONS.

    http://probcons.stanford.edu/

    Notes
    -----
    Last checked against version: 1.12

    References
    ----------
    Do, C.B., Mahabhashyam, M.S.P., Brudno, M., and Batzoglou, S. 2005.
    PROBCONS: Probabilistic Consistency-based Multiple Sequence Alignment.
    Genome Research 15: 330-340.

    Examples
    --------
    To align a FASTA file (unaligned.fasta) with the output in ClustalW
    format, and otherwise default settings, use:

    >>> from Bio.Align.Applications import ProbconsCommandline
    >>> probcons_cline = ProbconsCommandline(input="unaligned.fasta",
    ...                                      clustalw=True)
    >>> print(probcons_cline)
    probcons -clustalw unaligned.fasta

    You would typically run the command line with probcons_cline() or via
    the Python subprocess module, as described in the Biopython tutorial.

    Note that PROBCONS will write the alignment to stdout, which you may
    want to save to a file and then parse, e.g.::

        stdout, stderr = probcons_cline()
        with open("aligned.aln", "w") as handle:
            handle.write(stdout)
        from Bio import AlignIO
        align = AlignIO.read("aligned.fasta", "clustalw")

    Alternatively, to parse the output with AlignIO directly you can
    use StringIO to turn the string into a handle::

        stdout, stderr = probcons_cline()
        from io import StringIO
        from Bio import AlignIO
        align = AlignIO.read(StringIO(stdout), "clustalw")

    """

    def __init__(self, cmd="probcons", **kwargs):
        """Initialize the class."""
        self.parameters = [
            # Note that some options cannot be assigned via properties using the
            # original documented option (because hyphens are not valid for names in
            # python), e.g cmdline.pre-training = 3 will not work
            # In these cases the shortened option name should be used
            # cmdline.pre = 3
            _Switch(
                ["-clustalw", "clustalw"], "Use CLUSTALW output format instead of MFA"
            ),
            _Option(
                ["-c", "c", "--consistency", "consistency"],
                "Use 0 <= REPS <= 5 (default: 2) passes of consistency transformation",
                checker_function=lambda x: x in range(6),
                equate=False,
            ),
            _Option(
                ["-ir", "--iterative-refinement", "iterative-refinement", "ir"],
                "Use 0 <= REPS <= 1000 (default: 100) passes of iterative-refinement",
                checker_function=lambda x: x in range(1001),
                equate=False,
            ),
            _Option(
                ["-pre", "--pre-training", "pre-training", "pre"],
                "Use 0 <= REPS <= 20 (default: 0) rounds of pretraining",
                checker_function=lambda x: x in range(21),
                equate=False,
            ),
            _Switch(["-pairs", "pairs"], "Generate all-pairs pairwise alignments"),
            _Switch(
                ["-viterbi", "viterbi"],
                "Use Viterbi algorithm to generate all pairs "
                "(automatically enables -pairs)",
            ),
            _Switch(
                ["-verbose", "verbose"], "Report progress while aligning (default: off)"
            ),
            _Option(
                ["-annot", "annot"],
                "Write annotation for multiple alignment to FILENAME",
                equate=False,
            ),
            _Option(
                ["-t", "t", "--train", "train"],
                "Compute EM transition probabilities, store in FILENAME "
                "(default: no training)",
                equate=False,
            ),
            _Switch(
                ["-e", "e", "--emissions", "emissions"],
                "Also reestimate emission probabilities (default: off)",
            ),
            _Option(
                ["-p", "p", "--paramfile", "paramfile"],
                "Read parameters from FILENAME",
                equate=False,
            ),
            _Switch(
                ["-a", "--alignment-order", "alignment-order", "a"],
                "Print sequences in alignment order rather than input "
                "order (default: off)",
            ),
            # Input file name
            _Argument(
                ["input"],
                "Input file name. Must be multiple FASTA alignment (MFA) format",
                filename=True,
                is_required=True,
            ),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
