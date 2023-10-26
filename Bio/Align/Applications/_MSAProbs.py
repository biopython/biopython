# Copyright 2013 by Christian Brueffer. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Command line wrapper for the multiple sequence alignment program MSAProbs."""

from Bio.Application import _Argument, _Option, _Switch, AbstractCommandline


class MSAProbsCommandline(AbstractCommandline):
    """Command line wrapper for MSAProbs.

    http://msaprobs.sourceforge.net

    Notes
    -----
    Last checked against version: 0.9.7

    References
    ----------
    Yongchao Liu, Bertil Schmidt, Douglas L. Maskell: "MSAProbs: multiple
    sequence alignment based on pair hidden Markov models and partition
    function posterior probabilities". Bioinformatics, 2010, 26(16): 1958 -1964

    Examples
    --------
    >>> from Bio.Align.Applications import MSAProbsCommandline
    >>> in_file = "unaligned.fasta"
    >>> out_file = "aligned.cla"
    >>> cline = MSAProbsCommandline(infile=in_file, outfile=out_file, clustalw=True)
    >>> print(cline)
    msaprobs -o aligned.cla -clustalw unaligned.fasta

    You would typically run the command line with cline() or via
    the Python subprocess module, as described in the Biopython tutorial.

    """

    def __init__(self, cmd="msaprobs", **kwargs):
        """Initialize the class."""
        # order of parameters is the same as in msaprobs -help
        self.parameters = [
            _Option(
                ["-o", "--outfile", "outfile"],
                "specify the output file name (STDOUT by default)",
                filename=True,
                equate=False,
            ),
            _Option(
                ["-num_threads", "numthreads"],
                "specify the number of threads used, and otherwise detect automatically",
                checker_function=lambda x: isinstance(x, int),
            ),
            _Switch(
                ["-clustalw", "clustalw"],
                "use CLUSTALW output format instead of FASTA format",
            ),
            _Option(
                ["-c", "consistency"],
                "use 0 <= REPS <= 5 (default: 2) passes of consistency transformation",
                checker_function=lambda x: isinstance(x, int) and 0 <= x <= 5,
            ),
            _Option(
                ["-ir", "--iterative-refinement", "iterative_refinement"],
                "use 0 <= REPS <= 1000 (default: 10) passes of iterative-refinement",
                checker_function=lambda x: isinstance(x, int) and 0 <= x <= 1000,
            ),
            _Switch(["-v", "verbose"], "report progress while aligning (default: off)"),
            _Option(
                ["-annot", "annot"],
                "write annotation for multiple alignment to FILENAME",
                filename=True,
            ),
            _Switch(
                ["-a", "--alignment-order", "alignment_order"],
                "print sequences in alignment order rather than input order (default: off)",
            ),
            _Option(["-version", "version"], "print out version of MSAPROBS"),
            _Argument(["infile"], "Multiple sequence input file", filename=True),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
