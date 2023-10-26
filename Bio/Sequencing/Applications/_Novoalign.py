# Copyright 2009 by Osvaldo Zagordi.  All rights reserved.
# Revisions copyright 2010 by Peter Cock.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Command line wrapper for the short read aligner Novoalign by Novocraft."""


from Bio.Application import _Option, AbstractCommandline


class NovoalignCommandline(AbstractCommandline):
    """Command line wrapper for novoalign by Novocraft.

    See www.novocraft.com - novoalign is a short read alignment program.

    Examples
    --------
    >>> from Bio.Sequencing.Applications import NovoalignCommandline
    >>> novoalign_cline = NovoalignCommandline(database='some_db',
    ...                                        readfile='some_seq.txt')
    >>> print(novoalign_cline)
    novoalign -d some_db -f some_seq.txt

    As with all the Biopython application wrappers, you can also add or
    change options after creating the object:

    >>> novoalign_cline.format = 'PRBnSEQ'
    >>> novoalign_cline.r_method='0.99' # limited valid values
    >>> novoalign_cline.fragment = '250 20' # must be given as a string
    >>> novoalign_cline.miRNA = 100
    >>> print(novoalign_cline)
    novoalign -d some_db -f some_seq.txt -F PRBnSEQ -r 0.99 -i 250 20 -m 100

    You would typically run the command line with novoalign_cline() or via
    the Python subprocess module, as described in the Biopython tutorial.

    Last checked against version: 2.05.04

    """

    def __init__(self, cmd="novoalign", **kwargs):
        """Initialize the class."""
        READ_FORMAT = ["FA", "SLXFQ", "STDFQ", "ILMFQ", "PRB", "PRBnSEQ"]
        REPORT_FORMAT = ["Native", "Pairwise", "SAM"]
        REPEAT_METHOD = ["None", "Random", "All", "Exhaustive", "0.99"]

        self.parameters = [
            _Option(
                ["-d", "database"], "database filename", filename=True, equate=False
            ),
            _Option(["-f", "readfile"], "read file", filename=True, equate=False),
            _Option(
                ["-F", "format"],
                f"Format of read files.\n\nAllowed values: {', '.join(READ_FORMAT)}",
                checker_function=lambda x: x in READ_FORMAT,
                equate=False,
            ),
            # Alignment scoring options
            _Option(
                ["-t", "threshold"],
                "Threshold for alignment score",
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            _Option(
                ["-g", "gap_open"],
                "Gap opening penalty [default: 40]",
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            _Option(
                ["-x", "gap_extend"],
                "Gap extend penalty [default: 15]",
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            _Option(
                ["-u", "unconverted"],
                "Experimental: unconverted cytosines penalty in bisulfite mode\n\n"
                "Default: no penalty",
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            # Quality control and read filtering
            _Option(
                ["-l", "good_bases"],
                "Minimum number of good quality bases [default: log(N_g, 4) + 5]",
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            _Option(
                ["-h", "homopolymer"],
                "Homopolymer read filter [default: 20; disable: negative value]",
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            # Read preprocessing options
            _Option(
                ["-a", "adapter3"],
                "Strips a 3' adapter sequence prior to alignment.\n\n"
                "With paired ends two adapters can be specified",
                checker_function=lambda x: isinstance(x, str),
                equate=False,
            ),
            _Option(
                ["-n", "truncate"],
                "Truncate to specific length before alignment",
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            _Option(
                ["-s", "trimming"],
                "If fail to align, trim by s bases until they map or become shorter than l.\n\n"
                "Ddefault: 2",
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            _Option(
                ["-5", "adapter5"],
                "Strips a 5' adapter sequence.\n\n"
                "Similar to -a (adaptor3), but on the 5' end.",
                checker_function=lambda x: isinstance(x, str),
                equate=False,
            ),
            # Reporting options
            _Option(
                ["-o", "report"],
                "Specifies the report format.\n\nAllowed values: %s\nDefault: Native"
                % ", ".join(REPORT_FORMAT),
                checker_function=lambda x: x in REPORT_FORMAT,
                equate=False,
            ),
            _Option(
                ["-Q", "quality"],
                "Lower threshold for an alignment to be reported [default: 0]",
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            _Option(
                ["-R", "repeats"],
                "If score difference is higher, report repeats.\n\n"
                "Otherwise -r read method applies [default: 5]",
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            _Option(
                ["-r", "r_method"],
                "Methods to report reads with multiple matches.\n\n"
                "Allowed values: %s\n"
                "'All' and 'Exhaustive' accept limits." % ", ".join(REPEAT_METHOD),
                checker_function=lambda x: x.split()[0] in REPEAT_METHOD,
                equate=False,
            ),
            _Option(
                ["-e", "recorded"],
                "Alignments recorded with score equal to the best.\n\n"
                "Default: 1000 in default read method, otherwise no limit.",
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            _Option(
                ["-q", "qual_digits"],
                "Decimal digits for quality scores [default: 0]",
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            # Paired end options
            _Option(
                ["-i", "fragment"],
                "Fragment length (2 reads + insert) and standard deviation [default: 250 30]",
                checker_function=lambda x: len(x.split()) == 2,
                equate=False,
            ),
            _Option(
                ["-v", "variation"],
                "Structural variation penalty [default: 70]",
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            # miRNA mode
            _Option(
                ["-m", "miRNA"],
                "Sets miRNA mode and optionally sets a value for the region scanned [default: off]",
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            # Multithreading
            _Option(
                ["-c", "cores"],
                "Number of threads, disabled on free versions [default: number of cores]",
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            # Quality calibrations
            _Option(
                ["-k", "read_cal"],
                "Read quality calibration from file (mismatch counts)",
                checker_function=lambda x: isinstance(x, str),
                equate=False,
            ),
            _Option(
                ["-K", "write_cal"],
                "Accumulate mismatch counts and write to file",
                checker_function=lambda x: isinstance(x, str),
                equate=False,
            ),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
