# Copyright 2009 by Osvaldo Zagordi.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Command line wrapper for the short read aligner Novoalign by Novocraft."""
import types
from Bio.Application import _Option, AbstractCommandline

class NovoalignCommandline(AbstractCommandline):
    """Command line wrapper for novoalign by Novocraft.

    See www.novocraft.com - novoalign is a short read alignment program.

    Example:

    >>> from Bio.Sequencing.Applications import NovoalignCommandline
    >>> novoalign_cline = NovoalignCommandline(database='some_db',
    ...                                        readfile='some_seq.txt')
    >>> print novoalign_cline
    novoalign -d some_db -f some_seq.txt

    As will all the Biopython application wrappers, you can also add or
    change options after creating the object:

    >>> novoalign_cline.format = 'PRBnSEQ'
    >>> novoalign_cline.r_method='0.99' # limited valid values
    >>> novoalign_cline.fragment = '250 20' # must be given as a string
    >>> novoalign_cline.miRNA = 100
    >>> print novoalign_cline
    novoalign -d some_db -f some_seq.txt -F PRBnSEQ -r 0.99 -i 250 20 -m 100

    You would typically run the command line with novoalign_cline() or via
    the Python subprocess module, as described in the Biopython tutorial.

    Last checked against version: 2.05.04
    """
    def __init__(self, cmd="novoalign", **kwargs):
        
        READ_FORMAT = ['FA', 'SLXFQ', 'STDFQ', 'ILMFQ', 'PRB', 'PRBnSEQ']
        REPORT_FORMAT = ['Native', 'Pairwise', 'SAM']
        REPEAT_METHOD = ['None', 'Random', 'All', 'Exhaustive', '0.99']
        
        self.parameters = \
           [
            _Option(["-d", "database"], ["file"],
                    None, 0, "database filename",
                    0),
            _Option(["-f", "readfile"], ["file"],
                    None, 0, "read file",
                    0),
            _Option(["-F", "format"], [],
                    lambda x: x in READ_FORMAT,
                    0, "Format of read files.\n\nAllowed values: %s" % ", ".join(READ_FORMAT),
                    0),
            
            # Alignment scoring options
            _Option(["-t", "threshold"], [],
                    lambda x: isinstance(x, types.IntType),
                    0, "Threshold for alignment score",
                    0),
            _Option(["-g", "gap_open"], [],
                    lambda x: isinstance(x, types.IntType),
                    0, "Gap opening penalty [default: 40]",
                    0),
            _Option(["-x", "gap_extend"], [],
                    lambda x: isinstance(x, types.IntType),
                    0, "Gap extend penalty [default: 15]",
                    0),
            _Option(["-u", "unconverted"], [],
                    lambda x: isinstance(x, types.IntType), 0,
                    "Experimental: unconverted cytosines penalty in bisulfite mode\n\n"
                    "Default: no penalty",
                    0),
            
            # Quality control and read filtering
            _Option(["-l", "good_bases"], [],
                    lambda x: isinstance(x, types.IntType),
                    0, "Minimum number of good quality bases [default: log(N_g, 4) + 5]",
                    0),
            _Option(["-h", "homopolymer"], [],
                    lambda x: isinstance(x, types.IntType),
                    0, "Homopolymer read filter [default: 20; disable: negative value]",
                    0),
            
            # Read preprocessing options
            _Option(["-a", "adapter3"], [],
                    lambda x: isinstance(x, types.StringType),
                    0, "Strips a 3' adapter sequence prior to alignment.\n\n"
                    "With paired ends two adapters can be specified",
                    0),
            _Option(["-n", "truncate"], [],
                    lambda x: isinstance(x, types.IntType),
                    0, "Truncate to specific length before alignment",
                    0),
            _Option(["-s", "trimming"], [],
                    lambda x: isinstance(x, types.IntType),
                    0, "If fail to align, trim by s bases until they map or become shorter than l.\n\n"
                    "Ddefault: 2",
                    0),
            _Option(["-5", "adapter5"], [],
                    lambda x: isinstance(x, types.StringType),
                    0, "Strips a 5' adapter sequence.\n\n"
                    "Similar to -a (adaptor3), but on the 5' end.",
                    0),
            # Reporting options
            _Option(["-o", "report"], [],
                    lambda x: x in REPORT_FORMAT,
                    0, "Specifies the report format.\n\nAllowed values: %s\nDefault: Native" \
                    % ", ".join(REPORT_FORMAT),
                    0),
            _Option(["-Q", "quality"], [],
                    lambda x: isinstance(x, types.IntType),
                    0, "Lower threshold for an alignment to be reported [default: 0]",
                    0),
            _Option(["-R", "repeats"], [],
                    lambda x: isinstance(x, types.IntType),
                    0, "If score difference is higher, report repeats.\n\n"
                    "Otherwise -r read method applies [default: 5]",
                    0),
            _Option(["-r", "r_method"], [],
                    lambda x: x.split()[0] in REPEAT_METHOD,
                    0, "Methods to report reads with multiple matches.\n\n"
                    "Allowed values: %s\n"
                    "'All' and 'Exhaustive' accept limits." \
                    % ", ".join(REPEAT_METHOD),
                    0),
            _Option(["-e", "recorded"], [],
                    lambda x: isinstance(x, types.IntType),
                    0, "Alignments recorded with score equal to the best.\n\n"
                    "Default: 1000 in default read method, otherwise no limit.",
                    0),
            _Option(["-q", "qual_digits"], [],
                    lambda x: isinstance(x, types.IntType),
                    0, "Decimal digits for quality scores [default: 0]",
                    0),

            # Paired end options
            _Option(["-i", "fragment"], [],
                    lambda x: len(x.split()) == 2,
                    0, "Fragment length (2 reads + insert) and standard deviation [default: 250 30]",
                    0),
            _Option(["-v", "variation"], [],
                    lambda x: isinstance(x, types.IntType),
                    0, "Structural variation penalty [default: 70]",
                    0),
            
            # miRNA mode
            _Option(["-m", "miRNA"], [],
                    lambda x: isinstance(x, types.IntType),
                    0, "Sets miRNA mode and optionally sets a value for the region scanned [default: off]",
                    0),
            
            # Multithreading
            _Option(["-c", "cores"], [],
                    lambda x: isinstance(x, types.IntType),
                    0, "Number of threads, disabled on free versions [default: number of cores]",
                    0),
            
            # Quality calibrations
            _Option(["-k", "read_cal"], [],
                    lambda x: isinstance(x, types.StringType),
                    0, "Read quality calibration from file (mismatch counts)",
                    0),
            _Option(["-K", "write_cal"], [],
                    lambda x: isinstance(x, types.StringType),
                    0, "Accumulate mismatch counts and write to file",
                    0)
            ]
        AbstractCommandline.__init__(self, cmd, **kwargs)

def _test():
    """Run the module's doctests (PRIVATE)."""
    print "Runing Novoalign doctests..."
    import doctest
    doctest.testmod()
    print "Done"

if __name__ == "__main__":
    _test()
