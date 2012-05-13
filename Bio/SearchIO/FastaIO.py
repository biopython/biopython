# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO parser for Bill Pearson's FASTA tools.

This module adds support for parsing FASTA outputs. FASTA is a suite of
programs that finds regions of local or global similarity between protein
or nucleotide sequences, either by searching databases or identifying
local duplications.

The following FASTA output format is supported: format 10 (triggered by the
-m 10 flag).

The current iteration of the module supports the following FASTA program:
fasta.

More information are available through these links:
  - Website: http://fasta.bioch.virginia.edu/fasta_www2/fasta_list2.shtml
  - User guide: http://fasta.bioch.virginia.edu/fasta_www2/fasta_guide.pdf


"""

from Bio.SearchIO._objects import Result, Hit, HSP


def fasta_m10_iterator(handle):
    """Generator function to parse FASTA m10 output as Result objects.

    handle -- Handle to the file, or the filename as string.

    """



def _test():
    """Run the Bio.SearchIO.FastaIO module's doctests.

    This will try and locate the unit tests directory, and run the doctests
    from there in order that the relative paths used in the examples work.
    """
    import doctest
    import os

    test_dir = 'Tests'

    if os.path.isdir(os.path.join('..', '..', test_dir)):
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join('..', '..', test_dir))
        doctest.testmod()
        os.chdir(cur_dir)
        print "Done"


# if not used as a module, run the doctest
if __name__ == "__main__":
    _test()
