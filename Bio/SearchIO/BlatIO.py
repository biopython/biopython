# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO parser for BLAT output formats.

This module adds support for parsing BLAT outputs. BLAT (BLAST-Like Alignment
Tool) is a sequence homology search program initially built for annotating
the human genome.

Specifically, this module supports the following BLAT output formats:

  - PSL - 'blat-psl'

More information are available through these links:
  
  - Publication: http://genome.cshlp.org/content/12/4/656
  - User guide: http://genome.ucsc.edu/goldenPath/help/blatSpec.html
  - Source download: http://www.soe.ucsc.edu/~kent/src
  - Executable download: http://hgdownload.cse.ucsc.edu/admin/exe/

"""

from Bio.SearchIO._objects import Result, Hit, HSP
from Bio.SearchIO._index import SearchIndexer


def blat_psl_iterator(handle):
    """Generator function to parse BLAT PSL output as Result objects.

    handle -- Handle to the file, or the filename as string.

    """



class BlatPslIndexer(SearchIndexer):

    """Indexer class for BLAT PSL output."""

    def __init__(self, handle):
        pass



def _test():
    """Run the Bio.SearchIO.BlatIO module's doctests.

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
