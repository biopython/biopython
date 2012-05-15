# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO parser for HMMER output formats.

This module adds support for parsing HMMER outputs, from version 3.0 onwards.
HMMER is a suite of programs implementing the profile hidden Markov models
to find homology across protein sequences.

Specifically, this module supports the following HMMER output formats:

  - Plain text - 'hmmer-text'

And the following HMMER programs: hmmersearch, hmmerscan

More information are available through these links:
  - Web page: http://hmmer.janelia.org/
  - User guide: ftp://selab.janelia.org/pub/software/hmmer3/3.0/Userguide.pdf

"""

from Bio.SearchIO._objects import Result, Hit, HSP, SearchIndexer


def hmmer_text_iterator(handle):
    """Generator function to parse HMMER plain text output as Result objects.

    handle -- Handle to the file, or the filename as string.

    """



class HmmerTextIndexer(SearchIndexer):

    """Indexer class for HMMER plain text output."""

    def __init__(self, handle):
        pass



def _test():
    """Run the Bio.SearchIO.HmmerIO module's doctests.

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
