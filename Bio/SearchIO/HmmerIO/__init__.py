# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO parser for HMMER output formats.

This module adds support for parsing HMMER outputs, from version 3.0 onwards.
HMMER is a suite of programs implementing the profile hidden Markov models
to find homology across protein sequences.

Specifically, this module supports the following HMMER output formats:

  - Plain text - 'hmmer3-text'
  - Table - 'hmmer3-tab'
  - Domain table - 'hmmscan-domtab', 'hmmsearch-domtab', or 'phmmer3-domtab'

And the following HMMER programs: hmmersearch, hmmerscan, phmmer

More information are available through these links:
  - Web page: http://hmmer.janelia.org/
  - User guide: ftp://selab.janelia.org/pub/software/hmmer3/3.0/Userguide.pdf

"""

from hmmer3_domtab import *
from hmmer3_text import *
from hmmer3_tab import *


# if not used as a module, run the doctest
if __name__ == "__main__":
    from Bio.SearchIO._utils import run_doctest
    run_doctest()
