# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO parser for Exonerate output formats.

This module supports the following Exonerate output formats:

  - Plain text alignment - 'exonerate-text'
  - Vulgar line          - 'exonerate-vulgar'
  - Cigar line           - 'exonerate-cigar'

And the following BLAST+ programs: blastn, blastp, blastx, tblastn, tblastx

More information are available through these links:
  - Home page: www.ebi.ac.uk/~guy/exonerate/

"""

# Known issues & gotchas:
# - The cigar parser does not use the extended cigar string; only supports MID
# - Cigar and vulgar parsing results will most likely be different, due to the
#   different type of data stored by both formats

from exoneratetext import *
from exoneratevulgar import *
from exoneratecigar import *
