# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import warnings

from Bio import BiopythonDeprecationWarning


warnings.warn("The module Bio.ParserSupport is now obsolete, and will be "
        "deprecated and removed in a future release of Biopython.", 
        BiopythonDeprecationWarning)


from Bio.SearchIO._legacy.ParserSupport import *
