# Copyright 1999-2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
# Patches by Mike Poidinger to support multiple databases.
# Updated by Peter Cock in 2007 to do a better job on BLAST 2.2.15
# Migrated to Bio.SearchIO._legacy in 2013 by Wibowo Arindrarto for deprecation

import warnings

from Bio import BiopythonDeprecationWarning


warnings.warn("The module Bio.Blast.NCBIStandalone is now deprecated, and will "
        "be removed in a future Biopython release. To parse plain text BLAST "
        "output, please use Bio.SearchIO instead.", BiopythonDeprecationWarning)


from Bio.SearchIO._legacy.NCBIStandalone import *
