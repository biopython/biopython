# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Custom indexing for Bio.SearchIO objects (PRIVATE).

"""

try:
    from sqlite3 import dbapi2 as sqlite
    from sqlite3 import IntegrityError, OperationError
except ImportError:
    sqlite = None
