# Copyright 2002 by Jeffrey Chang.  All rights reserved.  This code is
# part of the Biopython distribution and governed by its license.
# Please see the LICENSE file that should have been included as part
# of this package.

# This module attempts to detect whether the internet is available.
# To use it, just insert "import requires_internet" into your Python
# code.  If the internet is available, then the import statement
# succeeds.  If it is not, then the statement will result in an
# MissingExternalDependencyError exception.

from Bio import MissingExternalDependencyError 

AVAILABLE = None  # Was the internet available the last time I checked?
TESTED = 0        # Have I checked before?  If so, just re-use the result.


# I'm going to check for internet availability 

RELIABLE_DOMAIN = "biopython.org"

if not TESTED:
    # Not thread-safe...
    TESTED = 1
    
    import socket
    try:
        socket.getaddrinfo(RELIABLE_DOMAIN, 80, socket.AF_UNSPEC, socket.SOCK_STREAM)
    except socket.gaierror, x:
        AVAILABLE = 0
    else:
        AVAILABLE = 1

if not AVAILABLE:
    raise MissingExternalDependencyError("internet not available")
