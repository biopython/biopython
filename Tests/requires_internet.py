# Copyright 2002 by Jeffrey Chang.  All rights reserved.  This code is
# part of the Biopython distribution and governed by its license.
# Please see the LICENSE file that should have been included as part
# of this package.

# This module attempts to detect whether the internet is available.
# To use it, import requires_internet into your Python code, and call
# requires_internet.check().  If the internet is available, then the
# import statement succeeds.  If it is not, then the statement will
# result in a MissingExternalDependencyError exception.

from Bio import MissingExternalDependencyError 

def check():
    try:
        check.available
    except AttributeError:
        # I'm going to check for internet availability 
        RELIABLE_DOMAIN = "biopython.org"
        import socket
        try:
            socket.getaddrinfo(RELIABLE_DOMAIN,
                               80,
                               socket.AF_UNSPEC,
                               socket.SOCK_STREAM)
        except socket.gaierror, x:
            check.available = False
        else:
            check.available = True
    if not check.available:
        raise MissingExternalDependencyError("internet not available")
