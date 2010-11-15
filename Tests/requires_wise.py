#!/usr/bin/env python

# Copyright 2004 by Michael Hoffman.  All rights reserved.  This code is
# part of the Biopython distribution and governed by its license.
# Please see the LICENSE file that should have been included as part
# of this package.

from Bio import MissingExternalDependencyError
import sys
if sys.platform=="win32":
    #Someone needs to find out if dnal works nicely on windows,
    #and if so where it is typically installed.
    raise MissingExternalDependencyError(\
        "Don't know how to find the Wise2 tool dnal on Windows.")

import commands
not_found_types = ["command not found", "dnal: not found", "not recognized"]
dnal_output = commands.getoutput("dnal")

for not_found in not_found_types:
    if dnal_output.find(not_found) != -1:
        #raise MissingExternalDependencyError(dnal_output)
        raise MissingExternalDependencyError(\
            "Install Wise2 (dnal) if you want to use Bio.Wise.")
