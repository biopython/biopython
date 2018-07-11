#!/usr/bin/env python

# Copyright 2004 by Michael Hoffman.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio import MissingExternalDependencyError
import sys
if sys.platform == "win32":
    # Someone needs to find out if dnal works nicely on windows,
    # and if so where it is typically installed.
    raise MissingExternalDependencyError(
        "Don't know how to find the Wise2 tool dnal on Windows.")

from Bio._py3k import getoutput
not_found_types = ["command not found", "dnal: not found", "not recognized"]
dnal_output = getoutput("dnal")

for not_found in not_found_types:
    if not_found in dnal_output:
        raise MissingExternalDependencyError(
            "Install Wise2 (dnal) if you want to use Bio.Wise.")
