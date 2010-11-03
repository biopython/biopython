#!/usr/bin/env python
"""Test the Bio.GFF module and dependencies
"""
import os
from Bio import MissingExternalDependencyError, BiopythonDeprecationWarning

# only do the test if we are set up to do it. We need to have MYSQLPASS
# set and have a GFF wormbase installed (see the code in Bio/GFF/__init_.py
if "MYSQLPASS" not in os.environ:
    raise MissingExternalDependencyError("Environment is not configured for this test (not important if you do not plan to use Bio.GFF).")

import warnings
warnings.filterwarnings("ignore", category=BiopythonDeprecationWarning)
import Bio.GFF
warnings.filters.pop()

print "Running Bio.GFF doctests..."
Bio.GFF._test()
print "Bio.GFF doctests complete."
