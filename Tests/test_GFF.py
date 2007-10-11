#!/usr/bin/env python
"""Test the Bio.GFF module and dependencies
"""
import os
from Bio import MissingExternalDependencyError

# only do the test if we are set up to do it. We need to have MYSQLPASS
# set and have a GFF wormbase installed (see the code in Bio/GFF/__init_.py
if not os.environ.has_key("MYSQLPASS"):
    raise MissingExternalDependencyError("Environment is not configured for this test (not important if you do not plan to use Bio.GFF).")

import Bio.GFF

"""
#Moved these tests to test_GFF2.py as they don't need the SQL database
import Bio.GFF.GenericTools
import Bio.GFF.easy

print "Running Bio.GFF.GenericTools doctests..."
Bio.GFF.GenericTools._test()
print "Bio.GFF.GenericTools doctests complete."

print "Running Bio.GFF.easy doctests..."
Bio.GFF.easy._test()
print "Bio.GFF.easy doctests complete."
"""

print "Running Bio.GFF doctests..."
Bio.GFF._test()
print "Bio.GFF doctests complete."
