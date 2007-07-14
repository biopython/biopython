#!/usr/bin/env python
"""Test the Bio.GFF module and dependencies
"""
import os
import MySQLdb
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
# only do the test if we are set up to do it. We need to have MYSQLPASS
# set and have a GFF wormbase installed (see the code in Bio/GFF/__init_.py
if os.environ.has_key("MYSQLPASS"):
    Bio.GFF._test()
else:
    raise ImportError("Environment not configured for GFF test")
print "Bio.GFF doctests complete."
