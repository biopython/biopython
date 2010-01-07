#!/usr/bin/env python
"""Test the Bio.GFF dependencies
"""

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
import Bio.GFF
import Bio.GFF.GenericTools
import Bio.GFF.easy
warnings.resetwarnings()

print "Running Bio.GFF.GenericTools doctests..."
Bio.GFF.GenericTools._test()
print "Bio.GFF.GenericTools doctests complete."

print "Running Bio.GFF.easy doctests..."
Bio.GFF.easy._test()
print "Bio.GFF.easy doctests complete."

