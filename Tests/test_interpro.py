# Copyright 2001 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

#!/usr/bin/env python
"""Test the InterPro parser and make sure everything is working smoothly.
"""
# standard library
import os
import copy

# Biopython
from Bio import File

from Bio import InterPro

test_files = ['IPR001064.htm', 'IPR001171.htm', 'IPR001391.htm',
              'IPR001442.htm', 'IPR001571.htm']



# test the parser
record_parser = InterPro.InterProParser()

for test in test_files:
    datafile = os.path.join( 'InterPro', test )
    src_handle = open( datafile )
    data = record_parser.parse( src_handle )
    print data
    print '\n'

