# Copyright 2001 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

#!/usr/bin/env python
"""Test the GenBank parser and make sure everything is working smoothly.
"""
# standard library
import os
import copy

# Biopython
from Bio import File

from Bio import Kabat

kabat_file_dir = os.path.join(os.getcwd(), 'Kabat')

test_files = [ 'k000245.txt', 'k000255.txt', 'k000347.txt', 'k000397.txt', \
              'k000456.txt']



# test the parser
record_parser = Kabat.RecordParser()

for test in test_files:
    datafile = os.path.join( 'Kabat', test )
    src_handle = open( datafile )
    iterator = Kabat.Iterator(src_handle, record_parser)
    data = iterator.next()
    data.print_kabat()
    print '\n'

