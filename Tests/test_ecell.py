# Copyright 2001 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

#!/usr/bin/env python
"""Test the ECell parser and make sure everything is working smoothly.
"""
# standard library
import os
import copy

# Biopython
from Bio import File
from Bio import FilteredReader

from Bio import ECell

test_files = [ 'sample.txt' ]



# test the parser
record_parser = ECell.RecordParser()

for test in test_files:
    datafile = os.path.join( 'ECell', test )
    src_handle = open( datafile )
    filtered_reader = FilteredReader.FilteredReader( src_handle )
    iterator = ECell.Iterator( filtered_reader, record_parser)
#    iterator = ECell.Iterator( filtered_reader, None)
    while 1:
        data = iterator.next()
        if not data: break
        print data

