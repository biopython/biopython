# Copyright 2001 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

#!/usr/bin/env python
"""Test the IntelliGenetics parser and make sure everything is working smoothly.
"""
# standard library
import os
import copy

# Biopython
from Bio import File

from Bio import IntelliGenetics

intelligenetics_file_dir = os.path.join(os.getcwd(), 'IntelliGenetics')

test_files = [ 'TAT_mase_nuc.txt', 'VIF_mase-pro.txt', 'vpu_nucaligned.txt' ]



# test the parser
record_parser = IntelliGenetics.RecordParser()

for test in test_files:
    datafile = os.path.join( 'IntelliGenetics', test )
    src_handle = open( datafile )
    iterator = IntelliGenetics.Iterator(src_handle, record_parser)
    while 1:
        data = iterator.next()
        if not data: break
        print data

