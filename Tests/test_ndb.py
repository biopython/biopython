# Copyright 2002 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

#!/usr/bin/env python
"""Test the Ndb parser and make sure everything is working smoothly.
"""
# standard library
import os
import copy

# Biopython
from Bio import File

from Bio import Ndb


test_files = ['AH0014.htm', 'PD0072.htm', 'PD0080.htm', 'PD0082.htm', 'PR0002.htm',
              'PR0004.htm', 'PR0007.htm', 'TRNA05.htm', 'URX035.htm', 'ZHF026.htm' ]



# test the parser
record_parser = Ndb.NdbParser()

for test in test_files:
    datafile = os.path.join( 'Ndb', test )
    src_handle = open( datafile )
    data = record_parser.parse( src_handle )
    print data
    print '\n'

