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

from Bio.LocusLink.web_parse import LocusLinkParser


test_files = ['Locus14789.htm', 'Locus16590.htm', 'Locus2523.htm', 'Locus5797.htm', 'Locus6373.htm' ]

 



# test the parser
record_parser = LocusLinkParser()

for test in test_files:
    datafile = os.path.join( 'LocusLink', test )
    print test
    src_handle = open( datafile )
    data = record_parser.parse( src_handle )
    print data
    print '\n'

