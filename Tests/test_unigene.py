# Copyright 1999 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
import types
import string
from types import *
from Bio import File
from Bio import ParserSupport
from Bio import UniGene


tests = [ 'Bt145.htm', 'Dr20.htm', 'Hs28199.htm', 'Rn35.htm', 'Mm421.htm', \
          'Rn35.htm', 'Rn20.htm', 'Hs13225.htm', 'Mm28919.htm', 'Hs227583.htm' ]
record_parser = UniGene.UniGeneParser()


for test in tests:
    print "testing %s" % test
    datafile = os.path.join( 'UniGene', test )
    src_handle = open( datafile )
    data = record_parser.parse( src_handle )
    print data
#    print '\n'