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
from Bio.Rebase import Rebase


tests = [ 'bamii.htm', 'cac81.htm', 'crei.htm', 'foki.htm', 'pvuii.htm', \
          'taqi.htm' ]
record_parser = Rebase.RebaseParser()


for test in tests:
    print "testing %s" % test
    datafile = os.path.join( 'Rebase', test )
    src_handle = open( datafile )
    data = record_parser.parse( src_handle )
    print data
    print '\n'