# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
from types import *
from Bio import File
from Bio import ParserSupport
from Bio.Rebase import Rebase


tests = [ 'cac81.htm', 'bamii.htm', 'pvuii.htm', 'taqi.htm' ]
#    'ps00488.txt', 'ps00546.txt' ]
record_parser = Rebase.RecordParser()

def print_list( list ):
    for item in list:
        print( '    ' + str( item ) )

for test in tests:
    print "testing %s" % test
    datafile = os.path.join( 'Rebase', test )
    src_handle = open( datafile )
    data = record_parser.parse( src_handle )
    print data.sequence
    print '!!! %s ' % data.enzyme_num
    print data.methylation
    print data.prototype
    print data.source
    print data.microorganism
    print data.temperature
    print data.date_entered
    print data.date_modified
    print data.num_Adeno2
    print data.num_Lambda
    print data.num_pBR322
    print data.num_PhiX174
    print data.num_SV40
