# Copyright 1999 by Cayte Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
from types import *
from Bio import File
from Bio import ParserSupport
from Bio.Prosite import Prodoc


#tests = [ 'ps00107.txt','ps00159.txt', 'ps00165.txt', 'ps00213.txt', 'ps00432.txt', \
#    'ps00488.txt', 'ps00546.txt' ]

tests = [ 'pdoc00100.txt', 'pdoc00113.txt', 'pdoc00144.txt', 'pdoc00149.txt', \
    'pdoc00340.txt', 'pdoc00424.txt', 'pdoc00472.txt', 'pdoc00640.txt', \
     'pdoc00787.txt', 'pdoc00933.txt' ]
record_parser = Prodoc.RecordParser()

def print_list( list ):
    for item in list:
        print( '    ' + str( item ) )

def print_references( list ):
    for item in list:
        print item.number
        print item.authors
        print item.citation

for test in tests:
#    print "testing %s" % test
    datafile = os.path.join( 'Prosite', 'Doc', test )
    src_handle = open( datafile )
    data = record_parser.parse( src_handle )
    print data.accession
    print 'prosite_refs'
    print_list( data.prosite_refs )
    print data.text
    print 'references'
    print_references( data.references )
    src_handle.close()
