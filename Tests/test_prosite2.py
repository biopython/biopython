# Copyright 1999 by Cayte Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
from types import *
from Bio import File
from Bio import ParserSupport
from Bio import Prosite


#tests = [ 'ps00107.txt', 'ps00159.txt', 'ps00165.txt', 'ps00213.txt', 'ps00432.txt', \
#    'ps00488.txt', 'ps00546.txt' ]

tests = [ 'ps00107.txt', 'ps00159.txt', 'ps00165.txt', 'ps00432.txt', \
    'ps00488.txt', 'ps00546.txt' ]
record_parser = Prosite.RecordParser()

def print_list( list ):
    for item in list:
        print( '    ' + str( item ) )

for test in tests:
    print "testing %s" % test
    datafile = os.path.join( 'Prosite', test )
    src_handle = open( datafile )
    data = record_parser.parse( src_handle )
    print data.name
    print data.type
    print data.accession
    print data.created
    print data.data_update
    print data.info_update
#    print data.data_info
    print data.pdoc
    print data.description
    print data.pattern
    print data.matrix
    print data.rules
    print data.nr_sp_release
    print data.nr_sp_seqs
    print data.cc_taxo_range
    print data.cc_max_repeat
#    print data.name
    print 'cc_site'
    print_list( data.cc_site )
    print 'dr_positive'
    print_list( data.dr_positive )
    print 'dr_false_neg'
    print_list( data.dr_false_neg )
    print 'dr_false_pos'
    print_list( data.dr_false_pos )
    print 'dr_potential'
    print_list( data.dr_potential )
    print 'dr_unknown'
    print_list( data.dr_unknown )
    print 'pdb_structs'
    print_list( data.pdb_structs )

