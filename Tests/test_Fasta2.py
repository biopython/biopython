# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
from types import *
from Bio import File
from Bio import ParserSupport
from Bio.Fasta import Fasta


tests = [ 'lupine.nu', 'elderberry.nu', 'phlox.nu', 'centaurea.nu', \
    'wisteria.nu', 'sweetpea.nu', 'lavender.nu' ]
record_parser = Fasta.RecordParser()
sequence_parser = Fasta.SequenceParser()

for test in tests:
    print "testing %s" % test
    datafile = os.path.join( 'Nucleic', test )
    src_handle = open( datafile )
    data = record_parser.parse( src_handle )
    print data


for test in tests:
    print "testing %s" % test
    datafile = os.path.join( 'Nucleic', test )
    src_handle = open( datafile )
    data = sequence_parser.parse( src_handle )
    print data.name
    print data.seq

tests = [ 'aster.pro', 'rosemary.pro', 'rose.pro', 'loveliesbleeding.pro' ]
for test in tests:
    print "testing %s" % test
    datafile = os.path.join( 'Amino', test )
    src_handle = open( datafile )
    data = record_parser.parse( src_handle )
    print data

for test in tests:
    print "testing %s" % test
    datafile = os.path.join( 'Amino', test )
    src_handle = open( datafile )
    data = sequence_parser.parse( src_handle )
    print data.name
    print data.seq[ :60]
