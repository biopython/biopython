"""Tests the basic functionality of the GEO parsers.
"""

import os

import Bio.Geo
from Bio.File import SGMLHandle

testfiles   = [  'GSE16.txt', 'GSM645.txt', 'GSM691.txt', 'GSM700.txt', 'GSM804.txt'  ]

for file in testfiles:
    fh = open(os.path.join("Geo", file))
    print "Testing Bio.Geo on " + file + "\n\n"
    records = Bio.Geo.Iterator( fh, Bio.Geo.RecordParser(debug_level=0))
    while 1:
#    for j in range( 0, 1 ):
        record = records.next()
        if record is not None:
            print record
            pass
        else:
            break
    print "\n"





