"""Tests the basic functionality of the SAF parsers.
"""

import os

import Bio.Saf
from Bio.File import SGMLHandle

testfiles   = [  'saf1.txt' ]

for file in testfiles:
    fh = open(os.path.join("Saf", file))
    print "Testing Bio.Saf on " + file + "\n\n"
    records = Bio.Saf.Iterator( fh, Bio.Saf.RecordParser(debug_level=0))
    while 1:
#    for j in range( 0, 1 ):
        record = records.next()
        if record is not None:
            print record
            pass
        else:
            break
    print "\n"





