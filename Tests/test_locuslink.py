"""Test the basic functionality of the LocusLink parsers.
"""

import os
import sys

import Bio.LocusLink
testfiles   = [  'LL_block', ]
for file in testfiles:
    fh = open(os.path.join("LocusLink", file ))
    print "Testing Bio.LocusLink on " + file + "\n\n"
    records = Bio.LocusLink.Iterator( fh, Bio.LocusLink.RecordParser(debug_level=0))
#    records = Bio.LocusLink.Iterator( fh, None)
#    while 1:
    for j in range( 0, 10 ):
        print 'record %d\n' % j
        record = records.next()
        if record is not None:
#            pass
            print record
        else:
            break
    print "\n"





