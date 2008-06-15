"""Tests the basic functionality of the CDD parsers.
"""

import os
import Bio.CDD

testfiles   = [  'pfam00035', 'pfam01356', 'pfam02903', 'smart00499', 'smart00505' ]
for file in testfiles:
    fh = open(os.path.join("CDD", file + '.htm' ))
    print "Testing Bio.CDD on " + file + "\n\n"
    records = Bio.CDD.Iterator( fh, Bio.CDD.RecordParser(debug_level=0))
#    records = Bio.CDD.Iterator( fh, None)
    while 1:
        record = records.next()
        if record is not None:
            print record
        else:
            break
    print "\n"





