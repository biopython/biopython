"""Tests the basic functionality of the CDD parsers.
"""

import os
import sys

import Bio.CDD
from Bio.FilteredReader import FilteredReader
from Bio.FilteredReader import remove_empty_line
from Bio.FilteredReader import remove_leading_whitespace
from Bio.SGMLExtractor import SGMLExtractorHandle

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





