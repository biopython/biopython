"""Tests the basic functionality of the GEO parsers.
"""

import os

import Bio.Geo
from Bio.File import SGMLHandle

testfiles   = [  'GSE16.txt', 'GSM645.txt', 'GSM691.txt', 'GSM700.txt', 'GSM804.txt'  ]

#Five additional files from the NCBI to document the GEO SOFT file format
#changes made in 2005.  Note the new table_begin and table_end lines.
testfiles.extend(['soft_ex_affy.txt',
                  'soft_ex_dual.txt',
                  'soft_ex_family.txt',
                  'soft_ex_platform.txt',
                  'soft_ex_series.txt'])

for file in testfiles:
    fh = open(os.path.join("Geo", file))
    print "Testing Bio.Geo on " + file + "\n\n"
    records = Bio.Geo.Iterator( fh, Bio.Geo.RecordParser(debug_level=0))
    while 1:
        record = records.next()
        if record is not None:
            print record
            pass
        else:
            break
    print "\n"





