"""Tests the basic functionality of the GEO parsers.
"""

import os, sys

import Bio.Geo

testfiles   = [  'GSE16.txt', 'GSM645.txt', 'GSM691.txt', 'GSM700.txt', 'GSM804.txt'  ]

#Five additional files from the NCBI to document the GEO SOFT file format
#changes made in 2005.  Note the new table_begin and table_end lines.
testfiles.extend(['soft_ex_affy.txt',
                  'soft_ex_affy_chp.txt',
                  'soft_ex_dual.txt',
                  'soft_ex_family.txt',
                  'soft_ex_platform.txt',
                  'soft_ex_series.txt'])

for file in testfiles:
    if sys.version_info[0] >= 3:
        #Python 3 problem: Can't use utf8 on Tests/Geo/soft_ex_*.txt
        #due to micro (\xb5) and degrees (\xb0) symbols
        fh = open(os.path.join("Geo", file), encoding="latin")
    else:
        fh = open(os.path.join("Geo", file))
    print "Testing Bio.Geo on " + file + "\n\n"
    records = Bio.Geo.parse(fh)
    for record in records:
        print record
    print "\n"





