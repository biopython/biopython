# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
from Bio import ParserSupport
from Bio.Blast import NCBIWWW


test_files = ['bt019', 'bt020', 'bt021', 'bt022', 'bt023',
              'bt024', 'bt025', 'bt026', 'bt027', 'bt028',
              'bt029', 'bt030', 'bt031', 'bt032', 'bt033',
              'bt034', 'bt035', 'bt036', 'bt037', 'bt038'
              ]



### _Scanner

print "Running tests on _Scanner"
    
scanner = NCBIWWW._Scanner()
for test in test_files:
    print "*" * 50, "TESTING %s" % test
    datafile = os.path.join("Blast", test)
    tc = ParserSupport.TaggingConsumer()
    scanner.feed(open(datafile), tc)


### BlastParser

print "Running tests on BlastParser"

for test in test_files:
    print "*" * 50, "TESTING %s" % test
    datafile = os.path.join("Blast", test)
    parser = NCBIWWW.BlastParser()
    rec = parser.parse(open(datafile))
    print "parsed without exception"
    # XXX should check this more thoroughly
    
