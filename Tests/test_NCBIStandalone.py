# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
from Bio import ParserSupport
from Bio.Blast import NCBIStandalone


### _Scanner

print "Running tests on _Scanner"

# BLAST v 2.0.10
test_2010 = ['bt001', 'bt002', 'bt003', 'bt004', 'bt005',
             'bt006', 'bt007', 'bt009', 'bt010', 'bt011',
             'bt012', 'bt013', 'bt014', 'bt015', 'bt016',
             'bt017', 'bt018', 'bt039', 'bt040'
             ]

# BLAST v 2.0.11
test_2011 = ['bt040', 'bt041', 'bt042', 'bt043', 'bt044',
             'bt045', 'bt046', 'bt047', 'bt048', 'bt049',
             'bt050', 'bt051', 'bt052', 'bt053', 'bt054',
             'bt055', 'bt056', 'bt057', 'bt058', 'bt059'
             ]

tests = test_2010 + test_2011

scanner = NCBIStandalone._Scanner()
for test in tests:
    print "*" * 50, "TESTING %s" % test
    datafile = os.path.join("Blast", test)
    tc = ParserSupport.TaggingConsumer()
    scanner.feed(open(datafile), tc)
