# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
from Bio import ParserSupport
from Bio.SwissProt import KeyWList


### _Scanner

print "Running tests on _Scanner"

tests = ['kw001', 'kw002']

class MyConsumer(ParserSupport.TaggingConsumer):
    def __init__(self, *args, **keywds):
        ParserSupport.TaggingConsumer.__init__( *(self,) + args, **keywds)
        self._keywd = 0

    # Only print the first keyword, so I don't generate a bunch of
    # output.
    def keyword(self, line):
        if not self._keywd:
            self._print_name("keyword", line)
            self._keywd = 1
        

scanner = KeyWList._Scanner()
for test in tests:
    print "testing %s" % test
    datafile = os.path.join("SwissProt", test)
    tc = MyConsumer()
    scanner.feed(open(datafile), tc)

