# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
from TestSupport import verbose, TestFailed
from Bio import ParserSupport
from Bio.SwissProt import KeyWList



### Scanner

if verbose:
    print "Running tests on Scanner"

tests = ['kw001', 'kw002']

class TestHandle:
    def __init__(self, h):
        self._h = h
    def write(self, s):
        assert self._h.readline() == s

scanner = KeyWList.Scanner()
for test in tests:
    datafile = os.path.join("SwissProt", test)
    modelfile = datafile + ".tagged"
    tc = ParserSupport.TaggingConsumer(handle=TestHandle(open(modelfile)))
    try:
        scanner.feed(open(datafile), tc)
    except:
        raise TestFailed, "Scanner (%s)" % test
