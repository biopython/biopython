# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
from TestSupport import verbose, TestFailed
from Bio import ParserSupport
from Bio.Blast import NCBIStandalone


### Scanner

if verbose:
    print "Running tests on Scanner"

tests = ['bt001', 'bt002', 'bt003', 'bt004', 'bt005',
         'bt006', 'bt007', 'bt009', 'bt010', 'bt011',
         'bt012', 'bt013', 'bt014', 'bt015', 'bt016',
         'bt017', 'bt018', 'bt039', 'bt040'
         ]

class TestHandle:
    def __init__(self, h):
        self._h = h
    def write(self, s):
        assert self._h.readline() == s

scanner = NCBIStandalone.Scanner()
for test in tests:
    datafile = os.path.join("Blast", test)
    modelfile = datafile + ".tagged"
    tc = ParserSupport.TaggingConsumer(handle=TestHandle(open(modelfile)))
    try:
        scanner.feed(open(datafile), tc)
    except:
        raise TestFailed, "Scanner (%s)" % test


