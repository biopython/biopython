# Copyright 1999 by Jeffrey Chang, Andrew Dalke.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
from types import *
from TestSupport import verbose, TestFailed
from Bio import ParserSupport
from Bio import Sequence


### Sequence

if verbose:
    print "Running tests on Sequence"

s = Sequence.Sequence()
try:
    assert type(s.seq) is StringType, "seq should be string"
    s.seq = "ANDREW"
    assert s[1:3] == "ND", "slicing error"
except Exception, x:
    raise TestFailed, "Sequence (%s)" % x
    

### SubSequence

if verbose:
    print "Running tests on SubSequence"

s = Sequence.Sequence("ANDREW")
try:
    assert s.subseq(1,3) == "AND", "__call__"
    assert s.subseq.omg(1, 3) == "AN", "omg"
    assert s.subseq.perl(1, 3) == "NDR", "perl"
    assert s.subseq.python(1, 3) == "ND", "python"
except Exception, x:
    raise TestFailed, "SubSequence (%s)" % x


### NamedSequence

if verbose:
    print "Running tests on NamedSequence"

s = Sequence.Sequence("ANDREW")
ns = Sequence.NamedSequence(s)
try:
    assert ns.name == '', "name"
    assert ns.uid == '', "uid"
    assert ns.dbid == '', "dbid"
    assert ns[1:3] == "ND", "delegation"
except Exception, x:
    raise TestFailed, "NamedSequence (%s)" % x
