# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
from types import *
from TestSupport import verbose, TestFailed
from Bio import File
from Bio import ParserSupport
from Bio.Fasta import Fasta


### Record

if verbose:
    print "Running tests on Record"

r = Fasta.Record()
try:
    assert type(r.title) is StringType, "title should be string"
    assert type(r.sequence) is StringType, "sequence should be string"
except Exception, x:
    raise TestFailed, "Record (%s)" % x
    

### Scanner

if verbose:
    print "Running tests on Scanner"

tests = ['f001', 'f002']

class TestHandle:
    def __init__(self, h):
        self._h = h
    def write(self, s):
        assert self._h.readline() == s

scanner = Fasta.Scanner()
for test in tests:
    datafile = os.path.join("Fasta", test)
    modelfile = datafile + ".tagged"
    tc = ParserSupport.TaggingConsumer(handle=TestHandle(open(modelfile)))
    try:
        scanner.feed(open(datafile), tc)
    except:
        raise TestFailed, "Scanner (%s)" % test


### StandardConsumer

if verbose:
    print "Running tests on StandardConsumer"

c = Fasta.StandardConsumer()
try:
    c.start_sequence()
    c.title('>This is a title\n')
    c.sequence('ABCD\n')
    c.sequence('EFG\n')
    c.end_sequence()

    assert c.data.title == "This is a title", "title is incorrect"
    assert c.data.sequence == "ABCDEFG", "sequence is incorrect"

    c.start_sequence()
    c.end_sequence()
    assert c.data.title == '', "record should be cleared"
except Exception, x:
    raise TestFailed, "StandardConsumer (%s)" % x


### SequenceConsumer

if verbose:
    print "Running tests on SequenceConsumer"

c = Fasta.SequenceConsumer()
try:
    c.start_sequence()
    c.title('>This is a title\n')
    c.sequence('ABCD\n')
    c.sequence('EFG\n')
    c.end_sequence()

    assert c.data.name == "This is a title", "title is incorrect"
    assert c.data.seq == "ABCDEFG", "sequence is incorrect"

    c.start_sequence()
    c.end_sequence()
    assert c.data.name == '', "record should be cleared"
except Exception, x:
    raise TestFailed, "SequenceConsumer (%s)" % x


### Iterator

if verbose:
    print "Running tests on Iterator"

i = Fasta.Iterator(open(os.path.join('Fasta', 'f002')))
try:
    assert i.next()[:11] == '>gi|1348912', "sequence 1"
    assert i.next()[:11] == '>gi|1348917', "sequence 2"
    assert i.next()[:11] == '>gi|1592936', "sequence 3"
    assert not i.next(), "no more sequences"
except Exception, x:
    raise TestFailed, "Iterator (%s)" % x
