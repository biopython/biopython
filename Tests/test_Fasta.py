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
        scanner.feed(File.open(datafile), tc)
    except:
        raise TestFailed, "Scanner (%s)" % test


### Consumer

if verbose:
    print "Running tests on Consumer"

c = Fasta.Consumer()
try:
    c.start_sequence()
    c.title('>This is a title\n')
    c.sequence('ABCD\n')
    c.sequence('EFG\n')
    c.end_sequence()

    assert len(c.records) == 1, "consumer should have 1 record"
    assert c.records[0].title == "This is a title", "title is incorrect"
    assert c.records[0].sequence == "ABCDEFG", "sequence is incorrect"

    c.start_sequence()
    c.end_sequence()
    assert len(c.records) == 2, "consumer should have 2 records"
except Exception, x:
    raise TestFailed, "Consumer (%s)" % x


### Iterator

if verbose:
    print "Running tests on Iterator"

i = Fasta.Iterator(File.open(os.path.join('Fasta', 'f002')))
try:
    assert i.next().title[:10] == 'gi|1348912', "sequence 1"
    assert i.next().title[:10] == 'gi|1348917', "sequence 2"
    assert i.next().title[:10] == 'gi|1592936', "sequence 3"
    assert i.next() is None, "no more sequences"
except Exception, x:
    raise TestFailed, "Iterator (%s)" % x


### parse

if verbose:
    print "Running tests on parse"

try:
    datafile = os.path.join("Fasta", 'f002')
    f = Fasta.parse(File.open(datafile))

    assert len(f) == 3, "Should have 3 records"
    assert str(f[2]) == '>gi|1592936|gb|G29385|G29385 human STS SHGC-32652\012GATCAAATCTGCACTGTGTCTACATATAGGAAAGGTCCTGGTGTGTGCTAATGTTCCCAATGCAGGACTTGAGGAAGAGC\012TCTGTTATATGTTTCCATTTCTCTTTATCAAAGATAACCAAACCTTATGGCCCTTATAACAATGGAGGCACTGGCTGCCT\012CTTAATTTTCAATCATGGACCTAAAGAAGTACTCTGAAGGGTCTCAACAATGCCAGGTGGGGACAGATATACTCAGAGAT\012TATCCAGGTCTGCCTCCCAGCGAGCCTGGAGTACACCAGACCCTCCTAGAGAAATCTGTTATAATTTACCACCCACTTAT\012CCACCTTTAAACTTGGGGAAGGNNGCNTTTCAAATTAAATTTAATCNTNGGGGGNTTTTAAACTTTAACCCTTTTNCCNT\012TNTNGGGGTNGGNANTTGNCCCCNTTAAAGGGGGNNCCCCTNCNNGGGGGAATAAAACAANTTNNTTTTTT', "__str__ isn't working correctly"
    
except Exception, x:
    raise TestFailed, "parse (%s)" % x
