# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
from types import *
from Bio import File
from Bio import ParserSupport
from Bio import Fasta


### Record

print "Running test on Record"

r = Fasta.Record()
print type(r.title)    # StringType
print type(r.sequence) #  StringType
    

### _Scanner

print "Running tests on _Scanner"

tests = ['f001', 'f002']

scanner = Fasta._Scanner()
for test in tests:
    print "testing %s" % test
    datafile = os.path.join("Fasta", test)
    tc = ParserSupport.TaggingConsumer()
    scanner.feed(open(datafile), tc)


### _RecordConsumer

print "Running tests on _RecordConsumer"

c = Fasta._RecordConsumer()
c.start_sequence()
c.title('>This is a title\n')
c.sequence('ABCD\n')
c.sequence('EFG\n')
c.end_sequence()

print c.data.title        # "This is a title"
print c.data.sequence     # "ABCDEFG"
# clear the record
c.start_sequence()   
c.end_sequence()
print repr(c.data.title)  # ''



### _SequenceConsumer

print "Running tests on _SequenceConsumer"

c = Fasta._SequenceConsumer()
c.start_sequence()
c.title('>This is a title\n')
c.sequence('ABCD\n')
c.sequence('EFG\n')
c.end_sequence()

print c.data.description  # "This is a title"
print c.data.seq          # "Seq('ABCDEFG', Alphabet())"
# clear the record
c.start_sequence()   
c.end_sequence()
print repr(c.data.name)   # '<unknown name>'



### Iterator

print "Running tests on Iterator"

i = Fasta.Iterator(open(os.path.join('Fasta', 'f002')))
print i.next()   # '>gi|1348912' [...]
print i.next()   # '>gi|1348917'
print i.next()   # '>gi|1592936'
print repr(i.next())   # None
