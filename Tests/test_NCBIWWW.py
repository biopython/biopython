# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
from Bio import ParserSupport
from Bio.Blast import NCBIWWW

# Biopython 1.40
# blast is deprecated now, and qblast only works for program=blastp or program=blastn
# Therefore, we only have two tests.
# all_tests = [
#    'bt019', 'bt020', 'bt021', 'bt022', 'bt023',
#    'bt024', 'bt025', 'bt026', 'bt027', 'bt028',
#    'bt029', 'bt030', 'bt031', 'bt032', 'bt033',
#    'bt034', 'bt035', 'bt036', 'bt037', 'bt038',
#    'bt061', 'bt064', 'bt066', 'bt069', 'bt070'
#    ]

all_tests = ['bt100', 'bt101']

# In order to keep the output file sizes reasonable, only generate
# a bunch of output for a few of the tests.
#detailed_tests = [
#    'bt019',    # 2.0.10 blastp
#    'bt023',    # 2.0.10 blastp master-slave
#    'bt026',    # 2.0.10 blastn
#    'bt028',    # 2.0.10 blastx
#    'bt030',    # 2.0.10 tblastn
#    'bt032',    # 2.0.10 tblastx
#    ]

detailed_tests = ['bt100', # blastn
                  'bt101'] # blastp
### _Scanner

print "Running tests on _Scanner"
    
scanner = NCBIWWW._Scanner()
for test in all_tests:
    print "*" * 50, "TESTING %s" % test
    datafile = os.path.join("Blast", test)
    scanner.feed(open(datafile), ParserSupport.AbstractConsumer())

for test in detailed_tests:
    print "*" * 50, "TESTING %s" % test
    datafile = os.path.join("Blast", test)
    scanner.feed(open(datafile), ParserSupport.TaggingConsumer())

### BlastParser

print "Running tests on BlastParser"

parser = NCBIWWW.BlastParser()
for test in all_tests:
    print "*" * 50, "TESTING %s" % test
    datafile = os.path.join("Blast", test)
    rec = parser.parse(open(datafile))
    # XXX should check this more thoroughly
    
