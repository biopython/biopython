# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
import string
from Bio import ParserSupport
from Bio.Blast import NCBIStandalone


all_tests = [
    'bt001', 'bt002', 'bt003', 'bt004', 'bt005',
    'bt006', 'bt007', 'bt009', 'bt010', 'bt011',
    'bt012', 'bt013', 'bt014', 'bt015', 'bt016',
    'bt017', 'bt018', 'bt039', 'bt040', 'bt041',
    'bt042', 'bt043', 'bt044', 'bt045', 'bt046',
    'bt047', 'bt048', 'bt049', 'bt050', 'bt051',
    'bt052', 'bt053', 'bt054', 'bt055', 'bt056',
    'bt057', 'bt058', 'bt059', 'bt060', 'bt062',
    'bt063', 'bt067', 'bt071'
    ]
# In order to keep the output file sizes reasonable, only generate
# a bunch of output for a few of the tests.
detailed_tests = [
    'bt001',   # 2.0.10 blastp
    'bt003',   # 2.0.10 blastp master-slave
    'bt006',   # 2.0.10 blastpgp
    'bt010',   # 2.0.10 blastn
    'bt014',   # 2.0.10 blastx
    'bt016',   # 2.0.10 tblastn
    'bt018',   # 2.0.10 tblastx
    'bt042',   # 2.0.11 blastp
    ]

### _Scanner

print "Running tests on _Scanner"

scanner = NCBIStandalone._Scanner()
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

parser = NCBIStandalone.BlastParser()
pb_parser = NCBIStandalone.PSIBlastParser()
for test in all_tests:
    print "*" * 50, "TESTING %s" % test
    datafile = os.path.join("Blast", test)
    try:
        # First, try parsing it with the normal parser.
        rec = parser.parse(open(datafile))
    except ValueError, x:
        # If it complains that the input is psiblast data, then
        # parse it with the psiblast parser.
        if string.find(str(x), 'PSI-BLAST data') >= 0:
            rec = pb_parser.parse(open(datafile))
        else:
            raise

### Blast Record
# - incredibly incomplete. Just testing what I'm adding -- Brad.

print "Running tests on Blast Records"

datafile = os.path.join("Blast", all_tests[4]) # bt005
rec = parser.parse(open(datafile))

print "\tTesting conversion of multiple alignments"
from Bio.Alphabet import IUPAC
generic_align = rec.multiple_alignment.to_generic(IUPAC.protein)
test_seq = generic_align.get_seq_by_num(0)
assert test_seq.alphabet == IUPAC.protein
assert test_seq.data[:60] == rec.multiple_alignment.alignment[0][2]
