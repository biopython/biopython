# Copyright 2004 by Jason A. Hackney.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
import string

#Bio.MEME is deprecated now, but I want to keep testing it until
#it is removed (without showing the deprecation warning)
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
from Bio.MEME import Parser
warnings.resetwarnings()

from Bio import ParserSupport
from Bio.File import UndoHandle

meme_tests = [
    "meme.dna.oops.txt",
    "meme.protein.oops.txt",
    "meme.protein.tcm.txt"
    ]
    
mast_tests = [
    "mast.dna.oops.txt",
    "mast.protein.oops.txt",
    "mast.protein.tcm.txt"
    ]


print "Testing MEME Scanner"

datafile = os.path.join("MEME", meme_tests[0])
uhandle = UndoHandle(open(datafile))
scanner = Parser._MEMEScanner()
consumer = ParserSupport.TaggingConsumer()
scanner.feed(uhandle, consumer)
    
print "Running tests on MEME parser"    

meme_parser = Parser.MEMEParser()

for test in meme_tests:
    print "*" * 50, "TESTING %s" % test
    datafile = os.path.join("MEME", test)
    rec = meme_parser.parse(open(datafile))

print "Testing MEME Scanner"

datafile = os.path.join("MEME", mast_tests[0])
uhandle = UndoHandle(open(datafile))
scanner = Parser._MASTScanner()
consumer = ParserSupport.TaggingConsumer()
scanner.feed(uhandle, consumer)    
   
print "Running tests on MAST parser"

mast_parser = Parser.MASTParser()

for test in mast_tests:
    print "*" * 50, "TESTING %s" % test
    datafile = os.path.join("MEME", test)
    rec = mast_parser.parse(open(datafile))
