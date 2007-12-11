# Copyright 2001 by Gavin E. Crooks.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


"""Unit test for Dom

This test requires the mini DOM file 'testDom.txt'
"""

import unittest

from Bio.SCOP import Dom
from Bio.SCOP.Residues import Residues

import sys

def run_tests(argv):
    test_suite = testing_suite()
    runner = unittest.TextTestRunner(sys.stdout, verbosity = 2)
    runner.run(test_suite)

def testing_suite():
    return test_suite()


class DomTests(unittest.TestCase):
    def setUp(self) :
        self.filename = './SCOP/testDom.txt'

    def testParse(self):
       f = open(self.filename)
       try: 
           count = 0
           i = Dom.Iterator(f, Dom.Parser())
           while 1 :
               rec = i.next() 
               if rec is None : break
               count +=1
           self.assertEquals(count,10)
       finally:
           f.close()
    
    def testStr(self):
       f = open(self.filename)
       try: 
           p = Dom.Parser()
           i = Dom.Iterator(f)
           while 1 :
               line = i.next() 
               if line is None : break
               rec = p.parse(line)
               #End of line is platform dependant. Strip it off
               self.assertEquals(str(rec).rstrip(),line.rstrip())
       finally:
           f.close()        

    def testError(self) :
        corruptDom = "49xxx268\tsp\tb.1.2.1\t-\n"
        self.assertRaises(ValueError, Dom.Parser().parse, corruptDom)


    def testRecord(self) :
        recLine = 'd7hbib_\t7hbi\tb:\t1.001.001.001.001.001'

        rec = Dom.Parser().parse(recLine)
        self.assertEquals(rec.sid, 'd7hbib_')
        self.assertEquals(rec.residues.pdbid,'7hbi')
        self.assertEquals(rec.residues.fragments,(('b','',''),) )        
        self.assertEquals(rec.hierarchy,'1.001.001.001.001.001')


def test_suite():
    return unittest.makeSuite(DomTests)

if __name__ == '__main__':
    unittest.main()







