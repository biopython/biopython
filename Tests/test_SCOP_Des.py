# Copyright 2001 by Gavin E. Crooks.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


"""Unit test for Des"""

import unittest

from Bio.SCOP import Des
from Bio.SCOP.Residues import Residues

import sys

def run_tests(argv):
    test_suite = testing_suite()
    runner = unittest.TextTestRunner(sys.stdout, verbosity = 2)
    runner.run(test_suite)

def testing_suite():
    return test_suite()


class DesTests(unittest.TestCase):

    def setUp(self) :
        self.filename = './SCOP/dir.des.scop.txt_test'

    def testParse(self):
       f = open(self.filename)
       try: 
           count = 0
           i = Des.Iterator(f, Des.Parser())
           while 1 :
               rec = i.next() 
               if rec is None : break
               count +=1
           assert count == 20, "Wrong number of records?!"
       finally:
           f.close()
    
    def testStr(self):
       f = open(self.filename)
       try: 
           p = Des.Parser()
           i = Des.Iterator(f)
           while 1 :
               line = i.next() 
               if line is None : break
               rec = p.parse(line)
               #End of line is plateform dependant. Strip it off
               assert str(rec).rstrip() == line.rstrip()
       finally:
           f.close()        

    def testError(self) :
        corruptRec = "49268\tsp\tb.1.2.1\t-\n"
        p = Des.Parser()

        try:
            rec = p.parse(corruptRec)
            assert 0, "Should never get here"
        except SyntaxError, e :
            pass

    def testRecord(self) :
        recLine = '49268\tsp\tb.1.2.1\t-\tHuman (Homo sapiens)    \n'
        recFields = (49268,'sp','b.1.2.1','','Human (Homo sapiens)')

        rec = Des.Parser().parse(recLine)
        assert rec.sunid == recFields[0]
        assert rec.nodetype == recFields[1]
        assert rec.sccs == recFields[2]
        assert rec.name == recFields[3]                
        assert rec.description == recFields[4]


def test_suite():
    return unittest.makeSuite(DesTests)

if __name__ == '__main__':
    unittest.main()




