# Copyright 2001 by Gavin E. Crooks.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


"""Unit test for Hie"""

import unittest

from Bio.SCOP import Hie
from Bio.SCOP.Residues import Residues
from Bio.SCOP.tests import findResource


class HieTests(unittest.TestCase):

    def setUp(self) :
        self.filename = findResource('Bio/SCOP/tests/hietest.txt')

    def testParse(self):
       f = open(self.filename)
       try: 
           count = 0
           i = Hie.Iterator(f, Hie.Parser())
           while 1 :
               rec = i.next() 
               if rec is None : break
               count +=1
           assert count == 21, "Wrong number of records?! "+str(count)
       finally:
           f.close()
    
    def testStr(self):
       f = open(self.filename)
       try: 
           p = Hie.Parser()
           i = Hie.Iterator(f)
           while 1 :
               line = i.next() 
               if line is None : break
               rec = p.parse(line)
               #End of line is plateform dependant. Strip it off
               assert str(rec).rstrip() == line.rstrip()
       finally:
           f.close()        

    def testError(self) :
        corruptRec = "4926sdfhjhfgyjdfyg"
        p = Hie.Parser()

        try:
            rec = p.parse(corruptRec)
            assert 0, "Should never get here"
        except SyntaxError, e :
            pass

def test_suite():
    return unittest.makeSuite(HieTests)

if __name__ == '__main__':
    unittest.main()








