# Copyright 2001 by Gavin E. Crooks.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


"""Unit test for Hie"""

import unittest

from Bio.SCOP import Hie
from Bio.SCOP.Residues import Residues




class HieTests(unittest.TestCase):

    def setUp(self):
        self.filename = './SCOP/dir.hie.scop.txt_test'

    def testParse(self):
       f = open(self.filename)
       try: 
           count = 0
           for record in Hie.parse(f):
               count +=1
           assert count == 21, "Wrong number of records?! "+str(count)
       finally:
           f.close()
    
    def testStr(self):
       f = open(self.filename)
       try: 
           for line in f:
               record = Hie.Record(line)
               #End of line is platform dependent. Strip it off
               assert str(record).rstrip() == line.rstrip()
       finally:
           f.close()        

    def testError(self):
        corruptRec = "4926sdfhjhfgyjdfyg"

        try:
            rec = Hie.Record(corruptRec)
            assert False, "Should never get here"
        except ValueError, e:
            pass


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
