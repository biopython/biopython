# Copyright 2001 by Gavin E. Crooks.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


"""Unit test for Hie"""

import unittest

from Bio.SCOP import Hie




class HieTests(unittest.TestCase):

    def setUp(self):
        self.filename = './SCOP/dir.hie.scop.txt_test'

    def testParse(self):
        """Test if all records in a HIE file are being read"""
        f = open(self.filename)
        try: 
            count = 0
            for record in Hie.parse(f):
                count +=1
            self.assertEqual(count, 21)
        finally:
            f.close()

    def testStr(self):
        """Test if we can convert each record to a string correctly"""
        f = open(self.filename)
        try: 
            for line in f:
                record = Hie.Record(line)
                #End of line is platform dependent. Strip it off
                self.assertEqual(str(record).rstrip(), line.rstrip())
        finally:
            f.close()        

    def testError(self):
        """Test if a corrupt record raises the appropriate exception"""
        corruptRec = "4926sdfhjhfgyjdfyg"
        self.assertRaises(ValueError, Hie.Record, corruptRec)


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
