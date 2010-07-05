# Copyright 2001 by Gavin E. Crooks.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


"""Unit test for Dom

This test requires the mini DOM file 'testDom.txt'
"""

import unittest

from Bio.SCOP import Dom




class DomTests(unittest.TestCase):
    def setUp(self):
        self.filename = './SCOP/testDom.txt'

    def testParse(self):
        """Test if all records in a DOM file are being read"""
        f = open(self.filename)
        try: 
            count = 0
            for record in Dom.parse(f):
                count +=1
            self.assertEqual(count,10)
        finally:
            f.close()
    
    def testStr(self):
        """Test if we can convert each record to a string correctly"""
        f = open(self.filename)
        try: 
            for line in f:
                record = Dom.Record(line)
                #End of line is platform dependent. Strip it off
                self.assertEqual(str(record).rstrip(),line.rstrip())
        finally:
            f.close()

    def testError(self):
        """Test if a corrupt record raises the appropriate exception"""
        corruptDom = "49xxx268\tsp\tb.1.2.1\t-\n"
        self.assertRaises(ValueError, Dom.Record, corruptDom)


    def testRecord(self):
        """Test one record in detail"""
        recLine = 'd7hbib_\t7hbi\tb:\t1.001.001.001.001.001'

        rec = Dom.Record(recLine)
        self.assertEqual(rec.sid, 'd7hbib_')
        self.assertEqual(rec.residues.pdbid,'7hbi')
        self.assertEqual(rec.residues.fragments,(('b','',''),) )        
        self.assertEqual(rec.hierarchy,'1.001.001.001.001.001')



if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
