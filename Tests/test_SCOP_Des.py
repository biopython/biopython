# Copyright 2001 by Gavin E. Crooks.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


"""Unit test for Des"""

import unittest

from Bio.SCOP import Des
from Bio.SCOP.Residues import Residues



class DesTests(unittest.TestCase):

    def setUp(self):
        self.filename = './SCOP/dir.des.scop.txt_test'

    def testParse(self):
       f = open(self.filename)
       try: 
           count = 0
           records = Des.parse(f)
           for record in records:
               count +=1
           assert count == 20, "Wrong number of records?!"
       finally:
           f.close()
    
    def testStr(self):
       f = open(self.filename)
       try: 
           for line in f:
               record = Des.Record(line)
               #End of line is plateform dependant. Strip it off
               assert str(record).rstrip() == line.rstrip()
       finally:
           f.close()        

    def testError(self):
        corruptRec = "49268\tsp\tb.1.2.1\t-\n"
        try:
            record = Des.Record(corruptRec)
            assert False, "Should never get here"
        except ValueError, e:
            pass

    def testRecord(self):
        recLine = '49268\tsp\tb.1.2.1\t-\tHuman (Homo sapiens)    \n'
        recFields = (49268,'sp','b.1.2.1','','Human (Homo sapiens)')

        record = Des.Record(recLine)
        assert record.sunid == recFields[0]
        assert record.nodetype == recFields[1]
        assert record.sccs == recFields[2]
        assert record.name == recFields[3]                
        assert record.description == recFields[4]



if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
