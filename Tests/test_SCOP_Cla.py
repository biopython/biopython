# Copyright 2001 by Gavin E. Crooks.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


"""Unit test for Cla"""

import unittest

from Bio.SCOP import Cla
from Bio.SCOP.Residues import Residues




class ClaTests(unittest.TestCase):

    def setUp(self):
        self.filename = './SCOP/dir.cla.scop.txt_test'

    def testParse(self):
        """Can we parse a CLA file?"""
        f=open(self.filename)
        try: 
            count = 0
            records = Cla.parse(f)
            for record in records:
                count +=1
            assert count == 14, "Wrong number of records?!"
        finally:
            f.close()
    
    def testStr(self):
        f = open(self.filename)
        try: 
            for line in f:
                record = Cla.Record(line)
                #End of line is platform dependent. Strip it off
                assert str(record).rstrip() == line.rstrip()
        finally:
            f.close()        

    def testError(self):
        corruptRec = "49268\tsp\tb.1.2.1\t-\n"

        try:
            record = Cla.Record(corruptRec)
            assert False, "Should never get here"
        except ValueError, e:
            pass

    def testRecord(self):
        recLine = 'd1dan.1\t1dan\tT:,U:91-106\tb.1.2.1\t21953\tcl=48724,cf=48725,sf=49265,fa=49266,dm=49267,sp=49268,px=21953'

        record = Cla.Record(recLine)
        assert record.sid =='d1dan.1'
        assert record.residues.pdbid =='1dan'
        assert record.residues.fragments ==(('T','',''),('U','91','106'))
        assert record.sccs == 'b.1.2.1'
        assert record.sunid == 21953
        assert record.hierarchy == [['cl',48724],['cf',48725],['sf',49265],
             ['fa',49266],['dm',49267],['sp',49268],
             ['px',21953]], record.hierarchy

    def testIndex(self):
        index = Cla.Index(self.filename)
        
        assert len(index)==14
        assert index.has_key('d4hbia_')

        rec = index['d1hbia_']
        assert rec.sunid == 14996



if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
