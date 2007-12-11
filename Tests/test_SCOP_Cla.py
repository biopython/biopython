# Copyright 2001 by Gavin E. Crooks.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


"""Unit test for Cla"""

import unittest

from Bio.SCOP import Cla
from Bio.SCOP.Residues import Residues

import sys

def run_tests(argv):
    test_suite = testing_suite()
    runner = unittest.TextTestRunner(sys.stdout, verbosity = 2)
    runner.run(test_suite)

def testing_suite():
    return test_suite()




class ClaTests(unittest.TestCase):

    def setUp(self) :
        self.filename = './SCOP/dir.cla.scop.txt_test'

    def testParse(self):
        """Can we parse a CLA file?"""
        f=open(self.filename)
        try: 
            count = 0
            i = Cla.Iterator(f, Cla.Parser())
            while 1 :
                rec = i.next() 
                if rec is None : break
                count +=1
            assert count == 14, "Wrong number of records?!"
        finally:
            f.close()
    
    def testStr(self):
        f = open(self.filename)
        try: 
            p = Cla.Parser()
            i = Cla.Iterator(f)
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
        p = Cla.Parser()

        try:
            rec = p.parse(corruptRec)
            assert False, "Should never get here"
        except ValueError, e :
            pass

    def testRecord(self) :
        recLine = 'd1dan.1\t1dan\tT:,U:91-106\tb.1.2.1\t21953\tcl=48724,cf=48725,sf=49265,fa=49266,dm=49267,sp=49268,px=21953'

        rec = Cla.Parser().parse(recLine)
        assert rec.sid =='d1dan.1'
        assert rec.residues.pdbid =='1dan'
        assert rec.residues.fragments ==(('T','',''),('U','91','106'))
        assert rec.sccs == 'b.1.2.1'
        assert rec.sunid == 21953
        assert rec.hierarchy == [['cl',48724],['cf',48725],['sf',49265],
             ['fa',49266],['dm',49267],['sp',49268],
             ['px',21953]], rec.hierarchy

    def testIndex(self) :
        index = Cla.Index(self.filename)
        
        assert len(index)==14
        assert index.has_key('d4hbia_')

        rec = index['d1hbia_']
        assert rec.sunid == 14996


def test_suite():
    return unittest.makeSuite(ClaTests)


if __name__ == '__main__':
    unittest.main()







