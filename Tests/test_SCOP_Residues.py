# Copyright 2001 by Gavin E. Crooks.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Unit test for Residues"""

import unittest
from Bio.SCOP.Residues import *




class ResiduesTests(unittest.TestCase):
    res = (
        ( "-",           () ),
        ( "A:",          (("A", "", ""),) ),
        ( "1:",          (("1", "", ""),) ),
        ( "1-100",       (("", "1", "100"),)  ),
        ( "B:1-101",     (("B",   "1" ,"101"),) ),
        ( "1:1a-100a",   (("1", "1a", "100a"),) ),
        ( "a:-100a--1a", (("a", "-100a", "-1a"),) ),
        ( "-1-100",      (("", "-1", "100"),) ),
        ( "-1-100",      (("", "-1", "100"),) ),
        ( "A:12-19,A:23-25", (("A","12","19"),("A","23","25")) ),
        ( "12-19,1:23-25", (("","12","19"),("1","23","25")) ),
        ( "0-1,1:-1a-25a,T:", (("","0","1"),("1","-1a","25a"),("T","","")) ),
        )


    def testParse(self):
        for loc in self.res:
            r = Residues(loc[0])
            self.assertEqual(r.fragments, loc[1])

    def testStr(self):
        for loc in self.res:
            r = Residues(loc[0])
            self.assertEqual(str(r), loc[0])

    def testAstralParse(self):
        """Test if we can parse residue subsets enclosed in brackets"""
        for loc in self.res:
            r = Residues("("+loc[0]+")")
            self.assertEqual(r.fragments, loc[1])

    def testPdbId(self):
        pdbid ="1ddf"
        for loc in self.res:
            r = Residues("\t 1ddf \t"+loc[0]+"\t\n\n\n")
            self.assertEqual(r.pdbid, pdbid)
            self.assertEqual(str(r), pdbid+" "+loc[0])

            r = Residues(pdbid+" "+loc[0])
            self.assertEqual(r.pdbid, pdbid)
            self.assertEqual(str(r), pdbid+" "+loc[0])


            r = Residues("104l A:112-113")
            self.assertEqual(r.pdbid, "104l")
            self.assertEqual(r.fragments, (('A', '112', '113'),))

            
    def testJustPdbId(self):
        r = Residues("1sds")
        self.assertEqual(r.pdbid, "1sds")
        self.assertFalse(r.fragments)


    def testParseError(self):
        self.assertRaises(ValueError, Residues, "09324923423hh./;,.389")



if __name__=='__main__':
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
