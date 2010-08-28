# Copyright 2001 by Gavin E. Crooks.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Unit test for Astral"""

import unittest
from Bio.SCOP import *



class AstralTests(unittest.TestCase):


    def setUp(self):
        self.scop = Scop(dir_path="SCOP", version="test")
        self.astral = Astral(scop=self.scop, dir_path="SCOP", version="test")
                

    def testGetSeq(self):
        self.assertEqual(self.astral.getSeqBySid('d3sdha_').tostring(), "AAAAA")
        self.assertEqual(self.astral.getSeqBySid('d4hbib_').tostring(), "KKKKK")

        dom = self.scop.getDomainBySid('d3sdha_')
        self.assertEqual(self.astral.getSeq(dom).tostring(), "AAAAA")

        
        

    def testConstructWithCustomFile(self):
        scop = Scop(dir_path="SCOP", version="test")
        astral = Astral(scop=scop, astral_file="SCOP/scopseq-test/astral-scopdom-seqres-all-test.fa")
        self.assertEqual(astral.getSeqBySid('d3sdha_').tostring(), "AAAAA")
        self.assertEqual(astral.getSeqBySid('d4hbib_').tostring(), "KKKKK")
                       
         
    def testGetDomainsFromFile(self):
        filename = "SCOP/scopseq-test/astral-scopdom-seqres-sel-gs-bib-20-test.id"
        domains = self.astral.getAstralDomainsFromFile(filename)

        self.assertEqual(len(domains), 3)
        self.assertEqual(domains[0].sid, "d3sdha_")
        self.assertEqual(domains[1].sid, "d4hbib_")
        self.assertEqual(domains[2].sid, "d5hbia_")

    def testGetDomainsClustered(self):
        domains1 = self.astral.domainsClusteredById(20)
        self.assertEqual(len(domains1), 3)
        self.assertEqual(domains1[0].sid, "d3sdha_")
        self.assertEqual(domains1[1].sid, "d4hbib_")
        self.assertEqual(domains1[2].sid, "d5hbia_")
                        
        domains2 = self.astral.domainsClusteredByEv(1e-15)
        self.assertEqual(len(domains2), 1)

        #d1 = scop.getDomainBySid("d3sdha_")
        #self.assertEqual(d1.isIn(astral.getHashedDomainsClusteredByPercentId(20))
        #self.assertEqual(d1.isIn(astral.getHashedDomainsClusteredByEv(-15))
        
        
        


if __name__=='__main__':
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
