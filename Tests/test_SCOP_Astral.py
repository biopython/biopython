# Copyright 2001 by Gavin E. Crooks.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Unit test for Astral"""

import unittest
from StringIO import *

from Bio.SCOP import *

import sys

def run_tests(argv):
    test_suite = testing_suite()
    runner = unittest.TextTestRunner(sys.stdout, verbosity = 2)
    runner.run(test_suite)

def testing_suite():
    test_suite = unittest.TestSuite()
    
    test_loader = unittest.TestLoader()
    test_loader.testMethodPrefix = 'test'
    tests = [AstralTests]
    
    for test in tests:
        cur_suite = test_loader.loadTestsFromTestCase(test)
        test_suite.addTest(cur_suite)
        
    return test_suite


class AstralTests(unittest.TestCase):


    def setUp(self):
        self.scop = Scop(dir_path="SCOP", version="test" )
        self.astral = Astral(scop=self.scop, dir_path="SCOP", version="test")
                

    def testGetSeq(self):
        assert self.astral.getSeqBySid('d3sdha_').data == "AAAAA"
        assert self.astral.getSeqBySid('d4hbib_').data == "KKKKK"

        dom = self.scop.getDomainBySid('d3sdha_')
        assert self.astral.getSeq(dom).data == "AAAAA"

        
        

    def testConstructWithCustomFile(self):
        scop = Scop(dir_path="SCOP", version="test" )
        astral = Astral(scop=scop, astral_file="SCOP/scopseq-test/astral-scopdom-seqres-all-test.fa")
        assert astral.getSeqBySid('d3sdha_').data == "AAAAA"
        assert astral.getSeqBySid('d4hbib_').data == "KKKKK"
                       
         
    def testGetDomainsFromFile(self):
        filename = "SCOP/scopseq-test/astral-scopdom-seqres-sel-gs-bib-20-test.id"
        domains = self.astral.getAstralDomainsFromFile(filename)

        assert len(domains)==3
        assert domains[0].sid == "d3sdha_"
        assert domains[1].sid == "d4hbib_"
        assert domains[2].sid == "d5hbia_"

    def testGetDomainsClustered(self):
        domains1 = self.astral.domainsClusteredById(20)
        assert len(domains1) == 3
        assert domains1[0].sid == "d3sdha_"
        assert domains1[1].sid == "d4hbib_"
        assert domains1[2].sid == "d5hbia_"
                        
        domains2 = self.astral.domainsClusteredByEv(1e-15)
        assert len(domains2) == 1

        #d1 = scop.getDomainBySid("d3sdha_")
        #assert d1.isIn(astral.getHashedDomainsClusteredByPercentId(20))
        #assert d1.isIn(astral.getHashedDomainsClusteredByEv(-15))
        
        
        


if __name__ == '__main__':
    unittest.main()





