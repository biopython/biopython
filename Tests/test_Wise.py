#!/usr/bin/env python2.3

__version__ = "$Revision: 1.1 $"

import cStringIO
import doctest, unittest
import sys

import requires_wise

from Bio import Wise

class TestWiseDryRun(unittest.TestCase):
    def setUp(self):
        self.old_stdout = sys.stdout
        sys.stdout = cStringIO.StringIO()
        
    def test_dnal(self):
        Wise.align(["dnal"], ("seq1.fna", "seq2.fna"), kbyte=100000, dry_run=True)
        self.assert_(sys.stdout.getvalue().startswith("dnal -kbyte 100000 seq1.fna seq2.fna > /tmp/tmp"))

    def test_psw(self):
        Wise.align(["psw"], ("seq1.faa", "seq2.faa"), dry_run=True)
        self.assert_(sys.stdout.getvalue().startswith("psw -kbyte 300000 seq1.faa seq2.faa > /tmp/tmp"))

    def tearDown(self):
        sys.stdout = self.old_stdout

class TestWise(unittest.TestCase):
    def test_align(self):
        temp_file = Wise.align(["dnal"], ("Wise/human_114_g01_exons.fna_01", "Wise/human_114_g02_exons.fna_01"), kbyte=100000, force_type="DNA", quiet=True)
        self.assertEqual(temp_file.readline().rstrip(), "ENSG00000172135   AGGGAAAGCCCCTAAGCTC--CTGATCTATGCTGCATCCAGTTTGCAAAGTGGGGTCCC")

def run_tests(argv):
    test_suite = testing_suite()
    runner = unittest.TextTestRunner(sys.stdout, verbosity = 2)
    runner.run(test_suite)

def testing_suite():
    """Generate the suite of tests.
    """
    unittest_suite = unittest.TestSuite()

    test_loader = unittest.TestLoader()
    test_loader.testMethodPrefix = 'test_'
    tests = [TestWiseDryRun, TestWise]
    
    for test in tests:
        cur_suite = test_loader.loadTestsFromTestCase(test)
        unittest_suite.addTest(cur_suite)

    doctest_suite = doctest.DocTestSuite(Wise)

    big_suite = unittest.TestSuite((unittest_suite, doctest_suite))

    return big_suite

if __name__ == "__main__":
    sys.exit(run_tests(sys.argv))
