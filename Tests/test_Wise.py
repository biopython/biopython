#!/usr/bin/env python
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

__version__ = "$Revision: 1.11 $"

import cStringIO
import doctest, unittest
import sys

if 'requires_wise' in sys.modules:
    del sys.modules['requires_wise']
import requires_wise

from Bio import Wise

class TestWiseDryRun(unittest.TestCase):
    def setUp(self):
        self.old_stdout = sys.stdout
        sys.stdout = cStringIO.StringIO()
        
    def test_dnal(self):
        """Call dnal, and do a trivial check on its output."""
        Wise.align(["dnal"], ("seq1.fna", "seq2.fna"), kbyte=100000, dry_run=True)
        #If test output is redirected to a file, the wrapper adds -quiet
        output = sys.stdout.getvalue().replace(" -quiet ", " ")
        self.assertTrue(output.startswith("dnal -kbyte 100000 seq1.fna seq2.fna"), output[:200])

    def test_psw(self):
        """Call psw, and do a trivial check on its output."""
        Wise.align(["psw"], ("seq1.faa", "seq2.faa"), dry_run=True, kbyte=4)
        #If test output is redirected to a file, the wrapper adds -quiet
        output = sys.stdout.getvalue().replace(" -quiet ", " ")
        self.assertTrue(output.startswith("psw -kbyte 4 seq1.faa seq2.faa"), output[:200])

    def tearDown(self):
        sys.stdout = self.old_stdout

class TestWise(unittest.TestCase):
    def test_align(self):
        """Call dnal with optional arguments, and do a trivial check on the output."""
        temp_file = Wise.align(["dnal"], ("Wise/human_114_g01_exons.fna_01", "Wise/human_114_g02_exons.fna_01"), kbyte=100000, force_type="DNA", quiet=True)
        line = temp_file.readline().rstrip()
        if line == "Score 114":
            #Wise 2.4.1 includes a score line, even in quiet mode, ignore this
            line = temp_file.readline().rstrip()
        if line == "ENSG00000172135   AGGGAAAGCCCCTAAGCTC--CTGATCTATGCTGCATCCAGTTTGCAAAGTGGGGTCCC":
            #This is what we expect from wise 2.2.0 (and earlier)
            pass
        elif line == "ENSG00000172135   AGGGAAAGCCCCTAAGCTC--CTGATCTATGCTGCATCCAGTTTGCAAAG-TGGGGTCC":
            #This is what we expect from wise 2.4.1
            pass
        else:
            #Bad!
            self.assertTrue(False, line)


if __name__ == "__main__":
    unittest_suite = unittest.TestLoader().loadTestsFromName("test_Wise")
    doctest_suite = doctest.DocTestSuite(Wise)
    suite = unittest.TestSuite((unittest_suite, doctest_suite))
    runner = unittest.TextTestRunner(sys.stdout, verbosity = 2)
    runner.run(suite)
