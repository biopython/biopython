"""Testing code for Bio-Registries.

Goals:
    Make sure that all retrieval is working as expected.
"""
import requires_internet

import os
import sys
import unittest
import StringIO

from Bio import GenBank

def run_tests(argv):
    test_suite = testing_suite()
    runner = unittest.TextTestRunner(sys.stdout, verbosity = 2)
    runner.run(test_suite)

def testing_suite():
    """Generate the suite of tests.
    """
    test_suite = unittest.TestSuite()

    test_loader = unittest.TestLoader()
    test_loader.testMethodPrefix = 't_'
    tests = [GenBankRetrievalTest]
    
    for test in tests:
        cur_suite = test_loader.loadTestsFromTestCase(test)
        test_suite.addTest(cur_suite)

    return test_suite

class GenBankRetrievalTest(unittest.TestCase):
    """Test convience functionality in Bio.GenBank for retrieval.
    """
    def setUp(self):
        self.gb_protein_ids = ["NP_058835", "8393944"]
        self.gb_nuc_ids = ["NM_017139", "8393943"]
        self.search_term = "Opuntia"

    def t_search_for(self):
        """Search for ids using the GenBank search_for function.
        """
        ids = GenBank.search_for(self.search_term, database = 'protein')
        assert len(ids) >= 9 # 9 in GenBank right now

    def t_download_many(self):
        """Retrieve sequences using the GenBank download_many function.
        """
        nuc_results = GenBank.download_many(self.gb_nuc_ids, 'nucleotide')
        iterator = GenBank.Iterator(nuc_results)
        num = 0
        for rec in iterator:
            num += 1
        assert num == 2

        prot_results = GenBank.download_many(self.gb_protein_ids, 'protein')
        iterator = GenBank.Iterator(prot_results)
        num = 0
        for rec in iterator:
            num += 1
        assert num == 2

if __name__ == "__main__":
    sys.exit(run_tests(sys.argv))
