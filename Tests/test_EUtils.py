"""Testing code for the EUtils interface for retrieving from Entrez.
"""
import requires_internet

import os
import sys
import unittest
import StringIO

from Bio import EUtils
from Bio.EUtils import HistoryClient

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
    tests = [HistoryClientTest]
    
    for test in tests:
        cur_suite = test_loader.loadTestsFromTestCase(test)
        test_suite.addTest(cur_suite)

    return test_suite

class HistoryClientTest(unittest.TestCase):
    """Tests for the EUtils HistoryClient interface.
    """
    def setUp(self):
        self.client = HistoryClient.HistoryClient()

    def tearDown(self):
        pass

    def t_simple_search(self):
        """Simple text search from PubMed.
        """
        search = self.client.search("narrow sheath")
        assert len(search) > 0

    def t_retrieve_search(self):
        """Search and retrieval of PubMed records.
        """
        num_to_fetch = 3
        search = self.client.search("opuntia ficus indica")
        results = search[:num_to_fetch].efetch(retmode = "text",
                                               rettype = "docsum").read()
        num_records = 0
        for line in results.split("\n"):
            if line.find("PMID:") == 0:
                num_records += 1
        assert num_records == num_to_fetch, results

    def t_retrieve_many(self):
        """Retrieving more than 500 sequences from nucleotide database.
        """
        search = self.client.search("Arabidopsis", db = "nucleotide",
                daterange = EUtils.DateRange("1995/01/01", "1995/04/15",
                "mdat"))
        results = search.efetch(retmode = "text", rettype = "seqid").read()
        num_records = 0
        for line in results.split("\n"):
            if line.find("::= gi ") > 0:
                num_records += 1
        assert num_records > 500, num_records

    def t_retrieve_by_accession(self):
        """Retrieve a single nucleotide record by accession.
        """
        search = self.client.search("AA759046", db="nucleotide")

if __name__ == "__main__":
    sys.exit(run_tests(sys.argv))
