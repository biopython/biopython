"""Testing code for Bio-Registries.

Goals:
    Make sure that all retrieval is working as expected.
"""
import requires_internet
import os
import sys
import unittest

from Bio import db

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
    tests = [BaseRetrievalTest]
    
    for test in tests:
        cur_suite = test_loader.loadTestsFromTestCase(test)
        test_suite.addTest(cur_suite)

    return test_suite

class BaseRetrievalTest(unittest.TestCase):
    def setUp(self):
        self.embl_id = "E33091"
        self.sp_id = "MTHC_DROME"
        self.gb_protein_id = "NP_058835"
        self.gb_nuc_id = "NM_017139"

    def tearDown(self):
        pass

    def t_swissprot(self):
        """Retrieval of Swissprot records from various sources.
        """
        sp_db = db["swissprot"]
        for swissprot_location in sp_db.objs:
            handle = swissprot_location[self.sp_id]
            first_line = handle.read(20)
            assert self.sp_id in first_line, \
                    (swissprot_location.name, first_line)

    def t_embl(self):
        """Retrieval of EMBL records from various sources.
        """
        embl_db = db["embl"]
        for embl_location in embl_db.objs:
            #if embl_location.name.find("xembl") == -1:
            handle = embl_location[self.embl_id]
            first_line = handle.read(20)
            assert self.embl_id in first_line, (embl_location.name, first_line)

    def t_embl_xml(self):
        """Retrieval of XML EMBL records in BSML format.
        """
        embl_db = db["embl-xml"]
        for embl_location in embl_db.objs:
            handle = embl_location[self.embl_id]
            # xml output from xembl
            first_line = handle.read(400)
            assert first_line.find("<?xml version") == 0, first_line
            assert self.embl_id in first_line, (embl_location.name, first_line)

    def t_genbank(self):
        """Retrieval of GenBank from various sources.
        """
        for gb_id, gb_location in \
                [(self.gb_protein_id, db["protein-genbank-eutils"]),
                 (self.gb_nuc_id, db["nucleotide-genbank-eutils"])]:
            handle = gb_location[gb_id]
            first_line = handle.read(25)
            assert gb_id in first_line, (gb_location.name, first_line)

if __name__ == "__main__":
    sys.exit(run_tests(sys.argv))
