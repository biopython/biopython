"""Testing code for Bio-Registries.

Goals:
    Make sure that all retrieval is working as expected.
"""
import requires_internet

import os
import sys
import unittest
import StringIO

from Bio import db
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
    tests = [BaseRetrievalTest, GenBankRetrievalTest]
    
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
        self.interpro_id = "IPR000836"
        self.pdb_id = "1BG2"
        self.prodoc_id = "PDOC00100"
        self.prosite_id = "PS00615"
        self.pubmed_id = "14734308"

    def tearDown(self):
        pass

    def t_swissprot(self):
        """Retrieval of Swissprot records from various sources.
        """
        sp_db = db["swissprot"]
        for swissprot_location in sp_db.objs:
            handle = swissprot_location[self.sp_id]
            first_line = handle.read(20)
            assert first_line.find(self.sp_id) >= 0, \
                    (swissprot_location.name, first_line)

    def t_embl(self):
        """Retrieval of EMBL records from various sources.
        """
        embl_db = db["embl"]
        for embl_location in embl_db.objs:
            #if embl_location.name.find("xembl") == -1:
            handle = embl_location[self.embl_id]
            first_line = handle.read(20)
            assert first_line.find(self.embl_id) >= 0, \
                    (embl_location.name, first_line)

    def t_embl_xml(self):
        """Retrieval of XML EMBL records in BSML format.
        """
        embl_db = db["embl-xml"]
        for embl_location in embl_db.objs:
            handle = embl_location[self.embl_id]
            # xml output from xembl
            first_line = handle.read(400)
            assert first_line.find("<?xml version") == 0, first_line
            assert first_line.find(self.embl_id) >= 0, \
                    (embl_location.name, first_line)

    def t_genbank(self):
        """Retrieval of GenBank from various sources.
        """
        for gb_id, gb_db in \
                [(self.gb_protein_id, db["genbank-protein"]),
                 (self.gb_nuc_id, db["genbank-nucleotide"])]:
            for gb_location in gb_db.objs:
                handle = gb_location[gb_id]
                first_line = handle.read(25)
                assert first_line.find(gb_id) >= 0,\
                        (gb_location.name, first_line)

    def t_interpro(self):
        """Retrieval of InterPro data from EBI.
        """
        interpro_db = db["interpro"]
        for interpro_loc in interpro_db.objs:
            handle = interpro_loc[self.interpro_id]
            data = handle.read()
            assert data.find(self.interpro_id) >= 0, \
                    (interpro_location.name, data)

    def t_pdb(self):
        """Retrieval of PDB data from various locations.
        """
        pdb_db = db["pdb"]
        for pdb_loc in pdb_db.objs:
            handle = pdb_loc[self.pdb_id]
            first_line = handle.read(80)
            assert first_line.find(self.pdb_id) >= 0, \
                    (pdb_loc.name, first_line)

    def t_prodoc(self):
        """Retrieval of Prodoc data from various locations.
        """
        prodoc_db = db["prodoc"]
        for prodoc_loc in prodoc_db.objs:
            handle = prodoc_loc[self.prodoc_id]
            first_line = handle.read(20)
            assert first_line.find(self.prodoc_id) >= 0, \
                    (prodoc_loc.name, first_line)

    def t_prosite(self):
        """Retrieval of Prosite data from various locations.
        """
        prosite_db = db["prosite"]
        for prosite_loc in prosite_db.objs:
            handle = prosite_loc[self.prosite_id]
            first_part = handle.read(80)
            assert first_part.find(self.prosite_id) >= 0, \
                    (prosite_loc.name, first_part)

    def t_medline(self):
        """Retrieval of Medline data from various locations.
        """
        medline_db = db["medline"]
        for medline_loc in medline_db.objs:
            handle = medline_loc[self.pubmed_id]
            first_part = handle.read(30)
            assert first_part.find(self.pubmed_id) >= 0, \
                    (medline_loc.name, first_part)

    def t_fasta(self):
        """Retrieval of FASTA protein and nucleotide sequences.
        """
        fasta_db = db["fasta"]
        for fasta_loc in fasta_db.objs:
            for fasta_id in [self.gb_protein_id, self.gb_nuc_id]:
                handle = fasta_loc[fasta_id]
                first_part = handle.read(80)
                assert first_part.find(fasta_id) >= 0, \
                        (fasta_loc.name, fasta_id, first_part)

class GenBankRetrievalTest(unittest.TestCase):
    """Test convience functionality in Bio.GenBank for retrieval.
    """
    def setUp(self):
        self.gb_protein_ids = ["NP_058835", "8393944"]
        self.gb_nuc_ids = ["NM_017139", "8393943"]
        self.search_term = "Opuntia"

    def t_ncbi_dictionary(self):
        """Retrieve sequences using a GenBank dictionary-like interface.
        """
        for ids, database_type in [(self.gb_protein_ids, "protein"),
                                   (self.gb_nuc_ids, "nucleotide")]:
            for format in ["genbank", "fasta"]:
                ncbi_dict = GenBank.NCBIDictionary(database_type, format)
                for gb_id in ids:
                    ncbi_dict[gb_id]

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
