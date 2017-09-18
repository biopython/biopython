# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for online functionality of EBI Search module."""

# Builtins
import unittest

from Bio.EBI.ebisearch import *
from Bio import SeqIO
import os
import sys

sys.path.insert(
    0,
    os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
try:
    import ebisearch
except:
    raise


class ebisearchTests(unittest.TestCase):
    """Tests for EBI Search API."""
    
    def cmp(la, lb):
        """Compare two lists"""
        return all(s in lb for s in la) and all(s in la for s in lb)
    
    
    def test_get_domain_details(self):
        """Test get_domain_details function"""
        domain_details = ebisearch.get_domain_details("metagenomics_runs")
        obs_keys = list(domain_details["domains"][0].keys())
        #exp_keys = ['description', 'id', 'fieldInfos', 'name', 'indexInfos']
        self.assertEqual(len(obs_keys), 5)
    
    
    def test_get_number_of_results(self):
        """Test get_number_of_results function"""
        nb_results = ebisearch.get_number_of_results(
            "metagenomics_runs",
            "experiment_type:(metagenomic)")
        self.assertTrue(nb_results >= 13227)
    
    
    def test_get_domains(self):
        """Test get_domains function"""
        domains = ebisearch.get_domains(verbose=False)
        self.assertTrue("uniref" in domains)
        self.assertTrue('sra-analysis' in domains)
    
    
    def test_print_domain_hierarchy(self):
        """Test print_domain_hierarchy function"""
        ebisearch.print_domain_hierarchy()
        #TODO: print and compare
    
    
    def test_get_fields(self):
        """Test get_fields function"""
        fields = ebisearch.get_fields("metagenomics_runs", verbose=False)
        obs_type = list(fields.keys())
        exp_type = ["searchable", "retrievable", "sortable", "facet", "topterms"]
        self.assertTrue(cmp(obs_type, exp_type))
        self.assertTrue("temperature" in fields["retrievable"])
    
    
    def test_get_domain_search_results(self):
        """Test get_domain_search_results function"""
        res = ebisearch.get_domain_search_results(
            domain="metagenomics_runs",
            query="experiment_type:(metagenomic) AND pipeline_version:(3.0)",
            fields="id,experiment_type",
            size=20)
        fields = list(res[0].keys())
        res_fields = list(res[0]["fields"].keys())
    
        exp_fields = ['source', 'id', 'fields']
        exp_res_fields = ['id', 'experiment_type']
        self.assertTrue(cmp(fields, exp_fields))
        self.assertTrue(cmp(res_fields, exp_res_fields))
        
    
    
    def test_get_all_domain_search_results(self):
        """Test get_all_domain_search_results function"""
        res = ebisearch.get_all_domain_search_results(
            domain="metagenomics_runs",
            query="experiment_type:(metagenomic) AND pipeline_version:(3.0)",
            fields="id,experiment_type")
        self.assertTrue(len(res) >= 2092)
    
    
    def test_get_entries(self):
        """Test get_entries function"""
        ent = ebisearch.get_entries(
            domain="metagenomics_runs",
            entryids="ERR1135279,SRR2135754",
            fields="id,experiment_type")
        fields = list(ent[0].keys())
        res_fields = list(ent[0]["fields"].keys())
        exp_fields = ['source', 'id', 'fields']
        exp_res_fields = ['id', 'experiment_type']
        self.assertTrue(len(ent) == 2)
        self.assertTrue(cmp(fields, exp_fields))
        self.assertTrue(cmp(res_fields, exp_res_fields))
    
    
if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
