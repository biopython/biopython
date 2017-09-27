# Copyright 2017 by Berenice Batut (berenice.batut@gmail.com). All rights reserved.
# Revision copyright 2017 by Francesco Gastaldello. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Tests for online functionality of EBI Search module."""

import unittest
import ebisearch


class ebisearchTests(unittest.TestCase):
    """Tests for EBI Search API."""

    def test_get_domain_details(self):
        """Test get_domain_details function"""
        domain_details = ebisearch.get_domain_details("metagenomics_runs")
        obs_keys = list(domain_details["domains"][0].keys())
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

    def test_get_fields(self):
        """Test get_fields function"""
        fields = ebisearch.get_fields("metagenomics_runs", verbose=False)
        obs_type = list(fields.keys())
        exp_type = ['facet', 'retrievable', 'searchable', 'sortable', 'topterms']
        self.assertEqual(sorted(obs_type), exp_type)
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
        exp_fields = ['fields', 'id', 'source']
        exp_res_fields = ['experiment_type', 'id']
        self.assertEqual(sorted(fields), exp_fields)
        self.assertEqual(sorted(res_fields), exp_res_fields)

    def test_check_domain(self):
        """Test check_domain function"""
        with self.assertRaises(ValueError) as context:
            ebisearch.check_domain("emm")
        self.assertTrue("The domain does not correspond to the id of a known domain in EBI. The list of EBI domains and their id can be accessed with get_domains" in str(context.exception))

    def test_get_searchable_fields(self):
        """Test get_searchable_fields function"""
        self.assertEqual(len(ebisearch.get_searchable_fields("metagenomics_runs", verbose=False)), 24)

    def test_get_sortable_fields(self):
        """Test get_Sortable_fields function"""
        self.assertEqual(len(ebisearch.get_sortable_fields("metagenomics_runs", verbose=False)), 0)

    def test_get_facet_fields(self):
        """Test get_facet_fields function"""
        self.assertEqual(len(ebisearch.get_facet_fields("metagenomics_runs", verbose=False)), 7)

    def test_get_topterms_fields(self):
        """Test get_topterms_fields function"""
        self.assertEqual(len(ebisearch.get_topterms_fields("metagenomics_runs", verbose=False)), 0)

    def test_check_retrievable_fields(self):
        """Test check_retrievable_fields function"""
        with self.assertRaises(ValueError) as context:
            ebisearch.check_retrievable_fields("protein_domain", "metagenomics_runs")
        self.assertTrue("The field protein_domain does not correspond to a retrievable field for the domain. The list of retrievable fields for a domain can be accessed with get_retrievable_fields" in str(context.exception))

    def test_check_sortable_field(self):
        """Test check_sortable_fields function"""
        with self.assertRaises(ValueError) as context:
            ebisearch.check_sortable_field("gene_name", "metagenomics_runs")
        self.assertTrue("The field gene_name is not a sortable field for the domain metagenomics_runsThe list of sortable fields for the domain can be accessed with get_sortable_fields" in str(context.exception))

    def test_check_facet_field(self):
        """Test check_facet_field function"""
        with self.assertRaises(ValueError) as context:
            ebisearch.check_facet_field("gene_name", "metagenomics_runs")
        self.assertTrue("The field gene_name is not a facet field for the domain metagenomics_runsThe list of sortable fields for the domain can be accessed with get_facet_fields" in str(context.exception))

    def test_check_topterms_field(self):
        """Test check_topterms_field function"""
        with self.assertRaises(ValueError) as context:
            ebisearch.check_topterms_field("gene_name", "metagenomics_runs")
        self.assertTrue("The field gene_name is not a topterms field for the domain metagenomics_runsThe list of sortable fields for the domain can be accessed with get_topterms_fields" in str(context.exception))

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
        exp_fields = ['fields', 'id', 'source']
        exp_res_fields = ['experiment_type', 'id']
        self.assertTrue(len(ent) == 2)
        self.assertEqual(sorted(fields), exp_fields)
        self.assertEqual(sorted(res_fields), exp_res_fields)


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
