"""Tests for the UniProt module."""

import unittest
from itertools import islice

import requires_internet

from Bio import UniProt

requires_internet.check()


class SearchTests(unittest.TestCase):
    def test_search_result_count(self):
        query = "(organism_id:2697049) AND (reviewed:true)"
        results = UniProt.search(query)
        search_result_count = len(results)

        self.assertIsInstance(search_result_count, int)
        self.assertGreaterEqual(search_result_count, 1)
        self.assertEqual(search_result_count, len(list(results)))

    def test_search_all_fields(self):
        query = "Insulin AND (reviewed:true)"
        results = list(islice(UniProt.search(query, batch_size=50), 50))
        self.assertEqual(len(results), 50)

        for result in results:
            self.assertIsInstance(result, dict)
            self.assertGreater(len(result), 0)

    def test_search_with_fields(self):
        query = "Insulin AND (reviewed:true)"
        fields = ["id", "accession"]
        results = list(islice(UniProt.search(query, fields=fields, batch_size=50), 50))
        self.assertEqual(len(results), 50)
        expected_keys = {"entryType", "primaryAccession", "uniProtkbId"}

        for result in results:
            self.assertIsInstance(result, dict)
            self.assertGreater(len(result), 0)
            self.assertEqual(set(result.keys()), expected_keys)

    def test_search_results_slicing(self):
        query = "Insulin AND (reviewed:true)"
        results = UniProt.search(query, batch_size=25)
        results_sliced = results[:50]
        self.assertEqual(50, len(results_sliced))
        results_list = list(islice(results, 50))
        self.assertEqual(50, len(results_list))
        self.assertEqual(results_sliced, results_list)

        # We need to fetch the results again to reset the iterator
        results = UniProt.search(query, batch_size=25)
        results_sliced = results[49:30:-1]
        self.assertEqual(19, len(results_sliced))
        results_list = list(islice(results, 50))[49:30:-1]
        self.assertEqual(19, len(results_list))
        self.assertEqual(results_sliced, results_list)

        results_sliced = results[0:0]
        self.assertEqual(0, len(results_sliced))

        # This query should return less than 500 results.
        query = "(organism_id:2697049) AND (reviewed:true)"
        results = UniProt.search(query)
        results_sliced = results[-10:-5]
        self.assertEqual(5, len(results_sliced))
        results_list = list(results)[-10:-5]
        self.assertEqual(5, len(results_list))
        self.assertEqual(results_sliced, results_list)

        results = UniProt.search(query)
        results_sliced = results[-1_000_001:1_000_000]
        self.assertEqual(len(results), len(results_sliced))
        results_list = list(results)[-1_000_001:1_000_000]
        self.assertEqual(len(results_list), len(results))
        self.assertEqual(results_sliced, results_list)

    def test_search_results_indexing(self):
        query = "(organism_id:2697049) AND (reviewed:true)"
        results = UniProt.search(query)
        self.assertIsInstance(results[0], dict)
        self.assertIsInstance(results[-1], dict)
        with self.assertRaises(IndexError):
            _ = results[1_000_000]
        with self.assertRaises(IndexError):
            _ = results[-1_000_000]
