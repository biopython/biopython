"""Tests for the alphafold_db module."""

import unittest

import requires_internet

from Bio.PDB import alphafold_db

requires_internet.check()


class AlphafoldDBTests(unittest.TestCase):
    def test_get_predictions(self):
        predictions = alphafold_db.get_predictions("P00520")

        for prediction in predictions:
            self.assertIsInstance(prediction, dict)
            self.assertGreater(len(prediction), 0)

    def test_get_mmcif_file_path_for(self):
        prediction = {
            "cifUrl": "https://alphafold.ebi.ac.uk/files/AF-P00520-F1-model_v4.cif",
        }
        file_path = alphafold_db._get_mmcif_file_path_for(prediction, "test")
        self.assertEqual(file_path, "test/AF-P00520-F1-model_v4.cif")
