"""Tests for the alphafold_db module."""

import os
import tempfile
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

    def test_get_predictions_invalid_accession(self):
        with self.assertRaises(ValueError):
            next(alphafold_db.get_predictions("ZZZZZZ"))

    def test_get_mmcif_file_path_for(self):
        prediction = {
            "cifUrl": "https://alphafold.ebi.ac.uk/files/AF-P00520-F1-model_v4.cif",
        }
        file_path = alphafold_db._get_mmcif_file_path_for(prediction, "test")
        self.assertEqual(file_path, os.path.join("test", "AF-P00520-F1-model_v4.cif"))

    def test_get_file_path_for(self):
        prediction = {
            "pdbUrl": "https://alphafold.ebi.ac.uk/files/AF-P00520-F1-model_v4.pdb",
        }
        file_path = alphafold_db._get_file_path_for(prediction, "pdbUrl", "test")
        self.assertEqual(file_path, os.path.join("test", "AF-P00520-F1-model_v4.pdb"))

    def test_download_pdb_for(self):
        prediction = next(alphafold_db.get_predictions("P00520"))
        with tempfile.TemporaryDirectory() as directory:
            file_path = alphafold_db.download_pdb_for(prediction, directory)
            self.assertTrue(file_path.endswith(".pdb"))
            self.assertTrue(os.path.isfile(file_path))

    def test_download_pae_for(self):
        prediction = next(alphafold_db.get_predictions("P00520"))
        with tempfile.TemporaryDirectory() as directory:
            file_path = alphafold_db.download_pae_for(prediction, directory)
            self.assertTrue(os.path.isfile(file_path))
