"""Tests for the alphafold_db module."""

import os
import tempfile
import unittest

import requires_internet

from Bio.PDB import alphafold_db
from Bio.PDB import AFDBList

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


class AFDBListTests(unittest.TestCase):
    def test_retrieve_pdb_file(self):
        """Download a single AlphaFold structure by UniProt ID."""
        with tempfile.TemporaryDirectory() as tmpdir:
            afdb = AFDBList(pdb=tmpdir, verbose=False)
            path = afdb.retrieve_pdb_file("P00520", pdir=tmpdir)
            self.assertIsNotNone(path)
            self.assertTrue(os.path.exists(path))
            self.assertTrue(path.endswith(".cif"))

    def test_retrieve_pdb_file_pdb_format(self):
        """Download as PDB format."""
        with tempfile.TemporaryDirectory() as tmpdir:
            afdb = AFDBList(pdb=tmpdir, verbose=False)
            path = afdb.retrieve_pdb_file("P00520", pdir=tmpdir, file_format="pdb")
            self.assertIsNotNone(path)
            self.assertTrue(os.path.exists(path))
            self.assertTrue(path.endswith(".pdb"))

    def test_retrieve_pdb_file_invalid_format(self):
        """Invalid file_format raises ValueError."""
        afdb = AFDBList(verbose=False)
        with self.assertRaises(ValueError):
            afdb.retrieve_pdb_file("P00520", file_format="xyz")
