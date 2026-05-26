"""Tests for the alphafold_db module."""

import os
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
        self.assertEqual(file_path, os.path.join("test", "AF-P00520-F1-model_v4.cif"))

    def test_list_proteomes(self):
        proteomes = alphafold_db.list_proteomes()

        self.assertIsInstance(proteomes, list)
        self.assertGreater(len(proteomes), 0)
        species_proteomes = [p for p in proteomes if p.get("type") in ("proteome", "global_health")]
        self.assertGreater(len(species_proteomes), 0)

        for proteome in species_proteomes:
            self.assertIn("species", proteome)
            self.assertIsInstance(proteome, dict)
            self.assertIn("species", proteome)
            self.assertIn("archive_name", proteome)
            self.assertIn("reference_proteome", proteome)
            self.assertIn("num_predicted_structures", proteome)
            self.assertIn("size_bytes", proteome)

    def test_download_proteome_invalid_species(self):
        with self.assertRaises(ValueError):
            alphafold_db.download_proteome("NOTASPECIES123")

    def test_download_proteome_match_by_common_name(self):
        proteomes = alphafold_db.list_proteomes()
        # find smallest proteome to keep test fast

        smallest = min(proteomes, key=lambda p: p["size_bytes"])
        # just test matching logic, not actual download

        match = None
        species_lower = smallest["common_name"].lower()

        for proteome in proteomes:
            if species_lower in (
                proteome["species"].lower(),
                proteome["common_name"].lower(),
                proteome["reference_proteome"].lower(),
            ):
                match = proteome
                break
        self.assertIsNotNone(match)
        self.assertEqual(match["reference_proteome"], smallest["reference_proteome"])

