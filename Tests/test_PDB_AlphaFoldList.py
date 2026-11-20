"""Tests for Bio.PDB.AlphaFoldList module."""

import os
import shutil
import tempfile
import unittest

import requires_internet

from Bio.PDB.AlphaFoldList import AlphaFoldList

requires_internet.check()


class AlphaFoldListTests(unittest.TestCase):
    """Test the AlphaFoldList module."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.afl = AlphaFoldList(local_dir=self.temp_dir, verbose=False)

    def tearDown(self):
        """Clean up test fixtures."""
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)

    def test_get_structure_metadata(self):
        """Test retrieving metadata from AlphaFold API."""
        metadata = self.afl.get_structure_metadata("P00520")

        self.assertIsNotNone(metadata)
        self.assertIsInstance(metadata, list)
        self.assertGreater(len(metadata), 0)

        # Check first entry has expected fields
        entry = metadata[0]
        self.assertIn("entryId", entry)
        self.assertIn("uniprotAccession", entry)
        self.assertIn("organismScientificName", entry)
        self.assertIn("latestVersion", entry)
        self.assertIn("pdbUrl", entry)
        self.assertIn("cifUrl", entry)

        # Verify UniProt ID matches
        self.assertTrue(entry["uniprotAccession"].startswith("P00520"))

    def test_get_structure_metadata_invalid(self):
        """Test handling of invalid UniProt ID."""
        metadata = self.afl.get_structure_metadata("INVALID123")
        self.assertIsNone(metadata)

    def test_retrieve_alphafold_structure_pdb(self):
        """Test downloading a structure in PDB format."""
        result = self.afl.retrieve_alphafold_structure(
            "P00520", file_format="pdb", overwrite=True
        )

        self.assertIsNotNone(result)
        self.assertTrue(os.path.exists(result))
        self.assertTrue(result.endswith(".pdb"))

        # Check file is not empty and has PDB content
        self.assertGreater(os.path.getsize(result), 1000)
        with open(result, "r") as f:
            content = f.read(200)
            # PDB files should have HEADER, ATOM, or MODEL lines
            self.assertTrue(
                any(keyword in content for keyword in ["HEADER", "ATOM", "MODEL"])
            )

    def test_retrieve_alphafold_structure_mmcif(self):
        """Test downloading a structure in mmCIF format."""
        result = self.afl.retrieve_alphafold_structure(
            "P00520", file_format="mmCif", overwrite=True
        )

        self.assertIsNotNone(result)
        self.assertTrue(os.path.exists(result))
        self.assertTrue(result.endswith(".cif"))

        # Check file is not empty
        self.assertGreater(os.path.getsize(result), 1000)

    def test_retrieve_with_pae(self):
        """Test downloading structure with PAE file."""
        result = self.afl.retrieve_alphafold_structure(
            "P00520", file_format="pdb", include_pae=True, overwrite=True
        )

        self.assertIsNotNone(result)
        self.assertTrue(os.path.exists(result))

        # Check if PAE JSON file exists
        base_dir = os.path.dirname(result)
        pae_files = [
            f for f in os.listdir(base_dir) if "predicted_aligned_error" in f
        ]
        self.assertGreater(len(pae_files), 0)

    def test_retrieve_with_confidence(self):
        """Test downloading structure with confidence file."""
        result = self.afl.retrieve_alphafold_structure(
            "P00520", file_format="pdb", include_confidence=True, overwrite=True
        )

        self.assertIsNotNone(result)
        self.assertTrue(os.path.exists(result))

        # Check if confidence JSON file exists
        base_dir = os.path.dirname(result)
        conf_files = [f for f in os.listdir(base_dir) if "confidence" in f]
        self.assertGreater(len(conf_files), 0)

    def test_overwrite_protection(self):
        """Test that existing files are not overwritten by default."""
        # Download once
        result1 = self.afl.retrieve_alphafold_structure(
            "P00520", file_format="pdb"
        )
        self.assertIsNotNone(result1)

        # Get modification time
        mtime1 = os.path.getmtime(result1)

        # Try to download again without overwrite
        result2 = self.afl.retrieve_alphafold_structure(
            "P00520", file_format="pdb", overwrite=False
        )

        # Should return same file without modifying it
        self.assertEqual(result1, result2)
        mtime2 = os.path.getmtime(result2)
        self.assertEqual(mtime1, mtime2)

    def test_get_available_proteomes(self):
        """Test retrieving list of available proteomes."""
        proteomes = self.afl.get_available_proteomes()

        self.assertIsNotNone(proteomes)
        self.assertIsInstance(proteomes, list)
        self.assertGreater(len(proteomes), 0)

        # Check structure of proteome entries
        for proteome in proteomes[:5]:  # Check first 5
            self.assertIn("filename", proteome)
            self.assertIn("organism", proteome)
            self.assertTrue(proteome["filename"].endswith(".tar"))

        # Verify some known organisms are present
        organism_codes = [p["organism"] for p in proteomes]
        common_organisms = ["HUMAN", "MOUSE", "ECOLI", "YEAST"]
        found = [org for org in common_organisms if org in organism_codes]
        self.assertGreater(len(found), 0, "Should find at least one common organism")

    def test_flat_tree_option(self):
        """Test flat directory tree option."""
        self.afl.flat_tree = True
        result = self.afl.retrieve_alphafold_structure(
            "P00520", file_format="pdb", overwrite=True
        )

        self.assertIsNotNone(result)
        # With flat_tree, file should be directly in temp_dir
        self.assertEqual(os.path.dirname(result), self.temp_dir)

    def test_directory_tree_structure(self):
        """Test hierarchical directory tree structure."""
        self.afl.flat_tree = False
        result = self.afl.retrieve_alphafold_structure(
            "P00520", file_format="pdb", overwrite=True
        )

        self.assertIsNotNone(result)
        # With tree structure, file should be in subdirectory
        self.assertNotEqual(os.path.dirname(result), self.temp_dir)
        # Directory should be based on entry ID characters
        self.assertTrue(os.path.exists(os.path.dirname(result)))

    def test_invalid_file_format(self):
        """Test handling of invalid file format."""
        with self.assertRaises(ValueError):
            self.afl.retrieve_alphafold_structure("P00520", file_format="invalid")

    def test_download_structures_by_accessions(self):
        """Test bulk download of multiple structures."""
        accessions = ["P00520", "P12345"]
        self.afl.download_structures_by_accessions(
            accessions, file_format="pdb", max_num_threads=2
        )

        # Check that files were created
        downloaded_files = []
        for root, dirs, files in os.walk(self.temp_dir):
            downloaded_files.extend([f for f in files if f.endswith(".pdb")])

        # Should have at least one file (some might fail, that's ok)
        self.assertGreater(len(downloaded_files), 0)


class AlphaFoldListOfflineTests(unittest.TestCase):
    """Tests that don't require internet."""

    def test_initialization(self):
        """Test AlphaFoldList initialization."""
        afl = AlphaFoldList()
        self.assertEqual(afl.api_server, "https://alphafold.ebi.ac.uk")
        self.assertEqual(
            afl.ftp_server, "https://ftp.ebi.ac.uk/pub/databases/alphafold"
        )
        self.assertTrue(afl._verbose)
        self.assertFalse(afl.flat_tree)

    def test_initialization_custom(self):
        """Test AlphaFoldList initialization with custom parameters."""
        custom_dir = "/custom/path"
        afl = AlphaFoldList(local_dir=custom_dir, verbose=False)
        self.assertEqual(afl.local_dir, custom_dir)
        self.assertFalse(afl._verbose)


if __name__ == "__main__":
    unittest.main()
