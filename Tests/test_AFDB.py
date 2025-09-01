import os
import shutil
import unittest
from unittest.mock import patch, Mock

# Import the class you want to test
from Bio.PDB.AFDB import AFDBList, BASE_URL


# Mock response for a successful download
class MockResponse:
    def __init__(self, content, status_code):
        self.content = content
        self.status_code = status_code


# The main test class, inheriting from unittest.TestCase
class TestAFDBList(unittest.TestCase):

    def setUp(self):
        """Set up a temporary directory for testing."""
        self.temp_dir = "test_cache"
        self.afdb = AFDBList(out_dir=self.temp_dir)

    def tearDown(self):
        """Clean up the temporary directory after each test."""
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)

    @patch('requests.get')
    def test_retrieve_success(self, mock_get):
        """Test a successful download of a PDB file."""
        # Configure the mock to return a successful response
        mock_get.return_value = MockResponse(
            b"ATOM      1  N   MET A   1      -9.011  15.938  20.408  1.00 97.43           N  ", 200)

        # Call the method
        uniprot_id = "P38398"
        filepath = self.afdb.retrieve(uniprot_id, fmt="pdb")

        # Assertions
        expected_path = os.path.join(self.temp_dir, f"AF-{uniprot_id}-F1-model_v4.pdb")
        self.assertTrue(os.path.exists(filepath))
        self.assertEqual(filepath, expected_path)
        mock_get.assert_called_once_with(f"{BASE_URL}/AF-{uniprot_id}-F1-model_v4.pdb")

    @patch('requests.get')
    def test_retrieve_not_found(self, mock_get):
        """Test retrieval for a non-existent UniProt ID."""
        # Configure the mock to return a 404 Not Found status
        mock_get.return_value = MockResponse(b"", 404)

        # Assert that the correct exception is raised
        with self.assertRaises(FileNotFoundError):
            self.afdb.retrieve("NONEXISTENT")

    def test_retrieve_invalid_format(self):
        """Test retrieval with an invalid format string."""
        with self.assertRaises(ValueError):
            self.afdb.retrieve("P38398", fmt="invalid")

    @patch('requests.get')
    def test_retrieve_cached(self, mock_get):
        """Test that a file is not downloaded if it already exists."""
        # Create a dummy cached file
        uniprot_id = "P38398"
        cached_path = os.path.join(self.temp_dir, f"AF-{uniprot_id}-F1-model_v4.pdb")
        os.makedirs(os.path.dirname(cached_path), exist_ok=True)
        with open(cached_path, "w") as f:
            f.write("DUMMY CONTENT")

        # Call the method
        filepath = self.afdb.retrieve(uniprot_id, fmt="pdb")

        # Assertions
        self.assertEqual(filepath, cached_path)
        # Verify that requests.get was NOT called
        mock_get.assert_not_called()

    @patch('requests.get')
    def test_retrieve_cif_format(self, mock_get):
        """Test successful download in CIF format."""
        # Configure the mock for a successful CIF download
        mock_get.return_value = MockResponse(b"data_AF-Q9Y6F2-F1-model_v4\nloop_\n_atom_site.group_PDB", 200)

        # Call the method
        uniprot_id = "Q9Y6F2"
        filepath = self.afdb.retrieve(uniprot_id, fmt="cif")

        # Assertions
        expected_path = os.path.join(self.temp_dir, f"AF-{uniprot_id}-F1-model_v4.cif")
        self.assertTrue(os.path.exists(filepath))
        self.assertEqual(filepath, expected_path)
        mock_get.assert_called_once_with(f"{BASE_URL}/AF-{uniprot_id}-F1-model_v4.cif")

if __name__ == '__main__':
    unittest.main()
