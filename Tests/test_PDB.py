# Copyright 2021 by Simon Duerr.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for PDB module."""

import gzip
import sys
import unittest
import warnings

from Bio import PDB

from Bio.PDB import PDBParser, MMCIFParser, MMTFParser


class PDBTestBaseClass(unittest.TestCase):
    """Base class for Bio.PDB unit tests."""

    def setUp(self):
        self.PDBParser = PDBParser()
        self.MMCIFParser = MMCIFParser()
        self.MMTFParser = MMTFParser()
        pass

    def test_PDB_read_pdb_pdbgz(self):
        """Test read function with supported formats.

        Can read: pdb, pdb.gz, cif, cif.gz, mmtf
        """
        structure = PDB.read("PDB/1A8O.pdb", format="pdb")
        structure_gzip = PDB.read("PDB/1A8O.pdb", format="pdb")

        with gzip.open("PDB/1A8O.pdb.gz", mode="rt") as handle:
            structure_gzip_handle = PDB.read(handle, format="pdb")

        with open("PDB/1A8O.pdb") as handle:
            structure_handle = PDB.read(handle, format="pdb")

        structure_parser = self.PDBParser.get_structure("1A8O", "PDB/1A8O.pdb")

        self.assertEqual(structure, structure_gzip)
        self.assertEqual(structure, structure_gzip_handle)
        self.assertEqual(structure, structure_handle)
        self.assertEqual(structure, structure_parser)

    def test_PDB_read_mmtf_handle(self):
        """Test mmtf error when using handle."""
        # MMTF reader does not accept handles
        with self.assertRaises(TypeError):
            with open("PDB/1A8O.mmtf") as handle:
                _ = PDB.read(handle, format="mmtf")

    @unittest.skip("warning not implemented yet")
    def test_PDB_read_mmtf_id(self):
        """Test mmtf error when setting id."""
        # MMTF reader does not accept ids, so a warning should show that it id is ignored
        _ = PDB.read("PDB/1A8O.mmtf", format="mmtf", id="test")

    def test_PDB_read_format(self):
        """Test format errors."""
        with self.assertRaises(TypeError):
            _ = structure = PDB.read("PDB/1A8O.pdb", format=0)

        with self.assertRaises(ValueError):
            _ = structure = PDB.read("PDB/1A8O.pdb", format="PDB")

        with self.assertRaises(ValueError):
            _ = structure = PDB.read("PDB/1A8O.pdb", format="")

        with self.assertRaises(ValueError):
            _ = structure = PDB.read("PDB/1A8O.pdb", format="blabla")

    def test_PDB_read_id(self):
        """Test reading from id."""
        # ID is read from PDB header
        structure = PDB.read("PDB/1A8O.pdb", format="pdb")

        structure_cif = PDB.read("PDB/1A8O.cif", format="cif")

        structure_mmtf = PDB.read("PDB/1A8O.mmtf", format="mmtf")

        with open("PDB/1A8O.pdb") as handle:
            structure_handle = PDB.read(handle, format="pdb")

        structure_parser = self.PDBParser.get_structure("1A8O", "PDB/1A8O.pdb")

        self.assertEqual(structure.id, "1A8O")
        self.assertEqual(structure_handle.id, "1A8O")
        self.assertEqual(structure_cif.id, "1A8O")
        self.assertEqual(structure_parser.id, "1A8O")
        self.assertEqual(structure_mmtf.id, "1A8O")

        # PDB file with no header, id is read from file name
        structure = PDB.read("PDB/ions.pdb", format="pdb")

        self.assertEqual(structure.id, "ions")

        # check if setting id works
        structure = PDB.read("PDB/1A8O.pdb", format="pdb", id="testid")
        self.assertEqual(structure.id, "testid")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
