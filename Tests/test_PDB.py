# Copyright 2021 by Simon Duerr.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for PDB module."""

import gzip
import sys
import unittest
import warnings

from pathlib import Path

from Bio import PDB
from Bio import BiopythonWarning

from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.mmtf import MMTFParser
from Bio.PDB.PDBExceptions import PDBConstructionException, PDBConstructionWarning


class PDBRead(unittest.TestCase):
    """Base class for Bio.PDB.read unit tests.

    This contains tests related to helper function PDB.read().
    All tests related to proper parsing are in test_PDB_MMTFParser, test_PDB_MMCIFParser, test_PDB_PDBParser
    """

    def setUp(self):
        self.pdb_parser = PDBParser()
        self.mmcif_parser = MMCIFParser()
        self.mmtf_parser = MMTFParser()
        pass

    def test_PDB_read_pathlike(self):
        """Test read function with Path input."""
        p = Path("PDB/1A8O.pdb")
        structure = PDB.read(p, fmt="pdb")

        pgzip = Path("PDB/1A8O.pdb.gz")
        structure = PDB.read(pgzip, fmt="pdb")

    def test_PDB_read_str(self):
        """Test read function with str input."""
        pgzip_str = "PDB/1A8O.pdb.gz"
        structure = PDB.read(pgzip_str, fmt="pdb")

        p_str = "PDB/1A8O.pdb"
        structure = PDB.read(p_str, fmt="pdb")

    def test_PDB_read_handle(self):
        """Test read function with handle input."""
        with open("PDB/1A8O.pdb", mode="rt") as fp:
            structure = PDB.read(fp, fmt="pdb")

        with gzip.open("PDB/1A8O.pdb.gz", mode="rt") as fp:
            structure = PDB.read(fp, fmt="pdb")

    def test_PDB_read_pdb(self):
        """Test read function with pdb files."""
        structure = PDB.read("PDB/1A8O.pdb", fmt="pdb")

        with open("PDB/1A8O.pdb") as handle:
            structure_handle = PDB.read(handle, fmt="pdb")
        structure_parser = self.pdb_parser.get_structure("1A8O", "PDB/1A8O.pdb")

        # Basic check that structure read from file and from handle are identical
        # more intricate checks are run in the unit tests of the Parser itself
        self.assertEqual(
            len(list(structure.get_residues())),
            len(list(structure_handle.get_residues())),
        )
        # Basic check that structure read directly using PDBParser and PDB.read are identical
        self.assertEqual(
            len(list(structure.get_residues())),
            len(list(structure_parser.get_residues())),
        )

    def test_PDB_read_pdbgz(self):
        """Test read function with gzipped pdb files."""
        structure_gzip = PDB.read("PDB/1A8O.pdb.gz", fmt="pdb")

        with gzip.open("PDB/1A8O.pdb.gz", mode="rt") as handle:
            structure_gzip_handle = PDB.read(handle, fmt="pdb")
        structure_parser = self.pdb_parser.get_structure("1A8O", "PDB/1A8O.pdb")

        # Basic check that structure read from file and from handle are identical
        self.assertEqual(
            len(list(structure_gzip.get_residues())),
            len(list(structure_gzip_handle.get_residues())),
        )
        # Basic check that non-gzipped structure read directly using PDBParser and PDB.read are identical
        self.assertEqual(
            len(list(structure_gzip.get_residues())),
            len(list(structure_parser.get_residues())),
        )

    def test_PDB_read_cif(self):
        """Test read function with cif files."""
        structure = PDB.read("PDB/1A8O.cif", fmt="cif")

        with open("PDB/1A8O.cif") as handle:
            structure_handle = PDB.read(handle, fmt="cif")
        structure_parser = self.mmcif_parser.get_structure("1A8O", "PDB/1A8O.cif")

        self.assertEqual(
            len(list(structure.get_residues())),
            len(list(structure_handle.get_residues())),
        )
        self.assertEqual(
            len(list(structure.get_residues())),
            len(list(structure_parser.get_residues())),
        )

    def test_PDB_read_cifgz(self):
        """Test read function with gzipped mmcif files."""
        structure_gzip = PDB.read("PDB/1A8O.cif.gz", fmt="cif")

        with gzip.open("PDB/1A8O.cif.gz", mode="rt") as handle:
            structure_gzip_handle = PDB.read(handle, fmt="cif")
        structure_parser = self.mmcif_parser.get_structure("1A8O", "PDB/1A8O.cif")

        self.assertEqual(
            len(list(structure_gzip.get_residues())),
            len(list(structure_gzip_handle.get_residues())),
        )
        self.assertEqual(
            len(list(structure_gzip.get_residues())),
            len(list(structure_parser.get_residues())),
        )

    def test_PDB_read_mmtf_handle(self):
        """Test mmtf error when using handle."""
        # MMTF reader does not accept handles
        with self.assertRaises(TypeError):
            with open("PDB/1A8O.mmtf") as handle:
                _ = PDB.read(handle, format="mmtf")

    def test_PDB_read_mmtf(self):
        """Test that read(fmt="mmtf") returns a valid structure."""
        structure_read = PDB.read("PDB/1A8O.mmtf", fmt="mmtf")
        structure_parser = self.mmtf_parser.get_structure("PDB/1A8O.mmtf")
        # basic comparison to see that they are equivalent
        self.assertEqual(
            len(list(structure_read.get_residues())),
            len(list(structure_parser.get_residues())),
        )

    def test_PDB_read_kwargs2Parser(self):
        """Check that kwargs are passed to parser."""
        # this structure can only be parsed in permissive mode as it contains errors
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", PDBConstructionWarning)
            structure = PDB.read("PDB/occupancy.pdb", fmt="pdb", PERMISSIVE=True)
            self.assertEqual(len(w), 3, w)
        # With quiet mode no warnings should be printed
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", PDBConstructionWarning)
            structure = PDB.read(
                "PDB/occupancy.pdb", fmt="pdb", QUIET=True, PERMISSIVE=True
            )
            self.assertEqual(len(w), 0, w)

    def test_PDB_read_format(self):
        """Test format errors."""
        with self.assertRaises(TypeError):
            _ = structure = PDB.read("PDB/1A8O.pdb", fmt=0)

        # uppercase format should not raise an error
        _ = structure = PDB.read("PDB/1A8O.pdb", fmt="PDB")

        # empty format
        with self.assertRaises(ValueError):
            _ = structure = PDB.read("PDB/1A8O.pdb", fmt="")
        # unknown format
        with self.assertRaises(ValueError):
            _ = structure = PDB.read("PDB/1A8O.pdb", fmt="blabla")

    def test_PDB_read_id(self):
        """Test reading structure from id if contained in structure."""
        # ID is read from PDB header
        structure = PDB.read("PDB/1A8O.pdb", fmt="pdb")

        with open("PDB/1A8O.pdb") as handle:
            structure_handle = PDB.read(handle, fmt="pdb")
        structure_parser = self.pdb_parser.get_structure("1A8O", "PDB/1A8O.pdb")

        self.assertEqual(structure.id, "1A8O")
        self.assertEqual(structure_handle.id, "1A8O")
        self.assertEqual(structure_parser.id, "1A8O")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
