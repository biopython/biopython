# Copyright 2019-2020 by Konstantinos Zisis and Artemi Bendandi.  All rights reserved.
#
# Converted by Konstantinos Zisis and Artemi Bendandi from an older unit test
# copyright 2009 by Eric Talevich.
#
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for the PQR parser in Bio.PDB module."""
from Bio.PDB.PDBParser import PDBParser
import unittest
import warnings
from io import StringIO
from Bio.PDB.PDBExceptions import PDBConstructionException, PDBConstructionWarning
from Bio.PDB.PDBIO import PDBIO

import tempfile
import os


class ParseSimplePQR(unittest.TestCase):
    """Parse a simple PQR entry and check behaviour for various problematic inputs."""

    def test_single_input(self):
        """Test if a single ATOM entry correctly parsed."""
        # Atom Entry
        data = "ATOM      1  N   PRO     1      000001  02.000 3.0000 -0.1000  1.0000       N\n"

        # Get sole atom of this structure
        parser = PDBParser(is_pqr=True)  # default initialization
        struct = parser.get_structure("test", StringIO(data))
        atom = next(struct.get_atoms())

        # Check that charge and radius are properly initallized
        self.assertEqual(atom.get_charge(), -0.1)
        self.assertEqual(atom.get_radius(), 1.0)
        self.assertIsNone(atom.get_occupancy())
        self.assertIsNone(atom.get_bfactor())

        # Coordinates
        for i in range(1, 3):
            self.assertEqual(atom.get_coord()[i], i + 1)

    def test_bad_xyz(self):
        """Test if bad coordinates exception is raised."""
        # Atom Entry
        data = "ATOM      1  N   PRO     1      00abc1  02.000 3.0000 -0.1000  1.0000       N\n"

        # Get sole atom of this structure
        parser = PDBParser(is_pqr=True)  # default initialization
        self.assertRaises(
            PDBConstructionException, parser.get_structure, "example", StringIO(data)
        )

    def test_bad_charge(self):
        """Test if missing or malformed charge case is handled correctly."""
        # Test Entries
        malformed = "ATOM      1  N   PRO     1      000001  02.000 3.0000 -0.W000  1.0000       N\n"
        missing = "ATOM      1  N   PRO     1      000001  02.000 3.0000          1.0000       N\n"

        # Malformed
        parser = PDBParser(PERMISSIVE=True, is_pqr=True)  # default initialization
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", PDBConstructionWarning)
            structure = parser.get_structure("test", StringIO(malformed))

        atom = next(structure.get_atoms())
        self.assertIsNone(atom.get_charge())

        # Missing
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", PDBConstructionWarning)
            structure = parser.get_structure("test", StringIO(missing))

        atom = next(structure.get_atoms())
        self.assertIsNone(atom.get_charge())

        # Test PERMISSIVE mode behaviour
        parser = PDBParser(PERMISSIVE=False, is_pqr=True)  # default initialization
        self.assertRaises(
            PDBConstructionException,
            parser.get_structure,
            "example",
            StringIO(malformed),
        )

    def test_bad_radius(self):
        """Test if missing, malformed or negative radius case is handled correctly."""
        # Test Entries
        malformed = "ATOM      1  N   PRO     1      000001  02.000 3.0000 -0.1000  1.a00f       N\n"
        missing = "ATOM      1  N   PRO     1      000001  02.000 3.0000 -0.1000               N\n"
        negative = "ATOM      1  N   PRO     1      000001  02.000 3.0000 -0.1000 -1.0000       N\n"

        # Malformed
        parser = PDBParser(PERMISSIVE=True, is_pqr=True)  # default initialization
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", PDBConstructionWarning)
            structure = parser.get_structure("test", StringIO(malformed))

        atom = next(structure.get_atoms())
        self.assertIsNone(atom.get_radius())

        # Missing
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", PDBConstructionWarning)
            structure = parser.get_structure("test", StringIO(missing))

        atom = next(structure.get_atoms())
        self.assertIsNone(atom.get_radius())

        # Negative
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", PDBConstructionWarning)
            structure = parser.get_structure("test", StringIO(negative))

        atom = next(structure.get_atoms())
        self.assertIsNone(atom.get_radius())

        # Test PERMISSIVE mode behaviour
        parser = PDBParser(PERMISSIVE=False, is_pqr=True)  # default initialization
        self.assertRaises(
            PDBConstructionException,
            parser.get_structure,
            "example",
            StringIO(malformed),
        )
        self.assertRaises(
            PDBConstructionException,
            parser.get_structure,
            "example",
            StringIO(negative),
        )
        self.assertRaises(
            PDBConstructionException, parser.get_structure, "example", StringIO(missing)
        )


class WriteTest(unittest.TestCase):
    """Test if the PDBIO module correctly exports .pqr files."""

    def setUp(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)

            # Open a parser in permissive mode and parse an example file
            self.pqr_parser = PDBParser(PERMISSIVE=1, is_pqr=True)
            self.example_structure = self.pqr_parser.get_structure(
                "example", "PQR/1A80.pqr"
            )

    def test_pdbio_write_pqr_structure(self):
        """Write a full structure using PDBIO."""
        # Create a PDBIO object in pqr mode with example_structure as an argument
        io = PDBIO(is_pqr=True)
        io.set_structure(self.example_structure)

        # Write to a temporary file
        filenumber, filename = tempfile.mkstemp()
        os.close(filenumber)
        try:
            # Export example_structure to a temp file
            io.save(filename)

            # Parse exported structure
            output_struct = self.pqr_parser.get_structure("1a8o", filename)

            # Comparisons
            self.assertEqual(
                len(output_struct), len(self.example_structure)
            )  # Structure Length

            original_residues = len(list(self.example_structure.get_residues()))
            parsed_residues = len(list(output_struct.get_residues()))
            self.assertEqual(parsed_residues, original_residues)  # Number of Residues

            # Atom-wise comparison
            original_atoms = self.example_structure.get_atoms()
            for atom in output_struct.get_atoms():
                self.assertEqual(atom, next(original_atoms))

        finally:
            os.remove(filename)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
