# Copyright 2017 by Bernhard Thiel.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Unit tests for PARTS of the parse_pdb_header module of Bio.PDB."""


import unittest

try:
    import numpy  # noqa F401
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError("Install NumPy if you want to use Bio.PDB.")

from Bio.PDB.ic_rebuild import structure_rebuild_test
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser


class Rebuild(unittest.TestCase):
    """Read PDB and mmCIF structures, convert to/from internal coordinates."""

    PDB_parser = PDBParser(PERMISSIVE=True, QUIET=True)
    CIF_parser = MMCIFParser(QUIET=True)

    def test_rebuild_multichain_missing_residues(self):
        """Convert multichain protein to internal coordinates and back."""
        struct = self.PDB_parser.get_structure("1LCD", "PDB/1LCD.pdb")
        r = structure_rebuild_test(struct, False)
        # print(r)
        self.assertEqual(r["residues"], 153)
        self.assertEqual(r["rCount"], 360)
        self.assertEqual(r["rMatchCount"], 360)
        self.assertEqual(r["aCount"], 1491)
        self.assertEqual(r["disAtmCount"], 0)
        self.assertEqual(r["aCoordMatchCount"], 1491)
        self.assertEqual(len(r["chains"]), 3)
        self.assertTrue(r["pass"])

    def test_rebuild_disordered_atoms_residues(self):
        """Convert disordered protein to internal coordinates and back."""
        struct = self.CIF_parser.get_structure("3JQH", "PDB/3JQH.cif")
        # 3jqh has both disordered residues
        # and disordered atoms in ordered residues
        r = structure_rebuild_test(struct, False)
        print(r)
        self.assertEqual(r["residues"], 26)
        self.assertEqual(r["rCount"], 47)
        self.assertEqual(r["rMatchCount"], 47)
        self.assertEqual(r["aCount"], 217)
        self.assertEqual(r["disAtmCount"], 50)
        self.assertEqual(r["aCoordMatchCount"], 217)
        self.assertEqual(len(r["chains"]), 1)
        self.assertTrue(r["pass"])


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
