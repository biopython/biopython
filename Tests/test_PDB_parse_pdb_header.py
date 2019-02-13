# Copyright 2017 by Bernhard Thiel.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Unit tests for PARTS of the parse_pdb_header module of Bio.PDB"""


import unittest
import warnings

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.PDB.")

from Bio.PDB.parse_pdb_header import parse_pdb_header, _parse_remark_465


class ParseReal(unittest.TestCase):
    """Testing with real PDB file(s)."""

    def test_parse_pdb_with_remark_465(self):
        """Tests that parse_pdb_header now can identify some REMARK 465 entries."""
        header = parse_pdb_header("PDB/2XHE.pdb")
        self.assertTrue(header["has_missing_residues"])
        self.assertEqual(len(header["missing_residues"]), 142)
        self.assertIn({"model": None, "res_name": "GLN", "chain": "B", "ssseq": 267, "insertion": None},
                      header["missing_residues"])
        header = parse_pdb_header("PDB/1A8O.pdb")
        self.assertFalse(header["has_missing_residues"])
        self.assertEqual(header["missing_residues"], [])

    def test_parse_remark_465(self):
        """A UNIT-test for the private function _parse_remark_465."""
        info = _parse_remark_465("GLU B   276")
        self.assertEqual(info, {"model": None, "res_name": "GLU",
                                "chain": "B", "ssseq": 276, "insertion": None})

        info = _parse_remark_465("2 GLU B   276B")
        self.assertEqual(info, {"model": 2, "res_name": "GLU",
                                "chain": "B", "ssseq": 276, "insertion": "B"})

        info = _parse_remark_465("A 2    11")
        self.assertEqual(info, {"model": None, "res_name": "A",
                                "chain": "2", "ssseq": 11, "insertion": None})


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
