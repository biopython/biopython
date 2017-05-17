"""Unit tests for PARTS of the parse_pdb_header module of Bio.PDB"""

import unittest
import warnings

try:
    import numpy
    from numpy import dot  # Missing on old PyPy's micronumpy
    del dot
    from numpy.linalg import svd, det  # Missing in PyPy 2.0 numpypy
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.PDB.")


from Bio.PDB.parse_pdb_header import parse_pdb_header, _parse_remark_465


class ParseReal(unittest.TestCase):
    """Testing with real PDB file(s)."""
    def test_parse_pdb_with_remark_465(self):
        """
        Tests that parse_pdb_header now can identify some REMARK 465 entries.
        """
        header = parse_pdb_header("PDB/2XHE.pdb")
        self.assertTrue(header["has_missing_residues"])
        self.assertEqual(len(header["missing_residues"]), 142)
        self.assertIn( {"model":None, "res_name":"GLN", "chain":"B", "ssseq":267, "insertion":None}  , header["missing_residues"])
        header = parse_pdb_header("PDB/1A8O.pdb")
        self.assertFalse(header["has_missing_residues"])
        self.assertEqual(header["missing_residues"], [])
    def test_parse_remark_465(self):
        "A UNIT-test for the private function _parse_remark_465"
        out_d = {"has_missing_residues": False, "missing_residues": []}
        _parse_remark_465("", out_d)
        self.assertFalse(out_d["has_missing_residues"])
        _parse_remark_465("GLU B   276", out_d)
        self.assertTrue(out_d["has_missing_residues"])
        self.assertEqual(len(out_d["missing_residues"]), 1)
        out_d = {"has_missing_residues": False, "missing_residues": []}
        _parse_remark_465("SOME OTHER TEXT THAT CAN'T BE PARSED WITH NUMBERS 1234", out_d)
        self.assertTrue(out_d["has_missing_residues"])
        self.assertEqual(out_d["missing_residues"], [])



if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
