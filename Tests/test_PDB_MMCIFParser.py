# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Unit tests for the Bio.PDB.MMCIFParser module."""

import unittest

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.PDB.")

from Bio.PDB import MMCIFParser, FastMMCIFParser


class MMCIFParserTests(unittest.TestCase):
    r"""Test module MMCIFParser with CIF files located in TEst\PDB."""

    def test_MMCIFParser(self):
        """Test method MMCIFParser when reading a CIF file."""
        filename = "PDB/1A8O.cif"
        p = MMCIFParser()
        structure = p.get_structure("test", filename)
        res_list = [res.get_resname() for res in structure.get_residues()]
        self.assertEqual(len(res_list), 158)
        self.assertEqual("-".join(res_list), 'MSE-ASP-ILE-ARG-GLN-GLY-PRO-LYS-GLU-PRO-PHE-ARG-ASP-TYR-VAL-ASP-ARG-PHE-TYR-LYS-THR-LEU-ARG-ALA-GLU-GLN-ALA-SER-GLN-GLU-VAL-LYS-ASN-TRP-MSE-THR-GLU-THR-LEU-LEU-VAL-GLN-ASN-ALA-ASN-PRO-ASP-CYS-LYS-THR-ILE-LEU-LYS-ALA-LEU-GLY-PRO-GLY-ALA-THR-LEU-GLU-GLU-MSE-MSE-THR-ALA-CYS-GLN-GLY-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH')

    def test_FastMMCIFParser(self):
        """Test alternative version of the previous method when reading a CIF file."""
        filename = "PDB/1A8O.cif"
        p = FastMMCIFParser()
        structure = p.get_structure("test", filename)
        res_list = [res.get_resname() for res in structure.get_residues()]
        self.assertEqual(len(res_list), 158)
        self.assertEqual("-".join(res_list), 'MSE-ASP-ILE-ARG-GLN-GLY-PRO-LYS-GLU-PRO-PHE-ARG-ASP-TYR-VAL-ASP-ARG-PHE-TYR-LYS-THR-LEU-ARG-ALA-GLU-GLN-ALA-SER-GLN-GLU-VAL-LYS-ASN-TRP-MSE-THR-GLU-THR-LEU-LEU-VAL-GLN-ASN-ALA-ASN-PRO-ASP-CYS-LYS-THR-ILE-LEU-LYS-ALA-LEU-GLY-PRO-GLY-ALA-THR-LEU-GLU-GLU-MSE-MSE-THR-ALA-CYS-GLN-GLY-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH-HOH')


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
