"""
Tests for the PDBML parser in the PDB package.

These tests rely on the principle that the structure returned by the PDBML parser should be the same as the structure
returned by the mmCIF parser for any PDB structure.
"""

import unittest
import warnings

from Bio.PDB import MMCIFParser
from Bio.PDB import PDBMLParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning


class TestPDBMLParser(unittest.TestCase):
    def test_get_structure(self):
        mmcif_parser = MMCIFParser()
        pdbml_parser = PDBMLParser()

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            for entry in ["1GBT", "6WG6", "3JQH"]:
                mmcif_structure = mmcif_parser.get_structure(entry, f"PDB/{entry}.cif")
                pdbml_structure = pdbml_parser.get_structure(f"PDB/{entry}.xml")
            self.assertEqual(mmcif_structure, pdbml_structure)

    def test_get_structure_filehandle(self):
        mmcif_parser = MMCIFParser()
        pdbml_parser = PDBMLParser()

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            for entry in ["1GBT"]:
                with (
                    open(f"PDB/{entry}.cif") as mmcif_file,
                    open(f"PDB/{entry}.xml") as pdbml_file,
                ):
                    mmcif_structure = mmcif_parser.get_structure(entry, mmcif_file)
                    pdbml_structure = pdbml_parser.get_structure(pdbml_file)
                self.assertEqual(mmcif_structure, pdbml_structure)


if __name__ == "__main__":
    unittest.main()
