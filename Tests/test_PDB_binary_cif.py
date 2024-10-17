"""
Tests for BinaryCIF code in the PDB package.
"""

import unittest

from Bio.PDB import MMCIFParser
from Bio.PDB.binary_cif import BinaryCIFParser


class TestBinaryCIFParser(unittest.TestCase):
    def test_get_structure(self):
        mmcif_parser = MMCIFParser(auth_chains=False)
        bcif_parser = BinaryCIFParser()

        for entry in ["1GBT", "6WG6", "3JQH"]:
            mmcif_structure = mmcif_parser.get_structure(entry, f"PDB/{entry}.cif")
            bcif_structure = bcif_parser.get_structure(
                entry, f"PDB/{entry.lower()}.bcif.gz"
            )
            self.assertTrue(
                mmcif_structure.strictly_equals(
                    bcif_structure, compare_coordinates=True
                )
            )
