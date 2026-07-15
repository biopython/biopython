"""
Tests for the PDBML parser in the PDB package.

These tests rely on the principle that the structure returned by the PDBML parser should be the same as the structure
returned by the mmCIF parser for any PDB structure.
"""

import io
import unittest
import warnings

from Bio.PDB import MMCIFParser
from Bio.PDB import PDBMLParser
from Bio.PDB.PDBExceptions import PDBConstructionException
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
                # `resolution` was previously always None because of the
                # childless-leaf `if element:` truthiness bug; it must now match
                # the mmCIF parser, i.e. _refine.ls_d_res_high
                # (1GBT: 2.0, 6WG6: 3.54, 3JQH: 2.201 A).
                self.assertAlmostEqual(
                    pdbml_structure.header["resolution"],
                    mmcif_structure.header["resolution"],
                )

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

    def test_missing_required_element_raises(self):
        """A missing required element raises a clear PDBConstructionException."""
        pdbml_parser = PDBMLParser()
        # Well-formed PDBML that declares the PDBx namespace but omits the
        # required structCategory/struct/title header element.
        incomplete_pdbml = (
            '<?xml version="1.0" encoding="UTF-8"?>\n'
            '<PDBx:datablock xmlns:PDBx="http://pdbml.pdb.org/schema/pdbx-v50.xsd">\n'
            "</PDBx:datablock>\n"
        )
        with self.assertRaisesRegex(
            PDBConstructionException, "Required PDBML element not found"
        ):
            pdbml_parser.get_structure(io.StringIO(incomplete_pdbml))

    def test_optional_header_text_may_be_absent(self):
        """A present-but-nil header text element yields None, not an error."""
        pdbml_parser = PDBMLParser()
        # struct_keywords/text is present but xsi:nil, as the wwPDB PDBML
        # generator emits missing scalar values. It should parse to None rather
        # than raising, matching the original parser and the mmCIF parser.
        pdbml = (
            '<?xml version="1.0" encoding="UTF-8"?>\n'
            '<PDBx:datablock xmlns:PDBx="http://pdbml.pdb.org/schema/pdbx-v50.xsd"'
            ' xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">\n'
            "<PDBx:structCategory><PDBx:struct>"
            "<PDBx:title>Test structure</PDBx:title>"
            "</PDBx:struct></PDBx:structCategory>\n"
            "<PDBx:struct_keywordsCategory><PDBx:struct_keywords>"
            '<PDBx:text xsi:nil="true"/>'
            "</PDBx:struct_keywords></PDBx:struct_keywordsCategory>\n"
            '<PDBx:entryCategory><PDBx:entry id="TEST"/></PDBx:entryCategory>\n'
            "<PDBx:pdbx_database_statusCategory><PDBx:pdbx_database_status>"
            "<PDBx:recvd_initial_deposition_date>2020-01-01"
            "</PDBx:recvd_initial_deposition_date>"
            "</PDBx:pdbx_database_status></PDBx:pdbx_database_statusCategory>\n"
            '<PDBx:exptlCategory><PDBx:exptl method="X-RAY DIFFRACTION"/>'
            "</PDBx:exptlCategory>\n"
            "<PDBx:atom_siteCategory/>\n"
            "</PDBx:datablock>\n"
        )
        structure = pdbml_parser.get_structure(io.StringIO(pdbml))
        self.assertEqual(structure.header["name"], "Test structure")
        self.assertIsNone(structure.header["head"])


if __name__ == "__main__":
    unittest.main()
