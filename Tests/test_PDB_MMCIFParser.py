# Copyright 2012 Lenna X. Peterson (arklenna@gmail.com).
# All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
#
# Tests adapted from test_PDB.py

"""Unit tests for the MMCIF portion of the Bio.PDB module."""

import tempfile
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
        "Install NumPy if you want to use Bio.PDB."
    ) from None


from Bio.Seq import Seq
from Bio.PDB.PDBExceptions import (
    PDBConstructionException,
    PDBConstructionWarning,
    PDBIOException,
)

from Bio.PDB import PPBuilder, CaPPBuilder
from Bio.PDB.MMCIFParser import MMCIFParser, FastMMCIFParser
from Bio.PDB import PDBParser, PDBIO


class ParseReal(unittest.TestCase):
    """Testing with real CIF file(s)."""

    def test_parsers(self):
        """Extract polypeptides from 1A80."""
        parser = MMCIFParser()
        fast_parser = FastMMCIFParser()

        structure = parser.get_structure("example", "PDB/1A8O.cif")
        f_structure = fast_parser.get_structure("example", "PDB/1A8O.cif")

        self.assertEqual(len(structure), 1)
        self.assertEqual(len(f_structure), 1)

        for ppbuild in [PPBuilder(), CaPPBuilder()]:
            # ==========================================================
            # Check that serial_num (model column) is stored properly
            self.assertEqual(structure[0].serial_num, 1)
            self.assertEqual(f_structure[0].serial_num, structure[0].serial_num)

            # First try allowing non-standard amino acids,
            polypeptides = ppbuild.build_peptides(structure[0], False)
            f_polypeptides = ppbuild.build_peptides(f_structure[0], False)

            self.assertEqual(len(polypeptides), 1)
            self.assertEqual(len(f_polypeptides), 1)

            pp = polypeptides[0]
            f_pp = f_polypeptides[0]

            # Check the start and end positions
            self.assertEqual(pp[0].get_id()[1], 151)
            self.assertEqual(pp[-1].get_id()[1], 220)

            self.assertEqual(f_pp[0].get_id()[1], 151)
            self.assertEqual(f_pp[-1].get_id()[1], 220)

            # Check the sequence
            s = pp.get_sequence()
            f_s = f_pp.get_sequence()

            self.assertEqual(s, f_s)  # enough to test this

            self.assertIsInstance(s, Seq)

            # Here non-standard MSE are shown as M
            self.assertEqual(
                "MDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQ"
                "NANPDCKTILKALGPGATLEEMMTACQG",
                s,
            )

            # ==========================================================
            # Now try strict version with only standard amino acids
            # Should ignore MSE 151 at start, and then break the chain
            # at MSE 185, and MSE 214,215
            polypeptides = ppbuild.build_peptides(structure[0], True)
            self.assertEqual(len(polypeptides), 3)

            # First fragment
            pp = polypeptides[0]
            self.assertEqual(pp[0].get_id()[1], 152)
            self.assertEqual(pp[-1].get_id()[1], 184)
            s = pp.get_sequence()
            self.assertIsInstance(s, Seq)
            self.assertEqual("DIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNW", s)

            # Second fragment
            pp = polypeptides[1]
            self.assertEqual(pp[0].get_id()[1], 186)
            self.assertEqual(pp[-1].get_id()[1], 213)
            s = pp.get_sequence()
            self.assertIsInstance(s, Seq)
            self.assertEqual("TETLLVQNANPDCKTILKALGPGATLEE", s)

            # Third fragment
            pp = polypeptides[2]
            self.assertEqual(pp[0].get_id()[1], 216)
            self.assertEqual(pp[-1].get_id()[1], 220)
            s = pp.get_sequence()
            self.assertIsInstance(s, Seq)
            self.assertEqual("TACQG", s)

        s_atoms = list(structure.get_atoms())
        f_atoms = list(f_structure.get_atoms())

        for atoms in [s_atoms, f_atoms]:
            self.assertEqual(len(atoms), 644)
            atom_names = ["N", "CA", "C", "O", "CB"]
            self.assertEqual([a.get_name() for a in atoms[:5]], atom_names)
            self.assertEqual([a.get_id() for a in atoms[:5]], atom_names)
            self.assertEqual([a.get_fullname() for a in atoms[:5]], atom_names)
            self.assertEqual(
                [a.get_occupancy() for a in atoms[:5]], [1.0, 1.0, 1.0, 1.0, 1.0]
            )
            self.assertIsInstance(atoms[0].get_coord(), numpy.ndarray)
            coord = numpy.array([19.594, 32.367, 28.012], dtype=numpy.float32)
            numpy.testing.assert_array_equal(atoms[0].get_coord(), coord)

            self.assertEqual(atoms[0].get_bfactor(), 18.03)
            for atom in atoms:
                self.assertIsNone(atom.get_anisou())

    def test_with_anisotrop(self):
        parser = MMCIFParser()
        fast_parser = FastMMCIFParser()

        structure = parser.get_structure("example", "PDB/4CUP.cif")
        f_structure = fast_parser.get_structure("example", "PDB/4CUP.cif")

        self.assertEqual(len(structure), 1)
        self.assertEqual(len(f_structure), 1)

        s_atoms = list(structure.get_atoms())
        f_atoms = list(f_structure.get_atoms())

        self.assertEqual(len(s_atoms), len(f_atoms))

        for atoms in [s_atoms, f_atoms]:
            atom_names = ["N", "CA", "C", "O", "CB"]
            self.assertEqual([a.get_name() for a in atoms[:5]], atom_names)
            self.assertEqual([a.get_id() for a in atoms[:5]], atom_names)
            self.assertEqual([a.get_fullname() for a in atoms[:5]], atom_names)
            self.assertEqual(
                [a.get_occupancy() for a in atoms[:5]], [1.0, 1.0, 1.0, 1.0, 1.0]
            )
            self.assertIsInstance(atoms[0].get_coord(), numpy.ndarray)
            coord = numpy.array([50.346, 19.287, 17.288], dtype=numpy.float32)
            numpy.testing.assert_array_equal(atoms[0].get_coord(), coord)
            self.assertEqual(atoms[0].get_bfactor(), 32.02)

            ansiou = numpy.array(
                [0.4738, -0.0309, -0.0231, 0.4524, 0.0036, 0.2904], dtype=numpy.float32
            )
            numpy.testing.assert_array_equal(atoms[0].get_anisou(), ansiou)
            ansiou = numpy.array(
                [1.1242, 0.2942, -0.0995, 1.1240, -0.1088, 0.8221], dtype=numpy.float32
            )
            atom_937 = list(f_structure[0]["A"])[114]["CB"]
            numpy.testing.assert_array_equal(atom_937.get_anisou(), ansiou)

    def testModels(self):
        """Test file with multiple models."""
        parser = MMCIFParser(QUIET=1)
        f_parser = FastMMCIFParser(QUIET=1)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            structure = parser.get_structure("example", "PDB/1LCD.cif")
            f_structure = f_parser.get_structure("example", "PDB/1LCD.cif")

        self.assertEqual(len(structure), 3)
        self.assertEqual(len(f_structure), 3)

        for ppbuild in [PPBuilder(), CaPPBuilder()]:
            # ==========================================================
            # Check that serial_num (model column) is stored properly
            self.assertEqual(structure[0].serial_num, 1)
            self.assertEqual(structure[1].serial_num, 2)
            self.assertEqual(structure[2].serial_num, 3)
            # First try allowing non-standard amino acids,
            polypeptides = ppbuild.build_peptides(structure[0], False)
            self.assertEqual(len(polypeptides), 1)
            pp = polypeptides[0]
            # Check the start and end positions
            self.assertEqual(pp[0].get_id()[1], 1)
            self.assertEqual(pp[-1].get_id()[1], 51)
            # Check the sequence
            s = pp.get_sequence()
            self.assertIsInstance(s, Seq)
            # Here non-standard MSE are shown as M
            self.assertEqual("MKPVTLYDVAEYAGVSYQTVSRVVNQASHVSAKTREKVEAAMAELNYIPNR", s)
            # ==========================================================
            # Now try strict version with only standard amino acids
            polypeptides = ppbuild.build_peptides(structure[0], True)
            self.assertEqual(len(polypeptides), 1)
            pp = polypeptides[0]
            # Check the start and end positions
            self.assertEqual(pp[0].get_id()[1], 1)
            self.assertEqual(pp[-1].get_id()[1], 51)
            # Check the sequence
            s = pp.get_sequence()
            self.assertIsInstance(s, Seq)
            self.assertEqual("MKPVTLYDVAEYAGVSYQTVSRVVNQASHVSAKTREKVEAAMAELNYIPNR", s)

        # This structure contains several models with multiple lengths.
        # The tests were failing.
        structure = parser.get_structure("example", "PDB/2OFG.cif")
        self.assertEqual(len(structure), 3)

    def test_insertions(self):
        """Test file with residue insertion codes."""
        parser = MMCIFParser(QUIET=1)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            structure = parser.get_structure("example", "PDB/4ZHL.cif")
        for ppbuild in [PPBuilder(), CaPPBuilder()]:
            # First try allowing non-standard amino acids,
            polypeptides = ppbuild.build_peptides(structure[0], False)
            self.assertEqual(len(polypeptides), 2)
            pp = polypeptides[0]
            # Check the start and end positions (first segment only)
            self.assertEqual(pp[0].get_id()[1], 16)
            self.assertEqual(pp[-1].get_id()[1], 244)
            # Check the sequence
            refseq = (
                "IIGGEFTTIENQPWFAAIYRRHRGGSVTYVCGGSLISPCWVISATHCFIDYPKKEDYIVYLGR"
                "SRLNSNTQGEMKFEVENLILHKDYSADTLAYHNDIALLKIRSKEGRCAQPSRTIQTIALPSMY"
                "NDPQFGTSCEITGFGKEQSTDYLYPEQLKMTVVKLISHRECQQPHYYGSEVTTKMLCAADPQW"
                "KTDSCQGDSGGPLVCSLQGRMTLTGIVSWGRGCALKDKPGVYTRVSHFLPWIRSHTKE"
            )
            s = pp.get_sequence()
            self.assertIsInstance(s, Seq)
            self.assertEqual(refseq, s)

    def test_filehandle(self):
        """Test if the parser can handle file handle as well as filename."""
        parser = MMCIFParser()
        structure = parser.get_structure("example", "PDB/1A8O.cif")
        self.assertEqual(len(structure), 1)

        with open("PDB/1A8O.cif") as handle:
            structure = parser.get_structure("example", handle)
        self.assertEqual(len(structure), 1)

    def test_point_mutations_main(self):
        """Test if MMCIFParser parse point mutations correctly."""
        self._run_point_mutation_tests(MMCIFParser(QUIET=True))

    def test_point_mutations_fast(self):
        """Test if FastMMCIFParser can parse point mutations correctly."""
        self._run_point_mutation_tests(FastMMCIFParser(QUIET=True))

    def _run_point_mutation_tests(self, parser):
        """Shared test code for testing point mutations."""
        structure = parser.get_structure("example", "PDB/3JQH.cif")

        # Residue 1 and 15 should be disordered.
        res_1 = structure[0]["A"][1]
        res_15 = structure[0]["A"][15]

        # Cursory check -- this would be true even if the residue just
        # contained some disordered atoms.
        self.assertTrue(res_1.is_disordered(), "Residue 1 is disordered")
        self.assertTrue(res_15.is_disordered(), "Residue 15 is disordered")

        # Check a non-mutated residue just to be sure we didn't break the
        # parser and cause everything to be disordered.
        self.assertFalse(
            structure[0]["A"][13].is_disordered(), "Residue 13 is not disordered"
        )

        # Check that the residue types were parsed correctly.
        self.assertEqual(
            set(res_1.disordered_get_id_list()),
            {"PRO", "SER"},
            "Residue 1 is proline/serine",
        )
        self.assertEqual(
            set(res_15.disordered_get_id_list()),
            {"ARG", "GLN", "GLU"},
            "Residue 15 is arginine/glutamine/glutamic acid",
        )

        # Quickly check that we can switch between residues and that the
        # correct set of residues was parsed.
        res_1.disordered_select("PRO")
        self.assertAlmostEqual(
            res_1["CA"].get_occupancy(), 0.83, 2, "Residue 1 proline occupancy correcy"
        )

        res_1.disordered_select("SER")
        self.assertAlmostEqual(
            res_1["CA"].get_occupancy(), 0.17, 2, "Residue 1 serine occupancy correcy"
        )

    def test_header(self):
        """Test if the parser populates header data."""
        parser = MMCIFParser(QUIET=1)

        # test default values
        structure = parser.get_structure("example", "PDB/a_structure.cif")
        self.assertEqual("", structure.header["idcode"])
        self.assertEqual("", structure.header["head"])
        self.assertEqual("", structure.header["deposition_date"])
        self.assertEqual("", structure.header["structure_method"])
        self.assertIsNone(structure.header["resolution"])

        # test extracting fields
        structure = parser.get_structure("example", "PDB/1A8O.cif")
        self.assertEqual("1A8O", structure.header["idcode"])
        self.assertEqual("Viral protein", structure.header["head"])
        self.assertEqual("", structure.header["deposition_date"])
        self.assertEqual("X-RAY DIFFRACTION", structure.header["structure_method"])
        self.assertEqual(1.7, structure.header["resolution"])

        # test not confused by '.'
        structure = parser.get_structure("example", "PDB/1SSU_mod.cif")
        self.assertIsNone(structure.header["resolution"])


class CIFtoPDB(unittest.TestCase):
    """Testing conversion between formats: CIF to PDB."""

    def test_conversion(self):
        """Parse 1LCD.cif, write 1LCD.pdb, parse again and compare."""
        cif_parser = MMCIFParser(QUIET=1)
        cif_struct = cif_parser.get_structure("example", "PDB/1LCD.cif")

        pdb_writer = PDBIO()
        pdb_writer.set_structure(cif_struct)
        filenumber, filename = tempfile.mkstemp()
        pdb_writer.save(filename)

        pdb_parser = PDBParser(QUIET=1)
        pdb_struct = pdb_parser.get_structure("example_pdb", filename)

        # comparisons
        self.assertEqual(len(pdb_struct), len(cif_struct))

        pdb_atom_names = [a.name for a in pdb_struct.get_atoms()]
        cif_atom_names = [a.name for a in cif_struct.get_atoms()]
        self.assertEqual(pdb_atom_names, cif_atom_names)

        pdb_atom_elems = [a.element for a in pdb_struct.get_atoms()]
        cif_atom_elems = [a.element for a in cif_struct.get_atoms()]
        self.assertEqual(pdb_atom_elems, cif_atom_elems)

    def test_conversion_not_preserve_numbering(self):
        """Convert mmCIF to PDB and renumber atom serials."""
        cif_parser = MMCIFParser(QUIET=1)
        cif_struct = cif_parser.get_structure("example", "PDB/a_structure.cif")

        pdb_writer = PDBIO()
        pdb_writer.set_structure(cif_struct)
        filenumber, filename = tempfile.mkstemp()

        pdb_writer.save(filename, preserve_atom_numbering=False)

    def test_conversion_preserve_numbering(self):
        """Convert mmCIF to PDB and preserve original serial numbering."""
        cif_parser = MMCIFParser(QUIET=1)
        cif_struct = cif_parser.get_structure("example", "PDB/a_structure.cif")

        pdb_writer = PDBIO()
        pdb_writer.set_structure(cif_struct)
        filenumber, filename = tempfile.mkstemp()

        with self.assertRaises(PDBIOException):
            pdb_writer.save(filename, preserve_atom_numbering=True)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
