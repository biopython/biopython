# Copyright 2009-2011 by Eric Talevich.  All rights reserved.
# Revisions copyright 2009-2013 by Peter Cock.  All rights reserved.
# Revisions copyright 2013 Lenna X. Peterson. All rights reserved.
# Revisions copyright 2020 Joao Rodrigues. All rights reserved.
#
# Converted by Eric Talevich from an older unit test copyright 2002
# by Thomas Hamelryck.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Unit tests for the Bio.PDB.PDBIO module."""

import os
import tempfile
import unittest
import warnings

from Bio import BiopythonWarning
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB import Atom, Residue
from Bio.PDB.PDBExceptions import (
    PDBConstructionException,
    PDBConstructionWarning,
    PDBIOException,
)


class WriteTest(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.io = PDBIO()
        self.parser = PDBParser(PERMISSIVE=1)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            self.structure = self.parser.get_structure("example", "PDB/1A8O.pdb")

    def test_pdbio_write_structure(self):
        """Write a full structure using PDBIO."""
        struct1 = self.structure
        # Ensure that set_structure doesn't alter parent
        parent = struct1.parent

        # Write full model to temp file
        self.io.set_structure(struct1)
        self.assertIs(parent, struct1.parent)

        filenumber, filename = tempfile.mkstemp()
        os.close(filenumber)

        try:
            self.io.save(filename)

            struct2 = self.parser.get_structure("1a8o", filename)
            nresidues = len(list(struct2.get_residues()))

            self.assertEqual(len(struct2), 1)
            self.assertEqual(nresidues, 158)
        finally:
            os.remove(filename)

    def test_pdbio_write_preserve_numbering(self):
        """Test writing PDB and preserve atom numbering."""
        self.io.set_structure(self.structure)

        filenumber, filename = tempfile.mkstemp()
        os.close(filenumber)

        try:
            self.io.save(filename)  # default preserve_atom_numbering=False

            struct = self.parser.get_structure("1a8o", filename)
            serials = [a.serial_number for a in struct.get_atoms()]
            og_serials = list(range(1, len(serials) + 1))
            self.assertEqual(og_serials, serials)
        finally:
            os.remove(filename)

    def test_pdbio_pdb_format_limits(self):
        """Test raising error when structure cannot meet PDB format limits."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            structure = self.parser.get_structure("example", "PDB/1A8O.pdb")

        # Modify structure and check if parser raises an error
        # Chain id
        structure[0]["A"].id = "AA"
        self.io.set_structure(structure)
        filenumber, filename = tempfile.mkstemp()
        os.close(filenumber)
        with self.assertRaises(PDBIOException):
            self.io.save(filename)
        structure[0]["AA"].id = "A"
        os.remove(filename)

        # Residue id
        het, ori, ins = structure[0]["A"][152].id
        structure[0]["A"][152].id = (het, 10000, ins)
        self.io.set_structure(structure)
        filenumber, filename = tempfile.mkstemp()
        os.close(filenumber)
        with self.assertRaises(PDBIOException):
            self.io.save(filename)
        structure[0]["A"][10000].id = (het, ori, ins)
        os.remove(filename)

        # Atom id
        structure[0]["A"][152]["CA"].serial_number = 1e6
        self.io.set_structure(structure)
        filenumber, filename = tempfile.mkstemp()
        os.close(filenumber)
        with self.assertRaises(PDBIOException):
            # perserve_... must be True for exception to trigger
            self.io.save(filename, preserve_atom_numbering=True)
        os.remove(filename)

    def test_pdbio_write_auto_numbering(self):
        """Test writing PDB and do not preserve atom numbering."""
        self.io.set_structure(self.structure)

        filenumber, filename = tempfile.mkstemp()
        os.close(filenumber)

        try:
            self.io.save(filename, preserve_atom_numbering=True)

            struct = self.parser.get_structure("1a8o", filename)
            serials = [a.serial_number for a in struct.get_atoms()]
            og_serials = [a.serial_number for a in self.structure.get_atoms()]

            self.assertEqual(og_serials, serials)
        finally:
            os.remove(filename)

    def test_pdbio_write_residue(self):
        """Write a single residue using PDBIO."""
        struct1 = self.structure
        residue1 = list(struct1.get_residues())[0]

        # Ensure that set_structure doesn't alter parent
        parent = residue1.parent

        # Write full model to temp file
        self.io.set_structure(residue1)
        self.assertIs(parent, residue1.parent)
        filenumber, filename = tempfile.mkstemp()
        os.close(filenumber)
        try:
            self.io.save(filename)
            struct2 = self.parser.get_structure("1a8o", filename)
            nresidues = len(list(struct2.get_residues()))
            self.assertEqual(nresidues, 1)
        finally:
            os.remove(filename)

    def test_pdbio_write_residue_w_chain(self):
        """Write a single residue (chain id == X) using PDBIO."""
        struct1 = self.structure.copy()  # make copy so we can change it
        residue1 = list(struct1.get_residues())[0]

        # Modify parent id
        parent = residue1.parent
        parent.id = "X"

        # Write full model to temp file
        self.io.set_structure(residue1)
        filenumber, filename = tempfile.mkstemp()
        os.close(filenumber)
        try:
            self.io.save(filename)
            struct2 = self.parser.get_structure("1a8o", filename)
            nresidues = len(list(struct2.get_residues()))
            self.assertEqual(nresidues, 1)

            # Assert chain remained the same
            chain_id = [c.id for c in struct2.get_chains()][0]
            self.assertEqual(chain_id, "X")
        finally:
            os.remove(filename)

    def test_pdbio_write_residue_wout_chain(self):
        """Write a single orphan residue using PDBIO."""
        struct1 = self.structure
        residue1 = list(struct1.get_residues())[0]

        residue1.parent = None  # detach residue

        # Write full model to temp file
        self.io.set_structure(residue1)

        filenumber, filename = tempfile.mkstemp()
        os.close(filenumber)
        try:
            self.io.save(filename)
            struct2 = self.parser.get_structure("1a8o", filename)
            nresidues = len(list(struct2.get_residues()))
            self.assertEqual(nresidues, 1)

            # Assert chain is default: "A"
            chain_id = [c.id for c in struct2.get_chains()][0]
            self.assertEqual(chain_id, "A")
        finally:
            os.remove(filename)

    def test_pdbio_write_custom_residue(self):
        """Write a chainless residue using PDBIO."""
        res = Residue.Residue((" ", 1, " "), "DUM", "")
        atm = Atom.Atom("CA", [0.1, 0.1, 0.1], 1.0, 1.0, " ", "CA", 1, "C")
        res.add(atm)

        # Ensure that set_structure doesn't alter parent
        parent = res.parent

        # Write full model to temp file
        self.io.set_structure(res)

        self.assertIs(parent, res.parent)
        filenumber, filename = tempfile.mkstemp()
        os.close(filenumber)
        try:
            self.io.save(filename)
            struct2 = self.parser.get_structure("res", filename)
            latoms = list(struct2.get_atoms())
            self.assertEqual(len(latoms), 1)
            self.assertEqual(latoms[0].name, "CA")
            self.assertEqual(latoms[0].parent.resname, "DUM")
            self.assertEqual(latoms[0].parent.parent.id, "A")
        finally:
            os.remove(filename)

    def test_pdbio_select(self):
        """Write a selection of the structure using a Select subclass."""
        # Selection class to filter all alpha carbons
        class CAonly(Select):
            """Accepts only CA residues."""

            def accept_atom(self, atom):
                if atom.name == "CA" and atom.element == "C":
                    return 1

        struct1 = self.structure
        # Ensure that set_structure doesn't alter parent
        parent = struct1.parent
        # Write to temp file
        self.io.set_structure(struct1)

        self.assertIs(parent, struct1.parent)
        filenumber, filename = tempfile.mkstemp()
        os.close(filenumber)
        try:
            self.io.save(filename, CAonly())
            struct2 = self.parser.get_structure("1a8o", filename)
            nresidues = len(list(struct2.get_residues()))
            self.assertEqual(nresidues, 70)
        finally:
            os.remove(filename)

    def test_pdbio_missing_occupancy(self):
        """Write PDB file with missing occupancy."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            structure = self.parser.get_structure("test", "PDB/occupancy.pdb")

        self.io.set_structure(structure)
        filenumber, filename = tempfile.mkstemp()
        os.close(filenumber)
        try:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always", BiopythonWarning)
                self.io.save(filename)
                self.assertEqual(len(w), 1, w)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", PDBConstructionWarning)
                struct2 = self.parser.get_structure("test", filename)
            atoms = struct2[0]["A"][(" ", 152, " ")]
            self.assertIsNone(atoms["N"].get_occupancy())
        finally:
            os.remove(filename)

    def test_pdbio_write_truncated(self):
        """Test parsing of truncated lines."""
        struct = self.structure

        # Write to temp file
        self.io.set_structure(struct)
        filenumber, filename = tempfile.mkstemp()
        os.close(filenumber)
        try:
            self.io.save(filename)
            # Check if there are lines besides 'ATOM', 'TER' and 'END'
            with open(filename) as handle:
                record_set = {line[0:6] for line in handle}
            record_set -= {
                "ATOM  ",
                "HETATM",
                "MODEL ",
                "ENDMDL",
                "TER\n",
                "TER   ",
                "END\n",
                "END   ",
            }
            self.assertEqual(len(record_set), 0)
        finally:
            os.remove(filename)

    def test_model_numbering(self):
        """Preserve model serial numbers during I/O."""

        def confirm_numbering(struct):
            self.assertEqual(len(struct), 3)
            for idx, model in enumerate(struct):
                self.assertEqual(model.serial_num, idx + 1)
                self.assertEqual(model.serial_num, model.id + 1)

        def confirm_single_end(fname):
            """Ensure there is only one END statement in multi-model files."""
            with open(fname) as handle:
                end_stment = []
                for iline, line in enumerate(handle):
                    if line.strip() == "END":
                        end_stment.append((line, iline))
            self.assertEqual(len(end_stment), 1)  # Only one?
            self.assertEqual(end_stment[0][1], iline)  # Last line of the file?

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            struct1 = self.parser.get_structure("1lcd", "PDB/1LCD.pdb")

        confirm_numbering(struct1)

        # Round trip: serialize and parse again
        self.io.set_structure(struct1)
        filenumber, filename = tempfile.mkstemp()
        os.close(filenumber)
        try:
            self.io.save(filename)
            struct2 = self.parser.get_structure("1lcd", filename)
            confirm_numbering(struct2)
            confirm_single_end(filename)
        finally:
            os.remove(filename)

    def test_pdbio_write_x_element(self):
        """Write a structure with atomic element X with PDBIO."""
        struct1 = self.structure

        # Change element of one atom
        atom = next(struct1.get_atoms())
        atom.element = "X"  # X is assigned in Atom.py as last resort

        self.io.set_structure(struct1)

        filenumber, filename = tempfile.mkstemp()
        os.close(filenumber)

        try:
            self.io.save(filename)
        finally:
            os.remove(filename)

    def test_pdbio_write_unk_element(self):
        """PDBIO raises PDBIOException when writing unrecognised atomic elements."""
        struct1 = self.structure

        atom = next(struct1.get_atoms())
        atom.element = "1"

        self.io.set_structure(struct1)

        filenumber, filename = tempfile.mkstemp()
        os.close(filenumber)

        with self.assertRaises(PDBIOException):
            self.io.save(filename)
        os.remove(filename)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
