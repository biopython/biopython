# Copyright 2016 by Anthony Bradley.  All rights reserved.
# Revisions copyright 2017 by Peter Cock.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

import unittest
import warnings
from Bio.PDB.mmtf import MMTFParser
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning


class ParseMMTF(unittest.TestCase):
    """Testing with real mmtf file(s)."""

    def check_atoms(self):
        """Check all atoms in self.mmtf_atoms and self.mmcif_atoms are equivalent"""
        self.assertEqual(len(self.mmcif_atoms), len(self.mmtf_atoms))
        for i, e in enumerate(self.mmcif_atoms):
            mmtf_atom = self.mmtf_atoms[i]
            mmcif_atom = self.mmcif_atoms[i]
            self.assertEqual(mmtf_atom.name, mmcif_atom.name)  # eg. CA, spaces are removed from atom name
            self.assertEqual(mmtf_atom.fullname, mmcif_atom.fullname)  # e.g. " CA ", spaces included
            self.assertAlmostEqual(mmtf_atom.coord[0], mmcif_atom.coord[0], places=3)
            self.assertAlmostEqual(mmtf_atom.coord[1], mmcif_atom.coord[1], places=3)
            self.assertAlmostEqual(mmtf_atom.coord[2], mmcif_atom.coord[2], places=3)
            self.assertEqual(mmtf_atom.bfactor, mmcif_atom.bfactor)
            self.assertEqual(mmtf_atom.occupancy, mmcif_atom.occupancy)
            self.assertEqual(mmtf_atom.altloc, mmcif_atom.altloc)
            self.assertEqual(mmtf_atom.full_id,
                             mmcif_atom.full_id)  # (structure id, model id, chain id, residue id, atom id)
            self.assertEqual(mmtf_atom.id, mmcif_atom.name)  # id of atom is the atom name (e.g. "CA")
            # self.assertEqual(mmtf_atom.serial_number,mmcif_atom.serial_number) # mmCIF serial number is none
            self.assertEqual(mmtf_atom - mmtf_atom, 0)
            self.assertEqual(mmtf_atom - mmcif_atom, 0)

    def check_residues(self):
        """Check all residues in self.mmcif_res and self.mmtf_res are equivalent"""
        self.assertEqual(len(self.mmcif_res), len(self.mmtf_res))
        for i, e in enumerate(self.mmcif_res):
            mmcif_r = self.mmcif_res[i]
            mmtf_r = self.mmtf_res[i]
            self.assertEqual(mmtf_r.level, mmcif_r.level)
            self.assertEqual(mmtf_r.disordered, mmcif_r.disordered)
            self.assertEqual(mmtf_r.resname, mmcif_r.resname)
            self.assertEqual(mmtf_r.segid, mmcif_r.segid)
            self.mmcif_atoms = [x for x in mmcif_r.get_atoms()]
            self.mmtf_atoms = [x for x in mmtf_r.get_atoms()]
            self.check_atoms()

    def check_mmtf_vs_cif(self, mmtf_filename, cif_filename):
        """Compare parsed structures for MMTF and CIF files."""
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', PDBConstructionWarning)
            mmtf_struct = MMTFParser.get_structure(mmtf_filename)
        mmcif_parser = MMCIFParser()
        mmcif_struct = mmcif_parser.get_structure("4CUP", cif_filename)
        self.mmcif_atoms = [x for x in mmcif_struct.get_atoms()]
        self.mmtf_atoms = [x for x in mmtf_struct.get_atoms()]
        self.check_atoms()
        mmcif_chains = [x for x in mmcif_struct.get_chains()]
        mmtf_chains = [x for x in mmtf_struct.get_chains()]
        self.assertEqual(len(mmcif_chains), len(mmtf_chains))
        for i, e in enumerate(mmcif_chains):
            self.mmcif_res = [x for x in mmcif_chains[i].get_residues()]
            self.mmtf_res = [x for x in mmtf_chains[i].get_residues()]
            self.check_residues()

        self.mmcif_res = [x for x in mmcif_struct.get_residues()]
        self.mmtf_res = [x for x in mmtf_struct.get_residues()]
        self.check_residues()
        self.assertEqual(len([x for x in mmcif_struct.get_models()]), len([x for x in mmtf_struct.get_models()]))

    def test_4CUP(self):
        """Compare parsing 4CUP.mmtf and 4CUP.cif"""
        self.check_mmtf_vs_cif("PDB/4CUP.mmtf", "PDB/4CUP.cif")

# TODO:
#    def test_1A8O(self):
#        """Compare parsing 1A8O.mmtf and 1A8O.cif"""
#        self.check_mmtf_vs_cif("PDB/1A8O.mmtf", "PDB/1A8O.cif")

# TODO:
#    def test_4ZHL(self):
#        """Compare parsing 4ZHL.mmtf and 4ZHL.cif"""
#        self.check_mmtf_vs_cif("PDB/4ZHL.mmtf", "PDB/4ZHL.cif")


class SimpleParseMMTF(unittest.TestCase):
    """Just parse some real mmtf files."""

    def test_4ZHL(self):
        """Parse 4ZHL.mmtf"""
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', PDBConstructionWarning)
            structure = MMTFParser.get_structure("PDB/4ZHL.mmtf")

    def test_1A80(self):
        """Parse 1A8O.mmtf"""
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', PDBConstructionWarning)
            structure = MMTFParser.get_structure("PDB/1A8O.mmtf")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
