import unittest
from Bio.PDB.MMTF import MMTFParser
from Bio.PDB.MMCIFParser import MMCIFParser

class ParseMMTF(unittest.TestCase):
    """Testing with real MMTF file(s)."""

    def test_parser(self):
        parser = MMTFParser()
        structure = parser.get_structure("PDB/4CUP.mmtf")

    def test_compare_to_mmcif(self):
        mmtf_parser = MMTFParser()
        mmtf_struct = mmtf_parser.get_structure("PDB/4CUP.mmtf")
        mmCIF_parser = MMCIFParser()
        mmCIF_struct = mmCIF_parser.get_structure("example", "PDB/4CUP.cif")
        self.mmCIF_atoms = [x for x in mmCIF_struct.get_atoms()]
        self.mmtf_atoms = [x for x in mmtf_struct.get_atoms()]
        self.check_atoms()
        mmCIF_chains = [x for x in mmCIF_struct.get_chains()]
        mmtf_chains = [x for x in mmtf_struct.get_chains()]
        self.assertEqual(len(mmCIF_chains), len(mmtf_chains))
        for i in range(len(mmCIF_chains)):
            self.mmCIF_res = [x for x in mmCIF_chains[i].get_residues()]
            self.mmtf_res = [x for x in mmtf_chains[i].get_residues()]
            self.check_residues()

        self.mmCIF_res = [x for x in mmCIF_struct.get_residues()]
        self.mmtf_res = [x for x in mmtf_struct.get_residues()]
        self.check_residues()


        self.assertEqual(len([x for x in mmCIF_struct.get_models()]), len([x for x in mmtf_struct.get_models()]))

    def check_residues(self):
        self.assertEqual(len(self.mmCIF_res), len(self.mmtf_res))
        for i in range(len(self.mmCIF_res)):
            mmCIF_r = self.mmCIF_res[i]
            mmtf_r=self.mmtf_res[i]
            self.assertEqual(mmtf_r.level,mmCIF_r.level)
            self.assertEqual(mmtf_r.disordered,mmCIF_r.disordered)
            self.assertEqual(mmtf_r.resname,mmCIF_r.resname)
            self.assertEqual(mmtf_r.segid,mmCIF_r.segid)
            self.mmCIF_atoms = [x for x in mmCIF_r.get_atom()]
            self.mmtf_atoms = [x for x in mmtf_r.get_atom()]
            self.check_atoms()

    def check_atoms(self):
        self.assertEqual(len(self.mmCIF_atoms),len(self.mmtf_atoms))
        for i in range(len(self.mmCIF_atoms)):
            mmtf_atom = self.mmtf_atoms[i]
            mmCIF_atom = self.mmCIF_atoms[i]
            self.assertEqual(mmtf_atom.name, mmCIF_atom.name)  # eg. CA, spaces are removed from atom name
            self.assertEqual(mmtf_atom.fullname, mmCIF_atom.fullname)  # e.g. " CA ", spaces included
            self.assertAlmostEqual(mmtf_atom.coord[0], mmCIF_atom.coord[0], places=3)
            self.assertAlmostEqual(mmtf_atom.coord[1], mmCIF_atom.coord[1], places=3)
            self.assertAlmostEqual(mmtf_atom.coord[2], mmCIF_atom.coord[2], places=3)
            self.assertEqual(mmtf_atom.bfactor, mmCIF_atom.bfactor)
            self.assertEqual(mmtf_atom.occupancy, mmCIF_atom.occupancy)
            self.assertEqual(mmtf_atom.altloc, mmCIF_atom.altloc)
            self.assertEqual(mmtf_atom.full_id,
                             mmCIF_atom.full_id)  # (structure id, model id, chain id, residue id, atom id)
            self.assertEqual(mmtf_atom.id, mmCIF_atom.name)  # id of atom is the atom name (e.g. "CA")
            # self.assertEqual(mmtf_atom.serial_number,mmCIF_atom.serial_number) # mmCIF serial number is none