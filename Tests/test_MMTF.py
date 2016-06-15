import unittest
import warnings
from Bio.PDB.mmtf import MMTFParser
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning


class ParseMMTF(unittest.TestCase):
    """Testing with real mmtf file(s)."""

    def test_compare_to_mmcif(self):
        """Compre the MMTF and mmCIF parsed structrues"""
        def test_atoms(parse_mmtf):
            """Test that all atoms in self.mmtf_atoms and self.mmcif_atoms are equivalent"""
            parse_mmtf.assertEqual(len(parse_mmtf.mmcif_atoms), len(parse_mmtf.mmtf_atoms))
            for i, e in enumerate(parse_mmtf.mmcif_atoms):
                mmtf_atom = parse_mmtf.mmtf_atoms[i]
                mmcif_atom = parse_mmtf.mmcif_atoms[i]
                parse_mmtf.assertEqual(mmtf_atom.name, mmcif_atom.name)  # eg. CA, spaces are removed from atom name
                parse_mmtf.assertEqual(mmtf_atom.fullname, mmcif_atom.fullname)  # e.g. " CA ", spaces included
                parse_mmtf.assertAlmostEqual(mmtf_atom.coord[0], mmcif_atom.coord[0], places=3)
                parse_mmtf.assertAlmostEqual(mmtf_atom.coord[1], mmcif_atom.coord[1], places=3)
                parse_mmtf.assertAlmostEqual(mmtf_atom.coord[2], mmcif_atom.coord[2], places=3)
                parse_mmtf.assertEqual(mmtf_atom.bfactor, mmcif_atom.bfactor)
                parse_mmtf.assertEqual(mmtf_atom.occupancy, mmcif_atom.occupancy)
                parse_mmtf.assertEqual(mmtf_atom.altloc, mmcif_atom.altloc)
                parse_mmtf.assertEqual(mmtf_atom.full_id,
                                       mmcif_atom.full_id)  # (structure id, model id, chain id, residue id, atom id)
                parse_mmtf.assertEqual(mmtf_atom.id, mmcif_atom.name)  # id of atom is the atom name (e.g. "CA")
                # self.assertEqual(mmtf_atom.serial_number,mmcif_atom.serial_number) # mmCIF serial number is none
        def test_residues(parse_mmtf):
            """Test that all residues in self.mmcif_res and self.mmtf_res are equivalent"""
            parse_mmtf.assertEqual(len(parse_mmtf.mmcif_res), len(parse_mmtf.mmtf_res))
            for i, e in enumerate(parse_mmtf.mmcif_res):
                mmcif_r = parse_mmtf.mmcif_res[i]
                mmtf_r = parse_mmtf.mmtf_res[i]
                parse_mmtf.assertEqual(mmtf_r.level, mmcif_r.level)
                parse_mmtf.assertEqual(mmtf_r.disordered, mmcif_r.disordered)
                parse_mmtf.assertEqual(mmtf_r.resname, mmcif_r.resname)
                parse_mmtf.assertEqual(mmtf_r.segid, mmcif_r.segid)
                parse_mmtf.mmcif_atoms = [x for x in mmcif_r.get_atom()]
                parse_mmtf.mmtf_atoms = [x for x in mmtf_r.get_atom()]
                test_atoms(parse_mmtf=parse_mmtf)

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', PDBConstructionWarning)
            mmtf_struct = MMTFParser.get_structure("PDB/4CUP.mmtf")
        mmcif_parser = MMCIFParser()
        mmcif_struct = mmcif_parser.get_structure("example", "PDB/4CUP.cif")
        self.mmcif_atoms = [x for x in mmcif_struct.get_atoms()]
        self.mmtf_atoms = [x for x in mmtf_struct.get_atoms()]
        test_atoms(self)
        mmcif_chains = [x for x in mmcif_struct.get_chains()]
        mmtf_chains = [x for x in mmtf_struct.get_chains()]
        self.assertEqual(len(mmcif_chains), len(mmtf_chains))
        for i, e in enumerate(mmcif_chains):
            self.mmcif_res = [x for x in mmcif_chains[i].get_residues()]
            self.mmtf_res = [x for x in mmtf_chains[i].get_residues()]
            test_residues(self)

        self.mmcif_res = [x for x in mmcif_struct.get_residues()]
        self.mmtf_res = [x for x in mmtf_struct.get_residues()]
        test_residues(self)
        self.assertEqual(len([x for x in mmcif_struct.get_models()]), len([x for x in mmtf_struct.get_models()]))



def test_parser():
    """Simply test that """
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', PDBConstructionWarning)
        structure = MMTFParser.get_structure("PDB/4CUP.mmtf")

if __name__ == '__main__':
    unittest.main()

