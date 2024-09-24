"""Tests for PDB Mutate module."""

import unittest
from Bio.PDB.Mutate import mutate
from Bio.PDB import PDBParser
from Bio.PDB import Mutate


class MutateTest(unittest.TestCase):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("my_structure", "Mutate/5e0m.pdb")
    mutate(structure, chain="A", mutate_to="TYR", res_num=44, mutation_type="best")

    def test_point_mutation(self):
        self.assertEqual(self.structure[0]["A"][44].get_resname(), "TYR")
        # self.assertRaises(ValueError, Mutate.mutate, self.structure, chain="A", mutate_to="TYR", res_num=44, mutation_type="best")

    def test_position(self):
        self.assertAlmostEqual(
            self.structure[0]["A"][44]["CD1"].coord[0], -27.298, places=3
        )
        self.assertAlmostEqual(
            self.structure[0]["A"][44]["CD1"].coord[1], 8.571, places=3
        )
        self.assertAlmostEqual(
            self.structure[0]["A"][44]["CD1"].coord[2], -5.160, places=3
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
