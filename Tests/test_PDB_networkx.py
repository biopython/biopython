import unittest
from Bio.PDB.Atom import Atom
from Bio.PDB.networks import to_networkx

class NetworkXExportTest(unittest.TestCase):
    def test_to_networkx_conversion(self):
        """Test conversion of atom pairs to weighted NetworkX graph."""
        a1 = Atom("CA", [0, 0, 0], 0, 0, " ", "CA", 1, "C")
        a2 = Atom("CA", [0, 0, 3], 0, 0, " ", "CA", 2, "C")
        neighbor_pairs = [(a1, a2)]

        G = to_networkx(neighbor_pairs)

        self.assertEqual(len(G.nodes), 2)
        self.assertEqual(len(G.edges), 1)
        self.assertAlmostEqual(G[a1][a2]['weight'], 3.0)

if __name__ == "__main__":
    unittest.main()
