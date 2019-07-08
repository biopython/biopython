# Copyright 2017 by Francesco Gastaldello. All rights reserved.
#
# Converted by Francesco Gastaldello from an older unit test copyright 2004
# by Thomas Hamelryck.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Unit tests for the Bio.PDB.Polypetide module."""

import unittest

from Bio.PDB import PDBParser, PPBuilder, CaPPBuilder


class PolypeptideTests(unittest.TestCase):
    """Test Polypetide module."""

    def test_polypeptide(self):
        """Tests on polypetide class and methods."""
        p = PDBParser(PERMISSIVE=True)
        pdb1 = "PDB/1A8O.pdb"
        s = p.get_structure("scr", pdb1)
        ppb = PPBuilder()
        pp = ppb.build_peptides(s)
        self.assertEqual(str(pp[0].get_sequence()), "DIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNW")
        self.assertEqual(str(pp[1].get_sequence()), "TETLLVQNANPDCKTILKALGPGATLEE")
        self.assertEqual(str(pp[2].get_sequence()), "TACQG")
        phi_psi = pp[0].get_phi_psi_list()
        self.assertEqual(phi_psi[0][0], None)
        self.assertAlmostEqual(phi_psi[0][1], -0.46297171497725553, places=3)
        self.assertAlmostEqual(phi_psi[1][0], -1.0873937604007962, places=3)
        self.assertAlmostEqual(phi_psi[1][1], 2.1337707832637109, places=3)
        self.assertAlmostEqual(phi_psi[2][0], -2.4052232743651878, places=3)
        self.assertAlmostEqual(phi_psi[2][1], 2.3807316946081554, places=3)
        phi_psi = pp[1].get_phi_psi_list()
        self.assertEqual(phi_psi[0][0], None)
        self.assertAlmostEqual(phi_psi[0][1], -0.6810077089092923, places=3)
        self.assertAlmostEqual(phi_psi[1][0], -1.2654003477656888, places=3)
        self.assertAlmostEqual(phi_psi[1][1], -0.58689987042756309, places=3)
        self.assertAlmostEqual(phi_psi[2][0], -1.7467679151684763, places=3)
        self.assertAlmostEqual(phi_psi[2][1], -1.5655066256698336, places=3)
        phi_psi = pp[2].get_phi_psi_list()
        self.assertEqual(phi_psi[0][0], None)
        self.assertAlmostEqual(phi_psi[0][1], -0.73222884210889716, places=3)
        self.assertAlmostEqual(phi_psi[1][0], -1.1044740234566259, places=3)
        self.assertAlmostEqual(phi_psi[1][1], -0.69681334592782884, places=3)
        self.assertAlmostEqual(phi_psi[2][0], -1.8497413300164958, places=3)
        self.assertAlmostEqual(phi_psi[2][1], 0.34762889834809058, places=3)
        ppb = CaPPBuilder()
        pp = ppb.build_peptides(s)
        self.assertEqual(str(pp[0].get_sequence()), "DIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNW")
        self.assertEqual(str(pp[1].get_sequence()), "TETLLVQNANPDCKTILKALGPGATLEE")
        self.assertEqual(str(pp[2].get_sequence()), "TACQG")
        self.assertEqual([ca.serial_number for ca in pp[0].get_ca_list()], [10, 18, 26, 37, 46, 50, 57, 66, 75, 82, 93, 104, 112, 124, 131, 139, 150, 161, 173, 182, 189, 197, 208, 213, 222, 231, 236, 242, 251, 260, 267, 276, 284])
        taus = pp[1].get_tau_list()
        self.assertAlmostEqual(taus[0], 0.3597907225123525, places=3)
        self.assertAlmostEqual(taus[1], 0.43239284636769254, places=3)
        self.assertAlmostEqual(taus[2], 0.99820157492712114, places=3)
        thetas = pp[2].get_theta_list()
        self.assertAlmostEqual(thetas[0], 1.6610069445335354, places=3)
        self.assertAlmostEqual(thetas[1], 1.7491703334817772, places=3)
        self.assertAlmostEqual(thetas[2], 2.0702447422720143, places=3)


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
