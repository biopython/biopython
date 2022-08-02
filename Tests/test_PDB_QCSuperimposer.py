# Copyright 2017 by Maximilian Greil.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for QCPSuperimposer module."""

import unittest

try:
    import numpy as np
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.QCPSuperimposer."
    ) from None


from Bio.PDB import PDBParser, Selection
from Bio.PDB.qcprot import QCPSuperimposer


class QCPSuperimposerTest(unittest.TestCase):
    def setUp(self):
        self.x = np.array(
            [
                [51.65, -1.90, 50.07],
                [50.40, -1.23, 50.65],
                [50.68, -0.04, 51.54],
                [50.22, -0.02, 52.85],
            ]
        )

        self.y = np.array(
            [
                [51.30, -2.99, 46.54],
                [51.09, -1.88, 47.58],
                [52.36, -1.20, 48.03],
                [52.71, -1.18, 49.38],
            ]
        )

        self.sup = QCPSuperimposer()
        self.sup.set(self.x, self.y)

    def _arr_to_list(self, a):
        """Return the array cast as a list, rounded to 3 decimals."""
        return np.around(a, decimals=3).tolist()

    # Public methods
    def test_set(self):
        """Test setting of initial parameters."""
        self.assertEqual(
            self._arr_to_list(self.sup.reference_coords), self._arr_to_list(self.x)
        )
        self.assertEqual(self._arr_to_list(self.sup.coords), self._arr_to_list(self.y))
        self.assertIsNone(self.sup.transformed_coords)
        self.assertIsNone(self.sup.rot)
        self.assertIsNone(self.sup.tran)
        self.assertIsNone(self.sup.rms)
        self.assertIsNone(self.sup.init_rms)

    def test_run(self):
        """Test QCP on dummy data."""
        self.sup.run()
        self.assertEqual(
            self._arr_to_list(self.sup.reference_coords), self._arr_to_list(self.x)
        )
        self.assertEqual(self._arr_to_list(self.sup.coords), self._arr_to_list(self.y))
        self.assertIsNone(self.sup.transformed_coords)
        calc_rot = [
            [0.683, 0.537, 0.495],
            [-0.523, 0.833, -0.181],
            [-0.510, -0.135, 0.849],
        ]
        self.assertEqual(self._arr_to_list(self.sup.rot), calc_rot)
        calc_tran = [38.786, -20.655, -15.422]
        self.assertEqual(self._arr_to_list(self.sup.tran), calc_tran)
        # We can reduce precision here since we do a similar calculation
        # for a full structure down below.
        self.assertAlmostEqual(self.sup.rms, 0.003, places=2)
        self.assertIsNone(self.sup.init_rms)

    def test_get_transformed(self):
        """Test transformation of coordinates after QCP."""
        self.sup.run()
        transformed_coords = [
            [51.652, -1.900, 50.071],
            [50.398, -1.229, 50.649],
            [50.680, -0.042, 51.537],
            [50.220, -0.019, 52.853],
        ]

        self.assertEqual(
            self._arr_to_list(self.sup.get_transformed()), transformed_coords
        )

    def test_get_init_rms(self):
        """Test initial RMS calculation."""
        x = np.array([[1.1, 1.2, 1.3], [1.4, 1.5, 1.6], [1.7, 1.8, 1.9]])
        y = np.array([[1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]])
        self.sup.set(x, y)
        self.assertIsNone(self.sup.init_rms)
        init_rms = [0.812, 0.90, 0.98]
        self.assertEqual(
            self._arr_to_list(self.sup.get_init_rms()),
            init_rms,
        )

    def test_on_pdb(self):
        """Align a PDB to itself."""
        pdb1 = "PDB/1A8O.pdb"
        p = PDBParser()
        s1 = p.get_structure("FIXED", pdb1)
        fixed = Selection.unfold_entities(s1, "A")
        s2 = p.get_structure("MOVING", pdb1)
        moving = Selection.unfold_entities(s2, "A")

        rot = np.eye(3, dtype=np.float64)
        tran = np.array([1.0, 2.0, 3.0], dtype=np.float64)
        for atom in moving:
            atom.transform(rot, tran)

        sup = QCPSuperimposer()
        sup.set_atoms(fixed, moving)
        self.assertEqual(self._arr_to_list(sup.rotran[0]), self._arr_to_list(rot))
        self.assertEqual(self._arr_to_list(sup.rotran[1]), self._arr_to_list(-tran))
        self.assertAlmostEqual(sup.rms, 0.0, places=6)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
