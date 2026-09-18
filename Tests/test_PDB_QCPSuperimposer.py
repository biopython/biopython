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


from Bio.PDB import PDBParser
from Bio.PDB import Selection
from Bio.PDB.qcprot import QCPSuperimposer
from Bio.SVDSuperimposer import SVDSuperimposer


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

    # Public methods
    def test_set(self):
        """Test setting of initial parameters."""
        sup = QCPSuperimposer()
        sup.set(self.x, self.y)

        self.assertTrue(np.allclose(sup.reference_coords, self.x, atol=1e-6))
        self.assertTrue(np.allclose(sup.coords, self.y, atol=1e-6))
        self.assertIsNone(sup.transformed_coords)
        self.assertIsNone(sup.rot)
        self.assertIsNone(sup.tran)
        self.assertIsNone(sup.rms)
        self.assertIsNone(sup.init_rms)

    def test_run(self):
        """Test QCP on dummy data."""
        sup = QCPSuperimposer()
        sup.set(self.x, self.y)

        sup.run()
        self.assertTrue(np.allclose(sup.reference_coords, self.x, atol=1e-6))
        self.assertTrue(np.allclose(sup.coords, self.y, atol=1e-6))
        self.assertIsNone(sup.transformed_coords)

        calc_rot = [
            [0.683, 0.537, 0.495],
            [-0.523, 0.833, -0.181],
            [-0.510, -0.135, 0.849],
        ]
        calc_tran = [38.786, -20.655, -15.422]

        self.assertTrue(np.allclose(np.array(calc_rot), sup.rot, atol=1e-3))
        self.assertTrue(np.allclose(np.array(calc_tran), sup.tran, atol=1e-3))

        # We can reduce precision here since we do a similar calculation
        # for a full structure down below.
        self.assertAlmostEqual(sup.rms, 0.003, places=3)
        self.assertIsNone(sup.init_rms)

    def test_compare_to_svd(self):
        """Compare results of QCP to SVD."""
        sup = QCPSuperimposer()
        sup.set(self.x, self.y)
        sup.run()

        svd_sup = SVDSuperimposer()
        svd_sup.set(self.x, self.y)
        svd_sup.run()

        self.assertAlmostEqual(svd_sup.get_rms(), sup.rms, places=3)
        self.assertTrue(np.allclose(svd_sup.rot, sup.rot, atol=1e-3))
        self.assertTrue(np.allclose(svd_sup.tran, sup.tran, atol=1e-3))
        self.assertTrue(
            np.allclose(svd_sup.get_transformed(), sup.get_transformed(), atol=1e-3)
        )

    def test_compare_to_svd_lines(self):
        """Compare results of QCP to SVD using simple lines."""
        ref = np.linspace([0, 0, 0], [10, 10, 10], 5)
        mob = np.linspace([20, 10, 30], [20, 20, 30], 5)

        sup = QCPSuperimposer()
        sup.set(ref, mob)
        sup.run()

        svd_sup = SVDSuperimposer()
        svd_sup.set(ref, mob)
        svd_sup.run()

        self.assertAlmostEqual(svd_sup.get_rms(), sup.rms, places=3)
        self.assertTrue(
            np.allclose(svd_sup.get_transformed(), sup.get_transformed(), atol=1e-3)
        )

    def test_get_transformed(self):
        """Test transformation of coordinates after QCP."""
        sup = QCPSuperimposer()
        sup.set(self.x, self.y)

        sup.run()
        transformed_coords = [
            [51.652, -1.900, 50.071],
            [50.398, -1.229, 50.649],
            [50.680, -0.042, 51.537],
            [50.220, -0.019, 52.853],
        ]

        self.assertTrue(
            np.allclose(sup.get_transformed(), np.array(transformed_coords), atol=1e-3)
        )

    def test_get_init_rms(self):
        """Test initial RMS calculation."""
        x = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        y = np.array([[2.0, 0.0, 0.0], [1.0, 1.0, 0.0], [1.0, 0.0, 1.0]])

        sup = QCPSuperimposer()
        sup.set(x, y)
        self.assertIsNone(sup.init_rms)

        expected_init_rms = 1.0
        self.assertAlmostEqual(sup.get_init_rms(), expected_init_rms, places=6)

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
        self.assertTrue(np.allclose(sup.rotran[0], rot, atol=1e-3))
        self.assertTrue(np.allclose(sup.rotran[1], -tran, atol=1e-3))
        self.assertAlmostEqual(sup.rms, 0.0, places=6)

    def test_compare_rmsd_to_transformed(self):
        """Compare RMSD from QCP algorithm to that from transformed"""
        ref = np.linspace([0, 0, 0], [10, 10, 10], 5)
        mob = np.linspace([20, 10, 30], [20, 20, 30], 5)

        sup = QCPSuperimposer()
        sup.set(ref, mob)
        sup.run()
        rms = sup.get_rms()
        mob_fitted = sup.get_transformed()
        rms_fitted = np.sqrt(((ref - mob_fitted) ** 2).sum() / ref.shape[0])
        self.assertAlmostEqual(rms, rms_fitted, places=6)

    # Tests for issue #5309: collinear / rank-deficient reference
    # coordinates produced a wrong rotation while reporting RMSD ~ 0.
    # The fix routes degenerate reference sets through SVDSuperimposer.

    def _random_rigid_transform(self, rng):
        """Return (rotation, translation) suitable for testing."""
        axis = rng.normal(size=3)
        axis /= np.linalg.norm(axis)
        theta = rng.uniform(0.0, np.pi)
        K = np.array(
            [
                [0.0, -axis[2], axis[1]],
                [axis[2], 0.0, -axis[0]],
                [-axis[1], axis[0], 0.0],
            ]
        )
        rot = np.eye(3) + np.sin(theta) * K + (1.0 - np.cos(theta)) * (K @ K)
        tran = rng.normal(scale=10.0, size=3)
        return rot, tran

    def test_collinear_reference(self):
        """Collinear reference must not produce a wrong rotation (#5309)."""
        rng = np.random.default_rng(2026)
        rot, tran = self._random_rigid_transform(rng)
        ref = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
        mob = ref @ rot.T + tran

        sup = QCPSuperimposer()
        sup.set(ref, mob)
        sup.run()

        # The reported RMSD must match the true residual RMSD of the
        # transformed coordinates -- the bug was RMSD ~ 0 while the
        # actual alignment was off.
        transformed = sup.get_transformed()
        true_rmsd = np.sqrt(((ref - transformed) ** 2).sum() / ref.shape[0])
        self.assertAlmostEqual(sup.rms, true_rmsd, places=10)
        self.assertLess(true_rmsd, 1e-9)

        # Should agree with SVDSuperimposer on this degenerate input.
        svd_sup = SVDSuperimposer()
        svd_sup.set(ref, mob)
        svd_sup.run()
        self.assertAlmostEqual(svd_sup.get_rms(), sup.rms, places=10)

    def test_collinear_reference_off_axis(self):
        """Collinear reference along a non-axis direction (#5309)."""
        rng = np.random.default_rng(2027)
        rot, tran = self._random_rigid_transform(rng)
        # Direction (1, 2, 3), normalised.
        d = np.array([1.0, 2.0, 3.0])
        d /= np.linalg.norm(d)
        ref = np.stack([np.zeros(3), d, 2.0 * d])
        mob = ref @ rot.T + tran

        sup = QCPSuperimposer()
        sup.set(ref, mob)
        sup.run()

        transformed = sup.get_transformed()
        true_rmsd = np.sqrt(((ref - transformed) ** 2).sum() / ref.shape[0])
        self.assertAlmostEqual(sup.rms, true_rmsd, places=10)
        self.assertLess(true_rmsd, 1e-9)

    def test_coplanar_reference_unchanged(self):
        """Coplanar (rank 2) reference must still align correctly."""
        rng = np.random.default_rng(2028)
        rot, tran = self._random_rigid_transform(rng)
        ref = np.array(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]]
        )
        mob = ref @ rot.T + tran

        sup = QCPSuperimposer()
        sup.set(ref, mob)
        sup.run()
        true_rmsd = np.sqrt(((ref - sup.get_transformed()) ** 2).sum() / 4)
        self.assertLess(true_rmsd, 1e-9)

    def test_coincident_points(self):
        """All reference points identical -> degenerate, no crash."""
        ref = np.zeros((3, 3))
        mob = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])

        sup = QCPSuperimposer()
        sup.set(ref, mob)
        # run() should not raise; result is degenerate but well-defined.
        sup.run()
        self.assertIsNotNone(sup.rot)
        self.assertIsNotNone(sup.tran)

    def test_full_rank_reference_unchanged(self):
        """Full-rank reference still uses the fast QCP path (no regression).

        We check that QCP and SVDSuperimposer agree on the existing
        well-conditioned 4-point fixture, and on 50 random rigid
        rotations to guard against any regression in the common path.
        """
        # Existing fixture: must match SVDSuperimposer to 1e-3.
        sup = QCPSuperimposer()
        sup.set(self.x, self.y)
        sup.run()
        svd_sup = SVDSuperimposer()
        svd_sup.set(self.x, self.y)
        svd_sup.run()
        self.assertAlmostEqual(svd_sup.get_rms(), sup.rms, places=3)
        self.assertTrue(np.allclose(svd_sup.rot, sup.rot, atol=1e-3))

        # Random regression sweep: QCP and SVD must agree on
        # well-conditioned inputs.  QCP is a closed-form quaternion
        # solver with ~1e-8 relative error per call, so we don't
        # tighten past 1e-6.
        rng = np.random.default_rng(2029)
        ref = rng.normal(size=(20, 3))
        for _ in range(50):
            rot, tran = self._random_rigid_transform(rng)
            mob = ref @ rot.T + tran

            q = QCPSuperimposer()
            q.set(ref, mob)
            q.run()
            s = SVDSuperimposer()
            s.set(ref, mob)
            s.run()
            self.assertAlmostEqual(q.rms, s.get_rms(), places=6)
            self.assertTrue(np.allclose(q.rot, s.rot, atol=1e-6))
            self.assertTrue(np.allclose(q.tran, s.tran, atol=1e-6))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
