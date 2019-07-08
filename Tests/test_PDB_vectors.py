# Copyright 2017 by Peter Cock.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Tests for the vector code in Bio.PDB."""

import unittest

try:
    import numpy
    from numpy.random import random
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.PDB.")

from Bio.PDB.vectors import Vector
from Bio.PDB import rotmat, refmat, calc_angle, calc_dihedral, rotaxis, m2rotaxis


class VectorTests(unittest.TestCase):
    """Tests for the Vector class."""

    def test_division(self):
        """Confirm division works."""
        v = Vector(1, 1, 1) / 2
        self.assertEqual(repr(v), "<Vector 0.50, 0.50, 0.50>")

    def test_Vector(self):
        """Test Vector object."""
        v1 = Vector(0, 0, 1)
        v2 = Vector(0, 0, 0)
        v3 = Vector(0, 1, 0)
        v4 = Vector(1, 1, 0)

        self.assertEqual(calc_angle(v1, v2, v3), 1.5707963267948966)
        self.assertEqual(calc_dihedral(v1, v2, v3, v4), 1.5707963267948966)
        self.assertTrue(numpy.array_equal((v1 - v2).get_array(), numpy.array([0.0, 0.0, 1.0])))
        self.assertTrue(numpy.array_equal((v1 - 1).get_array(), numpy.array([-1.0, -1.0, 0.0])))
        self.assertTrue(numpy.array_equal((v1 - (1, 2, 3)).get_array(), numpy.array([-1.0, -2.0, -2.0])))
        self.assertTrue(numpy.array_equal((v1 + v2).get_array(), numpy.array([0.0, 0.0, 1.0])))
        self.assertTrue(numpy.array_equal((v1 + 3).get_array(), numpy.array([3.0, 3.0, 4.0])))
        self.assertTrue(numpy.array_equal((v1 + (1, 2, 3)).get_array(), numpy.array([1.0, 2.0, 4.0])))
        self.assertTrue(numpy.array_equal(v1.get_array() / 2, numpy.array([0, 0, 0.5])))
        self.assertTrue(numpy.array_equal(v1.get_array() / 2, numpy.array([0, 0, 0.5])))
        self.assertEqual(v1 * v2, 0.0)
        self.assertTrue(numpy.array_equal((v1 ** v2).get_array(), numpy.array([0.0, -0.0, 0.0])))
        self.assertTrue(numpy.array_equal((v1 ** 2).get_array(), numpy.array([0.0, 0.0, 2.0])))
        self.assertTrue(numpy.array_equal((v1 ** (1, 2, 3)).get_array(), numpy.array([0.0, 0.0, 3.0])))
        self.assertEqual(v1.norm(), 1.0)
        self.assertEqual(v1.normsq(), 1.0)
        v1[2] = 10
        self.assertEqual(v1.__getitem__(2), 10)

    def test_normalization(self):
        """Test Vector normalization."""
        v1 = Vector([2, 0, 0])
        self.assertTrue(numpy.array_equal(v1.normalized().get_array(), numpy.array([1, 0, 0])))
        # State of v1 should not be affected by `normalized`
        self.assertTrue(numpy.array_equal(v1.get_array(), numpy.array([2, 0, 0])))
        v1.normalize()
        # State of v1 should be affected by `normalize`
        self.assertTrue(numpy.array_equal(v1.get_array(), numpy.array([1, 0, 0])))

    def test_refmat(self):
        v1 = Vector(0, 0, 1)
        v2 = Vector(0, 1, 0)
        ref = refmat(v1, v2)
        self.assertTrue(numpy.allclose(ref[0], [1.0, 0.0, 0.0]))
        self.assertTrue(numpy.allclose(ref[1], [0.0, 0.0, 1.0]))
        self.assertTrue(numpy.allclose(ref[2], [0.0, 1.0, 0.0]))
        self.assertTrue(numpy.allclose(v1.left_multiply(ref).get_array(), [0.0, 1.0, 0.0]))

    def test_rotmat_90(self):
        """Test regular 90 deg rotation."""
        v1 = Vector(0, 0, 1)
        v2 = Vector(0, 1, 0)
        rot = rotmat(v1, v2)
        self.assertTrue(numpy.allclose(rot[0], numpy.array([1.0, 0.0, 0.0])))
        self.assertTrue(numpy.allclose(rot[1], numpy.array([0.0, 0.0, 1.0])))
        self.assertTrue(numpy.allclose(rot[2], numpy.array([0.0, -1.0, 0.0])))
        self.assertTrue(numpy.allclose(v1.left_multiply(rot).get_array(), [0.0, 1.0, 0.0]))
        self.assertTrue(numpy.allclose(v1.right_multiply(numpy.transpose(rot)).get_array(), [0.0, 1.0, 0.0]))

    def test_rotmat_180(self):
        """Test rotmat when the rotation is 180 deg (singularity)."""
        v1 = Vector([1.0, 0.8, 0])
        v2 = Vector([-1.0, -0.8, 0])
        rot = rotmat(v1, v2)
        v3 = v1.left_multiply(rot)
        self.assertTrue(numpy.allclose(v2.get_array(), v3.get_array()))

    def test_rotmat_0(self):
        """Test rotmat when the rotation is 0 deg (singularity)."""
        v1 = Vector([1.0, 0.8, 0])
        v2 = Vector([1.0, 0.8, 0])
        rot = rotmat(v1, v2)
        v3 = v1.left_multiply(rot)
        self.assertTrue(numpy.allclose(v1.get_array(), v3.get_array()))

    def test_m2rotaxis_90(self):
        """Test 90 deg rotation."""
        v1 = Vector(0, 0, 1)
        v2 = Vector(0, 1, 0)
        rot = rotmat(v1, v2)
        angle, axis = m2rotaxis(rot)
        self.assertTrue(numpy.allclose(axis.get_array(), [-1.0, 0.0, 0.0]))
        self.assertTrue(abs(angle - numpy.pi / 2) < 1e-5)

    def test_m2rotaxis_180(self):
        """Test 180 deg rotation."""
        v1 = Vector([1.0, 0.8, 0])
        v2 = Vector([-1.0, -0.8, 0])
        rot = rotmat(v1, v2)
        angle, axis = m2rotaxis(rot)
        self.assertTrue(abs(axis * v1) < 1e-5)  # axis orthogonal to v1
        self.assertTrue(abs(angle - numpy.pi) < 1e-5)

    def test_m2rotaxis_0(self):
        """Test 0 deg rotation. Axis must be [1, 0, 0] as per Vector documentation."""
        v1 = Vector([1.0, 0.8, 0])
        v2 = Vector([1.0, 0.8, 0])
        rot = rotmat(v1, v2)
        angle, axis = m2rotaxis(rot)
        self.assertTrue(numpy.allclose(axis.get_array(), [1, 0, 0]))
        self.assertTrue(abs(angle) < 1e-5)

    def test_Vector_angles(self):
        """Test Vector angles."""
        angle = random() * numpy.pi
        axis = Vector(random(3) - random(3))
        axis.normalize()
        m = rotaxis(angle, axis)
        cangle, caxis = m2rotaxis(m)
        self.assertAlmostEqual(angle, cangle, places=3)
        self.assertTrue(numpy.allclose(list(map(int, (axis - caxis).get_array())), [0, 0, 0]),
                        "Want %r and %r to be almost equal" % (axis.get_array(), caxis.get_array()))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
