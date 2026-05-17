# Copyright 2024 by the Biopython team. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Unit tests for the Bio.Kinetics module."""

import math
import unittest

from Bio.Kinetics.Arrhenius import Arrenheius
from Bio.Kinetics.Arrhenius import firstOrder
from Bio.Kinetics.Arrhenius import secondOrder
from Bio.Kinetics.Arrhenius import zeroOrder


class TestArrenheius(unittest.TestCase):
    """Tests for the Arrenheius function."""

    def test_zero_activation_energy(self):
        """Returns A when Ea = 0."""
        self.assertEqual(Arrenheius(5e11, 0, 300), 5e11)

    def test_zero_prefactor(self):
        """Returns 0 when A = 0."""
        self.assertEqual(Arrenheius(0, 75000, 300), 0)

    def test_positive_values(self):
        """Returns positive result less than A."""
        result = Arrenheius(1e13, 75000, 298.15)
        self.assertGreater(result, 0)
        self.assertLess(result, 1e13)

    def test_matches_direct_formula(self):
        """Matches explicit formula A * exp(-Ea / (R * T))."""
        R = 8.314462618
        A, Ea, T = 1e10, 50000, 310.0
        expected = A * math.exp(-Ea / (R * T))
        self.assertAlmostEqual(Arrenheius(A, Ea, T), expected)

    def test_large_activation_energy(self):
        """Returns near zero for very large Ea."""
        result = Arrenheius(1e15, 500000, 300)
        self.assertAlmostEqual(result, 0.0)

    def test_high_temperature_increases_result(self):
        """Result increases with temperature."""
        A, Ea = 1e10, 75000
        low = Arrenheius(A, Ea, 500)
        high = Arrenheius(A, Ea, 5000)
        self.assertLess(low, high)
        self.assertGreater(high, low * 10)

    def test_negative_Ea_gives_greater_than_A(self):
        """Negative Ea gives result > A."""
        result = Arrenheius(1e10, -10000, 300)
        self.assertGreater(result, 1e10)


class TestFirstOrder(unittest.TestCase):
    """Tests for the firstOrder function (A = A0 * exp(-k*t))."""

    def test_zero_rate_constant(self):
        """Returns A0 when k = 0."""
        self.assertEqual(firstOrder(0, 100, 50), 100)

    def test_zero_time(self):
        """Returns A0 when t = 0."""
        self.assertEqual(firstOrder(0.5, 100, 0), 100)

    def test_zero_initial(self):
        """Returns 0 when A0 = 0."""
        self.assertEqual(firstOrder(0.5, 0, 50), 0)

    def test_decay(self):
        """Result decreases from A0 over time."""
        A0, k = 100, 0.1
        t1, t2 = 1, 10
        self.assertGreater(firstOrder(k, A0, t1), firstOrder(k, A0, t2))

    def test_matches_direct_formula(self):
        """Matches explicit formula A0 * exp(-k*t)."""
        k, A0, t = 0.05, 200, 30
        expected = A0 * math.exp(-k * t)
        self.assertAlmostEqual(firstOrder(k, A0, t), expected)

    def test_half_life(self):
        """At t = ln(2)/k, A = A0/2."""
        k = 0.1
        A0 = 100
        t_half = math.log(2) / k
        self.assertAlmostEqual(firstOrder(k, A0, t_half), A0 / 2)

    def test_large_t_approaches_zero(self):
        """Result approaches 0 for large t."""
        result = firstOrder(0.5, 100, 100)
        self.assertAlmostEqual(result, 0.0, places=10)

    def test_negative_k_increases(self):
        """Negative k gives growth instead of decay."""
        k, A0, t = -0.1, 100, 10
        self.assertGreater(firstOrder(k, A0, t), A0)


class TestSecondOrder(unittest.TestCase):
    """Tests for the secondOrder function (A = 1 / (1/A0 + k*t))."""

    def test_zero_rate_constant(self):
        """Returns A0 when k = 0."""
        self.assertEqual(secondOrder(0, 100, 50), 100)

    def test_zero_time(self):
        """Returns A0 when t = 0."""
        self.assertEqual(secondOrder(0.5, 100, 0), 100)

    def test_decay(self):
        """Result decreases from A0 over time."""
        A0, k = 100, 0.01
        t1, t2 = 1, 10
        self.assertGreater(secondOrder(k, A0, t1), secondOrder(k, A0, t2))

    def test_matches_direct_formula(self):
        """Matches explicit formula 1 / (1/A0 + k*t)."""
        k, A0, t = 0.02, 50, 5
        expected = 1 / (1 / A0 + k * t)
        self.assertAlmostEqual(secondOrder(k, A0, t), expected)

    def test_large_t_approaches_zero(self):
        """Result approaches 0 for large t."""
        result = secondOrder(10, 10, 1000000)
        self.assertAlmostEqual(result, 0.0, places=6)

    def test_zero_A0_raises(self):
        """Raises ZeroDivisionError when A0 = 0."""
        with self.assertRaises(ZeroDivisionError):
            secondOrder(0.1, 0, 5)


class TestZeroOrder(unittest.TestCase):
    """Tests for the zeroOrder function (A = A0 - k*t)."""

    def test_zero_rate_constant(self):
        """Returns A0 when k = 0."""
        self.assertEqual(zeroOrder(0, 100, 50), 100)

    def test_zero_time(self):
        """Returns A0 when t = 0."""
        self.assertEqual(zeroOrder(0.5, 100, 0), 100)

    def test_linear_decrease(self):
        """Result decreases linearly with time."""
        k, A0 = 2, 100
        self.assertEqual(zeroOrder(k, A0, 10), 80)
        self.assertEqual(zeroOrder(k, A0, 20), 60)
        self.assertEqual(zeroOrder(k, A0, 50), 0)

    def test_matches_direct_formula(self):
        """Matches explicit formula A0 - k*t."""
        k, A0, t = 0.5, 100, 30
        expected = A0 - k * t
        self.assertEqual(zeroOrder(k, A0, t), expected)

    def test_negative_after_depletion(self):
        """Returns negative when t > A0/k."""
        result = zeroOrder(10, 100, 20)
        self.assertLess(result, 0)

    def test_negative_k_increases(self):
        """Negative k gives linear increase."""
        k, A0, t = -2, 100, 10
        self.assertEqual(zeroOrder(k, A0, t), 120)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
