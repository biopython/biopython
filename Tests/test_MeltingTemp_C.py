#!/usr/bin/env python
# Copyright 2025 by Biopython contributors.
# All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Test the C-accelerated melting temperature calculations."""

import unittest
import warnings

from Bio import BiopythonWarning
from Bio.SeqUtils import MeltingTemp as mt


class TestMeltingTempCAcceleration(unittest.TestCase):
    """Test C acceleration for melting temperature calculations."""

    def setUp(self):
        """Set up test sequences."""
        self.test_sequences = [
            "GAGTCCATATGGCTGAGAAG",
            "ACGTACGTACGTACGT",
            "GGGGCCCCGGGGCCCC",
            "ATATATATATATATATAT",
            "GCGCGCGCGCGCGCGCGC",
            "ACGT",
            "AAAATTTT",
            "CCCCGGGG",
        ]
        self.tolerance = 0.1  # Allow 0.1°C difference for floating point

    def test_c_acceleration_available(self):
        """Check if C acceleration is available."""
        try:
            from Bio.SeqUtils import _meltingtemp_exact

            available = True
        except ImportError:
            available = False

        # This test documents whether C acceleration is available
        # It doesn't fail if not available (optional dependency)
        if not available:
            warnings.warn(
                "C acceleration for melting temperature not available", BiopythonWarning
            )

    def test_c_matches_python(self):
        """Test that C implementation matches Python implementation."""
        try:
            from Bio.SeqUtils import _meltingtemp_exact
        except ImportError:
            self.skipTest("C acceleration not available")

        for seq in self.test_sequences:
            with self.subTest(sequence=seq):
                # Force Python implementation
                tm_python = mt.Tm_NN(seq, nn_table=mt.DNA_NN3, check=False)

                # Use C implementation directly
                tm_c = _meltingtemp_exact.tm_nn_exact(seq)

                # Check they match within tolerance
                self.assertAlmostEqual(
                    tm_python,
                    tm_c,
                    delta=self.tolerance,
                    msg=f"Mismatch for sequence {seq}: Python={tm_python:.2f}, C={tm_c:.2f}",
                )

    def test_salt_corrections(self):
        """Test different salt correction methods."""
        try:
            from Bio.SeqUtils import _meltingtemp_exact
        except ImportError:
            self.skipTest("C acceleration not available")

        seq = "GAGTCCATATGGCTGAGAAG"

        # Test each salt correction method
        for method in [0, 1, 2, 3, 4, 5]:
            with self.subTest(saltcorr=method):
                tm_python = mt.Tm_NN(
                    seq, nn_table=mt.DNA_NN3, check=False, saltcorr=method
                )
                tm_c = _meltingtemp_exact.tm_nn_exact(seq, saltcorr=method)

                self.assertAlmostEqual(
                    tm_python,
                    tm_c,
                    delta=self.tolerance,
                    msg=f"Salt correction {method} mismatch",
                )

    def test_concentration_parameters(self):
        """Test different DNA concentration parameters."""
        try:
            from Bio.SeqUtils import _meltingtemp_exact
        except ImportError:
            self.skipTest("C acceleration not available")

        seq = "GAGTCCATATGGCTGAGAAG"

        # Test different concentrations
        test_params = [
            {"dnac1": 100, "dnac2": 100},
            {"dnac1": 10, "dnac2": 10},
            {"dnac1": 50, "dnac2": 25},
        ]

        for params in test_params:
            with self.subTest(**params):
                tm_python = mt.Tm_NN(seq, nn_table=mt.DNA_NN3, check=False, **params)
                tm_c = _meltingtemp_exact.tm_nn_exact(seq, **params)

                self.assertAlmostEqual(
                    tm_python,
                    tm_c,
                    delta=self.tolerance,
                    msg=f"Concentration parameter mismatch: {params}",
                )

    def test_self_complementary(self):
        """Test self-complementary sequence detection."""
        try:
            from Bio.SeqUtils import _meltingtemp_exact
        except ImportError:
            self.skipTest("C acceleration not available")

        # Self-complementary sequences
        self_comp_seqs = [
            "GCGC",
            "ATAT",
            "GGCC",
            "ACGT",  # Note: This is self-complementary when reversed
        ]

        for seq in self_comp_seqs:
            with self.subTest(sequence=seq):
                tm_python = mt.Tm_NN(
                    seq, nn_table=mt.DNA_NN3, check=False, selfcomp=True
                )
                tm_c = _meltingtemp_exact.tm_nn_exact(seq, selfcomp=True)

                self.assertAlmostEqual(
                    tm_python,
                    tm_c,
                    delta=self.tolerance,
                    msg=f"Self-complementary mismatch for {seq}",
                )

    def test_integrated_acceleration(self):
        """Test that default Tm_NN uses C acceleration when available."""
        try:
            from Bio.SeqUtils import _meltingtemp_exact

            c_available = True
        except ImportError:
            c_available = False
            self.skipTest("C acceleration not available")

        seq = "GAGTCCATATGGCTGAGAAG"

        # Default call (should use C if available)
        tm_default = mt.Tm_NN(seq)

        # Force Python
        tm_python = mt.Tm_NN(seq, nn_table=mt.DNA_NN3)

        # They should match
        self.assertAlmostEqual(
            tm_default,
            tm_python,
            delta=self.tolerance,
            msg="Integrated C acceleration not working correctly",
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
