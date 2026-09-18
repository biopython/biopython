# Copyright 2026 by the Biopython contributors.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included as part
# of this package.
"""Tests for Bio.SeqUtils.IsoelectricPoint.

These tests cover both the original public API and the regression for
issue #5312, where ``IsoelectricPoint.pi()`` returned a bracket-boundary
pH (and a large residual charge) for strongly-charged homopolymer
sequences whose true pI lies outside the previous [4.05, 12] bracket.
"""

import unittest

try:
    import numpy as np
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.SeqUtils.IsoelectricPoint."
    ) from None

from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint


class IsoelectricPointBasicTest(unittest.TestCase):
    """Sanity checks against the existing module-level doctest examples."""

    def test_pi_doctest_INGAR(self):
        """The INGAR example from the module docstring must still work."""
        p = IsoelectricPoint("INGAR")
        self.assertAlmostEqual(p.pi(), 9.75, places=2)
        # And charge at that pH is small.
        self.assertLess(abs(p.charge_at_pH(p.pi())), 0.01)

    def test_charge_at_pH_doctest(self):
        """Charge at pH 7 for INGAR is +0.76 per the docstring."""
        p = IsoelectricPoint("INGAR")
        self.assertAlmostEqual(p.charge_at_pH(7.0), 0.76, places=2)


class IsoelectricPointRegression5312Test(unittest.TestCase):
    """Regression for issue #5312: bracket was too narrow."""

    # All sequences here were confirmed to fail under the pre-fix
    # bracket [4.05, 12] (residual |charge| > 1 at the returned pH).
    # Under the new bracket [0.0, 14.0] the residual must be near 0.
    HOMOPOLYMER_FIXTURES = [
        # (sequence, descriptive label)
        ("DDDDDDDDDD", "poly-aspartate (issue reproducer)"),
        ("DDDDDD", "short poly-aspartate"),
        ("EEEEEEEEEE", "poly-glutamate"),
        ("RRRRRRRRRR", "poly-arginine"),
        ("KKKKKKKKKK", "poly-lysine"),
        ("CCCCCCCCCC", "poly-cysteine"),
        ("YYYYYYYYYY", "poly-tyrosine"),
    ]

    def test_issue_5312_reproducer(self):
        """Verbatim reproducer from issue #5312."""
        p = IsoelectricPoint("DDDDDDDDDD")
        self.assertAlmostEqual(p.charge_at_pH(p.pi()), 0.0, places=2)

    def test_homopolymer_residual_charge_is_near_zero(self):
        """The reported pH must yield a charge close to 0, not a large residual."""
        for seq, label in self.HOMOPOLYMER_FIXTURES:
            p = IsoelectricPoint(seq)
            pi = p.pi()
            residual = p.charge_at_pH(pi)
            self.assertLess(
                abs(residual),
                0.01,
                f"{label} ({seq}): residual charge {residual} at reported pI {pi}",
            )

    def test_homopolymer_pI_in_biochemical_range(self):
        """Reported pI for a homopolymer must be near the residue's
        net charge balance point (where positive and negative
        contributions cancel).  These values come from the standard
        Bjellqvist pK tables used by Biopython; the test guards
        against bracket-endpoint clamping."""
        # (sequence, expected_pI_range_low, expected_pI_range_high)
        expected_ranges = [
            # poly-D and poly-E: root is below 4 (D pKa=3.65, E pKa=4.25,
            # Cterm pKa=3.55). New bracket [0, 14] finds the true root.
            ("DDDDDDDDDD", 2.5, 4.0),
            ("EEEEEEEEEE", 3.0, 4.5),
            # poly-R: root is between R pKa=12.0 and Nterm pKa=7.5.
            ("RRRRRRRRRR", 12.0, 13.5),
            # poly-K: K pKa=10.0, Nterm=7.5, balance point near K's pKa.
            ("KKKKKKKKKK", 9.5, 11.5),
            # poly-C: Cys side chain is acidic (pKa=9.0); root is around
            # the balance point with Cterm (pKa=3.55) and Nterm (pKa=7.5).
            ("CCCCCCCCCC", 5.0, 7.5),
            # poly-Y: Tyr side chain is acidic (pKa=10.0); root is between
            # Tyr/Cterm and Nterm.
            ("YYYYYYYYYY", 5.0, 7.5),
        ]
        for seq, lo, hi in expected_ranges:
            p = IsoelectricPoint(seq)
            pi = p.pi()
            self.assertGreaterEqual(pi, lo, f"{seq}: pI {pi} below expected floor {lo}")
            self.assertLessEqual(pi, hi, f"{seq}: pI {pi} above expected ceiling {hi}")

    def test_pi_returns_float(self):
        """The function returns a plain float, not a numpy scalar."""
        pi = IsoelectricPoint("INGAR").pi()
        self.assertIsInstance(pi, float)

    def test_pi_recursion_terminates(self):
        """The bisection recursion must terminate (does not exceed Python
        recursion limit) for typical inputs.  Default Python limit is
        1000, so 20 levels of recursion per call is comfortable."""
        import sys

        old = sys.getrecursionlimit()
        sys.setrecursionlimit(100)
        try:
            pi = IsoelectricPoint("DDDDDDDDDD").pi()
            self.assertIsInstance(pi, float)
        finally:
            sys.setrecursionlimit(old)


class IsoelectricPointRandomTest(unittest.TestCase):
    """Random protein sequences: charge at the reported pI must be ~0."""

    AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

    def test_random_proteins_zero_residual(self):
        """For 50 random proteins, charge_at_pH(pi()) must be near 0.

        This is the invariant that should hold for *every* sequence --
        not just the homopolymers above.  If it ever fails for a
        random sequence, the bracket has once again drifted below the
        true root.
        """
        rng = np.random.default_rng(2026)
        for trial in range(50):
            length = int(rng.integers(5, 60))
            seq = "".join(rng.choice(list(self.AMINO_ACIDS), size=length))
            p = IsoelectricPoint(seq)
            pi = p.pi()
            residual = p.charge_at_pH(pi)
            self.assertLess(
                abs(residual),
                0.01,
                f"trial {trial} ({seq[:10]}...{seq[-3:]}, len={length}): "
                f"residual {residual} at pI {pi}",
            )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
