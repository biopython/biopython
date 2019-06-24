# Copyright 2017 by Peter Cock.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Tests for the vector code in Bio.PDB."""

import unittest

try:
    import numpy
    del numpy
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.PDB.")

from Bio.PDB.vectors import Vector


class VectorTests(unittest.TestCase):
    """Tests for the Vector class."""

    def test_division(self):
        """Confirm division works."""
        v = Vector(1, 1, 1) / 2
        self.assertEqual(repr(v), "<Vector 0.50, 0.50, 0.50>")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
