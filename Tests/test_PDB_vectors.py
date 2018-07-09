# Copyright 2017 by Peter Cock.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Tests for the vector code in Bio.PDB."""

import unittest
import warnings

try:
    import numpy
    del numpy
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.PDB.")

from Bio import BiopythonDeprecationWarning


class LegacyImportTests(unittest.TestCase):
    """Tests for the legacy Bio.PDB.Vector module."""

    def test_legacy_imports(self):
        """Confirm legacy (deprecated) vector imports still work."""
        # These now live under Bio.PDB.vectors instead.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonDeprecationWarning)
            from Bio.PDB.Vector import calc_angle, calc_dihedral
            from Bio.PDB.Vector import m2rotaxis, refmat, rotaxis2m, rotmat
            from Bio.PDB.Vector import vector_to_axis, Vector
            del calc_angle, calc_dihedral
            del m2rotaxis, refmat, rotaxis2m, rotmat
            del vector_to_axis, Vector


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
