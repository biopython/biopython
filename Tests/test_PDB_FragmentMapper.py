# Copyright 2017 by Sourav Singh. All rights reserved.
# Revisions copyright 2017 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Tests for PDB FragmentMapper module."""

import unittest

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.PDB."
    ) from None

from Bio.PDB import PDBParser
from Bio.PDB import FragmentMapper
from Bio.PDB import Selection


class FragmentMapperTests(unittest.TestCase):
    """Tests for FragmentMapper module."""

    def test_fragment_mapper(self):
        """Self test for FragmentMapper module."""
        p = PDBParser()
        pdb1 = "PDB/1A8O.pdb"
        s = p.get_structure("X", pdb1)
        m = s[0]
        fm = FragmentMapper(m, 10, 5, "PDB")
        for r in Selection.unfold_entities(m, "R"):
            if r in fm:
                self.assertTrue(str(fm[r]).startswith("<Fragment length=5 id="))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
