# Copyright 2009-2011 by Eric Talevich.  All rights reserved.
# Revisions copyright 2009-2013 by Peter Cock.  All rights reserved.
# Revisions copyright 2013 Lenna X. Peterson. All rights reserved.
#
# Converted by Eric Talevich from an older unit test copyright 2002
# by Thomas Hamelryck.
#
# Merged related test files into one, by Joao Rodrigues (2020)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Unit tests for the Bio.PDB exposure classes."""

from copy import deepcopy
import os
import sys
import tempfile
import unittest
import warnings
from io import StringIO

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.PDB."
    ) from None

from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB import HSExposureCA, HSExposureCB, ExposureCN
from Bio.PDB.PDBExceptions import PDBConstructionException, PDBConstructionWarning


class Exposure(unittest.TestCase):
    """Testing Bio.PDB.HSExposure."""

    def setUp(self):
        pdb_filename = "PDB/a_structure.pdb"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            structure = PDBParser(PERMISSIVE=True).get_structure("X", pdb_filename)
        structure[1].detach_child("B")
        self.model = structure[1]
        # Look at first chain only
        a_residues = list(self.model["A"].child_list)
        self.assertEqual(86, len(a_residues))
        self.assertEqual(a_residues[0].get_resname(), "CYS")
        self.assertEqual(a_residues[1].get_resname(), "ARG")
        self.assertEqual(a_residues[2].get_resname(), "CYS")
        self.assertEqual(a_residues[3].get_resname(), "GLY")
        # ...
        self.assertEqual(a_residues[-3].get_resname(), "TYR")
        self.assertEqual(a_residues[-2].get_resname(), "ARG")
        self.assertEqual(a_residues[-1].get_resname(), "CYS")
        self.a_residues = a_residues
        self.radius = 13.0

    def test_HSExposureCA(self):
        """HSExposureCA."""
        _ = HSExposureCA(self.model, self.radius)
        residues = self.a_residues
        self.assertEqual(0, len(residues[0].xtra))
        self.assertEqual(0, len(residues[1].xtra))
        self.assertEqual(3, len(residues[2].xtra))
        self.assertAlmostEqual(
            0.81250973133184456, residues[2].xtra["EXP_CB_PCB_ANGLE"]
        )
        self.assertEqual(14, residues[2].xtra["EXP_HSE_A_D"])
        self.assertEqual(14, residues[2].xtra["EXP_HSE_A_U"])
        self.assertEqual(3, len(residues[3].xtra))
        self.assertAlmostEqual(1.3383737, residues[3].xtra["EXP_CB_PCB_ANGLE"])
        self.assertEqual(13, residues[3].xtra["EXP_HSE_A_D"])
        self.assertEqual(16, residues[3].xtra["EXP_HSE_A_U"])
        # ...
        self.assertEqual(3, len(residues[-2].xtra))
        self.assertAlmostEqual(
            0.77124014456278489, residues[-2].xtra["EXP_CB_PCB_ANGLE"]
        )
        self.assertEqual(24, residues[-2].xtra["EXP_HSE_A_D"])
        self.assertEqual(24, residues[-2].xtra["EXP_HSE_A_U"])
        self.assertEqual(0, len(residues[-1].xtra))

    def test_HSExposureCB(self):
        """HSExposureCB."""
        _ = HSExposureCB(self.model, self.radius)
        residues = self.a_residues
        self.assertEqual(0, len(residues[0].xtra))
        self.assertEqual(2, len(residues[1].xtra))
        self.assertEqual(20, residues[1].xtra["EXP_HSE_B_D"])
        self.assertEqual(5, residues[1].xtra["EXP_HSE_B_U"])
        self.assertEqual(2, len(residues[2].xtra))
        self.assertEqual(10, residues[2].xtra["EXP_HSE_B_D"])
        self.assertEqual(18, residues[2].xtra["EXP_HSE_B_U"])
        self.assertEqual(2, len(residues[3].xtra))
        self.assertEqual(7, residues[3].xtra["EXP_HSE_B_D"])
        self.assertEqual(22, residues[3].xtra["EXP_HSE_B_U"])
        # ...
        self.assertEqual(2, len(residues[-2].xtra))
        self.assertEqual(14, residues[-2].xtra["EXP_HSE_B_D"])
        self.assertEqual(34, residues[-2].xtra["EXP_HSE_B_U"])
        self.assertEqual(2, len(residues[-1].xtra))
        self.assertEqual(23, residues[-1].xtra["EXP_HSE_B_D"])
        self.assertEqual(15, residues[-1].xtra["EXP_HSE_B_U"])

    def test_ExposureCN(self):
        """HSExposureCN."""
        _ = ExposureCN(self.model, self.radius)
        residues = self.a_residues
        self.assertEqual(0, len(residues[0].xtra))
        self.assertEqual(1, len(residues[1].xtra))
        self.assertEqual(25, residues[1].xtra["EXP_CN"])
        self.assertEqual(1, len(residues[2].xtra))
        self.assertEqual(28, residues[2].xtra["EXP_CN"])
        self.assertEqual(1, len(residues[3].xtra))
        self.assertEqual(29, residues[3].xtra["EXP_CN"])
        # ...
        self.assertEqual(1, len(residues[-2].xtra))
        self.assertEqual(48, residues[-2].xtra["EXP_CN"])
        self.assertEqual(1, len(residues[-1].xtra))
        self.assertEqual(38, residues[-1].xtra["EXP_CN"])


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
