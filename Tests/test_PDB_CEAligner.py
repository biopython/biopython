# Copyright 2022 by Joao Rodrigues. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Unit tests for the Bio.PDB.CEAligner module."""

import unittest

try:
    import numpy as np
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.PDB."
    ) from None

from Bio.PDB import CEAligner
from Bio.PDB import MMCIFParser


class CEAlignerTests(unittest.TestCase):
    """Test CEAligner class."""

    @staticmethod
    def _get_ca_coords_as_array(structure):
        xyz_list = [a.coord for a in structure.get_atoms() if a.name == "CA"]
        return np.asarray(xyz_list, dtype=np.float64)

    def test_cealigner(self):
        """Test aligning 7CFN on 6WQA."""
        ref = "PDB/6WQA.cif"
        mob = "PDB/7CFN.cif"
        result = "PDB/7CFN_aligned.cif"

        parser = MMCIFParser(QUIET=1)
        s1 = parser.get_structure("6wqa", ref)
        s2 = parser.get_structure("7cfn", mob)

        aligner = CEAligner()
        aligner.set_reference(s1)
        aligner.align(s2)

        self.assertAlmostEqual(aligner.rms, 3.83, places=2)

        # Assert the transformation was done right by comparing
        # the moved coordinates to a 'ground truth' reference.
        # Reference obtained with Pymol's CEAlign code.
        refe = parser.get_structure("7cfn_aligned", result)
        refe_coords = self._get_ca_coords_as_array(refe)
        s2_f_coords = self._get_ca_coords_as_array(s2)

        diff = refe_coords - s2_f_coords
        rmsd = np.sqrt((diff * diff).sum() / len(refe_coords))
        self.assertAlmostEqual(rmsd, 0.0, places=2)  # < 0.01

    def test_cealigner_no_transform(self):
        """Test aligning 7CFN on 6WQA without transforming 7CFN."""
        ref = "PDB/6WQA.cif"
        mob = "PDB/7CFN.cif"

        parser = MMCIFParser(QUIET=1)
        s1 = parser.get_structure("6wqa", ref)
        s2 = parser.get_structure("7cfn", mob)

        s2_original_coords = [list(a.coord) for a in s2.get_atoms()]

        aligner = CEAligner()
        aligner.set_reference(s1)
        aligner.align(s2, transform=False)
        s2_coords_final = [list(a.coord) for a in s2.get_atoms()]

        self.assertAlmostEqual(aligner.rms, 3.83, places=2)
        self.assertEqual(s2_original_coords, s2_coords_final)

    def test_cealigner_nucleic(self):
        """Test aligning 1LCD on 1LCD."""
        ref = "PDB/1LCD.cif"
        mob = "PDB/1LCD.cif"

        parser = MMCIFParser(QUIET=1)
        s1 = parser.get_structure("1lcd_ref", ref)
        s2 = parser.get_structure("1lcd_mob", mob)

        aligner = CEAligner()
        aligner.set_reference(s1)
        aligner.align(s2)

        self.assertAlmostEqual(aligner.rms, 0.0, places=3)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
