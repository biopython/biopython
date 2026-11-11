# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Unit tests for the Scripts/PDB/hsexpo.py example script."""

import os
import subprocess
import sys
import unittest
import warnings

from Bio.PDB import HSExposureCA
from Bio.PDB import HSExposureCB
from Bio.PDB import PDBParser
from Bio.PDB import Selection
from Bio.PDB.PDBExceptions import PDBConstructionWarning

SCRIPT = os.path.join(os.path.dirname(__file__), "..", "Scripts", "PDB", "hsexpo.py")
PDB_FILENAME = os.path.join("PDB", "a_structure.pdb")
RADIUS = 13.0  # default radius used by the script


def run_script(exp_type):
    """Run hsexpo.py for one exposure type, return its printed values."""
    out = subprocess.run(
        [sys.executable, SCRIPT, PDB_FILENAME, "-t", exp_type, "-m", "1"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    values = []
    for line in out.splitlines():
        resname, resseq, icode, value = line.split()
        values.append((resname, int(resseq), icode, int(value)))
    return values


class HSExpoScript(unittest.TestCase):
    """Check the script prints the exposure values its options ask for."""

    def setUp(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            structure = PDBParser().get_structure("X", PDB_FILENAME)
        # Model 0 of this fixture is a pathological parser test case with
        # no complete peptides; model 1 is the one used for exposure tests
        # (see test_PDB_Exposure.py).
        self.model = structure[1]

    def library_values(self, hse_class, key):
        """Compute the values the script is expected to print."""
        hse_class(self.model, RADIUS)
        values = []
        for r in Selection.unfold_entities(self.model, "R"):
            if key in r.xtra:
                hetflag, resseq, icode = r.get_id()
                if icode == " ":
                    icode = "_"  # the script prints blank icodes as underscores
                values.append((r.get_resname(), resseq, icode, r.xtra[key]))
        return values

    def check(self, exp_type, hse_class, key):
        self.assertEqual(run_script(exp_type), self.library_values(hse_class, key))

    def test_hsebeta_down(self):
        """-t HSEBD prints the down (EXP_HSE_B_D) values."""
        self.check("HSEBD", HSExposureCB, "EXP_HSE_B_D")

    def test_hsebeta_up(self):
        """-t HSEBU prints the up (EXP_HSE_B_U) values."""
        self.check("HSEBU", HSExposureCB, "EXP_HSE_B_U")

    def test_hsealpha_down(self):
        """-t HSEAD prints the down (EXP_HSE_A_D) values."""
        self.check("HSEAD", HSExposureCA, "EXP_HSE_A_D")

    def test_hsealpha_up(self):
        """-t HSEAU prints the up (EXP_HSE_A_U) values."""
        self.check("HSEAU", HSExposureCA, "EXP_HSE_A_U")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
