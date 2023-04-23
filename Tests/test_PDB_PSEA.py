# Copyright 2022 by Valentin Vareskic (valentin.vareskic@gmail.com).
# All rights reserved. This code is part of the Biopython distribution
# and governed by its license.  Please see the LICENSE file that should
# have been included as part of this package.

"""Tests for PDB PSEA."""

import io
import os
import unittest
from subprocess import getoutput
import sys

from Bio import MissingExternalDependencyError
from Bio.PDB import PDBParser
from Bio.PDB.PSEA import run_psea, psea, psea2HEC, PSEA


os.environ["LANG"] = "C"


cmd_output = getoutput("psea -h")
if not cmd_output.startswith("o---"):
    raise MissingExternalDependencyError(
        "Download and install psea from "
        "ftp://ftp.lmcp.jussieu.fr/pub/sincris/software/protein/p-sea/. "
        "Make sure that psea is on path"
    )


def remove_sea_files():
    for file in os.listdir():
        if file.endswith(".sea"):
            os.remove(file)


class TestPDBPSEA(unittest.TestCase):
    def tearDown(self):
        remove_sea_files()

    def test_run_psea_verbose(self):
        captured_ouput = io.StringIO()
        sys.stdout = captured_ouput
        psae_run = run_psea("PDB/1A8O.pdb", verbose=True)
        sys.stdout = sys.__stdout__
        self.assertEqual(psae_run, "1A8O.sea")
        self.assertTrue(captured_ouput.getvalue())

    def test_run_psea_quiet(self):
        captured_ouput = io.StringIO()
        sys.stdout = captured_ouput
        psae_run = run_psea("PDB/1A8O.pdb", verbose=False)
        sys.stdout = sys.__stdout__
        self.assertEqual(psae_run, "1A8O.sea")
        self.assertFalse(captured_ouput.getvalue())

    def test_psea(self):
        psae_run = psea("PDB/2BEG.pdb")
        self.assertEqual(psae_run, "ccccbbbbbbbccccbbbbbbbbbbc")

    def test_psea_2HEC(self):
        seq = psea("PDB/2BEG.pdb")
        psae_run = psea2HEC(seq)
        self.assertEqual(
            psae_run,
            [
                "C",
                "C",
                "C",
                "C",
                "E",
                "E",
                "E",
                "E",
                "E",
                "E",
                "E",
                "C",
                "C",
                "C",
                "C",
                "E",
                "E",
                "E",
                "E",
                "E",
                "E",
                "E",
                "E",
                "E",
                "E",
                "C",
            ],
        )


class TestPSEA(unittest.TestCase):
    def tearDown(self):
        remove_sea_files()

    def test_get_seq(self):
        p = PDBParser()
        s = p.get_structure("X", "PDB/2BEG.pdb")
        psea_class = PSEA(s[0], "PDB/2BEG.pdb")
        self.assertEqual(
            psea_class.get_seq(),
            [
                "C",
                "C",
                "C",
                "C",
                "E",
                "E",
                "E",
                "E",
                "E",
                "E",
                "E",
                "C",
                "C",
                "C",
                "C",
                "E",
                "E",
                "E",
                "E",
                "E",
                "E",
                "E",
                "E",
                "E",
                "E",
                "C",
            ],
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
