# Copyright 2022 by Valentin Vareskic (valentin.vareskic@gmail.com).
# All rights reserved. This code is part of the Biopython distribution
# and governed by its license.  Please see the LICENSE file that should
# have been included as part of this package.

"""Tests for PDB PSEA."""

import unittest
from Bio.PDB.PSEA import run_psea, psea, psea2HEC, PSEA
from Bio.PDB import PDBParser
from subprocess import getoutput
from Bio import MissingExternalDependencyError
import os

if "command not found" in getoutput("psea -h"):
    raise MissingExternalDependencyError(
        "Download and install psea from ftp://ftp.lmcp.jussieu.fr/pub/sincris/software/protein/p-sea/. Make sure that psea is on path"
    )


class TestPDBPSEA(unittest.TestCase):
    def test_run_psea(self):
        psae_run = run_psea("PDB/1A8O.pdb")
        self.assertEqual(psae_run, "1A8O.sea")

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


def clean_up():
    file_list: list = ["1A8O", "2BEG"]
    try:
        [os.remove(f"./{file}.sea") for file in file_list]
    except OSError as err:
        print("No eligible files found for deletion.")
