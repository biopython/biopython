# Copyright 2020 by Rob Miller.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Unit tests for PARTS of the parse_pdb_header module of Bio.PDB."""


import unittest
import re
import warnings

try:
    import numpy  # noqa F401
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError("Install NumPy if you want to use Bio.PDB.")

from Bio.PDB.ic_rebuild import structure_rebuild_test, write_PDB
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser
from io import StringIO
from Bio.PDB.SCADIO import write_SCAD
from Bio.PDB.PICIO import write_PIC
from Bio.File import as_handle
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue
from Bio.PDB.internal_coords import IC_Residue

from Bio.PDB.PDBExceptions import PDBConstructionWarning


class Rebuild(unittest.TestCase):
    """Read PDB and mmCIF structures, convert to/from internal coordinates."""

    PDB_parser = PDBParser(PERMISSIVE=True, QUIET=True)
    CIF_parser = MMCIFParser(QUIET=True)
    pdb_1LCD = PDB_parser.get_structure("1LCD", "PDB/1LCD.pdb")
    pdb_2XHE = PDB_parser.get_structure("2XHE", "PDB/2XHE.pdb")
    cif_3JQH = CIF_parser.get_structure("3JQH", "PDB/3JQH.cif")
    cif_4CUP = CIF_parser.get_structure("4CUP", "PDB/4CUP.cif")

    def test_rebuild_multichain_missing(self):
        """Convert multichain missing atom protein to internal coordinates and back."""
        # 2XHE has regions of missing chain, last residue has only N
        r = structure_rebuild_test(self.pdb_2XHE, False)
        self.assertEqual(r["residues"], 787)
        self.assertEqual(r["rCount"], 835)
        self.assertEqual(r["rMatchCount"], 835)
        self.assertEqual(r["aCount"], 6267)
        self.assertEqual(r["disAtmCount"], 0)
        self.assertEqual(r["aCoordMatchCount"], 6267)
        self.assertEqual(len(r["chains"]), 2)
        self.assertTrue(r["pass"])

    def test_rebuild_disordered_atoms_residues(self):
        """Convert disordered protein to internal coordinates and back."""
        # 3jqh has both disordered residues
        # and disordered atoms in ordered residues
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", PDBConstructionWarning)
            r = structure_rebuild_test(self.cif_3JQH, False)
        # print(r)
        self.assertEqual(r["residues"], 26)
        self.assertEqual(r["rCount"], 47)
        self.assertEqual(r["rMatchCount"], 47)
        self.assertEqual(r["aCount"], 217)
        self.assertEqual(r["disAtmCount"], 50)
        self.assertEqual(r["aCoordMatchCount"], 217)
        self.assertEqual(len(r["chains"]), 1)
        self.assertTrue(r["pass"])

    def test_model_change_internal_coords(self):
        """Get model internal coords, modify psi and chi1 values and check."""
        for mdl in self.pdb_1LCD:
            if mdl.serial_num == 2:
                break
        mdl.atom_to_internal_coordinates()
        # other tests show can build with arbitrary internal coords
        # build here so changes below trigger more comlicated
        # xAtoms_needs_update mask arrays
        mdl.internal_to_atom_coordinates()
        nvt = {}
        nvc1 = {}
        nvpsi = {}
        tcount = 0
        c1count = 0
        psicount = 0
        for r in mdl.get_residues():
            ric = r.internal_coord
            if ric:
                # hedra change
                tau = ric.get_angle("tau")
                if ric.rprev != [] and tau is not None:
                    tcount += 1
                    nv = tau + 0.5
                    ric.set_angle("tau", nv)
                    nvt[str(r)] = nv
                # sidechain dihedron change
                chi1 = ric.get_angle("chi1")
                if chi1 is not None:
                    c1count += 1
                    nv = chi1 + 90
                    if nv > 180.0:
                        nv -= 360.0
                    ric.set_angle("chi1", nv)
                    nvc1[str(r)] = nv
                # backbone dihedron change
                psi = ric.get_angle("psi")
                if psi is not None:
                    psicount += 1
                    nv = psi - 90
                    if nv < -180.0:
                        nv += 360.0
                    ric.set_angle("psi", nv)
                    nvpsi[str(r)] = nv
        mdl.internal_to_atom_coordinates()
        sf = StringIO()
        write_PDB(self.pdb_1LCD, sf)
        sf.seek(0)
        new_1LCD = self.PDB_parser.get_structure("1LCD", sf)
        for mdl in new_1LCD:
            if mdl.serial_num == 2:
                break
        mdl.atom_to_internal_coordinates()
        ttcount = 0
        c1tcount = 0
        psitcount = 0
        for r in mdl.get_residues():
            ric = r.internal_coord
            if ric:
                tau = ric.get_angle("tau")
                if ric.rprev != [] and tau is not None:
                    ttcount += 1
                    self.assertAlmostEqual(tau, nvt[str(r)], places=1)
                chi1 = ric.get_angle("chi1")
                if chi1 is not None:
                    c1tcount += 1
                    self.assertAlmostEqual(chi1, nvc1[str(r)], places=1)
                psi = ric.get_angle("psi")
                if psi is not None:
                    psitcount += 1
                    self.assertAlmostEqual(psi, nvpsi[str(r)], places=1)
        self.assertEqual(tcount, ttcount)
        self.assertEqual(c1count, c1tcount)
        self.assertEqual(psicount, psitcount)
        self.assertTrue(ttcount > 0)
        self.assertTrue(c1count > 0)
        self.assertTrue(psicount > 0)

    def test_write_SCAD(self):
        """Check SCAD output plus MaxPeptideBond and Gly CB.

        SCAD tests: scaling, transform mtx, extra bond created (allBonds)
        """
        sf = StringIO()
        write_SCAD(
            self.cif_4CUP, sf, 10.0, pdbid="4cup", backboneOnly=True, includeCode=False
        )
        sf.seek(0)
        next_one = False
        with as_handle(sf, mode="r") as handle:
            for aline in handle.readlines():
                if "// (1856_S_CB, 1856_S_CA, 1856_S_C)" in aline:
                    m = re.search(r"\[\s+(\d+\.\d+)\,", aline)
                    if m:
                        # test correctly scaled atom bond length
                        self.assertAlmostEqual(float(m.group(1)), 15.30582, places=3)
                    else:
                        self.fail("scaled atom bond length not found")
                elif '[ 1, "1857M",' in aline:
                    next_one = True
                elif next_one:
                    next_one = False
                    # test last residue transform looks roughly correct
                    # some differences due to sorting issues on different python
                    # versions
                    target = [-12.413, -3.303, 35.771, 1.0]
                    ms = re.findall(  # last column of each row
                        r"\s+(-?\d+\.\d+)\s+\]", aline
                    )
                    if ms:
                        for i in range(0, 3):
                            self.assertAlmostEqual(float(ms[i]), target[i], places=0)
                    else:
                        self.fail("transform not found")
        sf.seek(0)
        IC_Residue.gly_Cbeta = True
        write_SCAD(
            self.pdb_2XHE[0]["A"],
            sf,
            10.0,
            pdbid="2xhe",
            maxPeptideBond=100.0,
            includeCode=False,
        )
        sf.seek(0)
        allBondsPass = False
        maxPeptideBondPass = False
        glyCbetaFound = False
        with as_handle(sf, mode="r") as handle:
            for aline in handle.readlines():
                # test extra bond created in TRP (allBonds is True)
                if '"Cres", 0, 0, 1, 0, StdBond, "W", 24, "CD2CE3CZ3"' in aline:
                    allBondsPass = True
                # test 509_K-561_E long bond created
                if "509_K" in aline and "561_E" in aline:
                    maxPeptideBondPass = True
                if "(21_G_CB, 21_G_CA, 21_G_C)" in aline:
                    glyCbetaFound = True
                    target = [15.33630, 110.17513, 15.13861]
                    ms = re.findall(r"\s+(-?\d+\.\d+)", aline)
                    if ms:
                        for i in range(0, 3):
                            self.assertAlmostEqual(float(ms[i]), target[i], places=0)
                    else:
                        self.fail("Cbeta internal coords not found")

        self.assertTrue(allBondsPass)
        self.assertTrue(glyCbetaFound)
        self.assertTrue(maxPeptideBondPass)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
