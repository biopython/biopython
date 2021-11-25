# Copyright 2020-2021 by Robert T. Miller.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Unit tests for internal_coords module of Bio.PDB."""


import unittest
import re
import warnings
import copy

try:
    import numpy as np  # noqa F401
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError("Install NumPy if you want to use Bio.PDB.")

from Bio.PDB.ic_rebuild import (
    structure_rebuild_test,
    IC_duplicate,
    compare_residues,
)
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser
from io import StringIO
from Bio.PDB.SCADIO import write_SCAD
from Bio.PDB.PICIO import read_PIC_seq

from Bio.File import as_handle

from Bio.PDB.internal_coords import IC_Residue, IC_Chain, Dihedron

from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio import SeqIO


class Rebuild(unittest.TestCase):
    """Read PDB and mmCIF structures, convert to/from internal coordinates."""

    PDB_parser = PDBParser(PERMISSIVE=True, QUIET=True)
    CIF_parser = MMCIFParser(QUIET=True)
    pdb_1LCD = PDB_parser.get_structure("1LCD", "PDB/1LCD.pdb")
    # cif_1A7G = CIF_parser.get_structure("1A7G", "PDB/1A7G.cif")
    # cif_1A7G2 = CIF_parser.get_structure("1A7G", "PDB/1A7G.cif")
    pdb_2XHE = PDB_parser.get_structure("2XHE", "PDB/2XHE.pdb")
    pdb_2XHE2 = PDB_parser.get_structure("2XHE", "PDB/2XHE.pdb")
    cif_3JQH = CIF_parser.get_structure("3JQH", "PDB/3JQH.cif")
    cif_4CUP = CIF_parser.get_structure("4CUP", "PDB/4CUP.cif")
    cif_4CUP2 = CIF_parser.get_structure("4CUP", "PDB/4CUP.cif")
    cif_4ZHL = CIF_parser.get_structure("4ZHL", "PDB/4ZHL.cif")
    cif_4ZHL2 = CIF_parser.get_structure("4ZHL", "PDB/4ZHL.cif")

    def test_rebuild_multichain_missing(self):
        """Convert multichain missing atom struct to, from internal coords."""
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
        with warnings.catch_warnings(record=True):
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

    def test_no_crosstalk(self):
        """Deep copy, change few internal coords, test nothing else changes."""
        # IC_Chain.ParallelAssembleResidues = False
        self.cif_4CUP.atom_to_internal_coordinates()
        cpy4cup = copy.deepcopy(self.cif_4CUP)
        cic0 = self.cif_4CUP.child_list[0].child_list[0].internal_coord
        cic1 = cpy4cup.child_list[0].child_list[0].internal_coord
        alist = [
            "omg",
            "phi",
            "psi",
            "chi1",
            "chi2",
            "chi3",
            "chi4",
            "chi5",
            "tau",
        ]
        delta = 33  # degrees to change
        tdelta = delta / 10.0  # more realistic for bond angle
        targPos = 1
        for ang in alist:
            # skip by 2's with alist along original chain changing angle spec
            ricTarg = cic0.chain.child_list[targPos].internal_coord
            # print(targPos + 1, ricTarg.lc, ang)
            targPos += 2
            try:
                andx = ricTarg.pick_angle(ang).ndx
                if ang == "tau":
                    cic0.hedraAngle[andx] += tdelta
                    cic0.hAtoms_needs_update[andx] = True
                    cic0.atomArrayValid[cic0.h2aa[andx]] = False
                    cic0.hAtoms_needs_update[:] = True
                    cic0.atomArrayValid[:] = False
                    cic0.dAtoms_needs_update[:] = True
                else:
                    cic0.dihedraAngle[andx] += delta
                    if cic0.dihedraAngle[andx] > 180.0:
                        cic0.dihedraAngle[andx] -= 360.0
                    cic0.dihedraAngleRads[andx] = np.deg2rad(cic0.dihedraAngle[andx])
                    cic0.dAtoms_needs_update[andx] = True
                    cic0.atomArrayValid[cic0.d2aa[andx]] = False
            except AttributeError:
                pass  # skip if residue does not have e.g. chi5
        cic0.internal_to_atom_coordinates()  # move atoms
        cic0.atom_to_internal_coordinates()  # get new internal coords
        # generate hdelta and ddelta difference arrays so can look for what
        # changed
        hdelta = cic0.hedraAngle - cic1.hedraAngle
        hdelta[np.abs(hdelta) < 0.00001] = 0.0
        ddelta = cic0.dihedraAngle - cic1.dihedraAngle
        ddelta[np.abs(ddelta) < 0.00001] = 0.0
        ddelta[ddelta < -180.0] += 360.0  # wrap around circle values
        targPos = 1
        for ang in alist:
            # same skip along original chain looking at hdelta and ddelta
            # if change is as specified, set difference to 0 then we can test
            # for any remaining (spurious) changes
            ricTarg = cic0.chain.child_list[targPos].internal_coord
            # print(targPos + 1, ricTarg.lc, ang)
            targPos += 2
            try:
                andx = ricTarg.pick_angle(ang).ndx
                if ang == "tau":
                    self.assertAlmostEqual(hdelta[andx], tdelta, places=4)
                    hdelta[andx] = 0.0
                    # some other angle has to change to accommodate tau change
                    # N-Ca-Cb is artifact of choices in ic_data
                    # expected change so clear relevant hdelta here
                    adjAngNdx = ricTarg.pick_angle("N:CA:CB").ndx
                    self.assertNotAlmostEqual(hdelta[adjAngNdx], 0.0, places=1)
                    hdelta[adjAngNdx] = 0.0
                else:
                    self.assertAlmostEqual(ddelta[andx], delta, places=4)
                    ddelta[andx] = 0.0
            except AttributeError:
                pass  # if residue does not have e.g. chi5

        hsum = hdelta.sum()
        self.assertEqual(hsum, 0.0)
        dsum = ddelta.sum()
        self.assertEqual(dsum, 0.0)

    def test_model_change_internal_coords(self):
        """Get model internal coords, modify psi and chi1 values and check."""
        mdl = self.pdb_1LCD[1]
        mdl.atom_to_internal_coordinates()
        # other tests show can build with arbitrary internal coords
        # build here so changes below trigger more complicated
        # Atoms_needs_update mask arrays
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

        # prove not using stored results
        for chn in mdl.get_chains():
            if hasattr(chn, "hedraLen"):
                delattr(chn.internal_coord, "hedraLen")
                delattr(chn.internal_coord, "dihedraLen")
                delattr(chn.internal_coord, "hedraAngle")
                delattr(chn.internal_coord, "dihedraAngle")
                for r in chn.get_residues():
                    r.internal_coord.hedra = {}
                    r.internal_coord.dihedra = {}

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
                    # print(str(r), "tau", tau, nvt[str(r)])
                    self.assertAlmostEqual(tau, nvt[str(r)], places=3)
                chi1 = ric.get_angle("chi1")
                if chi1 is not None:
                    c1tcount += 1
                    # print(str(r), "chi1", chi1, nvc1[str(r)])
                    self.assertAlmostEqual(chi1, nvc1[str(r)], places=3)
                psi = ric.get_angle("psi")
                if psi is not None:
                    psitcount += 1
                    # print(str(r), "psi", psi, nvpsi[str(r)])
                    self.assertAlmostEqual(psi, nvpsi[str(r)], places=3)
        self.assertEqual(tcount, ttcount)
        self.assertEqual(c1count, c1tcount)
        self.assertEqual(psicount, psitcount)
        self.assertGreater(ttcount, 0)
        self.assertGreater(c1count, 0)
        self.assertGreater(psicount, 0)

    def test_write_SCAD(self):
        """Check SCAD output plus MaxPeptideBond and Gly CB.

        SCAD tests: scaling, transform mtx, extra bond created (allBonds)
        """
        sf = StringIO()
        write_SCAD(
            self.cif_4CUP2,
            sf,
            10.0,
            pdbid="4cup",
            backboneOnly=True,
            includeCode=False,
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
                    # some differences due to sorting issues on different
                    # python versions
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
            self.pdb_2XHE2[0]["A"],
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

    def test_i2a_start_fin(self):
        """Test assemble start/fin, default NCaC coordinates, IC_duplicate."""
        chn = self.pdb_1LCD[2]["A"]
        cpy = IC_duplicate(chn)[2]["A"]  # generates internal coords as needed
        cpy.internal_to_atom_coordinates(start=31, fin=45)
        cdict = compare_residues(chn, cpy, quick=True)
        self.assertFalse(cdict["pass"])
        # transform source coordinates to put res 31 tau at origin like
        # fragment
        res = chn[31]
        psi = res.internal_coord.pick_angle("psi")
        cst = np.transpose(psi.cst)
        chn.internal_coord.atomArray[:] = chn.internal_coord.atomArray.dot(cst)
        cdict = compare_residues(chn, cpy, rtol=1e-03, atol=1e-05)
        self.assertEqual(cdict["residues"], 51)
        self.assertEqual(cdict["rMatchCount"], 77)
        self.assertEqual(cdict["aCount"], 497)
        self.assertEqual(cdict["disAtmCount"], 0)
        self.assertEqual(cdict["aCoordMatchCount"], 140)
        self.assertEqual(cdict["aFullIdMatchCount"], 140)
        self.assertEqual(len(cdict["chains"]), 1)
        self.assertEqual(cdict["rCount"], 77)
        self.assertFalse(cdict["pass"])

    def test_distplot_rebuild(self):
        """Build identical structure from distplot and chirality data."""
        # load input chain
        for _chn1 in self.cif_4ZHL.get_chains():
            break
        # create atomArray and compute distplot and dihedral signs array
        _chn1.atom_to_internal_coordinates()
        dplot1 = _chn1.internal_coord.distance_plot()
        dsigns = _chn1.internal_coord.dihedral_signs()

        # load second copy (same again) input chain
        for _chn2 in self.cif_4ZHL2.get_chains():
            break
        # create internal coord structures but do not compute di/hedra
        cic2 = _chn2.internal_coord = IC_Chain(_chn2)
        cic2.init_edra()
        # load relevant interatomic distances from chn1 distance plot
        cic2.distplot_to_dh_arrays(dplot1)
        # compute di/hedra angles from dh_arrays
        cic2.distance_to_internal_coordinates(dsigns)

        # clear chn2 atom coordinates
        cic2.atomArrayValid[:] = False
        # initialize values but this is redundant to Valid=False above
        cic2.atomArray = np.zeros((cic2.AAsiz, 4), dtype=np.float64)
        cic2.atomArray[:, 3] = 1.0

        # 4zhl has chain breaks so copy initial coords of each segment
        cic2.copy_initNCaCs(_chn1.internal_coord)
        # compute chn2 atom coords from di/hedra data
        cic2.internal_to_atom_coordinates()

        # generate distance plot from second chain, confirm minimal distance
        # from original
        dp2 = cic2.distance_plot()
        dpdiff = np.abs(dplot1 - dp2)
        # print(np.amax(dpdiff))
        self.assertTrue(np.amax(dpdiff) < 0.000001)

    def test_seq_as_PIC(self):
        """Read protein sequence, generate default PIC data."""
        seqIter = SeqIO.parse("Fasta/f001", "fasta")
        for _record in seqIter:
            break
        pdb_structure = read_PIC_seq(_record)
        pdb_structure.internal_to_atom_coordinates()
        for _chn in pdb_structure.get_chains():
            break
        assert len(_chn.internal_coord.atomArrayValid) == 575

    def test_angle_fns(self):
        """Test angle_dif and angle_avg across +/-180 boundaries."""
        arr1 = np.array([179.0, 90.0, 88.0, 1.0])
        arr2 = np.array([-179.0, -90.0, -91.0, -1.0])
        assert (
            Dihedron.angle_dif(arr1, arr2) == np.array([2.0, 180.0, -179.0, -2.0])
        ).all()
        assert Dihedron.angle_avg(np.array([179.0, -179.0])) == 180.0
        assert Dihedron.angle_avg(np.array([1.0, -1.0])) == 0.0
        assert Dihedron.angle_avg(np.array([90.0, -90.0])) == 0.0
        assert Dihedron.angle_avg(np.array([91.0, -91.0])) == 180.0


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
