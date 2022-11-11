# Copyright 2009-2011 by Eric Talevich.  All rights reserved.
# Revisions copyright 2009-2013 by Peter Cock.  All rights reserved.
# Revisions copyright 2013 Lenna X. Peterson. All rights reserved.
# Revisions copyright 2020 Joao Rodrigues. All rights reserved.
#
# Converted by Eric Talevich from an older unit test copyright 2002
# by Thomas Hamelryck.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Unit tests for the Bio.PDB PDBParser module."""

from io import StringIO
import os
import tempfile
import unittest
import warnings

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.PDB."
    ) from None

from Bio import BiopythonWarning
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionException, PDBConstructionWarning


class FlawedPDB_tests(unittest.TestCase):
    """Errors and warnings while parsing flawed PDB files.

    These tests must be executed because of the way Python's warnings module
    works -- a warning is only logged the first time it is encountered.
    """

    def setUp(self):
        self.permissive = PDBParser(PERMISSIVE=True)
        self.strict = PDBParser(PERMISSIVE=False)

    def test_1_flawedpdb_permissive(self):
        """Parse a flawed PDB file in permissive mode: check warnings."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", PDBConstructionWarning)

            # Trigger warnings
            self.permissive.get_structure("example", "PDB/a_structure.pdb")

            self.assertEqual(len(w), 15)
            for wrn, msg in zip(
                w,
                [
                    # Expected warning messages:
                    "Used element 'N' for Atom (name=N) with given element ''",
                    "Used element 'C' for Atom (name=CA) with given element ''",
                    "Atom names ' CA ' and 'CA  ' differ only in spaces at line 18.",
                    "Used element 'CA' for Atom (name=CA  ) with given element ''",
                    "Atom N defined twice in residue <Residue ARG het=  resseq=2 icode= > at line 22.",
                    "disordered atom found with blank altloc before line 34.",
                    "Residue (' ', 4, ' ') redefined at line 44.",
                    "Blank altlocs in duplicate residue SER (' ', 4, ' ') at line 44.",
                    "Residue (' ', 10, ' ') redefined at line 76.",
                    "Residue (' ', 14, ' ') redefined at line 107.",
                    "Residue (' ', 16, ' ') redefined at line 136.",
                    "Residue (' ', 80, ' ') redefined at line 634.",
                    "Residue (' ', 81, ' ') redefined at line 647.",
                    "Ignoring unrecognized record 'ATOM 1' at line 777",
                    "Atom O defined twice in residue <Residue HOH het=W resseq=67 icode= > at line 904.",
                ],
            ):
                self.assertIn(msg, str(wrn))

    def test_2_flawedpdb_strict(self):
        """Parse a flawed PDB file in permissive mode: check errors."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", PDBConstructionWarning)
            self.assertRaises(
                PDBConstructionException,
                self.strict.get_structure,
                "example",
                "PDB/a_structure.pdb",
            )

            self.assertEqual(len(w), 4, w)

    def test_3_bad_xyz_permissive(self):
        """Parse an entry with bad x,y,z value with PERMISSIVE=True."""
        data = "ATOM      9  N   ASP A 152      21.554  34.953  27.691  1.00 19.26           N\n"
        _ = self.permissive.get_structure("example", StringIO(data))

    def test_4_bad_xyz_strict(self):
        """Parse an entry with bad x,y,z value with PERMISSIVE=False."""
        data = "ATOM      9  N   ASP A 152      21.ish  34.953  27.691  1.00 19.26           N\n"
        with self.assertRaises(PDBConstructionException):
            self.strict.get_structure("example", StringIO(data))

    def test_5_missing_occupancy_permissive(self):
        """Parse file with missing occupancy with PERMISSIVE=True."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", PDBConstructionWarning)
            structure = self.permissive.get_structure("test", "PDB/occupancy.pdb")
            self.assertEqual(len(w), 3, w)

        atoms = structure[0]["A"][(" ", 152, " ")]

        # Blank occupancy behavior set in Bio/PDB/PDBParser
        self.assertIsNone(atoms["N"].get_occupancy())
        self.assertEqual(atoms["CA"].get_occupancy(), 1.0)
        self.assertEqual(atoms["C"].get_occupancy(), 0.0)

    def test_6_missing_occupancy_strict(self):
        """Parse file with missing occupancy with PERMISSIVE=False."""
        with self.assertRaises(PDBConstructionException):
            _ = self.strict.get_structure("test", "PDB/occupancy.pdb")


class ParseDummyPDB_test(unittest.TestCase):
    """Tests for artificial/dummy PDB files."""

    @classmethod
    def setUpClass(cls):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            p = PDBParser(PERMISSIVE=1)
            cls.structure = p.get_structure("example", "PDB/a_structure.pdb")

    def test_structure_integrity(self):
        """Verify the structure of the parsed example PDB file."""
        # Structure contains 2 models
        self.assertEqual(len(self.structure), 2)
        # --- Checking model 0 ---
        m0 = self.structure[0]
        # Model 0 contains 1 chain
        self.assertEqual(len(m0), 1)
        # Chain 'A' contains 1 residue
        self.assertEqual(len(m0["A"]), 1)
        # Residue ('H_PCA', 1, ' ') contains 9 atoms.
        residue = m0["A"].get_list()[0]
        self.assertEqual(residue.get_id(), ("H_PCA", 1, " "))
        self.assertEqual(len(residue), 9)
        # --- Checking model 1 ---
        m1 = self.structure[1]
        # Model 1 contains 3 chains
        self.assertEqual(len(m1), 4)
        # Deconstruct this data structure to check each chain
        chain_data = [  # chain_id, chain_len, [(residue_id, residue_len), ...]
            (
                "A",
                86,
                [
                    ((" ", 0, " "), 1),
                    ((" ", 2, " "), 11),
                    ((" ", 3, " "), 6, 1),  # disordered
                    ((" ", 4, " "), 4),
                    ((" ", 5, " "), 6),
                    ((" ", 6, " "), 9),
                    ((" ", 7, " "), 4),
                    ((" ", 8, " "), 4),
                    ((" ", 9, " "), 4),
                    ((" ", 10, " "), 6, ["GLY", "SER"]),  # point mut
                    ((" ", 11, " "), 7),
                    ((" ", 12, " "), 6),
                    ((" ", 13, " "), 7),
                    ((" ", 14, " "), 4, ["ALA", "GLY"]),  # point mut
                    ((" ", 15, " "), 8, 3),  # disordered
                    ((" ", 16, " "), 11, ["ARG", "TRP"]),  # point mut
                    ((" ", 17, " "), 6),
                    ((" ", 18, " "), 6),
                    ((" ", 19, " "), 6),
                    ((" ", 20, " "), 8),
                    ((" ", 21, " "), 14),
                    ((" ", 22, " "), 4),
                    ((" ", 23, " "), 14),
                    ((" ", 24, " "), 6),
                    ((" ", 25, " "), 4),
                    ((" ", 26, " "), 8),
                    ((" ", 27, " "), 6),
                    ((" ", 28, " "), 9, 5),  # disordered
                    ((" ", 29, " "), 7),
                    ((" ", 30, " "), 12),
                    ((" ", 31, " "), 6),
                    ((" ", 32, " "), 4),
                    ((" ", 33, " "), 11),
                    ((" ", 34, " "), 7),
                    ((" ", 35, " "), 6),
                    ((" ", 36, " "), 9),
                    ((" ", 37, " "), 8),
                    ((" ", 38, " "), 9),
                    ((" ", 39, " "), 6),
                    ((" ", 40, " "), 14),
                    ((" ", 41, " "), 6),
                    ((" ", 42, " "), 4),
                    ((" ", 43, " "), 9),
                    ((" ", 44, " "), 11),
                    ((" ", 45, " "), 6, 1),  # disordered
                    ((" ", 46, " "), 8),
                    ((" ", 47, " "), 10),
                    ((" ", 48, " "), 11),
                    ((" ", 49, " "), 6),
                    ((" ", 50, " "), 4),
                    ((" ", 51, " "), 5),
                    ((" ", 52, " "), 5),
                    ((" ", 53, " "), 7),
                    ((" ", 54, " "), 4),
                    ((" ", 55, " "), 8),
                    ((" ", 56, " "), 7),
                    ((" ", 57, " "), 7),
                    ((" ", 58, " "), 6),
                    ((" ", 59, " "), 4),
                    ((" ", 60, " "), 9),
                    ((" ", 61, " "), 8),
                    ((" ", 62, " "), 11),
                    ((" ", 63, " "), 6),
                    ((" ", 64, " "), 6),
                    ((" ", 65, " "), 6),
                    ((" ", 66, " "), 7),
                    ((" ", 67, " "), 10),
                    ((" ", 68, " "), 4),
                    ((" ", 69, " "), 14),
                    ((" ", 70, " "), 6),
                    ((" ", 71, " "), 4),
                    ((" ", 72, " "), 4),
                    ((" ", 73, " "), 4),
                    ((" ", 74, " "), 8, 3),  # disordered
                    ((" ", 75, " "), 8),
                    ((" ", 76, " "), 12),
                    ((" ", 77, " "), 6),
                    ((" ", 78, " "), 6),
                    ((" ", 79, " "), 4, 4),  # disordered
                    ((" ", 80, " "), 4, ["GLY", "SER"]),  # point mut
                    ((" ", 81, " "), 8, ["ASN", "LYS"]),  # point mut
                    ((" ", 82, " "), 6),
                    ((" ", 83, " "), 9),
                    ((" ", 84, " "), 12),
                    ((" ", 85, " "), 11),
                    ((" ", 86, " "), 6),
                ],
            ),
            (
                "B",
                11,
                [
                    ((" ", 44, " "), 11),
                    (("H_SEP", 45, " "), 10),  # Phosphoserine
                    ((" ", 46, " "), 8),
                    ((" ", 47, " "), 10),
                    ((" ", 48, " "), 11),
                    ((" ", 49, " "), 6),
                    ((" ", 50, " "), 4),
                    ((" ", 51, " "), 5),
                    ((" ", 51, "A"), 5),
                    ((" ", 52, " "), 7),
                    (("W", 0, " "), 1),
                ],
            ),
            (
                "C",
                5,
                [
                    (("W", 0, " "), 1),
                    (("H_NAG", 1, " "), 14),
                    (("H_NAG", 2, " "), 14),
                    (("H_NAG", 4, " "), 14),
                    (("H_NAG", 3, " "), 14),
                ],
            ),
            (
                " ",
                76,
                [
                    (("W", 1, " "), 1),
                    (("W", 2, " "), 1),
                    (("W", 3, " "), 1),
                    (("W", 4, " "), 1),
                    (("W", 5, " "), 1),
                    (("W", 6, " "), 1),
                    (("W", 7, " "), 1),
                    (("W", 8, " "), 1),
                    (("W", 9, " "), 1),
                    (("W", 10, " "), 1),
                    (("W", 11, " "), 1),
                    (("W", 12, " "), 1),
                    (("W", 13, " "), 1),
                    (("W", 14, " "), 1),
                    (("W", 15, " "), 1),
                    (("W", 16, " "), 1),
                    (("W", 17, " "), 1),
                    (("W", 18, " "), 1),
                    (("W", 19, " "), 1),
                    (("W", 20, " "), 1),
                    (("W", 21, " "), 1),
                    (("W", 22, " "), 1),
                    (("W", 23, " "), 1),
                    (("W", 24, " "), 1),
                    (("W", 25, " "), 1),
                    (("W", 26, " "), 1),
                    (("W", 27, " "), 1),
                    (("W", 28, " "), 1),
                    (("W", 29, " "), 1),
                    (("W", 30, " "), 1),
                    (("W", 31, " "), 1),
                    (("W", 32, " "), 1),
                    (("W", 33, " "), 1),
                    (("W", 34, " "), 1),
                    (("W", 35, " "), 1),
                    (("W", 36, " "), 1),
                    (("W", 37, " "), 1),
                    (("W", 38, " "), 1),
                    (("W", 39, " "), 1),
                    (("W", 40, " "), 1),
                    (("W", 41, " "), 1),
                    (("W", 42, " "), 1),
                    (("W", 43, " "), 1),
                    (("W", 44, " "), 1),
                    (("W", 45, " "), 1),
                    (("W", 46, " "), 1),
                    (("W", 47, " "), 1),
                    (("W", 48, " "), 1),
                    (("W", 49, " "), 1),
                    (("W", 50, " "), 1),
                    (("W", 51, " "), 1),
                    (("W", 52, " "), 1),
                    (("W", 53, " "), 1),
                    (("W", 54, " "), 1),
                    (("W", 55, " "), 1),
                    (("W", 56, " "), 1),
                    (("W", 57, " "), 1),
                    (("W", 58, " "), 1),
                    (("W", 59, " "), 1),
                    (("W", 60, " "), 1),
                    (("W", 61, " "), 1),
                    (("W", 62, " "), 1),
                    (("W", 63, " "), 1),
                    (("W", 64, " "), 1),
                    (("W", 65, " "), 1),
                    (("W", 66, " "), 1),
                    (("W", 67, " "), 1),
                    (("W", 68, " "), 1),
                    (("W", 69, " "), 1),
                    (("W", 70, " "), 1),
                    (("W", 71, " "), 1),
                    (("W", 72, " "), 1),
                    (("W", 73, " "), 1),
                    (("W", 74, " "), 1),
                    (("W", 75, " "), 1),
                    (("W", 77, " "), 1),
                ],
            ),
        ]

        for c_idx, chn in enumerate(chain_data):
            # Check chain ID and length
            chain = m1.get_list()[c_idx]
            self.assertEqual(chain.get_id(), chn[0])
            self.assertEqual(len(chain), chn[1])
            for r_idx, res in enumerate(chn[2]):
                residue = chain.get_list()[r_idx]
                # Check residue ID and atom count
                self.assertEqual(residue.get_id(), res[0])
                self.assertEqual(len(residue), res[1])
                disorder_lvl = residue.is_disordered()
                if disorder_lvl == 1:
                    # Check the number of disordered atoms
                    disordered_count = sum(
                        1 for atom in residue if atom.is_disordered()
                    )
                    if disordered_count:
                        self.assertEqual(disordered_count, res[2])
                elif disorder_lvl == 2:
                    # Point mutation -- check residue names
                    self.assertEqual(residue.disordered_get_id_list(), res[2])

    def test_structure_details(self):
        """Verify details of the parsed example PDB file."""
        structure = self.structure
        self.assertEqual(len(structure), 2)

        # First model
        model = structure[0]
        self.assertEqual(model.id, 0)
        self.assertEqual(model.level, "M")
        self.assertEqual(len(model), 1)
        chain = model["A"]
        self.assertEqual(chain.id, "A")
        self.assertEqual(chain.level, "C")
        self.assertEqual(len(chain), 1)
        self.assertEqual(" ".join(residue.resname for residue in chain), "PCA")
        self.assertEqual(
            " ".join(atom.name for atom in chain.get_atoms()),
            "N CA CB CG DA OE C O CA  ",
        )
        self.assertEqual(
            " ".join(atom.element for atom in chain.get_atoms()), "N C C C D O C O CA"
        )
        # Second model
        model = structure[1]
        self.assertEqual(model.id, 1)
        self.assertEqual(model.level, "M")
        self.assertEqual(len(model), 4)
        chain = model["A"]
        self.assertEqual(chain.id, "A")
        self.assertEqual(chain.level, "C")
        self.assertEqual(len(chain), 86)
        self.assertEqual(
            " ".join(residue.resname for residue in chain),
            "CYS ARG CYS GLY SER GLN GLY GLY GLY SER THR CYS "
            "PRO GLY LEU ARG CYS CYS SER ILE TRP GLY TRP CYS "
            "GLY ASP SER GLU PRO TYR CYS GLY ARG THR CYS GLU "
            "ASN LYS CYS TRP SER GLY GLU ARG SER ASP HIS ARG "
            "CYS GLY ALA ALA VAL GLY ASN PRO PRO CYS GLY GLN "
            "ASP ARG CYS CYS SER VAL HIS GLY TRP CYS GLY GLY "
            "GLY ASN ASP TYR CYS SER GLY GLY ASN CYS GLN TYR "
            "ARG CYS",
        )

        self.assertEqual(
            " ".join(atom.name for atom in chain.get_atoms()),
            "C N CA C O CB CG CD NE CZ NH1 NH2 N CA C O CB SG N "
            "CA C O N CA C O CB OG N CA C O CB CG CD OE1 NE2 N CA "
            "C O N CA C O N CA C O N CA C O CB OG N CA C O CB OG1 "
            "CG2 N CA C O CB SG N CA C O CB CG CD N CA C O N CA C "
            "O CB CG CD1 CD2 N CA C O CB CG CD NE CZ NH1 NH2 N CA "
            "C O CB SG N CA C O CB SG N CA C O CB OG N CA C O CB "
            "CG1 CG2 CD1 N CA C O CB CG CD1 CD2 NE1 CE2 CE3 CZ2 "
            "CZ3 CH2 N CA C O N CA C O CB CG CD1 CD2 NE1 CE2 CE3 "
            "CZ2 CZ3 CH2 N CA C O CB SG N CA C O N CA C O CB CG "
            "OD1 OD2 N CA C O CB OG N CA C O CB CG CD OE1 OE2 N "
            "CA C O CB CG CD N CA C O CB CG CD1 CD2 CE1 CE2 CZ OH "
            "N CA C O CB SG N CA C O N CA C O CB CG CD NE CZ NH1 "
            "NH2 N CA C O CB OG1 CG2 N CA C O CB SG N CA C O CB "
            "CG CD OE1 OE2 N CA C O CB CG OD1 ND2 N CA C O CB CG "
            "CD CE NZ N CA C O CB SG N CA C O CB CG CD1 CD2 NE1 "
            "CE2 CE3 CZ2 CZ3 CH2 N CA C O CB OG N CA C O N CA C "
            "O CB CG CD OE1 OE2 N CA C O CB CG CD NE CZ NH1 NH2 "
            "N CA C O CB OG N CA C O CB CG OD1 OD2 N CA C O CB "
            "CG ND1 CD2 CE1 NE2 N CA C O CB CG CD NE CZ NH1 NH2 "
            "N CA C O CB SG N CA C O N CA C O CB N CA C O CB N "
            "CA C O CB CG1 CG2 N CA C O N CA C O CB CG OD1 ND2 "
            "N CA C O CB CG CD N CA C O CB CG CD N CA C O CB SG "
            "N CA C O N CA C O CB CG CD OE1 NE2 N CA C O CB CG "
            "OD1 OD2 N CA C O CB CG CD NE CZ NH1 NH2 N CA C O CB "
            "SG N CA C O CB SG N CA C O CB OG N CA C O CB CG1 CG2 "
            "N CA C O CB CG ND1 CD2 CE1 NE2 N CA C O N CA C O CB "
            "CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2 N CA C O CB SG N "
            "CA C O N CA C O N CA C O CA N C O CB CG OD1 ND2 N CA "
            "C O CB CG OD1 OD2 N CA C O CB CG CD1 CD2 CE1 CE2 CZ "
            "OH N CA C O CB SG N CA C O CB OG N CA C O N CA C O N "
            "CA C O CB CG OD1 ND2 N CA C O CB SG N CA C O CB CG "
            "CD OE1 NE2 N CA C O CB CG CD1 CD2 CE1 CE2 CZ OH N CA "
            "C O CB CG CD NE CZ NH1 NH2 N CA C O CB SG",
        )

        self.assertEqual(
            " ".join(atom.element for atom in chain.get_atoms()),
            "C N C C O C C C N C N N N C C O C S N C C O N C C O "
            "C O N C C O C C C O N N C C O N C C O N C C O N C C "
            "O C O N C C O C O C N C C O C S N C C O C C C N C C "
            "O N C C O C C C C N C C O C C C N C N N N C C O C S "
            "N C C O C S N C C O C O N C C O C C C C N C C O C C "
            "C C N C C C C C N C C O N C C O C C C C N C C C C C "
            "N C C O C S N C C O N C C O C C O O N C C O C O N C "
            "C O C C C O O N C C O C C C N C C O C C C C C C C O "
            "N C C O C S N C C O N C C O C C C N C N N N C C O C "
            "O C N C C O C S N C C O C C C O O N C C O C C O N N "
            "C C O C C C C N N C C O C S N C C O C C C C N C C C "
            "C C N C C O C O N C C O N C C O C C C O O N C C O C "
            "C C N C N N N C C O C O N C C O C C O O N C C O C C "
            "N C C N N C C O C C C N C N N N C C O C S N C C O N "
            "C C O C N C C O C N C C O C C C N C C O N C C O C C "
            "O N N C C O C C C N C C O C C C N C C O C S N C C O "
            "N C C O C C C O N N C C O C C O O N C C O C C C N C "
            "N N N C C O C S N C C O C S N C C O C O N C C O C C "
            "C N C C O C C N C C N N C C O N C C O C C C C N C C "
            "C C C N C C O C S N C C O N C C O N C C O C N C O C "
            "C O N N C C O C C O O N C C O C C C C C C C O N C C "
            "O C S N C C O C O N C C O N C C O N C C O C C O N N "
            "C C O C S N C C O C C C O N N C C O C C C C C C C O "
            "N C C O C C C N C N N N C C O C S",
        )


class ParseRealPDB_tests(unittest.TestCase):
    """Testing with real PDB files."""

    @classmethod
    def setUpClass(cls):
        cls.permissive = PDBParser()
        cls.strict = PDBParser(PERMISSIVE=False)

    def test_empty(self):
        """Parse an empty file."""
        handle = StringIO()
        with self.assertRaises(ValueError) as context_manager:
            _ = self.permissive.get_structure("MT", handle)
        self.assertEqual(str(context_manager.exception), "Empty file.")

    def test_SMCRA(self):
        """Walk down the structure hierarchy and test parser reliability."""
        s = self.permissive.get_structure("scr", "PDB/1A8O.pdb")
        for m in s:
            p = m.get_parent()
            self.assertEqual(s, p)
            for c in m:
                p = c.get_parent()
                self.assertEqual(m, p)
                for r in c:
                    p = r.get_parent()
                    self.assertEqual(c, p)
                    for a in r:
                        p = a.get_parent()
                        self.assertEqual(r.get_resname(), p.get_resname())

    def test_1A8O_strict(self):
        """Parse 1A8O.pdb file in strict mode."""
        structure = self.strict.get_structure("example", "PDB/1A8O.pdb")
        self.assertEqual(len(structure), 1)
        model = structure[0]
        self.assertEqual(model.id, 0)
        self.assertEqual(model.level, "M")
        self.assertEqual(len(model), 1)
        chain = model["A"]
        self.assertEqual(chain.id, "A")
        self.assertEqual(chain.level, "C")
        self.assertEqual(len(chain), 158)
        self.assertEqual(
            " ".join(residue.resname for residue in chain),
            "MSE ASP ILE ARG GLN GLY PRO LYS GLU PRO PHE ARG "
            "ASP TYR VAL ASP ARG PHE TYR LYS THR LEU ARG ALA "
            "GLU GLN ALA SER GLN GLU VAL LYS ASN TRP MSE THR "
            "GLU THR LEU LEU VAL GLN ASN ALA ASN PRO ASP CYS "
            "LYS THR ILE LEU LYS ALA LEU GLY PRO GLY ALA THR "
            "LEU GLU GLU MSE MSE THR ALA CYS GLN GLY HOH HOH "
            "HOH HOH HOH HOH HOH HOH HOH HOH HOH HOH HOH HOH "
            "HOH HOH HOH HOH HOH HOH HOH HOH HOH HOH HOH HOH "
            "HOH HOH HOH HOH HOH HOH HOH HOH HOH HOH HOH HOH "
            "HOH HOH HOH HOH HOH HOH HOH HOH HOH HOH HOH HOH "
            "HOH HOH HOH HOH HOH HOH HOH HOH HOH HOH HOH HOH "
            "HOH HOH HOH HOH HOH HOH HOH HOH HOH HOH HOH HOH "
            "HOH HOH HOH HOH HOH HOH HOH HOH HOH HOH HOH HOH "
            "HOH HOH",
        )
        self.assertEqual(
            " ".join(atom.name for atom in chain.get_atoms()),
            "N CA C O CB CG SE CE N CA C O CB CG OD1 OD2 N CA "
            "C O CB CG1 CG2 CD1 N CA C O CB CG CD NE CZ NH1 "
            "NH2 N CA C O CB CG CD OE1 NE2 N CA C O N CA C O "
            "CB CG CD N CA C O CB CG CD CE NZ N CA C O CB CG "
            "CD OE1 OE2 N CA C O CB CG CD N CA C O CB CG CD1 "
            "CD2 CE1 CE2 CZ N CA C O CB CG CD NE CZ NH1 NH2 N "
            "CA C O CB CG OD1 OD2 N CA C O CB CG CD1 CD2 CE1 "
            "CE2 CZ OH N CA C O CB CG1 CG2 N CA C O CB CG OD1 "
            "OD2 N CA C O CB CG CD NE CZ NH1 NH2 N CA C O CB "
            "CG CD1 CD2 CE1 CE2 CZ N CA C O CB CG CD1 CD2 CE1 "
            "CE2 CZ OH N CA C O CB CG CD CE NZ N CA C O CB "
            "OG1 CG2 N CA C O CB CG CD1 CD2 N CA C O CB CG CD "
            "NE CZ NH1 NH2 N CA C O CB N CA C O CB CG CD OE1 "
            "OE2 N CA C O CB CG CD OE1 NE2 N CA C O CB N CA C "
            "O CB OG N CA C O CB CG CD OE1 NE2 N CA C O CB CG "
            "CD OE1 OE2 N CA C O CB CG1 CG2 N CA C O CB CG CD "
            "CE NZ N CA C O CB CG OD1 ND2 N CA C O CB CG CD1 "
            "CD2 NE1 CE2 CE3 CZ2 CZ3 CH2 N CA C O CB CG SE CE "
            "N CA C O CB OG1 CG2 N CA C O CB CG CD OE1 OE2 N "
            "CA C O CB OG1 CG2 N CA C O CB CG CD1 CD2 N CA C "
            "O CB CG CD1 CD2 N CA C O CB CG1 CG2 N CA C O CB "
            "CG CD OE1 NE2 N CA C O CB CG OD1 ND2 N CA C O CB "
            "N CA C O CB CG OD1 ND2 N CA C O CB CG CD N CA C "
            "O CB CG OD1 OD2 N CA C O CB SG N CA C O CB CG CD "
            "CE NZ N CA C O CB OG1 CG2 N CA C O CB CG1 CG2 "
            "CD1 N CA C O CB CG CD1 CD2 N CA C O CB CG CD CE "
            "NZ N CA C O CB N CA C O CB CG CD1 CD2 N CA C O N "
            "CA C O CB CG CD N CA C O N CA C O CB N CA C O CB "
            "OG1 CG2 N CA C O CB CG CD1 CD2 N CA C O CB CG CD "
            "OE1 OE2 N CA C O CB CG CD OE1 OE2 N CA C O CB CG "
            "SE CE N CA C O CB CG SE CE N CA C O CB OG1 CG2 N "
            "CA C O CB N CA C O CB SG N CA C O CB CG CD OE1 "
            "NE2 N CA C O OXT O O O O O O O O O O O O O O O O "
            "O O O O O O O O O O O O O O O O O O O O O O O O "
            "O O O O O O O O O O O O O O O O O O O O O O O O "
            "O O O O O O O O O O O O O O O O O O O O O O O O",
        )
        self.assertEqual(
            " ".join(atom.element for atom in chain.get_atoms()),
            "N C C O C C SE C N C C O C C O O N C C O C C C C "
            "N C C O C C C N C N N N C C O C C C O N N C C O "
            "N C C O C C C N C C O C C C C N N C C O C C C O "
            "O N C C O C C C N C C O C C C C C C C N C C O C "
            "C C N C N N N C C O C C O O N C C O C C C C C C "
            "C O N C C O C C C N C C O C C O O N C C O C C C "
            "N C N N N C C O C C C C C C C N C C O C C C C C "
            "C C O N C C O C C C C N N C C O C O C N C C O C "
            "C C C N C C O C C C N C N N N C C O C N C C O C "
            "C C O O N C C O C C C O N N C C O C N C C O C O "
            "N C C O C C C O N N C C O C C C O O N C C O C C "
            "C N C C O C C C C N N C C O C C O N N C C O C C "
            "C C N C C C C C N C C O C C SE C N C C O C O C N "
            "C C O C C C O O N C C O C O C N C C O C C C C N "
            "C C O C C C C N C C O C C C N C C O C C C O N N "
            "C C O C C O N N C C O C N C C O C C O N N C C O "
            "C C C N C C O C C O O N C C O C S N C C O C C C "
            "C N N C C O C O C N C C O C C C C N C C O C C C "
            "C N C C O C C C C N N C C O C N C C O C C C C N "
            "C C O N C C O C C C N C C O N C C O C N C C O C "
            "O C N C C O C C C C N C C O C C C O O N C C O C "
            "C C O O N C C O C C SE C N C C O C C SE C N C C "
            "O C O C N C C O C N C C O C S N C C O C C C O N "
            "N C C O O O O O O O O O O O O O O O O O O O O O "
            "O O O O O O O O O O O O O O O O O O O O O O O O "
            "O O O O O O O O O O O O O O O O O O O O O O O O "
            "O O O O O O O O O O O O O O O O O O O O O",
        )

    def test_duplicated_residue_permissive(self):
        """Catch exception on duplicated residue."""
        data = (
            "HETATM 6289  O   HOH     5      28.182  -5.239  31.370  1.00 22.99           O\n"
            "HETATM 6513  O   HOH     6      21.829   3.361  14.003  1.00 14.25           O\n"
            "HETATM 6607  O   HOH     5      33.861  40.044  18.022  1.00 18.73           O\n"
            "END   \n"
        )

        with warnings.catch_warnings(record=True) as w:
            s = self.permissive.get_structure("example", StringIO(data))
            self.assertEqual(len(w), 1)

        reslist = list(s.get_residues())
        n_res = len(reslist)
        resids = [r.id[1] for r in reslist]
        self.assertEqual(n_res, 2)
        self.assertEqual(resids, [5, 6])

    def test_duplicated_residue_strict(self):
        """Throw exception on duplicated residue."""
        data = (
            "HETATM 6289  O   HOH     5      28.182  -5.239  31.370  1.00 22.99           O\n"
            "HETATM 6513  O   HOH     6      21.829   3.361  14.003  1.00 14.25           O\n"
            "HETATM 6607  O   HOH     5      33.861  40.044  18.022  1.00 18.73           O\n"
            "END   \n"
        )

        with self.assertRaises(PDBConstructionException):
            _ = self.strict.get_structure("example", StringIO(data))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
