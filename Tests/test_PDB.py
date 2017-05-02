# Copyright 2009-2011 by Eric Talevich.  All rights reserved.
# Revisions copyright 2009-2013 by Peter Cock.  All rights reserved.
# Revisions copyright 2013 Lenna X. Peterson. All rights reserved.
#
# Converted by Eric Talevich from an older unit test copyright 2002
# by Thomas Hamelryck.
#
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for the Bio.PDB module."""
from __future__ import print_function

from copy import deepcopy
import os
import sys
import tempfile
import unittest
import warnings
from Bio._py3k import StringIO

try:
    import numpy
    from numpy import dot  # Missing on old PyPy's micronumpy
    del dot
    from numpy.linalg import svd, det  # Missing in PyPy 2.0 numpypy
    from numpy.random import random
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.PDB.")

from Bio import BiopythonWarning
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein
from Bio.PDB import PDBParser, PPBuilder, CaPPBuilder, PDBIO, Select
from Bio.PDB import HSExposureCA, HSExposureCB, ExposureCN
from Bio.PDB.PDBExceptions import PDBConstructionException, PDBConstructionWarning
from Bio.PDB import rotmat, Vector, refmat, calc_angle, calc_dihedral, rotaxis, m2rotaxis
from Bio.PDB import Residue, Atom, StructureAlignment, Superimposer, Selection
from Bio.PDB import make_dssp_dict
from Bio.PDB import DSSP
from Bio.PDB.NACCESS import process_asa_data, process_rsa_data
from Bio.PDB.ResidueDepth import _get_atom_radius
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


# NB: the 'A_' prefix ensures this test case is run first
class A_ExceptionTest(unittest.TestCase):
    """Errors and warnings while parsing of flawed PDB files.

    These tests must be executed because of the way Python's warnings module
    works -- a warning is only logged the first time it is encountered.
    """
    def test_1_warnings(self):
        """Check warnings: Parse a flawed PDB file in permissive mode."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always', PDBConstructionWarning)

            # Trigger warnings
            p = PDBParser(PERMISSIVE=True)
            p.get_structure("example", "PDB/a_structure.pdb")

            self.assertEqual(len(w), 14)
            for wrn, msg in zip(w, [
                    # Expected warning messages:
                    "Used element 'N' for Atom (name=N) with given element ''",
                    "Used element 'C' for Atom (name=CA) with given element ''",
                    "Atom names ' CA ' and 'CA  ' differ only in spaces at line 17.",
                    "Used element 'CA' for Atom (name=CA  ) with given element ''",
                    'Atom N defined twice in residue <Residue ARG het=  resseq=2 icode= > at line 21.',
                    'disordered atom found with blank altloc before line 33.',
                    "Residue (' ', 4, ' ') redefined at line 43.",
                    "Blank altlocs in duplicate residue SER (' ', 4, ' ') at line 43.",
                    "Residue (' ', 10, ' ') redefined at line 75.",
                    "Residue (' ', 14, ' ') redefined at line 106.",
                    "Residue (' ', 16, ' ') redefined at line 135.",
                    "Residue (' ', 80, ' ') redefined at line 633.",
                    "Residue (' ', 81, ' ') redefined at line 646.",
                    'Atom O defined twice in residue <Residue HOH het=W resseq=67 icode= > at line 822.'
                    ]):
                self.assertIn(msg, str(wrn))

    def test_2_strict(self):
        """Check error: Parse a flawed PDB file in strict mode."""
        parser = PDBParser(PERMISSIVE=False)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", PDBConstructionWarning)
            self.assertRaises(PDBConstructionException,
                              parser.get_structure, "example", "PDB/a_structure.pdb")
            self.assertEqual(len(w), 4, w)

    def test_3_bad_xyz(self):
        """Check error: Parse an entry with bad x,y,z value."""
        data = "ATOM      9  N   ASP A 152      21.554  34.953  27.691  1.00 19.26           N\n"
        parser = PDBParser(PERMISSIVE=False)
        s = parser.get_structure("example", StringIO(data))
        data = "ATOM      9  N   ASP A 152      21.ish  34.953  27.691  1.00 19.26           N\n"
        self.assertRaises(PDBConstructionException,
                          parser.get_structure, "example", StringIO(data))

    def test_4_occupancy(self):
        """Parse file with missing occupancy"""
        permissive = PDBParser(PERMISSIVE=True)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", PDBConstructionWarning)
            structure = permissive.get_structure("test", "PDB/occupancy.pdb")
            self.assertEqual(len(w), 3, w)
        atoms = structure[0]['A'][(' ', 152, ' ')]
        # Blank occupancy behavior set in Bio/PDB/PDBParser
        self.assertEqual(atoms['N'].get_occupancy(), None)
        self.assertEqual(atoms['CA'].get_occupancy(), 1.0)
        self.assertEqual(atoms['C'].get_occupancy(), 0.0)

        strict = PDBParser(PERMISSIVE=False)
        self.assertRaises(PDBConstructionException,
                          strict.get_structure, "test", "PDB/occupancy.pdb")


class HeaderTests(unittest.TestCase):
    """Tests for parse_pdb_header."""

    def test_capsid(self):
        """Parse the header of a known PDB file (1A8O)."""
        parser = PDBParser()
        struct = parser.get_structure('1A8O', 'PDB/1A8O.pdb')
        self.assertAlmostEqual(struct.header['resolution'], 1.7)
        # Case-insensitive string comparisons
        known_strings = {
            'author': 'T.R.Gamble,S.Yoo,F.F.Vajdos,U.K.Von Schwedler,D.K.Worthylake,H.Wang,J.P.Mccutcheon,W.I.Sundquist,C.P.Hill',
            'deposition_date': '1998-03-27',
            'head': 'viral protein',
            'journal': 'AUTH   T.R.GAMBLE,S.YOO,F.F.VAJDOS,U.K.VON SCHWEDLER,AUTH 2 D.K.WORTHYLAKE,H.WANG,J.P.MCCUTCHEON,W.I.SUNDQUIST,AUTH 3 C.P.HILLTITL   STRUCTURE OF THE CARBOXYL-TERMINAL DIMERIZATIONTITL 2 DOMAIN OF THE HIV-1 CAPSID PROTEIN.REF    SCIENCE                       V. 278   849 1997REFN                   ISSN 0036-8075PMID   9346481DOI    10.1126/SCIENCE.278.5339.849',
            'journal_reference': 't.r.gamble,s.yoo,f.f.vajdos,u.k.von schwedler, d.k.worthylake,h.wang,j.p.mccutcheon,w.i.sundquist, c.p.hill structure of the carboxyl-terminal dimerization domain of the hiv-1 capsid protein. science v. 278 849 1997 issn 0036-8075 9346481 10.1126/science.278.5339.849 ',
            'keywords': 'capsid, core protein, hiv, c-terminal domain, viral protein',
            'name': ' hiv capsid c-terminal domain',
            'release_date': '1998-10-14',
            'structure_method': 'x-ray diffraction',
        }
        for key, expect in known_strings.items():
            self.assertEqual(struct.header[key].lower(), expect.lower())

    def test_fibril(self):
        """Parse the header of another PDB file (2BEG)."""
        parser = PDBParser()
        struct = parser.get_structure('2BEG', 'PDB/2BEG.pdb')
        known_strings = {
            'author': 'T.Luhrs,C.Ritter,M.Adrian,D.Riek-Loher,B.Bohrmann,H.Dobeli,D.Schubert,R.Riek',
            'deposition_date': '2005-10-24',
            'head': 'protein fibril',
            'journal': "AUTH   T.LUHRS,C.RITTER,M.ADRIAN,D.RIEK-LOHER,B.BOHRMANN,AUTH 2 H.DOBELI,D.SCHUBERT,R.RIEKTITL   3D STRUCTURE OF ALZHEIMER'S AMYLOID-{BETA}(1-42)TITL 2 FIBRILS.REF    PROC.NATL.ACAD.SCI.USA        V. 102 17342 2005REFN                   ISSN 0027-8424PMID   16293696DOI    10.1073/PNAS.0506723102",
            'journal_reference': "t.luhrs,c.ritter,m.adrian,d.riek-loher,b.bohrmann, h.dobeli,d.schubert,r.riek 3d structure of alzheimer's amyloid-{beta}(1-42) fibrils. proc.natl.acad.sci.usa v. 102 17342 2005 issn 0027-8424 16293696 10.1073/pnas.0506723102 ",
            'keywords': "alzheimer's, fibril, protofilament, beta-sandwich, quenched hydrogen/deuterium exchange, pairwise mutagenesis, protein fibril",
            'name': " 3d structure of alzheimer's abeta(1-42) fibrils",
            'release_date': '2005-11-22',
            'structure_method': 'solution nmr',
        }
        for key, expect in known_strings.items():
            self.assertEqual(struct.header[key].lower(), expect.lower())


class ParseTest(unittest.TestCase):
    def setUp(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            p = PDBParser(PERMISSIVE=1)
            self.structure = p.get_structure("example", "PDB/a_structure.pdb")

    def test_c_n(self):
        """Extract polypeptides using C-N."""
        ppbuild = PPBuilder()
        polypeptides = ppbuild.build_peptides(self.structure[1])
        self.assertEqual(len(polypeptides), 1)
        pp = polypeptides[0]
        # Check the start and end positions
        self.assertEqual(pp[0].get_id()[1], 2)
        self.assertEqual(pp[-1].get_id()[1], 86)
        # Check the sequence
        s = pp.get_sequence()
        self.assertTrue(isinstance(s, Seq))
        self.assertEqual(s.alphabet, generic_protein)
        self.assertEqual("RCGSQGGGSTCPGLRCCSIWGWCGDSEPYCGRTCENKCWSGER"
                         "SDHRCGAAVGNPPCGQDRCCSVHGWCGGGNDYCSGGNCQYRC",
                         str(s))

    def test_ca_ca(self):
        """Extract polypeptides using CA-CA."""
        ppbuild = CaPPBuilder()
        polypeptides = ppbuild.build_peptides(self.structure[1])
        self.assertEqual(len(polypeptides), 1)
        pp = polypeptides[0]
        # Check the start and end positions
        self.assertEqual(pp[0].get_id()[1], 2)
        self.assertEqual(pp[-1].get_id()[1], 86)
        # Check the sequence
        s = pp.get_sequence()
        self.assertTrue(isinstance(s, Seq))
        self.assertEqual(s.alphabet, generic_protein)
        self.assertEqual("RCGSQGGGSTCPGLRCCSIWGWCGDSEPYCGRTCENKCWSGER"
                         "SDHRCGAAVGNPPCGQDRCCSVHGWCGGGNDYCSGGNCQYRC",
                         str(s))

    def test_structure(self):
        """Verify the structure of the parsed example PDB file."""
        # Structure contains 2 models
        self.assertEqual(len(self.structure), 2)
        # --- Checking model 0 ---
        m0 = self.structure[0]
        # Model 0 contains 1 chain
        self.assertEqual(len(m0), 1)
        # Chain 'A' contains 1 residue
        self.assertEqual(len(m0['A']), 1)
        # Residue ('H_PCA', 1, ' ') contains 8 atoms.
        residue = m0['A'].get_list()[0]
        self.assertEqual(residue.get_id(), ('H_PCA', 1, ' '))
        self.assertEqual(len(residue), 9)
        # --- Checking model 1 ---
        m1 = self.structure[1]
        # Model 1 contains 3 chains
        self.assertEqual(len(m1), 3)
        # Deconstruct this data structure to check each chain
        chain_data = [  # chain_id, chain_len, [(residue_id, residue_len), ...]
            ('A', 86, [((' ', 0, ' '), 1),
                       ((' ', 2, ' '), 11),
                       ((' ', 3, ' '), 6, 1),  # disordered
                       ((' ', 4, ' '), 4),
                       ((' ', 5, ' '), 6),
                       ((' ', 6, ' '), 9),
                       ((' ', 7, ' '), 4),
                       ((' ', 8, ' '), 4),
                       ((' ', 9, ' '), 4),
                       ((' ', 10, ' '), 6, ['GLY', 'SER']),  # point mut
                       ((' ', 11, ' '), 7),
                       ((' ', 12, ' '), 6),
                       ((' ', 13, ' '), 7),
                       ((' ', 14, ' '), 4, ['ALA', 'GLY']),  # point mut
                       ((' ', 15, ' '), 8, 3),  # disordered
                       ((' ', 16, ' '), 11, ['ARG', 'TRP']),  # point mut
                       ((' ', 17, ' '), 6),
                       ((' ', 18, ' '), 6),
                       ((' ', 19, ' '), 6),
                       ((' ', 20, ' '), 8),
                       ((' ', 21, ' '), 14),
                       ((' ', 22, ' '), 4),
                       ((' ', 23, ' '), 14),
                       ((' ', 24, ' '), 6),
                       ((' ', 25, ' '), 4),
                       ((' ', 26, ' '), 8),
                       ((' ', 27, ' '), 6),
                       ((' ', 28, ' '), 9, 5),  # disordered
                       ((' ', 29, ' '), 7),
                       ((' ', 30, ' '), 12),
                       ((' ', 31, ' '), 6),
                       ((' ', 32, ' '), 4),
                       ((' ', 33, ' '), 11),
                       ((' ', 34, ' '), 7),
                       ((' ', 35, ' '), 6),
                       ((' ', 36, ' '), 9),
                       ((' ', 37, ' '), 8),
                       ((' ', 38, ' '), 9),
                       ((' ', 39, ' '), 6),
                       ((' ', 40, ' '), 14),
                       ((' ', 41, ' '), 6),
                       ((' ', 42, ' '), 4),
                       ((' ', 43, ' '), 9),
                       ((' ', 44, ' '), 11),
                       ((' ', 45, ' '), 6, 1),  # disordered
                       ((' ', 46, ' '), 8),
                       ((' ', 47, ' '), 10),
                       ((' ', 48, ' '), 11),
                       ((' ', 49, ' '), 6),
                       ((' ', 50, ' '), 4),
                       ((' ', 51, ' '), 5),
                       ((' ', 52, ' '), 5),
                       ((' ', 53, ' '), 7),
                       ((' ', 54, ' '), 4),
                       ((' ', 55, ' '), 8),
                       ((' ', 56, ' '), 7),
                       ((' ', 57, ' '), 7),
                       ((' ', 58, ' '), 6),
                       ((' ', 59, ' '), 4),
                       ((' ', 60, ' '), 9),
                       ((' ', 61, ' '), 8),
                       ((' ', 62, ' '), 11),
                       ((' ', 63, ' '), 6),
                       ((' ', 64, ' '), 6),
                       ((' ', 65, ' '), 6),
                       ((' ', 66, ' '), 7),
                       ((' ', 67, ' '), 10),
                       ((' ', 68, ' '), 4),
                       ((' ', 69, ' '), 14),
                       ((' ', 70, ' '), 6),
                       ((' ', 71, ' '), 4),
                       ((' ', 72, ' '), 4),
                       ((' ', 73, ' '), 4),
                       ((' ', 74, ' '), 8, 3),  # disordered
                       ((' ', 75, ' '), 8),
                       ((' ', 76, ' '), 12),
                       ((' ', 77, ' '), 6),
                       ((' ', 78, ' '), 6),
                       ((' ', 79, ' '), 4, 4),  # disordered
                       ((' ', 80, ' '), 4, ['GLY', 'SER']),  # point mut
                       ((' ', 81, ' '), 8, ['ASN', 'LYS']),  # point mut
                       ((' ', 82, ' '), 6),
                       ((' ', 83, ' '), 9),
                       ((' ', 84, ' '), 12),
                       ((' ', 85, ' '), 11),
                       ((' ', 86, ' '), 6),
                       ]),
            ('B', 4, [(('H_NAG', 1, ' '), 14),
                      (('H_NAG', 2, ' '), 14),
                      (('H_NAG', 3, ' '), 14),
                      (('H_NAG', 4, ' '), 14),
                      ]),
            (' ', 76, [(('W', 1, ' '), 1),
                       (('W', 2, ' '), 1),
                       (('W', 3, ' '), 1),
                       (('W', 4, ' '), 1),
                       (('W', 5, ' '), 1),
                       (('W', 6, ' '), 1),
                       (('W', 7, ' '), 1),
                       (('W', 8, ' '), 1),
                       (('W', 9, ' '), 1),
                       (('W', 10, ' '), 1),
                       (('W', 11, ' '), 1),
                       (('W', 12, ' '), 1),
                       (('W', 13, ' '), 1),
                       (('W', 14, ' '), 1),
                       (('W', 15, ' '), 1),
                       (('W', 16, ' '), 1),
                       (('W', 17, ' '), 1),
                       (('W', 18, ' '), 1),
                       (('W', 19, ' '), 1),
                       (('W', 20, ' '), 1),
                       (('W', 21, ' '), 1),
                       (('W', 22, ' '), 1),
                       (('W', 23, ' '), 1),
                       (('W', 24, ' '), 1),
                       (('W', 25, ' '), 1),
                       (('W', 26, ' '), 1),
                       (('W', 27, ' '), 1),
                       (('W', 28, ' '), 1),
                       (('W', 29, ' '), 1),
                       (('W', 30, ' '), 1),
                       (('W', 31, ' '), 1),
                       (('W', 32, ' '), 1),
                       (('W', 33, ' '), 1),
                       (('W', 34, ' '), 1),
                       (('W', 35, ' '), 1),
                       (('W', 36, ' '), 1),
                       (('W', 37, ' '), 1),
                       (('W', 38, ' '), 1),
                       (('W', 39, ' '), 1),
                       (('W', 40, ' '), 1),
                       (('W', 41, ' '), 1),
                       (('W', 42, ' '), 1),
                       (('W', 43, ' '), 1),
                       (('W', 44, ' '), 1),
                       (('W', 45, ' '), 1),
                       (('W', 46, ' '), 1),
                       (('W', 47, ' '), 1),
                       (('W', 48, ' '), 1),
                       (('W', 49, ' '), 1),
                       (('W', 50, ' '), 1),
                       (('W', 51, ' '), 1),
                       (('W', 52, ' '), 1),
                       (('W', 53, ' '), 1),
                       (('W', 54, ' '), 1),
                       (('W', 55, ' '), 1),
                       (('W', 56, ' '), 1),
                       (('W', 57, ' '), 1),
                       (('W', 58, ' '), 1),
                       (('W', 59, ' '), 1),
                       (('W', 60, ' '), 1),
                       (('W', 61, ' '), 1),
                       (('W', 62, ' '), 1),
                       (('W', 63, ' '), 1),
                       (('W', 64, ' '), 1),
                       (('W', 65, ' '), 1),
                       (('W', 66, ' '), 1),
                       (('W', 67, ' '), 1),
                       (('W', 68, ' '), 1),
                       (('W', 69, ' '), 1),
                       (('W', 70, ' '), 1),
                       (('W', 71, ' '), 1),
                       (('W', 72, ' '), 1),
                       (('W', 73, ' '), 1),
                       (('W', 74, ' '), 1),
                       (('W', 75, ' '), 1),
                       (('W', 77, ' '), 1),
                       ])
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
                    disordered_count = sum(1 for atom in residue
                                           if atom.is_disordered())
                    if disordered_count:
                        self.assertEqual(disordered_count, res[2])
                elif disorder_lvl == 2:
                    # Point mutation -- check residue names
                    self.assertEqual(residue.disordered_get_id_list(), res[2])

    def test_details(self):
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
        self.assertEqual(" ".join(atom.name for atom in chain.get_atoms()),
                         "N CA CB CG DA OE C O CA  ")
        self.assertEqual(" ".join(atom.element for atom in chain.get_atoms()),
                         "N C C C D O C O CA")
        # Second model
        model = structure[1]
        self.assertEqual(model.id, 1)
        self.assertEqual(model.level, "M")
        self.assertEqual(len(model), 3)
        chain = model["A"]
        self.assertEqual(chain.id, "A")
        self.assertEqual(chain.level, "C")
        self.assertEqual(len(chain), 86)
        self.assertEqual(" ".join(residue.resname for residue in chain),
                         "CYS ARG CYS GLY SER GLN GLY GLY GLY SER THR CYS "
                         "PRO GLY LEU ARG CYS CYS SER ILE TRP GLY TRP CYS "
                         "GLY ASP SER GLU PRO TYR CYS GLY ARG THR CYS GLU "
                         "ASN LYS CYS TRP SER GLY GLU ARG SER ASP HIS ARG "
                         "CYS GLY ALA ALA VAL GLY ASN PRO PRO CYS GLY GLN "
                         "ASP ARG CYS CYS SER VAL HIS GLY TRP CYS GLY GLY "
                         "GLY ASN ASP TYR CYS SER GLY GLY ASN CYS GLN TYR "
                         "ARG CYS")
        self.assertEqual(" ".join(atom.name for atom in chain.get_atoms()),
                         "C N CA C O CB CG CD NE CZ NH1 NH2 N CA C O CB SG "
                         "N CA C O N CA C O CB OG N CA C O CB CG CD OE1 NE2 "
                         "N CA C O N CA C O N CA C O N CA C O CB OG N CA C "
                         "O CB OG1 CG2 N CA C O CB SG N CA C O CB CG CD N "
                         "CA C O N CA C O CB CG CD1 CD2 N CA C O CB CG CD NE "
                         "CZ NH1 NH2 N CA C O CB SG N CA C O CB SG N CA C O "
                         "CB OG N CA C O CB CG1 CG2 CD1 N CA C O CB CG CD1 "
                         "CD2 NE1 CE2 CE3 CZ2 CZ3 CH2 N CA C O N CA C O CB "
                         "CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2 N CA C O CB SG "
                         "N CA C O N CA C O CB CG OD1 OD2 N CA C O CB OG N "
                         "CA C O CB CG CD OE1 OE2 N CA C O CB CG CD N CA C O "
                         "CB CG CD1 CD2 CE1 CE2 CZ OH N CA C O CB SG N CA C "
                         "O N CA C O CB CG CD NE CZ NH1 NH2 N CA C O CB OG1 "
                         "CG2 N CA C O CB SG N CA C O CB CG CD OE1 OE2 N CA "
                         "C O CB CG OD1 ND2 N CA C O CB CG CD CE NZ N CA C O "
                         "CB SG N CA C O CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 "
                         "CH2 N CA C O CB OG N CA C O N CA C O CB CG CD OE1 "
                         "OE2 N CA C O CB CG CD NE CZ NH1 NH2 N CA C O CB OG "
                         "N CA C O CB CG OD1 OD2 N CA C O CB CG ND1 CD2 CE1 "
                         "NE2 N CA C O CB CG CD NE CZ NH1 NH2 N CA C O CB SG "
                         "N CA C O N CA C O CB N CA C O CB N CA C O CB CG1 "
                         "CG2 N CA C O N CA C O CB CG OD1 ND2 N CA C O CB CG "
                         "CD N CA C O CB CG CD N CA C O CB SG N CA C O N CA "
                         "C O CB CG CD OE1 NE2 N CA C O CB CG OD1 OD2 N CA C "
                         "O CB CG CD NE CZ NH1 NH2 N CA C O CB SG N CA C O "
                         "CB SG N CA C O CB OG N CA C O CB CG1 CG2 N CA C O "
                         "CB CG ND1 CD2 CE1 NE2 N CA C O N CA C O CB CG CD1 "
                         "CD2 NE1 CE2 CE3 CZ2 CZ3 CH2 N CA C O CB SG N CA C "
                         "O N CA C O N CA C O N CA C O CB CG OD1 ND2 N CA C O "
                         "CB CG OD1 OD2 N CA C O CB CG CD1 CD2 CE1 CE2 CZ OH "
                         "N CA C O CB SG N CA C O CB OG N CA C O N CA C O N "
                         "CA C O CB CG OD1 ND2 N CA C O CB SG N CA C O CB CG "
                         "CD OE1 NE2 N CA C O CB CG CD1 CD2 CE1 CE2 CZ OH N "
                         "CA C O CB CG CD NE CZ NH1 NH2 N CA C O CB SG")
        self.assertEqual(" ".join(atom.element for atom in chain.get_atoms()),
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
                         "C C C N C C O C S N C C O N C C O N C C O N C C O C "
                         "C O N N C C O C C O O N C C O C C C C C C C O N C C "
                         "O C S N C C O C O N C C O N C C O N C C O C C O N N "
                         "C C O C S N C C O C C C O N N C C O C C C C C C C O "
                         "N C C O C C C N C N N N C C O C S")

    def test_pdbio_write_truncated(self):
        """Test parsing of truncated lines"""
        io = PDBIO()
        struct = self.structure
        # Write to temp file
        io.set_structure(struct)
        filenumber, filename = tempfile.mkstemp()
        os.close(filenumber)
        try:
            io.save(filename)
            # Check if there are lines besides 'ATOM', 'TER' and 'END'
            with open(filename, 'rU') as handle:
                record_set = set(l[0:6] for l in handle)
            record_set -= set(('ATOM  ', 'HETATM', 'MODEL ', 'ENDMDL', 'TER\n',  'TER   ', 'END\n', 'END   '))
            self.assertEqual(record_set, set())
        finally:
            os.remove(filename)

    def test_deepcopy_of_structure_with_disorder(self):
            """Test deepcopy of a structure with disordered atoms.
            Shouldn't cause recursion.
            """
            _ = deepcopy(self.structure)


class ParseReal(unittest.TestCase):
    """Testing with real PDB files."""

    def test_empty(self):
        """Parse an empty file."""
        parser = PDBParser()
        filenumber, filename = tempfile.mkstemp()
        os.close(filenumber)
        try:
            struct = parser.get_structure('MT', filename)
            # Structure has no children (models)
            self.assertFalse(len(struct))
        finally:
            os.remove(filename)

    def test_c_n(self):
        """Extract polypeptides from 1A80."""
        parser = PDBParser(PERMISSIVE=False)
        structure = parser.get_structure("example", "PDB/1A8O.pdb")
        self.assertEqual(len(structure), 1)
        for ppbuild in [PPBuilder(), CaPPBuilder()]:
            # ==========================================================
            # First try allowing non-standard amino acids,
            polypeptides = ppbuild.build_peptides(structure[0], False)
            self.assertEqual(len(polypeptides), 1)
            pp = polypeptides[0]
            # Check the start and end positions
            self.assertEqual(pp[0].get_id()[1], 151)
            self.assertEqual(pp[-1].get_id()[1], 220)
            # Check the sequence
            s = pp.get_sequence()
            self.assertTrue(isinstance(s, Seq))
            self.assertEqual(s.alphabet, generic_protein)
            # Here non-standard MSE are shown as M
            self.assertEqual("MDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQ"
                             "NANPDCKTILKALGPGATLEEMMTACQG", str(s))
            # ==========================================================
            # Now try strict version with only standard amino acids
            # Should ignore MSE 151 at start, and then break the chain
            # at MSE 185, and MSE 214,215
            polypeptides = ppbuild.build_peptides(structure[0], True)
            self.assertEqual(len(polypeptides), 3)
            # First fragment
            pp = polypeptides[0]
            self.assertEqual(pp[0].get_id()[1], 152)
            self.assertEqual(pp[-1].get_id()[1], 184)
            s = pp.get_sequence()
            self.assertTrue(isinstance(s, Seq))
            self.assertEqual(s.alphabet, generic_protein)
            self.assertEqual("DIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNW", str(s))
            # Second fragment
            pp = polypeptides[1]
            self.assertEqual(pp[0].get_id()[1], 186)
            self.assertEqual(pp[-1].get_id()[1], 213)
            s = pp.get_sequence()
            self.assertTrue(isinstance(s, Seq))
            self.assertEqual(s.alphabet, generic_protein)
            self.assertEqual("TETLLVQNANPDCKTILKALGPGATLEE", str(s))
            # Third fragment
            pp = polypeptides[2]
            self.assertEqual(pp[0].get_id()[1], 216)
            self.assertEqual(pp[-1].get_id()[1], 220)
            s = pp.get_sequence()
            self.assertTrue(isinstance(s, Seq))
            self.assertEqual(s.alphabet, generic_protein)
            self.assertEqual("TACQG", str(s))

    def test_strict(self):
        """Parse 1A8O.pdb file in strict mode."""
        parser = PDBParser(PERMISSIVE=False)
        structure = parser.get_structure("example", "PDB/1A8O.pdb")
        self.assertEqual(len(structure), 1)
        model = structure[0]
        self.assertEqual(model.id, 0)
        self.assertEqual(model.level, "M")
        self.assertEqual(len(model), 1)
        chain = model["A"]
        self.assertEqual(chain.id, "A")
        self.assertEqual(chain.level, "C")
        self.assertEqual(len(chain), 158)
        self.assertEqual(" ".join(residue.resname for residue in chain),
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
                         "HOH HOH")
        self.assertEqual(" ".join(atom.name for atom in chain.get_atoms()),
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
                         "O O O O O O O O O O O O O O O O O O O O O O O O")
        self.assertEqual(" ".join(atom.element for atom in chain.get_atoms()),
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
                         "O O O O O O O O O O O O O O O O O O O O O")

    def test_model_numbering(self):
        """Preserve model serial numbers during I/O."""
        def confirm_numbering(struct):
            self.assertEqual(len(struct), 3)
            for idx, model in enumerate(struct):
                self.assertEqual(model.serial_num, idx + 1)
                self.assertEqual(model.serial_num, model.id + 1)

        def confirm_single_end(fname):
            """Ensure there is only one END statement in multi-model files"""
            with open(fname) as handle:
                end_stment = []
                for iline, line in enumerate(handle):
                    if line.strip() == 'END':
                        end_stment.append((line, iline))
            self.assertEqual(len(end_stment), 1)  # Only one?
            self.assertEqual(end_stment[0][1], iline)  # Last line of the file?

        parser = PDBParser(QUIET=1)
        struct1 = parser.get_structure("1lcd", "PDB/1LCD.pdb")
        confirm_numbering(struct1)

        # Round trip: serialize and parse again
        io = PDBIO()
        io.set_structure(struct1)
        filenumber, filename = tempfile.mkstemp()
        os.close(filenumber)
        try:
            io.save(filename)
            struct2 = parser.get_structure("1lcd", filename)
            confirm_numbering(struct2)
            confirm_single_end(filename)
        finally:
            os.remove(filename)


class WriteTest(unittest.TestCase):
    def setUp(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            self.parser = PDBParser(PERMISSIVE=1)
            self.structure = self.parser.get_structure("example", "PDB/1A8O.pdb")

    def test_pdbio_write_structure(self):
        """Write a full structure using PDBIO"""
        io = PDBIO()
        struct1 = self.structure
        # Write full model to temp file
        io.set_structure(struct1)
        filenumber, filename = tempfile.mkstemp()
        os.close(filenumber)
        try:
            io.save(filename)
            struct2 = self.parser.get_structure("1a8o", filename)
            nresidues = len(list(struct2.get_residues()))
            self.assertEqual(len(struct2), 1)
            self.assertEqual(nresidues, 158)
        finally:
            os.remove(filename)

    def test_pdbio_write_residue(self):
        """Write a single residue using PDBIO"""
        io = PDBIO()
        struct1 = self.structure
        residue1 = list(struct1.get_residues())[0]
        # Write full model to temp file
        io.set_structure(residue1)
        filenumber, filename = tempfile.mkstemp()
        os.close(filenumber)
        try:
            io.save(filename)
            struct2 = self.parser.get_structure("1a8o", filename)
            nresidues = len(list(struct2.get_residues()))
            self.assertEqual(nresidues, 1)
        finally:
            os.remove(filename)

    def test_pdbio_write_custom_residue(self):
        """Write a chainless residue using PDBIO"""
        io = PDBIO()

        res = Residue.Residue((' ', 1, ' '), 'DUM', '')
        atm = Atom.Atom('CA', [0.1, 0.1, 0.1], 1.0, 1.0, ' ', 'CA', 1, 'C')
        res.add(atm)

        # Write full model to temp file
        io.set_structure(res)
        filenumber, filename = tempfile.mkstemp()
        os.close(filenumber)
        try:
            io.save(filename)
            struct2 = self.parser.get_structure("res", filename)
            latoms = list(struct2.get_atoms())
            self.assertEqual(len(latoms), 1)
            self.assertEqual(latoms[0].name, 'CA')
            self.assertEqual(latoms[0].parent.resname, 'DUM')
            self.assertEqual(latoms[0].parent.parent.id, 'A')
        finally:
            os.remove(filename)

    def test_pdbio_select(self):
        """Write a selection of the structure using a Select subclass"""
        # Selection class to filter all alpha carbons
        class CAonly(Select):
            """
            Accepts only CA residues
            """
            def accept_atom(self, atom):
                if atom.name == "CA" and atom.element == "C":
                    return 1

        io = PDBIO()
        struct1 = self.structure
        # Write to temp file
        io.set_structure(struct1)
        filenumber, filename = tempfile.mkstemp()
        os.close(filenumber)
        try:
            io.save(filename, CAonly())
            struct2 = self.parser.get_structure("1a8o", filename)
            nresidues = len(list(struct2.get_residues()))
            self.assertEqual(nresidues, 70)
        finally:
            os.remove(filename)

    def test_pdbio_missing_occupancy(self):
        """Write PDB file with missing occupancy"""
        io = PDBIO()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            structure = self.parser.get_structure("test", "PDB/occupancy.pdb")
        io.set_structure(structure)
        filenumber, filename = tempfile.mkstemp()
        os.close(filenumber)
        try:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always", BiopythonWarning)
                io.save(filename)
                self.assertEqual(len(w), 1, w)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", PDBConstructionWarning)
                struct2 = self.parser.get_structure("test", filename)
            atoms = struct2[0]['A'][(' ', 152, ' ')]
            self.assertEqual(atoms['N'].get_occupancy(), None)
        finally:
            os.remove(filename)


class Exposure(unittest.TestCase):
    """Testing Bio.PDB.HSExposure."""
    def setUp(self):
        pdb_filename = "PDB/a_structure.pdb"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            structure = PDBParser(PERMISSIVE=True).get_structure('X', pdb_filename)
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
        hse = HSExposureCA(self.model, self.radius)
        residues = self.a_residues
        self.assertEqual(0, len(residues[0].xtra))
        self.assertEqual(0, len(residues[1].xtra))
        self.assertEqual(3, len(residues[2].xtra))
        self.assertAlmostEqual(0.81250973133184456, residues[2].xtra["EXP_CB_PCB_ANGLE"])
        self.assertEqual(14, residues[2].xtra["EXP_HSE_A_D"])
        self.assertEqual(14, residues[2].xtra["EXP_HSE_A_U"])
        self.assertEqual(3, len(residues[3].xtra))
        self.assertAlmostEqual(1.3383737, residues[3].xtra["EXP_CB_PCB_ANGLE"])
        self.assertEqual(13, residues[3].xtra["EXP_HSE_A_D"])
        self.assertEqual(16, residues[3].xtra["EXP_HSE_A_U"])
        # ...
        self.assertEqual(3, len(residues[-2].xtra))
        self.assertAlmostEqual(0.77124014456278489, residues[-2].xtra["EXP_CB_PCB_ANGLE"])
        self.assertEqual(24, residues[-2].xtra["EXP_HSE_A_D"])
        self.assertEqual(24, residues[-2].xtra["EXP_HSE_A_U"])
        self.assertEqual(0, len(residues[-1].xtra))

    def test_HSExposureCB(self):
        """HSExposureCB."""
        hse = HSExposureCB(self.model, self.radius)
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
        hse = ExposureCN(self.model, self.radius)
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


class Atom_Element(unittest.TestCase):
    """induces Atom Element from Atom Name"""

    def setUp(self):
        pdb_filename = "PDB/a_structure.pdb"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            structure = PDBParser(PERMISSIVE=True).get_structure('X', pdb_filename)
        self.residue = structure[0]['A'][('H_PCA', 1, ' ')]

    def test_AtomElement(self):
        """ Atom Element """
        atoms = self.residue.child_list
        self.assertEqual('N', atoms[0].element)  # N
        self.assertEqual('C', atoms[1].element)  # Alpha Carbon
        self.assertEqual('CA', atoms[8].element)  # Calcium
        self.assertEqual('D', atoms[4].element)  # Deuterium

    def test_ions(self):
        """Element for magnesium is assigned correctly."""
        pdb_filename = "PDB/ions.pdb"
        structure = PDBParser(PERMISSIVE=True).get_structure('X', pdb_filename)
        # check magnesium atom
        atoms = structure[0]['A'][('H_ MG', 1, ' ')].child_list
        self.assertEqual('MG', atoms[0].element)

    def test_hydrogens(self):

        def quick_assign(fullname):
            return Atom.Atom(fullname.strip(), None, None, None, None,
                             fullname, None).element

        pdb_elements = dict(
            H=(' H  ', ' HA ', ' HB ', ' HD1', ' HD2', ' HE ', ' HE1', ' HE2',
               ' HE3', ' HG ', ' HG1', ' HH ', ' HH2', ' HZ ', ' HZ2', ' HZ3',
               '1H  ', '1HA ', '1HB ', '1HD ', '1HD1', '1HD2', '1HE ', '1HE2',
               '1HG ', '1HG1', '1HG2', '1HH1', '1HH2', '1HZ ', '2H  ', '2HA ',
               '2HB ', '2HD ', '2HD1', '2HD2', '2HE ', '2HE2', '2HG ', '2HG1',
               '2HG2', '2HH1', '2HH2', '2HZ ', '3H  ', '3HB ', '3HD1', '3HD2',
               '3HE ', '3HG1', '3HG2', '3HZ ', 'HE21'),
            O=(' OH ',),
            C=(' CH2',),
            N=(' NH1', ' NH2'),
        )

        for element, atom_names in pdb_elements.items():
            for fullname in atom_names:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", PDBConstructionWarning)
                    e = quick_assign(fullname)
                # warnings.warn("%s %s" % (fullname, e))
                self.assertEqual(e, element)


class IterationTests(unittest.TestCase):

    def setUp(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            self.struc = PDBParser(PERMISSIVE=True).get_structure('X', "PDB/a_structure.pdb")

    def test_get_chains(self):
        """Yields chains from different models separately."""
        chains = [chain.id for chain in self.struc.get_chains()]
        self.assertEqual(chains, ['A', 'A', 'B', ' '])

    def test_get_residues(self):
        """Yields all residues from all models."""
        residues = [resi.id for resi in self.struc.get_residues()]
        self.assertEqual(len(residues), 167)

    def test_get_atoms(self):
        """Yields all atoms from the structure, excluding duplicates and ALTLOCs which are not parsed."""
        atoms = ["%12s" % str((atom.id, atom.altloc)) for atom in self.struc.get_atoms()]
        self.assertEqual(len(atoms), 756)


class ChangingIdTests(unittest.TestCase):

    def setUp(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            self.struc = PDBParser(PERMISSIVE=True).get_structure(
                                                  'X', "PDB/a_structure.pdb")

    def test_change_model_id(self):
        """Change the id of a model"""
        for model in self.struc:
            break  # Get first model in structure
        model.id = 2
        self.assertEqual(model.id, 2)
        self.assertIn(2, self.struc)
        self.assertNotIn(0, self.struc)

    def test_change_model_id_raises(self):
        """Cannot change id to a value already in use by another child"""
        model = next(iter(self.struc))
        with self.assertRaises(ValueError):
            model.id = 1
        # Make sure nothing was changed
        self.assertEqual(model.id, 0)
        self.assertIn(0, self.struc)
        self.assertIn(1, self.struc)

    def test_change_chain_id(self):
        """Change the id of a model"""
        chain = next(iter(self.struc.get_chains()))
        chain.id = "R"
        self.assertEqual(chain.id, "R")
        model = next(iter(self.struc))
        self.assertIn("R", model)

    def test_change_residue_id(self):
        """Change the id of a residue"""
        chain = next(iter(self.struc.get_chains()))
        res = chain[('H_PCA', 1, ' ')]
        res.id = (' ', 1, ' ')

        self.assertEqual(res.id, (' ', 1, ' '))
        self.assertIn((' ', 1, ' '), chain)
        self.assertNotIn(('H_PCA', 1, ' '), chain)
        self.assertEqual(chain[(' ', 1, ' ')], res)

    def test_full_id_is_updated_residue(self):
        """
        Invalidate cached full_ids if an id is changed.
        """
        atom = next(iter(self.struc.get_atoms()))

        # Generate the original full id.
        original_id = atom.get_full_id()
        self.assertEqual(original_id,
                         ('X', 0, 'A', ('H_PCA', 1, ' '), ('N', ' ')))
        residue = next(iter(self.struc.get_residues()))

        # Make sure the full id was in fact cached,
        # so we need to invalidate it later.
        self.assertEqual(residue.full_id, ('X', 0, 'A', ('H_PCA', 1, ' ')))

        # Changing the residue's id should lead to an updated full id.
        residue.id = (' ', 1, ' ')
        new_id = atom.get_full_id()
        self.assertNotEqual(original_id, new_id)
        self.assertEqual(new_id, ('X', 0, 'A', (' ', 1, ' '), ('N', ' ')))

    def test_full_id_is_updated_chain(self):
        """
        Invalidate cached full_ids if an id is changed.
        """
        atom = next(iter(self.struc.get_atoms()))

        # Generate the original full id.
        original_id = atom.get_full_id()
        self.assertEqual(original_id,
                         ('X', 0, 'A', ('H_PCA', 1, ' '), ('N', ' ')))
        residue = next(iter(self.struc.get_residues()))

        # Make sure the full id was in fact cached,
        # so we need to invalidate it later.
        self.assertEqual(residue.full_id, ('X', 0, 'A', ('H_PCA', 1, ' ')))
        chain = next(iter(self.struc.get_chains()))

        # Changing the chain's id should lead to an updated full id.
        chain.id = 'Q'
        new_id = atom.get_full_id()
        self.assertNotEqual(original_id, new_id)
        self.assertEqual(new_id, ('X', 0, 'Q', ('H_PCA', 1, ' '), ('N', ' ')))

# class RenumberTests(unittest.TestCase):
#    """Tests renumbering of structures."""
#
#    def setUp(self):
#        pdb_filename = "PDB/1A8O.pdb"
#        self.structure=PDBParser(PERMISSIVE=True).get_structure('X', pdb_filename)
#
#    def test_renumber_residues(self):
#        """Residues in a structure are renumbered."""
#        self.structure.renumber_residues()
#        nums = [resi.id[1] for resi in self.structure[0]['A'].child_list]
#        print(nums)
#
# -------------------------------------------------------------


class TransformTests(unittest.TestCase):

    def setUp(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            self.s = PDBParser(PERMISSIVE=True).get_structure(
                'X', "PDB/a_structure.pdb")
        self.m = self.s.get_list()[0]
        self.c = self.m.get_list()[0]
        self.r = self.c.get_list()[0]
        self.a = self.r.get_list()[0]

    def get_total_pos(self, o):
        """
        Returns the sum of the positions of atoms in an entity along
        with the number of atoms.
        """
        if hasattr(o, "get_coord"):
            return o.get_coord(), 1
        total_pos = numpy.array((0.0, 0.0, 0.0))
        total_count = 0
        for p in o.get_list():
            pos, count = self.get_total_pos(p)
            total_pos += pos
            total_count += count
        return total_pos, total_count

    def get_pos(self, o):
        """
        Returns the average atom position in an entity.
        """
        pos, count = self.get_total_pos(o)
        return 1.0 * pos / count

    def test_transform(self):
        """Transform entities (rotation and translation)."""
        for o in (self.s, self.m, self.c, self.r, self.a):
            rotation = rotmat(Vector(1, 3, 5), Vector(1, 0, 0))
            translation = numpy.array((2.4, 0, 1), 'f')
            oldpos = self.get_pos(o)
            o.transform(rotation, translation)
            newpos = self.get_pos(o)
            newpos_check = numpy.dot(oldpos, rotation) + translation
            for i in range(0, 3):
                self.assertAlmostEqual(newpos[i], newpos_check[i])

    def test_Vector(self):
        """Test Vector object"""
        v1 = Vector(0, 0, 1)
        v2 = Vector(0, 0, 0)
        v3 = Vector(0, 1, 0)
        v4 = Vector(1, 1, 0)

        self.assertEqual(calc_angle(v1, v2, v3), 1.5707963267948966)
        self.assertEqual(calc_dihedral(v1, v2, v3, v4), 1.5707963267948966)
        self.assertTrue(numpy.array_equal((v1 - v2).get_array(), numpy.array([0.0, 0.0, 1.0])))
        self.assertTrue(numpy.array_equal((v1 - 1).get_array(), numpy.array([-1.0, -1.0, 0.0])))
        self.assertTrue(numpy.array_equal((v1 - (1, 2, 3)).get_array(), numpy.array([-1.0, -2.0, -2.0])))
        self.assertTrue(numpy.array_equal((v1 + v2).get_array(), numpy.array([0.0, 0.0, 1.0])))
        self.assertTrue(numpy.array_equal((v1 + 3).get_array(), numpy.array([3.0, 3.0, 4.0])))
        self.assertTrue(numpy.array_equal((v1 + (1, 2, 3)).get_array(), numpy.array([1.0, 2.0, 4.0])))
        self.assertTrue(numpy.array_equal(v1.get_array() / 2, numpy.array([0, 0, 0.5])))
        self.assertTrue(numpy.array_equal(v1.get_array() / 2, numpy.array([0, 0, 0.5])))
        self.assertEqual(v1 * v2, 0.0)
        self.assertTrue(numpy.array_equal((v1 ** v2).get_array(), numpy.array([0.0, -0.0, 0.0])))
        self.assertTrue(numpy.array_equal((v1 ** 2).get_array(), numpy.array([0.0, 0.0, 2.0])))
        self.assertTrue(numpy.array_equal((v1 ** (1, 2, 3)).get_array(), numpy.array([0.0, 0.0, 3.0])))
        self.assertEqual(v1.norm(), 1.0)
        self.assertEqual(v1.normsq(), 1.0)
        v1[2] = 10
        self.assertEqual(v1.__getitem__(2), 10)

        # Vector normalization
        v1 = Vector([2, 0, 0])
        self.assertTrue(numpy.array_equal(v1.normalized().get_array(), numpy.array([1, 0, 0])))
        # State of v1 should not be affected by `normalized`
        self.assertTrue(numpy.array_equal(v1.get_array(), numpy.array([2, 0, 0])))
        v1.normalize()
        # State of v1 should be affected by `normalize`
        self.assertTrue(numpy.array_equal(v1.get_array(), numpy.array([1, 0, 0])))

    def test_refmat(self):
        v1 = Vector(0, 0, 1)
        v2 = Vector(0, 1, 0)
        ref = refmat(v1, v2)
        self.assertTrue(numpy.allclose(ref[0], [1.0, 0.0, 0.0]))
        self.assertTrue(numpy.allclose(ref[1], [0.0, 0.0, 1.0]))
        self.assertTrue(numpy.allclose(ref[2], [0.0, 1.0, 0.0]))
        self.assertTrue(numpy.allclose(v1.left_multiply(ref).get_array(), [0.0, 1.0, 0.0]))

    def test_rotmat(self):
        # Regular 90 deg rotation
        v1 = Vector(0, 0, 1)
        v2 = Vector(0, 1, 0)
        rot = rotmat(v1, v2)
        self.assertTrue(numpy.allclose(rot[0], numpy.array([1.0, 0.0, 0.0])))
        self.assertTrue(numpy.allclose(rot[1], numpy.array([0.0, 0.0, 1.0])))
        self.assertTrue(numpy.allclose(rot[2], numpy.array([0.0, -1.0, 0.0])))
        self.assertTrue(numpy.allclose(v1.left_multiply(rot).get_array(), [0.0, 1.0, 0.0]))
        self.assertTrue(numpy.allclose(v1.right_multiply(numpy.transpose(rot)).get_array(), [0.0, 1.0, 0.0]))

        # Applying rotmat works when the rotation is 180 deg (singularity)
        v1 = Vector([1.0, 0.8, 0])
        v2 = Vector([-1.0, -0.8, 0])
        rot = rotmat(v1, v2)
        v3 = v1.left_multiply(rot)
        self.assertTrue(numpy.allclose(v2.get_array(), v3.get_array()))

        # Applying rotmat works when the rotation is 0 deg (singularity)
        v1 = Vector([1.0, 0.8, 0])
        v2 = Vector([1.0, 0.8, 0])
        rot = rotmat(v1, v2)
        v3 = v1.left_multiply(rot)
        self.assertTrue(numpy.allclose(v1.get_array(), v3.get_array()))

    def test_m2rotaxis(self):
        # Regular 90 deg rotation
        v1 = Vector(0, 0, 1)
        v2 = Vector(0, 1, 0)
        rot = rotmat(v1, v2)
        angle, axis = m2rotaxis(rot)
        self.assertTrue(numpy.allclose(axis.get_array(), [-1.0, 0.0, 0.0]))
        self.assertTrue(abs(angle - numpy.pi / 2) < 1e-5)

        # 180 deg rotation
        v1 = Vector([1.0, 0.8, 0])
        v2 = Vector([-1.0, -0.8, 0])
        rot = rotmat(v1, v2)
        angle, axis = m2rotaxis(rot)
        self.assertTrue(abs(axis * v1) < 1e-5)  # axis orthogonal to v1
        self.assertTrue(abs(angle - numpy.pi) < 1e-5)

        # 0 deg rotation. Axis must be [1, 0, 0] as per Vector documentation
        v1 = Vector([1.0, 0.8, 0])
        v2 = Vector([1.0, 0.8, 0])
        rot = rotmat(v1, v2)
        angle, axis = m2rotaxis(rot)
        self.assertTrue(numpy.allclose(axis.get_array(), [1, 0, 0]))
        self.assertTrue(abs(angle) < 1e-5)

    def test_Vector_angles(self):
        angle = random() * numpy.pi
        axis = Vector(random(3) - random(3))
        axis.normalize()
        m = rotaxis(angle, axis)
        cangle, caxis = m2rotaxis(m)
        self.assertAlmostEqual(angle, cangle, places=3)
        self.assertTrue(numpy.allclose(list(map(int, (axis - caxis).get_array())), [0, 0, 0]),
                        "Want %r and %r to be almost equal" % (axis.get_array(), caxis.get_array()))


class StructureAlignTests(unittest.TestCase):

    def test_StructAlign(self):
        """Tests on module to align two proteins according to a FASTA alignment file."""
        al_file = "PDB/alignment_file.fa"
        pdb2 = "PDB/1A8O.pdb"
        pdb1 = "PDB/2XHE.pdb"
        with open(al_file, 'r') as handle:
            records = AlignIO.read(handle, "fasta")
        p = PDBParser()
        s1 = p.get_structure('1', pdb1)
        p = PDBParser()
        s2 = p.get_structure('2', pdb2)
        m1 = s1[0]
        m2 = s2[0]
        al = StructureAlignment(records, m1, m2)
        self.assertFalse(al.map12 == al.map21)
        self.assertTrue(len(al.map12), 566)
        self.assertTrue(len(al.map21), 70)
        chain1_A = m1["A"]
        chain2_A = m2["A"]
        self.assertEqual(chain1_A[202].get_resname(), 'ILE')
        self.assertEqual(chain2_A[202].get_resname(), 'LEU')
        self.assertEqual(chain1_A[291].get_resname(), chain2_A[180].get_resname())
        self.assertNotEqual(chain1_A[291].get_resname(), chain2_A[181].get_resname())


class SuperimposerTests(unittest.TestCase):

    def test_Superimposer(self):
        """Test on module that superimpose two protein structures."""
        pdb1 = "PDB/1A8O.pdb"
        p = PDBParser()
        s1 = p.get_structure("FIXED", pdb1)
        fixed = Selection.unfold_entities(s1, "A")
        s2 = p.get_structure("MOVING", pdb1)
        moving = Selection.unfold_entities(s2, "A")
        rot = numpy.identity(3).astype('f')
        tran = numpy.array((1.0, 2.0, 3.0), 'f')
        for atom in moving:
            atom.transform(rot, tran)
        sup = Superimposer()
        sup.set_atoms(fixed, moving)
        self.assertTrue(numpy.allclose(sup.rotran[0], numpy.identity(3)))
        self.assertTrue(numpy.allclose(sup.rotran[1], numpy.array([-1.0, -2.0, -3.0])))
        self.assertAlmostEqual(sup.rms, 0.0, places=3)
        atom_list = ['N', 'C', 'C', 'O', 'C', 'C', 'SE', 'C', 'N', 'C', 'C',
                     'O', 'C', 'C', 'O', 'O', 'N', 'C', 'C', 'O', 'C', 'C',
                     'C', 'C', 'N', 'C', 'C', 'O', 'C', 'C', 'C', 'N', 'C',
                     'N', 'N', 'N', 'C', 'C', 'O', 'C', 'C', 'C', 'O', 'N',
                     'N', 'C', 'C', 'O', 'N', 'C', 'C', 'O', 'C', 'C', 'C',
                     'N', 'C', 'C', 'O', 'C', 'C', 'C', 'C', 'N', 'N', 'C',
                     'C', 'O', 'C', 'C', 'C', 'O', 'O', 'N', 'C', 'C', 'O',
                     'C', 'C', 'C', 'N', 'C', 'C', 'O', 'C', 'C', 'C', 'C',
                     'C', 'C', 'C', 'N', 'C', 'C', 'O', 'C', 'C', 'C', 'N',
                     'C', 'N', 'N', 'N', 'C', 'C', 'O', 'C', 'C', 'O', 'O',
                     'N', 'C', 'C', 'O', 'C', 'C', 'C', 'C', 'C', 'C', 'C',
                     'O', 'N', 'C', 'C', 'O', 'C', 'C', 'C', 'N', 'C', 'C',
                     'O', 'C', 'C', 'O', 'O', 'N', 'C', 'C', 'O', 'C', 'C',
                     'C', 'N', 'C', 'N', 'N', 'N', 'C', 'C', 'O', 'C', 'C',
                     'C', 'C', 'C', 'C', 'C', 'N', 'C', 'C', 'O', 'C', 'C',
                     'C', 'C', 'C', 'C', 'C', 'O', 'N', 'C', 'C', 'O', 'C',
                     'C', 'C', 'C', 'N', 'N', 'C', 'C', 'O', 'C', 'O', 'C',
                     'N', 'C', 'C', 'O', 'C', 'C', 'C', 'C', 'N', 'C', 'C',
                     'O', 'C', 'C', 'C', 'N', 'C', 'N', 'N', 'N', 'C', 'C',
                     'O', 'C', 'N', 'C', 'C', 'O', 'C', 'C', 'C', 'O', 'O',
                     'N', 'C', 'C', 'O', 'C', 'C', 'C', 'O', 'N', 'N', 'C',
                     'C', 'O', 'C', 'N', 'C', 'C', 'O', 'C', 'O', 'N', 'C',
                     'C', 'O', 'C', 'C', 'C', 'O', 'N', 'N', 'C', 'C', 'O',
                     'C', 'C', 'C', 'O', 'O', 'N', 'C', 'C', 'O', 'C', 'C',
                     'C', 'N', 'C', 'C', 'O', 'C', 'C', 'C', 'C', 'N', 'N',
                     'C', 'C', 'O', 'C', 'C', 'O', 'N', 'N', 'C', 'C', 'O',
                     'C', 'C', 'C', 'C', 'N', 'C', 'C', 'C', 'C', 'C', 'N',
                     'C', 'C', 'O', 'C', 'C', 'SE', 'C', 'N', 'C', 'C', 'O',
                     'C', 'O', 'C', 'N', 'C', 'C', 'O', 'C', 'C', 'C', 'O',
                     'O', 'N', 'C', 'C', 'O', 'C', 'O', 'C', 'N', 'C', 'C',
                     'O', 'C', 'C', 'C', 'C', 'N', 'C', 'C', 'O', 'C', 'C',
                     'C', 'C', 'N', 'C', 'C', 'O', 'C', 'C', 'C', 'N', 'C',
                     'C', 'O', 'C', 'C', 'C', 'O', 'N', 'N', 'C', 'C', 'O',
                     'C', 'C', 'O', 'N', 'N', 'C', 'C', 'O', 'C', 'N', 'C',
                     'C', 'O', 'C', 'C', 'O', 'N', 'N', 'C', 'C', 'O', 'C',
                     'C', 'C', 'N', 'C', 'C', 'O', 'C', 'C', 'O', 'O', 'N',
                     'C', 'C', 'O', 'C', 'S', 'N', 'C', 'C', 'O', 'C', 'C',
                     'C', 'C', 'N', 'N', 'C', 'C', 'O', 'C', 'O', 'C', 'N',
                     'C', 'C', 'O', 'C', 'C', 'C', 'C', 'N', 'C', 'C', 'O',
                     'C', 'C', 'C', 'C', 'N', 'C', 'C', 'O', 'C', 'C', 'C',
                     'C', 'N', 'N', 'C', 'C', 'O', 'C', 'N', 'C', 'C', 'O',
                     'C', 'C', 'C', 'C', 'N', 'C', 'C', 'O', 'N', 'C', 'C',
                     'O', 'C', 'C', 'C', 'N', 'C', 'C', 'O', 'N', 'C', 'C',
                     'O', 'C', 'N', 'C', 'C', 'O', 'C', 'O', 'C', 'N', 'C',
                     'C', 'O', 'C', 'C', 'C', 'C', 'N', 'C', 'C', 'O', 'C',
                     'C', 'C', 'O', 'O', 'N', 'C', 'C', 'O', 'C', 'C', 'C',
                     'O', 'O', 'N', 'C', 'C', 'O', 'C', 'C', 'SE', 'C', 'N',
                     'C', 'C', 'O', 'C', 'C', 'SE', 'C', 'N', 'C', 'C', 'O',
                     'C', 'O', 'C', 'N', 'C', 'C', 'O', 'C', 'N', 'C', 'C',
                     'O', 'C', 'S', 'N', 'C', 'C', 'O', 'C', 'C', 'C', 'O',
                     'N', 'N', 'C', 'C', 'O', 'O', 'O', 'O', 'O', 'O', 'O',
                     'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O',
                     'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O',
                     'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O',
                     'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O',
                     'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O',
                     'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O',
                     'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O',
                     'O', 'O', 'O', 'O', 'O', 'O']
        sup.apply(moving)
        atom_moved = []
        for aa in moving:
            atom_moved.append(aa.element)
        self.assertEqual(atom_moved, atom_list)


class PolypeptideTests(unittest.TestCase):

    def test_polypeptide(self):
        """Tests on polypetide class and methods."""
        p = PDBParser(PERMISSIVE=True)
        pdb1 = "PDB/1A8O.pdb"
        s = p.get_structure("scr", pdb1)
        ppb = PPBuilder()
        pp = ppb.build_peptides(s)
        self.assertEqual(str(pp[0].get_sequence()), "DIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNW")
        self.assertEqual(str(pp[1].get_sequence()), "TETLLVQNANPDCKTILKALGPGATLEE")
        self.assertEqual(str(pp[2].get_sequence()), "TACQG")
        phi_psi = pp[0].get_phi_psi_list()
        self.assertEqual(phi_psi[0][0], None)
        self.assertAlmostEqual(phi_psi[0][1], -0.46297171497725553, places=3)
        self.assertAlmostEqual(phi_psi[1][0], -1.0873937604007962, places=3)
        self.assertAlmostEqual(phi_psi[1][1], 2.1337707832637109, places=3)
        self.assertAlmostEqual(phi_psi[2][0], -2.4052232743651878, places=3)
        self.assertAlmostEqual(phi_psi[2][1], 2.3807316946081554, places=3)
        phi_psi = pp[1].get_phi_psi_list()
        self.assertEqual(phi_psi[0][0], None)
        self.assertAlmostEqual(phi_psi[0][1], -0.6810077089092923, places=3)
        self.assertAlmostEqual(phi_psi[1][0], -1.2654003477656888, places=3)
        self.assertAlmostEqual(phi_psi[1][1], -0.58689987042756309, places=3)
        self.assertAlmostEqual(phi_psi[2][0], -1.7467679151684763, places=3)
        self.assertAlmostEqual(phi_psi[2][1], -1.5655066256698336, places=3)
        phi_psi = pp[2].get_phi_psi_list()
        self.assertEqual(phi_psi[0][0], None)
        self.assertAlmostEqual(phi_psi[0][1], -0.73222884210889716, places=3)
        self.assertAlmostEqual(phi_psi[1][0], -1.1044740234566259, places=3)
        self.assertAlmostEqual(phi_psi[1][1], -0.69681334592782884, places=3)
        self.assertAlmostEqual(phi_psi[2][0], -1.8497413300164958, places=3)
        self.assertAlmostEqual(phi_psi[2][1], 0.34762889834809058, places=3)
        ppb = CaPPBuilder()
        pp = ppb.build_peptides(s)
        self.assertEqual(str(pp[0].get_sequence()), "DIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNW")
        self.assertEqual(str(pp[1].get_sequence()), "TETLLVQNANPDCKTILKALGPGATLEE")
        self.assertEqual(str(pp[2].get_sequence()), "TACQG")
        self.assertEqual([ca.serial_number for ca in pp[0].get_ca_list()], [10, 18, 26, 37, 46, 50, 57, 66, 75, 82, 93, 104, 112, 124, 131, 139, 150, 161, 173, 182, 189, 197, 208, 213, 222, 231, 236, 242, 251, 260, 267, 276, 284])
        taus = pp[1].get_tau_list()
        self.assertAlmostEqual(taus[0], 0.3597907225123525, places=3)
        self.assertAlmostEqual(taus[1], 0.43239284636769254, places=3)
        self.assertAlmostEqual(taus[2], 0.99820157492712114, places=3)
        thetas = pp[2].get_theta_list()
        self.assertAlmostEqual(thetas[0], 1.6610069445335354, places=3)
        self.assertAlmostEqual(thetas[1], 1.7491703334817772, places=3)
        self.assertAlmostEqual(thetas[2], 2.0702447422720143, places=3)


class MMCIF2dictTests(unittest.TestCase):

    def test_MMCIF2dict(self):
        filename = "PDB/1A8O.cif"
        mmcif = MMCIF2Dict(filename)
        self.assertEqual(len(mmcif.keys()), 575)
        self.assertEqual(mmcif['_entity_poly_seq.mon_id'], ['MSE', 'ASP', 'ILE', 'ARG', 'GLN', 'GLY', 'PRO', 'LYS', 'GLU', 'PRO', 'PHE', 'ARG', 'ASP', 'TYR', 'VAL', 'ASP', 'ARG', 'PHE', 'TYR', 'LYS', 'THR', 'LEU', 'ARG', 'ALA', 'GLU', 'GLN', 'ALA', 'SER', 'GLN', 'GLU', 'VAL', 'LYS', 'ASN', 'TRP', 'MSE', 'THR', 'GLU', 'THR', 'LEU', 'LEU', 'VAL', 'GLN', 'ASN', 'ALA', 'ASN', 'PRO', 'ASP', 'CYS', 'LYS', 'THR', 'ILE', 'LEU', 'LYS', 'ALA', 'LEU', 'GLY', 'PRO', 'GLY', 'ALA', 'THR', 'LEU', 'GLU', 'GLU', 'MSE', 'MSE', 'THR', 'ALA', 'CYS', 'GLN', 'GLY'])
        self.assertEqual(mmcif['_atom_site.Cartn_x'], ['19.594', '20.255', '20.351', '19.362', '19.457', '20.022', '21.718', '21.424', '21.554', '21.835', '21.947', '21.678', '23.126', '23.098', '23.433', '22.749', '22.322', '22.498', '21.220', '20.214', '23.062', '24.282', '23.423', '25.429', '21.280', '20.173', '20.766', '21.804', '19.444', '18.724', '18.011', '17.416', '16.221', '15.459', '15.824', '20.116', '20.613', '20.546', '19.488', '19.837', '20.385', '19.526', '18.365', '20.090', '21.675', '21.698', '20.859', '20.729', '20.260', '19.435', '20.158', '19.512', '18.993', '20.056', '20.300', '21.486', '22.285', '23.286', '24.155', '23.025', '22.117', '21.236', '20.159', '19.231', '23.152', '24.037', '23.563', '22.398', '24.086', '25.003', '24.858', '23.861', '25.748', '24.459', '24.089', '23.580', '24.111', '25.415', '26.116', '25.852', '22.544', '21.960', '22.965', '22.928', '20.793', '19.999', '19.234', '20.019', '18.495', '19.286', '18.523', '23.861', '24.870', '25.788', '26.158', '25.684', '26.777', '26.215', '27.235', '28.136', '28.155', '29.030', '26.137', '26.994', '26.279', '26.880', '27.408', '28.345', '28.814', '28.620', '24.992', '24.151', '24.025', '24.139', '22.787', '21.629', '21.657', '20.489', '20.571', '19.408', '19.450', '18.365', '23.839', '23.720', '24.962', '24.853', '23.502', '23.661', '22.120', '26.137', '27.387', '27.511', '27.925', '28.595', '28.723', '28.016', '29.545', '27.136', '27.202', '26.238', '26.585', '26.850', '27.835', '27.667', '26.352', '25.494', '25.797', '24.325', '25.037', '23.984', '24.456', '24.305', '22.761', '21.538', '21.301', '20.586', '20.130', '19.415', '19.186', '25.033', '25.526', '26.755', '27.015', '25.771', '24.608', '23.508', '24.583', '22.406', '23.490', '22.406', '21.326', '27.508', '28.691', '28.183', '28.705', '29.455', '30.787', '31.428', '32.618', '33.153', '27.116', '26.508', '25.826', '25.827', '25.475', '26.150', '24.741', '25.264', '24.587', '25.587', '25.302', '23.789', '22.707', '21.787', '21.910', '26.767', '27.806', '28.299', '28.656', '29.006', '28.944', '30.295', '30.744', '30.326', '29.441', '30.787', '28.332', '28.789', '27.943', '28.374', '28.803', '26.740', '25.833', '25.775', '24.998', '24.425', '24.354', '24.816', '24.535', '25.454', '26.601', '26.645', '25.240', '24.885', '27.391', '28.884', '29.200', '28.729', '29.998', '24.438', '23.066', '23.001', '23.824', '22.370', '22.035', '21.831', '21.174', '20.852', '20.917', '19.638', '20.949', '20.315', '18.908', '18.539', '20.262', '19.688', '20.414', '21.592', '19.714', '18.136', '16.775', '16.738', '15.875', '16.101', '15.478', '14.341', '13.247', '14.542', '17.668', '17.730', '18.064', '17.491', '18.754', '18.932', '18.279', '18.971', '19.343', '18.126', '17.905', '20.444', '21.777', '22.756', '24.069', '24.913', '17.344', '16.136', '15.146', '14.599', '15.468', '16.242', '17.164', '15.865', '14.932', '14.017', '14.495', '13.700', '13.904', '13.254', '12.332', '13.484', '11.975', '12.666', '14.303', '12.641', '14.280', '13.452', '15.793', '16.368', '16.285', '16.053', '17.815', '17.939', '17.221', '18.427', '16.438', '16.375', '14.950', '14.778', '16.869', '18.228', '16.791', '13.947', '12.529', '12.045', '11.151', '11.625', '11.950', '11.054', '11.086', '10.326', '12.589', '12.177', '13.076', '12.888', '11.978', '13.202', '10.883', '14.054', '14.963', '15.702', '15.846', '15.935', '15.286', '16.327', '14.580', '16.162', '16.876', '15.961', '16.391', '17.402', '18.238', '19.553', '18.506', '14.695', '13.703', '13.270', '13.262', '12.460', '11.372', '12.854', '12.954', '12.503', '13.541', '13.184', '12.008', '10.830', '10.505', '10.626', '10.093', '14.820', '15.887', '16.443', '17.416', '17.014', '16.627', '15.451', '17.619', '15.830', '16.248', '15.758', '14.809', '15.689', '16.404', '16.005', '14.639', '14.122', '17.109', '17.396', '16.559', '18.588', '14.018', '12.706', '12.516', '11.536', '12.617', '13.288', '14.522', '13.454', '13.383', '13.351', '12.406', '14.564', '14.482', '13.353', '15.552', '14.378', '14.488', '13.443', '12.968', '15.902', '16.144', '13.061', '12.087', '10.746', '10.157', '11.879', '11.014', '11.003', '10.171', '10.269', '10.273', '9.002', '9.101', '8.227', '8.612', '8.611', '7.224', '10.191', '10.458', '10.518', '9.916', '11.791', '11.677', '12.184', '12.967', '11.222', '11.377', '10.082', '9.885', '12.416', '13.824', '14.764', '14.287', '9.214', '7.937', '7.048', '6.294', '7.230', '7.828', '7.618', '8.090', '7.916', '7.189', '6.419', '6.871', '6.391', '6.449', '7.815', '8.305', '7.481', '7.371', '9.788', '10.832', '12.217', '10.789', '6.886', '6.080', '6.922', '8.149', '6.294', '7.024', '7.912', '7.680', '5.901', '4.734', '4.839', '8.952', '9.861', '10.886', '11.642', '10.910', '11.884', '13.285', '13.524', '11.599', '14.199', '15.563', '16.391', '16.022', '16.290', '16.498', '15.473', '17.509', '18.426', '18.875', '19.012', '19.645', '20.773', '20.264', '21.920', '19.082', '19.510', '18.471', '18.816', '19.784', '21.035', '20.954', '19.902', '21.955', '17.199', '16.109', '16.001', '15.690', '14.787', '14.776', '13.539', '13.220', '12.888', '16.301', '16.274', '17.413', '17.209', '16.429', '15.284', '15.332', '13.844', '18.606', '19.764', '19.548', '19.922', '21.047', '21.507', '23.105', '22.645', '18.915', '18.636', '17.640', '17.807', '18.050', '18.998', '17.730', '16.631', '15.593', '16.104', '15.685', '14.486', '17.033', '17.572', '18.985', '19.634', '17.525', '15.855', '19.451', '20.802', '21.001', '20.066', '21.152', '20.421', '20.725', '21.768', '19.817', '22.226', '22.536', '23.683', '24.328', '23.949', '15.165', '19.774', '22.152', '12.938', '23.499', '17.568', '13.544', '15.524', '31.249', '11.999', '14.511', '7.439', '19.303', '17.114', '21.867', '17.573', '26.151', '20.974', '20.796', '28.370', '29.565', '21.248', '25.744', '8.691', '30.789', '30.905', '28.623', '24.935', '23.462', '9.924', '28.729', '13.579', '23.652', '25.631', '17.799', '23.547', '16.363', '24.125', '33.063', '29.209', '10.391', '12.221', '18.997', '16.360', '27.915', '28.158', '21.975', '27.069', '30.148', '21.196', '8.864', '13.228', '18.577', '20.526', '25.758', '7.838', '20.569', '13.009', '19.229', '17.655', '30.445', '9.014', '3.398', '31.603', '16.543', '12.037', '7.261', '5.607', '23.532', '30.701', '32.300', '34.351', '9.450', '29.476', '13.681', '26.728', '10.004', '30.553', '23.569', '10.927', '17.983', '8.191', '32.095', '11.520', '13.249', '15.919', '11.187', '16.743'])
        self.assertEqual(mmcif['_struct_ref.pdbx_seq_one_letter_code'], 'GARASVLSGGELDKWEKIRLRPGGKKQYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYNTIAVLYCVHQRIDVKDTKEALDKIEEEQNKSKKKAQQAAADTGNNSQVSQNYPIVQNLQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEAAEWDRLHPVHAGPIAPGQMREPRGSDIAGTTSTLQEQIGWMTHNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPGATLEEMMTACQGVGGPGHKARVLAEAMSQVTNPATIMIQKGNFRNQRKTVKCFNCGKEGHIAKNCRAPRKKGCWKCGKEGHQMKDCTERQANFLGKIWPSHKGRPGNFLQSRPEPTAPPEESFRFGEETTTPSQKQEPIDKELYPLASLRSLFGSDPSSQ')


class CopyTests(unittest.TestCase):

    def setUp(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            self.s = PDBParser(PERMISSIVE=True).get_structure(
                'X', "PDB/a_structure.pdb")
        self.m = self.s.get_list()[0]
        self.c = self.m.get_list()[0]
        self.r = self.c.get_list()[0]
        self.a = self.r.get_list()[0]

    def test_atom_copy(self):
        aa = self.a.copy()
        self.assertFalse(self.a is aa)
        self.assertFalse(self.a.get_coord() is aa.get_coord())

    def test_entitity_copy(self):
        """Make a copy of a residue."""
        for e in (self.s, self.m, self.c, self.r):
            ee = e.copy()
            self.assertFalse(e is ee)
            self.assertFalse(e.get_list()[0] is ee.get_list()[0])


def eprint(*args, **kwargs):
    """Helper function that prints to stderr."""
    print(*args, file=sys.stderr, **kwargs)


def will_it_float(s):
    """ Helper function that converts the input into a float if it is a number.

    If the input is a string, the output does not change.
    """
    try:
        return float(s)
    except ValueError:
        return(s)


class DsspTests(unittest.TestCase):
    """Tests for DSSP parsing etc which don't need the binary tool.

    See also test_DSSP_tool.py for run time testing with the tool.
    """
    def test_DSSP_file(self):
        """Test parsing of pregenerated DSSP"""
        dssp, keys = make_dssp_dict("PDB/2BEG.dssp")
        self.assertEqual(len(dssp), 130)

    def test_DSSP_noheader_file(self):
        """Test parsing of pregenerated DSSP missing header information"""
        # New DSSP prints a line containing only whitespace and "."
        dssp, keys = make_dssp_dict("PDB/2BEG_noheader.dssp")
        self.assertEqual(len(dssp), 130)

    def test_DSSP_hbonds(self):
        """Test parsing of DSSP hydrogen bond information."""
        dssp, keys = make_dssp_dict("PDB/2BEG.dssp")

        dssp_indices = set(v[5] for v in dssp.values())
        hb_indices = set()

        # The integers preceding each hydrogen bond energy (kcal/mol) in the
        # "N-H-->O    O-->H-N    N-H-->O    O-->H-N" dssp output columns are
        # relative dssp indices. Therefore, "hb_indices" contains the absolute
        # dssp indices of residues participating in (provisional) h-bonds. Note
        # that actual h-bonds are typically determined by an energetic
        # threshold.
        for val in dssp.values():
            hb_indices |= set(
                (val[5] + x) for x in (val[6], val[8], val[10], val[12]))

        # Check if all h-bond partner indices were successfully parsed.
        self.assertEqual((dssp_indices & hb_indices), hb_indices)

    def test_DSSP_in_model_obj(self):
        """ Test that all the elements are added correctly to the xtra attribute of the input model object."""
        p = PDBParser()
        s = p.get_structure("example", "PDB/2BEG.pdb")
        m = s[0]
        # Read the DSSP data into the pdb object:
        trash_var = DSSP(m, "PDB/2BEG.dssp", 'dssp', 'Sander', 'DSSP')
        # Now compare the xtra attribute of the pdb object
        # residue by residue with the pre-computed values:
        i = 0
        with open("PDB/dssp_xtra_Sander.txt", 'r') as fh_ref:
            ref_lines = fh_ref.readlines()
            for chain in m:
                for res in chain:
                    # Split the pre-computed values into a list:
                    xtra_list_ref = ref_lines[i].rstrip().split('\t')
                    # Then convert each element to float where possible:
                    xtra_list_ref = list(map(will_it_float, xtra_list_ref))
                    # The xtra attribute is a dict.
                    # To compare with the pre-comouted values first sort according to keys:
                    xtra_itemts = sorted(res.xtra.items(), key=lambda s: s[0])
                    # Then extract the list of xtra values for the residue
                    # and convert to floats where possible:
                    xtra_list = [t[1] for t in xtra_itemts]
                    xtra_list = list(map(will_it_float, xtra_list))
                    # The reason for converting to float is, that casting a float to a string in python2.6
                    # will include fewer decimals than python3 and an assertion error will be thrown.
                    self.assertEqual(xtra_list, xtra_list_ref)
                    i += 1

    def test_DSSP_RSA(self):
        """Tests the usage of different ASA tables."""
        # Tests include Sander/default, Wilke and Miller
        p = PDBParser()
        # Sander/default:
        s = p.get_structure("example", "PDB/2BEG.pdb")
        m = s[0]
        # Read the DSSP data into the pdb object:
        trash_var = DSSP(m, "PDB/2BEG.dssp", 'dssp', 'Sander', 'DSSP')
        # Then compare the RASA values for each residue with the pre-computed values:
        i = 0
        with open("PDB/Sander_RASA.txt", 'r') as fh_ref:
            ref_lines = fh_ref.readlines()
            for chain in m:
                for res in chain:
                    rasa_ref = float(ref_lines[i].rstrip())
                    rasa = float(res.xtra['EXP_DSSP_RASA'])
                    self.assertAlmostEqual(rasa, rasa_ref)
                    i += 1

        # Wilke (procedure similar as for the Sander values above):
        s = p.get_structure("example", "PDB/2BEG.pdb")
        m = s[0]
        trash_var = DSSP(m, "PDB/2BEG.dssp", 'dssp', 'Wilke', 'DSSP')
        i = 0
        with open("PDB/Wilke_RASA.txt", 'r') as fh_ref:
            ref_lines = fh_ref.readlines()
            for chain in m:
                for res in chain:
                    rasa_ref = float(ref_lines[i].rstrip())
                    rasa = float(res.xtra['EXP_DSSP_RASA'])
                    self.assertAlmostEqual(rasa, rasa_ref)
                    i += 1

        # Miller (procedure similar as for the Sander values above):
        s = p.get_structure("example", "PDB/2BEG.pdb")
        m = s[0]
        trash_var = DSSP(m, "PDB/2BEG.dssp", 'dssp', 'Miller', 'DSSP')
        i = 0
        with open("PDB/Miller_RASA.txt", 'r') as fh_ref:
            ref_lines = fh_ref.readlines()
            for chain in m:
                for res in chain:
                    rasa_ref = float(ref_lines[i].rstrip())
                    rasa = float(res.xtra['EXP_DSSP_RASA'])
                    self.assertAlmostEqual(rasa, rasa_ref)
                    i += 1


class NACCESSTests(unittest.TestCase):
    """Tests for NACCESS parsing etc which don't need the binary tool.

    See also test_NACCESS_tool.py for run time testing with the tool.
    """
    def test_NACCESS_rsa_file(self):
        """Test parsing of pregenerated rsa NACCESS file"""
        with open("PDB/1A8O.rsa") as rsa:
            naccess = process_rsa_data(rsa)
        self.assertEqual(len(naccess), 66)

    def test_NACCESS_asa_file(self):
        """Test parsing of pregenerated asa NACCESS file"""
        with open("PDB/1A8O.asa") as asa:
            naccess = process_asa_data(asa)
        self.assertEqual(len(naccess), 524)


class ResidueDepthTests(unittest.TestCase):
    """Tests for ResidueDepth module, except for running MSMS itself."""
    def test_pdb_to_xyzr(self):
        """Test generation of xyzr (atomic radii) file"""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            p = PDBParser(PERMISSIVE=1)
            structure = p.get_structure("example", "PDB/1A8O.pdb")

        # Read radii produced with original shell script
        with open('PDB/1A8O.xyzr') as handle:
            msms_radii = []
            for line in handle:
                fields = line.split()
                radius = float(fields[3])
                msms_radii.append(radius)

        model = structure[0]
        biopy_radii = []
        for atom in model.get_atoms():
            biopy_radii.append(_get_atom_radius(atom, rtype='united'))

        assert len(msms_radii) == len(biopy_radii)
        self.assertSequenceEqual(msms_radii, biopy_radii)


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
