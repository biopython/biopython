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
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein
from Bio.PDB import PDBParser, PPBuilder, CaPPBuilder, PDBIO, Select, MMCIFParser, MMCIFIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB import HSExposureCA, HSExposureCB, ExposureCN
from Bio.PDB.PDBExceptions import PDBConstructionException, PDBConstructionWarning
from Bio.PDB import rotmat, Vector, refmat, calc_angle, calc_dihedral, rotaxis, m2rotaxis
from Bio.PDB import Residue, Atom
from Bio.PDB import make_dssp_dict
from Bio.PDB import DSSP
from Bio.PDB.NACCESS import process_asa_data, process_rsa_data
from Bio.PDB.ResidueDepth import _get_atom_radius


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
                    'Atom O defined twice in residue <Residue HOH het=W resseq=67 icode= > at line 902.'
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
        self.assertEqual(len(polypeptides), 2)
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
        self.assertEqual(len(polypeptides), 2)
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
        self.assertEqual(len(m1), 4)
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
            ('B', 11, [((' ', 44, ' '), 11),
                       (('H_SEP', 45, ' '), 10),  # Phosphoserine
                       ((' ', 46, ' '), 8),
                       ((' ', 47, ' '), 10),
                       ((' ', 48, ' '), 11),
                       ((' ', 49, ' '), 6),
                       ((' ', 50, ' '), 4),
                       ((' ', 51, ' '), 5),
                       ((' ', 51, 'A'), 5),
                       ((' ', 52, ' '), 7),
                       (('W', 0, ' '), 1),
                       ]),
            ('C', 5, [(('W', 0, ' '), 1),
                      (('H_NAG', 1, ' '), 14),
                      (('H_NAG', 2, ' '), 14),
                      (('H_NAG', 4, ' '), 14),
                      (('H_NAG', 3, ' '), 14),
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
        self.assertEqual(len(model), 4)
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
                         "C O CB CG CD NE CZ NH1 NH2 N CA C O CB SG")

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
                         "C C C N C C O C S N C C O N C C O N C C O C N C O C "
                         "C O N N C C O C C O O N C C O C C C C C C C O N C C "
                         "O C S N C C O C O N C C O N C C O N C C O C C O N N "
                         "C C O C S N C C O C C C O N N C C O C C C C C C C O "
                         "N C C O C C C N C N N N C C O C S")

    def test_pdbio_write_truncated(self):
        """Test parsing of truncated lines."""
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

    # Tests for sorting methods
    def test_comparison_entities(self):
        """Test comparing and sorting the several SMCRA objects"""

        struct = self.structure
        # Test deepcopy of a structure with disordered atoms
        struct2 = deepcopy(struct)

        # Sorting (<, >, <=, <=)
        # Chains (same code as models)
        model = struct[1]
        chains = [c.id for c in sorted(model)]
        self.assertEqual(chains, ['A', 'B', 'C', ' '])
        # Residues
        residues = [r.id[1] for r in sorted(struct[1]['C'])]
        self.assertEqual(residues, [1, 2, 3, 4, 0])
        # Atoms
        for residue in struct.get_residues():
            old = [a.name for a in residue]
            new = [a.name for a in sorted(residue)]

            special = [a for a in ('N', 'CA', 'C', 'O') if a in old]
            len_special = len(special)
            # Placed N, CA, C, O first?
            self.assertEqual(new[:len_special], special,
                             "Sorted residue did not place N, CA, C, "
                             "O first: %s" % new)
            # Placed everyone else alphabetically?
            self.assertEqual(new[len_special:], sorted(new[len_special:]),
                             "After N, CA, C, O order Should be "
                             "alphabetical: %s" % new)
        # DisorderedResidue
        residues = [r.id[1] for r in sorted(struct[1]['A'])][79:81]
        self.assertEqual(residues, [80, 81])
        # Insertion code + hetflag + chain
        residues = [r for r in struct[1]['B']] + [struct[1]['A'][44]]
        self.assertEqual([("{}" * 4).format(r.parent.id, *r.id) for r in sorted(residues)],
                         ['A 44 ', 'B 44 ', 'B 46 ', 'B 47 ', 'B 48 ', 'B 49 ', 'B 50 ',
                          'B 51 ', 'B 51A', 'B 52 ', 'BH_SEP45 ', 'BW0 '])
        # DisorderedAtom
        atoms = [a.altloc for a in sorted(struct[1]['A'][74]['OD1'])]
        self.assertEqual(atoms, ['A', 'B'])

        # Comparisons
        # Structure
        self.assertEqual(struct, struct2)
        self.assertLessEqual(struct, struct2)
        self.assertGreaterEqual(struct, struct2)
        struct2.id = 'new_id'
        self.assertNotEqual(struct, struct2)
        self.assertLess(struct, struct2)
        self.assertLessEqual(struct, struct2)
        self.assertGreater(struct2, struct)
        self.assertGreaterEqual(struct2, struct)

        # Model
        self.assertTrue(model == model)  # __eq__ same type
        self.assertFalse(struct[0] == struct[1])

        self.assertFalse(struct[0] == [])  # __eq__ diff. types
        self.assertFalse(struct == model)

        # residues with same ID string should not be equal if the parent is not equal
        res1, res2, res3 = residues[0], residues[-1], struct2[1]['A'][44]
        self.assertEqual(res1.id, res2.id)
        self.assertEqual(res2, res3)  # Equality of identical residues with different structure ID
        self.assertFalse(res1 == res2)
        self.assertGreater(res1, res2)
        self.assertGreaterEqual(res1, res2)
        self.assertLess(res2, res1)
        self.assertLessEqual(res2, res1)

        # atom should not be equal if the parent is not equal
        atom1, atom2, atom3 = res1['CA'], res2['CA'], res3['CA']
        self.assertEqual(atom2, atom3)  # Equality of identical atoms with different structure ID
        self.assertGreater(atom1, atom2)
        self.assertGreaterEqual(atom1, atom2)
        self.assertGreaterEqual(atom2, atom3)
        self.assertNotEqual(atom1, atom2)
        self.assertLess(atom2, atom1)
        self.assertLessEqual(atom2, atom1)
        self.assertLessEqual(atom2, atom3)

        # In Py2 this will be True/False, in Py3 it will raise a TypeError.
        try:
            self.assertTrue(atom1 < res1)  # __gt__ diff. types
        except TypeError:
            pass

        try:
            self.assertTrue(struct > model)  # __gt__ diff. types
        except TypeError:
            pass

        try:
            self.assertFalse(struct >= [])  # __le__ diff. types
        except TypeError:
            pass


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

    def test_residue_sort(self):
        """Sorting atoms in residues."""
        parser = PDBParser(PERMISSIVE=False)
        structure = parser.get_structure("example", "PDB/1A8O.pdb")
        for residue in structure.get_residues():
            old = [a.name for a in residue]
            new = [a.name for a in sorted(residue)]
            special = []
            for a in ['N', 'CA', 'C', 'O']:
                if a in old:
                    special.append(a)
            special_len = len(special)
            self.assertEqual(new[0:special_len], special,
                             "Sorted residue did not place N, CA, C, O first: %s" % new)
            self.assertEqual(new[special_len:], sorted(new[special_len:]),
                             "After N, CA, C, O should be alphabet: %s" % new)

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
            """Ensure there is only one END statement in multi-model files."""
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
            self.mmcif_parser = MMCIFParser()
            self.structure = self.parser.get_structure("example", "PDB/1A8O.pdb")
            self.mmcif_file = "PDB/1A8O.cif"
            self.mmcif_multimodel_pdb_file = "PDB/1SSU_mod.pdb"
            self.mmcif_multimodel_mmcif_file = "PDB/1SSU_mod.cif"

    def test_pdbio_write_structure(self):
        """Write a full structure using PDBIO."""
        io = PDBIO()
        struct1 = self.structure
        # Ensure that set_structure doesn't alter parent
        parent = struct1.parent
        # Write full model to temp file
        io.set_structure(struct1)
        self.assertIs(parent, struct1.parent)
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
        # Ensure that set_structure doesn't alter parent
        parent = residue1.parent
        # Write full model to temp file
        io.set_structure(residue1)
        self.assertIs(parent, residue1.parent)
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
        """Write a chainless residue using PDBIO."""
        io = PDBIO()

        res = Residue.Residue((' ', 1, ' '), 'DUM', '')
        atm = Atom.Atom('CA', [0.1, 0.1, 0.1], 1.0, 1.0, ' ', 'CA', 1, 'C')
        res.add(atm)

        # Ensure that set_structure doesn't alter parent
        parent = res.parent
        # Write full model to temp file
        io.set_structure(res)
        self.assertIs(parent, res.parent)
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
        """Write a selection of the structure using a Select subclass."""
        # Selection class to filter all alpha carbons
        class CAonly(Select):
            """Accepts only CA residues."""

            def accept_atom(self, atom):
                if atom.name == "CA" and atom.element == "C":
                    return 1

        io = PDBIO()
        struct1 = self.structure
        # Ensure that set_structure doesn't alter parent
        parent = struct1.parent
        # Write to temp file
        io.set_structure(struct1)
        self.assertIs(parent, struct1.parent)
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
        """Write PDB file with missing occupancy."""
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

    def test_mmcifio_write_structure(self):
        """Write a full structure using MMCIFIO."""
        io = MMCIFIO()
        struct1 = self.structure
        # Write full model to temp file
        io.set_structure(struct1)
        filenumber, filename = tempfile.mkstemp()
        os.close(filenumber)
        try:
            io.save(filename)
            struct2 = self.mmcif_parser.get_structure("1a8o", filename)
            nresidues = len(list(struct2.get_residues()))
            self.assertEqual(len(struct2), 1)
            self.assertEqual(nresidues, 158)
        finally:
            os.remove(filename)

    def test_mmcifio_write_residue(self):
        """Write a single residue using MMCIFIO."""
        io = MMCIFIO()
        struct1 = self.structure
        residue1 = list(struct1.get_residues())[0]
        # Write full model to temp file
        io.set_structure(residue1)
        filenumber, filename = tempfile.mkstemp()
        os.close(filenumber)
        try:
            io.save(filename)
            struct2 = self.mmcif_parser.get_structure("1a8o", filename)
            nresidues = len(list(struct2.get_residues()))
            self.assertEqual(nresidues, 1)
        finally:
            os.remove(filename)

    def test_mmcifio_select(self):
        """Write a selection of the structure using a Select subclass."""
        # Selection class to filter all alpha carbons
        class CAonly(Select):
            """Accepts only CA residues."""

            def accept_atom(self, atom):
                if atom.name == "CA" and atom.element == "C":
                    return 1

        io = MMCIFIO()
        struct1 = self.structure
        # Write to temp file
        io.set_structure(struct1)
        filenumber, filename = tempfile.mkstemp()
        os.close(filenumber)
        try:
            io.save(filename, CAonly())
            struct2 = self.mmcif_parser.get_structure("1a8o", filename)
            nresidues = len(list(struct2.get_residues()))
            self.assertEqual(nresidues, 70)
        finally:
            os.remove(filename)

    def test_mmcifio_write_dict(self):
        """Write an mmCIF dictionary out, read it in and compare them."""
        d1 = MMCIF2Dict(self.mmcif_file)
        io = MMCIFIO()
        # Write to temp file
        io.set_dict(d1)
        filenumber, filename = tempfile.mkstemp()
        os.close(filenumber)
        try:
            io.save(filename)
            d2 = MMCIF2Dict(filename)
            k1 = sorted(d1.keys())
            k2 = sorted(d2.keys())
            self.assertEqual(k1, k2)
            for key in k1:
                self.assertEqual(d1[key], d2[key])
        finally:
            os.remove(filename)

    def test_mmcifio_multimodel(self):
        """Write a multi-model, multi-chain mmCIF file."""
        pdb_struct = self.parser.get_structure("1SSU_mod_pdb", self.mmcif_multimodel_pdb_file)
        mmcif_struct = self.mmcif_parser.get_structure("1SSU_mod_mmcif", self.mmcif_multimodel_mmcif_file)
        io = MMCIFIO()
        for struct in [pdb_struct, mmcif_struct]:
            io.set_structure(struct)
            filenumber, filename = tempfile.mkstemp()
            os.close(filenumber)
            try:
                io.save(filename)
                struct_in = self.mmcif_parser.get_structure("1SSU_mod_in", filename)
                self.assertEqual(len(struct_in), 2)
                self.assertEqual(len(struct_in[1]), 2)
                self.assertEqual(round(float(struct_in[1]["B"][1]["N"].get_coord()[0]), 3), 6.259)
            finally:
                os.remove(filename)


class Exposure(unittest.TestCase):
    """Testing Bio.PDB.HSExposure."""

    def setUp(self):
        pdb_filename = "PDB/a_structure.pdb"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            structure = PDBParser(PERMISSIVE=True).get_structure('X', pdb_filename)
        structure[1].detach_child('B')
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
    """induces Atom Element from Atom Name."""

    def setUp(self):
        pdb_filename = "PDB/a_structure.pdb"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            structure = PDBParser(PERMISSIVE=True).get_structure('X', pdb_filename)
        self.residue = structure[0]['A'][('H_PCA', 1, ' ')]

    def test_AtomElement(self):
        """Atom Element."""
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
            O=(' OH ',),  # noqa: E741
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
        self.assertEqual(chains, ['A', 'A', 'B', 'C', ' '])

    def test_get_residues(self):
        """Yields all residues from all models."""
        residues = [resi.id for resi in self.struc.get_residues()]
        self.assertEqual(len(residues), 179)

    def test_get_atoms(self):
        """Yields all atoms from the structure, excluding duplicates and ALTLOCs which are not parsed."""
        atoms = ["%12s" % str((atom.id, atom.altloc)) for atom in self.struc.get_atoms()]
        self.assertEqual(len(atoms), 835)


class ChangingIdTests(unittest.TestCase):

    def setUp(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            self.struc = PDBParser(PERMISSIVE=True).get_structure(
                                                  'X', "PDB/a_structure.pdb")

    def test_change_model_id(self):
        """Change the id of a model."""
        for model in self.struc:
            break  # Get first model in structure
        model.id = 2
        self.assertEqual(model.id, 2)
        self.assertIn(2, self.struc)
        self.assertNotIn(0, self.struc)

    def test_change_model_id_raises(self):
        """Cannot change id to a value already in use by another child."""
        model = next(iter(self.struc))
        with self.assertRaises(ValueError):
            model.id = 1
        # Make sure nothing was changed
        self.assertEqual(model.id, 0)
        self.assertIn(0, self.struc)
        self.assertIn(1, self.struc)

    def test_change_chain_id(self):
        """Change the id of a model."""
        chain = next(iter(self.struc.get_chains()))
        chain.id = "R"
        self.assertEqual(chain.id, "R")
        model = next(iter(self.struc))
        self.assertIn("R", model)

    def test_change_id_to_self(self):
        """Changing the id to itself does nothing (does not raise)."""
        chain = next(iter(self.struc.get_chains()))
        chain_id = chain.id
        chain.id = chain_id
        self.assertEqual(chain.id, chain_id)

    def test_change_residue_id(self):
        """Change the id of a residue."""
        chain = next(iter(self.struc.get_chains()))
        res = chain[('H_PCA', 1, ' ')]
        res.id = (' ', 1, ' ')

        self.assertEqual(res.id, (' ', 1, ' '))
        self.assertIn((' ', 1, ' '), chain)
        self.assertNotIn(('H_PCA', 1, ' '), chain)
        self.assertEqual(chain[(' ', 1, ' ')], res)

    def test_full_id_is_updated_residue(self):
        """Invalidate cached full_ids if an id is changed."""
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
        """Invalidate cached full_ids if an id is changed."""
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
        """Sum of positions of atoms in an entity along with the number of atoms."""
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
        """Average atom position in an entity."""
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
        """Test Vector object."""
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


class PDBParserTests(unittest.TestCase):
    """Test PDBParser module."""

    def test_PDBParser(self):
        """Walk down the structure hierarchy and test parser reliability."""
        p = PDBParser(PERMISSIVE=True)
        filename = "PDB/1A8O.pdb"
        s = p.get_structure("scr", filename)
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

    def test_entity_copy(self):
        """Make a copy of a residue."""
        for e in (self.s, self.m, self.c, self.r):
            ee = e.copy()
            self.assertFalse(e is ee)
            self.assertFalse(e.get_list()[0] is ee.get_list()[0])


def eprint(*args, **kwargs):
    """Helper function that prints to stderr."""
    print(*args, file=sys.stderr, **kwargs)


def will_it_float(s):
    """Helper function that converts the input into a float if it is a number.

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
        """Test parsing of pregenerated DSSP."""
        dssp, keys = make_dssp_dict("PDB/2BEG.dssp")
        self.assertEqual(len(dssp), 130)

    def test_DSSP_noheader_file(self):
        """Test parsing of pregenerated DSSP missing header information."""
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
        """All elements correctly added to xtra attribute of input model object."""
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
