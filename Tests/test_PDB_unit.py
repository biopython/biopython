# Copyright 2009 by Eric Talevich.  All rights reserved.
# Revisions copyright 2009 by Peter Cock.  All rights reserved.
#
# Converted by Eric Talevich from an older unit test copyright 2002
# by Thomas Hamelryck.
# 
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for the Bio.PDB module."""
import unittest
import warnings

try:
    from numpy.random import random
except ImportError:
    from Bio import MissingExternalDependencyError
    raise MissingExternalDependencyError(\
        "Install NumPy if you want to use Bio.PDB.")
 
from Bio.PDB import PDBParser, PPBuilder, CaPPBuilder
from Bio.PDB import HSExposureCA, HSExposureCB, ExposureCN
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBExceptions import PDBConstructionException, PDBConstructionWarning

class PDBNeighborTest(unittest.TestCase):
    def setUp(self):
        warnings.resetwarnings()

    def test_neighbor_search(self):
        """NeighborSearch: Find nearby randomly generated coordinates.
         
        Based on the self test in Bio.PDB.NeighborSearch.
        """
        class RandomAtom:
            def __init__(self):
                self.coord = 100 * random(3)
            def get_coord(self):
                return self.coord
        for i in range(0, 20):
            atoms = [RandomAtom() for j in range(100)]
            ns = NeighborSearch(atoms)
            hits = ns.search_all(5.0)
            self.assert_(hits >= 0)
 
 
class PDBExceptionTest(unittest.TestCase):
    def test_strict(self):
        """Check error: Parse a flawed PDB file in strict mode."""
        warnings.resetwarnings()
        parser = PDBParser(PERMISSIVE=False)
        self.assertRaises(PDBConstructionException,
                parser.get_structure, "example", "PDB/a_structure.pdb")

    #TODO - check get expected warnings, may require Python 2.6+
    #See Bug 2820


class PDBParseTest(unittest.TestCase):
    def setUp(self):
        warnings.resetwarnings()
        warnings.simplefilter('ignore', PDBConstructionWarning)
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
 
    def test_ca_ca(self):
        """Extract polypeptides using CA-CA."""
        ppbuild = CaPPBuilder()
        polypeptides = ppbuild.build_peptides(self.structure[1])
        self.assertEqual(len(polypeptides), 1)
        pp = polypeptides[0]
        # Check the start and end positions
        self.assertEqual(pp[0].get_id()[1], 2)
        self.assertEqual(pp[-1].get_id()[1], 86)
 
    def test_structure(self):
        """Verify the structure of the parsed example PDB file."""
        # Structure contains 2 models
        self.assertEquals(len(self.structure), 2)
        # --- Checking model 0 ---
        m0 = self.structure[0]
        # Model 0 contains 1 chain
        self.assertEquals(len(m0), 1)
        # Chain 'A' contains 1 residue
        self.assertEquals(len(m0['A']), 1)
        # Residue ('H_PCA', 1, ' ') contains 8 atoms.
        residue = m0['A'].get_list()[0]
        self.assertEquals(residue.get_id(), ('H_PCA', 1, ' '))
        self.assertEquals(len(residue), 8)
        # --- Checking model 1 ---
        m1 = self.structure[1]
        # Model 1 contains 3 chains
        self.assertEquals(len(m1), 3)
        # Deconstruct this data structure to check each chain
        chain_data = [ # chain_id, chain_len, [(residue_id, residue_len), ...]
            ('A', 86, [ ((' ', 0, ' '), 1 ),
                        ((' ', 2, ' '), 11),
                        ((' ', 3, ' '), 6, 1), # disordered
                        ((' ', 4, ' '), 4 ),
                        ((' ', 5, ' '), 6 ),
                        ((' ', 6, ' '), 9 ),
                        ((' ', 7, ' '), 4 ),
                        ((' ', 8, ' '), 4 ),
                        ((' ', 9, ' '), 4 ),
                        ((' ', 10, ' '), 6, ['GLY', 'SER']), # point mut
                        ((' ', 11, ' '), 7 ),
                        ((' ', 12, ' '), 6 ),
                        ((' ', 13, ' '), 7 ),
                        ((' ', 14, ' '), 4, ['ALA', 'GLY']), # point mut
                        ((' ', 15, ' '), 8, 3), # disordered
                        ((' ', 16, ' '), 11, ['ARG', 'TRP']), # point mut
                        ((' ', 17, ' '), 6 ),
                        ((' ', 18, ' '), 6 ),
                        ((' ', 19, ' '), 6 ),
                        ((' ', 20, ' '), 8 ),
                        ((' ', 21, ' '), 14),
                        ((' ', 22, ' '), 4 ),
                        ((' ', 23, ' '), 14),
                        ((' ', 24, ' '), 6 ),
                        ((' ', 25, ' '), 4 ),
                        ((' ', 26, ' '), 8 ),
                        ((' ', 27, ' '), 6 ),
                        ((' ', 28, ' '), 9, 5), # disordered
                        ((' ', 29, ' '), 7 ),
                        ((' ', 30, ' '), 12),
                        ((' ', 31, ' '), 6 ),
                        ((' ', 32, ' '), 4 ),
                        ((' ', 33, ' '), 11),
                        ((' ', 34, ' '), 7 ),
                        ((' ', 35, ' '), 6 ),
                        ((' ', 36, ' '), 9 ),
                        ((' ', 37, ' '), 8 ),
                        ((' ', 38, ' '), 9 ),
                        ((' ', 39, ' '), 6 ),
                        ((' ', 40, ' '), 14),
                        ((' ', 41, ' '), 6 ),
                        ((' ', 42, ' '), 4 ),
                        ((' ', 43, ' '), 9 ),
                        ((' ', 44, ' '), 11),
                        ((' ', 45, ' '), 6, 1), # disordered
                        ((' ', 46, ' '), 8 ),
                        ((' ', 47, ' '), 10),
                        ((' ', 48, ' '), 11),
                        ((' ', 49, ' '), 6 ),
                        ((' ', 50, ' '), 4 ),
                        ((' ', 51, ' '), 5 ),
                        ((' ', 52, ' '), 5 ),
                        ((' ', 53, ' '), 7 ),
                        ((' ', 54, ' '), 4 ),
                        ((' ', 55, ' '), 8 ),
                        ((' ', 56, ' '), 7 ),
                        ((' ', 57, ' '), 7 ),
                        ((' ', 58, ' '), 6 ),
                        ((' ', 59, ' '), 4 ),
                        ((' ', 60, ' '), 9 ),
                        ((' ', 61, ' '), 8 ),
                        ((' ', 62, ' '), 11),
                        ((' ', 63, ' '), 6 ),
                        ((' ', 64, ' '), 6 ),
                        ((' ', 65, ' '), 6 ),
                        ((' ', 66, ' '), 7 ),
                        ((' ', 67, ' '), 10),
                        ((' ', 68, ' '), 4 ),
                        ((' ', 69, ' '), 14),
                        ((' ', 70, ' '), 6 ),
                        ((' ', 71, ' '), 4 ),
                        ((' ', 72, ' '), 4 ),
                        ((' ', 73, ' '), 4 ),
                        ((' ', 74, ' '), 8, 3), # disordered
                        ((' ', 75, ' '), 8 ),
                        ((' ', 76, ' '), 12),
                        ((' ', 77, ' '), 6 ),
                        ((' ', 78, ' '), 6 ),
                        ((' ', 79, ' '), 4, 4), # disordered
                        ((' ', 80, ' '), 4, ['GLY', 'SER']), # point mut
                        ((' ', 81, ' '), 8, ['ASN', 'LYS']), # point mut
                        ((' ', 82, ' '), 6 ),
                        ((' ', 83, ' '), 9 ),
                        ((' ', 84, ' '), 12),
                        ((' ', 85, ' '), 11),
                        ((' ', 86, ' '), 6 ),
                        ]),
            ('B', 4, [ (('H_NAG', 1, ' '), 14),
                        (('H_NAG', 2, ' '), 14),
                        (('H_NAG', 3, ' '), 14),
                        (('H_NAG', 4, ' '), 14),
                        ]),
            (' ', 76, [ (('W', 1, ' '), 1),
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
            self.assertEquals(chain.get_id(), chn[0])
            self.assertEquals(len(chain), chn[1])
            for r_idx, res in enumerate(chn[2]):
                residue = chain.get_list()[r_idx]
                # Check residue ID and atom count
                self.assertEquals(residue.get_id(), res[0])
                self.assertEquals(len(residue), res[1])
                disorder_lvl = residue.is_disordered()
                if disorder_lvl == 1:
                    # Check the number of disordered atoms
                    disordered_count = sum(1 for atom in residue
                                           if atom.is_disordered())
                    if disordered_count:
                        self.assertEquals(disordered_count, res[2])
                elif disorder_lvl == 2:
                    # Point mutation -- check residue names
                    self.assertEquals(residue.disordered_get_id_list(), res[2])

    def test_details(self):
        """Verify details of the parsed example PDB file."""
        structure = self.structure
        self.assertEqual(len(structure), 2)

        #First model
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
                         "N CA CB CG CD OE C O")
        self.assertEqual(" ".join(atom.element for atom in chain.get_atoms()),
                         "N C C C C O C O")
        #Second model
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

class Exposure(unittest.TestCase):
    "Testing Bio.PDB.HSExposure."
    def setUp(self):
        warnings.resetwarnings()
        warnings.simplefilter('ignore', PDBConstructionWarning)
        pdb_filename = "PDB/a_structure.pdb"
        structure=PDBParser(PERMISSIVE=True).get_structure('X', pdb_filename)
        self.model=structure[1]
        #Look at first chain only
        a_residues=list(self.model["A"].child_list)
        self.assertEqual(86, len(a_residues))
        self.assertEqual(a_residues[0].get_resname(), "CYS")
        self.assertEqual(a_residues[1].get_resname(), "ARG")
        self.assertEqual(a_residues[2].get_resname(), "CYS")
        self.assertEqual(a_residues[3].get_resname(), "GLY")
        #...
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
        #...
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
        #...
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
        #...
        self.assertEqual(1, len(residues[-2].xtra))
        self.assertEqual(48, residues[-2].xtra["EXP_CN"])
        self.assertEqual(1, len(residues[-1].xtra))
        self.assertEqual(38, residues[-1].xtra["EXP_CN"])

class AssortedMisc(unittest.TestCase):
    "Testing with real PDB files."

    def test_strict(self):
        """Parse 1A8O.pdb file in strict mode."""
        warnings.resetwarnings()
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


# -------------------------------------------------------------

if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
