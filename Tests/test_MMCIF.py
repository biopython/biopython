# Copyright 2012 Lenna X. Peterson (arklenna@gmail.com).
# All rights reserved.
#
# Tests adapted from test_PDB.py
#
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for the MMCIF portion of the Bio.PDB module."""

import unittest

try:
    import numpy
    from numpy import dot  # Missing on old PyPy's micronumpy
    del dot
    from numpy.linalg import svd, det # Missing in PyPy 2.0 numpypy
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.PDB.")


from Bio.Seq import Seq
from Bio.Alphabet import generic_protein
from Bio.PDB.PDBExceptions import PDBConstructionException, PDBConstructionWarning

from Bio.PDB import PPBuilder, CaPPBuilder
from Bio.PDB.MMCIFParser import MMCIFParser


class ParseReal(unittest.TestCase):
    """Testing with real CIF file(s)."""

    def test_parser(self):
        """Extract polypeptides from 1A80."""
        parser = MMCIFParser()
        structure = parser.get_structure("example", "PDB/1A8O.cif")
        self.assertEqual(len(structure), 1)
        for ppbuild in [PPBuilder(), CaPPBuilder()]:
            #==========================================================
            # Check that serial_num (model column) is stored properly
            self.assertEqual(structure[0].serial_num, 1)
            #First try allowing non-standard amino acids,
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
            #Here non-standard MSE are shown as M
            self.assertEqual("MDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQ"
                             "NANPDCKTILKALGPGATLEEMMTACQG", str(s))
            #==========================================================
            #Now try strict version with only standard amino acids
            #Should ignore MSE 151 at start, and then break the chain
            #at MSE 185, and MSE 214,215
            polypeptides = ppbuild.build_peptides(structure[0], True)
            self.assertEqual(len(polypeptides), 3)
            #First fragment
            pp = polypeptides[0]
            self.assertEqual(pp[0].get_id()[1], 152)
            self.assertEqual(pp[-1].get_id()[1], 184)
            s = pp.get_sequence()
            self.assertTrue(isinstance(s, Seq))
            self.assertEqual(s.alphabet, generic_protein)
            self.assertEqual("DIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNW", str(s))
            #Second fragment
            pp = polypeptides[1]
            self.assertEqual(pp[0].get_id()[1], 186)
            self.assertEqual(pp[-1].get_id()[1], 213)
            s = pp.get_sequence()
            self.assertTrue(isinstance(s, Seq))
            self.assertEqual(s.alphabet, generic_protein)
            self.assertEqual("TETLLVQNANPDCKTILKALGPGATLEE", str(s))
            #Third fragment
            pp = polypeptides[2]
            self.assertEqual(pp[0].get_id()[1], 216)
            self.assertEqual(pp[-1].get_id()[1], 220)
            s = pp.get_sequence()
            self.assertTrue(isinstance(s, Seq))
            self.assertEqual(s.alphabet, generic_protein)
            self.assertEqual("TACQG", str(s))

    def testModels(self):
        """Test file with multiple models"""
        parser = MMCIFParser()
        structure = parser.get_structure("example", "PDB/1LCD.cif")
        self.assertEqual(len(structure), 3)
        for ppbuild in [PPBuilder(), CaPPBuilder()]:
                #==========================================================
                # Check that serial_num (model column) is stored properly
                self.assertEqual(structure[0].serial_num, 1)
                self.assertEqual(structure[1].serial_num, 2)
                self.assertEqual(structure[2].serial_num, 3)
                #First try allowing non-standard amino acids,
                polypeptides = ppbuild.build_peptides(structure[0], False)
                self.assertEqual(len(polypeptides), 1)
                pp = polypeptides[0]
                # Check the start and end positions
                self.assertEqual(pp[0].get_id()[1], 1)
                self.assertEqual(pp[-1].get_id()[1], 51)
                # Check the sequence
                s = pp.get_sequence()
                self.assertTrue(isinstance(s, Seq))
                self.assertEqual(s.alphabet, generic_protein)
                #Here non-standard MSE are shown as M
                self.assertEqual("MKPVTLYDVAEYAGVSYQTVSRVVNQASHVSAKTREKVEAAMAELNYIPNR",
                                 str(s))
                #==========================================================
                #Now try strict version with only standard amino acids
                polypeptides = ppbuild.build_peptides(structure[0], True)
                self.assertEqual(len(polypeptides), 1)
                pp = polypeptides[0]
                # Check the start and end positions
                self.assertEqual(pp[0].get_id()[1], 1)
                self.assertEqual(pp[-1].get_id()[1], 51)
                # Check the sequence
                s = pp.get_sequence()
                self.assertTrue(isinstance(s, Seq))
                self.assertEqual(s.alphabet, generic_protein)
                self.assertEqual("MKPVTLYDVAEYAGVSYQTVSRVVNQASHVSAKTREKVEAAMAELNYIPNR",
                                 str(s))

if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
