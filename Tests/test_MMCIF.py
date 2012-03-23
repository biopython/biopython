# Written 2012 Lenna X. Peterson
# arklenna@gmail.com

"""Unit tests for the MMCIF portion of the Bio.PDB module."""

import os
import tempfile
import unittest
import warnings

try:
    import numpy
    from numpy import dot #Missing on PyPy's micronumpy
    del dot
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.PDB.")

#from Bio.Seq import Seq
#from Bio.Alphabet import generic_protein
#from Bio.PDB import PDBParser, PPBuilder, CaPPBuilder, PDBIO
#from Bio.PDB import HSExposureCA, HSExposureCB, ExposureCN
from Bio.PDB.PDBExceptions import PDBConstructionException, PDBConstructionWarning
#from Bio.PDB import rotmat, Vector

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB import PPBuilder

class ParseReal(unittest.TestCase):
    """Testing with real CIF file(s)."""

    def test_parser(self):
        parser = MMCIFParser()
        structure = parser.get_structure("example", "PDB/1A8O.cif")
        self.assertEqual(len(structure), 1)
        #polypeptides = 

if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
