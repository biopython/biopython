# Copyright 2016 Joao P. Rodrigues (j.p.g.l.m.rodrigues@gmail.com).
# All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests to compare objects parsed by different structure parsers."""

import os
import unittest
import warnings

try:
    import numpy
    from numpy import dot  # Missing on old PyPy's micronumpy

    del dot
    from numpy.linalg import svd, det  # Missing in PyPy 2.0 numpypy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.PDB."
    ) from None

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionException, PDBConstructionWarning


class CompareStructures(unittest.TestCase):
    """Tests for comparing the same structure parsed by PDB and MMCIF parsers."""

    def setUp(self):

        # Silence!
        warnings.simplefilter("ignore", PDBConstructionWarning)

        pdbparser = PDBParser(QUIET=1)
        cifparser = MMCIFParser(QUIET=1)

        modpath = os.path.abspath(os.path.dirname(__file__))

        pdb_file = os.path.join(modpath, "PDB", "1LCD.pdb")
        cif_file = os.path.join(modpath, "PDB", "1LCD.cif")

        self.pdbo = pdbparser.get_structure("pdb", pdb_file)
        self.cifo = cifparser.get_structure("pdb", cif_file)

    def test_compare_models(self):
        """Compared parsed models."""
        cif_models = [(m.id, len(m.child_list)) for m in self.cifo.get_models()]
        pdb_models = [(m.id, len(m.child_list)) for m in self.pdbo.get_models()]

        self.assertEqual(len(cif_models), len(pdb_models))
        self.assertEqual([i[0] for i in cif_models], [i[0] for i in pdb_models])
        self.assertEqual([i[1] for i in cif_models], [i[1] for i in pdb_models])

    def test_compare_chains(self):
        """Compare parsed chains."""
        cif_chains = [(c.id, len(c.child_list)) for c in self.cifo.get_chains()]
        pdb_chains = [(c.id, len(c.child_list)) for c in self.pdbo.get_chains()]

        self.assertEqual(len(cif_chains), len(pdb_chains))
        self.assertEqual([i[0] for i in cif_chains], [i[0] for i in pdb_chains])
        self.assertEqual([i[1] for i in cif_chains], [i[1] for i in pdb_chains])


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
