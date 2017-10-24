import os
import sys
import unittest

from Bio import MissingExternalDependencyError
from Bio.PDB import PDBParser, PSEA


os.environ['LANG'] = 'C'

psea = None
from Bio._py3k import getoutput
try:
    output = getoutput("p-sea -h")
    if output.startswith("displays this help"):
        psea = "p-sea"
except OSError:
    pass

if not psea:
    raise MissingExternalDependencyError(
        "Install p-sea if you want to use Bio.PSEA from Biopython.")


class PSEATests(unittest.TestCase):
    """Tests for PSEA module"""
    def test_psea(self):
        """Self-test for PSEA"""
        pdb1 = "PDB/1A8O.pdb"
        p = PDBParser()
        s = p.get_structure('X', pdb1)
        self.assertEqual(PSEA.PSEA(s[0], pdb1))


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
