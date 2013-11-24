# Copyright 2013 Gokcen Eraslan. All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit test for Bio.PDB.NACCESS which needs NACCESS tool.

See also test_PDB.py for dependency free NACCESS tests.
"""

import subprocess
import unittest


from Bio import MissingExternalDependencyError
from Bio.PDB import PDBParser
from Bio.PDB.NACCESS import NACCESS

# Check if NACCESS is installed
try:
    subprocess.check_call(["naccess", "-q"],
                          stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT)
except OSError:
    raise MissingExternalDependencyError(
        "Install naccess if you want to use it from Biopython.")


class NACCESS_test(unittest.TestCase):
    """Test NACCESS module"""

    def test_NACCESS(self):
        """Test NACCESS generation from PDB"""
        p = PDBParser()
        pdbfile = "PDB/1A8O.pdb"
        model = p.get_structure("1A8O", pdbfile)[0]
        naccess = NACCESS(model, pdbfile)
        self.assertEqual(len(naccess), 66)


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
