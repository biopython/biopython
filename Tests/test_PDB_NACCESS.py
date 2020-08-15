# Copyright 2009-2011 by Eric Talevich.  All rights reserved.
# Revisions copyright 2009-2013 by Peter Cock.  All rights reserved.
# Revisions copyright 2013 Lenna X. Peterson. All rights reserved.
# Revisions copyright 2013 Gokcen Eraslan. All rights reserved.
# Revisions copyright 2020 Joao Rodrigues. All rights reserved.
#
# Converted by Eric Talevich from an older unit test copyright 2002
# by Thomas Hamelryck.
#
# Merged related test files into one, by Joao Rodrigues (2020)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Unit tests for the Bio.PDB.NACCESS submodule."""

import subprocess
import unittest
import warnings

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.PDB."
    ) from None


from Bio import MissingExternalDependencyError
from Bio.PDB import PDBParser
from Bio.PDB.NACCESS import NACCESS, process_asa_data, process_rsa_data


class NACCESS_test(unittest.TestCase):
    """Tests for Bio.PDB.NACCESS and output parsing."""

    def test_NACCESS_rsa_file(self):
        """Test parsing of pregenerated rsa NACCESS file."""
        with open("PDB/1A8O.rsa") as rsa:
            naccess = process_rsa_data(rsa)
        self.assertEqual(len(naccess), 66)

    def test_NACCESS_asa_file(self):
        """Test parsing of pregenerated asa NACCESS file."""
        with open("PDB/1A8O.asa") as asa:
            naccess = process_asa_data(asa)
        self.assertEqual(len(naccess), 524)

    def test_NACCESS(self):
        """Test calling NACCESS from Bio.PDB."""
        # Check if NACCESS is available
        try:
            subprocess.check_call(
                ["naccess", "-q"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT
            )
        except OSError:
            raise self.skipTest("Install naccess if you want to use it from Biopython.")

        p = PDBParser()
        pdbfile = "PDB/1A8O.pdb"
        model = p.get_structure("1A8O", pdbfile)[0]
        naccess = NACCESS(model, pdbfile)
        self.assertEqual(len(naccess), 66)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
