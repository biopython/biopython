# Copyright 2013 Lenna X. Peterson. All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit test for Bio.PDB.DSSP which need the binary tool.

See also test_PDB.py for dependency free DSSP tests.
"""

from distutils.version import StrictVersion
import re
import subprocess
import unittest


from Bio import MissingExternalDependencyError
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB import DSSP

# DSSP version, if known, as DSSP <2.2.0 does not support mmcif files
dssp_version = '0.0.0'
# Check if DSSP is installed
quiet_kwargs = dict(stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
try:
    try:
        # Newer versions of DSSP
        version_string = subprocess.check_output(["dssp", "--version"],
                                                 universal_newlines=True)
        dssp_version = re.search(r'\s*([\d.]+)', version_string).group(1)
    except subprocess.CalledProcessError:
        # Older versions of DSSP
        subprocess.check_call(["dssp", "-h"], **quiet_kwargs)
except OSError:
    try:
        version_string = subprocess.check_output(["mkdssp", "--version"],
                                                 universal_newlines=True)
        dssp_version = re.search(r'\s*([\d.]+)', version_string).group(1)
    except OSError:
        raise MissingExternalDependencyError(
            "Install dssp if you want to use it from Biopython.")


class DSSP_test(unittest.TestCase):
    """Test DSSP module"""

    def test_dssp(self):
        """Test DSSP generation from PDB"""
        p = PDBParser()
        pdbfile = "PDB/2BEG.pdb"
        model = p.get_structure("2BEG", pdbfile)[0]
        dssp = DSSP(model, pdbfile)
        self.assertEqual(len(dssp), 130)

    # Only run mmCIF tests if DSSP version installed supports mmcif
    if StrictVersion(dssp_version) >= StrictVersion('2.2.0'):
        def test_dssp_with_mmcif_file(self):
            """Test DSSP generation from MMCIF"""
            p = MMCIFParser()
            pdbfile = "PDB/2BEG.cif"
            model = p.get_structure("2BEG", pdbfile)[0]
            dssp = DSSP(model, pdbfile)
            self.assertEqual(len(dssp), 130)

        def test_dssp_with_mmcif_file_and_nonstandard_residues(self):
            """Test DSSP generation from MMCIF with non-standard residues"""
            p = MMCIFParser()
            pdbfile = "PDB/1AS5.cif"
            model = p.get_structure("1AS5", pdbfile)[0]
            dssp = DSSP(model, pdbfile)
            self.assertEqual(len(dssp), 24)


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
