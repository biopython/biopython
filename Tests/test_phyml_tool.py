# Copyright (C) 2012 by Eric Talevich.
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for Bio.Phylo.Applications wrappers."""

import sys
import os
import unittest

from Bio import Phylo
from Bio.Phylo.Applications import PhymlCommandline
from Bio import MissingExternalDependencyError

#Try to avoid problems when the OS is in another language
os.environ['LANG'] = 'C'

phyml_exe = None
if sys.platform=="win32":
    raise MissingExternalDependencyError(
        "Testing PhyML on Windows not supported yet")
else:
    from Bio._py3k import getoutput
    output = getoutput("phyml --version")
    if "not found" not in output and "20" in output:
        phyml_exe = "phyml"

if not phyml_exe:
    raise MissingExternalDependencyError(
        "Install PhyML 3.0 if you want to use the Bio.Phylo.Applications wrapper.")


# Example Phylip file with 4 aligned protein sequences
EX_PHYLIP = 'Phylip/interlaced2.phy'


class AppTests(unittest.TestCase):
    """Tests for application wrappers."""

    def test_phyml(self):
        """Run PhyML using the wrapper."""
        cmd = PhymlCommandline(phyml_exe, input=EX_PHYLIP, datatype='aa')
        # Smoke test
        try:
            out, err = cmd()
            self.assertTrue(len(out) > 0)
            self.assertEqual(len(err), 0)
            # Check the output tree
            tree = Phylo.read(EX_PHYLIP + '_phyml_tree.txt', 'newick')
            self.assertEqual(tree.count_terminals(), 4)
        finally:
            # Clean up generated files
            for suffix in ['_phyml_tree.txt', '_phyml_stats.txt']:
                fname = EX_PHYLIP + suffix
                if os.path.isfile(fname):
                    os.remove(fname)


# ---------------------------------------------------------

if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
