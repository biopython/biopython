# Copyright (C) 2012 by Eric Talevich.
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for Bio.Phylo.Applications wrappers."""

import sys
import os
import unittest

from Bio import Phylo
from Bio.Phylo.Applications import RaxmlCommandline
from Bio import MissingExternalDependencyError

raxml_exe = None
if sys.platform=="win32":
    raise MissingExternalDependencyError(
          "Testing RAxML on Windows not supported yet")
else:
    from Bio._py3k import getoutput
    output = getoutput("raxmlHPC -v")
    if "not found" not in output and "This is RAxML" in output:
        raxml_exe = "raxmlHPC"
if not raxml_exe:
    raise MissingExternalDependencyError(
        "Install RAxML (binary raxmlHPC) if you want to test the Bio.Phylo.Applications wrapper.")

# Example Phylip file with 4 aligned protein sequences
EX_PHYLIP = 'Phylip/interlaced2.phy'


class AppTests(unittest.TestCase):
    """Tests for application wrappers."""

    def test_raxml(self):
        """Run RAxML using the wrapper."""
        cmd = RaxmlCommandline(raxml_exe,
                               sequences=EX_PHYLIP, model="PROTCATWAG",
                               name="test")
        # The parsimony seed should be set automatically
        self.assertTrue('-p' in str(cmd))
        # Smoke test
        try:
            out, err = cmd()
            self.assertTrue(len(out) > 0)
            self.assertEqual(len(err), 0)
            # Check the output tree
            tree = Phylo.read('RAxML_result.test', 'newick')
            self.assertEqual(tree.count_terminals(), 4)
        finally:
            # Remove RAxML-generated files, or RAxML will complain bitterly
            # during the next run
            for fname in ['RAxML_info.test',
                          'RAxML_log.test',
                          'RAxML_parsimonyTree.test',
                          'RAxML_result.test',
                          # Present in 7.2.X+  but not 7.0.4:
                          'RAxML_bestTree.test',
                         ]:
                if os.path.isfile(fname):
                    os.remove(fname)


# ---------------------------------------------------------

if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
