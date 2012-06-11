# Copyright (C) 2012 by Eric Talevich.
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for Bio.Phylo.Applications wrappers."""

import os
import unittest
# from cStringIO import StringIO

from Bio import Phylo
from Bio.Phylo.Applications import RaxmlCommandline

# Example Phylip file with 4 aligned protein sequences
EX_PHYLIP = 'Phylip/interlaced2.phy'

class AppTests(unittest.TestCase):
    """Tests for application wrappers."""

    def test_raxml(self):
        """Run RAxML using the wrapper."""
        cmd = RaxmlCommandline(sequences=EX_PHYLIP, model="PROTCATWAG",
                               name="test")
        # The parsimony seed should be set automatically
        self.assert_('-p' in str(cmd))
        # Smoke test
        try:
            out, err = cmd()
            self.assert_(len(out) > 0)
            self.assert_(len(err) == 0)
            # Check the output tree
            tree = Phylo.read('RAxML_bestTree.test', 'newick')
            self.assertEqual(tree.count_terminals(), 4)
        finally:
            # Remove RAxML-generated files, or RAxML will complain bitterly
            # during the next run
            for fname in ['RAxML_info.test',
                          'RAxML_log.test',
                          'RAxML_parsimonyTree.test',
                          'RAxML_bestTree.test',
                          'RAxML_result.test']:
                if os.path.isfile(fname):
                    os.remove(fname)


# ---------------------------------------------------------

if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
