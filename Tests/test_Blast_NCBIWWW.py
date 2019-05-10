# Copyright 2019 by Markus Piotrowski.
# All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Tests for the Bio.Blast.NCBIWWW module.

Be careful with adding more tests here, since these searches can take
a very long time.
"""

import unittest
import warnings
import time

from Bio import BiopythonWarning

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

import requires_internet

requires_internet.check()


class TestNCBIWWW(unittest.TestCase):
    """Test Blast online search via (formerly) QBLAST server."""

    def test_short_query(self):
        """Test SHORT_QUERY_ADJUST parameter."""
        # Should give no hits:
        my_search = NCBIWWW.qblast('blastp', 'nr', 'ICWENRM', hitlist_size=5)
        my_hits = NCBIXML.read(my_search)
        my_search.close()
        self.assertEqual(len(my_hits.alignments), 0)

        time.sleep(15)  # May prevent penalizing by NCBI?

        # Should give hits:
        my_search = NCBIWWW.qblast('blastp', 'nr', 'ICWENRM', hitlist_size=5,
                                   short_query=True)
        my_hits = NCBIXML.read(my_search)
        my_search.close()
        self.assertEqual(len(my_hits.alignments), 5)

        time.sleep(15)

        with warnings.catch_warnings(record=True) as w:
            # Cause all warnings to always be triggered.
            warnings.simplefilter('always')
            # Trigger a warning.
            my_search = NCBIWWW.qblast('blastn', 'nt', 'ATGTCAACTTCAGAA',
                                       hitlist_size=5, short_query=True)
            # Verify some things
            self.assertEqual(len(w), 1)
            self.assertEqual(w[-1].category, BiopythonWarning)
            self.assertIn('blastn', str(w[-1].message))
            my_hits = NCBIXML.read(my_search)
            my_search.close()
            self.assertEqual(len(my_hits.alignments), 5)

    def test_error_conditions(self):
        """Test if exceptions were properly handled."""
        self.assertRaises(ValueError, NCBIWWW.qblast, 'megablast', 'nt',
                          'ATGCGTACGCAGCTAAAGTAAACCTATCGCGTCTCCT')


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
