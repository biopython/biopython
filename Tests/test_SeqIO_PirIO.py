# Copyright 2026 by Contributors.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.SeqIO.PirIO module."""

import unittest
from io import StringIO

from Bio.SeqIO.PirIO import PirHeaderParser


class TestPirHeaderParser(unittest.TestCase):
    """Test PirHeaderParser directly."""

    ins_pir = [
        ">P1;HLA:HLA00489\ndescription\nsequence\n*",
        ">P1;HLA:HLA00490\ndesc\nseq\n*",
    ]
    outs_pir = [["P1;HLA:HLA00489"], ["P1;HLA:HLA00490"]]

    def test_regular_PirHeaderParser(self):
        for inp, out in zip(self.ins_pir, self.outs_pir):
            handle = StringIO(inp)
            self.assertEqual(list(PirHeaderParser(handle)), out)

    def test_edgecases(self):
        handle = StringIO(">P1;NoDesc\n\nseq\n*\n>P1;NoSeq\ndesc\n*")
        self.assertEqual(list(PirHeaderParser(handle)), ["P1;NoDesc", "P1;NoSeq"])


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
