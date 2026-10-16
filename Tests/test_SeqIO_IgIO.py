# Copyright 2026 by Contributors.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.SeqIO.IgIO module."""

import unittest
from io import StringIO

from Bio.SeqIO.IgIO import IgHeaderParser


class TestIgHeaderParser(unittest.TestCase):
    """Test IgHeaderParser directly."""

    ins_ig = [
        ";;Free text\n;Comment\nA_U455\nsequence1\n1\n",
        ";;Comment\n;;More\n;Comment\nB_HXB2R\nsequence2\n2\n",
        ";Comment\nC_UG268A\nseq\n1\n",
    ]
    outs_ig = [["A_U455"], ["B_HXB2R"], ["C_UG268A"]]

    def test_regular_IgHeaderParser(self):
        for inp, out in zip(self.ins_ig, self.outs_ig):
            handle = StringIO(inp)
            self.assertEqual(list(IgHeaderParser(handle)), out)

    def test_multiple(self):
        handle = StringIO(
            ";;head\n;comment\nA_U455\nseq\n1\n;comment\nB_HXB2R\nseq2\n2\n"
        )
        self.assertEqual(list(IgHeaderParser(handle)), ["A_U455", "B_HXB2R"])


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
