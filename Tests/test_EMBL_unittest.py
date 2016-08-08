# Copyright 2015 by Kai Blin.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import unittest
import warnings
from os import path

from Bio import SeqIO


class EMBLTests(unittest.TestCase):
    def test_embl_content_after_co(self):
        """Test an AssertionError is thrown by content after a CO line"""
        def parse_content_after_co():
            rec = SeqIO.read(path.join('EMBL', 'xx_after_co.embl'), 'embl')

        self.assertRaises(AssertionError, parse_content_after_co)

        try:
            parse_content_after_co()
        except AssertionError as e:
            self.assertEqual(str(e), "Unexpected content after SQ or CO line: 'XX'")
        else:
            self.assertTrue(False, "Error message without explanation raised by content after CO line")

    def test_embl_0_line(self):
        """Test an Assertion is thrown by SQ line with 0 length sequence"""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            SeqIO.read(path.join('EMBL', 'embl_with_0_line.embl'), 'embl')
            assert len(w) == 0, "Importing embl format error: sequence line with no seuqence but coordinate"


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
