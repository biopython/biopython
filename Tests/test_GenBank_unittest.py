# Copyright 2013 by Kai Blin.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import unittest
from os import path

from Bio import SeqIO


class GenBankTests(unittest.TestCase):
    def test_invalid_product_line_raises_value_error(self):
        "Test GenBank parsing invalid product line raises ValueError"
        def parse_invalid_product_line():
            rec = SeqIO.read(path.join('GenBank', 'invalid_product.gb'),
                             'genbank')
        self.assertRaises(ValueError, parse_invalid_product_line)



if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
