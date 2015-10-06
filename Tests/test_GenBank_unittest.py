# Copyright 2013 by Kai Blin.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import unittest
import warnings
from os import path

from Bio import BiopythonParserWarning
from Bio import GenBank
from Bio import SeqIO


class GenBankTests(unittest.TestCase):
    def test_invalid_product_line_raises_value_error(self):
        """Test GenBank parsing invalid product line raises ValueError"""
        def parse_invalid_product_line():
            rec = SeqIO.read(path.join('GenBank', 'invalid_product.gb'),
                             'genbank')
        self.assertRaises(ValueError, parse_invalid_product_line)

    def test_genbank_read(self):
        with open(path.join("GenBank", "NC_000932.gb")) as handle:
            record = GenBank.read(handle)
        self.assertEqual(['NC_000932'], record.accession)

    def test_genbank_read_multirecord(self):
        with open(path.join("GenBank", "cor6_6.gb")) as handle:
            self.assertRaises(ValueError, GenBank.read, handle)

    def test_genbank_read_invalid(self):
        with open(path.join("GenBank", "NC_000932.faa")) as handle:
            self.assertRaises(ValueError, GenBank.read, handle)

    def test_genbank_read_no_origin_no_end(self):
        with open(path.join("GenBank", "no_origin_no_end.gb")) as handle:
            self.assertRaises(ValueError, GenBank.read, handle)

    # Evil hack with 000 to manipulate sort order to ensure this is tested
    # first (otherwise something silences the warning)
    def test_000_genbank_bad_loc_wrap_warning(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error", BiopythonParserWarning)
            with open(path.join("GenBank", "bad_loc_wrap.gb")) as handle:
                # self.assertRaises(BiopythonParserWarning, GenBank.read, handle)
                try:
                    record = GenBank.read(handle)
                except BiopythonParserWarning as e:
                    self.assertEqual(str(e), "Non-standard feature line wrapping (didn't break on comma)?")
                else:
                    self.assertTrue(False, "Expected specified BiopythonParserWarning here.")

    # Similar hack as we also want to catch that warning here
    def test_001_negative_location_warning(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error", BiopythonParserWarning)
            try:
                SeqIO.read(path.join("GenBank", "negative_location.gb"), "genbank")
            except BiopythonParserWarning as e:
                self.assertEqual(str(e), "Couldn't parse feature location: '-2..492'")
            else:
                self.assertTrue(False, "Expected specified BiopythonParserWarning here.")

    def test_genbank_bad_loc_wrap_parsing(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonParserWarning)
            with open(path.join("GenBank", "bad_loc_wrap.gb")) as handle:
                record = GenBank.read(handle)
                self.assertEqual(1, len(record.features))
                loc = record.features[0].location
                self.assertEqual(loc, "join(3462..3615,3698..3978,4077..4307,4408..4797,4876..5028,5141..5332)")

    def test_negative_location(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonParserWarning)
            rec = SeqIO.read(path.join("GenBank", "negative_location.gb"), "genbank")
            self.assertEqual(None, rec.features[-1].location)

    def test_dot_lineage(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonParserWarning)
            rec = SeqIO.read("GenBank/bad_loc_wrap.gb", "genbank")
        self.assertEqual(rec.annotations["organism"], ".")
        self.assertEqual(rec.annotations["taxonomy"], [])

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
