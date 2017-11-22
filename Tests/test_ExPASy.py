# Copyright 2017 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Testing Bio.ExPASy online code."""

import unittest

import requires_internet
requires_internet.check()

# We want to test these:
from Bio import ExPASy

# In order to check any records returned
from Bio.ExPASy import Prodoc
from Bio.ExPASy import Prosite


# TODO - Use with statement when drop Python 2

class ExPASyOnlineTests(unittest.TestCase):
    """Test ExPASy online resources."""

    def test_prosite_raw(self):
        handle = ExPASy.get_prosite_raw('PS00001')
        record = Prosite.read(handle)
        handle.close()
        self.assertEqual(record.accession, 'PS00001')
        self.assertEqual(record.name, 'ASN_GLYCOSYLATION')

    def test_prodoc_raw(self):
        handle = ExPASy.get_prosite_raw('PDOC00001')
        record = Prodoc.read(handle)
        handle.close()
        self.assertEqual(record.accession, 'PDOC00001')

    def test_prosite_html(self):
        handle = ExPASy.get_prosite_entry('PS00001')
        html = handle.read()
        self.assertEqual(handle.url,
                         'https://prosite.expasy.org/cgi-bin/prosite/get-prosite-entry?PS00001')
        handle.close()
        self.assertTrue('<title>PROSITE: PS00001</title>' in html)

    def test_prodoc_html(self):
        handle = ExPASy.get_prodoc_entry('PDOC00001')
        html = handle.read()
        self.assertEqual(handle.url,
                         'https://prosite.expasy.org/cgi-bin/prosite/get-prodoc-entry?PDOC00001')
        handle.close()
        self.assertTrue('{PS00001; ASN_GLYCOSYLATION}' in html)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
