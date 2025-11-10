# Copyright 2017 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Testing Bio.ExPASy online code."""

import io
import unittest

import requires_internet

# We want to test these:
from Bio import ExPASy

# In order to check any records returned
from Bio.ExPASy import Prodoc
from Bio.ExPASy import Prosite
from Bio.ExPASy import ScanProsite

requires_internet.check()


class ExPASyOnlineTests(unittest.TestCase):
    """Test ExPASy online resources."""

    def test_prosite_raw(self):
        with ExPASy.get_prosite_raw("PS00001") as handle:
            record = Prosite.read(handle)
        self.assertEqual(record.accession, "PS00001")
        self.assertEqual(record.name, "ASN_GLYCOSYLATION")

    def test_prodoc_raw(self):
        with ExPASy.get_prosite_raw("PDOC00001") as handle:
            record = Prodoc.read(handle)
        self.assertEqual(record.accession, "PDOC00001")

    def test_prosite_html(self):
        with ExPASy.get_prosite_entry("PS00001") as handle:
            html = handle.read()
        self.assertEqual(
            handle.url,
            "https://prosite.expasy.org/cgi-bin/prosite/get-prosite-entry?PS00001",
        )
        self.assertIn("<title>PROSITE - PS00001</title>", html)

    def test_prodoc_html(self):
        with ExPASy.get_prodoc_entry("PDOC00001") as handle:
            html = handle.read()
        self.assertEqual(
            handle.url,
            "https://prosite.expasy.org/cgi-bin/prosite/get-prodoc-entry?PDOC00001",
        )
        self.assertIn("{PS00001; ASN_GLYCOSYLATION}", html)

    def test_scanprosite_swissprot(self):
        pattern = "P-x(2)-G-E-S-G(2)-[AS]"
        result = ScanProsite.scan(
            sig=pattern, mirror=ScanProsite.PROSITE_URL, output="xml", db="sp"
        )
        sequences = ScanProsite.read(result)
        assert len(sequences) > 0

    def test_scanprosite_trembl(self):
        pattern = "P-x(2)-G-E-S-G(2)-[AS]"
        result = ScanProsite.scan(
            sig=pattern, mirror=ScanProsite.PROSITE_URL, output="xml", db="tr"
        )
        sequences = ScanProsite.read(result)
        assert len(sequences) > 0

    def test_scanprosite_pdb(self):
        pattern = "P-x(2)-G-E-S-G(2)-[AS]"
        result = ScanProsite.scan(
            sig=pattern, mirror=ScanProsite.PROSITE_URL, output="xml", db="pdb"
        )
        sequences = ScanProsite.read(result)
        assert len(sequences) > 0

    def test_scanprosite_output_not_implemented(self):
        pattern = "P-x(2)-G-E-S-G(2)-[AS]"
        with self.assertRaises(NotImplementedError):
            ScanProsite.scan(
                sig=pattern, mirror=ScanProsite.PROSITE_URL, output="txt", db="sp"
            )


class ExPASyOfflineTests(unittest.TestCase):
    """Test ExPASy offline with pre-saved data."""

    def test_scanprosite_reads_xml(self):
        with open("ExPASy/scanprosite_response.xml", "rb") as f:
            example_xml = f.read()
            handle = io.BytesIO(example_xml)
            sequences = ScanProsite.read(handle)
            assert len(sequences) == 1081


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
