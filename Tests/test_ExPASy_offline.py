# Copyright 2025 by Timothy Dennis. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Testing Bio.ExPASy offline code."""

import io
import unittest

from Bio.ExPASy import ScanProsite


class ExPASyOfflineTests(unittest.TestCase):
    """Test ExPASy offline with pre-saved data."""

    def test_scanprosite_reads_xml(self):
        with open("ExPASy/scanprosite_response.xml", "rb") as f:
            example_xml = f.read()
            handle = io.BytesIO(example_xml)
            sequences = ScanProsite.read(handle)
            assert len(sequences) == 1081

    def test_scanprosite_mismatched_end_element(self):
        """Raise ValueError when an end element does not match its start.

        The check guards a defensive invariant in the SAX content handler;
        well-formed XML can never trip it via the parser, so we drive the
        handler directly with a mismatched end-element name.
        """
        content_handler = ScanProsite.ContentHandler()
        content_handler.startElement("scanprosite_response", {})
        content_handler.startElement("matchset", {"n_match": "0", "n_seq": "0"})
        with self.assertRaisesRegex(ValueError, "Unexpected XML end element"):
            content_handler.endElement("not_matchset")
