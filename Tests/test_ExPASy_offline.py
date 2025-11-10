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
