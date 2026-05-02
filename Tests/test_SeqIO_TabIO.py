# Copyright 2026 by the Biopython contributors. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for SeqIO TabIO module — error paths."""

import unittest

from Bio.Seq import Seq
from Bio.SeqIO.TabIO import TabWriter
from Bio.SeqRecord import SeqRecord


class TestTabIOInvalidInput(unittest.TestCase):
    """Tab format rejects records whose id or sequence contain reserved chars."""

    def test_tab_id_with_tab_raises(self):
        r"""Reject tab-format write when record.id contains a tab character.

        Note: ``_clean`` already replaces ``\n`` and ``\r`` in the id with
        spaces before this check runs, so only the tab branch of the id check
        is reachable from a normal write call.
        """
        record = SeqRecord(Seq("ACGT"), id="bad\tid", description="")
        with self.assertRaisesRegex(ValueError, "Record id contains"):
            TabWriter.to_string(record)

    def test_tab_seq_with_tab_raises(self):
        """Reject tab-format write when the sequence contains a tab character."""
        record = SeqRecord(Seq("AC\tGT"), id="ok_id", description="")
        with self.assertRaisesRegex(ValueError, "Record sequence contains"):
            TabWriter.to_string(record)

    def test_tab_seq_with_newline_raises(self):
        """Reject tab-format write when the sequence contains a newline."""
        record = SeqRecord(Seq("AC\nGT"), id="ok_id", description="")
        with self.assertRaisesRegex(ValueError, "Record sequence contains"):
            TabWriter.to_string(record)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
