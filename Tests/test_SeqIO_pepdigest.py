# Copyright 2010 by Thomas Schmitt.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for SeqIO SeqXML module."""
import sys
import unittest

from io import BytesIO

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class TestSimpleRead(unittest.TestCase):
    def test_check_SeqIO(self):
        """Files readable using parser via SeqIO."""
        records = list(SeqIO.parse("pepdigest/small.pepdigest", "pepdigest"))
        self.assertEqual(len(records), 41)
        records = list(SeqIO.parse("pepdigest/big.pepdigest", "pepdigest"))
        self.assertEqual(len(records), 969)


class TestDetailedRead(unittest.TestCase):

    records = {}

    def setUp(self):
        self.records = list(SeqIO.parse("pepdigest/small.pepdigest", "pepdigest"))

    def test_read_aux_information(self):
        self.assertEqual(
            self.records[0].annotations,
            {
                "Cterm": ".",
                "End": 35,
                "Mol_W": 3851.532,
                "Nterm": "S",
                "Start": 1,
                "digestion_enzyme": "Trypsin",
                "source_seq_from": 1,
                "source_seq_to": 243,
                "source_sequence": "GFAH01000435.1.p1",
            },
        )
        self.assertEqual(
            self.records[25].annotations,
            {
                "Cterm": "K",
                "End": 65,
                "Mol_W": 1997.481,
                "Nterm": ".",
                "Start": 49,
                "digestion_enzyme": "Trypsin",
                "source_seq_from": 1,
                "source_seq_to": 65,
                "source_sequence": "GFAH01000035.1.p2",
            },
        )
        self.assertEqual(
            self.records[39].annotations,
            {
                "Cterm": "K",
                "End": 53,
                "Mol_W": 261.278,
                "Nterm": "Q",
                "Start": 52,
                "digestion_enzyme": "Trypsin",
                "source_seq_from": 1,
                "source_seq_to": 141,
                "source_sequence": "GFAH01000730.1.p1",
            },
        )

    def test_read_sequence(self):
        self.assertEqual(
            self.records[0].seq, Seq("MGLIIVLVISVLSADAVLSMDNELYLNLEPHPSQR")
        )
        self.assertEqual(self.records[25].seq, Seq("CKRISSRKLICVCEQCK"))
        self.assertEqual(self.records[39].seq, Seq("DK"))


# class TestReadCorruptFiles(unittest.TestCase):
#     def test_for_errors(self):
#         """Handling of corrupt files."""
#         # SeqIO.parse reads the file in blocks until it finds the seqXML
#         # element with global information such as the source and sourceVersion.
#         # Since one block is likely large enough to cover the first few
#         # entries in the file, the ValueError may be raised after we call
#         # SeqIO.parse, before we start iterating over the file.
#         def f(path):
#             records = SeqIO.parse(path, "seqxml")
#             for record in records:
#                 pass

#         self.assertRaises(ValueError, f, "SeqXML/corrupt_example1.xml")
#         self.assertRaises(ValueError, f, "SeqXML/corrupt_example2.xml")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
