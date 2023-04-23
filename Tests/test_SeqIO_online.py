# Copyright 2007-2010 by Peter Cock.  All rights reserved.
# Revisions copyright 2007-2008 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Testing online code for fetching sequences, and parsing them.

Uses Bio.SeqIO to parse files downloaded with Bio.GenBank, Bio.WWW.NCBI,
Bio.ExPASy etc.

Goals:
    - Make sure that all retrieval is working as expected.
    - May catch some format changes early too.

"""
import unittest

import requires_internet

from Bio import Entrez  # Testing this
from Bio import ExPASy  # Testing this
from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid
from Bio.SwissProt import SwissProtParserError

requires_internet.check()

# This lets us set the email address to be sent to NCBI Entrez:
Entrez.email = "biopython@biopython.org"


class ExPASyTests(unittest.TestCase):
    """Tests for Bio.ExPASy module."""

    def test_get_sprot_raw(self):
        """Bio.ExPASy.get_sprot_raw("O23729")."""
        identifier = "O23729"
        handle = ExPASy.get_sprot_raw(identifier)
        try:
            record = SeqIO.read(handle, "swiss")
        except SwissProtParserError as e:
            # This is to catch an error page from our proxy
            if str(e) == "Failed to find ID in first line" and e.line.startswith(
                "<!DOCTYPE HTML"
            ):
                raise OSError from None
        handle.close()
        self.assertEqual(record.id, identifier)
        self.assertEqual(len(record), 394)
        self.assertEqual(seguid(record.seq), "5Y08l+HJRDIlhLKzFEfkcKd1dkM")


class EntrezTests(unittest.TestCase):
    def simple(self, database, formats, entry, length, checksum):
        for f in formats:
            handle = Entrez.efetch(db=database, id=entry, rettype=f, retmode="text")
            if f == "gbwithparts":
                f = "gb"
            record = SeqIO.read(handle, f)
            handle.close()
            # NCBI still takes GI on input, but phasing it out in output
            gi_to_acc = {"6273291": "AF191665.1", "16130152": "NP_416719.1"}
            if entry in gi_to_acc:
                entry = gi_to_acc[entry]
            self.assertTrue(
                (entry in record.name)
                or (entry in record.id)
                or ("gi" in record.annotations and record.annotations["gi"] == entry),
                f"{entry} got {record.name}, {record.id}",
            )
            self.assertEqual(len(record), length)
            self.assertEqual(seguid(record.seq), checksum)


for database, formats, entry, length, checksum in [
    ("nuccore", ["fasta", "gb"], "X52960", 248, "Ktxz0HgMlhQmrKTuZpOxPZJ6zGU"),
    ("nucleotide", ["fasta", "gb"], "6273291", 902, "bLhlq4mEFJOoS9PieOx4nhGnjAQ"),
    (
        "protein",
        ["fasta", "gbwithparts"],
        "16130152",
        367,
        "fCjcjMFeGIrilHAn6h+yju267lg",
    ),
]:

    def funct(d, f, e, l, c):  # noqa: E741
        method = lambda x: x.simple(d, f, e, l, c)  # noqa: E731
        method.__doc__ = f"Bio.Entrez.efetch({d!r}, id={e!r}, ...)"
        return method

    setattr(
        EntrezTests,
        f"test_{database}_{entry}",
        funct(database, formats, entry, length, checksum),
    )
    del funct
del database, formats, entry, length, checksum

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
