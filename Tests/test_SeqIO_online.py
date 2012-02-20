# Copyright 2007-2010 by Peter Cock.  All rights reserved.
# Revisions copyright 2007-2008 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Testing online code for fetching sequences, and parsing them

Uses Bio.SeqIO to parse files downloaded with Bio.GenBank, Bio.WWW.NCBI, 
Bio.ExPASy etc.

Goals:
    Make sure that all retrieval is working as expected.
    May catch some format changes early too.
"""
import unittest

import requires_internet
requires_internet.check()

from Bio import MissingExternalDependencyError

#We want to test these:
from Bio import Entrez
from Bio import ExPASy

#In order to check any sequences returned
from Bio import SeqIO
from StringIO import StringIO
from Bio.SeqUtils.CheckSum import seguid

from Bio.File import UndoHandle
from Bio._py3k import _as_string

#This lets us set the email address to be sent to NCBI Entrez:
Entrez.email = "biopython-dev@biopython.org"

class ExPASyTests(unittest.TestCase):
    """Tests for Bio.ExPASy module."""
    def test_get_sprot_raw(self):
        """Bio.ExPASy.get_sprot_raw("O23729")"""
        identifier = "O23729"
        try:
            #This is to catch an error page from our proxy:
            handle = UndoHandle(ExPASy.get_sprot_raw(identifier))
            if _as_string(handle.peekline()).startswith("<!DOCTYPE HTML"):
                raise IOError
            record = SeqIO.read(handle, "swiss")
            handle.close()
        except IOError:
            raise MissingExternalDependencyError(
                  "internet (or maybe just ExPASy) not available")
        self.assertEqual(record.id, identifier)
        self.assertEqual(len(record), 394)
        self.assertEqual(seguid(record.seq), "5Y08l+HJRDIlhLKzFEfkcKd1dkM")

class EntrezTests(unittest.TestCase):
    def simple(self, database, formats, entry, length, checksum):
        for f in formats:
            try:
                handle = Entrez.efetch(db=database, id=entry, rettype=f, retmode="text")
                record = SeqIO.read(handle, f)
                handle.close()
            except IOError:
                raise MissingExternalDependencyError(
                      "internet (or maybe just NCBI) not available")
            self.assertTrue((entry in record.name) or \
                         (entry in record.id) or \
                         ("gi" in record.annotations \
                          and record.annotations["gi"]==entry),
                         "%s got %s, %s" % (entry, record.name, record.id))
            self.assertEqual(len(record), length)
            self.assertEqual(seguid(record.seq), checksum)

for database, formats, entry, length, checksum in [
    ("nuccore", ["fasta", "gb"], "X52960", 248,
     "Ktxz0HgMlhQmrKTuZpOxPZJ6zGU"),
    ("nucleotide", ["fasta", "gb"], "6273291", 902,
     "bLhlq4mEFJOoS9PieOx4nhGnjAQ"),
    ("protein", ["fasta", "gb"], "16130152", 367,
     "fCjcjMFeGIrilHAn6h+yju267lg"),
    ]:
    def funct(d, f, e, l, c):
        method = lambda x : x.simple(d, f, e, l, c)
        method.__doc__ = "Bio.Entrez.efetch(%s, %s, ...)" % (d, e)
        return method
    setattr(EntrezTests, "test_%s_%s" % (database, entry),
            funct(database, formats, entry, length, checksum))
    del funct        
del database, formats, entry, length, checksum

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)

    
