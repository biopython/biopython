# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for online Entrez access.

This file include tests for accessing the online Entrez service and parsing the
returned results. Note that we are merely testing the access and whether the
results are parseable. Detailed tests on each Entrez service are not within the
scope of this file as they are already covered in test_Entrez.py.

"""
import unittest

import requires_internet
requires_internet.check()

from Bio import Entrez
from Bio import Medline
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


#This lets us set the email address to be sent to NCBI Entrez:
Entrez.email = "biopython-dev@biopython.org"


class EntrezOnlineCase(unittest.TestCase):

    def test_read_from_url(self):
        """Test Entrez.read from URL"""
        einfo = Entrez.einfo()
        rec = Entrez.read(einfo)
        self.assertTrue(isinstance(rec, dict))
        self.assertTrue('DbList' in rec)
        # arbitrary number, just to make sure that DbList has contents
        self.assertTrue(len(rec['DbList']) > 5)

    def test_parse_from_url(self):
        """Test Entrez.parse from URL"""
        efetch = Entrez.efetch(db='protein', id='15718680,157427902,119703751',
                retmode='xml')
        recs = Entrez.parse(efetch)
        recs = list(recs)
        self.assertEqual(3, len(recs))
        # arbitrary number, just to make sure the parser works
        self.assertTrue(all(len(rec).keys > 5) for rec in recs)

    def test_seqio_from_url(self):
        """Test Entrez into SeqIO.read from URL"""
        efetch = Entrez.efetch(db='nucleotide', id='186972394', rettype='gb',
                retmode='text')
        record = SeqIO.read(efetch, 'genbank')
        self.assertTrue(isinstance(record, SeqRecord))
        self.assertEqual('EU490707.1', record.id)
        self.assertEqual(1302, len(record))

    def test_medline_from_url(self):
        """Test Entrez into Medline.read from URL"""
        efetch = Entrez.efetch(db="pubmed", id='19304878', rettype="medline",
                retmode="text")
        record = Medline.read(efetch)
        self.assertTrue(isinstance(record, dict))
        self.assertEqual('19304878', record['PMID'])
        self.assertEqual('10.1093/bioinformatics/btp163 [doi]', record['LID'])


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
