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
import os
import unittest

import requires_internet
requires_internet.check()

from Bio import Entrez
from Bio import Medline
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

if os.name == 'java':
    try:
        from xml.parsers.expat import XML_PARAM_ENTITY_PARSING_ALWAYS
        del XML_PARAM_ENTITY_PARSING_ALWAYS
    except ImportError:
        from Bio import MissingPythonDependencyError
        raise MissingPythonDependencyError("The Bio.Entrez XML parser fails on "
                                  "Jython, see http://bugs.jython.org/issue1447")


# This lets us set the email address to be sent to NCBI Entrez:
Entrez.email = "biopython-dev@biopython.org"

URL_HEAD = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
URL_TOOL = "tool=biopython"
URL_EMAIL = "email=biopython-dev%40biopython.org"


class EntrezOnlineCase(unittest.TestCase):

    def test_read_from_url(self):
        """Test Entrez.read from URL"""
        handle = Entrez.einfo()
        self.assertTrue(handle.url.startswith(URL_HEAD + "einfo.fcgi?"), handle.url)
        self.assertTrue(URL_TOOL in handle.url)
        self.assertTrue(URL_EMAIL in handle.url)
        rec = Entrez.read(handle)
        handle.close()
        self.assertTrue(isinstance(rec, dict))
        self.assertTrue('DbList' in rec)
        # arbitrary number, just to make sure that DbList has contents
        self.assertTrue(len(rec['DbList']) > 5)

    def test_parse_from_url(self):
        """Test Entrez.parse from URL"""
        handle = Entrez.efetch(db='protein', id='15718680,157427902,119703751',
                               retmode='xml')
        self.assertTrue(handle.url.startswith(URL_HEAD + "efetch.fcgi?"), handle.url)
        self.assertTrue(URL_TOOL in handle.url)
        self.assertTrue(URL_EMAIL in handle.url)
        self.assertTrue("id=15718680%2C157427902%2C119703751" in handle.url, handle.url)
        recs = list(Entrez.parse(handle))
        handle.close()
        self.assertEqual(3, len(recs))
        # arbitrary number, just to make sure the parser works
        self.assertTrue(all(len(rec).keys > 5) for rec in recs)

    def test_webenv_search(self):
        """Test Entrez.search from link webenv history"""
        handle = Entrez.elink(db='nucleotide', dbfrom='protein',
                              id='22347800,48526535', webenv=None, query_key=None,
                              cmd='neighbor_history')
        self.assertTrue(handle.url.startswith(URL_HEAD + "elink.fcgi?"), handle.url)
        self.assertTrue(URL_TOOL in handle.url)
        self.assertTrue(URL_EMAIL in handle.url)
        self.assertTrue("id=22347800%2C48526535" in handle.url, handle.url)
        recs = Entrez.read(handle)
        handle.close()
        record = recs.pop()

        webenv = record['WebEnv']
        query_key = record['LinkSetDbHistory'][0]['QueryKey']
        handle = Entrez.esearch(db='nucleotide', term=None,
                                retstart=0, retmax=10,
                                webenv=webenv, query_key=query_key,
                                usehistory='y')
        self.assertTrue(handle.url.startswith(URL_HEAD + "esearch.fcgi?"), handle.url)
        self.assertTrue(URL_TOOL in handle.url)
        self.assertTrue(URL_EMAIL in handle.url)
        search_record = Entrez.read(handle)
        handle.close()
        self.assertEqual(2, len(search_record['IdList']))

    def test_seqio_from_url(self):
        """Test Entrez into SeqIO.read from URL"""
        handle = Entrez.efetch(db='nucleotide', id='186972394', rettype='gb',
                               retmode='text')
        self.assertTrue(handle.url.startswith(URL_HEAD + "efetch.fcgi?"), handle.url)
        self.assertTrue(URL_TOOL in handle.url)
        self.assertTrue(URL_EMAIL in handle.url)
        self.assertTrue("id=186972394" in handle.url)
        record = SeqIO.read(handle, 'genbank')
        handle.close()
        self.assertTrue(isinstance(record, SeqRecord))
        self.assertEqual('EU490707.1', record.id)
        self.assertEqual(1302, len(record))

    def test_medline_from_url(self):
        """Test Entrez into Medline.read from URL"""
        handle = Entrez.efetch(db="pubmed", id='19304878', rettype="medline",
                               retmode="text")
        self.assertTrue(handle.url.startswith(URL_HEAD + "efetch.fcgi?"), handle.url)
        self.assertTrue(URL_TOOL in handle.url)
        self.assertTrue(URL_EMAIL in handle.url)
        self.assertTrue("id=19304878" in handle.url)
        record = Medline.read(handle)
        handle.close()
        self.assertTrue(isinstance(record, dict))
        self.assertEqual('19304878', record['PMID'])
        self.assertEqual('10.1093/bioinformatics/btp163 [doi]', record['LID'])

    def test_elink(self):
        # Commas: Link from protein to gene
        handle = Entrez.elink(db="gene", dbfrom="protein",
                              id="15718680,157427902,119703751")
        self.assertTrue(handle.url.startswith(URL_HEAD + "elink.fcgi"), handle.url)
        self.assertTrue(URL_TOOL in handle.url)
        self.assertTrue(URL_EMAIL in handle.url)
        self.assertTrue("id=15718680%2C157427902%2C119703751" in handle.url, handle.url)
        handle.close()

        # Multiple ID entries: Find one-to-one links from protein to gene
        handle = Entrez.elink(db="gene", dbfrom="protein",
                              id=["15718680", "157427902", "119703751"])
        self.assertTrue(handle.url.startswith(URL_HEAD + "elink.fcgi"), handle.url)
        self.assertTrue(URL_TOOL in handle.url)
        self.assertTrue(URL_EMAIL in handle.url)
        self.assertTrue("id=15718680" in handle.url, handle.url)
        self.assertTrue("id=157427902" in handle.url, handle.url)
        self.assertTrue("id=119703751" in handle.url, handle.url)
        handle.close()

    def test_epost(self):
        handle = Entrez.epost("nuccore", id="186972394,160418")
        self.assertEqual(URL_HEAD + "epost.fcgi", handle.url)
        handle.close()
        handle = Entrez.epost("nuccore", id=["160418", "160351"])
        self.assertEqual(URL_HEAD + "epost.fcgi", handle.url)
        handle.close()


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
