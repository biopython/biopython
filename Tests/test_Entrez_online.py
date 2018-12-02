# -*- coding: utf-8 -*-
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
from Bio import Entrez
from Bio import Medline
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import doctest
import os
import locale
import sys
import unittest

import requires_internet
requires_internet.check()

if os.name == 'java':
    try:
        from xml.parsers.expat import XML_PARAM_ENTITY_PARSING_ALWAYS
        del XML_PARAM_ENTITY_PARSING_ALWAYS
    except ImportError:
        from Bio import MissingPythonDependencyError
        raise MissingPythonDependencyError("The Bio.Entrez XML parser fails on "
                                           "Jython, see http://bugs.jython.org/issue1447")


# This lets us set the email address to be sent to NCBI Entrez:
Entrez.email = "biopython@biopython.org"
Entrez.api_key = "5cfd4026f9df285d6cfc723c662d74bcbe09"

URL_HEAD = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
URL_TOOL = "tool=biopython"
URL_EMAIL = "email=biopython%40biopython.org"
URL_API_KEY = "api_key=5cfd4026f9df285d6cfc723c662d74bcbe09"


class EntrezOnlineCase(unittest.TestCase):

    def test_no_api_key(self):
        """Test Entrez.read without API key."""
        cached = Entrez.api_key
        Entrez.api_key = None  # default
        try:
            handle = Entrez.einfo()
        finally:
            # Do not want any failure here to break other tests
            Entrez.api_key = cached
        self.assertTrue(handle.url.startswith(URL_HEAD + "einfo.fcgi?"), handle.url)
        self.assertIn(URL_TOOL, handle.url)
        self.assertIn(URL_EMAIL, handle.url)
        self.assertNotIn("api_key=", handle.url)
        rec = Entrez.read(handle)
        handle.close()
        self.assertTrue(isinstance(rec, dict))
        self.assertIn('DbList', rec)
        # arbitrary number, just to make sure that DbList has contents
        self.assertTrue(len(rec['DbList']) > 5)

    def test_read_from_url(self):
        """Test Entrez.read from URL"""
        handle = Entrez.einfo()
        self.assertTrue(handle.url.startswith(URL_HEAD + "einfo.fcgi?"), handle.url)
        self.assertIn(URL_TOOL, handle.url)
        self.assertIn(URL_EMAIL, handle.url)
        self.assertIn(URL_API_KEY, handle.url)
        rec = Entrez.read(handle)
        handle.close()
        self.assertTrue(isinstance(rec, dict))
        self.assertIn('DbList', rec)
        # arbitrary number, just to make sure that DbList has contents
        self.assertTrue(len(rec['DbList']) > 5)

    def test_parse_from_url(self):
        """Test Entrez.parse from URL"""
        handle = Entrez.efetch(db='protein', id='15718680,157427902,119703751',
                               retmode='xml')
        self.assertTrue(handle.url.startswith(URL_HEAD + "efetch.fcgi?"), handle.url)
        self.assertIn(URL_TOOL, handle.url)
        self.assertIn(URL_EMAIL, handle.url)
        self.assertIn(URL_API_KEY, handle.url)
        self.assertIn("id=15718680%2C157427902%2C119703751", handle.url)
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
        self.assertIn(URL_TOOL, handle.url)
        self.assertIn(URL_EMAIL, handle.url)
        self.assertIn(URL_API_KEY, handle.url)
        self.assertIn("id=22347800%2C48526535", handle.url)
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
        self.assertIn(URL_TOOL, handle.url)
        self.assertIn(URL_EMAIL, handle.url)
        self.assertIn(URL_API_KEY, handle.url)
        search_record = Entrez.read(handle)
        handle.close()
        self.assertEqual(2, len(search_record['IdList']))

    def test_seqio_from_url(self):
        """Test Entrez into SeqIO.read from URL"""
        handle = Entrez.efetch(db='nucleotide', id='186972394', rettype='gb',
                               retmode='text')
        self.assertTrue(handle.url.startswith(URL_HEAD + "efetch.fcgi?"), handle.url)
        self.assertIn(URL_TOOL, handle.url)
        self.assertIn(URL_EMAIL, handle.url)
        self.assertIn(URL_API_KEY, handle.url)
        self.assertIn("id=186972394", handle.url)
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
        self.assertIn(URL_TOOL, handle.url)
        self.assertIn(URL_EMAIL, handle.url)
        self.assertIn(URL_API_KEY, handle.url)
        self.assertIn("id=19304878", handle.url)
        record = Medline.read(handle)
        handle.close()
        self.assertTrue(isinstance(record, dict))
        self.assertEqual('19304878', record['PMID'])
        self.assertEqual('10.1093/bioinformatics/btp163 [doi]', record['LID'])

    def test_efetch_biosystems_xml(self):
        """Test Entrez parser with XML from biosystems"""
        handle = Entrez.efetch(id="1134002", db="biosystems", retmode="xml")
        records = list(Entrez.parse(handle))
        handle.close()
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]['System_sysid']['Sys-id']['Sys-id_bsid'], '1134002')

    def test_efetch_taxonomy_xml(self):
        """Test Entrez using a integer id - like a taxon id"""
        handle = Entrez.efetch(db="taxonomy", id=3702, retmode="XML")
        taxon_record = Entrez.read(handle)
        self.assertTrue(1, len(taxon_record))
        self.assertIn('TaxId', taxon_record[0])
        self.assertTrue('3702', taxon_record[0]['TaxId'])
        handle.close()

    def test_elink(self):
        """Test Entrez.elink with multiple ids, both comma separated and as list"""
        # Commas: Link from protein to gene
        handle = Entrez.elink(db="gene", dbfrom="protein",
                              id="15718680,157427902,119703751")
        self.assertTrue(handle.url.startswith(URL_HEAD + "elink.fcgi"), handle.url)
        self.assertIn(URL_TOOL, handle.url)
        self.assertIn(URL_EMAIL, handle.url)
        self.assertIn(URL_API_KEY, handle.url)
        self.assertIn("id=15718680%2C157427902%2C119703751", handle.url)
        handle.close()

        # Multiple ID entries: Find one-to-one links from protein to gene
        handle = Entrez.elink(db="gene", dbfrom="protein",
                              id=["15718680", "157427902", "119703751"])
        self.assertTrue(handle.url.startswith(URL_HEAD + "elink.fcgi"), handle.url)
        self.assertIn(URL_TOOL, handle.url)
        self.assertIn(URL_EMAIL, handle.url)
        self.assertIn(URL_API_KEY, handle.url)
        self.assertIn("id=15718680", handle.url)
        self.assertIn("id=157427902", handle.url)
        self.assertIn("id=119703751", handle.url)
        handle.close()

    def test_epost(self):
        """Test Entrez.epost with multiple ids, both comma separated and as list"""
        handle = Entrez.epost("nuccore", id="186972394,160418")
        self.assertEqual(URL_HEAD + "epost.fcgi", handle.url)
        handle.close()
        handle = Entrez.epost("nuccore", id=["160418", "160351"])
        self.assertEqual(URL_HEAD + "epost.fcgi", handle.url)
        handle.close()

    def test_egquery(self):
        """Test Entrez.egquery which searches in all Entrez databases for a single text query"""
        handle = Entrez.egquery(term="biopython")
        record = Entrez.read(handle)
        handle.close()

        done = False
        for row in record["eGQueryResult"]:
            if "pmc" in row["DbName"]:
                self.assertTrue(int(row["Count"]) > 60)
                done = True
        self.assertTrue(done)

    def test_espell(self):
        """Test misspellings with Entrez.espell"""
        handle = Entrez.espell(term="biopythooon")
        record = Entrez.read(handle)
        handle.close()

        self.assertEqual(record["Query"], "biopythooon")
        self.assertEqual(record["CorrectedQuery"], "biopython")

    def test_ecitmatch(self):
        """Test Entrez.ecitmatch to search for a citation"""
        citation = {
            "journal_title": "proc natl acad sci u s a",
            "year": "1991", "volume": "88", "first_page": "3248",
            "author_name": "mann bj", "key": "citation_1"
        }
        handle = Entrez.ecitmatch(db="pubmed", bdata=[citation])
        self.assertIn("retmode=xml", handle.url)
        result = handle.read()
        expected_result = "proc natl acad sci u s a|1991|88|3248|mann bj|citation_1|2014248\n"
        self.assertEquals(result, expected_result)
        handle.close()

    def test_efetch_gds_utf8(self):
        """Test correct handling of encodings in Entrez.efetch"""
        # See issue #1402
        try:
            oldloc = locale.setlocale(locale.LC_ALL)
            locale.setlocale(locale.LC_ALL, 'C')
        except locale.Error as err:
            self.skipTest("Cannot set locale: {}".format(err))
        try:
            handle = Entrez.efetch(db="gds", id='200079209')
            self.assertTrue(handle.url.startswith(URL_HEAD + "efetch.fcgi?"), handle.url)
            self.assertIn(URL_TOOL, handle.url)
            self.assertIn(URL_EMAIL, handle.url)
            self.assertIn(URL_API_KEY, handle.url)
            self.assertIn("id=200079209", handle.url)
            result = handle.read()
            if sys.version_info[0] < 3:
                result = result.decode("UTF8")
            # Use of Unicode double quotation marks U+201C and U+201D
            expected_result = u'“field of injury”'
            self.assertEqual(result[342:359], expected_result)
            handle.close()
        finally:
            locale.setlocale(locale.LC_ALL, oldloc)

    def test_fetch_xml_schemas(self):
        handle = Entrez.efetch("protein", id="783730874", rettype="ipg", retmode="xml")
        records = list(Entrez.parse(handle, validate=False))
        handle.close()
        self.assertEqual(len(records), 1)
        self.assertIn("Product", records[0])
        self.assertIn("Statistics", records[0])
        self.assertIn("ProteinList", records[0])


if __name__ == "__main__":
    # When running test_Entrez.py directly, will also include the
    # Bio.Entrez doctests.
    # TODO: Include the doctests via run_tests.py when online.
    unittest_suite = unittest.TestLoader().loadTestsFromName("test_Entrez_online")
    doctest_suite = doctest.DocTestSuite(Entrez)
    suite = unittest.TestSuite((unittest_suite, doctest_suite))
    runner = unittest.TextTestRunner(sys.stdout, verbosity=2)
    runner.run(suite)
