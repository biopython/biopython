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

try:
    # Included in standard lib as of Python 3.3
    from unittest.mock import patch
except ImportError:
    try:
        from mock import patch
    except ImportError:
        from Bio import MissingPythonDependencyError
        raise MissingPythonDependencyError("Install the mock library backport on old versions of Python")

from Bio._py3k import _binary_to_string_handle
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
    @patch("Bio.Entrez._open", return_value=_binary_to_string_handle(open("Entrez/einfo1.xml", "rb")))
    def test_read_from_url(self, mock_open):
        """Test Entrez.read from URL"""
        handle = Entrez.einfo()
        rec = Entrez.read(handle)
        handle.close()
        self.assertTrue(isinstance(rec, dict))
        self.assertTrue('DbList' in rec)
        # arbitrary number, just to make sure that DbList has contents
        self.assertTrue(len(rec['DbList']) > 5)

    @patch("Bio.Entrez._open", return_value=_binary_to_string_handle(open("Entrez/protein2.xml", "rb")))
    def test_parse_from_url(self, mock_open):
        """Test Entrez.parse from URL"""
        handle = Entrez.efetch(db='protein', id='15718680,157427902,119703751',
                               retmode='xml')
        recs = list(Entrez.parse(handle))
        handle.close()
        self.assertEqual(3, len(recs))
        # arbitrary number, just to make sure the parser works
        self.assertTrue(all(len(rec).keys > 5) for rec in recs)

    @patch("Bio.Entrez._open", return_value=_binary_to_string_handle(open("Entrez/esearch9.xml", "rb")))
    def test_webenv_search(self, mock_open):
        """Test Entrez.search from link webenv history"""
        webenv = 'NCID_1_214573425_130.14.18.34_9001_1448030770_1520095147_0MetA0_S_MegaStore_F_1'
        query_key = 2
        handle = Entrez.esearch(db='nucleotide', term=None,
                                retstart=0, retmax=10,
                                webenv=webenv, query_key=query_key,
                                usehistory='y')
        search_record = Entrez.read(handle)
        handle.close()
        self.assertEqual(2, len(search_record['IdList']))

    @patch("Bio.Entrez._open", return_value=_binary_to_string_handle(open("Entrez/efetch1.txt", "rb")))
    def test_seqio_from_url(self, mock_open):
        """Test Entrez into SeqIO.read from URL"""
        handle = Entrez.efetch(db='nucleotide', id='186972394', rettype='gb',
                               retmode='text')
        record = SeqIO.read(handle, 'genbank')
        handle.close()
        self.assertTrue(isinstance(record, SeqRecord))
        self.assertEqual('EU490707.1', record.id)
        self.assertEqual(1302, len(record))

    @patch("Bio.Entrez._open", return_value=_binary_to_string_handle(open("Entrez/efetch2.txt", "rb")))
    def test_medline_from_url(self, mock_open):
        """Test Entrez into Medline.read from URL"""
        handle = Entrez.efetch(db="pubmed", id='19304878', rettype="medline",
                               retmode="text")
        record = Medline.read(handle)
        handle.close()
        self.assertTrue(isinstance(record, dict))
        self.assertEqual('19304878', record['PMID'])
        self.assertEqual('10.1093/bioinformatics/btp163 [doi]', record['LID'])

    @patch("Bio.Entrez._open", return_value=_binary_to_string_handle(open("Entrez/efetch3.xml", "rb")))
    def test_efetch_biosystems_xml(self, mock_open):
        """Test Entrez parser with XML from biosystems"""
        handle = Entrez.efetch(id="1134002", db="biosystems", retmode="xml")
        records = list(Entrez.parse(handle))
        handle.close()
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]['System_sysid']['Sys-id']['Sys-id_bsid'], '1134002')

    @patch("Bio.Entrez._open", return_value=_binary_to_string_handle(open("Entrez/ecitmatch.txt", "rb")))
    def test_ecitmatch(self, mock_open):
        citation = {
            "journal_title": "proc natl acad sci u s a",
            "year": "1991", "volume": "88", "first_page": "3248",
            "author_name": "mann bj", "key": "citation_1"
        }
        handle = Entrez.ecitmatch(db="pubmed", bdata=[citation])
        result = handle.read()
        expected_result = "proc natl acad sci u s a|1991|88|3248|mann bj|citation_1|2014248\n"
        self.assertEquals(result, expected_result)

    @patch("Bio.Entrez._open", return_value=_binary_to_string_handle(open("Entrez/efetch4.xml", "rb")))
    def test_fetch_xml_schemas(self, mock_open):
        handle = Entrez.efetch("protein", id="783730874", rettype="ipg", retmode="xml")
        records = list(Entrez.parse(handle))
        handle.close()
        self.assertEqual(len(records), 1)
        self.assertTrue("Product" in records[0])
        self.assertTrue("Statistics" in records[0])
        self.assertTrue("RedundantGiList" in records[0])


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
