# Copyright 2015 by Carlos Pena.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Offline tests for two Entrez features.

(1) the URL construction of NCBI's Entrez services.
(2) setting a custom directory for DTD and XSD downloads.
"""

import unittest
import warnings

from Bio import Entrez
from Bio.Entrez import Parser


# This lets us set the email address to be sent to NCBI Entrez:
Entrez.email = "biopython@biopython.org"
Entrez.api_key = "5cfd4026f9df285d6cfc723c662d74bcbe09"

URL_HEAD = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
URL_TOOL = "tool=biopython"
URL_EMAIL = "email=biopython%40biopython.org"
URL_API_KEY = "api_key=5cfd4026f9df285d6cfc723c662d74bcbe09"


class TestURLConstruction(unittest.TestCase):
    def test_email_warning(self):
        """Test issuing warning when user does not specify email address."""
        Entrez.email = None

        with warnings.catch_warnings(record=True) as w:
            Entrez._construct_params(params=None)
            self.assertEqual(len(w), 1)

    def test_construct_cgi_ecitmatch(self):
        citation = {
            "journal_title": "proc natl acad sci u s a",
            "year": "1991", "volume": "88", "first_page": "3248",
            "author_name": "mann bj", "key": "citation_1"
        }
        cgi = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/ecitmatch.cgi'
        variables = Entrez._update_ecitmatch_variables({'db': 'pubmed',
                                                        'bdata': [citation]})
        post = False

        params = Entrez._construct_params(variables)
        options = Entrez._encode_options(ecitmatch=True, params=params)
        result_url = Entrez._construct_cgi(cgi, post=post, options=options)
        self.assertIn("retmode=xml", result_url)
        self.assertIn(URL_API_KEY, result_url)

    def test_construct_cgi_einfo(self):
        """Test constructed url for request to Entrez."""
        cgi = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi'
        params = Entrez._construct_params(params=None)
        options = Entrez._encode_options(ecitmatch=False, params=params)
        result_url = Entrez._construct_cgi(cgi, post=False, options=options)
        self.assertTrue(result_url.startswith(URL_HEAD + "einfo.fcgi?"),
                        result_url)
        self.assertIn(URL_TOOL, result_url)
        self.assertIn(URL_EMAIL, result_url)

    def test_construct_cgi_epost1(self):
        cgi = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi'
        variables = {'db': 'nuccore', 'id': '186972394,160418'}
        post = True

        params = Entrez._construct_params(variables)
        options = Entrez._encode_options(ecitmatch=False, params=params)
        result_url = Entrez._construct_cgi(cgi, post=post, options=options)
        self.assertEqual(URL_HEAD + "epost.fcgi", result_url)

    def test_construct_cgi_epost2(self):
        cgi = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi'
        variables = {'db': 'nuccore', 'id': ["160418", "160351"]}
        post = True

        params = Entrez._construct_params(variables)
        options = Entrez._encode_options(ecitmatch=False, params=params)
        result_url = Entrez._construct_cgi(cgi, post=post, options=options)
        self.assertEqual(URL_HEAD + "epost.fcgi", result_url)

    def test_construct_cgi_elink1(self):
        cgi = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi'
        variables = {'cmd': 'neighbor_history', 'db': 'nucleotide',
                     'dbfrom': 'protein', 'id': '22347800,48526535',
                     'query_key': None, 'webenv': None}
        post = False

        params = Entrez._construct_params(variables)
        options = Entrez._encode_options(ecitmatch=False, params=params)
        result_url = Entrez._construct_cgi(cgi, post=post, options=options)
        self.assertTrue(result_url.startswith(URL_HEAD + "elink.fcgi?"),
                        result_url)
        self.assertIn(URL_TOOL, result_url)
        self.assertIn(URL_EMAIL, result_url)
        self.assertIn("id=22347800%2C48526535", result_url)
        self.assertIn(URL_API_KEY, result_url)

    def test_construct_cgi_elink2(self):
        """Commas: Link from protein to gene."""
        cgi = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi'
        variables = {'db': 'gene', 'dbfrom': 'protein',
                     'id': '15718680,157427902,119703751'}
        post = False

        params = Entrez._construct_params(variables)
        options = Entrez._encode_options(ecitmatch=False, params=params)
        result_url = Entrez._construct_cgi(cgi, post=post, options=options)
        self.assertTrue(result_url.startswith(URL_HEAD + "elink.fcgi"),
                        result_url)
        self.assertIn(URL_TOOL, result_url)
        self.assertIn(URL_EMAIL, result_url)
        self.assertTrue("id=15718680%2C157427902%2C119703751" in result_url,
                        result_url)
        self.assertIn(URL_API_KEY, result_url)

    def test_construct_cgi_elink3(self):
        """Multiple ID entries: Find one-to-one links from protein to gene."""
        cgi = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi'
        variables = {'db': 'gene', 'dbfrom': 'protein',
                     'id': ["15718680", "157427902", "119703751"]}
        post = False

        params = Entrez._construct_params(variables)
        options = Entrez._encode_options(ecitmatch=False, params=params)
        result_url = Entrez._construct_cgi(cgi, post=post, options=options)
        self.assertTrue(result_url.startswith(URL_HEAD + "elink.fcgi"),
                        result_url)
        self.assertIn(URL_TOOL, result_url)
        self.assertIn(URL_EMAIL, result_url)
        self.assertIn("id=15718680", result_url)
        self.assertIn("id=157427902", result_url)
        self.assertIn("id=119703751", result_url)
        self.assertIn(URL_API_KEY, result_url)

    def test_construct_cgi_efetch(self):
        cgi = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
        variables = {'db': 'protein', 'id': '15718680,157427902,119703751',
                     'retmode': 'xml'}
        post = False

        params = Entrez._construct_params(variables)
        options = Entrez._encode_options(ecitmatch=False, params=params)
        result_url = Entrez._construct_cgi(cgi, post=post, options=options)
        self.assertTrue(result_url.startswith(URL_HEAD + "efetch.fcgi?"),
                        result_url)
        self.assertIn(URL_TOOL, result_url)
        self.assertIn(URL_EMAIL, result_url)
        self.assertTrue("id=15718680%2C157427902%2C119703751" in result_url,
                        result_url)
        self.assertIn(URL_API_KEY, result_url)


class CustomDirectoryTest(unittest.TestCase):
    """Offline unit test for custom directory feature.

    Allow user to specify a custom directory for Entrez DTD/XSD files by setting Parser.DataHandler.directory.
    """
    def test_custom_directory(self):
        import tempfile
        import os
        import shutil

        handler = Parser.DataHandler(validate=False, escape=False)

        # Create a temporary directory
        tmpdir = tempfile.mkdtemp()
        # Set the custom directory to the temporary directory.
        # This assignment statement will also initialize the local DTD and XSD directories.
        handler.directory = tmpdir

        # Confirm that the two temp directories are named what we want.
        self.assertEqual(handler.local_dtd_dir, os.path.join(handler.directory, 'Bio', 'Entrez', 'DTDs'))
        self.assertEqual(handler.local_xsd_dir, os.path.join(handler.directory, 'Bio', 'Entrez', 'XSDs'))

        # And that they were created.
        self.assertTrue(os.path.isdir(handler.local_dtd_dir))
        self.assertTrue(os.path.isdir(handler.local_xsd_dir))
        shutil.rmtree(tmpdir)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
