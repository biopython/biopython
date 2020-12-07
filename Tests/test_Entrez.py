# Copyright 2015 by Carlos Pena.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Offline tests for two Entrez features.

(1) the URL construction of NCBI's Entrez services.
(2) setting a custom directory for DTD and XSD downloads.
"""

import unittest
from unittest import mock
import warnings
from http.client import HTTPMessage

from Bio import Entrez
from Bio.Entrez import Parser


# This lets us set the email address to be sent to NCBI Entrez:
Entrez.email = "biopython@biopython.org"
Entrez.api_key = "5cfd4026f9df285d6cfc723c662d74bcbe09"

URL_HEAD = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
URL_TOOL = "tool=biopython"
URL_EMAIL = "email=biopython%40biopython.org"
URL_API_KEY = "api_key=5cfd4026f9df285d6cfc723c662d74bcbe09"


def mock_httpresponse(code=200, content_type="/xml"):
    """Create a mocked version of a response object returned by urlopen().

    :param int code: Value of "code" attribute.
    :param str content_type: Used to set the "Content-Type" header in the "headers" attribute. This
        is checked in Entrez._open() to determine if the response data is plain text.
    """
    resp = mock.NonCallableMock()
    resp.code = code

    resp.headers = HTTPMessage()
    resp.headers.add_header("Content-Type", content_type + "; charset=UTF-8")

    return resp


def patch_urlopen(**kwargs):
    """Create a context manager which replaces Bio.Entrez.urlopen with a mocked version.

    Within the decorated function, Bio.Entrez.urlopen will be replaced with a unittest.mock.Mock
    object which when called simply records the arguments passed to it and returns a mocked response
    object. The actual urlopen function will not be called so no request will actually be made.
    """
    response = mock_httpresponse(**kwargs)
    return unittest.mock.patch("Bio.Entrez.urlopen", return_value=response)


def get_patched_get_url(patched_urlopen, testcase=None):
    """Get the URL of the GET request made to the patched urlopen() function.

    Expects that the patched function should have been called a single time with the url as the only
    positional argument and no keyword arguments.

    :param patched_urlopen: value returned when entering the context manager created by patch_urlopen.
    :type patched_urlopen: unittest.mock.Mock
    :param testcase: Test case currently being run, which is used to make asserts
    :type testcase: unittest.TestCase
    """
    args, kwargs = patched_urlopen.call_args

    if testcase is not None:
        testcase.assertEqual(patched_urlopen.call_count, 1)
        testcase.assertEqual(len(args), 1)
        testcase.assertEqual(len(kwargs), 0)

    return args[0]


def get_patched_post_args(patched_urlopen, testcase=None, decode=False):
    """Get the URL and content data of the POST request made to the patched urlopen() function.

    Expects that the patched function should have been called a single time with the url as the only
    positional argument and "data" as the only keyword argument. Returns a (url, data) tuple.

    :param patched_urlopen: value returned when entering the context manager created by patch_urlopen.
    :type patched_urlopen: unittest.mock.Mock
    :param testcase: Test case currently being run, which is used to make asserts
    :type testcase: unittest.TestCase
    :param bool decode: Decode the value of the "data" keyword argument before returning
    """
    args, kwargs = patched_urlopen.call_args

    if testcase is not None:
        testcase.assertEqual(patched_urlopen.call_count, 1)
        testcase.assertEqual(len(args), 1)
        testcase.assertEqual(list(kwargs), ["data"])

    data = kwargs["data"]
    if decode:
        data = data.decode("utf8")

    return args[0], data


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
            "year": "1991",
            "volume": "88",
            "first_page": "3248",
            "author_name": "mann bj",
            "key": "citation_1",
        }
        variables = Entrez._update_ecitmatch_variables(
            {"db": "pubmed", "bdata": [citation]}
        )

        with patch_urlopen() as patched:
            Entrez.ecitmatch(**variables)

        result_url = get_patched_get_url(patched, self)

        self.assertIn("retmode=xml", result_url)
        self.assertIn(URL_API_KEY, result_url)

    def test_construct_cgi_einfo(self):
        """Test constructed url for request to Entrez."""
        with patch_urlopen() as patched:
            Entrez.einfo()

        result_url = get_patched_get_url(patched, self)

        self.assertTrue(result_url.startswith(URL_HEAD + "einfo.fcgi?"), result_url)
        self.assertIn(URL_TOOL, result_url)
        self.assertIn(URL_EMAIL, result_url)

    def test_construct_cgi_epost1(self):
        variables = {"db": "nuccore", "id": "186972394,160418"}

        with patch_urlopen() as patched:
            Entrez.epost(**variables)

        result_url, options = get_patched_post_args(patched, self, decode=True)

        self.assertEqual(URL_HEAD + "epost.fcgi", result_url)

    def test_construct_cgi_epost2(self):
        variables = {"db": "nuccore", "id": ["160418", "160351"]}

        with patch_urlopen() as patched:
            Entrez.epost(**variables)

        result_url, options = get_patched_post_args(patched, self, decode=True)

        self.assertEqual(URL_HEAD + "epost.fcgi", result_url)

    def test_construct_cgi_elink1(self):
        variables = {
            "cmd": "neighbor_history",
            "db": "nucleotide",
            "dbfrom": "protein",
            "id": "22347800,48526535",
            "query_key": None,
            "webenv": None,
        }

        with patch_urlopen() as patched:
            Entrez.elink(**variables)

        result_url = get_patched_get_url(patched, self)

        self.assertTrue(result_url.startswith(URL_HEAD + "elink.fcgi?"), result_url)
        self.assertIn(URL_TOOL, result_url)
        self.assertIn(URL_EMAIL, result_url)
        self.assertIn("id=22347800%2C48526535", result_url)
        self.assertIn(URL_API_KEY, result_url)

    def test_construct_cgi_elink2(self):
        """Commas: Link from protein to gene."""
        variables = {
            "db": "gene",
            "dbfrom": "protein",
            "id": "15718680,157427902,119703751",
        }

        with patch_urlopen() as patched:
            Entrez.elink(**variables)

        result_url = get_patched_get_url(patched, self)

        self.assertTrue(result_url.startswith(URL_HEAD + "elink.fcgi"), result_url)
        self.assertIn(URL_TOOL, result_url)
        self.assertIn(URL_EMAIL, result_url)
        self.assertIn("id=15718680%2C157427902%2C119703751", result_url, result_url)
        self.assertIn(URL_API_KEY, result_url)

    def test_construct_cgi_elink3(self):
        """Multiple ID entries: Find one-to-one links from protein to gene."""
        variables = {
            "db": "gene",
            "dbfrom": "protein",
            "id": ["15718680", "157427902", "119703751"],
        }

        with patch_urlopen() as patched:
            Entrez.elink(**variables)

        result_url = get_patched_get_url(patched, self)

        self.assertTrue(result_url.startswith(URL_HEAD + "elink.fcgi"), result_url)
        self.assertIn(URL_TOOL, result_url)
        self.assertIn(URL_EMAIL, result_url)
        self.assertIn("id=15718680", result_url)
        self.assertIn("id=157427902", result_url)
        self.assertIn("id=119703751", result_url)
        self.assertIn(URL_API_KEY, result_url)

    def test_construct_cgi_efetch(self):
        variables = {
            "db": "protein",
            "id": "15718680,157427902,119703751",
            "retmode": "xml",
        }

        with patch_urlopen() as patched:
            Entrez.efetch(**variables)

        result_url = get_patched_get_url(patched, self)

        self.assertTrue(result_url.startswith(URL_HEAD + "efetch.fcgi?"), result_url)
        self.assertIn(URL_TOOL, result_url)
        self.assertIn(URL_EMAIL, result_url)
        self.assertIn("id=15718680%2C157427902%2C119703751", result_url, result_url)
        self.assertIn(URL_API_KEY, result_url)


class CustomDirectoryTest(unittest.TestCase):
    """Offline unit test for custom directory feature.

    Allow user to specify a custom directory for Entrez DTD/XSD files by setting
    Parser.DataHandler.directory.
    """

    def test_custom_directory(self):
        import tempfile
        import os
        import shutil

        handler = Parser.DataHandler(validate=False, escape=False)

        # Create a temporary directory
        tmpdir = tempfile.mkdtemp()
        # Set the custom directory to the temporary directory.
        # This assignment statement will also initialize the local DTD and XSD
        # directories.
        Parser.DataHandler.directory = tmpdir

        # Confirm that the two temp directories are named what we want.
        self.assertEqual(
            handler.local_dtd_dir, os.path.join(tmpdir, "Bio", "Entrez", "DTDs"),
        )
        self.assertEqual(
            handler.local_xsd_dir, os.path.join(tmpdir, "Bio", "Entrez", "XSDs"),
        )

        # And that they were created.
        self.assertTrue(os.path.isdir(handler.local_dtd_dir))
        self.assertTrue(os.path.isdir(handler.local_xsd_dir))
        shutil.rmtree(tmpdir)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
