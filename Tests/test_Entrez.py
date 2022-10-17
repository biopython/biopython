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
from urllib.parse import urlparse, parse_qs
from urllib.request import Request

from Bio import Entrez
from Bio.Entrez import Parser


# This lets us set the email address to be sent to NCBI Entrez:
Entrez.email = "biopython@biopython.org"
Entrez.api_key = "5cfd4026f9df285d6cfc723c662d74bcbe09"

URL_HEAD = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

# Default values of URL query string (or POST data) when parsed with urllib.parse.parse_qs
QUERY_DEFAULTS = {
    "tool": [Entrez.tool],
    "email": [Entrez.email],
    "api_key": [Entrez.api_key],
}


def get_base_url(parsed):
    """Convert a parsed URL back to string but only include scheme, netloc, and path, omitting query."""
    return parsed.scheme + "://" + parsed.netloc + parsed.path


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


def get_patched_request(patched_urlopen, testcase=None):
    """Get the Request object passed to the patched urlopen() function.

    Expects that the patched function should have been called a single time with a Request instance
    as the only positional argument and no keyword arguments.

    :param patched_urlopen: value returned when entering the context manager created by patch_urlopen.
    :type patched_urlopen: unittest.mock.Mock
    :param testcase: Test case currently being run, which is used to make asserts.
    :type testcase: unittest.TestCase
    :rtype: urllib.urlopen.Request
    """
    args, kwargs = patched_urlopen.call_args

    if testcase is not None:
        testcase.assertEqual(patched_urlopen.call_count, 1)
        testcase.assertEqual(len(args), 1)
        testcase.assertEqual(len(kwargs), 0)
        testcase.assertIsInstance(args[0], Request)

    return args[0]


def deconstruct_request(request, testcase=None):
    """Get the base URL and parsed parameters of a Request object.

    Method may be either GET or POST, POST data should be encoded query params.

    :param request: Request object passed to urlopen().
    :type request: urllib.request.Request
    :param testcase: Test case currently being run, which is used to make asserts.
    :type testcase: unittest.TestCase
    :returns: (base_url, params) tuple.
    """
    parsed = urlparse(request.full_url)

    if request.method == "GET":
        params = parse_qs(parsed.query)

    elif request.method == "POST":
        data = request.data.decode("utf8")
        params = parse_qs(data)

    else:
        raise ValueError(
            "Expected method to be either GET or POST, got %r" % request.method
        )

    return get_base_url(parsed), params


def check_request_ids(testcase, params, expected):
    """Check that the constructed request parameters contain the correct IDs.

    :param testcase: Test case currently being run, which is used to make asserts.
    :type testcase: unittest.TestCase
    :param params: Parsed parameter dictionary returned by `deconstruct_request`.
    :type params: dict
    :param expected: Expected set of IDs, as colleciton of strings.
    """
    testcase.assertEquals(len(params["id"]), 1)
    ids_str = params["id"][0]
    # Compare up to ordering
    testcase.assertCountEqual(ids_str.split(","), expected)


class TestURLConstruction(unittest.TestCase):
    def test_email_warning(self):
        """Test issuing warning when user does not specify email address."""
        email = Entrez.email
        Entrez.email = None

        try:
            with warnings.catch_warnings(record=True) as w:
                Entrez._construct_params(params=None)
                self.assertEqual(len(w), 1)
        finally:
            Entrez.email = email

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

        request = get_patched_request(patched, self)
        self.assertEqual(request.method, "GET")
        base_url, query = deconstruct_request(request, self)

        self.assertEqual(base_url, URL_HEAD + "ecitmatch.cgi")
        query.pop("bdata")  # TODO
        self.assertDictEqual(
            query, {"retmode": ["xml"], "db": [variables["db"]], **QUERY_DEFAULTS}
        )

    def test_construct_cgi_einfo(self):
        """Test constructed url for request to Entrez."""
        with patch_urlopen() as patched:
            Entrez.einfo()

        request = get_patched_request(patched, self)
        self.assertEqual(request.method, "GET")
        base_url, query = deconstruct_request(request, self)

        self.assertEqual(base_url, URL_HEAD + "einfo.fcgi")
        self.assertDictEqual(query, QUERY_DEFAULTS)

    def test_construct_cgi_epost1(self):
        variables = {"db": "nuccore", "id": "186972394,160418"}

        with patch_urlopen() as patched:
            Entrez.epost(**variables)

        request = get_patched_request(patched, self)
        self.assertEqual(request.method, "POST")
        base_url, query = deconstruct_request(request, self)

        self.assertEqual(base_url, URL_HEAD + "epost.fcgi")  # Params in POST data
        self.assertDictEqual(
            query, {"db": [variables["db"]], "id": [variables["id"]], **QUERY_DEFAULTS}
        )

    def test_construct_cgi_epost2(self):
        variables = {"db": "nuccore", "id": ["160418", "160351"]}

        with patch_urlopen() as patched:
            Entrez.epost(**variables)

        request = get_patched_request(patched, self)
        self.assertEqual(request.method, "POST")
        base_url, query = deconstruct_request(request, self)

        self.assertEqual(base_url, URL_HEAD + "epost.fcgi")  # Params in POST data
        check_request_ids(self, query, variables["id"])
        self.assertDictEqual(
            query,
            {
                "db": [variables["db"]],
                "id": query["id"],
                **QUERY_DEFAULTS,
            },
        )

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

        request = get_patched_request(patched, self)
        self.assertEqual(request.method, "GET")
        base_url, query = deconstruct_request(request, self)

        self.assertEqual(base_url, URL_HEAD + "elink.fcgi")
        self.assertDictEqual(
            query,
            {
                "cmd": [variables["cmd"]],
                "db": [variables["db"]],
                "dbfrom": [variables["dbfrom"]],
                "id": [variables["id"]],  # UIDs joined in single string
                **QUERY_DEFAULTS,
            },
        )

    def test_construct_cgi_elink2(self):
        """Commas: Link from protein to gene."""
        variables = {
            "db": "gene",
            "dbfrom": "protein",
            "id": "15718680,157427902,119703751",
        }

        with patch_urlopen() as patched:
            Entrez.elink(**variables)

        request = get_patched_request(patched, self)
        self.assertEqual(request.method, "GET")
        base_url, query = deconstruct_request(request, self)

        self.assertEqual(base_url, URL_HEAD + "elink.fcgi")
        self.assertDictEqual(
            query,
            {
                "db": [variables["db"]],
                "dbfrom": [variables["dbfrom"]],
                "id": [variables["id"]],  # UIDs joined in single string
                **QUERY_DEFAULTS,
            },
        )

    def test_construct_cgi_elink3(self):
        """Multiple ID entries: Find one-to-one links from protein to gene."""
        variables = {
            "db": "gene",
            "dbfrom": "protein",
            "id": ["15718680", "157427902", "119703751"],
        }

        with patch_urlopen() as patched:
            Entrez.elink(**variables)

        request = get_patched_request(patched, self)
        self.assertEqual(request.method, "GET")
        base_url, query = deconstruct_request(request, self)

        self.assertEqual(base_url, URL_HEAD + "elink.fcgi")
        self.assertDictEqual(
            query,
            {
                "db": [variables["db"]],
                "dbfrom": [variables["dbfrom"]],
                "id": query["id"],  # UIDs in multiple separate "id" parameters
                **QUERY_DEFAULTS,
            },
        )

    def test_construct_cgi_efetch(self):
        variables = {
            "db": "protein",
            "id": "15718680,157427902,119703751",
            "retmode": "xml",
        }

        with patch_urlopen() as patched:
            Entrez.efetch(**variables)

        request = get_patched_request(patched, self)
        self.assertEqual(request.method, "GET")
        base_url, query = deconstruct_request(request, self)

        self.assertEqual(base_url, URL_HEAD + "efetch.fcgi")
        self.assertDictEqual(
            query,
            {
                "db": [variables["db"]],
                "id": [variables["id"]],
                "retmode": [variables["retmode"]],
                **QUERY_DEFAULTS,
            },
        )

    def test_default_params(self):
        """Test overriding default values for the "email", "api_key", and "tool" parameters."""
        vars_base = {
            "db": "protein",
            "id": "15718680",
        }
        alt_values = {
            "tool": "mytool",
            "email": "example@example.com",
            "api_key": "test",
        }

        for param in alt_values.keys():
            # Try both an alternate value and None
            for alt_value in [alt_values[param], None]:
                # Try both altering global variable and also passing parameter directly
                for set_global in [False, True]:

                    variables = dict(vars_base)

                    with patch_urlopen() as patched:
                        if set_global:
                            with mock.patch("Bio.Entrez." + param, alt_value):
                                Entrez.efetch(**variables)
                        else:
                            variables[param] = alt_value
                            Entrez.efetch(**variables)

                    request = get_patched_request(patched, self)
                    base_url, query = deconstruct_request(request, self)

                    expected = {k: [v] for k, v in vars_base.items()}
                    expected.update(QUERY_DEFAULTS)
                    if alt_value is None:
                        del expected[param]
                    else:
                        expected[param] = [alt_value]

                    self.assertDictEqual(query, expected)

    def test_has_api_key(self):
        """Test checking whether a Request object specifies an API key.

        The _has_api_key() private function is used to set the delay in _open().
        """
        variables = {
            "db": "protein",
            "id": "15718680",
        }

        for etool in [Entrez.efetch, Entrez.epost]:  # Make both GET and POST requests

            with patch_urlopen() as patched:
                etool(**variables)
            assert Entrez._has_api_key(get_patched_request(patched, self))

            with patch_urlopen() as patched:
                etool(**variables, api_key=None)
            assert not Entrez._has_api_key(get_patched_request(patched, self))

            with patch_urlopen() as patched:
                with mock.patch("Bio.Entrez.api_key", None):
                    etool(**variables)
            assert not Entrez._has_api_key(get_patched_request(patched, self))

    def test_format_ids(self):
        ids = [
            15718680,
            157427902,
            119703751,
            "NP_001098858.1",  # Sequence databases accept accession #s as IDs
        ]
        ids_str = list(map(str, ids))
        ids_formatted = "15718680,157427902,119703751,NP_001098858.1"

        # Single integers or strings should just be converted to string
        for id_ in ids:
            self.assertEqual(Entrez._format_ids(id_), str(id_))

        # List:
        self.assertEqual(Entrez._format_ids(ids), ids_formatted)
        self.assertEqual(Entrez._format_ids(ids_str), ids_formatted)
        # Multiple IDs already joined by commas:
        self.assertEqual(Entrez._format_ids(ids_formatted), ids_formatted)
        # Other iterable types:
        self.assertEqual(Entrez._format_ids(tuple(ids)), ids_formatted)
        self.assertEqual(Entrez._format_ids(tuple(ids_str)), ids_formatted)
        # As set, compare up to reordering
        self.assertCountEqual(Entrez._format_ids(set(ids)).split(","), ids_str)
        self.assertCountEqual(Entrez._format_ids(set(ids_str)).split(","), ids_str)


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
            handler.local_dtd_dir, os.path.join(tmpdir, "Bio", "Entrez", "DTDs")
        )
        self.assertEqual(
            handler.local_xsd_dir, os.path.join(tmpdir, "Bio", "Entrez", "XSDs")
        )

        # And that they were created.
        self.assertTrue(os.path.isdir(handler.local_dtd_dir))
        self.assertTrue(os.path.isdir(handler.local_xsd_dir))
        shutil.rmtree(tmpdir)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
