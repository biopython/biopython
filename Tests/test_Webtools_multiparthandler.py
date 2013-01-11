# Copyright (c) 2012 Kevin Murray k.d.murray.91@gmail.com
# This test suite is adapted directly from multiparthandler
# multiparthandler is licensed under the LGPL v3

from __future__ import with_statement
import unittest
import tempfile
import os
try:  # Py3 hackery, needs to be tested w/ 2to3
    import urllib2 as a_urllib
except ImportError:
    import urllib.request as a_urllib
import re
from Bio.Webtools.multiparthandler import multiparthandler


class multiparthandler_t(unittest.TestCase):
    validator_url = "http://validator.w3.org/check"
    test_url = ("http://www.w3.org/History/19921103-hypertext/hypertext/WWW/"
        "TheProject.html")
    # Uncomment this, and the tests should fail
    # test_url = "http://www.google.com/"
    opener = a_urllib.build_opener(multiparthandler)
    validator_result_re = r"3 Errors, 4 warning\(s\)"  # expected re pattern

    def setUp(self):
        # Some inter-version portability here.
        # Python 2 returns a str from .read(), py3 returns bytes. if we convert
        # this to bytes, then we can reproducibly decode it to a string.
        try:
            self.test_html = bytes(self.opener.open(self.test_url).read()).\
                decode("UTF-8")
        except NameError:  # If the bytes() fn is not found, i.e. py2.5
            self.test_html = self.opener.open(self.test_url).read()

    def test_post_as_file(self):
        tmp_fd, tmp_fn = tempfile.mkstemp(suffix=".html")
        with open(tmp_fn, "w") as tmp_fh:
            tmp_fh.write(self.test_html)
        params = {
            "ss": "0",  # show source
            "doctype": "Inline",
            "uploaded_file": open(tmp_fn, "r")
            }
        # Portablity hack as above
        try:
            response_html = bytes(self.opener.open(self.validator_url, params).
                    read()).decode("UTF-8")
        except NameError:
            response_html = self.opener.open(self.validator_url, params)\
                    .read()
        os.remove(tmp_fn)
        re_match = re.search(self.validator_result_re, response_html)
        self.assertTrue(re_match is not None)

    def test_post_as_var(self):
        params = {
            "ss": "0",  # show source
            "doctype": "Inline",
            "fragment": self.test_html
            }
        # Portablity hack as above
        try:
            response_html = bytes(self.opener.open(self.validator_url, params).
                    read()).decode("UTF-8")
        except NameError:
            response_html = self.opener.open(self.validator_url, params)\
                    .read()
        re_match = re.search(self.validator_result_re, response_html)
        self.assertTrue(re_match is not None)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
