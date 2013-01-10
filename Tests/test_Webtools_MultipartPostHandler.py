import unittest
import tempfile
import os
import urllib2
import re
from Bio.Webtools.multiparthandler import multiparthandler

# Make sure test is skipped if not internet connected
import requires_internet
requires_internet.check()


class MultipartPost(unittest.TestCase):
    
    validator_url = "http://validator.w3.org/check"
    test_url = ("http://www.w3.org/History/19921103-hypertext/hypertext/WWW/"
        "TheProject.html")
    # Uncomment this, and the tests should fail
    # test_url = "http://www.google.com/"
    opener = urllib2.build_opener(multiparthandler)
    validator_result_re = "3 Errors, 4 warning\(s\)"  # expected re pattern

    def setUp(self):
        self.test_html = self.opener.open(self.test_url).read()

    def test_post_as_file(self):
        temp = tempfile.mkstemp(suffix=".html")
        os.write(temp[0], self.test_html )
        params = {
            "ss" : "0", # show source
            "doctype" : "Inline",
            "uploaded_file" : open(temp[1], "rb")
            }
        response_html = self.opener.open(self.validator_url, params).read()
        os.remove(temp[1])
        re_match = re.search(self.validator_result_re, response_html)
        self.assertTrue(re_match is not None)

    def test_post_as_var(self):
        params = {
            "ss" : "0", # show source
            "doctype" : "Inline",
            "fragment" : self.test_html
            }
        response_html = self.opener.open(self.validator_url, params).read()
        re_match = re.search(self.validator_result_re, response_html)
        self.assertTrue(re_match is not None)



if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
