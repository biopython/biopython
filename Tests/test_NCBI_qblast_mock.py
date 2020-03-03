# Copyright 2020 by Markus Piotrowski.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Testing online code for fetching NCBI qblast.

This test file imports the TestCases from ``NCBI_qblast_testCases.py``.
If you want to add more tests, do this over there.

Here, these tests are running *offline* with a mocked version of
``urllib.request.urlopen`` and ``time``.

"""
import unittest
from unittest import mock
from io import BytesIO


# Import all test classes deliberately:
from NCBI_qblast_testCases import *  # noqa: F401, F403
import NCBI_qblast_testCases


def mock_response(call_count):
    """Mimick an NCBI qblast response."""
    # Each use of NCBIWWW.qblast makes two urlopen calls with different responses:
    # a. the 'wait' page, and b. the result.
    if call_count % 2 == 0:
        # This mimicks the 'wait' page:
        return BytesIO(open("Blast/mock_wait.html", "rb").read())
    else:
        # These are the results:
        if call_count == 1:
            # test_blastp_nr_actin
            return BytesIO(open("Blast/mock_actin.xml", "rb").read())
        elif call_count == 3:
            # test_discomegablast
            return BytesIO(open("Blast/mock_disco.xml", "rb").read())
        elif call_count == 5:
            # test_orchid_est
            return BytesIO(open("Blast/mock_orchid.xml", "rb").read())
        elif call_count == 7:
            # test_pcr_primers
            return BytesIO(open("Blast/mock_pcr.xml", "rb").read())
        elif call_count == 9:
            # test_short_query #1: no results
            return BytesIO(open("Blast/mock_short_empty.xml", "rb").read())
        else:  # call_counts 11 & 13
            # test_short_query #2 & #3: 5 results
            return BytesIO(open("Blast/mock_short_result.xml", "rb").read())


NCBI_qblast_testCases.NCBIWWW.time.sleep = mock.Mock()

NCBI_qblast_testCases.NCBIWWW.urlopen = mock.Mock(
    side_effect=map(mock_response, range(14))
)

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
