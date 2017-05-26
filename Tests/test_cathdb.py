"""Testing Bio.cathdb online code."""

import time
import unittest

import requires_internet
requires_internet.check()

# We want to test these:

from Bio.cathdb import check_progress
from Bio.cathdb import retrieve_results
from Bio.cathdb import search_by_sequence



# TODO - Use with statement when drop Python 2

class cathdbOnlineTests(unittest.TestCase):
    """Test cathdb online resources."""

    def setUp(self):
        self.fasta = '''>tr|G4VGF5|G4VGF5_SCHMA
        MCSHYAQRNNFSCGGYGFIDFVSEDAANEALQQIKETHPSFTIKFAKENEKDKTNLYVTN'''

    def test_search(self):
        task_id = search_by_sequence(self.fasta)
        self.assertTrue(task_id)

    def test_check_progress(self):
        task_id = search_by_sequence(self.fasta)
        progress_response = check_progress(task_id)
        self.assertTrue(progress_response['message'])

    def test_retrieve_results(self):
        task_id = search_by_sequence(self.fasta)
        progress_response = check_progress(task_id)
        while progress_response['success'] != 1:
            progress_response = check_progress(task_id)
            time.sleep(2)
        results = retrieve_results(task_id)
        self.assertEqual(results['query_fasta'].strip(), self.fasta.strip().replace(' ',''))

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
