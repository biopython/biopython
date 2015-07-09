# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO HmmerIO hmmer3-pfamtab parsers."""


import os
import unittest

from Bio import BiopythonExperimentalWarning

import warnings
with warnings.catch_warnings():
   warnings.simplefilter('ignore', BiopythonExperimentalWarning)
   from Bio.SearchIO import parse

# test case files are in the Blast directory
TEST_DIR = 'Hmmer'


def get_file(filename):
    """Returns the path of a test file."""
    return os.path.join(TEST_DIR, filename)


class HmmscanCases(unittest.TestCase):

    fmt = 'hmmscan3-pfamtab'

    def test_pfamtab_31b2_hmmscan_001(self):
        "Test parsing hmmscan3-pfamtab, hmmscan 3.1b2, multiple queries (pfamtab_31b2_hmmscan_001)"

        tab_file = get_file('pfamtab_31b2_hmmscan_001.out')
        qresults = list(parse(tab_file, self.fmt))
        self.assertEqual(5, len(qresults))

        # first qresult
        qresult = qresults[0]
        self.assertEqual(0, len(qresult))
        self.assertEqual('query_0', qresult.id)

        # second qresult
        qresult = qresults[1]
        self.assertEqual(1, len(qresult))
        self.assertEqual('query_1', qresult.id)

        hit = qresult[0]
        self.assertEqual(1, len(hit))
        self.assertEqual(80.5, hit.bitscore)
        self.assertEqual(1e-22, hit.evalue)
        self.assertEqual(1.3, hit.domain_exp_num)
        self.assertEqual(0.3, hit.bias)
        self.assertEqual('query_1', hit.query_id)
        self.assertEqual('<unknown description>', hit.query_description)
        self.assertEqual('Globin', hit.id)
        self.assertEqual('Globin', hit.description)

        hsp = hit[0]
        self.assertEqual(79.8, hsp.bitscore)
        self.assertEqual(1.6e-22, hsp.evalue)
        self.assertEqual(1, hsp.domain_index)
        self.assertEqual(0.3, hsp.bias)
        self.assertEqual(6, hsp.env_start)
        self.assertEqual(113, hsp.env_end)
        self.assertEqual(6, hsp.query_start)
        self.assertEqual(112, hsp.query_end)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(109, hsp.hit_end)
        self.assertEqual('query_1', hit.query_id)
        self.assertEqual('<unknown description>', hit.query_description)
        self.assertEqual('Globin', hsp.hit_id)
        self.assertEqual('Globin', hsp.hit_description)

        # last qresult
        qresult = qresults[-1]
        self.assertEqual(5, len(qresult))
        self.assertEqual('query_4', qresult.id)

        # last hit with multiple hsps
        hit = qresult[2]
        self.assertEqual(2, len(hit))
        self.assertEqual(15.6, hit.bitscore)
        self.assertEqual(0.013, hit.evalue)
        self.assertEqual(2.2, hit.domain_exp_num)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual('query_4', hit.query_id)
        self.assertEqual('<unknown description>', hit.query_description)
        self.assertEqual('HTH_31', hit.id)
        self.assertEqual('Helix-turn-helix domain', hit.description)

        hsp = hit[-1]
        self.assertEqual(0.9, hsp.bitscore)
        self.assertEqual(5.3e+02, hsp.evalue)
        self.assertEqual(2, hsp.domain_index)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(242, hsp.env_start)
        self.assertEqual(270, hsp.env_end)
        self.assertEqual(244, hsp.query_start)
        self.assertEqual(268, hsp.query_end)
        self.assertEqual(38, hsp.hit_start)
        self.assertEqual(62, hsp.hit_end)
        self.assertEqual('query_4', hit.query_id)
        self.assertEqual('<unknown description>', hit.query_description)
        self.assertEqual('HTH_31', hsp.hit_id)
        self.assertEqual('Helix-turn-helix domain', hsp.hit_description)


class HmmersearchCases(unittest.TestCase):

    fmt = 'hmmsearch3-pfamtab'

    def test_pfamtab_31b2_hmmsearch_001(self):
        "Test parsing hmmsearch3-pfamtab, hmmsearch 3.1b2, multiple queries (pfamtab_31b2_hmmsearch_001)"

        tab_file = get_file('pfamtab_31b2_hmmsearch_001.out')
        qresults = list(parse(tab_file, self.fmt))

        self.assertEqual(2, len(qresults))

        # first qresult
        qresult = qresults[0]
        self.assertEqual(0, len(qresult))
        self.assertEqual('query_0', qresult.id)

        # last qresult
        qresult = qresults[-1]
        self.assertEqual(4, len(qresult))
        self.assertEqual('query_1', qresult.id)

        hit = qresult[-1]
        self.assertEqual(2, len(hit))
        self.assertEqual(487.5, hit.bitscore)
        self.assertEqual(2.6e-145, hit.evalue)
        self.assertEqual(2.1, hit.domain_exp_num)
        self.assertEqual(0.0, hit.bias)
        self.assertEqual('query_1', hit.query_id)
        self.assertEqual('<unknown description>', hit.query_description)
        self.assertEqual('sp|P18652|KS6AA_CHICK', hit.id)
        self.assertEqual('Ribosomal protein S6 kinase 2 alpha OS=Gallus '
                         'gallus GN=RPS6KA PE=2 SV=1', hit.description)

        # last qresult, last hit, last HSP (per table ordering)
        hsp = hit[-1]
        self.assertEqual(239.0, hsp.bitscore)
        self.assertEqual(1.6e-69, hsp.evalue)
        self.assertEqual(1, hsp.domain_index)
        self.assertEqual(0.0, hsp.bias)
        self.assertEqual(79, hsp.env_start)
        self.assertEqual(339, hsp.env_end)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(259, hsp.query_end)
        self.assertEqual(79, hsp.hit_start)
        self.assertEqual(338, hsp.hit_end)
        self.assertEqual('query_1', hsp.query_id)
        self.assertEqual('<unknown description>', hsp.query_description)
        self.assertEqual('sp|P18652|KS6AA_CHICK', hsp.hit_id)
        self.assertEqual('Ribosomal protein S6 kinase 2 alpha OS=Gallus '
                         'gallus GN=RPS6KA PE=2 SV=1', hsp.hit_description)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
