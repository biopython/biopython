# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO indexing.

Since the SearchIO indexing core shares a lot of similarity with SeqIO
indexing, here we are only testing the format-specific parsers and not
the SearchIO indexing core itself.

"""

import unittest
import warnings

from Bio import SearchIO

from search_tests_common import compare_qresult


class SearchIndexCases(unittest.TestCase):

    def check_index(self, filename, format):
        parsed = list(SearchIO.parse(filename, format))
        indexed = SearchIO.index(filename, format)
        db_indexed = SearchIO.index_db(':memory:', [filename], format)

        # check length of parsed and indexed
        self.assertEqual(len(parsed), len(indexed.keys()))
        self.assertEqual(len(parsed), len(db_indexed.keys()))

        for qres in parsed:
            idx_qres = indexed[qres.id]
            dbidx_qres = db_indexed[qres.id]
            # parsed and indexed qresult are different objects!
            self.assertNotEqual(id(qres), id(idx_qres))
            self.assertNotEqual(id(qres), id(dbidx_qres))
            # but they should have the same attribute values
            self.assertTrue(compare_qresult(qres, idx_qres, format))
            self.assertTrue(compare_qresult(qres, dbidx_qres, format))


class BlastXmlIndexCases(SearchIndexCases):

    fmt = 'blast-xml'

    def test_blastxml_blastp_2212(self):
        filename = 'Blast/xbt001.xml'
        self.check_index(filename, self.fmt)

    def test_blastxml_blastp_2218p(self):
        filename = 'Blast/xbt006.xml'
        self.check_index(filename, self.fmt)

    def test_blastxml_blastx_2222p(self):
        filename = 'Blast/xbt009.xml'
        self.check_index(filename, self.fmt)

    def test_blastxml_tblastn_2224p_mult(self):
        filename = 'Blast/xbt026.xml'
        self.check_index(filename, self.fmt)

    def test_blastxml_tblastn_2224p_none(self):
        filename = 'Blast/xbt027.xml'
        self.check_index(filename, self.fmt)

    def test_blastxml_tblastn_2224p_sing(self):
        filename = 'Blast/xbt029.xml'
        self.check_index(filename, self.fmt)


class BlastTabIndexCases(SearchIndexCases):

    fmt = 'blast-tab'

    def test_blasttab_mult(self):
        filename = 'Blast/tbt001.txt'
        self.check_index(filename, self.fmt)

    def test_blasttab_none(self):
        filename = 'Blast/tbt002.txt'
        self.check_index(filename, self.fmt)

    def test_blasttab_sing(self):
        filename = 'Blast/tbt004.txt'
        self.check_index(filename, self.fmt)

    def test_blasttab_comment_mult(self):
        filename = 'Blast/tbt005.txt'
        self.check_index(filename, self.fmt)

    def test_blasttab_comment_none(self):
        filename = 'Blast/tbt006.txt'
        self.check_index(filename, self.fmt)

    def test_blasttab_comment_sing(self):
        filename = 'Blast/tbt008.txt'
        self.check_index(filename, self.fmt)

    def test_blasttab_custom_cols(self):
        filename = 'Blast/tbt009.txt'
        self.assertRaises(AssertionError, self.check_index, filename, self.fmt)

    def test_blasttab_comment_custom_cols(self):
        filename = 'Blast/tbt010.txt'
        self.check_index(filename, self.fmt)

    def test_blasttab_comment_all_cols(self):
        filename = 'Blast/tbt011.txt'
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            self.check_index(filename, self.fmt)
            self.assertTrue(issubclass(w[-1].category, UserWarning))


class HmmerTextIndexCases(SearchIndexCases):

    fmt = 'hmmer-text'

    def test_hmmertext_hmmscan_mult(self):
        filename = 'Hmmer/text_hmm001.out'
        self.check_index(filename, self.fmt)

    def test_hmmertext_hmmscan_none(self):
        filename = 'Hmmer/text_hmm002.out'
        self.check_index(filename, self.fmt)

    def test_hmmertext_hmmscan_sing(self):
        filename = 'Hmmer/text_hmm006.out'
        self.check_index(filename, self.fmt)

    def test_hmmertext_hmmscan_sing_noali(self):
        filename = 'Hmmer/text_hmm007.out'
        self.check_index(filename, self.fmt)

    def test_hmmertext_hmmscan_sing_notextw(self):
        filename = 'Hmmer/text_hmm008.out'
        self.check_index(filename, self.fmt)

    def test_hmmertext_hmmsearch_mult(self):
        filename = 'Hmmer/text_hmm013.out'
        self.check_index(filename, self.fmt)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
