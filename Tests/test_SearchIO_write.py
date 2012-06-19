# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO writing."""

import os
import unittest

from Bio import SearchIO

from search_tests_common import compare_qresult


class WriteCases(unittest.TestCase):

    def tearDown(self):
        if os.path.exists(self.out):
            os.remove(self.out)

    def write_and_compare(self, source_file, source_format, out_file, out_format):
        source_qresults = list(SearchIO.parse(source_file, source_format))
        SearchIO.write(source_qresults, out_file, out_format)
        out_qresults = list(SearchIO.parse(out_file, out_format))
        for source, out in zip(source_qresults, out_qresults):
            self.assertTrue(compare_qresult(source, out, out_format))

class BlastXmlWriteCases(WriteCases):

    fmt = 'blast-xml'
    out = os.path.join('Blast', 'test_write.xml')

    def test_write_single_from_blastxml(self):
        """Test blast-xml writing from blast-xml, BLAST 2.2.26+, single query (xml_2226_blastp_004.xml)"""
        source = os.path.join('Blast', 'xml_2226_blastp_004.xml')
        self.write_and_compare(source, self.fmt, self.out, self.fmt)

    def test_write_multiple_from_blastxml(self):
        """Test blast-xml writing from blast-xml, BLAST 2.2.26+, multiple queries (xml_2226_blastp_001.xml)"""
        source = os.path.join('Blast', 'xml_2226_blastp_001.xml')
        self.write_and_compare(source, self.fmt, self.out, self.fmt)


class BlastTablWriteCases(WriteCases):

    fmt = 'blast-tab'
    out = os.path.join('Blast', 'test_write.txt')

    def test_write_single_from_blasttab(self):
        """Test blast-tab writing from blast-xml, BLAST 2.2.26+, single query (tab_2226_tblastn_004.txt)"""
        source = os.path.join('Blast', 'tab_2226_tblastn_004.txt')
        self.write_and_compare(source, self.fmt, self.out, self.fmt)

    def test_write_multiple_from_blasttab(self):
        """Test blast-tab writing from blast-xml, BLAST 2.2.26+, single query (tab_2226_tblastn_001.txt)"""
        source = os.path.join('Blast', 'tab_2226_tblastn_001.txt')
        self.write_and_compare(source, self.fmt, self.out, self.fmt)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
