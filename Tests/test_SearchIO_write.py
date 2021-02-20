# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO writing."""

import os
import unittest

from Bio import SearchIO

from search_tests_common import SearchTestBaseClass


class WriteCases(SearchTestBaseClass):
    def tearDown(self):
        if os.path.exists(self.out):
            os.remove(self.out)

    def parse_write_and_compare(
        self, source_file, source_format, out_file, out_format, **kwargs
    ):
        """Compare parsed QueryResults after they have been written to a file."""
        source_qresults = list(SearchIO.parse(source_file, source_format, **kwargs))
        SearchIO.write(source_qresults, out_file, out_format, **kwargs)
        out_qresults = list(SearchIO.parse(out_file, out_format, **kwargs))
        for source, out in zip(source_qresults, out_qresults):
            self.compare_search_obj(source, out)

    def read_write_and_compare(
        self, source_file, source_format, out_file, out_format, **kwargs
    ):
        """Compare read QueryResults after it has been written to a file."""
        source_qresult = SearchIO.read(source_file, source_format, **kwargs)
        SearchIO.write(source_qresult, out_file, out_format, **kwargs)
        out_qresult = SearchIO.read(out_file, out_format, **kwargs)
        self.compare_search_obj(source_qresult, out_qresult)


class BlastXmlWriteCases(WriteCases):

    fmt = "blast-xml"
    out = os.path.join("Blast", "test_write.xml")

    def test_write_single_from_blastxml(self):
        """Test blast-xml writing from blast-xml, BLAST 2.2.26+, single query (xml_2226_blastp_004.xml)."""
        source = os.path.join("Blast", "xml_2226_blastp_004.xml")
        self.parse_write_and_compare(source, self.fmt, self.out, self.fmt)
        self.read_write_and_compare(source, self.fmt, self.out, self.fmt)

    def test_write_multiple_from_blastxml(self):
        """Test blast-xml writing from blast-xml, BLAST 2.2.26+, multiple queries (xml_2226_blastp_001.xml)."""
        source = os.path.join("Blast", "xml_2226_blastp_001.xml")
        self.parse_write_and_compare(source, self.fmt, self.out, self.fmt)


class BlastTabWriteCases(WriteCases):

    fmt = "blast-tab"
    out = os.path.join("Blast", "test_write.txt")

    def test_write_single_from_blasttab(self):
        """Test blast-tab writing from blast-tab, BLAST 2.2.26+, single query (tab_2226_tblastn_004.txt)."""
        source = os.path.join("Blast", "tab_2226_tblastn_004.txt")
        self.parse_write_and_compare(source, self.fmt, self.out, self.fmt)
        self.read_write_and_compare(source, self.fmt, self.out, self.fmt)

    def test_write_multiple_from_blasttab(self):
        """Test blast-tab writing from blast-tab, BLAST 2.2.26+, multiple queries (tab_2226_tblastn_001.txt)."""
        source = os.path.join("Blast", "tab_2226_tblastn_001.txt")
        self.parse_write_and_compare(source, self.fmt, self.out, self.fmt)

    def test_write_single_from_blasttabc(self):
        """Test blast-tabc writing from blast-tabc, BLAST 2.2.26+, single query (tab_2226_tblastn_008.txt)."""
        source = os.path.join("Blast", "tab_2226_tblastn_008.txt")
        self.parse_write_and_compare(
            source, self.fmt, self.out, self.fmt, comments=True
        )
        self.read_write_and_compare(source, self.fmt, self.out, self.fmt, comments=True)

    def test_write_multiple_from_blasttabc(self):
        """Test blast-tabc writing from blast-tabc, BLAST 2.2.26+, multiple queries (tab_2226_tblastn_005.txt)."""
        source = os.path.join("Blast", "tab_2226_tblastn_005.txt")
        self.parse_write_and_compare(
            source, self.fmt, self.out, self.fmt, comments=True
        )

    def test_write_multiple_from_blasttabc_allfields(self):
        """Test blast-tabc writing from blast-tabc, BLAST 2.2.28+, multiple queries (tab_2228_tblastx_001.txt)."""
        source = os.path.join("Blast", "tab_2228_tblastx_001.txt")
        fields = [
            "qseqid",
            "qgi",
            "qacc",
            "qaccver",
            "qlen",
            "sseqid",
            "sallseqid",
            "sgi",
            "sallgi",
            "sacc",
            "saccver",
            "sallacc",
            "slen",
            "qstart",
            "qend",
            "sstart",
            "send",
            "qseq",
            "sseq",
            "evalue",
            "bitscore",
            "score",
            "length",
            "pident",
            "nident",
            "mismatch",
            "positive",
            "gapopen",
            "gaps",
            "ppos",
            "frames",
            "qframe",
            "sframe",
            "btop",
            "staxids",
            "sscinames",
            "scomnames",
            "sblastnames",
            "sskingdoms",
            "stitle",
            "salltitles",
            "sstrand",
            "qcovs",
            "qcovhsp",
        ]
        self.parse_write_and_compare(
            source, self.fmt, self.out, self.fmt, comments=True, fields=fields
        )


class HmmerTabWriteCases(WriteCases):

    fmt = "hmmer3-tab"
    out = os.path.join("Hmmer", "test_write.txt")

    def test_write_single_from_hmmertab(self):
        """Test hmmer3-tab writing from hmmer3-tab, HMMER 3.0, single query (tab_30_hmmscan_004.out)."""
        source = os.path.join("Hmmer", "tab_30_hmmscan_004.out")
        self.parse_write_and_compare(source, self.fmt, self.out, self.fmt)
        self.read_write_and_compare(source, self.fmt, self.out, self.fmt)

    def test_write_multiple_from_hmmertab(self):
        """Test hmmer3-tab writing from hmmer3-tab, HMMER 3.0, multiple queries (tab_30_hmmscan_001.out)."""
        source = os.path.join("Hmmer", "tab_30_hmmscan_001.out")
        self.parse_write_and_compare(source, self.fmt, self.out, self.fmt)


class HmmerDomtabWriteCases(WriteCases):

    out = os.path.join("Hmmer", "test_write.txt")

    def test_write_single_from_hmmscandomtab(self):
        """Test hmmscan-domtab writing from hmmscan-domtab, HMMER 3.0, single query (tab_30_hmmscan_004.out)."""
        source = os.path.join("Hmmer", "domtab_30_hmmscan_004.out")
        fmt = "hmmscan3-domtab"
        self.parse_write_and_compare(source, fmt, self.out, fmt)
        self.read_write_and_compare(source, fmt, self.out, fmt)

    def test_write_multiple_from_hmmscandomtab(self):
        """Test hmmscan-domtab writing from hmmscan-domtab, HMMER 3.0, multiple queries (tab_30_hmmscan_001.out)."""
        source = os.path.join("Hmmer", "domtab_30_hmmscan_001.out")
        fmt = "hmmscan3-domtab"
        self.parse_write_and_compare(source, fmt, self.out, fmt)

    def test_write_single_from_hmmsearchdomtab(self):
        """Test hmmsearch-domtab writing from hmmsearch-domtab, HMMER 3.0, single query (tab_30_hmmscan_004.out)."""
        source = os.path.join("Hmmer", "domtab_30_hmmsearch_001.out")
        fmt = "hmmsearch3-domtab"
        self.parse_write_and_compare(source, fmt, self.out, fmt)
        self.read_write_and_compare(source, fmt, self.out, fmt)


class BlatPslWriteCases(WriteCases):

    fmt = "blat-psl"
    out = os.path.join("Blat", "test_write.txt")

    def test_write_single_from_blatpsl(self):
        """Test blat-psl writing from blat-psl, single query (psl_34_004.psl)."""
        source = os.path.join("Blat", "psl_34_004.psl")
        self.parse_write_and_compare(source, self.fmt, self.out, self.fmt)
        self.read_write_and_compare(source, self.fmt, self.out, self.fmt)

    def test_write_single_from_blatpsl_protein_query(self):
        """Test blat-psl writing from blat-psl, single query (psl_35_002.psl)."""
        source = os.path.join("Blat", "psl_35_002.psl")
        self.parse_write_and_compare(source, self.fmt, self.out, self.fmt)
        self.read_write_and_compare(source, self.fmt, self.out, self.fmt)

    def test_write_multiple_from_blatpsl(self):
        """Test blat-psl writing from blat-psl, multiple queries (psl_34_001.psl)."""
        source = os.path.join("Blat", "psl_34_001.psl")
        self.parse_write_and_compare(source, self.fmt, self.out, self.fmt)

    def test_write_single_from_blatpslx(self):
        """Test blat-pslx writing from blat-pslx, single query (pslx_34_004.pslx)."""
        source = os.path.join("Blat", "pslx_34_004.pslx")
        self.parse_write_and_compare(source, self.fmt, self.out, self.fmt, pslx=True)
        self.read_write_and_compare(source, self.fmt, self.out, self.fmt, pslx=True)

    def test_write_single_from_blatpslx_protein_query(self):
        """Test blat-pslx writing from blat-pslx, single query (pslx_35_002.pslx)."""
        source = os.path.join("Blat", "pslx_35_002.pslx")
        self.parse_write_and_compare(source, self.fmt, self.out, self.fmt, pslx=True)
        self.read_write_and_compare(source, self.fmt, self.out, self.fmt, pslx=True)

    def test_write_multiple_from_blatpslx(self):
        """Test blat-pslx writing from blat-pslx, multiple queries (pslx_34_001.pslx)."""
        source = os.path.join("Blat", "pslx_34_001.pslx")
        self.parse_write_and_compare(source, self.fmt, self.out, self.fmt, pslx=True)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
