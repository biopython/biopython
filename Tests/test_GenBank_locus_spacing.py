"""Tests for GenBank LOCUS line whitespace fallback parsing (issue #4846)."""

import io
import unittest
import warnings

from Bio import SeqIO
from Bio import BiopythonParserWarning


class TestLocusSpacing(unittest.TestCase):
    """Tests for tolerance of malformed LOCUS line spacing."""

    def test_snapgene_style_missing_strand_spacing(self):
        """SnapGene-style LOCUS line should parse name with a warning."""
        gb = (
            "LOCUS       pAV-CAG-GFP-corr        2394 bp DNA"
            "                     01-JAN-1980\n"
            "DEFINITION  synthetic circular DNA.\n"
            "ACCESSION   .\n"
            "VERSION     .\n"
            "FEATURES             Location/Qualifiers\n"
            "ORIGIN\n"
            "        1 atgc\n"
            "//\n"
        )
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            record = SeqIO.read(io.StringIO(gb), "gb")
        self.assertEqual(record.name, "pAV-CAG-GFP-corr")
        self.assertTrue(any(issubclass(x.category, BiopythonParserWarning) for x in w))

    def test_correctly_formatted_locus_unaffected(self):
        """Correctly formatted LOCUS lines should not trigger the fallback warning."""
        gb = (
            "LOCUS       pAV-CAG-GFP-corr        2394 bp    DNA"
            "                  01-JAN-1980\n"
            "DEFINITION  synthetic circular DNA.\n"
            "ACCESSION   .\n"
            "VERSION     .\n"
            "FEATURES             Location/Qualifiers\n"
            "ORIGIN\n"
            "        1 atgc\n"
            "//\n"
        )
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            record = SeqIO.read(io.StringIO(gb), "gb")
        self.assertEqual(record.name, "pAV-CAG-GFP-corr")
        strand_warnings = [
            x
            for x in w
            if issubclass(x.category, BiopythonParserWarning)
            and "strand type" in str(x.message)
        ]
        self.assertEqual(len(strand_warnings), 0)


if __name__ == "__main__":
    unittest.main()
