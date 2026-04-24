"""Tests for GenBank LOCUS line whitespace fallback parsing (issue #4846)."""

import io
import unittest
import warnings

from Bio import SeqIO
from Bio import BiopythonParserWarning


class TestLocusSpacing(unittest.TestCase):

    def test_snapgene_style_missing_strand_spacing(self):
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
        self.assertEqual(str(record.annotations.get("molecule_type")), "DNA")
        self.assertEqual(record.annotations.get("date"), "01-JAN-1980")
        self.assertTrue(any(issubclass(x.category, BiopythonParserWarning) for x in w))

    def test_missing_strand_with_topology_and_division(self):
        gb = (
            "LOCUS       myplasmid               1000 bp DNA"
            "       circular  BCT 15-MAR-1991\n"
            "DEFINITION  test plasmid.\n"
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
        self.assertEqual(record.name, "myplasmid")
        self.assertEqual(str(record.annotations.get("molecule_type")), "DNA")
        self.assertEqual(record.annotations.get("topology"), "circular")
        self.assertEqual(record.annotations.get("data_file_division"), "BCT")
        self.assertEqual(record.annotations.get("date"), "15-MAR-1991")

    def test_correctly_formatted_locus_unaffected(self):
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
