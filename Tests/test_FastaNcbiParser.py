# Copyright 2020 by Tianyi Shi.  All rights reserved.
# Copyright 2026 by Martin Mokrejš.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Tests for NCBI FASTA header parsing in Bio.SeqIO.FastaIO."""

import os
import unittest
from io import StringIO

from Bio import SeqIO
from Bio.SeqIO.FastaIO import _parse_ncbi_fasta_header, FastaNcbiIterator

# Biopython tests expect to run from the Tests/ directory
if os.path.isdir("Tests"):
    os.chdir("Tests")


class TestParseNcbiFastaHeader(unittest.TestCase):
    """Tests for _parse_ncbi_fasta_header function."""

    def test_gi_genbank(self):
        """Parse combined gi and GenBank header."""
        id_, name, dbxrefs = _parse_ncbi_fasta_header(
            "gi|186972394|gb|EU490707.1| Selenipedium aequinoctiale matK"
        )
        self.assertEqual(id_, "gi|186972394|gb|EU490707.1|")
        self.assertEqual(name, "186972394")
        self.assertIn("GI:186972394", dbxrefs)
        self.assertIn("GenBank:EU490707.1", dbxrefs)

    def test_gi_ddbj(self):
        """Parse combined gi and DDBJ header."""
        id_, name, dbxrefs = _parse_ncbi_fasta_header(
            "gi|3298468|dbj|BAA31520.1| SAMIPF"
        )
        self.assertEqual(id_, "gi|3298468|dbj|BAA31520.1|")
        self.assertEqual(name, "3298468")
        self.assertIn("GI:3298468", dbxrefs)
        self.assertIn("DDBJ:BAA31520.1", dbxrefs)

    def test_swissprot(self):
        """Parse Swiss-Prot header."""
        id_, name, dbxrefs = _parse_ncbi_fasta_header(
            "sp|P05698|LYSC_HUMAN Lysozyme C OS=Homo sapiens"
        )
        self.assertEqual(id_, "sp|P05698|LYSC_HUMAN")
        self.assertEqual(name, "P05698")
        self.assertEqual(dbxrefs, ["UniProtKB/Swiss-Prot:P05698"])

    def test_trembl(self):
        """Parse TrEMBL header."""
        id_, name, dbxrefs = _parse_ncbi_fasta_header(
            "tr|A0A0C5B5G6|A0A0C5B5G6_9MICO Uncharacterized protein"
        )
        self.assertEqual(id_, "tr|A0A0C5B5G6|A0A0C5B5G6_9MICO")
        self.assertEqual(name, "A0A0C5B5G6")
        self.assertEqual(dbxrefs, ["UniProtKB/TrEMBL:A0A0C5B5G6"])

    def test_pdb(self):
        """Parse PDB header with chain."""
        id_, name, dbxrefs = _parse_ncbi_fasta_header("pdb|1A2B|C some PDB chain")
        self.assertEqual(id_, "pdb|1A2B|C")
        self.assertEqual(name, "1A2B")
        self.assertEqual(dbxrefs, ["PDB:1A2B"])

    def test_refseq(self):
        """Parse RefSeq header."""
        id_, name, dbxrefs = _parse_ncbi_fasta_header(
            "ref|NM_001301717.2| Homo sapiens TP53 mRNA"
        )
        self.assertEqual(id_, "ref|NM_001301717.2|")
        self.assertEqual(name, "NM_001301717.2")
        self.assertEqual(dbxrefs, ["RefSeq:NM_001301717.2"])

    def test_embl(self):
        """Parse EMBL header."""
        id_, name, dbxrefs = _parse_ncbi_fasta_header(
            "emb|CAA12345.6| some EMBL record"
        )
        self.assertEqual(id_, "emb|CAA12345.6|")
        self.assertEqual(name, "CAA12345.6")
        self.assertEqual(dbxrefs, ["EMBL:CAA12345.6"])

    def test_patent(self):
        """Parse patent header."""
        id_, name, dbxrefs = _parse_ncbi_fasta_header(
            "pat|US|RE33188|1 a patent sequence"
        )
        self.assertEqual(id_, "pat|US|RE33188|1")
        self.assertEqual(name, "US")
        self.assertEqual(dbxrefs, ["Patent:US"])

    def test_gnl(self):
        """Parse general database reference header."""
        id_, name, dbxrefs = _parse_ncbi_fasta_header(
            "gnl|taxon|9606 general database reference"
        )
        self.assertEqual(id_, "gnl|taxon|9606")
        self.assertEqual(name, "taxon")
        self.assertEqual(dbxrefs, ["General database reference:taxon"])

    def test_plain_header(self):
        """Parse plain header without pipes."""
        id_, name, dbxrefs = _parse_ncbi_fasta_header(
            "plain_id no pipe delimiters here"
        )
        self.assertEqual(id_, "plain_id")
        self.assertEqual(name, "plain_id")
        self.assertEqual(dbxrefs, [])

    def test_empty_title(self):
        """Parse empty title."""
        id_, name, dbxrefs = _parse_ncbi_fasta_header("")
        self.assertEqual(id_, "")
        self.assertEqual(name, "")
        self.assertEqual(dbxrefs, [])

    def test_unknown_prefix(self):
        """Unknown prefixes are silently skipped."""
        id_, name, dbxrefs = _parse_ncbi_fasta_header(
            "xyz|12345| some unknown database"
        )
        self.assertEqual(id_, "xyz|12345|")
        # No recognized prefix, falls back to first_word
        self.assertEqual(name, "xyz|12345|")
        self.assertEqual(dbxrefs, [])

    def test_trailing_pipe(self):
        """Trailing pipes don't cause crashes."""
        id_, name, dbxrefs = _parse_ncbi_fasta_header("gi|12345| some description")
        self.assertEqual(id_, "gi|12345|")
        self.assertEqual(name, "12345")
        self.assertEqual(dbxrefs, ["GI:12345"])

    def test_gi_only(self):
        """Parse header with only gi identifier."""
        id_, name, dbxrefs = _parse_ncbi_fasta_header("gi|12345 some protein")
        self.assertEqual(id_, "gi|12345")
        self.assertEqual(name, "12345")
        self.assertEqual(dbxrefs, ["GI:12345"])

    def test_multiple_identifiers(self):
        """Parse header with multiple chained identifiers."""
        id_, name, dbxrefs = _parse_ncbi_fasta_header("gi|10|pdb|1A2B|C fake protein")
        self.assertIn("GI:10", dbxrefs)
        self.assertIn("PDB:1A2B", dbxrefs)
        self.assertEqual(len(dbxrefs), 2)


class TestFastaNcbiIterator(unittest.TestCase):
    """Tests for FastaNcbiIterator via SeqIO.parse."""

    def setUp(self):
        """Load records from test file."""
        self.records = list(SeqIO.parse("Fasta/ncbi_headers.fasta", "fasta-ncbi"))

    def test_record_count(self):
        """Correct number of records parsed."""
        self.assertEqual(len(self.records), 10)

    def test_gi_ddbj_record(self):
        """First record: gi + DDBJ identifiers."""
        rec = self.records[0]
        self.assertEqual(rec.id, "gi|3298468|dbj|BAA31520.1|")
        self.assertEqual(rec.name, "3298468")
        self.assertIn("GI:3298468", rec.dbxrefs)
        self.assertIn("DDBJ:BAA31520.1", rec.dbxrefs)
        self.assertEqual(str(rec.seq), "GGHVNPAVTFGAFVGGNITLLRGIVYIIAQLLGSTVACLL")

    def test_gi_genbank_record(self):
        """Second record: gi + GenBank identifiers."""
        rec = self.records[1]
        self.assertEqual(rec.id, "gi|186972394|gb|EU490707.1|")
        self.assertIn("GenBank:EU490707.1", rec.dbxrefs)

    def test_swissprot_record(self):
        """Third record: Swiss-Prot identifier."""
        rec = self.records[2]
        self.assertEqual(rec.name, "P05698")
        self.assertEqual(rec.dbxrefs, ["UniProtKB/Swiss-Prot:P05698"])

    def test_pdb_record(self):
        """Fourth record: PDB identifier with chain."""
        rec = self.records[3]
        self.assertEqual(rec.name, "1A2B")
        self.assertEqual(rec.dbxrefs, ["PDB:1A2B"])

    def test_refseq_record(self):
        """Fifth record: RefSeq identifier."""
        rec = self.records[4]
        self.assertEqual(rec.name, "NM_001301717.2")
        self.assertEqual(rec.dbxrefs, ["RefSeq:NM_001301717.2"])

    def test_trembl_record(self):
        """Sixth record: TrEMBL identifier."""
        rec = self.records[5]
        self.assertEqual(rec.name, "A0A0C5B5G6")
        self.assertEqual(rec.dbxrefs, ["UniProtKB/TrEMBL:A0A0C5B5G6"])

    def test_plain_record(self):
        """Seventh record: no pipes, empty dbxrefs."""
        rec = self.records[6]
        self.assertEqual(rec.id, "plain_id")
        self.assertEqual(rec.name, "plain_id")
        self.assertEqual(rec.dbxrefs, [])

    def test_embl_record(self):
        """Eighth record: EMBL identifier."""
        rec = self.records[7]
        self.assertEqual(rec.name, "CAA12345.6")
        self.assertEqual(rec.dbxrefs, ["EMBL:CAA12345.6"])

    def test_description_preserved(self):
        """Full original title is preserved in description."""
        rec = self.records[1]
        self.assertEqual(
            rec.description,
            "gi|186972394|gb|EU490707.1| Selenipedium aequinoctiale matK",
        )


class TestFastaNcbiRoundtrip(unittest.TestCase):
    """Verify lossless FASTA roundtripping."""

    def test_roundtrip(self):
        """Writing fasta-ncbi records back as fasta preserves headers."""
        records = list(SeqIO.parse("Fasta/ncbi_headers.fasta", "fasta-ncbi"))
        output = StringIO()
        SeqIO.write(records, output, "fasta")
        output.seek(0)

        # Re-parse the written output with plain fasta
        reread = list(SeqIO.parse(output, "fasta"))
        self.assertEqual(len(records), len(reread))
        for orig, back in zip(records, reread):
            self.assertEqual(str(orig.seq), str(back.seq))
            self.assertEqual(orig.description, back.description)


class TestFastaNcbiVsPlain(unittest.TestCase):
    """Compare fasta-ncbi with plain fasta parser."""

    def test_same_sequences(self):
        """Both parsers yield identical sequences."""
        ncbi_records = list(SeqIO.parse("Fasta/ncbi_headers.fasta", "fasta-ncbi"))
        plain_records = list(SeqIO.parse("Fasta/ncbi_headers.fasta", "fasta"))
        self.assertEqual(len(ncbi_records), len(plain_records))
        for ncbi, plain in zip(ncbi_records, plain_records):
            self.assertEqual(str(ncbi.seq), str(plain.seq))
            self.assertEqual(ncbi.id, plain.id)
            self.assertEqual(ncbi.description, plain.description)

    def test_dbxrefs_populated(self):
        """fasta-ncbi populates dbxrefs, plain fasta does not."""
        ncbi_records = list(SeqIO.parse("Fasta/ncbi_headers.fasta", "fasta-ncbi"))
        plain_records = list(SeqIO.parse("Fasta/ncbi_headers.fasta", "fasta"))
        # At least some ncbi records should have dbxrefs
        ncbi_with_xrefs = [r for r in ncbi_records if r.dbxrefs]
        self.assertGreater(len(ncbi_with_xrefs), 0)
        # Plain parser never sets dbxrefs
        for r in plain_records:
            self.assertEqual(r.dbxrefs, [])


class TestExistingFastaFiles(unittest.TestCase):
    """Run fasta-ncbi parser against existing test files with NCBI headers."""

    def test_aster_pro(self):
        """Parse Fasta/aster.pro which has gi|...|dbj|...| headers."""
        records = list(SeqIO.parse("Fasta/aster.pro", "fasta-ncbi"))
        self.assertEqual(len(records), 1)
        rec = records[0]
        self.assertIn("GI:3298468", rec.dbxrefs)
        self.assertIn("DDBJ:BAA31520.1", rec.dbxrefs)
        self.assertTrue(str(rec.seq).startswith("GGHVNPAVTFG"))


class TestFastaNcbiEdgeCases(unittest.TestCase):
    """Edge cases for the NCBI FASTA parser."""

    def test_empty_file(self):
        """Empty file yields no records."""
        records = list(SeqIO.parse(StringIO(""), "fasta-ncbi"))
        self.assertEqual(len(records), 0)

    def test_single_record(self):
        """Parse a single NCBI-style record."""
        data = ">gi|12345|gb|AB000000.1| test\nACGT\n"
        records = list(SeqIO.parse(StringIO(data), "fasta-ncbi"))
        self.assertEqual(len(records), 1)
        self.assertIn("GI:12345", records[0].dbxrefs)
        self.assertIn("GenBank:AB000000.1", records[0].dbxrefs)

    def test_non_fasta_raises(self):
        """Non-FASTA input raises ValueError."""
        with self.assertRaises(ValueError):
            list(SeqIO.parse(StringIO("not fasta"), "fasta-ncbi"))

    def test_multiline_sequence(self):
        """Multi-line wrapped sequence is concatenated correctly."""
        data = ">sp|P12345|TEST_HUMAN test\nACGT\nTGCA\nAAAA\n"
        records = list(SeqIO.parse(StringIO(data), "fasta-ncbi"))
        self.assertEqual(str(records[0].seq), "ACGTTGCAAAAA")

    def test_alphabet_raises(self):
        """Passing alphabet raises ValueError."""
        with self.assertRaises(ValueError):
            FastaNcbiIterator(StringIO(">test\nACGT\n"), alphabet="DNA")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
