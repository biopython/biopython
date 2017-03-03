# Copyright 2013 by Kai Blin.
# Revisions copyright 2015-2016 by Peter Cock.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import unittest
import warnings
from os import path

from Bio import BiopythonParserWarning
from Bio import BiopythonWarning
from Bio import GenBank
from Bio.GenBank import Scanner
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio._py3k import StringIO


class GenBankTests(unittest.TestCase):
    """GenBank tests."""

    def test_invalid_product_line_raises_value_error(self):
        """Parsing invalid product line."""
        def parse_invalid_product_line():
            rec = SeqIO.read(path.join('GenBank', 'invalid_product.gb'),
                             'genbank')
        self.assertRaises(ValueError, parse_invalid_product_line)

    def test_genbank_read(self):
        """GenBank.read(...) simple test."""
        with open(path.join("GenBank", "NC_000932.gb")) as handle:
            record = GenBank.read(handle)
        self.assertEqual(['NC_000932'], record.accession)

    def test_genbank_read_multirecord(self):
        """GenBank.read(...) error on multiple record input."""
        with open(path.join("GenBank", "cor6_6.gb")) as handle:
            self.assertRaises(ValueError, GenBank.read, handle)

    def test_genbank_read_invalid(self):
        """GenBank.read(...) error on invalid file (e.g. FASTA file)."""
        with open(path.join("GenBank", "NC_000932.faa")) as handle:
            self.assertRaises(ValueError, GenBank.read, handle)

    def test_genbank_read_no_origin_no_end(self):
        """GenBank.read(...) error on malformed file."""
        with open(path.join("GenBank", "no_origin_no_end.gb")) as handle:
            self.assertRaises(ValueError, GenBank.read, handle)

    # Evil hack with 000 to manipulate sort order to ensure this is tested
    # first (otherwise something silences the warning)
    def test_000_genbank_bad_loc_wrap_warning(self):
        """Feature line wrapping warning."""
        with warnings.catch_warnings():
            warnings.simplefilter("error", BiopythonParserWarning)
            with open(path.join("GenBank", "bad_loc_wrap.gb")) as handle:
                # self.assertRaises(BiopythonParserWarning, GenBank.read, handle)
                try:
                    record = GenBank.read(handle)
                except BiopythonParserWarning as e:
                    self.assertEqual(str(e),
                                     "Non-standard feature line wrapping (didn't break on comma)?")
                else:
                    self.assertTrue(False, "Expected specified BiopythonParserWarning here.")

    # Similar hack as we also want to catch that warning here
    def test_001_negative_location_warning(self):
        """Un-parsable feature location warning."""
        with warnings.catch_warnings():
            warnings.simplefilter("error", BiopythonParserWarning)
            try:
                SeqIO.read(path.join("GenBank", "negative_location.gb"), "genbank")
            except BiopythonParserWarning as e:
                self.assertEqual(str(e), "Couldn't parse feature location: '-2..492'")
            else:
                self.assertTrue(False, "Expected specified BiopythonParserWarning here.")

    def test_genbank_bad_loc_wrap_parsing(self):
        """Bad location wrapping."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonParserWarning)
            with open(path.join("GenBank", "bad_loc_wrap.gb")) as handle:
                record = GenBank.read(handle)
                self.assertEqual(1, len(record.features))
                loc = record.features[0].location
                self.assertEqual(loc, "join(3462..3615,3698..3978,4077..4307,4408..4797,4876..5028,5141..5332)")

    def test_negative_location(self):
        """Negative feature locations."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonParserWarning)
            rec = SeqIO.read(path.join("GenBank", "negative_location.gb"), "genbank")
            self.assertEqual(None, rec.features[-1].location)

    def test_dot_lineage(self):
        """Missing taxonomy lineage."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonParserWarning)
            rec = SeqIO.read("GenBank/bad_loc_wrap.gb", "genbank")
        self.assertEqual(rec.annotations["organism"], ".")
        self.assertEqual(rec.annotations["taxonomy"], [])

    def test_dblink(self):
        """GenBank record with old DBLINK project entry."""
        record = SeqIO.read("GenBank/NC_005816.gb", "gb")
        self.assertEqual(record.dbxrefs, ["Project:58037"])
        gb = record.format("gb")
        self.assertIn("\nDBLINK      Project: 58037\n", gb)
        embl = record.format("embl")
        self.assertIn("XX\nPR   Project:58037;\nXX\n", embl)

    def test_dblink_two(self):
        """GenBank record with old and new DBLINK project entries."""
        record = SeqIO.read("GenBank/NP_416719.gbwithparts", "gb")
        self.assertEqual(record.dbxrefs,
                         ["Project:57779", "BioProject:PRJNA57779"])
        gb = record.format("gb")
        self.assertTrue("""
DBLINK      Project: 57779
            BioProject: PRJNA57779
KEYWORDS    """ in gb, gb)
        embl = record.format("embl")
        self.assertIn("XX\nPR   Project:PRJNA57779;\nXX\n", embl)

    def test_dbline_gb_embl(self):
        """GenBank / EMBL paired records with PR project entry: GenBank"""
        record = SeqIO.read("GenBank/DS830848.gb", "gb")
        self.assertIn("BioProject:PRJNA16232", record.dbxrefs)
        gb = record.format("gb")
        self.assertTrue("""
DBLINK      BioProject: PRJNA16232
            BioSample: SAMN03004382
KEYWORDS    """ in gb, gb)
        # Also check EMBL output
        embl = record.format("embl")
        self.assertIn("XX\nPR   Project:PRJNA16232;\nXX\n", embl)

    def test_dbline_embl_gb(self):
        """GenBank / EMBL paired records with PR project entry: EMBL"""
        record = SeqIO.read("EMBL/DS830848.embl", "embl")
        # TODO: Should we map this to BioProject:PRJNA16232
        self.assertIn("Project:PRJNA16232", record.dbxrefs)
        gb = record.format("gb")
        self.assertTrue("""
DBLINK      Project: PRJNA16232
            MD5: 387e72e4f7ae804780d06f875ab3bc41
            ENA: ABJB010000000
            ENA: ABJB000000000
            BioSample: SAMN03004382
KEYWORDS    """ in gb, gb)
        embl = record.format("embl")
        self.assertIn("XX\nPR   Project:PRJNA16232;\nXX\n", embl)

    def test_structured_comment_parsing(self):
        """Structued comment parsing."""
        # GISAID_EpiFlu(TM)Data, HM138502.gbk has both 'comment' and 'structured_comment'
        record = SeqIO.read(path.join('GenBank', 'HM138502.gbk'), 'genbank')
        self.assertEqual(record.annotations['comment'],
                         "Swine influenza A (H1N1) virus isolated during human swine flu\noutbreak of 2009.")
        self.assertEqual(record.annotations['structured_comment']['GISAID_EpiFlu(TM)Data']['Lineage'], 'swl')
        self.assertEqual(len(record.annotations['structured_comment']['GISAID_EpiFlu(TM)Data']), 3)
        with open(path.join('GenBank', 'HM138502_output.gbk'), "r") as ifile:
            self.assertEqual(record.format("gb"), ifile.read())
        # FluData structured comment
        record = SeqIO.read(path.join('GenBank', 'EU851978.gbk'), 'genbank')
        self.assertEqual(record.annotations['structured_comment']['FluData']['LabID'], '2008704957')
        self.assertEqual(len(record.annotations['structured_comment']['FluData']), 5)
        with open(path.join('GenBank', 'EU851978_output.gbk'), "r") as ifile:
            self.assertEqual(record.format("gb"), ifile.read())
        # Assembly-Data structured comment
        record = SeqIO.read(path.join('GenBank', 'KF527485.gbk'), 'genbank')
        self.assertEqual(record.annotations['structured_comment']['Assembly-Data']['Assembly Method'], 'Lasergene v. 10')
        self.assertEqual(len(record.annotations['structured_comment']['Assembly-Data']), 2)
        with open(path.join('GenBank', 'KF527485_output.gbk'), "r") as ifile:
            self.assertEqual(record.format("gb"), ifile.read())
        # No structured comment in NC_000932.gb, just a regular comment
        record = SeqIO.read(path.join('GenBank', 'NC_000932.gb'), 'genbank')
        self.assertFalse("structured_comment" in record.annotations)
        self.assertEqual(record.annotations['comment'],
                         'REVIEWED REFSEQ: This record has been curated by NCBI staff. The\n'
                         'reference sequence was derived from AP000423.\n'
                         'COMPLETENESS: full length.')

    def test_locus_line_topogoly(self):
        """Test if chromosome topology is conserved"""
        record = SeqIO.read('GenBank/DS830848.gb', 'genbank')
        self.assertEqual(record.annotations['topology'], 'linear')
        out_handle = StringIO()
        SeqIO.write([record], out_handle, 'genbank')
        first_line = out_handle.getvalue().split('\n')[0]
        self.assertIn('linear', first_line)
        with open('GenBank/DS830848.gb', 'r') as fh:
            orig_first_line = fh.readline().strip()
        self.assertEqual(first_line, orig_first_line)

    def test_qualifier_order(self):
        """Check the qualifier order is preserved."""
        record = SeqIO.read("GenBank/DS830848.gb", "gb")
        f = record.features[0]
        self.assertEqual(list(f.qualifiers),
                         ['organism', 'mol_type', 'strain', 'db_xref', 'dev_stage'])

    def test_long_names(self):
        """Various GenBank names which push the column based LOCUS line."""
        original = SeqIO.read("GenBank/iro.gb", "gb")
        self.assertEqual(len(original), 1326)
        for name, seq_len, ok in [
                ("short", 1, True),
                ("max_length_of_16", 1000, True),
                ("overly_long_at_17", 1000, True),
                ("excessively_long_at_22", 99999, True),
                ("excessively_long_at_22", 100000, False),
                ("pushing_the_limits_at_24", 999, True),
                ("pushing_the_limits_at_24", 1000, False),
                ("longest_possible_len_of_26", 10, False),  # 2 digits
                ("longest_possible_len_of_26", 9, True),  # 1 digit
                ]:
            # Make the length match the desired target
            record = original[:]
            # TODO - Implement Seq * int
            record.seq = Seq("N" * seq_len, original.seq.alphabet)
            # Set the identifer to the desired name
            record.id = record.name = name
            # Attempt to output the record...
            if not ok:
                # e.g. ValueError: Locus identifier 'excessively_long_at_22' is too long
                self.assertRaises(ValueError, record.format, "gb")
                continue
            with warnings.catch_warnings():
                # e.g. BiopythonWarning: Stealing space from length field to
                # allow long name in LOCUS line
                warnings.simplefilter("ignore", BiopythonWarning)
                # output = record.format("gb")
                handle = StringIO()
                self.assertEqual(1, SeqIO.write(record, handle, "gb"))
            handle.seek(0)
            line = handle.readline()
            self.assertIn(" %s " % name, line)
            self.assertIn(" %i bp " % seq_len, line)
            name_and_length = line[12:40]
            self.assertEqual(name_and_length.split(), [name, str(seq_len)], line)
            handle.seek(0)
            with warnings.catch_warnings():
                # e.g. BiopythonParserWarning: GenBank LOCUS line
                # identifier over 16 characters
                warnings.simplefilter("ignore", BiopythonWarning)
                new = SeqIO.read(handle, "gb")
            self.assertEqual(name, new.name)
            self.assertEqual(seq_len, len(new))


class LineOneTests(unittest.TestCase):
    """Check GenBank/EMBL topology / molecule_type parsing."""

    def test_topology_genbank(self):
        """Check GenBank LOCUS line parsing."""
        # This is a bit low level, but can test pasing the LOCUS line only
        tests = [
            ("LOCUS       U00096",
             None, None, None),
            # This example is actually fungal, accession U49845 from Saccharomyces cerevisiae:
            ("LOCUS       SCU49845     5028 bp    DNA             PLN       21-JUN-1999",
             None, "DNA", "PLN"),
            ("LOCUS       AB070938                6497 bp    DNA     linear   BCT 11-OCT-2001",
             "linear", "DNA", "BCT"),
            ("LOCUS       NC_005816               9609 bp    DNA     circular BCT 21-JUL-2008",
             "circular", "DNA", "BCT"),
            ("LOCUS       SCX3_BUTOC                64 aa            linear   INV 16-OCT-2001",
             "linear", None, "INV"),
        ]
        for (line, topo, mol_type, div) in tests:
            scanner = Scanner.GenBankScanner()
            consumer = GenBank._FeatureConsumer(1, GenBank.FeatureValueCleaner)
            scanner._feed_first_line(consumer, line)
            t = consumer.data.annotations.get('topology', None)
            self.assertEqual(t, topo,
                             "Wrong topology %r not %r from %r" % (t, topo, line))
            mt = consumer.data.annotations.get('molecule_type', None)
            self.assertEqual(mt, mol_type,
                             "Wrong molecule_type %r not %r from %r" %
                             (mt, mol_type, line))
            d = consumer.data.annotations.get('data_file_division', None)
            self.assertEqual(d, div,
                             "Wrong division %r not %r from %r" % (d, div, line))

    def test_topology_embl(self):
        """Check EMBL ID line parsing."""
        # This is a bit low level, but can test pasing the ID line only
        tests = [
            # Modern examples with sequence version
            ("ID   X56734; SV 1; linear; mRNA; STD; PLN; 1859 BP.",
             "linear", "mRNA", "PLN"),
            ("ID   CD789012; SV 4; linear; genomic DNA; HTG; MAM; 500 BP.",
             "linear", "genomic DNA", "MAM"),
            # Example to match GenBank example used above:
            ("ID   U49845; SV 1; linear; genomic DNA; STD; FUN; 5028 BP.",
             "linear", "genomic DNA", "FUN"),
            # Old examples:
            ("ID   BSUB9999   standard; circular DNA; PRO; 4214630 BP.",
             "circular", "DNA", "PRO"),
            ("ID   SC10H5 standard; DNA; PRO; 4870 BP.",
             None, "DNA", "PRO"),
            # Patent example from 2016-06-10
            # ftp://ftp.ebi.ac.uk/pub/databases/embl/patent/
            ("ID   A01679; SV 1; linear; unassigned DNA; PAT; MUS; 12 BP.",
             "linear", "unassigned DNA", "MUS"),
            # Old patent examples
            ("ID   NRP_AX000635; PRT; NR1; 15 SQ", None, None, "NR1"),
            ("ID   NRP0000016E; PRT; NR2; 5 SQ", None, None, "NR2"),
            # KIPO patent examples
            ("ID   DI500001       STANDARD;      PRT;   111 AA.", None, None, None),
            ("ID   DI644510   standard; PRT;  1852 AA.", None, None, None),
        ]
        for (line, topo, mol_type, div) in tests:
            scanner = Scanner.EmblScanner()
            consumer = GenBank._FeatureConsumer(1, GenBank.FeatureValueCleaner)
            scanner._feed_first_line(consumer, line)
            t = consumer.data.annotations.get('topology', None)
            self.assertEqual(t, topo,
                             "Wrong topology %r not %r from %r" % (t, topo, line))
            mt = consumer.data.annotations.get('molecule_type', None)
            self.assertEqual(mt, mol_type,
                             "Wrong molecule_type %r not %r from %r" %
                             (mt, mol_type, line))
            d = consumer.data.annotations.get('data_file_division', None)
            self.assertEqual(d, div,
                             "Wrong division %r not %r from %r" % (d, div, line))

    def test_first_line_imgt(self):
        """Check IMGT ID line parsing."""
        # This is a bit low level, but can test pasing the ID line only
        tests = [
            ("ID   HLA00001   standard; DNA; HUM; 3503 BP.",
             None, "DNA", "HUM"),
            ("ID   HLA00001; SV 1; standard; DNA; HUM; 3503 BP.",
             None, "DNA", "HUM"),
        ]
        for (line, topo, mol_type, div) in tests:
            scanner = Scanner._ImgtScanner()
            consumer = GenBank._FeatureConsumer(1, GenBank.FeatureValueCleaner)
            scanner._feed_first_line(consumer, line)
            t = consumer.data.annotations.get('topology', None)
            self.assertEqual(t, topo,
                             "Wrong topology %r not %r from %r" % (t, topo, line))
            mt = consumer.data.annotations.get('molecule_type', None)
            self.assertEqual(mt, mol_type,
                             "Wrong molecule_type %r not %r from %r" %
                             (mt, mol_type, line))
            d = consumer.data.annotations.get('data_file_division', None)
            self.assertEqual(d, div,
                             "Wrong division %r not %r from %r" % (d, div, line))


class OutputTests(unittest.TestCase):
    """GenBank output tests."""

    def test_mad_dots(self):
        """Writing and reading back accesssion.version variants."""
        for identifier in ["example",
                           "example.1a",
                           "example.1.2",
                           "example.1-2",
                           ]:
            old = SeqRecord(Seq("ACGT", generic_dna),
                            id=identifier,
                            name=identifier,
                            description="mad dots")
            new = SeqIO.read(StringIO(old.format("gb")), "gb")
            self.assertEqual(old.id, new.id)
            self.assertEqual(old.name, new.name)
            self.assertEqual(old.description, new.description)
            self.assertEqual(old.seq, new.seq)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
