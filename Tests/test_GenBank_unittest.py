# Copyright 2013 by Kai Blin.
# Revisions copyright 2015-2016 by Peter Cock.
# Revisions copyright 2019 by Sergio Valqui.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import unittest
from os import path
import warnings
from datetime import datetime
import sys

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

    def test_001_genbank_bad_origin_wrapping_location(self):
        """Bad origin wrapping."""
        with warnings.catch_warnings():
            warnings.simplefilter("error", BiopythonParserWarning)
            try:
                SeqIO.read(path.join("GenBank", "bad_origin_wrap_linear.gb"), "genbank")
            except BiopythonParserWarning as e:
                self.assertEqual(str(e), "Ignoring invalid location: '6801..100'")
            else:
                self.assertTrue(False, "Expected specified BiopythonParserWarning here.")

    def test_implicit_orign_wrap_fix(self):
        """Attempt to fix implied origin wrapping."""
        with warnings.catch_warnings():
            warnings.simplefilter("error", BiopythonParserWarning)
            try:
                SeqIO.read(path.join("GenBank", "bad_origin_wrap.gb"), "genbank")
            except BiopythonParserWarning as e:
                self.assertEqual(str(e), "Attempting to fix invalid location "
                                         "'6801..100' ""as it looks like "
                                         "incorrect origin wrapping. Please "
                                         "fix input file, this could have "
                                         "unintended behavior.")
            else:
                self.assertTrue(False, "Expected specified BiopythonParserWarning here.")

    def test_implicit_orign_wrap_extract_and_translate(self):
        """Test that features wrapped around origin give expected data."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonParserWarning)
            with open(path.join("GenBank", "bad_origin_wrap_CDS.gb")) as handle:
                seq_record = SeqIO.read(handle, "genbank")
                seq_features = seq_record.features
                self.assertEqual(str(seq_features[1].extract(seq_record).seq.lower()),
                                 "atgccctataaaacccagggctgccttggaaaaggcgcaacccc"
                                 "aaccccctcgagccgcggcatataa")
                self.assertEqual(str(seq_features[2].extract(seq_record).seq.lower()),
                                 "atgccgcggctcgagggggttggggttgcgccttttccaaggca"
                                 "gccctgggttttatag")
                self.assertEqual(str(seq_features[1].extract(seq_record).seq.translate()),
                                 "MPYKTQGCLGKGATPTPSSRGI*")
                self.assertEqual(str(seq_features[2].extract(seq_record).seq.translate()),
                                 "MPRLEGVGVAPFPRQPWVL*")

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
        """Parse GenBank record with old DBLINK project entry."""
        record = SeqIO.read("GenBank/NC_005816.gb", "gb")
        self.assertEqual(record.dbxrefs, ["Project:58037"])
        gb = record.format("gb")
        self.assertIn("\nDBLINK      Project: 58037\n", gb)
        embl = record.format("embl")
        self.assertIn("XX\nPR   Project:58037;\nXX\n", embl)

    def test_dblink_two(self):
        """Parse GenBank record with old and new DBLINK project entries."""
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
        """Parse GenBank/EMBL paired records with PR project entry: GenBank."""
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
        """Parse GenBank/EMBL paired records with PR project entry: EMBL."""
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
        """Test if chromosome topology is conserved."""
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

    def test_qualifier_escaping_read(self):
        """Check qualifier escaping is preserved when parsing."""
        # Make sure parsing improperly escaped qualifiers raises a warning
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            record = SeqIO.read("GenBank/qualifier_escaping_read.gb", "gb")
            self.assertEqual(len(caught), 4)
            self.assertEqual(caught[0].category, BiopythonParserWarning)
            self.assertEqual(str(caught[0].message), 'The NCBI states double-quote characters like " should be escaped'
                                                     ' as "" (two double - quotes), but here it was not: '
                                                     '%r' % 'One missing ""quotation mark" here')
        # Check records parsed as expected
        f1 = record.features[0]
        f2 = record.features[1]
        f3 = record.features[2]
        f4 = record.features[3]
        f5 = record.features[4]
        self.assertEqual(f1.qualifiers['note'][0], '"This" is "already" "escaped"')
        self.assertEqual(f2.qualifiers['note'][0], 'One missing "quotation mark" here')
        self.assertEqual(f3.qualifiers['note'][0], 'End not properly "escaped"')
        self.assertEqual(f4.qualifiers['note'][0], '"Start" not properly escaped')
        self.assertEqual(f5.qualifiers['note'][0], 'Middle not "properly" escaped')

    def test_qualifier_escaping_write(self):
        """Check qualifier escaping is preserved when writing."""
        # Write some properly escaped qualifiers and test
        genbank_out = "GenBank/qualifier_escaping_write.gb"
        record = SeqIO.read(genbank_out, "gb")
        f1 = record.features[0]
        f2 = record.features[1]
        f1.qualifiers['note'][0] = '"Should" now "be" escaped in "file"'
        f2.qualifiers['note'][0] = '"Should also be escaped in file"'
        SeqIO.write(record, genbank_out, "gb")
        # Read newly escaped qualifiers and test
        record = SeqIO.read(genbank_out, "gb")
        f1 = record.features[0]
        f2 = record.features[1]
        self.assertEqual(f1.qualifiers['note'][0], '"Should" now "be" escaped in "file"')
        self.assertEqual(f2.qualifiers['note'][0], '"Should also be escaped in file"')

    def test_long_names(self):
        """Various GenBank names which push the column based LOCUS line."""
        original = SeqIO.read("GenBank/iro.gb", "gb")
        self.assertEqual(len(original), 1326)
        # Acceptability of LOCUS line with length > 80 invalidates some of these tests
        for name, seq_len, ok in [
                ("short", 1, True),
                ("max_length_of_16", 1000, True),
                ("overly_long_at_17", 1000, True),
                ("excessively_long_at_22", 99999, True),
                ("excessively_long_at_22", 100000, True),
                ("pushing_the_limits_at_24", 999, True),
                ("pushing_the_limits_at_24", 1000, True),
                ("old_max_name_length_was_26", 10, True),  # 2 digits
                ("old_max_name_length_was_26", 9, True)]:  # 1 digit
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
            # Splitting based on whitespace rather than position due to
            # updated GenBank specification
            name_and_length = line.split()[1:3]
            self.assertEqual(name_and_length, [name, str(seq_len)], line)
            handle.seek(0)
            with warnings.catch_warnings():
                # e.g. BiopythonParserWarning: GenBank LOCUS line
                # identifier over 16 characters
                warnings.simplefilter("ignore", BiopythonWarning)
                new = SeqIO.read(handle, "gb")
            self.assertEqual(name, new.name)
            self.assertEqual(seq_len, len(new))

    def test_genbank_date_default(self):
        """Check if default date is handled correctly."""
        sequence_object = Seq("ATGC", generic_dna)
        # check if default value is inserted correctly
        record = SeqRecord(sequence_object,
                           id='123456789',
                           name='UnitTest',
                           description='Test case for date parsing')
        handle = StringIO()
        SeqIO.write(record, handle, 'genbank')
        handle.seek(0)
        gb = SeqIO.read(handle, "gb")
        self.assertEqual(gb.annotations["date"], "01-JAN-1980")

    def test_genbank_date_correct(self):
        """Check if user provided date is inserted correctly."""
        sequence_object = Seq("ATGC", generic_dna)
        record = SeqRecord(sequence_object,
                           id='123456789',
                           name='UnitTest',
                           description='Test case for date parsing')
        record.annotations["date"] = "24-DEC-2015"
        handle = StringIO()
        SeqIO.write(record, handle, 'genbank')
        handle.seek(0)
        gb = SeqIO.read(handle, "gb")
        self.assertEqual(gb.annotations["date"], "24-DEC-2015")

    def test_genbank_date_list(self):
        """Check if date lists are handled correctly."""
        sequence_object = Seq("ATGC", generic_dna)
        record = SeqRecord(sequence_object,
                           id='123456789',
                           name='UnitTest',
                           description='Test case for date parsing')
        record.annotations["date"] = ["24-DEC-2015"]
        handle = StringIO()
        SeqIO.write(record, handle, 'genbank')
        handle.seek(0)
        gb = SeqIO.read(handle, "gb")
        self.assertEqual(gb.annotations["date"], "24-DEC-2015")

        record = SeqRecord(sequence_object,
                           id='123456789',
                           name='UnitTest',
                           description='Test case for date parsing')
        record.annotations["date"] = ["24-DEC-2015", "25-JAN-2016"]
        handle = StringIO()
        SeqIO.write(record, handle, 'genbank')
        handle.seek(0)
        gb = SeqIO.read(handle, "gb")
        self.assertEqual(gb.annotations["date"], "01-JAN-1980")

    def test_genbank_date_datetime(self):
        """Check if datetime objects are handled correctly."""
        sequence_object = Seq("ATGC", generic_dna)
        record = SeqRecord(sequence_object,
                           id='123456789',
                           name='UnitTest',
                           description='Test case for date parsing')
        record.annotations["date"] = datetime(2000, 2, 2)
        handle = StringIO()
        SeqIO.write(record, handle, 'genbank')
        handle.seek(0)
        gb = SeqIO.read(handle, "gb")
        self.assertEqual(gb.annotations["date"], "02-FEB-2000")

    def test_genbank_date_invalid(self):
        """Check if invalid dates are treated as default."""
        invalid_dates = ("invalid date",
                         "29-2-1981",
                         "35-1-2018",
                         "1-1-80",
                         "1-9-99")

        sequence_object = Seq("ATGC", generic_dna)
        for invalid_date in invalid_dates:
            record = SeqRecord(sequence_object,
                               id='123456789',
                               name='UnitTest',
                               description='Test case for date parsing')

            record.annotations["date"] = invalid_date
            handle = StringIO()
            SeqIO.write(record, handle, 'genbank')
            handle.seek(0)
            gb = SeqIO.read(handle, "gb")
            self.assertEqual(gb.annotations["date"], "01-JAN-1980")

    def test_longer_locus_line(self):
        """Check that we can read and write files with longer locus lines."""
        # Create example file from existing file
        with open(path.join("GenBank", "DS830848.gb"), 'r') as inhandle:
            data = inhandle.readlines()
            data[0] = "LOCUS       AZZZAA021234567891234 2147483647 bp    DNA     linear   PRI 15-OCT-2018\n"

        # Create memory file from modified genbank file
        in_tmp = StringIO()
        in_tmp.writelines(data)
        in_tmp.seek(0)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            in_tmp.seek(0)
            record = SeqIO.read(in_tmp, 'genbank')

            # Create temporary output memory file
            out_tmp = StringIO()
            SeqIO.write(record, out_tmp, 'genbank')

            # Check that the written file can be read back in
            out_tmp.seek(0)
            record_in = SeqIO.read(out_tmp, 'genbank')
            self.assertEqual(record_in.id, "DS830848.1")
            self.assertEqual(record_in.name, "AZZZAA021234567891234")
            self.assertEqual(len(record_in.seq), 2147483647)

    if sys.maxsize > 2147483647:
        def test_extremely_long_sequence(self):
            """Tests if extremely long sequences can be read.

            This is only run if sys.maxsize > 2147483647.
            """
            # Create example file from existing file
            with open(path.join("GenBank", "DS830848.gb"), 'r') as inhandle:
                data = inhandle.readlines()
                data[0] = "LOCUS       AZZZAA02123456789 10000000000 bp    DNA     linear   PRI 15-OCT-2018\n"

            # Create memory file from modified genbank file
            in_tmp = StringIO()
            in_tmp.writelines(data)
            in_tmp.seek(0)

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                in_tmp.seek(0)
                record = SeqIO.read(in_tmp, 'genbank')

                # Create temporary output memory file
                out_tmp = StringIO()
                SeqIO.write(record, out_tmp, 'genbank')

                # Check that the written file can be read back in
                out_tmp.seek(0)
                record_in = SeqIO.read(out_tmp, 'genbank')
                self.assertEqual(record_in.id, "DS830848.1")
                self.assertEqual(record_in.name, "AZZZAA02123456789")
                self.assertEqual(len(record_in.seq), 10000000000)

            def read_longer_than_maxsize():
                with open(path.join("GenBank", "DS830848.gb"), 'r') as inhandle:
                    data2 = inhandle.readlines()
                    data2[0] = "LOCUS       AZZZAA02123456789 " + str(sys.maxsize + 1) + " bp    DNA     linear   PRI 15-OCT-2018\n"

                long_in_tmp = StringIO()
                long_in_tmp.writelines(data2)
                long_in_tmp.seek(0)
                record = SeqIO.read(long_in_tmp, 'genbank')

            self.assertRaises(ValueError, read_longer_than_maxsize)


class LineOneTests(unittest.TestCase):
    """Check GenBank/EMBL topology / molecule_type parsing."""

    def test_topology_genbank(self):
        """Check GenBank LOCUS line parsing."""
        # This is a bit low level, but can test pasing the LOCUS line only
        tests = [
            ("LOCUS       U00096",
             None, None, None, None),
            # This example is actually fungal, accession U49845 from Saccharomyces cerevisiae:
            ("LOCUS       SCU49845     5028 bp    DNA             PLN       21-JUN-1999",
             None, "DNA", "PLN", None),
            ("LOCUS       AB070938                6497 bp    DNA     linear   BCT 11-OCT-2001",
             "linear", "DNA", "BCT", None),
            ("LOCUS       NC_005816               9609 bp    DNA     circular BCT 21-JUL-2008",
             "circular", "DNA", "BCT", None),
            ("LOCUS       SCX3_BUTOC                64 aa            linear   INV 16-OCT-2001",
             "linear", None, "INV", None),
            ("LOCUS       pEH010                  5743 bp    DNA     circular",
             "circular", "DNA", None, [BiopythonParserWarning]),
            # This is a test of the format > 80 chars long
            ("LOCUS       AZZZAA02123456789 1000000000 bp    DNA     linear   PRI 15-OCT-2018",
             "linear", "DNA", "PRI", None)
        ]
        for (line, topo, mol_type, div, warning_list) in tests:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
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
                if warning_list is None:
                    self.assertEqual(len(caught), 0)
                else:
                    self.assertEqual(len(caught), len(warning_list))
                    for i, warning_class in enumerate(warning_list):
                        self.assertEqual(caught[i].category, warning_class)

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

    def test_seqrecord_default_description(self):
        """Read in file using SeqRecord default description."""
        old = SeqRecord(Seq("ACGT", generic_dna),
                        id="example",
                        name="short")
        self.assertEqual(old.description, "<unknown description>")
        txt = old.format("gb")
        self.assertIn("DEFINITION  .\n", txt)
        new = SeqIO.read(StringIO(txt), "gb")
        self.assertEqual(old.id, new.id)
        self.assertEqual(old.name, new.name)
        self.assertEqual("", new.description)
        self.assertEqual(old.seq, new.seq)

    # Evil hack with 000 to manipulate sort order to ensure this is tested
    # first (otherwise something silences the warning)
    def test_000_write_invalid_but_parsed_locus_line(self):
        """Make sure we survive writing slightly invalid LOCUS lines we could parse."""
        # grab a valid file
        with open(path.join('GenBank', 'NC_005816.gb'), 'r') as handle:
            lines = handle.readlines()

        # futz with the molecule type to make it lower case
        invalid_line = "LOCUS       NC_005816               9609 bp    dna     circular BCT 21-JUL-2008\n"
        lines[0] = invalid_line
        fake_handle = StringIO("".join(lines))

        # Make sure parsing this actually raises a warning
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            rec = SeqIO.read(fake_handle, 'genbank')
            self.assertEqual(len(caught), 1)
            self.assertEqual(caught[0].category, BiopythonParserWarning)
            self.assertEqual(str(caught[0].message), "Non-upper case molecule type in LOCUS line: dna")

        out_handle = StringIO()

        ret = SeqIO.write([rec], out_handle, 'genbank')
        self.assertEqual(ret, 1)

        out_handle.seek(0)
        out_lines = out_handle.readlines()
        self.assertEqual(out_lines[0], invalid_line)


class GenBankScannerTests(unittest.TestCase):
    """GenBank Scanner tests, test parsing gbk and embl files."""

    gb_s = Scanner.GenBankScanner()

    def gb_to_l_cds_f(self, filename, tags2id=None):
        """Gb file to Seq list parse CDS features."""
        with open(filename) as handle:
            if tags2id:
                l_cds_f = list(self.gb_s.parse_cds_features(handle, tags2id=tags2id))
            else:
                l_cds_f = list(self.gb_s.parse_cds_features(handle))
        return l_cds_f

    def gb_to_l_r(self, filename, do_features=False):
        """Gb file to Seq list parse records."""
        with open(filename) as handle:
            l_gb_r = list(self.gb_s.parse_records(handle, do_features=do_features))
        return l_gb_r

    def test_genbank_cds_interaction(self):
        """Test CDS interaction, parse CDS features on gb(k) files."""
        # Test parse CDS features on NC_000932.gb
        l_cds_f = self.gb_to_l_cds_f("GenBank/NC_000932.gb")
        # number of records, should be 85
        self.assertEqual(len(l_cds_f), 85)
        # Seq ID
        self.assertEqual(l_cds_f[0].id, 'NP_051037.1')
        self.assertEqual(l_cds_f[84].id, 'NP_051123.1')

        # Test parse CDS features on NC_005816.gb, Tag to ID
        l_cds_f = self.gb_to_l_cds_f("GenBank/NC_005816.gb",
                                     tags2id=('gene', 'locus_tag', 'product'))
        # number of records, should be 10
        self.assertEqual(len(l_cds_f), 10)
        # Seq ID
        self.assertEqual(l_cds_f[0].id, '<unknown id>')
        self.assertEqual(l_cds_f[0].name, 'YP_pPCP01')

        # Test parse CDS features on NC_000932.gb and NC_005816.gb combined
        l_cds_f1 = self.gb_to_l_cds_f("GenBank/NC_000932.gb",
                                      tags2id=('gene', 'locus_tag', 'product'))
        l_cds_f2 = self.gb_to_l_cds_f("GenBank/NC_005816.gb",
                                      tags2id=('gene', 'locus_tag', 'product'))
        l_cds_combined = l_cds_f1 + l_cds_f2
        # number of records combined, should be 95
        self.assertEqual(len(l_cds_combined), 95)
        # Seq ID
        self.assertEqual(l_cds_combined[0].id, 'rps12')
        self.assertEqual(l_cds_combined[0].description, 'ribosomal protein S12')
        self.assertEqual(l_cds_combined[94].id, '<unknown id>')
        self.assertEqual(l_cds_combined[94].description, 'hypothetical protein')

    def test_genbank_interaction(self):
        """Test GenBank records interaction on gbk files."""
        # Test parse records, on NC_005816, do_features False
        l_r = self.gb_to_l_r("GenBank/NC_005816.gb", do_features=False)
        # number of records, should be 1
        self.assertEqual(len(l_r), 1)
        self.assertEqual(l_r[0].id, 'NC_005816.1')
        self.assertEqual(l_r[0].name, 'NC_005816')
        self.assertEqual(l_r[0].description, 'Yersinia pestis biovar '
                                             'Microtus str. 91001 plasmid '
                                             'pPCP1, complete sequence')
        self.assertEqual(len(l_r[0].features), 0)

        # Test parse records on NC_005816, do_features True
        l_r = self.gb_to_l_r("GenBank/NC_005816.gb", do_features=True)
        # number of records, should be 1
        self.assertEqual(len(l_r), 1)
        self.assertEqual(l_r[0].id, 'NC_005816.1')
        self.assertEqual(l_r[0].name, 'NC_005816')
        self.assertEqual(l_r[0].description, 'Yersinia pestis biovar '
                                             'Microtus str. 91001 plasmid '
                                             'pPCP1, complete sequence')
        self.assertEqual(len(l_r[0].features), 41)

        # Test parse records on "GenBank/NC_000932.gb", do_features False
        l_r = self.gb_to_l_r("GenBank/NC_000932.gb", do_features=False)
        # number of records, should be 1
        self.assertEqual(len(l_r), 1)
        self.assertEqual(l_r[0].id, 'NC_000932.1')
        self.assertEqual(l_r[0].name, 'NC_000932')
        self.assertEqual(l_r[0].description, 'Arabidopsis thaliana chloroplast, '
                                             'complete genome')
        self.assertEqual(len(l_r[0].features), 0)

        # Test parse records on NC_000932, do_features True
        l_r = self.gb_to_l_r("GenBank/NC_000932.gb", do_features=True)
        # number of records, should be 1
        self.assertEqual(len(l_r), 1)
        self.assertEqual(l_r[0].id, 'NC_000932.1')
        self.assertEqual(l_r[0].name, 'NC_000932')
        self.assertEqual(l_r[0].description, 'Arabidopsis thaliana chloroplast, '
                                             'complete genome')
        self.assertEqual(len(l_r[0].features), 259)

    def test_embl_cds_interaction(self):
        """Test EMBL CDS interaction, parse CDS features on embl files."""
        embl_s = Scanner.EmblScanner()

        # Test parse CDS features on embl_file
        with open("EMBL/AE017046.embl") as handle_embl7046:
            l_cds_f = list(embl_s.parse_cds_features(handle_embl7046))
        # number of records, should be 10
        self.assertEqual(len(l_cds_f), 10)
        # Seq ID
        self.assertEqual(l_cds_f[0].id, 'AAS58758.1')
        self.assertEqual(l_cds_f[0].description, 'putative transposase')

    def test_embl_record_interaction(self):
        """Test EMBL Record interaction on embl files."""
        embl_s = Scanner.EmblScanner()

        #  Test parse records on embl_file
        with open("EMBL/AE017046.embl") as handle_embl7046:
            l_embl_r = list(embl_s.parse_records(handle_embl7046, do_features=True))
        # number of records, should be 1
        self.assertEqual(len(l_embl_r), 1)
        self.assertEqual(l_embl_r[0].id, 'AE017046.1')
        self.assertEqual(l_embl_r[0].description, 'Yersinia pestis biovar Microtus '
                                                  'str. 91001 plasmid pPCP1, complete '
                                                  'sequence.')
        self.assertEqual(len(l_embl_r[0].features), 29)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
