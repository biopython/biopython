# Copyright 2009-2013 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Additional unit tests for Bio.SeqIO.QualityIO (covering FASTQ and QUAL)."""


import os
import unittest
import warnings

from io import StringIO
from io import BytesIO

from Bio import BiopythonWarning, BiopythonParserWarning
from Bio.SeqIO import QualityIO
from Bio import SeqIO
from Bio.Seq import Seq, UnknownSeq, MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio.Data.IUPACData import ambiguous_dna_letters, ambiguous_rna_letters

from test_SeqIO import SeqIOTestBaseClass, SeqIOConverterTestBaseClass


class QualityIOTestBaseClass(SeqIOTestBaseClass):
    def compare_record(self, old, new, fmt=None, msg=None):
        """Quality aware SeqRecord comparison.

        This will check the mapping between Solexa and PHRED scores.
        It knows to ignore UnknownSeq objects for string matching (i.e. QUAL files).
        """
        super().compare_record(old, new, msg=None)
        if fmt in ["fastq-solexa", "fastq-illumina"]:
            truncate = 62
        elif fmt in ["fastq", "fastq-sanger"]:
            truncate = 93
        else:
            assert fmt in ["fasta", "qual", "phd", "sff", "tab", None]
            truncate = None
        for keyword in ("phred_quality", "solexa_quality"):
            q_old = old.letter_annotations.get(keyword)
            q_new = new.letter_annotations.get(keyword)
            if q_old is None or q_new is None:
                continue
            if truncate is not None and q_old != q_new:
                q_old = [min(q, truncate) for q in q_old]
                q_new = [min(q, truncate) for q in q_new]
            err_msg = "mismatch in %s" % keyword
            if msg is not None:
                err_msg = "%s: %s" % (msg, err_msg)
            self.assertEqual(q_old, q_new, msg=err_msg)

        q_old = old.letter_annotations.get("phred_quality")
        q_new = new.letter_annotations.get("solexa_quality")
        if q_old is not None and q_new is not None:
            # Mapping from Solexa to PHRED is lossy, but so is PHRED to Solexa.
            # Assume "old" is the original, and "new" has been converted.
            converted = [round(QualityIO.solexa_quality_from_phred(q)) for q in q_old]
            if truncate is not None:
                converted = [min(q, truncate) for q in converted]
            err_msg = "mismatch converting phred_quality %s to solexa_quality" % q_old
            if msg is not None:
                err_msg = "%s: %s" % (msg, err_msg)
            self.assertEqual(converted, q_new, msg=err_msg)

        q_old = old.letter_annotations.get("solexa_quality")
        q_new = new.letter_annotations.get("phred_quality")
        if q_old is not None and q_new is not None:
            # Mapping from Solexa to PHRED is lossy, but so is PHRED to Solexa.
            # Assume "old" is the original, and "new" has been converted.
            converted = [round(QualityIO.phred_quality_from_solexa(q)) for q in q_old]
            if truncate is not None:
                converted = [min(q, truncate) for q in converted]
            err_msg = "mismatch converting solexa_quality %s to phred_quality" % q_old
            if msg is not None:
                err_msg = "%s: %s" % (msg, err_msg)
            self.assertEqual(converted, q_new, msg=err_msg)


class TestFastqErrors(unittest.TestCase):
    """Test reject invalid FASTQ files."""

    def check_fails(self, filename, good_count, formats=None, raw=True):
        if not formats:
            formats = ["fastq-sanger", "fastq-solexa", "fastq-illumina"]
        msg = "SeqIO.parse failed to detect error in %s" % filename
        for fmt in formats:
            records = SeqIO.parse(filename, fmt)
            for i in range(good_count):
                record = next(records)  # Make sure no errors!
                self.assertIsInstance(record, SeqRecord)
            # Detect error in the next record:
            with self.assertRaises(ValueError, msg=msg) as cm:
                record = next(records)

    def check_general_fails(self, filename, good_count):
        tuples = QualityIO.FastqGeneralIterator(filename)
        msg = "FastqGeneralIterator failed to detect error in %s" % filename
        for i in range(good_count):
            title, seq, qual = next(tuples)  # Make sure no errors!
        # Detect error in the next record:
        with self.assertRaises(ValueError, msg=msg) as cm:
            title, seq, qual = next(tuples)

    def check_general_passes(self, filename, record_count):
        tuples = QualityIO.FastqGeneralIterator(filename)
        # This "raw" parser doesn't check the ASCII characters which means
        # certain invalid FASTQ files will get parsed without errors.
        msg = "FastqGeneralIterator failed to parse %s" % filename
        count = 0
        for title, seq, qual in tuples:
            self.assertEqual(len(seq), len(qual), msg=msg)
            count += 1
        self.assertEqual(count, record_count, msg=msg)

    def test_reject_high_and_low(self):
        # These FASTQ files will be rejected by both the low level parser AND
        # the high level SeqRecord parser:
        tests = [
            ("Quality/error_diff_ids.fastq", 2),
            ("Quality/error_no_qual.fastq", 0),
            ("Quality/error_long_qual.fastq", 3),
            ("Quality/error_short_qual.fastq", 2),
            ("Quality/error_double_seq.fastq", 3),
            ("Quality/error_double_qual.fastq", 2),
            ("Quality/error_tabs.fastq", 0),
            ("Quality/error_spaces.fastq", 0),
            ("Quality/error_trunc_in_title.fastq", 4),
            ("Quality/error_trunc_in_seq.fastq", 4),
            ("Quality/error_trunc_in_plus.fastq", 4),
            ("Quality/error_trunc_in_qual.fastq", 4),
            ("Quality/error_trunc_at_seq.fastq", 4),
            ("Quality/error_trunc_at_plus.fastq", 4),
            ("Quality/error_trunc_at_qual.fastq", 4),
        ]
        for path, count in tests:
            self.check_fails(path, count)
            self.check_general_fails(path, count)

    def test_reject_high_but_not_low(self):
        # These FASTQ files which will be rejected by the high level SeqRecord
        # parser, but will be accepted by the low level parser:
        tests = [
            ("Quality/error_qual_del.fastq", 3, 5),
            ("Quality/error_qual_space.fastq", 3, 5),
            ("Quality/error_qual_vtab.fastq", 0, 5),
            ("Quality/error_qual_escape.fastq", 4, 5),
            ("Quality/error_qual_unit_sep.fastq", 2, 5),
            ("Quality/error_qual_tab.fastq", 4, 5),
            ("Quality/error_qual_null.fastq", 0, 5),
        ]
        for path, good_count, full_count in tests:
            self.check_fails(path, good_count)
            self.check_general_passes(path, full_count)


class TestReferenceSffConversions(unittest.TestCase):
    def check(self, sff_name, sff_format, out_name, fmt):
        wanted = list(SeqIO.parse(out_name, fmt))
        data = StringIO()
        count = SeqIO.convert(sff_name, sff_format, data, fmt)
        self.assertEqual(count, len(wanted))
        data.seek(0)
        converted = list(SeqIO.parse(data, fmt))
        self.assertEqual(len(wanted), len(converted))
        for old, new in zip(wanted, converted):
            self.assertEqual(old.id, new.id)
            self.assertEqual(old.name, new.name)
            if fmt != "qual":
                self.assertEqual(str(old.seq), str(new.seq))
            elif fmt != "fasta":
                self.assertEqual(
                    old.letter_annotations["phred_quality"],
                    new.letter_annotations["phred_quality"],
                )

    def check_sff(self, sff_name):
        self.check(
            sff_name, "sff", "Roche/E3MFGYR02_random_10_reads_no_trim.fasta", "fasta"
        )
        self.check(
            sff_name, "sff", "Roche/E3MFGYR02_random_10_reads_no_trim.qual", "qual"
        )
        self.check(
            sff_name, "sff-trim", "Roche/E3MFGYR02_random_10_reads.fasta", "fasta"
        )
        self.check(sff_name, "sff-trim", "Roche/E3MFGYR02_random_10_reads.qual", "qual")

    def test_original(self):
        """Test converting E3MFGYR02_random_10_reads.sff into FASTA+QUAL."""
        self.check_sff("Roche/E3MFGYR02_random_10_reads.sff")

    def test_no_manifest(self):
        """Test converting E3MFGYR02_no_manifest.sff into FASTA+QUAL."""
        self.check_sff("Roche/E3MFGYR02_no_manifest.sff")

    def test_alt_index_at_start(self):
        """Test converting E3MFGYR02_alt_index_at_start into FASTA+QUAL."""
        self.check_sff("Roche/E3MFGYR02_alt_index_at_start.sff")

    def test_alt_index_in_middle(self):
        """Test converting E3MFGYR02_alt_index_in_middle into FASTA+QUAL."""
        self.check_sff("Roche/E3MFGYR02_alt_index_in_middle.sff")

    def test_alt_index_at_end(self):
        """Test converting E3MFGYR02_alt_index_at_end into FASTA+QUAL."""
        self.check_sff("Roche/E3MFGYR02_alt_index_at_end.sff")

    def test_index_at_start(self):
        """Test converting E3MFGYR02_index_at_start into FASTA+QUAL."""
        self.check_sff("Roche/E3MFGYR02_index_at_start.sff")

    def test_index_at_end(self):
        """Test converting E3MFGYR02_index_in_middle into FASTA+QUAL."""
        self.check_sff("Roche/E3MFGYR02_index_in_middle.sff")


class TestReferenceFastqConversions(unittest.TestCase):
    """Tests where we have reference output."""

    def simple_check(self, base_name, in_variant):
        for out_variant in ["sanger", "solexa", "illumina"]:
            in_filename = "Quality/%s_original_%s.fastq" % (base_name, in_variant)
            self.assertTrue(os.path.isfile(in_filename))
            # Load the reference output...
            with open("Quality/%s_as_%s.fastq" % (base_name, out_variant)) as handle:
                expected = handle.read()

            with warnings.catch_warnings():
                if out_variant != "sanger":
                    # Ignore data loss warnings from max qualities
                    warnings.simplefilter("ignore", BiopythonWarning)
                # Check matches using convert...
                handle = StringIO()
                SeqIO.convert(
                    in_filename, "fastq-" + in_variant, handle, "fastq-" + out_variant
                )
                self.assertEqual(expected, handle.getvalue())
                # Check matches using parse/write
                handle = StringIO()
                SeqIO.write(
                    SeqIO.parse(in_filename, "fastq-" + in_variant),
                    handle,
                    "fastq-" + out_variant,
                )
                self.assertEqual(expected, handle.getvalue())

    def test_reference_conversion(self):
        tests = [
            ("illumina_full_range", "illumina"),
            ("sanger_full_range", "sanger"),
            ("longreads", "sanger"),
            ("solexa_full_range", "solexa"),
            ("misc_dna", "sanger"),
            ("wrapping", "sanger"),
            ("misc_rna", "sanger"),
        ]
        for base_name, variant in tests:
            assert variant in ["sanger", "solexa", "illumina"]
            self.simple_check(base_name, variant)


class TestQual(QualityIOTestBaseClass):
    """Tests with QUAL files."""

    def test_paired(self):
        """Check FASTQ parsing matches FASTA+QUAL parsing."""
        with open("Quality/example.fasta") as f, open("Quality/example.qual") as q:
            records1 = list(QualityIO.PairedFastaQualIterator(f, q))
        records2 = list(SeqIO.parse("Quality/example.fastq", "fastq"))
        self.compare_records(records1, records2)

    def test_qual(self):
        """Check FASTQ parsing matches QUAL parsing."""
        records1 = list(SeqIO.parse("Quality/example.qual", "qual"))
        records2 = list(SeqIO.parse("Quality/example.fastq", "fastq"))
        # Will ignore the unknown sequences :)
        self.compare_records(records1, records2)

    def test_qual_out(self):
        """Check FASTQ to QUAL output."""
        records = SeqIO.parse("Quality/example.fastq", "fastq")
        h = StringIO()
        SeqIO.write(records, h, "qual")
        with open("Quality/example.qual") as expected:
            self.assertEqual(h.getvalue(), expected.read())

    def test_fasta(self):
        """Check FASTQ parsing matches FASTA parsing."""
        records1 = list(SeqIO.parse("Quality/example.fasta", "fasta"))
        records2 = list(SeqIO.parse("Quality/example.fastq", "fastq"))
        self.compare_records(records1, records2)

    def test_fasta_out(self):
        """Check FASTQ to FASTA output."""
        records = SeqIO.parse("Quality/example.fastq", "fastq")
        h = StringIO()
        SeqIO.write(records, h, "fasta")
        with open("Quality/example.fasta") as expected:
            self.assertEqual(h.getvalue(), expected.read())

    def test_qual_negative(self):
        """Check QUAL negative scores mapped to PHRED zero."""
        data = """>1117_10_107_F3
23 31 -1 -1 -1 29 -1 -1 20 32 -1 18 25 7 -1 6 -1 -1 -1 30 -1 20 13 7 -1 -1 21 30 -1 24 -1 22 -1 -1 22 14 -1 12 26 21 -1 5 -1 -1 -1 20 -1 -1 12 28
>1117_10_146_F3
20 33 -1 -1 -1 29 -1 -1 28 28 -1 7 16 5 -1 30 -1 -1 -1 14 -1 4 13 4 -1 -1 11 13 -1 5 -1 7 -1 -1 10 16 -1 4 12 15 -1 8 -1 -1 -1 16 -1 -1 10 4
>1117_10_1017_F3
33 33 -1 -1 -1 27 -1 -1 17 16 -1 28 24 11 -1 6 -1 -1 -1 29 -1 8 29 24 -1 -1 8 8 -1 20 -1 13 -1 -1 8 13 -1 28 10 24 -1 10 -1 -1 -1 4 -1 -1 7 6
>1117_11_136_F3
16 22 -1 -1 -1 33 -1 -1 30 27 -1 27 28 32 -1 29 -1 -1 -1 27 -1 18 9 6 -1 -1 23 16 -1 26 -1 5 7 -1 22 7 -1 18 14 8 -1 8 -1 -1 -1 11 -1 -1 4 24"""  # noqa : W291
        h = StringIO(data)
        h2 = StringIO()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonParserWarning)
            self.assertEqual(4, SeqIO.convert(h, "qual", h2, "fastq"))
        self.assertEqual(
            h2.getvalue(),
            """\
@1117_10_107_F3
??????????????????????????????????????????????????
+
8@!!!>!!5A!3:(!'!!!?!5.(!!6?!9!7!!7/!-;6!&!!!5!!-=
@1117_10_146_F3
??????????????????????????????????????????????????
+
5B!!!>!!==!(1&!?!!!/!%.%!!,.!&!(!!+1!%-0!)!!!1!!+%
@1117_10_1017_F3
??????????????????????????????????????????????????
+
BB!!!<!!21!=9,!'!!!>!)>9!!))!5!.!!).!=+9!+!!!%!!('
@1117_11_136_F3
??????????????????????????????????????????????????
+
17!!!B!!?<!<=A!>!!!<!3*'!!81!;!&(!7(!3/)!)!!!,!!%9
""",
        )


class TestReadWrite(unittest.TestCase):
    """Test can read and write back files."""

    def test_fastq_2000(self):
        """Read and write back simple example with upper case 2000bp read."""
        data = "@%s\n%s\n+\n%s\n" % ("id descr goes here", "ACGT" * 500, "!@a~" * 500)
        handle = StringIO()
        self.assertEqual(
            1, SeqIO.write(SeqIO.parse(StringIO(data), "fastq"), handle, "fastq")
        )
        self.assertEqual(data, handle.getvalue())

    def test_fastq_1000(self):
        """Read and write back simple example with mixed case 1000bp read."""
        data = "@%s\n%s\n+\n%s\n" % (
            "id descr goes here",
            "ACGTNncgta" * 100,
            "abcd!!efgh" * 100,
        )
        handle = StringIO()
        self.assertEqual(
            1, SeqIO.write(SeqIO.parse(StringIO(data), "fastq"), handle, "fastq")
        )
        self.assertEqual(data, handle.getvalue())

    def test_fastq_dna(self):
        """Read and write back simple example with ambiguous DNA."""
        # First in upper case...
        data = "@%s\n%s\n+\n%s\n" % (
            "id descr goes here",
            ambiguous_dna_letters.upper(),
            "".join(chr(33 + q) for q in range(len(ambiguous_dna_letters))),
        )
        handle = StringIO()
        self.assertEqual(
            1, SeqIO.write(SeqIO.parse(StringIO(data), "fastq"), handle, "fastq")
        )
        self.assertEqual(data, handle.getvalue())
        # Now in lower case...
        data = "@%s\n%s\n+\n%s\n" % (
            "id descr goes here",
            ambiguous_dna_letters.lower(),
            "".join(chr(33 + q) for q in range(len(ambiguous_dna_letters))),
        )
        handle = StringIO()
        self.assertEqual(
            1, SeqIO.write(SeqIO.parse(StringIO(data), "fastq"), handle, "fastq")
        )
        self.assertEqual(data, handle.getvalue())

    def test_fastq_rna(self):
        """Read and write back simple example with ambiguous RNA."""
        # First in upper case...
        data = "@%s\n%s\n+\n%s\n" % (
            "id descr goes here",
            ambiguous_rna_letters.upper(),
            "".join(chr(33 + q) for q in range(len(ambiguous_rna_letters))),
        )
        handle = StringIO()
        self.assertEqual(
            1, SeqIO.write(SeqIO.parse(StringIO(data), "fastq"), handle, "fastq")
        )
        self.assertEqual(data, handle.getvalue())
        # Now in lower case...
        data = "@%s\n%s\n+\n%s\n" % (
            "id descr goes here",
            ambiguous_rna_letters.lower(),
            "".join(chr(33 + q) for q in range(len(ambiguous_rna_letters))),
        )
        handle = StringIO()
        self.assertEqual(
            1, SeqIO.write(SeqIO.parse(StringIO(data), "fastq"), handle, "fastq")
        )
        self.assertEqual(data, handle.getvalue())


class TestWriteRead(QualityIOTestBaseClass):
    """Test can write and read back files."""

    def write_read(self, filename, in_format, out_format):
        records = list(SeqIO.parse(filename, in_format))
        mode = self.get_mode(out_format)
        if mode == "b":
            handle = BytesIO()
        else:
            handle = StringIO()
        SeqIO.write(records, handle, out_format)
        handle.seek(0)
        # Now load it back and check it agrees,
        records2 = list(SeqIO.parse(handle, out_format))
        self.compare_records(records, records2, out_format)

    def test_generated(self):
        """Write and read back odd SeqRecord objects."""
        record1 = SeqRecord(
            Seq("ACGT" * 500),
            id="Test",
            description="Long " * 500,
            letter_annotations={"phred_quality": [40, 30, 20, 10] * 500},
        )
        record2 = SeqRecord(
            MutableSeq("NGGC" * 1000),
            id="Mut",
            description="very " * 1000 + "long",
            letter_annotations={"phred_quality": [0, 5, 5, 10] * 1000},
        )
        record3 = SeqRecord(
            UnknownSeq(2000, character="N"),
            id="Unk",
            description="l" + ("o" * 1000) + "ng",
            letter_annotations={"phred_quality": [0, 1] * 1000},
        )
        record4 = SeqRecord(
            Seq("ACGT" * 500),
            id="no_descr",
            description="",
            name="",
            letter_annotations={"phred_quality": [40, 50, 60, 62] * 500},
        )
        record5 = SeqRecord(
            Seq(""),
            id="empty_p",
            description="(could have been trimmed lots)",
            letter_annotations={"phred_quality": []},
        )
        record6 = SeqRecord(
            Seq(""),
            id="empty_s",
            description="(could have been trimmed lots)",
            letter_annotations={"solexa_quality": []},
        )
        record7 = SeqRecord(
            Seq("ACNN" * 500),
            id="Test_Sol",
            description="Long " * 500,
            letter_annotations={"solexa_quality": [40, 30, 0, -5] * 500},
        )
        record8 = SeqRecord(
            Seq("ACGT"),
            id="HighQual",
            description="With very large qualities that even Sanger FASTQ can't hold!",
            letter_annotations={"solexa_quality": [0, 10, 100, 1000]},
        )
        # TODO - Record with no identifier?
        records = [
            record1,
            record2,
            record3,
            record4,
            record5,
            record6,
            record7,
            record8,
        ]
        for fmt in ["fasta", "fastq", "fastq-solexa", "fastq-illumina", "qual"]:
            handle = StringIO()
            with warnings.catch_warnings():
                # TODO - Have a Biopython defined "DataLossWarning?"
                warnings.simplefilter("ignore", BiopythonWarning)
                SeqIO.write(records, handle, fmt)
            handle.seek(0)
            self.compare_records(records, list(SeqIO.parse(handle, fmt)), fmt)

    def check(self, filename, fmt, out_formats):
        for f in out_formats:
            self.write_read(filename, fmt, f)

    def test_tricky(self):
        """Write and read back tricky.fastq."""
        self.check(
            os.path.join("Quality", "tricky.fastq"),
            "fastq",
            [
                "fastq",
                "fastq-sanger",
                "fastq-illumina",
                "fastq-solexa",
                "fasta",
                "qual",
                "phd",
            ],
        )

    def test_sanger_93(self):
        """Write and read back sanger_93.fastq."""
        self.check(
            os.path.join("Quality", "sanger_93.fastq"),
            "fastq",
            ["fastq", "fastq-sanger", "fasta", "qual", "phd"],
        )
        with warnings.catch_warnings():
            # TODO - Have a Biopython defined "DataLossWarning?"
            warnings.simplefilter("ignore", BiopythonWarning)
            self.check(
                os.path.join("Quality", "sanger_93.fastq"),
                "fastq",
                ["fastq-solexa", "fastq-illumina"],
            )

    def test_sanger_faked(self):
        """Write and read back sanger_faked.fastq."""
        self.check(
            os.path.join("Quality", "sanger_faked.fastq"),
            "fastq",
            [
                "fastq",
                "fastq-sanger",
                "fastq-illumina",
                "fastq-solexa",
                "fasta",
                "qual",
                "phd",
            ],
        )

    def test_example_fasta(self):
        """Write and read back example.fasta."""
        self.write_read(os.path.join("Quality", "example.fasta"), "fasta", "fasta")
        # TODO - tests to check can't write FASTQ or QUAL...

    def test_example_fastq(self):
        """Write and read back example.fastq."""
        self.check(
            os.path.join("Quality", "example.fastq"),
            "fastq",
            [
                "fastq",
                "fastq-sanger",
                "fastq-illumina",
                "fastq-solexa",
                "fasta",
                "qual",
                "phd",
            ],
        )

    def test_example_qual(self):
        """Write and read back example.qual."""
        self.check(
            os.path.join("Quality", "example.qual"),
            "qual",
            [
                "fastq",
                "fastq-sanger",
                "fastq-illumina",
                "fastq-solexa",
                "fasta",
                "qual",
                "phd",
            ],
        )

    def test_solexa_faked(self):
        """Write and read back solexa_faked.fastq."""
        self.check(
            os.path.join("Quality", "solexa_faked.fastq"),
            "fastq-solexa",
            [
                "fastq",
                "fastq-sanger",
                "fastq-illumina",
                "fastq-solexa",
                "fasta",
                "qual",
                "phd",
            ],
        )

    def test_solexa_example(self):
        """Write and read back solexa_example.fastq."""
        self.check(
            os.path.join("Quality", "solexa_example.fastq"),
            "fastq-solexa",
            [
                "fastq",
                "fastq-sanger",
                "fastq-illumina",
                "fastq-solexa",
                "fasta",
                "qual",
                "phd",
            ],
        )

    def test_illumina_faked(self):
        """Write and read back illumina_faked.fastq."""
        self.check(
            os.path.join("Quality", "illumina_faked.fastq"),
            "fastq-illumina",
            [
                "fastq",
                "fastq-sanger",
                "fastq-illumina",
                "fastq-solexa",
                "fasta",
                "qual",
                "phd",
            ],
        )

    def test_greek_sff(self):
        """Write and read back greek.sff."""
        self.check(
            os.path.join("Roche", "greek.sff"),
            "sff",
            [
                "fastq",
                "fastq-sanger",
                "fastq-illumina",
                "fastq-solexa",
                "fasta",
                "qual",
                "phd",
                "sff",
            ],
        )

    def test_paired_sff(self):
        """Write and read back paired.sff."""
        self.check(
            os.path.join("Roche", "paired.sff"),
            "sff",
            [
                "fastq",
                "fastq-sanger",
                "fastq-illumina",
                "fastq-solexa",
                "fasta",
                "qual",
                "phd",
                "sff",
            ],
        )

    def test_E3MFGYR02(self):
        """Write and read back E3MFGYR02_random_10_reads.sff."""
        self.check(
            os.path.join("Roche", "E3MFGYR02_random_10_reads.sff"),
            "sff",
            [
                "fastq",
                "fastq-sanger",
                "fastq-illumina",
                "fastq-solexa",
                "fasta",
                "qual",
                "phd",
                "sff",
            ],
        )

    def test_E3MFGYR02_no_manifest(self):
        """Write and read back E3MFGYR02_no_manifest.sff."""
        self.check(
            os.path.join("Roche", "E3MFGYR02_no_manifest.sff"),
            "sff",
            [
                "fastq",
                "fastq-sanger",
                "fastq-illumina",
                "fastq-solexa",
                "fasta",
                "qual",
                "phd",
                "sff",
            ],
        )

    def test_E3MFGYR02_index_at_start(self):
        """Write and read back E3MFGYR02_index_at_start.sff."""
        self.check(
            os.path.join("Roche", "E3MFGYR02_index_at_start.sff"),
            "sff",
            [
                "fastq",
                "fastq-sanger",
                "fastq-illumina",
                "fastq-solexa",
                "fasta",
                "qual",
                "phd",
                "sff",
            ],
        )

    def test_E3MFGYR02_index_in_middle(self):
        """Write and read back E3MFGYR02_index_in_middle.sff."""
        self.check(
            os.path.join("Roche", "E3MFGYR02_index_in_middle.sff"),
            "sff",
            [
                "fastq",
                "fastq-sanger",
                "fastq-illumina",
                "fastq-solexa",
                "fasta",
                "qual",
                "phd",
                "sff",
            ],
        )

    def test_E3MFGYR02_alt_index_at_start(self):
        """Write and read back E3MFGYR02_alt_index_at_start.sff."""
        self.check(
            os.path.join("Roche", "E3MFGYR02_alt_index_at_start.sff"),
            "sff",
            [
                "fastq",
                "fastq-sanger",
                "fastq-illumina",
                "fastq-solexa",
                "fasta",
                "qual",
                "phd",
                "sff",
            ],
        )

    def test_E3MFGYR02_alt_index_in_middle(self):
        """Write and read back E3MFGYR02_alt_index_in_middle.sff."""
        self.check(
            os.path.join("Roche", "E3MFGYR02_alt_index_in_middle.sff"),
            "sff",
            [
                "fastq",
                "fastq-sanger",
                "fastq-illumina",
                "fastq-solexa",
                "fasta",
                "qual",
                "phd",
                "sff",
            ],
        )

    def test_E3MFGYR02_alt_index_at_end(self):
        """Write and read back E3MFGYR02_alt_index_at_end.sff."""
        self.check(
            os.path.join("Roche", "E3MFGYR02_alt_index_at_end.sff"),
            "sff",
            [
                "fastq",
                "fastq-sanger",
                "fastq-illumina",
                "fastq-solexa",
                "fasta",
                "qual",
                "phd",
                "sff",
            ],
        )

    def test_E3MFGYR02_trimmed(self):
        """Write and read back E3MFGYR02_random_10_reads.sff (trimmed)."""
        self.check(
            os.path.join("Roche", "E3MFGYR02_random_10_reads.sff"),
            "sff-trim",
            [
                "fastq",
                "fastq-sanger",
                "fastq-illumina",
                "fastq-solexa",
                "fasta",
                "qual",
                "phd",
            ],
        )  # not sff as output


class MappingTests(unittest.TestCase):
    """Quality mapping tests."""

    def test_solexa_quality_from_phred(self):
        """Mapping check for function solexa_quality_from_phred."""
        self.assertEqual(-5, round(QualityIO.solexa_quality_from_phred(0)))
        self.assertEqual(-5, round(QualityIO.solexa_quality_from_phred(1)))
        self.assertEqual(-2, round(QualityIO.solexa_quality_from_phred(2)))
        self.assertEqual(0, round(QualityIO.solexa_quality_from_phred(3)))
        self.assertEqual(2, round(QualityIO.solexa_quality_from_phred(4)))
        self.assertEqual(3, round(QualityIO.solexa_quality_from_phred(5)))
        self.assertEqual(5, round(QualityIO.solexa_quality_from_phred(6)))
        self.assertEqual(6, round(QualityIO.solexa_quality_from_phred(7)))
        self.assertEqual(7, round(QualityIO.solexa_quality_from_phred(8)))
        self.assertEqual(8, round(QualityIO.solexa_quality_from_phred(9)))
        for i in range(10, 100):
            self.assertEqual(i, round(QualityIO.solexa_quality_from_phred(i)))

    def test_phred_quality_from_solexa(self):
        """Mapping check for function phred_quality_from_solexa."""
        self.assertEqual(1, round(QualityIO.phred_quality_from_solexa(-5)))
        self.assertEqual(1, round(QualityIO.phred_quality_from_solexa(-4)))
        self.assertEqual(2, round(QualityIO.phred_quality_from_solexa(-3)))
        self.assertEqual(2, round(QualityIO.phred_quality_from_solexa(-2)))
        self.assertEqual(3, round(QualityIO.phred_quality_from_solexa(-1)))
        self.assertEqual(3, round(QualityIO.phred_quality_from_solexa(0)))
        self.assertEqual(4, round(QualityIO.phred_quality_from_solexa(1)))
        self.assertEqual(4, round(QualityIO.phred_quality_from_solexa(2)))
        self.assertEqual(5, round(QualityIO.phred_quality_from_solexa(3)))
        self.assertEqual(5, round(QualityIO.phred_quality_from_solexa(4)))
        self.assertEqual(6, round(QualityIO.phred_quality_from_solexa(5)))
        self.assertEqual(7, round(QualityIO.phred_quality_from_solexa(6)))
        self.assertEqual(8, round(QualityIO.phred_quality_from_solexa(7)))
        self.assertEqual(9, round(QualityIO.phred_quality_from_solexa(8)))
        self.assertEqual(10, round(QualityIO.phred_quality_from_solexa(9)))
        for i in range(10, 100):
            self.assertEqual(i, round(QualityIO.phred_quality_from_solexa(i)))

    def test_sanger_to_solexa(self):
        """Mapping check for FASTQ Sanger (0 to 93) to Solexa (-5 to 62)."""
        # The point of this test is the writing code doesn't actually use the
        # solexa_quality_from_phred function directly. For speed it uses a
        # cached dictionary of the mappings.
        seq = "N" * 94
        qual = "".join(chr(33 + q) for q in range(0, 94))
        expected_sol = [
            min(62, int(round(QualityIO.solexa_quality_from_phred(q))))
            for q in range(0, 94)
        ]
        in_handle = StringIO("@Test\n%s\n+\n%s" % (seq, qual))
        out_handle = StringIO()
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", BiopythonWarning)
            SeqIO.write(
                SeqIO.parse(in_handle, "fastq-sanger"), out_handle, "fastq-solexa"
            )
            self.assertLessEqual(len(w), 1, w)
        out_handle.seek(0)
        record = SeqIO.read(out_handle, "fastq-solexa")
        self.assertEqual(str(record.seq), seq)
        self.assertEqual(record.letter_annotations["solexa_quality"], expected_sol)

    def test_solexa_to_sanger(self):
        """Mapping check for FASTQ Solexa (-5 to 62) to Sanger (0 to 62)."""
        # The point of this test is the writing code doesn't actually use the
        # solexa_quality_from_phred function directly. For speed it uses a
        # cached dictionary of the mappings.
        seq = "N" * 68
        qual = "".join(chr(64 + q) for q in range(-5, 63))
        expected_phred = [
            round(QualityIO.phred_quality_from_solexa(q)) for q in range(-5, 63)
        ]
        in_handle = StringIO("@Test\n%s\n+\n%s" % (seq, qual))
        out_handle = StringIO()
        SeqIO.write(SeqIO.parse(in_handle, "fastq-solexa"), out_handle, "fastq-sanger")
        out_handle.seek(0)
        record = SeqIO.read(out_handle, "fastq-sanger")
        self.assertEqual(str(record.seq), seq)
        self.assertEqual(record.letter_annotations["phred_quality"], expected_phred)

    def test_sanger_to_illumina(self):
        """Mapping check for FASTQ Sanger (0 to 93) to Illumina (0 to 62)."""
        seq = "N" * 94
        qual = "".join(chr(33 + q) for q in range(0, 94))
        expected_phred = [min(62, q) for q in range(0, 94)]
        in_handle = StringIO("@Test\n%s\n+\n%s" % (seq, qual))
        out_handle = StringIO()
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", BiopythonWarning)
            SeqIO.write(
                SeqIO.parse(in_handle, "fastq-sanger"), out_handle, "fastq-illumina"
            )
            self.assertLessEqual(len(w), 1, w)
        out_handle.seek(0)
        record = SeqIO.read(out_handle, "fastq-illumina")
        self.assertEqual(str(record.seq), seq)
        self.assertEqual(record.letter_annotations["phred_quality"], expected_phred)

    def test_illumina_to_sanger(self):
        """Mapping check for FASTQ Illumina (0 to 62) to Sanger (0 to 62)."""
        seq = "N" * 63
        qual = "".join(chr(64 + q) for q in range(0, 63))
        expected_phred = list(range(63))
        in_handle = StringIO("@Test\n%s\n+\n%s" % (seq, qual))
        out_handle = StringIO()
        SeqIO.write(
            SeqIO.parse(in_handle, "fastq-illumina"), out_handle, "fastq-sanger"
        )
        out_handle.seek(0)
        record = SeqIO.read(out_handle, "fastq-sanger")
        self.assertEqual(str(record.seq), seq)
        self.assertEqual(record.letter_annotations["phred_quality"], expected_phred)


class TestSFF(unittest.TestCase):
    """Test SFF specific details."""

    def test_overlapping_clip(self):
        record = next(SeqIO.parse("Roche/greek.sff", "sff"))
        self.assertEqual(len(record), 395)
        s = str(record.seq.lower())
        # Apply overlapping clipping
        record.annotations["clip_qual_left"] = 51
        record.annotations["clip_qual_right"] = 44
        record.annotations["clip_adapter_left"] = 50
        record.annotations["clip_adapter_right"] = 75
        self.assertEqual(len(record), 395)
        self.assertEqual(len(record.seq), 395)
        # Save the clipped record...
        h = BytesIO()
        count = SeqIO.write(record, h, "sff")
        self.assertEqual(count, 1)
        # Now reload it...
        h.seek(0)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", BiopythonParserWarning)
            record = SeqIO.read(h, "sff")
            self.assertEqual(len(w), 1, w)
        self.assertEqual(record.annotations["clip_qual_left"], 51)
        self.assertEqual(record.annotations["clip_qual_right"], 44)
        self.assertEqual(record.annotations["clip_adapter_left"], 50)
        self.assertEqual(record.annotations["clip_adapter_right"], 75)
        self.assertEqual(len(record), 395)
        self.assertEqual(s, str(record.seq.lower()))
        # And check with trimming applied...
        h.seek(0)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", BiopythonParserWarning)
            record = SeqIO.read(h, "sff-trim")
            self.assertEqual(len(w), 1, w)
        self.assertEqual(len(record), 0)

    def test_negative_clip(self):
        for clip in [
            "clip_qual_left",
            "clip_qual_right",
            "clip_adapter_left",
            "clip_adapter_right",
        ]:
            record = next(SeqIO.parse("Roche/greek.sff", "sff"))
            self.assertEqual(len(record), 395)
            self.assertLessEqual(0, record.annotations[clip])
            record.annotations[clip] = -1
            with BytesIO() as h:
                self.assertRaises(ValueError, SeqIO.write, record, h, "sff")


class NonFastqTests(unittest.TestCase):
    def check_wrong_format(self, filename):
        for f in ("fastq", "fastq-sanger", "fastq-solexa", "fastq-illumina"):
            generator = SeqIO.parse(filename, f)
            self.assertRaises(ValueError, next, generator)

    def test_fasta_as_fastq(self):
        self.check_wrong_format("Fasta/elderberry.nu")

    def test_sff_as_fastq(self):
        self.check_wrong_format("Roche/greek.sff")


class TestsConverter(SeqIOConverterTestBaseClass, QualityIOTestBaseClass):
    def check_conversion(self, filename, in_format, out_format):
        msg = "Convert %s from %s to %s" % (filename, in_format, out_format)
        records = list(SeqIO.parse(filename, in_format))
        # Write it out...
        handle = StringIO()
        with warnings.catch_warnings():
            if out_format in (
                "fastq",
                "fastq-sanger",
                "fastq-solexa",
                "fastq-illumina",
            ):
                warnings.simplefilter("ignore", BiopythonWarning)
            SeqIO.write(records, handle, out_format)
        handle.seek(0)
        # Now load it back and check it agrees,
        records2 = list(SeqIO.parse(handle, out_format))
        self.assertEqual(len(records), len(records2), msg=msg)
        for record1, record2 in zip(records, records2):
            self.compare_record(record1, record2, out_format, msg=msg)
        # Finally, use the convert function, and check that agrees:
        handle2 = StringIO()
        with warnings.catch_warnings():
            if out_format in (
                "fastq",
                "fastq-sanger",
                "fastq-solexa",
                "fastq-illumina",
            ):
                warnings.simplefilter("ignore", BiopythonWarning)
            SeqIO.convert(filename, in_format, handle2, out_format)
        # We could re-parse this, but it is simpler and stricter:
        self.assertEqual(handle.getvalue(), handle2.getvalue(), msg=msg)

    def failure_check(self, filename, in_format, out_format):
        msg = "Confirm failure detection converting %s from %s to %s" % (
            filename,
            in_format,
            out_format,
        )
        # We want the SAME error message from parse/write as convert!
        with self.assertRaises(ValueError, msg=msg) as cm:
            records = list(SeqIO.parse(filename, in_format))
            self.write_records(records, out_format)
        err1 = str(cm.exception)
        # Now do the conversion...
        with self.assertRaises(ValueError, msg=msg) as cm:
            handle = StringIO()
            SeqIO.convert(filename, in_format, handle, out_format)
        err2 = str(cm.exception)
        # Verify that parse/write and convert give the same failure
        err_msg = "%s: parse/write and convert gave different failures" % msg
        self.assertEqual(err1, err2, msg=err_msg)

    def test_conversion(self):
        tests = [
            ("Quality/example.fastq", "fastq"),
            ("Quality/example.fastq", "fastq-sanger"),
            ("Quality/tricky.fastq", "fastq"),
            ("Quality/sanger_93.fastq", "fastq-sanger"),
            ("Quality/sanger_faked.fastq", "fastq-sanger"),
            ("Quality/solexa_faked.fastq", "fastq-solexa"),
            ("Quality/illumina_faked.fastq", "fastq-illumina"),
        ]
        for filename, fmt in tests:
            for (in_format, out_format) in self.formats:
                if in_format != fmt:
                    continue
                self.check_conversion(filename, in_format, out_format)

    def test_failure_detection(self):
        tests = [
            ("Quality/error_diff_ids.fastq", "fastq"),
            ("Quality/error_long_qual.fastq", "fastq"),
            ("Quality/error_no_qual.fastq", "fastq"),
            ("Quality/error_qual_del.fastq", "fastq"),
            ("Quality/error_qual_escape.fastq", "fastq"),
            ("Quality/error_qual_null.fastq", "fastq"),
            ("Quality/error_qual_space.fastq", "fastq"),
            ("Quality/error_qual_tab.fastq", "fastq"),
            ("Quality/error_qual_unit_sep.fastq", "fastq"),
            ("Quality/error_qual_vtab.fastq", "fastq"),
            ("Quality/error_short_qual.fastq", "fastq"),
            ("Quality/error_spaces.fastq", "fastq"),
            ("Quality/error_tabs.fastq", "fastq"),
            ("Quality/error_trunc_at_plus.fastq", "fastq"),
            ("Quality/error_trunc_at_qual.fastq", "fastq"),
            ("Quality/error_trunc_at_seq.fastq", "fastq"),
            ("Quality/error_trunc_in_title.fastq", "fastq"),
            ("Quality/error_trunc_in_seq.fastq", "fastq"),
            ("Quality/error_trunc_in_plus.fastq", "fastq"),
            ("Quality/error_trunc_in_qual.fastq", "fastq"),
            ("Quality/error_double_seq.fastq", "fastq"),
            ("Quality/error_double_qual.fastq", "fastq"),
        ]
        for filename, fmt in tests:
            for (in_format, out_format) in self.formats:
                if in_format != fmt:
                    continue
                if (
                    in_format
                    in ["fastq", "fastq-sanger", "fastq-solexa", "fastq-illumina"]
                    and out_format in ["fasta", "tab"]
                    and filename.startswith("Quality/error_qual_")
                ):
                    # TODO? These conversions don't check for bad characters in the quality,
                    # and in order to pass this strict test they should.
                    continue
                self.failure_check(filename, in_format, out_format)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
