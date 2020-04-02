# Copyright 2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for Bio.SeqIO.convert(...) function."""

import unittest
import warnings
from Bio import BiopythonWarning
from Bio.Seq import UnknownSeq
from Bio import SeqIO
from Bio.SeqIO import QualityIO
from Bio.SeqIO._convert import _converter as converter_dict
from io import StringIO
from Bio.Alphabet import generic_nucleotide, generic_dna


# TODO - share this with the QualityIO tests...
def truncation_expected(format):
    if format in ["fastq-solexa", "fastq-illumina"]:
        return 62
    elif format in ["fastq", "fastq-sanger"]:
        return 93
    else:
        return None


class ConvertTests(unittest.TestCase):

    # TODO - move compare_record to a shared test module...
    def compare_record(self, old, new, truncate, msg):
        """Quality aware SeqRecord comparison.

        This will check the mapping between Solexa and PHRED scores.
        It knows to ignore UnknownSeq objects for string matching (i.e. QUAL files).
        """
        self.assertEqual(old.id, new.id, msg=msg)
        self.assertTrue(
            old.description == new.description
            or (old.id + " " + old.description).strip() == new.description
            or new.description == "<unknown description>"
            or new.description == "",
            msg=msg
        )
        self.assertEqual(len(old.seq), len(new.seq), msg=msg)
        if isinstance(old.seq, UnknownSeq) or isinstance(new.seq, UnknownSeq):
            pass
        else:
            if len(old.seq) < 200:
                err_msg = "%s: '%s' vs '%s'" % (msg, old.seq, new.seq)
            else:
                err_msg = "%s: '%s...' vs '%s...'" % (msg, old.seq[:100], new.seq[:100])
            # Don't use assertEqual here, as it would show the complete
            # sequence in case of a failure.
            self.assertTrue(str(old.seq) == str(new.seq), msg=err_msg)

        for keyword in ("phred_quality", "solexa_quality"):
            q_old = old.letter_annotations.get(keyword)
            q_new = new.letter_annotations.get(keyword)
            if q_old is None or q_new is None:
                continue
            if truncate and q_old != q_new:
                q_old = [min(q, truncate) for q in q_old]
                q_new = [min(q, truncate) for q in q_new]
            err_msg = "%s: mismatch in %s" % (msg, keyword)
            self.assertEqual(q_old, q_new, msg=err_msg)

        q_old = old.letter_annotations.get("phred_quality")
        q_new = new.letter_annotations.get("solexa_quality")
        if q_old is not None and q_new is not None:
            # Mapping from Solexa to PHRED is lossy, but so is PHRED to Solexa.
            # Assume "old" is the original, and "new" has been converted.
            converted = [round(QualityIO.solexa_quality_from_phred(q)) for q in q_old]
            if truncate:
                converted = [min(q, truncate) for q in converted]
            err_msg = "%s: mismatch in phred_quality vs solexa_quality" % msg
            self.assertEqual(converted, q_new, msg=err_msg)

        q_old = old.letter_annotations.get("solexa_quality")
        q_new = new.letter_annotations.get("phred_quality")
        if q_old is not None and q_new is not None:
            # Mapping from Solexa to PHRED is lossy, but so is PHRED to Solexa.
            # Assume "old" is the original, and "new" has been converted.
            converted = [round(QualityIO.phred_quality_from_solexa(q)) for q in q_old]
            if truncate:
                converted = [min(q, truncate) for q in converted]
            err_msg = "%s: mismatch in solexa_quality vs phred_quality" % msg
            self.assertEqual(converted, q_new, msg=err_msg)

    def check_conversion(self, filename, in_format, out_format, alphabet):
        msg = "Convert %s from %s to %s" % (filename, in_format, out_format)
        records = list(SeqIO.parse(filename, in_format, alphabet))
        # Write it out...
        handle = StringIO()
        qual_truncate = truncation_expected(out_format)
        with warnings.catch_warnings():
            if qual_truncate:
                warnings.simplefilter("ignore", BiopythonWarning)
            SeqIO.write(records, handle, out_format)
        handle.seek(0)
        # Now load it back and check it agrees,
        records2 = list(SeqIO.parse(handle, out_format, alphabet))
        self.assertEqual(len(records), len(records2), msg=msg)
        for record1, record2 in zip(records, records2):
            self.compare_record(record1, record2, qual_truncate, msg=msg)
        # Finally, use the convert function, and check that agrees:
        handle2 = StringIO()
        with warnings.catch_warnings():
            if qual_truncate:
                warnings.simplefilter("ignore", BiopythonWarning)
            SeqIO.convert(filename, in_format, handle2, out_format, alphabet)
        # We could re-parse this, but it is simpler and stricter:
        self.assertEqual(handle.getvalue(), handle2.getvalue(), msg=msg)

    def failure_check(self, filename, in_format, out_format, alphabet):
        msg = "Confirm failure detection converting %s from %s to %s" % (filename, in_format, out_format)
        qual_truncate = truncation_expected(out_format)
        # We want the SAME error message from parse/write as convert!
        with self.assertRaises(ValueError, msg=msg) as cm:
            records = list(SeqIO.parse(filename, in_format, alphabet))
            self.write_records(records, out_format)
        err1 = str(cm.exception)
        # Now do the conversion...
        with self.assertRaises(ValueError, msg=msg) as cm:
            handle = StringIO()
            SeqIO.convert(filename, in_format, handle, out_format, alphabet)
        err2 = str(cm.exception)
        # Verify that parse/write and convert give the same failure
        err_msg = "%s: parse/write and convert gave different failures" % msg
        self.assertEqual(err1, err2, msg=err_msg)

    def test_conversion(self):
        tests = [
            ("Quality/example.fastq", "fastq", None),
            ("Quality/example.fastq", "fastq-sanger", generic_dna),
            ("Quality/tricky.fastq", "fastq", generic_nucleotide),
            ("Quality/sanger_93.fastq", "fastq-sanger", None),
            ("Quality/sanger_faked.fastq", "fastq-sanger", generic_dna),
            ("Quality/solexa_faked.fastq", "fastq-solexa", generic_dna),
            ("Quality/illumina_faked.fastq", "fastq-illumina", generic_dna),
            ("EMBL/U87107.embl", "embl", None),
            ("EMBL/TRBG361.embl", "embl", None),
            ("GenBank/NC_005816.gb", "gb", None),
            ("GenBank/cor6_6.gb", "genbank", None),
        ]
        for filename, fmt, alphabet in tests:
            for (in_format, out_format) in converter_dict:
                if in_format != fmt:
                    continue
                self.check_conversion(filename, in_format, out_format, alphabet)

    def test_failure_detection(self):
        tests = [
            ("Quality/error_diff_ids.fastq", "fastq", None),
            ("Quality/error_long_qual.fastq", "fastq", None),
            ("Quality/error_no_qual.fastq", "fastq", None),
            ("Quality/error_qual_del.fastq", "fastq", None),
            ("Quality/error_qual_escape.fastq", "fastq", None),
            ("Quality/error_qual_null.fastq", "fastq", None),
            ("Quality/error_qual_space.fastq", "fastq", None),
            ("Quality/error_qual_tab.fastq", "fastq", None),
            ("Quality/error_qual_unit_sep.fastq", "fastq", None),
            ("Quality/error_qual_vtab.fastq", "fastq", None),
            ("Quality/error_short_qual.fastq", "fastq", None),
            ("Quality/error_spaces.fastq", "fastq", None),
            ("Quality/error_tabs.fastq", "fastq", None),
            ("Quality/error_trunc_at_plus.fastq", "fastq", None),
            ("Quality/error_trunc_at_qual.fastq", "fastq", None),
            ("Quality/error_trunc_at_seq.fastq", "fastq", None),
            ("Quality/error_trunc_in_title.fastq", "fastq", generic_dna),
            ("Quality/error_trunc_in_seq.fastq", "fastq", generic_nucleotide),
            ("Quality/error_trunc_in_plus.fastq", "fastq", None),
            ("Quality/error_trunc_in_qual.fastq", "fastq", generic_dna),
            ("Quality/error_double_seq.fastq", "fastq", generic_dna),
            ("Quality/error_double_qual.fastq", "fastq", generic_dna),
        ]
        for filename, fmt, alphabet in tests:
            for (in_format, out_format) in converter_dict:
                if in_format != fmt:
                    continue
                if (
                    in_format in ["fastq", "fastq-sanger", "fastq-solexa", "fastq-illumina"]
                    and out_format in ["fasta", "tab"]
                    and filename.startswith("Quality/error_qual_")
                ):
                    # TODO? These conversions don't check for bad characters in the quality,
                    # and in order to pass this strict test they should.
                    continue
                self.failure_check(filename, in_format, out_format, alphabet)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
