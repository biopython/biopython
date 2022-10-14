# Copyright 2007-2010 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for SeqIO write module."""
import os
import unittest
import warnings

from io import BytesIO
from io import StringIO

from Bio import AlignIO
from Bio import BiopythonWarning
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from test_SeqIO import SeqIOTestBaseClass


# List of formats including alignment only file formats we can read AND write.
# We don't care about the order
test_write_read_alignment_formats = sorted(SeqIO._FormatToWriter)
for fmt in sorted(AlignIO._FormatToWriter):
    if fmt not in test_write_read_alignment_formats:
        test_write_read_alignment_formats.append(fmt)
test_write_read_alignment_formats.remove("gb")  # an alias for genbank
test_write_read_alignment_formats.remove("fastq-sanger")  # an alias for fastq


# This is a list of three-tuples.  Each tuple contains a
# list of SeqRecord objects, a description (string), and
# a list of tuples for expected failures (each with a
# list of formats, exception type, exception message).
test_records = [
    ([], "zero records", {}),
    (
        [
            SeqRecord(
                Seq("CHSMAIKLSSEHNIPSGIANAL"),
                id="Alpha",
                annotations={"molecule_type": "protein"},
            ),
            SeqRecord(
                Seq("HNGFTALEGEIHHLTHGEKVAF"),
                id="Gamma",
                annotations={"molecule_type": "protein"},
            ),
            SeqRecord(
                Seq("DITHGVG"), id="delta", annotations={"molecule_type": "protein"}
            ),
        ],
        "three peptides of different lengths",
        [],
    ),
    (
        [
            SeqRecord(
                Seq("CHSMAIKLSSEHNIPSGIANAL"),
                id="Alpha",
                annotations={"molecule_type": "protein"},
            ),
            SeqRecord(
                Seq("VHGMAHPLGAFYNTPHGVANAI"),
                id="Beta",
                annotations={"molecule_type": "protein"},
            ),
            SeqRecord(
                Seq("HNGFTALEGEIHHLTHGEKVAF"),
                id="Gamma",
                annotations={"molecule_type": "protein"},
            ),
        ],
        "three proteins alignment",
        [],
    ),
    (
        [
            SeqRecord(
                Seq("AATAAACCTTGCTGGCCATTGTGATCCATCCA"),
                id="X",
                annotations={"molecule_type": "DNA"},
            ),
            SeqRecord(
                Seq("ACTCAACCTTGCTGGTCATTGTGACCCCAGCA"),
                id="Y",
                annotations={"molecule_type": "DNA"},
            ),
            SeqRecord(
                Seq("TTTCCTCGGAGGCCAATCTGGATCAAGACCAT"),
                id="Z",
                annotations={"molecule_type": "DNA"},
            ),
        ],
        "three DNA sequence alignment",
        [],
    ),
    (
        [
            SeqRecord(
                Seq("AATAAACCTTGCTGGCCATTGTGATCCATCCA"),
                id="X",
                name="The\nMystery\rSequece:\r\nX",
                annotations={"molecule_type": "DNA"},
            ),
            SeqRecord(
                Seq("ACTCAACCTTGCTGGTCATTGTGACCCCAGCA"),
                id="Y",
                description=f"an{os.linesep}evil\rdescription right\nhere",
                annotations={"molecule_type": "DNA"},
            ),
            SeqRecord(
                Seq("TTTCCTCGGAGGCCAATCTGGATCAAGACCAT"),
                id="Z",
                annotations={"molecule_type": "DNA"},
            ),
        ],
        "3 DNA seq alignment with CR/LF in name/descr",
        [
            (
                ["genbank"],
                ValueError,
                r"Invalid whitespace in 'The\nMystery\rSequece:\r\nX' for LOCUS line",
            )
        ],
    ),
    (
        [
            SeqRecord(
                Seq("CHSMAIKLSSEHNIPSGIANAL"),
                id="Alpha",
                annotations={"molecule_type": "protein"},
            ),
            SeqRecord(
                Seq("VHGMAHPLGAFYNTPHGVANAI"),
                id="Beta",
                annotations={"molecule_type": "protein"},
            ),
            SeqRecord(
                Seq("VHGMAHPLGAFYNTPHGVANAI"),
                id="Beta",
                annotations={"molecule_type": "protein"},
            ),
            SeqRecord(
                Seq("HNGFTALEGEIHHLTHGEKVAF"),
                id="Gamma",
                annotations={"molecule_type": "protein"},
            ),
        ],
        "alignment with repeated record",
        [
            (["stockholm"], ValueError, "Duplicate record identifier: Beta"),
            (
                ["phylip", "phylip-relaxed", "phylip-sequential"],
                ValueError,
                "Repeated name 'Beta' (originally 'Beta'), possibly due to truncation",
            ),
        ],
    ),
]
# Meddle with the annotation too:
assert test_records[4][1] == "3 DNA seq alignment with CR/LF in name/descr"
# Add a list of strings,
test_records[4][0][2].annotations["note"] = [
    f"Note{os.linesep}also\r\nhas\n evil line\rbreaks!",
    "Wow",
]
# Add a simple string
test_records[4][0][2].annotations[
    "comment"
] = f"More{os.linesep}of\r\nthese\n evil line\rbreaks!"
# Add a float too:
test_records[4][0][2].annotations["weight"] = 2.5


class WriterTests(SeqIOTestBaseClass):
    def check(self, records, fmt, descr):
        """General test function with with a little format specific information.

        This has some general expected exceptions hard coded!
        """
        # TODO - Check the exception messages?
        lengths = len({len(r) for r in records})
        dna = all(set(record.seq.upper()).issubset("ACGTN") for record in records)
        if not records and fmt in [
            "stockholm",
            "phylip",
            "phylip-relaxed",
            "phylip-sequential",
            "nexus",
            "clustal",
            "sff",
            "mauve",
        ]:
            self.check_write_fails(
                records, fmt, descr, ValueError, "Must have at least one sequence"
            )
        elif not records and fmt in ["nib", "xdna"]:
            self.check_write_fails(
                records, fmt, descr, ValueError, "Must have one sequence"
            )
        elif lengths > 1 and fmt in AlignIO._FormatToWriter:
            self.check_write_fails(
                records, fmt, descr, ValueError, "Sequences must all be the same length"
            )
        elif (not dna) and fmt == "nib":
            self.check_write_fails(
                records,
                fmt,
                descr,
                ValueError,
                "Sequence should contain A,C,G,T,N,a,c,g,t,n only",
            )
        elif len(records) > 1 and fmt in ["nib", "xdna"]:
            self.check_write_fails(
                records, fmt, descr, ValueError, "More than one sequence found"
            )
        elif records and fmt in [
            "fastq",
            "fastq-sanger",
            "fastq-solexa",
            "fastq-illumina",
            "qual",
            "phd",
        ]:
            self.check_write_fails(
                records,
                fmt,
                descr,
                ValueError,
                "No suitable quality scores found in "
                "letter_annotations of SeqRecord "
                "(id=%s)." % records[0].id,
            )
        elif records and fmt == "sff":
            self.check_write_fails(
                records, fmt, descr, ValueError, "Missing SFF flow information"
            )
        else:
            self.check_simple(records, fmt, descr)

    def check_simple(self, records, fmt, descr):
        msg = f"Test failure {fmt} for {descr}"
        mode = self.get_mode(fmt)
        if mode == "t":
            handle = StringIO()
        elif mode == "b":
            handle = BytesIO()
        count = SeqIO.write(records, handle, fmt)
        self.assertEqual(count, len(records), msg=msg)
        # Now read them back...
        handle.seek(0)
        new_records = list(SeqIO.parse(handle, fmt))
        self.assertEqual(len(new_records), len(records), msg=msg)
        for record, new_record in zip(records, new_records):
            # Using compare_record(record, new_record) is too strict
            if fmt == "nexus":
                # The nexus parser will dis-ambiguate repeated record ids.
                self.assertTrue(
                    record.id == new_record.id
                    or new_record.id.startswith(record.id + ".copy"),
                    msg=msg,
                )
            else:
                self.assertEqual(record.id, new_record.id, msg=msg)
            self.assertEqual(record.seq, new_record.seq, msg=msg)
        handle.close()

    def check_write_fails(self, records, fmt, descr, err_type, err_msg=""):
        msg = f"Test failure {fmt} for {descr}"
        mode = self.get_mode(fmt)
        if mode == "t":
            handle = StringIO()
        elif mode == "b":
            handle = BytesIO()
        if err_msg:
            with self.assertRaises(err_type, msg=msg) as cm:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", BiopythonWarning)
                    SeqIO.write(records, handle, fmt)
            self.assertEqual(str(cm.exception), err_msg, msg=msg)
        else:
            with self.assertRaises(err_type, msg=msg) as cm:
                SeqIO.write(records, handle, fmt)
        handle.close()

    def test_bad_handle(self):
        handle = os.devnull
        record = SeqRecord(Seq("CHSMAIKLSSEHNIPSGIANAL"), id="Alpha")
        records = [record]
        fmt = "fasta"
        # These deliberately mix up the handle and record order:
        self.assertRaises(TypeError, SeqIO.write, handle, record, fmt)
        self.assertRaises(TypeError, SeqIO.write, handle, records, fmt)
        self.assertEqual(1, SeqIO.write(records, handle, fmt))

    def test_alignment_formats(self):
        for (records, descr, errs) in test_records:
            for fmt in test_write_read_alignment_formats:
                for err_formats, err_type, err_msg in errs:
                    if fmt in err_formats:
                        self.check_write_fails(records, fmt, descr, err_type, err_msg)
                        break
                else:
                    self.check(records, fmt, descr)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
