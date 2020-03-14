# Copyright 2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for Bio.SeqIO.convert(...) function."""
import unittest
from io import StringIO

from Bio import AlignIO
from Bio.Alphabet import generic_protein, generic_nucleotide, generic_dna


class ConvertTests(unittest.TestCase):

    def check_convert(self, in_filename, in_format, out_format, alphabet=None):
        # Write it out using parse/write
        msg = "Failed converting %s from %s to %s" % (in_filename, in_format, out_format)
        handle = StringIO()
        aligns = list(AlignIO.parse(in_filename, in_format, None, alphabet))
        try:
            count = AlignIO.write(aligns, handle, out_format)
        except ValueError:
            count = 0
        # Write it out using convert passing filename and handle
        handle2 = StringIO()
        try:
            count2 = AlignIO.convert(in_filename, in_format, handle2, out_format, alphabet)
        except ValueError:
            count2 = 0
        self.assertEqual(count, count2, msg=msg)
        self.assertEqual(handle.getvalue(), handle2.getvalue(), msg=msg)
        # Write it out using convert passing handle and handle
        handle2 = StringIO()
        try:
            with open(in_filename) as handle1:
                count2 = AlignIO.convert(handle1, in_format, handle2, out_format, alphabet)
        except ValueError:
            count2 = 0
        self.assertEqual(count, count2, msg=msg)
        self.assertEqual(handle.getvalue(), handle2.getvalue(), msg=msg)
        # TODO - convert passing an output filename?

    def test_convert(self):
        tests = [("Clustalw/hedgehog.aln", "clustal", None),
                 ("Nexus/test_Nexus_input.nex", "nexus", None),
                 ("Stockholm/simple.sth", "stockholm", None),
                 ("GFF/multi.fna", "fasta", generic_nucleotide),
                 ("Quality/example.fastq", "fastq", None),
                 ("Quality/example.fastq", "fastq-sanger", generic_dna),
                 ("Fasta/output001.m10", "fasta-m10", None),
                 ("IntelliGenetics/VIF_mase-pro.txt", "ig", generic_protein),
                 ("NBRF/clustalw.pir", "pir", None),
                 ]
        output_formats = ["fasta"] + sorted(AlignIO._FormatToWriter)
        for filename, in_format, alphabet in tests:
            for out_format in output_formats:
                self.check_convert(filename, in_format, out_format, alphabet)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
