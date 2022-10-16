# Copyright 2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Unit tests for Bio.SeqIO.convert(...) function."""
import unittest

from io import StringIO

from Bio import AlignIO


class ConvertTests(unittest.TestCase):
    def check_convert(self, in_filename, in_format, out_format, molecule_type):
        # Write it out using parse/write
        msg = f"Failed converting {in_filename} from {in_format} to {out_format}"
        handle = StringIO()
        aligns = list(AlignIO.parse(in_filename, in_format, None))
        if molecule_type:
            # Applying the molecule type explicitly:
            for align in aligns:
                for record in align:
                    record.annotations["molecule_type"] = molecule_type
        try:
            count = AlignIO.write(aligns, handle, out_format)
        except ValueError:
            count = 0
        # Write it out using convert passing filename and handle
        handle2 = StringIO()
        try:
            count2 = AlignIO.convert(
                in_filename, in_format, handle2, out_format, molecule_type
            )
        except ValueError:
            count2 = 0
        self.assertEqual(count, count2, msg=msg)
        self.assertEqual(handle.getvalue(), handle2.getvalue(), msg=msg)
        # Write it out using convert passing handle and handle
        handle2 = StringIO()
        try:
            with open(in_filename) as handle1:
                count2 = AlignIO.convert(
                    handle1, in_format, handle2, out_format, molecule_type
                )
        except ValueError:
            count2 = 0
        self.assertEqual(count, count2, msg=msg)
        self.assertEqual(handle.getvalue(), handle2.getvalue(), msg=msg)
        # TODO - convert passing an output filename?

    def test_convert(self):
        tests = [
            ("Clustalw/hedgehog.aln", "clustal", None),
            ("Nexus/test_Nexus_input.nex", "nexus", None),
            ("Stockholm/simple.sth", "stockholm", None),
            ("GFF/multi.fna", "fasta", "DNA"),
            ("Quality/example.fastq", "fastq", None),
            ("Quality/example.fastq", "fastq-sanger", "DNA"),
            ("Fasta/output001.m10", "fasta-m10", None),
            ("IntelliGenetics/VIF_mase-pro.txt", "ig", "protein"),
            ("NBRF/clustalw.pir", "pir", None),
        ]
        output_formats = ["fasta"] + sorted(AlignIO._FormatToWriter)
        for filename, in_format, mol_type in tests:
            for out_format in output_formats:
                self.check_convert(filename, in_format, out_format, mol_type)

    def test_clustal_to_nexus_without_mol_type(self):
        """Converting Clustal to NEXUS without a molecule type."""
        handle = StringIO()
        self.assertRaises(
            ValueError,
            AlignIO.convert,
            "Clustalw/protein.aln",
            "clustal",
            handle,
            "nexus",
        )

    def test_clustal_to_nexus_with_mol_type(self):
        """Converting Clustal to NEXUS with a molecule type."""
        handle = StringIO()
        self.assertEqual(
            1,
            AlignIO.convert(
                "Clustalw/protein.aln", "clustal", handle, "nexus", "protein"
            ),
        )
        self.assertIn(" datatype=protein ", handle.getvalue())


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
