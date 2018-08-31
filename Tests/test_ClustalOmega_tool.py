# Copyright 2008-2011 by Peter Cock.  All rights reserved.
# Revisions copyright 2012 by Christian Brueffer.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio import MissingExternalDependencyError

import sys
import os
import unittest
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Application import ApplicationError

#################################################################

# Try to avoid problems when the OS is in another language
os.environ['LANG'] = 'C'

clustalo_exe = None
from Bio._py3k import getoutput
try:
    output = getoutput("clustalo --help")
    if output.startswith("Clustal Omega"):
        clustalo_exe = "clustalo"
except OSError:
    # TODO: Use FileNotFoundError once we drop Python 2
    pass

if not clustalo_exe:
    raise MissingExternalDependencyError(
        "Install clustalo if you want to use Clustal Omega from Biopython.")


class ClustalOmegaTestCase(unittest.TestCase):

    def setUp(self):
        self.files_to_clean = set()

    def tearDown(self):
        for filename in self.files_to_clean:
            if os.path.isfile(filename):
                os.remove(filename)

    def standard_test_procedure(self, cline):
        """Standard testing procedure used by all tests."""

        # Overwrite existing files.
        cline.force = True

        # Mark output files for later cleanup.
        self.add_file_to_clean(cline.outfile)
        if cline.guidetree_out:
            self.add_file_to_clean(cline.guidetree_out)

        input_records = SeqIO.to_dict(SeqIO.parse(cline.infile, "fasta"))
        self.assertEqual(str(eval(repr(cline))), str(cline))
        output, error = cline()
        self.assertTrue(not output or output.strip().startswith("CLUSTAL"))

        # Test if ClustalOmega executed successfully.
        self.assertTrue(error.strip() == "" or
                        error.startswith("WARNING: Sequence type is DNA.") or
                        error.startswith("WARNING: DNA alignment is still experimental."))

        # Check the output...
        align = AlignIO.read(cline.outfile, "clustal")
        output_records = SeqIO.to_dict(SeqIO.parse(cline.outfile, "clustal"))
        self.assertEqual(len(set(input_records.keys())), len(set(output_records.keys())))
        for record in align:
            self.assertEqual(str(record.seq), str(output_records[record.id].seq))

        # TODO - Try and parse this with Bio.Nexus?
        if cline.guidetree_out:
            self.assertTrue(os.path.isfile(cline.guidetree_out))

    def add_file_to_clean(self, filename):
        """Adds a file for deferred removal by the tearDown routine."""
        self.files_to_clean.add(filename)

#################################################################


class ClustalOmegaTestErrorConditions(ClustalOmegaTestCase):

    def test_empty_file(self):
        """Test an empty file."""
        input_file = "does_not_exist.fasta"
        self.assertFalse(os.path.isfile(input_file))
        cline = ClustalOmegaCommandline(clustalo_exe, infile=input_file)
        try:
            stdout, stderr = cline()
        except ApplicationError as err:
            self.assertTrue("Cannot open sequence file" in str(err) or
                            "Cannot open input file" in str(err) or
                            "Non-zero return code" in str(err), str(err))
        else:
            self.fail("Should have failed, returned:\n%s\n%s" % (stdout, stderr))

    def test_single_sequence(self):
        """Test an input file containing a single sequence."""
        input_file = "Fasta/f001"
        self.assertTrue(os.path.isfile(input_file))
        self.assertEqual(len(list(SeqIO.parse(input_file, "fasta"))), 1)
        cline = ClustalOmegaCommandline(clustalo_exe, infile=input_file)
        try:
            stdout, stderr = cline()
        except ApplicationError as err:
            self.assertIn("contains 1 sequence, nothing to align", str(err))
        else:
            self.fail("Should have failed, returned:\n%s\n%s" % (stdout, stderr))

    def test_invalid_format(self):
        """Test an input file in an invalid format."""
        input_file = "Medline/pubmed_result1.txt"
        self.assertTrue(os.path.isfile(input_file))
        cline = ClustalOmegaCommandline(clustalo_exe, infile=input_file)
        try:
            stdout, stderr = cline()
        except ApplicationError as err:
            # Ideally we'd catch the return code and raise the specific
            # error for "invalid format".
            self.assertIn("Can't determine format of sequence file", str(err))
        else:
            self.fail("Should have failed, returned:\n%s\n%s" % (stdout, stderr))

#################################################################


class ClustalOmegaTestNormalConditions(ClustalOmegaTestCase):

    def test_simple_fasta(self):
        """Test a simple fasta file."""
        input_file = "Registry/seqs.fasta"
        output_file = "temp_test.aln"

        cline = ClustalOmegaCommandline(clustalo_exe,
                                        infile=input_file,
                                        outfile=output_file,
                                        outfmt="clustal")

        self.standard_test_procedure(cline)

    def test_properties(self):
        """Test setting options via properties."""
        input_file = "Registry/seqs.fasta"
        output_file = "temp_test.aln"

        cline = ClustalOmegaCommandline(clustalo_exe)
        cline.infile = input_file
        cline.outfile = output_file
        cline.outfmt = "clustal"

        self.standard_test_procedure(cline)

    def test_input_filename_with_space(self):
        """Test an input filename containing a space."""
        input_file = "Clustalw/temp horses.fasta"
        handle = open(input_file, "w")
        SeqIO.write(SeqIO.parse("Phylip/hennigian.phy", "phylip"), handle, "fasta")
        handle.close()
        output_file = "temp_test.aln"

        cline = ClustalOmegaCommandline(clustalo_exe,
                                        infile=input_file,
                                        outfile=output_file,
                                        outfmt="clustal")

        self.add_file_to_clean(input_file)
        self.standard_test_procedure(cline)

    def test_output_filename_with_spaces(self):
        """Test an output filename containing spaces."""
        input_file = "Registry/seqs.fasta"
        output_file = "temp with spaces.aln"

        cline = ClustalOmegaCommandline(clustalo_exe,
                                        infile=input_file,
                                        outfile=output_file,
                                        outfmt="clustal")

        self.standard_test_procedure(cline)

    def test_large_fasta_file(self):
        """Test a large fasta input file."""
        # Create a large input file by converting another example file
        # (See Bug 2804, this will produce so much output on stdout that
        # subprocess could suffer a deadlock and hang).  Using all the
        # records should show the deadlock but is very slow - just thirty
        # seems to lockup on Mac OS X, even 20 on Linux (without the fix).
        input_file = "temp_cw_prot.fasta"
        handle = open(input_file, "w")
        records = list(SeqIO.parse("NBRF/Cw_prot.pir", "pir"))[:40]
        SeqIO.write(records, handle, "fasta")
        handle.close()
        del handle, records
        output_file = "temp_cw_prot.aln"

        cline = ClustalOmegaCommandline(clustalo_exe,
                                        infile=input_file,
                                        outfile=output_file,
                                        outfmt="clustal")

        self.add_file_to_clean(input_file)
        self.standard_test_procedure(cline)

    def test_newtree_files(self):
        """Test requesting a guide tree."""
        input_file = "Fasta/f002"
        output_file = "temp_test.aln"
        newtree_file = "temp_test.dnd"

        cline = ClustalOmegaCommandline(clustalo_exe,
                                        infile=input_file,
                                        outfile=output_file,
                                        guidetree_out=newtree_file,
                                        outfmt="clustal")

        self.standard_test_procedure(cline)
        cline.guidetree_out = "temp with space.dnd"
        self.standard_test_procedure(cline)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
