# Copyright 2013 by Christian Brueffer.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
import sys
import unittest
from Bio import AlignIO
from Bio import MissingExternalDependencyError
from Bio import SeqIO
from Bio.Align.Applications import MSAProbsCommandline
from Bio.Application import ApplicationError
from Bio._py3k import getoutput

#################################################################

#Try to avoid problems when the OS is in another language
os.environ['LANG'] = 'C'

msaprobs_exe = None
if sys.platform == "win32":
    #TODO
    raise MissingExternalDependencyError("Testing this on Windows is not implemented yet")
else:
    output = getoutput("msaprobs -version")
    if output.startswith("MSAPROBS version"):
        msaprobs_exe = "msaprobs"

if not msaprobs_exe:
    raise MissingExternalDependencyError(
        "Install msaprobs if you want to use MSAProbs from Biopython.")


class MSAProbsTestCase(unittest.TestCase):

    def setUp(self):
        self.files_to_clean = set()

    def tearDown(self):
        for filename in self.files_to_clean:
            if os.path.isfile(filename):
                os.remove(filename)

    def standard_test_procedure(self, cline):
        """Standard testing procedure used by all tests."""

        # Mark output files for later cleanup.
        self.add_file_to_clean(cline.outfile)

        input_records = SeqIO.to_dict(SeqIO.parse(cline.infile, "fasta"))
        self.assertEqual(str(eval(repr(cline))), str(cline))
        output, error = cline()

    def add_file_to_clean(self, filename):
        """Adds a file for deferred removal by the tearDown routine."""
        self.files_to_clean.add(filename)

#################################################################


class MSAProbsTestErrorConditions(MSAProbsTestCase):

    def test_empty_file(self):
        """Test an empty file."""
        input_file = "does_not_exist.fasta"
        self.assertFalse(os.path.isfile(input_file))
        cline = MSAProbsCommandline(msaprobs_exe, infile=input_file)
        try:
            stdout, stderr = cline()
        except ApplicationError as err:
            self.assertTrue("Cannot open sequence file" in str(err) or
                            "Cannot open input file" in str(err) or
                            "non-zero exit status" in str(err))
        else:
            self.fail("Should have failed, returned:\n%s\n%s" % (stdout, stderr))

    def test_single_sequence(self):
        """Test an input file containing a single sequence."""
        input_file = "Fasta/f001"
        self.assertTrue(os.path.isfile(input_file))
        self.assertEqual(len(list(SeqIO.parse(input_file, "fasta"))), 1)
        cline = MSAProbsCommandline(msaprobs_exe, infile=input_file)
        try:
            stdout, stderr = cline()
        except ApplicationError as err:
            self.assertEqual(err.returncode, 139)
        else:
            self.fail("Should have failed, returned:\n%s\n%s" % (stdout, stderr))

    def test_invalid_format(self):
        """Test an input file in an invalid format."""
        input_file = "Medline/pubmed_result1.txt"
        self.assertTrue(os.path.isfile(input_file))
        cline = MSAProbsCommandline(msaprobs_exe, infile=input_file)
        try:
            stdout, stderr = cline()
        except ApplicationError as err:
            self.assertEqual(err.returncode, 1)
        else:
            self.fail("Should have failed, returned:\n%s\n%s" % (stdout, stderr))

#################################################################


class MSAProbsTestNormalConditions(MSAProbsTestCase):

    def test_simple_fasta(self):
        """Test a simple fasta file."""
        input_file = "Registry/seqs.fasta"
        output_file = "temp_test.aln"

        cline = MSAProbsCommandline(msaprobs_exe,
                                        infile=input_file,
                                        outfile=output_file,
                                        clustalw=True)

        self.standard_test_procedure(cline)

    def test_properties(self):
        """Test setting options via properties."""
        input_file = "Registry/seqs.fasta"
        output_file = "temp_test.aln"

        cline = MSAProbsCommandline(msaprobs_exe)
        cline.infile = input_file
        cline.outfile = output_file
        cline.clustalw = True

        self.standard_test_procedure(cline)

    def test_input_filename_with_space(self):
        """Test an input filename containing a space."""
        input_file = "Clustalw/temp horses.fasta"
        handle = open(input_file, "w")
        SeqIO.write(SeqIO.parse("Phylip/hennigian.phy", "phylip"), handle, "fasta")
        handle.close()
        output_file = "temp_test.aln"

        cline = MSAProbsCommandline(msaprobs_exe,
                                        infile=input_file,
                                        outfile=output_file,
                                        clustalw=True)

        self.add_file_to_clean(input_file)
        self.standard_test_procedure(cline)

    def test_output_filename_with_spaces(self):
        """Test an output filename containing spaces."""
        input_file = "Registry/seqs.fasta"
        output_file = "temp with spaces.aln"

        cline = MSAProbsCommandline(msaprobs_exe,
                                        infile=input_file,
                                        outfile=output_file,
                                        clustalw=True)

        self.standard_test_procedure(cline)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
