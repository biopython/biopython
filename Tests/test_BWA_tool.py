# Copyright 2013 by Saket Choudhary.  Based on test_Clustalw_tool.py by Peter
# Cock .
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for calling BWA."""

from Bio import MissingExternalDependencyError
import sys
import os
import unittest

# TODO from Bio.Sequencing.Applications import BwaBwaswCommandline
from Bio.Sequencing.Applications import BwaIndexCommandline
from Bio.Sequencing.Applications import BwaAlignCommandline
from Bio.Sequencing.Applications import BwaSamseCommandline
from Bio.Sequencing.Applications import BwaSampeCommandline
from Bio.Sequencing.Applications import BwaMemCommandline


#################################################################

# Try to avoid problems when the OS is in another language
os.environ["LANG"] = "C"

bwa_exe = None
if sys.platform == "win32":
    # TODO - Check the path?
    try:
        # This can vary depending on the Windows language.
        prog_files = os.environ["PROGRAMFILES"]
    except KeyError:
        prog_files = r"C:\Program Files"
    likely_dirs = ["bwa", "bwa-0.6.2", ""]
    likely_exes = ["bwa"]
    for folder in likely_dirs:
        if os.path.isdir(os.path.join(prog_files, folder)):
            for filename in likely_exes:
                if os.path.isfile(os.path.join(prog_files, folder, filename)):
                    bwa_exe = os.path.join(prog_files, folder, filename)
                    break
            if bwa_exe:
                break
else:
    from subprocess import getoutput

    output = getoutput("bwa")

    # Since "not found" may be in another language, try and be sure this is
    # really the bwa tool's output
    bwa_found = False
    if "not found" not in output and "not recognized" not in output:
        if "bwa" in output and "alignment via Burrows-Wheeler transformation" in output:
            bwa_exe = "bwa"
    skip_aln_tests = False
    aln_output = getoutput("bwa aln")
    if "unrecognized" in aln_output:
        skip_aln_tests = True
        print("'bwa aln' is unrecognized, skipping aln/samse/sampe tests")

if not bwa_exe:
    raise MissingExternalDependencyError(
        "Install bwa and correctly set"
        " the file path to the program if"
        " you want to use it from Biopython"
    )


class BwaTestCase(unittest.TestCase):
    """Class for implementing BWA test cases."""

    def setUp(self):
        self.reference_file = "BWA/human_g1k_v37_truncated.fasta"
        self.reference_extensions = ["amb", "ann", "bwt", "pac", "sa"]
        self.infile1 = "BWA/HNSCC1_1_truncated.fastq"
        self.infile2 = "BWA/HNSCC1_2_truncated.fastq"
        self.saifile1 = "BWA/1.sai"
        self.saifile2 = "BWA/2.sai"
        self.samfile1 = "BWA/1.sam"
        self.samfile2 = "BWA/2.sam"
        self.samfile = "BWA/out.sam"
        self.files_to_clean = [
            self.saifile1,
            self.saifile2,
            self.samfile1,
            self.samfile2,
            self.samfile,
        ]

    def tearDown(self):
        for filename in self.files_to_clean:
            if os.path.isfile(filename):
                os.remove(filename)
        for extension in self.reference_extensions:
            index_file = self.reference_file + "." + extension
            if os.path.exists(index_file):
                os.remove(index_file)

    def test_index(self):
        """Test for creating index files for the reference genome fasta file."""
        cmdline = BwaIndexCommandline(bwa_exe)
        cmdline.set_parameter("infile", self.reference_file)
        cmdline.set_parameter("algorithm", "bwtsw")
        stdout, stderr = cmdline()
        for extension in self.reference_extensions:
            index_file = self.reference_file + "." + extension
            self.assertTrue(
                os.path.exists(index_file), f"Index File {index_file} not found"
            )
        self.assertIn(
            "Finished constructing BWT",
            str(stdout) + str(stderr),
            f"FASTA indexing failed:\n{cmdline}\nStdout:{stdout}\nStderr:{stderr}\n",
        )

    def do_aln(self, in_file, out_file):
        """Test for generating sai files given the reference and read file."""
        cmdline = BwaAlignCommandline(bwa_exe)
        cmdline.set_parameter("reference", self.reference_file)
        cmdline.read_file = in_file
        self.assertTrue(os.path.isfile(in_file))
        stdout, stderr = cmdline(stdout=out_file)
        self.assertNotIn(
            "fail to locate the index",
            str(stderr) + str(stdout),
            "Error aligning sequence to reference:\n%s\nStdout:%s\nStderr:%s\n"
            % (cmdline, stdout, stderr),
        )

    def create_fasta_index(self):
        """Test for generating index for fasta file.

        BWA requires an indexed fasta for each alignment operation.
        This should be called to create an index before any alignment
        operation.
        """
        cmdline = BwaIndexCommandline(bwa_exe)
        cmdline.set_parameter("infile", self.reference_file)
        cmdline.set_parameter("algorithm", "bwtsw")
        stdout, stderr = cmdline()

    if not skip_aln_tests:

        def test_samse(self):
            """Test for single end sequencing."""
            self.create_fasta_index()
            self.do_aln(self.infile1, self.saifile1)
            cmdline = BwaSamseCommandline(bwa_exe)
            cmdline.set_parameter("reference", self.reference_file)
            cmdline.set_parameter("read_file", self.infile1)
            cmdline.set_parameter("sai_file", self.saifile1)
            stdout, stderr = cmdline(stdout=self.samfile1)

            with open(self.samfile1) as handle:
                headline = handle.readline()
            self.assertTrue(
                headline.startswith("@SQ"),
                f"Error generating sam files:\n{cmdline}\nOutput starts:{headline}",
            )

        def test_sampe(self):
            """Test for generating samfile by paired end sequencing."""
            self.create_fasta_index()

            # Generate sai files from paired end data
            self.do_aln(self.infile1, self.saifile1)
            self.do_aln(self.infile2, self.saifile2)

            cmdline = BwaSampeCommandline(bwa_exe)
            cmdline.set_parameter("reference", self.reference_file)
            cmdline.set_parameter("sai_file1", self.saifile1)
            cmdline.set_parameter("sai_file2", self.saifile2)
            cmdline.set_parameter("read_file1", self.infile1)
            cmdline.set_parameter("read_file2", self.infile2)
            stdout, stderr = cmdline(stdout=self.samfile)

            with open(self.samfile) as handle:
                headline = handle.readline()
            self.assertTrue(
                headline.startswith("@SQ"),
                f"Error generating sam files:\n{cmdline}\nOutput starts:{headline}",
            )

        def test_mem(self):
            """Test for generating samfile by paired end sequencing using BWA-MEM."""
            self.create_fasta_index()

            cmdline = BwaMemCommandline(bwa_exe)
            cmdline.set_parameter("reference", self.reference_file)
            cmdline.set_parameter("read_file1", self.infile1)
            cmdline.set_parameter("read_file2", self.infile2)
            stdout, stderr = cmdline(stdout=self.samfile)

            with open(self.samfile) as handle:
                headline = handle.readline()
            self.assertTrue(
                headline.startswith("@SQ"),
                f"Error generating sam files:\n{cmdline}\nOutput starts:{headline}",
            )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
