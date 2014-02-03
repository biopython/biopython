# Copyright 2013 by Saket Choudhary.  Based on test_Clustalw_tool.py by Peter
# Cock .
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
from Bio.Sequencing.Applications import BwaIndexCommandline, BwaAlignCommandline
from Bio.Sequencing.Applications import BwaSamseCommandline, BwaSampeCommandline
from Bio.Sequencing.Applications import BwaBwaswCommandline
from Bio.Application import ApplicationError

#################################################################

#Try to avoid problems when the OS is in another language
os.environ['LANG'] = 'C'

bwa_exe = None
if sys.platform == "win32":
    #TODO - Check the path?
    try:
        #This can vary depending on the Windows language.
        prog_files = os.environ["PROGRAMFILES"]
    except KeyError:
        prog_files = r"C:\Program Files"
    #A default path of C:\Program Files\bwa.exe was chosen
    #but this path can be edited depending on where bwa is located

    likely_dirs = ["bwa", "BWA", "Bwa", "bwa-0.6.2"]
    likely_exes = ["bwa.exe", "BWA.exe", "Bwa.exe"]
    for folder in likely_dirs:
        if os.path.isdir(os.path.join(prog_files, folder)):
            for filename in likely_exes:
                if os.path.isfile(os.path.join(prog_files, folder, filename)):
                    bwa_exe = os.path.join(prog_files, folder, filename)
                    break
            if bwa_exe:
                break
else:
    from Bio._py3k import getoutput
    output = getoutput("bwa")

    #Since "not found" may be in another language, try and be sure this is
    #really the bwa tool's output
    bwa_found = False
    if "not found" not in output and "bwa" in output \
    and "alignment via Burrows-Wheeler transformation" in output:
        bwa_exe = "bwa"

if not bwa_exe:
    raise MissingExternalDependencyError(\
        "Install bwa and correctly set the file path to the program if you want to use it from Biopython")


class BwaTestCase(unittest.TestCase):
    """Class for implementing BWA test cases"""
    def setUp(self):
        self.reference_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "BWA", "human_g1k_v37_truncated.fasta")
        self.infile1 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "BWA", "HNSCC1_1_truncated.fastq")
        self.infile2 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "BWA", "HNSCC1_2_truncated.fastq")
        self.saifile1 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "BWA", "1.sai")
        self.saifile2 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "BWA", "2.sai")
        self.samfile1 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "BWA", "1.sam")
        self.samfile2 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "BWA", "2.sam")
        self.samfile = os.path.join(os.path.dirname(os.path.abspath(__file__)), "BWA", "out.sam")

    def test_index(self):
        """Test for creating index files for the reference genome fasta file"""
        cmdline = BwaIndexCommandline()
        cmdline.set_parameter("infile", self.reference_file)
        cmdline.set_parameter("algorithm", "bwtsw")
        stdout, stderr = cmdline()
        output = stdout.startswith("[bwt_gen]")
        self.assertTrue(stdout.startswith("[bwt_gen]"),
                        "FASTA indexing failed:\n%s\nStdout:%s" \
                        % (cmdline, stdout))

    def do_aln(self, in_file, out_file):
        """Test for generating sai files given the reference and read file"""
        cmdline = BwaAlignCommandline()
        cmdline.set_parameter("reference", self.reference_file)
        cmdline.read_file = in_file
        self.assertTrue(os.path.isfile(in_file))
        stdout, stderr = cmdline(stdout=out_file)

        self.assertTrue("fail to locate the index" not in stderr,
                        "Error aligning sequence to reference:\n%s\nStderr:%s" \
                        % (cmdline, stderr))

    def test_samse(self):
        """Test for single end sequencing """
        cmdline = BwaSamseCommandline()
        cmdline.set_parameter("reference", self.reference_file)
        cmdline.set_parameter("read_file", self.infile1)
        cmdline.set_parameter("sai_file", self.saifile1)
        stdout, stderr = cmdline(stdout=self.samfile1)

        with open(self.samfile1, "r") as handle:
            headline = handle.readline()
        self.assertTrue(headline.startswith("@SQ"),
                        "Error generating sam files:\n%s\nOutput starts:%s" \
                        % (cmdline, headline))

    def test_sampe(self):
        """Test for generating samfile by paired end sequencing"""
        ##Generate sai files from paired end data
        self.do_aln(self.infile1, self.saifile1)
        self.do_aln(self.infile2, self.saifile2)

        cmdline = BwaSampeCommandline()
        cmdline.set_parameter("reference", self.reference_file)
        cmdline.set_parameter("sai_file1", self.saifile1)
        cmdline.set_parameter("sai_file2", self.saifile2)
        cmdline.set_parameter("read_file1", self.infile1)
        cmdline.set_parameter("read_file2", self.infile2)
        stdout, stderr = cmdline(stdout=self.samfile)

        with open(self.samfile, "r") as handle:
            headline = handle.readline()
        self.assertTrue(headline.startswith("@SQ"),
                        "Error generating sam files:\n%s\nOutput starts:%s" \
                        % (cmdline, headline))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
