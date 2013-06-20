# Copyright 2013 by Saket Choudhary.  Based on test_Clustalw_tool.py by Peter
# Cock .
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from __future__ import with_statement

from Bio import MissingExternalDependencyError
import sys
import os
import unittest
from Bio import SeqIO
from Bio import AlignIO
from Bio.Sequencing.Applications import SamtoolsViewCommandline,SamtoolsCalmdCommandline
from Bio.Sequencing.Applications import SamtoolsCatCommandline, SamtoolsFaidxCommandline
from Bio.Sequencing.Applications import SamtoolsFixmateCommandline, SamtoolsIdxstatsCommandline
from Bio.Sequencing.Applications import SamtoolsIndexCommandline, SamtoolsMergeCommandline
from Bio.Sequencing.Applications import SamtoolsMpileupCommandline, SamtoolsPhaseCommandline
from Bio.Sequencing.Applications import SamtoolsReheaderCommandline, SamtoolsRmdupCommandline
from Bio.Sequencing.Applications import SamtoolsSortCommandline, SamtoolsTargetcutCommandline
from Bio.Application import ApplicationError

#################################################################

#Try to avoid problems when the OS is in another language
os.environ['LANG'] = 'C'

samtools_exe = None
if sys.platform == "win32":
    #TODO - Check the path?
    try:
        #This can vary depending on the Windows language.
        prog_files = os.environ["PROGRAMFILES"]
    except KeyError:
        prog_files = r"C:\Program Files"
    #A default path of C:\Program Files\samtools.exe was chosen
    #but this path can be edited depending on where samtools is located

    likely_dirs = ["samtools", "Samtools", "SAMtools"]
    likely_exes = ["samtools.exe", "Samtools.exe", "SAMtools.exe"]
    for folder in likely_dirs:
        if os.path.isdir(os.path.join(prog_files, folder)):
            for filename in likely_exes:
                if os.path.isfile(os.path.join(prog_files, folder, filename)):
                    samtools_exe = os.path.join(prog_files, folder, filename)
                    break
            if samtools_exe:
                break
else:
    import commands
    output = commands.getoutput("samtools")

    #Since "not found" may be in another language, try and be sure this is
    #really the samtools tool's output
    samtools_found = False
    if "not found" not in output and "samtools" in output \
    and "samtools (Tools for alignments in the SAM format)" in output:
        samtools_exe = "bwa"

if not samtools_exe:
    raise MissingExternalDependencyError(\
        "Install samtools and correctly set the file path to the program if you want to use it from Biopython")


class SamtoolsTestCase(unittest.TestCase):
    """Class for implementing Samtools test cases."""

    def setUp(self):
        self.files_to_clean = set()
        self.samfile1 = os.path.join("SamBam", "sam1.sam")
        self.reference = os.path.join("BWA", "human_g1k_v37_truncated.fasta")
        self.referenceindexfile = os.path.join("BWA", "human_g1k_v37_truncated.fasta.fai")
        self.samfile2 = os.path.join("SamBam", "sam2.sam")
        self.bamfile1 = os.path.join("SamBam", "bam1.bam")
        self.bamfile2 = os.path.join("SamBam", "bam2.bam")
        self.outsamfile = os.path.join("SamBam", "out.sam")
        self.outbamfile = os.path.join("SamBam", "out.bam")
        self.bamindexfile1 = os.path.join("SamBam", "bam1.bam.bai")

    def tearDown(self):
        for filename in self.files_to_clean:
            if os.path.isfile(filename):
                os.remove(filename)

    def add_files_to_clean(self,filename):
        self.files_to_clean.add(filename)



    def test_view(self):
        """Test for samtools view"""

        cmdline = SamtoolsViewCommandline()
        cmdline.set_parameter("input_file", self.bamfile1)
        stdout_bam,stderr_bam = cmdline()
        self.assertTrue(stderr_bam.startswith(""),"SAM file viewing failed: \n%s\nStdout:%s" \
                        % (cmdline, stdout_bam))
        cmdline.set_parameter("input_file", self.samfile1)
        cmdline.set_parameter("S", True)
        stdout_sam, stderr_sam = cmdline()
        self.assertTrue(stderr_sam.startswith("[samopen] SAM header is present:"),
                        "SAM file  viewing failed:\n%s\nStderr:%s" \
                        % (cmdline, stderr_sam))


    def test_faidx(self):
        cmdline = SamtoolsFaidxCommandline()
        cmdline.set_parameter("reference", self.reference)
        stdout, stderr = cmdline()
        self.assertFalse(stderr, "Samtools faidx failed:\n%s\nStderr:%s"\
                        % (cmdline, stderr))

        self.assertTrue(os.path.isfile(self.referenceindexfile))

    def test_calmd(self):
        """Test for samtools calmd"""
        cmdline = SamtoolsCalmdCommandline()
        cmdline.set_parameter("reference", self.reference)
        cmdline.set_parameter("input_bam", self.bamfile1)
        ## If there is no index file for the reference
        ## samtools calmd creates one at the time of calling

        if os.path.exists(self.referenceindexfile):
            print "exists"
            stderr_calmd_expected = ""
        else:
            print "doesnt exist"
            stderr_calmd_expected = "[fai_load] build FASTA index.\n"

        stdout, stderr = cmdline()
        if stderr == stderr_calmd_expected:
            print "SAME"
        else:
            print "NOT SAME"
            print stderr
            print stderr_calmd_expected
        self.assertTrue(stderr==stderr_calmd_expected, "Samtools calmd failed:\n%s\nStderr:%s" \
                        % (cmdline, stderr))

    def test_cat(self):
        cmdline = SamtoolsCatCommandline()
        cmdline.set_parameter("o", self.outbamfile)
        cmdline.set_parameter("input_bam", [self.bamfile1,self.bamfile2])
        stdout, stderr = cmdline()
        self.assertTrue(stderr=="", "Samtools cat failed:\n%s\nStderr:%s"\
                        % (cmdline, stderr))
        self.assertTrue(os.path.exists(self.outbamfile))
        self.add_files_to_clean(self.outbamfile)



    def test_fixmate(self):
        ##TODO : Needs a name-sorted alignment file
        pass


    def test_sort(self):
        cmdline = SamtoolsSortCommandline()
        cmdline.set_parameter("input_bam", self.bamfile1)
        cmdline.set_parameter("out_prefix", "bam1")
        stdout, stderr = cmdline()
        self.assertFalse(stderr, "Samtools sort failed:\n%s\nStderr:%s"\
                        % (cmdline, stderr))


    def test_index(self):
        cmdline = SamtoolsIndexCommandline()
        cmdline.set_parameter("input_bam", self.bamfile1)
        stdout, stderr = cmdline()
        self.assertFalse(stderr, "Samtools index failed:\n%s\nStderr:%s"\
                        % (cmdline, stderr))

        self.assertTrue(os.path.exists(self.bamindexfile1))

    def test_idxstats(self):
        self.test_index()
        cmdline = SamtoolsIdxstatsCommandline()
        cmdline.set_parameter("input_bam", self.bamfile1)
        stdout, stderr = cmdline()
        self.assertFalse(stderr, "Samtools idxstats failed:\n%s\nStderr:%s"\
                        % (cmdline, stderr))



    def test_merge(self):
        cmdline = SamtoolsMergeCommandline()
        cmdline.set_parameter("input_bam",[self.bamfile1, self.bamfile2])
        cmdline.set_parameter("out_bam", self.outbamfile)
        cmdline.set_parameter("f", True) ## Overwrite out.bam if it exists
        stdout, stderr = cmdline()

        self.assertFalse(stderr, "Samtools merge failed:\n%s\nStderr:%s"\
                        % (cmdline, stderr)    )

        self.assertTrue(os.path.exists(self.outbamfile))
        self.add_files_to_clean(self.outbamfile)


    def test_mpileup(self):
        cmdline = SamtoolsMpileupCommandline()
        cmdline.set_parameter("input_file", self.bamfile1)
        stdout, stderr = cmdline()
        self.assertFalse("[bam_pileup_core]" in stdout)



    def test_phase(self):
        pass

    def test_reheader(self):
        pass

    def test_rmdup(self):
        pass

    def test_targetcut(self):
        pass




if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
