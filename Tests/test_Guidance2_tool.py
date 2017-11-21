# Copyright 2017 by Rob Gilmore and Shaurita Hutchins. All rights reserved.
# Based on ClustalOmega wrapper copyright 2011 by Andreas Wilm.
#
# Wrapper for Guidance2 by Rob Gilmore (2017).  http://guidance.tau.ac.il/ver2/
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

from Bio import MissingExternalDependencyError

import sys
import os
import shutil
import unittest
import subprocess
from Bio import SeqIO
from Bio import AlignIO
from Bio._py3k import getoutput
from Bio.Align.Applications import Guidance2Commandline

#################################################################

# Try to avoid problems when the OS is in another language
os.environ['LANG'] = 'C'

guidance_exe = None
if sys.platform == "win32":
    raise MissingExternalDependencyError("Testing with Guidance2 not implemented on Windows yet")

try:
    output = getoutput("guidance")
    if output.startswith("USAGE: --seqFile <seqFile>"):
        guidance_exe = "guidance"
except FileNotFoundError:
    pass

if not guidance_exe:
    raise MissingExternalDependencyError(
        "Install guidance and put it in your PATH if you want to use Guidance2 from Biopython.")

##################################################################################################


class Guidance2TestCase(unittest.TestCase):

    def setUp(self):
        self.na_input = {
            "seqFile": "Guidance2/Seqs/HTR1A.ffn",
            "seqType": "nuc",
            "outDir": "Guidance2/Output"
                      }
        self.aa_input = {
            "seqFile": "Guidance2/Seqs/HTR1A.faa",
            "seqType": "aa",
            "outDir": "Guidance2/Output"
        }
        self.output = {
            "outDir": "Guidance2/Output",
            "outFASTA": "Guidance2/Output/Seqs.Orig.fas",
            "outALIGN": "Guidance2/Output/HTR1A.CLUSTALW.aln.orig"
        }

    def tearDown(self):
        outDir = self.output["outDir"]
        shutil.rmtree(outDir)

    def standard_test_procedure(self, cline=Guidance2Commandline):
        """Standard testing procedure used by all tests."""

        output_path = cline.outDir

        self.assertEqual(str(eval(repr(cline))), str(cline))
        guidance2 = subprocess.Popen([str(cline)], stderr=subprocess.PIPE,
                                     stdout=subprocess.PIPE, shell=True,
                                     encoding='utf-8')
        _error = guidance2.stderr.readlines()
        _out = guidance2.stdout.readlines()
        self.assertTrue(not _out or _out.strip().startswith("................."))
        input_records = SeqIO.to_dict(SeqIO.parse(cline.seqFile, "fasta"))

        # Check the output...
        align = AlignIO.read(self.output["outALIGN"], "clustal")
        output_records = SeqIO.to_dict(SeqIO.parse(self.output["outFASTA"], "clustal"))
        self.assertEqual(len(set(input_records.keys())), len(set(output_records.keys())))
        for record in align:
            self.assertEqual(str(record.seq), str(output_records[record.id].seq))


class Guidance2TestNormalConditions(Guidance2TestCase):

    def test_nucleic_acid_fasta(self):
        cline = Guidance2Commandline(msaProgram="CLUSTALW", **self.na_input)
        self.standard_test_procedure(cline)

    def test_amino_acid_fasta(self):
        cline = Guidance2Commandline(msaProgram="CLUSTALW", **self.aa_input)
        self.standard_test_procedure(cline)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
