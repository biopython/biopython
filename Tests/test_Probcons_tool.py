# Copyright 2009 by Cymon J. Cox.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Unittests for Bio.Align.Applications interface for PROBCONS."""

import sys
import os
import unittest
import subprocess
from cStringIO import StringIO
from Bio import AlignIO, SeqIO, MissingExternalDependencyError
from Bio.Align.Applications import ProbconsCommandline

probcons_exe = None
if sys.platform=="win32":
    raise MissingExternalDependencyError("PROBCONS not available on Windows")
else:
    import commands
    output = commands.getoutput("probcons")
    if "not found" not in output and "probcons" in output.lower():
        probcons_exe = "probcons"

if not probcons_exe:
    raise MissingExternalDependencyError(\
        "Install PROBCONS if you want to use the Bio.Align.Applications wrapper.")

class ProbconsApplication(unittest.TestCase):

    def setUp(self):
        self.infile1 = "Fasta/fa01"
        self.annotation_outfile = "Fasta/probcons_annot.out"

    def tearDown(self):
        if os.path.isfile(self.annotation_outfile):
            os.remove(self.annotation_outfile)

    def test_Probcons_alignment_fasta(self):
        """Round-trip through app and read fasta alignment from stdout
        """
        cmdline = ProbconsCommandline(probcons_exe, input=self.infile1)
        self.assertEqual(str(cmdline), probcons_exe + " Fasta/fa01")
        self.assertEqual(str(eval(repr(cmdline))), str(cmdline))
        child = subprocess.Popen(str(cmdline),
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True,
                                 shell=(sys.platform!="win32"))
        return_code = child.wait()
        self.assertEqual(return_code, 0)
        self.assertTrue(child.stderr.read().startswith("\nPROBCONS"))
        align = AlignIO.read(StringIO(child.stdout.read()), "fasta")
        records = list(SeqIO.parse(open(self.infile1),"fasta"))
        self.assertEqual(len(records),len(align))
        for old, new in zip(records, align):
            self.assertEqual(old.id, new.id)
            self.assertEqual(str(new.seq).replace("-",""), str(old.seq).replace("-",""))
        del child

    def test_Probcons_alignment_clustalw(self):
        """Round-trip through app and read clustalw alignment from stdout
        """
        cmdline = ProbconsCommandline(probcons_exe)
        cmdline.set_parameter("input", "Fasta/fa01")
        cmdline.clustalw = True
        self.assertEqual(str(cmdline), probcons_exe + " -clustalw Fasta/fa01")
        self.assertEqual(str(eval(repr(cmdline))), str(cmdline))
        child = subprocess.Popen(str(cmdline),
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True,
                                 shell=(sys.platform!="win32"))
        return_code = child.wait()
        self.assertEqual(return_code, 0)
        self.assertTrue(child.stderr.read().strip().startswith("PROBCONS"))
        #self.assertTrue(stdout.read().strip().startswith("PROBCONS"))
        align = AlignIO.read(StringIO(child.stdout.read()), "clustal")
        records = list(SeqIO.parse(open(self.infile1),"fasta"))
        self.assertEqual(len(records),len(align))
        for old, new in zip(records, align):
            self.assertEqual(old.id, new.id)
            self.assertEqual(str(new.seq).replace("-",""), str(old.seq).replace("-",""))
        del child

    def test_Probcons_complex_commandline(self):
        """Round-trip through app with complex command line and output file
        """
        cmdline = ProbconsCommandline(probcons_exe, pre=1)
        cmdline.set_parameter("input", "Fasta/fa01")
        cmdline.consistency = 4
        cmdline.set_parameter("--iterative-refinement", 222)
        cmdline.set_parameter("a", True)
        cmdline.annot = self.annotation_outfile
        self.assertEqual(str(cmdline), probcons_exe +
                " -c 4 -ir 222 -pre 1 -annot Fasta/probcons_annot.out "
                "-a Fasta/fa01")
        child = subprocess.Popen(str(cmdline),
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True,
                                 shell=(sys.platform!="win32"))
        return_code = child.wait()
        self.assertEqual(return_code, 0)
        self.assertTrue(child.stderr.read().startswith("\nPROBCONS"))
        self.assertTrue(child.stdout.read().startswith(">AK1H_ECOLI/1-378"))
        del child

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
