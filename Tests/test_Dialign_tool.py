# Copyright 2009 by Cymon J. Cox.  All rights reserved.
# Revisions 2009 copyright by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Unittests for Bio.Align.Applications interface for DIALIGN2-2."""

import sys
import os
import unittest
import subprocess
from Bio import MissingExternalDependencyError
from Bio.Align.Applications import DialignCommandline

dialign_exe = None
if sys.platform=="win32":
    raise MissingExternalDependencyError("DIALIGN2-2 not available on Windows")
else:
    import commands
    output = commands.getoutput("dialign2-2")
    if "not found" not in output and "dialign2-2" in output.lower():
        dialign_exe = "dialign2-2"
        if "DIALIGN2_DIR" not in os.environ:
            raise MissingExternalDependencyError(\
                "Environment variable DIALIGN2_DIR for DIALIGN2-2 missing.")
        if not os.path.isdir(os.environ["DIALIGN2_DIR"]):
            raise MissingExternalDependencyError(\
                "Environment variable DIALIGN2_DIR for DIALIGN2-2 is not a valid directory.")
        if not os.path.isfile(os.path.join(os.environ["DIALIGN2_DIR"], "BLOSUM")):
            raise MissingExternalDependencyError(\
                "Environment variable DIALIGN2_DIR directory missing BLOSUM file.")
        #TODO - check for tp400_dna, tp400_prot and tp400_trans too?
        
if not dialign_exe:
    raise MissingExternalDependencyError(\
        "Install DIALIGN2-2 if you want to use the Bio.Align.Applications wrapper.")

class DialignApplication(unittest.TestCase):

    def setUp(self):
        self.infile1 = "Fasta/f002" 
        #Standard output file
        self.outfile1 = "Fasta/f002.ali"
        #MSF output
        self.outfile2 = "Fasta/f002.ms"

    def tearDown(self):
        if os.path.isfile(self.outfile1):
            os.remove(self.outfile1)
        if os.path.isfile(self.outfile2):
            os.remove(self.outfile2)

    def test_Dialign_simple(self):
        """Simple round-trip through app with infile.
        """
        #Test using keyword arguments:
        cmdline = DialignCommandline(dialign_exe, input=self.infile1)
        self.assertEqual(str(cmdline), dialign_exe + " Fasta/f002")
        child = subprocess.Popen(str(cmdline),
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True,
                                 shell=(sys.platform!="win32"))
        return_code = child.wait()
        self.assertEqual(return_code, 0)
        self.assertEqual(child.stderr.read(), "")
        self.assertEqual(child.stdout.read(), "")
        self.assertTrue(os.path.exists(self.outfile1))
        del child

    def test_Dialign_simple_with_options(self):
        """Simple round-trip through app with infile and options
        """
        cmdline = DialignCommandline(dialign_exe)
        cmdline.set_parameter("input", self.infile1)
        cmdline.set_parameter("-max_link", True)
        cmdline.set_parameter("stars", 4)
        self.assertEqual(str(cmdline), dialign_exe + \
                         " -max_link -stars 4 Fasta/f002")
        child = subprocess.Popen(str(cmdline),
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True,
                                 shell=(sys.platform!="win32"))
        return_code = child.wait()
        self.assertEqual(return_code, 0)
        self.assertEqual(child.stderr.read(), "")
        self.assertEqual(child.stdout.read(), "")
        self.assertTrue(os.path.exists(self.outfile1))
        del child

    def test_Dialign_simple_with_MSF_output(self):
        """Simple round-trip through app with infile, output MSF
        """
        cmdline = DialignCommandline(dialign_exe)
        #Test with properties
        cmdline.input = self.infile1
        cmdline.msf = True
        self.assertEqual(str(cmdline), dialign_exe + " -msf Fasta/f002")
        child = subprocess.Popen(str(cmdline),
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True,
                                 shell=(sys.platform!="win32"))
        return_code = child.wait()
        self.assertEqual(return_code, 0)
        self.assertEqual(child.stdout.read(), "")
        self.assertEqual(child.stderr.read(), "")
        self.assertTrue(os.path.exists(self.outfile1))
        self.assertTrue(os.path.exists(self.outfile2))
        del child

    def test_Dialign_complex_command_line(self):
        """Round-trip through app with complex command line."""
        cmdline = DialignCommandline(dialign_exe)
        cmdline.set_parameter("input", self.infile1)
        cmdline.set_parameter("-nt", True)
        cmdline.set_parameter("-thr", 4)
        cmdline.set_parameter("stars", 9)
        cmdline.set_parameter("-ow", True)
        cmdline.set_parameter("mask", True)
        cmdline.set_parameter("-cs", True)
        self.assertEqual(str(cmdline), dialign_exe + \
                         " -cs -mask -nt -ow -stars 9 -thr 4 Fasta/f002")
        child = subprocess.Popen(str(cmdline),
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True,
                                 shell=(sys.platform!="win32"))
        return_code = child.wait()
        self.assertEqual(return_code, 0)
        self.assertEqual(child.stderr.read(), "")
        self.assertTrue(os.path.exists(self.outfile1))
        self.assertTrue(child.stdout.read().startswith(" e_len = 633"))
        del child

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
