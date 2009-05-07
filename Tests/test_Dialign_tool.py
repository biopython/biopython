# Copyright 2009 by Cymon J. Cox.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Unittests for Bio.Align.Applications interface for DIALIGN2-2."""

import sys
import os
import unittest
from Bio import Application
from Bio import MissingExternalDependencyError
from Bio.Align.Applications import DialignCommandline

dialign_exe = None
if sys.platform=="win32":
    raise MissingExternalDependencyError("DIALIGN2-2 not available on Windows")
else :
    import commands
    output = commands.getoutput("dialign2-2")
    if "not found" not in output and "dialign2-2" in output.lower():
        dialign_exe = "dialign2-2"
        #This check is currently not needed, as if the environment variable
        #is missing the tool outputs this:
        #
        #    Please set the environmentvariable DIALIGN2_DIR 
        #    as described in the README file 
        #
        #if "DIALIGN2_DIR" not in os.environ :
        #    raise MissingExternalDependencyError(\
        #        "Environment variable DIALIGN2_DIR for DIALIGN2-2 missing.")

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
        stdin, stdout, stderr = Application.generic_run(cmdline)
        
        self.assert_(stdin.return_code == 0)
        self.assert_(os.path.exists(self.outfile1))
        self.assert_(stdout.read() == "")
        self.assert_(stderr.read() == "")
        self.assert_(str(stdin._cl) == "dialign2-2 Fasta/f002 ")

    def test_Dialign_simple_with_options(self):
        """Simple round-trip through app with infile and options
        """
        cmdline = DialignCommandline(dialign_exe)
        cmdline.set_parameter("input", self.infile1)
        cmdline.set_parameter("-max_link", True)
        cmdline.set_parameter("stars", 4)
        stdin, stdout, stderr = Application.generic_run(cmdline)
        
        self.assert_(stdin.return_code == 0)
        self.assert_(os.path.exists(self.outfile1))
        self.assert_(stdout.read() == "")
        self.assert_(stderr.read() == "")
        self.assert_(str(stdin._cl) == "dialign2-2 -max_link -stars 4 Fasta/f002 ")

    def test_Dialign_simple_with_MSF_output(self):
        """Simple round-trip through app with infile, output MSF
        """
        cmdline = DialignCommandline(dialign_exe)
        #Test with properties
        cmdline.input = self.infile1
        cmdline.msf = True
        stdin, stdout, stderr = Application.generic_run(cmdline)
        
        self.assert_(stdin.return_code == 0)
        self.assert_(os.path.exists(self.outfile1))
        self.assert_(os.path.exists(self.outfile2))
        self.assert_(stdout.read() == "")
        self.assert_(stderr.read() == "")
        self.assert_(str(stdin._cl) == "dialign2-2 -msf Fasta/f002 ")

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

        stdin, stdout, stderr = Application.generic_run(cmdline)
        
        self.assert_(stdin.return_code == 0)
        self.assert_(os.path.exists(self.outfile1))
        self.assert_(stdout.read().startswith(" e_len = 633"))
        self.assert_(stderr.read() == "")
        self.assert_(str(stdin._cl) == "dialign2-2 -cs -mask -nt -ow -stars 9 -thr 4 Fasta/f002 ")

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
