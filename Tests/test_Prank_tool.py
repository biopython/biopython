"""
Unittests for Bio.Align.Applications interface for PRANK

This code is part of the Biopython distribution and governed by its
license.  Please see the LICENSE file that should have been included
as part of this package.
"""

import sys
import os
import unittest
from Bio import Application
from Bio import MissingExternalDependencyError
from Bio.Align.Applications import PrankCommandline

app_name = "prank"
if sys.platform=="win32" :
    raise MissingExternalDependencyError("Testing with PRANK not implemented on Windows yet")
else :
    import commands
    if "not found" in commands.getoutput("%s -help" % app_name):
        raise MissingExternalDependencyError(\
            "Alignment application PRANK not found.")

class PrankApplication(unittest.TestCase):
    
    def setUp(self):
        self.infile1 = "Fasta/f002"

    def tearDown(self):
        """
        output.1.dnd  output.1.fas  output.1.xml  output.2.dnd  output.2.fas  output.2.xml
        """
        if os.path.isfile("output.1.dnd"):
            os.remove("output.1.dnd")
        if os.path.isfile("output.1.fas"):
            os.remove("output.1.fas")
        if os.path.isfile("output.1.xml"):
            os.remove("output.1.xml")
        if os.path.isfile("output.2.dnd"):
            os.remove("output.2.dnd")
        if os.path.isfile("output.2.fas"):
            os.remove("output.2.fas")
        if os.path.isfile("output.2.xml"):
            os.remove("output.2.xml")
        if os.path.isfile("output.1.nex"):
            os.remove("output.1.nex")
        if os.path.isfile("output.2.nex"):
            os.remove("output.2.nex")

    def test_Prank_simple(self):
        """Simple round-trip through app with infile.
        output.?.??? files written to cwd - no way to redirect
        """
        cmdline = PrankCommandline()
        cmdline.set_parameter("d", self.infile1)
        self.assertEqual(str(cmdline), app_name + " -d=Fasta/f002 ")
        stdin, stdout, stderr = Application.generic_run(cmdline)
        self.assertEqual(stdin.return_code, 0)
        self.assert_("Total time" in stdout.read())
        self.assertEqual(stderr.read(), "")
        self.assertEqual(str(stdin._cl), str(cmdline))

    def test_Prank_simple_with_NEXUS_output(self):
        """Simple round-trip through app with infile, output in NEXUS
        output.?.??? files written to cwd - no way to redirect
        """
        cmdline = PrankCommandline()
        cmdline.set_parameter("d", self.infile1)
        cmdline.set_parameter("f", 17) #17. NEXUS FORMAT
        cmdline.set_parameter("-noxml")
        cmdline.set_parameter("notree")
        self.assertEqual(str(cmdline), app_name + " -d=Fasta/f002 -f=17 -noxml -notree ")
        stdin, stdout, stderr = Application.generic_run(cmdline)
        self.assertEqual(stdin.return_code, 0)
        self.assert_("Total time" in stdout.read())
        self.assertEqual(stderr.read(), "")
        self.assertEqual(str(stdin._cl), str(cmdline))

    def test_Prank_complex_command_line(self):
        """Round-trip with complex command line."""
        cmdline = PrankCommandline()
        cmdline.set_parameter("d", self.infile1)
        cmdline.set_parameter("-noxml")
        cmdline.set_parameter("notree")
        cmdline.set_parameter("-gaprate", 0.321)
        cmdline.set_parameter("gapext", 0.6)
        cmdline.set_parameter("-dots")
        cmdline.set_parameter("-kappa", 3)
        cmdline.set_parameter("-skipins")
        cmdline.set_parameter("-once")
        cmdline.set_parameter("realbranches")
        self.assertEqual(str(cmdline), app_name + " -d=Fasta/f002 -noxml" + \
                         " -notree -dots -gaprate=0.321 -gapext=0.6 -kappa=3" + \
                         " -once -skipins -realbranches ")
        stdin, stdout, stderr = Application.generic_run(cmdline)
        self.assertEqual(stdin.return_code, 0)
        self.assert_("Total time" in stdout.read())
        self.assertEqual(stderr.read(), "")
        self.assertEqual(str(stdin._cl), str(cmdline))

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
