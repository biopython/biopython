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

prank_exe = None
if sys.platform=="win32":
    try :
        #This can vary depending on the Windows language.
        prog_files = os.environ["PROGRAMFILES"]
    except KeyError :
        prog_files = r"C:\Program Files"
    #For Windows, PRANK just comes as a zip file which contains the
    #prank.exe file which the user could put anywhere.  We'll try a few
    #sensible locations under Program Files... and then the full path.
    likely_dirs = ["", #Current dir
                   prog_files,
                   os.path.join(prog_files,"Prank")] + sys.path
    for folder in likely_dirs :
        if os.path.isdir(folder) :
            if os.path.isfile(os.path.join(folder, "prank.exe")) :
                prank_exe = os.path.join(folder, "prank.exe")
                break
        if prank_exe : break
else :
    import commands
    output = commands.getoutput("prank")
    if "not found" not in output and "prank" in output.lower():
        prank_exe = "prank"
if not prank_exe:
    raise MissingExternalDependencyError(\
        "Install PRANK if you want to use the Bio.Align.Applications wrapper.")

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
        cmdline = PrankCommandline(prank_exe)
        cmdline.set_parameter("d", self.infile1)
        self.assertEqual(str(cmdline), prank_exe + " -d=Fasta/f002 ")
        stdin, stdout, stderr = Application.generic_run(cmdline)
        self.assertEqual(stdin.return_code, 0)
        self.assert_("Total time" in stdout.read())
        self.assertEqual(stderr.read(), "")
        self.assertEqual(str(stdin._cl), str(cmdline))

    def test_Prank_simple_with_NEXUS_output(self):
        """Simple round-trip through app with infile, output in NEXUS
        output.?.??? files written to cwd - no way to redirect
        """
        cmdline = PrankCommandline(prank_exe)
        cmdline.set_parameter("d", self.infile1)
        cmdline.set_parameter("f", 17) #17. NEXUS FORMAT
        cmdline.set_parameter("-noxml")
        cmdline.set_parameter("notree")
        self.assertEqual(str(cmdline), prank_exe + " -d=Fasta/f002 -f=17 -noxml -notree ")
        stdin, stdout, stderr = Application.generic_run(cmdline)
        self.assertEqual(stdin.return_code, 0)
        self.assert_("Total time" in stdout.read())
        self.assertEqual(stderr.read(), "")
        self.assertEqual(str(stdin._cl), str(cmdline))

    def test_Prank_complex_command_line(self):
        """Round-trip with complex command line."""
        cmdline = PrankCommandline(prank_exe)
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
        self.assertEqual(str(cmdline), prank_exe + " -d=Fasta/f002 -noxml" + \
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
