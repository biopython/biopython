# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
Unittests for Bio.Align.Applications interface for MAFFT

This code is part of the Biopython distribution and governed by its
license.  Please see the LICENSE file that should have been included
as part of this package.
"""

import sys
import os
import unittest
from Bio import Application
from Bio import MissingExternalDependencyError
from Bio.Align.Applications import MafftCommandline

mafft_exe = None
if sys.platform=="win32":
    raise MissingExternalDependencyError("Testing with MAFFT not implemented on Windows yet")
else:
    import commands
    output = commands.getoutput("mafft -help")
    if "not found" not in output and "MAFFT" in output:
        mafft_exe = "mafft"
if not mafft_exe:
    raise MissingExternalDependencyError(\
        "Install MAFFT if you want to use the Bio.Align.Applications wrapper.")

class MafftApplication(unittest.TestCase):

    def setUp(self):
        self.infile1  = "Fasta/f002"

    def tearDown(self):
        if os.path.isfile("Fasta/f002.tree"):
            os.remove("Fasta/f002.tree")

    def test_Mafft_simple(self):
        """Simple round-trip through app with infile.
        Result passed to stdout.
        """
        #Use a keyword argument at init,
        cmdline = MafftCommandline(mafft_exe, input=self.infile1)
        self.assertEqual(str(eval(repr(cmdline))), str(cmdline))
        result, stdout, stderr = Application.generic_run(cmdline)
        stderr_string = stderr.read()
        self.assertEqual(result.return_code, 0)
        self.assert_(stdout.read().startswith(">gi|1348912|gb|G26680|G26680"))
        self.assert_("STEP     2 / 2 d" in stderr_string)
        self.assert_("$#=0" not in stderr_string)
        self.assertEqual(str(result._cl), mafft_exe \
                         + " Fasta/f002")

    def test_Mafft_with_options(self):
        """Simple round-trip through app with infile and options.
        Result passed to stdout.
        """
        cmdline = MafftCommandline(mafft_exe)
        cmdline.set_parameter("input", self.infile1)
        cmdline.set_parameter("maxiterate", 100)
        cmdline.set_parameter("--localpair", True)
        self.assertEqual(str(eval(repr(cmdline))), str(cmdline))
        result, stdout, stderr = Application.generic_run(cmdline)
        self.assertEqual(result.return_code, 0)
        self.assert_(stdout.read().startswith(">gi|1348912|gb|G26680|G26680"))
        self.assert_("$#=0" not in stderr.read())
        self.assertEqual(str(result._cl), mafft_exe \
                         + " --localpair --maxiterate 100 Fasta/f002")

    def test_Mafft_with_Clustalw_output(self):
        """Simple round-trip through app with clustal output"""
        cmdline = MafftCommandline(mafft_exe)
        #Use some properties:
        cmdline.input = self.infile1
        cmdline.clustalout = True
        self.assertEqual(str(eval(repr(cmdline))), str(cmdline))
        result, stdout, stderr = Application.generic_run(cmdline)
        self.assertEqual(result.return_code, 0)
        output = stdout.read()
        #e.g. "CLUSTAL format alignment by MAFFT ..."
        #or "CLUSTAL (-like) formatted alignment by MAFFT FFT-NS-2 (v6.240)"
        self.assert_(output.startswith("CLUSTAL"), output)
        self.assert_("$#=0" not in stderr.read())
        self.assertEqual(str(result._cl), mafft_exe \
                         + " --clustalout Fasta/f002")

    def test_Mafft_with_complex_command_line(self):
        """Round-trip with complex command line."""
        cmdline = MafftCommandline(mafft_exe)
        cmdline.set_parameter("input", self.infile1)
        cmdline.set_parameter("--localpair", True)
        cmdline.set_parameter("--weighti", 4.2)
        cmdline.set_parameter("retree", 5)
        cmdline.set_parameter("maxiterate", 200)
        cmdline.set_parameter("--nofft", True)
        cmdline.set_parameter("op", 2.04)
        cmdline.set_parameter("--ep", 0.51)
        cmdline.set_parameter("--lop", 0.233)
        cmdline.set_parameter("lep", 0.2)
        cmdline.set_parameter("--reorder", True)
        cmdline.set_parameter("--treeout", True)
        cmdline.set_parameter("nuc", True)
        self.assertEqual(str(eval(repr(cmdline))), str(cmdline))
        result, stdout, stderr = Application.generic_run(cmdline)
        self.assertEqual(result.return_code, 0)
        self.assert_(stdout.read().startswith(">gi|1348912|gb|G26680|G26680"))
        self.assert_("$#=0" not in stderr.read())
        self.assertEqual(str(result._cl), mafft_exe \
                         + " --localpair --weighti 4.2 --retree 5 " \
                         + "--maxiterate 200 --nofft --op 2.04 --ep 0.51" \
                         + " --lop 0.233 --lep 0.2 --reorder --treeout" \
                         + " --nuc Fasta/f002")

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
