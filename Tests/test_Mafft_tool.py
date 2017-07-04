# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Unittests for Bio.Align.Applications interface for MAFFT."""

import sys
import os
import unittest
import subprocess
from Bio import MissingExternalDependencyError
from Bio.Align.Applications import MafftCommandline

# Try to avoid problems when the OS is in another language
os.environ['LANG'] = 'C'

mafft_exe = None
if sys.platform == "win32":
    raise MissingExternalDependencyError("Testing with MAFFT not implemented on Windows yet")
else:
    from Bio._py3k import getoutput
    output = getoutput("mafft -help")
    if "not found" not in output and "MAFFT" in output:
        mafft_exe = "mafft"
if not mafft_exe:
    raise MissingExternalDependencyError(
        "Install MAFFT if you want to use the Bio.Align.Applications wrapper.")


def check_mafft_version(mafft_exe):
    child = subprocess.Popen("%s --help" % mafft_exe,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             universal_newlines=True,
                             shell=(sys.platform != "win32"))
    stdoutdata, stderrdata = child.communicate()
    output = stdoutdata + "\n" + stderrdata
    return_code = child.returncode
    del child
    if "correctly installed?" in output \
    or "mafft binaries have to be installed" in output:
        raise MissingExternalDependencyError(
            "MAFFT does not seem to be correctly installed.")

    # e.g. "MAFFT version 5.732 (2005/09/14)\n"
    # e.g. "  MAFFT v6.717b (2009/12/03)\n"
    for marker in ["MAFFT version", "MAFFT v"]:
        index = output.find(marker)
        if index == -1:
            continue
        version = output[index + len(marker):].strip().split(None, 1)[0]
        major = int(version.split(".", 1)[0])
        if major < 6:
            raise MissingExternalDependencyError("Test requires MAFFT v6 or "
                                                 "later (found %s)." % version)
        return (major, version)
    raise MissingExternalDependencyError("Couldn't determine MAFFT version.")


# This also checks it actually runs!
version_major, version_string = check_mafft_version(mafft_exe)


class MafftApplication(unittest.TestCase):

    def setUp(self):
        self.infile1 = "Fasta/f002"

    def tearDown(self):
        if os.path.isfile("Fasta/f002.tree"):
            os.remove("Fasta/f002.tree")

    def test_Mafft_simple(self):
        """Simple round-trip through app with infile, result passed to stdout."""
        # Use a keyword argument at init,
        cmdline = MafftCommandline(mafft_exe, input=self.infile1)
        self.assertEqual(str(eval(repr(cmdline))), str(cmdline))
        stdoutdata, stderrdata = cmdline()
        self.assertTrue(stdoutdata.startswith(">gi|1348912|gb|G26680|G26680"))
        # Used to get "Progressive alignment ..." but in v7.245
        # became "Progressive alignment 1/2..." and "Progressive alignment 2/2..."
        self.assertTrue(("Progressive alignment ..." in stderrdata) or
                        ("Progressive alignment 1/" in stderrdata), stderrdata)
        self.assertNotIn("$#=0", stderrdata)

    def test_Mafft_with_options(self):
        """Simple round-trip through app with infile and options, result passed to stdout."""
        cmdline = MafftCommandline(mafft_exe)
        cmdline.set_parameter("input", self.infile1)
        cmdline.set_parameter("maxiterate", 100)
        cmdline.set_parameter("--localpair", True)
        self.assertEqual(str(eval(repr(cmdline))), str(cmdline))
        stdoutdata, stderrdata = cmdline()
        self.assertTrue(stdoutdata.startswith(">gi|1348912|gb|G26680|G26680"))
        self.assertNotIn("$#=0", stderrdata)

    def test_Mafft_with_Clustalw_output(self):
        """Simple round-trip through app with clustal output."""
        cmdline = MafftCommandline(mafft_exe)
        # Use some properties:
        cmdline.input = self.infile1
        cmdline.clustalout = True
        self.assertEqual(str(eval(repr(cmdline))), str(cmdline))
        stdoutdata, stderrdata = cmdline()
        # e.g. "CLUSTAL format alignment by MAFFT ..."
        # or "CLUSTAL (-like) formatted alignment by MAFFT FFT-NS-2 (v6.240)"
        self.assertTrue(stdoutdata.startswith("CLUSTAL"), stdoutdata)
        self.assertNotIn("$#=0", stderrdata)

    if version_major >= 7:
        def test_Mafft_with_PHYLIP_output(self):
            """Simple round-trip through app with PHYLIP output."""
            cmdline = MafftCommandline(mafft_exe, input=self.infile1,
                                       phylipout=True)
            self.assertEqual(str(eval(repr(cmdline))), str(cmdline))
            stdoutdata, stderrdata = cmdline()
            # e.g. " 3 706\n" or " 3 681" but allow some variation in the column count
            self.assertTrue(stdoutdata.startswith(" 3 68") or
                            stdoutdata.startswith(" 3 69") or
                            stdoutdata.startswith(" 3 70"), stdoutdata)
            self.assertTrue("gi|1348912 " in stdoutdata,
                            stdoutdata)
            self.assertTrue("gi|1348912|gb|G26680|G26680" not in stdoutdata,
                            stdoutdata)
            self.assertNotIn("$#=0", stderrdata)

        def test_Mafft_with_PHYLIP_namelength(self):
            """Check PHYLIP with --namelength."""
            cmdline = MafftCommandline(mafft_exe, input=self.infile1,
                                       phylipout=True, namelength=50)
            self.assertEqual(str(eval(repr(cmdline))), str(cmdline))
            stdoutdata, stderrdata = cmdline()
            # e.g. " 3 706\n" or " 3 681" but allow some variation in the column count
            self.assertTrue(stdoutdata.startswith(" 3 68") or
                            stdoutdata.startswith(" 3 69") or
                            stdoutdata.startswith(" 3 70"), stdoutdata)
            self.assertTrue("gi|1348912|gb|G26680|G26680" in stdoutdata,
                            stdoutdata)
            self.assertNotIn("$#=0", stderrdata)

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
        self.assertEqual(str(cmdline), mafft_exe +
                         " --localpair --weighti 4.2 --retree 5 " +
                         "--maxiterate 200 --nofft --op 2.04 --ep 0.51" +
                         " --lop 0.233 --lep 0.2 --reorder --treeout" +
                         " --nuc Fasta/f002")
        stdoutdata, stderrdata = cmdline()
        self.assertTrue(stdoutdata.startswith(">gi|1348912|gb|G26680|G26680"))
        self.assertNotIn("$#=0", stderrdata)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
