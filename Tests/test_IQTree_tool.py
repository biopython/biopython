
"""Test for IQTree tool."""

from Bio import MissingExternalDependencyError

import sys
import os
import unittest

from Bio import SeqIO
from Bio import AlignIO

from Bio.Align.Applications import IQTreeCommandline
from Bio.Application import ApplicationError

#################################################################

# Try to aboid problems when the OS is in another language
os.environ["LANG"] = "C"

iqtree_exe = None

if sys.platform == "Win32":

    try:
        prog_files = os.environ["PROGRAMFILES"]
    except KeyError:
        prog_files = r"C:\Program Files"


    likely_dirs = [
        "",
        "iqtree-1.6.12-Windows",
        r"iqtree-1.6.12-Windows\bin",
        ]
    likely_exes = ["iqtree.exe"]

    for folder in likely_dirs:
        if os.path.isdir(os.path.join(prog_files, folder)):
            for filename in likely_exes:
                if os.path.isfile(os.path.join(prog_files, folder, filename)):
                    iqtree_exe = os.path.join(prof_files, folder, filename)
                    break
                if iqtree_exe:
                    break
else:
    from subprocess import getoutput

    output = getoutput("iqtree -version")

    if "not found" not in output and "not recognized" not in output:
        if "IQ-TREE" in output:
            iqtree_exe = "iqtree"
    if not iqtree_exe:
        output = geoutput("iqtree2 -version")
        if "not found" not in output and "not recognized" not in output:
            if "IQ-TREE" in output:
                iqtree_exe = "iqtree2"

if not iqtree_exe:
    raise MissingExternalDependencyError(
        "Install IQTree or IQTree2 if you want to use the "
        "Bio.Phylo.Applications wrapper."
        )

#Example Phylip file provided on IQ-Tree website
Example_Phy = "example.phy"

class IQTreeTestCase(unittest.TestCase):
    """Tests for application wrappers"""
    cmd = IQTreeCommandline(iqtree_exe, s=Example_Phy)

    try:
        out, err = cmd()
        self.assertGreater(len(out), 0)
        self.assertEqual(len(err), 0)
        # Check the output tree
        outfname = Example_Phy + ".iqtree"
        if not os.path.isfile(outfname):
            outfname = outfname[:-4]

    except Exception as exc:
        self.fail("IQTree wrapper error: %s" % exc)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
