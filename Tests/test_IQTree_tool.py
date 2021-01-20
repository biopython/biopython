# Copyright 2021 by Cristian Ferrari.
# Based on test_phyml_tool.py by Eric Talevich.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Test for IQTree tool."""

from Bio import MissingExternalDependencyError

import sys
import os
import unittest

from subprocess import getoutput

from Bio import Phylo
from Bio.Phylo.Applications import IQTreeCommandline
from Bio.Application import ApplicationError

# Try to aboid problems when the OS is in another language
os.environ["LANG"] = "C"

iqtree_exe = None

if sys.platform == "Win32":

    try:
        prog_files = r"C:\Program Files\iqtree-1.6.12-Windows\bin"
    except KeyError:
        prog_files = r"C:\Program Files"


    likely_dirs = ["",
        "iqtree-1.6.12-Windows",
        r"iqtree-1.6.12-Windows\bin",
        "iqtree2-1.6.12-Windows",
        r"iqtree2-1.6.12-Windows\bin",]
    likely_exes = ["iqtree.exe", "iqtree2.exe"]

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
            iqtree_exe = "iqtree.exe"
    if not iqtree_exe:
        output = getoutput("iqtree2 -version")
        if "not found" not in output and "not recognized" not in output:
            if "IQ-TREE" in output:
                iqtree_exe = "iqtree2.exe"

if not iqtree_exe:
    raise MissingExternalDependencyError(
       "Install IQTree or IQTree2 if you want to use the Bio.Phylo.Applications wrapper."
    )

#Example Phylip file provided on IQ-Tree website
Example_Phy = "Phylip/example_iqtree.phy"

class IQTreeTestCase(unittest.TestCase):
    """Tests for application wrappers"""

    def test_iqtree(self):
        cmd = _IQtree.IQTreeCommandline(iqtree_exe, s=Example_Phy)

        try:
             out, err = cmd()
             self.assertGreater(len(out), 0)
             self.assertEqual(len(err), 0)
             # Check the output tree
             outfname = Example_Phy + ".iqtree"
             print(outfname)
             self.assertEqual(os.path.isfile(outfname), True)

        except Exception as exc:
             self.fail("IQTree wrapper error: %s" % exc)

        finally:
          # Clean up generated files
          for suffix in [".bionj",
              ".ckp.gz",
              ".iqtree",
              ".log",
              ".mldist",
              ".model.gz",
              ".treefile",]:
              fname = Example_Phy + suffix
              if os.path.isfile(fname):
                  os.remove(fname)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)