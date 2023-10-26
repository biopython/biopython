# Copyright (C) 2012 by Eric Talevich.
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for Bio.Phylo.Applications wrappers."""

import sys
import os
import unittest

from subprocess import getoutput

from Bio import Phylo
from Bio.Phylo.Applications import PhymlCommandline
from Bio import MissingExternalDependencyError

# Try to avoid problems when the OS is in another language
os.environ["LANG"] = "C"

phyml_exe = None
exe_name = "PhyML-3.1_win32.exe" if sys.platform == "win32" else "phyml"

output = getoutput(exe_name + " --version")
# Looks like this:
# . This is PhyML version 20120412.
if "20" in output and "PhyML" in output:
    phyml_exe = exe_name

if not phyml_exe:
    raise MissingExternalDependencyError(
        "Couldn't find the PhyML software. Install PhyML 3.0 or later if you want "
        "to use the Bio.Phylo.Applications wrapper."
    )


# Example Phylip file with 4 aligned protein sequences
EX_PHYLIP = "Phylip/interlaced2.phy"


class AppTests(unittest.TestCase):
    """Tests for application wrappers."""

    def test_phyml(self):
        """Run PhyML using the wrapper."""
        # Stabilize phyml tests by running in single threaded mode by default.
        # Note: PHYMLCPUS environment is specific to Debian and derivatives.
        if not os.getenv("PHYMLCPUS"):
            os.putenv("PHYMLCPUS", "1")
        cmd = PhymlCommandline(phyml_exe, input=EX_PHYLIP, datatype="aa")
        # Smoke test
        try:
            out, err = cmd()
            self.assertGreater(len(out), 0)
            self.assertEqual(len(err), 0)
            # Check the output tree
            outfname = EX_PHYLIP + "_phyml_tree.txt"
            if not os.path.isfile(outfname):
                # NB: Briefly, PhyML dropped the .txt suffix (#919)
                outfname = outfname[:-4]
            tree = Phylo.read(outfname, "newick")
            self.assertEqual(tree.count_terminals(), 4)
        except Exception as exc:
            self.fail(f"PhyML wrapper error: {exc}")
        finally:
            # Clean up generated files
            for suffix in [
                "_phyml_tree.txt",
                "_phyml_tree",
                "_phyml_stats.txt",
                "_phyml_stats",
            ]:
                fname = EX_PHYLIP + suffix
                if os.path.isfile(fname):
                    os.remove(fname)


# ---------------------------------------------------------

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
