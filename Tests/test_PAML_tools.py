# Copyright (C) 2011 by Brandon Invergo (b.invergo@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for PAML tools module."""

import unittest
import os
import sys
from Bio.Phylo.PAML import codeml, baseml, yn00
from Bio import MissingExternalDependencyError


def is_exe(filepath):
    """Test if a file is an executable."""
    return os.path.exists(filepath) and os.access(filepath, os.X_OK)


def which(program):
    """Find the path to an executable."""
    filepath, filename = os.path.split(program)
    os_path = os.environ["PATH"].split(os.pathsep)
    if sys.platform == "win32":
        try:
            # This can vary depending on the Windows language.
            prog_files = os.environ["PROGRAMFILES"]
        except KeyError:
            prog_files = r"C:\Program Files"
        # For Windows, the user is instructed to move the programs to a folder
        # and then to add the folder to the system path. Just in case they didn't
        # do that, we can check for it in Program Files.
        likely_dirs = [
            "",  # Current dir
            prog_files,
            os.path.join(prog_files, "paml41"),
            os.path.join(prog_files, "paml43"),
            os.path.join(prog_files, "paml44"),
            os.path.join(prog_files, "paml45"),
        ] + sys.path
        os_path.extend(likely_dirs)
    for path in os.environ["PATH"].split(os.pathsep):
        exe_file = os.path.join(path, program)
        if is_exe(exe_file):
            return exe_file
    return None


# Find the PAML binaries
if sys.platform == "win32":
    binaries = ["codeml.exe", "baseml.exe", "yn00.exe"]
else:
    binaries = ["codeml", "baseml", "yn00"]
for binary in binaries:
    if which(binary) is None:
        raise MissingExternalDependencyError(
            "Install PAML if you want to use the Bio.Phylo.PAML wrapper."
        )


class Common(unittest.TestCase):
    """Base class for PAML unit tests."""

    del_files = []

    def __del__(self):
        """Just in case tool creates some junk files, do a clean-up."""
        del_files = self.del_files
        for filename in del_files:
            if os.path.exists(filename):
                os.remove(filename)


class CodemlTest(Common):
    """Tests for PAML tool codeml."""

    def setUp(self):
        self.cml = codeml.Codeml()

    def testCodemlBinary(self):
        """Check codeml runs, generates correct output, and is the correct version."""
        ctl_file = os.path.join("PAML", "Control_files", "codeml", "codeml.ctl")
        self.cml.read_ctl_file(ctl_file)
        self.cml.alignment = os.path.join("PAML", "Alignments", "alignment.phylip")
        self.cml.tree = os.path.join("PAML", "Trees", "species.tree")
        self.cml.out_file = os.path.join("PAML", "temp.out")
        self.cml.working_dir = os.path.join("PAML", "codeml_test")
        results = self.cml.run()
        self.assertGreater(results["version"], "4.0")
        self.assertIn("NSsites", results)
        self.assertEqual(len(results["NSsites"]), 1)
        self.assertEqual(len(results["NSsites"][0]), 5)


class BasemlTest(Common):
    """Tests for PAML tool baseml."""

    def setUp(self):
        self.bml = baseml.Baseml()

    def testBasemlBinary(self):
        """Check baseml runs, generates correct output, and is the correct version."""
        ctl_file = os.path.join("PAML", "Control_files", "baseml", "baseml.ctl")
        self.bml.read_ctl_file(ctl_file)
        self.bml.alignment = os.path.join("PAML", "Alignments", "alignment.phylip")
        self.bml.tree = os.path.join("PAML", "Trees", "species.tree")
        self.bml.out_file = os.path.join("PAML", "temp.out")
        self.bml.working_dir = os.path.join("PAML", "baseml_test")
        results = self.bml.run()
        self.assertGreater(results["version"], "4.0")
        self.assertIn("parameters", results)
        self.assertEqual(len(results["parameters"]), 5)


class Yn00Test(Common):
    """Tests for PAML tool yn00."""

    def setUp(self):
        self.yn = yn00.Yn00()

    def testYn00Binary(self):
        """Check yn00 binary runs and generates correct output.

        yn00 output does not specify the version number.
        """
        ctl_file = os.path.join("PAML", "Control_files", "yn00", "yn00.ctl")
        self.yn.read_ctl_file(ctl_file)
        self.yn.alignment = os.path.join("PAML", "Alignments", "alignment.phylip")
        self.yn.out_file = os.path.join("PAML", "temp.out")
        self.yn.working_dir = os.path.join("PAML", "yn00_test")
        results = self.yn.run()
        self.assertEqual(len(results), 5)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
