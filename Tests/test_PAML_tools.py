# Copyright (C) 2011 by Brandon Invergo (b.invergo@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

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
            #This can vary depending on the Windows language.
            prog_files = os.environ["PROGRAMFILES"]
        except KeyError:
            prog_files = r"C:\Program Files"
        #For Windows, the user is instructed to move the programs to a folder
        #and then to add the folder to the system path. Just in case they didn't
        #do that, we can check for it in Program Files.
        likely_dirs = ["", #Current dir
                       prog_files,
                       os.path.join(prog_files,"paml41"),
                       os.path.join(prog_files,"paml43"),
                       os.path.join(prog_files,"paml44")] + sys.path
        os_path.extend(likely_dirs)
    for path in os.environ["PATH"].split(os.pathsep):
        exe_file = os.path.join(path, program)
        if is_exe(exe_file):
            return exe_file
    return None

#Find the PAML binaries
if sys.platform == "win32":
    binaries = ["codeml.exe", "baseml.exe", "yn00.exe"]
else:
    binaries = ["codeml", "baseml", "yn00"]
for binary in binaries:
    if which(binary) is None:
        raise MissingExternalDependencyError(\
            "Install PAML if you want to use the Bio.Phylo.PAML wrapper.")

class CodemlTest(unittest.TestCase):

    align_file = os.path.join("PAML", "alignment.phylip")
    tree_file = os.path.join("PAML", "species.tree")
    out_file = os.path.join("PAML", "temp.out")
    working_dir = os.path.join("PAML", "codeml_test")
    ctl_file = os.path.join("PAML", "codeml.ctl")

    def __del__(self):
        """Just in case CODEML creates some junk files, do a clean-up."""
        del_files = [self.out_file, "2NG.dN", "2NG.dS", "2NG.t", "codeml.ctl",
                     "lnf", "rst", "rst1", "rub"]
        for filename in del_files:
            if os.path.exists(filename):
                os.remove(filename)
        if os.path.exists(self.working_dir):
            for filename in os.listdir(self.working_dir):
                os.remove(filename)
            os.rmdir(self.working_dir)

    def setUp(self):
        self.cml = codeml.Codeml(working_dir=self.working_dir)


    def testCodemlBinary(self):
        """Test that the codeml binary runs and generates correct output
        and is the correct version.
        """
        self.cml.read_ctl_file(self.ctl_file)
        results = self.cml.run()
        self.assertTrue(results["version"] > "4.0")
        self.assertTrue("NSsites" in results)
        self.assertEqual(len(results["NSsites"]), 1)
        self.assertEqual(len(results["NSsites"][0]), 5)

class BasemlTest(unittest.TestCase):

    align_file = os.path.join("PAML", "alignment.phylip")
    tree_file = os.path.join("PAML", "species.tree")
    out_file = os.path.join("PAML", "temp.out")
    working_dir = os.path.join("PAML", "baseml_test")
    ctl_file = os.path.join("PAML", "baseml.ctl")

    def __del__(self):
        """Just in case BASEML creates some junk files, do a clean-up."""
        del_files = [self.out_file, "2base.t", "in.basemlg", "baseml.ctl",
                     "lnf", "rates", "rst", "rst1", "rub", "temp.out"]
        for filename in del_files:
            if os.path.exists(filename):
                os.remove(filename)
        if os.path.exists(self.working_dir):
            for filename in os.listdir(self.working_dir):
                os.remove(filename)
            os.rmdir(self.working_dir)

    def setUp(self):
        self.bml = baseml.Baseml(working_dir=self.working_dir)


    def testBasemlBinary(self):
        """Test that the baseml binary runs and generates correct output
        and is the correct version.
        """
        self.bml.read_ctl_file(self.ctl_file)
        results = self.bml.run()
        self.assertTrue(results["version"] > "4.0")
        self.assertTrue("parameters" in results)
        self.assertEqual(len(results["parameters"]), 5)

class Yn00Test(unittest.TestCase):

    align_file = os.path.join("PAML", "alignment.phylip")
    out_file = os.path.join("PAML", "temp.out")
    working_dir = os.path.join("PAML", "yn00_test")
    ctl_file = os.path.join("PAML", "yn00.ctl")

    def __del__(self):
        """Just in case BASEML creates some junk files, do a clean-up."""
        del_files = [self.out_file, "2YN.dN", "2YN.dS", "2YN.t", "yn00.ctl",
                     "rst", "rst1", "rub"]
        for filename in del_files:
            if os.path.exists(filename):
                os.remove(filename)
        if os.path.exists(self.working_dir):
            for filename in os.listdir(self.working_dir):
                os.remove(filename)
            os.rmdir(self.working_dir)

    def setUp(self):
        self.yn = yn00.Yn00(working_dir=self.working_dir)


    def testYn00Binary(self):
        """Test that the yn00 binary runs and generates correct output.
        yn00 output does not specify the version number.
        """
        self.yn.read_ctl_file(self.ctl_file)
        results = self.yn.run()
        self.assertEqual(len(results), 5)
