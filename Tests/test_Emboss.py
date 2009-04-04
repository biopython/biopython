# Copyright 2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Runs a few EMBOSS tools to check our wrappers and parsers."""

import os
import sys
import unittest

from Bio.Application import generic_run
from Bio.Emboss import Applications
from Bio import SeqIO
from Bio import AlignIO
from Bio import MissingExternalDependencyError

try :
    import subprocess
except ImportError :
    raise MissingExternalDependencyError(\
        "Python 2.3 not supported, this needs the subprocess module.")

#################################################################

exes_wanted = ["water", "needle"]
exes = dict() #Dictionary mapping from names to exe locations
if sys.platform=="win32" :
    #TODO - Find out where the default install goes, and/or
    #if we can use the registry to find it.
    raise MissingExternalDependencyError(\
        "Auto-detection of EMBOSS on Windows not supported (yet).")
else :
    import commands
    for name in exes_wanted :
        #This will "just work" if installed on the path as normal on Unix
        output = commands.getoutput("%s -help" % name)
        if "not found" not in output :
            exes[name] = name

if len(exes) < len(exes_wanted) :
    raise MissingExternalDependencyError(\
        "Install EMBOSS if you want to use Bio.EMBOSS.")

#################################################################

class PairwiseAlignmentTests(unittest.TestCase):
    """Run pairwise alignments with water and needle, and parse them."""
    def test_water_file(self):
        """Call water with the asis trick, output to a file."""
        #Setup,
        cline = Applications.WaterCommandline(cmd=exes["water"])
        cline.set_parameter("-asequence", "asis:ACCCGGGCGCGGT")
        cline.set_parameter("-bsequence", "asis:ACCCGAGCGCGGT")
        cline.set_parameter("-gapopen", "10")
        cline.set_parameter("-gapextend", "0.5")
        cline.set_parameter("-outfile", "temp_test.water")
        #Run the tool,
        result, out, err = generic_run(cline)
        #Check it worked,
        self.assertEqual(result.return_code, 0)
        self.assertEqual(out.read().strip(), "")
        self.assertEqual(err.read().strip(), "Smith-Waterman local alignment.")
        filename = result.get_result("-outfile")
        self.assertEqual(filename, "temp_test.water")
        assert os.path.isfile(filename)
        #Check we can parse the output...
        align = AlignIO.read(open(filename),"emboss")
        self.assertEqual(len(align), 2)
        self.assertEqual(str(align[0].seq), "ACCCGGGCGCGGT")
        self.assertEqual(str(align[1].seq), "ACCCGAGCGCGGT")

    def test_water_piped(self):
        """Call water with asis trick, output piped to stdout."""
        #TODO - Support -auto and -filter in Bio.Emboss.Applications
        cline = exes["water"]
        cline += " -asequence asis:ACCCGGGCGCGGT"
        cline += " -bsequence asis:ACCCGAGCGCGGT"
        cline += " -auto" #no prompting
        cline += " -filter" #use stdout
        #Run the tool,
        child = subprocess.Popen(str(cline),
                                 stdin=subprocess.PIPE,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 shell=(sys.platform!="win32"))
        child.stdin.close()
        #Check we could read it's output
        align = AlignIO.read(child.stdout, "emboss")
        self.assertEqual(len(align), 2)
        self.assertEqual(str(align[0].seq), "ACCCGGGCGCGGT")
        self.assertEqual(str(align[1].seq), "ACCCGAGCGCGGT")
        #Check no error output:
        assert child.stderr.read() == ""
        assert 0 == child.wait()

    def test_needle_piped(self):
        """Call needle with asis trick, output piped to stdout."""
        #TODO - Support needle in Bio.Emboss.Applications
        #(ideally with the -auto and -filter arguments)
        #Setup,
        cline = exes["needle"]
        cline += " -asequence asis:ACCCGGGCGCGGT"
        cline += " -bsequence asis:ACCCGAGCGCGGT"
        cline += " -auto" #no prompting
        cline += " -filter" #use stdout
        #Run the tool,
        child = subprocess.Popen(str(cline),
                                 stdin=subprocess.PIPE,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 shell=(sys.platform!="win32"))
        child.stdin.close()
        #Check we could read it's output
        align = AlignIO.read(child.stdout, "emboss")
        self.assertEqual(len(align), 2)
        self.assertEqual(str(align[0].seq), "ACCCGGGCGCGGT")
        self.assertEqual(str(align[1].seq), "ACCCGAGCGCGGT")
        #Check no error output:
        assert child.stderr.read() == ""
        assert 0 == child.wait()

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
