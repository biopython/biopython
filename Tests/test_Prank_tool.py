"""
Unittests for Bio.Align.Applications interface for PRANK

This code is part of the Biopython distribution and governed by its
license.  Please see the LICENSE file that should have been included
as part of this package.
"""

import sys
import os
import unittest
import subprocess
from Bio import AlignIO
from Bio import SeqIO
from Bio import MissingExternalDependencyError
from Bio.Align.Applications import PrankCommandline
from Bio.Nexus.Nexus import NexusError

#Try to avoid problems when the OS is in another language
os.environ['LANG'] = 'C'

prank_exe = None
if sys.platform=="win32":
    try:
        #This can vary depending on the Windows language.
        prog_files = os.environ["PROGRAMFILES"]
    except KeyError:
        prog_files = r"C:\Program Files"
    #For Windows, PRANK just comes as a zip file which contains the
    #prank.exe file which the user could put anywhere.  We'll try a few
    #sensible locations under Program Files... and then the full path.
    likely_dirs = ["", #Current dir
                   prog_files,
                   os.path.join(prog_files,"Prank")] + sys.path
    for folder in likely_dirs:
        if os.path.isdir(folder):
            if os.path.isfile(os.path.join(folder, "prank.exe")):
                prank_exe = os.path.join(folder, "prank.exe")
                break
        if prank_exe : break
else:
    import commands
    output = commands.getoutput("prank")
    if "not found" not in output and "prank" in output.lower():
        prank_exe = "prank"
if not prank_exe:
    raise MissingExternalDependencyError(\
        "Install PRANK if you want to use the Bio.Align.Applications wrapper.")

class PrankApplication(unittest.TestCase):
    
    def setUp(self):
        self.infile1 = "Fasta/fa01"

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
        self.assertEqual(str(cmdline), prank_exe + " -d=Fasta/fa01")
        self.assertEqual(str(eval(repr(cmdline))), str(cmdline))
        output, error = cmdline()
        self.assertEqual(error, "")
        self.assertTrue("Total time" in output)

    def test_Prank_simple_with_NEXUS_output(self):
        """Simple round-trip through app with infile, output in NEXUS
        output.?.??? files written to cwd - no way to redirect
        """
        records = list(SeqIO.parse(self.infile1,"fasta"))
        #Try using keyword argument,
        cmdline = PrankCommandline(prank_exe, d=self.infile1, noxml=True)
        #Try using a property,
        cmdline.d = self.infile1
        cmdline.f = 17 # NEXUS format
        cmdline.set_parameter("notree", True)
        self.assertEqual(str(cmdline), prank_exe + \
                         " -d=Fasta/fa01 -f=17 -noxml -notree")
        self.assertEqual(str(eval(repr(cmdline))), str(cmdline))
        stdout, stderr = cmdline()
        self.assertTrue("Total time" in stdout)
        self.assertEqual(stderr, "")
        try:
            align = AlignIO.read("output.2.nex", "nexus")
            for old, new in zip(records, align):
                #Old versions of Prank reduced name to 9 chars
                self.assertTrue(old.id==new.id or old.id[:9]==new.id)
                #infile1 has alignment gaps in it
                self.assertEqual(str(new.seq).replace("-",""),
                                 str(old.seq).replace("-",""))
        except NexusError:
            #See bug 3119,
            #Bio.Nexus can't parse output from prank v100701 (1 July 2010)
            pass

    def test_Prank_complex_command_line(self):
        """Round-trip with complex command line."""
        cmdline = PrankCommandline(prank_exe)
        cmdline.set_parameter("d", self.infile1)
        cmdline.set_parameter("-noxml", True)
        cmdline.set_parameter("notree", True)
        cmdline.set_parameter("-gaprate", 0.321)
        cmdline.set_parameter("gapext", 0.6)
        cmdline.set_parameter("-dots", 1) #i.e. True
        #Try using a property:
        cmdline.kappa = 3
        cmdline.skipins = True
        cmdline.set_parameter("-once", True)
        cmdline.realbranches = True
        self.assertEqual(str(cmdline), prank_exe + " -d=Fasta/fa01 -noxml" + \
                         " -notree -dots -gaprate=0.321 -gapext=0.6 -kappa=3" + \
                         " -once -skipins -realbranches")
        self.assertEqual(str(eval(repr(cmdline))), str(cmdline))
        stdout, stderr = cmdline()
        self.assertTrue("Total time" in stdout, stdout)


class PrankConversion(unittest.TestCase):
    def setUp(self):
        #As these reads are all 36, it can be seen as pre-aligned:
        self.input = "Quality/example.fasta"
        self.output = 'temp with space' #prefix, PRANK will pick extensions

    def conversion(self, prank_number, prank_ext, format):
        """Get PRANK to do a conversion, and check it with SeqIO."""
        filename = "%s.%s" % (self.output, prank_ext)
        if os.path.isfile(filename):
            os.remove(filename)
        cmdline = PrankCommandline(prank_exe, d=self.input,
                                   convert=True, f=prank_number,
                                   o='"%s"' % self.output)
        self.assertEqual(str(cmdline), prank_exe \
                         + ' -d=%s' % self.input \
                         + ' -o="%s"' % self.output \
                         + ' -f=%i' % prank_number \
                         + ' -convert')
        self.assertEqual(str(eval(repr(cmdline))), str(cmdline))
        message, error = cmdline()
        self.assertTrue(("PRANK: converting '%s' to '%s'" % (self.input, filename)) \
                        in message, message)
        self.assertEqual(error, "")
        self.assertTrue(os.path.isfile(filename))
        old = AlignIO.read(self.input, "fasta")
        #Hack...
        if format=="phylip":
            for record in old:
                record.id = record.id[:10]
        new = AlignIO.read(filename, format)
        assert len(old) == len(new)
        for old_r, new_r in zip(old, new):
            self.assertEqual(old_r.id, new_r.id)
            self.assertEqual(str(old_r.seq), str(new_r.seq))
        os.remove(filename)
        
    def test_convert_to_fasta(self):
        """Convert FASTA to FASTA format."""
        self.conversion(8, "fas", "fasta")

    #Prank v.100701 seems to output an invalid file here...
    #def test_convert_to_phylip32(self):
    #    """Convert FASTA to PHYLIP 3.2 format."""
    #    self.conversion(11, "phy", "phylip")

    def test_convert_to_phylip(self):
        """Convert FASTA to PHYLIP format."""
        self.conversion(12, "phy", "phylip")

    #PRANK truncated the record names in the matrix block. An error?
    #def test_convert_to_paup_nexus(self):
    #    """Convert FASTA to PAUP/NEXUS."""
    #    self.conversion(17, "nex", "nexus")

    #We don't support format 18, PAML


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
