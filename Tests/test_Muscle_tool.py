# Copyright 2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
import sys
import subprocess
import unittest

from Bio import MissingExternalDependencyError
from Bio.Align.Applications import MuscleCommandline
from Bio import SeqIO
from Bio import AlignIO

#################################################################

#Try to avoid problems when the OS is in another language
os.environ['LANG'] = 'C'

muscle_exe = None
if sys.platform=="win32":
    try:
        #This can vary depending on the Windows language.
        prog_files = os.environ["PROGRAMFILES"]
    except KeyError:
        prog_files = r"C:\Program Files"
    #For Windows, MUSCLE just comes as a zip file which contains the
    #a Muscle directory with the muscle.exe file plus a readme etc,
    #which the user could put anywhere.  We'll try a few sensible
    #locations under Program Files... and then the full path.
    likely_dirs = ["", #Current dir
                   prog_files,
                   os.path.join(prog_files,"Muscle3.6"),
                   os.path.join(prog_files,"Muscle3.7"),
                   os.path.join(prog_files,"Muscle3.8"),
                   os.path.join(prog_files,"Muscle3.9"),
                   os.path.join(prog_files,"Muscle")] + sys.path
    for folder in likely_dirs:
        if os.path.isdir(folder):
            if os.path.isfile(os.path.join(folder, "muscle.exe")):
                muscle_exe = os.path.join(folder, "muscle.exe")
                break
        if muscle_exe : break
else:
    import commands
    output = commands.getoutput("muscle -version")
    #Since "not found" may be in another language, try and be sure this is
    #really the MUSCLE tool's output
    if "not found" not in output and "MUSCLE" in output \
    and "Edgar" in output:
        muscle_exe = "muscle"

if not muscle_exe:
    raise MissingExternalDependencyError(\
        "Install MUSCLE if you want to use the Bio.Align.Applications wrapper.")

#################################################################

class MuscleApplication(unittest.TestCase):
    
    def setUp(self):
        self.infile1  = "Fasta/f002"
        self.infile2  = "Fasta/fa01"
        self.infile3  = "Fasta/f001"
        self.outfile1 = "Fasta/temp align out1.fa" #with spaces!
        self.outfile2 = "Fasta/temp_align_out2.fa"
        self.outfile3 = "Fasta/temp_align_out3.fa"
        self.outfile4 = "Fasta/temp_align_out4.fa"

    def tearDown(self):
        if os.path.isfile(self.outfile1):
            os.remove(self.outfile1)
        if os.path.isfile(self.outfile2):
            os.remove(self.outfile2)
        if os.path.isfile(self.outfile3):
            os.remove(self.outfile3)
        if os.path.isfile(self.outfile4):
            os.remove(self.outfile4)

    def test_Muscle_simple(self):
        """Simple round-trip through app just infile and outfile"""
        cmdline = MuscleCommandline(muscle_exe,
                                    input=self.infile1,
                                    out=self.outfile1)
        self.assertEqual(str(cmdline), muscle_exe \
                         + ' -in Fasta/f002 -out "Fasta/temp align out1.fa"')
        self.assertEqual(str(eval(repr(cmdline))), str(cmdline))
        output, error = cmdline()
        self.assertEqual(output, "")
        self.assertTrue("ERROR" not in error)

    def test_Muscle_with_options(self):
        """Round-trip through app with a switch and valued option"""
        cmdline = MuscleCommandline(muscle_exe)
        cmdline.set_parameter("input", self.infile1) #"input" is alias for "in"
        cmdline.set_parameter("out", self.outfile2)
        #Use property:
        cmdline.objscore = "sp"
        cmdline.noanchors = True
        self.assertEqual(str(cmdline), muscle_exe +\
                         " -in Fasta/f002" + \
                         " -out Fasta/temp_align_out2.fa" + \
                         " -objscore sp -noanchors")
        self.assertEqual(str(eval(repr(cmdline))), str(cmdline))
        output, error = cmdline()
        self.assertEqual(output, "")
        self.assertTrue("ERROR" not in error)
        self.assertTrue(error.strip().startswith("MUSCLE"), output)

    def test_Muscle_profile_simple(self):
        """Simple round-trip through app doing a profile alignment"""
        cmdline = MuscleCommandline(muscle_exe)
        cmdline.set_parameter("out", self.outfile3)
        cmdline.set_parameter("profile", True)
        cmdline.set_parameter("in1", self.infile2)
        cmdline.set_parameter("in2", self.infile3)
        self.assertEqual(str(cmdline), muscle_exe + \
                         " -out Fasta/temp_align_out3.fa" + \
                         " -profile -in1 Fasta/fa01 -in2 Fasta/f001")
        self.assertEqual(str(eval(repr(cmdline))), str(cmdline))
        output, error = cmdline()
        self.assertEqual(output, "")
        self.assertTrue("ERROR" not in error)
        self.assertTrue(error.strip().startswith("MUSCLE"), output)

    def test_Muscle_profile_with_options(self):
        """Profile alignment, and switch and valued options"""
        #Using some keyword arguments, note -stable isn't supported in v3.8
        cmdline = MuscleCommandline(muscle_exe, out=self.outfile4,
                                    in1=self.infile2, in2=self.infile3,
                                    profile=True, stable=True,
                                    cluster1="neighborjoining")
        self.assertEqual(str(cmdline), muscle_exe + \
                         " -out Fasta/temp_align_out4.fa" + \
                         " -profile -in1 Fasta/fa01 -in2 Fasta/f001" + \
                         " -cluster1 neighborjoining -stable")
        self.assertEqual(str(eval(repr(cmdline))), str(cmdline))
        """
        #TODO - Why doesn't this work with MUSCLE 3.6 on the Mac?
        #It may be another bug fixed in MUSCLE 3.7 ...
        result, stdout, stderr = generic_run(cmdline)
        #NOTE: generic_run has been removed from Biopython
        self.assertEqual(result.return_code, 0)
        self.assertEqual(stdout.read(), "")
        self.assertTrue("ERROR" not in stderr.read())
        self.assertEqual(str(result._cl), str(cmdline))
        """

class SimpleAlignTest(unittest.TestCase):
    """Simple MUSCLE tests"""

    """
    #FASTA output seems broken on Muscle 3.6 (on the Mac).
    def test_simple_fasta(self):
        input_file = "Fasta/f002"
        self.assertTrue(os.path.isfile(input_file))
        records = list(SeqIO.parse(input_file,"fasta"))
        #Prepare the command...
        cmdline = MuscleCommandline(muscle_exe)
        cmdline.set_parameter("in", input_file)
        #Preserve input record order (makes checking output easier)
        cmdline.set_parameter("stable")
        #Set some others options just to test them
        cmdline.set_parameter("maxiters", 2)
        self.assertEqual(str(cmdline).rstrip(), "muscle -in Fasta/f002 -maxiters 2 -stable")
        result, out_handle, err_handle = generic_run(cmdline)
        #NOTE: generic_run has been removed from Biopython
        print err_handle.read()
        print out_handle.read()
        align = AlignIO.read(out_handle, "fasta")
        self.assertEqual(len(records),len(align))
        for old, new in zip(records, align):
            self.assertEqual(old.id, new.id)
            self.assertEqual(str(new.seq).replace("-",""), str(old.seq))
    """

    def test_simple_clustal(self):
        """Simple muscle call using Clustal output with a MUSCLE header"""
        input_file = "Fasta/f002"
        self.assertTrue(os.path.isfile(input_file))
        records = list(SeqIO.parse(input_file,"fasta"))
        records.sort(key = lambda rec: rec.id)
        #Prepare the command... use Clustal output (with a MUSCLE header)
        cmdline = MuscleCommandline(muscle_exe, input=input_file, clw = True)
        self.assertEqual(str(cmdline).rstrip(), muscle_exe + \
                         " -in Fasta/f002 -clw")
        self.assertEqual(str(eval(repr(cmdline))), str(cmdline))
        child = subprocess.Popen(str(cmdline),
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True,
                                 shell=(sys.platform!="win32"))
        #Didn't use -quiet so there should be progress reports on stderr,
        align = AlignIO.read(child.stdout, "clustal")
        align.sort() #by record.id
        self.assertTrue(child.stderr.read().strip().startswith("MUSCLE"))
        return_code = child.wait()
        self.assertEqual(return_code, 0)
        child.stdout.close()
        child.stderr.close()
        del child
        self.assertEqual(len(records),len(align))
        for old, new in zip(records, align):
            self.assertEqual(old.id, new.id)
            self.assertEqual(str(new.seq).replace("-",""), str(old.seq))

    def test_simple_clustal_strict(self):
        """Simple muscle call using strict Clustal output"""
        input_file = "Fasta/f002"
        self.assertTrue(os.path.isfile(input_file))
        records = list(SeqIO.parse(input_file,"fasta"))
        records.sort(key = lambda rec: rec.id)
        #Prepare the command...
        cmdline = MuscleCommandline(muscle_exe)
        cmdline.set_parameter("in", input_file)
        #Use clustal output (with a CLUSTAL header)
        cmdline.set_parameter("clwstrict", True) #Default None treated as False!
        self.assertEqual(str(cmdline).rstrip(), muscle_exe + \
                         " -in Fasta/f002 -clwstrict")
        self.assertEqual(str(eval(repr(cmdline))), str(cmdline))
        child = subprocess.Popen(str(cmdline),
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True,
                                 shell=(sys.platform!="win32"))
        #Didn't use -quiet so there should be progress reports on stderr,
        align = AlignIO.read(child.stdout, "clustal")
        align.sort()
        self.assertTrue(child.stderr.read().strip().startswith("MUSCLE"))
        self.assertEqual(len(records),len(align))
        for old, new in zip(records, align):
            self.assertEqual(old.id, new.id)
            self.assertEqual(str(new.seq).replace("-",""), str(old.seq))
        return_code = child.wait()
        self.assertEqual(return_code, 0)
        child.stdout.close()
        child.stderr.close()
        del child

    def test_long(self):
        """Simple muscle call using long file"""
        #Create a large input file by converting some of another example file
        temp_large_fasta_file = "temp_cw_prot.fasta"
        records = list(SeqIO.parse("NBRF/Cw_prot.pir", "pir"))[:40]
        SeqIO.write(records, temp_large_fasta_file, "fasta")
        #Prepare the command...
        cmdline = MuscleCommandline(muscle_exe)
        cmdline.set_parameter("in", temp_large_fasta_file)
        #Use fast options
        cmdline.set_parameter("maxiters", 1)
        cmdline.set_parameter("diags", True) #Default None treated as False!
        #Use clustal output
        cmdline.set_parameter("clwstrict", True) #Default None treated as False!
        #Shoudn't need this, but just to make sure it is accepted
        cmdline.set_parameter("maxhours", 0.1)
        #No progress reports to stderr
        cmdline.set_parameter("quiet", True) #Default None treated as False!
        self.assertEqual(str(cmdline).rstrip(), muscle_exe + \
                         " -in temp_cw_prot.fasta -diags -maxhours 0.1" + \
                         " -maxiters 1 -clwstrict -quiet")
        self.assertEqual(str(eval(repr(cmdline))), str(cmdline))
        child = subprocess.Popen(str(cmdline),
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True,
                                 shell=(sys.platform!="win32"))
        align = AlignIO.read(child.stdout, "clustal")
        align.sort()
        records.sort(key = lambda rec: rec.id)
        self.assertEqual(len(records), len(align))
        for old, new in zip(records, align):
            self.assertEqual(old.id, new.id)
            self.assertEqual(str(new.seq).replace("-",""), str(old.seq))
        os.remove(temp_large_fasta_file)
        #See if quiet worked:
        self.assertEqual("", child.stderr.read().strip())
        return_code = child.wait()
        self.assertEqual(return_code, 0)
        child.stdout.close()
        child.stderr.close()
        del child

    def test_using_stdin(self):
        """Simple alignment using stdin"""
        input_file = "Fasta/f002"
        self.assertTrue(os.path.isfile(input_file))
        records = list(SeqIO.parse(input_file,"fasta"))
        #Prepare the command... use Clustal output (with a MUSCLE header)
        cline = MuscleCommandline(muscle_exe, clw=True)
        self.assertEqual(str(cline).rstrip(), muscle_exe + " -clw")
        self.assertEqual(str(eval(repr(cline))), str(cline))
        child = subprocess.Popen(str(cline),
                                 stdin=subprocess.PIPE,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True,
                                 shell=(sys.platform!="win32"))
        SeqIO.write(records, child.stdin, "fasta")
        child.stdin.close()
        #Alignment will now run...
        align = AlignIO.read(child.stdout, "clustal")
        align.sort()
        records.sort(key = lambda rec: rec.id)
        self.assertEqual(len(records),len(align))
        for old, new in zip(records, align):
            self.assertEqual(old.id, new.id)
            self.assertEqual(str(new.seq).replace("-",""), str(old.seq))
        self.assertEqual(0, child.wait())
        child.stdout.close()
        child.stderr.close()
        del child

    def test_with_multiple_output_formats(self):
        """Simple muscle call with multiple output formats"""
        input_file = "Fasta/f002"
        output_html = "temp_f002.html"
        output_clwstrict = "temp_f002.clw"
        self.assertTrue(os.path.isfile(input_file))
        records = list(SeqIO.parse(input_file,"fasta"))
        records.sort(key = lambda rec: rec.id)
        #Prepare the command... use Clustal output (with a MUSCLE header)
        cmdline = MuscleCommandline(muscle_exe, input=input_file,
                                    clw=True, htmlout = output_html,
                                    clwstrictout = output_clwstrict)
        self.assertEqual(str(cmdline).rstrip(), muscle_exe + \
                         " -in Fasta/f002 -clw -htmlout temp_f002.html" +\
                         " -clwstrictout temp_f002.clw")
        self.assertEqual(str(eval(repr(cmdline))), str(cmdline))
        child = subprocess.Popen(str(cmdline),
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True,
                                 shell=(sys.platform!="win32"))
        #Clustalw on stdout:
        align = AlignIO.read(child.stdout, "clustal")
        align.sort()
        #Didn't use -quiet so there should be progress reports on stderr,
        self.assertTrue(child.stderr.read().strip().startswith("MUSCLE"))
        return_code = child.wait()
        self.assertEqual(return_code, 0)
        self.assertEqual(len(records),len(align))
        for old, new in zip(records, align):
            self.assertEqual(old.id, new.id)
        child.stdout.close()
        child.stderr.close()
        del child
        handle = open(output_html,"rU")
        html = handle.read().strip().upper()
        handle.close()
        self.assertTrue(html.startswith("<HTML"))
        self.assertTrue(html.endswith("</HTML>"))
        #ClustalW strict:
        align = AlignIO.read(output_clwstrict, "clustal")
        align.sort()
        self.assertEqual(len(records),len(align))
        for old, new in zip(records, align):
            self.assertEqual(old.id, new.id)
        os.remove(output_html)
        os.remove(output_clwstrict)

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
