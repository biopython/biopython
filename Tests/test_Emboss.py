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

exes_wanted = ["water", "needle", "seqret"]
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

class SeqRetTests(unittest.TestCase):
    """Check EMBOSS seqret against Bio.SeqIO for converting files."""
    def emboss_convert(self, filename, old_format, new_format):
        """Run seqret, returns handle."""
        #TODO - Support seqret in Bio.Emboss.Applications
        #(ideally with the -auto and -filter arguments)
        #Setup, this assumes for all the format names used
        #Biopython and EMBOSS names are consistent!
        cline = exes["seqret"]
        cline += " -sequence " + filename
        cline += " -sformat " + old_format
        cline += " -osformat " + new_format
        cline += " -auto" #no prompting
        cline += " -filter" #use stdout
        #Run the tool,
        child = subprocess.Popen(str(cline),
                                 stdin=subprocess.PIPE,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 shell=(sys.platform!="win32"))
        child.stdin.close()
        return child.stdout

    def compare_records(self, old_list, new_list) :
        """Check two lists of SeqRecords agree."""
        self.assertEqual(len(old_list), len(new_list))
        for old, new in zip(old_list, new_list) :
            self.assertEqual(str(old.seq).upper(), str(new.seq).upper())
        return True
        
    def check_can_read_emboss_conversion(self, filename, old_format, new_format) :
        """Can Bio.SeqIO read seqret's conversion of the file?"""
        self.assert_(os.path.isfile(filename))
        old_records = list(SeqIO.parse(open(filename), old_format))
        handle = self.emboss_convert(filename, old_format, new_format)
        new_records = list(SeqIO.parse(handle, new_format))
        self.compare_records(old_records, new_records)

    def test_genbank(self) :
        """Can we read EMBOSS's conversions of a GenBank file?"""
        self.check_can_read_emboss_conversion("GenBank/cor6_6.gb", "genbank", "genbank")
        self.check_can_read_emboss_conversion("GenBank/cor6_6.gb", "genbank", "fasta")
        self.check_can_read_emboss_conversion("GenBank/cor6_6.gb", "genbank", "pir")
        self.check_can_read_emboss_conversion("GenBank/cor6_6.gb", "genbank", "embl")
        #TODO: Why can't we read EMBOSS's swiss output?
        #self.check_can_read_emboss_conversion("GenBank/cor6_6.gb", "genbank", "swiss")
        
    
        
class PairwiseAlignmentTests(unittest.TestCase):
    """Run pairwise alignments with water and needle, and parse them."""
    def pairwise_alignment_check(self, query_seq,
                                 targets, alignments,
                                 local=True) :
        """Check pairwise alignment data is sane."""
        #The datasets should be small, so making iterators into lists is OK
        targets = list(targets)
        alignments = list(alignments)
        self.assertEqual(len(targets), len(alignments))
        for target, alignment in zip(targets, alignments) :
            self.assertEqual(len(alignment), 2)
            #self.assertEqual(target.id, alignment[1].id) #too strict
            if alignment[1].id not in target.id \
            and alignment[1].id not in target.name :
                raise AssertionError("%s vs %s or %s" \
                                     % (alignment[1].id , target.id, target.name))
            if local :
                #Local alignment
                self.assert_(str(alignment[0].seq).replace("-","") \
                             in query_seq)
                self.assert_(str(alignment[1].seq).replace("-","").upper() \
                             in str(target.seq).upper())
            else :
                #Global alignment
                self.assertEqual(str(query_seq), str(alignment[0].seq).replace("-",""))
                self.assertEqual(str(target.seq).upper(), \
                                 str(alignment[1].seq).replace("-","").upper())
        return True

    def test_water_file(self):
        """water with the asis trick, output to a file."""
        #Setup,
        cline = Applications.WaterCommandline(cmd=exes["water"])
        cline.set_parameter("-asequence", "asis:ACCCGGGCGCGGT")
        cline.set_parameter("-bsequence", "asis:ACCCGAGCGCGGT")
        cline.set_parameter("-gapopen", "10")
        cline.set_parameter("-gapextend", "0.5")
        cline.set_parameter("-outfile", "Emboss/temp_test.water")
        #Run the tool,
        result, out, err = generic_run(cline)
        #Check it worked,
        self.assertEqual(result.return_code, 0)
        self.assertEqual(out.read().strip(), "")
        self.assertEqual(err.read().strip(), "Smith-Waterman local alignment.")
        filename = result.get_result("-outfile")
        self.assertEqual(filename, "Emboss/temp_test.water")
        assert os.path.isfile(filename)
        #Check we can parse the output...
        align = AlignIO.read(open(filename),"emboss")
        self.assertEqual(len(align), 2)
        self.assertEqual(str(align[0].seq), "ACCCGGGCGCGGT")
        self.assertEqual(str(align[1].seq), "ACCCGAGCGCGGT")
        #Clean up,
        os.remove(filename)            
        
    def test_water_piped(self):
        """water with asis trick, output piped to stdout."""
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
        """needle with asis trick, output piped to stdout."""
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

    def test_water_file2(self):
        """water with the asis trick and nucleotide FASTA file, output to a file."""
        #Setup,
        query = "ACACACTCACACACACTTGGTCAGAGATGCTGTGCTTCTTGGAAGCAAGGNCTCAAAGGCAAGGTGCACGCAGAGGGACGTTTGAGTCTGGGATGAAGCATGTNCGTATTATTTATATGATGGAATTTCACGTTTTTATG"
        out_file = "Emboss/temp_test2.water"
        in_file = "Fasta/f002"
        self.assert_(os.path.isfile(in_file))
        if os.path.isfile(out_file) :
            os.remove(out_file)
        cline = Applications.WaterCommandline(cmd=exes["water"])
        cline.set_parameter("-asequence", "asis:%s" % query)
        cline.set_parameter("-bsequence", in_file)
        cline.set_parameter("-gapopen", "10")
        cline.set_parameter("-gapextend", "0.5")
        cline.set_parameter("-outfile", out_file)
        #Run the tool,
        result, out, err = generic_run(cline)
        #Check it worked,
        if result.return_code != 0 : print cline
        self.assertEqual(result.return_code, 0)
        self.assertEqual(out.read().strip(), "")
        self.assertEqual(err.read().strip(), "Smith-Waterman local alignment.")
        self.assertEqual(result.get_result("-outfile"), out_file)
        assert os.path.isfile(out_file)
        #Check we can parse the output and it is sensible...
        self.pairwise_alignment_check(query,
                                      SeqIO.parse(open(in_file),"fasta"),
                                      AlignIO.parse(open(out_file),"emboss"),
                                      local=True)
        #Clean up,
        os.remove(out_file)

    def test_water_file3(self):
        """water with the asis trick and GenBank file, output to a file."""
        #Setup,
        query = "TGTTGTAATGTTTTAATGTTTCTTCTCCCTTTAGATGTACTACGTTTGGA"
        out_file = "Emboss/temp_test3.water"
        in_file = "GenBank/cor6_6.gb"
        self.assert_(os.path.isfile(in_file))
        if os.path.isfile(out_file) :
            os.remove(out_file)
        cline = Applications.WaterCommandline(cmd=exes["water"])
        cline.set_parameter("-asequence", "asis:%s" % query)
        cline.set_parameter("-bsequence", in_file)
        #TODO - Tell water this is a GenBank file!
        cline.set_parameter("-gapopen", "1")
        cline.set_parameter("-gapextend", "0.5")
        cline.set_parameter("-outfile", out_file)
        #Run the tool,
        result, out, err = generic_run(cline)
        #Check it worked,
        if result.return_code != 0 : print cline
        self.assertEqual(result.return_code, 0)
        self.assertEqual(out.read().strip(), "")
        self.assertEqual(err.read().strip(), "Smith-Waterman local alignment.")
        self.assertEqual(result.get_result("-outfile"), out_file)
        assert os.path.isfile(out_file)
        #Check we can parse the output and it is sensible...
        self.pairwise_alignment_check(query,
                                      SeqIO.parse(open(in_file),"genbank"),
                                      AlignIO.parse(open(out_file),"emboss"),
                                      local=True)
        #Clean up,
        os.remove(out_file)

    def test_water_file4(self):
        """water with the asis trick and SwissProt file, output to a file."""
        #Setup,
        query = "DVCTGKALCDPVTQNIKTYPVKIENLRVMI"
        out_file = "Emboss/temp_test4.water"
        in_file = "SwissProt/sp004"
        self.assert_(os.path.isfile(in_file))
        if os.path.isfile(out_file) :
            os.remove(out_file)
        cline = Applications.WaterCommandline(cmd=exes["water"])
        cline.set_parameter("-asequence", "asis:%s" % query)
        cline.set_parameter("-bsequence", in_file)
        #TODO - Tell water this is a SwissProt file!
        cline.set_parameter("-gapopen", "20")
        cline.set_parameter("-gapextend", "5")
        cline.set_parameter("-outfile", out_file)
        #Run the tool,
        result, out, err = generic_run(cline)
        #Check it worked,
        if result.return_code != 0 : print cline
        self.assertEqual(result.return_code, 0)
        self.assertEqual(out.read().strip(), "")
        self.assertEqual(err.read().strip(), "Smith-Waterman local alignment.")
        self.assertEqual(result.get_result("-outfile"), out_file)
        assert os.path.isfile(out_file)
        #Check we can parse the output and it is sensible...
        self.pairwise_alignment_check(query,
                                      SeqIO.parse(open(in_file),"swiss"),
                                      AlignIO.parse(open(out_file),"emboss"),
                                      local=True)
        #Clean up,
        os.remove(out_file)
        
    def test_needle_piped2(self):
        """needle with asis trick, and nucleotide FASTA file, output piped to stdout."""
        #TODO - Support needle in Bio.Emboss.Applications
        #(ideally with the -auto and -filter arguments)
        #Setup,
        query = "ACACACTCACACACACTTGGTCAGAGATGCTGTGCTTCTTGGAA"
        cline = exes["needle"]
        cline += " -asequence asis:" + query
        cline += " -bsequence Fasta/f002"
        cline += " -auto" #no prompting
        cline += " -filter" #use stdout
        #Run the tool,
        child = subprocess.Popen(str(cline),
                                 stdin=subprocess.PIPE,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 shell=(sys.platform!="win32"))
        child.stdin.close()
        #Check we can parse the output and it is sensible...
        self.pairwise_alignment_check(query,
                                      SeqIO.parse(open("Fasta/f002"),"fasta"),
                                      AlignIO.parse(child.stdout,"emboss"),
                                      local=False)
        #Check no error output:
        assert child.stderr.read() == ""
        assert 0 == child.wait()
        
def clean_up() :
    """Fallback clean up method to remove temp files."""
    for filename in os.listdir("Emboss") :
        if filename.startswith("temp_") :
            try :
                os.remove(filename)
            except :
                pass

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
    clean_up()
