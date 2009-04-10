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
from Bio.Alphabet import generic_protein, generic_dna

try :
    import subprocess
except ImportError :
    raise MissingExternalDependencyError(\
        "Python 2.3 not supported, this needs the subprocess module.")

#################################################################

exes_wanted = ["water", "needle", "seqret"]
exes = dict() #Dictionary mapping from names to exe locations
if sys.platform=="win32" :
    #The default installation path is C:\mEMBOSS which contains the exes.
    #EMBOSS also sets an environment variable which we will check for.
    try :
        path = os.environ["EMBOSS_ROOT"]
    except KeyError :
        #print >> sys.stderr, "Missing EMBOSS_ROOT environment variable!"
        raise MissingExternalDependencyError(\
            "Install EMBOSS if you want to use Bio.EMBOSS.")
    if os.path.isdir(path) :
        for name in exes_wanted :
            if os.path.isfile(os.path.join(path, name+".exe")) :
                exes[name] = os.path.join(path, name+".exe")
    del path, name
else :
    import commands
    for name in exes_wanted :
        #This will "just work" if installed on the path as normal on Unix
        if "not found" not in commands.getoutput("%s -help" % name) :
            exes[name] = name
    del name

if len(exes) < len(exes_wanted) :
    raise MissingExternalDependencyError(\
        "Install EMBOSS if you want to use Bio.EMBOSS.")

#################################################################

#Top level function as this makes it easier to use for debugging:
def emboss_convert(filename, old_format, new_format):
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

def compare_records(old_list, new_list) :
    """Check two lists of SeqRecords agree, raises a ValueError if mismatch."""
    if len(old_list) != len(new_list) :
        raise ValueError("%i vs %i records" % (len(old_list), len(new_list)))
    for old, new in zip(old_list, new_list) :
        #TODO - check annotation
        if old.id != new.id and old.name != new.name \
        and not new.id.endswith(old.id) and not old.id.endswith(new.id) :
            raise ValueError("'%s' or '%s' vs '%s' or '%s' records" \
                             % (old.id, old.name, new.id, new.name))

        if str(old.seq).upper() != str(new.seq).upper() :
            raise ValueError("'%s' vs '%s'" % (old.seq, new.seq))
        if old.features and new.features \
        and len(old.features) != len(new.features) :
            raise ValueError("%i vs %i features" \
                             % (len(old.features, len(new.features))))
    return True

def compare_alignments(old_list, new_list) :
    """Check two lists of Alignments agree, raises a ValueError if mismatch."""
    if len(old_list) != len(new_list) :
        raise ValueError("%i vs %i alignments" % (len(old_list), len(new_list)))
    for old, new in zip(old_list, new_list) :
        if len(old) != len(new) :
            raise ValueError("Alignment with %i vs %i records" \
                             % (len(old), len(new)))
        for r1, r2 in zip(old, new) :
            if str(r1.seq).upper() != str(r2.seq).upper() :
                raise ValueError("'%s' vs '%s'" % (r1.seq, r2.seq))
            #TODO - Be flexible on the names, PHYLIP truncates them.
    return True

class SeqRetTests(unittest.TestCase):
    """Check EMBOSS seqret against Bio.SeqIO for converting files."""

    def check_SeqIO_to_EMBOSS(self, in_filename, in_format, skip_formats=[],
                              alphabet=None) :
        """Can Bio.SeqIO write files seqret can read back?"""
        if alphabet :
            records = list(SeqIO.parse(open(in_filename), in_format, alphabet))
        else :
            records = list(SeqIO.parse(open(in_filename), in_format))
        for temp_format in ["genbank","fasta"] :
            if temp_format in skip_formats :
                continue
            #TODO - Handle this with a pipe?
            #i.e. can Bio.SeqIO write to the stdin of seqret?
            filename = "Emboss/temp_%s.txt" % temp_format
            temp_handle = open(filename,"w")
            SeqIO.write(records, temp_handle, temp_format)
            temp_handle.flush()
            temp_handle.close()
            
            handle = emboss_convert(filename, temp_format, "fasta")
            new_records = list(SeqIO.parse(handle, "fasta"))

            try :
                self.assert_(compare_records(records, new_records))
            except ValueError, err :
                raise ValueError("Disagree on file %s %s in %s format: %s" \
                                 % (in_format, in_filename, temp_format, err))
            os.remove(filename)
            
    def check_EMBOSS_to_SeqIO(self, filename, old_format,
                              skip_formats=[]) :
        """Can Bio.SeqIO read seqret's conversion of the file?"""
        #TODO: Why can't we read EMBOSS's swiss output?
        self.assert_(os.path.isfile(filename))
        old_records = list(SeqIO.parse(open(filename), old_format))
        for new_format in ["genbank","fasta","pir","embl", "ig"] :
            if new_format in skip_formats :
                continue
            handle = emboss_convert(filename, old_format, new_format)
            new_records = list(SeqIO.parse(handle, new_format))
            try :
                self.assert_(compare_records(old_records, new_records))
            except ValueError, err:
                raise ValueError("Disagree on %s file %s in %s format: %s" \
                                 % (old_format, filename, new_format, err))

    def check_SeqIO_with_EMBOSS(self, filename, old_format, skip_formats=[],
                                alphabet=None):
        #TODO - May need to supply the sequence alphabet for parsing.
        #Check EMBOSS can read Bio.SeqIO output...
        self.check_SeqIO_to_EMBOSS(filename, old_format, skip_formats,
                                   alphabet)
        #Check Bio.SeqIO can read EMBOSS seqret output...
        self.check_EMBOSS_to_SeqIO(filename, old_format, skip_formats)

    def test_genbank(self) :
        """SeqIO & EMBOSS reading each other's conversions of a GenBank file."""
        self.check_SeqIO_with_EMBOSS("GenBank/cor6_6.gb", "genbank")

    def test_embl(self) :
        """SeqIO & EMBOSS reading each other's conversions of an EMBL file."""
        self.check_SeqIO_with_EMBOSS("EMBL/U87107.embl", "embl")

    def test_ig(self) :
        """SeqIO & EMBOSS reading each other's conversions of an ig file."""
        self.check_SeqIO_to_EMBOSS("IntelliGenetics/VIF_mase-pro.txt", "ig",
                                   alphabet=generic_protein)
        #TODO - What does a % in an ig sequence mean?
        #e.g. "IntelliGenetics/vpu_nucaligned.txt"
        #and  "IntelliGenetics/TAT_mase_nuc.txt"
        #EMBOSS seems to ignore them.

    def test_pir(self) :
        """SeqIO & EMBOSS reading each other's conversions of a PIR file."""
        #Skip genbank here, EMBOSS mangles the LOCUS line:
        self.check_SeqIO_with_EMBOSS("NBRF/clustalw.pir", "pir",
                               skip_formats=["genbank"])
        #Skip EMBL here, EMBOSS mangles the ID line
        #Skip GenBank, EMBOSS 6.0.1 on Windows won't output proteins as GenBank
        self.check_SeqIO_with_EMBOSS("NBRF/DMB_prot.pir", "pir",
                               skip_formats=["embl","genbank"])
    def test_clustalw(self) :
        """SeqIO & EMBOSS reading each other's conversions of a Clustalw file."""
        self.check_SeqIO_with_EMBOSS("Clustalw/hedgehog.aln", "clustal",
                                   skip_formats=["embl","genbank"])
        self.check_SeqIO_with_EMBOSS("Clustalw/opuntia.aln", "clustal",
                                   skip_formats=["embl","genbank"])

    def check_EMBOSS_to_AlignIO(self, filename, old_format,
                              skip_formats=[]) :
        """Can AlignIO read seqret's conversion of the file?"""
        self.assert_(os.path.isfile(filename), filename)
        old_aligns = list(AlignIO.parse(open(filename), old_format))
        formats = ["clustal", "phylip", "ig"]
        if len(old_aligns) == 1 :
            formats.extend(["fasta","nexus"])
        for new_format in formats :
            if new_format in skip_formats :
                continue
            handle = emboss_convert(filename, old_format, new_format)
            try :
                new_aligns = list(AlignIO.parse(handle, new_format))
            except :
                raise ValueError("Can't parse %s file %s in %s format." \
                                 % (old_format, filename, new_format))
            try :
                self.assert_(compare_alignments(old_aligns, new_aligns))
            except ValueError, err :
                raise ValueError("Disagree on %s file %s in %s format: %s" \
                                 % (old_format, filename, new_format, err))

    def test_align_clustalw(self) :
        """Can AlignIO read EMBOSS's conversions of a Clustalw file?"""
        self.check_EMBOSS_to_AlignIO("Clustalw/hedgehog.aln", "clustal")
        self.check_EMBOSS_to_AlignIO("Clustalw/opuntia.aln", "clustal")
        self.check_EMBOSS_to_AlignIO("Clustalw/odd_consensus.aln", "clustal",
                               skip_formats=["nexus"]) #TODO - why not nexus?
        self.check_EMBOSS_to_AlignIO("Clustalw/protein.aln", "clustal")
        self.check_EMBOSS_to_AlignIO("Clustalw/promals3d.aln", "clustal")

    def test_clustalw(self) :
        """Can we read EMBOSS's conversions of a PHYLIP file?"""
        self.check_EMBOSS_to_AlignIO("Phylip/horses.phy", "phylip")
        self.check_EMBOSS_to_AlignIO("Phylip/hennigian.phy", "phylip")
        self.check_EMBOSS_to_AlignIO("Phylip/reference_dna.phy", "phylip")
        self.check_EMBOSS_to_AlignIO("Phylip/reference_dna2.phy", "phylip")
        self.check_EMBOSS_to_AlignIO("Phylip/interlaced.phy", "phylip")
        self.check_EMBOSS_to_AlignIO("Phylip/interlaced2.phy", "phylip")
        self.check_EMBOSS_to_AlignIO("Phylip/random.phy", "phylip")
        
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
        errors = err.read().strip()
        self.assert_(errors.startswith("Smith-Waterman local alignment"), errors)
        self.assertEqual(out.read().strip(), "")
        if result.return_code != 0 : print >> sys.stderr, "\n%s"%cline
        self.assertEqual(result.return_code, 0)
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
        errors = err.read().strip()
        self.assert_(errors.startswith("Smith-Waterman local alignment"), errors)
        self.assertEqual(out.read().strip(), "")
        if result.return_code != 0 : print >> sys.stderr, "\n%s"%cline
        self.assertEqual(result.return_code, 0)
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
        errors = err.read().strip()
        self.assert_(errors.startswith("Smith-Waterman local alignment"), errors)
        self.assertEqual(out.read().strip(), "")
        if result.return_code != 0 : print >> sys.stderr, "\n%s"%cline
        self.assertEqual(result.return_code, 0)
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
        errors = err.read().strip()
        self.assert_(errors.startswith("Smith-Waterman local alignment"), errors)
        self.assertEqual(out.read().strip(), "")
        if result.return_code != 0 : print >> sys.stderr, "\n%s"%cline
        self.assertEqual(result.return_code, 0)
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
