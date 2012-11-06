# Copyright 2008-2011 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

#TODO - Clean up the extra files created by clustalw?  e.g. *.dnd
#and *.aln where we have not requested an explicit name?

import sys
import os
import unittest
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.Application import ApplicationError

class ClustalWTestCase(unittest.TestCase):
    def setUp(self):
        if sys.platform=="win32":
            #TODO - Check the path?
            try:
                #This can vary depending on the Windows language.
                prog_files = os.environ["PROGRAMFILES"]
            except KeyError:
                prog_files = r"C:\Program Files"

            #Note that EBI's clustalw2 installer, e.g. clustalw-2.0.10-win.msi
            #uses C:\Program Files\ClustalW2\clustalw2.exe so we should check
            #for that.
            #
            #Some users doing a manual install have reported using
            #C:\Program Files\clustalw.exe
            #
            #Older installers might use something like this,
            #C:\Program Files\Clustalw\clustalw.exe
            #
            #One particular case is www.tc.cornell.edu currently provide a
            #clustalw1.83 installer which uses the following long location:
            #C:\Program Files\CTCBioApps\clustalw\v1.83\clustalw1.83.exe
            likely_dirs = ["ClustalW2", "",
                           "Clustal","Clustalw","Clustalw183","Clustalw1.83",
                           r"CTCBioApps\clustalw\v1.83"]
            likely_exes = ["clustalw2.exe",
                           "clustalw.exe", "clustalw1.83.exe"]
            for folder in likely_dirs:
                if os.path.isdir(os.path.join(prog_files, folder)):
                    for filename in likely_exes:
                        if os.path.isfile(os.path.join(prog_files, folder, filename)):
                            self.clustalw_exe = os.path.join(prog_files, folder, filename)
                            break
                    if self.clustalw_exe : break
        else:
            import commands
            #Note that clustalw 1.83 and clustalw 2.1 don't obey the --version
            #command, but this does cause them to quit cleanly.  Otherwise they prompt
            #the user for input (causing a lock up).
            output = commands.getoutput("clustalw2 --version")
            #Since "not found" may be in another language, try and be sure this is
            #really the clustalw tool's output
            if "not found" not in output and "CLUSTAL" in output \
            and "Multiple Sequence Alignments" in output:
                self.clustalw_exe = "clustalw2"
            if not self.clustalw_exe:
                output = commands.getoutput("clustalw --version")
                if "not found" not in output and "CLUSTAL" in output \
                and "Multiple Sequence Alignments" in output:
                    self.clustalw_exe = "clustalw"

        if not self.clustalw_exe:
            raise MissingExternalDependencyError(\
                "Install clustalw or clustalw2 if you want to use it from Biopython.")


class ClustalWTestErrorConditions(ClustalWTestCase):
    """Test general error conditions."""

    def test_empty_file(self):
        """Non-existing input file"""
        input_file = "does_not_exist.fasta"
        self.assertFalse(os.path.isfile(input_file))
        cline = ClustalwCommandline(self.clustalw_exe, infile=input_file)
	with self.assertRaises(ApplicationError) as cm:
            stdout, stderr = cline()
        err = cm.exception
        self.assertTrue("Cannot open sequence file" in str(err) or \
                        "Cannot open input file" in str(err) or \
                        "non-zero exit status" in str(err))

    def test_single_sequence(self):
        """Single sequence"""
        input_file = "Fasta/f001"
        self.assertTrue(os.path.isfile(input_file))
        self.assertTrue(len(list(SeqIO.parse(input_file, "fasta"))) == 1)
        cline = ClustalwCommandline(self.clustalw_exe, infile=input_file)
        with self.assertRaises(ApplicationError) as cm:
            stdout, stderr = cline()
            # Apparently some versions of clustal have a zero as return code
            # on error.  The following abuses ApplicationError to catch
            # both cases.
            if "cannot do multiple alignment" in (stdout + stderr):
                raise ApplicationError(-1, "", "", "cannot do multiple alignment")
        err = cm.exception
        self.assertTrue("No records found in handle" in str(err) or \
                        "cannot do multiple alignment" in str(err))

    def test_invalid_sequence(self):
        """Invalid sequence"""
        input_file = "Medline/pubmed_result1.txt"
        self.assertTrue(os.path.isfile(input_file))
        cline = ClustalwCommandline(self.clustalw_exe, infile=input_file)
        with self.assertRaises(ApplicationError) as cm:
            stdout, stderr = cline()
        err = cm.exception
        #Ideally we'd catch the return code and raise the specific
        #error for "invalid format", rather than just notice there
        #is not output file.
        #Note:
        #Python 2.3 on Windows gave (0, 'Error')
        #Python 2.5 on Windows gives [Errno 0] Error
        self.assertTrue("invalid format" in str(err) or \
                   "not produced" in str(err) or \
                   "No sequences in file" in str(err) or\
                   "non-zero exit status " in str(err))


class ClustalWTestNormalConditions(ClustalWTestCase):
    """Tests for normal situations."""

    def tearDown(self):
        if os.path.isfile("Fasta/f001.aln"):
            os.remove("Fasta/f001.aln")
        if os.path.isfile("Medline/pubmed_result1.aln"):
            os.remove("Medline/pubmed_result1.aln")
        if os.path.isfile(self.temp_filename_with_spaces):
            os.remove(self.temp_filename_with_spaces)
        if os.path.isfile(self.temp_large_fasta_file):
            os.remove(self.temp_large_fasta_file)

    def test_simple_fasta(self):
        input_file = "Fasta/f002"
        output_file = "temp_test.aln"
        statistics_file = "temp_stats.txt"

    def test_multiple_alignment(self):
        #Create a temp fasta file with a space in the name
        self.temp_filename_with_spaces = "Clustalw/temp horses.fasta"
        handle = open(self.temp_filename_with_spaces, "w")
        SeqIO.write(SeqIO.parse("Phylip/hennigian.phy","phylip"), handle, "fasta")
        handle.close()

        #Create a large input file by converting another example file
        #(See Bug 2804, this will produce so much output on stdout that
        #subprocess could suffer a deadlock and hang).  Using all the
        #records should show the deadlock but is very slow - just thirty
        #seems to lockup on Mac OS X, even 20 on Linux (without the fix).
        self.temp_large_fasta_file = "temp_cw_prot.fasta"
        handle = open(self.temp_large_fasta_file, "w")
        records = list(SeqIO.parse("NBRF/Cw_prot.pir", "pir"))[:40]
        SeqIO.write(records, handle, "fasta")
        handle.close()
        del handle, records

        for input_file, output_file, statistics_file, newtree_file in [
            ("Fasta/f002", "temp_test.aln", "temp_stats.txt", None),
            ("GFF/multi.fna", "temp with space.aln", "temp_stats.txt", None),
            ("Registry/seqs.fasta", "temp_test.aln", "temp_stats.txt", None),
            ("Registry/seqs.fasta", "temp_test.aln", "temp stats with space.txt", "temp_test.dnd"),
            ("Registry/seqs.fasta", "temp_test.aln", "temp_stats.txt", "temp with space.dnd"),
            (self.temp_filename_with_spaces, "temp_test.aln", "temp_stats.txt", None),
            (self.temp_filename_with_spaces, "temp with space.aln", "temp_stats", None),
            (self.temp_large_fasta_file, "temp_cw_prot.aln", "temp_stats.txt", None),
        ]:
            #Note that ClustalW will map ":" to "_" in it's output file
            input_records = SeqIO.to_dict(SeqIO.parse(input_file,"fasta"),
                                          lambda rec : rec.id.replace(":","_"))
            if os.path.isfile(output_file):
                os.remove(output_file)
            print "Calling clustalw on %s (with %i records)" \
                   % (repr(input_file), len(input_records))
            print "using output file %s" % repr(output_file)
            if newtree_file is not None:
                print "requesting output guide tree file %s" % repr(newtree_file)

            #Any filenames with spaces should get escaped with quotes automatically.
            #Using keyword arguments here.
            if self.clustalw_exe == "clustalw2":
                # By using the stats keyword, we require ClustalW 2.0.10 or higher.
                cline = ClustalwCommandline(self.clustalw_exe,
                                            infile=input_file,
                                            outfile=output_file,
                                            stats=statistics_file)
            else:
                cline = ClustalwCommandline(self.clustalw_exe,
                                            infile=input_file,
                                            outfile=output_file)
            self.assertTrue(str(eval(repr(cline))) == str(cline))
            if newtree_file is not None:
                #Test using a property:
                cline.newtree = newtree_file
                #I don't just want the tree, also want the alignment:
                cline.align = True
                self.assertTrue(str(eval(repr(cline))) == str(cline))
            output, error = cline()
            self.assertTrue(output.strip().startswith("CLUSTAL"))
            self.assertTrue(error.strip() == "")

            #Check the output...
            align = AlignIO.read(output_file, "clustal")
            #The length of the alignment will depend on the version of clustalw
            #(clustalw 2.0.10 and clustalw 1.83 are certainly different).
            print "Got an alignment, %i sequences" % (len(align))
            if self.clustalw_exe == "clustalw2":
                self.assertTrue(os.path.isfile(statistics_file))
                os.remove(statistics_file)
            output_records = SeqIO.to_dict(SeqIO.parse(output_file,"clustal"))
            self.assertTrue(set(input_records.keys()) == set(output_records.keys()))
            for record in align:
                self.assertTrue(str(record.seq) == str(output_records[record.id].seq))
                self.assertTrue(str(record.seq).replace("-","") == \
                       str(input_records[record.id].seq))

            #Clean up...
            os.remove(output_file)

            #Check the DND file was created.
            #TODO - Try and parse this with Bio.Nexus?
            if newtree_file is not None:
                tree_file = newtree_file
            else:
                #Clustalw will name it based on the input file
                tree_file = os.path.splitext(input_file)[0] + ".dnd"
            self.assertTrue(os.path.isfile(tree_file))
            os.remove(tree_file)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
