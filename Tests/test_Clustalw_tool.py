# Copyright 2008-2011 by Peter Cock.  All rights reserved.
# Revisions copyright 2012 by Christian Brueffer.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# TODO - Clean up the extra files created by clustalw?  e.g. *.dnd
# and *.aln where we have not requested an explicit name?
from __future__ import print_function

from Bio import MissingExternalDependencyError

import sys
import os
import unittest
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.Application import ApplicationError

#################################################################

# Try to avoid problems when the OS is in another language
os.environ['LANG'] = 'C'

clustalw_exe = None
if sys.platform == "win32":
    # TODO - Check the path?
    try:
        # This can vary depending on the Windows language.
        prog_files = os.environ["PROGRAMFILES"]
    except KeyError:
        prog_files = r"C:\Program Files"

    # Note that EBI's clustalw2 installer, e.g. clustalw-2.0.10-win.msi
    # uses C:\Program Files\ClustalW2\clustalw2.exe so we should check
    # for that.
    #
    # Some users doing a manual install have reported using
    # C:\Program Files\clustalw.exe
    #
    # Older installers might use something like this,
    # C:\Program Files\Clustalw\clustalw.exe
    #
    # One particular case is www.tc.cornell.edu currently provide a
    # clustalw1.83 installer which uses the following long location:
    # C:\Program Files\CTCBioApps\clustalw\v1.83\clustalw1.83.exe
    likely_dirs = ["ClustalW2", "",
                   "Clustal", "Clustalw", "Clustalw183", "Clustalw1.83",
                   r"CTCBioApps\clustalw\v1.83"]
    likely_exes = ["clustalw2.exe",
                   "clustalw.exe", "clustalw1.83.exe"]
    for folder in likely_dirs:
        if os.path.isdir(os.path.join(prog_files, folder)):
            for filename in likely_exes:
                if os.path.isfile(os.path.join(prog_files, folder, filename)):
                    clustalw_exe = os.path.join(prog_files, folder, filename)
                    break
            if clustalw_exe:
                break
else:
    from Bio._py3k import getoutput
    # Note that clustalw 1.83 and clustalw 2.1 don't obey the --version
    # command, but this does cause them to quit cleanly.  Otherwise they prompt
    # the user for input (causing a lock up).
    output = getoutput("clustalw2 --version")
    # Since "not found" may be in another language, try and be sure this is
    # really the clustalw tool's output
    if "not found" not in output and "CLUSTAL" in output \
    and "Multiple Sequence Alignments" in output:
        clustalw_exe = "clustalw2"
    if not clustalw_exe:
        output = getoutput("clustalw --version")
        if "not found" not in output and "CLUSTAL" in output \
        and "Multiple Sequence Alignments" in output:
            clustalw_exe = "clustalw"

if not clustalw_exe:
    raise MissingExternalDependencyError(
        "Install clustalw or clustalw2 if you want to use it from Biopython.")


class ClustalWTestCase(unittest.TestCase):
    """Class implementing common functions for ClustalW tests."""

    def setUp(self):
        self.files_to_clean = set()

    def tearDown(self):
        for filename in self.files_to_clean:
            if os.path.isfile(filename):
                os.remove(filename)

    def standard_test_procedure(self, cline):
        """Standard testing procedure used by all tests."""
        self.assertTrue(str(eval(repr(cline))) == str(cline))
        input_records = SeqIO.to_dict(SeqIO.parse(cline.infile, "fasta"),
                                      lambda rec: rec.id.replace(":", "_"))

        # Determine name of tree file
        if cline.newtree:
            tree_file = cline.newtree
        else:
            # Clustalw will name it based on the input file
            tree_file = os.path.splitext(cline.infile)[0] + ".dnd"

        # Mark generated files for later removal
        self.add_file_to_clean(cline.outfile)
        self.add_file_to_clean(tree_file)

        output, error = cline()
        self.assertTrue(output.strip().startswith("CLUSTAL"))
        self.assertTrue(error.strip() == "")

        # Check the output...
        align = AlignIO.read(cline.outfile, "clustal")
        # The length of the alignment will depend on the version of clustalw
        # (clustalw 2.1 and clustalw 1.83 are certainly different).
        output_records = SeqIO.to_dict(SeqIO.parse(cline.outfile, "clustal"))
        self.assertTrue(set(input_records.keys()) == set(output_records.keys()))
        for record in align:
            self.assertTrue(str(record.seq) == str(output_records[record.id].seq))
            self.assertTrue(str(record.seq).replace("-", "") ==
                   str(input_records[record.id].seq))

        # Check the DND file was created.
        # TODO - Try and parse this with Bio.Nexus?
        self.assertTrue(os.path.isfile(tree_file))

    def add_file_to_clean(self, filename):
        """Adds a file for deferred removal by the tearDown routine."""
        self.files_to_clean.add(filename)


class ClustalWTestErrorConditions(ClustalWTestCase):
    """Test general error conditions."""

    def test_empty_file(self):
        """Test a non-existing input file."""
        input_file = "does_not_exist.fasta"
        self.assertFalse(os.path.isfile(input_file))
        cline = ClustalwCommandline(clustalw_exe, infile=input_file)

        try:
            stdout, stderr = cline()
        except ApplicationError as err:
            self.assertTrue("Cannot open sequence file" in str(err) or
                            "Cannot open input file" in str(err) or
                            "Non-zero return code " in str(err), str(err))
        else:
            self.fail("expected an ApplicationError")

    def test_single_sequence(self):
        """Test an input file containing a single sequence."""
        input_file = "Fasta/f001"
        self.assertTrue(os.path.isfile(input_file))
        self.assertTrue(len(list(SeqIO.parse(input_file, "fasta"))) == 1)
        cline = ClustalwCommandline(clustalw_exe, infile=input_file)

        try:
            stdout, stderr = cline()
            # Zero return code is a possible bug in clustalw 2.1?
            self.assertIn("cannot do multiple alignment", (stdout + stderr))
        except ApplicationError as err:
            # Good, non-zero return code indicating an error in clustalw
            # e.g. Using clustalw 1.83 get:
            # Command 'clustalw -infile=Fasta/f001' returned non-zero exit status 4
            pass

        if os.path.isfile(input_file + ".aln"):
            # Clustalw 2.1 made an emtpy aln file, clustalw 1.83 did not
            self.add_file_to_clean(input_file + ".aln")

    def test_invalid_sequence(self):
        """Test an input file containing an invalid sequence."""
        input_file = "Medline/pubmed_result1.txt"
        self.assertTrue(os.path.isfile(input_file))
        cline = ClustalwCommandline(clustalw_exe, infile=input_file)

        try:
            stdout, stderr = cline()
        except ApplicationError as err:
            # Ideally we'd catch the return code and raise the specific
            # error for "invalid format", rather than just notice there
            # is not output file.
            # Note:
            # Python 2.3 on Windows gave (0, 'Error')
            # Python 2.5 on Windows gives [Errno 0] Error
            self.assertTrue("invalid format" in str(err) or
                            "not produced" in str(err) or
                            "No sequences in file" in str(err) or
                            "Non-zero return code " in str(err))
        else:
            self.fail("expected an ApplicationError")


class ClustalWTestNormalConditions(ClustalWTestCase):
    """Tests for normal conditions."""

    def test_properties(self):
        """Test passing options via properties."""
        cline = ClustalwCommandline(clustalw_exe)
        cline.infile = "Fasta/f002"
        cline.outfile = "temp_test.aln"
        cline.align = True

        self.standard_test_procedure(cline)

    def test_simple_fasta(self):
        """Test a simple fasta input file."""
        input_file = "Fasta/f002"
        output_file = "temp_test.aln"
        cline = ClustalwCommandline(clustalw_exe,
                                    infile=input_file,
                                    outfile=output_file)

        self.standard_test_procedure(cline)

    def test_newtree(self):
        """Test newtree files."""
        input_file = "Registry/seqs.fasta"
        output_file = "temp_test.aln"
        newtree_file = "temp_test.dnd"
        cline = ClustalwCommandline(clustalw_exe,
                                    infile=input_file,
                                    outfile=output_file,
                                    newtree=newtree_file,
                                    align=True)

        self.standard_test_procedure(cline)
        cline.newtree = "temp with space.dnd"
        self.standard_test_procedure(cline)

    def test_large_input_file(self):
        """Test a large input file."""

        # Create a large input file by converting another example file
        # (See Bug 2804, this will produce so much output on stdout that
        # subprocess could suffer a deadlock and hang).  Using all the
        # records should show the deadlock but is very slow - just thirty
        # seems to lockup on Mac OS X, even 20 on Linux (without the fix).
        input_file = "temp_cw_prot.fasta"
        handle = open(input_file, "w")
        records = list(SeqIO.parse("NBRF/Cw_prot.pir", "pir"))[:40]
        SeqIO.write(records, handle, "fasta")
        handle.close()
        del handle, records
        output_file = "temp_cw_prot.aln"

        cline = ClustalwCommandline(clustalw_exe,
                                    infile=input_file,
                                    outfile=output_file)

        self.add_file_to_clean(input_file)
        self.standard_test_procedure(cline)

    def test_input_filename_with_space(self):
        """Test an input filename containing a space."""
        input_file = "Clustalw/temp horses.fasta"
        handle = open(input_file, "w")
        SeqIO.write(SeqIO.parse("Phylip/hennigian.phy", "phylip"), handle, "fasta")
        handle.close()
        output_file = "temp with space.aln"

        cline = ClustalwCommandline(clustalw_exe,
                                    infile=input_file,
                                    outfile=output_file)

        self.add_file_to_clean(input_file)
        self.standard_test_procedure(cline)

    def test_output_filename_with_spaces(self):
        """Test an output filename containing spaces."""
        input_file = "GFF/multi.fna"
        output_file = "temp with space.aln"
        cline = ClustalwCommandline(clustalw_exe,
                                    infile=input_file,
                                    outfile=output_file)

        self.standard_test_procedure(cline)


class ClustalWTestVersionTwoSpecific(ClustalWTestCase):
    """Tests specific to ClustalW2."""

    def test_statistics(self):
        """Test a statistics file."""
        if clustalw_exe == "clustalw2":
            input_file = "Fasta/f002"
            output_file = "temp_test.aln"
            statistics_file = "temp_stats.txt"
            cline = ClustalwCommandline(clustalw_exe,
                                        infile=input_file,
                                        outfile=output_file,
                                        stats=statistics_file)

            self.add_file_to_clean(statistics_file)
            self.standard_test_procedure(cline)
            self.assertTrue(os.path.isfile(statistics_file))
        else:
            print("Skipping ClustalW2 specific test.")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
