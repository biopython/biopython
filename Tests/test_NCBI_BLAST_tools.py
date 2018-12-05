# Copyright 2009-2013 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
# This unit test attempts to locate the blastall executable and the nr
# database.

from __future__ import print_function

import os
import os.path
import sys
import subprocess
import unittest
import re

from Bio.Application import _escape_filename
from Bio import MissingExternalDependencyError
from Bio.Blast import Applications

# TODO - On windows, can we use the ncbi.ini file?
wanted = ["blastx", "blastp", "blastn", "tblastn", "tblastx",
          "rpsblast+",  # For Debian
          "rpsblast", "rpstblastn", "psiblast", "blast_formatter",
          "deltablast", "makeblastdb"]
exe_names = {}

if sys.platform == "win32":
    # The Windows 32 bit BLAST 2.2.22+ installer does add itself to the path,
    # and by default installs to C:\Program Files\NCBI\BLAST-2.2.22+\bin
    # To keep things simple, assume BLAST+ is on the path on Windows.
    #
    # On Windows the environment variable name isn't case sensitive,
    # but must split on ";" not ":"
    likely_dirs = os.environ.get("PATH", "").split(";")
else:
    likely_dirs = os.environ.get("PATH", "").split(":")

for folder in likely_dirs:
    if not os.path.isdir(folder):
        continue
    # Loop over copy as will remove entries from wanted:
    for name in wanted[:]:
        if sys.platform == "win32":
            exe_name = os.path.join(folder, name + ".exe")
        else:
            exe_name = os.path.join(folder, name)
        if not os.path.isfile(exe_name):
            continue
        # To tell the old and new rpsblast apart (since I have both on
        # my path and the old blast has priority), try -h as a parameter.
        # This should also reject WU-BLAST (since it doesn't like -h).
        child = subprocess.Popen(exe_name + " -h",
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True,
                                 shell=(sys.platform != "win32"))
        output, error = child.communicate()
        if child.returncode == 0 and "ERROR: Invalid argument: -h" not in output:
            # Special case, blast_formatter from BLAST 2.2.23+ (i.e. BLAST+)
            # has mandatory argument -rid, but no -archive. We don't support it.
            if name == "blast_formatter" and " -archive " not in output:
                continue
            exe_names[name] = exe_name
            wanted.remove(name)  # can stop search for this now
        # else:
        #    print("Rejecting %r" % exe_name)
        del exe_name, name

# To avoid the name clash with legacy BLAST, Debian introduced rpsblast+ alias
if "rpsblast+" in wanted:
    wanted.remove("rpsblast+")
if "rpsblast+" in exe_names:
    exe_names["rpsblast"] = exe_names["rpsblast+"]
    del exe_names["rpsblast+"]

# We can cope with blast_formatter being missing, only added in BLAST 2.2.24+
# We can cope with deltablast being missing, only added in BLAST 2.2.26+
optional = ["blast_formatter", "deltablast"]
if len(set(exe_names).difference(optional)) < len(set(wanted).difference(optional)):
    raise MissingExternalDependencyError("Install the NCBI BLAST+ command line "
                                         "tools if you want to use the "
                                         "Bio.Blast.Applications wrapper.")


class Pairwise(unittest.TestCase):
    def test_blastp(self):
        """Pairwise BLASTP search"""
        global exe_names
        cline = Applications.NcbiblastpCommandline(
            exe_names["blastp"],
            query="Fasta/rose.pro",
            subject="GenBank/NC_005816.faa",
            evalue=1)
        self.assertEqual(str(cline), _escape_filename(exe_names["blastp"]) +
                         " -query Fasta/rose.pro -evalue 1" +
                         " -subject GenBank/NC_005816.faa")
        child = subprocess.Popen(str(cline),
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True,
                                 shell=(sys.platform != "win32"))
        stdoutdata, stderrdata = child.communicate()
        return_code = child.returncode
        self.assertEqual(return_code, 0, "Got error code %i back from:\n%s"
                         % (return_code, cline))
        # Used to get 10 matches from 10 pairwise searches,
        # as of NCBI BLAST+ 2.3.0 only get 1 Query= line:
        if stdoutdata.count("Query= ") == 10:
            if stdoutdata.count("***** No hits found *****") == 7:
                # This happens with BLAST 2.2.26+ which is potentially a bug
                pass
            else:
                self.assertEqual(9, stdoutdata.count("***** No hits found *****"))
        else:
            # Assume this is NCBI BLAST+ 2.3.0 or later,
            self.assertEqual(1, stdoutdata.count("Query= "))
            self.assertEqual(0, stdoutdata.count("***** No hits found *****"))

    def test_blastn(self):
        """Pairwise BLASTN search"""
        global exe_names
        cline = Applications.NcbiblastnCommandline(
            exe_names["blastn"],
            query="GenBank/NC_005816.ffn",
            subject="GenBank/NC_005816.fna",
            evalue="0.000001")
        self.assertEqual(str(cline), _escape_filename(exe_names["blastn"]) +
                         " -query GenBank/NC_005816.ffn -evalue 0.000001" +
                         " -subject GenBank/NC_005816.fna")
        child = subprocess.Popen(str(cline),
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True,
                                 shell=(sys.platform != "win32"))
        stdoutdata, stderrdata = child.communicate()
        return_code = child.returncode
        self.assertEqual(return_code, 0, "Got error code %i back from:\n%s"
                         % (return_code, cline))
        self.assertEqual(10, stdoutdata.count("Query= "))
        self.assertEqual(0, stdoutdata.count("***** No hits found *****"))
        # TODO - Parse it?

    def test_tblastn(self):
        """Pairwise TBLASTN search"""
        global exe_names
        cline = Applications.NcbitblastnCommandline(
            exe_names["tblastn"],
            query="GenBank/NC_005816.faa",
            subject="GenBank/NC_005816.fna",
            evalue="1e-6")
        self.assertEqual(str(cline), _escape_filename(exe_names["tblastn"]) +
                         " -query GenBank/NC_005816.faa -evalue 1e-6" +
                         " -subject GenBank/NC_005816.fna")
        child = subprocess.Popen(str(cline),
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True,
                                 shell=(sys.platform != "win32"))
        stdoutdata, stderrdata = child.communicate()
        return_code = child.returncode
        self.assertEqual(return_code, 0, "Got error code %i back from:\n%s"
                         % (return_code, cline))
        self.assertEqual(10, stdoutdata.count("Query= "))
        self.assertEqual(0, stdoutdata.count("***** No hits found *****"))
        # TODO - Parse it?


class BlastDB(unittest.TestCase):
    def test_requires_dbtype(self):
        """Check that dbtype throws error if not set"""
        global exe_names
        cline = Applications.NcbimakeblastdbCommandline(
            exe_names["makeblastdb"],
            input_file="GenBank/NC_005816.faa")
        with self.assertRaises(ValueError):
            str(cline)

    def test_fasta_db_prot(self):
        """Test makeblastdb wrapper with protein database"""
        global exe_names
        cline = Applications.NcbimakeblastdbCommandline(
            exe_names["makeblastdb"],
            input_file="GenBank/NC_005816.faa",
            dbtype="prot",
            hash_index=True,
            max_file_sz="20MB",
            parse_seqids=True,
            taxid=10)

        self.assertEqual(str(cline),
                         _escape_filename(exe_names["makeblastdb"]) +
                         " -dbtype prot -in GenBank/NC_005816.faa"
                         " -parse_seqids -hash_index -max_file_sz 20MB"
                         " -taxid 10")

        child = subprocess.Popen(str(cline),
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True,
                                 shell=(sys.platform != "win32"))
        stdoutdata, stderrdata = child.communicate()
        return_code = child.returncode

        self.assertTrue(os.path.isfile("GenBank/NC_005816.faa.phd"))
        self.assertTrue(os.path.isfile("GenBank/NC_005816.faa.phi"))
        self.assertTrue(os.path.isfile("GenBank/NC_005816.faa.phr"))
        self.assertTrue(os.path.isfile("GenBank/NC_005816.faa.pin"))
        self.assertTrue(os.path.isfile("GenBank/NC_005816.faa.pog"))
        self.assertTrue(os.path.isfile("GenBank/NC_005816.faa.psd"))
        self.assertTrue(os.path.isfile("GenBank/NC_005816.faa.psi"))
        self.assertTrue(os.path.isfile("GenBank/NC_005816.faa.psq"))

    def test_fasta_db_nucl(self):
        """Test makeblastdb wrapper with nucleotide database"""
        global exe_names
        cline = Applications.NcbimakeblastdbCommandline(
            exe_names["makeblastdb"],
            input_file="GenBank/NC_005816.fna",
            dbtype="nucl",
            hash_index=True,
            max_file_sz="20MB",
            parse_seqids=True,
            taxid=10)

        self.assertEqual(str(cline),
                         _escape_filename(exe_names["makeblastdb"]) +
                         " -dbtype nucl -in GenBank/NC_005816.fna"
                         " -parse_seqids -hash_index -max_file_sz 20MB"
                         " -taxid 10")

        child = subprocess.Popen(str(cline),
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True,
                                 shell=(sys.platform != "win32"))
        stdoutdata, stderrdata = child.communicate()
        return_code = child.returncode

        self.assertTrue(os.path.isfile("GenBank/NC_005816.fna.nhd"))
        self.assertTrue(os.path.isfile("GenBank/NC_005816.fna.nhi"))
        self.assertTrue(os.path.isfile("GenBank/NC_005816.fna.nhr"))
        self.assertTrue(os.path.isfile("GenBank/NC_005816.fna.nin"))
        self.assertTrue(os.path.isfile("GenBank/NC_005816.fna.nog"))
        self.assertTrue(os.path.isfile("GenBank/NC_005816.fna.nsd"))
        self.assertTrue(os.path.isfile("GenBank/NC_005816.fna.nsi"))
        self.assertTrue(os.path.isfile("GenBank/NC_005816.fna.nsq"))

    # makeblastdb makes files in the same dir as the input, clean these up
    def tearDown(self):
        blastdb_matcher_prot = re.compile(r'NC_005816\.faa\.p.+')
        for file in os.listdir('GenBank/'):
            if blastdb_matcher_prot.match(file):
                path = os.path.join('GenBank/', file)
                os.remove(path)

        blastdb_matcher_nucl = re.compile(r'NC_005816\.fna\.n.+')
        for file in os.listdir('GenBank/'):
            if blastdb_matcher_nucl.match(file):
                path = os.path.join('GenBank/', file)
                os.remove(path)


class CheckCompleteArgList(unittest.TestCase):
    def check(self, exe_name, wrapper):
        global exe_names
        exe = exe_names[exe_name]
        # dbtype must be set to initialize NcbimakeblastdbCommandline
        if exe_name == 'makeblastdb':
            cline = wrapper(exe, h=True, dbtype='prot')
        else:
            cline = wrapper(exe, h=True)
        names = set(parameter.names[0]
                    for parameter in cline.parameters)

        child = subprocess.Popen(str(cline),
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True,
                                 shell=(sys.platform != "win32"))
        stdoutdata, stderrdata = child.communicate()
        self.assertEqual(stderrdata, "",
                         "%s\n%s" % (str(cline), stderrdata))
        names_in_tool = set()
        while stdoutdata:
            index = stdoutdata.find("[")
            if index == -1:
                break
            stdoutdata = stdoutdata[index + 1:]
            index = stdoutdata.find("]")
            assert index != -1
            name = stdoutdata[:index]
            if " " in name:
                name = name.split(None, 1)[0]
            names_in_tool.add(name)
            stdoutdata = stdoutdata[index + 1:]

        extra = names.difference(names_in_tool)
        missing = names_in_tool.difference(names)
        if "-verbose" in missing:
            # Known issue, seems to be present in some builds (Bug 3043)
            missing.remove("-verbose")
        if "-remote_verbose" in missing:
            # Known issue, seems to be present in some builds (Bug 3043)
            missing.remove("-remote_verbose")
        if "-use_test_remote_service" in missing:
            # Known issue, seems to be present in some builds (Bug 3043)
            missing.remove("-use_test_remote_service")
        if exe_name == "blastn" and "-off_diagonal_range" in extra:
            # Added in BLAST 2.2.23+
            extra.remove("-off_diagonal_range")
        if exe_name == "tblastx":
            # These appear to have been removed in BLAST 2.2.23+
            # (which seems a bit odd - TODO - check with NCBI?)
            extra = extra.difference(["-gapextend", "-gapopen",
                                      "-xdrop_gap", "-xdrop_gap_final"])
        if exe_name in ["rpsblast", "rpstblastn"]:
            # These appear to have been removed in BLAST 2.2.24+
            # (which seems a bit odd - TODO - check with NCBI?)
            extra = extra.difference(["-num_threads"])
        if exe_name in ["tblastn", "tblastx"]:
            # These appear to have been removed in BLAST 2.2.24+
            extra = extra.difference(["-db_soft_mask"])
        # This was added in BLAST 2.2.24+ to most/all the tools, so
        # will be seen as an extra argument on older versions:
        if "-seqidlist" in extra:
            extra.remove("-seqidlist")
        if "-db_hard_mask" in extra \
        and exe_name in ["blastn", "blastp", "blastx", "tblastx", "tblastn"]:
            # New in BLAST 2.2.25+ so will look like an extra arg on old BLAST
            extra.remove("-db_hard_mask")
        if "-msa_master_idx" in extra and exe_name == "psiblast":
            # New in BLAST 2.2.25+ so will look like an extra arg on old BLAST
            extra.remove("-msa_master_idx")
        if exe_name == "rpsblast":
            # New in BLAST 2.2.25+ so will look like an extra arg on old BLAST
            extra = extra.difference(["-best_hit_overhang",
                                      "-best_hit_score_edge",
                                      "-culling_limit"])
        if "-max_hsps_per_subject" in extra:
            # New in BLAST 2.2.26+ so will look like an extra arg on old BLAST
            extra.remove("-max_hsps_per_subject")
        if "-ignore_msa_master" in extra and exe_name == "psiblast":
            # New in BLAST 2.2.26+ so will look like an extra arg on old BLAST
            extra.remove("-ignore_msa_master")
        if exe_name == "blastx":
            # New in BLAST 2.2.27+ so will look like an extra arg on old BLAST
            extra = extra.difference(["-comp_based_stats",
                                      "-use_sw_tback"])
        if exe_name in ["blastx", "tblastn"]:
            # Removed in BLAST 2.2.27+ so will look like extra arg on new BLAST
            extra = extra.difference(["-frame_shift_penalty"])
        if exe_name in ["blastn", "blastp", "blastx", "tblastn", "tblastx",
                        "psiblast", "rpstblastn", "rpsblast"]:
            # New in BLAST 2.2.29+ so will look like extra args on old BLAST:
            extra = extra.difference(["-max_hsps", "-sum_statistics"])
        if exe_name in ["rpstblastn", "rpsblast"]:
            # Removed in BLAST 2.2.29+ so will look like extra args on new BLAST
            extra = extra.difference(["-gilist", "-negative_gilist"])
            # Removed in BLAST 2.2.30 so will look like extra args on new BLAST
            # Apparently -word_size should never have been added to these tools.
            extra = extra.difference(["-word_size"])
            # New in BLAST 2.2.28+ (for rpsblast) and BLAST 2.6+ (for rpstblastn)
            # so will look like extra args on old BLAST:
            extra = extra.difference(["-comp_based_stats", "-use_sw_tback"])
        if exe_name == "deltablast":
            # New in BLAST+ 2.2.29 so will look like extra args on BLAST+ 2.2.28
            extra = extra.difference(["-entrez_query", "-max_hsps", "-sum_statistics"])
            # New in BLAST+ 2.2.27 so looks like extra arg on BLAST+ 2.2.26
            extra = extra.difference(["-remote"])
        if exe_name in ["blastx", "tblastn"]:
            # New in BLAST+ 2.2.30 so will look like extra args on BLAST+ 2.2.29 etc
            extra = extra.difference(["-task"])
        if exe_name in ["blastn", "blastp", "blastx", "deltablast", "psiblast",
                        "rpstblastn", "rpsblast", "tblastn", "tblastx"]:
            # New in BLAST+ 2.2.30 so will look like extra args on BLAST+ 2.2.29 etc
            extra = extra.difference(["-line_length", "-qcov_hsp_perc", "-sum_stats"])
        if exe_name in ["deltablast", "psiblast"]:
            # New in BLAST+ 2.3.0 so will look like extra args on older verions
            extra = extra.difference(["-save_each_pssm", "-save_pssm_after_last_round"])
        if exe_name == "makeblastdb":
            # -input_type is new in BLAST+ 2.2.26+ so will look like extra
            # arg on older versions. -mask_desc and -mask_id were added after
            # 2.2.28+, but not documented in changelog.
            # -dbtype also ignored here, as it was set to initialize class.
            extra = extra.difference(["-input_type", "-mask_desc",
                                      "-mask_id", "-dbtype"])
        # This was added in BLAST+ 2.7.1 (or maybe 2.7.0) to most/all the tools,
        # so will be seen as an extra argument on older versions:
        if "-negative_seqidlist" in extra:
            extra.remove("-negative_seqidlist")

        if extra or missing:
            import warnings
            warnings.warn("NCBI BLAST+ %s and Biopython out sync. Please "
                          "update Biopython, or report this issue if you are "
                          "already using the latest version. (Extra args: %s; "
                          "Missing: %s)" % (exe_name,
                                            ",".join(sorted(extra)),
                                            ",".join(sorted(missing))))

        # An almost trivial example to test any validation
        if "-query" in names:
            cline = wrapper(exe, query="dummy")
        elif "-archive" in names:
            cline = wrapper(exe, archive="dummy")
        str(cline)

    def test_blastx(self):
        """Check all blastx arguments are supported"""
        self.check("blastx", Applications.NcbiblastxCommandline)

    def test_blastp(self):
        """Check all blastp arguments are supported"""
        self.check("blastp", Applications.NcbiblastpCommandline)

    def test_blastn(self):
        """Check all blastn arguments are supported"""
        self.check("blastn", Applications.NcbiblastnCommandline)

    def test_tblastx(self):
        """Check all tblastx arguments are supported"""
        self.check("tblastx", Applications.NcbitblastxCommandline)

    def test_tblastn(self):
        """Check all tblastn arguments are supported"""
        self.check("tblastn", Applications.NcbitblastnCommandline)

    def test_psiblast(self):
        """Check all psiblast arguments are supported"""
        self.check("psiblast", Applications.NcbipsiblastCommandline)

    def test_rpsblast(self):
        """Check all rpsblast arguments are supported"""
        self.check("rpsblast", Applications.NcbirpsblastCommandline)

    def test_rpstblastn(self):
        """Check all rpstblastn arguments are supported"""
        self.check("rpstblastn", Applications.NcbirpstblastnCommandline)

    def test_makeblastdb(self):
        """Check all makeblastdb arguments are supported"""
        self.check("makeblastdb", Applications.NcbimakeblastdbCommandline)

    if "blast_formatter" in exe_names:
        def test_blast_formatter(self):
            """Check all blast_formatter arguments are supported"""
            self.check("blast_formatter", Applications.NcbiblastformatterCommandline)

    if "deltablast" in exe_names:
        def test_deltablast(self):
            """Check all deltablast arguments are supported"""
            self.check("deltablast", Applications.NcbideltablastCommandline)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
