# Copyright 2009-2010 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
# This unit test attempts to locate the blastall executable and the nr
# database, and if it finds them then do some standalone blast searches
# using Bio.Blast.NCBIStandalone to call the command line tool.

import os, sys
import subprocess
import unittest

from Bio import MissingExternalDependencyError
from Bio.Blast import Applications

# TODO - On windows, can we use the ncbi.ini file?
wanted = ["blastx", "blastp", "blastn", "tblastn", "tblastx",
          "rpsblast", "rpstblastn", "psiblast", "blast_formatter"]
exe_names = {}

if sys.platform=="win32":
    #The Windows 32 bit BLAST 2.2.22+ installer does add itself to the path,
    #and by default installs to C:\Program Files\NCBI\BLAST-2.2.22+\bin
    #To keep things simple, assume BLAST+ is on the path on Windows.
    #
    #On Windows the environment variable name isn't case senstive,
    #but must split on ";" not ":"
    likely_dirs = os.environ.get("PATH", "").split(";")
else :
    likely_dirs = os.environ.get("PATH", "").split(":")

for folder in likely_dirs:
    if not os.path.isdir(folder): continue
    for name in wanted :
        if sys.platform=="win32":
            exe_name = os.path.join(folder, name+".exe")
        else:
            exe_name = os.path.join(folder, name)
        if not os.path.isfile(exe_name):
            continue
        #To tell the old and new rpsblast apart (since I have both on
        #my path and the old blast has priority), try -h as a parameter.
        #This should also reject WU-BLAST (since it doesn't like -h).
        child = subprocess.Popen(exe_name + " -h",
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True,
                                 shell=(sys.platform!="win32"))
        output, error = child.communicate()
        if child.returncode==0 and "ERROR: Invalid argument: -h" not in output:
            #Special case, blast_formatter from BLAST 2.2.23+ (i.e. BLAST+)
            #has mandatory argument -rid, but no -archive. We don't support it.
            if name == "blast_formatter" and " -archive " not in output:
                continue
            exe_names[name] = exe_name
        #else :
        #    print "Rejecting", exe_name
        del exe_name, name

#We can cope with blast_formatter being missing, only added in BLAST 2.2.24+
if len(set(exe_names).difference(["blast_formatter"])) < len(wanted)-1 :
    raise MissingExternalDependencyError("Install the NCBI BLAST+ command line "
                                         "tools if you want to use the "
                                         "Bio.Blast.Applications wrapper.")


class Pairwise(unittest.TestCase):
    def test_blastp(self):
        """Pairwise BLASTP search"""
        global exe_names
        cline = Applications.NcbiblastpCommandline(exe_names["blastp"],
                        query="Fasta/rose.pro",
                        subject="GenBank/NC_005816.faa",
                        evalue=1)
        self.assertEqual(str(cline), exe_names["blastp"] \
                         + " -query Fasta/rose.pro -evalue 1" \
                         + " -subject GenBank/NC_005816.faa")
        child = subprocess.Popen(str(cline),
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True,
                                 shell=(sys.platform!="win32"))
        stdoutdata, stderrdata = child.communicate()
        return_code = child.returncode
        self.assertEqual(return_code, 0, "Got error code %i back from:\n%s"
                         % (return_code, cline))
        self.assertEqual(10, stdoutdata.count("Query= "))
        self.assertEqual(9, stdoutdata.count("***** No hits found *****"))
        
        #TODO - Parse it? I think we'd need to update this obsole code :(
        #records = list(NCBIStandalone.Iterator(StringIO(stdoutdata),
        #                                       NCBIStandalone.BlastParser()))   

    def test_blastn(self):
        """Pairwise BLASTN search"""
        global exe_names
        cline = Applications.NcbiblastnCommandline(exe_names["blastn"],
                        query="GenBank/NC_005816.ffn",
                        subject="GenBank/NC_005816.fna",
                        evalue="0.000001")
        self.assertEqual(str(cline), exe_names["blastn"] \
                         + " -query GenBank/NC_005816.ffn -evalue 0.000001" \
                         + " -subject GenBank/NC_005816.fna")
        child = subprocess.Popen(str(cline),
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True,
                                 shell=(sys.platform!="win32"))
        stdoutdata, stderrdata = child.communicate()
        return_code = child.returncode
        self.assertEqual(return_code, 0, "Got error code %i back from:\n%s"
                         % (return_code, cline))
        self.assertEqual(10, stdoutdata.count("Query= "))
        self.assertEqual(0, stdoutdata.count("***** No hits found *****"))
        #TODO - Parse it?

    def test_tblastn(self):
        """Pairwise TBLASTN search"""
        global exe_names
        cline = Applications.NcbitblastnCommandline(exe_names["tblastn"],
                        query="GenBank/NC_005816.faa",
                        subject="GenBank/NC_005816.fna",
                        evalue="1e-6")
        self.assertEqual(str(cline), exe_names["tblastn"] \
                         + " -query GenBank/NC_005816.faa -evalue 1e-6" \
                         + " -subject GenBank/NC_005816.fna")
        child = subprocess.Popen(str(cline),
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True,
                                 shell=(sys.platform!="win32"))
        stdoutdata, stderrdata = child.communicate()
        return_code = child.returncode
        self.assertEqual(return_code, 0, "Got error code %i back from:\n%s"
                         % (return_code, cline))
        self.assertEqual(10, stdoutdata.count("Query= "))
        self.assertEqual(0, stdoutdata.count("***** No hits found *****"))
        #TODO - Parse it?

   
class CheckCompleteArgList(unittest.TestCase):
    def check(self, exe_name, wrapper) :
        global exe_names
        exe = exe_names[exe_name]
        cline = wrapper(exe, h=True)

        names = set(parameter.names[0] \
                    for parameter in cline.parameters)
        
        child = subprocess.Popen(str(cline),
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True,
                                 shell=(sys.platform!="win32"))
        stdoutdata, stderrdata = child.communicate()
        self.assertEqual(stderrdata, "",
                         "%s\n%s" % (str(cline), stderrdata))
        names_in_tool = set()
        while stdoutdata :
            index = stdoutdata.find("[")
            if index == -1 : break
            stdoutdata = stdoutdata[index+1:]
            index = stdoutdata.find("]")
            assert index != -1
            name = stdoutdata[:index]
            if " " in name : name = name.split(None,1)[0]
            names_in_tool.add(name)
            stdoutdata = stdoutdata[index+1:]
                
        extra = names.difference(names_in_tool)
        missing = names_in_tool.difference(names)
        if "-soft_masking" in missing :
            #Known issue, need to establish how this option works
            missing.remove("-soft_masking")
        if "-use_index" in missing :
            #Known issue, need to establish how this option works
            missing.remove("-use_index")
        if "-verbose" in missing :
            #Known issue, seems to be present in some builds (Bug 3043)
            missing.remove("-verbose")
        if "-remote_verbose" in missing :
            #Known issue, seems to be present in some builds (Bug 3043)
            missing.remove("-remote_verbose")
        if "-use_test_remote_service" in missing :
            #Known issue, seems to be present in some builds (Bug 3043)
            missing.remove("-use_test_remote_service")
        if exe_name == "blastn" and "-off_diagonal_range" in extra:
            #Added in BLAST 2.2.23+
            extra.remove("-off_diagonal_range")
        if exe_name == "tblastx":
            #These appear to have been removed in BLAST 2.2.23+
            #(which seems a bit odd - TODO - check with NCBI?)
            extra = extra.difference(["-gapextend","-gapopen",
                                      "-xdrop_gap","-xdrop_gap_final"])
        if exe_name in ["rpsblast", "rpstblastn"]:
            #These appear to have been removed in BLAST 2.2.24+
            #(which seems a bit odd - TODO - check with NCBI?)
            extra = extra.difference(["-num_threads"])
        if exe_name in ["tblastn", "tblastx"]:
            #These appear to have been removed in BLAST 2.2.24+
            extra = extra.difference(["-db_soft_mask"])
        #This was added in BLAST 2.2.24+ to most/all the tools, so
        #will be seen as an extra argument on older versions:
        if "-seqidlist" in extra:
            extra.remove("-seqidlist")
        if "-db_hard_mask" in extra \
        and exe_name in ["blastn", "blastp", "blastx", "tblastx", "tblastn"]:
            #New in BLAST 2.2.25+ so will look like an extra arg on old BLAST
            extra.remove("-db_hard_mask")
        if "-msa_master_idx" in extra and exe_name=="psiblast":
            #New in BLAST 2.2.25+ so will look like an extra arg on old BLAST
            extra.remove("-msa_master_idx")
        if exe_name=="rpsblast":
            #New in BLAST 2.2.25+ so will look like an extra arg on old BLAST
            extra = extra.difference(["-best_hit_overhang",
                                      "-best_hit_score_edge",
                                      "-culling_limit"])

        if extra or missing:
            import warnings
            warnings.warn("NCBI BLAST+ %s and Biopython out sync. Please "
                          "update Biopython, or report this issue if you are "
                          "already using the latest version. (Exta args: %s; "
                          "Missing: %s)" % (exe_name,
                          ",".join(sorted(extra)),
                          ",".join(sorted(missing))))

        #An almost trivial example to test any validation
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

    if "blast_formatter" in exe_names:
        def test_blast_formatter(self):
            """Check all blast_formatter arguments are supported"""
            self.check("blast_formatter", Applications.NcbiblastformatterCommandline)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
