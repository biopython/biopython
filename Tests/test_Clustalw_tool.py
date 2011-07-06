# Copyright 2008-2011 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

#TODO - Clean up the extra files created by clustalw?  e.g. *.dnd
#and *.aln where we have not requested an explicit name?
from Bio import MissingExternalDependencyError

import sys
import os
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.Application import ApplicationError

#################################################################

#Try to avoid problems when the OS is in another language
os.environ['LANG'] = 'C'

clustalw_exe = None
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
                    clustalw_exe = os.path.join(prog_files, folder, filename)
                    break
            if clustalw_exe : break
else:
    import commands
    #Note that clustalw 1.83 and clustalw 2.0.10 don't obey the --version
    #command, but this does cause them to quit cleanly.  Otherwise they prompt
    #the user for input (causing a lock up).
    output = commands.getoutput("clustalw2 --version")
    #Since "not found" may be in another language, try and be sure this is
    #really the clustalw tool's output
    if "not found" not in output and "CLUSTAL" in output \
    and "Multiple Sequence Alignments" in output:
        clustalw_exe = "clustalw2"
    if not clustalw_exe:
        output = commands.getoutput("clustalw --version")
        if "not found" not in output and "CLUSTAL" in output \
        and "Multiple Sequence Alignments" in output:
            clustalw_exe = "clustalw"

if not clustalw_exe:
    raise MissingExternalDependencyError(\
        "Install clustalw or clustalw2 if you want to use it from Biopython.")

#################################################################

print "Checking error conditions"
print "========================="

print "Empty file"
input_file = "does_not_exist.fasta"
assert not os.path.isfile(input_file)
cline = ClustalwCommandline(clustalw_exe, infile=input_file)
try:
    stdout, stderr = cline()
    assert False, "Should have failed, returned:\n%s\n%s" % (stdout, stderr)
except ApplicationError, err:
    print "Failed (good)"
    #Python 2.3 on Windows gave (0, 'Error')
    #Python 2.5 on Windows gives [Errno 0] Error
    assert "Cannot open sequence file" in str(err) or \
           "Cannot open input file" in str(err) or \
           "non-zero exit status" in str(err), str(err)

print
print "Single sequence"
input_file = "Fasta/f001"
assert os.path.isfile(input_file)
assert len(list(SeqIO.parse(input_file,"fasta")))==1
cline = ClustalwCommandline(clustalw_exe, infile=input_file)
try:
    stdout, stderr = cline()
    if "cannot do multiple alignment" in (stdout + stderr):
        #Zero return code is a possible bug in clustal?
        print "Failed (good)"
    else:
        assert False, "Should have failed, returned:\n%s\n%s" % (stdout, stderr)
except ApplicationError, err:
    print "Failed (good)"
    assert str(err) == "No records found in handle", str(err)

print
print "Invalid sequence"
input_file = "Medline/pubmed_result1.txt"
assert os.path.isfile(input_file)
cline = ClustalwCommandline(clustalw_exe, infile=input_file)
try:
    stdout, stderr = cline()
    assert False, "Should have failed, returned:\n%s\n%s" % (stdout, stderr)
except ApplicationError, err:
    print "Failed (good)"
    #Ideally we'd catch the return code and raise the specific
    #error for "invalid format", rather than just notice there
    #is not output file.
    #Note:
    #Python 2.3 on Windows gave (0, 'Error')
    #Python 2.5 on Windows gives [Errno 0] Error
    assert "invalid format" in str(err) \
           or "not produced" in str(err) \
           or "No sequences in file" in str(err) \
           or "non-zero exit status " in str(err), str(err)

#################################################################
print
print "Checking normal situations"
print "=========================="

#Create a temp fasta file with a space in the name
temp_filename_with_spaces = "Clustalw/temp horses.fasta"
handle = open(temp_filename_with_spaces, "w")
SeqIO.write(SeqIO.parse("Phylip/hennigian.phy","phylip"),handle, "fasta")
handle.close()

#Create a large input file by converting another example file
#(See Bug 2804, this will produce so much output on stdout that
#subprocess could suffer a deadlock and hang).  Using all the
#records should show the deadlock but is very slow - just thirty
#seems to lockup on Mac OS X, even 20 on Linux (without the fix).
temp_large_fasta_file = "temp_cw_prot.fasta"
handle = open(temp_large_fasta_file, "w")
records = list(SeqIO.parse("NBRF/Cw_prot.pir", "pir"))[:40]
SeqIO.write(records, handle, "fasta")
handle.close()
del handle, records

for input_file, output_file, newtree_file in [
    ("Fasta/f002", "temp_test.aln", None),
    ("GFF/multi.fna", "temp with space.aln", None),
    ("Registry/seqs.fasta", "temp_test.aln", None),
    ("Registry/seqs.fasta", "temp_test.aln", "temp_test.dnd"),
    ("Registry/seqs.fasta", "temp_test.aln", "temp with space.dnd"),
    (temp_filename_with_spaces, "temp_test.aln", None),
    (temp_filename_with_spaces, "temp with space.aln", None),
    (temp_large_fasta_file, "temp_cw_prot.aln", None),
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

    #Any filesnames with spaces should get escaped with quotes automatically.
    #Using keyword arguments here.
    cline = ClustalwCommandline(clustalw_exe,
                                infile=input_file,
                                outfile=output_file)
    assert str(eval(repr(cline)))==str(cline)
    if newtree_file is not None:
        #Test using a property:
        cline.newtree = newtree_file
        #I don't just want the tree, also want the alignment:
        cline.align = True
        assert str(eval(repr(cline)))==str(cline)
    #print cline
    output, error = cline()
    assert output.strip().startswith("CLUSTAL")
    assert error.strip() == ""
    #Check the output...
    align = AlignIO.read(output_file, "clustal")
    #The length of the alignment will depend on the version of clustalw
    #(clustalw 2.0.10 and clustalw 1.83 are certainly different).
    print "Got an alignment, %i sequences" % (len(align))
    output_records = SeqIO.to_dict(SeqIO.parse(output_file,"clustal"))
    assert set(input_records.keys()) == set(output_records.keys())
    for record in align:
        assert str(record.seq) == str(output_records[record.id].seq)
        assert str(record.seq).replace("-","") == \
               str(input_records[record.id].seq)

    #Clean up...
    os.remove(output_file)

    #Check the DND file was created.
    #TODO - Try and parse this with Bio.Nexus?
    if newtree_file is not None:
        tree_file = newtree_file
    else:
        #Clustalw will name it based on the input file
        tree_file = os.path.splitext(input_file)[0] + ".dnd"
    assert os.path.isfile(tree_file), \
           "Did not find tree file %s" % tree_file
    os.remove(tree_file)

#Clean up any stray temp files..
if os.path.isfile("Fasta/f001.aln"):
    os.remove("Fasta/f001.aln")
if os.path.isfile("Medline/pubmed_result1.aln"):
    os.remove("Medline/pubmed_result1.aln")
if os.path.isfile(temp_filename_with_spaces):
    os.remove(temp_filename_with_spaces)
if os.path.isfile(temp_large_fasta_file):
    os.remove(temp_large_fasta_file)

print "Done"
