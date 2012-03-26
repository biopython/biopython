# Copyright 2008-2011 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio import MissingExternalDependencyError

import sys
import os
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Application import ApplicationError

#################################################################

#Try to avoid problems when the OS is in another language
os.environ['LANG'] = 'C'

clustalo_exe = None
if sys.platform=="win32":
    #TODO
    raise MissingExternalDependencyError("Testing this on Windows not implemented yet")
else:
    import commands
    output = commands.getoutput("clustalo --help")
    if output.startswith("Clustal Omega"):
        clustalo_exe = "clustalo"

if not clustalo_exe:
    raise MissingExternalDependencyError(\
        "Install clustalo if you want to use Clustal Omega from Biopython.")

#################################################################

print "Checking error conditions"
print "========================="

print "Empty file"
input_file = "does_not_exist.fasta"
assert not os.path.isfile(input_file)
cline = ClustalOmegaCommandline(clustalo_exe, infile=input_file)
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
cline = ClustalOmegaCommandline(clustalo_exe, infile=input_file)
try:
    stdout, stderr = cline()
    assert False, "Should have failed, returned:\n%s\n%s" % (stdout, stderr)
except ApplicationError, err:
    print "Failed (good)"
    assert "contains 1 sequence, nothing to align" in str(err), str(err)

print
print "Invalid sequence"
input_file = "Medline/pubmed_result1.txt"
assert os.path.isfile(input_file)
cline = ClustalOmegaCommandline(clustalo_exe, infile=input_file)
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
    assert "Can't determine format of sequence file" in str(err), str(err)

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
    ("Registry/seqs.fasta", "temp with space.aln", None),
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
    #Using keyword arguments here, force=True to over-write the output file.
    cline = ClustalOmegaCommandline(clustalo_exe,
                                    infile=input_file,
                                    outfile=output_file,
                                    outfmt="clustal",
                                    force=True)
    assert str(eval(repr(cline)))==str(cline)
    if newtree_file is not None:
        #Test using a property:
        cline.guidetree_out = newtree_file
        assert str(eval(repr(cline)))==str(cline)
    #print cline
    output, error = cline()
    #assert output, "No output from: %s\n%s" % (cline, error)
    assert not output or output.strip().startswith("CLUSTAL"), output
    assert error.strip() == "" or \
           error.startswith("WARNING: Sequence type is DNA."), error
    #Check the output...
    align = AlignIO.read(output_file, "clustal")
    #The length of the alignment will depend on the version of clustalw
    #(clustalw 2.0.10 and clustalw 1.83 are certainly different).
    print "Got an alignment, %i sequences" % (len(align))
    output_records = SeqIO.to_dict(SeqIO.parse(output_file,"clustal"))
    assert len(set(input_records.keys())) == len(set(output_records.keys())), \
        "%r vs %r" %(sorted(input_records.keys()), sorted(output_records.keys()))
    #TODO - Investigate Clustal Omega underscore/semi-colon mapping
    #for record in align:
    #    assert str(record.seq) == str(output_records[record.id].seq)
    #    assert str(record.seq).replace("-","") == \
    #           str(input_records[record.id].seq)

    #Clean up...
    os.remove(output_file)

    #TODO Check the DND file was created.
    #TODO - Try and parse this with Bio.Nexus?
    if newtree_file is not None \
    and os.path.isfile(newtree_file):
        os.remove(newtree_file)

#Clean up any stray temp files..
if os.path.isfile(temp_filename_with_spaces):
    os.remove(temp_filename_with_spaces)
if os.path.isfile(temp_large_fasta_file):
    os.remove(temp_large_fasta_file)

print "Done"
