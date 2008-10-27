# Copyright 2008 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio import MissingExternalDependencyError
import sys
if sys.platform=="win32" :
    #Need to check the path, and then common places like under
    #C:\Program Files\Clustalw\Clustalw.exe etc.
    raise MissingExternalDependencyError(\
        "Don't know how to find the tool clustalw on Windows.")
    clustalw_exe = None
else :
    import commands
    #Note that clustalw 1.83 doesn't obey the --version command,
    #but this does cause it to quit cleanly.  Otherwise it prompts
    #the user for input!
    if "not found" in commands.getoutput("clustalw --version") :
        raise MissingExternalDependencyError(\
            "Install clustalw if you want to use Bio.Clustalw.")
    clustalw_exe = "clustalw"

#################################################################
import os
from Bio import Clustalw
from Bio.Clustalw import MultipleAlignCL
from Bio import SeqIO

#Create a temp fasta file with a space in the name
temp_filename_with_spaces = "Clustalw/temp horses.fasta"
handle = open(temp_filename_with_spaces, "w")
SeqIO.write(SeqIO.parse(open("Phylip/hennigian.phy"),"phylip"),handle, "fasta")
handle.close()

for input_file, output_file, newtree_file in [
    ("Fasta/f002", "temp_test.aln", None),
    ("GFF/multi.fna", "temp with space.aln", None),
    ("Registry/seqs.fasta", "temp_test.aln", None),
    ("Registry/seqs.fasta", "temp_test.aln", "temp_test.dnd"),
    ("Registry/seqs.fasta", "temp_test.aln", "temp with space.dnd"),
    (temp_filename_with_spaces, "temp_test.aln", None),
    (temp_filename_with_spaces, "temp with space.aln", None),
    ] :
    input_records = SeqIO.to_dict(SeqIO.parse(open(input_file),"fasta"))
    if os.path.isfile(output_file) :
        os.remove(output_file)
    print "Calling clustalw on %s (with %i records)" \
          % (repr(input_file), len(input_records))
    print "using output file %s" % repr(output_file)
    if newtree_file is not None :
        print "requesting output guide tree file %s" % repr(newtree_file)

    #Prepare the command...
    cline = MultipleAlignCL(input_file, command='clustalw')
    cline.set_output(output_file)
    if newtree_file is not None :
        cline.set_new_guide_tree(newtree_file)

    #Run the command...
    align = Clustalw.do_alignment(cline)

    #Check the output...
    print "Got an alignment, %i sequences of length %i" \
          % (len(align.get_all_seqs()), align.get_alignment_length())
    output_records = SeqIO.to_dict(SeqIO.parse(open(output_file),"clustal"))
    assert set(input_records.keys()) == set(output_records.keys())
    for record in align :
        assert str(record.seq) == str(output_records[record.id].seq)
        assert str(record.seq).replace("-","") == \
               str(input_records[record.id].seq)
    os.remove(output_file)
    if newtree_file is not None :
        assert os.path.isfile(newtree_file)
        #TODO - Try and parse this with Bio.Nexus?
        os.remove(newtree_file)
    print

#Clean up temp file..
os.remove(temp_filename_with_spaces)
print "Done"

