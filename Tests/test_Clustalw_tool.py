# Copyright 2008 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

#TODO - Clean up the extra files created by clustalw?  e.g. *.dnd
#and *.aln where we have not requested an explicit name?
from Bio import MissingExternalDependencyError

#TODO - Remove this work around once we drop python 2.3 support
try:
    set = set
except NameError:
    from sets import Set as set

import sys
import os
from Bio import Clustalw
from Bio.Clustalw import MultipleAlignCL
from Bio import SeqIO

#################################################################

clustalw_exe = None
if sys.platform=="win32" :
    #TODO - Check the path?
    try :
        #This can vary depending on the Windows language.
        prog_files = os.environ["PROGRAMFILES"]
    except KeyError :
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
    for folder in likely_dirs :
        if os.path.isdir(os.path.join(prog_files, folder)) :
            for filename in likely_exes :
                if os.path.isfile(os.path.join(prog_files, folder, filename)) :
                    clustalw_exe = os.path.join(prog_files, folder, filename)
                    break
            if clustalw_exe : break
else :
    import commands
    #Note that clustalw 1.83 and clustalw 2.0.10 don't obey the --version
    #command, but this does cause them to quit cleanly.  Otherwise they prompt
    #the user for input (causing a lock up).
    output = commands.getoutput("clustalw2 --version")
    if "not found" not in output and "CLUSTAL" in output.upper() :
        clustalw_exe = "clustalw2"
    if not clustalw_exe :
        output = commands.getoutput("clustalw --version")
        if "not found" not in output and "CLUSTAL" in output.upper() :
            clustalw_exe = "clustalw"

if not clustalw_exe :
    raise MissingExternalDependencyError(\
        "Install clustalw or clustalw2 if you want to use Bio.Clustalw.")

#################################################################

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
    cline = MultipleAlignCL(input_file, command=clustalw_exe)
    cline.set_output(output_file)
    if newtree_file is not None :
        cline.set_new_guide_tree(newtree_file)

    #Run the command...
    align = Clustalw.do_alignment(cline)

    #Check the output...
    print "Got an alignment, %i sequences" % (len(align.get_all_seqs()))
    #The length of the alignment will depend on the version of clustalw
    #(clustalw 2.0.10 and clustalw 1.83 are certainly different).
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

