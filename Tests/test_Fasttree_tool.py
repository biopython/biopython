# Copyright 2013 by Nate Sutton.  All rights reserved.
# Based on test_Clustalw_tool.py by Peter Cock.
# Example code used from Biopython's Phylo cookbook by Eric Talevich.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from __future__ import print_function

from Bio import MissingExternalDependencyError

import sys
import os
import itertools
from Bio._py3k import StringIO
from Bio._py3k import zip

from Bio import SeqIO
from Bio import Phylo
from Bio.Phylo.Applications import FastTreeCommandline
from Bio.Phylo.Applications import _Fasttree
from Bio.Application import ApplicationError

#################################################################

#Try to avoid problems when the OS is in another language
os.environ['LANG'] = 'C'

fasttree_exe = None
if sys.platform == "win32":
    try:
        #This can vary depending on the Windows language.
        prog_files = os.environ["PROGRAMFILES"]
    except KeyError:
        prog_files = r"C:\Program Files (x86)"

    # A default fasttree file path of "C:\Program Files (x86)\Fasttree.exe"
    # was chosen here but users can alter the path according to where
    # fasttree is located on their systems

    likely_dirs = ["Fasttree", "fasttree", "FastTree"]
    likely_exes = ["Fasttree.exe", "fasttree.exe", "FastTree.exe"]
    for folder in likely_dirs:
        if os.path.isdir(os.path.join(prog_files, folder)):
            for filename in likely_exes:
                if os.path.isfile(os.path.join(prog_files, folder, filename)):
                    fasttree_exe = os.path.join(prog_files, folder, filename)
                    break
            if fasttree_exe:
                break
else:
    from Bio._py3k import getoutput
    # Checking the -help argument
    output = getoutput("fasttree -help")
    # Since "is not recognized" may be in another language, try and be sure this
    # is really the fasttree tool's output
    fasttree_found = False
    if "is not recognized" not in output and "protein_alignment" in output \
    and "nucleotide_alignment" in output:
        fasttree_exe = "fasttree"

if not fasttree_exe:
    raise MissingExternalDependencyError(
        "Install fasttree and correctly set the file path to the program "
        "if you want to use it from Biopython.")

#################################################################

print("Checking error conditions")
print("=========================")

print("Empty file")
input_file = "does_not_exist.fasta"
assert not os.path.isfile(input_file)
cline = FastTreeCommandline(fasttree_exe, input=input_file)
try:
    stdout, stderr = cline()
    assert False, "Should have failed, returned:\n%s\n%s" % (stdout, stderr)
except ApplicationError as err:
    print("Failed (good)")
    #Python 2.3 on Windows gave (0, 'Error')
    #Python 2.5 on Windows gives [Errno 0] Error
    assert "Cannot open sequence file" in str(err) or \
           "Cannot open input file" in str(err) or \
           "non-zero exit status" in str(err), str(err)

print("")
print("Single sequence")
input_file = "Fasta/f001"
assert os.path.isfile(input_file)
assert len(list(SeqIO.parse(input_file, "fasta")))==1
cline = FastTreeCommandline(fasttree_exe, input=input_file)
try:
    stdout, stderr = cline()
    if "Unique: 1/1" in stderr:
        print("Failed (good)")
    else:
        assert False, "Should have failed, returned:\n%s\n%s" % (stdout, stderr)
except ApplicationError as err:
    print("Failed (good)")
    #assert str(err) == "No records found in handle", str(err)

print("")
print("Invalid sequence")
input_file = "Medline/pubmed_result1.txt"
assert os.path.isfile(input_file)
cline = FastTreeCommandline(fasttree_exe, input=input_file)
try:
    stdout, stderr = cline()
    assert False, "Should have failed, returned:\n%s\n%s" % (stdout, stderr)
except ApplicationError as err:
    print("Failed (good)")
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
print("")
print("Checking normal situations")
print("==========================")

#Create a temp fasta file with a space in the name
temp_filename_with_spaces = "Clustalw/temp horses.fasta"
handle = open(temp_filename_with_spaces, "w")
SeqIO.write(SeqIO.parse("Phylip/hennigian.phy", "phylip"), handle, "fasta")
handle.close()

for input_file in ["Quality/example.fasta", "Clustalw/temp horses.fasta"]:
    input_records = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))
    print("")
    print("Calling fasttree on %s (with %i records)" \
          % (repr(input_file), len(input_records)))

    #Any filesnames with spaces should get escaped with quotes automatically.
    #Using keyword arguments here.
    cline = _Fasttree.FastTreeCommandline(fasttree_exe, input=input_file)
    assert str(eval(repr(cline)))==str(cline)

    out, err = cline()
    assert err.strip().startswith("FastTree")

    print("")
    print("Checking generation of tree terminals")
    tree = Phylo.read(StringIO(out), 'newick')

    def lookup_by_names(tree):
        names = {}
        for clade in tree.find_clades():
            if clade.name:
                if clade.name in names:
                    raise ValueError("Duplicate key: %s" % clade.name)
                names[clade.name] = clade
        return names

    names = lookup_by_names(tree)

    assert len(names) > 0.0
    print("Success")

    print("")
    print("Checking distances between tree terminals")
    def terminal_neighbor_dists(self):
        """Return a list of distances between adjacent terminals."""
        def generate_pairs(self):
            pairs = itertools.tee(self)
            next(pairs[1]) # Advance second iterator one step
            return zip(pairs[0], pairs[1])
        return [self.distance(*i) for i in
                generate_pairs(self.find_clades(terminal=True))]

    for dist in terminal_neighbor_dists(tree):
        assert dist > 0.0

    print("Success")

print("")
print("Done")
