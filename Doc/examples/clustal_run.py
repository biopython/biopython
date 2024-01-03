#!/usr/bin/env python
# Copyright 2000 Brad Chapman.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Run clustalw and parse the output.

Example code to show how to create a clustalw command line, run clustalw
and parse the results into an object that can be dealt with easily.
"""
# standard library


import sys
import subprocess

# biopython
from Bio import Align
from Bio.motifs import Motif

# create the command line to run clustalw
# this assumes you've got clustalw somewhere on your path, otherwise
# you need to pass the full path of the executable to this via cmd="..."
command = "clustalw -infile=opuntia.fasta -outfile=test.aln"

# actually perform the alignment
return_code = subprocess.call(command, shell=(sys.platform != "win32"))
assert return_code == 0, "Calling ClustalW failed"

# Parse the output
alignment = Align.read("test.aln", "clustal")

print(alignment)

print("first description: %s" % alignment.sequences[0].description)
print("first sequence: %s" % alignment.sequences[0].seq)

# get the length of the alignment
print("length %i" % alignment.length)

motif = Motif("ACGT", alignment)
counts = motif.counts

# print out interesting information about the alignment
consensus = motif.counts.calculate_consensus(identity=0.7)
print("consensus %s" % consensus)

motif.background = {"A": 0.3, "G": 0.2, "T": 0.3, "C": 0.2}
relative_entropy = motif.relative_entropy
print("relative entropy for columns [5:30]: %f" % sum(relative_entropy[5:30]))
