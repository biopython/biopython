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
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from Bio.Align import AlignInfo

# create the command line to run clustalw
# this assumes you've got clustalw somewhere on your path, otherwise
# you need to pass the full path of the executable to this via cmd="..."
cline = ClustalwCommandline(infile="opuntia.fasta", outfile="test.aln")

# actually perform the alignment
return_code = subprocess.call(str(cline), shell=(sys.platform != "win32"))
assert return_code == 0, "Calling ClustalW failed"

# Parse the output
alignment = AlignIO.read("test.aln", "clustal")

print(alignment)

print("first description: %s" % alignment[0].description)
print("first sequence: %s" % alignment[0].seq)

# get the length of the alignment
print("length %i" % alignment.get_alignment_length())

print(alignment)

# print out interesting information about the alignment
summary_align = AlignInfo.SummaryInfo(alignment)

consensus = summary_align.dumb_consensus()
print("consensus %s" % consensus)

my_pssm = summary_align.pos_specific_score_matrix(consensus, chars_to_ignore=["N"])
print(my_pssm)

expect_freq = {"A": 0.3, "G": 0.2, "T": 0.3, "C": 0.2}

info_content = summary_align.information_content(
    5, 30, chars_to_ignore=["N"], e_freq_table=expect_freq
)

print("relative info content: %f" % info_content)
