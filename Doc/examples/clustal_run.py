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
from Bio.motifs import Motif

# create the command line to run clustalw
# this assumes you've got clustalw somewhere on your path, otherwise
# you need to pass the full path of the executable to this via cmd="..."
cline = ClustalwCommandline(infile="opuntia.fasta", outfile="test.aln")

# actually perform the alignment
return_code = subprocess.call(str(cline), shell=(sys.platform != "win32"))
assert return_code == 0, "Calling ClustalW failed"

# Parse the output
msa = AlignIO.read("test.aln", "clustal")

print(msa)

print("first description: %s" % msa[0].description)
print("first sequence: %s" % msa[0].seq)

# get the length of the alignment
print("length %i" % msa.get_alignment_length())

print(msa)

# print out interesting information about the alignment
summary_align = AlignInfo.SummaryInfo(msa)

consensus = summary_align.dumb_consensus()
print("consensus %s" % consensus)

alignment = msa.alignment
motif = Motif("ACGT", alignment)
counts = motif.counts
print(counts)

expect_freq = {"A": 0.3, "G": 0.2, "T": 0.3, "C": 0.2}

info_content = summary_align.information_content(
    5, 30, chars_to_ignore=["N"], e_freq_table=expect_freq
)

print("relative info content: %f" % info_content)
