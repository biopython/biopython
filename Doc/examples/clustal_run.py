#!/usr/bin/env python
"""clustal_run.py

Example code to show how to create a clustalw command line, run clustalw
and parse the results into an object that can be dealt with easily."""
# standard library
import os

# biopython
from Bio.Clustalw import Clustalw
from Bio.Clustalw.Clustalw import MultipleAlignCL

# create the command line to run clustalw
# this assumes you've got clustalw somewhere on your path, otherwise
# you need to pass a second argument to MultipleAlignCL with the complete
# path to clustalw
cline = MultipleAlignCL(os.path.join(os.curdir, 'opuntia.fasta'))
cline.set_output('test.aln')

# actually perform the alignment and get back an alignment object
alignment = Clustalw.do_alignment(cline)

# get the records in the alignment
all_records = alignment.get_all_seqs()

print 'description:', all_records[0].description
print 'sequence:', all_records[0].seq

# get the length of the alignment
print 'length', alignment.get_alignment_length()

# print out a consensus
print 'consensus', alignment.dumb_consensus()

print alignment
