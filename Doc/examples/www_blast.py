#!/usr/bin/env python
# Copyright 2000 Brad Chapman.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Example showing how to deal with internet BLAST from Biopython.

This code is described in great detail in the BLAST section of the Biopython
documentation.
"""
# standard library
from io import StringIO

# biopython
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# first get the sequence we want to parse from a FASTA file
f_record = next(SeqIO.parse("m_cold.fasta", "fasta"))

print("Doing the BLAST and retrieving the results...")
result_handle = NCBIWWW.qblast("blastn", "nr", f_record.format("fasta"))

# save the results for later, in case we want to look at it
with open("m_cold_blast.out", "w") as save_file:
    blast_results = result_handle.read()
    save_file.write(blast_results)

print("Parsing the results and extracting info...")

# option 1 -- open the saved file to parse it
# option 2 -- create a handle from the string and parse it
string_result_handle = StringIO(blast_results)
b_record = NCBIXML.read(string_result_handle)

# now get the alignment info for all e values greater than some threshold
E_VALUE_THRESH = 0.1

for alignment in b_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print("****Alignment****")
            print("sequence: %s" % alignment.title)
            print("length: %i" % alignment.length)
            print("e value: %f" % hsp.expect)
            print(hsp.query[0:75] + "...")
            print(hsp.match[0:75] + "...")
            print(hsp.sbjct[0:75] + "...")
