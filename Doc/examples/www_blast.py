#!/usr/bin/env python
"""Example showing how to deal with internet BLAST from Biopython.

This code is described in great detail in the BLAST section of the Biopython
documentation."""
# standard library
import copy

# biopython
from Bio.Blast import NCBIWWW
from Bio.Fasta import Fasta

# first get the sequence we want to parse from a FASTA file
file_for_blast = open('m_cold.fasta', 'r')
f_iterator = Fasta.Iterator(file_for_blast)

f_record = f_iterator.next()

print 'Doing the BLAST and retrieving the results...'
b_results = NCBIWWW.blast('blastn', 'nr', f_record)

# save the results for later
save_file = open('m_cold_blast.out', 'w')
# copy the handle so we can write it out and still have something for parsing
save_handle = copy.deepcopy(b_results)
save_file.write(save_handle.read())
save_file.close()

print 'Parsing the results and extracting info...'
b_parser = NCBIWWW.BlastParser()
b_record = b_parser.parse(b_results)

# now get the alignment info for all e values greater than some threshold
E_VALUE_THRESH = 0.04

for alignment in b_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print '****Alignment****'
            print 'sequence:', alignment.title
            print 'length:', alignment.length
            print 'e value:', hsp.expect
            print hsp.query[0:75] + '...'
            print hsp.match[0:75] + '...'
            print hsp.sbjct[0:75] + '...'
    
    
