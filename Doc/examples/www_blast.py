#!/usr/bin/env python
"""Example showing how to deal with internet BLAST from Biopython.

This code is described in great detail in the BLAST section of the Biopython
documentation.
"""
# standard library
import cStringIO

# biopython
from Bio.Blast import NCBIWWW
from Bio import Fasta

# first get the sequence we want to parse from a FASTA file
file_for_blast = open('m_cold.fasta', 'r')
f_iterator = Fasta.Iterator(file_for_blast)

f_record = f_iterator.next()

print 'Doing the BLAST and retrieving the results...'
result_handle = NCBIWWW.qblast('blastn', 'nr', f_record)

# save the results for later, in case we want to look at it
save_file = open('m_cold_blast.out', 'w')
blast_results = result_handle.read()
save_file.write(blast_results)
save_file.close()

print 'Parsing the results and extracting info...'
b_parser = NCBIWWW.BlastParser()

# option 1 -- parse the string directly
# b_record = b_parser.parse_str(blast_results)

# option 2 -- create a handle from the string and parse it
string_result_handle = cStringIO.StringIO(blast_results)
b_record = b_parser.parse(string_result_handle)

# now get the alignment info for all e values greater than some threshold
E_VALUE_THRESH = 0.1

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
    
    
