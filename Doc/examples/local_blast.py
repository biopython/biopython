#!/usr/bin/env python
"""Script demonstrating the ability to interact with local BLAST.

The contents of this script are described more fully in the available
documentation.
"""
# standard library
import os
import sys

# biopython
from Bio.Blast import NCBIStandalone

my_blast_db = os.path.join(os.getcwd(), 'at-est', 'a_cds-10-7.fasta')
my_blast_file = os.path.join(os.getcwd(), 'at-est', 'test_blast',
                             'sorghum_est-test.fasta')
my_blast_exe = os.path.join(os.getcwd(), 'blast', 'blastall')

print 'Running blastall...'
blast_out, error_info = NCBIStandalone.blastall(my_blast_exe, 'blastn',
                                                my_blast_db, my_blast_file)


b_parser = NCBIStandalone.BlastParser()

b_iterator = NCBIStandalone.Iterator(blast_out, b_parser)

while 1:
    b_record = b_iterator.next()

    if b_record is None:
        break
    
    E_VALUE_THRESH = 0.04
    for alignment in b_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                print '****Alignment****'
                print 'sequence:', alignment.title
                print 'length:', alignment.length
                print 'e value:', hsp.expect

                if len(hsp.query) > 75:
                    dots = '...'
                else:
                    dots = ''
                
                print hsp.query[0:75] + dots
                print hsp.match[0:75] + dots
                print hsp.sbjct[0:75] + dots


