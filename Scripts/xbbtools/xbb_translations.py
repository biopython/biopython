#!/usr/bin/env python
# Created: Wed Jun 21 15:53:22 2000
# Last changed: Time-stamp: <00/12/02 15:56:27 thomas>
# thomas@cbs.dtu.dk, http://www.cbs.dtu.dk/thomas
# File: xbb_translations.py

import string
import os, sys  # os.system, sys.argv
import time

sys.path.insert(0, '.')
from Tkinter import *

from Bio import Seq
from Bio import Alphabet
from Bio.Seq import reverse_complement, translate
from Bio.SeqUtils import GC

class xbb_translations:
    def __init__(self):
        ""

    def frame1(self, seq, translation_table = 1):
        return translate(seq, table=translation_table)

    def complement(self, seq):
        #TODO - use Seq methods instead of this hack:?
        return reverse_complement(seq)[::-1]

    def reverse(self, seq):
        return seq[::-1]

    def antiparallel(self, seq):
        return reverse_complement(seq)
    
    def frame(self, seq, frame, translation_table = 1):
        if not ((-3 <= frame <= -1) or (1 <= frame <= 3)):
            frame = 1
        if frame != 1 :
            raise NotImplementedError
            #TODO - Support the frame argument
            #The old code didn't, but I can guess from
            #the code the expected 1,2,3 for the forward
            #strands and -1,-2,-3 for the reverse.
        return translate(seq, table=translation_table)

    def header_nice(self, txt, seq):
        length = len(seq)
        if length > 20:
            short = '%s ... %s' % (seq[:10], seq[-10:])
        else:
            short = seq
            
        date = time.strftime('%y %b %d, %X', time.localtime(time.time()))
        res = '%s: %s, ' % (txt,date)

        for nt in ['a','t','g','c']:
            res = res + '%s:%d ' % (nt, string.count(seq, string.upper(nt)))

        res = res + '\nSequence: %s, %d nt, %0.2f %%GC\n' % (string.lower(short),length, self.gc(seq))       
        res = res + '\n\n'
        return res
        
    def frame_nice(self, seq, frame, translation_table = 1):
        length = len(seq)
        protein = self.frame(seq, frame, translation_table)
        res = self.header_nice('Plus one frame translation',seq)
        for i in range(0,length,60):
            subseq = seq[i:i+60]
            p = i/3
            res = res + '%d/%d\n' % (i+1, i/3+1)
            res = res + string.join(map(None,protein[p:p+20]),'  ') + '\n'
            # seq
            res = res + string.lower(subseq) + '%5d %%\n' % int(self.gc(subseq))

        return res
    
    def gc(self, seq):
        """Returns a float between 0 and 100."""
        return GC(seq)
    
    def gcframe(self, seq, translation_table = 1):
        # always use uppercase nt-sequence !!
        comp = self.complement(seq)
        anti = self.reverse(comp)
        length = len(seq)
        frames = {}
        for i in range(0,3):
            #print i+1, seq[i:]
            frames[i+1]  = self.frame1(seq[i:], translation_table)
            #print -(i+1), anti[i:]
            frames[-(i+1)] = self.reverse(self.frame1(anti[i:], translation_table))
            #print len(frames[i+1])

        res = self.header_nice('GCFrame', seq)
#         if length > 20:
#             short = '%s ... %s' % (seq[:10], seq[-10:])
#         else:
#             short = seq
            
#         date = time.strftime('%y %b %d, %X', time.localtime(time.time()))
#         res = 'GCFrame: %s, ' % date

#         for nt in ['a','t','g','c']:
#             res = res + '%s:%d ' % (nt, string.count(seq, string.upper(nt)))

#         res = res + '\nSequence: %s, %d nt, %0.2f %%GC\n' % (string.lower(short),length, self.gc(seq))       
#         res = res + '\n\n'

        for i in range(0,length,60):
            subseq = seq[i:i+60]
            csubseq = comp[i:i+60]
            p = i/3
#             print 3, frames[3][p:p+20]
#             print 2, frames[2][p:p+20]
#             print 1, frames[1][p:p+20]
#             print -1, frames[-1][p:p+20]
#             print -2, frames[-2][p:p+20]
#             print -3,frames[-3][p:p+20]
            # + frames
            res = res + '%d/%d\n' % (i+1, i/3+1)
            res = res + '  ' + string.join(map(None,frames[3][p:p+20]),'  ') + '\n'
            res = res + ' ' + string.join(map(None,frames[2][p:p+20]),'  ') + '\n'
            res = res + string.join(map(None,frames[1][p:p+20]),'  ') + '\n'
            # seq
            res = res + string.lower(subseq) + '%5d %%\n' % int(self.gc(subseq))
            res = res + string.lower(csubseq) + '\n'
            # - frames
            res = res + string.join(map(None,frames[-2][p:p+20]),'  ')  +' \n'
            res = res + ' ' + string.join(map(None,frames[-1][p:p+20]),'  ') + '\n'
            res = res + '  ' + string.join(map(None,frames[-3][p:p+20]),'  ') + '\n\n'
            
            
        return res
        
if __name__ == '__main__':
    #s = 'GCCCTTTCTTATTAGTGCTACCGCTAATAGGTAAATATGAAAAACCTTTG'
    s = 'ATTCCGGTTGATCCTGCCGGACCCGACCGCTATCGGGGTAGGGATAAGCCATGGGAGTCTTACACTCCCGGGTAAGGGAGTGTGGCGGACGGCTGAGTAACACGTGGCTAACCTACCCTCGGGACGGGGATAACCCCGGGAAACTGGGGATAATCCCCGATAGGGAAGGAGTCCTGGAATGGTTCCTTCCCTAAAGGGCTATAGGCTATTTCCCGTTTGTAGCCGCCCGAGGATGGGGCTACGGCCCATCAGGCTGTCGGTGGGGTAAAGGCCCACCGAACCTATAACGGGTAGGGGCCGTGGAAGCGGGAGCCTCCAGTTGGGCACTGAGACAAGGGCCCAGGCCCTACGGGGCGCACCAGGCGCGAAACGTCCCCAATGCGCGAAAGCGTGAGGGCGCTACCCCGAGTGCCTCCGCAAGGAGGCTTTTCCCCGCTCTAAAAAGGCGGGGGAATAAGCGGGGGGCAAGTCTGGTGTCAGCCGCCGCGGTAATACCAGCTCCGCGAGTGGTCGGGGTGATTACTGGGCCTAAAGCGCCTGTAGCCGGCCCACCAAGTCGCCCCTTAAAGTCCCCGGCTCAACCGGGGAACTGGGGGCGATACTGGTGGGCTAGGGGGCGGGAGAGGCGGGGGGTACTCCCGGAGTAGGGGCGAAATCCTTAGATACCGGGAGGACCACCAGTGGCGGAAGCGCCCCGCTA'

    test = xbb_translations()
#     for i in range(0,4):
#         print test.frame1(s[i:])
    #print s
    #print test.complement(s)
    print '============================================================'
    print test.gcframe(s)
    
#     for i in Translate.unambiguous_dna_by_id.keys():
#         print Translate.unambiguous_dna_by_id[i].table.names[0]
        
