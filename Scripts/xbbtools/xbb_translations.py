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
from Bio.Alphabet import IUPAC
from Bio import Translate
from Bio.Data import IUPACData


class xbb_translations:
    def __init__(self):
        ""

    def frame1(self, seq, translation_table = 1):
        dna = Seq.Seq(seq, IUPAC.unambiguous_dna)
        trans = Translate.unambiguous_dna_by_id[translation_table]
        protein = trans.translate(dna)
        return protein

    def complement(self, seq):
        return string.join(map(lambda x:IUPACData.ambiguous_dna_complement[x], map(None,seq)),'')

    def reverse(self, seq):
        r = map(None, seq)
        r.reverse()
        return string.join(r,'')

    def antiparallel(self, seq):
        s = self.complement(seq)
        s = self.reverse(s)
        return s
    
    def frame(self, seq, frame, translation_table = 1):
        dna = Seq.Seq(seq, IUPAC.unambiguous_dna)
        trans = Translate.unambiguous_dna_by_id[translation_table]
        if not ((-3 <= frame <= -1) or (1 <= frame <= 3)):
            frame = 1
        protein = trans.translate(dna)
        return protein

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
        ngc = string.count(seq,'G') + string.count(seq,'C')
        if ngc == 0: return 0.0
        gc = (100.0*ngc)/len(seq)
        return gc
    
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
        
