#!/usr/bin/env python
# Created: Sat Jul 21 08:01:13 2001
# Last changed: Time-stamp: <01/07/21 16:58:18 thomas>
# thomas@cbs.dtu.dk, # Cecilia.Alsmark@ebc.uu.se

# Copyright 2001 by Thomas Sicheritz-Ponten and Cecilia Alsmark.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
# File: sequtils.py

import os, sys
from Bio import Fasta
from Bio.Tools import Translate
from Bio.Seq import Seq
from Bio import Alphabet
from Bio.Alphabet import IUPAC
from Bio.Data import IUPACData, CodonTable
from PropertyManager import default_manager

"""

* sequtils is a temporary bucket for sequence utilities created by
  molecular biologs for molecular biologs.
  
* all functions can be used on multi FASTA files from the command line via the
  function multi_fasta_to.
  Usage: multi_fasta_to <name of FASTA file> <name of function to call>

* to get a list of functions type: sequtils.py describe

* functions starting with 'x' are graphical functions ... should definitely be   moved elsewhere ...

* most functions should probably be moved into better locations
    ... someday ... someone ...

* functions to add:
   - primer calculation and selection
   - melting point
   - pattern matcher
   more suggestions ?

"""

# should we provide a FASTA fasta reader for huge files ?
#  .... reading a complete genome file with the parser takes ages ...

# temporary hack for exception free translation of "dirty" DNA
# should be moved to ???

class ProteinX(Alphabet.ProteinAlphabet):
   letters = IUPACData.extended_protein_letters + "X"

proteinX = ProteinX()

class MissingTable:
  def __init__(self, table):
    self._table = table
  def get(self, codon, stop_symbol):
    try:
      return self._table.get(codon, stop_symbol)
    except CodonTable.TranslationError:
      return 'X'

def makeTableX(table):
  assert table.protein_alphabet == IUPAC.extended_protein
  return CodonTable.CodonTable(table.nucleotide_alphabet, proteinX,
                               MissingTable(table.forward_table),
                               table.back_table, table.start_codons,
                               table.stop_codons)



# end of hacks


def complement(seq):
   " returns the complementary sequence (NOT antiparallel) "
    return ''.join([IUPACData.ambiguous_dna_complement[x] for x in seq])

def reverse(seq):
   " reverse the sequence "
    r = map(None, seq)
    r.reverse()
    return ''.join(r)

def antiparallel(seq):
   " returns reversed complementary sequence ( = other strand ) "
    s = complement(seq)
    s = reverse(s)
    return s

def translate(seq, frame = 1, genetic_code = 1, translator = None):
   " translation of DNA in one of the six different reading frames "
    if frame not in [1,2,3,-1,-2,-3]:
        raise ValueError, 'invalid frame'

    if not translator:
        table = makeTableX(CodonTable.ambiguous_dna_by_id[genetic_code])
        translator = Translate.Translator(table)
        
    return translator.translate(Seq(seq[frame-1:], IUPAC.ambiguous_dna)).data

def GC(seq):
   " calculates G+C content "
    d = {}
    for nt in ['A','T','G','C']:
        d[nt] = seq.count(nt)
        gc = d.get('G',0) + d.get('C',0)
        
    if gc == 0: return 0
    return round(gc*100.0/(d['A'] +d['T'] + gc),1)
    
def GC123(seq):
   " calculates totla G+C content plus first, second and third position "
   l = len(seq)

   d= {}
   for nt in ['A','T','G','C']:
      d[nt] = [0,0,0]

   for i in range(0,l,3):
      codon = seq[i:i+3]
      if len(codon) <3: codon += '  '
      for pos in range(0,3):
         for nt in ['A','T','G','C']:
            if codon[pos] == nt: d[nt][pos] = d[nt][pos] +1


   gc = {}
   gcall = 0
   nall = 0
   for i in range(0,3):
      try:
         n = d['G'][i] + d['C'][i] +d['T'][i] + d['A'][i]
         gc[i] = (d['G'][i] + d['C'][i])*100.0/n
      except:
         gc[i] = 0

      gcall = gcall + d['G'][i] + d['C'][i]
      nall = nall + n

   gcall = 100.0*gcall/nall
   res = '%.1f%%, %.1f%%, %.1f%%, %.1f%%' % (gcall, gc[0], gc[1], gc[2])
   return res

def GC_skew(seq, window = 100):
   " calculates GC skew (G-C)/(G+C) "
   values = []
   for i in range(0, len(seq), window):
      s = seq[i: i + window]
      gc = GC(s)
      values.append(gc)
   return values

def Accumulated_GC_skew(seq, window = 1000):
   " calculates the accumulated GC skew (G-C)/(G+C) (easier to see)"
   values = []
   acc = 0
   for i in range(0, len(seq), window):
      s = seq[i: i + window]
      print GC(s)
      acc += GC(s)
      values.append(acc)
   return values

def xGC_skew(seq, window = 1000, zoom = 1,
                         r = 400, px = 100, py = 100):
   " calculates and plots the GC skew (GRAPHICS !!!) "
   
   from Tkinter import *
   from math import pi, sin, cos, log

   yscroll = Scrollbar(orient = VERTICAL)
   xscroll = Scrollbar(orient = HORIZONTAL)
   canvas = Canvas(yscrollcommand = yscroll.set,
                   xscrollcommand = xscroll.set, background = 'white')
   
   yscroll.config(command = canvas.yview)
   xscroll.config(command = canvas.xview)
   yscroll.pack(side = RIGHT, fill = Y)
   xscroll.pack(side = BOTTOM, fill = X)        
   canvas.pack(fill=BOTH, side = LEFT, expand = 1)
   canvas.update()

   X0, Y0  = r + px, r + py
   x1, x2, y1, y2 = X0 - r, X0 + r, Y0 -r, Y0 + r
   r1 = r
   
   canvas.create_oval(x1,y1, x2, y2)

   start = 0
   for gc in GC_skew(seq, window):
      alpha = pi - (2*pi*start)/len(seq)
      r2 = r1 - gc*zoom
      x1 = X0 + r1 * sin(alpha)
      y1 = Y0 + r1 * cos(alpha)
      x2 = X0 + r2 * sin(alpha)
      y2 = Y0 + r2 * cos(alpha)
      canvas.create_line(x1,y1,x2,y2, fill = 'green4')
      canvas.update()
      start = start + window
      

   canvas.configure(scrollregion = canvas.bbox(ALL))

def molecular_weight(seq):
   if type(seq) == type(''): seq = Seq(seq, IUPAC.unambiguous_dna)
   weight_table = IUPACData.unambiguous_dna_weights
   sum = 0
   for x in seq:
      sum += weight_table[x]
   return sum


def multi_fasta_to(file, function):
   " apply function on each sequence in a multiple FASTA file "
   try:
      f = globals()[function]
   except:
      raise NotImplementedError, "%s not implemented" % function
   
   parser = Fasta.RecordParser()
   handle = open(file, 'r')
   iter = Fasta.Iterator(handle, parser)

   while 1:
      record = iter.next()
      if not record: break
      print '>%s\n%s' % (record.title, f(record.sequence))



def test_all():
    string_seq = 'ATGACAAAGCTAATTATTCACTTGGTTTCAGACTCTTCTGTGCAAACTGCAAAACATGCAGCAAATTCTGCTCTTGCTCAATTTACTTCTATAAAACAAAAATTGTATCATTGGCCAATGATTAGAAATTGTGAATTACTAAATGAAGTATTAAGTAAAATAGAATCTAAACATGGAATAGTATTATACACAATTGCTGA'
    seq = Seq(string_seq, IUPAC.ambiguous_dna)
    
    print 'GC(Seq): ',GC(seq), 'GC(string):', GC(string_seq)
    print 'GC123(Seq):',GC123(seq), 'GC123(string):', GC123(string_seq)
    print 'seq:         ', seq.data
    print 'complement:  ', complement(seq)
    print 'reverse:     ', reverse(seq)
    print 'antiparallel:', antiparallel(seq)
    print 'comp+rev    :', complement(reverse(seq))




    
if __name__ == '__main__':

   if sys.argv[1] == 'describe':
      # get all new functions from this file
      mol_funcs = [x[0] for x in locals().items() if type(x[1]) == type(GC)]
      mol_funcs.sort()
      print 'available functions:'
      for f in mol_funcs: print '\t', f
      sys.exit(0)

   
   file = sys.argv[1]
   function = sys.argv[2]
   multi_fasta_to(file, function)

      
