#!/usr/bin/env python
# Created: Tue Aug  8 20:32:36 2000
# Last changed: Time-stamp: <00/08/10 19:03:53 thomas>
# thomas@cbs.dtu.dk, http://www.cbs.dtu.dk/thomas/index.html
# File: nextorf.py


import string, re, regsub
import posixpath, posix
import os, sys, commands
import getopt

from Bio.Fasta import Fasta
from Bio.Tools import Translate
from Bio.Seq import Seq
from Bio import Alphabet
from Bio.Alphabet import IUPAC
from Bio.Data import IUPACData, CodonTable

# From: "Andrew Dalke" <dalke@acm.org>
# To: <thomas@cbs.dtu.dk>, <biopython-dev@biopython.org>
# Subject: Re: [Biopython-dev] unambiguous DNA
# Date: Thu, 10 Aug 2000 01:43:07 -0600
# [..] There isn't an alphabet in biopython which supports output sequence
# uses the ambiguous protein alphabet along with '*' for stop codon
# symbol and 'X' for untranslatable protein residues , so we need to
# make a new one:

class ProteinX(Alphabet.ProteinAlphabet):
   letters = IUPACData.extended_protein_letters + "X"

proteinX = ProteinX()

# Forward translation table, mapping codon to protein
class MissingTable:
  def __init__(self, table):
    self._table = table
  def get(self, codon, stop_symbol):
    try:
      return self._table.get(codon, stop_symbol)
    except CodonTable.TranslationError:
      return 'X'

# Make the codon table given an existing table
def makeTableX(table):
  assert table.protein_alphabet == IUPAC.extended_protein
  return CodonTable.CodonTable(table.nucleotide_alphabet, proteinX,
                               MissingTable(table.forward_table),
                               table.back_table, table.start_codons,
                               table.stop_codons)

def complement(seq):
    return string.join(map(lambda x:IUPACData.ambiguous_dna_complement[x], map(None,seq)),'')

class NextORF:
    def __init__(self, file, options):
        self.code = 1
        self.file = file
        self.options = options
        table = makeTableX(CodonTable.ambiguous_dna_by_id[int(self.options['table'])])
        self.translator = Translate.Translator(table)
        
    def read_file(self):
        self.parser = Fasta.RecordParser()
        self.iter = Fasta.Iterator(handle = open(self.file), parser = self.parser)
        while 1:
            rec = self.iter.next()
            if not rec: break
            self.header = string.split(rec.title,',')[0]
            self.handle_record(rec)

    def toFasta(self, header, seq):
       seq = re.sub('(............................................................)','\\1\n',seq)
       return '>%s\n%s' % (header, seq)

    def gc(self, seq):
       d = {}
       for nt in ['A','T','G','C']:
          d[nt] = string.count(seq, nt)
       gc = d['G'] + d['C']
       if gc == 0: return 0
       return round(gc*100.0/(d['A'] +d['T'] + gc),1)

    def gc2(self,seq):
       l = len(seq)
       d= {}
       for nt in ['A','T','G','C']:
          d[nt] = [0,0,0]
          
       for i in range(0,l,3):
          codon = seq[i:i+3]
          for pos in range(0,3):
             for nt in ['A','T','G','C']:
                if codon[0] == nt: d[nt][pos] = d[nt][pos] +1


       gc = {}
       gcall = 0
       nall = 0
       for i in range(0,3):
          n = (d['G'][1] + d['C'][1] +d['T'][1] + d['A'][1])
          nall = nall + n
          if n == 0:
             gc[i] = 0
          else:
             gc[i] = (d['G'][1] + d['C'][1])*100.0/n
          gcall = gcall + gc[i]

       gcall = 100.0*gcall/nall
       return '%.1f%%, %.1f%%, %.1f%%, %.1f%%' % (gcall, gc[0], gc[1], gc[1])
          
   
    def handle_record(self, rec):
        plus, minus= 0,0
        if self.options['strand'] == 'both' or self.options['strand'] == 'plus': plus = 1
        if self.options['strand'] == 'both' or self.options['strand'] == 'minus': minus = 1         
        s = string.upper(rec.sequence[self.options['start']:self.options['stop']])
        seq = Seq(s,IUPAC.ambiguous_dna)
        if minus:
            r = map(None,s)
            r.reverse()
            rseq = Seq(complement(r), IUPAC.ambiguous_dna)
            length = len(s)

        n = 0
        if plus:
            for frame in self.options['frames']:
                orf = self.translator.translate(seq[frame-1:])
                orfs = string.split(orf.data,'*')
                start = 0
                for orf in orfs:
                    stop = start + 3*(len(orf) + 1)
                    subs = seq[frame-1:][start:stop]
                    _start = start
                    start = stop

                    # ORF too small ?
                    if len(orf) < int(self.options['minlength']): continue
                    # ORF too big ?
                    if self.options['maxlength'] and \
                       len(orf) > int(self.options['maxlength']): continue
                    
                    n = n + 1
                    # ORF just allright ...
                    out = self.options['output']
                    head = 'orf_%s:%s:+%d:%d:%d' % (n, self.header, frame, _start+1,stop+1)
                    if self.options['gc']: head = '%s:%s' % (head, self.gc2(subs.data))
                    if out == 'aa':
                        print self.toFasta(head, orf)
                    elif out == 'nt':
                        print self.toFasta(head, subs.data)
                    elif out == 'pos':
                        print head
                        
                        

        if minus:
            for frame in self.options['frames']:
                orf = self.translator.translate(rseq[frame-1:])
                orfs = string.split(orf.data,'*')
                start = 0
                for orf in orfs:
                    stop = start + 3*(len(orf) + 1)
                    subs = rseq[frame-1:][start:stop]
                    _start = start
                    start = stop

                    # ORF too small ?
                    if len(orf) < int(self.options['minlength']): continue
                    # ORF too big ?
                    if self.options['maxlength'] and \
                       len(orf) > int(self.options['maxlength']): continue
                    # ORF just allright ...
                    n = n + 1
                    head = 'orf_%s:%s:-%d:%d:%d' % (n, self.header, frame, length - stop +2, length - _start)
                    if self.options['gc']: head = '%s:%s' % (head, self.gc2(subs.data))
                    out = self.options['output']
                    if out == 'aa':
                        print self.toFasta(head, orf)
                    elif out == 'nt':
                        print self.toFasta(head, subs.data)
                    elif out == 'pos':
                        print head

            

def help():
    global options
    print 'Usage:', sys.argv[0], '<FASTA file> (<options>)'
    
    print 'Internal Options:'
    for a,b in options.items(): print '\t', a,b
    print ''
    print "NCBI's Codon Tables:"
    for key, table in CodonTable.ambiguous_dna_by_id.items():
        print '\t',key, table._codon_table.names[0]
        
    sys.exit(0)
    

options = {
    'start': 0,
    'stop': -1,
    'minlength': 100,
    'maxlength': None,
    'strand': 'both',
    'output': 'aa',
    'frames': [1,2,3],
    'gc': None,
    'cc': None,
    'nostart': None,
    'table': 1,
    }

if __name__ == '__main__':

    args = sys.argv[1:]

    show_help = len(sys.argv)<=1

    shorts = 'hv'
    longs = map(lambda x: x +'=', options.keys())
    optlist, args = getopt.getopt(args,shorts, longs)
    if show_help: help()
    for arg in optlist:
        if arg[0] == '-h' or arg[0] == '--help':
            help()
            sys.exit(0)
        for key in options.keys():
            if arg[0][2:] == key:
                if arg[1] == 'no': 
                    options[key] = None
                else:
                    options[key] = arg[1]
        if arg[0] == 'v':
            print options
    file = args[0]

    nextorf = NextORF(file, options)
    nextorf.read_file()

    
