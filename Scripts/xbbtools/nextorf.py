#!/usr/bin/env python
# Created: Thu Feb 15 14:22:12 2001
# Last changed: Time-stamp: <01/02/15 17:09:02 thomas>
# thomas@cbs.dtu.dk, http://www.cbs.dtu.dk/thomas
# File: nextorf.py

import string, re
import os, sys, commands
import getopt

from Bio import Fasta
from Bio.Tools import Translate
from Bio.Seq import Seq
from Bio import Alphabet
from Bio.Alphabet import IUPAC
from Bio.Data import IUPACData, CodonTable

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

# Make the codon table given an existing table
def makeTableX(table):
  assert table.protein_alphabet == IUPAC.extended_protein
  return CodonTable.CodonTable(table.nucleotide_alphabet, proteinX,
                               MissingTable(table.forward_table),
                               table.back_table, table.start_codons,
                               table.stop_codons)

def complement(seq):
    return string.join(map(lambda x:IUPACData.ambiguous_dna_complement[x], map(None,seq)),'')

def reverse(seq):
    r = map(None, seq)
    r.reverse()
    return string.join(r,'')

def antiparallel(seq):
    s = complement(seq)
    s = reverse(s)
    return s

class NextOrf:
    def __init__(self, file, options):
        self.options = options
        self.file = file
        self.genetic_code = int(self.options['table'])
        self.table = makeTableX(CodonTable.ambiguous_dna_by_id[self.genetic_code])
        self.translator = Translate.Translator(self.table)
        self.counter = 0
        self.ReadFile()
        
    def ReadFile(self):
        self.parser = Fasta.RecordParser()
        self.iter = Fasta.Iterator(handle = open(self.file), parser = self.parser)
        while 1:
            rec = self.iter.next()
            if not rec: break
            self.header = rec.title.split()[0].split(',')[0]
            self.HandleRecord(rec)

    def ToFasta(self, header, seq):
       seq = re.sub('(............................................................)','\\1\n',seq)
       return '>%s\n%s' % (header, seq)

    def Gc(self, seq):
       d = {}
       for nt in ['A','T','G','C']:
          d[nt] = string.count(seq, nt)
       gc = d['G'] + d['C']
       if gc == 0: return 0
       return round(gc*100.0/(d['A'] +d['T'] + gc),1)

    def Gc2(self,seq):
       l = len(seq)
       d= {}
       for nt in ['A','T','G','C']:
          d[nt] = [0,0,0]
          
       for i in range(0,l,3):
          codon = seq[i:i+3]
          if len(codon) <3: codon = codon + '  '
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
       print 'GC:', res
       return res
          

    def GetOrfCoordinates(self, seq):
        s = seq.data
        letters = []
        table = self.table
        get = self.table.forward_table.get
        n = len(seq)
        start_codons = self.table.start_codons
        stop_codons = self.table.stop_codons
        frame_coordinates = []
        for frame in range(0,3):
            coordinates = []
            for i in range(0+frame, n-n%3, 3):
                codon = s[i:i+3]
                if codon in start_codons: coordinates.append((i,1,codon))
                elif codon in stop_codons: coordinates.append((i,0,codon))


            frame_coordinates.append(coordinates)

        return frame_coordinates

    def HandleRecord(self, rec):
        frame_coordinates = ''
        dir = self.options['strand']
        plus = dir in ['both', 'plus']
        minus = dir in ['both', 'minus']

        start, stop = int(self.options['start']), int(self.options['stop'])
        s = string.upper(rec.sequence[start:stop])
        self.seq = Seq(s,IUPAC.ambiguous_dna)
        self.length = len(self.seq)
        self.rseq = None

        
        CDS = []
        if plus: CDS.extend(self.GetCDS(self.seq))
        if minus:
            self.rseq = Seq(antiparallel(s),IUPAC.ambiguous_dna)
            CDS.extend(self.GetCDS(self.rseq, strand = -1))

        self.Output(CDS)
            
        

    def GetCDS(self, seq, strand = 1):
        frame_coordinates = self.GetOrfCoordinates(seq)
        START, STOP = 1,0
        so = self.options
        nostart = so['nostart']
        minlength, maxlength = int(so['minlength']), int(so['maxlength'])
        CDS = []
        f = 0
        for frame in frame_coordinates:
            f+=1
            start_site = -1
            if nostart: start_site = 0
            frame.append((self.length, 0, 'XXX'))
            for pos, codon_type, codon in frame:
                if codon_type == START:
                    if start_site == -1:start_site = pos
                elif codon_type == STOP:
                    if start_site == -1: continue
                    stop = pos
                    length = stop - start_site +1
                    if length >= minlength and length <= maxlength:
                        s = seq[start_site:stop +1]
                        CDS.append((start_site, stop, length, s, strand*f))
                        start_site = -1
                        if nostart: start_site = pos + 3
                        del stop

        return CDS
    

    def Output(self, CDS):
        out = self.options['output']
        seqs = (self.seq, self.rseq)
        
        for start, stop, length, subs, strand in CDS:
            self.counter += 1
            head = 'orf_%s:%s:%d:%d:%d' % (self.counter, self.header, strand, start,stop+1)

            if self.options['gc']:
                head = '%s:%s' % (head, self.Gc2(subs.data))
                
            if out == 'aa':
                orf = self.translator.translate(subs)
                print self.ToFasta(head, orf.data)
            elif out == 'nt':
                print self.ToFasta(head, subs.data)
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

    print 'e.g.\n./nextorf.py --minlength 5 --strand plus --output nt --gc 1 testjan.fas'
    sys.exit(0)
    

options = {
    'start': 0,
    'stop': -1,
    'minlength': 100,
    'maxlength': 100000000,
    'strand': 'both',
    'output': 'aa',
    'frames': [1,2,3],
    'gc': 0,
    'nostart': 0,
    'table': 1,
    }

if __name__ == '__main__':
    args = sys.argv[1:]
    show_help = len(sys.argv)<=1

    shorts = 'hv'
    longs = map(lambda x: x +'=', options.keys()) + ['help']

    optlist, args = getopt.getopt(args,shorts, longs)
    if show_help: help()

    for arg in optlist:
        if arg[0] == '-h' or arg[0] == '--help':
            help()
            sys.exit(0)
        for key in options.keys():
            if arg[1].lower() == 'no': arg[1] = 0
            elif arg[1].lower() == 'yes': arg[1] = 1
            if arg[0][2:] == key: options[key] = arg[1]
            
        if arg[0] == '-v':print 'OPTIONS', options

    file = args[0]
    nextorf = NextOrf(file, options)
