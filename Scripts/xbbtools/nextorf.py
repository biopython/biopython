#!/usr/bin/env python
# Created: Tue Aug  8 20:32:36 2000
# Last changed: Time-stamp: <00/08/10 09:55:44 thomas>
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

    def gc(self,seq):
        gc = string.count(seq, 'G') + string.count(seq, 'C')
        if gc == 0: return 0
        return len(seq)*100.0/gc
        
    def handle_record(self, rec):
        plus, minus= 0,0
        if self.options['strand'] == 'both' or self.options['strand'] == 'plus': plus = 1
        if self.options['strand'] == 'both' or self.options['strand'] == 'minus': minus = 1         
        s = string.upper(rec.sequence[self.options['start']:self.options['stop']])
        seq = Seq(s,IUPAC.ambiguous_dna)
        if minus: rseq = Seq(complement(s), IUPAC.ambiguous_dna)


        n = 0
        if plus:
            for frame in self.options['frames']:
                orf = self.translator.translate(seq[frame-1:])
                orfs = string.split(orf.data,'*')
                l = int(self.options['minlength'])
                orfs = filter(lambda x,l=l: len(x) >= l, orfs)
                if self.options['maxlength']:
                    l = int(self.options['maxlength'])
                    orfs = filter(lambda x,l=l: len(x) <= l, orfs)
                for orf in orfs:
                    n = n + 1
                    print self.toFasta('orf_%s:%s:+%d' % (n, self.header,frame), orf)


        if minus:
            for frame in self.options['frames']:
                orf = self.translator.translate(rseq[frame-1:])
                orfs = string.split(orf.data,'*')

                l = int(self.options['minlength'])
                orfs = filter(lambda x,l=l: len(x) >= l, orfs)
                if self.options['maxlength']:
                    l = int(self.options['maxlength'])
                    orfs = filter(lambda x,l=l: len(x) <= l, orfs)
                for orf in orfs:
                    n = n + 1
                    print self.toFasta('orf_%s:%s:-%d' % (n, self.header,frame), orf)
            

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

    shorts = 'h'
    longs = map(lambda x: x +'=', options.keys())
    optlist, args = getopt.getopt(args,shorts, longs)
    if show_help: help()
    for arg in optlist:
        if arg[0] == '-h' or arg[0] == '--help':
            help()
            sys.exit(0)
        for key in options.keys():
            if arg[0][2:] == key:
                options[key] = arg[1]

    file = args[0]

    nextorf = NextORF(file, options)
    nextorf.read_file()

    
