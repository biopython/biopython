#!/usr/bin/env python
# Created: Tue Sep 11 17:21:54 2001
# Last changed: Time-stamp: <01/09/17 09:49:40 thomas>
# thomas@cbs.dtu.dk, http://www.cbs.dtu.dk/thomas/index.html
# File: generic.py

import sys
import os
sys.path.insert(0, os.path.expanduser('~thomas/cbs/python/biopython'))

import string
import Bio.Alphabet

from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord

class SeqRecord:
    def __init__(self, seq, id = "<unknown id>", name = "<unknown name>",
                 description = "<unknown description>"):
        self.seq = seq
        self.id = id
        self.name = name
        self.description = description
        # annotations about the whole sequence
        self.annotations = {}
        
        # annotations about parts of the sequence
        self.features = []
        
    def __str__(self):
        res = ''
        res += '%s %s' % (self.name, self.seq.data)
        return res
    

class GenericFormat:
    def __init__(self, instream=None, outstream=None,
                 alphabet = Bio.Alphabet.generic_alphabet,
                 start_indicator = None):
        self.instream = instream
        self.outstream = outstream
        self.alphabet = alphabet
        self._n = -1
        self._lookahead = None
        self.start_indicator = start_indicator
        
    def find_start(self):
        # find the start of data
        line = self.instream.readline()
        l = len(self.start_indicator)
        while line and line[:l] != self.start_indicator:
            line = self.instream.readline()
        self._lookahead = line
        self._n = 0
    
    def next(self):
        self._n = self._n + 1

        line = self._lookahead
        if not line: return None

        x = string.split(line[1:-1], None, 1)
        if len(x) == 1:
            id = x
            desc = ""
        else:
            id, desc = x
            
        lines = []
        line = self.instream.readline()
        l = len(self.start_indicator)
        while line:
            if line[:l] == self.start_indicator:
                break
            lines.append(line[:-1])
            line = self.instream.readline()
            
        self._lookahead = line

        return SeqRecord(Seq(string.join(lines, ""), self.alphabet),
                         id = id, name = id, description = desc)
        
    def __getitem__(self, i):
        # wrapper to the normal Python "for spam in list:" idiom
        assert i == self._n  # forward iteration only!
        x = self.next()
        if x is None:
            raise IndexError, i
        return x

    def write(self, record):
        pass
    
    def write_records(self, records):
        # In general, can assume homogenous records... useful?
        for record in records:
            self.write(record)

    def close(self):
        return self.outstream.close()

    def flush(self):
        return self.outstream.flush()

class FastaFormat(GenericFormat):
    def __init__(self, instream=None, outstream=None, alphabet = Bio.Alphabet.generic_alphabet):
        GenericFormat.__init__(self, instream, outstream, alphabet, '>')
        if instream: self.find_start()
        
    def write(self, record):
        id = record.id
        assert "\n" not in id

        description = record.description
        assert "\n" not in description
        
        self.outstream.write(">%s %s\n" % (id, description))

        data = record.seq.tostring()
        for i in range(0, len(data), 60):
            self.outstream.write(data[i:i+60] + "\n")
            
class PirFormat(GenericFormat):
    def __init__(self, instream=None, outstream=None, alphabet = Bio.Alphabet.generic_alphabet):
        GenericFormat.__init__(self, instream, outstream, alphabet, '>P1;')
        if instream: self.find_start()

    def write(self, record):
        id = record.id
        assert "\n" not in id

        description = record.description
        assert "\n" not in description
        
        self.outstream.write(">P1;%s %s\n" % (id, description))

        data = record.seq.tostring()
        for i in range(0, len(data), 60):
            self.outstream.write(data[i:i+60] + "\n")

        if data[-1] != '*':
            self.outstream.write("*\n")
            


if __name__ == '__main__':
    txt = """
>TM0001 hypothetical protein
MVYGKEGYGRSKNILLSECVCGIISLELNGFQYFLRGMETL
>TM0002 hypothetical protein
MSPEDWKRLICFHTSKEVLKQTLDDAQQNISDSVSIPLRKY
>TM0003 hypothetical protein
METVKAYEVEDIPAIGFNNSLEVWKLFPASSSRSTSSSFQ
>TM0004 hypothetical protein
MKDLYERFNNSLEVWKLVELFGTSIRIHLFQ
"""
    from StringIO import StringIO
    test = FastaFormat(instream = StringIO(txt))
    test2 = PirFormat(outstream = sys.stdout)
    while 1:
        r = test.next()
        if not r: break
        test2.write(r)
        
