#!/usr/bin/env python
# Created: Tue Sep 11 17:21:54 2001
# Last changed: Time-stamp: <01/09/19 13:28:22 thomas>
# thomas@cbs.dtu.dk, http://www.cbs.dtu.dk/thomas/index.html
# File: generic.py
# based on Brads's code


import sys
import os, re, time
sys.path.insert(0, os.path.expanduser('~thomas/cbs/python/biopython'))

import string
import Bio.Alphabet

from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord

class SeqRecord:
    # possible backwards incompatibility !
    # all id and descriptions are stripped - NO MORE '\n'
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

    def get_header(self, line):
        try:
            x = string.split(line[1:-1], None, 1)
            if len(x) == 1:
                id = x[0].strip()
                desc = ""
            else:
                id, desc = [x.strip() for x in x]

        except:
            print >> sys.stderr, 'Unable to get header !!!'
            print >> sys.stderr, 'Offending line:', line
            sys.exit(0)
                
        return (id, desc)
    
    def next(self):
        self._n = self._n + 1

        line = self._lookahead
        if not line: return None

        id, desc = self.get_header(line)
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
        description = record.description
        
        self.outstream.write(">%s %s\n" % (id, description))

        data = record.seq.tostring()
        for i in range(0, len(data), 60):
            self.outstream.write(data[i:i+60] + "\n")

class LargeFastaFormat(GenericFormat):
    def __init__(self, instream=None, outstream=None, alphabet = Bio.Alphabet.generic_alphabet):
        GenericFormat.__init__(self, instream, outstream, alphabet)
        self.entries = None
        
    def next(self):
        if not self.entries:
            txt = instream.read()
            self.entries = txt.split('>')[1:]
            self._n = -1

        self._n += 1
        if self._n >= len(self.entries): return None

        entry = self.entries[self._n]
        
        name,seq= entry.split('\n',1)
        name, desc = self.get_header(name)
        
        seq = seq.replace('\n','')
        return SeqRecord(Seq(seq, self.alphabet), id = name,
                         name = name, description = desc)

    
    def write(self, record):
        id = record.id
        description = record.description
        
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
            
class EMBLFormat(GenericFormat):
    order = ['AC', 'DT', 'DE', 'GN', 'OS', 'OC', 'DR']
    
    def __init__(self, instream=None, outstream=None, alphabet = Bio.Alphabet.generic_alphabet):
        GenericFormat.__init__(self, instream, outstream, alphabet, 'ID   ')
        if instream: self.find_start()

        
        
    def next(self):
        self._n = self._n + 1

        line = self._lookahead
        if not line: return None
        
        dict = {}
        while line:
            if line[:2] == '//': break
            feature = line[:2]
            if feature == '  ': feature = 'SQ'
            dict.setdefault(feature, [])
            dict[feature].append(line[5:].strip())

            line = self.instream.readline()
        assert 'ID' in dict.keys()
        
        self._lookahead = self.instream.readline()

        seq = Seq(string.join(dict['SQ'][1:], ''), self.alphabet)
        ID = dict['ID'][0].split()[0]

        rec = SeqRecord(seq, id = ID, name = ID, description = dict.get('DE',[])[0])
        rec.annotations = dict

        return rec
        
    def write(self, record):
        id = record.id

        description = record.description
        if description and not description[-1] == '\n':
            description = description + '\n'

        dataclass = 'STANDARD;'
        division = 'PRT;' # fix that to change for DNA sequence
        length = len(record.seq)

        dict = record.annotations
        put = self.outstream.write

        if dict.has_key('ID'):
            put('ID   %s' % dict['ID'][0])
        else:
            put('ID   %-12s%+12s%+10s% 6d AA.\n' % (id, dataclass, division, length))
            
        features = record.annotations.keys()
        if 'ID' in features: features.remove('ID')
        if 'SQ' in features: features.remove('SQ')

        for feature in self.order:
            if not feature in features: continue
            features.remove(feature)
            for line in dict[feature]:
                put('%s   %s\n' % (feature, line))

        for feature in features:
            if feature[0] == 'R': continue
            # TODO
            # fix the order of all R* features
            for line in dict[feature]:
                put('%s   %s\n' % (feature, line))

        if dict.has_key('SQ'):
            put('SQ   %s\n' % '\n     '.join(dict['SQ'][1:]))
        else:
            put('SQ   SEQUENCE%4d AA;\n' % length)
            data = record.seq.tostring()
            for i in range(0, len(data), 60):
                put(data[i:i+60] + "\n")
            
        put('//\n')
            
class GCGFormat(GenericFormat):
    def __init__(self, instream=None, outstream=None, alphabet = Bio.Alphabet.generic_alphabet):
        GenericFormat.__init__(self, instream, outstream, alphabet)

    def next(self):
        line = instream.readline()
        if not line: return None

        desc = ''
        while 1:
            if line.find('..') > -1: break
            desc += line.strip()
            line = instream.readline()

        id = line.split()[0]

        seq = ''
        while 1:
            line = instream.readline()
            if not line: break
            seq += re.sub('[^a-zA-Z-]','',line).upper()

        return SeqRecord(Seq(seq, self.alphabet),
                         id = id, name = id, description = desc)

        
    def write(self, record):
        id = record.id
        description = record.description

        put = self.outstream.write

        if not description: description = id
        put(description)
        if description[-1] != '\n': put('\n')

        timestamp = time.strftime('%B %d, %Y %H:%M', time.localtime(time.time()))
        put('%s Length: %d %s Type: P\n' % (id, len(record.seq), timestamp))
        data = record.seq.tostring()
        for i in range(0, len(data), 60):
            put('% 6d %s\n' % (i+1,data[i:i+60]))

        put('\n')
        
class ClustalFormat(GenericFormat):
    def __init__(self, instream=None, outstream=None, alphabet = Bio.Alphabet.generic_alphabet):
        GenericFormat.__init__(self, instream, outstream, alphabet)
        self.sequences = None
        self.ids = None

    def ParseAlignment(self):
        line = self.instream.readline()
        if not line: return
        if not line[:7] == 'CLUSTAL':
            print >> sys.stderr, 'Warning file does not start with CLUSTAL header'

        dict = {}
        self.ids = []
        while 1:
            line = self.instream.readline()
            if not line: break
            if line[0] == ' ': continue
            fields = line.split()
            if not len(fields): continue
            name, seq = fields
            if not name in self.ids: self.ids.append(name)
            dict.setdefault(name, '')
            dict[name] += seq.upper()

        self.sequences = dict
        self._n = -1
            
        
    def next(self):
        if not self.ids: self.ParseAlignment()

        self._n += 1
        if self._n >= len(self.ids): return None

        name = self.ids[self._n]
        seq = self.sequences[name]
        
        return SeqRecord(Seq(seq, self.alphabet),id = name, name = name,
                         description = 'Clustal Alignment')

class NexusFormat(GenericFormat):
    def __init__(self, instream=None, outstream=None, alphabet = Bio.Alphabet.generic_alphabet):
        GenericFormat.__init__(self, instream, outstream, alphabet)
        self.sequences = None
        self.ids = None

    def ParseNexus(self):
        line = self.instream.readline()
        if not line: return

        self.ids = []

        found_data, found_matrix = 0,0
        while 1:
            # search for the data block
            if line.lower().find('begin data;') > -1:
                found_data = 1
                break
            line = self.instream.readline()
            if not line: break

        while 1:
            # search for the matrix block
            if line.lower().find('matrix') > -1:
                found_matrix = 1
                break
            line = self.instream.readline()
            if not line: break


        dict = {}
        while 1:
            # read name, sequence pairs until first ';'
            line = self.instream.readline()
            if not line: break
            if line.find(';') > -1: break
            
            fields = line.split()
            if len(fields) <2: continue

            name = fields[0]
            if not name in self.ids: self.ids.append(name)
            dict.setdefault(name, '')
            dict[name] += ''.join(fields[1:])

        self.sequences = dict
        self._n = -1
            
        
    def next(self):
        if not self.ids: self.ParseNexus()

        self._n += 1
        if self._n >= len(self.ids): return None

        name = self.ids[self._n]
        seq = self.sequences[name]
        
        return SeqRecord(Seq(seq, self.alphabet),id = name, name = name,
                         description = '')

class ReadSeq:
    def __init__(self):
        self.fdict = {
            'fasta': FastaFormat,
            'largefasta': LargeFastaFormat,
            'embl': EMBLFormat,
            'pir': PirFormat,
            'gcg': GCGFormat,
            'clustal': ClustalFormat, # read only
            'nexus': NexusFormat,     # read only
            }
        

    def Convert(self, informat, outformat, instream=None, outstream=None):
        instream = instream or sys.stdin
        outstrem = outstream or sys.stdout

        if instream == '-': instream = sys.stdin
        if outstream == '-': outstream = sys.stdout
        if type(instream) != type(sys.stdin): instream = open(instream)
        if type(outstream) != type(sys.stdout): outstream = open(outstream, 'w+')

        try:
            reader = self.fdict[informat.lower()](instream=instream)
            writer = self.fdict[outformat.lower()](outstream = outstream)
        except:
            print >> sys.stderr, 'Unknown format: %s -> %s' % (informat, outformat)
            return

        while 1:
            rec = reader.next()
            if not rec: break
            writer.write(rec)


if __name__ == '__main__':

    readseq = ReadSeq()
    
    try:
        instream = sys.argv[1]
        informat = sys.argv[2]
        outstream = sys.argv[3]
        outformat = sys.argv[4]
    except IndexError: 
        p = os.path.basename(sys.argv[0])
        print >> sys.stderr, 'Usage: %s <instream> <informat> <outstream> <outformat>' % p
        print >> sys.stderr, '\twhere "-" can be used for stdin resp. stdout'
        print >> sys.stderr, '\tKnown formats: %s' % ', '.join(readseq.fdict.keys())

        print >> sys.stderr, '\n\te.g. %s eftu.fas fasta eftu.emb embl' % p 
        print >> sys.stderr, '\tor   zcat test.aln.gz | %s - clustal - fasta' % p
        sys.exit(0)

    readseq.Convert(informat, outformat, instream, outstream)
    
        
        
    
    
