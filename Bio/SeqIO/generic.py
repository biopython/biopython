#!/usr/bin/env python
# Created: Tue Sep 11 17:21:54 2001
# Last changed: Time-stamp: <01/09/19 10:47:20 thomas>
# thomas@cbs.dtu.dk, http://www.cbs.dtu.dk/thomas/index.html
# File: generic.py

import sys
import os, time
sys.path.insert(0, os.path.expanduser('~thomas/cbs/python/biopython'))

import string, re
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
        self.type = 'P' 
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
            id = x[0].strip()
            desc = ""
        else:
            id, desc = [x.strip() for x in x]

        lines = []
        line = self.instream.readline()
        l = len(self.start_indicator)
        print >> sys.stderr, id, desc
        while line:
            if line[:l] == self.start_indicator:
                break
            lines.append(line[:-1])
            line = self.instream.readline()
            
        self._lookahead = line
        print >> sys.stderr, 'record'
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
        if id == description: description = ''
        
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
            txt = self.instream.read()
            self.entries = txt.split('>')[1:]
            self._n = -1
            
        self._n = self._n + 1
        if self._n >= len(self.entries): return None
        
        entry = self.entries[self._n]
        line ,seq= entry.split('\n',1)
        seq = seq.replace('\n','')

        x = string.split(line[1:-1], None, 1)
        if len(x) == 1:
            id = x[0].strip()
            desc = ""
        else:
            id, desc = [x.strip() for x in x]

        return SeqRecord(Seq(seq, self.alphabet), id = id, name = id, description = desc)

    def write(self, record):
        id = record.id
        description = record.description
        if id == description: description = ''
        
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
        description = record.description
        
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

        seq = ''
        for line in dict['SQ'][1:]:
            seq += re.sub('[^a-zA-Z]', '', line).upper()

        seq = Seq(seq, self.alphabet)
        ID = dict['ID'][0].split()[0]

        rec = SeqRecord(seq, id = ID, name = ID, description = dict.get('DE',[])[0])
        rec.annotations = dict

        return rec
        
    def write(self, record):
        id = record.id

        description = record.description
        if description and not description[-1] == '\n': description = description + '\n'

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
        # GCG format allows only one sequence per file, so no need for find_start
        
    def next(self):
        self._n = self._n + 1
        line = self.instream.readline()
        if not line: return None

        desc = ''
        while line:
            # Get the descriptive info (anything before the line with '..')
            if line.find('..') > -1: break
            desc += line.strip()
            line = self.instream.readline()

        # Get ID and Type from the line containing '..'
        ID = line.split()[0]
        m = re.match('Type:\s(\w)\s', line)
        if m:
            Type = m.group(1)
        else:
            Type = 'P'

        # anything after '..' is sequence ...
        seq = ''
        line = self.instream.readline()
        while line:
            if line.find('..') > -1:
                # Arrgggghhh multiple GCG ???
                break
            
            # remove anything that is not alphabet char and make sequence uppercase
            seq += re.sub('[^a-zA-Z]', '', line).upper()
            line = self.instream.readline()
            
        
        seq = Seq(seq, self.alphabet)
        rec = SeqRecord(seq, id = ID, name = ID, description = desc)
        rec.type = Type
        return rec
        

    def write(self, record):
        id = record.id
        description = record.description

        put = self.outstream.write
        if description:
            put(description)
            if description[-1] != '\n': put('\n')
        else: put('%s\n' % id)
        
        timestamp = time.strftime('%B %d, %Y %H:%M', time.localtime(time.time()))
        put('%s Length: %d %s Type: %s ..\n' % (id, len(record.seq), timestamp, record.type))
        
        data = record.seq.tostring()
        for i in range(0, len(data), 60):
            put('% 8d %s\n' % (i+1, data[i:i+60]))
        put('\n')


class ClustalFormat(GenericFormat):
    def __init__(self, instream=None, outstream=None, alphabet = Bio.Alphabet.generic_alphabet):
        GenericFormat.__init__(self, instream, outstream, alphabet)
        self.sequences = None
        self.ids = None
        
    def ParseAlignment(self):
        line = self.instream.readline()
        if not line: return None

        if not line[:7] == 'CLUSTAL':
            sys.stderr.write('Warning, file does not start with a CLUSTAL header\n')
        else:
            line = self.instream.readline()
            
        self.sequences = {}
        dict = self.sequences
        self.ids = []
        while line:
            if not line: break

            if line[0] == ' ':
                line = self.instream.readline()
                continue
            
            fields = line.split()
            if not len(fields):
                line = self.instream.readline()
                continue
            
            name, seq = fields
            if not name in self.ids: self.ids.append(name)
            dict.setdefault(name,'')
            dict[name] += seq
            
            line = self.instream.readline()
            
        self._n = -1
        
    def next(self):
        if not self.ids: self.ParseAlignment()
        self._n = self._n + 1
        if self._n >= len(self.ids): return None
        
        name = self.ids[self._n]
        seq = self.sequences[name] 
        rec = SeqRecord(Seq(seq, self.alphabet), id = name, name = name, description = 'CLustal alignment')
        return rec
    
    def write(self, record):
        pass
    

    
class ReadSeq:
    def __init__(self):
        self.fdict = {'fasta': FastaFormat,
                      'largefasta': LargeFastaFormat,
                      'pir': PirFormat,
                      'embl': EMBLFormat,
                      'gcg': GCGFormat,
                      'clustal': ClustalFormat, # read only
                      }


    def Convert(self, from_format, to_format,
                instream = sys.stdin, outstream = sys.stdout):
        
        if instream == '-': instream = sys.stdin
        if not type(instream) == type(sys.stdin): instream = open(instream)
        if outstream == '-': outstream = sys.stdout
        if not type(outstream) == type(sys.stdout): outstream = open(outstream, 'w+')
        
        try:
            reader = self.fdict[from_format.lower()](instream)
            writer = self.fdict[to_format.lower()](outstream=outstream)
        except:
            import traceback
            traceback.print_exc()
            print >> sys.stderr, 'Unknown Format: %s->%s' % (from_format, to_format)
            return

        while 1:
            record = reader.next()
            if not record: break
            writer.write(record)
        
        

if __name__ == '__main__':
    readseq = ReadSeq()
    
    try:
        infile = sys.argv[1]
        informat = sys.argv[2]
        outfile = sys.argv[3]
        outformat = sys.argv[4]
    except:
        print >> sys.stderr, 'Usage: generic.py <infile> <informat> <outfile> <outformat>'
        print >> sys.stderr, '\tWhere "-" is stdin resp. stdout'
        print >> sys.stderr, '\tKnown formats: %s' % ', '.join(readseq.fdict.keys())
        print >> sys.stderr, '\n\te.g. generic.py myfile.fas fasta myfile.emb embl'
        print >> sys.stderr, '\tor   zcat eftu.aln.gz| generic.py - clustal - fasta'
        
        sys.exit(1)
        
        

    readseq.Convert(informat, outformat, infile, outfile)
    
