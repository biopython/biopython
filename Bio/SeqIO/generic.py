#!/usr/bin/env python
# Created: Tue Sep 11 17:21:54 2001
# Last changed: Time-stamp: <01/09/18 11:00:46 thomas>
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
    
    def next(self):
        self._n = self._n + 1

        line = self._lookahead
        if not line: return None

        x = string.split(line[1:-1], None, 1)
        if len(x) == 1:
            id = x.strip()
            desc = ""
        else:
            id, desc = [x.strip() for x in x]

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
        if not description[-1] == '\n': description = description + '\n'

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
    txt2 = """
ID   CAB59873    PRELIMINARY;      PRT;   188 AA.
AC   CAB59873;
DT   04-JUL-2000 (EMBLrel. 62, Created)
DT   04-JUL-2000 (EMBLrel. 62, Last sequence update)
DT   04-JUL-2000 (EMBLrel. 62, Last annotation update)
DE   60S ribosomal protein L11.
GN   P1421.04.
OS   Leishmania major.
OC   Eukaryota; Euglenozoa; Kinetoplastida; Trypanosomatidae; Leishmania.
OX   NCBI_TaxID=5664;
RN   [1]
RP   SEQUENCE FROM N.A.
RC   STRAIN=Friedlin;
RA   Ivens A.C., Lawson D., Murphy L., Quail M., Rajandream M.A.,
RA   Barrell B.G.;
RL   Submitted (DEC-1999) to the EMBL/GenBank/DDBJ databases.
RN   [2]
RP   SEQUENCE FROM N.A.
RC   STRAIN=Friedlin;
RA   Ivens A.C., Lewis S.M., Bagherzadeh A., Zhang L., Chan H.M.,
RA   Smith D.F.;
RT   \"A physical map of the Leishmania major Friedlin genome.\";
RL   Genome Res. 8:135-145(1998).
DR   EMBL; AL132764; CAB59873.1; -.
SQ   SEQUENCE   188 AA;  21645 MW;  9E70E090C1D0FA5C CRC64;
     MVAESKAANP MREIVVKKLC INICVGESGD RLTRASKVLE QLCEQTPVLS RARLTVRTFG
     IRRNEKIAVH CTVRGKKAEE LLEKGLKVKE FELKSYNFAD TGSFGFGIDE HIDLGIKYDP
     STGIYGMDFY VVLGRRGERV AHRKRKCSRV GHSHHVTKEE AMKWFEKVHD GIIFQAKKKK
     KMIRRRRR
//
"""    
    from StringIO import StringIO
    test = FastaFormat(instream = StringIO(txt))
    test2 = PirFormat(outstream = sys.stdout)
    test3 = EMBLFormat(instream = StringIO(txt2))
    test4 = EMBLFormat(outstream = sys.stdout)
    

    while 1:
        r = test.next()
        r2 = test3.next()
        if r:
            test2.write(r)
            test4.write(r)
        if r2:
            test2.write(r2)
            test4.write(r2)

        if not (r or r2): break
        
        
