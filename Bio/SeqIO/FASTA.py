import os, string
import Bio.Alphabet

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# This can be made a lot faster by using infile.readlines!
# (Alas, that's a bit complicated to write.)

class FastaReader:
    def __init__(self, infile, alphabet = Bio.Alphabet.generic_alphabet):
        self.infile = infile
        self.alphabet = alphabet

        # find the start of data
        line = infile.readline()
        while line and line[0] != ">":
            line = infile.readline()
        self._lookahead = line
        
        self._n = 0

    def next(self):
        self._n = self._n + 1

        line = self._lookahead
        if not line:
            return None

        # Like bioperl, I assume the first word is the name/id and the
        # rest of the line (after the first whitespace) is a
        # description.  If there's only one word, it's the id.
        x = string.split(line[1:].rstrip(), None, 1)
        if len(x) == 1:
            id = x[0]
            desc = ""
        else:
            id, desc = x
            
        lines = []
        line = self.infile.readline()
        while line:
            if line[0] == ">":
                break
            lines.append(line.rstrip())
            line = self.infile.readline()
            
        self._lookahead = line

        # Unlike bioperl, I assume whitespace is significant.
        return SeqRecord(Seq(string.join(lines, ""), self.alphabet),
                         id = id, name = id, description = desc)
        
    def __getitem__(self, i):
        # wrapper to the normal Python "for spam in list:" idiom
        assert i == self._n  # forward iteration only!
        x = self.next()
        if x is None:
            raise IndexError, i
        return x
        
class FastaWriter:
    def __init__(self, outfile):
        self.outfile = outfile

    def write(self, record):
        id = record.id
        assert os.linesep not in id

        description = record.description
        assert os.linesep not in description
        
        self.outfile.write(">%s %s%s" % (id, description,os.linesep))

        data = record.seq.tostring()
        for i in range(0, len(data), 60):
            self.outfile.write(data[i:i+60] + os.linesep)

    
    def write_records(self, records):
        # In general, can assume homogenous records... useful?
        for record in records:
            self.write(record)

    def close(self):
        return self.outfile.close()

    def flush(self):
        return self.outfile.flush()
