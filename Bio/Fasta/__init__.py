"""Utilities for working with FASTA-formatted sequences (DEPRECATED).

Classes:
Record             Holds FASTA sequence data.
Iterator           Iterates over sequence data in a FASTA file.
RecordParser       Parses FASTA sequence data into a Record object.
SequenceParser     Parses FASTA sequence data into a SeqRecord object.

For a long time this module was the most commonly used and best documented
FASTA parser in Biopython.  However, we now recommend using Bio.SeqIO instead.
After being declared obsolete, Bio.Fasta has now been officially deprecated
(with a warning message when imported) and will be removed in a future
release.

If you are already using Bio.Fasta with the SequenceParser to get SeqRecord
objects, then you should be able to switch to the more recent Bio.SeqIO module
very easily as that too uses SeqRecord objects.  For example,

from Bio import Fasta
handle = open("example.fas")
for seq_record in Fasta.Iterator(handle, Fasta.SequenceParser()):
    print seq_record.description
    print seq_record.seq
handle.close()

Using Bio.SeqIO instead this becomes:

from Bio import SeqIO
handle = open("example.fas")
for seq_record in SeqIO.parse(handle, "fasta"):
    print seq_record.description
    print seq_record.seq
handle.close()

Converting an existing code which uses the RecordParser is a little more
complicated as the Bio.Fasta.Record object differs from the SeqRecord.

from Bio import Fasta
handle = open("example.fas")
for record in Fasta.Iterator(handle, Fasta.RecordParser()):
    #record is a Bio.Fasta.Record object
    print record.title #The full title line as a string
    print record.sequence #The sequence as a string
handle.close()

Using Bio.SeqIO instead this becomes:

from Bio import SeqIO
handle = open("example.fas")
for seq_record in SeqIO.parse(handle, "fasta"):
    print seq_record.description #The full title line as a string
    print str(seq_record.seq) #The sequence as a string
handle.close()

Very old code may have used Bio.Fasta.index_file and Dictionary, which were
deprecated in Biopython 1.44 and removed in Biopython 1.46. These allowed
indexing of a FASTA file and access to the records with a dictionary like
interface. Currently using Bio.SeqIO.to_dict to create an in memory dictionary
of SeqRecord objects is the best replacement, but for very large files
additional indexing support for Bio.SeqIO is being considered.
"""
from Bio import Seq
from Bio import SeqRecord
from Bio import Alphabet

import warnings
warnings.warn('Bio.Fasta is deprecated. Please use the "fasta" support in '
              'Bio.SeqIO (or Bio.AlignIO) instead.', DeprecationWarning)

class Record:
    """Holds information from a FASTA record.

    Members:
    title       Title line ('>' character not included).
    sequence    The sequence.
    
    """
    def __init__(self, colwidth=60):
        """__init__(self, colwidth=60)

        Create a new Record.  colwidth specifies the number of residues
        to put on each line when generating FASTA format.

        """
        self.title = ''
        self.sequence = ''
        self._colwidth = colwidth
        
    def __str__(self):
        s = []
        s.append('>%s' % self.title)
        i = 0
        while i < len(self.sequence):
            s.append(self.sequence[i:i+self._colwidth])
            i = i + self._colwidth
        #Was having a problem getting the tests to pass on windows...
        #return os.linesep.join(s)
        return "\n".join(s)

class Iterator:
    """Returns one record at a time from a FASTA file.
    """
    def __init__(self, handle, parser = None, debug = 0):
        """Initialize a new iterator.
        """
        self.handle = handle
        self._parser = parser
        self._debug = debug

        #Skip any text before the first record (e.g. blank lines)
        while True:
            line = handle.readline()
            if not line or line[0] == ">":
                break
            if debug : print "Skipping: " + line
        self._lookahead = line

    def __iter__(self):
        return iter(self.next, None)

    def next(self):
        """Return the next record in the file"""
        line = self._lookahead
        if not line:
            return None
        assert line[0]==">", line
        lines = [line.rstrip()]
        line = self.handle.readline()
        while line:
            if line[0] == ">": break
            if line[0] == "#":
                if self._debug : print "Ignoring comment line"
                pass
            else:
                lines.append(line.rstrip())
            line = self.handle.readline()
        self._lookahead = line
        if self._debug : print "Debug: '%s' and '%s'" % (title, "".join(lines))
        if self._parser is None:
            return "\n".join(lines)
        else:
            return self._parser.parse_string("\n".join(lines))

class RecordParser:
    """Parses FASTA sequence data into a Fasta.Record object.
    """
    def __init__(self, debug = 0):
        pass

    def parse_string(self, text):
        text = text.replace("\r\n","\n") #Crude way of dealing with \r\n
        assert text[0] == ">", text
        text = text.split("\n>",1)[0] # Only do the first record if more than one
        title, sequence = text.split("\n", 1)
        title = title[1:]
        rec = Record()
        rec.title = title
        rec.sequence = sequence.replace("\n","")
        return rec
    
    def parse(self, handle):
        return self.parse_string(handle.read())

class SequenceParser:
    """Parses FASTA sequence data into a SeqRecord object.
    """
    def __init__(self, alphabet = Alphabet.generic_alphabet, title2ids = None,
            debug = 0):
        """Initialize a Scanner and Sequence Consumer.

        Arguments:
        o alphabet - The alphabet of the sequences to be parsed. If not
        passed, this will be set as generic_alphabet.
        o title2ids - A function that, when given the title of the FASTA
        file (without the beginning >), will return the id, name and
        description (in that order) for the record. If this is not given,
        then the entire title line will be used as the description.
        """
        self.alphabet = alphabet
        self.title2ids = title2ids
    
    def parse_string(self, text):
        text = text.replace("\r\n","\n") #Crude way of dealing with \r\n
        assert text[0] == ">", text
        text = text.split("\n>",1)[0] # Only do the first record if more than one
        title, sequence = text.split("\n", 1)
        title = title[1:]

        seq = Seq.Seq(sequence.replace("\n",""), self.alphabet)
        rec = SeqRecord.SeqRecord(seq)
        
        if self.title2ids:
            seq_id, name, descr = self.title2ids(title)
            rec.id = seq_id
            rec.name = name
            rec.description = descr
        else:
            rec.description = title

        return rec

    def parse(self, handle):
        return self.parse_string(handle.read())
