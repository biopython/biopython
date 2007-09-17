"""Utilities for working with FASTA-formatted sequences.

This module uses Martel-based parsing to speed up the parsing process.

Classes:
Record             Holds FASTA sequence data.
Iterator           Iterates over sequence data in a FASTA file.
Dictionary         Accesses a FASTA file using a dictionary interface.
RecordParser       Parses FASTA sequence data into a Record object.
SequenceParser     Parses FASTA sequence data into a Sequence object.

Functions:
index_file         Index a FASTA file for a Dictionary.
"""
from Bio import Seq
from Bio import SeqRecord
from Bio import Alphabet

#These imports are only used by the deprecated dictionary functions/classes
import cStringIO
from Bio import Mindy
from Bio.Mindy import SimpleSeqRecord

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
        while True :
            line = handle.readline()
            if line[0] == ">" :
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
            if line[0] == "#" :
                if self._debug : print "Ignoring comment line"
                pass
            else :
                lines.append(line.rstrip())
            line = self.handle.readline()
        self._lookahead = line
        if self._debug : print "Debug: '%s' and '%s'" % (title, "".join(lines))
        if self._parser is None:
            return "\n".join(lines)
        else :
            return self._parser.parse_string("\n".join(lines))

class RecordParser:
    """Parses FASTA sequence data into a Fasta.Record object.
    """
    def __init__(self, debug = 0):
        pass

    def parse_string(self, text) :
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
    
    def parse_string(self, text) :
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

class Dictionary(dict):
    """Accesses an indexed FASTA file using a dictionary interface. DEPRECATED
    """
    def __init__(self, indexname, parser=None, filename = None):
        """Open a Fasta Dictionary.  DEPRECATED
        
        indexname is the name of the index for the dictionary.  The index should 
        have been created using the index_file function.  
        
        parser is an optional Parser object to change the results into another 
        form.  If set to None, then the raw contents of the file will be returned.
        
        filename specifies the name of the file that this index references. 
        This is useful in cases where the file has been moved since indexing.
        If no filename is supplied (the default) the filename stored in the
        index will be used. XXX This is no longer supported -- use symbolic
        links in the filesystem.
        """
        
        import warnings
        warnings.warn("Bio.Fasta.index_file Bio.Fasta.Dictionary are deprecated." \
                      + " We hope an in memory dictionary, for example using the" \
                      + " Bio.SeqIO.to_dict() function, will be suitable for" \
                      + " most users.  Please get in touch on the mailing lists if" \
                      + " this (or its removal) causes any problems for you.",
                      DeprecationWarning)

        # we can't support finding the index file name if we want to follow
        # standard open-bio fetching protocols.
        if filename is not None:
            raise AttributeError("Specifying filenames is no longer supported")

        self._index = Mindy.open(indexname)
        self._parser = parser
        
        primary_key_retriever = self._index['id']
        for k in primary_key_retriever.keys():
            dict.__setitem__(self,k,None)


    def _load_seq(self,key):
        try:
            seqs = self._index.lookup(id = key)
            # if we can't do that, we have to try and fetch by alias
        except KeyError:
            seqs = self._index.lookup(aliases = key)

        if len(seqs) == 1:
            seq = seqs[0]
        else:
            raise KeyError("Multiple sequences found for %s" % key)

        if self._parser:
            handle = cStringIO.StringIO(seq.text)
            self[key] = self._parser.parse(handle)
        else:
            self[key] = seq.text
                                                                

    def __getitem__(self, key):
        if self.has_key(key) and dict.__getitem__(self,key) is None:
            self._load_seq(key)
        return dict.__getitem__(self,key)
            

def index_file(filename, indexname, rec2key = None, use_berkeley = 0):
    """Index a FASTA file. DEPRECATED

    filename is the name of the file to index.

    indexname is the name of the dictionary to be created. This can be
    just the name of the index, in which case the index information will
    be created in a directory of the given index name in the current
    directory, or a full pathname to a directory to save the indexing
    information.

    rec2key is an optional callback fuction that takes a Fasta Record and 
    generates a unique key (e.g. the accession number) for the record.
    Optionally, it can also return 3 items, to be used as the id (unique key)
    name, and aliases for the index. If not specified, the sequence title 
    will be used.

    use_berkeley specifies whether to use the BerkeleyDB indexer, which 
    uses the bsddb3 wrappers around the embedded database Berkeley DB. By
    default, the standard flat file (non-Berkeley) indexes are used.
    """

    import warnings
    warnings.warn("Bio.Fasta.index_file Bio.Fasta.Dictionary are deprecated." \
                  + " We hope an in memory dictionary, for example using the" \
                  + " Bio.SeqIO.to_dict() function, will be suitable for" \
                  + " most users.  Please get in touch on the mailing lists if" \
                  + " this (or its removal) causes any problems for you.",
                  DeprecationWarning)

    if rec2key:
        indexer = _FastaFunctionIndexer(rec2key)
    else:
        indexer = _FastaTitleIndexer()

    if use_berkeley:
        SimpleSeqRecord.create_berkeleydb([filename], indexname, indexer)
    else:
        SimpleSeqRecord.create_flatdb([filename], indexname, indexer)

class _FastaTitleIndexer(SimpleSeqRecord.BaseSeqRecordIndexer):
    """Simple indexer to index by the title of a FASTA record.

    This doesn't do anything fancy, just gets the title and uses that as the
    identifier.
    """
    def __init__(self):
        SimpleSeqRecord.BaseSeqRecordIndexer.__init__(self)

    def primary_key_name(self):
        return "id"

    def secondary_key_names(self):
        return ["name", "aliases"]

    def get_id_dictionary(self, seq_record):
        sequence_id = seq_record.description
            
        id_info = {"id" : [sequence_id],
                   "name" : [],
                   "aliases" : []}
        return id_info


class _FastaFunctionIndexer(SimpleSeqRecord.BaseSeqRecordIndexer):
    """Indexer to index based on values returned by a function.
    
    This class is passed a function to parse description titles from a Fasta
    title. It needs to return either one item, which is an id from the title,
    or three items which are (in order), the id, a list of names, and a list
    of aliases.

    This indexer allows indexing to be completely flexible based on passed
    functions.
    """
    def __init__(self, index_function):
        SimpleSeqRecord.BaseSeqRecordIndexer.__init__(self)
        self.index_function = index_function

    def primary_key_name(self):
        return "id"

    def secondary_key_names(self):
        return ["name", "aliases"]

    def get_id_dictionary(self, seq_record):
        # make a FASTA record to make this compatible with previous Biopython
        # code
        tmp_rec = Record()
        tmp_rec.title = seq_record.description
        tmp_rec.sequence = seq_record.seq.data
        items = self.index_function(tmp_rec)
        if type(items) is not type([]) and type(items) is not type(()):
            items = [items]
        if len(items) == 1:
            seq_id = items[0]
            name = []
            aliases = []
        elif len(items) == 3:
            seq_id, name, aliases = items
        else:
            raise ValueError("Unexpected items from index function: %s" %
                    (items))
        
        return {"id" : [seq_id],
                "name" : name,
                "aliases" : aliases}
