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
import os
import cStringIO

from Bio import Seq
from Bio import SeqRecord
from Bio import Alphabet

from Bio import Mindy
from Bio.Mindy import SimpleSeqRecord
from Bio.expressions import fasta

from Martel.LAX import LAX

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
        return os.linesep.join(s)

class Iterator:
    """Returns one record at a time from a FASTA file.
    """
    def __init__(self, handle, parser = None, debug = 0):
        """Initialize a new iterator.
        """
        it_builder = fasta.format.make_iterator("record", debug_level = debug)
        self._parser = parser
        self._iterator = it_builder.iterateFile(handle,
                LAX(fields = ['bioformat:sequence', 'bioformat:description']))

    def __iter__(self):
        return iter(self.next, None)

    def next(self):
        try:
            result = self._iterator.next()
        except StopIteration:
            return None

        if self._parser is None:
            # return a string formatted result comparable to the original
            parser = RecordParser()
            rec = parser.convert_lax(result)
            return str(rec) + "\n"
        else:
            return self._parser.convert_lax(result)

class _MartelBaseFastaParser:
    """Base class to provide a old-Biopython style to Martel parsing.
    """
    def __init__(self, debug):
        self._it_builder = fasta.format.make_iterator("record",
                debug_level = debug)

    def parse(self, handle):
        """Parse a single FASTA record in an file handle.

        Internally, this parses the record into a dictionary
        and then passes that on to the convert_lax function defined
        in a derived class.
        """
        iterator = self._it_builder.iterateFile(handle,
                LAX(fields = ['bioformat:sequence', 'bioformat:description']))
        return self.convert_lax(iterator.next())

    def parse_str(self, str):
        return self.parse(cStringIO.StringIO(str))

    def convert_lax(self, result):
        raise NotImplementedError("Derived class must implement")

class RecordParser(_MartelBaseFastaParser):
    """Parses FASTA sequence data into a Record object.
    """
    def __init__(self, debug = 0):
        _MartelBaseFastaParser.__init__(self, debug)

    def convert_lax(self, result):
        """Convert a dictionary LAX parsing result into a Record object.
        """
        rec = Record()
        rec.title = "".join(result['bioformat:description'])
        rec.sequence = "".join(result['bioformat:sequence'])
        return rec

class SequenceParser(_MartelBaseFastaParser):
    """Parses FASTA sequence data into a Sequence object.
    """
    def __init__(self, alphabet = Alphabet.generic_alphabet, title2ids = None,
            debug = 0):
        """Initialize a Scanner and Sequence Consumer.

        Arguments:
        o alphabet - The alphabet of the sequences to be parsed. If not
        passed, this will be set as generic_alphabet.
        o title2ids - A function that, when given the title of the FASTA
        file (without the beginning >), will return the name, id and
        description (in that order) for the record. If this is not given,
        then the entire title line will be used as the description.
        """
        _MartelBaseFastaParser.__init__(self, debug)
        self.alphabet = alphabet
        self.title2ids = title2ids
    
    def convert_lax(self, result):
        """Convert a dictionary LAX parsing result into a Sequence object.
        """
        seq = Seq.Seq("".join(result['bioformat:sequence']), self.alphabet)
        rec = SeqRecord.SeqRecord(seq)
        
        title = "".join(result['bioformat:description'])
        if self.title2ids:
            seq_id, name, descr = self.title2ids(title)
            rec.id = seq_id
            rec.name = name
            rec.description = descr
        else:
            rec.description = title

        return rec

class Dictionary(dict):
    """Accesses an indexed FASTA file using a dictionary interface.
    """
    def __init__(self, indexname, parser=None, filename = None):
        """Open a Fasta Dictionary.  
        
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
    """Index a FASTA file.

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
