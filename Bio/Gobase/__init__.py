# Copyright 2000 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Gobase

This module provides code to work with files from
http://megasun.bch.umontreal.ca/gobase/


Classes:
Record             Holds gobase sequence data.
Iterator           Iterates over sequence data in a gobase file.
Dictionary         Accesses a gobase file using a dictionary interface.
RecordParser       Parses gobase sequence data into a Record object.

_Scanner           Scans a gobase-format stream.
_RecordConsumer    Consumes gobase data to a Record object.


Functions:
index_file         Index a FASTA file for a Dictionary.

"""
from types import *
import string
import re
from Bio import File
from Bio import Index
from Bio.ParserSupport import *

class Record:
    """Holds information from a Gobase record.

    Members:
    species_name
    taxon_division
    gobase_id
    """
    def __init__(self, colwidth=60):
        """__init__(self, colwidth=60)

        Create a new Record.  colwidth specifies the number of residues
        to put on each line.

        """
        self.species_name = ''
        self.taxon_division = ''

class SequenceRecord( Record ):
    """Holds information from a Gobase record.

    Members:
    molecule_type
    is_plasmid
    shape
    submission_date
    update_date
    entrez_record
    genbank_accession
    """
    def __init__(self, colwidth=60):
        """__init__(self, colwidth=60)

        Create a new Record.  colwidth specifies the number of residues
        to put on each line.

        """
        Record.__init__( self )
        self.molecule_type = ''
        self.is_plasmid = ''
        self.shape = ''
        self.submission_date = ''
        self.update_date = ''
        self.entrez_record = ''
        self.genbank_accession = ''

class GeneRecord( Record ):
    """Holds information from a Gobase record.

    Members:
     """
    def __init__(self, colwidth=60):
        """__init__(self, colwidth=60)

        Create a new Record.  colwidth specifies the number of residues
        to put on each line.

        """
        Record.__init__( self )
        self.gene_class = ''
        self.plasmid_encoded = ''
        self.is_partial_gene = ''
        self.is_pseudo_gene = ''
        self.is_transpliced_gene = ''
        self.chloroplast_origin = ''
        self.contains_intron = ''
        self.orf = ''
        self.included_in_intron = ''
        self.published_info = ''
        self.genbank_accession = ''
        self.entrez_record = ''
        self.product_type = ''
        self.product_class = ''

class ProteinRecord( Record ):
    """Holds information from a Gobase record.

    Members:
    product_class
    gene_class
    is_partial_protein
    is_plasmid
    function
    entry_record
    """
    def __init__(self, colwidth=60):
        """__init__(self, colwidth=60)

        Create a new Record.  colwidth specifies the number of residues
        to put on each line.

        """
        Record.__init__( self )
        self.product_class = ''
        self.gene_class = ''
        self.is_partial_protein = ''
        self.is_plasmid = ''
        self.is_pseudo = ''
        self.function = ''
        self.entry_record = ''

class Iterator:
    """Returns one record at a time from a Gobase file.

    Methods:
    next   Return the next record from the stream, or None.

    """
    def __init__(self, handle, parser=None):
        """__init__(self, handle, parser=None)

        Create a new iterator.  handle is a file-like object.  parser
        is an optional Parser object to change the results into another form.
        If set to None, then the raw contents of the file will be returned.

        """
        if type(handle) is not FileType and type(handle) is not InstanceType:
            raise ValueError, "I expected a file handle or file-like object"
        self._uhandle = SGMLHandle( File.UndoHandle( handle ) )
        self._parser = parser

    def next(self):
        """next(self) -> object

        Return the next gobase record from the file.  If no more records,
        return None.

        """
        lines = []
        first_tag = 'Recognition Sequence'
        while 1:
            line = self._uhandle.readline()
            if not line:
                break
            if line[:len( first_tag )] == 'first_tag':
                self._uhandle.saveline(line)
                break

        if not line:
            return None

        if self._parser is not None:
            return self._parser.parse(File.StringHandle(data))
        return data

class Dictionary:
    """Accesses a gobase file using a dictionary interface.

    """
    __filename_key = '__filename'

    def __init__(self, indexname, parser=None):
        """__init__(self, indexname, parser=None)

        Open a Gobase Dictionary.  indexname is the name of the
        index for the dictionary.  The index should have been created
        using the index_file function.  parser is an optional Parser
        object to change the results into another form.  If set to None,
        then the raw contents of the file will be returned.

        """
        self._index = Index.Index(indexname)
        self._handle = open(self._index[Dictionary.__filename_key])
        self._parser = parser

    def __len__(self):
        return len(self._index)

    def __getitem__(self, key):
        start, len = self._index[key]
        self._handle.seek(start)
        data = self._handle.read(len)
        if self._parser is not None:
            return self._parser.parse(File.StringHandle(data))
        return data

    def __getattr__(self, name):
        return getattr(self._index, name)

class RecordParser:
    """Parses Gobase sequence data into a Record object.

    """
    def __init__(self):
        self._scanner = _Scanner()
        self._consumer = _RecordConsumer()

    def parse(self, handle):
        self._scanner.feed(handle, self._consumer)
        return self._consumer.data

class _Scanner:
    """Scans a gobase file.

    Methods:
    feed   Feed in one gobase record.

    """
    def feed(self, handle, consumer):
        """feed(self, handle, consumer)

        Feed in gobase data for scanning.  handle is a file-like object
        containing gobase data.  consumer is a Consumer object that will
        receive events as the gobase data is scanned.

        """
        if isinstance(handle, File.UndoHandle):
            uhandle = handle
        else:
            uhandle = File.UndoHandle(handle)
        uhandle = File.SGMLHandle( uhandle )

        if uhandle.peekline():
            self._scan_record(uhandle, consumer)

    def _scan_line(self, uhandle ):
        line = safe_readline( uhandle )
        line = string.join( string.split( line ), ' ' ) + ' '
        return line

    def _text_in( self, uhandle, text, count ):
        for j in range( count ):
            try:
                line = self._scan_line( uhandle )
                text = text + line
            except:
                if( line == '' ):
                    return text
        return text

    def _scan_sequence_record( self, text, consumer ):
        data = consumer.data
        next_item = self._scan_field( text, 'Molecule type:', 'Species name:' )
        data.molecule_type = consumer.text_field( next_item )

        next_item = self._scan_field( text, 'Shape of molecule:', 'Sequence length:' )
        data.shape = consumer.text_field( next_item )

        next_item = self._scan_field( text, 'Plasmid:', 'Complete genome:' )
        data.is_plasmid = consumer.text_field( next_item )

        next_item = self._scan_field( text, 'NCBI Entrez record:', 'Genbank accession:' )
        data.entrez_record = consumer.text_field( next_item )

        next_item = self._scan_field( text, 'Genbank accession:', 'Coding gene(s):' )
        data.genbank_accession = consumer.text_field( next_item )
        consumer.data = data

    def _scan_gene_record( self, text, consumer ):
        data = consumer.data
        next_item = self._scan_field( text, 'Gene Class:', 'Species name:' )
        data.gene_class = consumer.text_field( next_item )

        next_item = self._scan_field( text, 'Plasmid encoded:', 'Partial gene:' )
        data.is_plasmid = consumer.word_field( next_item )

        next_item = self._scan_field( text, 'Partial gene:', 'Pseudo:' )
        data.is_partial_gene = consumer.text_field( next_item )

        next_item = self._scan_field( text, 'Pseudo:', 'Transpliced gene:' )
        data.is_pseudo_gene = consumer.text_field( next_item )

        next_item = self._scan_field( text, 'Transpliced gene:', 'Chloroplast origin:' )
        data.is_transpliced_gene = consumer.text_field( next_item )

        next_item = self._scan_field( text, 'Chloroplast origin:', 'Contains intron(s):' )
        data.chloroplast_origin = consumer.word_field( next_item )

        next_item = self._scan_field( text, 'Contains intron(s):' )
        data.contains_intron = consumer.word_field( next_item )

        next_item = self._scan_field( text, 'Included in intron:' )
        data.included_in_intron = consumer.word_field( next_item )

        next_item = self._scan_field( text, 'ORF:' )
        data.orf = consumer.word_field( next_item )

        next_item = self._scan_field( text, 'NCBI Entrez record:' )
        data.entrez_record = consumer.word_field( next_item )

        next_item = self._scan_field( text, 'Genbank accession:', 'Product type:' )
        data.genbank_accession = consumer.word_field( next_item )

        next_item = self._scan_field( text, 'Product type:', 'Product Class:' )
        data.product_type = consumer.text_field( next_item )

        next_item = self._scan_field( text, 'Product Class:' )
        data.product_class = consumer.text_field( next_item )

        consumer.data = data

    def _scan_protein_record( self, text, consumer ):
        data = consumer.data
        next_item = self._scan_field( text, 'Product Class:', 'Species name:' )
        data.product_class = consumer.text_field( next_item )

        next_item = self._scan_field( text, 'Gene Class:', 'Partial protein:' )
        data.gene_class = consumer.text_field( next_item )

        next_item = self._scan_field( text, 'Partial protein:', 'Conflict:' )
        data.is_partial_protein = consumer.text_field( next_item )

        next_item = self._scan_field( text, 'Plasmid:', 'Sequence length:' )
        data.is_plasmid = consumer.text_field( next_item )

        next_item = self._scan_field( text, 'General function:' )
        data.function = consumer.text_field( next_item )

        next_item = self._scan_field( text, 'NCBI Entrez record:' )
        data.entrez_record = consumer.word_field( next_item )

        consumer.data = data

    def _scan_record(self, uhandle, consumer):
        text = ''
        text = self._text_in( uhandle, text, 100 )
        text = string.lstrip( text )

        if( string.find( text, 'Sequence' ) == 0 ):
            consumer.data = SequenceRecord()
            self._scan_sequence_record( text, consumer )
        elif( string.find( text, 'Gene' ) == 0 ):
            consumer.data = GeneRecord()
            self._scan_gene_record( text, consumer )
        elif( string.find( text, 'Protein' ) == 0 ):
            consumer.data = ProteinRecord()
            self._scan_protein_record( text, consumer )
        else:
            print 'UNKNOWN!!!!!!'

        data = consumer.data
        next_item = self._scan_field( text, 'Species name:', 'Taxon division' )
        data.species_name = consumer.text_field( next_item )

        next_item = self._scan_field( text, 'Taxon division:' )
        print next_item
        data.taxon_division = consumer.word_field( next_item )
        consumer.data = data

#        consumer.end_sequence()


    def _scan_field(self, text, field, next_field = None ):
        start = string.find( text, field )
        if( start == -1 ):
            return ''
        if( next_field == None ):
            pattern = re.compile( '[A-Z][a-z0-9 ]+:' )
            offset = start + len( field )
            match = pattern.search( text[ offset: ] )
            if match:
                end = offset + match.start()
            else:
                end = start + 40
        else:
            end = string.find( text, next_field )
            if( end == -1 ):
                return ''
        next_item = text[ start:end ]
        return( next_item )


class _RecordConsumer(AbstractConsumer):
    """Consumer that converts a gobase record to a Record object.

    Members:
    data    Record with gobase data.

    """
    def __init__(self):
        self.data = None

    def end_sequence(self):
        pass

    def text_field( self, line ):
        if( line == '' ):
            return ''
        cols = string.split( line, ': ' )
        return( cols[ 1 ] )

    def int_field( self, line ):
        if( line == '' ):
            return None
        cols = string.split( line, ': ' )
        return( int( cols[ 1 ] ) )

    def word_field( self, line ):
        if( line == '' ):
            return ''
        cols = string.split( line, ': ' )
        cols = string.split( cols[ 1 ] )
        return( cols[ 0 ] )

    def date_field( self, line ):
        if( line == '' ):
            return ''
        cols = string.split( line, ':' )
        cols = string.split( cols[ 1 ] )
        return( string.join( cols[ :3 ] ) )


def index_file(filename, indexname, rec2key=None):
    """index_file(filename, ind/exname, rec2key=None)

    Index a gobase file.  filename is the name of the file.
    indexname is the name of the dictionary.  rec2key is an
    optional callback that takes a Record and generates a unique key
    (e.g. the accession number) for the record.  If not specified,
    the sequence title will be used.

    """
    if not os.path.exists(filename):
        raise ValueError, "%s does not exist" % filename

    index = Index.Index(indexname, truncate=1)
    index[Dictionary._Dictionary__filename_key] = filename

    iter = Iterator(open(filename), parser=RecordParser())
    while 1:
        start = iter._uhandle.tell()
        rec = iter.next()
        length = iter._uhandle.tell() - start

        if rec is None:
            break
        if rec2key is not None:
            key = rec2key(rec)
        else:
            key = rec.title

        if not key:
            raise KeyError, "empty sequence key was produced"
        elif index.has_key(key):
            raise KeyError, "duplicate key %s found" % key

        index[key] = start, length
