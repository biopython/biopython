# Copyright 2000 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
This module provides code to work with files from Rebase.
http://rebase.neb.com/rebase/rebase.html


Classes:
Record             Holds rebase sequence data.
Iterator           Iterates over sequence data in a rebase file.
Dictionary         Accesses a rebase file using a dictionary interface.
RecordParser       Parses rebase sequence data into a Record object.

_Scanner           Scans a rebase-format stream.
_RecordConsumer    Consumes rebase data to a Record object.


Functions:
index_file         Index a FASTA file for a Dictionary.

"""
from types import *
import string
from Bio import File
from Bio import Index
from Bio.ParserSupport import *

class Record:
    """Holds information from a FASTA record.

    Members:
    seq_5_to_3       The sequence.
    seq_3_to_5
    enzyme_num       The enzyme number
    pos              Position of cleavage
    prototype        Prototype
    source
    microorganism
    temperature      Growth temperature
    misc             Miscellaneous information
    date_entered
    date_modified
    num_Adeno2
    num_Lambda
    num_pBR322
    num_PhiX174
    num_SV40

    """
    def __init__(self, colwidth=60):
        """__init__(self, colwidth=60)

        Create a new Record.  colwidth specifies the number of residues
        to put on each line.

        """
        self.seq_5_to_3 = ''
        self.seq_3_to_5 = ''
        self.methylation = ''
        self.enzyme_num = None
        self.prototype = ''
        self.source = ''
        self.microorganism = ''
        self.temperature = None
        self.misc = ''
        self.date_entered = ''
        self.date_modified = ''
        self._colwidth = colwidth
        self.num_Adeno2 = 0
        self.num_Lambda = 0
        self.num_pBR322 = 0
        self.num_PhiX174 = 0
        self.num_SV40 = 0

class Iterator:
    """Returns one record at a time from a Rebase file.

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

        Return the next rebase record from the file.  If no more records,
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
    
    def __iter__(self):
        return iter(self.next, None)

class Dictionary:
    """Accesses a rebase file using a dictionary interface.

    """
    __filename_key = '__filename'
    
    def __init__(self, indexname, parser=None):
        """__init__(self, indexname, parser=None)

        Open a Fasta Dictionary.  indexname is the name of the
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
    """Parses FASTA sequence data into a Record object.

    """
    def __init__(self):
        self._scanner = _Scanner()
        self._consumer = _RecordConsumer()

    def parse(self, handle):
        self._scanner.feed(handle, self._consumer)
        return self._consumer.data

class _Scanner:
    """Scans a rebase file.

    Methods:
    feed   Feed in one rebase record.

    """
    def feed(self, handle, consumer):
        """feed(self, handle, consumer)

        Feed in rebase data for scanning.  handle is a file-like object
        containing rebase data.  consumer is a Consumer object that will
        receive events as the rebase data is scanned.

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
            line = self._scan_line( uhandle )
            text = text + line
        return text

    def _scan_record(self, uhandle, consumer):
        consumer.start_sequence()
        text = ''
        text = self._text_in( uhandle, text, 100 )
        self._scan_sequence( text, consumer)
        self._scan_methylation( text, consumer)
        self._scan_enzyme_num( text, consumer )
        self._scan_prototype( text, consumer )
        self._scan_source( text, consumer )
        self._scan_microorganism( text, consumer )
        self._scan_temperature( text, consumer)
        self._scan_date_entered( text, consumer)
        self._scan_date_modified( text, consumer)
        self._scan_Adeno2( text, consumer)
        self._scan_Lambda( text, consumer)
        self._scan_pBR322( text, consumer)
        self._scan_PhiX174( text, consumer)
        self._scan_SV40( text, consumer)
#        consumer.end_sequence()


    def _scan_sequence(self, text, consumer ):
        start = string.find( text, 'Recognition Sequence:' )
        end = string.find( text, 'Base (Type of methylation):' )
        if( end == -1 ):
            end = string.find( text, 'REBASE enzyme #:' )
        next_item = text[ start:end ]
        consumer.sequence( next_item )

    def _scan_methylation(self, text, consumer ):
        start = string.find( text, 'Base (Type of methylation):' )
        if( start != -1 ):
            end = string.find( text, 'REBASE enzyme #:' )
            next_item = text[ start:end ]
            consumer.methylation( next_item )

    def _scan_enzyme_num(self, text, consumer ):
        start = string.find( text, 'REBASE enzyme #:' )
        end = string.find( text, 'Prototype:' )
        next_item = text[ start:end ]
        consumer.enzyme_num( next_item )

    def _scan_prototype(self,  text, consumer ):
        start = string.find( text, 'Prototype:' )
        end = string.find( text, 'Source:' )
        next_item = text[ start:end ]
        consumer.prototype( next_item )

    def _scan_source(self, text, consumer ):
        start = string.find( text, 'Source:' )
        end = string.find( text, 'Microorganism:' )
        next_item = text[ start:end ]
        consumer.source( next_item )


    def _scan_microorganism(self, text, consumer ):
        start = string.find( text, 'Microorganism:' )
        end = string.find( text, 'Growth Temperature:' )
        next_item = text[ start:end ]
        consumer.microorganism( next_item )

    def _scan_temperature(self, text, consumer):
        start = string.find( text, 'Growth Temperature:' )
        end = start + 30
        next_item = text[ start:end ]
        consumer.temperature( next_item )


    def _scan_date_entered(self, text, consumer):
        start = string.find( text, 'Entered:' )
        end = start + 30
        next_item = text[ start:end ]
        consumer.data_entered( next_item )

    def _scan_date_modified(self, text, consumer):
        start = string.find( text, 'Modified:' )
        if( start != -1 ):
            end = start + 30
            next_item = text[ start:end ]
            consumer.data_modified( next_item )

    def _scan_Adeno2( self, text, consumer ):
        start = string.find( text, 'Adeno2:' )
        end = string.find( text, 'Lambda:' )
        next_item = text[ start:end ]
        consumer.num_Adeno2( next_item )

    def _scan_Lambda( self, text, consumer ):
        start = string.find( text, 'Lambda:' )
        end = string.find( text, 'pBR322:' )
        next_item = text[ start:end ]
        consumer.num_Lambda( next_item )

    def _scan_pBR322(self, text, consumer ):
        start = string.find( text, 'pBR322:' )
        end = string.find( text, 'PhiX174:' )
        next_item = text[ start:end ]
        consumer.num_pBR322( next_item )

    def _scan_PhiX174(self, text, consumer ):
        start = string.find( text, 'PhiX174:' )
        end = string.find( text, 'SV40:' )
        next_item = text[ start:end ]
        consumer.num_PhiX174( next_item )

    def _scan_SV40(self, text, consumer ):
        start = string.find( text, 'SV40:' )
        end = start + 30
        next_item = text[ start:end ]
        consumer.num_SV40( next_item )


class _RecordConsumer(AbstractConsumer):
    """Consumer that converts a rebase record to a Record object.

    Members:
    data    Record with rebase data.

    """
    def __init__(self):
        self.data = None

    def start_sequence(self):
        self.data = Record()

    def end_sequence(self):
        pass

    def sequence( self, line ):
	cols = string.split( line, ': ' )
        sequence = cols[ 1 ]
        sequence = string.strip( sequence )
        if( string.find( sequence, ' ...' ) != -1 ):
            cols = string.split( sequence, '...' )
            self.data.seq_5_to_3 = cols[ 1 ]
        elif( string.lower( sequence ) != 'unknown' ):
            seq_len = len( sequence ) / 2
            self.data.seq_5_to_3 = string.strip( sequence[ :seq_len ] )
            self.data.seq_3_to_5 = string.strip( sequence[ seq_len: ] )

    def methylation( self, line ):
        cols = string.split( line, ': ' )
        self.data.methylation = cols[ 1 ]

    def enzyme_num( self, line ):
        cols = string.split( line, ': ' )
        self.data.enzyme_num = int( cols[ 1 ] )

    def prototype( self, line ):
        cols = string.split( line, ': ' )
        self.data.prototype = cols[ 1 ]

    def source( self, line ):
        cols = string.split( line, ': ' )
        self.data.source = cols[ 1 ]

    def microorganism( self, line ):
        cols = string.split( line, ': ' )
        self.data.microorganism = cols[ 1 ]

    def temperature( self, line ):
        cols = string.split( line, ':' )
        cols = string.split( cols[ 1 ], ' ' )
        self.data.temperature = cols[ 1 ]

    def data_entered( self, line ):
        cols = string.split( line, ':' )
        cols = string.split( cols[ 1 ] )
        self.data.date_entered = string.join( cols[ :3 ] )

    def data_modified( self, line ):
        cols = string.split( line, ':' )
        cols = string.split( cols[ 1 ] )
        self.data.date_modified = string.join( cols[ :3 ] )

    def num_Adeno2( self, line ):
        cols = string.split( line, ': ' )
        self.data.num_Adeno2 = int( cols[ 1 ] )

    def num_Lambda( self, line ):
        cols = string.split( line, ': ' )
        self.data.num_Lambda = int( cols[ 1 ] )

    def num_pBR322( self, line ):
        cols = string.split( line, ': ' )
        self.data.num_pBR322 = int( cols[ 1 ] )

    def num_PhiX174( self, line ):
        cols = string.split( line, ': ' )
        self.data.num_PhiX174 = int( cols[ 1 ] )

    def num_SV40( self, line ):
        cols = string.split( line, ':' )
        cols = string.split( cols[ 1 ], ' ' )
        self.data.num_SV40 = cols[ 1 ]

def index_file(filename, indexname, rec2key=None):
    """index_file(filename, ind/exname, rec2key=None)

    Index a rebase file.  filename is the name of the file.
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
