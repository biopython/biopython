
# Copyright 2000 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Rebase.py

This module provides code to work with files from
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
import urllib
import sgmllib
from Bio.ParserSupport import *
import Bio.File

class Record:
    """Holds information from a Rebase record.

    Members:
    seq_5_to_3       The sequence.
    seq_3_to_5

    """
    def __init__(self ):
        """__init__(self )

        Create a new Record.

        """
        self.seq_5_to_3 = ''
        self.seq_3_to_5 = ''

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
        first_tag = '<!DOCTYPE HTML>'
        while 1:
            line = self._uhandle.readline()
            if not line:
                break
            if line[:len( first_tag )] == first_tag:
                self._uhandle.saveline(line)
                break

        if not line:
            return None

        if self._parser is not None:
            return self._parser.parse(File.StringHandle(data))
        return data

class Dictionary:
    """Accesses a rebase file using a dictionary interface.

    """
    __filename_key = '__filename'

    def __init__(self, indexname, parser=None):
        """__init__(self, indexname, parser=None)

        Open a Rebase Dictionary.  indexname is the name of the
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

class EnzymeDict( dict ):
    def print_item( self, item, level = 1 ):
        indent = '    '
        for j in range( 0, level ):
            indent = indent + '    '
        out = ''
        if( type( item ) == type( '' ) ):
            if( item != '' ):
                out = out + '%s%s\n' % ( indent, item )
        elif( type( item ) == type([])):
            for subitem in item:
                out = out + self.print_item( subitem, level + 1 )
        elif( isinstance( item, dict ) ):
            keys = item.keys()
            keys.sort()
            for subitem in keys:
                out = out + '%s%s\n' % ( indent, subitem )
                out = out + self.print_item( item[ subitem ], level + 1 )
        elif( type( item ) == type( {} ) ):
            keys = item.keys()
            keys.sort()
            for subitem in keys:
                out = out +  '%s%s\n' % ( indent, subitem )
                out = out + self.print_item( item[ subitem ], level + 1 )
        else:
            out = out + item
            out = out + '\n'
        return out

    def __str__( self ):
        out = ''
        keys = self.keys()
        keys.sort()
        for key in keys:
            out = out +  '%s\n' % key
            out = out + self.print_item( self[ key ] )
        return out

class RebaseParser(  sgmllib.SGMLParser ):
    """Parses Rebase sequence data into a Record object.

    """
    def reset(self):
        sgmllib.SGMLParser.reset( self )
        self.text = ''
        self.enzyme_dict = EnzymeDict()
        self._state = 'title'
        self.key_waiting = ''
        self.enzyme_dict[ 'factoids' ] = []
        self.enzyme_dict[ 'number_of_sites' ] = dict()
        self._table_nest_level = 0
        self._table_nest_dir = ''
        self._factoid = ''

    def parse(self, handle):
        self.reset()
        self.feed(handle)
        return self.enzyme_dict

    def feed(self, handle):
        """feed(self, handle )

        Feed in rebase data for scanning.  handle is a file-like object
        containing rebase data.  consumer is a Consumer object that will
        receive events as the rebase data is scanned.

        """
        if isinstance(handle, File.UndoHandle):
            uhandle = handle
        else:
            uhandle = File.UndoHandle(handle)
        text = ''
        while 1:
            line = uhandle.readline()
            line = string.strip( line )
            if( is_blank_line( line, 0 ) ):
                break
            if( line[ -7: ] == '</HTML>' ):
                break
            text = text + ' ' + line

        sgmllib.SGMLParser.feed( self, text )


    def handle_data(self, newtext ):
        newtext = string.strip( newtext )
        self.text = self.text + newtext

    def start_table( self, attrs ):
        self._table_nest_level = self._table_nest_level + 1
        self._table_nest_dir = 'increasing'

    def end_table( self ):
        self._table_nest_level = self._table_nest_level - 1
        if( self._table_nest_dir == 'increasing' ):
            if( self._state == 'details' ):
                self.extract_details()
                self._state = 'number_of_sites'
            elif( self._state == 'number_of_sites' ):
                self.extract_number_of_sites()
                self._state = 'complete'
        self._table_nest_dir = 'decreasing'

    def start_b( self, attrs ):
        if( self._state == 'details' ):
            self.extract_details()
        else:
            self.flush_text()

    def end_b( self ):
        if( self._state == 'details' ):
            text = self.flush_text()
            if( text != '' ):
                if( text[ -1 ] == ':' ):
                    key = self.key_waiting
                    if( key != '' ):
                        self.enzyme_dict[ key ] = self._factoid
                        self._factoid = ''
                    self.key_waiting = text[ :-1 ]
                else:
                    self._factoid = self._factoid + ' ' + text

    def start_i( self, attrs ):
        if( self._state == 'number_of_sites' ):
            self.extract_number_of_sites()

    def end_i( self ):
        if( self._state == 'details' ):
            self.extract_details()
        elif( self._state == 'number_of_sites' ):
            text = self.flush_text()
            if( text != '' ):
                if( text[ -1 ] == ':' ):
                    self.key_waiting = text[ :-1 ]
                elif( self.key_waiting != '' ):
                    key = self.key_waiting
                    self.key_waiting = ''
                    self.enzyme_dict[ 'number_of_sites' ][ key ] = text

    def start_pre( self, attrs ):
        self.flush_text()
        if( self._state == 'sequence' ):
            self._state = '5_to_3_strand'

    def end_pre( self ):
        if( self._state == '5_to_3_strand' ):
            text = self.flush_text()
            cols = string.split( text, '...' )
            if( len( cols ) > 1 ):
                self.enzyme_dict[ '5_to_3_strand' ] = cols[ 1 ]
                self.enzyme_dict[ '3_to_5_strand' ] = ''

            self._state = 'details'

        elif( self._state == '3_to_5_strand' ):
            text = self.flush_text()
            self.enzyme_dict[ '3_to_5_strand' ] = text
            self._state = 'details'

    def do_br( self, attrs ):
        if( self._state == '5_to_3_strand' ):

            text = self.flush_text()
            if( text != '' ):
                self.enzyme_dict[ '5_to_3_strand' ] = text
                self._state = '3_to_5_strand'
        elif( self._state == 'details' ):
            self.extract_details()
            key = self.key_waiting
            if( self._factoid != '' ):
                if( key != '' ):
                    self.enzyme_dict[ key ] = self._factoid
                    self.key_waiting = ''
                else:
                    self.enzyme_dict[ 'factoids' ].append( self._factoid )
                self._factoid = ''

    def start_font( self, attrs ):
        self.flush_text()

    def end_font( self ):

        if( self._state == 'title' ):
            text = self.flush_text()
            if( text != '' ):
                self._state = 'descriptor'
                self.enzyme_dict[ 'title' ] = text
        elif( self._state == 'descriptor' ):
            text = self.flush_text()
            if( text != '' ):
                self._state = 'sequence'
                self.enzyme_dict[ 'descriptor' ] = text

    def extract_details( self ):
        text = self.flush_text()
        if( text != '' ):
            self._factoid = self._factoid + ' ' + text

    def extract_number_of_sites( self ):
        text = self.flush_text()
        if( text != '' ):
            key = self.key_waiting
            if( key != '' ):
                self.enzyme_dict[ 'number_of_sites' ][ key ] = text
                self.key_waiting = ''


    def flush_text( self ):
        text = string.strip( self.text[:] )
        self.text = ''
        return text


if( __name__ == '__main__' ):
    handle = open( 'bamii.htm')
    undo_handle = Bio.File.UndoHandle( handle )
    rebase_parser = RebaseParser()
    data = rebase_parser.parse( handle )
    print data
