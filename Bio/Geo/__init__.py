# Copyright 2001 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# standard library
import string
import array
import os
import re
import sgmllib
import urlparse
import copy

# XML from python 2.0
from xml.sax import handler

# Martel
import Martel
from Martel import RecordReader
from Martel import Dispatch

from Bio.ParserSupport import EventGenerator
from Bio.ParserSupport import AbstractConsumer
from Bio import File
from Bio.Align.Generic import Alignment
import Bio.Alphabet
import geo_format
import Record

__all__ = [
    'Record',
    'geo_format'
    ]

class Error( Exception ):
    """
    """
    def __init__( self ):
        pass

class GeoError( Error ):

    """
        message - description of error
    """

    def __init__( self, message ):
        self.message = message

class Iterator:
    """Iterator interface to move over a file of Geo entries one at a time.
    """
    def __init__(self, handle, parser = None):
        """Initialize the iterator.

        Arguments:
        o handle - A handle with GEO entries to iterate through.
        o parser - An optional parser to pass the entries through before
        returning them. If None, then the raw entry will be returned.
        """
        self.handle = File.UndoHandle( handle )
        self._reader = RecordReader.Everything( self.handle )
        self._parser = parser

    def next(self):
        """Return the next GEO record from the handle.

        Will return None if we ran out of records.
        """
        data = self._reader.next()

        if self._parser is not None:
            if data:
                dumpfile = open( 'dump', 'w' )
                dumpfile.write( data )
                dumpfile.close()
                return self._parser.parse(File.StringHandle(data))

        return data
    
    def __iter__(self):
        return iter(self.next, None)



class _Scanner:
    """Start up Martel to do the scanning of the file.

    This initialzes the Martel based parser and connects it to a handler
    that will generate events for a Feature Consumer.
    """
    def __init__(self, debug = 0):
        """Initialize the scanner by setting up our caches.

        Creating the parser takes a long time, so we want to cache it
        to reduce parsing time.

        Arguments:
        o debug - The level of debugging that the parser should
        display. Level 0 is no debugging, Level 2 displays the most
        debugging info (but is much slower). See Martel documentation
        for more info on this.
        """
        # a listing of all tags we are interested in scanning for
        # in the MartelParser
        self.interest_tags = [ 'entity_line', 'attribute_line', 'col_heading_line', \
            'row_line' ]

        # make a parser that returns only the tags we are interested in
        expression = Martel.select_names( geo_format.geo_record, self.interest_tags)
        self._parser = expression.make_parser(debug_level = debug)

    def feed(self, handle, consumer):
        """Feed a set of data into the scanner.

        Arguments:
        o handle - A handle with the information to parse.
        o consumer - The consumer that should be informed of events.
        """
        consumer.set_interest_tags( self.interest_tags )
        self._parser.setContentHandler( consumer )
#        self._parser.setErrorHandler(handle.ErrorHandler())

        self._parser.parseFile(handle)

class _RecordConsumer( Dispatch.Dispatcher ):
    """Create a Saf Record object from scanner generated information.
    """
    def __init__(self):
        Dispatch.Dispatcher.__init__( self )
        self.data = Record.Record()
        self._refresh()

    def _refresh( self ):
        self._col_headers = []

    def set_interest_tags( self, interest_tags ):
        self.interest_tags = interest_tags

    def startDocument(self):
        self.data = Record.Record()
        self._refresh()


    def start_entity_line(self, text, attrs):
        self.save_characters()

    def end_entity_line(self, text ):
        line = self.get_characters()
        cols = line.split( '=' )
        entity_type = ( cols[ 0 ] ).strip()
        entity_type = entity_type[ 1: ]
        self.data.entity_type = entity_type
        entity_id = cols[ 1 ].strip()
        self.data.entity_id = entity_id


    def start_attribute_line(self, text, attrs):
        self.save_characters()

    def end_attribute_line(self, text ):
        line = self.get_characters()
        cols = line.split( '=', 1 )
        key = cols[ 0 ].strip()
        key = key[ 1: ]
        val = cols[ 1 ]
        val = val.strip()

        if( self.data.entity_attributes.has_key( key ) ):
            contents = self.data.entity_attributes[ key ]
            if( type( contents ) == type( [] ) ):
                contents.append( val )
                lazy_list = contents
            else:
                lazy_list =  [ contents, val ]
            self.data.entity_attributes[ key ] = lazy_list
        else:
            self.data.entity_attributes[ key ] = val

    def start_col_heading_line(self, text, attrs):
        self.save_characters()

    def end_col_heading_line(self, text ):
        line = self.get_characters()
        line = line[ 1: ]
        cols = line.split( '=' )
        key = cols[ 0 ].strip()
        val = cols[ 1 ].strip()
        self.data.col_defs[ key ] = val

    def start_row_line(self, text, attrs):
        self.save_characters()

    def end_row_line(self, text ):
        line = self.get_characters()
        contents = line.strip()
        if( contents != '' ):
            cols = line.split( '\t' )
            cols = map( string.strip, cols )
            self.data.table_rows.append( cols )



class RecordParser:
    """Parse GEO files into Record objects
    """
    def __init__(self, debug_level = 0):
        """Initialize the parser.

        Arguments:
        o debug_level - An optional argument that specifies the amount of
        debugging information Martel should spit out. By default we have
        no debugging info (the fastest way to do things), but if you want
        you can set this as high as two and see exactly where a parse fails.
        """
        self._scanner = _Scanner(debug_level)

    def parse(self, handle):
        """Parse the specified handle into an GEO record.
        """
        self._consumer = _RecordConsumer()
        self._scanner.feed(handle, self._consumer)
        return self._consumer.data

