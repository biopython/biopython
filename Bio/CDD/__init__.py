# Copyright 2002 by Katharine Lindner.  All rights reserved.
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


# XML from python 2.0
from xml.sax import handler

# Martel
import Martel
from Martel import RecordReader

from Bio.RecordFile import RecordFile
from Bio.FilteredReader import FilteredReader
from Bio.FilteredReader import remove_empty_line
from Bio.FilteredReader import remove_leading_whitespace
from Bio.SGMLExtractor import SGMLExtractorHandle
from Bio import File
from Bio.Seq import Seq
from Martel.Dispatch import Dispatcher
import cdd_format
import Record

__all__ = [
    'Record',
    'cdd_format'
    'Iterator'
    ]



class Iterator:
    """Iterator interface to move over a file of CDD entries one at a time.
       Iterator expects a handle to an sgml file.  It extracts data bracketed
       by specified tag pairs, then removes blank lines and leading white space.
       The parser operates on the filtered data.
    """
    def __init__(self, handle, parser = None):
        """Initialize the iterator.

        Arguments:
        o handle - A handle with CDD entries to iterate through.
        o parser - An optional parser to pass the entries through before
        returning them. If None, then the raw entry will be returned.
        """
        record_handle = SGMLExtractorHandle( handle, [ 'title', 'table', ] )
        filtered_handle = FilteredReader( record_handle )
        filtered_handle.filter_chain = [ remove_empty_line, remove_leading_whitespace ]
        self.handle = File.UndoHandle( filtered_handle )
        self._reader = RecordReader.Everything( self.handle  )
        self._parser = parser

    def next(self):
        """Return the next CDD record from the handle.

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

class _Scanner:
    """Start up Martel to do the scanning of the file.

    This initialzes the Martel based parser and connects it to a handler
    that will generate events for a Feature Consumer.
    """
    def __init__(self, debug_level = 0):
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
        self.interest_tags = [ "cd_tag", \
            "description_tag", \
            "status_tag", \
            "source_tag", \
            "date_tag", \
            "taxonomy_tag", \
            "aligned_tag", \
            "representative_tag", \
            "range_tag", \
            "sequence_tag", \
            "description_contents_multiline", \
            "status_contents_multiline", \
            "source_contents_multiline", \
            "date_contents_multiline", \
            "reference_contents_multiline", \
            "taxonomy_contents_multiline", \
            "aligned_contents_multiline", \
            "representative_contents_multiline", \
            "range_contents_multiline", \
            "cd_contents_multiline", \
            "sequence_contents_multiline", \
            "table_entry" ]

        # make a parser that returns only the tags we are interested in
        expression = Martel.select_names( cdd_format.cdd_record, self.interest_tags)
        self._parser = expression.make_parser(debug_level )

    def feed(self, handle, consumer):
        """Feeed a set of data into the scanner.

        Arguments:
        o handle - A handle with the information to parse.
        o consumer - The consumer that should be informed of events.
        """
        consumer.set_interest_tags( self.interest_tags )
        self._parser.setContentHandler( consumer )
#        self._parser.setErrorHandler(handle.ErrorHandler())

        self._parser.parseFile(handle)

class _RecordConsumer( Dispatcher ):
    """Create a CDD Record object from scanner generated information.
    """
    def __init__(self):
        Dispatcher.__init__( self )
        self.data = Record.Record()
        self._pending_key = ''


    def set_interest_tags( self, interest_tags ):
        self.interest_tags = interest_tags

    def start_cd_tag( self, line, attrs ):
        self.save_characters()

    def end_cd_tag( self, cdd_record ):
        key = self.save_key()

    def start_cd_contents_multiline( self, text, attrs ):
        self.save_characters()

    def end_cd_contents_multiline( self, cdd_record ):
        self.add_entry()

    def start_description_tag( self, text, attrs ):
        self.save_characters()

    def end_description_tag( self, cdd_record ):
        key = self.save_key()

    def start_description_contents_multiline( self, text, attrs ):
        self.save_characters()

    def end_description_contents_multiline( self, cdd_record ):
        self.add_entry()

    def start_status_tag( self, text, attrs ):
        self.save_characters()

    def end_status_tag( self, cdd_record ):
        key = self.save_key()

    def start_status_contents_multiline( self, text, attrs ):
        self.save_characters()

    def end_status_contents_multiline( self, cdd_record ):
        self.add_entry()

    def start_source_tag( self, text, attrs ):
        self.save_characters()

    def end_source_tag( self, cdd_record ):
        key = self.save_key()

    def start_source_contents_multiline( self, text, attrs ):
        self.save_characters()

    def end_source_contents_multiline( self, cdd_record ):
        self.add_entry()

    def start_date_tag( self, text, attrs ):
        self.save_characters()

    def end_date_tag( self, cdd_record ):
        key = self.save_key()

    def start_date_contents_multiline( self, text, attrs ):
        self.save_characters()

    def end_date_contents_multiline( self, cdd_record ):
        self.add_entry()

    def start_reference_contents_multiline( self, text, attrs ):
        self.save_characters()

    def end_reference_contents_multiline( self, cdd_record ):
        reference = self.get_characters()
        self.data[ 'references' ].append( reference )

    def start_taxonomy_tag( self, text, attrs ):
        self.save_characters()

    def end_taxonomy_tag( self, cdd_record ):
        key = self.save_key()

    def start_taxonomy_contents_multiline( self, text, attrs ):
        self.save_characters()

    def end_taxonomy_contents_multiline( self, cdd_record ):
        self.add_entry()

    def start_aligned_tag( self, text, attrs ):
        self.save_characters()

    def end_aligned_tag( self, cdd_record ):
        key = self.save_key()

    def start_aligned_contents_multiline( self, text, attrs ):
        self.save_characters()

    def end_aligned_contents_multiline( self, cdd_record ):
        self.add_entry()

    def start_representative_tag( self, text, attrs ):
        self.save_characters()

    def end_representative_tag( self, cdd_record ):
        key = self.save_key()

    def start_representative_contents_multiline( self, text, attrs ):
        self.save_characters()

    def end_representative_contents_multiline( self, cdd_record ):
        self.add_entry()

    def start_range_tag( self, text, attrs ):
        self.save_characters()

    def end_range_tag( self, cdd_record ):
        key = self.save_key()

    def start_range_contents_multiline( self, text, attrs ):
        self.save_characters()

    def end_range_contents_multiline( self, cdd_record ):
        self.add_entry()

    def start_sequence_tag( self, text, attrs ):
        self.save_characters()

    def end_sequence_tag( self, cdd_record ):
        key = self.save_key()

    def start_sequence_contents_multiline( self, text, attrs ):
        self.save_characters()

    def end_sequence_contents_multiline( self, cdd_record ):
        line = self.get_characters()
        ( lines )  = line.splitlines()
        key = self._pending_key
        val = ''
        for line in lines:
            line = line.strip()
            val = val + line
        self.data[ key ] = Seq( val )

    def start_table_entry( self, text, attrs ):
        self.save_characters()

    def end_table_entry( self, cdd_record ):
        line = self.get_characters()
        ( lines )  = line.splitlines()
        key = ''
        val = ''
        state = 'key'
        for line in lines:
            line = line.strip()
            upper_line = line.upper()
            if( upper_line.endswith( '[CD]' ) ):
                line = line[ :-4 ]
                state = 'val'
            elif( len( line ) > 60 ):
                state = 'val'
            else:
                state = 'key'
            if( state == 'key' ):
                key = key + line
            else:
                val = val + line
        self.data[ 'alignment_lookup' ][ key ] = val

    def save_key( self ):
        key = self.get_characters()
        self._pending_key = key[ : -1 ]

    def add_entry( self ):
        key = self._pending_key
        self._pending_key = ""
        self.data[ key ] = self.get_characters()

class RecordParser:
    """Parse CDD files into Record objects
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
        """Parse the specified handle into an NBRF record.
        """
        self._consumer = _RecordConsumer()
        self._scanner.feed(handle, self._consumer)
        return self._consumer.data

