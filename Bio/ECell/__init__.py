# Copyright 2001 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""For reading the ECell spreadsheet format from Ecell2 (DEPRECATED).

Ecell converts the ECell input from spreadsheet format to an intermediate
format, described in http://www.e-cell.org/manual/chapter2E.html#3.2.  It
provides an alternative to the perl script supplied with the Ecell2
distribution at http://bioinformatics.org/project/?group_id=49.

ECell expects a spreadsheet exported in delimited text format. The file should
be read with FilteredReader using the default filter chain to remove extraneous
characters.
"""
import warnings
warnings.warn("Bio.ECell was deprecated, as it does not seem to have any" \
              + " users. If you do use this module, please contact the"\
              + " Biopython developers at biopython-dev@biopython.org to"\
              + " avoid permanent removal of this module", DeprecationWarning)

# standard library
import sys
import string
import copy
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

from Bio.ParserSupport import EventGenerator
from Bio.ParserSupport import AbstractConsumer
from Bio import File
from Bio.Align.Generic import Alignment
import Bio.Alphabet
import ecell_format
import Record

class Error( Exception ):
    """
    """
    def __init__( self ):
        pass

class ECellError( Error ):

    """
        message - description of error
    """

    def __init__( self, message ):
        self.message = message



class Iterator:
    """Iterator interface to move over a file of ecell entries one at a time.
    """
    def __init__(self, handle, parser = None):
        """Initialize the iterator.

        Arguments:
        o handle - A handle with ECell entries to iterate through.
        o parser - An optional parser to pass the entries through before
        returning them. If None, then the raw entry will be returned.
        """
        self.handle = File.UndoHandle( handle )
        self._reader = RecordReader.Everything( self.handle )
        self._parser = parser

    def next(self):
        """Return the next ecell record from the handle.

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
        self.interest_tags = [ 'header_line', 'system_line', 'substance_multiline', \
            'reactor_multiline', 'include_line' ]

        # make a parser that returns only the tags we are interested in
        expression = Martel.select_names( ecell_format.ecell_record, self.interest_tags)
        self._parser = expression.make_parser(debug_level = debug)

    def feed(self, handle, consumer):
        """Feed a set of data into the scanner.

        Arguments:
        o handle - A handle with the information to parse.
        o consumer - The consumer that should be informed of events.
        """
        self._parser.setContentHandler( EventGenerator(consumer,
                                                       self.interest_tags))
#        self._parser.setErrorHandler(handle.ErrorHandler())

        self._parser.parseFile(handle)

class _RecordConsumer:
    """Create an ECell Record object from scanner generated information.
    """
    def __init__(self):
        self.data = Record.Record()
        self._header = []
        self._database = {}
        self._state = ''

    def include_line( self, line ):
        self.data.include_buf = self.data.include_buf + line

    def header_line( self, lines ):
        for line in lines:
            items = line.split( '\t')
            items[ 0 ] = items[ 0 ].lower()
            self._header = []
            self._state = items[ 0 ]
            for item in items:
                item = item.strip()
                self._header.append( item.lower() )


    def system_line( self, lines ):
        for line in lines:
            line_dict = self._make_line_dict( line )
            if( not self._check_missing_header( line_dict ) ):
                raise EcellError( "invalid header" )
            self.data.num_systems = self.data.num_systems + 1
            _set_defaults( line_dict )
            self._build_system_entry( line_dict )


    def substance_multiline( self, multiline ):
        for line in multiline:
            self.parse_substance_lines( line )

    def parse_substance_lines( self, multiline ):
        lines = multiline.splitlines()
        line_no = 0
        for line in lines:
            line_dict = self._make_line_dict( line )
            try:
                if( not _is_valid_substance( line_dict ) ):
                    raise ECellError( "quantity and concentration are mutually exclusive" )
            except ECellError, e:
                print sys.stderr, e.message

            qty = Record.get_entry( line_dict, 'qty' )
            conc = Record.get_entry( line_dict, 'conc' )
            if( ( qty.lower() != 'fix' ) and ( conc.lower() != 'fix' ) ):
                self.data.num_substances = self.data.num_substances + 1
            else:
                line_no = line_no + 1
            if( line.lower().startswith( 'substance' ) ):
                _set_defaults( line_dict )
                self._convert_conc( line_dict )

            self._build_substance_entry( line_dict, line_no )

    def reactor_multiline( self, multiline ):
        for line in multiline:
            self.parse_reactor_lines( line )

    def parse_reactor_lines( self, multiline ):
        lines = multiline.splitlines()
        for line in lines:
            line_dict = self._make_line_dict( line )
            if( line.lower().startswith( 'reactor' ) ):
                if( not self._check_missing_header( line_dict ) ):
                    raise ECellError( "invalid header" )
            try:
                if( not is_only_digits( line_dict[ 's_coeff' ] ) ):
                    raise ECellError( 's_coeff must contain only digits' )
                if( not is_only_digits( line_dict[ 'p_coeff' ] ) ):
                    raise ECellError( 'p_coeff must contain only digits' )
            except KeyError:
                pass
            if( line.lower().startswith( 'reactor' ) ):
                _set_reactor_defaults( line_dict )
            line_dict = self._remove_if_inconsistent( line_dict )

            if( line_dict.has_key( 'class' ) ):
                self.data.num_reactors = self.data.num_reactors + 1
                num_substrates = 0
                num_products = 0
                num_catalysts = 0
                num_effectors = 0
                num_options = 0
                num_args = 0
            if( line_dict.has_key( 's_id' ) ): num_substrates = num_substrates + 1
            if( line_dict.has_key( 'p_id' ) ): num_products = num_products + 1
            if( line_dict.has_key( 'c_id' ) ): num_catalysts = num_catalysts + 1
            if( line_dict.has_key( 'e_id' ) ): num_effectors = num_effectors + 1
            if( line_dict.has_key( 'o_type' ) ): num_options = num_options + 1
            if( line_dict.has_key( 'arg_tag' ) ): num_args = num_args + 1
            counter_dict = { \
                's_' : num_substrates, \
                'p_' : num_products, \
                'c_' : num_catalysts, \
                'e_' : num_effectors, \
                'o_' : num_options, \
                'arg_tag' : num_args
            }
            self._set_max( counter_dict )
            self._build_reactor_entry( line_dict, counter_dict )


    def _set_max( self, counter_dict ):
        num_reactors = self.data.num_reactors
        for key in counter_dict.keys():
            composite_key = key + str( num_reactors )
            self.data._max_dict[ composite_key ] = counter_dict[ key ]

    def _build_system_entry( self, line_dict ):
        for key in line_dict.keys():
            item = line_dict[ key ]
            composite_key = 'system' + str( self.data.num_systems ) + key + '0'

            if( not self.data.cell_dict.has_key( composite_key ) ):
                self.data.cell_dict[ composite_key ] = item

    def _build_substance_entry( self, line_dict, line_no ):
        for key in line_dict.keys():
            item = line_dict[ key ]
            composite_key = 'substance' + str( self.data.num_substances ) + key + \
                str( line_no )
            if( not self.data.cell_dict.has_key( composite_key ) ):
                self.data.cell_dict[ composite_key ] = item

    def _convert_conc( self, line_dict ):
        if( line_dict.has_key( 'conc' ) ):
            if( not line_dict.has_key( 'qty' ) ):
                contents = 'QTY(%s,%s)' % ( line_dict[ 'conc' ], line_dict[ 'path' ] )
                composite_key = 'substance' + str( self.data.num_substances ) + 'qty' + '0'
                self.data.cell_dict[ composite_key ] = contents
                self.data.contains_concentration = 1

    def _build_reactor_entry( self, line_dict, counter_dict ):
        for key in line_dict.keys():
            item = line_dict[ key ]
            prefix = key[ :2 ]
            if( key.startswith( 'arg_' ) ):
                index = counter_dict[ 'arg_tag' ]
            elif( counter_dict.has_key( prefix ) ):
                index = counter_dict[ prefix ]
            else:
                index = '0'
            composite_key = 'reactor' + str( self.data.num_reactors ) + str( key ) + str( index )
            if( not self.data.cell_dict.has_key( composite_key ) ):
                self.data.cell_dict[ composite_key ] = item


    def _check_missing_header( self, line_dict ):
        ok = 1
        items = [ 'id', 'path', 'class' ]
        for item in items:
            if( line_dict.has_key( item ) == 0 ):
                others = copy.deepcopy( items )
                others.remove( item )
                for other in others:
                    if( line_dict.has_key( other ) ):
                        if( item.lower() != 'class' ):
                            ok = 0
                            break
        return ok

    def _remove_if_inconsistent( self, list_dict ):
        valid_keys = list_dict.keys()
        for label in [ 'id', 'path', 'type' ]:
            for prefix in [ 's_', 'p_', 'c_', 'e_' ]:
                node = prefix + label
                valid_keys = self._consistency_filter( prefix, node, valid_keys )
        for key in list_dict.keys():
            if( not key in valid_keys ):
                del list_dict[ key ]
        return list_dict

    def _consistency_filter( self, prefix, tag, valid_keys ):
        block = []
        for suffix in [ 'id', 'path', 'coeff', 'type' ]:
            node = prefix + suffix
            block.append( node )
        for node in block:
            if( ( not tag in valid_keys ) and ( node in valid_keys ) ):
                if( ( prefix == 'o_' ) or ( not tag.endswith( 'type' ) ) ):
                    valid_keys.remove( node )
        return valid_keys

    def _make_line_dict( self, line ):
        line_dict = {}
        items = line.split( '\t' )
        num = 0
        for item in items:
            item = item.strip()
            if( item != '' ):
                line_dict[ self._header[ num ] ] = item
            num = num + 1
        return line_dict

def _clear_bad_block( block, items ):
    for label in block:
        items = items.remove( items.index( label ) )
    return items

def _is_valid_substance( line_dict ):
    ok = 1
    if( line_dict.has_key( 'qty' ) and line_dict.has_key( 'conc' ) ):
        if( not ( line_dict[ 'qty' ] == 'QTY' ) ):
            ok = 0
    return ok

def is_only_digits( line ):
    ok = 1
    text = line.strip()
    if( text != '' ):
        if( not text.isdigit() ):
            ok = 0
    return ok

def _set_reactor_defaults( line_dict ):
    line_dict = _set_defaults( line_dict )
    for item in [ 's_', 'p_', 'c_', 'e_' ]:
        id = item + 'id'
        coeff = item + 'coeff'
        path = item + 'path'
        if( line_dict.has_key( id ) ):
            if( not line_dict.has_key( coeff ) ):
                line_dict[ coeff ] = 1
        if( not line_dict.has_key( path ) ):
            line_dict[ path ] = line_dict[ 'path' ]

    return( line_dict )

def _set_defaults( line_dict ):
    if( not line_dict.has_key( 'name' ) ):
        line_dict[ 'name' ] = line_dict[ 'id' ]
    if( line_dict.has_key( 'arg_tag' ) ):
        if( not line_dict.has_key( 'arg_coeff' ) ):
            line_dict[ 'arg_coeff' ] = 0

    return( line_dict )







class RecordParser:
    """Parse ECell files into Record objects
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
        """Parse the specified handle into an ECell record.
        """
        self._consumer = _RecordConsumer()
        self._scanner.feed(handle, self._consumer)
        return self._consumer.data

