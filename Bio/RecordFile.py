# Copyright 2001 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Code for more fancy file handles.

Classes:
RecordFile is a decorator for File that allows the user to skip over 
boilerplate and read the contents of a record.  The initializer requires the
starting tag and the ending tag of the record.  RecordFile processes multiple
records, provided they all have the same starting and ending tags.

The implementation is based on a state machine. and assumes sequential access.
RecordRead.read provides the same interface as read.  It has an optional
parameter, size.  

RecordFile.read, with no parameters, searches for the next record and returns
it in its entirety.  A subsequent call returns an empty string to signal the
end of the record.  However, the next call will advance to the next record
and return it, it it exists.  Otherwise RecordFile.read will return a second
 empty string.

If a calling program passes a size parameter RecordFile will check its own
position with respect to record data. If it just past the end of record, it
will return an empty string.  It it is between records it will advance to
the next record and retuurn the specified number of bytes from the record,
if they are available.  Otherwise, it will return all the data remaining in
the current record.  If it is already within the record, RecordFile.read
will return data from the current position.

"""

import os
import string
from File import UndoHandle



"""Used for debugging"""
def dump_saved( name, text, j ):
    dump_file = open( name + '%d' % j, "w" )
    k = 0
    for i in range ( 0, len( text ), 80 ):
        dump_file.write(  '%s\n' % text[ i : i + 80 ] )
    dump_file.close()




class RecordFile:
    def __init__(self, handle, start_tag, end_tag ):
        self._handle = handle
        self._start_tag = start_tag
        self._start_len = len( self._start_tag )
        self._end_tag = end_tag
        self._end_len = len( self._end_tag )
        self._set_state_info()
        self._valid_delims = [ chr( 10 ), chr( 13 ), chr( 10 ) + chr( 13 ), chr( 13 ) + chr( 10 ) ]
        self._debug_count = 0

    def _set_state_info( self ):
        self._file_state = 'CLOSED'
        self._record_state = 'SEARCHING'
        self._saved_text = ''
        self._look_back = ''
        self._tag_chars_pending = 0


    def restart( self ):
        self._set_state_info()


    def close(self, *args, **keywds ):
        return apply(self._handle.close, args, keywds)

    def read( self, *args, **keywds ):
        line = ''
        if( self._record_state == 'SEARCHING' ):
            line = self._search_start()
            self._saved_text = line + self._saved_text
        if( self._record_state == 'IN_RECORD' ):
            text = self._in_record_state( args, keywds )

        elif( self._record_state == 'SCANNING_END_TAG' ):
            len_expected = self._get_len_expected( args, keywds )
            if( len_expected ):
                text = self._scan_tag( self._saved_text[ : len_expected ] )
            else:
                text = self._scan_tag( self._saved_text[ : ] )

        elif( self._record_state == 'AT_END_RECORD' ):
            self._record_state = 'SEARCHING'
            text =  ''
        elif( self._record_state == 'EXIT' ):
            text = ''
        else:
            text = ''
        self._debug_count = self._debug_count + 1
        return text

    def _in_record_state( self, args, keywds ):
        saved_text = self._saved_text
        lookahead_len = self._end_len + 1
        adjustment = lookahead_len - len( self._saved_text )
        len_expected = self._get_len_expected( args, keywds )
        if( len_expected ):
            len_to_read = len_expected + adjustment
            if( len_to_read > 0 ):
                text_read = self._handle.read( len_to_read )
                text = saved_text + text_read
                if( len( text_read ) < len_to_read ):
                    self._file_state = 'AT_END_FILE'
            else:
                len_to_retrieve = len_expected + lookahead_len
                text = saved_text[ : len_to_retrieve ]
                requested_text = text[ :len_expected ]
        else:
            self._file_state = 'AT_END_FILE'
            requested_text = text
        if( text == '' ):
            self._file_state = 'AT_END_FILE'
            self._record_state = 'EXIT'
        else:
            requested_text = self._search_end( text )
            self._saved_text = text[ len( requested_text ) : ]

        return requested_text

    def _get_len_expected( self, args, keywds ):

        if( len( args) > 0 ):

            len_expected = args[ 0 ]
        elif( keywds.has_key( 'size' ) ):
            len_expected = keywds[ 'size' ]
        else:
            len_expected = None
        return len_expected

    def _scan_tag( self, text ):
        text_len = len( text )
        tag_chars_pending = self._tag_chars_pending
        if( text_len >= tag_chars_pending ):
            if( self._record_state == 'SCANNING_START_TAG' ):
                self._record_state = 'IN_RECORD'
            elif( self._record_state == 'SCANNING_END_TAG' ):
                self._record_state = 'AT_END_RECORD'
            self._tag_chars_pending = 0
            self._saved_text = text[ tag_chars_pending : ]
            requested_text = text[ : tag_chars_pending ]
        else:
            self._scan_tag_index = self_scan_tag_index + text_len
            requested_text = text
        return requested_text


    def _search_start( self ):
        while(  self._record_state == 'SEARCHING' ):
            line = self.extract_saved_line( self._valid_delims )
            if( line == '' ):
                line = self._handle.readline()
                line = self._saved_text + line
                self._saved_text = ''
            if( line == '' ):
                self_file_state = 'AT_END_FILE'
                self._record_state = 'EXIT'
            elif( line.startswith( self._start_tag ) ):
                self._record_state = 'IN_RECORD'
        return( line )

    def _search_end( self, text ):
        pos = text.find( self._end_tag )
        lookahead_len = self._end_len + 1
        requested_text_len = len( text )
        if( self._file_state != 'AT_END_FILE' ):
            requested_text_len = requested_text_len - lookahead_len
        requested_text = text[ : requested_text_len ]
        len_to_end = self.find_record_boundary( text, requested_text_len )
        requested_text = text[ : len_to_end ]
        requested_text_len = len( requested_text )
        self._look_back = text[ requested_text_len -2 : requested_text_len ]
        return requested_text

    def find_record_boundary( self, text, text_len ):
        len_to_end  = text_len
        look_back = self._look_back
        pos = 0
        accum_pos = 0
        newpos = 0
        while( pos >= 0  and pos < text_len ):
            acum_pos = accum_pos + pos
            pos = text.find( self._end_tag, newpos )
            current_pos = accum_pos + pos
            newpos = pos + self._end_len
            suffix = text[ current_pos + self._end_len :  ]
            if( len( suffix ) == 0 ):
                len_to_end = current_pos + self._end_len
                self._record_state = 'AT_END_RECORD'
                break
            elif( ( pos >= 0 ) and ( current_pos < text_len ) ):
                ( delim_pos, delim_len )= find_delim( suffix, self._valid_delims )
                if( delim_len > 0 ):
                    prefix = look_back + text[ : current_pos ]
                    if( self.at_line_boundary( prefix, delim_pos ) ):
                        len_to_end =  current_pos + self._end_len + delim_len
                        len_to_end = min( len_to_end, text_len )
                        self._saved_text = text[ len_to_end : ] + self._saved_text
                        if( len_to_end <= text_len ):
                            self._record_state = 'AT_END_RECORD'
                        else:
                            self._record_state = 'SCANNING_END_TAG'
                            scan_tag_index = len_to_end - text_len
                            self._tag_chars_pending = self._end_len + delim_len - scan_tag_index
                        break
        return len_to_end

    def at_line_boundary( self, prefix, delim_pos ):
        if( is_prefix_in_set( prefix, self._valid_delims ) or \
            (  delim_pos == 0 ) ):
            return 1
        else:
            return 0

    def extract_saved_line( self, delims ):
        text = self._saved_text
        if( text == '' ):
            return text
        ( delim_pos, delim_len )= find_delim( text, delims )
        if( delim_pos >= 0 and delim_pos < len( text ) ):
            line_len = delim_pos + delim_len
            line = text[ : line_len ]
            self._saved_text = text[ line_len : ]
        else:
            line = ''
        return line

def find_delim( text, delims ):
    pos = len( text )
    delim_len = 0
    for delim in delims:
        first = text.find( delim )
        if( first >= 0 ):
            if( first < pos ):
                pos = first
                delim_len = len( delim )
            elif( first == pos ):
                delim_len = max( delim_len, len( delim ) )
    return ( pos, delim_len )

def is_prefix_in_set( text, delims ):
    pos = len( text )
    for delim in delims:
        last = text.rfind( delim )
        if( ( last + len( delim ) ) == len( text ) ):
            pos = min( pos, last )
    if( pos < len( text ) ):
        return 1
    else:
        return 0


