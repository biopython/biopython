# Copyright 2001 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Code for more fancy file handles.

Classes:
Filtered is a decorator for File that allows the user to filter the output
on a line by line basis.
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


def remove_empty_line( line ):
    stripped_line = line.strip()
    if( stripped_line ):
        return line[ : ]
    else:
        return ''

def remove_useless_dot( line ):
    before = line
    while( 1 ):
        after = before.replace( "\t.\t", "\t\t" )
        if( len( before ) == len( after ) ):
            break
        before = after
    if( after.endswith( '.' ) ):
        after = after[ :-1 ]
    return after

def fix_punctuation( line ):
    line = line.replace( "'", '' )
    line = line.replace( '"', '' )
    line = line.replace( ';', '\t' )
    line = line.replace( 'entryname', 'id' )
#    line = line.lower( )
    if( line ):
        return line[ : ]
    else:
        return ''



class FilteredReader:
    def __init__(self, handle ):
        self._handle = handle
        self._start_line = ''
        self._debug_count = 0

    def __getattr__(self, attr):
        return getattr(self._handle, attr)



    def close(self, *args, **keywds ):
        return apply(self._handle.close, args, keywds)

    def read( self, *args, **keywds ):
        line = ''
        len_expected = self._get_len_expected( args, keywds )
        if( len_expected ):
            filtered_text = self.read_block( len_expected )
        else:
            filtered_text = self.read_to_end()
        return filtered_text

    def read_block( self, len_expected ):

        len_filtered = 0
        len_adjusted = len_expected - len( self._start_line )
        filtered_text = ''
        while( len_filtered < len_expected ):

            text_read = self._handle.read( len_adjusted )
            full_text = self._start_line + text_read
            filtered_text = filtered_text + self._filter( full_text )
            if( text_read == '' ):
                break
            len_filtered_text = len( filtered_text )
            len_adjusted = len_adjusted - len_filtered_text
        return filtered_text[ : ]

    def read_to_end( self ):
        filtered_text = ''
        text_read = self._handle.read()
        full_text = self._start_line + text_read
        filtered_text = filtered_text + self._filter( full_text )
        filtered_text = filtered_text + self._start_line
        return filtered_text[ : ]






    def _get_len_expected( self, args, keywds ):

        if( len( args) > 0 ):
            len_expected = args[ 0 ]
            if( len_expected < 0 ):
                len_expected = None
        elif( keywds.has_key( 'size' ) ):
            len_expected = keywds[ 'size' ]
        else:
            len_expected = None
        return len_expected

    def _filter( self, text, filter_chain = [ remove_empty_line, remove_useless_dot, fix_punctuation ] ):
        lines = text.splitlines( 1 )
        filtered_text = ''
        all_but_last_line = lines[ :-1 ]
        for line in all_but_last_line:
            for filter in filter_chain:
                line = apply( filter, ( line, ) )
                print line
            filtered_text = filtered_text + line
        last_line = lines[ -1 ]
        if( has_trailing_linefeed( last_line ) ):
            for filter in filter_chain:
                line = apply( filter, ( line, ) )
            filtered_text = filtered_text + line
        else:
            self._start_line = last_line
        return filtered_text

def has_trailing_linefeed( line ):
    if( line.endswith( chr( 13 ) ) or \
        line.endswith( chr( 10 ) ) ):
        return 1
    else:
        return 0