# Copyright 2001 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Code for more fancy file handles (DEPRECATED).

This module is now deprecated, and will be removed in a future release of
Biopython.

Classes:
Filtered is a decorator for File that allows the user to filter the output
on a line by line basis.

The FilteredReader module reads a file and applies a sequence of filters to the input
The constructor sets a default filter chain, but the user can select another filter by setting
Bio.FilteredReader.filter_chain.

handle = open( "filename" )
filtered_reader = Bio.FilteredReader( handle )
filtered_reader.filter_chain = [ remove_asterisks, replace_dot_with_dash ]
filtered_reasder.read()

All filters in the chain must provide the same interface with a line of text as the single
input parameter and altered text as the return value.
"""

import warnings
warnings.warn("Bio.FilteredReader is deprecated, and will be removed in a"\
              " future release of Biopython.  If you want to continue to use"\
              " this code, please get in contact with the developers"\
              " via the mailing lists to avoid its permanent removal from"\
              " Biopython.", DeprecationWarning)

def dump_saved( name, text, j ):
    """Used for debugging."""
    dump_file = open( name + '%d' % j, "w" )
    k = 0
    for i in range ( 0, len( text ), 80 ):
        dump_file.write(  '%s\n' % text[ i : i + 80 ] )
    dump_file.close()

def remove_leading_whitespace( line ):
    return line.lstrip()


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
        self.filter_chain = [ remove_empty_line, remove_useless_dot, fix_punctuation ]

    def __getattr__(self, attr):
        return getattr(self._handle, attr)



    def close(self, *args, **keywds ):
        return self._handle.close( *args, **keywds)

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
        len_adjusted -= len( self._start_line )
        filtered_text = ''
        while( len_filtered < len_expected ):

            text_read = self._handle.read( len_adjusted )
            full_text = self._start_line + text_read
            lines = full_text.splitlines( 1 )
            if( text_read == '' ):
                filtered_text = filtered_text + self.filter( lines )
                break
            else:
                all_but_last_line = lines[ :-1 ]
                self._start_line = lines[ -1 ]
                filtered_text = filtered_text + self.filter( all_but_last_line )
            len_filtered_text = len( filtered_text )
            len_adjusted = len_adjusted - len_filtered_text
        return filtered_text[ : ]

    def read_to_end( self ):
        filtered_text = ''
        text_read = self._handle.read()
        full_text = self._start_line + text_read
        lines = full_text.splitlines( 1 )
        filtered_text += self.filter( lines[:] )
        return filtered_text[ : ]

    def _get_len_expected( self, args, keywds ):

        if( len( args) > 0 ):
            len_expected = args[ 0 ]
            if( len_expected < 0 ):
                len_expected = None
        elif 'size' in keywds:
            len_expected = keywds['size']
        else:
            len_expected = None
        return len_expected

    def filter( self, lines  ):
        filter_chain = self.filter_chain
        filtered_text = ''
        for line in lines:
            for filter in filter_chain:
                line = filter( *( line, ) )
            filtered_text += line

        return filtered_text

def has_trailing_linefeed( line ):
    if( line.endswith( chr( 13 ) ) or \
        line.endswith( chr( 10 ) ) ):
        return 1
    else:
        return 0
