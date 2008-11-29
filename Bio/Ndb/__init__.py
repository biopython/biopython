# Copyright 2002 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
This module provided code to parse HTML files from NDB (DEPRECATED).

This module provides an HTML parser designed for the NDB website
http://ndbserver.rutgers.edu/ as it was circa 2002.  The site has since
been redesigned, breaking the parser.  Bio.Ndb is therefore deprecated,
and will be removed in a future release of Biopython.

Classes:
Record             Holds NDB sequence data.
NdbParser          Parses NDB sequence data into a Record object.

The algorithm is based on a state machine because the record has multiple 
sections and the handling of tags varies depending on the section.  
Citations have their own state machine.
"""
import warnings
warnings.warn("Bio.Ndb has been deprecated as the NDB website it used to"\
              " parse has been redesigned.", DeprecationWarning)

from types import *
from Bio import File
from Bio import Index
from Bio.Crystal import Hetero
from Bio.Crystal import Chain
from Bio.Crystal import Crystal
from Bio.SeqFeature import Reference
import urllib
import sgmllib
from Bio.ParserSupport import *
from Bio.SeqFeature import Reference


class Record( dict ):

    def __init__( self ):
        self[ 'Id' ] = ''
        self[ 'Features' ] = ''
        self[ 'Name' ] = ''
        self[ 'Sequence' ] = Crystal( {} )
        self[ 'Citation' ] = Reference()
        self[ 'Space Group' ] = ''
        self[ 'Cell Constants' ] = {}
        self[ 'Crystallization Conditions' ] = []
        self[ 'Refinement' ] = ''
        self[ 'Coordinates' ] = ''

    def __str__( self ):
        keys = self.keys()
        keys.sort()
        out = ''
        for key in keys:
            val = self[ key ]
            if( type( val ) == type( [] ) ):
                out = out + '\n%s\n' % key
                for item in val:
                    out = out + '%s\n' % item

            elif( type( val ) == type( {} ) ):
                out = out + '\n%s\n' % key
                subkeys = val.keys()
                subkeys.sort()
                for item in subkeys:
                    out = out + '%s : %s\n' % ( item, val[ item ] )
            elif( isinstance( val, dict ) ):
                out = out + '\n%s\n' % key
                subkeys = val.keys()
                subkeys.sort()
                for item in subkeys:
                    out = out + '%s : %s\n' % ( item, val[ item ] )

            else:
                out = out + '%s: %s\n' % ( key, self[ key ] )
        return out

def _parse_constants( text ):
    items = text.split( '=' )
    constants = {}
    key = ''
    for i in range( 0, ( len( items ) - 1 ) ):
        item = items[ i ]
        item = item.strip()
        separator = item.rfind( ' ' )
        if( separator < 0 ):
            separator = 0
        val = item[ :separator ]
        val = val.strip()
        if( key != '' ):
            constants[ key ] = val
        key = item[ separator: ]
        key = key.strip()
    constants[ key ] = items[ -1 ]
    return constants





class NdbParser(  sgmllib.SGMLParser ):
    """Parses Ndb sequence data into a Record object.
    data available at: http://ndbserver.rutgers.edu/NDB/NDBATLAS/index.html
    """
    def reset(self):
        sgmllib.SGMLParser.reset( self )
        self.ndb_dict = Record()
        self.text = ''
        self._space_group = ''
        self._state = 'id'
        self._reference_state = 'authors'
        self._current_reference = Reference()

    def parse(self, handle):
        self.reset()
        self.feed(handle)
        return self.ndb_dict

    def feed(self, handle):
        """feed(self, handle )

        Feed in ndb data for scanning.  handle is a file-like object
        containing ndb data.  consumer is a Consumer object that will
        receive events as the ndb data is scanned.

        """
        if isinstance(handle, File.UndoHandle):
            uhandle = handle
        else:
            uhandle = File.UndoHandle(handle)
        text = ''
        while 1:
            line = uhandle.readline()
            if( not line ):
                break
            line = line.strip()
            if( line[ -7: ] == '</HTML>' ):
                break
            text = text + ' ' + line

        sgmllib.SGMLParser.feed( self, text )


    def handle_data(self, newtext ):
        newtext = newtext.strip()
        self.text = self.text + newtext

    def start_h1( self, attrs ):
        self._flush_text()

    def end_h1( self ):
        text = self._flush_text()
        if( self._state == 'id' ):
            cols = text.split( ':' )
            self.ndb_dict[ 'Id' ] = ( cols[ 1 ] ).upper()
            self._state = 'id_found'

    def start_h2( self, attrs ):
        text = self._flush_text()
        if( self._state == 'features' ):
            self.ndb_dict[ 'Features' ] = text
        elif( self._state == 'name' ):
            self.ndb_dict[ 'Name' ] = text
        elif( self._state == 'sequence' ):
            pass
        elif( self._state == 'citation' ):
            if( self._reference_state == 'journal' ):
                self._current_reference.journal = text
            self.ndb_dict[ 'Citation' ] = self._current_reference
        elif( self._state == 'space' ):
            self._space_group = self._space_group + text
            self.ndb_dict[ 'Space Group' ] = self._space_group
        elif( self._state == 'constants' ):
            self.ndb_dict[ 'Cell Constants' ] = _parse_constants( text )
        elif( self._state == 'crystallization' ):
            pass
        elif( self._state == 'refinement' ):
            self.ndb_dict[ 'Refinement' ] = text
        elif( self._state == 'coordinates' ):
            self.ndb_dict[ 'Coordinates' ] = text

    def end_h2( self ):
        text = self._flush_text()
        text = text.lower()
        if( self._state == 'id' ):
            if( text.find( 'id' ) >= 0 ):
                cols = text.split( ':' )
                self.ndb_dict[ 'Id' ] = ( cols[ 1 ] ).upper()
                self._state = 'id_found'
        elif( text.find( 'feature' ) >= 0 ):
            self._state = 'features'
        elif( text.find( 'name' ) >= 0 ):
            self._state = 'name'
        elif( text.find( 'sequence' ) >= 0 ):
            self._state = 'sequence'
        elif( text.find( 'citation' ) >= 0 ):
            self._state = 'citation'
        elif( text.find( 'space' ) >= 0 ):
            self._state = 'space'
        elif( text.find( 'constants' ) >= 0 ):
            self._state = 'constants'
        elif( text.find( 'crystallization' ) >= 0 ):
            self._state = 'crystallization'
        elif( text.find( 'refinement' ) >= 0 ):
            self._state = 'refinement'
        elif( text.find( 'coordinates' ) >= 0 ):
            self._state = 'coordinates'


    def start_ul( self, attrs ):
        if( self._state == 'sequence' ):
            self._flush_text()

        elif( self._state == 'crystallization' ):
            self._flush_text()

    def end_ul( self ):
        if( self._state == 'sequence' ):
            self._parse_chain()
        elif( self._state == 'crystallization' ):
            text = self._flush_text()
            ( self.ndb_dict[ 'Crystallization Conditions' ] ).append( text )
        elif( self._state == 'citation' ):
            if( self._reference_state == 'journal' ):
                self._current_reference.journal = self._flush_text()
                self._reference_state = 'done'

    def start_sub( self, attrs ):
        if( self._state == 'space' ):
            self._space_group = self._space_group + self._flush_text()

    def end_sub( self ):
        if( self._state == 'space' ):
            self._space_group = self._space_group + '(%s) ' % self._flush_text()

    def start_li( self, attrs ):
        if( self._state == 'sequence' ):
            self._parse_chain()
        elif( self._state == 'crystallization' ):
            text = self._flush_text()
            ( self.ndb_dict[ 'Crystallization Conditions' ] ).append( text )

    def end_li( self ):
        if( self._state == 'sequence' ):
            self._parse_chain()
        elif( self._state == 'crystallization' ):
            text = self._flush_text()
            ( self.ndb_dict[ 'Crystallization Conditions' ] ).append( text )

    def do_br( self, attrs ):
        if( self._state == 'citation' ):
            if( self._reference_state == 'authors' ):
                self._current_reference.authors = self._flush_text()
                self._reference_state = 'title'
            elif( self._reference_state == 'title' ):
                self._current_reference.title = self._flush_text()
                self._reference_state = 'journal'

    def start_i( self, attrs ):
        pass

    def end_i( self ):
        if( self._state == 'references' ):
            if( self._reference_state == 'title' ):
                text = self._flush_text()
                self._current_reference.title = text
                self._reference_state = 'journal'


    def _parse_chain(  self ):
        text = self._flush_text()
        text = text.strip()
        if( text.lower().startswith( 'chain' ) ):
            fields = text.split( ':' )
            words = fields[ 0 ].split()
            key = words[ 1 ]
            val = fields[ 1 ]
            self.ndb_dict[ 'Sequence' ][ key ] = val



    def _flush_text( self ):
        text = self.text.strip()
        self.text = ''
        return text[:]


if( __name__ == '__main__' ):
    handle = open( 'PR0004.htm')
    undo_handle = File.UndoHandle( handle )
    ndb_parser = NdbParser()
    record = ndb_parser.parse( handle )
    print str( record )
