
# Copyright 2001 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""InterPro

This module provides code to work with html files from
http://www.ebi.ac.uk/interpro/


Classes:
Record             Holds interpro sequence data.
InterProParser     Parses interpro sequence data into a Record object.




"""

__all__ = [
    'Record',
    'InterProParser',
    ]

from types import *
import string
from Bio import File
from Bio import Index
import urllib
import sgmllib
from Bio.ParserSupport import *
from Bio.SeqFeature import Reference
import UserDict


class Record( UserDict.UserDict ):

    def __str__( self ):
        keys = self.keys()
        keys.sort()
        out = ''
        for key in keys:
            val = self[ key ]
            if( key == 'References' ):
                out = out + '\n%s\n' % key
                for reference in val:
                    out = out + '%s\n' % str( reference )
                out = out + '\n'
            elif( key == 'Examples' ):
                out = out + '\n%s\n' % key
                for example in val:
                    out = out + '%s\n' % example
            elif( key == 'Abstract' ):
                out = out + '\n%s\n' % key
                out = out + '%s...\n' % val[ : 80 ]
            elif( type( self[ key ] ) == type( [] ) ):
                out = out + '\n%s\n' % key
                for item in val:
                    out = out + '%s\n' % item

            else:
                out = out + '%s: %s\n' % ( key, self[ key ] )
        return out



class InterProParser(  sgmllib.SGMLParser ):
    """Parses InterPro sequence data into a Record object.

    """
    def reset(self):
        sgmllib.SGMLParser.reset( self )
        self.text = ''
        self.inter_pro_dict = Record()
        self.inter_pro_dict[ 'Database' ] = ''
        self.inter_pro_dict[ 'Accession' ] = ''
        self.inter_pro_dict[ 'Name' ] = ''
        self.inter_pro_dict[ 'Dates' ] = ''
        self.inter_pro_dict[ 'Type' ] = ''
        self.inter_pro_dict[ 'Parent' ] = ''
        self.inter_pro_dict[ 'Process' ] = ''
        self.inter_pro_dict[ 'Function' ] = ''
        self.inter_pro_dict[ 'Component' ] = ''
        self.inter_pro_dict[ 'Signatures' ] = []
        self.inter_pro_dict[ 'Abstract' ] = ''
        self.inter_pro_dict[ 'Examples' ] = []
        self.inter_pro_dict[ 'References' ] = []
        self.inter_pro_dict[ 'Database links' ] = []
        self._state = 'title'
        self._reference_state = ''
        self._key_waiting = ''
        self._current_reference = ''

    def parse(self, handle):
        self.reset()
        self.feed(handle)
        return self.inter_pro_dict

    def feed(self, handle):
        """feed(self, handle )

        Feed in interpro data for scanning.  handle is a file-like object
        containing interpro data.  consumer is a Consumer object that will
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
            line = string.strip( line )
            if( line[ -7: ] == '</HTML>' ):
                break
            text = text + ' ' + line

        sgmllib.SGMLParser.feed( self, text )


    def handle_data(self, newtext ):
        newtext = string.strip( newtext )
        self.text = self.text + newtext

    def start_table( self, attrs ):
        dict = pairlist_to_dict( attrs )
        for key in dict.keys():
            val = dict[ key ]

    def start_h2( self, attrs ):
        pass

    def end_h2( self ):
        self._state = 'chugging_along'

    def start_td( self, attrs ):
        dict = pairlist_to_dict( attrs )
        if( self._state == 'chugging_along' ):
            if( dict.has_key( 'class' ) ):
                if( dict[ 'class' ] == 'tag' ):
                    self._state = 'waiting_tag'
                    self._flush_text()
                elif( dict[ 'class' ] == 'inf' ):
                    self._state = 'waiting_inf'
                    self._flush_text()

    def end_td( self ):
        if( self._state == 'waiting_tag' ):
            self._key_waiting = self._flush_text()
            self._state = 'chugging_along'
        elif( self._state == 'waiting_inf' ):
            key = self._key_waiting
            if( self.inter_pro_dict.has_key( key ) ):
                val = self._flush_text()
                if( key == 'Signatures' ):
                    pass
                elif( key == 'Database links' ):
                    pass
                else:
                    self.inter_pro_dict[ key ] = val
            self._key_waiting = ''
            self._state = 'chugging_along'


    def start_ul( self, attrs ):
        if( self._key_waiting == 'Examples' ):
            self._state = 'examples'
            self._flush_text()

    def end_ul( self ):
        self._key_waiting = ''
        self._state = 'chugging_along'

    def start_ol( self, attrs ):
        if( self._key_waiting == 'References' ):
            self._state = 'references'
            self._reference_state = 'pubmed_id'
            self._flush_text()
            self._references = []

    def end_ol( self ):
        if( self._state == 'references' ):
            self._references.append( self._current_reference )
            self.inter_pro_dict[ 'References' ] = self._references
        self._state = 'chugging_along'

    def start_li( self, attrs ):
        if( self._state == 'references' ):
            self._reference_state = 'pubmed_id'
            self._flush_text()
            if( self._current_reference != '' ):
                self._references.append( self._current_reference )
            self._current_reference = Reference()

    def end_li( self ):
        if( self._state == 'examples' ):
            text = self._flush_text()
            ( self.inter_pro_dict[ 'Examples' ] ).append( text )

    def start_a( self, attrs ):
        dict = pairlist_to_dict( attrs )
        if( self._state == 'references' ):
            if( self._reference_state == 'pubmed_id' ):
                if( dict.has_key( 'name' ) ):
                    self._current_reference.pubmed_id = dict[ 'name' ]
                    self._reference_state = 'authors'
            elif( self._reference_state == 'journal' ):
                self._current_reference.journal = self._flush_text()
                self._reference_state = 'medline_id'

    def end_a( self ):
        if( self._state == 'references' ):
            if( self._reference_state == 'medline_id' ):
                text = self._flush_text()
                cols = text.split( ':' )
                try:
                    medline_id = cols[ 1 ]
                except IndexError:
                    medline_id = None
                else:
                    medline_id = medline_id[ : -1 ]
                self._current_reference.medline_id = medline_id

    def do_br( self, attrs ):
        if( self._state == 'references' ):
            if( self._reference_state == 'authors' ):
                self._current_reference.authors = self._flush_text()
                self._reference_state = 'title'
        elif( self._key_waiting == 'Signatures' ):
            self.inter_pro_dict[ 'Signatures' ].append( self._flush_text() )
        elif( self._key_waiting == 'Database links' ):
            self.inter_pro_dict[ 'Database links' ].append( self._flush_text() )

    def start_i( self, attrs ):
        pass

    def end_i( self ):
        if( self._state == 'references' ):
            if( self._reference_state == 'title' ):
                text = self._flush_text()
                self._current_reference.title = text
                self._reference_state = 'journal'


    def handle_starttag(self, tag, method, attrs):
        if( self._state == 'references' ):
            if( tag == 'li' ):
                self.stack.pop()
            elif( tag == 'a' ):
                if( self._reference_state == 'pubmed_id' ):
                    self.stack.pop()
        method(attrs)


    def _flush_text( self ):
        text = string.strip( self.text )
        self.text = ''
        return text[:]

def pairlist_to_dict( pairs ):
    dict = {}
    for pair in pairs:
        key = pair[ 0 ]
        val = pair[ 1 ]
        dict[key ] = val
    return dict



if( __name__ == '__main__' ):
    handle = open( 'IPR001064.htm')
    undo_handle = Bio.File.UndoHandle( handle )
    interpro_parser = InterProParser()
    record = interpro_parser.parse( handle )
    print str( record )
