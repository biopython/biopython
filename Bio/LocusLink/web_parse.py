import sgmllib
import Bio.File

"""
The LocusLink site is:
http://www.ncbi.nlm.nih.gov/LocusLink/
Parses a Locus web page.
"""

import warnings
warnings.warn("Bio.LocusLink was deprecated, as NCBI's LocusLink was superceded by Entrez Gene. If you still need this module, please get in touch with the Biopython developers (biopython-dev@biopython.org) to avoid permanent removal of this module", DeprecationWarning)


def is_empty_container( item ):
    response = 0
    if is_container( item ):
        if len( item ) == 0:
            response = 1
    return response

def is_container( item ):
    response = 0
    if type( item ) in [ type( [] ), type( {} ) ]:
        response = 1
    return response

def is_substring( a, b ):
    if( a.find( b ) < 0 ):
        return 0
    else:
        return 1

def print_params( params ):
    print '%s!!!!!\n' % 'PARAMS'
    for item in params:

        print 'param ' + str( item )
    print '-----------------'

def process_list( params ):
    len_params = len( params )
    container = []
    while 1:
        try:
            element = params.pop()
        except:
            break
        if is_close_token( element ): break
        elif is_open_token( element ):
            break
        else:
            container.append( element )
    return container

def put( dict, key, val ):
    if dict.has_key( key ):
        element = dict[ key ]
        dict[ key ] = [ element, val ]
    else:
        dict[ key ] = val


def process_dict( params ):
    container = {}
    while len( params ) > 0:
        element = params.pop()
        if type( element ) == type( {} ):
            for key, val in element.items():
                put( container, key, val )
        elif is_close_token( element ): break
        elif is_open_token( element ):
            params.append( element )
        else:
            val = params.pop()
            if type( val ) == type( [] ):
                if len( val ) == 1:
                    val = val[ 0 ]
                try:
                    put( container, element, val )
                except:
                    print 'Element'
                    print element
                    params.append( element )

            elif(  not is_close_token( val ) ):
                try:
                    put( container, element, val )
                except:
                    print 'Element'
                    print element
                    params.append( element )
            else:
                break
    return container

class Token:
    def __init__( self, token  ):
        self.token = token

    def __eq__( self, other ):
        if not isinstance( other, self.__class__ ):
            return 0
        if self.token == other.token:
            return 1
        return 0

    def __ne__( self, other ):
        if not isinstance( other, Token ):
            return 1
        if self.token != other.token:
            return 1
        return 0

    def __str__( self ):
        output = 'token_%s\n' % self.token
        return output


open_list = Token( 'open_list' )
close_list = Token( 'close_list' )
open_dict = Token( 'open_dict' )
close_dict = Token( 'close_dict' )

def is_open_token( target ):
    answer = 0
    if isinstance( target, Token ):
        if ( open_list.__eq__( target ) ) or ( open_dict.__eq__(
target ) ):
            answer = 1
    return answer

def is_close_token( target ):
    answer = 0
    if isinstance( target, Token ):
        if ( close_list.__eq__( target ) ) or ( close_dict.__eq__(
target ) ):
            answer = 1
    return answer

def is_token( target ):
    return is_open_token( target ) or is_close_token( target )

class Url:

    def __init__( self, url, label = '', description = '' ):
        self.url = url
        self.label = label
        self.description = description

    def __str__( self ):
        output = '%s\n' % self.label
        output = output + 'url = %s\n' % self.url
        output = output + '%s\n' % self.description
        return output


class Record(dict):

    def __init__( self ):
        dict.__init__( self )

    def __str__( self ):
        queue_keys = self.keys()
        queue_keys.sort()
        out = ''
        for key in queue_keys:
            out = out +  '%s:\n' % key.upper()
            out = out + self.print_item( self[ key ] )
        out = out + '\n'

        return out

    def print_item( self, item, level = 1 ):
        indent = '    '
        out = ''
        for j in range( 0, level ):
            indent = indent + '    '
        if( type( item ) == type( '' ) ):
            if( item != '' ):
                out = out + '%s%s\n' % ( indent, item )
        elif( type( item ) == type([])):
            for subitem in item:
                out = out + self.print_item( subitem, level + 1 )
            out = out + '----------------------------------------------\n'
        elif( type( item ) == type ( {} ) ):
            keys = item.keys()
            keys.sort()
            for subitem in keys:
                out = out + '%skey is %s\n' % ( indent, subitem )
                out = out + self.print_item( item[ subitem ], level + 1 )
        elif( isinstance( item, dict ) ):
            keys = item.keys()
            keys.sort()
            for subitem in keys:
                out = out + '%skey is %s\n' % ( indent, subitem )
                out = out + self.print_item( item[ subitem ], level + 1 )
        else:
            out = out + '%s%s\n' % ( indent, str( item ) )
        return out


class LocusLinkParser( sgmllib.SGMLParser ):

    def reset( self ):
        sgmllib.SGMLParser.reset( self )
        self.text = ''
        self.record = Record()
        self.open_tag_stack = []
        self.open_tag = 'open_html'
        self.outer_state = 'undefined'
        self.section_state = 'undefined'
        self.local_title = ''
        self.structure_stack = []
        self.category = ''
        self.context_chain = []
        self.outer_state_dict = { 'nomenclature' : 'nomenclature', 'overview' : 'overview', \
            'function' : 'function', \
            'relationships' : 'relationships', \
            'locus' : 'locus', \
            'map' : 'map', \
            'refseq' : 'refseq', \
            'genbank' : 'genbank', \
            'external' : 'external_annotation', \
           'additional' : 'additional_links' \
        }


    def parse( self, handle ):
        self.reset()
        self.feed( handle )
        return self.record

#
# Assumes an empty line between records
#
    def feed( self, handle ):
        if isinstance(handle, Bio.File.UndoHandle):
            uhandle = handle
        else:
            uhandle = Bio.File.UndoHandle(handle)
        text = ''
        while 1:
            line = uhandle.readline()
            if not line:
                break
            text = text + ' ' + line

        sgmllib.SGMLParser.feed( self, text )

    def get_text( self ):
        text = self.text
        self.text = ''
        return text

    def handle_comment( self, comment ):
        while comment.startswith( '-' ):
            comment = comment[ 1: ]
        comment = comment.strip()
        comment = comment.lower()

        keys = self.outer_state_dict.keys()
        for key in keys:
            if comment.startswith( key ):
                if key in [ 'nomenclature', 'overview', 'function',
'relationships', 'map', 'locus', 'external' ]:
                    self.structure_stack.append( open_dict )
                elif key in [ 'genbank', 'additional' ]:
                    self.structure_stack.append( open_list )
                elif key in [ 'refseq' ]:
                    self.structure_stack.append( open_list )
                self.outer_state = key
                self.section_state = 'local_title'
                self.detail_state = 'undefined'
                if( key == 'refseq' ):
                    self.detail_state = 'waiting_category'
                else:
                    self.detail_state = 'waiting_key'
                break
        if comment.startswith( 'end' ):
            if is_substring( comment.lower(), self.outer_state ):
                if self.outer_state == 'refseq':
                    self.structure_stack.append( close_list )
                elif self.outer_state == 'function':
                    self.structure_stack.append( close_list )
                    self.structure_stack.append( close_dict )
                self.process_structure_stack()
                while 1:
                    try:
                        item = self.structure_stack.pop()
                    except:
                        item = 'Not Available'
                    if not is_token( item ) : break
                key = self.outer_state
                self.record[ self.outer_state_dict[ key ] ] = item
                self.outer_state = 'undefined'


    def handle_data(self, newtext ):
        newtext = newtext.strip()
        self.text = self.text + newtext

    def start_a( self, attrs ):
        self.open_tag_stack.append( self.open_tag )
        self.open_tag = 'open_a'
        attr_dict = {}
        for key, val in attrs:
            attr_dict[ key ] = val
        outer_state = self.outer_state
        if( outer_state in [ 'nomenclature', 'overview', 'relationships', 'locus', 'map', 'genbank', 'refseq', 'additional', 'external' ] ):
            if self.section_state == 'local_contents':
                if self.detail_state in [ 'scan_val', 'unpaired_key' ]:
                    if attr_dict.has_key( 'href' ):
                        href = attr_dict[ 'href' ]
                        self.text = ''
                        self.structure_stack.append( Url( href, '' ) )
        elif outer_state == 'function':
            if self.section_state == 'local_contents':
                if self.detail_state in [ 'scan_val', 'unpaired_key', 'may_be_val' ]:
                    if attr_dict.has_key( 'href' ):
                        href = attr_dict[ 'href' ]
                        self.text = ''
                        self.structure_stack.append( Url( href, '' ) )


    def end_a( self ):
        try:
            self.open_tag = self.open_tag_stack.pop()
        except:
            self.open_tag = 'open_html'
        outer_state = self.outer_state
        if( outer_state in [ 'nomenclature', 'overview', 'relationships', 'locus', 'map', 'refseq', 'genbank', 'additional', 'external' ] ):
            if self.section_state == 'local_contents':
                if self.detail_state in [ 'scan_val', 'unpaired_key' ]:
                    text = self.get_text()
                    url = self.structure_stack.pop()
                    if isinstance( url, Url ):
                        url.label = text
                    self.structure_stack.append( url )

        elif outer_state == 'function':
            if self.section_state == 'local_contents':
                if self.detail_state in [ 'scan_val', 'unpaired_key',
'may_be_val' ]:
                    text = self.get_text()
                    url = self.structure_stack.pop()
                    if isinstance( url, Url ):
                        url.label = text
                    self.structure_stack.append( url )

    def start_b( self, attrs ):

        self.open_tag_stack.append( self.open_tag )
        self.open_tag = 'open_b'
        outer_state = self.outer_state
        if( outer_state in [ 'nomenclature', 'overview', 'function', 'relationships', 'locus', 'map', 'refseq', 'genbank', 'additional', 'external' ] ):
            self.text = ''



    def end_b( self ):
        try:
            self.open_tag = self.open_tag_stack.pop()
        except:
            self.open_tag = 'open_html'
        outer_state = self.outer_state
        if( outer_state in [ 'nomenclature', 'overview', 'function', 'relationships', 'locus', 'map', 'refseq', 'genbank', 'additional', 'external' ] ):
            if self.section_state == 'local_contents':
                text = self.get_text()
                cols = text.split( ':', 1 )
                key = cols[ 0 ]
                if( outer_state == 'refseq' ):
                    self.structure_stack.append( cols[ 1 ] )
                    self.structure_stack.append( open_dict )
                    self.detail_state = 'waiting_key'
                elif outer_state == 'relationships':
                    self.structure_stack.append( key )
                    self.structure_stack.append( open_list )
                    self.detail_state = 'skip'
                elif outer_state == 'additional':
                    self.structure_stack.append( open_dict )
                    self.structure_stack.append( key )
                    self.structure_stack.append( open_list )
                    self.detail_state = 'unpaired_key'
                elif outer_state == 'function':
                    if self.detail_state != 'waiting_key':
                        self.structure_stack.append( close_list )
                    self.structure_stack.append( key )
                    self.detail_state = 'unpaired_key'
                    self.structure_stack.append( open_list )
                    self.structure_stack.append( open_list )
                    try:
                        val = cols[ 1 ]
                        if val.strip() != '':
                            self.structure_stack.append( val )
                            self.detail_state = 'unpaired_key'

                    except IndexError:
                        pass
                else:
                    if self.detail_state != 'waiting_key':
                        self.structure_stack.append( close_list )
                    self.detail_state = 'scan_val'
                    self.structure_stack.append( key )
                    self.structure_stack.append( open_list )
                    self.structure_stack.append( open_list )
                    try:
                        val = cols[ 1 ]
                        if val.strip() != '':
                            self.structure_stack.append( val )
                    except IndexError:
                        pass


    def start_th( self, attrs ):

        self.open_tag_stack.append( self.open_tag )
        self.open_tag = 'open_th'
        outer_state = self.outer_state
        self.text = ''
        if outer_state in [ 'function', 'relationships', 'map', 'locus', 'genbank', 'additional', 'external' ]:
            if self.section_state == 'local_contents':
                self.detail_state = 'scan_headings'



    def end_th( self ):
        try:
            self.open_tag = self.open_tag_stack.pop()
        except:
            self.open_tag = 'open_html'
        outer_state = self.outer_state
        if outer_state == 'refseq':
            if self.section_state == 'local_contents':
                text = self.get_text()
                cols = text.strip().split( ':', 1 )
                if text.strip().lower().startswith( 'category' ):
                    self.structure_stack.append( open_dict )
                    self.structure_stack.append( cols[ 1 ] )
                    self.structure_stack.append( open_list )
                    self.structure_stack.append( open_dict )
                    self.detail_state = 'found_category'

                elif self.detail_state in [ 'found_category', 'may_be_val' ]:
                    if  text.strip() != '':
                        if self.detail_state != 'found_category':
                            self.structure_stack.append( close_list )
                        cols  = text.split( ':' )
                        self.structure_stack.append( cols[ 0 ] )
                        self.structure_stack.append( open_list )
                        try:
                            val = cols[ 1 ]
                            self.structure_stack.append(  open_list )
                            self.structure_stack.append( val )
                            self.detail_state = 'scan_val'
                        except IndexError:
                            self.detail_state = 'may_be_val'





    def start_table( self, attrs ):
        self.open_tag_stack.append( self.open_tag )
        self.open_tag = 'open_table'
        self.text = ''
        if self.outer_state == 'genbank':
            if self.section_state == 'local_contents':
                self.detail_state = 'skip'
        elif( self.outer_state in [ 'nomenclature', 'overview', 'relationships', 'locus', 'map', 'genbank', 'additional', 'external' ] ):

            if self.section_state == 'local_contents':
                self.detail_state = 'waiting_key'

    def end_table( self ):
        try:
            self.open_tag = self.open_tag_stack.pop()
        except:
            self.open_tag = 'open_html'
        if( self.section_state == 'local_title' ):
            if self.outer_state == 'refseq':
                self.section_state = 'local_contents'
            elif self.outer_state == 'additional':
                self.section_state = 'local_contents'
                self.detail_state = 'scan_val'
            else:
                self.section_state = 'local_contents'
                self.detail_state = 'waiting_key'
        elif self.section_state == 'local_contents':
            if( self.outer_state in [  'nomenclature', 'relationships', 'locus', 'map', 'external' ] ):
                self.structure_stack.append( close_list )
            elif ( self.outer_state in [ 'genbank', 'additional' ] ):
                if self.detail_state == 'scan_val':
                    self.structure_stack.append( close_list )

            elif self.outer_state == 'refseq':
                if self.detail_state in ['may_be_val', 'scan_val' ]:
                    self.structure_stack.append( close_list )
                    self.structure_stack.append( close_dict )
                    self.structure_stack.append( close_list )
                    self.structure_stack.append( close_dict )
                    self.detail_state = 'scan_category'


    def start_tr( self, attrs ):
        top = self.open_tag
        self.open_tag_stack.append( self.open_tag )
        if top == 'open_table_row':
            if self.outer_state == 'refseq':
                if self.section_state == 'local_contents':
                    if self.detail_state in [ 'scan_val', ]:
                        self.structure_stack.append( close_list )
                        self.detail_state = 'may_be_val'
                        self.open_tag_stack.pop()
        self.open_tag = 'open_table_row'
        self.text = ''
        outer_state = self.outer_state
        if( outer_state in [   'relationships', 'locus', 'function', 'genbank', 'external'
] ):
            if self.section_state == 'local_contents':
                if self.detail_state == 'scan_val':
                    self.structure_stack.append( open_list )
        elif outer_state == 'map':
            if self.section_state == 'local_contents':
                if self.detail_state == 'scan_val':
                    self.structure_stack.append( open_list )

        elif outer_state == 'additional':
            if self.section_state == 'local_contents':
                self.detail_state = 'scan_val'
                self.structure_stack.append( open_list )


    def end_tr( self ):
        try:
            self.open_tag = self.open_tag_stack.pop()
        except:
            self.open_tag = 'open_html'
        if self.section_state == 'local_contents':
            if( self.outer_state in [  'overview', 'nomenclature', 'relationships',
'locus', 'genbank', 'external' ] ):
                if self.detail_state == 'scan_val':
                    self.structure_stack.append( close_list )
                elif self.detail_state == 'unpaired_key':
                    self.structure_stack.append( close_list )
                elif self.detail_state == 'skip':
                    self.detail_state = 'scan_val'
                elif self.detail_state == 'scan_headings':
                    self.detail_state = 'scan_val'
            elif self.outer_state in [ 'additional', ]:
                if self.detail_state == 'unpaired_key':
                    self.structure_stack.append( close_list )
                    self.structure_stack.append( close_dict )
                    self.structure_stack.append( close_list )
                elif self.detail_state == 'scan_val':
                    self.structure_stack.append( close_list )
            elif self.outer_state in [ 'function', ]:
                if self.detail_state == 'scan_headings':
                    self.detail_state = 'scan_val'
                elif self.detail_state == 'unpaired_key':
                    self.detail_state = 'may_be_val'
                    self.structure_stack.append( close_list )
                elif self.detail_state == 'scan_val':
                    self.detail_state = 'may_be_val'
                    self.structure_stack.append( close_list )
            elif self.outer_state in [ 'refseq', ]:
                if self.section_state == 'local_contents':
                    if self.detail_state == 'scan_val':
                        self.structure_stack.append( close_list )
                        self.detail_state = 'may_be_val'
            elif self.outer_state == 'map':
                if self.section_state == 'local_contents':
                    if self.detail_state == 'scan_val':
                        self.structure_stack.append( close_list )
                    self.detail_state = 'may_be_val'


    def start_td( self, attrs ):
        self.open_tag_stack.append( self.open_tag )
        self.open_tag = 'open_table_data'
        if self.outer_state in [ 'nomenclature', 'overview', 'relationships', 'map', 'locus', 'genbank', 'additional', 'external' ]:
            if( self.section_state == 'local_contents' ):
                self.text = ''
        elif self.outer_state == 'refseq':
            if self.section_state == 'local_contents':
                self.text = ''
                if self.detail_state == 'may_be_val':
                    self.structure_stack.append( open_list )
                    self.detail_state = 'scan_val'

    def end_td( self ):
        try:
            self.open_tag = self.open_tag_stack.pop()
        except:
            self.open_tag = 'open_html'



        if self.outer_state in [ 'nomenclature', 'overview',  'relationships', 'locus', 'genbank', 'additional', 'external' ]:
            if( self.section_state == 'local_contents' ):
                if self.detail_state == 'scan_val':
                    text = self.get_text()
                    if( text != '' ):
                        self.structure_stack.append( text )
        elif self.outer_state == 'function':
            if self.section_state == 'local_contents':
                text = self.get_text()
                if( text != '' ):
                    if self.detail_state == 'may_be_val':
                        if text.strip() != '':
                            self.structure_stack.append( open_list )
                            self.detail_state = 'scan_val'
                    if self.detail_state in [ 'unpaired_key', 'scan_val' ]:
                        self.structure_stack.append( text )
        elif self.outer_state == 'map':
            if self.section_state == 'local_contents':
                text = self.get_text()
                if( text != '' ):
                    if self.detail_state == 'may_be_val':
                        if text.strip() != '':
                            self.structure_stack.append( open_list )
                            self.detail_state = 'scan_val'
                    if self.detail_state == 'scan_val':
                        self.structure_stack.append( text )
        elif self.outer_state == 'refseq':
            if self.section_state == 'local_contents':
                if self.detail_state == 'scan_val':
                    text = self.get_text()
                    if( text != '' ):
                        self.add_text_to_object( text )

    def do_br( self, attrs ):
        if self.outer_state in [ 'nomenclature', 'overview', 'function', 'relationships', 'map', 'locus', 'genbank', 'additional', 'external' ]:
            if( self.section_state == 'local_contents' ):
                if self.detail_state == 'scan_val':
                    if self.is_contained_by( 'open_table_data' ):
                        text = self.get_text()
                        if( text != '' ):
                            self.structure_stack.append( text )


    def add_text_to_object( self, text ):
        stack_item = self.structure_stack.pop()
        if isinstance( stack_item, Url ):
            if stack_item.description == '':
                stack_item.description = text
            self.structure_stack.append( stack_item )
        else:
            self.structure_stack.append( stack_item )
            self.structure_stack.append( text )



    def is_contained_by( self, tag ):
        return tag in self.open_tag_stack

    def process_structure_stack( self ):
        params = []
        outer_state = self.outer_state
        if outer_state in [ 'nomenclature', 'overview', 'function', 'relationships', 'refseq', 'locus', 'map', 'genbank', 'additional', 'external' ]:
            while len( self.structure_stack ) > 1:
                len_stack = len( self.structure_stack )
#                self.print_stack()
                for i in range ( 0, len_stack ):
                    item = self.structure_stack.pop()
                    if not is_open_token( item ):
                        params.append( item )
                    else: break
                if( open_list.__eq__( item ) ):
                    container = process_list( params )
                    params.append( container )
                else:
                    container = process_dict( params )
                    if len( container ) > 0:
                        params.append( container )
                if ( len( self.structure_stack ) == 0  ) or is_open_token(
self.structure_stack[ -1 ] ):
                    for j in range( 0, len( params ) ):
                        item = params.pop()
                        self.structure_stack.append( item )
                    params = []


    def print_stack( self ):
        print '%s!!!!!\n' % self.outer_state.upper()
        for stack_item in self.structure_stack:
            print 'stack has ' + str( stack_item )
        print '-----------------'




if( __name__ == '__main__' ):
    handle = open( 'Hs13225.htm')
    undo_handle = Bio.File.UndoHandle( handle )
    locuslink_parser = LocusLinkParser()
    record = locuslink_parser.parse( handle )
    print record
