import string
import operator
import urllib
import sgmllib
import UserDict
import Bio.File
import Martel
from mx import TextTools
import unigene_format



class UniGeneParser( sgmllib.SGMLParser ):

    def reset( self ):
        sgmllib.SGMLParser.reset( self )
        self.text = ''
        self.queue = UserDict.UserDict()
        self.taglist = []
        self.tag = 'html'
        self.nextkey = ''
        self.table = ''
        self.context = 'general_info'

    def parse( self, handle ):
        self.reset()
        self.feed( handle )
        for key in self.queue.keys():
            if( self.queue[ key ] == {} ):
                if( key[ :15 ] == 'UniGene Cluster' ):
                    self.queue[ 'UniGene Cluster' ] = key[ 16: ]
                del self.queue[ key ]
        return self.queue

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
            if( string.strip( line ) == '' ):
                break
            text = text + line
        sgmllib.SGMLParser.feed( self, text )



    def handle_data(self, newtext ):
        newtext = string.strip( newtext )
        self.text = self.text + newtext

    def start_a( self, attrs ):
        if( self.context == 'seq_info' ):
            self.text = ''

#        self.queue.append( attrs )

    def end_a( self ):
        if( self.context == 'seq_info' ):
            self.nextkey = self.text
            self.text = ''

    def start_b( self, attrs ):

        self.taglist.append( self.tag )
        self.tag = 'label'
        self.text = ''

    def end_b( self ):
        key = string.strip( self.text )
        self.text = ''

        try:
            self.tag = self.taglist.pop()
        except:
            self.tag = 'html'
        if( self.tag == 'table_data' ):
            if( self.context == 'general_info' ):
                self.nextkey = key
                self.text = ''
            elif( self.context == 'seq_info' ):
                if( key == 'Key to Symbols' ):
                    self.context = 'legend'
                    self.table = key
        elif( self.context == 'general_info' ):
            self.table = key
            if( string.find( key, 'SEQUENCE' ) != -1 ):
                self.context = 'seq_info'
            self.queue[ key ] = UserDict.UserDict()
        elif( self.context == 'seq_info' ):
            self.queue[ key ] = UserDict.UserDict()
            self.table = key



    def start_table( self, attrs ):
        self.taglist.append( self.tag )
        self.tag = 'table'

    def end_table( self ):
        try:
            self.tag = self.taglist.pop()
        except:
            self.tag = 'html'

    def start_tr( self, attrs ):
        self.taglist.append( self.tag )
        self.tag = 'table_row'
        self.text = ''

    def end_tr( self ):
        try:
            self.tag = self.taglist.pop()
        except:
            self.tag = 'html'
        text = self.text
        self.text = ''
        if( text[ 0 ] == ':' ):
            text = text[ 1: ]
        if( self.context == 'general_info' ):
            self.queue[ self.table ][ self.nextkey ] = text
        elif( self.context == 'seq_info' ):
            self.queue[ self.table ][ self.nextkey ] = text



    def start_td( self, attrs ):
        self.taglist.append( self.tag )
        self.tag = 'table_data'

    def end_td( self ):
        try:
            self.tag = self.taglist.pop()
        except:
            self.tag = 'html'
        if( self.context == 'seq_info' ):
            self.text = self.text + ' '

    def print_item( self, item ):
        if( type( item ) == type( '' ) ):
            if( item != '' ):
                print item
        elif( type( item ) == type( [] ) ):
            for subitem in item:
                self.print_item( subitem )
        elif( type( item ) == type( {} ) ):
            for subitem in item.keys():
                self.print_item( subitem )
                self.print_item( item[ subitem ] )
        else:
            print item

    def print_tags( self ):
        print '\nTAGS\n'
        for key in self.queue.keys():
            print key
            self.print_item( self.queue[ key ] )

    def join_tags( self ):
        self.data = '\n'.join( self.queue ) + '\n'
        print self.data


if( __name__ == '__main__' ):
    handle = urllib.urlopen( 'http://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Hs&CID=222015&OPT=text')
    undo_handle = Bio.File.UndoHandle( handle )
    unigene_parser = UniGeneParser()
    unigene_parser.parse( handle )
    unigene_parser.print_tags()
#    unigene_parser.print_data()


