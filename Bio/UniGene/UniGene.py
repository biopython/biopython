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
            if( key == 'SEQUENCE INFORMATION' ):
                self.context = 'seq_info'
            else:
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
    text = ''
    unigene_parser = UniGeneParser()
    while 1:
        line = undo_handle.readline()
        if( line == '' ):
            break
        text = text + line
    unigene_parser.feed( text )
    unigene_parser.print_tags()
#    unigene_parser.print_data()


