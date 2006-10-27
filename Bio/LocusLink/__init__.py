import string
import operator
from Bio import File
import Martel
from Martel.Dispatch import Dispatcher
from Martel import RecordReader
from mx import TextTools
from locus_format import locus_record

"""Parser for NCBI's LocusLink, curated sequence and descriptive information 
about genetic loci.

The LocusLink site is:
http://www.ncbi.nlm.nih.gov/LocusLink/
"""

class Record( dict):

    def __init__( self ):
        dict.__init__( self )

    def __str__( self ):
        queue_keys = self.keys()
        queue_keys.sort()
        out = ''
        for key in queue_keys:
            out = out +  '%s:\n' % key
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
            out = out + '\n'
        elif( isinstance( item, dict ) ):
            keys = item.keys()
            keys.sort()
            for subitem in keys:
                out = out + '%s %s:\n' % ( indent, subitem )
                out = out + self.print_item( item[ subitem ], level + 1 )
            out = out + '\n'
        elif( type( item ) == type( {} ) ):
            keys = item.keys()
            keys.sort()
            for subitem in keys:
                out = out + '%s %s:\n' % ( indent, subitem )
                out = out + self.print_item( item[ subitem ], level + 1 )
            out = out + '\n'
        else:
            out = out + '%s\n' % str( item )
        return out

class Iterator:
    """Iterator interface to move over a file of LocusLink entries one at a time.

    """
    def __init__(self, handle, parser = None):
        """Initialize the iterator.

        Arguments:
        o handle - A handle with LocusLink entries to iterate through.
        o parser - An optional parser to pass the entries through before
        returning them. If None, then the raw entry will be returned.
        """
        self.handle = File.UndoHandle( handle )
        self._reader = RecordReader.StartsWith( self.handle, '>>'  )
        self._parser = parser

    def next(self):
        """Return the next LocusLink record from the handle.

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
        self.interest_tags = [ "locus_line", "accnum_block", "phenotype_block", "db_block" ]

        # make a parser that returns only the tags we are interested in
        expression = Martel.select_names( locus_format.locus_record, self.interest_tags)
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
    """Create a LocusLink Record object from scanner generated information.
    """
    def __init__(self):
        Dispatcher.__init__( self )

    def startDocument( self ):
        self.data = Record()        


    def set_interest_tags( self, interest_tags ):
        self.interest_tags = interest_tags

    def start_locus_line( self, line, attrs ):
        self.save_characters()

    def end_locus_line( self, locus_record ):
        line = self.get_characters()
        cols = line.split( ':', 1 )

        key = cols[ 0 ]
        key = key.strip()
        newval = cols[ 1 ]
        newval = newval.strip()
        if key == 'BUTTON':
            pass
        elif not self.data.has_key( key ):
            self.data[ key ] = newval
        else:
 
            val = self.data[ key ]
            if( type( val ) == type( '' ) ):
                self.data[ key ] = [ val, newval ]
            elif( type( val ) == type( [] ) ):
                val.append( newval )
                self.data[ key ] = val

    def start_accnum_block( self, line, attrs ):
        self.save_characters()

    def end_accnum_block( self, locus_record ):
        block = self.get_characters()
        self.parse_block( block, 'ACCNUM' )
            
    def start_phenotype_block( self, line, attrs ):
        self.save_characters()

    def end_phenotype_block( self, locus_record ):
        block = self.get_characters()
        self.parse_block( block, 'PHENOTYPE' )

    def start_db_block( self, line, attrs ):
        self.save_characters()

    def end_db_block( self, locus_record ):
        block = self.get_characters()
        self.parse_block( block, 'DB' )

    def parse_block( self, block, block_key ):
        lines = block.splitlines()
        entry = {}
        for line in lines:
            cols = line.split( ':', 1 )

            key = cols[ 0 ]
            key = key.strip()
            newval = cols[ 1 ]
            newval = newval.strip()
            entry[ key ] =  newval

        if not self.data.has_key( block_key ):
            self.data[ block_key ] = [ entry, ]
        else:
 
            val = self.data[ block_key ]
            val.append( entry )
            self.data[ block_key ] = val

class RecordParser:
    """Parse LocusLink files into Record objects
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
