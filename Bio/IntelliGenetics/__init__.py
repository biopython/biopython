# Copyright 2001 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Parser for the MASE/IntelliGenetics alignment file format (DEPRECATED).

Please use Bio.SeqIO with the "ig" format instead."""

import warnings
warnings.warn("Bio.IntelliGenetics is deprecated." \
              + " We hope the new 'ig' format support in Bio.SeqIO will be" \
              + " suitable for most users.  Please get in touch on the " \
              + " mailing lists if this (or its removal) causes any problems "\
              + "for you.",
              DeprecationWarning)



# Martel
import Martel
from Martel import RecordReader

from Bio.ParserSupport import EventGenerator
from Bio import File
import intelligenetics_format
import Record
class Iterator:
    """Iterator interface to move over a file of IntelliGenetics entries one at a time.
    """
    def __init__(self, handle, parser = None):
        """Initialize the iterator.

        Arguments:
        o handle - A handle with IntelliGenetics entries to iterate through.
        o parser - An optional parser to pass the entries through before
        returning them. If None, then the raw entry will be returned.
        """
        self.handle = File.UndoHandle( handle )
        self._reader = IntelliGeneticsReader( self.handle )
        self._parser = parser

    def next(self):
        """Return the next IntelliGenetics record from the handle.

        Will return None if we ran out of records.
        """
        data = self._reader.next()

        if self._parser is not None:
            if data:
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
        self.interest_tags = ["comment", "title_line", "sequence" ]

        # make a parser that returns only the tags we are interested in
        expression = Martel.select_names(intelligenetics_format.intelligenetics_record, self.interest_tags)
        self._parser = expression.make_parser(debug_level = debug)

    def feed(self, handle, consumer):
        """Feeed a set of data into the scanner.

        Arguments:
        o handle - A handle with the information to parse.
        o consumer - The consumer that should be informed of events.
        """
        self._parser.setContentHandler( EventGenerator(consumer,
                                                       self.interest_tags))
#        self._parser.setErrorHandler(handle.ErrorHandler())

        self._parser.parseFile(handle)

class _RecordConsumer:
    """Create an IntelliGenetics Record object from scanner generated information.
    """
    def __init__(self):
        self.data = Record.Record()


    def title_line(self, title):
        self.data.title = title[ 0 ]

    def comment(self, comments ):
        for comment in comments:
            self.data.comments.append( comment )

    def sequence( self, sequences ):
        for sequence in sequences:
            self.data.sequence.data = self.data.sequence.data + sequence.strip()

class RecordParser:
    """Parse IntelliGenetics files into Record objects
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
        """Parse the specified handle into a GenBank record.
        """
        self._consumer = _RecordConsumer()
        self._scanner.feed(handle, self._consumer)
        return self._consumer.data

class IntelliGeneticsReader( Martel.RecordReader.RecordReader ):

    def __init__( self, infile ):
        Martel.RecordReader.RecordReader.__init__( self, infile )

    def next( self ):
        infile = self.infile
        state = 'COMMENT_STATE'
        record = ''
        while( state != 'DONE' ):
            line = infile.readline()
            if( line == '' ):
                state = 'DONE'
                break
            if( line[ 0 ] == ';' ):
                if( state == 'SEQUENCE_STATE' ):
                    state = 'DONE'
                    infile.saveline( line )
                elif( state == 'COMMENT_STATE' ):
                    record = record + line
            else:
                if( state == 'COMMENT_STATE' ):
                    record = record + line
                    state = 'SEQUENCE_STATE'
                elif( state == 'SEQUENCE_STATE' ):
                    record = record + line
        return record

