# Copyright 2001 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Parser for SAF (Simple Alignment Format).

This is a fairly liberal multiple sequence alignment format, where
record names may contain up to 14 characters and no blanks.  Lines
beginging with a hash (#) are ignored.

SAF has been described as a simplified subset of MSF, dropping the
checksum and with more flexibility in terms of line length.

A current URL describing this file format is:
http://www.predictprotein.org/Dexa/optin_saf.html

This appears to replace the old URL of:
http://www.embl-heidelberg.de/predictprotein/Dexa/optin_safDes.html
"""

import warnings
warnings.warn("Bio.Saf has been deprecated, due to problems with Martel"\
              +" and recent versions of mxTextTools. If you want to"\
              +" continue to use this module (or read this file format),"\
              +" please get touch to avoid permanent removal of this"\
              +" module from Biopython.", DeprecationWarning)

# Martel
import Martel
from Martel import RecordReader
from Martel import Dispatch


from Bio import File
import saf_format
import Record


class Iterator:
    """Iterator interface to move over a file of Saf entries one at a time.
    """
    def __init__(self, handle, parser = None):
        """Initialize the iterator.

        Arguments:
        o handle - A handle with Saf entries to iterate through.
        o parser - An optional parser to pass the entries through before
        returning them. If None, then the raw entry will be returned.
        """
        self.handle = File.UndoHandle( handle )
        self._reader = RecordReader.Everything( self.handle )
        self._parser = parser

    def next(self):
        """Return the next Saf record from the handle.

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
        self.interest_tags = [ 'candidate_line', 'saf_record' ]

        # make a parser that returns only the tags we are interested in
        expression = Martel.select_names( saf_format.saf_record, self.interest_tags)
        self._parser = expression.make_parser(debug_level = debug)

    def feed(self, handle, consumer):
        """Feed a set of data into the scanner.

        Arguments:
        o handle - A handle with the information to parse.
        o consumer - The consumer that should be informed of events.
        """
        consumer.set_interest_tags( self.interest_tags )
        self._parser.setContentHandler( consumer )
#        self._parser.setErrorHandler(handle.ErrorHandler())

        self._parser.parseFile(handle)

class _RecordConsumer( Dispatch.Dispatcher ):
    """Create a Saf Record object from scanner generated information.
    """
    def __init__(self ):
        Dispatch.Dispatcher.__init__( self )
        self.data = Record.Record()
        self._refresh()

    def _refresh( self ):
        self._sequences = {}
        self._names = {}
        self._history = []
        self._guide = ''
        self._ref_length = 0
        self._ordinal = 0

    def set_interest_tags( self, interest_tags ):
        self.interest_tags = interest_tags

    def startDocument(self):
        self.data = Record.Record()
        self._refresh()


    def start_candidate_line(self, name, attrs):
        self.save_characters()

    def end_candidate_line(self, candidate_lines ):
        candidate_line = self.get_characters()
        name = candidate_line.split( ' ' )[ 0 ]
        sequence = candidate_line[ len( name ): ]
        name = name.strip()
        sequence = sequence.replace( " ", "" )
        if( self._guide == '' ):
            self._guide = name
            self._ref_length = len( sequence )
        elif( name == self._guide ):
            history = []
            self._ref_length = len( sequence )
        try:
            self._history.index( name )
        except ValueError:
            self._names[ self._ordinal ] = name
            self._ordinal = self._ordinal + 1
            self._history.append( name )
        sequence = sequence.strip()
        try:
            sequence = self._sequences[ name ] + sequence
        except KeyError:
            pass
        self._sequences[ name ] = sequence

    def start_saf_record( self, sequence, attrs ):
        self._sequences = {}

    def end_saf_record( self, saf_record ):
        ordinals = self._names.keys()
        ordinals.sort()
        for ordinal in ordinals:
            name = self._names[ ordinal ]
            sequence = self._sequences[ name ]
            self.data.alignment.add_sequence( name, sequence )
        self._refresh()

class RecordParser:
    """Parse Saf files into Record objects.
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
        """Parse the specified handle into a SAF record.
        """
        self._consumer = _RecordConsumer()
        self._scanner.feed(handle, self._consumer)
        return self._consumer.data

