# Copyright 2001 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Parser for the NBRF/PIR file format (DEPRECATED).

Please use Bio.SeqIO with the "pir" format instead."""

import warnings
warnings.warn("Bio.NBRF is deprecated." \
              + " We hope the new 'pir' format support in Bio.SeqIO will be" \
              + " suitable for most users.  Please get in touch on the " \
              + " mailing lists if this (or its removal) causes any problems "\
              + "for you.",
              DeprecationWarning)


# Martel
import Martel
from Martel import RecordReader

from Bio.ParserSupport import EventGenerator
from Bio import File
import nbrf_format
import Record

class Iterator:
    """Iterator interface to move over a file of NBRF entries one at a time.
    """
    def __init__(self, handle, parser = None):
        """Initialize the iterator.

        Arguments:
        o handle - A handle with NBRF entries to iterate through.
        o parser - An optional parser to pass the entries through before
        returning them. If None, then the raw entry will be returned.
        """
        self.handle = File.UndoHandle( handle )
        self._reader = RecordReader.StartsWith( self.handle, ">"  )
        self._parser = parser

    def next(self):
        """Return the next NBRF record from the handle.

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
        self.interest_tags = [ "sequence_type", \
            "sequence_name", \
            "comment", \
            "sequence_line", \
            "sequence_final_text" ]

        # make a parser that returns only the tags we are interested in
        expression = Martel.select_names( nbrf_format.nbrf_record, self.interest_tags)
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
    """Create an NBRF Record object from scanner generated information.
    """
    def __init__(self):
        self.data = Record.Record()
        self._sequences = []

    def sequence_type(self, sequence_type ):
        self.data.sequence_type = "".join(sequence_type)

    def sequence_name(self, sequence_name ):
        self.data.sequence_name = "".join(sequence_name)

    def comment(self, comment ):
        self.data.comment = "".join(comment)

    def sequence_line( self, sequences ):
        new_seq = "".join(sequences)
        parts = new_seq.split()
        self._sequences.append("".join(parts))

    def sequence_final_text( self, sequences ):
        new_seq = "".join(sequences)
        parts = new_seq.split()
        self._sequences.append("".join(parts))
        
        self.data.sequence.data = "".join(self._sequences)

class RecordParser:
    """Parse NBRF files into Record objects
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

