# Copyright 2000-2001, Dalke Scientific Software, LLC
# Distributed under the Biopython License Agreement (see the LICENSE file).

"""Implement Martel parsers.

The classes in this module are used by other Martel modules and not
typically by external users.

There are two major parsers, 'Parser' and 'RecordParser.'  The first
is the standard one, which parses the file as one string in memory
then generates the SAX events.  The other reads a record at a time
using a RecordReader and generates events after each read.  The
generated event callbacks are identical.

At some level, both parsers use "_do_callback" to convert mxTextTools
tags into SAX events.

XXX finish this documentation

XXX need a better way to get closer to the likely error position when
parsing.

XXX need to implement Locator

"""
import urllib, traceback, sys
from xml.sax import handler, saxutils
import Parser, RecordReader

try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO


class IterRecords:
    def __init__(self, record_parser, make_reader, reader_args, marker_tag):
        self.record_parser = record_parser
        self.make_reader = make_reader
        self.reader_args = reader_args
        self.marker_tag = marker_tag

    def copy(self):
        return IterRecords(self.record_parser.copy(),
                           self.make_reader,
                           self.reader_args,
                           self.marker_tag)

    def iterate(self, source, cont_handler = None):
        source = saxutils.prepare_input_source(source)
        file = source.getCharacterStream() or source.getByteStream()
        return self.iterateFile(file, cont_handler)

    def iterateString(self, s, cont_handler = None):
        return self.iterateFile(StringIO(s), cont_handler)

    def iterateFile(self, fileobj, cont_handler = None):
        self.start_position = 0
        if cont_handler is None:
            import LAX
            cont_handler = LAX.LAX()
        self.record_parser.setContentHandler(cont_handler)

        reader = self.make_reader(fileobj, *self.reader_args)
        while 1:
            try:
                rec = reader.next()
            except RecordReader.ReaderError:
                raise Parser.ParserPositionException(self.start_position)
            if rec is None:
                break
            self.end_position = self.start_position + len(rec)
            try:
                self.record_parser.parseString(rec)
            except Parser.ParserPositionException, exc:
                exc += self.start_position
                raise
            
            yield cont_handler
            self.start_position = self.end_position

        fileobj, lookahead = reader.remainder()
        if lookahead or fileobj.read(1):
            raise Parser.ParserPositionException(self.start_position)
        

class IterHeaderFooter:
    def __init__(self,
                 header_parser, make_header_reader, header_args,
                 record_parser, make_record_reader, record_args,
                 footer_parser, make_footer_reader, footer_args,
                 marker_tag):
        self.header_parser = header_parser
        self.make_header_reader = make_header_reader
        self.header_args = header_args
        
        self.record_parser = record_parser
        self.make_record_reader = make_record_reader
        self.record_args = record_args
        
        self.footer_parser = footer_parser
        self.make_footer_reader = make_footer_reader
        self.footer_args = footer_args
        
        self.marker_tag = marker_tag

    def copy(self):
        header_parser = self.header_parser
        if header_parser is not None:
            header_parser = header_parser.copy()
        record_parser = self.record_parser.copy()
        footer_parser = self.footer_parser
        if footer_parser is not None:
            footer_parser = footer_parser.copy()
            
        return IterHeaderFooter(
            header_parser, self.make_header_reader, self.header_args,
            record_parser, self.make_record_reader, self.record_args,
            footer_parser, self.make_footer_reader, self.footer_args,
            self.marker_tag)

    def iterate(self, source, cont_handler = None):
        """parse using the URL or file handle"""
        source = saxutils.prepare_input_source(source)
        file = source.getCharacterStream() or source.getByteStream()
        return self.iterateFile(file, cont_handler)

    def iterateString(self, s, cont_handler = None):
        return self.iterateFile(StringIO(s), cont_handler)

    def iterateFile(self, fileobj, cont_handler = None):
        self.start_position = self.end_position = 0
        if cont_handler is None:
            import LAX
            cont_handler = LAX.LAX()
        self.record_parser.setContentHandler(cont_handler)

        lookahead = ""

        # By construction, we never need events from the header
        # nor from the footer
        if self.header_parser is not None:
            reader = self.make_header_reader(fileobj, *self.header_args,
                                             **{"lookahead": lookahead})
            try:
                rec = reader.next()
            except RecordReader.ReaderError:
                raise Parser.ParserPositionException(self.start_position)
            self.end_position = self.start_position + len(rec)
            self.header_parser.parseString(rec)
            self.start_position = self.end_position
            fileobj, lookahead = reader.remainder()

        reader = self.make_record_reader(fileobj, *self.record_args,
                                         **{"lookahead": lookahead})

        if not self.footer_parser:
            while 1:
                try:
                    rec = reader.next()
                except RecordReader.ReaderError:
                    raise Parser.ParserPositionException(self.start_position)
                if rec is None:
                    break
                self.end_position = self.start_position + len(rec)
                try:
                    self.record_parser.parseString(rec)
                except Parser.ParserPositionException, exc:
                    exc += self.start_position
                    raise
                yield cont_handler
                self.start_position = self.end_position
            return
            
        # This one is tedious
        while 1:
            try:
                rec = reader.next()
            except RecordReader.ReaderError:
                # we may have stumbled into the footer
                fileobj, lookahead = reader.remainder()
                break

            if not rec:
                # maybe there's a footer left
                fileobj, lookahead = reader.remainder()
                break

            try:
                self.record_parser.parseString(rec)
            except Parser.ParserException:
                # we may have tried to parse the footer
                fileobj, lookahead = reader.remainder()
                lookahead = rec + lookahead
                break
            self.end_position = self.start_position + len(rec)
            yield cont_handler
            self.start_position = self.end_position

        # Try to read the footer
        reader = self.make_footer_reader(fileobj, *self.footer_args,
                                         **{"lookahead": lookahead})
        try:
            rec = reader.next()
        except RecordReader.ReaderError:
            raise Parser.ParserPositionException(self.start_position)

        if rec is None:
            # Could read any footer
            raise Parser.ParserPositionException(self.start_position)

        try:
            self.footer_parser.parseString(rec)
        except Parser.ParserPositionException, exc:
            exc += self.start_position
            raise
        self.end_position = self.start_position + len(rec)
        self.start_position = self.end_position
        
        fileobj, lookahead = reader.remainder()
        if lookahead or fileobj.read(1):
            raise Parser.ParserIncompleteException(self.start_position)
