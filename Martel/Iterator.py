# Copyright 2000-2001, Dalke Scientific Software, LLC
# Distributed under the Biopython License Agreement (see the LICENSE file).

"""Iterate over records of a XML parse tree.

The standard parser is callback based over all the elements of a file.
If the file contains records, many people would like to be able to
iterate over each record and only use the callback parser to analyze
the record.

If the expression is a 'ParseRecords', then the code to do this is
easy; use its make_reader to grab records and its record_expression to
parse them.  However, this isn't general enough.  The use of a
ParseRecords in the format definition should be strictly a
implementation decision for better memory use.  So there needs to be
an API which allows both full and record oriented parsers.

Here's an example use of the API:
>>> import sys
>>> import swissprot38  # one is in Martel/test/testformats
>>> from xml.dom import pulldom
>>> iterator = swissprot38.format.make_iterator("swissprot38_record")
>>> text = open("sample.swissprot").read()
>>> for record in iterator.iterateString(text, pulldom.SAX2DOM()):
..      print "Read a record with the following AC numbers:"
...     for acc in record.document.getElementsByTagName("ac_number"):
...         acc.writexml(sys.stdout)
...         sys.stdout.write("\n")
...


There are several parts to this API.  First is the 'Iterator

There are two parts to the API.  One is the EventStream.  This
contains a single method called "next()" which returns a list of SAX
events in the 2-ple (event_name, args).  It is called multiple times
to return successive event lists and returns None if no events are
available.

The other is the Iterator

Sean McGrath has a RAX parser (Record API for XML) which uses a
concept similar to this.
"""


import sys, urllib, traceback
from xml.sax import saxutils
import Parser
try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO

class StoreEvents:
    def __init__(self):
        self.events = []
        self.has_error = 0
        self.characters = lambda ch, append = self.events.append: \
                          append( ("characters", ch) )

    def startDocument(self):
        pass
    def endDocument(self):
        pass

    def startElement(self, *args):
        self.events.append( ("startElement", args) )
##    def characters(self, s):
##        # Note: This doesn't store the args as a tuple!
##        self.events.append( ("characters", s) )
    def endElement(self, *args):
        self.events.append( ("endElement", args) )

    def error(self, *args):
        self.has_error = 1
        self.events.append( ("error", args) )
    def fatalError(self, *args):
        self.has_error = 1
        self.events.append( ("fatalError", args) )

class EventStream:
    def __init__(self, event_list):
        self.events = event_list
    def next(self):
        if self.events:
            x = self.events
            self.events = None
            return x
        return None

class Iterator:
    def __init__(self, parser, tag):
        self.parser = parser
        self.tag = tag

    def iterateString(self, s, cont_handler = None):
        """create an iterator over a string"""
        events = StoreEvents()
        self.parser.setContentHandler(events)
        self.parser.setErrorHandler(events)
        self.parser.parseString(s)
        return Iterate(self, EventStream(events.events), self.tag,
                       cont_handler)

    def iterateFile(self, fileobj, cont_handler = None):
        return self.iterateString(fileobj.read(), cont_handler)
        
    def iterate(self, source, cont_handler = None):
        """parse using the URL or file handle"""
        source = saxutils.prepare_input_source(source)
        file = source.getCharacterStream() or source.getByteStream()
        return self.iterateFile(file, cont_handler)
    
class RecordEventStream:
    def __init__(self, reader, parser):
        self.reader = reader
        self.parser = parser
    def next(self):
        text = self.reader.next()
        if text is None:
            return None
        events = StoreEvents()
        self.parser.setContentHandler(events)
        self.parser.setErrorHandler(events)
        self.parser.parseString(text)
        return events.events

class IteratorRecords:
    def __init__(self, record_parser, make_reader, reader_args, marker_tag):
        self.record_parser = record_parser
        self.make_reader = make_reader
        self.reader_args = reader_args
        self.marker_tag = marker_tag

    def copy(self):
        return self  # This is (so far) an immutable object

    def iterateString(self, s, cont_handler = None):
        return self.iterateFile(StringIO(s), cont_handler)

    def iterateFile(self, fileobj, cont_handler = None):
        record_reader = self.make_reader(
                              *(fileobj,) + self.reader_args)
        return Iterate(self,
                       RecordEventStream(record_reader, self.record_parser),
                       self.marker_tag, cont_handler)

    def iterate(self, source, cont_handler = None):
        """parse using the URL or file handle"""
        source = saxutils.prepare_input_source(source)
        file = source.getCharacterStream() or source.getByteStream()
        return self.iterateFile(file, cont_handler)

def _get_next_text(reader):
    try:
        return reader.next(), None
    except (KeyboardInterrupt, SystemExit):
        raise
    except:
        # Something unusual happened (couldn't find a record?)
        # so call it a fatal error and stop
        outfile = StringIO()
        traceback.print_exc(file=outfile)
        exc = Parser.ParserRecordException(
            outfile.getvalue(), sys.exc_info()[1])
        events = [ ("fatalError", (exc,)) ]
        return None, events


class HeaderFooterEventStream:
    def __init__(self, fileobj, 
                 header_parser, make_header_reader, header_args,
                 record_parser, make_record_reader, record_args,
                 footer_parser, make_footer_reader, footer_args):
        self.fileobj = fileobj
        
        self.header_parser = header_parser
        self.make_header_reader = make_header_reader
        self.header_args = header_args
        
        self.record_parser = record_parser
        self.make_record_reader = make_record_reader
        self.record_args = record_args
        
        self.footer_parser = footer_parser
        self.make_footer_reader = make_footer_reader
        self.footer_args = footer_args

        self._state = "header"
        self._reader = None
        self._lookahead = ""

    def next(self):
        if self._state == "header":
            x = self._header_next()
            self._state = "record"
            if x is not None:
                return x
            
        if self._state == "record":
            x = self._record_next()
            if x is not None:
                return x
            self._state = "footer"

        if self._state == "footer":
            x = self._footer_next()
            self._state = "end"
            if x is not None:
                return x

        if self._state == "end":
            if self._lookahead:
                return [ ("fatalError",
                          (Parser.ParserIncompleteException(0),)) ]
            return None
        
        raise AssertionError("Should not get here")

    def _header_next(self):
        assert self._reader is None
        if self.header_parser is None:
            return None
        reader = self.make_header_reader(
                       *(self.fileobj,) + self.header_args,
                       **{"lookahead": self._lookahead})
        text, errors = _get_next_text(reader)
        self.fileobj, self._lookahead = reader.remainder()
        if text is None:
            return errors
        events = StoreEvents()
        self.header_parser.setContentHandler(events)
        self.header_parser.setErrorHandler(events)
        self.header_parser.parseString(text)
        return events.events
        
    def _record_next(self):
        if self._reader is None:
            assert self.record_parser is not None
            reader = self.make_record_reader(
                           *(self.fileobj,) + self.record_args,
                           **{"lookahead": self._lookahead})
            self._lookahead = None
            self._reader = reader
        else:
            reader = self._reader
        text, errors = _get_next_text(reader)
        if text is None:
            self.fileobj, self._lookahead = reader.remainder()
            self._reader = None
            return errors
        
        events = StoreEvents()
        self.record_parser.setContentHandler(events)
        self.record_parser.setErrorHandler(events)
        self.record_parser.parseString(text)

        if events.has_error:
            # Couldn't parse the record.
            if self.footer_parser is not None:
                # perhaps there's a footer here?
                # We'll need to try reading that
                self.fileobj, self._lookahead = reader.remainder()
                self._lookahead = text + self._lookahead
                self._reader = None
                return None
            # If no footer is possible, go on and pass
            # back the error as normal
            
        return events.events
            
    def _footer_next(self):
        assert self._reader is None
        if self.footer_parser is None:
            return None
        reader = self.make_footer_reader(
                       *(self.fileobj,) + self.footer_args,
                       **{"lookahead": self._lookahead})
        text, errors = _get_next_text(reader)
        self.fileobj, self._lookahead = reader.remainder()
        if text is None:
            return errors
        events = StoreEvents()
        self.footer_parser.setContentHandler(events)
        self.footer_parser.setErrorHandler(events)
        self.footer_parser.parseString(text)
        return events.events
        
            

class IteratorHeaderFooter:
    def __init__(self,
                 header_parser, make_header_reader, header_args,
                 record_parser, make_record_reader, record_args,
                 footer_parser, make_footer_reader, footer_args,
                 marker_tag):

        self.args = header_parser, make_header_reader, header_args, \
                    record_parser, make_record_reader, record_args, \
                    footer_parser, make_footer_reader, footer_args
        self.marker_tag = marker_tag

    def iterateString(self, s, cont_handler = None):
        return self.iterateFile(StringIO(s), cont_handler)
        
    def iterateFile(self, fileobj, cont_handler = None):
        args = (fileobj,) + self.args
        return Iterate(self, HeaderFooterEventStream(*args),
                       self.marker_tag, cont_handler)
    
    def iterate(self, source, cont_handler = None):
        """parse using the URL or file handle"""
        source = saxutils.prepare_input_source(source)
        file = source.getCharacterStream() or source.getByteStream()
        return self.iterateFile(file, cont_handler)
    

class Iterate:
    def __init__(self, parent, event_stream, tag, cont_handler = None):
        self.parent = parent
        if cont_handler is None:
            import LAX
            cont_handler = LAX.LAX()
        self.event_stream = event_stream
        self.events = None
        self.tag = tag
        self.cont_handler = cont_handler
        self._n = 0
        self.parent.start_position = 0
        self.parent.end_position = 0
        self.current_position = 0
        
    def next(self):
        events = self.events
        if not events:
            events = self.event_stream.next()
            if events is None:
                return None
            self.events = events

        i = 0
        n = len(events)
        # Look for the start of the next record
        while 1:
            if i == n:
                new_events = self.event_stream.next()
                if new_events is None:
                    break
                events.extend(new_events)
                n = len(events)

            name, args = events[i]
            if name == "error" or name == "fatalError":
                # at this level the error is unrecoverable
                self.events = None
                if isinstance(args[0], Parser.ParserPositionException):
                    exc = args[0]
                    exc.pos = 0
                    exc += self.current_position
                raise args[0]

            if name == "startElement" and args[0] == self.tag:
                self.parent.start_position = self.current_position
                cont_handler = self.cont_handler
                cont_handler.startDocument()
                while i < n:
                    name, args = events[i]
                    if name == "characters":
                        # This is the most common case.
                        # Recall, args is not a tuple
                        cont_handler.characters(args)
                        self.current_position += len(args)
                        i = i + 1
                    elif name == "error":
                        # in theory this is recoverable, so scan forward
                        # until there's an endElement
                        exc = args[0]
                        while i < n:
                            name, args = events[i]
                            if name == "endElement" and args[0] == self.tag:
                                del self.events[:i+1]
                                if isinstance(exc, Parser.ParserPositionException):
                                    exc.pos = 0
                                    exc += self.current_position
                                raise exc
                            elif name == "characters":
                                self.current_position += len(args)
                            i = i + 1
                        # no end found, so not recoverable
                        self.events = None
                        if isinstance(exc, Parser.ParserPositionException):
                            exc.pos = 0
                            exc += self.parent.start_position
                        raise exc
                    elif name == "fatalError":
                        # not recoverable
                        self.events = None
                        if isinstance(args[0], Parser.ParserPositionException):
                            exc = args[0]
                            exc = 0
                            exc += self.parent.start_position
                        raise args[0]
                    else:
                        getattr(cont_handler, name)(*args)
                        if name == "endElement" and args[0] == self.tag:
                            self.parent.end_position = self.current_position
                            del self.events[:i+1]
                            cont_handler.endDocument()
                            self._n = self._n + 1
                            return cont_handler
                        i = i + 1

                # Got here without an endElement?  Not supposed to happen!
                raise AssertionError, "no endElement(%s) and no errors?" % \
                      repr(self.tag)
            else:
                if name == "characters":
                    self.current_position += len(args)
                i = i + 1

        # Went through the document and no more records were found
        self.events = None
        return None

    def __getitem__(self, n):
        assert n == self._n, "forward iteration only"
        x = self.next()
        if x is None:
            raise IndexError, n
        return x

    def __iter__(self):
        return iter(self.next, None)
