# Copyright 2000-2001, Dalke Scientific Software, LLC
# Distributed under the Biopython License Agreement (see the LICENSE file).

"""Iterate over records of a XML parse tree

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
>>> from Martel.formats import swissprot38
>>> from xml.dom import pulldom
>>> iterator = swissprot38.format.make_iterator("swissprot38_record")
>>> text = open("sample.swissprot").read()
>>> for record in iterator.iterateString(text, pulldom.SAX2DOM):
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


import urllib
import Parser
try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO

class StoreEvents:
    def __init__(self):
        self.events = []

    def startDocument(self):
        pass
    def endDocument(self):
        pass

    def startElement(self, *args):
        #print "startElement", args
        self.events.append( ("startElement", args) )
    def characters(self, *args):
        self.events.append( ("characters", args) )
    def endElement(self, *args):
        self.events.append( ("endElement", args) )

    def error(self, *args):
        self.events.append( ("error", args) )
    def fatalError(self, *args):
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

    def iterateString(self, s, make_cont_handler):
        """create an iterator over a string"""
        events = StoreEvents()
        self.parser.setContentHandler(events)
        self.parser.setErrorHandler(events)
        self.parser.parseString(s)
        return Iterate(EventStream(events.events), self.tag, make_cont_handler)

    def iterateFile(self, fileobj, make_cont_handler):
        return self.iterateString(self, fileobj.read(), make_cont_handler)
        
    def iterate(self, systemId, doc_handler):
        return self.iterateFile(urllib.urlopen(systemId), make_cont_handler)
    
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

    def iterateString(self, s, make_cont_handler):
        return self.iterateFile(StringIO(s), make_cont_handler)

    def iterateFile(self, fileobj, make_cont_handler):
        record_reader = apply(self.make_reader,
                              (fileobj,) + self.reader_args)
        return Iterate(RecordEventStream(record_reader, self.record_parser),
                       self.marker_tag, make_cont_handler)

    def iterate(self, systemId, make_cont_handler):
        return self.iterateFile(urllib.urlopen(systemId), make_cont_handler)

class Iterate:
    def __init__(self, event_stream, tag, make_cont_handler):
        self.event_stream = event_stream
        self.events = None
        self.tag = tag
        self.make_cont_handler = make_cont_handler
        self._n = 0
        self.start_position = 0
        self.end_position = 0
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
                raise args[0]

            if name == "startElement" and args[0] == self.tag:
                self.start_position = self.current_position
                cont_handler = self.make_cont_handler()
                cont_handler.startDocument()
                while i < n:
                    name, args = events[i]
                    if name == "error":
                        # in theory this is recoverable, so scan forward
                        # until there's an endElement
                        exc = args[0]
                        while i < n:
                            name, args = events[i]
                            if name == "endElement" and args[0] == self.tag:
                                del self.events[:i+1]
                                raise exc
                            elif name == "characters":
                                self.current_position += len(args[0])
                            i = i + 1
                        # no end found, so not recoverable
                        self.events = None
                        raise exc
                    elif name == "fatalError":
                        # not recoverable
                        self.events = None
                        raise args[0]
                    else:
                        apply(getattr(cont_handler, name), args)
                        if name == "endElement" and args[0] == self.tag:
                            self.end_position = self.current_position
                            del self.events[:i+1]
                            cont_handler.endDocument()
                            self._n = self._n + 1
                            return cont_handler
                        elif name == "characters":
                            self.current_position += len(args[0])
                        i = i + 1
                # Got here without an endElement?  Not supposed to happen!
                raise AssertionError, "no endElement(%s) and no errors?" % \
                      repr(self.tag)
            else:
                if name == "characters":
                    self.current_position += len(args[0])
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

def test1():
    from Martel.formats import swissprot38
    from Martel.test import support
    iterator = swissprot38.format.make_iterator("swissprot38_record")
    infile = open("/home/dalke/ftps/swissprot/smaller_sprot38.dat")
    output = support.Dump()
    from xml.sax import saxlib
    #output = saxlib.HandlerBase()
    for record in iterator.iterateFile(infile, output):
        pass

def test2():
    record_format = Martel.Group("letter", Martel.Is("a"))
    format = Group("all", Rep(record_format))
    
def test():
    test1()
    test2()

if __name__ == "__main__":
    test2()
    
