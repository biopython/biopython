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

import urllib, pprint, traceback, sys, string
from xml.sax import xmlreader, _exceptions, handler, saxutils

try:
    from mx import TextTools
except ImportError:
    import TextTools

try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO

import Dispatch

# These exceptions are liable to change in the future
class ParserException(_exceptions.SAXException):
    """used when a parse cannot be done"""
    def setLocation(self, text):
        self._msg += "; in %s" % repr(text)

class ParserPositionException(ParserException):
    def __init__(self, pos):
        ParserException.__init__(self,
                    "error parsing at or beyond character %d" % pos,
                                     None)
        self.pos = pos
    def __iadd__(self, offset):
        self.pos += offset
        self._msg = "error parsing at or beyond character %d" % self.pos
        return self

class ParserIncompleteException(ParserPositionException):
    def __init__(self, pos):
        ParserPositionException.__init__(self, pos)
        self._msg += " (unparsed text remains)"
    def __iadd__(self, offset):
        ParserPositionException.__iadd__(self, offset)
        self._msg += " (unparsed text remains)"

class ParserRecordException(ParserException):
    """used by the RecordParser when it can't read a record"""
    pass


# Uses a hack to support back references in mxTextTools!

# THIS MEANS SINGLE THREADED SUPPORT for anything using
# backreferences!  There is a much more complicated solution where the
# first element of any taglist is defined to contain the _match_group
# for that parse session.  I don't want to do that, since another
# approach is to modify mxTextTools to pass around an extra state
# object, or to write my own code.  (Did someone say NIH syndrome? :)
_match_group = {}



# The SAX startElements take an AttributeList as the second argument.
# Martel's attributes can be empty, so make a simple class which
# doesn't do anything and which I can guarantee won't be modified.
class MartelAttributeList(xmlreader.AttributesImpl):
    def getLength(self):
        return 0
    def getName(self, i):
        raise IndexError, i
    def getType(self, i):
        raise IndexError, i
    def getValue(self, i):
        raise IndexError, i
    def __len__(self):
        return 0
    def __getitem__(self, key):
        if type(key) == type(0):
            raise IndexError, key
        else:
            raise KeyError, key
    def keys(self):
        return []
    def values(self):
        return []
    def items(self):
        return []
    def has_key(self, key):
        return 0
    def get(self, key, alternative):
        return alternative
    def __repr__(self):
        return "{}"
    def __str__(self):
        return "{}"

# singleton object shared amoung all startElement calls
_attribute_list = MartelAttributeList([])


def _do_callback(s, begin, end, taglist, cont_handler, attrlookup):
    """internal function to convert the tagtable into ContentHandler events

    's' is the input text
    'begin' is the current position in the text
    'end' is 1 past the last position of the text allowed to be parsed
    'taglist' is the tag list from mxTextTools.parse
    'cont_handler' is the SAX ContentHandler
    'attrlookup' is a dict mapping the encoded tag name to the element info
    """
    # bind functions to local names for a slight speedup
    characters   = cont_handler.characters
    startElement = cont_handler.startElement
    endElement   = cont_handler.endElement

    for tag, l, r, subtags in taglist:
        # If the tag's beginning is after the current position, then
        # the text from here to the tag's beginning are characters()
        assert begin <= l, "begin = %d and l = %d" % (begin, l)
        if begin < l:
            characters(s[begin:l])

        if tag.startswith(">"):
            # Named groups doesn't create ">ignore" tags, so pass them on
            # to the ContentHandler.  Unnamed groups still need a name so
            # mxTextTools can create subtags for them.  I named them
            # ">ignore" - don't create events for them.
            if not tag == ">ignore":
                assert tag.startswith(">G"),"Unknown special tag %s" % repr(tag)
                # This needs a lookup to get the full attrs
                realtag, attrs = attrlookup[tag]
                startElement(realtag, attrs)
        
        else:
            # Normal tags
            startElement(tag, _attribute_list)
        
        # Recurse if it has any children
        if subtags:
            _do_callback(s, l, r, subtags, cont_handler, attrlookup)
        else:
            characters(s[l:r])
        begin = r

        if tag.startswith(">"):
            if tag.startswith(">G"):
                realtag, attrs = attrlookup[tag]
                endElement(realtag)
        else:
            endElement(tag)

    # anything after the last tag and before the end of the current
    # range are characters
    if begin < end:
        characters(s[begin:end])

def _do_dispatch_callback(s, begin, end, taglist,
                          start_table_get, cont_handler, save_stack,
                          end_table_get,
                          attrlookup):
    """internal function to convert the tagtable into ContentHandler events

    THIS IS A SPECIAL CASE FOR Dispatch.Dispatcher objects

    's' is the input text
    'begin' is the current position in the text
    'end' is 1 past the last position of the text allowed to be parsed
    'taglist' is the tag list from mxTextTools.parse
    'start_table_get' is the Dispatcher._start_table
    'cont_handler' is the Dispatcher
    'end_table_get' is the Dispatcher._end_table
    'cont_handler' is the SAX ContentHandler
    'attrlookup' is a dict mapping the encoded tag name to the element info
    """
    for tag, l, r, subtags in taglist:
        # If the tag's beginning is after the current position, then
        # the text from here to the tag's beginning are characters()
        assert begin <= l, "begin = %d and l = %d" % (begin, l)
        if begin < l and save_stack:
            cont_handler._save_text += s[begin:l]

        # Normal tags, see if the start function exists and call it
        #  ** This is a bit of a hack, in that this check also occurs
        #     with special tags.  But those begin with a '>' so will
        #     always fail.  This makes the logic a bit simpler and
        #     faster than checking the '>G' and '>ignore' terms.
        #     However, it is possible that specially constructed
        #     handlers could mess things up.  That cannot happen by
        #     accident, so I won't worry about it.
        # Yes, this reaches into the implementation of the Dispatcher.
        f = start_table_get(tag)
        if f is not None:
            f(tag, _attribute_list)
        else:
            # Tags with attributes
            x = attrlookup.get(tag)
            if x is not None:
                realtag, attrs = x
                # Does this function exist?
                f = start_table_get(realtag)
                if f is not None:
                    f(realtag, attrs)
        
        # Recurse if it has any children
        if subtags:
            _do_dispatch_callback(s, l, r, subtags,
                                  start_table_get,
                                  cont_handler, save_stack,
                                  end_table_get,
                                  attrlookup)
        elif save_stack:
            # Yes, this reaches into the implementation of the Dispatcher.
            cont_handler._save_text += s[l:r]
        begin = r

        # See if theres' a function for the normal tag
        f = end_table_get(tag)
        if f is not None:
            f(tag)
        else:
            # See if the special attribute tag exists
            x = attrlookup.get(tag)
            if x is not None:
                realtag, attrs = x
                # Yes, this reaches into the implementation of the Dispatcher.
                f = end_table_get(realtag)
                if f is not None:
                    f(realtag)

    # anything after the last tag and before the end of the current
    # range are characters
    if begin < end and save_stack:
        cont_handler._save_text += s[begin:end]

def _parse_elements(s, tagtable, cont_handler, debug_level, attrlookup):
    """parse the string with the tagtable and send the ContentHandler events

    Specifically, it sends the startElement, endElement and characters
    events but not startDocument and endDocument.
    """
    if debug_level:
        import Generate
        Generate._position = 0

    result, taglist, pos = TextTools.tag(s, tagtable, 0, len(s))

    # Special case test for the base ContentHandler since I know that
    # object does nothing and I want to test the method call overhead.
    if isinstance(cont_handler, Dispatch.Dispatcher):
        _do_dispatch_callback(s, 0, pos, taglist,
                              cont_handler._start_table.get,
                              cont_handler, cont_handler._save_stack,
                              cont_handler._end_table.get,
                              attrlookup)
    elif cont_handler.__class__ != handler.ContentHandler:
        # Send any tags to the client (there can be some even if there
        _do_callback(s, 0, pos, taglist, cont_handler, attrlookup)

    if not result:
        if debug_level:
            return ParserPositionException(Generate._position)
        else:
            return ParserPositionException(pos)
    elif pos != len(s):
        return pos
    else:
        return None

# This needs an interface like the standard XML parser
class Parser(xmlreader.XMLReader):
    """Parse the input data all in memory"""

    def __init__(self, tagtable, (want_groupref_names, debug_level, attrlookup) = (0, 1, {})):
        xmlreader.XMLReader.__init__(self)

        assert type(tagtable) == type( () ), "mxTextTools only allows a tuple tagtable"
        self.tagtable = tagtable

        # WARNING: This attribute is set directly by Generate - it bypasses
        # the value used in __init__.
        # Used to tell if the global "match_group" dict needs to be cleared.
        self.want_groupref_names = want_groupref_names

        self.debug_level = debug_level
        self.attrlookup = attrlookup

    def copy(self):
        parser = Parser(self.tagtable, (self.want_groupref_names,
                                        self.debug_level, self.attrlookup))
        parser.setContentHandler(self.getContentHandler())
        parser.setErrorHandler(self.getErrorHandler())
        parser.setDTDHandler(self.getDTDHandler())
        return parser
            
    def __str__(self):
        x = StringIO()
        pprint.pprint(self.tagtable, x)
        return x.getvalue()

    def parseFile(self, fileobj):
        """parse using the input file object

        XXX will be removed with the switch to Python 2.0, where parse()
        takes an 'InputSource'
        """
        # Just parse as a string
        self.parseString(fileobj.read())
    
    def parse(self, source):
        """parse using the URL or file handle"""
        source = saxutils.prepare_input_source(source)
        self.parseFile(source.getCharacterStream() or source.getByteStream())
        
    def parseString(self, s):
        """parse using the given string

        XXX will be removed with the switch to Python 2.0, where parse()
        takes an 'InputSource'
        """
        self._cont_handler.startDocument()

        if self.want_groupref_names:
            _match_group.clear()

        # parse the text and send the SAX events
        result = _parse_elements(s, self.tagtable, self._cont_handler,
                                 self.debug_level, self.attrlookup)

        if result is None:
            # Successful parse
            pass

        elif isinstance(result, _exceptions.SAXException):
            # could not parse record, and wasn't EOF
            self._err_handler.fatalError(result)
        
        else:
            # Parsed a record, but extra text remains
            pos = result
            self._err_handler.fatalError(ParserIncompleteException(pos))
        
        # Send an endDocument event even after errors
        self._cont_handler.endDocument()

    def close(self):
        pass

class RecordParser(xmlreader.XMLReader):
    """Parse the input data a record at a time"""
    def __init__(self, format_name, attrs, record_tagtable,
                 (want_groupref_names, debug_level, attrlookup),
                 make_reader, reader_args = ()):
        """parse the input data a record at a time

        format_name - XML tag name for the whole data file
        record_tagtable - mxTexTools tag table for each record
        want_groupref_names - flag to say if the match_group table needs to
              be reset (will disappear with better support from mxTextTools)

        make_reader - callable object which creates a RecordReader; first
              parameter will be an input file object
        reader_args - optional arguments to pass to make_reader after the
              input file object
        """
        xmlreader.XMLReader.__init__(self)
        
        self.format_name = format_name
        self.attrs = attrs
        assert type(record_tagtable) == type( () ), \
               "mxTextTools only allows a tuple tagtable"
        self.tagtable = record_tagtable
        self.want_groupref_names = want_groupref_names
        self.debug_level = debug_level
        self.attrlookup = attrlookup
        self.make_reader = make_reader
        self.reader_args = reader_args

    def copy(self):
        parser =  RecordParser(self.format_name, self.attrs, self.tagtable,
                               (self.want_groupref_names, self.debug_level,
                                self.attrlookup),
                               self.make_reader, self.reader_args)
        parser.setContentHandler(self.getContentHandler())
        parser.setErrorHandler(self.getErrorHandler())
        parser.setDTDHandler(self.getDTDHandler())
        return parser

    
    def __str__(self):
        x = StringIO()
        pprint.pprint(self.tagtable, x)
        return "parse records: " + x.getvalue()
    
    def parseFile(self, fileobj):
        """parse using the input file object

        XXX will be removed with the switch to Python 2.0, where parse()
        takes an 'InputSource'
        """
        self._cont_handler.startDocument()
        
        try:
            reader = self.make_reader( *(fileobj,) + self.reader_args)
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            # something unexpected happened
            # so call it a fatal error and stop
            outfile = StringIO()
            traceback.print_exc(file=outfile)
            self._err_handler.fatalError(ParserRecordException(
                outfile.getvalue(), sys.exc_info()[1]))
            self._cont_handler.endDocument()
            return

        if self.want_groupref_names:
            _match_group.clear()

        self._cont_handler.startElement(self.format_name, self.attrs)
        filepos = 0  # can get mixed up with DOS style "\r\n"
        while 1:
            try:
                record = reader.next()  
            except (KeyboardInterrupt, SystemExit):
                raise
            except:
                # something unexpected happened (couldn't find a record?)
                # so call it a fatal error and stop
                outfile = StringIO()
                traceback.print_exc(file=outfile)
                self._err_handler.fatalError(ParserRecordException(
                    outfile.getvalue(), sys.exc_info()[1]))
                self._cont_handler.endDocument()
                return
            
            if record is None:
                break
            result = _parse_elements(record, self.tagtable, self._cont_handler,
                                     self.debug_level, self.attrlookup)

            if result is None:
                # Successfully read the record
                pass
            elif isinstance(result, _exceptions.SAXException):
                # Wrong format or a SAX problem, but this is recoverable
                result += filepos
                self._err_handler.error(result)
            else:
                # Did not reach end of string, but this is recoverable
                pos = filepos + result
                self._err_handler.error(ParserPositionException(pos))

            filepos = filepos + len(record)

        self._cont_handler.endElement(self.format_name)
        self._cont_handler.endDocument()

    def parse(self, source):
        """parse using the URL or file handle"""
        source = saxutils.prepare_input_source(source)
        self.parseFile(source.getCharacterStream() or source.getByteStream())
        
    def parseString(self, s):
        """parse using the given string

        XXX will be removed with the switch to Python 2.0, where parse()
        takes an 'InputSource'
        """
        # Just parse it as a file
        strfile = StringIO(s)
        self.parseFile(strfile)

    def close(self):
        pass

class HeaderFooterParser(xmlreader.XMLReader):
    """Header followed by 0 or more records followed by a footer"""
    def __init__(self, format_name, attrs,
                 make_header_reader, header_reader_args, header_tagtable,
                 make_reader, reader_args, record_tagtable,
                 make_footer_reader, footer_reader_args, footer_tagtable,
                 (want_groupref_names, debug_level, attrlookup)):
        xmlreader.XMLReader.__init__(self)

        self.format_name = format_name
        self.attrs = attrs

        self.make_header_reader = make_header_reader
        self.header_reader_args = header_reader_args
        self.header_tagtable = header_tagtable

        self.make_reader = make_reader
        self.reader_args = reader_args
        self.record_tagtable = record_tagtable
        
        self.make_footer_reader = make_footer_reader
        self.footer_reader_args = footer_reader_args
        self.footer_tagtable = footer_tagtable
        
        self.want_groupref_names = want_groupref_names
        self.debug_level = debug_level
        self.attrlookup = attrlookup

    def __str__(self):
        x = StringIO()
        pprint.pprint( (self.header_tagtable, self.record_tagtable,
                        self.footer_tagtable), x)
        return "header footer records: " + x.getvalue()

    def copy(self):
        parser = HeaderFooterParser(self.format_name, self.attrs,
    self.make_header_reader, self.header_reader_args, self.header_tagtable,
    self.make_reader, self.reader_args, self.record_tagtable,
    self.make_footer_reader, self.footer_reader_args, self.footer_tagtable,
    (self.want_groupref_names, self.debug_level, self.attrlookup))

        parser.setContentHandler(self.getContentHandler())
        parser.setErrorHandler(self.getErrorHandler())
        parser.setDTDHandler(self.getDTDHandler())
        return parser


    def parseString(self, s):
        strfile = StringIO(s)
        self.parseFile(strfile)

    def parse(self, source):
        """parse using the URL or file handle"""
        source = saxutils.prepare_input_source(source)
        self.parseFile(source.getCharacterStream() or source.getByteStream())

    def parseFile(self, fileobj):
        self._cont_handler.startDocument()
        self._cont_handler.startElement(self.format_name, self.attrs)

        if self.want_groupref_names:
            _match_group.clear()

        # Read the header
        filepos = 0
        lookahead = ""
        if self.make_header_reader is not None:
            try:
                header_reader = self.make_header_reader(
                                     *(fileobj,) + self.header_reader_args)
                header = header_reader.next()
            except (KeyboardInterrupt, SystemExit):
                raise
            except:
                # Something unexpected happend so call it a fatal error
                outfile = StringIO()
                traceback.print_exc(file=outfile)
                exc = ParserRecordException(outfile.getvalue(),
                                            sys.exc_info()[1])
                self._err_handler.fatalError(exc)
                self._cont_handler.endDocument()
                return

            # Parse the text (if any) and send the SAX events
            if header is None:
                header = ""
            filepos += len(header)

            result = _parse_elements(header, self.header_tagtable,
                                     self._cont_handler, self.debug_level,
                                     self.attrlookup)
            if result is None:
                # Successful parse
                pass
            elif isinstance(result, _exceptions.SAXException):
                # Could not parse header and wasn't EOF
                self._err_handler.fatalError(result)
                self._cont_handler.endDocument()
                return
            else:
                # Reached EOF
                pos = result
                self._err_handler.fatalError(ParserPositionException(pos))
                self._cont_handler.endDocument()
                return

        # We've successfully parsed the header, now parse the records

        # Get any forward data from the header reader
        if self.make_header_reader is None:
            x, lookahead = fileobj, ""
        else:
            x, lookahead = header_reader.remainder()

        if self.make_footer_reader is None:
            # Only records - no footer
            try:
                reader = self.make_reader( *(fileobj,) + self.reader_args,
                             **{"lookahead": lookahead})
            except (KeyboardInterrupt, SystemExit):
                raise
            except:
                # Something unexpected happened so call it a fatal
                # error and stop
                outfile = StringIO()
                traceback.print_exc(file=outfile)
                exc = ParserRecordException(outfile.getvalue(),
                                            sys.exc_info()[1])
                self._err_handler.fatalError(exc)
                self._cont_handler.endDocument()
                return

            while 1:
                try:
                    record = reader.next()
                except (KeyboardInterrupt, SystemExit):
                    raise
                except:
                    # Something unexpected happened and I cannot recover
                    outfile = StringIO()
                    traceback.print_exc(file=outfile)
                    exc = ParserRecordException(outfile.getvalue(),
                                                sys.exc_info()[1])
                    self._err_handler.fatalError(exc)
                    self._cont_handler.endDocument()
                    return

                if record is None:
                    # Reached EOF, so that's it (since there's no footer)
                    self._cont_handler.endElement(self.format_name)
                    self._cont_handler.endDocument()
                    return

                result = _parse_elements(record, self.record_tagtable,
                                         self._cont_handler, self.debug_level,
                                         self.attrlookup)
                if result is None:
                    # Successfully parsed the record
                    pass
                else:
                    # Failed to parse the record, but can recover
                    if isinstance(result, _exceptions.SAXException):
                        result += filepos
                    else:
                        result = ParserPositionException(filepos + result)
                    self._err_handler.error(result)

                filepos += len(record)
                    
        assert self.make_footer_reader is not None, "internal error"
            
        # This gets to be quite complicated :(

        # If the record fails, try the footer.  If that fails,
        # skip the record and try again
        record_exc = None
        try:
            reader = self.make_reader( *(fileobj,) + self.reader_args,
                         **{"lookahead": lookahead})
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            # Something unexpected happened - could be that there was
            # no record and only a footer?  Save the current exception.
            outfile = StringIO()
            traceback.print_exc(file=outfile)
            record_exc = ParserRecordException(outfile.getvalue(),
                                               sys.exc_info()[1])

        while record_exc is None:
            try:
                record = reader.next()
            except (KeyboardInterrupt, SystemExit):
                raise
            except:
                # Something unexpected happened.  Could be the footer,
                # but save the current exception in case it isn't
                outfile = StringIO()
                traceback.print_exc(file=outfile)
                record_exc = ParserRecordException(outfile.getvalue(),
                                                   sys.exc_info()[1])
                break

            if record is None:
                # Reached EOF, but there should have been a footer
                record_exc = ParserPositionException(filepos)
                break

            result = _parse_elements(record, self.record_tagtable,
                                     self._cont_handler, self.debug_level,
                                     self.attrlookup)
            if result is None:
                # Successfully parsed the record
                pass
            else:
                # Failed to parse the record, but may recover of it
                # isn't the footer
                if isinstance(result, _exceptions.SAXException):
                    result += filepos
                else:
                    result = ParserPositionException(filepos + result)
                record_exc = result

                # Is there a valid footer?
                try:
                    footer = ""
                    x, lookahead = reader.remainder()
                    footer_reader = self.make_footer_reader(
                                         *(fileobj,) + self.footer_reader_args,
                                        **{"lookahead": record + lookahead})
                    footer = footer_reader.next()
                except (KeyboardInterrupt, SystemExit):
                    raise
                except:
                    # Not a footer either, so call this an error and
                    # attempt the next record
                    self._err_handler.error(record_exc)
                    record_exc = None

                    # But that means I need to reset the record reader(!)
                    x, lookahead = footer_reader.remainder()
                    try:
                        reader = self.make_reader(
                                      *(fileobj,) + self.reader_args,
                                     **{"lookahead": footer + lookahead})
                    except (KeyboardInterrupt, SystemExit):
                        raise
                    except:
                        # Something unexpected happened.  Save the
                        # current exception and stop reading
                        outfile = StringIO()
                        traceback.print_exc(file=outfile)
                        record_exc = ParserRecordException(outfile.getvalue(),
                                                           sys.exc_info()[1])
                        break
                    
                    

                # Hmm, okay, it was a valid footer, but can be it be
                # parsed?
                result = _parse_elements(footer, self.footer_tagtable,
                                         self._cont_handler, self.debug_level,
                                         self.attrlookup)

                if result is None:
                    # parsed the footer, but need to check that it's
                    # at EOF
                    x, remainder = footer_reader.remainder()
                    if remainder or x.read(1):
                        # Acck, there's data left over
                        record_exc = ParserPositionException(filepos +
                                                             len(footer))
                        self._err_handler.fatalError(record_exc)
                        self._cont_handler.endDocument()
                        return
                    # Success!
                    self._cont_handler.endElement(self.format_name)
                    self._cont_handler.endDocument()
                    return
                else:
                    # Wasn't a footer, so reset the reader stream and skip
                    # past the record which I know I can read.
                    x, remainder = footer_reader.remainder()
                    reader = self.make_reader(
                                  *(fileobj, ) + self.reader_args,
                                 **{"lookahead": footer + remainder})
                    record = reader.next()
                    self._err_handler.error(record_exc)
                    record_exc = None

            filepos += len(record)
                    
        # Could not read a record or reached EOF.  Try to parse the
        # trailer
        x, remainder = reader.remainder()
        try:
            footer_reader = self.make_footer_reader(
                                 *(fileobj,) + self.footer_reader_args,
                                **{"lookahead": remainder})
            footer = footer_reader.next()
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            # Cannot read the record, so use the older error
            self._err_handler.fatalError(record_exc)
            self._cont_handler.endDocument()
            return

        if footer is None:
            footer = ""
        result = _parse_elements(footer, self.footer_tagtable,
                                 self._cont_handler, self.debug_level,
                                 self.attrlookup)
        if result is None:
            # parsed the footer, but need to check that it's
            # at EOF
            x, remainder = footer_reader.remainder()
            if remainder or x.read(1):
                # Acck, there's data left over
                record_exc = ParserPositionException(filepos +
                                                     len(footer))
                self._err_handler.fatalError(record_exc)
                self._cont_handler.endDocument()
                return
            # Success!
            self._cont_handler.endElement(self.format_name)
            self._cont_handler.endDocument()
            return
        else:
            # Okay, use the old error
            self._err_handler.fatalError(record_exc)
            self._cont_handler.endDocument()
            return
