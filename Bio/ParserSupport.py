# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Code to support writing parsers (OBSOLETE).



Classes:
AbstractParser         Base class for parsers.
AbstractConsumer       Base class of all Consumers.
TaggingConsumer        Consumer that tags output with its event.  For debugging
SGMLStrippingConsumer  Consumer that strips SGML tags from output.
EventGenerator         Generate Biopython Events from Martel XML output
                       (note that Martel is now DEPRECATED)

Functions:
safe_readline          Read a line from a handle, with check for EOF.
safe_peekline          Peek at next line, with check for EOF.
read_and_call          Read a line from a handle and pass it to a method.
read_and_call_while    Read many lines, as long as a condition is met.
read_and_call_until    Read many lines, until a condition is met.
attempt_read_and_call  Like read_and_call, but forgiving of errors.
is_blank_line          Test whether a line is blank.

"""


import warnings
warnings.warn("The module Bio.ParserSupport is now obsolete, and will be deprecated and removed in a future release of Biopython.", PendingDeprecationWarning)


import sys
import traceback
from types import *
import StringIO

from Bio import File

# XML from python 2.0
try:
    from xml.sax import handler
    xml_support = 1
except ImportError:
    sys.stderr.write("Warning: Could not import SAX for dealing with XML.\n" +
                     "This causes problems with some ParserSupport modules\n")
    xml_support = 0

class AbstractParser(object):
    """Base class for other parsers.

    """
    def parse(self, handle):
        raise NotImplementedError("Please implement in a derived class")

    def parse_str(self, string):
        return self.parse(StringIO.StringIO(string))

    def parse_file(self, filename):
        h = open(filename)
        try:
            retval = self.parse(h)
        finally:
            h.close()
        return retval

class AbstractConsumer(object):
    """Base class for other Consumers.

    Derive Consumers from this class and implement appropriate
    methods for each event that you want to receive.
    
    """
    def _unhandled_section(self):
        pass
    def _unhandled(self, data):
        pass
    def __getattr__(self, attr):
        if attr[:6] == 'start_' or attr[:4] == 'end_':
            method = self._unhandled_section
        else:
            method = self._unhandled
        return method

class TaggingConsumer(AbstractConsumer):
    """A Consumer that tags the data stream with the event and
    prints it to a handle.  Useful for debugging.

    """
    def __init__(self, handle=None, colwidth=15, maxwidth=80):
        """TaggingConsumer(handle=sys.stdout, colwidth=15, maxwidth=80)"""
        # I can't assign sys.stdout to handle in the argument list.
        # If I do that, handle will be assigned the value of sys.stdout
        # the first time this function is called.  This will fail if
        # the user has assigned sys.stdout to some other file, which may
        # be closed or invalid at a later time.
        if handle is None:
            handle = sys.stdout
        self._handle = handle
        self._colwidth = colwidth
        self._maxwidth = maxwidth

    def unhandled_section(self):
        self._print_name('unhandled_section')

    def unhandled(self, data):
        self._print_name('unhandled', data)

    def _print_name(self, name, data=None):
        if data is None:
            # Write the name of a section.
            self._handle.write("%s %s\n" % ("*"*self._colwidth, name))
        else:
            # Write the tag and line.
            self._handle.write("%-*s: %s\n" % (
                self._colwidth, name[:self._colwidth],
                data[:self._maxwidth-self._colwidth-2].rstrip()))

    def __getattr__(self, attr):
        if attr[:6] == 'start_' or attr[:4] == 'end_':
            method = lambda a=attr, s=self: s._print_name(a)
        else:
            method = lambda x, a=attr, s=self: s._print_name(a, x)
        return method

class SGMLStrippingConsumer(object):
    """A consumer that strips off SGML tags.

    This is meant to be used as a decorator for other consumers.

    """
    def __init__(self, consumer):
        import Bio
        warnings.warn("SGMLStrippingConsumer is deprecated, and is likely to be removed in a future version of Biopython", Bio.BiopythonDeprecationWarning)
        if type(consumer) is not InstanceType:
            raise ValueError("consumer should be an instance")
        self._consumer = consumer
        self._prev_attr = None
        self._stripper = File.SGMLStripper()

    def _apply_clean_data(self, data):
        clean = self._stripper.strip(data)
        self._prev_attr(clean)

    def __getattr__(self, name):
        if name in ['_prev_attr', '_stripper']:
            return getattr(self, name)
        attr = getattr(self._consumer, name)
        # If this is not a method, then return it as is.
        if type(attr) is not MethodType:
            return attr
        # If it's a section method, then return it.
        if name[:6] == 'start_' or name[:4] == 'end_':
            return attr
        # Otherwise, it's an info event, and return my method.
        self._prev_attr = attr
        return self._apply_clean_data

# onle use the Event Generator if XML handling is okay
if xml_support:
    class EventGenerator(handler.ContentHandler):
        """Handler to generate events associated with a Martel parsed file.

        This acts like a normal SAX handler, and accepts XML generated by
        Martel during parsing. These events are then converted into
        'Biopython events', which can then be caught by a standard
        biopython consumer.

        Note that Martel is now DEPRECATED.
        """
        def __init__(self, consumer, interest_tags, callback_finalizer = None,
                     exempt_tags = []):
            """Initialize to begin catching and firing off events.

            Arguments:
            o consumer - The consumer that we'll send Biopython events to.
            
            o interest_tags - A listing of all the tags we are interested in.

            o callback_finalizer - A function to deal with the collected
            information before passing it on to the consumer. By default
            the collected information is a list of all of the lines read
            for a particular tag -- if there are multiple tags in a row
            like:

            <some_info>Spam<some_info>
            <some_info>More Spam<some_info>

            In this case the list of information would be:

            ['Spam', 'More Spam']
            
            This list of lines will be passed to the callback finalizer if
            it is present. Otherwise the consumer will be called with the
            list of content information.

            o exempt_tags - A listing of particular tags that are exempt from
            being processed by the callback_finalizer. This allows you to
            use a finalizer to deal with most tags, but leave those you don't
            want touched.
            """
            self._consumer = consumer
            self.interest_tags = interest_tags
            self._finalizer = callback_finalizer
            self._exempt_tags = exempt_tags

            # a dictionary of content for each tag of interest
            # the information for each tag is held as a list of the lines.
            # This allows us to collect information from multiple tags
            # in a row, and return it all at once.
            self.info = {}
            for tag in self.interest_tags:
                self.info[tag] = []

            # the previous tag we were collecting information for.
            # We set a delay in sending info to the consumer so that we can
            # collect a bunch of tags in a row and append all of the info
            # together.
            self._previous_tag = ''

            # the current character information for a tag
            self._cur_content = []
            # whether we should be collecting information
            self._collect_characters = 0

        def startElement(self, name, attrs):
            """Determine if we should collect characters from this tag.
            """
            if name in self.interest_tags:
                self._collect_characters = 1

        def characters(self, content):
            """Extract the information if we are interested in it.
            """
            if self._collect_characters:
                self._cur_content.append(content)

        def endElement(self, name):
            """Send the information to the consumer.

            Once we've got the end element we've collected up all of the
            character information we need, and we need to send this on to
            the consumer to do something with it.

            We have a delay of one tag on doing this, so that we can collect
            all of the info from multiple calls to the same element at once.
            """
            # only deal with the tag if it is something we are
            # interested in and potentially have information for
            if self._collect_characters:
                # add all of the information collected inside this tag
                self.info[name].append("".join(self._cur_content))
                # reset our information and flags
                self._cur_content = []
                self._collect_characters = 0
                
                # if we are at a new tag, pass on the info from the last tag
                if self._previous_tag and self._previous_tag != name:
                    self._make_callback(self._previous_tag)

                # set this tag as the next to be passed
                self._previous_tag = name

        def _make_callback(self, name):
            """Call the callback function with the info with the given name.
            """
            # strip off whitespace and call the consumer
            callback_function = getattr(self._consumer, name)

            # --- pass back the information
            # if there is a finalizer, use that
            if self._finalizer is not None and name not in self._exempt_tags:
                info_to_pass = self._finalizer(self.info[name])
            # otherwise pass back the entire list of information
            else:
                info_to_pass = self.info[name]
            
            callback_function(info_to_pass)

            # reset the information for the tag
            self.info[name] = []

        def endDocument(self):
            """Make sure all of our information has been passed.

            This just flushes out any stored tags that need to be passed.
            """
            if self._previous_tag:
                self._make_callback(self._previous_tag)

def read_and_call(uhandle, method, **keywds):
    """read_and_call(uhandle, method[, start][, end][, contains][, blank][, has_re])

    Read a line from uhandle, check it, and pass it to the method.
    Raises a ValueError if the line does not pass the checks.

    start, end, contains, blank, and has_re specify optional conditions
    that the line must pass.  start and end specifies what the line must
    begin or end with (not counting EOL characters).  contains
    specifies a substring that must be found in the line.  If blank
    is a true value, then the line must be blank.  has_re should be
    a regular expression object with a pattern that the line must match
    somewhere.

    """
    line = safe_readline(uhandle)
    errmsg = _fails_conditions(*(line,), **keywds)
    if errmsg is not None:
        raise ValueError(errmsg)
    method(line)

def read_and_call_while(uhandle, method, **keywds):
    """read_and_call_while(uhandle, method[, start][, end][, contains][, blank][, has_re]) -> number of lines

    Read a line from uhandle and pass it to the method as long as
    some condition is true.  Returns the number of lines that were read.

    See the docstring for read_and_call for a description of the parameters.
    
    """
    nlines = 0
    while 1:
        line = safe_readline(uhandle)
        # If I've failed the condition, then stop reading the line.
        if _fails_conditions(*(line,), **keywds):
            uhandle.saveline(line)
            break
        method(line)
        nlines = nlines + 1
    return nlines

def read_and_call_until(uhandle, method, **keywds):
    """read_and_call_until(uhandle, method, 
    start=None, end=None, contains=None, blank=None) -> number of lines

    Read a line from uhandle and pass it to the method until
    some condition is true.  Returns the number of lines that were read.

    See the docstring for read_and_call for a description of the parameters.
    
    """
    nlines = 0
    while 1:
        line = safe_readline(uhandle)
        # If I've met the condition, then stop reading the line.
        if not _fails_conditions(*(line,), **keywds):
            uhandle.saveline(line)
            break
        method(line)
        nlines = nlines + 1
    return nlines

def attempt_read_and_call(uhandle, method, **keywds):
    """attempt_read_and_call(uhandle, method, **keywds) -> boolean

    Similar to read_and_call, but returns a boolean specifying
    whether the line has passed the checks.  Does not raise
    exceptions.

    See docs for read_and_call for a description of the function
    arguments.

    """
    line = safe_readline(uhandle)
    passed = not _fails_conditions(*(line,), **keywds)
    if passed:
        method(line)
    else:
        uhandle.saveline(line)
    return passed

def _fails_conditions(line, start=None, end=None, contains=None, blank=None,
                      has_re=None):
    if start is not None:
        if line[:len(start)] != start:
            return "Line does not start with '%s':\n%s" % (start, line)
    if end is not None:
        if line.rstrip()[-len(end):] != end:
            return "Line does not end with '%s':\n%s" % (end, line)
    if contains is not None:
        if line.find(contains) == -1:
            return "Line does not contain '%s':\n%s" % (contains, line)
    if blank is not None:
        if blank:
            if not is_blank_line(line):
                return "Expected blank line, but got:\n%s" % line
        else:
            if is_blank_line(line):
                return "Expected non-blank line, but got a blank one"
    if has_re is not None:
        if has_re.search(line) is None:
            return "Line does not match regex '%s':\n%s" % (
                has_re.pattern, line)
    return None

def is_blank_line(line, allow_spaces=0):
    """is_blank_line(line, allow_spaces=0) -> boolean

    Return whether a line is blank.  allow_spaces specifies whether to
    allow whitespaces in a blank line.  A true value signifies that a
    line containing whitespaces as well as end-of-line characters
    should be considered blank.

    """
    if not line:
        return 1
    if allow_spaces:
        return line.rstrip() == ''
    return line[0] == '\n' or line[0] == '\r'

def safe_readline(handle):
    """safe_readline(handle) -> line

    Read a line from an UndoHandle and return it.  If there are no more
    lines to read, I will raise a ValueError.

    """
    line = handle.readline()
    if not line:
        raise ValueError("Unexpected end of stream.")
    return line

def safe_peekline(handle):
    """safe_peekline(handle) -> line

    Peek at the next line in an UndoHandle and return it.  If there are no
    more lines to peek, I will raise a ValueError.
    
    """
    line = handle.peekline()
    if not line:
        raise ValueError("Unexpected end of stream.")
    return line
