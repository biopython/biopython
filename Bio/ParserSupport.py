# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Code to support writing parsers.



Classes:
AbstractConsumer Base class of all Consumers.
TaggingConsumer  Consumer that tags output with its event.

Functions:
safe_readline    Read a line from a handle, with check for EOF.
safe_peekline    Peek at next line, with check for EOF.
read_and_call    Read a line from a handle and pass it to a method.
attempt_read_and_call  Like read_and_call, but forgiving of errors.
is_blank_line    Test whether a line is blank.

"""

import sys
import string
import traceback

from Bio import File


class AbstractConsumer:
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
    def __init__(self, handle=sys.stdout, colwidth=15, maxwidth=80):
        """__init__(self, handle=sys.stdout)"""
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
                string.rstrip(data[:self._maxwidth-self._colwidth-2])))

    def __getattr__(self, attr):
        if attr[:6] == 'start_' or attr[:4] == 'end_':
            method = lambda a=attr, s=self: s._print_name(a)
        else:
            method = lambda x, a=attr, s=self: s._print_name(a, x)
        return method

def read_and_call(uhandle, method,
                  start=None, end=None, contains=None, blank=None):
    """_read_and_call(uhandle, method,
    start=None, end=None, contains=None, blank=None)

    Read a line from uhandle, check it, and pass it to the method.
    Raises a SyntaxError if the line does not pass the checks.

    start, end, contains, and blank specify optional conditions that
    the line must pass.  start and end specifies what the line must
    begin or end with (not counting EOL characters).  contains
    specifies a substring that must be found in the line.  If blank
    is a true value, then the line must be blank.  Set these parameters
    to None if the check is not necessary.

    """
    line = uhandle.readline()
    if not line:
        raise SyntaxError, "Unexpected end of stream."
    if start:
        if line[:len(start)] != start:
            raise SyntaxError, "Line does not start with '%s': %s" % \
                  (start, line)
    if end:
        if string.rstrip(line)[-len(end):] != end:
            raise SyntaxError, "Line does not end with '%s': %s" % \
                  (end, line)
    if contains:
        if string.find(line, contains) == -1:
            raise SyntaxError, "Line does not contain '%s': %s" % \
                  (contains, line)
    if blank is not None:
        if blank:
            if not is_blank_line(line):
                raise SyntaxError, "Expected blank line, but got: %s" % \
                      line
        else:
            if is_blank_line(line):
                raise SyntaxError, "Expected non-blank line, but got a blank one"
    apply(method, (line,))

def attempt_read_and_call(uhandle, method, **keywds):
    """attempt_read_and_call(uhandle, method, **keywds) -> boolean

    Similar to read_and_call, but returns a boolean specifying
    whether the line has passed the checks.  Does not raise
    exceptions.

    See docs for read_and_call for a description of the function
    arguments.

    """
    # Delegate as much work as possible to read_and_call.  If the
    # line fails the test inside that function, I will need to catch
    # the SyntaxError exception and restore the state of the
    # buffer.  Thus, I will need to first take a peek at the next line
    # so I can restore it if necessary.
    line = uhandle.peekline()
        
    try:
        apply(read_and_call, (uhandle, method), keywds)
    except SyntaxError:
        # I only want to catch the exception if it was raised by
        # the read_and_call method.  Thus, I will examine the traceback.
        # If it contains more than 2 stack frames, then pass the exception on.
        try:
            tb = sys.exc_info()[2]
            if len(traceback.extract_tb(tb)) > 2:
                raise
        finally:
            del tb  # prevent circular reference to tb
        
        # read_and_call has read out a line and failed.
        # Put the line back in the buffer.
        uhandle.saveline(line)
        return 0
    return 1

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
        return string.rstrip(line) == ''
    return line[0] == '\n' or line[0] == '\r'

def safe_readline(handle):
    """safe_readline(handle) -> line

    Read a line from an UndoHandle and return it.  If there are no more
    lines to read, I will raise a SyntaxError.

    """
    line = handle.readline()
    if not line:
        raise SyntaxError, "Unexpected end of stream."
    return line

def safe_peekline(handle):
    """safe_peekline(handle) -> line

    Peek at the next line in an UndoHandle and return it.  If there are no
    more lines to peek, I will raise a SyntaxError.
    
    """
    line = handle.peekline()
    if not line:
        raise SyntaxError, "Unexpected end of stream."
    return line
