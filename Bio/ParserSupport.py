# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Code to support writing parsers.



Classes:
Consumer         Base class of all Consumers.
TaggingConsumer  Consumer that tags output with its event.
OopsHandle       File object decorator with support for undo-like things.

Functions:
read_and_call    Read a line from a handle and pass it to a method.
attempt_read_and_call  Like read_and_call, but forgiving of errors.
is_blank_line    Test whether a line is blank.
"""

import sys
import string


class Consumer:
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

class TaggingConsumer(Consumer):
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
                string.rstrip(data[:self._maxwidth-self._colwidth])))

    def __getattr__(self, attr):
        if attr[:6] == 'start_' or attr[:4] == 'end_':
            method = lambda a=attr, s=self: s._print_name(a)
        else:
            method = lambda x, a=attr, s=self: s._print_name(a, x)
        return method

class OopsHandle:
    """A Python handle that adds functionality for saving lines.

    XXX Bug: I don't handle the arguments to readline or readlines
    correctly, since their behavior isn't documented anywhere.
    I'll just pass the arguments along, and hope I don't do anything
    too bad...

    Added methods:
    saveline    Save a line to be returned next time.
    peekline    Peek at the next line without consuming it.

    """
    def __init__(self, handle):
        self._handle = handle
        self._saved = []

    def readlines(self, *args, **keywds):
        if self._saved:
            lines = self._saved + apply(self._handle.readlines, args, keywds)
            self._saved = []
            return lines
        return apply(self._handle.readlines, args, keywds)

    def readline(self, *args, **keywds):
        if self._saved:
            line = self._saved.pop(0)
        else:
            line = apply(self._handle.readline, args, keywds)
        return line

    def saveline(self, line):
        self._saved.append(line)

    def peekline(self):
        if self._saved:
            line = self._saved[0]
        else:
            line = self._handle.readline()
            self._saved.append(line)
        return line

    def __getattr__(self, attr):
        return getattr(self._handle, attr)

def read_and_call(ohandle, method,
                  start=None, end=None, contains=None, blank=None):
    """_read_and_call(ohandle, method,
    start=None, end=None, contains=None, blank=None)

    Read a line from ohandle, check it, and pass it to the method.
    Raises a SyntaxError if the line does not pass the checks.

    start, end, contains, and blank specify optional conditions that
    the line must pass.  start and end specifies what the line must
    begin or end with (not counting EOL characters).  contains
    specifies a substring that must be found in the line.  If blank
    is a true value, then the line must be blank.  Set these parameters
    to None if the check is not necessary.

    """
    line = ohandle.readline()
    if not line:
        raise SyntaxError, "Unexpected end of blast report."
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
    if blank:
        if not is_blank_line(line):
            raise SyntaxError, "Expected blank line, but got: %s" % \
                  line
    apply(method, (line,))

def attempt_read_and_call(ohandle, method, **keywds):
    """attempt_read_and_call(ohandle, method, **keywds) -> boolean

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
    #
    # XXX Bug here: If the method raises a SyntaxError, I will erroneously
    # catch that too.  
    line = ohandle.peekline()
        
    try:
        apply(read_and_call, (ohandle, method), keywds)
    except SyntaxError:
        # read_and_call has read out a line and failed.
        # Put the line back in the buffer.
        ohandle.saveline(line)
        return 0
    return 1

def is_blank_line(line, allow_spaces=0):
    """is_blank_line(line, allow_spaces=0) -> boolean

    Return whether a line is blank.  allow_spaces specifies whether to
    allow whitespaces in a blank line.  A true value signifies that a
    line containing whitespaces as well as end-of-line characters
    should be considered blank.

    """
    if allow_spaces:
        return string.rstrip(line) == ''
    return line[0] == '\n' or line[0] == '\r'
