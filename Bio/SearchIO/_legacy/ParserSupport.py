# Copyright 1999 by Jeffrey Chang.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Code to support writing parsers (DEPRECATED).

Classes:
 - UndoHandle             File object decorator with support for undo-like operations.
 - AbstractParser         Base class for parsers.
 - AbstractConsumer       Base class of all Consumers.
 - TaggingConsumer        Consumer that tags output with its event.  For debugging

Functions:
 - safe_readline          Read a line from a handle, with check for EOF.
 - safe_peekline          Peek at next line, with check for EOF.
 - read_and_call          Read a line from a handle and pass it to a method.
 - read_and_call_while    Read many lines, as long as a condition is met.
 - read_and_call_until    Read many lines, until a condition is met.
 - attempt_read_and_call  Like read_and_call, but forgiving of errors.
 - is_blank_line          Test whether a line is blank.

"""

import sys
from io import StringIO

from abc import ABC, abstractmethod


class UndoHandle:
    """A Python handle that adds functionality for saving lines.

    Saves lines in a LIFO fashion.
    """

    def __init__(self, handle):
        """Initialize the class."""
        self._handle = handle
        self._saved = []
        try:
            # If wrapping an online handle, this this is nice to have:
            self.url = handle.url
        except AttributeError:
            pass

    def __iter__(self):
        """Iterate over the lines in the File."""
        return self

    def __next__(self):
        """Return the next line."""
        next = self.readline()
        if not next:
            raise StopIteration
        return next

    def readlines(self, *args, **keywds):
        """Read all the lines from the file as a list of strings."""
        lines = self._saved + self._handle.readlines(*args, **keywds)
        self._saved = []
        return lines

    def readline(self, *args, **keywds):
        """Read the next line from the file as string."""
        if self._saved:
            line = self._saved.pop(0)
        else:
            line = self._handle.readline(*args, **keywds)
        return line

    def read(self, size=-1):
        """Read the File."""
        if size == -1:
            saved = "".join(self._saved)
            self._saved[:] = []
        else:
            saved = ""
            while size > 0 and self._saved:
                if len(self._saved[0]) <= size:
                    size = size - len(self._saved[0])
                    saved = saved + self._saved.pop(0)
                else:
                    saved = saved + self._saved[0][:size]
                    self._saved[0] = self._saved[0][size:]
                    size = 0
        return saved + self._handle.read(size)

    def saveline(self, line):
        """Store a line in the cache memory for later use.

        This acts to undo a readline, reflecting the name of the class: UndoHandle.
        """
        if line:
            self._saved = [line] + self._saved

    def peekline(self):
        """Return the next line in the file, but do not move forward though the file."""
        if self._saved:
            line = self._saved[0]
        else:
            line = self._handle.readline()
            self.saveline(line)
        return line

    def tell(self):
        """Return the current position of the file read/write pointer within the File."""
        return self._handle.tell() - sum(len(line) for line in self._saved)

    def seek(self, *args):
        """Set the current position at the offset specified."""
        self._saved = []
        self._handle.seek(*args)

    def __getattr__(self, attr):
        """Return File attribute."""
        return getattr(self._handle, attr)

    def __enter__(self):
        """Call special method when opening the file using a with-statement."""
        return self

    def __exit__(self, type, value, traceback):
        """Call special method when closing the file using a with-statement."""
        self._handle.close()


class AbstractParser(ABC):
    """Abstract base class for other parsers."""

    @abstractmethod
    def parse(self, handle):
        """Provision for parsing a file handle."""
        raise NotImplementedError

    def parse_str(self, string):
        """Make string a handle, so it can be taken by parse."""
        return self.parse(StringIO(string))

    def parse_file(self, filename):
        """Parse a file, open the file as handle so it can be taken by parse."""
        with open(filename) as h:
            retval = self.parse(h)
        return retval


class AbstractConsumer:
    """Base class for other Consumers.

    Derive Consumers from this class and implement appropriate
    methods for each event that you want to receive.

    """

    # Optionally implement in the sub-class
    def _unhandled_section(self):
        pass

    # Optionally implement in the sub-class
    def _unhandled(self, data):
        pass

    def __getattr__(self, attr):
        if attr[:6] == "start_" or attr[:4] == "end_":
            method = self._unhandled_section
        else:
            method = self._unhandled
        return method


class TaggingConsumer(AbstractConsumer):
    """Debugging consumer which tags data with the event and logs it.

    This is a Consumer that tags the data stream with the event and
    prints it to a handle.  Useful for debugging.

    """

    def __init__(self, handle=None, colwidth=15, maxwidth=80):
        """Initialize.

        Arguments:
         - handle to log to, defaults to ``sys.stdout``
         - colwidth for logging to the handle
         - maxwidth for truncation when logging

        """
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
        """Tag an unhandled section."""
        self._print_name("unhandled_section")

    def unhandled(self, data):
        """Tag unhandled data."""
        self._print_name("unhandled", data)

    def _print_name(self, name, data=None):
        if data is None:
            # Write the name of a section.
            self._handle.write("%s %s\n" % ("*" * self._colwidth, name))
        else:
            # Write the tag and line.
            self._handle.write(
                "%-*s: %s\n"
                % (
                    self._colwidth,
                    name[: self._colwidth],
                    data[: self._maxwidth - self._colwidth - 2].rstrip(),
                )
            )

    def __getattr__(self, attr):
        if attr[:6] == "start_" or attr[:4] == "end_":
            method = lambda a=attr, s=self: s._print_name(a)  # noqa: E731
        else:
            method = lambda x, a=attr, s=self: s._print_name(a, x)  # noqa: E731
        return method


def read_and_call(uhandle, method, **keywds):
    """Read line and pass it to the method.

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
    """Read line and pass it to the method while condition is true.

    Read a line from uhandle and pass it to the method as long as
    some condition is true.  Returns the number of lines that were read.

    See the docstring for read_and_call for a description of the parameters.

    """
    nlines = 0
    while True:
        line = safe_readline(uhandle)
        # If I've failed the condition, then stop reading the line.
        if _fails_conditions(*(line,), **keywds):
            uhandle.saveline(line)
            break
        method(line)
        nlines = nlines + 1
    return nlines


def read_and_call_until(uhandle, method, **keywds):
    """Read line and pass it to the method until condition is true.

    Read a line from uhandle and pass it to the method until
    some condition is true.  Returns the number of lines that were read.

    See the docstring for read_and_call for a description of the parameters.

    """
    nlines = 0
    while True:
        line = safe_readline(uhandle)
        # If I've met the condition, then stop reading the line.
        if not _fails_conditions(*(line,), **keywds):
            uhandle.saveline(line)
            break
        method(line)
        nlines = nlines + 1
    return nlines


def attempt_read_and_call(uhandle, method, **keywds):
    """Attempt read line and call method.

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


def _fails_conditions(
    line, start=None, end=None, contains=None, blank=None, has_re=None
):
    if start is not None:
        if line[: len(start)] != start:
            return "Line does not start with '%s':\n%s" % (start, line)
    if end is not None:
        if line.rstrip()[-len(end) :] != end:
            return "Line does not end with '%s':\n%s" % (end, line)
    if contains is not None:
        if contains not in line:
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
            return "Line does not match regex '%s':\n%s" % (has_re.pattern, line)
    return None


def is_blank_line(line, allow_spaces=0):
    """Check if a line is blank.

    Return whether a line is blank.  allow_spaces specifies whether to
    allow whitespaces in a blank line.  A true value signifies that a
    line containing whitespaces as well as end-of-line characters
    should be considered blank.

    """
    if not line:
        return 1
    if allow_spaces:
        return line.rstrip() == ""
    return line[0] == "\n" or line[0] == "\r"


def safe_readline(handle):
    """Read a line, otherwise raises ValueError.

    Read a line from an UndoHandle and return it.  If there are no more
    lines to read, I will raise a ValueError.

    """
    line = handle.readline()
    if not line:
        raise ValueError("Unexpected end of stream.")
    return line


def safe_peekline(handle):
    """Peek at the next line if present, otherwise raises ValueError.

    Peek at the next line in an UndoHandle and return it.  If there are no
    more lines to peek, I will raise a ValueError.

    """
    line = handle.peekline()
    if not line:
        raise ValueError("Unexpected end of stream.")
    return line
