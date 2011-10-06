# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Code for more fancy file handles.


Classes:

UndoHandle     File object decorator with support for undo-like operations.

StringHandle   Wraps a file object around a string.  This is now DEPRECATED,
               and is likely to be removed in a future release of Biopython.

SGMLStripper   Object that strips SGML.  This is now DEPRECATED, and is likely
               to be removed in a future release of Biopython.

"""
# For with statement in Python 2.5
from __future__ import with_statement
import contextlib
import StringIO


@contextlib.contextmanager
def as_handle(handleish, mode='r', **kwargs):
    """
    Context manager for arguments that can be passed to
    SeqIO and AlignIO read, write, and parse methods: either file objects or strings.

    When given a string, returns a file handle open to handleish with provided
    mode which will be closed when the manager exits.

    All other inputs are returned, and are *not* closed

    - handleish  - Either a string or file handle
    - mode       - Mode to open handleish (used only if handleish is a string)
    - kwargs     - Further arguments to pass to open(...)

    Example:

    >>> with as_handle('seqs.fasta', 'w') as fp:
    ...     fp.write('>test\nACGT')
    >>> fp.closed
    True

    >>> handle = open('seqs.fasta', 'w')
    >>> with as_handle(handle) as fp:
    ...     fp.write('>test\nACGT')
    >>> fp.closed
    False
    >>> fp.close()
    """
    if isinstance(handleish, basestring):
        with open(handleish, mode, **kwargs) as fp:
            yield fp
    else:
        yield handleish


class UndoHandle(object):
    """A Python handle that adds functionality for saving lines.

    Saves lines in a LIFO fashion.

    Added methods:
    saveline    Save a line to be returned next time.
    peekline    Peek at the next line without consuming it.

    """
    def __init__(self, handle):
        self._handle = handle
        self._saved = []

    def __iter__(self):
        return self

    def next(self):
        next = self.readline()
        if not next:
            raise StopIteration
        return next

    def readlines(self, *args, **keywds):
        lines = self._saved + self._handle.readlines(*args,**keywds)
        self._saved = []
        return lines

    def readline(self, *args, **keywds):
        if self._saved:
            line = self._saved.pop(0)
        else:
            line = self._handle.readline(*args,**keywds)
        return line

    def read(self, size=-1):
        if size == -1:
            saved = "".join(self._saved)
            self._saved[:] = []
        else:
            saved = ''
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
        if line:
            self._saved = [line] + self._saved

    def peekline(self):
        if self._saved:
            line = self._saved[0]
        else:
            line = self._handle.readline()
            self.saveline(line)
        return line

    def tell(self):
        lengths = map(len, self._saved)
        sum = reduce(lambda x, y: x+y, lengths, 0)
        return self._handle.tell() - sum

    def seek(self, *args):
        self._saved = []
        self._handle.seek(*args)

    def __getattr__(self, attr):
        return getattr(self._handle, attr)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self._handle.close()


# I could make this faster by using cStringIO.
# However, cStringIO (in v1.52) does not implement the
# readlines method.
class StringHandle(StringIO.StringIO):
    def __init__(self, buffer=''):
        import warnings
        import Bio
        warnings.warn("This class is deprecated, and is likely to be removed in a future version of Biopython. Please use the class StringIO in the module StringIO in the Python standard library instead", Bio.BiopythonDeprecationWarning)
        StringIO.StringIO.__init__(self, buffer)

try:
    import sgmllib
except ImportError:
    #This isn't available on Python 3, but we don't care much as SGMLStripper
    #is obsolete
    pass
else:
    class SGMLStripper(object):
        """Object to strip SGML tags (OBSOLETE)."""
        class MyParser(sgmllib.SGMLParser):
            def __init__(self):
                sgmllib.SGMLParser.__init__(self)
                self.data = ''
            def handle_data(self, data):
                self.data = self.data + data

        def __init__(self):
            import warnings
            import Bio
            warnings.warn("This class is deprecated, and is likely to be removed in a future version of Biopython", Bio.BiopythonDeprecationWarning)
            self._parser = SGMLStripper.MyParser()

        def strip(self, str):
            """S.strip(str) -> string

            Strip the SGML tags from str.

            """
            if not str:  # empty string, don't do anything.
                return ''
            # I need to make sure that I don't return an empty string if
            # the buffer is not empty.  This can happen if there's a newline
            # character embedded within a tag.  Thus, I'll first check to
            # see if the last character is a newline.  If it is, and it's stripped
            # away, I'll add it back.
            is_newline = str[-1] in ['\n', '\r']

            self._parser.data = ''    # clear the parser's data (don't reset)
            self._parser.feed(str)
            if self._parser.data:
                str = self._parser.data
            elif is_newline:
                str = '\n'
            else:
                str = ''
            return str


