# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Code for more fancy file handles.


Classes:

UndoHandle     File object decorator with support for undo-like operations.

StringHandle   Wraps a file object around a string.

SGMLHandle     File object that automatically strips SGML tags from data
               (DEPRECATED).

SGMLStripper   Object that strips SGML.  This is now considered OBSOLETE, and
               is likely to be deprecated in a future release of Biopython,
               and later removed.

"""
import os
import StringIO
import sgmllib

class UndoHandle:
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
StringHandle = StringIO.StringIO

class SGMLHandle:
    """A Python handle that automatically strips SGML tags from data (DEPRECATED).

    This module is now considered to be obsolete, and is likely to be
    deprecated in a future release of Biopython, and later removed.
    """
    def __init__(self, handle):
        """SGMLStripper(handle)

        handle is a file handle to SGML-formatted data.
        
        """
        import warnings
        warnings.warn("Bio.File.SGMLHandle is deprecated, and will be removed"\
                      " in a future release of Biopython.  If you want to"\
                      " continue to use this code, please get in contact via"\
                      " the mailing lists to avoid its permanent removal from"\
                      " Biopython.", DeprecationWarning)
        self._handle = handle
        self._stripper = SGMLStripper()

    def read(self, *args, **keywds):
        data = self._handle.read(*args, **keywds)
        return self._stripper.strip(data)

    def readline(self, *args, **keywds):
        line = self._handle.readline(*args, **keywds)
        return self._stripper.strip(line)

    def readlines(self, *args, **keywds):
        lines = self._handle.readlines(*args, **keywds)
        for i in range(len(lines)):
            lines[i] = self._stripper.strip(str)
        return lines

    def __getattr__(self, attr):
        return getattr(self._handle, attr)


class SGMLStripper:
    class MyParser(sgmllib.SGMLParser):
        def __init__(self):
            sgmllib.SGMLParser.__init__(self)
            self.data = ''
        def handle_data(self, data):
            self.data = self.data + data

    def __init__(self):
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

