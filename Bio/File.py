# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Code for more fancy file handles.


Classes:
UndoHandle     File object decorator with support for undo-like operations.
StringHandle   Wraps a file object around a string.

"""
import os
import string
import tempfile
import cStringIO

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

    def readlines(self, *args, **keywds):
        lines = self._saved + apply(self._handle.readlines, args, keywds)
        self._saved = []
        return lines

    def readline(self, *args, **keywds):
        if self._saved:
            line = self._saved.pop(0)
        else:
            line = apply(self._handle.readline, args, keywds)
        return line

    def read(self, *args, **keywds):
        saved = string.join(self._saved, '')
        return saved + apply(self._handle.read, args, keywds)

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
        apply(self._handle.seek, args)

    def __getattr__(self, attr):
        return getattr(self._handle, attr)

StringHandle = cStringIO.StringIO
