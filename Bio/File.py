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


# Old implementation of StringHandle based on pipes.
# This fails if str fills up the pipe.
# Maybe I should use a pipe if the string is short, but
# a "real" file if it's long?

#class StringHandle:
#    def __init__(self, str):
#        r, w = os.pipe()
#        os.fdopen(w, 'w').write(str)
#        self._handle = os.fdopen(r, 'r')
#
#    def __getattr__(self, name):
#        return getattr(self._handle, name)
#        
#    def __setattr__(self, name, value):
#        if name == '_handle':
#            self.__dict__[name] = value
#        else:
#            setattr(self._handle, name, value)

# XXX need to replace this with StringIO
class StringHandle:
    def __init__(self, str):
        self._filename = tempfile.mktemp()
        open(self._filename, 'w').write(str)
        self._handle = open(self._filename)

    def __del__(self):
        if self.__dict__.has_key('_handle'):
            self._handle.close()
        if os.path.exists(self._filename):
            os.unlink(self._filename)

    def __getattr__(self, name):
        return getattr(self._handle, name)
        
    def __setattr__(self, name, value):
        if name in ['_handle', '_filename']:
            self.__dict__[name] = value
        else:
            setattr(self._handle, name, value)
