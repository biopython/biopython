# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Code for more fancy file handles.


Classes:
UndoHandle     File object decorator with support for undo-like operations.
StringHandle   Wraps a file object around a string.
"""

import urllib  # XXX remove this
import os

class UndoHandle:
    """A Python handle that adds functionality for saving lines.

    Saves lines in a LIFO fashion.

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
        lines = self._saved + apply(self._handle.readlines, args, keywds)
        self._saved = []
        return lines

    def readline(self, *args, **keywds):
        if self._saved:
            line = self._saved.pop(0)
        else:
            line = apply(self._handle.readline, args, keywds)
        return line

    def saveline(self, line):
        self._saved = [line] + self._saved

    def peekline(self):
        if self._saved:
            line = self._saved[0]
        else:
            line = self._handle.readline()
            self.saveline(line)
        return line

    def __getattr__(self, attr):
        return getattr(self._handle, attr)

class StringHandle:
    def __init__(self, str):
        r, w = os.pipe()
        os.fdopen(w, 'w').write(str)
        self._handle = os.fdopen(r, 'r')

    def __getattr__(self, name):
        if name == '_handle':
            return self.__dict__[name]
        return getattr(self._handle, name)
        
    def __setattr__(self, name, value):
        if name == '_handle':
            self.__dict__[name] = value
        else:
            setattr(self._handle, name, value)




### THESE ARE DEPRECATED

# XXX to be removed
fileopen = open
def _open(*args):
    """open(*args) -> UndoHandle

    Open a file.  The arguments are the same as the standard open function.
    
    """
    return UndoHandle(apply(fileopen, args))

# XXX to be removed
def _urlopen(*args):
    """urlopen(*args) -> UndoHandle

    Open a URL.  The arguments are the same as urllib.urlopen.

    """
    return UndoHandle(apply(urllib.urlopen, args))
