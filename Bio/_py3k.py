# Copyright 2010 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Python 3 compatibility tools (PRIVATE)."""

import sys

if sys.version_info[0] >= 3:
    #Python 3 code (which will be converted using 2to3 script)

    _bytes_to_string = lambda b: b.decode() # bytes to unicode string
    _string_to_bytes = lambda s: s.encode() # unicode string to bytes

    def _as_unicode(s):
        """Turn byte string or unicode string into a unicode string."""
        if isinstance(s, str):
            return s
        #Assume it is a bytes string
        return s.decode()


    def _as_bytes(s):
        """Turn byte string or unicode string into a bytes string.
        
        The Python 2 version returns a (byte) string.
        """
        if isinstance(s, bytes):
            return s
        #Assume it is a unicode string
        return s.encode()
    
    _as_string = _as_unicode

    def _is_int_or_long(i):
        """Check if the value is an integer.

        Note there are no longs on Python 3.
        """
        return isinstance(i, int)

else:
    #Python 2 code

    _bytes_to_string = lambda b: b # bytes to string, i.e. do nothing
    _string_to_bytes = lambda s: str(s) # str (or unicode) to bytes string

    def _as_unicode(s):
        """Turn a (byte) string or a unicode string into a (byte) string."""
        #Will be changed by 2to3 to "isinstance(s, str)" but doesn't matter:
        if isinstance(s, unicode):
            return s
        return s.decode()
    
    def _as_bytes(s):
        """Turn a (byte) string or a unicode string into a (byte) string."""
        return str(s)
    
    _as_string = _as_bytes

    def _is_int_or_long(i):
        """Check if the value is an integer or long."""
        #If the 2to3 long fixer is enabled (which it is by default), this
        #will be changed to "isinstance(i, int) or isinstance(i, int)"
        #but that doesn't matter.
        return isinstance(i, int) or isinstance(i, long)
