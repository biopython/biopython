# Copyright 2010 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Python 3 compatibility tools (PRIVATE)."""

import sys

if sys.version_info[0] >= 3:
    #Python 3 code (which will be converted using 2to3 script)
    import codecs

    _bytes_to_string = lambda b: b.decode() # bytes to unicode string
    _string_to_bytes = lambda s: s.encode() # unicode string to bytes

    def _as_unicode(s):
        """Turn byte string or unicode string into a unicode string."""
        if isinstance(s, str):
            return s
        #Assume it is a bytes string
        #Note ISO-8859-1 aka Latin-1 preserves first 256 chars
        return codecs.latin_1_decode(s)[0]


    def _as_bytes(s):
        """Turn byte string or unicode string into a bytes string.
        
        The Python 2 version returns a (byte) string.
        """
        if isinstance(s, bytes):
            return s
        #Assume it is a unicode string
        #Note ISO-8859-1 aka Latin-1 preserves first 256 chars
        return codecs.latin_1_encode(s)[0]

    
    _as_string = _as_unicode

    def _is_int_or_long(i):
        """Check if the value is an integer.

        Note there are no longs on Python 3.
        """
        return isinstance(i, int)

    import io
    def _binary_to_string_handle(handle):
        """Treat a binary (bytes) handle like a text (unicode) handle."""
        #See also http://bugs.python.org/issue5628
        #and http://bugs.python.org/issue13541
        #and http://bugs.python.org/issue13464 which should be fixed in Python 3.3
        #return io.TextIOWrapper(io.BufferedReader(handle))
        #TODO - Re-evaluate this workaround under Python 3.3
        #(perhaps we will only need it on Python 3.1 and 3.2?)
        class EvilHandleHack(object):
            def __init__(self, handle):
                self._handle = handle
            def read(self, length=None):
                return _as_string(self._handle.read(length))
            def readline(self):
                return _as_string(self._handle.readline())
            def __iter__(self):
                for line in self._handle:
                    yield _as_string(line)
            def close(self):
                return self._handle.close()
            def seek(self, pos):
                return self._handle.seek(pos)
            def tell(self):
                return self._handle.tell(pos)
        return EvilHandleHack(handle)

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

    def _binary_to_string_handle(handle):
        """Treat a binary handle like a text handle."""
        return handle


###########################################################
# Custom OrderedDict, to enable Python <2.7 compatibility #
###########################################################

# Copyright (c) 2009 Raymond Hettinger
#
# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation files
# (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
#     The above copyright notice and this permission notice shall be
#     included in all copies or substantial portions of the Software.
#
#     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
#     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
#     OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#     NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
#     WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#     FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
#     OTHER DEALINGS IN THE SOFTWARE.


try:

    from collections import OrderedDict

except ImportError:

    from UserDict import DictMixin

    class OrderedDict(dict, DictMixin):

        def __init__(self, *args, **kwds):
            if len(args) > 1:
                raise TypeError('expected at most 1 arguments, got %d' % len(args))
            try:
                self.__end
            except AttributeError:
                self.clear()
            self.update(*args, **kwds)

        def clear(self):
            self.__end = end = []
            end += [None, end, end]         # sentinel node for doubly linked list
            self.__map = {}                 # key --> [key, prev, next]
            dict.clear(self)

        def __setitem__(self, key, value):
            if key not in self:
                end = self.__end
                curr = end[1]
                curr[2] = end[1] = self.__map[key] = [key, curr, end]
            dict.__setitem__(self, key, value)

        def __delitem__(self, key):
            dict.__delitem__(self, key)
            key, prev, next = self.__map.pop(key)
            prev[2] = next
            next[1] = prev

        def __iter__(self):
            end = self.__end
            curr = end[2]
            while curr is not end:
                yield curr[0]
                curr = curr[2]

        def __reversed__(self):
            end = self.__end
            curr = end[1]
            while curr is not end:
                yield curr[0]
                curr = curr[1]

        def popitem(self, last=True):
            if not self:
                raise KeyError('dictionary is empty')
            if last:
                key = reversed(self).next()
            else:
                key = iter(self).next()
            value = self.pop(key)
            return key, value

        def __reduce__(self):
            items = [[k, self[k]] for k in self]
            tmp = self.__map, self.__end
            del self.__map, self.__end
            inst_dict = vars(self).copy()
            self.__map, self.__end = tmp
            if inst_dict:
                return (self.__class__, (items,), inst_dict)
            return self.__class__, (items,)

        def keys(self):
            return list(self)

        setdefault = DictMixin.setdefault
        update = DictMixin.update
        pop = DictMixin.pop
        values = DictMixin.values
        items = DictMixin.items
        iterkeys = DictMixin.iterkeys
        itervalues = DictMixin.itervalues
        iteritems = DictMixin.iteritems

        def __repr__(self):
            if not self:
                return '%s()' % (self.__class__.__name__,)
            return '%s(%r)' % (self.__class__.__name__, self.items())

        def copy(self):
            return self.__class__(self)

        @classmethod
        def fromkeys(cls, iterable, value=None):
            d = cls()
            for key in iterable:
                d[key] = value
            return d

        def __eq__(self, other):
            if isinstance(other, OrderedDict):
                if len(self) != len(other):
                    return False
                for p, q in  zip(self.items(), other.items()):
                    if p != q:
                        return False
                return True
            return dict.__eq__(self, other)

        def __ne__(self, other):
            return not self == other
