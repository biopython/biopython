# Copyright 2010-2013 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Python 3 compatibility tools (PRIVATE).

We currently have lines like this under Python 2 in order
to use iterator based zip, map and filter:

    from future_builtins import zip

There is no similar option for range yet, other than:

    range = xrange
    input = raw_input

or:

    from __builtin__ import xrange as range
    from __builtin__ import raw_input as input

Under Python 3 this imports need to be removed. Also, deliberate
importing of built in functions like open changes from Python 2:

    from __builtin__ import open

to this under Python 3:

    from builtins import open

Instead, we can do this under either Python 2 or 3:

    from Bio._py3k import open
    from Bio._py3k import zip

Once we drop support for Python 2, the whole of Bio._py3k will
go away.
"""
import sys

if sys.version_info[0] >= 3:
    #Code for Python 3
    from builtins import open, zip, map, filter, range, input

    import codecs

    #Lots of our Python 2 code uses isinstance(x, basestring)
    #which after 2to3 becomes isinstance(x, str)
    basestring = str
    unicode = str

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
                return self._handle.tell()

        return EvilHandleHack(handle)

    #On Python 3, can depend on OrderedDict being present:
    from collections import OrderedDict

    #On Python 3, this will be a unicode StringIO
    from io import StringIO

    #On Python 3 urllib, urllib2, and urlparse were merged:
    from urllib.request import urlopen, Request, urlretrieve, urlparse
    from urllib.parse import urlencode, quote
    from urllib.error import HTTPError

else:
    #Python 2 code
    from __builtin__ import open, basestring, unicode

    #Import Python3 like iterator functions:
    from future_builtins import zip, map, filter
    from __builtin__ import xrange as range
    from __builtin__ import raw_input as input

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
        return isinstance(i, (int, long))

    def _binary_to_string_handle(handle):
        """Treat a binary handle like a text handle."""
        return handle

    try:
        #Present on Python 2.7
        from collections import OrderedDict
    except ImportError:
        try:
            #Raymond Hettinger's backport available on PyPI
            from ordereddict import OrderedDict
        except ImportError:
            #Use our bundled copy instead
            from ._ordereddict import OrderedDict

    # On Python 2 this will be a (bytes) string based handle.
    # Note this doesn't work as it is unicode based:
    # from io import StringIO
    try:
        from cStringIO import StringIO
    except ImportError:
        from StringIO import StringIO

    #Under urllib.request on Python 3:
    from urllib2 import urlopen, Request
    from urllib import urlretrieve
    from urlparse import urlparse

    #Under urllib.parse on Python 3:
    from urllib import urlencode, quote

    #Under urllib.error on Python 3:
    from urllib2 import HTTPError


if sys.platform == "win32":
    # Can't use commands.getoutput on Python 2, Unix only/broken:
    # http://bugs.python.org/issue15073
    # Can't use subprocess.getoutput on Python 3, Unix only/broken:
    # http://bugs.python.org/issue10197
    def getoutput(cmd):
        import subprocess
        child = subprocess.Popen(cmd,
                                 stdin=subprocess.PIPE,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 universal_newlines=True,
                                 shell=False)
        stdout, stderr = child.communicate()
        # Remove trailing \n to match the Unix function,
        return stdout.rstrip("\n")
elif sys.version_info[0] >= 3:
    # Use subprocess.getoutput on Python 3,
    from subprocess import getoutput
else:
    # Use commands.getoutput on Python 2,
    from commands import getoutput
