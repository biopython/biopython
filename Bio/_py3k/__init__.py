# Copyright 2010-2018 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Python 3 compatibility tools (PRIVATE).

Once we drop support for Python 2, the whole of Bio._py3k will
go away.
"""

# From the point of view of pep8 and flake8, there are lots of issues with
# this file. This line tells flake8 to ignore it for quality assurance:
# flake8: noqa

import sys

import codecs


def _bytes_bytearray_to_str(s):
    """If s is bytes or bytearray, convert to a unicode string (PRIVATE)."""
    if isinstance(s, (bytes, bytearray)):
        return s.decode()
    return s


def _as_unicode(s):
    """Turn byte string or unicode string into a unicode string (PRIVATE)."""
    if isinstance(s, str):
        return s
    # Assume it is a bytes string
    # Note ISO-8859-1 aka Latin-1 preserves first 256 chars
    return codecs.latin_1_decode(s)[0]


def _as_bytes(s):
    """Turn byte string or unicode string into a bytes string (PRIVATE).

    The Python 2 version returns a (byte) string.
    """
    if isinstance(s, bytes):
        return s
    # Assume it is a unicode string
    # Note ISO-8859-1 aka Latin-1 preserves first 256 chars
    return codecs.latin_1_encode(s)[0]


_as_string = _as_unicode


def _is_int_or_long(i):
    """Check if the value is an integer (PRIVATE).

    Note there are no longs on Python 3.
    """
    return isinstance(i, int)


import io
import locale

# Python 3.4 onwards, the standard library wrappers should work:
def _binary_to_string_handle(handle):
    """Treat a binary (bytes) handle like a text (unicode) handle (PRIVATE)."""
    try:
        # If this is a network handle from urllib,
        # the HTTP headers may tell us the encoding.
        encoding = handle.headers.get_content_charset()
    except AttributeError:
        encoding = None
    if encoding is None:
        # The W3C recommendation is:
        # When no explicit charset parameter is provided by the sender,
        # media subtypes of the "text" type are defined to have a default
        # charset value of "ISO-8859-1" when received via HTTP.
        # "ISO-8859-1" is also known as 'latin-1'
        # See the following for more detail:
        # https://www.w3.org/Protocols/rfc2616/rfc2616-sec3.html#sec3.7.1
        encoding = "latin-1"
    wrapped = io.TextIOWrapper(io.BufferedReader(handle), encoding=encoding)
    try:
        # If wrapping an online handle, this is nice to have:
        wrapped.url = handle.url
    except AttributeError:
        pass
    return wrapped


# This is to avoid the deprecation warning from open(filename, "rU")
_universal_read_mode = "r"  # text mode does universal new lines

# On Python 3 urllib, urllib2, and urlparse were merged:
from urllib.request import urlopen, Request, urlparse, urlcleanup
from urllib.parse import urlencode, quote
from urllib.error import URLError, HTTPError
