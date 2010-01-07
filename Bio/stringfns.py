# Copyright 2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""This provides useful general functions for working with strings (DEPRECATED).

This module is considered to be deprecated, and is likely to be removed in a
future release of Biopython.  Its C code implementation has already been
removed. Please get in touch via the mailing list if this will affect you.

Functions:
splitany       Split a string using many delimiters.
find_anychar   Find one of a list of characters in a string.
rfind_anychar  Find one of a list of characters in a string, from end to start.

"""
import warnings
warnings.warn("Bio.stringfns and its C code equivalent Bio.cstringfns are" \
              +" deprecated, and will be removed in a future release of"\
              +" Biopython.  If you want to continue to use this code,"\
              +" please get in contact with the Biopython developers via"\
              +" the mailing lists to avoid its permanent removal from"\
              +" Biopython.", \
              DeprecationWarning)

def splitany(s, sep=" \011\012\013\014\015", maxsplit=None, negate=0):
    """splitany(s [,sep [,maxsplit [,negate]]]) -> list of strings

    Split a string.  Similar to string.split, except that this considers
    any one of the characters in sep to be a delimiter.  If negate is
    true, then everything but sep will be a separator.

    """
    strlist = []
    prev = 0
    for i in range(len(s)):
        if maxsplit is not None and len(strlist) >= maxsplit:
            break
        if (s[i] in sep) == (not negate):
            strlist.append(s[prev:i])
            prev = i+1
    strlist.append(s[prev:])
    return strlist

def find_anychar(string, chars, index=None, negate=0):
    """find_anychar(string, chars[, index]) -> index of a character or -1

    Find a character in string.  chars is a list of characters to look
    for.  Return the index of the first occurrence of any of the
    characters, or -1 if not found.  index is the index where the
    search should start.  By default, I search from the beginning of
    the string.
    
    """
    if index is None:
        index = 0
    while index < len(string) and \
          ((not negate and string[index] not in chars) or
           (negate and string[index] in chars)):
        index += 1
    if index == len(string):
        return -1
    return index

def rfind_anychar(string, chars, index=None, negate=0):
    """rfind_anychar(string, chars[, index]) -> index of a character or -1

    Find a character in string, looking from the end to the start.
    chars is a list of characters to look for.  Return the index of
    the first occurrence of any of the characters, or -1 if not found.
    index is the index where the search should start.  By default, I
    search from the end of the string.
    
    """
    if index is None:
        index = len(string)-1
    while index >= 0 and \
          ((not negate and string[index] not in chars) or
           (negate and string[index] in chars)):
        index -= 1
    # If not found, index will already be -1.
    return index
