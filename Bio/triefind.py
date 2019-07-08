# Copyright 2002-2003 Jeff Chang.  All rights reserved.
# Revisions copyright 2012 by Christian Brueffer.  All rights reserved.
# Revisions copyright 2012-2017 by Peter Cock.  All rights reserved.
# Revisions copyright 2015 by Brian Osborne.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
#

"""Given a trie, find all occurrences of a word in the trie in a string.

Like searching a string for a substring, except that the substring is
any word in a trie.

Functions:
 - match         Find longest key in a trie matching the beginning of the string.
 - match_all     Find all keys in a trie matching the beginning of the string.
 - find          Find keys in a trie matching anywhere in a string.
 - find_words    Find keys in a trie matching whole words in a string.

This module is DEPRECATED. We encourage users to switch to alternative libraries
implementing a trie data structure, for example pygtrie.
"""

import string
import re

from Bio import BiopythonDeprecationWarning
import warnings
warnings.warn("This module has been deprecated. We encourage users to switch "
              "to alternative libraries implementing a trie data structure, "
              "for example pygtrie.", BiopythonDeprecationWarning)


def match(string, trie):
    """Find longest key, or return None.

    Find the longest key in the trie that matches the beginning of the
    string.
    """
    longest = None
    for i in range(len(string)):
        substr = string[:i + 1]
        if not trie.has_prefix(substr):
            break
        if substr in trie:
            longest = substr
    return longest


def match_all(string, trie):
    """Find and return a list of keys.

    Find all the keys in the trie that matches the beginning of the
    string.
    """
    matches = []
    for i in range(len(string)):
        substr = string[:i + 1]
        if not trie.has_prefix(substr):
            break
        if substr in trie:
            matches.append(substr)
    return matches


def find(string, trie):
    """Find all the keys in the trie that match anywhere in the string.

    Returns a list of tuples (key, start, end).
    """
    results = []
    start = 0     # index to start the search
    while start < len(string):
        # Look for a match.
        keys = match_all(string[start:], trie)
        for key in keys:
            results.append((key, start, start + len(key)))
        start += 1
    return results


DEFAULT_BOUNDARY_CHARS = string.punctuation + string.whitespace


def find_words(string, trie):
    """Find all the keys in the trie that match full words in the string.

    Find all the keys in the trie that match full words in the string.
    Word boundaries are defined as any punctuation or whitespace.

    Returns a list of tuples (key, start, end).
    """
    _boundary_re = re.compile(r"[%s]+" % re.escape(DEFAULT_BOUNDARY_CHARS))

    results = []
    start = 0     # index of word boundary
    while start < len(string):
        # Look for a match.
        keys = match_all(string[start:], trie)
        for key in keys:
            length = len(key)
            # Make sure it ends at a boundary.
            if start + length == len(string) or \
               _boundary_re.match(string[start + length]):
                results.append((key, start, start + length))
        # Move forward to the next boundary.
        m = _boundary_re.search(string, start)
        if m is None:
            break
        start = m.end()
    return results
