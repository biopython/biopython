"""
Given a trie, find all occurrences of a word in the trie in a string.

Like searching a string for a substring, except that the substring is
any word in a trie.

Functions:
match         Find longest key in a trie matching the beginning of the string.
match_all     Find all keys in a trie matching the beginning of the string.
find          Find keys in a trie matching anywhere in a string.
find_words    Find keys in a trie matching whole words in a string.

"""
import string
import re

def match(string, trie):
    """match(string, trie) -> longest key or None

    Find the longest key in the trie that matches the beginning of the
    string.

    """
    longest = None
    for i in range(len(string)):
        substr = string[:i+1]
        if not trie.has_prefix(substr):
            break
        if trie.has_key(substr):
            longest = substr
    return longest

def match_all(string, trie):
    """match_all(string, trie) -> list of keys

    Find all the keys in the trie that matches the beginning of the
    string.

    """
    matches = []
    for i in range(len(string)):
        substr = string[:i+1]
        if not trie.has_prefix(substr):
            break
        if trie.has_key(substr):
            matches.append(substr)
    return matches

def find(string, trie):
    """find(string, trie) -> list of tuples (key, start, end)

    Find all the keys in the trie that match anywhere in the string.

    """
    results = []
    start = 0     # index to start the search
    while start < len(string):
        # Look for a match.
        keys = match_all(string[start:], trie)
        for key in keys:
            results.append((key, start, start+len(key)))
        start += 1
    return results

DEFAULT_BOUNDARY_CHARS = string.punctuation + string.whitespace

def find_words(string, trie):
    """find_words(string, trie) -> list of tuples (key, start, end)

    Find all the keys in the trie that match full words in the string.
    Word boundaries are defined as any punctuation or whitespace.

    """
    _boundary_re = re.compile(r"[%s]+" % re.escape(DEFAULT_BOUNDARY_CHARS))
        
    results = []
    start = 0     # index of word boundary
    while start < len(string):
        # Look for a match.
        keys = match_all(string[start:], trie)
        for key in keys:
            l = len(key)
            # Make sure it ends at a boundary.
            if start+l == len(string) or \
               _boundary_re.match(string[start+l]):
                results.append((key, start, start+l))
        # Move forward to the next boundary.
        m = _boundary_re.search(string, start)
        if m is None:
            break
        start = m.end()
    return results
