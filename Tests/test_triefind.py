#!/usr/bin/env python

import StringIO
from operator import truth

from Bio import trie
from Bio import triefind

trieobj = trie.trie()

trieobj["hello"] = 5
trieobj["he"] = 7
trieobj["hej"] = 9
trieobj["foo"] = "bar"
trieobj["wor"] = "ld"

print triefind.match("hello world!", trieobj)    # "hello"
k = triefind.match_all("hello world!", trieobj)
k.sort()
print k     # ["he", "hello"]

k = triefind.find("hello world!", trieobj)
k.sort()
print k     # [("he", 0, 2), ("hello", 0, 5), ("wor", 6, 9)]

k = triefind.find_words("hello world!", trieobj)
k.sort()
print k     # [("hello", 0, 5)]

trieobj["world"] = "full"
k = triefind.find("hello world!", trieobj)
k.sort()
print k     # [("he", 0, 2), ("hello", 0, 5), ("wor", 6, 9), ("world", 6, 11)]

k = triefind.find_words("hello world!", trieobj)
k.sort()
print k     # [("hello", 0, 5), ("world", 6, 11)]
