#!/usr/bin/env python

import StringIO
from operator import truth

try :
    from Bio import trie
except ImportError :
    import os
    from Bio import MissingExternalDependencyError
    if os.name=="java" :
        message = "Not available on Jython, Bio.trie requires compiled C code."
    else :
        message = "Could not import Bio.trie, check C code was compiled."
    raise MissingExternalDependencyError(message)

trieobj = trie.trie()

trieobj["hello"] = 5
trieobj["he"] = 7
trieobj["hej"] = 9
trieobj["foo"] = "bar"

k = trieobj.keys()
k.sort()
print k                          # ["foo", "he", "hej", "hello"]
print trieobj["hello"]           # 5
print trieobj.get("bye")         # None

print trieobj.has_key("hello")   # 1
print trieobj.has_key("he")      # 1
print trieobj.has_key("bye")     # 0

print trieobj.has_prefix("h")    # 1
print trieobj.has_prefix("hel")  # 1
print trieobj.has_prefix("foa")  # 0
print trieobj.has_prefix("hello world")   # 0

print len(trieobj)               # 4

k = trieobj.with_prefix("he")
k.sort()
print k                          # ["he", "hej", "hello"]
k = trieobj.with_prefix("l")
k.sort()
print k                          # []
k = trieobj.with_prefix("hej")
k.sort()
print k                          # ["hej"]
k = trieobj.with_prefix("hejk")
k.sort()
print k                          # []

trieobj2 = trie.trie()
trieobj2["foo"] = 1
k = trieobj2.keys()
k.sort()
print k                          # ["foo"]
v = trieobj2.values()
v.sort()
print v                          # [1]

print trieobj2.get("bar", 99)    # 99

trieobj2["hello"] = '55a'

print trieobj2.get_approximate("foo", 0)    # [("foo", 1, 0)]
print trieobj2.get_approximate("foo", 1)    # [("foo", 1, 0)]
print trieobj2.get_approximate("foa", 0)    # []
print trieobj2.get_approximate("foa", 1)    # [("foo", 1, 1)]
x = trieobj2.get_approximate("foa", 2)
print "found %d matches" % len(x)           # 3
x.sort()
print x                       # [("foo", 1, 1), ("foo", 1, 2), ("foo", 1, 2)]
# foo  foo-  foo-
# foa  f-oa  fo-a


#import sys; sys.exit(0)

# mismatch a->o
# insertion after f, deletion of o
# insertion after o, deletion of o

x = trieobj2.get_approximate("foo", 4)
y = {}
for z in x:
    y[z] = y.get(z, 0) + 1
x = y.items()
x.sort()
print x                       # [(('foo', 1, 0), 1), (('hello', '55a', 4), 6)]

h = StringIO.StringIO()
trie.save(h, trieobj2)
h.seek(0)
trieobj3 = trie.load(h)
k = trieobj3.keys()
k.sort()
for m in k:                       # foo 1
    print m, repr(trieobj3[m])    # hello '55a'


# Found bug, doesn't handle insertions and deletions at end properly.
trieobj = trie.trie()
trieobj["hello"] = 1
print trieobj.get_approximate('he', 2)        # []
print trieobj.get_approximate('he', 3)        # [('hello', 1, 3)]
print trieobj.get_approximate('hello me!', 3) # []
print trieobj.get_approximate('hello me!', 4) # [('hello', 1, 4)]
print trieobj.get_approximate('hello me!', 5) # [('hello', 1, 4)]
