#!/usr/bin/env python

# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import unittest
from io import BytesIO

try:
    from Bio import trie
except ImportError:
    import os
    from Bio import MissingPythonDependencyError
    if os.name == "java":
        message = "Not available on Jython, Bio.trie requires compiled C code."
    else:
        message = "Could not import Bio.trie, check C code was compiled."
    raise MissingPythonDependencyError(message)


class TestTrie(unittest.TestCase):

    def test_get_set(self):
        trieobj = trie.trie()
        trieobj["hello world"] = "s1"
        trieobj["bye"] = "s2"
        trieobj["hell sucks"] = "s3"
        trieobj["hebee"] = "s4"
        self.assertEqual(trieobj["hello world"], "s1")
        self.assertEqual(trieobj["bye"], "s2")
        self.assertEqual(trieobj["hell sucks"], "s3")
        self.assertEqual(trieobj["hebee"], "s4")
        trieobj["blah"] = "s5"
        self.assertEqual(trieobj["blah"], "s5")
        self.assertEqual(trieobj.get("foobar"), None)
        self.assertEqual(len(trieobj), 5)
        trieobj["blah"] = "snew"
        self.assertEqual(trieobj["blah"], "snew")

    def test_prefix(self):
        trieobj = trie.trie()
        trieobj["hello"] = 5
        trieobj["he"] = 7
        trieobj["hej"] = 9
        trieobj["foo"] = "bar"
        k = sorted(trieobj.keys())
        self.assertEqual(k, ["foo", "he", "hej", "hello"])
        self.assertEqual(trieobj["hello"], 5)
        self.assertEqual(trieobj.get("bye"), None)
        self.assertIn("hello", trieobj)
        self.assertIn("he", trieobj)
        self.assertFalse("bye" in trieobj)
        self.assertTrue(trieobj.has_prefix("h"))
        self.assertTrue(trieobj.has_prefix("hel"))
        self.assertFalse(trieobj.has_prefix("foa"))
        self.assertFalse(trieobj.has_prefix("hello world"))
        self.assertEqual(len(trieobj), 4)
        k = sorted(trieobj.with_prefix("he"))
        self.assertEqual(k, ["he", "hej", "hello"])
        k = trieobj.with_prefix("l")
        self.assertEqual(k, [])
        k = trieobj.with_prefix("hej")
        self.assertEqual(k, ["hej"])
        k = trieobj.with_prefix("hejk")
        self.assertEqual(k, [])

    def test_save(self):
        trieobj = trie.trie()
        trieobj["foo"] = 1
        k = list(trieobj.keys())
        self.assertEqual(k, ["foo"])
        v = list(trieobj.values())
        self.assertEqual(v, [1])
        self.assertEqual(trieobj.get("bar", 99), 99)
        trieobj["hello"] = '55a'
        self.assertEqual(trieobj.get_approximate("foo", 0), [("foo", 1, 0)])
        self.assertEqual(trieobj.get_approximate("foo", 1), [("foo", 1, 0)])
        self.assertEqual(trieobj.get_approximate("foa", 0), [])
        self.assertEqual(trieobj.get_approximate("foa", 1), [("foo", 1, 1)])
        x = sorted(trieobj.get_approximate("foa", 2))
        self.assertEqual(x, [("foo", 1, 1), ("foo", 1, 2), ("foo", 1, 2)])
        # foo  foo-  foo-
        # foa  f-oa  fo-a
        # mismatch a->o
        # insertion after f, deletion of o
        # insertion after o, deletion of o
        x = trieobj.get_approximate("foo", 4)
        y = {}
        for z in x:
            y[z] = y.get(z, 0) + 1
        x = sorted(y.items())
        self.assertEqual(x, [(('foo', 1, 0), 1), (('hello', '55a', 4), 6)])
        h = BytesIO()
        trie.save(h, trieobj)
        h.seek(0)
        trieobj = trie.load(h)
        k = list(trieobj.keys())
        self.assertIn("foo", k)
        self.assertIn("hello", k)
        self.assertEqual(repr(trieobj["foo"]), '1')
        self.assertEqual(repr(trieobj["hello"]), "'55a'")

    def test_get_approximate(self):
        # Found bug, doesn't handle insertions and deletions at end properly.
        trieobj = trie.trie()
        trieobj["hello"] = 1
        self.assertEqual(trieobj.get_approximate('he', 2), [])
        self.assertEqual(trieobj.get_approximate('he', 3), [('hello', 1, 3)])
        self.assertEqual(trieobj.get_approximate('hello me!', 3), [])
        self.assertEqual(trieobj.get_approximate('hello me!', 4), [('hello', 1, 4)])
        self.assertEqual(trieobj.get_approximate('hello me!', 5), [('hello', 1, 4)])

    def test_with_prefix(self):
        trieobj = trie.trie()
        s = "BANANA"
        for i in range(len(s)):  # insert all suffixes into trie
            trieobj[s[i:]] = i
            self.assertEqual(trieobj[s[i:]], i)
        self.assertEqual(set(trieobj.values()), set(range(6)))
        self.assertEqual(set(['A', 'ANA', 'ANANA', 'BANANA', 'NA', 'NANA']),
                         set(trieobj.keys()))
        self.assertEqual(set(['NA', 'NANA']),
                         set(trieobj.with_prefix("N")))
        self.assertEqual(set(['NA', 'NANA']),
                         set(trieobj.with_prefix("NA")))
        self.assertEqual(set(['A', 'ANA', 'ANANA']),
                         set(trieobj.with_prefix("A")))
        self.assertEqual(set(['ANA', 'ANANA']),
                         set(trieobj.with_prefix("AN")))


class TestTrieFind(unittest.TestCase):

    def test_find(self):
        from Bio import triefind
        trieobj = trie.trie()
        trieobj["hello"] = 5
        trieobj["he"] = 7
        trieobj["hej"] = 9
        trieobj["foo"] = "bar"
        trieobj["wor"] = "ld"
        self.assertEqual(triefind.match("hello world!", trieobj), "hello")
        k = sorted(triefind.match_all("hello world!", trieobj))
        self.assertEqual(k, ["he", "hello"])
        k = sorted(triefind.find("hello world!", trieobj))
        self.assertEqual(k, [("he", 0, 2), ("hello", 0, 5), ("wor", 6, 9)])
        k = sorted(triefind.find_words("hello world!", trieobj))
        self.assertEqual(k, [("hello", 0, 5)])
        trieobj["world"] = "full"
        k = sorted(triefind.find("hello world!", trieobj))
        self.assertEqual(k, [("he", 0, 2), ("hello", 0, 5), ("wor", 6, 9), ("world", 6, 11)])
        k = sorted(triefind.find_words("hello world!", trieobj))
        self.assertEqual(k, [("hello", 0, 5), ("world", 6, 11)])


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
