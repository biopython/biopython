#!/usr/bin/env python

# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import unittest

try:
    from Bio import trie
except ImportError:
    import os
    from Bio import MissingPythonDependencyError
    if os.name=="java":
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
        k = trieobj.keys()
        k.sort()
        self.assertEqual(k, ["foo", "he", "hej", "hello"])
        self.assertEqual(trieobj["hello"], 5)
        self.assertEqual(trieobj.get("bye"), None)
        self.assertEqual(trieobj.has_key("hello"), True)
        self.assertEqual(trieobj.has_key("he"), True)
        self.assertEqual(trieobj.has_key("bye"), False)
        self.assertEqual(trieobj.has_prefix("h"), True)
        self.assertEqual(trieobj.has_prefix("hel"), True)
        self.assertEqual(trieobj.has_prefix("foa"), False)
        self.assertEqual(trieobj.has_prefix("hello world"), False)
        self.assertEqual(len(trieobj), 4)
        k = trieobj.with_prefix("he")
        k.sort()
        self.assertEqual(k, ["he", "hej", "hello"])
        k = trieobj.with_prefix("l")
        self.assertEqual(k, [])
        k = trieobj.with_prefix("hej")
        self.assertEqual(k, ["hej"])
        k = trieobj.with_prefix("hejk")
        self.assertEqual(k, [])

    def test_save(self):
        import StringIO
        trieobj = trie.trie()
        trieobj["foo"] = 1
        k = trieobj.keys()
        self.assertEqual(k, ["foo"])
        v = trieobj.values()
        self.assertEqual(v, [1])
        self.assertEqual(trieobj.get("bar", 99), 99)
        trieobj["hello"] = '55a'
        self.assertEqual(trieobj.get_approximate("foo", 0), [("foo", 1, 0)])
        self.assertEqual(trieobj.get_approximate("foo", 1), [("foo", 1, 0)])
        self.assertEqual(trieobj.get_approximate("foa", 0), [])
        self.assertEqual(trieobj.get_approximate("foa", 1), [("foo", 1, 1)])
        x = trieobj.get_approximate("foa", 2)
        x.sort()
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
        x = y.items()
        x.sort()
        self.assertEqual(x,[(('foo', 1, 0), 1), (('hello', '55a', 4), 6)])
        h = StringIO.StringIO()
        trie.save(h, trieobj)
        h.seek(0)
        trieobj = trie.load(h)
        k = trieobj.keys()
        self.assertTrue("foo" in k)
        self.assertTrue("hello" in k)
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
        k = triefind.match_all("hello world!", trieobj)
        k.sort()
        self.assertEqual(k, ["he", "hello"])
        k = triefind.find("hello world!", trieobj)
        k.sort()
        self.assertEqual(k, [("he", 0, 2), ("hello", 0, 5), ("wor", 6, 9)])
        k = triefind.find_words("hello world!", trieobj)
        k.sort()
        self.assertEqual(k, [("hello", 0, 5)])
        trieobj["world"] = "full"
        k = triefind.find("hello world!", trieobj)
        k.sort()
        self.assertEqual(k, [("he", 0, 2), ("hello", 0, 5), ("wor", 6, 9), ("world", 6, 11)])
        k = triefind.find_words("hello world!", trieobj)
        k.sort()
        self.assertEqual(k, [("hello", 0, 5), ("world", 6, 11)])


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
