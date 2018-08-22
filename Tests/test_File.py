# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from __future__ import print_function

import os.path
import shutil
import sys
import tempfile
import unittest

from Bio import bgzf
from Bio import File
from Bio._py3k import StringIO

data = """This
is
a multi-line
file"""


class UndoHandleTests(unittest.TestCase):

    def test_one(self):
        h = File.UndoHandle(StringIO(data))
        self.assertEqual(h.readline(), "This\n")
        self.assertEqual(h.peekline(), "is\n")
        self.assertEqual(h.readline(), "is\n")
        # TODO - Meaning of saveline lacking \n?
        h.saveline("saved\n")
        self.assertEqual(h.peekline(), "saved\n")
        h.saveline("another\n")
        self.assertEqual(h.readline(), "another\n")
        self.assertEqual(h.readline(), "saved\n")
        # Test readlines after saveline
        h.saveline("saved again\n")
        lines = h.readlines()
        self.assertEqual(len(lines), 3)
        self.assertEqual(lines[0], "saved again\n")
        self.assertEqual(lines[1], "a multi-line\n")
        self.assertEqual(lines[2], "file")  # no trailing \n
        # should be empty now
        self.assertEqual(h.readline(), "")
        h.saveline("save after empty\n")
        self.assertEqual(h.readline(), "save after empty\n")
        self.assertEqual(h.readline(), "")

    def test_read(self):
        """Test read method"""
        h = File.UndoHandle(StringIO("some text"))
        h.saveline("more text")
        self.assertEqual(h.read(), 'more textsome text')

    def test_undohandle_read_block(self):
        for block in [1, 2, 10]:
            s = StringIO(data)
            h = File.UndoHandle(s)
            h.peekline()
            new = ""
            while True:
                tmp = h.read(block)
                if not tmp:
                    break
                new += tmp
            self.assertEqual(data, new)
            h.close()


class RandomAccess(unittest.TestCase):

    def test_plain(self):
        with File._open_for_random_access("Quality/example.fastq") as handle:
            self.assertTrue("r" in handle.mode)
            self.assertTrue("b" in handle.mode)

    def test_bgzf(self):
        with File._open_for_random_access("Quality/example.fastq.bgz") as handle:
            self.assertIsInstance(handle, bgzf.BgzfReader)

    def test_gzip(self):
        self.assertRaises(ValueError,
                          File._open_for_random_access,
                          "Quality/example.fastq.gz")


class AsHandleTestCase(unittest.TestCase):

    def setUp(self):
        # Create a directory to work in
        self.temp_dir = tempfile.mkdtemp(prefix='biopython-test')

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def _path(self, *args):
        return os.path.join(self.temp_dir, *args)

    def test_handle(self):
        "Test as_handle with a file-like object argument"
        p = self._path('test_file.fasta')
        with open(p, 'wb') as fp:
            with File.as_handle(fp) as handle:
                self.assertEqual(fp, handle, "as_handle should "
                                 "return argument when given a "
                                 "file-like object")
                self.assertFalse(handle.closed)

            self.assertFalse(handle.closed,
                             "Exiting as_handle given a file-like object "
                             "should not close the file")

    def test_string_path(self):
        "Test as_handle with a string path argument"
        p = self._path('test_file.fasta')
        mode = 'wb'
        with File.as_handle(p, mode=mode) as handle:
            self.assertEqual(p, handle.name)
            self.assertEqual(mode, handle.mode)
            self.assertFalse(handle.closed)
        self.assertTrue(handle.closed)

    @unittest.skipIf(
        sys.version_info < (3, 6),
        'Passing Path objects to File.as_handle requires Python >= 3.6',
    )
    def test_path_object(self):
        "Test as_handle with a pathlib.Path object"
        from pathlib import Path
        p = Path(self._path('test_file.fasta'))
        mode = 'wb'
        with File.as_handle(p, mode=mode) as handle:
            self.assertEqual(str(p.absolute()), handle.name)
            self.assertEqual(mode, handle.mode)
            self.assertFalse(handle.closed)
        self.assertTrue(handle.closed)

    @unittest.skipIf(
        sys.version_info < (3, 6),
        'Passing path-like objects to File.as_handle requires Python >= 3.6',
    )
    def test_custom_path_like_object(self):
        "Test as_handle with a custom path-like object"
        class CustomPathLike:
            def __init__(self, path):
                self.path = path

            def __fspath__(self):
                return self.path

        p = CustomPathLike(self._path('test_file.fasta'))
        mode = 'wb'
        with File.as_handle(p, mode=mode) as handle:
            self.assertEqual(p.path, handle.name)
            self.assertEqual(mode, handle.mode)
            self.assertFalse(handle.closed)
        self.assertTrue(handle.closed)

    def test_stringio(self):
        s = StringIO()
        with File.as_handle(s) as handle:
            self.assertIs(s, handle)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
