# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.File module."""

import os.path
import shutil
import tempfile
import unittest
from io import StringIO

from Bio import bgzf
from Bio import File
from Bio import StreamModeError


class RandomAccess(unittest.TestCase):
    """Random access tests."""

    def test_plain(self):
        """Test plain text file."""
        with File._open_for_random_access("Quality/example.fastq") as handle:
            self.assertIn("r", handle.mode)
            self.assertIn("b", handle.mode)

    def test_bgzf(self):
        """Test BGZF compressed file."""
        with File._open_for_random_access("Quality/example.fastq.bgz") as handle:
            self.assertIsInstance(handle, bgzf.BgzfReader)

    def test_gzip(self):
        """Test gzip compressed file."""
        self.assertRaises(
            ValueError, File._open_for_random_access, "Quality/example.fastq.gz"
        )


class AsHandleTestCase(unittest.TestCase):
    """Tests for as_handle function."""

    def setUp(self):
        """Initialise temporary directory."""
        # Create a directory to work in
        self.temp_dir = tempfile.mkdtemp(prefix="biopython-test")

    def tearDown(self):
        """Remove temporary directory."""
        shutil.rmtree(self.temp_dir)

    def _path(self, *args):
        return os.path.join(self.temp_dir, *args)

    def test_handle(self):
        """Test as_handle with a file-like object argument."""
        p = self._path("test_file.fasta")
        with open(p, "wb") as fp:
            with File.as_handle(fp) as handle:
                self.assertEqual(
                    fp,
                    handle,
                    "as_handle should "
                    "return argument when given a "
                    "file-like object",
                )
                self.assertFalse(handle.closed)

            self.assertFalse(
                handle.closed,
                "Exiting as_handle given a file-like object "
                "should not close the file",
            )

    def test_string_path(self):
        """Test as_handle with a string path argument."""
        p = self._path("test_file.fasta")
        mode = "wb"
        with File.as_handle(p, mode=mode) as handle:
            self.assertEqual(p, handle.name)
            self.assertEqual(mode, handle.mode)
            self.assertFalse(handle.closed)
        self.assertTrue(handle.closed)

    def test_path_object(self):
        """Test as_handle with a pathlib.Path object."""
        from pathlib import Path

        p = Path(self._path("test_file.fasta"))
        mode = "wb"
        with File.as_handle(p, mode=mode) as handle:
            self.assertEqual(str(p.absolute()), handle.name)
            self.assertEqual(mode, handle.mode)
            self.assertFalse(handle.closed)
        self.assertTrue(handle.closed)

    def test_custom_path_like_object(self):
        """Test as_handle with a custom path-like object."""

        class CustomPathLike:
            def __init__(self, path):
                self.path = path

            def __fspath__(self):
                return self.path

        p = CustomPathLike(self._path("test_file.fasta"))
        mode = "wb"
        with File.as_handle(p, mode=mode) as handle:
            self.assertEqual(p.path, handle.name)
            self.assertEqual(mode, handle.mode)
            self.assertFalse(handle.closed)
        self.assertTrue(handle.closed)

    def test_stringio(self):
        """Testing passing StringIO handles."""
        s = StringIO()
        with File.as_handle(s) as handle:
            self.assertIs(s, handle)

    def test_as_handle_and_flag(self):
        p = self._path("test_file.fasta")
        fp, shouldclose = File.as_handle_and_flag(p, "w")
        self.assertTrue(shouldclose)
        fp.write(">test\nACGT")
        fp.close()
        fp, shouldclose = File.as_handle_and_flag(p, "r")
        self.assertTrue(shouldclose)
        self.assertEqual(fp.read(), ">test\nACGT")
        fp.close()
        fp, shouldclose = File.as_handle_and_flag(p, "wb")
        self.assertTrue(shouldclose)
        fp.write(b">protein\nACDEFG")
        fp.close()
        fp, shouldclose = File.as_handle_and_flag(p, "rb")
        self.assertTrue(shouldclose)
        self.assertEqual(fp.read(), b">protein\nACDEFG")
        fp.close()
        with open(p, "w") as handle:
            fp, shouldclose = File.as_handle_and_flag(handle, "w")
            self.assertFalse(shouldclose)
            fp, shouldclose = File.as_handle_and_flag(handle, "a")
            self.assertFalse(shouldclose)
            fp.write(">test\nACGT")
            self.assertRaises(StreamModeError, File.as_handle_and_flag, handle, "wb")
            self.assertRaises(StreamModeError, File.as_handle_and_flag, handle, "r")
            self.assertRaises(StreamModeError, File.as_handle_and_flag, handle, "r+")
            self.assertRaises(StreamModeError, File.as_handle_and_flag, handle, "rb")
        with open(p, "rt") as handle:
            fp, shouldclose = File.as_handle_and_flag(handle, "r")
            self.assertFalse(shouldclose)
            fp, shouldclose = File.as_handle_and_flag(handle, "rt")
            self.assertFalse(shouldclose)
            self.assertEqual(fp.read(), ">test\nACGT")
            self.assertRaises(StreamModeError, File.as_handle_and_flag, handle, "rb")
            self.assertRaises(StreamModeError, File.as_handle_and_flag, handle, "r+")
            self.assertRaises(StreamModeError, File.as_handle_and_flag, handle, "w")
            self.assertRaises(StreamModeError, File.as_handle_and_flag, handle, "wb")
        with open(p, "r+") as handle:
            self.assertFalse(File.as_handle_and_flag(handle, "r")[1])
            self.assertFalse(File.as_handle_and_flag(handle, "rt")[1])
            self.assertFalse(File.as_handle_and_flag(handle, "r+")[1])
            self.assertFalse(File.as_handle_and_flag(handle, "w")[1])
            self.assertFalse(File.as_handle_and_flag(handle, "a")[1])
            self.assertRaises(StreamModeError, File.as_handle_and_flag, handle, "wb")
            self.assertRaises(StreamModeError, File.as_handle_and_flag, handle, "rb")
            self.assertRaises(StreamModeError, File.as_handle_and_flag, handle, "r+b")
        with open(p, "wb") as handle:
            fp, shouldclose = File.as_handle_and_flag(handle, "wb")
            self.assertFalse(shouldclose)
            fp, shouldclose = File.as_handle_and_flag(handle, "ab")
            self.assertFalse(shouldclose)
            fp.write(b">protein\nACDEFG")
            self.assertRaises(StreamModeError, File.as_handle_and_flag, handle, "w")
            self.assertRaises(StreamModeError, File.as_handle_and_flag, handle, "r")
            self.assertRaises(StreamModeError, File.as_handle_and_flag, handle, "r+")
            self.assertRaises(StreamModeError, File.as_handle_and_flag, handle, "rb")
        with open(p, "rb") as handle:
            fp, shouldclose = File.as_handle_and_flag(handle, "rb")
            self.assertFalse(shouldclose)
            self.assertEqual(fp.read(), b">protein\nACDEFG")
            self.assertRaises(StreamModeError, File.as_handle_and_flag, handle, "r")
            self.assertRaises(StreamModeError, File.as_handle_and_flag, handle, "r+b")
            self.assertRaises(StreamModeError, File.as_handle_and_flag, handle, "w")
            self.assertRaises(StreamModeError, File.as_handle_and_flag, handle, "wb")
        with open(p, "r+b") as handle:
            self.assertFalse(File.as_handle_and_flag(handle, "rb")[1])
            self.assertFalse(File.as_handle_and_flag(handle, "r+b")[1])
            self.assertFalse(File.as_handle_and_flag(handle, "wb")[1])
            self.assertFalse(File.as_handle_and_flag(handle, "ab")[1])
            self.assertRaises(StreamModeError, File.as_handle_and_flag, handle, "w")
            self.assertRaises(StreamModeError, File.as_handle_and_flag, handle, "r")
            self.assertRaises(StreamModeError, File.as_handle_and_flag, handle, "r+")
            self.assertRaises(StreamModeError, File.as_handle_and_flag, handle, "r+t")


class BaseClassTests(unittest.TestCase):
    """Tests for _IndexedSeqFileProxy base class."""

    def test_instance_exception(self):
        self.assertRaises(TypeError, File._IndexedSeqFileProxy)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
