# Copyright 2013 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.Application related tests for command line application wrappers.

This is intended to check generic things like argument parsing, and
stdin/stdout/stderr handling.
"""

import os
import unittest

from Bio.Application import AbstractCommandline, _Argument


class EchoApp(AbstractCommandline):
    """Minimal command line wrapper for echo command."""

    def __init__(self, cmd="echo", **kwargs):
        """Initialize wrapper for echo command."""
        self.parameters = [_Argument(["text"], "Text to echo")]
        AbstractCommandline.__init__(self, cmd, **kwargs)


class TestApp(unittest.TestCase):
    def test_echo(self):
        cline = EchoApp(text="Hello World")
        stdout, stderr = cline()
        self.assertEqual(stderr, "")
        self.assertEqual(stdout, "Hello World\n")

    def test_echo_capture_both(self):
        cline = EchoApp(text="Hello World")
        stdout, stderr = cline(stdout=True, stderr=True)
        self.assertEqual(stderr, "")
        self.assertEqual(stdout, "Hello World\n")

    def test_echo_capture_stdout(self):
        cline = EchoApp(text="Hello World")
        stdout, stderr = cline(stdout=True, stderr=False)
        self.assertEqual(stderr, None)
        self.assertEqual(stdout, "Hello World\n")

    def test_echo_capture_stderr(self):
        cline = EchoApp(text="Hello World")
        stdout, stderr = cline(stdout=False, stderr=True)
        self.assertEqual(stderr, "")
        self.assertEqual(stdout, None)

    def test_echo_capture_neither(self):
        cline = EchoApp(text="Hello World")
        stdout, stderr = cline(stdout=False, stderr=False)
        self.assertEqual(stderr, None)
        self.assertEqual(stdout, None)

    def test_echo_file_stdout(self):
        cline = EchoApp(text="Hello World")
        tmp = "echo_stdout.tmp"
        if os.path.isfile(tmp):
            os.remove(tmp)
        stdout, stderr = cline(stdout=tmp)
        self.assertEqual(stderr, "")
        self.assertEqual(stdout, None)
        self.assertTrue(os.path.isfile(tmp))
        with open(tmp) as h:
            contents = h.read()
        self.assertEqual(contents, "Hello World\n")
        os.remove(tmp)

    def test_echo_file_stderr(self):
        cline = EchoApp(text="Hello World")
        tmp = "echo_stderr.tmp"
        if os.path.isfile(tmp):
            os.remove(tmp)
        stdout, stderr = cline(stderr=tmp)
        self.assertEqual(stderr, None)
        self.assertEqual(stdout, "Hello World\n")
        self.assertTrue(os.path.isfile(tmp))
        with open(tmp) as h:
            contents = h.read()
            self.assertEqual(contents, "")
        os.remove(tmp)

    def test_echo_file_same(self):
        cline = EchoApp(text="Hello World")
        tmp = "echo_stdout_stderr.tmp"
        if os.path.isfile(tmp):
            os.remove(tmp)
        stdout, stderr = cline(stdout=tmp, stderr=tmp)
        self.assertEqual(stderr, None)
        self.assertEqual(stdout, None)
        self.assertTrue(os.path.isfile(tmp))
        with open(tmp) as h:
            contents = h.read()
        self.assertEqual(contents, "Hello World\n")  # stdout + stderr
        os.remove(tmp)

    def test_echo_file_both(self):
        cline = EchoApp(text="Hello World")
        tmp = "echo_stdout.tmp"
        if os.path.isfile(tmp):
            os.remove(tmp)
        tmp2 = "echo_stderr.tmp"
        if os.path.isfile(tmp2):
            os.remove(tmp2)
        stdout, stderr = cline(stdout=tmp, stderr=tmp2)
        self.assertEqual(stderr, None)
        self.assertEqual(stdout, None)
        self.assertTrue(os.path.isfile(tmp), tmp)
        with open(tmp) as h:
            contents = h.read()
        self.assertEqual(contents, "Hello World\n")  # stdout
        os.remove(tmp)
        self.assertTrue(os.path.isfile(tmp2), tmp2)
        with open(tmp2) as h:
            contents = h.read()
        self.assertEqual(contents, "")  # stderr
        os.remove(tmp2)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
