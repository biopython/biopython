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
    def __init__(self, cmd="echo", **kwargs):
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


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
