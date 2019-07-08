# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO legacy module."""


import string
import unittest

from Bio._py3k import StringIO
from Bio import File


import warnings
from Bio import BiopythonWarning
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonWarning)
    from Bio.SearchIO._legacy import ParserSupport


class TestParserSupport(unittest.TestCase):

    def test_TaggingConsumer(self):

        h = StringIO()
        tc = ParserSupport.TaggingConsumer(handle=h, colwidth=5)
        tc.start_section()
        self.assertEqual(h.getvalue(), '***** start_section\n')
        h.seek(0)
        h.truncate(0)
        tc.test1('myline')
        self.assertEqual(h.getvalue(), 'test1: myline\n')
        h.seek(0)
        h.truncate(0)
        tc.end_section()
        self.assertEqual(h.getvalue(), '***** end_section\n')

    def test_is_blank_line(self):
        is_blank_line = ParserSupport.is_blank_line
        self.assertTrue(is_blank_line('\n'))
        self.assertTrue(is_blank_line('\r\n'))
        self.assertTrue(is_blank_line('\r'))
        self.assertTrue(is_blank_line(''))
        self.assertTrue(is_blank_line('', allow_spaces=1))
        self.assertTrue(is_blank_line('', allow_spaces=0))
        self.assertTrue(is_blank_line(string.whitespace, allow_spaces=1))
        self.assertFalse(is_blank_line('hello'))
        self.assertFalse(is_blank_line('hello', allow_spaces=1))
        self.assertFalse(is_blank_line('hello', allow_spaces=0))
        self.assertFalse(is_blank_line(string.whitespace, allow_spaces=0))

    def test_safe_readline(self):
        data = """\
This
file"""
        h = File.UndoHandle(StringIO(data))
        safe_readline = ParserSupport.safe_readline
        self.assertEqual(safe_readline(h), "This\n")
        self.assertEqual(safe_readline(h), "file")
        self.assertRaises(ValueError, safe_readline, h)

    def test_safe_peekline(self):
        safe_peekline = ParserSupport.safe_peekline
        data = """\
This
file"""
        h = File.UndoHandle(StringIO(data))
        self.assertEqual(safe_peekline(h), "This\n")
        h.readline()
        self.assertEqual(safe_peekline(h), "file")
        h.readline()
        self.assertRaises(ValueError, safe_peekline, h)
        h.saveline('hello')
        self.assertEqual(safe_peekline(h), 'hello')

    def test_read_and_call(self):
        data = """\
>gi|132871|sp|P19947|RL30_BACSU 50S RIBOSOMAL PROTEIN L30 (BL27)
MAKLEITLKRSVIGRPEDQRVTVRTLGLKKTNQTVVHEDNAAIRGMINKVSHLVSVKEQ
>gi|132679|sp|P19946|RL15_BACSU 50S RIBOSOMAL PROTEIN L15
MKLHELKPSEGSRKTRNRVGRGIGSGNGKTAGKGHKGQNARSGGGVRPGFEGGQMPLFQRLPKRGFTNIN
RKEYAVVNLDKLNGFAEGTEVTPELLLETGVISKLNAGVKILGNGKLEKKLTVKANKFSASAKEAVEAAG
GTAEVI


"""
        h = File.UndoHandle(StringIO(data))
        lines = []
        rac = ParserSupport.read_and_call
        rac(h, lines.append)
        self.assertEqual(lines[-1][:10], '>gi|132871')
        rac(h, lines.append, start='MAKLE', end='KEQ', contains='SVIG')
        self.assertRaises(ValueError, rac, h, lines.append, blank=1)
        self.assertRaises(ValueError, rac, h, lines.append, start='foobar')
        self.assertRaises(ValueError, rac, h, lines.append, end='foobar')
        self.assertRaises(ValueError, rac, h, lines.append, contains='foobar')
        self.assertRaises(ValueError, rac, h, lines.append, blank=0)

    def test_attempt_read_and_call(self):
        data = """\
>gi|132871|sp|P19947|RL30_BACSU 50S RIBOSOMAL PROTEIN L30 (BL27)
MAKLEITLKRSVIGRPEDQRVTVRTLGLKKTNQTVVHEDNAAIRGMINKVSHLVSVKEQ
>gi|132679|sp|P19946|RL15_BACSU 50S RIBOSOMAL PROTEIN L15
MKLHELKPSEGSRKTRNRVGRGIGSGNGKTAGKGHKGQNARSGGGVRPGFEGGQMPLFQRLPKRGFTNIN
RKEYAVVNLDKLNGFAEGTEVTPELLLETGVISKLNAGVKILGNGKLEKKLTVKANKFSASAKEAVEAAG
GTAEVI"""
        h = File.UndoHandle(StringIO(data))
        lines = []
        arac = ParserSupport.attempt_read_and_call
        self.assertTrue(arac(h, lines.append, contains="RIBOSOMAL PROTEIN"))
        self.assertFalse(arac(h, lines.append, start="foobar"))
        self.assertFalse(arac(h, lines.append, blank=1))
        self.assertTrue(arac(h, lines.append, end="LVSVKEQ"))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
