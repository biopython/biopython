# Copyright 2006-2014 by Peter Cock.  All rights reserved.
# Revisions copyright 2011 Brandon Invergo. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for Bio.AlignIO.MauveIO"""

import unittest

from Bio._py3k import StringIO

from Bio.AlignIO.MauveIO import MauveIterator, MauveWriter


class TestMauveIO(unittest.TestCase):

    def test_one(self):
        handle = open('Mauve/simple.xmfa')
        ids = []
        for alignment in MauveIterator(handle):
            for record in alignment:
                ids.append(record.id)
        self.assertEqual(ids, ['1', '2', '3', '1', '2', '3'])

        expected = """TTCGGTACCCTCCATGACCCACGAAATGAGGGCCCAGGGTATGCTT"""
        self.assertEqual(str(record.seq).replace("-", ""), expected)

    def test_write_read(self):
        handle = open('Mauve/simple.xmfa')
        aln_list = list(MauveIterator(handle))
        handle.close()

        handle = StringIO()
        MauveWriter(handle).write_file(aln_list)
        handle.seek(0)
        aln_list_out = list(MauveIterator(handle))

        for a1, a2 in zip(aln_list, aln_list_out):
            self.assertEqual(len(a1), len(a2))
            for r1, r2 in zip(a1, a2):
                self.assertEqual(r1.id, r2.id)
                self.assertEqual(str(r1.seq), str(r2.seq))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
