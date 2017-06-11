# Copyright 2017 by Adil Iqbal.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
import sys
import unittest

from Bio.Alphabet import IUPAC, NucleotideAlphabet
from Bio.SeqUtils import GC

# Find file path to 'testseq' module.
file_path = os.getcwd()
cwd_path = file_path
file_path = os.path.split(file_path)
if file_path[1] == "Tests":
    file_path = os.path.join(file_path[0], "Scripts")
    file_path = os.path.join(file_path, "testseq.py")
else:
    raise ImportError("File path could not be derived from: " + cwd_path)

# Import module manually
if sys.version_info[0] < 3:
    # Python version 2.7
    import imp
    foo = imp.load_source('testseq', file_path)
else:
    if sys.version_info[1] <= 4:
        # Python version 3.3 and 3.4
        from importlib.machinery import SourceFileLoader
        foo = SourceFileLoader('testseq', file_path).load_module()
    else:
        # Python version 3.5+
        from importlib.util import spec_from_file_location, module_from_spec
        spec = spec_from_file_location('testseq', file_path)
        foo = module_from_spec(spec)
        spec.loader.exec_module(foo)


class TestTestseq(unittest.TestCase):
    """Perform unit test for 'testseq' function."""

    def test_alphabet(self):
        """Testing 'alphabet' argument..."""
        alphabets = [IUPAC.unambiguous_dna, IUPAC.ambiguous_dna, IUPAC.extended_dna,
                     IUPAC.unambiguous_rna, IUPAC.ambiguous_rna, IUPAC.protein, IUPAC.extended_protein]

        for i in alphabets:
            seq = foo.testseq(alphabet=i)._data
            for j, k in enumerate(seq):
                self.assertTrue(k in i.letters)

        with self.assertRaises(TypeError):
            foo.testseq(alphabet=NucleotideAlphabet)

    def test_size(self):
        """Testing 'size' and 'truncate' arguments..."""
        seq = foo.testseq()
        self.assertEqual(len(seq), 30)

        seq = foo.testseq(100)
        self.assertEqual(len(seq), 99)

        seq = foo.testseq(100, alphabet=IUPAC.protein)
        self.assertEqual(len(seq), 100)

        seq = foo.testseq(100, truncate=False)
        self.assertEqual(len(seq), 100)

    def test_gc_target(self):
        """Testing 'gc_target' argument..."""
        seq = foo.testseq(1000, gc_target=90)
        error = 90 - GC(seq)
        self.assertTrue(-5 < error < 5)

    def test_codon_tables(self):
        """Testing codon sets..."""
        seq = foo.testseq(table=5)
        seq1 = seq.translate(table=5)._data
        seq2 = seq.translate(table=6)._data

        self.assertTrue(seq1[0] == "M")
        self.assertFalse("*" in seq1[1:-1])
        self.assertTrue(seq1[-1] == "*")
        self.assertFalse(seq2[-1] == "*")

        seq = foo.testseq(size=1000, from_start=False, to_stop=False, persistent=False)
        seq = seq.translate()._data

        self.assertFalse(seq[0] == "M")
        self.assertFalse(seq[-1] == "*")
        self.assertTrue("*" in seq[1:-1])

    def test_messenger(self):
        """Testing 'messenger' argument..."""
        seq = foo.testseq(alphabet=IUPAC.unambiguous_rna, messenger=True)._data
        self.assertTrue("A" * 20 == seq[-20:])

    def test_seeding(self):
        """Testing 'rand_seed' argument..."""
        seq1 = foo.testseq(rand_seed=None)
        seq2 = foo.testseq(rand_seed=None)
        self.assertFalse(seq1 == seq2)

        seq1 = foo.testseq(rand_seed=50)
        seq2 = foo.testseq(rand_seed=50)
        self.assertTrue(seq1 == seq2)


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
