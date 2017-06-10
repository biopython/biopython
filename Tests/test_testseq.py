# Copyright 2017 by Adil Iqbal.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import unittest

if __name__ == '__main__' and __package__ is None:
    __package__ = 'test_testseq'

from Bio.Alphabet import IUPAC, NucleotideAlphabet
from Bio.SeqUtils import GC

# Find file path to 'testseq' module.
file_path = os.getcwd()
cwd_path = file_path
file_path = file_path.split("\\")
for i in reversed(range(len(file_path))):
    if file_path[i] == "Tests":
        file_path[i] = "Scripts\\testseq.py"
        break
    else:
        del file_path[i]
    if len(file_path) == 0:
        raise ImportError("Could not find file path from", cwd_path)
file_path = "\\".join(file_path)
del cwd_path

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


class TestTestseq(unittest.TestCase):

    def test_alphabet(self):
        """Testing 'alphabet' argument..."""
        seq = foo.testseq()
        seq = seq._data
        self.assertEqual(seq, "ATGTCCTCTAATAGTATGGTCGTCTACTGA")

        seq = foo.testseq(alphabet=IUPAC.ambiguous_dna)
        seq = seq._data
        self.assertEqual(seq, "ATGNWGRMSWNVDRSYKNNCMTRTVAKTRA")

        seq = foo.testseq(alphabet=IUPAC.extended_dna)
        seq = seq._data
        self.assertEqual(seq, "ATGBWSBWDCTBTABTBAADWADSDCWTGA")

        seq = foo.testseq(alphabet=IUPAC.unambiguous_rna)
        seq = seq._data
        self.assertEqual(seq, "AUGUCCUCUAAUAGUAUGGUCGUCUACUGA")

        seq = foo.testseq(alphabet=IUPAC.ambiguous_rna)
        seq = seq._data
        self.assertEqual(seq, "AUGNWGRMSWNVDRSYKNNCMURUVAKURA")

        seq = foo.testseq(alphabet=IUPAC.protein)
        seq = seq._data
        self.assertEqual(seq, "MQCKTSPLSNWHTFLFEYKVYFLEDMSVE*")

        seq = foo.testseq(alphabet=IUPAC.extended_protein)
        seq = seq._data
        self.assertEqual(seq, "MUQCKTSPOLSNWHTFLFUEYOKVZOYFL*")

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
        self.assertEqual(len(sfoo.eq), 100)

    def test_gc_target(self):
        """Testing 'gc_target'foo. argument..."""
        seq = foo.testseq(1000, gc_target=90)
        self.assertEqual(GC(seq), 89.48948948948949)

    def test_codon_tables(self):
        """Testing codon sets.foo..."""
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
