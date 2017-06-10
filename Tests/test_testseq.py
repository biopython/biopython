# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
import unittest

from Bio.Alphabet import IUPAC, NucleotideAlphabet
from Bio.SeqUtils import GC

file_path = os.getcwd()
file_path = file_path.split("/")
for i in reversed(range(len(file_path))):
    if file_path[i] == 'Tests':
        file_path[i] = 'Scripts/testseq.py'
        break
    else:
        del file_path[i]
    if len(file_path) == 0:
        raise ImportError('File path not found.')
file_path = "/".join(file_path)
run_testseq = "python " + file_path + " unittest"
os.system(run_testseq)


class TestTestseq(unittest.TestCase):

    def test_alphabet(self):
        """Testing 'alphabet' argument..."""
        seq = testseq()
        seq = seq._data
        self.assertEqual(seq, "ATGTCCTCTAATAGTATGGTCGTCTACTGA")

        seq = testseq(alphabet=IUPAC.ambiguous_dna)
        seq = seq._data
        self.assertEqual(seq, "ATGNWGRMSWNVDRSYKNNCMTRTVAKTRA")

        seq = testseq(alphabet=IUPAC.extended_dna)
        seq = seq._data
        self.assertEqual(seq, "ATGBWSBWDCTBTABTBAADWADSDCWTGA")

        seq = testseq(alphabet=IUPAC.unambiguous_rna)
        seq = seq._data
        self.assertEqual(seq, "AUGUCCUCUAAUAGUAUGGUCGUCUACUGA")

        seq = testseq(alphabet=IUPAC.ambiguous_rna)
        seq = seq._data
        self.assertEqual(seq, "AUGNWGRMSWNVDRSYKNNCMURUVAKURA")

        seq = testseq(alphabet=IUPAC.protein)
        seq = seq._data
        self.assertEqual(seq, "MQCKTSPLSNWHTFLFEYKVYFLEDMSVE*")

        seq = testseq(alphabet=IUPAC.extended_protein)
        seq = seq._data
        self.assertEqual(seq, "MUQCKTSPOLSNWHTFLFUEYOKVZOYFL*")

        with self.assertRaises(TypeError):
            testseq(alphabet=NucleotideAlphabet)

    def test_size(self):
        """Testing 'size' and 'truncate' arguments..."""
        seq = testseq()
        self.assertEqual(len(seq), 30)

        seq = testseq(100)
        self.assertEqual(len(seq), 99)

        seq = testseq(100, alphabet=IUPAC.protein)
        self.assertEqual(len(seq), 100)

        seq = testseq(100, truncate=False)
        self.assertEqual(len(seq), 100)

    def test_gc_target(self):
        """Testing 'gc_target' argument..."""
        seq = testseq(1000, gc_target=90)
        self.assertEqual(GC(seq), 89.48948948948949)

    def test_codon_tables(self):
        """Testing codon sets..."""
        seq = testseq(table=5)
        seq1 = seq.translate(table=5)._data
        seq2 = seq.translate(table=6)._data

        self.assertTrue(seq1[0] == "M")
        self.assertFalse("*" in seq1[1:-1])
        self.assertTrue(seq1[-1] == "*")
        self.assertFalse(seq2[-1] == "*")

        seq = testseq(size=1000, from_start=False, to_stop=False, persistent=False)
        seq = seq.translate()._data

        self.assertFalse(seq[0] == "M")
        self.assertFalse(seq[-1] == "*")
        self.assertTrue("*" in seq[1:-1])

    def test_messenger(self):
        """Testing 'messenger' argument..."""
        seq = testseq(alphabet=IUPAC.unambiguous_rna, messenger=True)._data
        self.assertTrue("A" * 20 == seq[-20:])

    def test_seeding(self):
        """Testing 'rand_seed' argument..."""
        seq1 = testseq(rand_seed=None)
        seq2 = testseq(rand_seed=None)
        self.assertFalse(seq1 == seq2)

        seq1 = testseq(rand_seed=50)
        seq2 = testseq(rand_seed=50)
        self.assertTrue(seq1 == seq2)


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
