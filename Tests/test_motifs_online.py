# Copyright 2015 by Peter Cock.  All rights reserved.
# Copyright 2015 by Antony Lee. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Testing online code for Bio.motifs (weblogo etc)."""
import os
import unittest

# We want to test these:
from Bio import motifs
from Bio.Alphabet.IUPAC import extended_dna, unambiguous_rna
from Bio.Alphabet.IUPAC import extended_protein

# In order to check any sequences returned
from Bio.Seq import Seq

import requires_internet
requires_internet.check()


class TestotifWeblogo(unittest.TestCase):
    """Tests Bio.motifs online code."""

    def check(self, seqs_as_strs, alpha):
        # Using Seq objects and passing exactly the same alphabet:
        m = motifs.create([Seq(s, alpha) for s in seqs_as_strs], alpha)
        m.weblogo(os.devnull)
        # Using Seq objects but not passing alphabet:
        m = motifs.create([Seq(s, alpha) for s in seqs_as_strs])
        m.weblogo(os.devnull)
        # Using strings and passing alphabet:
        m = motifs.create(seqs_as_strs, alpha)
        m.weblogo(os.devnull)

    def test_dna(self):
        """Test Bio.Motif.weblogo with a DNA sequence."""
        self.check(["TACAA", "TACGC", "TACAC", "TACCC",
                    "AACCC", "AATGC", "AATGC"], extended_dna)

    def test_rna(self):
        """Test Bio.Motif.weblogo with an RNA sequence."""
        self.check(["UACAA", "UACGC", "UACAC", "UACCC",
                    "AACCC", "AAUGC", "AAUGC"], unambiguous_rna)

    def test_protein(self):
        """Test Bio.Motif.weblogo with a protein sequence."""
        self.check(["ACDEG", "AYCRN", "HYLID", "AYHEL",
                    "ACDEH", "AYYRN", "HYIID"], extended_protein)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
