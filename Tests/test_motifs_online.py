# Copyright 2015 by Peter Cock.  All rights reserved.
# Copyright 2015 by Antony Lee. All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Testing online code for Bio.motifs (weblogo etc)."""
import os
import unittest

# We want to test these:
from Bio import motifs

# In order to check any sequences returned
from Bio.Seq import Seq

import requires_internet

requires_internet.check()


class TestMotifWeblogo(unittest.TestCase):
    """Tests Bio.motifs online code."""

    def check(self, seqs_as_strs, alpha):
        # Using Seq objects:
        m = motifs.create([Seq(s) for s in seqs_as_strs], alpha)
        m.weblogo(os.devnull)
        # Using strings:
        m = motifs.create(seqs_as_strs, alpha)
        m.weblogo(os.devnull)

    def test_dna(self):
        """Test Bio.motifs.weblogo with a DNA sequence."""
        self.check(
            ["TACAA", "TACGC", "TACAC", "TACCC", "AACCC", "AATGC", "AATGC"], "GATCBDSW"
        )

    def test_rna(self):
        """Test Bio.motifs.weblogo with an RNA sequence."""
        self.check(
            ["UACAA", "UACGC", "UACAC", "UACCC", "AACCC", "AAUGC", "AAUGC"], "GAUC"
        )

    def test_protein(self):
        """Test Bio.motifs.weblogo with a protein sequence."""
        self.check(
            ["ACDEG", "AYCRN", "HYLID", "AYHEL", "ACDEH", "AYYRN", "HYIID"],
            "ACDEFGHIKLMNPQRSTVWYBXZJUO",
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
