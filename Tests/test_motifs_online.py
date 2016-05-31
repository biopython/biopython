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

    def test_dna(self):
        """Test Bio.Motif.weblogo with a DNA sequence."""
        a = extended_dna
        m = motifs.create([Seq("TACAA", a), Seq("TACGC", a),
                           Seq("TACAC", a), Seq("TACCC", a),
                           Seq("AACCC", a), Seq("AATGC", a),
                           Seq("AATGC", a)], alphabet=a)
        m.weblogo(os.devnull)

    def test_rna(self):
        """Test Bio.Motif.weblogo with an RNA sequence."""
        a = unambiguous_rna
        m = motifs.create([Seq("UACAA", a), Seq("UACGC", a),
                           Seq("UACAC", a), Seq("UACCC", a),
                           Seq("AACCC", a), Seq("AAUGC", a),
                           Seq("AAUGC", a)], alphabet=a)
        m.weblogo(os.devnull)

    def test_protein(self):
        """Test Bio.Motif.weblogo with a protein sequence."""
        a = extended_protein
        m = motifs.create([Seq("ACDEG", a), Seq("AYCRN", a),
                           Seq("HYLID", a), Seq("AYHEL", a),
                           Seq("ACDEH", a), Seq("AYYRN", a),
                           Seq("HYIID", a)], alphabet=a)
        m.weblogo(os.devnull)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
