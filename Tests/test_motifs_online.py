# Copyright 2015 by Peter Cock.  All rights reserved.
# Copyright 2015 by Antony Lee. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Testing online code for Bio.motifs (weblogo etc)."""
import os
import unittest

import requires_internet
requires_internet.check()

# We want to test these:
from Bio import motifs
from Bio.Alphabet import IUPAC

# In order to check any sequences returned
from Bio.Seq import Seq


class TestDNAMotifWeblogo(unittest.TestCase):
    """Tests Bio.motifs online code with DNA sequences."""

    def setUp(self):
        """Create DNA motif for testing."""
        self.m = motifs.create(
            [
                Seq("TACAA"), Seq("TACGC"), Seq("TACAC"), Seq("TACCC"),
                Seq("AACCC"), Seq("AATGC"), Seq("AATGC")
            ],
            alphabet=IUPAC.extended_dna,
        )

    def test_weblogo(self):
        """Test Bio.Motif.weblogo with a DNA sequence."""
        self.m.weblogo(os.devnull)


class TestRNAMotifWeblogo(unittest.TestCase):
    """Tests Bio.motifs online code with RNA sequences."""

    def setUp(self):
        """Create RNA motif for testing."""
        self.m = motifs.create(
            [
                Seq("UACAA"), Seq("UACGC"), Seq("UACAC"), Seq("UACCC"),
                Seq("AACCC"), Seq("AAUGC"), Seq("AAUGC")
            ],
            alphabet=IUPAC.unambiguous_rna,
        )

    def test_weblogo(self):
        """Test Bio.Motif.weblogo with an RNA sequence."""
        self.m.weblogo(os.devnull)


class TestProteinMotifWeblogo(unittest.TestCase):
    """Tests Bio.motifs online code with protein sequences."""

    def setUp(self):
        """Create protein motif for testing."""
        self.m = motifs.create(
            [
                Seq("ACDEG"), Seq("AYCRN"), Seq("HYLID"), Seq("AYHEL"),
                Seq("ACDEH"), Seq("AYYRN"), Seq("HYIID")
            ],
            alphabet=IUPAC.extended_protein,
        )

    def test_weblogo(self):
        """Test Bio.Motif.weblogo with a protein sequence."""
        self.m.weblogo(os.devnull)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
