# Copyright 2016 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.Align.AlignInfo related tests."""
import unittest

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align.AlignInfo import SummaryInfo
from Bio.Data import IUPACData
import math


class AlignInfoTests(unittest.TestCase):
    """Test basic usage."""

    def assertAlmostEqualList(self, list1, list2, **kwargs):
        self.assertEqual(len(list1), len(list2))
        for (v1, v2) in zip(list1, list2):
            self.assertAlmostEqual(v1, v2, **kwargs)

    def test_nucleotides(self):
        filename = "GFF/multi.fna"
        fmt = "fasta"
        alignment = AlignIO.read(filename, fmt)
        summary = SummaryInfo(alignment)

        c = summary.dumb_consensus(ambiguous="N")
        self.assertEqual(c, "NNNNNNNN")

        c = summary.gap_consensus(ambiguous="N")
        self.assertEqual(c, "NNNNNNNN")

        expected = {"A": 0.25, "G": 0.25, "T": 0.25, "C": 0.25}

        m = summary.pos_specific_score_matrix(chars_to_ignore=["-"], axis_seq=c)
        self.assertEqual(
            str(m),
            """    A   C   G   T
N  2.0 0.0 1.0 0.0
N  1.0 1.0 1.0 0.0
N  1.0 0.0 2.0 0.0
N  0.0 1.0 1.0 1.0
N  1.0 2.0 0.0 0.0
N  0.0 2.0 1.0 0.0
N  1.0 2.0 0.0 0.0
N  0.0 2.0 1.0 0.0
""",
        )

        # provide the frequencies and chars to ignore explicitly.
        ic = summary.information_content(e_freq_table=expected, chars_to_ignore=["-"])
        self.assertAlmostEqual(ic, 7.32029999423075, places=6)

    def test_proteins(self):
        a = MultipleSeqAlignment(
            [
                SeqRecord(Seq("MHQAIFIYQIGYP*LKSGYIQSIRSPEYDNW-"), id="ID001"),
                SeqRecord(Seq("MH--IFIYQIGYAYLKSGYIQSIRSPEY-NW*"), id="ID002"),
                SeqRecord(Seq("MHQAIFIYQIGYPYLKSGYIQSIRSPEYDNW*"), id="ID003"),
            ]
        )
        self.assertEqual(32, a.get_alignment_length())

        s = SummaryInfo(a)

        c = s.dumb_consensus(ambiguous="X")
        self.assertEqual(c, "MHQAIFIYQIGYXXLKSGYIQSIRSPEYDNW*")

        c = s.gap_consensus(ambiguous="X")
        self.assertEqual(c, "MHXXIFIYQIGYXXLKSGYIQSIRSPEYXNWX")

        m = s.pos_specific_score_matrix(chars_to_ignore=["-", "*"], axis_seq=c)
        self.assertEqual(
            str(m),
            """    A   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   W   Y
M  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 3.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
H  0.0 0.0 0.0 0.0 0.0 3.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
X  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 2.0 0.0 0.0 0.0 0.0
X  2.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
I  0.0 0.0 0.0 0.0 0.0 0.0 3.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
F  0.0 0.0 0.0 3.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
I  0.0 0.0 0.0 0.0 0.0 0.0 3.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
Y  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 3.0
Q  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 3.0 0.0 0.0 0.0 0.0
I  0.0 0.0 0.0 0.0 0.0 0.0 3.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
G  0.0 0.0 0.0 0.0 3.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
Y  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 3.0
X  1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 2.0 0.0 0.0 0.0 0.0 0.0
X  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 2.0
L  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 3.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
K  0.0 0.0 0.0 0.0 0.0 0.0 0.0 3.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
S  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 3.0 0.0 0.0
G  0.0 0.0 0.0 0.0 3.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
Y  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 3.0
I  0.0 0.0 0.0 0.0 0.0 0.0 3.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
Q  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 3.0 0.0 0.0 0.0 0.0
S  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 3.0 0.0 0.0
I  0.0 0.0 0.0 0.0 0.0 0.0 3.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
R  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 3.0 0.0 0.0 0.0
S  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 3.0 0.0 0.0
P  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 3.0 0.0 0.0 0.0 0.0 0.0
E  0.0 0.0 3.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
Y  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 3.0
X  0.0 2.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
N  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 3.0 0.0 0.0 0.0 0.0 0.0 0.0
W  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 3.0 0.0
X  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
""",
        )

        letters = IUPACData.protein_letters
        base_freq = 1.0 / len(letters)
        e_freq_table = {letter: base_freq for letter in letters}
        ic = s.information_content(
            e_freq_table=e_freq_table, chars_to_ignore=["-", "*"]
        )
        self.assertAlmostEqual(ic, 133.061475107, places=6)

    def test_pseudo_count(self):
        # use example from
        # http://biologie.univ-mrs.fr/upload/p202/01.4.PSSM_theory.pdf
        dna_align = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AACCACGTTTAA"), id="ID001"),
                SeqRecord(Seq("CACCACGTGGGT"), id="ID002"),
                SeqRecord(Seq("CACCACGTTCGC"), id="ID003"),
                SeqRecord(Seq("GCGCACGTGGGG"), id="ID004"),
                SeqRecord(Seq("TCGCACGTTGTG"), id="ID005"),
                SeqRecord(Seq("TGGCACGTGTTT"), id="ID006"),
                SeqRecord(Seq("TGACACGTGGGA"), id="ID007"),
                SeqRecord(Seq("TTACACGTGCGC"), id="ID008"),
            ]
        )

        summary = SummaryInfo(dna_align)
        expected = {"A": 0.325, "G": 0.175, "T": 0.325, "C": 0.175}
        ic = summary.information_content(
            e_freq_table=expected, log_base=math.exp(1), pseudo_count=1
        )
        self.assertAlmostEqualList(
            summary.ic_vector,
            [
                0.110,
                0.090,
                0.360,
                1.290,
                0.800,
                1.290,
                1.290,
                0.80,
                0.610,
                0.390,
                0.470,
                0.040,
            ],
            places=2,
        )
        self.assertAlmostEqual(ic, 7.546, places=3)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
