# Copyright 2016 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.Align.AlignInfo related tests."""
import unittest

from Bio import BiopythonDeprecationWarning
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align.AlignInfo import SummaryInfo
from Bio.Data import IUPACData
from Bio.motifs import Motif
import math


class AlignInfoTests(unittest.TestCase):
    """Test basic usage."""

    def assertAlmostEqualList(self, list1, list2, **kwargs):
        self.assertEqual(len(list1), len(list2))
        for v1, v2 in zip(list1, list2):
            self.assertAlmostEqual(v1, v2, **kwargs)

    def test_nucleotides(self):
        filename = "GFF/multi.fna"
        fmt = "fasta"
        msa = AlignIO.read(filename, fmt)
        summary = SummaryInfo(msa)
        alignment = msa.alignment
        motif = Motif("ACGT", alignment)

        with self.assertWarns(BiopythonDeprecationWarning):
            c = summary.dumb_consensus(threshold=0.1, ambiguous="N")
        # dumb_consensus uses ambiguous if multiple letters have the same score
        self.assertEqual(c, "ANGNCCCC")
        c = motif.counts.calculate_consensus(identity=0.1)
        # Instead, EMBOSS uses the first letter it encounters
        self.assertEqual(c, "AaGcCCCC")
        with self.assertWarns(BiopythonDeprecationWarning):
            c = summary.dumb_consensus(ambiguous="N")
        self.assertEqual(c, "NNNNNNNN")
        c = motif.counts.calculate_consensus(identity=0.7)
        self.assertEqual(c, "NNNNNNNN")

        with self.assertWarns(BiopythonDeprecationWarning):
            c = summary.gap_consensus(ambiguous="N")
        self.assertEqual(c, "NNNNNNNN")

        expected = {"A": 0.25, "G": 0.25, "T": 0.25, "C": 0.25}

        with self.assertWarns(BiopythonDeprecationWarning):
            m = summary.pos_specific_score_matrix(chars_to_ignore=["-"], axis_seq=c)

        counts = motif.counts

        for i in range(alignment.length):
            for letter in "ACGT":
                self.assertAlmostEqual(counts[letter][i], m[i][letter])

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
        with self.assertWarns(BiopythonDeprecationWarning):
            ic = summary.information_content(
                e_freq_table=expected, chars_to_ignore=["-"]
            )
        self.assertAlmostEqual(ic, 7.32029999423075)
        ic = sum(motif.relative_entropy)
        self.assertAlmostEqual(ic, 7.32029999423075)

    def test_proteins(self):
        letters = IUPACData.protein_letters
        a = MultipleSeqAlignment(
            [
                SeqRecord(Seq("MHQAIFIYQIGYP*LKSGYIQSIRSPEYDNW-"), id="ID001"),
                SeqRecord(Seq("MH--IFIYQIGYAYLKSGYIQSIRSPEY-NW*"), id="ID002"),
                SeqRecord(Seq("MHQAIFIYQIGYPYLKSGYIQSIRSPEYDNW*"), id="ID003"),
            ]
        )
        self.assertEqual(32, a.get_alignment_length())

        s = SummaryInfo(a)

        alignment = a.alignment
        motif = Motif(letters + "*", alignment)
        counts = motif.counts

        with self.assertWarns(BiopythonDeprecationWarning):
            dumb_consensus = s.dumb_consensus()
        self.assertEqual(dumb_consensus, "MHQAIFIYQIGYXXLKSGYIQSIRSPEYDNW*")
        consensus = counts.calculate_consensus(identity=0.7)
        self.assertEqual(consensus, dumb_consensus)

        with self.assertWarns(BiopythonDeprecationWarning):
            c = s.gap_consensus(ambiguous="X")
        self.assertEqual(c, "MHXXIFIYQIGYXXLKSGYIQSIRSPEYXNWX")

        with self.assertWarns(BiopythonDeprecationWarning):
            m = s.pos_specific_score_matrix(chars_to_ignore=["-", "*"], axis_seq=c)
        j = 0
        all_letters = s._get_all_letters()
        for i in range(alignment.length):
            for letter in letters:
                count = counts[letter][i]
                if letter in all_letters:
                    self.assertAlmostEqual(count, m[j][letter])
                else:
                    self.assertAlmostEqual(count, 0.0)
            j += 1
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

        base_freq = 1.0 / len(letters)
        e_freq_table = {letter: base_freq for letter in letters}
        with self.assertWarns(BiopythonDeprecationWarning):
            ic = s.information_content(
                e_freq_table=e_freq_table, chars_to_ignore=["-", "*"]
            )
        self.assertAlmostEqual(ic, 133.061475107)
        motif = Motif(letters, alignment)
        ic = sum(motif.relative_entropy)
        self.assertAlmostEqual(ic, 133.061475107)

    def test_pseudo_count(self):
        # use example from
        # http://biologie.univ-mrs.fr/upload/p202/01.4.PSSM_theory.pdf
        msa = MultipleSeqAlignment(
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

        summary = SummaryInfo(msa)
        expected = {"A": 0.325, "G": 0.175, "T": 0.325, "C": 0.175}
        with self.assertWarns(BiopythonDeprecationWarning):
            ic = summary.information_content(
                e_freq_table=expected, log_base=math.exp(1), pseudo_count=1
            )
        self.assertAlmostEqual(ic, 7.546369561463767)
        ic_vector = [
            0.11112361,
            0.08677812,
            0.35598044,
            1.29445419,
            0.80272907,
            1.29445419,
            1.29445419,
            0.80272907,
            0.60929642,
            0.39157892,
            0.46539767,
            0.03739368,
        ]
        self.assertAlmostEqualList(summary.ic_vector, ic_vector)
        # One more time, now using a new-style Alignment object:
        alignment = msa.alignment
        motif = Motif("ACGT", alignment)
        motif.background = expected
        motif.pseudocounts = expected
        self.assertAlmostEqualList(motif.relative_entropy * math.log(2), ic_vector)
        ic = sum(ic_vector)
        self.assertAlmostEqual(ic, 7.546369561463767)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
