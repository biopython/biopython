# Copyright 2021 by Michiel de Hoon.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.Align.nexus module."""
import unittest
from io import StringIO
from tempfile import NamedTemporaryFile

from Bio import Align

try:
    import numpy as np
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align.nexus."
    ) from None


class TestNexusReading(unittest.TestCase):
    def check_reading_writing(self, path):
        alignments = Align.parse(path, "nexus")
        stream = StringIO()
        n = Align.write(alignments, stream, "nexus")
        self.assertEqual(n, 1)
        alignments = Align.parse(path, "nexus")
        alignments = list(alignments)
        alignment = alignments[0]
        stream.seek(0)
        saved_alignments = Align.parse(stream, "nexus")
        saved_alignments = list(saved_alignments)
        self.assertEqual(len(alignments), len(saved_alignments))
        saved_alignment = saved_alignments[0]
        for i, (sequence, saved_sequence) in enumerate(
            zip(alignment.sequences, saved_alignment.sequences)
        ):
            self.assertEqual(sequence.id, saved_sequence.id)
            self.assertEqual(sequence.seq, saved_sequence.seq)
            self.assertEqual(sequence.annotations, saved_sequence.annotations)
            self.assertEqual(alignment[i], saved_alignment[i])
            self.assertTrue(
                np.array_equal(alignment.coordinates, saved_alignment.coordinates)
            )

    def test_nexus1(self):
        path = "Nexus/test_Nexus_input.nex"
        alignments = Align.parse(path, "nexus")
        self.check_nexus1(alignments)
        alignments = iter(alignments)
        self.check_nexus1(alignments)
        with Align.parse(path, "nexus") as alignments:
            self.check_nexus1(alignments)
        with self.assertRaises(AttributeError):
            alignments._stream
        with Align.parse(path, "nexus") as alignments:
            pass
        with self.assertRaises(AttributeError):
            alignments._stream
        self.check_reading_writing(path)
        with open(path) as stream:
            data = stream.read()
        stream = NamedTemporaryFile("w+t")
        stream.write(data)
        stream.seek(0)
        alignments = Align.parse(stream, "nexus")
        self.check_nexus1(alignments)

    def check_nexus1(self, alignments):
        alignment = next(alignments)
        self.assertEqual(len(alignment), 9)
        self.assertEqual(alignment.shape, (9, 46))
        self.assertEqual(alignment.sequences[0].id, "t1")
        self.assertEqual(alignment.sequences[1].id, "t2 the name")
        self.assertEqual(alignment.sequences[2].id, "isn'that [a] strange name?")
        self.assertEqual(
            alignment.sequences[3].id, "one should be punished, for (that)!"
        )
        self.assertEqual(alignment.sequences[4].id, "t5")
        self.assertEqual(alignment.sequences[5].id, "t6")
        self.assertEqual(alignment.sequences[6].id, "t7")
        self.assertEqual(alignment.sequences[7].id, "t8")
        self.assertEqual(alignment.sequences[8].id, "t9")
        self.assertEqual(alignment.sequences[0].annotations, {"molecule_type": "DNA"})
        self.assertEqual(alignment.sequences[1].annotations, {"molecule_type": "DNA"})
        self.assertEqual(alignment.sequences[2].annotations, {"molecule_type": "DNA"})
        self.assertEqual(alignment.sequences[3].annotations, {"molecule_type": "DNA"})
        self.assertEqual(alignment.sequences[4].annotations, {"molecule_type": "DNA"})
        self.assertEqual(alignment.sequences[5].annotations, {"molecule_type": "DNA"})
        self.assertEqual(alignment.sequences[6].annotations, {"molecule_type": "DNA"})
        self.assertEqual(alignment.sequences[7].annotations, {"molecule_type": "DNA"})
        self.assertEqual(alignment.sequences[8].annotations, {"molecule_type": "DNA"})
        self.assertEqual(
            alignment.sequences[0].seq, "ACGTcgtgtgtgctctttacgtgtgtgctcttt"
        )
        self.assertEqual(alignment.sequences[1].seq, "ACGcTcgtgtctttacacgtgtcttt")
        self.assertEqual(alignment.sequences[2].seq, "ACcGcTcgtgtgtgctacacacgtgtgtgct")
        self.assertEqual(alignment.sequences[3].seq, "ACGT")
        self.assertEqual(
            alignment.sequences[4].seq, "AC?GT?acgt???????????acgt????????"
        )
        self.assertEqual(
            alignment.sequences[5].seq, "AcCaGtTc?aaaaaaaaaaacgactac?aaaaaaaaaa"
        )
        self.assertEqual(
            alignment.sequences[6].seq, "A?CGgTgggggggggggggg???gggggggggggggggg"
        )
        self.assertEqual(
            alignment.sequences[7].seq, "AtCtGtTtttttttttttt??ttttttttttttttttttt??"
        )
        self.assertEqual(
            alignment.sequences[8].seq, "cccccccccccccccccccNcccccccccccccccccccccNcc"
        )
        self.assertEqual(
            str(alignment),
            """\
t1                0 A-C-G-Tcgtgtgtgctct-t-t------acgtgtgtgctct-t-t 33
t2 the na         0 A-C-GcTcgtg-----tct-t-t----acacgtg-----tct-t-t 26
isn'that          0 A-CcGcTcgtgtgtgct--------acacacgtgtgtgct------ 31
one shoul         0 A-C-G-T---------------------------------------  4
t5                0 A-C?G-T?-acgt??-???-???--??---?-acgt??-???-??? 33
t6                0 AcCaGtTc?--aaaaaaaa-a-aacgactac?--aaaaaaaa-a-a 38
t7                0 A?C-GgTgggggggggggg-g-g??--?gggggggggggggg-g-g 39
t8                0 AtCtGtTtttttttttttt-?-?ttttttttttttttttttt-?-? 42
t9                0 cccccccccccccccccccNc-ccccccccccccccccccccNc-c 44
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[ 0,  1,  1,  2,  2,  3,  3,  4,  5,  6,  8, 12,
                           13, 14, 16, 16, 17, 17, 18, 18, 18, 18, 19, 20,
                           21, 23, 27, 28, 29, 31, 31, 32, 32, 33],
                          [ 0,  1,  1,  2,  2,  3,  4,  5,  6,  7,  9,  9,
                            9, 10, 12, 12, 13, 13, 14, 14, 14, 16, 17, 18,
                           19, 21, 21, 21, 22, 24, 24, 25, 25, 26],
                          [ 0,  1,  1,  2,  3,  4,  5,  6,  7,  8, 10, 14,
                           15, 16, 16, 16, 16, 16, 16, 16, 18, 20, 21, 22,
                           23, 25, 29, 30, 31, 31, 31, 31, 31, 31],
                          [ 0,  1,  1,  2,  2,  3,  3,  4,  4,  4,  4,  4,
                            4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
                            4,  4,  4,  4,  4,  4,  4,  4,  4,  4],
                          [ 0,  1,  1,  2,  3,  4,  4,  5,  6,  6,  8, 12,
                           12, 13, 15, 15, 16, 17, 18, 18, 20, 20, 20, 21,
                           21, 23, 27, 27, 28, 30, 30, 31, 32, 33],
                          [ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  9, 13,
                           14, 15, 17, 17, 18, 18, 19, 21, 23, 25, 26, 27,
                           28, 28, 32, 33, 34, 36, 36, 37, 37, 38],
                          [ 0,  1,  2,  3,  3,  4,  5,  6,  7,  8, 10, 14,
                           15, 16, 18, 18, 19, 19, 20, 22, 22, 24, 25, 26,
                           27, 29, 33, 34, 35, 37, 37, 38, 38, 39],
                          [ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 11, 15,
                           16, 17, 19, 19, 20, 20, 21, 23, 25, 27, 28, 29,
                           30, 32, 36, 37, 38, 40, 40, 41, 41, 42],
                          [ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 11, 15,
                           16, 17, 19, 20, 21, 21, 22, 24, 26, 28, 29, 30,
                           31, 33, 37, 38, 39, 41, 42, 43, 43, 44]])
                # fmt: on
            )
        )
        self.assertEqual(
            alignment[0],
            "A-C-G-Tcgtgtgtgctct-t-t------acgtgtgtgctct-t-t",
        )
        self.assertEqual(
            alignment[1],
            "A-C-GcTcgtg-----tct-t-t----acacgtg-----tct-t-t",
        )
        self.assertEqual(alignment[2], "A-CcGcTcgtgtgtgct--------acacacgtgtgtgct------")
        self.assertEqual(alignment[3], "A-C-G-T---------------------------------------")
        self.assertEqual(alignment[4], "A-C?G-T?-acgt??-???-???--??---?-acgt??-???-???")
        self.assertEqual(alignment[5], "AcCaGtTc?--aaaaaaaa-a-aacgactac?--aaaaaaaa-a-a")
        self.assertEqual(alignment[6], "A?C-GgTgggggggggggg-g-g??--?gggggggggggggg-g-g")
        self.assertEqual(alignment[7], "AtCtGtTtttttttttttt-?-?ttttttttttttttttttt-?-?")
        self.assertEqual(alignment[8], "cccccccccccccccccccNc-ccccccccccccccccccccNc-c")
        self.assertEqual(
            format(alignment, "nexus"),
            """\
#NEXUS
begin data;
dimensions ntax=9 nchar=46;
format datatype=dna missing=? gap=-;
matrix
t1                                    A-C-G-Tcgtgtgtgctct-t-t------acgtgtgtgctct-t-t
't2 the name'                         A-C-GcTcgtg-----tct-t-t----acacgtg-----tct-t-t
'isn''that [a] strange name?'         A-CcGcTcgtgtgtgct--------acacacgtgtgtgct------
'one should be punished, for (that)!' A-C-G-T---------------------------------------
t5                                    A-C?G-T?-acgt??-???-???--??---?-acgt??-???-???
t6                                    AcCaGtTc?--aaaaaaaa-a-aacgactac?--aaaaaaaa-a-a
t7                                    A?C-GgTgggggggggggg-g-g??--?gggggggggggggg-g-g
t8                                    AtCtGtTtttttttttttt-?-?ttttttttttttttttttt-?-?
t9                                    cccccccccccccccccccNc-ccccccccccccccccccccNc-c
;
end;
""",
        )
        self.assertTrue(
            np.array_equal(
                np.array(alignment, "U"),
                # fmt: off
np.array([['A', '-', 'C', '-', 'G', '-', 'T', 'c', 'g', 't', 'g', 't', 'g',
           't', 'g', 'c', 't', 'c', 't', '-', 't', '-', 't', '-', '-', '-',
           '-', '-', '-', 'a', 'c', 'g', 't', 'g', 't', 'g', 't', 'g', 'c',
           't', 'c', 't', '-', 't', '-', 't'],
          ['A', '-', 'C', '-', 'G', 'c', 'T', 'c', 'g', 't', 'g', '-', '-',
           '-', '-', '-', 't', 'c', 't', '-', 't', '-', 't', '-', '-', '-',
           '-', 'a', 'c', 'a', 'c', 'g', 't', 'g', '-', '-', '-', '-', '-',
           't', 'c', 't', '-', 't', '-', 't'],
          ['A', '-', 'C', 'c', 'G', 'c', 'T', 'c', 'g', 't', 'g', 't', 'g',
           't', 'g', 'c', 't', '-', '-', '-', '-', '-', '-', '-', '-', 'a',
           'c', 'a', 'c', 'a', 'c', 'g', 't', 'g', 't', 'g', 't', 'g', 'c',
           't', '-', '-', '-', '-', '-', '-'],
          ['A', '-', 'C', '-', 'G', '-', 'T', '-', '-', '-', '-', '-', '-',
           '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
           '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
           '-', '-', '-', '-', '-', '-', '-'],
          ['A', '-', 'C', '?', 'G', '-', 'T', '?', '-', 'a', 'c', 'g', 't',
           '?', '?', '-', '?', '?', '?', '-', '?', '?', '?', '-', '-', '?',
           '?', '-', '-', '-', '?', '-', 'a', 'c', 'g', 't', '?', '?', '-',
           '?', '?', '?', '-', '?', '?', '?'],
          ['A', 'c', 'C', 'a', 'G', 't', 'T', 'c', '?', '-', '-', 'a', 'a',
           'a', 'a', 'a', 'a', 'a', 'a', '-', 'a', '-', 'a', 'a', 'c', 'g',
           'a', 'c', 't', 'a', 'c', '?', '-', '-', 'a', 'a', 'a', 'a', 'a',
           'a', 'a', 'a', '-', 'a', '-', 'a'],
          ['A', '?', 'C', '-', 'G', 'g', 'T', 'g', 'g', 'g', 'g', 'g', 'g',
           'g', 'g', 'g', 'g', 'g', 'g', '-', 'g', '-', 'g', '?', '?', '-',
           '-', '?', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g',
           'g', 'g', 'g', '-', 'g', '-', 'g'],
          ['A', 't', 'C', 't', 'G', 't', 'T', 't', 't', 't', 't', 't', 't',
           't', 't', 't', 't', 't', 't', '-', '?', '-', '?', 't', 't', 't',
           't', 't', 't', 't', 't', 't', 't', 't', 't', 't', 't', 't', 't',
           't', 't', 't', '-', '?', '-', '?'],
          ['c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c',
           'c', 'c', 'c', 'c', 'c', 'c', 'N', 'c', '-', 'c', 'c', 'c', 'c',
           'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c',
           'c', 'c', 'c', 'N', 'c', '-', 'c']], dtype='U')
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (862 aligned letters; 256 identities; 606 mismatches; 596 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 862:
        identities = 256,
        mismatches = 606.
    gaps = 596:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 327:
            internal_insertions = 254:
                open_internal_insertions = 128,
                extend_internal_insertions = 126;
            internal_deletions = 73:
                open_internal_deletions = 44,
                extend_internal_deletions = 29;
        right_gaps = 269:
            right_insertions = 186:
                open_right_insertions = 10,
                extend_right_insertions = 176;
            right_deletions = 83:
                open_right_deletions = 5,
                extend_right_deletions = 78.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 186)
        self.assertEqual(counts.right_deletions, 83)
        self.assertEqual(counts.internal_insertions, 254)
        self.assertEqual(counts.internal_deletions, 73)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 269)
        self.assertEqual(counts.internal_gaps, 327)
        self.assertEqual(counts.insertions, 440)
        self.assertEqual(counts.deletions, 156)
        self.assertEqual(counts.gaps, 596)
        self.assertEqual(counts.aligned, 862)
        self.assertEqual(counts.identities, 256)
        self.assertEqual(counts.mismatches, 606)
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_nexus2(self):
        path = "Nexus/codonposset.nex"
        alignments = Align.parse(path, "nexus")
        alignment = next(alignments)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 22))
        self.assertEqual(alignment.sequences[0].id, "Aegotheles")
        self.assertEqual(alignment.sequences[1].id, "Aerodramus")
        self.assertEqual(alignment.sequences[0].annotations, {"molecule_type": "DNA"})
        self.assertEqual(alignment.sequences[1].annotations, {"molecule_type": "DNA"})
        self.assertEqual(alignment.sequences[0].seq, "AAAAAGGCATTGTGGTGGGAAT")
        self.assertEqual(alignment.sequences[1].seq, "?????????TTGTGGTGGGAAT")
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 22], [0, 22]]))
        )
        self.assertEqual(alignment[0], "AAAAAGGCATTGTGGTGGGAAT")
        self.assertEqual(alignment[1], "?????????TTGTGGTGGGAAT")
        self.assertEqual(
            str(alignment),
            """\
Aegothele         0 AAAAAGGCATTGTGGTGGGAAT 22
                  0 .........||||||||||||| 22
Aerodramu         0 ?????????TTGTGGTGGGAAT 22
""",
        )
        self.assertEqual(
            format(alignment, "nexus"),
            """\
#NEXUS
begin data;
dimensions ntax=2 nchar=22;
format datatype=dna missing=? gap=-;
matrix
Aegotheles AAAAAGGCATTGTGGTGGGAAT
Aerodramus ?????????TTGTGGTGGGAAT
;
end;
""",
        )
        self.assertTrue(
            np.array_equal(
                np.array(alignment, "U"),
                # fmt: off
np.array([['A', 'A', 'A', 'A', 'A', 'G', 'G', 'C', 'A', 'T', 'T', 'G', 'T',
           'G', 'G', 'T', 'G', 'G', 'G', 'A', 'A', 'T'],
          ['?', '?', '?', '?', '?', '?', '?', '?', '?', 'T', 'T', 'G', 'T',
           'G', 'G', 'T', 'G', 'G', 'G', 'A', 'A', 'T']], dtype='U')
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (22 aligned letters; 13 identities; 9 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 22:
        identities = 13,
        mismatches = 9.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 22)
        self.assertEqual(counts.identities, 13)
        self.assertEqual(counts.mismatches, 9)
        with self.assertRaises(StopIteration):
            next(alignments)
        self.check_reading_writing(path)


class TestNexusBasic(unittest.TestCase):
    def test_empty(self):
        import io

        stream = io.StringIO()
        with self.assertRaisesRegex(ValueError, "Empty file."):
            Align.parse(stream, "nexus")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
