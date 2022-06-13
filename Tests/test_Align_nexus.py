# Copyright 2021 by Michiel de Hoon.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.Align.nexus module."""
import unittest
import warnings
from io import StringIO

from Bio import BiopythonExperimentalWarning

with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonExperimentalWarning)
    from Bio.Align.nexus import AlignmentIterator, AlignmentWriter

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align.nexus."
    ) from None


class TestNexusReading(unittest.TestCase):
    def check_reading_writing(self, path):
        alignments = AlignmentIterator(path)
        stream = StringIO()
        writer = AlignmentWriter(stream)
        n = writer.write_file(alignments)
        self.assertEqual(n, 1)
        alignments = AlignmentIterator(path)
        alignments = list(alignments)
        alignment = alignments[0]
        stream.seek(0)
        saved_alignments = AlignmentIterator(stream)
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
                numpy.array_equal(alignment.coordinates, saved_alignment.coordinates)
            )

    def test_nexus1(self):
        path = "Nexus/test_Nexus_input.nex"
        with open(path) as stream:
            alignments = AlignmentIterator(stream)
            alignments = list(alignments)
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
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
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [
                        [
                            0,
                            1,
                            1,
                            2,
                            2,
                            3,
                            3,
                            4,
                            5,
                            6,
                            8,
                            12,
                            13,
                            14,
                            16,
                            16,
                            17,
                            17,
                            18,
                            18,
                            18,
                            18,
                            19,
                            20,
                            21,
                            23,
                            27,
                            28,
                            29,
                            31,
                            31,
                            32,
                            32,
                            33,
                        ],
                        [
                            0,
                            1,
                            1,
                            2,
                            2,
                            3,
                            4,
                            5,
                            6,
                            7,
                            9,
                            9,
                            9,
                            10,
                            12,
                            12,
                            13,
                            13,
                            14,
                            14,
                            14,
                            16,
                            17,
                            18,
                            19,
                            21,
                            21,
                            21,
                            22,
                            24,
                            24,
                            25,
                            25,
                            26,
                        ],
                        [
                            0,
                            1,
                            1,
                            2,
                            3,
                            4,
                            5,
                            6,
                            7,
                            8,
                            10,
                            14,
                            15,
                            16,
                            16,
                            16,
                            16,
                            16,
                            16,
                            16,
                            18,
                            20,
                            21,
                            22,
                            23,
                            25,
                            29,
                            30,
                            31,
                            31,
                            31,
                            31,
                            31,
                            31,
                        ],
                        [
                            0,
                            1,
                            1,
                            2,
                            2,
                            3,
                            3,
                            4,
                            4,
                            4,
                            4,
                            4,
                            4,
                            4,
                            4,
                            4,
                            4,
                            4,
                            4,
                            4,
                            4,
                            4,
                            4,
                            4,
                            4,
                            4,
                            4,
                            4,
                            4,
                            4,
                            4,
                            4,
                            4,
                            4,
                        ],
                        [
                            0,
                            1,
                            1,
                            2,
                            3,
                            4,
                            4,
                            5,
                            6,
                            6,
                            8,
                            12,
                            12,
                            13,
                            15,
                            15,
                            16,
                            17,
                            18,
                            18,
                            20,
                            20,
                            20,
                            21,
                            21,
                            23,
                            27,
                            27,
                            28,
                            30,
                            30,
                            31,
                            32,
                            33,
                        ],
                        [
                            0,
                            1,
                            2,
                            3,
                            4,
                            5,
                            6,
                            7,
                            8,
                            9,
                            9,
                            13,
                            14,
                            15,
                            17,
                            17,
                            18,
                            18,
                            19,
                            21,
                            23,
                            25,
                            26,
                            27,
                            28,
                            28,
                            32,
                            33,
                            34,
                            36,
                            36,
                            37,
                            37,
                            38,
                        ],
                        [
                            0,
                            1,
                            2,
                            3,
                            3,
                            4,
                            5,
                            6,
                            7,
                            8,
                            10,
                            14,
                            15,
                            16,
                            18,
                            18,
                            19,
                            19,
                            20,
                            22,
                            22,
                            24,
                            25,
                            26,
                            27,
                            29,
                            33,
                            34,
                            35,
                            37,
                            37,
                            38,
                            38,
                            39,
                        ],
                        [
                            0,
                            1,
                            2,
                            3,
                            4,
                            5,
                            6,
                            7,
                            8,
                            9,
                            11,
                            15,
                            16,
                            17,
                            19,
                            19,
                            20,
                            20,
                            21,
                            23,
                            25,
                            27,
                            28,
                            29,
                            30,
                            32,
                            36,
                            37,
                            38,
                            40,
                            40,
                            41,
                            41,
                            42,
                        ],
                        [
                            0,
                            1,
                            2,
                            3,
                            4,
                            5,
                            6,
                            7,
                            8,
                            9,
                            11,
                            15,
                            16,
                            17,
                            19,
                            20,
                            21,
                            21,
                            22,
                            24,
                            26,
                            28,
                            29,
                            30,
                            31,
                            33,
                            37,
                            38,
                            39,
                            41,
                            42,
                            43,
                            43,
                            44,
                        ],
                    ]
                ),
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
        self.check_reading_writing(path)

    def test_nexus2(self):
        path = "Nexus/codonposset.nex"
        with open(path) as stream:
            alignments = AlignmentIterator(stream)
            alignments = list(alignments)
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 22))
        self.assertEqual(alignment.sequences[0].id, "Aegotheles")
        self.assertEqual(alignment.sequences[1].id, "Aerodramus")
        self.assertEqual(alignment.sequences[0].annotations, {"molecule_type": "DNA"})
        self.assertEqual(alignment.sequences[1].annotations, {"molecule_type": "DNA"})
        self.assertEqual(alignment.sequences[0].seq, "AAAAAGGCATTGTGGTGGGAAT")
        self.assertEqual(alignment.sequences[1].seq, "?????????TTGTGGTGGGAAT")
        self.assertTrue(
            numpy.array_equal(alignment.coordinates, numpy.array([[0, 22], [0, 22]]))
        )
        self.assertEqual(alignment[0], "AAAAAGGCATTGTGGTGGGAAT")
        self.assertEqual(alignment[1], "?????????TTGTGGTGGGAAT")
        self.check_reading_writing(path)


class TestNexusBasic(unittest.TestCase):
    def test_empty(self):
        import io

        stream = io.StringIO()
        with self.assertRaisesRegex(ValueError, "Empty file."):
            AlignmentIterator(stream)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
