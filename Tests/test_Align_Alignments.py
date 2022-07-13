# Copyright 2022 by Michiel de Hoon.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for the Alignments class in Bio.Align.interfaces."""
import unittest
import warnings

from io import StringIO

from Bio import BiopythonExperimentalWarning

with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonExperimentalWarning)
    from Bio.Align.emboss import AlignmentIterator


try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align.emboss."
    ) from None


class TestEmboss(unittest.TestCase):
    def setUp(self):
        self.path = "Emboss/needle.txt"
        self.alignments = AlignmentIterator(self.path)

    def test_repr(self):
        representation = "<Bio.Align.emboss.AlignmentIterator object at %s>" % hex(id(self.alignments))
        self.assertFalse(self.alignments._loaded)
        self.assertEqual(repr(self.alignments), representation)
        alignments = self.alignments[::1]
        self.assertTrue(self.alignments._loaded)
        self.assertEqual(repr(self.alignments), representation)
        self.assertEqual(repr(alignments), representation)

    def test_lt_false(self):
        same_alignments = AlignmentIterator(self.path)
        self.assertFalse(self.alignments._loaded)
        self.assertFalse(same_alignments._loaded)
        result = (self.alignments < same_alignments)
        self.assertTrue(self.alignments._loaded)
        self.assertTrue(same_alignments._loaded)
        self.assertFalse(result)

    def test_lt_true(self):
        other_alignments = AlignmentIterator("Emboss/water.txt")
        self.assertFalse(self.alignments._loaded)
        self.assertFalse(other_alignments._loaded)
        result = (self.alignments < other_alignments)
        self.assertTrue(self.alignments._loaded)
        self.assertTrue(other_alignments._loaded)
        self.assertTrue(result)

    def test_le_false(self):
        other_alignments = AlignmentIterator("Emboss/water.txt")
        self.assertFalse(self.alignments._loaded)
        self.assertFalse(other_alignments._loaded)
        result = (other_alignments <= self.alignments)
        self.assertTrue(self.alignments._loaded)
        self.assertTrue(other_alignments._loaded)
        self.assertFalse(result)

    def test_le_true(self):
        same_alignments = AlignmentIterator(self.path)
        self.assertFalse(self.alignments._loaded)
        self.assertFalse(same_alignments._loaded)
        result = (self.alignments <= same_alignments)
        self.assertTrue(self.alignments._loaded)
        self.assertTrue(same_alignments._loaded)
        self.assertTrue(result)

    def test_eq_false(self):
        other_alignments = AlignmentIterator("Emboss/water.txt")
        self.assertFalse(self.alignments._loaded)
        self.assertFalse(other_alignments._loaded)
        result = (other_alignments == self.alignments)
        self.assertTrue(self.alignments._loaded)
        self.assertTrue(other_alignments._loaded)
        self.assertFalse(result)

    def test_eq_true(self):
        same_alignments = AlignmentIterator(self.path)
        self.assertFalse(self.alignments._loaded)
        self.assertFalse(same_alignments._loaded)
        result = (self.alignments == same_alignments)
        self.assertTrue(self.alignments._loaded)
        self.assertTrue(same_alignments._loaded)
        self.assertTrue(result)

    def test_gt_false(self):
        same_alignments = AlignmentIterator(self.path)
        self.assertFalse(self.alignments._loaded)
        self.assertFalse(same_alignments._loaded)
        result = (same_alignments > self.alignments)
        self.assertTrue(self.alignments._loaded)
        self.assertTrue(same_alignments._loaded)
        self.assertFalse(result)

    def test_gt_true(self):
        other_alignments = AlignmentIterator("Emboss/water.txt")
        self.assertFalse(self.alignments._loaded)
        self.assertFalse(other_alignments._loaded)
        result = (other_alignments > self.alignments)
        self.assertTrue(self.alignments._loaded)
        self.assertTrue(other_alignments._loaded)
        self.assertTrue(result)

    def test_ge_false(self):
        other_alignments = AlignmentIterator("Emboss/water.txt")
        self.assertFalse(self.alignments._loaded)
        self.assertFalse(other_alignments._loaded)
        result = (self.alignments >= other_alignments)
        self.assertTrue(self.alignments._loaded)
        self.assertTrue(other_alignments._loaded)
        self.assertFalse(result)

    def test_ge_true(self):
        same_alignments = AlignmentIterator(self.path)
        self.assertFalse(self.alignments._loaded)
        self.assertFalse(same_alignments._loaded)
        result = (self.alignments >= same_alignments)
        self.assertTrue(self.alignments._loaded)
        self.assertTrue(same_alignments._loaded)
        self.assertTrue(result)

    def test_contains_false(self):
        other_alignments = AlignmentIterator("Emboss/water.txt")
        alignment = next(other_alignments)
        self.assertFalse(self.alignments._loaded)
        result = (alignment in self.alignments)
        self.assertTrue(self.alignments._loaded)
        self.assertFalse(result)

    def test_contains_true(self):
        same_alignments = AlignmentIterator(self.path)
        alignment = next(same_alignments)
        self.assertFalse(self.alignments._loaded)
        result = (alignment in self.alignments)
        self.assertTrue(self.alignments._loaded)
        self.assertTrue(result)

    def test_len(self):
        self.assertFalse(self.alignments._loaded)
        length = len(self.alignments)
        self.assertTrue(self.alignments._loaded)
        self.assertEqual(length, 5)

    def test_getitem(self):
        self.assertFalse(self.alignments._loaded)
        alignment = self.alignments[3]
        self.assertTrue(self.alignments._loaded)
        same_alignments = AlignmentIterator(self.path)
        zeroth_alignment = next(same_alignments)
        self.assertNotEqual(zeroth_alignment, alignment)
        first_alignment = next(same_alignments)
        self.assertNotEqual(first_alignment, alignment)
        second_alignment = next(same_alignments)
        self.assertNotEqual(second_alignment, alignment)
        third_alignment = next(same_alignments)
        self.assertEqual(third_alignment, alignment)
        fourth_alignment = next(same_alignments)
        self.assertNotEqual(fourth_alignment, alignment)
        with self.assertRaises(StopIteration):
            next(same_alignments)
        self.assertFalse(same_alignments._loaded)

    def test_setitem(self):
        same_alignments = AlignmentIterator(self.path)
        zeroth_alignment = next(same_alignments)
        first_alignment = next(same_alignments)
        second_alignment = next(same_alignments)
        third_alignment = next(same_alignments)
        fourth_alignment = next(same_alignments)
        self.assertFalse(same_alignments._loaded)
        self.assertFalse(self.alignments._loaded)
        self.alignments[3] = first_alignment
        self.assertTrue(self.alignments._loaded)
        alignment = next(self.alignments)
        self.assertEqual(zeroth_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(first_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(second_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(first_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(fourth_alignment, alignment)
        with self.assertRaises(StopIteration):
            next(self.alignments)

    def test_delitem(self):
        same_alignments = AlignmentIterator(self.path)
        zeroth_alignment = next(same_alignments)
        first_alignment = next(same_alignments)
        second_alignment = next(same_alignments)
        third_alignment = next(same_alignments)
        fourth_alignment = next(same_alignments)
        self.assertFalse(same_alignments._loaded)
        self.assertFalse(self.alignments._loaded)
        del self.alignments[3]
        self.assertTrue(self.alignments._loaded)
        alignment = next(self.alignments)
        self.assertEqual(zeroth_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(first_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(second_alignment, alignment)
        alignment = next(self.alignments)
        self.assertEqual(fourth_alignment, alignment)
        with self.assertRaises(StopIteration):
            next(self.alignments)

    def test_add(self):
        same_repr = ("<Bio.Align.Alignment object (2 rows x 124 columns) at ",
                     "<Bio.Align.Alignment object (2 rows x 119 columns) at ",
                     "<Bio.Align.Alignment object (2 rows x 120 columns) at ",
                     "<Bio.Align.Alignment object (2 rows x 118 columns) at ",
                     "<Bio.Align.Alignment object (2 rows x 125 columns) at ",
                    )
        other_repr = ("<Bio.Align.Alignment object (2 rows x 131 columns) at ",
                     )
        same_alignments = AlignmentIterator(self.path)
        self.assertFalse(self.alignments._loaded)
        self.assertFalse(same_alignments._loaded)
        alignments = self.alignments + same_alignments
        self.assertTrue(self.alignments._loaded)
        self.assertTrue(same_alignments._loaded)
        for alignment, representation in zip(alignments, 2 * same_repr):
            self.assertTrue(repr(alignment).startswith(representation))
        self.assertEqual(len(alignments), 10)
        other_alignments = AlignmentIterator("Emboss/water.txt")
        self.assertFalse(other_alignments._loaded)
        alignments = self.alignments + other_alignments
        self.assertTrue(other_alignments._loaded)
        for alignment, representation in zip(alignments, same_repr + other_repr):
            self.assertTrue(repr(alignment).startswith(representation))
        self.assertEqual(len(alignments), 6)

    def test_needle(self):
        alignments = self.alignments
        self.assertEqual(alignments.metadata["Program"], "needle")
        self.assertEqual(alignments.metadata["Rundate"], "Sun 27 Apr 2007 17:20:35")
        self.assertEqual(
            alignments.metadata["Command line"],
            "needle [-asequence] Spo0F.faa [-bsequence] paired_r.faa -sformat2 pearson",
        )
        self.assertEqual(alignments.metadata["Align_format"], "srspair")
        self.assertEqual(alignments.metadata["Report_file"], "ref_rec .needle")
        alignment = next(alignments)
        self.assertEqual(alignment.annotations["matrix"], "EBLOSUM62")
        self.assertAlmostEqual(alignment.annotations["gap_penalty"], 10.0)
        self.assertAlmostEqual(alignment.annotations["extend_penalty"], 0.5)
        self.assertEqual(alignment.annotations["identity"], 32)
        self.assertEqual(alignment.annotations["similarity"], 64)
        self.assertEqual(alignment.annotations["gaps"], 17)
        self.assertAlmostEqual(alignment.annotations["score"], 112.0)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 124))
        self.assertEqual(alignment.sequences[0].id, "ref_rec")
        self.assertEqual(alignment.sequences[1].id, "gi|94968718|receiver")
        self.assertEqual(
            alignment.sequences[0].seq,
            "KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIVTKERPDLVLLDMKIPGMDGIEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHFAKPFDIDEIRDAV",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "VLLADDHALVRRGFRLMLEDDPEIEIVAEAGDGAQAVKLAGELHPRVVVMDCAMPGMSGMDATKQIRTQWPDIAVLMLTMHSEDTWVRLALEAGANGYILKSAIDLDLIQAVRRVANGET",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [
                        [0, 1, 7, 7, 17, 19, 100, 100, 108, 109, 111, 111],
                        [0, 0, 6, 10, 20, 20, 101, 102, 110, 110, 112, 120],
                    ]
                ),
            )
        )
        self.assertEqual(
            alignment[0],
            "KILIVDD----QYGIRILLNEVFNKEGYQTFQAANGLQALDIVTKERPDLVLLDMKIPGMDGIEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHFAK-PFDIDEIRDAV--------",
        )
        self.assertEqual(
            alignment[1],
            "-VLLADDHALVRRGFRLMLED--DPEIEIVAEAGDGAQAVKLAGELHPRVVVMDCAMPGMSGMDATKQIRTQWPDIAVLMLTMHSEDTWVRLALEAGANGYILKSAIDLDLIQ-AVRRVANGET",
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            " :|:.||    :.|.|::|.:  :.|.....:|.:|.||:.:..:..|.:|::|..:|||.|::..|:::....:|.|:::|.:.|...::.:.|.||..:..| ..|:|.|: ||        ",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.annotations["matrix"], "EBLOSUM62")
        self.assertAlmostEqual(alignment.annotations["gap_penalty"], 10.0)
        self.assertAlmostEqual(alignment.annotations["extend_penalty"], 0.5)
        self.assertEqual(alignment.annotations["identity"], 34)
        self.assertEqual(alignment.annotations["similarity"], 58)
        self.assertEqual(alignment.annotations["gaps"], 9)
        self.assertAlmostEqual(alignment.annotations["score"], 154.0)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 119))
        self.assertEqual(alignment.sequences[0].id, "ref_rec")
        self.assertEqual(alignment.sequences[1].id, "gi|94968761|receiver")
        self.assertEqual(
            alignment.sequences[0].seq,
            "KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIVTKERPDLVLLDMKIPGMDGIEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHFAKPFDIDEIRDAV",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "ILIVDDEANTLASLSRAFRLAGHEATVCDNAVRALEIAKSKPFDLILSDVVMPGRDGLTLLEDLKTAGVQAPVVMMSGQAHIEMAVKATRLGALDFLEKPLSTDKLLLTVENALKLKR",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[0, 1, 111, 111], [0, 0, 110, 118]])
            )
        )
        self.assertEqual(
            alignment[0],
            "KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIVTKERPDLVLLDMKIPGMDGIEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHFAKPFDIDEIRDAV--------",
        )
        self.assertEqual(
            alignment[1],
            "-ILIVDDEANTLASLSRAFRLAGHEATVCDNAVRALEIAKSKPFDLILSDVVMPGRDGLTLLEDLKTAGVQAPVVMMSGQAHIEMAVKATRLGALDFLEKPLSTDKLLLTVENALKLKR",
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            " ||||||:......|:..|...|::.....|.::||:|...:..||:|.|:.:||.||:.:|:.:|.......|::|:....::|..::..||||....||...|::...|        ",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.annotations["matrix"], "EBLOSUM62")
        self.assertAlmostEqual(alignment.annotations["gap_penalty"], 10.0)
        self.assertAlmostEqual(alignment.annotations["extend_penalty"], 0.5)
        self.assertEqual(alignment.annotations["identity"], 29)
        self.assertEqual(alignment.annotations["similarity"], 53)
        self.assertEqual(alignment.annotations["gaps"], 9)
        self.assertAlmostEqual(alignment.annotations["score"], 121.0)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 120))
        self.assertEqual(alignment.sequences[0].id, "ref_rec")
        self.assertEqual(alignment.sequences[1].id, "gi|94967506|receiver")
        self.assertEqual(
            alignment.sequences[0].seq,
            "KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIVTKERPDLVLLDMKIPGMDGIEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHFAKPFDIDEIRDAV",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "LHIVVVDDDPGTCVYIESVFAELGHTCKSFVRPEAAEEYILTHPVDLAIVDVYLGSTTGVEVLRRCRVHRPKLYAVIITGQISLEMAARSIAEGAVDYIQKPIDIDALLNIAERALEHKE",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[0, 0, 111, 111], [0, 1, 112, 120]])
            )
        )
        self.assertEqual(
            alignment[0],
            "-KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIVTKERPDLVLLDMKIPGMDGIEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHFAKPFDIDEIRDAV--------",
        )
        self.assertEqual(
            alignment[1],
            "LHIVVVDDDPGTCVYIESVFAELGHTCKSFVRPEAAEEYILTHPVDLAIVDVYLGSTTGVEVLRRCRVHRPKLYAVIITGQISLEMAARSIAEGAVDYIQKPIDIDALLNIAERALEHKE",
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            " .|::|||..|..:.:..||.:.|:..........|.:.:.....||.::|:.:....|:|:|:|.:|....:..:|:|....|:|...|...||:.:..||.|||.:.:..        ",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.annotations["matrix"], "EBLOSUM62")
        self.assertAlmostEqual(alignment.annotations["gap_penalty"], 10.0)
        self.assertAlmostEqual(alignment.annotations["extend_penalty"], 0.5)
        self.assertEqual(alignment.annotations["identity"], 30)
        self.assertEqual(alignment.annotations["similarity"], 64)
        self.assertEqual(alignment.annotations["gaps"], 9)
        self.assertAlmostEqual(alignment.annotations["score"], 126.0)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 118))
        self.assertEqual(alignment.sequences[0].id, "ref_rec")
        self.assertEqual(alignment.sequences[1].id, "gi|94970045|receiver")
        self.assertEqual(
            alignment.sequences[0].seq,
            "KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIVTKERPDLVLLDMKIPGMDGIEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHFAKPFDIDEIRDAV",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "VLLVEDEEALRAAAGDFLETRGYKIMTARDGTEALSMASKFAERIDVLITDLVMPGISGRVLAQELVKIHPETKVMYMSGYDDETVMVNGEIDSSSAFLRKPFRMDALSAKIREVL",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [
                        [0, 1, 41, 41, 82, 83, 98, 98, 105, 105, 111],
                        [0, 0, 40, 42, 83, 83, 98, 99, 106, 110, 116],
                    ]
                ),
            )
        )
        self.assertEqual(
            alignment[0],
            "KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIVTK--ERPDLVLLDMKIPGMDGIEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHF-AKPFDID----EIRDAV",
        )
        self.assertEqual(
            alignment[1],
            "-VLLVEDEEALRAAAGDFLETRGYKIMTARDGTEALSMASKFAERIDVLITDLVMPGISGRVLAQELVKIHPETKVMYMSGYDD-ETVMVNGEIDSSSAFLRKPFRMDALSAKIREVL",
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            " :|:|:|:..:|....:.....||:...|.:|.:||.:.:|  ||.|:::.|:.:||:.|..:.:.:..|....:|:.|:.|.: :.:..:.|:.:.:.| .|||.:|    :||:.:",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.annotations["matrix"], "EBLOSUM62")
        self.assertAlmostEqual(alignment.annotations["gap_penalty"], 10.0)
        self.assertAlmostEqual(alignment.annotations["extend_penalty"], 0.5)
        self.assertEqual(alignment.annotations["identity"], 35)
        self.assertEqual(alignment.annotations["similarity"], 70)
        self.assertEqual(alignment.annotations["gaps"], 18)
        self.assertAlmostEqual(alignment.annotations["score"], 156.5)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 125))
        self.assertEqual(alignment.sequences[0].id, "ref_rec")
        self.assertEqual(alignment.sequences[1].id, "gi|94970041|receiver")
        self.assertEqual(
            alignment.sequences[0].seq,
            "KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIVTKERPDLVLLDMKIPGMDGIEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHFAKPFDIDEIRDAV",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "TVLLVEDEEGVRKLVRGILSRQGYHVLEATSGEEALEIVRESTQKIDMLLSDVVLVGMSGRELSERLRIQMPSLKVIYMSGYTDDAIVRHGVLTESAEFLQKPFTSDSLLRKVRAVLQKRQ",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [
                        [0, 39, 39, 88, 92, 99, 99, 111, 111],
                        [0, 39, 41, 90, 90, 97, 101, 113, 121],
                    ]
                ),
            )
        )
        self.assertEqual(
            alignment[0],
            "KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIV--TKERPDLVLLDMKIPGMDGIEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHFA----KPFDIDEIRDAV--------",
        )
        self.assertEqual(
            alignment[1],
            "TVLLVEDEEGVRKLVRGILSRQGYHVLEATSGEEALEIVRESTQKIDMLLSDVVLVGMSGRELSERLRIQMPSLKVIYMSGYTDDAIVRH----GVLTESAEFLQKPFTSDSLLRKVRAVLQKRQ",
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            ".:|:|:|:.|:|.|:..:.:::||...:|.:|.:||:||  :.::.|::|.|:.:.||.|.|:.:|:::...:::||.|:.|.:..:::.    |.||..|    |||..|.:...|        ",
        )
        with self.assertRaises(StopIteration):
            next(alignments)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
