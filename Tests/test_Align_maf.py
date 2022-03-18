# Copyright 2021 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Align.maf module."""
import unittest
import warnings
from io import StringIO


from Bio import BiopythonExperimentalWarning

with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonExperimentalWarning)
    from Bio.Align import maf


try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align.maf."
    ) from None


class TestAlign_reading(unittest.TestCase):
    def test_reading_bundle_without_target(self):
        """Test parsing bundle_without_target.maf."""
        path = "MAF/bundle_without_target.maf"
        alignments = maf.AlignmentIterator(path)
        self.assertEqual(alignments.metadata["version"], "1")
        self.assertEqual(alignments.metadata["scoring"], "autoMZ.v1")
        alignment = next(alignments)
        self.assertRaises(StopIteration, next, alignments)
        self.assertEqual(alignment.score, 6441)
        self.assertEqual(len(alignment.sequences), 2)
        self.assertEqual(alignment.sequences[0].id, "mm8.chr10")
        self.assertEqual(alignment.sequences[1].id, "oryCun1.scaffold_133159")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(len(alignment.sequences[1].seq), 13221)
        self.assertEqual(
            alignment.sequences[0].seq[3009319 : 3009319 + 162],
            "TCATAGGTATTTATTTTTAAATATGGTTTGCTTTATGGCTAGAACACACCGATTACTTAAAATAGGATTAACCCCCATACACTTTAAAAATGATTAAACAACATTTCTGCTGCTCGCTCACATTCTTCATAGAAGATGACATAATGTATTTTCCTTTTGGTT",
        )
        self.assertEqual(
            alignment.sequences[1].seq[11087 : 11087 + 164],
            "TCACAGATATTTACTATTAAATATGGTTTGTTATATGGTTACGGTTCATAGGTTACTTGGAATTGGATTAACCTTCTTATTCATTGCAGAATTGGTTACACTGTGTTCTTGACCTTTGCTTGTTTTCTCCATGGAAACTGATGTCAAATACTTTCCCTTTGGTT",
        )
        self.assertEqual(
            alignment[0],
            "TCATAGGTATTTATTTTTAAATATGGTTTGCTTTATGGCTAGAACACACCGATTACTTAAAATAGGATTAACC--CCCATACACTTTAAAAATGATTAAACAACATTTCTGCTGCTCGCTCACATTCTTCATAGAAGATGACATAATGTATTTTCCTTTTGGTT",
        )
        self.assertEqual(
            alignment[1],
            "TCACAGATATTTACTATTAAATATGGTTTGTTATATGGTTACGGTTCATAGGTTACTTGGAATTGGATTAACCTTCTTATTCATTGCAGAATTGGTTACACTGTGTTCTTGACCTTTGCTTGTTTTCTCCATGGAAACTGATGTCAAATACTTTCCCTTTGGTT",
        )
        self.assertEqual(
            alignment.sequences[1].annotations["quality"],
            "99569899999998999999999999999999999999999999999999999999999999999999999757878999975999999999999999979999999999997899999999999997997999999869999996999988997997999999",
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "N")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "N")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 0)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[3009319, 3009392, 3009392, 3009481],
                             [  11087,   11160,   11162,   11251],
                            ])
                # fmt: on
            )
        )

    def test_reading_ucsc_mm9_chr10(self):
        """Test parsing MAF file ucsc_mm9_chr10.maf."""
        path = "MAF/ucsc_mm9_chr10.maf"
        alignments = maf.AlignmentIterator(path)
        self.assertEqual(alignments.metadata["version"], "1")
        self.assertEqual(alignments.metadata["scoring"], "autoMZ.v1")
        alignment = next(alignments)
        self.assertEqual(alignment.score, 6441)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3009319 : 3009319 + 162],
            "TCATAGGTATTTATTTTTAAATATGGTTTGCTTTATGGCTAGAACACACCGATTACTTAAAATAGGATTAACCCCCATACACTTTAAAAATGATTAAACAACATTTCTGCTGCTCGCTCACATTCTTCATAGAAGATGACATAATGTATTTTCCTTTTGGTT",
        )
        self.assertEqual(
            alignment[0],
            "TCATAGGTATTTATTTTTAAATATGGTTTGCTTTATGGCTAGAACACACCGATTACTTAAAATAGGATTAACC--CCCATACACTTTAAAAATGATTAAACAACATTTCTGCTGCTCGCTCACATTCTTCATAGAAGATGACATAATGTATTTTCCTTTTGGTT",
        )
        self.assertEqual(alignment.sequences[1].id, "oryCun1.scaffold_133159")
        self.assertEqual(len(alignment.sequences[1].seq), 13221)
        self.assertEqual(
            alignment.sequences[1].seq[11087 : 11087 + 164],
            "TCACAGATATTTACTATTAAATATGGTTTGTTATATGGTTACGGTTCATAGGTTACTTGGAATTGGATTAACCTTCTTATTCATTGCAGAATTGGTTACACTGTGTTCTTGACCTTTGCTTGTTTTCTCCATGGAAACTGATGTCAAATACTTTCCCTTTGGTT",
        )
        self.assertEqual(
            alignment[1],
            "TCACAGATATTTACTATTAAATATGGTTTGTTATATGGTTACGGTTCATAGGTTACTTGGAATTGGATTAACCTTCTTATTCATTGCAGAATTGGTTACACTGTGTTCTTGACCTTTGCTTGTTTTCTCCATGGAAACTGATGTCAAATACTTTCCCTTTGGTT",
        )
        self.assertEqual(
            alignment.sequences[1].annotations["quality"],
            "99569899999998999999999999999999999999999999999999999999999999999999999757878999975999999999999999979999999999997899999999999997997999999869999996999988997997999999",
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "N")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "N")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 0)
        self.assertEqual(len(alignment.sequences), 2)
        self.assertNotIn("empty", alignment.annotations)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
                numpy.array([[3009319, 3009392, 3009392, 3009481],
                             [  11087,   11160,   11162,   11251],
                            ])
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 103072)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3012076 : 3012076 + 365],
            "AGTCTTTCCAATGGGACCTGTGAGTCCTAACTATGCCAGCACTCCCAACAGCAAGACACTAAGTTCACTCATCCTTGGTGGATGGGATTTTGCTCCTGGAGTGTCACCAAATTAAATAACCAGTGAGCAGAGTTGTGACGAGCATCAGGCTCTGGATTTAGGTGAGAGACCTTAGTGTATGTCTCCTGTAGGTCGCAGCTCCCTATGGATGAGTCAAGTGAAGGTCCTGAGACAACAAGTCCTCGGCTATGTGGGGGTGAGGGATGCAGCTGGAACCTCAGGGATCTCTGTAAGCAGTGGCATAAATGCTTGGCGGGAGAGAGCATGTTAGAGCTCACACGACATAGGAAGCCACTGAGACACTG",
        )
        self.assertEqual(
            alignment[0],
            "AGTCTTTCCAATGGGACCTGTGAGTCCTAACTATGCCAGC-----ACTCCCAACAGCAAGACACTAAGTT---------CACTCATCCTTGGTGGATGGGATTTTGCTCCTGGAGTGTCAC-----CAAATTAAATAACCAGTGAGCAGAGTTG--TGACGAGCATCAGGCTCTGGATTTAGGTGAGAGACCTTAGTGTATGTCTCCTGTAGGTCGCAGCTCCCTATGGAT--------------------------GAGTCAAGTGAAGGTCCTGAGACAA-------------------CAAGTCCTC----GGCTATGTGGGGGTGAGGG-------------ATGC----AG--------CTGGAACCTCAGGGA-TCTCTGT-AAGCAGTGGCATAAATGCTTGGCGG--GAGAGAGCATGTTAGAGCTCACACGACATAGGAAGCCACTGA--GACACTG",
        )
        self.assertEqual(alignment.sequences[1].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[1].seq), 174210431)
        self.assertEqual(
            alignment.sequences[1].seq[158049785 : 158049785 + 443],
            "GAACACAAGTCCATGGCCTTTTCAACCACGTTAGGCGTTGGAGGCCACGTCAGCCCCCACAGATCCACTCTACAACTTTAGAAAGAATTTGTTAGGCACCAGTGGCTTGCTTCAGAGTTCTAGTTCCTACGCTGTCAGGGAGTCTACAGCAGCCCTATCCACAACTTACCTCCAGGATTTTGCTTAATTTCTGTGGAAAGCGATCTCAACTTCCCAGTTCACACACAGATGCAACCAGCATGAGACATTCCTAGAAATCAACAACCTTTGGTGACTGGTCTATCTCATGCTATTCTAAATACTGCTTGCTCACTTCATAAGTACGACATTCCAAGAGCACAAGGCTATCTGGTAAGGATGAGTATGGAGATAAATATTAACATCTTTTTACTGTGATGATTTGGCTGGAATAATTAAAACTTATATTTCCACTTATGAAGACT",
        )
        self.assertEqual(
            alignment[1],
            "AGTCTTCATAAGTGGAAATATAAGTTTTAATTATTCCAGCCAAATCATCACAGTAAAAAGATGTTAATATTTATCTCCATACTCATCCTTACCAGATAGCCTTGTGCTCTTGGAATGTCGTACTTATGAAGTGAGCAAGCAGTATTTAGAATAGCATGAGATAGACCAGTCACCAA----AGGTTGTTGATTTCTAGGAATGTCTCATGCTGGTTGCATCTGTGTGTGAACTGGGAAGTTGAGATCGCTTTCCACAGAAATTAAGCAAAATCCTGGAGGTAAGTTGTGGA-----------TAGGGCTGCTGTAGACTCCCTGACAGCGTAGGAACT--------AGAACTCTGAAGCAAGCCACTGGTGCCTAACAAATTCTTTCTAAAGTTGTAGAGTGGATCTGTGGGGGCTGACGTGGCCTCCAACGCCTAACGTGGTTGAAAAGGCCATGGACTTGTGTTC",
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "N")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[2].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[2].seq), 173908612)
        self.assertEqual(
            alignment.sequences[2].seq[157529165 : 157529165 + 443],
            "GAACACAAGTCCATGGCCTTTTCAACCATGTTAGGCGTTGGAGGCCACGTCAGCCCCCACAGATCCACTCTACAACTTTATAAAGAATTTGTTAGGCACCAGTGGCTTGCTTCAGAGTTCTAGTTCCTATGCTGTCAGGGAGTCTACAGCATCCCTATCCACAACTTACCTCCAGGATTTTGCTTAATTTCTGTGGAAAGCGATCTCAACTTCCCAGTTCACACACAGATGCAACCAGCATGAGACATTCCTAGAAATCAACAACCTTTGGTGACTGGTCTATCTCATGCTATTCTACATACTGCTTGCTCACTTCACAAATACGACATTCCAAGAGCACAAGGCTATCTGGTAAGGATGAGTATGGAGATAAATATTAACATCTTTTTACTGTGATGATTTAGCTGGAATAATTAAAACTTATATTTCCACTTATGAAGACT",
        )
        self.assertEqual(
            alignment[2],
            "AGTCTTCATAAGTGGAAATATAAGTTTTAATTATTCCAGCTAAATCATCACAGTAAAAAGATGTTAATATTTATCTCCATACTCATCCTTACCAGATAGCCTTGTGCTCTTGGAATGTCGTATTTGTGAAGTGAGCAAGCAGTATGTAGAATAGCATGAGATAGACCAGTCACCAA----AGGTTGTTGATTTCTAGGAATGTCTCATGCTGGTTGCATCTGTGTGTGAACTGGGAAGTTGAGATCGCTTTCCACAGAAATTAAGCAAAATCCTGGAGGTAAGTTGTGGA-----------TAGGGATGCTGTAGACTCCCTGACAGCATAGGAACTAGAACTCTGAAGC----AAGCCA----CTGGTGCCTAACAAATTCTTTATAAAGTTGTAGAGTGGATCTGTGGGGGCTGACGTGGCCTCCAACGCCTAACATGGTTGAAAAGGCCATGGACTTGTGTTC",
        )
        self.assertEqual(
            alignment.sequences[2].annotations["quality"],
            "99999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[2].annotations["leftStatus"], "N")
        self.assertEqual(alignment.sequences[2].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[2].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[3].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[3].seq), 170899992)
        self.assertEqual(
            alignment.sequences[3].seq[155039093 : 155039093 + 443],
            "GAACACAAGTCCATGGCCTTTTCAACCATGTTAGGCGTTGGAGGCCACGTCAGCCCCCATGGATCCACTCTACAACTTTATAAAGAATTTGTTAGGCACCAGTGGCTTGCTTCAGAGTTCTAGTTCCTATGCTGTCAGGGAGTCTACAGCATCCCTATCCACAACTTACCTCCAGGATTTTGCTTAATTTCTGTGGAAAGCGATCTCAACTTCCCAGTTCACACACAGATGCAACCAGCATGAGACATTCCTAGAAATCAACAACCTTTGGTGACTGGTCTATCTCATGCTATTCTACATACTGCTTGCTCACTTCACAAATACGACATTCCAAGAGCACAAGGCTATCTGGTAAGGATGAGTATGGAGATAAATATTAACATCTTTTTACTGTGATGATTTAGCTGGAATAATTAAAACTTATATTTCCACTTATGAAGACT",
        )
        self.assertEqual(
            alignment[3],
            "AGTCTTCATAAGTGGAAATATAAGTTTTAATTATTCCAGCTAAATCATCACAGTAAAAAGATGTTAATATTTATCTCCATACTCATCCTTACCAGATAGCCTTGTGCTCTTGGAATGTCGTATTTGTGAAGTGAGCAAGCAGTATGTAGAATAGCATGAGATAGACCAGTCACCAA----AGGTTGTTGATTTCTAGGAATGTCTCATGCTGGTTGCATCTGTGTGTGAACTGGGAAGTTGAGATCGCTTTCCACAGAAATTAAGCAAAATCCTGGAGGTAAGTTGTGGATAGGGATGCTGTAGACTCCCT---GACAGCATAGGAAC-TAGAACTCTGAAGC---AAGC----CA--------CTGGTGCCTAACAAATTCTTTATAAAGTTGTAGAGTGGATCCATGGGGGCTGACGTGGCCTCCAACGCCTAACATGGTTGAAAAGGCCATGGACTTGTGTTC",
        )
        self.assertEqual(alignment.sequences[3].annotations["leftStatus"], "N")
        self.assertEqual(alignment.sequences[3].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[3].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[3].annotations["rightCount"], 0)
        self.assertEqual(len(alignment.sequences), 4)
        self.assertNotIn("empty", alignment.annotations)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
            [[  3012076,   3012116,   3012116,   3012141,   3012141,   3012183,
                3012183,   3012211,   3012211,   3012231,   3012235,   3012286,
                3012286,   3012311,   3012311,   3012311,   3012320,   3012320,
                3012320,   3012334,   3012335,   3012339,   3012339,   3012339,
                3012339,   3012339,   3012343,   3012343,   3012345,   3012345,
                3012345,   3012360,   3012360,   3012367,   3012367,   3012392,
                3012392,   3012434,   3012434,   3012441],
             [158050228, 158050188, 158050183, 158050158, 158050149, 158050107,
              158050102, 158050074, 158050072, 158050052, 158050052, 158050001,
              158049975, 158049950, 158049942, 158049942, 158049933, 158049932,
              158049929, 158049915, 158049914, 158049910, 158049906, 158049906,
              158049906, 158049905, 158049901, 158049897, 158049895, 158049891,
              158049887, 158049872, 158049871, 158049864, 158049863, 158049838,
              158049836, 158049794, 158049792, 158049785],
             [157529608, 157529568, 157529563, 157529538, 157529529, 157529487,
              157529482, 157529454, 157529452, 157529432, 157529432, 157529381,
              157529355, 157529330, 157529322, 157529322, 157529313, 157529312,
              157529309, 157529295, 157529294, 157529290, 157529286, 157529280,
              157529278, 157529277, 157529273, 157529273, 157529271, 157529267,
              157529267, 157529252, 157529251, 157529244, 157529243, 157529218,
              157529216, 157529174, 157529172, 157529165],
             [155039536, 155039496, 155039491, 155039466, 155039457, 155039415,
              155039410, 155039382, 155039380, 155039360, 155039360, 155039309,
              155039283, 155039258, 155039250, 155039239, 155039230, 155039229,
              155039229, 155039215, 155039215, 155039211, 155039207, 155039201,
              155039201, 155039201, 155039197, 155039197, 155039195, 155039195,
              155039195, 155039180, 155039179, 155039172, 155039171, 155039146,
              155039144, 155039102, 155039100, 155039093],
             ]
                    # fmt: on
                ),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 49128)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3012441 : 3012441 + 125],
            "TGGGTCCCCTTGGCACATCCAGATCTCCCCAGTTAACCTGTCCTGCTTAGACCACTTACCTGAATTGAATTGGGAGGAGAGAAAGAAGCCAGTTTCCCAGAGAGGGAAAAGGAAAAGCTCGACAC",
        )
        self.assertEqual(
            alignment[0],
            "TGGGTCCCCTTGGCACATCCAGATCTCCCCAGTTAACCTGTCCTGCTTAGACCACTTACCTGAATTG--AATTGGGAGGAGAGAAAGAAGCCAGTTTCCCAGAGAGGGAAAAGGAAAAGCTCGACAC",
        )
        self.assertEqual(alignment.sequences[1].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[1].seq), 170899992)
        self.assertEqual(
            alignment.sequences[1].seq[155038979 : 155038979 + 114],
            "GCTTCCAGCATTTTCCATTGGCCCTCAGCAAGGTGCTCCCCAAGATCAGCGGTCAAGAGAGGTGGTCAATACACACCAGGTTACGTGAGGACTTGGTTATTCTAGAGGAACCCA",
        )
        self.assertEqual(
            alignment[1],
            "TGGGTTCCTCTAGAATAACCAAG--TCCTCACGTAACCTGGTGTGTATTGACCACCTCTCTTGACCGCTGATCTTGGGGAG----------CACCTTGCT-GAGGGCCAATGGAAAATGCTGGAAGC",
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[2].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[2].seq), 173908612)
        self.assertEqual(
            alignment.sequences[2].seq[157529051 : 157529051 + 114],
            "GCTTCCAGCATTTTCCATTGGCCCTCAGCAAGGTGCTCCCCAAGATCAGCGGTCAAGAGAGGTGGTCAATACACACCAGGTTACGTGAGGACTTGGTTATTCTAGAGGAACCCA",
        )
        self.assertEqual(
            alignment[2],
            "TGGGTTCCTCTAGAATAACCAAG--TCCTCACGTAACCTGGTGTGTATTGACCACCTCTCTTGACCGCTGATCTTGGGGAG----------CACCTTGCT-GAGGGCCAATGGAAAATGCTGGAAGC",
        )
        self.assertEqual(
            alignment.sequences[2].annotations["quality"],
            "999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[2].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[2].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[3].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[3].seq), 174210431)
        self.assertEqual(
            alignment.sequences[3].seq[158049671 : 158049671 + 114],
            "GCTTCCAGCATTTTCCATTGGCCCTCAGCAAGGTGCTCCCCAAGATCAACAGTCAAGACAGGTGGTCAATATAGACCAGGTTACGTGAGGACTTGGTTATTCTAGAGGAACCCA",
        )
        self.assertEqual(
            alignment[3],
            "TGGGTTCCTCTAGAATAACCAAG--TCCTCACGTAACCTGGTCTATATTGACCACCTGTCTTGACTGTTGATCTTGGGGAG----------CACCTTGCT-GAGGGCCAATGGAAAATGCTGGAAGC",
        )
        self.assertEqual(alignment.sequences[3].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[3].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[3].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[3].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[4].id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(alignment.sequences[4].seq), 359464)
        self.assertEqual(
            alignment.sequences[4].seq[178838 : 178838 + 101],
            "GCTTCCCGTGTTTGCTGTTGCCCCTCAGCAAAGGGCTCCCCAAGTTTGGTGGTCAGTGTAAACCCGGGTGAGTGAGGACTTGGCCAGGTCAGGGGGGCCCA",
        )
        self.assertEqual(
            alignment[4],
            "TGGGCCCCCCTGACCTGGCCAAG--TCCTCACTCACCCGGGTTTACACTGACCACC-------------AAACTTGGGGAG----------CCCTTTGCT-GAGGGGCAACAGCAAACACGGGAAGC",
        )
        self.assertEqual(
            alignment.sequences[4].annotations["quality"],
            "99999999999999999999999999999999999999999999999999999999999999999999999999999999979999999999999999999",
        )
        self.assertEqual(alignment.sequences[4].annotations["leftStatus"], "N")
        self.assertEqual(alignment.sequences[4].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[4].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[4].annotations["rightCount"], 0)
        self.assertEqual(len(alignment.sequences), 5)
        self.assertNotIn("empty", alignment.annotations)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
            [[  3012441,   3012464,   3012466,   3012497,   3012508,   3012508,
                3012520,   3012530,   3012539,   3012540,   3012566],
             [155039093, 155039070, 155039070, 155039039, 155039028, 155039026,
              155039014, 155039014, 155039005, 155039005, 155038979],
             [157529165, 157529142, 157529142, 157529111, 157529100, 157529098,
              157529086, 157529086, 157529077, 157529077, 157529051],
             [158049785, 158049762, 158049762, 158049731, 158049720, 158049718,
              158049706, 158049706, 158049697, 158049697, 158049671],
             [   178939,    178916,    178916,    178885,    178885,    178885,
                 178873,    178873,    178864,    178864,    178838],
            ]
                    # fmt: on
                ),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 117109)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3012566 : 3012566 + 262],
            "TGTGGGCTCCTCACTTTCCTGTCTCAGGTGTGTCTGTGAGTTTCGGTGAGTGTCGTACAGGAAAGAGGGTGAAAACTCAGTCTGAGCTGTCATTCTTGCCAGCTATGTTGCTTTCCTGTCCTCTTTAGCTTATCTCAGGCAACCTATCTTATTTTGTTTGCTTTCAGAAGGCAAGCGAtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtatgtgtgtgtgcgtgCGCGCGCGCGAGCACATGTGCATGCATGCGCACTCGTG",
        )
        self.assertEqual(
            alignment[0],
            "--TGTGGGCTCCTCACTTTCCTG-TCTCAGGTGTGTCTGTGAGTTTCGGTGAGTGTCGTACAGGAAAGAGGGTGAAAACTCAGTCTGAGCTGTCATTCTTGCCAGCTATGTTGCTTTCCTGTCCTCTTTAGC-------TTATCTCAGGCAACCTATCTTATTTTGTT-TGCTTTC--AGAAGGCAAG---CGAtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtatgtgtgtgtgcgtgCGCGCGCGCGAGCACATGTGCATGCATGCGCACTCGTG",
        )
        self.assertEqual(alignment.sequences[1].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[1].seq), 170899992)
        self.assertEqual(
            alignment.sequences[1].seq[155038780 : 155038780 + 199],
            "CCAAAACTCATACACAACTCATCTTCTGGGAAAACAGAATAAGTAAGATGTGTTGCTAGAAATAAGGTAGTGAAGGAAAAAAATAAGAGGATTGATTGACAGGAATGAAAACTCAGACTGGATTTTGAGGCTGTTTCTGGATGGCTTCCACCAAAAACCATAGACCCATCTAAGATCAGGGATACAAGGAGGGCTCATG",
        )
        self.assertEqual(
            alignment[1],
            "CATGAGCCCTCCTTGTATCCCTGATCTTAGATGGGTCTATGGTTTTTGGTGGAAGCCATCCAG-AAACAGCCTCAAAATCCAGTCTGAGTTTTCATTCCTGTCAATCAATCCTCTTATTTTTTTCCTTCACTACC----TTATTTCTAGCAACACATCTTAC-TTATTCTGTTTTCCCAGAAGA-------TGAGTTGTGTATGAGTTTTGG------------------------------------------------------------------",
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[2].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[2].seq), 173908612)
        self.assertEqual(
            alignment.sequences[2].seq[157528852 : 157528852 + 199],
            "CCAAAACTCATACACAACTCATCTTCTGGGAAAACAAAATAAGTAAGATGTGTTGCTAGAAATAAGGTAGTGAAGGAAAAAAATAAGAGGATTGATTGACAGGAATGAAAACTCAGACTGGATTTTGAGGCTGTTTCTGGATGGCTTCCACCAAAAACCATAGACCCATCTAAGATCAGGGATACAAGGAGGGCTCATG",
        )
        self.assertEqual(
            alignment[2],
            "CATGAGCCCTCCTTGTATCCCTGATCTTAGATGGGTCTATGGTTTTTGGTGGAAGCCATCCAG-AAACAGCCTCAAAATCCAGTCTGAGTTTTCATTCCTGTCAATCAATCCTCTTATTTTTTTCCTTCACTACC----TTATTTCTAGCAACACATCTTAC-TTATTTTGTTTTCCCAGAAGA-------TGAGTTGTGTATGAGTTTTGG------------------------------------------------------------------",
        )
        self.assertEqual(
            alignment.sequences[2].annotations["quality"],
            "9999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[2].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[2].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[3].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[3].seq), 174210431)
        self.assertEqual(
            alignment.sequences[3].seq[158049472 : 158049472 + 199],
            "CCAAAACTCATACACAACTCATCTTCTGGGAAAACAAAATAAGTAAGATATGTTGCTAGAAATAAGGTAATGAAGGAAAAAAACAAGAGGATTGATTGACAGGAATGACAACTCAGACTGGATTTTGAGGCTGTTTCTGGATGGCTTCCACCAAAAACCATAGACCCATCTAAGATCAGGGATACAAGGAGGGCTCATG",
        )
        self.assertEqual(
            alignment[3],
            "CATGAGCCCTCCTTGTATCCCTGATCTTAGATGGGTCTATGGTTTTTGGTGGAAGCCATCCAG-AAACAGCCTCAAAATCCAGTCTGAGTTGTCATTCCTGTCAATCAATCCTCTTGTTTTTTTCCTTCATTACC----TTATTTCTAGCAACATATCTTAC-TTATTTTGTTTTCCCAGAAGA-------TGAGTTGTGTATGAGTTTTGG------------------------------------------------------------------",
        )
        self.assertEqual(alignment.sequences[3].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[3].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[3].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[3].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[4].id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(alignment.sequences[4].seq), 359464)
        self.assertEqual(
            alignment.sequences[4].seq[178634 : 178634 + 204],
            "TCAAAACTCATACACAAGTCACCACTTGTCTTCTGGGAAAACAAAACAAAGTAAGCTATGTCACCTGAAATAACAAGGCAACCTAAAAATAAGGGACTGATTGCCAGCAATGACCACAGAGACTGGGTTTTCAGGTGGTTTTCTGGACAGCTTTCCCCAGAAACCAGGGCTACGCCTAAGACAGAGACACGATGAGTGCATGTG",
        )
        self.assertEqual(
            alignment[4],
            "CACATGCACTCATCGTGTCTCTG-TCTTAGGCGTAGCCCTGGTTTCTGGGGAAAGCTGTCCAGAAAACCACCTGAAAACCCAGTCTCTGTGGTCATTGCTGGCAATCAGTCC-CTTATTTTTAGGTTGCCTTG------TTATTTCAGGTGACATAGCTTACTTTGTTTTGTTTTCCCAGAAGACAAGTGGTGACTTGTGTATGAGTTTTGA------------------------------------------------------------------",
        )
        self.assertEqual(
            alignment.sequences[4].annotations["quality"],
            "999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[4].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[4].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[4].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[4].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[5].id, "cavPor2.scaffold_290371")
        self.assertEqual(len(alignment.sequences[5].seq), 39932)
        self.assertEqual(
            alignment.sequences[5].seq[39417 : 39417 + 205],
            "CTCAAACtcacacaaaattcatttatcttccaggaaaacagaacaaaataagatATGTTTTGCAATAACGAGGTCACTATGGGGTAAAAACAAGAGGATTGGCAGGAGGGAATAACAATTTAGACTAGGTTTTGAAGTTCTTTTCTGGATGGCTTTCACCAAAAACCACAGCCACACCTCGGATCAGGATCGCAAGGGATACGCA",
        )
        self.assertEqual(
            alignment[5],
            "--TGCGTATCCCTTGCGATCCTGATCCGAGGTGTGGCTGTGGTTTTTGGTGAAAGCCATCCAGAAAAGAACTTCAAAACCTAGTCTAAATTGTTATTCCCTCCTGCCAATCCTCTTGTTTTTACCCCATAGTGACCTCGTTATTGCAA--AACATatcttattttgttctgttttcctggaagataaa---tgaattttgtgtgaGTTTGAG------------------------------------------------------------------",
        )
        self.assertEqual(
            alignment.sequences[5].annotations["quality"],
            "9999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999799999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[5].annotations["leftStatus"], "N")
        self.assertEqual(alignment.sequences[5].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[5].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[5].annotations["rightCount"], 0)
        self.assertEqual(len(alignment.sequences), 6)
        self.assertNotIn("empty", alignment.annotations)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
            [[  3012566,   3012566,   3012587,   3012587,   3012626,   3012627,
                3012675,   3012676,   3012695,   3012695,   3012695,   3012695,
                3012704,   3012706,   3012718,   3012719,   3012724,   3012724,
                3012731,   3012731,   3012737,   3012741,   3012741,   3012762,
                3012828],
             [155038979, 155038977, 155038956, 155038955, 155038916, 155038916,
              155038868, 155038867, 155038848, 155038847, 155038845, 155038845,
              155038836, 155038834, 155038822, 155038822, 155038817, 155038816,
              155038809, 155038807, 155038801, 155038801, 155038801, 155038780,
              155038780],
             [157529051, 157529049, 157529028, 157529027, 157528988, 157528988,
              157528940, 157528939, 157528920, 157528919, 157528917, 157528917,
              157528908, 157528906, 157528894, 157528894, 157528889, 157528888,
              157528881, 157528879, 157528873, 157528873, 157528873, 157528852,
              157528852],
             [158049671, 158049669, 158049648, 158049647, 158049608, 158049608,
              158049560, 158049559, 158049540, 158049539, 158049537, 158049537,
              158049528, 158049526, 158049514, 158049514, 158049509, 158049508,
              158049501, 158049499, 158049493, 158049493, 158049493, 158049472,
              158049472],
             [   178838,    178836,    178815,    178815,    178776,    178775,
                 178727,    178727,    178708,    178707,    178707,    178707,
                 178698,    178696,    178684,    178683,    178678,    178677,
                 178670,    178668,    178662,    178658,    178655,    178634,
                 178634],
             [    39622,     39622,     39601,     39600,     39561,     39560,
                  39512,     39511,     39492,     39491,     39489,     39485,
                  39476,     39476,     39464,     39463,     39458,     39457,
                  39450,     39448,     39442,     39438,     39438,     39417,
                  39417],
            ]
                    # fmt: on
                ),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 128047)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3012828 : 3012828 + 168],
            "TTTGCATAGACTCTCTTGGCAACAAAATAACGTTATATTTAAACATCCATTAAAATAATGCACTTAGCACAGCCTGCCCTGAGGGATGAACACTATTGTTAAAGAACTATTCCGCTAAGGCAGCAACCTCTGGATCTTCAGCATTCTGGCGCCATCTGCTGGTCATAT",
        )
        self.assertEqual(
            alignment[0],
            "TTTGCATAGACTCTCTTGGCAACAAAATAACGTTATATTTAAACATCCATTAAAATAATGCACTTAGCACAGCCTGCCCTGAGGGAT----GAACACT--ATTGTTAA-AGAACTATTCCGCTAAGGCAGCAACCTCTGGATCTTCAGCATTCTGGCGCCATCTGCTGGTCATAT",
        )
        self.assertEqual(alignment.sequences[1].id, "cavPor2.scaffold_290371")
        self.assertEqual(len(alignment.sequences[1].seq), 39932)
        self.assertEqual(
            alignment.sequences[1].seq[39251 : 39251 + 166],
            "ATATGACCAGCAGATGGCACTAGAATCCCACGGAACCTGAAGTTGCCTCCTGCTAAACAGCTTTGCTGACAATAAAGGTCTCAGAGCTACTAGAGCGAACCATGATAACCTCATTATTTTCataaatgcttaaatataagACTATTTAGCTGTCAAGACATCCCAC",
        )
        self.assertEqual(
            alignment[1],
            "-----GTGGGATGTCTTGACAGCTAAATAGTcttatatttaagcatttatGAAAATAATGAGGTTATCATGGTTCGCTCTAGTAGCTCTGAGACCTTT--ATTGTCAGCAAAGCTGTTTAGC--AGGAGGCAACTTCAGGTTCCGTGGGATTCTAGTGCCATCTGCTGGTCATAT",
        )
        self.assertEqual(
            alignment.sequences[1].annotations["quality"],
            "9999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[2].id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(alignment.sequences[2].seq), 359464)
        self.assertEqual(
            alignment.sequences[2].seq[178473 : 178473 + 161],
            "CTGTGACCAGCAGATGGCGCTAAAATCCCACAGATCCTGAGGCTGCTGCTCCGCAAAACCGCTTCCCTGACAGTAGCCTCTCACTTCTAAGGGGCTGAGCTAAGTTTATCATTTCAACAAATGTTTAAATAGAACACTATTTGCCTGTCACGAATGTCCAC",
        )
        self.assertEqual(
            alignment[2],
            "-----GTGGACATTCGTGACAGGCAAATAGTGTTCTATTTAAACATTTGTTGAAATGATAAACTTAGCTCAGCC--CCTT---AGAAGTG-AGAGGCT--ACTGTCAGGGAAGCGGTTTTGCG-GAGCAGCAGCCTCAGGATCTGTGGGATTTTAGCGCCATCTGCTGGTCACAG",
        )
        self.assertEqual(
            alignment.sequences[2].annotations["quality"],
            "99999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999989996999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[2].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[2].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[3].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[3].seq), 174210431)
        self.assertEqual(
            alignment.sequences[3].seq[158049311 : 158049311 + 161],
            "ATGTGACCAGCAGATGGCATTAGAGTCCCACAGATCCTGAGATTGAGCCTTAGCTAAACAACTTCACTAGGCATACCACATCTGTCAAGGCAGGTTATGCTAAGCTTATTATTTTCATAAATGTTTAAATATAACGCTATTTGGTTGTCAGGATAGTCCAC",
        )
        self.assertEqual(
            alignment[3],
            "-----GTGGACTATCCTGACAACCAAATAGCGTTATATTTAAACATTTATGAAAATAATAAGCTTAGCATAACCTGCCTTGACAGATGTG-GTATGC-------CTAGTGAAGTTGTTTAGCT-AAGGCTCAATCTCAGGATCTGTGGGACTCTAATGCCATCTGCTGGTCACAT",
        )
        self.assertEqual(alignment.sequences[3].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[3].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[3].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[3].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[4].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[4].seq), 173908612)
        self.assertEqual(
            alignment.sequences[4].seq[157528691 : 157528691 + 161],
            "GTGTGACCAGCAGATGGCATTAGAGTCCCACAGATCCTGAGATTGAGCCATAGCTAAGCAACTTCACTAGGCATACCACATCTGTCAAGGCAGGTTATGCTAAGCTTATTATTTTCATAAACGTTTAAATATAACGCTATTTGGTTGTCAAGATAGTCCGC",
        )
        self.assertEqual(
            alignment[4],
            "-----GCGGACTATCTTGACAACCAAATAGCGTTATATTTAAACGTTTATGAAAATAATAAGCTTAGCATAACCTGCCTTGACAGATGTG-GTATGC-------CTAGTGAAGTTGCTTAGCT-ATGGCTCAATCTCAGGATCTGTGGGACTCTAATGCCATCTGCTGGTCACAC",
        )
        self.assertEqual(
            alignment.sequences[4].annotations["quality"],
            "99999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[4].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[4].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[4].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[4].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[5].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[5].seq), 170899992)
        self.assertEqual(
            alignment.sequences[5].seq[155038619 : 155038619 + 161],
            "GTGTGACCAGCAGATGGCATTAGAGTCCCACAGATCCTGAGATTGAGCCATAGCTAAGCAACTTCACTAGGCATACCACATCTGTCAAGGCAGGTTATGCTAAGCTTATTATTTTCACAAACGTTTAAATATAACGCTATTTGGTTGTCAAGATAGTCCAC",
        )
        self.assertEqual(
            alignment[5],
            "-----GTGGACTATCTTGACAACCAAATAGCGTTATATTTAAACGTTTGTGAAAATAATAAGCTTAGCATAACCTGCCTTGACAGATGTG-GTATGC-------CTAGTGAAGTTGCTTAGCT-ATGGCTCAATCTCAGGATCTGTGGGACTCTAATGCCATCTGCTGGTCACAC",
        )
        self.assertEqual(alignment.sequences[5].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[5].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[5].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[5].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[6].id, "echTel1.scaffold_288249")
        self.assertEqual(len(alignment.sequences[6].seq), 100002)
        self.assertEqual(
            alignment.sequences[6].seq[87492 : 87492 + 169],
            "TTGGGGTGGATGCTCTTGGCAGTCACACAGTGCTCTATTTTAGGATTTACTAGAACAATGAGTTTGTCATAACTGCCTCCTCCCAAGTGGGAAGCTGAAGCACCAAGTACTGACTCAGCAGTCCCTGACCTCACATCCATGGAAATCTAGTGACATCTGCTGGACACAT",
        )
        self.assertEqual(
            alignment[6],
            "TTGGGGTGGATGCTCTTGGCAGTCACACAGTGCTCTATTTTAGGATTTACTAGAACAATGAGTTTGTCATAACT-GCCTCCTCCCAAGTG-GGAAGCTGAAGCACCAA-GTACTGACTCAGC--AGTCCCTGACCTCA-CATCCATGGAAATCTAGTGACATCTGCTGGACACAT",
        )
        self.assertEqual(
            alignment.sequences[6].annotations["quality"],
            "9998999999897966589999999999967689989799789997987889999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[6].annotations["leftStatus"], "N")
        self.assertEqual(alignment.sequences[6].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[6].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[6].annotations["rightCount"], 7564)
        self.assertEqual(len(alignment.sequences), 7)
        self.assertNotIn("empty", alignment.annotations)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
            [[  3012828,   3012833,   3012902,   3012903,   3012904,   3012908,
                3012911,   3012915,   3012915,   3012915,   3012921,   3012922,
                3012922,   3012926,   3012930,   3012930,   3012943,   3012944,
                3012945,   3012959,   3012960,   3012996],
             [    39417,     39417,     39348,     39347,     39346,     39342,
                  39339,     39335,     39332,     39331,     39325,     39324,
                  39324,     39320,     39316,     39315,     39302,     39302,
                  39302,     39288,     39287,     39251],
             [   178634,    178634,    178565,    178565,    178565,    178561,
                 178561,    178557,    178554,    178554,    178548,    178547,
                 178547,    178543,    178539,    178538,    178525,    178524,
                 178524,    178510,    178509,    178473],
             [158049472, 158049472, 158049403, 158049402, 158049401, 158049397,
              158049394, 158049390, 158049387, 158049387, 158049381, 158049381,
              158049381, 158049381, 158049377, 158049376, 158049363, 158049362,
              158049362, 158049348, 158049347, 158049311],
             [157528852, 157528852, 157528783, 157528782, 157528781, 157528777,
              157528774, 157528770, 157528767, 157528767, 157528761, 157528761,
              157528761, 157528761, 157528757, 157528756, 157528743, 157528742,
              157528742, 157528728, 157528727, 157528691],
             [155038780, 155038780, 155038711, 155038710, 155038709, 155038705,
              155038702, 155038698, 155038695, 155038695, 155038689, 155038689,
              155038689, 155038689, 155038685, 155038684, 155038671, 155038670,
              155038670, 155038656, 155038655, 155038619],
             [    87492,     87497,     87566,     87566,     87567,     87571,
                  87574,     87578,     87581,     87581,     87587,     87588,
                  87590,     87594,     87598,     87598,     87611,     87611,
                  87611,     87625,     87625,     87661],
            ]
                    # fmt: on
                ),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 98097)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3012996 : 3012996 + 222],
            "AGATGTCTGCTGTGGAGACCTGGCCAACTTTGCTTTCTTCAAAAAGGCAACAGAAGGTAATCAGTTGAATGCCCACCATTAGGAAGGCGACCTCTAGTGCACAAACCTTGACATTTTCCCTTTTAATGGAATTTAACAGAAGTTCAGGATGTTCTTTGGGTAATTTACAATTAGGGGGCAAAAATCAAAAGTATTTCGAGCATATCAAAACTGTTAGCTATG",
        )
        self.assertEqual(
            alignment[0],
            "AGATGTCTGCTGTGGAGA-------CCTGGCCAACTTTG----CTT--TCTTC-----AAAAAGGCAACAGAAGGTAATCAGTTGAATGCCCACCA-----TTAGGAAGGCGACCTCTAGTGCACAAACCTTGAC-ATTTTCCCTTTTAATGGAA-TTTAACAGAAGTTCAGGATGTTCTTTGGGTAATTTACAATT---A----GGGGGCAAAAATCAAAAGTATTTCGAGCATATCAAAACTGTTAGCTATG",
        )
        self.assertEqual(alignment.sequences[1].id, "cavPor2.scaffold_290371")
        self.assertEqual(len(alignment.sequences[1].seq), 39932)
        self.assertEqual(
            alignment.sequences[1].seq[39074 : 39074 + 177],
            "AAATGcctctgcttttccttttctagctgtAAATTGTTTAAATGGAATTCTGAACCTGACAGAGTcgaagaaaaaaattgcaaagctTTTCCGCTAGAGGTCACTCTCCCCACTATGTGGTGTGGTCGGCTTTAGCGCGGAAGACTGGAGTTAGCCAGCCCGTGGCTGCAGACGTCA",
        )
        self.assertEqual(
            alignment[1],
            "TGACGTCTGCAGCCACGG-------GCTGGCTAACTCCA--GTCTT--CCGCG-----CTAAAGCCGAC--------------------CACACCACATAGTGGGGAGAGTGACCTCTAGCGGAAAagctttgca-atttttttcttc----gAC-TCTGTCAG--GTTCAGAATTCCATTTAAACAATTTacagct---a----gaaaaggaaaagcagaggCATTT--------------------------",
        )
        self.assertEqual(
            alignment.sequences[1].annotations["quality"],
            "999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999997999999999995699999999999336991774687",
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "N")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[2].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[2].seq), 170899992)
        self.assertEqual(
            alignment.sequences[2].seq[155038435 : 155038435 + 184],
            "ATACTCTTACCTCCCCCCACCATTACAAATTGTCTAAATGGCATTCTAAGCCTTTATGAATTCAATTAAAAAGAAAAATGCCAAGGCCTGCTCGCTAGAGGTCGCTCTCCTAACTATGGGGTCCTGTTAGTCTTATAGAGTAAGACAGCAGAGTAGGCCAAGAAAGTGATGGCAGCAGCCATCA",
        )
        self.assertEqual(
            alignment[2],
            "TGATGGCTGCTGCCATCA---CTTTCTTGGCCTACTCTGCTGTCTT--ACTCT-----ATAAGACTAAC--------------------AGGACCCCATAGTTAGGAGAGCGACCTCTAGCGAGCAGGCCTTGGC-ATTTTTCTTTTTAATTGAA-TTCATAAA-GGCTTAGAATGCCATTTAGACAATTTGTAATG---GTGGGGGGAGGTAAGAGTAT----------------------------------",
        )
        self.assertEqual(alignment.sequences[2].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[2].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[2].annotations["rightCount"], 14)
        self.assertEqual(alignment.sequences[3].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[3].seq), 173908612)
        self.assertEqual(
            alignment.sequences[3].seq[157528507 : 157528507 + 184],
            "ATACTCTTACCTCCCCCCACCATTACAAATTGTCTAAATGGCATTCTAAGCCTTTATGAATTTAATTAAAAAGAAAAATTCCAAGGCCTGCTCGCTAGAGGTCGCTCTCCTAACCATGGGGTCCTGTTAGTCTTATAGAATAAGACAGCAGAGTAGGCCAAGAAAGTGATGGCAGCAGCCATCA",
        )
        self.assertEqual(
            alignment[3],
            "TGATGGCTGCTGCCATCA---CTTTCTTGGCCTACTCTGCTGTCTT--ATTCT-----ATAAGACTAAC--------------------AGGACCCCATGGTTAGGAGAGCGACCTCTAGCGAGCAGGCCTTGGA-ATTTTTCTTTTTAATTAAA-TTCATAAA-GGCTTAGAATGCCATTTAGACAATTTGTAATG---GTGGGGGGAGGTAAGAGTAT----------------------------------",
        )
        self.assertEqual(
            alignment.sequences[3].annotations["quality"],
            "9999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[3].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[3].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[3].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[3].annotations["rightCount"], 14)
        self.assertEqual(alignment.sequences[4].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[4].seq), 174210431)
        self.assertEqual(
            alignment.sequences[4].seq[158049094 : 158049094 + 217],
            "CATCATTGACAGGTTTGAATACATTCAATATTTATACTCTTACCTCCCCCCACCATTACAAATTGTCTAAATGGCATTCTAAGCCTTTATGAATTCAATTAAAAAGAAAAATTCCAAGGCCTGCTCGCTAGAGGTCGCTCTCCTAACTATGGGGTCCTGTTAGTCTTATAGAGTAAGACGGCCTAGTAGGCCAAGAAAGTGATGGCAGCAGCCATCA",
        )
        self.assertEqual(
            alignment[4],
            "TGATGGCTGCTGCCATCA---CTTTCTTGGCCTACTAGGCCGTCTT--ACTCT-----ATAAGACTAAC--------------------AGGACCCCATAGTTAGGAGAGCGACCTCTAGCGAGCAGGCCTTGGA-ATTTTTCTTTTTAATTGAA-TTCATAAA-GGCTTAGAATGCCATTTAGACAATTTGTAATG---GTGGGGGGAGGTAAGAGTATAAATATTGAATGTAT-TCAAACCTGTCAATGATG",
        )
        self.assertEqual(alignment.sequences[4].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[4].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[4].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[4].annotations["rightCount"], 2)
        self.assertEqual(alignment.sequences[5].id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(alignment.sequences[5].seq), 359464)
        self.assertEqual(
            alignment.sequences[5].seq[178247 : 178247 + 226],
            "TAACAAGGCTGAATGCCAGCAATATTTGCACTCTCAGGTCCCCACACATACACTATGAATTGTCTAACTTGCATTTTGAGCCGTTGTAAAATTCAATTAAAAAGAAAAAAATCCAAGGCTTGTTCACCAGAGGTCGCTCTAGTAACTGCAGGGTCCGATTTTTCTTTTGTTGCAGAATTTAAGACGGTGGAGTTGGCAGGGAAAGCGGTGGTGGCGGCGGCCGTCC",
        )
        self.assertEqual(
            alignment[5],
            "GGACGGCCGCCGCCACCACCGCTTTCCCTGCCAACTCCACCGTCTTAAATTCTGCAACAAAAGAAAAAT--------------------CGGACCCTGCAGTTACTAGAGCGACCTCTGGTGAACAAGCCTTGGATTTTTTTCTTTTTAATTGAATTTTACAAC-GGCTCAAAATGCAAGTTAGACAATTCATAGTGTATGTGTGGGGACCTGAGAGTGCAAATATTGCTGGCAT-TCAGCCTTGTTA------",
        )
        self.assertEqual(
            alignment.sequences[5].annotations["quality"],
            "9999989999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999699999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[5].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[5].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[5].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[5].annotations["rightCount"], 2931)
        self.assertEqual(alignment.sequences[6].id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(alignment.sequences[6].seq), 498454)
        self.assertEqual(
            alignment.sequences[6].seq[331078 : 331078 + 159],
            "TAATTATAAGTTGTCTAAATGGCATTTGAACCTTTATAAATTCTGTTTTAAAGAAAAATTCCAAGGTTTGCTCACTAGAGGTCTCCTCTCCCGGAGGATCCAGTTTGTCTTTCGTTGAGAGAGCTGAGTGGTCCAAGAAAGTGACAGCTGTAGCCATCA",
        )
        self.assertEqual(
            alignment[6],
            "TGATGGCTACAGCTGTCA---CTTTCTTGGACCACTCAGCTCTCTC--AAC-------GAAAGACAAAC--------------------TGGATCCTCCGG---GAGAGGAGACCTCTAGTGAGCAAACCTTGGA-ATTTTTCTTTAAAACAGAA-TTTATAAA-GGTTCA-AATGCCATTTAGACAACTTATAATT---A-----------------------------------------------------",
        )
        self.assertEqual(
            alignment.sequences[6].annotations["quality"],
            "234233433332122232158211222021213444332433213323732111754365002326236241111233524253535324593652222413766453735782535545832457354545484445655854554657999679999",
        )
        self.assertEqual(alignment.sequences[6].annotations["leftStatus"], "N")
        self.assertEqual(alignment.sequences[6].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[6].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[6].annotations["rightCount"], 4145)
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "echTel1.scaffold_288249")
        self.assertEqual(len(record.seq), 100002)
        self.assertEqual(segment, (87661, 95225))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 7)
        self.assertEqual(len(alignment.annotations["empty"]), 1)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
            [[  3012996,   3013014,   3013014,   3013014,   3013028,   3013028,
                3013028,   3013031,   3013031,   3013034,   3013036,   3013036,
                3013047,   3013067,   3013074,   3013074,   3013077,   3013108,
                3013108,   3013120,   3013124,   3013127,   3013127,   3013135,
                3013136,   3013137,   3013142,   3013143,   3013168,   3013168,
                3013169,   3013169,   3013184,   3013192,   3013199,   3013200,
                3013212,   3013218],
             [    39251,     39233,     39233,     39233,     39219,     39219,
                  39217,     39214,     39214,     39211,     39209,     39209,
                  39198,     39198,     39191,     39186,     39183,     39152,
                  39152,     39140,     39140,     39137,     39137,     39129,
                  39129,     39129,     39124,     39123,     39098,     39098,
                  39097,     39097,     39082,     39074,     39074,     39074,
                  39074,     39074],
             [155038619, 155038601, 155038601, 155038597, 155038583, 155038581,
              155038579, 155038576, 155038576, 155038573, 155038571, 155038571,
              155038560, 155038560, 155038553, 155038548, 155038545, 155038514,
              155038514, 155038502, 155038498, 155038495, 155038495, 155038487,
              155038487, 155038486, 155038481, 155038480, 155038455, 155038455,
              155038454, 155038450, 155038435, 155038435, 155038435, 155038435,
              155038435, 155038435],
             [157528691, 157528673, 157528673, 157528669, 157528655, 157528653,
              157528651, 157528648, 157528648, 157528645, 157528643, 157528643,
              157528632, 157528632, 157528625, 157528620, 157528617, 157528586,
              157528586, 157528574, 157528570, 157528567, 157528567, 157528559,
              157528559, 157528558, 157528553, 157528552, 157528527, 157528527,
              157528526, 157528522, 157528507, 157528507, 157528507, 157528507,
              157528507, 157528507],
             [158049311, 158049293, 158049293, 158049289, 158049275, 158049273,
              158049271, 158049268, 158049268, 158049265, 158049263, 158049263,
              158049252, 158049252, 158049245, 158049240, 158049237, 158049206,
              158049206, 158049194, 158049190, 158049187, 158049187, 158049179,
              158049179, 158049178, 158049173, 158049172, 158049147, 158049147,
              158049146, 158049142, 158049127, 158049119, 158049112, 158049112,
              158049100, 158049094],
             [   178473,    178455,    178452,    178448,    178434,    178432,
                 178430,    178427,    178425,    178422,    178420,    178415,
                 178404,    178404,    178397,    178392,    178389,    178358,
                 178357,    178345,    178341,    178338,    178337,    178329,
                 178329,    178328,    178323,    178322,    178297,    178294,
                 178293,    178289,    178274,    178266,    178259,    178259,
                 178247,    178247],
             [   331237,    331219,    331219,    331215,    331201,    331199,
                 331197,    331194,    331194,    331191,    331191,    331191,
                 331180,    331180,    331173,    331168,    331168,    331137,
                 331137,    331125,    331121,    331118,    331118,    331110,
                 331110,    331109,    331104,    331104,    331079,    331079,
                 331078,    331078,    331078,    331078,    331078,    331078,
                 331078,    331078],
            ]
                    # fmt: on
                ),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3013218 : 3013218 + 219],
            "agccaggcgtggtggcacacacctttactcccagcatttggggggcagaggcaggtggatctgtgagtttgaggccagcctggtctacagagggagtctcaggacagccagagctacacagaaataacctgcctagaaaaacaaaacaaaacaaaacatcaaaactcaaaacaaaTAAAAAAAATAAAAAACCCAACCTAAACCAAATAACAAAACACT",
        )
        self.assertEqual(
            alignment[0],
            "agccaggcgtggtggcacacacctttactcccagcatttggggggcagaggcaggtggatctgtgagtttgaggccagcctggtctacagagggagtctcaggacagccagagctacacagaaataacctgcctagaaaaacaaaacaaaacaaaacatcaaaactcaaaacaaaTAAAAAAAATAAAAAACCCAACCTAAACCAAATAACAAAACACT",
        )
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(record.seq), 498454)
        self.assertEqual(segment, (331078, 326933))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "echTel1.scaffold_288249")
        self.assertEqual(len(record.seq), 100002)
        self.assertEqual(segment, (87661, 95225))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][2]
        (record, segment, status) = empty
        self.assertEqual(record.id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(record.seq), 359464)
        self.assertEqual(segment, (178247, 175316))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][3]
        (record, segment, status) = empty
        self.assertEqual(record.id, "hg18.chr6")
        self.assertEqual(len(record.seq), 170899992)
        self.assertEqual(segment, (155038435, 155038421))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][4]
        (record, segment, status) = empty
        self.assertEqual(record.id, "panTro2.chr6")
        self.assertEqual(len(record.seq), 173908612)
        self.assertEqual(segment, (157528507, 157528493))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][5]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ponAbe2.chr6")
        self.assertEqual(len(record.seq), 174210431)
        self.assertEqual(segment, (158049094, 158049092))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 1)
        self.assertEqual(len(alignment.annotations["empty"]), 6)
        self.assertTrue(
            numpy.array_equal(alignment.coordinates, numpy.array([[3013218, 3013437]]))
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 40604)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3013437 : 3013437 + 166],
            "TCCAAAATGGTTAGCTATGCCCAACTCCTTTCACTCCAAGAAAATATCCTAACCATGTAAGAGAGCTAGCCTGTTGGTGGCAGCCAAGCCTGATGGTGGCAGACTAGATTGATGGTGCCAGACTACTTTATGGCTGTATCATTTTCCATTCATGTGTTGTGTTATA",
        )
        self.assertEqual(
            alignment[0],
            "TCCAAAATGGTTAGCTATGCCCAACTCCTTTCACTCCAAGAAAATATCCTAACCATGTAAGAGAGCTAGCCTGTTGGTGGCAGCCAAGCCTGATGGTGGCAGACTAGATTGATGGTGCCAGACTACTTTATGGCTGTATCATTTTCCATTCATGTGTTGTGTTATA",
        )
        self.assertEqual(alignment.sequences[1].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[1].seq), 173908612)
        self.assertEqual(
            alignment.sequences[1].seq[157528363 : 157528363 + 130],
            "TATCAAACAGCCCATGAATGCAAAACGATACACCCATAAAATAGGATGTAATCGCCAAAATCTAATTCTATTGCATCTAGAAGATCTTTATCCACAGTGAAAAGAGTTGAACATCATTGACAGGTTTGAA",
        )
        self.assertEqual(
            alignment[1],
            "TTCAAACCTGTCAATGATGTTCAACTCTTTTCACTGTGGATAAAGATCTTCTAGATGCAATAGAATTAGATT------------------------------------TTGGCGATTACATCCTATTTTATGGGTGTATCGTTTTGCATTCATGGGCTGTTTGATA",
        )
        self.assertEqual(
            alignment.sequences[1].annotations["quality"],
            "9999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999099999999999999999999",
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 14)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 9106)
        self.assertEqual(alignment.sequences[2].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[2].seq), 170899992)
        self.assertEqual(
            alignment.sequences[2].seq[155038291 : 155038291 + 130],
            "TATCAAACAGCCCATGAATGCAAAACGATACACCCATAAAATAGGATGTAATCGCCAAAATCTAATTCTATTGCATCTAGAAGATCTTTATCCACAGTGAAAAGAGTTGAACATCATTGACAGGTTTGAA",
        )
        self.assertEqual(
            alignment[2],
            "TTCAAACCTGTCAATGATGTTCAACTCTTTTCACTGTGGATAAAGATCTTCTAGATGCAATAGAATTAGATT------------------------------------TTGGCGATTACATCCTATTTTATGGGTGTATCGTTTTGCATTCATGGGCTGTTTGATA",
        )
        self.assertEqual(alignment.sequences[2].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[2].annotations["leftCount"], 14)
        self.assertEqual(alignment.sequences[2].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[2].annotations["rightCount"], 9085)
        self.assertEqual(alignment.sequences[3].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[3].seq), 174210431)
        self.assertEqual(
            alignment.sequences[3].seq[158048983 : 158048983 + 109],
            "TATCAAACAGCCCATGAATGCAAAACGATACACCCATAAAATAGGATGTAATCGCCAAAATCTAATTCTATTGCATCTAGAAGATCTTTATCCACAGTGAAAAGAGTTG",
        )
        self.assertEqual(
            alignment[3],
            "---------------------CAACTCTTTTCACTGTGGATAAAGATCTTCTAGATGCAATAGAATTAGATT------------------------------------TTGGCGATTACATCCTATTTTATGGGTGTATCGTTTTGCATTCATGGGCTGTTTGATA",
        )
        self.assertEqual(alignment.sequences[3].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[3].annotations["leftCount"], 2)
        self.assertEqual(alignment.sequences[3].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[3].annotations["rightCount"], 8044)
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(record.seq), 498454)
        self.assertEqual(segment, (331078, 326933))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "echTel1.scaffold_288249")
        self.assertEqual(len(record.seq), 100002)
        self.assertEqual(segment, (87661, 95225))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][2]
        (record, segment, status) = empty
        self.assertEqual(record.id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(record.seq), 359464)
        self.assertEqual(segment, (178247, 175316))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 4)
        self.assertEqual(len(alignment.annotations["empty"]), 3)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
            [[  3013437,   3013458,   3013509,   3013545,   3013603],
             [157528493, 157528472, 157528421, 157528421, 157528363],
             [155038421, 155038400, 155038349, 155038349, 155038291],
             [158049092, 158049092, 158049041, 158049041, 158048983],
            ]
                    # fmt: on
                ),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3013603 : 3013603 + 1041],
            "CCTTCTTAAGAACACTAGACTCAggactggggagatggctcagcagttaagaatcggtgctgttaagagtgggagacaagtttggttcccagctcccacattggtcagctcacagccacccgtaactctaagatggtacacacctttaatcccaggagacagaggcaatcagatctgagttcaagattcagcctgagacagagcatgttccaaattcaggcatggtgggtcatacctttaatatgggacataccttctgctggaggcctacctaaggacaacggagaaaggaagtattcgttcttctcctgcttgcacttacttgccagcgcatctactggaacccacttcttcaggattccagcttatacaggagaccagctgaaatatccagcctctcgggactgaacaagtactagagtctcagacttcccattcacagctgcccattgttggttggttgtactacagactgtaagtcattgtaataatttcccttaatatatagagacattatataagttctgtgactctagagaaccctgactagtacaCGTGGCTAACTAGAAAGctctggtatgtgcttacttaatgctgaggttttaggcatggccacggtgctctgcttcttatgtgggtgctgggaatgcagactcaggtcctcatgtgtatgcagcaaacacttcatacactcagctgcttccctaacccTATGCTTGTGTCTTATTACTAACTTGTGAAAAGCTTTGAGTTTATTTTCTATGTTTTCAACCACTTTCTTGAGTATGCTCAGCTCGTGGCTTTAAACTGGATTTCCCCCTAATATGTAATGACTATAAGTATTCCTTAAATAGGACACACTTTTGTTATACTTTTTGTTATCatataaaatatttcaaaaaaatttttttGCTATTTTTATCTTTGAGCCATTGGTCATTTTGACGTGTATCTCTTGATTTTTATAGATGGTAATATTTTATGTATTGCTAGCCAATCTCGTTTTCTTGTTTGCTTGCTTGTTTGTTTTGGTCAATGCAG",
        )
        self.assertEqual(
            alignment[0],
            "CCTTCTTAAGAACACTAGACTCAggactggggagatggctcagcagttaagaatcggtgctgttaagagtgggagacaagtttggttcccagctcccacattggtcagctcacagccacccgtaactctaagatggtacacacctttaatcccaggagacagaggcaatcagatctgagttcaagattcagcctgagacagagcatgttccaaattcaggcatggtgggtcatacctttaatatgggacataccttctgctggaggcctacctaaggacaacggagaaaggaagtattcgttcttctcctgcttgcacttacttgccagcgcatctactggaacccacttcttcaggattccagcttatacaggagaccagctgaaatatccagcctctcgggactgaacaagtactagagtctcagacttcccattcacagctgcccattgttggttggttgtactacagactgtaagtcattgtaataatttcccttaatatatagagacattatataagttctgtgactctagagaaccctgactagtacaCGTGGCTAACTAGAAAGctctggtatgtgcttacttaatgctgaggttttaggcatggccacggtgctctgcttcttatgtgggtgctgggaatgcagactcaggtcctcatgtgtatgcagcaaacacttcatacactcagctgcttccctaacccTATGCTTGTGTCTTATTACTAACTTGTGAAAAGCTTTGAGTTTATTTTCTATGTTTTCAACCACTTTCTTGAGTATGCTCAGCTCGTGGCTTTAAACTGGATTTCCCCCTAATATGTAATGACTATAAGTATTCCTTAAATAGGACACACTTTTGTTATACTTTTTGTTATCatataaaatatttcaaaaaaatttttttGCTATTTTTATCTTTGAGCCATTGGTCATTTTGACGTGTATCTCTTGATTTTTATAGATGGTAATATTTTATGTATTGCTAGCCAATCTCGTTTTCTTGTTTGCTTGCTTGTTTGTTTTGGTCAATGCAG",
        )
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(record.seq), 498454)
        self.assertEqual(segment, (331078, 326933))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "echTel1.scaffold_288249")
        self.assertEqual(len(record.seq), 100002)
        self.assertEqual(segment, (87661, 95225))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][2]
        (record, segment, status) = empty
        self.assertEqual(record.id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(record.seq), 359464)
        self.assertEqual(segment, (178247, 175316))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][3]
        (record, segment, status) = empty
        self.assertEqual(record.id, "hg18.chr6")
        self.assertEqual(len(record.seq), 170899992)
        self.assertEqual(segment, (155038291, 155029206))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][4]
        (record, segment, status) = empty
        self.assertEqual(record.id, "panTro2.chr6")
        self.assertEqual(len(record.seq), 173908612)
        self.assertEqual(segment, (157528363, 157519257))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][5]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ponAbe2.chr6")
        self.assertEqual(len(record.seq), 174210431)
        self.assertEqual(segment, (158048983, 158040939))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 1)
        self.assertEqual(len(alignment.annotations["empty"]), 6)
        self.assertTrue(
            numpy.array_equal(alignment.coordinates, numpy.array([[3013603, 3014644]]))
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 19159)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3014644 : 3014644 + 45],
            "CCTGTACCCTTTGGTGAGAATTTTTGTTTCAGTGTTAAAAGTTTG",
        )
        self.assertEqual(
            alignment[0], "CCTGTACC---CTTTGGTGAGAATTTTTGTTTCAGTGTTAAAAGTTTG"
        )
        self.assertEqual(alignment.sequences[1].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[1].seq), 170899992)
        self.assertEqual(
            alignment.sequences[1].seq[155029160 : 155029160 + 46],
            "AAAAGTTTAGGATTAAAACAAAATTCTCATAAAAGAAAGGTATAGG",
        )
        self.assertEqual(
            alignment[1], "CCTATACCTTTCTTTTATGAGAA-TTTTGTTTTAATCCTAAAC-TTTT"
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 9085)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[2].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[2].seq), 173908612)
        self.assertEqual(
            alignment.sequences[2].seq[157519211 : 157519211 + 46],
            "AAAAGTTTAGGATTAAAACAAAATTCTCATAAAAGAAAGGTATAGG",
        )
        self.assertEqual(
            alignment[2], "CCTATACCTTTCTTTTATGAGAA-TTTTGTTTTAATCCTAAAC-TTTT"
        )
        self.assertEqual(
            alignment.sequences[2].annotations["quality"],
            "9999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[2].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[2].annotations["leftCount"], 9106)
        self.assertEqual(alignment.sequences[2].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[3].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[3].seq), 133105)
        self.assertEqual(
            alignment.sequences[3].seq[6182 : 6182 + 46],
            "CCTATACCTTTCTTTCATGAGAATTTTGTTTGAATCCTAAACTTTT",
        )
        self.assertEqual(
            alignment[3], "CCTATACCTTTCTTTCATGAGAA-TTTTGTTTGAATCCTAAAC-TTTT"
        )
        self.assertEqual(alignment.sequences[3].annotations["leftStatus"], "N")
        self.assertEqual(alignment.sequences[3].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[3].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[3].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[4].id, "loxAfr1.scaffold_75566")
        self.assertEqual(len(alignment.sequences[4].seq), 10574)
        self.assertEqual(
            alignment.sequences[4].seq[9373 : 9373 + 34],
            "GGAAGTTTTGAATTAAAGCATAATTCTAACCAAA",
        )
        self.assertEqual(
            alignment[4], "------------TTTGGTTAGAA-TTATGCTTTAATTCAAAAC-TTCC"
        )
        self.assertEqual(
            alignment.sequences[4].annotations["quality"],
            "9999969989999999999999998699989997",
        )
        self.assertEqual(alignment.sequences[4].annotations["leftStatus"], "N")
        self.assertEqual(alignment.sequences[4].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[4].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[4].annotations["rightCount"], 0)
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(record.seq), 498454)
        self.assertEqual(segment, (331078, 326933))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "echTel1.scaffold_288249")
        self.assertEqual(len(record.seq), 100002)
        self.assertEqual(segment, (87661, 95225))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][2]
        (record, segment, status) = empty
        self.assertEqual(record.id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(record.seq), 359464)
        self.assertEqual(segment, (178247, 175316))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][3]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ponAbe2.chr6")
        self.assertEqual(len(record.seq), 174210431)
        self.assertEqual(segment, (158048983, 158040939))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 5)
        self.assertEqual(len(alignment.annotations["empty"]), 4)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
            [[  3014644,   3014652,   3014652,   3014653,   3014664,   3014665,
                3014684,   3014685,   3014689],
             [155029206, 155029198, 155029195, 155029194, 155029183, 155029183,
              155029164, 155029164, 155029160],
             [157519257, 157519249, 157519246, 157519245, 157519234, 157519234,
              157519215, 157519215, 157519211],
             [     6182,      6190,      6193,      6194,      6205,      6205,
                   6224,      6224,      6228],
             [     9407,      9407,      9407,      9407,      9396,      9396,
                   9377,      9377,      9373],
            ]
                    # fmt: on
                ),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 40840)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3014689 : 3014689 + 53],
            "GGGAGCATAAAACTCTAAATCTGCTAAATGTCTTGTCCCTTTGGAAAGAGTTG",
        )
        self.assertEqual(
            alignment[0], "GGGAGCATAAAACTCTAAATCTGCTAAATGTCTTGTCCCT-TTGGAAAGAGTTG"
        )
        self.assertEqual(alignment.sequences[1].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[1].seq), 170899992)
        self.assertEqual(
            alignment.sequences[1].seq[155029107 : 155029107 + 53],
            "CCACTATTTCCCAAAAGATTAGATATTTCACAGATTAAATGGTTTATGATCCC",
        )
        self.assertEqual(
            alignment[1], "GGGATCATAAACCATTTAATCTGTGAAATATCTAATCTTT-TGGGAAATAGTGG"
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 401)
        self.assertEqual(alignment.sequences[2].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[2].seq), 173908612)
        self.assertEqual(
            alignment.sequences[2].seq[157519158 : 157519158 + 53],
            "CCACTATTTCCCAAAAGATTAGATATTTCACAGATTAAATGGTTTATGATCCC",
        )
        self.assertEqual(
            alignment[2], "GGGATCATAAACCATTTAATCTGTGAAATATCTAATCTTT-TGGGAAATAGTGG"
        )
        self.assertEqual(
            alignment.sequences[2].annotations["quality"],
            "99999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[2].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[2].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[2].annotations["rightCount"], 400)
        self.assertEqual(alignment.sequences[3].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[3].seq), 133105)
        self.assertEqual(
            alignment.sequences[3].seq[6228 : 6228 + 53],
            "GGGATCATAAGCCATTTAATCTGTGAAATGTGAAATCTTTTGGGAAACAGTGG",
        )
        self.assertEqual(
            alignment[3], "GGGATCATAAGCCATTTAATCTGTGAAATGTGAAATCTTT-TGGGAAACAGTGG"
        )
        self.assertEqual(alignment.sequences[3].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[3].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[3].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[3].annotations["rightCount"], 2)
        self.assertEqual(alignment.sequences[4].id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(alignment.sequences[4].seq), 359464)
        self.assertEqual(
            alignment.sequences[4].seq[175264 : 175264 + 52],
            "CAGCTATTGCCCAAGTGATTTGATATTTCATAGATTAAAAGTTTATGCTTCC",
        )
        self.assertEqual(
            alignment[4], "GGAAGCATAAACT-TTTAATCTATGAAATATCAAATCACT-TGGGCAATAGCTG"
        )
        self.assertEqual(
            alignment.sequences[4].annotations["quality"],
            "7455455669566996656997698955556899975999984787795599",
        )
        self.assertEqual(alignment.sequences[4].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[4].annotations["leftCount"], 2931)
        self.assertEqual(alignment.sequences[4].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[4].annotations["rightCount"], 2)
        self.assertEqual(alignment.sequences[5].id, "loxAfr1.scaffold_75566")
        self.assertEqual(len(alignment.sequences[5].seq), 10574)
        self.assertEqual(
            alignment.sequences[5].seq[9319 : 9319 + 54],
            "CAGCTTTTTCCCCTGAAGATTTGGCATTTCGCAGACTAAATGGTTTATACTCCC",
        )
        self.assertEqual(
            alignment[5], "GGGAGTATAAACCATTTAGTCTGCGAAATGCCAAATCTTCAGGGGAAAAAGCTG"
        )
        self.assertEqual(
            alignment.sequences[5].annotations["quality"],
            "899989799999979999999999999999797999999999999999999999",
        )
        self.assertEqual(alignment.sequences[5].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[5].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[5].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[5].annotations["rightCount"], 2)
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(record.seq), 498454)
        self.assertEqual(segment, (331078, 326933))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "echTel1.scaffold_288249")
        self.assertEqual(len(record.seq), 100002)
        self.assertEqual(segment, (87661, 95225))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][2]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ponAbe2.chr6")
        self.assertEqual(len(record.seq), 174210431)
        self.assertEqual(segment, (158048983, 158040939))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 6)
        self.assertEqual(len(alignment.annotations["empty"]), 3)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
            [[  3014689,   3014702,   3014703,   3014729,   3014729,   3014742],
             [155029160, 155029147, 155029146, 155029120, 155029120, 155029107],
             [157519211, 157519198, 157519197, 157519171, 157519171, 157519158],
             [     6228,      6241,      6242,      6268,      6268,      6281],
             [   175316,    175303,    175303,    175277,    175277,    175264],
             [     9373,      9360,      9359,      9333,      9332,      9319],
            ]
                    # fmt: on
                ),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 411)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3014742 : 3014742 + 36],
            "AAGTTCCCTCCATAATTCCTTCCTCCCACCCCCACA",
        )
        self.assertEqual(alignment[0], "AAGTTCCCTCCATAATTCCTTCCTCCCACCCCCACA")
        self.assertEqual(alignment.sequences[1].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[1].seq), 133105)
        self.assertEqual(
            alignment.sequences[1].seq[6283 : 6283 + 28], "AAATGTATGATCTCCCCATCCTGCCCTG"
        )
        self.assertEqual(alignment[1], "AAATGTA-----TGATCTCCCCATCCTGCCCTG---")
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 2)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 54)
        self.assertEqual(alignment.sequences[2].id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(alignment.sequences[2].seq), 359464)
        self.assertEqual(
            alignment.sequences[2].seq[175231 : 175231 + 31],
            "TGCACGGAGGGGGTGAGGGCATCAGAAATCT",
        )
        self.assertEqual(alignment[2], "AGATTTC-----TGATGCCCTCACCCCCTCCGTGCA")
        self.assertEqual(
            alignment.sequences[2].annotations["quality"],
            "9996999965974999999999999999999",
        )
        self.assertEqual(alignment.sequences[2].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[2].annotations["leftCount"], 2)
        self.assertEqual(alignment.sequences[2].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[2].annotations["rightCount"], 24)
        self.assertEqual(alignment.sequences[3].id, "loxAfr1.scaffold_75566")
        self.assertEqual(len(alignment.sequences[3].seq), 10574)
        self.assertEqual(
            alignment.sequences[3].seq[9290 : 9290 + 27], "TGTGGGGGTGGGGGGTGGCATAAGCCT"
        )
        self.assertEqual(alignment[3], "AGGCTTA-----TG----CCACCCCCCACCCCCACA")
        self.assertEqual(
            alignment.sequences[3].annotations["quality"], "999999999997999999999999999"
        )
        self.assertEqual(alignment.sequences[3].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[3].annotations["leftCount"], 2)
        self.assertEqual(alignment.sequences[3].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[3].annotations["rightCount"], 25)
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(record.seq), 498454)
        self.assertEqual(segment, (331078, 326933))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "echTel1.scaffold_288249")
        self.assertEqual(len(record.seq), 100002)
        self.assertEqual(segment, (87661, 95225))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][2]
        (record, segment, status) = empty
        self.assertEqual(record.id, "hg18.chr6")
        self.assertEqual(len(record.seq), 170899992)
        self.assertEqual(segment, (155029107, 155028706))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][3]
        (record, segment, status) = empty
        self.assertEqual(record.id, "panTro2.chr6")
        self.assertEqual(len(record.seq), 173908612)
        self.assertEqual(segment, (157519158, 157518758))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][4]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ponAbe2.chr6")
        self.assertEqual(len(record.seq), 174210431)
        self.assertEqual(segment, (158048983, 158040939))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 4)
        self.assertEqual(len(alignment.annotations["empty"]), 5)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
            [[3014742, 3014749, 3014754, 3014756, 3014760, 3014775, 3014778],
             [   6283,    6290,    6290,    6292,    6296,    6311,    6311],
             [ 175262,  175255,  175255,  175253,  175249,  175234,  175231],
             [   9317,    9310,    9310,    9308,    9308,    9293,    9290],
            ]
                    # fmt: on
                ),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3014778 : 3014778 + 17], "TCCCATGTCCACCCTGA"
        )
        self.assertEqual(alignment[0], "TCCCATGTCCACCCTGA")
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "loxAfr1.scaffold_75566")
        self.assertEqual(len(record.seq), 10574)
        self.assertEqual(segment, (9290, 9265))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "calJac1.Contig6394")
        self.assertEqual(len(record.seq), 133105)
        self.assertEqual(segment, (6311, 6365))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][2]
        (record, segment, status) = empty
        self.assertEqual(record.id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(record.seq), 498454)
        self.assertEqual(segment, (331078, 326933))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][3]
        (record, segment, status) = empty
        self.assertEqual(record.id, "echTel1.scaffold_288249")
        self.assertEqual(len(record.seq), 100002)
        self.assertEqual(segment, (87661, 95225))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][4]
        (record, segment, status) = empty
        self.assertEqual(record.id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(record.seq), 359464)
        self.assertEqual(segment, (175231, 175207))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][5]
        (record, segment, status) = empty
        self.assertEqual(record.id, "hg18.chr6")
        self.assertEqual(len(record.seq), 170899992)
        self.assertEqual(segment, (155029107, 155028706))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][6]
        (record, segment, status) = empty
        self.assertEqual(record.id, "panTro2.chr6")
        self.assertEqual(len(record.seq), 173908612)
        self.assertEqual(segment, (157519158, 157518758))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][7]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ponAbe2.chr6")
        self.assertEqual(len(record.seq), 174210431)
        self.assertEqual(segment, (158048983, 158040939))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 1)
        self.assertEqual(len(alignment.annotations["empty"]), 8)
        self.assertTrue(
            numpy.array_equal(alignment.coordinates, numpy.array([[3014778, 3014795]]))
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, -12243)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3014795 : 3014795 + 47],
            "GTTTCAGGGGCAGCTCGCTGTTAGCAGCTAAGGCATGGTGTCTCTCA",
        )
        self.assertEqual(
            alignment[0],
            "GTTTCAGGGGCAGCTCGCTG----------------TTAGCAG-CTAAGGCATGGTGTCTCTCA",
        )
        self.assertEqual(alignment.sequences[1].id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(alignment.sequences[1].seq), 359464)
        self.assertEqual(
            alignment.sequences[1].seq[175147 : 175147 + 60],
            "AGAAATACACCCCACCTAAACCATCTAAACCAAATATCCATGACTGCAACTTGTTCCCGA",
        )
        self.assertEqual(
            alignment[1],
            "---TCGGGAACAAGTTGCAGTCATGGATAT-TTGGTTTAGATGGTTTAGGTGGGGTGTATTTCT",
        )
        self.assertEqual(
            alignment.sequences[1].annotations["quality"],
            "899999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 24)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[2].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[2].seq), 133105)
        self.assertEqual(
            alignment.sequences[2].seq[6365 : 6365 + 39],
            "GCCATGAATATTTTAGACATGCAGGTGTGGCGTGTTTCT",
        )
        self.assertEqual(
            alignment[2],
            "-------------------GCCATGAATAT-----TTTAGAC-ATGCAGGTGTGGCGTGTTTCT",
        )
        self.assertEqual(alignment.sequences[2].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[2].annotations["leftCount"], 54)
        self.assertEqual(alignment.sequences[2].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[3].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[3].seq), 170899992)
        self.assertEqual(
            alignment.sequences[3].seq[155028686 : 155028686 + 20],
            "AGAAACATGCCACACCTGCt",
        )
        self.assertEqual(
            alignment[3],
            "--------------------------------------------aGCAGGTGTGGCATGTTTCT",
        )
        self.assertEqual(alignment.sequences[3].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[3].annotations["leftCount"], 401)
        self.assertEqual(alignment.sequences[3].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[3].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[4].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[4].seq), 173908612)
        self.assertEqual(
            alignment.sequences[4].seq[157518738 : 157518738 + 20],
            "AGAAACATGCCACACCTGCt",
        )
        self.assertEqual(
            alignment[4],
            "--------------------------------------------aGCAGGTGTGGCATGTTTCT",
        )
        self.assertEqual(
            alignment.sequences[4].annotations["quality"], "99999999999999999999"
        )
        self.assertEqual(alignment.sequences[4].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[4].annotations["leftCount"], 400)
        self.assertEqual(alignment.sequences[4].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[4].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[5].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[5].seq), 174210431)
        self.assertEqual(
            alignment.sequences[5].seq[158040919 : 158040919 + 20],
            "AGAAACATGCCACACCTGCt",
        )
        self.assertEqual(
            alignment[5],
            "--------------------------------------------aGCAGGTGTGGCATGTTTCT",
        )
        self.assertEqual(alignment.sequences[5].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[5].annotations["leftCount"], 8044)
        self.assertEqual(alignment.sequences[5].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[5].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[6].id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(alignment.sequences[6].seq), 498454)
        self.assertEqual(
            alignment.sequences[6].seq[326906 : 326906 + 27],
            "AGCGATGCGCCACACCCACATTTCTAA",
        )
        self.assertEqual(
            alignment[6],
            "------------------------------------TTAGAAA-TGTGGGTGTGGCGCATCGCT",
        )
        self.assertEqual(
            alignment.sequences[6].annotations["quality"], "999999999999989998899999699"
        )
        self.assertEqual(alignment.sequences[6].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[6].annotations["leftCount"], 4145)
        self.assertEqual(alignment.sequences[6].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[6].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[7].id, "cavPor2.scaffold_216473")
        self.assertEqual(len(alignment.sequences[7].seq), 10026)
        self.assertEqual(
            alignment.sequences[7].seq[2184 : 2184 + 26], "GAACCATGCCACACCTAAATTTCTAA"
        )
        self.assertEqual(
            alignment[7],
            "------------------------------------TTAGAAA-TTTAGGTGTGGCATGGTTC-",
        )
        self.assertEqual(
            alignment.sequences[7].annotations["quality"], "42558311324566557465575854"
        )
        self.assertEqual(alignment.sequences[7].annotations["leftStatus"], "N")
        self.assertEqual(alignment.sequences[7].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[7].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[7].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[8].id, "loxAfr1.scaffold_75566")
        self.assertEqual(len(alignment.sequences[8].seq), 10574)
        self.assertEqual(
            alignment.sequences[8].seq[9203 : 9203 + 62],
            "aGAGATACAGCAAACCTCAATTTCTGAAACAACGTATTGATGTCCTAAACTTGCTCTCAAAC",
        )
        self.assertEqual(
            alignment[8],
            "GTT-TGAGAGCAAGTTTAGGACATCAATACGTTGTTTCAGAAA-TTGAGGTTTGCTGTATCTCt",
        )
        self.assertEqual(
            alignment.sequences[8].annotations["quality"],
            "99999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[8].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[8].annotations["leftCount"], 25)
        self.assertEqual(alignment.sequences[8].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[8].annotations["rightCount"], 0)
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "echTel1.scaffold_288249")
        self.assertEqual(len(record.seq), 100002)
        self.assertEqual(segment, (87661, 95225))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 9)
        self.assertEqual(len(alignment.annotations["empty"]), 1)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
            [[  3014795,   3014798,   3014799,   3014814,   3014815,   3014815,
                3014815,   3014815,   3014815,   3014821,   3014822,   3014822,
                3014841,   3014842],
             [   175207,    175207,    175206,    175191,    175190,    175180,
                 175180,    175176,    175175,    175169,    175168,    175167,
                 175148,    175147],
             [     6365,      6365,      6365,      6365,      6366,      6376,
                   6376,      6376,      6377,      6383,      6383,      6384,
                   6403,      6404],
             [155028706, 155028706, 155028706, 155028706, 155028706, 155028706,
              155028706, 155028706, 155028706, 155028706, 155028706, 155028706,
              155028687, 155028686],
             [157518758, 157518758, 157518758, 157518758, 157518758, 157518758,
              157518758, 157518758, 157518758, 157518758, 157518758, 157518758,
              157518739, 157518738],
             [158040939, 158040939, 158040939, 158040939, 158040939, 158040939,
              158040939, 158040939, 158040939, 158040939, 158040939, 158040939,
              158040920, 158040919],
             [   326933,    326933,    326933,    326933,    326933,    326933,
                 326933,    326933,    326933,    326927,    326926,    326926,
                 326907,    326906],
             [     2210,      2210,      2210,      2210,      2210,      2210,
                   2210,      2210,      2210,      2204,      2203,      2203,
                   2184,      2184],
             [     9265,      9262,      9262,      9247,      9246,      9236,
                   9235,      9231,      9230,      9224,      9223,      9223,
                   9204,      9203],
            ]
                    # fmt: on
                ),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 320596)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3014842 : 3014842 + 186],
            "CTTGGGATGCTTTATAGTGGAAATGGAAAGCAATTTATTTAGATCTTAAATCATTTTGAAGGTTAATAAAATGACCATATTAATATTCCCATGAACAAAGCCTTCATTTTTAAAATATTGCATCCTATAATACACATAAATCTTGTTCTCGtttttatttttttatttatttatttttttttcttt",
        )
        self.assertEqual(
            alignment[0],
            "C--TTGGGA---------TGCTTTATAGTGGAAATGGAAAGCA----A-TTTATTTAGATCTTAAATCATTTT-GAAGGTTAATAAAATGACCATATTAATATTCCCATGAACAAAGCCTTCATTT----TTAAAATATTGCATCCTATAATACACATAA-ATCTTGT-----TCTCGtttttatttttt----tatt-tat-----------------------ttattttttttt------------cttt",
        )
        self.assertEqual(alignment.sequences[1].id, "loxAfr1.scaffold_75566")
        self.assertEqual(len(alignment.sequences[1].seq), 10574)
        self.assertEqual(
            alignment.sequences[1].seq[8958 : 8958 + 245],
            "AGATTTGGCCAAAGTCAATGAAAAAAGAGAGAGAGAGTGAACTTGCTAAGACACCTGCTTAAAAAGGGAATGAGATTTTGAAAAGATGCTGTGTGTGGTATAAAACCCAATTTTTTTTTTTAAATGAGGGTATTGTTCACAGGAATATTAAAGTGAAAATTTCATTATACTTCAAAGGGATTTATGGCCAAAAGAAACAGTGTGACTTTCACTTCAGCtttaaaaaaaaaaaaaaatcaaaaata",
        )
        self.assertEqual(
            alignment[1],
            "tatttttgatttttttttttttttaaaGCTGAAGTGAAAGTCACACTG-TTTCTTTTGGCCATAAATCCCTTT-GAAGTATAATGAAATTTTCACTTTAATATTCCTGTGAACAATACCCTCATTT-AAAAAAAAAAATTGGGTTTTATACCACACACAGCATCTTTTCAAAATCTCATTCCC-TTTTTAAGCAGGTG-TCT---TAGCAAGTTCACTCTCTCTCTCTTTTTTCATTGACTTTGGCCAAATCT",
        )
        self.assertEqual(
            alignment.sequences[1].annotations["quality"],
            "99999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[2].id, "cavPor2.scaffold_216473")
        self.assertEqual(len(alignment.sequences[2].seq), 10026)
        self.assertEqual(
            alignment.sequences[2].seq[1986 : 1986 + 198],
            "agagaaagagagactaaTAAACCTGCCAGAAGAAGGAAATGAGATGTTGGAAATATGCTGTATGTGGTACAAAATCCAAACTTTTTTAAATGAGGGCATTGTTAGTAGGAATATTAATGTGATCATTTTATTATACTTCAAAAGGATTTAAGTCCTAAAGAAACAGGGTGCCTTTCATTCCAACTGTGAATATCCAAA",
        )
        self.assertEqual(
            alignment[2],
            "---TTTGGA---------TA-TTCACAGTTGGAATGAAAGGCACCCTG-TTTCTTTAGGACTTAAATCCTTTT-GAAGTATAATAAAATGATCACATTAATATTCCTACTAACAATGCCCTCATTT----AAAAAAGTTTGGATTTTGTACCACATACAGCATATTTCCAACATCTCATTTCCTTCTTCTGGCAGGTT-TAt-----------------------tagtctctcttt------------ctct",
        )
        self.assertEqual(
            alignment.sequences[2].annotations["quality"],
            "610137772001955312362668764253688587789879568878689568989568988778987788768588885664786777656586678636299978766899797899369899566878676899958889788869976598977898999989967788999979899987999997779899",
        )
        self.assertEqual(alignment.sequences[2].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[2].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[3].id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(alignment.sequences[3].seq), 498454)
        self.assertEqual(
            alignment.sequences[3].seq[326691 : 326691 + 215],
            "AAGTCAATAGACAAAAAGGGAACTTGCTAATAAACCTGCCAAAAAAGGAAATGAGATTTTGGAAATATGCTGTACGTGGTATAAAATCCCAGATTTTGTAAAATGAGGGCATTGTTCACAGGAATAGTAAAGTGATCATTTTATTATACTTCAAAAGGATTTAAGACCTACAGACACAGTGTGCTTTTTATTTCAGCTGTAAAAAATCCAAAAGG",
        )
        self.assertEqual(
            alignment[3],
            "CCTTTTGGA---------TTTTTTACAGCTGAAATAAAAAGCACACTG-TGTCTGTAGGTCTTAAATCCTTTT-GAAGTATAATAAAATGATCACTTTACTATTCCTGTGAACAATGCCCTCATTT---TACAAAATCTGGGATTTTATACCACGTACAGCATATTTCCAAAATCTCATTTCC-TTTTTTGGCAGGTT-TAT---TAGCAAGTTCCCT-------TTTTGTCTATTG------------ACTT",
        )
        self.assertEqual(
            alignment.sequences[3].annotations["quality"],
            "98999999989999999899999999999999999999999999999999999999999999999999989999999999999999999999999998999999999999999999988999999999739999999989999999999999999799999999999999999999769984999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[3].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[3].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[3].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[3].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[4].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[4].seq), 174210431)
        self.assertEqual(
            alignment.sequences[4].seq[158040700 : 158040700 + 219],
            "AAGTAAATAAGAGAGAGAGATAACTTGCTATAAATAAACCTGCCAAGAAAGGAAATGAGATTTCAGAAATATGCTGCATGTGGTATAAAATCCAGGATTTTTTTAAAAGAGGACATTGCTCACAGGAGTATTAAAGTGATCGTTTTATTATACTTCAAAAGTATTTAAGACCTAAAGAAACAGTGTGCCTTCCATTTCAGCTATAAAAAGTCCAAAAGG",
        )
        self.assertEqual(
            alignment[4],
            "CCTTTTGGA---------CTTTTTATAGCTGAAATGGAAGGCACACTG-TTTCTTTAGGTCTTAAATACTTTT-GAAGTATAATAAAACGATCACTTTAATACTCCTGTGAGCAATGTCCTCTTTT---AAAAAAATCCTGGATTTTATACCACATGCAGCATATTTCTGAAATCTCATTTCC-TTTCTTGGCAGGTT-TATTTATAGCAAGTTATCTC------TCTCTCTTATTT------------ACTT",
        )
        self.assertEqual(alignment.sequences[4].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[4].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[4].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[4].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[5].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[5].seq), 173908612)
        self.assertEqual(
            alignment.sequences[5].seq[157518570 : 157518570 + 168],
            "AAATGGGATTTCAGAAATATGCTGCATGTGGTATAAAATCCAGGATTTTTTTAAAAAGAGGGCATTGTTCACAGGAGTATTAAAGTGATCATTTTTTTATACTTCAAAAGTATTTAAGACCTAAAGAAATATTGTGCCTTCCATTTCAGCTATAAAAAGTCCAAAAGG",
        )
        self.assertEqual(
            alignment[5],
            "CCTTTTGGA---------CTTTTTATAGCTGAAATGGAAGGCACAATA-TTTCTTTAGGTCTTAAATACTTTT-GAAGTATAAAAAAATGATCACTTTAATACTCCTGTGAACAATGCCCTCTTTTT--AAAAAAATCCTGGATTTTATACCACATGCAGCATATTTCTGAAATCCCATTT------------------------------------------------------------------------",
        )
        self.assertEqual(
            alignment.sequences[5].annotations["quality"],
            "999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[5].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[5].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[5].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[5].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[6].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[6].seq), 170899992)
        self.assertEqual(
            alignment.sequences[6].seq[155028517 : 155028517 + 169],
            "AAATGAGATTTCAGAAATATGCTGCATGTGGTATAAAATCCAGGATTTTTTTTAAAAAGAGGGCATTGTTCACAGGAGTATTAAAGTGATCATTTTTTTATACTTCAAAAGTATTTAAGACCTAAAGAAATATTGTGCCTTCCATTTCAGCTATAAAAAGTCCAAAAGG",
        )
        self.assertEqual(
            alignment[6],
            "CCTTTTGGA---------CTTTTTATAGCTGAAATGGAAGGCACAATA-TTTCTTTAGGTCTTAAATACTTTT-GAAGTATAAAAAAATGATCACTTTAATACTCCTGTGAACAATGCCCTCTTTTT-AAAAAAAATCCTGGATTTTATACCACATGCAGCATATTTCTGAAATCTCATTT------------------------------------------------------------------------",
        )
        self.assertEqual(alignment.sequences[6].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[6].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[6].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[6].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[7].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[7].seq), 133105)
        self.assertEqual(
            alignment.sequences[7].seq[6404 : 6404 + 223],
            "TCTTTTGTATTTTTTTATAGCTGAAAAGGAAGGCACACTGTCTCTTTAGGTCTTAAATACGTTTGAAGCATAATAAAATGATCACTTTCATACTCCTGTGAATAATGCCCTCCTTTTAAAAAGAAATCTTGGATTTTATATCACACGCAGCATATTTCTGAAATCTCATTTCCTTTCCTGGCAGGTTTATTAATAGCAAGTTCTCTCTCTTATTTATTTATTT",
        )
        self.assertEqual(
            alignment[7],
            "TCTTTTGTA--------TTTTTTTATAGCTGAAAAGGAAGGCACACTG-TCTCTTTAGGTCTTAAATACGTTT-GAAGCATAATAAAATGATCACTTTCATACTCCTGTGAATAATGCCCTCCTTTTAAAAAGAAATCTTGGATTTTATATCACACGCAGCATATTTCTGAAATCTCATTTCC-TTTCCTGGCAGGTT-TATTAATAGCAAGTTCTCTC------TCTTATTTATTT------------ATTT",
        )
        self.assertEqual(alignment.sequences[7].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[7].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[7].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[7].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[8].id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(alignment.sequences[8].seq), 359464)
        self.assertEqual(
            alignment.sequences[8].seq[174925 : 174925 + 222],
            "GAGAAAAAGAAGAAGAAAGCCGGAGAGAACTTGCTAACAAACCCGCCAAAAAGGAAATGAGATTTTGGAAATATGCTGGATCTGGTATAAAATCCAAGATTTTTTTAAATGAGGGCATTGTTCACAGGAACAGTAAAGTAATCGTTTTATTATACTTCAAAAGGATTTAAGTCCTAAAAAAACAGTGTGCCTTTCATTTCACCCATAAAAAAGCCCAAAATA",
        )
        self.assertEqual(
            alignment[8],
            "TATTTTGGG--------CTTTTTTATGGGTGAAATGAAAGGCACACTG-TTTTTTTAGGACTTAAATCCTTTT-GAAGTATAATAAAACGATTACTTTACTGTTCCTGTGAACAATGCCCTCATTT---AAAAAAATCTTGGATTTTATACCAGATCCAGCATATTTCCAAAATCTCATTT-C-CTTTTTGGCGGGTT-TGT---TAGCAAGTTCTCTCCGGCTTTCTTCTTCTTTT------------TCTC",
        )
        self.assertEqual(
            alignment.sequences[8].annotations["quality"],
            "999999999999999999999999999999999799999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[8].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[8].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[8].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[8].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[9].id, "ornAna1.chr2")
        self.assertEqual(len(alignment.sequences[9].seq), 54797317)
        self.assertEqual(
            alignment.sequences[9].seq[40046122 : 40046122 + 201],
            "AAAGAGAGAACGGTTAATACAATGCACCCAAAAAGGAACTAAGATTTTGGAAATGGGCTGTATGTCGTTAAAATACAAAGATATTTTAATGAGGGCTTTGTTAGCGGGGAAATGAAGATCATAATTACAATACACTTTTAAAAAGGCTTAAGATAGAAAGAAAACAATGAGCCTTTCACTTTTGCAGTAACATTTCCAAGG",
        )
        self.assertEqual(
            alignment[9],
            "--CCTTGGA---------AATGTTACTGCAAAAGTGAAAGGCTCATTGTTTTCTTTCTATCTTAAGCCTTTTTAAAAGTGTATTGTAATTATGATCTTCATTTCCCCGCTAACAAAGCCCTCATTA----AAATATCTTTGTATTTTA-ACGACATACAGCCCATTTCCAAAATCTTAGTTCC-TTTTTGGGTGCATTGTAT---TAAC----------------CGTTCTCTCTTT----------------",
        )
        self.assertEqual(alignment.sequences[9].annotations["leftStatus"], "N")
        self.assertEqual(alignment.sequences[9].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[9].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[9].annotations["rightCount"], 5690)
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "echTel1.scaffold_288249")
        self.assertEqual(len(record.seq), 100002)
        self.assertEqual(segment, (87661, 95225))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 10)
        self.assertEqual(len(alignment.annotations["empty"]), 1)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
            [[  3014842,   3014843,   3014843,   3014843,   3014849,   3014849,
                3014849,   3014851,   3014852,   3014874,   3014874,   3014875,
                3014875,   3014899,   3014899,   3014951,   3014951,   3014951,
                3014951,   3014951,   3014969,   3014970,   3014981,   3014981,
                3014988,   3014988,   3014996,   3014997,   3014998,   3014999,
                3015005,   3015005,   3015009,   3015009,   3015012,   3015012,
                3015012,   3015012,   3015012,   3015012,   3015024,   3015024,
                3015028],
             [     9203,      9202,      9201,      9200,      9194,      9186,
                   9185,      9183,      9182,      9160,      9156,      9155,
                   9155,      9131,      9131,      9079,      9079,      9078,
                   9077,      9076,      9058,      9057,      9046,      9045,
                   9038,      9033,      9025,      9024,      9023,      9023,
                   9017,      9013,      9009,      9009,      9006,      9006,
                   9002,      8993,      8992,      8986,      8974,      8962,
                   8958],
             [     2184,      2184,      2184,      2184,      2178,      2178,
                   2178,      2176,      2176,      2154,      2150,      2149,
                   2149,      2125,      2125,      2073,      2073,      2073,
                   2073,      2073,      2055,      2054,      2043,      2042,
                   2035,      2030,      2022,      2021,      2020,      2019,
                   2013,      2009,      2005,      2005,      2002,      2002,
                   2002,      2002,      2002,      2002,      1990,      1990,
                   1986],
             [   326906,    326905,    326904,    326903,    326897,    326897,
                 326897,    326895,    326894,    326872,    326868,    326867,
                 326867,    326843,    326843,    326791,    326791,    326791,
                 326791,    326790,    326772,    326771,    326760,    326759,
                 326752,    326747,    326739,    326738,    326737,    326737,
                 326731,    326727,    326723,    326723,    326720,    326720,
                 326716,    326707,    326707,    326707,    326695,    326695,
                 326691],
             [158040919, 158040918, 158040917, 158040916, 158040910, 158040910,
              158040910, 158040908, 158040907, 158040885, 158040881, 158040880,
              158040880, 158040856, 158040856, 158040804, 158040804, 158040804,
              158040804, 158040803, 158040785, 158040784, 158040773, 158040772,
              158040765, 158040760, 158040752, 158040751, 158040750, 158040750,
              158040744, 158040740, 158040736, 158040736, 158040733, 158040730,
              158040726, 158040717, 158040716, 158040716, 158040704, 158040704,
              158040700],
             [157518738, 157518737, 157518736, 157518735, 157518729, 157518729,
              157518729, 157518727, 157518726, 157518704, 157518700, 157518699,
              157518699, 157518675, 157518675, 157518623, 157518622, 157518622,
              157518622, 157518621, 157518603, 157518602, 157518591, 157518590,
              157518583, 157518578, 157518570, 157518570, 157518570, 157518570,
              157518570, 157518570, 157518570, 157518570, 157518570, 157518570,
              157518570, 157518570, 157518570, 157518570, 157518570, 157518570,
              157518570],
             [155028686, 155028685, 155028684, 155028683, 155028677, 155028677,
              155028677, 155028675, 155028674, 155028652, 155028648, 155028647,
              155028647, 155028623, 155028623, 155028571, 155028570, 155028570,
              155028569, 155028568, 155028550, 155028549, 155028538, 155028537,
              155028530, 155028525, 155028517, 155028517, 155028517, 155028517,
              155028517, 155028517, 155028517, 155028517, 155028517, 155028517,
              155028517, 155028517, 155028517, 155028517, 155028517, 155028517,
              155028517],
             [     6404,      6405,      6406,      6407,      6413,      6413,
                   6414,      6416,      6417,      6439,      6443,      6444,
                   6444,      6468,      6468,      6520,      6521,      6522,
                   6523,      6524,      6542,      6543,      6554,      6555,
                   6562,      6567,      6575,      6576,      6577,      6577,
                   6583,      6587,      6591,      6591,      6594,      6597,
                   6601,      6610,      6611,      6611,      6623,      6623,
                   6627],
             [   175147,    175146,    175145,    175144,    175138,    175138,
                 175137,    175135,    175134,    175112,    175108,    175107,
                 175107,    175083,    175083,    175031,    175031,    175031,
                 175031,    175030,    175012,    175011,    175000,    174999,
                 174992,    174987,    174979,    174979,    174978,    174978,
                 174972,    174968,    174964,    174964,    174961,    174961,
                 174957,    174948,    174947,    174941,    174929,    174929,
                 174925],
             [ 40046323,  40046323,  40046323,  40046322,  40046316,  40046316,
               40046316,  40046314,  40046313,  40046291,  40046287,  40046286,
               40046285,  40046261,  40046260,  40046208,  40046208,  40046208,
               40046208,  40046208,  40046190,  40046190,  40046179,  40046178,
               40046171,  40046166,  40046158,  40046157,  40046156,  40046156,
               40046150,  40046146,  40046142,  40046141,  40046138,  40046138,
               40046134,  40046134,  40046134,  40046134,  40046122,  40046122,
               40046122],
            ]
                    # fmt: on
                ),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, -36127)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3015028 : 3015028 + 58],
            "ccattttttattaggtatttagctcatttacatttccaatgctataccaaaagtcccc",
        )
        self.assertEqual(
            alignment[0],
            "ccatttt----------ttattaggtatttagctcatttacatttccaatgctatac----caaaagtcccc",
        )
        self.assertEqual(alignment.sequences[1].id, "loxAfr1.scaffold_75566")
        self.assertEqual(len(alignment.sequences[1].seq), 10574)
        self.assertEqual(
            alignment.sequences[1].seq[8925 : 8925 + 33],
            "GAGAACTTTTGTAAGGAATGGAGGTAGAAGTGA",
        )
        self.assertEqual(
            alignment[1],
            "TCACTTCTA---------------------------------------CCTCCATTCCTTACAAAAGTTCTC",
        )
        self.assertEqual(
            alignment.sequences[1].annotations["quality"],
            "999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "N")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[2].id, "cavPor2.scaffold_216473")
        self.assertEqual(len(alignment.sequences[2].seq), 10026)
        self.assertEqual(alignment.sequences[2].seq[1978 : 1978 + 8], "agacacag")
        self.assertEqual(
            alignment[2],
            "ctgtgtc----------t------------------------------------------------------",
        )
        self.assertEqual(alignment.sequences[2].annotations["quality"], "67889899")
        self.assertEqual(alignment.sequences[2].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[2].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[2].annotations["rightCount"], 1372)
        self.assertEqual(alignment.sequences[3].id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(alignment.sequences[3].seq), 498454)
        self.assertEqual(alignment.sequences[3].seq[326690 : 326690 + 1], "A")
        self.assertEqual(
            alignment[3],
            "T-----------------------------------------------------------------------",
        )
        self.assertEqual(alignment.sequences[3].annotations["quality"], "9")
        self.assertEqual(alignment.sequences[3].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[3].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[3].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[3].annotations["rightCount"], 2374)
        self.assertEqual(alignment.sequences[4].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[4].seq), 174210431)
        self.assertEqual(
            alignment.sequences[4].seq[158040688 : 158040688 + 12], "GAAAAGGCAAAA"
        )
        self.assertEqual(
            alignment[4],
            "TTTTGCCTTTTC------------------------------------------------------------",
        )
        self.assertEqual(alignment.sequences[4].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[4].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[4].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[4].annotations["rightCount"], 75)
        self.assertEqual(alignment.sequences[5].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[5].seq), 133105)
        self.assertEqual(
            alignment.sequences[5].seq[6627 : 6627 + 22], "TGCTTTTTCTCAAATCTCCACT"
        )
        self.assertEqual(
            alignment[5],
            "TGCTTTTTCTCAAATCTCCACT--------------------------------------------------",
        )
        self.assertEqual(alignment.sequences[5].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[5].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[5].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[5].annotations["rightCount"], 3479)
        self.assertEqual(alignment.sequences[6].id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(alignment.sequences[6].seq), 359464)
        self.assertEqual(
            alignment.sequences[6].seq[174889 : 174889 + 36],
            "TGGAAGTGAATATTTGGCAAAAGTCAATAGAAAAAA",
        )
        self.assertEqual(
            alignment[6],
            "TTTTTTCTATT------------GACTTTTGCCAAATATTCACTTCCA------------------------",
        )
        self.assertEqual(
            alignment.sequences[6].annotations["quality"],
            "999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[6].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[6].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[6].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[6].annotations["rightCount"], 137)
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ornAna1.chr2")
        self.assertEqual(len(record.seq), 54797317)
        self.assertEqual(segment, (40046122, 40040432))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "echTel1.scaffold_288249")
        self.assertEqual(len(record.seq), 100002)
        self.assertEqual(segment, (87661, 95225))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][2]
        (record, segment, status) = empty
        self.assertEqual(record.id, "hg18.chr6")
        self.assertEqual(len(record.seq), 170899992)
        self.assertEqual(segment, (155028517, 155028517))
        self.assertEqual(status, "C")
        empty = alignment.annotations["empty"][3]
        (record, segment, status) = empty
        self.assertEqual(record.id, "panTro2.chr6")
        self.assertEqual(len(record.seq), 173908612)
        self.assertEqual(segment, (157518570, 157518570))
        self.assertEqual(status, "C")
        self.assertEqual(len(alignment.sequences), 7)
        self.assertEqual(len(alignment.annotations["empty"]), 4)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
            [[  3015028,   3015029,   3015035,   3015035,   3015035,   3015035,
                3015035,   3015036,   3015040,   3015041,   3015066,   3015075,
                3015075,   3015086],
             [     8958,      8957,      8951,      8949,      8949,      8949,
                   8949,      8949,      8949,      8949,      8949,      8940,
                   8936,      8925],
             [     1986,      1985,      1979,      1979,      1979,      1979,
                   1979,      1978,      1978,      1978,      1978,      1978,
                   1978,      1978],
             [   326691,    326690,    326690,    326690,    326690,    326690,
                 326690,    326690,    326690,    326690,    326690,    326690,
                 326690,    326690],
             [158040700, 158040699, 158040693, 158040691, 158040689, 158040688,
              158040688, 158040688, 158040688, 158040688, 158040688, 158040688,
              158040688, 158040688],
             [     6627,      6628,      6634,      6636,      6638,      6639,
                   6644,      6645,      6649,      6649,      6649,      6649,
                   6649,      6649],
             [   174925,    174924,    174918,    174916,    174914,    174914,
                 174914,    174914,    174914,    174914,    174889,    174889,
                 174889,    174889],
            ]
                    # fmt: on
                ),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3015086 : 3015086 + 2572],
            "catacccacccacccccactcccctacccgcccactccccctttttggccctggcgttcccctgttctggggcatataaagtttgtgtgtccaatgggcctctctttccagtgatggccgactaggccatcttttgatacatatgcagctagagtcaagagctccggggtactggttagttcataatgttgatccacctatagggttgcagatccctttagctccttgggtactttctctagctcccccattgggagccctgtgatccatccattagctgactgtgggcatccacttctgtgtttgctaggccccggcatagtctcacaagagacagctacatctgggtcctttcgataagatcttgctagtgtatgcaatggtgtcagcgtttggatgctgattatggggtagatccctggataaggcagtctctacatggtccatcctttcatctcagctccaaactttgtctctgtaactccttccaagggtgttttgttcccacttctaaggaggggcatagtgtccacacttcagtcttcatttttcttgagtttcatgtgtttaggaaattgtatcttatatcttgggtatcctaggttttgggctaatatccacttatcagtgagtacatattgtgtgagttcctttgtgaatgtgttacctcactcaggatgatgccctccaggtccatccatttggctaggaatttcataaattaattctttttaatagctgagtagtactccattgtgtagatgtaccacattttctgtatccattcctctgttgaggggcatctgggttctttccagcttctggctattataaataaggctgctatgaacatagtggagcatgtgtccttcttaccagttggggcatcttttggatatatgcccaggagaggtattgctggatcctccggtagtactatgtccaattttctgaggaaccgccagacggatttccagagtggttgtacaagcctgcaatcccaccaacaatggaggagtgttcctctttctccacatcctcgccagcatctgctgtcacctgaatttttgatcttagccattctgactggtgtgaggtggaatctcagggttgttttgatttgcatttctctgatgattaaggatgttgaacatgttttcaggtgcttctctgccattcggtattcctcaggtgagaattctttgttcagttctgagccccattttttaatggggttatttgattttctgaagtccaccttcttgagttctttatatatgttggatattagtcccctatctgatttaggataggtaaagatcctttcccaatctgttggtggtctctttgtgttattgacggtgtcttttgccttgcagaaactttggagtttcattaggtcccatttgtcaattctcgatcttacagcacaagccattgctgttctgttcaggaatttttcccctgtgcccatatcttcaaggcttttccccactttctcctctataagtttcagtgtctctggttttatgtggagttctttgatccatttagatttgaccttagtacaaggagataagtatggatcgattcgcattcttctacatgataacaaccagttgtgccagcaccaattgttgaaaatgctgtctttcttccactggatggttttagctcccttgtcgaagatcaagtgaccataggtgtgtgggttcatttctgggtcttcaattctattccattggtctacttgtctgtctctataccagtaccatgcagtttttaccacaattgctctgtagtaaagctttaggtcaggcatggtgattccaccagaggttcttttatccttgagaagagtttttgctatcctaggttttttgttattccagatgaatttgcaaattgctccttctaattcgttgaagaattgagttggaattgtgatggggattgcattgaatctgtagattgcttttggcaagatagccatttttacaatgttgatcctgccaatccatgagcatgggagagctttccatcttctgagatcttctttaatttctttcttcagagacttgaagtttttatcatacagatctttcacttccttagttagagtcacgccgagatattttatattatttgtgactattgagaagggtgttgtttccctaatttctttctcagcctgtttattctttgtgtagagaaaggccattgacttgtttgagttaattttatatccagctacttcaccgaagctgtttatcaggtttaggagttctctggtggaatttttagggtcacttatatatactatcatatcatctgcaaaaagtgatattttgacttcctcctttccaatttgtatccccttgatctccttttgttgtcgaattgctctggctaatacttcaagtactatgttgaaaaggtagggagaaagtgggcagccttgtctagtccctgattttagtgagattgcttccagcttctctccatttactttgatgttggctactggtttgctgtagattgcttttatcatgtttaggtatgggTGTTCTCG",
        )
        self.assertEqual(
            alignment[0],
            "catacccacccacccccactcccctacccgcccactccccctttttggccctggcgttcccctgttctggggcatataaagtttgtgtgtccaatgggcctctctttccagtgatggccgactaggccatcttttgatacatatgcagctagagtcaagagctccggggtactggttagttcataatgttgatccacctatagggttgcagatccctttagctccttgggtactttctctagctcccccattgggagccctgtgatccatccattagctgactgtgggcatccacttctgtgtttgctaggccccggcatagtctcacaagagacagctacatctgggtcctttcgataagatcttgctagtgtatgcaatggtgtcagcgtttggatgctgattatggggtagatccctggataaggcagtctctacatggtccatcctttcatctcagctccaaactttgtctctgtaactccttccaagggtgttttgttcccacttctaaggaggggcatagtgtccacacttcagtcttcatttttcttgagtttcatgtgtttaggaaattgtatcttatatcttgggtatcctaggttttgggctaatatccacttatcagtgagtacatattgtgtgagttcctttgtgaatgtgttacctcactcaggatgatgccctccaggtccatccatttggctaggaatttcataaattaattctttttaatagctgagtagtactccattgtgtagatgtaccacattttctgtatccattcctctgttgaggggcatctgggttctttccagcttctggctattataaataaggctgctatgaacatagtggagcatgtgtccttcttaccagttggggcatcttttggatatatgcccaggagaggtattgctggatcctccggtagtactatgtccaattttctgaggaaccgccagacggatttccagagtggttgtacaagcctgcaatcccaccaacaatggaggagtgttcctctttctccacatcctcgccagcatctgctgtcacctgaatttttgatcttagccattctgactggtgtgaggtggaatctcagggttgttttgatttgcatttctctgatgattaaggatgttgaacatgttttcaggtgcttctctgccattcggtattcctcaggtgagaattctttgttcagttctgagccccattttttaatggggttatttgattttctgaagtccaccttcttgagttctttatatatgttggatattagtcccctatctgatttaggataggtaaagatcctttcccaatctgttggtggtctctttgtgttattgacggtgtcttttgccttgcagaaactttggagtttcattaggtcccatttgtcaattctcgatcttacagcacaagccattgctgttctgttcaggaatttttcccctgtgcccatatcttcaaggcttttccccactttctcctctataagtttcagtgtctctggttttatgtggagttctttgatccatttagatttgaccttagtacaaggagataagtatggatcgattcgcattcttctacatgataacaaccagttgtgccagcaccaattgttgaaaatgctgtctttcttccactggatggttttagctcccttgtcgaagatcaagtgaccataggtgtgtgggttcatttctgggtcttcaattctattccattggtctacttgtctgtctctataccagtaccatgcagtttttaccacaattgctctgtagtaaagctttaggtcaggcatggtgattccaccagaggttcttttatccttgagaagagtttttgctatcctaggttttttgttattccagatgaatttgcaaattgctccttctaattcgttgaagaattgagttggaattgtgatggggattgcattgaatctgtagattgcttttggcaagatagccatttttacaatgttgatcctgccaatccatgagcatgggagagctttccatcttctgagatcttctttaatttctttcttcagagacttgaagtttttatcatacagatctttcacttccttagttagagtcacgccgagatattttatattatttgtgactattgagaagggtgttgtttccctaatttctttctcagcctgtttattctttgtgtagagaaaggccattgacttgtttgagttaattttatatccagctacttcaccgaagctgtttatcaggtttaggagttctctggtggaatttttagggtcacttatatatactatcatatcatctgcaaaaagtgatattttgacttcctcctttccaatttgtatccccttgatctccttttgttgtcgaattgctctggctaatacttcaagtactatgttgaaaaggtagggagaaagtgggcagccttgtctagtccctgattttagtgagattgcttccagcttctctccatttactttgatgttggctactggtttgctgtagattgcttttatcatgtttaggtatgggTGTTCTCG",
        )
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ornAna1.chr2")
        self.assertEqual(len(record.seq), 54797317)
        self.assertEqual(segment, (40046122, 40040432))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "calJac1.Contig6394")
        self.assertEqual(len(record.seq), 133105)
        self.assertEqual(segment, (6649, 10128))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][2]
        (record, segment, status) = empty
        self.assertEqual(record.id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(record.seq), 498454)
        self.assertEqual(segment, (326690, 324316))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][3]
        (record, segment, status) = empty
        self.assertEqual(record.id, "echTel1.scaffold_288249")
        self.assertEqual(len(record.seq), 100002)
        self.assertEqual(segment, (87661, 95225))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][4]
        (record, segment, status) = empty
        self.assertEqual(record.id, "cavPor2.scaffold_216473")
        self.assertEqual(len(record.seq), 10026)
        self.assertEqual(segment, (1978, 606))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][5]
        (record, segment, status) = empty
        self.assertEqual(record.id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(record.seq), 359464)
        self.assertEqual(segment, (174889, 174752))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][6]
        (record, segment, status) = empty
        self.assertEqual(record.id, "hg18.chr6")
        self.assertEqual(len(record.seq), 170899992)
        self.assertEqual(segment, (155028517, 155028517))
        self.assertEqual(status, "C")
        empty = alignment.annotations["empty"][7]
        (record, segment, status) = empty
        self.assertEqual(record.id, "panTro2.chr6")
        self.assertEqual(len(record.seq), 173908612)
        self.assertEqual(segment, (157518570, 157518570))
        self.assertEqual(status, "C")
        empty = alignment.annotations["empty"][8]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ponAbe2.chr6")
        self.assertEqual(len(record.seq), 174210431)
        self.assertEqual(segment, (158040688, 158040613))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 1)
        self.assertEqual(len(alignment.annotations["empty"]), 9)
        self.assertTrue(
            numpy.array_equal(alignment.coordinates, numpy.array([[3015086, 3017658]]))
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 12170)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3017658 : 3017658 + 85],
            "TTTTTATTTGCAGGTTTCTTTACAGTTCTCTTTCATTCTTCTCCTCTTTTCTTCTGTTGACCTTTATCAGATTTCTGCTTTAACC",
        )
        self.assertEqual(
            alignment[0],
            "TTTTTATTTGCAGGTTTCTTTAC----AGTTCTCTTTCATTCTTCTCCTCTTTTCTTCTGTTGACCTTTATCAGATTTCTGCTTTAACC",
        )
        self.assertEqual(alignment.sequences[1].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[1].seq), 170899992)
        self.assertEqual(
            alignment.sequences[1].seq[155028434 : 155028434 + 83],
            "GGTGGAAGTGGAGATTTGGGAAAAGGCAAAAAAATAAATAAGAGAGAGAGATAACTTGCTATAAATAACCCTGCCAAGAAAGG",
        )
        self.assertEqual(
            alignment[1],
            "CCTTTCTTGGCAGGGTTATTTATAGCAAGTTATCTCTCTCTCTTA------TTTATTTTTTTGCCTTTTCCCAAATCTCCACTTCCACC",
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 53)
        self.assertEqual(alignment.sequences[2].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[2].seq), 173908612)
        self.assertEqual(
            alignment.sequences[2].seq[157518487 : 157518487 + 83],
            "GGTGGAAGTGGAGATTTGGGAAAAGGCAAAAAAATAAATAAGAGAGAGAGATAACTTGCTATAAATAACCCTGCCAAGAAAGG",
        )
        self.assertEqual(
            alignment[2],
            "CCTTTCTTGGCAGGGTTATTTATAGCAAGTTATCTCTCTCTCTTA------TTTATTTTTTTGCCTTTTCCCAAATCTCCACTTCCACC",
        )
        self.assertEqual(
            alignment.sequences[2].annotations["quality"],
            "99999999999999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[2].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[2].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[2].annotations["rightCount"], 53)
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ornAna1.chr2")
        self.assertEqual(len(record.seq), 54797317)
        self.assertEqual(segment, (40046122, 40040432))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "calJac1.Contig6394")
        self.assertEqual(len(record.seq), 133105)
        self.assertEqual(segment, (6649, 10128))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][2]
        (record, segment, status) = empty
        self.assertEqual(record.id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(record.seq), 498454)
        self.assertEqual(segment, (326690, 324316))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][3]
        (record, segment, status) = empty
        self.assertEqual(record.id, "echTel1.scaffold_288249")
        self.assertEqual(len(record.seq), 100002)
        self.assertEqual(segment, (87661, 95225))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][4]
        (record, segment, status) = empty
        self.assertEqual(record.id, "cavPor2.scaffold_216473")
        self.assertEqual(len(record.seq), 10026)
        self.assertEqual(segment, (1978, 606))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][5]
        (record, segment, status) = empty
        self.assertEqual(record.id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(record.seq), 359464)
        self.assertEqual(segment, (174889, 174752))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][6]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ponAbe2.chr6")
        self.assertEqual(len(record.seq), 174210431)
        self.assertEqual(segment, (158040688, 158040613))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 3)
        self.assertEqual(len(alignment.annotations["empty"]), 7)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
            [[  3017658,   3017681,   3017681,   3017699,   3017705,   3017743],
             [155028517, 155028494, 155028490, 155028472, 155028472, 155028434],
             [157518570, 157518547, 157518543, 157518525, 157518525, 157518487],
            ]
                    # fmt: on
                ),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3017743 : 3017743 + 418],
            "ACCACAGACCTTCTGTTTAGTCCAAAGGACGCAAATTATGTATCCACTTtagtaggaggctgacccgcagcctacatgaaccaggtatttctggaaggcaggctggggttgaaagagaaattagatggtgagaaaagaataatgaggccaagacaaatttttctcttatcaaggcccaagagagtttactaagagactatgcttaaaagggggaaggcccatcccccccccccctcgcgccagtctatccttggtgctttgtcaccatgccatcagcacttggtcggcaggtagcagaatctcagggcagttgacacttcaaaagaaaccagccaagtcagaaagctgcactgcaggagacctgcactcagtggtgacaaggtctgtaccagcctgcttcaggctgggggaggctaca",
        )
        self.assertEqual(
            alignment[0],
            "ACCACAGACCTTCTGTTTAGTCCAAAGGACGCAAATTATGTATCCACTTtagtaggaggctgacccgcagcctacatgaaccaggtatttctggaaggcaggctggggttgaaagagaaattagatggtgagaaaagaataatgaggccaagacaaatttttctcttatcaaggcccaagagagtttactaagagactatgcttaaaagggggaaggcccatcccccccccccctcgcgccagtctatccttggtgctttgtcaccatgccatcagcacttggtcggcaggtagcagaatctcagggcagttgacacttcaaaagaaaccagccaagtcagaaagctgcactgcaggagacctgcactcagtggtgacaaggtctgtaccagcctgcttcaggctgggggaggctaca",
        )
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ornAna1.chr2")
        self.assertEqual(len(record.seq), 54797317)
        self.assertEqual(segment, (40046122, 40040432))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "calJac1.Contig6394")
        self.assertEqual(len(record.seq), 133105)
        self.assertEqual(segment, (6649, 10128))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][2]
        (record, segment, status) = empty
        self.assertEqual(record.id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(record.seq), 498454)
        self.assertEqual(segment, (326690, 324316))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][3]
        (record, segment, status) = empty
        self.assertEqual(record.id, "echTel1.scaffold_288249")
        self.assertEqual(len(record.seq), 100002)
        self.assertEqual(segment, (87661, 95225))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][4]
        (record, segment, status) = empty
        self.assertEqual(record.id, "cavPor2.scaffold_216473")
        self.assertEqual(len(record.seq), 10026)
        self.assertEqual(segment, (1978, 606))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][5]
        (record, segment, status) = empty
        self.assertEqual(record.id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(record.seq), 359464)
        self.assertEqual(segment, (174889, 174752))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][6]
        (record, segment, status) = empty
        self.assertEqual(record.id, "hg18.chr6")
        self.assertEqual(len(record.seq), 170899992)
        self.assertEqual(segment, (155028434, 155028381))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][7]
        (record, segment, status) = empty
        self.assertEqual(record.id, "panTro2.chr6")
        self.assertEqual(len(record.seq), 173908612)
        self.assertEqual(segment, (157518487, 157518434))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][8]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ponAbe2.chr6")
        self.assertEqual(len(record.seq), 174210431)
        self.assertEqual(segment, (158040688, 158040613))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 1)
        self.assertEqual(len(alignment.annotations["empty"]), 9)
        self.assertTrue(
            numpy.array_equal(alignment.coordinates, numpy.array([[3017743, 3018161]]))
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 22499)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3018161 : 3018161 + 69],
            "ATCCACAAAAGAGACAAAGAAGAAAACCAAAAGAAAAGATTGTAGCTTAAAACAATTCCATTTTATTGA",
        )
        self.assertEqual(
            alignment[0],
            "ATCCACAAAAGAGAC-----AAAGAAGAAAACCAAAAGAAAAGATTGTAGCTTAAAACAATTCCATTTTATTGA",
        )
        self.assertEqual(alignment.sequences[1].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[1].seq), 173908612)
        self.assertEqual(
            alignment.sequences[1].seq[157518369 : 157518369 + 65],
            "TCAATAAAATAAAATTGTTTTATGGTACAAACTTTTCTTATGGTTAATTGTTTTCCTTTGTAGGT",
        )
        self.assertEqual(
            alignment[1],
            "ACCTACAAAGG---------AAAACAATTAACCATAAGAAAAGTTTGTACCATAAAACAATTTTATTTTATTGA",
        )
        self.assertEqual(
            alignment.sequences[1].annotations["quality"],
            "99999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 53)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[2].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[2].seq), 170899992)
        self.assertEqual(
            alignment.sequences[2].seq[155028316 : 155028316 + 65],
            "TCAATAAAATAAAATTGTTTTATGGTACAAACTTTTCTTGTGGTTAATTGTTTTCCTTTGTAGGT",
        )
        self.assertEqual(
            alignment[2],
            "ACCTACAAAGG---------AAAACAATTAACCACAAGAAAAGTTTGTACCATAAAACAATTTTATTTTATTGA",
        )
        self.assertEqual(alignment.sequences[2].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[2].annotations["leftCount"], 53)
        self.assertEqual(alignment.sequences[2].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[3].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[3].seq), 174210431)
        self.assertEqual(
            alignment.sequences[3].seq[158040552 : 158040552 + 61],
            "TCAATAAAATAAAATTGTTTTATGGTACAAACTTTTCTTATGGTTAATTGTTTTCCTTTGT",
        )
        self.assertEqual(
            alignment[3],
            "-------------ACAAAGGAAAACAATTAACCATAAGAAAAGTTTGTACCATAAAACAATTTTATTTTATTGA",
        )
        self.assertEqual(alignment.sequences[3].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[3].annotations["leftCount"], 75)
        self.assertEqual(alignment.sequences[3].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[3].annotations["rightCount"], 97)
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ornAna1.chr2")
        self.assertEqual(len(record.seq), 54797317)
        self.assertEqual(segment, (40046122, 40040432))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "calJac1.Contig6394")
        self.assertEqual(len(record.seq), 133105)
        self.assertEqual(segment, (6649, 10128))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][2]
        (record, segment, status) = empty
        self.assertEqual(record.id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(record.seq), 498454)
        self.assertEqual(segment, (326690, 324316))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][3]
        (record, segment, status) = empty
        self.assertEqual(record.id, "echTel1.scaffold_288249")
        self.assertEqual(len(record.seq), 100002)
        self.assertEqual(segment, (87661, 95225))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][4]
        (record, segment, status) = empty
        self.assertEqual(record.id, "cavPor2.scaffold_216473")
        self.assertEqual(len(record.seq), 10026)
        self.assertEqual(segment, (1978, 606))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][5]
        (record, segment, status) = empty
        self.assertEqual(record.id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(record.seq), 359464)
        self.assertEqual(segment, (174889, 174752))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 4)
        self.assertEqual(len(alignment.annotations["empty"]), 6)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
            [[  3018161,   3018172,   3018174,   3018176,   3018176,   3018230],
             [157518434, 157518423, 157518423, 157518423, 157518423, 157518369],
             [155028381, 155028370, 155028370, 155028370, 155028370, 155028316],
             [158040613, 158040613, 158040613, 158040611, 158040606, 158040552],
            ]
                    # fmt: on
                ),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 4781)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3018230 : 3018230 + 129],
            "AGGACAAAATAATACAGAtttttttttttttttttttttGCAGTACTGGAAATGGAATGAATGTCCCTCACAATCACTATCAAGGTCCCTATCAAGGCAATCACTCTGTCACCGAGCTACAGCCCCAGC",
        )
        self.assertEqual(
            alignment[0],
            "AGGA-CAAAATAATACAGAtttttttttttttttttttttGCAGTACTGGAAATGGAATGAATGTCCCTCACAATCACTATCAAGGTCCCTATCAAGGCAATCACTCTGTCACCGAGCTA-CAGCCCCAGC",
        )
        self.assertEqual(alignment.sequences[1].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[1].seq), 170899992)
        self.assertEqual(
            alignment.sequences[1].seq[155028221 : 155028221 + 95],
            "CATATACCTAATAAATCTATATATGTGTAAATTCTTTTTCATAATGACTATCAGAACATTGGGAGCCAGGTTCTGATATATATTGATTGGATTAT",
        )
        self.assertEqual(
            alignment[1],
            "ATAATCCAATCAATATATAT---------------------CAGAACCTGGCTCCCAATG-----TTCTGATAGTCATTATGAA----------AAAGAATTTACACATATATAGATTTATTAGGTATATG",
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[2].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[2].seq), 173908612)
        self.assertEqual(
            alignment.sequences[2].seq[157518274 : 157518274 + 95],
            "catatacctaataaatctatatatgtgtaaatTCTTTTTCACAATGACTATCAGAACATTGGGAGCCAGGTTCTGATATATATTGATTGGATTAT",
        )
        self.assertEqual(
            alignment[2],
            "ATAATCCAATCAATATATAT---------------------CAGAACCTGGCTCCCAATG-----TTCTGATAGTCATTGTGAA----------AAAGAatttacacatatatagatttattaggtatatg",
        )
        self.assertEqual(
            alignment.sequences[2].annotations["quality"],
            "99999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[2].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[2].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["rightCount"], 0)
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ornAna1.chr2")
        self.assertEqual(len(record.seq), 54797317)
        self.assertEqual(segment, (40046122, 40040432))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "calJac1.Contig6394")
        self.assertEqual(len(record.seq), 133105)
        self.assertEqual(segment, (6649, 10128))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][2]
        (record, segment, status) = empty
        self.assertEqual(record.id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(record.seq), 498454)
        self.assertEqual(segment, (326690, 324316))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][3]
        (record, segment, status) = empty
        self.assertEqual(record.id, "echTel1.scaffold_288249")
        self.assertEqual(len(record.seq), 100002)
        self.assertEqual(segment, (87661, 95225))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][4]
        (record, segment, status) = empty
        self.assertEqual(record.id, "cavPor2.scaffold_216473")
        self.assertEqual(len(record.seq), 10026)
        self.assertEqual(segment, (1978, 606))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][5]
        (record, segment, status) = empty
        self.assertEqual(record.id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(record.seq), 359464)
        self.assertEqual(segment, (174889, 174752))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][6]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ponAbe2.chr6")
        self.assertEqual(len(record.seq), 174210431)
        self.assertEqual(segment, (158040552, 158040455))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 3)
        self.assertEqual(len(alignment.annotations["empty"]), 7)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
            [[  3018230,   3018234,   3018234,   3018249,   3018270,   3018289,
                3018294,   3018313,   3018323,   3018349,   3018349,   3018359],
             [155028316, 155028312, 155028311, 155028296, 155028296, 155028277,
              155028277, 155028258, 155028258, 155028232, 155028231, 155028221],
             [157518369, 157518365, 157518364, 157518349, 157518349, 157518330,
              157518330, 157518311, 157518311, 157518285, 157518284, 157518274],
            ]
                    # fmt: on
                ),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 61520)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3018359 : 3018359 + 123],
            "TTCAAACATGCATACATGCATTCATGTCTCATAATAATTATTAACATTGTCTTAGGCCAGAGGCTCGACTGCCCCAAAGCAATCCACTTAAACTGTCCCTGAGAAAGTCAttcctctccctaa",
        )
        self.assertEqual(
            alignment[0],
            "TT-CAAACATGCATACATGCATTCATGTCTCATAA-TAATTATTAACA-TTGTCTTAGGCCAGAGGCTCGACTGCCCCAAAGCAATCCACT-------TAAACTGTCCCTGAGAA-AGTCAttcctctccctaa",
        )
        self.assertEqual(alignment.sequences[1].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[1].seq), 173908612)
        self.assertEqual(
            alignment.sequences[1].seq[157518143 : 157518143 + 131],
            "AAAGAAAAGGATGACCCTTCCCGGAGAGAATATAAAATATGGATGGCTTAGTTTGTATCACTCAGGCCACTGACCTAAGATGACTGTTAATAATTATTTATGAGacatacacatatatttgtacatttata",
        )
        self.assertEqual(
            alignment[1],
            "-tataaatgtacaaatatatgtgtatgtCTCATAAATAATTATTAACAGTCATCTTAGGTCAGTGGCCTGAGTGATACAAACTAAGCCATCCATATTTTATATTCTCTCCGGGAAGGGTCATCCTTTTCTTT--",
        )
        self.assertEqual(
            alignment.sequences[1].annotations["quality"],
            "99999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 2400)
        self.assertEqual(alignment.sequences[2].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[2].seq), 170899992)
        self.assertEqual(
            alignment.sequences[2].seq[155028090 : 155028090 + 131],
            "AAAGAAAAGGATGACCCTTCCCGGAGAGAATATAAAATATGGATGGCTTAGTTTGTATCATTCAGGCCACTGACCTAAGATGACTGTTAATAATTATTTATGAGACATACACATATATTTGTACGTTTATA",
        )
        self.assertEqual(
            alignment[2],
            "-TATAAACGTACAAATATATGTGTATGTCTCATAAATAATTATTAACAGTCATCTTAGGTCAGTGGCCTGAATGATACAAACTAAGCCATCCATATTTTATATTCTCTCCGGGAAGGGTCATCCTTTTCTTT--",
        )
        self.assertEqual(alignment.sequences[2].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[2].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[2].annotations["rightCount"], 2402)
        self.assertEqual(alignment.sequences[3].id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(alignment.sequences[3].seq), 359464)
        self.assertEqual(
            alignment.sequences[3].seq[174625 : 174625 + 127],
            "TCAGAAGAAAGGATGACCCTTCACACAGGGAGGATAGGATGGGTTAGCTTGTGTCACTCGGGCCTCTAACCTAAGACAAATGTTAATAATTACTTACAGGACCTACATTTACACATACATGTATAAA",
        )
        self.assertEqual(
            alignment[3],
            "TT-TATACATGTATGTGTAAATGTAGGTCCTGTAAGTAATTATTAACATTTGTCTTAGGTTAGAGGCCCGAGTGACACAAGCTAACCCATCC------TATCCTCCCTGTGTGAAGGGTCATCCTTTCTTCTGA",
        )
        self.assertEqual(
            alignment.sequences[3].annotations["quality"],
            "9999999999999999999999999999999999999999999999999999667736999999999999999999999999999995666677755798899998999967967999999999589",
        )
        self.assertEqual(alignment.sequences[3].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[3].annotations["leftCount"], 137)
        self.assertEqual(alignment.sequences[3].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[3].annotations["rightCount"], 2731)
        self.assertEqual(alignment.sequences[4].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[4].seq), 174210431)
        self.assertEqual(
            alignment.sequences[4].seq[158040326 : 158040326 + 129],
            "AAAGAAAAGGATGACCCTTCCCAGAGAGAATATAAAATATGGATGGCTTAGTTTCTATCACTCAGGCCACTGACCTAAGATGACTGTTAATAATtatttatgacacatacacatatatttgtacattta",
        )
        self.assertEqual(
            alignment[4],
            "---taaatgtacaaatatatgtgtatgtgtcataaataATTATTAACAGTCATCTTAGGTCAGTGGCCTGAGTGATAGAAACTAAGCCATCCATATTTTATATTCTCTCTGGGAAGGGTCATCCTTTTCTTT--",
        )
        self.assertEqual(alignment.sequences[4].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[4].annotations["leftCount"], 97)
        self.assertEqual(alignment.sequences[4].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[4].annotations["rightCount"], 2523)
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ornAna1.chr2")
        self.assertEqual(len(record.seq), 54797317)
        self.assertEqual(segment, (40046122, 40040432))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "calJac1.Contig6394")
        self.assertEqual(len(record.seq), 133105)
        self.assertEqual(segment, (6649, 10128))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][2]
        (record, segment, status) = empty
        self.assertEqual(record.id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(record.seq), 498454)
        self.assertEqual(segment, (326690, 324316))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][3]
        (record, segment, status) = empty
        self.assertEqual(record.id, "echTel1.scaffold_288249")
        self.assertEqual(len(record.seq), 100002)
        self.assertEqual(segment, (87661, 95225))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][4]
        (record, segment, status) = empty
        self.assertEqual(record.id, "cavPor2.scaffold_216473")
        self.assertEqual(len(record.seq), 10026)
        self.assertEqual(segment, (1978, 606))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 5)
        self.assertEqual(len(alignment.annotations["empty"]), 5)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
            [[  3018359,   3018360,   3018361,   3018361,   3018393,   3018393,
                3018405,   3018405,   3018447,   3018447,   3018447,   3018464,
                3018464,   3018480,   3018482],
             [157518274, 157518274, 157518273, 157518272, 157518240, 157518239,
              157518227, 157518226, 157518184, 157518183, 157518177, 157518160,
              157518159, 157518143, 157518143],
             [155028221, 155028221, 155028220, 155028219, 155028187, 155028186,
              155028174, 155028173, 155028131, 155028130, 155028124, 155028107,
              155028106, 155028090, 155028090],
             [   174752,    174751,    174750,    174750,    174718,    174717,
                 174705,    174704,    174662,    174661,    174661,    174644,
                 174643,    174627,    174625],
             [158040455, 158040455, 158040455, 158040455, 158040423, 158040422,
              158040410, 158040409, 158040367, 158040366, 158040360, 158040343,
              158040342, 158040326, 158040326],
            ]
                    # fmt: on
                ),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3018482 : 3018482 + 162],
            "tcttcatctcctcttttcctccttttttttttctcatttctctttctctttcttttgtccttttccttTATAGCAAGCAAGGCAAGTAGTCTCTATTTAGAAGGCATggagagaatggggagaggaggaaaggaggagaggggaggagaggaggggagGTAT",
        )
        self.assertEqual(
            alignment[0],
            "tcttcatctcctcttttcctccttttttttttctcatttctctttctctttcttttgtccttttccttTATAGCAAGCAAGGCAAGTAGTCTCTATTTAGAAGGCATggagagaatggggagaggaggaaaggaggagaggggaggagaggaggggagGTAT",
        )
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ornAna1.chr2")
        self.assertEqual(len(record.seq), 54797317)
        self.assertEqual(segment, (40046122, 40040432))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "calJac1.Contig6394")
        self.assertEqual(len(record.seq), 133105)
        self.assertEqual(segment, (6649, 10128))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][2]
        (record, segment, status) = empty
        self.assertEqual(record.id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(record.seq), 498454)
        self.assertEqual(segment, (326690, 324316))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][3]
        (record, segment, status) = empty
        self.assertEqual(record.id, "echTel1.scaffold_288249")
        self.assertEqual(len(record.seq), 100002)
        self.assertEqual(segment, (87661, 95225))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][4]
        (record, segment, status) = empty
        self.assertEqual(record.id, "cavPor2.scaffold_216473")
        self.assertEqual(len(record.seq), 10026)
        self.assertEqual(segment, (1978, 606))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][5]
        (record, segment, status) = empty
        self.assertEqual(record.id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(record.seq), 359464)
        self.assertEqual(segment, (174625, 171894))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][6]
        (record, segment, status) = empty
        self.assertEqual(record.id, "hg18.chr6")
        self.assertEqual(len(record.seq), 170899992)
        self.assertEqual(segment, (155028090, 155025688))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][7]
        (record, segment, status) = empty
        self.assertEqual(record.id, "panTro2.chr6")
        self.assertEqual(len(record.seq), 173908612)
        self.assertEqual(segment, (157518143, 157515743))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][8]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ponAbe2.chr6")
        self.assertEqual(len(record.seq), 174210431)
        self.assertEqual(segment, (158040326, 158037803))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 1)
        self.assertEqual(len(alignment.annotations["empty"]), 9)
        self.assertTrue(
            numpy.array_equal(alignment.coordinates, numpy.array([[3018482, 3018644]]))
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1520)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3018644 : 3018644 + 178],
            "AGGGTAGGCCAAGTGCCTTGGGAAGTAGTTGTTGGTAGACTGAAAGTGTGTTCTGAGTGTCAGTGATGTTCATGAGATTATCACCAGCAAGGATGGCTGACGGGAACTGCAAGAGGCATAGCCCTGAGTTCTAAAGGAGAGGGAAACGTCACAGAAAGGATGCACTGTTTCAGCATCT",
        )
        self.assertEqual(
            alignment[0],
            "AGGGTAGGCCAAGTGCCTTGGGAAGTAGTTGTTGGTAGACTGAAAGTGTGTTC---TGAGTGTCAGTGATGTTCA-TGAGATTATCACCAGCAAGGATG--GCTGACGGGAACTG---CAAGAGGCATAGCCCTGAGTTCTAAAGGAGAGGGAAACGTCACAGAAAGGATG--------------------------CACTGTTTCAGCATCT",
        )
        self.assertEqual(alignment.sequences[1].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[1].seq), 125616256)
        self.assertEqual(
            alignment.sequences[1].seq[47545632 : 47545632 + 204],
            "AAATTCTTCAATAGGGAAATTTTATATTTCATATTACAATTCTATTCTATCCATGATCTTTCCTCTCATTTTAGAGTTCTGGCTGGTGGTATTTGCTTCTGATTCCATCAACCGCATCCTTAGGAGAGGATCTCAGTATCCATCACTTACACGCATGTGAATAAACTCAATGTGCTATCCCTCCCAGGGGCATTTGCCTTTTCT",
        )
        self.assertEqual(
            alignment[1],
            "AGAAAAGGCAAATGCCCCTGGGAGGGA-----TAGCACATTGA--GTTTATTCACATGCGTGTAAGTGATGGATACTGAGATCCTCTCC--TAAGGATGCGGTTGATGGAATCAGAAGCAAATACCACCAGCCAGAACTCTAAAATGAGAGGAAAGATCATGGATAGAATAGAATTGTAATATGAAATATAAAATTTCCCTATTGAAGAATTT",
        )
        self.assertEqual(
            alignment.sequences[1].annotations["quality"],
            "999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "N")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 160)
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ornAna1.chr2")
        self.assertEqual(len(record.seq), 54797317)
        self.assertEqual(segment, (40046122, 40040432))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "calJac1.Contig6394")
        self.assertEqual(len(record.seq), 133105)
        self.assertEqual(segment, (6649, 10128))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][2]
        (record, segment, status) = empty
        self.assertEqual(record.id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(record.seq), 498454)
        self.assertEqual(segment, (326690, 324316))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][3]
        (record, segment, status) = empty
        self.assertEqual(record.id, "echTel1.scaffold_288249")
        self.assertEqual(len(record.seq), 100002)
        self.assertEqual(segment, (87661, 95225))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][4]
        (record, segment, status) = empty
        self.assertEqual(record.id, "cavPor2.scaffold_216473")
        self.assertEqual(len(record.seq), 10026)
        self.assertEqual(segment, (1978, 606))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][5]
        (record, segment, status) = empty
        self.assertEqual(record.id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(record.seq), 359464)
        self.assertEqual(segment, (174625, 171894))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][6]
        (record, segment, status) = empty
        self.assertEqual(record.id, "hg18.chr6")
        self.assertEqual(len(record.seq), 170899992)
        self.assertEqual(segment, (155028090, 155025688))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][7]
        (record, segment, status) = empty
        self.assertEqual(record.id, "panTro2.chr6")
        self.assertEqual(len(record.seq), 173908612)
        self.assertEqual(segment, (157518143, 157515743))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][8]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ponAbe2.chr6")
        self.assertEqual(len(record.seq), 174210431)
        self.assertEqual(segment, (158040326, 158037803))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 2)
        self.assertEqual(len(alignment.annotations["empty"]), 9)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
            [[ 3018644,  3018671,  3018676,  3018687,  3018689,  3018697,
               3018697,  3018716,  3018716,  3018729,  3018731,  3018739,
               3018739,  3018753,  3018753,  3018806,  3018806,  3018822],
             [47545836, 47545809, 47545809, 47545798, 47545798, 47545790,
              47545787, 47545768, 47545767, 47545754, 47545754, 47545746,
              47545744, 47545730, 47545727, 47545674, 47545648, 47545632],
            ]
                    # fmt: on
                ),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1986)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3018822 : 3018822 + 110],
            "CTGCCTTCCATTACGATTTACTGATCACTTACAACCCTCCCACAGAAGAGAACCTAACTTGCTTAGGAGCATATGTACAGTTAATCAAGACAAAAATAAGAATGGAGACt",
        )
        self.assertEqual(
            alignment[0],
            "CTGCCTTCCATTACGATTTACTGATCACTTACAACCCTCCCACA----GAAGAGAACCTAACTTG-CTTAGGAGCATATGTACAGTTAATCAAGAC-----AAAAATAAGAATGGAGACt",
        )
        self.assertEqual(alignment.sequences[1].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[1].seq), 125616256)
        self.assertEqual(
            alignment.sequences[1].seq[47545353 : 47545353 + 119],
            "ATTTTCCATACTTATTTTTAGTTGATCTTGCTAAATTGTGTGCATCTCTAGAAGTCGATTTGAATTCCCTTTACAATGTGGAGAATTTTAAGTGAGTAATAATGGGCAGGGAAAAAAAG",
        )
        self.assertEqual(
            alignment[1],
            "CTTTTTTTCCCTGCCCATTATTACTCACTTAAAATTCTCC-ACATTGTAAAGGGAATTCAAATCGACTTCTAGAGATGCACACAATTTAGCAAGATCAACTAAAAATAAGTATGGAAAAT",
        )
        self.assertEqual(
            alignment.sequences[1].annotations["quality"],
            "99999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 160)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 0)
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ornAna1.chr2")
        self.assertEqual(len(record.seq), 54797317)
        self.assertEqual(segment, (40046122, 40040432))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "calJac1.Contig6394")
        self.assertEqual(len(record.seq), 133105)
        self.assertEqual(segment, (6649, 10128))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][2]
        (record, segment, status) = empty
        self.assertEqual(record.id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(record.seq), 498454)
        self.assertEqual(segment, (326690, 324316))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][3]
        (record, segment, status) = empty
        self.assertEqual(record.id, "echTel1.scaffold_288249")
        self.assertEqual(len(record.seq), 100002)
        self.assertEqual(segment, (87661, 95225))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][4]
        (record, segment, status) = empty
        self.assertEqual(record.id, "cavPor2.scaffold_216473")
        self.assertEqual(len(record.seq), 10026)
        self.assertEqual(segment, (1978, 606))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][5]
        (record, segment, status) = empty
        self.assertEqual(record.id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(record.seq), 359464)
        self.assertEqual(segment, (174625, 171894))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][6]
        (record, segment, status) = empty
        self.assertEqual(record.id, "hg18.chr6")
        self.assertEqual(len(record.seq), 170899992)
        self.assertEqual(segment, (155028090, 155025688))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][7]
        (record, segment, status) = empty
        self.assertEqual(record.id, "panTro2.chr6")
        self.assertEqual(len(record.seq), 173908612)
        self.assertEqual(segment, (157518143, 157515743))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][8]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ponAbe2.chr6")
        self.assertEqual(len(record.seq), 174210431)
        self.assertEqual(segment, (158040326, 158037803))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 2)
        self.assertEqual(len(alignment.annotations["empty"]), 9)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
            [[ 3018822,  3018862,  3018863,  3018866,  3018866,  3018883,
               3018883,  3018913,  3018913,  3018932],
             [47545472, 47545432, 47545432, 47545429, 47545425, 47545408,
              47545407, 47545377, 47545372, 47545353],
            ]
                    # fmt: on
                ),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3018932 : 3018932 + 339],
            "gtctgagttagggttttactgctgtgaacagacaccatgaccaaggcatgtcttataaaaaaaatttaattagggctggcttacagattcagaggttcagtgggagcatcaaggtgggggcatggcagcatccaggcaggcatggtgcaggcagagctgagagttctacatcttcatccaaaggcttctagtggaagactgacttccaggcacctagggtgagggtcttaagcccacacccacagtgacacacctattccaaccaggtcacacctattccaacaaggccatacctccaaatggcaccactcctggtccaagaatatacaaaccatgaca",
        )
        self.assertEqual(
            alignment[0],
            "gtctgagttagggttttactgctgtgaacagacaccatgaccaaggcatgtcttataaaaaaaatttaattagggctggcttacagattcagaggttcagtgggagcatcaaggtgggggcatggcagcatccaggcaggcatggtgcaggcagagctgagagttctacatcttcatccaaaggcttctagtggaagactgacttccaggcacctagggtgagggtcttaagcccacacccacagtgacacacctattccaaccaggtcacacctattccaacaaggccatacctccaaatggcaccactcctggtccaagaatatacaaaccatgaca",
        )
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "canFam2.chr1")
        self.assertEqual(len(record.seq), 125616256)
        self.assertEqual(segment, (47545353, 47545353))
        self.assertEqual(status, "C")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ornAna1.chr2")
        self.assertEqual(len(record.seq), 54797317)
        self.assertEqual(segment, (40046122, 40040432))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][2]
        (record, segment, status) = empty
        self.assertEqual(record.id, "calJac1.Contig6394")
        self.assertEqual(len(record.seq), 133105)
        self.assertEqual(segment, (6649, 10128))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][3]
        (record, segment, status) = empty
        self.assertEqual(record.id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(record.seq), 498454)
        self.assertEqual(segment, (326690, 324316))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][4]
        (record, segment, status) = empty
        self.assertEqual(record.id, "echTel1.scaffold_288249")
        self.assertEqual(len(record.seq), 100002)
        self.assertEqual(segment, (87661, 95225))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][5]
        (record, segment, status) = empty
        self.assertEqual(record.id, "cavPor2.scaffold_216473")
        self.assertEqual(len(record.seq), 10026)
        self.assertEqual(segment, (1978, 606))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][6]
        (record, segment, status) = empty
        self.assertEqual(record.id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(record.seq), 359464)
        self.assertEqual(segment, (174625, 171894))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][7]
        (record, segment, status) = empty
        self.assertEqual(record.id, "hg18.chr6")
        self.assertEqual(len(record.seq), 170899992)
        self.assertEqual(segment, (155028090, 155025688))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][8]
        (record, segment, status) = empty
        self.assertEqual(record.id, "panTro2.chr6")
        self.assertEqual(len(record.seq), 173908612)
        self.assertEqual(segment, (157518143, 157515743))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][9]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ponAbe2.chr6")
        self.assertEqual(len(record.seq), 174210431)
        self.assertEqual(segment, (158040326, 158037803))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 1)
        self.assertEqual(len(alignment.annotations["empty"]), 10)
        self.assertTrue(
            numpy.array_equal(alignment.coordinates, numpy.array([[3018932, 3019271]]))
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 228)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3019271 : 3019271 + 106],
            "GAGACCAAATGTGGCGCTCACGTGAGGCCAGGAGTAAATCGCACACACAGCCCATGCTTTCACCATCTGCTAGGGTGCTCTGGAGCAGGGCAGGCTTCTAACCTGG",
        )
        self.assertEqual(
            alignment[0],
            "GAGACCAAATG------------TGGCGCTCACG-TGAGGCCAGGAGTAAATCGCACACACAGCCCATGCTTTCACCATCTGCTAGGGTGCTCTGGAGCAGGGCAGGCTTCTAACCTGG",
        )
        self.assertEqual(alignment.sequences[1].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[1].seq), 125616256)
        self.assertEqual(
            alignment.sequences[1].seq[47545247 : 47545247 + 106],
            "CCAGGTTTGGGACCTATGTGTGTCAAGTGAATAGAAGCAGGTGGGACAGCATGTGCGTCACACATGGACACAGCCTGGACATCGTGGGGCATTTTCCATTTGCCTC",
        )
        self.assertEqual(
            alignment[1],
            "GAGGC-AAATGGAAAATGCCCCACGATGTCCAGGCTGTGTCCATGTGTGA--CGCACATGCTGTCC----------CACCTGCTTCTATTCACTTGACACACATAGGTCCCAAACCTGG",
        )
        self.assertEqual(
            alignment.sequences[1].annotations["quality"],
            "9999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 0)
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ornAna1.chr2")
        self.assertEqual(len(record.seq), 54797317)
        self.assertEqual(segment, (40046122, 40040432))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "calJac1.Contig6394")
        self.assertEqual(len(record.seq), 133105)
        self.assertEqual(segment, (6649, 10128))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][2]
        (record, segment, status) = empty
        self.assertEqual(record.id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(record.seq), 498454)
        self.assertEqual(segment, (326690, 324316))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][3]
        (record, segment, status) = empty
        self.assertEqual(record.id, "echTel1.scaffold_288249")
        self.assertEqual(len(record.seq), 100002)
        self.assertEqual(segment, (87661, 95225))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][4]
        (record, segment, status) = empty
        self.assertEqual(record.id, "cavPor2.scaffold_216473")
        self.assertEqual(len(record.seq), 10026)
        self.assertEqual(segment, (1978, 606))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][5]
        (record, segment, status) = empty
        self.assertEqual(record.id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(record.seq), 359464)
        self.assertEqual(segment, (174625, 171894))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][6]
        (record, segment, status) = empty
        self.assertEqual(record.id, "hg18.chr6")
        self.assertEqual(len(record.seq), 170899992)
        self.assertEqual(segment, (155028090, 155025688))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][7]
        (record, segment, status) = empty
        self.assertEqual(record.id, "panTro2.chr6")
        self.assertEqual(len(record.seq), 173908612)
        self.assertEqual(segment, (157518143, 157515743))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][8]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ponAbe2.chr6")
        self.assertEqual(len(record.seq), 174210431)
        self.assertEqual(segment, (158040326, 158037803))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 2)
        self.assertEqual(len(alignment.annotations["empty"]), 9)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
            [[ 3019271,  3019276,  3019277,  3019282,  3019282,  3019293,
               3019293,  3019308,  3019310,  3019324,  3019334,  3019377],
             [47545353, 47545348, 47545348, 47545343, 47545331, 47545320,
              47545319, 47545304, 47545304, 47545290, 47545290, 47545247],
            ]
                    # fmt: on
                ),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 10938)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3019377 : 3019377 + 88],
            "CCCCAGCATTCTGGCAGACACAGTGAAAAGAGACAGATGGTCACTAATAAAATCTGTATAAATTAGATCTCAGAGGATGGATGGACCA",
        )
        self.assertEqual(
            alignment[0],
            "CCCCAGCATTCTGGCAGACACAGTG-AAAAGAGACAGATGGTCACTAATAAAATCTGT-ATAAATTAG-ATCTCAGAGGATGGATGGACCA",
        )
        self.assertEqual(alignment.sequences[1].id, "felCat3.scaffold_205680")
        self.assertEqual(len(alignment.sequences[1].seq), 119354)
        self.assertEqual(
            alignment.sequences[1].seq[46757 : 46757 + 88],
            "GAGTCCCATCCTCTGAGATTCTAGTTCACCACAAAGCTTACTGGTGGGCACATGCTTCTTTTTCACATTAGCTATCAGAACACTTGGG",
        )
        self.assertEqual(
            alignment[1],
            "CCCAAGTGTTCTGATAGCTAATGTGAAAAAGAAGCATGTGCCCACCAGTAAGCTTTGTGGTGAACTAGAATCTCAGAGGATG---GGACTC",
        )
        self.assertEqual(
            alignment.sequences[1].annotations["quality"],
            "9999999999999999999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "N")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[2].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[2].seq), 125616256)
        self.assertEqual(
            alignment.sequences[2].seq[47545159 : 47545159 + 88],
            "gagtcccatcctctgagattctaggtcattgcaaatctTATTAGCGGGCCCATGTTTCTTTTTCACAGAGGCAATCAGAACACTTGGG",
        )
        self.assertEqual(
            alignment[2],
            "CCCAAGTGTTCTGATTGCCTCTGTGAAAAAGAAACATGGGCCCGCTAATAagatttgcaatgacctagaatctcagaggatg---ggactc",
        )
        self.assertEqual(
            alignment.sequences[2].annotations["quality"],
            "9999999999999999999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[2].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[2].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["rightCount"], 0)
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ornAna1.chr2")
        self.assertEqual(len(record.seq), 54797317)
        self.assertEqual(segment, (40046122, 40040432))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "calJac1.Contig6394")
        self.assertEqual(len(record.seq), 133105)
        self.assertEqual(segment, (6649, 10128))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][2]
        (record, segment, status) = empty
        self.assertEqual(record.id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(record.seq), 498454)
        self.assertEqual(segment, (326690, 324316))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][3]
        (record, segment, status) = empty
        self.assertEqual(record.id, "echTel1.scaffold_288249")
        self.assertEqual(len(record.seq), 100002)
        self.assertEqual(segment, (87661, 95225))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][4]
        (record, segment, status) = empty
        self.assertEqual(record.id, "cavPor2.scaffold_216473")
        self.assertEqual(len(record.seq), 10026)
        self.assertEqual(segment, (1978, 606))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][5]
        (record, segment, status) = empty
        self.assertEqual(record.id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(record.seq), 359464)
        self.assertEqual(segment, (174625, 171894))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][6]
        (record, segment, status) = empty
        self.assertEqual(record.id, "hg18.chr6")
        self.assertEqual(len(record.seq), 170899992)
        self.assertEqual(segment, (155028090, 155025688))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][7]
        (record, segment, status) = empty
        self.assertEqual(record.id, "panTro2.chr6")
        self.assertEqual(len(record.seq), 173908612)
        self.assertEqual(segment, (157518143, 157515743))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][8]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ponAbe2.chr6")
        self.assertEqual(len(record.seq), 174210431)
        self.assertEqual(segment, (158040326, 158037803))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 3)
        self.assertEqual(len(alignment.annotations["empty"]), 9)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
            [[ 3019377,  3019402,  3019402,  3019434,  3019434,  3019443,
               3019443,  3019456,  3019459,  3019465],
             [   46845,    46820,    46819,    46787,    46786,    46777,
                 46776,    46763,    46763,    46757],
             [47545247, 47545222, 47545221, 47545189, 47545188, 47545179,
              47545178, 47545165, 47545165, 47545159],
            ]
                    # fmt: on
                ),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 36924)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3019465 : 3019465 + 47],
            "AAGATAGATATTTAGAAGTAGCTTTTTATGTTTTTCTGATGTGTGTT",
        )
        self.assertEqual(
            alignment[0], "AAGATAGATATTTAGAAGTAGCTTTTTATGTTTTTCTGATGTGTGTT"
        )
        self.assertEqual(alignment.sequences[1].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[1].seq), 133105)
        self.assertEqual(
            alignment.sequences[1].seq[10128 : 10128 + 47],
            "aacaTCTATATTTTGAAATGGCTTTTCATGTTACTCTGATGTGTGTC",
        )
        self.assertEqual(
            alignment[1], "aacaTCTATATTTTGAAATGGCTTTTCATGTTACTCTGATGTGTGTC"
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 3479)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[2].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[2].seq), 173908612)
        self.assertEqual(
            alignment.sequences[2].seq[157515703 : 157515703 + 40],
            "AAAACAcatcagaataacatgaaaagccatttcaaaatat",
        )
        self.assertEqual(
            alignment[2], "-------atattttgaaatggcttttcatgttattctgatgTGTTTT"
        )
        self.assertEqual(
            alignment.sequences[2].annotations["quality"],
            "9999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[2].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[2].annotations["leftCount"], 2400)
        self.assertEqual(alignment.sequences[2].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[3].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[3].seq), 170899992)
        self.assertEqual(
            alignment.sequences[3].seq[155025648 : 155025648 + 40],
            "AAAACAcatcagaataacatgaaaagccatttcaaaatat",
        )
        self.assertEqual(
            alignment[3], "-------atattttgaaatggcttttcatgttattctgatgTGTTTT"
        )
        self.assertEqual(alignment.sequences[3].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[3].annotations["leftCount"], 2402)
        self.assertEqual(alignment.sequences[3].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[3].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[4].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[4].seq), 125616256)
        self.assertEqual(
            alignment.sequences[4].seq[47545115 : 47545115 + 44],
            "aacacacattagaacaacatgaagagctatttttaaatatacct",
        )
        self.assertEqual(
            alignment[4], "---aggtatatttaaaaatagctcttcatgttgttctaatgtgtgtt"
        )
        self.assertEqual(
            alignment.sequences[4].annotations["quality"],
            "99999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[4].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[4].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[4].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[4].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[5].id, "felCat3.scaffold_205680")
        self.assertEqual(len(alignment.sequences[5].seq), 119354)
        self.assertEqual(
            alignment.sequences[5].seq[46714 : 46714 + 43],
            "AACACACATCAGAACAACATGAAAACTATTTTTAAAAACACGT",
        )
        self.assertEqual(
            alignment[5], "---ACGTGTTTTTAAAAATAG-TTTTCATGTTGTTCTGATGTGTGTT"
        )
        self.assertEqual(
            alignment.sequences[5].annotations["quality"],
            "9999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[5].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[5].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[5].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[5].annotations["rightCount"], 193)
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ornAna1.chr2")
        self.assertEqual(len(record.seq), 54797317)
        self.assertEqual(segment, (40046122, 40040432))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(record.seq), 498454)
        self.assertEqual(segment, (326690, 324316))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][2]
        (record, segment, status) = empty
        self.assertEqual(record.id, "echTel1.scaffold_288249")
        self.assertEqual(len(record.seq), 100002)
        self.assertEqual(segment, (87661, 95225))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][3]
        (record, segment, status) = empty
        self.assertEqual(record.id, "cavPor2.scaffold_216473")
        self.assertEqual(len(record.seq), 10026)
        self.assertEqual(segment, (1978, 606))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][4]
        (record, segment, status) = empty
        self.assertEqual(record.id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(record.seq), 359464)
        self.assertEqual(segment, (174625, 171894))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][5]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ponAbe2.chr6")
        self.assertEqual(len(record.seq), 174210431)
        self.assertEqual(segment, (158040326, 158037803))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 6)
        self.assertEqual(len(alignment.annotations["empty"]), 6)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
            [[  3019465,   3019468,   3019472,   3019486,   3019487,   3019512],
             [    10128,     10131,     10135,     10149,     10150,     10175],
             [157515743, 157515743, 157515743, 157515729, 157515728, 157515703],
             [155025688, 155025688, 155025688, 155025674, 155025673, 155025648],
             [ 47545159,  47545159,  47545155,  47545141,  47545140,  47545115],
             [    46757,     46757,     46753,     46739,     46739,     46714],
            ]
                    # fmt: on
                ),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 20303)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3019512 : 3019512 + 92],
            "TGCATCATTAAGACTAGAGTTCCTTTCTGTCTTTGCTTTCTTGACAGGGCCATGCTCGGCAGTCATTCTTAGACTGCTTTTTGTTTgtttgg",
        )
        self.assertEqual(
            alignment[0],
            "TGCATCATTAAGACTAGAGTTCCT---------------TTCTGTCTT---TGCTTTCTTG--------ACAGGGCCATGCTCGGCAGTCATTCTTAGACTGCTTTTTGTTTgtttgg",
        )
        self.assertEqual(alignment.sequences[1].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[1].seq), 133105)
        self.assertEqual(
            alignment.sequences[1].seq[10175 : 10175 + 86],
            "CCTCATGAAGGCCCAAGTTCCTAAAACATTAATTCTCTTCCTATTTCCTAGCTGTCTCGCCTAGGCCTTGCCCGCCAGCAATTCCC",
        )
        self.assertEqual(
            alignment[1],
            "--CCTCATGAAGGCCCAAGTTCCTAAA-----ACATTAATTCTCTTCC---TATTTCCTAGCTGTCTCGCCTAGGCCTTGCCCGCCAGCAATTCCC----------------------",
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[2].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[2].seq), 173908612)
        self.assertEqual(
            alignment.sequences[2].seq[157515610 : 157515610 + 93],
            "CCAATACCCATGTGGGAATTGCTGGCGGGCAAGACCTCTTTTCTAGAAAAGAGGAAGAGAATGAATGTTTTAGGAACTTGGGCCTTAATGAAG",
        )
        self.assertEqual(
            alignment[2],
            "--CTTCATTAAGGCCCAAGTTCCTAAA-----ACATTCATTCTCTTCC---TCTTTTCTAG------AAAAGAGGTCTTGCCCGCCAGCAATTCCCACATGGGTATTGG---------",
        )
        self.assertEqual(
            alignment.sequences[2].annotations["quality"],
            "999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[2].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[2].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[3].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[3].seq), 170899992)
        self.assertEqual(
            alignment.sequences[3].seq[155025555 : 155025555 + 93],
            "CCAATACCCATGTGGGAATTGCTGGCGGGCAAGACCTCTTTTCTAGAAAAGAGGAAGAGAATGAATGTTTTAGGAACTTGGGCCTTAATGAAG",
        )
        self.assertEqual(
            alignment[3],
            "--CTTCATTAAGGCCCAAGTTCCTAAA-----ACATTCATTCTCTTCC---TCTTTTCTAG------AAAAGAGGTCTTGCCCGCCAGCAATTCCCACATGGGTATTGG---------",
        )
        self.assertEqual(alignment.sequences[3].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[3].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[3].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[3].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[4].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[4].seq), 125616256)
        self.assertEqual(
            alignment.sequences[4].seq[47545018 : 47545018 + 97],
            "CCGAACAAAGTAGGAATGACTGGCAGGCAAGACACTGATGAGGAAGGGTTGAAGAAAGAAAAAATATAAGCCTGTAAGAACCTGGGTTTAGATGAGG",
        )
        self.assertEqual(
            alignment[4],
            "--CCTCATCTAAACCCAGGTTCTTACAGGCTTATATTTTTTCTTTCTTCAACCCTTCCTCA--------TCAGTGTCTTGCCTGCCAGTCATTCCTAC-----------TTTGTTCGG",
        )
        self.assertEqual(
            alignment.sequences[4].annotations["quality"],
            "9999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[4].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[4].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[4].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[4].annotations["rightCount"], 0)
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "felCat3.scaffold_205680")
        self.assertEqual(len(record.seq), 119354)
        self.assertEqual(segment, (46714, 46521))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ornAna1.chr2")
        self.assertEqual(len(record.seq), 54797317)
        self.assertEqual(segment, (40046122, 40040432))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][2]
        (record, segment, status) = empty
        self.assertEqual(record.id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(record.seq), 498454)
        self.assertEqual(segment, (326690, 324316))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][3]
        (record, segment, status) = empty
        self.assertEqual(record.id, "echTel1.scaffold_288249")
        self.assertEqual(len(record.seq), 100002)
        self.assertEqual(segment, (87661, 95225))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][4]
        (record, segment, status) = empty
        self.assertEqual(record.id, "cavPor2.scaffold_216473")
        self.assertEqual(len(record.seq), 10026)
        self.assertEqual(segment, (1978, 606))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][5]
        (record, segment, status) = empty
        self.assertEqual(record.id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(record.seq), 359464)
        self.assertEqual(segment, (174625, 171894))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][6]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ponAbe2.chr6")
        self.assertEqual(len(record.seq), 174210431)
        self.assertEqual(segment, (158040326, 158037803))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 5)
        self.assertEqual(len(alignment.annotations["empty"]), 7)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
            [[  3019512,   3019514,   3019536,   3019536,   3019536,   3019536,
                3019545,   3019545,   3019555,   3019555,   3019555,   3019582,
                3019584,   3019595,   3019604],
             [    10175,     10175,     10197,     10200,     10200,     10207,
                  10216,     10216,     10226,     10232,     10234,     10261,
                  10261,     10261,     10261],
             [157515703, 157515703, 157515681, 157515678, 157515678, 157515671,
              157515662, 157515662, 157515652, 157515652, 157515650, 157515623,
              157515621, 157515610, 157515610],
             [155025648, 155025648, 155025626, 155025623, 155025623, 155025616,
              155025607, 155025607, 155025597, 155025597, 155025595, 155025568,
              155025566, 155025555, 155025555],
             [ 47545115,  47545115,  47545093,  47545090,  47545085,  47545078,
               47545069,  47545066,  47545056,  47545056,  47545056,  47545029,
               47545027,  47545027,  47545018],
            ]
                    # fmt: on
                ),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3019604 : 3019604 + 98],
            "tttggtttggtttggttttttcaagacagggtttctttgtatagtcctagctgtcctggaactcactttgtagaccagactggccttgaactcagaaa",
        )
        self.assertEqual(
            alignment[0],
            "tttggtttggtttggttttttcaagacagggtttctttgtatagtcctagctgtcctggaactcactttgtagaccagactggccttgaactcagaaa",
        )
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "felCat3.scaffold_205680")
        self.assertEqual(len(record.seq), 119354)
        self.assertEqual(segment, (46714, 46521))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "canFam2.chr1")
        self.assertEqual(len(record.seq), 125616256)
        self.assertEqual(segment, (47545018, 47545018))
        self.assertEqual(status, "C")
        empty = alignment.annotations["empty"][2]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ornAna1.chr2")
        self.assertEqual(len(record.seq), 54797317)
        self.assertEqual(segment, (40046122, 40040432))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][3]
        (record, segment, status) = empty
        self.assertEqual(record.id, "calJac1.Contig6394")
        self.assertEqual(len(record.seq), 133105)
        self.assertEqual(segment, (10261, 10261))
        self.assertEqual(status, "C")
        empty = alignment.annotations["empty"][4]
        (record, segment, status) = empty
        self.assertEqual(record.id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(record.seq), 498454)
        self.assertEqual(segment, (326690, 324316))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][5]
        (record, segment, status) = empty
        self.assertEqual(record.id, "echTel1.scaffold_288249")
        self.assertEqual(len(record.seq), 100002)
        self.assertEqual(segment, (87661, 95225))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][6]
        (record, segment, status) = empty
        self.assertEqual(record.id, "cavPor2.scaffold_216473")
        self.assertEqual(len(record.seq), 10026)
        self.assertEqual(segment, (1978, 606))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][7]
        (record, segment, status) = empty
        self.assertEqual(record.id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(record.seq), 359464)
        self.assertEqual(segment, (174625, 171894))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][8]
        (record, segment, status) = empty
        self.assertEqual(record.id, "hg18.chr6")
        self.assertEqual(len(record.seq), 170899992)
        self.assertEqual(segment, (155025555, 155025555))
        self.assertEqual(status, "C")
        empty = alignment.annotations["empty"][9]
        (record, segment, status) = empty
        self.assertEqual(record.id, "panTro2.chr6")
        self.assertEqual(len(record.seq), 173908612)
        self.assertEqual(segment, (157515610, 157515610))
        self.assertEqual(status, "C")
        empty = alignment.annotations["empty"][10]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ponAbe2.chr6")
        self.assertEqual(len(record.seq), 174210431)
        self.assertEqual(segment, (158040326, 158037803))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 1)
        self.assertEqual(len(alignment.annotations["empty"]), 11)
        self.assertTrue(
            numpy.array_equal(alignment.coordinates, numpy.array([[3019604, 3019702]]))
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 45)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3019702 : 3019702 + 42],
            "tctgcctgcctctgcctcccaagtcctgggattaaaggcgtg",
        )
        self.assertEqual(
            alignment[0],
            "tctgcctgcctctgcctcccaag--------------------------------tcctgggattaaaggcgtg",
        )
        self.assertEqual(alignment.sequences[1].id, "felCat3.scaffold_205680")
        self.assertEqual(len(alignment.sequences[1].seq), 119354)
        self.assertEqual(
            alignment.sequences[1].seq[46447 : 46447 + 74],
            "TAAACCTTTAAGAACTTGAGTTTGTTTGTTTGGttttttttaatgtttatttttgagagagagagagagacaga",
        )
        self.assertEqual(
            alignment[1],
            "tctgtctctctctctctctcaaaaataaacattaaaaaaaaCCAAACAAACAAACTCAAGTTCTTAAAGGTTTA",
        )
        self.assertEqual(
            alignment.sequences[1].annotations["quality"],
            "99999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 193)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 0)
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "canFam2.chr1")
        self.assertEqual(len(record.seq), 125616256)
        self.assertEqual(segment, (47545018, 47545018))
        self.assertEqual(status, "C")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ornAna1.chr2")
        self.assertEqual(len(record.seq), 54797317)
        self.assertEqual(segment, (40046122, 40040432))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][2]
        (record, segment, status) = empty
        self.assertEqual(record.id, "calJac1.Contig6394")
        self.assertEqual(len(record.seq), 133105)
        self.assertEqual(segment, (10261, 10261))
        self.assertEqual(status, "C")
        empty = alignment.annotations["empty"][3]
        (record, segment, status) = empty
        self.assertEqual(record.id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(record.seq), 498454)
        self.assertEqual(segment, (326690, 324316))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][4]
        (record, segment, status) = empty
        self.assertEqual(record.id, "echTel1.scaffold_288249")
        self.assertEqual(len(record.seq), 100002)
        self.assertEqual(segment, (87661, 95225))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][5]
        (record, segment, status) = empty
        self.assertEqual(record.id, "cavPor2.scaffold_216473")
        self.assertEqual(len(record.seq), 10026)
        self.assertEqual(segment, (1978, 606))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][6]
        (record, segment, status) = empty
        self.assertEqual(record.id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(record.seq), 359464)
        self.assertEqual(segment, (174625, 171894))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][7]
        (record, segment, status) = empty
        self.assertEqual(record.id, "hg18.chr6")
        self.assertEqual(len(record.seq), 170899992)
        self.assertEqual(segment, (155025555, 155025555))
        self.assertEqual(status, "C")
        empty = alignment.annotations["empty"][8]
        (record, segment, status) = empty
        self.assertEqual(record.id, "panTro2.chr6")
        self.assertEqual(len(record.seq), 173908612)
        self.assertEqual(segment, (157515610, 157515610))
        self.assertEqual(status, "C")
        empty = alignment.annotations["empty"][9]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ponAbe2.chr6")
        self.assertEqual(len(record.seq), 174210431)
        self.assertEqual(segment, (158040326, 158037803))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 2)
        self.assertEqual(len(alignment.annotations["empty"]), 10)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
                numpy.array([[3019702, 3019725, 3019725, 3019744],
                             [  46521,   46498,   46466,   46447],
                            ])
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, -16865)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3019744 : 3019744 + 33],
            "cgccaccactgccctgcCTTAAACTGCTCTTAA",
        )
        self.assertEqual(
            alignment[0],
            "-----------------cgccaccactgccctgcCT------------TAAACTGCTCTTAA",
        )
        self.assertEqual(alignment.sequences[1].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[1].seq), 174210431)
        self.assertEqual(
            alignment.sequences[1].seq[158037776 : 158037776 + 27],
            "CCAATACCCATGTGGGAATTGCTGGCG",
        )
        self.assertEqual(
            alignment[1],
            "-----------------CGCCAGCAATTCC------------------CACATGGGTATTGG",
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 2523)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[2].id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(alignment.sequences[2].seq), 498454)
        self.assertEqual(alignment.sequences[2].seq[324312 : 324312 + 4], "TCAA")
        self.assertEqual(
            alignment[2],
            "----------------------------------------------------------TTGA",
        )
        self.assertEqual(alignment.sequences[2].annotations["quality"], "9999")
        self.assertEqual(alignment.sequences[2].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[2].annotations["leftCount"], 2374)
        self.assertEqual(alignment.sequences[2].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[3].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[3].seq), 133105)
        self.assertEqual(
            alignment.sequences[3].seq[10261 : 10261 + 13], "ACATGGCTACTGG"
        )
        self.assertEqual(
            alignment[3],
            "-------------------------------------------------ACATGGCTACTGG",
        )
        self.assertEqual(alignment.sequences[3].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[3].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[3].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[3].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[4].id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(alignment.sequences[4].seq), 359464)
        self.assertEqual(
            alignment.sequences[4].seq[171883 : 171883 + 11], "TCAACACCAAT"
        )
        self.assertEqual(
            alignment[4],
            "---------------------------------------------------ATTGGTGTTGA",
        )
        self.assertEqual(alignment.sequences[4].annotations["quality"], "87784564678")
        self.assertEqual(alignment.sequences[4].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[4].annotations["leftCount"], 2731)
        self.assertEqual(alignment.sequences[4].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[4].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[5].id, "oryCun1.scaffold_156751")
        self.assertEqual(len(alignment.sequences[5].seq), 4726)
        self.assertEqual(
            alignment.sequences[5].seq[4493 : 4493 + 41],
            "TCAACACCCCTATGGGACCTGAGACAGGCAAGACAATGGTG",
        )
        self.assertEqual(
            alignment[5],
            "--------------------CACCATTGTCTTGCCTGTC-TCAGGTCCCATAGGGGTGTTGA",
        )
        self.assertEqual(
            alignment.sequences[5].annotations["quality"],
            "99999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[5].annotations["leftStatus"], "N")
        self.assertEqual(alignment.sequences[5].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[5].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[5].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[6].id, "felCat3.scaffold_205680")
        self.assertEqual(len(alignment.sequences[6].seq), 119354)
        self.assertEqual(
            alignment.sequences[6].seq[46385 : 46385 + 62],
            "TCTACACCAATGTGGGAATGACTGAAAGACAAGACACTGATGAGTAAGAGAAAGAGAAATCA",
        )
        self.assertEqual(
            alignment[6],
            "TGATTTCTCTTTCTCTTACTCATCAGTGTCTTGTCTTTCAGTCATTCCCACATTGGTGTAGA",
        )
        self.assertEqual(
            alignment.sequences[6].annotations["quality"],
            "99999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[6].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[6].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[6].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[6].annotations["rightCount"], 0)
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "canFam2.chr1")
        self.assertEqual(len(record.seq), 125616256)
        self.assertEqual(segment, (47545018, 47545018))
        self.assertEqual(status, "C")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ornAna1.chr2")
        self.assertEqual(len(record.seq), 54797317)
        self.assertEqual(segment, (40046122, 40040432))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][2]
        (record, segment, status) = empty
        self.assertEqual(record.id, "echTel1.scaffold_288249")
        self.assertEqual(len(record.seq), 100002)
        self.assertEqual(segment, (87661, 95225))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][3]
        (record, segment, status) = empty
        self.assertEqual(record.id, "cavPor2.scaffold_216473")
        self.assertEqual(len(record.seq), 10026)
        self.assertEqual(segment, (1978, 606))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][4]
        (record, segment, status) = empty
        self.assertEqual(record.id, "hg18.chr6")
        self.assertEqual(len(record.seq), 170899992)
        self.assertEqual(segment, (155025555, 155025555))
        self.assertEqual(status, "C")
        empty = alignment.annotations["empty"][5]
        (record, segment, status) = empty
        self.assertEqual(record.id, "panTro2.chr6")
        self.assertEqual(len(record.seq), 173908612)
        self.assertEqual(segment, (157515610, 157515610))
        self.assertEqual(status, "C")
        self.assertEqual(len(alignment.sequences), 7)
        self.assertEqual(len(alignment.annotations["empty"]), 6)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
            [[  3019744,   3019744,   3019747,   3019757,   3019763,   3019763,
                3019763,   3019763,   3019764,   3019766,   3019773,   3019777],
             [158037803, 158037803, 158037800, 158037790, 158037790, 158037790,
              158037790, 158037790, 158037789, 158037787, 158037780, 158037776],
             [   324316,    324316,    324316,    324316,    324316,    324316,
                 324316,    324316,    324316,    324316,    324316,    324312],
             [    10261,     10261,     10261,     10261,     10261,     10261,
                  10261,     10261,     10261,     10263,     10270,     10274],
             [   171894,    171894,    171894,    171894,    171894,    171894,
                 171894,    171894,    171894,    171894,    171887,    171883],
             [     4534,      4534,      4534,      4524,      4518,      4515,
                   4515,      4507,      4506,      4504,      4497,      4493],
             [    46447,     46430,     46427,     46417,     46411,     46408,
                  46407,     46399,     46398,     46396,     46389,     46385],
            ]
                    # fmt: on
                ),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 367532)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3019777 : 3019777 + 183],
            "GGCAATGTCAGGCTATGCGTTCTAGACAGGGCACAAGAAAAGCTTTTAGCAGCAGAATAAACTTTTAAAGTAAATTACTTTCCTTGATAGCAACTAGACGACCCAATTGATACAGTGGAAAGAGGCCTTTGAGAATGCATGAGAGAATATTTCCTGTAAGAGTTGAACAATTTAGAATTTACc",
        )
        self.assertEqual(
            alignment[0],
            "GGCAATGTCA-GGCTATGCGT------TCTAGACAGGGCACAAGAAAAGCTTTTAGCAGCAGAATAAACTTTT-AAAGTAAATTACTTTCCTTGATAGCAACTAGACGACCCAATTGA-TACAGT------GGAAAG-----A----------GGCCTTTGAGAAT---GCATGAGAGAATAT---TTCCTGTA------AGAGTTGAACAATTTAGAATTTACc",
        )
        self.assertEqual(alignment.sequences[1].id, "oryCun1.scaffold_156751")
        self.assertEqual(len(alignment.sequences[1].seq), 4726)
        self.assertEqual(
            alignment.sequences[1].seq[4300 : 4300 + 193],
            "AATCCTAAACTATTTAGCTCTGAGAAATAGGGGAAGTAATATTATTCCATCCTGATTCCCAAAAGGCTTTTTTTTTCTTTGTTGTATCACTTTGAACACTTGATTACTAGCAAGGGAAGCAGTTTCCTCTAAAGGCTTATTCTACTGGAAAAAGCTTTTCCTGTCCCTTAGAACATAAGGTCATGATGTTGCC",
        )
        self.assertEqual(
            alignment[1],
            "GGCAACATCATGACCTTATGT------TCTAA----GGGACAGGAAAAGCTTTTTCCAGTAGAATAAGCCTTT-AGAGGAAACTGCTTCCCTTGCTAGTAATCAAGTGTTCAAAGTGA-TACAACAAAGAAAAAAAA-----A----------GCCTTTTGGGAAT-CAGGATGGAATAATATTACTTCCCCTATTTCTCAGAGCTAAATAGTTTAGGATT----",
        )
        self.assertEqual(
            alignment.sequences[1].annotations["quality"],
            "9999999999999999999999999999999999999899999999989999999999999999999999999999999999999999999999999999898999999989979999999999999997999999998978998999999999999999978999999979689999999999999999979",
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 37)
        self.assertEqual(alignment.sequences[2].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[2].seq), 170899992)
        self.assertEqual(
            alignment.sequences[2].seq[155025365 : 155025365 + 190],
            "AATCCTAAATTGTTTAACTCTTAGAAATAAGGGAAATAATATTATTCTATCCTGCATTCTCAAAAAAATTTTTCTTCATCTCACTTTGGTGACCTGATTGCTATCAAGGAAAGCAATTCCCTTTAACATTTTATTCCACTGCAAAAAGCTTTTCCTGTACTTTATCTTGAATATGTATCATGGTGTTGCC",
        )
        self.assertEqual(
            alignment[2],
            "GGCAACACCATGA-TACATAT------TCAAGATAAAGTACAGGAAAAGCTTTTTGCAGTGGAATAAAATGTT-AAAGGGAATTGCTTTCCTTGATAGCAATCAGGTCACCAAAGTGA-GATGA-------AGAAAA-----A----------TTTTTTTGAGAATGCAGGATAGAATAATATTATTTCCCTTATTTCTAAGAGTTAAACAATTTAGGATT----",
        )
        self.assertEqual(alignment.sequences[2].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[2].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[3].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[3].seq), 173908612)
        self.assertEqual(
            alignment.sequences[3].seq[157515419 : 157515419 + 191],
            "AATCCTAAATTGTTTAACTCTTAGAAATAAGGGAAATAATATTATTCTATCCTGCATTCTCAAAAAAAATTTTTCTTCATCTCACTTTGGTGACCTGATTGCTATCAAGGAAAGCAATTTCCTTTAACATTTTATTCCACTGCAAAAAGCTTTTCCTGTACTTTATCTTGAATATGTATCATGGTGTTGCC",
        )
        self.assertEqual(
            alignment[3],
            "GGCAACACCATGA-TACATAT------TCAAGATAAAGTACAGGAAAAGCTTTTTGCAGTGGAATAAAATGTT-AAAGGAAATTGCTTTCCTTGATAGCAATCAGGTCACCAAAGTGA-GATGA-------AGAAAA-----A---------TTTTTTTTGAGAATGCAGGATAGAATAATATTATTTCCCTTATTTCTAAGAGTTAAACAATTTAGGATT----",
        )
        self.assertEqual(
            alignment.sequences[3].annotations["quality"],
            "99999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[3].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[3].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[3].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[3].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[4].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[4].seq), 174210431)
        self.assertEqual(
            alignment.sequences[4].seq[158037586 : 158037586 + 190],
            "AATCCTAAATTGTTTAACTCTAGAAATAAGGGAAATAATATTATTCTGTGCTGCATTCTCAAAAAGAATTTTTCTTCATCTCACTTTGGTGACCTGATTGCTATCAAGGAAAGCAATTTCCTTTAACATTTTATTCCACCGCAAAAAGCTTTTCCTGTACTTTATCTTGAATATGTATCATGGTGTTGCC",
        )
        self.assertEqual(
            alignment[4],
            "GGCAACACCATGA-TACATAT------TCAAGATAAAGTACAGGAAAAGCTTTTTGCGGTGGAATAAAATGTT-AAAGGAAATTGCTTTCCTTGATAGCAATCAGGTCACCAAAGTGA-GATGA-------AGAAAA-----A---------TTCTTTTTGAGAATGCAGCACAGAATAATATTATTTCCCTTATTTCT-AGAGTTAAACAATTTAGGATT----",
        )
        self.assertEqual(alignment.sequences[4].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[4].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[4].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[4].annotations["rightCount"], 21)
        self.assertEqual(alignment.sequences[5].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[5].seq), 133105)
        self.assertEqual(
            alignment.sequences[5].seq[10274 : 10274 + 187],
            "GGCAACACCGTGATACATATCCGAGATAAAGTACAGGAAAAGCTTTTGCAGTGGACGAAAATGTTGAAGGAAATTGCTTTTCTTGATAGCAATCAGGTCACCAAAGTGAGATGAAGAAAGTTTTTTGAGAATGTGGGATAGAATAATATTGTTTTCGTTATTTCTAAGAGTTAAACAATTTAGGATT",
        )
        self.assertEqual(
            alignment[5],
            "GGCAACACCGTGA-TACATAT------CCGAGATAAAGTACAGGAAAAGC-TTTTGCAGTGGACGAAAATGTT-GAAGGAAATTGCTTTTCTTGATAGCAATCAGGTCACCAAAGTGA-GATGA-------AGAAA-----------------GTTTTTTGAGAATGTGGGATAGAATAATATTGTTTTCGTTATTTCTAAGAGTTAAACAATTTAGGATT----",
        )
        self.assertEqual(alignment.sequences[5].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[5].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[5].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[5].annotations["rightCount"], 21)
        self.assertEqual(alignment.sequences[6].id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(alignment.sequences[6].seq), 359464)
        self.assertEqual(
            alignment.sequences[6].seq[171686 : 171686 + 197],
            "AATCCTGAATTATTTAACTCTGAAATGGAGGGTGGATAAAATTATTCTTTCCAGCATTCTCAAAACATTTTTTCTTTATAATCACTTTGGTCACCTGATTGCTATCAAGAAAAGCAATTTCCTTTAACCTTTTATTCTACCGTAAAGGGCTTTTCCTGCACTTTGTCTTTATCTACAGCAAGGATCACGGCGTTGCC",
        )
        self.assertEqual(
            alignment[6],
            "GGCAACGCCGTGA-TCCTTGCTGTAGATAAAGACAAAGTGCAGGAAAAGCCCTTTACGGTAGAATAAAAGGTT-AAAGGAAATTGCTTTTCTTGATAGCAATCAGGTGACCAAAGTGATTATAA-------AGAAAA-----A----------ATGTTTTGAGAATGCTGGAAAGAATAATTTTATCCACCCTCCATTTCAGAGTTAAATAATTCAGGATT----",
        )
        self.assertEqual(
            alignment.sequences[6].annotations["quality"],
            "54667653455648997776896979788699699955776666776867789776453646557877552667697777697779896767576676677696675376877786786858478887756858859779666569356677477759968759378767656876769557597655596536566",
        )
        self.assertEqual(alignment.sequences[6].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[6].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[6].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[6].annotations["rightCount"], 358)
        self.assertEqual(alignment.sequences[7].id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(alignment.sequences[7].seq), 498454)
        self.assertEqual(
            alignment.sequences[7].seq[324116 : 324116 + 196],
            "AAGCCCAAGTTATTTAACTTAGACATAGGGAAAATAATATGCTTCTATACTGCATCCACAAAAGTTTGTTTTTGTTTTTTTTTTCAGTATGACTTTGGTCACCTGATTGCTGTCAAGGAAGGCAATTTCCTTTAAAAGTTTATCCTACTGCAAAAAGCTTTTCCTGTACTTTATCTAGAACATTCATGACATTGCC",
        )
        self.assertEqual(
            alignment[7],
            "GGCAATGTCATGA----ATGT------TCTAGATAAAGTACAGGAAAAGCTTTTTGCAGTAGGATAAACTTTT-AAAGGAAATTGCCTTCCTTGACAGCAATCAGGTGACCAAAGTCATACTGA-------AAAAAA-----AAACAAAAACAAACTTTTGTGGATGCAGTATAGAAGCATATTATTTTCCCTATGTCT--AAGTTAAATAACTTGGGCTT----",
        )
        self.assertEqual(
            alignment.sequences[7].annotations["quality"],
            "9999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[7].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[7].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[7].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[7].annotations["rightCount"], 623)
        self.assertEqual(alignment.sequences[8].id, "felCat3.scaffold_205680")
        self.assertEqual(len(alignment.sequences[8].seq), 119354)
        self.assertEqual(
            alignment.sequences[8].seq[46199 : 46199 + 186],
            "aattttaaattatttaaCTCTTAGATGTAGGGGAAATAATATTATTGTATCCTGTATTCTCAAAAAAAATTCCCTTTATATTACTTTAATCACCTAATGTCTCTTGAAGAAAGCAATTTCCTACAACATTTTATTCTACTCAAAATGCTTTTCCTGCACTTTATCTAGAATGTACTTTATGATGTC",
        )
        self.assertEqual(
            alignment[8],
            "GACATC---ATAA-AGTACAT------TCTAGATAAAGTGCAGGAAAAGCATTTTG-AGTAGAATAAAATGTT-GTAGGAAATTGCTTTCTTCAAGAGACATTAGGTGATTAAAGTAA-TATAA-------AGGGAA-----T----------TTTTTTTGAGAATACAGGATACAATAATATTATTTCCCCTACATCTAAGAGttaaataatttaaaatt----",
        )
        self.assertEqual(
            alignment.sequences[8].annotations["quality"],
            "999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[8].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[8].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[8].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[8].annotations["rightCount"], 27)
        self.assertEqual(alignment.sequences[9].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[9].seq), 125616256)
        self.assertEqual(
            alignment.sequences[9].seq[47544825 : 47544825 + 193],
            "attaaattttaaattatttaaCTCTTAGAAGTAAGGGAAACAATATTATTCTAACTTGGATTCACTGAACTTTGTTTTTTCCCTTTATATTCCTTCAATCACCTGATTGCTCTCAAAGAAAGTAATTTCCTATTAACATTTTATCCTGTTCAAAATCCTTTTGTTGTACTTTACCCGGAATTTACATCATGAT",
        )
        self.assertEqual(
            alignment[9],
            "------ATCATGA-TGTAAAT------TCCGGGTAAAGTACAACAAAAGGATTTTG-AACAGGATAAAATGTTAATAGGAAATTACTTTCTTTGAGAGCAATCAGGTGATTGAAGGAA-TATAA-------AGGGAAAAAACA----------AAGTTCAGTGAATCCAAGTTAGAATAATATTGTTTCCCTTACTTCTAAGAGttaaataatttaaaatttaat",
        )
        self.assertEqual(
            alignment.sequences[9].annotations["quality"],
            "9999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[9].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[9].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[9].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[9].annotations["rightCount"], 0)
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ornAna1.chr2")
        self.assertEqual(len(record.seq), 54797317)
        self.assertEqual(segment, (40046122, 40040432))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "echTel1.scaffold_288249")
        self.assertEqual(len(record.seq), 100002)
        self.assertEqual(segment, (87661, 95225))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][2]
        (record, segment, status) = empty
        self.assertEqual(record.id, "cavPor2.scaffold_216473")
        self.assertEqual(len(record.seq), 10026)
        self.assertEqual(segment, (1978, 606))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 10)
        self.assertEqual(len(alignment.annotations["empty"]), 3)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
            [[  3019777,   3019783,   3019786,   3019787,   3019787,   3019789,
                3019790,   3019793,   3019797,   3019797,   3019802,   3019806,
                3019820,   3019821,   3019826,   3019827,   3019843,   3019843,
                3019887,   3019887,   3019892,   3019893,   3019893,   3019898,
                3019899,   3019899,   3019900,   3019900,   3019900,   3019913,
                3019913,   3019913,   3019927,   3019927,   3019935,   3019935,
                3019935,   3019936,   3019956,   3019960],
             [     4493,      4487,      4484,      4483,      4482,      4480,
                   4479,      4476,      4472,      4472,      4467,      4467,
                   4453,      4452,      4447,      4446,      4430,      4430,
                   4386,      4386,      4381,      4380,      4374,      4369,
                   4368,      4368,      4367,      4367,      4367,      4354,
                   4354,      4352,      4338,      4335,      4327,      4322,
                   4321,      4320,      4300,      4300],
             [155025555, 155025549, 155025546, 155025545, 155025544, 155025542,
              155025542, 155025539, 155025535, 155025535, 155025530, 155025526,
              155025512, 155025511, 155025506, 155025505, 155025489, 155025489,
              155025445, 155025445, 155025440, 155025440, 155025440, 155025435,
              155025434, 155025434, 155025433, 155025433, 155025433, 155025420,
              155025419, 155025417, 155025403, 155025400, 155025392, 155025387,
              155025386, 155025385, 155025365, 155025365],
             [157515610, 157515604, 157515601, 157515600, 157515599, 157515597,
              157515597, 157515594, 157515590, 157515590, 157515585, 157515581,
              157515567, 157515566, 157515561, 157515560, 157515544, 157515544,
              157515500, 157515500, 157515495, 157515495, 157515495, 157515490,
              157515489, 157515489, 157515488, 157515488, 157515487, 157515474,
              157515473, 157515471, 157515457, 157515454, 157515446, 157515441,
              157515440, 157515439, 157515419, 157515419],
             [158037776, 158037770, 158037767, 158037766, 158037765, 158037763,
              158037763, 158037760, 158037756, 158037756, 158037751, 158037747,
              158037733, 158037732, 158037727, 158037726, 158037710, 158037710,
              158037666, 158037666, 158037661, 158037661, 158037661, 158037656,
              158037655, 158037655, 158037654, 158037654, 158037653, 158037640,
              158037639, 158037637, 158037623, 158037620, 158037612, 158037607,
              158037607, 158037606, 158037586, 158037586],
             [    10274,     10280,     10283,     10284,     10285,     10287,
                  10287,     10290,     10294,     10294,     10299,     10303,
                  10317,     10317,     10322,     10323,     10339,     10339,
                  10383,     10383,     10388,     10388,     10388,     10393,
                  10393,     10393,     10393,     10393,     10393,     10406,
                  10407,     10409,     10423,     10426,     10434,     10439,
                  10440,     10441,     10461,     10461],
             [   171883,    171877,    171874,    171873,    171872,    171870,
                 171870,    171867,    171863,    171857,    171852,    171848,
                 171834,    171833,    171828,    171827,    171811,    171811,
                 171767,    171766,    171761,    171761,    171761,    171756,
                 171755,    171755,    171754,    171754,    171754,    171741,
                 171740,    171738,    171724,    171721,    171713,    171708,
                 171707,    171706,    171686,    171686],
             [   324312,    324306,    324303,    324302,    324301,    324299,
                 324299,    324299,    324295,    324295,    324290,    324286,
                 324272,    324271,    324266,    324265,    324249,    324249,
                 324205,    324204,    324199,    324199,    324199,    324194,
                 324193,    324193,    324192,    324183,    324182,    324169,
                 324168,    324166,    324152,    324149,    324141,    324136,
                 324136,    324136,    324116,    324116],
             [    46385,     46379,     46379,     46378,     46377,     46375,
                  46375,     46372,     46368,     46368,     46363,     46359,
                  46345,     46344,     46339,     46339,     46323,     46323,
                  46279,     46279,     46274,     46274,     46274,     46269,
                  46268,     46268,     46267,     46267,     46267,     46254,
                  46253,     46251,     46237,     46234,     46226,     46221,
                  46220,     46219,     46199,     46199],
             [ 47545018,  47545018,  47545015,  47545014,  47545013,  47545011,
               47545011,  47545008,  47545004,  47545004,  47544999,  47544995,
               47544981,  47544980,  47544975,  47544975,  47544959,  47544958,
               47544914,  47544914,  47544909,  47544909,  47544909,  47544904,
               47544903,  47544898,  47544897,  47544897,  47544897,  47544884,
               47544883,  47544881,  47544867,  47544864,  47544856,  47544851,
               47544850,  47544849,  47544829,  47544825],
            ]
                    # fmt: on
                ),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3019960 : 3019960 + 757],
            "actagggatgggagaggctcccagaacccagtaatgatgacattaagaaatacacaacagttgggaaatggaacccaaagagaacacctccagtagataagcatgacccccagttgagggatgggcccatgcacccatcttaaaattttggacccagaattattcttctcaaaaggaaatgcagggatgaaaatggagcagagactggaagaaaggccaaccagagactgccctaactcaggatccatcgcatgtgcaggcaccaaccccaacactattgctgatgccatgttgtacttgctgatggaagcctggcatggctgtcctctgagagtctcaaatgaggcacctgacagatgcagatacttacagccaaccaatggactgagccccgggacctcaataaaagaatgaggggatggcaaccccataggaagaacaacagtatcaactccctggactcctcagagctcccggggactaagccaccaactaaagagcatacataggctgctctgaggccccagatacatatgtagcagaggactgcctcagtgggaggggatgtgcttggtcttgtgaaggcttgatgctccagagaaggaggatgctagaggggtgaggtgggagtggatgggtgggtgggcaggggagcaccctcttagaggacaagggctctggggtgggggagctcatggagggggaactgggaaggagggagaacatttgaaatgtaaataaataaaataataaaaaa",
        )
        self.assertEqual(
            alignment[0],
            "actagggatgggagaggctcccagaacccagtaatgatgacattaagaaatacacaacagttgggaaatggaacccaaagagaacacctccagtagataagcatgacccccagttgagggatgggcccatgcacccatcttaaaattttggacccagaattattcttctcaaaaggaaatgcagggatgaaaatggagcagagactggaagaaaggccaaccagagactgccctaactcaggatccatcgcatgtgcaggcaccaaccccaacactattgctgatgccatgttgtacttgctgatggaagcctggcatggctgtcctctgagagtctcaaatgaggcacctgacagatgcagatacttacagccaaccaatggactgagccccgggacctcaataaaagaatgaggggatggcaaccccataggaagaacaacagtatcaactccctggactcctcagagctcccggggactaagccaccaactaaagagcatacataggctgctctgaggccccagatacatatgtagcagaggactgcctcagtgggaggggatgtgcttggtcttgtgaaggcttgatgctccagagaaggaggatgctagaggggtgaggtgggagtggatgggtgggtgggcaggggagcaccctcttagaggacaagggctctggggtgggggagctcatggagggggaactgggaaggagggagaacatttgaaatgtaaataaataaaataataaaaaa",
        )
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "felCat3.scaffold_205680")
        self.assertEqual(len(record.seq), 119354)
        self.assertEqual(segment, (46199, 46172))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "canFam2.chr1")
        self.assertEqual(len(record.seq), 125616256)
        self.assertEqual(segment, (47544825, 47544825))
        self.assertEqual(status, "C")
        empty = alignment.annotations["empty"][2]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ornAna1.chr2")
        self.assertEqual(len(record.seq), 54797317)
        self.assertEqual(segment, (40046122, 40040432))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][3]
        (record, segment, status) = empty
        self.assertEqual(record.id, "calJac1.Contig6394")
        self.assertEqual(len(record.seq), 133105)
        self.assertEqual(segment, (10461, 10482))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][4]
        (record, segment, status) = empty
        self.assertEqual(record.id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(record.seq), 498454)
        self.assertEqual(segment, (324116, 323493))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][5]
        (record, segment, status) = empty
        self.assertEqual(record.id, "echTel1.scaffold_288249")
        self.assertEqual(len(record.seq), 100002)
        self.assertEqual(segment, (87661, 95225))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][6]
        (record, segment, status) = empty
        self.assertEqual(record.id, "cavPor2.scaffold_216473")
        self.assertEqual(len(record.seq), 10026)
        self.assertEqual(segment, (1978, 606))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][7]
        (record, segment, status) = empty
        self.assertEqual(record.id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(record.seq), 359464)
        self.assertEqual(segment, (171686, 171328))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][8]
        (record, segment, status) = empty
        self.assertEqual(record.id, "hg18.chr6")
        self.assertEqual(len(record.seq), 170899992)
        self.assertEqual(segment, (155025365, 155025365))
        self.assertEqual(status, "C")
        empty = alignment.annotations["empty"][9]
        (record, segment, status) = empty
        self.assertEqual(record.id, "panTro2.chr6")
        self.assertEqual(len(record.seq), 173908612)
        self.assertEqual(segment, (157515419, 157515419))
        self.assertEqual(status, "C")
        empty = alignment.annotations["empty"][10]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ponAbe2.chr6")
        self.assertEqual(len(record.seq), 174210431)
        self.assertEqual(segment, (158037586, 158037565))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][11]
        (record, segment, status) = empty
        self.assertEqual(record.id, "oryCun1.scaffold_156751")
        self.assertEqual(len(record.seq), 4726)
        self.assertEqual(segment, (4300, 4263))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 1)
        self.assertEqual(len(alignment.annotations["empty"]), 12)
        self.assertTrue(
            numpy.array_equal(alignment.coordinates, numpy.array([[3019960, 3020717]]))
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 8951)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3020717 : 3020717 + 44],
            "TGTCAAACATGCATAAAGATATACTGAGGAGCCCATGAATTTTA",
        )
        self.assertEqual(alignment[0], "TGTCAAACATGCATAAAGATATACT-GAGGAGCCCATGAATTTTA")
        self.assertEqual(alignment.sequences[1].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[1].seq), 125616256)
        self.assertEqual(
            alignment.sequences[1].seq[47544792 : 47544792 + 33],
            "TAAAATTCATGGGCCCCTCTATTATGTTAAACa",
        )
        self.assertEqual(alignment[1], "tGTT------------TAACATAATAGAGGGGCCCATGAATTTTA")
        self.assertEqual(
            alignment.sequences[1].annotations["quality"],
            "999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 9)
        self.assertEqual(alignment.sequences[2].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[2].seq), 173908612)
        self.assertEqual(
            alignment.sequences[2].seq[157515381 : 157515381 + 38],
            "TAAAATTCATGGACCCCTCTAGTATATTTAAAATTTTT",
        )
        self.assertEqual(alignment[2], "----AAAAAT---TTTAAATATACTAGAGGGGTCCATGAATTTTA")
        self.assertEqual(
            alignment.sequences[2].annotations["quality"],
            "99999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[2].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[2].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[2].annotations["rightCount"], 11)
        self.assertEqual(alignment.sequences[3].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[3].seq), 170899992)
        self.assertEqual(
            alignment.sequences[3].seq[155025327 : 155025327 + 38],
            "TAAAATTCATGGACCCCTCTAGTATATTTAAAATTTTT",
        )
        self.assertEqual(alignment[3], "----AAAAAT---TTTAAATATACTAGAGGGGTCCATGAATTTTA")
        self.assertEqual(alignment.sequences[3].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[3].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[3].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[3].annotations["rightCount"], 11)
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "felCat3.scaffold_205680")
        self.assertEqual(len(record.seq), 119354)
        self.assertEqual(segment, (46199, 46172))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ornAna1.chr2")
        self.assertEqual(len(record.seq), 54797317)
        self.assertEqual(segment, (40046122, 40040432))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][2]
        (record, segment, status) = empty
        self.assertEqual(record.id, "calJac1.Contig6394")
        self.assertEqual(len(record.seq), 133105)
        self.assertEqual(segment, (10461, 10482))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][3]
        (record, segment, status) = empty
        self.assertEqual(record.id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(record.seq), 498454)
        self.assertEqual(segment, (324116, 323493))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][4]
        (record, segment, status) = empty
        self.assertEqual(record.id, "echTel1.scaffold_288249")
        self.assertEqual(len(record.seq), 100002)
        self.assertEqual(segment, (87661, 95225))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][5]
        (record, segment, status) = empty
        self.assertEqual(record.id, "cavPor2.scaffold_216473")
        self.assertEqual(len(record.seq), 10026)
        self.assertEqual(segment, (1978, 606))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][6]
        (record, segment, status) = empty
        self.assertEqual(record.id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(record.seq), 359464)
        self.assertEqual(segment, (171686, 171328))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][7]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ponAbe2.chr6")
        self.assertEqual(len(record.seq), 174210431)
        self.assertEqual(segment, (158037586, 158037565))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][8]
        (record, segment, status) = empty
        self.assertEqual(record.id, "oryCun1.scaffold_156751")
        self.assertEqual(len(record.seq), 4726)
        self.assertEqual(segment, (4300, 4263))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 4)
        self.assertEqual(len(alignment.annotations["empty"]), 9)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
            [[  3020717,   3020721,   3020727,   3020730,   3020733,   3020742,
                3020742,   3020761],
             [ 47544825,  47544821,  47544821,  47544821,  47544821,  47544812,
               47544811,  47544792],
             [157515419, 157515419, 157515413, 157515413, 157515410, 157515401,
              157515400, 157515381],
             [155025365, 155025365, 155025359, 155025359, 155025356, 155025347,
              155025346, 155025327],
            ]
                    # fmt: on
                ),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3020761 : 3020761 + 157],
            "TATATATGCTATCCGTGTGCTGTGATTTTTGTTTTAAATGTTATTTTATGTATATGcaagattttgcattgtagcagaaggtggcttcaaactcacgatcctcctgcctcagccttccaagtgctgagatcatacctctgcaccatcctgcccACCT",
        )
        self.assertEqual(
            alignment[0],
            "TATATATGCTATCCGTGTGCTGTGATTTTTGTTTTAAATGTTATTTTATGTATATGcaagattttgcattgtagcagaaggtggcttcaaactcacgatcctcctgcctcagccttccaagtgctgagatcatacctctgcaccatcctgcccACCT",
        )
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "felCat3.scaffold_205680")
        self.assertEqual(len(record.seq), 119354)
        self.assertEqual(segment, (46199, 46172))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "canFam2.chr1")
        self.assertEqual(len(record.seq), 125616256)
        self.assertEqual(segment, (47544792, 47544783))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][2]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ornAna1.chr2")
        self.assertEqual(len(record.seq), 54797317)
        self.assertEqual(segment, (40046122, 40040432))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][3]
        (record, segment, status) = empty
        self.assertEqual(record.id, "calJac1.Contig6394")
        self.assertEqual(len(record.seq), 133105)
        self.assertEqual(segment, (10461, 10482))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][4]
        (record, segment, status) = empty
        self.assertEqual(record.id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(record.seq), 498454)
        self.assertEqual(segment, (324116, 323493))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][5]
        (record, segment, status) = empty
        self.assertEqual(record.id, "echTel1.scaffold_288249")
        self.assertEqual(len(record.seq), 100002)
        self.assertEqual(segment, (87661, 95225))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][6]
        (record, segment, status) = empty
        self.assertEqual(record.id, "cavPor2.scaffold_216473")
        self.assertEqual(len(record.seq), 10026)
        self.assertEqual(segment, (1978, 606))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][7]
        (record, segment, status) = empty
        self.assertEqual(record.id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(record.seq), 359464)
        self.assertEqual(segment, (171686, 171328))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][8]
        (record, segment, status) = empty
        self.assertEqual(record.id, "hg18.chr6")
        self.assertEqual(len(record.seq), 170899992)
        self.assertEqual(segment, (155025327, 155025316))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][9]
        (record, segment, status) = empty
        self.assertEqual(record.id, "panTro2.chr6")
        self.assertEqual(len(record.seq), 173908612)
        self.assertEqual(segment, (157515381, 157515370))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][10]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ponAbe2.chr6")
        self.assertEqual(len(record.seq), 174210431)
        self.assertEqual(segment, (158037586, 158037565))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][11]
        (record, segment, status) = empty
        self.assertEqual(record.id, "oryCun1.scaffold_156751")
        self.assertEqual(len(record.seq), 4726)
        self.assertEqual(segment, (4300, 4263))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 1)
        self.assertEqual(len(alignment.annotations["empty"]), 12)
        self.assertTrue(
            numpy.array_equal(alignment.coordinates, numpy.array([[3020761, 3020918]]))
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 85471)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3020918 : 3020918 + 96],
            "GAGGTTTGTGACTTTTAATACTGATTGTTATCTAACATCACAGAATTCTCAGTTCTTAAGGAAACAATTGTTCTGTGTGTTATTTGTCTAGGAGGA",
        )
        self.assertEqual(
            alignment[0],
            "GAGGTTTGTGACTTTTAATA----------CTGATTGTTATCTAACATCACAGAATTCTCAGTTCTTAAGGAAACAATTGTTCTGTGTGTTATTTGTCTAGGAGGA",
        )
        self.assertEqual(alignment.sequences[1].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[1].seq), 170899992)
        self.assertEqual(
            alignment.sequences[1].seq[155025243 : 155025243 + 73],
            "TTGTCCTAAATAATTAATAAGTCAAACATGTTTTTCCTTAAAAGCTGAGGATTGTGCAGTATTAAATAACCAT",
        )
        self.assertEqual(
            alignment[1],
            "---------------------------------ATGGTTATTTAATACTGCACAATCCTCAGCTTTTAAGGAAAAACATGTTTGACTTATTAATTATTTAGGACAA",
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 11)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[2].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[2].seq), 173908612)
        self.assertEqual(
            alignment.sequences[2].seq[157515297 : 157515297 + 73],
            "TTGTCCTAAATAATTAATAAGTCAAACATGTTTTTCCTTAAAAGCTGAGGATTGTGCAGTATTAAATAACCAT",
        )
        self.assertEqual(
            alignment[2],
            "---------------------------------ATGGTTATTTAATACTGCACAATCCTCAGCTTTTAAGGAAAAACATGTTTGACTTATTAATTATTTAGGACAA",
        )
        self.assertEqual(
            alignment.sequences[2].annotations["quality"],
            "9999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[2].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[2].annotations["leftCount"], 11)
        self.assertEqual(alignment.sequences[2].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[3].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[3].seq), 125616256)
        self.assertEqual(
            alignment.sequences[3].seq[47544708 : 47544708 + 75],
            "TTGTCCTAAGTAATTAACAAATCACTTTTTCCTTAAGAGCTGAGAACTTCATAATGGTACAGAATTATTTTATTA",
        )
        self.assertEqual(
            alignment[3],
            "--------------------------TAATAAAATAATTCTGTACCATTATGAAGTTCTCAGCTCTTAAGGAAAAA-----GTGATTTGTTAATTACTTAGGACAA",
        )
        self.assertEqual(
            alignment.sequences[3].annotations["quality"],
            "999999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[3].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[3].annotations["leftCount"], 9)
        self.assertEqual(alignment.sequences[3].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[3].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[4].id, "felCat3.scaffold_205680")
        self.assertEqual(len(alignment.sequences[4].seq), 119354)
        self.assertEqual(
            alignment.sequences[4].seq[46096 : 46096 + 76],
            "TTGTCCTAAATAATTAATAAATTACTTTGTGATATTAAAGAATTATTTCATTACTTCACCATTAAAATTCATGGAC",
        )
        self.assertEqual(
            alignment[4],
            "---GTCCATGAATTTTAATGGTGAAGTAATGAAATAATTCTTTAATATCAC----------------------AAA-----GTAATTTATTAATTATTTAGGACAA",
        )
        self.assertEqual(
            alignment.sequences[4].annotations["quality"],
            "9999899999999999999999999999999999999999999999999999999999999999999999999769",
        )
        self.assertEqual(alignment.sequences[4].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[4].annotations["leftCount"], 27)
        self.assertEqual(alignment.sequences[4].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[4].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[5].id, "cavPor2.scaffold_216473")
        self.assertEqual(len(alignment.sequences[5].seq), 10026)
        self.assertEqual(
            alignment.sequences[5].seq[524 : 524 + 82],
            "tttaaaaaaatcacactttTTCCATCAGAACTGAAAACTTTGTAATATTAatcttcttttcctattaaaATTTGCAGAACTC",
        )
        self.assertEqual(
            alignment[5],
            "GAGTTCTGCAAATtttaata----------ggaaaagaagatTAATATTACAAAGTTTTCAGTTCTGATGGAAaaa----------gtgtgatttttttaaa----",
        )
        self.assertEqual(
            alignment.sequences[5].annotations["quality"],
            "9967996679966856646585678288383465687882688656636765583677657766676965686997776684",
        )
        self.assertEqual(alignment.sequences[5].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[5].annotations["leftCount"], 1372)
        self.assertEqual(alignment.sequences[5].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[5].annotations["rightCount"], 28)
        self.assertEqual(alignment.sequences[6].id, "oryCun1.scaffold_156751")
        self.assertEqual(len(alignment.sequences[6].seq), 4726)
        self.assertEqual(
            alignment.sequences[6].seq[4161 : 4161 + 102],
            "TCATCCCAAGTCATTAACAAGTCAAACACTTTTTTCCTTAAGACTTGAGGATTTTCCAACATTAAATAATTATCTCATTGCTTTATTAAAATTCATGGATCT",
        )
        self.assertEqual(
            alignment[6],
            "-AGATCCATGAATTTTAATA---AAGCAATGAGATAATTATTTAATGTTGGAAAATCCTCAAGTCTTAAGGAAAAAAGTGTTTGACTTGTTAATGACTTGGGATGA",
        )
        self.assertEqual(
            alignment.sequences[6].annotations["quality"],
            "999999799999999999999999999999999999999999999999999998677999969999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[6].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[6].annotations["leftCount"], 37)
        self.assertEqual(alignment.sequences[6].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[6].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[7].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[7].seq), 133105)
        self.assertEqual(
            alignment.sequences[7].seq[10482 : 10482 + 100],
            "GGGGTCCATGAATTTTAATAGTAACAAAATGGTTATTCAATATTGCAAAATCCTCAGCTTTAAGGAAAAACATATTTGATTTGTTAATTATTTAGGACAA",
        )
        self.assertEqual(
            alignment[7],
            "GGGGTCCATGAATTTTAATA-----GTAACAAAATGGTTATTCAATATTGCAAAATCCTCAGC-TTTAAGGAAAAACATATTTGATTTGTTAATTATTTAGGACAA",
        )
        self.assertEqual(alignment.sequences[7].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[7].annotations["leftCount"], 21)
        self.assertEqual(alignment.sequences[7].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[7].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[8].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[8].seq), 174210431)
        self.assertEqual(
            alignment.sequences[8].seq[158037464 : 158037464 + 101],
            "TTGTCCTAAATAATTAATAAGTCAAACACGTTTTTCCTTAAAAGCTGAGGATTTTGCAATATTAAATAACCATTTTATTACTATTAAAATTCATGGACCCC",
        )
        self.assertEqual(
            alignment[8],
            "GGGGTCCATGAATTTTAATA-----GTAATAAAATGGTTATTTAATATTGCAAAATCCTCAGCTTTTAAGGAAAAACGTGTTTGACTTATTAATTATTTAGGACAA",
        )
        self.assertEqual(alignment.sequences[8].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[8].annotations["leftCount"], 21)
        self.assertEqual(alignment.sequences[8].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[8].annotations["rightCount"], 0)
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ornAna1.chr2")
        self.assertEqual(len(record.seq), 54797317)
        self.assertEqual(segment, (40046122, 40040432))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(record.seq), 498454)
        self.assertEqual(segment, (324116, 323493))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][2]
        (record, segment, status) = empty
        self.assertEqual(record.id, "echTel1.scaffold_288249")
        self.assertEqual(len(record.seq), 100002)
        self.assertEqual(segment, (87661, 95225))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][3]
        (record, segment, status) = empty
        self.assertEqual(record.id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(record.seq), 359464)
        self.assertEqual(segment, (171686, 171328))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 9)
        self.assertEqual(len(alignment.annotations["empty"]), 4)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
            [[  3020918,   3020919,   3020921,   3020938,   3020938,   3020938,
                3020938,   3020938,   3020941,   3020959,   3020971,   3020972,
                3020981,   3020984,   3020989,   3020994,   3021010,   3021014],
             [155025316, 155025316, 155025316, 155025316, 155025316, 155025316,
              155025316, 155025316, 155025316, 155025298, 155025286, 155025285,
              155025276, 155025273, 155025268, 155025263, 155025247, 155025243],
             [157515370, 157515370, 157515370, 157515370, 157515370, 157515370,
              157515370, 157515370, 157515370, 157515352, 157515340, 157515339,
              157515330, 157515327, 157515322, 157515317, 157515301, 157515297],
             [ 47544783,  47544783,  47544783,  47544783,  47544783,  47544783,
               47544783,  47544779,  47544776,  47544758,  47544746,  47544745,
               47544736,  47544733,  47544733,  47544728,  47544712,  47544708],
             [    46172,     46172,     46172,     46155,     46152,     46150,
                  46149,     46145,     46142,     46124,     46124,     46124,
                  46124,     46121,     46121,     46116,     46100,     46096],
             [      606,       605,       603,       586,       586,       586,
                    586,       586,       583,       565,       553,       552,
                    543,       540,       540,       540,       524,       524],
             [     4263,      4263,      4261,      4244,      4244,      4242,
                   4241,      4237,      4234,      4216,      4204,      4203,
                   4194,      4191,      4186,      4181,      4165,      4161],
             [    10482,     10483,     10485,     10502,     10502,     10502,
                  10503,     10507,     10510,     10528,     10540,     10540,
                  10549,     10552,     10557,     10562,     10578,     10582],
             [158037565, 158037564, 158037562, 158037545, 158037545, 158037545,
              158037544, 158037540, 158037537, 158037519, 158037507, 158037506,
              158037497, 158037494, 158037489, 158037484, 158037468, 158037464],
            ]
                    # fmt: on
                ),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 105724)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3021014 : 3021014 + 40],
            "ACCTTGGTGACGCCACTGGATTTTGTATGACTGAATACTG",
        )
        self.assertEqual(alignment[0], "ACCTTGGTGACGCCACTGGATTTTGTATGACTGAATACTG")
        self.assertEqual(alignment.sequences[1].id, "oryCun1.scaffold_156751")
        self.assertEqual(len(alignment.sequences[1].seq), 4726)
        self.assertEqual(
            alignment.sequences[1].seq[4121 : 4121 + 40],
            "AATTCCTAAATCATCCAAGTTGGAATGACTTCATCAAGAT",
        )
        self.assertEqual(alignment[1], "ATCTTGATGAAGTCATTCCAACTTGGATGATTTAGGAATT")
        self.assertEqual(
            alignment.sequences[1].annotations["quality"],
            "9999999999999969999699999999999999999999",
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[2].id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(alignment.sequences[2].seq), 359464)
        self.assertEqual(
            alignment.sequences[2].seq[171288 : 171288 + 40],
            "CATTCCTAAATCACATAAATTGGAATAACTTCACCAAGAT",
        )
        self.assertEqual(alignment[2], "ATCTTGGTGAAGTTATTCCAATTTATGTGATTTAGGAATG")
        self.assertEqual(
            alignment.sequences[2].annotations["quality"],
            "9999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[2].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[2].annotations["leftCount"], 358)
        self.assertEqual(alignment.sequences[2].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[3].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[3].seq), 133105)
        self.assertEqual(
            alignment.sequences[3].seq[10582 : 10582 + 40],
            "ATTTTGGTGAAGTTATTCCAACTTGTGTGGCTTAGGAATG",
        )
        self.assertEqual(alignment[3], "ATTTTGGTGAAGTTATTCCAACTTGTGTGGCTTAGGAATG")
        self.assertEqual(alignment.sequences[3].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[3].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[3].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[3].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[4].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[4].seq), 170899992)
        self.assertEqual(
            alignment.sequences[4].seq[155025203 : 155025203 + 40],
            "CATTCCTAAGCCATGCAAGTTGGAATAACTTCACCAAAAT",
        )
        self.assertEqual(alignment[4], "ATTTTGGTGAAGTTATTCCAACTTGCATGGCTTAGGAATG")
        self.assertEqual(alignment.sequences[4].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[4].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[4].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[4].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[5].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[5].seq), 173908612)
        self.assertEqual(
            alignment.sequences[5].seq[157515257 : 157515257 + 40],
            "CATTCCTAAGCCATGCAAGTTGGAATAACTTCACCAAAAT",
        )
        self.assertEqual(alignment[5], "ATTTTGGTGAAGTTATTCCAACTTGCATGGCTTAGGAATG")
        self.assertEqual(
            alignment.sequences[5].annotations["quality"],
            "9999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[5].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[5].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[5].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[5].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[6].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[6].seq), 174210431)
        self.assertEqual(
            alignment.sequences[6].seq[158037424 : 158037424 + 40],
            "CATTCCTAAGCCATGCAAGTTGGAATAACTTCACCAAAAT",
        )
        self.assertEqual(alignment[6], "ATTTTGGTGAAGTTATTCCAACTTGCATGGCTTAGGAATG")
        self.assertEqual(alignment.sequences[6].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[6].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[6].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[6].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[7].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[7].seq), 125616256)
        self.assertEqual(
            alignment.sequences[7].seq[47544668 : 47544668 + 40],
            "CATTCCTAAATCATGCGAGTCAGAATGACTTCACTGAGAT",
        )
        self.assertEqual(alignment[7], "ATCTCAGTGAAGTCATTCTGACTCGCATGATTTAGGAATG")
        self.assertEqual(
            alignment.sequences[7].annotations["quality"],
            "9999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[7].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[7].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[7].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[7].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[8].id, "felCat3.scaffold_205680")
        self.assertEqual(len(alignment.sequences[8].seq), 119354)
        self.assertEqual(
            alignment.sequences[8].seq[46056 : 46056 + 40],
            "CATACCTAAATCATGCAAGTCAGAATAACTTCACTGAGAT",
        )
        self.assertEqual(alignment[8], "ATCTCAGTGAAGTTATTCTGACTTGCATGATTTAGGTATG")
        self.assertEqual(
            alignment.sequences[8].annotations["quality"],
            "9999989999989988999997999979997996167779",
        )
        self.assertEqual(alignment.sequences[8].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[8].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[8].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[8].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[9].id, "dasNov1.scaffold_56749")
        self.assertEqual(len(alignment.sequences[9].seq), 10470)
        self.assertEqual(
            alignment.sequences[9].seq[7314 : 7314 + 40],
            "CATCGCTCAGTCATACGAGTCGGAATGATTTCACTGATGT",
        )
        self.assertEqual(alignment[9], "ACATCAGTGAAATCATTCCGACTCGTATGACTGAGCGATG")
        self.assertEqual(
            alignment.sequences[9].annotations["quality"],
            "9759855999977756667495765475885678385647",
        )
        self.assertEqual(alignment.sequences[9].annotations["leftStatus"], "N")
        self.assertEqual(alignment.sequences[9].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[9].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[9].annotations["rightCount"], 0)
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ornAna1.chr2")
        self.assertEqual(len(record.seq), 54797317)
        self.assertEqual(segment, (40046122, 40040432))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(record.seq), 498454)
        self.assertEqual(segment, (324116, 323493))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][2]
        (record, segment, status) = empty
        self.assertEqual(record.id, "echTel1.scaffold_288249")
        self.assertEqual(len(record.seq), 100002)
        self.assertEqual(segment, (87661, 95225))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][3]
        (record, segment, status) = empty
        self.assertEqual(record.id, "cavPor2.scaffold_216473")
        self.assertEqual(len(record.seq), 10026)
        self.assertEqual(segment, (524, 496))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 10)
        self.assertEqual(len(alignment.annotations["empty"]), 4)
        self.assertTrue(
            numpy.array_equal(
                # fmt: off
                alignment.coordinates, numpy.array([[  3021014,   3021054],
                                                    [     4161,      4121],
                                                    [   171328,    171288],
                                                    [    10582,     10622],
                                                    [155025243, 155025203],
                                                    [157515297, 157515257],
                                                    [158037464, 158037424],
                                                    [ 47544708,  47544668],
                                                    [    46096,     46056],
                                                    [     7354,      7314],
                                                   ])
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 115790)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3021054 : 3021054 + 50],
            "CTCATTTGGGAACTTACAGGTCAGCAAAGGCTTCCAGGACTTACATGCAG",
        )
        self.assertEqual(
            alignment[0],
            "CTCATTTGGGAACTTACAGGTCAGCAAAGGCTTCCAG--------------------GACTTACATGCAG",
        )
        self.assertEqual(alignment.sequences[1].id, "cavPor2.scaffold_216473")
        self.assertEqual(len(alignment.sequences[1].seq), 10026)
        self.assertEqual(
            alignment.sequences[1].seq[458 : 458 + 38],
            "TCACAAGTATTTATAGCAACTCTGAACTCCCAAATGAG",
        )
        self.assertEqual(
            alignment[1],
            "CTCATTTGGGAGTTCAGAGTT--------GCTATAAA--------------------TACTTGTGA----",
        )
        self.assertEqual(
            alignment.sequences[1].annotations["quality"],
            "78877766789884666566698766677876665669",
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 28)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[2].id, "oryCun1.scaffold_156751")
        self.assertEqual(len(alignment.sequences[2].seq), 4726)
        self.assertEqual(
            alignment.sequences[2].seq[4052 : 4052 + 69],
            "ATACGTGTAATTAAGGAGGAAAAAAGAAGACTCCATTAGGCAGTTATAAAATGTAAGAGCCCAAATTAG",
        )
        self.assertEqual(
            alignment[2],
            "CTAATTTGGGCTCTTACATTTTATAACTGCCTAATGGAGTCTTCTTTTTTCCT-CCTTAATTACACGTAT",
        )
        self.assertEqual(
            alignment.sequences[2].annotations["quality"],
            "999999999999979997999999999999999999999997997999999979979999999495999",
        )
        self.assertEqual(alignment.sequences[2].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[2].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[3].id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(alignment.sequences[3].seq), 359464)
        self.assertEqual(
            alignment.sequences[3].seq[171219 : 171219 + 69],
            "GTACATAGAGTCTGGGGAAGGGACCAGTGGGCTCCATTAAGCCTTTATGAACCTGTGTCCCCAAGTTAG",
        )
        self.assertEqual(
            alignment[3],
            "CTAACTTGGGGACACAGG-TTCATAAAGGCTTAATGGAGCCCACTGGTCCCTTCCCCAGACTCTATGTAC",
        )
        self.assertEqual(
            alignment.sequences[3].annotations["quality"],
            "999999999999999999999999999999999999998999999689999999999999999999999",
        )
        self.assertEqual(alignment.sequences[3].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[3].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[3].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[3].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[4].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[4].seq), 133105)
        self.assertEqual(
            alignment.sequences[4].seq[10622 : 10622 + 70],
            "CCAGTTTGGGGACTTAGATTTTCTAACTGCCTAATGAAGTCTGCTCTTTCCTTCCCTAGCCTATATGTAT",
        )
        self.assertEqual(
            alignment[4],
            "CCAGTTTGGGGACTTAGATTTTCTAACTGCCTAATGAAGTCTGCTCTTTCCTTCCCTAGCCTATATGTAT",
        )
        self.assertEqual(alignment.sequences[4].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[4].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[4].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[4].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[5].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[5].seq), 170899992)
        self.assertEqual(
            alignment.sequences[5].seq[155025133 : 155025133 + 70],
            "ATACATATAGTCTAGGGAAGGAAGGAGCAGACTTCATTAGGCAGTTAGGAAATGTAAATCCCCAAATTGT",
        )
        self.assertEqual(
            alignment[5],
            "ACAATTTGGGGATTTACATTTCCTAACTGCCTAATGAAGTCTGCTCCTTCCTTCCCTAGACTATATGTAT",
        )
        self.assertEqual(alignment.sequences[5].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[5].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[5].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[5].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[6].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[6].seq), 173908612)
        self.assertEqual(
            alignment.sequences[6].seq[157515187 : 157515187 + 70],
            "ATACATATAGTCTAGGGAAGGAAAGAGCAGACTTCATTAGGCAGTTAGGAAATGTAAATCCCCAAATTGG",
        )
        self.assertEqual(
            alignment[6],
            "CCAATTTGGGGATTTACATTTCCTAACTGCCTAATGAAGTCTGCTCTTTCCTTCCCTAGACTATATGTAT",
        )
        self.assertEqual(
            alignment.sequences[6].annotations["quality"],
            "9999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[6].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[6].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[6].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[6].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[7].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[7].seq), 174210431)
        self.assertEqual(
            alignment.sequences[7].seq[158037354 : 158037354 + 70],
            "ATACATATAGTCTAGGGAAGGAAAGAGCAGACTTCATTAGGCAGTTAGGAAATGTAAATCCCCAAATTGG",
        )
        self.assertEqual(
            alignment[7],
            "CCAATTTGGGGATTTACATTTCCTAACTGCCTAATGAAGTCTGCTCTTTCCTTCCCTAGACTATATGTAT",
        )
        self.assertEqual(alignment.sequences[7].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[7].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[7].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[7].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[8].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[8].seq), 125616256)
        self.assertEqual(
            alignment.sequences[8].seq[47544598 : 47544598 + 70],
            "ATACATATGGTGTAGGGAAGAAAACATCAGCTTTCATTAAGCAATTATGAAATTTGAGTCACCAAATTAG",
        )
        self.assertEqual(
            alignment[8],
            "CTAATTTGGTGACTCAAATTTCATAATTGCTTAATGAAAGCTGATGTTTTCTTCCCTACACCATATGTAT",
        )
        self.assertEqual(
            alignment.sequences[8].annotations["quality"],
            "9999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[8].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[8].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[8].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[8].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[9].id, "felCat3.scaffold_205680")
        self.assertEqual(len(alignment.sequences[9].seq), 119354)
        self.assertEqual(
            alignment.sequences[9].seq[45986 : 45986 + 70],
            "ATGCATATAGTCTTGGGAAGGAAACAGCAGGCTTCCTTAAGCAATTATGAAATTCAAGTCACCAAATTAG",
        )
        self.assertEqual(
            alignment[9],
            "CTAATTTGGTGACTTGAATTTCATAATTGCTTAAGGAAGCCTGCTGTTTCCTTCCCAAGACTATATGCAT",
        )
        self.assertEqual(
            alignment.sequences[9].annotations["quality"],
            "9769999999975699868977669777966666596959759669595666758736585676666655",
        )
        self.assertEqual(alignment.sequences[9].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[9].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[9].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[9].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[10].id, "dasNov1.scaffold_56749")
        self.assertEqual(len(alignment.sequences[10].seq), 10470)
        self.assertEqual(
            alignment.sequences[10].seq[7244 : 7244 + 70],
            "ATACTTATAGTCTAGAGAAGAAAACAGCAGGCTTCGGTTAGCAATTTTAAAATGTGAGTCCCCCAATTAG",
        )
        self.assertEqual(
            alignment[10],
            "CTAATTGGGGGACTCACATTTTAAAATTGCTAACCGAAGCCTGCTGTTTTCTTCTCTAGACTATAAGTAT",
        )
        self.assertEqual(
            alignment.sequences[10].annotations["quality"],
            "5556576999664654656985688667655565647767537567688856666555556565555656",
        )
        self.assertEqual(alignment.sequences[10].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[10].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[10].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[10].annotations["rightCount"], 0)
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ornAna1.chr2")
        self.assertEqual(len(record.seq), 54797317)
        self.assertEqual(segment, (40046122, 40040432))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(record.seq), 498454)
        self.assertEqual(segment, (324116, 323493))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][2]
        (record, segment, status) = empty
        self.assertEqual(record.id, "echTel1.scaffold_288249")
        self.assertEqual(len(record.seq), 100002)
        self.assertEqual(segment, (87661, 95225))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 11)
        self.assertEqual(len(alignment.annotations["empty"]), 3)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
            [[  3021054,   3021072,   3021073,   3021075,   3021083,   3021091,
                3021091,   3021091,   3021091,   3021100,   3021104],
             [      496,       478,       477,       475,       475,       467,
                    467,       467,       467,       458,       458],
             [     4121,      4103,      4102,      4100,      4092,      4084,
                   4068,      4068,      4065,      4056,      4052],
             [   171288,    171270,    171270,    171268,    171260,    171252,
                 171236,    171235,    171232,    171223,    171219],
             [    10622,     10640,     10641,     10643,     10651,     10659,
                  10675,     10676,     10679,     10688,     10692],
             [155025203, 155025185, 155025184, 155025182, 155025174, 155025166,
              155025150, 155025149, 155025146, 155025137, 155025133],
             [157515257, 157515239, 157515238, 157515236, 157515228, 157515220,
              157515204, 157515203, 157515200, 157515191, 157515187],
             [158037424, 158037406, 158037405, 158037403, 158037395, 158037387,
              158037371, 158037370, 158037367, 158037358, 158037354],
             [ 47544668,  47544650,  47544649,  47544647,  47544639,  47544631,
               47544615,  47544614,  47544611,  47544602,  47544598],
             [    46056,     46038,     46037,     46035,     46027,     46019,
                  46003,     46002,     45999,     45990,     45986],
             [     7314,      7296,      7295,      7293,      7285,      7277,
                   7261,      7260,      7257,      7248,      7244],
            ]
                    # fmt: on
                ),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 44222)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3021104 : 3021104 + 32],
            "CTGTTAGTGCTGTTTTAATGTACCTCGCAGTA",
        )
        self.assertEqual(alignment[0], "CTGTTAGTGCTGTTTT---AATGTACCTCGCAGTA")
        self.assertEqual(alignment.sequences[1].id, "cavPor2.scaffold_216473")
        self.assertEqual(len(alignment.sequences[1].seq), 10026)
        self.assertEqual(
            alignment.sequences[1].seq[431 : 431 + 27], "TATAAAAATTTACATTAAGAAAGTAAT"
        )
        self.assertEqual(alignment[1], "-----ATTACTTTCTT---AATGTAAATTTTTATA")
        self.assertEqual(
            alignment.sequences[1].annotations["quality"], "554558687467957999989884575"
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[2].id, "oryCun1.scaffold_156751")
        self.assertEqual(len(alignment.sequences[2].seq), 4726)
        self.assertEqual(
            alignment.sequences[2].seq[4022 : 4022 + 30],
            "TATAAAATATTACATTGAGAATACTAGCAA",
        )
        self.assertEqual(alignment[2], "TTGCTAGTA--TTCTC---AATGTAATATTTTATA")
        self.assertEqual(
            alignment.sequences[2].annotations["quality"],
            "999979999966656999999999999999",
        )
        self.assertEqual(alignment.sequences[2].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[2].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[3].id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(alignment.sequences[3].seq), 359464)
        self.assertEqual(
            alignment.sequences[3].seq[171187 : 171187 + 32],
            "TATAAAACATTATATTAAGGGAACACCAGCAA",
        )
        self.assertEqual(alignment[3], "TTGCTGGTGTTCCCTT---AATATAATGTTTTATA")
        self.assertEqual(
            alignment.sequences[3].annotations["quality"],
            "99999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[3].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[3].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[3].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[3].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[4].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[4].seq), 133105)
        self.assertEqual(
            alignment.sequences[4].seq[10692 : 10692 + 30],
            "TTACTCGTGCCCTTAATATAGCATTTTATA",
        )
        self.assertEqual(alignment[4], "TTACTCGTG--CCCTT---AATATAGCATTTTATA")
        self.assertEqual(alignment.sequences[4].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[4].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[4].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[4].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[5].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[5].seq), 170899992)
        self.assertEqual(
            alignment.sequences[5].seq[155025103 : 155025103 + 30],
            "TATAAAATGTGATATTAAGAGTACCAGCAA",
        )
        self.assertEqual(alignment[5], "TTGCTGGTA--CTCTT---AATATCACATTTTATA")
        self.assertEqual(alignment.sequences[5].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[5].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[5].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[5].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[6].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[6].seq), 173908612)
        self.assertEqual(
            alignment.sequences[6].seq[157515157 : 157515157 + 30],
            "TATAAAATGTGATATTAAGAGTACCAGCAA",
        )
        self.assertEqual(alignment[6], "TTGCTGGTA--CTCTT---AATATCACATTTTATA")
        self.assertEqual(
            alignment.sequences[6].annotations["quality"],
            "999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[6].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[6].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[6].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[6].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[7].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[7].seq), 174210431)
        self.assertEqual(
            alignment.sequences[7].seq[158037324 : 158037324 + 30],
            "TATAAAATGTTATATTAAGAGCACCAGCAA",
        )
        self.assertEqual(alignment[7], "TTGCTGGTG--CTCTT---AATATAACATTTTATA")
        self.assertEqual(alignment.sequences[7].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[7].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[7].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[7].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[8].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[8].seq), 125616256)
        self.assertEqual(
            alignment.sequences[8].seq[47544573 : 47544573 + 25],
            "TGCTACATTTTGAGAGCACCAGCAA",
        )
        self.assertEqual(alignment[8], "TTGCTGGTGCTCTCAA---AATGTAGCA-------")
        self.assertEqual(
            alignment.sequences[8].annotations["quality"], "9999999999999999999999999"
        )
        self.assertEqual(alignment.sequences[8].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[8].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[8].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[8].annotations["rightCount"], 196)
        self.assertEqual(alignment.sequences[9].id, "felCat3.scaffold_205680")
        self.assertEqual(len(alignment.sequences[9].seq), 119354)
        self.assertEqual(
            alignment.sequences[9].seq[45961 : 45961 + 25], "TGTTATATTTTGAGAGCACCAGCAA"
        )
        self.assertEqual(alignment[9], "TTGCTGGTGCTCTCAA---AATATAACA-------")
        self.assertEqual(
            alignment.sequences[9].annotations["quality"], "8677668566555658876555655"
        )
        self.assertEqual(alignment.sequences[9].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[9].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[9].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[9].annotations["rightCount"], 6)
        self.assertEqual(alignment.sequences[10].id, "dasNov1.scaffold_56749")
        self.assertEqual(len(alignment.sequences[10].seq), 10470)
        self.assertEqual(
            alignment.sequences[10].seq[7211 : 7211 + 33],
            "CATAAAATGTTACATTAATCAGATCACCAGCAA",
        )
        self.assertEqual(alignment[10], "TTGCTGGTG--ATCTGATTAATGTAACATTTTATG")
        self.assertEqual(
            alignment.sequences[10].annotations["quality"],
            "856647736775356546747663745776545",
        )
        self.assertEqual(alignment.sequences[10].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[10].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[10].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[10].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[11].id, "echTel1.scaffold_288249")
        self.assertEqual(len(alignment.sequences[11].seq), 100002)
        self.assertEqual(
            alignment.sequences[11].seq[95225 : 95225 + 21], "CTGTTAATGCTCTGTTTTATG"
        )
        self.assertEqual(alignment[11], "CTGTTAATG--CTCTG------------TTTTATG")
        self.assertEqual(
            alignment.sequences[11].annotations["quality"], "999999999999999999999"
        )
        self.assertEqual(alignment.sequences[11].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[11].annotations["leftCount"], 7564)
        self.assertEqual(alignment.sequences[11].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[11].annotations["rightCount"], 0)
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ornAna1.chr2")
        self.assertEqual(len(record.seq), 54797317)
        self.assertEqual(segment, (40046122, 40040432))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(record.seq), 498454)
        self.assertEqual(segment, (324116, 323493))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 12)
        self.assertEqual(len(alignment.annotations["empty"]), 2)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
            [[  3021104,   3021109,   3021113,   3021115,   3021120,   3021120,
                3021129,   3021136],
             [      458,       458,       454,       452,       447,       447,
                    438,       431],
             [     4052,      4047,      4043,      4043,      4038,      4038,
                   4029,      4022],
             [   171219,    171214,    171210,    171208,    171203,    171203,
                 171194,    171187],
             [    10692,     10697,     10701,     10701,     10706,     10706,
                  10715,     10722],
             [155025133, 155025128, 155025124, 155025124, 155025119, 155025119,
              155025110, 155025103],
             [157515187, 157515182, 157515178, 157515178, 157515173, 157515173,
              157515164, 157515157],
             [158037354, 158037349, 158037345, 158037345, 158037340, 158037340,
              158037331, 158037324],
             [ 47544598,  47544593,  47544589,  47544587,  47544582,  47544582,
               47544573,  47544573],
             [    45986,     45981,     45977,     45975,     45970,     45970,
                  45961,     45961],
             [     7244,      7239,      7235,      7235,      7230,      7227,
                   7218,      7211],
             [    95225,     95230,     95234,     95234,     95239,     95239,
                  95239,     95246],
            ]
                    # fmt: on
                ),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 43757)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3021136 : 3021136 + 44],
            "AGGCAAATGAGGTGATAAGATTGTGTTTACTCCCTCTGTGCTTG",
        )
        self.assertEqual(
            alignment[0],
            "AGGCAAATGAGGTGATAAGA-------TTGTGTT-----TAC----TCCCTCTGTGC----------TTG",
        )
        self.assertEqual(alignment.sequences[1].id, "cavPor2.scaffold_216473")
        self.assertEqual(len(alignment.sequences[1].seq), 10026)
        self.assertEqual(
            alignment.sequences[1].seq[388 : 388 + 43],
            "catacacacagagaatgAATATATCACTGTTATCTCATCTGCT",
        )
        self.assertEqual(
            alignment[1],
            "AG-CAGATGAGATAACAGTG-------ATATATT-----cat----tctctgtgtgt----------atg",
        )
        self.assertEqual(
            alignment.sequences[1].annotations["quality"],
            "4996766988786798889867956675666896967579888",
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[2].id, "oryCun1.scaffold_156751")
        self.assertEqual(len(alignment.sequences[2].seq), 4726)
        self.assertEqual(
            alignment.sequences[2].seq[3980 : 3980 + 42],
            "CACACAAATAAGAAATAGATGCCCTCATCACCTCATTTGCTT",
        )
        self.assertEqual(
            alignment[2],
            "AAGCAAATGAGGTGATGAGG---------GCATC-----TAT----TTCTTATTTGT----------GTG",
        )
        self.assertEqual(
            alignment.sequences[2].annotations["quality"],
            "999989999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[2].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[2].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[3].id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(alignment.sequences[3].seq), 359464)
        self.assertEqual(
            alignment.sequences[3].seq[171142 : 171142 + 45],
            "CACACGTGGGGGACAAGTACATATGTCCTGTGATCTCATTTGTTT",
        )
        self.assertEqual(
            alignment[3],
            "AAACAAATGAGATCACA-GG-------ACATATG-----TA--CTTGTCCCCCACGT----------GTG",
        )
        self.assertEqual(
            alignment.sequences[3].annotations["quality"],
            "999999999999999999999999999999939999999999999",
        )
        self.assertEqual(alignment.sequences[3].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[3].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[3].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[3].annotations["rightCount"], 15)
        self.assertEqual(alignment.sequences[4].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[4].seq), 133105)
        self.assertEqual(
            alignment.sequences[4].seq[10722 : 10722 + 45],
            "AAGCAAATGAGATCACAGCACATGTATATTTTTTCTCCGTGTGTG",
        )
        self.assertEqual(
            alignment[4],
            "AAGCAAATGAGATCACA----------GCACATG-----TATATTTTTTCTCCGTGT----------GTG",
        )
        self.assertEqual(alignment.sequences[4].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[4].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[4].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[4].annotations["rightCount"], 15)
        self.assertEqual(alignment.sequences[5].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[5].seq), 170899992)
        self.assertEqual(
            alignment.sequences[5].seq[155025056 : 155025056 + 47],
            "CACACACAGAGAAAAAAACATATATGCCCTTGTGATCTCATTTGTTT",
        )
        self.assertEqual(
            alignment[5],
            "AAACAAATGAGATCACAAGG-------GCATATA-----TGT-TTTTTTCTCTGTGT----------GTG",
        )
        self.assertEqual(alignment.sequences[5].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[5].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[5].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[5].annotations["rightCount"], 15)
        self.assertEqual(alignment.sequences[6].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[6].seq), 173908612)
        self.assertEqual(
            alignment.sequences[6].seq[157515110 : 157515110 + 47],
            "CACACACAGAGAAAAAAACATATATGCCCTTGTGATCTCATTTGTTT",
        )
        self.assertEqual(
            alignment[6],
            "AAACAAATGAGATCACAAGG-------GCATATA-----TGT-TTTTTTCTCTGTGT----------GTG",
        )
        self.assertEqual(
            alignment.sequences[6].annotations["quality"],
            "99999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[6].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[6].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[6].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[6].annotations["rightCount"], 15)
        self.assertEqual(alignment.sequences[7].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[7].seq), 174210431)
        self.assertEqual(
            alignment.sequences[7].seq[158037277 : 158037277 + 47],
            "CACACACAGAGAAAAAAATACATATGCCCTTGTGATCTCATTTGTTT",
        )
        self.assertEqual(
            alignment[7],
            "AAACAAATGAGATCACAAGG-------GCATATG-----TAT-TTTTTTCTCTGTGT----------GTG",
        )
        self.assertEqual(alignment.sequences[7].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[7].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[7].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[7].annotations["rightCount"], 15)
        self.assertEqual(alignment.sequences[8].id, "eriEur1.scaffold_266115")
        self.assertEqual(len(alignment.sequences[8].seq), 4589)
        self.assertEqual(
            alignment.sequences[8].seq[358 : 358 + 65],
            "taaataataagataagatgaaagCATAGCATGTATTTTCTtgccctctccttctctgtctctgtc",
        )
        self.assertEqual(
            alignment[8],
            "taaataataagataagatgaaagCATAGCATGTA-----TTTTCTtgccctctccttctctgtctctgtc",
        )
        self.assertEqual(
            alignment.sequences[8].annotations["quality"],
            "99999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[8].annotations["leftStatus"], "N")
        self.assertEqual(alignment.sequences[8].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[8].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[8].annotations["rightCount"], 9)
        self.assertEqual(alignment.sequences[9].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[9].seq), 125616256)
        self.assertEqual(
            alignment.sequences[9].seq[47544326 : 47544326 + 51],
            "CACAGACCCCACACAGAGAGAAAGCATATACACACGGTATCTCATTTTCTa",
        )
        self.assertEqual(
            alignment[9],
            "tAGAAAATGAGATACC-----------GTGTGTA-----TATGCTTTCTCTCTGTGT---GGGGTCTGTG",
        )
        self.assertEqual(
            alignment.sequences[9].annotations["quality"],
            "999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[9].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[9].annotations["leftCount"], 196)
        self.assertEqual(alignment.sequences[9].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[9].annotations["rightCount"], 9)
        self.assertEqual(alignment.sequences[10].id, "felCat3.scaffold_205680")
        self.assertEqual(len(alignment.sequences[10].seq), 119354)
        self.assertEqual(
            alignment.sequences[10].seq[45904 : 45904 + 51],
            "CACAGACACCACACAGAGAGGAAGAATACACGTGTCTTATCTCATTTGCTT",
        )
        self.assertEqual(
            alignment[10],
            "AAGCAAATGAGATAAG-----------ACACGTG-----TATTCTTCCTCTCTGTGT---GGTGTCTGTG",
        )
        self.assertEqual(
            alignment.sequences[10].annotations["quality"],
            "975559665645435353463242434515353544222333635339999",
        )
        self.assertEqual(alignment.sequences[10].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[10].annotations["leftCount"], 6)
        self.assertEqual(alignment.sequences[10].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[10].annotations["rightCount"], 9)
        self.assertEqual(alignment.sequences[11].id, "dasNov1.scaffold_56749")
        self.assertEqual(len(alignment.sequences[11].seq), 10470)
        self.assertEqual(
            alignment.sequences[11].seq[7165 : 7165 + 46],
            "TACACAGAGAAAGAATATTTATTTGCCCTGGTCCTCTCACTTTCCT",
        )
        self.assertEqual(
            alignment[11],
            "AGGAAAGTGAGAGGACCAGG-------GCAAATA-----AATATTCTTTCTCTGTGT----------A--",
        )
        self.assertEqual(
            alignment.sequences[11].annotations["quality"],
            "5556976455665765856867685558586864578595356565",
        )
        self.assertEqual(alignment.sequences[11].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[11].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[11].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[11].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[12].id, "echTel1.scaffold_288249")
        self.assertEqual(len(alignment.sequences[12].seq), 100002)
        self.assertEqual(
            alignment.sequences[12].seq[95246 : 95246 + 51],
            "AAGAAGGTGAGATGACAAGGGTGTATAGATAGGATATTCTTGCTTTGGGTG",
        )
        self.assertEqual(
            alignment[12],
            "AAGAAGGTGAGATGACAAGG-------GTGTATAGATAGGATATTCTTGCTTTGGGT----------G--",
        )
        self.assertEqual(
            alignment.sequences[12].annotations["quality"],
            "999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[12].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[12].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[12].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[12].annotations["rightCount"], 0)
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ornAna1.chr2")
        self.assertEqual(len(record.seq), 54797317)
        self.assertEqual(segment, (40046122, 40040432))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(record.seq), 498454)
        self.assertEqual(segment, (324116, 323493))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 13)
        self.assertEqual(len(alignment.annotations["empty"]), 2)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
            [[  3021136,   3021138,   3021139,   3021152,   3021153,   3021154,
                3021156,   3021156,   3021158,   3021163,   3021163,   3021165,
                3021166,   3021166,   3021166,   3021177,   3021177,   3021177,
                3021178,   3021180],
             [      431,       429,       429,       416,       415,       414,
                    412,       412,       410,       405,       405,       403,
                    402,       402,       402,       391,       391,       391,
                    390,       388],
             [     4022,      4020,      4019,      4006,      4005,      4004,
                   4002,      4002,      4002,      3997,      3997,      3995,
                   3994,      3994,      3994,      3983,      3983,      3983,
                   3982,      3980],
             [   171187,    171185,    171184,    171171,    171170,    171170,
                 171168,    171168,    171166,    171161,    171161,    171159,
                 171159,    171159,    171156,    171145,    171145,    171145,
                 171144,    171142],
             [    10722,     10724,     10725,     10738,     10739,     10739,
                  10739,     10739,     10741,     10746,     10746,     10748,
                  10749,     10750,     10753,     10764,     10764,     10764,
                  10765,     10767],
             [155025103, 155025101, 155025100, 155025087, 155025086, 155025085,
              155025083, 155025083, 155025081, 155025076, 155025076, 155025074,
              155025073, 155025073, 155025070, 155025059, 155025059, 155025059,
              155025058, 155025056],
             [157515157, 157515155, 157515154, 157515141, 157515140, 157515139,
              157515137, 157515137, 157515135, 157515130, 157515130, 157515128,
              157515127, 157515127, 157515124, 157515113, 157515113, 157515113,
              157515112, 157515110],
             [158037324, 158037322, 158037321, 158037308, 158037307, 158037306,
              158037304, 158037304, 158037302, 158037297, 158037297, 158037295,
              158037294, 158037294, 158037291, 158037280, 158037280, 158037280,
              158037279, 158037277],
             [      358,       360,       361,       374,       375,       376,
                    378,       385,       387,       392,       392,       394,
                    395,       396,       399,       410,       413,       420,
                    421,       423],
             [ 47544377,  47544375,  47544374,  47544361,  47544361,  47544361,
               47544361,  47544361,  47544359,  47544354,  47544354,  47544352,
               47544351,  47544350,  47544347,  47544336,  47544336,  47544329,
               47544328,  47544326],
             [    45955,     45953,     45952,     45939,     45939,     45939,
                  45939,     45939,     45937,     45932,     45932,     45930,
                  45929,     45928,     45925,     45914,     45914,     45907,
                  45906,     45904],
             [     7211,      7209,      7208,      7195,      7194,      7193,
                   7191,      7191,      7189,      7184,      7184,      7182,
                   7181,      7180,      7177,      7166,      7166,      7166,
                   7165,      7165],
             [    95246,     95248,     95249,     95262,     95263,     95264,
                  95266,     95266,     95268,     95273,     95278,     95280,
                  95281,     95282,     95285,     95296,     95296,     95296,
                  95297,     95297],
            ]
                    # fmt: on
                ),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 32886)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3021180 : 3021180 + 24],
            "TCCCAGAGAGTCTGATAGGAGGAG",
        )
        self.assertEqual(
            alignment[0], "-------------------TCCC-------AGAGAGTCTGA-TAGGAGGAG"
        )
        self.assertEqual(alignment.sequences[1].id, "cavPor2.scaffold_216473")
        self.assertEqual(len(alignment.sequences[1].seq), 10026)
        self.assertEqual(
            alignment.sequences[1].seq[356 : 356 + 32],
            "TACTTTATACTCAGAGCCACTATACAAaggca",
        )
        self.assertEqual(
            alignment[1], "-------------------tgcctTTGTATAGTGGCTCTGAGTATAAAGTA"
        )
        self.assertEqual(
            alignment.sequences[1].annotations["quality"],
            "67576649966655666885655548785776",
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[2].id, "oryCun1.scaffold_156751")
        self.assertEqual(len(alignment.sequences[2].seq), 4726)
        self.assertEqual(
            alignment.sequences[2].seq[3953 : 3953 + 27], "TATATTCAGAAATGCTATACACAGGCA"
        )
        self.assertEqual(
            alignment[2], "-------------------TGCCTGTGTATAGCATTTCTGAATATA-----"
        )
        self.assertEqual(
            alignment.sequences[2].annotations["quality"], "999999999999999999999999999"
        )
        self.assertEqual(alignment.sequences[2].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[2].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[3].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[3].seq), 174210431)
        self.assertEqual(
            alignment.sequences[3].seq[158037237 : 158037237 + 25],
            "TACTTCATATTCATACACACTCAGA",
        )
        self.assertEqual(
            alignment[3], "-----------------------TCTG---AGTGTGTATGAATATGAAGTA"
        )
        self.assertEqual(alignment.sequences[3].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[3].annotations["leftCount"], 15)
        self.assertEqual(alignment.sequences[3].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[3].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[4].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[4].seq), 173908612)
        self.assertEqual(
            alignment.sequences[4].seq[157515070 : 157515070 + 25],
            "TACTTCATATTCATACATACTCAGA",
        )
        self.assertEqual(
            alignment[4], "-----------------------TCTG---AGTATGTATGAATATGAAGTA"
        )
        self.assertEqual(
            alignment.sequences[4].annotations["quality"], "9999999999999999999999999"
        )
        self.assertEqual(alignment.sequences[4].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[4].annotations["leftCount"], 15)
        self.assertEqual(alignment.sequences[4].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[4].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[5].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[5].seq), 170899992)
        self.assertEqual(
            alignment.sequences[5].seq[155025016 : 155025016 + 25],
            "TACTTCATATTCATACATACTCAGA",
        )
        self.assertEqual(
            alignment[5], "-----------------------TCTG---AGTATGTATGAATATGAAGTA"
        )
        self.assertEqual(alignment.sequences[5].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[5].annotations["leftCount"], 15)
        self.assertEqual(alignment.sequences[5].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[5].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[6].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[6].seq), 133105)
        self.assertEqual(
            alignment.sequences[6].seq[10782 : 10782 + 25], "TCTGAGTATGTCTGAATATGAAGTG"
        )
        self.assertEqual(
            alignment[6], "-----------------------TCTG---AGTATGTCTGAATATGAAGTG"
        )
        self.assertEqual(alignment.sequences[6].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[6].annotations["leftCount"], 15)
        self.assertEqual(alignment.sequences[6].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[6].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[7].id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(alignment.sequences[7].seq), 359464)
        self.assertEqual(
            alignment.sequences[7].seq[171102 : 171102 + 25],
            "CACATCATGTTCAGACACACTTAGA",
        )
        self.assertEqual(
            alignment[7], "-----------------------TCTA---AGTGTGTCTGAACATGATGTG"
        )
        self.assertEqual(
            alignment.sequences[7].annotations["quality"], "9999999999999999999999999"
        )
        self.assertEqual(alignment.sequences[7].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[7].annotations["leftCount"], 15)
        self.assertEqual(alignment.sequences[7].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[7].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[8].id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(alignment.sequences[8].seq), 498454)
        self.assertEqual(
            alignment.sequences[8].seq[323468 : 323468 + 25],
            "TACCTCAAGTTCAGACACTCAGAGA",
        )
        self.assertEqual(
            alignment[8], "-----------------------TCTC---TGAGTGTCTGAACTTGAGGTA"
        )
        self.assertEqual(
            alignment.sequences[8].annotations["quality"], "9999999999999999999999999"
        )
        self.assertEqual(alignment.sequences[8].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[8].annotations["leftCount"], 623)
        self.assertEqual(alignment.sequences[8].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[8].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[9].id, "eriEur1.scaffold_266115")
        self.assertEqual(len(alignment.sequences[9].seq), 4589)
        self.assertEqual(
            alignment.sequences[9].seq[432 : 432 + 19], "TCTCAGTGTGTCTGACCAG"
        )
        self.assertEqual(
            alignment[9], "-----------------------TCTC---AGTGTGTCTGACCAG------"
        )
        self.assertEqual(
            alignment.sequences[9].annotations["quality"], "9999999999999999999"
        )
        self.assertEqual(alignment.sequences[9].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[9].annotations["leftCount"], 9)
        self.assertEqual(alignment.sequences[9].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[9].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[10].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[10].seq), 125616256)
        self.assertEqual(
            alignment.sequences[10].seq[47544292 : 47544292 + 25],
            "TACCTCACATTCAGACACACTCAGA",
        )
        self.assertEqual(
            alignment[10], "-----------------------TCTG---AGTGTGTCTGAATGTGAGGTA"
        )
        self.assertEqual(
            alignment.sequences[10].annotations["quality"], "9999999999999999999999999"
        )
        self.assertEqual(alignment.sequences[10].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[10].annotations["leftCount"], 9)
        self.assertEqual(alignment.sequences[10].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[10].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[11].id, "felCat3.scaffold_205680")
        self.assertEqual(len(alignment.sequences[11].seq), 119354)
        self.assertEqual(
            alignment.sequences[11].seq[45870 : 45870 + 25], "CACCTCACATTCAGGAACACTCAGA"
        )
        self.assertEqual(
            alignment[11], "-----------------------TCTG---AGTGTTCCTGAATGTGAGGTG"
        )
        self.assertEqual(
            alignment.sequences[11].annotations["quality"], "9999999999999999999999999"
        )
        self.assertEqual(alignment.sequences[11].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[11].annotations["leftCount"], 9)
        self.assertEqual(alignment.sequences[11].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[11].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[12].id, "dasNov1.scaffold_56749")
        self.assertEqual(len(alignment.sequences[12].seq), 10470)
        self.assertEqual(
            alignment.sequences[12].seq[7118 : 7118 + 47],
            "TACCTCAAGTTCAGACACGCTCAGAAATGCTCCGCAAAGGCACACAA",
        )
        self.assertEqual(
            alignment[12], "T-TGTGTGCCTTTGCGGAGCATTTCTG---AGCGTGTCTGAACTTGAGGTA"
        )
        self.assertEqual(
            alignment.sequences[12].annotations["quality"],
            "73659557766555777595699547965955937797755656587",
        )
        self.assertEqual(alignment.sequences[12].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[12].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[12].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[12].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[13].id, "echTel1.scaffold_288249")
        self.assertEqual(len(alignment.sequences[13].seq), 100002)
        self.assertEqual(
            alignment.sequences[13].seq[95297 : 95297 + 48],
            "TGCCCAGGACTGCGCATGGTATTTCTTGGTGTGTCTGAAGGTGAGATA",
        )
        self.assertEqual(
            alignment[13], "TGCCCAGGACTGCGCATGGTATTTCTT---GGTGTGTCTGAAGGTGAGATA"
        )
        self.assertEqual(
            alignment.sequences[13].annotations["quality"],
            "999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[13].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[13].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[13].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[13].annotations["rightCount"], 0)
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ornAna1.chr2")
        self.assertEqual(len(record.seq), 54797317)
        self.assertEqual(segment, (40046122, 40040432))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 14)
        self.assertEqual(len(alignment.annotations["empty"]), 1)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
            [[  3021180,   3021180,   3021180,   3021180,   3021184,   3021184,
                3021184,   3021195,   3021195,   3021198,   3021199,   3021204],
             [      388,       388,       388,       388,       384,       380,
                    377,       366,       365,       362,       361,       356],
             [     3980,      3980,      3980,      3980,      3976,      3972,
                   3969,      3958,      3957,      3954,      3953,      3953],
             [158037262, 158037262, 158037262, 158037262, 158037262, 158037258,
              158037258, 158037247, 158037246, 158037243, 158037242, 158037237],
             [157515095, 157515095, 157515095, 157515095, 157515095, 157515091,
              157515091, 157515080, 157515079, 157515076, 157515075, 157515070],
             [155025041, 155025041, 155025041, 155025041, 155025041, 155025037,
              155025037, 155025026, 155025025, 155025022, 155025021, 155025016],
             [    10782,     10782,     10782,     10782,     10782,     10786,
                  10786,     10797,     10798,     10801,     10802,     10807],
             [   171127,    171127,    171127,    171127,    171127,    171123,
                 171123,    171112,    171111,    171108,    171107,    171102],
             [   323493,    323493,    323493,    323493,    323493,    323489,
                 323489,    323478,    323477,    323474,    323473,    323468],
             [      432,       432,       432,       432,       432,       436,
                    436,       447,       448,       451,       451,       451],
             [ 47544317,  47544317,  47544317,  47544317,  47544317,  47544313,
               47544313,  47544302,  47544301,  47544298,  47544297,  47544292],
             [    45895,     45895,     45895,     45895,     45895,     45891,
                  45891,     45880,     45879,     45876,     45875,     45870],
             [     7165,      7164,      7164,      7147,      7143,      7139,
                   7139,      7128,      7127,      7124,      7123,      7118],
             [    95297,     95298,     95299,     95316,     95320,     95324,
                  95324,     95335,     95336,     95339,     95340,     95345],
            ]
                    # fmt: on
                ),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 309116)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3021204 : 3021204 + 71],
            "TGCACTGGTTTTCCTGCAGTGGTTCTCAGTAATAGGAAGACAACAGAATTTGAAGTATCCGGCTTTGGCCA",
        )
        self.assertEqual(
            alignment[0],
            "TGCACTGGTTTTCC-TGCAGTGGTTCTCAGTAATAGGAAGACA-ACAGAATTTGAAGTATCCGGCTTTGGCCA",
        )
        self.assertEqual(alignment.sequences[1].id, "echTel1.scaffold_288249")
        self.assertEqual(len(alignment.sequences[1].seq), 100002)
        self.assertEqual(
            alignment.sequences[1].seq[95345 : 95345 + 73],
            "GGCATTGGTTTTTAGAGAGAGAACCCACATAAGTAGGAAAACATTTTGAATTTATAGTAAATATTCTTGGCTA",
        )
        self.assertEqual(
            alignment[1],
            "GGCATTGGTTTTTAGAGAGAGAACCCACATAAGTAGGAAAACATTTTGAATTTATAGTAAATATTCTTGGCTA",
        )
        self.assertEqual(
            alignment.sequences[1].annotations["quality"],
            "9999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[2].id, "dasNov1.scaffold_56749")
        self.assertEqual(len(alignment.sequences[2].seq), 10470)
        self.assertEqual(
            alignment.sequences[2].seq[7048 : 7048 + 70],
            "TAGCCAAGAGCATTTATTATAAATTCAAAATGCTCTCTTAGTGTGATTTCCCTCACTGAAAATCAATGCA",
        )
        self.assertEqual(
            alignment[2],
            "TGCATTGATTTTCAGTGAGGGAAATCACAC---TAAGAGAGCATTTTGAATTTATAATAAATGCTCTTGGCTA",
        )
        self.assertEqual(
            alignment.sequences[2].annotations["quality"],
            "5796635799979575966667835948658458696597898258979678997999677999997699",
        )
        self.assertEqual(alignment.sequences[2].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[2].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[3].id, "felCat3.scaffold_205680")
        self.assertEqual(len(alignment.sequences[3].seq), 119354)
        self.assertEqual(
            alignment.sequences[3].seq[45798 : 45798 + 72],
            "CAGCCAAAAACATTTGCTATAATTTCAAAATGCATTCCTATTCATGTGAATCAATGCATGAAAATCAATGCA",
        )
        self.assertEqual(
            alignment[3],
            "TGCATTGATTTTCA-TGCATTGATTCACATGAATAGGAATGCATTTTGAAATTATAGCAAATGTTTTTGGCTG",
        )
        self.assertEqual(
            alignment.sequences[3].annotations["quality"],
            "999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[3].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[3].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[3].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[3].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[4].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[4].seq), 125616256)
        self.assertEqual(
            alignment.sequences[4].seq[47544221 : 47544221 + 71],
            "CAGCCAAAAGCATTTCCTATAAATTCAAATGTATTCCTATTCATGTGAATCAACGCATGAAAATTAATGCA",
        )
        self.assertEqual(
            alignment[4],
            "TGCATTAATTTTCA-TGCGTTGATTCACATGAATAGGAATACA-TTTGAATTTATAGGAAATGCTTTTGGCTG",
        )
        self.assertEqual(
            alignment.sequences[4].annotations["quality"],
            "99999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[4].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[4].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[4].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[4].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[5].id, "eriEur1.scaffold_266115")
        self.assertEqual(len(alignment.sequences[5].seq), 4589)
        self.assertEqual(
            alignment.sequences[5].seq[451 : 451 + 70],
            "GCACTGATTCTCGAGGGTTGATTCCCAGTAAGAGGAAACTGCGTGAGTTTACAGTACATGGGCTTGGCTG",
        )
        self.assertEqual(
            alignment[5],
            "-GCACTGATTCTCG-AGGGTTGATTCCCAGTAAGAGGAAACTG-CGTGAGTTTACAGTACATGGGCTTGGCTG",
        )
        self.assertEqual(
            alignment.sequences[5].annotations["quality"],
            "9999999999999999999999999999999999999999999999999999999999999999999799",
        )
        self.assertEqual(alignment.sequences[5].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[5].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[5].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[5].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[6].id, "sorAra1.scaffold_2476")
        self.assertEqual(len(alignment.sequences[6].seq), 4997)
        self.assertEqual(
            alignment.sequences[6].seq[1615 : 1615 + 70],
            "TGGCCAAATCATTTATCATAAATTGAAGTGCTTTCCTTTTCATGTGAATCAATACTTAAAAATCAATGCA",
        )
        self.assertEqual(
            alignment[6],
            "TGCATTGATTTTTA-AGTATTGATTCACATGAAAAGGAAAGCA-CTTCAATTTATGATAAAT-GATTTGGCCA",
        )
        self.assertEqual(
            alignment.sequences[6].annotations["quality"],
            "8999999999998999999999999999999999999999999999999997999999999999999999",
        )
        self.assertEqual(alignment.sequences[6].annotations["leftStatus"], "N")
        self.assertEqual(alignment.sequences[6].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[6].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[6].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[7].id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(alignment.sequences[7].seq), 498454)
        self.assertEqual(
            alignment.sequences[7].seq[323397 : 323397 + 71],
            "TGGTCAGAAGCATTTCCTATCAATTAGAATCGTTTCCTATTAATGTGACATAACACATGGAATTCAATGCA",
        )
        self.assertEqual(
            alignment[7],
            "TGCATTGAATTCCA-TGTGTTATGTCACATTAATAGGAAACGA-TTCTAATTGATAGGAAATGCTTCTGACCA",
        )
        self.assertEqual(
            alignment.sequences[7].annotations["quality"],
            "99999899899999999999999999999999999999999999999999999999999999999999997",
        )
        self.assertEqual(alignment.sequences[7].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[7].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[7].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[7].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[8].id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(alignment.sequences[8].seq), 359464)
        self.assertEqual(
            alignment.sequences[8].seq[171035 : 171035 + 67],
            "TGGCCAAATGCCTTTCCTACAAATTCACATGTTTTCCCGTTAAGGTGAATCAATGTGTAAAAATCCA",
        )
        self.assertEqual(
            alignment[8],
            "TG----GATTTTTA-CACATTGATTCACCTTAACGGGAAAACA-TGTGAATTTGTAGGAAAGGCATTTGGCCA",
        )
        self.assertEqual(
            alignment.sequences[8].annotations["quality"],
            "8899999999989999989994384888899999999966359999569699999923799351281",
        )
        self.assertEqual(alignment.sequences[8].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[8].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[8].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[8].annotations["rightCount"], 6280)
        self.assertEqual(alignment.sequences[9].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[9].seq), 133105)
        self.assertEqual(
            alignment.sequences[9].seq[10807 : 10807 + 71],
            "TGCATTGGTTTTTATGCCTTGATTCACATGAATAGGAAAACGTTTGAATTTATAGGAAATGGTTTTGGCCA",
        )
        self.assertEqual(
            alignment[9],
            "TGCATTGGTTTTTA-TGCCTTGATTCACATGAATAGGAAAACG-TTTGAATTTATAGGAAATGGTTTTGGCCA",
        )
        self.assertEqual(alignment.sequences[9].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[9].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[9].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[9].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[10].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[10].seq), 170899992)
        self.assertEqual(
            alignment.sequences[10].seq[155024945 : 155024945 + 71],
            "TGGCCAAAACCATTTACTATAAATTCAAACATTTTCCTATTCATGTAAATCAATGCATAAAAATCAAAGCA",
        )
        self.assertEqual(
            alignment[10],
            "TGCTTTGATTTTTA-TGCATTGATTTACATGAATAGGAAAATG-TTTGAATTTATAGTAAATGGTTTTGGCCA",
        )
        self.assertEqual(alignment.sequences[10].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[10].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[10].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[10].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[11].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[11].seq), 173908612)
        self.assertEqual(
            alignment.sequences[11].seq[157514999 : 157514999 + 71],
            "TGGCCAAAACCATTTACTATAAATTCAAATGTTTTCCTATTCATGTAAATCAATGCATAAAAATCAAAGCA",
        )
        self.assertEqual(
            alignment[11],
            "TGCTTTGATTTTTA-TGCATTGATTTACATGAATAGGAAAACA-TTTGAATTTATAGTAAATGGTTTTGGCCA",
        )
        self.assertEqual(
            alignment.sequences[11].annotations["quality"],
            "99999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[11].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[11].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[11].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[11].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[12].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[12].seq), 174210431)
        self.assertEqual(
            alignment.sequences[12].seq[158037166 : 158037166 + 71],
            "TGGCCAAAACCATTTACTATAAATTCAAACGTTTTCCTATCCATGTAAATCAATGCATAAAAATCAAAGCA",
        )
        self.assertEqual(
            alignment[12],
            "TGCTTTGATTTTTA-TGCATTGATTTACATGGATAGGAAAACG-TTTGAATTTATAGTAAATGGTTTTGGCCA",
        )
        self.assertEqual(alignment.sequences[12].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[12].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[12].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[12].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[13].id, "oryCun1.scaffold_156751")
        self.assertEqual(len(alignment.sequences[13].seq), 4726)
        self.assertEqual(
            alignment.sequences[13].seq[3888 : 3888 + 65],
            "TGACCAAAAGCCATTGCTATAAATGCAAATGTTTTCCCGTCGATGTAAATCAAAGTATGTAAAGA",
        )
        self.assertEqual(
            alignment[13],
            "------TCTTTACA-TACTTTGATTTACATCGACGGGAAAACA-TTTGCATTTATAGCAATGGCTTTTGGTCA",
        )
        self.assertEqual(
            alignment.sequences[13].annotations["quality"],
            "99999994999999999999989999999999999999999999999999899999999999999",
        )
        self.assertEqual(alignment.sequences[13].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[13].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[13].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[13].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[14].id, "cavPor2.scaffold_216473")
        self.assertEqual(len(alignment.sequences[14].seq), 10026)
        self.assertEqual(
            alignment.sequences[14].seq[285 : 285 + 71],
            "TGGCCAAAACCATTTACTATAAATTGGACTGTTTTCCTATTAATGTGAATTAATGCATGGAAATCAATGCC",
        )
        self.assertEqual(
            alignment[14],
            "GGCATTGATTTCCA-TGCATTAATTCACATTAATAGGAAAACA-GTCCAATTTATAGTAAATGGTTTTGGCCA",
        )
        self.assertEqual(
            alignment.sequences[14].annotations["quality"],
            "77888768786695388675879644655668865666547868687676669688688666687574686",
        )
        self.assertEqual(alignment.sequences[14].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[14].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[14].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[14].annotations["rightCount"], 0)
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "ornAna1.chr2")
        self.assertEqual(len(record.seq), 54797317)
        self.assertEqual(segment, (40046122, 40040432))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 15)
        self.assertEqual(len(alignment.annotations["empty"]), 1)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
            [[  3021204,   3021205,   3021206,   3021210,   3021218,   3021218,
                3021233,   3021236,   3021246,   3021246,   3021264,   3021265,
                3021275],
             [    95345,     95346,     95347,     95351,     95359,     95360,
                  95375,     95378,     95388,     95389,     95407,     95408,
                  95418],
             [     7118,      7117,      7116,      7112,      7104,      7103,
                   7088,      7088,      7078,      7077,      7059,      7058,
                   7048],
             [    45870,     45869,     45868,     45864,     45856,     45856,
                  45841,     45838,     45828,     45827,     45809,     45808,
                  45798],
             [ 47544292,  47544291,  47544290,  47544286,  47544278,  47544278,
               47544263,  47544260,  47544250,  47544250,  47544232,  47544231,
               47544221],
             [      451,       451,       452,       456,       464,       464,
                    479,       482,       492,       492,       510,       511,
                    521],
             [     1685,      1684,      1683,      1679,      1671,      1671,
                   1656,      1653,      1643,      1643,      1625,      1625,
                   1615],
             [   323468,    323467,    323466,    323462,    323454,    323454,
                 323439,    323436,    323426,    323426,    323408,    323407,
                 323397],
             [   171102,    171101,    171100,    171100,    171092,    171092,
                 171077,    171074,    171064,    171064,    171046,    171045,
                 171035],
             [    10807,     10808,     10809,     10813,     10821,     10821,
                  10836,     10839,     10849,     10849,     10867,     10868,
                  10878],
             [155025016, 155025015, 155025014, 155025010, 155025002, 155025002,
              155024987, 155024984, 155024974, 155024974, 155024956, 155024955,
              155024945],
             [157515070, 157515069, 157515068, 157515064, 157515056, 157515056,
              157515041, 157515038, 157515028, 157515028, 157515010, 157515009,
              157514999],
             [158037237, 158037236, 158037235, 158037231, 158037223, 158037223,
              158037208, 158037205, 158037195, 158037195, 158037177, 158037176,
              158037166],
             [     3953,      3953,      3953,      3953,      3945,      3945,
                   3930,      3927,      3917,      3917,      3899,      3898,
                   3888],
             [      356,       355,       354,       350,       342,       342,
                    327,       324,       314,       314,       296,       295,
                    285],
            ]
                    # fmt: on
                ),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 891219)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3021275 : 3021275 + 146],
            "CTTTCTCTGACATTACTGTAACTGAAGTAGCTCAGAAGCACAAAAGGTCACATCATGCATCCATGCAGAATCCACTGAAGCTGTTTGGAAAGGCCACGTGTCTTCCCAGAAGGCCAGTTACACCATCATTTCCTTCCATGTTTCAG",
        )
        self.assertEqual(
            alignment[0],
            "CTTTCTCTGACATTACTGTAACTGAAGTAGCTC-AGAAGCACAAAAGGTCACATCATGCATCCATGCAGAATCCACTGAAGCTGTTTGGAAAGGC-----------------------CACGTGTCTTCCCAGAAGGCCAGTTACACCATCATTTCCTTCCATGTTTCAG",
        )
        self.assertEqual(alignment.sequences[1].id, "cavPor2.scaffold_216473")
        self.assertEqual(len(alignment.sequences[1].seq), 10026)
        self.assertEqual(
            alignment.sequences[1].seq[116 : 116 + 169],
            "TTGAAATATGGAAGGAAATGATGGCGTGACTGGATTTCCTGCGAGACACATGAACTAAGTAATAAAATGCGAGGTGCCTTTATGAACGGAACTAGTAAATTTTGCATAAATGCCTGACATGACCTTTTGTGCTTTTCAGCTAGTTCGCTTACAGTAACATCAGCAAAGG",
        )
        self.assertEqual(
            alignment[1],
            "CCTTTGCTGATGTTACTGTAAGCGAACTAGCTG-AAAAGCACAAAAGGTCATGTCAGGCATTTATGCAAAATTTACTAGTTCCGTTCATAAAGGCACCTCGCATTTTATTACTTAGTTCATGTGTCTCGCAGGAAATCCAGTCACGCCATCATTTCCTTCCATATTTCAA",
        )
        self.assertEqual(
            alignment.sequences[1].annotations["quality"],
            "6685365455476645398666666666666688787886799997788875886655536666688786768476688854689668876558688874778656675665668746666788474786888875667457786687996867886548798845555",
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[2].id, "oryCun1.scaffold_156751")
        self.assertEqual(len(alignment.sequences[2].seq), 4726)
        self.assertEqual(
            alignment.sequences[2].seq[3730 : 3730 + 158],
            "TTGAAATGTGGAAGGAAATCGCGGCGTGACTGGATTTCCTGAACAAAGTAATACAACACAAGGTGCTGTTATAAACAGTGCTAGTACATTTGACATAAATGCCTGCCATGGCCTTTTGTGCTTTGGAGCGCTTTCAGCTAGAGTAACATCAGAGACGA",
        )
        self.assertEqual(
            alignment[2],
            "TCGTCTCTGATGTTACTCTAGCTGAAAGCGCTC-CAAAGCACAAAAGGCCATGGCAGGCATTTATGTCAAATGTACTAGCACTGTTTATAACAGCACCTTGTGTTGTATTACTTTGTTCA-----------GGAAATCCAGTCACGCCGCGATTTCCTTCCACATTTCAA",
        )
        self.assertEqual(
            alignment.sequences[2].annotations["quality"],
            "99999987999999999999999999999799999999799999999898999999988999999988999999989996999999999988897996999999999999979999987899998987999799998998999999999777999899",
        )
        self.assertEqual(alignment.sequences[2].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[2].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[3].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[3].seq), 174210431)
        self.assertEqual(
            alignment.sequences[3].seq[158036997 : 158036997 + 169],
            "ATGAAATATGGAAGGAAATGATGGCGTGACTGAATTTCCTGAAAGATACATTAACAAAGTAATAAAACACAAGGTACCCTTATAAACAGCACTAGTAAGTTTTACATGAATGCCTCCCATGACCTTTTGTGCTTTTGAGCTATTTCAGTTACAGTAACATCAGAGAAGG",
        )
        self.assertEqual(
            alignment[3],
            "CCTTCTCTGATGTTACTGTAACTGAAATAGCTC-AAAAGCACAAAAGGTCATGGGAGGCATTCATGTAAAACTTACTAGTGCTGTTTATAAGGGTACCTTGTGTTTTATTACTTTGTTAATGTATCTTTCAGGAAATTCAGTCACGCCATCATTTCCTTCCATATTTCAT",
        )
        self.assertEqual(alignment.sequences[3].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[3].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[3].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[3].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[4].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[4].seq), 173908612)
        self.assertEqual(
            alignment.sequences[4].seq[157514830 : 157514830 + 169],
            "ATGAAATATGGAAGGAAATGATGGCATGACTGGATTTCCTGAAAGGTACATTAACAAAGTAATAAAACACAAGGTACCCTTATAAACAGCACTAGTAAGTTTTACATAAATGCCTCCCATGACCTTTTGTGCTTTTGAGCTATTTCAGTTACAGTAACATCAGAGAAGG",
        )
        self.assertEqual(
            alignment[4],
            "CCTTCTCTGATGTTACTGTAACTGAAATAGCTC-AAAAGCACAAAAGGTCATGGGAGGCATTTATGTAAAACTTACTAGTGCTGTTTATAAGGGTACCTTGTGTTTTATTACTTTGTTAATGTACCTTTCAGGAAATCCAGTCATGCCATCATTTCCTTCCATATTTCAT",
        )
        self.assertEqual(
            alignment.sequences[4].annotations["quality"],
            "9999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[4].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[4].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[4].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[4].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[5].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[5].seq), 170899992)
        self.assertEqual(
            alignment.sequences[5].seq[155024776 : 155024776 + 169],
            "ATGAAATATGGAAGGAAATGATGGCATGACTGGATTTCCTGAAAGGTACATTAACAAAGTAATAAAACACAAGATACCCTTATAAACAGCACTAGTAAGTTTTACATAAATGCCTCCCATGACCTTTTGTGCTTTTGAGCTATTTCAGTTACAGTAACATCAGAGAAGG",
        )
        self.assertEqual(
            alignment[5],
            "CCTTCTCTGATGTTACTGTAACTGAAATAGCTC-AAAAGCACAAAAGGTCATGGGAGGCATTTATGTAAAACTTACTAGTGCTGTTTATAAGGGTATCTTGTGTTTTATTACTTTGTTAATGTACCTTTCAGGAAATCCAGTCATGCCATCATTTCCTTCCATATTTCAT",
        )
        self.assertEqual(alignment.sequences[5].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[5].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[5].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[5].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[6].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[6].seq), 133105)
        self.assertEqual(
            alignment.sequences[6].seq[10878 : 10878 + 169],
            "CCTTCTCTGATATTACTGTAACTGAAATAGCTCAAAAGCACAAAAGGTCATGGCAGGCATTTATGTAAAACTTACTAGTGCTGTTTATAAAGGCACCTTGTGTTTTATTACTTTGTTAATGTACCTTTCAGGAAATCCAGTCATGCCATCATTTCCTTCCATATTTCAT",
        )
        self.assertEqual(
            alignment[6],
            "CCTTCTCTGATATTACTGTAACTGAAATAGCTC-AAAAGCACAAAAGGTCATGGCAGGCATTTATGTAAAACTTACTAGTGCTGTTTATAAAGGCACCTTGTGTTTTATTACTTTGTTAATGTACCTTTCAGGAAATCCAGTCATGCCATCATTTCCTTCCATATTTCAT",
        )
        self.assertEqual(alignment.sequences[6].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[6].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[6].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[6].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[7].id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(alignment.sequences[7].seq), 498454)
        self.assertEqual(
            alignment.sequences[7].seq[323241 : 323241 + 156],
            "CTGAGATATGGAAGGAAATGATGGCCCGACTGGATTTCCTGTAAGGCACATGAACAAAGTAATAAAACACGAGATGCCTTCAGAAACTTTACATCAATGCCTGGCATGACCTTTTGTGCTTCTGAGCTACTTTGGTTACAGGAACATCTGAGAGGA",
        )
        self.assertEqual(
            alignment[7],
            "TCCTCTCAGATGTTCCTGTAACCAAAGTAGCTC-AGAAGCACAAAAGGTCATGCCAGGCATTGATGTAAA-------------GTTTCTGAAGGCATCTCGTGTTTTATTACTTTGTTCATGTGCCTTACAGGAAATCCAGTCGGGCCATCATTTCCTTCCATATCTCAG",
        )
        self.assertEqual(
            alignment.sequences[7].annotations["quality"],
            "999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[7].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[7].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[7].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[7].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[8].id, "sorAra1.scaffold_2476")
        self.assertEqual(len(alignment.sequences[8].seq), 4997)
        self.assertEqual(
            alignment.sequences[8].seq[1445 : 1445 + 170],
            "CTGGAATATGGGAGGAAATGATGAGATTACAGGATTTCCGGAAAAGTTCATGATCAACAACATTATATGCAAGGGGTCTTTATAAACAGCACTGGTATATTTTACATCAATGCCTGACATGACCTTTTGTGCTTTTTGAGCTACCTCAGTTACAGTAACATCAGAGAAGG",
        )
        self.assertEqual(
            alignment[8],
            "CCTTCTCTGATGTTACTGTAACTGAGGTAGCTCAAAAAGCACAAAAGGTCATGTCAGGCATTGATGTAAAATATACCAGTGCTGTTTATAAAGACCCCTTGCATATAATGTTGTTGATCATGAACTTTTCCGGAAATCCTGTAATCTCATCATTTCCTCCCATATTCCAG",
        )
        self.assertEqual(
            alignment.sequences[8].annotations["quality"],
            "99999999999999997999999999999999999999999999999999999999999999999999999999999979999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[8].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[8].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[8].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[8].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[9].id, "eriEur1.scaffold_266115")
        self.assertEqual(len(alignment.sequences[9].seq), 4589)
        self.assertEqual(
            alignment.sequences[9].seq[521 : 521 + 156],
            "CCTTCTCTGATCTGACTGTAACCGAAGTCTCTCGAAAGCACAAAAGCTCACGGCAGGCCTTTCTGTAAAATACAGCGGCGCTGCTTCCAAAGGCaccttgcagatgtctctcacgcATTTGGCAGGAGTCCCTGTCACTCGGCCAGTTCCTTCCTG",
        )
        self.assertEqual(
            alignment[9],
            "CCTTCTCTGATCTGACTGTAACCGAAGTCTCTC-GAAAGCACAAAAGCTCACGGCAGGCCTTTCTGTAAAATACAGCGGCGCTGCTTCCAAAGGCaccttgca------gatgtctctcacgcATTTGGCAGGAGTCCCTGTCACTCGGCCAGTTCC-------TTCCTG",
        )
        self.assertEqual(
            alignment.sequences[9].annotations["quality"],
            "999999999999999999999999999999999999999999999999999999999799999999999999999999998987999899999999999999999999999999999999999999999999999999999999999998999989",
        )
        self.assertEqual(alignment.sequences[9].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[9].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[9].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[9].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[10].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[10].seq), 125616256)
        self.assertEqual(
            alignment.sequences[10].seq[47544053 : 47544053 + 168],
            "TTGAAATATGGAAGGAAGTAACAGGATGACAGGATTTCCTGTAAAGTACATGAACAAAGCAATATTACACGAGGTGCCTTTATAAACAGCTCTAGTAAATTTTACATAAAGCCTGACATGACCTTTTGTGCTTTTGAGCTACTTCAGTTACAGTAACATCAGAGACAG",
        )
        self.assertEqual(
            alignment[10],
            "CTGTCTCTGATGTTACTGTAACTGAAGTAGCTC-AAAAGCACAAAAGGTCATGTCAGGC-TTTATGTAAAATTTACTAGAGCTGTTTATAAAGGCACCTCGTGTAATATTGCTTTGTTCATGTACTTTACAGGAAATCCTGTCATCCTGTTACTTCCTTCCATATTTCAA",
        )
        self.assertEqual(
            alignment.sequences[10].annotations["quality"],
            "999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[10].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[10].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[10].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[10].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[11].id, "felCat3.scaffold_205680")
        self.assertEqual(len(alignment.sequences[11].seq), 119354)
        self.assertEqual(
            alignment.sequences[11].seq[45629 : 45629 + 169],
            "CTGAACTATGGAAGGAAATGACAGGATGACAGGATTTCCTATAAAGTGCATGAACAAAGCAATAATACACAAGGTACCTTTATAAACAGCACTAGTAAATTTTACATAAATGCCTGACATCACCTTTTGTGCTTTTGAACTATTTCAGTTACAGTAACATCAGAGACAG",
        )
        self.assertEqual(
            alignment[11],
            "CTGTCTCTGATGTTACTGTAACTGAAATAGTTC-AAAAGCACAAAAGGTGATGTCAGGCATTTATGTAAAATTTACTAGTGCTGTTTATAAAGGTACCTTGTGTATTATTGCTTTGTTCATGCACTTTATAGGAAATCCTGTCATCCTGTCATTTCCTTCCATAGTTCAG",
        )
        self.assertEqual(
            alignment.sequences[11].annotations["quality"],
            "9999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[11].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[11].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[11].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[11].annotations["rightCount"], 34165)
        self.assertEqual(alignment.sequences[12].id, "dasNov1.scaffold_56749")
        self.assertEqual(len(alignment.sequences[12].seq), 10470)
        self.assertEqual(
            alignment.sequences[12].seq[6880 : 6880 + 168],
            "GTGACTTATGGAAGGAAATGATGGTATGACAGGATTTCCTGTAAGGGGGATAGAAAAAGTAATAGAACACAAGGTCCTTTATAAACAGCACTAGTAAATTTTACATCAATGCCCGACATGACCTTTTGTGCTTTTGAGCCTCTTCAGTTACAGTAACATCAGAGAAGG",
        )
        self.assertEqual(
            alignment[12],
            "CCTTCTCTGATGTTACTGTAACTGAAGAGGCTC-AAAAGCACAAAAGGTCATGTCGGGCATTGATGTAAAATTTACTAGTGCTG-TTTATAAAGGACCTTGTGTTCTATTACTTTTTCTATCCCCCTTACAGGAAATCCTGTCATACCATCATTTCCTTCCATAAGTCAC",
        )
        self.assertEqual(
            alignment.sequences[12].annotations["quality"],
            "989999999999999937699999999999999799999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[12].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[12].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[12].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[12].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[13].id, "echTel1.scaffold_288249")
        self.assertEqual(len(alignment.sequences[13].seq), 100002)
        self.assertEqual(
            alignment.sequences[13].seq[95418 : 95418 + 169],
            "CCATCTCTGATATTACTGTAATGGAATTACTTTGAAAGCACAAAAGGTCAGGACAAGCGTGTATGTGAAATTTCCTAGAGCTGTTTTCCTGCACACCTTGGATTTTATTGCTTAGTTCATTTGCTTTCCCAGAAATCCCGCCATGCCATCATTTCCTTCCACATCTCAG",
        )
        self.assertEqual(
            alignment[13],
            "CCATCTCTGATATTACTGTAATGGAATTACTTT-GAAAGCACAAAAGGTCAGGACAAGCGTGTATGTGAAATTTCCTAGAGCTGTTTTCCTGCACACCTTGGATTTTATTGCTTAGTTCATTTGCTTTCCCAGAAATCCCGCCATGCCATCATTTCCTTCCACATCTCAG",
        )
        self.assertEqual(
            alignment.sequences[13].annotations["quality"],
            "9999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[13].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[13].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[13].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[13].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[14].id, "ornAna1.chr2")
        self.assertEqual(len(alignment.sequences[14].seq), 54797317)
        self.assertEqual(
            alignment.sequences[14].seq[40040263 : 40040263 + 169],
            "TGAAGATATGGAAGGAAATGATGGCTTGACAGGATTTCCCCTAAGGTGAATGAATAGCCTCATAAAACACAAAGCAGCTTTATGAACAGGGCTATTAAGTTGCACATGATTGGCTGATATGACCTTTCCTATCTTCTGGGTAGCTGAGTTACAGTGATATCAGAGGTGG",
        )
        self.assertEqual(
            alignment[14],
            "CCACCTCTGATATCACTGTAACTCAGCTACCCA-GAAGATAGGAAAGGTCATATCAGCCAATCATGTGCAACTTAATAGCCCTGTTCATAAAGCTGCTTTGTGTTTTATGAGGCTATTCATTCACCTTAGGGGAAATCCTGTCAAGCCATCATTTCCTTCCATATCTTCA",
        )
        self.assertEqual(alignment.sequences[14].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[14].annotations["leftCount"], 5690)
        self.assertEqual(alignment.sequences[14].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[14].annotations["rightCount"], 0)
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(record.seq), 359464)
        self.assertEqual(segment, (171035, 164755))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 15)
        self.assertEqual(len(alignment.annotations["empty"]), 1)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
            [[  3021275,   3021308,   3021308,   3021333,   3021334,   3021344,
                3021357,   3021358,   3021359,   3021369,   3021369,   3021369,
                3021369,   3021371,   3021382,   3021408,   3021415,   3021421],
             [      285,       252,       252,       227,       226,       216,
                    203,       202,       201,       191,       183,       177,
                    168,       166,       155,       129,       122,       116],
             [     3888,      3855,      3855,      3830,      3829,      3819,
                   3806,      3805,      3804,      3794,      3786,      3780,
                   3771,      3769,      3769,      3743,      3736,      3730],
             [158037166, 158037133, 158037133, 158037108, 158037107, 158037097,
              158037084, 158037083, 158037082, 158037072, 158037064, 158037058,
              158037049, 158037047, 158037036, 158037010, 158037003, 158036997],
             [157514999, 157514966, 157514966, 157514941, 157514940, 157514930,
              157514917, 157514916, 157514915, 157514905, 157514897, 157514891,
              157514882, 157514880, 157514869, 157514843, 157514836, 157514830],
             [155024945, 155024912, 155024912, 155024887, 155024886, 155024876,
              155024863, 155024862, 155024861, 155024851, 155024843, 155024837,
              155024828, 155024826, 155024815, 155024789, 155024782, 155024776],
             [    10878,     10911,     10911,     10936,     10937,     10947,
                  10960,     10961,     10962,     10972,     10980,     10986,
                  10995,     10997,     11008,     11034,     11041,     11047],
             [   323397,    323364,    323364,    323339,    323338,    323328,
                 323328,    323327,    323326,    323316,    323308,    323302,
                 323293,    323291,    323280,    323254,    323247,    323241],
             [     1615,      1582,      1581,      1556,      1555,      1545,
                   1532,      1531,      1530,      1520,      1512,      1506,
                   1497,      1495,      1484,      1458,      1451,      1445],
             [      521,       554,       554,       579,       580,       590,
                    603,       604,       605,       615,       623,       623,
                    632,       634,       645,       671,       671,       677],
             [ 47544221,  47544188,  47544188,  47544163,  47544163,  47544153,
               47544140,  47544139,  47544138,  47544128,  47544120,  47544114,
               47544105,  47544103,  47544092,  47544066,  47544059,  47544053],
             [    45798,     45765,     45765,     45740,     45739,     45729,
                  45716,     45715,     45714,     45704,     45696,     45690,
                  45681,     45679,     45668,     45642,     45635,     45629],
             [     7048,      7015,      7015,      6990,      6989,      6979,
                   6966,      6965,      6965,      6955,      6947,      6941,
                   6932,      6930,      6919,      6893,      6886,      6880],
             [    95418,     95451,     95451,     95476,     95477,     95487,
                  95500,     95501,     95502,     95512,     95520,     95526,
                  95535,     95537,     95548,     95574,     95581,     95587],
             [ 40040432,  40040399,  40040399,  40040374,  40040373,  40040363,
               40040350,  40040349,  40040348,  40040338,  40040330,  40040324,
               40040315,  40040313,  40040302,  40040276,  40040269,  40040263],
            ]
                    # fmt: on
                ),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 30254)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3021421 : 3021421 + 44],
            "ACGTTGCTCATTGTAATTGAAGCATTTATTACCAATGCCTTCCC",
        )
        self.assertEqual(
            alignment[0],
            "-ACGTTGCTCATTGT-----AATTGAAGCATTTATTACCAA--------TG--------------CCTTCCC",
        )
        self.assertEqual(alignment.sequences[1].id, "cavPor2.scaffold_216473")
        self.assertEqual(len(alignment.sequences[1].seq), 10026)
        self.assertEqual(
            alignment.sequences[1].seq[69 : 69 + 47],
            "CAGAGGGGCCGTCAGCAAAACATACTTTTATTTGTAACAAGGAACAG",
        )
        self.assertEqual(
            alignment[1],
            "-CTGTTCCTTGTTACA----AATAAAAGTATGTTTTGCTGA--------CG-------G-----CCCCTCTG",
        )
        self.assertEqual(
            alignment.sequences[1].annotations["quality"],
            "47675464566778886867867744466523545576436669356",
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[2].id, "oryCun1.scaffold_156751")
        self.assertEqual(len(alignment.sequences[2].seq), 4726)
        self.assertEqual(
            alignment.sequences[2].seq[3693 : 3693 + 37],
            "ATTAGCAGTGAATGCTTTAATCCATAATGAGGAACCA",
        )
        self.assertEqual(
            alignment[2],
            "-TGGTTCCTCATTATG----GATTAAAGCATTCACTGCTAA--------T----------------------",
        )
        self.assertEqual(
            alignment.sequences[2].annotations["quality"],
            "9999999998999998879999999999999999999",
        )
        self.assertEqual(alignment.sequences[2].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[2].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[2].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[2].annotations["rightCount"], 2345)
        self.assertEqual(alignment.sequences[3].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[3].seq), 174210431)
        self.assertEqual(
            alignment.sequences[3].seq[158036952 : 158036952 + 45],
            "GAGAAAGCATTAGCAATAAACACTTTTATTTATAATGAGCAACAG",
        )
        self.assertEqual(
            alignment[3],
            "-CTGTTGCTCATTATA----AATAAAAGTGTTTATTGCTAA--------T--------------GCTTTCTC",
        )
        self.assertEqual(alignment.sequences[3].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[3].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[3].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[3].annotations["rightCount"], 4)
        self.assertEqual(alignment.sequences[4].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[4].seq), 173908612)
        self.assertEqual(
            alignment.sequences[4].seq[157514785 : 157514785 + 45],
            "GAGAAAGCGTTAGCAATAAACACTTTTATTTATAATGAGGAACAG",
        )
        self.assertEqual(
            alignment[4],
            "-CTGTTCCTCATTATA----AATAAAAGTGTTTATTGCTAA--------C--------------GCTTTCTC",
        )
        self.assertEqual(
            alignment.sequences[4].annotations["quality"],
            "999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[4].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[4].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[4].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[4].annotations["rightCount"], 4)
        self.assertEqual(alignment.sequences[5].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[5].seq), 170899992)
        self.assertEqual(
            alignment.sequences[5].seq[155024731 : 155024731 + 45],
            "GAGAAAGCGTTAGCAATAAACACTTTTATTTATAATGAGGAACAG",
        )
        self.assertEqual(
            alignment[5],
            "-CTGTTCCTCATTATA----AATAAAAGTGTTTATTGCTAA--------C--------------GCTTTCTC",
        )
        self.assertEqual(alignment.sequences[5].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[5].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[5].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[5].annotations["rightCount"], 4)
        self.assertEqual(alignment.sequences[6].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[6].seq), 133105)
        self.assertEqual(
            alignment.sequences[6].seq[11047 : 11047 + 43],
            "CTGTTTCTCATTATAAATAAGTGTTTATTGCTAACGCTTTCTC",
        )
        self.assertEqual(
            alignment[6],
            "-CTGTTTCTCATTATA----AAT--AAGTGTTTATTGCTAA--------C--------------GCTTTCTC",
        )
        self.assertEqual(alignment.sequences[6].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[6].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[6].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[6].annotations["rightCount"], 701)
        self.assertEqual(alignment.sequences[7].id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(alignment.sequences[7].seq), 498454)
        self.assertEqual(
            alignment.sequences[7].seq[323206 : 323206 + 35],
            "TTAGCAACAAACGCTATATTTATCATGAGGAGCAG",
        )
        self.assertEqual(
            alignment[7],
            "-CTGCTCCTCATGATA----AAT-ATAGCGTTTGTTGCTAA-------------------------------",
        )
        self.assertEqual(
            alignment.sequences[7].annotations["quality"],
            "99999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[7].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[7].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[7].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[7].annotations["rightCount"], 10695)
        self.assertEqual(alignment.sequences[8].id, "sorAra1.scaffold_2476")
        self.assertEqual(len(alignment.sequences[8].seq), 4997)
        self.assertEqual(
            alignment.sequences[8].seq[1397 : 1397 + 48],
            "GGAACCTGCTGGTCATAAACGCTCtttttttAATCTAAGAAGGAACAG",
        )
        self.assertEqual(
            alignment[8],
            "-CTGTTCCTTCTTAGATTaaaaaaaGAGCGTTTATGACCAG--------CA-------GGTTCC--------",
        )
        self.assertEqual(
            alignment.sequences[8].annotations["quality"],
            "999999999999999999999999999999999766975599999999",
        )
        self.assertEqual(alignment.sequences[8].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[8].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[8].annotations["rightStatus"], "N")
        self.assertEqual(alignment.sequences[8].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[9].id, "eriEur1.scaffold_266115")
        self.assertEqual(len(alignment.sequences[9].seq), 4589)
        self.assertEqual(
            alignment.sequences[9].seq[677 : 677 + 35],
            "CTCCTTGTTCTTGAAGCCGTTTTTTATTGTCAGTG",
        )
        self.assertEqual(
            alignment[9],
            "-CTCCTTGTTCTTGAAGC----------CGTTTTTTATTGT--------CA-------GTG-----------",
        )
        self.assertEqual(
            alignment.sequences[9].annotations["quality"],
            "98959997997999999999999999999989999",
        )
        self.assertEqual(alignment.sequences[9].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[9].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[9].annotations["rightStatus"], "N")
        self.assertEqual(alignment.sequences[9].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[10].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[10].seq), 125616256)
        self.assertEqual(
            alignment.sequences[10].seq[47544013 : 47544013 + 40],
            "GACATTAGCAATAAATGCTTTTATTTATAATGAGAAACAG",
        )
        self.assertEqual(
            alignment[10],
            "-CTGTTTCTCATTATA----AATAAAAGCATTTATTGCTAA--------TG-------T------------C",
        )
        self.assertEqual(
            alignment.sequences[10].annotations["quality"],
            "9999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[10].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[10].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[10].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[10].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[11].id, "dasNov1.scaffold_56749")
        self.assertEqual(len(alignment.sequences[11].seq), 10470)
        self.assertEqual(
            alignment.sequences[11].seq[6836 : 6836 + 44],
            "GAGGAGCCGTTAGCAGTAAATGCTTTTATTATAAAGAGGAACAG",
        )
        self.assertEqual(
            alignment[11],
            "-CTGTTCCTCTTTAT-----AATAAAAGCATTTACTGCTAA--------CGGCTCCTC--------------",
        )
        self.assertEqual(
            alignment.sequences[11].annotations["quality"],
            "99999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[11].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[11].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[11].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[11].annotations["rightCount"], 904)
        self.assertEqual(alignment.sequences[12].id, "echTel1.scaffold_288249")
        self.assertEqual(len(alignment.sequences[12].seq), 100002)
        self.assertEqual(
            alignment.sequences[12].seq[95587 : 95587 + 38],
            "CTGATCCTCATTATAAATAAAAGTGTTTGTTACTAATG",
        )
        self.assertEqual(
            alignment[12],
            "-CTGATCCTCATTATA----AATAAAAGTGTTTGTTACTAA--------TG---------------------",
        )
        self.assertEqual(
            alignment.sequences[12].annotations["quality"],
            "99999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[12].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[12].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[12].annotations["rightStatus"], "N")
        self.assertEqual(alignment.sequences[12].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[13].id, "ornAna1.chr2")
        self.assertEqual(len(alignment.sequences[13].seq), 54797317)
        self.assertEqual(
            alignment.sequences[13].seq[40040218 : 40040218 + 45],
            "GGGGGATCTTTAGTAATAAAAGAGTCCCTGTACAATGAGGACAGT",
        )
        self.assertEqual(
            alignment[13],
            "ACTG-TCCTCATTGTA----CAGGGACTCTTTTATTACTAAAGATCCCCC----------------------",
        )
        self.assertEqual(alignment.sequences[13].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[13].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[13].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[13].annotations["rightCount"], 0)
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "felCat3.scaffold_205680")
        self.assertEqual(len(record.seq), 119354)
        self.assertEqual(segment, (45629, 11464))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(record.seq), 359464)
        self.assertEqual(segment, (171035, 164755))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 14)
        self.assertEqual(len(alignment.annotations["empty"]), 2)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
            [[  3021421,   3021421,   3021424,   3021425,   3021435,   3021435,
                3021435,   3021435,   3021438,   3021439,   3021440,   3021443,
                3021456,   3021456,   3021457,   3021458,   3021458,   3021458,
                3021458,   3021458,   3021458,   3021464,   3021465],
             [      116,       116,       113,       112,       102,       101,
                    101,       101,        98,        97,        96,        93,
                     80,        80,        79,        78,        78,        77,
                     77,        77,        76,        70,        69],
             [     3730,      3730,      3727,      3726,      3716,      3715,
                   3715,      3715,      3712,      3711,      3710,      3707,
                   3694,      3694,      3693,      3693,      3693,      3693,
                   3693,      3693,      3693,      3693,      3693],
             [158036997, 158036997, 158036994, 158036993, 158036983, 158036982,
              158036982, 158036982, 158036979, 158036978, 158036977, 158036974,
              158036961, 158036961, 158036960, 158036960, 158036960, 158036960,
              158036960, 158036960, 158036959, 158036953, 158036952],
             [157514830, 157514830, 157514827, 157514826, 157514816, 157514815,
              157514815, 157514815, 157514812, 157514811, 157514810, 157514807,
              157514794, 157514794, 157514793, 157514793, 157514793, 157514793,
              157514793, 157514793, 157514792, 157514786, 157514785],
             [155024776, 155024776, 155024773, 155024772, 155024762, 155024761,
              155024761, 155024761, 155024758, 155024757, 155024756, 155024753,
              155024740, 155024740, 155024739, 155024739, 155024739, 155024739,
              155024739, 155024739, 155024738, 155024732, 155024731],
             [    11047,     11047,     11050,     11051,     11061,     11062,
                  11062,     11062,     11065,     11065,     11065,     11068,
                  11081,     11081,     11082,     11082,     11082,     11082,
                  11082,     11082,     11083,     11089,     11090],
             [   323241,    323241,    323238,    323237,    323227,    323226,
                 323226,    323226,    323223,    323223,    323222,    323219,
                 323206,    323206,    323206,    323206,    323206,    323206,
                 323206,    323206,    323206,    323206,    323206],
             [     1445,      1445,      1442,      1441,      1431,      1430,
                   1428,      1426,      1423,      1422,      1421,      1418,
                   1405,      1405,      1404,      1403,      1403,      1402,
                   1400,      1397,      1397,      1397,      1397],
             [      677,       677,       680,       681,       691,       692,
                    694,       694,       694,       694,       694,       694,
                    707,       707,       708,       709,       709,       710,
                    712,       712,       712,       712,       712],
             [ 47544053,  47544053,  47544050,  47544049,  47544039,  47544038,
               47544038,  47544038,  47544035,  47544034,  47544033,  47544030,
               47544017,  47544017,  47544016,  47544015,  47544015,  47544014,
               47544014,  47544014,  47544014,  47544014,  47544013],
             [     6880,      6880,      6877,      6876,      6866,      6866,
                   6866,      6866,      6863,      6862,      6861,      6858,
                   6845,      6845,      6844,      6843,      6836,      6836,
                   6836,      6836,      6836,      6836,      6836],
             [    95587,     95587,     95590,     95591,     95601,     95602,
                  95602,     95602,     95605,     95606,     95607,     95610,
                  95623,     95623,     95624,     95625,     95625,     95625,
                  95625,     95625,     95625,     95625,     95625],
             [ 40040263,  40040262,  40040259,  40040259,  40040249,  40040248,
               40040248,  40040248,  40040245,  40040244,  40040243,  40040240,
               40040227,  40040219,  40040218,  40040218,  40040218,  40040218,
               40040218,  40040218,  40040218,  40040218,  40040218],
            ]
                    # fmt: on
                ),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, -9167)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3021465 : 3021465 + 29],
            "CCCTACACTGTCAAGTGGGAGGAGACAGT",
        )
        self.assertEqual(
            alignment[0], "CCCT--ACACTGTC----AAGTGGGAGGAGACAGT--------------------"
        )
        self.assertEqual(alignment.sequences[1].id, "cavPor2.scaffold_216473")
        self.assertEqual(len(alignment.sequences[1].seq), 10026)
        self.assertEqual(
            alignment.sequences[1].seq[41 : 41 + 28], "accatcccaccccacccccagtgtGGCT"
        )
        self.assertEqual(
            alignment[1], "AGCC--acactgg-----gggtggggtgggatggt--------------------"
        )
        self.assertEqual(
            alignment.sequences[1].annotations["quality"],
            "7667687856666544895554554677",
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "N")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 0)
        self.assertEqual(alignment.sequences[2].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[2].seq), 174210431)
        self.assertEqual(
            alignment.sequences[2].seq[158036924 : 158036924 + 24],
            "GCTGCTCCTATCGCCCCCACAGGG",
        )
        self.assertEqual(
            alignment[2], "-------CCCTGTG----GGGGCGATAGGAGCAGC--------------------"
        )
        self.assertEqual(alignment.sequences[2].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[2].annotations["leftCount"], 4)
        self.assertEqual(alignment.sequences[2].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[2].annotations["rightCount"], 9)
        self.assertEqual(alignment.sequences[3].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[3].seq), 173908612)
        self.assertEqual(
            alignment.sequences[3].seq[157514757 : 157514757 + 24],
            "GCTGCCCCTATCGCCCCCACATGG",
        )
        self.assertEqual(
            alignment[3], "-------CCATGTG----GGGGCGATAGGGGCAGC--------------------"
        )
        self.assertEqual(
            alignment.sequences[3].annotations["quality"], "999999999999999999999999"
        )
        self.assertEqual(alignment.sequences[3].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[3].annotations["leftCount"], 4)
        self.assertEqual(alignment.sequences[3].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[3].annotations["rightCount"], 9)
        self.assertEqual(alignment.sequences[4].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[4].seq), 170899992)
        self.assertEqual(
            alignment.sequences[4].seq[155024703 : 155024703 + 24],
            "GCTGCCCCTATCGCCCCCACATGG",
        )
        self.assertEqual(
            alignment[4], "-------CCATGTG----GGGGCGATAGGGGCAGC--------------------"
        )
        self.assertEqual(alignment.sequences[4].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[4].annotations["leftCount"], 4)
        self.assertEqual(alignment.sequences[4].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[4].annotations["rightCount"], 9)
        self.assertEqual(alignment.sequences[5].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[5].seq), 125616256)
        self.assertEqual(
            alignment.sequences[5].seq[47543982 : 47543982 + 31],
            "CACACACACACACCCCCCATGTGGTTCAGGG",
        )
        self.assertEqual(
            alignment[5], "CCCTGAACCACATG----GGGGGTGTGTGTGTGTG--------------------"
        )
        self.assertEqual(
            alignment.sequences[5].annotations["quality"],
            "9999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[5].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[5].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[5].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[5].annotations["rightCount"], 13)
        self.assertEqual(alignment.sequences[6].id, "ornAna1.chr2")
        self.assertEqual(len(alignment.sequences[6].seq), 54797317)
        self.assertEqual(
            alignment.sequences[6].seq[40040173 : 40040173 + 45],
            "GCCAGCAAAGATGAAAGAGGGCTACATCCAAACTCCTATGACACA",
        )
        self.assertEqual(
            alignment[6], "----------TGTGTCATAGGAGTTTGGATGTAGCCCTCTTTCATCTTTGCTGGC"
        )
        self.assertEqual(alignment.sequences[6].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[6].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[6].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[6].annotations["rightCount"], 0)
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "dasNov1.scaffold_56749")
        self.assertEqual(len(record.seq), 10470)
        self.assertEqual(segment, (6836, 5932))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "felCat3.scaffold_205680")
        self.assertEqual(len(record.seq), 119354)
        self.assertEqual(segment, (45629, 11464))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][2]
        (record, segment, status) = empty
        self.assertEqual(record.id, "calJac1.Contig6394")
        self.assertEqual(len(record.seq), 133105)
        self.assertEqual(segment, (11090, 11791))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][3]
        (record, segment, status) = empty
        self.assertEqual(record.id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(record.seq), 498454)
        self.assertEqual(segment, (323206, 312511))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][4]
        (record, segment, status) = empty
        self.assertEqual(record.id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(record.seq), 359464)
        self.assertEqual(segment, (171035, 164755))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][5]
        (record, segment, status) = empty
        self.assertEqual(record.id, "oryCun1.scaffold_156751")
        self.assertEqual(len(record.seq), 4726)
        self.assertEqual(segment, (3693, 1348))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 7)
        self.assertEqual(len(alignment.annotations["empty"]), 6)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
            [[  3021465,   3021469,   3021469,   3021470,   3021473,   3021476,
                3021477,   3021477,   3021494,   3021494],
             [       69,        65,        65,        64,        61,        58,
                     58,        58,        41,        41],
             [158036948, 158036948, 158036948, 158036948, 158036945, 158036942,
              158036941, 158036941, 158036924, 158036924],
             [157514781, 157514781, 157514781, 157514781, 157514778, 157514775,
              157514774, 157514774, 157514757, 157514757],
             [155024727, 155024727, 155024727, 155024727, 155024724, 155024721,
              155024720, 155024720, 155024703, 155024703],
             [ 47544013,  47544009,  47544007,  47544006,  47544003,  47544000,
               47543999,  47543999,  47543982,  47543982],
             [ 40040218,  40040218,  40040218,  40040218,  40040218,  40040215,
               40040214,  40040210,  40040193,  40040173],
            ]
                    # fmt: on
                ),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 15763)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3021494 : 3021494 + 42],
            "TGTTTAGTACCATGCTTAGGAATGATAAACTCACTTAGTGtt",
        )
        self.assertEqual(alignment[0], "TGTTTAGTACC----ATGCTTAGGAATGATAAACTCACTTAGTGtt")
        self.assertEqual(alignment.sequences[1].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[1].seq), 174210431)
        self.assertEqual(
            alignment.sequences[1].seq[158036869 : 158036869 + 46],
            "AAGATTGGGTGAGCCTATCACGCCAAAGAATAAAGGACATGCAACA",
        )
        self.assertEqual(alignment[1], "TGTTGCATGTCCTTTATTCTTTGGCGTGATAGGCTCACCCAATCTT")
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 9)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 943)
        self.assertEqual(alignment.sequences[2].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[2].seq), 173908612)
        self.assertEqual(
            alignment.sequences[2].seq[157514702 : 157514702 + 46],
            "AAGATTGGGTGAGCCTATCACGCCAAAGAATAAAGGATATGCAACA",
        )
        self.assertEqual(alignment[2], "TGTTGCATATCCTTTATTCTTTGGCGTGATAGGCTCACCCAATCTT")
        self.assertEqual(
            alignment.sequences[2].annotations["quality"],
            "9999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[2].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[2].annotations["leftCount"], 9)
        self.assertEqual(alignment.sequences[2].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[2].annotations["rightCount"], 10)
        self.assertEqual(alignment.sequences[3].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[3].seq), 170899992)
        self.assertEqual(
            alignment.sequences[3].seq[155024648 : 155024648 + 46],
            "AAGATTGGGTGAGCCTATCACGCCAAAGAATAAACGACATGCAACA",
        )
        self.assertEqual(alignment[3], "TGTTGCATGTCGTTTATTCTTTGGCGTGATAGGCTCACCCAATCTT")
        self.assertEqual(alignment.sequences[3].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[3].annotations["leftCount"], 9)
        self.assertEqual(alignment.sequences[3].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[3].annotations["rightCount"], 931)
        self.assertEqual(alignment.sequences[4].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[4].seq), 125616256)
        self.assertEqual(
            alignment.sequences[4].seq[47543923 : 47543923 + 46],
            "ATGATGGAGTGAAGCTATCACTTTGAACAGCAAGTGAGACTTAACA",
        )
        self.assertEqual(alignment[4], "TGTTAAGTCTCACTTGCTGTTCAAAGTGATAGCTTCACTCCATCAT")
        self.assertEqual(
            alignment.sequences[4].annotations["quality"],
            "9999999999999999999999999999999999999999999999",
        )
        self.assertEqual(alignment.sequences[4].annotations["leftStatus"], "I")
        self.assertEqual(alignment.sequences[4].annotations["leftCount"], 13)
        self.assertEqual(alignment.sequences[4].annotations["rightStatus"], "I")
        self.assertEqual(alignment.sequences[4].annotations["rightCount"], 1)
        self.assertEqual(alignment.sequences[5].id, "ornAna1.chr2")
        self.assertEqual(len(alignment.sequences[5].seq), 54797317)
        self.assertEqual(
            alignment.sequences[5].seq[40040137 : 40040137 + 36],
            "TCCAGTGAGTAGAAGTTCTAGCAATCATTTTAAACA",
        )
        self.assertEqual(alignment[5], "TGTTTAAAATG----ATTGCTAGAACTTCTA--CTCACTGGA----")
        self.assertEqual(alignment.sequences[5].annotations["leftStatus"], "C")
        self.assertEqual(alignment.sequences[5].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[5].annotations["rightStatus"], "C")
        self.assertEqual(alignment.sequences[5].annotations["rightCount"], 0)
        empty = alignment.annotations["empty"][0]
        (record, segment, status) = empty
        self.assertEqual(record.id, "dasNov1.scaffold_56749")
        self.assertEqual(len(record.seq), 10470)
        self.assertEqual(segment, (6836, 5932))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][1]
        (record, segment, status) = empty
        self.assertEqual(record.id, "felCat3.scaffold_205680")
        self.assertEqual(len(record.seq), 119354)
        self.assertEqual(segment, (45629, 11464))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][2]
        (record, segment, status) = empty
        self.assertEqual(record.id, "calJac1.Contig6394")
        self.assertEqual(len(record.seq), 133105)
        self.assertEqual(segment, (11090, 11791))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][3]
        (record, segment, status) = empty
        self.assertEqual(record.id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(record.seq), 498454)
        self.assertEqual(segment, (323206, 312511))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][4]
        (record, segment, status) = empty
        self.assertEqual(record.id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(record.seq), 359464)
        self.assertEqual(segment, (171035, 164755))
        self.assertEqual(status, "I")
        empty = alignment.annotations["empty"][5]
        (record, segment, status) = empty
        self.assertEqual(record.id, "oryCun1.scaffold_156751")
        self.assertEqual(len(record.seq), 4726)
        self.assertEqual(segment, (3693, 1348))
        self.assertEqual(status, "I")
        self.assertEqual(len(alignment.sequences), 6)
        self.assertEqual(len(alignment.annotations["empty"]), 6)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    # fmt: off
            [[  3021494,   3021505,   3021505,   3021521,   3021523,   3021532,
                3021536],
             [158036915, 158036904, 158036900, 158036884, 158036882, 158036873,
              158036869],
             [157514748, 157514737, 157514733, 157514717, 157514715, 157514706,
              157514702],
             [155024694, 155024683, 155024679, 155024663, 155024661, 155024652,
              155024648],
             [ 47543969,  47543958,  47543954,  47543938,  47543936,  47543927,
               47543923],
             [ 40040173,  40040162,  40040162,  40040146,  40040146,  40040137,
               40040137],
            ]
                    # fmt: on
                ),
            )
        )

    def test_reading_missing_signature(self):
        """Test parsing MAF file ucsc_mm9_chr10_big.maf with missing signature."""
        path = "MAF/ucsc_mm9_chr10_big.maf"
        with self.assertRaises(ValueError) as cm:
            maf.AlignmentIterator(path)
        self.assertEqual(str(cm.exception), "header line does not start with ##maf")

    def test_reading_ucsc_mm9_chr10_bad(self):
        """Test parsing MAF file ucsc_mm9_chr10_bad.maf with incorrect sequence size."""
        path = "MAF/ucsc_mm9_chr10_bad.maf"
        alignments = maf.AlignmentIterator(path)
        self.assertEqual(alignments.metadata["version"], "1")
        self.assertEqual(alignments.metadata["scoring"], "autoMZ.v1")
        next(alignments)
        next(alignments)
        next(alignments)
        next(alignments)
        next(alignments)
        next(alignments)
        with self.assertRaises(ValueError) as cm:
            next(alignments)
        self.assertEqual(
            str(cm.exception), "sequence size is incorrect (found 219, expected 319)"
        )

    def test_reading_length_coords_mismatch(self):
        """Test parsing inconsistent MAF file length_coords_mismatch.maf."""
        path = "MAF/length_coords_mismatch.maf"
        alignments = maf.AlignmentIterator(path)
        self.assertEqual(alignments.metadata["version"], "1")
        self.assertEqual(alignments.metadata["scoring"], "autoMZ.v1")
        alignment = next(alignments)
        self.assertEqual(alignment.score, 6441)
        self.assertEqual(len(alignment.sequences), 2)
        self.assertEqual(alignment.sequences[0].id, "mm8.chr10")
        self.assertEqual(len(alignment.sequences[0]), 129993255)
        self.assertEqual(
            alignment.sequences[0].seq[3009319 : 3009319 + 162],
            "TCATAGGTATTTATTTTTAAATATGGTTTGCTTTATGGCTAGAACACACCGATTACTTAAAATAGGATTAACCCCCATACACTTTAAAAATGATTAAACAACATTTCTGCTGCTCGCTCACATTCTTCATAGAAGATGACATAATGTATTTTCCTTTTGGTT",
        )
        self.assertEqual(
            alignment[0],
            "TCATAGGTATTTATTTTTAAATATGGTTTGCTTTATGGCTAGAACACACCGATTACTTAAAATAGGATTAACC--CCCATACACTTTAAAAATGATTAAACAACATTTCTGCTGCTCGCTCACATTCTTCATAGAAGATGACATAATGTATTTTCCTTTTGGTT",
        )
        self.assertEqual(alignment.sequences[1].id, "oryCun1.scaffold_133159")
        self.assertEqual(len(alignment.sequences[1]), 13221)
        self.assertEqual(
            alignment.sequences[1].seq[11087 : 11087 + 164],
            "TCACAGATATTTACTATTAAATATGGTTTGTTATATGGTTACGGTTCATAGGTTACTTGGAATTGGATTAACCTTCTTATTCATTGCAGAATTGGTTACACTGTGTTCTTGACCTTTGCTTGTTTTCTCCATGGAAACTGATGTCAAATACTTTCCCTTTGGTT",
        )
        self.assertEqual(
            alignment[1],
            "TCACAGATATTTACTATTAAATATGGTTTGTTATATGGTTACGGTTCATAGGTTACTTGGAATTGGATTAACCTTCTTATTCATTGCAGAATTGGTTACACTGTGTTCTTGACCTTTGCTTGTTTTCTCCATGGAAACTGATGTCAAATACTTTCCCTTTGGTT",
        )
        self.assertEqual(
            alignment.sequences[1].annotations["quality"],
            "99569899999998999999999999999999999999999999999999999999999999999999999757878999975999999999999999979999999999997899999999999997997999999869999996999988997997999999",
        )
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "N")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "N")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 0)
        with self.assertRaises(ValueError) as cm:
            next(alignments)
        self.assertEqual(
            str(cm.exception), "sequence size is incorrect (found 219, expected 319)"
        )

    def test_reading_bug2453(self):
        """Test parsing bug2453.maf."""
        path = "MAF/bug2453.maf"
        alignments = maf.AlignmentIterator(path)
        self.assertEqual(len(alignments.metadata), 3)
        self.check_ucsc_test(alignments)

    def test_reading_ucsc_test(self):
        """Test parsing ucsc_test.maf."""
        path = "MAF/ucsc_test.maf"
        alignments = maf.AlignmentIterator(path)
        self.assertEqual(len(alignments.metadata), 9)
        self.assertEqual(alignments.metadata["name"], "euArc")
        self.assertEqual(alignments.metadata["visibility"], "pack")
        self.assertEqual(alignments.metadata["mafDot"], "off")
        self.assertEqual(alignments.metadata["frames"], "multiz28wayFrames")
        self.assertEqual(
            alignments.metadata["speciesOrder"],
            ["hg16", "panTro1", "baboon", "mm4", "rn3"],
        )
        self.assertEqual(alignments.metadata["description"], "A sample alignment")
        self.check_ucsc_test(alignments)

    def check_ucsc_test(self, alignments):
        self.assertEqual(alignments.metadata["version"], "1")
        self.assertEqual(alignments.metadata["scoring"], "tba.v8")
        self.assertEqual(
            alignments.metadata["comments"],
            [
                "tba.v8 (((human chimp) baboon) (mouse rat))",
                "multiz.v7",
                "maf_project.v5 _tba_right.maf3 mouse _tba_C",
                "single_cov2.v4 single_cov2 /dev/stdin",
            ],
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 23262)
        self.assertEqual(len(alignment.sequences), 5)
        self.assertEqual(alignment.sequences[0].id, "hg16.chr7")
        self.assertEqual(len(alignment.sequences[0]), 158545518)
        self.assertEqual(
            alignment.sequences[0].seq[27578828 : 27578828 + 38],
            "AAAGGGAATGTTAACCAAATGAATTGTCTCTTACGGTG",
        )
        self.assertEqual(alignment[0], "AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG")
        self.assertEqual(alignment.sequences[1].id, "panTro1.chr6")
        self.assertEqual(len(alignment.sequences[1]), 161576975)
        self.assertEqual(
            alignment.sequences[1].seq[28741140 : 28741140 + 38],
            "AAAGGGAATGTTAACCAAATGAATTGTCTCTTACGGTG",
        )
        self.assertEqual(alignment[1], "AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG")
        self.assertEqual(alignment.sequences[2].id, "baboon")
        self.assertEqual(len(alignment.sequences[2]), 4622798)
        self.assertEqual(
            alignment.sequences[2].seq[116834 : 116834 + 38],
            "AAAGGGAATGTTAACCAAATGAGTTGTCTCTTATGGTG",
        )
        self.assertEqual(alignment[2], "AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG")
        self.assertEqual(alignment.sequences[3].id, "mm4.chr6")
        self.assertEqual(len(alignment.sequences[3]), 151104725)
        self.assertEqual(
            alignment.sequences[3].seq[53215344 : 53215344 + 38],
            "AATGGGAATGTTAAGCAAACGAATTGTCTCTCAGTGTG",
        )
        self.assertEqual(alignment[3], "-AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG")
        self.assertEqual(alignment.sequences[4].id, "rn3.chr4")
        self.assertEqual(len(alignment.sequences[4]), 187371129)
        self.assertEqual(
            alignment.sequences[4].seq[81344243 : 81344243 + 40],
            "AAGGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG",
        )
        self.assertEqual(alignment[4], "-AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [
                        # fmt: off
        [27578828, 27578829, 27578831, 27578831, 27578850, 27578850, 27578866],
        [28741140, 28741141, 28741143, 28741143, 28741162, 28741162, 28741178],
        [  116834,   116835,   116837,   116837,   116856,   116856,   116872],
        [53215344, 53215344, 53215346, 53215347, 53215366, 53215366, 53215382],
        [81344243, 81344243, 81344245, 81344245, 81344264, 81344267, 81344283],
                        # fmt: on
                    ]
                ),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 5062.0)
        self.assertEqual(len(alignment.sequences), 5)
        self.assertEqual(alignment.sequences[0].id, "hg16.chr7")
        self.assertEqual(len(alignment.sequences[0]), 158545518)
        self.assertEqual(alignment.sequences[0].seq[27699739 : 27699739 + 6], "TAAAGA")
        self.assertEqual(alignment[0], "TAAAGA")
        self.assertEqual(alignment.sequences[1].id, "panTro1.chr6")
        self.assertEqual(len(alignment.sequences[1]), 161576975)
        self.assertEqual(alignment.sequences[1].seq[28862317 : 28862317 + 6], "TAAAGA")
        self.assertEqual(alignment[1], "TAAAGA")
        self.assertEqual(alignment.sequences[2].id, "baboon")
        self.assertEqual(len(alignment.sequences[2]), 4622798)
        self.assertEqual(alignment.sequences[2].seq[241163 : 241163 + 6], "TAAAGA")
        self.assertEqual(alignment[2], "TAAAGA")
        self.assertEqual(alignment.sequences[3].id, "mm4.chr6")
        self.assertEqual(len(alignment.sequences[3]), 151104725)
        self.assertEqual(alignment.sequences[3].seq[53303881 : 53303881 + 6], "TAAAGA")
        self.assertEqual(alignment[3], "TAAAGA")
        self.assertEqual(alignment.sequences[4].id, "rn3.chr4")
        self.assertEqual(len(alignment.sequences[4]), 187371129)
        self.assertEqual(alignment.sequences[4].seq[81444246 : 81444246 + 6], "taagga")
        self.assertEqual(alignment[4], "taagga")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [
                        # fmt: off
                             [27699739, 27699745],
                             [28862317, 28862323],
                             [  241163,   241169],
                             [53303881, 53303887],
                             [81444246, 81444252],
                        # fmt: on
                    ]
                ),
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 6636.0)
        self.assertEqual(len(alignment.sequences), 4)
        self.assertEqual(alignment.sequences[0].id, "hg16.chr7")
        self.assertEqual(len(alignment.sequences[0]), 158545518)
        self.assertEqual(
            alignment.sequences[0].seq[27707221 : 27707221 + 13], "gcagctgaaaaca"
        )
        self.assertEqual(alignment[0], "gcagctgaaaaca")
        self.assertEqual(alignment.sequences[1].id, "panTro1.chr6")
        self.assertEqual(len(alignment.sequences[1]), 161576975)
        self.assertEqual(
            alignment.sequences[1].seq[28869787 : 28869787 + 13], "gcagctgaaaaca"
        )
        self.assertEqual(alignment[1], "gcagctgaaaaca")
        self.assertEqual(alignment.sequences[2].id, "baboon")
        self.assertEqual(len(alignment.sequences[2]), 4622798)
        self.assertEqual(
            alignment.sequences[2].seq[249182 : 249182 + 13], "gcagctgaaaaca"
        )
        self.assertEqual(alignment[2], "gcagctgaaaaca")
        self.assertEqual(alignment.sequences[3].id, "mm4.chr6")
        self.assertEqual(len(alignment.sequences[3]), 151104725)
        self.assertEqual(
            alignment.sequences[3].seq[53310102 : 53310102 + 13], "ACAGCTGAAAATA"
        )
        self.assertEqual(alignment[3], "ACAGCTGAAAATA")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [
                        # fmt: off
                    [27707221, 27707234],
                    [28869787, 28869800],
                    [  249182,   249195],
                    [53310102, 53310115],
                        # fmt: on
                    ]
                ),
            )
        )
        self.assertRaises(StopIteration, next, alignments)


class TestAlign_writing(unittest.TestCase):
    def test_writing_ucsc_test(self):
        """Test reading and writing ucsc_test.maf."""
        path = "MAF/ucsc_test.maf"
        alignments = maf.AlignmentIterator(path)
        output = StringIO()
        writer = maf.AlignmentWriter(output)
        writer.write_file(alignments)
        output.seek(0)
        with open(path) as stream:
            line1 = next(output)
            line2 = next(stream)
            self.assertEqual(
                line1,
                'track name=euArc visibility=pack mafDot=off frames=multiz28wayFrames speciesOrder="hg16 panTro1 baboon mm4 rn3" description="A sample alignment"\n',
            )
            self.assertEqual(
                line2,
                'track name=euArc visibility=pack mafDot=off frames="multiz28wayFrames" speciesOrder="hg16 panTro1 baboon mm4 rn3" description="A sample alignment"\n',
            )
            for line1, line2 in zip(output, stream):
                if line1.startswith("a score="):
                    prefix1, score1 = line1.split("=")
                    prefix2, score2 = line2.split("=")
                    self.assertEqual(prefix1, prefix2)
                    self.assertAlmostEqual(float(score1), float(score2))
                else:
                    self.assertEqual(line1.strip(), line2.strip())

    def test_writing_bug2453(self):
        """Test reading and writing bug2453.maf."""
        path = "MAF/bug2453.maf"
        alignments = maf.AlignmentIterator(path)
        output = StringIO()
        writer = maf.AlignmentWriter(output)
        writer.write_file(alignments)
        output.seek(0)
        with open(path) as stream:
            for line1, line2 in zip(output, stream):
                if line1.startswith("a score="):
                    prefix1, score1 = line1.split("=")
                    prefix2, score2 = line2.split("=")
                    self.assertEqual(prefix1, prefix2)
                    self.assertAlmostEqual(float(score1), float(score2))
                else:
                    self.assertEqual(line1.strip(), line2.strip())

    def test_writing_bundle_without_target(self):
        """Test reading and writing bundle_without_target.maf."""
        path = "MAF/bundle_without_target.maf"
        alignments = maf.AlignmentIterator(path)
        output = StringIO()
        writer = maf.AlignmentWriter(output)
        writer.write_file(alignments)
        output.seek(0)
        with open(path) as stream:
            for line1, line2 in zip(output, stream):
                self.assertEqual(line1, line2)

    def test_writing_ucsc_mm9_chr10(self):
        """Test reading and writing ucsc_mm9_chr10.maf."""
        path = "MAF/ucsc_mm9_chr10.maf"
        alignments = maf.AlignmentIterator(path)
        output = StringIO()
        writer = maf.AlignmentWriter(output)
        writer.write_file(alignments)
        output.seek(0)
        with open(path) as stream:
            for line1, line2 in zip(output, stream):
                words1 = line1.split()
                words2 = line2.split()
                if line1.startswith("e "):
                    (
                        prefix1,
                        name1,
                        start1,
                        size1,
                        strand1,
                        length1,
                        status1,
                    ) = line1.split()
                    (
                        prefix2,
                        name2,
                        start2,
                        size2,
                        strand2,
                        length2,
                        status2,
                    ) = line2.split()
                    size1 = int(size1)
                    size2 = int(size2)
                    if size1 == 0 and strand1 != strand2:
                        start1 = int(start1)
                        start2 = int(start2)
                        size1 = int(size1)
                        size2 = int(size2)
                        length1 = int(length1)
                        length2 = int(length2)
                        self.assertEqual(prefix1, "e")
                        self.assertEqual(prefix2, "e")
                        self.assertEqual(name1, name2)
                        self.assertEqual(status1, status2)
                        self.assertEqual(size1, size2)
                        self.assertEqual(length1, length2)
                        self.assertEqual(start1, length2 - start2 - size2)
                        self.assertEqual(start2, length1 - start1 - size1)
                        continue
                self.assertEqual(words1, words2)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
