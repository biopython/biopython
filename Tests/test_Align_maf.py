# Copyright 2021 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Align.maf module."""
import unittest


from Bio.Align import maf


try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align.maf."
    ) from None


class TestAlignIO_reading(unittest.TestCase):
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
        self.assertEqual(alignment.sequences[0].seq[3009319: 3009319+162], "TCATAGGTATTTATTTTTAAATATGGTTTGCTTTATGGCTAGAACACACCGATTACTTAAAATAGGATTAACCCCCATACACTTTAAAAATGATTAAACAACATTTCTGCTGCTCGCTCACATTCTTCATAGAAGATGACATAATGTATTTTCCTTTTGGTT")
        self.assertEqual(alignment.sequences[1].seq[11087: 11087+164], "TCACAGATATTTACTATTAAATATGGTTTGTTATATGGTTACGGTTCATAGGTTACTTGGAATTGGATTAACCTTCTTATTCATTGCAGAATTGGTTACACTGTGTTCTTGACCTTTGCTTGTTTTCTCCATGGAAACTGATGTCAAATACTTTCCCTTTGGTT")
        self.assertEqual(alignment[0], "TCATAGGTATTTATTTTTAAATATGGTTTGCTTTATGGCTAGAACACACCGATTACTTAAAATAGGATTAACC--CCCATACACTTTAAAAATGATTAAACAACATTTCTGCTGCTCGCTCACATTCTTCATAGAAGATGACATAATGTATTTTCCTTTTGGTT")
        self.assertEqual(alignment[1], "TCACAGATATTTACTATTAAATATGGTTTGTTATATGGTTACGGTTCATAGGTTACTTGGAATTGGATTAACCTTCTTATTCATTGCAGAATTGGTTACACTGTGTTCTTGACCTTTGCTTGTTTTCTCCATGGAAACTGATGTCAAATACTTTCCCTTTGGTT")
        self.assertEqual(alignment.column_annotations["oryCun1.scaffold_133159"], "99569899999998999999999999999999999999999999999999999999999999999999999757878999975999999999999999979999999999997899999999999997997999999869999996999988997997999999")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'N')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'N')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 0)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[3009319, 3009392, 3009392, 3009481],
                             [  11087,   11160,   11162,   11251],
                            ])
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
        self.assertEqual(alignment.sequences[0].seq[3009319:3009319+162], "TCATAGGTATTTATTTTTAAATATGGTTTGCTTTATGGCTAGAACACACCGATTACTTAAAATAGGATTAACCCCCATACACTTTAAAAATGATTAAACAACATTTCTGCTGCTCGCTCACATTCTTCATAGAAGATGACATAATGTATTTTCCTTTTGGTT")
        self.assertEqual(alignment[0], "TCATAGGTATTTATTTTTAAATATGGTTTGCTTTATGGCTAGAACACACCGATTACTTAAAATAGGATTAACC--CCCATACACTTTAAAAATGATTAAACAACATTTCTGCTGCTCGCTCACATTCTTCATAGAAGATGACATAATGTATTTTCCTTTTGGTT")
        self.assertEqual(alignment.sequences[1].id, "oryCun1.scaffold_133159")
        self.assertEqual(len(alignment.sequences[1].seq), 13221)
        self.assertEqual(alignment.sequences[1].seq[11087:11087+164], "TCACAGATATTTACTATTAAATATGGTTTGTTATATGGTTACGGTTCATAGGTTACTTGGAATTGGATTAACCTTCTTATTCATTGCAGAATTGGTTACACTGTGTTCTTGACCTTTGCTTGTTTTCTCCATGGAAACTGATGTCAAATACTTTCCCTTTGGTT")
        self.assertEqual(alignment[1], "TCACAGATATTTACTATTAAATATGGTTTGTTATATGGTTACGGTTCATAGGTTACTTGGAATTGGATTAACCTTCTTATTCATTGCAGAATTGGTTACACTGTGTTCTTGACCTTTGCTTGTTTTCTCCATGGAAACTGATGTCAAATACTTTCCCTTTGGTT")
        self.assertEqual(alignment.column_annotations["oryCun1.scaffold_133159"], "99569899999998999999999999999999999999999999999999999999999999999999999757878999975999999999999999979999999999997899999999999997997999999869999996999988997997999999")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'N')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'N')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 0)
        self.assertEqual(len(alignment.sequences), 2)
        self.assertNotIn("empty", alignment.annotations)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[3009319, 3009392, 3009392, 3009481],
                             [  11087,   11160,   11162,   11251]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 103072)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3012076:3012076+365], "AGTCTTTCCAATGGGACCTGTGAGTCCTAACTATGCCAGCACTCCCAACAGCAAGACACTAAGTTCACTCATCCTTGGTGGATGGGATTTTGCTCCTGGAGTGTCACCAAATTAAATAACCAGTGAGCAGAGTTGTGACGAGCATCAGGCTCTGGATTTAGGTGAGAGACCTTAGTGTATGTCTCCTGTAGGTCGCAGCTCCCTATGGATGAGTCAAGTGAAGGTCCTGAGACAACAAGTCCTCGGCTATGTGGGGGTGAGGGATGCAGCTGGAACCTCAGGGATCTCTGTAAGCAGTGGCATAAATGCTTGGCGGGAGAGAGCATGTTAGAGCTCACACGACATAGGAAGCCACTGAGACACTG")
        self.assertEqual(alignment[0], "AGTCTTTCCAATGGGACCTGTGAGTCCTAACTATGCCAGC-----ACTCCCAACAGCAAGACACTAAGTT---------CACTCATCCTTGGTGGATGGGATTTTGCTCCTGGAGTGTCAC-----CAAATTAAATAACCAGTGAGCAGAGTTG--TGACGAGCATCAGGCTCTGGATTTAGGTGAGAGACCTTAGTGTATGTCTCCTGTAGGTCGCAGCTCCCTATGGAT--------------------------GAGTCAAGTGAAGGTCCTGAGACAA-------------------CAAGTCCTC----GGCTATGTGGGGGTGAGGG-------------ATGC----AG--------CTGGAACCTCAGGGA-TCTCTGT-AAGCAGTGGCATAAATGCTTGGCGG--GAGAGAGCATGTTAGAGCTCACACGACATAGGAAGCCACTGA--GACACTG")
        self.assertEqual(alignment.sequences[1].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[1].seq), 174210431)
        self.assertEqual(alignment.sequences[1].seq[16160203:16160203+443], "AGTCTTCATAAGTGGAAATATAAGTTTTAATTATTCCAGCCAAATCATCACAGTAAAAAGATGTTAATATTTATCTCCATACTCATCCTTACCAGATAGCCTTGTGCTCTTGGAATGTCGTACTTATGAAGTGAGCAAGCAGTATTTAGAATAGCATGAGATAGACCAGTCACCAAAGGTTGTTGATTTCTAGGAATGTCTCATGCTGGTTGCATCTGTGTGTGAACTGGGAAGTTGAGATCGCTTTCCACAGAAATTAAGCAAAATCCTGGAGGTAAGTTGTGGATAGGGCTGCTGTAGACTCCCTGACAGCGTAGGAACTAGAACTCTGAAGCAAGCCACTGGTGCCTAACAAATTCTTTCTAAAGTTGTAGAGTGGATCTGTGGGGGCTGACGTGGCCTCCAACGCCTAACGTGGTTGAAAAGGCCATGGACTTGTGTTC")
        self.assertEqual(alignment[1], "AGTCTTCATAAGTGGAAATATAAGTTTTAATTATTCCAGCCAAATCATCACAGTAAAAAGATGTTAATATTTATCTCCATACTCATCCTTACCAGATAGCCTTGTGCTCTTGGAATGTCGTACTTATGAAGTGAGCAAGCAGTATTTAGAATAGCATGAGATAGACCAGTCACCAA----AGGTTGTTGATTTCTAGGAATGTCTCATGCTGGTTGCATCTGTGTGTGAACTGGGAAGTTGAGATCGCTTTCCACAGAAATTAAGCAAAATCCTGGAGGTAAGTTGTGGA-----------TAGGGCTGCTGTAGACTCCCTGACAGCGTAGGAACT--------AGAACTCTGAAGCAAGCCACTGGTGCCTAACAAATTCTTTCTAAAGTTGTAGAGTGGATCTGTGGGGGCTGACGTGGCCTCCAACGCCTAACGTGGTTGAAAAGGCCATGGACTTGTGTTC")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'N')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[2].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[2].seq), 173908612)
        self.assertEqual(alignment.sequences[2].seq[16379004:16379004+443], "AGTCTTCATAAGTGGAAATATAAGTTTTAATTATTCCAGCTAAATCATCACAGTAAAAAGATGTTAATATTTATCTCCATACTCATCCTTACCAGATAGCCTTGTGCTCTTGGAATGTCGTATTTGTGAAGTGAGCAAGCAGTATGTAGAATAGCATGAGATAGACCAGTCACCAAAGGTTGTTGATTTCTAGGAATGTCTCATGCTGGTTGCATCTGTGTGTGAACTGGGAAGTTGAGATCGCTTTCCACAGAAATTAAGCAAAATCCTGGAGGTAAGTTGTGGATAGGGATGCTGTAGACTCCCTGACAGCATAGGAACTAGAACTCTGAAGCAAGCCACTGGTGCCTAACAAATTCTTTATAAAGTTGTAGAGTGGATCTGTGGGGGCTGACGTGGCCTCCAACGCCTAACATGGTTGAAAAGGCCATGGACTTGTGTTC")
        self.assertEqual(alignment[2], "AGTCTTCATAAGTGGAAATATAAGTTTTAATTATTCCAGCTAAATCATCACAGTAAAAAGATGTTAATATTTATCTCCATACTCATCCTTACCAGATAGCCTTGTGCTCTTGGAATGTCGTATTTGTGAAGTGAGCAAGCAGTATGTAGAATAGCATGAGATAGACCAGTCACCAA----AGGTTGTTGATTTCTAGGAATGTCTCATGCTGGTTGCATCTGTGTGTGAACTGGGAAGTTGAGATCGCTTTCCACAGAAATTAAGCAAAATCCTGGAGGTAAGTTGTGGA-----------TAGGGATGCTGTAGACTCCCTGACAGCATAGGAACTAGAACTCTGAAGC----AAGCCA----CTGGTGCCTAACAAATTCTTTATAAAGTTGTAGAGTGGATCTGTGGGGGCTGACGTGGCCTCCAACGCCTAACATGGTTGAAAAGGCCATGGACTTGTGTTC")
        self.assertEqual(alignment.column_annotations["panTro2.chr6"], "99999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999----99999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999-----------9999999999999999999999999999999999999999999999999----999999----999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999")
        self.assertEqual(alignment.sequences[2].annotations['leftStatus'], 'N')
        self.assertEqual(alignment.sequences[2].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[2].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[3].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[3].seq), 170899992)
        self.assertEqual(alignment.sequences[3].seq[15860456:15860456+443], "AGTCTTCATAAGTGGAAATATAAGTTTTAATTATTCCAGCTAAATCATCACAGTAAAAAGATGTTAATATTTATCTCCATACTCATCCTTACCAGATAGCCTTGTGCTCTTGGAATGTCGTATTTGTGAAGTGAGCAAGCAGTATGTAGAATAGCATGAGATAGACCAGTCACCAAAGGTTGTTGATTTCTAGGAATGTCTCATGCTGGTTGCATCTGTGTGTGAACTGGGAAGTTGAGATCGCTTTCCACAGAAATTAAGCAAAATCCTGGAGGTAAGTTGTGGATAGGGATGCTGTAGACTCCCTGACAGCATAGGAACTAGAACTCTGAAGCAAGCCACTGGTGCCTAACAAATTCTTTATAAAGTTGTAGAGTGGATCCATGGGGGCTGACGTGGCCTCCAACGCCTAACATGGTTGAAAAGGCCATGGACTTGTGTTC")
        self.assertEqual(alignment[3], "AGTCTTCATAAGTGGAAATATAAGTTTTAATTATTCCAGCTAAATCATCACAGTAAAAAGATGTTAATATTTATCTCCATACTCATCCTTACCAGATAGCCTTGTGCTCTTGGAATGTCGTATTTGTGAAGTGAGCAAGCAGTATGTAGAATAGCATGAGATAGACCAGTCACCAA----AGGTTGTTGATTTCTAGGAATGTCTCATGCTGGTTGCATCTGTGTGTGAACTGGGAAGTTGAGATCGCTTTCCACAGAAATTAAGCAAAATCCTGGAGGTAAGTTGTGGATAGGGATGCTGTAGACTCCCT---GACAGCATAGGAAC-TAGAACTCTGAAGC---AAGC----CA--------CTGGTGCCTAACAAATTCTTTATAAAGTTGTAGAGTGGATCCATGGGGGCTGACGTGGCCTCCAACGCCTAACATGGTTGAAAAGGCCATGGACTTGTGTTC")
        self.assertEqual(alignment.sequences[3].annotations['leftStatus'], 'N')
        self.assertEqual(alignment.sequences[3].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[3].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[3].annotations['rightCount'], 0)
        self.assertEqual(len(alignment.sequences), 4)
        self.assertNotIn("empty", alignment.annotations)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[ 3012076,  3012116,  3012116,  3012141,  3012141,  3012183,
                               3012183,  3012211,  3012211,  3012231,  3012235,  3012286,
                               3012286,  3012311,  3012311,  3012311,  3012320,  3012320,
                               3012320,  3012334,  3012335,  3012339,  3012339,  3012339,
                               3012339,  3012339,  3012343,  3012343,  3012345,  3012345,
                               3012345,  3012360,  3012360,  3012367,  3012367,  3012392,
                               3012392,  3012434,  3012434,  3012441],
                             [16160203, 16160243, 16160248, 16160273, 16160282, 16160324,
                              16160329, 16160357, 16160359, 16160379, 16160379, 16160430,
                              16160456, 16160481, 16160489, 16160489, 16160498, 16160499,
                              16160502, 16160516, 16160517, 16160521, 16160525, 16160525,
                              16160525, 16160526, 16160530, 16160534, 16160536, 16160540,
                              16160544, 16160559, 16160560, 16160567, 16160568, 16160593,
                              16160595, 16160637, 16160639, 16160646],
                             [16379004, 16379044, 16379049, 16379074, 16379083, 16379125,
                              16379130, 16379158, 16379160, 16379180, 16379180, 16379231,
                              16379257, 16379282, 16379290, 16379290, 16379299, 16379300,
                              16379303, 16379317, 16379318, 16379322, 16379326, 16379332,
                              16379334, 16379335, 16379339, 16379339, 16379341, 16379345,
                              16379345, 16379360, 16379361, 16379368, 16379369, 16379394,
                              16379396, 16379438, 16379440, 16379447],
                             [15860456, 15860496, 15860501, 15860526, 15860535, 15860577,
                              15860582, 15860610, 15860612, 15860632, 15860632, 15860683,
                              15860709, 15860734, 15860742, 15860753, 15860762, 15860763,
                              15860763, 15860777, 15860777, 15860781, 15860785, 15860791,
                              15860791, 15860791, 15860795, 15860795, 15860797, 15860797,
                              15860797, 15860812, 15860813, 15860820, 15860821, 15860846,
                              15860848, 15860890, 15860892, 15860899]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 49128)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3012441:3012441+125], "TGGGTCCCCTTGGCACATCCAGATCTCCCCAGTTAACCTGTCCTGCTTAGACCACTTACCTGAATTGAATTGGGAGGAGAGAAAGAAGCCAGTTTCCCAGAGAGGGAAAAGGAAAAGCTCGACAC")
        self.assertEqual(alignment[0], "TGGGTCCCCTTGGCACATCCAGATCTCCCCAGTTAACCTGTCCTGCTTAGACCACTTACCTGAATTG--AATTGGGAGGAGAGAAAGAAGCCAGTTTCCCAGAGAGGGAAAAGGAAAAGCTCGACAC")
        self.assertEqual(alignment.sequences[1].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[1].seq), 170899992)
        self.assertEqual(alignment.sequences[1].seq[15860899:15860899+114], "TGGGTTCCTCTAGAATAACCAAGTCCTCACGTAACCTGGTGTGTATTGACCACCTCTCTTGACCGCTGATCTTGGGGAGCACCTTGCTGAGGGCCAATGGAAAATGCTGGAAGC")
        self.assertEqual(alignment[1], "TGGGTTCCTCTAGAATAACCAAG--TCCTCACGTAACCTGGTGTGTATTGACCACCTCTCTTGACCGCTGATCTTGGGGAG----------CACCTTGCT-GAGGGCCAATGGAAAATGCTGGAAGC")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[2].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[2].seq), 173908612)
        self.assertEqual(alignment.sequences[2].seq[16379447:16379447+114], "TGGGTTCCTCTAGAATAACCAAGTCCTCACGTAACCTGGTGTGTATTGACCACCTCTCTTGACCGCTGATCTTGGGGAGCACCTTGCTGAGGGCCAATGGAAAATGCTGGAAGC")
        self.assertEqual(alignment[2], "TGGGTTCCTCTAGAATAACCAAG--TCCTCACGTAACCTGGTGTGTATTGACCACCTCTCTTGACCGCTGATCTTGGGGAG----------CACCTTGCT-GAGGGCCAATGGAAAATGCTGGAAGC")
        self.assertEqual(alignment.column_annotations["panTro2.chr6"], "99999999999999999999999--99999999999999999999999999999999999999999999999999999999----------999999999-99999999999999999999999999")
        self.assertEqual(alignment.sequences[2].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[2].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[3].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[3].seq), 174210431)
        self.assertEqual(alignment.sequences[3].seq[16160646:16160646+114], "TGGGTTCCTCTAGAATAACCAAGTCCTCACGTAACCTGGTCTATATTGACCACCTGTCTTGACTGTTGATCTTGGGGAGCACCTTGCTGAGGGCCAATGGAAAATGCTGGAAGC")
        self.assertEqual(alignment[3], "TGGGTTCCTCTAGAATAACCAAG--TCCTCACGTAACCTGGTCTATATTGACCACCTGTCTTGACTGTTGATCTTGGGGAG----------CACCTTGCT-GAGGGCCAATGGAAAATGCTGGAAGC")
        self.assertEqual(alignment.sequences[3].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[3].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[3].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[3].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[4].id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(alignment.sequences[4].seq), 359464)
        self.assertEqual(alignment.sequences[4].seq[180525:180525+101], "TGGGCCCCCCTGACCTGGCCAAGTCCTCACTCACCCGGGTTTACACTGACCACCAAACTTGGGGAGCCCTTTGCTGAGGGGCAACAGCAAACACGGGAAGC")
        self.assertEqual(alignment[4], "TGGGCCCCCCTGACCTGGCCAAG--TCCTCACTCACCCGGGTTTACACTGACCACC-------------AAACTTGGGGAG----------CCCTTTGCT-GAGGGGCAACAGCAAACACGGGAAGC")
        self.assertEqual(alignment.column_annotations["otoGar1.scaffold_334.1-359464"], "99999999999999999999999--9999999999999999999999999999999-------------999999999999----------999999999-99999979999999999999999999")
        self.assertEqual(alignment.sequences[4].annotations['leftStatus'], 'N')
        self.assertEqual(alignment.sequences[4].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[4].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[4].annotations['rightCount'], 0)
        self.assertEqual(len(alignment.sequences), 5)
        self.assertNotIn("empty", alignment.annotations)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[ 3012441,  3012464,  3012466,  3012497,  3012508,  3012508,
                               3012520,  3012530,  3012539,  3012540,  3012566],
                             [15860899, 15860922, 15860922, 15860953, 15860964, 15860966,
                              15860978, 15860978, 15860987, 15860987, 15861013],
                             [16379447, 16379470, 16379470, 16379501, 16379512, 16379514,
                              16379526, 16379526, 16379535, 16379535, 16379561],
                             [16160646, 16160669, 16160669, 16160700, 16160711, 16160713,
                              16160725, 16160725, 16160734, 16160734, 16160760],
                             [  180525,   180548,   180548,   180579,   180579,   180579,
                                180591,   180591,   180600,   180600,   180626]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 117109)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3012566:3012566+262], "TGTGGGCTCCTCACTTTCCTGTCTCAGGTGTGTCTGTGAGTTTCGGTGAGTGTCGTACAGGAAAGAGGGTGAAAACTCAGTCTGAGCTGTCATTCTTGCCAGCTATGTTGCTTTCCTGTCCTCTTTAGCTTATCTCAGGCAACCTATCTTATTTTGTTTGCTTTCAGAAGGCAAGCGAtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtatgtgtgtgtgcgtgCGCGCGCGCGAGCACATGTGCATGCATGCGCACTCGTG")
        self.assertEqual(alignment[0], "--TGTGGGCTCCTCACTTTCCTG-TCTCAGGTGTGTCTGTGAGTTTCGGTGAGTGTCGTACAGGAAAGAGGGTGAAAACTCAGTCTGAGCTGTCATTCTTGCCAGCTATGTTGCTTTCCTGTCCTCTTTAGC-------TTATCTCAGGCAACCTATCTTATTTTGTT-TGCTTTC--AGAAGGCAAG---CGAtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtatgtgtgtgtgcgtgCGCGCGCGCGAGCACATGTGCATGCATGCGCACTCGTG")
        self.assertEqual(alignment.sequences[1].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[1].seq), 170899992)
        self.assertEqual(alignment.sequences[1].seq[15861013:15861013+199], "CATGAGCCCTCCTTGTATCCCTGATCTTAGATGGGTCTATGGTTTTTGGTGGAAGCCATCCAGAAACAGCCTCAAAATCCAGTCTGAGTTTTCATTCCTGTCAATCAATCCTCTTATTTTTTTCCTTCACTACCTTATTTCTAGCAACACATCTTACTTATTCTGTTTTCCCAGAAGATGAGTTGTGTATGAGTTTTGG")
        self.assertEqual(alignment[1], "CATGAGCCCTCCTTGTATCCCTGATCTTAGATGGGTCTATGGTTTTTGGTGGAAGCCATCCAG-AAACAGCCTCAAAATCCAGTCTGAGTTTTCATTCCTGTCAATCAATCCTCTTATTTTTTTCCTTCACTACC----TTATTTCTAGCAACACATCTTAC-TTATTCTGTTTTCCCAGAAGA-------TGAGTTGTGTATGAGTTTTGG------------------------------------------------------------------")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[2].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[2].seq), 173908612)
        self.assertEqual(alignment.sequences[2].seq[16379561:16379561+199], "CATGAGCCCTCCTTGTATCCCTGATCTTAGATGGGTCTATGGTTTTTGGTGGAAGCCATCCAGAAACAGCCTCAAAATCCAGTCTGAGTTTTCATTCCTGTCAATCAATCCTCTTATTTTTTTCCTTCACTACCTTATTTCTAGCAACACATCTTACTTATTTTGTTTTCCCAGAAGATGAGTTGTGTATGAGTTTTGG")
        self.assertEqual(alignment[2], "CATGAGCCCTCCTTGTATCCCTGATCTTAGATGGGTCTATGGTTTTTGGTGGAAGCCATCCAG-AAACAGCCTCAAAATCCAGTCTGAGTTTTCATTCCTGTCAATCAATCCTCTTATTTTTTTCCTTCACTACC----TTATTTCTAGCAACACATCTTAC-TTATTTTGTTTTCCCAGAAGA-------TGAGTTGTGTATGAGTTTTGG------------------------------------------------------------------")
        self.assertEqual(alignment.column_annotations["panTro2.chr6"], "999999999999999999999999999999999999999999999999999999999999999-99999999999999999999999999999999999999999999999999999999999999999999999----99999999999999999999999-999999999999999999999-------999999999999999999999------------------------------------------------------------------")
        self.assertEqual(alignment.sequences[2].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[2].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[3].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[3].seq), 174210431)
        self.assertEqual(alignment.sequences[3].seq[16160760:16160760+199], "CATGAGCCCTCCTTGTATCCCTGATCTTAGATGGGTCTATGGTTTTTGGTGGAAGCCATCCAGAAACAGCCTCAAAATCCAGTCTGAGTTGTCATTCCTGTCAATCAATCCTCTTGTTTTTTTCCTTCATTACCTTATTTCTAGCAACATATCTTACTTATTTTGTTTTCCCAGAAGATGAGTTGTGTATGAGTTTTGG")
        self.assertEqual(alignment[3], "CATGAGCCCTCCTTGTATCCCTGATCTTAGATGGGTCTATGGTTTTTGGTGGAAGCCATCCAG-AAACAGCCTCAAAATCCAGTCTGAGTTGTCATTCCTGTCAATCAATCCTCTTGTTTTTTTCCTTCATTACC----TTATTTCTAGCAACATATCTTAC-TTATTTTGTTTTCCCAGAAGA-------TGAGTTGTGTATGAGTTTTGG------------------------------------------------------------------")
        self.assertEqual(alignment.sequences[3].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[3].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[3].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[3].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[4].id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(alignment.sequences[4].seq), 359464)
        self.assertEqual(alignment.sequences[4].seq[180626:180626+204], "CACATGCACTCATCGTGTCTCTGTCTTAGGCGTAGCCCTGGTTTCTGGGGAAAGCTGTCCAGAAAACCACCTGAAAACCCAGTCTCTGTGGTCATTGCTGGCAATCAGTCCCTTATTTTTAGGTTGCCTTGTTATTTCAGGTGACATAGCTTACTTTGTTTTGTTTTCCCAGAAGACAAGTGGTGACTTGTGTATGAGTTTTGA")
        self.assertEqual(alignment[4], "CACATGCACTCATCGTGTCTCTG-TCTTAGGCGTAGCCCTGGTTTCTGGGGAAAGCTGTCCAGAAAACCACCTGAAAACCCAGTCTCTGTGGTCATTGCTGGCAATCAGTCC-CTTATTTTTAGGTTGCCTTG------TTATTTCAGGTGACATAGCTTACTTTGTTTTGTTTTCCCAGAAGACAAGTGGTGACTTGTGTATGAGTTTTGA------------------------------------------------------------------")
        self.assertEqual(alignment.column_annotations["otoGar1.scaffold_334.1-359464"], "99999999999999999999999-9999999999999999999999999999999999999999999999999999999999999999999999999999999999999999-99999999999999999999------9999999999999999999999999999999999999999999999999999999999999999999999999------------------------------------------------------------------")
        self.assertEqual(alignment.sequences[4].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[4].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[4].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[4].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[5].id, "cavPor2.scaffold_290371")
        self.assertEqual(len(alignment.sequences[5].seq), 39932)
        self.assertEqual(alignment.sequences[5].seq[310:310+205], "TGCGTATCCCTTGCGATCCTGATCCGAGGTGTGGCTGTGGTTTTTGGTGAAAGCCATCCAGAAAAGAACTTCAAAACCTAGTCTAAATTGTTATTCCCTCCTGCCAATCCTCTTGTTTTTACCCCATAGTGACCTCGTTATTGCAAAACATatcttattttgttctgttttcctggaagataaatgaattttgtgtgaGTTTGAG")
        self.assertEqual(alignment[5], "--TGCGTATCCCTTGCGATCCTGATCCGAGGTGTGGCTGTGGTTTTTGGTGAAAGCCATCCAGAAAAGAACTTCAAAACCTAGTCTAAATTGTTATTCCCTCCTGCCAATCCTCTTGTTTTTACCCCATAGTGACCTCGTTATTGCAA--AACATatcttattttgttctgttttcctggaagataaa---tgaattttgtgtgaGTTTGAG------------------------------------------------------------------")
        self.assertEqual(alignment.column_annotations["cavPor2.scaffold_290371"], "--99999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999997999999999999999999999999999999999999999--99999999999999999999999999999999999999---999999999999999999999------------------------------------------------------------------")
        self.assertEqual(alignment.sequences[5].annotations['leftStatus'], 'N')
        self.assertEqual(alignment.sequences[5].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[5].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[5].annotations['rightCount'], 0)
        self.assertEqual(len(alignment.sequences), 6)
        self.assertNotIn("empty", alignment.annotations)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[ 3012566,  3012566,  3012587,  3012587,  3012626,  3012627,
                               3012675,  3012676,  3012695,  3012695,  3012695,  3012695,
                               3012704,  3012706,  3012718,  3012719,  3012724,  3012724,
                               3012731,  3012731,  3012737,  3012741,  3012741,  3012762,
                               3012828],
                             [15861013, 15861015, 15861036, 15861037, 15861076, 15861076,
                              15861124, 15861125, 15861144, 15861145, 15861147, 15861147,
                              15861156, 15861158, 15861170, 15861170, 15861175, 15861176,
                              15861183, 15861185, 15861191, 15861191, 15861191, 15861212,
                              15861212],
                             [16379561, 16379563, 16379584, 16379585, 16379624, 16379624,
                              16379672, 16379673, 16379692, 16379693, 16379695, 16379695,
                              16379704, 16379706, 16379718, 16379718, 16379723, 16379724,
                              16379731, 16379733, 16379739, 16379739, 16379739, 16379760,
                              16379760],
                             [16160760, 16160762, 16160783, 16160784, 16160823, 16160823,
                              16160871, 16160872, 16160891, 16160892, 16160894, 16160894,
                              16160903, 16160905, 16160917, 16160917, 16160922, 16160923,
                              16160930, 16160932, 16160938, 16160938, 16160938, 16160959,
                              16160959],
                             [  180626,   180628,   180649,   180649,   180688,   180689,
                                180737,   180737,   180756,   180757,   180757,   180757,
                                180766,   180768,   180780,   180781,   180786,   180787,
                                180794,   180796,   180802,   180806,   180809,   180830,
                                180830],
                             [     310,      310,      331,      332,      371,      372,
                                   420,      421,      440,      441,      443,      447,
                                   456,      456,      468,      469,      474,      475,
                                   482,      484,      490,      494,      494,      515,
                                   515]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 128047)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3012828:3012828+168], "TTTGCATAGACTCTCTTGGCAACAAAATAACGTTATATTTAAACATCCATTAAAATAATGCACTTAGCACAGCCTGCCCTGAGGGATGAACACTATTGTTAAAGAACTATTCCGCTAAGGCAGCAACCTCTGGATCTTCAGCATTCTGGCGCCATCTGCTGGTCATAT")
        self.assertEqual(alignment[0], "TTTGCATAGACTCTCTTGGCAACAAAATAACGTTATATTTAAACATCCATTAAAATAATGCACTTAGCACAGCCTGCCCTGAGGGAT----GAACACT--ATTGTTAA-AGAACTATTCCGCTAAGGCAGCAACCTCTGGATCTTCAGCATTCTGGCGCCATCTGCTGGTCATAT")
        self.assertEqual(alignment.sequences[1].id, "cavPor2.scaffold_290371")
        self.assertEqual(len(alignment.sequences[1].seq), 39932)
        self.assertEqual(alignment.sequences[1].seq[515:515+166], "GTGGGATGTCTTGACAGCTAAATAGTcttatatttaagcatttatGAAAATAATGAGGTTATCATGGTTCGCTCTAGTAGCTCTGAGACCTTTATTGTCAGCAAAGCTGTTTAGCAGGAGGCAACTTCAGGTTCCGTGGGATTCTAGTGCCATCTGCTGGTCATAT")
        self.assertEqual(alignment[1], "-----GTGGGATGTCTTGACAGCTAAATAGTcttatatttaagcatttatGAAAATAATGAGGTTATCATGGTTCGCTCTAGTAGCTCTGAGACCTTT--ATTGTCAGCAAAGCTGTTTAGC--AGGAGGCAACTTCAGGTTCCGTGGGATTCTAGTGCCATCTGCTGGTCATAT")
        self.assertEqual(alignment.column_annotations["cavPor2.scaffold_290371"], "-----999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999--9999999999999999999999--999999999999999999999999999999999999999999999999999")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[2].id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(alignment.sequences[2].seq), 359464)
        self.assertEqual(alignment.sequences[2].seq[180830:180830+161], "GTGGACATTCGTGACAGGCAAATAGTGTTCTATTTAAACATTTGTTGAAATGATAAACTTAGCTCAGCCCCTTAGAAGTGAGAGGCTACTGTCAGGGAAGCGGTTTTGCGGAGCAGCAGCCTCAGGATCTGTGGGATTTTAGCGCCATCTGCTGGTCACAG")
        self.assertEqual(alignment[2], "-----GTGGACATTCGTGACAGGCAAATAGTGTTCTATTTAAACATTTGTTGAAATGATAAACTTAGCTCAGCC--CCTT---AGAAGTG-AGAGGCT--ACTGTCAGGGAAGCGGTTTTGCG-GAGCAGCAGCCTCAGGATCTGTGGGATTTTAGCGCCATCTGCTGGTCACAG")
        self.assertEqual(alignment.column_annotations["otoGar1.scaffold_334.1-359464"], "-----999999999999999999999999999999999999999999999999999999999999999999999--9999---9999999-9999999--99999999999999999999989-996999999999999999999999999999999999999999999999999")
        self.assertEqual(alignment.sequences[2].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[2].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[3].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[3].seq), 174210431)
        self.assertEqual(alignment.sequences[3].seq[16160959:16160959+161], "GTGGACTATCCTGACAACCAAATAGCGTTATATTTAAACATTTATGAAAATAATAAGCTTAGCATAACCTGCCTTGACAGATGTGGTATGCCTAGTGAAGTTGTTTAGCTAAGGCTCAATCTCAGGATCTGTGGGACTCTAATGCCATCTGCTGGTCACAT")
        self.assertEqual(alignment[3], "-----GTGGACTATCCTGACAACCAAATAGCGTTATATTTAAACATTTATGAAAATAATAAGCTTAGCATAACCTGCCTTGACAGATGTG-GTATGC-------CTAGTGAAGTTGTTTAGCT-AAGGCTCAATCTCAGGATCTGTGGGACTCTAATGCCATCTGCTGGTCACAT")
        self.assertEqual(alignment.sequences[3].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[3].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[3].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[3].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[4].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[4].seq), 173908612)
        self.assertEqual(alignment.sequences[4].seq[16379760:16379760+161], "GCGGACTATCTTGACAACCAAATAGCGTTATATTTAAACGTTTATGAAAATAATAAGCTTAGCATAACCTGCCTTGACAGATGTGGTATGCCTAGTGAAGTTGCTTAGCTATGGCTCAATCTCAGGATCTGTGGGACTCTAATGCCATCTGCTGGTCACAC")
        self.assertEqual(alignment[4], "-----GCGGACTATCTTGACAACCAAATAGCGTTATATTTAAACGTTTATGAAAATAATAAGCTTAGCATAACCTGCCTTGACAGATGTG-GTATGC-------CTAGTGAAGTTGCTTAGCT-ATGGCTCAATCTCAGGATCTGTGGGACTCTAATGCCATCTGCTGGTCACAC")
        self.assertEqual(alignment.column_annotations["panTro2.chr6"], "-----9999999999999999999999999999999999999999999999999999999999999999999999999999999999999-999999-------9999999999999999999-999999999999999999999999999999999999999999999999999")
        self.assertEqual(alignment.sequences[4].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[4].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[4].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[4].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[5].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[5].seq), 170899992)
        self.assertEqual(alignment.sequences[5].seq[15861212:15861212+161], "GTGGACTATCTTGACAACCAAATAGCGTTATATTTAAACGTTTGTGAAAATAATAAGCTTAGCATAACCTGCCTTGACAGATGTGGTATGCCTAGTGAAGTTGCTTAGCTATGGCTCAATCTCAGGATCTGTGGGACTCTAATGCCATCTGCTGGTCACAC")
        self.assertEqual(alignment[5], "-----GTGGACTATCTTGACAACCAAATAGCGTTATATTTAAACGTTTGTGAAAATAATAAGCTTAGCATAACCTGCCTTGACAGATGTG-GTATGC-------CTAGTGAAGTTGCTTAGCT-ATGGCTCAATCTCAGGATCTGTGGGACTCTAATGCCATCTGCTGGTCACAC")
        self.assertEqual(alignment.sequences[5].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[5].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[5].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[5].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[6].id, "echTel1.scaffold_288249")
        self.assertEqual(len(alignment.sequences[6].seq), 100002)
        self.assertEqual(alignment.sequences[6].seq[87492:87492+169], "TTGGGGTGGATGCTCTTGGCAGTCACACAGTGCTCTATTTTAGGATTTACTAGAACAATGAGTTTGTCATAACTGCCTCCTCCCAAGTGGGAAGCTGAAGCACCAAGTACTGACTCAGCAGTCCCTGACCTCACATCCATGGAAATCTAGTGACATCTGCTGGACACAT")
        self.assertEqual(alignment[6], "TTGGGGTGGATGCTCTTGGCAGTCACACAGTGCTCTATTTTAGGATTTACTAGAACAATGAGTTTGTCATAACT-GCCTCCTCCCAAGTG-GGAAGCTGAAGCACCAA-GTACTGACTCAGC--AGTCCCTGACCTCA-CATCCATGGAAATCTAGTGACATCTGCTGGACACAT")
        self.assertEqual(alignment.column_annotations["echTel1.scaffold_288249"], "99989999998979665899999999999676899897997899979878899999999999999999999999-999999999999999-99999999999999999-9999999999999--99999999999999-999999999999999999999999999999999999")
        self.assertEqual(alignment.sequences[6].annotations['leftStatus'], 'N')
        self.assertEqual(alignment.sequences[6].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[6].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[6].annotations['rightCount'], 7564)
        self.assertEqual(len(alignment.sequences), 7)
        self.assertNotIn("empty", alignment.annotations)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[ 3012828,  3012833,  3012902,  3012903,  3012904,  3012908,
                               3012911,  3012915,  3012915,  3012915,  3012921,  3012922,
                               3012922,  3012926,  3012930,  3012930,  3012943,  3012944,
                               3012945,  3012959,  3012960,  3012996],
                             [     515,      515,      584,      585,      586,      590,
                                   593,      597,      600,      601,      607,      608,
                                   608,      612,      616,      617,      630,      630,
                                   630,      644,      645,      681],
                             [  180830,   180830,   180899,   180899,   180899,   180903,
                                180903,   180907,   180910,   180910,   180916,   180917,
                                180917,   180921,   180925,   180926,   180939,   180940,
                                180940,   180954,   180955,   180991],
                             [16160959, 16160959, 16161028, 16161029, 16161030, 16161034,
                              16161037, 16161041, 16161044, 16161044, 16161050, 16161050,
                              16161050, 16161050, 16161054, 16161055, 16161068, 16161069,
                              16161069, 16161083, 16161084, 16161120],
                             [16379760, 16379760, 16379829, 16379830, 16379831, 16379835,
                              16379838, 16379842, 16379845, 16379845, 16379851, 16379851,
                              16379851, 16379851, 16379855, 16379856, 16379869, 16379870,
                              16379870, 16379884, 16379885, 16379921],
                             [15861212, 15861212, 15861281, 15861282, 15861283, 15861287,
                              15861290, 15861294, 15861297, 15861297, 15861303, 15861303,
                              15861303, 15861303, 15861307, 15861308, 15861321, 15861322,
                              15861322, 15861336, 15861337, 15861373],
                             [   87492,    87497,    87566,    87566,    87567,    87571,
                                 87574,    87578,    87581,    87581,    87587,    87588,
                                 87590,    87594,    87598,    87598,    87611,    87611,
                                 87611,    87625,    87625,    87661]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 98097)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3012996:3012996+222], "AGATGTCTGCTGTGGAGACCTGGCCAACTTTGCTTTCTTCAAAAAGGCAACAGAAGGTAATCAGTTGAATGCCCACCATTAGGAAGGCGACCTCTAGTGCACAAACCTTGACATTTTCCCTTTTAATGGAATTTAACAGAAGTTCAGGATGTTCTTTGGGTAATTTACAATTAGGGGGCAAAAATCAAAAGTATTTCGAGCATATCAAAACTGTTAGCTATG")
        self.assertEqual(alignment[0], "AGATGTCTGCTGTGGAGA-------CCTGGCCAACTTTG----CTT--TCTTC-----AAAAAGGCAACAGAAGGTAATCAGTTGAATGCCCACCA-----TTAGGAAGGCGACCTCTAGTGCACAAACCTTGAC-ATTTTCCCTTTTAATGGAA-TTTAACAGAAGTTCAGGATGTTCTTTGGGTAATTTACAATT---A----GGGGGCAAAAATCAAAAGTATTTCGAGCATATCAAAACTGTTAGCTATG")
        self.assertEqual(alignment.sequences[1].id, "cavPor2.scaffold_290371")
        self.assertEqual(len(alignment.sequences[1].seq), 39932)
        self.assertEqual(alignment.sequences[1].seq[681:681+177], "TGACGTCTGCAGCCACGGGCTGGCTAACTCCAGTCTTCCGCGCTAAAGCCGACCACACCACATAGTGGGGAGAGTGACCTCTAGCGGAAAagctttgcaatttttttcttcgACTCTGTCAGGTTCAGAATTCCATTTAAACAATTTacagctagaaaaggaaaagcagaggCATTT")
        self.assertEqual(alignment[1], "TGACGTCTGCAGCCACGG-------GCTGGCTAACTCCA--GTCTT--CCGCG-----CTAAAGCCGAC--------------------CACACCACATAGTGGGGAGAGTGACCTCTAGCGGAAAagctttgca-atttttttcttc----gAC-TCTGTCAG--GTTCAGAATTCCATTTAAACAATTTacagct---a----gaaaaggaaaagcagaggCATTT--------------------------")
        self.assertEqual(alignment.column_annotations["cavPor2.scaffold_290371"], "999999999999999999-------99999999999999--99999--99999-----99999999999--------------------9999999999999999999999999999999999999999999999-999999999999----999-99999999--9999999999999999997999999999995---6----99999999999336991774687--------------------------")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'N')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[2].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[2].seq), 170899992)
        self.assertEqual(alignment.sequences[2].seq[15861373:15861373+184], "TGATGGCTGCTGCCATCACTTTCTTGGCCTACTCTGCTGTCTTACTCTATAAGACTAACAGGACCCCATAGTTAGGAGAGCGACCTCTAGCGAGCAGGCCTTGGCATTTTTCTTTTTAATTGAATTCATAAAGGCTTAGAATGCCATTTAGACAATTTGTAATGGTGGGGGGAGGTAAGAGTAT")
        self.assertEqual(alignment[2], "TGATGGCTGCTGCCATCA---CTTTCTTGGCCTACTCTGCTGTCTT--ACTCT-----ATAAGACTAAC--------------------AGGACCCCATAGTTAGGAGAGCGACCTCTAGCGAGCAGGCCTTGGC-ATTTTTCTTTTTAATTGAA-TTCATAAA-GGCTTAGAATGCCATTTAGACAATTTGTAATG---GTGGGGGGAGGTAAGAGTAT----------------------------------")
        self.assertEqual(alignment.sequences[2].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[2].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[2].annotations['rightCount'], 14)
        self.assertEqual(alignment.sequences[3].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[3].seq), 173908612)
        self.assertEqual(alignment.sequences[3].seq[16379921:16379921+184], "TGATGGCTGCTGCCATCACTTTCTTGGCCTACTCTGCTGTCTTATTCTATAAGACTAACAGGACCCCATGGTTAGGAGAGCGACCTCTAGCGAGCAGGCCTTGGAATTTTTCTTTTTAATTAAATTCATAAAGGCTTAGAATGCCATTTAGACAATTTGTAATGGTGGGGGGAGGTAAGAGTAT")
        self.assertEqual(alignment[3], "TGATGGCTGCTGCCATCA---CTTTCTTGGCCTACTCTGCTGTCTT--ATTCT-----ATAAGACTAAC--------------------AGGACCCCATGGTTAGGAGAGCGACCTCTAGCGAGCAGGCCTTGGA-ATTTTTCTTTTTAATTAAA-TTCATAAA-GGCTTAGAATGCCATTTAGACAATTTGTAATG---GTGGGGGGAGGTAAGAGTAT----------------------------------")
        self.assertEqual(alignment.column_annotations["panTro2.chr6"], "999999999999999999---9999999999999999999999999--99999-----99999999999--------------------9999999999999999999999999999999999999999999999-9999999999999999999-99999999-99999999999999999999999999999999---99999999999999999999----------------------------------")
        self.assertEqual(alignment.sequences[3].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[3].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[3].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[3].annotations['rightCount'], 14)
        self.assertEqual(alignment.sequences[4].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[4].seq), 174210431)
        self.assertEqual(alignment.sequences[4].seq[16161120:16161120+217], "TGATGGCTGCTGCCATCACTTTCTTGGCCTACTAGGCCGTCTTACTCTATAAGACTAACAGGACCCCATAGTTAGGAGAGCGACCTCTAGCGAGCAGGCCTTGGAATTTTTCTTTTTAATTGAATTCATAAAGGCTTAGAATGCCATTTAGACAATTTGTAATGGTGGGGGGAGGTAAGAGTATAAATATTGAATGTATTCAAACCTGTCAATGATG")
        self.assertEqual(alignment[4], "TGATGGCTGCTGCCATCA---CTTTCTTGGCCTACTAGGCCGTCTT--ACTCT-----ATAAGACTAAC--------------------AGGACCCCATAGTTAGGAGAGCGACCTCTAGCGAGCAGGCCTTGGA-ATTTTTCTTTTTAATTGAA-TTCATAAA-GGCTTAGAATGCCATTTAGACAATTTGTAATG---GTGGGGGGAGGTAAGAGTATAAATATTGAATGTAT-TCAAACCTGTCAATGATG")
        self.assertEqual(alignment.sequences[4].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[4].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[4].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[4].annotations['rightCount'], 2)
        self.assertEqual(alignment.sequences[5].id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(alignment.sequences[5].seq), 359464)
        self.assertEqual(alignment.sequences[5].seq[180991:180991+226], "GGACGGCCGCCGCCACCACCGCTTTCCCTGCCAACTCCACCGTCTTAAATTCTGCAACAAAAGAAAAATCGGACCCTGCAGTTACTAGAGCGACCTCTGGTGAACAAGCCTTGGATTTTTTTCTTTTTAATTGAATTTTACAACGGCTCAAAATGCAAGTTAGACAATTCATAGTGTATGTGTGGGGACCTGAGAGTGCAAATATTGCTGGCATTCAGCCTTGTTA")
        self.assertEqual(alignment[5], "GGACGGCCGCCGCCACCACCGCTTTCCCTGCCAACTCCACCGTCTTAAATTCTGCAACAAAAGAAAAAT--------------------CGGACCCTGCAGTTACTAGAGCGACCTCTGGTGAACAAGCCTTGGATTTTTTTCTTTTTAATTGAATTTTACAAC-GGCTCAAAATGCAAGTTAGACAATTCATAGTGTATGTGTGGGGACCTGAGAGTGCAAATATTGCTGGCAT-TCAGCCTTGTTA------")
        self.assertEqual(alignment.column_annotations["otoGar1.scaffold_334.1-359464"], "999998999999999999999999999999999999999999999999999999999999999999999--------------------999999999999999999999999999999999999999999999999999999999999999999999999999-9999999999999999999999999999999999999999999699999999999999999999999999-999999999999------")
        self.assertEqual(alignment.sequences[5].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[5].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[5].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[5].annotations['rightCount'], 2931)
        self.assertEqual(alignment.sequences[6].id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(alignment.sequences[6].seq), 498454)
        self.assertEqual(alignment.sequences[6].seq[167217:167217+159], "TGATGGCTACAGCTGTCACTTTCTTGGACCACTCAGCTCTCTCAACGAAAGACAAACTGGATCCTCCGGGAGAGGAGACCTCTAGTGAGCAAACCTTGGAATTTTTCTTTAAAACAGAATTTATAAAGGTTCAAATGCCATTTAGACAACTTATAATTA")
        self.assertEqual(alignment[6], "TGATGGCTACAGCTGTCA---CTTTCTTGGACCACTCAGCTCTCTC--AAC-------GAAAGACAAAC--------------------TGGATCCTCCGG---GAGAGGAGACCTCTAGTGAGCAAACCTTGGA-ATTTTTCTTTAAAACAGAA-TTTATAAA-GGTTCA-AATGCCATTTAGACAACTTATAATT---A-----------------------------------------------------")
        self.assertEqual(alignment.column_annotations["tupBel1.scaffold_114895.1-498454"], "234233433332122232---1582112220212134443324332--133-------23732111754--------------------365002326236---2411112335242535353245936522224-1376645373578253554-58324573-545454-8444565585455465799967999---9-----------------------------------------------------")
        self.assertEqual(alignment.sequences[6].annotations['leftStatus'], 'N')
        self.assertEqual(alignment.sequences[6].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[6].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[6].annotations['rightCount'], 4145)
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
                numpy.array([[ 3012996,  3013014,  3013014,  3013014,  3013028,  3013028,
                               3013028,  3013031,  3013031,  3013034,  3013036,  3013036,
                               3013047,  3013067,  3013074,  3013074,  3013077,  3013108,
                               3013108,  3013120,  3013124,  3013127,  3013127,  3013135,
                               3013136,  3013137,  3013142,  3013143,  3013168,  3013168,
                               3013169,  3013169,  3013184,  3013192,  3013199,  3013200,
                               3013212,  3013218],
                             [     681,      699,      699,      699,      713,      713,
                                   715,      718,      718,      721,      723,      723,
                                   734,      734,      741,      746,      749,      780,
                                   780,      792,      792,      795,      795,      803,
                                   803,      803,      808,      809,      834,      834,
                                   835,      835,      850,      858,      858,      858,
                                   858,      858],
                             [15861373, 15861391, 15861391, 15861395, 15861409, 15861411,
                              15861413, 15861416, 15861416, 15861419, 15861421, 15861421,
                              15861432, 15861432, 15861439, 15861444, 15861447, 15861478,
                              15861478, 15861490, 15861494, 15861497, 15861497, 15861505,
                              15861505, 15861506, 15861511, 15861512, 15861537, 15861537,
                              15861538, 15861542, 15861557, 15861557, 15861557, 15861557,
                              15861557, 15861557],
                             [16379921, 16379939, 16379939, 16379943, 16379957, 16379959,
                              16379961, 16379964, 16379964, 16379967, 16379969, 16379969,
                              16379980, 16379980, 16379987, 16379992, 16379995, 16380026,
                              16380026, 16380038, 16380042, 16380045, 16380045, 16380053,
                              16380053, 16380054, 16380059, 16380060, 16380085, 16380085,
                              16380086, 16380090, 16380105, 16380105, 16380105, 16380105,
                              16380105, 16380105],
                             [16161120, 16161138, 16161138, 16161142, 16161156, 16161158,
                              16161160, 16161163, 16161163, 16161166, 16161168, 16161168,
                              16161179, 16161179, 16161186, 16161191, 16161194, 16161225,
                              16161225, 16161237, 16161241, 16161244, 16161244, 16161252,
                              16161252, 16161253, 16161258, 16161259, 16161284, 16161284,
                              16161285, 16161289, 16161304, 16161312, 16161319, 16161319,
                              16161331, 16161337],
                             [  180991,   181009,   181012,   181016,   181030,   181032,
                                181034,   181037,   181039,   181042,   181044,   181049,
                                181060,   181060,   181067,   181072,   181075,   181106,
                                181107,   181119,   181123,   181126,   181127,   181135,
                                181135,   181136,   181141,   181142,   181167,   181170,
                                181171,   181175,   181190,   181198,   181205,   181205,
                                181217,   181217],
                             [  167217,   167235,   167235,   167239,   167253,   167255,
                                167257,   167260,   167260,   167263,   167263,   167263,
                                167274,   167274,   167281,   167286,   167286,   167317,
                                167317,   167329,   167333,   167336,   167336,   167344,
                                167344,   167345,   167350,   167350,   167375,   167375,
                                167376,   167376,   167376,   167376,   167376,   167376,
                                167376,   167376]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3013218:3013218+219], "agccaggcgtggtggcacacacctttactcccagcatttggggggcagaggcaggtggatctgtgagtttgaggccagcctggtctacagagggagtctcaggacagccagagctacacagaaataacctgcctagaaaaacaaaacaaaacaaaacatcaaaactcaaaacaaaTAAAAAAAATAAAAAACCCAACCTAAACCAAATAACAAAACACT")
        self.assertEqual(alignment[0], "agccaggcgtggtggcacacacctttactcccagcatttggggggcagaggcaggtggatctgtgagtttgaggccagcctggtctacagagggagtctcaggacagccagagctacacagaaataacctgcctagaaaaacaaaacaaaacaaaacatcaaaactcaaaacaaaTAAAAAAAATAAAAAACCCAACCTAAACCAAATAACAAAACACT")
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
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[3013218, 3013437]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 40604)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3013437:3013437+166], "TCCAAAATGGTTAGCTATGCCCAACTCCTTTCACTCCAAGAAAATATCCTAACCATGTAAGAGAGCTAGCCTGTTGGTGGCAGCCAAGCCTGATGGTGGCAGACTAGATTGATGGTGCCAGACTACTTTATGGCTGTATCATTTTCCATTCATGTGTTGTGTTATA")
        self.assertEqual(alignment[0], "TCCAAAATGGTTAGCTATGCCCAACTCCTTTCACTCCAAGAAAATATCCTAACCATGTAAGAGAGCTAGCCTGTTGGTGGCAGCCAAGCCTGATGGTGGCAGACTAGATTGATGGTGCCAGACTACTTTATGGCTGTATCATTTTCCATTCATGTGTTGTGTTATA")
        self.assertEqual(alignment.sequences[1].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[1].seq), 173908612)
        self.assertEqual(alignment.sequences[1].seq[16380119:16380119+130], "TTCAAACCTGTCAATGATGTTCAACTCTTTTCACTGTGGATAAAGATCTTCTAGATGCAATAGAATTAGATTTTGGCGATTACATCCTATTTTATGGGTGTATCGTTTTGCATTCATGGGCTGTTTGATA")
        self.assertEqual(alignment[1], "TTCAAACCTGTCAATGATGTTCAACTCTTTTCACTGTGGATAAAGATCTTCTAGATGCAATAGAATTAGATT------------------------------------TTGGCGATTACATCCTATTTTATGGGTGTATCGTTTTGCATTCATGGGCTGTTTGATA")
        self.assertEqual(alignment.column_annotations["panTro2.chr6"], "999999999999999999999999999999999999999999999999999999999999999999999999------------------------------------9999999999999999999999999999999999999099999999999999999999")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 14)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 9106)
        self.assertEqual(alignment.sequences[2].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[2].seq), 170899992)
        self.assertEqual(alignment.sequences[2].seq[15861571:15861571+130], "TTCAAACCTGTCAATGATGTTCAACTCTTTTCACTGTGGATAAAGATCTTCTAGATGCAATAGAATTAGATTTTGGCGATTACATCCTATTTTATGGGTGTATCGTTTTGCATTCATGGGCTGTTTGATA")
        self.assertEqual(alignment[2], "TTCAAACCTGTCAATGATGTTCAACTCTTTTCACTGTGGATAAAGATCTTCTAGATGCAATAGAATTAGATT------------------------------------TTGGCGATTACATCCTATTTTATGGGTGTATCGTTTTGCATTCATGGGCTGTTTGATA")
        self.assertEqual(alignment.sequences[2].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[2].annotations['leftCount'], 14)
        self.assertEqual(alignment.sequences[2].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[2].annotations['rightCount'], 9085)
        self.assertEqual(alignment.sequences[3].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[3].seq), 174210431)
        self.assertEqual(alignment.sequences[3].seq[16161339:16161339+109], "CAACTCTTTTCACTGTGGATAAAGATCTTCTAGATGCAATAGAATTAGATTTTGGCGATTACATCCTATTTTATGGGTGTATCGTTTTGCATTCATGGGCTGTTTGATA")
        self.assertEqual(alignment[3], "---------------------CAACTCTTTTCACTGTGGATAAAGATCTTCTAGATGCAATAGAATTAGATT------------------------------------TTGGCGATTACATCCTATTTTATGGGTGTATCGTTTTGCATTCATGGGCTGTTTGATA")
        self.assertEqual(alignment.sequences[3].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[3].annotations['leftCount'], 2)
        self.assertEqual(alignment.sequences[3].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[3].annotations['rightCount'], 8044)
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
                numpy.array([[ 3013437,  3013458,  3013509,  3013545,  3013603],
                             [16380119, 16380140, 16380191, 16380191, 16380249],
                             [15861571, 15861592, 15861643, 15861643, 15861701],
                             [16161339, 16161339, 16161390, 16161390, 16161448]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3013603:3013603+1041], "CCTTCTTAAGAACACTAGACTCAggactggggagatggctcagcagttaagaatcggtgctgttaagagtgggagacaagtttggttcccagctcccacattggtcagctcacagccacccgtaactctaagatggtacacacctttaatcccaggagacagaggcaatcagatctgagttcaagattcagcctgagacagagcatgttccaaattcaggcatggtgggtcatacctttaatatgggacataccttctgctggaggcctacctaaggacaacggagaaaggaagtattcgttcttctcctgcttgcacttacttgccagcgcatctactggaacccacttcttcaggattccagcttatacaggagaccagctgaaatatccagcctctcgggactgaacaagtactagagtctcagacttcccattcacagctgcccattgttggttggttgtactacagactgtaagtcattgtaataatttcccttaatatatagagacattatataagttctgtgactctagagaaccctgactagtacaCGTGGCTAACTAGAAAGctctggtatgtgcttacttaatgctgaggttttaggcatggccacggtgctctgcttcttatgtgggtgctgggaatgcagactcaggtcctcatgtgtatgcagcaaacacttcatacactcagctgcttccctaacccTATGCTTGTGTCTTATTACTAACTTGTGAAAAGCTTTGAGTTTATTTTCTATGTTTTCAACCACTTTCTTGAGTATGCTCAGCTCGTGGCTTTAAACTGGATTTCCCCCTAATATGTAATGACTATAAGTATTCCTTAAATAGGACACACTTTTGTTATACTTTTTGTTATCatataaaatatttcaaaaaaatttttttGCTATTTTTATCTTTGAGCCATTGGTCATTTTGACGTGTATCTCTTGATTTTTATAGATGGTAATATTTTATGTATTGCTAGCCAATCTCGTTTTCTTGTTTGCTTGCTTGTTTGTTTTGGTCAATGCAG")
        self.assertEqual(alignment[0], "CCTTCTTAAGAACACTAGACTCAggactggggagatggctcagcagttaagaatcggtgctgttaagagtgggagacaagtttggttcccagctcccacattggtcagctcacagccacccgtaactctaagatggtacacacctttaatcccaggagacagaggcaatcagatctgagttcaagattcagcctgagacagagcatgttccaaattcaggcatggtgggtcatacctttaatatgggacataccttctgctggaggcctacctaaggacaacggagaaaggaagtattcgttcttctcctgcttgcacttacttgccagcgcatctactggaacccacttcttcaggattccagcttatacaggagaccagctgaaatatccagcctctcgggactgaacaagtactagagtctcagacttcccattcacagctgcccattgttggttggttgtactacagactgtaagtcattgtaataatttcccttaatatatagagacattatataagttctgtgactctagagaaccctgactagtacaCGTGGCTAACTAGAAAGctctggtatgtgcttacttaatgctgaggttttaggcatggccacggtgctctgcttcttatgtgggtgctgggaatgcagactcaggtcctcatgtgtatgcagcaaacacttcatacactcagctgcttccctaacccTATGCTTGTGTCTTATTACTAACTTGTGAAAAGCTTTGAGTTTATTTTCTATGTTTTCAACCACTTTCTTGAGTATGCTCAGCTCGTGGCTTTAAACTGGATTTCCCCCTAATATGTAATGACTATAAGTATTCCTTAAATAGGACACACTTTTGTTATACTTTTTGTTATCatataaaatatttcaaaaaaatttttttGCTATTTTTATCTTTGAGCCATTGGTCATTTTGACGTGTATCTCTTGATTTTTATAGATGGTAATATTTTATGTATTGCTAGCCAATCTCGTTTTCTTGTTTGCTTGCTTGTTTGTTTTGGTCAATGCAG")
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
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[3013603, 3014644]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 19159)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3014644:3014644+45], "CCTGTACCCTTTGGTGAGAATTTTTGTTTCAGTGTTAAAAGTTTG")
        self.assertEqual(alignment[0], "CCTGTACC---CTTTGGTGAGAATTTTTGTTTCAGTGTTAAAAGTTTG")
        self.assertEqual(alignment.sequences[1].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[1].seq), 170899992)
        self.assertEqual(alignment.sequences[1].seq[15870786:15870786+46], "CCTATACCTTTCTTTTATGAGAATTTTGTTTTAATCCTAAACTTTT")
        self.assertEqual(alignment[1], "CCTATACCTTTCTTTTATGAGAA-TTTTGTTTTAATCCTAAAC-TTTT")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 9085)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[2].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[2].seq), 173908612)
        self.assertEqual(alignment.sequences[2].seq[16389355:16389355+46], "CCTATACCTTTCTTTTATGAGAATTTTGTTTTAATCCTAAACTTTT")
        self.assertEqual(alignment[2], "CCTATACCTTTCTTTTATGAGAA-TTTTGTTTTAATCCTAAAC-TTTT")
        self.assertEqual(alignment.column_annotations["panTro2.chr6"], "99999999999999999999999-9999999999999999999-9999")
        self.assertEqual(alignment.sequences[2].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[2].annotations['leftCount'], 9106)
        self.assertEqual(alignment.sequences[2].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[3].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[3].seq), 133105)
        self.assertEqual(alignment.sequences[3].seq[6182:6182+46], "CCTATACCTTTCTTTCATGAGAATTTTGTTTGAATCCTAAACTTTT")
        self.assertEqual(alignment[3], "CCTATACCTTTCTTTCATGAGAA-TTTTGTTTGAATCCTAAAC-TTTT")
        self.assertEqual(alignment.sequences[3].annotations['leftStatus'], 'N')
        self.assertEqual(alignment.sequences[3].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[3].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[3].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[4].id, "loxAfr1.scaffold_75566")
        self.assertEqual(len(alignment.sequences[4].seq), 10574)
        self.assertEqual(alignment.sequences[4].seq[1167:1167+34], "TTTGGTTAGAATTATGCTTTAATTCAAAACTTCC")
        self.assertEqual(alignment[4], "------------TTTGGTTAGAA-TTATGCTTTAATTCAAAAC-TTCC")
        self.assertEqual(alignment.column_annotations["loxAfr1.scaffold_75566"], "------------99999699899-9999999999999869998-9997")
        self.assertEqual(alignment.sequences[4].annotations['leftStatus'], 'N')
        self.assertEqual(alignment.sequences[4].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[4].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[4].annotations['rightCount'], 0)
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
                numpy.array([[ 3014644,  3014652,  3014652,  3014653,  3014664,  3014665,
                               3014684,  3014685,  3014689],
                             [15870786, 15870794, 15870797, 15870798, 15870809, 15870809,
                              15870828, 15870828, 15870832],
                             [16389355, 16389363, 16389366, 16389367, 16389378, 16389378,
                              16389397, 16389397, 16389401],
                             [    6182,     6190,     6193,     6194,     6205,     6205,
                                  6224,     6224,     6228],
                             [    1167,     1167,     1167,     1167,     1178,     1178,
                                  1197,     1197,     1201]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 40840)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3014689:3014689+53], "GGGAGCATAAAACTCTAAATCTGCTAAATGTCTTGTCCCTTTGGAAAGAGTTG")
        self.assertEqual(alignment[0], "GGGAGCATAAAACTCTAAATCTGCTAAATGTCTTGTCCCT-TTGGAAAGAGTTG")
        self.assertEqual(alignment.sequences[1].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[1].seq), 170899992)
        self.assertEqual(alignment.sequences[1].seq[15870832:15870832+53], "GGGATCATAAACCATTTAATCTGTGAAATATCTAATCTTTTGGGAAATAGTGG")
        self.assertEqual(alignment[1], "GGGATCATAAACCATTTAATCTGTGAAATATCTAATCTTT-TGGGAAATAGTGG")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 401)
        self.assertEqual(alignment.sequences[2].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[2].seq), 173908612)
        self.assertEqual(alignment.sequences[2].seq[16389401:16389401+53], "GGGATCATAAACCATTTAATCTGTGAAATATCTAATCTTTTGGGAAATAGTGG")
        self.assertEqual(alignment[2], "GGGATCATAAACCATTTAATCTGTGAAATATCTAATCTTT-TGGGAAATAGTGG")
        self.assertEqual(alignment.column_annotations["panTro2.chr6"], "9999999999999999999999999999999999999999-9999999999999")
        self.assertEqual(alignment.sequences[2].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[2].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[2].annotations['rightCount'], 400)
        self.assertEqual(alignment.sequences[3].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[3].seq), 133105)
        self.assertEqual(alignment.sequences[3].seq[6228:6228+53], "GGGATCATAAGCCATTTAATCTGTGAAATGTGAAATCTTTTGGGAAACAGTGG")
        self.assertEqual(alignment[3], "GGGATCATAAGCCATTTAATCTGTGAAATGTGAAATCTTT-TGGGAAACAGTGG")
        self.assertEqual(alignment.sequences[3].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[3].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[3].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[3].annotations['rightCount'], 2)
        self.assertEqual(alignment.sequences[4].id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(alignment.sequences[4].seq), 359464)
        self.assertEqual(alignment.sequences[4].seq[184148:184148+52], "GGAAGCATAAACTTTTAATCTATGAAATATCAAATCACTTGGGCAATAGCTG")
        self.assertEqual(alignment[4], "GGAAGCATAAACT-TTTAATCTATGAAATATCAAATCACT-TGGGCAATAGCTG")
        self.assertEqual(alignment.column_annotations["otoGar1.scaffold_334.1-359464"], "7455455669566-99665699769895555689997599-9984787795599")
        self.assertEqual(alignment.sequences[4].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[4].annotations['leftCount'], 2931)
        self.assertEqual(alignment.sequences[4].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[4].annotations['rightCount'], 2)
        self.assertEqual(alignment.sequences[5].id, "loxAfr1.scaffold_75566")
        self.assertEqual(len(alignment.sequences[5].seq), 10574)
        self.assertEqual(alignment.sequences[5].seq[1201:1201+54], "GGGAGTATAAACCATTTAGTCTGCGAAATGCCAAATCTTCAGGGGAAAAAGCTG")
        self.assertEqual(alignment[5], "GGGAGTATAAACCATTTAGTCTGCGAAATGCCAAATCTTCAGGGGAAAAAGCTG")
        self.assertEqual(alignment.column_annotations["loxAfr1.scaffold_75566"], "899989799999979999999999999999797999999999999999999999")
        self.assertEqual(alignment.sequences[5].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[5].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[5].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[5].annotations['rightCount'], 2)
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
                numpy.array([[ 3014689,  3014702,  3014703,  3014729,  3014729,  3014742],
                             [15870832, 15870845, 15870846, 15870872, 15870872, 15870885],
                             [16389401, 16389414, 16389415, 16389441, 16389441, 16389454],
                             [    6228,     6241,     6242,     6268,     6268,     6281],
                             [  184148,   184161,   184161,   184187,   184187,   184200],
                             [    1201,     1214,     1215,     1241,     1242,     1255]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 411)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3014742:3014742+36], "AAGTTCCCTCCATAATTCCTTCCTCCCACCCCCACA")
        self.assertEqual(alignment[0], "AAGTTCCCTCCATAATTCCTTCCTCCCACCCCCACA")
        self.assertEqual(alignment.sequences[1].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[1].seq), 133105)
        self.assertEqual(alignment.sequences[1].seq[6283:6283+28], "AAATGTATGATCTCCCCATCCTGCCCTG")
        self.assertEqual(alignment[1], "AAATGTA-----TGATCTCCCCATCCTGCCCTG---")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 2)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 54)
        self.assertEqual(alignment.sequences[2].id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(alignment.sequences[2].seq), 359464)
        self.assertEqual(alignment.sequences[2].seq[184202:184202+31], "AGATTTCTGATGCCCTCACCCCCTCCGTGCA")
        self.assertEqual(alignment[2], "AGATTTC-----TGATGCCCTCACCCCCTCCGTGCA")
        self.assertEqual(alignment.column_annotations["otoGar1.scaffold_334.1-359464"], "9996999-----965974999999999999999999")
        self.assertEqual(alignment.sequences[2].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[2].annotations['leftCount'], 2)
        self.assertEqual(alignment.sequences[2].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[2].annotations['rightCount'], 24)
        self.assertEqual(alignment.sequences[3].id, "loxAfr1.scaffold_75566")
        self.assertEqual(len(alignment.sequences[3].seq), 10574)
        self.assertEqual(alignment.sequences[3].seq[1257:1257+27], "AGGCTTATGCCACCCCCCACCCCCACA")
        self.assertEqual(alignment[3], "AGGCTTA-----TG----CCACCCCCCACCCCCACA")
        self.assertEqual(alignment.column_annotations["loxAfr1.scaffold_75566"], "9999999-----99----997999999999999999")
        self.assertEqual(alignment.sequences[3].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[3].annotations['leftCount'], 2)
        self.assertEqual(alignment.sequences[3].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[3].annotations['rightCount'], 25)
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
                numpy.array([[3014742, 3014749, 3014754, 3014756, 3014760, 3014775, 3014778],
                             [   6283,    6290,    6290,    6292,    6296,    6311,    6311],
                             [ 184202,  184209,  184209,  184211,  184215,  184230,  184233],
                             [   1257,    1264,    1264,    1266,    1266,    1281,    1284]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3014778:3014778+17], "TCCCATGTCCACCCTGA")
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
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[3014778, 3014795]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, -12243)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3014795:3014795+47], "GTTTCAGGGGCAGCTCGCTGTTAGCAGCTAAGGCATGGTGTCTCTCA")
        self.assertEqual(alignment[0], "GTTTCAGGGGCAGCTCGCTG----------------TTAGCAG-CTAAGGCATGGTGTCTCTCA")
        self.assertEqual(alignment.sequences[1].id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(alignment.sequences[1].seq), 359464)
        self.assertEqual(alignment.sequences[1].seq[184257:184257+60], "TCGGGAACAAGTTGCAGTCATGGATATTTGGTTTAGATGGTTTAGGTGGGGTGTATTTCT")
        self.assertEqual(alignment[1], "---TCGGGAACAAGTTGCAGTCATGGATAT-TTGGTTTAGATGGTTTAGGTGGGGTGTATTTCT")
        self.assertEqual(alignment.column_annotations["otoGar1.scaffold_334.1-359464"], "---899999999999999999999999999-999999999999999999999999999999999")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 24)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[2].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[2].seq), 133105)
        self.assertEqual(alignment.sequences[2].seq[6365:6365+39], "GCCATGAATATTTTAGACATGCAGGTGTGGCGTGTTTCT")
        self.assertEqual(alignment[2], "-------------------GCCATGAATAT-----TTTAGAC-ATGCAGGTGTGGCGTGTTTCT")
        self.assertEqual(alignment.sequences[2].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[2].annotations['leftCount'], 54)
        self.assertEqual(alignment.sequences[2].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[3].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[3].seq), 170899992)
        self.assertEqual(alignment.sequences[3].seq[15871286:15871286+20], "aGCAGGTGTGGCATGTTTCT")
        self.assertEqual(alignment[3], "--------------------------------------------aGCAGGTGTGGCATGTTTCT")
        self.assertEqual(alignment.sequences[3].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[3].annotations['leftCount'], 401)
        self.assertEqual(alignment.sequences[3].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[3].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[4].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[4].seq), 173908612)
        self.assertEqual(alignment.sequences[4].seq[16389854:16389854+20], "aGCAGGTGTGGCATGTTTCT")
        self.assertEqual(alignment[4], "--------------------------------------------aGCAGGTGTGGCATGTTTCT")
        self.assertEqual(alignment.column_annotations["panTro2.chr6"], "--------------------------------------------99999999999999999999")
        self.assertEqual(alignment.sequences[4].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[4].annotations['leftCount'], 400)
        self.assertEqual(alignment.sequences[4].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[4].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[5].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[5].seq), 174210431)
        self.assertEqual(alignment.sequences[5].seq[16169492:16169492+20], "aGCAGGTGTGGCATGTTTCT")
        self.assertEqual(alignment[5], "--------------------------------------------aGCAGGTGTGGCATGTTTCT")
        self.assertEqual(alignment.sequences[5].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[5].annotations['leftCount'], 8044)
        self.assertEqual(alignment.sequences[5].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[5].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[6].id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(alignment.sequences[6].seq), 498454)
        self.assertEqual(alignment.sequences[6].seq[171521:171521+27], "TTAGAAATGTGGGTGTGGCGCATCGCT")
        self.assertEqual(alignment[6], "------------------------------------TTAGAAA-TGTGGGTGTGGCGCATCGCT")
        self.assertEqual(alignment.column_annotations["tupBel1.scaffold_114895.1-498454"], "------------------------------------9999999-99999989998899999699")
        self.assertEqual(alignment.sequences[6].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[6].annotations['leftCount'], 4145)
        self.assertEqual(alignment.sequences[6].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[6].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[7].id, "cavPor2.scaffold_216473")
        self.assertEqual(len(alignment.sequences[7].seq), 10026)
        self.assertEqual(alignment.sequences[7].seq[7816:7816+26], "TTAGAAATTTAGGTGTGGCATGGTTC")
        self.assertEqual(alignment[7], "------------------------------------TTAGAAA-TTTAGGTGTGGCATGGTTC-")
        self.assertEqual(alignment.column_annotations["cavPor2.scaffold_216473"], "------------------------------------4255831-1324566557465575854-")
        self.assertEqual(alignment.sequences[7].annotations['leftStatus'], 'N')
        self.assertEqual(alignment.sequences[7].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[7].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[7].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[8].id, "loxAfr1.scaffold_75566")
        self.assertEqual(len(alignment.sequences[8].seq), 10574)
        self.assertEqual(alignment.sequences[8].seq[1309:1309+62], "GTTTGAGAGCAAGTTTAGGACATCAATACGTTGTTTCAGAAATTGAGGTTTGCTGTATCTCt")
        self.assertEqual(alignment[8], "GTT-TGAGAGCAAGTTTAGGACATCAATACGTTGTTTCAGAAA-TTGAGGTTTGCTGTATCTCt")
        self.assertEqual(alignment.column_annotations["loxAfr1.scaffold_75566"], "999-999999999999999999999999999999999999999-99999999999999999999")
        self.assertEqual(alignment.sequences[8].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[8].annotations['leftCount'], 25)
        self.assertEqual(alignment.sequences[8].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[8].annotations['rightCount'], 0)
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
                numpy.array([[ 3014795,  3014798,  3014799,  3014814,  3014815,  3014815,
                               3014815,  3014815,  3014815,  3014821,  3014822,  3014822,
                               3014841,  3014842],
                             [  184257,   184257,   184258,   184273,   184274,   184284,
                                184284,   184288,   184289,   184295,   184296,   184297,
                                184316,   184317],
                             [    6365,     6365,     6365,     6365,     6366,     6376,
                                  6376,     6376,     6377,     6383,     6383,     6384,
                                  6403,     6404],
                             [15871286, 15871286, 15871286, 15871286, 15871286, 15871286,
                              15871286, 15871286, 15871286, 15871286, 15871286, 15871286,
                              15871305, 15871306],
                             [16389854, 16389854, 16389854, 16389854, 16389854, 16389854,
                              16389854, 16389854, 16389854, 16389854, 16389854, 16389854,
                              16389873, 16389874],
                             [16169492, 16169492, 16169492, 16169492, 16169492, 16169492,
                              16169492, 16169492, 16169492, 16169492, 16169492, 16169492,
                              16169511, 16169512],
                             [  171521,   171521,   171521,   171521,   171521,   171521,
                                171521,   171521,   171521,   171527,   171528,   171528,
                                171547,   171548],
                             [    7816,     7816,     7816,     7816,     7816,     7816,
                                  7816,     7816,     7816,     7822,     7823,     7823,
                                  7842,     7842],
                             [    1309,     1312,     1312,     1327,     1328,     1338,
                                  1339,     1343,     1344,     1350,     1351,     1351,
                                  1370,     1371]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 320596)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3014842:3014842+186], "CTTGGGATGCTTTATAGTGGAAATGGAAAGCAATTTATTTAGATCTTAAATCATTTTGAAGGTTAATAAAATGACCATATTAATATTCCCATGAACAAAGCCTTCATTTTTAAAATATTGCATCCTATAATACACATAAATCTTGTTCTCGtttttatttttttatttatttatttttttttcttt")
        self.assertEqual(alignment[0], "C--TTGGGA---------TGCTTTATAGTGGAAATGGAAAGCA----A-TTTATTTAGATCTTAAATCATTTT-GAAGGTTAATAAAATGACCATATTAATATTCCCATGAACAAAGCCTTCATTT----TTAAAATATTGCATCCTATAATACACATAA-ATCTTGT-----TCTCGtttttatttttt----tatt-tat-----------------------ttattttttttt------------cttt")
        self.assertEqual(alignment.sequences[1].id, "loxAfr1.scaffold_75566")
        self.assertEqual(len(alignment.sequences[1].seq), 10574)
        self.assertEqual(alignment.sequences[1].seq[1371:1371+245], "tatttttgatttttttttttttttaaaGCTGAAGTGAAAGTCACACTGTTTCTTTTGGCCATAAATCCCTTTGAAGTATAATGAAATTTTCACTTTAATATTCCTGTGAACAATACCCTCATTTAAAAAAAAAAATTGGGTTTTATACCACACACAGCATCTTTTCAAAATCTCATTCCCTTTTTAAGCAGGTGTCTTAGCAAGTTCACTCTCTCTCTCTTTTTTCATTGACTTTGGCCAAATCT")
        self.assertEqual(alignment[1], "tatttttgatttttttttttttttaaaGCTGAAGTGAAAGTCACACTG-TTTCTTTTGGCCATAAATCCCTTT-GAAGTATAATGAAATTTTCACTTTAATATTCCTGTGAACAATACCCTCATTT-AAAAAAAAAAATTGGGTTTTATACCACACACAGCATCTTTTCAAAATCTCATTCCC-TTTTTAAGCAGGTG-TCT---TAGCAAGTTCACTCTCTCTCTCTTTTTTCATTGACTTTGGCCAAATCT")
        self.assertEqual(alignment.column_annotations["loxAfr1.scaffold_75566"], "999999999999999999999999999999999999999999999999-999999999999999999999999-9999999999999999999999999999999999999999999999999999-99999999999999999999999999999999999999999999999999999999-99999999999999-999---999999999999999999999999999999999999999999999999")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[2].id, "cavPor2.scaffold_216473")
        self.assertEqual(len(alignment.sequences[2].seq), 10026)
        self.assertEqual(alignment.sequences[2].seq[7842:7842+198], "TTTGGATATTCACAGTTGGAATGAAAGGCACCCTGTTTCTTTAGGACTTAAATCCTTTTGAAGTATAATAAAATGATCACATTAATATTCCTACTAACAATGCCCTCATTTAAAAAAGTTTGGATTTTGTACCACATACAGCATATTTCCAACATCTCATTTCCTTCTTCTGGCAGGTTTAttagtctctctttctct")
        self.assertEqual(alignment[2], "---TTTGGA---------TA-TTCACAGTTGGAATGAAAGGCACCCTG-TTTCTTTAGGACTTAAATCCTTTT-GAAGTATAATAAAATGATCACATTAATATTCCTACTAACAATGCCCTCATTT----AAAAAAGTTTGGATTTTGTACCACATACAGCATATTTCCAACATCTCATTTCCTTCTTCTGGCAGGTT-TAt-----------------------tagtctctcttt------------ctct")
        self.assertEqual(alignment.column_annotations["cavPor2.scaffold_216473"], "---610137---------77-200195531236266876425368858-778987956887868956898956-8988778987788768588885664786777656586678636299978766----89979789936989956687867689995888978886997659897789899998996778899997-989-----------------------998799999777------------9899")
        self.assertEqual(alignment.sequences[2].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[2].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[3].id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(alignment.sequences[3].seq), 498454)
        self.assertEqual(alignment.sequences[3].seq[171548:171548+215], "CCTTTTGGATTTTTTACAGCTGAAATAAAAAGCACACTGTGTCTGTAGGTCTTAAATCCTTTTGAAGTATAATAAAATGATCACTTTACTATTCCTGTGAACAATGCCCTCATTTTACAAAATCTGGGATTTTATACCACGTACAGCATATTTCCAAAATCTCATTTCCTTTTTTGGCAGGTTTATTAGCAAGTTCCCTTTTTGTCTATTGACTT")
        self.assertEqual(alignment[3], "CCTTTTGGA---------TTTTTTACAGCTGAAATAAAAAGCACACTG-TGTCTGTAGGTCTTAAATCCTTTT-GAAGTATAATAAAATGATCACTTTACTATTCCTGTGAACAATGCCCTCATTT---TACAAAATCTGGGATTTTATACCACGTACAGCATATTTCCAAAATCTCATTTCC-TTTTTTGGCAGGTT-TAT---TAGCAAGTTCCCT-------TTTTGTCTATTG------------ACTT")
        self.assertEqual(alignment.column_annotations["tupBel1.scaffold_114895.1-498454"], "989999999---------899999998999999999999999999999-999999999999999999999999-9999998999999999999999999999999999899999999999999999---998899999999973999999998999999999999999979999999999999-99999997699849-999---9999999999999-------999999999999------------9999")
        self.assertEqual(alignment.sequences[3].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[3].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[3].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[3].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[4].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[4].seq), 174210431)
        self.assertEqual(alignment.sequences[4].seq[16169512:16169512+219], "CCTTTTGGACTTTTTATAGCTGAAATGGAAGGCACACTGTTTCTTTAGGTCTTAAATACTTTTGAAGTATAATAAAACGATCACTTTAATACTCCTGTGAGCAATGTCCTCTTTTAAAAAAATCCTGGATTTTATACCACATGCAGCATATTTCTGAAATCTCATTTCCTTTCTTGGCAGGTTTATTTATAGCAAGTTATCTCTCTCTCTTATTTACTT")
        self.assertEqual(alignment[4], "CCTTTTGGA---------CTTTTTATAGCTGAAATGGAAGGCACACTG-TTTCTTTAGGTCTTAAATACTTTT-GAAGTATAATAAAACGATCACTTTAATACTCCTGTGAGCAATGTCCTCTTTT---AAAAAAATCCTGGATTTTATACCACATGCAGCATATTTCTGAAATCTCATTTCC-TTTCTTGGCAGGTT-TATTTATAGCAAGTTATCTC------TCTCTCTTATTT------------ACTT")
        self.assertEqual(alignment.sequences[4].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[4].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[4].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[4].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[5].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[5].seq), 173908612)
        self.assertEqual(alignment.sequences[5].seq[16389874:16389874+168], "CCTTTTGGACTTTTTATAGCTGAAATGGAAGGCACAATATTTCTTTAGGTCTTAAATACTTTTGAAGTATAAAAAAATGATCACTTTAATACTCCTGTGAACAATGCCCTCTTTTTAAAAAAATCCTGGATTTTATACCACATGCAGCATATTTCTGAAATCCCATTT")
        self.assertEqual(alignment[5], "CCTTTTGGA---------CTTTTTATAGCTGAAATGGAAGGCACAATA-TTTCTTTAGGTCTTAAATACTTTT-GAAGTATAAAAAAATGATCACTTTAATACTCCTGTGAACAATGCCCTCTTTTT--AAAAAAATCCTGGATTTTATACCACATGCAGCATATTTCTGAAATCCCATTT------------------------------------------------------------------------")
        self.assertEqual(alignment.column_annotations["panTro2.chr6"], "999999999---------999999999999999999999999999999-999999999999999999999999-99999999999999999999999999999999999999999999999999999--9999999999999999999999999999999999999999999999999999------------------------------------------------------------------------")
        self.assertEqual(alignment.sequences[5].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[5].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[5].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[5].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[6].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[6].seq), 170899992)
        self.assertEqual(alignment.sequences[6].seq[15871306:15871306+169], "CCTTTTGGACTTTTTATAGCTGAAATGGAAGGCACAATATTTCTTTAGGTCTTAAATACTTTTGAAGTATAAAAAAATGATCACTTTAATACTCCTGTGAACAATGCCCTCTTTTTAAAAAAAATCCTGGATTTTATACCACATGCAGCATATTTCTGAAATCTCATTT")
        self.assertEqual(alignment[6], "CCTTTTGGA---------CTTTTTATAGCTGAAATGGAAGGCACAATA-TTTCTTTAGGTCTTAAATACTTTT-GAAGTATAAAAAAATGATCACTTTAATACTCCTGTGAACAATGCCCTCTTTTT-AAAAAAAATCCTGGATTTTATACCACATGCAGCATATTTCTGAAATCTCATTT------------------------------------------------------------------------")
        self.assertEqual(alignment.sequences[6].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[6].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[6].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[6].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[7].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[7].seq), 133105)
        self.assertEqual(alignment.sequences[7].seq[6404:6404+223], "TCTTTTGTATTTTTTTATAGCTGAAAAGGAAGGCACACTGTCTCTTTAGGTCTTAAATACGTTTGAAGCATAATAAAATGATCACTTTCATACTCCTGTGAATAATGCCCTCCTTTTAAAAAGAAATCTTGGATTTTATATCACACGCAGCATATTTCTGAAATCTCATTTCCTTTCCTGGCAGGTTTATTAATAGCAAGTTCTCTCTCTTATTTATTTATTT")
        self.assertEqual(alignment[7], "TCTTTTGTA--------TTTTTTTATAGCTGAAAAGGAAGGCACACTG-TCTCTTTAGGTCTTAAATACGTTT-GAAGCATAATAAAATGATCACTTTCATACTCCTGTGAATAATGCCCTCCTTTTAAAAAGAAATCTTGGATTTTATATCACACGCAGCATATTTCTGAAATCTCATTTCC-TTTCCTGGCAGGTT-TATTAATAGCAAGTTCTCTC------TCTTATTTATTT------------ATTT")
        self.assertEqual(alignment.sequences[7].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[7].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[7].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[7].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[8].id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(alignment.sequences[8].seq), 359464)
        self.assertEqual(alignment.sequences[8].seq[184317:184317+222], "TATTTTGGGCTTTTTTATGGGTGAAATGAAAGGCACACTGTTTTTTTAGGACTTAAATCCTTTTGAAGTATAATAAAACGATTACTTTACTGTTCCTGTGAACAATGCCCTCATTTAAAAAAATCTTGGATTTTATACCAGATCCAGCATATTTCCAAAATCTCATTTCCTTTTTGGCGGGTTTGTTAGCAAGTTCTCTCCGGCTTTCTTCTTCTTTTTCTC")
        self.assertEqual(alignment[8], "TATTTTGGG--------CTTTTTTATGGGTGAAATGAAAGGCACACTG-TTTTTTTAGGACTTAAATCCTTTT-GAAGTATAATAAAACGATTACTTTACTGTTCCTGTGAACAATGCCCTCATTT---AAAAAAATCTTGGATTTTATACCAGATCCAGCATATTTCCAAAATCTCATTT-C-CTTTTTGGCGGGTT-TGT---TAGCAAGTTCTCTCCGGCTTTCTTCTTCTTTT------------TCTC")
        self.assertEqual(alignment.column_annotations["otoGar1.scaffold_334.1-359464"], "999999999--------9999999999999999999999997999999-999999999999999999999999-9999999999999999999999999999999999999999999999999999---9999999999999999999999999999999999999999999999999999-9-99999999999999-999---99999999999999999999999999999999------------9999")
        self.assertEqual(alignment.sequences[8].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[8].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[8].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[8].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[9].id, "ornAna1.chr2")
        self.assertEqual(len(alignment.sequences[9].seq), 54797317)
        self.assertEqual(alignment.sequences[9].seq[14750994:14750994+201], "CCTTGGAAATGTTACTGCAAAAGTGAAAGGCTCATTGTTTTCTTTCTATCTTAAGCCTTTTTAAAAGTGTATTGTAATTATGATCTTCATTTCCCCGCTAACAAAGCCCTCATTAAAATATCTTTGTATTTTAACGACATACAGCCCATTTCCAAAATCTTAGTTCCTTTTTGGGTGCATTGTATTAACCGTTCTCTCTTT")
        self.assertEqual(alignment[9], "--CCTTGGA---------AATGTTACTGCAAAAGTGAAAGGCTCATTGTTTTCTTTCTATCTTAAGCCTTTTTAAAAGTGTATTGTAATTATGATCTTCATTTCCCCGCTAACAAAGCCCTCATTA----AAATATCTTTGTATTTTA-ACGACATACAGCCCATTTCCAAAATCTTAGTTCC-TTTTTGGGTGCATTGTAT---TAAC----------------CGTTCTCTCTTT----------------")
        self.assertEqual(alignment.sequences[9].annotations['leftStatus'], 'N')
        self.assertEqual(alignment.sequences[9].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[9].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[9].annotations['rightCount'], 5690)
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
                numpy.array([[ 3014842,  3014843,  3014843,  3014843,  3014849,  3014849,
                               3014849,  3014851,  3014852,  3014874,  3014874,  3014875,
                               3014875,  3014899,  3014899,  3014951,  3014951,  3014951,
                               3014951,  3014951,  3014969,  3014970,  3014981,  3014981,
                               3014988,  3014988,  3014996,  3014997,  3014998,  3014999,
                               3015005,  3015005,  3015009,  3015009,  3015012,  3015012,
                               3015012,  3015012,  3015012,  3015012,  3015024,  3015024,
                               3015028],
                             [    1371,     1372,     1373,     1374,     1380,     1388,
                                  1389,     1391,     1392,     1414,     1418,     1419,
                                  1419,     1443,     1443,     1495,     1495,     1496,
                                  1497,     1498,     1516,     1517,     1528,     1529,
                                  1536,     1541,     1549,     1550,     1551,     1551,
                                  1557,     1561,     1565,     1565,     1568,     1568,
                                  1572,     1581,     1582,     1588,     1600,     1612,
                                  1616],
                             [    7842,     7842,     7842,     7842,     7848,     7848,
                                  7848,     7850,     7850,     7872,     7876,     7877,
                                  7877,     7901,     7901,     7953,     7953,     7953,
                                  7953,     7953,     7971,     7972,     7983,     7984,
                                  7991,     7996,     8004,     8005,     8006,     8007,
                                  8013,     8017,     8021,     8021,     8024,     8024,
                                  8024,     8024,     8024,     8024,     8036,     8036,
                                  8040],
                             [  171548,   171549,   171550,   171551,   171557,   171557,
                                171557,   171559,   171560,   171582,   171586,   171587,
                                171587,   171611,   171611,   171663,   171663,   171663,
                                171663,   171664,   171682,   171683,   171694,   171695,
                                171702,   171707,   171715,   171716,   171717,   171717,
                                171723,   171727,   171731,   171731,   171734,   171734,
                                171738,   171747,   171747,   171747,   171759,   171759,
                                171763],
                             [16169512, 16169513, 16169514, 16169515, 16169521, 16169521,
                              16169521, 16169523, 16169524, 16169546, 16169550, 16169551,
                              16169551, 16169575, 16169575, 16169627, 16169627, 16169627,
                              16169627, 16169628, 16169646, 16169647, 16169658, 16169659,
                              16169666, 16169671, 16169679, 16169680, 16169681, 16169681,
                              16169687, 16169691, 16169695, 16169695, 16169698, 16169701,
                              16169705, 16169714, 16169715, 16169715, 16169727, 16169727,
                              16169731],
                             [16389874, 16389875, 16389876, 16389877, 16389883, 16389883,
                              16389883, 16389885, 16389886, 16389908, 16389912, 16389913,
                              16389913, 16389937, 16389937, 16389989, 16389990, 16389990,
                              16389990, 16389991, 16390009, 16390010, 16390021, 16390022,
                              16390029, 16390034, 16390042, 16390042, 16390042, 16390042,
                              16390042, 16390042, 16390042, 16390042, 16390042, 16390042,
                              16390042, 16390042, 16390042, 16390042, 16390042, 16390042,
                              16390042],
                             [15871306, 15871307, 15871308, 15871309, 15871315, 15871315,
                              15871315, 15871317, 15871318, 15871340, 15871344, 15871345,
                              15871345, 15871369, 15871369, 15871421, 15871422, 15871422,
                              15871423, 15871424, 15871442, 15871443, 15871454, 15871455,
                              15871462, 15871467, 15871475, 15871475, 15871475, 15871475,
                              15871475, 15871475, 15871475, 15871475, 15871475, 15871475,
                              15871475, 15871475, 15871475, 15871475, 15871475, 15871475,
                              15871475],
                             [    6404,     6405,     6406,     6407,     6413,     6413,
                                  6414,     6416,     6417,     6439,     6443,     6444,
                                  6444,     6468,     6468,     6520,     6521,     6522,
                                  6523,     6524,     6542,     6543,     6554,     6555,
                                  6562,     6567,     6575,     6576,     6577,     6577,
                                  6583,     6587,     6591,     6591,     6594,     6597,
                                  6601,     6610,     6611,     6611,     6623,     6623,
                                  6627],
                             [  184317,   184318,   184319,   184320,   184326,   184326,
                                184327,   184329,   184330,   184352,   184356,   184357,
                                184357,   184381,   184381,   184433,   184433,   184433,
                                184433,   184434,   184452,   184453,   184464,   184465,
                                184472,   184477,   184485,   184485,   184486,   184486,
                                184492,   184496,   184500,   184500,   184503,   184503,
                                184507,   184516,   184517,   184523,   184535,   184535,
                                184539],
                             [14750994, 14750994, 14750994, 14750995, 14751001, 14751001,
                              14751001, 14751003, 14751004, 14751026, 14751030, 14751031,
                              14751032, 14751056, 14751057, 14751109, 14751109, 14751109,
                              14751109, 14751109, 14751127, 14751127, 14751138, 14751139,
                              14751146, 14751151, 14751159, 14751160, 14751161, 14751161,
                              14751167, 14751171, 14751175, 14751176, 14751179, 14751179,
                              14751183, 14751183, 14751183, 14751183, 14751195, 14751195,
                              14751195]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, -36127)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3015028:3015028+58], "ccattttttattaggtatttagctcatttacatttccaatgctataccaaaagtcccc")
        self.assertEqual(alignment[0], "ccatttt----------ttattaggtatttagctcatttacatttccaatgctatac----caaaagtcccc")
        self.assertEqual(alignment.sequences[1].id, "loxAfr1.scaffold_75566")
        self.assertEqual(len(alignment.sequences[1].seq), 10574)
        self.assertEqual(alignment.sequences[1].seq[1616:1616+33], "TCACTTCTACCTCCATTCCTTACAAAAGTTCTC")
        self.assertEqual(alignment[1], "TCACTTCTA---------------------------------------CCTCCATTCCTTACAAAAGTTCTC")
        self.assertEqual(alignment.column_annotations["loxAfr1.scaffold_75566"], "999999999---------------------------------------999999999999999999999999")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'N')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[2].id, "cavPor2.scaffold_216473")
        self.assertEqual(len(alignment.sequences[2].seq), 10026)
        self.assertEqual(alignment.sequences[2].seq[8040:8040+8], "ctgtgtct")
        self.assertEqual(alignment[2], "ctgtgtc----------t------------------------------------------------------")
        self.assertEqual(alignment.column_annotations["cavPor2.scaffold_216473"], "6788989----------9------------------------------------------------------")
        self.assertEqual(alignment.sequences[2].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[2].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[2].annotations['rightCount'], 1372)
        self.assertEqual(alignment.sequences[3].id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(alignment.sequences[3].seq), 498454)
        self.assertEqual(alignment.sequences[3].seq[171763:171763+1], "T")
        self.assertEqual(alignment[3], "T-----------------------------------------------------------------------")
        self.assertEqual(alignment.column_annotations["tupBel1.scaffold_114895.1-498454"], "9-----------------------------------------------------------------------")
        self.assertEqual(alignment.sequences[3].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[3].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[3].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[3].annotations['rightCount'], 2374)
        self.assertEqual(alignment.sequences[4].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[4].seq), 174210431)
        self.assertEqual(alignment.sequences[4].seq[16169731:16169731+12], "TTTTGCCTTTTC")
        self.assertEqual(alignment[4], "TTTTGCCTTTTC------------------------------------------------------------")
        self.assertEqual(alignment.sequences[4].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[4].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[4].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[4].annotations['rightCount'], 75)
        self.assertEqual(alignment.sequences[5].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[5].seq), 133105)
        self.assertEqual(alignment.sequences[5].seq[6627:6627+22], "TGCTTTTTCTCAAATCTCCACT")
        self.assertEqual(alignment[5], "TGCTTTTTCTCAAATCTCCACT--------------------------------------------------")
        self.assertEqual(alignment.sequences[5].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[5].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[5].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[5].annotations['rightCount'], 3479)
        self.assertEqual(alignment.sequences[6].id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(alignment.sequences[6].seq), 359464)
        self.assertEqual(alignment.sequences[6].seq[184539:184539+36], "TTTTTTCTATTGACTTTTGCCAAATATTCACTTCCA")
        self.assertEqual(alignment[6], "TTTTTTCTATT------------GACTTTTGCCAAATATTCACTTCCA------------------------")
        self.assertEqual(alignment.column_annotations["otoGar1.scaffold_334.1-359464"], "99999999999------------9999999999999999999999999------------------------")
        self.assertEqual(alignment.sequences[6].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[6].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[6].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[6].annotations['rightCount'], 137)
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
                numpy.array([[ 3015028,  3015029,  3015035,  3015035,  3015035,  3015035,
                               3015035,  3015036,  3015040,  3015041,  3015066,  3015075,
                               3015075,  3015086],
                             [    1616,     1617,     1623,     1625,     1625,     1625,
                                  1625,     1625,     1625,     1625,     1625,     1634,
                                  1638,     1649],
                             [    8040,     8041,     8047,     8047,     8047,     8047,
                                  8047,     8048,     8048,     8048,     8048,     8048,
                                  8048,     8048],
                             [  171763,   171764,   171764,   171764,   171764,   171764,
                                171764,   171764,   171764,   171764,   171764,   171764,
                                171764,   171764],
                             [16169731, 16169732, 16169738, 16169740, 16169742, 16169743,
                              16169743, 16169743, 16169743, 16169743, 16169743, 16169743,
                              16169743, 16169743],
                             [    6627,     6628,     6634,     6636,     6638,     6639,
                                  6644,     6645,     6649,     6649,     6649,     6649,
                                  6649,     6649],
                             [  184539,   184540,   184546,   184548,   184550,   184550,
                                184550,   184550,   184550,   184550,   184575,   184575,
                                184575,   184575]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3015086:3015086+2572], "catacccacccacccccactcccctacccgcccactccccctttttggccctggcgttcccctgttctggggcatataaagtttgtgtgtccaatgggcctctctttccagtgatggccgactaggccatcttttgatacatatgcagctagagtcaagagctccggggtactggttagttcataatgttgatccacctatagggttgcagatccctttagctccttgggtactttctctagctcccccattgggagccctgtgatccatccattagctgactgtgggcatccacttctgtgtttgctaggccccggcatagtctcacaagagacagctacatctgggtcctttcgataagatcttgctagtgtatgcaatggtgtcagcgtttggatgctgattatggggtagatccctggataaggcagtctctacatggtccatcctttcatctcagctccaaactttgtctctgtaactccttccaagggtgttttgttcccacttctaaggaggggcatagtgtccacacttcagtcttcatttttcttgagtttcatgtgtttaggaaattgtatcttatatcttgggtatcctaggttttgggctaatatccacttatcagtgagtacatattgtgtgagttcctttgtgaatgtgttacctcactcaggatgatgccctccaggtccatccatttggctaggaatttcataaattaattctttttaatagctgagtagtactccattgtgtagatgtaccacattttctgtatccattcctctgttgaggggcatctgggttctttccagcttctggctattataaataaggctgctatgaacatagtggagcatgtgtccttcttaccagttggggcatcttttggatatatgcccaggagaggtattgctggatcctccggtagtactatgtccaattttctgaggaaccgccagacggatttccagagtggttgtacaagcctgcaatcccaccaacaatggaggagtgttcctctttctccacatcctcgccagcatctgctgtcacctgaatttttgatcttagccattctgactggtgtgaggtggaatctcagggttgttttgatttgcatttctctgatgattaaggatgttgaacatgttttcaggtgcttctctgccattcggtattcctcaggtgagaattctttgttcagttctgagccccattttttaatggggttatttgattttctgaagtccaccttcttgagttctttatatatgttggatattagtcccctatctgatttaggataggtaaagatcctttcccaatctgttggtggtctctttgtgttattgacggtgtcttttgccttgcagaaactttggagtttcattaggtcccatttgtcaattctcgatcttacagcacaagccattgctgttctgttcaggaatttttcccctgtgcccatatcttcaaggcttttccccactttctcctctataagtttcagtgtctctggttttatgtggagttctttgatccatttagatttgaccttagtacaaggagataagtatggatcgattcgcattcttctacatgataacaaccagttgtgccagcaccaattgttgaaaatgctgtctttcttccactggatggttttagctcccttgtcgaagatcaagtgaccataggtgtgtgggttcatttctgggtcttcaattctattccattggtctacttgtctgtctctataccagtaccatgcagtttttaccacaattgctctgtagtaaagctttaggtcaggcatggtgattccaccagaggttcttttatccttgagaagagtttttgctatcctaggttttttgttattccagatgaatttgcaaattgctccttctaattcgttgaagaattgagttggaattgtgatggggattgcattgaatctgtagattgcttttggcaagatagccatttttacaatgttgatcctgccaatccatgagcatgggagagctttccatcttctgagatcttctttaatttctttcttcagagacttgaagtttttatcatacagatctttcacttccttagttagagtcacgccgagatattttatattatttgtgactattgagaagggtgttgtttccctaatttctttctcagcctgtttattctttgtgtagagaaaggccattgacttgtttgagttaattttatatccagctacttcaccgaagctgtttatcaggtttaggagttctctggtggaatttttagggtcacttatatatactatcatatcatctgcaaaaagtgatattttgacttcctcctttccaatttgtatccccttgatctccttttgttgtcgaattgctctggctaatacttcaagtactatgttgaaaaggtagggagaaagtgggcagccttgtctagtccctgattttagtgagattgcttccagcttctctccatttactttgatgttggctactggtttgctgtagattgcttttatcatgtttaggtatgggTGTTCTCG")
        self.assertEqual(alignment[0], "catacccacccacccccactcccctacccgcccactccccctttttggccctggcgttcccctgttctggggcatataaagtttgtgtgtccaatgggcctctctttccagtgatggccgactaggccatcttttgatacatatgcagctagagtcaagagctccggggtactggttagttcataatgttgatccacctatagggttgcagatccctttagctccttgggtactttctctagctcccccattgggagccctgtgatccatccattagctgactgtgggcatccacttctgtgtttgctaggccccggcatagtctcacaagagacagctacatctgggtcctttcgataagatcttgctagtgtatgcaatggtgtcagcgtttggatgctgattatggggtagatccctggataaggcagtctctacatggtccatcctttcatctcagctccaaactttgtctctgtaactccttccaagggtgttttgttcccacttctaaggaggggcatagtgtccacacttcagtcttcatttttcttgagtttcatgtgtttaggaaattgtatcttatatcttgggtatcctaggttttgggctaatatccacttatcagtgagtacatattgtgtgagttcctttgtgaatgtgttacctcactcaggatgatgccctccaggtccatccatttggctaggaatttcataaattaattctttttaatagctgagtagtactccattgtgtagatgtaccacattttctgtatccattcctctgttgaggggcatctgggttctttccagcttctggctattataaataaggctgctatgaacatagtggagcatgtgtccttcttaccagttggggcatcttttggatatatgcccaggagaggtattgctggatcctccggtagtactatgtccaattttctgaggaaccgccagacggatttccagagtggttgtacaagcctgcaatcccaccaacaatggaggagtgttcctctttctccacatcctcgccagcatctgctgtcacctgaatttttgatcttagccattctgactggtgtgaggtggaatctcagggttgttttgatttgcatttctctgatgattaaggatgttgaacatgttttcaggtgcttctctgccattcggtattcctcaggtgagaattctttgttcagttctgagccccattttttaatggggttatttgattttctgaagtccaccttcttgagttctttatatatgttggatattagtcccctatctgatttaggataggtaaagatcctttcccaatctgttggtggtctctttgtgttattgacggtgtcttttgccttgcagaaactttggagtttcattaggtcccatttgtcaattctcgatcttacagcacaagccattgctgttctgttcaggaatttttcccctgtgcccatatcttcaaggcttttccccactttctcctctataagtttcagtgtctctggttttatgtggagttctttgatccatttagatttgaccttagtacaaggagataagtatggatcgattcgcattcttctacatgataacaaccagttgtgccagcaccaattgttgaaaatgctgtctttcttccactggatggttttagctcccttgtcgaagatcaagtgaccataggtgtgtgggttcatttctgggtcttcaattctattccattggtctacttgtctgtctctataccagtaccatgcagtttttaccacaattgctctgtagtaaagctttaggtcaggcatggtgattccaccagaggttcttttatccttgagaagagtttttgctatcctaggttttttgttattccagatgaatttgcaaattgctccttctaattcgttgaagaattgagttggaattgtgatggggattgcattgaatctgtagattgcttttggcaagatagccatttttacaatgttgatcctgccaatccatgagcatgggagagctttccatcttctgagatcttctttaatttctttcttcagagacttgaagtttttatcatacagatctttcacttccttagttagagtcacgccgagatattttatattatttgtgactattgagaagggtgttgtttccctaatttctttctcagcctgtttattctttgtgtagagaaaggccattgacttgtttgagttaattttatatccagctacttcaccgaagctgtttatcaggtttaggagttctctggtggaatttttagggtcacttatatatactatcatatcatctgcaaaaagtgatattttgacttcctcctttccaatttgtatccccttgatctccttttgttgtcgaattgctctggctaatacttcaagtactatgttgaaaaggtagggagaaagtgggcagccttgtctagtccctgattttagtgagattgcttccagcttctctccatttactttgatgttggctactggtttgctgtagattgcttttatcatgtttaggtatgggTGTTCTCG")
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
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[3015086, 3017658]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 12170)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3017658:3017658+85], "TTTTTATTTGCAGGTTTCTTTACAGTTCTCTTTCATTCTTCTCCTCTTTTCTTCTGTTGACCTTTATCAGATTTCTGCTTTAACC")
        self.assertEqual(alignment[0], "TTTTTATTTGCAGGTTTCTTTAC----AGTTCTCTTTCATTCTTCTCCTCTTTTCTTCTGTTGACCTTTATCAGATTTCTGCTTTAACC")
        self.assertEqual(alignment.sequences[1].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[1].seq), 170899992)
        self.assertEqual(alignment.sequences[1].seq[15871475:15871475+83], "CCTTTCTTGGCAGGGTTATTTATAGCAAGTTATCTCTCTCTCTTATTTATTTTTTTGCCTTTTCCCAAATCTCCACTTCCACC")
        self.assertEqual(alignment[1], "CCTTTCTTGGCAGGGTTATTTATAGCAAGTTATCTCTCTCTCTTA------TTTATTTTTTTGCCTTTTCCCAAATCTCCACTTCCACC")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 53)
        self.assertEqual(alignment.sequences[2].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[2].seq), 173908612)
        self.assertEqual(alignment.sequences[2].seq[16390042:16390042+83], "CCTTTCTTGGCAGGGTTATTTATAGCAAGTTATCTCTCTCTCTTATTTATTTTTTTGCCTTTTCCCAAATCTCCACTTCCACC")
        self.assertEqual(alignment[2], "CCTTTCTTGGCAGGGTTATTTATAGCAAGTTATCTCTCTCTCTTA------TTTATTTTTTTGCCTTTTCCCAAATCTCCACTTCCACC")
        self.assertEqual(alignment.column_annotations["panTro2.chr6"], "999999999999999999999999999999999999999999999------99999999999999999999999999999999999999")
        self.assertEqual(alignment.sequences[2].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[2].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[2].annotations['rightCount'], 53)
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
                numpy.array([[ 3017658,  3017681,  3017681,  3017699,  3017705,  3017743],
                             [15871475, 15871498, 15871502, 15871520, 15871520, 15871558],
                             [16390042, 16390065, 16390069, 16390087, 16390087, 16390125]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3017743:3017743+418], "ACCACAGACCTTCTGTTTAGTCCAAAGGACGCAAATTATGTATCCACTTtagtaggaggctgacccgcagcctacatgaaccaggtatttctggaaggcaggctggggttgaaagagaaattagatggtgagaaaagaataatgaggccaagacaaatttttctcttatcaaggcccaagagagtttactaagagactatgcttaaaagggggaaggcccatcccccccccccctcgcgccagtctatccttggtgctttgtcaccatgccatcagcacttggtcggcaggtagcagaatctcagggcagttgacacttcaaaagaaaccagccaagtcagaaagctgcactgcaggagacctgcactcagtggtgacaaggtctgtaccagcctgcttcaggctgggggaggctaca")
        self.assertEqual(alignment[0], "ACCACAGACCTTCTGTTTAGTCCAAAGGACGCAAATTATGTATCCACTTtagtaggaggctgacccgcagcctacatgaaccaggtatttctggaaggcaggctggggttgaaagagaaattagatggtgagaaaagaataatgaggccaagacaaatttttctcttatcaaggcccaagagagtttactaagagactatgcttaaaagggggaaggcccatcccccccccccctcgcgccagtctatccttggtgctttgtcaccatgccatcagcacttggtcggcaggtagcagaatctcagggcagttgacacttcaaaagaaaccagccaagtcagaaagctgcactgcaggagacctgcactcagtggtgacaaggtctgtaccagcctgcttcaggctgggggaggctaca")
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
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[3017743, 3018161]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 22499)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3018161:3018161+69], "ATCCACAAAAGAGACAAAGAAGAAAACCAAAAGAAAAGATTGTAGCTTAAAACAATTCCATTTTATTGA")
        self.assertEqual(alignment[0], "ATCCACAAAAGAGAC-----AAAGAAGAAAACCAAAAGAAAAGATTGTAGCTTAAAACAATTCCATTTTATTGA")
        self.assertEqual(alignment.sequences[1].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[1].seq), 173908612)
        self.assertEqual(alignment.sequences[1].seq[16390178:16390178+65], "ACCTACAAAGGAAAACAATTAACCATAAGAAAAGTTTGTACCATAAAACAATTTTATTTTATTGA")
        self.assertEqual(alignment[1], "ACCTACAAAGG---------AAAACAATTAACCATAAGAAAAGTTTGTACCATAAAACAATTTTATTTTATTGA")
        self.assertEqual(alignment.column_annotations["panTro2.chr6"], "99999999999---------999999999999999999999999999999999999999999999999999999")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 53)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[2].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[2].seq), 170899992)
        self.assertEqual(alignment.sequences[2].seq[15871611:15871611+65], "ACCTACAAAGGAAAACAATTAACCACAAGAAAAGTTTGTACCATAAAACAATTTTATTTTATTGA")
        self.assertEqual(alignment[2], "ACCTACAAAGG---------AAAACAATTAACCACAAGAAAAGTTTGTACCATAAAACAATTTTATTTTATTGA")
        self.assertEqual(alignment.sequences[2].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[2].annotations['leftCount'], 53)
        self.assertEqual(alignment.sequences[2].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[3].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[3].seq), 174210431)
        self.assertEqual(alignment.sequences[3].seq[16169818:16169818+61], "ACAAAGGAAAACAATTAACCATAAGAAAAGTTTGTACCATAAAACAATTTTATTTTATTGA")
        self.assertEqual(alignment[3], "-------------ACAAAGGAAAACAATTAACCATAAGAAAAGTTTGTACCATAAAACAATTTTATTTTATTGA")
        self.assertEqual(alignment.sequences[3].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[3].annotations['leftCount'], 75)
        self.assertEqual(alignment.sequences[3].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[3].annotations['rightCount'], 97)
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
                numpy.array([[ 3018161,  3018172,  3018174,  3018176,  3018176,  3018230],
                             [16390178, 16390189, 16390189, 16390189, 16390189, 16390243],
                             [15871611, 15871622, 15871622, 15871622, 15871622, 15871676],
                             [16169818, 16169818, 16169818, 16169820, 16169825, 16169879]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 4781)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3018230:3018230+129], "AGGACAAAATAATACAGAtttttttttttttttttttttGCAGTACTGGAAATGGAATGAATGTCCCTCACAATCACTATCAAGGTCCCTATCAAGGCAATCACTCTGTCACCGAGCTACAGCCCCAGC")
        self.assertEqual(alignment[0], "AGGA-CAAAATAATACAGAtttttttttttttttttttttGCAGTACTGGAAATGGAATGAATGTCCCTCACAATCACTATCAAGGTCCCTATCAAGGCAATCACTCTGTCACCGAGCTA-CAGCCCCAGC")
        self.assertEqual(alignment.sequences[1].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[1].seq), 170899992)
        self.assertEqual(alignment.sequences[1].seq[15871676:15871676+95], "ATAATCCAATCAATATATATCAGAACCTGGCTCCCAATGTTCTGATAGTCATTATGAAAAAGAATTTACACATATATAGATTTATTAGGTATATG")
        self.assertEqual(alignment[1], "ATAATCCAATCAATATATAT---------------------CAGAACCTGGCTCCCAATG-----TTCTGATAGTCATTATGAA----------AAAGAATTTACACATATATAGATTTATTAGGTATATG")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[2].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[2].seq), 173908612)
        self.assertEqual(alignment.sequences[2].seq[16390243:16390243+95], "ATAATCCAATCAATATATATCAGAACCTGGCTCCCAATGTTCTGATAGTCATTGTGAAAAAGAatttacacatatatagatttattaggtatatg")
        self.assertEqual(alignment[2], "ATAATCCAATCAATATATAT---------------------CAGAACCTGGCTCCCAATG-----TTCTGATAGTCATTGTGAA----------AAAGAatttacacatatatagatttattaggtatatg")
        self.assertEqual(alignment.column_annotations["panTro2.chr6"], "99999999999999999999---------------------9999999999999999999-----9999999999999999999----------9999999999999999999999999999999999999")
        self.assertEqual(alignment.sequences[2].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[2].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['rightCount'], 0)
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
                numpy.array([[ 3018230,  3018234,  3018234,  3018249,  3018270,  3018289,
                               3018294,  3018313,  3018323,  3018349,  3018349,  3018359],
                             [15871676, 15871680, 15871681, 15871696, 15871696, 15871715,
                              15871715, 15871734, 15871734, 15871760, 15871761, 15871771],
                             [16390243, 16390247, 16390248, 16390263, 16390263, 16390282,
                              16390282, 16390301, 16390301, 16390327, 16390328, 16390338]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 61520)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3018359:3018359+123], "TTCAAACATGCATACATGCATTCATGTCTCATAATAATTATTAACATTGTCTTAGGCCAGAGGCTCGACTGCCCCAAAGCAATCCACTTAAACTGTCCCTGAGAAAGTCAttcctctccctaa")
        self.assertEqual(alignment[0], "TT-CAAACATGCATACATGCATTCATGTCTCATAA-TAATTATTAACA-TTGTCTTAGGCCAGAGGCTCGACTGCCCCAAAGCAATCCACT-------TAAACTGTCCCTGAGAA-AGTCAttcctctccctaa")
        self.assertEqual(alignment.sequences[1].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[1].seq), 173908612)
        self.assertEqual(alignment.sequences[1].seq[16390338:16390338+131], "tataaatgtacaaatatatgtgtatgtCTCATAAATAATTATTAACAGTCATCTTAGGTCAGTGGCCTGAGTGATACAAACTAAGCCATCCATATTTTATATTCTCTCCGGGAAGGGTCATCCTTTTCTTT")
        self.assertEqual(alignment[1], "-tataaatgtacaaatatatgtgtatgtCTCATAAATAATTATTAACAGTCATCTTAGGTCAGTGGCCTGAGTGATACAAACTAAGCCATCCATATTTTATATTCTCTCCGGGAAGGGTCATCCTTTTCTTT--")
        self.assertEqual(alignment.column_annotations["panTro2.chr6"], "-99999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999--")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 2400)
        self.assertEqual(alignment.sequences[2].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[2].seq), 170899992)
        self.assertEqual(alignment.sequences[2].seq[15871771:15871771+131], "TATAAACGTACAAATATATGTGTATGTCTCATAAATAATTATTAACAGTCATCTTAGGTCAGTGGCCTGAATGATACAAACTAAGCCATCCATATTTTATATTCTCTCCGGGAAGGGTCATCCTTTTCTTT")
        self.assertEqual(alignment[2], "-TATAAACGTACAAATATATGTGTATGTCTCATAAATAATTATTAACAGTCATCTTAGGTCAGTGGCCTGAATGATACAAACTAAGCCATCCATATTTTATATTCTCTCCGGGAAGGGTCATCCTTTTCTTT--")
        self.assertEqual(alignment.sequences[2].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[2].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[2].annotations['rightCount'], 2402)
        self.assertEqual(alignment.sequences[3].id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(alignment.sequences[3].seq), 359464)
        self.assertEqual(alignment.sequences[3].seq[184712:184712+127], "TTTATACATGTATGTGTAAATGTAGGTCCTGTAAGTAATTATTAACATTTGTCTTAGGTTAGAGGCCCGAGTGACACAAGCTAACCCATCCTATCCTCCCTGTGTGAAGGGTCATCCTTTCTTCTGA")
        self.assertEqual(alignment[3], "TT-TATACATGTATGTGTAAATGTAGGTCCTGTAAGTAATTATTAACATTTGTCTTAGGTTAGAGGCCCGAGTGACACAAGCTAACCCATCC------TATCCTCCCTGTGTGAAGGGTCATCCTTTCTTCTGA")
        self.assertEqual(alignment.column_annotations["otoGar1.scaffold_334.1-359464"], "99-99999999999999999999999999999999999999999999999999667736999999999999999999999999999995666------677755798899998999967967999999999589")
        self.assertEqual(alignment.sequences[3].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[3].annotations['leftCount'], 137)
        self.assertEqual(alignment.sequences[3].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[3].annotations['rightCount'], 2731)
        self.assertEqual(alignment.sequences[4].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[4].seq), 174210431)
        self.assertEqual(alignment.sequences[4].seq[16169976:16169976+129], "taaatgtacaaatatatgtgtatgtgtcataaataATTATTAACAGTCATCTTAGGTCAGTGGCCTGAGTGATAGAAACTAAGCCATCCATATTTTATATTCTCTCTGGGAAGGGTCATCCTTTTCTTT")
        self.assertEqual(alignment[4], "---taaatgtacaaatatatgtgtatgtgtcataaataATTATTAACAGTCATCTTAGGTCAGTGGCCTGAGTGATAGAAACTAAGCCATCCATATTTTATATTCTCTCTGGGAAGGGTCATCCTTTTCTTT--")
        self.assertEqual(alignment.sequences[4].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[4].annotations['leftCount'], 97)
        self.assertEqual(alignment.sequences[4].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[4].annotations['rightCount'], 2523)
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
                numpy.array([[ 3018359,  3018360,  3018361,  3018361,  3018393,  3018393,
                               3018405,  3018405,  3018447,  3018447,  3018447,  3018464,
                               3018464,  3018480,  3018482],
                             [16390338, 16390338, 16390339, 16390340, 16390372, 16390373,
                              16390385, 16390386, 16390428, 16390429, 16390435, 16390452,
                              16390453, 16390469, 16390469],
                             [15871771, 15871771, 15871772, 15871773, 15871805, 15871806,
                              15871818, 15871819, 15871861, 15871862, 15871868, 15871885,
                              15871886, 15871902, 15871902],
                             [  184712,   184713,   184714,   184714,   184746,   184747,
                                184759,   184760,   184802,   184803,   184803,   184820,
                                184821,   184837,   184839],
                             [16169976, 16169976, 16169976, 16169976, 16170008, 16170009,
                              16170021, 16170022, 16170064, 16170065, 16170071, 16170088,
                              16170089, 16170105, 16170105]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3018482:3018482+162], "tcttcatctcctcttttcctccttttttttttctcatttctctttctctttcttttgtccttttccttTATAGCAAGCAAGGCAAGTAGTCTCTATTTAGAAGGCATggagagaatggggagaggaggaaaggaggagaggggaggagaggaggggagGTAT")
        self.assertEqual(alignment[0], "tcttcatctcctcttttcctccttttttttttctcatttctctttctctttcttttgtccttttccttTATAGCAAGCAAGGCAAGTAGTCTCTATTTAGAAGGCATggagagaatggggagaggaggaaaggaggagaggggaggagaggaggggagGTAT")
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
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[3018482, 3018644]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1520)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3018644:3018644+178], "AGGGTAGGCCAAGTGCCTTGGGAAGTAGTTGTTGGTAGACTGAAAGTGTGTTCTGAGTGTCAGTGATGTTCATGAGATTATCACCAGCAAGGATGGCTGACGGGAACTGCAAGAGGCATAGCCCTGAGTTCTAAAGGAGAGGGAAACGTCACAGAAAGGATGCACTGTTTCAGCATCT")
        self.assertEqual(alignment[0], "AGGGTAGGCCAAGTGCCTTGGGAAGTAGTTGTTGGTAGACTGAAAGTGTGTTC---TGAGTGTCAGTGATGTTCA-TGAGATTATCACCAGCAAGGATG--GCTGACGGGAACTG---CAAGAGGCATAGCCCTGAGTTCTAAAGGAGAGGGAAACGTCACAGAAAGGATG--------------------------CACTGTTTCAGCATCT")
        self.assertEqual(alignment.sequences[1].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[1].seq), 125616256)
        self.assertEqual(alignment.sequences[1].seq[78070420:78070420+204], "AGAAAAGGCAAATGCCCCTGGGAGGGATAGCACATTGAGTTTATTCACATGCGTGTAAGTGATGGATACTGAGATCCTCTCCTAAGGATGCGGTTGATGGAATCAGAAGCAAATACCACCAGCCAGAACTCTAAAATGAGAGGAAAGATCATGGATAGAATAGAATTGTAATATGAAATATAAAATTTCCCTATTGAAGAATTT")
        self.assertEqual(alignment[1], "AGAAAAGGCAAATGCCCCTGGGAGGGA-----TAGCACATTGA--GTTTATTCACATGCGTGTAAGTGATGGATACTGAGATCCTCTCC--TAAGGATGCGGTTGATGGAATCAGAAGCAAATACCACCAGCCAGAACTCTAAAATGAGAGGAAAGATCATGGATAGAATAGAATTGTAATATGAAATATAAAATTTCCCTATTGAAGAATTT")
        self.assertEqual(alignment.column_annotations["canFam2.chr1"], "999999999999999999999999999-----99999999999--99999999999999999999999999999999999999999999--99999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'N')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 160)
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
                numpy.array([[ 3018644,  3018671,  3018676,  3018687,  3018689,  3018697,
                               3018697,  3018716,  3018716,  3018729,  3018731,  3018739,
                               3018739,  3018753,  3018753,  3018806,  3018806,  3018822],
                             [78070420, 78070447, 78070447, 78070458, 78070458, 78070466,
                              78070469, 78070488, 78070489, 78070502, 78070502, 78070510,
                              78070512, 78070526, 78070529, 78070582, 78070608, 78070624]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1986)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3018822:3018822+110], "CTGCCTTCCATTACGATTTACTGATCACTTACAACCCTCCCACAGAAGAGAACCTAACTTGCTTAGGAGCATATGTACAGTTAATCAAGACAAAAATAAGAATGGAGACt")
        self.assertEqual(alignment[0], "CTGCCTTCCATTACGATTTACTGATCACTTACAACCCTCCCACA----GAAGAGAACCTAACTTG-CTTAGGAGCATATGTACAGTTAATCAAGAC-----AAAAATAAGAATGGAGACt")
        self.assertEqual(alignment.sequences[1].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[1].seq), 125616256)
        self.assertEqual(alignment.sequences[1].seq[78070784:78070784+119], "CTTTTTTTCCCTGCCCATTATTACTCACTTAAAATTCTCCACATTGTAAAGGGAATTCAAATCGACTTCTAGAGATGCACACAATTTAGCAAGATCAACTAAAAATAAGTATGGAAAAT")
        self.assertEqual(alignment[1], "CTTTTTTTCCCTGCCCATTATTACTCACTTAAAATTCTCC-ACATTGTAAAGGGAATTCAAATCGACTTCTAGAGATGCACACAATTTAGCAAGATCAACTAAAAATAAGTATGGAAAAT")
        self.assertEqual(alignment.column_annotations["canFam2.chr1"], "9999999999999999999999999999999999999999-9999999999999999999999999999999999999999999999999999999999999999999999999999999")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 160)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 0)
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
                numpy.array([[ 3018822,  3018862,  3018863,  3018866,  3018866,  3018883,
                               3018883,  3018913,  3018913,  3018932],
                             [78070784, 78070824, 78070824, 78070827, 78070831, 78070848,
                              78070849, 78070879, 78070884, 78070903]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3018932:3018932+339], "gtctgagttagggttttactgctgtgaacagacaccatgaccaaggcatgtcttataaaaaaaatttaattagggctggcttacagattcagaggttcagtgggagcatcaaggtgggggcatggcagcatccaggcaggcatggtgcaggcagagctgagagttctacatcttcatccaaaggcttctagtggaagactgacttccaggcacctagggtgagggtcttaagcccacacccacagtgacacacctattccaaccaggtcacacctattccaacaaggccatacctccaaatggcaccactcctggtccaagaatatacaaaccatgaca")
        self.assertEqual(alignment[0], "gtctgagttagggttttactgctgtgaacagacaccatgaccaaggcatgtcttataaaaaaaatttaattagggctggcttacagattcagaggttcagtgggagcatcaaggtgggggcatggcagcatccaggcaggcatggtgcaggcagagctgagagttctacatcttcatccaaaggcttctagtggaagactgacttccaggcacctagggtgagggtcttaagcccacacccacagtgacacacctattccaaccaggtcacacctattccaacaaggccatacctccaaatggcaccactcctggtccaagaatatacaaaccatgaca")
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
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[3018932, 3019271]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 228)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3019271:3019271+106], "GAGACCAAATGTGGCGCTCACGTGAGGCCAGGAGTAAATCGCACACACAGCCCATGCTTTCACCATCTGCTAGGGTGCTCTGGAGCAGGGCAGGCTTCTAACCTGG")
        self.assertEqual(alignment[0], "GAGACCAAATG------------TGGCGCTCACG-TGAGGCCAGGAGTAAATCGCACACACAGCCCATGCTTTCACCATCTGCTAGGGTGCTCTGGAGCAGGGCAGGCTTCTAACCTGG")
        self.assertEqual(alignment.sequences[1].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[1].seq), 125616256)
        self.assertEqual(alignment.sequences[1].seq[78070903:78070903+106], "GAGGCAAATGGAAAATGCCCCACGATGTCCAGGCTGTGTCCATGTGTGACGCACATGCTGTCCCACCTGCTTCTATTCACTTGACACACATAGGTCCCAAACCTGG")
        self.assertEqual(alignment[1], "GAGGC-AAATGGAAAATGCCCCACGATGTCCAGGCTGTGTCCATGTGTGA--CGCACATGCTGTCC----------CACCTGCTTCTATTCACTTGACACACATAGGTCCCAAACCTGG")
        self.assertEqual(alignment.column_annotations["canFam2.chr1"], "99999-99999999999999999999999999999999999999999999--99999999999999----------9999999999999999999999999999999999999999999")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 0)
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
                numpy.array([[ 3019271,  3019276,  3019277,  3019282,  3019282,  3019293,
                               3019293,  3019308,  3019310,  3019324,  3019334,  3019377],
                             [78070903, 78070908, 78070908, 78070913, 78070925, 78070936,
                              78070937, 78070952, 78070952, 78070966, 78070966, 78071009]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 10938)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3019377:3019377+88], "CCCCAGCATTCTGGCAGACACAGTGAAAAGAGACAGATGGTCACTAATAAAATCTGTATAAATTAGATCTCAGAGGATGGATGGACCA")
        self.assertEqual(alignment[0], "CCCCAGCATTCTGGCAGACACAGTG-AAAAGAGACAGATGGTCACTAATAAAATCTGT-ATAAATTAG-ATCTCAGAGGATGGATGGACCA")
        self.assertEqual(alignment.sequences[1].id, "felCat3.scaffold_205680")
        self.assertEqual(len(alignment.sequences[1].seq), 119354)
        self.assertEqual(alignment.sequences[1].seq[72509:72509+88], "CCCAAGTGTTCTGATAGCTAATGTGAAAAAGAAGCATGTGCCCACCAGTAAGCTTTGTGGTGAACTAGAATCTCAGAGGATGGGACTC")
        self.assertEqual(alignment[1], "CCCAAGTGTTCTGATAGCTAATGTGAAAAAGAAGCATGTGCCCACCAGTAAGCTTTGTGGTGAACTAGAATCTCAGAGGATG---GGACTC")
        self.assertEqual(alignment.column_annotations["felCat3.scaffold_205680"], "9999999999999999999999999999999999999999999999999999999999999999999999999999999999---999999")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'N')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[2].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[2].seq), 125616256)
        self.assertEqual(alignment.sequences[2].seq[78071009:78071009+88], "CCCAAGTGTTCTGATTGCCTCTGTGAAAAAGAAACATGGGCCCGCTAATAagatttgcaatgacctagaatctcagaggatgggactc")
        self.assertEqual(alignment[2], "CCCAAGTGTTCTGATTGCCTCTGTGAAAAAGAAACATGGGCCCGCTAATAagatttgcaatgacctagaatctcagaggatg---ggactc")
        self.assertEqual(alignment.column_annotations["canFam2.chr1"], "9999999999999999999999999999999999999999999999999999999999999999999999999999999999---999999")
        self.assertEqual(alignment.sequences[2].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[2].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['rightCount'], 0)
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
                numpy.array([[ 3019377,  3019402,  3019402,  3019434,  3019434,  3019443,
                               3019443,  3019456,  3019459,  3019465],
                             [   72509,    72534,    72535,    72567,    72568,    72577,
                                 72578,    72591,    72591,    72597],
                             [78071009, 78071034, 78071035, 78071067, 78071068, 78071077,
                              78071078, 78071091, 78071091, 78071097]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 36924)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3019465:3019465+47], "AAGATAGATATTTAGAAGTAGCTTTTTATGTTTTTCTGATGTGTGTT")
        self.assertEqual(alignment[0], "AAGATAGATATTTAGAAGTAGCTTTTTATGTTTTTCTGATGTGTGTT")
        self.assertEqual(alignment.sequences[1].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[1].seq), 133105)
        self.assertEqual(alignment.sequences[1].seq[10128:10128+47], "aacaTCTATATTTTGAAATGGCTTTTCATGTTACTCTGATGTGTGTC")
        self.assertEqual(alignment[1], "aacaTCTATATTTTGAAATGGCTTTTCATGTTACTCTGATGTGTGTC")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 3479)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[2].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[2].seq), 173908612)
        self.assertEqual(alignment.sequences[2].seq[16392869:16392869+40], "atattttgaaatggcttttcatgttattctgatgTGTTTT")
        self.assertEqual(alignment[2], "-------atattttgaaatggcttttcatgttattctgatgTGTTTT")
        self.assertEqual(alignment.column_annotations["panTro2.chr6"], "-------9999999999999999999999999999999999999999")
        self.assertEqual(alignment.sequences[2].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[2].annotations['leftCount'], 2400)
        self.assertEqual(alignment.sequences[2].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[3].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[3].seq), 170899992)
        self.assertEqual(alignment.sequences[3].seq[15874304:15874304+40], "atattttgaaatggcttttcatgttattctgatgTGTTTT")
        self.assertEqual(alignment[3], "-------atattttgaaatggcttttcatgttattctgatgTGTTTT")
        self.assertEqual(alignment.sequences[3].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[3].annotations['leftCount'], 2402)
        self.assertEqual(alignment.sequences[3].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[3].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[4].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[4].seq), 125616256)
        self.assertEqual(alignment.sequences[4].seq[78071097:78071097+44], "aggtatatttaaaaatagctcttcatgttgttctaatgtgtgtt")
        self.assertEqual(alignment[4], "---aggtatatttaaaaatagctcttcatgttgttctaatgtgtgtt")
        self.assertEqual(alignment.column_annotations["canFam2.chr1"], "---99999999999999999999999999999999999999999999")
        self.assertEqual(alignment.sequences[4].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[4].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[4].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[4].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[5].id, "felCat3.scaffold_205680")
        self.assertEqual(len(alignment.sequences[5].seq), 119354)
        self.assertEqual(alignment.sequences[5].seq[72597:72597+43], "ACGTGTTTTTAAAAATAGTTTTCATGTTGTTCTGATGTGTGTT")
        self.assertEqual(alignment[5], "---ACGTGTTTTTAAAAATAG-TTTTCATGTTGTTCTGATGTGTGTT")
        self.assertEqual(alignment.column_annotations["felCat3.scaffold_205680"], "---999999999999999999-9999999999999999999999999")
        self.assertEqual(alignment.sequences[5].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[5].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[5].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[5].annotations['rightCount'], 193)
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
                numpy.array([[ 3019465,  3019468,  3019472,  3019486,  3019487,  3019512],
                             [   10128,    10131,    10135,    10149,    10150,    10175],
                             [16392869, 16392869, 16392869, 16392883, 16392884, 16392909],
                             [15874304, 15874304, 15874304, 15874318, 15874319, 15874344],
                             [78071097, 78071097, 78071101, 78071115, 78071116, 78071141],
                             [   72597,    72597,    72601,    72615,    72615,    72640]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 20303)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3019512:3019512+92], "TGCATCATTAAGACTAGAGTTCCTTTCTGTCTTTGCTTTCTTGACAGGGCCATGCTCGGCAGTCATTCTTAGACTGCTTTTTGTTTgtttgg")
        self.assertEqual(alignment[0], "TGCATCATTAAGACTAGAGTTCCT---------------TTCTGTCTT---TGCTTTCTTG--------ACAGGGCCATGCTCGGCAGTCATTCTTAGACTGCTTTTTGTTTgtttgg")
        self.assertEqual(alignment.sequences[1].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[1].seq), 133105)
        self.assertEqual(alignment.sequences[1].seq[10175:10175+86], "CCTCATGAAGGCCCAAGTTCCTAAAACATTAATTCTCTTCCTATTTCCTAGCTGTCTCGCCTAGGCCTTGCCCGCCAGCAATTCCC")
        self.assertEqual(alignment[1], "--CCTCATGAAGGCCCAAGTTCCTAAA-----ACATTAATTCTCTTCC---TATTTCCTAGCTGTCTCGCCTAGGCCTTGCCCGCCAGCAATTCCC----------------------")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[2].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[2].seq), 173908612)
        self.assertEqual(alignment.sequences[2].seq[16392909:16392909+93], "CTTCATTAAGGCCCAAGTTCCTAAAACATTCATTCTCTTCCTCTTTTCTAGAAAAGAGGTCTTGCCCGCCAGCAATTCCCACATGGGTATTGG")
        self.assertEqual(alignment[2], "--CTTCATTAAGGCCCAAGTTCCTAAA-----ACATTCATTCTCTTCC---TCTTTTCTAG------AAAAGAGGTCTTGCCCGCCAGCAATTCCCACATGGGTATTGG---------")
        self.assertEqual(alignment.column_annotations["panTro2.chr6"], "--9999999999999999999999999-----9999999999999999---9999999999------999999999999999999999999999999999999999999---------")
        self.assertEqual(alignment.sequences[2].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[2].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[3].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[3].seq), 170899992)
        self.assertEqual(alignment.sequences[3].seq[15874344:15874344+93], "CTTCATTAAGGCCCAAGTTCCTAAAACATTCATTCTCTTCCTCTTTTCTAGAAAAGAGGTCTTGCCCGCCAGCAATTCCCACATGGGTATTGG")
        self.assertEqual(alignment[3], "--CTTCATTAAGGCCCAAGTTCCTAAA-----ACATTCATTCTCTTCC---TCTTTTCTAG------AAAAGAGGTCTTGCCCGCCAGCAATTCCCACATGGGTATTGG---------")
        self.assertEqual(alignment.sequences[3].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[3].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[3].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[3].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[4].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[4].seq), 125616256)
        self.assertEqual(alignment.sequences[4].seq[78071141:78071141+97], "CCTCATCTAAACCCAGGTTCTTACAGGCTTATATTTTTTCTTTCTTCAACCCTTCCTCATCAGTGTCTTGCCTGCCAGTCATTCCTACTTTGTTCGG")
        self.assertEqual(alignment[4], "--CCTCATCTAAACCCAGGTTCTTACAGGCTTATATTTTTTCTTTCTTCAACCCTTCCTCA--------TCAGTGTCTTGCCTGCCAGTCATTCCTAC-----------TTTGTTCGG")
        self.assertEqual(alignment.column_annotations["canFam2.chr1"], "--99999999999999999999999999999999999999999999999999999999999--------99999999999999999999999999999-----------999999999")
        self.assertEqual(alignment.sequences[4].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[4].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[4].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[4].annotations['rightCount'], 0)
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
                numpy.array([[ 3019512,  3019514,  3019536,  3019536,  3019536,  3019536,
                               3019545,  3019545,  3019555,  3019555,  3019555,  3019582,
                               3019584,  3019595,  3019604],
                             [   10175,    10175,    10197,    10200,    10200,    10207,
                                 10216,    10216,    10226,    10232,    10234,    10261,
                                 10261,    10261,    10261],
                             [16392909, 16392909, 16392931, 16392934, 16392934, 16392941,
                              16392950, 16392950, 16392960, 16392960, 16392962, 16392989,
                              16392991, 16393002, 16393002],
                             [15874344, 15874344, 15874366, 15874369, 15874369, 15874376,
                              15874385, 15874385, 15874395, 15874395, 15874397, 15874424,
                              15874426, 15874437, 15874437],
                             [78071141, 78071141, 78071163, 78071166, 78071171, 78071178,
                              78071187, 78071190, 78071200, 78071200, 78071200, 78071227,
                              78071229, 78071229, 78071238]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3019604:3019604+98], "tttggtttggtttggttttttcaagacagggtttctttgtatagtcctagctgtcctggaactcactttgtagaccagactggccttgaactcagaaa")
        self.assertEqual(alignment[0], "tttggtttggtttggttttttcaagacagggtttctttgtatagtcctagctgtcctggaactcactttgtagaccagactggccttgaactcagaaa")
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
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[3019604, 3019702]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 45)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3019702:3019702+42], "tctgcctgcctctgcctcccaagtcctgggattaaaggcgtg")
        self.assertEqual(alignment[0], "tctgcctgcctctgcctcccaag--------------------------------tcctgggattaaaggcgtg")
        self.assertEqual(alignment.sequences[1].id, "felCat3.scaffold_205680")
        self.assertEqual(len(alignment.sequences[1].seq), 119354)
        self.assertEqual(alignment.sequences[1].seq[72833:72833+74], "tctgtctctctctctctctcaaaaataaacattaaaaaaaaCCAAACAAACAAACTCAAGTTCTTAAAGGTTTA")
        self.assertEqual(alignment[1], "tctgtctctctctctctctcaaaaataaacattaaaaaaaaCCAAACAAACAAACTCAAGTTCTTAAAGGTTTA")
        self.assertEqual(alignment.column_annotations["felCat3.scaffold_205680"], "99999999999999999999999999999999999999999999999999999999999999999999999999")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 193)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 0)
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
                numpy.array([[3019702, 3019725, 3019725, 3019744],
                             [  72833,   72856,   72888,   72907]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, -16865)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3019744:3019744+33], "cgccaccactgccctgcCTTAAACTGCTCTTAA")
        self.assertEqual(alignment[0], "-----------------cgccaccactgccctgcCT------------TAAACTGCTCTTAA")
        self.assertEqual(alignment.sequences[1].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[1].seq), 174210431)
        self.assertEqual(alignment.sequences[1].seq[16172628:16172628+27], "CGCCAGCAATTCCCACATGGGTATTGG")
        self.assertEqual(alignment[1], "-----------------CGCCAGCAATTCC------------------CACATGGGTATTGG")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 2523)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[2].id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(alignment.sequences[2].seq), 498454)
        self.assertEqual(alignment.sequences[2].seq[174138:174138+4], "TTGA")
        self.assertEqual(alignment[2], "----------------------------------------------------------TTGA")
        self.assertEqual(alignment.column_annotations["tupBel1.scaffold_114895.1-498454"], "----------------------------------------------------------9999")
        self.assertEqual(alignment.sequences[2].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[2].annotations['leftCount'], 2374)
        self.assertEqual(alignment.sequences[2].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[3].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[3].seq), 133105)
        self.assertEqual(alignment.sequences[3].seq[10261:10261+13], "ACATGGCTACTGG")
        self.assertEqual(alignment[3], "-------------------------------------------------ACATGGCTACTGG")
        self.assertEqual(alignment.sequences[3].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[3].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[3].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[3].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[4].id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(alignment.sequences[4].seq), 359464)
        self.assertEqual(alignment.sequences[4].seq[187570:187570+11], "ATTGGTGTTGA")
        self.assertEqual(alignment[4], "---------------------------------------------------ATTGGTGTTGA")
        self.assertEqual(alignment.column_annotations["otoGar1.scaffold_334.1-359464"], "---------------------------------------------------87784564678")
        self.assertEqual(alignment.sequences[4].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[4].annotations['leftCount'], 2731)
        self.assertEqual(alignment.sequences[4].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[4].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[5].id, "oryCun1.scaffold_156751")
        self.assertEqual(len(alignment.sequences[5].seq), 4726)
        self.assertEqual(alignment.sequences[5].seq[192:192+41], "CACCATTGTCTTGCCTGTCTCAGGTCCCATAGGGGTGTTGA")
        self.assertEqual(alignment[5], "--------------------CACCATTGTCTTGCCTGTC-TCAGGTCCCATAGGGGTGTTGA")
        self.assertEqual(alignment.column_annotations["oryCun1.scaffold_156751"], "--------------------9999999999999999999-9999999999999999999999")
        self.assertEqual(alignment.sequences[5].annotations['leftStatus'], 'N')
        self.assertEqual(alignment.sequences[5].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[5].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[5].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[6].id, "felCat3.scaffold_205680")
        self.assertEqual(len(alignment.sequences[6].seq), 119354)
        self.assertEqual(alignment.sequences[6].seq[72907:72907+62], "TGATTTCTCTTTCTCTTACTCATCAGTGTCTTGTCTTTCAGTCATTCCCACATTGGTGTAGA")
        self.assertEqual(alignment[6], "TGATTTCTCTTTCTCTTACTCATCAGTGTCTTGTCTTTCAGTCATTCCCACATTGGTGTAGA")
        self.assertEqual(alignment.column_annotations["felCat3.scaffold_205680"], "99999999999999999999999999999999999999999999999999999999999999")
        self.assertEqual(alignment.sequences[6].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[6].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[6].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[6].annotations['rightCount'], 0)
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
                numpy.array([[ 3019744,  3019744,  3019747,  3019757,  3019763,  3019763,
                               3019763,  3019763,  3019764,  3019766,  3019773,  3019777],
                             [16172628, 16172628, 16172631, 16172641, 16172641, 16172641,
                              16172641, 16172641, 16172642, 16172644, 16172651, 16172655],
                             [  174138,   174138,   174138,   174138,   174138,   174138,
                                174138,   174138,   174138,   174138,   174138,   174142],
                             [   10261,    10261,    10261,    10261,    10261,    10261,
                                 10261,    10261,    10261,    10263,    10270,    10274],
                             [  187570,   187570,   187570,   187570,   187570,   187570,
                                187570,   187570,   187570,   187570,   187577,   187581],
                             [     192,      192,      192,      202,      208,      211,
                                   211,      219,      220,      222,      229,      233],
                             [   72907,    72924,    72927,    72937,    72943,    72946,
                                 72947,    72955,    72956,    72958,    72965,    72969]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 367532)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3019777:3019777+183], "GGCAATGTCAGGCTATGCGTTCTAGACAGGGCACAAGAAAAGCTTTTAGCAGCAGAATAAACTTTTAAAGTAAATTACTTTCCTTGATAGCAACTAGACGACCCAATTGATACAGTGGAAAGAGGCCTTTGAGAATGCATGAGAGAATATTTCCTGTAAGAGTTGAACAATTTAGAATTTACc")
        self.assertEqual(alignment[0], "GGCAATGTCA-GGCTATGCGT------TCTAGACAGGGCACAAGAAAAGCTTTTAGCAGCAGAATAAACTTTT-AAAGTAAATTACTTTCCTTGATAGCAACTAGACGACCCAATTGA-TACAGT------GGAAAG-----A----------GGCCTTTGAGAAT---GCATGAGAGAATAT---TTCCTGTA------AGAGTTGAACAATTTAGAATTTACc")
        self.assertEqual(alignment.sequences[1].id, "oryCun1.scaffold_156751")
        self.assertEqual(len(alignment.sequences[1].seq), 4726)
        self.assertEqual(alignment.sequences[1].seq[233:233+193], "GGCAACATCATGACCTTATGTTCTAAGGGACAGGAAAAGCTTTTTCCAGTAGAATAAGCCTTTAGAGGAAACTGCTTCCCTTGCTAGTAATCAAGTGTTCAAAGTGATACAACAAAGAAAAAAAAAGCCTTTTGGGAATCAGGATGGAATAATATTACTTCCCCTATTTCTCAGAGCTAAATAGTTTAGGATT")
        self.assertEqual(alignment[1], "GGCAACATCATGACCTTATGT------TCTAA----GGGACAGGAAAAGCTTTTTCCAGTAGAATAAGCCTTT-AGAGGAAACTGCTTCCCTTGCTAGTAATCAAGTGTTCAAAGTGA-TACAACAAAGAAAAAAAA-----A----------GCCTTTTGGGAAT-CAGGATGGAATAATATTACTTCCCCTATTTCTCAGAGCTAAATAGTTTAGGATT----")
        self.assertEqual(alignment.column_annotations["oryCun1.scaffold_156751"], "999999999999999999999------99999----9999999999989999999998999999999999999-99999999999999999999999999999999999998989999-999899799999999999-----9----------9997999999998-978998999999999999999978999999979689999999999999999979----")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 37)
        self.assertEqual(alignment.sequences[2].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[2].seq), 170899992)
        self.assertEqual(alignment.sequences[2].seq[15874437:15874437+190], "GGCAACACCATGATACATATTCAAGATAAAGTACAGGAAAAGCTTTTTGCAGTGGAATAAAATGTTAAAGGGAATTGCTTTCCTTGATAGCAATCAGGTCACCAAAGTGAGATGAAGAAAAATTTTTTTGAGAATGCAGGATAGAATAATATTATTTCCCTTATTTCTAAGAGTTAAACAATTTAGGATT")
        self.assertEqual(alignment[2], "GGCAACACCATGA-TACATAT------TCAAGATAAAGTACAGGAAAAGCTTTTTGCAGTGGAATAAAATGTT-AAAGGGAATTGCTTTCCTTGATAGCAATCAGGTCACCAAAGTGA-GATGA-------AGAAAA-----A----------TTTTTTTGAGAATGCAGGATAGAATAATATTATTTCCCTTATTTCTAAGAGTTAAACAATTTAGGATT----")
        self.assertEqual(alignment.sequences[2].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[2].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[3].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[3].seq), 173908612)
        self.assertEqual(alignment.sequences[3].seq[16393002:16393002+191], "GGCAACACCATGATACATATTCAAGATAAAGTACAGGAAAAGCTTTTTGCAGTGGAATAAAATGTTAAAGGAAATTGCTTTCCTTGATAGCAATCAGGTCACCAAAGTGAGATGAAGAAAAATTTTTTTTGAGAATGCAGGATAGAATAATATTATTTCCCTTATTTCTAAGAGTTAAACAATTTAGGATT")
        self.assertEqual(alignment[3], "GGCAACACCATGA-TACATAT------TCAAGATAAAGTACAGGAAAAGCTTTTTGCAGTGGAATAAAATGTT-AAAGGAAATTGCTTTCCTTGATAGCAATCAGGTCACCAAAGTGA-GATGA-------AGAAAA-----A---------TTTTTTTTGAGAATGCAGGATAGAATAATATTATTTCCCTTATTTCTAAGAGTTAAACAATTTAGGATT----")
        self.assertEqual(alignment.column_annotations["panTro2.chr6"], "9999999999999-9999999------9999999999999999999999999999999999999999999999-99999999999999999999999999999999999999999999-99999-------999999-----9---------999999999999999999999999999999999999999999999999999999999999999999999----")
        self.assertEqual(alignment.sequences[3].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[3].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[3].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[3].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[4].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[4].seq), 174210431)
        self.assertEqual(alignment.sequences[4].seq[16172655:16172655+190], "GGCAACACCATGATACATATTCAAGATAAAGTACAGGAAAAGCTTTTTGCGGTGGAATAAAATGTTAAAGGAAATTGCTTTCCTTGATAGCAATCAGGTCACCAAAGTGAGATGAAGAAAAATTCTTTTTGAGAATGCAGCACAGAATAATATTATTTCCCTTATTTCTAGAGTTAAACAATTTAGGATT")
        self.assertEqual(alignment[4], "GGCAACACCATGA-TACATAT------TCAAGATAAAGTACAGGAAAAGCTTTTTGCGGTGGAATAAAATGTT-AAAGGAAATTGCTTTCCTTGATAGCAATCAGGTCACCAAAGTGA-GATGA-------AGAAAA-----A---------TTCTTTTTGAGAATGCAGCACAGAATAATATTATTTCCCTTATTTCT-AGAGTTAAACAATTTAGGATT----")
        self.assertEqual(alignment.sequences[4].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[4].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[4].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[4].annotations['rightCount'], 21)
        self.assertEqual(alignment.sequences[5].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[5].seq), 133105)
        self.assertEqual(alignment.sequences[5].seq[10274:10274+187], "GGCAACACCGTGATACATATCCGAGATAAAGTACAGGAAAAGCTTTTGCAGTGGACGAAAATGTTGAAGGAAATTGCTTTTCTTGATAGCAATCAGGTCACCAAAGTGAGATGAAGAAAGTTTTTTGAGAATGTGGGATAGAATAATATTGTTTTCGTTATTTCTAAGAGTTAAACAATTTAGGATT")
        self.assertEqual(alignment[5], "GGCAACACCGTGA-TACATAT------CCGAGATAAAGTACAGGAAAAGC-TTTTGCAGTGGACGAAAATGTT-GAAGGAAATTGCTTTTCTTGATAGCAATCAGGTCACCAAAGTGA-GATGA-------AGAAA-----------------GTTTTTTGAGAATGTGGGATAGAATAATATTGTTTTCGTTATTTCTAAGAGTTAAACAATTTAGGATT----")
        self.assertEqual(alignment.sequences[5].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[5].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[5].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[5].annotations['rightCount'], 21)
        self.assertEqual(alignment.sequences[6].id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(alignment.sequences[6].seq), 359464)
        self.assertEqual(alignment.sequences[6].seq[187581:187581+197], "GGCAACGCCGTGATCCTTGCTGTAGATAAAGACAAAGTGCAGGAAAAGCCCTTTACGGTAGAATAAAAGGTTAAAGGAAATTGCTTTTCTTGATAGCAATCAGGTGACCAAAGTGATTATAAAGAAAAAATGTTTTGAGAATGCTGGAAAGAATAATTTTATCCACCCTCCATTTCAGAGTTAAATAATTCAGGATT")
        self.assertEqual(alignment[6], "GGCAACGCCGTGA-TCCTTGCTGTAGATAAAGACAAAGTGCAGGAAAAGCCCTTTACGGTAGAATAAAAGGTT-AAAGGAAATTGCTTTTCTTGATAGCAATCAGGTGACCAAAGTGATTATAA-------AGAAAA-----A----------ATGTTTTGAGAATGCTGGAAAGAATAATTTTATCCACCCTCCATTTCAGAGTTAAATAATTCAGGATT----")
        self.assertEqual(alignment.column_annotations["otoGar1.scaffold_334.1-359464"], "5466765345564-89977768969797886996999557766667768677897764536465578775526-67697777697779896767576676677696675376877786786858-------478887-----7----------56858859779666569356677477759968759378767656876769557597655596536566----")
        self.assertEqual(alignment.sequences[6].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[6].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[6].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[6].annotations['rightCount'], 358)
        self.assertEqual(alignment.sequences[7].id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(alignment.sequences[7].seq), 498454)
        self.assertEqual(alignment.sequences[7].seq[174142:174142+196], "GGCAATGTCATGAATGTTCTAGATAAAGTACAGGAAAAGCTTTTTGCAGTAGGATAAACTTTTAAAGGAAATTGCCTTCCTTGACAGCAATCAGGTGACCAAAGTCATACTGAAAAAAAAAACAAAAACAAACTTTTGTGGATGCAGTATAGAAGCATATTATTTTCCCTATGTCTAAGTTAAATAACTTGGGCTT")
        self.assertEqual(alignment[7], "GGCAATGTCATGA----ATGT------TCTAGATAAAGTACAGGAAAAGCTTTTTGCAGTAGGATAAACTTTT-AAAGGAAATTGCCTTCCTTGACAGCAATCAGGTGACCAAAGTCATACTGA-------AAAAAA-----AAACAAAAACAAACTTTTGTGGATGCAGTATAGAAGCATATTATTTTCCCTATGTCT--AAGTTAAATAACTTGGGCTT----")
        self.assertEqual(alignment.column_annotations["tupBel1.scaffold_114895.1-498454"], "9999999999999----9999------9999999999999999999999999999999999999999999999-99999999999999999999999999999999999999999999999999-------999999-----999999999999999999999999999999999999999999999999999999999--99999999999999999999----")
        self.assertEqual(alignment.sequences[7].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[7].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[7].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[7].annotations['rightCount'], 623)
        self.assertEqual(alignment.sequences[8].id, "felCat3.scaffold_205680")
        self.assertEqual(len(alignment.sequences[8].seq), 119354)
        self.assertEqual(alignment.sequences[8].seq[72969:72969+186], "GACATCATAAAGTACATTCTAGATAAAGTGCAGGAAAAGCATTTTGAGTAGAATAAAATGTTGTAGGAAATTGCTTTCTTCAAGAGACATTAGGTGATTAAAGTAATATAAAGGGAATTTTTTTTGAGAATACAGGATACAATAATATTATTTCCCCTACATCTAAGAGttaaataatttaaaatt")
        self.assertEqual(alignment[8], "GACATC---ATAA-AGTACAT------TCTAGATAAAGTGCAGGAAAAGCATTTTG-AGTAGAATAAAATGTT-GTAGGAAATTGCTTTCTTCAAGAGACATTAGGTGATTAAAGTAA-TATAA-------AGGGAA-----T----------TTTTTTTGAGAATACAGGATACAATAATATTATTTCCCCTACATCTAAGAGttaaataatttaaaatt----")
        self.assertEqual(alignment.column_annotations["felCat3.scaffold_205680"], "999999---9999-9999999------99999999999999999999999999999-9999999999999999-99999999999999999999999999999999999999999999-99999-------999999-----9----------99999999999999999999999999999999999999999999999999999999999999999999----")
        self.assertEqual(alignment.sequences[8].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[8].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[8].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[8].annotations['rightCount'], 27)
        self.assertEqual(alignment.sequences[9].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[9].seq), 125616256)
        self.assertEqual(alignment.sequences[9].seq[78071238:78071238+193], "ATCATGATGTAAATTCCGGGTAAAGTACAACAAAAGGATTTTGAACAGGATAAAATGTTAATAGGAAATTACTTTCTTTGAGAGCAATCAGGTGATTGAAGGAATATAAAGGGAAAAAACAAAGTTCAGTGAATCCAAGTTAGAATAATATTGTTTCCCTTACTTCTAAGAGttaaataatttaaaatttaat")
        self.assertEqual(alignment[9], "------ATCATGA-TGTAAAT------TCCGGGTAAAGTACAACAAAAGGATTTTG-AACAGGATAAAATGTTAATAGGAAATTACTTTCTTTGAGAGCAATCAGGTGATTGAAGGAA-TATAA-------AGGGAAAAAACA----------AAGTTCAGTGAATCCAAGTTAGAATAATATTGTTTCCCTTACTTCTAAGAGttaaataatttaaaatttaat")
        self.assertEqual(alignment.column_annotations["canFam2.chr1"], "------9999999-9999999------99999999999999999999999999999-9999999999999999999999999999999999999999999999999999999999999-99999-------999999999999----------999999999999999999999999999999999999999999999999999999999999999999999999")
        self.assertEqual(alignment.sequences[9].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[9].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[9].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[9].annotations['rightCount'], 0)
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
                numpy.array([[ 3019777,  3019783,  3019786,  3019787,  3019787,  3019789,
                               3019790,  3019793,  3019797,  3019797,  3019802,  3019806,
                               3019820,  3019821,  3019826,  3019827,  3019843,  3019843,
                               3019887,  3019887,  3019892,  3019893,  3019893,  3019898,
                               3019899,  3019899,  3019900,  3019900,  3019900,  3019913,
                               3019913,  3019913,  3019927,  3019927,  3019935,  3019935,
                               3019935,  3019936,  3019956,  3019960],
                             [     233,      239,      242,      243,      244,      246,
                                   247,      250,      254,      254,      259,      259,
                                   273,      274,      279,      280,      296,      296,
                                   340,      340,      345,      346,      352,      357,
                                   358,      358,      359,      359,      359,      372,
                                   372,      374,      388,      391,      399,      404,
                                   405,      406,      426,      426],
                             [15874437, 15874443, 15874446, 15874447, 15874448, 15874450,
                              15874450, 15874453, 15874457, 15874457, 15874462, 15874466,
                              15874480, 15874481, 15874486, 15874487, 15874503, 15874503,
                              15874547, 15874547, 15874552, 15874552, 15874552, 15874557,
                              15874558, 15874558, 15874559, 15874559, 15874559, 15874572,
                              15874573, 15874575, 15874589, 15874592, 15874600, 15874605,
                              15874606, 15874607, 15874627, 15874627],
                             [16393002, 16393008, 16393011, 16393012, 16393013, 16393015,
                              16393015, 16393018, 16393022, 16393022, 16393027, 16393031,
                              16393045, 16393046, 16393051, 16393052, 16393068, 16393068,
                              16393112, 16393112, 16393117, 16393117, 16393117, 16393122,
                              16393123, 16393123, 16393124, 16393124, 16393125, 16393138,
                              16393139, 16393141, 16393155, 16393158, 16393166, 16393171,
                              16393172, 16393173, 16393193, 16393193],
                             [16172655, 16172661, 16172664, 16172665, 16172666, 16172668,
                              16172668, 16172671, 16172675, 16172675, 16172680, 16172684,
                              16172698, 16172699, 16172704, 16172705, 16172721, 16172721,
                              16172765, 16172765, 16172770, 16172770, 16172770, 16172775,
                              16172776, 16172776, 16172777, 16172777, 16172778, 16172791,
                              16172792, 16172794, 16172808, 16172811, 16172819, 16172824,
                              16172824, 16172825, 16172845, 16172845],
                             [   10274,    10280,    10283,    10284,    10285,    10287,
                                 10287,    10290,    10294,    10294,    10299,    10303,
                                 10317,    10317,    10322,    10323,    10339,    10339,
                                 10383,    10383,    10388,    10388,    10388,    10393,
                                 10393,    10393,    10393,    10393,    10393,    10406,
                                 10407,    10409,    10423,    10426,    10434,    10439,
                                 10440,    10441,    10461,    10461],
                             [  187581,   187587,   187590,   187591,   187592,   187594,
                                187594,   187597,   187601,   187607,   187612,   187616,
                                187630,   187631,   187636,   187637,   187653,   187653,
                                187697,   187698,   187703,   187703,   187703,   187708,
                                187709,   187709,   187710,   187710,   187710,   187723,
                                187724,   187726,   187740,   187743,   187751,   187756,
                                187757,   187758,   187778,   187778],
                             [  174142,   174148,   174151,   174152,   174153,   174155,
                                174155,   174155,   174159,   174159,   174164,   174168,
                                174182,   174183,   174188,   174189,   174205,   174205,
                                174249,   174250,   174255,   174255,   174255,   174260,
                                174261,   174261,   174262,   174271,   174272,   174285,
                                174286,   174288,   174302,   174305,   174313,   174318,
                                174318,   174318,   174338,   174338],
                             [   72969,    72975,    72975,    72976,    72977,    72979,
                                 72979,    72982,    72986,    72986,    72991,    72995,
                                 73009,    73010,    73015,    73015,    73031,    73031,
                                 73075,    73075,    73080,    73080,    73080,    73085,
                                 73086,    73086,    73087,    73087,    73087,    73100,
                                 73101,    73103,    73117,    73120,    73128,    73133,
                                 73134,    73135,    73155,    73155],
                             [78071238, 78071238, 78071241, 78071242, 78071243, 78071245,
                              78071245, 78071248, 78071252, 78071252, 78071257, 78071261,
                              78071275, 78071276, 78071281, 78071281, 78071297, 78071298,
                              78071342, 78071342, 78071347, 78071347, 78071347, 78071352,
                              78071353, 78071358, 78071359, 78071359, 78071359, 78071372,
                              78071373, 78071375, 78071389, 78071392, 78071400, 78071405,
                              78071406, 78071407, 78071427, 78071431]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3019960:3019960+757], "actagggatgggagaggctcccagaacccagtaatgatgacattaagaaatacacaacagttgggaaatggaacccaaagagaacacctccagtagataagcatgacccccagttgagggatgggcccatgcacccatcttaaaattttggacccagaattattcttctcaaaaggaaatgcagggatgaaaatggagcagagactggaagaaaggccaaccagagactgccctaactcaggatccatcgcatgtgcaggcaccaaccccaacactattgctgatgccatgttgtacttgctgatggaagcctggcatggctgtcctctgagagtctcaaatgaggcacctgacagatgcagatacttacagccaaccaatggactgagccccgggacctcaataaaagaatgaggggatggcaaccccataggaagaacaacagtatcaactccctggactcctcagagctcccggggactaagccaccaactaaagagcatacataggctgctctgaggccccagatacatatgtagcagaggactgcctcagtgggaggggatgtgcttggtcttgtgaaggcttgatgctccagagaaggaggatgctagaggggtgaggtgggagtggatgggtgggtgggcaggggagcaccctcttagaggacaagggctctggggtgggggagctcatggagggggaactgggaaggagggagaacatttgaaatgtaaataaataaaataataaaaaa")
        self.assertEqual(alignment[0], "actagggatgggagaggctcccagaacccagtaatgatgacattaagaaatacacaacagttgggaaatggaacccaaagagaacacctccagtagataagcatgacccccagttgagggatgggcccatgcacccatcttaaaattttggacccagaattattcttctcaaaaggaaatgcagggatgaaaatggagcagagactggaagaaaggccaaccagagactgccctaactcaggatccatcgcatgtgcaggcaccaaccccaacactattgctgatgccatgttgtacttgctgatggaagcctggcatggctgtcctctgagagtctcaaatgaggcacctgacagatgcagatacttacagccaaccaatggactgagccccgggacctcaataaaagaatgaggggatggcaaccccataggaagaacaacagtatcaactccctggactcctcagagctcccggggactaagccaccaactaaagagcatacataggctgctctgaggccccagatacatatgtagcagaggactgcctcagtgggaggggatgtgcttggtcttgtgaaggcttgatgctccagagaaggaggatgctagaggggtgaggtgggagtggatgggtgggtgggcaggggagcaccctcttagaggacaagggctctggggtgggggagctcatggagggggaactgggaaggagggagaacatttgaaatgtaaataaataaaataataaaaaa")
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
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[3019960, 3020717]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 8951)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3020717:3020717+44], "TGTCAAACATGCATAAAGATATACTGAGGAGCCCATGAATTTTA")
        self.assertEqual(alignment[0], "TGTCAAACATGCATAAAGATATACT-GAGGAGCCCATGAATTTTA")
        self.assertEqual(alignment.sequences[1].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[1].seq), 125616256)
        self.assertEqual(alignment.sequences[1].seq[78071431:78071431+33], "tGTTTAACATAATAGAGGGGCCCATGAATTTTA")
        self.assertEqual(alignment[1], "tGTT------------TAACATAATAGAGGGGCCCATGAATTTTA")
        self.assertEqual(alignment.column_annotations["canFam2.chr1"], "9999------------99999999999999999999999999999")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 9)
        self.assertEqual(alignment.sequences[2].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[2].seq), 173908612)
        self.assertEqual(alignment.sequences[2].seq[16393193:16393193+38], "AAAAATTTTAAATATACTAGAGGGGTCCATGAATTTTA")
        self.assertEqual(alignment[2], "----AAAAAT---TTTAAATATACTAGAGGGGTCCATGAATTTTA")
        self.assertEqual(alignment.column_annotations["panTro2.chr6"], "----999999---99999999999999999999999999999999")
        self.assertEqual(alignment.sequences[2].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[2].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[2].annotations['rightCount'], 11)
        self.assertEqual(alignment.sequences[3].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[3].seq), 170899992)
        self.assertEqual(alignment.sequences[3].seq[15874627:15874627+38], "AAAAATTTTAAATATACTAGAGGGGTCCATGAATTTTA")
        self.assertEqual(alignment[3], "----AAAAAT---TTTAAATATACTAGAGGGGTCCATGAATTTTA")
        self.assertEqual(alignment.sequences[3].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[3].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[3].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[3].annotations['rightCount'], 11)
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
                numpy.array([[ 3020717,  3020721,  3020727,  3020730,  3020733,  3020742,
                               3020742,  3020761],
                             [78071431, 78071435, 78071435, 78071435, 78071435, 78071444,
                              78071445, 78071464],
                             [16393193, 16393193, 16393199, 16393199, 16393202, 16393211,
                              16393212, 16393231],
                             [15874627, 15874627, 15874633, 15874633, 15874636, 15874645,
                              15874646, 15874665]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3020761:3020761+157], "TATATATGCTATCCGTGTGCTGTGATTTTTGTTTTAAATGTTATTTTATGTATATGcaagattttgcattgtagcagaaggtggcttcaaactcacgatcctcctgcctcagccttccaagtgctgagatcatacctctgcaccatcctgcccACCT")
        self.assertEqual(alignment[0], "TATATATGCTATCCGTGTGCTGTGATTTTTGTTTTAAATGTTATTTTATGTATATGcaagattttgcattgtagcagaaggtggcttcaaactcacgatcctcctgcctcagccttccaagtgctgagatcatacctctgcaccatcctgcccACCT")
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
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[3020761, 3020918]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 85471)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3020918:3020918+96], "GAGGTTTGTGACTTTTAATACTGATTGTTATCTAACATCACAGAATTCTCAGTTCTTAAGGAAACAATTGTTCTGTGTGTTATTTGTCTAGGAGGA")
        self.assertEqual(alignment[0], "GAGGTTTGTGACTTTTAATA----------CTGATTGTTATCTAACATCACAGAATTCTCAGTTCTTAAGGAAACAATTGTTCTGTGTGTTATTTGTCTAGGAGGA")
        self.assertEqual(alignment.sequences[1].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[1].seq), 170899992)
        self.assertEqual(alignment.sequences[1].seq[15874676:15874676+73], "ATGGTTATTTAATACTGCACAATCCTCAGCTTTTAAGGAAAAACATGTTTGACTTATTAATTATTTAGGACAA")
        self.assertEqual(alignment[1], "---------------------------------ATGGTTATTTAATACTGCACAATCCTCAGCTTTTAAGGAAAAACATGTTTGACTTATTAATTATTTAGGACAA")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 11)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[2].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[2].seq), 173908612)
        self.assertEqual(alignment.sequences[2].seq[16393242:16393242+73], "ATGGTTATTTAATACTGCACAATCCTCAGCTTTTAAGGAAAAACATGTTTGACTTATTAATTATTTAGGACAA")
        self.assertEqual(alignment[2], "---------------------------------ATGGTTATTTAATACTGCACAATCCTCAGCTTTTAAGGAAAAACATGTTTGACTTATTAATTATTTAGGACAA")
        self.assertEqual(alignment.column_annotations["panTro2.chr6"], "---------------------------------9999999999999999999999999999999999999999999999999999999999999999999999999")
        self.assertEqual(alignment.sequences[2].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[2].annotations['leftCount'], 11)
        self.assertEqual(alignment.sequences[2].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[3].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[3].seq), 125616256)
        self.assertEqual(alignment.sequences[3].seq[78071473:78071473+75], "TAATAAAATAATTCTGTACCATTATGAAGTTCTCAGCTCTTAAGGAAAAAGTGATTTGTTAATTACTTAGGACAA")
        self.assertEqual(alignment[3], "--------------------------TAATAAAATAATTCTGTACCATTATGAAGTTCTCAGCTCTTAAGGAAAAA-----GTGATTTGTTAATTACTTAGGACAA")
        self.assertEqual(alignment.column_annotations["canFam2.chr1"], "--------------------------99999999999999999999999999999999999999999999999999-----9999999999999999999999999")
        self.assertEqual(alignment.sequences[3].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[3].annotations['leftCount'], 9)
        self.assertEqual(alignment.sequences[3].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[3].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[4].id, "felCat3.scaffold_205680")
        self.assertEqual(len(alignment.sequences[4].seq), 119354)
        self.assertEqual(alignment.sequences[4].seq[73182:73182+76], "GTCCATGAATTTTAATGGTGAAGTAATGAAATAATTCTTTAATATCACAAAGTAATTTATTAATTATTTAGGACAA")
        self.assertEqual(alignment[4], "---GTCCATGAATTTTAATGGTGAAGTAATGAAATAATTCTTTAATATCAC----------------------AAA-----GTAATTTATTAATTATTTAGGACAA")
        self.assertEqual(alignment.column_annotations["felCat3.scaffold_205680"], "---999989999999999999999999999999999999999999999999----------------------999-----9999999999999999999999769")
        self.assertEqual(alignment.sequences[4].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[4].annotations['leftCount'], 27)
        self.assertEqual(alignment.sequences[4].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[4].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[5].id, "cavPor2.scaffold_216473")
        self.assertEqual(len(alignment.sequences[5].seq), 10026)
        self.assertEqual(alignment.sequences[5].seq[9420:9420+82], "GAGTTCTGCAAATtttaataggaaaagaagatTAATATTACAAAGTTTTCAGTTCTGATGGAAaaagtgtgatttttttaaa")
        self.assertEqual(alignment[5], "GAGTTCTGCAAATtttaata----------ggaaaagaagatTAATATTACAAAGTTTTCAGTTCTGATGGAAaaa----------gtgtgatttttttaaa----")
        self.assertEqual(alignment.column_annotations["cavPor2.scaffold_216473"], "99679966799668566465----------8567828838346568788268865663676558367765776667----------6965686997776684----")
        self.assertEqual(alignment.sequences[5].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[5].annotations['leftCount'], 1372)
        self.assertEqual(alignment.sequences[5].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[5].annotations['rightCount'], 28)
        self.assertEqual(alignment.sequences[6].id, "oryCun1.scaffold_156751")
        self.assertEqual(len(alignment.sequences[6].seq), 4726)
        self.assertEqual(alignment.sequences[6].seq[463:463+102], "AGATCCATGAATTTTAATAAAGCAATGAGATAATTATTTAATGTTGGAAAATCCTCAAGTCTTAAGGAAAAAAGTGTTTGACTTGTTAATGACTTGGGATGA")
        self.assertEqual(alignment[6], "-AGATCCATGAATTTTAATA---AAGCAATGAGATAATTATTTAATGTTGGAAAATCCTCAAGTCTTAAGGAAAAAAGTGTTTGACTTGTTAATGACTTGGGATGA")
        self.assertEqual(alignment.column_annotations["oryCun1.scaffold_156751"], "-9999997999999999999---99999999999999999999999999999999998677999969999999999999999999999999999999999999999")
        self.assertEqual(alignment.sequences[6].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[6].annotations['leftCount'], 37)
        self.assertEqual(alignment.sequences[6].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[6].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[7].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[7].seq), 133105)
        self.assertEqual(alignment.sequences[7].seq[10482:10482+100], "GGGGTCCATGAATTTTAATAGTAACAAAATGGTTATTCAATATTGCAAAATCCTCAGCTTTAAGGAAAAACATATTTGATTTGTTAATTATTTAGGACAA")
        self.assertEqual(alignment[7], "GGGGTCCATGAATTTTAATA-----GTAACAAAATGGTTATTCAATATTGCAAAATCCTCAGC-TTTAAGGAAAAACATATTTGATTTGTTAATTATTTAGGACAA")
        self.assertEqual(alignment.sequences[7].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[7].annotations['leftCount'], 21)
        self.assertEqual(alignment.sequences[7].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[7].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[8].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[8].seq), 174210431)
        self.assertEqual(alignment.sequences[8].seq[16172866:16172866+101], "GGGGTCCATGAATTTTAATAGTAATAAAATGGTTATTTAATATTGCAAAATCCTCAGCTTTTAAGGAAAAACGTGTTTGACTTATTAATTATTTAGGACAA")
        self.assertEqual(alignment[8], "GGGGTCCATGAATTTTAATA-----GTAATAAAATGGTTATTTAATATTGCAAAATCCTCAGCTTTTAAGGAAAAACGTGTTTGACTTATTAATTATTTAGGACAA")
        self.assertEqual(alignment.sequences[8].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[8].annotations['leftCount'], 21)
        self.assertEqual(alignment.sequences[8].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[8].annotations['rightCount'], 0)
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
                numpy.array([[ 3020918,  3020919,  3020921,  3020938,  3020938,  3020938,
                               3020938,  3020938,  3020941,  3020959,  3020971,  3020972,
                               3020981,  3020984,  3020989,  3020994,  3021010,  3021014],
                             [15874676, 15874676, 15874676, 15874676, 15874676, 15874676,
                              15874676, 15874676, 15874676, 15874694, 15874706, 15874707,
                              15874716, 15874719, 15874724, 15874729, 15874745, 15874749],
                             [16393242, 16393242, 16393242, 16393242, 16393242, 16393242,
                              16393242, 16393242, 16393242, 16393260, 16393272, 16393273,
                              16393282, 16393285, 16393290, 16393295, 16393311, 16393315],
                             [78071473, 78071473, 78071473, 78071473, 78071473, 78071473,
                              78071473, 78071477, 78071480, 78071498, 78071510, 78071511,
                              78071520, 78071523, 78071523, 78071528, 78071544, 78071548],
                             [   73182,    73182,    73182,    73199,    73202,    73204,
                                 73205,    73209,    73212,    73230,    73230,    73230,
                                 73230,    73233,    73233,    73238,    73254,    73258],
                             [    9420,     9421,     9423,     9440,     9440,     9440,
                                  9440,     9440,     9443,     9461,     9473,     9474,
                                  9483,     9486,     9486,     9486,     9502,     9502],
                             [     463,      463,      465,      482,      482,      484,
                                   485,      489,      492,      510,      522,      523,
                                   532,      535,      540,      545,      561,      565],
                             [   10482,    10483,    10485,    10502,    10502,    10502,
                                 10503,    10507,    10510,    10528,    10540,    10540,
                                 10549,    10552,    10557,    10562,    10578,    10582],
                             [16172866, 16172867, 16172869, 16172886, 16172886, 16172886,
                              16172887, 16172891, 16172894, 16172912, 16172924, 16172925,
                              16172934, 16172937, 16172942, 16172947, 16172963, 16172967]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 105724)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3021014:3021014+40], "ACCTTGGTGACGCCACTGGATTTTGTATGACTGAATACTG")
        self.assertEqual(alignment[0], "ACCTTGGTGACGCCACTGGATTTTGTATGACTGAATACTG")
        self.assertEqual(alignment.sequences[1].id, "oryCun1.scaffold_156751")
        self.assertEqual(len(alignment.sequences[1].seq), 4726)
        self.assertEqual(alignment.sequences[1].seq[565:565+40], "ATCTTGATGAAGTCATTCCAACTTGGATGATTTAGGAATT")
        self.assertEqual(alignment[1], "ATCTTGATGAAGTCATTCCAACTTGGATGATTTAGGAATT")
        self.assertEqual(alignment.column_annotations["oryCun1.scaffold_156751"], "9999999999999969999699999999999999999999")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[2].id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(alignment.sequences[2].seq), 359464)
        self.assertEqual(alignment.sequences[2].seq[188136:188136+40], "ATCTTGGTGAAGTTATTCCAATTTATGTGATTTAGGAATG")
        self.assertEqual(alignment[2], "ATCTTGGTGAAGTTATTCCAATTTATGTGATTTAGGAATG")
        self.assertEqual(alignment.column_annotations["otoGar1.scaffold_334.1-359464"], "9999999999999999999999999999999999999999")
        self.assertEqual(alignment.sequences[2].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[2].annotations['leftCount'], 358)
        self.assertEqual(alignment.sequences[2].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[3].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[3].seq), 133105)
        self.assertEqual(alignment.sequences[3].seq[10582:10582+40], "ATTTTGGTGAAGTTATTCCAACTTGTGTGGCTTAGGAATG")
        self.assertEqual(alignment[3], "ATTTTGGTGAAGTTATTCCAACTTGTGTGGCTTAGGAATG")
        self.assertEqual(alignment.sequences[3].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[3].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[3].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[3].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[4].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[4].seq), 170899992)
        self.assertEqual(alignment.sequences[4].seq[15874749:15874749+40], "ATTTTGGTGAAGTTATTCCAACTTGCATGGCTTAGGAATG")
        self.assertEqual(alignment[4], "ATTTTGGTGAAGTTATTCCAACTTGCATGGCTTAGGAATG")
        self.assertEqual(alignment.sequences[4].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[4].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[4].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[4].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[5].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[5].seq), 173908612)
        self.assertEqual(alignment.sequences[5].seq[16393315:16393315+40], "ATTTTGGTGAAGTTATTCCAACTTGCATGGCTTAGGAATG")
        self.assertEqual(alignment[5], "ATTTTGGTGAAGTTATTCCAACTTGCATGGCTTAGGAATG")
        self.assertEqual(alignment.column_annotations["panTro2.chr6"], "9999999999999999999999999999999999999999")
        self.assertEqual(alignment.sequences[5].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[5].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[5].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[5].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[6].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[6].seq), 174210431)
        self.assertEqual(alignment.sequences[6].seq[16172967:16172967+40], "ATTTTGGTGAAGTTATTCCAACTTGCATGGCTTAGGAATG")
        self.assertEqual(alignment[6], "ATTTTGGTGAAGTTATTCCAACTTGCATGGCTTAGGAATG")
        self.assertEqual(alignment.sequences[6].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[6].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[6].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[6].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[7].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[7].seq), 125616256)
        self.assertEqual(alignment.sequences[7].seq[78071548:78071548+40], "ATCTCAGTGAAGTCATTCTGACTCGCATGATTTAGGAATG")
        self.assertEqual(alignment[7], "ATCTCAGTGAAGTCATTCTGACTCGCATGATTTAGGAATG")
        self.assertEqual(alignment.column_annotations["canFam2.chr1"], "9999999999999999999999999999999999999999")
        self.assertEqual(alignment.sequences[7].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[7].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[7].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[7].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[8].id, "felCat3.scaffold_205680")
        self.assertEqual(len(alignment.sequences[8].seq), 119354)
        self.assertEqual(alignment.sequences[8].seq[73258:73258+40], "ATCTCAGTGAAGTTATTCTGACTTGCATGATTTAGGTATG")
        self.assertEqual(alignment[8], "ATCTCAGTGAAGTTATTCTGACTTGCATGATTTAGGTATG")
        self.assertEqual(alignment.column_annotations["felCat3.scaffold_205680"], "9999989999989988999997999979997996167779")
        self.assertEqual(alignment.sequences[8].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[8].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[8].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[8].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[9].id, "dasNov1.scaffold_56749")
        self.assertEqual(len(alignment.sequences[9].seq), 10470)
        self.assertEqual(alignment.sequences[9].seq[3116:3116+40], "ACATCAGTGAAATCATTCCGACTCGTATGACTGAGCGATG")
        self.assertEqual(alignment[9], "ACATCAGTGAAATCATTCCGACTCGTATGACTGAGCGATG")
        self.assertEqual(alignment.column_annotations["dasNov1.scaffold_56749"], "9759855999977756667495765475885678385647")
        self.assertEqual(alignment.sequences[9].annotations['leftStatus'], 'N')
        self.assertEqual(alignment.sequences[9].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[9].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[9].annotations['rightCount'], 0)
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
                alignment.coordinates,
                numpy.array([[ 3021014,  3021054],
                             [     565,      605],
                             [  188136,   188176],
                             [   10582,    10622],
                             [15874749, 15874789],
                             [16393315, 16393355],
                             [16172967, 16173007],
                             [78071548, 78071588],
                             [   73258,    73298],
                             [    3116,     3156]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 115790)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3021054:3021054+50], "CTCATTTGGGAACTTACAGGTCAGCAAAGGCTTCCAGGACTTACATGCAG")
        self.assertEqual(alignment[0], "CTCATTTGGGAACTTACAGGTCAGCAAAGGCTTCCAG--------------------GACTTACATGCAG")
        self.assertEqual(alignment.sequences[1].id, "cavPor2.scaffold_216473")
        self.assertEqual(len(alignment.sequences[1].seq), 10026)
        self.assertEqual(alignment.sequences[1].seq[9530:9530+38], "CTCATTTGGGAGTTCAGAGTTGCTATAAATACTTGTGA")
        self.assertEqual(alignment[1], "CTCATTTGGGAGTTCAGAGTT--------GCTATAAA--------------------TACTTGTGA----")
        self.assertEqual(alignment.column_annotations["cavPor2.scaffold_216473"], "788777667898846665666--------98766677--------------------876665669----")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 28)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[2].id, "oryCun1.scaffold_156751")
        self.assertEqual(len(alignment.sequences[2].seq), 4726)
        self.assertEqual(alignment.sequences[2].seq[605:605+69], "CTAATTTGGGCTCTTACATTTTATAACTGCCTAATGGAGTCTTCTTTTTTCCTCCTTAATTACACGTAT")
        self.assertEqual(alignment[2], "CTAATTTGGGCTCTTACATTTTATAACTGCCTAATGGAGTCTTCTTTTTTCCT-CCTTAATTACACGTAT")
        self.assertEqual(alignment.column_annotations["oryCun1.scaffold_156751"], "99999999999997999799999999999999999999999799799999997-9979999999495999")
        self.assertEqual(alignment.sequences[2].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[2].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[3].id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(alignment.sequences[3].seq), 359464)
        self.assertEqual(alignment.sequences[3].seq[188176:188176+69], "CTAACTTGGGGACACAGGTTCATAAAGGCTTAATGGAGCCCACTGGTCCCTTCCCCAGACTCTATGTAC")
        self.assertEqual(alignment[3], "CTAACTTGGGGACACAGG-TTCATAAAGGCTTAATGGAGCCCACTGGTCCCTTCCCCAGACTCTATGTAC")
        self.assertEqual(alignment.column_annotations["otoGar1.scaffold_334.1-359464"], "999999999999999999-999999999999999999998999999689999999999999999999999")
        self.assertEqual(alignment.sequences[3].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[3].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[3].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[3].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[4].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[4].seq), 133105)
        self.assertEqual(alignment.sequences[4].seq[10622:10622+70], "CCAGTTTGGGGACTTAGATTTTCTAACTGCCTAATGAAGTCTGCTCTTTCCTTCCCTAGCCTATATGTAT")
        self.assertEqual(alignment[4], "CCAGTTTGGGGACTTAGATTTTCTAACTGCCTAATGAAGTCTGCTCTTTCCTTCCCTAGCCTATATGTAT")
        self.assertEqual(alignment.sequences[4].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[4].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[4].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[4].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[5].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[5].seq), 170899992)
        self.assertEqual(alignment.sequences[5].seq[15874789:15874789+70], "ACAATTTGGGGATTTACATTTCCTAACTGCCTAATGAAGTCTGCTCCTTCCTTCCCTAGACTATATGTAT")
        self.assertEqual(alignment[5], "ACAATTTGGGGATTTACATTTCCTAACTGCCTAATGAAGTCTGCTCCTTCCTTCCCTAGACTATATGTAT")
        self.assertEqual(alignment.sequences[5].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[5].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[5].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[5].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[6].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[6].seq), 173908612)
        self.assertEqual(alignment.sequences[6].seq[16393355:16393355+70], "CCAATTTGGGGATTTACATTTCCTAACTGCCTAATGAAGTCTGCTCTTTCCTTCCCTAGACTATATGTAT")
        self.assertEqual(alignment[6], "CCAATTTGGGGATTTACATTTCCTAACTGCCTAATGAAGTCTGCTCTTTCCTTCCCTAGACTATATGTAT")
        self.assertEqual(alignment.column_annotations["panTro2.chr6"], "9999999999999999999999999999999999999999999999999999999999999999999999")
        self.assertEqual(alignment.sequences[6].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[6].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[6].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[6].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[7].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[7].seq), 174210431)
        self.assertEqual(alignment.sequences[7].seq[16173007:16173007+70], "CCAATTTGGGGATTTACATTTCCTAACTGCCTAATGAAGTCTGCTCTTTCCTTCCCTAGACTATATGTAT")
        self.assertEqual(alignment[7], "CCAATTTGGGGATTTACATTTCCTAACTGCCTAATGAAGTCTGCTCTTTCCTTCCCTAGACTATATGTAT")
        self.assertEqual(alignment.sequences[7].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[7].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[7].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[7].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[8].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[8].seq), 125616256)
        self.assertEqual(alignment.sequences[8].seq[78071588:78071588+70], "CTAATTTGGTGACTCAAATTTCATAATTGCTTAATGAAAGCTGATGTTTTCTTCCCTACACCATATGTAT")
        self.assertEqual(alignment[8], "CTAATTTGGTGACTCAAATTTCATAATTGCTTAATGAAAGCTGATGTTTTCTTCCCTACACCATATGTAT")
        self.assertEqual(alignment.column_annotations["canFam2.chr1"], "9999999999999999999999999999999999999999999999999999999999999999999999")
        self.assertEqual(alignment.sequences[8].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[8].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[8].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[8].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[9].id, "felCat3.scaffold_205680")
        self.assertEqual(len(alignment.sequences[9].seq), 119354)
        self.assertEqual(alignment.sequences[9].seq[73298:73298+70], "CTAATTTGGTGACTTGAATTTCATAATTGCTTAAGGAAGCCTGCTGTTTCCTTCCCAAGACTATATGCAT")
        self.assertEqual(alignment[9], "CTAATTTGGTGACTTGAATTTCATAATTGCTTAAGGAAGCCTGCTGTTTCCTTCCCAAGACTATATGCAT")
        self.assertEqual(alignment.column_annotations["felCat3.scaffold_205680"], "9769999999975699868977669777966666596959759669595666758736585676666655")
        self.assertEqual(alignment.sequences[9].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[9].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[9].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[9].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[10].id, "dasNov1.scaffold_56749")
        self.assertEqual(len(alignment.sequences[10].seq), 10470)
        self.assertEqual(alignment.sequences[10].seq[3156:3156+70], "CTAATTGGGGGACTCACATTTTAAAATTGCTAACCGAAGCCTGCTGTTTTCTTCTCTAGACTATAAGTAT")
        self.assertEqual(alignment[10], "CTAATTGGGGGACTCACATTTTAAAATTGCTAACCGAAGCCTGCTGTTTTCTTCTCTAGACTATAAGTAT")
        self.assertEqual(alignment.column_annotations["dasNov1.scaffold_56749"], "5556576999664654656985688667655565647767537567688856666555556565555656")
        self.assertEqual(alignment.sequences[10].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[10].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[10].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[10].annotations['rightCount'], 0)
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
                numpy.array([[ 3021054,  3021072,  3021073,  3021075,  3021083,  3021091,
                               3021091,  3021091,  3021091,  3021100,  3021104],
                             [    9530,     9548,     9549,     9551,     9551,     9559,
                                  9559,     9559,     9559,     9568,     9568],
                             [     605,      623,      624,      626,      634,      642,
                                   658,      658,      661,      670,      674],
                             [  188176,   188194,   188194,   188196,   188204,   188212,
                                188228,   188229,   188232,   188241,   188245],
                             [   10622,    10640,    10641,    10643,    10651,    10659,
                                 10675,    10676,    10679,    10688,    10692],
                             [15874789, 15874807, 15874808, 15874810, 15874818, 15874826,
                              15874842, 15874843, 15874846, 15874855, 15874859],
                             [16393355, 16393373, 16393374, 16393376, 16393384, 16393392,
                              16393408, 16393409, 16393412, 16393421, 16393425],
                             [16173007, 16173025, 16173026, 16173028, 16173036, 16173044,
                              16173060, 16173061, 16173064, 16173073, 16173077],
                             [78071588, 78071606, 78071607, 78071609, 78071617, 78071625,
                              78071641, 78071642, 78071645, 78071654, 78071658],
                             [   73298,    73316,    73317,    73319,    73327,    73335,
                                 73351,    73352,    73355,    73364,    73368],
                             [    3156,     3174,     3175,     3177,     3185,     3193,
                                  3209,     3210,     3213,     3222,     3226]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 44222)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3021104:3021104+32], "CTGTTAGTGCTGTTTTAATGTACCTCGCAGTA")
        self.assertEqual(alignment[0], "CTGTTAGTGCTGTTTT---AATGTACCTCGCAGTA")
        self.assertEqual(alignment.sequences[1].id, "cavPor2.scaffold_216473")
        self.assertEqual(len(alignment.sequences[1].seq), 10026)
        self.assertEqual(alignment.sequences[1].seq[9568:9568+27], "ATTACTTTCTTAATGTAAATTTTTATA")
        self.assertEqual(alignment[1], "-----ATTACTTTCTT---AATGTAAATTTTTATA")
        self.assertEqual(alignment.column_annotations["cavPor2.scaffold_216473"], "-----55455868746---7957999989884575")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[2].id, "oryCun1.scaffold_156751")
        self.assertEqual(len(alignment.sequences[2].seq), 4726)
        self.assertEqual(alignment.sequences[2].seq[674:674+30], "TTGCTAGTATTCTCAATGTAATATTTTATA")
        self.assertEqual(alignment[2], "TTGCTAGTA--TTCTC---AATGTAATATTTTATA")
        self.assertEqual(alignment.column_annotations["oryCun1.scaffold_156751"], "999979999--96665---6999999999999999")
        self.assertEqual(alignment.sequences[2].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[2].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[3].id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(alignment.sequences[3].seq), 359464)
        self.assertEqual(alignment.sequences[3].seq[188245:188245+32], "TTGCTGGTGTTCCCTTAATATAATGTTTTATA")
        self.assertEqual(alignment[3], "TTGCTGGTGTTCCCTT---AATATAATGTTTTATA")
        self.assertEqual(alignment.column_annotations["otoGar1.scaffold_334.1-359464"], "9999999999999999---9999999999999999")
        self.assertEqual(alignment.sequences[3].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[3].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[3].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[3].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[4].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[4].seq), 133105)
        self.assertEqual(alignment.sequences[4].seq[10692:10692+30], "TTACTCGTGCCCTTAATATAGCATTTTATA")
        self.assertEqual(alignment[4], "TTACTCGTG--CCCTT---AATATAGCATTTTATA")
        self.assertEqual(alignment.sequences[4].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[4].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[4].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[4].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[5].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[5].seq), 170899992)
        self.assertEqual(alignment.sequences[5].seq[15874859:15874859+30], "TTGCTGGTACTCTTAATATCACATTTTATA")
        self.assertEqual(alignment[5], "TTGCTGGTA--CTCTT---AATATCACATTTTATA")
        self.assertEqual(alignment.sequences[5].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[5].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[5].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[5].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[6].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[6].seq), 173908612)
        self.assertEqual(alignment.sequences[6].seq[16393425:16393425+30], "TTGCTGGTACTCTTAATATCACATTTTATA")
        self.assertEqual(alignment[6], "TTGCTGGTA--CTCTT---AATATCACATTTTATA")
        self.assertEqual(alignment.column_annotations["panTro2.chr6"], "999999999--99999---9999999999999999")
        self.assertEqual(alignment.sequences[6].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[6].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[6].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[6].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[7].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[7].seq), 174210431)
        self.assertEqual(alignment.sequences[7].seq[16173077:16173077+30], "TTGCTGGTGCTCTTAATATAACATTTTATA")
        self.assertEqual(alignment[7], "TTGCTGGTG--CTCTT---AATATAACATTTTATA")
        self.assertEqual(alignment.sequences[7].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[7].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[7].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[7].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[8].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[8].seq), 125616256)
        self.assertEqual(alignment.sequences[8].seq[78071658:78071658+25], "TTGCTGGTGCTCTCAAAATGTAGCA")
        self.assertEqual(alignment[8], "TTGCTGGTGCTCTCAA---AATGTAGCA-------")
        self.assertEqual(alignment.column_annotations["canFam2.chr1"], "9999999999999999---999999999-------")
        self.assertEqual(alignment.sequences[8].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[8].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[8].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[8].annotations['rightCount'], 196)
        self.assertEqual(alignment.sequences[9].id, "felCat3.scaffold_205680")
        self.assertEqual(len(alignment.sequences[9].seq), 119354)
        self.assertEqual(alignment.sequences[9].seq[73368:73368+25], "TTGCTGGTGCTCTCAAAATATAACA")
        self.assertEqual(alignment[9], "TTGCTGGTGCTCTCAA---AATATAACA-------")
        self.assertEqual(alignment.column_annotations["felCat3.scaffold_205680"], "8677668566555658---876555655-------")
        self.assertEqual(alignment.sequences[9].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[9].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[9].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[9].annotations['rightCount'], 6)
        self.assertEqual(alignment.sequences[10].id, "dasNov1.scaffold_56749")
        self.assertEqual(len(alignment.sequences[10].seq), 10470)
        self.assertEqual(alignment.sequences[10].seq[3226:3226+33], "TTGCTGGTGATCTGATTAATGTAACATTTTATG")
        self.assertEqual(alignment[10], "TTGCTGGTG--ATCTGATTAATGTAACATTTTATG")
        self.assertEqual(alignment.column_annotations["dasNov1.scaffold_56749"], "856647736--775356546747663745776545")
        self.assertEqual(alignment.sequences[10].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[10].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[10].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[10].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[11].id, "echTel1.scaffold_288249")
        self.assertEqual(len(alignment.sequences[11].seq), 100002)
        self.assertEqual(alignment.sequences[11].seq[95225:95225+21], "CTGTTAATGCTCTGTTTTATG")
        self.assertEqual(alignment[11], "CTGTTAATG--CTCTG------------TTTTATG")
        self.assertEqual(alignment.column_annotations["echTel1.scaffold_288249"], "999999999--99999------------9999999")
        self.assertEqual(alignment.sequences[11].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[11].annotations['leftCount'], 7564)
        self.assertEqual(alignment.sequences[11].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[11].annotations['rightCount'], 0)
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
                numpy.array([[ 3021104,  3021109,  3021113,  3021115,  3021120,  3021120,
         3021129,  3021136],
                             [    9568,     9568,     9572,     9574,     9579,     9579,
                                  9588,     9595],
                             [     674,      679,      683,      683,      688,      688,
                                   697,      704],
                             [  188245,   188250,   188254,   188256,   188261,   188261,
                                188270,   188277],
                             [   10692,    10697,    10701,    10701,    10706,    10706,
                                 10715,    10722],
                             [15874859, 15874864, 15874868, 15874868, 15874873, 15874873,
                              15874882, 15874889],
                             [16393425, 16393430, 16393434, 16393434, 16393439, 16393439,
                              16393448, 16393455],
                             [16173077, 16173082, 16173086, 16173086, 16173091, 16173091,
                              16173100, 16173107],
                             [78071658, 78071663, 78071667, 78071669, 78071674, 78071674,
                              78071683, 78071683],
                             [   73368,    73373,    73377,    73379,    73384,    73384,
                                 73393,    73393],
                             [    3226,     3231,     3235,     3235,     3240,     3243,
                                  3252,     3259],
                             [   95225,    95230,    95234,    95234,    95239,    95239,
                                 95239,    95246]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 43757)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3021136:3021136+44], "AGGCAAATGAGGTGATAAGATTGTGTTTACTCCCTCTGTGCTTG")
        self.assertEqual(alignment[0], "AGGCAAATGAGGTGATAAGA-------TTGTGTT-----TAC----TCCCTCTGTGC----------TTG")
        self.assertEqual(alignment.sequences[1].id, "cavPor2.scaffold_216473")
        self.assertEqual(len(alignment.sequences[1].seq), 10026)
        self.assertEqual(alignment.sequences[1].seq[9595:9595+43], "AGCAGATGAGATAACAGTGATATATTcattctctgtgtgtatg")
        self.assertEqual(alignment[1], "AG-CAGATGAGATAACAGTG-------ATATATT-----cat----tctctgtgtgt----------atg")
        self.assertEqual(alignment.column_annotations["cavPor2.scaffold_216473"], "49-96766988786798889-------8679566-----756----66896967579----------888")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[2].id, "oryCun1.scaffold_156751")
        self.assertEqual(len(alignment.sequences[2].seq), 4726)
        self.assertEqual(alignment.sequences[2].seq[704:704+42], "AAGCAAATGAGGTGATGAGGGCATCTATTTCTTATTTGTGTG")
        self.assertEqual(alignment[2], "AAGCAAATGAGGTGATGAGG---------GCATC-----TAT----TTCTTATTTGT----------GTG")
        self.assertEqual(alignment.column_annotations["oryCun1.scaffold_156751"], "99998999999999999999---------99999-----999----99999999999----------999")
        self.assertEqual(alignment.sequences[2].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[2].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[3].id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(alignment.sequences[3].seq), 359464)
        self.assertEqual(alignment.sequences[3].seq[188277:188277+45], "AAACAAATGAGATCACAGGACATATGTACTTGTCCCCCACGTGTG")
        self.assertEqual(alignment[3], "AAACAAATGAGATCACA-GG-------ACATATG-----TA--CTTGTCCCCCACGT----------GTG")
        self.assertEqual(alignment.column_annotations["otoGar1.scaffold_334.1-359464"], "99999999999999999-99-------9999999-----99--99939999999999----------999")
        self.assertEqual(alignment.sequences[3].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[3].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[3].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[3].annotations['rightCount'], 15)
        self.assertEqual(alignment.sequences[4].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[4].seq), 133105)
        self.assertEqual(alignment.sequences[4].seq[10722:10722+45], "AAGCAAATGAGATCACAGCACATGTATATTTTTTCTCCGTGTGTG")
        self.assertEqual(alignment[4], "AAGCAAATGAGATCACA----------GCACATG-----TATATTTTTTCTCCGTGT----------GTG")
        self.assertEqual(alignment.sequences[4].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[4].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[4].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[4].annotations['rightCount'], 15)
        self.assertEqual(alignment.sequences[5].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[5].seq), 170899992)
        self.assertEqual(alignment.sequences[5].seq[15874889:15874889+47], "AAACAAATGAGATCACAAGGGCATATATGTTTTTTTCTCTGTGTGTG")
        self.assertEqual(alignment[5], "AAACAAATGAGATCACAAGG-------GCATATA-----TGT-TTTTTTCTCTGTGT----------GTG")
        self.assertEqual(alignment.sequences[5].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[5].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[5].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[5].annotations['rightCount'], 15)
        self.assertEqual(alignment.sequences[6].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[6].seq), 173908612)
        self.assertEqual(alignment.sequences[6].seq[16393455:16393455+47], "AAACAAATGAGATCACAAGGGCATATATGTTTTTTTCTCTGTGTGTG")
        self.assertEqual(alignment[6], "AAACAAATGAGATCACAAGG-------GCATATA-----TGT-TTTTTTCTCTGTGT----------GTG")
        self.assertEqual(alignment.column_annotations["panTro2.chr6"], "99999999999999999999-------9999999-----999-99999999999999----------999")
        self.assertEqual(alignment.sequences[6].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[6].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[6].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[6].annotations['rightCount'], 15)
        self.assertEqual(alignment.sequences[7].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[7].seq), 174210431)
        self.assertEqual(alignment.sequences[7].seq[16173107:16173107+47], "AAACAAATGAGATCACAAGGGCATATGTATTTTTTTCTCTGTGTGTG")
        self.assertEqual(alignment[7], "AAACAAATGAGATCACAAGG-------GCATATG-----TAT-TTTTTTCTCTGTGT----------GTG")
        self.assertEqual(alignment.sequences[7].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[7].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[7].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[7].annotations['rightCount'], 15)
        self.assertEqual(alignment.sequences[8].id, "eriEur1.scaffold_266115")
        self.assertEqual(len(alignment.sequences[8].seq), 4589)
        self.assertEqual(alignment.sequences[8].seq[358:358+65], "taaataataagataagatgaaagCATAGCATGTATTTTCTtgccctctccttctctgtctctgtc")
        self.assertEqual(alignment[8], "taaataataagataagatgaaagCATAGCATGTA-----TTTTCTtgccctctccttctctgtctctgtc")
        self.assertEqual(alignment.column_annotations["eriEur1.scaffold_266115"], "9999999999999999999999999999999999-----9999999999999999999999999999999")
        self.assertEqual(alignment.sequences[8].annotations['leftStatus'], 'N')
        self.assertEqual(alignment.sequences[8].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[8].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[8].annotations['rightCount'], 9)
        self.assertEqual(alignment.sequences[9].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[9].seq), 125616256)
        self.assertEqual(alignment.sequences[9].seq[78071879:78071879+51], "tAGAAAATGAGATACCGTGTGTATATGCTTTCTCTCTGTGTGGGGTCTGTG")
        self.assertEqual(alignment[9], "tAGAAAATGAGATACC-----------GTGTGTA-----TATGCTTTCTCTCTGTGT---GGGGTCTGTG")
        self.assertEqual(alignment.column_annotations["canFam2.chr1"], "9999999999999999-----------9999999-----999999999999999999---9999999999")
        self.assertEqual(alignment.sequences[9].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[9].annotations['leftCount'], 196)
        self.assertEqual(alignment.sequences[9].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[9].annotations['rightCount'], 9)
        self.assertEqual(alignment.sequences[10].id, "felCat3.scaffold_205680")
        self.assertEqual(len(alignment.sequences[10].seq), 119354)
        self.assertEqual(alignment.sequences[10].seq[73399:73399+51], "AAGCAAATGAGATAAGACACGTGTATTCTTCCTCTCTGTGTGGTGTCTGTG")
        self.assertEqual(alignment[10], "AAGCAAATGAGATAAG-----------ACACGTG-----TATTCTTCCTCTCTGTGT---GGTGTCTGTG")
        self.assertEqual(alignment.column_annotations["felCat3.scaffold_205680"], "9755596656454353-----------5346324-----243451535354422233---3635339999")
        self.assertEqual(alignment.sequences[10].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[10].annotations['leftCount'], 6)
        self.assertEqual(alignment.sequences[10].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[10].annotations['rightCount'], 9)
        self.assertEqual(alignment.sequences[11].id, "dasNov1.scaffold_56749")
        self.assertEqual(len(alignment.sequences[11].seq), 10470)
        self.assertEqual(alignment.sequences[11].seq[3259:3259+46], "AGGAAAGTGAGAGGACCAGGGCAAATAAATATTCTTTCTCTGTGTA")
        self.assertEqual(alignment[11], "AGGAAAGTGAGAGGACCAGG-------GCAAATA-----AATATTCTTTCTCTGTGT----------A--")
        self.assertEqual(alignment.column_annotations["dasNov1.scaffold_56749"], "55569764556657658568-------6768555-----858686457859535656----------5--")
        self.assertEqual(alignment.sequences[11].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[11].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[11].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[11].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[12].id, "echTel1.scaffold_288249")
        self.assertEqual(len(alignment.sequences[12].seq), 100002)
        self.assertEqual(alignment.sequences[12].seq[95246:95246+51], "AAGAAGGTGAGATGACAAGGGTGTATAGATAGGATATTCTTGCTTTGGGTG")
        self.assertEqual(alignment[12], "AAGAAGGTGAGATGACAAGG-------GTGTATAGATAGGATATTCTTGCTTTGGGT----------G--")
        self.assertEqual(alignment.column_annotations["echTel1.scaffold_288249"], "99999999999999999999-------999999999999999999999999999999----------9--")
        self.assertEqual(alignment.sequences[12].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[12].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[12].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[12].annotations['rightCount'], 0)
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
                numpy.array([[ 3021136,  3021138,  3021139,  3021152,  3021153,  3021154,
                               3021156,  3021156,  3021158,  3021163,  3021163,  3021165,
                               3021166,  3021166,  3021166,  3021177,  3021177,  3021177,
                               3021178,  3021180],
                             [    9595,     9597,     9597,     9610,     9611,     9612,
                                  9614,     9614,     9616,     9621,     9621,     9623,
                                  9624,     9624,     9624,     9635,     9635,     9635,
                                  9636,     9638],
                             [     704,      706,      707,      720,      721,      722,
                                   724,      724,      724,      729,      729,      731,
                                   732,      732,      732,      743,      743,      743,
                                   744,      746],
                             [  188277,   188279,   188280,   188293,   188294,   188294,
                                188296,   188296,   188298,   188303,   188303,   188305,
                                188305,   188305,   188308,   188319,   188319,   188319,
                                188320,   188322],
                             [   10722,    10724,    10725,    10738,    10739,    10739,
                                 10739,    10739,    10741,    10746,    10746,    10748,
                                 10749,    10750,    10753,    10764,    10764,    10764,
                                 10765,    10767],
                             [15874889, 15874891, 15874892, 15874905, 15874906, 15874907,
                              15874909, 15874909, 15874911, 15874916, 15874916, 15874918,
                              15874919, 15874919, 15874922, 15874933, 15874933, 15874933,
                              15874934, 15874936],
                             [16393455, 16393457, 16393458, 16393471, 16393472, 16393473,
                              16393475, 16393475, 16393477, 16393482, 16393482, 16393484,
                              16393485, 16393485, 16393488, 16393499, 16393499, 16393499,
                              16393500, 16393502],
                             [16173107, 16173109, 16173110, 16173123, 16173124, 16173125,
                              16173127, 16173127, 16173129, 16173134, 16173134, 16173136,
                              16173137, 16173137, 16173140, 16173151, 16173151, 16173151,
                              16173152, 16173154],
                             [     358,      360,      361,      374,      375,      376,
                                   378,      385,      387,      392,      392,      394,
                                   395,      396,      399,      410,      413,      420,
                                   421,      423],
                             [78071879, 78071881, 78071882, 78071895, 78071895, 78071895,
                              78071895, 78071895, 78071897, 78071902, 78071902, 78071904,
                              78071905, 78071906, 78071909, 78071920, 78071920, 78071927,
                              78071928, 78071930],
                             [   73399,    73401,    73402,    73415,    73415,    73415,
                                 73415,    73415,    73417,    73422,    73422,    73424,
                                 73425,    73426,    73429,    73440,    73440,    73447,
                                 73448,    73450],
                             [    3259,     3261,     3262,     3275,     3276,     3277,
                                  3279,     3279,     3281,     3286,     3286,     3288,
                                  3289,     3290,     3293,     3304,     3304,     3304,
                                  3305,     3305],
                             [   95246,    95248,    95249,    95262,    95263,    95264,
                                 95266,    95266,    95268,    95273,    95278,    95280,
                                 95281,    95282,    95285,    95296,    95296,    95296,
                                 95297,    95297]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 32886)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3021180:3021180+24], "TCCCAGAGAGTCTGATAGGAGGAG")
        self.assertEqual(alignment[0], "-------------------TCCC-------AGAGAGTCTGA-TAGGAGGAG")
        self.assertEqual(alignment.sequences[1].id, "cavPor2.scaffold_216473")
        self.assertEqual(len(alignment.sequences[1].seq), 10026)
        self.assertEqual(alignment.sequences[1].seq[9638:9638+32], "tgcctTTGTATAGTGGCTCTGAGTATAAAGTA")
        self.assertEqual(alignment[1], "-------------------tgcctTTGTATAGTGGCTCTGAGTATAAAGTA")
        self.assertEqual(alignment.column_annotations["cavPor2.scaffold_216473"], "-------------------67576649966655666885655548785776")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[2].id, "oryCun1.scaffold_156751")
        self.assertEqual(len(alignment.sequences[2].seq), 4726)
        self.assertEqual(alignment.sequences[2].seq[746:746+27], "TGCCTGTGTATAGCATTTCTGAATATA")
        self.assertEqual(alignment[2], "-------------------TGCCTGTGTATAGCATTTCTGAATATA-----")
        self.assertEqual(alignment.column_annotations["oryCun1.scaffold_156751"], "-------------------999999999999999999999999999-----")
        self.assertEqual(alignment.sequences[2].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[2].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[3].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[3].seq), 174210431)
        self.assertEqual(alignment.sequences[3].seq[16173169:16173169+25], "TCTGAGTGTGTATGAATATGAAGTA")
        self.assertEqual(alignment[3], "-----------------------TCTG---AGTGTGTATGAATATGAAGTA")
        self.assertEqual(alignment.sequences[3].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[3].annotations['leftCount'], 15)
        self.assertEqual(alignment.sequences[3].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[3].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[4].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[4].seq), 173908612)
        self.assertEqual(alignment.sequences[4].seq[16393517:16393517+25], "TCTGAGTATGTATGAATATGAAGTA")
        self.assertEqual(alignment[4], "-----------------------TCTG---AGTATGTATGAATATGAAGTA")
        self.assertEqual(alignment.column_annotations["panTro2.chr6"], "-----------------------9999---999999999999999999999")
        self.assertEqual(alignment.sequences[4].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[4].annotations['leftCount'], 15)
        self.assertEqual(alignment.sequences[4].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[4].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[5].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[5].seq), 170899992)
        self.assertEqual(alignment.sequences[5].seq[15874951:15874951+25], "TCTGAGTATGTATGAATATGAAGTA")
        self.assertEqual(alignment[5], "-----------------------TCTG---AGTATGTATGAATATGAAGTA")
        self.assertEqual(alignment.sequences[5].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[5].annotations['leftCount'], 15)
        self.assertEqual(alignment.sequences[5].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[5].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[6].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[6].seq), 133105)
        self.assertEqual(alignment.sequences[6].seq[10782:10782+25], "TCTGAGTATGTCTGAATATGAAGTG")
        self.assertEqual(alignment[6], "-----------------------TCTG---AGTATGTCTGAATATGAAGTG")
        self.assertEqual(alignment.sequences[6].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[6].annotations['leftCount'], 15)
        self.assertEqual(alignment.sequences[6].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[6].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[7].id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(alignment.sequences[7].seq), 359464)
        self.assertEqual(alignment.sequences[7].seq[188337:188337+25], "TCTAAGTGTGTCTGAACATGATGTG")
        self.assertEqual(alignment[7], "-----------------------TCTA---AGTGTGTCTGAACATGATGTG")
        self.assertEqual(alignment.column_annotations["otoGar1.scaffold_334.1-359464"], "-----------------------9999---999999999999999999999")
        self.assertEqual(alignment.sequences[7].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[7].annotations['leftCount'], 15)
        self.assertEqual(alignment.sequences[7].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[7].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[8].id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(alignment.sequences[8].seq), 498454)
        self.assertEqual(alignment.sequences[8].seq[174961:174961+25], "TCTCTGAGTGTCTGAACTTGAGGTA")
        self.assertEqual(alignment[8], "-----------------------TCTC---TGAGTGTCTGAACTTGAGGTA")
        self.assertEqual(alignment.column_annotations["tupBel1.scaffold_114895.1-498454"], "-----------------------9999---999999999999999999999")
        self.assertEqual(alignment.sequences[8].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[8].annotations['leftCount'], 623)
        self.assertEqual(alignment.sequences[8].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[8].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[9].id, "eriEur1.scaffold_266115")
        self.assertEqual(len(alignment.sequences[9].seq), 4589)
        self.assertEqual(alignment.sequences[9].seq[432:432+19], "TCTCAGTGTGTCTGACCAG")
        self.assertEqual(alignment[9], "-----------------------TCTC---AGTGTGTCTGACCAG------")
        self.assertEqual(alignment.column_annotations["eriEur1.scaffold_266115"], "-----------------------9999---999999999999999------")
        self.assertEqual(alignment.sequences[9].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[9].annotations['leftCount'], 9)
        self.assertEqual(alignment.sequences[9].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[9].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[10].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[10].seq), 125616256)
        self.assertEqual(alignment.sequences[10].seq[78071939:78071939+25], "TCTGAGTGTGTCTGAATGTGAGGTA")
        self.assertEqual(alignment[10], "-----------------------TCTG---AGTGTGTCTGAATGTGAGGTA")
        self.assertEqual(alignment.column_annotations["canFam2.chr1"], "-----------------------9999---999999999999999999999")
        self.assertEqual(alignment.sequences[10].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[10].annotations['leftCount'], 9)
        self.assertEqual(alignment.sequences[10].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[10].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[11].id, "felCat3.scaffold_205680")
        self.assertEqual(len(alignment.sequences[11].seq), 119354)
        self.assertEqual(alignment.sequences[11].seq[73459:73459+25], "TCTGAGTGTTCCTGAATGTGAGGTG")
        self.assertEqual(alignment[11], "-----------------------TCTG---AGTGTTCCTGAATGTGAGGTG")
        self.assertEqual(alignment.column_annotations["felCat3.scaffold_205680"], "-----------------------9999---999999999999999999999")
        self.assertEqual(alignment.sequences[11].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[11].annotations['leftCount'], 9)
        self.assertEqual(alignment.sequences[11].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[11].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[12].id, "dasNov1.scaffold_56749")
        self.assertEqual(len(alignment.sequences[12].seq), 10470)
        self.assertEqual(alignment.sequences[12].seq[3305:3305+47], "TTGTGTGCCTTTGCGGAGCATTTCTGAGCGTGTCTGAACTTGAGGTA")
        self.assertEqual(alignment[12], "T-TGTGTGCCTTTGCGGAGCATTTCTG---AGCGTGTCTGAACTTGAGGTA")
        self.assertEqual(alignment.column_annotations["dasNov1.scaffold_56749"], "7-3659557766555777595699547---965955937797755656587")
        self.assertEqual(alignment.sequences[12].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[12].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[12].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[12].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[13].id, "echTel1.scaffold_288249")
        self.assertEqual(len(alignment.sequences[13].seq), 100002)
        self.assertEqual(alignment.sequences[13].seq[95297:95297+48], "TGCCCAGGACTGCGCATGGTATTTCTTGGTGTGTCTGAAGGTGAGATA")
        self.assertEqual(alignment[13], "TGCCCAGGACTGCGCATGGTATTTCTT---GGTGTGTCTGAAGGTGAGATA")
        self.assertEqual(alignment.column_annotations["echTel1.scaffold_288249"], "999999999999999999999999999---999999999999999999999")
        self.assertEqual(alignment.sequences[13].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[13].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[13].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[13].annotations['rightCount'], 0)
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
                numpy.array([[ 3021180,  3021180,  3021180,  3021180,  3021184,  3021184,
                               3021184,  3021195,  3021195,  3021198,  3021199,  3021204],
                             [    9638,     9638,     9638,     9638,     9642,     9646,
                                  9649,     9660,     9661,     9664,     9665,     9670],
                             [     746,      746,      746,      746,      750,      754,
                                   757,      768,      769,      772,      773,      773],
                             [16173169, 16173169, 16173169, 16173169, 16173169, 16173173,
                              16173173, 16173184, 16173185, 16173188, 16173189, 16173194],
                             [16393517, 16393517, 16393517, 16393517, 16393517, 16393521,
                              16393521, 16393532, 16393533, 16393536, 16393537, 16393542],
                             [15874951, 15874951, 15874951, 15874951, 15874951, 15874955,
                              15874955, 15874966, 15874967, 15874970, 15874971, 15874976],
                             [   10782,    10782,    10782,    10782,    10782,    10786,
                                 10786,    10797,    10798,    10801,    10802,    10807],
                             [  188337,   188337,   188337,   188337,   188337,   188341,
                                188341,   188352,   188353,   188356,   188357,   188362],
                             [  174961,   174961,   174961,   174961,   174961,   174965,
                                174965,   174976,   174977,   174980,   174981,   174986],
                             [     432,      432,      432,      432,      432,      436,
                                   436,      447,      448,      451,      451,      451],
                             [78071939, 78071939, 78071939, 78071939, 78071939, 78071943,
                              78071943, 78071954, 78071955, 78071958, 78071959, 78071964],
                             [   73459,    73459,    73459,    73459,    73459,    73463,
                                 73463,    73474,    73475,    73478,    73479,    73484],
                             [    3305,     3306,     3306,     3323,     3327,     3331,
                                  3331,     3342,     3343,     3346,     3347,     3352],
                             [   95297,    95298,    95299,    95316,    95320,    95324,
                                 95324,    95335,    95336,    95339,    95340,    95345]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 309116)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3021204:3021204+71], "TGCACTGGTTTTCCTGCAGTGGTTCTCAGTAATAGGAAGACAACAGAATTTGAAGTATCCGGCTTTGGCCA")
        self.assertEqual(alignment[0], "TGCACTGGTTTTCC-TGCAGTGGTTCTCAGTAATAGGAAGACA-ACAGAATTTGAAGTATCCGGCTTTGGCCA")
        self.assertEqual(alignment.sequences[1].id, "echTel1.scaffold_288249")
        self.assertEqual(len(alignment.sequences[1].seq), 100002)
        self.assertEqual(alignment.sequences[1].seq[95345:95345+73], "GGCATTGGTTTTTAGAGAGAGAACCCACATAAGTAGGAAAACATTTTGAATTTATAGTAAATATTCTTGGCTA")
        self.assertEqual(alignment[1], "GGCATTGGTTTTTAGAGAGAGAACCCACATAAGTAGGAAAACATTTTGAATTTATAGTAAATATTCTTGGCTA")
        self.assertEqual(alignment.column_annotations["echTel1.scaffold_288249"], "9999999999999999999999999999999999999999999999999999999999999999999999999")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[2].id, "dasNov1.scaffold_56749")
        self.assertEqual(len(alignment.sequences[2].seq), 10470)
        self.assertEqual(alignment.sequences[2].seq[3352:3352+70], "TGCATTGATTTTCAGTGAGGGAAATCACACTAAGAGAGCATTTTGAATTTATAATAAATGCTCTTGGCTA")
        self.assertEqual(alignment[2], "TGCATTGATTTTCAGTGAGGGAAATCACAC---TAAGAGAGCATTTTGAATTTATAATAAATGCTCTTGGCTA")
        self.assertEqual(alignment.column_annotations["dasNov1.scaffold_56749"], "579663579997957596666783594865---8458696597898258979678997999677999997699")
        self.assertEqual(alignment.sequences[2].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[2].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[3].id, "felCat3.scaffold_205680")
        self.assertEqual(len(alignment.sequences[3].seq), 119354)
        self.assertEqual(alignment.sequences[3].seq[73484:73484+72], "TGCATTGATTTTCATGCATTGATTCACATGAATAGGAATGCATTTTGAAATTATAGCAAATGTTTTTGGCTG")
        self.assertEqual(alignment[3], "TGCATTGATTTTCA-TGCATTGATTCACATGAATAGGAATGCATTTTGAAATTATAGCAAATGTTTTTGGCTG")
        self.assertEqual(alignment.column_annotations["felCat3.scaffold_205680"], "99999999999999-9999999999999999999999999999999999999999999999999999999999")
        self.assertEqual(alignment.sequences[3].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[3].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[3].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[3].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[4].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[4].seq), 125616256)
        self.assertEqual(alignment.sequences[4].seq[78071964:78071964+71], "TGCATTAATTTTCATGCGTTGATTCACATGAATAGGAATACATTTGAATTTATAGGAAATGCTTTTGGCTG")
        self.assertEqual(alignment[4], "TGCATTAATTTTCA-TGCGTTGATTCACATGAATAGGAATACA-TTTGAATTTATAGGAAATGCTTTTGGCTG")
        self.assertEqual(alignment.column_annotations["canFam2.chr1"], "99999999999999-9999999999999999999999999999-99999999999999999999999999999")
        self.assertEqual(alignment.sequences[4].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[4].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[4].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[4].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[5].id, "eriEur1.scaffold_266115")
        self.assertEqual(len(alignment.sequences[5].seq), 4589)
        self.assertEqual(alignment.sequences[5].seq[451:451+70], "GCACTGATTCTCGAGGGTTGATTCCCAGTAAGAGGAAACTGCGTGAGTTTACAGTACATGGGCTTGGCTG")
        self.assertEqual(alignment[5], "-GCACTGATTCTCG-AGGGTTGATTCCCAGTAAGAGGAAACTG-CGTGAGTTTACAGTACATGGGCTTGGCTG")
        self.assertEqual(alignment.column_annotations["eriEur1.scaffold_266115"], "-9999999999999-9999999999999999999999999999-99999999999999999999999999799")
        self.assertEqual(alignment.sequences[5].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[5].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[5].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[5].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[6].id, "sorAra1.scaffold_2476")
        self.assertEqual(len(alignment.sequences[6].seq), 4997)
        self.assertEqual(alignment.sequences[6].seq[3312:3312+70], "TGCATTGATTTTTAAGTATTGATTCACATGAAAAGGAAAGCACTTCAATTTATGATAAATGATTTGGCCA")
        self.assertEqual(alignment[6], "TGCATTGATTTTTA-AGTATTGATTCACATGAAAAGGAAAGCA-CTTCAATTTATGATAAAT-GATTTGGCCA")
        self.assertEqual(alignment.column_annotations["sorAra1.scaffold_2476"], "89999999999989-9999999999999999999999999999-999999999799999999-9999999999")
        self.assertEqual(alignment.sequences[6].annotations['leftStatus'], 'N')
        self.assertEqual(alignment.sequences[6].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[6].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[6].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[7].id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(alignment.sequences[7].seq), 498454)
        self.assertEqual(alignment.sequences[7].seq[174986:174986+71], "TGCATTGAATTCCATGTGTTATGTCACATTAATAGGAAACGATTCTAATTGATAGGAAATGCTTCTGACCA")
        self.assertEqual(alignment[7], "TGCATTGAATTCCA-TGTGTTATGTCACATTAATAGGAAACGA-TTCTAATTGATAGGAAATGCTTCTGACCA")
        self.assertEqual(alignment.column_annotations["tupBel1.scaffold_114895.1-498454"], "99999899899999-9999999999999999999999999999-99999999999999999999999999997")
        self.assertEqual(alignment.sequences[7].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[7].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[7].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[7].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[8].id, "otoGar1.scaffold_334.1-359464")
        self.assertEqual(len(alignment.sequences[8].seq), 359464)
        self.assertEqual(alignment.sequences[8].seq[188362:188362+67], "TGGATTTTTACACATTGATTCACCTTAACGGGAAAACATGTGAATTTGTAGGAAAGGCATTTGGCCA")
        self.assertEqual(alignment[8], "TG----GATTTTTA-CACATTGATTCACCTTAACGGGAAAACA-TGTGAATTTGTAGGAAAGGCATTTGGCCA")
        self.assertEqual(alignment.column_annotations["otoGar1.scaffold_334.1-359464"], "88----99999999-9899999899943848888999999999-66359999569699999923799351281")
        self.assertEqual(alignment.sequences[8].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[8].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[8].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[8].annotations['rightCount'], 6280)
        self.assertEqual(alignment.sequences[9].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[9].seq), 133105)
        self.assertEqual(alignment.sequences[9].seq[10807:10807+71], "TGCATTGGTTTTTATGCCTTGATTCACATGAATAGGAAAACGTTTGAATTTATAGGAAATGGTTTTGGCCA")
        self.assertEqual(alignment[9], "TGCATTGGTTTTTA-TGCCTTGATTCACATGAATAGGAAAACG-TTTGAATTTATAGGAAATGGTTTTGGCCA")
        self.assertEqual(alignment.sequences[9].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[9].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[9].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[9].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[10].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[10].seq), 170899992)
        self.assertEqual(alignment.sequences[10].seq[15874976:15874976+71], "TGCTTTGATTTTTATGCATTGATTTACATGAATAGGAAAATGTTTGAATTTATAGTAAATGGTTTTGGCCA")
        self.assertEqual(alignment[10], "TGCTTTGATTTTTA-TGCATTGATTTACATGAATAGGAAAATG-TTTGAATTTATAGTAAATGGTTTTGGCCA")
        self.assertEqual(alignment.sequences[10].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[10].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[10].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[10].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[11].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[11].seq), 173908612)
        self.assertEqual(alignment.sequences[11].seq[16393542:16393542+71], "TGCTTTGATTTTTATGCATTGATTTACATGAATAGGAAAACATTTGAATTTATAGTAAATGGTTTTGGCCA")
        self.assertEqual(alignment[11], "TGCTTTGATTTTTA-TGCATTGATTTACATGAATAGGAAAACA-TTTGAATTTATAGTAAATGGTTTTGGCCA")
        self.assertEqual(alignment.column_annotations["panTro2.chr6"], "99999999999999-9999999999999999999999999999-99999999999999999999999999999")
        self.assertEqual(alignment.sequences[11].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[11].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[11].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[11].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[12].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[12].seq), 174210431)
        self.assertEqual(alignment.sequences[12].seq[16173194:16173194+71], "TGCTTTGATTTTTATGCATTGATTTACATGGATAGGAAAACGTTTGAATTTATAGTAAATGGTTTTGGCCA")
        self.assertEqual(alignment[12], "TGCTTTGATTTTTA-TGCATTGATTTACATGGATAGGAAAACG-TTTGAATTTATAGTAAATGGTTTTGGCCA")
        self.assertEqual(alignment.sequences[12].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[12].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[12].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[12].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[13].id, "oryCun1.scaffold_156751")
        self.assertEqual(len(alignment.sequences[13].seq), 4726)
        self.assertEqual(alignment.sequences[13].seq[773:773+65], "TCTTTACATACTTTGATTTACATCGACGGGAAAACATTTGCATTTATAGCAATGGCTTTTGGTCA")
        self.assertEqual(alignment[13], "------TCTTTACA-TACTTTGATTTACATCGACGGGAAAACA-TTTGCATTTATAGCAATGGCTTTTGGTCA")
        self.assertEqual(alignment.column_annotations["oryCun1.scaffold_156751"], "------99999994-9999999999999899999999999999-99999999999999899999999999999")
        self.assertEqual(alignment.sequences[13].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[13].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[13].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[13].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[14].id, "cavPor2.scaffold_216473")
        self.assertEqual(len(alignment.sequences[14].seq), 10026)
        self.assertEqual(alignment.sequences[14].seq[9670:9670+71], "GGCATTGATTTCCATGCATTAATTCACATTAATAGGAAAACAGTCCAATTTATAGTAAATGGTTTTGGCCA")
        self.assertEqual(alignment[14], "GGCATTGATTTCCA-TGCATTAATTCACATTAATAGGAAAACA-GTCCAATTTATAGTAAATGGTTTTGGCCA")
        self.assertEqual(alignment.column_annotations["cavPor2.scaffold_216473"], "77888768786695-3886758796446556688656665478-68687676669688688666687574686")
        self.assertEqual(alignment.sequences[14].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[14].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[14].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[14].annotations['rightCount'], 0)
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
                numpy.array([[ 3021204,  3021205,  3021206,  3021210,  3021218,  3021218,
                               3021233,  3021236,  3021246,  3021246,  3021264,  3021265,
                               3021275],
                             [   95345,    95346,    95347,    95351,    95359,    95360,
                                 95375,    95378,    95388,    95389,    95407,    95408,
                                 95418],
                             [    3352,     3353,     3354,     3358,     3366,     3367,
                                  3382,     3382,     3392,     3393,     3411,     3412,
                                  3422],
                             [   73484,    73485,    73486,    73490,    73498,    73498,
                                 73513,    73516,    73526,    73527,    73545,    73546,
                                 73556],
                             [78071964, 78071965, 78071966, 78071970, 78071978, 78071978,
                              78071993, 78071996, 78072006, 78072006, 78072024, 78072025,
                              78072035],
                             [     451,      451,      452,      456,      464,      464,
                                   479,      482,      492,      492,      510,      511,
                                   521],
                             [    3312,     3313,     3314,     3318,     3326,     3326,
                                  3341,     3344,     3354,     3354,     3372,     3372,
                                  3382],
                             [  174986,   174987,   174988,   174992,   175000,   175000,
                                175015,   175018,   175028,   175028,   175046,   175047,
                                175057],
                             [  188362,   188363,   188364,   188364,   188372,   188372,
                                188387,   188390,   188400,   188400,   188418,   188419,
                                188429],
                             [   10807,    10808,    10809,    10813,    10821,    10821,
                                 10836,    10839,    10849,    10849,    10867,    10868,
                                 10878],
                             [15874976, 15874977, 15874978, 15874982, 15874990, 15874990,
                              15875005, 15875008, 15875018, 15875018, 15875036, 15875037,
                              15875047],
                             [16393542, 16393543, 16393544, 16393548, 16393556, 16393556,
                              16393571, 16393574, 16393584, 16393584, 16393602, 16393603,
                              16393613],
                             [16173194, 16173195, 16173196, 16173200, 16173208, 16173208,
                              16173223, 16173226, 16173236, 16173236, 16173254, 16173255,
                              16173265],
                             [     773,      773,      773,      773,      781,      781,
                                   796,      799,      809,      809,      827,      828,
                                   838],
                             [    9670,     9671,     9672,     9676,     9684,     9684,
                                  9699,     9702,     9712,     9712,     9730,     9731,
                                  9741]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 891219)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3021275:3021275+146], "CTTTCTCTGACATTACTGTAACTGAAGTAGCTCAGAAGCACAAAAGGTCACATCATGCATCCATGCAGAATCCACTGAAGCTGTTTGGAAAGGCCACGTGTCTTCCCAGAAGGCCAGTTACACCATCATTTCCTTCCATGTTTCAG")
        self.assertEqual(alignment[0], "CTTTCTCTGACATTACTGTAACTGAAGTAGCTC-AGAAGCACAAAAGGTCACATCATGCATCCATGCAGAATCCACTGAAGCTGTTTGGAAAGGC-----------------------CACGTGTCTTCCCAGAAGGCCAGTTACACCATCATTTCCTTCCATGTTTCAG")
        self.assertEqual(alignment.sequences[1].id, "cavPor2.scaffold_216473")
        self.assertEqual(len(alignment.sequences[1].seq), 10026)
        self.assertEqual(alignment.sequences[1].seq[9741:9741+169], "CCTTTGCTGATGTTACTGTAAGCGAACTAGCTGAAAAGCACAAAAGGTCATGTCAGGCATTTATGCAAAATTTACTAGTTCCGTTCATAAAGGCACCTCGCATTTTATTACTTAGTTCATGTGTCTCGCAGGAAATCCAGTCACGCCATCATTTCCTTCCATATTTCAA")
        self.assertEqual(alignment[1], "CCTTTGCTGATGTTACTGTAAGCGAACTAGCTG-AAAAGCACAAAAGGTCATGTCAGGCATTTATGCAAAATTTACTAGTTCCGTTCATAAAGGCACCTCGCATTTTATTACTTAGTTCATGTGTCTCGCAGGAAATCCAGTCACGCCATCATTTCCTTCCATATTTCAA")
        self.assertEqual(alignment.column_annotations["cavPor2.scaffold_216473"], "668536545547664539866666666666668-8787886799997788875886655536666688786768476688854689668876558688874778656675665668746666788474786888875667457786687996867886548798845555")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[2].id, "oryCun1.scaffold_156751")
        self.assertEqual(len(alignment.sequences[2].seq), 4726)
        self.assertEqual(alignment.sequences[2].seq[838:838+158], "TCGTCTCTGATGTTACTCTAGCTGAAAGCGCTCCAAAGCACAAAAGGCCATGGCAGGCATTTATGTCAAATGTACTAGCACTGTTTATAACAGCACCTTGTGTTGTATTACTTTGTTCAGGAAATCCAGTCACGCCGCGATTTCCTTCCACATTTCAA")
        self.assertEqual(alignment[2], "TCGTCTCTGATGTTACTCTAGCTGAAAGCGCTC-CAAAGCACAAAAGGCCATGGCAGGCATTTATGTCAAATGTACTAGCACTGTTTATAACAGCACCTTGTGTTGTATTACTTTGTTCA-----------GGAAATCCAGTCACGCCGCGATTTCCTTCCACATTTCAA")
        self.assertEqual(alignment.column_annotations["oryCun1.scaffold_156751"], "999999879999999999999999999997999-99999799999999898999999988999999988999999989996999999999988897996999999999999979999987-----------899998987999799998998999999999777999899")
        self.assertEqual(alignment.sequences[2].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[2].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[3].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[3].seq), 174210431)
        self.assertEqual(alignment.sequences[3].seq[16173265:16173265+169], "CCTTCTCTGATGTTACTGTAACTGAAATAGCTCAAAAGCACAAAAGGTCATGGGAGGCATTCATGTAAAACTTACTAGTGCTGTTTATAAGGGTACCTTGTGTTTTATTACTTTGTTAATGTATCTTTCAGGAAATTCAGTCACGCCATCATTTCCTTCCATATTTCAT")
        self.assertEqual(alignment[3], "CCTTCTCTGATGTTACTGTAACTGAAATAGCTC-AAAAGCACAAAAGGTCATGGGAGGCATTCATGTAAAACTTACTAGTGCTGTTTATAAGGGTACCTTGTGTTTTATTACTTTGTTAATGTATCTTTCAGGAAATTCAGTCACGCCATCATTTCCTTCCATATTTCAT")
        self.assertEqual(alignment.sequences[3].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[3].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[3].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[3].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[4].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[4].seq), 173908612)
        self.assertEqual(alignment.sequences[4].seq[16393613:16393613+169], "CCTTCTCTGATGTTACTGTAACTGAAATAGCTCAAAAGCACAAAAGGTCATGGGAGGCATTTATGTAAAACTTACTAGTGCTGTTTATAAGGGTACCTTGTGTTTTATTACTTTGTTAATGTACCTTTCAGGAAATCCAGTCATGCCATCATTTCCTTCCATATTTCAT")
        self.assertEqual(alignment[4], "CCTTCTCTGATGTTACTGTAACTGAAATAGCTC-AAAAGCACAAAAGGTCATGGGAGGCATTTATGTAAAACTTACTAGTGCTGTTTATAAGGGTACCTTGTGTTTTATTACTTTGTTAATGTACCTTTCAGGAAATCCAGTCATGCCATCATTTCCTTCCATATTTCAT")
        self.assertEqual(alignment.column_annotations["panTro2.chr6"], "999999999999999999999999999999999-9999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999")
        self.assertEqual(alignment.sequences[4].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[4].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[4].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[4].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[5].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[5].seq), 170899992)
        self.assertEqual(alignment.sequences[5].seq[15875047:15875047+169], "CCTTCTCTGATGTTACTGTAACTGAAATAGCTCAAAAGCACAAAAGGTCATGGGAGGCATTTATGTAAAACTTACTAGTGCTGTTTATAAGGGTATCTTGTGTTTTATTACTTTGTTAATGTACCTTTCAGGAAATCCAGTCATGCCATCATTTCCTTCCATATTTCAT")
        self.assertEqual(alignment[5], "CCTTCTCTGATGTTACTGTAACTGAAATAGCTC-AAAAGCACAAAAGGTCATGGGAGGCATTTATGTAAAACTTACTAGTGCTGTTTATAAGGGTATCTTGTGTTTTATTACTTTGTTAATGTACCTTTCAGGAAATCCAGTCATGCCATCATTTCCTTCCATATTTCAT")
        self.assertEqual(alignment.sequences[5].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[5].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[5].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[5].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[6].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[6].seq), 133105)
        self.assertEqual(alignment.sequences[6].seq[10878:10878+169], "CCTTCTCTGATATTACTGTAACTGAAATAGCTCAAAAGCACAAAAGGTCATGGCAGGCATTTATGTAAAACTTACTAGTGCTGTTTATAAAGGCACCTTGTGTTTTATTACTTTGTTAATGTACCTTTCAGGAAATCCAGTCATGCCATCATTTCCTTCCATATTTCAT")
        self.assertEqual(alignment[6], "CCTTCTCTGATATTACTGTAACTGAAATAGCTC-AAAAGCACAAAAGGTCATGGCAGGCATTTATGTAAAACTTACTAGTGCTGTTTATAAAGGCACCTTGTGTTTTATTACTTTGTTAATGTACCTTTCAGGAAATCCAGTCATGCCATCATTTCCTTCCATATTTCAT")
        self.assertEqual(alignment.sequences[6].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[6].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[6].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[6].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[7].id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(alignment.sequences[7].seq), 498454)
        self.assertEqual(alignment.sequences[7].seq[175057:175057+156], "TCCTCTCAGATGTTCCTGTAACCAAAGTAGCTCAGAAGCACAAAAGGTCATGCCAGGCATTGATGTAAAGTTTCTGAAGGCATCTCGTGTTTTATTACTTTGTTCATGTGCCTTACAGGAAATCCAGTCGGGCCATCATTTCCTTCCATATCTCAG")
        self.assertEqual(alignment[7], "TCCTCTCAGATGTTCCTGTAACCAAAGTAGCTC-AGAAGCACAAAAGGTCATGCCAGGCATTGATGTAAA-------------GTTTCTGAAGGCATCTCGTGTTTTATTACTTTGTTCATGTGCCTTACAGGAAATCCAGTCGGGCCATCATTTCCTTCCATATCTCAG")
        self.assertEqual(alignment.column_annotations["tupBel1.scaffold_114895.1-498454"], "999999999999999999999999999999999-999999999999999999999999999999999999-------------999999999999999999999999999999999999999999999999999999999999999999999999999999999999999")
        self.assertEqual(alignment.sequences[7].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[7].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[7].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[7].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[8].id, "sorAra1.scaffold_2476")
        self.assertEqual(len(alignment.sequences[8].seq), 4997)
        self.assertEqual(alignment.sequences[8].seq[3382:3382+170], "CCTTCTCTGATGTTACTGTAACTGAGGTAGCTCAAAAAGCACAAAAGGTCATGTCAGGCATTGATGTAAAATATACCAGTGCTGTTTATAAAGACCCCTTGCATATAATGTTGTTGATCATGAACTTTTCCGGAAATCCTGTAATCTCATCATTTCCTCCCATATTCCAG")
        self.assertEqual(alignment[8], "CCTTCTCTGATGTTACTGTAACTGAGGTAGCTCAAAAAGCACAAAAGGTCATGTCAGGCATTGATGTAAAATATACCAGTGCTGTTTATAAAGACCCCTTGCATATAATGTTGTTGATCATGAACTTTTCCGGAAATCCTGTAATCTCATCATTTCCTCCCATATTCCAG")
        self.assertEqual(alignment.column_annotations["sorAra1.scaffold_2476"], "99999999999999997999999999999999999999999999999999999999999999999999999999999979999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999")
        self.assertEqual(alignment.sequences[8].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[8].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[8].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[8].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[9].id, "eriEur1.scaffold_266115")
        self.assertEqual(len(alignment.sequences[9].seq), 4589)
        self.assertEqual(alignment.sequences[9].seq[521:521+156], "CCTTCTCTGATCTGACTGTAACCGAAGTCTCTCGAAAGCACAAAAGCTCACGGCAGGCCTTTCTGTAAAATACAGCGGCGCTGCTTCCAAAGGCaccttgcagatgtctctcacgcATTTGGCAGGAGTCCCTGTCACTCGGCCAGTTCCTTCCTG")
        self.assertEqual(alignment[9], "CCTTCTCTGATCTGACTGTAACCGAAGTCTCTC-GAAAGCACAAAAGCTCACGGCAGGCCTTTCTGTAAAATACAGCGGCGCTGCTTCCAAAGGCaccttgca------gatgtctctcacgcATTTGGCAGGAGTCCCTGTCACTCGGCCAGTTCC-------TTCCTG")
        self.assertEqual(alignment.column_annotations["eriEur1.scaffold_266115"], "999999999999999999999999999999999-999999999999999999999999799999999999999999999998987999899999999999999------999999999999999999999999999999999999999999999998-------999989")
        self.assertEqual(alignment.sequences[9].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[9].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[9].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[9].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[10].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[10].seq), 125616256)
        self.assertEqual(alignment.sequences[10].seq[78072035:78072035+168], "CTGTCTCTGATGTTACTGTAACTGAAGTAGCTCAAAAGCACAAAAGGTCATGTCAGGCTTTATGTAAAATTTACTAGAGCTGTTTATAAAGGCACCTCGTGTAATATTGCTTTGTTCATGTACTTTACAGGAAATCCTGTCATCCTGTTACTTCCTTCCATATTTCAA")
        self.assertEqual(alignment[10], "CTGTCTCTGATGTTACTGTAACTGAAGTAGCTC-AAAAGCACAAAAGGTCATGTCAGGC-TTTATGTAAAATTTACTAGAGCTGTTTATAAAGGCACCTCGTGTAATATTGCTTTGTTCATGTACTTTACAGGAAATCCTGTCATCCTGTTACTTCCTTCCATATTTCAA")
        self.assertEqual(alignment.column_annotations["canFam2.chr1"], "999999999999999999999999999999999-9999999999999999999999999-99999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999")
        self.assertEqual(alignment.sequences[10].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[10].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[10].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[10].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[11].id, "felCat3.scaffold_205680")
        self.assertEqual(len(alignment.sequences[11].seq), 119354)
        self.assertEqual(alignment.sequences[11].seq[73556:73556+169], "CTGTCTCTGATGTTACTGTAACTGAAATAGTTCAAAAGCACAAAAGGTGATGTCAGGCATTTATGTAAAATTTACTAGTGCTGTTTATAAAGGTACCTTGTGTATTATTGCTTTGTTCATGCACTTTATAGGAAATCCTGTCATCCTGTCATTTCCTTCCATAGTTCAG")
        self.assertEqual(alignment[11], "CTGTCTCTGATGTTACTGTAACTGAAATAGTTC-AAAAGCACAAAAGGTGATGTCAGGCATTTATGTAAAATTTACTAGTGCTGTTTATAAAGGTACCTTGTGTATTATTGCTTTGTTCATGCACTTTATAGGAAATCCTGTCATCCTGTCATTTCCTTCCATAGTTCAG")
        self.assertEqual(alignment.column_annotations["felCat3.scaffold_205680"], "999999999999999999999999999999999-9999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999")
        self.assertEqual(alignment.sequences[11].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[11].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[11].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[11].annotations['rightCount'], 34165)
        self.assertEqual(alignment.sequences[12].id, "dasNov1.scaffold_56749")
        self.assertEqual(len(alignment.sequences[12].seq), 10470)
        self.assertEqual(alignment.sequences[12].seq[3422:3422+168], "CCTTCTCTGATGTTACTGTAACTGAAGAGGCTCAAAAGCACAAAAGGTCATGTCGGGCATTGATGTAAAATTTACTAGTGCTGTTTATAAAGGACCTTGTGTTCTATTACTTTTTCTATCCCCCTTACAGGAAATCCTGTCATACCATCATTTCCTTCCATAAGTCAC")
        self.assertEqual(alignment[12], "CCTTCTCTGATGTTACTGTAACTGAAGAGGCTC-AAAAGCACAAAAGGTCATGTCGGGCATTGATGTAAAATTTACTAGTGCTG-TTTATAAAGGACCTTGTGTTCTATTACTTTTTCTATCCCCCTTACAGGAAATCCTGTCATACCATCATTTCCTTCCATAAGTCAC")
        self.assertEqual(alignment.column_annotations["dasNov1.scaffold_56749"], "989999999999999937699999999999999-79999999999999999999999999999999999999999999999999-9999999999999999999999999999999999999999999999999999999999999999999999999999999999999")
        self.assertEqual(alignment.sequences[12].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[12].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[12].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[12].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[13].id, "echTel1.scaffold_288249")
        self.assertEqual(len(alignment.sequences[13].seq), 100002)
        self.assertEqual(alignment.sequences[13].seq[95418:95418+169], "CCATCTCTGATATTACTGTAATGGAATTACTTTGAAAGCACAAAAGGTCAGGACAAGCGTGTATGTGAAATTTCCTAGAGCTGTTTTCCTGCACACCTTGGATTTTATTGCTTAGTTCATTTGCTTTCCCAGAAATCCCGCCATGCCATCATTTCCTTCCACATCTCAG")
        self.assertEqual(alignment[13], "CCATCTCTGATATTACTGTAATGGAATTACTTT-GAAAGCACAAAAGGTCAGGACAAGCGTGTATGTGAAATTTCCTAGAGCTGTTTTCCTGCACACCTTGGATTTTATTGCTTAGTTCATTTGCTTTCCCAGAAATCCCGCCATGCCATCATTTCCTTCCACATCTCAG")
        self.assertEqual(alignment.column_annotations["echTel1.scaffold_288249"], "999999999999999999999999999999999-9999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999")
        self.assertEqual(alignment.sequences[13].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[13].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[13].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[13].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[14].id, "ornAna1.chr2")
        self.assertEqual(len(alignment.sequences[14].seq), 54797317)
        self.assertEqual(alignment.sequences[14].seq[14756885:14756885+169], "CCACCTCTGATATCACTGTAACTCAGCTACCCAGAAGATAGGAAAGGTCATATCAGCCAATCATGTGCAACTTAATAGCCCTGTTCATAAAGCTGCTTTGTGTTTTATGAGGCTATTCATTCACCTTAGGGGAAATCCTGTCAAGCCATCATTTCCTTCCATATCTTCA")
        self.assertEqual(alignment[14], "CCACCTCTGATATCACTGTAACTCAGCTACCCA-GAAGATAGGAAAGGTCATATCAGCCAATCATGTGCAACTTAATAGCCCTGTTCATAAAGCTGCTTTGTGTTTTATGAGGCTATTCATTCACCTTAGGGGAAATCCTGTCAAGCCATCATTTCCTTCCATATCTTCA")
        self.assertEqual(alignment.sequences[14].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[14].annotations['leftCount'], 5690)
        self.assertEqual(alignment.sequences[14].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[14].annotations['rightCount'], 0)
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
                numpy.array([[ 3021275,  3021308,  3021308,  3021333,  3021334,  3021344,
                               3021357,  3021358,  3021359,  3021369,  3021369,  3021369,
                               3021369,  3021371,  3021382,  3021408,  3021415,  3021421],
                             [    9741,     9774,     9774,     9799,     9800,     9810,
                                  9823,     9824,     9825,     9835,     9843,     9849,
                                  9858,     9860,     9871,     9897,     9904,     9910],
                             [     838,      871,      871,      896,      897,      907,
                                   920,      921,      922,      932,      940,      946,
                                   955,      957,      957,      983,      990,      996],
                             [16173265, 16173298, 16173298, 16173323, 16173324, 16173334,
                              16173347, 16173348, 16173349, 16173359, 16173367, 16173373,
                              16173382, 16173384, 16173395, 16173421, 16173428, 16173434],
                             [16393613, 16393646, 16393646, 16393671, 16393672, 16393682,
                              16393695, 16393696, 16393697, 16393707, 16393715, 16393721,
                              16393730, 16393732, 16393743, 16393769, 16393776, 16393782],
                             [15875047, 15875080, 15875080, 15875105, 15875106, 15875116,
                              15875129, 15875130, 15875131, 15875141, 15875149, 15875155,
                              15875164, 15875166, 15875177, 15875203, 15875210, 15875216],
                             [   10878,    10911,    10911,    10936,    10937,    10947,
                                 10960,    10961,    10962,    10972,    10980,    10986,
                                 10995,    10997,    11008,    11034,    11041,    11047],
                             [  175057,   175090,   175090,   175115,   175116,   175126,
                                175126,   175127,   175128,   175138,   175146,   175152,
                                175161,   175163,   175174,   175200,   175207,   175213],
                             [    3382,     3415,     3416,     3441,     3442,     3452,
                                  3465,     3466,     3467,     3477,     3485,     3491,
                                  3500,     3502,     3513,     3539,     3546,     3552],
                             [     521,      554,      554,      579,      580,      590,
                                   603,      604,      605,      615,      623,      623,
                                   632,      634,      645,      671,      671,      677],
                             [78072035, 78072068, 78072068, 78072093, 78072093, 78072103,
                              78072116, 78072117, 78072118, 78072128, 78072136, 78072142,
                              78072151, 78072153, 78072164, 78072190, 78072197, 78072203],
                             [   73556,    73589,    73589,    73614,    73615,    73625,
                                 73638,    73639,    73640,    73650,    73658,    73664,
                                 73673,    73675,    73686,    73712,    73719,    73725],
                             [    3422,     3455,     3455,     3480,     3481,     3491,
                                  3504,     3505,     3505,     3515,     3523,     3529,
                                  3538,     3540,     3551,     3577,     3584,     3590],
                             [   95418,    95451,    95451,    95476,    95477,    95487,
                                 95500,    95501,    95502,    95512,    95520,    95526,
                                 95535,    95537,    95548,    95574,    95581,    95587],
                             [14756885, 14756918, 14756918, 14756943, 14756944, 14756954,
                              14756967, 14756968, 14756969, 14756979, 14756987, 14756993,
                              14757002, 14757004, 14757015, 14757041, 14757048, 14757054]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 30254)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3021421:3021421+44], "ACGTTGCTCATTGTAATTGAAGCATTTATTACCAATGCCTTCCC")
        self.assertEqual(alignment[0], "-ACGTTGCTCATTGT-----AATTGAAGCATTTATTACCAA--------TG--------------CCTTCCC")
        self.assertEqual(alignment.sequences[1].id, "cavPor2.scaffold_216473")
        self.assertEqual(len(alignment.sequences[1].seq), 10026)
        self.assertEqual(alignment.sequences[1].seq[9910:9910+47], "CTGTTCCTTGTTACAAATAAAAGTATGTTTTGCTGACGGCCCCTCTG")
        self.assertEqual(alignment[1], "-CTGTTCCTTGTTACA----AATAAAAGTATGTTTTGCTGA--------CG-------G-----CCCCTCTG")
        self.assertEqual(alignment.column_annotations["cavPor2.scaffold_216473"], "-476754645667788----868678677444665235455--------76-------4-----36669356")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[2].id, "oryCun1.scaffold_156751")
        self.assertEqual(len(alignment.sequences[2].seq), 4726)
        self.assertEqual(alignment.sequences[2].seq[996:996+37], "TGGTTCCTCATTATGGATTAAAGCATTCACTGCTAAT")
        self.assertEqual(alignment[2], "-TGGTTCCTCATTATG----GATTAAAGCATTCACTGCTAA--------T----------------------")
        self.assertEqual(alignment.column_annotations["oryCun1.scaffold_156751"], "-999999999899999----887999999999999999999--------9----------------------")
        self.assertEqual(alignment.sequences[2].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[2].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[2].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[2].annotations['rightCount'], 2345)
        self.assertEqual(alignment.sequences[3].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[3].seq), 174210431)
        self.assertEqual(alignment.sequences[3].seq[16173434:16173434+45], "CTGTTGCTCATTATAAATAAAAGTGTTTATTGCTAATGCTTTCTC")
        self.assertEqual(alignment[3], "-CTGTTGCTCATTATA----AATAAAAGTGTTTATTGCTAA--------T--------------GCTTTCTC")
        self.assertEqual(alignment.sequences[3].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[3].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[3].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[3].annotations['rightCount'], 4)
        self.assertEqual(alignment.sequences[4].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[4].seq), 173908612)
        self.assertEqual(alignment.sequences[4].seq[16393782:16393782+45], "CTGTTCCTCATTATAAATAAAAGTGTTTATTGCTAACGCTTTCTC")
        self.assertEqual(alignment[4], "-CTGTTCCTCATTATA----AATAAAAGTGTTTATTGCTAA--------C--------------GCTTTCTC")
        self.assertEqual(alignment.column_annotations["panTro2.chr6"], "-999999999999999----999999999999999999999--------9--------------99999999")
        self.assertEqual(alignment.sequences[4].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[4].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[4].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[4].annotations['rightCount'], 4)
        self.assertEqual(alignment.sequences[5].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[5].seq), 170899992)
        self.assertEqual(alignment.sequences[5].seq[15875216:15875216+45], "CTGTTCCTCATTATAAATAAAAGTGTTTATTGCTAACGCTTTCTC")
        self.assertEqual(alignment[5], "-CTGTTCCTCATTATA----AATAAAAGTGTTTATTGCTAA--------C--------------GCTTTCTC")
        self.assertEqual(alignment.sequences[5].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[5].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[5].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[5].annotations['rightCount'], 4)
        self.assertEqual(alignment.sequences[6].id, "calJac1.Contig6394")
        self.assertEqual(len(alignment.sequences[6].seq), 133105)
        self.assertEqual(alignment.sequences[6].seq[11047:11047+43], "CTGTTTCTCATTATAAATAAGTGTTTATTGCTAACGCTTTCTC")
        self.assertEqual(alignment[6], "-CTGTTTCTCATTATA----AAT--AAGTGTTTATTGCTAA--------C--------------GCTTTCTC")
        self.assertEqual(alignment.sequences[6].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[6].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[6].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[6].annotations['rightCount'], 701)
        self.assertEqual(alignment.sequences[7].id, "tupBel1.scaffold_114895.1-498454")
        self.assertEqual(len(alignment.sequences[7].seq), 498454)
        self.assertEqual(alignment.sequences[7].seq[175213:175213+35], "CTGCTCCTCATGATAAATATAGCGTTTGTTGCTAA")
        self.assertEqual(alignment[7], "-CTGCTCCTCATGATA----AAT-ATAGCGTTTGTTGCTAA-------------------------------")
        self.assertEqual(alignment.column_annotations["tupBel1.scaffold_114895.1-498454"], "-999999999999999----999-99999999999999999-------------------------------")
        self.assertEqual(alignment.sequences[7].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[7].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[7].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[7].annotations['rightCount'], 10695)
        self.assertEqual(alignment.sequences[8].id, "sorAra1.scaffold_2476")
        self.assertEqual(len(alignment.sequences[8].seq), 4997)
        self.assertEqual(alignment.sequences[8].seq[3552:3552+48], "CTGTTCCTTCTTAGATTaaaaaaaGAGCGTTTATGACCAGCAGGTTCC")
        self.assertEqual(alignment[8], "-CTGTTCCTTCTTAGATTaaaaaaaGAGCGTTTATGACCAG--------CA-------GGTTCC--------")
        self.assertEqual(alignment.column_annotations["sorAra1.scaffold_2476"], "-9999999999999999999999999999999997669755--------99-------999999--------")
        self.assertEqual(alignment.sequences[8].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[8].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[8].annotations['rightStatus'], 'N')
        self.assertEqual(alignment.sequences[8].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[9].id, "eriEur1.scaffold_266115")
        self.assertEqual(len(alignment.sequences[9].seq), 4589)
        self.assertEqual(alignment.sequences[9].seq[677:677+35], "CTCCTTGTTCTTGAAGCCGTTTTTTATTGTCAGTG")
        self.assertEqual(alignment[9], "-CTCCTTGTTCTTGAAGC----------CGTTTTTTATTGT--------CA-------GTG-----------")
        self.assertEqual(alignment.column_annotations["eriEur1.scaffold_266115"], "-98959997997999999----------9999999999999--------89-------999-----------")
        self.assertEqual(alignment.sequences[9].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[9].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[9].annotations['rightStatus'], 'N')
        self.assertEqual(alignment.sequences[9].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[10].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[10].seq), 125616256)
        self.assertEqual(alignment.sequences[10].seq[78072203:78072203+40], "CTGTTTCTCATTATAAATAAAAGCATTTATTGCTAATGTC")
        self.assertEqual(alignment[10], "-CTGTTTCTCATTATA----AATAAAAGCATTTATTGCTAA--------TG-------T------------C")
        self.assertEqual(alignment.column_annotations["canFam2.chr1"], "-999999999999999----999999999999999999999--------99-------9------------9")
        self.assertEqual(alignment.sequences[10].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[10].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[10].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[10].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[11].id, "dasNov1.scaffold_56749")
        self.assertEqual(len(alignment.sequences[11].seq), 10470)
        self.assertEqual(alignment.sequences[11].seq[3590:3590+44], "CTGTTCCTCTTTATAATAAAAGCATTTACTGCTAACGGCTCCTC")
        self.assertEqual(alignment[11], "-CTGTTCCTCTTTAT-----AATAAAAGCATTTACTGCTAA--------CGGCTCCTC--------------")
        self.assertEqual(alignment.column_annotations["dasNov1.scaffold_56749"], "-99999999999999-----999999999999999999999--------999999999--------------")
        self.assertEqual(alignment.sequences[11].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[11].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[11].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[11].annotations['rightCount'], 904)
        self.assertEqual(alignment.sequences[12].id, "echTel1.scaffold_288249")
        self.assertEqual(len(alignment.sequences[12].seq), 100002)
        self.assertEqual(alignment.sequences[12].seq[95587:95587+38], "CTGATCCTCATTATAAATAAAAGTGTTTGTTACTAATG")
        self.assertEqual(alignment[12], "-CTGATCCTCATTATA----AATAAAAGTGTTTGTTACTAA--------TG---------------------")
        self.assertEqual(alignment.column_annotations["echTel1.scaffold_288249"], "-999999999999999----999999999999999999999--------99---------------------")
        self.assertEqual(alignment.sequences[12].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[12].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[12].annotations['rightStatus'], 'N')
        self.assertEqual(alignment.sequences[12].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[13].id, "ornAna1.chr2")
        self.assertEqual(len(alignment.sequences[13].seq), 54797317)
        self.assertEqual(alignment.sequences[13].seq[14757054:14757054+45], "ACTGTCCTCATTGTACAGGGACTCTTTTATTACTAAAGATCCCCC")
        self.assertEqual(alignment[13], "ACTG-TCCTCATTGTA----CAGGGACTCTTTTATTACTAAAGATCCCCC----------------------")
        self.assertEqual(alignment.sequences[13].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[13].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[13].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[13].annotations['rightCount'], 0)
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
                numpy.array([[ 3021421,  3021421,  3021424,  3021425,  3021435,  3021435,
                               3021435,  3021435,  3021438,  3021439,  3021440,  3021443,
                               3021456,  3021456,  3021457,  3021458,  3021458,  3021458,
                               3021458,  3021458,  3021458,  3021464,  3021465],
                             [    9910,     9910,     9913,     9914,     9924,     9925,
                                  9925,     9925,     9928,     9929,     9930,     9933,
                                  9946,     9946,     9947,     9948,     9948,     9949,
                                  9949,     9949,     9950,     9956,     9957],
                             [     996,      996,      999,     1000,     1010,     1011,
                                  1011,     1011,     1014,     1015,     1016,     1019,
                                  1032,     1032,     1033,     1033,     1033,     1033,
                                  1033,     1033,     1033,     1033,     1033],
                             [16173434, 16173434, 16173437, 16173438, 16173448, 16173449,
                              16173449, 16173449, 16173452, 16173453, 16173454, 16173457,
                              16173470, 16173470, 16173471, 16173471, 16173471, 16173471,
                              16173471, 16173471, 16173472, 16173478, 16173479],
                             [16393782, 16393782, 16393785, 16393786, 16393796, 16393797,
                              16393797, 16393797, 16393800, 16393801, 16393802, 16393805,
                              16393818, 16393818, 16393819, 16393819, 16393819, 16393819,
                              16393819, 16393819, 16393820, 16393826, 16393827],
                             [15875216, 15875216, 15875219, 15875220, 15875230, 15875231,
                              15875231, 15875231, 15875234, 15875235, 15875236, 15875239,
                              15875252, 15875252, 15875253, 15875253, 15875253, 15875253,
                              15875253, 15875253, 15875254, 15875260, 15875261],
                             [   11047,    11047,    11050,    11051,    11061,    11062,
                                 11062,    11062,    11065,    11065,    11065,    11068,
                                 11081,    11081,    11082,    11082,    11082,    11082,
                                 11082,    11082,    11083,    11089,    11090],
                             [  175213,   175213,   175216,   175217,   175227,   175228,
                                175228,   175228,   175231,   175231,   175232,   175235,
                                175248,   175248,   175248,   175248,   175248,   175248,
                                175248,   175248,   175248,   175248,   175248],
                             [    3552,     3552,     3555,     3556,     3566,     3567,
                                  3569,     3571,     3574,     3575,     3576,     3579,
                                  3592,     3592,     3593,     3594,     3594,     3595,
                                  3597,     3600,     3600,     3600,     3600],
                             [     677,      677,      680,      681,      691,      692,
                                   694,      694,      694,      694,      694,      694,
                                   707,      707,      708,      709,      709,      710,
                                   712,      712,      712,      712,      712],
                             [78072203, 78072203, 78072206, 78072207, 78072217, 78072218,
                              78072218, 78072218, 78072221, 78072222, 78072223, 78072226,
                              78072239, 78072239, 78072240, 78072241, 78072241, 78072242,
                              78072242, 78072242, 78072242, 78072242, 78072243],
                             [    3590,     3590,     3593,     3594,     3604,     3604,
                                  3604,     3604,     3607,     3608,     3609,     3612,
                                  3625,     3625,     3626,     3627,     3634,     3634,
                                  3634,     3634,     3634,     3634,     3634],
                             [   95587,    95587,    95590,    95591,    95601,    95602,
                                 95602,    95602,    95605,    95606,    95607,    95610,
                                 95623,    95623,    95624,    95625,    95625,    95625,
                                 95625,    95625,    95625,    95625,    95625],
                             [14757054, 14757055, 14757058, 14757058, 14757068, 14757069,
                              14757069, 14757069, 14757072, 14757073, 14757074, 14757077,
                              14757090, 14757098, 14757099, 14757099, 14757099, 14757099,
                              14757099, 14757099, 14757099, 14757099, 14757099]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, -9167)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3021465:3021465+29], "CCCTACACTGTCAAGTGGGAGGAGACAGT")
        self.assertEqual(alignment[0], "CCCT--ACACTGTC----AAGTGGGAGGAGACAGT--------------------")
        self.assertEqual(alignment.sequences[1].id, "cavPor2.scaffold_216473")
        self.assertEqual(len(alignment.sequences[1].seq), 10026)
        self.assertEqual(alignment.sequences[1].seq[9957:9957+28], "AGCCacactgggggtggggtgggatggt")
        self.assertEqual(alignment[1], "AGCC--acactgg-----gggtggggtgggatggt--------------------")
        self.assertEqual(alignment.column_annotations["cavPor2.scaffold_216473"], "7667--6878566-----66544895554554677--------------------")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'N')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 0)
        self.assertEqual(alignment.sequences[2].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[2].seq), 174210431)
        self.assertEqual(alignment.sequences[2].seq[16173483:16173483+24], "CCCTGTGGGGGCGATAGGAGCAGC")
        self.assertEqual(alignment[2], "-------CCCTGTG----GGGGCGATAGGAGCAGC--------------------")
        self.assertEqual(alignment.sequences[2].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[2].annotations['leftCount'], 4)
        self.assertEqual(alignment.sequences[2].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[2].annotations['rightCount'], 9)
        self.assertEqual(alignment.sequences[3].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[3].seq), 173908612)
        self.assertEqual(alignment.sequences[3].seq[16393831:16393831+24], "CCATGTGGGGGCGATAGGGGCAGC")
        self.assertEqual(alignment[3], "-------CCATGTG----GGGGCGATAGGGGCAGC--------------------")
        self.assertEqual(alignment.column_annotations["panTro2.chr6"], "-------9999999----99999999999999999--------------------")
        self.assertEqual(alignment.sequences[3].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[3].annotations['leftCount'], 4)
        self.assertEqual(alignment.sequences[3].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[3].annotations['rightCount'], 9)
        self.assertEqual(alignment.sequences[4].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[4].seq), 170899992)
        self.assertEqual(alignment.sequences[4].seq[15875265:15875265+24], "CCATGTGGGGGCGATAGGGGCAGC")
        self.assertEqual(alignment[4], "-------CCATGTG----GGGGCGATAGGGGCAGC--------------------")
        self.assertEqual(alignment.sequences[4].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[4].annotations['leftCount'], 4)
        self.assertEqual(alignment.sequences[4].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[4].annotations['rightCount'], 9)
        self.assertEqual(alignment.sequences[5].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[5].seq), 125616256)
        self.assertEqual(alignment.sequences[5].seq[78072243:78072243+31], "CCCTGAACCACATGGGGGGTGTGTGTGTGTG")
        self.assertEqual(alignment[5], "CCCTGAACCACATG----GGGGGTGTGTGTGTGTG--------------------")
        self.assertEqual(alignment.column_annotations["canFam2.chr1"], "99999999999999----99999999999999999--------------------")
        self.assertEqual(alignment.sequences[5].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[5].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[5].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[5].annotations['rightCount'], 13)
        self.assertEqual(alignment.sequences[6].id, "ornAna1.chr2")
        self.assertEqual(len(alignment.sequences[6].seq), 54797317)
        self.assertEqual(alignment.sequences[6].seq[14757099:14757099+45], "TGTGTCATAGGAGTTTGGATGTAGCCCTCTTTCATCTTTGCTGGC")
        self.assertEqual(alignment[6], "----------TGTGTCATAGGAGTTTGGATGTAGCCCTCTTTCATCTTTGCTGGC")
        self.assertEqual(alignment.sequences[6].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[6].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[6].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[6].annotations['rightCount'], 0)
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
                numpy.array([[ 3021465,  3021469,  3021469,  3021470,  3021473,  3021476,
                               3021477,  3021477,  3021494,  3021494],
                             [    9957,     9961,     9961,     9962,     9965,     9968,
                                  9968,     9968,     9985,     9985],
                             [16173483, 16173483, 16173483, 16173483, 16173486, 16173489,
                              16173490, 16173490, 16173507, 16173507],
                             [16393831, 16393831, 16393831, 16393831, 16393834, 16393837,
                              16393838, 16393838, 16393855, 16393855],
                             [15875265, 15875265, 15875265, 15875265, 15875268, 15875271,
                              15875272, 15875272, 15875289, 15875289],
                             [78072243, 78072247, 78072249, 78072250, 78072253, 78072256,
                              78072257, 78072257, 78072274, 78072274],
                             [14757099, 14757099, 14757099, 14757099, 14757099, 14757102,
                              14757103, 14757107, 14757124, 14757144]])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 15763)
        self.assertEqual(alignment.sequences[0].id, "mm9.chr10")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(alignment.sequences[0].seq[3021494:3021494+42], "TGTTTAGTACCATGCTTAGGAATGATAAACTCACTTAGTGtt")
        self.assertEqual(alignment[0], "TGTTTAGTACC----ATGCTTAGGAATGATAAACTCACTTAGTGtt")
        self.assertEqual(alignment.sequences[1].id, "ponAbe2.chr6")
        self.assertEqual(len(alignment.sequences[1].seq), 174210431)
        self.assertEqual(alignment.sequences[1].seq[16173516:16173516+46], "TGTTGCATGTCCTTTATTCTTTGGCGTGATAGGCTCACCCAATCTT")
        self.assertEqual(alignment[1], "TGTTGCATGTCCTTTATTCTTTGGCGTGATAGGCTCACCCAATCTT")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 9)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 943)
        self.assertEqual(alignment.sequences[2].id, "panTro2.chr6")
        self.assertEqual(len(alignment.sequences[2].seq), 173908612)
        self.assertEqual(alignment.sequences[2].seq[16393864:16393864+46], "TGTTGCATATCCTTTATTCTTTGGCGTGATAGGCTCACCCAATCTT")
        self.assertEqual(alignment[2], "TGTTGCATATCCTTTATTCTTTGGCGTGATAGGCTCACCCAATCTT")
        self.assertEqual(alignment.column_annotations["panTro2.chr6"], "9999999999999999999999999999999999999999999999")
        self.assertEqual(alignment.sequences[2].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[2].annotations['leftCount'], 9)
        self.assertEqual(alignment.sequences[2].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[2].annotations['rightCount'], 10)
        self.assertEqual(alignment.sequences[3].id, "hg18.chr6")
        self.assertEqual(len(alignment.sequences[3].seq), 170899992)
        self.assertEqual(alignment.sequences[3].seq[15875298:15875298+46], "TGTTGCATGTCGTTTATTCTTTGGCGTGATAGGCTCACCCAATCTT")
        self.assertEqual(alignment[3], "TGTTGCATGTCGTTTATTCTTTGGCGTGATAGGCTCACCCAATCTT")
        self.assertEqual(alignment.sequences[3].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[3].annotations['leftCount'], 9)
        self.assertEqual(alignment.sequences[3].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[3].annotations['rightCount'], 931)
        self.assertEqual(alignment.sequences[4].id, "canFam2.chr1")
        self.assertEqual(len(alignment.sequences[4].seq), 125616256)
        self.assertEqual(alignment.sequences[4].seq[78072287:78072287+46], "TGTTAAGTCTCACTTGCTGTTCAAAGTGATAGCTTCACTCCATCAT")
        self.assertEqual(alignment[4], "TGTTAAGTCTCACTTGCTGTTCAAAGTGATAGCTTCACTCCATCAT")
        self.assertEqual(alignment.column_annotations["canFam2.chr1"], "9999999999999999999999999999999999999999999999")
        self.assertEqual(alignment.sequences[4].annotations['leftStatus'], 'I')
        self.assertEqual(alignment.sequences[4].annotations['leftCount'], 13)
        self.assertEqual(alignment.sequences[4].annotations['rightStatus'], 'I')
        self.assertEqual(alignment.sequences[4].annotations['rightCount'], 1)
        self.assertEqual(alignment.sequences[5].id, "ornAna1.chr2")
        self.assertEqual(len(alignment.sequences[5].seq), 54797317)
        self.assertEqual(alignment.sequences[5].seq[14757144:14757144+36], "TGTTTAAAATGATTGCTAGAACTTCTACTCACTGGA")
        self.assertEqual(alignment[5], "TGTTTAAAATG----ATTGCTAGAACTTCTA--CTCACTGGA----")
        self.assertEqual(alignment.sequences[5].annotations['leftStatus'], 'C')
        self.assertEqual(alignment.sequences[5].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[5].annotations['rightStatus'], 'C')
        self.assertEqual(alignment.sequences[5].annotations['rightCount'], 0)
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
                numpy.array([[ 3021494,  3021505,  3021505,  3021521,  3021523,  3021532,
                               3021536],
                             [16173516, 16173527, 16173531, 16173547, 16173549, 16173558,
                              16173562],
                             [16393864, 16393875, 16393879, 16393895, 16393897, 16393906,
                              16393910],
                             [15875298, 15875309, 15875313, 15875329, 15875331, 15875340,
                              15875344],
                             [78072287, 78072298, 78072302, 78072318, 78072320, 78072329,
                              78072333],
                             [14757144, 14757155, 14757155, 14757171, 14757171, 14757180,
                              14757180]])
                )
            )


    def test_reading_missing_signature(self):
        """Test parsing MAF file ucsc_mm9_chr10_big.maf with missing signature."""
        path = "MAF/ucsc_mm9_chr10_big.maf"
        with self.assertRaises(ValueError) as cm:
            alignments = maf.AlignmentIterator(path)
        self.assertEqual(str(cm.exception), "header line does not start with ##maf")

    def test_reading_ucsc_mm9_chr10_bad(self):
        """Test parsing MAF file ucsc_mm9_chr10_bad.maf with incorrect sequence size."""
        path = "MAF/ucsc_mm9_chr10_bad.maf"
        alignments = maf.AlignmentIterator(path)
        self.assertEqual(alignments.metadata["version"], "1")
        self.assertEqual(alignments.metadata["scoring"], "autoMZ.v1")
        alignment = next(alignments)
        alignment = next(alignments)
        alignment = next(alignments)
        alignment = next(alignments)
        alignment = next(alignments)
        alignment = next(alignments)
        with self.assertRaises(ValueError) as cm:
            alignment = next(alignments)
        self.assertEqual(str(cm.exception), "sequence size is incorrect (found 219, expected 319)")

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
        self.assertEqual(alignment.sequences[0].seq[3009319:3009319+162], "TCATAGGTATTTATTTTTAAATATGGTTTGCTTTATGGCTAGAACACACCGATTACTTAAAATAGGATTAACCCCCATACACTTTAAAAATGATTAAACAACATTTCTGCTGCTCGCTCACATTCTTCATAGAAGATGACATAATGTATTTTCCTTTTGGTT")
        self.assertEqual(alignment[0], "TCATAGGTATTTATTTTTAAATATGGTTTGCTTTATGGCTAGAACACACCGATTACTTAAAATAGGATTAACC--CCCATACACTTTAAAAATGATTAAACAACATTTCTGCTGCTCGCTCACATTCTTCATAGAAGATGACATAATGTATTTTCCTTTTGGTT")
        self.assertEqual(alignment.sequences[1].id, "oryCun1.scaffold_133159")
        self.assertEqual(len(alignment.sequences[1]), 13221)
        self.assertEqual(alignment.sequences[1].seq[11087:11087+164], "TCACAGATATTTACTATTAAATATGGTTTGTTATATGGTTACGGTTCATAGGTTACTTGGAATTGGATTAACCTTCTTATTCATTGCAGAATTGGTTACACTGTGTTCTTGACCTTTGCTTGTTTTCTCCATGGAAACTGATGTCAAATACTTTCCCTTTGGTT")
        self.assertEqual(alignment[1], "TCACAGATATTTACTATTAAATATGGTTTGTTATATGGTTACGGTTCATAGGTTACTTGGAATTGGATTAACCTTCTTATTCATTGCAGAATTGGTTACACTGTGTTCTTGACCTTTGCTTGTTTTCTCCATGGAAACTGATGTCAAATACTTTCCCTTTGGTT")
        self.assertEqual(alignment.column_annotations["oryCun1.scaffold_133159"], "99569899999998999999999999999999999999999999999999999999999999999999999757878999975999999999999999979999999999997899999999999997997999999869999996999988997997999999")
        self.assertEqual(alignment.sequences[1].annotations["leftStatus"], "N")
        self.assertEqual(alignment.sequences[1].annotations["leftCount"], 0)
        self.assertEqual(alignment.sequences[1].annotations["rightStatus"], "N")
        self.assertEqual(alignment.sequences[1].annotations["rightCount"], 0)
        with self.assertRaises(ValueError) as cm:
            next(alignments)
        self.assertEqual(str(cm.exception), "sequence size is incorrect (found 219, expected 319)")

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
        self.assertEqual(alignments.metadata["speciesOrder"], ["hg16", "panTro1", "baboon", "mm4", "rn3"])
        self.assertEqual(alignments.metadata["description"], "A sample alignment")
        self.check_ucsc_test(alignments)

    def check_ucsc_test(self, alignments):
        self.assertEqual(alignments.metadata["version"], "1")
        self.assertEqual(alignments.metadata["scoring"], "tba.v8")
        self.assertEqual(alignments.metadata["comments"], [
                 "tba.v8 (((human chimp) baboon) (mouse rat))",
                 "multiz.v7",
                 "maf_project.v5 _tba_right.maf3 mouse _tba_C",
                 "single_cov2.v4 single_cov2 /dev/stdin",
         ])
        alignment = next(alignments)
        self.assertEqual(alignment.score, 23262)
        self.assertEqual(len(alignment.sequences), 5)
        self.assertEqual(alignment.sequences[0].id, "hg16.chr7")
        self.assertEqual(len(alignment.sequences[0]), 158545518)
        self.assertEqual(alignment.sequences[0].seq[27578828:27578828+38], "AAAGGGAATGTTAACCAAATGAATTGTCTCTTACGGTG")
        self.assertEqual(alignment[0], "AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG")
        self.assertEqual(alignment.sequences[1].id, "panTro1.chr6")
        self.assertEqual(len(alignment.sequences[1]), 161576975)
        self.assertEqual(alignment.sequences[1].seq[28741140:28741140+38], "AAAGGGAATGTTAACCAAATGAATTGTCTCTTACGGTG")
        self.assertEqual(alignment[1], "AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG")
        self.assertEqual(alignment.sequences[2].id, "baboon")
        self.assertEqual(len(alignment.sequences[2]), 4622798)
        self.assertEqual(alignment.sequences[2].seq[116834:116834+38], "AAAGGGAATGTTAACCAAATGAGTTGTCTCTTATGGTG")
        self.assertEqual(alignment[2], "AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG")
        self.assertEqual(alignment.sequences[3].id, "mm4.chr6")
        self.assertEqual(len(alignment.sequences[3]), 151104725)
        self.assertEqual(alignment.sequences[3].seq[53215344:53215344+38], "AATGGGAATGTTAAGCAAACGAATTGTCTCTCAGTGTG")
        self.assertEqual(alignment[3], "-AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG")
        self.assertEqual(alignment.sequences[4].id, "rn3.chr4")
        self.assertEqual(len(alignment.sequences[4]), 187371129)
        self.assertEqual(alignment.sequences[4].seq[81344243:81344243+40], "AAGGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG")
        self.assertEqual(alignment[4], "-AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([
        [27578828, 27578829, 27578831, 27578831, 27578850, 27578850, 27578866],
        [28741140, 28741141, 28741143, 28741143, 28741162, 28741162, 28741178],
        [  116834,   116835,   116837,   116837,   116856,   116856,   116872],
        [53215344, 53215344, 53215346, 53215347, 53215366, 53215366, 53215382],
        [81344243, 81344243, 81344245, 81344245, 81344264, 81344267, 81344283],
                            ])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 5062.0 )
        self.assertEqual(len(alignment.sequences), 5)
        self.assertEqual(alignment.sequences[0].id, "hg16.chr7")
        self.assertEqual(len(alignment.sequences[0]), 158545518)
        self.assertEqual(alignment.sequences[0].seq[27699739:27699739+6], "TAAAGA")
        self.assertEqual(alignment[0], "TAAAGA")
        self.assertEqual(alignment.sequences[1].id, "panTro1.chr6")
        self.assertEqual(len(alignment.sequences[1]), 161576975)
        self.assertEqual(alignment.sequences[1].seq[28862317:28862317+6], "TAAAGA")
        self.assertEqual(alignment[1], "TAAAGA")
        self.assertEqual(alignment.sequences[2].id, "baboon")
        self.assertEqual(len(alignment.sequences[2]), 4622798)
        self.assertEqual(alignment.sequences[2].seq[241163:241163+6], "TAAAGA")
        self.assertEqual(alignment[2], "TAAAGA")
        self.assertEqual(alignment.sequences[3].id, "mm4.chr6")
        self.assertEqual(len(alignment.sequences[3]), 151104725)
        self.assertEqual(alignment.sequences[3].seq[53303881:53303881+6], "TAAAGA")
        self.assertEqual(alignment[3], "TAAAGA")
        self.assertEqual(alignment.sequences[4].id, "rn3.chr4")
        self.assertEqual(len(alignment.sequences[4]), 187371129)
        self.assertEqual(alignment.sequences[4].seq[81444246:81444246+6], "taagga")
        self.assertEqual(alignment[4], "taagga")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([
                             [27699739, 27699745],
                             [28862317, 28862323],
                             [  241163,   241169],
                             [53303881, 53303887],
                             [81444246, 81444252],
                            ])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 6636.0)
        self.assertEqual(len(alignment.sequences), 4)
        self.assertEqual(alignment.sequences[0].id, "hg16.chr7")
        self.assertEqual(len(alignment.sequences[0]), 158545518)
        self.assertEqual(alignment.sequences[0].seq[27707221:27707221+13], "gcagctgaaaaca")
        self.assertEqual(alignment[0], "gcagctgaaaaca")
        self.assertEqual(alignment.sequences[1].id, "panTro1.chr6")
        self.assertEqual(len(alignment.sequences[1]), 161576975)
        self.assertEqual(alignment.sequences[1].seq[28869787:28869787+13], "gcagctgaaaaca")
        self.assertEqual(alignment[1], "gcagctgaaaaca")
        self.assertEqual(alignment.sequences[2].id, "baboon")
        self.assertEqual(len(alignment.sequences[2]), 4622798)
        self.assertEqual(alignment.sequences[2].seq[249182:249182+13], "gcagctgaaaaca")
        self.assertEqual(alignment[2], "gcagctgaaaaca")
        self.assertEqual(alignment.sequences[3].id, "mm4.chr6")
        self.assertEqual(len(alignment.sequences[3]), 151104725)
        self.assertEqual(alignment.sequences[3].seq[53310102:53310102+13], "ACAGCTGAAAATA")
        self.assertEqual(alignment[3], "ACAGCTGAAAATA")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([
                    [27707221, 27707234],
                    [28869787, 28869800],
                    [  249182,   249195],
                    [53310102, 53310115],
                            ])
                )
            )
        self.assertRaises(StopIteration, next, alignments)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
