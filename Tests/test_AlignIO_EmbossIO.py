# Copyright 2008-2014 by Peter Cock.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.AlignIO.EmbossIO module."""
import unittest

from io import StringIO

from Bio.AlignIO.EmbossIO import EmbossIterator

# http://emboss.sourceforge.net/docs/themes/alnformats/align.simple
simple_example = """\
########################################
# Program:  alignret
# Rundate:  Wed Jan 16 17:16:13 2002
# Report_file: stdout
########################################
#=======================================
#
# Aligned_sequences: 4
# 1: IXI_234
# 2: IXI_235
# 3: IXI_236
# 4: IXI_237
# Matrix: EBLOSUM62
# Gap_penalty: 10.0
# Extend_penalty: 0.5
#
# Length: 131
# Identity:      95/131 (72.5%)
# Similarity:   127/131 (96.9%)
# Gaps:          25/131 (19.1%)
# Score: 100.0
#
#
#=======================================

IXI_234            1 TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQAT     50
IXI_235            1 TSPASIRPPAGPSSR---------RPSPPGPRRPTGRPCCSAAPRRPQAT     41
IXI_236            1 TSPASIRPPAGPSSRPAMVSSR--RPSPPPPRRPPGRPCCSAAPPRPQAT     48
IXI_237            1 TSPASLRPPAGPSSRPAMVSSRR-RPSPPGPRRPT----CSAAPRRPQAT     45
                     |||||:|||||||||:::::::  |||||:||||:::::|||||:|||||

IXI_234           51 GGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAG    100
IXI_235           42 GGWKTCSGTCTTSTSTRHRGRSGW----------RASRKSMRAACSRSAG     81
IXI_236           49 GGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSR--G     96
IXI_237           46 GGYKTCSGTCTTSTSTRHRGRSGYSARTTTAACLRASRKSMRAACSR--G     93
                     ||:||||||||||||||||||||:::::::::::|||||||||||||  |

IXI_234          101 SRPNRFAPTLMSSCITSTTGPPAWAGDRSHE    131
IXI_235           82 SRPNRFAPTLMSSCITSTTGPPAWAGDRSHE    112
IXI_236           97 SRPPRFAPPLMSSCITSTTGPPPPAGDRSHE    127
IXI_237           94 SRPNRFAPTLMSSCLTSTTGPPAYAGDRSHE    124
                     |||:||||:|||||:|||||||::|||||||


#---------------------------------------
#---------------------------------------

"""  # noqa: E122 not clear to me, why this comes up here

# http://emboss.sourceforge.net/docs/themes/alnformats/align.pair
pair_example = """\
########################################
# Program:  water
# Rundate:  Wed Jan 16 17:23:19 2002
# Report_file: stdout
########################################
#=======================================
#
# Aligned_sequences: 2
# 1: IXI_234
# 2: IXI_235
# Matrix: EBLOSUM62
# Gap_penalty: 10.0
# Extend_penalty: 0.5
#
# Length: 131
# Identity:     112/131 (85.5%)
# Similarity:   112/131 (85.5%)
# Gaps:          19/131 (14.5%)
# Score: 591.5
#
#
#=======================================

IXI_234            1 TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQAT     50
                     |||||||||||||||         ||||||||||||||||||||||||||
IXI_235            1 TSPASIRPPAGPSSR---------RPSPPGPRRPTGRPCCSAAPRRPQAT     41

IXI_234           51 GGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAG    100
                     ||||||||||||||||||||||||          ||||||||||||||||
IXI_235           42 GGWKTCSGTCTTSTSTRHRGRSGW----------RASRKSMRAACSRSAG     81

IXI_234          101 SRPNRFAPTLMSSCITSTTGPPAWAGDRSHE    131
                     |||||||||||||||||||||||||||||||
IXI_235           82 SRPNRFAPTLMSSCITSTTGPPAWAGDRSHE    112


#---------------------------------------
#---------------------------------------       


"""  # noqa : W291


with open("Emboss/needle.txt") as handle:
    pair_example2 = handle.read()


pair_example3 = """########################################
# Program: needle
# Rundate: Mon 14 Jul 2008 11:45:42
# Commandline: needle
#    [-asequence] asis:TGTGGTTAGGTTTGGTTTTATTGGGGGCTTGGTTTGGGCCCACCCCAAATAGGGAGTGGGGGTATGACCTCAGATAGACGAGCTTATTTTAGGGCGGCGACTATAATTATTTCGTTTCCTACAAGGATTAAAGTTTTTTCTTTTACTGTGGGAGGGGGTTTGGTATTAAGAAACGCTAGTCCGGATGTGGCTCTCCATGATACTTATTGTGTAGTAGCTCATTTTCATTATGTTCTTCGAATGGGAGCAGTCATTGGTATTTTTTTGGTTTTTTTTTGAAATTTTTAGGTTATTTAGACCATTTTTTTTTGTTTCGCTAATTAGAATTTTATTAGCCTTTGGTTTTTTTTTATTTTTTGGGGTTAAGACAAGGTGTCGTTGAATTAGTTTAGCAAAATACTGCTTAAGGTAGGCTATAGGATCTACCTTTTATCTTTCTAATCTTTTGTTTTAGTATAATTGGTCTTCGATTCAACAATTTTTAGTCTTCAGTCTTTTTTTTTATTTTGAAAAGGTTTTAACACTCTTGGTTTTGGAGGCTTTGGCTTTCTTCTTACTCTTAGGAGGATGGGCGCTAGAAAGAGTTTTAAGAGGGTGTGAAAGGGGGTTAATAGC
#    [-bsequence] asis:TTATTAATCTTATGGTTTTGCCGTAAAATTTCTTTCTTTATTTTTTATTGTTAGGATTTTGTTGATTTTATTTTTCTCAAGAATTTTTAGGTCAATTAGACCGGCTTATTTTTTTGTCAGTGTTTAAAGTTTTATTAATTTTTGGGGGGGGGGGGAGACGGGGTGTTATCTGAATTAGTTTTTGGGAGTCTCTAGACATCTCATGGGTTGGCCGGGGGCCTGCCGTCTATAGTTCTTATTCCTTTTAAGGGAGTAAGAATTTCGATTCAGCAACTTTAGTTCACAGTCTTTTTTTTTATTAAGAAAGGTTT
#    -filter
# Align_format: srspair
# Report_file: stdout
########################################

#=======================================
#
# Aligned_sequences: 2
# 1: asis
# 2: asis
# Matrix: EDNAFULL
# Gap_penalty: 10.0
# Extend_penalty: 0.5
#
# Length: 667
# Identity:     210/667 (31.5%)
# Similarity:   210/667 (31.5%)
# Gaps:         408/667 (61.2%)
# Score: 561.0
# 
#
#=======================================

asis               1 TGTGGTTAGGTTTGGTTTTATTGGGGGCTTGGTTTGGGCCCACCCCAAAT     50
                                                                       
asis               0 --------------------------------------------------      0

asis              51 AGGGAGTGGGGGTATGACCTCAGATAGACGAGCTTATTTTAGGGCGGCGA    100
                                                                       
asis               0 --------------------------------------------------      0

asis             101 CTATAATTATTTCGTTTCCTACAAGGATTAAAGTTTTTTCTTTTACTGTG    150
                                                                       
asis               0 --------------------------------------------------      0

asis             151 GGAGGGGGTTTGGTATTAAGAAACGCTAGTCCGGATGTGGCTCTCCATGA    200
                                 .||||||                               
asis               1 ------------TTATTAA-------------------------------      7

asis             201 TACTTATTGT------GTAGTAGCTCATTTTCATTATGTTCTTCGAATGG    244
                      .|||||.||      |||..|..||  ||||.||||.||.|    ||.|
asis               8 -TCTTATGGTTTTGCCGTAAAATTTC--TTTCTTTATTTTTT----ATTG     50

asis             245 GAGCAGTCATTGGTATTTTTTTGGTTTTTTTTT------GAAATTTTTAG    288
                              ||.|.|||||.|||.||||.||||      | |||||||||
asis              51 ---------TTAGGATTTTGTTGATTTTATTTTTCTCAAG-AATTTTTAG     90

asis             289 GTTATTTAGACC-----ATTTTTTTTT--GTTTCGCTAATTAGAATTTTA    331
                     ||.|.|||||||     ||||||||.|  ||.|      |||.|.|||||
asis              91 GTCAATTAGACCGGCTTATTTTTTTGTCAGTGT------TTAAAGTTTTA    134

asis             332 TTAGCCTTTGGTTTTTTTTTATTTTT----TGGGGTTAAGACAAGGTGTC    377
                     |||                 ||||||    .||||...||||..|||||.
asis             135 TTA-----------------ATTTTTGGGGGGGGGGGGAGACGGGGTGTT    167

asis             378 GT-TGAATTAGTTTAGCAAAATACTGCTTAAGGTAGGCTATA--------    418
                     .| |||||||||||             ||  ||.||.||.||        
asis             168 ATCTGAATTAGTTT-------------TT--GGGAGTCTCTAGACATCTC    202

asis             419 -------------GGATCTACCTTTTATCTTTCTAAT--CTTTT----GT    449
                                  ||..||.||.|.|||..||||.||  |||||    | 
asis             203 ATGGGTTGGCCGGGGGCCTGCCGTCTATAGTTCTTATTCCTTTTAAGGG-    251

asis             450 TTTAGT-ATAATTGGTCTTCGATTCAACAATTTTTAGTCTTCAGTCTTTT    498
                        ||| |.|||     |||||||||.||| .||||||...|||||||||
asis             252 ---AGTAAGAAT-----TTCGATTCAGCAA-CTTTAGTTCACAGTCTTTT    292

asis             499 TTTTTATTTTGAAAAGGTTTTAACACTCTTGGTTTTGGAGGCTTTGGCTT    548
                     ||||||||..| ||||||||                              
asis             293 TTTTTATTAAG-AAAGGTTT------------------------------    311

asis             549 TCTTCTTACTCTTAGGAGGATGGGCGCTAGAAAGAGTTTTAAGAGGGTGT    598
                                                                       
asis             311 --------------------------------------------------    311

asis             599 GAAAGGGGGTTAATAGC    615
                                      
asis             311 -----------------    311


#---------------------------------------
#---------------------------------------"""  # noqa : W291


class TestEmbossIO(unittest.TestCase):
    def test_pair_example(self):
        alignments = list(EmbossIterator(StringIO(pair_example)))
        self.assertEqual(len(alignments), 1)
        self.assertEqual(len(alignments[0]), 2)
        self.assertEqual([r.id for r in alignments[0]], ["IXI_234", "IXI_235"])

    def test_simple_example(self):
        alignments = list(EmbossIterator(StringIO(simple_example)))
        self.assertEqual(len(alignments), 1)
        self.assertEqual(len(alignments[0]), 4)
        self.assertEqual(
            [r.id for r in alignments[0]], ["IXI_234", "IXI_235", "IXI_236", "IXI_237"]
        )

    def test_pair_plus_simple(self):
        alignments = list(EmbossIterator(StringIO(pair_example + simple_example)))
        self.assertEqual(len(alignments), 2)
        self.assertEqual(len(alignments[0]), 2)
        self.assertEqual(len(alignments[1]), 4)
        self.assertEqual([r.id for r in alignments[0]], ["IXI_234", "IXI_235"])
        self.assertEqual(
            [r.id for r in alignments[1]], ["IXI_234", "IXI_235", "IXI_236", "IXI_237"]
        )

    def test_pair_example2(self):
        alignments = list(EmbossIterator(StringIO(pair_example2)))
        self.assertEqual(len(alignments), 5)
        self.assertEqual(len(alignments[0]), 2)
        self.assertEqual(
            [r.id for r in alignments[0]], ["ref_rec", "gi|94968718|receiver"]
        )
        self.assertEqual(
            [r.id for r in alignments[4]], ["ref_rec", "gi|94970041|receiver"]
        )

    def test_pair_example3(self):
        alignments = list(EmbossIterator(StringIO(pair_example3)))
        self.assertEqual(len(alignments), 1)
        self.assertEqual(len(alignments[0]), 2)
        self.assertEqual([r.id for r in alignments[0]], ["asis", "asis"])


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
