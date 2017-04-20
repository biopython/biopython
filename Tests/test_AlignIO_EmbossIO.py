# Copyright 2008-2014 by Peter Cock.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.AlignIO.EmbossIO"""

import unittest

from Bio._py3k import StringIO

from Bio.AlignIO.EmbossIO import EmbossWriter, EmbossIterator

# http://emboss.sourceforge.net/docs/themes/alnformats/align.simple
simple_example = \
"""########################################
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

"""

# http://emboss.sourceforge.net/docs/themes/alnformats/align.pair
pair_example = \
"""########################################
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

pair_example2 = \
"""########################################
# Program: needle
# Rundate: Sun 27 Apr 2007 17:20:35
# Commandline: needle
#    [-asequence] Spo0F.faa
#    [-bsequence] paired_r.faa
#    -sformat2 pearson
# Align_format: srspair
# Report_file: ref_rec .needle
########################################

#=======================================
#
# Aligned_sequences: 2
# 1: ref_rec
# 2: gi|94968718|receiver
# Matrix: EBLOSUM62
# Gap_penalty: 10.0
# Extend_penalty: 0.5
#
# Length: 124
# Identity:      32/124 (25.8%)
# Similarity:    64/124 (51.6%)
# Gaps:          17/124 (13.7%)
# Score: 112.0
# 
#
#=======================================

ref_rec            1 KILIVDD----QYGIRILLNEVFNKEGYQTFQAANGLQALDIVTKERPDL     46
                      :|:.||    :.|.|::|.:  :.|.....:|.:|.||:.:..:..|.:
gi|94968718|r      1 -VLLADDHALVRRGFRLMLED--DPEIEIVAEAGDGAQAVKLAGELHPRV     47

ref_rec           47 VLLDMKIPGMDGIEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALT     96
                     |::|..:|||.|::..|:::....:|.|:::|.:.|...::.:.|.||..
gi|94968718|r     48 VVMDCAMPGMSGMDATKQIRTQWPDIAVLMLTMHSEDTWVRLALEAGANG     97

ref_rec           97 HFAK-PFDIDEIRDAV--------    111
                     :..| ..|:|.|: ||        
gi|94968718|r     98 YILKSAIDLDLIQ-AVRRVANGET    120


#=======================================
#
# Aligned_sequences: 2
# 1: ref_rec
# 2: gi|94968761|receiver
# Matrix: EBLOSUM62
# Gap_penalty: 10.0
# Extend_penalty: 0.5
#
# Length: 119
# Identity:      34/119 (28.6%)
# Similarity:    58/119 (48.7%)
# Gaps:           9/119 ( 7.6%)
# Score: 154.0
# 
#
#=======================================

ref_rec            1 KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIVTKERPDLVLLD     50
                      ||||||:......|:..|...|::.....|.::||:|...:..||:|.|
gi|94968761|r      1 -ILIVDDEANTLASLSRAFRLAGHEATVCDNAVRALEIAKSKPFDLILSD     49

ref_rec           51 MKIPGMDGIEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHFAK    100
                     :.:||.||:.:|:.:|.......|::|:....::|..::..||||....|
gi|94968761|r     50 VVMPGRDGLTLLEDLKTAGVQAPVVMMSGQAHIEMAVKATRLGALDFLEK     99

ref_rec          101 PFDIDEIRDAV--------    111
                     |...|::...|        
gi|94968761|r    100 PLSTDKLLLTVENALKLKR    118


#=======================================
#
# Aligned_sequences: 2
# 1: ref_rec
# 2: gi|94967506|receiver
# Matrix: EBLOSUM62
# Gap_penalty: 10.0
# Extend_penalty: 0.5
#
# Length: 120
# Identity:      29/120 (24.2%)
# Similarity:    53/120 (44.2%)
# Gaps:           9/120 ( 7.5%)
# Score: 121.0
# 
#
#=======================================

ref_rec            1 -KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIVTKERPDLVLL     49
                      .|::|||..|..:.:..||.:.|:..........|.:.:.....||.::
gi|94967506|r      1 LHIVVVDDDPGTCVYIESVFAELGHTCKSFVRPEAAEEYILTHPVDLAIV     50

ref_rec           50 DMKIPGMDGIEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHFA     99
                     |:.:....|:|:|:|.:|....:..:|:|....|:|...|...||:.:..
gi|94967506|r     51 DVYLGSTTGVEVLRRCRVHRPKLYAVIITGQISLEMAARSIAEGAVDYIQ    100

ref_rec          100 KPFDIDEIRDAV--------    111
                     ||.|||.:.:..        
gi|94967506|r    101 KPIDIDALLNIAERALEHKE    120


#=======================================
#
# Aligned_sequences: 2
# 1: ref_rec
# 2: gi|94970045|receiver
# Matrix: EBLOSUM62
# Gap_penalty: 10.0
# Extend_penalty: 0.5
#
# Length: 118
# Identity:      30/118 (25.4%)
# Similarity:    64/118 (54.2%)
# Gaps:           9/118 ( 7.6%)
# Score: 126.0
# 
#
#=======================================

ref_rec            1 KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIVTK--ERPDLVL     48
                      :|:|:|:..:|....:.....||:...|.:|.:||.:.:|  ||.|:::
gi|94970045|r      1 -VLLVEDEEALRAAAGDFLETRGYKIMTARDGTEALSMASKFAERIDVLI     49

ref_rec           49 LDMKIPGMDGIEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHF     98
                     .|:.:||:.|..:.:.:..|....:|:.|:.|.: :.:..:.|:.:.:.|
gi|94970045|r     50 TDLVMPGISGRVLAQELVKIHPETKVMYMSGYDD-ETVMVNGEIDSSSAF     98

ref_rec           99 -AKPFDID----EIRDAV    111
                      .|||.:|    :||:.:
gi|94970045|r     99 LRKPFRMDALSAKIREVL    116


#=======================================
#
# Aligned_sequences: 2
# 1: ref_rec
# 2: gi|94970041|receiver
# Matrix: EBLOSUM62
# Gap_penalty: 10.0
# Extend_penalty: 0.5
#
# Length: 125
# Identity:      35/125 (28.0%)
# Similarity:    70/125 (56.0%)
# Gaps:          18/125 (14.4%)
# Score: 156.5
# 
#
#=======================================

ref_rec            1 KILIVDDQYGIRILLNEVFNKEGYQTFQAANGLQALDIV--TKERPDLVL     48
                     .:|:|:|:.|:|.|:..:.:::||...:|.:|.:||:||  :.::.|::|
gi|94970041|r      1 TVLLVEDEEGVRKLVRGILSRQGYHVLEATSGEEALEIVRESTQKIDMLL     50

ref_rec           49 LDMKIPGMDGIEILKRMKVIDENIRVIIMTAYGELDMIQESKELGALTHF     98
                     .|:.:.||.|.|:.:|:::...:::||.|:.|.:..:::.    |.||..
gi|94970041|r     51 SDVVLVGMSGRELSERLRIQMPSLKVIYMSGYTDDAIVRH----GVLTES     96

ref_rec           99 A----KPFDIDEIRDAV--------    111
                     |    |||..|.:...|        
gi|94970041|r     97 AEFLQKPFTSDSLLRKVRAVLQKRQ    121


#---------------------------------------
#---------------------------------------

"""  # noqa : W291

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
        self.assertEqual([r.id for r in alignments[0]], ["IXI_234", "IXI_235", "IXI_236", "IXI_237"])

    def test_pair_plus_simple(self):
        alignments = list(EmbossIterator(StringIO(pair_example + simple_example)))
        self.assertEqual(len(alignments), 2)
        self.assertEqual(len(alignments[0]), 2)
        self.assertEqual(len(alignments[1]), 4)
        self.assertEqual([r.id for r in alignments[0]], ["IXI_234", "IXI_235"])
        self.assertEqual([r.id for r in alignments[1]], ["IXI_234", "IXI_235", "IXI_236", "IXI_237"])

    def test_pair_example2(self):
        alignments = list(EmbossIterator(StringIO(pair_example2)))
        self.assertEqual(len(alignments), 5)
        self.assertEqual(len(alignments[0]), 2)
        self.assertEqual([r.id for r in alignments[0]], ["ref_rec", "gi|94968718|receiver"])
        self.assertEqual([r.id for r in alignments[4]], ["ref_rec", "gi|94970041|receiver"])

    def test_pair_example3(self):
        alignments = list(EmbossIterator(StringIO(pair_example3)))
        self.assertEqual(len(alignments), 1)
        self.assertEqual(len(alignments[0]), 2)
        self.assertEqual([r.id for r in alignments[0]], ["asis", "asis"])


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
