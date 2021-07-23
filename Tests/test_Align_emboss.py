# Copyright 2008-2014 by Peter Cock.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.Align.emboss module."""
import unittest

from io import StringIO

from Bio.Align.emboss import AlignmentIterator

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align.emboss."
    ) from None


pair_example2 = """\
########################################
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


class TestEmboss(unittest.TestCase):
    def test_pair_example(self):
        # Alignment file obtained from EMBOSS:
        # http://emboss.sourceforge.net/docs/themes/alnformats/align.pair
        path = "Emboss/align.pair"
        with open(path) as stream:
            alignments = AlignmentIterator(stream)
            self.assertEqual(alignments.program, "water")
            self.assertEqual(alignments.rundate, "Wed Jan 16 17:23:19 2002")
            self.assertEqual(alignments.report_file, "stdout")
            alignments = list(alignments)
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertEqual(alignment.matrix, "EBLOSUM62")
        self.assertAlmostEqual(alignment.gap_penalty, 10.0)
        self.assertAlmostEqual(alignment.extend_penalty, 0.5)
        self.assertEqual(alignment.identity, 112)
        self.assertEqual(alignment.similarity, 112)
        self.assertEqual(alignment.gaps, 19)
        self.assertAlmostEqual(alignment.score, 591.5)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 131))
        self.assertEqual(alignment.sequences[0].id, "IXI_234")
        self.assertEqual(alignment.sequences[1].id, "IXI_235")
        self.assertEqual(
            alignment.sequences[0].seq,
            "TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "TSPASIRPPAGPSSRRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[0, 15, 24, 74, 84, 131], [0, 15, 15, 65, 65, 112]]),
            )
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            "|||||||||||||||         ||||||||||||||||||||||||||||||||||||||||||||||||||          |||||||||||||||||||||||||||||||||||||||||||||||",
        )

    def test_pair_example_nobrief(self):
        # Variation on the alignment file obtained from EMBOSS
        # (http://emboss.sourceforge.net/docs/themes/alnformats/align.pair)
        # if we include 3 sequences to align against, and we use the -nobrief
        # command line option.
        path = "Emboss/needle_nobrief_multiple.pair"
        with open(path) as stream:
            alignments = AlignmentIterator(stream)
            self.assertEqual(alignments.program, "needle")
            self.assertEqual(alignments.rundate, "Fri 23 Jul 2021 22:45:41")
            self.assertEqual(
                alignments.commandline,
                "needle -asequence seqa.fa -bsequence seqb.fa -datafile EBLOSUM62 -gapopen 10 -gapextend 0.5 -nobrief -outfile stdout",
            )
            self.assertEqual(alignments.align_format, "srspair")
            self.assertEqual(alignments.report_file, "stdout")
            alignments = list(alignments)
        self.assertEqual(len(alignments), 3)
        alignment = alignments[0]
        self.assertEqual(alignment.matrix, "EBLOSUM62")
        self.assertAlmostEqual(alignment.gap_penalty, 10.0)
        self.assertAlmostEqual(alignment.extend_penalty, 0.5)
        self.assertEqual(alignment.identity, 112)
        self.assertEqual(alignment.similarity, 112)
        self.assertEqual(alignment.gaps, 19)
        self.assertAlmostEqual(alignment.score, 591.5)
        self.assertEqual(alignment.longest_identity, "100.00%")
        self.assertEqual(alignment.longest_similarity, "100.00%")
        self.assertEqual(alignment.shortest_identity, "85.50%")
        self.assertEqual(alignment.shortest_similarity, "85.50%")
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 131))
        self.assertEqual(alignment.sequences[0].id, "IXI_234")
        self.assertEqual(alignment.sequences[1].id, "IXI_235")
        self.assertEqual(
            alignment.sequences[0].seq,
            "TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "TSPASIRPPAGPSSRRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[0, 15, 24, 74, 84, 131], [0, 15, 15, 65, 65, 112]]),
            )
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            "|||||||||||||||         ||||||||||||||||||||||||||||||||||||||||||||||||||          |||||||||||||||||||||||||||||||||||||||||||||||",
        )
        alignment = alignments[1]
        self.assertEqual(alignment.matrix, "EBLOSUM62")
        self.assertAlmostEqual(alignment.gap_penalty, 10.0)
        self.assertAlmostEqual(alignment.extend_penalty, 0.5)
        self.assertEqual(alignment.identity, 120)
        self.assertEqual(alignment.similarity, 120)
        self.assertEqual(alignment.gaps, 4)
        self.assertAlmostEqual(alignment.score, 618.0)
        self.assertEqual(alignment.longest_identity, "94.49%")
        self.assertEqual(alignment.longest_similarity, "94.49%")
        self.assertEqual(alignment.shortest_identity, "91.60%")
        self.assertEqual(alignment.shortest_similarity, "91.60%")
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 131))
        self.assertEqual(alignment.sequences[0].id, "IXI_234")
        self.assertEqual(alignment.sequences[1].id, "IXI_236")
        self.assertEqual(
            alignment.sequences[0].seq,
            "TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "TSPASIRPPAGPSSRPAMVSSRRPSPPPPRRPPGRPCCSAAPPRPQATGGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRGSRPPRFAPPLMSSCITSTTGPPPPAGDRSHE",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[0, 22, 24, 97, 99, 131], [0, 22, 22, 95, 95, 127]]),
            )
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            "||||||||||||||||||||||  |||||.||||.|||||||||.||||||||||||||||||||||||||||||||||||||||||||||||||||  ||||.||||.|||||||||||||..|||||||",
        )
        alignment = alignments[2]
        self.assertEqual(alignment.matrix, "EBLOSUM62")
        self.assertAlmostEqual(alignment.gap_penalty, 10.0)
        self.assertAlmostEqual(alignment.extend_penalty, 0.5)
        self.assertEqual(alignment.identity, 119)
        self.assertEqual(alignment.similarity, 124)
        self.assertEqual(alignment.gaps, 7)
        self.assertAlmostEqual(alignment.score, 609.0)
        self.assertEqual(alignment.longest_identity, "95.97%")
        self.assertEqual(alignment.longest_similarity, "100.00%")
        self.assertEqual(alignment.shortest_identity, "90.84%")
        self.assertEqual(alignment.shortest_similarity, "94.66%")
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 131))
        self.assertEqual(alignment.sequences[0].id, "IXI_234")
        self.assertEqual(alignment.sequences[1].id, "IXI_237")
        self.assertEqual(
            alignment.sequences[0].seq,
            "TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "TSPASLRPPAGPSSRPAMVSSRRRPSPPGPRRPTCSAAPRRPQATGGYKTCSGTCTTSTSTRHRGRSGYSARTTTAACLRASRKSMRAACSRGSRPNRFAPTLMSSCLTSTTGPPAYAGDRSHE",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [[0, 23, 24, 35, 39, 97, 99, 131], [0, 23, 23, 34, 34, 92, 92, 124]]
                ),
            )
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            "|||||:||||||||||||||||| |||||||||||    |||||||||||||:||||||||||||||||||||:|||||||||||||||||||||||  |||||||||||||||:||||||||:|||||||",
        )

    def test_simple_example(self):
        # Alignment file obtained from EMBOSS:
        # http://emboss.sourceforge.net/docs/themes/alnformats/align.simple
        path = "Emboss/align.simple"
        with open(path) as stream:
            alignments = AlignmentIterator(stream)
            self.assertEqual(alignments.program, "alignret")
            self.assertEqual(alignments.rundate, "Wed Jan 16 17:16:13 2002")
            self.assertEqual(alignments.report_file, "stdout")
            alignments = list(alignments)
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertEqual(alignment.matrix, "EBLOSUM62")
        self.assertAlmostEqual(alignment.gap_penalty, 10.0)
        self.assertAlmostEqual(alignment.extend_penalty, 0.5)
        self.assertEqual(alignment.identity, 95)
        self.assertEqual(alignment.similarity, 127)
        self.assertEqual(alignment.gaps, 25)
        self.assertAlmostEqual(alignment.score, 100.0)
        self.assertEqual(len(alignment), 4)
        self.assertEqual(alignment.shape, (4, 131))
        self.assertEqual(alignment.sequences[0].id, "IXI_234")
        self.assertEqual(alignment.sequences[1].id, "IXI_235")
        self.assertEqual(alignment.sequences[2].id, "IXI_236")
        self.assertEqual(alignment.sequences[3].id, "IXI_237")
        self.assertEqual(
            alignment.sequences[0].seq,
            "TSPASIRPPAGPSSRPAMVSSRRTRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "TSPASIRPPAGPSSRRPSPPGPRRPTGRPCCSAAPRRPQATGGWKTCSGTCTTSTSTRHRGRSGWRASRKSMRAACSRSAGSRPNRFAPTLMSSCITSTTGPPAWAGDRSHE",
        )
        self.assertEqual(
            alignment.sequences[2].seq,
            "TSPASIRPPAGPSSRPAMVSSRRPSPPPPRRPPGRPCCSAAPPRPQATGGWKTCSGTCTTSTSTRHRGRSGWSARTTTAACLRASRKSMRAACSRGSRPPRFAPPLMSSCITSTTGPPPPAGDRSHE",
        )
        self.assertEqual(
            alignment.sequences[3].seq,
            "TSPASLRPPAGPSSRPAMVSSRRRPSPPGPRRPTCSAAPRRPQATGGYKTCSGTCTTSTSTRHRGRSGYSARTTTAACLRASRKSMRAACSRGSRPNRFAPTLMSSCLTSTTGPPAYAGDRSHE",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [
                        [0, 15, 22, 23, 24, 35, 39, 74, 84, 97, 99, 131],
                        [0, 15, 15, 15, 15, 26, 30, 65, 65, 78, 80, 112],
                        [0, 15, 22, 22, 22, 33, 37, 72, 82, 95, 95, 127],
                        [0, 15, 22, 23, 23, 34, 34, 69, 79, 92, 92, 124],
                    ]
                ),
            )
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            "|||||:|||||||||:::::::  |||||:||||:::::|||||:|||||||:||||||||||||||||||||:::::::::::|||||||||||||  ||||:||||:|||||:|||||||::|||||||",
        )

    def test_pair_example2(self):
        alignments = AlignmentIterator(StringIO(pair_example2))
        self.assertEqual(alignments.program, "needle")
        self.assertEqual(alignments.rundate, "Sun 27 Apr 2007 17:20:35")
        self.assertEqual(
            alignments.commandline,
            "needle [-asequence] Spo0F.faa [-bsequence] paired_r.faa -sformat2 pearson",
        )
        self.assertEqual(alignments.align_format, "srspair")
        self.assertEqual(alignments.report_file, "ref_rec .needle")
        alignments = list(alignments)
        self.assertEqual(len(alignments), 5)
        alignment = alignments[0]
        self.assertEqual(alignment.matrix, "EBLOSUM62")
        self.assertAlmostEqual(alignment.gap_penalty, 10.0)
        self.assertAlmostEqual(alignment.extend_penalty, 0.5)
        self.assertEqual(alignment.identity, 32)
        self.assertEqual(alignment.similarity, 64)
        self.assertEqual(alignment.gaps, 17)
        self.assertAlmostEqual(alignment.score, 112.0)
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
            alignment.column_annotations["emboss_consensus"],
            " :|:.||    :.|.|::|.:  :.|.....:|.:|.||:.:..:..|.:|::|..:|||.|::..|:::....:|.|:::|.:.|...::.:.|.||..:..| ..|:|.|: ||        ",
        )
        alignment = alignments[1]
        self.assertEqual(alignment.matrix, "EBLOSUM62")
        self.assertAlmostEqual(alignment.gap_penalty, 10.0)
        self.assertAlmostEqual(alignment.extend_penalty, 0.5)
        self.assertEqual(alignment.identity, 34)
        self.assertEqual(alignment.similarity, 58)
        self.assertEqual(alignment.gaps, 9)
        self.assertAlmostEqual(alignment.score, 154.0)
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
            alignment.column_annotations["emboss_consensus"],
            " ||||||:......|:..|...|::.....|.::||:|...:..||:|.|:.:||.||:.:|:.:|.......|::|:....::|..::..||||....||...|::...|        ",
        )
        alignment = alignments[2]
        self.assertEqual(alignment.matrix, "EBLOSUM62")
        self.assertAlmostEqual(alignment.gap_penalty, 10.0)
        self.assertAlmostEqual(alignment.extend_penalty, 0.5)
        self.assertEqual(alignment.identity, 29)
        self.assertEqual(alignment.similarity, 53)
        self.assertEqual(alignment.gaps, 9)
        self.assertAlmostEqual(alignment.score, 121.0)
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
            alignment.column_annotations["emboss_consensus"],
            " .|::|||..|..:.:..||.:.|:..........|.:.:.....||.::|:.:....|:|:|:|.:|....:..:|:|....|:|...|...||:.:..||.|||.:.:..        ",
        )
        alignment = alignments[3]
        self.assertEqual(alignment.matrix, "EBLOSUM62")
        self.assertAlmostEqual(alignment.gap_penalty, 10.0)
        self.assertAlmostEqual(alignment.extend_penalty, 0.5)
        self.assertEqual(alignment.identity, 30)
        self.assertEqual(alignment.similarity, 64)
        self.assertEqual(alignment.gaps, 9)
        self.assertAlmostEqual(alignment.score, 126.0)
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
            alignment.column_annotations["emboss_consensus"],
            " :|:|:|:..:|....:.....||:...|.:|.:||.:.:|  ||.|:::.|:.:||:.|..:.:.:..|....:|:.|:.|.: :.:..:.|:.:.:.| .|||.:|    :||:.:",
        )
        alignment = alignments[4]
        self.assertEqual(alignment.matrix, "EBLOSUM62")
        self.assertAlmostEqual(alignment.gap_penalty, 10.0)
        self.assertAlmostEqual(alignment.extend_penalty, 0.5)
        self.assertEqual(alignment.identity, 35)
        self.assertEqual(alignment.similarity, 70)
        self.assertEqual(alignment.gaps, 18)
        self.assertAlmostEqual(alignment.score, 156.5)
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
            alignment.column_annotations["emboss_consensus"],
            ".:|:|:|:.|:|.|:..:.:::||...:|.:|.:||:||  :.::.|::|.|:.:.||.|.|:.:|:::...:::||.|:.|.:..:::.    |.||..|    |||..|.:...|        ",
        )

    def test_pair_example3(self):
        alignments = AlignmentIterator(StringIO(pair_example3))
        self.assertEqual(alignments.program, "needle")
        self.assertEqual(alignments.rundate, "Mon 14 Jul 2008 11:45:42")
        self.assertEqual(
            alignments.commandline,
            "needle [-asequence] asis:TGTGGTTAGGTTTGGTTTTATTGGGGGCTTGGTTTGGGCCCACCCCAAATAGGGAGTGGGGGTATGACCTCAGATAGACGAGCTTATTTTAGGGCGGCGACTATAATTATTTCGTTTCCTACAAGGATTAAAGTTTTTTCTTTTACTGTGGGAGGGGGTTTGGTATTAAGAAACGCTAGTCCGGATGTGGCTCTCCATGATACTTATTGTGTAGTAGCTCATTTTCATTATGTTCTTCGAATGGGAGCAGTCATTGGTATTTTTTTGGTTTTTTTTTGAAATTTTTAGGTTATTTAGACCATTTTTTTTTGTTTCGCTAATTAGAATTTTATTAGCCTTTGGTTTTTTTTTATTTTTTGGGGTTAAGACAAGGTGTCGTTGAATTAGTTTAGCAAAATACTGCTTAAGGTAGGCTATAGGATCTACCTTTTATCTTTCTAATCTTTTGTTTTAGTATAATTGGTCTTCGATTCAACAATTTTTAGTCTTCAGTCTTTTTTTTTATTTTGAAAAGGTTTTAACACTCTTGGTTTTGGAGGCTTTGGCTTTCTTCTTACTCTTAGGAGGATGGGCGCTAGAAAGAGTTTTAAGAGGGTGTGAAAGGGGGTTAATAGC [-bsequence] asis:TTATTAATCTTATGGTTTTGCCGTAAAATTTCTTTCTTTATTTTTTATTGTTAGGATTTTGTTGATTTTATTTTTCTCAAGAATTTTTAGGTCAATTAGACCGGCTTATTTTTTTGTCAGTGTTTAAAGTTTTATTAATTTTTGGGGGGGGGGGGAGACGGGGTGTTATCTGAATTAGTTTTTGGGAGTCTCTAGACATCTCATGGGTTGGCCGGGGGCCTGCCGTCTATAGTTCTTATTCCTTTTAAGGGAGTAAGAATTTCGATTCAGCAACTTTAGTTCACAGTCTTTTTTTTTATTAAGAAAGGTTT -filter",
        )
        self.assertEqual(alignments.align_format, "srspair")
        self.assertEqual(alignments.report_file, "stdout")
        alignments = list(alignments)
        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertEqual(alignment.matrix, "EDNAFULL")
        self.assertAlmostEqual(alignment.gap_penalty, 10.0)
        self.assertAlmostEqual(alignment.extend_penalty, 0.5)
        self.assertEqual(alignment.identity, 210)
        self.assertEqual(alignment.similarity, 210)
        self.assertEqual(alignment.gaps, 408)
        self.assertAlmostEqual(alignment.score, 561.0)
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.shape, (2, 667))
        self.assertEqual(alignment.sequences[0].id, "asis")
        self.assertEqual(alignment.sequences[1].id, "asis")
        self.assertEqual(
            alignment.sequences[0].seq,
            "TGTGGTTAGGTTTGGTTTTATTGGGGGCTTGGTTTGGGCCCACCCCAAATAGGGAGTGGGGGTATGACCTCAGATAGACGAGCTTATTTTAGGGCGGCGACTATAATTATTTCGTTTCCTACAAGGATTAAAGTTTTTTCTTTTACTGTGGGAGGGGGTTTGGTATTAAGAAACGCTAGTCCGGATGTGGCTCTCCATGATACTTATTGTGTAGTAGCTCATTTTCATTATGTTCTTCGAATGGGAGCAGTCATTGGTATTTTTTTGGTTTTTTTTTGAAATTTTTAGGTTATTTAGACCATTTTTTTTTGTTTCGCTAATTAGAATTTTATTAGCCTTTGGTTTTTTTTTATTTTTTGGGGTTAAGACAAGGTGTCGTTGAATTAGTTTAGCAAAATACTGCTTAAGGTAGGCTATAGGATCTACCTTTTATCTTTCTAATCTTTTGTTTTAGTATAATTGGTCTTCGATTCAACAATTTTTAGTCTTCAGTCTTTTTTTTTATTTTGAAAAGGTTTTAACACTCTTGGTTTTGGAGGCTTTGGCTTTCTTCTTACTCTTAGGAGGATGGGCGCTAGAAAGAGTTTTAAGAGGGTGTGAAAGGGGGTTAATAGC",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "TTATTAATCTTATGGTTTTGCCGTAAAATTTCTTTCTTTATTTTTTATTGTTAGGATTTTGTTGATTTTATTTTTCTCAAGAATTTTTAGGTCAATTAGACCGGCTTATTTTTTTGTCAGTGTTTAAAGTTTTATTAATTTTTGGGGGGGGGGGGAGACGGGGTGTTATCTGAATTAGTTTTTGGGAGTCTCTAGACATCTCATGGGTTGGCCGGGGGCCTGCCGTCTATAGTTCTTATTCCTTTTAAGGGAGTAAGAATTTCGATTCAGCAACTTTAGTTCACAGTCTTTTTTTTTATTAAGAAAGGTTT",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [
                        [
                            0,
                            162,
                            169,
                            201,
                            210,
                            210,
                            220,
                            222,
                            236,
                            240,
                            244,
                            253,
                            277,
                            277,
                            278,
                            279,
                            300,
                            300,
                            310,
                            310,
                            314,
                            320,
                            334,
                            351,
                            357,
                            357,
                            379,
                            379,
                            390,
                            403,
                            405,
                            407,
                            418,
                            418,
                            442,
                            442,
                            447,
                            447,
                            448,
                            452,
                            455,
                            455,
                            460,
                            465,
                            478,
                            479,
                            509,
                            510,
                            518,
                            615,
                        ],
                        [
                            0,
                            0,
                            7,
                            7,
                            16,
                            22,
                            32,
                            32,
                            46,
                            46,
                            50,
                            50,
                            74,
                            80,
                            81,
                            81,
                            102,
                            107,
                            117,
                            119,
                            123,
                            123,
                            137,
                            137,
                            143,
                            147,
                            169,
                            170,
                            181,
                            181,
                            183,
                            183,
                            194,
                            215,
                            239,
                            241,
                            246,
                            250,
                            251,
                            251,
                            254,
                            255,
                            260,
                            260,
                            273,
                            273,
                            303,
                            303,
                            311,
                            311,
                        ],
                    ]
                ),
            )
        )
        self.assertEqual(
            alignment.column_annotations["emboss_consensus"],
            "                                                                                                                                                                  .||||||                                .|||||.||      |||..|..||  ||||.||||.||.|    ||.|         ||.|.|||||.|||.||||.||||      | |||||||||||.|.|||||||     ||||||||.|  ||.|      |||.|.||||||||                 ||||||    .||||...||||..|||||..| |||||||||||             ||  ||.||.||.||                     ||..||.||.|.|||..||||.||  |||||    |    ||| |.|||     |||||||||.||| .||||||...|||||||||||||||||..| ||||||||                                                                                                 ",
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
