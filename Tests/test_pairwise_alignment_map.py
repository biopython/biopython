# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


"""Tests for mapping pairwise alignments."""

import os
import random
import unittest

from Bio.Seq import Seq
from Bio.Align import PairwiseAligner


class TestSimple(unittest.TestCase):
    def setUp(self):
        aligner = PairwiseAligner()
        aligner.internal_open_gap_score = -1
        aligner.internal_extend_gap_score = -0.0
        aligner.match_score = +1
        aligner.mismatch_score = -1
        aligner.mode = "local"
        self.aligner = aligner

    def test_internal(self):
        aligner = self.aligner
        chromosome = Seq("AAAAAAAAAAAAGGGGGGGCCCCCGGGGGGAAAAAAAAAA")
        chromosome.id = "chromosome"
        transcript = Seq("GGGGGGGCCCCCGGGGGGA")
        transcript.id = "transcript"
        sequence = Seq("GGCCCCCGGG")
        sequence.id = "sequence"
        alignments1 = aligner.align(chromosome, transcript)
        self.assertEqual(len(alignments1), 1)
        alignment1 = alignments1[0]
        self.assertEqual(alignment1.path, ((12, 0), (31, 19)))
        self.assertEqual(
            str(alignment1),
            """\
AAAAAAAAAAAAGGGGGGGCCCCCGGGGGGAAAAAAAAAA
            |||||||||||||||||||         
            GGGGGGGCCCCCGGGGGGA         
""",  # noqa: W291
        )
        alignments2 = aligner.align(transcript, sequence)
        self.assertEqual(len(alignments2), 1)
        alignment2 = alignments2[0]
        self.assertEqual(alignment2.path, ((5, 0), (15, 10)))
        self.assertEqual(
            str(alignment2),
            """\
GGGGGGGCCCCCGGGGGGA
     ||||||||||    
     GGCCCCCGGG    
""",  # noqa: W291
        )
        alignment = alignment1.map(alignment2)
        self.assertEqual(len(alignment.path), 2)
        self.assertSequenceEqual(alignment.path[0], [17, 0])
        self.assertSequenceEqual(alignment.path[1], [27, 10])
        self.assertEqual(
            str(alignment),
            """\
AAAAAAAAAAAAGGGGGGGCCCCCGGGGGGAAAAAAAAAA
                 ||||||||||             
                 GGCCCCCGGG             
""",  # noqa: W291
        )
        psl = format(alignment, "psl")
        self.assertEqual(
            psl,
            """\
10	0	0	0	0	0	0	0	+	sequence	10	0	10	chromosome	40	17	27	1	10,	0,	17,
""",
        )

    def test_left_overhang(self):
        aligner = self.aligner
        chromosome = Seq("GGGCCCCCGGGGGGAAAAAAAAAA")
        chromosome.id = "chromosome"
        transcript = Seq("AGGGGGCCCCCGGGGGGA")
        transcript.id = "transcript"
        sequence = Seq("GGGGGCCCCCGGG")
        sequence.id = "sequence"
        alignments1 = aligner.align(chromosome, transcript)
        self.assertEqual(len(alignments1), 1)
        alignment1 = alignments1[0]
        self.assertEqual(
            str(alignment1),
            """\
   GGGCCCCCGGGGGGAAAAAAAAAA
   |||||||||||||||         
AGGGGGCCCCCGGGGGGA         
""",  # noqa: W291
        )
        alignments2 = aligner.align(transcript, sequence)
        self.assertEqual(len(alignments2), 1)
        alignment2 = alignments2[0]
        self.assertEqual(
            str(alignment2),
            """\
AGGGGGCCCCCGGGGGGA
 |||||||||||||    
 GGGGGCCCCCGGG    
""",  # noqa: W291
        )
        alignment = alignment1.map(alignment2)
        self.assertEqual(len(alignment.path), 2)
        self.assertSequenceEqual(alignment.path[0], [0, 2])
        self.assertSequenceEqual(alignment.path[1], [11, 13])
        self.assertEqual(
            str(alignment),
            """\
  GGGCCCCCGGGGGGAAAAAAAAAA
  |||||||||||             
GGGGGCCCCCGGG             
""",  # noqa: W291
        )
        psl = format(alignment, "psl")
        self.assertEqual(
            psl,
            """\
11	0	0	0	0	0	0	0	+	sequence	13	2	13	chromosome	24	0	11	1	11,	2,	0,
""",
        )

    def test_right_overhang(self):
        aligner = self.aligner
        chromosome = Seq("AAAAAAAAAAAAGGGGGGGCCCCCGGG")
        chromosome.id = "chromosome"
        transcript = Seq("GGGGGGGCCCCCGGGGGGA")
        transcript.id = "transcript"
        sequence = Seq("GGCCCCCGGGGG")
        sequence.id = "sequence"
        alignments1 = aligner.align(chromosome, transcript)
        self.assertEqual(len(alignments1), 1)
        alignment1 = alignments1[0]
        self.assertEqual(
            str(alignment1),
            """\
AAAAAAAAAAAAGGGGGGGCCCCCGGG    
            |||||||||||||||    
            GGGGGGGCCCCCGGGGGGA
""",  # noqa: W291
        )
        alignments2 = aligner.align(transcript, sequence)
        self.assertEqual(len(alignments2), 1)
        alignment2 = alignments2[0]
        self.assertEqual(
            str(alignment2),
            """\
GGGGGGGCCCCCGGGGGGA
     ||||||||||||  
     GGCCCCCGGGGG  
""",  # noqa: W291
        )
        alignment = alignment1.map(alignment2)
        self.assertEqual(len(alignment.path), 2)
        self.assertSequenceEqual(alignment.path[0], [17, 0])
        self.assertSequenceEqual(alignment.path[1], [27, 10])
        self.assertEqual(
            str(alignment),
            """\
AAAAAAAAAAAAGGGGGGGCCCCCGGG  
                 ||||||||||  
                 GGCCCCCGGGGG
""",  # noqa: W291
        )
        psl = format(alignment, "psl")
        self.assertEqual(
            psl,
            """\
10	0	0	0	0	0	0	0	+	sequence	12	0	10	chromosome	27	17	27	1	10,	0,	17,
""",
        )

    def test_reverse_transcript(self):
        aligner = self.aligner
        chromosome = Seq("AAAAAAAAAAAAGGGGGGGCCCCCGGGGGGAAAAAAAAAA")
        chromosome.id = "chromosome"
        transcript = Seq("TCCCCCCGGGGGCCCCCCC")
        transcript.id = "transcript"
        sequence = Seq("GGCCCCCGGG")
        sequence.id = "sequence"
        alignments1 = aligner.align(chromosome, transcript, strand="-")
        self.assertEqual(len(alignments1), 1)
        alignment1 = alignments1[0]
        self.assertEqual(alignment1.path, ((12, 19), (31, 0)))
        self.assertEqual(
            str(alignment1),
            """\
AAAAAAAAAAAAGGGGGGGCCCCCGGGGGGAAAAAAAAAA
            |||||||||||||||||||         
            GGGGGGGCCCCCGGGGGGA         
""",  # noqa: W291
        )
        alignments2 = aligner.align(transcript, sequence, strand="-")
        self.assertEqual(len(alignments2), 1)
        alignment2 = alignments2[0]
        self.assertEqual(alignment2.path, ((4, 10), (14, 0)))
        self.assertEqual(
            str(alignment2),
            """\
TCCCCCCGGGGGCCCCCCC
    ||||||||||     
    CCCGGGGGCC     
""",  # noqa: W291
        )
        alignment = alignment1.map(alignment2)
        self.assertEqual(len(alignment.path), 2)
        self.assertSequenceEqual(alignment.path[0], [17, 0])
        self.assertSequenceEqual(alignment.path[1], [27, 10])
        self.assertEqual(
            str(alignment),
            """\
AAAAAAAAAAAAGGGGGGGCCCCCGGGGGGAAAAAAAAAA
                 ||||||||||             
                 GGCCCCCGGG             
""",  # noqa: W291
        )
        psl = format(alignment, "psl")
        self.assertEqual(
            psl,
            """\
10	0	0	0	0	0	0	0	+	sequence	10	0	10	chromosome	40	17	27	1	10,	0,	17,
""",
        )

    def test_reverse_sequence(self):
        aligner = self.aligner
        chromosome = Seq("AAAAAAAAAAAAGGGGGGGCCCCCGGGGGGAAAAAAAAAA")
        chromosome.id = "chromosome"
        transcript = Seq("GGGGGGGCCCCCGGGGGGA")
        transcript.id = "transcript"
        sequence = Seq("CCCGGGGGCC")
        sequence.id = "sequence"
        alignments1 = aligner.align(chromosome, transcript)
        self.assertEqual(len(alignments1), 1)
        alignment1 = alignments1[0]
        self.assertEqual(alignment1.path, ((12, 0), (31, 19)))
        self.assertEqual(
            str(alignment1),
            """\
AAAAAAAAAAAAGGGGGGGCCCCCGGGGGGAAAAAAAAAA
            |||||||||||||||||||         
            GGGGGGGCCCCCGGGGGGA         
""",  # noqa: W291
        )
        alignments2 = aligner.align(transcript, sequence, "-")
        self.assertEqual(len(alignments2), 1)
        alignment2 = alignments2[0]
        self.assertEqual(alignment2.path, ((5, 10), (15, 0)))
        self.assertEqual(
            str(alignment2),
            """\
GGGGGGGCCCCCGGGGGGA
     ||||||||||    
     GGCCCCCGGG    
""",  # noqa: W291
        )
        alignment = alignment1.map(alignment2)
        self.assertEqual(len(alignment.path), 2)
        self.assertSequenceEqual(alignment.path[0], [17, 10])
        self.assertSequenceEqual(alignment.path[1], [27, 0])
        self.assertEqual(
            str(alignment),
            """\
AAAAAAAAAAAAGGGGGGGCCCCCGGGGGGAAAAAAAAAA
                 ||||||||||             
                 GGCCCCCGGG             
""",  # noqa: W291
        )
        psl = format(alignment, "psl")
        self.assertEqual(
            psl,
            """\
10	0	0	0	0	0	0	0	-	sequence	10	0	10	chromosome	40	17	27	1	10,	0,	17,
""",
        )

    def test_reverse_transcript_sequence(self):
        aligner = self.aligner
        chromosome = Seq("AAAAAAAAAAAAGGGGGGGCCCCCGGGGGGAAAAAAAAAA")
        chromosome.id = "chromosome"
        transcript = Seq("TCCCCCCGGGGGCCCCCCC")
        transcript.id = "transcript"
        sequence = Seq("CCCGGGGGCC")
        sequence.id = "sequence"
        alignments1 = aligner.align(chromosome, transcript, "-")
        self.assertEqual(len(alignments1), 1)
        alignment1 = alignments1[0]
        self.assertEqual(alignment1.path, ((12, 19), (31, 0)))
        self.assertEqual(
            str(alignment1),
            """\
AAAAAAAAAAAAGGGGGGGCCCCCGGGGGGAAAAAAAAAA
            |||||||||||||||||||         
            GGGGGGGCCCCCGGGGGGA         
""",  # noqa: W291
        )
        alignments2 = aligner.align(transcript, sequence)
        self.assertEqual(len(alignments2), 1)
        alignment2 = alignments2[0]
        self.assertEqual(alignment2.path, ((4, 0), (14, 10)))
        self.assertEqual(
            str(alignment2),
            """\
TCCCCCCGGGGGCCCCCCC
    ||||||||||     
    CCCGGGGGCC     
""",  # noqa: W291
        )
        alignment = alignment1.map(alignment2)
        self.assertEqual(len(alignment.path), 2)
        self.assertSequenceEqual(alignment.path[0], [17, 10])
        self.assertSequenceEqual(alignment.path[1], [27, 0])
        self.assertEqual(
            str(alignment),
            """\
AAAAAAAAAAAAGGGGGGGCCCCCGGGGGGAAAAAAAAAA
                 ||||||||||             
                 GGCCCCCGGG             
""",  # noqa: W291
        )
        psl = format(alignment, "psl")
        self.assertEqual(
            psl,
            """\
10	0	0	0	0	0	0	0	-	sequence	10	0	10	chromosome	40	17	27	1	10,	0,	17,
""",
        )


class TestComplex(unittest.TestCase):
    def setUp(self):
        aligner = PairwiseAligner()
        aligner.internal_open_gap_score = -1
        aligner.internal_extend_gap_score = -0.0
        aligner.match_score = +1
        aligner.mismatch_score = -1
        aligner.mode = "local"
        self.aligner = aligner

    def test1(self):
        aligner = self.aligner
        chromosome = Seq(
            "GCCTACCGTATAACAATGGTTATAATACAAGGCGGTCATAATTAAAGGGAGTGCAGCAACGGCCTGCTCTCCAAAAAAACAGGTTTTATGAAAAGAAAGTGCATTAACTGTTAAAGCCGTCATATCGGTGGGTTCTGCCAGTCACCGGCATACGTCCTGGGACAAAGACTTTTTACTACAATGCCAGGCGGGAGAGTCACCCGCCGCGGTGTCGACCCAGGGGACAGCGGGAAGATGTCGTGGTTTCCTTGTCATTAACCAACTCCATCTTAAAAGCTCCTCTAGCCATGGCATGGTACGTTGCGCGCACCCTTTTATCGGTAAGGCGCGGTGACTCTCTCCCAAAACAGTGCCATAATGGTTCGCTTCCTACCTAAGGCACTTACGGCCAATTAATGCGCAAGCGAGCGGAAGGTCTAACAGGGCACCGAATTCGATTA"
        )
        chromosome.id = "chromosome"
        transcript = Seq(
            "GGAATTTTAGCAGCCAAAGGACGGATCCTCCAAGGGGCCCCAGCACAGCACATTTTTAACGCGAACTAAGCGGGAGCGCATGTGGGACAGTTGATCCCATCCGCCTCAAAATTTCTCGCAATATCGGTTGGGGCACAGGTCCACTTTACGAATTCATACCGTGGTAGAGACCTTTATTAGATAGATATGACTGTTTGATTGCGGCATAGTACGACGAAGCAAGGGGATGGACGTTTCGGTTGCATTCGACCGGGTTGGGTCGAAAAACAGGTTTTATGAAAAGAAAGTGCATTAACTGTTAAAGCCGTCATATCGGTGGGTTC"
        )
        transcript.id = "transcript"
        sequence = Seq(
            "TCCAAGGGGCCCCAGCACAGCACATTTTTAACGCGGGGACAGTTGATCCCATCCGCCTTTTACGAATTCATACCGTGGTAGGCGGCATAGTACGACGAAGCGGTTGGGTCGAAAAACAGGTTGCCGTCATATCGGTGGGTTC"
        )
        sequence.id = "sequence"
        alignments1 = aligner.align(chromosome, transcript)
        alignment1 = alignments1[0]
        self.assertEqual(len(alignment1.path), 164)
        self.assertEqual(
            str(alignment1),
            """\
GCCTACCGTATAACAATGGTTATA------ATACAAGG-CGG----TCATAATTAAAGGGAGTG---CAGCAACGGCCTGCTCTCCAAAAAAACAGGTTTTATGAAAAGAAAGTGCATTAACTGTTAAAGC-----CGTCATATCGGTGG----GTTCTGCCAGTCACCGGCATACGTCCTGGGACAAAGACTTTTTACT-ACAATGCCAGGCGGGAGAGTCACCCGCCGCGGTGTCGACCCAGGGG-ACAGCGGGAAGATGTCGTGGTTTC-CTT---G---TCATTAACC-------A-ACTCCATCTTA--AAAGCTCCTCTAGCCATGGCATG---GT---ACGTTGCGCGCACCCTTTTA-T----CG--GTAAGG-------CG---CGGT-------GACTCTC--------TCCCAAAACAGTGCCATAATGGTTCGCTTCCTACCT-------AAG-GCACTT-ACGGCCAATTAATGCGCAAGCGAGCGGAAGGTC-TAACAG-GGCACCGAATTCGATTA
              |||--||-||------|---||||-|||----||------.|||||---|---|||||-----|.||-----------|||--||||-|------||.|.|----||||----||||-----||-|||----||||----||--||--|-|--||--|||.||-|||----||||-|---|||-||-.||||------------|-|---------||||-|-------||||-||||---------|||-------|-|||---|---||||--|||-------|-|--||-|-|||--|.|-------|||--||---|||---||---|--|||||-|||------||-|----||--|.||||-------||---||||-------|||---|--------||..||||||----------|||----||--||--|-------|||-|||-||-||.|----|||------||||---|-----|||-||.|.|-||----|--|||     
            GGAAT--TT-TAGCAGCCA---AAGGACGGATCCTC------CAAGGG---GCCCCAGCA-----CAGC-----------ACA--TTTT-T------AACGCG----AACT----AAGCGGGAGCG-CAT----GTGGGACAGT--TG--A-T--CC--CATCCG-CCT----CAAA-A---TTT-CTCGCAAT------------A-T---------CGGT-T-------GGGGCACAG---------GTC-------CACTTTACGAATTCAT--ACCGTGGTAGAGA--CC-T-TTATTAGA-------TAG--AT---ATGACTGTTTGA--TTGCG-GCA------TAGTACGACGAAGCAAGGGGATGGACGTTTCGGTTGCATTCGAC---CGGGTTGGGTCGAAAAACA----------GGT----TT--TA--TGAAAAGAAAGTGCA-TTAACTG----TTA------AAGC---C-----GTCATATCGGTGG----G--TTC     
""",  # noqa: W291
        )
        alignments2 = aligner.align(transcript, sequence)
        alignment2 = alignments2[0]
        self.assertEqual(len(alignment2.path), 12)
        self.assertEqual(
            str(alignment2),
            """\
GGAATTTTAGCAGCCAAAGGACGGATCCTCCAAGGGGCCCCAGCACAGCACATTTTTAACGCGAACTAAGCGGGAGCGCATGTGGGACAGTTGATCCCATCCGCCTCAAAATTTCTCGCAATATCGGTTGGGGCACAGGTCCACTTTACGAATTCATACCGTGGTAGAGACCTTTATTAGATAGATATGACTGTTTGATTGCGGCATAGTACGACGAAGCAAGGGGATGGACGTTTCGGTTGCATTCGACCGGGTTGGGTCGAAAAACAGGTTTTATGAAAAGAAAGTGCATTAACTGTTAAAGCCGTCATATCGGTGGGTTC
                            |||||||||||||||||||||||||||||||||||--------------------|||||||||||||||||||||||--------------------------------------|||||||||||||||||||||||---------------------------------||||||||||||||||||||--------------------------------|||||||||||||||||||||------------------------------||||||||||||||||||||
                            TCCAAGGGGCCCCAGCACAGCACATTTTTAACGCG--------------------GGGACAGTTGATCCCATCCGCCT--------------------------------------TTTACGAATTCATACCGTGGTAG---------------------------------GCGGCATAGTACGACGAAGC--------------------------------GGTTGGGTCGAAAAACAGGTT------------------------------GCCGTCATATCGGTGGGTTC
""",
        )
        alignment = alignment1.map(alignment2)
        self.assertEqual(len(alignment.path), 76)
        self.assertEqual(
            str(alignment),
            """\
GCCTACCGTATAACAATGGTTATAATACAAGGCGGTCATAATTAAAGGGAGTG---CAGCAACGGCCTGCTCTCCAAAAAAACAGGTTTTATGAAAAGAAAGTGCATTAACTGTTAAAGCCGTCATATCGGTGG----GTTCTGCCAGTCACCGGCATACGTCCTGGGACAAAGACTTTTTACTACAATGCCAGGCGGGAGAGTCACCCGCCGCGGTGTCGACCCAGGGGACAGCGGGAAGATGTCGTGGTTTCCTT---G---TCATTAACCAACTCCATCTTAAAAGCTCCTCTAGCCATGGCATGGTACGTT-------GCGCGCACCCTTTTA-T----CG--GTAAGGCGCGGTGACTCTC-------TCCCAAAACAGTGCCATAATGGTTCGCTTCCTACCTAAGGCACTTACGGCCAATTAATGCGCAAGCGAGCGGAAGGTC-TAACAG-GGCACCGAATTCGATTA
                                   ||------.|||||---|---|||||-----|.||-----------|||--||||-|------||.|.|----------------------------||----||--||--|-|--||--|||.||-|||------------------------------------------------------------------------------------------||---|---||||--|||-------------------------------------------------|||-|||------||-|----||--|.------------------------||..||||||----------|||----|------------------------------------||---|-----|||-||.|.|-||----|--|||     
                                   TC------CAAGGG---GCCCCAGCA-----CAGC-----------ACA--TTTT-T------AACGCG----------------------------GGGACAGT--TG--A-T--CC--CATCCG-CCT------------------------------------------------------------------------------------------TTTACGAATTCAT--ACC------------------------------------------GTGGTAGGCG-GCA------TAGTACGACGAAGC-----------------GGTTGGGTCGAAAAACA----------GGT----T------------------------------------GC---C-----GTCATATCGGTGG----G--TTC     
""",  # noqa: W291
        )
        psl = format(alignment, "psl")
        self.assertEqual(
            psl,
            """\
96	10	0	0	11	36	27	294	+	sequence	142	0	142	chromosome	440	35	435	37	2,6,1,5,4,3,4,1,6,2,2,2,1,1,2,6,3,2,1,4,3,3,3,2,1,2,2,10,3,1,2,1,3,6,2,1,3,	0,2,8,12,17,21,24,28,29,35,41,43,45,46,47,49,55,58,63,67,71,81,84,87,90,95,99,108,118,121,122,124,125,129,136,138,139,	35,43,52,53,63,78,83,88,95,129,131,135,139,141,144,148,155,248,250,251,257,302,306,315,317,318,320,339,359,366,403,408,414,417,423,429,432,
""",
        )

    def test2(self):
        aligner = self.aligner
        chromosome = Seq(
            "CTAATGCGCCTTGGTTTTGGCTTAACTAGAAGCAACCTGTAAGATTGCCAATTCTTCAGTCGAAGTAAATCTTCAATGTTTTGGACTCTTAGCGGATATGCGGCTGAGAAGTACGACATGTGTACATTCATACCTGCGTGACGGTCAGCCTCCCCCGGGACCTCATTGGGCGAATCTAGGTGTGATAATTGACACACTCTTGGTAAGAAGCACTCTTTACCCGATCTCCAAGTACCGACGCCAAGGCCAAGCTCTGCGATCTAAAGCTGCCGATCGTAGATCCAAGTCCTCAGCAAGCTCGCACGAATACGCAGTTCGAAGGCTGGGTGTTGTACGACGGTACGGTTGCTATAGCACTTTCGCGGTCTCGCTATTTTCAGTTTGACTCACCAGTCAGTATTGTCATCGACCAACTTGGAATAGTGTAACGCAGCGCTTGA"
        )
        chromosome.id = "chromosome"
        transcript = Seq(
            "CACCGGCGTCGGTACCAGAGGGCGTGAGTACCTTGTACTAGTACTCATTGGAATAATGCTCTTAGAAGTCATCTAAAAGTGACAACGCCTGTTTGGTTATGACGTTCACGACGCGTCTTAACAGACTAGCATTAGACCGACGGGTTGAGGCGTCTGGGTTGATACAGCCGTTTGCATCAGTGTATCTAACACTCTGAGGGATAATTGATGAACCGTGTTTTCCGATAGGTATGTACAGTACCACCACGCACGACTAAGGACCATTTTCTGCGTGCGACGGTTAAAATAACCTCAATCACT"
        )
        transcript.id = "transcript"
        sequence = Seq(
            "TCCCCTTCTAATGGAATCCCCCTCCGAAGGTCGCAGAAGCGGCCACGCCGGAGATACCAGTTCCACGCCTCAGGTTGGACTTGTCACACTTGTACGCGAT"
        )
        sequence.id = "sequence"
        alignments1 = aligner.align(chromosome, transcript)
        alignment1 = alignments1[0]
        self.assertEqual(len(alignment1.path), 126)
        self.assertEqual(
            str(alignment1),
            """\
CTAATGCGCCTTGGTTTTGGCTTAACTAGA-------AGCAACC-TGTAAGATTGCCAATTCTTCAGTCGAAGTAAATCTTCAATGTTTTGGA------CTCTTAG----CGGATATGCGGCTGAGAAGTACGACA-----TGT---GT----ACATTCATAC--CTGCGT-------GACGGTCAGCCT----CCCCCGGGACCTCATTG-GGCGAATCTAGGTGTGATA-A-----TTGACA-CA----CTCTTGGTAAGAAGCACTCT---------TTACCCGATCTCCAAGTACCGACGCCAAGGCCAAGCTCTG-----CGATCTAAAGCTGCCGATCGTAGATCCAAGTCCTCAGCAAGCTCGCACGAATACGCAG-------TTCGAAGGCTGGGTGTTGTACGACGGTACGGTTGCTATAGCACTTTCGCGGTCTCGCTATTTTCAGTTTGACTCACCAGTCAGTATTGTCATCGACCAACTTGGAATAGTGTAACGCAGCGCTTGA
     |||--|.|||---------||.|||-------||-.|||-||||------------||--------||||---|-|||-----|||||------|||||||----|----||----||-|.||||--||||-----|||---||----||.|||--||--|-||||-------|||--|-|||.|----||..||||------|||-||||--|||.|||-|||||-|-----|||-||-||----.|||-------||-||||||---------||----|||------|.||||-----------------||-----||||----||.|----||-|||----|-|||.|-||.||----|||||||.||---||-------||------|||.|||-----|||||||--------||-|--|----------------||----------|---|||--|||--------|||-----|||                         
CACCGGCG--TCGGT---------ACCAGAGGGCGTGAG-TACCTTGTA------------CT--------AGTA---C-TCA-----TTGGAATAATGCTCTTAGAAGTC----AT----CT-AAAAGT--GACAACGCCTGTTTGGTTATGACGTTC--ACGAC-GCGTCTTAACAGAC--T-AGCATTAGACCGACGGG------TTGAGGCG--TCTGGGT-TGATACAGCCGTTTG-CATCAGTGTATCT-------AA-CACTCTGAGGGATAATT----GAT------GAACCG-----------------TGTTTTCCGAT----AGGT----AT-GTA----C-AGTAC-CACCA----CGCACGACTA---AGGACCATTTT------CTGCGTG-----CGACGGT--------TA-A--A----------------AT----------A---ACC--TCA--------ATC-----ACT                         
""",  # noqa: W291
        )
        alignments2 = aligner.align(transcript, sequence)
        alignment2 = alignments2[0]
        self.assertEqual(len(alignment2.path), 66)
        self.assertEqual(
            str(alignment2),
            """\
CACCGGCGTCGGTACCAGAGGGCGTGAGTACCTTGTACTAGTACTCATTGGAATAATGCTCTTAGAAGTCATCTAAAAGTGACAACGCCTGTTTGGTTATGACGTTCACGACGCGTCTTAACAGACTAGCATTAGACCGACG--GGTTGAGGCGTCTGGGTTGATACAGCCGTTTGCATCAGTGTATCTAACA---CTCTGAGGGATAATTGATGAACCGTGTTTTCCGATAGGTATGTACAGTACCACCACGCACGACTAAGGACCATTTTCTG--CGTGCGACGGTTAAAATAACCTCAATCACT
        ||------------|-------||||---|||------|-||||||-------------------------------|.||-------------|--||-|||.|-|||---.||||--|||----|.||-|||--|---||------------|||||------------||||---||---||---|||--|||-----|||--||--|-----||----------------||---------|||-|||-------------||--|--||||                       
        TC------------C-------CCTT---CTA------A-TGGAAT-------------------------------CCCC-------------C--TC-CGAAG-GTC---GCAGA--AGC----GGCC-ACGCCG---GA------------GATAC------------CAGT---TC---CACGCCTC--AGG-----TTG--GA--C-----TT----------------GT---------CAC-ACT-------------TGTAC--GCGAT                      
""",  # noqa: W291
        )
        alignment = alignment1.map(alignment2)
        self.assertEqual(len(alignment.path), 78)
        self.assertEqual(
            str(alignment),
            """\
CTAATGCGCCTTGGTTTTGGCTTAACTAGAAGCAA-CC-TGTAAGATTGCCAATTCTTCAGTCGAAGTAAATCTTCAATGTTTTGGACTCTTAGCGGATATGCGGCTGAGAAGTACGACATGTGTA------CATTCATAC--CTGCGT----GACGGTCAGCCT--CCCCCG--GGACCTCATTGGGCGAATCTAGGTGT-GATAATTGACA-CAC--TCTTGGTAAGAAGCA---CTCT---TTACCCGATCTCCAAGTACCGACGCCAAGGCCAAGCTCTGCGATCTAAAGCTGCCGATCGTAGATCCAA--GTCCTCAGCAAGCTCGCACGAATACGCAGTTCGAAGGCTG--GGTGTTGTACGACGGTACGGTTGCTATAGCACTTTCGCGGTCTCGCTATTTTCAGTTTGACTCACCAGTCAGTATTGTCATCGACCAACTTGGAATAGTGTAACGCAGCGCTTGA
          |.------------------------||-|---------------||--------|----------|------||||---------------------------------------------|--||---|--.-|-||----||-----|||----||-.||--|---------|----------------||||--------||---||-----------||---|||----||----|--------|.--|---------------------------------------------------||--------------|||-|.|---------------||--.--|-----|||                                                                                                       
          TC-----------------------CCCTT---------------CT--------A----------A------TGGA---------------------------------------ATCCCCC--TC---CGAA-G-GTCGCAGA-----AGC--GGCC-ACGCCG---------G---------------AGATA-------CCA-GTTC-----------CACGCCTC-AGGTT----G--------GA--C-------------------------------------------------TTGT--------------CAC-ACT---------------TGTAC--G-----CGAT                                                                                                      
""",  # noqa: W291
        )
        psl = format(alignment, "psl")
        self.assertEqual(
            psl,
            """\
61	6	0	0	14	32	28	260	+	sequence	100	0	99	chromosome	440	10	337	35	2,2,1,2,1,1,4,1,2,1,1,1,2,2,3,2,3,1,1,4,2,2,2,3,2,1,2,1,2,3,3,2,1,1,3,	0,3,6,7,9,10,11,21,22,24,27,28,29,35,37,42,44,49,50,52,57,61,63,68,74,76,77,79,82,84,87,90,94,95,96,	10,35,37,53,63,74,81,124,127,132,133,135,137,139,146,151,154,157,167,183,194,197,210,212,216,222,231,235,285,301,305,323,325,328,334,
""",
        )


def map_check(alignment1, alignment2):
    psl1 = format(alignment1, "psl")
    handle = open("transcript.psl", "w")
    handle.write(psl1)
    handle.close()
    psl2 = format(alignment2, "psl")
    handle = open("sequence.psl", "w")
    handle.write(psl2)
    handle.close()
    stdout = os.popen("pslMap sequence.psl transcript.psl stdout")
    psl = stdout.read()
    os.remove("transcript.psl")
    os.remove("sequence.psl")
    return psl


def test_random(aligner, nBlocks1=1, nBlocks2=1, strand1="+", strand2="+"):
    chromosome = "".join(["ACGT"[random.randint(0, 3)] for i in range(1000)])
    nBlocks = nBlocks1
    transcript = ""
    position = 0
    for i in range(nBlocks):
        position += random.randint(60, 80)
        blockSize = random.randint(60, 80)
        transcript += chromosome[position : position + blockSize]
        position += blockSize
    nBlocks = nBlocks2
    sequence = ""
    position = 0
    for i in range(nBlocks):
        position += random.randint(20, 40)
        blockSize = random.randint(20, 40)
        sequence += transcript[position : position + blockSize]
        position += blockSize
    chromosome = Seq(chromosome)
    transcript = Seq(transcript)
    sequence = Seq(sequence)
    if strand1 == "-":
        chromosome = chromosome.reverse_complement()
    if strand2 == "-":
        sequence = sequence.reverse_complement()
    chromosome.id = "chromosome"
    transcript.id = "transcript"
    sequence.id = "sequence"
    alignments1 = aligner.align(chromosome, transcript, strand=strand1)
    alignment1 = alignments1[0]
    alignments2 = aligner.align(transcript, sequence, strand=strand2)
    alignment2 = alignments2[0]
    alignment = alignment1.map(alignment2)
    psl_check = map_check(alignment1, alignment2)
    psl = format(alignment, "psl")
    assert psl == psl_check
    print("Randomized test %d, %d, %s, %s OK" % (nBlocks1, nBlocks2, strand1, strand2))


def test_random_sequences(aligner, strand1="+", strand2="+"):
    chromosome = "".join(["ACGT"[random.randint(0, 3)] for i in range(1000)])
    transcript = "".join(["ACGT"[random.randint(0, 3)] for i in range(300)])
    sequence = "".join(["ACGT"[random.randint(0, 3)] for i in range(100)])
    chromosome = Seq(chromosome)
    transcript = Seq(transcript)
    sequence = Seq(sequence)
    chromosome.id = "chromosome"
    transcript.id = "transcript"
    sequence.id = "sequence"
    alignments = aligner.align(chromosome, transcript, strand=strand1)
    alignment1 = alignments[0]
    alignments = aligner.align(transcript, sequence, strand=strand2)
    alignment2 = alignments[0]
    psl_check = map_check(alignment1, alignment2)
    alignment = alignment1.map(alignment2)
    psl_check = psl_check.split()
    psl = format(alignment, "psl")
    psl = psl.split()
    assert psl[8:] == psl_check[8:]
    psl1 = format(alignment1, "psl")
    words = psl1.split()
    nBlocks1 = int(words[17])
    psl2 = format(alignment2, "psl")
    words = psl2.split()
    nBlocks2 = int(words[17])
    print(
        "Randomized sequence test %d, %d, %s, %s OK"
        % (nBlocks1, nBlocks2, strand1, strand2)
    )


def perform_randomized_tests(n=1000):
    """Perform randomized tests and compare to pslMap.

    Run this function to perform 8 x n mappings for alignments of randomly
    generated sequences, get the alignment in PSL format, and compare the
    result to that of pslMap.
    """
    aligner = PairwiseAligner()
    aligner.internal_open_gap_score = -1
    aligner.internal_extend_gap_score = -0.0
    aligner.match_score = +1
    aligner.mismatch_score = -1
    aligner.mode = "local"
    for i in range(n):
        nBlocks1 = random.randint(1, 10)
        nBlocks2 = random.randint(1, 10)
        test_random(aligner, nBlocks1, nBlocks2, "+", "+")
        test_random(aligner, nBlocks1, nBlocks2, "+", "-")
        test_random(aligner, nBlocks1, nBlocks2, "-", "+")
        test_random(aligner, nBlocks1, nBlocks2, "-", "-")
        test_random_sequences("+", "+")
        test_random_sequences("+", "-")
        test_random_sequences("-", "+")
        test_random_sequences("-", "-")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
