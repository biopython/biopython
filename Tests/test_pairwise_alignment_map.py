# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


"""Tests for mapping pairwise alignments."""

import os
import random
import unittest

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align.Alignment.map."
    ) from None

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
        self.assertTrue(
            numpy.array_equal(alignment1.coordinates, numpy.array([[12, 31], [0, 19]]))
        )
        self.assertEqual(
            str(alignment1),
            """\
chromosom        12 GGGGGGGCCCCCGGGGGGA 31
                  0 ||||||||||||||||||| 19
transcrip         0 GGGGGGGCCCCCGGGGGGA 19
""",
        )
        alignments2 = aligner.align(transcript, sequence)
        self.assertEqual(len(alignments2), 1)
        alignment2 = alignments2[0]
        self.assertTrue(
            numpy.array_equal(alignment2.coordinates, numpy.array([[5, 15], [0, 10]]))
        )
        self.assertEqual(
            str(alignment2),
            """\
transcrip         5 GGCCCCCGGG 15
                  0 |||||||||| 10
sequence          0 GGCCCCCGGG 10
""",
        )
        alignment = alignment1.map(alignment2)
        self.assertTrue(
            numpy.array_equal(alignment.coordinates, numpy.array([[17, 27], [0, 10]]))
        )
        self.assertEqual(
            str(alignment),
            """\
chromosom        17 GGCCCCCGGG 27
                  0 |||||||||| 10
sequence          0 GGCCCCCGGG 10
""",
        )
        line = format(alignment, "psl")
        self.assertEqual(
            line,
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
chromosom         0 GGGCCCCCGGGGGGA 15
                  0 ||||||||||||||| 15
transcrip         3 GGGCCCCCGGGGGGA 18
""",
        )
        alignments2 = aligner.align(transcript, sequence)
        self.assertEqual(len(alignments2), 1)
        alignment2 = alignments2[0]
        self.assertEqual(
            str(alignment2),
            """\
transcrip         1 GGGGGCCCCCGGG 14
                  0 ||||||||||||| 13
sequence          0 GGGGGCCCCCGGG 13
""",
        )
        alignment = alignment1.map(alignment2)
        self.assertTrue(
            numpy.array_equal(alignment.coordinates, numpy.array([[0, 11], [2, 13]]))
        )
        self.assertEqual(
            str(alignment),
            """\
chromosom         0 GGGCCCCCGGG 11
                  0 ||||||||||| 11
sequence          2 GGGCCCCCGGG 13
""",
        )
        line = format(alignment, "psl")
        self.assertEqual(
            line,
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
chromosom        12 GGGGGGGCCCCCGGG 27
                  0 ||||||||||||||| 15
transcrip         0 GGGGGGGCCCCCGGG 15
""",
        )
        alignments2 = aligner.align(transcript, sequence)
        self.assertEqual(len(alignments2), 1)
        alignment2 = alignments2[0]
        self.assertEqual(
            str(alignment2),
            """\
transcrip         5 GGCCCCCGGGGG 17
                  0 |||||||||||| 12
sequence          0 GGCCCCCGGGGG 12
""",
        )
        alignment = alignment1.map(alignment2)
        self.assertTrue(
            numpy.array_equal(alignment.coordinates, numpy.array([[17, 27], [0, 10]]))
        )
        self.assertEqual(
            str(alignment),
            """\
chromosom        17 GGCCCCCGGG 27
                  0 |||||||||| 10
sequence          0 GGCCCCCGGG 10
""",
        )
        line = format(alignment, "psl")
        self.assertEqual(
            line,
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
        self.assertTrue(
            numpy.array_equal(alignment1.coordinates, numpy.array([[12, 31], [19, 0]]))
        )
        self.assertEqual(
            str(alignment1),
            """\
chromosom        12 GGGGGGGCCCCCGGGGGGA 31
                  0 ||||||||||||||||||| 19
transcrip        19 GGGGGGGCCCCCGGGGGGA  0
""",
        )
        alignments2 = aligner.align(transcript, sequence, strand="-")
        self.assertEqual(len(alignments2), 1)
        alignment2 = alignments2[0]
        self.assertTrue(
            numpy.array_equal(alignment2.coordinates, numpy.array([[4, 14], [10, 0]]))
        )
        self.assertEqual(
            str(alignment2),
            """\
transcrip         4 CCCGGGGGCC 14
                  0 |||||||||| 10
sequence         10 CCCGGGGGCC  0
""",
        )
        alignment = alignment1.map(alignment2)
        self.assertTrue(
            numpy.array_equal(alignment.coordinates, numpy.array([[17, 27], [0, 10]]))
        )
        self.assertEqual(
            str(alignment),
            """\
chromosom        17 GGCCCCCGGG 27
                  0 |||||||||| 10
sequence          0 GGCCCCCGGG 10
""",
        )
        line = format(alignment, "psl")
        self.assertEqual(
            line,
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
        self.assertTrue(
            numpy.array_equal(alignment1.coordinates, numpy.array([[12, 31], [0, 19]]))
        )
        self.assertEqual(
            str(alignment1),
            """\
chromosom        12 GGGGGGGCCCCCGGGGGGA 31
                  0 ||||||||||||||||||| 19
transcrip         0 GGGGGGGCCCCCGGGGGGA 19
""",
        )
        alignments2 = aligner.align(transcript, sequence, "-")
        self.assertEqual(len(alignments2), 1)
        alignment2 = alignments2[0]
        self.assertTrue(
            numpy.array_equal(alignment2.coordinates, numpy.array([[5, 15], [10, 0]]))
        )
        self.assertEqual(
            str(alignment2),
            """\
transcrip         5 GGCCCCCGGG 15
                  0 |||||||||| 10
sequence         10 GGCCCCCGGG  0
""",
        )
        alignment = alignment1.map(alignment2)
        self.assertTrue(
            numpy.array_equal(alignment.coordinates, numpy.array([[17, 27], [10, 0]]))
        )
        self.assertEqual(
            str(alignment),
            """\
chromosom        17 GGCCCCCGGG 27
                  0 |||||||||| 10
sequence         10 GGCCCCCGGG  0
""",
        )
        line = format(alignment, "psl")
        self.assertEqual(
            line,
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
        self.assertTrue(
            numpy.array_equal(alignment1.coordinates, numpy.array([[12, 31], [19, 0]]))
        )
        self.assertEqual(
            str(alignment1),
            """\
chromosom        12 GGGGGGGCCCCCGGGGGGA 31
                  0 ||||||||||||||||||| 19
transcrip        19 GGGGGGGCCCCCGGGGGGA  0
""",
        )
        alignments2 = aligner.align(transcript, sequence)
        self.assertEqual(len(alignments2), 1)
        alignment2 = alignments2[0]
        self.assertTrue(
            numpy.array_equal(alignment2.coordinates, numpy.array([[4, 14], [0, 10]]))
        )
        self.assertEqual(
            str(alignment2),
            """\
transcrip         4 CCCGGGGGCC 14
                  0 |||||||||| 10
sequence          0 CCCGGGGGCC 10
""",
        )
        alignment = alignment1.map(alignment2)
        self.assertTrue(
            numpy.array_equal(alignment.coordinates, numpy.array([[17, 27], [10, 0]]))
        )
        self.assertEqual(
            str(alignment),
            """\
chromosom        17 GGCCCCCGGG 27
                  0 |||||||||| 10
sequence         10 GGCCCCCGGG  0
""",
        )
        line = format(alignment, "psl")
        self.assertEqual(
            line,
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
        self.assertEqual(alignment1.coordinates.shape[1], 164)
        self.assertEqual(
            str(alignment1),
            """\
chromosom        14 AATGGTTATA------ATACAAGG-CGG----TCATAATTAAAGGGAGTG---CAGCAAC
                  0 |||--||-||------|---||||-|||----||------.|||||---|---|||||--
transcrip         2 AAT--TT-TAGCAGCCA---AAGGACGGATCCTC------CAAGGG---GCCCCAGCA--

chromosom        60 GGCCTGCTCTCCAAAAAAACAGGTTTTATGAAAAGAAAGTGCATTAACTGTTAAAGC---
                 60 ---|.||-----------|||--||||-|------||.|.|----||||----||||---
transcrip        45 ---CAGC-----------ACA--TTTT-T------AACGCG----AACT----AAGCGGG

chromosom       117 --CGTCATATCGGTGG----GTTCTGCCAGTCACCGGCATACGTCCTGGGACAAAGACTT
                120 --||-|||----||||----||--||--|-|--||--|||.||-|||----||||-|---
transcrip        74 AGCG-CAT----GTGGGACAGT--TG--A-T--CC--CATCCG-CCT----CAAA-A---

chromosom       171 TTTACT-ACAATGCCAGGCGGGAGAGTCACCCGCCGCGGTGTCGACCCAGGGG-ACAGCG
                180 |||-||-.||||------------|-|---------||||-|-------||||-||||--
transcrip       111 TTT-CTCGCAAT------------A-T---------CGGT-T-------GGGGCACAG--

chromosom       229 GGAAGATGTCGTGGTTTC-CTT---G---TCATTAACC-------A-ACTCCATCTTA--
                240 -------|||-------|-|||---|---||||--|||-------|-|--||-|-|||--
transcrip       138 -------GTC-------CACTTTACGAATTCAT--ACCGTGGTAGAGA--CC-T-TTATT

chromosom       272 AAAGCTCCTCTAGCCATGGCATG---GT---ACGTTGCGCGCACCCTTTTA-T----CG-
                300 |.|-------|||--||---|||---||---|--|||||-|||------||-|----||-
transcrip       178 AGA-------TAG--AT---ATGACTGTTTGA--TTGCG-GCA------TAGTACGACGA

chromosom       320 -GTAAGG-------CG---CGGT-------GACTCTC--------TCCCAAAACAGTGCC
                360 -|.||||-------||---||||-------|||---|--------||..||||||-----
transcrip       217 AGCAAGGGGATGGACGTTTCGGTTGCATTCGAC---CGGGTTGGGTCGAAAAACA-----

chromosom       354 ATAATGGTTCGCTTCCTACCT-------AAG-GCACTT-ACGGCCAATTAATGCGCAAGC
                420 -----|||----||--||--|-------|||-|||-||-||.|----|||------||||
transcrip       269 -----GGT----TT--TA--TGAAAAGAAAGTGCA-TTAACTG----TTA------AAGC

chromosom       405 GAGCGGAAGGTC-TAACAG-GGCACCGAATTC 435
                480 ---|-----|||-||.|.|-||----|--||| 512
transcrip       305 ---C-----GTCATATCGGTGG----G--TTC 323
""",
        )
        alignments2 = aligner.align(transcript, sequence)
        alignment2 = alignments2[0]
        self.assertEqual(alignment2.coordinates.shape[1], 12)
        self.assertEqual(
            str(alignment2),
            """\
transcrip        28 TCCAAGGGGCCCCAGCACAGCACATTTTTAACGCGAACTAAGCGGGAGCGCATGTGGGAC
                  0 |||||||||||||||||||||||||||||||||||--------------------|||||
sequence          0 TCCAAGGGGCCCCAGCACAGCACATTTTTAACGCG--------------------GGGAC

transcrip        88 AGTTGATCCCATCCGCCTCAAAATTTCTCGCAATATCGGTTGGGGCACAGGTCCACTTTA
                 60 ||||||||||||||||||--------------------------------------||||
sequence         40 AGTTGATCCCATCCGCCT--------------------------------------TTTA

transcrip       148 CGAATTCATACCGTGGTAGAGACCTTTATTAGATAGATATGACTGTTTGATTGCGGCATA
                120 |||||||||||||||||||---------------------------------||||||||
sequence         62 CGAATTCATACCGTGGTAG---------------------------------GCGGCATA

transcrip       208 GTACGACGAAGCAAGGGGATGGACGTTTCGGTTGCATTCGACCGGGTTGGGTCGAAAAAC
                180 ||||||||||||--------------------------------||||||||||||||||
sequence         89 GTACGACGAAGC--------------------------------GGTTGGGTCGAAAAAC

transcrip       268 AGGTTTTATGAAAAGAAAGTGCATTAACTGTTAAAGCCGTCATATCGGTGGGTTC 323
                240 |||||------------------------------|||||||||||||||||||| 295
sequence        117 AGGTT------------------------------GCCGTCATATCGGTGGGTTC 142
""",
        )
        alignment = alignment1.map(alignment2)
        self.assertEqual(alignment.coordinates.shape[1], 76)
        self.assertEqual(
            str(alignment),
            """\
chromosom        35 TCATAATTAAAGGGAGTG---CAGCAACGGCCTGCTCTCCAAAAAAACAGGTTTTATGAA
                  0 ||------.|||||---|---|||||-----|.||-----------|||--||||-|---
sequence          0 TC------CAAGGG---GCCCCAGCA-----CAGC-----------ACA--TTTT-T---

chromosom        92 AAGAAAGTGCATTAACTGTTAAAGCCGTCATATCGGTGG----GTTCTGCCAGTCACCGG
                 60 ---||.|.|----------------------------||----||--||--|-|--||--
sequence         29 ---AACGCG----------------------------GGGACAGT--TG--A-T--CC--

chromosom       148 CATACGTCCTGGGACAAAGACTTTTTACTACAATGCCAGGCGGGAGAGTCACCCGCCGCG
                120 |||.||-|||--------------------------------------------------
sequence         49 CATCCG-CCT--------------------------------------------------

chromosom       208 GTGTCGACCCAGGGGACAGCGGGAAGATGTCGTGGTTTCCTT---G---TCATTAACCAA
                180 ----------------------------------------||---|---||||--|||--
sequence         58 ----------------------------------------TTTACGAATTCAT--ACC--

chromosom       262 CTCCATCTTAAAAGCTCCTCTAGCCATGGCATGGTACGTT-------GCGCGCACCCTTT
                240 -----------------------------------------------|||-|||------
sequence         74 ----------------------------------------GTGGTAGGCG-GCA------

chromosom       315 TA-T----CG--GTAAGGCGCGGTGACTCTC-------TCCCAAAACAGTGCCATAATGG
                300 ||-|----||--|.------------------------||..||||||----------||
sequence         87 TAGTACGACGAAGC-----------------GGTTGGGTCGAAAAACA----------GG

chromosom       361 TTCGCTTCCTACCTAAGGCACTTACGGCCAATTAATGCGCAAGCGAGCGGAAGGTC-TAA
                360 |----|------------------------------------||---|-----|||-||.
sequence        120 T----T------------------------------------GC---C-----GTCATAT

chromosom       420 CAG-GGCACCGAATTC 435
                420 |.|-||----|--||| 436
sequence        132 CGGTGG----G--TTC 142
""",
        )
        line = format(alignment, "psl")
        self.assertEqual(
            line,
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
        self.assertEqual(alignment1.coordinates.shape[1], 126)
        self.assertEqual(
            str(alignment1),
            """\
chromosom         5 GCGCCTTGGTTTTGGCTTAACTAGA-------AGCAACC-TGTAAGATTGCCAATTCTTC
                  0 |||--|.|||---------||.|||-------||-.|||-||||------------||--
transcrip         5 GCG--TCGGT---------ACCAGAGGGCGTGAG-TACCTTGTA------------CT--

chromosom        57 AGTCGAAGTAAATCTTCAATGTTTTGGA------CTCTTAG----CGGATATGCGGCTGA
                 60 ------||||---|-|||-----|||||------|||||||----|----||----||-|
transcrip        39 ------AGTA---C-TCA-----TTGGAATAATGCTCTTAGAAGTC----AT----CT-A

chromosom       107 GAAGTACGACA-----TGT---GT----ACATTCATAC--CTGCGT-------GACGGTC
                120 .||||--||||-----|||---||----||.|||--||--|-||||-------|||--|-
transcrip        75 AAAGT--GACAACGCCTGTTTGGTTATGACGTTC--ACGAC-GCGTCTTAACAGAC--T-

chromosom       146 AGCCT----CCCCCGGGACCTCATTG-GGCGAATCTAGGTGTGATA-A-----TTGACA-
                180 |||.|----||..||||------|||-||||--|||.|||-|||||-|-----|||-||-
transcrip       127 AGCATTAGACCGACGGG------TTGAGGCG--TCTGGGT-TGATACAGCCGTTTG-CAT

chromosom       194 CA----CTCTTGGTAAGAAGCACTCT---------TTACCCGATCTCCAAGTACCGACGC
                240 ||----.|||-------||-||||||---------||----|||------|.||||----
transcrip       177 CAGTGTATCT-------AA-CACTCTGAGGGATAATT----GAT------GAACCG----

chromosom       241 CAAGGCCAAGCTCTG-----CGATCTAAAGCTGCCGATCGTAGATCCAAGTCCTCAGCAA
                300 -------------||-----||||----||.|----||-|||----|-|||.|-||.||-
transcrip       215 -------------TGTTTTCCGAT----AGGT----AT-GTA----C-AGTAC-CACCA-

chromosom       296 GCTCGCACGAATACGCAG-------TTCGAAGGCTGGGTGTTGTACGACGGTACGGTTGC
                360 ---|||||||.||---||-------||------|||.|||-----|||||||--------
transcrip       246 ---CGCACGACTA---AGGACCATTTT------CTGCGTG-----CGACGGT--------

chromosom       349 TATAGCACTTTCGCGGTCTCGCTATTTTCAGTTTGACTCACCAGTCAGTATTGTCATCGA
                420 ||-|--|----------------||----------|---|||--|||--------|||--
transcrip       281 TA-A--A----------------AT----------A---ACC--TCA--------ATC--

chromosom       409 CCAACT 415
                480 ---||| 486
transcrip       297 ---ACT 300
""",
        )
        alignments2 = aligner.align(transcript, sequence)
        alignment2 = alignments2[0]
        self.assertEqual(alignment2.coordinates.shape[1], 66)
        self.assertEqual(
            str(alignment2),
            """\
transcrip         8 TCGGTACCAGAGGGCGTGAGTACCTTGTACTAGTACTCATTGGAATAATGCTCTTAGAAG
                  0 ||------------|-------||||---|||------|-||||||--------------
sequence          0 TC------------C-------CCTT---CTA------A-TGGAAT--------------

transcrip        68 TCATCTAAAAGTGACAACGCCTGTTTGGTTATGACGTTCACGACGCGTCTTAACAGACTA
                 60 -----------------|.||-------------|--||-|||.|-|||---.||||--|
sequence         17 -----------------CCCC-------------C--TC-CGAAG-GTC---GCAGA--A

transcrip       128 GCATTAGACCGACG--GGTTGAGGCGTCTGGGTTGATACAGCCGTTTGCATCAGTGTATC
                120 ||----|.||-|||--|---||------------|||||------------||||---||
sequence         38 GC----GGCC-ACGCCG---GA------------GATAC------------CAGT---TC

transcrip       186 TAACA---CTCTGAGGGATAATTGATGAACCGTGTTTTCCGATAGGTATGTACAGTACCA
                180 ---||---|||--|||-----|||--||--|-----||----------------||----
sequence         63 ---CACGCCTC--AGG-----TTG--GA--C-----TT----------------GT----

transcrip       243 CCACGCACGACTAAGGACCATTTTCTG--CGTGCGA 277
                240 -----|||-|||-------------||--|--|||| 276
sequence         84 -----CAC-ACT-------------TGTAC--GCGA  99
""",
        )
        alignment = alignment1.map(alignment2)
        self.assertEqual(alignment.coordinates.shape[1], 78)
        self.assertEqual(
            str(alignment),
            """\
chromosom        10 TTGGTTTTGGCTTAACTAGAAGCAA-CC-TGTAAGATTGCCAATTCTTCAGTCGAAGTAA
                  0 |.------------------------||-|---------------||--------|----
sequence          0 TC-----------------------CCCTT---------------CT--------A----

chromosom        68 ATCTTCAATGTTTTGGACTCTTAGCGGATATGCGGCTGAGAAGTACGACATGTGTA----
                 60 ------|------||||-------------------------------------------
sequence         10 ------A------TGGA---------------------------------------ATCC

chromosom       124 --CATTCATAC--CTGCGT----GACGGTCAGCCT--CCCCCG--GGACCTCATTGGGCG
                120 --|--||---|--.-|-||----||-----|||----||-.||--|---------|----
sequence         19 CCC--TC---CGAA-G-GTCGCAGA-----AGC--GGCC-ACGCCG---------G----

chromosom       172 AATCTAGGTGT-GATAATTGACA-CAC--TCTTGGTAAGAAGCA---CTCT---TTACCC
                180 ------------||||--------||---||-----------||---|||----||----
sequence         51 -----------AGATA-------CCA-GTTC-----------CACGCCTC-AGGTT----

chromosom       222 GATCTCCAAGTACCGACGCCAAGGCCAAGCTCTGCGATCTAAAGCTGCCGATCGTAGATC
                240 |--------|.--|----------------------------------------------
sequence         76 G--------GA--C----------------------------------------------

chromosom       282 CAA--GTCCTCAGCAAGCTCGCACGAATACGCAGTTCGAAGGCTG--GGTGTTGTACGA
                300 -----||--------------|||-|.|---------------||--.--|-----|||
sequence         80 ---TTGT--------------CAC-ACT---------------TGTAC--G-----CGA

chromosom       337
                359
sequence         99
""",
        )
        line = format(alignment, "psl")
        self.assertEqual(
            line,
            """\
61	6	0	0	14	32	28	260	+	sequence	100	0	99	chromosome	440	10	337	35	2,2,1,2,1,1,4,1,2,1,1,1,2,2,3,2,3,1,1,4,2,2,2,3,2,1,2,1,2,3,3,2,1,1,3,	0,3,6,7,9,10,11,21,22,24,27,28,29,35,37,42,44,49,50,52,57,61,63,68,74,76,77,79,82,84,87,90,94,95,96,	10,35,37,53,63,74,81,124,127,132,133,135,137,139,146,151,154,157,167,183,194,197,210,212,216,222,231,235,285,301,305,323,325,328,334,
""",
        )


def map_check(alignment1, alignment2):
    line1 = format(alignment1, "psl")
    handle = open("transcript.psl", "w")
    handle.write(line1)
    handle.close()
    line2 = format(alignment2, "psl")
    handle = open("sequence.psl", "w")
    handle.write(line2)
    handle.close()
    stdout = os.popen("pslMap sequence.psl transcript.psl stdout")
    line = stdout.read()
    os.remove("transcript.psl")
    os.remove("sequence.psl")
    return line


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
    line_check = map_check(alignment1, alignment2)
    line = format(alignment, "psl")
    assert line == line_check
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
    line_check = map_check(alignment1, alignment2)
    alignment = alignment1.map(alignment2)
    line_check = line_check.split()
    line = format(alignment, "psl")
    line = line.split()
    assert line[8:] == line_check[8:]
    line1 = format(alignment1, "psl")
    words = line1.split()
    nBlocks1 = int(words[17])
    line2 = format(alignment2, "psl")
    words = line2.split()
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
        test_random_sequences(aligner, "+", "+")
        test_random_sequences(aligner, "+", "-")
        test_random_sequences(aligner, "-", "+")
        test_random_sequences(aligner, "-", "-")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
