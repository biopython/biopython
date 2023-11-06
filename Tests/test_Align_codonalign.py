# Copyright (C) 2013 by Zheng Ruan (zruan1991@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for CodonAligner and Bio.Align.analysis."""

import unittest

import numpy as np

from Bio import SeqIO
from Bio import Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable
from Bio.Align import Alignment
from Bio.Align import CodonAligner
from Bio.Align.analysis import calculate_dn_ds, calculate_dn_ds_matrix, mktest


class TestBasic(unittest.TestCase):
    def test_aligner(self):
        aligner = CodonAligner()
        self.assertEqual(
            str(aligner),
            """\
Codon aligner with parameters
  wildcard: 'X'
  match_score: 1.0
  mismatch_score: 0.0
  frameshift_minus_two_score: -3.0
  frameshift_minus_one_score: -3.0
  frameshift_plus_one_score: -3.0
  frameshift_plus_two_score: -3.0
""",
        )
        aligner.wildcard = "Y"
        aligner.match_score = 2.0
        aligner.mismatch_score = -1.0
        aligner.frameshift_score = -5.0
        aligner.frameshift_two_score = -2.0
        aligner.frameshift_minus_score = -4.0
        self.assertEqual(
            str(aligner),
            """\
Codon aligner with parameters
  wildcard: 'Y'
  match_score: 2.0
  mismatch_score: -1.0
  frameshift_minus_two_score: -4.0
  frameshift_minus_one_score: -4.0
  frameshift_plus_one_score: -5.0
  frameshift_plus_two_score: -2.0
""",
        )
        self.assertEqual(aligner.wildcard, "Y")
        self.assertAlmostEqual(aligner.match_score, 2.0)
        self.assertAlmostEqual(aligner.mismatch_score, -1.0)
        self.assertAlmostEqual(aligner.frameshift_minus_two_score, -4.0)
        self.assertAlmostEqual(aligner.frameshift_minus_one_score, -4.0)
        self.assertAlmostEqual(aligner.frameshift_plus_one_score, -5.0)
        self.assertAlmostEqual(aligner.frameshift_plus_two_score, -2.0)

    def test_alignments(self):
        aligner = CodonAligner()
        aligner.frameshift_score = -1.0
        dna = SeqRecord(Seq("TTTAAAAAAAAATTT"), id="dna")
        pro = SeqRecord(Seq("FKKKF"), id="pro")
        alignments = aligner.align(pro, dna)
        self.assertEqual(len(alignments), 1)
        self.assertAlmostEqual(alignments.score, 5.0)
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
pro               0 F  K  K  K  F    5
dna               0 TTTAAAAAAAAATTT 15
""",
        )
        dna = SeqRecord(Seq("TTTAAAAAAAATTT"), id="dna")
        score = aligner.score(pro, dna)
        self.assertAlmostEqual(score, 4.0)
        alignments = aligner.align(pro, dna)
        self.assertEqual(len(alignments), 2)
        self.assertAlmostEqual(alignments.score, 4.0)
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
pro               0 F  K  K   3
dna               0 TTTAAAAAA 9

pro               3 K  F    5
dna               8 AAATTT 14
""",
        )
        alignment = alignments[1]
        self.assertEqual(
            str(alignment),
            """\
pro               0 F  K   2
dna               0 TTTAAA 6

pro               2 K  K  F    5
dna               5 AAAAAATTT 14
""",
        )
        dna = SeqRecord(Seq("TTTAAAAAAATTT"), id="dna")
        score = aligner.score(pro, dna)
        self.assertAlmostEqual(score, 4.0)
        alignments = aligner.align(pro, dna)
        self.assertEqual(len(alignments), 2)
        self.assertAlmostEqual(alignments.score, 4.0)
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
pro               0 F  K  K   3
dna               0 TTTAAAAAA 9

pro               3 K  F    5
dna               7 AAATTT 13
""",
        )
        alignment = alignments[1]
        self.assertEqual(
            str(alignment),
            """\
pro               0 F  K   2
dna               0 TTTAAA 6

pro               2 K  K  F    5
dna               4 AAAAAATTT 13
""",
        )
        dna = SeqRecord(Seq("TTTAAAAAATTT"), id="dna")
        score = aligner.score(pro, dna)
        self.assertAlmostEqual(score, 3.0)
        alignments = aligner.align(pro, dna)
        self.assertEqual(len(alignments), 2)
        self.assertAlmostEqual(alignments.score, 3.0)
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
pro               0 F  K   2
dna               0 TTTAAA 6

pro               2 K   3
dna               5 AAA 8

pro               3 K  F    5
dna               6 AAATTT 12
""",
        )
        alignment = alignments[1]
        self.assertEqual(
            str(alignment),
            """\
pro               0 F  K   2
dna               0 TTTAAA 6

pro               2 K   3
dna               4 AAA 7

pro               3 K  F    5
dna               6 AAATTT 12
""",
        )
        dna = SeqRecord(Seq("TTTAAAAATTT"), id="dna")
        score = aligner.score(pro, dna)
        self.assertAlmostEqual(score, 3.0)
        alignments = aligner.align(pro, dna)
        self.assertEqual(len(alignments), 1)
        self.assertAlmostEqual(alignments.score, 3.0)
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
pro               0 F  K   2
dna               0 TTTAAA 6

pro               2 K   3
dna               4 AAA 7

pro               3 K  F    5
dna               5 AAATTT 11
""",
        )
        dna = SeqRecord(Seq("TTTAAAAAAAAAATTT"), id="dna")
        alignments = aligner.align(pro, dna)
        self.assertEqual(len(alignments), 6)
        self.assertAlmostEqual(alignments.score, 4.0)
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
pro               0 F  K  K  K  F    5
dna               0 TTTAAAAAAAAAATT 15
""",
        )
        alignment = alignments[1]
        self.assertEqual(
            str(alignment),
            """\
pro               0 F  K  K  K  F    5
dna               1 TTAAAAAAAAAATTT 16
""",
        )
        alignment = alignments[2]
        self.assertEqual(
            str(alignment),
            """\
pro               0 F  -K  K  K  F    5
dna               0 TTTAAAAAAAAAATTT 16
""",
        )
        alignment = alignments[3]
        self.assertEqual(
            str(alignment),
            """\
pro               0 F  K  -K  K  F    5
dna               0 TTTAAAAAAAAAATTT 16
""",
        )
        alignment = alignments[4]
        self.assertEqual(
            str(alignment),
            """\
pro               0 F  K  K  -K  F    5
dna               0 TTTAAAAAAAAAATTT 16
""",
        )
        alignment = alignments[5]
        self.assertEqual(
            str(alignment),
            """\
pro               0 F  K  K  K  -F    5
dna               0 TTTAAAAAAAAAATTT 16
""",
        )
        dna = SeqRecord(Seq("TTTAAAAAAAAAAATTT"), id="dna")
        alignments = aligner.align(pro, dna)
        self.assertEqual(len(alignments), 6)
        self.assertAlmostEqual(alignments.score, 4.0)
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
pro               0 F  K  K  K  F    5
dna               0 TTTAAAAAAAAAAAT 15
""",
        )
        alignment = alignments[1]
        self.assertEqual(
            str(alignment),
            """\
pro               0 F  K  K  K  F    5
dna               2 TAAAAAAAAAAATTT 17
""",
        )
        alignment = alignments[2]
        self.assertEqual(
            str(alignment),
            """\
pro               0 F  --K  K  K  F    5
dna               0 TTTAAAAAAAAAAATTT 17
""",
        )
        alignment = alignments[3]
        self.assertEqual(
            str(alignment),
            """\
pro               0 F  K  --K  K  F    5
dna               0 TTTAAAAAAAAAAATTT 17
""",
        )
        alignment = alignments[4]
        self.assertEqual(
            str(alignment),
            """\
pro               0 F  K  K  --K  F    5
dna               0 TTTAAAAAAAAAAATTT 17
""",
        )
        alignment = alignments[5]
        self.assertEqual(
            str(alignment),
            """\
pro               0 F  K  K  K  --F    5
dna               0 TTTAAAAAAAAAAATTT 17
""",
        )
        dna = SeqRecord(Seq("TTTAAAAAAAAAAAATTT"), id="dna")
        alignments = aligner.align(pro, dna)
        self.assertEqual(len(alignments), 2)
        self.assertAlmostEqual(alignments.score, 4.0)
        alignment = alignments[0]
        self.assertEqual(
            str(alignment),
            """\
pro               0 F  K  K  K  F    5
dna               0 TTTAAAAAAAAAAAA 15
""",
        )
        alignment = alignments[1]
        self.assertEqual(
            str(alignment),
            """\
pro               0 F  K  K  K  F    5
dna               3 AAAAAAAAAAAATTT 18
""",
        )


class TestBuildAndIO(unittest.TestCase):
    def test1(self):
        aligner = CodonAligner()
        nucleotide_records = SeqIO.parse("codonalign/nucl1.fa", "fasta")
        protein_alignment = Align.read("codonalign/pro1.aln", "clustal")
        self.assertEqual(len(protein_alignment.sequences), 3)
        codon_alignments = []
        nucleotide_record = next(nucleotide_records)
        protein_record = protein_alignment.sequences[0]
        self.assertEqual(nucleotide_record.id, protein_record.id)
        alignments = aligner.align(protein_record, nucleotide_record)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 183], [0, 549]]))
        )
        self.assertEqual(
            str(alignment),
            """\
isotig697         0 R  G  D  Q  R  S  N  F  Q  L  S  P  S  T  M  Q  I  S  T  G  
isotig697         0 AGAGGCGATCAACGCAGCAACTTCCAGCTGTCTCCCTCCACCATGCAGATCTCCACAGGG

isotig697        20 L  L  C  L  L  L  V  A  T  G  F  T  S  Q  V  L  A  H  P  G  
isotig697        60 CTTCTGTGCcTGCTGCTTGTGGCCACTGGCTTCACTTCCCAGGTGCTGGCTCACCCAGGC

isotig697        40 S  I  P  S  T  Y  C  F  V  M  T  S  K  K  I  P  K  S  L  L  
isotig697       120 TCTATCCCATCTACCTaCTGCTTTGTTATGACCAGTAAGAaGATCCCCAAATCACTACTG

isotig697        60 K  S  Y  K  R  I  S  N  S  R  C  T  L  K  A  I  L  F  K  T  
isotig697       180 AaGAGCTACAAAaGAATCTCCAACAGCAGATGCACCcTGAAAGCCATACTCTTCAAGACC

isotig697        80 K  S  G  K  E  I  C  A  D  P  K  K  K  W  V  Q  D  A  T  K  
isotig697       240 AAGTCGGGCAAAGAGATCTGTGCTGACCCCAAGAAGAAGTGGGTCcAGGATGCCACAAAG

isotig697       100 H  L  D  Q  I  L  Q  T  P  K  P  T  I  P  S  F  E  T  H  P  
isotig697       300 CACCTGGACCAAATCCTTCAAACTCCAAAACCGACAATCCCCTCTTTTGAGACTCACCCA

isotig697       120 E  T  K  K  C  F  I  H  S  P  F  L  R  R  A  P  R  S  T  Q  
isotig697       360 GAGACTAAGAAATGCTTCATTCATTCTCCATTCCTAAGACGTGCTCCAAGGTCAACTCAG

isotig697       140 H  H  S  P  R  T  W  L  H  L  V  M  D  R  T  E  S  H  Y  V  
isotig697       420 CACCATTCCCCAAGGACTTGGCTTCATTTAGTTATGGATAGAACTGAAAGTCATTATGTT

isotig697       160 Q  N  K  P  D  L  K  R  L  C  N  F  L  N  M  Q  N  L  K  R  
isotig697       480 CAGAATAAGCCAGACTTGAAGAGGTTGTGTAATTTCTTGAATATGCAAAATCTTAAAAGG

isotig697       180 G  A  C   183
isotig697       540 GGGGCATGC 549
""",
        )
        nucleotide_record = next(nucleotide_records)
        protein_record = protein_alignment.sequences[1]
        self.assertEqual(nucleotide_record.id, protein_record.id)
        alignments = aligner.align(protein_record, nucleotide_record)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 65], [0, 195]]))
        )
        self.assertEqual(
            str(alignment),
            """\
ENSG00000         0 M  K  V  S  A  A  L  L  C  L  L  L  I  A  A  T  F  I  P  Q  
ENSG00000         0 ATGAAAGTCTCTGCCGCCCTTCTGTGCCTGCTGCTCATAGCAGCCACCTTCATTCCCCAA

ENSG00000        20 G  L  A  Q  P  D  A  I  N  A  P  V  T  C  C  Y  N  F  T  N  
ENSG00000        60 GGGCTCGCTCAGCCAGATGCAATCAATGCCCCAGTCACCTGCTGTTATAACTTCACCAAT

ENSG00000        40 R  K  I  S  V  Q  R  L  A  S  Y  R  R  I  T  S  S  K  C  P  
ENSG00000       120 AGGAAGATCTCAGTGCAGAGGCTCGCGAGCTATAGAAGAATCACCAGCAGCAAGTGTCCC

ENSG00000        60 K  E  A  V  M    65
ENSG00000       180 AAAGAAGCTGTGATG 195
""",
        )
        nucleotide_record = next(nucleotide_records)
        protein_record = protein_alignment.sequences[2]
        self.assertEqual(nucleotide_record.id, protein_record.id)
        alignments = aligner.align(protein_record, nucleotide_record)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 99], [9, 306]]))
        )
        self.assertEqual(
            str(alignment),
            """\
ENSG00000         0 M  K  V  S  A  A  L  L  C  L  L  L  I  A  A  T  F  I  P  Q  
ENSG00000         9 ATGAAAGTCTCTGCCGCCCTTCTGTGCCTGCTGCTCATAGCAGCCACCTTCATTCCCCAA

ENSG00000        20 G  L  A  Q  P  D  A  I  N  A  P  V  T  C  C  Y  N  F  T  N  
ENSG00000        69 GGGCTCGCTCAGCCAGATGCAATCAATGCCCCAGTCACCTGCTGTTATAACTTCACCAAT

ENSG00000        40 R  K  I  S  V  Q  R  L  A  S  Y  R  R  I  T  S  S  K  C  P  
ENSG00000       129 AGGAAGATCTCAGTGCAGAGGCTCGCGAGCTATAGAAGAATCACCAGCAGCAAGTGTCCC

ENSG00000        60 K  E  A  V  I  F  K  T  I  V  A  K  E  I  C  A  D  P  K  Q  
ENSG00000       189 AAAGAAGCTGTGATCTTCAAGACCATTGTGGCCAAGGAGATCTGTGCTGACCCCAAGCAG

ENSG00000        80 K  W  V  Q  D  S  M  D  H  L  D  K  Q  T  Q  T  P  K  T  
ENSG00000       249 AAGTGGGTTCAGGATTCCATGGACCACCTGGACAAGCAAACCCAAACTCCGAAGACT

ENSG00000        99
ENSG00000       306
""",
        )
        alignment = protein_alignment.mapall(codon_alignments)
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[0, 42, 126, 126, 231, 333, 549],
                          [0,  0,  84,  90, 195, 195, 195],
                          [9,  9,  93,  99, 204, 306, 306]])
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "clustal"),
            """\
isotig69710                         AGAGGCGATCAACGCAGCAACTTCCAGCTGTCTCCCTCCACCATGCAGAT
ENSG00000108691:ENST0000058090      ------------------------------------------ATGAAAGT
ENSG00000108691:ENST0000022583      ------------------------------------------ATGAAAGT

isotig69710                         CTCCACAGGGCTTCTGTGCcTGCTGCTTGTGGCCACTGGCTTCACTTCCC
ENSG00000108691:ENST0000058090      CTCTGCCGCCCTTCTGTGCCTGCTGCTCATAGCAGCCACCTTCATTCCCC
ENSG00000108691:ENST0000022583      CTCTGCCGCCCTTCTGTGCCTGCTGCTCATAGCAGCCACCTTCATTCCCC

isotig69710                         AGGTGCTGGCTCACCCAGGCTCTATC------CCATCTACCTaCTGCTTT
ENSG00000108691:ENST0000058090      AAGGGCTCGCTCAGCCAGATGCAATCAATGCCCCAGTCACCTGCTGTTAT
ENSG00000108691:ENST0000022583      AAGGGCTCGCTCAGCCAGATGCAATCAATGCCCCAGTCACCTGCTGTTAT

isotig69710                         GTTATGACCAGTAAGAaGATCCCCAAATCACTACTGAaGAGCTACAAAaG
ENSG00000108691:ENST0000058090      AACTTCACCAATAGGAAGATCTCAGTGCAGAGGCTCGCGAGCTATAGAAG
ENSG00000108691:ENST0000022583      AACTTCACCAATAGGAAGATCTCAGTGCAGAGGCTCGCGAGCTATAGAAG

isotig69710                         AATCTCCAACAGCAGATGCACCcTGAAAGCCATACTCTTCAAGACCAAGT
ENSG00000108691:ENST0000058090      AATCACCAGCAGCAAGTGTCCCAAAGAAGCTGTGATG-------------
ENSG00000108691:ENST0000022583      AATCACCAGCAGCAAGTGTCCCAAAGAAGCTGTGATCTTCAAGACCATTG

isotig69710                         CGGGCAAAGAGATCTGTGCTGACCCCAAGAAGAAGTGGGTCcAGGATGCC
ENSG00000108691:ENST0000058090      --------------------------------------------------
ENSG00000108691:ENST0000022583      TGGCCAAGGAGATCTGTGCTGACCCCAAGCAGAAGTGGGTTCAGGATTCC

isotig69710                         ACAAAGCACCTGGACCAAATCCTTCAAACTCCAAAACCGACAATCCCCTC
ENSG00000108691:ENST0000058090      --------------------------------------------------
ENSG00000108691:ENST0000022583      ATGGACCACCTGGACAAGCAAACCCAAACTCCGAAGACT-----------

isotig69710                         TTTTGAGACTCACCCAGAGACTAAGAAATGCTTCATTCATTCTCCATTCC
ENSG00000108691:ENST0000058090      --------------------------------------------------
ENSG00000108691:ENST0000022583      --------------------------------------------------

isotig69710                         TAAGACGTGCTCCAAGGTCAACTCAGCACCATTCCCCAAGGACTTGGCTT
ENSG00000108691:ENST0000058090      --------------------------------------------------
ENSG00000108691:ENST0000022583      --------------------------------------------------

isotig69710                         CATTTAGTTATGGATAGAACTGAAAGTCATTATGTTCAGAATAAGCCAGA
ENSG00000108691:ENST0000058090      --------------------------------------------------
ENSG00000108691:ENST0000022583      --------------------------------------------------

isotig69710                         CTTGAAGAGGTTGTGTAATTTCTTGAATATGCAAAATCTTAAAAGGGGGG
ENSG00000108691:ENST0000058090      --------------------------------------------------
ENSG00000108691:ENST0000022583      --------------------------------------------------

isotig69710                         CATGC
ENSG00000108691:ENST0000058090      -----
ENSG00000108691:ENST0000022583      -----


""",
        )

    def test2(self):
        aligner = CodonAligner()
        nucleotide_records = SeqIO.parse("codonalign/nucl2.fa", "fasta")
        protein_alignment = Align.read("codonalign/pro2.aln", "clustal")
        self.assertEqual(len(protein_alignment.sequences), 3)
        codon_alignments = []
        nucleotide_record = next(nucleotide_records)
        protein_record = protein_alignment.sequences[0]
        self.assertEqual(nucleotide_record.id, protein_record.id)
        alignments = aligner.align(protein_record, nucleotide_record)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 1094], [0, 3282]]))
        )
        self.assertEqual(
            str(alignment),
            """\
isotig351         0 E  R  Q  G  R  W  C  V  P  G  E  A  V  E  A  R  V  S  R  S  
isotig351         0 GAAAGGCAGGGTCGctGGTGCGTGCCCGGCGAGGCTGTGGAGGCCcgTGTGTCTAGAAGC

isotig351        20 C  V  R  E  M  A  E  P  G  R  R  R  G  P  R  S  R  G  G  G  
isotig351        60 TGTGTGAGAGAGATGGCGGAACCCGGGAGGAGACGGGGTCCGAGGTCCCGCGGTGGCGGC

isotig351        40 A  G  R  G  A  R  R  A  R  V  A  R  G  R  R  P  R  A  P  Q  
isotig351       120 GCCGGCCGAGGCgCTCGAAGAGCCCGGGTCGCCCGTGGCCGGCGTCCTCGCGCCCCGCAG

isotig351        60 S  L  S  R  L  I  P  D  T  V  L  V  D  L  V  S  D  S  D  E  
isotig351       180 TCTCTGTCCCGGCTCATTCCAGACACGGTGCTTGTGGACTTGGTCAGTGACAGCGACGAG

isotig351        80 E  I  L  E  V  V  A  D  P  V  E  A  P  A  A  R  A  P  A  P  
isotig351       240 GAGATCCTGGAAGTCGTCGCGGACCCGGTAGAAGCGCCCGCCGCCCGGGCcCCcGCGCCG

isotig351       100 A  A  H  G  Q  D  S  D  S  D  S  A  G  A  D  E  G  P  A  G  
isotig351       300 GCCGCACATGGGCAGGACAGCGACAGCGACAGTGCAGGGGCGgACGAGGGGCCTGCAGGA

isotig351       120 A  P  Q  T  L  V  R  R  R  R  R  R  L  L  D  P  G  E  A  P  
isotig351       360 GCCCCTCAGACCTTGGTCCGGCGGCGGCGCCGGCGGCTGCTGGATCCCGGCGAGGCACCG

isotig351       140 V  V  P  V  Y  S  G  K  V  Q  S  S  L  N  L  I  P  D  N  S  
isotig351       420 GTGGTTCCTGTGTACTCCGGGAAGGTACAAAGCAGCCTCAACCTCATCCCAGATAATTCA

isotig351       160 S  L  L  K  L  C  P  S  E  P  E  D  E  A  D  V  T  D  C  G  
isotig351       480 TCCCTCTTGAAACTTTGCCCCTCAGAGCCTGAAGATGAGGCAGATGTGACAGATTGTGGC

isotig351       180 S  P  P  P  E  D  A  L  I  P  G  S  P  W  K  K  K  L  R  N  
isotig351       540 AGTCCTCCTCCTGAGGATGCCCTAATTCCAGGTTCTCCCTGGAAGAAGAAGCTGAGGAAT

isotig351       200 K  H  E  K  E  E  M  K  M  E  E  F  P  D  Q  D  I  S  P  L  
isotig351       600 AAGCATGAAAAaGAAGAGATGAAGATGGAAGAGTTtCCGGACCAGGACATCTCTCCTTTG

isotig351       220 P  R  P  S  S  R  N  K  S  R  K  H  T  E  A  L  Q  K  L  R  
isotig351       660 CCCCGACCTTCATCAAGAAaCAAAAGCAGAAAGCATACCGAGGCACTCCAGAAGTTGAGG

isotig351       240 E  V  N  K  R  L  Q  D  L  R  S  C  L  S  P  K  Q  H  Q  S  
isotig351       720 GAAGTGAACAAGCGCCTCCAAGATCTCCGATCCTGCCTGAGCCCCAAGCAGCACCAGAGT

isotig351       260 P  A  L  Q  N  P  D  D  E  V  V  L  V  D  G  P  V  L  S  Q  
isotig351       780 CCAGCCCTTCAGAACCCAGATGATGAGGTGGTCCTCGTGGACGGGCCTGtCTTGTCACAG

isotig351       280 S  P  R  L  F  T  L  K  I  R  C  R  A  D  L  V  R  L  P  V  
isotig351       840 AGCCCGAGACTCTTCACCCTCAAGATCCGGTGCCGGGCTGACCTAGTCAGATTGCCCGTC

isotig351       300 M  T  S  E  P  L  Q  N  V  V  D  Y  M  A  N  H  L  G  V  S  
isotig351       900 ATGACATCGGAACCCCTTCAGAATGTGGTGGATTACATGGCCAATCATCTTGGGGTGTCT

isotig351       320 P  S  R  I  L  L  L  F  G  E  T  E  L  S  P  T  A  T  P  R  
isotig351       960 CCAAGCAGGATTCTTTTACTCTTTGGAGAGACAGAACTGTCCCCTACTGCCACCCCTAGG

isotig351       340 T  L  K  L  G  V  A  D  I  I  D  C  V  V  L  T  S  S  S  E  
isotig351      1020 ACCCTAAAGCTTGGTGTGGCTGACATCATTGATTGTGTGGTGCTGACAAGTTCTTCAGAG

isotig351       360 A  T  E  T  T  Q  Q  L  C  L  R  V  Q  G  K  E  K  H  Q  M  
isotig351      1080 GCCACAGAGACAACCCAGCAGCTCTGCCTCCGGGTGCAGGGGAAGGAGAAGCACCAGATG

isotig351       380 L  E  I  S  L  S  P  D  S  P  L  E  V  L  M  A  H  Y  E  E  
isotig351      1140 TTGGAGATCTCACTCTCTCCTGACTCtCCTCTTGAGGTCCTCATGGCGCACTATGAGGAG

isotig351       400 A  M  G  L  S  G  H  K  L  S  F  F  F  D  G  T  K  L  S  G  
isotig351      1200 GCCATGGGACTCTCCGGACACAAGCTCTCCTTCTTCTTCGATGGGACAAAGCTGTCGGGC

isotig351       420 K  E  L  P  A  D  L  G  M  E  S  G  D  L  I  E  V  W  G  S  
isotig351      1260 AAGGAGCTGCCAGCTGATCTGGGCATGGAATCCGGGGATCTCATTGAAGTTTGGGGCAGC

isotig351       440 F  L  L  L  F  G  C  R  A  K  T  W  G  Q  Q  L  P  L  L  L  
isotig351      1320 TTCCTCCTCCTGTTTGGATGCAGAGCCAAGACTTGGGGACAACAGCTCCCACTTTtATTA

isotig351       460 L  F  F  A  P  G  L  T  E  T  E  L  E  L  V  Y  L  F  P  A  
isotig351      1380 TTATTTTTtGCCCcAGGGCTAACAGAAACCGAATTAGAACTCGTTTATTTATTTCCGGCA

isotig351       480 L  G  I  E  P  R  A  V  H  M  L  R  M  C  E  V  E  V  K  P  
isotig351      1440 CTGGGGATTGAACCCAGGGCTGTGCATATGCTAAGGATGTGTGAAGTTGAGGTAAAaCCA

isotig351       500 R  H  D  L  C  P  V  S  L  T  V  V  S  S  G  L  I  G  N  C  
isotig351      1500 AGGCATGACCTTTGCCcTGTCTCGTTGACCGTAGTCTCAAGCGGTCTGATTGGTAATTGT

isotig351       520 V  T  V  A  A  L  G  G  L  W  L  C  V  L  V  L  C  T  A  A  
isotig351      1560 GTGACTGTGGCTGCCCTGGGTGgCCTGTGGCTGTGTGTGTTGGTGCTGTGTACAGCAGCT

isotig351       540 P  G  A  W  R  R  G  I  G  F  P  T  I  S  F  S  S  A  Q  M  
isotig351      1620 CCTGGGGCATGGAGAAGGGGTATTGGCTTCCCTACCATTTCGTTCAGTAGTGCACAAATG

isotig351       560 S  L  A  F  G  C  P  G  L  F  M  P  H  W  G  P  T  G  I  L  
isotig351      1680 AGCCTTGCATTTGGGTGCCCAGGCTTGTTTATGCCACATTGGGGACCAACAGGGATTTTA

isotig351       580 I  L  I  W  G  L  R  C  R  Y  L  P  W  W  G  S  V  A  G  K  
isotig351      1740 ATTCTCATTTGGGGGCTGAGATGCAGGTACCTTCCCTGGTGGGGATCGGTTGCAGGGAAA

isotig351       600 T  S  N  D  I  Y  A  E  A  E  D  V  D  S  G  F  L  K  A  N  
isotig351      1800 ACAAGCAATGACATCTATGCTGAAGCTGAGGACGTAGACAGTGGGTTTTTAAAGGCTAAC

isotig351       620 R  S  G  V  F  S  P  G  H  P  L  S  Q  D  S  A  S  I  W  Y  
isotig351      1860 AGGAGTGGTGTCTTCAGCCCTGGGCATCCACTCTCCCAGGACTCGGCCAGCATCTGGTAC

isotig351       640 I  L  P  G  M  H  I  P  H  L  S  T  C  S  C  V  C  L  L  F  
isotig351      1920 ATACTTCCTGGCATGCACATCCCACATCTGAGCACATGCAGCTGTGTTTGTTTACTCTTC

isotig351       660 R  P  P  A  L  S  C  Q  H  P  L  G  I  F  V  S  E  C  L  Y  
isotig351      1980 cGTCCTCCCGCCCTCTCCTGTCAGCACCCACTTGGTATATTTGTATCTGAATGTCTTTAT

isotig351       680 K  M  I  F  M  C  V  I  L  V  Y  V  H  G  T  V  Q  T  F  C  
isotig351      2040 AAGATGATTTTCATGTGTGTGATTTTAGTGTATGTACATGGTACTGTGCAAACATTCTGT

isotig351       700 L  H  S  V  V  T  G  P  C  H  L  L  D  L  L  L  F  R  S  A  
isotig351      2100 CTTCACTCAGTGGTGACAGGCCCATGTCACCTTCTAGATCTGTTGTTGTTCCGATCTGCT

isotig351       720 V  A  A  P  W  H  M  I  L  L  T  H  C  L  N  W  R  D  T  V  
isotig351      2160 GTGGCTGCCCCATGGCATATGATATTACTCACACATTGCCTTAATTGGAGGGATACCGTC

isotig351       740 A  S  S  E  P  N  S  R  P  Q  Q  T  S  F  C  V  F  P  S  G  
isotig351      2220 GCCAGCTCCGAACCAAATTCCAGACCACAGCAAACATCTTTTTGTGTGTTCCCTTCTGGG

isotig351       760 P  G  V  I  V  W  H  T  C  P  G  W  D  C  G  L  R  R  F  I  
isotig351      2280 CCAGGGGTAATTGTCTGGCATACTTGCCCAGGATGGGACTgTGGGTTACGCAGGTTCATT

isotig351       780 S  L  V  F  L  G  F  S  S  E  C  G  L  L  F  L  E  L  L  F  
isotig351      2340 TcTCTGGTTTTTCTAGGgTTTTCTTCAGAATGTGGACTGCTTTTCTTGGAGCTTCTcTTT

isotig351       800 P  S  L  D  Q  Y  L  V  L  L  Y  F  L  I  C  Q  S  T  V  Y  
isotig351      2400 CCTTcTCTTGACCAGTATTTAGTGTTGttgTACTTTCTGATTTGCCAGTCTACTGTGTAT

isotig351       820 K  A  L  H  I  L  I  S  F  I  P  P  F  L  I  N  V  H  L  L  
isotig351      2460 AAAGCACTGCATATCTTAATATCCTTTATCCCACCATTTCTTATCAATGTCCACCTTCTG

isotig351       840 V  F  L  C  L  L  V  F  L  G  F  L  L  C  V  F  S  V  R  S  
isotig351      2520 GTTTTTCTTTGTCTGCTTGTTTTCTTGGGCTTCCTTTTGTGTGTTTTTtCTGTGAGGAGT

isotig351       860 C  L  C  S  R  Q  L  Q  I  T  F  R  I  A  L  C  V  P  R  G  
isotig351      2580 TGCTTGTGTTCTaGgCAATTGCAAATAACCTTCAGAATAGCACTTTGTGTACCCAGGGGG

isotig351       880 I  F  T  K  Q  T  F  L  T  L  M  Q  T  A  G  F  C  K  G  F  
isotig351      2640 ATTTTCACCAAACAGACATTCTTAACTTTAATGCAGACAGCTGGCTTTTGTAAAGGGTTT

isotig351       900 F  V  L  C  F  F  F  F  C  F  C  F  C  F  L  R  E  N  E  K  
isotig351      2700 TTTGTTTTGTGTTTtTTTTTtttttGCTTTTGTTTTtGcTTTTTAAGGGAAAATGAGAAA

isotig351       920 E  L  P  E  W  Q  G  F  R  G  V  S  G  L  L  L  L  D  L  T  
isotig351      2760 GAACTTCCAGAATGGCAGGGGTTCCGGGGAGTCTCAGGATTGCTACTTTTGGATCTCACT

isotig351       940 S  E  T  F  I  V  F  C  L  A  S  F  L  F  F  C  F  Y  L  S  
isotig351      2820 AGTGAAACTTTCATTGTCTTCTGTCTTGCTAGTTTTTTGTTCTTTTGTTTTTATCTCTCA

isotig351       960 Y  D  H  M  R  I  E  F  L  I  C  C  G  V  D  T  D  C  S  E  
isotig351      2880 TATGACCACATGAGAATAGAGTTCCTGATCTGCTGTGGTGTGGATaCTGACTGTTCTGAa

isotig351       980 L  S  T  L  L  L  T  G  Y  S  T  V  L  V  P  D  Y  C  M  L  
isotig351      2940 CTTTCCACTCTGCTTCTCACTGGTTACAGCACTGTACTGGTTCCTGACTATTGTATGTtA

isotig351      1000 R  L  N  H  L  K  C  T  F  V  W  E  S  K  F  Y  C  V  L  I  
isotig351      3000 CGTCTAAACCATCTGAAATGCACCTTTGTGTGGGAGAGCAAATTTTATTGTGTTTTAATT

isotig351      1020 I  M  S  Q  F  S  K  L  C  L  L  K  N  Y  V  V  L  S  Q  E  
isotig351      3060 ATAATGAGCCAGTTTTCtAAGCTTTGTCTaCTAAAAAATTATGTAGTCCTTTCACAAGAA

isotig351      1040 V  C  Q  L  L  K  N  F  F  C  E  N  M  L  Y  L  I  L  I  F  
isotig351      3120 GTATGTCAACTTCTtAAAAATTTTTTTtGtGAGAATATGTTATACCTTATTTTAATATTT

isotig351      1060 K  N  S  M  T  Y  S  F  Q  V  S  Q  L  L  Y  D  N  I  T  H  
isotig351      3180 AAAAATTCCATGACATATAGCTTCCAGGTTTCTCAGCTCTTGTATGACAATATCACACAT

isotig351      1080 Y  Y  Y  V  F  N  L  S  Q  R  E  M  P  L   1094
isotig351      3240 TaCTATTATGTATTCAATCTCAGCCAAAGGGAGATGCCTTTA 3282
""",
        )
        nucleotide_record = next(nucleotide_records)
        protein_record = protein_alignment.sequences[1]
        self.assertEqual(nucleotide_record.id, protein_record.id)
        alignments = aligner.align(protein_record, nucleotide_record)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 1104], [0, 3312]]))
        )
        self.assertEqual(
            str(alignment),
            """\
isotig351         0 E  R  Q  G  R  W  C  V  P  G  E  A  V  E  A  R  V  S  R  S  
isotig351         0 GAAAGGCAGGGTCGctGGTGCGTGCCCGGCGAGGCTGTGGAGGCCcgTGTGTCTAGAAGC

isotig351        20 C  V  R  E  M  A  E  P  G  R  R  R  G  P  R  S  R  G  G  G  
isotig351        60 TGTGTGAGAGAGATGGCGGAACCCGGGAGGAGACGGGGTCCGAGGTCCCGCGGTGGCGGC

isotig351        40 A  G  R  G  A  R  R  A  R  V  A  R  G  R  R  P  R  A  P  Q  
isotig351       120 GCCGGCCGAGGCgCTCGAAGAGCCCGGGTCGCCCGTGGCCGGCGTCCTCGCGCCCCGCAG

isotig351        60 S  L  S  R  L  I  P  D  T  V  L  V  D  L  V  S  D  S  D  E  
isotig351       180 TCTCTGTCCCGGCTCATTCCAGACACGGTGCTTGTGGACTTGGTCAGTGACAGCGACGAG

isotig351        80 E  I  L  E  V  V  A  D  P  V  E  A  P  A  A  R  A  P  A  P  
isotig351       240 GAGATCCTGGAAGTCGTCGCGGACCCGGTAGAAGCGCCCGCCGCCCGGGCcCCcGCGCCG

isotig351       100 A  A  H  G  Q  D  S  D  S  D  S  A  G  A  D  E  G  P  A  G  
isotig351       300 GCCGCACATGGGCAGGACAGCGACAGCGACAGTGCAGGGGCGgACGAGGGGCCTGCAGGA

isotig351       120 A  P  Q  T  L  V  R  R  R  R  R  R  L  L  D  P  G  E  A  P  
isotig351       360 GCCCCTCAGACCTTGGTCCGGCGGCGGCGCCGGCGGCTGCTGGATCCCGGCGAGGCACCG

isotig351       140 V  V  P  V  Y  S  G  K  V  Q  S  S  L  N  L  I  P  D  N  S  
isotig351       420 GTGGTTCCTGTGTACTCCGGGAAGGTACAAAGCAGCCTCAACCTCATCCCAGATAATTCA

isotig351       160 S  L  L  K  L  C  P  S  E  P  E  D  E  A  D  V  T  D  C  G  
isotig351       480 TCCCTCTTGAAACTTTGCCCCTCAGAGCCTGAAGATGAGGCAGATGTGACAGATTGTGGC

isotig351       180 S  P  P  P  E  D  A  L  I  P  G  S  P  W  K  K  K  L  R  N  
isotig351       540 AGTCCTCCTCCTGAGGATGCCCTAATTCCAGGTTCTCCCTGGAAGAAGAAGCTGAGGAAT

isotig351       200 K  H  E  K  E  E  M  K  M  E  E  F  P  D  Q  D  I  S  P  L  
isotig351       600 AAGCATGAAAAaGAAGAGATGAAGATGGAAGAGTTtCCGGACCAGGACATCTCTCCTTTG

isotig351       220 P  R  P  S  S  R  N  K  S  R  K  H  T  E  A  L  Q  K  L  R  
isotig351       660 CCCCGACCTTCATCAAGAAaCAAAAGCAGAAAGCATACCGAGGCACTCCAGAAGTTGAGG

isotig351       240 E  V  N  K  R  L  Q  D  L  R  S  C  L  S  P  K  Q  H  Q  S  
isotig351       720 GAAGTGAACAAGCGCCTCCAAGATCTCCGATCCTGCCTGAGCCCCAAGCAGCACCAGAGT

isotig351       260 P  A  L  Q  N  P  D  D  E  V  V  L  V  D  G  P  V  L  S  Q  
isotig351       780 CCAGCCCTTCAGAACCCAGATGATGAGGTGGTCCTCGTGGACGGGCCTGtCTTGTCACAG

isotig351       280 S  P  R  L  F  T  L  K  I  R  C  R  A  D  L  V  R  L  P  V  
isotig351       840 AGCCCGAGACTCTTCACCCTCAAGATCCGGTGCCGGGCTGACCTAGTCAGATTGCCCGTC

isotig351       300 M  T  S  E  P  L  Q  N  V  V  D  Y  M  A  N  H  L  G  V  S  
isotig351       900 ATGACATCGGAACCCCTTCAGAATGTGGTGGATTACATGGCCAATCATCTTGGGGTGTCT

isotig351       320 P  S  R  I  L  L  L  F  G  E  T  E  L  S  P  T  A  T  P  R  
isotig351       960 CCAAGCAGGATTCTTTTACTCTTTGGAGAGACAGAACTGTCCCCTACTGCCACCCCTAGG

isotig351       340 T  L  K  L  G  V  A  D  I  I  D  C  V  V  L  T  S  S  S  E  
isotig351      1020 ACCCTAAAGCTTGGTGTGGCTGACATCATTGATTGTGTGGTGCTGACAAGTTCTTCAGAG

isotig351       360 A  T  E  T  T  Q  Q  L  C  L  R  V  Q  G  K  E  K  H  Q  M  
isotig351      1080 GCCACAGAGACAACCCAGCAGCTCTGCCTCCGGGTGCAGGGGAAGGAGAAGCACCAGATG

isotig351       380 L  E  I  S  L  S  P  D  S  P  L  E  V  L  M  A  H  Y  E  E  
isotig351      1140 TTGGAGATCTCACTCTCTCCTGACTCgCCTCTTGAGGTCCTCATGGCGCACTATGAGGAG

isotig351       400 A  M  G  L  S  G  H  K  L  S  F  F  F  D  G  T  K  L  S  G  
isotig351      1200 GCCATGGGACTCTCCGGACACAAGCTCTCCTTCTTCTTCGATGGGACAAAGCTGTCGGGC

isotig351       420 K  E  L  P  A  D  L  G  M  E  S  G  D  L  I  E  V  W  G  S  
isotig351      1260 AAGGAACTGCCAGCTGATCTGGGCATGGAATCCGGGGATCTCATTGAAGTTTGGGGCAGC

isotig351       440 F  L  L  L  F  G  C  R  A  K  T  W  G  Q  Q  L  P  L  L  L  
isotig351      1320 TTCCTCCTCCTGTTTGGATGCAGAGCCAAGACTTGGGGACAACAGCTCCCACTTTTATTA

isotig351       460 L  F  F  A  P  G  L  T  E  T  E  L  E  L  V  Y  L  F  P  A  
isotig351      1380 TTATTTTTTGCCCCAGGGCTAACAGAAACCGAATTAGAACTCGTTTATTTATTTCCGGCA

isotig351       480 L  G  I  E  P  W  A  V  H  M  L  R  M  C  E  V  E  V  K  P  
isotig351      1440 CTGGGGATTGAACCCTGGGCTGTGCATATGCTAAGGATGTGTGAAGTTGAGGTAAAACCA

isotig351       500 R  H  D  L  C  P  V  S  L  T  D  W  L  C  D  C  G  C  P  G  
isotig351      1500 AGGCATGACCTTTGCCCTGTCTCGTTGACCGATTGGTTGTGTGACTGTGGCTGCCCTGGG

isotig351       520 W  P  V  A  V  C  V  G  A  V  D  S  S  S  W  G  M  E  K  G  
isotig351      1560 TGGCCTGTGGCTGTGTGTGTTGGTGCTGTGGACAGCAGCTCCTGGGGCATGGAGAAGGGG

isotig351       540 Y  W  L  P  Y  H  F  V  Q  C  T  N  E  P  C  I  W  V  P  R  
isotig351      1620 TATTGGCTTCCCTACCATTTCGTTCAGTGCACAAATGAGCCTTGCATTTGGGTGCCCAGG

isotig351       560 L  V  I  Y  A  T  L  G  T  N  R  D  F  N  S  H  L  G  A  E  
isotig351      1680 CTTGTGATTTATGCCACATTGGGGACCAACAGGGATTTTAATTCTCATTTGGGGGCTGAG

isotig351       580 M  Q  V  P  S  L  V  G  I  V  G  C  R  E  I  N  K  Q  H  L  
isotig351      1740 ATGCAGGTACCTTCCCTGGTGGGGATAGTCGGTTGCAGGGAAATAAACAAGCAACATCTA

isotig351       600 C  S  G  L  G  R  Q  W  V  I  F  K  G  Q  E  W  C  L  Q  P  
isotig351      1800 TGCAGCGGACTAGGTAGACAGTGGGTAATTTTTAAAGGCCAGGAGTGGTGTCTTCAGCCC

isotig351       620 W  A  S  T  L  P  G  L  G  Q  H  L  V  H  T  S  W  H  A  H  
isotig351      1860 TGGGCATCCACTCTCCCAGGACTCGGCCAGCATCTGGTACATACTTCCTGGCATGCACAT

isotig351       640 L  A  T  S  E  H  M  Q  L  L  S  L  M  F  T  L  V  P  S  S  
isotig351      1920 CTAGCCACATCTGAGCACATGCAGCTGTTGAGTTTGATGTTTACTCTAGTTCcGTCCTCC

isotig351       660 R  P  L  L  S  A  P  T  W  Y  I  C  I  M  S  L  D  D  F  L  
isotig351      1980 CGCCCTCTCCTGTCAGCACCCACTTGGTATATTTGTATCATGTCTTTAGATGATTTTCTG

isotig351       680 N  V  C  D  F  S  V  C  T  W  Y  C  A  I  D  I  L  S  S  L  
isotig351      2040 AATGTGTGTGATTTTAGTGTATGTACATGGTACTGTGCAATAGACATTCTGTCTTCACTC

isotig351       700 S  V  S  D  R  L  A  M  L  A  P  S  R  S  V  V  V  P  I  C  
isotig351      2100 AGTGTGAGTGACAGGCTAGCCATGTTAGCACCTTCTAGATCTGTTGTTGTTCCGATCTGC

isotig351       720 C  G  C  P  I  V  A  L  D  D  I  T  H  T  L  P  L  E  G  Y  
isotig351      2160 TGTGGCTGCCCCATAGTGGCATTAGATGATATTACTCACACATTGCCTTTGGAGGGATAC

isotig351       740 L  S  R  Q  L  L  G  T  K  F  Q  T  T  A  N  I  F  L  C  V  
isotig351      2220 CTGAGTCGCCAGCTCCTAGGAACCAAATTCCAGACCACAGCAAACATCTTTTTGTGTGTT

isotig351       760 P  F  W  A  R  V  G  N  C  L  A  Y  L  P  R  M  G  L  W  V  
isotig351      2280 CCCTTCTGGGCCAGGGTAGGTAATTGTCTGGCATACTTGCCCAGGATGGGACTgTGGGTT

isotig351       780 T  Q  V  H  F  S  G  F  S  R  V  F  F  R  M  W  T  A  F  L  
isotig351      2340 ACGCAGGTTCATTTcTCTGGTTTTTCTAGGgTTTTCTTCAGAATGTGGACTGCTTTTCTT

isotig351       800 G  A  S  L  S  F  S  P  V  F  S  V  V  V  L  S  D  L  P  V  
isotig351      2400 GGAGCTTCTcTTTCCTTcTCTCCAGTATTTAGTGTTGttgTACTTTCTGATTTGCCAGTC

isotig351       820 Y  C  V  S  T  A  Y  L  N  I  L  Y  L  T  T  I  S  Y  Q  C  
isotig351      2460 TACTGTGTAAGCACTGCATATCTTAATATCCTTTATCTGACCACCATTTCTTATCAATGT

isotig351       840 P  P  S  G  F  S  L  L  A  A  C  F  L  V  R  L  P  F  V  C  
isotig351      2520 CCACCTTCTGGTTTTTCTTTGTTAGCTGCTTGTTTTCTTGTGAGGCTTCCTTTTGTGTGT

isotig351       860 F  F  C  E  E  L  L  V  F  A  I  A  N  N  L  Q  N  S  T  L  
isotig351      2580 TTTTtCTGTGAGGAGTTGCTTGTGTTCgCAATTGCAAATAACCTTCAGAATAGCACTTTG

isotig351       880 C  T  Q  G  D  F  H  Q  T  D  I  L  N  F  N  A  D  S  W  L  
isotig351      2640 TGTACCCAGGGGGATTTTCACCAAACAGACATTCTTAACTTTAATGCAGACAGCTGGCTT

isotig351       900 L  L  E  R  V  F  C  F  V  F  F  F  F  L  L  L  F  L  L  F  
isotig351      2700 TTGTTAGAAAGGGTTTTTTGTTTTGTGTTTtTTTTTtttttGCTTTTGTTTTtGcTTTTT

isotig351       920 K  G  K  E  R  T  S  R  M  A  G  V  L  A  G  S  L  R  I  A  
isotig351      2760 AAGGGAAAAGAAAGAACTTCCAGAATGGCAGGGGTTCTAGCGGGGAGTCTCAGGATTGCT

isotig351       940 T  F  G  S  L  N  L  R  N  F  H  C  L  L  L  A  C  F  F  V  
isotig351      2820 ACTTTTGGATCTCTAAACTTAAGAAACTTTCATTGTCTTCTGTTAGCTTGCTTTTTTGTT

isotig351       960 L  L  F  L  S  L  I  L  R  P  L  N  E  N  R  V  L  T  D  L  
isotig351      2880 CTTTTGTTTTTATCTCTCATATTGAGACCACTAAATGAGAATAGAGTTCTAACTGATCTG

isotig351       980 L  W  C  G  Y  L  F  T  F  H  S  A  S  H  W  L  Q  L  D  C  
isotig351      2940 CTGTGGTGTGGATaCCTGTTCaCTTTCCACTCTGCTTCTCACTGGTTACAGCTAGACTGT

isotig351      1000 T  G  S  L  L  Y  V  T  S  K  P  S  E  M  H  L  C  V  G  E  
isotig351      3000 ACTGGTTCCCTATTGTATGTtACGTCTAAACCATCTGAAATGCACCTTTGTGTGGGAGAG

isotig351      1020 Q  I  L  L  C  F  N  Y  N  E  P  V  F  A  L  S  T  K  K  L  
isotig351      3060 CAAATTTTATTGTGTTTTAATTATAATGAGCCAGTTTTCGCTTTGTCTaCTAAAAAATTA

isotig351      1040 I  C  S  L  T  F  T  R  S  M  S  T  S  K  F  F  L  E  Y  V  
isotig351      3120 ATATGTAGTCTGACTTTCACAAGAAGTATGTCAACTTCTAAATTTTTTTtGGAATATGTT

isotig351      1060 I  P  L  N  F  N  I  K  F  H  D  I  V  L  P  G  F  S  A  L  
isotig351      3180 ATACCTTTGAATTTTAATATTAAATTCCATGACATAGTACTTCCAGGTTTCTCAGCTCTT

isotig351      1080 V  Q  L  M  N  L  N  T  L  L  L  C  I  Q  S  Q  P  I  K  G  
isotig351      3240 GTACAATTAATGAATCTGAACACATTaCTATTATGTATTCAATCTCAGCCAATAAAGGGA

isotig351      1100 D  A  F  N   1104
isotig351      3300 GATGCCTTTAAC 3312
""",
        )
        nucleotide_record = next(nucleotide_records)
        protein_record = protein_alignment.sequences[2]
        self.assertEqual(nucleotide_record.id, protein_record.id)
        alignments = aligner.align(protein_record, nucleotide_record)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 419], [0, 1257]]))
        )
        self.assertEqual(
            str(alignment),
            """\
ENSG00000         0 M  A  E  P  V  G  K  R  G  R  W  S  G  G  S  G  A  G  R  G  
ENSG00000         0 ATGGCGGAGCCTGTGGGGAAGCGGGGCCGCTGGTCCGGAGGTAGCGGTGCCGGCCGAGGG

ENSG00000        20 G  R  G  G  W  G  G  R  G  R  R  P  R  A  Q  R  S  P  S  R  
ENSG00000        60 GGTCGGGGCGGCTGGGGCGGTCGGGGCCGGCGTCCTCGGGCCCAGCGGTCTCCATCCCGG

ENSG00000        40 G  T  L  D  V  V  S  V  D  L  V  T  D  S  D  E  E  I  L  E  
ENSG00000       120 GGCACGCTGGACGTAGTGTCTGTGGACTTGGTCACCGACAGCGATGAGGAAATTCTGGAG

ENSG00000        60 V  A  T  A  R  G  A  A  D  E  V  E  V  E  P  P  E  P  P  G  
ENSG00000       180 GTCGCCACCGCTCGCGGTGCCGCGGACGAGGTTGAGGTGGAGCCCCCGGAGCCCCCGGGG

ENSG00000        80 P  V  A  S  R  D  N  S  N  S  D  S  E  G  E  D  R  R  P  A  
ENSG00000       240 CCGGTCGCGTCCCGGGATAACAGCAACAGTGACAGCGAAGGGGAGGACAGGCGGCCCGCA

ENSG00000       100 G  P  P  R  E  P  V  R  R  R  R  R  L  V  L  D  P  G  E  A  
ENSG00000       300 GGACCCCCGCGGGAGCCGGTCAGGCGGCGGCGGCGGCTGGTGCTGGATCCGGGGGAGGCG

ENSG00000       120 P  L  V  P  V  Y  S  G  K  V  K  S  S  L  R  L  I  P  D  D  
ENSG00000       360 CCGCTGGTTCCGGTGTACTCGGGGAAGGTTAAAAGCAGCCTTCGCCTTATCCCAGATGAT

ENSG00000       140 L  S  L  L  K  L  Y  P  P  G  D  E  E  E  A  E  L  A  D  S  
ENSG00000       420 CTATCCCTCCTGAAACTCTACCCTCCAGGGGATGAGGAAGAGGCAGAGCTGGCAGATTCG

ENSG00000       160 S  G  L  Y  H  E  G  S  P  S  P  G  S  P  W  K  T  K  L  R  
ENSG00000       480 AGTGGTCTCTACCATGAGGGCTCCCCATCACCAGGCTCTCCCTGGAAGACAAAGCTGAGG

ENSG00000       180 T  K  D  K  E  E  K  K  K  T  E  F  L  D  L  D  N  S  P  L  
ENSG00000       540 ACTAAGGATAAAGAAGAGAAGAAAAAGACAGAGTTTCTGGATCTGGACAACTCTCCTCTG

ENSG00000       200 S  P  P  S  P  R  T  K  S  R  T  H  T  R  A  L  K  K  L  S  
ENSG00000       600 TCCCCACCTTCACCAAGGACCAAAAGCAGAACGCATACTCGGGCACTCAAGAAGTTAAGT

ENSG00000       220 E  V  N  K  R  L  Q  D  L  R  S  C  L  S  P  K  P  P  Q  G  
ENSG00000       660 GAGGTGAACAAGCGCCTCCAGGATCTCCGTTCCTGTCTGAGCCCCAAGCCACCTCAGGGT

ENSG00000       240 Q  E  Q  Q  G  Q  E  D  E  V  V  L  V  E  G  P  T  L  P  E  
ENSG00000       720 CAAGAGCAACAGGGCCAAGAGGATGAAGTGGTCTTGGTGGAAGGGCCCACCCTCCCAGAG

ENSG00000       260 T  P  R  L  F  P  L  K  I  R  C  R  A  D  L  V  R  L  P  L  
ENSG00000       780 ACCCCCCGACTCTTCCCACTCAAAATCCGTTGCCGGGCTGACCTGGTCAGATTGCCCCTC

ENSG00000       280 R  M  S  E  P  L  Q  S  V  V  D  H  M  A  T  H  L  G  V  S  
ENSG00000       840 AGGATGTCGGAGCCCCTGCAGAGTGTGGTGGACCACATGGCCACCCACCTTGGGGTGTCC

ENSG00000       300 P  S  R  I  L  L  L  F  G  E  T  E  L  S  P  T  A  T  P  R  
ENSG00000       900 CCAAGCAGGATCCTTTTGCTTTTTGGAGAGACAGAGCTATCACCTACTGCCACTCCCAGG

ENSG00000       320 T  L  K  L  G  V  A  D  I  I  D  C  V  V  L  T  S  S  P  E  
ENSG00000       960 ACCCTAAAGCTCGGAGTGGCTGACATCATTGACTGTGTGGTACTAACAAGTTCTCCAGAG

ENSG00000       340 A  T  E  T  S  Q  Q  L  Q  L  R  V  Q  G  K  E  K  H  Q  T  
ENSG00000      1020 GCCACAGAGACGTCCCAACAGCTCCAGCTCCGGGTGCAGGGAAAGGAGAAACACCAGACA

ENSG00000       360 L  E  V  S  L  S  R  D  S  P  L  K  T  L  M  S  H  Y  E  E  
ENSG00000      1080 CTGGAAGTCTCACTGTCTCGAGATTCCCCTCTAAAGACCCTCATGTCCCACTATGAGGAG

ENSG00000       380 A  M  G  L  S  G  R  K  L  S  F  F  F  D  G  T  K  L  S  G  
ENSG00000      1140 GCCATGGGACTGTCGGGACGGAAGCTCTCCTTCTTCTTTGATGGGACAAAGCTTTCAGGC

ENSG00000       400 R  E  L  P  A  D  L  G  M  E  S  G  D  L  I  E  V  W  G  
ENSG00000      1200 AGGGAGCTGCCAGCTGACCTGGGCATGGAATCTGGGGACCTCATTGAGGTCTGGGGC

ENSG00000       419
ENSG00000      1257
""",
        )
        alignment = protein_alignment.mapall(codon_alignments)
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[   0,   72,  255,  255,  603,  606, 1317, 1530,
                           1530, 1737, 1737, 1902, 1902, 1926, 1926, 2031,
                           2088, 2343, 2346, 2418, 2442, 2556, 2562, 2763,
                           2784, 2853, 2859, 3165, 3165, 3282, 3282],
                          [   0,   72,  255,  255,  603,  606, 1317, 1530,
                           1569, 1776, 1821, 1986, 1992, 2016, 2031, 2136,
                           2136, 2391, 2391, 2463, 2463, 2577, 2577, 2778,
                           2778, 2847, 2847, 3153, 3162, 3279, 3312],
                          [   0,    0,  183,  198,  546,  546, 1257, 1257,
                           1257, 1257, 1257, 1257, 1257, 1257, 1257, 1257,
                           1257, 1257, 1257, 1257, 1257, 1257, 1257, 1257,
                           1257, 1257, 1257, 1257, 1257, 1257, 1257]])
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "clustal"),
            """\
isotig35100                         GAAAGGCAGGGTCGctGGTGCGTGCCCGGCGAGGCTGTGGAGGCCcgTGT
isotig35101                         GAAAGGCAGGGTCGctGGTGCGTGCCCGGCGAGGCTGTGGAGGCCcgTGT
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         GTCTAGAAGCTGTGTGAGAGAGATGGCGGAACCCGGGAGGAGACGGGGTC
isotig35101                         GTCTAGAAGCTGTGTGAGAGAGATGGCGGAACCCGGGAGGAGACGGGGTC
ENSG00000176953:ENST0000032080      ----------------------ATGGCGGAGCCTGTGGGGAAGCGGGGCC

isotig35100                         CGAGGTCCCGCGGTGGCGGCGCCGGCCGAGGCgCTCGAAGAGCCCGGGTC
isotig35101                         CGAGGTCCCGCGGTGGCGGCGCCGGCCGAGGCgCTCGAAGAGCCCGGGTC
ENSG00000176953:ENST0000032080      GCTGGTCCGGAGGTAGCGGTGCCGGCCGAGGGGGTCGGGGCGGCTGGGGC

isotig35100                         GCCCGTGGCCGGCGTCCTCGCGCCCCGCAGTCTCTGTCCCGGCTCATTCC
isotig35101                         GCCCGTGGCCGGCGTCCTCGCGCCCCGCAGTCTCTGTCCCGGCTCATTCC
ENSG00000176953:ENST0000032080      GGTCGGGGCCGGCGTCCTCGGGCCCAGCGGTCTCCATCCCGGGGCACGCT

isotig35100                         AGACACGGTGCTTGTGGACTTGGTCAGTGACAGCGACGAGGAGATCCTGG
isotig35101                         AGACACGGTGCTTGTGGACTTGGTCAGTGACAGCGACGAGGAGATCCTGG
ENSG00000176953:ENST0000032080      GGACGTAGTGTCTGTGGACTTGGTCACCGACAGCGATGAGGAAATTCTGG

isotig35100                         AAGTC---------------GTCGCGGACCCGGTAGAAGCGCCCGCCGCC
isotig35101                         AAGTC---------------GTCGCGGACCCGGTAGAAGCGCCCGCCGCC
ENSG00000176953:ENST0000032080      AGGTCGCCACCGCTCGCGGTGCCGCGGACGAGGTTGAGGTGGAGCCCCCG

isotig35100                         CGGGCcCCcGCGCCGGCCGCACATGGGCAGGACAGCGACAGCGACAGTGC
isotig35101                         CGGGCcCCcGCGCCGGCCGCACATGGGCAGGACAGCGACAGCGACAGTGC
ENSG00000176953:ENST0000032080      GAGCCCCCGGGGCCGGTCGCGTCCCGGGATAACAGCAACAGTGACAGCGA

isotig35100                         AGGGGCGgACGAGGGGCCTGCAGGAGCCCCTCAGACCTTGGTCCGGCGGC
isotig35101                         AGGGGCGgACGAGGGGCCTGCAGGAGCCCCTCAGACCTTGGTCCGGCGGC
ENSG00000176953:ENST0000032080      AGGGGAGGACAGGCGGCCCGCAGGACCCCCGCGGGAGCCGGTCAGGCGGC

isotig35100                         GGCGCCGGCGGCTGCTGGATCCCGGCGAGGCACCGGTGGTTCCTGTGTAC
isotig35101                         GGCGCCGGCGGCTGCTGGATCCCGGCGAGGCACCGGTGGTTCCTGTGTAC
ENSG00000176953:ENST0000032080      GGCGGCGGCTGGTGCTGGATCCGGGGGAGGCGCCGCTGGTTCCGGTGTAC

isotig35100                         TCCGGGAAGGTACAAAGCAGCCTCAACCTCATCCCAGATAATTCATCCCT
isotig35101                         TCCGGGAAGGTACAAAGCAGCCTCAACCTCATCCCAGATAATTCATCCCT
ENSG00000176953:ENST0000032080      TCGGGGAAGGTTAAAAGCAGCCTTCGCCTTATCCCAGATGATCTATCCCT

isotig35100                         CTTGAAACTTTGCCCCTCAGAGCCTGAAGATGAGGCAGATGTGACAGATT
isotig35101                         CTTGAAACTTTGCCCCTCAGAGCCTGAAGATGAGGCAGATGTGACAGATT
ENSG00000176953:ENST0000032080      CCTGAAACTCTACCCTCCAGGGGATGAGGAAGAGGCAGAGCTGGCAGATT

isotig35100                         GTGGCAGTCCTCCTCCTGAGGATGCCCTAATTCCAGGTTCTCCCTGGAAG
isotig35101                         GTGGCAGTCCTCCTCCTGAGGATGCCCTAATTCCAGGTTCTCCCTGGAAG
ENSG00000176953:ENST0000032080      CGAGTGGTCTCTACCATGAGGGCTCCCCATCACCAGGCTCTCCCTGGAAG

isotig35100                         AAGAAGCTGAGGAATAAGCATGAAAAaGAAGAGATGAAGATGGAAGAGTT
isotig35101                         AAGAAGCTGAGGAATAAGCATGAAAAaGAAGAGATGAAGATGGAAGAGTT
ENSG00000176953:ENST0000032080      ACAAAGCTGAGGACTAAG---GATAAAGAAGAGAAGAAAAAGACAGAGTT

isotig35100                         tCCGGACCAGGACATCTCTCCTTTGCCCCGACCTTCATCAAGAAaCAAAA
isotig35101                         tCCGGACCAGGACATCTCTCCTTTGCCCCGACCTTCATCAAGAAaCAAAA
ENSG00000176953:ENST0000032080      TCTGGATCTGGACAACTCTCCTCTGTCCCCACCTTCACCAAGGACCAAAA

isotig35100                         GCAGAAAGCATACCGAGGCACTCCAGAAGTTGAGGGAAGTGAACAAGCGC
isotig35101                         GCAGAAAGCATACCGAGGCACTCCAGAAGTTGAGGGAAGTGAACAAGCGC
ENSG00000176953:ENST0000032080      GCAGAACGCATACTCGGGCACTCAAGAAGTTAAGTGAGGTGAACAAGCGC

isotig35100                         CTCCAAGATCTCCGATCCTGCCTGAGCCCCAAGCAGCACCAGAGTCCAGC
isotig35101                         CTCCAAGATCTCCGATCCTGCCTGAGCCCCAAGCAGCACCAGAGTCCAGC
ENSG00000176953:ENST0000032080      CTCCAGGATCTCCGTTCCTGTCTGAGCCCCAAGCCACCTCAGGGTCAAGA

isotig35100                         CCTTCAGAACCCAGATGATGAGGTGGTCCTCGTGGACGGGCCTGtCTTGT
isotig35101                         CCTTCAGAACCCAGATGATGAGGTGGTCCTCGTGGACGGGCCTGtCTTGT
ENSG00000176953:ENST0000032080      GCAACAGGGCCAAGAGGATGAAGTGGTCTTGGTGGAAGGGCCCACCCTCC

isotig35100                         CACAGAGCCCGAGACTCTTCACCCTCAAGATCCGGTGCCGGGCTGACCTA
isotig35101                         CACAGAGCCCGAGACTCTTCACCCTCAAGATCCGGTGCCGGGCTGACCTA
ENSG00000176953:ENST0000032080      CAGAGACCCCCCGACTCTTCCCACTCAAAATCCGTTGCCGGGCTGACCTG

isotig35100                         GTCAGATTGCCCGTCATGACATCGGAACCCCTTCAGAATGTGGTGGATTA
isotig35101                         GTCAGATTGCCCGTCATGACATCGGAACCCCTTCAGAATGTGGTGGATTA
ENSG00000176953:ENST0000032080      GTCAGATTGCCCCTCAGGATGTCGGAGCCCCTGCAGAGTGTGGTGGACCA

isotig35100                         CATGGCCAATCATCTTGGGGTGTCTCCAAGCAGGATTCTTTTACTCTTTG
isotig35101                         CATGGCCAATCATCTTGGGGTGTCTCCAAGCAGGATTCTTTTACTCTTTG
ENSG00000176953:ENST0000032080      CATGGCCACCCACCTTGGGGTGTCCCCAAGCAGGATCCTTTTGCTTTTTG

isotig35100                         GAGAGACAGAACTGTCCCCTACTGCCACCCCTAGGACCCTAAAGCTTGGT
isotig35101                         GAGAGACAGAACTGTCCCCTACTGCCACCCCTAGGACCCTAAAGCTTGGT
ENSG00000176953:ENST0000032080      GAGAGACAGAGCTATCACCTACTGCCACTCCCAGGACCCTAAAGCTCGGA

isotig35100                         GTGGCTGACATCATTGATTGTGTGGTGCTGACAAGTTCTTCAGAGGCCAC
isotig35101                         GTGGCTGACATCATTGATTGTGTGGTGCTGACAAGTTCTTCAGAGGCCAC
ENSG00000176953:ENST0000032080      GTGGCTGACATCATTGACTGTGTGGTACTAACAAGTTCTCCAGAGGCCAC

isotig35100                         AGAGACAACCCAGCAGCTCTGCCTCCGGGTGCAGGGGAAGGAGAAGCACC
isotig35101                         AGAGACAACCCAGCAGCTCTGCCTCCGGGTGCAGGGGAAGGAGAAGCACC
ENSG00000176953:ENST0000032080      AGAGACGTCCCAACAGCTCCAGCTCCGGGTGCAGGGAAAGGAGAAACACC

isotig35100                         AGATGTTGGAGATCTCACTCTCTCCTGACTCtCCTCTTGAGGTCCTCATG
isotig35101                         AGATGTTGGAGATCTCACTCTCTCCTGACTCgCCTCTTGAGGTCCTCATG
ENSG00000176953:ENST0000032080      AGACACTGGAAGTCTCACTGTCTCGAGATTCCCCTCTAAAGACCCTCATG

isotig35100                         GCGCACTATGAGGAGGCCATGGGACTCTCCGGACACAAGCTCTCCTTCTT
isotig35101                         GCGCACTATGAGGAGGCCATGGGACTCTCCGGACACAAGCTCTCCTTCTT
ENSG00000176953:ENST0000032080      TCCCACTATGAGGAGGCCATGGGACTGTCGGGACGGAAGCTCTCCTTCTT

isotig35100                         CTTCGATGGGACAAAGCTGTCGGGCAAGGAGCTGCCAGCTGATCTGGGCA
isotig35101                         CTTCGATGGGACAAAGCTGTCGGGCAAGGAACTGCCAGCTGATCTGGGCA
ENSG00000176953:ENST0000032080      CTTTGATGGGACAAAGCTTTCAGGCAGGGAGCTGCCAGCTGACCTGGGCA

isotig35100                         TGGAATCCGGGGATCTCATTGAAGTTTGGGGCAGCTTCCTCCTCCTGTTT
isotig35101                         TGGAATCCGGGGATCTCATTGAAGTTTGGGGCAGCTTCCTCCTCCTGTTT
ENSG00000176953:ENST0000032080      TGGAATCTGGGGACCTCATTGAGGTCTGGGGC------------------

isotig35100                         GGATGCAGAGCCAAGACTTGGGGACAACAGCTCCCACTTTtATTATTATT
isotig35101                         GGATGCAGAGCCAAGACTTGGGGACAACAGCTCCCACTTTTATTATTATT
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         TTTtGCCCcAGGGCTAACAGAAACCGAATTAGAACTCGTTTATTTATTTC
isotig35101                         TTTTGCCCCAGGGCTAACAGAAACCGAATTAGAACTCGTTTATTTATTTC
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         CGGCACTGGGGATTGAACCCAGGGCTGTGCATATGCTAAGGATGTGTGAA
isotig35101                         CGGCACTGGGGATTGAACCCTGGGCTGTGCATATGCTAAGGATGTGTGAA
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         GTTGAGGTAAAaCCAAGGCATGACCTTTGCCcTGTCTCGTTGACC-----
isotig35101                         GTTGAGGTAAAACCAAGGCATGACCTTTGCCCTGTCTCGTTGACCGATTG
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         ----------------------------------GTAGTCTCAAGCGGTC
isotig35101                         GTTGTGTGACTGTGGCTGCCCTGGGTGGCCTGTGGCTGTGTGTGTTGGTG
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         TGATTGGTAATTGTGTGACTGTGGCTGCCCTGGGTGgCCTGTGGCTGTGT
isotig35101                         CTGTGGACAGCAGCTCCTGGGGCATGGAGAAGGGGTATTGGCTTCCCTAC
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         GTGTTGGTGCTGTGTACAGCAGCTCCTGGGGCATGGAGAAGGGGTATTGG
isotig35101                         CATTTCGTTCAGTGCACAAATGAGCCTTGCATTTGGGTGCCCAGGCTTGT
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         CTTCCCTACCATTTCGTTCAGTAGTGCACAAATGAGCCTTGCATTTGGGT
isotig35101                         GATTTATGCCACATTGGGGACCAACAGGGATTTTAATTCTCATTTGGGGG
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         GCCCAGGCTTGTTTATGCCACATTGGGGACCAACAGGGATT---------
isotig35101                         CTGAGATGCAGGTACCTTCCCTGGTGGGGATAGTCGGTTGCAGGGAAATA
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         ------------------------------------TTAATTCTCATTTG
isotig35101                         AACAAGCAACATCTATGCAGCGGACTAGGTAGACAGTGGGTAATTTTTAA
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         GGGGCTGAGATGCAGGTACCTTCCCTGGTGGGGATCGGTTGCAGGGAAAA
isotig35101                         AGGCCAGGAGTGGTGTCTTCAGCCCTGGGCATCCACTCTCCCAGGACTCG
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         CAAGCAATGACATCTATGCTGAAGCTGAGGACGTAGACAGTGGGTTTTTA
isotig35101                         GCCAGCATCTGGTACATACTTCCTGGCATGCACATCTAGCCACATCTGAG
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         AAGGCTAACAGGAGTGGTGTCTTCAGCCCTGGGCATCCACTCTCCCAGGA
isotig35101                         CACATGCAGCTGTTGAGTTTGATGTTTACTCTAGTTCcGTCCTCCCGCCC
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         C------TCGGCCAGCATCTGGTACATACTT---------------CCTG
isotig35101                         TCTCCTGTCAGCACCCACTTGGTATATTTGTATCATGTCTTTAGATGATT
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         GCATGCACATCCCACATCTGAGCACATGCAGCTGTGTTTGTTTACTCTTC
isotig35101                         TTCTGAATGTGTGTGATTTTAGTGTATGTACATGGTACTGTGCAATAGAC
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         cGTCCTCCCGCCCTCTCCTGTCAGCACCCACTTGGTATATTTGTATCTGA
isotig35101                         ATTCTGTCTTCACTCAGTGTGAGTGACAGGCTAGCCATGTTAGCACCTTC
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         ATGTCTTTATAAGATGATTTTCATGTGTGTGATTTTAGTGTATGTACATG
isotig35101                         T-------------------------------------------------
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         GTACTGTGCAAACATTCTGTCTTCACTCAGTGGTGACAGGCCCATGTCAC
isotig35101                         --------AGATCTGTTGTTGTTCCGATCTGCTGTGGCTGCCCCATAGTG
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         CTTCTAGATCTGTTGTTGTTCCGATCTGCTGTGGCTGCCCCATGGCATAT
isotig35101                         GCATTAGATGATATTACTCACACATTGCCTTTGGAGGGATACCTGAGTCG
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         GATATTACTCACACATTGCCTTAATTGGAGGGATACCGTCGCCAGCTCCG
isotig35101                         CCAGCTCCTAGGAACCAAATTCCAGACCACAGCAAACATCTTTTTGTGTG
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         AACCAAATTCCAGACCACAGCAAACATCTTTTTGTGTGTTCCCTTCTGGG
isotig35101                         TTCCCTTCTGGGCCAGGGTAGGTAATTGTCTGGCATACTTGCCCAGGATG
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         CCAGGGGTAATTGTCTGGCATACTTGCCCAGGATGGGACTgTGGGTTACG
isotig35101                         GGACTgTGGGTTACGCAGGTTCATTTcTCTGGTTTTTCTAGGgTTTTCTT
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         CAGGTTCATTTcTCTGGTTTTTCTAGGgTTTTCTTCAGAATGTGGACTGC
isotig35101                         CAGAATGTGGACT---GCTTTTCTTGGAGCTTCTcTTTCCTTcTCTCCAG
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         TTTTCTTGGAGCTTCTcTTTCCTTcTCTTGACCAGTATTTAGTGTTGttg
isotig35101                         TATTTAGTGTTGttgTACTTTCTGATTTGCCAGTCTAC------------
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         TACTTTCTGATTTGCCAGTCTACTGTGTATAAAGCACTGCATATCTTAAT
isotig35101                         ------------TGTGTAAGCACTGCATATCTTAATATCCTTTATCTGAC
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         ATCCTTTATCCCACCATTTCTTATCAATGTCCACCTTCTGGTTTTTCTTT
isotig35101                         CACCATTTCTTATCAATGTCCACCTTCTGGTTTTTCTTTGTTAGCTGCTT
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         GTCTGCTTGTTTTCTTGGGCTTCCTTTTGTGTGTTTTTtCTGTGAGGAGT
isotig35101                         GTTTTCTTGTGAGGCTTCCTTTTGTG------TGTTTTTtCTGTGAGGAG
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         TGCTTGTGTTCTaGgCAATTGCAAATAACCTTCAGAATAGCACTTTGTGT
isotig35101                         TTGCTTGTGTTCgCAATTGCAAATAACCTTCAGAATAGCACTTTGTGTAC
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         ACCCAGGGGGATTTTCACCAAACAGACATTCTTAACTTTAATGCAGACAG
isotig35101                         CCAGGGGGATTTTCACCAAACAGACATTCTTAACTTTAATGCAGACAGCT
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         CTGGCTTTTGTAAAGGGTTTTTTGTTTTGTGTTTtTTTTTtttttGCTTT
isotig35101                         GGCTTTTGTTAGAAAGGGTTTTTTGTTTTGTGTTTtTTTTTtttttGCTT
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         TGTTTTtGcTTTTTAAGGGAAAATGAGAAAGAACTTCCAGAATGGCAGGG
isotig35101                         TTGTTTTtGcTTTTTAAGGGAAAAGAAAGAACT-----------------
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         GTTCCGGGGAGTCTCAGGATTGCTACTTTTGGATCTCACTAGTGAAACTT
isotig35101                         ----TCCAGAATGGCAGGGGTTCTAGCGGGGAGTCTCAGGATTGCTACTT
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         TCATTGTCTTCTGTCTTGCTAGTTTTTTGTTCTTTTGTTTTTATCTCTCA
isotig35101                         TTGGATCTCTAAACTTAAGAAAC------TTTCATTGTCTTCTGTTAGCT
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         TATGACCACATGAGAATAGAGTTCCTGATCTGCTGTGGTGTGGATaCTGA
isotig35101                         TGCTTTTTTGTTCTTTTGTTTTTATCTCTCATATTGAGACCACTAAATGA
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         CTGTTCTGAaCTTTCCACTCTGCTTCTCACTGGTTACAGCACTGTACTGG
isotig35101                         GAATAGAGTTCTAACTGATCTGCTGTGGTGTGGATaCCTGTTCaCTTTCC
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         TTCCTGACTATTGTATGTtACGTCTAAACCATCTGAAATGCACCTTTGTG
isotig35101                         ACTCTGCTTCTCACTGGTTACAGCTAGACTGTACTGGTTCCCTATTGTAT
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         TGGGAGAGCAAATTTTATTGTGTTTTAATTATAATGAGCCAGTTTTCtAA
isotig35101                         GTtACGTCTAAACCATCTGAAATGCACCTTTGTGTGGGAGAGCAAATTTT
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         GCTTTGTCTaCTAAAAAATTATGTAGTCCTTTCACAAGAAGTATGTCAAC
isotig35101                         ATTGTGTTTTAATTATAATGAGCCAGTTTTCGCTTTGTCTaCTAAAAAAT
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         TTCTtAAAAATTTTTTTtGtGAGAATATGTTATAC---------CTTATT
isotig35101                         TAATATGTAGTCTGACTTTCACAAGAAGTATGTCAACTTCTAAATTTTTT
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         TTAATATTTAAAAATTCCATGACATATAGCTTCCAGGTTTCTCAGCTCTT
isotig35101                         TtGGAATATGTTATACCTTTGAATTTTAATATTAAATTCCATGACATAGT
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         GTATGACAATATCACACATTaCTATTATGTATTCAATCTCAGCCAAAGGG
isotig35101                         ACTTCCAGGTTTCTCAGCTCTTGTACAATTAATGAATCTGAACACATTaC
ENSG00000176953:ENST0000032080      --------------------------------------------------

isotig35100                         AGATGCCTTTA---------------------------------
isotig35101                         TATTATGTATTCAATCTCAGCCAATAAAGGGAGATGCCTTTAAC
ENSG00000176953:ENST0000032080      --------------------------------------------


""",
        )

    def test3(self):
        aligner = CodonAligner()
        nucleotide_records = SeqIO.index("codonalign/nucl3.fa", "fasta")
        protein_alignment = Align.read("codonalign/pro3.aln", "clustal")
        self.assertEqual(len(protein_alignment.sequences), 10)
        codon_alignments = []
        protein_record = protein_alignment.sequences[0]
        nucleotide_record = nucleotide_records[protein_record.id]
        self.assertEqual(nucleotide_record.id, protein_record.id)
        alignments = aligner.align(protein_record, nucleotide_record)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 847], [0, 2541]]))
        )
        self.assertEqual(
            str(alignment),
            """\
ENSG00000         0 M  S  L  W  G  L  V  S  K  M  P  P  E  K  V  Q  R  L  Y  V  
ENSG00000         0 ATGTCTCTGTGGGGTCTGGTCTCCAAGATGCCCCCAGAAAAAGTGCAGCGGCTCTATGTC

ENSG00000        20 D  F  P  Q  H  L  R  H  L  L  G  D  W  L  E  S  Q  P  W  E  
ENSG00000        60 GACTTTCCCCAACACCTGCGGCATCTTCTGGGTGACTGGCTGGAGAGCCAGCCCTGGGAG

ENSG00000        40 F  L  V  G  S  D  A  F  C  C  N  L  A  S  A  L  L  S  D  T  
ENSG00000       120 TTCCTGGTCGGCTCCGACGCCTTCTGCTGCAACTTGGCTAGTGCCCTACTTTCAGACACT

ENSG00000        60 V  Q  H  L  Q  A  S  V  G  E  Q  G  E  G  S  T  I  L  Q  H  
ENSG00000       180 GTCCAGCACCTTCAGGCCTCGGTGGGAGAGCAGGGGGAGGGGAGCACCATCTTGCAACAC

ENSG00000        80 I  S  T  L  E  S  I  Y  Q  R  D  P  L  K  L  V  A  T  F  R  
ENSG00000       240 ATCAGCACCCTTGAGAGCATATATCAGAGGGACCCCCTGAAGCTGGTGGCCACTTTCAGA

ENSG00000       100 Q  I  L  Q  G  E  K  K  A  V  M  E  Q  F  R  H  L  P  M  P  
ENSG00000       300 CAAATACTTCAAGGAGAGAAAAAAGCTGTTATGGAACAGTTCCGCCACTTGCCAATGCCT

ENSG00000       120 F  H  W  K  Q  E  E  L  K  F  K  T  G  L  R  R  L  Q  H  R  
ENSG00000       360 TTCCACTGGAAGCAGGAAGAACTCAAGTTTAAGACAGGCTTGCGGAGGCTGCAGCACCGA

ENSG00000       140 V  G  E  I  H  L  L  R  E  A  L  Q  K  G  A  E  A  G  Q  V  
ENSG00000       420 GTAGGGGAGATCCACCTTCTCCGAGAAGCCCTGCAGAAGGGGGCTGAGGCTGGCCAAGTG

ENSG00000       160 S  L  H  S  L  I  E  T  P  A  N  G  T  G  P  S  E  A  L  A  
ENSG00000       480 TCTCTGCACAGCTTGATAGAAACTCCTGCTAATGGGACTGGGCCAAGTGAGGCCCTGGCC

ENSG00000       180 M  L  L  Q  E  T  T  G  E  L  E  A  A  K  A  L  V  L  K  R  
ENSG00000       540 ATGCTACTGCAGGAGACCACTGGAGAGCTAGAGGCAGCCAAAGCCCTAGTGCTGAAGAGG

ENSG00000       200 I  Q  I  W  K  R  Q  Q  Q  L  A  G  N  G  A  P  F  E  E  S  
ENSG00000       600 ATCCAGATTTGGAAACGGCAGCAGCAGCTGGCAGGGAATGGCGCACCGTTTGAGGAGAGC

ENSG00000       220 L  A  P  L  Q  E  R  C  E  S  L  V  D  I  Y  S  Q  L  Q  Q  
ENSG00000       660 CTGGCCCCACTCCAGGAGAGGTGTGAAAGCCTGGTGGACATTTATTCCCAGCTACAGCAG

ENSG00000       240 E  V  G  A  A  G  G  E  L  E  P  K  T  R  A  S  L  T  G  R  
ENSG00000       720 GAGGTAGGGGCGGCTGGTGGGGAGCTTGAGCCCAAGACCCGGGCATCGCTGACTGGCCGG

ENSG00000       260 L  D  E  V  L  R  T  L  V  T  S  C  F  L  V  E  K  Q  P  P  
ENSG00000       780 CTGGATGAAGTCCTGAGAACCCTCGTCACCAGTTGCTTCCTGGTGGAGAAGCAGCCCCCC

ENSG00000       280 Q  V  L  K  T  Q  T  K  F  Q  A  G  V  R  F  L  L  G  L  R  
ENSG00000       840 CAGGTACTGAAGACTCAGACCAAGTTCCAGGCTGGAGTTCGATTCCTGTTGGGCTTGAGG

ENSG00000       300 F  L  G  A  P  A  K  P  P  L  V  R  A  D  M  V  T  E  K  Q  
ENSG00000       900 TTCCTGGGGGCCCCAGCCAAGCCTCCGCTGGTCAGGGCCGACATGGTGACAGAGAAGCAG

ENSG00000       320 A  R  E  L  S  V  P  Q  G  P  G  A  G  A  E  S  T  G  E  I  
ENSG00000       960 GCGCGGGAGCTGAGTGTGCCTCAGGGTCCTGGGGCTGGAGCAGAAAGCACTGGAGAAATC

ENSG00000       340 I  N  N  T  V  P  L  E  N  S  I  P  G  N  C  C  S  A  L  F  
ENSG00000      1020 ATCAACAACACTGTGCCCTTGGAGAACAGCATTCCTGGGAACTGCTGCTCTGCCCTGTTC

ENSG00000       360 K  N  L  L  L  K  K  I  K  R  C  E  R  K  G  T  E  S  V  T  
ENSG00000      1080 AAGAACCTGCTTCTCAAGAAGATCAAGCGGTGTGAGCGGAAGGGCACTGAGTCTGTCACA

ENSG00000       380 E  E  K  C  A  V  L  F  S  A  S  F  T  L  G  P  G  K  L  P  
ENSG00000      1140 GAGGAGAAGTGCGCTGTGCTCTTCTCTGCCAGCTTCACACTTGGCCCCGGCAAACTCCCC

ENSG00000       400 I  Q  L  Q  A  L  S  L  P  L  V  V  I  V  H  G  N  Q  D  N  
ENSG00000      1200 ATCCAGCTCCAGGCCCTGTCTCTGCCCCTGGTGGTCATCGTCCATGGCAACCAAGACAAC

ENSG00000       420 N  A  K  A  T  I  L  W  D  N  A  F  S  E  M  D  R  V  P  F  
ENSG00000      1260 AATGCCAAAGCCACTATCCTGTGGGACAATGCCTTCTCTGAGATGGACCGCGTGCCCTTT

ENSG00000       440 V  V  A  E  R  V  P  W  E  K  M  C  E  T  L  N  L  K  F  M  
ENSG00000      1320 GTGGTGGCTGAGCGGGTGCCCTGGGAGAAGATGTGTGAAACTCTGAACCTGAAGTTCATG

ENSG00000       460 A  E  V  G  T  N  R  G  L  L  P  E  H  F  L  F  L  A  Q  K  
ENSG00000      1380 GCTGAGGTGGGGACCAACCGGGGGCTGCTCCCAGAGCACTTCCTCTTCCTGGCCCAGAAG

ENSG00000       480 I  F  N  D  N  S  L  S  M  E  A  F  Q  H  R  S  V  S  W  S  
ENSG00000      1440 ATCTTCAATGACAACAGCCTCAGTATGGAGGCCTTCCAGCACCGTTCTGTGTCCTGGTCG

ENSG00000       500 Q  F  N  K  E  I  L  L  G  R  G  F  T  F  W  Q  W  F  D  G  
ENSG00000      1500 CAGTTCAACAAGGAGATCCTGCTGGGCCGTGGCTTCACCTTTTGGCAGTGGTTTGATGGT

ENSG00000       520 V  L  D  L  T  K  R  C  L  R  S  Y  W  S  D  R  L  I  I  G  
ENSG00000      1560 GTCCTGGACCTCACCAAACGCTGTCTCCGGAGCTACTGGTCTGACCGGCTGATCATTGGC

ENSG00000       540 F  I  S  K  Q  Y  V  T  S  L  L  L  N  E  P  D  G  T  F  L  
ENSG00000      1620 TTCATCAGCAAACAGTACGTTACTAGCCTTCTTCTCAATGAGCCCGACGGAACCTTTCTC

ENSG00000       560 L  R  F  S  D  S  E  I  G  G  I  T  I  A  H  V  I  R  G  Q  
ENSG00000      1680 CTCCGCTTCAGCGACTCAGAGATTGGGGGCATCACCATTGCCCATGTCATCCGGGGCCAG

ENSG00000       580 D  G  S  P  Q  I  E  N  I  Q  P  F  S  A  K  D  L  S  I  R  
ENSG00000      1740 GATGGCTCTCCACAGATAGAGAACATCCAGCCATTCTCTGCCAAAGACCTGTCCATTCGC

ENSG00000       600 S  L  G  D  R  I  R  D  L  A  Q  L  K  N  L  Y  P  K  K  P  
ENSG00000      1800 TCACTGGGGGACCGAATCCGGGATCTTGCTCAGCTCAAAAATCTCTATCCCAAGAAGCCC

ENSG00000       620 K  D  E  A  F  R  S  H  Y  K  P  E  Q  M  G  K  D  G  R  G  
ENSG00000      1860 AAGGATGAGGCTTTCCGGAGCCACTACAAGCCTGAACAGATGGGTAAGGATGGCAGGGGT

ENSG00000       640 Y  V  P  A  T  I  K  M  T  V  E  R  D  Q  P  L  P  T  P  E  
ENSG00000      1920 TATGTCCCAGCTACCATCAAGATGACCGTGGAAAGGGACCAACCACTTCCTACCCCAGAG

ENSG00000       660 L  Q  M  P  T  M  V  P  S  Y  D  L  G  M  A  P  D  S  S  M  
ENSG00000      1980 CTCCAGATGCCTACCATGGTGCCTTCTTATGACCTTGGAATGGCCCCTGATTCCTCCATG

ENSG00000       680 S  M  Q  L  G  P  D  M  V  P  Q  V  Y  P  P  H  S  H  S  I  
ENSG00000      2040 AGCATGCAGCTTGGCCCAGATATGGTGCCCCAGGTGTACCCACCACACTCTCACTCCATC

ENSG00000       700 P  P  Y  Q  G  L  S  P  E  E  S  V  N  V  L  S  A  F  Q  E  
ENSG00000      2100 CCCCCGTATCAAGGCCTCTCCCCAGAAGAATCAGTCAACGTGTTGTCAGCCTTCCAGGAG

ENSG00000       720 P  H  L  Q  M  P  P  S  L  G  Q  M  S  L  P  F  D  Q  P  H  
ENSG00000      2160 CCTCACCTGCAGATGCCCCCCAGCCTGGGCCAGATGAGCCTGCCCTTTGACCAGCCTCAC

ENSG00000       740 P  Q  G  L  L  P  C  Q  P  Q  E  H  A  V  S  S  P  D  P  L  
ENSG00000      2220 CCCCAGGGCCTGCTGCCGTGCCAGCCTCAGGAGCATGCTGTGTCCAGCCCTGACCCCCTG

ENSG00000       760 L  C  S  D  V  T  M  V  E  D  S  C  L  S  Q  P  V  T  A  F  
ENSG00000      2280 CTCTGCTCAGATGTGACCATGGTGGAAGACAGCTGCCTGAGCCAGCCAGTGACAGCGTTT

ENSG00000       780 P  Q  G  T  W  I  G  E  D  I  F  P  P  L  L  P  P  T  E  Q  
ENSG00000      2340 CCTCAGGGCACTTGGATTGGTGAAGACATATTCCCTCCTCTGCTGCCTCCCACTGAACAG

ENSG00000       800 D  L  T  K  L  L  L  E  G  Q  G  E  S  G  G  G  S  L  G  A  
ENSG00000      2400 GACCTCACTAAGCTTCTCCTGGAGGGGCAAGGGGAGTCGGGGGGAGGGTCCTTGGGGGCA

ENSG00000       820 Q  P  L  L  Q  P  S  H  Y  G  Q  S  G  I  S  M  S  H  M  D  
ENSG00000      2460 CAGCCCCTCCTGCAGCCCTCCCACTATGGGCAATCTGGGATCTCAATGTCCCACATGGAC

ENSG00000       840 L  R  A  N  P  S  W    847
ENSG00000      2520 CTAAGGGCCAACCCCAGTTGG 2541
""",
        )
        protein_record = protein_alignment.sequences[1]
        nucleotide_record = nucleotide_records[protein_record.id]
        self.assertEqual(nucleotide_record.id, protein_record.id)
        alignments = aligner.align(protein_record, nucleotide_record)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 847], [0, 2541]]))
        )
        self.assertEqual(
            str(alignment),
            """\
ENSG00000         0 M  S  L  W  G  L  V  S  K  M  P  P  E  K  V  Q  R  L  Y  V  
ENSG00000         0 ATGTCTCTGTGGGGTCTGGTCTCCAAGATGCCCCCAGAAAAAGTGCAGCGGCTCTATGTC

ENSG00000        20 D  F  P  Q  H  L  R  H  L  L  G  D  W  L  E  S  Q  P  W  E  
ENSG00000        60 GACTTTCCCCAACACCTGCGGCATCTTCTGGGTGACTGGCTGGAGAGCCAGCCCTGGGAG

ENSG00000        40 F  L  V  G  S  D  A  F  C  C  N  L  A  S  A  L  L  S  D  T  
ENSG00000       120 TTCCTGGTCGGCTCCGACGCCTTCTGCTGCAACTTGGCTAGTGCCCTACTTTCAGACACT

ENSG00000        60 V  Q  H  L  Q  A  S  V  G  E  Q  G  E  G  S  T  I  L  Q  H  
ENSG00000       180 GTCCAGCACCTTCAGGCCTCGGTGGGAGAGCAGGGGGAGGGGAGCACCATCTTGCAACAC

ENSG00000        80 I  S  T  L  E  S  I  Y  Q  R  D  P  L  K  L  V  A  T  F  R  
ENSG00000       240 ATCAGCACCCTTGAGAGCATATATCAGAGGGACCCCCTGAAGCTGGTGGCCACTTTCAGA

ENSG00000       100 Q  I  L  Q  G  E  K  K  A  V  M  E  Q  F  R  H  L  P  M  P  
ENSG00000       300 CAAATACTTCAAGGAGAGAAAAAAGCTGTTATGGAACAGTTCCGCCACTTGCCAATGCCT

ENSG00000       120 F  H  W  K  Q  E  E  L  K  F  K  T  G  L  R  R  L  Q  H  R  
ENSG00000       360 TTCCACTGGAAGCAGGAAGAACTCAAGTTTAAGACAGGCTTGCGGAGGCTGCAGCACCGA

ENSG00000       140 V  G  E  I  H  L  L  R  E  A  L  Q  K  G  A  E  A  G  Q  V  
ENSG00000       420 GTAGGGGAGATCCACCTTCTCCGAGAAGCCCTGCAGAAGGGGGCTGAGGCTGGCCAAGTG

ENSG00000       160 S  L  H  S  L  I  E  T  P  A  N  G  T  G  P  S  E  A  L  A  
ENSG00000       480 TCTCTGCACAGCTTGATAGAAACTCCTGCTAATGGGACTGGGCCAAGTGAGGCCCTGGCC

ENSG00000       180 M  L  L  Q  E  T  T  G  E  L  E  A  A  K  A  L  V  L  K  R  
ENSG00000       540 ATGCTACTGCAGGAGACCACTGGAGAGCTAGAGGCAGCCAAAGCCCTAGTGCTGAAGAGG

ENSG00000       200 I  Q  I  W  K  R  Q  Q  Q  L  A  G  N  G  A  P  F  E  E  S  
ENSG00000       600 ATCCAGATTTGGAAACGGCAGCAGCAGCTGGCAGGGAATGGCGCACCGTTTGAGGAGAGC

ENSG00000       220 L  A  P  L  Q  E  R  C  E  S  L  V  D  I  Y  S  Q  L  Q  Q  
ENSG00000       660 CTGGCCCCACTCCAGGAGAGGTGTGAAAGCCTGGTGGACATTTATTCCCAGCTACAGCAG

ENSG00000       240 E  V  G  A  A  G  G  E  L  E  P  K  T  R  A  S  L  T  G  R  
ENSG00000       720 GAGGTAGGGGCGGCTGGTGGGGAGCTTGAGCCCAAGACCCGGGCATCGCTGACTGGCCGG

ENSG00000       260 L  D  E  V  L  R  T  L  V  T  S  C  F  L  V  E  K  Q  P  P  
ENSG00000       780 CTGGATGAAGTCCTGAGAACCCTCGTCACCAGTTGCTTCCTGGTGGAGAAGCAGCCCCCC

ENSG00000       280 Q  V  L  K  T  Q  T  K  F  Q  A  G  V  R  F  L  L  G  L  R  
ENSG00000       840 CAGGTACTGAAGACTCAGACCAAGTTCCAGGCTGGAGTTCGATTCCTGTTGGGCTTGAGG

ENSG00000       300 F  L  G  A  P  A  K  P  P  L  V  R  A  D  M  V  T  E  K  Q  
ENSG00000       900 TTCCTGGGGGCCCCAGCCAAGCCTCCGCTGGTCAGGGCCGACATGGTGACAGAGAAGCAG

ENSG00000       320 A  R  E  L  S  V  P  Q  G  P  G  A  G  A  E  S  T  G  E  I  
ENSG00000       960 GCGCGGGAGCTGAGTGTGCCTCAGGGTCCTGGGGCTGGAGCAGAAAGCACTGGAGAAATC

ENSG00000       340 I  N  N  T  V  P  L  E  N  S  I  P  G  N  C  C  S  A  L  F  
ENSG00000      1020 ATCAACAACACTGTGCCCTTGGAGAACAGCATTCCTGGGAACTGCTGCTCTGCCCTGTTC

ENSG00000       360 K  N  L  L  L  K  K  I  K  R  C  E  R  K  G  T  E  S  V  T  
ENSG00000      1080 AAGAACCTGCTTCTCAAGAAGATCAAGCGGTGTGAGCGGAAGGGCACTGAGTCTGTCACA

ENSG00000       380 E  E  K  C  A  V  L  F  S  A  S  F  T  L  G  P  G  K  L  P  
ENSG00000      1140 GAGGAGAAGTGCGCTGTGCTCTTCTCTGCCAGCTTCACACTTGGCCCCGGCAAACTCCCC

ENSG00000       400 I  Q  L  Q  A  L  S  L  P  L  V  V  I  V  H  G  N  Q  D  N  
ENSG00000      1200 ATCCAGCTCCAGGCCCTGTCTCTGCCCCTGGTGGTCATCGTCCATGGCAACCAAGACAAC

ENSG00000       420 N  A  K  A  T  I  L  W  D  N  A  F  S  E  M  D  R  V  P  F  
ENSG00000      1260 AATGCCAAAGCCACTATCCTGTGGGACAATGCCTTCTCTGAGATGGACCGCGTGCCCTTT

ENSG00000       440 V  V  A  E  R  V  P  W  E  K  M  C  E  T  L  N  L  K  F  M  
ENSG00000      1320 GTGGTGGCTGAGCGGGTGCCCTGGGAGAAGATGTGTGAAACTCTGAACCTGAAGTTCATG

ENSG00000       460 A  E  V  G  T  N  R  G  L  L  P  E  H  F  L  F  L  A  Q  K  
ENSG00000      1380 GCTGAGGTGGGGACCAACCGGGGGCTGCTCCCAGAGCACTTCCTCTTCCTGGCCCAGAAG

ENSG00000       480 I  F  N  D  N  S  L  S  M  E  A  F  Q  H  R  S  V  S  W  S  
ENSG00000      1440 ATCTTCAATGACAACAGCCTCAGTATGGAGGCCTTCCAGCACCGTTCTGTGTCCTGGTCG

ENSG00000       500 Q  F  N  K  E  I  L  L  G  R  G  F  T  F  W  Q  W  F  D  G  
ENSG00000      1500 CAGTTCAACAAGGAGATCCTGCTGGGCCGTGGCTTCACCTTTTGGCAGTGGTTTGATGGT

ENSG00000       520 V  L  D  L  T  K  R  C  L  R  S  Y  W  S  D  R  L  I  I  G  
ENSG00000      1560 GTCCTGGACCTCACCAAACGCTGTCTCCGGAGCTACTGGTCTGACCGGCTGATCATTGGC

ENSG00000       540 F  I  S  K  Q  Y  V  T  S  L  L  L  N  E  P  D  G  T  F  L  
ENSG00000      1620 TTCATCAGCAAACAGTACGTTACTAGCCTTCTTCTCAATGAGCCCGACGGAACCTTTCTC

ENSG00000       560 L  R  F  S  D  S  E  I  G  G  I  T  I  A  H  V  I  R  G  Q  
ENSG00000      1680 CTCCGCTTCAGCGACTCAGAGATTGGGGGCATCACCATTGCCCATGTCATCCGGGGCCAG

ENSG00000       580 D  G  S  P  Q  I  E  N  I  Q  P  F  S  A  K  D  L  S  I  R  
ENSG00000      1740 GATGGCTCTCCACAGATAGAGAACATCCAGCCATTCTCTGCCAAAGACCTGTCCATTCGC

ENSG00000       600 S  L  G  D  R  I  R  D  L  A  Q  L  K  N  L  Y  P  K  K  P  
ENSG00000      1800 TCACTGGGGGACCGAATCCGGGATCTTGCTCAGCTCAAAAATCTCTATCCCAAGAAGCCC

ENSG00000       620 K  D  E  A  F  R  S  H  Y  K  P  E  Q  M  G  K  D  G  R  G  
ENSG00000      1860 AAGGATGAGGCTTTCCGGAGCCACTACAAGCCTGAACAGATGGGTAAGGATGGCAGGGGT

ENSG00000       640 Y  V  P  A  T  I  K  M  T  V  E  R  D  Q  P  L  P  T  P  E  
ENSG00000      1920 TATGTCCCAGCTACCATCAAGATGACCGTGGAAAGGGACCAACCACTTCCTACCCCAGAG

ENSG00000       660 L  Q  M  P  T  M  V  P  S  Y  D  L  G  M  A  P  D  S  S  M  
ENSG00000      1980 CTCCAGATGCCTACCATGGTGCCTTCTTATGACCTTGGAATGGCCCCTGATTCCTCCATG

ENSG00000       680 S  M  Q  L  G  P  D  M  V  P  Q  V  Y  P  P  H  S  H  S  I  
ENSG00000      2040 AGCATGCAGCTTGGCCCAGATATGGTGCCCCAGGTGTACCCACCACACTCTCACTCCATC

ENSG00000       700 P  P  Y  Q  G  L  S  P  E  E  S  V  N  V  L  S  A  F  Q  E  
ENSG00000      2100 CCCCCGTATCAAGGCCTCTCCCCAGAAGAATCAGTCAACGTGTTGTCAGCCTTCCAGGAG

ENSG00000       720 P  H  L  Q  M  P  P  S  L  G  Q  M  S  L  P  F  D  Q  P  H  
ENSG00000      2160 CCTCACCTGCAGATGCCCCCCAGCCTGGGCCAGATGAGCCTGCCCTTTGACCAGCCTCAC

ENSG00000       740 P  Q  G  L  L  P  C  Q  P  Q  E  H  A  V  S  S  P  D  P  L  
ENSG00000      2220 CCCCAGGGCCTGCTGCCGTGCCAGCCTCAGGAGCATGCTGTGTCCAGCCCTGACCCCCTG

ENSG00000       760 L  C  S  D  V  T  M  V  E  D  S  C  L  S  Q  P  V  T  A  F  
ENSG00000      2280 CTCTGCTCAGATGTGACCATGGTGGAAGACAGCTGCCTGAGCCAGCCAGTGACAGCGTTT

ENSG00000       780 P  Q  G  T  W  I  G  E  D  I  F  P  P  L  L  P  P  T  E  Q  
ENSG00000      2340 CCTCAGGGCACTTGGATTGGTGAAGACATATTCCCTCCTCTGCTGCCTCCCACTGAACAG

ENSG00000       800 D  L  T  K  L  L  L  E  G  Q  G  E  S  G  G  G  S  L  G  A  
ENSG00000      2400 GACCTCACTAAGCTTCTCCTGGAGGGGCAAGGGGAGTCGGGGGGAGGGTCCTTGGGGGCA

ENSG00000       820 Q  P  L  L  Q  P  S  H  Y  G  Q  S  G  I  S  M  S  H  M  D  
ENSG00000      2460 CAGCCCCTCCTGCAGCCCTCCCACTATGGGCAATCTGGGATCTCAATGTCCCACATGGAC

ENSG00000       840 L  R  A  N  P  S  W    847
ENSG00000      2520 CTAAGGGCCAACCCCAGTTGG 2541
""",
        )
        protein_record = protein_alignment.sequences[2]
        nucleotide_record = nucleotide_records[protein_record.id]
        self.assertEqual(nucleotide_record.id, protein_record.id)
        alignments = aligner.align(protein_record, nucleotide_record)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 847], [0, 2541]]))
        )
        self.assertEqual(
            str(alignment),
            """\
ENSG00000         0 M  S  L  W  G  L  V  S  K  M  P  P  E  K  V  Q  R  L  Y  V  
ENSG00000         0 ATGTCTCTGTGGGGTCTGGTCTCCAAGATGCCCCCAGAAAAAGTGCAGCGGCTCTATGTC

ENSG00000        20 D  F  P  Q  H  L  R  H  L  L  G  D  W  L  E  S  Q  P  W  E  
ENSG00000        60 GACTTTCCCCAACACCTGCGGCATCTTCTGGGTGACTGGCTGGAGAGCCAGCCCTGGGAG

ENSG00000        40 F  L  V  G  S  D  A  F  C  C  N  L  A  S  A  L  L  S  D  T  
ENSG00000       120 TTCCTGGTCGGCTCCGACGCCTTCTGCTGCAACTTGGCTAGTGCCCTACTTTCAGACACT

ENSG00000        60 V  Q  H  L  Q  A  S  V  G  E  Q  G  E  G  S  T  I  L  Q  H  
ENSG00000       180 GTCCAGCACCTTCAGGCCTCGGTGGGAGAGCAGGGGGAGGGGAGCACCATCTTGCAACAC

ENSG00000        80 I  S  T  L  E  S  I  Y  Q  R  D  P  L  K  L  V  A  T  F  R  
ENSG00000       240 ATCAGCACCCTTGAGAGCATATATCAGAGGGACCCCCTGAAGCTGGTGGCCACTTTCAGA

ENSG00000       100 Q  I  L  Q  G  E  K  K  A  V  M  E  Q  F  R  H  L  P  M  P  
ENSG00000       300 CAAATACTTCAAGGAGAGAAAAAAGCTGTTATGGAACAGTTCCGCCACTTGCCAATGCCT

ENSG00000       120 F  H  W  K  Q  E  E  L  K  F  K  T  G  L  R  R  L  Q  H  R  
ENSG00000       360 TTCCACTGGAAGCAGGAAGAACTCAAGTTTAAGACAGGCTTGCGGAGGCTGCAGCACCGA

ENSG00000       140 V  G  E  I  H  L  L  R  E  A  L  Q  K  G  A  E  A  G  Q  V  
ENSG00000       420 GTAGGGGAGATCCACCTTCTCCGAGAAGCCCTGCAGAAGGGGGCTGAGGCTGGCCAAGTG

ENSG00000       160 S  L  H  S  L  I  E  T  P  A  N  G  T  G  P  S  E  A  L  A  
ENSG00000       480 TCTCTGCACAGCTTGATAGAAACTCCTGCTAATGGGACTGGGCCAAGTGAGGCCCTGGCC

ENSG00000       180 M  L  L  Q  E  T  T  G  E  L  E  A  A  K  A  L  V  L  K  R  
ENSG00000       540 ATGCTACTGCAGGAGACCACTGGAGAGCTAGAGGCAGCCAAAGCCCTAGTGCTGAAGAGG

ENSG00000       200 I  Q  I  W  K  R  Q  Q  Q  L  A  G  N  G  A  P  F  E  E  S  
ENSG00000       600 ATCCAGATTTGGAAACGGCAGCAGCAGCTGGCAGGGAATGGCGCACCGTTTGAGGAGAGC

ENSG00000       220 L  A  P  L  Q  E  R  C  E  S  L  V  D  I  Y  S  Q  L  Q  Q  
ENSG00000       660 CTGGCCCCACTCCAGGAGAGGTGTGAAAGCCTGGTGGACATTTATTCCCAGCTACAGCAG

ENSG00000       240 E  V  G  A  A  G  G  E  L  E  P  K  T  R  A  S  L  T  G  R  
ENSG00000       720 GAGGTAGGGGCGGCTGGTGGGGAGCTTGAGCCCAAGACCCGGGCATCGCTGACTGGCCGG

ENSG00000       260 L  D  E  V  L  R  T  L  V  T  S  C  F  L  V  E  K  Q  P  P  
ENSG00000       780 CTGGATGAAGTCCTGAGAACCCTCGTCACCAGTTGCTTCCTGGTGGAGAAGCAGCCCCCC

ENSG00000       280 Q  V  L  K  T  Q  T  K  F  Q  A  G  V  R  F  L  L  G  L  R  
ENSG00000       840 CAGGTACTGAAGACTCAGACCAAGTTCCAGGCTGGAGTTCGATTCCTGTTGGGCTTGAGG

ENSG00000       300 F  L  G  A  P  A  K  P  P  L  V  R  A  D  M  V  T  E  K  Q  
ENSG00000       900 TTCCTGGGGGCCCCAGCCAAGCCTCCGCTGGTCAGGGCCGACATGGTGACAGAGAAGCAG

ENSG00000       320 A  R  E  L  S  V  P  Q  G  P  G  A  G  A  E  S  T  G  E  I  
ENSG00000       960 GCGCGGGAGCTGAGTGTGCCTCAGGGTCCTGGGGCTGGAGCAGAAAGCACTGGAGAAATC

ENSG00000       340 I  N  N  T  V  P  L  E  N  S  I  P  G  N  C  C  S  A  L  F  
ENSG00000      1020 ATCAACAACACTGTGCCCTTGGAGAACAGCATTCCTGGGAACTGCTGCTCTGCCCTGTTC

ENSG00000       360 K  N  L  L  L  K  K  I  K  R  C  E  R  K  G  T  E  S  V  T  
ENSG00000      1080 AAGAACCTGCTTCTCAAGAAGATCAAGCGGTGTGAGCGGAAGGGCACTGAGTCTGTCACA

ENSG00000       380 E  E  K  C  A  V  L  F  S  A  S  F  T  L  G  P  G  K  L  P  
ENSG00000      1140 GAGGAGAAGTGCGCTGTGCTCTTCTCTGCCAGCTTCACACTTGGCCCCGGCAAACTCCCC

ENSG00000       400 I  Q  L  Q  A  L  S  L  P  L  V  V  I  V  H  G  N  Q  D  N  
ENSG00000      1200 ATCCAGCTCCAGGCCCTGTCTCTGCCCCTGGTGGTCATCGTCCATGGCAACCAAGACAAC

ENSG00000       420 N  A  K  A  T  I  L  W  D  N  A  F  S  E  M  D  R  V  P  F  
ENSG00000      1260 AATGCCAAAGCCACTATCCTGTGGGACAATGCCTTCTCTGAGATGGACCGCGTGCCCTTT

ENSG00000       440 V  V  A  E  R  V  P  W  E  K  M  C  E  T  L  N  L  K  F  M  
ENSG00000      1320 GTGGTGGCTGAGCGGGTGCCCTGGGAGAAGATGTGTGAAACTCTGAACCTGAAGTTCATG

ENSG00000       460 A  E  V  G  T  N  R  G  L  L  P  E  H  F  L  F  L  A  Q  K  
ENSG00000      1380 GCTGAGGTGGGGACCAACCGGGGGCTGCTCCCAGAGCACTTCCTCTTCCTGGCCCAGAAG

ENSG00000       480 I  F  N  D  N  S  L  S  M  E  A  F  Q  H  R  S  V  S  W  S  
ENSG00000      1440 ATCTTCAATGACAACAGCCTCAGTATGGAGGCCTTCCAGCACCGTTCTGTGTCCTGGTCG

ENSG00000       500 Q  F  N  K  E  I  L  L  G  R  G  F  T  F  W  Q  W  F  D  G  
ENSG00000      1500 CAGTTCAACAAGGAGATCCTGCTGGGCCGTGGCTTCACCTTTTGGCAGTGGTTTGATGGT

ENSG00000       520 V  L  D  L  T  K  R  C  L  R  S  Y  W  S  D  R  L  I  I  G  
ENSG00000      1560 GTCCTGGACCTCACCAAACGCTGTCTCCGGAGCTACTGGTCTGACCGGCTGATCATTGGC

ENSG00000       540 F  I  S  K  Q  Y  V  T  S  L  L  L  N  E  P  D  G  T  F  L  
ENSG00000      1620 TTCATCAGCAAACAGTACGTTACTAGCCTTCTTCTCAATGAGCCCGACGGAACCTTTCTC

ENSG00000       560 L  R  F  S  D  S  E  I  G  G  I  T  I  A  H  V  I  R  G  Q  
ENSG00000      1680 CTCCGCTTCAGCGACTCAGAGATTGGGGGCATCACCATTGCCCATGTCATCCGGGGCCAG

ENSG00000       580 D  G  S  P  Q  I  E  N  I  Q  P  F  S  A  K  D  L  S  I  R  
ENSG00000      1740 GATGGCTCTCCACAGATAGAGAACATCCAGCCATTCTCTGCCAAAGACCTGTCCATTCGC

ENSG00000       600 S  L  G  D  R  I  R  D  L  A  Q  L  K  N  L  Y  P  K  K  P  
ENSG00000      1800 TCACTGGGGGACCGAATCCGGGATCTTGCTCAGCTCAAAAATCTCTATCCCAAGAAGCCC

ENSG00000       620 K  D  E  A  F  R  S  H  Y  K  P  E  Q  M  G  K  D  G  R  G  
ENSG00000      1860 AAGGATGAGGCTTTCCGGAGCCACTACAAGCCTGAACAGATGGGTAAGGATGGCAGGGGT

ENSG00000       640 Y  V  P  A  T  I  K  M  T  V  E  R  D  Q  P  L  P  T  P  E  
ENSG00000      1920 TATGTCCCAGCTACCATCAAGATGACCGTGGAAAGGGACCAACCACTTCCTACCCCAGAG

ENSG00000       660 L  Q  M  P  T  M  V  P  S  Y  D  L  G  M  A  P  D  S  S  M  
ENSG00000      1980 CTCCAGATGCCTACCATGGTGCCTTCTTATGACCTTGGAATGGCCCCTGATTCCTCCATG

ENSG00000       680 S  M  Q  L  G  P  D  M  V  P  Q  V  Y  P  P  H  S  H  S  I  
ENSG00000      2040 AGCATGCAGCTTGGCCCAGATATGGTGCCCCAGGTGTACCCACCACACTCTCACTCCATC

ENSG00000       700 P  P  Y  Q  G  L  S  P  E  E  S  V  N  V  L  S  A  F  Q  E  
ENSG00000      2100 CCCCCGTATCAAGGCCTCTCCCCAGAAGAATCAGTCAACGTGTTGTCAGCCTTCCAGGAG

ENSG00000       720 P  H  L  Q  M  P  P  S  L  G  Q  M  S  L  P  F  D  Q  P  H  
ENSG00000      2160 CCTCACCTGCAGATGCCCCCCAGCCTGGGCCAGATGAGCCTGCCCTTTGACCAGCCTCAC

ENSG00000       740 P  Q  G  L  L  P  C  Q  P  Q  E  H  A  V  S  S  P  D  P  L  
ENSG00000      2220 CCCCAGGGCCTGCTGCCGTGCCAGCCTCAGGAGCATGCTGTGTCCAGCCCTGACCCCCTG

ENSG00000       760 L  C  S  D  V  T  M  V  E  D  S  C  L  S  Q  P  V  T  A  F  
ENSG00000      2280 CTCTGCTCAGATGTGACCATGGTGGAAGACAGCTGCCTGAGCCAGCCAGTGACAGCGTTT

ENSG00000       780 P  Q  G  T  W  I  G  E  D  I  F  P  P  L  L  P  P  T  E  Q  
ENSG00000      2340 CCTCAGGGCACTTGGATTGGTGAAGACATATTCCCTCCTCTGCTGCCTCCCACTGAACAG

ENSG00000       800 D  L  T  K  L  L  L  E  G  Q  G  E  S  G  G  G  S  L  G  A  
ENSG00000      2400 GACCTCACTAAGCTTCTCCTGGAGGGGCAAGGGGAGTCGGGGGGAGGGTCCTTGGGGGCA

ENSG00000       820 Q  P  L  L  Q  P  S  H  Y  G  Q  S  G  I  S  M  S  H  M  D  
ENSG00000      2460 CAGCCCCTCCTGCAGCCCTCCCACTATGGGCAATCTGGGATCTCAATGTCCCACATGGAC

ENSG00000       840 L  R  A  N  P  S  W    847
ENSG00000      2520 CTAAGGGCCAACCCCAGTTGG 2541
""",
        )
        protein_record = protein_alignment.sequences[3]
        nucleotide_record = nucleotide_records[protein_record.id]
        self.assertEqual(nucleotide_record.id, protein_record.id)
        alignments = aligner.align(protein_record, nucleotide_record)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 847], [0, 2541]]))
        )
        self.assertEqual(
            str(alignment),
            """\
ENSG00000         0 M  S  L  W  G  L  V  S  K  M  P  P  E  K  V  Q  R  L  Y  V  
ENSG00000         0 ATGTCTCTGTGGGGTCTGGTCTCCAAGATGCCCCCAGAAAAAGTGCAGCGGCTCTATGTC

ENSG00000        20 D  F  P  Q  H  L  R  H  L  L  G  D  W  L  E  S  Q  P  W  E  
ENSG00000        60 GACTTTCCCCAACACCTGCGGCATCTTCTGGGTGACTGGCTGGAGAGCCAGCCCTGGGAG

ENSG00000        40 F  L  V  G  S  D  A  F  C  C  N  L  A  S  A  L  L  S  D  T  
ENSG00000       120 TTCCTGGTCGGCTCCGACGCCTTCTGCTGCAACTTGGCTAGTGCCCTACTTTCAGACACT

ENSG00000        60 V  Q  H  L  Q  A  S  V  G  E  Q  G  E  G  S  T  I  L  Q  H  
ENSG00000       180 GTCCAGCACCTTCAGGCCTCGGTGGGAGAGCAGGGGGAGGGGAGCACCATCTTGCAACAC

ENSG00000        80 I  S  T  L  E  S  I  Y  Q  R  D  P  L  K  L  V  A  T  F  R  
ENSG00000       240 ATCAGCACCCTTGAGAGCATATATCAGAGGGACCCCCTGAAGCTGGTGGCCACTTTCAGA

ENSG00000       100 Q  I  L  Q  G  E  K  K  A  V  M  E  Q  F  R  H  L  P  M  P  
ENSG00000       300 CAAATACTTCAAGGAGAGAAAAAAGCTGTTATGGAACAGTTCCGCCACTTGCCAATGCCT

ENSG00000       120 F  H  W  K  Q  E  E  L  K  F  K  T  G  L  R  R  L  Q  H  R  
ENSG00000       360 TTCCACTGGAAGCAGGAAGAACTCAAGTTTAAGACAGGCTTGCGGAGGCTGCAGCACCGA

ENSG00000       140 V  G  E  I  H  L  L  R  E  A  L  Q  K  G  A  E  A  G  Q  V  
ENSG00000       420 GTAGGGGAGATCCACCTTCTCCGAGAAGCCCTGCAGAAGGGGGCTGAGGCTGGCCAAGTG

ENSG00000       160 S  L  H  S  L  I  E  T  P  A  N  G  T  G  P  S  E  A  L  A  
ENSG00000       480 TCTCTGCACAGCTTGATAGAAACTCCTGCTAATGGGACTGGGCCAAGTGAGGCCCTGGCC

ENSG00000       180 M  L  L  Q  E  T  T  G  E  L  E  A  A  K  A  L  V  L  K  R  
ENSG00000       540 ATGCTACTGCAGGAGACCACTGGAGAGCTAGAGGCAGCCAAAGCCCTAGTGCTGAAGAGG

ENSG00000       200 I  Q  I  W  K  R  Q  Q  Q  L  A  G  N  G  A  P  F  E  E  S  
ENSG00000       600 ATCCAGATTTGGAAACGGCAGCAGCAGCTGGCAGGGAATGGCGCACCGTTTGAGGAGAGC

ENSG00000       220 L  A  P  L  Q  E  R  C  E  S  L  V  D  I  Y  S  Q  L  Q  Q  
ENSG00000       660 CTGGCCCCACTCCAGGAGAGGTGTGAAAGCCTGGTGGACATTTATTCCCAGCTACAGCAG

ENSG00000       240 E  V  G  A  A  G  G  E  L  E  P  K  T  R  A  S  L  T  G  R  
ENSG00000       720 GAGGTAGGGGCGGCTGGTGGGGAGCTTGAGCCCAAGACCCGGGCATCGCTGACTGGCCGG

ENSG00000       260 L  D  E  V  L  R  T  L  V  T  S  C  F  L  V  E  K  Q  P  P  
ENSG00000       780 CTGGATGAAGTCCTGAGAACCCTCGTCACCAGTTGCTTCCTGGTGGAGAAGCAGCCCCCC

ENSG00000       280 Q  V  L  K  T  Q  T  K  F  Q  A  G  V  R  F  L  L  G  L  R  
ENSG00000       840 CAGGTACTGAAGACTCAGACCAAGTTCCAGGCTGGAGTTCGATTCCTGTTGGGCTTGAGG

ENSG00000       300 F  L  G  A  P  A  K  P  P  L  V  R  A  D  M  V  T  E  K  Q  
ENSG00000       900 TTCCTGGGGGCCCCAGCCAAGCCTCCGCTGGTCAGGGCCGACATGGTGACAGAGAAGCAG

ENSG00000       320 A  R  E  L  S  V  P  Q  G  P  G  A  G  A  E  S  T  G  E  I  
ENSG00000       960 GCGCGGGAGCTGAGTGTGCCTCAGGGTCCTGGGGCTGGAGCAGAAAGCACTGGAGAAATC

ENSG00000       340 I  N  N  T  V  P  L  E  N  S  I  P  G  N  C  C  S  A  L  F  
ENSG00000      1020 ATCAACAACACTGTGCCCTTGGAGAACAGCATTCCTGGGAACTGCTGCTCTGCCCTGTTC

ENSG00000       360 K  N  L  L  L  K  K  I  K  R  C  E  R  K  G  T  E  S  V  T  
ENSG00000      1080 AAGAACCTGCTTCTCAAGAAGATCAAGCGGTGTGAGCGGAAGGGCACTGAGTCTGTCACA

ENSG00000       380 E  E  K  C  A  V  L  F  S  A  S  F  T  L  G  P  G  K  L  P  
ENSG00000      1140 GAGGAGAAGTGCGCTGTGCTCTTCTCTGCCAGCTTCACACTTGGCCCCGGCAAACTCCCC

ENSG00000       400 I  Q  L  Q  A  L  S  L  P  L  V  V  I  V  H  G  N  Q  D  N  
ENSG00000      1200 ATCCAGCTCCAGGCCCTGTCTCTGCCCCTGGTGGTCATCGTCCATGGCAACCAAGACAAC

ENSG00000       420 N  A  K  A  T  I  L  W  D  N  A  F  S  E  M  D  R  V  P  F  
ENSG00000      1260 AATGCCAAAGCCACTATCCTGTGGGACAATGCCTTCTCTGAGATGGACCGCGTGCCCTTT

ENSG00000       440 V  V  A  E  R  V  P  W  E  K  M  C  E  T  L  N  L  K  F  M  
ENSG00000      1320 GTGGTGGCTGAGCGGGTGCCCTGGGAGAAGATGTGTGAAACTCTGAACCTGAAGTTCATG

ENSG00000       460 A  E  V  G  T  N  R  G  L  L  P  E  H  F  L  F  L  A  Q  K  
ENSG00000      1380 GCTGAGGTGGGGACCAACCGGGGGCTGCTCCCAGAGCACTTCCTCTTCCTGGCCCAGAAG

ENSG00000       480 I  F  N  D  N  S  L  S  M  E  A  F  Q  H  R  S  V  S  W  S  
ENSG00000      1440 ATCTTCAATGACAACAGCCTCAGTATGGAGGCCTTCCAGCACCGTTCTGTGTCCTGGTCG

ENSG00000       500 Q  F  N  K  E  I  L  L  G  R  G  F  T  F  W  Q  W  F  D  G  
ENSG00000      1500 CAGTTCAACAAGGAGATCCTGCTGGGCCGTGGCTTCACCTTTTGGCAGTGGTTTGATGGT

ENSG00000       520 V  L  D  L  T  K  R  C  L  R  S  Y  W  S  D  R  L  I  I  G  
ENSG00000      1560 GTCCTGGACCTCACCAAACGCTGTCTCCGGAGCTACTGGTCTGACCGGCTGATCATTGGC

ENSG00000       540 F  I  S  K  Q  Y  V  T  S  L  L  L  N  E  P  D  G  T  F  L  
ENSG00000      1620 TTCATCAGCAAACAGTACGTTACTAGCCTTCTTCTCAATGAGCCCGACGGAACCTTTCTC

ENSG00000       560 L  R  F  S  D  S  E  I  G  G  I  T  I  A  H  V  I  R  G  Q  
ENSG00000      1680 CTCCGCTTCAGCGACTCAGAGATTGGGGGCATCACCATTGCCCATGTCATCCGGGGCCAG

ENSG00000       580 D  G  S  P  Q  I  E  N  I  Q  P  F  S  A  K  D  L  S  I  R  
ENSG00000      1740 GATGGCTCTCCACAGATAGAGAACATCCAGCCATTCTCTGCCAAAGACCTGTCCATTCGC

ENSG00000       600 S  L  G  D  R  I  R  D  L  A  Q  L  K  N  L  Y  P  K  K  P  
ENSG00000      1800 TCACTGGGGGACCGAATCCGGGATCTTGCTCAGCTCAAAAATCTCTATCCCAAGAAGCCC

ENSG00000       620 K  D  E  A  F  R  S  H  Y  K  P  E  Q  M  G  K  D  G  R  G  
ENSG00000      1860 AAGGATGAGGCTTTCCGGAGCCACTACAAGCCTGAACAGATGGGTAAGGATGGCAGGGGT

ENSG00000       640 Y  V  P  A  T  I  K  M  T  V  E  R  D  Q  P  L  P  T  P  E  
ENSG00000      1920 TATGTCCCAGCTACCATCAAGATGACCGTGGAAAGGGACCAACCACTTCCTACCCCAGAG

ENSG00000       660 L  Q  M  P  T  M  V  P  S  Y  D  L  G  M  A  P  D  S  S  M  
ENSG00000      1980 CTCCAGATGCCTACCATGGTGCCTTCTTATGACCTTGGAATGGCCCCTGATTCCTCCATG

ENSG00000       680 S  M  Q  L  G  P  D  M  V  P  Q  V  Y  P  P  H  S  H  S  I  
ENSG00000      2040 AGCATGCAGCTTGGCCCAGATATGGTGCCCCAGGTGTACCCACCACACTCTCACTCCATC

ENSG00000       700 P  P  Y  Q  G  L  S  P  E  E  S  V  N  V  L  S  A  F  Q  E  
ENSG00000      2100 CCCCCGTATCAAGGCCTCTCCCCAGAAGAATCAGTCAACGTGTTGTCAGCCTTCCAGGAG

ENSG00000       720 P  H  L  Q  M  P  P  S  L  G  Q  M  S  L  P  F  D  Q  P  H  
ENSG00000      2160 CCTCACCTGCAGATGCCCCCCAGCCTGGGCCAGATGAGCCTGCCCTTTGACCAGCCTCAC

ENSG00000       740 P  Q  G  L  L  P  C  Q  P  Q  E  H  A  V  S  S  P  D  P  L  
ENSG00000      2220 CCCCAGGGCCTGCTGCCGTGCCAGCCTCAGGAGCATGCTGTGTCCAGCCCTGACCCCCTG

ENSG00000       760 L  C  S  D  V  T  M  V  E  D  S  C  L  S  Q  P  V  T  A  F  
ENSG00000      2280 CTCTGCTCAGATGTGACCATGGTGGAAGACAGCTGCCTGAGCCAGCCAGTGACAGCGTTT

ENSG00000       780 P  Q  G  T  W  I  G  E  D  I  F  P  P  L  L  P  P  T  E  Q  
ENSG00000      2340 CCTCAGGGCACTTGGATTGGTGAAGACATATTCCCTCCTCTGCTGCCTCCCACTGAACAG

ENSG00000       800 D  L  T  K  L  L  L  E  G  Q  G  E  S  G  G  G  S  L  G  A  
ENSG00000      2400 GACCTCACTAAGCTTCTCCTGGAGGGGCAAGGGGAGTCGGGGGGAGGGTCCTTGGGGGCA

ENSG00000       820 Q  P  L  L  Q  P  S  H  Y  G  Q  S  G  I  S  M  S  H  M  D  
ENSG00000      2460 CAGCCCCTCCTGCAGCCCTCCCACTATGGGCAATCTGGGATCTCAATGTCCCACATGGAC

ENSG00000       840 L  R  A  N  P  S  W    847
ENSG00000      2520 CTAAGGGCCAACCCCAGTTGG 2541
""",
        )
        protein_record = protein_alignment.sequences[4]
        nucleotide_record = nucleotide_records[protein_record.id]
        self.assertEqual(nucleotide_record.id, protein_record.id)
        alignments = aligner.align(protein_record, nucleotide_record)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 737], [0, 2211]]))
        )
        self.assertEqual(
            str(alignment),
            """\
ENSG00000         0 M  E  Q  F  R  H  L  P  M  P  F  H  W  K  Q  E  E  L  K  F  
ENSG00000         0 ATGGAACAGTTCCGCCACTTGCCAATGCCTTTCCACTGGAAGCAGGAAGAACTCAAGTTT

ENSG00000        20 K  T  G  L  R  R  L  Q  H  R  V  G  E  I  H  L  L  R  E  A  
ENSG00000        60 AAGACAGGCTTGCGGAGGCTGCAGCACCGAGTAGGGGAGATCCACCTTCTCCGAGAAGCC

ENSG00000        40 L  Q  K  G  A  E  A  G  Q  V  S  L  H  S  L  I  E  T  P  A  
ENSG00000       120 CTGCAGAAGGGGGCTGAGGCTGGCCAAGTGTCTCTGCACAGCTTGATAGAAACTCCTGCT

ENSG00000        60 N  G  T  G  P  S  E  A  L  A  M  L  L  Q  E  T  T  G  E  L  
ENSG00000       180 AATGGGACTGGGCCAAGTGAGGCCCTGGCCATGCTACTGCAGGAGACCACTGGAGAGCTA

ENSG00000        80 E  A  A  K  A  L  V  L  K  R  I  Q  I  W  K  R  Q  Q  Q  L  
ENSG00000       240 GAGGCAGCCAAAGCCCTAGTGCTGAAGAGGATCCAGATTTGGAAACGGCAGCAGCAGCTG

ENSG00000       100 A  G  N  G  A  P  F  E  E  S  L  A  P  L  Q  E  R  C  E  S  
ENSG00000       300 GCAGGGAATGGCGCACCGTTTGAGGAGAGCCTGGCCCCACTCCAGGAGAGGTGTGAAAGC

ENSG00000       120 L  V  D  I  Y  S  Q  L  Q  Q  E  V  G  A  A  G  G  E  L  E  
ENSG00000       360 CTGGTGGACATTTATTCCCAGCTACAGCAGGAGGTAGGGGCGGCTGGTGGGGAGCTTGAG

ENSG00000       140 P  K  T  R  A  S  L  T  G  R  L  D  E  V  L  R  T  L  V  T  
ENSG00000       420 CCCAAGACCCGGGCATCGCTGACTGGCCGGCTGGATGAAGTCCTGAGAACCCTCGTCACC

ENSG00000       160 S  C  F  L  V  E  K  Q  P  P  Q  V  L  K  T  Q  T  K  F  Q  
ENSG00000       480 AGTTGCTTCCTGGTGGAGAAGCAGCCCCCCCAGGTACTGAAGACTCAGACCAAGTTCCAG

ENSG00000       180 A  G  V  R  F  L  L  G  L  R  F  L  G  A  P  A  K  P  P  L  
ENSG00000       540 GCTGGAGTTCGATTCCTGTTGGGCTTGAGGTTCCTGGGGGCCCCAGCCAAGCCTCCGCTG

ENSG00000       200 V  R  A  D  M  V  T  E  K  Q  A  R  E  L  S  V  P  Q  G  P  
ENSG00000       600 GTCAGGGCCGACATGGTGACAGAGAAGCAGGCGCGGGAGCTGAGTGTGCCTCAGGGTCCT

ENSG00000       220 G  A  G  A  E  S  T  G  E  I  I  N  N  T  V  P  L  E  N  S  
ENSG00000       660 GGGGCTGGAGCAGAAAGCACTGGAGAAATCATCAACAACACTGTGCCCTTGGAGAACAGC

ENSG00000       240 I  P  G  N  C  C  S  A  L  F  K  N  L  L  L  K  K  I  K  R  
ENSG00000       720 ATTCCTGGGAACTGCTGCTCTGCCCTGTTCAAGAACCTGCTTCTCAAGAAGATCAAGCGG

ENSG00000       260 C  E  R  K  G  T  E  S  V  T  E  E  K  C  A  V  L  F  S  A  
ENSG00000       780 TGTGAGCGGAAGGGCACTGAGTCTGTCACAGAGGAGAAGTGCGCTGTGCTCTTCTCTGCC

ENSG00000       280 S  F  T  L  G  P  G  K  L  P  I  Q  L  Q  A  L  S  L  P  L  
ENSG00000       840 AGCTTCACACTTGGCCCCGGCAAACTCCCCATCCAGCTCCAGGCCCTGTCTCTGCCCCTG

ENSG00000       300 V  V  I  V  H  G  N  Q  D  N  N  A  K  A  T  I  L  W  D  N  
ENSG00000       900 GTGGTCATCGTCCATGGCAACCAAGACAACAATGCCAAAGCCACTATCCTGTGGGACAAT

ENSG00000       320 A  F  S  E  M  D  R  V  P  F  V  V  A  E  R  V  P  W  E  K  
ENSG00000       960 GCCTTCTCTGAGATGGACCGCGTGCCCTTTGTGGTGGCTGAGCGGGTGCCCTGGGAGAAG

ENSG00000       340 M  C  E  T  L  N  L  K  F  M  A  E  V  G  T  N  R  G  L  L  
ENSG00000      1020 ATGTGTGAAACTCTGAACCTGAAGTTCATGGCTGAGGTGGGGACCAACCGGGGGCTGCTC

ENSG00000       360 P  E  H  F  L  F  L  A  Q  K  I  F  N  D  N  S  L  S  M  E  
ENSG00000      1080 CCAGAGCACTTCCTCTTCCTGGCCCAGAAGATCTTCAATGACAACAGCCTCAGTATGGAG

ENSG00000       380 A  F  Q  H  R  S  V  S  W  S  Q  F  N  K  E  I  L  L  G  R  
ENSG00000      1140 GCCTTCCAGCACCGTTCTGTGTCCTGGTCGCAGTTCAACAAGGAGATCCTGCTGGGCCGT

ENSG00000       400 G  F  T  F  W  Q  W  F  D  G  V  L  D  L  T  K  R  C  L  R  
ENSG00000      1200 GGCTTCACCTTTTGGCAGTGGTTTGATGGTGTCCTGGACCTCACCAAACGCTGTCTCCGG

ENSG00000       420 S  Y  W  S  D  R  L  I  I  G  F  I  S  K  Q  Y  V  T  S  L  
ENSG00000      1260 AGCTACTGGTCTGACCGGCTGATCATTGGCTTCATCAGCAAACAGTACGTTACTAGCCTT

ENSG00000       440 L  L  N  E  P  D  G  T  F  L  L  R  F  S  D  S  E  I  G  G  
ENSG00000      1320 CTTCTCAATGAGCCCGACGGAACCTTTCTCCTCCGCTTCAGCGACTCAGAGATTGGGGGC

ENSG00000       460 I  T  I  A  H  V  I  R  G  Q  D  G  S  P  Q  I  E  N  I  Q  
ENSG00000      1380 ATCACCATTGCCCATGTCATCCGGGGCCAGGATGGCTCTCCACAGATAGAGAACATCCAG

ENSG00000       480 P  F  S  A  K  D  L  S  I  R  S  L  G  D  R  I  R  D  L  A  
ENSG00000      1440 CCATTCTCTGCCAAAGACCTGTCCATTCGCTCACTGGGGGACCGAATCCGGGATCTTGCT

ENSG00000       500 Q  L  K  N  L  Y  P  K  K  P  K  D  E  A  F  R  S  H  Y  K  
ENSG00000      1500 CAGCTCAAAAATCTCTATCCCAAGAAGCCCAAGGATGAGGCTTTCCGGAGCCACTACAAG

ENSG00000       520 P  E  Q  M  G  K  D  G  R  G  Y  V  P  A  T  I  K  M  T  V  
ENSG00000      1560 CCTGAACAGATGGGTAAGGATGGCAGGGGTTATGTCCCAGCTACCATCAAGATGACCGTG

ENSG00000       540 E  R  D  Q  P  L  P  T  P  E  L  Q  M  P  T  M  V  P  S  Y  
ENSG00000      1620 GAAAGGGACCAACCACTTCCTACCCCAGAGCTCCAGATGCCTACCATGGTGCCTTCTTAT

ENSG00000       560 D  L  G  M  A  P  D  S  S  M  S  M  Q  L  G  P  D  M  V  P  
ENSG00000      1680 GACCTTGGAATGGCCCCTGATTCCTCCATGAGCATGCAGCTTGGCCCAGATATGGTGCCC

ENSG00000       580 Q  V  Y  P  P  H  S  H  S  I  P  P  Y  Q  G  L  S  P  E  E  
ENSG00000      1740 CAGGTGTACCCACCACACTCTCACTCCATCCCCCCGTATCAAGGCCTCTCCCCAGAAGAA

ENSG00000       600 S  V  N  V  L  S  A  F  Q  E  P  H  L  Q  M  P  P  S  L  G  
ENSG00000      1800 TCAGTCAACGTGTTGTCAGCCTTCCAGGAGCCTCACCTGCAGATGCCCCCCAGCCTGGGC

ENSG00000       620 Q  M  S  L  P  F  D  Q  P  H  P  Q  G  L  L  P  C  Q  P  Q  
ENSG00000      1860 CAGATGAGCCTGCCCTTTGACCAGCCTCACCCCCAGGGCCTGCTGCCGTGCCAGCCTCAG

ENSG00000       640 E  H  A  V  S  S  P  D  P  L  L  C  S  D  V  T  M  V  E  D  
ENSG00000      1920 GAGCATGCTGTGTCCAGCCCTGACCCCCTGCTCTGCTCAGATGTGACCATGGTGGAAGAC

ENSG00000       660 S  C  L  S  Q  P  V  T  A  F  P  Q  G  T  W  I  G  E  D  I  
ENSG00000      1980 AGCTGCCTGAGCCAGCCAGTGACAGCGTTTCCTCAGGGCACTTGGATTGGTGAAGACATA

ENSG00000       680 F  P  P  L  L  P  P  T  E  Q  D  L  T  K  L  L  L  E  G  Q  
ENSG00000      2040 TTCCCTCCTCTGCTGCCTCCCACTGAACAGGACCTCACTAAGCTTCTCCTGGAGGGGCAA

ENSG00000       700 G  E  S  G  G  G  S  L  G  A  Q  P  L  L  Q  P  S  H  Y  G  
ENSG00000      2100 GGGGAGTCGGGGGGAGGGTCCTTGGGGGCACAGCCCCTCCTGCAGCCCTCCCACTATGGG

ENSG00000       720 Q  S  G  I  S  M  S  H  M  D  L  R  A  N  P  S  W    737
ENSG00000      2160 CAATCTGGGATCTCAATGTCCCACATGGACCTAAGGGCCAACCCCAGTTGG 2211
""",
        )
        protein_record = protein_alignment.sequences[5]
        nucleotide_record = nucleotide_records[protein_record.id]
        self.assertEqual(nucleotide_record.id, protein_record.id)
        alignments = aligner.align(protein_record, nucleotide_record)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 737], [0, 2211]]))
        )
        self.assertEqual(
            str(alignment),
            """\
ENSG00000         0 M  E  Q  F  R  H  L  P  M  P  F  H  W  K  Q  E  E  L  K  F  
ENSG00000         0 ATGGAACAGTTCCGCCACTTGCCAATGCCTTTCCACTGGAAGCAGGAAGAACTCAAGTTT

ENSG00000        20 K  T  G  L  R  R  L  Q  H  R  V  G  E  I  H  L  L  R  E  A  
ENSG00000        60 AAGACAGGCTTGCGGAGGCTGCAGCACCGAGTAGGGGAGATCCACCTTCTCCGAGAAGCC

ENSG00000        40 L  Q  K  G  A  E  A  G  Q  V  S  L  H  S  L  I  E  T  P  A  
ENSG00000       120 CTGCAGAAGGGGGCTGAGGCTGGCCAAGTGTCTCTGCACAGCTTGATAGAAACTCCTGCT

ENSG00000        60 N  G  T  G  P  S  E  A  L  A  M  L  L  Q  E  T  T  G  E  L  
ENSG00000       180 AATGGGACTGGGCCAAGTGAGGCCCTGGCCATGCTACTGCAGGAGACCACTGGAGAGCTA

ENSG00000        80 E  A  A  K  A  L  V  L  K  R  I  Q  I  W  K  R  Q  Q  Q  L  
ENSG00000       240 GAGGCAGCCAAAGCCCTAGTGCTGAAGAGGATCCAGATTTGGAAACGGCAGCAGCAGCTG

ENSG00000       100 A  G  N  G  A  P  F  E  E  S  L  A  P  L  Q  E  R  C  E  S  
ENSG00000       300 GCAGGGAATGGCGCACCGTTTGAGGAGAGCCTGGCCCCACTCCAGGAGAGGTGTGAAAGC

ENSG00000       120 L  V  D  I  Y  S  Q  L  Q  Q  E  V  G  A  A  G  G  E  L  E  
ENSG00000       360 CTGGTGGACATTTATTCCCAGCTACAGCAGGAGGTAGGGGCGGCTGGTGGGGAGCTTGAG

ENSG00000       140 P  K  T  R  A  S  L  T  G  R  L  D  E  V  L  R  T  L  V  T  
ENSG00000       420 CCCAAGACCCGGGCATCGCTGACTGGCCGGCTGGATGAAGTCCTGAGAACCCTCGTCACC

ENSG00000       160 S  C  F  L  V  E  K  Q  P  P  Q  V  L  K  T  Q  T  K  F  Q  
ENSG00000       480 AGTTGCTTCCTGGTGGAGAAGCAGCCCCCCCAGGTACTGAAGACTCAGACCAAGTTCCAG

ENSG00000       180 A  G  V  R  F  L  L  G  L  R  F  L  G  A  P  A  K  P  P  L  
ENSG00000       540 GCTGGAGTTCGATTCCTGTTGGGCTTGAGGTTCCTGGGGGCCCCAGCCAAGCCTCCGCTG

ENSG00000       200 V  R  A  D  M  V  T  E  K  Q  A  R  E  L  S  V  P  Q  G  P  
ENSG00000       600 GTCAGGGCCGACATGGTGACAGAGAAGCAGGCGCGGGAGCTGAGTGTGCCTCAGGGTCCT

ENSG00000       220 G  A  G  A  E  S  T  G  E  I  I  N  N  T  V  P  L  E  N  S  
ENSG00000       660 GGGGCTGGAGCAGAAAGCACTGGAGAAATCATCAACAACACTGTGCCCTTGGAGAACAGC

ENSG00000       240 I  P  G  N  C  C  S  A  L  F  K  N  L  L  L  K  K  I  K  R  
ENSG00000       720 ATTCCTGGGAACTGCTGCTCTGCCCTGTTCAAGAACCTGCTTCTCAAGAAGATCAAGCGG

ENSG00000       260 C  E  R  K  G  T  E  S  V  T  E  E  K  C  A  V  L  F  S  A  
ENSG00000       780 TGTGAGCGGAAGGGCACTGAGTCTGTCACAGAGGAGAAGTGCGCTGTGCTCTTCTCTGCC

ENSG00000       280 S  F  T  L  G  P  G  K  L  P  I  Q  L  Q  A  L  S  L  P  L  
ENSG00000       840 AGCTTCACACTTGGCCCCGGCAAACTCCCCATCCAGCTCCAGGCCCTGTCTCTGCCCCTG

ENSG00000       300 V  V  I  V  H  G  N  Q  D  N  N  A  K  A  T  I  L  W  D  N  
ENSG00000       900 GTGGTCATCGTCCATGGCAACCAAGACAACAATGCCAAAGCCACTATCCTGTGGGACAAT

ENSG00000       320 A  F  S  E  M  D  R  V  P  F  V  V  A  E  R  V  P  W  E  K  
ENSG00000       960 GCCTTCTCTGAGATGGACCGCGTGCCCTTTGTGGTGGCTGAGCGGGTGCCCTGGGAGAAG

ENSG00000       340 M  C  E  T  L  N  L  K  F  M  A  E  V  G  T  N  R  G  L  L  
ENSG00000      1020 ATGTGTGAAACTCTGAACCTGAAGTTCATGGCTGAGGTGGGGACCAACCGGGGGCTGCTC

ENSG00000       360 P  E  H  F  L  F  L  A  Q  K  I  F  N  D  N  S  L  S  M  E  
ENSG00000      1080 CCAGAGCACTTCCTCTTCCTGGCCCAGAAGATCTTCAATGACAACAGCCTCAGTATGGAG

ENSG00000       380 A  F  Q  H  R  S  V  S  W  S  Q  F  N  K  E  I  L  L  G  R  
ENSG00000      1140 GCCTTCCAGCACCGTTCTGTGTCCTGGTCGCAGTTCAACAAGGAGATCCTGCTGGGCCGT

ENSG00000       400 G  F  T  F  W  Q  W  F  D  G  V  L  D  L  T  K  R  C  L  R  
ENSG00000      1200 GGCTTCACCTTTTGGCAGTGGTTTGATGGTGTCCTGGACCTCACCAAACGCTGTCTCCGG

ENSG00000       420 S  Y  W  S  D  R  L  I  I  G  F  I  S  K  Q  Y  V  T  S  L  
ENSG00000      1260 AGCTACTGGTCTGACCGGCTGATCATTGGCTTCATCAGCAAACAGTACGTTACTAGCCTT

ENSG00000       440 L  L  N  E  P  D  G  T  F  L  L  R  F  S  D  S  E  I  G  G  
ENSG00000      1320 CTTCTCAATGAGCCCGACGGAACCTTTCTCCTCCGCTTCAGCGACTCAGAGATTGGGGGC

ENSG00000       460 I  T  I  A  H  V  I  R  G  Q  D  G  S  P  Q  I  E  N  I  Q  
ENSG00000      1380 ATCACCATTGCCCATGTCATCCGGGGCCAGGATGGCTCTCCACAGATAGAGAACATCCAG

ENSG00000       480 P  F  S  A  K  D  L  S  I  R  S  L  G  D  R  I  R  D  L  A  
ENSG00000      1440 CCATTCTCTGCCAAAGACCTGTCCATTCGCTCACTGGGGGACCGAATCCGGGATCTTGCT

ENSG00000       500 Q  L  K  N  L  Y  P  K  K  P  K  D  E  A  F  R  S  H  Y  K  
ENSG00000      1500 CAGCTCAAAAATCTCTATCCCAAGAAGCCCAAGGATGAGGCTTTCCGGAGCCACTACAAG

ENSG00000       520 P  E  Q  M  G  K  D  G  R  G  Y  V  P  A  T  I  K  M  T  V  
ENSG00000      1560 CCTGAACAGATGGGTAAGGATGGCAGGGGTTATGTCCCAGCTACCATCAAGATGACCGTG

ENSG00000       540 E  R  D  Q  P  L  P  T  P  E  L  Q  M  P  T  M  V  P  S  Y  
ENSG00000      1620 GAAAGGGACCAACCACTTCCTACCCCAGAGCTCCAGATGCCTACCATGGTGCCTTCTTAT

ENSG00000       560 D  L  G  M  A  P  D  S  S  M  S  M  Q  L  G  P  D  M  V  P  
ENSG00000      1680 GACCTTGGAATGGCCCCTGATTCCTCCATGAGCATGCAGCTTGGCCCAGATATGGTGCCC

ENSG00000       580 Q  V  Y  P  P  H  S  H  S  I  P  P  Y  Q  G  L  S  P  E  E  
ENSG00000      1740 CAGGTGTACCCACCACACTCTCACTCCATCCCCCCGTATCAAGGCCTCTCCCCAGAAGAA

ENSG00000       600 S  V  N  V  L  S  A  F  Q  E  P  H  L  Q  M  P  P  S  L  G  
ENSG00000      1800 TCAGTCAACGTGTTGTCAGCCTTCCAGGAGCCTCACCTGCAGATGCCCCCCAGCCTGGGC

ENSG00000       620 Q  M  S  L  P  F  D  Q  P  H  P  Q  G  L  L  P  C  Q  P  Q  
ENSG00000      1860 CAGATGAGCCTGCCCTTTGACCAGCCTCACCCCCAGGGCCTGCTGCCGTGCCAGCCTCAG

ENSG00000       640 E  H  A  V  S  S  P  D  P  L  L  C  S  D  V  T  M  V  E  D  
ENSG00000      1920 GAGCATGCTGTGTCCAGCCCTGACCCCCTGCTCTGCTCAGATGTGACCATGGTGGAAGAC

ENSG00000       660 S  C  L  S  Q  P  V  T  A  F  P  Q  G  T  W  I  G  E  D  I  
ENSG00000      1980 AGCTGCCTGAGCCAGCCAGTGACAGCGTTTCCTCAGGGCACTTGGATTGGTGAAGACATA

ENSG00000       680 F  P  P  L  L  P  P  T  E  Q  D  L  T  K  L  L  L  E  G  Q  
ENSG00000      2040 TTCCCTCCTCTGCTGCCTCCCACTGAACAGGACCTCACTAAGCTTCTCCTGGAGGGGCAA

ENSG00000       700 G  E  S  G  G  G  S  L  G  A  Q  P  L  L  Q  P  S  H  Y  G  
ENSG00000      2100 GGGGAGTCGGGGGGAGGGTCCTTGGGGGCACAGCCCCTCCTGCAGCCCTCCCACTATGGG

ENSG00000       720 Q  S  G  I  S  M  S  H  M  D  L  R  A  N  P  S  W    737
ENSG00000      2160 CAATCTGGGATCTCAATGTCCCACATGGACCTAAGGGCCAACCCCAGTTGG 2211
""",
        )
        protein_record = protein_alignment.sequences[6]
        nucleotide_record = nucleotide_records[protein_record.id]
        self.assertEqual(nucleotide_record.id, protein_record.id)
        alignments = aligner.align(protein_record, nucleotide_record)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 737], [0, 2211]]))
        )
        self.assertEqual(
            str(alignment),
            """\
ENSG00000         0 M  E  Q  F  R  H  L  P  M  P  F  H  W  K  Q  E  E  L  K  F  
ENSG00000         0 ATGGAACAGTTCCGCCACTTGCCAATGCCTTTCCACTGGAAGCAGGAAGAACTCAAGTTT

ENSG00000        20 K  T  G  L  R  R  L  Q  H  R  V  G  E  I  H  L  L  R  E  A  
ENSG00000        60 AAGACAGGCTTGCGGAGGCTGCAGCACCGAGTAGGGGAGATCCACCTTCTCCGAGAAGCC

ENSG00000        40 L  Q  K  G  A  E  A  G  Q  V  S  L  H  S  L  I  E  T  P  A  
ENSG00000       120 CTGCAGAAGGGGGCTGAGGCTGGCCAAGTGTCTCTGCACAGCTTGATAGAAACTCCTGCT

ENSG00000        60 N  G  T  G  P  S  E  A  L  A  M  L  L  Q  E  T  T  G  E  L  
ENSG00000       180 AATGGGACTGGGCCAAGTGAGGCCCTGGCCATGCTACTGCAGGAGACCACTGGAGAGCTA

ENSG00000        80 E  A  A  K  A  L  V  L  K  R  I  Q  I  W  K  R  Q  Q  Q  L  
ENSG00000       240 GAGGCAGCCAAAGCCCTAGTGCTGAAGAGGATCCAGATTTGGAAACGGCAGCAGCAGCTG

ENSG00000       100 A  G  N  G  A  P  F  E  E  S  L  A  P  L  Q  E  R  C  E  S  
ENSG00000       300 GCAGGGAATGGCGCACCGTTTGAGGAGAGCCTGGCCCCACTCCAGGAGAGGTGTGAAAGC

ENSG00000       120 L  V  D  I  Y  S  Q  L  Q  Q  E  V  G  A  A  G  G  E  L  E  
ENSG00000       360 CTGGTGGACATTTATTCCCAGCTACAGCAGGAGGTAGGGGCGGCTGGTGGGGAGCTTGAG

ENSG00000       140 P  K  T  R  A  S  L  T  G  R  L  D  E  V  L  R  T  L  V  T  
ENSG00000       420 CCCAAGACCCGGGCATCGCTGACTGGCCGGCTGGATGAAGTCCTGAGAACCCTCGTCACC

ENSG00000       160 S  C  F  L  V  E  K  Q  P  P  Q  V  L  K  T  Q  T  K  F  Q  
ENSG00000       480 AGTTGCTTCCTGGTGGAGAAGCAGCCCCCCCAGGTACTGAAGACTCAGACCAAGTTCCAG

ENSG00000       180 A  G  V  R  F  L  L  G  L  R  F  L  G  A  P  A  K  P  P  L  
ENSG00000       540 GCTGGAGTTCGATTCCTGTTGGGCTTGAGGTTCCTGGGGGCCCCAGCCAAGCCTCCGCTG

ENSG00000       200 V  R  A  D  M  V  T  E  K  Q  A  R  E  L  S  V  P  Q  G  P  
ENSG00000       600 GTCAGGGCCGACATGGTGACAGAGAAGCAGGCGCGGGAGCTGAGTGTGCCTCAGGGTCCT

ENSG00000       220 G  A  G  A  E  S  T  G  E  I  I  N  N  T  V  P  L  E  N  S  
ENSG00000       660 GGGGCTGGAGCAGAAAGCACTGGAGAAATCATCAACAACACTGTGCCCTTGGAGAACAGC

ENSG00000       240 I  P  G  N  C  C  S  A  L  F  K  N  L  L  L  K  K  I  K  R  
ENSG00000       720 ATTCCTGGGAACTGCTGCTCTGCCCTGTTCAAGAACCTGCTTCTCAAGAAGATCAAGCGG

ENSG00000       260 C  E  R  K  G  T  E  S  V  T  E  E  K  C  A  V  L  F  S  A  
ENSG00000       780 TGTGAGCGGAAGGGCACTGAGTCTGTCACAGAGGAGAAGTGCGCTGTGCTCTTCTCTGCC

ENSG00000       280 S  F  T  L  G  P  G  K  L  P  I  Q  L  Q  A  L  S  L  P  L  
ENSG00000       840 AGCTTCACACTTGGCCCCGGCAAACTCCCCATCCAGCTCCAGGCCCTGTCTCTGCCCCTG

ENSG00000       300 V  V  I  V  H  G  N  Q  D  N  N  A  K  A  T  I  L  W  D  N  
ENSG00000       900 GTGGTCATCGTCCATGGCAACCAAGACAACAATGCCAAAGCCACTATCCTGTGGGACAAT

ENSG00000       320 A  F  S  E  M  D  R  V  P  F  V  V  A  E  R  V  P  W  E  K  
ENSG00000       960 GCCTTCTCTGAGATGGACCGCGTGCCCTTTGTGGTGGCTGAGCGGGTGCCCTGGGAGAAG

ENSG00000       340 M  C  E  T  L  N  L  K  F  M  A  E  V  G  T  N  R  G  L  L  
ENSG00000      1020 ATGTGTGAAACTCTGAACCTGAAGTTCATGGCTGAGGTGGGGACCAACCGGGGGCTGCTC

ENSG00000       360 P  E  H  F  L  F  L  A  Q  K  I  F  N  D  N  S  L  S  M  E  
ENSG00000      1080 CCAGAGCACTTCCTCTTCCTGGCCCAGAAGATCTTCAATGACAACAGCCTCAGTATGGAG

ENSG00000       380 A  F  Q  H  R  S  V  S  W  S  Q  F  N  K  E  I  L  L  G  R  
ENSG00000      1140 GCCTTCCAGCACCGTTCTGTGTCCTGGTCGCAGTTCAACAAGGAGATCCTGCTGGGCCGT

ENSG00000       400 G  F  T  F  W  Q  W  F  D  G  V  L  D  L  T  K  R  C  L  R  
ENSG00000      1200 GGCTTCACCTTTTGGCAGTGGTTTGATGGTGTCCTGGACCTCACCAAACGCTGTCTCCGG

ENSG00000       420 S  Y  W  S  D  R  L  I  I  G  F  I  S  K  Q  Y  V  T  S  L  
ENSG00000      1260 AGCTACTGGTCTGACCGGCTGATCATTGGCTTCATCAGCAAACAGTACGTTACTAGCCTT

ENSG00000       440 L  L  N  E  P  D  G  T  F  L  L  R  F  S  D  S  E  I  G  G  
ENSG00000      1320 CTTCTCAATGAGCCCGACGGAACCTTTCTCCTCCGCTTCAGCGACTCAGAGATTGGGGGC

ENSG00000       460 I  T  I  A  H  V  I  R  G  Q  D  G  S  P  Q  I  E  N  I  Q  
ENSG00000      1380 ATCACCATTGCCCATGTCATCCGGGGCCAGGATGGCTCTCCACAGATAGAGAACATCCAG

ENSG00000       480 P  F  S  A  K  D  L  S  I  R  S  L  G  D  R  I  R  D  L  A  
ENSG00000      1440 CCATTCTCTGCCAAAGACCTGTCCATTCGCTCACTGGGGGACCGAATCCGGGATCTTGCT

ENSG00000       500 Q  L  K  N  L  Y  P  K  K  P  K  D  E  A  F  R  S  H  Y  K  
ENSG00000      1500 CAGCTCAAAAATCTCTATCCCAAGAAGCCCAAGGATGAGGCTTTCCGGAGCCACTACAAG

ENSG00000       520 P  E  Q  M  G  K  D  G  R  G  Y  V  P  A  T  I  K  M  T  V  
ENSG00000      1560 CCTGAACAGATGGGTAAGGATGGCAGGGGTTATGTCCCAGCTACCATCAAGATGACCGTG

ENSG00000       540 E  R  D  Q  P  L  P  T  P  E  L  Q  M  P  T  M  V  P  S  Y  
ENSG00000      1620 GAAAGGGACCAACCACTTCCTACCCCAGAGCTCCAGATGCCTACCATGGTGCCTTCTTAT

ENSG00000       560 D  L  G  M  A  P  D  S  S  M  S  M  Q  L  G  P  D  M  V  P  
ENSG00000      1680 GACCTTGGAATGGCCCCTGATTCCTCCATGAGCATGCAGCTTGGCCCAGATATGGTGCCC

ENSG00000       580 Q  V  Y  P  P  H  S  H  S  I  P  P  Y  Q  G  L  S  P  E  E  
ENSG00000      1740 CAGGTGTACCCACCACACTCTCACTCCATCCCCCCGTATCAAGGCCTCTCCCCAGAAGAA

ENSG00000       600 S  V  N  V  L  S  A  F  Q  E  P  H  L  Q  M  P  P  S  L  G  
ENSG00000      1800 TCAGTCAACGTGTTGTCAGCCTTCCAGGAGCCTCACCTGCAGATGCCCCCCAGCCTGGGC

ENSG00000       620 Q  M  S  L  P  F  D  Q  P  H  P  Q  G  L  L  P  C  Q  P  Q  
ENSG00000      1860 CAGATGAGCCTGCCCTTTGACCAGCCTCACCCCCAGGGCCTGCTGCCGTGCCAGCCTCAG

ENSG00000       640 E  H  A  V  S  S  P  D  P  L  L  C  S  D  V  T  M  V  E  D  
ENSG00000      1920 GAGCATGCTGTGTCCAGCCCTGACCCCCTGCTCTGCTCAGATGTGACCATGGTGGAAGAC

ENSG00000       660 S  C  L  S  Q  P  V  T  A  F  P  Q  G  T  W  I  G  E  D  I  
ENSG00000      1980 AGCTGCCTGAGCCAGCCAGTGACAGCGTTTCCTCAGGGCACTTGGATTGGTGAAGACATA

ENSG00000       680 F  P  P  L  L  P  P  T  E  Q  D  L  T  K  L  L  L  E  G  Q  
ENSG00000      2040 TTCCCTCCTCTGCTGCCTCCCACTGAACAGGACCTCACTAAGCTTCTCCTGGAGGGGCAA

ENSG00000       700 G  E  S  G  G  G  S  L  G  A  Q  P  L  L  Q  P  S  H  Y  G  
ENSG00000      2100 GGGGAGTCGGGGGGAGGGTCCTTGGGGGCACAGCCCCTCCTGCAGCCCTCCCACTATGGG

ENSG00000       720 Q  S  G  I  S  M  S  H  M  D  L  R  A  N  P  S  W    737
ENSG00000      2160 CAATCTGGGATCTCAATGTCCCACATGGACCTAAGGGCCAACCCCAGTTGG 2211
""",
        )
        protein_record = protein_alignment.sequences[7]
        nucleotide_record = nucleotide_records[protein_record.id]
        self.assertEqual(nucleotide_record.id, protein_record.id)
        nucleotide_record = nucleotide_record.upper()
        alignments = aligner.align(protein_record, nucleotide_record)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 1021], [0, 3063]]))
        )
        self.assertEqual(
            str(alignment),
            """\
isotig466         0 E  V  T  Q  S  R  R  K  P  V  R  E  G  R  P  W  E  P  S  Q  
isotig466         0 GAGGTTACTCAGAGTAGGAGGAAGCCGGTCAGAGAGGGCAGACCCTGGGAACCTTCGCAG

isotig466        20 E  L  A  Q  T  P  E  G  G  C  P  R  Q  E  S  R  P  E  G  G  
isotig466        60 GAGTTGGCCCAAACCCCTGAGGGAGGCTGTCCTCGGCAGGAGAGCAGACCTGAGGGAGGT

isotig466        40 W  C  A  G  P  G  C  P  A  L  P  S  H  P  R  A  G  S  P  E  
isotig466       120 TGGTGTGCCGGTCCTGGTTGTCCAGCCCTCCCGTCCCATCCCAGGGCTGGCTCACCAGAA

isotig466        60 K  S  T  G  I  P  G  P  L  G  R  E  Q  R  V  S  S  C  P  G  
isotig466       180 AAAAGTACAGGCATCCCAGGCCCTCTGGGGCGGGAGCAGAGGGTCTCCTCTTGTCCGGGG

isotig466        80 E  D  K  N  K  N  K  K  R  V  C  S  P  H  W  K  S  S  G  A  
isotig466       240 GAAGACAAAAACAAAAACAAAAAACGAGTGTGTAGCCCTCACTGGAAGTCTTCTGGTGCT

isotig466       100 L  G  P  F  A  L  G  S  L  A  S  G  H  G  G  R  A  P  P  G  
isotig466       300 CTGGGGCCGTTTGCACTTGGGAGCCTGGCTTCTGGGCATGGTGGCCGGGCTCCGCCGGGG

isotig466       120 A  S  A  L  G  K  A  F  L  E  Q  S  G  W  E  P  T  E  V  P  
isotig466       360 GCTTCAGCCCTCGGCAAAGCGTTTCTAGAACAGAGTGGGTGGGAGCCGACTGAAGTCCCA

isotig466       140 E  P  R  T  L  T  H  R  K  A  S  S  G  W  R  N  A  N  P  P  
isotig466       420 GAACCGCGAACACTGACGCACAGGAAAGCCTCGAGTGGGTGGAGAAATGCAAATCCCCCA

isotig466       160 V  A  W  P  Q  P  L  R  T  L  T  A  D  S  D  F  G  L  E  R  
isotig466       480 GTGGCATGGCCTCAGCCGCTGCGGACCCTGACGGCCGATTCTGACTTCGGACTTGAGCGG

isotig466       180 D  Q  H  P  D  P  G  R  A  S  M  S  Q  W  N  Q  V  Q  Q  L  
isotig466       540 GACCAGCACCCGGACCCAGGGCGAGCCAGCATGTCTCAGTGGAATCAAGTCCAACAGTTA

isotig466       200 E  I  K  F  L  E  Q  V  D  Q  F  Y  D  D  N  F  P  M  E  I  
isotig466       600 GAAATCAAGTTTTTGGAGCAAGTTGATCAGTTCTATGATGACAACTTTCCCATGGAAATC

isotig466       220 R  H  L  L  A  Q  W  I  E  H  Q  D  W  E  V  A  S  N  N  E  
isotig466       660 CGACATCTGCTGGCCCAGTGGATTGAGCATCAAGACTGGGAGGTGGCCTCTAACAATGAA

isotig466       240 T  M  A  T  I  L  L  Q  N  L  L  I  Q  L  D  E  Q  L  G  R  
isotig466       720 ACTATGGCAACAATTCTTCTTCAAAACTTATTAATACAATTGGATGAACAGTTAGGTCGT

isotig466       260 V  S  K  E  K  N  L  L  L  I  H  N  L  K  R  I  R  K  V  L  
isotig466       780 GTTTCCAAAGAGAAAAACCTGCTATTGATCCACAATCTAAAGAGAATTAGAAAAGTACTT

isotig466       280 Q  G  K  F  H  G  N  P  M  H  V  A  V  V  I  S  N  C  L  R  
isotig466       840 CAGGGGAAGTTTCATGGAAATCCAATGCATGTAGCCGTGGTAATCTCAAATTGTTTAAGG

isotig466       300 E  E  R  R  I  L  A  A  A  N  M  P  I  Q  G  P  L  E  K  S  
isotig466       900 GAAGAGAGGAGAATACTGGCTGCAGCGAACATGCCTATCCAGGGACCTCTGGAGAAATCC

isotig466       320 L  Q  S  S  S  V  S  E  R  Q  R  N  V  E  H  K  V  A  A  I  
isotig466       960 TTACAAAGTTCGTCGGTTTCAGAAAGACAGAGAAATGTGGAACACAAAGTGGCTGCCATT

isotig466       340 K  N  S  V  Q  M  T  E  Q  D  T  K  Y  L  E  D  L  Q  D  E  
isotig466      1020 AAAAACAGTGTGCAGATGACAGAACAAGACACCAAATACTTGGAAGATCTGCAAGATGAA

isotig466       360 F  D  Y  R  Y  K  T  I  Q  T  M  D  Q  G  D  K  N  S  I  L  
isotig466      1080 TTTGACTACAGGTATAAAACAATTCAGACAATGGACCAGGGTGACAAGAATAGCATCCTA

isotig466       380 M  N  Q  E  V  L  T  L  Q  E  M  L  N  S  L  D  F  K  R  K  
isotig466      1140 ATGAACCAGGAGGTTTTGACACTCCAAGAAATGCTTAATAGCCTGGACTTCAAGAGAAAG

isotig466       400 E  A  L  T  K  M  T  Q  I  V  N  E  S  D  L  L  M  S  S  M  
isotig466      1200 GAAGCACTCACTAAGATGACACAGATAGTGAACGAGTCGGACCTGCTGATGAGCAGCATG

isotig466       420 L  I  E  E  L  Q  D  W  K  R  R  Q  Q  I  A  C  I  G  G  P  
isotig466      1260 CTCATAGAAGAGCTGCAGGACTGGAAGAGGAGGCAGCAGATCGCCTGCATCGGTGGCCCA

isotig466       440 L  H  N  G  L  D  Q  L  Q  N  C  F  T  L  L  A  E  S  L  F  
isotig466      1320 CTCCACAACGGGCTGGACCAGCTTCAGAACTGCTTTACCCTGTTGGCAGAAAGTCTTTTC

isotig466       460 Q  L  R  R  Q  L  E  K  L  E  E  Q  S  S  K  M  T  Y  E  G  
isotig466      1380 CAACTCAGACGACAGCTGGAGAAATTAGAGGAGCAGTCTTCCAAGATGACTTACGAAGGA

isotig466       480 D  P  I  P  T  Q  R  A  H  L  L  E  R  A  T  F  L  I  Y  N  
isotig466      1440 GACCCCATCCCCACGCAGAGAGCACACCTGCTGGAGAGAGCCACCTTCCTGATCTACAAC

isotig466       500 L  F  K  N  S  F  V  V  E  R  Q  P  C  M  P  T  H  P  Q  R  
isotig466      1500 CTTTTCAAGAACTCATTTGTGGTTGAGCGACAGCCCTGCATGCCAACACACCCTCAGAGG

isotig466       520 P  L  V  L  K  T  L  I  Q  F  T  A  K  L  R  L  L  I  K  L  
isotig466      1560 CCGCTGGTACTCAAAACCCTCATTCAGTTCACCGCGAAACTGAGACTACTAATAAAATTG

isotig466       540 P  E  L  N  Y  Q  V  K  V  K  A  S  I  D  K  N  V  S  T  L  
isotig466      1620 CCGGAACTCAACTATCAGGTGAAAGTAAAGGCATCGATCGACAAGAATGTTTCAACGCTA

isotig466       560 S  N  R  R  F  V  L  C  G  T  Q  V  K  A  M  S  I  E  E  S  
isotig466      1680 AGCAATAGAAGATTTGTGCTTTGTGGAACTCAAGTCAAAGCCATGTCCATCGAGGAATCC

isotig466       580 S  N  G  S  L  S  V  E  F  R  H  L  Q  P  K  E  M  K  S  S  
isotig466      1740 TCCAATGGGAGCCTCTCAGTAGAATTTAGACATTTGCAACCGAAGGAAATGAAATCCAGT

isotig466       600 A  G  S  K  G  N  E  G  C  H  M  V  T  E  E  L  H  S  I  A  
isotig466      1800 GCCGGAAGTAAAGGAAATGAGGGCTGCCACATGGTGACGGAAGAGCTGCATTCCATAGCC

isotig466       620 F  E  T  Q  I  C  L  Y  G  L  T  I  D  L  E  T  S  S  L  P  
isotig466      1860 TTTGAGACCCAGATCTGCCTCTATGGCCTCACCATCGACTTGGAGACAAGCTCATTACCT

isotig466       640 V  V  M  I  S  N  V  S  Q  L  P  N  A  W  A  S  I  I  W  Y  
isotig466      1920 GTGGTGATGATTTCTAATGTCAGCCAACTGCCTAATGCTTGGGCATCCATCATTTGGTAC

isotig466       660 N  V  S  T  N  D  C  Q  N  L  V  F  F  N  N  P  P  P  V  T  
isotig466      1980 AATGTGTCAACCAACGATTGCCAGAACTTGGTTTTCTTTAATAATCCTCCGCCTGTCACT

isotig466       680 L  S  Q  L  L  E  V  M  S  W  Q  F  S  S  Y  V  G  R  G  L  
isotig466      2040 TTGAGTCAACTCCTGGAAGTGATGAGCTGGCAGTTTTCATCCTATGTTGGTCGTGGCCTT

isotig466       700 N  S  D  Q  L  N  M  L  A  E  K  L  T  V  Q  S  N  Y  S  D  
isotig466      2100 AATTCAGACCAGCTCAACATGCTGGCAGAGAAGCTCACAGTTCAGTCTAACTACAGCGAT

isotig466       720 G  H  L  T  W  A  K  F  C  K  E  H  L  P  G  K  P  F  T  F  
isotig466      2160 GGTCACCTCACCTGGGCCAAGTTCTGCAAGGAACACTTGCCTGGCAAACCATTTACCTTC

isotig466       740 W  T  W  L  E  A  I  L  D  L  I  K  K  H  I  L  P  L  W  I  
isotig466      2220 TGGACCTGGCTTGAAGCAATATTGGACCTAATTAAAAAACACATTCTTCCCCTCTGGATT

isotig466       760 D  G  Y  I  M  G  F  V  S  K  E  K  E  R  F  L  L  K  D  K  
isotig466      2280 GATGGGTACATCATGGGCTTCGTGAGCAAAGAGAAGGAGAGGTTTCTGCTCAAGGATAAA

isotig466       780 M  P  G  T  F  L  L  R  F  S  E  S  H  L  G  G  I  T  F  T  
isotig466      2340 ATGCCCGGGACATTTTTGTTACGATTCAGTGAGAGCCATCTCGGAGGGATCACCTTCACC

isotig466       800 W  V  D  H  S  E  N  G  E  V  R  F  H  S  V  E  P  Y  N  K  
isotig466      2400 TGGGTGGACCACTCTGAAAACGGAGAAGTGAGATTCCACTCCGTAGAACCCTACAACAAA

isotig466       820 G  R  L  S  A  L  P  F  A  D  I  L  R  D  Y  K  V  I  M  A  
isotig466      2460 GGGCGTCTGTCGGCCCTGCCATTTGCTGACATCCTGCGGGACTACAAGGTCATCATGGCT

isotig466       840 E  N  I  P  E  N  P  L  K  Y  L  Y  P  D  I  P  K  D  K  A  
isotig466      2520 GAGAACATTCCCGAGAACCCTCTCAAGTACCTCTACCCCGACATCCCCAAAGACAAAGCC

isotig466       860 F  G  K  H  Y  S  S  Q  P  C  E  V  S  R  P  T  E  R  G  D  
isotig466      2580 TTCGGTAAACACTACAGCTCCCAGCCTTGCGAAGTTTCAAGGCCAACAGAACGGGGAGAC

isotig466       880 K  G  Y  V  P  S  V  F  I  P  I  S  T  I  R  S  D  A  M  E  
isotig466      2640 AAAGGTTATGTTCCTTCAGTTTTTATCCCTATTTCAACAATCCGCAGCGACGCCATGGAG

isotig466       900 P  Q  S  P  S  D  L  L  P  M  S  P  S  V  Y  A  V  L  R  E  
isotig466      2700 CCGCAGTCTCCTTCAGACCTTCTCCCCATGTCTCCGAGTGTATACGCTGTGCTGAGAGAA

isotig466       920 N  L  S  P  T  T  I  E  T  A  M  K  S  P  Y  S  E  R  Y  K  
isotig466      2760 AACCTGAGCCCTACCACAATTGAAACAGCAATGAAGTCTCCATATTCTGAGCGGTACAAA

isotig466       940 A  T  L  Q  G  R  E  Q  M  K  T  E  T  A  L  C  Q  S  P  Q  
isotig466      2820 GCGACTCTTCAAGGAAGAGAGCAGATGAAAACGGAGACTGCTCTTTGCCAAAGTCCACAA

isotig466       960 F  I  S  S  A  L  I  L  V  S  R  K  W  H  K  S  E  A  F  L  
isotig466      2880 TTCATTTCTTCAGCTTTGATACTGGTTTCTAGAAAATGGCACAAATCCGAAGCTTTCCTC

isotig466       980 S  L  G  D  I  P  Q  L  G  V  L  L  K  C  K  P  K  L  Q  I  
isotig466      2940 TCACTAGGTGACATTCCCCAACTGGGAGTGCTGCTGAAATGCAAACCAAAGCTTCAGATA

isotig466      1000 N  T  Q  E  K  T  A  S  R  N  L  C  S  Q  Y  N  R  R  L  L  
isotig466      3000 AACACGCAGGAAAAGACAGCTTCGAGAAACCTATGTTCGCAATATAACAGAAGGCTGCTT

isotig466      1020 C   1021
isotig466      3060 TGC 3063
""",
        )
        protein_record = protein_alignment.sequences[8]
        nucleotide_record = nucleotide_records[protein_record.id]
        self.assertEqual(nucleotide_record.id, protein_record.id)
        nucleotide_record = nucleotide_record.upper()
        alignments = aligner.align(protein_record, nucleotide_record)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 806], [0, 2418]]))
        )
        self.assertEqual(
            str(alignment),
            """\
isotig125         0 A  R  R  G  Q  A  A  L  G  S  P  A  A  R  T  W  S  Q  R  S  
isotig125         0 GCTAGGAGAGGCCAGGCGGCCCTCGGGAGCCCAGCTGCTCGCACCTGGAGCCAGCGCAGC

isotig125        20 P  A  S  R  A  S  A  R  E  T  V  T  P  P  D  C  G  R  M  A  
isotig125        60 CCGGCCAGTCGGGCCTCAGCCCGGGAGACAGTTACGCCCCCTGATTGCGGCAGGATGGCC

isotig125        40 Q  W  N  Q  L  Q  Q  L  D  T  R  Y  L  E  Q  L  H  Q  L  Y  
isotig125       120 CAGTGGAACCAGCTGCAGCAGCTGGACACTCGGTACCTGGAGCAGCTGCACCAGCTGTAC

isotig125        60 S  D  S  F  P  M  E  L  R  Q  F  L  A  P  W  I  E  S  Q  D  
isotig125       180 AGCGACAGCTTCCCCATGGAGCTGCGGCAGTTCCTGGCACCTTGGATTGAGAGTCAAGAT

isotig125        80 W  A  Y  A  A  S  K  E  S  H  A  T  L  V  F  H  N  L  L  G  
isotig125       240 TGGGCATATGCAGCCAGCAAAGAGTCACATGCCACACTGGTGTTTCATAATCTCTTGGGT

isotig125       100 E  I  D  Q  Q  Y  S  R  F  L  Q  E  S  N  V  L  Y  Q  H  N  
isotig125       300 GAGATTGACCAGCAGTACAGCCGATTCCTGCAAGAGTCCAACGTCCTCTATCAGCACAAC

isotig125       120 L  R  R  I  K  Q  F  L  Q  S  R  Y  L  E  K  P  M  E  I  A  
isotig125       360 CTTCGGAGGATCAAGCAGTTCCTACAGAGCAGGTATCTTGAGAAGCCGATGGAAATCGCC

isotig125       140 R  I  V  A  R  C  L  W  E  E  S  R  L  L  Q  T  A  A  T  A  
isotig125       420 CGCATCGTGGCCCGATGCCTGTGGGAAGAGTCTCGCCTCCTCCAGACGGCAGCCACTGCA

isotig125       160 A  Q  Q  G  G  Q  A  N  H  P  T  A  A  V  V  T  E  K  Q  Q  
isotig125       480 GCCCAGCAAGGAGGCCAGGCCAACCACCCAACAGCTGCTGTGGTGACGGAGAAACAGCAG

isotig125       180 M  L  E  Q  H  L  Q  D  V  R  K  R  V  Q  D  L  E  Q  K  M  
isotig125       540 ATGCTGGAGCAGCATCTTCAGGATGTCCGGAAACGTGTGCAGGATCTAGAACAGAAAATG

isotig125       200 K  V  V  E  N  L  Q  D  D  F  D  F  N  Y  K  T  L  K  S  Q  
isotig125       600 AAAGTGGTAGAGAATCTCCAGGATGACTTTGATTTCAACTATAAAACCCTCAAGAGTCAA

isotig125       220 G  D  M  Q  D  L  N  G  N  N  Q  S  V  T  R  Q  K  M  Q  Q  
isotig125       660 GGAGACATGCAGGATCTGAATGGAAACAACCAGTCTGTGACCAGGCAGAAGATGCAGCAG

isotig125       240 L  E  Q  M  L  T  A  L  D  Q  M  R  R  S  I  V  S  E  L  A  
isotig125       720 CTGGAACAGATGCTCACGGCGCTGGACCAGATGCGGAGAAGCATTGTGAGTGAGCTGGCG

isotig125       260 G  L  L  S  A  M  E  Y  V  Q  K  T  L  T  D  E  E  L  A  D  
isotig125       780 GGGCTTTTGTCGGCAATGGAGTACGTGCAGAAAACACTCACAGACGAGGAGCTGGCTGAC

isotig125       280 W  K  R  R  Q  Q  I  A  C  I  G  G  P  P  N  I  C  L  D  R  
isotig125       840 TGGAAGAGGCGGCAGCAGATCGCGTGCATTGGAGGCCCTCCCAACATCTGCCTGGATCGC

isotig125       300 L  E  N  W  I  T  S  L  A  E  S  Q  L  Q  T  R  Q  Q  I  K  
isotig125       900 CTGGAAAACTGGATAACTTCGTTAGCAGAATCTCAACTTCAGACCCGCCAACAAATTAAG

isotig125       320 K  L  E  E  L  Q  Q  K  V  S  Y  K  G  D  P  I  V  Q  H  R  
isotig125       960 AAACTGGAGGAGCTACAGCAGAAGGTGTCCTACAAGGGGGACCCCATTGTGCAGCACCGG

isotig125       340 P  M  L  E  E  R  I  V  E  L  F  R  N  L  M  K  S  A  F  V  
isotig125      1020 CCGATGCTGGAGGAGCGGATCGTGGAGCTGTTCAGAAACTTGATGAAGAGTGCCTTCGTG

isotig125       360 V  E  R  Q  P  C  M  P  M  H  P  D  R  P  L  V  I  K  T  G  
isotig125      1080 GTGGAGCGACAGCCCTGCATGCCGATGCACCCCGACCGGCCCTTGGTCATCAAGACTGGT

isotig125       380 V  Q  F  T  T  K  V  R  L  L  V  K  F  P  E  L  N  Y  Q  L  
isotig125      1140 GTCCAGTTCACTACTAAAGTCAGGTTGTTGGTCAAGTTTCCCGAGTTGAATTATCAGCTT

isotig125       400 K  I  K  V  C  I  D  K  D  S  G  D  V  A  A  L  R  G  S  R  
isotig125      1200 AAAATTAAAGTGTGCATTGACAAAGATTCTGGGGACGTTGCTGCTCTCAGAGGATCTCGG

isotig125       420 K  F  N  I  L  G  T  N  T  K  V  M  N  M  E  E  S  N  N  G  
isotig125      1260 AAATTTAACATTCTGGGCACAAACACGAAGGTGATGAACATGGAAGAATCCAACAACGGC

isotig125       440 S  L  S  A  E  F  K  H  L  T  L  R  E  Q  R  C  G  N  G  G  
isotig125      1320 AGCCTGTCTGCGGAGTTCAAGCACTTGACCCTGAGGGAGCAGAGATGTGGGAATGGAGGC

isotig125       460 R  A  N  C  D  A  S  L  I  V  T  E  E  L  H  L  I  T  F  E  
isotig125      1380 CGTGCCAATTGTGATGCCTCCTTGATTGTGACCGAGGAGCTGCATCTGATCACCTTCGAG

isotig125       480 T  E  V  Y  H  Q  G  L  K  I  D  L  E  T  H  S  L  P  V  V  
isotig125      1440 ACTGAGGTGTACCACCAAGGCCTCAAGATTGACCTGGAGACCCATTCTTTGCCAGTTGTG

isotig125       500 V  I  S  N  I  C  Q  M  P  N  A  W  A  S  I  L  W  Y  N  M  
isotig125      1500 GTGATCTCCAACATCTGTCAGATGCCAAATGCCTGGGCATCCATCCTGTGGTATAACATG

isotig125       520 L  T  N  N  P  K  N  V  N  F  F  T  K  P  P  I  G  T  W  D  
isotig125      1560 CTGACCAACAACCCCAAGAACGTGAACTTCTTCACCAAGCCACCAATCGGAACCTGGGAC

isotig125       540 Q  V  A  E  V  L  S  W  Q  F  S  S  T  T  K  R  G  L  S  I  
isotig125      1620 CAGGTGGCCGAGGTGCTCAGCTGGCAGTTCTCATCCACCACAAAGCGAGGGCTGAGCATC

isotig125       560 E  Q  L  T  T  L  A  E  K  L  L  G  P  G  V  N  Y  S  G  C  
isotig125      1680 GAGCAGCTGACTACGCTGGCCGAGAAGCTCCTAGGACCTGGTGTCAACTACTCCGGGTGT

isotig125       580 Q  I  T  W  A  K  F  C  K  E  N  M  A  G  K  G  F  S  F  W  
isotig125      1740 CAGATCACATGGGCTAAATTTTGCAAAGAAAACATGGCTGGCAAGGGCTTCTCCTTCTGG

isotig125       600 V  W  L  D  N  I  I  D  L  V  K  K  Y  I  L  A  L  W  N  E  
isotig125      1800 GTGTGGCTAGACAATATCATTGACCTTGTGAAAAAGTATATCTTGGCCCTCTGGAATGAA

isotig125       620 G  Y  I  M  G  F  I  S  K  E  R  E  R  A  I  L  S  T  K  P  
isotig125      1860 GGGTACATCATGGGCTTCATTAGCAAGGAGCGGGAGCGGGCGATCCTGAGCACGAAACCC

isotig125       640 P  G  T  F  L  L  R  F  S  E  S  S  K  E  G  G  V  T  F  T  
isotig125      1920 CCGGGCACCTTCCTGCTGAGATTCAGCGAGAGCAGCAAAGAAGGAGGGGTCACTTTCACT

isotig125       660 W  V  E  K  D  I  S  G  K  T  Q  I  Q  S  V  E  P  Y  T  K  
isotig125      1980 TGGGTGGAAAAGGACATCAGTGGCAAGACCCAGATCCAGTCTGTAGAGCCGTACACCAAG

isotig125       680 Q  Q  L  N  N  M  S  F  A  E  I  I  M  G  Y  K  I  M  D  A  
isotig125      2040 CAGCAGCTGAACAACATGTCCTTTGCTGAAATCATCATGGGCTACAAGATCATGGATGCC

isotig125       700 T  N  I  L  V  S  P  L  V  Y  L  Y  P  D  I  P  K  E  E  A  
isotig125      2100 ACCAACATCCTGGTGTCCCCATTGGTCTACCTCTACCCTGACATTCCCAAGGAGGAGGCG

isotig125       720 F  G  K  Y  C  R  P  E  S  Q  E  H  P  E  A  D  P  G  S  A  
isotig125      2160 TTCGGGAAGTACTGTCGACCAGAGAGCCAGGAGCATCCTGAAGCTGACCCCGGTAGTGCC

isotig125       740 A  P  Y  L  K  T  K  F  I  C  V  T  P  T  T  C  S  N  T  I  
isotig125      2220 GCCCCTTACCTGAAGACCAAGTTCATCTGTGTGACACCAACGACCTGCAGCAATACCATT

isotig125       760 D  L  P  M  S  P  P  H  F  R  F  I  D  A  V  W  K  R  R  C  
isotig125      2280 GACCTGCCGATGTCCCCCCCGCACTTTAGATTCATTGATGCAGTTTGGAAACGGAGGTGC

isotig125       780 A  L  G  R  R  A  V  V  T  H  V  H  G  S  D  F  G  V  R  Y  
isotig125      2340 GCCCTCGGCAGGAGGGCAGTTGTCACTCACGTTCATGGATCTGACTTCGGAGTGCGCTAC

isotig125       800 L  P  H  V  R  S    806
isotig125      2400 CTCCCCCATGTGAGGAGC 2418
""",
        )
        protein_record = protein_alignment.sequences[9]
        nucleotide_record = nucleotide_records[protein_record.id]
        self.assertEqual(nucleotide_record.id, protein_record.id)
        nucleotide_record = nucleotide_record.upper()
        alignments = aligner.align(protein_record, nucleotide_record)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 796], [0, 2388]]))
        )
        self.assertEqual(
            str(alignment),
            """\
isotig125         0 A  R  R  G  Q  A  A  L  G  S  P  A  A  R  T  W  S  Q  R  S  
isotig125         0 GCTAGGAGAGGCCAGGCGGCCCTCGGGAGCCCAGCTGCTCGCACCTGGAGCCAGCGCAGC

isotig125        20 P  A  S  R  A  S  A  R  E  T  V  T  P  P  D  C  G  R  M  A  
isotig125        60 CCGGCCAGTCGGGCCTCAGCCCGGGAGACAGTTACGCCCCCTGATTGCGGCAGGATGGCC

isotig125        40 Q  W  N  Q  L  Q  Q  L  D  T  R  Y  L  E  Q  L  H  Q  L  Y  
isotig125       120 CAGTGGAACCAGCTGCAGCAGCTGGACACTCGGTACCTGGAGCAGCTGCACCAGCTGTAC

isotig125        60 S  D  S  F  P  M  E  L  R  Q  F  L  A  P  W  I  E  S  Q  D  
isotig125       180 AGCGACAGCTTCCCCATGGAGCTGCGGCAGTTCCTGGCACCTTGGATTGAGAGTCAAGAT

isotig125        80 W  A  Y  A  A  S  K  E  S  H  A  T  L  V  F  H  N  L  L  G  
isotig125       240 TGGGCATATGCAGCCAGCAAAGAGTCACATGCCACACTGGTGTTTCATAATCTCTTGGGT

isotig125       100 E  I  D  Q  Q  Y  S  R  F  L  Q  E  S  N  V  L  Y  Q  H  N  
isotig125       300 GAGATTGACCAGCAGTACAGCCGATTCCTGCAAGAGTCCAACGTCCTCTATCAGCACAAC

isotig125       120 L  R  R  I  K  Q  F  L  Q  S  R  Y  L  E  K  P  M  E  I  A  
isotig125       360 CTTCGGAGGATCAAGCAGTTCCTACAGAGCAGGTATCTTGAGAAGCCGATGGAAATCGCC

isotig125       140 R  I  V  A  R  C  L  W  E  E  S  R  L  L  Q  T  A  A  T  A  
isotig125       420 CGCATCGTGGCCCGATGCCTGTGGGAAGAGTCTCGCCTCCTCCAGACGGCAGCCACTGCA

isotig125       160 A  Q  Q  G  G  Q  A  N  H  P  T  A  A  V  V  T  E  K  Q  Q  
isotig125       480 GCCCAGCAAGGAGGCCAGGCCAACCACCCAACAGCTGCTGTGGTGACGGAGAAACAGCAG

isotig125       180 M  L  E  Q  H  L  Q  D  V  R  K  R  V  Q  D  L  E  Q  K  M  
isotig125       540 ATGCTGGAGCAGCATCTTCAGGATGTCCGGAAACGTGTGCAGGATCTAGAACAGAAAATG

isotig125       200 K  V  V  E  N  L  Q  D  D  F  D  F  N  Y  K  T  L  K  S  Q  
isotig125       600 AAAGTGGTAGAGAATCTCCAGGATGACTTTGATTTCAACTATAAAACCCTCAAGAGTCAA

isotig125       220 G  D  M  Q  D  L  N  G  N  N  Q  S  V  T  R  Q  K  M  Q  Q  
isotig125       660 GGAGACATGCAGGATCTGAATGGAAACAACCAGTCTGTGACCAGGCAGAAGATGCAGCAG

isotig125       240 L  E  Q  M  L  T  A  L  D  Q  M  R  R  S  I  V  S  E  L  A  
isotig125       720 CTGGAACAGATGCTCACGGCGCTGGACCAGATGCGGAGAAGCATTGTGAGTGAGCTGGCG

isotig125       260 G  L  L  S  A  M  E  Y  V  Q  K  T  L  T  D  E  E  L  A  D  
isotig125       780 GGGCTTTTGTCGGCAATGGAGTACGTGCAGAAAACACTCACAGACGAGGAGCTGGCTGAC

isotig125       280 W  K  R  R  Q  Q  I  A  C  I  G  G  P  P  N  I  C  L  D  R  
isotig125       840 TGGAAGAGGCGGCAGCAGATCGCGTGCATTGGAGGCCCTCCCAACATCTGCCTGGATCGC

isotig125       300 L  E  N  W  I  T  S  L  A  E  S  Q  L  Q  T  R  Q  Q  I  K  
isotig125       900 CTGGAAAACTGGATAACTTCGTTAGCAGAATCTCAACTTCAGACCCGCCAACAAATTAAG

isotig125       320 K  L  E  E  L  Q  Q  K  V  S  Y  K  G  D  P  I  V  Q  H  R  
isotig125       960 AAACTGGAGGAGCTACAGCAGAAGGTGTCCTACAAGGGGGACCCCATTGTGCAGCACCGG

isotig125       340 P  M  L  E  E  R  I  V  E  L  F  R  N  L  M  K  S  A  F  V  
isotig125      1020 CCGATGCTGGAGGAGCGGATCGTGGAGCTGTTCAGAAACTTGATGAAGAGTGCCTTCGTG

isotig125       360 V  E  R  Q  P  C  M  P  M  H  P  D  R  P  L  V  I  K  T  G  
isotig125      1080 GTGGAGCGACAGCCCTGCATGCCGATGCACCCCGACCGGCCCTTGGTCATCAAGACTGGT

isotig125       380 V  Q  F  T  T  K  V  R  L  L  V  K  F  P  E  L  N  Y  Q  L  
isotig125      1140 GTCCAGTTCACTACTAAAGTCAGGTTGTTGGTCAAGTTTCCCGAGTTGAATTATCAGCTT

isotig125       400 K  I  K  V  C  I  D  K  D  S  G  D  V  A  A  L  R  G  S  R  
isotig125      1200 AAAATTAAAGTGTGCATTGACAAAGATTCTGGGGACGTTGCTGCTCTCAGAGGATCTCGG

isotig125       420 K  F  N  I  L  G  T  N  T  K  V  M  N  M  E  E  S  N  N  G  
isotig125      1260 AAATTTAACATTCTGGGCACAAACACGAAGGTGATGAACATGGAAGAATCCAACAACGGC

isotig125       440 S  L  S  A  E  F  K  H  L  T  L  R  E  Q  R  C  G  N  G  G  
isotig125      1320 AGCCTGTCTGCGGAGTTCAAGCACTTGACCCTGAGGGAGCAGAGATGTGGGAATGGAGGC

isotig125       460 R  A  N  C  D  A  S  L  I  V  T  E  E  L  H  L  I  T  F  E  
isotig125      1380 CGTGCCAATTGTGATGCCTCCTTGATTGTGACCGAGGAGCTGCATCTGATCACCTTCGAG

isotig125       480 T  E  V  Y  H  Q  G  L  K  I  D  L  E  T  H  S  L  P  V  V  
isotig125      1440 ACTGAGGTGTACCACCAAGGCCTCAAGATTGACCTGGAGACCCATTCTTTGCCAGTTGTG

isotig125       500 V  I  S  N  I  C  Q  M  P  N  A  W  A  S  I  L  W  Y  N  M  
isotig125      1500 GTGATCTCCAACATCTGTCAGATGCCAAATGCCTGGGCATCCATCCTGTGGTATAACATG

isotig125       520 L  T  N  N  P  K  N  V  N  F  F  T  K  P  P  I  G  T  W  D  
isotig125      1560 CTGACCAACAACCCCAAGAACGTGAACTTCTTCACCAAGCCACCAATCGGAACCTGGGAC

isotig125       540 Q  V  A  E  V  L  S  W  Q  F  S  S  T  T  K  R  G  L  S  I  
isotig125      1620 CAGGTGGCCGAGGTGCTCAGCTGGCAGTTCTCATCCACCACAAAGCGAGGGCTGAGCATC

isotig125       560 E  Q  L  T  T  L  A  E  K  L  L  G  P  G  V  N  Y  S  G  C  
isotig125      1680 GAGCAGCTGACTACGCTGGCCGAGAAGCTCCTAGGACCTGGTGTCAACTACTCCGGGTGT

isotig125       580 Q  I  T  W  A  K  F  C  K  E  N  M  A  G  K  G  F  S  F  W  
isotig125      1740 CAGATCACATGGGCTAAATTTTGCAAAGAAAACATGGCTGGCAAGGGCTTCTCCTTCTGG

isotig125       600 V  W  L  D  N  I  I  D  L  V  K  K  Y  I  L  A  L  W  N  E  
isotig125      1800 GTGTGGCTAGACAATATCATTGACCTTGTGAAAAAGTATATCTTGGCCCTCTGGAATGAA

isotig125       620 G  Y  I  M  G  F  I  S  K  E  R  E  R  A  I  L  S  T  K  P  
isotig125      1860 GGGTACATCATGGGCTTCATTAGCAAGGAGCGGGAGCGGGCGATCCTGAGCACGAAACCC

isotig125       640 P  G  T  F  L  L  R  F  S  E  S  S  K  E  G  G  V  T  F  T  
isotig125      1920 CCGGGCACCTTCCTGCTGAGATTCAGCGAGAGCAGCAAAGAAGGAGGGGTCACTTTCACT

isotig125       660 W  V  E  K  D  I  S  G  K  T  Q  I  Q  S  V  E  P  Y  T  K  
isotig125      1980 TGGGTGGAAAAGGACATCAGTGGCAAGACCCAGATCCAGTCTGTAGAGCCGTACACCAAG

isotig125       680 Q  Q  L  N  N  M  S  F  A  E  I  I  M  G  Y  K  I  M  D  A  
isotig125      2040 CAGCAGCTGAACAACATGTCCTTTGCTGAAATCATCATGGGCTACAAGATCATGGATGCC

isotig125       700 T  N  I  L  V  S  P  L  V  Y  L  Y  P  D  I  P  K  E  E  A  
isotig125      2100 ACCAACATCCTGGTGTCCCCATTGGTCTACCTCTACCCTGACATTCCCAAGGAGGAGGCG

isotig125       720 F  G  K  Y  C  R  P  E  S  Q  E  H  P  E  A  D  P  G  S  C  
isotig125      2160 TTCGGGAAGTACTGTCGACCAGAGAGCCAGGAGCATCCTGAAGCTGACCCCGGTAGTTGT

isotig125       740 F  S  M  V  L  V  S  L  L  G  K  G  G  Q  C  R  S  L  E  E  
isotig125      2220 TTTTCCATGGTTCTGGTTTCGCTGTTAGGGAAAGGGGGACAGTGCAGGTCCTTGGAGGAG

isotig125       760 R  Q  G  H  D  R  V  S  G  G  E  S  C  Y  G  R  A  V  Y  W  
isotig125      2280 AGACAAGGACATGACCGGGTGTCTGGTGGTGAGTCCTGCTATGGAAGAGCTGTTTATTGG

isotig125       780 V  L  Q  G  D  R  D  S  R  E  D  Q  N  Q  A  S    796
isotig125      2340 GTACTTCAGGGTGACCGGGATTCAAGAGAAGACCAGAATCAGGCCTCA 2388
""",
        )
        alignment = protein_alignment.mapall(codon_alignments)
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[   0,    0,    0,    0,    0,    0,    0,    0,
                              0,    0,    0,    0,   36,   36,   63,   63,
                            213,  222,  240,  240,  330,  330,  330,  405,
                            405,  438,  450,  579,  579,  729,  729,  837,
                            837,  897,  912,  975,  987, 1020, 1026, 1074,
                           1077, 1107, 1107, 1116, 1119, 1182, 1188, 1302,
                           1302, 1392, 1398, 1455, 1458, 1470, 1695, 1695,
                           1824, 1824, 1887, 1887, 1896, 1902, 1956, 1980,
                           2022, 2025, 2043, 2046, 2061, 2106, 2187, 2220,
                           2232, 2253, 2271, 2331, 2334, 2352, 2364, 2373,
                           2382, 2424, 2463, 2484, 2517, 2541],
                          [   0,    0,    0,    0,    0,    0,    0,    0,
                              0,    0,    0,    0,   36,   36,   63,   63,
                            213,  222,  240,  240,  330,  330,  330,  405,
                            405,  438,  450,  579,  579,  729,  729,  837,
                            837,  897,  912,  975,  987, 1020, 1026, 1074,
                           1077, 1107, 1107, 1116, 1119, 1182, 1188, 1302,
                           1302, 1392, 1398, 1455, 1458, 1470, 1695, 1695,
                           1824, 1824, 1887, 1887, 1896, 1902, 1956, 1980,
                           2022, 2025, 2043, 2046, 2061, 2106, 2187, 2220,
                           2232, 2253, 2271, 2331, 2334, 2352, 2364, 2373,
                           2382, 2424, 2463, 2484, 2517, 2541],
                          [   0,    0,    0,    0,    0,    0,    0,    0,
                              0,    0,    0,    0,   36,   36,   63,   63,
                            213,  222,  240,  240,  330,  330,  330,  405,
                            405,  438,  450,  579,  579,  729,  729,  837,
                            837,  897,  912,  975,  987, 1020, 1026, 1074,
                           1077, 1107, 1107, 1116, 1119, 1182, 1188, 1302,
                           1302, 1392, 1398, 1455, 1458, 1470, 1695, 1695,
                           1824, 1824, 1887, 1887, 1896, 1902, 1956, 1980,
                           2022, 2025, 2043, 2046, 2061, 2106, 2187, 2220,
                           2232, 2253, 2271, 2331, 2334, 2352, 2364, 2373,
                           2382, 2424, 2463, 2484, 2517, 2541],
                          [   0,    0,    0,    0,    0,    0,    0,    0,
                              0,    0,    0,    0,   36,   36,   63,   63,
                            213,  222,  240,  240,  330,  330,  330,  405,
                            405,  438,  450,  579,  579,  729,  729,  837,
                            837,  897,  912,  975,  987, 1020, 1026, 1074,
                           1077, 1107, 1107, 1116, 1119, 1182, 1188, 1302,
                           1302, 1392, 1398, 1455, 1458, 1470, 1695, 1695,
                           1824, 1824, 1887, 1887, 1896, 1902, 1956, 1980,
                           2022, 2025, 2043, 2046, 2061, 2106, 2187, 2220,
                           2232, 2253, 2271, 2331, 2334, 2352, 2364, 2373,
                           2382, 2424, 2463, 2484, 2517, 2541],
                          [   0,    0,    0,    0,    0,    0,    0,    0,
                              0,    0,    0,    0,    0,    0,    0,    0,
                              0,    0,    0,    0,    0,    0,    0,   75,
                             75,  108,  120,  249,  249,  399,  399,  507,
                            507,  567,  582,  645,  657,  690,  696,  744,
                            747,  777,  777,  786,  789,  852,  858,  972,
                            972, 1062, 1068, 1125, 1128, 1140, 1365, 1365,
                           1494, 1494, 1557, 1557, 1566, 1572, 1626, 1650,
                           1692, 1695, 1713, 1716, 1731, 1776, 1857, 1890,
                           1902, 1923, 1941, 2001, 2004, 2022, 2034, 2043,
                           2052, 2094, 2133, 2154, 2187, 2211],
                          [   0,    0,    0,    0,    0,    0,    0,    0,
                              0,    0,    0,    0,    0,    0,    0,    0,
                              0,    0,    0,    0,    0,    0,    0,   75,
                             75,  108,  120,  249,  249,  399,  399,  507,
                            507,  567,  582,  645,  657,  690,  696,  744,
                            747,  777,  777,  786,  789,  852,  858,  972,
                            972, 1062, 1068, 1125, 1128, 1140, 1365, 1365,
                           1494, 1494, 1557, 1557, 1566, 1572, 1626, 1650,
                           1692, 1695, 1713, 1716, 1731, 1776, 1857, 1890,
                           1902, 1923, 1941, 2001, 2004, 2022, 2034, 2043,
                           2052, 2094, 2133, 2154, 2187, 2211],
                          [   0,    0,    0,    0,    0,    0,    0,    0,
                              0,    0,    0,    0,    0,    0,    0,    0,
                              0,    0,    0,    0,    0,    0,    0,   75,
                             75,  108,  120,  249,  249,  399,  399,  507,
                            507,  567,  582,  645,  657,  690,  696,  744,
                            747,  777,  777,  786,  789,  852,  858,  972,
                            972, 1062, 1068, 1125, 1128, 1140, 1365, 1365,
                           1494, 1494, 1557, 1557, 1566, 1572, 1626, 1650,
                           1692, 1695, 1713, 1716, 1731, 1776, 1857, 1890,
                           1902, 1923, 1941, 2001, 2004, 2022, 2034, 2043,
                           2052, 2094, 2133, 2154, 2187, 2211],
                          [   0,   12,   21,  357,  378,  417,  444,  474,
                            516,  549,  564,  570,  606,  615,  642,  645,
                            795,  795,  813,  828,  918,  918,  999, 1074,
                           1077, 1110, 1110, 1239, 1251, 1401, 1428, 1536,
                           1560, 1620, 1620, 1683, 1683, 1716, 1716, 1764,
                           1764, 1794, 1809, 1818, 1818, 1881, 1881, 1995,
                           2001, 2091, 2091, 2148, 2148, 2148, 2373, 2373,
                           2502, 2532, 2595, 2619, 2628, 2628, 2682, 2706,
                           2706, 2709, 2727, 2727, 2727, 2772, 2772, 2805,
                           2817, 2817, 2835, 2895, 2898, 2916, 2928, 2937,
                           2946, 2946, 2985, 3006, 3039, 3063],
                          [   0,    0,    9,    9,   30,   30,   57,   57,
                             99,   99,  114,  114,  150,  159,  186,  189,
                            339,  339,  357,  372,  462,  465,  546,  621,
                            624,  657,  669,  798,  810,  960,  987, 1095,
                           1119, 1179, 1179, 1242, 1254, 1287, 1287, 1335,
                           1335, 1365, 1380, 1389, 1392, 1455, 1455, 1569,
                           1575, 1665, 1665, 1722, 1725, 1725, 1950, 1953,
                           2082, 2112, 2175, 2175, 2184, 2184, 2184, 2208,
                           2208, 2208, 2226, 2226, 2226, 2271, 2271, 2304,
                           2304, 2304, 2304, 2304, 2307, 2325, 2337, 2346,
                           2346, 2346, 2385, 2385, 2418, 2418],
                          [   0,    0,    9,    9,   30,   30,   57,   57,
                             99,   99,  114,  114,  150,  159,  186,  189,
                            339,  339,  357,  372,  462,  465,  546,  621,
                            624,  657,  669,  798,  810,  960,  987, 1095,
                           1119, 1179, 1179, 1242, 1254, 1287, 1287, 1335,
                           1335, 1365, 1380, 1389, 1392, 1455, 1455, 1569,
                           1575, 1665, 1665, 1722, 1725, 1725, 1950, 1953,
                           2082, 2112, 2175, 2175, 2184, 2184, 2184, 2208,
                           2208, 2208, 2226, 2229, 2229, 2229, 2229, 2262,
                           2274, 2274, 2274, 2334, 2337, 2337, 2349, 2349,
                           2349, 2349, 2388, 2388, 2388, 2388]])
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "clustal"),
            """\
ENSG00000166888:ENST0000030013      --------------------------------------------------
ENSG00000166888:ENST0000054387      --------------------------------------------------
ENSG00000166888:ENST0000055615      --------------------------------------------------
ENSG00000166888:ENST0000045407      --------------------------------------------------
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         GAGGTTACTCAGAGTAGGAGGAAGCCGGTCAGAGAGGGCAGACCCTGGGA
isotig12565                         ------------GCTAGGAGA-----------------------------
isotig12566                         ------------GCTAGGAGA-----------------------------

ENSG00000166888:ENST0000030013      --------------------------------------------------
ENSG00000166888:ENST0000054387      --------------------------------------------------
ENSG00000166888:ENST0000055615      --------------------------------------------------
ENSG00000166888:ENST0000045407      --------------------------------------------------
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         ACCTTCGCAGGAGTTGGCCCAAACCCCTGAGGGAGGCTGTCCTCGGCAGG
isotig12565                         --------------------------------------------------
isotig12566                         --------------------------------------------------

ENSG00000166888:ENST0000030013      --------------------------------------------------
ENSG00000166888:ENST0000054387      --------------------------------------------------
ENSG00000166888:ENST0000055615      --------------------------------------------------
ENSG00000166888:ENST0000045407      --------------------------------------------------
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         AGAGCAGACCTGAGGGAGGTTGGTGTGCCGGTCCTGGTTGTCCAGCCCTC
isotig12565                         --------------------------------------------------
isotig12566                         --------------------------------------------------

ENSG00000166888:ENST0000030013      --------------------------------------------------
ENSG00000166888:ENST0000054387      --------------------------------------------------
ENSG00000166888:ENST0000055615      --------------------------------------------------
ENSG00000166888:ENST0000045407      --------------------------------------------------
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         CCGTCCCATCCCAGGGCTGGCTCACCAGAAAAAAGTACAGGCATCCCAGG
isotig12565                         --------------------------------------------------
isotig12566                         --------------------------------------------------

ENSG00000166888:ENST0000030013      --------------------------------------------------
ENSG00000166888:ENST0000054387      --------------------------------------------------
ENSG00000166888:ENST0000055615      --------------------------------------------------
ENSG00000166888:ENST0000045407      --------------------------------------------------
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         CCCTCTGGGGCGGGAGCAGAGGGTCTCCTCTTGTCCGGGGGAAGACAAAA
isotig12565                         --------------------------------------------------
isotig12566                         --------------------------------------------------

ENSG00000166888:ENST0000030013      --------------------------------------------------
ENSG00000166888:ENST0000054387      --------------------------------------------------
ENSG00000166888:ENST0000055615      --------------------------------------------------
ENSG00000166888:ENST0000045407      --------------------------------------------------
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         ACAAAAACAAAAAACGAGTGTGTAGCCCTCACTGGAAGTCTTCTGGTGCT
isotig12565                         --------------------------------------------------
isotig12566                         --------------------------------------------------

ENSG00000166888:ENST0000030013      --------------------------------------------------
ENSG00000166888:ENST0000054387      --------------------------------------------------
ENSG00000166888:ENST0000055615      --------------------------------------------------
ENSG00000166888:ENST0000045407      --------------------------------------------------
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         CTGGGGCCGTTTGCACTTGGGAGCCTGGCTTCTGGGCATGGTGGCCGGGC
isotig12565                         --------------------------------------------------
isotig12566                         --------------------------------------------------

ENSG00000166888:ENST0000030013      --------------------------------------------------
ENSG00000166888:ENST0000054387      --------------------------------------------------
ENSG00000166888:ENST0000055615      --------------------------------------------------
ENSG00000166888:ENST0000045407      --------------------------------------------------
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         TCCGCCGGGGGCTTCAGCCCTCGGCAAAGCGTTTCTAGAACAGAGTGGGT
isotig12565                         -------GGCCAGGCGGCCCTCGGGAGC----------------------
isotig12566                         -------GGCCAGGCGGCCCTCGGGAGC----------------------

ENSG00000166888:ENST0000030013      --------------------------------------------------
ENSG00000166888:ENST0000054387      --------------------------------------------------
ENSG00000166888:ENST0000055615      --------------------------------------------------
ENSG00000166888:ENST0000045407      --------------------------------------------------
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         GGGAGCCGACTGAAGTCCCAGAACCGCGAACACTGACGCACAGGAAAGCC
isotig12565                         -----------------CCAGCTGCTCGCACCTGGAGCCAGCGC------
isotig12566                         -----------------CCAGCTGCTCGCACCTGGAGCCAGCGC------

ENSG00000166888:ENST0000030013      --------------------------------------------------
ENSG00000166888:ENST0000054387      --------------------------------------------------
ENSG00000166888:ENST0000055615      --------------------------------------------------
ENSG00000166888:ENST0000045407      --------------------------------------------------
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         TCGAGTGGGTGGAGAAATGCAAATCCCCCAGTGGCATGGCCTCAGCCGCT
isotig12565                         ------------------------AGCCCGGCCAGTCGGGCCTCAGCCCG
isotig12566                         ------------------------AGCCCGGCCAGTCGGGCCTCAGCCCG

ENSG00000166888:ENST0000030013      --------------------------------------------------
ENSG00000166888:ENST0000054387      --------------------------------------------------
ENSG00000166888:ENST0000055615      --------------------------------------------------
ENSG00000166888:ENST0000045407      --------------------------------------------------
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         GCGGACCCTGACGGCCGATTCTGACTTCGGACTTGAGCGGGACCAGCACC
isotig12565                         GGAGACAGTTACGCCC---------------------------------C
isotig12566                         GGAGACAGTTACGCCC---------------------------------C

ENSG00000166888:ENST0000030013      --------------------ATGTCTCTGTGGGGTCTGGTCTCCAAGATG
ENSG00000166888:ENST0000054387      --------------------ATGTCTCTGTGGGGTCTGGTCTCCAAGATG
ENSG00000166888:ENST0000055615      --------------------ATGTCTCTGTGGGGTCTGGTCTCCAAGATG
ENSG00000166888:ENST0000045407      --------------------ATGTCTCTGTGGGGTCTGGTCTCCAAGATG
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         CGGACCCAGGGCGAGCCAGCATGTCTCAGTGGAATCAAGTCCAACAGTTA
isotig12565                         CTGATTGCGGCAGG------ATGGCCCAGTGGAACCAGCTGCAGCAGCTG
isotig12566                         CTGATTGCGGCAGG------ATGGCCCAGTGGAACCAGCTGCAGCAGCTG

ENSG00000166888:ENST0000030013      CCCCCA---------GAAAAAGTGCAGCGGCTCTATGTCGAC---TTTCC
ENSG00000166888:ENST0000054387      CCCCCA---------GAAAAAGTGCAGCGGCTCTATGTCGAC---TTTCC
ENSG00000166888:ENST0000055615      CCCCCA---------GAAAAAGTGCAGCGGCTCTATGTCGAC---TTTCC
ENSG00000166888:ENST0000045407      CCCCCA---------GAAAAAGTGCAGCGGCTCTATGTCGAC---TTTCC
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         GAAATCAAGTTTTTGGAGCAAGTTGATCAGTTCTATGATGACAACTTTCC
isotig12565                         GACACTCGGTACCTGGAGCAGCTGCACCAGCTGTACAGCGACAGCTTCCC
isotig12566                         GACACTCGGTACCTGGAGCAGCTGCACCAGCTGTACAGCGACAGCTTCCC

ENSG00000166888:ENST0000030013      CCAACACCTGCGGCATCTTCTGGGTGACTGGCTGGAGAGCCAGCCCTGGG
ENSG00000166888:ENST0000054387      CCAACACCTGCGGCATCTTCTGGGTGACTGGCTGGAGAGCCAGCCCTGGG
ENSG00000166888:ENST0000055615      CCAACACCTGCGGCATCTTCTGGGTGACTGGCTGGAGAGCCAGCCCTGGG
ENSG00000166888:ENST0000045407      CCAACACCTGCGGCATCTTCTGGGTGACTGGCTGGAGAGCCAGCCCTGGG
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         CATGGAAATCCGACATCTGCTGGCCCAGTGGATTGAGCATCAAGACTGGG
isotig12565                         CATGGAGCTGCGGCAGTTCCTGGCACCTTGGATTGAGAGTCAAGATTGGG
isotig12566                         CATGGAGCTGCGGCAGTTCCTGGCACCTTGGATTGAGAGTCAAGATTGGG

ENSG00000166888:ENST0000030013      AGTTCCTGGTCGGCTCCGACGCCTTCTGCTGCAACTTGGCTAGTGCCCTA
ENSG00000166888:ENST0000054387      AGTTCCTGGTCGGCTCCGACGCCTTCTGCTGCAACTTGGCTAGTGCCCTA
ENSG00000166888:ENST0000055615      AGTTCCTGGTCGGCTCCGACGCCTTCTGCTGCAACTTGGCTAGTGCCCTA
ENSG00000166888:ENST0000045407      AGTTCCTGGTCGGCTCCGACGCCTTCTGCTGCAACTTGGCTAGTGCCCTA
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         AGGTGGCCTCTAACAATGAAACTATGGCAACAATTCTTCTTCAAAACTTA
isotig12565                         CATATGCAGCCAGCAAAGAGTCACATGCCACACTGGTGTTTCATAATCTC
isotig12566                         CATATGCAGCCAGCAAAGAGTCACATGCCACACTGGTGTTTCATAATCTC

ENSG00000166888:ENST0000030013      CTTTCAGACACTGTCCAGCACCTTCAGGCCTCGGTGGGAGAGCAGGGGGA
ENSG00000166888:ENST0000054387      CTTTCAGACACTGTCCAGCACCTTCAGGCCTCGGTGGGAGAGCAGGGGGA
ENSG00000166888:ENST0000055615      CTTTCAGACACTGTCCAGCACCTTCAGGCCTCGGTGGGAGAGCAGGGGGA
ENSG00000166888:ENST0000045407      CTTTCAGACACTGTCCAGCACCTTCAGGCCTCGGTGGGAGAGCAGGGGGA
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         TTAATACAATTGGATGAACAGTTAGGTCGTGTTTCCAAAGAGAAA-----
isotig12565                         TTGGGTGAGATTGACCAGCAGTACAGCCGATTCCTGCAAGAGTCC-----
isotig12566                         TTGGGTGAGATTGACCAGCAGTACAGCCGATTCCTGCAAGAGTCC-----

ENSG00000166888:ENST0000030013      GGGGAGCACCATCTTGCAACAC---------------ATCAGCACCCTTG
ENSG00000166888:ENST0000054387      GGGGAGCACCATCTTGCAACAC---------------ATCAGCACCCTTG
ENSG00000166888:ENST0000055615      GGGGAGCACCATCTTGCAACAC---------------ATCAGCACCCTTG
ENSG00000166888:ENST0000045407      GGGGAGCACCATCTTGCAACAC---------------ATCAGCACCCTTG
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         ----AACCTGCTATTGATCCACAATCTAAAGAGAATTAGAAAAGTACTTC
isotig12565                         ----AACGTCCTCTATCAGCACAACCTTCGGAGGATCAAGCAGTTCCTAC
isotig12566                         ----AACGTCCTCTATCAGCACAACCTTCGGAGGATCAAGCAGTTCCTAC

ENSG00000166888:ENST0000030013      AGAGCATATATCAGAGGGACCCCCTGAAGCTGGTGGCCACTTTCAGACAA
ENSG00000166888:ENST0000054387      AGAGCATATATCAGAGGGACCCCCTGAAGCTGGTGGCCACTTTCAGACAA
ENSG00000166888:ENST0000055615      AGAGCATATATCAGAGGGACCCCCTGAAGCTGGTGGCCACTTTCAGACAA
ENSG00000166888:ENST0000045407      AGAGCATATATCAGAGGGACCCCCTGAAGCTGGTGGCCACTTTCAGACAA
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         AGGGGAAGTTTCATGGAAATCCAATGCATGTAGCCGTGGTAATCTCAAAT
isotig12565                         AGAGCAGGTATCTTGAGAAGCCGATGGAAATCGCCCGCATCGTGGCCCGA
isotig12566                         AGAGCAGGTATCTTGAGAAGCCGATGGAAATCGCCCGCATCGTGGCCCGA

ENSG00000166888:ENST0000030013      ATACTTCAAGGAGAGAAAAAAGCTGTT-----------------------
ENSG00000166888:ENST0000054387      ATACTTCAAGGAGAGAAAAAAGCTGTT-----------------------
ENSG00000166888:ENST0000055615      ATACTTCAAGGAGAGAAAAAAGCTGTT-----------------------
ENSG00000166888:ENST0000045407      ATACTTCAAGGAGAGAAAAAAGCTGTT-----------------------
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         TGTTTAAGGGAAGAGAGGAGAATACTG---GCTGCAGCGAACATGCCTAT
isotig12565                         TGCCTGTGGGAAGAGTCTCGCCTCCTCCAGACGGCAGCCACTGCAGCCCA
isotig12566                         TGCCTGTGGGAAGAGTCTCGCCTCCTCCAGACGGCAGCCACTGCAGCCCA

ENSG00000166888:ENST0000030013      --------------------------------------------------
ENSG00000166888:ENST0000054387      --------------------------------------------------
ENSG00000166888:ENST0000055615      --------------------------------------------------
ENSG00000166888:ENST0000045407      --------------------------------------------------
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         CCAGGGACCTCTGGAGAAATCCTTACAAAGTTCGTCGGTTTCAGAAAGAC
isotig12565                         GCAAGGAGGCCAGGCCAACCACCCAACAGCTGCTGTGGTGACGGAGAAAC
isotig12566                         GCAAGGAGGCCAGGCCAACCACCCAACAGCTGCTGTGGTGACGGAGAAAC

ENSG00000166888:ENST0000030013      -----------ATGGAACAGTTCCGCCACTTGCCAATGCCTTTCCACTGG
ENSG00000166888:ENST0000054387      -----------ATGGAACAGTTCCGCCACTTGCCAATGCCTTTCCACTGG
ENSG00000166888:ENST0000055615      -----------ATGGAACAGTTCCGCCACTTGCCAATGCCTTTCCACTGG
ENSG00000166888:ENST0000045407      -----------ATGGAACAGTTCCGCCACTTGCCAATGCCTTTCCACTGG
ENSG00000166888:ENST0000053891      -----------ATGGAACAGTTCCGCCACTTGCCAATGCCTTTCCACTGG
ENSG00000166888:ENST0000053721      -----------ATGGAACAGTTCCGCCACTTGCCAATGCCTTTCCACTGG
ENSG00000166888:ENST0000053520      -----------ATGGAACAGTTCCGCCACTTGCCAATGCCTTTCCACTGG
isotig46679                         AGAGAAATGTGGAACACAAAGTGGCTGCCATTAAAAACAGTGTGCAGATG
isotig12565                         AGCAGATGCTGGAGCAGCATCTTCAGGATGTCCGGAAACGTGTGCAGGAT
isotig12566                         AGCAGATGCTGGAGCAGCATCTTCAGGATGTCCGGAAACGTGTGCAGGAT

ENSG00000166888:ENST0000030013      AAGCAGGAAGAACTCAAGTTTAAGACAGGCTTGCGG---AGGCTGCAGCA
ENSG00000166888:ENST0000054387      AAGCAGGAAGAACTCAAGTTTAAGACAGGCTTGCGG---AGGCTGCAGCA
ENSG00000166888:ENST0000055615      AAGCAGGAAGAACTCAAGTTTAAGACAGGCTTGCGG---AGGCTGCAGCA
ENSG00000166888:ENST0000045407      AAGCAGGAAGAACTCAAGTTTAAGACAGGCTTGCGG---AGGCTGCAGCA
ENSG00000166888:ENST0000053891      AAGCAGGAAGAACTCAAGTTTAAGACAGGCTTGCGG---AGGCTGCAGCA
ENSG00000166888:ENST0000053721      AAGCAGGAAGAACTCAAGTTTAAGACAGGCTTGCGG---AGGCTGCAGCA
ENSG00000166888:ENST0000053520      AAGCAGGAAGAACTCAAGTTTAAGACAGGCTTGCGG---AGGCTGCAGCA
isotig46679                         ACAGAACAAGACACCAAATACTTGGAAGATCTGCAAGATGAATTTGACTA
isotig12565                         CTAGAACAGAAAATGAAAGTGGTAGAGAATCTCCAGGATGACTTTGATTT
isotig12566                         CTAGAACAGAAAATGAAAGTGGTAGAGAATCTCCAGGATGACTTTGATTT

ENSG00000166888:ENST0000030013      CCGAGTAGGGGAGATCCACCTTCTCCGAGAAGCCCTGCAGAAGGGGGCTG
ENSG00000166888:ENST0000054387      CCGAGTAGGGGAGATCCACCTTCTCCGAGAAGCCCTGCAGAAGGGGGCTG
ENSG00000166888:ENST0000055615      CCGAGTAGGGGAGATCCACCTTCTCCGAGAAGCCCTGCAGAAGGGGGCTG
ENSG00000166888:ENST0000045407      CCGAGTAGGGGAGATCCACCTTCTCCGAGAAGCCCTGCAGAAGGGGGCTG
ENSG00000166888:ENST0000053891      CCGAGTAGGGGAGATCCACCTTCTCCGAGAAGCCCTGCAGAAGGGGGCTG
ENSG00000166888:ENST0000053721      CCGAGTAGGGGAGATCCACCTTCTCCGAGAAGCCCTGCAGAAGGGGGCTG
ENSG00000166888:ENST0000053520      CCGAGTAGGGGAGATCCACCTTCTCCGAGAAGCCCTGCAGAAGGGGGCTG
isotig46679                         CAGGTATAAAACAATTCAGACA------------ATGGACCAGGGTGACA
isotig12565                         CAACTATAAAACCCTCAAGAGTCAAGGAGACATGCAGGATCTGAATGGAA
isotig12566                         CAACTATAAAACCCTCAAGAGTCAAGGAGACATGCAGGATCTGAATGGAA

ENSG00000166888:ENST0000030013      AGGCTGGCCAAGTGTCTCTGCACAGCTTGATAGAAACTCCTGCTAATGGG
ENSG00000166888:ENST0000054387      AGGCTGGCCAAGTGTCTCTGCACAGCTTGATAGAAACTCCTGCTAATGGG
ENSG00000166888:ENST0000055615      AGGCTGGCCAAGTGTCTCTGCACAGCTTGATAGAAACTCCTGCTAATGGG
ENSG00000166888:ENST0000045407      AGGCTGGCCAAGTGTCTCTGCACAGCTTGATAGAAACTCCTGCTAATGGG
ENSG00000166888:ENST0000053891      AGGCTGGCCAAGTGTCTCTGCACAGCTTGATAGAAACTCCTGCTAATGGG
ENSG00000166888:ENST0000053721      AGGCTGGCCAAGTGTCTCTGCACAGCTTGATAGAAACTCCTGCTAATGGG
ENSG00000166888:ENST0000053520      AGGCTGGCCAAGTGTCTCTGCACAGCTTGATAGAAACTCCTGCTAATGGG
isotig46679                         AGAATAGCATCCTAATGAACCAGGAGGTTTTGACACTCCAAGAAATGCTT
isotig12565                         ACAACCAGTCTGTGACCAGGCAGAAGATGCAGCAGCTGGAACAGATGCTC
isotig12566                         ACAACCAGTCTGTGACCAGGCAGAAGATGCAGCAGCTGGAACAGATGCTC

ENSG00000166888:ENST0000030013      ACTGGGCCAAGTGAGGCCCTGGCCATGCTACTGCAGGAGACCACTGGAGA
ENSG00000166888:ENST0000054387      ACTGGGCCAAGTGAGGCCCTGGCCATGCTACTGCAGGAGACCACTGGAGA
ENSG00000166888:ENST0000055615      ACTGGGCCAAGTGAGGCCCTGGCCATGCTACTGCAGGAGACCACTGGAGA
ENSG00000166888:ENST0000045407      ACTGGGCCAAGTGAGGCCCTGGCCATGCTACTGCAGGAGACCACTGGAGA
ENSG00000166888:ENST0000053891      ACTGGGCCAAGTGAGGCCCTGGCCATGCTACTGCAGGAGACCACTGGAGA
ENSG00000166888:ENST0000053721      ACTGGGCCAAGTGAGGCCCTGGCCATGCTACTGCAGGAGACCACTGGAGA
ENSG00000166888:ENST0000053520      ACTGGGCCAAGTGAGGCCCTGGCCATGCTACTGCAGGAGACCACTGGAGA
isotig46679                         AATAGCCTGGACTTCAAGAGAAAGGAAGCACTCACTAAGATGACACAGAT
isotig12565                         ACGGCGCTGGACCAGATGCGGAGAAGCATTGTGAGTGAGCTGGCGGGGCT
isotig12566                         ACGGCGCTGGACCAGATGCGGAGAAGCATTGTGAGTGAGCTGGCGGGGCT

ENSG00000166888:ENST0000030013      GCTAGAGGCAGCC------------AAAGCCCTAGTGCTGAAGAGGATCC
ENSG00000166888:ENST0000054387      GCTAGAGGCAGCC------------AAAGCCCTAGTGCTGAAGAGGATCC
ENSG00000166888:ENST0000055615      GCTAGAGGCAGCC------------AAAGCCCTAGTGCTGAAGAGGATCC
ENSG00000166888:ENST0000045407      GCTAGAGGCAGCC------------AAAGCCCTAGTGCTGAAGAGGATCC
ENSG00000166888:ENST0000053891      GCTAGAGGCAGCC------------AAAGCCCTAGTGCTGAAGAGGATCC
ENSG00000166888:ENST0000053721      GCTAGAGGCAGCC------------AAAGCCCTAGTGCTGAAGAGGATCC
ENSG00000166888:ENST0000053520      GCTAGAGGCAGCC------------AAAGCCCTAGTGCTGAAGAGGATCC
isotig46679                         AGTGAACGAGTCGGACCTGCTGATGAGCAGCATGCTCATAGAAGAGCTGC
isotig12565                         TTTGTCGGCAATGGAGTACGTGCAGAAAACACTCACAGACGAGGAGCTGG
isotig12566                         TTTGTCGGCAATGGAGTACGTGCAGAAAACACTCACAGACGAGGAGCTGG

ENSG00000166888:ENST0000030013      AGATTTGGAAACGGCAGCAGCAGCTGGCAGGGAATGGCGCACCGTTTGAG
ENSG00000166888:ENST0000054387      AGATTTGGAAACGGCAGCAGCAGCTGGCAGGGAATGGCGCACCGTTTGAG
ENSG00000166888:ENST0000055615      AGATTTGGAAACGGCAGCAGCAGCTGGCAGGGAATGGCGCACCGTTTGAG
ENSG00000166888:ENST0000045407      AGATTTGGAAACGGCAGCAGCAGCTGGCAGGGAATGGCGCACCGTTTGAG
ENSG00000166888:ENST0000053891      AGATTTGGAAACGGCAGCAGCAGCTGGCAGGGAATGGCGCACCGTTTGAG
ENSG00000166888:ENST0000053721      AGATTTGGAAACGGCAGCAGCAGCTGGCAGGGAATGGCGCACCGTTTGAG
ENSG00000166888:ENST0000053520      AGATTTGGAAACGGCAGCAGCAGCTGGCAGGGAATGGCGCACCGTTTGAG
isotig46679                         AGGACTGGAAGAGGAGGCAGCAGATCGCCTGCATCGGTGGCCCACTCCAC
isotig12565                         CTGACTGGAAGAGGCGGCAGCAGATCGCGTGCATTGGAGGCCCTCCCAAC
isotig12566                         CTGACTGGAAGAGGCGGCAGCAGATCGCGTGCATTGGAGGCCCTCCCAAC

ENSG00000166888:ENST0000030013      GAGAGCCTGGCCCCACTCCAGGAGAGGTGTGAAAGCCTGGTGGACATTTA
ENSG00000166888:ENST0000054387      GAGAGCCTGGCCCCACTCCAGGAGAGGTGTGAAAGCCTGGTGGACATTTA
ENSG00000166888:ENST0000055615      GAGAGCCTGGCCCCACTCCAGGAGAGGTGTGAAAGCCTGGTGGACATTTA
ENSG00000166888:ENST0000045407      GAGAGCCTGGCCCCACTCCAGGAGAGGTGTGAAAGCCTGGTGGACATTTA
ENSG00000166888:ENST0000053891      GAGAGCCTGGCCCCACTCCAGGAGAGGTGTGAAAGCCTGGTGGACATTTA
ENSG00000166888:ENST0000053721      GAGAGCCTGGCCCCACTCCAGGAGAGGTGTGAAAGCCTGGTGGACATTTA
ENSG00000166888:ENST0000053520      GAGAGCCTGGCCCCACTCCAGGAGAGGTGTGAAAGCCTGGTGGACATTTA
isotig46679                         AACGGGCTGGACCAGCTTCAGAACTGCTTTACCCTGTTGGCAGAAAGTCT
isotig12565                         ATCTGCCTGGATCGCCTGGAAAACTGGATAACTTCGTTAGCAGAATCTCA
isotig12566                         ATCTGCCTGGATCGCCTGGAAAACTGGATAACTTCGTTAGCAGAATCTCA

ENSG00000166888:ENST0000030013      TTCCCAGCTACAGCAGGAGGTAGGG-------------------------
ENSG00000166888:ENST0000054387      TTCCCAGCTACAGCAGGAGGTAGGG-------------------------
ENSG00000166888:ENST0000055615      TTCCCAGCTACAGCAGGAGGTAGGG-------------------------
ENSG00000166888:ENST0000045407      TTCCCAGCTACAGCAGGAGGTAGGG-------------------------
ENSG00000166888:ENST0000053891      TTCCCAGCTACAGCAGGAGGTAGGG-------------------------
ENSG00000166888:ENST0000053721      TTCCCAGCTACAGCAGGAGGTAGGG-------------------------
ENSG00000166888:ENST0000053520      TTCCCAGCTACAGCAGGAGGTAGGG-------------------------
isotig46679                         TTTCCAACTCAGACGACAGCTGGAGAAATTAGAGGAGCAGTCTTCCAAGA
isotig12565                         ACTTCAGACCCGCCAACAAATTAAGAAACTGGAGGAGCTACAGCAGAAGG
isotig12566                         ACTTCAGACCCGCCAACAAATTAAGAAACTGGAGGAGCTACAGCAGAAGG

ENSG00000166888:ENST0000030013      --GCGGCTGGTGGGGAGCTTGAGCCCAAGACCCGGGCATCGCTGACTGGC
ENSG00000166888:ENST0000054387      --GCGGCTGGTGGGGAGCTTGAGCCCAAGACCCGGGCATCGCTGACTGGC
ENSG00000166888:ENST0000055615      --GCGGCTGGTGGGGAGCTTGAGCCCAAGACCCGGGCATCGCTGACTGGC
ENSG00000166888:ENST0000045407      --GCGGCTGGTGGGGAGCTTGAGCCCAAGACCCGGGCATCGCTGACTGGC
ENSG00000166888:ENST0000053891      --GCGGCTGGTGGGGAGCTTGAGCCCAAGACCCGGGCATCGCTGACTGGC
ENSG00000166888:ENST0000053721      --GCGGCTGGTGGGGAGCTTGAGCCCAAGACCCGGGCATCGCTGACTGGC
ENSG00000166888:ENST0000053520      --GCGGCTGGTGGGGAGCTTGAGCCCAAGACCCGGGCATCGCTGACTGGC
isotig46679                         TGACTTACGAAGGAGACCCCATCCCCACGCAGAGAGCACACCTGCTGGAG
isotig12565                         TGTCCTACAAGGGGGACCCCATTGTGCAGCACCGGCCGATGCTGGAGGAG
isotig12566                         TGTCCTACAAGGGGGACCCCATTGTGCAGCACCGGCCGATGCTGGAGGAG

ENSG00000166888:ENST0000030013      CGGCTGGATGAAGTCCTGAGAACCCTCGTCACCAGTTGCTTCCTGGTGGA
ENSG00000166888:ENST0000054387      CGGCTGGATGAAGTCCTGAGAACCCTCGTCACCAGTTGCTTCCTGGTGGA
ENSG00000166888:ENST0000055615      CGGCTGGATGAAGTCCTGAGAACCCTCGTCACCAGTTGCTTCCTGGTGGA
ENSG00000166888:ENST0000045407      CGGCTGGATGAAGTCCTGAGAACCCTCGTCACCAGTTGCTTCCTGGTGGA
ENSG00000166888:ENST0000053891      CGGCTGGATGAAGTCCTGAGAACCCTCGTCACCAGTTGCTTCCTGGTGGA
ENSG00000166888:ENST0000053721      CGGCTGGATGAAGTCCTGAGAACCCTCGTCACCAGTTGCTTCCTGGTGGA
ENSG00000166888:ENST0000053520      CGGCTGGATGAAGTCCTGAGAACCCTCGTCACCAGTTGCTTCCTGGTGGA
isotig46679                         AGAGCCACCTTCCTGATCTACAACCTTTTCAAGAACTCATTTGTGGTTGA
isotig12565                         CGGATCGTGGAGCTGTTCAGAAACTTGATGAAGAGTGCCTTCGTGGTGGA
isotig12566                         CGGATCGTGGAGCTGTTCAGAAACTTGATGAAGAGTGCCTTCGTGGTGGA

ENSG00000166888:ENST0000030013      GAAGCAGCCC------------------------CCCCAGGTACTGAAGA
ENSG00000166888:ENST0000054387      GAAGCAGCCC------------------------CCCCAGGTACTGAAGA
ENSG00000166888:ENST0000055615      GAAGCAGCCC------------------------CCCCAGGTACTGAAGA
ENSG00000166888:ENST0000045407      GAAGCAGCCC------------------------CCCCAGGTACTGAAGA
ENSG00000166888:ENST0000053891      GAAGCAGCCC------------------------CCCCAGGTACTGAAGA
ENSG00000166888:ENST0000053721      GAAGCAGCCC------------------------CCCCAGGTACTGAAGA
ENSG00000166888:ENST0000053520      GAAGCAGCCC------------------------CCCCAGGTACTGAAGA
isotig46679                         GCGACAGCCCTGCATGCCAACACACCCTCAGAGGCCGCTGGTACTCAAAA
isotig12565                         GCGACAGCCCTGCATGCCGATGCACCCCGACCGGCCCTTGGTCATCAAGA
isotig12566                         GCGACAGCCCTGCATGCCGATGCACCCCGACCGGCCCTTGGTCATCAAGA

ENSG00000166888:ENST0000030013      CTCAGACCAAGTTCCAGGCTGGAGTTCGATTCCTGTTGGGCTTGAGGTTC
ENSG00000166888:ENST0000054387      CTCAGACCAAGTTCCAGGCTGGAGTTCGATTCCTGTTGGGCTTGAGGTTC
ENSG00000166888:ENST0000055615      CTCAGACCAAGTTCCAGGCTGGAGTTCGATTCCTGTTGGGCTTGAGGTTC
ENSG00000166888:ENST0000045407      CTCAGACCAAGTTCCAGGCTGGAGTTCGATTCCTGTTGGGCTTGAGGTTC
ENSG00000166888:ENST0000053891      CTCAGACCAAGTTCCAGGCTGGAGTTCGATTCCTGTTGGGCTTGAGGTTC
ENSG00000166888:ENST0000053721      CTCAGACCAAGTTCCAGGCTGGAGTTCGATTCCTGTTGGGCTTGAGGTTC
ENSG00000166888:ENST0000053520      CTCAGACCAAGTTCCAGGCTGGAGTTCGATTCCTGTTGGGCTTGAGGTTC
isotig46679                         CCCTCATTCAGTTCACCGCGAAACTGAGACTACTAATAAAATTG------
isotig12565                         CTGGTGTCCAGTTCACTACTAAAGTCAGGTTGTTGGTCAAGTTT------
isotig12566                         CTGGTGTCCAGTTCACTACTAAAGTCAGGTTGTTGGTCAAGTTT------

ENSG00000166888:ENST0000030013      CTGGGGGCCCCAGCCAAGCCTCCGCTGGTCAGGGCCGACATGGTGACAGA
ENSG00000166888:ENST0000054387      CTGGGGGCCCCAGCCAAGCCTCCGCTGGTCAGGGCCGACATGGTGACAGA
ENSG00000166888:ENST0000055615      CTGGGGGCCCCAGCCAAGCCTCCGCTGGTCAGGGCCGACATGGTGACAGA
ENSG00000166888:ENST0000045407      CTGGGGGCCCCAGCCAAGCCTCCGCTGGTCAGGGCCGACATGGTGACAGA
ENSG00000166888:ENST0000053891      CTGGGGGCCCCAGCCAAGCCTCCGCTGGTCAGGGCCGACATGGTGACAGA
ENSG00000166888:ENST0000053721      CTGGGGGCCCCAGCCAAGCCTCCGCTGGTCAGGGCCGACATGGTGACAGA
ENSG00000166888:ENST0000053520      CTGGGGGCCCCAGCCAAGCCTCCGCTGGTCAGGGCCGACATGGTGACAGA
isotig46679                         ---------CCGGAACTCAACTATCAGGTGAAAGTAAAGGCATCGATCGA
isotig12565                         ---------CCCGAGTTGAATTATCAGCTTAAAATTAAAGTGTGCATTGA
isotig12566                         ---------CCCGAGTTGAATTATCAGCTTAAAATTAAAGTGTGCATTGA

ENSG00000166888:ENST0000030013      GAAGCAGGCGCGGGAGCTGAGTGTGCCTCAGGGTCCTGGGGCTGGAGCAG
ENSG00000166888:ENST0000054387      GAAGCAGGCGCGGGAGCTGAGTGTGCCTCAGGGTCCTGGGGCTGGAGCAG
ENSG00000166888:ENST0000055615      GAAGCAGGCGCGGGAGCTGAGTGTGCCTCAGGGTCCTGGGGCTGGAGCAG
ENSG00000166888:ENST0000045407      GAAGCAGGCGCGGGAGCTGAGTGTGCCTCAGGGTCCTGGGGCTGGAGCAG
ENSG00000166888:ENST0000053891      GAAGCAGGCGCGGGAGCTGAGTGTGCCTCAGGGTCCTGGGGCTGGAGCAG
ENSG00000166888:ENST0000053721      GAAGCAGGCGCGGGAGCTGAGTGTGCCTCAGGGTCCTGGGGCTGGAGCAG
ENSG00000166888:ENST0000053520      GAAGCAGGCGCGGGAGCTGAGTGTGCCTCAGGGTCCTGGGGCTGGAGCAG
isotig46679                         CAAGAATGTTTCAACGCTAAGC------------AATAGAAGATTTGTGC
isotig12565                         CAAAGATTCTGGGGACGTTGCTGCTCTCAGAGGATCTCGGAAATTTAACA
isotig12566                         CAAAGATTCTGGGGACGTTGCTGCTCTCAGAGGATCTCGGAAATTTAACA

ENSG00000166888:ENST0000030013      AAAGCACTGGAGAAATCATCAACAACACTGTGCCCTTGGAGAACAGCATT
ENSG00000166888:ENST0000054387      AAAGCACTGGAGAAATCATCAACAACACTGTGCCCTTGGAGAACAGCATT
ENSG00000166888:ENST0000055615      AAAGCACTGGAGAAATCATCAACAACACTGTGCCCTTGGAGAACAGCATT
ENSG00000166888:ENST0000045407      AAAGCACTGGAGAAATCATCAACAACACTGTGCCCTTGGAGAACAGCATT
ENSG00000166888:ENST0000053891      AAAGCACTGGAGAAATCATCAACAACACTGTGCCCTTGGAGAACAGCATT
ENSG00000166888:ENST0000053721      AAAGCACTGGAGAAATCATCAACAACACTGTGCCCTTGGAGAACAGCATT
ENSG00000166888:ENST0000053520      AAAGCACTGGAGAAATCATCAACAACACTGTGCCCTTGGAGAACAGCATT
isotig46679                         TTTGTGGAACTCAAGTC------AAAGCCATGTCCATCGAGGAATCCTCC
isotig12565                         TTCTGGGCACAAACACG------AAGGTGATGAACATGGAAGAATCCAAC
isotig12566                         TTCTGGGCACAAACACG------AAGGTGATGAACATGGAAGAATCCAAC

ENSG00000166888:ENST0000030013      CCTGGGAACTGCTGCTCTGCCCTGTTCAAGAACCTGCTTCTCAAGAAGAT
ENSG00000166888:ENST0000054387      CCTGGGAACTGCTGCTCTGCCCTGTTCAAGAACCTGCTTCTCAAGAAGAT
ENSG00000166888:ENST0000055615      CCTGGGAACTGCTGCTCTGCCCTGTTCAAGAACCTGCTTCTCAAGAAGAT
ENSG00000166888:ENST0000045407      CCTGGGAACTGCTGCTCTGCCCTGTTCAAGAACCTGCTTCTCAAGAAGAT
ENSG00000166888:ENST0000053891      CCTGGGAACTGCTGCTCTGCCCTGTTCAAGAACCTGCTTCTCAAGAAGAT
ENSG00000166888:ENST0000053721      CCTGGGAACTGCTGCTCTGCCCTGTTCAAGAACCTGCTTCTCAAGAAGAT
ENSG00000166888:ENST0000053520      CCTGGGAACTGCTGCTCTGCCCTGTTCAAGAACCTGCTTCTCAAGAAGAT
isotig46679                         AATGGGAGCCTCTCAGTAGAA---TTTAGACATTTGCAACCGAAGGAAAT
isotig12565                         AACGGCAGCCTGTCTGCGGAG---TTCAAGCACTTGACCCTGAGGGAGCA
isotig12566                         AACGGCAGCCTGTCTGCGGAG---TTCAAGCACTTGACCCTGAGGGAGCA

ENSG00000166888:ENST0000030013      CAAG---------------CGGTGTGAGCGGAAGGGCACTGAGTCTGTCA
ENSG00000166888:ENST0000054387      CAAG---------------CGGTGTGAGCGGAAGGGCACTGAGTCTGTCA
ENSG00000166888:ENST0000055615      CAAG---------------CGGTGTGAGCGGAAGGGCACTGAGTCTGTCA
ENSG00000166888:ENST0000045407      CAAG---------------CGGTGTGAGCGGAAGGGCACTGAGTCTGTCA
ENSG00000166888:ENST0000053891      CAAG---------------CGGTGTGAGCGGAAGGGCACTGAGTCTGTCA
ENSG00000166888:ENST0000053721      CAAG---------------CGGTGTGAGCGGAAGGGCACTGAGTCTGTCA
ENSG00000166888:ENST0000053520      CAAG---------------CGGTGTGAGCGGAAGGGCACTGAGTCTGTCA
isotig46679                         GAAATCCAGTGCCGGAAGTAAAGGAAAT---GAGGGCTGCCACATGGTGA
isotig12565                         GAGATGTGGGAATGGAGGCCGTGCCAATTGTGATGCCTCCTTGATTGTGA
isotig12566                         GAGATGTGGGAATGGAGGCCGTGCCAATTGTGATGCCTCCTTGATTGTGA

ENSG00000166888:ENST0000030013      CAGAGGAGAAGTGCGCTGTGCTCTTCTCTGCCAGCTTCACACTTGGCCCC
ENSG00000166888:ENST0000054387      CAGAGGAGAAGTGCGCTGTGCTCTTCTCTGCCAGCTTCACACTTGGCCCC
ENSG00000166888:ENST0000055615      CAGAGGAGAAGTGCGCTGTGCTCTTCTCTGCCAGCTTCACACTTGGCCCC
ENSG00000166888:ENST0000045407      CAGAGGAGAAGTGCGCTGTGCTCTTCTCTGCCAGCTTCACACTTGGCCCC
ENSG00000166888:ENST0000053891      CAGAGGAGAAGTGCGCTGTGCTCTTCTCTGCCAGCTTCACACTTGGCCCC
ENSG00000166888:ENST0000053721      CAGAGGAGAAGTGCGCTGTGCTCTTCTCTGCCAGCTTCACACTTGGCCCC
ENSG00000166888:ENST0000053520      CAGAGGAGAAGTGCGCTGTGCTCTTCTCTGCCAGCTTCACACTTGGCCCC
isotig46679                         CGGAAGAGCTGCATTCCATAGCCTTTGAGACCCAGATCTGCCTC------
isotig12565                         CCGAGGAGCTGCATCTGATCACCTTCGAGACTGAGGTGTACCAC------
isotig12566                         CCGAGGAGCTGCATCTGATCACCTTCGAGACTGAGGTGTACCAC------

ENSG00000166888:ENST0000030013      GGCAAACTCCCCATCCAGCTCCAGGCCCTGTCTCTGCCCCTGGTGGTCAT
ENSG00000166888:ENST0000054387      GGCAAACTCCCCATCCAGCTCCAGGCCCTGTCTCTGCCCCTGGTGGTCAT
ENSG00000166888:ENST0000055615      GGCAAACTCCCCATCCAGCTCCAGGCCCTGTCTCTGCCCCTGGTGGTCAT
ENSG00000166888:ENST0000045407      GGCAAACTCCCCATCCAGCTCCAGGCCCTGTCTCTGCCCCTGGTGGTCAT
ENSG00000166888:ENST0000053891      GGCAAACTCCCCATCCAGCTCCAGGCCCTGTCTCTGCCCCTGGTGGTCAT
ENSG00000166888:ENST0000053721      GGCAAACTCCCCATCCAGCTCCAGGCCCTGTCTCTGCCCCTGGTGGTCAT
ENSG00000166888:ENST0000053520      GGCAAACTCCCCATCCAGCTCCAGGCCCTGTCTCTGCCCCTGGTGGTCAT
isotig46679                         TATGGCCTCACCATCGACTTGGAGACAAGCTCATTACCTGTGGTGATGAT
isotig12565                         CAAGGCCTCAAGATTGACCTGGAGACCCATTCTTTGCCAGTTGTGGTGAT
isotig12566                         CAAGGCCTCAAGATTGACCTGGAGACCCATTCTTTGCCAGTTGTGGTGAT

ENSG00000166888:ENST0000030013      CGTCCATGGCAACCAAGACAACAATGCCAAAGCCACTATCCTGTGGGACA
ENSG00000166888:ENST0000054387      CGTCCATGGCAACCAAGACAACAATGCCAAAGCCACTATCCTGTGGGACA
ENSG00000166888:ENST0000055615      CGTCCATGGCAACCAAGACAACAATGCCAAAGCCACTATCCTGTGGGACA
ENSG00000166888:ENST0000045407      CGTCCATGGCAACCAAGACAACAATGCCAAAGCCACTATCCTGTGGGACA
ENSG00000166888:ENST0000053891      CGTCCATGGCAACCAAGACAACAATGCCAAAGCCACTATCCTGTGGGACA
ENSG00000166888:ENST0000053721      CGTCCATGGCAACCAAGACAACAATGCCAAAGCCACTATCCTGTGGGACA
ENSG00000166888:ENST0000053520      CGTCCATGGCAACCAAGACAACAATGCCAAAGCCACTATCCTGTGGGACA
isotig46679                         TTCTAATGTCAGCCAACTGCCTAATGCTTGGGCATCCATCATTTGGTACA
isotig12565                         CTCCAACATCTGTCAGATGCCAAATGCCTGGGCATCCATCCTGTGGTATA
isotig12566                         CTCCAACATCTGTCAGATGCCAAATGCCTGGGCATCCATCCTGTGGTATA

ENSG00000166888:ENST0000030013      ATGCCTTCTCTGAG------ATGGACCGCGTGCCCTTTGTGGTGGCTGAG
ENSG00000166888:ENST0000054387      ATGCCTTCTCTGAG------ATGGACCGCGTGCCCTTTGTGGTGGCTGAG
ENSG00000166888:ENST0000055615      ATGCCTTCTCTGAG------ATGGACCGCGTGCCCTTTGTGGTGGCTGAG
ENSG00000166888:ENST0000045407      ATGCCTTCTCTGAG------ATGGACCGCGTGCCCTTTGTGGTGGCTGAG
ENSG00000166888:ENST0000053891      ATGCCTTCTCTGAG------ATGGACCGCGTGCCCTTTGTGGTGGCTGAG
ENSG00000166888:ENST0000053721      ATGCCTTCTCTGAG------ATGGACCGCGTGCCCTTTGTGGTGGCTGAG
ENSG00000166888:ENST0000053520      ATGCCTTCTCTGAG------ATGGACCGCGTGCCCTTTGTGGTGGCTGAG
isotig46679                         ATGTGTCAACCAACGATTGCCAGAACTTGGTTTTCTTTAATAATCCTCCG
isotig12565                         ACATGCTGACCAACAACCCCAAGAACGTGAACTTCTTCACCAAGCCACCA
isotig12566                         ACATGCTGACCAACAACCCCAAGAACGTGAACTTCTTCACCAAGCCACCA

ENSG00000166888:ENST0000030013      CGGGTGCCCTGGGAGAAGATGTGTGAAACTCTGAACCTGAAGTTCATGGC
ENSG00000166888:ENST0000054387      CGGGTGCCCTGGGAGAAGATGTGTGAAACTCTGAACCTGAAGTTCATGGC
ENSG00000166888:ENST0000055615      CGGGTGCCCTGGGAGAAGATGTGTGAAACTCTGAACCTGAAGTTCATGGC
ENSG00000166888:ENST0000045407      CGGGTGCCCTGGGAGAAGATGTGTGAAACTCTGAACCTGAAGTTCATGGC
ENSG00000166888:ENST0000053891      CGGGTGCCCTGGGAGAAGATGTGTGAAACTCTGAACCTGAAGTTCATGGC
ENSG00000166888:ENST0000053721      CGGGTGCCCTGGGAGAAGATGTGTGAAACTCTGAACCTGAAGTTCATGGC
ENSG00000166888:ENST0000053520      CGGGTGCCCTGGGAGAAGATGTGTGAAACTCTGAACCTGAAGTTCATGGC
isotig46679                         CCTGTCACTTTGAGTCAACTCCTGGAAGTGATGAGCTGGCAGTTTTCATC
isotig12565                         ATCGGAACCTGGGACCAGGTGGCCGAGGTGCTCAGCTGGCAGTTCTCATC
isotig12566                         ATCGGAACCTGGGACCAGGTGGCCGAGGTGCTCAGCTGGCAGTTCTCATC

ENSG00000166888:ENST0000030013      TGAGGTGGGGACCAACCGGGGGCTGCTCCCAGAGCACTTCCTCTTCCTGG
ENSG00000166888:ENST0000054387      TGAGGTGGGGACCAACCGGGGGCTGCTCCCAGAGCACTTCCTCTTCCTGG
ENSG00000166888:ENST0000055615      TGAGGTGGGGACCAACCGGGGGCTGCTCCCAGAGCACTTCCTCTTCCTGG
ENSG00000166888:ENST0000045407      TGAGGTGGGGACCAACCGGGGGCTGCTCCCAGAGCACTTCCTCTTCCTGG
ENSG00000166888:ENST0000053891      TGAGGTGGGGACCAACCGGGGGCTGCTCCCAGAGCACTTCCTCTTCCTGG
ENSG00000166888:ENST0000053721      TGAGGTGGGGACCAACCGGGGGCTGCTCCCAGAGCACTTCCTCTTCCTGG
ENSG00000166888:ENST0000053520      TGAGGTGGGGACCAACCGGGGGCTGCTCCCAGAGCACTTCCTCTTCCTGG
isotig46679                         CTATGTTGGT------CGTGGCCTTAATTCAGACCAGCTCAACATGCTGG
isotig12565                         CACCACAAAG------CGAGGGCTGAGCATCGAGCAGCTGACTACGCTGG
isotig12566                         CACCACAAAG------CGAGGGCTGAGCATCGAGCAGCTGACTACGCTGG

ENSG00000166888:ENST0000030013      CCCAGAAGATCTTCAATGACAACAGCCTCAGTATGGAGGCCTTCCAGCAC
ENSG00000166888:ENST0000054387      CCCAGAAGATCTTCAATGACAACAGCCTCAGTATGGAGGCCTTCCAGCAC
ENSG00000166888:ENST0000055615      CCCAGAAGATCTTCAATGACAACAGCCTCAGTATGGAGGCCTTCCAGCAC
ENSG00000166888:ENST0000045407      CCCAGAAGATCTTCAATGACAACAGCCTCAGTATGGAGGCCTTCCAGCAC
ENSG00000166888:ENST0000053891      CCCAGAAGATCTTCAATGACAACAGCCTCAGTATGGAGGCCTTCCAGCAC
ENSG00000166888:ENST0000053721      CCCAGAAGATCTTCAATGACAACAGCCTCAGTATGGAGGCCTTCCAGCAC
ENSG00000166888:ENST0000053520      CCCAGAAGATCTTCAATGACAACAGCCTCAGTATGGAGGCCTTCCAGCAC
isotig46679                         CAGAGAAGCTCACAGTTCAGTCT---------------AACTACAGCGAT
isotig12565                         CCGAGAAGCTCCTAGGACCTGGTGTC------------AACTACTCCGGG
isotig12566                         CCGAGAAGCTCCTAGGACCTGGTGTC------------AACTACTCCGGG

ENSG00000166888:ENST0000030013      CGTTCTGTGTCCTGGTCGCAGTTCAACAAGGAGATCCTGCTGGGCCGTGG
ENSG00000166888:ENST0000054387      CGTTCTGTGTCCTGGTCGCAGTTCAACAAGGAGATCCTGCTGGGCCGTGG
ENSG00000166888:ENST0000055615      CGTTCTGTGTCCTGGTCGCAGTTCAACAAGGAGATCCTGCTGGGCCGTGG
ENSG00000166888:ENST0000045407      CGTTCTGTGTCCTGGTCGCAGTTCAACAAGGAGATCCTGCTGGGCCGTGG
ENSG00000166888:ENST0000053891      CGTTCTGTGTCCTGGTCGCAGTTCAACAAGGAGATCCTGCTGGGCCGTGG
ENSG00000166888:ENST0000053721      CGTTCTGTGTCCTGGTCGCAGTTCAACAAGGAGATCCTGCTGGGCCGTGG
ENSG00000166888:ENST0000053520      CGTTCTGTGTCCTGGTCGCAGTTCAACAAGGAGATCCTGCTGGGCCGTGG
isotig46679                         GGTCACCTCACCTGGGCCAAGTTCTGCAAGGAACACTTGCCTGGCAAACC
isotig12565                         TGTCAGATCACATGGGCTAAATTTTGCAAAGAAAACATGGCTGGCAAGGG
isotig12566                         TGTCAGATCACATGGGCTAAATTTTGCAAAGAAAACATGGCTGGCAAGGG

ENSG00000166888:ENST0000030013      CTTCACCTTTTGGCAGTGGTTTGATGGTGTCCTGGACCTCACCAAACGCT
ENSG00000166888:ENST0000054387      CTTCACCTTTTGGCAGTGGTTTGATGGTGTCCTGGACCTCACCAAACGCT
ENSG00000166888:ENST0000055615      CTTCACCTTTTGGCAGTGGTTTGATGGTGTCCTGGACCTCACCAAACGCT
ENSG00000166888:ENST0000045407      CTTCACCTTTTGGCAGTGGTTTGATGGTGTCCTGGACCTCACCAAACGCT
ENSG00000166888:ENST0000053891      CTTCACCTTTTGGCAGTGGTTTGATGGTGTCCTGGACCTCACCAAACGCT
ENSG00000166888:ENST0000053721      CTTCACCTTTTGGCAGTGGTTTGATGGTGTCCTGGACCTCACCAAACGCT
ENSG00000166888:ENST0000053520      CTTCACCTTTTGGCAGTGGTTTGATGGTGTCCTGGACCTCACCAAACGCT
isotig46679                         ATTTACCTTCTGGACCTGGCTTGAAGCAATATTGGACCTAATTAAAAAAC
isotig12565                         CTTCTCCTTCTGGGTGTGGCTAGACAATATCATTGACCTTGTGAAAAAGT
isotig12566                         CTTCTCCTTCTGGGTGTGGCTAGACAATATCATTGACCTTGTGAAAAAGT

ENSG00000166888:ENST0000030013      GTCTCCGGAGCTACTGGTCTGACCGGCTGATCATTGGCTTCATCAGCAAA
ENSG00000166888:ENST0000054387      GTCTCCGGAGCTACTGGTCTGACCGGCTGATCATTGGCTTCATCAGCAAA
ENSG00000166888:ENST0000055615      GTCTCCGGAGCTACTGGTCTGACCGGCTGATCATTGGCTTCATCAGCAAA
ENSG00000166888:ENST0000045407      GTCTCCGGAGCTACTGGTCTGACCGGCTGATCATTGGCTTCATCAGCAAA
ENSG00000166888:ENST0000053891      GTCTCCGGAGCTACTGGTCTGACCGGCTGATCATTGGCTTCATCAGCAAA
ENSG00000166888:ENST0000053721      GTCTCCGGAGCTACTGGTCTGACCGGCTGATCATTGGCTTCATCAGCAAA
ENSG00000166888:ENST0000053520      GTCTCCGGAGCTACTGGTCTGACCGGCTGATCATTGGCTTCATCAGCAAA
isotig46679                         ACATTCTTCCCCTCTGGATTGATGGGTACATCATGGGCTTCGTGAGCAAA
isotig12565                         ATATCTTGGCCCTCTGGAATGAAGGGTACATCATGGGCTTCATTAGCAAG
isotig12566                         ATATCTTGGCCCTCTGGAATGAAGGGTACATCATGGGCTTCATTAGCAAG

ENSG00000166888:ENST0000030013      CAGTACGTTACTAGCCTTCTTCTCAATGAGCCCGACGGAACCTTTCTCCT
ENSG00000166888:ENST0000054387      CAGTACGTTACTAGCCTTCTTCTCAATGAGCCCGACGGAACCTTTCTCCT
ENSG00000166888:ENST0000055615      CAGTACGTTACTAGCCTTCTTCTCAATGAGCCCGACGGAACCTTTCTCCT
ENSG00000166888:ENST0000045407      CAGTACGTTACTAGCCTTCTTCTCAATGAGCCCGACGGAACCTTTCTCCT
ENSG00000166888:ENST0000053891      CAGTACGTTACTAGCCTTCTTCTCAATGAGCCCGACGGAACCTTTCTCCT
ENSG00000166888:ENST0000053721      CAGTACGTTACTAGCCTTCTTCTCAATGAGCCCGACGGAACCTTTCTCCT
ENSG00000166888:ENST0000053520      CAGTACGTTACTAGCCTTCTTCTCAATGAGCCCGACGGAACCTTTCTCCT
isotig46679                         GAGAAGGAGAGGTTTCTGCTCAAGGATAAAATGCCCGGGACATTTTTGTT
isotig12565                         GAGCGGGAGCGGGCGATCCTGAGCACGAAACCCCCGGGCACCTTCCTGCT
isotig12566                         GAGCGGGAGCGGGCGATCCTGAGCACGAAACCCCCGGGCACCTTCCTGCT

ENSG00000166888:ENST0000030013      CCGCTTCAGCGAC---TCAGAGATTGGGGGCATCACCATTGCCCATGTCA
ENSG00000166888:ENST0000054387      CCGCTTCAGCGAC---TCAGAGATTGGGGGCATCACCATTGCCCATGTCA
ENSG00000166888:ENST0000055615      CCGCTTCAGCGAC---TCAGAGATTGGGGGCATCACCATTGCCCATGTCA
ENSG00000166888:ENST0000045407      CCGCTTCAGCGAC---TCAGAGATTGGGGGCATCACCATTGCCCATGTCA
ENSG00000166888:ENST0000053891      CCGCTTCAGCGAC---TCAGAGATTGGGGGCATCACCATTGCCCATGTCA
ENSG00000166888:ENST0000053721      CCGCTTCAGCGAC---TCAGAGATTGGGGGCATCACCATTGCCCATGTCA
ENSG00000166888:ENST0000053520      CCGCTTCAGCGAC---TCAGAGATTGGGGGCATCACCATTGCCCATGTCA
isotig46679                         ACGATTCAGTGAG---AGCCATCTCGGAGGGATCACCTTCACCTGGGTGG
isotig12565                         GAGATTCAGCGAGAGCAGCAAAGAAGGAGGGGTCACTTTCACTTGGGTGG
isotig12566                         GAGATTCAGCGAGAGCAGCAAAGAAGGAGGGGTCACTTTCACTTGGGTGG

ENSG00000166888:ENST0000030013      TCCGGGGCCAGGATGGCTCTCCACAGATAGAGAACATCCAGCCATTCTCT
ENSG00000166888:ENST0000054387      TCCGGGGCCAGGATGGCTCTCCACAGATAGAGAACATCCAGCCATTCTCT
ENSG00000166888:ENST0000055615      TCCGGGGCCAGGATGGCTCTCCACAGATAGAGAACATCCAGCCATTCTCT
ENSG00000166888:ENST0000045407      TCCGGGGCCAGGATGGCTCTCCACAGATAGAGAACATCCAGCCATTCTCT
ENSG00000166888:ENST0000053891      TCCGGGGCCAGGATGGCTCTCCACAGATAGAGAACATCCAGCCATTCTCT
ENSG00000166888:ENST0000053721      TCCGGGGCCAGGATGGCTCTCCACAGATAGAGAACATCCAGCCATTCTCT
ENSG00000166888:ENST0000053520      TCCGGGGCCAGGATGGCTCTCCACAGATAGAGAACATCCAGCCATTCTCT
isotig46679                         ACCACTCTGAAAACGGAGAAGTGAGATTCCACTCCGTAGAACCCTACAAC
isotig12565                         AAAAGGACATCAGTGGCAAGACCCAGATCCAGTCTGTAGAGCCGTACACC
isotig12566                         AAAAGGACATCAGTGGCAAGACCCAGATCCAGTCTGTAGAGCCGTACACC

ENSG00000166888:ENST0000030013      GCCAAAGACCTGTCCATTCGCTCACTGGGGGACCGAATCCGGGAT-----
ENSG00000166888:ENST0000054387      GCCAAAGACCTGTCCATTCGCTCACTGGGGGACCGAATCCGGGAT-----
ENSG00000166888:ENST0000055615      GCCAAAGACCTGTCCATTCGCTCACTGGGGGACCGAATCCGGGAT-----
ENSG00000166888:ENST0000045407      GCCAAAGACCTGTCCATTCGCTCACTGGGGGACCGAATCCGGGAT-----
ENSG00000166888:ENST0000053891      GCCAAAGACCTGTCCATTCGCTCACTGGGGGACCGAATCCGGGAT-----
ENSG00000166888:ENST0000053721      GCCAAAGACCTGTCCATTCGCTCACTGGGGGACCGAATCCGGGAT-----
ENSG00000166888:ENST0000053520      GCCAAAGACCTGTCCATTCGCTCACTGGGGGACCGAATCCGGGAT-----
isotig46679                         AAAGGGCGTCTGTCGGCCCTGCCATTTGCTGACATCCTGCGGGACTACAA
isotig12565                         AAGCAGCAGCTGAACAACATGTCCTTTGCTGAAATCATCATGGGCTACAA
isotig12566                         AAGCAGCAGCTGAACAACATGTCCTTTGCTGAAATCATCATGGGCTACAA

ENSG00000166888:ENST0000030013      -------------------------CTTGCTCAGCTCAAAAATCTCTATC
ENSG00000166888:ENST0000054387      -------------------------CTTGCTCAGCTCAAAAATCTCTATC
ENSG00000166888:ENST0000055615      -------------------------CTTGCTCAGCTCAAAAATCTCTATC
ENSG00000166888:ENST0000045407      -------------------------CTTGCTCAGCTCAAAAATCTCTATC
ENSG00000166888:ENST0000053891      -------------------------CTTGCTCAGCTCAAAAATCTCTATC
ENSG00000166888:ENST0000053721      -------------------------CTTGCTCAGCTCAAAAATCTCTATC
ENSG00000166888:ENST0000053520      -------------------------CTTGCTCAGCTCAAAAATCTCTATC
isotig46679                         GGTCATCATGGCTGAGAACATTCCCGAGAACCCTCTCAAGTACCTCTACC
isotig12565                         GATCATGGATGCCACCAACATCCTGGTGTCCCCATTGGTCTACCTCTACC
isotig12566                         GATCATGGATGCCACCAACATCCTGGTGTCCCCATTGGTCTACCTCTACC

ENSG00000166888:ENST0000030013      CCAAGAAGCCCAAGGATGAGGCTTTCCGGAGCCACTAC------------
ENSG00000166888:ENST0000054387      CCAAGAAGCCCAAGGATGAGGCTTTCCGGAGCCACTAC------------
ENSG00000166888:ENST0000055615      CCAAGAAGCCCAAGGATGAGGCTTTCCGGAGCCACTAC------------
ENSG00000166888:ENST0000045407      CCAAGAAGCCCAAGGATGAGGCTTTCCGGAGCCACTAC------------
ENSG00000166888:ENST0000053891      CCAAGAAGCCCAAGGATGAGGCTTTCCGGAGCCACTAC------------
ENSG00000166888:ENST0000053721      CCAAGAAGCCCAAGGATGAGGCTTTCCGGAGCCACTAC------------
ENSG00000166888:ENST0000053520      CCAAGAAGCCCAAGGATGAGGCTTTCCGGAGCCACTAC------------
isotig46679                         CCGACATCCCCAAAGACAAAGCCTTCGGTAAACACTACAGCTCCCAGCCT
isotig12565                         CTGACATTCCCAAGGAGGAGGCGTTCGGGAAGTACTGT------------
isotig12566                         CTGACATTCCCAAGGAGGAGGCGTTCGGGAAGTACTGT------------

ENSG00000166888:ENST0000030013      ------------AAGCCTGAACAGATGGGTAAGGATGGCAGGGGTTATGT
ENSG00000166888:ENST0000054387      ------------AAGCCTGAACAGATGGGTAAGGATGGCAGGGGTTATGT
ENSG00000166888:ENST0000055615      ------------AAGCCTGAACAGATGGGTAAGGATGGCAGGGGTTATGT
ENSG00000166888:ENST0000045407      ------------AAGCCTGAACAGATGGGTAAGGATGGCAGGGGTTATGT
ENSG00000166888:ENST0000053891      ------------AAGCCTGAACAGATGGGTAAGGATGGCAGGGGTTATGT
ENSG00000166888:ENST0000053721      ------------AAGCCTGAACAGATGGGTAAGGATGGCAGGGGTTATGT
ENSG00000166888:ENST0000053520      ------------AAGCCTGAACAGATGGGTAAGGATGGCAGGGGTTATGT
isotig46679                         TGCGAAGTTTCAAGGCCAACA------GAACGGGGAGACAAAGGTTATGT
isotig12565                         ------------CGACCAGAG-----------------------------
isotig12566                         ------------CGACCAGAG-----------------------------

ENSG00000166888:ENST0000030013      CCCAGCTACCATCAAGATGACCGTGGAAAGGGACCAACCACTTCCTACCC
ENSG00000166888:ENST0000054387      CCCAGCTACCATCAAGATGACCGTGGAAAGGGACCAACCACTTCCTACCC
ENSG00000166888:ENST0000055615      CCCAGCTACCATCAAGATGACCGTGGAAAGGGACCAACCACTTCCTACCC
ENSG00000166888:ENST0000045407      CCCAGCTACCATCAAGATGACCGTGGAAAGGGACCAACCACTTCCTACCC
ENSG00000166888:ENST0000053891      CCCAGCTACCATCAAGATGACCGTGGAAAGGGACCAACCACTTCCTACCC
ENSG00000166888:ENST0000053721      CCCAGCTACCATCAAGATGACCGTGGAAAGGGACCAACCACTTCCTACCC
ENSG00000166888:ENST0000053520      CCCAGCTACCATCAAGATGACCGTGGAAAGGGACCAACCACTTCCTACCC
isotig46679                         TCCTTCAGTTTTTATCCCTATTTCAACAATCCGCAGCGACGCCATGGAGC
isotig12565                         -------------------------------AGCCAGGAGCATCCTGAAG
isotig12566                         -------------------------------AGCCAGGAGCATCCTGAAG

ENSG00000166888:ENST0000030013      CAGAGCTCCAGATGCCTACCATGGTGCCTTCTTATGACCTTGGAATGGCC
ENSG00000166888:ENST0000054387      CAGAGCTCCAGATGCCTACCATGGTGCCTTCTTATGACCTTGGAATGGCC
ENSG00000166888:ENST0000055615      CAGAGCTCCAGATGCCTACCATGGTGCCTTCTTATGACCTTGGAATGGCC
ENSG00000166888:ENST0000045407      CAGAGCTCCAGATGCCTACCATGGTGCCTTCTTATGACCTTGGAATGGCC
ENSG00000166888:ENST0000053891      CAGAGCTCCAGATGCCTACCATGGTGCCTTCTTATGACCTTGGAATGGCC
ENSG00000166888:ENST0000053721      CAGAGCTCCAGATGCCTACCATGGTGCCTTCTTATGACCTTGGAATGGCC
ENSG00000166888:ENST0000053520      CAGAGCTCCAGATGCCTACCATGGTGCCTTCTTATGACCTTGGAATGGCC
isotig46679                         CGCAG------------------------------------------TCT
isotig12565                         CTGAC---------------------------------------------
isotig12566                         CTGAC---------------------------------------------

ENSG00000166888:ENST0000030013      CCTGATTCCTCCATGAGCATGCAGCTTGGCCCAGATATGGTGCCCCAGGT
ENSG00000166888:ENST0000054387      CCTGATTCCTCCATGAGCATGCAGCTTGGCCCAGATATGGTGCCCCAGGT
ENSG00000166888:ENST0000055615      CCTGATTCCTCCATGAGCATGCAGCTTGGCCCAGATATGGTGCCCCAGGT
ENSG00000166888:ENST0000045407      CCTGATTCCTCCATGAGCATGCAGCTTGGCCCAGATATGGTGCCCCAGGT
ENSG00000166888:ENST0000053891      CCTGATTCCTCCATGAGCATGCAGCTTGGCCCAGATATGGTGCCCCAGGT
ENSG00000166888:ENST0000053721      CCTGATTCCTCCATGAGCATGCAGCTTGGCCCAGATATGGTGCCCCAGGT
ENSG00000166888:ENST0000053520      CCTGATTCCTCCATGAGCATGCAGCTTGGCCCAGATATGGTGCCCCAGGT
isotig46679                         CCTTCAGACCTTCTCCCC------------------ATGTCTCCGAGTGT
isotig12565                         CCCGGTAGTGCCGCCCCT------------------TACCTGAAGACCAA
isotig12566                         CCCGGTAGTTGTTTTTCCATG-----------------------------

ENSG00000166888:ENST0000030013      GTACCCACCACACTCTCACTCCATCCCCCCGTATCAAGGCCTCTCCCCAG
ENSG00000166888:ENST0000054387      GTACCCACCACACTCTCACTCCATCCCCCCGTATCAAGGCCTCTCCCCAG
ENSG00000166888:ENST0000055615      GTACCCACCACACTCTCACTCCATCCCCCCGTATCAAGGCCTCTCCCCAG
ENSG00000166888:ENST0000045407      GTACCCACCACACTCTCACTCCATCCCCCCGTATCAAGGCCTCTCCCCAG
ENSG00000166888:ENST0000053891      GTACCCACCACACTCTCACTCCATCCCCCCGTATCAAGGCCTCTCCCCAG
ENSG00000166888:ENST0000053721      GTACCCACCACACTCTCACTCCATCCCCCCGTATCAAGGCCTCTCCCCAG
ENSG00000166888:ENST0000053520      GTACCCACCACACTCTCACTCCATCCCCCCGTATCAAGGCCTCTCCCCAG
isotig46679                         ATACGCTGTGCTGAGAGAAAACCTGAGCCCT-------------------
isotig12565                         GTTCATCTGTGTGACACCAACGACCTGCAGC-------------------
isotig12566                         --------------------------------------------------

ENSG00000166888:ENST0000030013      AAGAATCAGTCAACGTGTTGTCAGCCTTCCAGGAGCCTCACCTGCAGATG
ENSG00000166888:ENST0000054387      AAGAATCAGTCAACGTGTTGTCAGCCTTCCAGGAGCCTCACCTGCAGATG
ENSG00000166888:ENST0000055615      AAGAATCAGTCAACGTGTTGTCAGCCTTCCAGGAGCCTCACCTGCAGATG
ENSG00000166888:ENST0000045407      AAGAATCAGTCAACGTGTTGTCAGCCTTCCAGGAGCCTCACCTGCAGATG
ENSG00000166888:ENST0000053891      AAGAATCAGTCAACGTGTTGTCAGCCTTCCAGGAGCCTCACCTGCAGATG
ENSG00000166888:ENST0000053721      AAGAATCAGTCAACGTGTTGTCAGCCTTCCAGGAGCCTCACCTGCAGATG
ENSG00000166888:ENST0000053520      AAGAATCAGTCAACGTGTTGTCAGCCTTCCAGGAGCCTCACCTGCAGATG
isotig46679                         --------------------------------------------------
isotig12565                         --------------------------------------------------
isotig12566                         --------------------------------------------------

ENSG00000166888:ENST0000030013      CCCCCCAGCCTGGGCCAGATGAGCCTGCCCTTTGACCAGCCTCACCCCCA
ENSG00000166888:ENST0000054387      CCCCCCAGCCTGGGCCAGATGAGCCTGCCCTTTGACCAGCCTCACCCCCA
ENSG00000166888:ENST0000055615      CCCCCCAGCCTGGGCCAGATGAGCCTGCCCTTTGACCAGCCTCACCCCCA
ENSG00000166888:ENST0000045407      CCCCCCAGCCTGGGCCAGATGAGCCTGCCCTTTGACCAGCCTCACCCCCA
ENSG00000166888:ENST0000053891      CCCCCCAGCCTGGGCCAGATGAGCCTGCCCTTTGACCAGCCTCACCCCCA
ENSG00000166888:ENST0000053721      CCCCCCAGCCTGGGCCAGATGAGCCTGCCCTTTGACCAGCCTCACCCCCA
ENSG00000166888:ENST0000053520      CCCCCCAGCCTGGGCCAGATGAGCCTGCCCTTTGACCAGCCTCACCCCCA
isotig46679                         ------------ACCACAATTGAAACAGCAATGAAGTCTCCATATTCTGA
isotig12565                         ------------AATACCATTGACCTGCCGATGTCCCCCCCGCAC-----
isotig12566                         ------------GTTCTGGTTTCGCTGTTAGGGAAAGGGGGACAGTGCAG

ENSG00000166888:ENST0000030013      GGGCCTGCTGCCGTGCCAGCCTCAGGAGCATGCTGTGTCCAGCCCTGACC
ENSG00000166888:ENST0000054387      GGGCCTGCTGCCGTGCCAGCCTCAGGAGCATGCTGTGTCCAGCCCTGACC
ENSG00000166888:ENST0000055615      GGGCCTGCTGCCGTGCCAGCCTCAGGAGCATGCTGTGTCCAGCCCTGACC
ENSG00000166888:ENST0000045407      GGGCCTGCTGCCGTGCCAGCCTCAGGAGCATGCTGTGTCCAGCCCTGACC
ENSG00000166888:ENST0000053891      GGGCCTGCTGCCGTGCCAGCCTCAGGAGCATGCTGTGTCCAGCCCTGACC
ENSG00000166888:ENST0000053721      GGGCCTGCTGCCGTGCCAGCCTCAGGAGCATGCTGTGTCCAGCCCTGACC
ENSG00000166888:ENST0000053520      GGGCCTGCTGCCGTGCCAGCCTCAGGAGCATGCTGTGTCCAGCCCTGACC
isotig46679                         GCGGTAC---------------------AAAGCGACTCTTCAAGGAAGAG
isotig12565                         --------------------------------------------------
isotig12566                         GTCCTTG---------------------------------------GAGG

ENSG00000166888:ENST0000030013      CCCTGCTCTGCTCAGATGTGACCATGGTGGAAGACAGCTGCCTGAGCCAG
ENSG00000166888:ENST0000054387      CCCTGCTCTGCTCAGATGTGACCATGGTGGAAGACAGCTGCCTGAGCCAG
ENSG00000166888:ENST0000055615      CCCTGCTCTGCTCAGATGTGACCATGGTGGAAGACAGCTGCCTGAGCCAG
ENSG00000166888:ENST0000045407      CCCTGCTCTGCTCAGATGTGACCATGGTGGAAGACAGCTGCCTGAGCCAG
ENSG00000166888:ENST0000053891      CCCTGCTCTGCTCAGATGTGACCATGGTGGAAGACAGCTGCCTGAGCCAG
ENSG00000166888:ENST0000053721      CCCTGCTCTGCTCAGATGTGACCATGGTGGAAGACAGCTGCCTGAGCCAG
ENSG00000166888:ENST0000053520      CCCTGCTCTGCTCAGATGTGACCATGGTGGAAGACAGCTGCCTGAGCCAG
isotig46679                         AGCAGATGAAAACGGAGACTGCTCTTTGCCAAAGTCCACAATTCATTTCT
isotig12565                         --------------------------------------------------
isotig12566                         AGAGACAAGGACATGACCGGGTGTCTGGTGGTGAGTCCTGCTATGGAAGA

ENSG00000166888:ENST0000030013      CCAGTGACAGCGTTTCCTCAGGGCACTTGGATTGGTGAAGACATATTCCC
ENSG00000166888:ENST0000054387      CCAGTGACAGCGTTTCCTCAGGGCACTTGGATTGGTGAAGACATATTCCC
ENSG00000166888:ENST0000055615      CCAGTGACAGCGTTTCCTCAGGGCACTTGGATTGGTGAAGACATATTCCC
ENSG00000166888:ENST0000045407      CCAGTGACAGCGTTTCCTCAGGGCACTTGGATTGGTGAAGACATATTCCC
ENSG00000166888:ENST0000053891      CCAGTGACAGCGTTTCCTCAGGGCACTTGGATTGGTGAAGACATATTCCC
ENSG00000166888:ENST0000053721      CCAGTGACAGCGTTTCCTCAGGGCACTTGGATTGGTGAAGACATATTCCC
ENSG00000166888:ENST0000053520      CCAGTGACAGCGTTTCCTCAGGGCACTTGGATTGGTGAAGACATATTCCC
isotig46679                         TCAGCTTTGATACTGGTTTCTAGAAAATGGCACAAATCCGAAGCTTTCCT
isotig12565                         ------TTTAGATTCATTGATGCAGTTTGGAAACGGAGGTGCGCCCTC--
isotig12566                         GCTGTTTAT------------------TGGGTACTTCAG-----------

ENSG00000166888:ENST0000030013      TCCTCTGCTGCCTCCCACTGAACAGGACCTCACTAAGCTTCTCCTGGAGG
ENSG00000166888:ENST0000054387      TCCTCTGCTGCCTCCCACTGAACAGGACCTCACTAAGCTTCTCCTGGAGG
ENSG00000166888:ENST0000055615      TCCTCTGCTGCCTCCCACTGAACAGGACCTCACTAAGCTTCTCCTGGAGG
ENSG00000166888:ENST0000045407      TCCTCTGCTGCCTCCCACTGAACAGGACCTCACTAAGCTTCTCCTGGAGG
ENSG00000166888:ENST0000053891      TCCTCTGCTGCCTCCCACTGAACAGGACCTCACTAAGCTTCTCCTGGAGG
ENSG00000166888:ENST0000053721      TCCTCTGCTGCCTCCCACTGAACAGGACCTCACTAAGCTTCTCCTGGAGG
ENSG00000166888:ENST0000053520      TCCTCTGCTGCCTCCCACTGAACAGGACCTCACTAAGCTTCTCCTGGAGG
isotig46679                         CTCACTA------------------------------------------G
isotig12565                         -------------------------------------------------G
isotig12566                         -------------------------------------------------G

ENSG00000166888:ENST0000030013      GGCAAGGGGAGTCGGGGGGAGGGTCCTTGGGGGCACAGCCCCTCCTGCAG
ENSG00000166888:ENST0000054387      GGCAAGGGGAGTCGGGGGGAGGGTCCTTGGGGGCACAGCCCCTCCTGCAG
ENSG00000166888:ENST0000055615      GGCAAGGGGAGTCGGGGGGAGGGTCCTTGGGGGCACAGCCCCTCCTGCAG
ENSG00000166888:ENST0000045407      GGCAAGGGGAGTCGGGGGGAGGGTCCTTGGGGGCACAGCCCCTCCTGCAG
ENSG00000166888:ENST0000053891      GGCAAGGGGAGTCGGGGGGAGGGTCCTTGGGGGCACAGCCCCTCCTGCAG
ENSG00000166888:ENST0000053721      GGCAAGGGGAGTCGGGGGGAGGGTCCTTGGGGGCACAGCCCCTCCTGCAG
ENSG00000166888:ENST0000053520      GGCAAGGGGAGTCGGGGGGAGGGTCCTTGGGGGCACAGCCCCTCCTGCAG
isotig46679                         GTGACATTCCCCAACTGGGAGTGCTGCTGAAATGCAAACCAAAGCTTCAG
isotig12565                         GCAGGAGGGCAGTTGTCACTCACGTTCATGGATCTGAC------------
isotig12566                         GTGACCGGGATTCAAGAGAAGACCAGAATCAGGCCTCA------------

ENSG00000166888:ENST0000030013      CCCTCCCACTATGGGCAATCTGGGATCTCAATGTCCCACATGGACCTAAG
ENSG00000166888:ENST0000054387      CCCTCCCACTATGGGCAATCTGGGATCTCAATGTCCCACATGGACCTAAG
ENSG00000166888:ENST0000055615      CCCTCCCACTATGGGCAATCTGGGATCTCAATGTCCCACATGGACCTAAG
ENSG00000166888:ENST0000045407      CCCTCCCACTATGGGCAATCTGGGATCTCAATGTCCCACATGGACCTAAG
ENSG00000166888:ENST0000053891      CCCTCCCACTATGGGCAATCTGGGATCTCAATGTCCCACATGGACCTAAG
ENSG00000166888:ENST0000053721      CCCTCCCACTATGGGCAATCTGGGATCTCAATGTCCCACATGGACCTAAG
ENSG00000166888:ENST0000053520      CCCTCCCACTATGGGCAATCTGGGATCTCAATGTCCCACATGGACCTAAG
isotig46679                         ATAAACACGCAGGAAAAGACAGCTTCGAGAAACCTATGTTCGCAATATAA
isotig12565                         ---------TTCGGAGTGCGCTACCTCCCCCATGTGAGGAGC--------
isotig12566                         --------------------------------------------------

ENSG00000166888:ENST0000030013      GGCCAACCCCAGTTGG
ENSG00000166888:ENST0000054387      GGCCAACCCCAGTTGG
ENSG00000166888:ENST0000055615      GGCCAACCCCAGTTGG
ENSG00000166888:ENST0000045407      GGCCAACCCCAGTTGG
ENSG00000166888:ENST0000053891      GGCCAACCCCAGTTGG
ENSG00000166888:ENST0000053721      GGCCAACCCCAGTTGG
ENSG00000166888:ENST0000053520      GGCCAACCCCAGTTGG
isotig46679                         CAGAAGGCTGCTTTGC
isotig12565                         ----------------
isotig12566                         ----------------


""",
        )

    def test4(self):
        aligner = CodonAligner()
        nucleotide_records = SeqIO.index("codonalign/nucl4.fa", "fasta")
        protein_alignment = Align.read("codonalign/pro4.aln", "clustal")
        self.assertEqual(len(protein_alignment.sequences), 10)
        codon_alignments = []
        protein_record = protein_alignment.sequences[0]
        nucleotide_record = nucleotide_records[protein_record.id]
        self.assertEqual(nucleotide_record.id, protein_record.id)
        alignments = aligner.align(protein_record, nucleotide_record)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 847], [0, 2541]]))
        )
        self.assertEqual(
            str(alignment),
            """\
ENSG00000         0 M  S  L  W  G  L  V  S  K  M  P  P  E  K  V  Q  R  L  Y  V  
ENSG00000         0 ATGTCTCTGTGGGGTCTGGTCTCCAAGATGCCCCCAGAAAAAGTGCAGCGGCTCTATGTC

ENSG00000        20 D  F  P  Q  H  L  R  H  L  L  G  D  W  L  E  S  Q  P  W  E  
ENSG00000        60 GACTTTCCCCAACACCTGCGGCATCTTCTGGGTGACTGGCTGGAGAGCCAGCCCTGGGAG

ENSG00000        40 F  L  V  G  S  D  A  F  C  C  N  L  A  S  A  L  L  S  D  T  
ENSG00000       120 TTCCTGGTCGGCTCCGACGCCTTCTGCTGCAACTTGGCTAGTGCCCTACTTTCAGACACT

ENSG00000        60 V  Q  H  L  Q  A  S  V  G  E  Q  G  E  G  S  T  I  L  Q  H  
ENSG00000       180 GTCCAGCACCTTCAGGCCTCGGTGGGAGAGCAGGGGGAGGGGAGCACCATCTTGCAACAC

ENSG00000        80 I  S  T  L  E  S  I  Y  Q  R  D  P  L  K  L  V  A  T  F  R  
ENSG00000       240 ATCAGCACCCTTGAGAGCATATATCAGAGGGACCCCCTGAAGCTGGTGGCCACTTTCAGA

ENSG00000       100 Q  I  L  Q  G  E  K  K  A  V  M  E  Q  F  R  H  L  P  M  P  
ENSG00000       300 CAAATACTTCAAGGAGAGAAAAAAGCTGTTATGGAACAGTTCCGCCACTTGCCAATGCCT

ENSG00000       120 F  H  W  K  Q  E  E  L  K  F  K  T  G  L  R  R  L  Q  H  R  
ENSG00000       360 TTCCACTGGAAGCAGGAAGAACTCAAGTTTAAGACAGGCTTGCGGAGGCTGCAGCACCGA

ENSG00000       140 V  G  E  I  H  L  L  R  E  A  L  Q  K  G  A  E  A  G  Q  V  
ENSG00000       420 GTAGGGGAGATCCACCTTCTCCGAGAAGCCCTGCAGAAGGGGGCTGAGGCTGGCCAAGTG

ENSG00000       160 S  L  H  S  L  I  E  T  P  A  N  G  T  G  P  S  E  A  L  A  
ENSG00000       480 TCTCTGCACAGCTTGATAGAAACTCCTGCTAATGGGACTGGGCCAAGTGAGGCCCTGGCC

ENSG00000       180 M  L  L  Q  E  T  T  G  E  L  E  A  A  K  A  L  V  L  K  R  
ENSG00000       540 ATGCTACTGCAGGAGACCACTGGAGAGCTAGAGGCAGCCAAAGCCCTAGTGCTGAAGAGG

ENSG00000       200 I  Q  I  W  K  R  Q  Q  Q  L  A  G  N  G  A  P  F  E  E  S  
ENSG00000       600 ATCCAGATTTGGAAACGGCAGCAGCAGCTGGCAGGGAATGGCGCACCGTTTGAGGAGAGC

ENSG00000       220 L  A  P  L  Q  E  R  C  E  S  L  V  D  I  Y  S  Q  L  Q  Q  
ENSG00000       660 CTGGCCCCACTCCAGGAGAGGTGTGAAAGCCTGGTGGACATTTATTCCCAGCTACAGCAG

ENSG00000       240 E  V  G  A  A  G  G  E  L  E  P  K  T  R  A  S  L  T  G  R  
ENSG00000       720 GAGGTAGGGGCGGCTGGTGGGGAGCTTGAGCCCAAGACCCGGGCATCGCTGACTGGCCGG

ENSG00000       260 L  D  E  V  L  R  T  L  V  T  S  C  F  L  V  E  K  Q  P  P  
ENSG00000       780 CTGGATGAAGTCCTGAGAACCCTCGTCACCAGTTGCTTCCTGGTGGAGAAGCAGCCCCCC

ENSG00000       280 Q  V  L  K  T  Q  T  K  F  Q  A  G  V  R  F  L  L  G  L  R  
ENSG00000       840 CAGGTACTGAAGACTCAGACCAAGTTCCAGGCTGGAGTTCGATTCCTGTTGGGCTTGAGG

ENSG00000       300 F  L  G  A  P  A  K  P  P  L  V  R  A  D  M  V  T  E  K  Q  
ENSG00000       900 TTCCTGGGGGCCCCAGCCAAGCCTCCGCTGGTCAGGGCCGACATGGTGACAGAGAAGCAG

ENSG00000       320 A  R  E  L  S  V  P  Q  G  P  G  A  G  A  E  S  T  G  E  I  
ENSG00000       960 GCGCGGGAGCTGAGTGTGCCTCAGGGTCCTGGGGCTGGAGCAGAAAGCACTGGAGAAATC

ENSG00000       340 I  N  N  T  V  P  L  E  N  S  I  P  G  N  C  C  S  A  L  F  
ENSG00000      1020 ATCAACAACACTGTGCCCTTGGAGAACAGCATTCCTGGGAACTGCTGCTCTGCCCTGTTC

ENSG00000       360 K  N  L  L  L  K  K  I  K  R  C  E  R  K  G  T  E  S  V  T  
ENSG00000      1080 AAGAACCTGCTTCTCAAGAAGATCAAGCGGTGTGAGCGGAAGGGCACTGAGTCTGTCACA

ENSG00000       380 E  E  K  C  A  V  L  F  S  A  S  F  T  L  G  P  G  K  L  P  
ENSG00000      1140 GAGGAGAAGTGCGCTGTGCTCTTCTCTGCCAGCTTCACACTTGGCCCCGGCAAACTCCCC

ENSG00000       400 I  Q  L  Q  A  L  S  L  P  L  V  V  I  V  H  G  N  Q  D  N  
ENSG00000      1200 ATCCAGCTCCAGGCCCTGTCTCTGCCCCTGGTGGTCATCGTCCATGGCAACCAAGACAAC

ENSG00000       420 N  A  K  A  T  I  L  W  D  N  A  F  S  E  M  D  R  V  P  F  
ENSG00000      1260 AATGCCAAAGCCACTATCCTGTGGGACAATGCCTTCTCTGAGATGGACCGCGTGCCCTTT

ENSG00000       440 V  V  A  E  R  V  P  W  E  K  M  C  E  T  L  N  L  K  F  M  
ENSG00000      1320 GTGGTGGCTGAGCGGGTGCCCTGGGAGAAGATGTGTGAAACTCTGAACCTGAAGTTCATG

ENSG00000       460 A  E  V  G  T  N  R  G  L  L  P  E  H  F  L  F  L  A  Q  K  
ENSG00000      1380 GCTGAGGTGGGGACCAACCGGGGGCTGCTCCCAGAGCACTTCCTCTTCCTGGCCCAGAAG

ENSG00000       480 I  F  N  D  N  S  L  S  M  E  A  F  Q  H  R  S  V  S  W  S  
ENSG00000      1440 ATCTTCAATGACAACAGCCTCAGTATGGAGGCCTTCCAGCACCGTTCTGTGTCCTGGTCG

ENSG00000       500 Q  F  N  K  E  I  L  L  G  R  G  F  T  F  W  Q  W  F  D  G  
ENSG00000      1500 CAGTTCAACAAGGAGATCCTGCTGGGCCGTGGCTTCACCTTTTGGCAGTGGTTTGATGGT

ENSG00000       520 V  L  D  L  T  K  R  C  L  R  S  Y  W  S  D  R  L  I  I  G  
ENSG00000      1560 GTCCTGGACCTCACCAAACGCTGTCTCCGGAGCTACTGGTCTGACCGGCTGATCATTGGC

ENSG00000       540 F  I  S  K  Q  Y  V  T  S  L  L  L  N  E  P  D  G  T  F  L  
ENSG00000      1620 TTCATCAGCAAACAGTACGTTACTAGCCTTCTTCTCAATGAGCCCGACGGAACCTTTCTC

ENSG00000       560 L  R  F  S  D  S  E  I  G  G  I  T  I  A  H  V  I  R  G  Q  
ENSG00000      1680 CTCCGCTTCAGCGACTCAGAGATTGGGGGCATCACCATTGCCCATGTCATCCGGGGCCAG

ENSG00000       580 D  G  S  P  Q  I  E  N  I  Q  P  F  S  A  K  D  L  S  I  R  
ENSG00000      1740 GATGGCTCTCCACAGATAGAGAACATCCAGCCATTCTCTGCCAAAGACCTGTCCATTCGC

ENSG00000       600 S  L  G  D  R  I  R  D  L  A  Q  L  K  N  L  Y  P  K  K  P  
ENSG00000      1800 TCACTGGGGGACCGAATCCGGGATCTTGCTCAGCTCAAAAATCTCTATCCCAAGAAGCCC

ENSG00000       620 K  D  E  A  F  R  S  H  Y  K  P  E  Q  M  G  K  E  G  R  G  
ENSG00000      1860 AAGGATGAGGCTTTCCGGAGCCACTACAAGCCTGAACAGATGGGTAAGGATGGCAGGGGT

ENSG00000       640 Y  V  P  A  T  I  K  M  T  V  E  R  D  Q  P  L  P  T  P  E  
ENSG00000      1920 TATGTCCCAGCTACCATCAAGATGACCGTGGAAAGGGACCAACCACTTCCTACCCCAGAG

ENSG00000       660 L  Q  M  P  T  M  V  P  S  Y  D  L  G  M  A  P  D  S  S  M  
ENSG00000      1980 CTCCAGATGCCTACCATGGTGCCTTCTTATGACCTTGGAATGGCCCCTGATTCCTCCATG

ENSG00000       680 S  M  Q  L  G  P  D  M  V  P  Q  V  Y  P  P  H  S  H  S  I  
ENSG00000      2040 AGCATGCAGCTTGGCCCAGATATGGTGCCCCAGGTGTACCCACCACACTCTCACTCCATC

ENSG00000       700 P  P  Y  Q  G  L  S  P  E  E  S  V  N  V  L  S  A  F  Q  E  
ENSG00000      2100 CCCCCGTATCAAGGCCTCTCCCCAGAAGAATCAGTCAACGTGTTGTCAGCCTTCCAGGAG

ENSG00000       720 P  H  L  Q  M  P  P  S  L  G  Q  M  S  L  P  F  D  Q  P  H  
ENSG00000      2160 CCTCACCTGCAGATGCCCCCCAGCCTGGGCCAGATGAGCCTGCCCTTTGACCAGCCTCAC

ENSG00000       740 P  Q  G  L  L  P  C  Q  P  Q  E  H  A  V  S  S  P  D  P  L  
ENSG00000      2220 CCCCAGGGCCTGCTGCCGTGCCAGCCTCAGGAGCATGCTGTGTCCAGCCCTGACCCCCTG

ENSG00000       760 L  C  S  D  V  T  M  V  E  D  S  C  L  S  Q  P  V  T  A  F  
ENSG00000      2280 CTCTGCTCAGATGTGACCATGGTGGAAGACAGCTGCCTGAGCCAGCCAGTGACAGCGTTT

ENSG00000       780 P  Q  G  T  W  I  G  E  D  I  F  P  P  L  L  P  P  T  E  Q  
ENSG00000      2340 CCTCAGGGCACTTGGATTGGTGAAGACATATTCCCTCCTCTGCTGCCTCCCACTGAACAG

ENSG00000       800 D  L  T  K  L  L  L  E  G  Q  G  E  S  G  G  G  S  L  G  A  
ENSG00000      2400 GACCTCACTAAGCTTCTCCTGGAGGGGCAAGGGGAGTCGGGGGGAGGGTCCTTGGGGGCA

ENSG00000       820 Q  P  L  L  Q  P  S  H  Y  G  Q  S  G  I  S  M  S  H  M  D  
ENSG00000      2460 CAGCCCCTCCTGCAGCCCTCCCACTATGGGCAATCTGGGATCTCAATGTCCCACATGGAC

ENSG00000       840 L  R  A  N  P  S  W    847
ENSG00000      2520 CTAAGGGCCAACCCCAGTTGG 2541
""",
        )
        protein_record = protein_alignment.sequences[1]
        nucleotide_record = nucleotide_records[protein_record.id]
        self.assertEqual(nucleotide_record.id, protein_record.id)
        alignments = aligner.align(protein_record, nucleotide_record)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 847], [0, 2541]]))
        )
        self.assertEqual(
            str(alignment),
            """\
ENSG00000         0 M  S  L  W  G  L  V  S  K  M  P  P  E  K  V  Q  R  L  Y  V  
ENSG00000         0 ATGTCTCTGTGGGGTCTGGTCTCCAAGATGCCCCCAGAAAAAGTGCAGCGGCTCTATGTC

ENSG00000        20 D  F  P  Q  H  L  R  H  L  L  G  D  W  L  E  S  Q  P  W  E  
ENSG00000        60 GACTTTCCCCAACACCTGCGGCATCTTCTGGGTGACTGGCTGGAGAGCCAGCCCTGGGAG

ENSG00000        40 F  L  V  G  S  D  A  F  C  C  N  L  A  S  A  L  L  S  D  T  
ENSG00000       120 TTCCTGGTCGGCTCCGACGCCTTCTGCTGCAACTTGGCTAGTGCCCTACTTTCAGACACT

ENSG00000        60 V  Q  H  L  Q  A  S  V  G  E  Q  G  E  G  S  T  I  L  Q  H  
ENSG00000       180 GTCCAGCACCTTCAGGCCTCGGTGGGAGAGCAGGGGGAGGGGAGCACCATCTTGCAACAC

ENSG00000        80 I  S  T  L  E  S  I  Y  Q  R  D  P  L  K  L  V  A  T  F  R  
ENSG00000       240 ATCAGCACCCTTGAGAGCATATATCAGAGGGACCCCCTGAAGCTGGTGGCCACTTTCAGA

ENSG00000       100 Q  I  L  Q  G  E  K  K  A  V  M  E  Q  F  R  H  L  P  M  P  
ENSG00000       300 CAAATACTTCAAGGAGAGAAAAAAGCTGTTATGGAACAGTTCCGCCACTTGCCAATGCCT

ENSG00000       120 F  H  W  K  Q  E  E  L  K  F  K  T  G  L  R  R  L  Q  H  R  
ENSG00000       360 TTCCACTGGAAGCAGGAAGAACTCAAGTTTAAGACAGGCTTGCGGAGGCTGCAGCACCGA

ENSG00000       140 V  G  E  I  H  L  L  R  E  A  L  Q  K  G  A  E  A  G  Q  V  
ENSG00000       420 GTAGGGGAGATCCACCTTCTCCGAGAAGCCCTGCAGAAGGGGGCTGAGGCTGGCCAAGTG

ENSG00000       160 S  L  H  S  L  I  E  T  P  A  N  G  T  G  P  S  E  A  L  A  
ENSG00000       480 TCTCTGCACAGCTTGATAGAAACTCCTGCTAATGGGACTGGGCCAAGTGAGGCCCTGGCC

ENSG00000       180 M  L  L  Q  E  T  T  G  E  L  E  A  A  K  A  L  V  L  K  R  
ENSG00000       540 ATGCTACTGCAGGAGACCACTGGAGAGCTAGAGGCAGCCAAAGCCCTAGTGCTGAAGAGG

ENSG00000       200 I  Q  I  W  K  R  Q  Q  Q  L  A  G  N  G  A  P  F  E  E  S  
ENSG00000       600 ATCCAGATTTGGAAACGGCAGCAGCAGCTGGCAGGGAATGGCGCACCGTTTGAGGAGAGC

ENSG00000       220 L  A  P  L  Q  E  R  C  E  S  L  V  D  I  Y  S  Q  L  Q  Q  
ENSG00000       660 CTGGCCCCACTCCAGGAGAGGTGTGAAAGCCTGGTGGACATTTATTCCCAGCTACAGCAG

ENSG00000       240 E  V  G  A  A  G  G  E  L  E  P  K  T  R  A  S  L  T  G  R  
ENSG00000       720 GAGGTAGGGGCGGCTGGTGGGGAGCTTGAGCCCAAGACCCGGGCATCGCTGACTGGCCGG

ENSG00000       260 L  D  E  V  L  R  T  L  V  T  S  C  F  L  V  E  K  Q  P  P  
ENSG00000       780 CTGGATGAAGTCCTGAGAACCCTCGTCACCAGTTGCTTCCTGGTGGAGAAGCAGCCCCCC

ENSG00000       280 Q  V  L  K  T  Q  T  K  F  Q  A  G  V  R  F  L  L  G  L  R  
ENSG00000       840 CAGGTACTGAAGACTCAGACCAAGTTCCAGGCTGGAGTTCGATTCCTGTTGGGCTTGAGG

ENSG00000       300 F  L  G  A  P  A  K  P  P  L  V  R  A  D  M  V  T  E  K  Q  
ENSG00000       900 TTCCTGGGGGCCCCAGCCAAGCCTCCGCTGGTCAGGGCCGACATGGTGACAGAGAAGCAG

ENSG00000       320 A  R  E  L  S  V  P  Q  G  P  G  A  G  A  E  S  T  G  E  I  
ENSG00000       960 GCGCGGGAGCTGAGTGTGCCTCAGGGTCCTGGGGCTGGAGCAGAAAGCACTGGAGAAATC

ENSG00000       340 I  N  N  T  V  P  L  E  N  S  I  P  G  N  C  C  S  A  L  F  
ENSG00000      1020 ATCAACAACACTGTGCCCTTGGAGAACAGCATTCCTGGGAACTGCTGCTCTGCCCTGTTC

ENSG00000       360 K  N  L  L  L  K  K  I  K  R  C  E  R  K  G  T  E  S  V  T  
ENSG00000      1080 AAGAACCTGCTTCTCAAGAAGATCAAGCGGTGTGAGCGGAAGGGCACTGAGTCTGTCACA

ENSG00000       380 E  E  K  C  A  V  L  F  S  A  S  F  T  L  G  P  G  K  L  P  
ENSG00000      1140 GAGGAGAAGTGCGCTGTGCTCTTCTCTGCCAGCTTCACACTTGGCCCCGGCAAACTCCCC

ENSG00000       400 I  Q  L  Q  A  L  S  L  P  L  V  V  I  V  H  G  N  Q  D  N  
ENSG00000      1200 ATCCAGCTCCAGGCCCTGTCTCTGCCCCTGGTGGTCATCGTCCATGGCAACCAAGACAAC

ENSG00000       420 N  A  K  A  T  I  L  W  D  N  A  F  S  E  M  D  R  V  P  F  
ENSG00000      1260 AATGCCAAAGCCACTATCCTGTGGGACAATGCCTTCTCTGAGATGGACCGCGTGCCCTTT

ENSG00000       440 V  V  A  E  R  V  P  W  E  K  M  C  E  T  L  N  L  K  F  M  
ENSG00000      1320 GTGGTGGCTGAGCGGGTGCCCTGGGAGAAGATGTGTGAAACTCTGAACCTGAAGTTCATG

ENSG00000       460 A  E  V  G  T  N  R  G  L  L  P  E  H  F  L  F  L  A  Q  K  
ENSG00000      1380 GCTGAGGTGGGGACCAACCGGGGGCTGCTCCCAGAGCACTTCCTCTTCCTGGCCCAGAAG

ENSG00000       480 I  F  N  D  N  S  L  S  M  E  A  F  Q  H  R  S  V  S  W  S  
ENSG00000      1440 ATCTTCAATGACAACAGCCTCAGTATGGAGGCCTTCCAGCACCGTTCTGTGTCCTGGTCG

ENSG00000       500 Q  F  N  K  E  I  L  L  G  R  G  F  T  F  W  Q  W  F  D  G  
ENSG00000      1500 CAGTTCAACAAGGAGATCCTGCTGGGCCGTGGCTTCACCTTTTGGCAGTGGTTTGATGGT

ENSG00000       520 V  L  D  L  T  K  R  C  L  R  S  Y  W  S  D  R  L  I  I  G  
ENSG00000      1560 GTCCTGGACCTCACCAAACGCTGTCTCCGGAGCTACTGGTCTGACCGGCTGATCATTGGC

ENSG00000       540 F  I  S  K  Q  Y  V  T  S  L  L  L  N  E  P  D  G  T  F  L  
ENSG00000      1620 TTCATCAGCAAACAGTACGTTACTAGCCTTCTTCTCAATGAGCCCGACGGAACCTTTCTC

ENSG00000       560 L  R  F  S  D  S  E  I  G  G  I  T  I  A  H  V  I  R  G  Q  
ENSG00000      1680 CTCCGCTTCAGCGACTCAGAGATTGGGGGCATCACCATTGCCCATGTCATCCGGGGCCAG

ENSG00000       580 D  G  S  P  Q  I  E  N  I  Q  P  F  S  A  K  D  L  S  I  R  
ENSG00000      1740 GATGGCTCTCCACAGATAGAGAACATCCAGCCATTCTCTGCCAAAGACCTGTCCATTCGC

ENSG00000       600 S  L  G  D  R  I  R  D  L  A  Q  L  K  N  L  Y  P  K  K  P  
ENSG00000      1800 TCACTGGGGGACCGAATCCGGGATCTTGCTCAGCTCAAAAATCTCTATCCCAAGAAGCCC

ENSG00000       620 K  D  E  A  F  R  S  H  Y  K  P  E  Q  M  G  K  D  G  R  G  
ENSG00000      1860 AAGGATGAGGCTTTCCGGAGCCACTACAAGCCTGAACAGATGGGTAAGGATGGCAGGGGT

ENSG00000       640 Y  V  P  A  T  I  K  M  T  V  E  R  D  Q  P  L  P  T  P  E  
ENSG00000      1920 TATGTCCCAGCTACCATCAAGATGACCGTGGAAAGGGACCAACCACTTCCTACCCCAGAG

ENSG00000       660 L  Q  M  P  T  M  V  P  S  Y  D  L  G  M  A  P  D  S  S  M  
ENSG00000      1980 CTCCAGATGCCTACCATGGTGCCTTCTTATGACCTTGGAATGGCCCCTGATTCCTCCATG

ENSG00000       680 S  M  Q  L  G  P  D  M  V  P  Q  V  Y  P  P  H  S  H  S  I  
ENSG00000      2040 AGCATGCAGCTTGGCCCAGATATGGTGCCCCAGGTGTACCCACCACACTCTCACTCCATC

ENSG00000       700 P  P  Y  Q  G  L  S  P  E  E  S  V  N  V  L  S  A  F  Q  E  
ENSG00000      2100 CCCCCGTATCAAGGCCTCTCCCCAGAAGAATCAGTCAACGTGTTGTCAGCCTTCCAGGAG

ENSG00000       720 P  H  L  Q  M  P  P  S  L  G  Q  M  S  L  P  F  D  Q  P  H  
ENSG00000      2160 CCTCACCTGCAGATGCCCCCCAGCCTGGGCCAGATGAGCCTGCCCTTTGACCAGCCTCAC

ENSG00000       740 P  Q  G  L  L  P  C  Q  P  Q  E  H  A  V  S  S  P  D  P  L  
ENSG00000      2220 CCCCAGGGCCTGCTGCCGTGCCAGCCTCAGGAGCATGCTGTGTCCAGCCCTGACCCCCTG

ENSG00000       760 L  C  S  D  V  T  M  V  E  D  S  C  L  S  Q  P  V  T  A  F  
ENSG00000      2280 CTCTGCTCAGATGTGACCATGGTGGAAGACAGCTGCCTGAGCCAGCCAGTGACAGCGTTT

ENSG00000       780 P  Q  G  T  W  I  G  E  D  I  F  P  P  L  L  P  P  T  E  Q  
ENSG00000      2340 CCTCAGGGCACTTGGATTGGTGAAGACATATTCCCTCCTCTGCTGCCTCCCACTGAACAG

ENSG00000       800 D  L  T  K  L  L  L  E  G  Q  G  E  S  G  G  G  S  L  G  A  
ENSG00000      2400 GACCTCACTAAGCTTCTCCTGGAGGGGCAAGGGGAGTCGGGGGGAGGGTCCTTGGGGGCA

ENSG00000       820 Q  P  L  L  Q  P  S  H  Y  G  Q  S  G  I  S  M  S  H  M  D  
ENSG00000      2460 CAGCCCCTCCTGCAGCCCTCCCACTATGGGCAATCTGGGATCTCAATGTCCCACATGGAC

ENSG00000       840 L  R  A  N  P  S  W    847
ENSG00000      2520 CTAAGGGCCAACCCCAGTTGG 2541
""",
        )
        protein_record = protein_alignment.sequences[2]
        nucleotide_record = nucleotide_records[protein_record.id]
        self.assertEqual(nucleotide_record.id, protein_record.id)
        alignments = aligner.align(protein_record, nucleotide_record)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 847], [0, 2541]]))
        )
        self.assertEqual(
            str(alignment),
            """\
ENSG00000         0 M  S  L  W  G  L  V  S  K  M  P  P  E  K  V  Q  R  L  Y  V  
ENSG00000         0 ATGTCTCTGTGGGGTCTGGTCTCCAAGATGCCCCCAGAAAAAGTGCAGCGGCTCTATGTC

ENSG00000        20 D  F  P  Q  H  L  R  H  L  L  G  D  W  L  E  S  Q  P  W  E  
ENSG00000        60 GACTTTCCCCAACACCTGCGGCATCTTCTGGGTGACTGGCTGGAGAGCCAGCCCTGGGAG

ENSG00000        40 F  L  V  G  S  D  A  F  C  C  N  L  A  S  A  L  L  S  D  T  
ENSG00000       120 TTCCTGGTCGGCTCCGACGCCTTCTGCTGCAACTTGGCTAGTGCCCTACTTTCAGACACT

ENSG00000        60 V  Q  H  L  Q  A  S  V  G  E  Q  G  E  G  S  T  I  L  Q  H  
ENSG00000       180 GTCCAGCACCTTCAGGCCTCGGTGGGAGAGCAGGGGGAGGGGAGCACCATCTTGCAACAC

ENSG00000        80 I  S  T  L  E  S  I  Y  Q  R  D  P  L  K  L  V  A  T  F  R  
ENSG00000       240 ATCAGCACCCTTGAGAGCATATATCAGAGGGACCCCCTGAAGCTGGTGGCCACTTTCAGA

ENSG00000       100 Q  I  L  Q  G  E  K  K  A  V  M  E  Q  F  R  H  L  P  M  P  
ENSG00000       300 CAAATACTTCAAGGAGAGAAAAAAGCTGTTATGGAACAGTTCCGCCACTTGCCAATGCCT

ENSG00000       120 F  H  W  K  Q  E  E  L  K  F  K  T  G  L  R  R  L  Q  H  R  
ENSG00000       360 TTCCACTGGAAGCAGGAAGAACTCAAGTTTAAGACAGGCTTGCGGAGGCTGCAGCACCGA

ENSG00000       140 V  G  E  I  H  L  L  R  E  A  L  Q  K  G  A  E  A  G  Q  V  
ENSG00000       420 GTAGGGGAGATCCACCTTCTCCGAGAAGCCCTGCAGAAGGGGGCTGAGGCTGGCCAAGTG

ENSG00000       160 S  L  H  S  L  I  E  T  P  A  N  G  T  G  P  S  E  A  L  A  
ENSG00000       480 TCTCTGCACAGCTTGATAGAAACTCCTGCTAATGGGACTGGGCCAAGTGAGGCCCTGGCC

ENSG00000       180 M  L  L  Q  E  T  T  G  E  L  E  A  A  K  A  L  V  L  K  R  
ENSG00000       540 ATGCTACTGCAGGAGACCACTGGAGAGCTAGAGGCAGCCAAAGCCCTAGTGCTGAAGAGG

ENSG00000       200 I  Q  I  W  K  R  Q  Q  Q  L  A  G  N  G  A  P  F  E  E  S  
ENSG00000       600 ATCCAGATTTGGAAACGGCAGCAGCAGCTGGCAGGGAATGGCGCACCGTTTGAGGAGAGC

ENSG00000       220 L  A  P  L  Q  E  R  C  E  S  L  V  D  I  Y  S  Q  L  Q  Q  
ENSG00000       660 CTGGCCCCACTCCAGGAGAGGTGTGAAAGCCTGGTGGACATTTATTCCCAGCTACAGCAG

ENSG00000       240 E  V  G  A  A  G  G  E  L  E  P  K  T  R  A  S  L  T  G  R  
ENSG00000       720 GAGGTAGGGGCGGCTGGTGGGGAGCTTGAGCCCAAGACCCGGGCATCGCTGACTGGCCGG

ENSG00000       260 L  D  E  V  L  R  T  L  V  T  S  C  F  L  V  E  K  Q  P  P  
ENSG00000       780 CTGGATGAAGTCCTGAGAACCCTCGTCACCAGTTGCTTCCTGGTGGAGAAGCAGCCCCCC

ENSG00000       280 Q  V  L  K  T  Q  T  K  F  Q  A  G  V  R  F  L  L  G  L  R  
ENSG00000       840 CAGGTACTGAAGACTCAGACCAAGTTCCAGGCTGGAGTTCGATTCCTGTTGGGCTTGAGG

ENSG00000       300 F  L  G  A  P  A  K  P  P  L  V  R  A  D  M  V  T  E  K  Q  
ENSG00000       900 TTCCTGGGGGCCCCAGCCAAGCCTCCGCTGGTCAGGGCCGACATGGTGACAGAGAAGCAG

ENSG00000       320 A  R  E  L  S  V  P  Q  G  P  G  A  G  A  E  S  T  G  E  I  
ENSG00000       960 GCGCGGGAGCTGAGTGTGCCTCAGGGTCCTGGGGCTGGAGCAGAAAGCACTGGAGAAATC

ENSG00000       340 I  N  N  T  V  P  L  E  N  S  I  P  G  N  C  C  S  A  L  F  
ENSG00000      1020 ATCAACAACACTGTGCCCTTGGAGAACAGCATTCCTGGGAACTGCTGCTCTGCCCTGTTC

ENSG00000       360 K  N  L  L  L  K  K  I  K  R  C  E  R  K  G  T  E  S  V  T  
ENSG00000      1080 AAGAACCTGCTTCTCAAGAAGATCAAGCGGTGTGAGCGGAAGGGCACTGAGTCTGTCACA

ENSG00000       380 E  E  K  C  A  V  L  F  S  A  S  F  T  L  G  P  G  K  L  P  
ENSG00000      1140 GAGGAGAAGTGCGCTGTGCTCTTCTCTGCCAGCTTCACACTTGGCCCCGGCAAACTCCCC

ENSG00000       400 I  Q  L  Q  A  L  S  L  P  L  V  V  I  V  H  G  N  Q  D  N  
ENSG00000      1200 ATCCAGCTCCAGGCCCTGTCTCTGCCCCTGGTGGTCATCGTCCATGGCAACCAAGACAAC

ENSG00000       420 N  A  K  A  T  I  L  W  D  N  A  F  S  E  M  D  R  V  P  F  
ENSG00000      1260 AATGCCAAAGCCACTATCCTGTGGGACAATGCCTTCTCTGAGATGGACCGCGTGCCCTTT

ENSG00000       440 V  V  A  E  R  V  P  W  E  K  M  C  E  T  L  N  L  K  F  M  
ENSG00000      1320 GTGGTGGCTGAGCGGGTGCCCTGGGAGAAGATGTGTGAAACTCTGAACCTGAAGTTCATG

ENSG00000       460 A  E  V  G  T  N  R  G  L  L  P  E  H  F  L  F  L  A  Q  K  
ENSG00000      1380 GCTGAGGTGGGGACCAACCGGGGGCTGCTCCCAGAGCACTTCCTCTTCCTGGCCCAGAAG

ENSG00000       480 I  F  N  D  N  S  L  S  M  E  A  F  Q  H  R  S  V  S  W  S  
ENSG00000      1440 ATCTTCAATGACAACAGCCTCAGTATGGAGGCCTTCCAGCACCGTTCTGTGTCCTGGTCG

ENSG00000       500 Q  F  N  K  E  I  L  L  G  R  G  F  T  F  W  Q  W  F  D  G  
ENSG00000      1500 CAGTTCAACAAGGAGATCCTGCTGGGCCGTGGCTTCACCTTTTGGCAGTGGTTTGATGGT

ENSG00000       520 V  L  D  L  T  K  R  C  L  R  S  Y  W  S  D  R  L  I  I  G  
ENSG00000      1560 GTCCTGGACCTCACCAAACGCTGTCTCCGGAGCTACTGGTCTGACCGGCTGATCATTGGC

ENSG00000       540 F  I  S  K  Q  Y  V  T  S  L  L  L  N  E  P  D  G  T  F  L  
ENSG00000      1620 TTCATCAGCAAACAGTACGTTACTAGCCTTCTTCTCAATGAGCCCGACGGAACCTTTCTC

ENSG00000       560 L  R  F  S  D  S  E  I  G  G  I  T  I  A  H  V  I  R  G  Q  
ENSG00000      1680 CTCCGCTTCAGCGACTCAGAGATTGGGGGCATCACCATTGCCCATGTCATCCGGGGCCAG

ENSG00000       580 D  G  S  P  Q  I  E  N  I  Q  P  F  S  A  K  D  L  S  I  R  
ENSG00000      1740 GATGGCTCTCCACAGATAGAGAACATCCAGCCATTCTCTGCCAAAGACCTGTCCATTCGC

ENSG00000       600 S  L  G  D  R  I  R  D  L  A  Q  L  K  N  L  Y  P  K  K  P  
ENSG00000      1800 TCACTGGGGGACCGAATCCGGGATCTTGCTCAGCTCAAAAATCTCTATCCCAAGAAGCCC

ENSG00000       620 K  D  E  A  F  R  S  H  Y  K  P  E  Q  M  G  K  D  G  R  G  
ENSG00000      1860 AAGGATGAGGCTTTCCGGAGCCACTACAAGCCTGAACAGATGGGTAAGGATGGCAGGGGT

ENSG00000       640 Y  V  P  A  T  I  K  M  T  V  E  R  D  Q  P  L  P  T  P  E  
ENSG00000      1920 TATGTCCCAGCTACCATCAAGATGACCGTGGAAAGGGACCAACCACTTCCTACCCCAGAG

ENSG00000       660 L  Q  M  P  T  M  V  P  S  Y  D  L  G  M  A  P  D  S  S  M  
ENSG00000      1980 CTCCAGATGCCTACCATGGTGCCTTCTTATGACCTTGGAATGGCCCCTGATTCCTCCATG

ENSG00000       680 S  M  Q  L  G  P  D  M  V  P  Q  V  Y  P  P  H  S  H  S  I  
ENSG00000      2040 AGCATGCAGCTTGGCCCAGATATGGTGCCCCAGGTGTACCCACCACACTCTCACTCCATC

ENSG00000       700 P  P  Y  Q  G  L  S  P  E  E  S  V  N  V  L  S  A  F  Q  E  
ENSG00000      2100 CCCCCGTATCAAGGCCTCTCCCCAGAAGAATCAGTCAACGTGTTGTCAGCCTTCCAGGAG

ENSG00000       720 P  H  L  Q  M  P  P  S  L  G  Q  M  S  L  P  F  D  Q  P  H  
ENSG00000      2160 CCTCACCTGCAGATGCCCCCCAGCCTGGGCCAGATGAGCCTGCCCTTTGACCAGCCTCAC

ENSG00000       740 P  Q  G  L  L  P  C  Q  P  Q  E  H  A  V  S  S  P  D  P  L  
ENSG00000      2220 CCCCAGGGCCTGCTGCCGTGCCAGCCTCAGGAGCATGCTGTGTCCAGCCCTGACCCCCTG

ENSG00000       760 L  C  S  D  V  T  M  V  E  D  S  C  L  S  Q  P  V  T  A  F  
ENSG00000      2280 CTCTGCTCAGATGTGACCATGGTGGAAGACAGCTGCCTGAGCCAGCCAGTGACAGCGTTT

ENSG00000       780 P  Q  G  T  W  I  G  E  D  I  F  P  P  L  L  P  P  T  E  Q  
ENSG00000      2340 CCTCAGGGCACTTGGATTGGTGAAGACATATTCCCTCCTCTGCTGCCTCCCACTGAACAG

ENSG00000       800 D  L  T  K  L  L  L  E  G  Q  G  E  S  G  G  G  S  L  G  A  
ENSG00000      2400 GACCTCACTAAGCTTCTCCTGGAGGGGCAAGGGGAGTCGGGGGGAGGGTCCTTGGGGGCA

ENSG00000       820 Q  P  L  L  Q  P  S  H  Y  G  Q  S  G  I  S  M  S  H  M  D  
ENSG00000      2460 CAGCCCCTCCTGCAGCCCTCCCACTATGGGCAATCTGGGATCTCAATGTCCCACATGGAC

ENSG00000       840 L  R  A  N  P  S  W    847
ENSG00000      2520 CTAAGGGCCAACCCCAGTTGG 2541
""",
        )
        protein_record = protein_alignment.sequences[3]
        nucleotide_record = nucleotide_records[protein_record.id]
        self.assertEqual(nucleotide_record.id, protein_record.id)
        alignments = aligner.align(protein_record, nucleotide_record)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 847], [0, 2541]]))
        )
        self.assertEqual(
            str(alignment),
            """\
ENSG00000         0 M  S  L  W  G  L  V  S  K  M  P  P  E  K  V  Q  R  L  Y  V  
ENSG00000         0 ATGTCTCTGTGGGGTCTGGTCTCCAAGATGCCCCCAGAAAAAGTGCAGCGGCTCTATGTC

ENSG00000        20 D  F  P  Q  H  L  R  H  L  L  G  D  W  L  E  S  Q  P  W  E  
ENSG00000        60 GACTTTCCCCAACACCTGCGGCATCTTCTGGGTGACTGGCTGGAGAGCCAGCCCTGGGAG

ENSG00000        40 F  L  V  G  S  D  A  F  C  C  N  L  A  S  A  L  L  S  D  T  
ENSG00000       120 TTCCTGGTCGGCTCCGACGCCTTCTGCTGCAACTTGGCTAGTGCCCTACTTTCAGACACT

ENSG00000        60 V  Q  H  L  Q  A  S  V  G  E  Q  G  E  G  S  T  I  L  Q  H  
ENSG00000       180 GTCCAGCACCTTCAGGCCTCGGTGGGAGAGCAGGGGGAGGGGAGCACCATCTTGCAACAC

ENSG00000        80 I  S  T  L  E  S  I  Y  Q  R  D  P  L  K  L  V  A  T  F  R  
ENSG00000       240 ATCAGCACCCTTGAGAGCATATATCAGAGGGACCCCCTGAAGCTGGTGGCCACTTTCAGA

ENSG00000       100 Q  I  L  Q  G  E  K  K  A  V  M  E  Q  F  R  H  L  P  M  P  
ENSG00000       300 CAAATACTTCAAGGAGAGAAAAAAGCTGTTATGGAACAGTTCCGCCACTTGCCAATGCCT

ENSG00000       120 F  H  W  K  Q  E  E  L  K  F  K  T  G  L  R  R  L  Q  H  R  
ENSG00000       360 TTCCACTGGAAGCAGGAAGAACTCAAGTTTAAGACAGGCTTGCGGAGGCTGCAGCACCGA

ENSG00000       140 V  G  E  I  H  L  L  R  E  A  L  Q  K  G  A  E  A  G  Q  V  
ENSG00000       420 GTAGGGGAGATCCACCTTCTCCGAGAAGCCCTGCAGAAGGGGGCTGAGGCTGGCCAAGTG

ENSG00000       160 S  L  H  S  L  I  E  T  P  A  N  G  T  G  P  S  E  A  L  A  
ENSG00000       480 TCTCTGCACAGCTTGATAGAAACTCCTGCTAATGGGACTGGGCCAAGTGAGGCCCTGGCC

ENSG00000       180 M  L  L  Q  E  T  T  G  E  L  E  A  A  K  A  L  V  L  K  R  
ENSG00000       540 ATGCTACTGCAGGAGACCACTGGAGAGCTAGAGGCAGCCAAAGCCCTAGTGCTGAAGAGG

ENSG00000       200 I  Q  I  W  K  R  Q  Q  Q  L  A  G  N  G  A  P  F  E  E  S  
ENSG00000       600 ATCCAGATTTGGAAACGGCAGCAGCAGCTGGCAGGGAATGGCGCACCGTTTGAGGAGAGC

ENSG00000       220 L  A  P  L  Q  E  R  C  E  S  L  V  D  I  Y  S  Q  L  Q  Q  
ENSG00000       660 CTGGCCCCACTCCAGGAGAGGTGTGAAAGCCTGGTGGACATTTATTCCCAGCTACAGCAG

ENSG00000       240 E  V  G  A  A  G  G  E  L  E  P  K  T  R  A  S  L  T  G  R  
ENSG00000       720 GAGGTAGGGGCGGCTGGTGGGGAGCTTGAGCCCAAGACCCGGGCATCGCTGACTGGCCGG

ENSG00000       260 L  D  E  V  L  R  T  L  V  T  S  C  F  L  V  E  K  Q  P  P  
ENSG00000       780 CTGGATGAAGTCCTGAGAACCCTCGTCACCAGTTGCTTCCTGGTGGAGAAGCAGCCCCCC

ENSG00000       280 Q  V  L  K  T  Q  T  K  F  Q  A  G  V  R  F  L  L  G  L  R  
ENSG00000       840 CAGGTACTGAAGACTCAGACCAAGTTCCAGGCTGGAGTTCGATTCCTGTTGGGCTTGAGG

ENSG00000       300 F  L  G  A  P  A  K  P  P  L  V  R  A  D  M  V  T  E  K  Q  
ENSG00000       900 TTCCTGGGGGCCCCAGCCAAGCCTCCGCTGGTCAGGGCCGACATGGTGACAGAGAAGCAG

ENSG00000       320 A  R  E  L  S  V  P  Q  G  P  G  A  G  A  E  S  T  G  E  I  
ENSG00000       960 GCGCGGGAGCTGAGTGTGCCTCAGGGTCCTGGGGCTGGAGCAGAAAGCACTGGAGAAATC

ENSG00000       340 I  N  N  T  V  P  L  E  N  S  I  P  G  N  C  C  S  A  L  F  
ENSG00000      1020 ATCAACAACACTGTGCCCTTGGAGAACAGCATTCCTGGGAACTGCTGCTCTGCCCTGTTC

ENSG00000       360 K  N  L  L  L  K  K  I  K  R  C  E  R  K  G  T  E  S  V  T  
ENSG00000      1080 AAGAACCTGCTTCTCAAGAAGATCAAGCGGTGTGAGCGGAAGGGCACTGAGTCTGTCACA

ENSG00000       380 E  E  K  C  A  V  L  F  S  A  S  F  T  L  G  P  G  K  L  P  
ENSG00000      1140 GAGGAGAAGTGCGCTGTGCTCTTCTCTGCCAGCTTCACACTTGGCCCCGGCAAACTCCCC

ENSG00000       400 I  Q  L  Q  A  L  S  L  P  L  V  V  I  V  H  G  N  Q  D  N  
ENSG00000      1200 ATCCAGCTCCAGGCCCTGTCTCTGCCCCTGGTGGTCATCGTCCATGGCAACCAAGACAAC

ENSG00000       420 N  A  K  A  T  I  L  W  D  N  A  F  S  E  M  D  R  V  P  F  
ENSG00000      1260 AATGCCAAAGCCACTATCCTGTGGGACAATGCCTTCTCTGAGATGGACCGCGTGCCCTTT

ENSG00000       440 V  V  A  E  R  V  P  W  E  K  M  C  E  T  L  N  L  K  F  M  
ENSG00000      1320 GTGGTGGCTGAGCGGGTGCCCTGGGAGAAGATGTGTGAAACTCTGAACCTGAAGTTCATG

ENSG00000       460 A  E  V  G  T  N  R  G  L  L  P  E  H  F  L  F  L  A  Q  K  
ENSG00000      1380 GCTGAGGTGGGGACCAACCGGGGGCTGCTCCCAGAGCACTTCCTCTTCCTGGCCCAGAAG

ENSG00000       480 I  F  N  D  N  S  L  S  M  E  A  F  Q  H  R  S  V  S  W  S  
ENSG00000      1440 ATCTTCAATGACAACAGCCTCAGTATGGAGGCCTTCCAGCACCGTTCTGTGTCCTGGTCG

ENSG00000       500 Q  F  N  K  E  I  L  L  G  R  G  F  T  F  W  Q  W  F  D  G  
ENSG00000      1500 CAGTTCAACAAGGAGATCCTGCTGGGCCGTGGCTTCACCTTTTGGCAGTGGTTTGATGGT

ENSG00000       520 V  L  D  L  T  K  R  C  L  R  S  Y  W  S  D  R  L  I  I  G  
ENSG00000      1560 GTCCTGGACCTCACCAAACGCTGTCTCCGGAGCTACTGGTCTGACCGGCTGATCATTGGC

ENSG00000       540 F  I  S  K  Q  Y  V  T  S  L  L  L  N  E  P  D  G  T  F  L  
ENSG00000      1620 TTCATCAGCAAACAGTACGTTACTAGCCTTCTTCTCAATGAGCCCGACGGAACCTTTCTC

ENSG00000       560 L  R  F  S  D  S  E  I  G  G  I  T  I  A  H  V  I  R  G  Q  
ENSG00000      1680 CTCCGCTTCAGCGACTCAGAGATTGGGGGCATCACCATTGCCCATGTCATCCGGGGCCAG

ENSG00000       580 D  G  S  P  Q  I  E  N  I  Q  P  F  S  A  K  D  L  S  I  R  
ENSG00000      1740 GATGGCTCTCCACAGATAGAGAACATCCAGCCATTCTCTGCCAAAGACCTGTCCATTCGC

ENSG00000       600 S  L  G  D  R  I  R  D  L  A  Q  L  K  N  L  Y  P  K  K  P  
ENSG00000      1800 TCACTGGGGGACCGAATCCGGGATCTTGCTCAGCTCAAAAATCTCTATCCCAAGAAGCCC

ENSG00000       620 K  D  E  A  F  R  S  H  Y  K  P  E  Q  M  G  K  D  G  R  G  
ENSG00000      1860 AAGGATGAGGCTTTCCGGAGCCACTACAAGCCTGAACAGATGGGTAAGGATGGCAGGGGT

ENSG00000       640 Y  V  P  A  T  I  K  M  T  V  E  R  D  Q  P  L  P  T  P  E  
ENSG00000      1920 TATGTCCCAGCTACCATCAAGATGACCGTGGAAAGGGACCAACCACTTCCTACCCCAGAG

ENSG00000       660 L  Q  M  P  T  M  V  P  S  Y  D  L  G  M  A  P  D  S  S  M  
ENSG00000      1980 CTCCAGATGCCTACCATGGTGCCTTCTTATGACCTTGGAATGGCCCCTGATTCCTCCATG

ENSG00000       680 S  M  Q  L  G  P  D  M  V  P  Q  V  Y  P  P  H  S  H  S  I  
ENSG00000      2040 AGCATGCAGCTTGGCCCAGATATGGTGCCCCAGGTGTACCCACCACACTCTCACTCCATC

ENSG00000       700 P  P  Y  Q  G  L  S  P  E  E  S  V  N  V  L  S  A  F  Q  E  
ENSG00000      2100 CCCCCGTATCAAGGCCTCTCCCCAGAAGAATCAGTCAACGTGTTGTCAGCCTTCCAGGAG

ENSG00000       720 P  H  L  Q  M  P  P  S  L  G  Q  M  S  L  P  F  D  Q  P  H  
ENSG00000      2160 CCTCACCTGCAGATGCCCCCCAGCCTGGGCCAGATGAGCCTGCCCTTTGACCAGCCTCAC

ENSG00000       740 P  Q  G  L  L  P  C  Q  P  Q  E  H  A  V  S  S  P  D  P  L  
ENSG00000      2220 CCCCAGGGCCTGCTGCCGTGCCAGCCTCAGGAGCATGCTGTGTCCAGCCCTGACCCCCTG

ENSG00000       760 L  C  S  D  V  T  M  V  E  D  S  C  L  S  Q  P  V  T  A  F  
ENSG00000      2280 CTCTGCTCAGATGTGACCATGGTGGAAGACAGCTGCCTGAGCCAGCCAGTGACAGCGTTT

ENSG00000       780 P  Q  G  T  W  I  G  E  D  I  F  P  P  L  L  P  P  T  E  Q  
ENSG00000      2340 CCTCAGGGCACTTGGATTGGTGAAGACATATTCCCTCCTCTGCTGCCTCCCACTGAACAG

ENSG00000       800 D  L  T  K  L  L  L  E  G  Q  G  E  S  G  G  G  S  L  G  A  
ENSG00000      2400 GACCTCACTAAGCTTCTCCTGGAGGGGCAAGGGGAGTCGGGGGGAGGGTCCTTGGGGGCA

ENSG00000       820 Q  P  L  L  Q  P  S  H  Y  G  Q  S  G  I  S  M  S  H  M  D  
ENSG00000      2460 CAGCCCCTCCTGCAGCCCTCCCACTATGGGCAATCTGGGATCTCAATGTCCCACATGGAC

ENSG00000       840 L  R  A  N  P  S  W    847
ENSG00000      2520 CTAAGGGCCAACCCCAGTTGG 2541
""",
        )
        protein_record = protein_alignment.sequences[4]
        nucleotide_record = nucleotide_records[protein_record.id]
        self.assertEqual(nucleotide_record.id, protein_record.id)
        alignments = aligner.align(protein_record, nucleotide_record)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 737], [0, 2211]]))
        )
        self.assertEqual(
            str(alignment),
            """\
ENSG00000         0 M  E  Q  F  R  H  L  P  M  P  F  H  W  K  Q  E  E  L  K  F  
ENSG00000         0 ATGGAACAGTTCCGCCACTTGCCAATGCCTTTCCACTGGAAGCAGGAAGAACTCAAGTTT

ENSG00000        20 K  T  G  L  R  R  L  Q  H  R  V  G  E  I  H  L  L  R  E  A  
ENSG00000        60 AAGACAGGCTTGCGGAGGCTGCAGCACCGAGTAGGGGAGATCCACCTTCTCCGAGAAGCC

ENSG00000        40 L  Q  K  G  A  E  A  G  Q  V  S  L  H  S  L  I  E  T  P  A  
ENSG00000       120 CTGCAGAAGGGGGCTGAGGCTGGCCAAGTGTCTCTGCACAGCTTGATAGAAACTCCTGCT

ENSG00000        60 N  G  T  G  P  S  E  A  L  A  M  L  L  Q  E  T  T  G  E  L  
ENSG00000       180 AATGGGACTGGGCCAAGTGAGGCCCTGGCCATGCTACTGCAGGAGACCACTGGAGAGCTA

ENSG00000        80 E  A  A  K  A  L  V  L  K  R  I  Q  I  W  K  R  Q  Q  Q  L  
ENSG00000       240 GAGGCAGCCAAAGCCCTAGTGCTGAAGAGGATCCAGATTTGGAAACGGCAGCAGCAGCTG

ENSG00000       100 A  G  N  G  A  P  F  E  E  S  L  A  P  L  Q  E  R  C  E  S  
ENSG00000       300 GCAGGGAATGGCGCACCGTTTGAGGAGAGCCTGGCCCCACTCCAGGAGAGGTGTGAAAGC

ENSG00000       120 L  V  D  I  Y  S  Q  L  Q  Q  E  V  G  A  A  G  G  E  L  E  
ENSG00000       360 CTGGTGGACATTTATTCCCAGCTACAGCAGGAGGTAGGGGCGGCTGGTGGGGAGCTTGAG

ENSG00000       140 P  K  T  R  A  S  L  T  G  R  L  D  E  V  L  R  T  L  V  T  
ENSG00000       420 CCCAAGACCCGGGCATCGCTGACTGGCCGGCTGGATGAAGTCCTGAGAACCCTCGTCACC

ENSG00000       160 S  C  F  L  V  E  K  Q  P  P  Q  V  L  K  T  Q  T  K  F  Q  
ENSG00000       480 AGTTGCTTCCTGGTGGAGAAGCAGCCCCCCCAGGTACTGAAGACTCAGACCAAGTTCCAG

ENSG00000       180 A  G  V  R  F  L  L  G  L  R  F  L  G  A  P  A  K  P  P  L  
ENSG00000       540 GCTGGAGTTCGATTCCTGTTGGGCTTGAGGTTCCTGGGGGCCCCAGCCAAGCCTCCGCTG

ENSG00000       200 V  R  A  D  M  V  T  E  K  Q  A  R  E  L  S  V  P  Q  G  P  
ENSG00000       600 GTCAGGGCCGACATGGTGACAGAGAAGCAGGCGCGGGAGCTGAGTGTGCCTCAGGGTCCT

ENSG00000       220 G  A  G  A  E  S  T  G  E  I  I  N  N  T  V  P  L  E  N  S  
ENSG00000       660 GGGGCTGGAGCAGAAAGCACTGGAGAAATCATCAACAACACTGTGCCCTTGGAGAACAGC

ENSG00000       240 I  P  G  N  C  C  S  A  L  F  K  N  L  L  L  K  K  I  K  R  
ENSG00000       720 ATTCCTGGGAACTGCTGCTCTGCCCTGTTCAAGAACCTGCTTCTCAAGAAGATCAAGCGG

ENSG00000       260 C  E  R  K  G  T  E  S  V  T  E  E  K  C  A  V  L  F  S  A  
ENSG00000       780 TGTGAGCGGAAGGGCACTGAGTCTGTCACAGAGGAGAAGTGCGCTGTGCTCTTCTCTGCC

ENSG00000       280 S  F  T  L  G  P  G  K  L  P  I  Q  L  Q  A  L  S  L  P  L  
ENSG00000       840 AGCTTCACACTTGGCCCCGGCAAACTCCCCATCCAGCTCCAGGCCCTGTCTCTGCCCCTG

ENSG00000       300 V  V  I  V  H  G  N  Q  D  N  N  A  K  A  T  I  L  W  D  N  
ENSG00000       900 GTGGTCATCGTCCATGGCAACCAAGACAACAATGCCAAAGCCACTATCCTGTGGGACAAT

ENSG00000       320 A  F  S  E  M  D  R  V  P  F  V  V  A  E  R  V  P  W  E  K  
ENSG00000       960 GCCTTCTCTGAGATGGACCGCGTGCCCTTTGTGGTGGCTGAGCGGGTGCCCTGGGAGAAG

ENSG00000       340 M  C  E  T  L  N  L  K  F  M  A  E  V  G  T  N  R  G  L  L  
ENSG00000      1020 ATGTGTGAAACTCTGAACCTGAAGTTCATGGCTGAGGTGGGGACCAACCGGGGGCTGCTC

ENSG00000       360 P  E  H  F  L  F  L  A  Q  K  I  F  N  D  N  S  L  S  M  E  
ENSG00000      1080 CCAGAGCACTTCCTCTTCCTGGCCCAGAAGATCTTCAATGACAACAGCCTCAGTATGGAG

ENSG00000       380 A  F  Q  H  R  S  V  S  W  S  Q  F  N  K  E  I  L  L  G  R  
ENSG00000      1140 GCCTTCCAGCACCGTTCTGTGTCCTGGTCGCAGTTCAACAAGGAGATCCTGCTGGGCCGT

ENSG00000       400 G  F  T  F  W  Q  W  F  D  G  V  L  D  L  T  K  R  C  L  R  
ENSG00000      1200 GGCTTCACCTTTTGGCAGTGGTTTGATGGTGTCCTGGACCTCACCAAACGCTGTCTCCGG

ENSG00000       420 S  Y  W  S  D  R  L  I  I  G  F  I  S  K  Q  Y  V  T  S  L  
ENSG00000      1260 AGCTACTGGTCTGACCGGCTGATCATTGGCTTCATCAGCAAACAGTACGTTACTAGCCTT

ENSG00000       440 L  L  N  E  P  D  G  T  F  L  L  R  F  S  D  S  E  I  G  G  
ENSG00000      1320 CTTCTCAATGAGCCCGACGGAACCTTTCTCCTCCGCTTCAGCGACTCAGAGATTGGGGGC

ENSG00000       460 I  T  I  A  H  V  I  R  G  Q  D  G  S  P  Q  I  E  N  I  Q  
ENSG00000      1380 ATCACCATTGCCCATGTCATCCGGGGCCAGGATGGCTCTCCACAGATAGAGAACATCCAG

ENSG00000       480 P  F  S  A  K  D  L  S  I  R  S  L  G  D  R  I  R  D  L  A  
ENSG00000      1440 CCATTCTCTGCCAAAGACCTGTCCATTCGCTCACTGGGGGACCGAATCCGGGATCTTGCT

ENSG00000       500 Q  L  K  N  L  Y  P  K  K  P  K  D  E  A  F  R  S  H  Y  K  
ENSG00000      1500 CAGCTCAAAAATCTCTATCCCAAGAAGCCCAAGGATGAGGCTTTCCGGAGCCACTACAAG

ENSG00000       520 P  E  Q  M  G  K  D  G  R  G  Y  V  P  A  T  I  K  M  T  V  
ENSG00000      1560 CCTGAACAGATGGGTAAGGATGGCAGGGGTTATGTCCCAGCTACCATCAAGATGACCGTG

ENSG00000       540 E  R  D  Q  P  L  P  T  P  E  L  Q  M  P  T  M  V  P  S  Y  
ENSG00000      1620 GAAAGGGACCAACCACTTCCTACCCCAGAGCTCCAGATGCCTACCATGGTGCCTTCTTAT

ENSG00000       560 D  L  G  M  A  P  D  S  S  M  S  M  Q  L  G  P  D  M  V  P  
ENSG00000      1680 GACCTTGGAATGGCCCCTGATTCCTCCATGAGCATGCAGCTTGGCCCAGATATGGTGCCC

ENSG00000       580 Q  V  Y  P  P  H  S  H  S  I  P  P  Y  Q  G  L  S  P  E  E  
ENSG00000      1740 CAGGTGTACCCACCACACTCTCACTCCATCCCCCCGTATCAAGGCCTCTCCCCAGAAGAA

ENSG00000       600 S  V  N  V  L  S  A  F  Q  E  P  H  L  Q  M  P  P  S  L  G  
ENSG00000      1800 TCAGTCAACGTGTTGTCAGCCTTCCAGGAGCCTCACCTGCAGATGCCCCCCAGCCTGGGC

ENSG00000       620 Q  M  S  L  P  F  D  Q  P  H  P  Q  G  L  L  P  C  Q  P  Q  
ENSG00000      1860 CAGATGAGCCTGCCCTTTGACCAGCCTCACCCCCAGGGCCTGCTGCCGTGCCAGCCTCAG

ENSG00000       640 E  H  A  V  S  S  P  D  P  L  L  C  S  D  V  T  M  V  E  D  
ENSG00000      1920 GAGCATGCTGTGTCCAGCCCTGACCCCCTGCTCTGCTCAGATGTGACCATGGTGGAAGAC

ENSG00000       660 S  C  L  S  Q  P  V  T  A  F  P  Q  G  T  W  I  G  E  D  I  
ENSG00000      1980 AGCTGCCTGAGCCAGCCAGTGACAGCGTTTCCTCAGGGCACTTGGATTGGTGAAGACATA

ENSG00000       680 F  P  P  L  L  P  P  T  E  Q  D  L  T  K  L  L  L  E  G  Q  
ENSG00000      2040 TTCCCTCCTCTGCTGCCTCCCACTGAACAGGACCTCACTAAGCTTCTCCTGGAGGGGCAA

ENSG00000       700 G  E  S  G  G  G  S  L  G  A  Q  P  L  L  Q  P  S  H  Y  G  
ENSG00000      2100 GGGGAGTCGGGGGGAGGGTCCTTGGGGGCACAGCCCCTCCTGCAGCCCTCCCACTATGGG

ENSG00000       720 Q  S  G  I  S  M  S  H  M  D  L  R  A  N  P  S  W    737
ENSG00000      2160 CAATCTGGGATCTCAATGTCCCACATGGACCTAAGGGCCAACCCCAGTTGG 2211
""",
        )
        protein_record = protein_alignment.sequences[5]
        nucleotide_record = nucleotide_records[protein_record.id]
        self.assertEqual(nucleotide_record.id, protein_record.id)
        alignments = aligner.align(protein_record, nucleotide_record)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 737], [0, 2211]]))
        )
        self.assertEqual(
            str(alignment),
            """\
ENSG00000         0 M  E  Q  F  R  H  L  P  M  P  F  H  W  K  Q  E  E  L  K  F  
ENSG00000         0 ATGGAACAGTTCCGCCACTTGCCAATGCCTTTCCACTGGAAGCAGGAAGAACTCAAGTTT

ENSG00000        20 K  T  G  L  R  R  L  Q  H  R  V  G  E  I  H  L  L  R  E  A  
ENSG00000        60 AAGACAGGCTTGCGGAGGCTGCAGCACCGAGTAGGGGAGATCCACCTTCTCCGAGAAGCC

ENSG00000        40 L  Q  K  G  A  E  A  G  Q  V  S  L  H  S  L  I  E  T  P  A  
ENSG00000       120 CTGCAGAAGGGGGCTGAGGCTGGCCAAGTGTCTCTGCACAGCTTGATAGAAACTCCTGCT

ENSG00000        60 N  G  T  G  P  S  E  A  L  A  M  L  L  Q  E  T  T  G  E  L  
ENSG00000       180 AATGGGACTGGGCCAAGTGAGGCCCTGGCCATGCTACTGCAGGAGACCACTGGAGAGCTA

ENSG00000        80 E  A  A  K  A  L  V  L  K  R  I  Q  I  W  K  R  Q  Q  Q  L  
ENSG00000       240 GAGGCAGCCAAAGCCCTAGTGCTGAAGAGGATCCAGATTTGGAAACGGCAGCAGCAGCTG

ENSG00000       100 A  G  N  G  A  P  F  E  E  S  L  A  P  L  Q  E  R  C  E  S  
ENSG00000       300 GCAGGGAATGGCGCACCGTTTGAGGAGAGCCTGGCCCCACTCCAGGAGAGGTGTGAAAGC

ENSG00000       120 L  V  D  I  Y  S  Q  L  Q  Q  E  V  G  A  A  G  G  E  L  E  
ENSG00000       360 CTGGTGGACATTTATTCCCAGCTACAGCAGGAGGTAGGGGCGGCTGGTGGGGAGCTTGAG

ENSG00000       140 P  K  T  R  A  S  L  T  G  R  L  D  E  V  L  R  T  L  V  T  
ENSG00000       420 CCCAAGACCCGGGCATCGCTGACTGGCCGGCTGGATGAAGTCCTGAGAACCCTCGTCACC

ENSG00000       160 S  C  F  L  V  E  K  Q  P  P  Q  V  L  K  T  Q  T  K  F  Q  
ENSG00000       480 AGTTGCTTCCTGGTGGAGAAGCAGCCCCCCCAGGTACTGAAGACTCAGACCAAGTTCCAG

ENSG00000       180 A  G  V  R  F  L  L  G  L  R  F  L  G  A  P  A  K  P  P  L  
ENSG00000       540 GCTGGAGTTCGATTCCTGTTGGGCTTGAGGTTCCTGGGGGCCCCAGCCAAGCCTCCGCTG

ENSG00000       200 V  R  A  D  M  V  T  E  K  Q  A  R  E  L  S  V  P  Q  G  P  
ENSG00000       600 GTCAGGGCCGACATGGTGACAGAGAAGCAGGCGCGGGAGCTGAGTGTGCCTCAGGGTCCT

ENSG00000       220 G  A  G  A  E  S  T  G  E  I  I  N  N  T  V  P  L  E  N  S  
ENSG00000       660 GGGGCTGGAGCAGAAAGCACTGGAGAAATCATCAACAACACTGTGCCCTTGGAGAACAGC

ENSG00000       240 I  P  G  N  C  C  S  A  L  F  K  N  L  L  L  K  K  I  K  R  
ENSG00000       720 ATTCCTGGGAACTGCTGCTCTGCCCTGTTCAAGAACCTGCTTCTCAAGAAGATCAAGCGG

ENSG00000       260 C  E  R  K  G  T  E  S  V  T  E  E  K  C  A  V  L  F  S  A  
ENSG00000       780 TGTGAGCGGAAGGGCACTGAGTCTGTCACAGAGGAGAAGTGCGCTGTGCTCTTCTCTGCC

ENSG00000       280 S  F  T  L  G  P  G  K  L  P  I  Q  L  Q  A  L  S  L  P  L  
ENSG00000       840 AGCTTCACACTTGGCCCCGGCAAACTCCCCATCCAGCTCCAGGCCCTGTCTCTGCCCCTG

ENSG00000       300 V  V  I  V  H  G  N  Q  D  N  N  A  K  A  T  I  L  W  D  N  
ENSG00000       900 GTGGTCATCGTCCATGGCAACCAAGACAACAATGCCAAAGCCACTATCCTGTGGGACAAT

ENSG00000       320 A  F  S  E  M  D  R  V  P  F  V  V  A  E  R  V  P  W  E  K  
ENSG00000       960 GCCTTCTCTGAGATGGACCGCGTGCCCTTTGTGGTGGCTGAGCGGGTGCCCTGGGAGAAG

ENSG00000       340 M  C  E  T  L  N  L  K  F  M  A  E  V  G  T  N  R  G  L  L  
ENSG00000      1020 ATGTGTGAAACTCTGAACCTGAAGTTCATGGCTGAGGTGGGGACCAACCGGGGGCTGCTC

ENSG00000       360 P  E  H  F  L  F  L  A  Q  K  I  F  N  D  N  S  L  S  M  E  
ENSG00000      1080 CCAGAGCACTTCCTCTTCCTGGCCCAGAAGATCTTCAATGACAACAGCCTCAGTATGGAG

ENSG00000       380 A  F  Q  H  R  S  V  S  W  S  Q  F  N  K  E  I  L  L  G  R  
ENSG00000      1140 GCCTTCCAGCACCGTTCTGTGTCCTGGTCGCAGTTCAACAAGGAGATCCTGCTGGGCCGT

ENSG00000       400 G  F  T  F  W  Q  W  F  D  G  V  L  D  L  T  K  R  C  L  R  
ENSG00000      1200 GGCTTCACCTTTTGGCAGTGGTTTGATGGTGTCCTGGACCTCACCAAACGCTGTCTCCGG

ENSG00000       420 S  Y  W  S  D  R  L  I  I  G  F  I  S  K  Q  Y  V  T  S  L  
ENSG00000      1260 AGCTACTGGTCTGACCGGCTGATCATTGGCTTCATCAGCAAACAGTACGTTACTAGCCTT

ENSG00000       440 L  L  N  E  P  D  G  T  F  L  L  R  F  S  D  S  E  I  G  G  
ENSG00000      1320 CTTCTCAATGAGCCCGACGGAACCTTTCTCCTCCGCTTCAGCGACTCAGAGATTGGGGGC

ENSG00000       460 I  T  I  A  H  V  I  R  G  Q  D  G  S  P  Q  I  E  N  I  Q  
ENSG00000      1380 ATCACCATTGCCCATGTCATCCGGGGCCAGGATGGCTCTCCACAGATAGAGAACATCCAG

ENSG00000       480 P  F  S  A  K  D  L  S  I  R  S  L  G  D  R  I  R  D  L  A  
ENSG00000      1440 CCATTCTCTGCCAAAGACCTGTCCATTCGCTCACTGGGGGACCGAATCCGGGATCTTGCT

ENSG00000       500 Q  L  K  N  L  Y  P  K  K  P  K  D  E  A  F  R  S  H  Y  K  
ENSG00000      1500 CAGCTCAAAAATCTCTATCCCAAGAAGCCCAAGGATGAGGCTTTCCGGAGCCACTACAAG

ENSG00000       520 P  E  Q  M  G  K  D  G  R  G  Y  V  P  A  T  I  K  M  T  V  
ENSG00000      1560 CCTGAACAGATGGGTAAGGATGGCAGGGGTTATGTCCCAGCTACCATCAAGATGACCGTG

ENSG00000       540 E  R  D  Q  P  L  P  T  P  E  L  Q  M  P  T  M  V  P  S  Y  
ENSG00000      1620 GAAAGGGACCAACCACTTCCTACCCCAGAGCTCCAGATGCCTACCATGGTGCCTTCTTAT

ENSG00000       560 D  L  G  M  A  P  D  S  S  M  S  M  Q  L  G  P  D  M  V  P  
ENSG00000      1680 GACCTTGGAATGGCCCCTGATTCCTCCATGAGCATGCAGCTTGGCCCAGATATGGTGCCC

ENSG00000       580 Q  V  Y  P  P  H  S  H  S  I  P  P  Y  Q  G  L  S  P  E  E  
ENSG00000      1740 CAGGTGTACCCACCACACTCTCACTCCATCCCCCCGTATCAAGGCCTCTCCCCAGAAGAA

ENSG00000       600 S  V  N  V  L  S  A  F  Q  E  P  H  L  Q  M  P  P  S  L  G  
ENSG00000      1800 TCAGTCAACGTGTTGTCAGCCTTCCAGGAGCCTCACCTGCAGATGCCCCCCAGCCTGGGC

ENSG00000       620 Q  M  S  L  P  F  D  Q  P  H  P  Q  G  L  L  P  C  Q  P  Q  
ENSG00000      1860 CAGATGAGCCTGCCCTTTGACCAGCCTCACCCCCAGGGCCTGCTGCCGTGCCAGCCTCAG

ENSG00000       640 E  H  A  V  S  S  P  D  P  L  L  C  S  D  V  T  M  V  E  D  
ENSG00000      1920 GAGCATGCTGTGTCCAGCCCTGACCCCCTGCTCTGCTCAGATGTGACCATGGTGGAAGAC

ENSG00000       660 S  C  L  S  Q  P  V  T  A  F  P  Q  G  T  W  I  G  E  D  I  
ENSG00000      1980 AGCTGCCTGAGCCAGCCAGTGACAGCGTTTCCTCAGGGCACTTGGATTGGTGAAGACATA

ENSG00000       680 F  P  P  L  L  P  P  T  E  Q  D  L  T  K  L  L  L  E  G  Q  
ENSG00000      2040 TTCCCTCCTCTGCTGCCTCCCACTGAACAGGACCTCACTAAGCTTCTCCTGGAGGGGCAA

ENSG00000       700 G  E  S  G  G  G  S  L  G  A  Q  P  L  L  Q  P  S  H  Y  G  
ENSG00000      2100 GGGGAGTCGGGGGGAGGGTCCTTGGGGGCACAGCCCCTCCTGCAGCCCTCCCACTATGGG

ENSG00000       720 Q  S  G  I  S  M  S  H  M  D  L  R  A  N  P  S  W    737
ENSG00000      2160 CAATCTGGGATCTCAATGTCCCACATGGACCTAAGGGCCAACCCCAGTTGG 2211
""",
        )
        protein_record = protein_alignment.sequences[6]
        nucleotide_record = nucleotide_records[protein_record.id]
        self.assertEqual(nucleotide_record.id, protein_record.id)
        alignments = aligner.align(protein_record, nucleotide_record)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 737], [0, 2211]]))
        )
        self.assertEqual(
            str(alignment),
            """\
ENSG00000         0 M  E  Q  F  R  H  L  P  M  P  F  H  W  K  Q  E  E  L  K  F  
ENSG00000         0 ATGGAACAGTTCCGCCACTTGCCAATGCCTTTCCACTGGAAGCAGGAAGAACTCAAGTTT

ENSG00000        20 K  T  G  L  R  R  L  Q  H  R  V  G  E  I  H  L  L  R  E  A  
ENSG00000        60 AAGACAGGCTTGCGGAGGCTGCAGCACCGAGTAGGGGAGATCCACCTTCTCCGAGAAGCC

ENSG00000        40 L  Q  K  G  A  E  A  G  Q  V  S  L  H  S  L  I  E  T  P  A  
ENSG00000       120 CTGCAGAAGGGGGCTGAGGCTGGCCAAGTGTCTCTGCACAGCTTGATAGAAACTCCTGCT

ENSG00000        60 N  G  T  G  P  S  E  A  L  A  M  L  L  Q  E  T  T  G  E  L  
ENSG00000       180 AATGGGACTGGGCCAAGTGAGGCCCTGGCCATGCTACTGCAGGAGACCACTGGAGAGCTA

ENSG00000        80 E  A  A  K  A  L  V  L  K  R  I  Q  I  W  K  R  Q  Q  Q  L  
ENSG00000       240 GAGGCAGCCAAAGCCCTAGTGCTGAAGAGGATCCAGATTTGGAAACGGCAGCAGCAGCTG

ENSG00000       100 A  G  N  G  A  P  F  E  E  S  L  A  P  L  Q  E  R  C  E  S  
ENSG00000       300 GCAGGGAATGGCGCACCGTTTGAGGAGAGCCTGGCCCCACTCCAGGAGAGGTGTGAAAGC

ENSG00000       120 L  V  D  I  Y  S  Q  L  Q  Q  E  V  G  A  A  G  G  E  L  E  
ENSG00000       360 CTGGTGGACATTTATTCCCAGCTACAGCAGGAGGTAGGGGCGGCTGGTGGGGAGCTTGAG

ENSG00000       140 P  K  T  R  A  S  L  T  G  R  L  D  E  V  L  R  T  L  V  T  
ENSG00000       420 CCCAAGACCCGGGCATCGCTGACTGGCCGGCTGGATGAAGTCCTGAGAACCCTCGTCACC

ENSG00000       160 S  C  F  L  V  E  K  Q  P  P  Q  V  L  K  T  Q  T  K  F  Q  
ENSG00000       480 AGTTGCTTCCTGGTGGAGAAGCAGCCCCCCCAGGTACTGAAGACTCAGACCAAGTTCCAG

ENSG00000       180 A  G  V  R  F  L  L  G  L  R  F  L  G  A  P  A  K  P  P  L  
ENSG00000       540 GCTGGAGTTCGATTCCTGTTGGGCTTGAGGTTCCTGGGGGCCCCAGCCAAGCCTCCGCTG

ENSG00000       200 V  R  A  D  M  V  T  E  K  Q  A  R  E  L  S  V  P  Q  G  P  
ENSG00000       600 GTCAGGGCCGACATGGTGACAGAGAAGCAGGCGCGGGAGCTGAGTGTGCCTCAGGGTCCT

ENSG00000       220 G  A  G  A  E  S  T  G  E  I  I  N  N  T  V  P  L  E  N  S  
ENSG00000       660 GGGGCTGGAGCAGAAAGCACTGGAGAAATCATCAACAACACTGTGCCCTTGGAGAACAGC

ENSG00000       240 I  P  G  N  C  C  S  A  L  F  K  N  L  L  L  K  K  I  K  R  
ENSG00000       720 ATTCCTGGGAACTGCTGCTCTGCCCTGTTCAAGAACCTGCTTCTCAAGAAGATCAAGCGG

ENSG00000       260 C  E  R  K  G  T  E  S  V  T  E  E  K  C  A  V  L  F  S  A  
ENSG00000       780 TGTGAGCGGAAGGGCACTGAGTCTGTCACAGAGGAGAAGTGCGCTGTGCTCTTCTCTGCC

ENSG00000       280 S  F  T  L  G  P  G  K  L  P  I  Q  L  Q  A  L  S  L  P  L  
ENSG00000       840 AGCTTCACACTTGGCCCCGGCAAACTCCCCATCCAGCTCCAGGCCCTGTCTCTGCCCCTG

ENSG00000       300 V  V  I  V  H  G  N  Q  D  N  N  A  K  A  T  I  L  W  D  N  
ENSG00000       900 GTGGTCATCGTCCATGGCAACCAAGACAACAATGCCAAAGCCACTATCCTGTGGGACAAT

ENSG00000       320 A  F  S  E  M  D  R  V  P  F  V  V  A  E  R  V  P  W  E  K  
ENSG00000       960 GCCTTCTCTGAGATGGACCGCGTGCCCTTTGTGGTGGCTGAGCGGGTGCCCTGGGAGAAG

ENSG00000       340 M  C  E  T  L  N  L  K  F  M  A  E  V  G  T  N  R  G  L  L  
ENSG00000      1020 ATGTGTGAAACTCTGAACCTGAAGTTCATGGCTGAGGTGGGGACCAACCGGGGGCTGCTC

ENSG00000       360 P  E  H  F  L  F  L  A  Q  K  I  F  N  D  N  S  L  S  M  E  
ENSG00000      1080 CCAGAGCACTTCCTCTTCCTGGCCCAGAAGATCTTCAATGACAACAGCCTCAGTATGGAG

ENSG00000       380 A  F  Q  H  R  S  V  S  W  S  Q  F  N  K  E  I  L  L  G  R  
ENSG00000      1140 GCCTTCCAGCACCGTTCTGTGTCCTGGTCGCAGTTCAACAAGGAGATCCTGCTGGGCCGT

ENSG00000       400 G  F  T  F  W  Q  W  F  D  G  V  L  D  L  T  K  R  C  L  R  
ENSG00000      1200 GGCTTCACCTTTTGGCAGTGGTTTGATGGTGTCCTGGACCTCACCAAACGCTGTCTCCGG

ENSG00000       420 S  Y  W  S  D  R  L  I  I  G  F  I  S  K  Q  Y  V  T  S  L  
ENSG00000      1260 AGCTACTGGTCTGACCGGCTGATCATTGGCTTCATCAGCAAACAGTACGTTACTAGCCTT

ENSG00000       440 L  L  N  E  P  D  G  T  F  L  L  R  F  S  D  S  E  I  G  G  
ENSG00000      1320 CTTCTCAATGAGCCCGACGGAACCTTTCTCCTCCGCTTCAGCGACTCAGAGATTGGGGGC

ENSG00000       460 I  T  I  A  H  V  I  R  G  Q  D  G  S  P  Q  I  E  N  I  Q  
ENSG00000      1380 ATCACCATTGCCCATGTCATCCGGGGCCAGGATGGCTCTCCACAGATAGAGAACATCCAG

ENSG00000       480 P  F  S  A  K  D  L  S  I  R  S  L  G  D  R  I  R  D  L  A  
ENSG00000      1440 CCATTCTCTGCCAAAGACCTGTCCATTCGCTCACTGGGGGACCGAATCCGGGATCTTGCT

ENSG00000       500 Q  L  K  N  L  Y  P  K  K  P  K  D  E  A  F  R  S  H  Y  K  
ENSG00000      1500 CAGCTCAAAAATCTCTATCCCAAGAAGCCCAAGGATGAGGCTTTCCGGAGCCACTACAAG

ENSG00000       520 P  E  Q  M  G  K  D  G  R  G  Y  V  P  A  T  I  K  M  T  V  
ENSG00000      1560 CCTGAACAGATGGGTAAGGATGGCAGGGGTTATGTCCCAGCTACCATCAAGATGACCGTG

ENSG00000       540 E  R  D  Q  P  L  P  T  P  E  L  Q  M  P  T  M  V  P  S  Y  
ENSG00000      1620 GAAAGGGACCAACCACTTCCTACCCCAGAGCTCCAGATGCCTACCATGGTGCCTTCTTAT

ENSG00000       560 D  L  G  M  A  P  D  S  S  M  S  M  Q  L  G  P  D  M  V  P  
ENSG00000      1680 GACCTTGGAATGGCCCCTGATTCCTCCATGAGCATGCAGCTTGGCCCAGATATGGTGCCC

ENSG00000       580 Q  V  Y  P  P  H  S  H  S  I  P  P  Y  Q  G  L  S  P  E  E  
ENSG00000      1740 CAGGTGTACCCACCACACTCTCACTCCATCCCCCCGTATCAAGGCCTCTCCCCAGAAGAA

ENSG00000       600 S  V  N  V  L  S  A  F  Q  E  P  H  L  Q  M  P  P  S  L  G  
ENSG00000      1800 TCAGTCAACGTGTTGTCAGCCTTCCAGGAGCCTCACCTGCAGATGCCCCCCAGCCTGGGC

ENSG00000       620 Q  M  S  L  P  F  D  Q  P  H  P  Q  G  L  L  P  C  Q  P  Q  
ENSG00000      1860 CAGATGAGCCTGCCCTTTGACCAGCCTCACCCCCAGGGCCTGCTGCCGTGCCAGCCTCAG

ENSG00000       640 E  H  A  V  S  S  P  D  P  L  L  C  S  D  V  T  M  V  E  D  
ENSG00000      1920 GAGCATGCTGTGTCCAGCCCTGACCCCCTGCTCTGCTCAGATGTGACCATGGTGGAAGAC

ENSG00000       660 S  C  L  S  Q  P  V  T  A  F  P  Q  G  T  W  I  G  E  D  I  
ENSG00000      1980 AGCTGCCTGAGCCAGCCAGTGACAGCGTTTCCTCAGGGCACTTGGATTGGTGAAGACATA

ENSG00000       680 F  P  P  L  L  P  P  T  E  Q  D  L  T  K  L  L  L  E  G  Q  
ENSG00000      2040 TTCCCTCCTCTGCTGCCTCCCACTGAACAGGACCTCACTAAGCTTCTCCTGGAGGGGCAA

ENSG00000       700 G  E  S  G  G  G  S  L  G  A  Q  P  L  L  Q  P  S  H  Y  G  
ENSG00000      2100 GGGGAGTCGGGGGGAGGGTCCTTGGGGGCACAGCCCCTCCTGCAGCCCTCCCACTATGGG

ENSG00000       720 Q  S  G  I  S  M  S  H  M  D  L  R  A  N  P  S  W    737
ENSG00000      2160 CAATCTGGGATCTCAATGTCCCACATGGACCTAAGGGCCAACCCCAGTTGG 2211
""",
        )
        protein_record = protein_alignment.sequences[7]
        nucleotide_record = nucleotide_records[protein_record.id]
        self.assertEqual(nucleotide_record.id, protein_record.id)
        nucleotide_record = nucleotide_record.upper()
        alignments = aligner.align(protein_record, nucleotide_record)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 1021], [0, 3063]]))
        )
        self.assertEqual(
            str(alignment),
            """\
isotig466         0 E  V  T  Q  S  R  R  K  P  V  R  E  G  R  P  W  E  P  S  Q  
isotig466         0 AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

isotig466        20 E  L  A  Q  T  P  E  G  G  C  P  R  Q  E  S  R  P  E  G  G  
isotig466        60 GAGTTGGCCCAAACCCCTGAGGGAGGCTGTCCTCGGCAGGAGAGCAGACCTGAGGGAGGT

isotig466        40 W  C  A  G  P  G  C  P  A  L  P  S  H  P  R  A  G  S  P  E  
isotig466       120 TGGTGTGCCGGTCCTGGTTGTCCAGCCCTCCCGTCCCATCCCAGGGCTGGCTCACCAGAA

isotig466        60 K  S  T  G  I  P  G  P  L  G  R  E  Q  R  V  S  S  C  P  G  
isotig466       180 AAAAGTACAGGCATCCCAGGCCCTCTGGGGCGGGAGCAGAGGGTCTCCTCTTGTCCGGGG

isotig466        80 E  D  K  N  K  N  K  K  R  V  C  S  P  H  W  K  S  S  G  A  
isotig466       240 GAAGACAAAAACAAAAACAAAAAACGAGTGTGTAGCCCTCACTGGAAGTCTTCTGGTGCT

isotig466       100 L  G  P  F  A  L  G  S  L  A  S  G  H  G  G  R  A  P  P  G  
isotig466       300 CTGGGGCCGTTTGCACTTGGGAGCCTGGCTTCTGGGCATGGTGGCCGGGCTCCGCCGGGG

isotig466       120 A  S  A  L  G  K  A  F  L  E  Q  S  G  W  E  P  T  E  V  P  
isotig466       360 GCTTCAGCCCTCGGCAAAGCGTTTCTAGAACAGAGTGGGTGGGAGCCGACTGAAGTCCCA

isotig466       140 E  P  R  T  L  T  H  R  K  A  S  S  G  W  R  N  A  N  P  P  
isotig466       420 GAACCGCGAACACTGACGCACAGGAAAGCCTCGAGTGGGTGGAGAAATGCAAATCCCCCA

isotig466       160 V  A  W  P  Q  P  L  R  T  L  T  A  D  S  D  F  G  L  E  R  
isotig466       480 GTGGCATGGCCTCAGCCGCTGCGGACCCTGACGGCCGATTCTGACTTCGGACTTGAGCGG

isotig466       180 D  Q  H  P  D  P  G  R  A  S  M  S  Q  W  N  Q  V  Q  Q  L  
isotig466       540 GACCAGCACCCGGACCCAGGGCGAGCCAGCATGTCTCAGTGGAATCAAGTCCAACAGTTA

isotig466       200 E  I  K  F  L  E  Q  V  D  Q  F  Y  D  D  N  F  P  M  E  I  
isotig466       600 GAAATCAAGTTTTTGGAGCAAGTTGATCAGTTCTATGATGACAACTTTCCCATGGAAATC

isotig466       220 R  H  L  L  A  Q  W  I  E  H  Q  D  W  E  V  A  S  N  N  E  
isotig466       660 CGACATCTGCTGGCCCAGTGGATTGAGCATCAAGACTGGGAGGTGGCCTCTAACAATGAA

isotig466       240 T  M  A  T  I  L  L  Q  N  L  L  I  Q  L  D  E  Q  L  G  R  
isotig466       720 ACTATGGCAACAATTCTTCTTCAAAACTTATTAATACAATTGGATGAACAGTTAGGTCGT

isotig466       260 V  S  K  E  K  N  L  L  L  I  H  N  L  K  R  I  R  K  V  L  
isotig466       780 GTTTCCAAAGAGAAAAACCTGCTATTGATCCACAATCTAAAGAGAATTAGAAAAGTACTT

isotig466       280 Q  G  K  F  H  G  N  P  M  H  V  A  V  V  I  S  N  C  L  R  
isotig466       840 CAGGGGAAGTTTCATGGAAATCCAATGCATGTAGCCGTGGTAATCTCAAATTGTTTAAGG

isotig466       300 E  E  R  R  I  L  A  A  A  N  M  P  I  Q  G  P  L  E  K  S  
isotig466       900 GAAGAGAGGAGAATACTGGCTGCAGCGAACATGCCTATCCAGGGACCTCTGGAGAAATCC

isotig466       320 L  Q  S  S  S  V  S  E  R  Q  R  N  V  E  H  K  V  A  A  I  
isotig466       960 TTACAAAGTTCGTCGGTTTCAGAAAGACAGAGAAATGTGGAACACAAAGTGGCTGCCATT

isotig466       340 K  N  S  V  Q  M  T  E  Q  D  T  K  Y  L  E  D  L  Q  D  E  
isotig466      1020 AAAAACAGTGTGCAGATGACAGAACAAGACACCAAATACTTGGAAGATCTGCAAGATGAA

isotig466       360 F  D  Y  R  Y  K  T  I  Q  T  M  D  Q  G  D  K  N  S  I  L  
isotig466      1080 TTTGACTACAGGTATAAAACAATTCAGACAATGGACCAGGGTGACAAGAATAGCATCCTA

isotig466       380 M  N  Q  E  V  L  T  L  Q  E  M  L  N  S  L  D  F  K  R  K  
isotig466      1140 ATGAACCAGGAGGTTTTGACACTCCAAGAAATGCTTAATAGCCTGGACTTCAAGAGAAAG

isotig466       400 E  A  L  T  K  M  T  Q  I  V  N  E  S  D  L  L  M  S  S  M  
isotig466      1200 GAAGCACTCACTAAGATGACACAGATAGTGAACGAGTCGGACCTGCTGATGAGCAGCATG

isotig466       420 L  I  E  E  L  Q  D  W  K  R  R  Q  Q  I  A  C  I  G  G  P  
isotig466      1260 CTCATAGAAGAGCTGCAGGACTGGAAGAGGAGGCAGCAGATCGCCTGCATCGGTGGCCCA

isotig466       440 L  H  N  G  L  D  Q  L  Q  N  C  F  T  L  L  A  E  S  L  F  
isotig466      1320 CTCCACAACGGGCTGGACCAGCTTCAGAACTGCTTTACCCTGTTGGCAGAAAGTCTTTTC

isotig466       460 Q  L  R  R  Q  L  E  K  L  E  E  Q  S  S  K  M  T  Y  E  G  
isotig466      1380 CAACTCAGACGACAGCTGGAGAAATTAGAGGAGCAGTCTTCCAAGATGACTTACGAAGGA

isotig466       480 D  P  I  P  T  Q  R  A  H  L  L  E  R  A  T  F  L  I  Y  N  
isotig466      1440 GACCCCATCCCCACGCAGAGAGCACACCTGCTGGAGAGAGCCACCTTCCTGATCTACAAC

isotig466       500 L  F  K  N  S  F  V  V  E  R  Q  P  C  M  P  T  H  P  Q  R  
isotig466      1500 CTTTTCAAGAACTCATTTGTGGTTGAGCGACAGCCCTGCATGCCAACACACCCTCAGAGG

isotig466       520 P  L  V  L  K  T  L  I  Q  F  T  A  K  L  R  L  L  I  K  L  
isotig466      1560 CCGCTGGTACTCAAAACCCTCATTCAGTTCACCGCGAAACTGAGACTACTAATAAAATTG

isotig466       540 P  E  L  N  Y  Q  V  K  V  K  A  S  I  D  K  N  V  S  T  L  
isotig466      1620 CCGGAACTCAACTATCAGGTGAAAGTAAAGGCATCGATCGACAAGAATGTTTCAACGCTA

isotig466       560 S  N  R  R  F  V  L  C  G  T  Q  V  K  A  M  S  I  E  E  S  
isotig466      1680 AGCAATAGAAGATTTGTGCTTTGTGGAACTCAAGTCAAAGCCATGTCCATCGAGGAATCC

isotig466       580 S  N  G  S  L  S  V  E  F  R  H  L  Q  P  K  E  M  K  S  S  
isotig466      1740 TCCAATGGGAGCCTCTCAGTAGAATTTAGACATTTGCAACCGAAGGAAATGAAATCCAGT

isotig466       600 A  G  S  K  G  N  E  G  C  H  M  V  T  E  E  L  H  S  I  A  
isotig466      1800 GCCGGAAGTAAAGGAAATGAGGGCTGCCACATGGTGACGGAAGAGCTGCATTCCATAGCC

isotig466       620 F  E  T  Q  I  C  L  Y  G  L  T  I  D  L  E  T  S  S  L  P  
isotig466      1860 TTTGAGACCCAGATCTGCCTCTATGGCCTCACCATCGACTTGGAGACAAGCTCATTACCT

isotig466       640 V  V  M  I  S  N  V  S  Q  L  P  N  A  W  A  S  I  I  W  Y  
isotig466      1920 GTGGTGATGATTTCTAATGTCAGCCAACTGCCTAATGCTTGGGCATCCATCATTTGGTAC

isotig466       660 N  V  S  T  N  D  C  Q  N  L  V  F  F  N  N  P  P  P  V  T  
isotig466      1980 AATGTGTCAACCAACGATTGCCAGAACTTGGTTTTCTTTAATAATCCTCCGCCTGTCACT

isotig466       680 L  S  Q  L  L  E  V  M  S  W  Q  F  S  S  Y  V  G  R  G  L  
isotig466      2040 TTGAGTCAACTCCTGGAAGTGATGAGCTGGCAGTTTTCATCCTATGTTGGTCGTGGCCTT

isotig466       700 N  S  D  Q  L  N  M  L  A  E  K  L  T  V  Q  S  N  Y  S  D  
isotig466      2100 AATTCAGACCAGCTCAACATGCTGGCAGAGAAGCTCACAGTTCAGTCTAACTACAGCGAT

isotig466       720 G  H  L  T  W  A  K  F  C  K  E  H  L  P  G  K  P  F  T  F  
isotig466      2160 GGTCACCTCACCTGGGCCAAGTTCTGCAAGGAACACTTGCCTGGCAAACCATTTACCTTC

isotig466       740 W  T  W  L  E  A  I  L  D  L  I  K  K  H  I  L  P  L  W  I  
isotig466      2220 TGGACCTGGCTTGAAGCAATATTGGACCTAATTAAAAAACACATTCTTCCCCTCTGGATT

isotig466       760 D  G  Y  I  M  G  F  V  S  K  E  K  E  R  F  L  L  K  D  K  
isotig466      2280 GATGGGTACATCATGGGCTTCGTGAGCAAAGAGAAGGAGAGGTTTCTGCTCAAGGATAAA

isotig466       780 M  P  G  T  F  L  L  R  F  S  E  S  H  L  G  G  I  T  F  T  
isotig466      2340 ATGCCCGGGACATTTTTGTTACGATTCAGTGAGAGCCATCTCGGAGGGATCACCTTCACC

isotig466       800 W  V  D  H  S  E  N  G  E  V  R  F  H  S  V  E  P  Y  N  K  
isotig466      2400 TGGGTGGACCACTCTGAAAACGGAGAAGTGAGATTCCACTCCGTAGAACCCTACAACAAA

isotig466       820 G  R  L  S  A  L  P  F  A  D  I  L  R  D  Y  K  V  I  M  A  
isotig466      2460 GGGCGTCTGTCGGCCCTGCCATTTGCTGACATCCTGCGGGACTACAAGGTCATCATGGCT

isotig466       840 E  N  I  P  E  N  P  L  K  Y  L  Y  P  D  I  P  K  D  K  A  
isotig466      2520 GAGAACATTCCCGAGAACCCTCTCAAGTACCTCTACCCCGACATCCCCAAAGACAAAGCC

isotig466       860 F  G  K  H  Y  S  S  Q  P  C  E  V  S  R  P  T  E  R  G  D  
isotig466      2580 TTCGGTAAACACTACAGCTCCCAGCCTTGCGAAGTTTCAAGGCCAACAGAACGGGGAGAC

isotig466       880 K  G  Y  V  P  S  V  F  I  P  I  S  T  I  R  S  D  A  M  E  
isotig466      2640 AAAGGTTATGTTCCTTCAGTTTTTATCCCTATTTCAACAATCCGCAGCGACGCCATGGAG

isotig466       900 P  Q  S  P  S  D  L  L  P  M  S  P  S  V  Y  A  V  L  R  E  
isotig466      2700 CCGCAGTCTCCTTCAGACCTTCTCCCCATGTCTCCGAGTGTATACGCTGTGCTGAGAGAA

isotig466       920 N  L  S  P  T  T  I  E  T  A  M  K  S  P  Y  S  E  R  Y  K  
isotig466      2760 AACCTGAGCCCTACCACAATTGAAACAGCAATGAAGTCTCCATATTCTGAGCGGTACAAA

isotig466       940 A  T  L  Q  G  R  E  Q  M  K  T  E  T  A  L  C  Q  S  P  Q  
isotig466      2820 GCGACTCTTCAAGGAAGAGAGCAGATGAAAACGGAGACTGCTCTTTGCCAAAGTCCACAA

isotig466       960 F  I  S  S  A  L  I  L  V  S  R  K  W  H  K  S  E  A  F  L  
isotig466      2880 TTCATTTCTTCAGCTTTGATACTGGTTTCTAGAAAATGGCACAAATCCGAAGCTTTCCTC

isotig466       980 S  L  G  D  I  P  Q  L  G  V  L  L  K  C  K  P  K  L  Q  I  
isotig466      2940 TCACTAGGTGACATTCCCCAACTGGGAGTGCTGCTGAAATGCAAACCAAAGCTTCAGATA

isotig466      1000 N  T  Q  E  K  T  A  S  R  N  L  C  S  Q  Y  N  R  R  L  L  
isotig466      3000 AACACGCAGGAAAAGACAGCTTCGAGAAACCTATGTTCGCAATATAACAGAAGGCTGCTT

isotig466      1020 C   1021
isotig466      3060 TGC 3063
""",
        )
        protein_record = protein_alignment.sequences[8]
        nucleotide_record = nucleotide_records[protein_record.id]
        self.assertEqual(nucleotide_record.id, protein_record.id)
        nucleotide_record = nucleotide_record.upper()
        alignments = aligner.align(protein_record, nucleotide_record)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 806], [0, 2418]]))
        )
        self.assertEqual(
            str(alignment),
            """\
isotig125         0 A  R  R  G  Q  A  A  L  G  S  P  A  A  R  T  W  S  Q  R  S  
isotig125         0 GCTAGGAGAGGCCAGGCGGCCCTCGGGAGCCCAGCTGCTCGCACCTGGAGCCAGCGCAGC

isotig125        20 P  A  S  R  A  S  A  R  E  T  V  T  P  P  D  C  G  R  M  A  
isotig125        60 CCGGCCAGTCGGGCCTCAGCCCGGGAGACAGTTACGCCCCCTGATTGCGGCAGGATGGCC

isotig125        40 Q  W  N  Q  L  Q  Q  L  D  T  R  Y  L  E  Q  L  H  Q  L  Y  
isotig125       120 CAGTGGAACCAGCTGCAGCAGCTGGACACTCGGTACCTGGAGCAGCTGCACCAGCTGTAC

isotig125        60 S  D  S  F  P  M  E  L  R  Q  F  L  A  P  W  I  E  S  Q  D  
isotig125       180 AGCGACAGCTTCCCCATGGAGCTGCGGCAGTTCCTGGCACCTTGGATTGAGAGTCAAGAT

isotig125        80 W  A  Y  A  A  S  K  E  S  H  A  T  L  V  F  H  N  L  L  G  
isotig125       240 TGGGCATATGCAGCCAGCAAAGAGTCACATGCCACACTGGTGTTTCATAATCTCTTGGGT

isotig125       100 E  I  D  Q  Q  Y  S  R  F  L  Q  E  S  N  V  L  Y  Q  H  N  
isotig125       300 GAGATTGACCAGCAGTACAGCCGATTCCTGCAAGAGTCCAACGTCCTCTATCAGCACAAC

isotig125       120 L  R  R  I  G  Q  F  V  Q  S  R  Y  L  E  K  P  M  E  I  A  
isotig125       360 CTTCGGAGGATCAAGCAGTTCCTACAGAGCAGGTATCTTGAGAAGCCGATGGAAATCGCC

isotig125       140 R  I  V  A  R  C  L  W  E  E  S  R  L  L  Q  T  A  A  T  A  
isotig125       420 CGCATCGTGGCCCGATGCCTGTGGGAAGAGTCTCGCCTCCTCCAGACGGCAGCCACTGCA

isotig125       160 A  Q  Q  G  G  Q  A  N  H  P  T  A  A  V  V  T  E  K  Q  Q  
isotig125       480 GCCCAGCAAGGAGGCCAGGCCAACCACCCAACAGCTGCTGTGGTGACGGAGAAACAGCAG

isotig125       180 M  L  E  Q  H  L  Q  D  V  R  K  R  V  Q  D  L  E  Q  K  M  
isotig125       540 ATGCTGGAGCAGCATCTTCAGGATGTCCGGAAACGTGTGCAGGATCTAGAACAGAAAATG

isotig125       200 K  V  V  E  N  L  Q  D  D  F  D  F  N  Y  K  T  L  K  S  Q  
isotig125       600 AAAGTGGTAGAGAATCTCCAGGATGACTTTGATTTCAACTATAAAACCCTCAAGAGTCAA

isotig125       220 G  D  M  Q  D  L  N  G  N  N  Q  S  V  T  R  Q  K  M  Q  Q  
isotig125       660 GGAGACATGCAGGATCTGAATGGAAACAACCAGTCTGTGACCAGGCAGAAGATGCAGCAG

isotig125       240 L  E  Q  M  L  T  A  L  D  Q  M  R  R  S  I  V  S  E  L  A  
isotig125       720 CTGGAACAGATGCTCACGGCGCTGGACCAGATGCGGAGAAGCATTGTGAGTGAGCTGGCG

isotig125       260 G  L  L  S  A  M  E  Y  V  Q  K  T  L  T  D  E  E  L  A  D  
isotig125       780 GGGCTTTTGTCGGCAATGGAGTACGTGCAGAAAACACTCACAGACGAGGAGCTGGCTGAC

isotig125       280 W  K  R  R  Q  Q  I  A  C  I  G  G  P  P  N  I  C  L  D  R  
isotig125       840 TGGAAGAGGCGGCAGCAGATCGCGTGCATTGGAGGCCCTCCCAACATCTGCCTGGATCGC

isotig125       300 L  E  N  W  I  T  S  L  A  E  S  Q  L  Q  T  R  Q  Q  I  K  
isotig125       900 CTGGAAAACTGGATAACTTCGTTAGCAGAATCTCAACTTCAGACCCGCCAACAAATTAAG

isotig125       320 K  L  E  E  L  Q  Q  K  V  S  Y  K  G  D  P  I  V  Q  H  R  
isotig125       960 AAACTGGAGGAGCTACAGCAGAAGGTGTCCTACAAGGGGGACCCCATTGTGCAGCACCGG

isotig125       340 P  M  L  E  E  R  I  V  E  L  F  R  N  L  M  K  S  A  F  V  
isotig125      1020 CCGATGCTGGAGGAGCGGATCGTGGAGCTGTTCAGAAACTTGATGAAGAGTGCCTTCGTG

isotig125       360 V  E  R  Q  P  C  M  P  M  H  P  D  R  P  L  V  I  K  T  G  
isotig125      1080 GTGGAGCGACAGCCCTGCATGCCGATGCACCCCGACCGGCCCTTGGTCATCAAGACTGGT

isotig125       380 V  Q  F  T  T  K  V  R  L  L  V  K  F  P  E  L  N  Y  Q  L  
isotig125      1140 GTCCAGTTCACTACTAAAGTCAGGTTGTTGGTCAAGTTTCCCGAGTTGAATTATCAGCTT

isotig125       400 K  I  K  V  C  I  D  K  D  S  G  D  V  A  A  L  R  G  S  R  
isotig125      1200 AAAATTAAAGTGTGCATTGACAAAGATTCTGGGGACGTTGCTGCTCTCAGAGGATCTCGG

isotig125       420 K  F  N  I  L  G  T  N  T  K  V  M  N  M  E  E  S  N  N  G  
isotig125      1260 AAATTTAACATTCTGGGCACAAACACGAAGGTGATGAACATGGAAGAATCCAACAACGGC

isotig125       440 S  L  S  A  E  F  K  H  L  T  L  R  E  Q  R  C  G  N  G  G  
isotig125      1320 AGCCTGTCTGCGGAGTTCAAGCACTTGACCCTGAGGGAGCAGAGATGTGGGAATGGAGGC

isotig125       460 R  A  N  C  D  A  S  L  I  V  T  E  E  L  H  L  I  T  F  E  
isotig125      1380 CGTGCCAATTGTGATGCCTCCTTGATTGTGACCGAGGAGCTGCATCTGATCACCTTCGAG

isotig125       480 T  E  V  Y  H  Q  G  L  K  I  D  L  E  T  H  S  L  P  V  V  
isotig125      1440 ACTGAGGTGTACCACCAAGGCCTCAAGATTGACCTGGAGACCCATTCTTTGCCAGTTGTG

isotig125       500 V  I  S  N  I  C  Q  M  P  N  A  W  A  S  I  L  W  Y  N  M  
isotig125      1500 GTGATCTCCAACATCTGTCAGATGCCAAATGCCTGGGCATCCATCCTGTGGTATAACATG

isotig125       520 L  T  N  N  P  K  N  V  N  F  F  T  K  P  P  I  G  T  W  D  
isotig125      1560 CTGACCAACAACCCCAAGAACGTGAACTTCTTCACCAAGCCACCAATCGGAACCTGGGAC

isotig125       540 Q  V  A  E  V  L  S  W  Q  F  S  S  T  T  K  R  G  L  S  I  
isotig125      1620 CAGGTGGCCGAGGTGCTCAGCTGGCAGTTCTCATCCACCACAAAGCGAGGGCTGAGCATC

isotig125       560 E  Q  L  T  T  L  A  E  K  L  L  G  P  G  V  N  Y  S  G  C  
isotig125      1680 GAGCAGCTGACTACGCTGGCCGAGAAGCTCCTAGGACCTGGTGTCAACTACTCCGGGTGT

isotig125       580 Q  I  T  W  A  K  F  C  K  E  N  M  A  G  K  G  F  S  F  W  
isotig125      1740 CAGATCACATGGGCTAAATTTTGCAAAGAAAACATGGCTGGCAAGGGCTTCTCCTTCTGG

isotig125       600 V  W  L  D  N  I  I  D  L  V  K  K  Y  I  L  A  L  W  N  E  
isotig125      1800 GTGTGGCTAGACAATATCATTGACCTTGTGAAAAAGTATATCTTGGCCCTCTGGAATGAA

isotig125       620 G  Y  I  M  G  F  I  S  K  E  R  E  R  A  I  L  S  T  K  P  
isotig125      1860 GGGTACATCATGGGCTTCATTAGCAAGGAGCGGGAGCGGGCGATCCTGAGCACGAAACCC

isotig125       640 P  G  T  F  L  L  R  F  S  E  S  S  K  E  G  G  V  T  F  T  
isotig125      1920 CCGGGCACCTTCCTGCTGAGATTCAGCGAGAGCAGCAAAGAAGGAGGGGTCACTTTCACT

isotig125       660 W  V  E  K  D  I  S  G  K  T  Q  I  Q  S  V  E  P  Y  T  K  
isotig125      1980 TGGGTGGAAAAGGACATCAGTGGCAAGACCCAGATCCAGTCTGTAGAGCCGTACACCAAG

isotig125       680 Q  Q  L  N  N  M  S  F  A  E  I  I  M  G  Y  K  I  M  D  A  
isotig125      2040 CAGCAGCTGAACAACATGTCCTTTGCTGAAATCATCATGGGCTACAAGATCATGGATGCC

isotig125       700 T  N  I  L  V  S  P  L  V  Y  L  Y  P  D  I  P  K  E  E  A  
isotig125      2100 ACCAACATCCTGGTGTCCCCATTGGTCTACCTCTACCCTGACATTCCCAAGGAGGAGGCG

isotig125       720 F  G  K  Y  C  R  P  E  S  Q  E  H  P  E  A  D  P  G  S  A  
isotig125      2160 TTCGGGAAGTACTGTCGACCAGAGAGCCAGGAGCATCCTGAAGCTGACCCCGGTAGTGCC

isotig125       740 A  P  Y  L  K  T  K  F  I  C  V  T  P  T  T  C  S  N  T  I  
isotig125      2220 GCCCCTTACCTGAAGACCAAGTTCATCTGTGTGACACCAACGACCTGCAGCAATACCATT

isotig125       760 D  L  P  M  S  P  P  H  F  R  F  I  D  A  V  W  K  R  R  C  
isotig125      2280 GACCTGCCGATGTCCCCCCCGCACTTTAGATTCATTGATGCAGTTTGGAAACGGAGGTGC

isotig125       780 A  L  G  R  R  A  V  V  T  H  V  H  G  S  D  F  G  V  R  Y  
isotig125      2340 GCCCTCGGCAGGAGGGCAGTTGTCACTCACGTTCATGGATCTGACTTCGGAGTGCGCTAC

isotig125       800 L  P  H  V  R  S    806
isotig125      2400 CTCCCCCATGTGAGGAGC 2418
""",
        )
        protein_record = protein_alignment.sequences[9]
        nucleotide_record = nucleotide_records[protein_record.id]
        self.assertEqual(nucleotide_record.id, protein_record.id)
        nucleotide_record = nucleotide_record.upper()
        alignments = aligner.align(protein_record, nucleotide_record)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 796], [0, 2388]]))
        )
        self.assertEqual(
            str(alignment),
            """\
isotig125         0 A  R  R  G  Q  A  A  L  G  S  P  A  A  R  T  W  S  Q  R  S  
isotig125         0 GCTAGGAGAGGCCAGGCGGCCCTCGGGAGCCCAGCTGCTCGCACCTGGAGCCAGCGCAGC

isotig125        20 P  A  S  R  A  S  A  R  E  T  V  T  P  P  D  C  G  R  M  A  
isotig125        60 CCGGCCAGTCGGGCCTCAGCCCGGGAGACAGTTACGCCCCCTGATTGCGGCAGGATGGCC

isotig125        40 Q  W  N  Q  L  Q  Q  L  D  T  R  Y  L  E  Q  L  H  Q  L  Y  
isotig125       120 CAGTGGAACCAGCTGCAGCAGCTGGACACTCGGTACCTGGAGCAGCTGCACCAGCTGTAC

isotig125        60 S  D  S  F  P  M  E  L  R  Q  F  L  A  P  W  I  E  S  Q  D  
isotig125       180 AGCGACAGCTTCCCCATGGAGCTGCGGCAGTTCCTGGCACCTTGGATTGAGAGTCAAGAT

isotig125        80 W  A  Y  A  A  S  K  E  S  H  A  T  L  V  F  H  N  L  L  G  
isotig125       240 TGGGCATATGCAGCCAGCAAAGAGTCACATGCCACACTGGTGTTTCATAATCTCTTGGGT

isotig125       100 E  I  D  Q  Q  Y  S  R  F  L  Q  E  S  N  V  L  Y  Q  H  N  
isotig125       300 GAGATTGACCAGCAGTACAGCCGATTCCTGCAAGAGTCCAACGTCCTCTATCAGCACAAC

isotig125       120 L  R  R  I  K  Q  F  L  Q  S  R  Y  L  E  K  P  M  E  I  A  
isotig125       360 CTTCGGAGGATCAAGCAGTTCCTACAGAGCAGGTATCTTGAGAAGCCGATGGAAATCGCC

isotig125       140 R  I  V  A  R  C  L  W  E  E  S  R  L  L  Q  T  A  A  T  A  
isotig125       420 CGCATCGTGGCCCGATGCCTGTGGGAAGAGTCTCGCCTCCTCCAGACGGCAGCCACTGCA

isotig125       160 A  Q  Q  G  G  Q  A  N  H  P  T  A  A  V  V  T  E  K  Q  Q  
isotig125       480 GCCCAGCAAGGAGGCCAGGCCAACCACCCAACAGCTGCTGTGGTGACGGAGAAACAGCAG

isotig125       180 M  L  E  Q  H  L  Q  D  V  R  K  R  V  Q  D  L  E  Q  K  M  
isotig125       540 ATGCTGGAGCAGCATCTTCAGGATGTCCGGAAACGTGTGCAGGATCTAGAACAGAAAATG

isotig125       200 K  V  V  E  N  L  Q  D  D  F  D  F  N  Y  K  T  L  K  S  Q  
isotig125       600 AAAGTGGTAGAGAATCTCCAGGATGACTTTGATTTCAACTATAAAACCCTCAAGAGTCAA

isotig125       220 G  D  M  Q  D  L  N  G  N  N  Q  S  V  T  R  Q  K  M  Q  Q  
isotig125       660 GGAGACATGCAGGATCTGAATGGAAACAACCAGTCTGTGACCAGGCAGAAGATGCAGCAG

isotig125       240 L  E  Q  M  L  T  A  L  D  Q  M  R  R  S  I  V  S  E  L  A  
isotig125       720 CTGGAACAGATGCTCACGGCGCTGGACCAGATGCGGAGAAGCATTGTGAGTGAGCTGGCG

isotig125       260 G  L  L  S  A  M  E  Y  V  Q  K  T  L  T  D  E  E  L  A  D  
isotig125       780 GGGCTTTTGTCGGCAATGGAGTACGTGCAGAAAACACTCACAGACGAGGAGCTGGCTGAC

isotig125       280 W  K  R  R  Q  Q  I  A  C  I  G  G  P  P  N  I  C  L  D  R  
isotig125       840 TGGAAGAGGCGGCAGCAGATCGCGTGCATTGGAGGCCCTCCCAACATCTGCCTGGATCGC

isotig125       300 L  E  N  W  I  T  S  L  A  E  S  Q  L  Q  T  R  Q  Q  I  K  
isotig125       900 CTGGAAAACTGGATAACTTCGTTAGCAGAATCTCAACTTCAGACCCGCCAACAAATTAAG

isotig125       320 K  L  E  E  L  Q  Q  K  V  S  Y  K  G  D  P  I  V  Q  H  R  
isotig125       960 AAACTGGAGGAGCTACAGCAGAAGGTGTCCTACAAGGGGGACCCCATTGTGCAGCACCGG

isotig125       340 P  M  L  E  E  R  I  V  E  L  F  R  N  L  M  K  S  A  F  V  
isotig125      1020 CCGATGCTGGAGGAGCGGATCGTGGAGCTGTTCAGAAACTTGATGAAGAGTGCCTTCGTG

isotig125       360 V  E  R  Q  P  C  M  P  M  H  P  D  R  P  L  V  I  K  T  G  
isotig125      1080 GTGGAGCGACAGCCCTGCATGCCGATGCACCCCGACCGGCCCTTGGTCATCAAGACTGGT

isotig125       380 V  Q  F  T  T  K  V  R  L  L  V  K  F  P  E  L  N  Y  Q  L  
isotig125      1140 GTCCAGTTCACTACTAAAGTCAGGTTGTTGGTCAAGTTTCCCGAGTTGAATTATCAGCTT

isotig125       400 K  I  K  V  C  I  D  K  D  S  G  D  V  A  A  L  R  G  S  R  
isotig125      1200 AAAATTAAAGTGTGCATTGACAAAGATTCTGGGGACGTTGCTGCTCTCAGAGGATCTCGG

isotig125       420 K  F  N  I  L  G  T  N  T  K  V  M  N  M  E  E  S  N  N  G  
isotig125      1260 AAATTTAACATTCTGGGCACAAACACGAAGGTGATGAACATGGAAGAATCCAACAACGGC

isotig125       440 S  L  S  A  E  F  K  H  L  T  L  R  E  Q  R  C  G  N  G  G  
isotig125      1320 AGCCTGTCTGCGGAGTTCAAGCACTTGACCCTGAGGGAGCAGAGATGTGGGAATGGAGGC

isotig125       460 R  A  N  C  D  A  S  L  I  V  T  E  E  L  H  L  I  T  F  E  
isotig125      1380 CGTGCCAATTGTGATGCCTCCTTGATTGTGACCGAGGAGCTGCATCTGATCACCTTCGAG

isotig125       480 T  E  V  Y  H  Q  G  L  K  I  D  L  E  T  H  S  L  P  V  V  
isotig125      1440 ACTGAGGTGTACCACCAAGGCCTCAAGATTGACCTGGAGACCCATTCTTTGCCAGTTGTG

isotig125       500 V  I  S  N  I  C  Q  M  P  N  A  W  A  S  I  L  W  Y  N  M  
isotig125      1500 GTGATCTCCAACATCTGTCAGATGCCAAATGCCTGGGCATCCATCCTGTGGTATAACATG

isotig125       520 L  T  N  N  P  K  N  V  N  F  F  T  K  P  P  I  G  T  W  D  
isotig125      1560 CTGACCAACAACCCCAAGAACGTGAACTTCTTCACCAAGCCACCAATCGGAACCTGGGAC

isotig125       540 Q  V  A  E  V  L  S  W  Q  F  S  S  T  T  K  R  G  L  S  I  
isotig125      1620 CAGGTGGCCGAGGTGCTCAGCTGGCAGTTCTCATCCACCACAAAGCGAGGGCTGAGCATC

isotig125       560 E  Q  L  T  T  L  A  E  K  L  L  G  P  G  V  N  Y  S  G  C  
isotig125      1680 GAGCAGCTGACTACGCTGGCCGAGAAGCTCCTAGGACCTGGTGTCAACTACTCCGGGTGT

isotig125       580 Q  I  T  W  A  K  F  C  K  E  N  M  A  G  K  G  F  S  F  W  
isotig125      1740 CAGATCACATGGGCTAAATTTTGCAAAGAAAACATGGCTGGCAAGGGCTTCTCCTTCTGG

isotig125       600 V  W  L  D  N  I  I  D  L  V  K  K  Y  I  L  A  L  W  N  E  
isotig125      1800 GTGTGGCTAGACAATATCATTGACCTTGTGAAAAAGTATATCTTGGCCCTCTGGAATGAA

isotig125       620 G  Y  I  M  G  F  I  S  K  E  R  E  R  A  I  L  S  T  K  P  
isotig125      1860 GGGTACATCATGGGCTTCATTAGCAAGGAGCGGGAGCGGGCGATCCTGAGCACGAAACCC

isotig125       640 P  G  T  F  L  L  R  F  S  E  S  S  K  E  G  G  V  T  F  T  
isotig125      1920 CCGGGCACCTTCCTGCTGAGATTCAGCGAGAGCAGCAAAGAAGGAGGGGTCACTTTCACT

isotig125       660 W  V  E  K  D  I  S  G  K  T  Q  I  Q  S  V  E  P  Y  T  K  
isotig125      1980 TGGGTGGAAAAGGACATCAGTGGCAAGACCCAGATCCAGTCTGTAGAGCCGTACACCAAG

isotig125       680 Q  Q  L  N  N  M  S  F  A  E  I  I  M  G  Y  K  I  M  D  A  
isotig125      2040 CAGCAGCTGAACAACATGTCCTTTGCTGAAATCATCATGGGCTACAAGATCATGGATGCC

isotig125       700 T  N  I  L  V  S  P  L  V  Y  L  Y  P  D  I  P  K  E  E  A  
isotig125      2100 ACCAACATCCTGGTGTCCCCATTGGTCTACCTCTACCCTGACATTCCCAAGGAGGAGGCG

isotig125       720 F  G  K  Y  C  R  P  E  S  Q  E  H  P  E  A  D  P  G  S  C  
isotig125      2160 TTCGGGAAGTACTGTCGACCAGAGAGCCAGGAGCATCCTGAAGCTGACCCCGGTAGTTGT

isotig125       740 F  S  M  V  L  V  S  L  L  G  K  G  G  Q  C  R  S  L  E  E  
isotig125      2220 TTTTCCATGGTTCTGGTTTCGCTGTTAGGGAAAGGGGGACAGTGCAGGTCCTTGGAGGAG

isotig125       760 R  Q  G  H  D  R  V  S  G  G  E  S  C  Y  G  R  A  V  Y  W  
isotig125      2280 AGACAAGGACATGACCGGGTGTCTGGTGGTGAGTCCTGCTATGGAAGAGCTGTTTATTGG

isotig125       780 V  L  Q  G  D  R  D  S  R  E  D  Q  N  Q  A  S    796
isotig125      2340 GTACTTCAGGGTGACCGGGATTCAAGAGAAGACCAGAATCAGGCCTCA 2388
""",
        )
        alignment = protein_alignment.mapall(codon_alignments)
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[   0,    0,    0,    0,    0,    0,    0,    0,
                              0,    0,    0,    0,   36,   36,   63,   63,
                            213,  222,  240,  240,  330,  330,  330,  405,
                            405,  438,  450,  579,  579,  729,  729,  837,
                            837,  897,  912,  975,  987, 1020, 1026, 1074,
                           1077, 1107, 1107, 1116, 1119, 1182, 1188, 1302,
                           1302, 1392, 1398, 1455, 1458, 1470, 1695, 1695,
                           1824, 1824, 1887, 1887, 1896, 1902, 1956, 1980,
                           2022, 2025, 2043, 2046, 2061, 2106, 2187, 2220,
                           2232, 2253, 2271, 2331, 2334, 2352, 2364, 2373,
                           2382, 2424, 2463, 2484, 2517, 2541],
                          [   0,    0,    0,    0,    0,    0,    0,    0,
                              0,    0,    0,    0,   36,   36,   63,   63,
                            213,  222,  240,  240,  330,  330,  330,  405,
                            405,  438,  450,  579,  579,  729,  729,  837,
                            837,  897,  912,  975,  987, 1020, 1026, 1074,
                           1077, 1107, 1107, 1116, 1119, 1182, 1188, 1302,
                           1302, 1392, 1398, 1455, 1458, 1470, 1695, 1695,
                           1824, 1824, 1887, 1887, 1896, 1902, 1956, 1980,
                           2022, 2025, 2043, 2046, 2061, 2106, 2187, 2220,
                           2232, 2253, 2271, 2331, 2334, 2352, 2364, 2373,
                           2382, 2424, 2463, 2484, 2517, 2541],
                          [   0,    0,    0,    0,    0,    0,    0,    0,
                              0,    0,    0,    0,   36,   36,   63,   63,
                            213,  222,  240,  240,  330,  330,  330,  405,
                            405,  438,  450,  579,  579,  729,  729,  837,
                            837,  897,  912,  975,  987, 1020, 1026, 1074,
                           1077, 1107, 1107, 1116, 1119, 1182, 1188, 1302,
                           1302, 1392, 1398, 1455, 1458, 1470, 1695, 1695,
                           1824, 1824, 1887, 1887, 1896, 1902, 1956, 1980,
                           2022, 2025, 2043, 2046, 2061, 2106, 2187, 2220,
                           2232, 2253, 2271, 2331, 2334, 2352, 2364, 2373,
                           2382, 2424, 2463, 2484, 2517, 2541],
                          [   0,    0,    0,    0,    0,    0,    0,    0,
                              0,    0,    0,    0,   36,   36,   63,   63,
                            213,  222,  240,  240,  330,  330,  330,  405,
                            405,  438,  450,  579,  579,  729,  729,  837,
                            837,  897,  912,  975,  987, 1020, 1026, 1074,
                           1077, 1107, 1107, 1116, 1119, 1182, 1188, 1302,
                           1302, 1392, 1398, 1455, 1458, 1470, 1695, 1695,
                           1824, 1824, 1887, 1887, 1896, 1902, 1956, 1980,
                           2022, 2025, 2043, 2046, 2061, 2106, 2187, 2220,
                           2232, 2253, 2271, 2331, 2334, 2352, 2364, 2373,
                           2382, 2424, 2463, 2484, 2517, 2541],
                          [   0,    0,    0,    0,    0,    0,    0,    0,
                              0,    0,    0,    0,    0,    0,    0,    0,
                              0,    0,    0,    0,    0,    0,    0,   75,
                             75,  108,  120,  249,  249,  399,  399,  507,
                            507,  567,  582,  645,  657,  690,  696,  744,
                            747,  777,  777,  786,  789,  852,  858,  972,
                            972, 1062, 1068, 1125, 1128, 1140, 1365, 1365,
                           1494, 1494, 1557, 1557, 1566, 1572, 1626, 1650,
                           1692, 1695, 1713, 1716, 1731, 1776, 1857, 1890,
                           1902, 1923, 1941, 2001, 2004, 2022, 2034, 2043,
                           2052, 2094, 2133, 2154, 2187, 2211],
                          [   0,    0,    0,    0,    0,    0,    0,    0,
                              0,    0,    0,    0,    0,    0,    0,    0,
                              0,    0,    0,    0,    0,    0,    0,   75,
                             75,  108,  120,  249,  249,  399,  399,  507,
                            507,  567,  582,  645,  657,  690,  696,  744,
                            747,  777,  777,  786,  789,  852,  858,  972,
                            972, 1062, 1068, 1125, 1128, 1140, 1365, 1365,
                           1494, 1494, 1557, 1557, 1566, 1572, 1626, 1650,
                           1692, 1695, 1713, 1716, 1731, 1776, 1857, 1890,
                           1902, 1923, 1941, 2001, 2004, 2022, 2034, 2043,
                           2052, 2094, 2133, 2154, 2187, 2211],
                          [   0,    0,    0,    0,    0,    0,    0,    0,
                              0,    0,    0,    0,    0,    0,    0,    0,
                              0,    0,    0,    0,    0,    0,    0,   75,
                             75,  108,  120,  249,  249,  399,  399,  507,
                            507,  567,  582,  645,  657,  690,  696,  744,
                            747,  777,  777,  786,  789,  852,  858,  972,
                            972, 1062, 1068, 1125, 1128, 1140, 1365, 1365,
                           1494, 1494, 1557, 1557, 1566, 1572, 1626, 1650,
                           1692, 1695, 1713, 1716, 1731, 1776, 1857, 1890,
                           1902, 1923, 1941, 2001, 2004, 2022, 2034, 2043,
                           2052, 2094, 2133, 2154, 2187, 2211],
                          [   0,   12,   21,  357,  378,  417,  444,  474,
                            516,  549,  564,  570,  606,  615,  642,  645,
                            795,  795,  813,  828,  918,  918,  999, 1074,
                           1077, 1110, 1110, 1239, 1251, 1401, 1428, 1536,
                           1560, 1620, 1620, 1683, 1683, 1716, 1716, 1764,
                           1764, 1794, 1809, 1818, 1818, 1881, 1881, 1995,
                           2001, 2091, 2091, 2148, 2148, 2148, 2373, 2373,
                           2502, 2532, 2595, 2619, 2628, 2628, 2682, 2706,
                           2706, 2709, 2727, 2727, 2727, 2772, 2772, 2805,
                           2817, 2817, 2835, 2895, 2898, 2916, 2928, 2937,
                           2946, 2946, 2985, 3006, 3039, 3063],
                          [   0,    0,    9,    9,   30,   30,   57,   57,
                             99,   99,  114,  114,  150,  159,  186,  189,
                            339,  339,  357,  372,  462,  465,  546,  621,
                            624,  657,  669,  798,  810,  960,  987, 1095,
                           1119, 1179, 1179, 1242, 1254, 1287, 1287, 1335,
                           1335, 1365, 1380, 1389, 1392, 1455, 1455, 1569,
                           1575, 1665, 1665, 1722, 1725, 1725, 1950, 1953,
                           2082, 2112, 2175, 2175, 2184, 2184, 2184, 2208,
                           2208, 2208, 2226, 2226, 2226, 2271, 2271, 2304,
                           2304, 2304, 2304, 2304, 2307, 2325, 2337, 2346,
                           2346, 2346, 2385, 2385, 2418, 2418],
                          [   0,    0,    9,    9,   30,   30,   57,   57,
                             99,   99,  114,  114,  150,  159,  186,  189,
                            339,  339,  357,  372,  462,  465,  546,  621,
                            624,  657,  669,  798,  810,  960,  987, 1095,
                           1119, 1179, 1179, 1242, 1254, 1287, 1287, 1335,
                           1335, 1365, 1380, 1389, 1392, 1455, 1455, 1569,
                           1575, 1665, 1665, 1722, 1725, 1725, 1950, 1953,
                           2082, 2112, 2175, 2175, 2184, 2184, 2184, 2208,
                           2208, 2208, 2226, 2229, 2229, 2229, 2229, 2262,
                           2274, 2274, 2274, 2334, 2337, 2337, 2349, 2349,
                           2349, 2349, 2388, 2388, 2388, 2388]])
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "clustal"),
            """\
ENSG00000166888:ENST0000030013      --------------------------------------------------
ENSG00000166888:ENST0000054387      --------------------------------------------------
ENSG00000166888:ENST0000055615      --------------------------------------------------
ENSG00000166888:ENST0000045407      --------------------------------------------------
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
isotig12565                         ------------GCTAGGAGA-----------------------------
isotig12566                         ------------GCTAGGAGA-----------------------------

ENSG00000166888:ENST0000030013      --------------------------------------------------
ENSG00000166888:ENST0000054387      --------------------------------------------------
ENSG00000166888:ENST0000055615      --------------------------------------------------
ENSG00000166888:ENST0000045407      --------------------------------------------------
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         AAAAAAAAAAGAGTTGGCCCAAACCCCTGAGGGAGGCTGTCCTCGGCAGG
isotig12565                         --------------------------------------------------
isotig12566                         --------------------------------------------------

ENSG00000166888:ENST0000030013      --------------------------------------------------
ENSG00000166888:ENST0000054387      --------------------------------------------------
ENSG00000166888:ENST0000055615      --------------------------------------------------
ENSG00000166888:ENST0000045407      --------------------------------------------------
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         AGAGCAGACCTGAGGGAGGTTGGTGTGCCGGTCCTGGTTGTCCAGCCCTC
isotig12565                         --------------------------------------------------
isotig12566                         --------------------------------------------------

ENSG00000166888:ENST0000030013      --------------------------------------------------
ENSG00000166888:ENST0000054387      --------------------------------------------------
ENSG00000166888:ENST0000055615      --------------------------------------------------
ENSG00000166888:ENST0000045407      --------------------------------------------------
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         CCGTCCCATCCCAGGGCTGGCTCACCAGAAAAAAGTACAGGCATCCCAGG
isotig12565                         --------------------------------------------------
isotig12566                         --------------------------------------------------

ENSG00000166888:ENST0000030013      --------------------------------------------------
ENSG00000166888:ENST0000054387      --------------------------------------------------
ENSG00000166888:ENST0000055615      --------------------------------------------------
ENSG00000166888:ENST0000045407      --------------------------------------------------
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         CCCTCTGGGGCGGGAGCAGAGGGTCTCCTCTTGTCCGGGGGAAGACAAAA
isotig12565                         --------------------------------------------------
isotig12566                         --------------------------------------------------

ENSG00000166888:ENST0000030013      --------------------------------------------------
ENSG00000166888:ENST0000054387      --------------------------------------------------
ENSG00000166888:ENST0000055615      --------------------------------------------------
ENSG00000166888:ENST0000045407      --------------------------------------------------
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         ACAAAAACAAAAAACGAGTGTGTAGCCCTCACTGGAAGTCTTCTGGTGCT
isotig12565                         --------------------------------------------------
isotig12566                         --------------------------------------------------

ENSG00000166888:ENST0000030013      --------------------------------------------------
ENSG00000166888:ENST0000054387      --------------------------------------------------
ENSG00000166888:ENST0000055615      --------------------------------------------------
ENSG00000166888:ENST0000045407      --------------------------------------------------
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         CTGGGGCCGTTTGCACTTGGGAGCCTGGCTTCTGGGCATGGTGGCCGGGC
isotig12565                         --------------------------------------------------
isotig12566                         --------------------------------------------------

ENSG00000166888:ENST0000030013      --------------------------------------------------
ENSG00000166888:ENST0000054387      --------------------------------------------------
ENSG00000166888:ENST0000055615      --------------------------------------------------
ENSG00000166888:ENST0000045407      --------------------------------------------------
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         TCCGCCGGGGGCTTCAGCCCTCGGCAAAGCGTTTCTAGAACAGAGTGGGT
isotig12565                         -------GGCCAGGCGGCCCTCGGGAGC----------------------
isotig12566                         -------GGCCAGGCGGCCCTCGGGAGC----------------------

ENSG00000166888:ENST0000030013      --------------------------------------------------
ENSG00000166888:ENST0000054387      --------------------------------------------------
ENSG00000166888:ENST0000055615      --------------------------------------------------
ENSG00000166888:ENST0000045407      --------------------------------------------------
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         GGGAGCCGACTGAAGTCCCAGAACCGCGAACACTGACGCACAGGAAAGCC
isotig12565                         -----------------CCAGCTGCTCGCACCTGGAGCCAGCGC------
isotig12566                         -----------------CCAGCTGCTCGCACCTGGAGCCAGCGC------

ENSG00000166888:ENST0000030013      --------------------------------------------------
ENSG00000166888:ENST0000054387      --------------------------------------------------
ENSG00000166888:ENST0000055615      --------------------------------------------------
ENSG00000166888:ENST0000045407      --------------------------------------------------
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         TCGAGTGGGTGGAGAAATGCAAATCCCCCAGTGGCATGGCCTCAGCCGCT
isotig12565                         ------------------------AGCCCGGCCAGTCGGGCCTCAGCCCG
isotig12566                         ------------------------AGCCCGGCCAGTCGGGCCTCAGCCCG

ENSG00000166888:ENST0000030013      --------------------------------------------------
ENSG00000166888:ENST0000054387      --------------------------------------------------
ENSG00000166888:ENST0000055615      --------------------------------------------------
ENSG00000166888:ENST0000045407      --------------------------------------------------
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         GCGGACCCTGACGGCCGATTCTGACTTCGGACTTGAGCGGGACCAGCACC
isotig12565                         GGAGACAGTTACGCCC---------------------------------C
isotig12566                         GGAGACAGTTACGCCC---------------------------------C

ENSG00000166888:ENST0000030013      --------------------ATGTCTCTGTGGGGTCTGGTCTCCAAGATG
ENSG00000166888:ENST0000054387      --------------------ATGTCTCTGTGGGGTCTGGTCTCCAAGATG
ENSG00000166888:ENST0000055615      --------------------ATGTCTCTGTGGGGTCTGGTCTCCAAGATG
ENSG00000166888:ENST0000045407      --------------------ATGTCTCTGTGGGGTCTGGTCTCCAAGATG
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         CGGACCCAGGGCGAGCCAGCATGTCTCAGTGGAATCAAGTCCAACAGTTA
isotig12565                         CTGATTGCGGCAGG------ATGGCCCAGTGGAACCAGCTGCAGCAGCTG
isotig12566                         CTGATTGCGGCAGG------ATGGCCCAGTGGAACCAGCTGCAGCAGCTG

ENSG00000166888:ENST0000030013      CCCCCA---------GAAAAAGTGCAGCGGCTCTATGTCGAC---TTTCC
ENSG00000166888:ENST0000054387      CCCCCA---------GAAAAAGTGCAGCGGCTCTATGTCGAC---TTTCC
ENSG00000166888:ENST0000055615      CCCCCA---------GAAAAAGTGCAGCGGCTCTATGTCGAC---TTTCC
ENSG00000166888:ENST0000045407      CCCCCA---------GAAAAAGTGCAGCGGCTCTATGTCGAC---TTTCC
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         GAAATCAAGTTTTTGGAGCAAGTTGATCAGTTCTATGATGACAACTTTCC
isotig12565                         GACACTCGGTACCTGGAGCAGCTGCACCAGCTGTACAGCGACAGCTTCCC
isotig12566                         GACACTCGGTACCTGGAGCAGCTGCACCAGCTGTACAGCGACAGCTTCCC

ENSG00000166888:ENST0000030013      CCAACACCTGCGGCATCTTCTGGGTGACTGGCTGGAGAGCCAGCCCTGGG
ENSG00000166888:ENST0000054387      CCAACACCTGCGGCATCTTCTGGGTGACTGGCTGGAGAGCCAGCCCTGGG
ENSG00000166888:ENST0000055615      CCAACACCTGCGGCATCTTCTGGGTGACTGGCTGGAGAGCCAGCCCTGGG
ENSG00000166888:ENST0000045407      CCAACACCTGCGGCATCTTCTGGGTGACTGGCTGGAGAGCCAGCCCTGGG
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         CATGGAAATCCGACATCTGCTGGCCCAGTGGATTGAGCATCAAGACTGGG
isotig12565                         CATGGAGCTGCGGCAGTTCCTGGCACCTTGGATTGAGAGTCAAGATTGGG
isotig12566                         CATGGAGCTGCGGCAGTTCCTGGCACCTTGGATTGAGAGTCAAGATTGGG

ENSG00000166888:ENST0000030013      AGTTCCTGGTCGGCTCCGACGCCTTCTGCTGCAACTTGGCTAGTGCCCTA
ENSG00000166888:ENST0000054387      AGTTCCTGGTCGGCTCCGACGCCTTCTGCTGCAACTTGGCTAGTGCCCTA
ENSG00000166888:ENST0000055615      AGTTCCTGGTCGGCTCCGACGCCTTCTGCTGCAACTTGGCTAGTGCCCTA
ENSG00000166888:ENST0000045407      AGTTCCTGGTCGGCTCCGACGCCTTCTGCTGCAACTTGGCTAGTGCCCTA
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         AGGTGGCCTCTAACAATGAAACTATGGCAACAATTCTTCTTCAAAACTTA
isotig12565                         CATATGCAGCCAGCAAAGAGTCACATGCCACACTGGTGTTTCATAATCTC
isotig12566                         CATATGCAGCCAGCAAAGAGTCACATGCCACACTGGTGTTTCATAATCTC

ENSG00000166888:ENST0000030013      CTTTCAGACACTGTCCAGCACCTTCAGGCCTCGGTGGGAGAGCAGGGGGA
ENSG00000166888:ENST0000054387      CTTTCAGACACTGTCCAGCACCTTCAGGCCTCGGTGGGAGAGCAGGGGGA
ENSG00000166888:ENST0000055615      CTTTCAGACACTGTCCAGCACCTTCAGGCCTCGGTGGGAGAGCAGGGGGA
ENSG00000166888:ENST0000045407      CTTTCAGACACTGTCCAGCACCTTCAGGCCTCGGTGGGAGAGCAGGGGGA
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         TTAATACAATTGGATGAACAGTTAGGTCGTGTTTCCAAAGAGAAA-----
isotig12565                         TTGGGTGAGATTGACCAGCAGTACAGCCGATTCCTGCAAGAGTCC-----
isotig12566                         TTGGGTGAGATTGACCAGCAGTACAGCCGATTCCTGCAAGAGTCC-----

ENSG00000166888:ENST0000030013      GGGGAGCACCATCTTGCAACAC---------------ATCAGCACCCTTG
ENSG00000166888:ENST0000054387      GGGGAGCACCATCTTGCAACAC---------------ATCAGCACCCTTG
ENSG00000166888:ENST0000055615      GGGGAGCACCATCTTGCAACAC---------------ATCAGCACCCTTG
ENSG00000166888:ENST0000045407      GGGGAGCACCATCTTGCAACAC---------------ATCAGCACCCTTG
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         ----AACCTGCTATTGATCCACAATCTAAAGAGAATTAGAAAAGTACTTC
isotig12565                         ----AACGTCCTCTATCAGCACAACCTTCGGAGGATCAAGCAGTTCCTAC
isotig12566                         ----AACGTCCTCTATCAGCACAACCTTCGGAGGATCAAGCAGTTCCTAC

ENSG00000166888:ENST0000030013      AGAGCATATATCAGAGGGACCCCCTGAAGCTGGTGGCCACTTTCAGACAA
ENSG00000166888:ENST0000054387      AGAGCATATATCAGAGGGACCCCCTGAAGCTGGTGGCCACTTTCAGACAA
ENSG00000166888:ENST0000055615      AGAGCATATATCAGAGGGACCCCCTGAAGCTGGTGGCCACTTTCAGACAA
ENSG00000166888:ENST0000045407      AGAGCATATATCAGAGGGACCCCCTGAAGCTGGTGGCCACTTTCAGACAA
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         AGGGGAAGTTTCATGGAAATCCAATGCATGTAGCCGTGGTAATCTCAAAT
isotig12565                         AGAGCAGGTATCTTGAGAAGCCGATGGAAATCGCCCGCATCGTGGCCCGA
isotig12566                         AGAGCAGGTATCTTGAGAAGCCGATGGAAATCGCCCGCATCGTGGCCCGA

ENSG00000166888:ENST0000030013      ATACTTCAAGGAGAGAAAAAAGCTGTT-----------------------
ENSG00000166888:ENST0000054387      ATACTTCAAGGAGAGAAAAAAGCTGTT-----------------------
ENSG00000166888:ENST0000055615      ATACTTCAAGGAGAGAAAAAAGCTGTT-----------------------
ENSG00000166888:ENST0000045407      ATACTTCAAGGAGAGAAAAAAGCTGTT-----------------------
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         TGTTTAAGGGAAGAGAGGAGAATACTG---GCTGCAGCGAACATGCCTAT
isotig12565                         TGCCTGTGGGAAGAGTCTCGCCTCCTCCAGACGGCAGCCACTGCAGCCCA
isotig12566                         TGCCTGTGGGAAGAGTCTCGCCTCCTCCAGACGGCAGCCACTGCAGCCCA

ENSG00000166888:ENST0000030013      --------------------------------------------------
ENSG00000166888:ENST0000054387      --------------------------------------------------
ENSG00000166888:ENST0000055615      --------------------------------------------------
ENSG00000166888:ENST0000045407      --------------------------------------------------
ENSG00000166888:ENST0000053891      --------------------------------------------------
ENSG00000166888:ENST0000053721      --------------------------------------------------
ENSG00000166888:ENST0000053520      --------------------------------------------------
isotig46679                         CCAGGGACCTCTGGAGAAATCCTTACAAAGTTCGTCGGTTTCAGAAAGAC
isotig12565                         GCAAGGAGGCCAGGCCAACCACCCAACAGCTGCTGTGGTGACGGAGAAAC
isotig12566                         GCAAGGAGGCCAGGCCAACCACCCAACAGCTGCTGTGGTGACGGAGAAAC

ENSG00000166888:ENST0000030013      -----------ATGGAACAGTTCCGCCACTTGCCAATGCCTTTCCACTGG
ENSG00000166888:ENST0000054387      -----------ATGGAACAGTTCCGCCACTTGCCAATGCCTTTCCACTGG
ENSG00000166888:ENST0000055615      -----------ATGGAACAGTTCCGCCACTTGCCAATGCCTTTCCACTGG
ENSG00000166888:ENST0000045407      -----------ATGGAACAGTTCCGCCACTTGCCAATGCCTTTCCACTGG
ENSG00000166888:ENST0000053891      -----------ATGGAACAGTTCCGCCACTTGCCAATGCCTTTCCACTGG
ENSG00000166888:ENST0000053721      -----------ATGGAACAGTTCCGCCACTTGCCAATGCCTTTCCACTGG
ENSG00000166888:ENST0000053520      -----------ATGGAACAGTTCCGCCACTTGCCAATGCCTTTCCACTGG
isotig46679                         AGAGAAATGTGGAACACAAAGTGGCTGCCATTAAAAACAGTGTGCAGATG
isotig12565                         AGCAGATGCTGGAGCAGCATCTTCAGGATGTCCGGAAACGTGTGCAGGAT
isotig12566                         AGCAGATGCTGGAGCAGCATCTTCAGGATGTCCGGAAACGTGTGCAGGAT

ENSG00000166888:ENST0000030013      AAGCAGGAAGAACTCAAGTTTAAGACAGGCTTGCGG---AGGCTGCAGCA
ENSG00000166888:ENST0000054387      AAGCAGGAAGAACTCAAGTTTAAGACAGGCTTGCGG---AGGCTGCAGCA
ENSG00000166888:ENST0000055615      AAGCAGGAAGAACTCAAGTTTAAGACAGGCTTGCGG---AGGCTGCAGCA
ENSG00000166888:ENST0000045407      AAGCAGGAAGAACTCAAGTTTAAGACAGGCTTGCGG---AGGCTGCAGCA
ENSG00000166888:ENST0000053891      AAGCAGGAAGAACTCAAGTTTAAGACAGGCTTGCGG---AGGCTGCAGCA
ENSG00000166888:ENST0000053721      AAGCAGGAAGAACTCAAGTTTAAGACAGGCTTGCGG---AGGCTGCAGCA
ENSG00000166888:ENST0000053520      AAGCAGGAAGAACTCAAGTTTAAGACAGGCTTGCGG---AGGCTGCAGCA
isotig46679                         ACAGAACAAGACACCAAATACTTGGAAGATCTGCAAGATGAATTTGACTA
isotig12565                         CTAGAACAGAAAATGAAAGTGGTAGAGAATCTCCAGGATGACTTTGATTT
isotig12566                         CTAGAACAGAAAATGAAAGTGGTAGAGAATCTCCAGGATGACTTTGATTT

ENSG00000166888:ENST0000030013      CCGAGTAGGGGAGATCCACCTTCTCCGAGAAGCCCTGCAGAAGGGGGCTG
ENSG00000166888:ENST0000054387      CCGAGTAGGGGAGATCCACCTTCTCCGAGAAGCCCTGCAGAAGGGGGCTG
ENSG00000166888:ENST0000055615      CCGAGTAGGGGAGATCCACCTTCTCCGAGAAGCCCTGCAGAAGGGGGCTG
ENSG00000166888:ENST0000045407      CCGAGTAGGGGAGATCCACCTTCTCCGAGAAGCCCTGCAGAAGGGGGCTG
ENSG00000166888:ENST0000053891      CCGAGTAGGGGAGATCCACCTTCTCCGAGAAGCCCTGCAGAAGGGGGCTG
ENSG00000166888:ENST0000053721      CCGAGTAGGGGAGATCCACCTTCTCCGAGAAGCCCTGCAGAAGGGGGCTG
ENSG00000166888:ENST0000053520      CCGAGTAGGGGAGATCCACCTTCTCCGAGAAGCCCTGCAGAAGGGGGCTG
isotig46679                         CAGGTATAAAACAATTCAGACA------------ATGGACCAGGGTGACA
isotig12565                         CAACTATAAAACCCTCAAGAGTCAAGGAGACATGCAGGATCTGAATGGAA
isotig12566                         CAACTATAAAACCCTCAAGAGTCAAGGAGACATGCAGGATCTGAATGGAA

ENSG00000166888:ENST0000030013      AGGCTGGCCAAGTGTCTCTGCACAGCTTGATAGAAACTCCTGCTAATGGG
ENSG00000166888:ENST0000054387      AGGCTGGCCAAGTGTCTCTGCACAGCTTGATAGAAACTCCTGCTAATGGG
ENSG00000166888:ENST0000055615      AGGCTGGCCAAGTGTCTCTGCACAGCTTGATAGAAACTCCTGCTAATGGG
ENSG00000166888:ENST0000045407      AGGCTGGCCAAGTGTCTCTGCACAGCTTGATAGAAACTCCTGCTAATGGG
ENSG00000166888:ENST0000053891      AGGCTGGCCAAGTGTCTCTGCACAGCTTGATAGAAACTCCTGCTAATGGG
ENSG00000166888:ENST0000053721      AGGCTGGCCAAGTGTCTCTGCACAGCTTGATAGAAACTCCTGCTAATGGG
ENSG00000166888:ENST0000053520      AGGCTGGCCAAGTGTCTCTGCACAGCTTGATAGAAACTCCTGCTAATGGG
isotig46679                         AGAATAGCATCCTAATGAACCAGGAGGTTTTGACACTCCAAGAAATGCTT
isotig12565                         ACAACCAGTCTGTGACCAGGCAGAAGATGCAGCAGCTGGAACAGATGCTC
isotig12566                         ACAACCAGTCTGTGACCAGGCAGAAGATGCAGCAGCTGGAACAGATGCTC

ENSG00000166888:ENST0000030013      ACTGGGCCAAGTGAGGCCCTGGCCATGCTACTGCAGGAGACCACTGGAGA
ENSG00000166888:ENST0000054387      ACTGGGCCAAGTGAGGCCCTGGCCATGCTACTGCAGGAGACCACTGGAGA
ENSG00000166888:ENST0000055615      ACTGGGCCAAGTGAGGCCCTGGCCATGCTACTGCAGGAGACCACTGGAGA
ENSG00000166888:ENST0000045407      ACTGGGCCAAGTGAGGCCCTGGCCATGCTACTGCAGGAGACCACTGGAGA
ENSG00000166888:ENST0000053891      ACTGGGCCAAGTGAGGCCCTGGCCATGCTACTGCAGGAGACCACTGGAGA
ENSG00000166888:ENST0000053721      ACTGGGCCAAGTGAGGCCCTGGCCATGCTACTGCAGGAGACCACTGGAGA
ENSG00000166888:ENST0000053520      ACTGGGCCAAGTGAGGCCCTGGCCATGCTACTGCAGGAGACCACTGGAGA
isotig46679                         AATAGCCTGGACTTCAAGAGAAAGGAAGCACTCACTAAGATGACACAGAT
isotig12565                         ACGGCGCTGGACCAGATGCGGAGAAGCATTGTGAGTGAGCTGGCGGGGCT
isotig12566                         ACGGCGCTGGACCAGATGCGGAGAAGCATTGTGAGTGAGCTGGCGGGGCT

ENSG00000166888:ENST0000030013      GCTAGAGGCAGCC------------AAAGCCCTAGTGCTGAAGAGGATCC
ENSG00000166888:ENST0000054387      GCTAGAGGCAGCC------------AAAGCCCTAGTGCTGAAGAGGATCC
ENSG00000166888:ENST0000055615      GCTAGAGGCAGCC------------AAAGCCCTAGTGCTGAAGAGGATCC
ENSG00000166888:ENST0000045407      GCTAGAGGCAGCC------------AAAGCCCTAGTGCTGAAGAGGATCC
ENSG00000166888:ENST0000053891      GCTAGAGGCAGCC------------AAAGCCCTAGTGCTGAAGAGGATCC
ENSG00000166888:ENST0000053721      GCTAGAGGCAGCC------------AAAGCCCTAGTGCTGAAGAGGATCC
ENSG00000166888:ENST0000053520      GCTAGAGGCAGCC------------AAAGCCCTAGTGCTGAAGAGGATCC
isotig46679                         AGTGAACGAGTCGGACCTGCTGATGAGCAGCATGCTCATAGAAGAGCTGC
isotig12565                         TTTGTCGGCAATGGAGTACGTGCAGAAAACACTCACAGACGAGGAGCTGG
isotig12566                         TTTGTCGGCAATGGAGTACGTGCAGAAAACACTCACAGACGAGGAGCTGG

ENSG00000166888:ENST0000030013      AGATTTGGAAACGGCAGCAGCAGCTGGCAGGGAATGGCGCACCGTTTGAG
ENSG00000166888:ENST0000054387      AGATTTGGAAACGGCAGCAGCAGCTGGCAGGGAATGGCGCACCGTTTGAG
ENSG00000166888:ENST0000055615      AGATTTGGAAACGGCAGCAGCAGCTGGCAGGGAATGGCGCACCGTTTGAG
ENSG00000166888:ENST0000045407      AGATTTGGAAACGGCAGCAGCAGCTGGCAGGGAATGGCGCACCGTTTGAG
ENSG00000166888:ENST0000053891      AGATTTGGAAACGGCAGCAGCAGCTGGCAGGGAATGGCGCACCGTTTGAG
ENSG00000166888:ENST0000053721      AGATTTGGAAACGGCAGCAGCAGCTGGCAGGGAATGGCGCACCGTTTGAG
ENSG00000166888:ENST0000053520      AGATTTGGAAACGGCAGCAGCAGCTGGCAGGGAATGGCGCACCGTTTGAG
isotig46679                         AGGACTGGAAGAGGAGGCAGCAGATCGCCTGCATCGGTGGCCCACTCCAC
isotig12565                         CTGACTGGAAGAGGCGGCAGCAGATCGCGTGCATTGGAGGCCCTCCCAAC
isotig12566                         CTGACTGGAAGAGGCGGCAGCAGATCGCGTGCATTGGAGGCCCTCCCAAC

ENSG00000166888:ENST0000030013      GAGAGCCTGGCCCCACTCCAGGAGAGGTGTGAAAGCCTGGTGGACATTTA
ENSG00000166888:ENST0000054387      GAGAGCCTGGCCCCACTCCAGGAGAGGTGTGAAAGCCTGGTGGACATTTA
ENSG00000166888:ENST0000055615      GAGAGCCTGGCCCCACTCCAGGAGAGGTGTGAAAGCCTGGTGGACATTTA
ENSG00000166888:ENST0000045407      GAGAGCCTGGCCCCACTCCAGGAGAGGTGTGAAAGCCTGGTGGACATTTA
ENSG00000166888:ENST0000053891      GAGAGCCTGGCCCCACTCCAGGAGAGGTGTGAAAGCCTGGTGGACATTTA
ENSG00000166888:ENST0000053721      GAGAGCCTGGCCCCACTCCAGGAGAGGTGTGAAAGCCTGGTGGACATTTA
ENSG00000166888:ENST0000053520      GAGAGCCTGGCCCCACTCCAGGAGAGGTGTGAAAGCCTGGTGGACATTTA
isotig46679                         AACGGGCTGGACCAGCTTCAGAACTGCTTTACCCTGTTGGCAGAAAGTCT
isotig12565                         ATCTGCCTGGATCGCCTGGAAAACTGGATAACTTCGTTAGCAGAATCTCA
isotig12566                         ATCTGCCTGGATCGCCTGGAAAACTGGATAACTTCGTTAGCAGAATCTCA

ENSG00000166888:ENST0000030013      TTCCCAGCTACAGCAGGAGGTAGGG-------------------------
ENSG00000166888:ENST0000054387      TTCCCAGCTACAGCAGGAGGTAGGG-------------------------
ENSG00000166888:ENST0000055615      TTCCCAGCTACAGCAGGAGGTAGGG-------------------------
ENSG00000166888:ENST0000045407      TTCCCAGCTACAGCAGGAGGTAGGG-------------------------
ENSG00000166888:ENST0000053891      TTCCCAGCTACAGCAGGAGGTAGGG-------------------------
ENSG00000166888:ENST0000053721      TTCCCAGCTACAGCAGGAGGTAGGG-------------------------
ENSG00000166888:ENST0000053520      TTCCCAGCTACAGCAGGAGGTAGGG-------------------------
isotig46679                         TTTCCAACTCAGACGACAGCTGGAGAAATTAGAGGAGCAGTCTTCCAAGA
isotig12565                         ACTTCAGACCCGCCAACAAATTAAGAAACTGGAGGAGCTACAGCAGAAGG
isotig12566                         ACTTCAGACCCGCCAACAAATTAAGAAACTGGAGGAGCTACAGCAGAAGG

ENSG00000166888:ENST0000030013      --GCGGCTGGTGGGGAGCTTGAGCCCAAGACCCGGGCATCGCTGACTGGC
ENSG00000166888:ENST0000054387      --GCGGCTGGTGGGGAGCTTGAGCCCAAGACCCGGGCATCGCTGACTGGC
ENSG00000166888:ENST0000055615      --GCGGCTGGTGGGGAGCTTGAGCCCAAGACCCGGGCATCGCTGACTGGC
ENSG00000166888:ENST0000045407      --GCGGCTGGTGGGGAGCTTGAGCCCAAGACCCGGGCATCGCTGACTGGC
ENSG00000166888:ENST0000053891      --GCGGCTGGTGGGGAGCTTGAGCCCAAGACCCGGGCATCGCTGACTGGC
ENSG00000166888:ENST0000053721      --GCGGCTGGTGGGGAGCTTGAGCCCAAGACCCGGGCATCGCTGACTGGC
ENSG00000166888:ENST0000053520      --GCGGCTGGTGGGGAGCTTGAGCCCAAGACCCGGGCATCGCTGACTGGC
isotig46679                         TGACTTACGAAGGAGACCCCATCCCCACGCAGAGAGCACACCTGCTGGAG
isotig12565                         TGTCCTACAAGGGGGACCCCATTGTGCAGCACCGGCCGATGCTGGAGGAG
isotig12566                         TGTCCTACAAGGGGGACCCCATTGTGCAGCACCGGCCGATGCTGGAGGAG

ENSG00000166888:ENST0000030013      CGGCTGGATGAAGTCCTGAGAACCCTCGTCACCAGTTGCTTCCTGGTGGA
ENSG00000166888:ENST0000054387      CGGCTGGATGAAGTCCTGAGAACCCTCGTCACCAGTTGCTTCCTGGTGGA
ENSG00000166888:ENST0000055615      CGGCTGGATGAAGTCCTGAGAACCCTCGTCACCAGTTGCTTCCTGGTGGA
ENSG00000166888:ENST0000045407      CGGCTGGATGAAGTCCTGAGAACCCTCGTCACCAGTTGCTTCCTGGTGGA
ENSG00000166888:ENST0000053891      CGGCTGGATGAAGTCCTGAGAACCCTCGTCACCAGTTGCTTCCTGGTGGA
ENSG00000166888:ENST0000053721      CGGCTGGATGAAGTCCTGAGAACCCTCGTCACCAGTTGCTTCCTGGTGGA
ENSG00000166888:ENST0000053520      CGGCTGGATGAAGTCCTGAGAACCCTCGTCACCAGTTGCTTCCTGGTGGA
isotig46679                         AGAGCCACCTTCCTGATCTACAACCTTTTCAAGAACTCATTTGTGGTTGA
isotig12565                         CGGATCGTGGAGCTGTTCAGAAACTTGATGAAGAGTGCCTTCGTGGTGGA
isotig12566                         CGGATCGTGGAGCTGTTCAGAAACTTGATGAAGAGTGCCTTCGTGGTGGA

ENSG00000166888:ENST0000030013      GAAGCAGCCC------------------------CCCCAGGTACTGAAGA
ENSG00000166888:ENST0000054387      GAAGCAGCCC------------------------CCCCAGGTACTGAAGA
ENSG00000166888:ENST0000055615      GAAGCAGCCC------------------------CCCCAGGTACTGAAGA
ENSG00000166888:ENST0000045407      GAAGCAGCCC------------------------CCCCAGGTACTGAAGA
ENSG00000166888:ENST0000053891      GAAGCAGCCC------------------------CCCCAGGTACTGAAGA
ENSG00000166888:ENST0000053721      GAAGCAGCCC------------------------CCCCAGGTACTGAAGA
ENSG00000166888:ENST0000053520      GAAGCAGCCC------------------------CCCCAGGTACTGAAGA
isotig46679                         GCGACAGCCCTGCATGCCAACACACCCTCAGAGGCCGCTGGTACTCAAAA
isotig12565                         GCGACAGCCCTGCATGCCGATGCACCCCGACCGGCCCTTGGTCATCAAGA
isotig12566                         GCGACAGCCCTGCATGCCGATGCACCCCGACCGGCCCTTGGTCATCAAGA

ENSG00000166888:ENST0000030013      CTCAGACCAAGTTCCAGGCTGGAGTTCGATTCCTGTTGGGCTTGAGGTTC
ENSG00000166888:ENST0000054387      CTCAGACCAAGTTCCAGGCTGGAGTTCGATTCCTGTTGGGCTTGAGGTTC
ENSG00000166888:ENST0000055615      CTCAGACCAAGTTCCAGGCTGGAGTTCGATTCCTGTTGGGCTTGAGGTTC
ENSG00000166888:ENST0000045407      CTCAGACCAAGTTCCAGGCTGGAGTTCGATTCCTGTTGGGCTTGAGGTTC
ENSG00000166888:ENST0000053891      CTCAGACCAAGTTCCAGGCTGGAGTTCGATTCCTGTTGGGCTTGAGGTTC
ENSG00000166888:ENST0000053721      CTCAGACCAAGTTCCAGGCTGGAGTTCGATTCCTGTTGGGCTTGAGGTTC
ENSG00000166888:ENST0000053520      CTCAGACCAAGTTCCAGGCTGGAGTTCGATTCCTGTTGGGCTTGAGGTTC
isotig46679                         CCCTCATTCAGTTCACCGCGAAACTGAGACTACTAATAAAATTG------
isotig12565                         CTGGTGTCCAGTTCACTACTAAAGTCAGGTTGTTGGTCAAGTTT------
isotig12566                         CTGGTGTCCAGTTCACTACTAAAGTCAGGTTGTTGGTCAAGTTT------

ENSG00000166888:ENST0000030013      CTGGGGGCCCCAGCCAAGCCTCCGCTGGTCAGGGCCGACATGGTGACAGA
ENSG00000166888:ENST0000054387      CTGGGGGCCCCAGCCAAGCCTCCGCTGGTCAGGGCCGACATGGTGACAGA
ENSG00000166888:ENST0000055615      CTGGGGGCCCCAGCCAAGCCTCCGCTGGTCAGGGCCGACATGGTGACAGA
ENSG00000166888:ENST0000045407      CTGGGGGCCCCAGCCAAGCCTCCGCTGGTCAGGGCCGACATGGTGACAGA
ENSG00000166888:ENST0000053891      CTGGGGGCCCCAGCCAAGCCTCCGCTGGTCAGGGCCGACATGGTGACAGA
ENSG00000166888:ENST0000053721      CTGGGGGCCCCAGCCAAGCCTCCGCTGGTCAGGGCCGACATGGTGACAGA
ENSG00000166888:ENST0000053520      CTGGGGGCCCCAGCCAAGCCTCCGCTGGTCAGGGCCGACATGGTGACAGA
isotig46679                         ---------CCGGAACTCAACTATCAGGTGAAAGTAAAGGCATCGATCGA
isotig12565                         ---------CCCGAGTTGAATTATCAGCTTAAAATTAAAGTGTGCATTGA
isotig12566                         ---------CCCGAGTTGAATTATCAGCTTAAAATTAAAGTGTGCATTGA

ENSG00000166888:ENST0000030013      GAAGCAGGCGCGGGAGCTGAGTGTGCCTCAGGGTCCTGGGGCTGGAGCAG
ENSG00000166888:ENST0000054387      GAAGCAGGCGCGGGAGCTGAGTGTGCCTCAGGGTCCTGGGGCTGGAGCAG
ENSG00000166888:ENST0000055615      GAAGCAGGCGCGGGAGCTGAGTGTGCCTCAGGGTCCTGGGGCTGGAGCAG
ENSG00000166888:ENST0000045407      GAAGCAGGCGCGGGAGCTGAGTGTGCCTCAGGGTCCTGGGGCTGGAGCAG
ENSG00000166888:ENST0000053891      GAAGCAGGCGCGGGAGCTGAGTGTGCCTCAGGGTCCTGGGGCTGGAGCAG
ENSG00000166888:ENST0000053721      GAAGCAGGCGCGGGAGCTGAGTGTGCCTCAGGGTCCTGGGGCTGGAGCAG
ENSG00000166888:ENST0000053520      GAAGCAGGCGCGGGAGCTGAGTGTGCCTCAGGGTCCTGGGGCTGGAGCAG
isotig46679                         CAAGAATGTTTCAACGCTAAGC------------AATAGAAGATTTGTGC
isotig12565                         CAAAGATTCTGGGGACGTTGCTGCTCTCAGAGGATCTCGGAAATTTAACA
isotig12566                         CAAAGATTCTGGGGACGTTGCTGCTCTCAGAGGATCTCGGAAATTTAACA

ENSG00000166888:ENST0000030013      AAAGCACTGGAGAAATCATCAACAACACTGTGCCCTTGGAGAACAGCATT
ENSG00000166888:ENST0000054387      AAAGCACTGGAGAAATCATCAACAACACTGTGCCCTTGGAGAACAGCATT
ENSG00000166888:ENST0000055615      AAAGCACTGGAGAAATCATCAACAACACTGTGCCCTTGGAGAACAGCATT
ENSG00000166888:ENST0000045407      AAAGCACTGGAGAAATCATCAACAACACTGTGCCCTTGGAGAACAGCATT
ENSG00000166888:ENST0000053891      AAAGCACTGGAGAAATCATCAACAACACTGTGCCCTTGGAGAACAGCATT
ENSG00000166888:ENST0000053721      AAAGCACTGGAGAAATCATCAACAACACTGTGCCCTTGGAGAACAGCATT
ENSG00000166888:ENST0000053520      AAAGCACTGGAGAAATCATCAACAACACTGTGCCCTTGGAGAACAGCATT
isotig46679                         TTTGTGGAACTCAAGTC------AAAGCCATGTCCATCGAGGAATCCTCC
isotig12565                         TTCTGGGCACAAACACG------AAGGTGATGAACATGGAAGAATCCAAC
isotig12566                         TTCTGGGCACAAACACG------AAGGTGATGAACATGGAAGAATCCAAC

ENSG00000166888:ENST0000030013      CCTGGGAACTGCTGCTCTGCCCTGTTCAAGAACCTGCTTCTCAAGAAGAT
ENSG00000166888:ENST0000054387      CCTGGGAACTGCTGCTCTGCCCTGTTCAAGAACCTGCTTCTCAAGAAGAT
ENSG00000166888:ENST0000055615      CCTGGGAACTGCTGCTCTGCCCTGTTCAAGAACCTGCTTCTCAAGAAGAT
ENSG00000166888:ENST0000045407      CCTGGGAACTGCTGCTCTGCCCTGTTCAAGAACCTGCTTCTCAAGAAGAT
ENSG00000166888:ENST0000053891      CCTGGGAACTGCTGCTCTGCCCTGTTCAAGAACCTGCTTCTCAAGAAGAT
ENSG00000166888:ENST0000053721      CCTGGGAACTGCTGCTCTGCCCTGTTCAAGAACCTGCTTCTCAAGAAGAT
ENSG00000166888:ENST0000053520      CCTGGGAACTGCTGCTCTGCCCTGTTCAAGAACCTGCTTCTCAAGAAGAT
isotig46679                         AATGGGAGCCTCTCAGTAGAA---TTTAGACATTTGCAACCGAAGGAAAT
isotig12565                         AACGGCAGCCTGTCTGCGGAG---TTCAAGCACTTGACCCTGAGGGAGCA
isotig12566                         AACGGCAGCCTGTCTGCGGAG---TTCAAGCACTTGACCCTGAGGGAGCA

ENSG00000166888:ENST0000030013      CAAG---------------CGGTGTGAGCGGAAGGGCACTGAGTCTGTCA
ENSG00000166888:ENST0000054387      CAAG---------------CGGTGTGAGCGGAAGGGCACTGAGTCTGTCA
ENSG00000166888:ENST0000055615      CAAG---------------CGGTGTGAGCGGAAGGGCACTGAGTCTGTCA
ENSG00000166888:ENST0000045407      CAAG---------------CGGTGTGAGCGGAAGGGCACTGAGTCTGTCA
ENSG00000166888:ENST0000053891      CAAG---------------CGGTGTGAGCGGAAGGGCACTGAGTCTGTCA
ENSG00000166888:ENST0000053721      CAAG---------------CGGTGTGAGCGGAAGGGCACTGAGTCTGTCA
ENSG00000166888:ENST0000053520      CAAG---------------CGGTGTGAGCGGAAGGGCACTGAGTCTGTCA
isotig46679                         GAAATCCAGTGCCGGAAGTAAAGGAAAT---GAGGGCTGCCACATGGTGA
isotig12565                         GAGATGTGGGAATGGAGGCCGTGCCAATTGTGATGCCTCCTTGATTGTGA
isotig12566                         GAGATGTGGGAATGGAGGCCGTGCCAATTGTGATGCCTCCTTGATTGTGA

ENSG00000166888:ENST0000030013      CAGAGGAGAAGTGCGCTGTGCTCTTCTCTGCCAGCTTCACACTTGGCCCC
ENSG00000166888:ENST0000054387      CAGAGGAGAAGTGCGCTGTGCTCTTCTCTGCCAGCTTCACACTTGGCCCC
ENSG00000166888:ENST0000055615      CAGAGGAGAAGTGCGCTGTGCTCTTCTCTGCCAGCTTCACACTTGGCCCC
ENSG00000166888:ENST0000045407      CAGAGGAGAAGTGCGCTGTGCTCTTCTCTGCCAGCTTCACACTTGGCCCC
ENSG00000166888:ENST0000053891      CAGAGGAGAAGTGCGCTGTGCTCTTCTCTGCCAGCTTCACACTTGGCCCC
ENSG00000166888:ENST0000053721      CAGAGGAGAAGTGCGCTGTGCTCTTCTCTGCCAGCTTCACACTTGGCCCC
ENSG00000166888:ENST0000053520      CAGAGGAGAAGTGCGCTGTGCTCTTCTCTGCCAGCTTCACACTTGGCCCC
isotig46679                         CGGAAGAGCTGCATTCCATAGCCTTTGAGACCCAGATCTGCCTC------
isotig12565                         CCGAGGAGCTGCATCTGATCACCTTCGAGACTGAGGTGTACCAC------
isotig12566                         CCGAGGAGCTGCATCTGATCACCTTCGAGACTGAGGTGTACCAC------

ENSG00000166888:ENST0000030013      GGCAAACTCCCCATCCAGCTCCAGGCCCTGTCTCTGCCCCTGGTGGTCAT
ENSG00000166888:ENST0000054387      GGCAAACTCCCCATCCAGCTCCAGGCCCTGTCTCTGCCCCTGGTGGTCAT
ENSG00000166888:ENST0000055615      GGCAAACTCCCCATCCAGCTCCAGGCCCTGTCTCTGCCCCTGGTGGTCAT
ENSG00000166888:ENST0000045407      GGCAAACTCCCCATCCAGCTCCAGGCCCTGTCTCTGCCCCTGGTGGTCAT
ENSG00000166888:ENST0000053891      GGCAAACTCCCCATCCAGCTCCAGGCCCTGTCTCTGCCCCTGGTGGTCAT
ENSG00000166888:ENST0000053721      GGCAAACTCCCCATCCAGCTCCAGGCCCTGTCTCTGCCCCTGGTGGTCAT
ENSG00000166888:ENST0000053520      GGCAAACTCCCCATCCAGCTCCAGGCCCTGTCTCTGCCCCTGGTGGTCAT
isotig46679                         TATGGCCTCACCATCGACTTGGAGACAAGCTCATTACCTGTGGTGATGAT
isotig12565                         CAAGGCCTCAAGATTGACCTGGAGACCCATTCTTTGCCAGTTGTGGTGAT
isotig12566                         CAAGGCCTCAAGATTGACCTGGAGACCCATTCTTTGCCAGTTGTGGTGAT

ENSG00000166888:ENST0000030013      CGTCCATGGCAACCAAGACAACAATGCCAAAGCCACTATCCTGTGGGACA
ENSG00000166888:ENST0000054387      CGTCCATGGCAACCAAGACAACAATGCCAAAGCCACTATCCTGTGGGACA
ENSG00000166888:ENST0000055615      CGTCCATGGCAACCAAGACAACAATGCCAAAGCCACTATCCTGTGGGACA
ENSG00000166888:ENST0000045407      CGTCCATGGCAACCAAGACAACAATGCCAAAGCCACTATCCTGTGGGACA
ENSG00000166888:ENST0000053891      CGTCCATGGCAACCAAGACAACAATGCCAAAGCCACTATCCTGTGGGACA
ENSG00000166888:ENST0000053721      CGTCCATGGCAACCAAGACAACAATGCCAAAGCCACTATCCTGTGGGACA
ENSG00000166888:ENST0000053520      CGTCCATGGCAACCAAGACAACAATGCCAAAGCCACTATCCTGTGGGACA
isotig46679                         TTCTAATGTCAGCCAACTGCCTAATGCTTGGGCATCCATCATTTGGTACA
isotig12565                         CTCCAACATCTGTCAGATGCCAAATGCCTGGGCATCCATCCTGTGGTATA
isotig12566                         CTCCAACATCTGTCAGATGCCAAATGCCTGGGCATCCATCCTGTGGTATA

ENSG00000166888:ENST0000030013      ATGCCTTCTCTGAG------ATGGACCGCGTGCCCTTTGTGGTGGCTGAG
ENSG00000166888:ENST0000054387      ATGCCTTCTCTGAG------ATGGACCGCGTGCCCTTTGTGGTGGCTGAG
ENSG00000166888:ENST0000055615      ATGCCTTCTCTGAG------ATGGACCGCGTGCCCTTTGTGGTGGCTGAG
ENSG00000166888:ENST0000045407      ATGCCTTCTCTGAG------ATGGACCGCGTGCCCTTTGTGGTGGCTGAG
ENSG00000166888:ENST0000053891      ATGCCTTCTCTGAG------ATGGACCGCGTGCCCTTTGTGGTGGCTGAG
ENSG00000166888:ENST0000053721      ATGCCTTCTCTGAG------ATGGACCGCGTGCCCTTTGTGGTGGCTGAG
ENSG00000166888:ENST0000053520      ATGCCTTCTCTGAG------ATGGACCGCGTGCCCTTTGTGGTGGCTGAG
isotig46679                         ATGTGTCAACCAACGATTGCCAGAACTTGGTTTTCTTTAATAATCCTCCG
isotig12565                         ACATGCTGACCAACAACCCCAAGAACGTGAACTTCTTCACCAAGCCACCA
isotig12566                         ACATGCTGACCAACAACCCCAAGAACGTGAACTTCTTCACCAAGCCACCA

ENSG00000166888:ENST0000030013      CGGGTGCCCTGGGAGAAGATGTGTGAAACTCTGAACCTGAAGTTCATGGC
ENSG00000166888:ENST0000054387      CGGGTGCCCTGGGAGAAGATGTGTGAAACTCTGAACCTGAAGTTCATGGC
ENSG00000166888:ENST0000055615      CGGGTGCCCTGGGAGAAGATGTGTGAAACTCTGAACCTGAAGTTCATGGC
ENSG00000166888:ENST0000045407      CGGGTGCCCTGGGAGAAGATGTGTGAAACTCTGAACCTGAAGTTCATGGC
ENSG00000166888:ENST0000053891      CGGGTGCCCTGGGAGAAGATGTGTGAAACTCTGAACCTGAAGTTCATGGC
ENSG00000166888:ENST0000053721      CGGGTGCCCTGGGAGAAGATGTGTGAAACTCTGAACCTGAAGTTCATGGC
ENSG00000166888:ENST0000053520      CGGGTGCCCTGGGAGAAGATGTGTGAAACTCTGAACCTGAAGTTCATGGC
isotig46679                         CCTGTCACTTTGAGTCAACTCCTGGAAGTGATGAGCTGGCAGTTTTCATC
isotig12565                         ATCGGAACCTGGGACCAGGTGGCCGAGGTGCTCAGCTGGCAGTTCTCATC
isotig12566                         ATCGGAACCTGGGACCAGGTGGCCGAGGTGCTCAGCTGGCAGTTCTCATC

ENSG00000166888:ENST0000030013      TGAGGTGGGGACCAACCGGGGGCTGCTCCCAGAGCACTTCCTCTTCCTGG
ENSG00000166888:ENST0000054387      TGAGGTGGGGACCAACCGGGGGCTGCTCCCAGAGCACTTCCTCTTCCTGG
ENSG00000166888:ENST0000055615      TGAGGTGGGGACCAACCGGGGGCTGCTCCCAGAGCACTTCCTCTTCCTGG
ENSG00000166888:ENST0000045407      TGAGGTGGGGACCAACCGGGGGCTGCTCCCAGAGCACTTCCTCTTCCTGG
ENSG00000166888:ENST0000053891      TGAGGTGGGGACCAACCGGGGGCTGCTCCCAGAGCACTTCCTCTTCCTGG
ENSG00000166888:ENST0000053721      TGAGGTGGGGACCAACCGGGGGCTGCTCCCAGAGCACTTCCTCTTCCTGG
ENSG00000166888:ENST0000053520      TGAGGTGGGGACCAACCGGGGGCTGCTCCCAGAGCACTTCCTCTTCCTGG
isotig46679                         CTATGTTGGT------CGTGGCCTTAATTCAGACCAGCTCAACATGCTGG
isotig12565                         CACCACAAAG------CGAGGGCTGAGCATCGAGCAGCTGACTACGCTGG
isotig12566                         CACCACAAAG------CGAGGGCTGAGCATCGAGCAGCTGACTACGCTGG

ENSG00000166888:ENST0000030013      CCCAGAAGATCTTCAATGACAACAGCCTCAGTATGGAGGCCTTCCAGCAC
ENSG00000166888:ENST0000054387      CCCAGAAGATCTTCAATGACAACAGCCTCAGTATGGAGGCCTTCCAGCAC
ENSG00000166888:ENST0000055615      CCCAGAAGATCTTCAATGACAACAGCCTCAGTATGGAGGCCTTCCAGCAC
ENSG00000166888:ENST0000045407      CCCAGAAGATCTTCAATGACAACAGCCTCAGTATGGAGGCCTTCCAGCAC
ENSG00000166888:ENST0000053891      CCCAGAAGATCTTCAATGACAACAGCCTCAGTATGGAGGCCTTCCAGCAC
ENSG00000166888:ENST0000053721      CCCAGAAGATCTTCAATGACAACAGCCTCAGTATGGAGGCCTTCCAGCAC
ENSG00000166888:ENST0000053520      CCCAGAAGATCTTCAATGACAACAGCCTCAGTATGGAGGCCTTCCAGCAC
isotig46679                         CAGAGAAGCTCACAGTTCAGTCT---------------AACTACAGCGAT
isotig12565                         CCGAGAAGCTCCTAGGACCTGGTGTC------------AACTACTCCGGG
isotig12566                         CCGAGAAGCTCCTAGGACCTGGTGTC------------AACTACTCCGGG

ENSG00000166888:ENST0000030013      CGTTCTGTGTCCTGGTCGCAGTTCAACAAGGAGATCCTGCTGGGCCGTGG
ENSG00000166888:ENST0000054387      CGTTCTGTGTCCTGGTCGCAGTTCAACAAGGAGATCCTGCTGGGCCGTGG
ENSG00000166888:ENST0000055615      CGTTCTGTGTCCTGGTCGCAGTTCAACAAGGAGATCCTGCTGGGCCGTGG
ENSG00000166888:ENST0000045407      CGTTCTGTGTCCTGGTCGCAGTTCAACAAGGAGATCCTGCTGGGCCGTGG
ENSG00000166888:ENST0000053891      CGTTCTGTGTCCTGGTCGCAGTTCAACAAGGAGATCCTGCTGGGCCGTGG
ENSG00000166888:ENST0000053721      CGTTCTGTGTCCTGGTCGCAGTTCAACAAGGAGATCCTGCTGGGCCGTGG
ENSG00000166888:ENST0000053520      CGTTCTGTGTCCTGGTCGCAGTTCAACAAGGAGATCCTGCTGGGCCGTGG
isotig46679                         GGTCACCTCACCTGGGCCAAGTTCTGCAAGGAACACTTGCCTGGCAAACC
isotig12565                         TGTCAGATCACATGGGCTAAATTTTGCAAAGAAAACATGGCTGGCAAGGG
isotig12566                         TGTCAGATCACATGGGCTAAATTTTGCAAAGAAAACATGGCTGGCAAGGG

ENSG00000166888:ENST0000030013      CTTCACCTTTTGGCAGTGGTTTGATGGTGTCCTGGACCTCACCAAACGCT
ENSG00000166888:ENST0000054387      CTTCACCTTTTGGCAGTGGTTTGATGGTGTCCTGGACCTCACCAAACGCT
ENSG00000166888:ENST0000055615      CTTCACCTTTTGGCAGTGGTTTGATGGTGTCCTGGACCTCACCAAACGCT
ENSG00000166888:ENST0000045407      CTTCACCTTTTGGCAGTGGTTTGATGGTGTCCTGGACCTCACCAAACGCT
ENSG00000166888:ENST0000053891      CTTCACCTTTTGGCAGTGGTTTGATGGTGTCCTGGACCTCACCAAACGCT
ENSG00000166888:ENST0000053721      CTTCACCTTTTGGCAGTGGTTTGATGGTGTCCTGGACCTCACCAAACGCT
ENSG00000166888:ENST0000053520      CTTCACCTTTTGGCAGTGGTTTGATGGTGTCCTGGACCTCACCAAACGCT
isotig46679                         ATTTACCTTCTGGACCTGGCTTGAAGCAATATTGGACCTAATTAAAAAAC
isotig12565                         CTTCTCCTTCTGGGTGTGGCTAGACAATATCATTGACCTTGTGAAAAAGT
isotig12566                         CTTCTCCTTCTGGGTGTGGCTAGACAATATCATTGACCTTGTGAAAAAGT

ENSG00000166888:ENST0000030013      GTCTCCGGAGCTACTGGTCTGACCGGCTGATCATTGGCTTCATCAGCAAA
ENSG00000166888:ENST0000054387      GTCTCCGGAGCTACTGGTCTGACCGGCTGATCATTGGCTTCATCAGCAAA
ENSG00000166888:ENST0000055615      GTCTCCGGAGCTACTGGTCTGACCGGCTGATCATTGGCTTCATCAGCAAA
ENSG00000166888:ENST0000045407      GTCTCCGGAGCTACTGGTCTGACCGGCTGATCATTGGCTTCATCAGCAAA
ENSG00000166888:ENST0000053891      GTCTCCGGAGCTACTGGTCTGACCGGCTGATCATTGGCTTCATCAGCAAA
ENSG00000166888:ENST0000053721      GTCTCCGGAGCTACTGGTCTGACCGGCTGATCATTGGCTTCATCAGCAAA
ENSG00000166888:ENST0000053520      GTCTCCGGAGCTACTGGTCTGACCGGCTGATCATTGGCTTCATCAGCAAA
isotig46679                         ACATTCTTCCCCTCTGGATTGATGGGTACATCATGGGCTTCGTGAGCAAA
isotig12565                         ATATCTTGGCCCTCTGGAATGAAGGGTACATCATGGGCTTCATTAGCAAG
isotig12566                         ATATCTTGGCCCTCTGGAATGAAGGGTACATCATGGGCTTCATTAGCAAG

ENSG00000166888:ENST0000030013      CAGTACGTTACTAGCCTTCTTCTCAATGAGCCCGACGGAACCTTTCTCCT
ENSG00000166888:ENST0000054387      CAGTACGTTACTAGCCTTCTTCTCAATGAGCCCGACGGAACCTTTCTCCT
ENSG00000166888:ENST0000055615      CAGTACGTTACTAGCCTTCTTCTCAATGAGCCCGACGGAACCTTTCTCCT
ENSG00000166888:ENST0000045407      CAGTACGTTACTAGCCTTCTTCTCAATGAGCCCGACGGAACCTTTCTCCT
ENSG00000166888:ENST0000053891      CAGTACGTTACTAGCCTTCTTCTCAATGAGCCCGACGGAACCTTTCTCCT
ENSG00000166888:ENST0000053721      CAGTACGTTACTAGCCTTCTTCTCAATGAGCCCGACGGAACCTTTCTCCT
ENSG00000166888:ENST0000053520      CAGTACGTTACTAGCCTTCTTCTCAATGAGCCCGACGGAACCTTTCTCCT
isotig46679                         GAGAAGGAGAGGTTTCTGCTCAAGGATAAAATGCCCGGGACATTTTTGTT
isotig12565                         GAGCGGGAGCGGGCGATCCTGAGCACGAAACCCCCGGGCACCTTCCTGCT
isotig12566                         GAGCGGGAGCGGGCGATCCTGAGCACGAAACCCCCGGGCACCTTCCTGCT

ENSG00000166888:ENST0000030013      CCGCTTCAGCGAC---TCAGAGATTGGGGGCATCACCATTGCCCATGTCA
ENSG00000166888:ENST0000054387      CCGCTTCAGCGAC---TCAGAGATTGGGGGCATCACCATTGCCCATGTCA
ENSG00000166888:ENST0000055615      CCGCTTCAGCGAC---TCAGAGATTGGGGGCATCACCATTGCCCATGTCA
ENSG00000166888:ENST0000045407      CCGCTTCAGCGAC---TCAGAGATTGGGGGCATCACCATTGCCCATGTCA
ENSG00000166888:ENST0000053891      CCGCTTCAGCGAC---TCAGAGATTGGGGGCATCACCATTGCCCATGTCA
ENSG00000166888:ENST0000053721      CCGCTTCAGCGAC---TCAGAGATTGGGGGCATCACCATTGCCCATGTCA
ENSG00000166888:ENST0000053520      CCGCTTCAGCGAC---TCAGAGATTGGGGGCATCACCATTGCCCATGTCA
isotig46679                         ACGATTCAGTGAG---AGCCATCTCGGAGGGATCACCTTCACCTGGGTGG
isotig12565                         GAGATTCAGCGAGAGCAGCAAAGAAGGAGGGGTCACTTTCACTTGGGTGG
isotig12566                         GAGATTCAGCGAGAGCAGCAAAGAAGGAGGGGTCACTTTCACTTGGGTGG

ENSG00000166888:ENST0000030013      TCCGGGGCCAGGATGGCTCTCCACAGATAGAGAACATCCAGCCATTCTCT
ENSG00000166888:ENST0000054387      TCCGGGGCCAGGATGGCTCTCCACAGATAGAGAACATCCAGCCATTCTCT
ENSG00000166888:ENST0000055615      TCCGGGGCCAGGATGGCTCTCCACAGATAGAGAACATCCAGCCATTCTCT
ENSG00000166888:ENST0000045407      TCCGGGGCCAGGATGGCTCTCCACAGATAGAGAACATCCAGCCATTCTCT
ENSG00000166888:ENST0000053891      TCCGGGGCCAGGATGGCTCTCCACAGATAGAGAACATCCAGCCATTCTCT
ENSG00000166888:ENST0000053721      TCCGGGGCCAGGATGGCTCTCCACAGATAGAGAACATCCAGCCATTCTCT
ENSG00000166888:ENST0000053520      TCCGGGGCCAGGATGGCTCTCCACAGATAGAGAACATCCAGCCATTCTCT
isotig46679                         ACCACTCTGAAAACGGAGAAGTGAGATTCCACTCCGTAGAACCCTACAAC
isotig12565                         AAAAGGACATCAGTGGCAAGACCCAGATCCAGTCTGTAGAGCCGTACACC
isotig12566                         AAAAGGACATCAGTGGCAAGACCCAGATCCAGTCTGTAGAGCCGTACACC

ENSG00000166888:ENST0000030013      GCCAAAGACCTGTCCATTCGCTCACTGGGGGACCGAATCCGGGAT-----
ENSG00000166888:ENST0000054387      GCCAAAGACCTGTCCATTCGCTCACTGGGGGACCGAATCCGGGAT-----
ENSG00000166888:ENST0000055615      GCCAAAGACCTGTCCATTCGCTCACTGGGGGACCGAATCCGGGAT-----
ENSG00000166888:ENST0000045407      GCCAAAGACCTGTCCATTCGCTCACTGGGGGACCGAATCCGGGAT-----
ENSG00000166888:ENST0000053891      GCCAAAGACCTGTCCATTCGCTCACTGGGGGACCGAATCCGGGAT-----
ENSG00000166888:ENST0000053721      GCCAAAGACCTGTCCATTCGCTCACTGGGGGACCGAATCCGGGAT-----
ENSG00000166888:ENST0000053520      GCCAAAGACCTGTCCATTCGCTCACTGGGGGACCGAATCCGGGAT-----
isotig46679                         AAAGGGCGTCTGTCGGCCCTGCCATTTGCTGACATCCTGCGGGACTACAA
isotig12565                         AAGCAGCAGCTGAACAACATGTCCTTTGCTGAAATCATCATGGGCTACAA
isotig12566                         AAGCAGCAGCTGAACAACATGTCCTTTGCTGAAATCATCATGGGCTACAA

ENSG00000166888:ENST0000030013      -------------------------CTTGCTCAGCTCAAAAATCTCTATC
ENSG00000166888:ENST0000054387      -------------------------CTTGCTCAGCTCAAAAATCTCTATC
ENSG00000166888:ENST0000055615      -------------------------CTTGCTCAGCTCAAAAATCTCTATC
ENSG00000166888:ENST0000045407      -------------------------CTTGCTCAGCTCAAAAATCTCTATC
ENSG00000166888:ENST0000053891      -------------------------CTTGCTCAGCTCAAAAATCTCTATC
ENSG00000166888:ENST0000053721      -------------------------CTTGCTCAGCTCAAAAATCTCTATC
ENSG00000166888:ENST0000053520      -------------------------CTTGCTCAGCTCAAAAATCTCTATC
isotig46679                         GGTCATCATGGCTGAGAACATTCCCGAGAACCCTCTCAAGTACCTCTACC
isotig12565                         GATCATGGATGCCACCAACATCCTGGTGTCCCCATTGGTCTACCTCTACC
isotig12566                         GATCATGGATGCCACCAACATCCTGGTGTCCCCATTGGTCTACCTCTACC

ENSG00000166888:ENST0000030013      CCAAGAAGCCCAAGGATGAGGCTTTCCGGAGCCACTAC------------
ENSG00000166888:ENST0000054387      CCAAGAAGCCCAAGGATGAGGCTTTCCGGAGCCACTAC------------
ENSG00000166888:ENST0000055615      CCAAGAAGCCCAAGGATGAGGCTTTCCGGAGCCACTAC------------
ENSG00000166888:ENST0000045407      CCAAGAAGCCCAAGGATGAGGCTTTCCGGAGCCACTAC------------
ENSG00000166888:ENST0000053891      CCAAGAAGCCCAAGGATGAGGCTTTCCGGAGCCACTAC------------
ENSG00000166888:ENST0000053721      CCAAGAAGCCCAAGGATGAGGCTTTCCGGAGCCACTAC------------
ENSG00000166888:ENST0000053520      CCAAGAAGCCCAAGGATGAGGCTTTCCGGAGCCACTAC------------
isotig46679                         CCGACATCCCCAAAGACAAAGCCTTCGGTAAACACTACAGCTCCCAGCCT
isotig12565                         CTGACATTCCCAAGGAGGAGGCGTTCGGGAAGTACTGT------------
isotig12566                         CTGACATTCCCAAGGAGGAGGCGTTCGGGAAGTACTGT------------

ENSG00000166888:ENST0000030013      ------------AAGCCTGAACAGATGGGTAAGGATGGCAGGGGTTATGT
ENSG00000166888:ENST0000054387      ------------AAGCCTGAACAGATGGGTAAGGATGGCAGGGGTTATGT
ENSG00000166888:ENST0000055615      ------------AAGCCTGAACAGATGGGTAAGGATGGCAGGGGTTATGT
ENSG00000166888:ENST0000045407      ------------AAGCCTGAACAGATGGGTAAGGATGGCAGGGGTTATGT
ENSG00000166888:ENST0000053891      ------------AAGCCTGAACAGATGGGTAAGGATGGCAGGGGTTATGT
ENSG00000166888:ENST0000053721      ------------AAGCCTGAACAGATGGGTAAGGATGGCAGGGGTTATGT
ENSG00000166888:ENST0000053520      ------------AAGCCTGAACAGATGGGTAAGGATGGCAGGGGTTATGT
isotig46679                         TGCGAAGTTTCAAGGCCAACA------GAACGGGGAGACAAAGGTTATGT
isotig12565                         ------------CGACCAGAG-----------------------------
isotig12566                         ------------CGACCAGAG-----------------------------

ENSG00000166888:ENST0000030013      CCCAGCTACCATCAAGATGACCGTGGAAAGGGACCAACCACTTCCTACCC
ENSG00000166888:ENST0000054387      CCCAGCTACCATCAAGATGACCGTGGAAAGGGACCAACCACTTCCTACCC
ENSG00000166888:ENST0000055615      CCCAGCTACCATCAAGATGACCGTGGAAAGGGACCAACCACTTCCTACCC
ENSG00000166888:ENST0000045407      CCCAGCTACCATCAAGATGACCGTGGAAAGGGACCAACCACTTCCTACCC
ENSG00000166888:ENST0000053891      CCCAGCTACCATCAAGATGACCGTGGAAAGGGACCAACCACTTCCTACCC
ENSG00000166888:ENST0000053721      CCCAGCTACCATCAAGATGACCGTGGAAAGGGACCAACCACTTCCTACCC
ENSG00000166888:ENST0000053520      CCCAGCTACCATCAAGATGACCGTGGAAAGGGACCAACCACTTCCTACCC
isotig46679                         TCCTTCAGTTTTTATCCCTATTTCAACAATCCGCAGCGACGCCATGGAGC
isotig12565                         -------------------------------AGCCAGGAGCATCCTGAAG
isotig12566                         -------------------------------AGCCAGGAGCATCCTGAAG

ENSG00000166888:ENST0000030013      CAGAGCTCCAGATGCCTACCATGGTGCCTTCTTATGACCTTGGAATGGCC
ENSG00000166888:ENST0000054387      CAGAGCTCCAGATGCCTACCATGGTGCCTTCTTATGACCTTGGAATGGCC
ENSG00000166888:ENST0000055615      CAGAGCTCCAGATGCCTACCATGGTGCCTTCTTATGACCTTGGAATGGCC
ENSG00000166888:ENST0000045407      CAGAGCTCCAGATGCCTACCATGGTGCCTTCTTATGACCTTGGAATGGCC
ENSG00000166888:ENST0000053891      CAGAGCTCCAGATGCCTACCATGGTGCCTTCTTATGACCTTGGAATGGCC
ENSG00000166888:ENST0000053721      CAGAGCTCCAGATGCCTACCATGGTGCCTTCTTATGACCTTGGAATGGCC
ENSG00000166888:ENST0000053520      CAGAGCTCCAGATGCCTACCATGGTGCCTTCTTATGACCTTGGAATGGCC
isotig46679                         CGCAG------------------------------------------TCT
isotig12565                         CTGAC---------------------------------------------
isotig12566                         CTGAC---------------------------------------------

ENSG00000166888:ENST0000030013      CCTGATTCCTCCATGAGCATGCAGCTTGGCCCAGATATGGTGCCCCAGGT
ENSG00000166888:ENST0000054387      CCTGATTCCTCCATGAGCATGCAGCTTGGCCCAGATATGGTGCCCCAGGT
ENSG00000166888:ENST0000055615      CCTGATTCCTCCATGAGCATGCAGCTTGGCCCAGATATGGTGCCCCAGGT
ENSG00000166888:ENST0000045407      CCTGATTCCTCCATGAGCATGCAGCTTGGCCCAGATATGGTGCCCCAGGT
ENSG00000166888:ENST0000053891      CCTGATTCCTCCATGAGCATGCAGCTTGGCCCAGATATGGTGCCCCAGGT
ENSG00000166888:ENST0000053721      CCTGATTCCTCCATGAGCATGCAGCTTGGCCCAGATATGGTGCCCCAGGT
ENSG00000166888:ENST0000053520      CCTGATTCCTCCATGAGCATGCAGCTTGGCCCAGATATGGTGCCCCAGGT
isotig46679                         CCTTCAGACCTTCTCCCC------------------ATGTCTCCGAGTGT
isotig12565                         CCCGGTAGTGCCGCCCCT------------------TACCTGAAGACCAA
isotig12566                         CCCGGTAGTTGTTTTTCCATG-----------------------------

ENSG00000166888:ENST0000030013      GTACCCACCACACTCTCACTCCATCCCCCCGTATCAAGGCCTCTCCCCAG
ENSG00000166888:ENST0000054387      GTACCCACCACACTCTCACTCCATCCCCCCGTATCAAGGCCTCTCCCCAG
ENSG00000166888:ENST0000055615      GTACCCACCACACTCTCACTCCATCCCCCCGTATCAAGGCCTCTCCCCAG
ENSG00000166888:ENST0000045407      GTACCCACCACACTCTCACTCCATCCCCCCGTATCAAGGCCTCTCCCCAG
ENSG00000166888:ENST0000053891      GTACCCACCACACTCTCACTCCATCCCCCCGTATCAAGGCCTCTCCCCAG
ENSG00000166888:ENST0000053721      GTACCCACCACACTCTCACTCCATCCCCCCGTATCAAGGCCTCTCCCCAG
ENSG00000166888:ENST0000053520      GTACCCACCACACTCTCACTCCATCCCCCCGTATCAAGGCCTCTCCCCAG
isotig46679                         ATACGCTGTGCTGAGAGAAAACCTGAGCCCT-------------------
isotig12565                         GTTCATCTGTGTGACACCAACGACCTGCAGC-------------------
isotig12566                         --------------------------------------------------

ENSG00000166888:ENST0000030013      AAGAATCAGTCAACGTGTTGTCAGCCTTCCAGGAGCCTCACCTGCAGATG
ENSG00000166888:ENST0000054387      AAGAATCAGTCAACGTGTTGTCAGCCTTCCAGGAGCCTCACCTGCAGATG
ENSG00000166888:ENST0000055615      AAGAATCAGTCAACGTGTTGTCAGCCTTCCAGGAGCCTCACCTGCAGATG
ENSG00000166888:ENST0000045407      AAGAATCAGTCAACGTGTTGTCAGCCTTCCAGGAGCCTCACCTGCAGATG
ENSG00000166888:ENST0000053891      AAGAATCAGTCAACGTGTTGTCAGCCTTCCAGGAGCCTCACCTGCAGATG
ENSG00000166888:ENST0000053721      AAGAATCAGTCAACGTGTTGTCAGCCTTCCAGGAGCCTCACCTGCAGATG
ENSG00000166888:ENST0000053520      AAGAATCAGTCAACGTGTTGTCAGCCTTCCAGGAGCCTCACCTGCAGATG
isotig46679                         --------------------------------------------------
isotig12565                         --------------------------------------------------
isotig12566                         --------------------------------------------------

ENSG00000166888:ENST0000030013      CCCCCCAGCCTGGGCCAGATGAGCCTGCCCTTTGACCAGCCTCACCCCCA
ENSG00000166888:ENST0000054387      CCCCCCAGCCTGGGCCAGATGAGCCTGCCCTTTGACCAGCCTCACCCCCA
ENSG00000166888:ENST0000055615      CCCCCCAGCCTGGGCCAGATGAGCCTGCCCTTTGACCAGCCTCACCCCCA
ENSG00000166888:ENST0000045407      CCCCCCAGCCTGGGCCAGATGAGCCTGCCCTTTGACCAGCCTCACCCCCA
ENSG00000166888:ENST0000053891      CCCCCCAGCCTGGGCCAGATGAGCCTGCCCTTTGACCAGCCTCACCCCCA
ENSG00000166888:ENST0000053721      CCCCCCAGCCTGGGCCAGATGAGCCTGCCCTTTGACCAGCCTCACCCCCA
ENSG00000166888:ENST0000053520      CCCCCCAGCCTGGGCCAGATGAGCCTGCCCTTTGACCAGCCTCACCCCCA
isotig46679                         ------------ACCACAATTGAAACAGCAATGAAGTCTCCATATTCTGA
isotig12565                         ------------AATACCATTGACCTGCCGATGTCCCCCCCGCAC-----
isotig12566                         ------------GTTCTGGTTTCGCTGTTAGGGAAAGGGGGACAGTGCAG

ENSG00000166888:ENST0000030013      GGGCCTGCTGCCGTGCCAGCCTCAGGAGCATGCTGTGTCCAGCCCTGACC
ENSG00000166888:ENST0000054387      GGGCCTGCTGCCGTGCCAGCCTCAGGAGCATGCTGTGTCCAGCCCTGACC
ENSG00000166888:ENST0000055615      GGGCCTGCTGCCGTGCCAGCCTCAGGAGCATGCTGTGTCCAGCCCTGACC
ENSG00000166888:ENST0000045407      GGGCCTGCTGCCGTGCCAGCCTCAGGAGCATGCTGTGTCCAGCCCTGACC
ENSG00000166888:ENST0000053891      GGGCCTGCTGCCGTGCCAGCCTCAGGAGCATGCTGTGTCCAGCCCTGACC
ENSG00000166888:ENST0000053721      GGGCCTGCTGCCGTGCCAGCCTCAGGAGCATGCTGTGTCCAGCCCTGACC
ENSG00000166888:ENST0000053520      GGGCCTGCTGCCGTGCCAGCCTCAGGAGCATGCTGTGTCCAGCCCTGACC
isotig46679                         GCGGTAC---------------------AAAGCGACTCTTCAAGGAAGAG
isotig12565                         --------------------------------------------------
isotig12566                         GTCCTTG---------------------------------------GAGG

ENSG00000166888:ENST0000030013      CCCTGCTCTGCTCAGATGTGACCATGGTGGAAGACAGCTGCCTGAGCCAG
ENSG00000166888:ENST0000054387      CCCTGCTCTGCTCAGATGTGACCATGGTGGAAGACAGCTGCCTGAGCCAG
ENSG00000166888:ENST0000055615      CCCTGCTCTGCTCAGATGTGACCATGGTGGAAGACAGCTGCCTGAGCCAG
ENSG00000166888:ENST0000045407      CCCTGCTCTGCTCAGATGTGACCATGGTGGAAGACAGCTGCCTGAGCCAG
ENSG00000166888:ENST0000053891      CCCTGCTCTGCTCAGATGTGACCATGGTGGAAGACAGCTGCCTGAGCCAG
ENSG00000166888:ENST0000053721      CCCTGCTCTGCTCAGATGTGACCATGGTGGAAGACAGCTGCCTGAGCCAG
ENSG00000166888:ENST0000053520      CCCTGCTCTGCTCAGATGTGACCATGGTGGAAGACAGCTGCCTGAGCCAG
isotig46679                         AGCAGATGAAAACGGAGACTGCTCTTTGCCAAAGTCCACAATTCATTTCT
isotig12565                         --------------------------------------------------
isotig12566                         AGAGACAAGGACATGACCGGGTGTCTGGTGGTGAGTCCTGCTATGGAAGA

ENSG00000166888:ENST0000030013      CCAGTGACAGCGTTTCCTCAGGGCACTTGGATTGGTGAAGACATATTCCC
ENSG00000166888:ENST0000054387      CCAGTGACAGCGTTTCCTCAGGGCACTTGGATTGGTGAAGACATATTCCC
ENSG00000166888:ENST0000055615      CCAGTGACAGCGTTTCCTCAGGGCACTTGGATTGGTGAAGACATATTCCC
ENSG00000166888:ENST0000045407      CCAGTGACAGCGTTTCCTCAGGGCACTTGGATTGGTGAAGACATATTCCC
ENSG00000166888:ENST0000053891      CCAGTGACAGCGTTTCCTCAGGGCACTTGGATTGGTGAAGACATATTCCC
ENSG00000166888:ENST0000053721      CCAGTGACAGCGTTTCCTCAGGGCACTTGGATTGGTGAAGACATATTCCC
ENSG00000166888:ENST0000053520      CCAGTGACAGCGTTTCCTCAGGGCACTTGGATTGGTGAAGACATATTCCC
isotig46679                         TCAGCTTTGATACTGGTTTCTAGAAAATGGCACAAATCCGAAGCTTTCCT
isotig12565                         ------TTTAGATTCATTGATGCAGTTTGGAAACGGAGGTGCGCCCTC--
isotig12566                         GCTGTTTAT------------------TGGGTACTTCAG-----------

ENSG00000166888:ENST0000030013      TCCTCTGCTGCCTCCCACTGAACAGGACCTCACTAAGCTTCTCCTGGAGG
ENSG00000166888:ENST0000054387      TCCTCTGCTGCCTCCCACTGAACAGGACCTCACTAAGCTTCTCCTGGAGG
ENSG00000166888:ENST0000055615      TCCTCTGCTGCCTCCCACTGAACAGGACCTCACTAAGCTTCTCCTGGAGG
ENSG00000166888:ENST0000045407      TCCTCTGCTGCCTCCCACTGAACAGGACCTCACTAAGCTTCTCCTGGAGG
ENSG00000166888:ENST0000053891      TCCTCTGCTGCCTCCCACTGAACAGGACCTCACTAAGCTTCTCCTGGAGG
ENSG00000166888:ENST0000053721      TCCTCTGCTGCCTCCCACTGAACAGGACCTCACTAAGCTTCTCCTGGAGG
ENSG00000166888:ENST0000053520      TCCTCTGCTGCCTCCCACTGAACAGGACCTCACTAAGCTTCTCCTGGAGG
isotig46679                         CTCACTA------------------------------------------G
isotig12565                         -------------------------------------------------G
isotig12566                         -------------------------------------------------G

ENSG00000166888:ENST0000030013      GGCAAGGGGAGTCGGGGGGAGGGTCCTTGGGGGCACAGCCCCTCCTGCAG
ENSG00000166888:ENST0000054387      GGCAAGGGGAGTCGGGGGGAGGGTCCTTGGGGGCACAGCCCCTCCTGCAG
ENSG00000166888:ENST0000055615      GGCAAGGGGAGTCGGGGGGAGGGTCCTTGGGGGCACAGCCCCTCCTGCAG
ENSG00000166888:ENST0000045407      GGCAAGGGGAGTCGGGGGGAGGGTCCTTGGGGGCACAGCCCCTCCTGCAG
ENSG00000166888:ENST0000053891      GGCAAGGGGAGTCGGGGGGAGGGTCCTTGGGGGCACAGCCCCTCCTGCAG
ENSG00000166888:ENST0000053721      GGCAAGGGGAGTCGGGGGGAGGGTCCTTGGGGGCACAGCCCCTCCTGCAG
ENSG00000166888:ENST0000053520      GGCAAGGGGAGTCGGGGGGAGGGTCCTTGGGGGCACAGCCCCTCCTGCAG
isotig46679                         GTGACATTCCCCAACTGGGAGTGCTGCTGAAATGCAAACCAAAGCTTCAG
isotig12565                         GCAGGAGGGCAGTTGTCACTCACGTTCATGGATCTGAC------------
isotig12566                         GTGACCGGGATTCAAGAGAAGACCAGAATCAGGCCTCA------------

ENSG00000166888:ENST0000030013      CCCTCCCACTATGGGCAATCTGGGATCTCAATGTCCCACATGGACCTAAG
ENSG00000166888:ENST0000054387      CCCTCCCACTATGGGCAATCTGGGATCTCAATGTCCCACATGGACCTAAG
ENSG00000166888:ENST0000055615      CCCTCCCACTATGGGCAATCTGGGATCTCAATGTCCCACATGGACCTAAG
ENSG00000166888:ENST0000045407      CCCTCCCACTATGGGCAATCTGGGATCTCAATGTCCCACATGGACCTAAG
ENSG00000166888:ENST0000053891      CCCTCCCACTATGGGCAATCTGGGATCTCAATGTCCCACATGGACCTAAG
ENSG00000166888:ENST0000053721      CCCTCCCACTATGGGCAATCTGGGATCTCAATGTCCCACATGGACCTAAG
ENSG00000166888:ENST0000053520      CCCTCCCACTATGGGCAATCTGGGATCTCAATGTCCCACATGGACCTAAG
isotig46679                         ATAAACACGCAGGAAAAGACAGCTTCGAGAAACCTATGTTCGCAATATAA
isotig12565                         ---------TTCGGAGTGCGCTACCTCCCCCATGTGAGGAGC--------
isotig12566                         --------------------------------------------------

ENSG00000166888:ENST0000030013      GGCCAACCCCAGTTGG
ENSG00000166888:ENST0000054387      GGCCAACCCCAGTTGG
ENSG00000166888:ENST0000055615      GGCCAACCCCAGTTGG
ENSG00000166888:ENST0000045407      GGCCAACCCCAGTTGG
ENSG00000166888:ENST0000053891      GGCCAACCCCAGTTGG
ENSG00000166888:ENST0000053721      GGCCAACCCCAGTTGG
ENSG00000166888:ENST0000053520      GGCCAACCCCAGTTGG
isotig46679                         CAGAAGGCTGCTTTGC
isotig12565                         ----------------
isotig12566                         ----------------


""",
        )

    def test5(self):
        aligner = CodonAligner()
        # aligner.frameshift_score = -10.0
        nucleotide_records = SeqIO.parse("codonalign/nucl5.fa", "fasta")
        protein_alignment = Align.read("codonalign/pro5.aln", "clustal")
        self.assertEqual(len(protein_alignment.sequences), 3)
        codon_alignments = []
        nucleotide_record = next(nucleotide_records)
        protein_record = protein_alignment.sequences[0]
        self.assertEqual(nucleotide_record.id, protein_record.id)
        alignments = aligner.align(protein_record, nucleotide_record)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 183], [0, 549]]))
        )
        self.assertEqual(
            str(alignment),
            """\
isotig697         0 R  G  D  Q  R  S  N  F  Q  L  S  P  S  T  M  Q  I  S  T  G  
isotig697         0 TGAGGCGATCAACGCAGCAACTTCCAGCTGTCTCCCTCCACCATGCAGATCTCCACAGGG

isotig697        20 L  L  C  L  L  L  V  A  T  G  F  T  S  Q  V  L  A  H  P  G  
isotig697        60 CTTCTGTGCcTGCTGCTTGTGGCCACTGGCTTCACTTCCCAGGTGCTGGCTCACCCAGGC

isotig697        40 S  I  P  S  T  Y  C  F  V  M  T  S  K  K  I  P  K  S  L  L  
isotig697       120 TCTATCCCATCTACCTaCTGCTTTGTTATGACCAGTAAGAaGATCCCCAAATCACTACTG

isotig697        60 K  S  Y  K  R  I  S  N  S  R  C  T  L  K  A  I  L  F  K  T  
isotig697       180 AaGAGCTACAAAaGAATCTCCAACAGCAGATGCACCcTGAAAGCCATACTCTTCAAGACC

isotig697        80 K  S  G  K  E  I  C  A  D  P  K  K  K  W  V  Q  D  A  T  K  
isotig697       240 AAGTCGGGCAAAGAGATCTGTGCTGACCCCAAGAAGAAGTGGGTCcAGGATGCCACAAAG

isotig697       100 H  L  D  Q  I  L  Q  T  P  K  P  T  I  P  S  F  E  T  H  P  
isotig697       300 CACCTGGACCAAATCCTTCAAACTCCAAAACCGACAATCCCCTCTTTTGAGACTCACCCA

isotig697       120 E  T  K  K  C  F  I  H  S  P  F  L  R  R  A  P  R  S  T  Q  
isotig697       360 GAGACTAAGAAATGCTTCATTCATTCTCCATTCCTAAGACGTGCTCCAAGGTCAACTCAG

isotig697       140 H  H  S  P  R  T  W  L  H  L  V  M  D  R  T  E  S  H  Y  V  
isotig697       420 CACCATTCCCCAAGGACTTGGCTTCATTTAGTTATGGATAGAACTGAAAGTCATTATGTT

isotig697       160 Q  N  K  P  D  L  K  R  L  C  N  F  L  N  M  Q  N  L  K  R  
isotig697       480 CAGAATAAGCCAGACTTGAAGAGGTTGTGTAATTTCTTGAATATGCAAAATCTTAAAAGG

isotig697       180 G  A  C   183
isotig697       540 GGGGCATGC 549
""",
        )
        nucleotide_record = next(nucleotide_records)
        protein_record = protein_alignment.sequences[1]
        self.assertEqual(nucleotide_record.id, protein_record.id)
        alignments = aligner.align(protein_record, nucleotide_record)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 65], [0, 195]]))
        )
        self.assertEqual(
            str(alignment),
            """\
ENSG00000         0 M  K  V  S  A  A  L  L  C  L  L  L  I  A  A  T  F  I  P  Q  
ENSG00000         0 ATGAAAGTCTCTGCCGCCCTTCTGTGCCTGCTGCTCATAGCAGCCACCTTCATTCCCCAA

ENSG00000        20 G  L  A  Q  P  D  A  I  N  A  P  V  T  C  C  Y  N  F  T  N  
ENSG00000        60 GGGCTCGCTCAGCCAGATGCAATCAATGCCCCAGTCACCTGCTGTTATAACTTCACCAAT

ENSG00000        40 R  K  I  S  V  Q  R  L  A  S  Y  R  R  I  T  S  S  K  C  P  
ENSG00000       120 AGGAAGATCTCAGTGCAGAGGCTCGCGAGCTATAGAAGAATCACCAGCAGCAAGTGTCCC

ENSG00000        60 K  E  A  V  M    65
ENSG00000       180 AAAGAAGCTGTGATG 195
""",
        )
        nucleotide_record = next(nucleotide_records)
        protein_record = protein_alignment.sequences[2]
        self.assertEqual(nucleotide_record.id, protein_record.id)
        alignments = aligner.align(protein_record, nucleotide_record)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 99], [9, 306]]))
        )
        self.assertEqual(
            str(alignment),
            """\
ENSG00000         0 M  K  V  S  A  A  L  L  C  L  L  L  I  A  A  T  F  I  P  Q  
ENSG00000         9 ATGAAAGTCTCTGCCGCCCTTCTGTGCCTGCTGCTCATAGCAGCCACCTTCATTCCCCAA

ENSG00000        20 G  L  A  Q  P  D  A  I  N  A  P  V  T  C  C  Y  N  F  T  N  
ENSG00000        69 GGGCTCGCTCAGCCAGATGCAATCAATGCCCCAGTCACCTGCTGTTATAACTTCACCAAT

ENSG00000        40 R  K  I  S  V  Q  R  L  A  S  Y  R  R  I  T  S  S  K  C  P  
ENSG00000       129 AGGAAGATCTCAGTGCAGAGGCTCGCGAGCTATAGAAGAATCACCAGCAGCAAGTGTCCC

ENSG00000        60 K  E  A  V  I  F  K  T  I  V  A  K  E  I  C  A  D  P  K  Q  
ENSG00000       189 AAAGAAGCTGTGATCTTCAAGACCATTGTGGCCAAGGAGATCTGTGCTGACCCCAAGCAG

ENSG00000        80 K  W  V  Q  D  S  M  D  H  L  D  K  Q  T  Q  T  P  K  T  
ENSG00000       249 AAGTGGGTTCAGGATTCCATGGACCACCTGGACAAGCAAACCCAAACTCCGAAGACT

ENSG00000        99
ENSG00000       306
""",
        )
        alignment = protein_alignment.mapall(codon_alignments)
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[0, 42, 126, 126, 231, 333, 549],
                          [0,  0,  84,  90, 195, 195, 195],
                          [9,  9,  93,  99, 204, 306, 306]])
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "clustal"),
            """\
isotig69710                         TGAGGCGATCAACGCAGCAACTTCCAGCTGTCTCCCTCCACCATGCAGAT
ENSG00000108691:ENST0000058090      ------------------------------------------ATGAAAGT
ENSG00000108691:ENST0000022583      ------------------------------------------ATGAAAGT

isotig69710                         CTCCACAGGGCTTCTGTGCcTGCTGCTTGTGGCCACTGGCTTCACTTCCC
ENSG00000108691:ENST0000058090      CTCTGCCGCCCTTCTGTGCCTGCTGCTCATAGCAGCCACCTTCATTCCCC
ENSG00000108691:ENST0000022583      CTCTGCCGCCCTTCTGTGCCTGCTGCTCATAGCAGCCACCTTCATTCCCC

isotig69710                         AGGTGCTGGCTCACCCAGGCTCTATC------CCATCTACCTaCTGCTTT
ENSG00000108691:ENST0000058090      AAGGGCTCGCTCAGCCAGATGCAATCAATGCCCCAGTCACCTGCTGTTAT
ENSG00000108691:ENST0000022583      AAGGGCTCGCTCAGCCAGATGCAATCAATGCCCCAGTCACCTGCTGTTAT

isotig69710                         GTTATGACCAGTAAGAaGATCCCCAAATCACTACTGAaGAGCTACAAAaG
ENSG00000108691:ENST0000058090      AACTTCACCAATAGGAAGATCTCAGTGCAGAGGCTCGCGAGCTATAGAAG
ENSG00000108691:ENST0000022583      AACTTCACCAATAGGAAGATCTCAGTGCAGAGGCTCGCGAGCTATAGAAG

isotig69710                         AATCTCCAACAGCAGATGCACCcTGAAAGCCATACTCTTCAAGACCAAGT
ENSG00000108691:ENST0000058090      AATCACCAGCAGCAAGTGTCCCAAAGAAGCTGTGATG-------------
ENSG00000108691:ENST0000022583      AATCACCAGCAGCAAGTGTCCCAAAGAAGCTGTGATCTTCAAGACCATTG

isotig69710                         CGGGCAAAGAGATCTGTGCTGACCCCAAGAAGAAGTGGGTCcAGGATGCC
ENSG00000108691:ENST0000058090      --------------------------------------------------
ENSG00000108691:ENST0000022583      TGGCCAAGGAGATCTGTGCTGACCCCAAGCAGAAGTGGGTTCAGGATTCC

isotig69710                         ACAAAGCACCTGGACCAAATCCTTCAAACTCCAAAACCGACAATCCCCTC
ENSG00000108691:ENST0000058090      --------------------------------------------------
ENSG00000108691:ENST0000022583      ATGGACCACCTGGACAAGCAAACCCAAACTCCGAAGACT-----------

isotig69710                         TTTTGAGACTCACCCAGAGACTAAGAAATGCTTCATTCATTCTCCATTCC
ENSG00000108691:ENST0000058090      --------------------------------------------------
ENSG00000108691:ENST0000022583      --------------------------------------------------

isotig69710                         TAAGACGTGCTCCAAGGTCAACTCAGCACCATTCCCCAAGGACTTGGCTT
ENSG00000108691:ENST0000058090      --------------------------------------------------
ENSG00000108691:ENST0000022583      --------------------------------------------------

isotig69710                         CATTTAGTTATGGATAGAACTGAAAGTCATTATGTTCAGAATAAGCCAGA
ENSG00000108691:ENST0000058090      --------------------------------------------------
ENSG00000108691:ENST0000022583      --------------------------------------------------

isotig69710                         CTTGAAGAGGTTGTGTAATTTCTTGAATATGCAAAATCTTAAAAGGGGGG
ENSG00000108691:ENST0000058090      --------------------------------------------------
ENSG00000108691:ENST0000022583      --------------------------------------------------

isotig69710                         CATGC
ENSG00000108691:ENST0000058090      -----
ENSG00000108691:ENST0000022583      -----


""",
        )


class Test_build(unittest.TestCase):
    def test_build1(self):
        aligner = CodonAligner()
        codon_alignments = []
        seq1 = SeqRecord(
            Seq(
                "TCAGGGACTGCGAGAACCAAGCTACTGCTGCTGCTGGCTGCGCTCTGCGCCGCAGGTGGGGCGCTGGAG"
            ),
            id="pro1",
        )
        seq2 = SeqRecord(
            Seq("TCAGGGACTTCGAGAACCAAGCGCTCCTGCTGCTGGCTGCGCTCGGCGCCGCAGGTGGAGCACTGGAG"),
            id="pro2",
        )
        pro1 = SeqRecord(Seq("SGTARTKLLLLLAALCAAGGALE"), id="pro1")
        pro2 = SeqRecord(Seq("SGTSRTKRLLLLAALGAAGGALE"), id="pro2")
        alignments = aligner.align(pro1, seq1)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 23], [0, 69]]))
        )
        self.assertEqual(
            str(alignment),
            """\
pro1              0 S  G  T  A  R  T  K  L  L  L  L  L  A  A  L  C  A  A  G  G  
pro1              0 TCAGGGACTGCGAGAACCAAGCTACTGCTGCTGCTGGCTGCGCTCTGCGCCGCAGGTGGG

pro1             20 A  L  E   23
pro1             60 GCGCTGGAG 69
""",
        )
        alignments = aligner.align(pro2, seq2)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(
                alignment.coordinates, np.array([[0, 8, 8, 23], [0, 24, 23, 68]])
            )
        )
        self.assertEqual(
            str(alignment),
            """\
pro2              0 S  G  T  S  R  T  K  R   8
pro2              0 TCAGGGACTTCGAGAACCAAGCGC 24

pro2              8 L  L  L  L  A  A  L  G  A  A  G  G  A  L  E   23
pro2             23 CTCCTGCTGCTGGCTGCGCTCGGCGCCGCAGGTGGAGCACTGGAG 68
""",
        )
        alignment = Alignment([pro1, pro2])
        alignment = alignment.mapall(codon_alignments)
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[0, 24, 24, 69],
                          [0, 24, 23, 68]])
                # fmt: on
            )
        )
        self.assertEqual(
            str(alignment),
            """\
pro1              0 TCAGGGACTGCGAGAACCAAGCTA 24
                  0 |||||||||.||||||||||||..
pro2              0 TCAGGGACTTCGAGAACCAAGCGC 24

pro1             24 CTGCTGCTGCTGGCTGCGCTCTGCGCCGCAGGTGGGGCGCTGGAG 69
                 24 ||.||||||||||||||||||.|||||||||||||.||.|||||| 69
pro2             23 CTCCTGCTGCTGGCTGCGCTCGGCGCCGCAGGTGGAGCACTGGAG 68
""",
        )

    def test_build2(self):
        aligner = CodonAligner()
        codon_alignments = []
        seq1 = SeqRecord(
            Seq(
                "ATGAAAAAGCACGAGTTACTTTGCCAAGGGACAAGTAACAAGCTCACCCAGTTGGGCACTTTTGAAGACCACTTTCTGAGCCTACAGAGGATGTTCAACAACTGTGAGGTGGTCCTTGGGAATTTGGAAATTACCTACATGCAGAGTAGTTACAACCTTTCTTTTCTCAAGACCATCCAGGAGGTTGCCGGCTATGTACTCATTGCCCTC"
            ),
            id="pro1",
        )
        seq2 = SeqRecord(
            Seq(
                "ATGAAAAAGCACGAGTTCTTTGCCAAGGGACAAGTAACAAGCTCACCCAGTTGGGCACTTTTGAAGACCACTTTCTGAGCCTACAGAGGATGTTCAACAATGTGAGGTGGTCCTTGGGAATTTGGAAATTACCTACATGCAGAGTAGTTACAACCTTTCTTTTCTCAAGACCATCCAGGAGGTTGCCGGCTATGTACTCATTGCCCTC"
            ),
            id="pro2",
        )
        seq3 = SeqRecord(
            Seq(
                "ATGAAAAAGCACGAGTTACTTTGCCAAGGGACAAGTAACAAGCTCACCCTTGGGCACTTTTGAAGACCACTTTCTGAGCCTACAGAGGATGTTCAACAACTGTGAGGTGGTCCTTGGGAATTTGGAAATTACCTACATGCAGAGTAGTTACAACCTTTCTTTTCTCAAGACCATCCAGGAGGTTGCCGGCTATGTACTCATTGCCCTC"
            ),
            id="pro3",
        )
        pro1 = SeqRecord(
            Seq(
                "MKKHELLCQGTSNKLTQLGTFEDHFLSLQRMFNNCEVVLGNLEITYMQSSYNLSFLKTIQEVAGYVLIAL"
            ),
            id="pro1",
        )
        pro2 = SeqRecord(
            Seq(
                "MKKHEFLCQGTSNKLTQLGTFEDHFLSLQRMFNNCEVVLGNLEITYMQSSYNLSFLKTIQEVAGYVLIAL"
            ),
            id="pro2",
        )
        pro3 = SeqRecord(
            Seq(
                "MKKHELLCQGTSNKLTLLGTFEDHFLSLQRMFNNCEVVLGNLEITYMQSSYNLSFLKTIQEVAGYVLIAL"
            ),
            id="pro3",
        )
        alignments = aligner.align(pro1, seq1)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 70], [0, 210]]))
        )
        self.assertEqual(
            str(alignment),
            """\
pro1              0 M  K  K  H  E  L  L  C  Q  G  T  S  N  K  L  T  Q  L  G  T  
pro1              0 ATGAAAAAGCACGAGTTACTTTGCCAAGGGACAAGTAACAAGCTCACCCAGTTGGGCACT

pro1             20 F  E  D  H  F  L  S  L  Q  R  M  F  N  N  C  E  V  V  L  G  
pro1             60 TTTGAAGACCACTTTCTGAGCCTACAGAGGATGTTCAACAACTGTGAGGTGGTCCTTGGG

pro1             40 N  L  E  I  T  Y  M  Q  S  S  Y  N  L  S  F  L  K  T  I  Q  
pro1            120 AATTTGGAAATTACCTACATGCAGAGTAGTTACAACCTTTCTTTTCTCAAGACCATCCAG

pro1             60 E  V  A  G  Y  V  L  I  A  L    70
pro1            180 GAGGTTGCCGGCTATGTACTCATTGCCCTC 210
""",
        )
        alignments = aligner.align(pro2, seq2)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                np.array([[0, 6, 6, 34, 34, 70], [0, 18, 17, 101, 100, 208]]),
            )
        )
        self.assertEqual(
            str(alignment),
            """\
pro2              0 M  K  K  H  E  F   6
pro2              0 ATGAAAAAGCACGAGTTC 18

pro2              6 L  C  Q  G  T  S  N  K  L  T  Q  L  G  T  F  E  D  H  F  L  
pro2             17 CTTTGCCAAGGGACAAGTAACAAGCTCACCCAGTTGGGCACTTTTGAAGACCACTTTCTG

pro2             26 S  L  Q  R  M  F  N  N   34
pro2             77 AGCCTACAGAGGATGTTCAACAAT 101

pro2             34 C  E  V  V  L  G  N  L  E  I  T  Y  M  Q  S  S  Y  N  L  S  
pro2            100 TGTGAGGTGGTCCTTGGGAATTTGGAAATTACCTACATGCAGAGTAGTTACAACCTTTCT

pro2             54 F  L  K  T  I  Q  E  V  A  G  Y  V  L  I  A  L    70
pro2            160 TTTCTCAAGACCATCCAGGAGGTTGCCGGCTATGTACTCATTGCCCTC 208
""",
        )
        alignments = aligner.align(pro3, seq3)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(
                alignment.coordinates, np.array([[0, 17, 17, 70], [0, 51, 49, 208]])
            )
        )
        self.assertEqual(
            str(alignment),
            """\
pro3              0 M  K  K  H  E  L  L  C  Q  G  T  S  N  K  L  T  L   17
pro3              0 ATGAAAAAGCACGAGTTACTTTGCCAAGGGACAAGTAACAAGCTCACCCTT 51

pro3             17 L  G  T  F  E  D  H  F  L  S  L  Q  R  M  F  N  N  C  E  V  
pro3             49 TTGGGCACTTTTGAAGACCACTTTCTGAGCCTACAGAGGATGTTCAACAACTGTGAGGTG

pro3             37 V  L  G  N  L  E  I  T  Y  M  Q  S  S  Y  N  L  S  F  L  K  
pro3            109 GTCCTTGGGAATTTGGAAATTACCTACATGCAGAGTAGTTACAACCTTTCTTTTCTCAAG

pro3             57 T  I  Q  E  V  A  G  Y  V  L  I  A  L    70
pro3            169 ACCATCCAGGAGGTTGCCGGCTATGTACTCATTGCCCTC 208
""",
        )
        alignment = Alignment([pro1, pro2, pro3])
        alignment = alignment.mapall(codon_alignments)
        self.assertEqual(
            str(alignment),
            """\
pro1              0 ATGAAAAAGCACGAGTTA 18
pro2              0 ATGAAAAAGCACGAGTTC 18
pro3              0 ATGAAAAAGCACGAGTTA 18

pro1             18 CTTTGCCAAGGGACAAGTAACAAGCTCACCCAG 51
pro2             17 CTTTGCCAAGGGACAAGTAACAAGCTCACCCAG 50
pro3             18 CTTTGCCAAGGGACAAGTAACAAGCTCACCCTT 51

pro1             51 TTGGGCACTTTTGAAGACCACTTTCTGAGCCTACAGAGGATGTTCAACAAC 102
pro2             50 TTGGGCACTTTTGAAGACCACTTTCTGAGCCTACAGAGGATGTTCAACAAT 101
pro3             49 TTGGGCACTTTTGAAGACCACTTTCTGAGCCTACAGAGGATGTTCAACAAC 100

pro1            102 TGTGAGGTGGTCCTTGGGAATTTGGAAATTACCTACATGCAGAGTAGTTACAACCTTTCT
pro2            100 TGTGAGGTGGTCCTTGGGAATTTGGAAATTACCTACATGCAGAGTAGTTACAACCTTTCT
pro3            100 TGTGAGGTGGTCCTTGGGAATTTGGAAATTACCTACATGCAGAGTAGTTACAACCTTTCT

pro1            162 TTTCTCAAGACCATCCAGGAGGTTGCCGGCTATGTACTCATTGCCCTC 210
pro2            160 TTTCTCAAGACCATCCAGGAGGTTGCCGGCTATGTACTCATTGCCCTC 208
pro3            160 TTTCTCAAGACCATCCAGGAGGTTGCCGGCTATGTACTCATTGCCCTC 208
""",
        )
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[0, 18, 18, 51, 51, 102, 102, 210],
                          [0, 18, 17, 50, 50, 101, 100, 208],
                          [0, 18, 18, 51, 49, 100, 100, 208]])
                # fmt: on
            )
        )

    def test_build3(self):
        # use Yeast mitochondrial codon table
        codon_table = CodonTable.unambiguous_dna_by_id[3]
        aligner = CodonAligner(codon_table=codon_table)
        codon_alignments = []
        seq1 = SeqRecord(
            Seq(
                "ATGGCAAGGGACCACCCAGTTGGGCACTGATATGATCGGGTGTATTTGCAGAGTAGTAACCTTTCTTTTCTCAAGACCATCCAG"
            ),
            id="pro1",
        )
        seq2 = SeqRecord(
            Seq(
                "ATGGCAAGGCACCATCCAGTTGAGCACTGATATGATCGGGTGTATTTGCAGAGTAGTAACGTGTCTCTGCTCAAGACCATCCAG"
            ),
            id="pro2",
        )
        seq3 = SeqRecord(
            Seq(
                "ATGGCAGGGGACCACCCAGTTGGGCACTGATATGATCGTGTGTATCTGCAGAGTAGTAACCACTCTTTTCTCATGACCATCCAG"
            ),
            id="pro3",
        )
        pro1 = SeqRecord(Seq("MARDHPVGHWYDRVYLQSSNTSFTKTIQ"), id="pro1")
        pro2 = SeqRecord(Seq("MARHHPVEHWYDRVYLQSSNVSTTKTIQ"), id="pro2")
        pro3 = SeqRecord(Seq("MAGDHPVGHWYDRVYTQSSNHSFTMTIQ"), id="pro3")
        alignments = aligner.align(pro1, seq1)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 28], [0, 84]]))
        )
        self.assertEqual(
            str(alignment),
            """\
pro1              0 M  A  R  D  H  P  V  G  H  W  Y  D  R  V  Y  L  Q  S  S  N  
pro1              0 ATGGCAAGGGACCACCCAGTTGGGCACTGATATGATCGGGTGTATTTGCAGAGTAGTAAC

pro1             20 T  S  F  T  K  T  I  Q   28
pro1             60 CTTTCTTTTCTCAAGACCATCCAG 84
""",
        )
        alignments = aligner.align(pro2, seq2)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 28], [0, 84]]))
        )
        self.assertEqual(
            str(alignment),
            """\
pro2              0 M  A  R  H  H  P  V  E  H  W  Y  D  R  V  Y  L  Q  S  S  N  
pro2              0 ATGGCAAGGCACCATCCAGTTGAGCACTGATATGATCGGGTGTATTTGCAGAGTAGTAAC

pro2             20 V  S  T  T  K  T  I  Q   28
pro2             60 GTGTCTCTGCTCAAGACCATCCAG 84
""",
        )
        alignments = aligner.align(pro3, seq3)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 28], [0, 84]]))
        )
        self.assertEqual(
            str(alignment),
            """\
pro3              0 M  A  G  D  H  P  V  G  H  W  Y  D  R  V  Y  T  Q  S  S  N  
pro3              0 ATGGCAGGGGACCACCCAGTTGGGCACTGATATGATCGTGTGTATCTGCAGAGTAGTAAC

pro3             20 H  S  F  T  M  T  I  Q   28
pro3             60 CACTCTTTTCTCATGACCATCCAG 84
""",
        )
        alignment = Alignment([pro1, pro2, pro3])
        alignment = alignment.mapall(codon_alignments)
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[0, 84],
                          [0, 84],
                          [0, 84]])
                # fmt: on
            )
        )
        self.assertEqual(
            str(alignment),
            """\
pro1              0 ATGGCAAGGGACCACCCAGTTGGGCACTGATATGATCGGGTGTATTTGCAGAGTAGTAAC
pro2              0 ATGGCAAGGCACCATCCAGTTGAGCACTGATATGATCGGGTGTATTTGCAGAGTAGTAAC
pro3              0 ATGGCAGGGGACCACCCAGTTGGGCACTGATATGATCGTGTGTATCTGCAGAGTAGTAAC

pro1             60 CTTTCTTTTCTCAAGACCATCCAG 84
pro2             60 GTGTCTCTGCTCAAGACCATCCAG 84
pro3             60 CACTCTTTTCTCATGACCATCCAG 84
""",
        )


class Test_dn_ds(unittest.TestCase):
    def test_dn_ds(self):
        aligner = CodonAligner()
        nucleotide_records = SeqIO.index("codonalign/egfr_nucl.fa", "fasta")
        protein_alignment = Align.read("codonalign/egfr_pro.aln", "clustal")
        self.assertEqual(len(protein_alignment.sequences), 6)
        codon_alignments = []
        protein_record = protein_alignment.sequences[0]
        self.assertEqual(protein_record.id, "gi|17136534|ref|NP_476758.1|")
        nucleotide_record = nucleotide_records["gi|24657088|ref|NM_057410.3|"]
        alignments = aligner.align(protein_record, nucleotide_record)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 1377], [84, 4215]]))
        )
        self.assertEqual(
            str(alignment),
            """\
gi|171365         0 M  M  I  I  S  M  W  M  S  I  S  R  G  L  W  D  S  S  S  I  
gi|246570        84 ATGATGATTATCAGCATGTGGATGAGCATATCGCGAGGATTGTGGGACAGCAGCTCCATC

gi|171365        20 W  S  V  L  L  I  L  A  C  M  A  S  I  T  T  S  S  S  V  S  
gi|246570       144 TGGTCGGTGCTGCTGATCCTCGCCTGCATGGCATCCATCACCACAAGCTCATCGGTCAGC

gi|171365        40 N  A  G  Y  V  D  N  G  N  M  K  V  C  I  G  T  K  S  R  L  
gi|246570       204 AATGCCGGCTATGTGGATAATGGCAATATGAAAGTCTGCATCGGCACTAAATCTCGGCTC

gi|171365        60 S  V  P  S  N  K  E  H  H  Y  R  N  L  R  D  R  Y  T  N  C  
gi|246570       264 TCCGTGCCCTCCAACAAGGAACATCATTACCGGAACCTCAGAGATCGGTACACGAACTGT

gi|171365        80 T  Y  V  D  G  N  L  E  L  T  W  L  P  N  E  N  L  D  L  S  
gi|246570       324 ACGTATGTGGATGGCAACCTGGAGCTGACCTGGCTGCCCAACGAGAATTTGGACCTCAGC

gi|171365       100 F  L  D  N  I  R  E  V  T  G  Y  I  L  I  S  H  V  D  V  K  
gi|246570       384 TTCCTGGACAACATACGGGAGGTCACCGGCTATATTCTGATCAGTCATGTGGACGTTAAG

gi|171365       120 K  V  V  F  P  K  L  Q  I  I  R  G  R  T  L  F  S  L  S  V  
gi|246570       444 AAAGTGGTATTTCCCAAACTACAAATCATTCGCGGACGCACGCTGTTCAGCTTATCCGTG

gi|171365       140 E  E  E  K  Y  A  L  F  V  T  Y  S  K  M  Y  T  L  E  I  P  
gi|246570       504 GAGGAGGAGAAGTATGCCTTGTTCGTCACTTATTCCAAAATGTACACGCTGGAGATTCCC

gi|171365       160 D  L  R  D  V  L  N  G  Q  V  G  F  H  N  N  Y  N  L  C  H  
gi|246570       564 GATCTACGCGATGTCTTAAATGGCCAAGTGGGCTTCCACAACAACTACAATCTCTGCCAC

gi|171365       180 M  R  T  I  Q  W  S  E  I  V  S  N  G  T  D  A  Y  Y  N  Y  
gi|246570       624 ATGCGAACGATCCAGTGGTCGGAGATTGTATCCAACGGCACGGATGCATACTACAACTAC

gi|171365       200 D  F  T  A  P  E  R  E  C  P  K  C  H  E  S  C  T  H  G  C  
gi|246570       684 GACTTTACTGCTCCGGAGCGCGAGTGTCCCAAGTGCCACGAGAGCTGCACGCACGGATGT

gi|171365       220 W  G  E  G  P  K  N  C  Q  K  F  S  K  L  T  C  S  P  Q  C  
gi|246570       744 TGGGGCGAGGGTCCCAAGAATTGCCAGAAGTTCAGCAAGCTCACCTGCTCGCCACAGTGT

gi|171365       240 A  G  G  R  C  Y  G  P  K  P  R  E  C  C  H  L  F  C  A  G  
gi|246570       804 GCCGGAGGTCGTTGCTATGGACCAAAGCCGCGGGAGTGTTGTCACCTCTTCTGCGCCGGA

gi|171365       260 G  C  T  G  P  T  Q  K  D  C  I  A  C  K  N  F  F  D  E  G  
gi|246570       864 GGATGCACTGGTCCCACGCAAAAGGATTGCATCGCCTGCAAGAACTTCTTCGACGAGGGC

gi|171365       280 V  C  K  E  E  C  P  P  M  R  K  Y  N  P  T  T  Y  V  L  E  
gi|246570       924 GTATGCAAGGAGGAATGCCCGCCCATGCGCAAGTACAATCCCACCACCTATGTTCTTGAA

gi|171365       300 T  N  P  E  G  K  Y  A  Y  G  A  T  C  V  K  E  C  P  G  H  
gi|246570       984 ACGAATCCTGAGGGAAAGTATGCCTATGGTGCCACCTGCGTCAAGGAGTGTCCCGGTCAT

gi|171365       320 L  L  R  D  N  G  A  C  V  R  S  C  P  Q  D  K  M  D  K  G  
gi|246570      1044 CTGTTGCGTGATAATGGCGCCTGCGTGCGCAGCTGTCCCCAGGACAAGATGGACAAGGGG

gi|171365       340 G  E  C  V  P  C  N  G  P  C  P  K  T  C  P  G  V  T  V  L  
gi|246570      1104 GGCGAGTGTGTGCCCTGCAATGGACCGTGCCCCAAAACCTGCCCGGGCGTTACTGTCCTG

gi|171365       360 H  A  G  N  I  D  S  F  R  N  C  T  V  I  D  G  N  I  R  I  
gi|246570      1164 CATGCCGGCAACATTGACTCGTTCCGGAATTGTACGGTGATCGATGGCAACATTCGCATT

gi|171365       380 L  D  Q  T  F  S  G  F  Q  D  V  Y  A  N  Y  T  M  G  P  R  
gi|246570      1224 TTGGATCAGACCTTCTCGGGCTTCCAGGATGTCTATGCCAACTACACGATGGGACCACGA

gi|171365       400 Y  I  P  L  D  P  E  R  L  E  V  F  S  T  V  K  E  I  T  G  
gi|246570      1284 TACATACCGCTGGATCCCGAGCGACTGGAGGTGTTCTCCACGGTGAAGGAGATCACCGGG

gi|171365       420 Y  L  N  I  E  G  T  H  P  Q  F  R  N  L  S  Y  F  R  N  L  
gi|246570      1344 TATCTGAATATCGAGGGAACCCACCCGCAGTTCCGGAATCTGTCGTACTTCCGCAATCTG

gi|171365       440 E  T  I  H  G  R  Q  L  M  E  S  M  F  A  A  L  A  I  V  K  
gi|246570      1404 GAAACAATTCATGGCCGCCAGCTGATGGAGAGCATGTTTGCCGCTTTGGCGATCGTTAAG

gi|171365       460 S  S  L  Y  S  L  E  M  R  N  L  K  Q  I  S  S  G  S  V  V  
gi|246570      1464 TCATCCCTGTACAGCCTGGAGATGCGCAATCTGAAGCAGATTAGTTCCGGCAGTGTGGTC

gi|171365       480 I  Q  H  N  R  D  L  C  Y  V  S  N  I  R  W  P  A  I  Q  K  
gi|246570      1524 ATCCAGCATAATAGAGACCTCTGCTACGTAAGCAATATCCGTTGGCCGGCCATTCAGAAG

gi|171365       500 E  P  E  Q  K  V  W  V  N  E  N  L  R  A  D  L  C  E  K  N  
gi|246570      1584 GAGCCCGAACAGAAGGTGTGGGTCAACGAGAATCTCAGGGCGGATCTATGCGAGAAAAAT

gi|171365       520 G  T  I  C  S  D  Q  C  N  E  D  G  C  W  G  A  G  T  D  Q  
gi|246570      1644 GGAACCATTTGCTCGGATCAGTGCAACGAGGACGGCTGCTGGGGAGCTGGCACGGATCAG

gi|171365       540 C  L  T  C  K  N  F  N  F  N  G  T  C  I  A  D  C  G  Y  I  
gi|246570      1704 TGCCTTACCTGCAAGAACTTCAATTTCAATGGCACCTGCATCGCCGACTGTGGTTATATA

gi|171365       560 S  N  A  Y  K  F  D  N  R  T  C  K  I  C  H  P  E  C  R  T  
gi|246570      1764 TCCAATGCCTACAAGTTTGACAATAGAACGTGCAAGATATGCCATCCAGAGTGCCGGACT

gi|171365       580 C  N  G  A  G  A  D  H  C  Q  E  C  V  H  V  R  D  G  Q  H  
gi|246570      1824 TGCAATGGAGCTGGAGCAGATCACTGCCAGGAGTGCGTCCATGTGAGGGACGGTCAGCAC

gi|171365       600 C  V  S  E  C  P  K  N  K  Y  N  D  R  G  V  C  R  E  C  H  
gi|246570      1884 TGTGTGTCCGAGTGCCCGAAGAACAAGTACAACGATCGTGGTGTCTGCCGAGAGTGCCAC

gi|171365       620 A  T  C  D  G  C  T  G  P  K  D  T  I  G  I  G  A  C  T  T  
gi|246570      1944 GCCACCTGCGATGGATGCACTGGGCCCAAGGACACCATCGGCATTGGAGCGTGTACAACG

gi|171365       640 C  N  L  A  I  I  N  N  D  A  T  V  K  R  C  L  L  K  D  D  
gi|246570      2004 TGCAATTTGGCCATTATCAACAATGACGCCACAGTAAAACGCTGCCTGCTGAAGGACGAC

gi|171365       660 K  C  P  D  G  Y  F  W  E  Y  V  H  P  Q  E  Q  G  S  L  K  
gi|246570      2064 AAGTGCCCCGATGGGTACTTCTGGGAGTATGTGCATCCACAAGAGCAGGGATCGCTAAAG

gi|171365       680 P  L  A  G  R  A  V  C  R  K  C  H  P  L  C  E  L  C  T  N  
gi|246570      2124 CCATTGGCCGGCAGAGCAGTTTGCCGAAAGTGCCATCCCCTTTGCGAGCTGTGCACCAAC

gi|171365       700 Y  G  Y  H  E  Q  V  C  S  K  C  T  H  Y  K  R  R  E  Q  C  
gi|246570      2184 TACGGATACCATGAACAGGTGTGCTCCAAGTGCACCCACTACAAGCGACGAGAGCAGTGC

gi|171365       720 E  T  E  C  P  A  D  H  Y  T  D  E  E  Q  R  E  C  F  Q  C  
gi|246570      2244 GAGACCGAGTGTCCGGCCGATCACTACACGGATGAGGAGCAGCGCGAGTGCTTCCAGTGC

gi|171365       740 H  P  E  C  N  G  C  T  G  P  G  A  D  D  C  K  S  C  R  N  
gi|246570      2304 CACCCAGAATGCAACGGTTGCACTGGTCCGGGTGCCGACGATTGCAAGTCTTGTCGCAAC

gi|171365       760 F  K  L  F  D  A  N  E  T  G  P  Y  V  N  S  T  M  F  N  C  
gi|246570      2364 TTCAAGTTGTTCGACGCGAATGAGACGGGTCCCTATGTGAACTCCACGATGTTCAATTGC

gi|171365       780 T  S  K  C  P  L  E  M  R  H  V  N  Y  Q  Y  T  A  I  G  P  
gi|246570      2424 ACCTCGAAGTGTCCCTTGGAGATGCGACATGTGAACTATCAGTACACGGCCATTGGACCC

gi|171365       800 Y  C  A  A  S  P  P  R  S  S  K  I  T  A  N  L  D  V  N  M  
gi|246570      2484 TACTGTGCAGCTAGTCCGCCGAGGAGCAGCAAGATAACTGCCAATCTGGATGTGAACATG

gi|171365       820 I  F  I  I  T  G  A  V  L  V  P  T  I  C  I  L  C  V  V  T  
gi|246570      2544 ATCTTCATTATCACTGGTGCTGTTCTGGTGCCGACGATCTGCATCCTCTGCGTGGTCACA

gi|171365       840 Y  I  C  R  Q  K  Q  K  A  K  K  E  T  V  K  M  T  M  A  L  
gi|246570      2604 TACATTTGTCGGCAAAAGCAAAAGGCCAAGAAAGAAACAGTGAAGATGACCATGGCTCTG

gi|171365       860 S  G  C  E  D  S  E  P  L  R  P  S  N  I  G  A  N  L  C  K  
gi|246570      2664 TCCGGCTGTGAGGATTCCGAGCCGCTGCGTCCCTCGAACATTGGAGCCAATCTATGCAAG

gi|171365       880 L  R  I  V  K  D  A  E  L  R  K  G  G  V  L  G  M  G  A  F  
gi|246570      2724 TTGCGCATTGTCAAGGACGCCGAGTTGCGCAAGGGCGGAGTCCTCGGAATGGGAGCCTTT

gi|171365       900 G  R  V  Y  K  G  V  W  V  P  E  G  E  N  V  K  I  P  V  A  
gi|246570      2784 GGACGAGTGTACAAGGGCGTTTGGGTGCCGGAGGGTGAGAACGTCAAGATTCCAGTGGCC

gi|171365       920 I  K  E  L  L  K  S  T  G  A  E  S  S  E  E  F  L  R  E  A  
gi|246570      2844 ATTAAGGAGCTGCTCAAGTCCACAGGCGCCGAGTCAAGCGAAGAGTTCCTCCGCGAAGCC

gi|171365       940 Y  I  M  A  S  V  E  H  V  N  L  L  K  L  L  A  V  C  M  S  
gi|246570      2904 TACATCATGGCCTCTGTGGAGCACGTTAATCTGCTGAAGCTCCTGGCCGTCTGCATGTCC

gi|171365       960 S  Q  M  M  L  I  T  Q  L  M  P  L  G  C  L  L  D  Y  V  R  
gi|246570      2964 TCACAAATGATGCTAATCACGCAACTGATGCCGCTTGGCTGCCTGTTGGACTATGTGCGA

gi|171365       980 N  N  R  D  K  I  G  S  K  A  L  L  N  W  S  T  Q  I  A  K  
gi|246570      3024 AATAACCGGGACAAGATCGGCTCTAAGGCTCTGCTCAACTGGAGCACGCAAATCGCCAAG

gi|171365      1000 G  M  S  Y  L  E  E  K  R  L  V  H  R  D  L  A  A  R  N  V  
gi|246570      3084 GGCATGTCGTATCTGGAGGAGAAGCGACTGGTCCACAGAGACTTGGCTGCCCGCAATGTC

gi|171365      1020 L  V  Q  T  P  S  L  V  K  I  T  D  F  G  L  A  K  L  L  S  
gi|246570      3144 CTGGTGCAGACTCCCTCGCTGGTGAAGATCACCGACTTTGGGCTGGCCAAGTTGCTGAGC

gi|171365      1040 S  D  S  N  E  Y  K  A  A  G  G  K  M  P  I  K  W  L  A  L  
gi|246570      3204 AGCGATTCCAATGAGTACAAGGCTGCTGGCGGCAAGATGCCCATCAAGTGGTTGGCACTG

gi|171365      1060 E  C  I  R  N  R  V  F  T  S  K  S  D  V  W  A  F  G  V  T  
gi|246570      3264 GAGTGCATTCGCAATCGTGTATTCACCAGCAAGTCCGATGTCTGGGCCTTTGGTGTGACA

gi|171365      1080 I  W  E  L  L  T  F  G  Q  R  P  H  E  N  I  P  A  K  D  I  
gi|246570      3324 ATTTGGGAACTGCTGACCTTTGGCCAGCGTCCACACGAGAACATCCCCGCTAAGGATATT

gi|171365      1100 P  D  L  I  E  V  G  L  K  L  E  Q  P  E  I  C  S  L  D  I  
gi|246570      3384 CCCGATCTTATTGAAGTCGGTCTGAAGCTGGAGCAGCCGGAGATTTGTTCGCTGGACATT

gi|171365      1120 Y  C  T  L  L  S  C  W  H  L  D  A  A  M  R  P  T  F  K  Q  
gi|246570      3444 TACTGCACACTTCTCTCGTGCTGGCACTTGGATGCCGCCATGCGTCCAACCTTCAAGCAG

gi|171365      1140 L  T  T  V  F  A  E  F  A  R  D  P  G  R  Y  L  A  I  P  G  
gi|246570      3504 CTGACTACGGTCTTTGCTGAGTTCGCCAGAGATCCGGGTCGCTATCTGGCCATTCCCGGG

gi|171365      1160 D  K  F  T  R  L  P  A  Y  T  S  Q  D  E  K  D  L  I  R  K  
gi|246570      3564 GATAAGTTCACCCGGCTGCCGGCCTACACGAGTCAGGATGAGAAGGATCTCATCCGAAAA

gi|171365      1180 L  A  P  T  T  D  G  S  E  A  I  A  E  P  D  D  Y  L  Q  P  
gi|246570      3624 TTGGCTCCCACCACCGATGGGTCCGAAGCCATTGCGGAACCCGATGACTACCTGCAACCC

gi|171365      1200 K  A  A  P  G  P  S  H  R  T  D  C  T  D  E  I  P  K  L  N  
gi|246570      3684 AAGGCAGCACCTGGTCCTAGTCACAGAACCGACTGCACGGATGAGATACCCAAGCTGAAC

gi|171365      1220 R  Y  C  K  D  P  S  N  K  N  S  S  T  G  D  D  E  T  D  S  
gi|246570      3744 CGCTACTGCAAGGATCCTAGCAACAAGAATTCGAGTACCGGAGACGATGAGACGGATTCG

gi|171365      1240 S  A  R  E  V  G  V  G  N  L  R  L  D  L  P  V  D  E  D  D  
gi|246570      3804 AGTGCCCGGGAAGTGGGCGTGGGTAATCTGCGCCTCGATCTACCAGTCGATGAGGATGAT

gi|171365      1260 Y  L  M  P  T  C  Q  P  G  P  N  N  N  N  N  I  N  N  P  N  
gi|246570      3864 TACCTGATGCCCACATGCCAACCGGGGCCCAACAACAACAACAACATAAATAATCCCAAT

gi|171365      1280 Q  N  N  M  A  A  V  G  V  A  A  G  Y  M  D  L  I  G  V  P  
gi|246570      3924 CAAAACAATATGGCAGCTGTGGGCGTGGCTGCCGGCTACATGGATCTCATCGGAGTGCCC

gi|171365      1300 V  S  V  D  N  P  E  Y  L  L  N  A  Q  T  L  G  V  G  E  S  
gi|246570      3984 GTTAGTGTGGACAATCCGGAGTATCTGCTAAACGCGCAGACACTGGGTGTTGGGGAGTCG

gi|171365      1320 P  I  P  T  Q  T  I  G  I  P  V  M  G  V  P  G  T  M  E  V  
gi|246570      4044 CCGATACCCACCCAGACCATCGGGATACCGGTGATGGGAGTCCCGGGCACCATGGAGGTC

gi|171365      1340 K  V  P  M  P  G  S  E  P  T  S  S  D  H  E  Y  Y  N  D  T  
gi|246570      4104 AAGGTGCCAATGCCAGGCAGTGAGCCAACGAGCTCCGATCACGAGTACTACAATGATACC

gi|171365      1360 Q  R  E  L  Q  P  L  H  R  N  R  N  T  E  T  R  V   1377
gi|246570      4164 CAACGGGAGTTGCAGCCACTGCATCGAAACCGCAACACGGAGACGAGGGTG 4215
""",
        )
        protein_record = protein_alignment.sequences[1]
        self.assertEqual(protein_record.id, "gi|17136536|ref|NP_476759.1|")
        nucleotide_record = nucleotide_records["gi|24657104|ref|NM_057411.3|"]
        alignments = aligner.align(protein_record, nucleotide_record)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 1426], [22, 4300]]))
        )
        self.assertEqual(
            str(alignment),
            """\
gi|171365         0 M  L  L  R  R  R  N  G  P  C  P  F  P  L  L  L  L  L  L  A  
gi|246571        22 ATGCTGCTGCGACGGCGCAACGGCCCCTGCCCCTTCCCCCTGCTGCTCCTGCTCCTGGCC

gi|171365        20 H  C  I  C  I  W  P  A  S  A  A  R  D  R  Y  A  R  Q  N  N  
gi|246571        82 CACTGCATTTGCATTTGGCCCGCGTCGGCGGCCCGCGATCGCTACGCCCGCCAGAACAAT

gi|171365        40 R  Q  R  H  Q  D  I  D  R  D  R  D  R  D  R  F  L  Y  R  S  
gi|246571       142 CGCCAGCGCCATCAGGATATAGATCGCGATCGGGATCGAGATCGATTCCTATACCGCAGC

gi|171365        60 S  S  A  Q  N  R  Q  R  G  G  A  N  F  A  L  G  L  G  A  N  
gi|246571       202 AGTTCGGCCCAAAATCGACAGAGGGGCGGGGCCAACTTCGCCCTGGGACTGGGAGCCAAC

gi|171365        80 G  V  T  I  P  T  S  L  E  D  K  N  K  N  E  F  V  K  G  K  
gi|246571       262 GGAGTCACCATTCCCACCAGTCTGGAGGATAAGAACAAGAACGAGTTCGTCAAGGGGAAA

gi|171365       100 I  C  I  G  T  K  S  R  L  S  V  P  S  N  K  E  H  H  Y  R  
gi|246571       322 ATCTGCATCGGCACTAAATCTCGGCTCTCCGTGCCCTCCAACAAGGAACATCATTACCGG

gi|171365       120 N  L  R  D  R  Y  T  N  C  T  Y  V  D  G  N  L  E  L  T  W  
gi|246571       382 AACCTCAGAGATCGGTACACGAACTGTACGTATGTGGATGGCAACCTGGAGCTGACCTGG

gi|171365       140 L  P  N  E  N  L  D  L  S  F  L  D  N  I  R  E  V  T  G  Y  
gi|246571       442 CTGCCCAACGAGAATTTGGACCTCAGCTTCCTGGACAACATACGGGAGGTCACCGGCTAT

gi|171365       160 I  L  I  S  H  V  D  V  K  K  V  V  F  P  K  L  Q  I  I  R  
gi|246571       502 ATTCTGATCAGTCATGTGGACGTTAAGAAAGTGGTATTTCCCAAACTACAAATCATTCGC

gi|171365       180 G  R  T  L  F  S  L  S  V  E  E  E  K  Y  A  L  F  V  T  Y  
gi|246571       562 GGACGCACGCTGTTCAGCTTATCCGTGGAGGAGGAGAAGTATGCCTTGTTCGTCACTTAT

gi|171365       200 S  K  M  Y  T  L  E  I  P  D  L  R  D  V  L  N  G  Q  V  G  
gi|246571       622 TCCAAAATGTACACGCTGGAGATTCCCGATCTACGCGATGTCTTAAATGGCCAAGTGGGC

gi|171365       220 F  H  N  N  Y  N  L  C  H  M  R  T  I  Q  W  S  E  I  V  S  
gi|246571       682 TTCCACAACAACTACAATCTCTGCCACATGCGAACGATCCAGTGGTCGGAGATTGTATCC

gi|171365       240 N  G  T  D  A  Y  Y  N  Y  D  F  T  A  P  E  R  E  C  P  K  
gi|246571       742 AACGGCACGGATGCATACTACAACTACGACTTTACTGCTCCGGAGCGCGAGTGTCCCAAG

gi|171365       260 C  H  E  S  C  T  H  G  C  W  G  E  G  P  K  N  C  Q  K  F  
gi|246571       802 TGCCACGAGAGCTGCACGCACGGATGTTGGGGCGAGGGTCCCAAGAATTGCCAGAAGTTC

gi|171365       280 S  K  L  T  C  S  P  Q  C  A  G  G  R  C  Y  G  P  K  P  R  
gi|246571       862 AGCAAGCTCACCTGCTCGCCACAGTGTGCCGGAGGTCGTTGCTATGGACCAAAGCCGCGG

gi|171365       300 E  C  C  H  L  F  C  A  G  G  C  T  G  P  T  Q  K  D  C  I  
gi|246571       922 GAGTGTTGTCACCTCTTCTGCGCCGGAGGATGCACTGGTCCCACGCAAAAGGATTGCATC

gi|171365       320 A  C  K  N  F  F  D  E  G  V  C  K  E  E  C  P  P  M  R  K  
gi|246571       982 GCCTGCAAGAACTTCTTCGACGAGGGCGTATGCAAGGAGGAATGCCCGCCCATGCGCAAG

gi|171365       340 Y  N  P  T  T  Y  V  L  E  T  N  P  E  G  K  Y  A  Y  G  A  
gi|246571      1042 TACAATCCCACCACCTATGTTCTTGAAACGAATCCTGAGGGAAAGTATGCCTATGGTGCC

gi|171365       360 T  C  V  K  E  C  P  G  H  L  L  R  D  N  G  A  C  V  R  S  
gi|246571      1102 ACCTGCGTCAAGGAGTGTCCCGGTCATCTGTTGCGTGATAATGGCGCCTGCGTGCGCAGC

gi|171365       380 C  P  Q  D  K  M  D  K  G  G  E  C  V  P  C  N  G  P  C  P  
gi|246571      1162 TGTCCCCAGGACAAGATGGACAAGGGGGGCGAGTGTGTGCCCTGCAATGGACCGTGCCCC

gi|171365       400 K  T  C  P  G  V  T  V  L  H  A  G  N  I  D  S  F  R  N  C  
gi|246571      1222 AAAACCTGCCCGGGCGTTACTGTCCTGCATGCCGGCAACATTGACTCGTTCCGGAATTGT

gi|171365       420 T  V  I  D  G  N  I  R  I  L  D  Q  T  F  S  G  F  Q  D  V  
gi|246571      1282 ACGGTGATCGATGGCAACATTCGCATTTTGGATCAGACCTTCTCGGGCTTCCAGGATGTC

gi|171365       440 Y  A  N  Y  T  M  G  P  R  Y  I  P  L  D  P  E  R  L  E  V  
gi|246571      1342 TATGCCAACTACACGATGGGACCACGATACATACCGCTGGATCCCGAGCGACTGGAGGTG

gi|171365       460 F  S  T  V  K  E  I  T  G  Y  L  N  I  E  G  T  H  P  Q  F  
gi|246571      1402 TTCTCCACGGTGAAGGAGATCACCGGGTATCTGAATATCGAGGGAACCCACCCGCAGTTC

gi|171365       480 R  N  L  S  Y  F  R  N  L  E  T  I  H  G  R  Q  L  M  E  S  
gi|246571      1462 CGGAATCTGTCGTACTTCCGCAATCTGGAAACAATTCATGGCCGCCAGCTGATGGAGAGC

gi|171365       500 M  F  A  A  L  A  I  V  K  S  S  L  Y  S  L  E  M  R  N  L  
gi|246571      1522 ATGTTTGCCGCTTTGGCGATCGTTAAGTCATCCCTGTACAGCCTGGAGATGCGCAATCTG

gi|171365       520 K  Q  I  S  S  G  S  V  V  I  Q  H  N  R  D  L  C  Y  V  S  
gi|246571      1582 AAGCAGATTAGTTCCGGCAGTGTGGTCATCCAGCATAATAGAGACCTCTGCTACGTAAGC

gi|171365       540 N  I  R  W  P  A  I  Q  K  E  P  E  Q  K  V  W  V  N  E  N  
gi|246571      1642 AATATCCGTTGGCCGGCCATTCAGAAGGAGCCCGAACAGAAGGTGTGGGTCAACGAGAAT

gi|171365       560 L  R  A  D  L  C  E  K  N  G  T  I  C  S  D  Q  C  N  E  D  
gi|246571      1702 CTCAGGGCGGATCTATGCGAGAAAAATGGAACCATTTGCTCGGATCAGTGCAACGAGGAC

gi|171365       580 G  C  W  G  A  G  T  D  Q  C  L  T  C  K  N  F  N  F  N  G  
gi|246571      1762 GGCTGCTGGGGAGCTGGCACGGATCAGTGCCTTACCTGCAAGAACTTCAATTTCAATGGC

gi|171365       600 T  C  I  A  D  C  G  Y  I  S  N  A  Y  K  F  D  N  R  T  C  
gi|246571      1822 ACCTGCATCGCCGACTGTGGTTATATATCCAATGCCTACAAGTTTGACAATAGAACGTGC

gi|171365       620 K  I  C  H  P  E  C  R  T  C  N  G  A  G  A  D  H  C  Q  E  
gi|246571      1882 AAGATATGCCATCCAGAGTGCCGGACTTGCAATGGAGCTGGAGCAGATCACTGCCAGGAG

gi|171365       640 C  V  H  V  R  D  G  Q  H  C  V  S  E  C  P  K  N  K  Y  N  
gi|246571      1942 TGCGTCCATGTGAGGGACGGTCAGCACTGTGTGTCCGAGTGCCCGAAGAACAAGTACAAC

gi|171365       660 D  R  G  V  C  R  E  C  H  A  T  C  D  G  C  T  G  P  K  D  
gi|246571      2002 GATCGTGGTGTCTGCCGAGAGTGCCACGCCACCTGCGATGGATGCACTGGGCCCAAGGAC

gi|171365       680 T  I  G  I  G  A  C  T  T  C  N  L  A  I  I  N  N  D  A  T  
gi|246571      2062 ACCATCGGCATTGGAGCGTGTACAACGTGCAATTTGGCCATTATCAACAATGACGCCACA

gi|171365       700 V  K  R  C  L  L  K  D  D  K  C  P  D  G  Y  F  W  E  Y  V  
gi|246571      2122 GTAAAACGCTGCCTGCTGAAGGACGACAAGTGCCCCGATGGGTACTTCTGGGAGTATGTG

gi|171365       720 H  P  Q  E  Q  G  S  L  K  P  L  A  G  R  A  V  C  R  K  C  
gi|246571      2182 CATCCACAAGAGCAGGGATCGCTAAAGCCATTGGCCGGCAGAGCAGTTTGCCGAAAGTGC

gi|171365       740 H  P  L  C  E  L  C  T  N  Y  G  Y  H  E  Q  V  C  S  K  C  
gi|246571      2242 CATCCCCTTTGCGAGCTGTGCACCAACTACGGATACCATGAACAGGTGTGCTCCAAGTGC

gi|171365       760 T  H  Y  K  R  R  E  Q  C  E  T  E  C  P  A  D  H  Y  T  D  
gi|246571      2302 ACCCACTACAAGCGACGAGAGCAGTGCGAGACCGAGTGTCCGGCCGATCACTACACGGAT

gi|171365       780 E  E  Q  R  E  C  F  Q  C  H  P  E  C  N  G  C  T  G  P  G  
gi|246571      2362 GAGGAGCAGCGCGAGTGCTTCCAGTGCCACCCAGAATGCAACGGTTGCACTGGTCCGGGT

gi|171365       800 A  D  D  C  K  S  C  R  N  F  K  L  F  D  A  N  E  T  G  P  
gi|246571      2422 GCCGACGATTGCAAGTCTTGTCGCAACTTCAAGTTGTTCGACGCGAATGAGACGGGTCCC

gi|171365       820 Y  V  N  S  T  M  F  N  C  T  S  K  C  P  L  E  M  R  H  V  
gi|246571      2482 TATGTGAACTCCACGATGTTCAATTGCACCTCGAAGTGTCCCTTGGAGATGCGACATGTG

gi|171365       840 N  Y  Q  Y  T  A  I  G  P  Y  C  A  A  S  P  P  R  S  S  K  
gi|246571      2542 AACTATCAGTACACGGCCATTGGACCCTACTGTGCAGCTAGTCCGCCGAGGAGCAGCAAG

gi|171365       860 I  T  A  N  L  D  V  N  M  I  F  I  I  T  G  A  V  L  V  P  
gi|246571      2602 ATAACTGCCAATCTGGATGTGAACATGATCTTCATTATCACTGGTGCTGTTCTGGTGCCG

gi|171365       880 T  I  C  I  L  C  V  V  T  Y  I  C  R  Q  K  Q  K  A  K  K  
gi|246571      2662 ACGATCTGCATCCTCTGCGTGGTCACATACATTTGTCGGCAAAAGCAAAAGGCCAAGAAA

gi|171365       900 E  T  V  K  M  T  M  A  L  S  G  C  E  D  S  E  P  L  R  P  
gi|246571      2722 GAAACAGTGAAGATGACCATGGCTCTGTCCGGCTGTGAGGATTCCGAGCCGCTGCGTCCC

gi|171365       920 S  N  I  G  A  N  L  C  K  L  R  I  V  K  D  A  E  L  R  K  
gi|246571      2782 TCGAACATTGGAGCCAATCTATGCAAGTTGCGCATTGTCAAGGACGCCGAGTTGCGCAAG

gi|171365       940 G  G  V  L  G  M  G  A  F  G  R  V  Y  K  G  V  W  V  P  E  
gi|246571      2842 GGCGGAGTCCTCGGAATGGGAGCCTTTGGACGAGTGTACAAGGGCGTTTGGGTGCCGGAG

gi|171365       960 G  E  N  V  K  I  P  V  A  I  K  E  L  L  K  S  T  G  A  E  
gi|246571      2902 GGTGAGAACGTCAAGATTCCAGTGGCCATTAAGGAGCTGCTCAAGTCCACAGGCGCCGAG

gi|171365       980 S  S  E  E  F  L  R  E  A  Y  I  M  A  S  V  E  H  V  N  L  
gi|246571      2962 TCAAGCGAAGAGTTCCTCCGCGAAGCCTACATCATGGCCTCTGTGGAGCACGTTAATCTG

gi|171365      1000 L  K  L  L  A  V  C  M  S  S  Q  M  M  L  I  T  Q  L  M  P  
gi|246571      3022 CTGAAGCTCCTGGCCGTCTGCATGTCCTCACAAATGATGCTAATCACGCAACTGATGCCG

gi|171365      1020 L  G  C  L  L  D  Y  V  R  N  N  R  D  K  I  G  S  K  A  L  
gi|246571      3082 CTTGGCTGCCTGTTGGACTATGTGCGAAATAACCGGGACAAGATCGGCTCTAAGGCTCTG

gi|171365      1040 L  N  W  S  T  Q  I  A  K  G  M  S  Y  L  E  E  K  R  L  V  
gi|246571      3142 CTCAACTGGAGCACGCAAATCGCCAAGGGCATGTCGTATCTGGAGGAGAAGCGACTGGTC

gi|171365      1060 H  R  D  L  A  A  R  N  V  L  V  Q  T  P  S  L  V  K  I  T  
gi|246571      3202 CACAGAGACTTGGCTGCCCGCAATGTCCTGGTGCAGACTCCCTCGCTGGTGAAGATCACC

gi|171365      1080 D  F  G  L  A  K  L  L  S  S  D  S  N  E  Y  K  A  A  G  G  
gi|246571      3262 GACTTTGGGCTGGCCAAGTTGCTGAGCAGCGATTCCAATGAGTACAAGGCTGCTGGCGGC

gi|171365      1100 K  M  P  I  K  W  L  A  L  E  C  I  R  N  R  V  F  T  S  K  
gi|246571      3322 AAGATGCCCATCAAGTGGTTGGCACTGGAGTGCATTCGCAATCGTGTATTCACCAGCAAG

gi|171365      1120 S  D  V  W  A  F  G  V  T  I  W  E  L  L  T  F  G  Q  R  P  
gi|246571      3382 TCCGATGTCTGGGCCTTTGGTGTGACAATTTGGGAACTGCTGACCTTTGGCCAGCGTCCA

gi|171365      1140 H  E  N  I  P  A  K  D  I  P  D  L  I  E  V  G  L  K  L  E  
gi|246571      3442 CACGAGAACATCCCCGCTAAGGATATTCCCGATCTTATTGAAGTCGGTCTGAAGCTGGAG

gi|171365      1160 Q  P  E  I  C  S  L  D  I  Y  C  T  L  L  S  C  W  H  L  D  
gi|246571      3502 CAGCCGGAGATTTGTTCGCTGGACATTTACTGCACACTTCTCTCGTGCTGGCACTTGGAT

gi|171365      1180 A  A  M  R  P  T  F  K  Q  L  T  T  V  F  A  E  F  A  R  D  
gi|246571      3562 GCCGCCATGCGTCCAACCTTCAAGCAGCTGACTACGGTCTTTGCTGAGTTCGCCAGAGAT

gi|171365      1200 P  G  R  Y  L  A  I  P  G  D  K  F  T  R  L  P  A  Y  T  S  
gi|246571      3622 CCGGGTCGCTATCTGGCCATTCCCGGGGATAAGTTCACCCGGCTGCCGGCCTACACGAGT

gi|171365      1220 Q  D  E  K  D  L  I  R  K  L  A  P  T  T  D  G  S  E  A  I  
gi|246571      3682 CAGGATGAGAAGGATCTCATCCGAAAATTGGCTCCCACCACCGATGGGTCCGAAGCCATT

gi|171365      1240 A  E  P  D  D  Y  L  Q  P  K  A  A  P  G  P  S  H  R  T  D  
gi|246571      3742 GCGGAACCCGATGACTACCTGCAACCCAAGGCAGCACCTGGTCCTAGTCACAGAACCGAC

gi|171365      1260 C  T  D  E  I  P  K  L  N  R  Y  C  K  D  P  S  N  K  N  S  
gi|246571      3802 TGCACGGATGAGATACCCAAGCTGAACCGCTACTGCAAGGATCCTAGCAACAAGAATTCG

gi|171365      1280 S  T  G  D  D  E  T  D  S  S  A  R  E  V  G  V  G  N  L  R  
gi|246571      3862 AGTACCGGAGACGATGAGACGGATTCGAGTGCCCGGGAAGTGGGCGTGGGTAATCTGCGC

gi|171365      1300 L  D  L  P  V  D  E  D  D  Y  L  M  P  T  C  Q  P  G  P  N  
gi|246571      3922 CTCGATCTACCAGTCGATGAGGATGATTACCTGATGCCCACATGCCAACCGGGGCCCAAC

gi|171365      1320 N  N  N  N  I  N  N  P  N  Q  N  N  M  A  A  V  G  V  A  A  
gi|246571      3982 AACAACAACAACATAAATAATCCCAATCAAAACAATATGGCAGCTGTGGGCGTGGCTGCC

gi|171365      1340 G  Y  M  D  L  I  G  V  P  V  S  V  D  N  P  E  Y  L  L  N  
gi|246571      4042 GGCTACATGGATCTCATCGGAGTGCCCGTTAGTGTGGACAATCCGGAGTATCTGCTAAAC

gi|171365      1360 A  Q  T  L  G  V  G  E  S  P  I  P  T  Q  T  I  G  I  P  V  
gi|246571      4102 GCGCAGACACTGGGTGTTGGGGAGTCGCCGATACCCACCCAGACCATCGGGATACCGGTG

gi|171365      1380 M  G  V  P  G  T  M  E  V  K  V  P  M  P  G  S  E  P  T  S  
gi|246571      4162 ATGGGAGTCCCGGGCACCATGGAGGTCAAGGTGCCAATGCCAGGCAGTGAGCCAACGAGC

gi|171365      1400 S  D  H  E  Y  Y  N  D  T  Q  R  E  L  Q  P  L  H  R  N  R  
gi|246571      4222 TCCGATCACGAGTACTACAATGATACCCAACGGGAGTTGCAGCCACTGCATCGAAACCGC

gi|171365      1420 N  T  E  T  R  V   1426
gi|246571      4282 AACACGGAGACGAGGGTG 4300
""",
        )
        protein_record = protein_alignment.sequences[2]
        self.assertEqual(protein_record.id, "gi|302179501|gb|ADK98534.1|")
        nucleotide_record = nucleotide_records["gi|302179500|gb|HM749883.1|"]
        alignments = aligner.align(protein_record, nucleotide_record)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 1185], [0, 3555]]))
        )
        self.assertEqual(
            str(alignment),
            """\
gi|302179         0 M  K  K  H  E  L  L  C  Q  G  T  S  N  K  L  T  Q  L  G  T  
gi|302179         0 ATGAAAAAGCACGAGTTACTTTGCCAAGGGACAAGTAACAAGCTCACCCAGTTGGGCACT

gi|302179        20 F  E  D  H  F  L  S  L  Q  R  M  F  N  N  C  E  V  V  L  G  
gi|302179        60 TTTGAAGACCACTTTCTGAGCCTACAGAGGATGTTCAACAACTGTGAGGTGGTCCTTGGG

gi|302179        40 N  L  E  I  T  Y  M  Q  S  S  Y  N  L  S  F  L  K  T  I  Q  
gi|302179       120 AATTTGGAAATTACCTACATGCAGAGTAGTTACAACCTTTCTTTTCTCAAGACCATCCAG

gi|302179        60 E  V  A  G  Y  V  L  I  A  L  N  T  V  E  K  I  P  L  E  N  
gi|302179       180 GAGGTTGCCGGCTATGTACTCATTGCCCTCAACACAGTGGAGAAGATTCCGCTGGAAAAC

gi|302179        80 L  Q  I  I  R  G  N  V  L  Y  E  N  T  H  A  L  A  V  L  S  
gi|302179       240 CTGCAGATCATCCGAGGAAATGTGCTTTATGAAAACACCCATGCCTTAGCCGTCTTATCC

gi|302179       100 N  Y  G  A  N  K  T  G  L  R  E  L  P  L  R  N  L  Q  E  I  
gi|302179       300 AACTATGGAGCAAACAAAACCGGACTGAGGGAGCTGCCCTTGAGAAACTTACAGGAAATT

gi|302179       120 L  Q  G  A  V  R  F  S  N  N  P  V  L  C  N  V  E  T  I  Q  
gi|302179       360 CTGCAAGGTGCCGTGAGATTCAGCAACAACCCTGTCCTCTGCAACGTGGAGACCATCCAG

gi|302179       140 W  R  D  I  V  N  P  D  F  L  S  N  M  T  G  D  F  Q  N  Q  
gi|302179       420 TGGCGGGACATCGTCAACCCTGATTTTCTAAGCAACATGACAGGGGACTTTCAGAACCAG

gi|302179       160 Q  G  N  C  P  K  C  D  P  A  C  L  N  R  S  C  W  G  A  G  
gi|302179       480 CAGGGCAACTGCCCAAAGTGTGATCCAGCCTGTCTCAACAGAAGCTGCTGGGGTGCCGGG

gi|302179       180 E  E  N  C  Q  K  L  T  K  I  I  C  A  Q  Q  C  S  G  R  C  
gi|302179       540 GAGGAGAACTGTCAGAAATTGACCAAAATCATCTGTGCCCAGCAGTGTTCCGGGCGCTGC

gi|302179       200 R  G  R  S  P  S  D  C  C  H  N  Q  C  A  A  G  C  T  G  P  
gi|302179       600 CGTGGCAGGTCCCCCAGTGACTGCTGCCACAACCAGTGTGCCGCTGGCTGCACAGGGCCA

gi|302179       220 R  E  S  D  C  L  V  C  R  R  F  R  D  E  A  T  C  K  D  T  
gi|302179       660 CGGGAGAGCGACTGCCTGGTCTGCCGCAGGTTCCGTGATGAAGCCACCTGCAAGGACACG

gi|302179       240 C  P  P  L  M  L  Y  D  P  T  T  Y  E  M  K  V  N  P  L  G  
gi|302179       720 TGTCCGCCACTCATGCTCTATGACCCTACCACCTACGAAATGAAGGTCAACCCGCTGGGG

gi|302179       260 K  Y  S  F  G  A  T  C  V  K  K  C  P  R  N  Y  V  V  T  D  
gi|302179       780 AAGTACAGCTTTGGCGCCACCTGTGTCAAGAAGTGTCCCCGTAACTACGTGGTGACAGAC

gi|302179       280 H  G  S  C  V  R  A  C  S  S  D  S  Q  E  V  E  E  D  G  V  
gi|302179       840 CACGGCTCCTGCGTCCGCGCCTGCAGTTCTGACAGCCAGGAGGTAGAGGAAGACGGTGTC

gi|302179       300 R  K  C  K  K  C  D  G  P  C  G  K  V  C  N  G  I  G  I  G  
gi|302179       900 CGCAAGTGTAAAAAGTGTGACGGGCCTTGTGGCAAAGTTTGTAACGGAATAGGAATCGGT

gi|302179       320 E  F  K  D  T  L  S  I  N  A  T  N  I  K  H  F  R  N  C  T  
gi|302179       960 GAGTTTAAAGACACACTTTCCATAAATGCTACAAACATTAAACACTTCAGAAACTGCACA

gi|302179       340 S  I  S  G  D  L  H  I  L  P  V  A  F  R  G  D  S  F  T  R  
gi|302179      1020 TCCATCAGTGGAGATCTTCATATCCTGCCAGTAGCATTTAGGGGTGACTCCTTCACACGT

gi|302179       360 T  A  P  L  D  P  K  E  L  D  I  L  R  T  V  K  E  I  T  G  
gi|302179      1080 ACTGCACCTCTGGACCCGAAAGAACTGGACATTCTAAGAACTGTAAAAGAAATAACAGGG

gi|302179       380 F  L  L  I  Q  A  W  P  E  N  R  T  D  L  H  A  F  E  N  L  
gi|302179      1140 TTTTTGCTGATTCAGGCCTGGCCCGAAAACAGGACTGACCTCCATGCTTTTGAGAACCTG

gi|302179       400 E  I  I  R  G  R  T  K  Q  H  G  Q  F  S  L  A  V  V  G  L  
gi|302179      1200 GAAATCATACGTGGCAGAACGAAGCAGCATGGCCAGTTTTCTCTTGCGGTTGTCGGCCTG

gi|302179       420 D  I  T  S  L  G  L  R  S  L  K  E  I  S  D  G  D  V  I  I  
gi|302179      1260 GATATAACATCTTTGGGATTACGCTCCCTCAAGGAGATAAGTGATGGTGATGTGATAATT

gi|302179       440 S  G  N  R  N  L  C  Y  A  D  T  I  R  W  K  K  L  F  G  T  
gi|302179      1320 TCAGGAAATCGAAACTTGTGCTATGCAGATACAATACGCTGGAAAAAACTTTTTGGGACC

gi|302179       460 S  T  Q  K  T  K  I  L  N  N  R  S  E  K  Q  C  K  A  A  G  
gi|302179      1380 TCAACTCAGAAAACCAAAATTTTAAACAACAGGAGTGAAAAACAGTGCAAGGCCGCAGGC

gi|302179       480 H  I  C  H  P  L  C  S  S  E  G  C  W  G  P  G  P  K  Y  C  
gi|302179      1440 CACATCTGTCACCCGCTGTGCTCATCAGAGGGCTGCTGGGGACCGGGACCCAAATACTGC

gi|302179       500 M  S  C  Q  N  F  S  R  G  K  E  C  V  G  K  C  N  I  L  E  
gi|302179      1500 ATGTCCTGCCAGAACTTCAGTCGTGGCAAGGAGTGTGTGGGAAAGTGCAACATTCTAGAG

gi|302179       520 G  E  P  R  E  F  V  E  N  S  E  C  V  Q  C  H  P  E  C  L  
gi|302179      1560 GGAGAGCCCAGAGAATTCGTGGAGAACTCCGAGTGTGTGCAGTGCCATCCAGAATGCCTG

gi|302179       540 P  Q  A  M  N  V  T  C  T  G  R  G  P  G  N  C  V  K  C  A  
gi|302179      1620 CCCCAGGCCATGAACGTGACCTGCACTGGACGCGGACCAGGCAACTGTGTAAAGTGCGCC

gi|302179       560 H  Y  I  D  G  P  H  C  V  K  T  C  P  A  G  V  A  G  E  N  
gi|302179      1680 CACTACATTGATGGCCCTCACTGCGTCAAGACCTGCCCTGCTGGAGTCGCGGGAGAGAAT

gi|302179       580 G  T  L  I  W  K  F  A  D  A  N  H  V  C  L  L  C  H  P  N  
gi|302179      1740 GGCACCCTGATCTGGAAGTTTGCAGATGCCAACCACGTGTGTCTCCTGTGCCACCCCAAC

gi|302179       600 C  T  Y  G  C  E  G  P  G  L  E  G  C  P  Q  K  G  P  K  I  
gi|302179      1800 TGCACCTATGGCTGTGAAGGGCCAGGTCTCGAAGGCTGTCCACAAAAAGGGCCCAAGATC

gi|302179       620 P  S  I  A  T  G  I  V  G  G  L  L  L  V  V  V  L  A  L  S  
gi|302179      1860 CCGTCCATTGCCACGGGCATCGTGGGCGGCCTGCTGCTGGTGGTGGTGCTGGCCCTGAGC

gi|302179       640 V  G  L  F  M  R  R  R  H  I  V  R  K  R  T  L  R  R  L  L  
gi|302179      1920 GTCGGCCTCTTCATGCGCAGGCGCCACATCGTGCGCAAGCGCACACTGCGCCGGCTGCTG

gi|302179       660 Q  E  R  E  L  V  E  P  L  T  P  S  G  E  A  P  N  Q  A  L  
gi|302179      1980 CAGGAGCGTGAGCTCGTGGAGCCTCTGACGCCCAGCGGAGAAGCTCCCAACCAAGCTCTC

gi|302179       680 L  R  I  L  K  E  T  E  F  K  K  V  K  V  L  G  S  G  A  F  
gi|302179      2040 TTGAGGATCCTAAAGGAAACAGAATTCAAGAAGGTCAAGGTGCTGGGCTCGGGAGCATTT

gi|302179       700 G  T  V  Y  K  G  L  W  I  P  E  G  E  K  V  K  I  P  V  A  
gi|302179      2100 GGCACCGTGTACAAGGGACTCTGGATCCCAGAAGGCGAGAAGGTTAAAATTCCTGTAGCT

gi|302179       720 I  K  E  L  R  E  A  T  S  P  K  A  N  K  E  I  L  D  E  A  
gi|302179      2160 ATCAAGGAATTAAGAGAAGCCACATCTCCAAAAGCCAACAAGGAAATTCTTGATGAGGCC

gi|302179       740 Y  V  M  A  S  V  D  N  P  H  V  C  R  L  L  G  I  C  L  T  
gi|302179      2220 TACGTGATGGCCAGTGTGGACAACCCCCATGTGTGCCGCCTCCTGGGCATCTGCCTGACC

gi|302179       760 S  T  V  Q  L  I  T  Q  L  M  P  F  G  C  L  L  D  Y  V  R  
gi|302179      2280 TCCACCGTGCAGCTCATCACACAGCTCATGCCCTTCGGCTGCCTGCTGGACTACGTCCGC

gi|302179       780 E  H  K  D  N  V  G  S  Q  Y  L  L  N  W  C  V  Q  I  A  K  
gi|302179      2340 GAGCACAAGGACAATGTCGGCTCCCAGTACCTGCTCAACTGGTGTGTGCAGATCGCAAAG

gi|302179       800 G  M  N  Y  L  E  D  R  R  L  V  H  R  D  L  A  A  R  N  V  
gi|302179      2400 GGCATGAATTACCTGGAAGACCGGCGCTTGGTGCATAGGGACCTGGCAGCCAGGAACGTG

gi|302179       820 L  V  K  T  P  Q  H  V  K  I  T  D  F  G  L  A  K  L  L  G  
gi|302179      2460 CTGGTGAAGACGCCGCAGCACGTGAAGATCACAGACTTCGGGCTGGCCAAGCTGCTGGGT

gi|302179       840 A  E  E  K  E  Y  H  A  E  G  G  K  V  P  I  K  W  M  A  L  
gi|302179      2520 GCCGAGGAGAAGGAGTATCATGCAGAAGGAGGCAAGGTCCCTATCAAATGGATGGCTTTG

gi|302179       860 E  S  I  L  H  R  I  Y  T  H  Q  S  D  V  W  S  Y  G  V  T  
gi|302179      2580 GAATCAATTTTACACCGAATTTATACCCATCAGAGTGATGTCTGGAGCTATGGAGTCACT

gi|302179       880 V  W  E  L  M  T  F  G  S  K  P  Y  D  G  I  P  A  S  E  I  
gi|302179      2640 GTTTGGGAGTTGATGACCTTTGGATCCAAGCCTTACGATGGAATCCCTGCGAGTGAGATC

gi|302179       900 S  T  V  L  E  K  G  E  R  L  P  Q  P  P  I  C  T  I  D  V  
gi|302179      2700 TCGACTGTCCTGGAGAAAGGAGAGCGCCTCCCACAGCCACCCATCTGCACCATCGACGTC

gi|302179       920 Y  M  I  M  V  K  C  W  M  I  D  A  D  S  R  P  K  F  R  E  
gi|302179      2760 TACATGATCATGGTCAAGTGCTGGATGATAGATGCAGACAGTCGCCCAAAGTTCCGTGAG

gi|302179       940 L  I  L  E  F  S  K  M  A  R  D  P  Q  R  Y  L  V  I  Q  G  
gi|302179      2820 TTGATCCTTGAATTCTCCAAGATGGCCCGAGACCCCCAGCGCTACCTTGTCATCCAGGGG

gi|302179       960 D  E  R  M  H  L  P  S  P  T  D  S  N  F  Y  R  A  L  M  D  
gi|302179      2880 GACGAGAGAATGCATTTGCCAAGCCCTACGGACTCCAACTTCTACCGCGCCCTGATGGAT

gi|302179       980 E  E  D  M  E  D  V  V  D  A  D  E  Y  L  V  P  Q  Q  G  F  
gi|302179      2940 GAGGAGGACATGGAGGATGTTGTGGATGCCGATGAGTACCTCGTCCCCCAGCAGGGCTTC

gi|302179      1000 F  H  S  P  T  T  S  R  T  P  L  L  S  S  L  S  T  S  S  N  
gi|302179      3000 TTCCACAGCCCCACCACCTCCCGGACACCCCTCCTCAGCTCGCTGAGCACCTCCAGCAAC

gi|302179      1020 T  P  T  V  T  C  V  D  R  N  G  S  Y  P  L  K  E  D  S  F  
gi|302179      3060 ACTCCCACTGTGACTTGCGTTGATAGAAATGGGAGCTACCCTCTCAAGGAAGACAGCTTC

gi|302179      1040 L  Q  R  Y  S  S  D  P  T  G  A  L  I  E  D  S  M  D  D  A  
gi|302179      3120 CTGCAGCGCTACAGCTCAGACCCCACTGGTGCCCTCATCGAGGACAGCATGGACGACGCT

gi|302179      1060 F  L  P  V  P  E  Y  V  N  Q  S  V  P  K  R  P  A  G  S  V  
gi|302179      3180 TTCCTCCCAGTACCCGAATATGTAAACCAATCTGTTCCCAAAAGACCCGCAGGCTCTGTC

gi|302179      1080 Q  N  P  V  Y  H  N  Q  P  L  Y  P  A  P  G  R  D  P  Q  Y  
gi|302179      3240 CAGAACCCTGTCTATCACAATCAGCCTCTATATCCAGCTCCTGGCAGAGACCCTCAGTAC

gi|302179      1100 Q  N  S  L  S  N  A  V  D  N  P  E  Y  L  N  T  T  H  P  A  
gi|302179      3300 CAAAATTCACTCAGCAACGCCGTGGACAACCCTGAGTATCTCAACACCACCCATCCTGCC

gi|302179      1120 C  I  N  G  V  L  D  G  P  A  L  W  A  Q  K  G  S  H  Q  F  
gi|302179      3360 TGTATCAATGGTGTGCTCGACGGCCCTGCCCTCTGGGCTCAGAAGGGCAGTCACCAATTT

gi|302179      1140 S  L  D  N  P  D  Y  Q  Q  A  F  F  P  K  E  A  K  S  N  G  
gi|302179      3420 AGCCTAGACAACCCTGACTACCAGCAGGCCTTCTTTCCCAAGGAAGCCAAGTCGAATGGC

gi|302179      1160 I  F  K  G  P  A  A  E  N  A  E  Y  L  R  A  A  P  A  G  S  
gi|302179      3480 ATCTTTAAGGGGCCTGCAGCTGAAAATGCAGAATACCTGCGGGCAGCACCAGCAGGCAGT

gi|302179      1180 D  F  T  G  A   1185
gi|302179      3540 GACTTTACTGGGGCC 3555
""",
        )
        protein_record = protein_alignment.sequences[3]
        self.assertEqual(protein_record.id, "gi|47522840|ref|NP_999172.1|")
        nucleotide_record = nucleotide_records["gi|47522839|ref|NM_214007.1|"]
        alignments = aligner.align(protein_record, nucleotide_record)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 1209], [126, 3753]]))
        )
        self.assertEqual(
            str(alignment),
            """\
gi|475228         0 M  R  R  S  W  A  G  G  A  A  L  L  A  L  L  A  A  H  F  Q  
gi|475228       126 ATGCGACGCTCCTGGGCGGGCGGCGCCGCGCTCCTGGCGCTGCTGGCCGCGCACTTCCAG

gi|475228        20 A  S  P  A  L  E  E  K  K  V  C  Q  G  T  S  N  K  L  T  Q  
gi|475228       186 GCGAGTCCGGCGCTGGAGGAGAAGAAAGTTTGCCAAGGTACAAGTAACAAGCTCACCCAG

gi|475228        40 L  G  T  F  E  D  H  F  L  S  L  Q  R  M  F  N  N  C  E  V  
gi|475228       246 CTGGGCACTTTCGAAGACCACTTTCTGAGCCTCCAGAGGATGTTCAATAACTGCGAGGTG

gi|475228        60 V  L  G  N  L  E  I  T  Y  M  Q  N  S  Y  N  L  S  F  L  K  
gi|475228       306 GTCCTTGGGAACTTGGAGATCACCTACATGCAGAACAGCTACAACCTGTCTTTCCTAAAG

gi|475228        80 T  I  Q  E  V  A  G  Y  V  L  I  A  L  N  T  V  E  K  I  P  
gi|475228       366 ACCATTCAGGAGGTCGCCGGCTACGTGCTCATCGCCCTCAACACCGTGGAGAAGATCCCT

gi|475228       100 L  E  N  L  Q  I  I  R  G  N  V  L  Y  E  N  T  H  A  L  A  
gi|475228       426 TTGGAAAACCTGCAGATCATCCGAGGAAATGTACTGTATGAAAACACCCATGCCTTAGCC

gi|475228       120 V  L  S  N  Y  G  A  N  K  T  G  L  R  E  L  P  M  R  N  L  
gi|475228       486 GTCTTATCCAACTACGGGGCCAATAAAACCGGCCTGAGGGAGCTGCCCATGAGGAACTTA

gi|475228       140 Q  E  I  L  Q  G  A  V  R  F  S  N  N  P  A  L  C  H  A  E  
gi|475228       546 CAAGAGATCCTGCAAGGCGCCGTGCGCTTCAGCAACAACCCTGCCCTCTGTCACGCGGAG

gi|475228       160 S  I  Q  W  R  D  I  V  N  S  D  F  L  S  N  M  S  M  D  F  
gi|475228       606 TCCATCCAGTGGAGGGACATTGTCAACAGCGACTTTCTAAGCAACATGTCCATGGACTTT

gi|475228       180 Q  S  Q  L  G  S  C  P  K  C  D  P  G  C  L  N  G  S  C  W  
gi|475228       666 CAGAGCCAGCTGGGCAGCTGCCCGAAGTGTGATCCAGGCTGTCTCAATGGGAGCTGCTGG

gi|475228       200 G  A  G  K  E  N  C  Q  K  L  T  K  V  I  C  A  Q  Q  C  S  
gi|475228       726 GGTGCTGGGAAGGAGAACTGCCAGAAATTGACCAAAGTCATCTGTGCCCAGCAGTGCTCC

gi|475228       220 G  R  C  R  G  R  S  P  S  D  C  C  H  N  Q  C  A  A  G  C  
gi|475228       786 GGGCGCTGCCGCGGCCGGTCGCCCAGTGACTGCTGCCACAACCAGTGCGCCGCTGGCTGC

gi|475228       240 T  G  P  R  E  S  D  C  L  V  C  R  R  F  R  D  E  A  T  C  
gi|475228       846 ACGGGGCCGCGGGAGAGCGACTGCCTGGTTTGCCGCAGATTCCGTGACGAGGCCACCTGC

gi|475228       260 K  D  T  C  P  P  L  M  L  Y  N  P  T  T  Y  Q  M  D  V  N  
gi|475228       906 AAGGACACATGCCCGCCGCTCATGCTCTACAACCCCACCACCTACCAGATGGACGTCAAC

gi|475228       280 P  L  G  K  Y  S  F  G  A  T  C  V  K  K  C  P  R  N  Y  V  
gi|475228       966 CCGCTGGGGAAGTACAGCTTTGGCGCCACCTGTGTCAAGAAGTGCCCTCGTAACTACGTG

gi|475228       300 V  T  D  H  G  S  C  V  R  A  C  S  S  D  S  Y  E  V  E  E  
gi|475228      1026 GTGACAGACCATGGCTCCTGTGTCCGTGCCTGCAGCTCCGACAGCTACGAGGTGGAGGAG

gi|475228       320 D  G  V  R  K  C  K  K  C  D  G  P  C  G  K  V  C  N  G  I  
gi|475228      1086 GACGGCGTCCGCAAGTGTAAAAAGTGTGACGGGCCCTGCGGCAAAGTTTGTAACGGGATA

gi|475228       340 G  I  G  E  F  K  D  T  L  S  I  N  A  T  N  I  K  H  F  R  
gi|475228      1146 GGGATTGGCGAGTTTAAAGACACACTTTCCATAAATGCTACGAATATCAAGCACTTCAGG

gi|475228       360 N  C  T  S  I  S  G  D  L  H  I  L  P  V  A  F  R  G  D  S  
gi|475228      1206 AACTGCACCTCGATCAGCGGAGATCTTCATATCCTGCCGGTAGCATTTAGGGGTGACTCC

gi|475228       380 F  T  R  T  P  P  L  D  P  K  E  L  D  I  L  K  T  V  K  E  
gi|475228      1266 TTCACACGCACGCCGCCTCTGGACCCCAAGGAACTGGACATCCTGAAAACCGTGAAGGAA

gi|475228       400 I  T  G  F  L  L  I  Q  A  W  P  E  N  R  T  G  L  H  A  F  
gi|475228      1326 ATAACAGGGTTTTTACTGATTCAGGCCTGGCCTGAAAACAGGACTGGCCTCCATGCTTTT

gi|475228       420 E  N  L  E  I  I  R  G  R  T  K  Q  H  G  Q  F  S  L  A  V  
gi|475228      1386 GAGAACCTGGAAATCATACGTGGCAGGACGAAGCAACATGGTCAGTTTTCCCTCGCGGTT

gi|475228       440 V  G  L  D  I  A  S  L  G  L  R  S  L  K  E  I  S  D  G  D  
gi|475228      1446 GTTGGCCTGGACATAGCGTCCTTGGGGATGCGCTCCCTCAAGGAGATCAGCGACGGAGAC

gi|475228       460 V  I  V  S  G  N  R  N  L  C  Y  A  N  T  I  S  W  K  K  L  
gi|475228      1506 GTGATCGTCTCAGGAAACCGAAACCTGTGCTATGCAAATACAATCAGCTGGAAAAAACTA

gi|475228       480 F  G  T  A  S  Q  K  T  K  I  I  N  N  R  S  E  K  E  C  K  
gi|475228      1566 TTTGGGACCGCAAGTCAGAAAACCAAAATTATAAACAACAGGAGCGAAAAAGAGTGCAAA

gi|475228       500 A  M  G  H  I  C  N  P  L  C  S  S  E  G  C  W  G  P  E  P  
gi|475228      1626 GCCATGGGCCACATCTGTAACCCGCTGTGCTCATCAGAGGGCTGCTGGGGCCCTGAACCC

gi|475228       520 R  D  C  M  S  C  R  N  F  S  R  G  K  E  C  V  E  K  C  N  
gi|475228      1686 AGAGACTGCATGTCCTGTCGAAACTTTAGCCGCGGCAAGGAATGTGTGGAGAAGTGCAAC

gi|475228       540 V  L  E  G  E  P  R  E  F  V  E  N  A  E  C  V  Q  C  H  P  
gi|475228      1746 GTTCTGGAGGGGGAGCCGAGAGAGTTCGTGGAGAATGCCGAGTGTGTGCAGTGCCACCCG

gi|475228       560 E  C  L  P  Q  A  K  N  V  T  C  M  G  R  G  P  D  S  C  V  
gi|475228      1806 GAGTGCCTGCCCCAGGCCAAGAACGTGACCTGCATGGGACGCGGACCGGACAGCTGTGTC

gi|475228       580 R  C  A  H  Y  I  D  G  P  H  C  V  K  T  C  P  A  G  I  A  
gi|475228      1866 CGGTGTGCTCACTACATCGACGGCCCTCACTGTGTCAAGACCTGCCCCGCGGGAATCGCA

gi|475228       600 G  E  N  S  T  L  I  W  K  F  A  D  A  N  H  V  C  H  L  C  
gi|475228      1926 GGAGAAAACAGCACCCTCATCTGGAAGTTTGCGGATGCCAACCACGTGTGTCACCTGTGC

gi|475228       620 H  P  N  C  T  Y  G  C  V  G  P  G  L  E  G  C  A  V  D  R  
gi|475228      1986 CACCCCAACTGCACCTACGGCTGTGTCGGACCAGGTCTCGAGGGCTGTGCGGTGGACAGG

gi|475228       640 P  K  I  P  S  I  A  T  G  I  V  G  G  L  L  L  A  V  V  L  
gi|475228      2046 CCCAAGATCCCGTCCATCGCCACCGGGATAGTGGGGGGCCTGCTTCTGGCCGTGGTGCTG

gi|475228       660 A  L  G  V  G  L  F  L  R  R  R  H  I  V  R  K  R  T  L  R  
gi|475228      2106 GCCCTGGGGGTCGGCCTCTTTCTGCGCAGGCGCCACATCGTCCGCAAGCGCACGCTGCGC

gi|475228       680 R  L  L  Q  E  R  E  L  V  E  P  L  T  P  S  G  E  A  P  N  
gi|475228      2166 CGGCTGCTGCAGGAGCGGGAGCTGGTTGAGCCTCTCACACCCAGTGGAGAAGCTCCCAAC

gi|475228       700 Q  A  L  L  R  I  L  K  E  T  E  F  K  K  V  K  V  L  G  S  
gi|475228      2226 CAAGCTCTCTTGAGGATCCTGAAGGAGACGGAATTCAAAAAGGTCAAGGTGCTGGGCTCC

gi|475228       720 G  A  F  G  T  V  Y  K  G  L  W  I  P  E  G  E  K  V  K  I  
gi|475228      2286 GGCGCGTTCGGCACGGTGTACAAGGGCCTCTGGATCCCAGAAGGTGAGAAGGTGAAAATT

gi|475228       740 P  V  A  I  K  E  L  R  E  A  T  S  P  K  A  N  K  E  I  L  
gi|475228      2346 CCTGTGGCTATCAAGGAATTAAGAGAAGCCACTTCTCCAAAAGCCAACAAGGAAATTCTT

gi|475228       760 D  E  A  Y  V  M  A  S  V  D  N  P  H  V  C  R  L  L  G  I  
gi|475228      2406 GACGAAGCCTACGTGATGGCCAGTGTGGACAATCCTCATGTGTGCCGCCTCCTGGGCATC

gi|475228       780 C  L  T  S  T  V  Q  L  I  T  Q  L  M  P  F  G  C  L  L  D  
gi|475228      2466 TGCCTGACCTCCACGGTGCAGCTCATCACGCAGCTCATGCCCTTCGGCTGCCTCCTGGAC

gi|475228       800 Y  V  R  E  H  K  D  N  I  G  S  Q  H  L  L  N  W  C  V  Q  
gi|475228      2526 TACGTCCGCGAGCACAAGGACAACATCGGCTCCCAGCACCTGCTCAACTGGTGTGTGCAG

gi|475228       820 I  A  K  G  M  N  Y  L  E  D  R  R  L  V  H  R  D  L  A  A  
gi|475228      2586 ATCGCAAAGGGCATGAACTATCTGGAAGACCGGCGCTTGGTGCACCGAGACCTGGCGGCC

gi|475228       840 R  N  V  L  V  K  T  P  Q  H  V  K  I  T  D  F  G  L  A  K  
gi|475228      2646 AGGAATGTGCTGGTGAAGACACCGCAGCATGTCAAGATCACTGACTTTGGGCTGGCCAAG

gi|475228       860 L  L  G  A  E  E  K  E  Y  H  A  E  G  G  K  V  P  I  K  W  
gi|475228      2706 CTGCTGGGCGCCGAGGAGAAAGAGTACCACGCGGAAGGAGGCAAAGTGCCCATCAAGTGG

gi|475228       880 L  A  L  E  S  I  L  H  R  V  Y  T  H  Q  S  D  V  W  S  Y  
gi|475228      2766 CTGGCTCTAGAGTCAATCCTGCACCGTGTATACACCCACCAGAGTGACGTCTGGAGCTAC

gi|475228       900 G  V  T  V  W  E  L  M  T  F  G  S  K  P  Y  D  G  I  P  A  
gi|475228      2826 GGAGTCACCGTTTGGGAGCTGATGACCTTTGGGTCCAAGCCTTATGACGGGATCCCCGCG

gi|475228       920 S  E  I  S  T  V  L  E  K  G  E  R  L  P  Q  P  P  I  C  T  
gi|475228      2886 AGTGAGATCTCGACCGTCCTGGAGAAGGGAGAGCGCCTCCCGCAGCCCCCCATCTGCACC

gi|475228       940 I  D  V  Y  M  I  M  V  K  C  W  M  I  D  A  D  S  R  P  K  
gi|475228      2946 ATTGATGTCTACATGATCATGGTCAAGTGCTGGATGATAGATGCTGATAGTCGCCCAAAG

gi|475228       960 F  R  E  L  I  I  E  F  S  K  M  A  R  D  P  Q  R  Y  L  V  
gi|475228      3006 TTCCGTGAGCTGATCATCGAATTCTCCAAAATGGCCCGAGACCCCCAGCGCTACCTTGTC

gi|475228       980 I  Q  G  D  E  R  M  H  L  P  S  P  T  D  S  N  F  Y  R  A  
gi|475228      3066 ATCCAGGGAGACGAGCGAATGCACTTGCCAAGCCCTACGGACTCCAACTTCTACCGCGCC

gi|475228      1000 L  M  D  E  E  D  M  E  D  V  V  D  A  D  E  Y  L  V  P  Q  
gi|475228      3126 CTGATGGACGAGGAGGACATGGAGGATGTGGTGGACGCCGACGAGTACCTCGTCCCCCAG

gi|475228      1020 Q  G  F  F  H  S  P  A  T  S  R  T  P  L  L  S  S  L  S  A  
gi|475228      3186 CAGGGCTTCTTCCACAGCCCCGCCACCTCCCGGACGCCGCTGCTCAGCTCTCTGAGCGCC

gi|475228      1040 T  S  S  T  P  A  V  A  C  V  D  R  N  G  Q  S  Y  P  L  K  
gi|475228      3246 ACCAGCAGCACCCCCGCTGTGGCTTGCGTTGACAGAAACGGGCAGAGTTATCCCCTCAAG

gi|475228      1060 E  D  S  F  L  Q  R  Y  S  S  D  P  T  G  A  L  T  E  D  S  
gi|475228      3306 GAAGACAGCTTCCTGCAGCGGTACAGCTCCGACCCCACTGGCGCCCTGACCGAGGACAGC

gi|475228      1080 L  D  D  T  F  L  P  A  P  E  Y  V  N  Q  S  V  P  K  R  P  
gi|475228      3366 CTAGACGACACTTTTCTCCCAGCACCCGAATATGTAAACCAGTCTGTTCCCAAGAGGCCC

gi|475228      1100 A  G  S  V  Q  N  P  V  Y  H  N  Q  P  L  S  A  A  P  G  R  
gi|475228      3426 GCGGGCTCCGTCCAGAACCCTGTCTACCACAATCAGCCTCTCAGTGCAGCTCCTGGCCGG

gi|475228      1120 D  P  H  Y  Q  N  S  H  S  N  A  V  G  N  P  E  Y  L  N  T  
gi|475228      3486 GACCCCCACTACCAGAACTCCCACAGCAATGCCGTGGGCAACCCTGAGTATCTCAACACC

gi|475228      1140 P  R  P  A  C  I  N  G  G  L  D  G  P  A  F  W  A  Q  T  G  
gi|475228      3546 CCCCGCCCCGCCTGCATCAACGGAGGACTGGACGGCCCTGCCTTCTGGGCACAGACAGGC

gi|475228      1160 S  H  Q  I  N  L  D  N  P  D  Y  Q  Q  A  F  F  P  K  E  A  
gi|475228      3606 AGCCACCAGATTAATCTGGACAACCCAGACTACCAGCAGGCCTTCTTCCCCAAGGAAGCC

gi|475228      1180 K  S  N  G  I  C  K  G  P  A  A  E  N  A  E  Y  L  R  A  A  
gi|475228      3666 AAGTCAAACGGCATCTGTAAGGGTCCCGCCGCCGAAAACGCAGAGTACCTAAGGGCGGCA

gi|475228      1200 P  A  S  S  D  L  T  G  A   1209
gi|475228      3726 CCAGCCAGCAGTGACCTTACTGGGGCA 3753
""",
        )
        protein_record = protein_alignment.sequences[4]
        self.assertEqual(protein_record.id, "gi|29725609|ref|NP_005219.2|")
        nucleotide_record = nucleotide_records["gi|41327737|ref|NM_005228.3|"]
        alignments = aligner.align(protein_record, nucleotide_record)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 1210], [246, 3876]]))
        )
        self.assertEqual(
            str(alignment),
            """\
gi|297256         0 M  R  P  S  G  T  A  G  A  A  L  L  A  L  L  A  A  L  C  P  
gi|413277       246 ATGCGACCCTCCGGGACGGCCGGGGCAGCGCTCCTGGCGCTGCTGGCTGCGCTCTGCCCG

gi|297256        20 A  S  R  A  L  E  E  K  K  V  C  Q  G  T  S  N  K  L  T  Q  
gi|413277       306 GCGAGTCGGGCTCTGGAGGAAAAGAAAGTTTGCCAAGGCACGAGTAACAAGCTCACGCAG

gi|297256        40 L  G  T  F  E  D  H  F  L  S  L  Q  R  M  F  N  N  C  E  V  
gi|413277       366 TTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGATGTTCAATAACTGTGAGGTG

gi|297256        60 V  L  G  N  L  E  I  T  Y  V  Q  R  N  Y  D  L  S  F  L  K  
gi|413277       426 GTCCTTGGGAATTTGGAAATTACCTATGTGCAGAGGAATTATGATCTTTCCTTCTTAAAG

gi|297256        80 T  I  Q  E  V  A  G  Y  V  L  I  A  L  N  T  V  E  R  I  P  
gi|413277       486 ACCATCCAGGAGGTGGCTGGTTATGTCCTCATTGCCCTCAACACAGTGGAGCGAATTCCT

gi|297256       100 L  E  N  L  Q  I  I  R  G  N  M  Y  Y  E  N  S  Y  A  L  A  
gi|413277       546 TTGGAAAACCTGCAGATCATCAGAGGAAATATGTACTACGAAAATTCCTATGCCTTAGCA

gi|297256       120 V  L  S  N  Y  D  A  N  K  T  G  L  K  E  L  P  M  R  N  L  
gi|413277       606 GTCTTATCTAACTATGATGCAAATAAAACCGGACTGAAGGAGCTGCCCATGAGAAATTTA

gi|297256       140 Q  E  I  L  H  G  A  V  R  F  S  N  N  P  A  L  C  N  V  E  
gi|413277       666 CAGGAAATCCTGCATGGCGCCGTGCGGTTCAGCAACAACCCTGCCCTGTGCAACGTGGAG

gi|297256       160 S  I  Q  W  R  D  I  V  S  S  D  F  L  S  N  M  S  M  D  F  
gi|413277       726 AGCATCCAGTGGCGGGACATAGTCAGCAGTGACTTTCTCAGCAACATGTCGATGGACTTC

gi|297256       180 Q  N  H  L  G  S  C  Q  K  C  D  P  S  C  P  N  G  S  C  W  
gi|413277       786 CAGAACCACCTGGGCAGCTGCCAAAAGTGTGATCCAAGCTGTCCCAATGGGAGCTGCTGG

gi|297256       200 G  A  G  E  E  N  C  Q  K  L  T  K  I  I  C  A  Q  Q  C  S  
gi|413277       846 GGTGCAGGAGAGGAGAACTGCCAGAAACTGACCAAAATCATCTGTGCCCAGCAGTGCTCC

gi|297256       220 G  R  C  R  G  K  S  P  S  D  C  C  H  N  Q  C  A  A  G  C  
gi|413277       906 GGGCGCTGCCGTGGCAAGTCCCCCAGTGACTGCTGCCACAACCAGTGTGCTGCAGGCTGC

gi|297256       240 T  G  P  R  E  S  D  C  L  V  C  R  K  F  R  D  E  A  T  C  
gi|413277       966 ACAGGCCCCCGGGAGAGCGACTGCCTGGTCTGCCGCAAATTCCGAGACGAAGCCACGTGC

gi|297256       260 K  D  T  C  P  P  L  M  L  Y  N  P  T  T  Y  Q  M  D  V  N  
gi|413277      1026 AAGGACACCTGCCCCCCACTCATGCTCTACAACCCCACCACGTACCAGATGGATGTGAAC

gi|297256       280 P  E  G  K  Y  S  F  G  A  T  C  V  K  K  C  P  R  N  Y  V  
gi|413277      1086 CCCGAGGGCAAATACAGCTTTGGTGCCACCTGCGTGAAGAAGTGTCCCCGTAATTATGTG

gi|297256       300 V  T  D  H  G  S  C  V  R  A  C  G  A  D  S  Y  E  M  E  E  
gi|413277      1146 GTGACAGATCACGGCTCGTGCGTCCGAGCCTGTGGGGCCGACAGCTATGAGATGGAGGAA

gi|297256       320 D  G  V  R  K  C  K  K  C  E  G  P  C  R  K  V  C  N  G  I  
gi|413277      1206 GACGGCGTCCGCAAGTGTAAGAAGTGCGAAGGGCCTTGCCGCAAAGTGTGTAACGGAATA

gi|297256       340 G  I  G  E  F  K  D  S  L  S  I  N  A  T  N  I  K  H  F  K  
gi|413277      1266 GGTATTGGTGAATTTAAAGACTCACTCTCCATAAATGCTACGAATATTAAACACTTCAAA

gi|297256       360 N  C  T  S  I  S  G  D  L  H  I  L  P  V  A  F  R  G  D  S  
gi|413277      1326 AACTGCACCTCCATCAGTGGCGATCTCCACATCCTGCCGGTGGCATTTAGGGGTGACTCC

gi|297256       380 F  T  H  T  P  P  L  D  P  Q  E  L  D  I  L  K  T  V  K  E  
gi|413277      1386 TTCACACATACTCCTCCTCTGGATCCACAGGAACTGGATATTCTGAAAACCGTAAAGGAA

gi|297256       400 I  T  G  F  L  L  I  Q  A  W  P  E  N  R  T  D  L  H  A  F  
gi|413277      1446 ATCACAGGGTTTTTGCTGATTCAGGCTTGGCCTGAAAACAGGACGGACCTCCATGCCTTT

gi|297256       420 E  N  L  E  I  I  R  G  R  T  K  Q  H  G  Q  F  S  L  A  V  
gi|413277      1506 GAGAACCTAGAAATCATACGCGGCAGGACCAAGCAACATGGTCAGTTTTCTCTTGCAGTC

gi|297256       440 V  S  L  N  I  T  S  L  G  L  R  S  L  K  E  I  S  D  G  D  
gi|413277      1566 GTCAGCCTGAACATAACATCCTTGGGATTACGCTCCCTCAAGGAGATAAGTGATGGAGAT

gi|297256       460 V  I  I  S  G  N  K  N  L  C  Y  A  N  T  I  N  W  K  K  L  
gi|413277      1626 GTGATAATTTCAGGAAACAAAAATTTGTGCTATGCAAATACAATAAACTGGAAAAAACTG

gi|297256       480 F  G  T  S  G  Q  K  T  K  I  I  S  N  R  G  E  N  S  C  K  
gi|413277      1686 TTTGGGACCTCCGGTCAGAAAACCAAAATTATAAGCAACAGAGGTGAAAACAGCTGCAAG

gi|297256       500 A  T  G  Q  V  C  H  A  L  C  S  P  E  G  C  W  G  P  E  P  
gi|413277      1746 GCCACAGGCCAGGTCTGCCATGCCTTGTGCTCCCCCGAGGGCTGCTGGGGCCCGGAGCCC

gi|297256       520 R  D  C  V  S  C  R  N  V  S  R  G  R  E  C  V  D  K  C  N  
gi|413277      1806 AGGGACTGCGTCTCTTGCCGGAATGTCAGCCGAGGCAGGGAATGCGTGGACAAGTGCAAC

gi|297256       540 L  L  E  G  E  P  R  E  F  V  E  N  S  E  C  I  Q  C  H  P  
gi|413277      1866 CTTCTGGAGGGTGAGCCAAGGGAGTTTGTGGAGAACTCTGAGTGCATACAGTGCCACCCA

gi|297256       560 E  C  L  P  Q  A  M  N  I  T  C  T  G  R  G  P  D  N  C  I  
gi|413277      1926 GAGTGCCTGCCTCAGGCCATGAACATCACCTGCACAGGACGGGGACCAGACAACTGTATC

gi|297256       580 Q  C  A  H  Y  I  D  G  P  H  C  V  K  T  C  P  A  G  V  M  
gi|413277      1986 CAGTGTGCCCACTACATTGACGGCCCCCACTGCGTCAAGACCTGCCCGGCAGGAGTCATG

gi|297256       600 G  E  N  N  T  L  V  W  K  Y  A  D  A  G  H  V  C  H  L  C  
gi|413277      2046 GGAGAAAACAACACCCTGGTCTGGAAGTACGCAGACGCCGGCCATGTGTGCCACCTGTGC

gi|297256       620 H  P  N  C  T  Y  G  C  T  G  P  G  L  E  G  C  P  T  N  G  
gi|413277      2106 CATCCAAACTGCACCTACGGATGCACTGGGCCAGGTCTTGAAGGCTGTCCAACGAATGGG

gi|297256       640 P  K  I  P  S  I  A  T  G  M  V  G  A  L  L  L  L  L  V  V  
gi|413277      2166 CCTAAGATCCCGTCCATCGCCACTGGGATGGTGGGGGCCCTCCTCTTGCTGCTGGTGGTG

gi|297256       660 A  L  G  I  G  L  F  M  R  R  R  H  I  V  R  K  R  T  L  R  
gi|413277      2226 GCCCTGGGGATCGGCCTCTTCATGCGAAGGCGCCACATCGTTCGGAAGCGCACGCTGCGG

gi|297256       680 R  L  L  Q  E  R  E  L  V  E  P  L  T  P  S  G  E  A  P  N  
gi|413277      2286 AGGCTGCTGCAGGAGAGGGAGCTTGTGGAGCCTCTTACACCCAGTGGAGAAGCTCCCAAC

gi|297256       700 Q  A  L  L  R  I  L  K  E  T  E  F  K  K  I  K  V  L  G  S  
gi|413277      2346 CAAGCTCTCTTGAGGATCTTGAAGGAAACTGAATTCAAAAAGATCAAAGTGCTGGGCTCC

gi|297256       720 G  A  F  G  T  V  Y  K  G  L  W  I  P  E  G  E  K  V  K  I  
gi|413277      2406 GGTGCGTTCGGCACGGTGTATAAGGGACTCTGGATCCCAGAAGGTGAGAAAGTTAAAATT

gi|297256       740 P  V  A  I  K  E  L  R  E  A  T  S  P  K  A  N  K  E  I  L  
gi|413277      2466 CCCGTCGCTATCAAGGAATTAAGAGAAGCAACATCTCCGAAAGCCAACAAGGAAATCCTC

gi|297256       760 D  E  A  Y  V  M  A  S  V  D  N  P  H  V  C  R  L  L  G  I  
gi|413277      2526 GATGAAGCCTACGTGATGGCCAGCGTGGACAACCCCCACGTGTGCCGCCTGCTGGGCATC

gi|297256       780 C  L  T  S  T  V  Q  L  I  T  Q  L  M  P  F  G  C  L  L  D  
gi|413277      2586 TGCCTCACCTCCACCGTGCAGCTCATCACGCAGCTCATGCCCTTCGGCTGCCTCCTGGAC

gi|297256       800 Y  V  R  E  H  K  D  N  I  G  S  Q  Y  L  L  N  W  C  V  Q  
gi|413277      2646 TATGTCCGGGAACACAAAGACAATATTGGCTCCCAGTACCTGCTCAACTGGTGTGTGCAG

gi|297256       820 I  A  K  G  M  N  Y  L  E  D  R  R  L  V  H  R  D  L  A  A  
gi|413277      2706 ATCGCAAAGGGCATGAACTACTTGGAGGACCGTCGCTTGGTGCACCGCGACCTGGCAGCC

gi|297256       840 R  N  V  L  V  K  T  P  Q  H  V  K  I  T  D  F  G  L  A  K  
gi|413277      2766 AGGAACGTACTGGTGAAAACACCGCAGCATGTCAAGATCACAGATTTTGGGCTGGCCAAA

gi|297256       860 L  L  G  A  E  E  K  E  Y  H  A  E  G  G  K  V  P  I  K  W  
gi|413277      2826 CTGCTGGGTGCGGAAGAGAAAGAATACCATGCAGAAGGAGGCAAAGTGCCTATCAAGTGG

gi|297256       880 M  A  L  E  S  I  L  H  R  I  Y  T  H  Q  S  D  V  W  S  Y  
gi|413277      2886 ATGGCATTGGAATCAATTTTACACAGAATCTATACCCACCAGAGTGATGTCTGGAGCTAC

gi|297256       900 G  V  T  V  W  E  L  M  T  F  G  S  K  P  Y  D  G  I  P  A  
gi|413277      2946 GGGGTGACCGTTTGGGAGTTGATGACCTTTGGATCCAAGCCATATGACGGAATCCCTGCC

gi|297256       920 S  E  I  S  S  I  L  E  K  G  E  R  L  P  Q  P  P  I  C  T  
gi|413277      3006 AGCGAGATCTCCTCCATCCTGGAGAAAGGAGAACGCCTCCCTCAGCCACCCATATGTACC

gi|297256       940 I  D  V  Y  M  I  M  V  K  C  W  M  I  D  A  D  S  R  P  K  
gi|413277      3066 ATCGATGTCTACATGATCATGGTCAAGTGCTGGATGATAGACGCAGATAGTCGCCCAAAG

gi|297256       960 F  R  E  L  I  I  E  F  S  K  M  A  R  D  P  Q  R  Y  L  V  
gi|413277      3126 TTCCGTGAGTTGATCATCGAATTCTCCAAAATGGCCCGAGACCCCCAGCGCTACCTTGTC

gi|297256       980 I  Q  G  D  E  R  M  H  L  P  S  P  T  D  S  N  F  Y  R  A  
gi|413277      3186 ATTCAGGGGGATGAAAGAATGCATTTGCCAAGTCCTACAGACTCCAACTTCTACCGTGCC

gi|297256      1000 L  M  D  E  E  D  M  D  D  V  V  D  A  D  E  Y  L  I  P  Q  
gi|413277      3246 CTGATGGATGAAGAAGACATGGACGACGTGGTGGATGCCGACGAGTACCTCATCCCACAG

gi|297256      1020 Q  G  F  F  S  S  P  S  T  S  R  T  P  L  L  S  S  L  S  A  
gi|413277      3306 CAGGGCTTCTTCAGCAGCCCCTCCACGTCACGGACTCCCCTCCTGAGCTCTCTGAGTGCA

gi|297256      1040 T  S  N  N  S  T  V  A  C  I  D  R  N  G  L  Q  S  C  P  I  
gi|413277      3366 ACCAGCAACAATTCCACCGTGGCTTGCATTGATAGAAATGGGCTGCAAAGCTGTCCCATC

gi|297256      1060 K  E  D  S  F  L  Q  R  Y  S  S  D  P  T  G  A  L  T  E  D  
gi|413277      3426 AAGGAAGACAGCTTCTTGCAGCGATACAGCTCAGACCCCACAGGCGCCTTGACTGAGGAC

gi|297256      1080 S  I  D  D  T  F  L  P  V  P  E  Y  I  N  Q  S  V  P  K  R  
gi|413277      3486 AGCATAGACGACACCTTCCTCCCAGTGCCTGAATACATAAACCAGTCCGTTCCCAAAAGG

gi|297256      1100 P  A  G  S  V  Q  N  P  V  Y  H  N  Q  P  L  N  P  A  P  S  
gi|413277      3546 CCCGCTGGCTCTGTGCAGAATCCTGTCTATCACAATCAGCCTCTGAACCCCGCGCCCAGC

gi|297256      1120 R  D  P  H  Y  Q  D  P  H  S  T  A  V  G  N  P  E  Y  L  N  
gi|413277      3606 AGAGACCCACACTACCAGGACCCCCACAGCACTGCAGTGGGCAACCCCGAGTATCTCAAC

gi|297256      1140 T  V  Q  P  T  C  V  N  S  T  F  D  S  P  A  H  W  A  Q  K  
gi|413277      3666 ACTGTCCAGCCCACCTGTGTCAACAGCACATTCGACAGCCCTGCCCACTGGGCCCAGAAA

gi|297256      1160 G  S  H  Q  I  S  L  D  N  P  D  Y  Q  Q  D  F  F  P  K  E  
gi|413277      3726 GGCAGCCACCAAATTAGCCTGGACAACCCTGACTACCAGCAGGACTTCTTTCCCAAGGAA

gi|297256      1180 A  K  P  N  G  I  F  K  G  S  T  A  E  N  A  E  Y  L  R  V  
gi|413277      3786 GCCAAGCCAAATGGCATCTTTAAGGGCTCCACAGCTGAAAATGCAGAATACCTAAGGGTC

gi|297256      1200 A  P  Q  S  S  E  F  I  G  A   1210
gi|413277      3846 GCGCCACAAAGCAGTGAATTTATTGGAGCA 3876
""",
        )
        protein_record = protein_alignment.sequences[5]
        self.assertEqual(protein_record.id, "gi|6478868|gb|AAF14008.1|")
        nucleotide_record = nucleotide_records["gi|6478867|gb|M37394.2|RATEGFR"]
        alignments = aligner.align(protein_record, nucleotide_record)
        self.assertEqual(len(alignments), 1)
        alignment = next(alignments)
        codon_alignments.append(alignment)
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 1209], [153, 3780]]))
        )
        self.assertEqual(
            str(alignment),
            """\
gi|647886         0 M  R  P  S  G  T  A  R  T  K  L  L  L  L  L  A  A  L  C  A  
gi|647886       153 ATGCGACCCTCAGGGACTGCGAGAACCAAGCTACTGCTGCTGCTGGCTGCGCTCTGCGCC

gi|647886        20 A  G  G  A  L  E  E  K  K  V  C  Q  G  T  S  N  R  L  T  Q  
gi|647886       213 GCAGGTGGGGCGCTGGAGGAAAAGAAAGTTTGCCAAGGCACAAGTAACAGGCTCACCCAA

gi|647886        40 L  G  T  F  E  D  H  F  L  S  L  Q  R  M  F  N  N  C  E  V  
gi|647886       273 CTAGGCACCTTTGAAGACCACTTTCTGAGCCTCCAGAGGATGTTCAACAACTGTGAAGTG

gi|647886        60 V  L  G  N  L  E  I  T  Y  V  Q  R  N  Y  D  L  S  F  L  K  
gi|647886       333 GTCCTTGGAAACTTGGAAATCACCTATGTGCAAAGGAATTATGACCTTTCCTTCTTAAAG

gi|647886        80 T  I  Q  E  V  A  G  Y  V  L  I  A  L  N  T  V  E  R  I  P  
gi|647886       393 ACCATCCAGGAGGTGGCTGGCTATGTTCTCATTGCCCTGAACACCGTGGAGAGAATCCCT

gi|647886       100 L  E  N  L  Q  I  I  R  G  N  A  L  Y  E  N  T  Y  A  L  A  
gi|647886       453 TTGGAGAACCTGCAGATCATCAGGGGAAATGCTCTCTACGAAAACACCTACGCCTTAGCC

gi|647886       120 V  L  S  N  Y  G  T  N  K  T  G  L  R  E  L  P  M  R  N  L  
gi|647886       513 GTCCTGTCCAACTATGGAACCAACAAAACTGGGCTTAGGGAACTGCCCATGCGGAACTTA

gi|647886       140 Q  E  I  L  I  G  A  V  R  F  S  N  N  P  I  L  C  N  M  E  
gi|647886       573 CAGGAAATTCTGATCGGTGCTGTGCGATTTAGCAACAACCCCATCCTCTGCAATATGGAG

gi|647886       160 T  I  Q  W  R  D  I  V  Q  D  V  F  L  S  N  M  S  M  D  V  
gi|647886       633 ACCATCCAGTGGAGGGACATCGTCCAAGATGTCTTTCTGAGCAACATGTCAATGGACGTA

gi|647886       180 Q  R  H  L  T  G  C  P  K  C  D  P  S  C  P  N  G  S  C  W  
gi|647886       693 CAGCGCCACCTGACGGGCTGCCCGAAATGTGATCCGAGCTGTCCCAATGGAAGCTGCTGG

gi|647886       200 G  R  G  E  E  N  C  Q  K  L  T  K  I  I  C  A  Q  Q  C  S  
gi|647886       753 GGAAGAGGAGAGGAGAACTGCCAGAAATTGACCAAAATCATCTGCGCCCAGCAATGTTCC

gi|647886       220 R  R  C  R  G  R  S  P  S  D  C  C  H  N  Q  C  A  A  G  C  
gi|647886       813 CGGCGTTGTCGTGGCAGGTCCCCTAGCGACTGCTGCCACAACCAGTGTGCCGCAGGGTGT

gi|647886       240 T  G  P  R  E  S  D  C  L  V  C  H  R  F  R  D  E  A  T  C  
gi|647886       873 ACAGGGCCCAGAGAGAGTGACTGTCTGGTCTGCCACAGGTTCCGAGATGAAGCCACGTGC

gi|647886       260 K  D  T  C  P  P  L  M  L  Y  N  P  T  T  Y  Q  M  D  V  N  
gi|647886       933 AAAGACACCTGCCCACCACTCATGCTGTACAACCCCACCACGTACCAGATGGATGTCAAC

gi|647886       280 P  E  G  K  Y  S  F  G  A  T  C  V  K  K  C  P  R  N  Y  V  
gi|647886       993 CCTGAGGGGAAGTACAGCTTTGGTGCCACCTGTGTGAAGAAATGCCCCAGAAACTACGTG

gi|647886       300 V  T  D  H  G  S  C  V  R  A  C  G  P  D  Y  Y  E  V  E  E  
gi|647886      1053 GTGACAGATCACGGCTCGTGTGTCCGGGCCTGTGGGCCAGACTACTATGAAGTAGAAGAA

gi|647886       320 D  G  V  S  K  C  K  K  C  D  G  P  C  R  K  V  C  N  G  I  
gi|647886      1113 GATGGAGTCAGCAAGTGTAAAAAATGTGACGGGCCCTGCCGCAAAGTTTGCAATGGCATA

gi|647886       340 G  I  G  E  F  K  D  T  L  S  I  N  A  T  N  I  K  H  F  K  
gi|647886      1173 GGCATTGGTGAATTTAAAGACACACTCTCCATAAATGCTACAAACATCAAACACTTCAAG

gi|647886       360 Y  C  T  A  I  S  G  D  L  H  I  L  P  V  A  F  K  G  D  S  
gi|647886      1233 TACTGCACTGCCATCAGTGGGGACCTCCACATCCTGCCAGTGGCCTTTAAGGGGGATTCT

gi|647886       380 F  T  R  T  P  P  L  D  P  R  E  L  E  I  L  K  T  V  K  E  
gi|647886      1293 TTCACCCGCACTCCTCCTCTAGACCCACGGGAACTAGAAATTCTCAAAACTGTGAAGGAA

gi|647886       400 I  T  G  F  L  L  I  Q  A  W  P  E  N  W  T  D  L  H  A  F  
gi|647886      1353 ATAACAGGGTTTTTGCTGATTCAGGCTTGGCCTGAAAACTGGACTGACCTCCATGCTTTT

gi|647886       420 E  N  L  E  I  I  R  G  R  T  K  Q  H  G  Q  F  S  L  A  V  
gi|647886      1413 GAGAACCTAGAAATAATTCGTGGCAGAACAAAGCAACATGGTCAGTTTTCTCTGGCGGTT

gi|647886       440 V  G  L  N  I  T  S  L  G  L  R  S  L  K  E  I  S  D  G  D  
gi|647886      1473 GTCGGCCTGAACATAACATCGCTGGGGTTGCGTTCCCTCAAGGAGATCAGTGATGGGGAT

gi|647886       460 V  I  I  S  G  N  R  N  L  C  Y  A  N  T  I  N  W  K  K  L  
gi|647886      1533 GTGATTATTTCTGGGAACCGAAATTTGTGCTACGCAAACACTATAAACTGGAAAAAACTC

gi|647886       480 F  G  T  P  N  Q  K  T  K  I  M  N  N  R  A  E  K  D  C  K  
gi|647886      1593 TTCGGGACGCCCAATCAAAAGACCAAAATCATGAACAACAGAGCTGAAAAGGACTGCAAG

gi|647886       500 A  T  N  H  V  C  N  P  L  C  S  S  E  G  C  W  G  P  E  P  
gi|647886      1653 GCCACGAACCACGTCTGTAATCCTTTATGCTCCTCGGAAGGCTGCTGGGGCCCTGAGCCC

gi|647886       520 T  D  C  V  S  C  Q  N  V  S  R  G  R  E  C  V  D  K  C  N  
gi|647886      1713 ACGGACTGTGTCTCCTGCCAGAATGTGAGCAGAGGCAGGGAGTGCGTGGACAAGTGCAAC

gi|647886       540 I  L  E  G  E  P  R  E  F  V  E  N  S  E  C  I  Q  C  H  P  
gi|647886      1773 ATCCTGGAGGGGGAACCGAGGGAGTTTGTGGAAAATTCTGAATGCATCCAGTGCCATCCA

gi|647886       560 E  C  L  P  Q  T  M  N  I  T  C  T  G  R  G  P  D  N  C  I  
gi|647886      1833 GAATGTCTGCCCCAGACCATGAACATCACCTGTACAGGCCGGGGGCCAGACAACTGCATC

gi|647886       580 K  C  A  H  Y  V  D  G  P  H  C  V  K  T  C  P  S  G  I  M  
gi|647886      1893 AAGTGTGCCCACTATGTTGATGGTCCCCACTGTGTCAAGACCTGCCCTTCGGGCATCATG

gi|647886       600 G  E  N  N  T  L  V  W  K  F  A  D  A  N  N  V  C  H  L  C  
gi|647886      1953 GGGGAGAACAACACCCTGGTCTGGAAGTTTGCAGATGCCAATAACGTCTGCCACCTCTGC

gi|647886       620 H  A  N  C  T  Y  G  C  A  G  P  G  L  K  G  C  Q  Q  P  E  
gi|647886      2013 CATGCAAACTGTACCTATGGATGTGCTGGGCCAGGCCTTAAAGGATGTCAACAACCAGAA

gi|647886       640 G  P  K  I  P  S  I  A  T  G  I  V  G  G  L  L  F  I  V  V  
gi|647886      2073 GGGCCAAAGATCCCATCCATCGCCACTGGGATTGTGGGTGGCCTCCTCTTCATAGTAGTG

gi|647886       660 V  A  L  G  I  G  L  F  M  R  R  R  Q  L  V  R  K  R  T  L  
gi|647886      2133 GTGGCCCTTGGGATCGGCCTCTTCATGCGTCGACGTCAGCTTGTCCGAAAACGTACACTA

gi|647886       680 R  R  L  L  Q  E  R  E  L  V  E  P  L  T  P  S  G  E  A  P  
gi|647886      2193 CGCCGCCTGCTTCAAGAGAGAGAGCTCGTGGAACCTCTCACACCCAGCGGAGAAGCTCCG

gi|647886       700 N  Q  A  H  L  R  I  L  K  E  T  E  F  K  K  I  K  V  L  G  
gi|647886      2253 AACCAAGCCCACTTGAGGATATTAAAGGAAACAGAATTCAAAAAGATCAAAGTTCTGGGT

gi|647886       720 S  G  A  F  G  T  V  Y  K  G  L  W  I  P  E  G  E  K  V  K  
gi|647886      2313 TCAGGAGCATTTGGCACAGTGTATAAGGGTCTCTGGATCCCAGAAGGCGAGAAAGTGAAA

gi|647886       740 I  P  V  A  I  K  E  L  R  E  A  T  S  P  K  A  N  K  E  I  
gi|647886      2373 ATCCCTGTGGCCATCAAGGAGTTAAGAGAAGCCACATCTCCCAAAGCCAACAAGGAAATC

gi|647886       760 L  D  E  A  Y  V  M  A  S  V  D  N  P  H  V  C  R  L  L  G  
gi|647886      2433 CTTGATGAAGCCTACGTGATGGCCAGTGTGGACAACCCTCATGTATGCCGCCTCCTGGGC

gi|647886       780 I  C  L  T  S  T  V  Q  L  I  T  Q  L  M  P  Y  G  C  L  L  
gi|647886      2493 ATCTGTCTGACCTCCACTGTCCAGCTCATTACACAACTCATGCCCTATGGTTGCCTCCTG

gi|647886       800 D  Y  V  R  E  H  K  D  N  I  G  S  Q  Y  L  L  N  W  C  V  
gi|647886      2553 GACTATGTCCGAGAACATAAGGACAACATTGGCTCCCAGTACCTACTCAACTGGTGTGTG

gi|647886       820 Q  I  A  K  G  M  N  Y  L  E  D  R  R  L  V  H  R  D  L  A  
gi|647886      2613 CAGATTGCAAAGGGCATGAACTACCTGGAAGACCGGCGTTTGGTACACCGTGACTTGGCA

gi|647886       840 A  R  N  V  L  V  K  T  P  Q  H  V  K  I  T  D  F  G  L  A  
gi|647886      2673 GCCAGGAATGTACTGGTAAAGACACCACAGCATGTCAAGATCACAGATTTTGGACTGGCC

gi|647886       860 K  L  L  G  A  E  E  K  E  Y  H  A  E  G  G  K  V  P  I  K  
gi|647886      2733 AAACTGCTTGGTGCTGAGGAGAAAGAATACCATGCAGAGGGGGGCAAAGTGCCTATCAAG

gi|647886       880 W  M  A  L  E  S  I  L  H  R  I  Y  T  H  Q  S  D  V  W  S  
gi|647886      2793 TGGATGGCTTTGGAATCAATTTTACACCGAATTTATACACACCAAAGCGACGTCTGGAGC

gi|647886       900 Y  G  V  T  V  W  E  L  M  T  F  G  S  K  P  Y  D  G  I  P  
gi|647886      2853 TATGGAGTCACCGTGTGGGAACTGATGACCTTTGGGTCCAAGCCTTATGATGGGATCCCT

gi|647886       920 A  S  E  I  S  S  I  L  E  K  G  E  R  L  P  Q  P  P  I  C  
gi|647886      2913 GCAAGTGAGATCTCATCCATCCTAGAGAAAGGAGAGCGCCTTCCACAGCCACCTATCTGC

gi|647886       940 T  I  D  V  Y  M  I  M  V  K  C  W  M  I  D  A  D  S  R  P  
gi|647886      2973 ACCATCGACGTCTACATGATCATGGTCAAGTGCTGGATGATAGATGCTGATAGCCGCCCA

gi|647886       960 K  F  R  E  L  I  L  E  F  S  K  M  A  R  D  P  Q  R  Y  L  
gi|647886      3033 AAGTTCCGAGAGTTGATTCTCGAATTCTCCAAAATGGCCAGAGACCCACAGCGCTACCTT

gi|647886       980 V  I  Q  G  D  E  R  M  H  L  P  S  P  T  D  S  N  F  Y  R  
gi|647886      3093 GTTATCCAGGGGGATGAAAGGATGCATTTGCCGAGCCCTACAGACTCCAACTTTTACCGA

gi|647886      1000 A  L  M  E  E  E  D  M  E  D  V  V  D  A  D  E  Y  L  I  P  
gi|647886      3153 GCCCTGATGGAGGAGGAGGACATGGAAGACGTAGTTGATGCTGATGAATACCTCATCCCA

gi|647886      1020 Q  Q  G  F  F  N  S  P  S  T  S  R  T  P  L  L  S  S  L  S  
gi|647886      3213 CAGCAAGGCTTCTTCAACAGCCCATCCACGTCACGGACTCCACTCTTGAGCTCTCTGAGT

gi|647886      1040 A  N  S  N  S  S  T  V  A  C  I  N  R  N  G  S  C  R  V  K  
gi|647886      3273 GCAAATAGCAACAGTTCCACTGTGGCTTGCATTAATAGAAATGGGAGCTGCCGTGTCAAA

gi|647886      1060 E  D  A  F  L  Q  R  Y  S  S  D  P  T  S  V  L  T  E  D  N  
gi|647886      3333 GAAGACGCCTTCTTGCAACGGTATAGCTCCGATCCCACCAGCGTCCTGACAGAGGACAAC

gi|647886      1080 I  D  D  T  F  L  P  V  P  E  Y  I  N  Q  S  V  P  K  R  P  
gi|647886      3393 ATAGATGACACATTCCTTCCCGTGCCTGAATATATAAACCAATCTGTTCCCAAGAGGCCG

gi|647886      1100 A  G  S  V  Q  N  P  V  Y  H  N  Q  P  L  H  P  A  P  G  R  
gi|647886      3453 GCTGGCTCTGTGCAGAACCCAGTCTATCACAATCAGCCCCTGCATCCAGCTCCTGGAAGA

gi|647886      1120 D  L  H  Y  Q  N  P  H  S  N  A  V  S  N  P  E  Y  L  N  T  
gi|647886      3513 GACCTGCATTATCAAAATCCCCATAGCAATGCGGTGAGCAACCCTGAGTATCTCAACACT

gi|647886      1140 A  Q  P  T  C  L  S  S  G  F  D  S  S  A  L  W  I  Q  K  G  
gi|647886      3573 GCCCAGCCGACCTGCCTCAGTAGTGGGTTTGACAGCTCTGCCCTCTGGATCCAGAAAGGC

gi|647886      1160 S  H  Q  M  S  L  D  N  P  D  Y  Q  Q  D  F  F  P  K  E  A  
gi|647886      3633 AGCCACCAAATGAGCCTGGACAACCCTGACTACCAGCAGGACTTCTTTCCCAAAGAAGCC

gi|647886      1180 K  P  N  G  I  F  K  G  P  T  A  E  N  A  E  Y  L  R  V  A  
gi|647886      3693 AAGCCGAATGGCATCTTTAAGGGCCCCACAGCTGAAAATGCAGAGTACCTGCGGGTGGCA

gi|647886      1200 P  P  S  S  E  F  S  G  A   1209
gi|647886      3753 CCGCCAAGCAGTGAGTTTAGTGGAGCA 3780
""",
        )
        alignment = protein_alignment.mapall(codon_alignments)
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
            np.array([[  84,  105,  114,  129,  171,  183,  198,  198,  210,
                        234,  366,  369,  492,  507,  531,  531,  741,  741,
                        807,  810, 1038, 1038, 1089, 1089, 1107, 1107, 1161,
                       1161, 1245, 1260, 1281, 1290, 1443, 1446, 1779, 1779,
                       1821, 1821, 1902, 1953, 1986, 2085, 2097, 2151, 2175,
                       2301, 2319, 2319, 2355, 2358, 2373, 2394, 2418, 2559,
                       2610, 2613, 2655, 2661, 3585, 3594, 3642, 3645, 3765,
                       3765, 3777, 3777, 3777, 3777, 3843, 3843, 3858, 3864,
                       3894, 3900, 3921, 3921, 3969, 3975, 4038, 4038, 4215],
                      [  22,   22,   31,   46,   88,  100,  115,  283,  295,
                        319,  451,  454,  577,  592,  616,  616,  826,  826,
                        892,  895, 1123, 1123, 1174, 1174, 1192, 1192, 1246,
                       1246, 1330, 1345, 1366, 1375, 1528, 1531, 1864, 1864,
                       1906, 1906, 1987, 2038, 2071, 2170, 2182, 2236, 2260,
                       2386, 2404, 2404, 2440, 2443, 2458, 2479, 2503, 2644,
                       2695, 2698, 2740, 2746, 3670, 3679, 3727, 3730, 3850,
                       3850, 3862, 3862, 3862, 3862, 3928, 3928, 3943, 3949,
                       3979, 3985, 4006, 4006, 4054, 4060, 4123, 4123, 4300],
                      [   0,    0,    0,   15,   15,   15,   15,   15,   15,
                         15,  147,  147,  270,  270,  294,  312,  522,  525,
                        591,  591,  819,  822,  873,  879,  897,  903,  957,
                        981, 1065, 1065, 1086, 1086, 1239, 1239, 1572, 1578,
                       1620, 1638, 1719, 1719, 1752, 1752, 1764, 1764, 1788,
                       1788, 1806, 1809, 1845, 1845, 1860, 1860, 1884, 1884,
                       1935, 1935, 1977, 1977, 2901, 2901, 2949, 2949, 3069,
                       3081, 3093, 3093, 3093, 3099, 3165, 3183, 3198, 3198,
                       3228, 3228, 3249, 3258, 3306, 3306, 3369, 3378, 3555],
                      [ 126,  126,  126,  141,  183,  183,  198,  198,  210,
                        210,  342,  342,  465,  465,  489,  507,  717,  720,
                        786,  786, 1014, 1017, 1068, 1074, 1092, 1098, 1152,
                       1176, 1260, 1260, 1281, 1281, 1434, 1434, 1767, 1773,
                       1815, 1833, 1914, 1914, 1947, 1947, 1959, 1959, 1983,
                       1983, 2001, 2004, 2040, 2040, 2055, 2055, 2079, 2079,
                       2130, 2130, 2172, 2172, 3096, 3096, 3144, 3144, 3264,
                       3276, 3288, 3288, 3291, 3297, 3363, 3381, 3396, 3396,
                       3426, 3426, 3447, 3456, 3504, 3504, 3567, 3576, 3753],
                      [ 246,  246,  246,  261,  303,  303,  318,  318,  330,
                        330,  462,  462,  585,  585,  609,  627,  837,  840,
                        906,  906, 1134, 1137, 1188, 1194, 1212, 1218, 1272,
                       1296, 1380, 1380, 1401, 1401, 1554, 1554, 1887, 1893,
                       1935, 1953, 2034, 2034, 2067, 2067, 2079, 2079, 2103,
                       2103, 2121, 2124, 2160, 2160, 2175, 2175, 2199, 2199,
                       2250, 2250, 2292, 2292, 3216, 3216, 3264, 3264, 3384,
                       3396, 3408, 3411, 3414, 3420, 3486, 3504, 3519, 3519,
                       3549, 3549, 3570, 3579, 3627, 3627, 3690, 3699, 3876],
                      [ 153,  153,  153,  168,  210,  210,  225,  225,  237,
                        237,  369,  369,  492,  492,  516,  534,  744,  747,
                        813,  813, 1041, 1044, 1095, 1101, 1119, 1125, 1179,
                       1203, 1287, 1287, 1308, 1308, 1461, 1461, 1794, 1800,
                       1842, 1860, 1941, 1941, 1974, 1974, 1986, 1986, 2010,
                       2010, 2028, 2031, 2067, 2070, 2085, 2085, 2109, 2109,
                       2160, 2160, 2202, 2202, 3126, 3126, 3174, 3174, 3294,
                       3306, 3318, 3318, 3318, 3324, 3390, 3408, 3423, 3423,
                       3453, 3453, 3474, 3483, 3531, 3531, 3594, 3603, 3780]])
                # fmt: on
            )
        )
        self.assertEqual(
            format(alignment, "clustal"),
            """\
gi|24657088|ref|NM_057410.3|        ATGATGATTATCAGCATGTGGATGAGCATATCGCGAGGATTGTGGGACAG
gi|24657104|ref|NM_057411.3|        ---------------------ATGCTGCTGCGACGGCGCAACGGCCCCTG
gi|302179500|gb|HM749883.1|         ------------------------------ATGAAAAAGCACGAG-----
gi|47522839|ref|NM_214007.1|        ------------------------------ATGCGACGCTCCTGGGCGGG
gi|41327737|ref|NM_005228.3|        ------------------------------ATGCGACCCTCCGGGACGGC
gi|6478867|gb|M37394.2|RATEGFR      ------------------------------ATGCGACCCTCAGGGACTGC

gi|24657088|ref|NM_057410.3|        CAGCTCCATCTGGTCGGTGCTGCTGATCCTCGCCTGCATGGCATCCATCA
gi|24657104|ref|NM_057411.3|        CCCCTTCCCCCTGCTGCTCCTGCTCCTGGCCCACTGCATTTGCATTTGGC
gi|302179500|gb|HM749883.1|         --------------------------------------------------
gi|47522839|ref|NM_214007.1|        CGGCGCCGCGCTCCTGGCGCTGCTGGCCGCGCACTTC------------C
gi|41327737|ref|NM_005228.3|        CGGGGCAGCGCTCCTGGCGCTGCTGGCTGCGCTCTGC------------C
gi|6478867|gb|M37394.2|RATEGFR      GAGAACCAAGCTACTGCTGCTGCTGGCTGCGCTCTGC------------G

gi|24657088|ref|NM_057410.3|        CCACAAGCTCATCG------------------------------------
gi|24657104|ref|NM_057411.3|        CCGCGTCGGCGGCCCGCGATCGCTACGCCCGCCAGAACAATCGCCAGCGC
gi|302179500|gb|HM749883.1|         --------------------------------------------------
gi|47522839|ref|NM_214007.1|        AGGCGAGTCCGGCG------------------------------------
gi|41327737|ref|NM_005228.3|        CGGCGAGTCGGGCT------------------------------------
gi|6478867|gb|M37394.2|RATEGFR      CCGCAGGTGGGGCG------------------------------------

gi|24657088|ref|NM_057410.3|        --------------------------------------------------
gi|24657104|ref|NM_057411.3|        CATCAGGATATAGATCGCGATCGGGATCGAGATCGATTCCTATACCGCAG
gi|302179500|gb|HM749883.1|         --------------------------------------------------
gi|47522839|ref|NM_214007.1|        --------------------------------------------------
gi|41327737|ref|NM_005228.3|        --------------------------------------------------
gi|6478867|gb|M37394.2|RATEGFR      --------------------------------------------------

gi|24657088|ref|NM_057410.3|        --------------------------------------------------
gi|24657104|ref|NM_057411.3|        CAGTTCGGCCCAAAATCGACAGAGGGGCGGGGCCAACTTCGCCCTGGGAC
gi|302179500|gb|HM749883.1|         --------------------------------------------------
gi|47522839|ref|NM_214007.1|        --------------------------------------------------
gi|41327737|ref|NM_005228.3|        --------------------------------------------------
gi|6478867|gb|M37394.2|RATEGFR      --------------------------------------------------

gi|24657088|ref|NM_057410.3|        --------------------------------GTCAGCAATGCCGGCTAT
gi|24657104|ref|NM_057411.3|        TGGGAGCCAACGGAGTCACCATTCCCACCAGTCTGGAGGATAAGAACAAG
gi|302179500|gb|HM749883.1|         --------------------------------------------------
gi|47522839|ref|NM_214007.1|        --------------------------------CTGGAGGAGAAG------
gi|41327737|ref|NM_005228.3|        --------------------------------CTGGAGGAAAAG------
gi|6478867|gb|M37394.2|RATEGFR      --------------------------------CTGGAGGAAAAG------

gi|24657088|ref|NM_057410.3|        GTGGATAATGGCAATATGAAAGTCTGCATCGGCACTAAATCTCGGCTCTC
gi|24657104|ref|NM_057411.3|        AACGAGTTCGTCAAGGGGAAAATCTGCATCGGCACTAAATCTCGGCTCTC
gi|302179500|gb|HM749883.1|         ------------------TTACTTTGCCAAGGGACAAGTAACAAGCTCAC
gi|47522839|ref|NM_214007.1|        ------------------AAAGTTTGCCAAGGTACAAGTAACAAGCTCAC
gi|41327737|ref|NM_005228.3|        ------------------AAAGTTTGCCAAGGCACGAGTAACAAGCTCAC
gi|6478867|gb|M37394.2|RATEGFR      ------------------AAAGTTTGCCAAGGCACAAGTAACAGGCTCAC

gi|24657088|ref|NM_057410.3|        CGTGCCCTCCAACAAGGAACATCATTACCGGAACCTCAGAGATCGGTACA
gi|24657104|ref|NM_057411.3|        CGTGCCCTCCAACAAGGAACATCATTACCGGAACCTCAGAGATCGGTACA
gi|302179500|gb|HM749883.1|         CCAGTTGGGCACTTTTGAAGACCACTTTCTGAGCCTACAGAGGATGTTCA
gi|47522839|ref|NM_214007.1|        CCAGCTGGGCACTTTCGAAGACCACTTTCTGAGCCTCCAGAGGATGTTCA
gi|41327737|ref|NM_005228.3|        GCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGATGTTCA
gi|6478867|gb|M37394.2|RATEGFR      CCAACTAGGCACCTTTGAAGACCACTTTCTGAGCCTCCAGAGGATGTTCA

gi|24657088|ref|NM_057410.3|        CGAACTGTACGTATGTGGATGGCAACCTGGAGCTGACCTGGCTGCCCAAC
gi|24657104|ref|NM_057411.3|        CGAACTGTACGTATGTGGATGGCAACCTGGAGCTGACCTGGCTGCCCAAC
gi|302179500|gb|HM749883.1|         ACAACTGTGAGGTGGTCCTTGGGAATTTGGAAATTACCTACATGCAGAGT
gi|47522839|ref|NM_214007.1|        ATAACTGCGAGGTGGTCCTTGGGAACTTGGAGATCACCTACATGCAGAAC
gi|41327737|ref|NM_005228.3|        ATAACTGTGAGGTGGTCCTTGGGAATTTGGAAATTACCTATGTGCAGAGG
gi|6478867|gb|M37394.2|RATEGFR      ACAACTGTGAAGTGGTCCTTGGAAACTTGGAAATCACCTATGTGCAAAGG

gi|24657088|ref|NM_057410.3|        GAGAATTTGGACCTCAGCTTCCTGGACAACATACGGGAGGTCACCGGCTA
gi|24657104|ref|NM_057411.3|        GAGAATTTGGACCTCAGCTTCCTGGACAACATACGGGAGGTCACCGGCTA
gi|302179500|gb|HM749883.1|         ---AGTTACAACCTTTCTTTTCTCAAGACCATCCAGGAGGTTGCCGGCTA
gi|47522839|ref|NM_214007.1|        ---AGCTACAACCTGTCTTTCCTAAAGACCATTCAGGAGGTCGCCGGCTA
gi|41327737|ref|NM_005228.3|        ---AATTATGATCTTTCCTTCTTAAAGACCATCCAGGAGGTGGCTGGTTA
gi|6478867|gb|M37394.2|RATEGFR      ---AATTATGACCTTTCCTTCTTAAAGACCATCCAGGAGGTGGCTGGCTA

gi|24657088|ref|NM_057410.3|        TATTCTGATCAGTCATGTGGACGTTAAGAAAGTGGTATTTCCCAAACTAC
gi|24657104|ref|NM_057411.3|        TATTCTGATCAGTCATGTGGACGTTAAGAAAGTGGTATTTCCCAAACTAC
gi|302179500|gb|HM749883.1|         TGTACTCATTGCCCTCAACACAGTGGAGAAGATTCCGCTGGAAAACCTGC
gi|47522839|ref|NM_214007.1|        CGTGCTCATCGCCCTCAACACCGTGGAGAAGATCCCTTTGGAAAACCTGC
gi|41327737|ref|NM_005228.3|        TGTCCTCATTGCCCTCAACACAGTGGAGCGAATTCCTTTGGAAAACCTGC
gi|6478867|gb|M37394.2|RATEGFR      TGTTCTCATTGCCCTGAACACCGTGGAGAGAATCCCTTTGGAGAACCTGC

gi|24657088|ref|NM_057410.3|        AAATCATTCGCGGACGCACGCTGTTCAGCTTATCCGTGGAGGAGGAGAAG
gi|24657104|ref|NM_057411.3|        AAATCATTCGCGGACGCACGCTGTTCAGCTTATCCGTGGAGGAGGAGAAG
gi|302179500|gb|HM749883.1|         AGATCATCCGAGGAAATGTGCTTTAT---------------GAAAACACC
gi|47522839|ref|NM_214007.1|        AGATCATCCGAGGAAATGTACTGTAT---------------GAAAACACC
gi|41327737|ref|NM_005228.3|        AGATCATCAGAGGAAATATGTACTAC---------------GAAAATTCC
gi|6478867|gb|M37394.2|RATEGFR      AGATCATCAGGGGAAATGCTCTCTAC---------------GAAAACACC

gi|24657088|ref|NM_057410.3|        TATGCCTTGTTCGTC------------------ACTTATTCCAAAATGTA
gi|24657104|ref|NM_057411.3|        TATGCCTTGTTCGTC------------------ACTTATTCCAAAATGTA
gi|302179500|gb|HM749883.1|         CATGCCTTAGCCGTCTTATCCAACTATGGAGCAAACAAAACCGGACTGAG
gi|47522839|ref|NM_214007.1|        CATGCCTTAGCCGTCTTATCCAACTACGGGGCCAATAAAACCGGCCTGAG
gi|41327737|ref|NM_005228.3|        TATGCCTTAGCAGTCTTATCTAACTATGATGCAAATAAAACCGGACTGAA
gi|6478867|gb|M37394.2|RATEGFR      TACGCCTTAGCCGTCCTGTCCAACTATGGAACCAACAAAACTGGGCTTAG

gi|24657088|ref|NM_057410.3|        CACGCTGGAGATTCCCGATCTACGCGATGTCTTAAATGGCCAAGTGGGCT
gi|24657104|ref|NM_057411.3|        CACGCTGGAGATTCCCGATCTACGCGATGTCTTAAATGGCCAAGTGGGCT
gi|302179500|gb|HM749883.1|         GGAGCTGCCCTTGAGAAACTTACAGGAAATTCTGCAAGGTGCCGTGAGAT
gi|47522839|ref|NM_214007.1|        GGAGCTGCCCATGAGGAACTTACAAGAGATCCTGCAAGGCGCCGTGCGCT
gi|41327737|ref|NM_005228.3|        GGAGCTGCCCATGAGAAATTTACAGGAAATCCTGCATGGCGCCGTGCGGT
gi|6478867|gb|M37394.2|RATEGFR      GGAACTGCCCATGCGGAACTTACAGGAAATTCTGATCGGTGCTGTGCGAT

gi|24657088|ref|NM_057410.3|        TCCACAACAACTACAATCTCTGCCACATGCGAACGATCCAGTGGTCGGAG
gi|24657104|ref|NM_057411.3|        TCCACAACAACTACAATCTCTGCCACATGCGAACGATCCAGTGGTCGGAG
gi|302179500|gb|HM749883.1|         TCAGCAACAACCCTGTCCTCTGCAACGTGGAGACCATCCAGTGGCGGGAC
gi|47522839|ref|NM_214007.1|        TCAGCAACAACCCTGCCCTCTGTCACGCGGAGTCCATCCAGTGGAGGGAC
gi|41327737|ref|NM_005228.3|        TCAGCAACAACCCTGCCCTGTGCAACGTGGAGAGCATCCAGTGGCGGGAC
gi|6478867|gb|M37394.2|RATEGFR      TTAGCAACAACCCCATCCTCTGCAATATGGAGACCATCCAGTGGAGGGAC

gi|24657088|ref|NM_057410.3|        ATTGTATCCAACGGCACGGATGCATACTACAACTACGACTTTACTGCTCC
gi|24657104|ref|NM_057411.3|        ATTGTATCCAACGGCACGGATGCATACTACAACTACGACTTTACTGCTCC
gi|302179500|gb|HM749883.1|         ATCGTCAACCCTGATTTTCTAAGCAACATGACAGGGGACTTTCAGAACCA
gi|47522839|ref|NM_214007.1|        ATTGTCAACAGCGACTTTCTAAGCAACATGTCCATGGACTTTCAGAGCCA
gi|41327737|ref|NM_005228.3|        ATAGTCAGCAGTGACTTTCTCAGCAACATGTCGATGGACTTCCAGAACCA
gi|6478867|gb|M37394.2|RATEGFR      ATCGTCCAAGATGTCTTTCTGAGCAACATGTCAATGGACGTACAGCGCCA

gi|24657088|ref|NM_057410.3|        GGAGCGCGAGTGTCCCAAGTGCCACGAGAGCTGCACGCACGGA---TGTT
gi|24657104|ref|NM_057411.3|        GGAGCGCGAGTGTCCCAAGTGCCACGAGAGCTGCACGCACGGA---TGTT
gi|302179500|gb|HM749883.1|         GCAGGGCAACTGCCCAAAGTGTGATCCAGCCTGTCTCAACAGAAGCTGCT
gi|47522839|ref|NM_214007.1|        GCTGGGCAGCTGCCCGAAGTGTGATCCAGGCTGTCTCAATGGGAGCTGCT
gi|41327737|ref|NM_005228.3|        CCTGGGCAGCTGCCAAAAGTGTGATCCAAGCTGTCCCAATGGGAGCTGCT
gi|6478867|gb|M37394.2|RATEGFR      CCTGACGGGCTGCCCGAAATGTGATCCGAGCTGTCCCAATGGAAGCTGCT

gi|24657088|ref|NM_057410.3|        GGGGCGAGGGTCCCAAGAATTGCCAGAAGTTCAGCAAGCTCACCTGCTCG
gi|24657104|ref|NM_057411.3|        GGGGCGAGGGTCCCAAGAATTGCCAGAAGTTCAGCAAGCTCACCTGCTCG
gi|302179500|gb|HM749883.1|         GGGGTGCCGGGGAGGAGAACTGTCAGAAATTGACCAAAATCATCTGTGCC
gi|47522839|ref|NM_214007.1|        GGGGTGCTGGGAAGGAGAACTGCCAGAAATTGACCAAAGTCATCTGTGCC
gi|41327737|ref|NM_005228.3|        GGGGTGCAGGAGAGGAGAACTGCCAGAAACTGACCAAAATCATCTGTGCC
gi|6478867|gb|M37394.2|RATEGFR      GGGGAAGAGGAGAGGAGAACTGCCAGAAATTGACCAAAATCATCTGCGCC

gi|24657088|ref|NM_057410.3|        CCACAGTGTGCCGGAGGTCGTTGCTATGGACCAAAGCCGCGGGAGTGTTG
gi|24657104|ref|NM_057411.3|        CCACAGTGTGCCGGAGGTCGTTGCTATGGACCAAAGCCGCGGGAGTGTTG
gi|302179500|gb|HM749883.1|         CAGCAGTGTTCC---GGGCGCTGCCGTGGCAGGTCCCCCAGTGACTGCTG
gi|47522839|ref|NM_214007.1|        CAGCAGTGCTCC---GGGCGCTGCCGCGGCCGGTCGCCCAGTGACTGCTG
gi|41327737|ref|NM_005228.3|        CAGCAGTGCTCC---GGGCGCTGCCGTGGCAAGTCCCCCAGTGACTGCTG
gi|6478867|gb|M37394.2|RATEGFR      CAGCAATGTTCC---CGGCGTTGTCGTGGCAGGTCCCCTAGCGACTGCTG

gi|24657088|ref|NM_057410.3|        TCACCTCTTCTGCGCCGGAGGATGCACTGGTCCCACGCAAAAGGATTGCA
gi|24657104|ref|NM_057411.3|        TCACCTCTTCTGCGCCGGAGGATGCACTGGTCCCACGCAAAAGGATTGCA
gi|302179500|gb|HM749883.1|         CCACAACCAGTGTGCCGCTGGCTGCACAGGGCCACGGGAGAGCGACTGCC
gi|47522839|ref|NM_214007.1|        CCACAACCAGTGCGCCGCTGGCTGCACGGGGCCGCGGGAGAGCGACTGCC
gi|41327737|ref|NM_005228.3|        CCACAACCAGTGTGCTGCAGGCTGCACAGGCCCCCGGGAGAGCGACTGCC
gi|6478867|gb|M37394.2|RATEGFR      CCACAACCAGTGTGCCGCAGGGTGTACAGGGCCCAGAGAGAGTGACTGTC

gi|24657088|ref|NM_057410.3|        TCGCCTGCAAGAACTTCTTCGACGAGGGCGTATGCAAGGAGGAATGCCCG
gi|24657104|ref|NM_057411.3|        TCGCCTGCAAGAACTTCTTCGACGAGGGCGTATGCAAGGAGGAATGCCCG
gi|302179500|gb|HM749883.1|         TGGTCTGCCGCAGGTTCCGTGATGAAGCCACCTGCAAGGACACGTGTCCG
gi|47522839|ref|NM_214007.1|        TGGTTTGCCGCAGATTCCGTGACGAGGCCACCTGCAAGGACACATGCCCG
gi|41327737|ref|NM_005228.3|        TGGTCTGCCGCAAATTCCGAGACGAAGCCACGTGCAAGGACACCTGCCCC
gi|6478867|gb|M37394.2|RATEGFR      TGGTCTGCCACAGGTTCCGAGATGAAGCCACGTGCAAAGACACCTGCCCA

gi|24657088|ref|NM_057410.3|        CCCATGCGCAAGTACAATCCCACCACCTATGTTCTTGAAACGAATCCTGA
gi|24657104|ref|NM_057411.3|        CCCATGCGCAAGTACAATCCCACCACCTATGTTCTTGAAACGAATCCTGA
gi|302179500|gb|HM749883.1|         CCACTCATGCTCTATGACCCTACCACCTACGAAATGAAGGTCAACCCGCT
gi|47522839|ref|NM_214007.1|        CCGCTCATGCTCTACAACCCCACCACCTACCAGATGGACGTCAACCCGCT
gi|41327737|ref|NM_005228.3|        CCACTCATGCTCTACAACCCCACCACGTACCAGATGGATGTGAACCCCGA
gi|6478867|gb|M37394.2|RATEGFR      CCACTCATGCTGTACAACCCCACCACGTACCAGATGGATGTCAACCCTGA

gi|24657088|ref|NM_057410.3|        GGGAAAGTATGCCTATGGTGCCACCTGCGTCAAGGAGTGTCCC---GGTC
gi|24657104|ref|NM_057411.3|        GGGAAAGTATGCCTATGGTGCCACCTGCGTCAAGGAGTGTCCC---GGTC
gi|302179500|gb|HM749883.1|         GGGGAAGTACAGCTTTGGCGCCACCTGTGTCAAGAAGTGTCCCCGTAACT
gi|47522839|ref|NM_214007.1|        GGGGAAGTACAGCTTTGGCGCCACCTGTGTCAAGAAGTGCCCTCGTAACT
gi|41327737|ref|NM_005228.3|        GGGCAAATACAGCTTTGGTGCCACCTGCGTGAAGAAGTGTCCCCGTAATT
gi|6478867|gb|M37394.2|RATEGFR      GGGGAAGTACAGCTTTGGTGCCACCTGTGTGAAGAAATGCCCCAGAAACT

gi|24657088|ref|NM_057410.3|        ATCTGTTGCGTGATAATGGCGCCTGCGTGCGCAGCTGTCCCCAGGAC---
gi|24657104|ref|NM_057411.3|        ATCTGTTGCGTGATAATGGCGCCTGCGTGCGCAGCTGTCCCCAGGAC---
gi|302179500|gb|HM749883.1|         ACGTGGTGACAGACCACGGCTCCTGCGTCCGCGCCTGCAGTTCTGACAGC
gi|47522839|ref|NM_214007.1|        ACGTGGTGACAGACCATGGCTCCTGTGTCCGTGCCTGCAGCTCCGACAGC
gi|41327737|ref|NM_005228.3|        ATGTGGTGACAGATCACGGCTCGTGCGTCCGAGCCTGTGGGGCCGACAGC
gi|6478867|gb|M37394.2|RATEGFR      ACGTGGTGACAGATCACGGCTCGTGTGTCCGGGCCTGTGGGCCAGACTAC

gi|24657088|ref|NM_057410.3|        ---AAGATGGACAAGGGGGGC------GAGTGTGTGCCCTGCAATGGACC
gi|24657104|ref|NM_057411.3|        ---AAGATGGACAAGGGGGGC------GAGTGTGTGCCCTGCAATGGACC
gi|302179500|gb|HM749883.1|         CAGGAGGTAGAGGAAGACGGTGTCCGCAAGTGTAAAAAGTGTGACGGGCC
gi|47522839|ref|NM_214007.1|        TACGAGGTGGAGGAGGACGGCGTCCGCAAGTGTAAAAAGTGTGACGGGCC
gi|41327737|ref|NM_005228.3|        TATGAGATGGAGGAAGACGGCGTCCGCAAGTGTAAGAAGTGCGAAGGGCC
gi|6478867|gb|M37394.2|RATEGFR      TATGAAGTAGAAGAAGATGGAGTCAGCAAGTGTAAAAAATGTGACGGGCC

gi|24657088|ref|NM_057410.3|        GTGCCCCAAAACCTGCCCGGGCGTTACTGTC-------------------
gi|24657104|ref|NM_057411.3|        GTGCCCCAAAACCTGCCCGGGCGTTACTGTC-------------------
gi|302179500|gb|HM749883.1|         TTGTGGCAAAGTTTGTAACGGAATAGGAATCGGTGAGTTTAAAGACACAC
gi|47522839|ref|NM_214007.1|        CTGCGGCAAAGTTTGTAACGGGATAGGGATTGGCGAGTTTAAAGACACAC
gi|41327737|ref|NM_005228.3|        TTGCCGCAAAGTGTGTAACGGAATAGGTATTGGTGAATTTAAAGACTCAC
gi|6478867|gb|M37394.2|RATEGFR      CTGCCGCAAAGTTTGCAATGGCATAGGCATTGGTGAATTTAAAGACACAC

gi|24657088|ref|NM_057410.3|        -----CTGCATGCCGGCAACATTGACTCGTTCCGGAATTGTACGGTGATC
gi|24657104|ref|NM_057411.3|        -----CTGCATGCCGGCAACATTGACTCGTTCCGGAATTGTACGGTGATC
gi|302179500|gb|HM749883.1|         TTTCCATAAATGCTACAAACATTAAACACTTCAGAAACTGCACATCCATC
gi|47522839|ref|NM_214007.1|        TTTCCATAAATGCTACGAATATCAAGCACTTCAGGAACTGCACCTCGATC
gi|41327737|ref|NM_005228.3|        TCTCCATAAATGCTACGAATATTAAACACTTCAAAAACTGCACCTCCATC
gi|6478867|gb|M37394.2|RATEGFR      TCTCCATAAATGCTACAAACATCAAACACTTCAAGTACTGCACTGCCATC

gi|24657088|ref|NM_057410.3|        GATGGCAACATTCGCATTTTGGATCAGACCTTCTCGGGCTTCCAGGATGT
gi|24657104|ref|NM_057411.3|        GATGGCAACATTCGCATTTTGGATCAGACCTTCTCGGGCTTCCAGGATGT
gi|302179500|gb|HM749883.1|         AGTGGAGATCTTCATATCCTGCCAGTAGCATTTAGGGGT-----------
gi|47522839|ref|NM_214007.1|        AGCGGAGATCTTCATATCCTGCCGGTAGCATTTAGGGGT-----------
gi|41327737|ref|NM_005228.3|        AGTGGCGATCTCCACATCCTGCCGGTGGCATTTAGGGGT-----------
gi|6478867|gb|M37394.2|RATEGFR      AGTGGGGACCTCCACATCCTGCCAGTGGCCTTTAAGGGG-----------

gi|24657088|ref|NM_057410.3|        CTATGCCAACTACACGATGGGACCACGATACATACCGCTGGATCCCGAGC
gi|24657104|ref|NM_057411.3|        CTATGCCAACTACACGATGGGACCACGATACATACCGCTGGATCCCGAGC
gi|302179500|gb|HM749883.1|         ----GACTCCTTCACACGTACTGCA---------CCTCTGGACCCGAAAG
gi|47522839|ref|NM_214007.1|        ----GACTCCTTCACACGCACGCCG---------CCTCTGGACCCCAAGG
gi|41327737|ref|NM_005228.3|        ----GACTCCTTCACACATACTCCT---------CCTCTGGATCCACAGG
gi|6478867|gb|M37394.2|RATEGFR      ----GATTCTTTCACCCGCACTCCT---------CCTCTAGACCCACGGG

gi|24657088|ref|NM_057410.3|        GACTGGAGGTGTTCTCCACGGTGAAGGAGATCACCGGGTATCTGAATATC
gi|24657104|ref|NM_057411.3|        GACTGGAGGTGTTCTCCACGGTGAAGGAGATCACCGGGTATCTGAATATC
gi|302179500|gb|HM749883.1|         AACTGGACATTCTAAGAACTGTAAAAGAAATAACAGGGTTTTTGCTGATT
gi|47522839|ref|NM_214007.1|        AACTGGACATCCTGAAAACCGTGAAGGAAATAACAGGGTTTTTACTGATT
gi|41327737|ref|NM_005228.3|        AACTGGATATTCTGAAAACCGTAAAGGAAATCACAGGGTTTTTGCTGATT
gi|6478867|gb|M37394.2|RATEGFR      AACTAGAAATTCTCAAAACTGTGAAGGAAATAACAGGGTTTTTGCTGATT

gi|24657088|ref|NM_057410.3|        GAGGGAACCCACCCGCAGTTCCGGAATCTGTCGTACTTCCGCAATCTGGA
gi|24657104|ref|NM_057411.3|        GAGGGAACCCACCCGCAGTTCCGGAATCTGTCGTACTTCCGCAATCTGGA
gi|302179500|gb|HM749883.1|         CAGGCCTGGCCCGAAAACAGGACTGACCTCCATGCTTTTGAGAACCTGGA
gi|47522839|ref|NM_214007.1|        CAGGCCTGGCCTGAAAACAGGACTGGCCTCCATGCTTTTGAGAACCTGGA
gi|41327737|ref|NM_005228.3|        CAGGCTTGGCCTGAAAACAGGACGGACCTCCATGCCTTTGAGAACCTAGA
gi|6478867|gb|M37394.2|RATEGFR      CAGGCTTGGCCTGAAAACTGGACTGACCTCCATGCTTTTGAGAACCTAGA

gi|24657088|ref|NM_057410.3|        AACAATTCATGGCCGCCAGCTGATGGAGAGCATGTTTGCCGCTTTGGCGA
gi|24657104|ref|NM_057411.3|        AACAATTCATGGCCGCCAGCTGATGGAGAGCATGTTTGCCGCTTTGGCGA
gi|302179500|gb|HM749883.1|         AATCATACGTGGCAGAACGAAGCAGCATGGCCAGTTT---TCTCTTGCGG
gi|47522839|ref|NM_214007.1|        AATCATACGTGGCAGGACGAAGCAACATGGTCAGTTT---TCCCTCGCGG
gi|41327737|ref|NM_005228.3|        AATCATACGCGGCAGGACCAAGCAACATGGTCAGTTT---TCTCTTGCAG
gi|6478867|gb|M37394.2|RATEGFR      AATAATTCGTGGCAGAACAAAGCAACATGGTCAGTTT---TCTCTGGCGG

gi|24657088|ref|NM_057410.3|        TCGTTAAGTCATCCCTGTACAGCCTGGAGATGCGCAATCTGAAGCAGATT
gi|24657104|ref|NM_057411.3|        TCGTTAAGTCATCCCTGTACAGCCTGGAGATGCGCAATCTGAAGCAGATT
gi|302179500|gb|HM749883.1|         TTGTCGGCCTGGATATAACATCTTTGGGATTACGCTCCCTCAAGGAGATA
gi|47522839|ref|NM_214007.1|        TTGTTGGCCTGGACATAGCGTCCTTGGGGATGCGCTCCCTCAAGGAGATC
gi|41327737|ref|NM_005228.3|        TCGTCAGCCTGAACATAACATCCTTGGGATTACGCTCCCTCAAGGAGATA
gi|6478867|gb|M37394.2|RATEGFR      TTGTCGGCCTGAACATAACATCGCTGGGGTTGCGTTCCCTCAAGGAGATC

gi|24657088|ref|NM_057410.3|        AGTTCCGGCAGTGTGGTCATCCAGCATAATAGAGACCTCTGCTACGTAAG
gi|24657104|ref|NM_057411.3|        AGTTCCGGCAGTGTGGTCATCCAGCATAATAGAGACCTCTGCTACGTAAG
gi|302179500|gb|HM749883.1|         AGTGATGGTGATGTGATAATTTCAGGAAATCGAAACTTGTGCTATGCAGA
gi|47522839|ref|NM_214007.1|        AGCGACGGAGACGTGATCGTCTCAGGAAACCGAAACCTGTGCTATGCAAA
gi|41327737|ref|NM_005228.3|        AGTGATGGAGATGTGATAATTTCAGGAAACAAAAATTTGTGCTATGCAAA
gi|6478867|gb|M37394.2|RATEGFR      AGTGATGGGGATGTGATTATTTCTGGGAACCGAAATTTGTGCTACGCAAA

gi|24657088|ref|NM_057410.3|        CAATATCCGTTGGCCGGCCATTCAGAAGGAGCCCGAACAGAAGGTGTGGG
gi|24657104|ref|NM_057411.3|        CAATATCCGTTGGCCGGCCATTCAGAAGGAGCCCGAACAGAAGGTGTGGG
gi|302179500|gb|HM749883.1|         TACAATACGCTGGAAAAAACTTTTTGGGACCTCAACTCAGAAAACCAAAA
gi|47522839|ref|NM_214007.1|        TACAATCAGCTGGAAAAAACTATTTGGGACCGCAAGTCAGAAAACCAAAA
gi|41327737|ref|NM_005228.3|        TACAATAAACTGGAAAAAACTGTTTGGGACCTCCGGTCAGAAAACCAAAA
gi|6478867|gb|M37394.2|RATEGFR      CACTATAAACTGGAAAAAACTCTTCGGGACGCCCAATCAAAAGACCAAAA

gi|24657088|ref|NM_057410.3|        TCAACGAGAATCTCAGGGCGGATCTATGCGAGAAAAATGGAACCATTTGC
gi|24657104|ref|NM_057411.3|        TCAACGAGAATCTCAGGGCGGATCTATGCGAGAAAAATGGAACCATTTGC
gi|302179500|gb|HM749883.1|         TTTTAAACAACAGGAGTGAAAAACAGTGCAAGGCCGCAGGCCACATCTGT
gi|47522839|ref|NM_214007.1|        TTATAAACAACAGGAGCGAAAAAGAGTGCAAAGCCATGGGCCACATCTGT
gi|41327737|ref|NM_005228.3|        TTATAAGCAACAGAGGTGAAAACAGCTGCAAGGCCACAGGCCAGGTCTGC
gi|6478867|gb|M37394.2|RATEGFR      TCATGAACAACAGAGCTGAAAAGGACTGCAAGGCCACGAACCACGTCTGT

gi|24657088|ref|NM_057410.3|        TCGGATCAGTGCAACGAGGACGGCTGCTGGGGAGCTGGCACGGATCAGTG
gi|24657104|ref|NM_057411.3|        TCGGATCAGTGCAACGAGGACGGCTGCTGGGGAGCTGGCACGGATCAGTG
gi|302179500|gb|HM749883.1|         CACCCGCTGTGCTCATCAGAGGGCTGCTGGGGACCGGGACCCAAATACTG
gi|47522839|ref|NM_214007.1|        AACCCGCTGTGCTCATCAGAGGGCTGCTGGGGCCCTGAACCCAGAGACTG
gi|41327737|ref|NM_005228.3|        CATGCCTTGTGCTCCCCCGAGGGCTGCTGGGGCCCGGAGCCCAGGGACTG
gi|6478867|gb|M37394.2|RATEGFR      AATCCTTTATGCTCCTCGGAAGGCTGCTGGGGCCCTGAGCCCACGGACTG

gi|24657088|ref|NM_057410.3|        CCTTACCTGCAAGAACTTCAATTTCAATGGCACCTGCATCGCCGACTGTG
gi|24657104|ref|NM_057411.3|        CCTTACCTGCAAGAACTTCAATTTCAATGGCACCTGCATCGCCGACTGTG
gi|302179500|gb|HM749883.1|         CATGTCCTGCCAGAACTTCAGTCGTGGCAAGGAGTGTGTGGGAAAGTGCA
gi|47522839|ref|NM_214007.1|        CATGTCCTGTCGAAACTTTAGCCGCGGCAAGGAATGTGTGGAGAAGTGCA
gi|41327737|ref|NM_005228.3|        CGTCTCTTGCCGGAATGTCAGCCGAGGCAGGGAATGCGTGGACAAGTGCA
gi|6478867|gb|M37394.2|RATEGFR      TGTCTCCTGCCAGAATGTGAGCAGAGGCAGGGAGTGCGTGGACAAGTGCA

gi|24657088|ref|NM_057410.3|        GTTATATATCCAATGCCTACAAG------TTTGACAATAGAACGTGCAAG
gi|24657104|ref|NM_057411.3|        GTTATATATCCAATGCCTACAAG------TTTGACAATAGAACGTGCAAG
gi|302179500|gb|HM749883.1|         ACATTCTAGAGGGAGAGCCCAGAGAATTCGTGGAGAACTCCGAGTGTGTG
gi|47522839|ref|NM_214007.1|        ACGTTCTGGAGGGGGAGCCGAGAGAGTTCGTGGAGAATGCCGAGTGTGTG
gi|41327737|ref|NM_005228.3|        ACCTTCTGGAGGGTGAGCCAAGGGAGTTTGTGGAGAACTCTGAGTGCATA
gi|6478867|gb|M37394.2|RATEGFR      ACATCCTGGAGGGGGAACCGAGGGAGTTTGTGGAAAATTCTGAATGCATC

gi|24657088|ref|NM_057410.3|        ATATGCCATCCAGAGTGCCGG------------------ACTTGCAATGG
gi|24657104|ref|NM_057411.3|        ATATGCCATCCAGAGTGCCGG------------------ACTTGCAATGG
gi|302179500|gb|HM749883.1|         CAGTGCCATCCAGAATGCCTGCCCCAGGCCATGAACGTGACCTGCACTGG
gi|47522839|ref|NM_214007.1|        CAGTGCCACCCGGAGTGCCTGCCCCAGGCCAAGAACGTGACCTGCATGGG
gi|41327737|ref|NM_005228.3|        CAGTGCCACCCAGAGTGCCTGCCTCAGGCCATGAACATCACCTGCACAGG
gi|6478867|gb|M37394.2|RATEGFR      CAGTGCCATCCAGAATGTCTGCCCCAGACCATGAACATCACCTGTACAGG

gi|24657088|ref|NM_057410.3|        AGCTGGAGCAGATCACTGCCAGGAGTGCGTCCATGTGAGGGACGGTCAGC
gi|24657104|ref|NM_057411.3|        AGCTGGAGCAGATCACTGCCAGGAGTGCGTCCATGTGAGGGACGGTCAGC
gi|302179500|gb|HM749883.1|         ACGCGGACCAGGCAACTGTGTAAAGTGCGCCCACTACATTGATGGCCCTC
gi|47522839|ref|NM_214007.1|        ACGCGGACCGGACAGCTGTGTCCGGTGTGCTCACTACATCGACGGCCCTC
gi|41327737|ref|NM_005228.3|        ACGGGGACCAGACAACTGTATCCAGTGTGCCCACTACATTGACGGCCCCC
gi|6478867|gb|M37394.2|RATEGFR      CCGGGGGCCAGACAACTGCATCAAGTGTGCCCACTATGTTGATGGTCCCC

gi|24657088|ref|NM_057410.3|        ACTGTGTGTCCGAGTGCCCGAAGAACAAGTACAACGATCGTGGTGTCTGC
gi|24657104|ref|NM_057411.3|        ACTGTGTGTCCGAGTGCCCGAAGAACAAGTACAACGATCGTGGTGTCTGC
gi|302179500|gb|HM749883.1|         ACTGCGTCAAGACCTGCCCT------------------------------
gi|47522839|ref|NM_214007.1|        ACTGTGTCAAGACCTGCCCC------------------------------
gi|41327737|ref|NM_005228.3|        ACTGCGTCAAGACCTGCCCG------------------------------
gi|6478867|gb|M37394.2|RATEGFR      ACTGTGTCAAGACCTGCCCT------------------------------

gi|24657088|ref|NM_057410.3|        CGAGAGTGCCACGCCACCTGCGATGGATGCACTGGGCCCAAGGACACCAT
gi|24657104|ref|NM_057411.3|        CGAGAGTGCCACGCCACCTGCGATGGATGCACTGGGCCCAAGGACACCAT
gi|302179500|gb|HM749883.1|         ---------------------GCTGGAGTCGCGGGAGAGAATGGCACCCT
gi|47522839|ref|NM_214007.1|        ---------------------GCGGGAATCGCAGGAGAAAACAGCACCCT
gi|41327737|ref|NM_005228.3|        ---------------------GCAGGAGTCATGGGAGAAAACAACACCCT
gi|6478867|gb|M37394.2|RATEGFR      ---------------------TCGGGCATCATGGGGGAGAACAACACCCT

gi|24657088|ref|NM_057410.3|        CGGCATTGGAGCGTGTACAACGTGCAATTTGGCCATTATCAACAATGACG
gi|24657104|ref|NM_057411.3|        CGGCATTGGAGCGTGTACAACGTGCAATTTGGCCATTATCAACAATGACG
gi|302179500|gb|HM749883.1|         GATC----------------------------------------------
gi|47522839|ref|NM_214007.1|        CATC----------------------------------------------
gi|41327737|ref|NM_005228.3|        GGTC----------------------------------------------
gi|6478867|gb|M37394.2|RATEGFR      GGTC----------------------------------------------

gi|24657088|ref|NM_057410.3|        CCACAGTAAAACGCTGCCTGCTGAAGGACGACAAGTGCCCCGATGGGTAC
gi|24657104|ref|NM_057411.3|        CCACAGTAAAACGCTGCCTGCTGAAGGACGACAAGTGCCCCGATGGGTAC
gi|302179500|gb|HM749883.1|         --------------------------------------------------
gi|47522839|ref|NM_214007.1|        --------------------------------------------------
gi|41327737|ref|NM_005228.3|        --------------------------------------------------
gi|6478867|gb|M37394.2|RATEGFR      --------------------------------------------------

gi|24657088|ref|NM_057410.3|        TTCTGGGAGTATGTGCATCCACAAGAGCAGGGATCGCTAAAGCCATTGGC
gi|24657104|ref|NM_057411.3|        TTCTGGGAGTATGTGCATCCACAAGAGCAGGGATCGCTAAAGCCATTGGC
gi|302179500|gb|HM749883.1|         ---TGGAAGTTTGCA-----------------------------------
gi|47522839|ref|NM_214007.1|        ---TGGAAGTTTGCG-----------------------------------
gi|41327737|ref|NM_005228.3|        ---TGGAAGTACGCA-----------------------------------
gi|6478867|gb|M37394.2|RATEGFR      ---TGGAAGTTTGCA-----------------------------------

gi|24657088|ref|NM_057410.3|        CGGCAGAGCAGTTTGCCGAAAGTGCCATCCCCTTTGCGAGCTGTGCACCA
gi|24657104|ref|NM_057411.3|        CGGCAGAGCAGTTTGCCGAAAGTGCCATCCCCTTTGCGAGCTGTGCACCA
gi|302179500|gb|HM749883.1|         -------------------GATGCCAACCACGTGTGTCTCCTG-------
gi|47522839|ref|NM_214007.1|        -------------------GATGCCAACCACGTGTGTCACCTG-------
gi|41327737|ref|NM_005228.3|        -------------------GACGCCGGCCATGTGTGCCACCTG-------
gi|6478867|gb|M37394.2|RATEGFR      -------------------GATGCCAATAACGTCTGCCACCTC-------

gi|24657088|ref|NM_057410.3|        ACTACGGATACCATGAACAGGTGTGCTCCAAGTGCACCCACTACAAGCGA
gi|24657104|ref|NM_057411.3|        ACTACGGATACCATGAACAGGTGTGCTCCAAGTGCACCCACTACAAGCGA
gi|302179500|gb|HM749883.1|         --------------------------------------------------
gi|47522839|ref|NM_214007.1|        --------------------------------------------------
gi|41327737|ref|NM_005228.3|        --------------------------------------------------
gi|6478867|gb|M37394.2|RATEGFR      --------------------------------------------------

gi|24657088|ref|NM_057410.3|        CGAGAGCAGTGCGAGACCGAGTGTCCGGCCGATCACTACACGGATGAGGA
gi|24657104|ref|NM_057411.3|        CGAGAGCAGTGCGAGACCGAGTGTCCGGCCGATCACTACACGGATGAGGA
gi|302179500|gb|HM749883.1|         --------------------------------------------------
gi|47522839|ref|NM_214007.1|        --------------------------------------------------
gi|41327737|ref|NM_005228.3|        --------------------------------------------------
gi|6478867|gb|M37394.2|RATEGFR      --------------------------------------------------

gi|24657088|ref|NM_057410.3|        GCAGCGCGAGTGCTTCCAGTGCCACCCAGAATGCAAC---GGTTGCACTG
gi|24657104|ref|NM_057411.3|        GCAGCGCGAGTGCTTCCAGTGCCACCCAGAATGCAAC---GGTTGCACTG
gi|302179500|gb|HM749883.1|         -------------------TGCCACCCCAACTGCACCTATGGCTGTGAAG
gi|47522839|ref|NM_214007.1|        -------------------TGCCACCCCAACTGCACCTACGGCTGTGTCG
gi|41327737|ref|NM_005228.3|        -------------------TGCCATCCAAACTGCACCTACGGATGCACTG
gi|6478867|gb|M37394.2|RATEGFR      -------------------TGCCATGCAAACTGTACCTATGGATGTGCTG

gi|24657088|ref|NM_057410.3|        GTCCGGGTGCCGACGATTGCAAGTCTTGTCGCAACTTCAAGTTGTTCGAC
gi|24657104|ref|NM_057411.3|        GTCCGGGTGCCGACGATTGCAAGTCTTGTCGCAACTTCAAGTTGTTCGAC
gi|302179500|gb|HM749883.1|         GGCCAGGTCTCGAAGGCTGTCCACAA---AAAGGGCCCAAGATC------
gi|47522839|ref|NM_214007.1|        GACCAGGTCTCGAGGGCTGTGCGGTG---GACAGGCCCAAGATC------
gi|41327737|ref|NM_005228.3|        GGCCAGGTCTTGAAGGCTGTCCAACG---AATGGGCCTAAGATC------
gi|6478867|gb|M37394.2|RATEGFR      GGCCAGGCCTTAAAGGATGTCAACAACCAGAAGGGCCAAAGATC------

gi|24657088|ref|NM_057410.3|        GCGAATGAGACGGGTCCCTATGTGAACTCCACGATGTTCAATTGCACCTC
gi|24657104|ref|NM_057411.3|        GCGAATGAGACGGGTCCCTATGTGAACTCCACGATGTTCAATTGCACCTC
gi|302179500|gb|HM749883.1|         ---------------CCGTCCATTGCCACGGGCATCGTG-----------
gi|47522839|ref|NM_214007.1|        ---------------CCGTCCATCGCCACCGGGATAGTG-----------
gi|41327737|ref|NM_005228.3|        ---------------CCGTCCATCGCCACTGGGATGGTG-----------
gi|6478867|gb|M37394.2|RATEGFR      ---------------CCATCCATCGCCACTGGGATTGTG-----------

gi|24657088|ref|NM_057410.3|        GAAGTGTCCCTTGGAGATGCGACATGTGAACTATCAGTACACGGCCATTG
gi|24657104|ref|NM_057411.3|        GAAGTGTCCCTTGGAGATGCGACATGTGAACTATCAGTACACGGCCATTG
gi|302179500|gb|HM749883.1|         --------------------------------------------------
gi|47522839|ref|NM_214007.1|        --------------------------------------------------
gi|41327737|ref|NM_005228.3|        --------------------------------------------------
gi|6478867|gb|M37394.2|RATEGFR      --------------------------------------------------

gi|24657088|ref|NM_057410.3|        GACCCTACTGTGCAGCTAGTCCGCCGAGGAGCAGCAAGATAACTGCCAAT
gi|24657104|ref|NM_057411.3|        GACCCTACTGTGCAGCTAGTCCGCCGAGGAGCAGCAAGATAACTGCCAAT
gi|302179500|gb|HM749883.1|         --------------------------------------------------
gi|47522839|ref|NM_214007.1|        --------------------------------------------------
gi|41327737|ref|NM_005228.3|        --------------------------------------------------
gi|6478867|gb|M37394.2|RATEGFR      --------------------------------------------------

gi|24657088|ref|NM_057410.3|        CTGGATGTGAACATGATCTTCATTATCACTGGTGCTGTTCTGGTGCCGAC
gi|24657104|ref|NM_057411.3|        CTGGATGTGAACATGATCTTCATTATCACTGGTGCTGTTCTGGTGCCGAC
gi|302179500|gb|HM749883.1|         ------------------------------GGCGGCCTGCTGCTGGTGGT
gi|47522839|ref|NM_214007.1|        ------------------------------GGGGGCCTGCTTCTGGCCGT
gi|41327737|ref|NM_005228.3|        ------------------------------GGGGCCCTCCTCTTGCTGCT
gi|6478867|gb|M37394.2|RATEGFR      ------------------------------GGTGGCCTCCTCTTCATAGT

gi|24657088|ref|NM_057410.3|        GATCTGCATCCTCTGCGTGGTCACATACATTTGTCGGCAAAAGCAAAAGG
gi|24657104|ref|NM_057411.3|        GATCTGCATCCTCTGCGTGGTCACATACATTTGTCGGCAAAAGCAAAAGG
gi|302179500|gb|HM749883.1|         GGTGCTGGCCCTGAGCGTCGGCCTCTTCATG---CGCAGGCGCCACATCG
gi|47522839|ref|NM_214007.1|        GGTGCTGGCCCTGGGGGTCGGCCTCTTTCTG---CGCAGGCGCCACATCG
gi|41327737|ref|NM_005228.3|        GGTGGTGGCCCTGGGGATCGGCCTCTTCATG---CGAAGGCGCCACATCG
gi|6478867|gb|M37394.2|RATEGFR      AGTGGTGGCCCTTGGGATCGGCCTCTTCATG---CGTCGACGTCAGCTTG

gi|24657088|ref|NM_057410.3|        CCAAGAAAGAAACAGTGAAGATGACCATGGCTCTGTCCGGCTGTGAGGAT
gi|24657104|ref|NM_057411.3|        CCAAGAAAGAAACAGTGAAGATGACCATGGCTCTGTCCGGCTGTGAGGAT
gi|302179500|gb|HM749883.1|         TGCGCAAGCGCACACTGCGCCGGCTG------CTGCAGGAGCGTGAGCTC
gi|47522839|ref|NM_214007.1|        TCCGCAAGCGCACGCTGCGCCGGCTG------CTGCAGGAGCGGGAGCTG
gi|41327737|ref|NM_005228.3|        TTCGGAAGCGCACGCTGCGGAGGCTG------CTGCAGGAGAGGGAGCTT
gi|6478867|gb|M37394.2|RATEGFR      TCCGAAAACGTACACTACGCCGCCTG------CTTCAAGAGAGAGAGCTC

gi|24657088|ref|NM_057410.3|        TCCGAGCCGCTGCGTCCCTCGAACATTGGAGCCAATCTATGCAAGTTGCG
gi|24657104|ref|NM_057411.3|        TCCGAGCCGCTGCGTCCCTCGAACATTGGAGCCAATCTATGCAAGTTGCG
gi|302179500|gb|HM749883.1|         GTGGAGCCTCTGACGCCCAGCGGAGAAGCTCCCAACCAAGCTCTCTTGAG
gi|47522839|ref|NM_214007.1|        GTTGAGCCTCTCACACCCAGTGGAGAAGCTCCCAACCAAGCTCTCTTGAG
gi|41327737|ref|NM_005228.3|        GTGGAGCCTCTTACACCCAGTGGAGAAGCTCCCAACCAAGCTCTCTTGAG
gi|6478867|gb|M37394.2|RATEGFR      GTGGAACCTCTCACACCCAGCGGAGAAGCTCCGAACCAAGCCCACTTGAG

gi|24657088|ref|NM_057410.3|        CATTGTCAAGGACGCCGAGTTGCGCAAGGGCGGAGTCCTCGGAATGGGAG
gi|24657104|ref|NM_057411.3|        CATTGTCAAGGACGCCGAGTTGCGCAAGGGCGGAGTCCTCGGAATGGGAG
gi|302179500|gb|HM749883.1|         GATCCTAAAGGAAACAGAATTCAAGAAGGTCAAGGTGCTGGGCTCGGGAG
gi|47522839|ref|NM_214007.1|        GATCCTGAAGGAGACGGAATTCAAAAAGGTCAAGGTGCTGGGCTCCGGCG
gi|41327737|ref|NM_005228.3|        GATCTTGAAGGAAACTGAATTCAAAAAGATCAAAGTGCTGGGCTCCGGTG
gi|6478867|gb|M37394.2|RATEGFR      GATATTAAAGGAAACAGAATTCAAAAAGATCAAAGTTCTGGGTTCAGGAG

gi|24657088|ref|NM_057410.3|        CCTTTGGACGAGTGTACAAGGGCGTTTGGGTGCCGGAGGGTGAGAACGTC
gi|24657104|ref|NM_057411.3|        CCTTTGGACGAGTGTACAAGGGCGTTTGGGTGCCGGAGGGTGAGAACGTC
gi|302179500|gb|HM749883.1|         CATTTGGCACCGTGTACAAGGGACTCTGGATCCCAGAAGGCGAGAAGGTT
gi|47522839|ref|NM_214007.1|        CGTTCGGCACGGTGTACAAGGGCCTCTGGATCCCAGAAGGTGAGAAGGTG
gi|41327737|ref|NM_005228.3|        CGTTCGGCACGGTGTATAAGGGACTCTGGATCCCAGAAGGTGAGAAAGTT
gi|6478867|gb|M37394.2|RATEGFR      CATTTGGCACAGTGTATAAGGGTCTCTGGATCCCAGAAGGCGAGAAAGTG

gi|24657088|ref|NM_057410.3|        AAGATTCCAGTGGCCATTAAGGAGCTGCTCAAGTCCACAGGCGCCGAGTC
gi|24657104|ref|NM_057411.3|        AAGATTCCAGTGGCCATTAAGGAGCTGCTCAAGTCCACAGGCGCCGAGTC
gi|302179500|gb|HM749883.1|         AAAATTCCTGTAGCTATCAAGGAATTAAGAGAAGCCACATCTCCAAAAGC
gi|47522839|ref|NM_214007.1|        AAAATTCCTGTGGCTATCAAGGAATTAAGAGAAGCCACTTCTCCAAAAGC
gi|41327737|ref|NM_005228.3|        AAAATTCCCGTCGCTATCAAGGAATTAAGAGAAGCAACATCTCCGAAAGC
gi|6478867|gb|M37394.2|RATEGFR      AAAATCCCTGTGGCCATCAAGGAGTTAAGAGAAGCCACATCTCCCAAAGC

gi|24657088|ref|NM_057410.3|        AAGCGAAGAGTTCCTCCGCGAAGCCTACATCATGGCCTCTGTGGAGCACG
gi|24657104|ref|NM_057411.3|        AAGCGAAGAGTTCCTCCGCGAAGCCTACATCATGGCCTCTGTGGAGCACG
gi|302179500|gb|HM749883.1|         CAACAAGGAAATTCTTGATGAGGCCTACGTGATGGCCAGTGTGGACAACC
gi|47522839|ref|NM_214007.1|        CAACAAGGAAATTCTTGACGAAGCCTACGTGATGGCCAGTGTGGACAATC
gi|41327737|ref|NM_005228.3|        CAACAAGGAAATCCTCGATGAAGCCTACGTGATGGCCAGCGTGGACAACC
gi|6478867|gb|M37394.2|RATEGFR      CAACAAGGAAATCCTTGATGAAGCCTACGTGATGGCCAGTGTGGACAACC

gi|24657088|ref|NM_057410.3|        TTAATCTGCTGAAGCTCCTGGCCGTCTGCATGTCCTCACAAATGATGCTA
gi|24657104|ref|NM_057411.3|        TTAATCTGCTGAAGCTCCTGGCCGTCTGCATGTCCTCACAAATGATGCTA
gi|302179500|gb|HM749883.1|         CCCATGTGTGCCGCCTCCTGGGCATCTGCCTGACCTCCACCGTGCAGCTC
gi|47522839|ref|NM_214007.1|        CTCATGTGTGCCGCCTCCTGGGCATCTGCCTGACCTCCACGGTGCAGCTC
gi|41327737|ref|NM_005228.3|        CCCACGTGTGCCGCCTGCTGGGCATCTGCCTCACCTCCACCGTGCAGCTC
gi|6478867|gb|M37394.2|RATEGFR      CTCATGTATGCCGCCTCCTGGGCATCTGTCTGACCTCCACTGTCCAGCTC

gi|24657088|ref|NM_057410.3|        ATCACGCAACTGATGCCGCTTGGCTGCCTGTTGGACTATGTGCGAAATAA
gi|24657104|ref|NM_057411.3|        ATCACGCAACTGATGCCGCTTGGCTGCCTGTTGGACTATGTGCGAAATAA
gi|302179500|gb|HM749883.1|         ATCACACAGCTCATGCCCTTCGGCTGCCTGCTGGACTACGTCCGCGAGCA
gi|47522839|ref|NM_214007.1|        ATCACGCAGCTCATGCCCTTCGGCTGCCTCCTGGACTACGTCCGCGAGCA
gi|41327737|ref|NM_005228.3|        ATCACGCAGCTCATGCCCTTCGGCTGCCTCCTGGACTATGTCCGGGAACA
gi|6478867|gb|M37394.2|RATEGFR      ATTACACAACTCATGCCCTATGGTTGCCTCCTGGACTATGTCCGAGAACA

gi|24657088|ref|NM_057410.3|        CCGGGACAAGATCGGCTCTAAGGCTCTGCTCAACTGGAGCACGCAAATCG
gi|24657104|ref|NM_057411.3|        CCGGGACAAGATCGGCTCTAAGGCTCTGCTCAACTGGAGCACGCAAATCG
gi|302179500|gb|HM749883.1|         CAAGGACAATGTCGGCTCCCAGTACCTGCTCAACTGGTGTGTGCAGATCG
gi|47522839|ref|NM_214007.1|        CAAGGACAACATCGGCTCCCAGCACCTGCTCAACTGGTGTGTGCAGATCG
gi|41327737|ref|NM_005228.3|        CAAAGACAATATTGGCTCCCAGTACCTGCTCAACTGGTGTGTGCAGATCG
gi|6478867|gb|M37394.2|RATEGFR      TAAGGACAACATTGGCTCCCAGTACCTACTCAACTGGTGTGTGCAGATTG

gi|24657088|ref|NM_057410.3|        CCAAGGGCATGTCGTATCTGGAGGAGAAGCGACTGGTCCACAGAGACTTG
gi|24657104|ref|NM_057411.3|        CCAAGGGCATGTCGTATCTGGAGGAGAAGCGACTGGTCCACAGAGACTTG
gi|302179500|gb|HM749883.1|         CAAAGGGCATGAATTACCTGGAAGACCGGCGCTTGGTGCATAGGGACCTG
gi|47522839|ref|NM_214007.1|        CAAAGGGCATGAACTATCTGGAAGACCGGCGCTTGGTGCACCGAGACCTG
gi|41327737|ref|NM_005228.3|        CAAAGGGCATGAACTACTTGGAGGACCGTCGCTTGGTGCACCGCGACCTG
gi|6478867|gb|M37394.2|RATEGFR      CAAAGGGCATGAACTACCTGGAAGACCGGCGTTTGGTACACCGTGACTTG

gi|24657088|ref|NM_057410.3|        GCTGCCCGCAATGTCCTGGTGCAGACTCCCTCGCTGGTGAAGATCACCGA
gi|24657104|ref|NM_057411.3|        GCTGCCCGCAATGTCCTGGTGCAGACTCCCTCGCTGGTGAAGATCACCGA
gi|302179500|gb|HM749883.1|         GCAGCCAGGAACGTGCTGGTGAAGACGCCGCAGCACGTGAAGATCACAGA
gi|47522839|ref|NM_214007.1|        GCGGCCAGGAATGTGCTGGTGAAGACACCGCAGCATGTCAAGATCACTGA
gi|41327737|ref|NM_005228.3|        GCAGCCAGGAACGTACTGGTGAAAACACCGCAGCATGTCAAGATCACAGA
gi|6478867|gb|M37394.2|RATEGFR      GCAGCCAGGAATGTACTGGTAAAGACACCACAGCATGTCAAGATCACAGA

gi|24657088|ref|NM_057410.3|        CTTTGGGCTGGCCAAGTTGCTGAGCAGCGATTCCAATGAGTACAAGGCTG
gi|24657104|ref|NM_057411.3|        CTTTGGGCTGGCCAAGTTGCTGAGCAGCGATTCCAATGAGTACAAGGCTG
gi|302179500|gb|HM749883.1|         CTTCGGGCTGGCCAAGCTGCTGGGTGCCGAGGAGAAGGAGTATCATGCAG
gi|47522839|ref|NM_214007.1|        CTTTGGGCTGGCCAAGCTGCTGGGCGCCGAGGAGAAAGAGTACCACGCGG
gi|41327737|ref|NM_005228.3|        TTTTGGGCTGGCCAAACTGCTGGGTGCGGAAGAGAAAGAATACCATGCAG
gi|6478867|gb|M37394.2|RATEGFR      TTTTGGACTGGCCAAACTGCTTGGTGCTGAGGAGAAAGAATACCATGCAG

gi|24657088|ref|NM_057410.3|        CTGGCGGCAAGATGCCCATCAAGTGGTTGGCACTGGAGTGCATTCGCAAT
gi|24657104|ref|NM_057411.3|        CTGGCGGCAAGATGCCCATCAAGTGGTTGGCACTGGAGTGCATTCGCAAT
gi|302179500|gb|HM749883.1|         AAGGAGGCAAGGTCCCTATCAAATGGATGGCTTTGGAATCAATTTTACAC
gi|47522839|ref|NM_214007.1|        AAGGAGGCAAAGTGCCCATCAAGTGGCTGGCTCTAGAGTCAATCCTGCAC
gi|41327737|ref|NM_005228.3|        AAGGAGGCAAAGTGCCTATCAAGTGGATGGCATTGGAATCAATTTTACAC
gi|6478867|gb|M37394.2|RATEGFR      AGGGGGGCAAAGTGCCTATCAAGTGGATGGCTTTGGAATCAATTTTACAC

gi|24657088|ref|NM_057410.3|        CGTGTATTCACCAGCAAGTCCGATGTCTGGGCCTTTGGTGTGACAATTTG
gi|24657104|ref|NM_057411.3|        CGTGTATTCACCAGCAAGTCCGATGTCTGGGCCTTTGGTGTGACAATTTG
gi|302179500|gb|HM749883.1|         CGAATTTATACCCATCAGAGTGATGTCTGGAGCTATGGAGTCACTGTTTG
gi|47522839|ref|NM_214007.1|        CGTGTATACACCCACCAGAGTGACGTCTGGAGCTACGGAGTCACCGTTTG
gi|41327737|ref|NM_005228.3|        AGAATCTATACCCACCAGAGTGATGTCTGGAGCTACGGGGTGACCGTTTG
gi|6478867|gb|M37394.2|RATEGFR      CGAATTTATACACACCAAAGCGACGTCTGGAGCTATGGAGTCACCGTGTG

gi|24657088|ref|NM_057410.3|        GGAACTGCTGACCTTTGGCCAGCGTCCACACGAGAACATCCCCGCTAAGG
gi|24657104|ref|NM_057411.3|        GGAACTGCTGACCTTTGGCCAGCGTCCACACGAGAACATCCCCGCTAAGG
gi|302179500|gb|HM749883.1|         GGAGTTGATGACCTTTGGATCCAAGCCTTACGATGGAATCCCTGCGAGTG
gi|47522839|ref|NM_214007.1|        GGAGCTGATGACCTTTGGGTCCAAGCCTTATGACGGGATCCCCGCGAGTG
gi|41327737|ref|NM_005228.3|        GGAGTTGATGACCTTTGGATCCAAGCCATATGACGGAATCCCTGCCAGCG
gi|6478867|gb|M37394.2|RATEGFR      GGAACTGATGACCTTTGGGTCCAAGCCTTATGATGGGATCCCTGCAAGTG

gi|24657088|ref|NM_057410.3|        ATATTCCCGATCTTATTGAAGTCGGTCTGAAGCTGGAGCAGCCGGAGATT
gi|24657104|ref|NM_057411.3|        ATATTCCCGATCTTATTGAAGTCGGTCTGAAGCTGGAGCAGCCGGAGATT
gi|302179500|gb|HM749883.1|         AGATCTCGACTGTCCTGGAGAAAGGAGAGCGCCTCCCACAGCCACCCATC
gi|47522839|ref|NM_214007.1|        AGATCTCGACCGTCCTGGAGAAGGGAGAGCGCCTCCCGCAGCCCCCCATC
gi|41327737|ref|NM_005228.3|        AGATCTCCTCCATCCTGGAGAAAGGAGAACGCCTCCCTCAGCCACCCATA
gi|6478867|gb|M37394.2|RATEGFR      AGATCTCATCCATCCTAGAGAAAGGAGAGCGCCTTCCACAGCCACCTATC

gi|24657088|ref|NM_057410.3|        TGTTCGCTGGACATTTACTGCACACTTCTCTCGTGCTGGCACTTGGATGC
gi|24657104|ref|NM_057411.3|        TGTTCGCTGGACATTTACTGCACACTTCTCTCGTGCTGGCACTTGGATGC
gi|302179500|gb|HM749883.1|         TGCACCATCGACGTCTACATGATCATGGTCAAGTGCTGGATGATAGATGC
gi|47522839|ref|NM_214007.1|        TGCACCATTGATGTCTACATGATCATGGTCAAGTGCTGGATGATAGATGC
gi|41327737|ref|NM_005228.3|        TGTACCATCGATGTCTACATGATCATGGTCAAGTGCTGGATGATAGACGC
gi|6478867|gb|M37394.2|RATEGFR      TGCACCATCGACGTCTACATGATCATGGTCAAGTGCTGGATGATAGATGC

gi|24657088|ref|NM_057410.3|        CGCCATGCGTCCAACCTTCAAGCAGCTGACTACGGTCTTTGCTGAGTTCG
gi|24657104|ref|NM_057411.3|        CGCCATGCGTCCAACCTTCAAGCAGCTGACTACGGTCTTTGCTGAGTTCG
gi|302179500|gb|HM749883.1|         AGACAGTCGCCCAAAGTTCCGTGAGTTGATCCTTGAATTCTCCAAGATGG
gi|47522839|ref|NM_214007.1|        TGATAGTCGCCCAAAGTTCCGTGAGCTGATCATCGAATTCTCCAAAATGG
gi|41327737|ref|NM_005228.3|        AGATAGTCGCCCAAAGTTCCGTGAGTTGATCATCGAATTCTCCAAAATGG
gi|6478867|gb|M37394.2|RATEGFR      TGATAGCCGCCCAAAGTTCCGAGAGTTGATTCTCGAATTCTCCAAAATGG

gi|24657088|ref|NM_057410.3|        CCAGAGATCCGGGTCGCTATCTGGCCATTCCCGGGGATAAGTTCACCCGG
gi|24657104|ref|NM_057411.3|        CCAGAGATCCGGGTCGCTATCTGGCCATTCCCGGGGATAAGTTCACCCGG
gi|302179500|gb|HM749883.1|         CCCGAGACCCCCAGCGCTACCTTGTCATCCAGGGGGACGAGAGAATGCAT
gi|47522839|ref|NM_214007.1|        CCCGAGACCCCCAGCGCTACCTTGTCATCCAGGGAGACGAGCGAATGCAC
gi|41327737|ref|NM_005228.3|        CCCGAGACCCCCAGCGCTACCTTGTCATTCAGGGGGATGAAAGAATGCAT
gi|6478867|gb|M37394.2|RATEGFR      CCAGAGACCCACAGCGCTACCTTGTTATCCAGGGGGATGAAAGGATGCAT

gi|24657088|ref|NM_057410.3|        CTGCCGGCCTACACGAGTCAGGATGAGAAGGATCTCATCCGAAAATTGGC
gi|24657104|ref|NM_057411.3|        CTGCCGGCCTACACGAGTCAGGATGAGAAGGATCTCATCCGAAAATTGGC
gi|302179500|gb|HM749883.1|         TTGCCA---------AGCCCTACGGACTCCAACTTCTACCGCGCCCTGAT
gi|47522839|ref|NM_214007.1|        TTGCCA---------AGCCCTACGGACTCCAACTTCTACCGCGCCCTGAT
gi|41327737|ref|NM_005228.3|        TTGCCA---------AGTCCTACAGACTCCAACTTCTACCGTGCCCTGAT
gi|6478867|gb|M37394.2|RATEGFR      TTGCCG---------AGCCCTACAGACTCCAACTTTTACCGAGCCCTGAT

gi|24657088|ref|NM_057410.3|        TCCCACCACCGATGGGTCCGAAGCCATTGCGGAACCCGATGACTACCTGC
gi|24657104|ref|NM_057411.3|        TCCCACCACCGATGGGTCCGAAGCCATTGCGGAACCCGATGACTACCTGC
gi|302179500|gb|HM749883.1|         GGATGAGGAGGAC---ATGGAGGATGTTGTGGATGCCGATGAGTACCTCG
gi|47522839|ref|NM_214007.1|        GGACGAGGAGGAC---ATGGAGGATGTGGTGGACGCCGACGAGTACCTCG
gi|41327737|ref|NM_005228.3|        GGATGAAGAAGAC---ATGGACGACGTGGTGGATGCCGACGAGTACCTCA
gi|6478867|gb|M37394.2|RATEGFR      GGAGGAGGAGGAC---ATGGAAGACGTAGTTGATGCTGATGAATACCTCA

gi|24657088|ref|NM_057410.3|        AACCCAAGGCAGCACCTGGTCCTAGTCACAGAACCGACTGCACGGATGAG
gi|24657104|ref|NM_057411.3|        AACCCAAGGCAGCACCTGGTCCTAGTCACAGAACCGACTGCACGGATGAG
gi|302179500|gb|HM749883.1|         TCCCCCAGCAGGGCTTCTTCCACAGCCCCACCACCTCCCGGACACCCCTC
gi|47522839|ref|NM_214007.1|        TCCCCCAGCAGGGCTTCTTCCACAGCCCCGCCACCTCCCGGACGCCGCTG
gi|41327737|ref|NM_005228.3|        TCCCACAGCAGGGCTTCTTCAGCAGCCCCTCCACGTCACGGACTCCCCTC
gi|6478867|gb|M37394.2|RATEGFR      TCCCACAGCAAGGCTTCTTCAACAGCCCATCCACGTCACGGACTCCACTC

gi|24657088|ref|NM_057410.3|        ATACCCAAGCTGAACCGCTACTGCAAGGATCCTAGC------------AA
gi|24657104|ref|NM_057411.3|        ATACCCAAGCTGAACCGCTACTGCAAGGATCCTAGC------------AA
gi|302179500|gb|HM749883.1|         CTCAGCTCGCTGAGCACCTCCAGCAACACTCCCACTGTGACTTGCGTTGA
gi|47522839|ref|NM_214007.1|        CTCAGCTCTCTGAGCGCCACCAGCAGCACCCCCGCTGTGGCTTGCGTTGA
gi|41327737|ref|NM_005228.3|        CTGAGCTCTCTGAGTGCAACCAGCAACAATTCCACCGTGGCTTGCATTGA
gi|6478867|gb|M37394.2|RATEGFR      TTGAGCTCTCTGAGTGCAAATAGCAACAGTTCCACTGTGGCTTGCATTAA

gi|24657088|ref|NM_057410.3|        CAAGAATTCG------------AGTACCGGAGACGATGAGACGGATTCGA
gi|24657104|ref|NM_057411.3|        CAAGAATTCG------------AGTACCGGAGACGATGAGACGGATTCGA
gi|302179500|gb|HM749883.1|         TAGAAATGGG------AGCTACCCTCTCAAGGAAGACAGCTTCCTGCAGC
gi|47522839|ref|NM_214007.1|        CAGAAACGGG---CAGAGTTATCCCCTCAAGGAAGACAGCTTCCTGCAGC
gi|41327737|ref|NM_005228.3|        TAGAAATGGGCTGCAAAGCTGTCCCATCAAGGAAGACAGCTTCTTGCAGC
gi|6478867|gb|M37394.2|RATEGFR      TAGAAATGGG------AGCTGCCGTGTCAAAGAAGACGCCTTCTTGCAAC

gi|24657088|ref|NM_057410.3|        GTGCCCGGGAAGTGGGCGTGGGTAATCTGCGCCTCGAT------------
gi|24657104|ref|NM_057411.3|        GTGCCCGGGAAGTGGGCGTGGGTAATCTGCGCCTCGAT------------
gi|302179500|gb|HM749883.1|         GCTACAGCTCAGACCCCACTGGTGCCCTCATCGAGGACAGCATGGACGAC
gi|47522839|ref|NM_214007.1|        GGTACAGCTCCGACCCCACTGGCGCCCTGACCGAGGACAGCCTAGACGAC
gi|41327737|ref|NM_005228.3|        GATACAGCTCAGACCCCACAGGCGCCTTGACTGAGGACAGCATAGACGAC
gi|6478867|gb|M37394.2|RATEGFR      GGTATAGCTCCGATCCCACCAGCGTCCTGACAGAGGACAACATAGATGAC

gi|24657088|ref|NM_057410.3|        ------CTACCAGTCGATGAGGATGATTACCTGATGCCCACATGCCAACC
gi|24657104|ref|NM_057411.3|        ------CTACCAGTCGATGAGGATGATTACCTGATGCCCACATGCCAACC
gi|302179500|gb|HM749883.1|         GCTTTCCTCCCAGTACCCGAA------TATGTAAACCAATCTGTTCCCAA
gi|47522839|ref|NM_214007.1|        ACTTTTCTCCCAGCACCCGAA------TATGTAAACCAGTCTGTTCCCAA
gi|41327737|ref|NM_005228.3|        ACCTTCCTCCCAGTGCCTGAA------TACATAAACCAGTCCGTTCCCAA
gi|6478867|gb|M37394.2|RATEGFR      ACATTCCTTCCCGTGCCTGAA------TATATAAACCAATCTGTTCCCAA

gi|24657088|ref|NM_057410.3|        GGGGCCCAACAACAACAACAACATAAATAATCCC---------AATCAAA
gi|24657104|ref|NM_057411.3|        GGGGCCCAACAACAACAACAACATAAATAATCCC---------AATCAAA
gi|302179500|gb|HM749883.1|         AAGACCC------GCAGGCTCTGTCCAGAACCCTGTCTATCACAATCAGC
gi|47522839|ref|NM_214007.1|        GAGGCCC------GCGGGCTCCGTCCAGAACCCTGTCTACCACAATCAGC
gi|41327737|ref|NM_005228.3|        AAGGCCC------GCTGGCTCTGTGCAGAATCCTGTCTATCACAATCAGC
gi|6478867|gb|M37394.2|RATEGFR      GAGGCCG------GCTGGCTCTGTGCAGAACCCAGTCTATCACAATCAGC

gi|24657088|ref|NM_057410.3|        ACAATATGGCAGCTGTGGGCGTGGCTGCCGGCTACATGGATCTCATCGGA
gi|24657104|ref|NM_057411.3|        ACAATATGGCAGCTGTGGGCGTGGCTGCCGGCTACATGGATCTCATCGGA
gi|302179500|gb|HM749883.1|         CTCTATATCCAGCTCCTGGCAGAGACCCTCAGTACCAAAAT------TCA
gi|47522839|ref|NM_214007.1|        CTCTCAGTGCAGCTCCTGGCCGGGACCCCCACTACCAGAAC------TCC
gi|41327737|ref|NM_005228.3|        CTCTGAACCCCGCGCCCAGCAGAGACCCACACTACCAGGAC------CCC
gi|6478867|gb|M37394.2|RATEGFR      CCCTGCATCCAGCTCCTGGAAGAGACCTGCATTATCAAAAT------CCC

gi|24657088|ref|NM_057410.3|        GTGCCCGTTAGTGTGGACAATCCGGAGTATCTGCTAAACGCGCAGACACT
gi|24657104|ref|NM_057411.3|        GTGCCCGTTAGTGTGGACAATCCGGAGTATCTGCTAAACGCGCAGACACT
gi|302179500|gb|HM749883.1|         CTCAGCAACGCCGTGGACAACCCTGAGTATCTCAACACCACCCATCCTGC
gi|47522839|ref|NM_214007.1|        CACAGCAATGCCGTGGGCAACCCTGAGTATCTCAACACCCCCCGCCCCGC
gi|41327737|ref|NM_005228.3|        CACAGCACTGCAGTGGGCAACCCCGAGTATCTCAACACTGTCCAGCCCAC
gi|6478867|gb|M37394.2|RATEGFR      CATAGCAATGCGGTGAGCAACCCTGAGTATCTCAACACTGCCCAGCCGAC

gi|24657088|ref|NM_057410.3|        GGGTGTTGGG---------GAGTCGCCGATACCCACCCAGACCATCGGGA
gi|24657104|ref|NM_057411.3|        GGGTGTTGGG---------GAGTCGCCGATACCCACCCAGACCATCGGGA
gi|302179500|gb|HM749883.1|         CTGTATCAATGGTGTGCTCGACGGCCCTGCCCTCTGGGCTCAGAAGGGCA
gi|47522839|ref|NM_214007.1|        CTGCATCAACGGAGGACTGGACGGCCCTGCCTTCTGGGCACAGACAGGCA
gi|41327737|ref|NM_005228.3|        CTGTGTCAACAGCACATTCGACAGCCCTGCCCACTGGGCCCAGAAAGGCA
gi|6478867|gb|M37394.2|RATEGFR      CTGCCTCAGTAGTGGGTTTGACAGCTCTGCCCTCTGGATCCAGAAAGGCA

gi|24657088|ref|NM_057410.3|        TACCGGTGATGGGAGTCCCGGGCACCATGGAGGTCAAGGTGCCAATGCCA
gi|24657104|ref|NM_057411.3|        TACCGGTGATGGGAGTCCCGGGCACCATGGAGGTCAAGGTGCCAATGCCA
gi|302179500|gb|HM749883.1|         GTCACCAATTTAGCCTAGACAACCCTGACTACCAGCAGGCCTTCTTTCCC
gi|47522839|ref|NM_214007.1|        GCCACCAGATTAATCTGGACAACCCAGACTACCAGCAGGCCTTCTTCCCC
gi|41327737|ref|NM_005228.3|        GCCACCAAATTAGCCTGGACAACCCTGACTACCAGCAGGACTTCTTTCCC
gi|6478867|gb|M37394.2|RATEGFR      GCCACCAAATGAGCCTGGACAACCCTGACTACCAGCAGGACTTCTTTCCC

gi|24657088|ref|NM_057410.3|        GGCAGTGAGCCAACGAGCTCCGATCACGAGTACTACAATGATACCCAACG
gi|24657104|ref|NM_057411.3|        GGCAGTGAGCCAACGAGCTCCGATCACGAGTACTACAATGATACCCAACG
gi|302179500|gb|HM749883.1|         AAGGAAGCCAAGTCGAATGGCATCTTTAAGGGGCCTGCAGCTGAAAATGC
gi|47522839|ref|NM_214007.1|        AAGGAAGCCAAGTCAAACGGCATCTGTAAGGGTCCCGCCGCCGAAAACGC
gi|41327737|ref|NM_005228.3|        AAGGAAGCCAAGCCAAATGGCATCTTTAAGGGCTCCACAGCTGAAAATGC
gi|6478867|gb|M37394.2|RATEGFR      AAAGAAGCCAAGCCGAATGGCATCTTTAAGGGCCCCACAGCTGAAAATGC

gi|24657088|ref|NM_057410.3|        GGAGTTGCAGCCACTGCATCGAAACCGCAACACGGAGACGAGGGTG
gi|24657104|ref|NM_057411.3|        GGAGTTGCAGCCACTGCATCGAAACCGCAACACGGAGACGAGGGTG
gi|302179500|gb|HM749883.1|         AGAATACCTGCGGGCAGCACCAGCAGGCAGTGACTTTACTGGGGCC
gi|47522839|ref|NM_214007.1|        AGAGTACCTAAGGGCGGCACCAGCCAGCAGTGACCTTACTGGGGCA
gi|41327737|ref|NM_005228.3|        AGAATACCTAAGGGTCGCGCCACAAAGCAGTGAATTTATTGGAGCA
gi|6478867|gb|M37394.2|RATEGFR      AGAGTACCTGCGGGTGGCACCGCCAAGCAGTGAGTTTAGTGGAGCA


""",
        )

        pairwise_alignment = alignment[:2]
        dN, dS = calculate_dn_ds(pairwise_alignment, method="NG86")
        self.assertAlmostEqual(dN, 0.0209, places=4)
        self.assertAlmostEqual(dS, 0.0178, places=4)
        dN, dS = calculate_dn_ds(pairwise_alignment, method="LWL85")
        self.assertAlmostEqual(dN, 0.0203, places=4)
        self.assertAlmostEqual(dS, 0.0164, places=4)

        try:
            import scipy
        except ImportError:
            # Silently skip the rest of the test
            return

        # This should be present:
        from scipy.linalg import expm

        dN, dS = calculate_dn_ds(pairwise_alignment, method="YN00")
        self.assertAlmostEqual(dN, 0.0198, places=4)
        self.assertAlmostEqual(dS, 0.0222, places=4)

        try:
            # New in scipy v0.11
            from scipy.optimize import minimize

            dN, dS = calculate_dn_ds(pairwise_alignment, method="ML")
            self.assertAlmostEqual(dN, 0.0194, places=4)
            self.assertAlmostEqual(dS, 0.0217, places=4)
        except ImportError:
            pass

        # NG86 method with default codon table
        dn_correct = [
            0,
            0.02090783050583131,
            0,
            0.6115239249238438,
            0.6102203266798018,
            0,
            0.6140350835631757,
            0.6040168621204747,
            0.041180350405913294,
            0,
            0.6141532531400524,
            0.6018263135601294,
            0.06701051445629494,
            0.061470360954086874,
            0,
            0.6187088340904762,
            0.6068687248870475,
            0.07386903034833081,
            0.07357890927918581,
            0.05179847072570129,
            0,
        ]
        ds_correct = [
            0,
            0.01783718763890243,
            0,
            2.9382055377913687,
            3.0375115405379267,
            0,
            2.008913071877126,
            2.0182088023715616,
            0.5638033197005285,
            0,
            2.771425931736778,
            2.7353083173058295,
            0.6374483799734671,
            0.723542095485497,
            0,
            -1,
            -1,
            0.953865978141643,
            1.182154857347706,
            0.843182957978177,
            0,
        ]
        dn, ds = calculate_dn_ds_matrix(alignment)
        dn_list = []
        for i in dn.matrix:
            dn_list.extend(i)
        for dn_cal, dn_corr in zip(dn_list, dn_correct):
            self.assertAlmostEqual(dn_cal, dn_corr, places=4)
        ds_list = []
        for i in ds.matrix:
            ds_list.extend(i)
        for ds_cal, ds_corr in zip(ds_list, ds_correct):
            self.assertAlmostEqual(ds_cal, ds_corr, places=4)
        # YN00 method with user specified codon table
        dn_correct = [
            0,
            0.019701773284646867,
            0,
            0.6109649819852769,
            0.6099903856901369,
            0,
            0.6114499930666559,
            0.6028068208599121,
            0.045158286242251426,
            0,
            0.6151835071687592,
            0.6053227393422296,
            0.07034397741651377,
            0.06956967795096626,
            0,
            0.6103850655769698,
            0.5988716898831496,
            0.07905930042150053,
            0.08203052937107111,
            0.05659346894088538,
            0,
        ]
        ds_correct = [
            0,
            0.01881718550096053,
            0,
            1.814457265482046,
            1.8417575124882066,
            0,
            1.5627041719628896,
            1.563930819079887,
            0.4748890153032888,
            0,
            1.6754828466084355,
            1.6531212012501901,
            0.5130923627791538,
            0.5599667707191436,
            0,
            2.0796114236540943,
            2.1452591651827304,
            0.7243066372971764,
            0.8536617406770075,
            0.6509203399899367,
            0,
        ]
        dn, ds = calculate_dn_ds_matrix(
            alignment, method="LWL85", codon_table=CodonTable.unambiguous_dna_by_id[3]
        )
        dn_list = []
        for i in dn.matrix:
            dn_list.extend(i)
        for dn_cal, dn_corr in zip(dn_list, dn_correct):
            self.assertAlmostEqual(dn_cal, dn_corr, places=4)
        ds_list = []
        for i in ds.matrix:
            ds_list.extend(i)
        for ds_cal, ds_corr in zip(ds_list, ds_correct):
            self.assertAlmostEqual(ds_cal, ds_corr, places=4)


try:
    import numpy
except ImportError:
    numpy = None

if numpy:

    class Test_MK(unittest.TestCase):
        def test_mk(self):
            aligner = CodonAligner()
            nucleotide_records = SeqIO.index("codonalign/drosophila.fasta", "fasta")
            protein_alignment = Align.read("codonalign/adh.aln", "clustal")
            self.assertEqual(len(protein_alignment.sequences), 27)
            codon_alignments = []
            protein_record = protein_alignment.sequences[0]
            nucleotide_record = nucleotide_records[protein_record.id]
            self.assertEqual(nucleotide_record.id, protein_record.id)
            alignments = aligner.align(protein_record, nucleotide_record)
            self.assertEqual(len(alignments), 1)
            alignment = next(alignments)
            codon_alignments.append(alignment)
            self.assertTrue(
                np.array_equal(alignment.coordinates, np.array([[0, 256], [0, 768]]))
            )
            self.assertEqual(
                str(alignment),
                """\
gi|9217|e         0 M  A  F  T  L  T  N  K  N  V  V  F  V  A  G  L  G  G  I  G  
gi|9217|e         0 ATGGCGTTTACCTTGACCAACAAGAACGTGGTTTTCGTGGCCGGTCTGGGAGGCATTGGT

gi|9217|e        20 L  D  T  S  K  E  L  V  K  R  D  L  K  N  L  V  I  L  D  R  
gi|9217|e        60 CTGGACACCAGCAAGGAGCTGGTCAAGCGGGACCTGAAGAACCTGGTGATCCTCGACCGC

gi|9217|e        40 I  E  N  P  A  A  I  A  E  L  K  A  I  N  P  K  V  T  V  T  
gi|9217|e       120 ATTGAGAACCCGGCTGCCATCGCCGAGCTGAAGGCAATCAATCCAAAGGTGACCGTCACC

gi|9217|e        60 F  Y  P  Y  D  V  T  V  P  I  A  E  T  T  K  L  L  K  T  I  
gi|9217|e       180 TTCTACCCCTATGATGTGACCGTGCCCATTGCCGAGACCACCAAGCTGCTGAAGACCATC

gi|9217|e        80 F  A  Q  L  K  T  I  D  V  L  I  N  G  A  G  I  L  D  D  H  
gi|9217|e       240 TTCGCCCAGCTGAAGACCATCGATGTCCTGATCAACGGAGCTGGCATCCTGGACGATCAC

gi|9217|e       100 Q  I  E  R  T  I  A  V  N  Y  T  G  L  V  N  T  T  T  A  I  
gi|9217|e       300 CAGATCGAGCGCACCATCGCCGTCAACTACACCGGCCTGGTGAACACCACGACTGCCATC

gi|9217|e       120 L  D  F  W  D  K  R  K  G  G  P  G  G  I  I  C  N  I  G  S  
gi|9217|e       360 CTGGACTTCTGGGACAAGCGCAAGGGTGGACCCGGTGGTATCATCTGCAACATTGGATCC

gi|9217|e       140 V  T  G  F  N  A  I  Y  Q  V  P  V  Y  S  G  T  K  A  A  V  
gi|9217|e       420 GTGACTGGATTCAACGCCATCTACCAGGTGCCCGTTTACTCCGGCACCAAGGCTGCCGTG

gi|9217|e       160 V  N  F  T  S  S  L  A  K  L  A  P  I  T  G  V  T  A  Y  T  
gi|9217|e       480 GTCAACTTCACCAGCTCCCTGGCGAAACTGGCCCCCATCACCGGCGTGACCGCTTACACC

gi|9217|e       180 V  N  P  G  I  T  R  T  T  L  V  H  K  F  N  S  W  L  D  V  
gi|9217|e       540 GTGAACCCCGGCATCACCCGCACCACCCTGGTGCACAAGTTCAACTCCTGGCTGGATGTT

gi|9217|e       200 E  P  Q  V  A  E  K  L  L  A  H  P  T  Q  P  S  L  A  C  A  
gi|9217|e       600 GAGCCCCAGGTGGCCGAGAAGCTCCTGGCTCACCCCACCCAGCCCTCGTTGGCCTGCGCC

gi|9217|e       220 Q  N  F  V  K  A  I  E  L  N  Q  N  G  A  I  W  K  L  D  L  
gi|9217|e       660 CAGAACTTTGTGAAGGCCATCGAGCTGAACCAGAACGGTGCCATCTGGAAACTGGACTTG

gi|9217|e       240 G  T  L  E  A  I  Q  W  S  K  H  W  D  S  G  I   256
gi|9217|e       720 GGCACCCTGGAGGCCATCCAGTGGTCCAAGCACTGGGACTCCGGCATC 768
""",
            )
            protein_record = protein_alignment.sequences[1]
            nucleotide_record = nucleotide_records[protein_record.id]
            self.assertEqual(nucleotide_record.id, protein_record.id)
            alignments = aligner.align(protein_record, nucleotide_record)
            self.assertEqual(len(alignments), 1)
            alignment = next(alignments)
            codon_alignments.append(alignment)
            self.assertTrue(
                np.array_equal(alignment.coordinates, np.array([[0, 256], [0, 768]]))
            )
            self.assertEqual(
                str(alignment),
                """\
gi|9219|e         0 M  A  F  T  L  T  N  K  N  V  V  F  V  A  G  L  G  G  I  G  
gi|9219|e         0 ATGGCGTTTACCTTGACCAACAAGAACGTGGTTTTCGTGGCCGGTCTGGGAGGCATTGGT

gi|9219|e        20 L  D  T  S  K  E  L  V  K  R  D  L  K  N  L  V  I  L  D  R  
gi|9219|e        60 CTGGACACCAGCAAGGAGCTGGTCAAGCGGGACCTGAAGAACCTGGTGATCCTCGACCGC

gi|9219|e        40 I  E  N  P  A  A  I  A  E  L  K  A  I  N  P  K  V  T  V  T  
gi|9219|e       120 ATTGAGAACCCGGCTGCCATCGCCGAGCTGAAGGCAATCAATCCAAAGGTGACCGTCACC

gi|9219|e        60 F  Y  P  Y  D  V  T  V  P  I  A  E  T  T  K  L  L  K  T  I  
gi|9219|e       180 TTCTACCCCTATGATGTGACCGTGCCCATTGCCGAGACCACCAAGCTGCTGAAGACCATC

gi|9219|e        80 F  A  Q  L  K  T  I  D  V  L  I  N  G  A  G  I  L  D  D  H  
gi|9219|e       240 TTCGCCCAGCTGAAGACCATCGATGTCCTGATCAACGGAGCTGGCATCCTGGACGATCAC

gi|9219|e       100 Q  I  E  R  T  I  A  V  N  Y  T  G  L  V  N  T  T  T  A  I  
gi|9219|e       300 CAGATCGAGCGCACCATCGCCGTCAACTACACCGGCCTGGTGAACACCACGACTGCCATC

gi|9219|e       120 L  D  F  W  D  K  R  K  G  G  P  G  G  I  I  C  N  I  G  S  
gi|9219|e       360 CTGGACTTCTGGGACAAGCGCAAGGGTGGACCCGGTGGTATCATCTGCAACATTGGATCC

gi|9219|e       140 V  T  G  F  N  A  I  Y  Q  V  P  V  Y  S  G  T  K  A  A  V  
gi|9219|e       420 GTGACTGGATTCAACGCCATCTACCAGGTGCCCGTTTACTCCGGCACCAAGGCTGCCGTG

gi|9219|e       160 V  N  F  T  S  S  L  A  K  L  A  P  I  T  G  V  T  A  Y  T  
gi|9219|e       480 GTCAACTTCACCAGCTCCCTGGCGAAACTGGCCCCCATCACCGGCGTGACCGCTTACACC

gi|9219|e       180 V  N  P  G  I  T  R  T  T  L  V  H  K  F  N  S  W  L  D  V  
gi|9219|e       540 GTGAACCCCGGCATCACCCGCACCACCCTGGTGCACAAGTTCAACTCCTGGCTGGATGTT

gi|9219|e       200 E  P  Q  V  A  E  K  L  L  A  H  P  T  Q  P  S  L  A  C  A  
gi|9219|e       600 GAGCCCCAGGTGGCCGAGAAGCTCCTGGCTCACCCCACCCAGCCCTCGTTGGCCTGCGCC

gi|9219|e       220 Q  N  F  V  K  A  I  E  L  N  Q  N  G  A  I  W  K  L  D  L  
gi|9219|e       660 CAGAACTTTGTGAAGGCCATCGAGCTGAACCAGAACGGTGCCATCTGGAAACTGGACTTG

gi|9219|e       240 G  T  L  E  A  I  Q  W  S  K  H  W  D  S  G  I   256
gi|9219|e       720 GGCACCCTGGAGGCCATCCAGTGGTCCAAGCACTGGGACTCCGGCATC 768
""",
            )
            protein_record = protein_alignment.sequences[2]
            nucleotide_record = nucleotide_records[protein_record.id]
            self.assertEqual(nucleotide_record.id, protein_record.id)
            alignments = aligner.align(protein_record, nucleotide_record)
            self.assertEqual(len(alignments), 1)
            alignment = next(alignments)
            codon_alignments.append(alignment)
            self.assertTrue(
                np.array_equal(alignment.coordinates, np.array([[0, 256], [0, 768]]))
            )
            self.assertEqual(
                str(alignment),
                """\
gi|9221|e         0 M  A  F  T  L  T  N  K  N  V  V  F  V  A  G  L  G  G  I  G  
gi|9221|e         0 ATGGCGTTTACCTTGACCAACAAGAACGTGGTTTTCGTGGCCGGTCTGGGAGGCATTGGT

gi|9221|e        20 L  D  T  S  K  E  L  V  K  R  D  L  K  N  L  V  I  L  D  R  
gi|9221|e        60 CTGGACACCAGCAAGGAGCTGGTCAAGCGCGACCTGAAGAACCTGGTGATCCTCGACCGC

gi|9221|e        40 I  E  N  P  A  A  I  A  E  L  K  A  I  N  P  K  V  T  V  T  
gi|9221|e       120 ATTGAGAACCCGGCTGCCATCGCCGAGCTGAAGGCAATCAATCCAAAGGTGACCGTCACC

gi|9221|e        60 F  Y  P  Y  D  V  T  V  P  I  A  E  T  T  K  L  L  K  T  I  
gi|9221|e       180 TTCTACCCCTATGATGTGACCGTGCCCATTGCCGAGACCACCAAGCTGCTGAAGACCATC

gi|9221|e        80 F  A  Q  L  K  T  I  D  V  L  I  N  G  A  G  I  L  D  D  H  
gi|9221|e       240 TTCGCCCAGCTGAAGACCATCGATGTCCTGATCAACGGAGCTGGCATCCTGGACGATCAC

gi|9221|e       100 Q  I  E  R  T  I  A  V  N  Y  T  G  L  V  N  T  T  T  A  I  
gi|9221|e       300 CAGATCGAGCGCACCATCGCCGTCAACTACACCGGCCTGGTGAACACAACGACGGCCATC

gi|9221|e       120 L  D  F  W  D  K  R  K  G  G  P  G  G  I  I  C  N  I  G  S  
gi|9221|e       360 CTGGACTTCTGGGACAAGCGCAAGGGTGGACCCGGTGGTATAATCTGCAACATTGGATCC

gi|9221|e       140 V  T  G  F  N  A  I  Y  Q  V  P  V  Y  S  G  T  K  A  A  V  
gi|9221|e       420 GTGACTGGATTCAACGCCATCTACCAGGTGCCCGTTTACTCCGGCACCAAGGCTGCCGTG

gi|9221|e       160 V  N  F  T  S  S  L  A  K  L  A  P  I  T  G  V  T  A  Y  T  
gi|9221|e       480 GTCAACTTCACCAGCTCCCTGGCGAAACTGGCCCCCATCACCGGCGTGACCGCTTACACC

gi|9221|e       180 V  N  P  G  I  T  R  T  T  L  V  H  K  F  N  S  W  L  D  V  
gi|9221|e       540 GTGAACCCCGGCATCACCCGCACCACCCTGGTGCACAAGTTCAACTCCTGGCTGGATGTT

gi|9221|e       200 E  P  Q  V  A  E  K  L  L  A  H  P  T  Q  P  S  L  A  C  A  
gi|9221|e       600 GAGCCCCAGGTGGCCGAGAAGCTCCTGGCTCACCCCACCCAGCCCTCGTTGGCCTGCGCC

gi|9221|e       220 Q  N  F  V  K  A  I  E  L  N  Q  N  G  A  I  W  K  L  D  L  
gi|9221|e       660 CAGAACTTTGTCAAGGCCATCGAGCTGAACCAGAACGGTGCCATCTGGAAACTGGACTTG

gi|9221|e       240 G  T  L  E  A  I  Q  W  S  K  H  W  D  S  G  I   256
gi|9221|e       720 GGTACCCTGGAGGCCATCCAGTGGTCCAAGCACTGGGACTCCGGCATC 768
""",
            )
            protein_record = protein_alignment.sequences[3]
            nucleotide_record = nucleotide_records[protein_record.id]
            self.assertEqual(nucleotide_record.id, protein_record.id)
            alignments = aligner.align(protein_record, nucleotide_record)
            self.assertEqual(len(alignments), 1)
            alignment = next(alignments)
            codon_alignments.append(alignment)
            self.assertTrue(
                np.array_equal(alignment.coordinates, np.array([[0, 256], [0, 768]]))
            )
            self.assertEqual(
                str(alignment),
                """\
gi|9223|e         0 M  A  F  T  L  T  N  K  N  V  V  F  V  A  G  L  G  G  I  G  
gi|9223|e         0 ATGGCGTTTACCTTGACCAACAAGAACGTGGTTTTCGTGGCCGGTCTGGGAGGCATTGGT

gi|9223|e        20 L  D  T  S  K  E  L  V  K  R  D  L  K  N  L  V  I  L  D  R  
gi|9223|e        60 CTGGACACCAGCAAGGAGCTGGTCAAGCGGGACCTGAAGAACCTGGTGATCCTCGACCGC

gi|9223|e        40 I  E  N  P  A  A  I  A  E  L  K  A  I  N  P  K  V  T  V  T  
gi|9223|e       120 ATTGAGAACCCGGCTGCCATCGCCGAGCTGAAGGCAATCAATCCAAAGGTGACCGTCACC

gi|9223|e        60 F  Y  P  Y  D  V  T  V  P  I  A  E  T  T  K  L  L  K  T  I  
gi|9223|e       180 TTCTACCCCTATGATGTGACCGTGCCCATTGCCGAGACCACCAAGCTGCTGAAGACCATC

gi|9223|e        80 F  A  Q  L  K  T  I  D  V  L  I  N  G  A  G  I  L  D  D  H  
gi|9223|e       240 TTCGCCCAGCTGAAGACCATCGATGTCCTGATCAACGGAGCTGGCATCCTGGACGATCAC

gi|9223|e       100 Q  I  E  R  T  I  A  V  N  Y  T  G  L  V  N  T  T  T  A  I  
gi|9223|e       300 CAGATCGAGCGCACCATCGCCGTCAACTACACCGGCCTGGTGAACACCACGACTGCCATC

gi|9223|e       120 L  D  F  W  D  K  R  K  G  G  P  G  G  I  I  C  N  I  G  S  
gi|9223|e       360 CTGGACTTCTGGGACAAGCGCAAGGGTGGACCCGGTGGTATCATCTGCAACATTGGATCC

gi|9223|e       140 V  T  G  F  N  A  I  Y  Q  V  P  V  Y  S  G  T  K  A  A  V  
gi|9223|e       420 GTGACTGGATTCAACGCCATCTACCAGGTGCCCGTTTACTCCGGCACCAAGGCTGCCGTG

gi|9223|e       160 V  N  F  T  S  S  L  A  K  L  A  P  I  T  G  V  T  A  Y  T  
gi|9223|e       480 GTCAACTTCACCAGCTCCCTGGCGAAACTGGCCCCCATCACCGGCGTGACCGCTTACACC

gi|9223|e       180 V  N  P  G  I  T  R  T  T  L  V  H  K  F  N  S  W  L  D  V  
gi|9223|e       540 GTGAACCCCGGCATCACCCGCACCACCCTGGTGCACAAGTTCAACTCCTGGCTGGATGTT

gi|9223|e       200 E  P  Q  V  A  E  K  L  L  A  H  P  T  Q  P  S  L  A  C  A  
gi|9223|e       600 GAGCCCCAGGTGGCCGAGAAGCTCCTGGCTCACCCCACCCAGCCCTCGTTGGCCTGCGCC

gi|9223|e       220 Q  N  F  V  K  A  I  E  L  N  Q  N  G  A  I  W  K  L  D  L  
gi|9223|e       660 CAGAACTTTGTGAAGGCCATCGAGCTGAACCAGAACGGTGCCATCTGGAAACTGGACTTG

gi|9223|e       240 G  T  L  E  A  I  Q  W  S  K  H  W  D  S  G  I   256
gi|9223|e       720 GGCACCCTGGAGGCCATCCAGTGGTCCAAGCACTGGGACTCCGGCATC 768
""",
            )
            protein_record = protein_alignment.sequences[4]
            nucleotide_record = nucleotide_records[protein_record.id]
            self.assertEqual(nucleotide_record.id, protein_record.id)
            alignments = aligner.align(protein_record, nucleotide_record)
            self.assertEqual(len(alignments), 1)
            alignment = next(alignments)
            codon_alignments.append(alignment)
            self.assertTrue(
                np.array_equal(alignment.coordinates, np.array([[0, 256], [0, 768]]))
            )
            self.assertEqual(
                str(alignment),
                """\
gi|9225|e         0 M  A  F  T  L  T  N  K  N  V  V  F  V  A  G  L  G  G  I  G  
gi|9225|e         0 ATGGCGTTTACCTTGACCAACAAGAACGTGGTTTTCGTGGCCGGTCTGGGAGGCATTGGT

gi|9225|e        20 L  D  T  S  K  E  L  V  K  R  D  L  K  N  L  V  I  L  D  R  
gi|9225|e        60 CTGGACACCAGCAAGGAGCTGGTCAAGCGGGACCTGAAGAACCTGGTGATCCTCGACCGC

gi|9225|e        40 I  E  N  P  A  A  I  A  E  L  K  A  I  N  P  K  V  T  V  T  
gi|9225|e       120 ATTGAGAACCCGGCTGCCATCGCCGAGCTGAAGGCAATCAATCCAAAGGTGACCGTAACC

gi|9225|e        60 F  Y  P  Y  D  V  T  V  P  I  A  E  T  T  K  L  L  K  T  I  
gi|9225|e       180 TTCTACCCCTATGATGTGACAGTGCCCATTGCCGAGACCACCAAGCTGCTGAAGACCATC

gi|9225|e        80 F  A  Q  L  K  T  I  D  V  L  I  N  G  A  G  I  L  D  D  H  
gi|9225|e       240 TTCGCCCAGCTGAAGACCATCGATGTCCTGATCAACGGAGCTGGCATCTTGGACGATCAC

gi|9225|e       100 Q  I  E  R  T  I  A  V  N  Y  T  G  L  V  N  T  T  T  A  I  
gi|9225|e       300 CAGATCGAGCGCACCATCGCCGTCAACTACACCGGCCTGGTGAACACCACGACGGCCATC

gi|9225|e       120 L  D  F  W  D  K  R  K  G  G  P  G  G  I  I  C  N  I  G  S  
gi|9225|e       360 CTGGACTTCTGGGACAAGCGCAAGGGTGGACCCGGTGGTATCATCTGCAACATTGGATCC

gi|9225|e       140 V  T  G  F  N  A  I  Y  Q  V  P  V  Y  S  G  T  K  A  A  V  
gi|9225|e       420 GTGACTGGATTCAACGCCATCTACCAGGTGCCCGTTTACTCCGGCACCAAGGCTGCCGTG

gi|9225|e       160 V  N  F  T  S  S  L  A  K  L  A  P  I  T  G  V  T  A  Y  T  
gi|9225|e       480 GTCAACTTCACCAGCTCCCTGGCGAAACTGGCCCCCATCACCGGCGTGACCGCTTACACC

gi|9225|e       180 V  N  P  G  I  T  R  T  T  L  V  H  K  F  N  S  W  L  D  V  
gi|9225|e       540 GTGAACCCCGGCATCACCCGCACCACCCTGGTGCACAAGTTCAACTCCTGGCTGGATGTT

gi|9225|e       200 E  P  Q  V  A  E  K  L  L  A  H  P  T  Q  P  S  L  A  C  A  
gi|9225|e       600 GAGCCCCAGGTGGCTGAGAAGCTCCTGGCTCACCCAACCCAGCCCTCGTTGGCCTGCGCC

gi|9225|e       220 Q  N  F  V  K  A  I  E  L  N  Q  N  G  A  I  W  K  L  D  L  
gi|9225|e       660 CAGAACTTTGTCAAGGCCATCGAGCTGAACCAGAACGGTGCCATCTGGAAACTGGACTTG

gi|9225|e       240 G  T  L  E  A  I  Q  W  S  K  H  W  D  S  G  I   256
gi|9225|e       720 GGTACCCTGGAGGCCATCCAGTGGTCCAAGCACTGGGACTCCGGCATC 768
""",
            )
            protein_record = protein_alignment.sequences[5]
            nucleotide_record = nucleotide_records[protein_record.id]
            self.assertEqual(nucleotide_record.id, protein_record.id)
            alignments = aligner.align(protein_record, nucleotide_record)
            self.assertEqual(len(alignments), 1)
            alignment = next(alignments)
            codon_alignments.append(alignment)
            self.assertTrue(
                np.array_equal(alignment.coordinates, np.array([[0, 256], [0, 768]]))
            )
            self.assertEqual(
                str(alignment),
                """\
gi|9227|e         0 M  A  F  T  L  T  N  K  N  V  V  F  V  A  G  L  G  G  I  G  
gi|9227|e         0 ATGGCGTTTACCTTGACCAACAAGAACGTGGTTTTCGTGGCCGGTCTGGGAGGCATTGGT

gi|9227|e        20 L  D  T  S  K  E  L  V  K  R  D  L  K  N  L  V  I  L  D  R  
gi|9227|e        60 CTGGACACCAGCAAGGAGCTGGTCAAGCGAGACCTGAAGAACCTGGTGATCCTCGACCGC

gi|9227|e        40 I  E  N  P  A  A  I  A  E  L  K  A  I  N  P  K  V  T  V  T  
gi|9227|e       120 ATTGAGAACCCGGCTGCCATCGCCGAGCTGAAGGCAATCAATCCAAAGGTGACCGTCACC

gi|9227|e        60 F  Y  P  Y  D  V  T  V  P  I  A  E  T  T  K  L  L  K  T  I  
gi|9227|e       180 TTCTACCCCTATGATGTGACCGTGCCCATTGCCGAGACCACCAAGCTGCTGAAGACCATC

gi|9227|e        80 F  A  Q  L  K  T  I  D  V  L  I  N  G  A  G  I  L  D  D  H  
gi|9227|e       240 TTCGCCCAGCTGAAGACCATCGATGTCCTGATCAACGGAGCTGGCATCCTGGACGATCAC

gi|9227|e       100 Q  I  E  R  T  I  A  V  N  Y  T  G  L  V  N  T  T  T  A  I  
gi|9227|e       300 CAGATCGAGCGCACCATCGCCGTCAACTACACCGGCCTGGTGAACACCACGACGGCCATC

gi|9227|e       120 L  D  F  W  D  K  R  K  G  G  P  G  G  I  I  C  N  I  G  S  
gi|9227|e       360 CTGGACTTCTGGGACAAGCGCAAGGGTGGACCCGGTGGTATCATCTGCAACATTGGATCC

gi|9227|e       140 V  T  G  F  N  A  I  Y  Q  V  P  V  Y  S  G  T  K  A  A  V  
gi|9227|e       420 GTGACTGGATTCAACGCCATCTACCAGGTGCCCGTTTACTCCGGCACCAAGGCTGCCGTG

gi|9227|e       160 V  N  F  T  S  S  L  A  K  L  A  P  I  T  G  V  T  A  Y  T  
gi|9227|e       480 GTCAACTTCACCAGCTCCCTGGCGAAACTGGCCCCCATTACCGGCGTGACCGCTTACACC

gi|9227|e       180 V  N  P  G  I  T  R  T  T  L  V  H  K  F  N  S  W  L  D  V  
gi|9227|e       540 GTGAACCCCGGCATCACCCGCACCACCCTGGTGCACAAGTTCAACTCCTGGCTGGATGTT

gi|9227|e       200 E  P  Q  V  A  E  K  L  L  A  H  P  T  Q  P  S  L  A  C  A  
gi|9227|e       600 GAGCCCCAGGTGGCCGAGAAGCTCCTGGCTCACCCCACCCAGCCCTCGTTGGCCTGCGCC

gi|9227|e       220 Q  N  F  V  K  A  I  E  L  N  Q  N  G  A  I  W  K  L  D  L  
gi|9227|e       660 CAGAACTTTGTGAAGGCCATCGAGCTGAACCAGAACGGTGCCATCTGGAAACTGGACTTG

gi|9227|e       240 G  T  L  E  A  I  Q  W  S  K  H  W  D  S  G  I   256
gi|9227|e       720 GGCACCCTGGAGGCCATCCAGTGGTCCAAGCACTGGGACTCCGGCATC 768
""",
            )
            protein_record = protein_alignment.sequences[6]
            nucleotide_record = nucleotide_records[protein_record.id]
            self.assertEqual(nucleotide_record.id, protein_record.id)
            alignments = aligner.align(protein_record, nucleotide_record)
            self.assertEqual(len(alignments), 1)
            alignment = next(alignments)
            codon_alignments.append(alignment)
            self.assertTrue(
                np.array_equal(alignment.coordinates, np.array([[0, 256], [0, 768]]))
            )
            self.assertEqual(
                str(alignment),
                """\
gi|9229|e         0 M  A  F  T  L  T  N  K  N  V  V  F  V  A  G  L  G  G  I  G  
gi|9229|e         0 ATGGCGTTTACCTTGACCAACAAGAACGTGGTTTTCGTGGCCGGTCTGGGAGGCATTGGT

gi|9229|e        20 L  D  T  S  K  E  L  V  K  R  D  L  K  N  L  V  I  L  D  R  
gi|9229|e        60 CTGGACACCAGCAAGGAGCTGGTCAAGCGGGACCTGAAGAACCTGGTGATCCTCGACCGC

gi|9229|e        40 I  E  N  P  A  A  I  A  E  L  K  A  I  N  P  K  V  T  V  T  
gi|9229|e       120 ATTGAGAACCCGGCTGCCATCGCCGAGCTGAAGGCAATCAATCCAAAGGTGACCGTCACC

gi|9229|e        60 F  Y  P  Y  D  V  T  V  P  I  A  E  T  T  K  L  L  K  T  I  
gi|9229|e       180 TTCTACCCCTATGATGTGACCGTGCCCATTGCCGAGACCACCAAGCTGCTGAAGACCATC

gi|9229|e        80 F  A  Q  L  K  T  I  D  V  L  I  N  G  A  G  I  L  D  D  H  
gi|9229|e       240 TTCGCCCAGCTGAAGACCATCGATGTCCTGATCAACGGAGCTGGCATCCTGGACGATCAC

gi|9229|e       100 Q  I  E  R  T  I  A  V  N  Y  T  G  L  V  N  T  T  T  A  I  
gi|9229|e       300 CAGATCGAGCGCACCATCGCCGTCAACTACACCGGCCTGGTGAACACCACGACTGCCATC

gi|9229|e       120 L  D  F  W  D  K  R  K  G  G  P  G  G  I  I  C  N  I  G  S  
gi|9229|e       360 CTGGACTTCTGGGACAAGCGCAAGGGTGGACCCGGTGGTATCATCTGCAACATTGGATCC

gi|9229|e       140 V  T  G  F  N  A  I  Y  Q  V  P  V  Y  S  G  T  K  A  A  V  
gi|9229|e       420 GTGACTGGATTCAACGCCATCTACCAGGTGCCCGTTTACTCCGGCACCAAGGCTGCCGTG

gi|9229|e       160 V  N  F  T  S  S  L  A  K  L  A  P  I  T  G  V  T  A  Y  T  
gi|9229|e       480 GTCAACTTCACCAGCTCCCTGGCGAAACTGGCCCCCATCACCGGCGTGACCGCTTACACC

gi|9229|e       180 V  N  P  G  I  T  R  T  T  L  V  H  K  F  N  S  W  L  D  V  
gi|9229|e       540 GTGAACCCCGGCATCACCCGCACCACCCTGGTGCACAAGTTCAACTCCTGGCTGGATGTT

gi|9229|e       200 E  P  Q  V  A  E  K  L  L  A  H  P  T  Q  P  S  L  A  C  A  
gi|9229|e       600 GAGCCCCAGGTGGCCGAGAAGCTCCTGGCTCACCCCACCCAGCCCTCGTTGGCCTGCGCC

gi|9229|e       220 Q  N  F  V  K  A  I  E  L  N  Q  N  G  A  I  W  K  L  D  L  
gi|9229|e       660 CAGAACTTTGTGAAGGCCATCGAGCTGAACCAGAACGGTGCCATCTGGAAACTGGACTTG

gi|9229|e       240 G  T  L  E  A  I  Q  W  S  K  H  W  D  S  G  I   256
gi|9229|e       720 GGCACCCTGGAGGCCATCCAGTGGTCCAAGCACTGGGACTCCGGCATC 768
""",
            )
            protein_record = protein_alignment.sequences[7]
            nucleotide_record = nucleotide_records[protein_record.id]
            self.assertEqual(nucleotide_record.id, protein_record.id)
            alignments = aligner.align(protein_record, nucleotide_record)
            self.assertEqual(len(alignments), 1)
            alignment = next(alignments)
            codon_alignments.append(alignment)
            self.assertTrue(
                np.array_equal(alignment.coordinates, np.array([[0, 256], [0, 768]]))
            )
            self.assertEqual(
                str(alignment),
                """\
gi|9231|e         0 M  A  F  T  L  T  N  K  N  V  V  F  V  A  G  L  G  G  I  G  
gi|9231|e         0 ATGGCGTTTACCTTGACCAACAAGAACGTGGTTTTCGTGGCCGGTCTGGGAGGCATTGGT

gi|9231|e        20 L  D  T  S  K  E  L  V  K  R  D  L  K  N  L  V  I  L  D  R  
gi|9231|e        60 CTGGACACCAGCAAGGAGCTGGTCAAGCGCGACCTGAAGAACCTGGTGATCCTCGACCGC

gi|9231|e        40 I  E  N  P  A  A  I  A  E  L  K  A  I  N  P  K  V  T  V  T  
gi|9231|e       120 ATTGAGAACCCGGCTGCCATCGCCGAGCTGAAGGCAATCAATCCAAAGGTGACCGTCACC

gi|9231|e        60 F  Y  P  Y  D  V  T  V  P  I  A  E  T  T  K  L  L  K  T  I  
gi|9231|e       180 TTCTACCCCTATGATGTGACCGTGCCCATTGCCGAGACCACCAAGCTGCTGAAGACCATC

gi|9231|e        80 F  A  Q  L  K  T  I  D  V  L  I  N  G  A  G  I  L  D  D  H  
gi|9231|e       240 TTCGCCCAGCTGAAGACCATCGATGTCCTGATCAACGGAGCTGGCATCCTGGACGATCAC

gi|9231|e       100 Q  I  E  R  T  I  A  V  N  Y  T  G  L  V  N  T  T  T  A  I  
gi|9231|e       300 CAGATCGAGCGCACCATCGCCGTCAACTACACCGGCCTGGTGAACACCACGACTGCCATC

gi|9231|e       120 L  D  F  W  D  K  R  K  G  G  P  G  G  I  I  C  N  I  G  S  
gi|9231|e       360 CTGGACTTCTGGGACAAGCGCAAGGGTGGACCCGGTGGTATCATCTGCAACATTGGATCC

gi|9231|e       140 V  T  G  F  N  A  I  Y  Q  V  P  V  Y  S  G  T  K  A  A  V  
gi|9231|e       420 GTGACTGGATTCAACGCCATCTACCAGGTGCCCGTTTACTCCGGCACCAAGGCTGCCGTG

gi|9231|e       160 V  N  F  T  S  S  L  A  K  L  A  P  I  T  G  V  T  A  Y  T  
gi|9231|e       480 GTCAACTTCACCAGCTCCCTGGCGAAACTGGCCCCCATCACCGGCGTGACCGCTTACACC

gi|9231|e       180 V  N  P  G  I  T  R  T  T  L  V  H  K  F  N  S  W  L  D  V  
gi|9231|e       540 GTGAACCCCGGCATCACCCGCACCACCCTGGTGCACAAGTTCAACTCCTGGCTGGATGTT

gi|9231|e       200 E  P  Q  V  A  E  K  L  L  A  H  P  T  Q  P  S  L  A  C  A  
gi|9231|e       600 GAGCCCCAGGTGGCCGAGAAGCTCCTGGCTCACCCCACCCAGCCCTCGTTGGCCTGCGCC

gi|9231|e       220 Q  N  F  V  K  A  I  E  L  N  Q  N  G  A  I  W  K  L  D  L  
gi|9231|e       660 CAGAACTTTGTGAAGGCCATCGAGCTGAACCAGAACGGTGCTATCTGGAAACTGGACTTG

gi|9231|e       240 G  T  L  E  A  I  Q  W  S  K  H  W  D  S  G  I   256
gi|9231|e       720 GGCACCCTGGAGGCCATCCAGTGGTCCAAGCACTGGGACTCCGGCATC 768
""",
            )
            protein_record = protein_alignment.sequences[8]
            nucleotide_record = nucleotide_records[protein_record.id]
            self.assertEqual(nucleotide_record.id, protein_record.id)
            alignments = aligner.align(protein_record, nucleotide_record)
            self.assertEqual(len(alignments), 1)
            alignment = next(alignments)
            codon_alignments.append(alignment)
            self.assertTrue(
                np.array_equal(alignment.coordinates, np.array([[0, 256], [0, 768]]))
            )
            self.assertEqual(
                str(alignment),
                """\
gi|9233|e         0 M  A  F  T  L  T  N  K  N  V  V  F  V  A  G  L  G  G  I  G  
gi|9233|e         0 ATGGCGTTTACCTTGACCAACAAGAACGTGGTTTTCGTGGCCGGTCTGGGAGGCATTGGT

gi|9233|e        20 L  D  T  S  K  E  L  V  K  R  D  L  K  N  L  V  I  L  D  R  
gi|9233|e        60 CTGGACACCAGCAAGGAGCTGGTCAAGCGGGACCTGAAGAACCTGGTGATCCTCGACCGC

gi|9233|e        40 I  E  N  P  A  A  I  A  E  L  K  A  I  N  P  K  V  T  V  T  
gi|9233|e       120 ATTGAGAACCCGGCTGCCATCGCCGAGCTGAAGGCAATCAATCCAAAGGTGACCGTCACC

gi|9233|e        60 F  Y  P  Y  D  V  T  V  P  I  A  E  T  T  K  L  L  K  T  I  
gi|9233|e       180 TTCTACCCATACGATGTGACCGTGCCCATTGCCGAGACCACCAAGCTGCTGAAGACCATC

gi|9233|e        80 F  A  Q  L  K  T  I  D  V  L  I  N  G  A  G  I  L  D  D  H  
gi|9233|e       240 TTCGCCCAGCTGAAGACCATCGATGTCCTGATCAACGGAGCTGGCATCCTGGACGATCAC

gi|9233|e       100 Q  I  E  R  T  I  A  V  N  Y  T  G  L  V  N  T  T  T  A  I  
gi|9233|e       300 CAGATCGAGCGCACCATCGCCGTCAACTACACCGGCCTGGTGAACACCACGACGGCCATC

gi|9233|e       120 L  D  F  W  D  K  R  K  G  G  P  G  G  I  I  C  N  I  G  S  
gi|9233|e       360 CTGGACTTCTGGGACAAGCGCAAGGGTGGACCCGGTGGTATCATCTGCAACATTGGATCC

gi|9233|e       140 V  T  G  F  N  A  I  Y  Q  V  P  V  Y  S  G  T  K  A  A  V  
gi|9233|e       420 GTGACTGGATTCAACGCCATCTACCAGGTGCCCGTTTACTCCGGCACCAAGGCTGCCGTG

gi|9233|e       160 V  N  F  T  S  S  L  A  K  L  A  P  I  T  G  V  T  A  Y  T  
gi|9233|e       480 GTCAACTTCACCAGCTCCCTGGCGAAACTGGCCCCCATCACCGGCGTGACCGCTTACACC

gi|9233|e       180 V  N  P  G  I  T  R  T  T  L  V  H  K  F  N  S  W  L  D  V  
gi|9233|e       540 GTGAACCCCGGCATCACCCGCACCACCCTGGTGCACAAGTTCAACTCCTGGCTGGATGTT

gi|9233|e       200 E  P  Q  V  A  E  K  L  L  A  H  P  T  Q  P  S  L  A  C  A  
gi|9233|e       600 GAGCCCCAGGTGGCCGAGAAGCTCCTGGCTCACCCCACCCAGCCCTCGTTGGCCTGCGCC

gi|9233|e       220 Q  N  F  V  K  A  I  E  L  N  Q  N  G  A  I  W  K  L  D  L  
gi|9233|e       660 CAGAACTTTGTGAAGGCCATCGAGCTGAACCAGAACGGTGCCATCTGGAAACTGGACTTG

gi|9233|e       240 G  T  L  E  A  I  Q  W  S  K  H  W  D  S  G  I   256
gi|9233|e       720 GGCACCCTGGAGGCCATCCAGTGGTCCAAGCACTGGGACTCCGGCATC 768
""",
            )
            protein_record = protein_alignment.sequences[9]
            nucleotide_record = nucleotide_records[protein_record.id]
            self.assertEqual(nucleotide_record.id, protein_record.id)
            alignments = aligner.align(protein_record, nucleotide_record)
            self.assertEqual(len(alignments), 1)
            alignment = next(alignments)
            codon_alignments.append(alignment)
            self.assertTrue(
                np.array_equal(alignment.coordinates, np.array([[0, 256], [0, 768]]))
            )
            self.assertEqual(
                str(alignment),
                """\
gi|9235|e         0 M  A  F  T  L  T  N  K  N  V  V  F  V  A  G  L  G  G  I  G  
gi|9235|e         0 ATGGCGTTTACCTTGACCAACAAGAACGTGGTTTTCGTGGCCGGTCTGGGAGGCATTGGT

gi|9235|e        20 L  D  T  S  K  E  L  V  K  R  D  L  K  N  L  V  I  L  D  R  
gi|9235|e        60 CTGGACACCAGCAAGGAGCTGGTCAAGCGGGACCTGAAGAACCTGGTGATCCTCGACCGC

gi|9235|e        40 I  E  N  P  A  A  I  A  E  L  K  A  I  N  P  K  V  T  V  T  
gi|9235|e       120 ATTGAGAACCCGGCTGCCATCGCCGAGCTGAAGGCAATCAATCCAAAGGTGACCGTCACC

gi|9235|e        60 F  Y  P  Y  D  V  T  V  P  I  A  E  T  T  K  L  L  K  T  I  
gi|9235|e       180 TTCTACCCCTATGATGTGACCGTGCCCATTGCCGAGACCACCAAGCTGCTGAAGACCATC

gi|9235|e        80 F  A  Q  L  K  T  I  D  V  L  I  N  G  A  G  I  L  D  D  H  
gi|9235|e       240 TTCGCCCAGCTGAAGACCATCGATGTCCTGATCAACGGAGCTGGCATCCTGGACGATCAC

gi|9235|e       100 Q  I  E  R  T  I  A  V  N  Y  T  G  L  V  N  T  T  T  A  I  
gi|9235|e       300 CAGATCGAGCGCACCATCGCCGTCAACTACACCGGCCTGGTGAACACCACGACTGCCATC

gi|9235|e       120 L  D  F  W  D  K  R  K  G  G  P  G  G  I  I  C  N  I  G  S  
gi|9235|e       360 CTGGACTTCTGGGACAAGCGCAAGGGTGGACCCGGTGGTATCATCTGCAACATTGGATCC

gi|9235|e       140 V  T  G  F  N  A  I  Y  Q  V  P  V  Y  S  G  T  K  A  A  V  
gi|9235|e       420 GTGACTGGATTCAACGCCATCTACCAGGTGCCCGTTTACTCCGGCACCAAGGCTGCCGTG

gi|9235|e       160 V  N  F  T  S  S  L  A  K  L  A  P  I  T  G  V  T  A  Y  T  
gi|9235|e       480 GTCAACTTCACCAGCTCCCTGGCGAAACTGGCCCCCATCACCGGCGTGACCGCTTACACC

gi|9235|e       180 V  N  P  G  I  T  R  T  T  L  V  H  K  F  N  S  W  L  D  V  
gi|9235|e       540 GTGAACCCCGGCATCACCCGCACCACCCTGGTGCACAAGTTCAACTCCTGGCTGGATGTT

gi|9235|e       200 E  P  Q  V  A  E  K  L  L  A  H  P  T  Q  P  S  L  A  C  A  
gi|9235|e       600 GAGCCCCAGGTGGCCGAGAAGCTCCTGGCTCACCCCACCCAGCCCTCGTTGGCCTGCGCC

gi|9235|e       220 Q  N  F  V  K  A  I  E  L  N  Q  N  G  A  I  W  K  L  D  L  
gi|9235|e       660 CAGAACTTTGTGAAGGCCATCGAGCTGAACCAGAACGGTGCCATCTGGAAACTGGACTTG

gi|9235|e       240 G  T  L  E  A  I  Q  W  S  K  H  W  D  S  G  I   256
gi|9235|e       720 GGCACCCTGGAGGCCATCCAGTGGTCCAAGCATTGGGACTCCGGCATC 768
""",
            )
            protein_record = protein_alignment.sequences[10]
            nucleotide_record = nucleotide_records[protein_record.id]
            self.assertEqual(nucleotide_record.id, protein_record.id)
            alignments = aligner.align(protein_record, nucleotide_record)
            self.assertEqual(len(alignments), 1)
            alignment = next(alignments)
            codon_alignments.append(alignment)
            self.assertTrue(
                np.array_equal(alignment.coordinates, np.array([[0, 256], [0, 768]]))
            )
            self.assertEqual(
                str(alignment),
                """\
gi|9237|e         0 M  A  F  T  L  T  N  K  N  V  V  F  V  A  G  L  G  G  I  G  
gi|9237|e         0 ATGGCGTTTACCTTGACCAACAAGAACGTGGTTTTCGTGGCCGGTCTGGGAGGCATTGGT

gi|9237|e        20 L  D  T  S  K  E  L  V  K  R  D  L  K  N  L  V  I  L  D  R  
gi|9237|e        60 CTGGACACCAGCAAGGAGCTGGTCAAGCGGGACCTGAAGAACCTGGTGATCCTCGACCGC

gi|9237|e        40 I  E  N  P  A  A  I  A  E  L  K  A  I  N  P  K  V  T  V  T  
gi|9237|e       120 ATTGAGAACCCGGCTGCCATCGCCGAGCTGAAGGCAATCAATCCAAAGGTGACCGTCACC

gi|9237|e        60 F  Y  P  Y  D  V  T  V  P  I  A  E  T  T  K  L  L  K  T  I  
gi|9237|e       180 TTCTACCCCTACGATGTGACCGTGCCCATTGCCGAGACCACCAAGCTGCTGAAGACCATC

gi|9237|e        80 F  A  Q  L  K  T  I  D  V  L  I  N  G  A  G  I  L  D  D  H  
gi|9237|e       240 TTCGCCCAGCTGAAGACCATCGATGTCCTGATCAACGGAGCTGGCATCCTGGACGATCAC

gi|9237|e       100 Q  I  E  R  T  I  A  V  N  Y  T  G  L  V  N  T  T  T  A  I  
gi|9237|e       300 CAGATCGAGCGCACCATCGCCGTCAACTACACCGGACTGGTGAACACCACGACGGCCATC

gi|9237|e       120 L  D  F  W  D  K  R  K  G  G  P  G  G  I  I  C  N  I  G  S  
gi|9237|e       360 CTGGACTTCTGGGACAAGCGCAAGGGTGGACCCGGTGGTATCATCTGCAACATTGGATCC

gi|9237|e       140 V  T  G  F  N  A  I  Y  Q  V  P  V  Y  S  G  T  K  A  A  V  
gi|9237|e       420 GTGACTGGATTCAACGCCATCTACCAGGTGCCCGTTTACTCCGGCACCAAGGCTGCCGTG

gi|9237|e       160 V  N  F  T  S  S  L  A  K  L  A  P  I  T  G  V  T  A  Y  T  
gi|9237|e       480 GTCAACTTCACCAGCTCCCTGGCGAAACTGGCCCCCATCACCGGCGTGACCGCTTACACC

gi|9237|e       180 V  N  P  G  I  T  R  T  T  L  V  H  K  F  N  S  W  L  D  V  
gi|9237|e       540 GTGAACCCCGGCATCACCCGCACCACCCTGGTGCACAAGTTCAACTCCTGGCTGGATGTT

gi|9237|e       200 E  P  Q  V  A  E  K  L  L  A  H  P  T  Q  P  S  L  A  C  A  
gi|9237|e       600 GAGCCCCAGGTGGCCGAGAAGCTCCTGGCTCACCCCACCCAGCCCTCGTTGGCCTGCGCC

gi|9237|e       220 Q  N  F  V  K  A  I  E  L  N  Q  N  G  A  I  W  K  L  D  L  
gi|9237|e       660 CAGAACTTTGTGAAGGCCATCGAGCTGAACCAGAACGGTGCCATCTGGAAACTGGACTTG

gi|9237|e       240 G  T  L  E  A  I  Q  W  S  K  H  W  D  S  G  I   256
gi|9237|e       720 GGCACCCTGGAGGCCATCCAGTGGTCCAAGCACTGGGACTCCGGCATC 768
""",
            )
            protein_record = protein_alignment.sequences[11]
            nucleotide_record = nucleotide_records[protein_record.id]
            self.assertEqual(nucleotide_record.id, protein_record.id)
            alignments = aligner.align(protein_record, nucleotide_record)
            self.assertEqual(len(alignments), 1)
            alignment = next(alignments)
            codon_alignments.append(alignment)
            self.assertTrue(
                np.array_equal(alignment.coordinates, np.array([[0, 256], [0, 768]]))
            )
            self.assertEqual(
                str(alignment),
                """\
gi|9239|e         0 M  A  F  T  L  T  N  K  N  V  V  F  V  A  G  L  G  G  I  G  
gi|9239|e         0 ATGGCGTTTACCTTGACCAACAAGAACGTGGTTTTCGTGGCCGGTCTGGGAGGCATTGGT

gi|9239|e        20 L  D  T  S  K  E  L  V  K  R  D  L  K  N  L  V  I  L  D  R  
gi|9239|e        60 CTGGACACCAGCAAGGAGCTGGTCAAGCGGGACCTGAAGAACCTGGTGATCCTCGACCGC

gi|9239|e        40 I  E  N  P  A  A  I  A  E  L  K  A  I  N  P  K  V  T  V  T  
gi|9239|e       120 ATTGAGAACCCGGCTGCCATCGCCGAGCTGAAGGCAATCAATCCAAAGGTGACCGTCACC

gi|9239|e        60 F  Y  P  Y  D  V  T  V  P  I  A  E  T  T  K  L  L  K  T  I  
gi|9239|e       180 TTCTACCCCTACGATGTGACCGTGCCCATTGCCGAGACCACCAAGCTGCTGAAGACCATC

gi|9239|e        80 F  A  Q  L  K  T  I  D  V  L  I  N  G  A  G  I  L  D  D  H  
gi|9239|e       240 TTCGCCCAGCTGAAGACCATCGATGTCCTGATCAACGGAGCTGGCATCCTGGACGATCAC

gi|9239|e       100 Q  I  E  R  T  I  A  V  N  Y  T  G  L  V  N  T  T  T  A  I  
gi|9239|e       300 CAGATCGAGCGCACCATCGCCGTCAACTACACCGGCCTGGTGAACACCACGACGGCCATC

gi|9239|e       120 L  D  F  W  D  K  R  K  G  G  P  G  G  I  I  C  N  I  G  S  
gi|9239|e       360 CTGGACTTCTGGGACAAGCGCAAGGGTGGACCCGGTGGTATCATCTGCAACATTGGATCC

gi|9239|e       140 V  T  G  F  N  A  I  Y  Q  V  P  V  Y  S  G  T  K  A  A  V  
gi|9239|e       420 GTGACTGGATTCAACGCCATCTACCAGGTGCCCGTTTACTCCGGCACCAAGGCTGCCGTG

gi|9239|e       160 V  N  F  T  S  S  L  A  K  L  A  P  I  T  G  V  T  A  Y  T  
gi|9239|e       480 GTCAACTTCACCAGCTCCCTGGCGAAACTGGCCCCCATCACCGGCGTGACCGCTTACACC

gi|9239|e       180 V  N  P  G  I  T  R  T  T  L  V  H  K  F  N  S  W  L  D  V  
gi|9239|e       540 GTGAACCCCGGCATCACCCGCACCACCCTGGTGCACAAGTTCAACTCCTGGCTGGATGTT

gi|9239|e       200 E  P  Q  V  A  E  K  L  L  A  H  P  T  Q  P  S  L  A  C  A  
gi|9239|e       600 GAGCCCCAGGTGGCCGAGAAGCTCCTGGCTCACCCCACCCAGCCCTCGTTGGCCTGCGCC

gi|9239|e       220 Q  N  F  V  K  A  I  E  L  N  Q  N  G  A  I  W  K  L  D  L  
gi|9239|e       660 CAGAACTTTGTCAAGGCCATCGAGCTGAACCAGAACGGTGCCATCTGGAAACTGGACTTG

gi|9239|e       240 G  T  L  E  A  I  Q  W  S  K  H  W  D  S  G  I   256
gi|9239|e       720 GGCACCCTGGAGGCCATCCAGTGGTCCAAGCACTGGGACTCCGGCATC 768
""",
            )
            protein_record = protein_alignment.sequences[12]
            nucleotide_record = nucleotide_records[protein_record.id]
            self.assertEqual(nucleotide_record.id, protein_record.id)
            alignments = aligner.align(protein_record, nucleotide_record)
            self.assertEqual(len(alignments), 1)
            alignment = next(alignments)
            codon_alignments.append(alignment)
            self.assertTrue(
                np.array_equal(alignment.coordinates, np.array([[0, 256], [0, 768]]))
            )
            self.assertEqual(
                str(alignment),
                """\
gi|9097|e         0 M  A  F  T  L  T  N  K  N  V  I  F  V  A  G  L  G  G  I  G  
gi|9097|e         0 ATGGCGTTTACTTTGACCAACAAGAACGTGATTTTCGTTGCCGGTCTGGGAGGCATTGGT

gi|9097|e        20 L  D  T  S  K  E  L  L  K  R  D  L  K  N  L  V  I  L  D  R  
gi|9097|e        60 CTGGACACCAGCAAGGAGCTGCTCAAGCGCGACCTGAAGAACCTGGTGATCCTCGACCGC

gi|9097|e        40 I  E  N  P  A  A  I  A  E  L  K  A  I  N  P  K  V  T  V  T  
gi|9097|e       120 ATTGAGAACCCTGCTGCCATTGCCGAGCTGAAGGCAATCAATCCAAAGGTGACCGTCACC

gi|9097|e        60 F  Y  P  Y  D  V  T  V  P  I  A  E  T  T  K  L  L  K  T  I  
gi|9097|e       180 TTCTACCCCTATGATGTGACCGTGCCCATTGCCGAGACCACCAAGCTGCTGAAGACCATC

gi|9097|e        80 F  A  K  L  K  T  V  D  V  L  I  N  G  A  G  I  L  D  D  H  
gi|9097|e       240 TTCGCCAAGCTGAAGACCGTCGATGTCCTGATCAACGGAGCTGGTATCCTGGACGATCAC

gi|9097|e       100 Q  I  E  R  T  I  A  V  N  Y  T  G  L  V  N  T  T  T  A  I  
gi|9097|e       300 CAGATCGAGCGCACCATTGCCGTCAACTACACTGGCCTGGTCAACACCACGACGGCCATT

gi|9097|e       120 L  D  F  W  D  K  R  K  G  G  P  G  G  I  I  C  N  I  G  S  
gi|9097|e       360 TTGGACTTCTGGGACAAGCGCAAGGGTGGTCCCGGTGGTATCATCTGCAACATTGGATCC

gi|9097|e       140 V  T  G  F  N  A  I  Y  Q  V  P  V  Y  S  G  T  K  A  A  V  
gi|9097|e       420 GTCACTGGATTCAATGCCATCTACCAGGTGCCCGTCTACTCCGGCACCAAGGCTGCCGTG

gi|9097|e       160 V  N  F  T  S  S  L  A  K  L  A  P  I  T  G  V  T  A  Y  T  
gi|9097|e       480 GTCAACTTCACCAGCTCCCTGGCGAAACTGGCCCCCATTACCGGCGTGACCGCTTACACC

gi|9097|e       180 V  N  P  G  I  T  R  T  T  L  V  H  K  F  N  S  W  L  D  V  
gi|9097|e       540 GTGAACCCCGGCATCACCCGCACCACCCTGGTGCACAAGTTCAACTCCTGGCTGGATGTT

gi|9097|e       200 E  P  Q  V  A  E  K  L  L  A  H  P  T  Q  P  S  L  A  C  A  
gi|9097|e       600 GAGCCCCAGGTTGCCGAGAAGCTCCTGGCTCATCCCACCCAGCCCTCGTTGGCCTGCGCC

gi|9097|e       220 E  N  F  V  K  A  I  E  L  N  Q  N  G  A  I  W  K  L  D  L  
gi|9097|e       660 GAGAACTTCGTCAAGGCTATCGAGCTGAACCAGAACGGAGCCATCTGGAAACTGGACTTG

gi|9097|e       240 G  T  L  E  A  I  Q  W  T  K  H  W  D  S  G  I   256
gi|9097|e       720 GGCACCCTGGAGGCCATCCAGTGGACCAAGCACTGGGACTCCGGCATC 768
""",
            )
            protein_record = protein_alignment.sequences[13]
            nucleotide_record = nucleotide_records[protein_record.id]
            self.assertEqual(nucleotide_record.id, protein_record.id)
            alignments = aligner.align(protein_record, nucleotide_record)
            self.assertEqual(len(alignments), 1)
            alignment = next(alignments)
            codon_alignments.append(alignment)
            self.assertTrue(
                np.array_equal(alignment.coordinates, np.array([[0, 256], [0, 768]]))
            )
            self.assertEqual(
                str(alignment),
                """\
gi|9099|e         0 M  A  F  T  L  T  N  K  N  V  I  F  V  A  G  L  G  G  I  G  
gi|9099|e         0 ATGGCGTTTACTTTGACCAACAAGAACGTGATTTTCGTTGCCGGTCTGGGAGGCATTGGT

gi|9099|e        20 L  D  T  S  K  E  L  L  K  R  D  L  K  N  L  V  I  L  D  R  
gi|9099|e        60 CTGGACACCAGCAAGGAGCTGCTCAAGCGCGACCTGAAGAACCTGGTGATCCTCGACCGC

gi|9099|e        40 I  E  N  P  A  A  I  A  E  L  K  A  I  N  P  K  V  T  V  T  
gi|9099|e       120 ATTGAGAACCCTGCTGCCATTGCCGAGCTGAAGGCAATCAATCCAAAGGTGACCGTCACC

gi|9099|e        60 F  Y  P  Y  D  V  T  V  P  I  A  E  T  T  K  L  L  K  T  I  
gi|9099|e       180 TTCTACCCCTATGATGTGACCGTGCCCATTGCCGAGACCACCAAGCTGCTGAAGACCATC

gi|9099|e        80 F  A  K  L  K  T  V  D  V  L  I  N  G  A  G  I  L  D  D  H  
gi|9099|e       240 TTCGCCAAGCTGAAGACCGTCGATGTCCTGATCAACGGAGCTGGTATCCTGGACGATCAC

gi|9099|e       100 Q  I  E  R  T  I  A  V  N  Y  T  G  L  V  N  T  T  T  A  I  
gi|9099|e       300 CAGATCGAGCGCACCATTGCCGTCAACTACACTGGCCTGGTCAACACCACGACGGCCATT

gi|9099|e       120 L  D  F  W  D  K  R  K  G  G  P  G  G  I  I  C  N  I  G  S  
gi|9099|e       360 CTGGACTTCTGGGACAAGCGCAAGGGTGGTCCCGGTGGTATCATCTGCAACATTGGATCC

gi|9099|e       140 V  T  G  F  N  A  I  Y  Q  V  P  V  Y  S  G  T  K  A  A  V  
gi|9099|e       420 GTCACTGGTTTCAATGCCATCTACCAGGTGCCCGTCTACTCCGGCACCAAGGCTGCCGTG

gi|9099|e       160 V  N  F  T  S  S  L  A  K  L  A  P  I  T  G  V  T  A  Y  T  
gi|9099|e       480 GTCAACTTCACCAGCTCCCTGGCGAAACTGGCCCCCATTACCGGCGTGACCGCTTACACC

gi|9099|e       180 V  N  P  G  I  T  R  T  T  L  V  H  K  F  N  S  W  L  D  V  
gi|9099|e       540 GTGAACCCCGGCATCACCCGCACCACCCTGGTGCACAAGTTCAACTCCTGGCTGGATGTT

gi|9099|e       200 E  P  Q  V  A  E  K  L  L  A  H  P  T  Q  P  S  L  A  C  A  
gi|9099|e       600 GAGCCCCAGGTTGCCGAGAAGCTCCTGGCTCATCCCACTCAGCCCTCATTGGCCTGCGCC

gi|9099|e       220 E  N  F  V  K  A  I  E  L  N  Q  N  G  A  I  W  K  L  D  L  
gi|9099|e       660 GAGAACTTCGTCAAGGCCATCGAGCTGAACCAGAACGGAGCCATCTGGAAACTGGACTTG

gi|9099|e       240 G  T  L  E  A  I  Q  W  T  K  H  W  D  S  G  I   256
gi|9099|e       720 GGCACCCTGGAGGCCATCCAGTGGACCAAGCACTGGGACTCCGGCATC 768
""",
            )
            protein_record = protein_alignment.sequences[14]
            nucleotide_record = nucleotide_records[protein_record.id]
            self.assertEqual(nucleotide_record.id, protein_record.id)
            alignments = aligner.align(protein_record, nucleotide_record)
            self.assertEqual(len(alignments), 1)
            alignment = next(alignments)
            codon_alignments.append(alignment)
            self.assertTrue(
                np.array_equal(alignment.coordinates, np.array([[0, 256], [0, 768]]))
            )
            self.assertEqual(
                str(alignment),
                """\
gi|9101|e         0 M  A  F  T  L  T  N  K  N  V  I  F  V  A  G  L  G  G  I  G  
gi|9101|e         0 ATGGCGTTTACTTTGACCAACAAGAACGTGATTTTCGTTGCCGGTCTGGGAGGCATTGGT

gi|9101|e        20 L  D  T  S  K  E  L  L  K  R  D  L  K  N  L  V  I  L  D  R  
gi|9101|e        60 CTGGACACCAGCAAGGAGCTGCTCAAGCGCGACCTGAAGAACCTGGTGATCCTCGACCGC

gi|9101|e        40 I  E  N  P  A  A  I  A  E  L  K  A  I  N  P  K  V  T  V  T  
gi|9101|e       120 ATTGAGAACCCTGCTGCCATTGCCGAGCTGAAGGCAATCAATCCAAAGGTGACCGTCACC

gi|9101|e        60 F  Y  P  Y  D  V  T  V  P  I  A  E  T  T  K  L  L  K  T  I  
gi|9101|e       180 TTCTACCCCTATGATGTGACCGTGCCCATTGCCGAGACCACCAAGCTGCTGAAGACCATC

gi|9101|e        80 F  A  K  L  K  T  V  D  V  L  I  N  G  A  G  I  L  D  D  H  
gi|9101|e       240 TTCGCCAAGCTGAAGACCGTCGATGTCCTGATCAACGGAGCTGGTATCCTGGACGATCAC

gi|9101|e       100 Q  I  E  R  T  I  A  V  N  Y  T  G  L  V  N  T  T  T  A  I  
gi|9101|e       300 CAGATCGAGCGCACCATTGCCGTCAACTACACTGGCCTGGTCAACACCACGACGGCCATT

gi|9101|e       120 L  D  F  W  D  K  R  K  G  G  P  G  G  I  I  C  N  I  G  S  
gi|9101|e       360 CTGGACTTCTGGGACAAGCGCAAGGGTGGTCCCGGTGGTATCATCTGCAACATTGGATCC

gi|9101|e       140 V  T  G  F  N  A  I  Y  Q  V  P  V  Y  S  G  T  K  A  A  V  
gi|9101|e       420 GTCACTGGATTCAATGCCATCTACCAGGTGCCCGTCTACTCTGGCACCAAGGCCGCCGTG

gi|9101|e       160 V  N  F  T  S  S  L  A  K  L  A  P  I  T  G  V  T  A  Y  T  
gi|9101|e       480 GTCAACTTCACCAGCTCCCTGGCGAAACTGGCCCCCATTACCGGCGTGACCGCTTACACC

gi|9101|e       180 V  N  P  G  I  T  R  T  T  L  V  H  K  F  N  S  W  L  D  V  
gi|9101|e       540 GTGAACCCCGGCATCACCCGCACCACCCTGGTGCACAAGTTCAACTCCTGGCTGGATGTT

gi|9101|e       200 E  P  Q  V  A  E  K  L  L  A  H  P  T  Q  P  S  L  A  C  A  
gi|9101|e       600 GAGCCCCAGGTTGCCGAGAAGCTCCTGGCTCATCCCACCCAGCCCTCGTTGGCCTGCGCC

gi|9101|e       220 E  N  F  V  K  A  I  E  L  N  Q  N  G  A  I  W  K  L  D  L  
gi|9101|e       660 GAGAACTTCGTCAAGGCCATCGAGCTGAACCAGAACGGAGCCATCTGGAAACTGGACTTG

gi|9101|e       240 G  T  L  E  A  I  Q  W  T  K  H  W  D  S  G  I   256
gi|9101|e       720 GGCACCCTGGAGGCCATCCAGTGGACCAAGCACTGGGACTCCGGCATC 768
""",
            )
            protein_record = protein_alignment.sequences[15]
            nucleotide_record = nucleotide_records[protein_record.id]
            self.assertEqual(nucleotide_record.id, protein_record.id)
            alignments = aligner.align(protein_record, nucleotide_record)
            self.assertEqual(len(alignments), 1)
            alignment = next(alignments)
            codon_alignments.append(alignment)
            self.assertTrue(
                np.array_equal(alignment.coordinates, np.array([[0, 256], [0, 768]]))
            )
            self.assertEqual(
                str(alignment),
                """\
gi|9103|e         0 M  A  F  T  L  T  N  K  N  V  I  F  V  A  G  L  G  G  I  G  
gi|9103|e         0 ATGGCGTTTACTTTGACCAACAAGAACGTGATTTTCGTTGCCGGTCTGGGAGGCATCGGT

gi|9103|e        20 L  D  T  S  K  E  L  L  K  R  D  L  K  N  L  V  I  L  D  R  
gi|9103|e        60 CTGGACACCAGCAAGGAGCTGCTCAAGCGCGACCTGAAGAACCTGGTGATCCTCGACCGC

gi|9103|e        40 I  E  N  P  A  A  I  A  E  L  K  A  I  N  P  K  V  T  V  T  
gi|9103|e       120 ATTGAGAACCCTGCTGCCATTGCCGAGCTGAAGGCAATCAATCCAAAGGTGACCGTCACC

gi|9103|e        60 F  Y  P  Y  D  V  T  V  P  I  A  E  T  T  K  L  L  K  T  I  
gi|9103|e       180 TTCTACCCCTATGATGTGACCGTGCCCATTGCCGAGACCACCAAGCTGCTGAAGACCATC

gi|9103|e        80 F  A  K  L  K  T  V  D  V  L  I  N  G  A  G  I  L  D  D  H  
gi|9103|e       240 TTCGCCAAGCTGAAGACCGTCGATGTCCTGATCAACGGAGCTGGTATCCTGGACGATCAC

gi|9103|e       100 Q  I  E  R  T  I  A  V  N  Y  T  G  L  V  N  T  T  T  A  I  
gi|9103|e       300 CAGATCGAGCGCACCATTGCCGTCAACTACACTGGCCTGGTCAACACCACGACGGCCATT

gi|9103|e       120 L  D  F  W  D  K  R  K  G  G  P  G  G  I  I  C  N  I  G  S  
gi|9103|e       360 CTGGACTTCTGGGACAAGCGCAAGGGTGGTCCCGGTGGTATCATCTGCAACATTGGATCC

gi|9103|e       140 V  T  G  F  N  A  I  Y  Q  V  P  V  Y  S  G  T  K  A  A  V  
gi|9103|e       420 GTCACTGGATTCAATGCCATCTACCAGGTGCCCGTCTACTCCGGCACCAAGGCCGCCGTG

gi|9103|e       160 V  N  F  T  S  S  L  A  K  L  A  P  I  T  G  V  T  A  Y  T  
gi|9103|e       480 GTCAACTTCACCAGCTCCCTGGCGAAACTGGCCCCCATTACCGGCGTGACCGCTTACACC

gi|9103|e       180 V  N  P  G  I  T  R  T  T  L  V  H  K  F  N  S  W  L  D  V  
gi|9103|e       540 GTGAACCCCGGCATCACCCGCACCACCCTGGTGCACAAGTTCAACTCCTGGCTGGATGTT

gi|9103|e       200 E  P  Q  V  A  E  K  L  L  A  H  P  T  Q  P  S  L  A  C  A  
gi|9103|e       600 GAGCCCCAGGTTGCCGAGAAGCTCCTGGCTCATCCCACCCAGCCCTCGTTGGCCTGCGCC

gi|9103|e       220 E  N  F  V  K  A  I  E  L  N  Q  N  G  A  I  W  K  L  D  L  
gi|9103|e       660 GAGAACTTCGTCAAGGCCATCGAGCTGAACCAGAACGGAGCCATCTGGAAACTGGACTTG

gi|9103|e       240 G  T  L  E  A  I  Q  W  T  K  H  W  D  S  G  I   256
gi|9103|e       720 GGCACCCTGGAGGCCATCCAGTGGACCAAGCACTGGGACTCCGGCATC 768
""",
            )
            protein_record = protein_alignment.sequences[16]
            nucleotide_record = nucleotide_records[protein_record.id]
            self.assertEqual(nucleotide_record.id, protein_record.id)
            alignments = aligner.align(protein_record, nucleotide_record)
            self.assertEqual(len(alignments), 1)
            alignment = next(alignments)
            codon_alignments.append(alignment)
            self.assertTrue(
                np.array_equal(alignment.coordinates, np.array([[0, 256], [0, 768]]))
            )
            self.assertEqual(
                str(alignment),
                """\
gi|156879         0 M  S  F  T  L  T  N  K  N  V  I  F  V  A  G  L  G  G  I  G  
gi|156879         0 ATGTCGTTTACTTTGACCAACAAGAACGTGATTTTCGTGGCCGGTCTGGGAGGCATTGGT

gi|156879        20 L  D  T  S  K  E  L  L  K  R  D  L  K  N  L  V  I  L  D  R  
gi|156879        60 CTGGACACCAGCAAGGAGCTGCTCAAGCGCGATCTGAAGAACCTGGTGATCCTCGACCGC

gi|156879        40 I  E  N  P  A  A  I  A  E  L  K  A  I  N  P  K  V  T  V  T  
gi|156879       120 ATTGAGAACCCGGCTGCCATTGCCGAGCTGAAGGCAATCAATCCAAAGGTGACCGTCACC

gi|156879        60 F  Y  P  Y  D  V  T  V  P  I  A  E  T  T  K  L  L  K  T  I  
gi|156879       180 TTCTACCCCTATGATGTGACCGTGCCCATTGCCGAGACCACCAAGCTGCTGAAGACCATC

gi|156879        80 F  A  Q  L  K  T  V  D  V  L  I  N  G  A  G  I  L  D  D  H  
gi|156879       240 TTCGCCCAGCTGAAGACCGTCGATGTCCTGATCAACGGAGCTGGTATCCTGGACGATCAC

gi|156879       100 Q  I  E  R  T  I  A  V  N  Y  T  G  L  V  N  T  T  T  A  I  
gi|156879       300 CAGATCGAGCGCACCATTGCCGTCAACTACACTGGCCTGGTCAACACCACGACGGCCATT

gi|156879       120 L  D  F  W  D  K  R  K  G  G  P  G  G  I  I  C  N  I  G  S  
gi|156879       360 CTGGACTTCTGGGACAAGCGCAAGGGCGGTCCAGGTGGTATCATCTGCAACATTGGATCC

gi|156879       140 V  T  G  F  N  A  I  Y  Q  V  P  V  Y  S  G  T  K  A  A  V  
gi|156879       420 GTCACTGGATTCAATGCCATCTACCAGGTGCCCGTCTACTCCGGCACCAAGGCCGCCGTG

gi|156879       160 V  N  F  T  S  S  L  A  K  L  A  P  I  T  G  V  T  A  Y  T  
gi|156879       480 GTCAACTTCACCAGCTCCCTGGCGAAACTGGCCCCCATTACCGGCGTGACGGCTTACACT

gi|156879       180 V  N  P  G  I  T  R  T  T  L  V  H  T  F  N  S  W  L  D  V  
gi|156879       540 GTGAACCCCGGCATCACCCGCACCACCCTGGTGCACACGTTCAACTCCTGGTTGGATGTT

gi|156879       200 E  P  Q  V  A  E  K  L  L  A  H  P  T  Q  P  S  L  A  C  A  
gi|156879       600 GAGCCTCAGGTTGCCGAGAAGCTCCTGGCTCATCCCACCCAGCCCTCGTTGGCCTGCGCC

gi|156879       220 E  N  F  V  K  A  I  E  L  N  Q  N  G  A  I  W  K  L  D  L  
gi|156879       660 GAGAACTTCGTCAAGGCTATCGAGCTGAACCAGAACGGAGCCATCTGGAAACTGGACTTG

gi|156879       240 G  T  L  E  A  I  Q  W  T  K  H  W  D  S  G  I   256
gi|156879       720 GGCACCCTGGAGGCCATCCAGTGGACCAAGCACTGGGACTCCGGCATC 768
""",
            )
            protein_record = protein_alignment.sequences[17]
            nucleotide_record = nucleotide_records[protein_record.id]
            self.assertEqual(nucleotide_record.id, protein_record.id)
            alignments = aligner.align(protein_record, nucleotide_record)
            self.assertEqual(len(alignments), 1)
            alignment = next(alignments)
            codon_alignments.append(alignment)
            self.assertTrue(
                np.array_equal(alignment.coordinates, np.array([[0, 256], [0, 768]]))
            )
            self.assertEqual(
                str(alignment),
                """\
gi|156877         0 M  S  F  T  L  T  N  K  N  V  I  F  V  A  G  L  G  G  I  G  
gi|156877         0 ATGTCGTTTACTTTGACCAACAAGAACGTGATTTTCGTGGCCGGTCTGGGAGGCATTGGT

gi|156877        20 L  D  T  S  K  E  L  L  K  R  D  L  K  N  L  V  I  L  D  R  
gi|156877        60 CTGGACACCAGCAAGGAGCTGCTCAAGCGCGATCTGAAGAACCTGGTGATCCTCGACCGC

gi|156877        40 I  E  N  P  A  A  I  A  E  L  K  A  I  N  P  K  V  T  V  T  
gi|156877       120 ATTGAGAACCCGGCTGCCATTGCCGAGCTGAAGGCAATCAATCCAAAGGTGACCGTCACC

gi|156877        60 F  Y  P  Y  D  V  T  V  P  I  A  E  T  T  K  L  L  K  T  I  
gi|156877       180 TTCTACCCCTATGATGTGACCGTGCCCATTGCCGAGACCACCAAGCTGCTGAAGACCATC

gi|156877        80 F  A  Q  L  K  T  V  D  V  L  I  N  G  A  G  I  L  D  D  H  
gi|156877       240 TTCGCCCAGCTGAAGACCGTCGATGTCCTGATCAACGGAGCTGGTATCCTGGACGATCAC

gi|156877       100 Q  I  E  R  T  I  A  V  N  Y  T  G  L  V  N  T  T  T  A  I  
gi|156877       300 CAGATCGAGCGCACCATTGCCGTCAACTACACTGGCCTGGTCAACACCACGACGGCCATT

gi|156877       120 L  D  F  W  D  K  R  K  G  G  P  G  G  I  I  C  N  I  G  S  
gi|156877       360 CTGGACTTCTGGGACAAGCGCAAGGGCGGTCCCGGTGGTATCATCTGCAACATTGGATCC

gi|156877       140 V  T  G  F  N  A  I  Y  Q  V  P  V  Y  S  G  T  K  A  A  V  
gi|156877       420 GTCACTGGATTCAATGCCATCTACCAGGTGCCCGTCTACTCCGGCACCAAGGCCGCCGTG

gi|156877       160 V  N  F  T  S  S  L  A  K  L  A  P  I  T  G  V  T  A  Y  T  
gi|156877       480 GTCAACTTCACCAGCTCCCTGGCGAAACTGGCCCCCATTACCGGCGTGACGGCTTACACT

gi|156877       180 V  N  P  G  I  T  R  T  T  L  V  H  T  F  N  S  W  L  D  V  
gi|156877       540 GTGAACCCCGGCATCACCCGCACCACCCTGGTGCACACGTTCAACTCCTGGTTGGATGTT

gi|156877       200 E  P  Q  V  A  E  K  L  L  A  H  P  T  Q  P  S  L  A  C  A  
gi|156877       600 GAGCCTCAGGTTGCCGAGAAGCTCCTGGCTCATCCCACCCAGCCCTCGTTGGCCTGCGCC

gi|156877       220 E  N  F  V  K  A  I  E  L  N  Q  N  G  A  I  W  K  L  D  L  
gi|156877       660 GAGAACTTCGTCAAGGCTATCGAGCTGAACCAGAACGGAGCCATCTGGAAACTGGACTTG

gi|156877       240 G  T  L  E  A  I  Q  W  T  K  H  W  D  S  G  I   256
gi|156877       720 GGCACCCTGGAGGCCATCCAGTGGACCAAGCACTGGGACTCCGGCATC 768
""",
            )
            protein_record = protein_alignment.sequences[18]
            nucleotide_record = nucleotide_records[protein_record.id]
            self.assertEqual(nucleotide_record.id, protein_record.id)
            alignments = aligner.align(protein_record, nucleotide_record)
            self.assertEqual(len(alignments), 1)
            alignment = next(alignments)
            codon_alignments.append(alignment)
            self.assertTrue(
                np.array_equal(alignment.coordinates, np.array([[0, 256], [0, 768]]))
            )
            self.assertEqual(
                str(alignment),
                """\
gi|156875         0 M  S  F  T  L  T  N  K  N  V  I  F  V  A  G  L  G  G  I  G  
gi|156875         0 ATGTCGTTTACTTTGACCAACAAGAACGTGATTTTCGTGGCCGGTCTGGGAGGCATTGGT

gi|156875        20 L  D  T  S  K  E  L  L  K  R  D  L  K  N  L  V  I  L  D  R  
gi|156875        60 CTGGACACCAGCAAGGAGCTGCTCAAGCGCGATCTGAAGAACCTGGTGATCCTCGACCGC

gi|156875        40 I  E  N  P  A  A  I  A  E  L  K  A  I  N  P  K  V  T  V  T  
gi|156875       120 ATTGAGAACCCGGCTGCCATTGCCGAGCTGAAGGCAATCAATCCAAAGGTGACCGTCACC

gi|156875        60 F  Y  P  Y  D  V  T  V  P  I  A  E  T  T  K  L  L  K  T  I  
gi|156875       180 TTCTACCCCTATGATGTGACCGTGCCCATTGCCGAGACCACCAAGCTGCTGAAGACCATC

gi|156875        80 F  A  Q  L  K  T  V  D  V  L  I  N  G  A  G  I  L  D  D  H  
gi|156875       240 TTCGCCCAGCTGAAGACCGTCGATGTCCTGATCAACGGAGCTGGTATCCTGGACGATCAC

gi|156875       100 Q  I  E  R  T  I  A  V  N  Y  T  G  L  V  N  T  T  T  A  I  
gi|156875       300 CAGATCGAGCGCACCATTGCCGTCAACTACACTGGCCTGGTCAACACCACGACGGCCATT

gi|156875       120 L  D  F  W  D  K  R  K  G  G  P  G  G  I  I  C  N  I  G  S  
gi|156875       360 CTGGACTTCTGGGACAAGCGCAAGGGCGGTCCCGGTGGTATCATCTGCAACATTGGATCC

gi|156875       140 V  T  G  F  N  A  I  Y  Q  V  P  V  Y  S  G  T  K  A  A  V  
gi|156875       420 GTCACTGGATTCAATGCCATCTACCAGGTGCCCGTCTACTCCGGCACCAAGGCCGCCGTG

gi|156875       160 V  N  F  T  S  S  L  A  K  L  A  P  I  T  G  V  T  A  Y  T  
gi|156875       480 GTCAACTTCACCAGCTCCCTGGCGAAACTGGCCCCCATTACCGGCGTGACGGCTTACACT

gi|156875       180 V  N  P  G  I  T  R  T  T  L  V  H  T  F  N  S  W  L  D  V  
gi|156875       540 GTGAACCCCGGCATCACCCGCACCACCCTGGTGCACACGTTCAACTCCTGGTTGGATGTT

gi|156875       200 E  P  Q  V  A  E  K  L  L  A  H  P  T  Q  P  S  L  A  C  A  
gi|156875       600 GAGCCTCAGGTTGCCGAGAAGCTCCTGGCTCATCCCACCCAGCCCTCGTTGGCCTGCGCC

gi|156875       220 E  N  F  V  K  A  I  E  L  N  Q  N  G  A  I  W  K  L  D  L  
gi|156875       660 GAGAACTTCGTCAAGGCTATCGAGCTGAACCAGAACGGAGCCATCTGGAAACTGGACTTG

gi|156875       240 G  T  L  E  A  I  Q  W  T  K  H  W  D  S  G  I   256
gi|156875       720 GGCACCCTGGAGGCCATCCAGTGGACCAAGCACTGGGACTCCGGCATC 768
""",
            )
            protein_record = protein_alignment.sequences[19]
            nucleotide_record = nucleotide_records[protein_record.id]
            self.assertEqual(nucleotide_record.id, protein_record.id)
            alignments = aligner.align(protein_record, nucleotide_record)
            self.assertEqual(len(alignments), 1)
            alignment = next(alignments)
            codon_alignments.append(alignment)
            self.assertTrue(
                np.array_equal(alignment.coordinates, np.array([[0, 256], [0, 768]]))
            )
            self.assertEqual(
                str(alignment),
                """\
gi|156873         0 M  S  F  T  L  T  N  K  N  V  I  F  V  A  G  L  G  G  I  G  
gi|156873         0 ATGTCGTTTACTTTGACCAACAAGAACGTGATTTTCGTGGCCGGTCTGGGAGGCATTGGT

gi|156873        20 L  D  T  S  K  E  L  L  K  R  D  L  K  N  L  V  I  L  D  R  
gi|156873        60 CTGGACACCAGCAAGGAGCTGCTCAAGCGCGATCTGAAGAACCTGGTGATCCTCGACCGC

gi|156873        40 I  E  N  P  A  A  I  A  E  L  K  A  I  N  P  K  V  T  V  T  
gi|156873       120 ATTGAGAACCCGGCTGCCATTGCCGAGCTGAAGGCAATCAATCCAAAGGTGACCGTCACC

gi|156873        60 F  Y  P  Y  D  V  T  V  P  I  A  E  T  T  K  L  L  K  T  I  
gi|156873       180 TTCTACCCCTATGATGTGACCGTGCCCATTGCCGAGACCACCAAGCTGCTGAAGACCATC

gi|156873        80 F  A  Q  L  K  T  V  D  V  L  I  N  G  A  G  I  L  D  D  H  
gi|156873       240 TTCGCCCAGCTGAAGACCGTCGATGTCCTGATCAACGGAGCTGGTATCCTGGACGATCAC

gi|156873       100 Q  I  E  R  T  I  A  V  N  Y  T  G  L  V  N  T  T  T  A  I  
gi|156873       300 CAGATCGAGCGCACCATTGCCGTCAACTACACTGGCCTGGTCAACACCACGACGGCCATT

gi|156873       120 L  D  F  W  D  K  R  K  G  G  P  G  G  I  I  C  N  I  G  S  
gi|156873       360 CTGGACTTCTGGGACAAGCGCAAGGGCGGTCCCGGTGGTATCATCTGCAACATTGGATCC

gi|156873       140 V  T  G  F  N  A  I  Y  Q  V  P  V  Y  S  G  T  K  A  A  V  
gi|156873       420 GTCACTGGATTCAATGCCATCTACCAGGTGCCCGTCTACTCCGGCACCAAGGCCGCCGTG

gi|156873       160 V  N  F  T  S  S  L  A  K  L  A  P  I  T  G  V  T  A  Y  T  
gi|156873       480 GTCAACTTCACCAGCTCCCTGGCGAAACTGGCCCCCATTACCGGCGTGACGGCTTACACT

gi|156873       180 V  N  P  G  I  T  R  T  T  L  V  H  T  F  N  S  W  L  D  V  
gi|156873       540 GTGAACCCCGGCATCACCCGCACCACCCTGGTGCACACGTTCAACTCCTGGTTGGATGTT

gi|156873       200 E  P  Q  V  A  E  K  L  L  A  H  P  T  Q  P  S  L  A  C  A  
gi|156873       600 GAGCCTCAGGTTGCCGAGAAGCTCCTGGCTCATCCCACCCAGCCCTCGTTGGCCTGCGCC

gi|156873       220 E  N  F  V  K  A  I  E  L  N  Q  N  G  A  I  W  K  L  D  L  
gi|156873       660 GAGAACTTCGTCAAGGCTATCGAGCTGAACCAGAACGGAGCCATCTGGAAACTGGACTTG

gi|156873       240 G  T  L  E  A  I  Q  W  T  K  H  W  D  S  G  I   256
gi|156873       720 GGCACCCTGGAGGCCATCCAGTGGACCAAGCACTGGGACTCCGGCATC 768
""",
            )
            protein_record = protein_alignment.sequences[20]
            nucleotide_record = nucleotide_records[protein_record.id]
            self.assertEqual(nucleotide_record.id, protein_record.id)
            alignments = aligner.align(protein_record, nucleotide_record)
            self.assertEqual(len(alignments), 1)
            alignment = next(alignments)
            codon_alignments.append(alignment)
            self.assertTrue(
                np.array_equal(alignment.coordinates, np.array([[0, 256], [0, 768]]))
            )
            self.assertEqual(
                str(alignment),
                """\
gi|156871         0 M  S  F  T  L  T  N  K  N  V  I  F  V  A  G  L  G  G  I  G  
gi|156871         0 ATGTCGTTTACTTTGACCAACAAGAACGTGATTTTCGTGGCCGGTCTGGGAGGCATTGGT

gi|156871        20 L  D  T  S  K  E  L  L  K  R  D  L  K  N  L  V  I  L  D  R  
gi|156871        60 CTGGACACCAGCAAGGAGCTGCTCAAGCGCGATCTGAAGAACCTGGTGATCCTCGACCGC

gi|156871        40 I  E  N  P  A  A  I  A  E  L  K  A  I  N  P  K  V  T  V  T  
gi|156871       120 ATTGAGAACCCGGCTGCCATTGCCGAGCTGAAGGCAATCAATCCAAAGGTGACCGTCACC

gi|156871        60 F  Y  P  Y  D  V  T  V  P  I  A  E  T  T  K  L  L  K  T  I  
gi|156871       180 TTCTACCCCTATGATGTGACCGTGCCCATTGCCGAGACCACCAAGCTGCTGAAGACCATC

gi|156871        80 F  A  Q  L  K  T  V  D  V  L  I  N  G  A  G  I  L  D  D  H  
gi|156871       240 TTCGCCCAGCTGAAGACCGTCGATGTCCTGATCAACGGAGCTGGTATCCTGGACGATCAC

gi|156871       100 Q  I  E  R  T  I  A  V  N  Y  T  G  L  V  N  T  T  T  A  I  
gi|156871       300 CAGATCGAGCGCACCATTGCCGTCAACTACACTGGCCTGGTCAACACCACGACGGCCATT

gi|156871       120 L  D  F  W  D  K  R  K  G  G  P  G  G  I  I  C  N  I  G  S  
gi|156871       360 CTGGACTTCTGGGACAAGCGCAAGGGCGGTCCCGGTGGTATCATCTGCAACATTGGATCC

gi|156871       140 V  T  G  F  N  A  I  Y  Q  V  P  V  Y  S  G  T  K  A  A  V  
gi|156871       420 GTCACTGGATTCAATGCCATCTACCAGGTGCCCGTCTACTCCGGCACCAAGGCCGCCGTG

gi|156871       160 V  N  F  T  S  S  L  A  K  L  A  P  I  T  G  V  T  A  Y  T  
gi|156871       480 GTCAACTTCACCAGCTCCCTGGCGAAACTGGCCCCCATTACCGGCGTGACGGCTTACACT

gi|156871       180 V  N  P  G  I  T  R  T  T  L  V  H  T  F  N  S  W  L  D  V  
gi|156871       540 GTGAACCCCGGCATCACCCGCACCACCCTGGTGCACACGTTCAACTCCTGGTTGGATGTT

gi|156871       200 E  P  Q  V  A  E  K  L  L  A  H  P  T  Q  P  S  L  A  C  A  
gi|156871       600 GAGCCTCAGGTTGCCGAGAAGCTCCTGGCTCATCCCACCCAGCCCTCGTTGGCCTGCGCC

gi|156871       220 E  N  F  V  K  A  I  E  L  N  Q  N  G  A  I  W  K  L  D  L  
gi|156871       660 GAGAACTTCGTCAAGGCTATCGAGCTGAACCAGAACGGAGCCATCTGGAAACTGGACTTG

gi|156871       240 G  T  L  E  A  I  Q  W  T  K  H  W  D  S  G  I   256
gi|156871       720 GGCACCCTGGAGGCCATCCAGTGGACCAAGCACTGGGACTCCGGCATC 768
""",
            )
            protein_record = protein_alignment.sequences[21]
            nucleotide_record = nucleotide_records[protein_record.id]
            self.assertEqual(nucleotide_record.id, protein_record.id)
            alignments = aligner.align(protein_record, nucleotide_record)
            self.assertEqual(len(alignments), 1)
            alignment = next(alignments)
            codon_alignments.append(alignment)
            self.assertTrue(
                np.array_equal(alignment.coordinates, np.array([[0, 256], [0, 768]]))
            )
            self.assertEqual(
                str(alignment),
                """\
gi|156863         0 M  S  F  T  L  T  N  K  N  V  I  F  V  A  G  L  G  G  I  G  
gi|156863         0 ATGTCGTTTACTTTGACCAACAAGAACGTGATTTTCGTTGCCGGTCTGGGAGGCATTGGT

gi|156863        20 L  D  T  S  K  E  L  L  K  R  D  L  K  N  L  V  I  L  D  R  
gi|156863        60 CTGGACACCAGCAAGGAGCTGCTCAAGCGCGATCTGAAGAACCTGGTGATCCTCGACCGC

gi|156863        40 I  E  N  P  A  A  I  A  E  L  K  A  I  N  P  K  V  T  V  T  
gi|156863       120 ATTGAGAACCCGGCTGCCATTGCCGAGCTGAAGGCAATCAATCCAAAGGTGACCGTCACC

gi|156863        60 F  Y  P  Y  D  V  T  V  P  I  A  E  T  T  K  L  L  K  T  I  
gi|156863       180 TTCTACCCCTATGATGTGACCGTGCCCATTGCCGAGACCACCAAGCTGCTGAAGACCATC

gi|156863        80 F  A  Q  L  K  T  V  D  V  L  I  N  G  A  G  I  L  D  D  H  
gi|156863       240 TTCGCCCAGCTGAAGACCGTCGATGTCCTGATCAACGGAGCTGGTATCCTGGACGATCAC

gi|156863       100 Q  I  E  R  T  I  A  V  N  Y  T  G  L  V  N  T  T  T  A  I  
gi|156863       300 CAGATCGAGCGCACCATTGCCGTCAACTACACTGGCCTGGTCAACACCACGACGGCCATT

gi|156863       120 L  D  F  W  D  K  R  K  G  G  P  G  G  I  I  C  N  I  G  S  
gi|156863       360 CTGGACTTCTGGGACAAGCGCAAGGGCGGTCCCGGTGGTATCATCTGCAACATTGGATCC

gi|156863       140 V  T  G  F  N  A  I  Y  Q  V  P  V  Y  S  G  T  K  A  A  V  
gi|156863       420 GTCACTGGATTCAATGCCATCTACCAGGTGCCCGTCTACTCCGGCACCAAGGCCGCCGTG

gi|156863       160 V  N  F  T  S  S  L  A  K  L  A  P  I  T  G  V  T  A  Y  T  
gi|156863       480 GTCAACTTCACCAGCTCCCTGGCGAAACTGGCCCCCATTACCGGCGTGACCGCTTACACC

gi|156863       180 V  N  P  G  I  T  R  T  T  L  V  H  K  F  N  S  W  L  D  V  
gi|156863       540 GTGAACCCCGGCATCACCCGCACCACCCTGGTGCACAAGTTCAACTCCTGGTTGGATGTT

gi|156863       200 E  P  Q  V  A  E  K  L  L  A  H  P  T  Q  P  S  L  A  C  A  
gi|156863       600 GAGCCCCAGGTTGCTGAGAAGCTCCTGGCTCATCCCACCCAGCCATCGTTGGCCTGCGCC

gi|156863       220 E  N  F  V  K  A  I  E  L  N  Q  N  G  A  I  W  K  L  D  L  
gi|156863       660 GAGAACTTCGTCAAGGCTATCGAACTGAACCAGAACGGAGCCATCTGGAAACTGGACTTG

gi|156863       240 G  T  L  E  A  I  Q  W  T  K  H  W  D  S  G  I   256
gi|156863       720 GGCACCCTGGAGGCCATCCAGTGGACCAAGCACTGGGACTCCGGCATC 768
""",
            )
            protein_record = protein_alignment.sequences[22]
            nucleotide_record = nucleotide_records[protein_record.id]
            self.assertEqual(nucleotide_record.id, protein_record.id)
            alignments = aligner.align(protein_record, nucleotide_record)
            self.assertEqual(len(alignments), 1)
            alignment = next(alignments)
            codon_alignments.append(alignment)
            self.assertTrue(
                np.array_equal(alignment.coordinates, np.array([[0, 256], [0, 768]]))
            )
            self.assertEqual(
                str(alignment),
                """\
gi|156869         0 M  S  F  T  L  T  N  K  N  V  I  F  V  A  G  L  G  G  I  G  
gi|156869         0 ATGTCGTTTACTTTGACCAACAAGAACGTGATTTTCGTGGCCGGTCTGGGAGGCATTGGT

gi|156869        20 L  D  T  S  K  E  L  L  K  R  D  L  K  N  L  V  I  L  D  R  
gi|156869        60 CTGGACACCAGCAAGGAGCTGCTCAAGCGCGATCTGAAGAACCTGGTGATCCTCGACCGC

gi|156869        40 I  E  N  P  A  A  I  A  E  L  K  A  I  N  P  K  V  T  V  T  
gi|156869       120 ATTGAGAACCCGGCTGCCATTGCCGAGCTGAAGGCAATCAATCCAAAGGTGACCGTCACC

gi|156869        60 F  Y  P  Y  D  V  T  V  P  I  A  E  T  T  K  L  L  K  T  I  
gi|156869       180 TTCTACCCCTATGATGTGACCGTGCCCATTGCCGAGACCACCAAGCTGCTGAAGACCATC

gi|156869        80 F  A  Q  L  K  T  V  D  V  L  I  N  G  A  G  I  L  D  D  H  
gi|156869       240 TTCGCCCAGCTGAAGACCGTCGATGTCCTGATCAACGGAGCTGGTATCCTGGACGATCAC

gi|156869       100 Q  I  E  R  T  I  A  V  N  Y  T  G  L  V  N  T  T  T  A  I  
gi|156869       300 CAGATCGAGCGCACCATTGCCGTCAACTACACTGGCCTGGTCAACACCACGACGGCCATT

gi|156869       120 L  D  F  W  D  K  R  K  G  G  P  G  G  I  I  C  N  I  G  S  
gi|156869       360 CTGGACTTCTGGGACAAGCGCAAGGGCGGTCCCGGTGGTATCATCTGCAACATTGGATCC

gi|156869       140 V  T  G  F  N  A  I  Y  Q  V  P  V  Y  S  G  T  K  A  A  V  
gi|156869       420 GTCACTGGATTCAATGCCATCTACCAGGTGCCCGTCTACTCCGGCACCAAGGCCGCCGTG

gi|156869       160 V  N  F  T  S  S  L  A  K  L  A  P  I  T  G  V  T  A  Y  T  
gi|156869       480 GTCAACTTCACCAGCTCCCTGGCGAAACTGGCCCCCATTACCGGCGTGACCGCTTACACT

gi|156869       180 V  N  P  G  I  T  R  T  T  L  V  H  K  F  N  S  W  L  D  V  
gi|156869       540 GTGAACCCCGGCATCACCCGCACCACCCTGGTGCACAAGTTCAACTCCTGGTTGGATGTT

gi|156869       200 E  P  Q  V  A  E  K  L  L  A  H  P  T  Q  P  S  L  A  C  A  
gi|156869       600 GAGCCTCAGGTTGCTGAGAAGCTCCTGGCTCATCCCACCCAGCCCTCGTTGGCCTGCGCC

gi|156869       220 E  N  F  V  K  A  I  E  L  N  Q  N  G  A  I  W  K  L  D  L  
gi|156869       660 GAGAACTTCGTCAAGGCTATCGAACTGAACCAGAACGGAGCCATCTGGAAACTGGACTTG

gi|156869       240 G  T  L  E  A  I  Q  W  T  K  H  W  D  S  G  I   256
gi|156869       720 GGCACCCTGGAGGCCATCCAGTGGACCAAGCACTGGGACTCCGGCATC 768
""",
            )
            protein_record = protein_alignment.sequences[23]
            nucleotide_record = nucleotide_records[protein_record.id]
            self.assertEqual(nucleotide_record.id, protein_record.id)
            alignments = aligner.align(protein_record, nucleotide_record)
            self.assertEqual(len(alignments), 1)
            alignment = next(alignments)
            codon_alignments.append(alignment)
            self.assertTrue(
                np.array_equal(alignment.coordinates, np.array([[0, 256], [0, 768]]))
            )
            self.assertEqual(
                str(alignment),
                """\
gi|156867         0 M  S  F  T  L  T  N  K  N  V  I  F  V  A  G  L  G  G  I  G  
gi|156867         0 ATGTCGTTTACTTTGACCAACAAGAACGTGATTTTCGTGGCCGGTCTGGGAGGCATTGGT

gi|156867        20 L  D  T  S  K  E  L  L  K  R  D  L  K  N  L  V  I  L  D  R  
gi|156867        60 CTGGACACCAGCAAGGAGCTGCTCAAGCGCGATCTGAAGAACCTGGTGATCCTCGACCGC

gi|156867        40 I  E  N  P  A  A  I  A  E  L  K  A  I  N  P  K  V  T  V  T  
gi|156867       120 ATTGAGAACCCGGCTGCCATTGCCGAGCTGAAGGCAATCAATCCAAAGGTGACCGTCACC

gi|156867        60 F  Y  P  Y  D  V  T  V  P  I  A  E  T  T  K  L  L  K  T  I  
gi|156867       180 TTCTACCCCTATGATGTGACCGTGCCCATTGCCGAGACCACCAAGCTGCTGAAGACCATC

gi|156867        80 F  A  Q  L  K  T  V  D  V  L  I  N  G  A  G  I  L  D  D  H  
gi|156867       240 TTCGCCCAGCTGAAGACCGTCGATGTCCTGATCAACGGAGCTGGTATCCTGGACGATCAC

gi|156867       100 Q  I  E  R  T  I  A  V  N  Y  T  G  L  V  N  T  T  T  A  I  
gi|156867       300 CAGATCGAGCGCACCATTGCCGTCAACTACACTGGCCTGGTCAACACCACGACGGCCATT

gi|156867       120 L  D  F  W  D  K  R  K  G  G  P  G  G  I  I  C  N  I  G  S  
gi|156867       360 CTGGACTTCTGGGACAAGCGCAAGGGCGGTCCCGGTGGTATCATCTGCAACATTGGATCC

gi|156867       140 V  T  G  F  N  A  I  Y  Q  V  P  V  Y  S  G  T  K  A  A  V  
gi|156867       420 GTCACTGGATTCAATGCCATCTACCAGGTGCCCGTCTACTCCGGCACCAAGGCCGCCGTG

gi|156867       160 V  N  F  T  S  S  L  A  K  L  A  P  I  T  G  V  T  A  Y  T  
gi|156867       480 GTCAACTTCACCAGCTCCCTGGCGAAACTGGCCCCCATTACCGGCGTGACCGCTTACACC

gi|156867       180 V  N  P  G  I  T  R  T  T  L  V  H  K  F  N  S  W  L  D  V  
gi|156867       540 GTGAACCCCGGCATCACCCGCACCACCCTGGTGCACAAGTTCAACTCCTGGTTGGATGTT

gi|156867       200 E  P  Q  V  A  E  K  L  L  A  H  P  T  Q  P  S  L  A  C  A  
gi|156867       600 GAGCCCCAGGTTGCTGAGAAGCTCCTGGCTCATCCCACCCAGCCATCGTTGGCCTGCGCC

gi|156867       220 E  N  F  V  K  A  I  E  L  N  Q  N  G  A  I  W  K  L  D  L  
gi|156867       660 GAGAACTTCGTCAAGGCTATCGAGCTGAACCAGAACGGAGCCATCTGGAAACTGGACTTG

gi|156867       240 G  T  L  E  A  I  Q  W  T  K  H  W  D  S  G  I   256
gi|156867       720 GGCACCCTGGAGGCCATCCAGTGGACCAAGCACTGGGACTCCGGCATC 768
""",
            )
            protein_record = protein_alignment.sequences[24]
            nucleotide_record = nucleotide_records[protein_record.id]
            self.assertEqual(nucleotide_record.id, protein_record.id)
            alignments = aligner.align(protein_record, nucleotide_record)
            self.assertEqual(len(alignments), 1)
            alignment = next(alignments)
            codon_alignments.append(alignment)
            self.assertTrue(
                np.array_equal(alignment.coordinates, np.array([[0, 256], [0, 768]]))
            )
            self.assertEqual(
                str(alignment),
                """\
gi|156865         0 M  S  F  T  L  T  N  K  N  V  I  F  V  A  G  L  G  G  I  G  
gi|156865         0 ATGTCGTTTACTTTGACCAACAAGAACGTGATTTTCGTTGCCGGTCTGGGAGGCATTGGT

gi|156865        20 L  D  T  S  K  E  L  L  K  R  D  L  K  N  L  V  I  L  D  R  
gi|156865        60 CTGGACACCAGCAAGGAGCTGCTCAAGCGCGATCTGAAGAACCTGGTGATCCTCGACCGC

gi|156865        40 I  E  N  P  A  A  I  A  E  L  K  A  I  N  P  K  V  T  V  T  
gi|156865       120 ATTGAGAACCCGGCTGCCATTGCCGAGCTGAAGGCAATCAATCCAAAGGTGACCGTCACC

gi|156865        60 F  Y  P  Y  D  V  T  V  P  I  A  E  T  T  K  L  L  K  T  I  
gi|156865       180 TTCTACCCCTATGATGTGACCGTGCCCATTGCCGAGACCACCAAGCTGCTGAAGACCATC

gi|156865        80 F  A  Q  L  K  T  V  D  V  L  I  N  G  A  G  I  L  D  D  H  
gi|156865       240 TTCGCCCAGCTGAAGACCGTCGATGTCCTGATCAACGGAGCTGGTATCCTGGACGATCAC

gi|156865       100 Q  I  E  R  T  I  A  V  N  Y  T  G  L  V  N  T  T  T  A  I  
gi|156865       300 CAGATCGAGCGCACCATTGCCGTCAACTACACTGGCCTGGTCAACACCACGACGGCCATT

gi|156865       120 L  D  F  W  D  K  R  K  G  G  P  G  G  I  I  C  N  I  G  S  
gi|156865       360 CTGGACTTCTGGGACAAGCGCAAGGGCGGTCCCGGTGGTATCATCTGCAACATTGGATCC

gi|156865       140 V  T  G  F  N  A  I  Y  Q  V  P  V  Y  S  G  T  K  A  A  V  
gi|156865       420 GTCACTGGATTCAATGCCATCTACCAGGTGCCCGTCTACTCCGGCACCAAGGCCGCCGTG

gi|156865       160 V  N  F  T  S  S  L  A  K  L  A  P  I  T  G  V  T  A  Y  T  
gi|156865       480 GTCAACTTCACCAGCTCCCTGGCGAAACTGGCCCCCATTACCGGCGTGACCGCTTACACC

gi|156865       180 V  N  P  G  I  T  R  T  T  L  V  H  K  F  N  S  W  L  D  V  
gi|156865       540 GTGAACCCCGGCATCACCCGCACCACCCTGGTGCACAAGTTCAACTCCTGGTTGGATGTT

gi|156865       200 E  P  Q  V  A  E  K  L  L  A  H  P  T  Q  P  S  L  A  C  A  
gi|156865       600 GAGCCCCAGGTTGCTGAGAAGCTCCTGGCTCATCCCACCCAGCCATCGTTGGCCTGCGCC

gi|156865       220 E  N  F  V  K  A  I  E  L  N  Q  N  G  A  I  W  K  L  D  L  
gi|156865       660 GAGAACTTCGTCAAGGCTATCGAACTGAACCAGAACGGAGCCATCTGGAAACTGGACTTG

gi|156865       240 G  T  L  E  A  I  Q  W  T  K  H  W  D  S  G  I   256
gi|156865       720 GGCACCCTGGAGGCCATCCAGTGGACCAAGCACTGGGACTCCGGCATC 768
""",
            )
            protein_record = protein_alignment.sequences[25]
            nucleotide_record = nucleotide_records[protein_record.id]
            self.assertEqual(nucleotide_record.id, protein_record.id)
            alignments = aligner.align(protein_record, nucleotide_record)
            self.assertEqual(len(alignments), 1)
            alignment = next(alignments)
            codon_alignments.append(alignment)
            self.assertTrue(
                np.array_equal(alignment.coordinates, np.array([[0, 256], [0, 768]]))
            )
            self.assertEqual(
                str(alignment),
                """\
gi|156861         0 M  S  F  T  L  T  N  K  N  V  I  F  V  A  G  L  G  G  I  G  
gi|156861         0 ATGTCGTTTACTTTGACCAACAAGAACGTGATTTTCGTTGCCGGTCTGGGAGGCATTGGT

gi|156861        20 L  D  T  S  K  E  L  L  K  R  D  L  K  N  L  V  I  L  D  R  
gi|156861        60 CTGGACACCAGCAAGGAGCTGCTCAAGCGCGATCTGAAGAACCTGGTGATCCTCGACCGC

gi|156861        40 I  E  N  P  A  A  I  A  E  L  K  A  I  N  P  K  V  T  V  T  
gi|156861       120 ATTGAGAACCCGGCTGCCATTGCCGAGCTGAAGGCAATCAATCCAAAGGTGACCGTCACC

gi|156861        60 F  Y  P  Y  D  V  T  V  P  I  A  E  T  T  K  L  L  K  T  I  
gi|156861       180 TTCTACCCCTATGATGTGACCGTGCCCATTGCCGAGACCACCAAGTTGCTGAAGACCATC

gi|156861        80 F  A  Q  L  K  T  V  D  V  L  I  N  G  A  G  I  L  D  D  H  
gi|156861       240 TTCGCCCAGCTGAAGACCGTCGATGTCCTGATCAACGGAGCTGGTATCCTGGACGATCAC

gi|156861       100 Q  I  E  R  T  I  A  V  N  Y  T  G  L  V  N  T  T  T  A  I  
gi|156861       300 CAGATCGAGCGCACCATTGCCGTCAACTACACTGGCCTGGTCAACACCACGACGGCTATT

gi|156861       120 L  D  F  W  D  K  R  K  G  G  P  G  G  I  I  C  N  I  G  S  
gi|156861       360 CTGGACTTCTGGGACAAGCGCAAGGGTGGTCCCGGTGGTATCATCTGCAACATTGGATCC

gi|156861       140 V  T  G  F  N  A  I  Y  Q  V  P  V  Y  S  G  T  K  A  A  V  
gi|156861       420 GTCACTGGATTCAATGCCATATACCAGGTGCCCGTCTACTCCGGCACCAAGGCCGCCGTG

gi|156861       160 V  N  F  T  S  S  L  A  K  L  A  P  I  T  G  V  T  A  Y  T  
gi|156861       480 GTCAACTTCACCAGCTCCCTGGCGAAACTGGCACCCATCACCGGCGTGACCGCTTACACC

gi|156861       180 V  N  P  G  I  T  R  T  T  L  V  H  K  F  N  S  W  L  D  V  
gi|156861       540 GTGAACCCCGGCATCACCCGCACCACCCTGGTGCACAAGTTCAACTCCTGGTTGGATGTT

gi|156861       200 E  P  Q  V  A  E  K  L  L  A  H  P  T  Q  P  S  L  A  C  A  
gi|156861       600 GAGCCCCAGGTTGCTGAGAAGCTCCTGGCTCATCCCACCCAGCCATCGTTGGCCTGCGCC

gi|156861       220 E  N  F  V  K  A  I  E  L  N  Q  N  G  A  I  W  K  L  D  L  
gi|156861       660 GAGAACTTCGTCAAGGCTATCGAGCTGAACCAGAACGGAGCCATCTGGAAACTGGACTTG

gi|156861       240 G  T  L  E  A  I  Q  W  T  K  H  W  D  S  G  I   256
gi|156861       720 GGCACCCTGGAGGCCATCCAGTGGACCAAGCACTGGGACTCCGGCATC 768
""",
            )
            protein_record = protein_alignment.sequences[26]
            nucleotide_record = nucleotide_records[protein_record.id]
            self.assertEqual(nucleotide_record.id, protein_record.id)
            alignments = aligner.align(protein_record, nucleotide_record)
            self.assertEqual(len(alignments), 1)
            alignment = next(alignments)
            codon_alignments.append(alignment)
            self.assertTrue(
                np.array_equal(alignment.coordinates, np.array([[0, 256], [0, 768]]))
            )
            self.assertEqual(
                str(alignment),
                """\
gi|156859         0 M  S  F  T  L  T  N  K  N  V  I  F  V  A  G  L  G  G  I  G  
gi|156859         0 ATGTCGTTTACTTTGACCAACAAGAACGTGATTTTCGTTGCCGGTCTGGGAGGCATTGGT

gi|156859        20 L  D  T  S  K  E  L  L  K  R  D  L  K  N  L  V  I  L  D  R  
gi|156859        60 CTGGACACCAGCAAGGAGCTGCTCAAGCGCGATCTGAAGAACCTGGTGATCCTCGACCGC

gi|156859        40 I  E  N  P  A  A  I  A  E  L  K  A  I  N  P  K  V  T  V  T  
gi|156859       120 ATTGAGAACCCGGCTGCCATTGCCGAGCTGAAGGCAATCAATCCAAAGGTGACCGTCACC

gi|156859        60 F  Y  P  Y  D  V  T  V  P  I  A  E  T  T  K  L  L  K  T  I  
gi|156859       180 TTCTACCCCTATGATGTGACCGTGCCCATTGCCGAGACCACCAAGTTGCTGAAGACCATC

gi|156859        80 F  A  Q  L  K  T  V  D  V  L  I  N  G  A  G  I  L  D  D  H  
gi|156859       240 TTCGCCCAGCTGAAGACCGTCGATGTCCTGATCAACGGAGCTGGTATCCTGGACGATCAC

gi|156859       100 Q  I  E  R  T  I  A  V  N  Y  T  G  L  V  N  T  T  T  A  I  
gi|156859       300 CAGATCGAGCGCACCATTGCCGTCAACTACACTGGCCTGGTCAACACCACGACGGCCATT

gi|156859       120 L  D  F  W  D  K  R  K  G  G  P  G  G  I  I  C  N  I  G  S  
gi|156859       360 CTGGACTTCTGGGACAAGCGCAAGGGTGGTCCCGGTGGTATCATCTGCAACATTGGATCC

gi|156859       140 V  T  G  F  N  A  I  Y  Q  V  P  V  Y  S  G  T  K  A  A  V  
gi|156859       420 GTCACTGGATTCAATGCCATATACCAGGTGCCCGTCTACTCCGGCACCAAGGCCGCCGTG

gi|156859       160 V  N  F  T  S  S  L  A  K  L  A  P  I  T  G  V  T  A  Y  T  
gi|156859       480 GTCAACTTCACCAGCTCCCTGGCGAAACTGGCACCCATCACCGGCGTGACCGCTTACACC

gi|156859       180 V  N  P  G  I  T  R  T  T  L  V  H  K  F  N  S  W  L  D  V  
gi|156859       540 GTGAACCCCGGCATCACCCGCACCACCCTGGTGCACAAGTTCAACTCCTGGTTGGATGTT

gi|156859       200 E  P  Q  V  A  E  K  L  L  A  H  P  T  Q  P  S  L  A  C  A  
gi|156859       600 GAGCCCCAGGTTGCTGAGAAGCTCCTGGCTCATCCCACCCAGCCATCGTTGGCCTGCGCC

gi|156859       220 E  N  F  V  K  A  I  E  L  N  Q  N  G  A  I  W  K  L  D  L  
gi|156859       660 GAGAACTTCGTCAAGGCTATCGAGCTGAACCAGAACGGAGCCATCTGGAAACTGGACTTG

gi|156859       240 G  T  L  E  A  I  Q  W  T  K  H  W  D  S  G  I   256
gi|156859       720 GGCACCCTGGAGGCCATCCAGTGGACCAAGCACTGGGACTCCGGCATC 768
""",
            )
            nucleotide_records.close()  # Close indexed FASTA file
            alignment = protein_alignment.mapall(codon_alignments)
            unique_species = [
                "Drosophila simulans",
                "Drosophila yakuba",
                "D.melanogaster",
            ]
            species = []
            for record in alignment.sequences:
                description = record.description
                for s in unique_species:
                    if s in description:
                        break
                else:
                    raise Exception(f"Failed to find species for {description}")
                species.append(s)
            pvalue = mktest(alignment, species)
            self.assertAlmostEqual(pvalue, 0.0020645725725430097)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
