# Copyright (C) 2013 by Zheng Ruan (zruan1991@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for the CodonAlign modules.
"""
import sys
import warnings
import tempfile
import platform
import unittest

from Bio import BiopythonExperimentalWarning

import warnings
with warnings.catch_warnings():
   warnings.simplefilter('ignore', BiopythonExperimentalWarning)
   from Bio import CodonAlign, SeqIO, AlignIO

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Align import MultipleSeqAlignment


TEST_ALIGN_FILE1 = [('CodonAlign/nucl1.fa', 'CodonAlign/pro1.aln'), 'parse']
TEST_ALIGN_FILE2 = [('CodonAlign/nucl2.fa', 'CodonAlign/pro2.aln'), 'parse']
TEST_ALIGN_FILE3 = [('CodonAlign/nucl3.fa', 'CodonAlign/pro3.aln'), 'index']
TEST_ALIGN_FILE4 = [('CodonAlign/nucl4.fa', 'CodonAlign/pro4.aln'), 'index']
TEST_ALIGN_FILE5 = [('CodonAlign/nucl5.fa', 'CodonAlign/pro5.aln'), 'parse']
TEST_ALIGN_FILE6 = [('CodonAlign/egfr_nucl.fa', 'CodonAlign/egfr_pro.aln', 'CodonAlign/egfr_id'), 'id']
TEST_ALIGN_FILE7 = [('CodonAlign/drosophilla.fasta', 'CodonAlign/adh.aln'), 'index']

temp_dir = tempfile.mkdtemp()

class TestCodonSeq(unittest.TestCase):
    def test_seq(self):
        codonseq1 = CodonAlign.CodonSeq('AAATTT---TTTGGACCC', rf_table=[0,3,6,9,12])
        self.assertEqual(len(codonseq1), 18)
        self.assertEqual(codonseq1.get_codon_num(), 5)
        self.assertEqual(str(codonseq1.get_codon(0)), 'AAA')
        self.assertEqual(str(codonseq1.get_codon(-1)), 'CCC')
        self.assertEqual(str(codonseq1.get_codon(slice(1,3))), 'TTT---')
        self.assertEqual(str(codonseq1.get_codon(slice(None,None,-1))), 'CCCGGATTT---TTTAAA')

        self.assertRaises(ValueError, CodonAlign.CodonSeq, 'AAA-TT')
        self.assertRaises(AssertionError, CodonAlign.CodonSeq, 'AAA-T')
        self.assertRaises(ValueError, CodonAlign.CodonSeq, 'YVVRRDQQQ')
        self.assertTrue(isinstance(codonseq1.toSeq(), Seq))

class TestCodonAlignment(unittest.TestCase):
    def setUp(self):
        codonseq1 = CodonAlign.CodonSeq('AAATTT---TTTGGACCC', CodonAlign.default_codon_alphabet)
        codonseq2 = CodonAlign.CodonSeq('AAGTTT---TTTGGGCCC', CodonAlign.default_codon_alphabet)
        codonseq3 = CodonAlign.CodonSeq('AAGTAT---TTTGGACCC', CodonAlign.default_codon_alphabet)
        codonseq4 = CodonAlign.CodonSeq('AACTTT---TTTGGACGC', CodonAlign.default_codon_alphabet)

        self.seqrec = [SeqRecord(codonseq1, id="alpha"), 
                       SeqRecord(codonseq2, id="beta" ),
                       SeqRecord(codonseq3, id="gamma"),
                       SeqRecord(codonseq4, id="delta")]

    def test_align(self):
        codonAlign = CodonAlign.CodonAlignment(self.seqrec)
        self.assertEqual(codonAlign.get_aln_length(), 6)
        self.assertTrue(isinstance(codonAlign.toMultipleSeqAlignment(), MultipleSeqAlignment))

class TestBuildAndIO(unittest.TestCase):
    def setUp(self):
        self.aln_file = [TEST_ALIGN_FILE1,
                         TEST_ALIGN_FILE2,
                         TEST_ALIGN_FILE3,
                         TEST_ALIGN_FILE4,
                         TEST_ALIGN_FILE5,
                         TEST_ALIGN_FILE6]
        alns = []
        for i in self.aln_file:
            if i[1] == 'parse':
                nucl = SeqIO.parse(i[0][0], 'fasta', alphabet=IUPAC.IUPACUnambiguousDNA())
                prot = AlignIO.read(i[0][1], 'clustal', alphabet=IUPAC.protein)
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    caln = CodonAlign.build(prot, nucl, alphabet=CodonAlign.default_codon_alphabet)
            elif i[1] == 'index':
                nucl = SeqIO.index(i[0][0], 'fasta', alphabet=IUPAC.IUPACUnambiguousDNA())
                prot = AlignIO.read(i[0][1], 'clustal', alphabet=IUPAC.protein)
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    caln = CodonAlign.build(prot, nucl, alphabet=CodonAlign.default_codon_alphabet, max_score=20)
            elif i[1] == 'id':
                nucl = SeqIO.parse(i[0][0], 'fasta', alphabet=IUPAC.IUPACUnambiguousDNA())
                prot = AlignIO.read(i[0][1], 'clustal', alphabet=IUPAC.protein)
                with open(i[0][2]) as handle:
                    id = dict((i.split()[0], i.split()[1]) for i in handle)
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    caln = CodonAlign.build(prot, nucl, corr_dict=id, alphabet=CodonAlign.default_codon_alphabet)
            alns.append(caln)
            nucl.close() # Close the indexed FASTA file
        self.alns = alns

    def test_IO(self):
        self.assertEqual(len(self.alns), 6)
        #print temp_dir
        for n, i in enumerate(self.alns):
            aln = i.toMultipleSeqAlignment()
            AlignIO.write(aln, temp_dir + '/aln' + str(n) + '.clw', 'clustal')

class Test_build(unittest.TestCase):
    def setUp(self):
        # Test set 1
        seq1 = SeqRecord(Seq('TCAGGGACTGCGAGAACCAAGCTACTGCTGCTGCTGGCTGCGCTCTGCGCCGCAGGTGGGGCGCTGGAG', \
                alphabet=IUPAC.IUPACUnambiguousDNA()), id='pro1')
        seq2 = SeqRecord(Seq('TCAGGGACTTCGAGAACCAAGCGCTCCTGCTGCTGGCTGCGCTCGGCGCCGCAGGTGGAGCACTGGAG', \
                alphabet=IUPAC.IUPACUnambiguousDNA()), id='pro2')
        pro1 = SeqRecord(Seq('SGTARTKLLLLLAALCAAGGALE', alphabet=IUPAC.protein),id='pro1')
        pro2 = SeqRecord(Seq('SGTSRTKRLLLLAALGAAGGALE', alphabet=IUPAC.protein),id='pro2')
        aln1 = MultipleSeqAlignment([pro1, pro2])
        self.aln1 = aln1
        self.seqlist1 = [seq1, seq2]
        # Test set 2
        #                      M  K  K  H  E L(F)L  C  Q  G  T  S  N  K  L  T  Q(L)L  G  T  F  E  D  H  F  L  S  L  Q  R  M  F  N  N  C  E  V  V  
        seq3 = SeqRecord(Seq('ATGAAAAAGCACGAGTTACTTTGCCAAGGGACAAGTAACAAGCTCACCCAGTTGGGCACTTTTGAAGACCACTTTCTGAGCCTACAGAGGATGTTCAACAACTGTGAGGTGGTCCTTGGGAATTTGGAAATTACCTACATGCAGAGTAGTTACAACCTTTCTTTTCTCAAGACCATCCAGGAGGTTGCCGGCTATGTACTCATTGCCCTC', alphabet=IUPAC.IUPACUnambiguousDNA()), id='pro1')
        #seq4 =SeqRecord(Seq('ATGAAAAAGCACGAGTT CTTTGCCAAGGGACAAGTAACAAGCTCACCCAGTTGGGCACTTTTGAAGACCACTTTCTGAGCCTACAGAGGATGTTCAACAA TGTGAGGTGGTCCTTGGGAATTTGGAAATTACCTACATGCAGAGTAGTTACAACCTTTCTTTTCTCAAGACCATCCAGGAGGTTGCCGGCTATGTACTCATTGCCCTC', alphabet=IUPAC.IUPACUnambiguousDNA()), id='pro2')
        seq4 = SeqRecord(Seq('ATGAAAAAGCACGAGTTCTTTGCCAAGGGACAAGTAACAAGCTCACCCAGTTGGGCACTTTTGAAGACCACTTTCTGAGCCTACAGAGGATGTTCAACAATGTGAGGTGGTCCTTGGGAATTTGGAAATTACCTACATGCAGAGTAGTTACAACCTTTCTTTTCTCAAGACCATCCAGGAGGTTGCCGGCTATGTACTCATTGCCCTC', alphabet=IUPAC.IUPACUnambiguousDNA()), id='pro2')
        #seq5 =SeqRecord(Seq('ATGAAAAAGCACGAGTT CTTTGCCAAGGGACAAGTAACAAGCTCACCC  TTGGGCACTTTTGAAGACCACTTTCTGAGCCTACAGAGGATGTTCAACAACTGTGAGGTGGTCCTTGGGAATTTGGAAATTACCTACATGCAGAGTAGTTACAACCTTTCTTTTCTCAAGACCATCCAGGAGGTTGCCGGCTATGTACTCATTGCCCTC', alphabet=IUPAC.IUPACUnambiguousDNA()), id='pro3')
        seq5 = SeqRecord(Seq('ATGAAAAAGCACGAGTTACTTTGCCAAGGGACAAGTAACAAGCTCACCCTTGGGCACTTTTGAAGACCACTTTCTGAGCCTACAGAGGATGTTCAACAACTGTGAGGTGGTCCTTGGGAATTTGGAAATTACCTACATGCAGAGTAGTTACAACCTTTCTTTTCTCAAGACCATCCAGGAGGTTGCCGGCTATGTACTCATTGCCCTC', alphabet=IUPAC.IUPACUnambiguousDNA()), id='pro3')
        pro3 = SeqRecord(Seq('MKKHELLCQGTSNKLTQLGTFEDHFLSLQRMFNNCEVVLGNLEITYMQSSYNLSFLKTIQEVAGYVLIAL', alphabet=IUPAC.protein), id='pro1')
        pro4 = SeqRecord(Seq('MKKHEFLCQGTSNKLTQLGTFEDHFLSLQRMFNNCEVVLGNLEITYMQSSYNLSFLKTIQEVAGYVLIAL', alphabet=IUPAC.protein), id='pro2')
        pro5 = SeqRecord(Seq('MKKHELLCQGTSNKLTLLGTFEDHFLSLQRMFNNCEVVLGNLEITYMQSSYNLSFLKTIQEVAGYVLIAL', alphabet=IUPAC.protein), id='pro3')
        aln2 = MultipleSeqAlignment([pro3, pro4, pro5])
        self.aln2 = aln2
        self.seqlist2 = [seq3, seq4, seq5]
    
    def test_build(self):
        codon_aln1 = CodonAlign.build(self.aln1, self.seqlist1)
        codon_aln2 = CodonAlign.build(self.aln2, self.seqlist2)

class Test_dn_ds(unittest.TestCase):
    def setUp(self):
        nucl = SeqIO.parse(TEST_ALIGN_FILE6[0][0], 'fasta', alphabet=IUPAC.IUPACUnambiguousDNA())
        prot = AlignIO.read(TEST_ALIGN_FILE6[0][1], 'clustal', alphabet=IUPAC.protein)
        with open(TEST_ALIGN_FILE6[0][2]) as handle:
            id_corr = dict((i.split()[0], i.split()[1]) for i in handle)
        aln = CodonAlign.build(prot, nucl, corr_dict=id_corr, alphabet=CodonAlign.default_codon_alphabet)
        self.aln = aln

    def test_dn_ds(self):
        from Bio.CodonAlign.CodonSeq import cal_dn_ds
        codon_seq1 = self.aln[0]
        codon_seq2 = self.aln[1]
        dN, dS = cal_dn_ds(codon_seq1, codon_seq2, method='NG86')
        self.assertAlmostEqual(round(dN, 4), 0.0209, places=4)
        self.assertAlmostEqual(round(dS, 4), 0.0178, places=4)
        dN, dS = cal_dn_ds(codon_seq1, codon_seq2, method='LWL85')
        self.assertAlmostEqual(round(dN, 4), 0.0203, places=4)
        self.assertAlmostEqual(round(dS, 4), 0.0164, places=4)

        try:
            import scipy
        except ImportError:
            # Silently skip the rest of the test
            return

        # This should be present:
        from scipy.linalg import expm
        dN, dS = cal_dn_ds(codon_seq1, codon_seq2, method='YN00')
        self.assertAlmostEqual(round(dN, 4), 0.0198, places=4)
        self.assertAlmostEqual(round(dS, 4), 0.0222, places=4)

        try:
            # New in scipy v0.11
            from scipy.optimize import minimize
            dN, dS = cal_dn_ds(codon_seq1, codon_seq2, method='ML')
            self.assertAlmostEqual(round(dN, 4), 0.0194, places=4)
            self.assertAlmostEqual(round(dS, 4), 0.0217, places=4)
        except ImportError:
            # TODO - Show a warning?
            pass


from run_tests import is_numpy
try:
    from math import lgamma # New in Python 2.7
except ImportError:
    lgamma = None
if is_numpy() and lgamma:
    class Test_MK(unittest.TestCase):
        def test_mk(self):
            p = SeqIO.index(TEST_ALIGN_FILE7[0][0], 'fasta', alphabet=IUPAC.IUPACUnambiguousDNA())
            pro_aln = AlignIO.read(TEST_ALIGN_FILE7[0][1], 'clustal', alphabet=IUPAC.protein)
            codon_aln = CodonAlign.build(pro_aln, p)
            p.close() # Close indexed FASTA file
            self.assertAlmostEqual(round(CodonAlign.mktest([codon_aln[1:12], codon_aln[12:16], codon_aln[16:]]), 4), 0.0021, places=4)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
