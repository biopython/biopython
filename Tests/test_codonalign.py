# Copyright (C) 2013 by Zheng Ruan (zruan1991@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for the Bio.codonalign modules."""

import warnings
import tempfile
import unittest

from Bio import BiopythonWarning, BiopythonExperimentalWarning
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Data import CodonTable


with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonExperimentalWarning)
    from Bio import codonalign


TEST_ALIGN_FILE1 = [("codonalign/nucl1.fa", "codonalign/pro1.aln"), "parse"]
TEST_ALIGN_FILE2 = [("codonalign/nucl2.fa", "codonalign/pro2.aln"), "parse"]
TEST_ALIGN_FILE3 = [("codonalign/nucl3.fa", "codonalign/pro3.aln"), "index"]
TEST_ALIGN_FILE4 = [("codonalign/nucl4.fa", "codonalign/pro4.aln"), "index"]
TEST_ALIGN_FILE5 = [("codonalign/nucl5.fa", "codonalign/pro5.aln"), "parse"]
TEST_ALIGN_FILE6 = [
    ("codonalign/egfr_nucl.fa", "codonalign/egfr_pro.aln", "codonalign/egfr_id"),
    "id",
]
TEST_ALIGN_FILE7 = [("codonalign/drosophilla.fasta", "codonalign/adh.aln"), "index"]

temp_dir = tempfile.mkdtemp()


class TestCodonSeq(unittest.TestCase):
    def test_seq(self):
        codonseq1 = codonalign.CodonSeq("AAATTT---TTTGGACCC", rf_table=[0, 3, 6, 9, 12])
        self.assertEqual(len(codonseq1), 18)
        self.assertEqual(codonseq1.get_codon_num(), 5)
        self.assertEqual(str(codonseq1.get_codon(0)), "AAA")
        self.assertEqual(str(codonseq1.get_codon(-1)), "CCC")
        self.assertEqual(str(codonseq1.get_codon(slice(1, 3))), "TTT---")
        self.assertEqual(
            str(codonseq1.get_codon(slice(None, None, -1))), "CCCGGATTT---TTTAAA"
        )

        self.assertRaises(ValueError, codonalign.CodonSeq, "AAA-T")
        self.assertIsInstance(codonseq1.toSeq(), Seq)


class TestCodonAlignment(unittest.TestCase):
    def setUp(self):
        codonseq1 = codonalign.CodonSeq("AAATTT---TTTGGACCC")
        codonseq2 = codonalign.CodonSeq("AAGTTT---TTTGGGCCC")
        codonseq3 = codonalign.CodonSeq("AAGTAT---TTTGGACCC")
        codonseq4 = codonalign.CodonSeq("AACTTT---TTTGGACGC")

        self.seqrec = [
            SeqRecord(codonseq1, id="alpha"),
            SeqRecord(codonseq2, id="beta"),
            SeqRecord(codonseq3, id="gamma"),
            SeqRecord(codonseq4, id="delta"),
        ]

    def test_align(self):
        codonAlign = codonalign.CodonAlignment(self.seqrec)
        self.assertEqual(codonAlign.get_aln_length(), 6)
        self.assertIsInstance(codonAlign.toMultipleSeqAlignment(), MultipleSeqAlignment)


class TestAddition(unittest.TestCase):
    def setUp(self):
        self.seq1 = SeqRecord(Seq("ATGTCTCGT"), id="pro1")
        self.seq2 = SeqRecord(Seq("ATGCGT"), id="pro2")
        self.pro1 = SeqRecord(Seq("MSR"), id="pro1")
        self.pro2 = SeqRecord(Seq("M-R"), id="pro2")
        self.aln = MultipleSeqAlignment([self.pro1, self.pro2])
        self.codon_aln = codonalign.build(self.aln, [self.seq1, self.seq2])

        tail1 = SeqRecord(Seq("AAA"), id="pro1")
        tail2 = SeqRecord(Seq("AAA"), id="pro2")
        self.multi_aln = MultipleSeqAlignment([tail1, tail2])

    def test_addition_MultipleSeqAlignment(self):
        """Check addition of CodonAlignment and MultipleSeqAlignment."""
        new_aln1 = self.codon_aln + self.multi_aln

        self.assertIsInstance(new_aln1, MultipleSeqAlignment)
        for x in range(len(self.codon_aln)):
            self.assertEqual(
                new_aln1[x].seq, self.codon_aln[x].seq + self.multi_aln[x].seq
            )

        new_aln2 = self.multi_aln + self.codon_aln

        self.assertIsInstance(new_aln2, MultipleSeqAlignment)
        for x in range(len(self.codon_aln)):
            self.assertEqual(
                new_aln2[x].seq, self.multi_aln[x].seq + self.codon_aln[x].seq
            )

    def test_addition_CodonAlignment(self):
        """Check addition of CodonAlignment and CodonAlignment."""
        new_aln = self.codon_aln + self.codon_aln

        self.assertIsInstance(new_aln, codonalign.CodonAlignment)
        for x in range(len(self.codon_aln)):
            self.assertEqual(
                new_aln[x].seq, self.codon_aln[x].seq + self.codon_aln[x].seq
            )

    def test_ValueError(self):
        """Check that ValueError is thrown for Alignments of different lengths."""
        # original len(self.aln) = 2 , len(aln) = 3
        aln = MultipleSeqAlignment(
            [self.pro1, self.pro2, SeqRecord(Seq("M--"), id="pro3")]
        )
        triple_codon = codonalign.build(
            aln, [self.seq1, self.seq2, SeqRecord(Seq("ATG"), id="pro3")]
        )
        with self.assertRaises(ValueError):
            triple_codon + self.multi_aln
        with self.assertRaises(ValueError):
            triple_codon + self.codon_aln

    def test_TypeError(self):
        """Check that TypeError is thrown for non CodonAlignment/MultipleSequenceAlignment objects."""
        for obj in [0, "string", ["str1", "str2"], Seq("ATGTCTCGT")]:
            with self.assertRaises(TypeError):
                self.codon_aln + obj


class TestBuildAndIO(unittest.TestCase):
    def setUp(self):
        self.aln_file = [
            TEST_ALIGN_FILE1,
            TEST_ALIGN_FILE2,
            TEST_ALIGN_FILE3,
            TEST_ALIGN_FILE4,
            TEST_ALIGN_FILE5,
            TEST_ALIGN_FILE6,
        ]
        alns = []
        for i in self.aln_file:
            if i[1] == "parse":
                nucl = SeqIO.parse(i[0][0], "fasta")
                prot = AlignIO.read(i[0][1], "clustal")
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    caln = codonalign.build(prot, nucl)
            elif i[1] == "index":
                nucl = SeqIO.index(i[0][0], "fasta")
                prot = AlignIO.read(i[0][1], "clustal")
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    caln = codonalign.build(prot, nucl, max_score=20)
                nucl.close()  # Close the indexed FASTA file
            elif i[1] == "id":
                nucl = SeqIO.parse(i[0][0], "fasta")
                prot = AlignIO.read(i[0][1], "clustal")
                with open(i[0][2]) as handle:
                    id = {i.split()[0]: i.split()[1] for i in handle}
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    caln = codonalign.build(prot, nucl, corr_dict=id)
            alns.append(caln)
        self.alns = alns

    def test_IO(self):
        self.assertEqual(len(self.alns), 6)
        # print(temp_dir)
        for n, i in enumerate(self.alns):
            aln = i.toMultipleSeqAlignment()
            AlignIO.write(aln, temp_dir + "/aln" + str(n) + ".clw", "clustal")


class Test_build(unittest.TestCase):
    def setUp(self):
        # Test set 1
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
        aln1 = MultipleSeqAlignment([pro1, pro2])
        self.aln1 = aln1
        self.seqlist1 = [seq1, seq2]
        # Test set 2
        #                      M  K  K  H  E L(F)L  C  Q  G  T  S  N  K  L  T  Q(L)L  G  T  F  E  D  H  F  L  S  L  Q  R  M  F  N  N  C  E  V  V
        seq3 = SeqRecord(
            Seq(
                "ATGAAAAAGCACGAGTTACTTTGCCAAGGGACAAGTAACAAGCTCACCCAGTTGGGCACTTTTGAAGACCACTTTCTGAGCCTACAGAGGATGTTCAACAACTGTGAGGTGGTCCTTGGGAATTTGGAAATTACCTACATGCAGAGTAGTTACAACCTTTCTTTTCTCAAGACCATCCAGGAGGTTGCCGGCTATGTACTCATTGCCCTC"
            ),
            id="pro1",
        )
        # seq4 =SeqRecord(Seq('ATGAAAAAGCACGAGTT CTTTGCCAAGGGACAAGTAACAAGCTCACCCAGTTGGGCACTTTTGAAGACCACTTTCTGAGCCTACAGAGGATGTTCAACAA TGTGAGGTGGTCCTTGGGAATTTGGAAATTACCTACATGCAGAGTAGTTACAACCTTTCTTTTCTCAAGACCATCCAGGAGGTTGCCGGCTATGTACTCATTGCCCTC'), id='pro2')
        seq4 = SeqRecord(
            Seq(
                "ATGAAAAAGCACGAGTTCTTTGCCAAGGGACAAGTAACAAGCTCACCCAGTTGGGCACTTTTGAAGACCACTTTCTGAGCCTACAGAGGATGTTCAACAATGTGAGGTGGTCCTTGGGAATTTGGAAATTACCTACATGCAGAGTAGTTACAACCTTTCTTTTCTCAAGACCATCCAGGAGGTTGCCGGCTATGTACTCATTGCCCTC"
            ),
            id="pro2",
        )
        # seq5 =SeqRecord(Seq('ATGAAAAAGCACGAGTT CTTTGCCAAGGGACAAGTAACAAGCTCACCC  TTGGGCACTTTTGAAGACCACTTTCTGAGCCTACAGAGGATGTTCAACAACTGTGAGGTGGTCCTTGGGAATTTGGAAATTACCTACATGCAGAGTAGTTACAACCTTTCTTTTCTCAAGACCATCCAGGAGGTTGCCGGCTATGTACTCATTGCCCTC'), id='pro3')
        seq5 = SeqRecord(
            Seq(
                "ATGAAAAAGCACGAGTTACTTTGCCAAGGGACAAGTAACAAGCTCACCCTTGGGCACTTTTGAAGACCACTTTCTGAGCCTACAGAGGATGTTCAACAACTGTGAGGTGGTCCTTGGGAATTTGGAAATTACCTACATGCAGAGTAGTTACAACCTTTCTTTTCTCAAGACCATCCAGGAGGTTGCCGGCTATGTACTCATTGCCCTC"
            ),
            id="pro3",
        )
        pro3 = SeqRecord(
            Seq(
                "MKKHELLCQGTSNKLTQLGTFEDHFLSLQRMFNNCEVVLGNLEITYMQSSYNLSFLKTIQEVAGYVLIAL"
            ),
            id="pro1",
        )
        pro4 = SeqRecord(
            Seq(
                "MKKHEFLCQGTSNKLTQLGTFEDHFLSLQRMFNNCEVVLGNLEITYMQSSYNLSFLKTIQEVAGYVLIAL"
            ),
            id="pro2",
        )
        pro5 = SeqRecord(
            Seq(
                "MKKHELLCQGTSNKLTLLGTFEDHFLSLQRMFNNCEVVLGNLEITYMQSSYNLSFLKTIQEVAGYVLIAL"
            ),
            id="pro3",
        )
        aln2 = MultipleSeqAlignment([pro3, pro4, pro5])
        self.aln2 = aln2
        self.seqlist2 = [seq3, seq4, seq5]

        # Test set 3
        # use Yeast mitochondrial codon table
        seq6 = SeqRecord(
            Seq(
                "ATGGCAAGGGACCACCCAGTTGGGCACTGATATGATCGGGTGTATTTGCAGAGTAGTAACCTTTCTTTTCTCAAGACCATCCAG"
            ),
            id="pro6",
        )
        seq7 = SeqRecord(
            Seq(
                "ATGGCAAGGCACCATCCAGTTGAGCACTGATATGATCGGGTGTATTTGCAGAGTAGTAACGTGTCTCTGCTCAAGACCATCCAG"
            ),
            id="pro7",
        )
        seq8 = SeqRecord(
            Seq(
                "ATGGCAGGGGACCACCCAGTTGGGCACTGATATGATCGTGTGTATCTGCAGAGTAGTAACCACTCTTTTCTCATGACCATCCAG"
            ),
            id="pro8",
        )
        pro6 = SeqRecord(Seq("MARDHPVGHWYDRVYLQSSNTSFTKTIQ"), id="pro6")
        pro7 = SeqRecord(Seq("MARHHPVEHWYDRVYLQSSNVSTTKTIQ"), id="pro7")
        pro8 = SeqRecord(Seq("MAGDHPVGHWYDRVYTQSSNHSFTMTIQ"), id="pro8")
        aln3 = MultipleSeqAlignment([pro6, pro7, pro8])
        self.aln3 = aln3
        self.seqlist3 = [seq6, seq7, seq8]
        self.codontable3 = CodonTable.unambiguous_dna_by_id[3]

    def test_build(self):
        codon_aln1 = codonalign.build(self.aln1, self.seqlist1)
        codon_aln2 = codonalign.build(self.aln2, self.seqlist2)
        codon_aln3 = codonalign.build(
            self.aln3, self.seqlist3, codon_table=self.codontable3
        )
        codon_aln4 = codonalign.build(self.aln1, self.seqlist1, complete_protein=True)


class Test_dn_ds(unittest.TestCase):
    def setUp(self):
        nucl = SeqIO.parse(TEST_ALIGN_FILE6[0][0], "fasta")
        prot = AlignIO.read(TEST_ALIGN_FILE6[0][1], "clustal")
        with open(TEST_ALIGN_FILE6[0][2]) as handle:
            id_corr = {i.split()[0]: i.split()[1] for i in handle}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonWarning)
            aln = codonalign.build(prot, nucl, corr_dict=id_corr)
        self.aln = aln

    def test_dn_ds(self):
        from Bio.codonalign.codonseq import cal_dn_ds

        codon_seq1 = self.aln[0]
        codon_seq2 = self.aln[1]
        dN, dS = cal_dn_ds(codon_seq1, codon_seq2, method="NG86")
        self.assertAlmostEqual(dN, 0.0209, places=4)
        self.assertAlmostEqual(dS, 0.0178, places=4)
        dN, dS = cal_dn_ds(codon_seq1, codon_seq2, method="LWL85")
        self.assertAlmostEqual(dN, 0.0203, places=4)
        self.assertAlmostEqual(dS, 0.0164, places=4)

        try:
            import scipy
        except ImportError:
            # Silently skip the rest of the test
            return

        # This should be present:
        from scipy.linalg import expm

        dN, dS = cal_dn_ds(codon_seq1, codon_seq2, method="YN00")
        self.assertAlmostEqual(dN, 0.0198, places=4)
        self.assertAlmostEqual(dS, 0.0222, places=4)

        try:
            # New in scipy v0.11
            from scipy.optimize import minimize

            dN, dS = cal_dn_ds(codon_seq1, codon_seq2, method="ML")
            self.assertAlmostEqual(dN, 0.0194, places=4)
            self.assertAlmostEqual(dS, 0.0217, places=4)
        except ImportError:
            # TODO - Show a warning?
            pass

    def test_dn_ds_matrix(self):
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
        dn, ds = self.aln.get_dn_ds_matrix()
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
        dn, ds = self.aln.get_dn_ds_matrix(
            method="LWL85", codon_table=CodonTable.unambiguous_dna_by_id[3]
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
            p = SeqIO.index(TEST_ALIGN_FILE7[0][0], "fasta")
            pro_aln = AlignIO.read(TEST_ALIGN_FILE7[0][1], "clustal")
            codon_aln = codonalign.build(pro_aln, p)
            p.close()  # Close indexed FASTA file
            self.assertAlmostEqual(
                codonalign.mktest([codon_aln[1:12], codon_aln[12:16], codon_aln[16:]]),
                0.0021,
                places=4,
            )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
