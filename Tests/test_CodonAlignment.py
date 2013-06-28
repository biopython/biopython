from Bio import CodonAlign, SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Align import MultipleSeqAlignment
import tempfile
import unittest

TEST_ALIGN_FILE1 = [('CodonAlign/nucl1.fa', 'CodonAlign/pro1.aln'), 'parse']
TEST_ALIGN_FILE2 = [('CodonAlign/nucl2.fa', 'CodonAlign/pro2.aln'), 'parse']
TEST_ALIGN_FILE3 = [('CodonAlign/nucl3.fa', 'CodonAlign/pro3.aln'), 'index']
TEST_ALIGN_FILE4 = [('CodonAlign/nucl4.fa', 'CodonAlign/pro4.aln'), 'index']
TEST_ALIGN_FILE5 = [('CodonAlign/nucl5.fa', 'CodonAlign/pro5.aln'), 'parse']
TEST_ALIGN_FILE6 = [('CodonAlign/egfr_nucl.fa', 'CodonAlign/egfr_pro.aln', 'CodonAlign/egfr_id'), 'id']

temp_dir = tempfile.mkdtemp()

class TestCodonSeq(unittest.TestCase):

    def test_seq(self):
        codonseq1 = CodonAlign.CodonSeq('AAATTT---TTTGGACCC', CodonAlign.default_codon_alphabet)
        self.assertEqual(len(codonseq1), 18)
        self.assertEqual(codonseq1.get_codon_num(), 5)
        self.assertEqual(str(codonseq1[0]), 'AAA')
        self.assertEqual(str(codonseq1[-1]), 'CCC')
        self.assertEqual(str(codonseq1[1:3]), 'TTT---')
        self.assertEqual(str(codonseq1[::-1]), 'CCCGGATTT---TTTAAA')

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
                caln = CodonAlign.build(prot, nucl, alphabet=CodonAlign.default_codon_alphabet)
            elif i[1] == 'index':
                nucl = SeqIO.index(i[0][0], 'fasta', alphabet=IUPAC.IUPACUnambiguousDNA())
                prot = AlignIO.read(i[0][1], 'clustal', alphabet=IUPAC.protein)
                caln = CodonAlign.build(prot, nucl, alphabet=CodonAlign.default_codon_alphabet)
            elif i[1] == 'id':
                nucl = SeqIO.parse(i[0][0], 'fasta', alphabet=IUPAC.IUPACUnambiguousDNA())
                prot = AlignIO.read(i[0][1], 'clustal', alphabet=IUPAC.protein)
                id = {i.split()[0]: i.split()[1] for i in open(i[0][2]).readlines()}
                caln = CodonAlign.build(prot, nucl, corr_dict=id, alphabet=CodonAlign.default_codon_alphabet)

            alns.append(caln)
        self.alns = alns

    def test_IO(self):
        self.assertEqual(len(self.alns), 6)
        #print temp_dir
        for n, i in enumerate(self.alns):
            aln = i.toMultipleSeqAlignment()
            AlignIO.write(aln, temp_dir + '/aln' + str(n) + '.clw', 'clustal')


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
        
