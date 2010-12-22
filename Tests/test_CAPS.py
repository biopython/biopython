import unittest

from Bio import CAPS
from Bio.Restriction import EcoRI, AluI
from Bio import Alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

def createAlignment(sequences, alphabet):
    """Create an Alignment object from a list of sequences"""
    return MultipleSeqAlignment((SeqRecord(Seq(s,alphabet), id="sequence%i"%(i+1)) \
                                 for (i,s) in enumerate(sequences)),
                                alphabet)
    
class TestCAPS(unittest.TestCase):

    def test_trivial(self):
        enzymes = [EcoRI]
        alignment = ["gaattc",
                     "gaactc",
                    ]
        align = createAlignment(alignment, Alphabet.generic_dna)
        map = CAPS.CAPSMap(align, enzymes)

        self.assertEqual(len(map.dcuts), 1)
        self.assertEqual(map.dcuts[0].enzyme, EcoRI)
        self.assertEqual(map.dcuts[0].start, 1)
        self.assertEqual(map.dcuts[0].cuts_in, [0])
        self.assertEqual(map.dcuts[0].blocked_in, [1])


    def test(self):
        alignment = [
"""\
AAAagaattcTAGATATACCAAACCAGAGAAAACAAATACATAATCGGAGAAATACAGAT
AGAGAGCGAGAGAGATCGACGGCGAAGCTCTTTACCCGGAAACCATTGAAATCGGACGGT
TTAGTGAAAATGGAGGATCAAGTagAtTTTGGGTTCCGTCCGAACGACGAGGAGCTCGTT
GGTCACTATCTCCGTAACAAAATCGAAGGAAACACTAGCCGCGACGTTGAAGTAGCCATC
AGCGAGGTCAACATCTGTAGCTACGATCCTTGGAACTTGCGCTGTAAGTTCCGAATTTTC
""",
"""\
AAAagaTttcTAGATATACCAAACCAGAGAAAACAAATACATAATCGGAGAAATACAGAT
AGAGAGCGAGAGAGATCGACGGCGAAGCTCTTTACCCGGAAACCATTGAAATCGGACGGT
TTAGTGAAAATGGAGGATCAAGTagctTTTGGGTTCCGTCCGAACGACGAGGAGCTCGTT
GGTCACTATCTCCGTAACAAAATCGAAGGAAACACTAGCCGCGACGTTGAAGTAGCCATC
AGCGAGGTCAACATCTGTAGCTACGATCCTTGGAACTTGCGCTGTAAGTTCCGAATTTTC
""",
"""\
AAAagaTttcTAGATATACCAAACCAGAGAAAACAAATACATAATCGGAGAAATACAGAT
AGAGAGCGAGAGAGATCGACGGCGAAGCTCTTTACCCGGAAACCATTGAAATCGGACGGT
TTAGTGAAAATGGAGGATCAAGTagctTTTGGGTTCCGTCCGAACGACGAGGAGCTCGTT
GGTCACTATCTCCGTAACAAAATCGAAGGAAACACTAGCCGCGACGTTGAAGTAGCCATC
AGCGAGGTCAACATCTGTAGCTACGATCCTTGGAACTTGCGCTGTAAGTTCCGAATTTTC
""",
                    ]
        enzymes = [EcoRI, AluI]
        align = createAlignment(alignment, Alphabet.generic_dna)
        map = CAPS.CAPSMap(align, enzymes)

        self.assertEqual(len(map.dcuts), 2)
        self.assertEqual(map.dcuts[0].enzyme, EcoRI)
        self.assertEqual(map.dcuts[0].start, 5)
        self.assertEqual(map.dcuts[0].cuts_in, [0])
        self.assertEqual(map.dcuts[0].blocked_in, [1,2])
        self.assertEqual(map.dcuts[1].enzyme, AluI)
        self.assertEqual(map.dcuts[1].start, 144)
        self.assertEqual(map.dcuts[1].cuts_in, [1,2])
        self.assertEqual(map.dcuts[1].blocked_in, [0])


    def testNoCAPS(self):
        alignment = ["aaaaaaaaaaaaaaaaaaaa",
                     "aaaaaaaaaaaaaaaaaaaa",
                    ]
        enzymes = []
        align = createAlignment(alignment, Alphabet.generic_nucleotide)
        map = CAPS.CAPSMap(align, enzymes)
        self.assertEqual(map.dcuts, [])


    def test_uneven(self):
        alignment = ["aaaaaaaaaaaaaa",
                     "aaaaaaaaaaaaaa", #we'll change this below
                     "aaaaaaaaaaaaaa",
                    ]
        align = createAlignment(alignment, Alphabet.generic_nucleotide)
        align[1].seq = align[1].seq[:8] #evil
        self.assertRaises(CAPS.AlignmentHasDifferentLengthsError,
                          CAPS.CAPSMap,
                          align)

 

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
