import unittest
from Bio import TAIR

test_agis = [
    "AT4G36450.1",
    "AT4G36900.1",
    "AT5G63980.1"
    ]


class TAIRDirect(unittest.TestCase):
    """
    """
    test_datasets = [
        "protein",
        "upstream_500",
        "downstream_500",
        "intergenic",
        "intron",
        "3prime_utr",
        "5prime_utr"
        ]

    def test_transcript(self):
        test_seqs = {
            "AT4G36450.1": ("ATGGCGATGC", 1086),
            "AT4G36900.1": ("AGTGTCGGTG", 1024),
            "AT5G63980.1": ("GACATATATT", 1383)
            ]
        for agi, tests in test_agis:
            seq = TAIR.get(agi, dataset, "rep_gene")
            self.assertEqual(seq.id, agi)
            self.assertEqual(seq.seq[:10], test_values[0])
            self.assertEqual(len(seq.seq), test_values[1])

    def test_cds(self):
        test_seqs = {
            "AT4G36450.1": ("ATGGCGATGC", 1086),
            "AT4G36900.1": ("ATGGAGACGG", 591),
            "AT5G63980.1": ("ATGATGTCTA", 1224)
            ]
        for agi, tests in test_agis:
            seq = TAIR.get(agi, dataset, "rep_gene")
            self.assertEqual(seq.id, agi)
            self.assertEqual(seq.seq[:10], test_values[0])
            self.assertEqual(len(seq.seq), test_values[1])

    def test_gene(self):
        test_seqs = {
            "AT4G36450.1": ("ATGGCGATGC", 1169),
            "AT4G36900.1": ("AGTGTCGGTG", 1024),
            "AT5G63980.1": ("GACATATATT", 2122)
            ]
        for agi, tests in test_agis:
            seq = TAIR.get(agi, dataset, "rep_gene")
            self.assertEqual(seq.id, agi)
            self.assertEqual(seq.seq[:10], test_values[0])
            self.assertEqual(len(seq.seq), test_values[1])

    def test_intergenic(self):
        test_seqs = {
            "AT4G36450.1": ("MAMLVDPPNG", 361),
            "AT4G36900.1": ("METATEVATV", 196),
            "AT5G63980.1": ("MMSINCFRTA", 407)
            ]
        for agi, tests in test_agis:
            seq = TAIR.get(agi, dataset, "rep_gene")
            self.assertEqual(seq.id, agi)
            self.assertEqual(seq.seq[:10], test_values[0])
            self.assertEqual(len(seq.seq), test_values[1])

    def test_intergenic(self):
        test_seqs = {
            "AT4G36450.1": ("MAMLVDPPNG", 361),
            "AT4G36900.1": ("METATEVATV", 196),
            "AT5G63980.1": ("MMSINCFRTA", 407)
            ]
        for agi, tests in test_agis:
            seqs = TAIR.get(agi, dataset, "rep_gene")
            upstream_seq = seqs[0]
            downstream_seq = seqs[1]
            self.assertEqual(upstream_seq.id, agi)
            self.assertEqual(seq.seq[:10], test_values[0])
            self.assertEqual(len(seq.seq), test_values[1])
            self.assertEqual(seq.id, agi)
            self.assertEqual(seq.seq[:10], test_values[0])
            self.assertEqual(len(seq.seq), test_values[1])
    
    def test_intron(self):
        test_seqs = {
            "AT4G36450.1": ("MAMLVDPPNG", 361),
            "AT4G36900.1": ("METATEVATV", 196),
            "AT5G63980.1": ("MMSINCFRTA", 407)
            ]
        for agi, tests in test_agis:
            seq = TAIR.get(agi, dataset, "rep_gene")
            self.assertEqual(seq.id, agi)
            self.assertEqual(seq.seq[:10], test_values[0])
            self.assertEqual(len(seq.seq), test_values[1])

    def test_3prime_utr(self):
        test_seqs = {
            "AT4G36450.1": ("MAMLVDPPNG", 361),
            "AT4G36900.1": ("METATEVATV", 196),
            "AT5G63980.1": ("MMSINCFRTA", 407)
            ]
        for agi, tests in test_agis:
            seq = TAIR.get(agi, dataset, "rep_gene")
            self.assertEqual(seq.id, agi)
            self.assertEqual(seq.seq[:10], test_values[0])
            self.assertEqual(len(seq.seq), test_values[1])

    def test_5prime_utr(self):
        test_seqs = {
            "AT4G36450.1": ("MAMLVDPPNG", 361),
            "AT4G36900.1": ("METATEVATV", 196),
            "AT5G63980.1": ("MMSINCFRTA", 407)
            ]
        for agi, tests in test_agis:
            seq = TAIR.get(agi, dataset, "rep_gene")
            self.assertEqual(seq.id, agi)
            self.assertEqual(seq.seq[:10], test_values[0])
            self.assertEqual(len(seq.seq), test_values[1])


class TAIRNCBI(unittest.TestCase):
    pass


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)

