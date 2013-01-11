from __future__ import with_statement
import unittest
from Bio.Webtools import TAIR

# Make sure test is skipped if not internet connected
import requires_internet
requires_internet.check()


class TAIRDirect(unittest.TestCase):
    test_agis = [
        "AT4G36450.1",
        "AT4G36900.1",
        "AT5G63980.1"
        ]

    def _perform_test(self, returned, expected):
        for seq, expected_seq in zip(returned, expected):
            self.assertEqual(seq.id, expected_seq[0])
            self.assertEqual(str(seq.seq[:10]), expected_seq[1])
            self.assertEqual(len(seq.seq), expected_seq[2])

    def test_transcript(self):
        expected = [
            ("AT4G36450.1", "ATGGCGATGC", 1086),
            ("AT4G36900.1", "AGTGTCGGTG", 1024),
            ("AT5G63980.1", "GACATATATT", 1383)
            ]
        returned = TAIR.get(self.test_agis, "transcript", "rep_gene")
        self._perform_test(returned, expected)

    def test_cds(self):
        expected = [
            ("AT4G36450.1", "ATGGCGATGC", 1086),
            ("AT4G36900.1", "ATGGAGACGG", 591),
            ("AT5G63980.1", "ATGATGTCTA", 1224)
            ]
        returned = TAIR.get(self.test_agis, "cds", "rep_gene")
        self._perform_test(returned, expected)

    def test_gene(self):
        expected = [
            ("AT4G36450.1", "ATGGCGATGC", 1169),
            ("AT4G36900.1", "AGTGTCGGTG", 1024),
            ("AT5G63980.1", "GACATATATT", 2122)
            ]
        returned = TAIR.get(self.test_agis, "gene", "rep_gene")
        self._perform_test(returned, expected)

    def test_protein(self):
        expected = [
            ("AT4G36450.1", "MAMLVDPPNG", 361),
            ("AT4G36900.1", "METATEVATV", 196),
            ("AT5G63980.1", "MMSINCFRTA", 407)
            ]
        returned = TAIR.get(self.test_agis, "protein", "rep_gene")
        self._perform_test(returned, expected)

    def test_intergenic(self):
        expected = [
            ("AT4G36440-AT4G36450", "AGAGTCAAAA", 237),
            ("AT4G36450-AT4G36460", "TCCTCTTCTT", 1248),
            ("AT4G36890-AT4G36900", "ACCAATCTCT", 7033),
            ("AT4G36900-AT4G36910", "AATCCCCCTC", 788),
            ("AT5G63970-AT5G63980", "TTGCGCCGGC", 564),
            ("AT5G63980-AT5G63990", "GTATAAAGGG", 1324)
            ]
        returned = TAIR.get(self.test_agis, "intergenic", "rep_gene")
        self._perform_test(returned, expected)

    def test_intron(self):
        expected = [
            ("AT4G36450.1-1", "GTACGTTGGT", 83),
            ("AT5G63980.1-1", "GTTAGGGTTT", 226),
            ("AT5G63980.1-2", "GTTAGTTTGT", 106),
            ("AT5G63980.1-3", "GTGAAACTGC", 97),
            ("AT5G63980.1-4", "GTACGTTTTA", 119),
            ("AT5G63980.1-5", "GTAAATTGCT", 91),
            ("AT5G63980.1-6", "GTAACATTAA", 100)
            ]
        returned = TAIR.get(self.test_agis, "intron", "rep_gene")
        self._perform_test(returned, expected)

    def test_3prime_utr(self):
        expected = [
            ("AT4G36900.1", "GAAAGCAAAA", 232),
            ("AT5G63980.1", "TTTGTTTTTT", 131)
            ]
        returned = TAIR.get(self.test_agis, "3prime_utr", "rep_gene")
        self._perform_test(returned, expected)

    def test_5prime_utr(self):
        expected = [
            ("AT4G36900.1", "AGTGTCGGTG", 201),
            ("AT5G63980.1", "GACATATATT", 28)
            ]
        returned = TAIR.get(self.test_agis, "5prime_utr", "rep_gene")
        self._perform_test(returned, expected)

    def test_upstream500(self):
        expected = [
            ("AT4G36450", "TCTCAAAATA", 500),
            ("AT4G36900", "CTTAAATTAA", 500),
            ("AT5G63980", "TTTTTAGTCT", 500)
            ]
        returned = TAIR.get(self.test_agis, "upstream_500", "rep_gene")
        self._perform_test(returned, expected)

    def test_upstream1000(self):
        expected = [
            ("AT4G36450", "TTTAAATTTA", 1000),
            ("AT4G36900", "AAATATCCTC", 1000),
            ("AT5G63980", "CTAAAATCAA", 1000)
            ]
        returned = TAIR.get(self.test_agis, "upstream_1000", "rep_gene")
        self._perform_test(returned, expected)

    def test_upstream3000(self):
        expected = [
            ("AT4G36450", "TAGAAAATTA", 3000),
            ("AT4G36900", "GGTTTAGTTA", 3000),
            ("AT5G63980", "TGATTAATAC", 3000)
            ]
        returned = TAIR.get(self.test_agis, "upstream_3000", "rep_gene")
        self._perform_test(returned, expected)

    def test_downstream500(self):
        expected = [
            ("AT4G36450", "ATCAATATTT", 500),
            ("AT4G36900", "AATCCCCCTC", 500),
            ("AT5G63980", "GTATAAAGGG", 500)
            ]
        returned = TAIR.get(self.test_agis, "downstream_500", "rep_gene")
        self._perform_test(returned, expected)

    def test_downstream1000(self):
        expected = [
            ("AT4G36450", "ATCAATATTT", 1000),
            ("AT4G36900", "AATCCCCCTC", 1000),
            ("AT5G63980", "GTATAAAGGG", 1000)
            ]
        returned = TAIR.get(self.test_agis, "downstream_1000", "rep_gene")
        self._perform_test(returned, expected)

    def test_downstream3000(self):
        expected = [
            ("AT4G36450", "ATCAATATTT", 3000),
            ("AT4G36900", "AATCCCCCTC", 3000),
            ("AT5G63980", "GTATAAAGGG", 3000)
            ]
        returned = TAIR.get(self.test_agis, "downstream_3000", "rep_gene")
        self._perform_test(returned, expected)

    def test_empty_return(self):
        # This checks to see that no errors are raised when a invalid AGI is
        # given.
        returned = TAIR.get(["AT5G66483.1"], "downstream_3000", "rep_gene")
        # self.assertTrue(len(returned)) doesn't work with generators, so do:
        try:
            returned.next()
        except StopIteration:
            # This should happen, it means generator is empty
            self.assertTrue(True)
            return
        self.assertTrue(False)  # This shouldnt happen

    def test_agiregex(self):
        agis = ["AT5G63980", "AT5G63980.1", "ATCG12345", "ATCG12345.6"]
        self.assertEqual(agis, TAIR._sanitise_agis(agis))
        not_agis = ["notanagi", "not_an_agi", "\as23sd1321\ '.XS"]
        self.assertEqual([], TAIR._sanitise_agis(not_agis))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
