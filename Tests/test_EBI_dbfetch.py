# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for online functionality of EBI dbfetch module."""

# Builtins
import unittest

from Bio.EBI.debfetch import db_info, dbfetch
from Bio import SeqIO


class dbfetchTests(unittest.TestCase):
    """Tests for EBI dbfetch REST API."""

    def test_dbinfo(self):
        """Test the meta-info retriever method."""
        h = db_info()
        self.assertEqual(h.url, "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/dbfetch.databases?style=json")

    def test_dbfetch_edam(self):
        """Query EDAM database via module method."""
        h1 = dbfetch(db="edam", ids=["0338", "1929"], format="obo")
        h2 = dbfetch(db="edam", ids="0338,1929", format="obo")
        self.assertEqual(h1.url, "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/edam/0338,1929/obo")
        self.assertEqual(h1.read(), h2.read())

    def test_dbfetch_ena_coding(self):
        """Query ENA Coding database via module method."""
        h1 = dbfetch(db="ena_coding", ids="AAA59452", format="annot")
        self.assertEqual(h1.url, "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/ena_coding/AAA59452/annot")

    def test_dbfetch_ensembl_gene(self):
        """Query Ensembl Gene database via module method."""
        h = dbfetch("ensembl_gene", "ENSBTAG00000000988", "fasta")
        data = SeqIO.read(h, "fasta")
        self.assertEqual(len(str(data.seq)), 50988)
        self.assertEqual(str(data.seq)[:100], 'ATGCCGATTGGATGCAAAGAGAGGCCAACTTTTTTTGACATTTTTAAGGCGCGATGCAACAAAGCAGGTATTAACAGATTTTTTATACGAATTGATCGTT')

        h = dbfetch("ensembl_gene", ["ENSG00000139618", "ENSMUSG00000041147"], "fasta")
        data = list(SeqIO.parse(h, "fasta"))
        self.assertEqual(len(str(data[0].seq)), 84793)
        self.assertEqual(str(data[0].seq)[:100], 'GGGCTTGTGGCGCGAGCTTCTGAAACTAGGCGGCAGAGGCGGAGCCGCTGTGGCACTGCTGCGCCTCTGCTGCGCCTCGGGTGTCTTTTGCGGCGGTGGG')
        self.assertEqual(len(str(data[1].seq)), 47700)
        self.assertEqual(str(data[1].seq)[:100], 'GCGGGAGCGGGAGCCGTGAGGCGTTGCCGTCAGTCAGCTACCGCTGCGGGAGCGGAGCGGGTCGGTGCGGCCGGGTGTGGCGGGCGTGCGCTCCGGGGTC')

    def test_dbfetch_uniprotkb(self):
        """Query UniProtKB database via module method."""
        h = dbfetch("uniprotkb", "P06213", "fasta")
        data = SeqIO.read(h, "fasta")
        self.assertEqual(len(str(data.seq)), 1382)
        self.assertEqual(str(data.seq)[:100], 'MATGGRRGAAAAPLLVAVAALLLGAAGHLYPGEVCPGMDIRNNLTRLHELENCSVIEGHLQILLMFKTRPEDFRDLSFPKLIMITDYLLLFRVYGLESLK')

        h = dbfetch("uniprotkb", ["1433X_MAIZE", "1433T_RAT"], "fasta")
        data = list(SeqIO.parse(h, "fasta"))
        self.assertEqual(len(str(data[0].seq)), 61)
        self.assertEqual(str(data[0].seq)[:100], 'ILNSPDRACNLAKQAFDEAISELDSLGEESYKDSTLIMQLLXDNLTLWTSDTNEDGGDEIK')
        self.assertEqual(len(str(data[1].seq)), 245)
        self.assertEqual(str(data[1].seq)[:100], 'MEKTELIQKAKLAEQAERYDDMATCMKAVTEQGAELSNEERNLLSVAYKNVVGGRRSAWRVISSIEQKTDTSDKKLQLIKDYREKVESELRSICTTVLEL')

    def test_dbfetch_pdb(self):
        """Query PDB database via module method."""
        h = dbfetch("pdb", "1GAG", "pdb")
        data = list(SeqIO.parse(h, "pdb-atom"))
        self.assertEqual(len(str(data[0].seq)), 303)
        self.assertEqual(str(data[0].seq)[:100], 'SSVFVPDEWEVSREKITLLRELGQGSFGMVYEGNARDIIKGEAETRVAVKTVNESASLRERIEFLNEASVMKGFTCHHVVRLLGVVSKGQPTLVVMELMA')
        self.assertEqual(len(str(data[1].seq)), 13)
        self.assertEqual(str(data[1].seq)[:100], 'PATGDFMNMSPVG')

    def test_dbfetch_interpro(self):
        """Query InterPro database via module method."""
        h = dbfetch("interpro", "IPR006212", "tab", "raw")
        for line in h.readlines():
            if line[0] != '#':
                cols = line.split('\t')
                self.assertEqual(cols[0], "IPR006212")
                self.assertEqual(cols[1], "Repeat")
                self.assertEqual(cols[2], "Furin_repeat")
                self.assertEqual(cols[3][:-1], "Furin-like repeat")
        h = dbfetch("interpro", "IPR008958,IPR009030", "tab", "raw")
        data = h.readlines()
        cols = data[2].split('\t')
        self.assertEqual(cols[0], "IPR008958")
        self.assertEqual(cols[1], "Domain")
        self.assertEqual(cols[2], "Transglutaminase_C")
        self.assertEqual(cols[3][:-1], "Transglutaminase, C-terminal")
        cols = data[3].split('\t')
        self.assertEqual(cols[0], "IPR009030")
        self.assertEqual(cols[1], "Domain")
        self.assertEqual(cols[2], "Growth_fac_rcpt_")
        self.assertEqual(cols[3][:-1], "Growth factor receptor cysteine-rich domain")


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
