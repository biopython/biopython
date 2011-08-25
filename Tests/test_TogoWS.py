# Copyright 2010-2011 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Testing Bio.TogoWS online code.
"""
import sys
import unittest
import urllib

#import requires_internet
#requires_internet.check()

#We want to test these:
from Bio import TogoWS

#In order to check any sequences returned
from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid
from Bio import Medline

#####################################################################

class TogoFields(unittest.TestCase):
    def test_databases(self):
        """Check supported databases"""
        dbs = set(TogoWS._get_entry_dbs())
        self.assert_(dbs.issuperset(['nuccore', 'nucest', 'nucgss',
                                     'nucleotide', 'protein', 'gene',
                                     'omim', 'homologene', 'snp',
                                     'mesh', 'pubmed', 'embl',
                                     'uniprot', 'uniparc', 'uniref100',
                                     'uniref90', 'uniref50', 'ddbj',
                                     'dad', 'pdb', 'compound', 'drug',
                                     'enzyme', 'genes', 'glycan',
                                     'orthology', 'reaction', 'module',
                                     'pathway']), dbs)

    def test_pubmed(self):
        """Check supported fields for pubmed database"""
        fields = set(TogoWS._get_entry_fields("pubmed"))
        self.assert_(fields.issuperset(["abstract", "au", "authors", "doi",
                                        "mesh", "so"]), fields)

    def test_ncbi_protein(self):
        """Check supported fields for NCBI protein database"""
        fields = set(TogoWS._get_entry_fields("ncbi-protein"))
        self.assert_(fields.issuperset(['entry_id', 'length', 'strand',
                                        'moltype', 'linearity', 'division',
                                        'date', 'definition', 'accession',
                                        'accessions', 'version', 'versions',
                                        'acc_version', 'gi', 'keywords',
                                        'organism', 'common_name',
                                        'taxonomy', 'comment', 'seq']),
                                        fields)

    def test_ddbj(self):
        """Check supported fields for ddbj database"""
        fields = set(TogoWS._get_entry_fields("ddbj"))
        self.assert_(fields.issuperset(['entry_id', 'length', 'strand',
                                        'moltype', 'linearity', 'division',
                                        'date', 'definition', 'accession',
                                        'accessions', 'version', 'versions',
                                        'acc_version', 'gi', 'keywords',
                                        'organism', 'common_name',
                                        'taxonomy', 'comment', 'seq']),
                                        fields)

    def test_embl(self):
        """Check supported fields for embl database"""
        fields = set(TogoWS._get_entry_fields("embl"))
        self.assert_(fields.issuperset(["definition", "entry_id", "seq"]),
                     fields)

    def test_uniprot(self):
        """Check supported fields for uniprot database"""
        fields = set(TogoWS._get_entry_fields("uniprot"))
        self.assert_(fields.issuperset(["definition", "entry_id", "seq"]),
                     fields)

    def test_pdb(self):
        """Check supported fields for pdb database"""
        fields = set(TogoWS._get_entry_fields("pdb"))
        self.assert_(fields.issuperset(["accession", "chains", "keywords",
                                        "models"]), fields)

class TogoTests(unittest.TestCase):
    def test_pubmed_16381885(self):
        """Bio.TogoWS.entry("pubmed", "16381885")"""
        #Gives Medline plain text
        handle = TogoWS.entry("pubmed", "16381885")
        data = Medline.read(handle)
        handle.close()
        self.assertEqual(data["TI"],
             'From genomics to chemical genomics: new developments in KEGG.')
        self.assertEqual(data["AU"], ['Kanehisa M', 'Goto S', 'Hattori M',
                                      'Aoki-Kinoshita KF', 'Itoh M',
                                      'Kawashima S', 'Katayama T', 'Araki M',
                                      'Hirakawa M'])

    def test_pubmed_16381885_ti(self):
        """Bio.TogoWS.entry("pubmed", "16381885", field="ti")"""
        handle = TogoWS.entry("pubmed", "16381885", field="ti")
        data = handle.read().strip()
        handle.close()
        self.assertEqual(data,
             'From genomics to chemical genomics: new developments in KEGG.')
    
    def test_pubmed_16381885_title(self):
        """Bio.TogoWS.entry("pubmed", "16381885", field="title")"""
        handle = TogoWS.entry("pubmed", "16381885", field="title")
        data = handle.read().strip()
        handle.close()
        self.assertEqual(data,
             'From genomics to chemical genomics: new developments in KEGG.')
    
    def test_pubmed_16381885_au(self):
        """Bio.TogoWS.entry("pubmed", "16381885", field="au")"""
        #Gives one name per line (i.e. \n separated), no dots
        handle = TogoWS.entry("pubmed", "16381885", field="au")
        data = handle.read().strip().split("\n")
        handle.close()
        self.assertEqual(data, ['Kanehisa M', 'Goto S', 'Hattori M',
                                'Aoki-Kinoshita KF', 'Itoh M',
                                'Kawashima S', 'Katayama T', 'Araki M',
                                'Hirakawa M'])

    def test_pubmed_16381885_authors(self):
        """Bio.TogoWS.entry("pubmed", "16381885", field="authors")"""
        #Gives names tab separated (i.e. \t separated)
        handle = TogoWS.entry("pubmed", "16381885", field="authors")
        data = handle.read().strip().split("\t")
        handle.close()
        self.assertEqual(data, ['Kanehisa, M.', 'Goto, S.', 'Hattori, M.',
                                'Aoki-Kinoshita, K. F.', 'Itoh, M.',
                                'Kawashima, S.', 'Katayama, T.', 'Araki, M.',
                                'Hirakawa, M.'])

    def test_pubmed_16381885_invalid_field(self):
        """Bio.TogoWS.entry("pubmed", "16381885", field="invalid_for_testing")"""
        self.assertRaises(IOError, TogoWS.entry,
                          "pubmed", "16381885", field="invalid_for_testing")

    def test_pubmed_16381885_invalid_format(self):
        """Bio.TogoWS.entry("pubmed", "16381885", format="invalid_for_testing")"""
        self.assertRaises(ValueError, TogoWS.entry,
                          "pubmed", "16381885", format="invalid_for_testing")

    def test_pubmed_invalid_id(self):
        """Bio.TogoWS.entry("pubmed", "invalid_for_testing")"""
        self.assertRaises(IOError, TogoWS.entry,
                          "pubmed", "invalid_for_testing")

    def test_invalid_db(self):
        """Bio.TogoWS.entry("invalid_db", "invalid_id")"""
        self.assertRaises(ValueError, TogoWS.entry,
                          "invalid_db", "invalid_id")

    def test_ddbj_genbank_length(self):
        """Bio.TogoWS.entry("genbank", "NC_000913.2", field="length")"""
        handle = TogoWS.entry("genbank", "NC_000913.2", field="length")
        data = handle.read().strip() #ignore trailing \n
        handle.close()
        self.assertEqual(data, "4639675")

    def test_ddbj_genbank(self):
        """Bio.TogoWS.entry("ddbj", "X52960")"""
        handle = TogoWS.entry("ddbj", "X52960") #Returns "genbank" format
        record  = SeqIO.read(handle, "gb")
        handle.close()
        self.assertEqual(record.id, "X52960.1")
        self.assertEqual(record.name, "X52960")
        self.assertEqual(len(record), 248)
        self.assertEqual(seguid(record.seq), "Ktxz0HgMlhQmrKTuZpOxPZJ6zGU")

    def test_ddbj_genbank_length(self):
        """Bio.TogoWS.entry("ddbj", "X52960", field="length")"""
        handle = TogoWS.entry("ddbj", "X52960", field="length")
        data = handle.read().strip() #ignore trailing \n
        handle.close()
        self.assertEqual(data, "248")

    def test_ddbj_genbank_seq(self):
        """Bio.TogoWS.entry("ddbj", "X52960", field="seq")"""
        handle = TogoWS.entry("ddbj", "X52960", field="seq")
        data = handle.read().strip() #ignore trailing \n
        handle.close()
        self.assertEqual(seguid(data), "Ktxz0HgMlhQmrKTuZpOxPZJ6zGU")

    def test_ddbj_genbank_definition(self):
        """Bio.TogoWS.entry("ddbj", "X52960", field="definition")"""
        handle = TogoWS.entry("ddbj", "X52960", field="definition")
        data = handle.read().strip() #ignore trailing \n
        handle.close()
        self.assertEqual(data, "Coleus blumei viroid 1 (CbVd) RNA.")

    def test_ddbj_genbank_accession(self):
        """Bio.TogoWS.entry("ddbj", "X52960", field="accession")"""
        handle = TogoWS.entry("ddbj", "X52960", field="accession")
        data = handle.read().strip() #ignore trailing \n
        handle.close()
        self.assertEqual(data, "X52960")

    def test_ddbj_genbank_accession(self):
        """Bio.TogoWS.entry("ddbj", "X52960", field="version")"""
        handle = TogoWS.entry("ddbj", "X52960", field="version")
        data = handle.read().strip() #ignore trailing \n
        handle.close()
        self.assertEqual(data, "1")

    def test_ddbj_genbank_acc_version(self):
        """Bio.TogoWS.entry("ddbj", "X52960", field="acc_version")"""
        handle = TogoWS.entry("ddbj", "X52960", field="acc_version")
        data = handle.read().strip() #ignore trailing \n
        handle.close()
        self.assertEqual(data, "X52960.1")

    def test_ddbj_genbank_organism(self):
        """Bio.TogoWS.entry("ddbj", "X52960", field="organism")"""
        handle = TogoWS.entry("ddbj", "X52960", field="organism")
        data = handle.read().strip() #ignore trailing \n
        handle.close()
        self.assertEqual(data, "Coleus blumei viroid 1")
        
    def test_ddbj_genbank_invalid_field(self):
        """Bio.TogoWS.entry("ddbj", "X52960", field="invalid_for_testing")"""
        self.assertRaises(IOError, TogoWS.entry,
                          "ddbj", "X52960", field="invalid_for_testing")

    def test_ddbj_fasta(self):
        """Bio.TogoWS.entry("ddbj", "X52960", "fasta")"""
        handle = TogoWS.entry("ddbj", "X52960", "fasta")
        record  = SeqIO.read(handle, "fasta")
        handle.close()
        self.assert_("X52960" in record.id, record.id)
        self.assert_("X52960" in record.name, record.name)
        self.assertEqual(len(record), 248)
        self.assertEqual(seguid(record.seq), "Ktxz0HgMlhQmrKTuZpOxPZJ6zGU")

    def test_uniprot_swiss(self):
        """Bio.TogoWS.entry("uniprot", ["A1AG1_HUMAN","A1AG1_MOUSE"])"""
        #Returns "swiss" format:
        handle = TogoWS.entry("uniprot", ["A1AG1_HUMAN","A1AG1_MOUSE"])
        record1, record2  = SeqIO.parse(handle, "swiss")
        handle.close()

        self.assertEqual(record1.id, "P02763")
        self.assertEqual(record1.name, "A1AG1_HUMAN")
        self.assertEqual(len(record1), 201)
        self.assertEqual(seguid(record1.seq), "LHDJJ6oC7gUXo8CC7Xn6EUeA8Gk")

        self.assertEqual(record2.id, "Q60590")
        self.assertEqual(record2.name, "A1AG1_MOUSE")
        self.assertEqual(len(record2), 207)
        self.assertEqual(seguid(record2.seq), "FGcj+RFQhP2gRusCmwPFty5PJT0")

    def test_nucleotide_fasta(self):
        """Bio.TogoWS.entry("nucleotide", "6273291", "fasta")"""
        handle = TogoWS.entry("nucleotide", "6273291", "fasta")
        record  = SeqIO.read(handle, "fasta")
        handle.close()
        self.assert_("6273291" in record.id, record.id)
        self.assert_("6273291" in record.name, record.name)
        self.assertEqual(len(record), 902)
        self.assertEqual(seguid(record.seq), "bLhlq4mEFJOoS9PieOx4nhGnjAQ")

    def test_protein_fasta(self):
        """Bio.TogoWS.entry("protein", "16130152", "fasta")"""
        handle = TogoWS.entry("protein", "16130152", "fasta")
        record  = SeqIO.read(handle, "fasta")
        handle.close()
        #Could use assertIn but requires Python 2.7+
        self.assert_("16130152" in record.id, record.id)
        self.assert_("16130152" in record.name, record.name)
        self.assert_("porin protein" in record.description, record.description)
        self.assertEqual(len(record), 367)
        self.assertEqual(seguid(record.seq), "fCjcjMFeGIrilHAn6h+yju267lg")

class TogoConvert(unittest.TestCase):
    """Conversion tests."""
    
    def test_genbank_to_fasta(self):
        """Conversion of GenBank to FASTA."""
        filename = "GenBank/NC_005816.gb"
        old = SeqIO.read(filename, "gb")
        new = SeqIO.read(TogoWS.tconvert(open(filename), "genbank", "fasta"), "fasta")
        self.assertEqual(str(old.seq), str(new.seq))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)

    
