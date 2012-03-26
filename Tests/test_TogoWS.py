# Copyright 2010-2011 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Testing Bio.TogoWS online code.
"""
import sys
import unittest
import urllib
from StringIO import StringIO

import requires_internet
requires_internet.check()

#We want to test these:
from Bio import TogoWS

#In order to check any sequences returned
from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid
from Bio import Medline

#####################################################################

class TogoFields(unittest.TestCase):
    def test_invalid_database(self):
        """Check asking for fields of invalid database fails"""
        self.assertRaises(IOError, TogoWS._get_fields,
                          "http://togows.dbcls.jp/entry/invalid?fields")

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
        self.assert_(fields.issuperset(['abstract', 'au', 'authors',
                                        'doi', 'mesh', 'so', 'ti',
                                        'title']), fields)

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


class TogoEntry(unittest.TestCase):
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
        self.assertRaises(ValueError, TogoWS.entry,
                          "pubmed", "16381885", field="invalid_for_testing")

    def test_pubmed_16381885_invalid_format(self):
        """Bio.TogoWS.entry("pubmed", "16381885", format="invalid_for_testing")"""
        self.assertRaises(ValueError, TogoWS.entry,
                          "pubmed", "16381885", format="invalid_for_testing")

    def test_pubmed_invalid_id(self):
        """Bio.TogoWS.entry("pubmed", "invalid_for_testing")"""
        self.assertRaises(IOError, TogoWS.entry,
                          "pubmed", "invalid_for_testing")

    def test_pubmed_16381885_and_19850725(self):
        """Bio.TogoWS.entry("pubmed", "16381885,19850725")"""
        handle = TogoWS.entry("pubmed", "16381885,19850725")
        records = list(Medline.parse(handle))
        handle.close()
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0]["TI"],
             'From genomics to chemical genomics: new developments in KEGG.')
        self.assertEqual(records[0]["AU"], ['Kanehisa M', 'Goto S',
                                            'Hattori M', 'Aoki-Kinoshita KF',
                                            'Itoh M', 'Kawashima S',
                                            'Katayama T', 'Araki M',
                                            'Hirakawa M'])
        self.assertEqual(records[1]["TI"],
             'DDBJ launches a new archive database with analytical tools ' + \
             'for next-generation sequence data.')
        self.assertEqual(records[1]["AU"], ['Kaminuma E', 'Mashima J',
                                            'Kodama Y', 'Gojobori T',
                                            'Ogasawara O', 'Okubo K',
                                            'Takagi T', 'Nakamura Y'])

    def test_pubmed_16381885_and_19850725_authors(self):
        """Bio.TogoWS.entry("pubmed", "16381885,19850725", field="authors")"""
        handle = TogoWS.entry("pubmed", "16381885,19850725", field="authors")
        #Little hack to remove blank lines...
        #names = handle.read().replace("\n\n", "\n").strip().split("\n")
        names = handle.read().strip().split("\n")
        handle.close()
        self.assertEqual(2, len(names))
        names1, names2 = names
        self.assertEqual(names1.split("\t"),
                         ['Kanehisa, M.', 'Goto, S.', 'Hattori, M.',
                          'Aoki-Kinoshita, K. F.', 'Itoh, M.',
                          'Kawashima, S.', 'Katayama, T.',
                          'Araki, M.', 'Hirakawa, M.'])
        self.assertEqual(names2.split("\t"),
                         ['Kaminuma, E.', 'Mashima, J.', 'Kodama, Y.',
                          'Gojobori, T.', 'Ogasawara, O.', 'Okubo, K.',
                          'Takagi, T.', 'Nakamura, Y.'])

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
        self.assertRaises(ValueError, TogoWS.entry,
                          "ddbj", "X52960", field="invalid_for_testing")

    def test_ddbj_invalid_format(self):
        """Bio.TogoWS.entry("ddbj", "X52960", format="invalid_for_testing")"""
        self.assertRaises(ValueError, TogoWS.entry,
                          "ddbj", "X52960", format="invalid_for_testing")

    def test_ddbj_gff3(self):
        """Bio.TogoWS.entry("ddbj", "X52960", format="gff")"""
        handle = TogoWS.entry("ddbj", "X52960", format="gff")
        data = handle.read()
        handle.close()
        self.assert_(data.startswith("##gff-version 3\nX52960\tDDBJ\t"), data)

    def test_genbank_gff3(self):
        """Bio.TogoWS.entry("nucleotide", "X52960", format="gff")"""
        #Note - Using manual URL with genbank instead of nucleotide works
        handle = TogoWS.entry("nucleotide", "X52960", format="gff")
        data = handle.read()
        handle.close()
        self.assert_(data.startswith("##gff-version 3\nX52960\tGenbank\t"), data)

    def test_embl_AM905444_gff3(self):
        """Bio.TogoWS.entry("embl", "AM905444", format="gff")"""
        handle = TogoWS.entry("embl", "AM905444", format="gff")
        data = handle.read()
        handle.close()
        self.assert_(data.startswith("##gff-version 3\nAM905444\tembl\t"), data)

    def test_embl_AM905444_seq(self):
        """Bio.TogoWS.entry("embl", "AM905444", field="seq")"""
        handle = TogoWS.entry("embl", "AM905444", field="seq")
        data = handle.read().strip() #ignore any trailing \n
        handle.close()
        self.assertEqual(seguid(data), "G0HtLpwF7i4FXUaUjDUPTjok79c")

    def test_embl_AM905444_definition(self):
        """Bio.TogoWS.entry("embl", "AM905444", field="definition")"""
        handle = TogoWS.entry("embl", "AM905444", field="definition")
        data = handle.read().strip() #ignore any trailing \n
        handle.close()
        self.assertEqual(data, "Herbaspirillum seropedicae locus tag HS193.0074 for porin")

    def test_embl_AM905444(self):
        """Bio.TogoWS.entry("embl", "AM905444")"""
        handle = TogoWS.entry("embl", "AM905444")
        record = SeqIO.read(handle, "embl")
        handle.close()
        self.assert_("AM905444" in record.id, record.id)
        self.assert_("AM905444" in record.name, record.name)
        self.assert_("porin" in record.description, record.description)
        self.assertEqual(len(record), 1164)
        self.assertEqual(seguid(record.seq), "G0HtLpwF7i4FXUaUjDUPTjok79c")

    def test_ddbj_fasta(self):
        """Bio.TogoWS.entry("ddbj", "X52960", "fasta")"""
        handle = TogoWS.entry("ddbj", "X52960", "fasta")
        record = SeqIO.read(handle, "fasta")
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

class TogoSearch(unittest.TestCase):
    """Search tests."""

    def test_bad_args_just_limit(self):
        """Reject Bio.TogoWS.search(...) with just limit"""
        self.assertRaises(ValueError, TogoWS.search,
                          "pubmed", "lung+cancer", limit=10)

    def test_bad_args_just_offset(self):
        """Reject Bio.TogoWS.search(...) with just offset"""
        self.assertRaises(ValueError, TogoWS.search,
                          "pubmed", "lung+cancer", offset=10)

    def test_bad_args_zero_limit(self):
        """Reject Bio.TogoWS.search(...) with zero limit"""
        self.assertRaises(ValueError, TogoWS.search,
                          "pubmed", "lung+cancer", offset=1, limit=0)

    def test_bad_args_zero_offset(self):
        """Reject Bio.TogoWS.search(...) with zero offset"""
        self.assertRaises(ValueError, TogoWS.search,
                          "pubmed", "lung+cancer", offset=0, limit=10)

    def test_bad_args_non_int_offset(self):
        """Reject Bio.TogoWS.search(...) with non-integer offset"""
        self.assertRaises(ValueError, TogoWS.search,
                          "pubmed", "lung+cancer", offset="test", limit=10)

    def test_bad_args_non_int_limit(self):
        """Reject Bio.TogoWS.search(...) with non-integer limit"""
        self.assertRaises(ValueError, TogoWS.search,
                          "pubmed", "lung+cancer", offset=1, limit="lots")

    def test_pubmed_search_togows(self):
        """Bio.TogoWS.search_iter("pubmed", "TogoWS") etc"""
        self.check("pubmed", "TogoWS", ["20472643"])

    def test_pubmed_search_bioruby(self):
        """Bio.TogoWS.search_iter("pubmed", "BioRuby") etc"""
        self.check("pubmed", "BioRuby", ["20739307", "20015970", "14693808"])

    def test_pubmed_search_porin(self):
        """Bio.TogoWS.search_iter("pubmed", "human porin") etc

        Count was 339 at time of writing, this was choosen to
        be larger than the default chunk size for iteration,
        but still not too big to download the full list.
        """
        self.check("pubmed", "human porin", ["21189321", "21835183"])

    def test_pdb_search_porin(self):
        """Bio.TogoWS.search_iter("pdb", "porin") etc

        Count was about 130 at time of writing.
        """
        self.check("pdb", "porin", ["2j1n", "2vqg", "3m8b", "2k0l"])

    def test_embl_search_porin(self):
        """Bio.TogoWS.search_iter("embl", "human pore", limit=200) etc

        Count was about 255 at time of writing.
        """
        self.check("embl", "human pore", limit=200)

    def test_uniprot_search_lung_cancer(self):
        """Bio.TogoWS.search_iter("uniprot", "lung+cancer", limit=150) etc

        Search count was 1327 at time of writing, a bit large to
        download all the results in a unit test. Want to use a limit
        larger than the batch size (100) to ensure at least two
        batches.
        """
        self.check("uniprot", "lung+cancer", limit=150)

    def check(self, database, search_term, expected_matches=[], limit=None):
        if expected_matches and limit:
            raise ValueError("Bad test - TogoWS makes no promises about order")
        search_count = TogoWS.search_count(database, search_term)
        if expected_matches and search_count < len(expected_matches):
            raise ValueError("Only %i matches, expected at least %i" \
                             % (search_count, len(expected_matches)))
        if search_count > 5000 and not limit:
            print "%i results, skipping" % search_count
            return
        if limit:
            count = min(search_count, limit)
        else:
            count = search_count

        #Iteration should find everything... unless a limit is used
        search_iter = list(TogoWS.search_iter(database, search_term, limit))
        self.assertEqual(count, len(search_iter))
        for match in expected_matches:
            self.assert_(match in search_iter,
                         "Expected %s in results but not" % match)


class TogoConvert(unittest.TestCase):
    """Conversion tests."""

    def test_invalid_format(self):
        """Check convert file format checking."""
        self.assertRaises(ValueError, TogoWS.convert,
                          StringIO("PLACEHOLDER"),
                          "genbank", "invalid_for_testing")
        self.assertRaises(ValueError, TogoWS.convert,
                          StringIO("PLACEHOLDER"),
                          "invalid_for_testing", "fasta")

    def test_genbank_to_fasta(self):
        """Conversion of GenBank to FASTA."""
        filename = "GenBank/NC_005816.gb"
        old = SeqIO.read(filename, "gb")
        new = SeqIO.read(TogoWS.convert(open(filename), "genbank", "fasta"), "fasta")
        self.assertEqual(str(old.seq), str(new.seq))

    def test_genbank_to_embl(self):
        """Conversion of GenBank to EMBL."""
        filename = "GenBank/NC_005816.gb"
        old = SeqIO.read(filename, "gb")
        new = SeqIO.read(TogoWS.convert(open(filename), "genbank", "embl"), "embl")
        self.assertEqual(str(old.seq), str(new.seq))

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)

    
