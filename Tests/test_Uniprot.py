#!/usr/bin/env python
"""Test for the Uniprot parser on Uniprot XML files.
"""
import os
import copy
import unittest

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from seq_tests_common import compare_reference

class TestUniprot(unittest.TestCase):

    def test_uni001(self):
        "Parsing Uniprot file uni001"
        filename = 'uni001'
        # test the record parser

        datafile = os.path.join('SwissProt', filename)

        test_handle = open(datafile)
        seq_record = SeqIO.read(test_handle, "uniprot")
        test_handle.close()

        self.assertTrue(isinstance(seq_record, SeqRecord))

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(seq_record.id, "Q91G55")
        self.assertEqual(seq_record.name, "043L_IIV6")
        self.assertEqual(seq_record.description, "Uncharacterized protein 043L")
        self.assertEqual(repr(seq_record.seq), "Seq('MDLINNKLNIEIQKFCLDLEKKYNINYNNLIDLWFNKESTERLIKCEVNLENKI...IPI', ProteinAlphabet())")

        # self.assertEqual(seq_record.accessions, ['Q91G55']) #seq_record.accessions does not exist
        # self.assertEqual(seq_record.organism_classification, ['Eukaryota', 'Metazoa', 'Chordata', 'Craniata', 'Vertebrata', 'Mammalia', 'Eutheria', 'Primates', 'Catarrhini', 'Hominidae', 'Homo'])
        # self.assertEqual(record.seqinfo, (348, 39676, '75818910'))
    
        self.assertEqual(len(seq_record.features), 1)
        self.assertEqual(repr(seq_record.features[0]), "SeqFeature(FeatureLocation(ExactPosition(0),ExactPosition(115)), type='chain', id='PRO_0000377969')")

        self.assertEqual(len(seq_record.annotations['references']), 2)
        self.assertEqual(seq_record.annotations['references'][0].authors, 'Jakob N.J., Mueller K., Bahr U., Darai G.')
        self.assertEqual(seq_record.annotations['references'][0].title, 'Analysis of the first complete DNA sequence of an invertebrate iridovirus: coding strategy of the genome of Chilo iridescent virus.')
        self.assertEqual(seq_record.annotations['references'][0].journal, 'Virology 286:182-196(2001)')
        self.assertEqual(seq_record.annotations['references'][0].comment, 'journal article | 2001 | Scope: NUCLEOTIDE SEQUENCE [LARGE SCALE GENOMIC DNA] | ')

        self.assertEqual(len(seq_record.dbxrefs), 11)
        self.assertEqual(seq_record.dbxrefs[0], 'DOI:10.1006/viro.2001.0963')

        self.assertEqual(seq_record.annotations['sequence_length'], 116)
        self.assertEqual(seq_record.annotations['sequence_checksum'], '4A29B35FB716523C')
        self.assertEqual(seq_record.annotations['modified'], '2009-07-07')
        self.assertEqual(seq_record.annotations['accessions'], ['Q91G55'])
        self.assertEqual(seq_record.annotations['taxonomy'], ['Viruses', 'dsDNA viruses, no RNA stage', 'Iridoviridae', 'Iridovirus'])
        self.assertEqual(seq_record.annotations['sequence_mass'], 13673)
        self.assertEqual(seq_record.annotations['dataset'], 'Swiss-Prot')
        self.assertEqual(seq_record.annotations['gene_name_ORF'], ['IIV6-043L'])
        self.assertEqual(seq_record.annotations['version'], 21)
        self.assertEqual(seq_record.annotations['sequence_modified'], '2001-12-01')
        self.assertEqual(seq_record.annotations['keywords'], ['Complete proteome', 'Virus reference strain'])
        self.assertEqual(seq_record.annotations['organismHost_name'], ['Acheta domesticus', 'House cricket', 'Chilo suppressalis', 'striped riceborer', 'Gryllus bimaculatus', 'Two-spotted cricket', 'Gryllus campestris', 'Spodoptera frugiperda', 'Fall armyworm'])
        self.assertEqual(seq_record.annotations['created'], '2009-06-16')
        self.assertEqual(seq_record.annotations['organism_name'], ['Chilo iridescent virus'])
        self.assertEqual(seq_record.annotations['organism'], 'Invertebrate iridescent virus 6 (IIV-6)')
        self.assertEqual(seq_record.annotations['recommendedName_fullName'], ['Uncharacterized protein 043L'])
        self.assertEqual(seq_record.annotations['sequence_version'], 1)
        self.assertEqual(seq_record.annotations['proteinExistence'], ['Predicted'])

    def compare_txt_xml(self, old, new):
        self.assertEqual(old.id, new.id)
        self.assertEqual(old.name, new.name)
        self.assertEqual(len(old), len(new))
        self.assertEqual(str(old.seq), str(new.seq))
        for key in set(old.annotations).intersection(new.annotations):
            if key == "references":
                self.assertEqual(len(old.annotations[key]),
                                 len(new.annotations[key]))
                for r1, r2 in zip(old.annotations[key], new.annotations[key]):
                    #Tweak for line breaks in plain text SwissProt
                    r1.title = r1.title.replace("- ", "-")
                    r2.title = r2.title.replace("- ", "-")
                    r1.journal = r1.journal.rstrip(".") #Should parser do this?
                    r1.medline_id = "" #Missing in UniPort MXL? TODO - check
                    #Lots of extra comments in UniProt XML
                    r1.comment = ""
                    r2.comment = ""
                    if not r2.journal: r1.journal = ""
                    compare_reference(r1, r2)
	    elif old.annotations[key] == new.annotations[key]:
		pass
	    elif key in ["date"]:
		#TODO - Why is this a list vs str?
		pass
	    elif type(old.annotations[key]) != type(new.annotations[key]):
		raise TypeError("%s gives %s vs %s" % \
				 (key, old.annotations[key], new.annotations[key]))
	    elif key in ["organism"]:
		if old.annotations[key] == new.annotations[key]:
		    pass
		elif old.annotations[key].startswith(new.annotations[key]+" "):
		    pass
		else:
		    raise ValueError(key)
	    elif isinstance(old.annotations[key], list) \
	    and sorted(old.annotations[key]) == sorted(new.annotations[key]):
		pass
            else:
		raise ValueError("%s gives %s vs %s" % \
				 (key, old.annotations[key], new.annotations[key]))
        #TODO - Parse features in plain text, and compare those

    def test_Q13639(self):
	"""Compare SwissProt text and uniprot XML versions of Q13639."""
	old = SeqIO.read("SwissProt/Q13639.txt", "swiss")
        new = SeqIO.read("SwissProt/Q13639.xml", "uniprot")
	self.compare_txt_xml(old, new)
    
    def test_swiss4(self):
	"""Compare SwissProt text and uniprot XML versions of 4 entries."""
	old = list(SeqIO.parse("SwissProt/swiss4.txt", "swiss"))
	new = list(SeqIO.parse("SwissProt/swiss4.xml", "uniprot"))
	fasta = list(SeqIO.parse("SwissProt/swiss4.fasta", "fasta"))
	ids = [x.strip() for x in open("SwissProt/swiss4.list")]
	self.assertEqual(len(old), len(fasta))
	self.assertEqual(len(old), len(ids))
	self.assertEqual(len(old), len(new))
	for txt, xml, fas, id in zip(old, new, fasta, ids):
	    self.assertEqual(txt.id, id)
	    self.assertTrue(txt.id in fas.id.split("|"))
	    self.assertEqual(str(txt.seq), str(fas.seq))
	    self.compare_txt_xml(txt, xml)

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
