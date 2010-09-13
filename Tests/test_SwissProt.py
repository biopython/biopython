#!/usr/bin/env python
"""Test for the SwissProt parser on SwissProt files.
"""
import os
import unittest

from Bio import SeqIO
from Bio import SwissProt
from Bio.SeqRecord import SeqRecord



class TestSwissProt(unittest.TestCase):

    def test_sp001(self):
        "Parsing SwissProt file sp001"
        filename = 'sp001'
        # test the record parser

        datafile = os.path.join('SwissProt', filename)

        test_handle = open(datafile)
        seq_record = SeqIO.read(test_handle, "swiss")
        test_handle.close()

        self.assertTrue(isinstance(seq_record, SeqRecord))

        self.assertEqual(seq_record.id, "Q13454")
        self.assertEqual(seq_record.name, "N33_HUMAN")
        self.assertEqual(seq_record.description, "N33 PROTEIN.")
        self.assertEqual(repr(seq_record.seq), "Seq('MGARGAPSRRRQAGRRLRYLPTGSFPFLLLLLLLCIQLGGGQKKKENLLAEKVE...DFE', ProteinAlphabet())")

        test_handle = open(datafile)
        record = SwissProt.read(test_handle)
        test_handle.close()

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(record.entry_name, "N33_HUMAN")
        self.assertEqual(record.accessions, ['Q13454', 'Q14911', 'Q14912'])
        self.assertEqual(record.organism_classification, ['Eukaryota', 'Metazoa', 'Chordata', 'Craniata', 'Vertebrata', 'Mammalia', 'Eutheria', 'Primates', 'Catarrhini', 'Hominidae', 'Homo'])
        self.assertEqual(record.seqinfo, (348, 39676, '75818910'))
    
        self.assertEqual(len(record.features), 6)
        self.assertEqual(record.features[0], ('TRANSMEM', 20, 40, 'POTENTIAL.', ''))
        self.assertEqual(record.features[1], ('TRANSMEM', 197, 217, 'POTENTIAL.', ''))
        self.assertEqual(record.features[2], ('TRANSMEM', 222, 242, 'POTENTIAL.', ''))
        self.assertEqual(record.features[3], ('TRANSMEM', 277, 297, 'POTENTIAL.', ''))
        self.assertEqual(record.features[4], ('TRANSMEM', 313, 333, 'POTENTIAL.', ''))
        self.assertEqual(record.features[5], ('VARSPLIC', 344, 348, 'DLDFE -> FLIK (IN FORM 2).', ''))

        self.assertEqual(len(record.references), 1)
        self.assertEqual(record.references[0].authors, "MACGROGAN D., LEVY A., BOVA G.S., ISAACS W.B., BOOKSTEIN R.")
        self.assertEqual(record.references[0].title, "Structure and methylation-associated silencing of a gene within a homozygously deleted region of human chromosome band 8p22.")
        self.assertEqual(len(record.references[0].references), 1)
        self.assertEqual(record.references[0].references[0], ('MEDLINE', '96299740'))

        #Check the two parsers agree on the essentials
        self.assertEqual(seq_record.seq.tostring(), record.sequence)
        self.assertEqual(seq_record.description, record.description)
        self.assertEqual(seq_record.name, record.entry_name)
        self.assertTrue(seq_record.id in record.accessions)

        #Now try using the iterator - note that all these
        #test cases have only one record.

        # With the SequenceParser
        test_handle = open(datafile)
        records = list(SeqIO.parse(test_handle, "swiss"))
        test_handle.close()

        self.assertEqual(len(records), 1)
        self.assertTrue(isinstance(records[0], SeqRecord))

        #Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].seq.tostring(), seq_record.seq.tostring())
        self.assertEqual(records[0].description, seq_record.description)
        self.assertEqual(records[0].name, seq_record.name)
        self.assertEqual(records[0].id, seq_record.id)

        # With the RecordParser
        test_handle = open(datafile)
        records = list(SwissProt.parse(test_handle))
        test_handle.close()

        self.assertEqual(len(records), 1)
        self.assertTrue(isinstance(records[0], SwissProt.Record))

        #Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].sequence, record.sequence)
        self.assertEqual(records[0].description, record.description)
        self.assertEqual(records[0].entry_name, record.entry_name)
        self.assertEqual(records[0].accessions, record.accessions)


    def test_sp002(self):
        "Parsing SwissProt file sp002"

        filename = 'sp002'
        # test the record parser

        datafile = os.path.join('SwissProt', filename)

        test_handle = open(datafile)
        seq_record = SeqIO.read(test_handle, "swiss")
        test_handle.close()

        self.assertTrue(isinstance(seq_record, SeqRecord))

        self.assertEqual(seq_record.id, "P54101")
        self.assertEqual(seq_record.name, "CSP_MOUSE")
        self.assertEqual(seq_record.description, "CYSTEINE STRING PROTEIN (CSP).")
        self.assertEqual(repr(seq_record.seq), "Seq('MADQRQRSLSTSGESLYHVLGLDKNATSDDIKKSYRKLALKYHPDKNPDNPEAA...GFN', ProteinAlphabet())")

        test_handle = open(datafile)
        record = SwissProt.read(test_handle)
        test_handle.close()

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(record.entry_name, "CSP_MOUSE")
        self.assertEqual(record.accessions, ['P54101'])
        self.assertEqual(record.organism_classification, ['Eukaryota', 'Metazoa', 'Chordata', 'Craniata', 'Vertebrata', 'Mammalia', 'Eutheria', 'Rodentia', 'Sciurognathi', 'Muridae', 'Murinae', 'Mus'])
        self.assertEqual(record.seqinfo, (198, 22100, '9DF0142B'))
    
        self.assertEqual(len(record.features), 2) 
        self.assertEqual(record.features[0], ('DOMAIN', 13, 82, 'DNAJ-LIKE.', ''))
        self.assertEqual(record.features[1], ('DOMAIN', 118, 128, 'POLY-CYS.', ''))

        self.assertEqual(len(record.references), 3)
        self.assertEqual(record.references[0].authors, "QIN N., LIN T., BIRNBAUMER L.")
        self.assertEqual(record.references[0].title, "")
        self.assertEqual(len(record.references[0].references), 0)
        self.assertEqual(record.references[1].authors, "MASTROGIACOMO A., GUNDERSEN C.B.")
        self.assertEqual(record.references[1].title, "The nucleotide and deduced amino acid sequence of a rat cysteine string protein.")
        self.assertEqual(len(record.references[1].references), 1)
        self.assertEqual(record.references[1].references[0], ('MEDLINE', '95223109'))
        self.assertEqual(record.references[2].authors, "BRAUN J.E., SCHELLER R.H.")
        self.assertEqual(record.references[2].title, "Cysteine string protein, a DnaJ family member, is present on diverse secretory vesicles.")
        self.assertEqual(len(record.references[2].references), 1)
        self.assertEqual(record.references[2].references[0], ('MEDLINE', '96188189'))

        #Check the two parsers agree on the essentials
        self.assertEqual(seq_record.seq.tostring(), record.sequence)
        self.assertEqual(seq_record.description, record.description)
        self.assertEqual(seq_record.name, record.entry_name)
        self.assertTrue(seq_record.id in record.accessions)

        #Now try using the iterator - note that all these
        #test cases have only one record.

        # With the SequenceParser
        test_handle = open(datafile)
        records = list(SeqIO.parse(test_handle, "swiss"))
        test_handle.close()

        self.assertEqual(len(records), 1)
        self.assertTrue(isinstance(records[0], SeqRecord))

        #Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].seq.tostring(), seq_record.seq.tostring())
        self.assertEqual(records[0].description, seq_record.description)
        self.assertEqual(records[0].name, seq_record.name)
        self.assertEqual(records[0].id, seq_record.id)

        # With the RecordParser
        test_handle = open(datafile)
        records = list(SwissProt.parse(test_handle))
        test_handle.close()

        self.assertEqual(len(records), 1)
        self.assertTrue(isinstance(records[0], SwissProt.Record))

        #Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].sequence, record.sequence)
        self.assertEqual(records[0].description, record.description)
        self.assertEqual(records[0].entry_name, record.entry_name)
        self.assertEqual(records[0].accessions, record.accessions)

    def test_sp003(self):
        "Parsing SwissProt file sp003"

        filename = 'sp003'
        # test the record parser

        datafile = os.path.join('SwissProt', filename)

        test_handle = open(datafile)
        seq_record = SeqIO.read(test_handle, "swiss")
        test_handle.close()

        self.assertTrue(isinstance(seq_record, SeqRecord))

        self.assertEqual(seq_record.id, "P42655")
        self.assertEqual(seq_record.name, "143E_HUMAN")
        self.assertEqual(seq_record.description, "14-3-3 PROTEIN EPSILON (MITOCHONDRIAL IMPORT STIMULATION FACTOR L SUBUNIT) (PROTEIN KINASE C INHIBITOR PROTEIN-1) (KCIP-1) (14-3-3E).")
        self.assertEqual(repr(seq_record.seq), "Seq('MDDREDLVYQAKLAEQAERYDEMVESMKKVAGMDVELTVEERNLLSVAYKNVIG...ENQ', ProteinAlphabet())")

        test_handle = open(datafile)
        record = SwissProt.read(test_handle)
        test_handle.close()

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(record.entry_name, "143E_HUMAN")
        self.assertEqual(record.accessions, ['P42655', 'P29360', 'Q63631'])
        self.assertEqual(record.organism_classification, ['EUKARYOTA', 'METAZOA', 'CHORDATA', 'VERTEBRATA', 'MAMMALIA', 'EUTHERIA', 'PRIMATES', 'CATARRHINI', 'HOMINIDAE', 'HOMO'])
        self.assertEqual(record.seqinfo, (255, 29174, '40A43E62'))

        self.assertEqual(len(record.features), 5)
        self.assertEqual(record.features[0], ('MOD_RES', 1, 1, 'ACETYLATION.', ''))
        self.assertEqual(record.features[1], ('CONFLICT', 73, 73, 'K -> T (IN REF. 8).', ''))
        self.assertEqual(record.features[2], ('CONFLICT', 120, 120, 'F -> S (IN REF. 8).', ''))
        self.assertEqual(record.features[3], ('CONFLICT', 123, 123, 'K -> Y (IN REF. 8).', ''))
        self.assertEqual(record.features[4], ('CONFLICT', 129, 129, 'H -> Y (IN REF. 13).', ''))

        self.assertEqual(len(record.references), 13)
        self.assertEqual(record.references[0].authors, "CONKLIN D.S., GALAKTIONOV K., BEACH D.")
        self.assertEqual(record.references[0].title, "14-3-3 proteins associate with cdc25 phosphatases.")
        self.assertEqual(len(record.references[0].references), 1)
        self.assertEqual(record.references[0].references[0], ('MEDLINE', '95372385'))
        self.assertEqual(record.references[1].authors, "LUK S.C.W., LEE C.Y., WAYE M.M.Y.")
        self.assertEqual(record.references[1].title, "")
        self.assertEqual(len(record.references[1].references), 0)
        self.assertEqual(record.references[2].authors, "JIN D.Y., LYU M.S., KOZAK C.A., JEANG K.T.")
        self.assertEqual(record.references[2].title, "Function of 14-3-3 proteins.")
        self.assertEqual(len(record.references[2].references), 1)
        self.assertEqual(record.references[2].references[0], ('MEDLINE', '96300316'))
        self.assertEqual(record.references[3].authors, "CHONG S.S., TANIGAMI A., ROSCHKE A.V., LEDBETTER D.H.")
        self.assertEqual(record.references[3].title, "14-3-3 epsilon has no homology to LIS1 and lies telomeric to it on chromosome 17p13.3 outside the Miller-Dieker syndrome chromosome region.")
        self.assertEqual(len(record.references[3].references), 1)
        self.assertEqual(record.references[3].references[0], ('MEDLINE', '97011338'))
        self.assertEqual(record.references[4].authors, "TANIGAMI A., CHONG S.S., LEDBETTER D.H.")
        self.assertEqual(record.references[4].title, "14-3-3 epsilon genomic sequence.")
        self.assertEqual(len(record.references[4].references), 0)
        self.assertEqual(record.references[5].authors, "ROSEBOOM P.H., WELLER J.L., BABILA T., AITKEN A., SELLERS L.A., MOFFET J.R., NAMBOODIRI M.A., KLEIN D.C.")
        self.assertEqual(record.references[5].title, "Cloning and characterization of the epsilon and zeta isoforms of the 14-3-3 proteins.")
        self.assertEqual(len(record.references[5].references), 1)
        self.assertEqual(record.references[5].references[0], ('MEDLINE', '94296566'))
        self.assertEqual(record.references[6].authors, "ALAM R., HACHIYA N., SAKAGUCHI M., SHUN-ICHIRO K., IWANAGA S., KITAJIMA M., MIHARA K., OMURA T.")
        self.assertEqual(record.references[6].title, "cDNA cloning and characterization of mitochondrial import stimulation factor (MSF) purified from rat liver cytosol.")
        self.assertEqual(len(record.references[6].references), 1)
        self.assertEqual(record.references[6].references[0], ('MEDLINE', '95122474'))
        self.assertEqual(record.references[7].authors, "GAO L., GU X.B., YU D.S., YU R.K., ZENG G.")
        self.assertEqual(record.references[7].title, "Association of a 14-3-3 protein with CMP-NeuAc:GM1 alpha 2,3- sialyltransferase.")
        self.assertEqual(len(record.references[7].references), 1)
        self.assertEqual(record.references[7].references[0], ('MEDLINE', '96280718'))
        self.assertEqual(record.references[8].authors, "MCCONNELL J.E., ARMSTRONG J.F., BARD J.B.")
        self.assertEqual(record.references[8].title, "The mouse 14-3-3 epsilon isoform, a kinase regulator whose expression pattern is modulated in mesenchyme and neuronal differentiation.")
        self.assertEqual(len(record.references[8].references), 1)
        self.assertEqual(record.references[8].references[0], ('MEDLINE', '95269876'))
        self.assertEqual(record.references[9].authors, "TAKIHARA Y., IRIE K., NOMURA M., MOTALEB M., MATSUMOTO K., SHIMADA K.")
        self.assertEqual(record.references[9].title, "")
        self.assertEqual(len(record.references[9].references), 0)
        self.assertEqual(record.references[10].authors, "JONES J.M., NIIKURA T., PINKE R.M., GUO W., MOLDAY L., LEYKAM J., MCCONNELL D.G.")
        self.assertEqual(record.references[10].title, "Expression of 14-3-3 proteins in bovine retinal photoreceptors.")
        self.assertEqual(len(record.references[10].references), 0)
        self.assertEqual(record.references[11].authors, "TOKER A., SELLERS L.A., AMESS B., PATEL Y., HARRIS A., AITKEN A.")
        self.assertEqual(record.references[11].title, "Multiple isoforms of a protein kinase C inhibitor (KCIP-1/14-3-3) from sheep brain. Amino acid sequence of phosphorylated forms.")
        self.assertEqual(len(record.references[11].references), 1)
        self.assertEqual(record.references[11].references[0], ('MEDLINE', '92283271'))
        self.assertEqual(record.references[12].authors, "TOKER A., ELLIS C.A., SELLERS L.A., AITKEN A.")
        self.assertEqual(record.references[12].title, "Protein kinase C inhibitor proteins. Purification from sheep brain and sequence similarity to lipocortins and 14-3-3 protein.")
        self.assertEqual(len(record.references[12].references), 1)
        self.assertEqual(record.references[12].references[0], ('MEDLINE', '90345949'))

        #Check the two parsers agree on the essentials
        self.assertEqual(seq_record.seq.tostring(), record.sequence)
        self.assertEqual(seq_record.description, record.description)
        self.assertEqual(seq_record.name, record.entry_name)
        self.assertTrue(seq_record.id in record.accessions)

        #Now try using the iterator - note that all these
        #test cases have only one record.

        # With the SequenceParser
        test_handle = open(datafile)
        records = list(SeqIO.parse(test_handle, "swiss"))
        test_handle.close()

        self.assertEqual(len(records), 1)
        self.assertTrue(isinstance(records[0], SeqRecord))

        #Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].seq.tostring(), seq_record.seq.tostring())
        self.assertEqual(records[0].description, seq_record.description)
        self.assertEqual(records[0].name, seq_record.name)
        self.assertEqual(records[0].id, seq_record.id)

        # With the RecordParser
        test_handle = open(datafile)
        records = list(SwissProt.parse(test_handle))
        test_handle.close()

        self.assertEqual(len(records), 1)
        self.assertTrue(isinstance(records[0], SwissProt.Record))

        #Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].sequence, record.sequence)
        self.assertEqual(records[0].description, record.description)
        self.assertEqual(records[0].entry_name, record.entry_name)
        self.assertEqual(records[0].accessions, record.accessions)


    def test_sp004(self):
        "Parsing SwissProt file sp004"

        filename = 'sp004'
        # test the record parser

        datafile = os.path.join('SwissProt', filename)

        test_handle = open(datafile)
        seq_record = SeqIO.read(test_handle, "swiss")
        test_handle.close()

        self.assertTrue(isinstance(seq_record, SeqRecord))

        self.assertEqual(seq_record.id, "P23082")
        self.assertEqual(seq_record.name, "NDOA_PSEPU")
        self.assertEqual(seq_record.description, "NAPHTHALENE 1,2-DIOXYGENASE SYSTEM FERREDOXIN COMPONENT.")
        self.assertEqual(repr(seq_record.seq), "Seq('TVKWIEAVALSDILEGDVLGVTVEGKELALYEVEGEIYATDNLCTHGSARMSDG...DLS', ProteinAlphabet())")

        test_handle = open(datafile)
        record = SwissProt.read(test_handle)
        test_handle.close()

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(record.entry_name, "NDOA_PSEPU")
        self.assertEqual(record.accessions, ['P23082', 'Q52123', 'O07829'])
        self.assertEqual(record.organism_classification, ['Bacteria', 'Proteobacteria', 'gamma subdivision', 'Pseudomonas group', 'Pseudomonas'])
        self.assertEqual(record.seqinfo, (103, 11315, '9F91B3C8'))
    
        self.assertEqual(len(record.features), 12)
        self.assertEqual(record.features[0], ('INIT_MET', 0, 0, '', ''))
        self.assertEqual(record.features[1], ('METAL', 44, 44, 'IRON-SULFUR (2FE-2S) (POTENTIAL).', ''))
        self.assertEqual(record.features[2], ('METAL', 46, 46, 'IRON-SULFUR (2FE-2S) (POTENTIAL).', ''))
        self.assertEqual(record.features[3], ('METAL', 63, 63, 'IRON-SULFUR (2FE-2S) (POTENTIAL).', ''))
        self.assertEqual(record.features[4], ('METAL', 66, 66, 'IRON-SULFUR (2FE-2S) (POTENTIAL).', ''))
        self.assertEqual(record.features[5], ('VARIANT', 2, 2, 'V -> E (IN STRAIN G7).', ''))
        self.assertEqual(record.features[6], ('VARIANT', 14, 14, 'L -> P (IN STRAIN G7).', ''))
        self.assertEqual(record.features[7], ('VARIANT', 48, 48, 'S -> A (IN STRAIN G7).', ''))
        self.assertEqual(record.features[8], ('VARIANT', 76, 76, 'K -> R (IN STRAIN G7).', ''))
        self.assertEqual(record.features[9], ('VARIANT', 84, 84, 'Q -> E (IN STRAIN G7).', ''))
        self.assertEqual(record.features[10], ('VARIANT', 90, 90, 'P -> A (IN STRAIN G7).', ''))
        self.assertEqual(record.features[11], ('VARIANT', 103, 103, 'S -> GEF (IN STRAIN G7).', ''))

        self.assertEqual(len(record.references), 4) 
        self.assertEqual(record.references[0].authors, "KURKELA S., LEHVAESLAIHO H., PALVA E.T., TEERI T.H.")
        self.assertEqual(record.references[0].title, "Cloning, nucleotide sequence and characterization of genes encoding naphthalene dioxygenase of Pseudomonas putida strain NCIB9816.")
        self.assertEqual(len(record.references[0].references), 1)
        self.assertEqual(record.references[0].references[0], ('MEDLINE', '89211973'))
        self.assertEqual(record.references[1].authors, "SIMON M.J., OSSLUND T.D., SAUNDERS R., ENSLEY B.D., SUGGS S., HARCOURT A.A., SUEN W.-C., CRUDEN D.L., GIBSON D.T., ZYLSTRA G.J.")
        self.assertEqual(record.references[1].title, "Sequences of genes encoding naphthalene dioxygenase in Pseudomonas putida strains G7 and NCIB 9816-4.")
        self.assertEqual(len(record.references[1].references), 1)
        self.assertEqual(record.references[1].references[0], ('MEDLINE', '93252277'))
        self.assertEqual(record.references[2].authors, "DENOME S.A., STANLEY D.C., OLSON E.S., YOUNG K.D.")
        self.assertEqual(record.references[2].title, "Metabolism of dibenzothiophene and naphthalene in Pseudomonas strains: complete DNA sequence of an upper naphthalene catabolic pathway.")
        self.assertEqual(len(record.references[2].references), 1)
        self.assertEqual(record.references[2].references[0], ('MEDLINE', '94042852'))
        self.assertEqual(record.references[3].authors, "HAMANN C.")
        self.assertEqual(record.references[3].title, "")
        self.assertEqual(len(record.references[3].references), 0)

        #Check the two parsers agree on the essentials
        self.assertEqual(seq_record.seq.tostring(), record.sequence)
        self.assertEqual(seq_record.description, record.description)
        self.assertEqual(seq_record.name, record.entry_name)
        self.assertTrue(seq_record.id in record.accessions)

        #Now try using the iterator - note that all these
        #test cases have only one record.

        # With the SequenceParser
        test_handle = open(datafile)
        records = list(SeqIO.parse(test_handle, "swiss"))
        test_handle.close()

        self.assertEqual(len(records), 1)
        self.assertTrue(isinstance(records[0], SeqRecord))

        #Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].seq.tostring(), seq_record.seq.tostring())
        self.assertEqual(records[0].description, seq_record.description)
        self.assertEqual(records[0].name, seq_record.name)
        self.assertEqual(records[0].id, seq_record.id)

        # With the RecordParser
        test_handle = open(datafile)
        records = list(SwissProt.parse(test_handle))
        test_handle.close()

        self.assertEqual(len(records), 1)
        self.assertTrue(isinstance(records[0], SwissProt.Record))

        #Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].sequence, record.sequence)
        self.assertEqual(records[0].description, record.description)
        self.assertEqual(records[0].entry_name, record.entry_name)
        self.assertEqual(records[0].accessions, record.accessions)

    def test_sp005(self):
        "Parsing SwissProt file sp005"

        filename = 'sp005'
        # test the record parser

        datafile = os.path.join('SwissProt', filename)

        test_handle = open(datafile)
        seq_record = SeqIO.read(test_handle, "swiss")
        test_handle.close()

        self.assertTrue(isinstance(seq_record, SeqRecord))

        self.assertEqual(seq_record.id, "P24973")
        self.assertEqual(seq_record.name, "NU3M_BALPH")
        self.assertEqual(seq_record.description, "NADH-UBIQUINONE OXIDOREDUCTASE CHAIN 3 (EC 1.6.5.3).")
        self.assertEqual(repr(seq_record.seq), "Seq('MNLLLTLLTNTTLALLLVFIAFWLPQLNVYAEKTSPYECGFDPMGSARLPFSMK...WAE', ProteinAlphabet())")

        test_handle = open(datafile)
        record = SwissProt.read(test_handle)
        test_handle.close()

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(record.entry_name, "NU3M_BALPH")
        self.assertEqual(record.accessions, ['P24973'])
        self.assertEqual(record.organism_classification, ['Eukaryota', 'Metazoa', 'Chordata', 'Craniata', 'Vertebrata', 'Mammalia', 'Eutheria', 'Cetartiodactyla', 'Cetacea', 'Mysticeti', 'Balaenopteridae', 'Balaenoptera'])
        self.assertEqual(record.seqinfo, (115, 13022, 'ACF02965'))
    
        self.assertEqual(len(record.features), 0)

        self.assertEqual(len(record.references), 2)
        self.assertEqual(record.references[0].authors, "ARNASON U., GULLBERG A., WIDEGREN B.")
        self.assertEqual(record.references[0].title, "The complete nucleotide sequence of the mitochondrial DNA of the fin whale, Balaenoptera physalus.")
        self.assertEqual(len(record.references[0].references), 1)
        self.assertEqual(record.references[0].references[0], ('MEDLINE', '92139449'))
        self.assertEqual(record.references[1].authors, "ARNASON U., GULLBERG A.")
        self.assertEqual(record.references[1].title, "Comparison between the complete mtDNA sequences of the blue and the fin whale, two species that can hybridize in nature.")
        self.assertEqual(len(record.references[1].references), 1)
        self.assertEqual(record.references[1].references[0], ('MEDLINE', '94141932'))

        #Check the two parsers agree on the essentials
        self.assertEqual(seq_record.seq.tostring(), record.sequence)
        self.assertEqual(seq_record.description, record.description)
        self.assertEqual(seq_record.name, record.entry_name)
        self.assertTrue(seq_record.id in record.accessions)

        #Now try using the iterator - note that all these
        #test cases have only one record.

        # With the SequenceParser
        test_handle = open(datafile)
        records = list(SeqIO.parse(test_handle, "swiss"))
        test_handle.close()

        self.assertEqual(len(records), 1)
        self.assertTrue(isinstance(records[0], SeqRecord))

        #Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].seq.tostring(), seq_record.seq.tostring())
        self.assertEqual(records[0].description, seq_record.description)
        self.assertEqual(records[0].name, seq_record.name)
        self.assertEqual(records[0].id, seq_record.id)

        # With the RecordParser
        test_handle = open(datafile)
        records = list(SwissProt.parse(test_handle))
        test_handle.close()

        self.assertEqual(len(records), 1)
        self.assertTrue(isinstance(records[0], SwissProt.Record))

        #Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].sequence, record.sequence)
        self.assertEqual(records[0].description, record.description)
        self.assertEqual(records[0].entry_name, record.entry_name)
        self.assertEqual(records[0].accessions, record.accessions)


    def test_sp006(self):
        "Parsing SwissProt file sp006"

        filename = 'sp006'
        # test the record parser

        datafile = os.path.join('SwissProt', filename)

        test_handle = open(datafile)
        seq_record = SeqIO.read(test_handle, "swiss")
        test_handle.close()

        self.assertTrue(isinstance(seq_record, SeqRecord))

        self.assertEqual(seq_record.id, "P39896")
        self.assertEqual(seq_record.name, "TCMO_STRGA")
        self.assertEqual(seq_record.description, "TETRACENOMYCIN POLYKETIDE SYNTHESIS 8-O-METHYL TRANSFERASE TCMO (EC 2.1.1.-).")
        self.assertEqual(repr(seq_record.seq), "Seq('MTPHTHVRGPGDILQLTMAFYGSRALISAVELDLFTLLAGKPLPLGELCERAGI...KPR', ProteinAlphabet())")

        test_handle = open(datafile)
        record = SwissProt.read(test_handle)
        test_handle.close()

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(record.entry_name, "TCMO_STRGA")
        self.assertEqual(record.accessions, ['P39896'])
        self.assertEqual(record.organism_classification, ['BACTERIA', 'FIRMICUTES', 'ACTINOBACTERIA', 'ACTINOBACTERIDAE', 'ACTINOMYCETALES', 'STREPTOMYCINEAE', 'STREPTOMYCETACEAE', 'STREPTOMYCES'])
        self.assertEqual(record.seqinfo, (339, 37035, '848B7337'))
    
        self.assertEqual(len(record.features), 0)

        self.assertEqual(len(record.references), 1)
        self.assertEqual(record.references[0].authors, "SUMMERS R.G., WENDT-PIENKOWSKI E., MOTAMEDI H., HUTCHINSON C.R.")
        self.assertEqual(record.references[0].title, "Nucleotide sequence of the tcmII-tcmIV region of the tetracenomycin C biosynthetic gene cluster of Streptomyces glaucescens and evidence that the tcmN gene encodes a multifunctional cyclase-dehydratase-O-methyl transferase.")
        self.assertEqual(len(record.references[0].references), 1)
        self.assertEqual(record.references[0].references[0], ('MEDLINE', '92193265'))

        #Check the two parsers agree on the essentials
        self.assertEqual(seq_record.seq.tostring(), record.sequence)
        self.assertEqual(seq_record.description, record.description)
        self.assertEqual(seq_record.name, record.entry_name)
        self.assertTrue(seq_record.id in record.accessions)

        #Now try using the iterator - note that all these
        #test cases have only one record.

        # With the SequenceParser
        test_handle = open(datafile)
        records = list(SeqIO.parse(test_handle, "swiss"))
        test_handle.close()

        self.assertEqual(len(records), 1)
        self.assertTrue(isinstance(records[0], SeqRecord))

        #Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].seq.tostring(), seq_record.seq.tostring())
        self.assertEqual(records[0].description, seq_record.description)
        self.assertEqual(records[0].name, seq_record.name)
        self.assertEqual(records[0].id, seq_record.id)

        # With the RecordParser
        test_handle = open(datafile)
        records = list(SwissProt.parse(test_handle))
        test_handle.close()

        self.assertEqual(len(records), 1)
        self.assertTrue(isinstance(records[0], SwissProt.Record))

        #Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].sequence, record.sequence)
        self.assertEqual(records[0].description, record.description)
        self.assertEqual(records[0].entry_name, record.entry_name)
        self.assertEqual(records[0].accessions, record.accessions)


    def test_sp007(self):
        "Parsing SwissProt file sp007"

        filename = 'sp007'
        # test the record parser

        datafile = os.path.join('SwissProt', filename)

        test_handle = open(datafile)
        seq_record = SeqIO.read(test_handle, "swiss")
        test_handle.close()

        self.assertTrue(isinstance(seq_record, SeqRecord))

        self.assertEqual(seq_record.id, "O95832")
        self.assertEqual(seq_record.name, "CLD1_HUMAN")
        self.assertEqual(seq_record.description, "CLAUDIN-1 (SENESCENCE-ASSOCIATED EPITHELIAL MEMBRANE PROTEIN).")
        self.assertEqual(repr(seq_record.seq), "Seq('MANAGLQLLGFILAFLGWIGAIVSTALPQWRIYSYAGDNIVTAQAMYEGLWMSC...DYV', ProteinAlphabet())")

        test_handle = open(datafile)
        record = SwissProt.read(test_handle)
        test_handle.close()

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(record.entry_name, "CLD1_HUMAN")
        self.assertEqual(record.accessions, ['O95832'])
        self.assertEqual(record.organism_classification, ['Eukaryota', 'Metazoa', 'Chordata', 'Craniata', 'Vertebrata', 'Mammalia', 'Eutheria', 'Primates', 'Catarrhini', 'Hominidae', 'Homo'])
        self.assertEqual(record.seqinfo, (211, 22744, '07269000E6C214F0'))
    
        self.assertEqual(len(record.features), 6)
        self.assertEqual(record.features[0], ('TRANSMEM', 8, 28, 'POTENTIAL.', ''))
        self.assertEqual(record.features[1], ('TRANSMEM', 82, 102, 'POTENTIAL.', ''))
        self.assertEqual(record.features[2], ('TRANSMEM', 116, 136, 'POTENTIAL.', ''))
        self.assertEqual(record.features[3], ('TRANSMEM', 164, 184, 'POTENTIAL.', ''))
        self.assertEqual(record.features[4], ('CONFLICT', 62, 62, 'I -> V (IN REF. 2).', ''))
        self.assertEqual(record.features[5], ('CONFLICT', 135, 135, 'V -> A (IN REF. 2).', ''))

        self.assertEqual(len(record.references), 2)
        self.assertEqual(record.references[0].authors, "Swisshelm K.L., Machl A., Planitzer S., Robertson R., Kubbies M., Hosier S.")
        self.assertEqual(record.references[0].title, "SEMP1, a senescence-associated cDNA isolated from human mammary epithelial cells, is a member of an epithelial membrane protein superfamily.")
        self.assertEqual(len(record.references[0].references), 1)
        self.assertEqual(record.references[0].references[0], ('MEDLINE', '99132301'))
        self.assertEqual(record.references[1].authors, "Mitic L.M., Anderson J.M.")
        self.assertEqual(record.references[1].title, "Human claudin-1 isolated from Caco-2 mRNA.")
        self.assertEqual(len(record.references[1].references), 0)

        #Check the two parsers agree on the essentials
        self.assertEqual(seq_record.seq.tostring(), record.sequence)
        self.assertEqual(seq_record.description, record.description)
        self.assertEqual(seq_record.name, record.entry_name)
        self.assertTrue(seq_record.id in record.accessions)

        #Now try using the iterator - note that all these
        #test cases have only one record.

        # With the SequenceParser
        test_handle = open(datafile)
        records = list(SeqIO.parse(test_handle, "swiss"))
        test_handle.close()

        self.assertEqual(len(records), 1)
        self.assertTrue(isinstance(records[0], SeqRecord))

        #Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].seq.tostring(), seq_record.seq.tostring())
        self.assertEqual(records[0].description, seq_record.description)
        self.assertEqual(records[0].name, seq_record.name)
        self.assertEqual(records[0].id, seq_record.id)

        # With the RecordParser
        test_handle = open(datafile)
        records = list(SwissProt.parse(test_handle))
        test_handle.close()

        self.assertEqual(len(records), 1)
        self.assertTrue(isinstance(records[0], SwissProt.Record))

        #Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].sequence, record.sequence)
        self.assertEqual(records[0].description, record.description)
        self.assertEqual(records[0].entry_name, record.entry_name)
        self.assertEqual(records[0].accessions, record.accessions)


    def test_sp008(self):
        "Parsing SwissProt file sp008"

        filename = 'sp008'
        # test the record parser

        datafile = os.path.join('SwissProt', filename)

        test_handle = open(datafile)
        seq_record = SeqIO.read(test_handle, "swiss")
        test_handle.close()

        self.assertTrue(isinstance(seq_record, SeqRecord))

        self.assertEqual(seq_record.id, "P01892")
        self.assertEqual(seq_record.name, "1A02_HUMAN")
        self.assertEqual(seq_record.description, "HLA CLASS I HISTOCOMPATIBILITY ANTIGEN, A-2 ALPHA CHAIN PRECURSOR.")
        self.assertEqual(repr(seq_record.seq), "Seq('MAVMAPRTLVLLLSGALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGYVDD...CKV', ProteinAlphabet())")

        test_handle = open(datafile)
        record = SwissProt.read(test_handle)
        test_handle.close()

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(record.entry_name, "1A02_HUMAN")
        self.assertEqual(record.accessions, ['P01892', 'P06338', 'P30514', 'P30444', 'P30445', 'P30446', 'Q29680', 'Q29899', 'Q95352', 'Q29837', 'Q95380'])
        self.assertEqual(record.organism_classification, ['Eukaryota', 'Metazoa', 'Chordata', 'Craniata', 'Vertebrata', 'Mammalia', 'Eutheria', 'Primates', 'Catarrhini', 'Hominidae', 'Homo'])
        self.assertEqual(record.seqinfo, (365, 40922, 'B54A97B24B337C08'))
    
        self.assertEqual(len(record.features), 71)
        self.assertEqual(record.features[0], ('SIGNAL', 1, 24, '', ''))
        self.assertEqual(record.features[1], ('CHAIN', 25, 365, 'HLA CLASS I HISTOCOMPATIBILITY ANTIGEN, A-2 ALPHA CHAIN.', ''))
        self.assertEqual(record.features[2], ('DOMAIN', 25, 114, 'EXTRACELLULAR ALPHA-1.', ''))
        self.assertEqual(record.features[3], ('DOMAIN', 115, 206, 'EXTRACELLULAR ALPHA-2.', ''))
        self.assertEqual(record.features[4], ('DOMAIN', 207, 298, 'EXTRACELLULAR ALPHA-3.', ''))
        self.assertEqual(record.features[5], ('DOMAIN', 299, 308, 'CONNECTING PEPTIDE.', ''))
        self.assertEqual(record.features[6], ('TRANSMEM', 309, 332, '', ''))
        self.assertEqual(record.features[7], ('DOMAIN', 333, 365, 'CYTOPLASMIC TAIL.', ''))
        self.assertEqual(record.features[8], ('CARBOHYD', 110, 110, '', ''))
        self.assertEqual(record.features[9], ('DISULFID', 125, 188, '', ''))
        self.assertEqual(record.features[10], ('DISULFID', 227, 283, '', ''))
        self.assertEqual(record.features[11], ('STRAND', 27, 36, '', ''))
        self.assertEqual(record.features[12], ('STRAND', 45, 52, '', ''))
        self.assertEqual(record.features[13], ('TURN', 53, 54, '', ''))
        self.assertEqual(record.features[14], ('STRAND', 55, 61, '', ''))
        self.assertEqual(record.features[15], ('TURN', 62, 63, '', ''))
        self.assertEqual(record.features[16], ('STRAND', 70, 71, '', ''))
        self.assertEqual(record.features[17], ('HELIX', 74, 76, '', ''))
        self.assertEqual(record.features[18], ('TURN', 77, 78, '', ''))
        self.assertEqual(record.features[19], ('HELIX', 81, 108, '', ''))
        self.assertEqual(record.features[20], ('TURN', 109, 110, '', ''))
        self.assertEqual(record.features[21], ('TURN', 113, 114, '', ''))
        self.assertEqual(record.features[22], ('STRAND', 118, 127, '', ''))
        self.assertEqual(record.features[23], ('TURN', 129, 130, '', ''))
        self.assertEqual(record.features[24], ('STRAND', 133, 142, '', ''))
        self.assertEqual(record.features[25], ('TURN', 143, 144, '', ''))
        self.assertEqual(record.features[26], ('STRAND', 145, 150, '', ''))
        self.assertEqual(record.features[27], ('TURN', 152, 153, '', ''))
        self.assertEqual(record.features[28], ('STRAND', 157, 159, '', ''))
        self.assertEqual(record.features[29], ('TURN', 163, 163, '', ''))
        self.assertEqual(record.features[30], ('HELIX', 164, 173, '', ''))
        self.assertEqual(record.features[31], ('TURN', 174, 175, '', ''))
        self.assertEqual(record.features[32], ('HELIX', 176, 185, '', ''))
        self.assertEqual(record.features[33], ('TURN', 186, 186, '', ''))
        self.assertEqual(record.features[34], ('HELIX', 187, 198, '', ''))
        self.assertEqual(record.features[35], ('TURN', 199, 199, '', ''))
        self.assertEqual(record.features[36], ('HELIX', 200, 203, '', ''))
        self.assertEqual(record.features[37], ('TURN', 204, 204, '', ''))
        self.assertEqual(record.features[38], ('STRAND', 207, 207, '', ''))
        self.assertEqual(record.features[39], ('STRAND', 210, 219, '', ''))
        self.assertEqual(record.features[40], ('TURN', 220, 221, '', ''))
        self.assertEqual(record.features[41], ('STRAND', 222, 233, '', ''))
        self.assertEqual(record.features[42], ('STRAND', 238, 243, '', ''))
        self.assertEqual(record.features[43], ('TURN', 244, 245, '', ''))
        self.assertEqual(record.features[44], ('STRAND', 246, 247, '', ''))
        self.assertEqual(record.features[45], ('HELIX', 249, 251, '', ''))
        self.assertEqual(record.features[46], ('STRAND', 253, 254, '', ''))
        self.assertEqual(record.features[47], ('STRAND', 258, 259, '', ''))
        self.assertEqual(record.features[48], ('STRAND', 265, 274, '', ''))
        self.assertEqual(record.features[49], ('TURN', 275, 276, '', ''))
        self.assertEqual(record.features[50], ('HELIX', 278, 280, '', ''))
        self.assertEqual(record.features[51], ('STRAND', 281, 286, '', ''))
        self.assertEqual(record.features[52], ('TURN', 288, 289, '', ''))
        self.assertEqual(record.features[53], ('STRAND', 294, 297, '', ''))
        self.assertEqual(record.features[54], ('VARIANT', 33, 33, 'F -> Y (IN A*0205, A*0206, A*0208, A*0210 AND A*0221).', 'VAR_004334'))
        self.assertEqual(record.features[55], ('VARIANT', 54, 54, 'D -> N (IN A*0221).', 'VAR_004335'))
        self.assertEqual(record.features[56], ('VARIANT', 67, 67, 'Q -> R (IN A*0202, A*0205, AND A*0208).', 'VAR_004336'))
        self.assertEqual(record.features[57], ('VARIANT', 90, 90, 'K -> N (IN A*0208 AND A*0220).', 'VAR_004337'))
        self.assertEqual(record.features[58], ('VARIANT', 97, 98, 'TH -> ID (IN A*0211).', 'VAR_004338'))
        self.assertEqual(record.features[59], ('VARIANT', 119, 119, 'V -> L (IN A*0202, A*0205, A*0208 AND A*0217).', 'VAR_004339'))
        self.assertEqual(record.features[60], ('VARIANT', 121, 121, 'R -> M (IN A*0204 AND A*0217).', 'VAR_004340'))
        self.assertEqual(record.features[61], ('VARIANT', 123, 123, 'Y -> C (IN A*0207 AND A*0218).', 'VAR_004341'))
        self.assertEqual(record.features[62], ('VARIANT', 123, 123, 'Y -> F (IN A*0210 AND A*0217).', 'VAR_004342'))
        self.assertEqual(record.features[63], ('VARIANT', 131, 131, 'W -> G (IN A*0210).', 'VAR_004343'))
        self.assertEqual(record.features[64], ('VARIANT', 162, 162, 'M -> K (IN A*0218).', 'VAR_004344'))
        self.assertEqual(record.features[65], ('VARIANT', 173, 173, 'A -> T (IN A*0203).', 'VAR_004345'))
        self.assertEqual(record.features[66], ('VARIANT', 176, 176, 'V -> E (IN A*0203 AND A*0213).', 'VAR_004346'))
        self.assertEqual(record.features[67], ('VARIANT', 180, 180, 'L -> W (IN A*0202, A*0203, A*0205 AND A*0208).', 'VAR_004347'))
        self.assertEqual(record.features[68], ('VARIANT', 180, 180, 'L -> Q (IN A*0212 AND A*0213).', 'VAR_004348'))
        self.assertEqual(record.features[69], ('VARIANT', 187, 187, 'T -> E (IN A*0216).', 'VAR_004349'))
        self.assertEqual(record.features[70], ('VARIANT', 260, 260, 'A -> E (IN A*0209).', 'VAR_004350'))

        self.assertEqual(len(record.references), 27)
        self.assertEqual(record.references[0].authors, "Koller B.H., Orr H.T.")
        self.assertEqual(record.references[0].title, "Cloning and complete sequence of an HLA-A2 gene: analysis of two HLA-A alleles at the nucleotide level.")
        self.assertEqual(len(record.references[0].references), 1)
        self.assertEqual(record.references[0].references[0], ('MEDLINE', '85132727'))
        self.assertEqual(record.references[1].authors, "Cianetti L., Testa U., Scotto L., la Valle R., Simeone A., Boccoli G., Giannella G., Peschle C., Boncinelli E.")
        self.assertEqual(record.references[1].title, "Three new class I HLA alleles: structure of mRNAs and alternative mechanisms of processing.")
        self.assertEqual(len(record.references[1].references), 1)
        self.assertEqual(record.references[1].references[0], ('MEDLINE', '89122144'))
        self.assertEqual(record.references[2].authors, "Ennis P.D., Zemmour J., Salter R.D., Parham P.")
        self.assertEqual(record.references[2].title, "Rapid cloning of HLA-A,B cDNA by using the polymerase chain reaction: frequency and nature of errors produced in amplification.")
        self.assertEqual(len(record.references[2].references), 1)
        self.assertEqual(record.references[2].references[0], ('MEDLINE', '90207291'))
        self.assertEqual(record.references[3].authors, "Belich M.P., Madrigal J.A., Hildebrand W.H., Zemmour J., Williams R.C., Luz R., Petzl-Erler M.L., Parham P.")
        self.assertEqual(record.references[3].title, "Unusual HLA-B alleles in two tribes of Brazilian Indians.")
        self.assertEqual(len(record.references[3].references), 1)
        self.assertEqual(record.references[3].references[0], ('MEDLINE', '92269955'))
        self.assertEqual(record.references[4].authors, "Krangel M.S.")
        self.assertEqual(record.references[4].title, "Unusual RNA splicing generates a secreted form of HLA-A2 in a mutagenized B lymphoblastoid cell line.")
        self.assertEqual(len(record.references[4].references), 1)
        self.assertEqual(record.references[4].references[0], ('MEDLINE', '85230571'))
        self.assertEqual(record.references[5].authors, "Orr H.T., Lopez de Castro J.A., Parham P., Ploegh H.L., Strominger J.L.")
        self.assertEqual(record.references[5].title, "Comparison of amino acid sequences of two human histocompatibility antigens, HLA-A2 and HLA-B7: location of putative alloantigenic sites.")
        self.assertEqual(len(record.references[5].references), 1)
        self.assertEqual(record.references[5].references[0], ('MEDLINE', '80056745'))
        self.assertEqual(record.references[6].authors, "Lopez de Castro J.A., Strominger J.L., Strong D.M., Orr H.T.")
        self.assertEqual(record.references[6].title, "Structure of crossreactive human histocompatibility antigens HLA-A28 and HLA-A2: possible implications for the generation of HLA polymorphism.")
        self.assertEqual(len(record.references[6].references), 1)
        self.assertEqual(record.references[6].references[0], ('MEDLINE', '82247941'))
        self.assertEqual(record.references[7].authors, "Mattson D.H., Handy D.E., Bradley D.A., Coligan J.E., Cowan E.P., Biddison W.E.")
        self.assertEqual(record.references[7].title, "DNA sequences of the genes that encode the CTL-defined HLA-A2 variants M7 and DK1.")
        self.assertEqual(len(record.references[7].references), 1)
        self.assertEqual(record.references[7].references[0], ('MEDLINE', '87306734'))
        self.assertEqual(record.references[8].authors, "Holmes N., Ennis P., Wan A.M., Denney D.W., Parham P.")
        self.assertEqual(record.references[8].title, "Multiple genetic mechanisms have contributed to the generation of the HLA-A2/A28 family of class I MHC molecules.")
        self.assertEqual(len(record.references[8].references), 1)
        self.assertEqual(record.references[8].references[0], ('MEDLINE', '87252273'))
        self.assertEqual(record.references[9].authors, "Domena J.D.")
        self.assertEqual(record.references[9].title, "")
        self.assertEqual(len(record.references[9].references), 0)
        self.assertEqual(record.references[10].authors, "Castano A.R., Lopez de Castro J.A.")
        self.assertEqual(record.references[10].title, "Structure of the HLA-A*0204 antigen, found in South American Indians. Spatial clustering of HLA-A2 subtype polymorphism.")
        self.assertEqual(len(record.references[10].references), 1)
        self.assertEqual(record.references[10].references[0], ('MEDLINE', '92039809'))
        self.assertEqual(record.references[11].authors, "Watkins D.I., McAdam S.N., Liu X., Stang C.R., Milford E.L., Levine C.G., Garber T.L., Dogon A.L., Lord C.I., Ghim S.H., Troup G.M., Hughes A.L., Letvin N.L.")
        self.assertEqual(record.references[11].title, "New recombinant HLA-B alleles in a tribe of South American Amerindians indicate rapid evolution of MHC class I loci.")
        self.assertEqual(len(record.references[11].references), 1)
        self.assertEqual(record.references[11].references[0], ('MEDLINE', '92269956'))
        self.assertEqual(record.references[12].authors, "Parham P., Lawlor D.A., Lomen C.E., Ennis P.D.")
        self.assertEqual(record.references[12].title, "Diversity and diversification of HLA-A,B,C alleles.")
        self.assertEqual(len(record.references[12].references), 1)
        self.assertEqual(record.references[12].references[0], ('MEDLINE', '89235215'))
        self.assertEqual(record.references[13].authors, "Ezquerra A., Domenech N., van der Poel J., Strominger J.L., Vega M.A., Lopez de Castro J.A.")
        self.assertEqual(record.references[13].title, "Molecular analysis of an HLA-A2 functional variant CLA defined by cytolytic T lymphocytes.")
        self.assertEqual(len(record.references[13].references), 1)
        self.assertEqual(record.references[13].references[0], ('MEDLINE', '86305811'))
        self.assertEqual(record.references[14].authors, "Domenech N., Ezquerra A., Castano R., Lopez de Castro J.A.")
        self.assertEqual(record.references[14].title, "Structural analysis of HLA-A2.4 functional variant KNE. Implications for the mapping of HLA-A2-specific T-cell epitopes.")
        self.assertEqual(len(record.references[14].references), 1)
        self.assertEqual(record.references[14].references[0], ('MEDLINE', '88113844'))
        self.assertEqual(record.references[15].authors, "Domenech N., Castano R., Goulmy E., Lopez de Castro J.A.")
        self.assertEqual(record.references[15].title, "Molecular analysis of HLA-A2.4 functional variant KLO: close structural and evolutionary relatedness to the HLA-A2.2 subtype.")
        self.assertEqual(len(record.references[15].references), 1)
        self.assertEqual(record.references[15].references[0], ('MEDLINE', '88314183'))
        self.assertEqual(record.references[16].authors, "Castano R., Ezquerra A., Domenech N., Lopez de Castro J.A.")
        self.assertEqual(record.references[16].title, "An HLA-A2 population variant with structural polymorphism in the alpha 3 region.")
        self.assertEqual(len(record.references[16].references), 1)
        self.assertEqual(record.references[16].references[0], ('MEDLINE', '88186100'))
        self.assertEqual(record.references[17].authors, "Epstein H., Kennedy L., Holmes N.")
        self.assertEqual(record.references[17].title, "An Oriental HLA-A2 subtype is closely related to a subset of Caucasoid HLA-A2 alleles.")
        self.assertEqual(len(record.references[17].references), 1)
        self.assertEqual(record.references[17].references[0], ('MEDLINE', '89122133'))
        self.assertEqual(record.references[18].authors, "Castano A.R., Lopez de Castro J.A.")
        self.assertEqual(record.references[18].title, "Structure of the HLA-A*0211 (A2.5) subtype: further evidence for selection-driven diversification of HLA-A2 antigens.")
        self.assertEqual(len(record.references[18].references), 1)
        self.assertEqual(record.references[18].references[0], ('MEDLINE', '92218010'))
        self.assertEqual(record.references[19].authors, "Barber D.F., Fernandez J.M., Lopez de Castro J.A.")
        self.assertEqual(record.references[19].title, "Primary structure of a new HLA-A2 subtype: HLA-A*0213.")
        self.assertEqual(len(record.references[19].references), 1)
        self.assertEqual(record.references[19].references[0], ('MEDLINE', '94222455'))
        self.assertEqual(record.references[20].authors, "Barouch D., Krausa P., Bodmer J., Browning M.J., McMichael A.J.")
        self.assertEqual(record.references[20].title, "Identification of a novel HLA-A2 subtype, HLA-A*0216.")
        self.assertEqual(len(record.references[20].references), 1)
        self.assertEqual(record.references[20].references[0], ('MEDLINE', '95278976'))
        self.assertEqual(record.references[21].authors, "Selvakumar A., Granja C.B., Salazar M., Alosco S.M., Yunis E.J., Dupont B.")
        self.assertEqual(record.references[21].title, "A novel subtype of A2 (A*0217) isolated from the South American Indian B-cell line AMALA.")
        self.assertEqual(len(record.references[21].references), 1)
        self.assertEqual(record.references[21].references[0], ('MEDLINE', '95381236'))
        self.assertEqual(record.references[22].authors, "Kashiwase K., Tokunaga K., Ishikawa Y., Oohashi H., Hashimoto M., Akaza T., Tadokoro K., Juji T.")
        self.assertEqual(record.references[22].title, "A new A2 sequence HLA-A2K from Japanese.")
        self.assertEqual(len(record.references[22].references), 0)
        self.assertEqual(record.references[23].authors, "Fleischhauer K., Zino E., Mazzi B., Severini G.M., Benazzi E., Bordignon C.")
        self.assertEqual(record.references[23].title, "HLA-A*02 subtype distribution in Caucasians from northern Italy: identification of A*0220.")
        self.assertEqual(len(record.references[23].references), 1)
        self.assertEqual(record.references[23].references[0], ('MEDLINE', '97161038'))
        self.assertEqual(record.references[24].authors, "Szmania S., Baxter-Lowe L.A.")
        self.assertEqual(record.references[24].title, "Nucleotide sequence of a novel HLA-A2 gene.")
        self.assertEqual(len(record.references[24].references), 0)
        self.assertEqual(record.references[25].authors, "Bjorkman P.J., Saper M.A., Samraoui B., Bennett W.S., Strominger J.L., Wiley D.C.")
        self.assertEqual(record.references[25].title, "Structure of the human class I histocompatibility antigen, HLA-A2.")
        self.assertEqual(len(record.references[25].references), 1)
        self.assertEqual(record.references[25].references[0], ('MEDLINE', '88014204'))
        self.assertEqual(record.references[26].authors, "Saper M.A., Bjorkman P.J., Wiley D.C.")
        self.assertEqual(record.references[26].title, "Refined structure of the human histocompatibility antigen HLA-A2 at 2.6-A resolution.")
        self.assertEqual(len(record.references[26].references), 1)
        self.assertEqual(record.references[26].references[0], ('MEDLINE', '91245570'))

        #Check the two parsers agree on the essentials
        self.assertEqual(seq_record.seq.tostring(), record.sequence)
        self.assertEqual(seq_record.description, record.description)
        self.assertEqual(seq_record.name, record.entry_name)
        self.assertTrue(seq_record.id in record.accessions)

        #Now try using the iterator - note that all these
        #test cases have only one record.

        # With the SequenceParser
        test_handle = open(datafile)
        records = list(SeqIO.parse(test_handle, "swiss"))
        test_handle.close()

        self.assertEqual(len(records), 1)
        self.assertTrue(isinstance(records[0], SeqRecord))

        #Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].seq.tostring(), seq_record.seq.tostring())
        self.assertEqual(records[0].description, seq_record.description)
        self.assertEqual(records[0].name, seq_record.name)
        self.assertEqual(records[0].id, seq_record.id)

        # With the RecordParser
        test_handle = open(datafile)
        records = list(SwissProt.parse(test_handle))
        test_handle.close()

        self.assertEqual(len(records), 1)
        self.assertTrue(isinstance(records[0], SwissProt.Record))

        #Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].sequence, record.sequence)
        self.assertEqual(records[0].description, record.description)
        self.assertEqual(records[0].entry_name, record.entry_name)
        self.assertEqual(records[0].accessions, record.accessions)

    def test_sp009(self):
        "Parsing SwissProt file sp009"

        filename = 'sp009'
        # test the record parser

        datafile = os.path.join('SwissProt', filename)

        test_handle = open(datafile)
        seq_record = SeqIO.read(test_handle, "swiss")
        test_handle.close()

        self.assertTrue(isinstance(seq_record, SeqRecord))

        self.assertEqual(seq_record.id, "O23729")
        self.assertEqual(seq_record.name, "CHS3_BROFI")
        self.assertEqual(seq_record.description, "CHALCONE SYNTHASE 3 (EC 2.3.1.74) (NARINGENIN-CHALCONE SYNTHASE 3).")
        self.assertEqual(repr(seq_record.seq), "Seq('MAPAMEEIRQAQRAEGPAAVLAIGTSTPPNALYQADYPDYYFRITKSEHLTELK...GAE', ProteinAlphabet())")

        test_handle = open(datafile)
        record = SwissProt.read(test_handle)
        test_handle.close()

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(record.entry_name, "CHS3_BROFI")
        self.assertEqual(record.accessions, ['O23729'])
        self.assertEqual(record.organism_classification, ['Eukaryota', 'Viridiplantae', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Bromheadia'])
        self.assertEqual(record.seqinfo, (394, 42941, '2F8D14AF4870BBB2'))
    
        self.assertEqual(len(record.features), 1)
        self.assertEqual(record.features[0], ('ACT_SITE', 165, 165, 'BY SIMILARITY.', ''))

        self.assertEqual(len(record.references), 1)
        self.assertEqual(record.references[0].authors, "Liew C.F., Lim S.H., Loh C.S., Goh C.J.")
        self.assertEqual(record.references[0].title, "Molecular cloning and sequence analysis of chalcone synthase cDNAs of Bromheadia finlaysoniana.")
        self.assertEqual(len(record.references[0].references), 0)


        #Check the two parsers agree on the essentials
        self.assertEqual(seq_record.seq.tostring(), record.sequence)
        self.assertEqual(seq_record.description, record.description)
        self.assertEqual(seq_record.name, record.entry_name)
        self.assertTrue(seq_record.id in record.accessions)

        #Now try using the iterator - note that all these
        #test cases have only one record.

        # With the SequenceParser
        test_handle = open(datafile)
        records = list(SeqIO.parse(test_handle, "swiss"))
        test_handle.close()

        self.assertEqual(len(records), 1)
        self.assertTrue(isinstance(records[0], SeqRecord))

        #Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].seq.tostring(), seq_record.seq.tostring())
        self.assertEqual(records[0].description, seq_record.description)
        self.assertEqual(records[0].name, seq_record.name)
        self.assertEqual(records[0].id, seq_record.id)

        # With the RecordParser
        test_handle = open(datafile)
        records = list(SwissProt.parse(test_handle))
        test_handle.close()

        self.assertEqual(len(records), 1)
        self.assertTrue(isinstance(records[0], SwissProt.Record))

        #Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].sequence, record.sequence)
        self.assertEqual(records[0].description, record.description)
        self.assertEqual(records[0].entry_name, record.entry_name)
        self.assertEqual(records[0].accessions, record.accessions)

    def test_sp010(self):
        "Parsing SwissProt file sp010"

        filename = 'sp010'
        # test the record parser

        datafile = os.path.join('SwissProt', filename)

        test_handle = open(datafile)
        seq_record = SeqIO.read(test_handle, "swiss")
        test_handle.close()

        self.assertTrue(isinstance(seq_record, SeqRecord))

        self.assertEqual(seq_record.id, "Q13639")
        self.assertEqual(seq_record.name, "5H4_HUMAN")
        self.assertEqual(seq_record.description, "5-HYDROXYTRYPTAMINE 4 RECEPTOR (5-HT-4) (SEROTONIN RECEPTOR) (5-HT4).")
        self.assertEqual(repr(seq_record.seq), "Seq('MDKLDANVSSEEGFGSVEKVVLLTFLSTVILMAILGNLLVMVAVCWDRQLRKIK...SDT', ProteinAlphabet())")

        test_handle = open(datafile)
        record = SwissProt.read(test_handle)
        test_handle.close()

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(record.entry_name, "5H4_HUMAN")
        self.assertEqual(record.accessions, ['Q13639', 'Q9UBM6', 'Q9UQR6', 'Q9UE22', 'Q9UE23', 'Q9UBT4', 'Q9NY73'])
        self.assertEqual(record.organism_classification, ['Eukaryota', 'Metazoa', 'Chordata', 'Craniata', 'Vertebrata', 'Euteleostomi', 'Mammalia', 'Eutheria', 'Primates', 'Catarrhini', 'Hominidae', 'Homo'])
        self.assertEqual(record.seqinfo, (388, 43761, '7FCFEC60E7BDF560'))
    
        self.assertEqual(len(record.features), 23)
        self.assertEqual(record.features[0], ('DOMAIN', 1, 19, 'EXTRACELLULAR (POTENTIAL).', ''))
        self.assertEqual(record.features[1], ('TRANSMEM', 20, 40, '1 (POTENTIAL).', ''))
        self.assertEqual(record.features[2], ('DOMAIN', 41, 58, 'CYTOPLASMIC (POTENTIAL).', ''))
        self.assertEqual(record.features[3], ('TRANSMEM', 59, 79, '2 (POTENTIAL).', ''))
        self.assertEqual(record.features[4], ('DOMAIN', 80, 93, 'EXTRACELLULAR (POTENTIAL).', ''))
        self.assertEqual(record.features[5], ('TRANSMEM', 94, 116, '3 (POTENTIAL).', ''))
        self.assertEqual(record.features[6], ('DOMAIN', 117, 137, 'CYTOPLASMIC (POTENTIAL).', ''))
        self.assertEqual(record.features[7], ('TRANSMEM', 138, 158, '4 (POTENTIAL).', ''))
        self.assertEqual(record.features[8], ('DOMAIN', 159, 192, 'EXTRACELLULAR (POTENTIAL).', ''))
        self.assertEqual(record.features[9], ('TRANSMEM', 193, 213, '5 (POTENTIAL).', ''))
        self.assertEqual(record.features[10], ('DOMAIN', 214, 260, 'CYTOPLASMIC (POTENTIAL).', ''))
        self.assertEqual(record.features[11], ('TRANSMEM', 261, 281, '6 (POTENTIAL).', ''))
        self.assertEqual(record.features[12], ('DOMAIN', 282, 294, 'EXTRACELLULAR (POTENTIAL).', ''))
        self.assertEqual(record.features[13], ('TRANSMEM', 295, 315, '7 (POTENTIAL).', ''))
        self.assertEqual(record.features[14], ('DOMAIN', 316, 388, 'CYTOPLASMIC (POTENTIAL).', ''))
        self.assertEqual(record.features[15], ('CARBOHYD', 7, 7, 'N-LINKED (GLCNAC...) (POTENTIAL).', ''))
        self.assertEqual(record.features[16], ('DISULFID', 93, 184, 'BY SIMILARITY.', ''))
        self.assertEqual(record.features[17], ('LIPID', 329, 329, 'PALMITATE (BY SIMILARITY).', ''))
        self.assertEqual(record.features[18], ('VARSPLIC', 169, 169, 'L -> LERSLNQGLGQDFHA (IN ISOFORM 5- HT4(F)).', ''))
        self.assertEqual(record.features[19], ('VARSPLIC', 359, 388, 'RDAVECGGQWESQCHPPATSPLVAAQPSDT -> SGCSPVSSFLLLFCNRPVPV (IN ISOFORM 5-HT4(E)).', ''))
        self.assertEqual(record.features[20], ('VARSPLIC', 359, 388, 'RDAVECGGQWESQCHPPATSPLVAAQPSDT -> SSGTETDRRNFGIRKRRLTKPS (IN ISOFORM 5-HT4(D)).', ''))
        self.assertEqual(record.features[21], ('VARSPLIC', 360, 388, 'DAVECGGQWESQCHPPATSPLVAAQPSDT -> F (IN ISOFORM 5-HT4(C)).', ''))
        self.assertEqual(record.features[22], ('VARSPLIC', 360, 388, 'DAVECGGQWESQCHPPATSPLVAAQPSDT -> YTVLHRGHHQELEKLPIHNDPESLESCF (IN ISOFORM 5- HT4(A)).', ''))
        self.assertEqual(len(record.references), 6)

        self.assertEqual(record.references[0].authors, "Blondel O., Gastineau M., Dahmoune Y., Langlois M., Fischmeister R.")
        self.assertEqual(record.references[0].title, "Cloning, expression, and pharmacology of four human 5- hydroxytryptamine receptor isoforms produced by alternative splicing in the carboxyl terminus.")
        self.assertEqual(len(record.references[0].references), 1)
        self.assertEqual(record.references[0].references[0], ('PubMed', '9603189'))
        self.assertEqual(record.references[1].authors, "Van den Wyngaert I., Gommeren W., Jurzak M., Verhasselt P., Gordon R., Leysen J., Luyten W., Bender E.")
        self.assertEqual(record.references[1].title, "Cloning and expression of 5-HT4 receptor species and splice variants.")
        self.assertEqual(len(record.references[1].references), 0)
        self.assertEqual(record.references[2].authors, "Claeysen S., Faye P., Sebben M., Lemaire S., Bockaert J., Dumuis A.")
        self.assertEqual(record.references[2].title, "Cloning and expression of human 5-HT4S receptors. Effect of receptor density on their coupling to adenylyl cyclase.")
        self.assertEqual(len(record.references[2].references), 1)
        self.assertEqual(record.references[2].references[0], ('PubMed', '9351641'))
        self.assertEqual(record.references[3].authors, "Claeysen S., Sebben M., Becamel C., Bockaert J., Dumuis A.")
        self.assertEqual(record.references[3].title, "Novel brain-specific 5-HT4 receptors splice variants show marked constitutive activity: role of the c-terminal intracellular domain.")
        self.assertEqual(len(record.references[3].references), 0)
        self.assertEqual(record.references[4].authors, "Bender E., Pindon A., van Oers I., Zhang Y.B., Gommeren W., Verhasselt P., Jurzak M., Leysen J., Luyten W.")
        self.assertEqual(record.references[4].title, "Structure of the human serotonin 5-HT4 receptor gene and cloning of a novel 5-HT4 splice variant.")
        self.assertEqual(len(record.references[4].references), 1)
        self.assertEqual(record.references[4].references[0], ('PubMed', '10646498'))
        self.assertEqual(record.references[5].authors, "Ullmer C., Schmuck K., Kalkman H.O., Lubbert H.")
        self.assertEqual(record.references[5].title, "Expression of serotonin receptor mRNAs in blood vessels.")
        self.assertEqual(len(record.references[5].references), 2)
        self.assertEqual(record.references[5].references[0], ('MEDLINE', '95385798'))
        self.assertEqual(record.references[5].references[1], ('PubMed', '7656980'))

        #Check the two parsers agree on the essentials
        self.assertEqual(seq_record.seq.tostring(), record.sequence)
        self.assertEqual(seq_record.description, record.description)
        self.assertEqual(seq_record.name, record.entry_name)
        self.assertTrue(seq_record.id in record.accessions)

        #Now try using the iterator - note that all these
        #test cases have only one record.

        # With the SequenceParser
        test_handle = open(datafile)
        records = list(SeqIO.parse(test_handle, "swiss"))
        test_handle.close()

        self.assertEqual(len(records), 1)
        self.assertTrue(isinstance(records[0], SeqRecord))

        #Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].seq.tostring(), seq_record.seq.tostring())
        self.assertEqual(records[0].description, seq_record.description)
        self.assertEqual(records[0].name, seq_record.name)
        self.assertEqual(records[0].id, seq_record.id)

        # With the RecordParser
        test_handle = open(datafile)
        records = list(SwissProt.parse(test_handle))
        test_handle.close()

        self.assertEqual(len(records), 1)
        self.assertTrue(isinstance(records[0], SwissProt.Record))

        #Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].sequence, record.sequence)
        self.assertEqual(records[0].description, record.description)
        self.assertEqual(records[0].entry_name, record.entry_name)
        self.assertEqual(records[0].accessions, record.accessions)

    def test_sp011(self):
        "Parsing SwissProt file sp011"

        filename = 'sp011'
        # test the record parser

        datafile = os.path.join('SwissProt', filename)

        test_handle = open(datafile)
        seq_record = SeqIO.read(test_handle, "swiss")
        test_handle.close()

        self.assertTrue(isinstance(seq_record, SeqRecord))

        self.assertEqual(seq_record.id, "P16235")
        self.assertEqual(seq_record.name, "LSHR_RAT")
        self.assertEqual(seq_record.description, "LUTROPIN-CHORIOGONADOTROPIC HORMONE RECEPTOR PRECURSOR (LH/CG-R) (LSH-R) (LUTEINIZING HORMONE RECEPTOR).")
        self.assertEqual(repr(seq_record.seq), "Seq('MGRRVPALRQLLVLAVLLLKPSQLQSRELSGSRCPEPCDCAPDGALRCPGPRAG...LTH', ProteinAlphabet())")

        test_handle = open(datafile)
        record = SwissProt.read(test_handle)
        test_handle.close()

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(record.entry_name, "LSHR_RAT")
        self.assertEqual(record.accessions, ['P16235', 'P70646', 'Q63807', 'Q63808', 'Q63809'])
        self.assertEqual(record.organism_classification, ['Eukaryota', 'Metazoa', 'Chordata', 'Craniata', 'Vertebrata', 'Euteleostomi', 'Mammalia', 'Eutheria', 'Rodentia', 'Sciurognathi', 'Muridae', 'Murinae', 'Rattus'])
        self.assertEqual(record.seqinfo, (700, 78035, '31807E73BAC94F1F'))
    
        self.assertEqual(len(record.features), 52)
        self.assertEqual(record.features[0], ('SIGNAL', 1, 26, '', ''))
        self.assertEqual(record.features[1], ('CHAIN', 27, 700, 'LUTROPIN-CHORIOGONADOTROPIC HORMONE RECEPTOR.', ''))
        self.assertEqual(record.features[2], ('DOMAIN', 27, 362, 'EXTRACELLULAR (POTENTIAL).', ''))
        self.assertEqual(record.features[3], ('TRANSMEM', 363, 390, '1 (POTENTIAL).', ''))
        self.assertEqual(record.features[4], ('DOMAIN', 391, 399, 'CYTOPLASMIC (POTENTIAL).', ''))
        self.assertEqual(record.features[5], ('TRANSMEM', 400, 422, '2 (POTENTIAL).', ''))
        self.assertEqual(record.features[6], ('DOMAIN', 423, 443, 'EXTRACELLULAR (POTENTIAL).', ''))
        self.assertEqual(record.features[7], ('TRANSMEM', 444, 466, '3 (POTENTIAL).', ''))
        self.assertEqual(record.features[8], ('DOMAIN', 467, 486, 'CYTOPLASMIC (POTENTIAL).', ''))
        self.assertEqual(record.features[9], ('TRANSMEM', 487, 509, '4 (POTENTIAL).', ''))
        self.assertEqual(record.features[10], ('DOMAIN', 510, 529, 'EXTRACELLULAR (POTENTIAL).', ''))
        self.assertEqual(record.features[11], ('TRANSMEM', 530, 551, '5 (POTENTIAL).', ''))
        self.assertEqual(record.features[12], ('DOMAIN', 552, 574, 'CYTOPLASMIC (POTENTIAL).', ''))
        self.assertEqual(record.features[13], ('TRANSMEM', 575, 598, '6 (POTENTIAL).', ''))
        self.assertEqual(record.features[14], ('DOMAIN', 599, 609, 'EXTRACELLULAR (POTENTIAL).', ''))
        self.assertEqual(record.features[15], ('TRANSMEM', 610, 631, '7 (POTENTIAL).', ''))
        self.assertEqual(record.features[16], ('DOMAIN', 632, 700, 'CYTOPLASMIC (POTENTIAL).', ''))
        self.assertEqual(record.features[17], ('REPEAT', 52, 75, 'LRR 1.', ''))
        self.assertEqual(record.features[18], ('REPEAT', 126, 150, 'LRR 2.', ''))
        self.assertEqual(record.features[19], ('REPEAT', 152, 175, 'LRR 3.', ''))
        self.assertEqual(record.features[20], ('REPEAT', 176, 200, 'LRR 4.', ''))
        self.assertEqual(record.features[21], ('REPEAT', 202, 224, 'LRR 5.', ''))
        self.assertEqual(record.features[22], ('REPEAT', 225, 248, 'LRR 6.', ''))
        self.assertEqual(record.features[23], ('REPEAT', 250, 271, 'LRR 7.', ''))
        self.assertEqual(record.features[24], ('DISULFID', 443, 518, 'BY SIMILARITY.', ''))
        self.assertEqual(record.features[25], ('CARBOHYD', 103, 103, 'N-LINKED (GLCNAC...) (POTENTIAL).', ''))
        self.assertEqual(record.features[26], ('CARBOHYD', 178, 178, 'N-LINKED (GLCNAC...) (POTENTIAL).', ''))
        self.assertEqual(record.features[27], ('CARBOHYD', 199, 199, 'N-LINKED (GLCNAC...) (POTENTIAL).', ''))
        self.assertEqual(record.features[28], ('CARBOHYD', 295, 295, 'N-LINKED (GLCNAC...) (POTENTIAL).', ''))
        self.assertEqual(record.features[29], ('CARBOHYD', 303, 303, 'N-LINKED (GLCNAC...) (POTENTIAL).', ''))
        self.assertEqual(record.features[30], ('CARBOHYD', 317, 317, 'N-LINKED (GLCNAC...) (POTENTIAL).', ''))
        self.assertEqual(record.features[31], ('VARSPLIC', 83, 132, 'MISSING (IN ISOFORM 1950).', ''))
        self.assertEqual(record.features[32], ('VARSPLIC', 133, 157, 'MISSING (IN ISOFORM 1759).', ''))
        self.assertEqual(record.features[33], ('VARSPLIC', 184, 700, 'MISSING (IN ISOFORM C2).', ''))
        self.assertEqual(record.features[34], ('VARSPLIC', 232, 251, 'DISSTKLQALPSHGLESIQT -> PCRATGWSPFRRSSPCLPTH (IN ISOFORM 2075).', ''))
        self.assertEqual(record.features[35], ('VARSPLIC', 232, 293, 'MISSING (IN ISOFORM E/A2, ISOFORM EB AND ISOFORM B1).', ''))
        self.assertEqual(record.features[36], ('VARSPLIC', 252, 700, 'MISSING (IN ISOFORM 2075).', ''))
        self.assertEqual(record.features[37], ('VARSPLIC', 294, 367, 'QNFSFSIFENFSKQCESTVRKADNETLYSAIFEENELSGWDYDYGFCSPKTLQCAPEPDAFNPCEDIMGYAFLR -> IFHFPFLKTSPNNAKAQLEKQITRRFIPPSLRRMNSVAGIMIMASVHPRHSNVLQNQMLSTPVKILWAMPSLGS (IN ISOFORM B1 AND ISOFORM B3).', ''))
        self.assertEqual(record.features[38], ('VARSPLIC', 294, 294, 'Q -> P (IN ISOFORM C1).', ''))
        self.assertEqual(record.features[39], ('VARSPLIC', 295, 700, 'MISSING (IN ISOFORM C1).', ''))
        self.assertEqual(record.features[40], ('VARSPLIC', 321, 342, 'YSAIFEENELSGWDYDYGFCSP -> LHGALPAAHCLRGLPNKRPVL (IN ISOFORM 1834, ISOFORM 1759 AND ISOFORM EB).', ''))
        self.assertEqual(record.features[41], ('VARSPLIC', 343, 700, 'MISSING (IN ISOFORMS 1834, ISOFORM 1759 AND ISOFORM EB).', ''))
        self.assertEqual(record.features[42], ('VARSPLIC', 368, 700, 'MISSING (IN ISOFORM B1 AND ISOFORM B3).', ''))
        self.assertEqual(record.features[43], ('VARIANT', 82, 82, 'I -> M (IN ISOFORM 1950).', ''))
        self.assertEqual(record.features[44], ('VARIANT', 179, 179, 'E -> G (IN ISOFORM 1759).', ''))
        self.assertEqual(record.features[45], ('VARIANT', 233, 233, 'I -> T (IN ISOFORM 1950).', ''))
        self.assertEqual(record.features[46], ('VARIANT', 646, 646, 'G -> S (IN ISOFORM 1950).', ''))
        self.assertEqual(record.features[47], ('MUTAGEN', 409, 409, 'D->N: SIGNIFICANT REDUCTION OF BINDING.', ''))
        self.assertEqual(record.features[48], ('MUTAGEN', 436, 436, 'D->N: NO CHANGE IN BINDING OR CAMP PROD.', ''))
        self.assertEqual(record.features[49], ('MUTAGEN', 455, 455, 'E->Q: NO CHANGE IN BINDING OR CAMP PROD.', ''))
        self.assertEqual(record.features[50], ('MUTAGEN', 582, 582, 'D->N: NO CHANGE IN BINDING OR CAMP PROD.', ''))
        self.assertEqual(record.features[51], ('CONFLICT', 33, 33, 'R -> L (IN REF. 7).', ''))

        self.assertEqual(len(record.references), 8)
        self.assertEqual(record.references[0].authors, "McFarland K.C., Sprengel R., Phillips H.S., Koehler M., Rosemblit N., Nikolics K., Segaloff D.L., Seeburg P.H.")
        self.assertEqual(record.references[0].title, "Lutropin-choriogonadotropin receptor: an unusual member of the G protein-coupled receptor family.")
        self.assertEqual(len(record.references[0].references), 2)
        self.assertEqual(record.references[0].references[0], ('MEDLINE', '89332512'))
        self.assertEqual(record.references[0].references[1], ('PubMed', '2502842'))
        self.assertEqual(record.references[1].authors, "Aatsinki J.T., Pietila E.M., Lakkakorpi J.T., Rajaniemi H.J.")
        self.assertEqual(record.references[1].title, "Expression of the LH/CG receptor gene in rat ovarian tissue is regulated by an extensive alternative splicing of the primary transcript.")
        self.assertEqual(len(record.references[1].references), 2)
        self.assertEqual(record.references[1].references[0], ('MEDLINE', '92347604'))
        self.assertEqual(record.references[1].references[1], ('PubMed', '1353463'))
        self.assertEqual(record.references[2].authors, "Koo Y.B., Slaughter R.G., Ji T.H.")
        self.assertEqual(record.references[2].title, "Structure of the luteinizing hormone receptor gene and multiple exons of the coding sequence.")
        self.assertEqual(len(record.references[2].references), 2)
        self.assertEqual(record.references[2].references[0], ('MEDLINE', '91209270'))
        self.assertEqual(record.references[2].references[1], ('PubMed', '2019252'))
        self.assertEqual(record.references[3].authors, "Bernard M.P., Myers R.V., Moyle W.R.")
        self.assertEqual(record.references[3].title, "Cloning of rat lutropin (LH) receptor analogs lacking the soybean lectin domain.")
        self.assertEqual(len(record.references[3].references), 2)
        self.assertEqual(record.references[3].references[0], ('MEDLINE', '91006819'))
        self.assertEqual(record.references[3].references[1], ('PubMed', '1976554'))
        self.assertEqual(record.references[4].authors, "Segaloff D.L., Sprengel R., Nikolics K., Ascoli M.")
        self.assertEqual(record.references[4].title, "Structure of the lutropin/choriogonadotropin receptor.")
        self.assertEqual(len(record.references[4].references), 2)
        self.assertEqual(record.references[4].references[0], ('MEDLINE', '91126285'))
        self.assertEqual(record.references[4].references[1], ('PubMed', '2281186'))
        self.assertEqual(record.references[5].authors, "Tsai-Morris C.H., Buczko E., Wang W., Dufau M.L.")
        self.assertEqual(record.references[5].title, "Intronic nature of the rat luteinizing hormone receptor gene defines a soluble receptor subspecies with hormone binding activity.")
        self.assertEqual(len(record.references[5].references), 2)
        self.assertEqual(record.references[5].references[0], ('MEDLINE', '91060531'))
        self.assertEqual(record.references[5].references[1], ('PubMed', '2174034'))
        self.assertEqual(record.references[6].authors, "Roche P.C., Ryan R.J.")
        self.assertEqual(record.references[6].title, "Purification, characterization, and amino-terminal sequence of rat ovarian receptor for luteinizing hormone/human choriogonadotropin.")
        self.assertEqual(len(record.references[6].references), 2)
        self.assertEqual(record.references[6].references[0], ('MEDLINE', '89174723'))
        self.assertEqual(record.references[6].references[1], ('PubMed', '2925659'))
        self.assertEqual(record.references[7].authors, "Ji I., Ji T.H.")
        self.assertEqual(record.references[7].title, "Asp383 in the second transmembrane domain of the lutropin receptor is important for high affinity hormone binding and cAMP production.")
        self.assertEqual(len(record.references[7].references), 2)
        self.assertEqual(record.references[7].references[0], ('MEDLINE', '91332007'))
        self.assertEqual(record.references[7].references[1], ('PubMed', '1714448'))

        #Check the two parsers agree on the essentials
        self.assertEqual(seq_record.seq.tostring(), record.sequence)
        self.assertEqual(seq_record.description, record.description)
        self.assertEqual(seq_record.name, record.entry_name)
        self.assertTrue(seq_record.id in record.accessions)

        #Now try using the iterator - note that all these
        #test cases have only one record.

        # With the SequenceParser
        test_handle = open(datafile)
        records = list(SeqIO.parse(test_handle, "swiss"))
        test_handle.close()

        self.assertEqual(len(records), 1)
        self.assertTrue(isinstance(records[0], SeqRecord))

        #Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].seq.tostring(), seq_record.seq.tostring())
        self.assertEqual(records[0].description, seq_record.description)
        self.assertEqual(records[0].name, seq_record.name)
        self.assertEqual(records[0].id, seq_record.id)

        # With the RecordParser
        test_handle = open(datafile)
        records = list(SwissProt.parse(test_handle))
        test_handle.close()

        self.assertEqual(len(records), 1)
        self.assertTrue(isinstance(records[0], SwissProt.Record))

        #Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].sequence, record.sequence)
        self.assertEqual(records[0].description, record.description)
        self.assertEqual(records[0].entry_name, record.entry_name)
        self.assertEqual(records[0].accessions, record.accessions)


    def test_sp012(self):
        "Parsing SwissProt file sp012"

        filename = 'sp012'
        # test the record parser

        datafile = os.path.join('SwissProt', filename)

        test_handle = open(datafile)
        seq_record = SeqIO.read(test_handle, "swiss")
        test_handle.close()

        self.assertTrue(isinstance(seq_record, SeqRecord))

        self.assertEqual(seq_record.id, "Q9Y736")
        self.assertEqual(seq_record.name, "Q9Y736")
        self.assertEqual(seq_record.description, "UBIQUITIN.")
        self.assertEqual(repr(seq_record.seq), "Seq('MQIFVKTLTGKTITLEVESSDTIDNVKTKIQDKEGIPPDQQRLIFAGKQLEDGR...GGN', ProteinAlphabet())")

        test_handle = open(datafile)
        record = SwissProt.read(test_handle)
        test_handle.close()

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(record.entry_name, "Q9Y736")
        self.assertEqual(record.accessions, ['Q9Y736'])
        self.assertEqual(record.organism_classification, ['Eukaryota', 'Fungi', 'Ascomycota', 'Pezizomycotina', 'Eurotiomycetes', 'Onygenales', 'Arthrodermataceae', 'mitosporic Arthrodermataceae', 'Trichophyton'])
        self.assertEqual(record.seqinfo, (153, 17238, '01153CF30C2DEDFF'))
    
        self.assertEqual(len(record.features), 0)

        self.assertEqual(len(record.references), 2)
        self.assertEqual(record.references[0].authors, "Kano R., Nakamura Y., Watanabe S., Hasegawa A.")
        self.assertEqual(record.references[0].title, "Trichophyton mentagrophytes mRNA for ubiquitin.")
        self.assertEqual(len(record.references[0].references), 0)
        self.assertEqual(record.references[1].authors, "Kano R.")
        self.assertEqual(record.references[1].title, "Microsporum canis mRNA for ubiquitin, complete cds.")
        self.assertEqual(len(record.references[1].references), 0)

        #Check the two parsers agree on the essentials
        self.assertEqual(seq_record.seq.tostring(), record.sequence)
        self.assertEqual(seq_record.description, record.description)
        self.assertEqual(seq_record.name, record.entry_name)
        self.assertTrue(seq_record.id in record.accessions)

        #Now try using the iterator - note that all these
        #test cases have only one record.

        # With the SequenceParser
        test_handle = open(datafile)
        records = list(SeqIO.parse(test_handle, "swiss"))
        test_handle.close()

        self.assertEqual(len(records), 1)
        self.assertTrue(isinstance(records[0], SeqRecord))

        #Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].seq.tostring(), seq_record.seq.tostring())
        self.assertEqual(records[0].description, seq_record.description)
        self.assertEqual(records[0].name, seq_record.name)
        self.assertEqual(records[0].id, seq_record.id)

        # With the RecordParser
        test_handle = open(datafile)
        records = list(SwissProt.parse(test_handle))
        test_handle.close()

        self.assertEqual(len(records), 1)
        self.assertTrue(isinstance(records[0], SwissProt.Record))

        #Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].sequence, record.sequence)
        self.assertEqual(records[0].description, record.description)
        self.assertEqual(records[0].entry_name, record.entry_name)
        self.assertEqual(records[0].accessions, record.accessions)


    def test_sp013(self):
        "Parsing SwissProt file sp013"

        filename = 'sp013'
        # test the record parser

        datafile = os.path.join('SwissProt', filename)

        test_handle = open(datafile)
        seq_record = SeqIO.read(test_handle, "swiss")
        test_handle.close()

        self.assertTrue(isinstance(seq_record, SeqRecord))

        self.assertEqual(seq_record.id, "P82909")
        self.assertEqual(seq_record.name, "P82909")
        self.assertEqual(seq_record.description, "MITOCHONDRIAL 28S RIBOSOMAL PROTEIN S36 (MRP-S36).")
        self.assertEqual(repr(seq_record.seq), "Seq('MGSKMASASRVVQVVKPHTPLIRFPDRRDNPKPNVSEALRSAGLPSHSSVISQH...GPE', ProteinAlphabet())")

        test_handle = open(datafile)
        record = SwissProt.read(test_handle)
        test_handle.close()

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(record.entry_name, "P82909")
        self.assertEqual(record.accessions, ['P82909'])
        self.assertEqual(record.organism_classification, ['Eukaryota', 'Metazoa', 'Chordata', 'Craniata', 'Vertebrata', 'Euteleostomi', 'Mammalia', 'Eutheria', 'Primates', 'Catarrhini', 'Hominidae', 'Homo'])
        self.assertEqual(record.seqinfo, (102, 11335, '83EF107B42E2FCFD'))
    
        self.assertEqual(len(record.features), 0)

        self.assertEqual(len(record.references), 2)
        self.assertEqual(record.references[0].authors, "Strausberg R.")
        self.assertEqual(record.references[0].title, "")
        self.assertEqual(len(record.references[0].references), 0)
        self.assertEqual(record.references[1].authors, "Koc E.C., Burkhart W., Blackburn K., Moseley A., Spremulli L.L.")
        self.assertEqual(record.references[1].title, "The small subunit of the mammalian mitochondrial ribosome. Identification of the full complement ribosomal proteins present.")
        self.assertEqual(len(record.references[1].references), 0)

        #Check the two parsers agree on the essentials
        self.assertEqual(seq_record.seq.tostring(), record.sequence)
        self.assertEqual(seq_record.description, record.description)
        self.assertEqual(seq_record.name, record.entry_name)
        self.assertTrue(seq_record.id in record.accessions)

        #Now try using the iterator - note that all these
        #test cases have only one record.

        # With the SequenceParser
        test_handle = open(datafile)
        records = list(SeqIO.parse(test_handle, "swiss"))
        test_handle.close()

        self.assertEqual(len(records), 1)
        self.assertTrue(isinstance(records[0], SeqRecord))

        #Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].seq.tostring(), seq_record.seq.tostring())
        self.assertEqual(records[0].description, seq_record.description)
        self.assertEqual(records[0].name, seq_record.name)
        self.assertEqual(records[0].id, seq_record.id)

        # With the RecordParser
        test_handle = open(datafile)
        records = list(SwissProt.parse(test_handle))
        test_handle.close()

        self.assertEqual(len(records), 1)
        self.assertTrue(isinstance(records[0], SwissProt.Record))

        #Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].sequence, record.sequence)
        self.assertEqual(records[0].description, record.description)
        self.assertEqual(records[0].entry_name, record.entry_name)
        self.assertEqual(records[0].accessions, record.accessions)


    def test_sp014(self):
        "Parsing SwissProt file sp014"

        filename = 'sp014'
        # test the record parser

        datafile = os.path.join('SwissProt', filename)

        test_handle = open(datafile)
        seq_record = SeqIO.read(test_handle, "swiss")
        test_handle.close()

        self.assertTrue(isinstance(seq_record, SeqRecord))

        self.assertEqual(seq_record.id, "P12166")
        self.assertEqual(seq_record.name, "PSBL_ORYSA")
        self.assertEqual(seq_record.description, "PHOTOSYSTEM II REACTION CENTER L PROTEIN (PSII 5 KDA PROTEIN).")
        self.assertEqual(repr(seq_record.seq), "Seq('TQSNPNEQNVELNRTSLYWGLLLIFVLAVLFSNYFFN', ProteinAlphabet())")

        test_handle = open(datafile)
        record = SwissProt.read(test_handle)
        test_handle.close()

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(record.entry_name, "PSBL_ORYSA")
        self.assertEqual(record.accessions, ['P12166', 'P12167', 'Q34007'])
        self.assertEqual(record.organism_classification, ['Eukaryota', 'Viridiplantae', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Poales', 'Poaceae', 'Ehrhartoideae', 'Oryzeae', 'Oryza'])
        self.assertEqual(record.seqinfo, (37, 4366, 'CC537AEC50B2C784'))
    
        self.assertEqual(len(record.features), 1)
        self.assertEqual(record.features[0], ('INIT_MET', 0, 0, '', ''))

        self.assertEqual(len(record.references), 22)
        self.assertEqual(record.references[0].authors, "Sugiura M.")
        self.assertEqual(record.references[0].title, "")
        self.assertEqual(len(record.references[0].references), 0)
        self.assertEqual(record.references[1].authors, "Hiratsuka J., Shimada H., Whittier R., Ishibashi T., Sakamoto M., Mori M., Kondo C., Honji Y., Sun C.-R., Meng B.-Y., Li Y.-Q., Kanno A., Nishizawa Y., Hirai A., Shinozaki K., Sugiura M.")
        self.assertEqual(record.references[1].title, "The complete sequence of the rice (Oryza sativa) chloroplast genome: intermolecular recombination between distinct tRNA genes accounts for a major plastid DNA inversion during the evolution of the cereals.")
        self.assertEqual(len(record.references[1].references), 2)
        self.assertEqual(record.references[1].references[0], ('MEDLINE', '89364698'))
        self.assertEqual(record.references[1].references[1], ('PubMed', '2770692'))
        self.assertEqual(record.references[2].authors, "Sugiura M.")
        self.assertEqual(record.references[2].title, "")
        self.assertEqual(len(record.references[2].references), 0)
        self.assertEqual(record.references[3].authors, "Shinozaki K., Ohme M., Tanaka M., Wakasugi T., Hayashida N., Matsubayashi T., Zaita N., Chunwongse J., Obokata J., Yamaguchi-Shinozaki K., Ohto C., Torazawa K., Meng B.Y., Sugita M., Deno H., Kamogashira T., Yamada K., Kusuda J., Takaiwa F., Kato A., Tohdoh N., Shimada H., Sugiura M.")
        self.assertEqual(record.references[3].title, "The complete nucleotide sequence of the tobacco chloroplast genome: its gene organization and expression.")
        self.assertEqual(len(record.references[3].references), 0)
        self.assertEqual(record.references[4].authors, "Chaudhuri S., Maliga P.")
        self.assertEqual(record.references[4].title, "Sequences directing C to U editing of the plastid psbL mRNA are located within a 22 nucleotide segment spanning the editing site.")
        self.assertEqual(len(record.references[4].references), 2)
        self.assertEqual(record.references[4].references[0], ('MEDLINE', '97076156'))
        self.assertEqual(record.references[4].references[1], ('PubMed', '8918473'))
        self.assertEqual(record.references[5].authors, "Chakhmakhcheva O.G., Andreeva A.V., Buryakova A.A., Reverdatto S.V., Efimov V.A.")
        self.assertEqual(record.references[5].title, "Nucleotide sequence of the barley chloroplast psbE, psbF genes and flanking regions.")
        self.assertEqual(len(record.references[5].references), 2)
        self.assertEqual(record.references[5].references[0], ('MEDLINE', '89240046'))
        self.assertEqual(record.references[5].references[1], ('PubMed', '2654886'))
        self.assertEqual(record.references[6].authors, "Efimov V.A., Andreeva A.V., Reverdatto S.V., Chakhmakhcheva O.G.")
        self.assertEqual(record.references[6].title, "Photosystem II of rye. Nucleotide sequence of the psbB, psbC, psbE, psbF, psbH genes of rye and chloroplast DNA regions adjacent to them.")
        self.assertEqual(len(record.references[6].references), 2)
        self.assertEqual(record.references[6].references[0], ('MEDLINE', '92207253'))
        self.assertEqual(record.references[6].references[1], ('PubMed', '1804121'))
        self.assertEqual(record.references[7].authors, "Webber A.N., Hird S.M., Packman L.C., Dyer T.A., Gray J.C.")
        self.assertEqual(record.references[7].title, "A photosystem II polypeptide is encoded by an open reading frame co-transcribed with genes for cytochrome b-559 in wheat chloroplast DNA.")
        self.assertEqual(len(record.references[7].references), 0)
        self.assertEqual(record.references[8].authors, "Kudla J., Igloi G.L., Metzlaff M., Hagemann R., Koessel H.")
        self.assertEqual(record.references[8].title, "RNA editing in tobacco chloroplasts leads to the formation of a translatable psbL mRNA by a C to U substitution within the initiation codon.")
        self.assertEqual(len(record.references[8].references), 2)
        self.assertEqual(record.references[8].references[0], ('MEDLINE', '92191997'))
        self.assertEqual(record.references[8].references[1], ('PubMed', '1547774'))
        self.assertEqual(record.references[9].authors, "Zolotarev A.S., Kolosov V.L.")
        self.assertEqual(record.references[9].title, "Nucleotide sequence of the rye chloroplast DNA fragment, comprising psbE and psbF genes.")
        self.assertEqual(len(record.references[9].references), 2)
        self.assertEqual(record.references[9].references[0], ('MEDLINE', '89160331'))
        self.assertEqual(record.references[9].references[1], ('PubMed', '2646599'))
        self.assertEqual(record.references[10].authors, "Kolosov V.L., Klezovich O.N., Abdulaev N.G., Zolotarev A.S.")
        self.assertEqual(record.references[10].title, "Photosystem II of rye. Nucleotide sequence of genes psbE, psbF, psbL and OPC40 of chloroplast DNA.")
        self.assertEqual(len(record.references[10].references), 2)
        self.assertEqual(record.references[10].references[0], ('MEDLINE', '90073796'))
        self.assertEqual(record.references[10].references[1], ('PubMed', '2686655'))
        self.assertEqual(record.references[11].authors, "Haley J., Bogorad L.")
        self.assertEqual(record.references[11].title, "")
        self.assertEqual(len(record.references[11].references), 0)
        self.assertEqual(record.references[12].authors, "Maier R.M., Neckermann K., Igloi G.L., Koessel H.")
        self.assertEqual(record.references[12].title, "Complete sequence of the maize chloroplast genome: gene content, hotspots of divergence and fine tuning of genetic information by transcript editing.")
        self.assertEqual(len(record.references[12].references), 2)
        self.assertEqual(record.references[12].references[0], ('MEDLINE', '95395841'))
        self.assertEqual(record.references[12].references[1], ('PubMed', '7666415'))
        self.assertEqual(record.references[13].authors, "Willey D.L., Gray J.C.")
        self.assertEqual(record.references[13].title, "Two small open reading frames are co-transcribed with the pea chloroplast genes for the polypeptides of cytochrome b-559.")
        self.assertEqual(len(record.references[13].references), 2)
        self.assertEqual(record.references[13].references[0], ('MEDLINE', '89354671'))
        self.assertEqual(record.references[13].references[1], ('PubMed', '2766383'))
        self.assertEqual(record.references[14].authors, "Bock R., Hagemann R., Koessel H., Kudla J.")
        self.assertEqual(record.references[14].title, "Tissue- and stage-specific modulation of RNA editing of the psbF and psbL transcript from spinach plastids -- a new regulatory mechanism?")
        self.assertEqual(len(record.references[14].references), 2)
        self.assertEqual(record.references[14].references[0], ('MEDLINE', '93360903'))
        self.assertEqual(record.references[14].references[1], ('PubMed', '8355656'))
        self.assertEqual(record.references[15].authors, "Hermann R.G., Alt J., Schiller B., Widger W.R., Cramer W.A.")
        self.assertEqual(record.references[15].title, "Nucleotide sequence of the gene for apocytochrome b-559 on the spinach plastid chromosome: implications for the structure of the membrane protein.")
        self.assertEqual(len(record.references[15].references), 0)
        self.assertEqual(record.references[16].authors, "Kuntz M., Camara B., Weil J.-H., Schantz R.")
        self.assertEqual(record.references[16].title, "The psbL gene from bell pepper (Capsicum annuum): plastid RNA editing also occurs in non-photosynthetic chromoplasts.")
        self.assertEqual(len(record.references[16].references), 2)
        self.assertEqual(record.references[16].references[0], ('MEDLINE', '93099270'))
        self.assertEqual(record.references[16].references[1], ('PubMed', '1463853'))
        self.assertEqual(record.references[17].authors, "Forsthoefel N.R., Cushman J.C.")
        self.assertEqual(record.references[17].title, "Characterization and expression of photosystem II genes (psbE, psbF, and psbL) from the facultative crassulacean acid metabolism plant Mesembryanthemum crystallinum.")
        self.assertEqual(len(record.references[17].references), 2)
        self.assertEqual(record.references[17].references[0], ('MEDLINE', '94345017'))
        self.assertEqual(record.references[17].references[1], ('PubMed', '8066140'))
        self.assertEqual(record.references[18].authors, "Kubo T., Yanai Y., Kinoshita T., Mikami T.")
        self.assertEqual(record.references[18].title, "The chloroplast trnP-trnW-petG gene cluster in the mitochondrial genomes of Beta vulgaris, B. trigyna and B. webbiana: evolutionary aspects.")
        self.assertEqual(len(record.references[18].references), 2)
        self.assertEqual(record.references[18].references[0], ('MEDLINE', '95254673'))
        self.assertEqual(record.references[18].references[1], ('PubMed', '7736615'))
        self.assertEqual(record.references[19].authors, "Naithani S.")
        self.assertEqual(record.references[19].title, "")
        self.assertEqual(len(record.references[19].references), 0)
        self.assertEqual(record.references[20].authors, "Ikeuchi M., Takio K., Inoue Y.")
        self.assertEqual(record.references[20].title, "N-terminal sequencing of photosystem II low-molecular-mass proteins. 5 and 4.1 kDa components of the O2-evolving core complex from higher plants.")
        self.assertEqual(len(record.references[20].references), 2)
        self.assertEqual(record.references[20].references[0], ('MEDLINE', '89121082'))
        self.assertEqual(record.references[20].references[1], ('PubMed', '2644131'))
        self.assertEqual(record.references[21].authors, "Zheleva D., Sharma J., Panico M., Morris H.R., Barber J.")
        self.assertEqual(record.references[21].title, "Isolation and characterization of monomeric and dimeric CP47-reaction center photosystem II complexes.")
        self.assertEqual(len(record.references[21].references), 2)
        self.assertEqual(record.references[21].references[0], ('MEDLINE', '98298118'))
        self.assertEqual(record.references[21].references[1], ('PubMed', '9632665'))

        #Check the two parsers agree on the essentials
        self.assertEqual(seq_record.seq.tostring(), record.sequence)
        self.assertEqual(seq_record.description, record.description)
        self.assertEqual(seq_record.name, record.entry_name)
        self.assertTrue(seq_record.id in record.accessions)

        #Now try using the iterator - note that all these
        #test cases have only one record.

        # With the SequenceParser
        test_handle = open(datafile)
        records = list(SeqIO.parse(test_handle, "swiss"))
        test_handle.close()

        self.assertEqual(len(records), 1)
        self.assertTrue(isinstance(records[0], SeqRecord))

        #Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].seq.tostring(), seq_record.seq.tostring())
        self.assertEqual(records[0].description, seq_record.description)
        self.assertEqual(records[0].name, seq_record.name)
        self.assertEqual(records[0].id, seq_record.id)

        # With the RecordParser
        test_handle = open(datafile)
        records = list(SwissProt.parse(test_handle))
        test_handle.close()

        self.assertEqual(len(records), 1)
        self.assertTrue(isinstance(records[0], SwissProt.Record))

        #Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].sequence, record.sequence)
        self.assertEqual(records[0].description, record.description)
        self.assertEqual(records[0].entry_name, record.entry_name)
        self.assertEqual(records[0].accessions, record.accessions)


    def test_sp015(self):
        "Parsing SwissProt file sp015"

        filename = 'sp015'
        # test the record parser

        datafile = os.path.join('SwissProt', filename)

        test_handle = open(datafile)
        seq_record = SeqIO.read(test_handle, "swiss")
        test_handle.close()

        self.assertTrue(isinstance(seq_record, SeqRecord))

        self.assertEqual(seq_record.id, "IPI00383150")
        self.assertEqual(seq_record.name, "IPI00383150.2")
        self.assertEqual(seq_record.description, "")
        self.assertEqual(repr(seq_record.seq), "Seq('MSFQAPRRLLELAGQSLLRDQALAISVLDELPRELFPRLFVEAFTSRRCEVLKV...TPC', ProteinAlphabet())")

        test_handle = open(datafile)
        record = SwissProt.read(test_handle)
        test_handle.close()

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(record.entry_name, "IPI00383150.2")
        self.assertEqual(record.accessions, ['IPI00383150'])
        self.assertEqual(record.organism_classification, ['Eukaryota', 'Metazoa', 'Chordata', 'Craniata', 'Vertebrata', 'Euteleostomi', 'Mammalia', 'Eutheria', 'Primates', 'Catarrhini', 'Hominidae', 'Homo'])
        self.assertEqual(record.seqinfo, (457, 52856, '5C3151AAADBDE232'))
    
        self.assertEqual(len(record.features), 0)
        self.assertEqual(len(record.references), 0)

        #Check the two parsers agree on the essentials
        self.assertEqual(seq_record.seq.tostring(), record.sequence)
        self.assertEqual(seq_record.description, record.description)
        self.assertEqual(seq_record.name, record.entry_name)
        self.assertTrue(seq_record.id in record.accessions)

        #Now try using the iterator - note that all these
        #test cases have only one record.

        # With the SequenceParser
        test_handle = open(datafile)
        records = list(SeqIO.parse(test_handle, "swiss"))
        test_handle.close()

        self.assertEqual(len(records), 1)
        self.assertTrue(isinstance(records[0], SeqRecord))

        #Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].seq.tostring(), seq_record.seq.tostring())
        self.assertEqual(records[0].description, seq_record.description)
        self.assertEqual(records[0].name, seq_record.name)
        self.assertEqual(records[0].id, seq_record.id)

        # With the RecordParser
        test_handle = open(datafile)
        records = list(SwissProt.parse(test_handle))
        test_handle.close()

        self.assertEqual(len(records), 1)
        self.assertTrue(isinstance(records[0], SwissProt.Record))

        #Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].sequence, record.sequence)
        self.assertEqual(records[0].description, record.description)
        self.assertEqual(records[0].entry_name, record.entry_name)
        self.assertEqual(records[0].accessions, record.accessions)



if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
