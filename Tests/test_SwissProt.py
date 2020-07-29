# Copyright 2009 by Michiel de Hoon.  All rights reserved.
# Revisions copyright 2010 by Peter Cock. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Test for the SwissProt parser on SwissProt files."""
import os
import unittest

from Bio import SeqIO
from Bio import SwissProt
from Bio.SeqRecord import SeqRecord


class TestSwissProt(unittest.TestCase):
    def test_sp001(self):
        """Parsing SwissProt file sp001."""
        filename = "sp001"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)

        with open(datafile) as test_handle:
            seq_record = SeqIO.read(test_handle, "swiss")

        self.assertIsInstance(seq_record, SeqRecord)

        self.assertEqual(seq_record.id, "Q13454")
        self.assertEqual(seq_record.name, "N33_HUMAN")
        self.assertEqual(seq_record.description, "N33 PROTEIN.")
        self.assertEqual(
            repr(seq_record.seq),
            "Seq('MGARGAPSRRRQAGRRLRYLPTGSFPFLLLLLLLCIQLGGGQKKKENLLAEKVE...DFE')",
        )

        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(record.entry_name, "N33_HUMAN")
        self.assertEqual(record.accessions, ["Q13454", "Q14911", "Q14912"])
        self.assertEqual(
            record.organism_classification,
            [
                "Eukaryota",
                "Metazoa",
                "Chordata",
                "Craniata",
                "Vertebrata",
                "Mammalia",
                "Eutheria",
                "Primates",
                "Catarrhini",
                "Hominidae",
                "Homo",
            ],
        )
        self.assertEqual(record.seqinfo, (348, 39676, "75818910"))

        self.assertEqual(len(record.features), 6)
        feature = record.features[0]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 19)
        self.assertEqual(feature.location.end, 40)
        self.assertEqual(feature.qualifiers["description"], "POTENTIAL.")
        self.assertIsNone(feature.id)
        feature = record.features[1]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 196)
        self.assertEqual(feature.location.end, 217)
        self.assertEqual(feature.qualifiers["description"], "POTENTIAL.")
        self.assertIsNone(feature.id)
        feature = record.features[2]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 221)
        self.assertEqual(feature.location.end, 242)
        self.assertEqual(feature.qualifiers["description"], "POTENTIAL.")
        self.assertIsNone(feature.id)
        feature = record.features[3]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 276)
        self.assertEqual(feature.location.end, 297)
        self.assertEqual(feature.qualifiers["description"], "POTENTIAL.")
        self.assertIsNone(feature.id)
        feature = record.features[4]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 312)
        self.assertEqual(feature.location.end, 333)
        self.assertEqual(feature.qualifiers["description"], "POTENTIAL.")
        self.assertIsNone(feature.id)
        feature = record.features[5]
        self.assertEqual(feature.type, "VARSPLIC")
        self.assertEqual(feature.location.start, 343)
        self.assertEqual(feature.location.end, 348)
        self.assertEqual(
            feature.qualifiers["description"], "DLDFE -> FLIK (IN FORM 2)."
        )
        self.assertIsNone(feature.id)

        self.assertEqual(len(record.references), 1)
        self.assertEqual(
            record.references[0].authors,
            "MACGROGAN D., LEVY A., BOVA G.S., ISAACS W.B., BOOKSTEIN R.",
        )
        self.assertEqual(
            record.references[0].title,
            "Structure and methylation-associated silencing of a gene within a homozygously deleted region of human chromosome band 8p22.",
        )
        self.assertEqual(len(record.references[0].references), 1)
        self.assertEqual(record.references[0].references[0], ("MEDLINE", "96299740"))

        # Check that the two parsers agree on the essentials
        self.assertEqual(str(seq_record.seq), record.sequence)
        self.assertEqual(seq_record.description, record.description)
        self.assertEqual(seq_record.name, record.entry_name)
        self.assertIn(seq_record.id, record.accessions)

        # Now try using the iterator - note that all these
        # test cases have only one record.

        # With the SequenceParser
        with open(datafile) as test_handle:
            records = list(SeqIO.parse(test_handle, "swiss"))

        self.assertEqual(len(records), 1)
        self.assertIsInstance(records[0], SeqRecord)

        # Check matches what we got earlier without the iterator:
        self.assertEqual(str(records[0].seq), str(seq_record.seq))
        self.assertEqual(records[0].description, seq_record.description)
        self.assertEqual(records[0].name, seq_record.name)
        self.assertEqual(records[0].id, seq_record.id)

        # With the RecordParser
        with open(datafile) as test_handle:
            records = list(SwissProt.parse(test_handle))

        self.assertEqual(len(records), 1)
        self.assertIsInstance(records[0], SwissProt.Record)

        # Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].sequence, record.sequence)
        self.assertEqual(records[0].description, record.description)
        self.assertEqual(records[0].entry_name, record.entry_name)
        self.assertEqual(records[0].accessions, record.accessions)

    def test_sp002(self):
        """Parsing SwissProt file sp002."""
        filename = "sp002"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)

        with open(datafile) as test_handle:
            seq_record = SeqIO.read(test_handle, "swiss")

        self.assertIsInstance(seq_record, SeqRecord)

        self.assertEqual(seq_record.id, "P54101")
        self.assertEqual(seq_record.name, "CSP_MOUSE")
        self.assertEqual(seq_record.description, "CYSTEINE STRING PROTEIN (CSP).")
        self.assertEqual(
            repr(seq_record.seq),
            "Seq('MADQRQRSLSTSGESLYHVLGLDKNATSDDIKKSYRKLALKYHPDKNPDNPEAA...GFN')",
        )

        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(record.entry_name, "CSP_MOUSE")
        self.assertEqual(record.accessions, ["P54101"])
        self.assertEqual(
            record.organism_classification,
            [
                "Eukaryota",
                "Metazoa",
                "Chordata",
                "Craniata",
                "Vertebrata",
                "Mammalia",
                "Eutheria",
                "Rodentia",
                "Sciurognathi",
                "Muridae",
                "Murinae",
                "Mus",
            ],
        )
        self.assertEqual(record.seqinfo, (198, 22100, "9DF0142B"))

        self.assertEqual(len(record.features), 2)
        feature = record.features[0]
        self.assertEqual(feature.type, "DOMAIN")
        self.assertEqual(feature.location.start, 12)
        self.assertEqual(feature.location.end, 82)
        self.assertEqual(feature.qualifiers["description"], "DNAJ-LIKE.")
        self.assertIsNone(feature.id)
        feature = record.features[1]
        self.assertEqual(feature.type, "DOMAIN")
        self.assertEqual(feature.location.start, 117)
        self.assertEqual(feature.location.end, 128)
        self.assertEqual(feature.qualifiers["description"], "POLY-CYS.")
        self.assertIsNone(feature.id)

        self.assertEqual(len(record.references), 3)
        self.assertEqual(record.references[0].authors, "QIN N., LIN T., BIRNBAUMER L.")
        self.assertEqual(record.references[0].title, "")
        self.assertEqual(len(record.references[0].references), 0)
        self.assertEqual(
            record.references[1].authors, "MASTROGIACOMO A., GUNDERSEN C.B."
        )
        self.assertEqual(
            record.references[1].title,
            "The nucleotide and deduced amino acid sequence of a rat cysteine string protein.",
        )
        self.assertEqual(len(record.references[1].references), 1)
        self.assertEqual(record.references[1].references[0], ("MEDLINE", "95223109"))
        self.assertEqual(record.references[2].authors, "BRAUN J.E., SCHELLER R.H.")
        self.assertEqual(
            record.references[2].title,
            "Cysteine string protein, a DnaJ family member, is present on diverse secretory vesicles.",
        )
        self.assertEqual(len(record.references[2].references), 1)
        self.assertEqual(record.references[2].references[0], ("MEDLINE", "96188189"))

        # Check the two parsers agree on the essentials
        self.assertEqual(str(seq_record.seq), record.sequence)
        self.assertEqual(seq_record.description, record.description)
        self.assertEqual(seq_record.name, record.entry_name)
        self.assertIn(seq_record.id, record.accessions)

        # Now try using the iterator - note that all these
        # test cases have only one record.

        # With the SequenceParser
        with open(datafile) as test_handle:
            records = list(SeqIO.parse(test_handle, "swiss"))

        self.assertEqual(len(records), 1)
        self.assertIsInstance(records[0], SeqRecord)

        # Check matches what we got earlier without the iterator:
        self.assertEqual(str(records[0].seq), str(seq_record.seq))
        self.assertEqual(records[0].description, seq_record.description)
        self.assertEqual(records[0].name, seq_record.name)
        self.assertEqual(records[0].id, seq_record.id)

        # With the RecordParser
        with open(datafile) as test_handle:
            records = list(SwissProt.parse(test_handle))

        self.assertEqual(len(records), 1)
        self.assertIsInstance(records[0], SwissProt.Record)

        # Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].sequence, record.sequence)
        self.assertEqual(records[0].description, record.description)
        self.assertEqual(records[0].entry_name, record.entry_name)
        self.assertEqual(records[0].accessions, record.accessions)

    def test_sp003(self):
        """Parsing SwissProt file sp003."""
        filename = "sp003"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)

        with open(datafile) as test_handle:
            seq_record = SeqIO.read(test_handle, "swiss")

        self.assertIsInstance(seq_record, SeqRecord)

        self.assertEqual(seq_record.id, "P42655")
        self.assertEqual(seq_record.name, "143E_HUMAN")
        self.assertEqual(
            seq_record.description,
            "14-3-3 PROTEIN EPSILON (MITOCHONDRIAL IMPORT STIMULATION FACTOR L SUBUNIT) (PROTEIN KINASE C INHIBITOR PROTEIN-1) (KCIP-1) (14-3-3E).",
        )
        self.assertEqual(
            repr(seq_record.seq),
            "Seq('MDDREDLVYQAKLAEQAERYDEMVESMKKVAGMDVELTVEERNLLSVAYKNVIG...ENQ')",
        )

        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(record.entry_name, "143E_HUMAN")
        self.assertEqual(record.accessions, ["P42655", "P29360", "Q63631"])
        self.assertEqual(
            record.organism_classification,
            [
                "EUKARYOTA",
                "METAZOA",
                "CHORDATA",
                "VERTEBRATA",
                "MAMMALIA",
                "EUTHERIA",
                "PRIMATES",
                "CATARRHINI",
                "HOMINIDAE",
                "HOMO",
            ],
        )
        self.assertEqual(record.seqinfo, (255, 29174, "40A43E62"))

        self.assertEqual(len(record.features), 5)
        feature = record.features[0]
        self.assertEqual(feature.type, "MOD_RES")
        self.assertEqual(feature.location.start, 0)
        self.assertEqual(feature.location.end, 1)
        self.assertEqual(feature.qualifiers["description"], "ACETYLATION.")
        self.assertIsNone(feature.id)
        feature = record.features[1]
        self.assertEqual(feature.type, "CONFLICT")
        self.assertEqual(feature.location.start, 72)
        self.assertEqual(feature.location.end, 73)
        self.assertEqual(feature.qualifiers["description"], "K -> T (IN REF. 8).")
        self.assertIsNone(feature.id)
        feature = record.features[2]
        self.assertEqual(feature.type, "CONFLICT")
        self.assertEqual(feature.location.start, 119)
        self.assertEqual(feature.location.end, 120)
        self.assertEqual(feature.qualifiers["description"], "F -> S (IN REF. 8).")
        self.assertIsNone(feature.id)
        feature = record.features[3]
        self.assertEqual(feature.type, "CONFLICT")
        self.assertEqual(feature.location.start, 122)
        self.assertEqual(feature.location.end, 123)
        self.assertEqual(feature.qualifiers["description"], "K -> Y (IN REF. 8).")
        self.assertIsNone(feature.id)
        feature = record.features[4]
        self.assertEqual(feature.type, "CONFLICT")
        self.assertEqual(feature.location.start, 128)
        self.assertEqual(feature.location.end, 129)
        self.assertEqual(feature.qualifiers["description"], "H -> Y (IN REF. 13).")
        self.assertIsNone(feature.id)

        self.assertEqual(len(record.references), 13)
        self.assertEqual(
            record.references[0].authors, "CONKLIN D.S., GALAKTIONOV K., BEACH D."
        )
        self.assertEqual(
            record.references[0].title,
            "14-3-3 proteins associate with cdc25 phosphatases.",
        )
        self.assertEqual(len(record.references[0].references), 1)
        self.assertEqual(record.references[0].references[0], ("MEDLINE", "95372385"))
        self.assertEqual(
            record.references[1].authors, "LUK S.C.W., LEE C.Y., WAYE M.M.Y."
        )
        self.assertEqual(record.references[1].title, "")
        self.assertEqual(len(record.references[1].references), 0)
        self.assertEqual(
            record.references[2].authors, "JIN D.Y., LYU M.S., KOZAK C.A., JEANG K.T."
        )
        self.assertEqual(record.references[2].title, "Function of 14-3-3 proteins.")
        self.assertEqual(len(record.references[2].references), 1)
        self.assertEqual(record.references[2].references[0], ("MEDLINE", "96300316"))
        self.assertEqual(
            record.references[3].authors,
            "CHONG S.S., TANIGAMI A., ROSCHKE A.V., LEDBETTER D.H.",
        )
        self.assertEqual(
            record.references[3].title,
            "14-3-3 epsilon has no homology to LIS1 and lies telomeric to it on chromosome 17p13.3 outside the Miller-Dieker syndrome chromosome region.",
        )
        self.assertEqual(len(record.references[3].references), 1)
        self.assertEqual(record.references[3].references[0], ("MEDLINE", "97011338"))
        self.assertEqual(
            record.references[4].authors, "TANIGAMI A., CHONG S.S., LEDBETTER D.H."
        )
        self.assertEqual(record.references[4].title, "14-3-3 epsilon genomic sequence.")
        self.assertEqual(len(record.references[4].references), 0)
        self.assertEqual(
            record.references[5].authors,
            "ROSEBOOM P.H., WELLER J.L., BABILA T., AITKEN A., SELLERS L.A., MOFFET J.R., NAMBOODIRI M.A., KLEIN D.C.",
        )
        self.assertEqual(
            record.references[5].title,
            "Cloning and characterization of the epsilon and zeta isoforms of the 14-3-3 proteins.",
        )
        self.assertEqual(len(record.references[5].references), 1)
        self.assertEqual(record.references[5].references[0], ("MEDLINE", "94296566"))
        self.assertEqual(
            record.references[6].authors,
            "ALAM R., HACHIYA N., SAKAGUCHI M., SHUN-ICHIRO K., IWANAGA S., KITAJIMA M., MIHARA K., OMURA T.",
        )
        self.assertEqual(
            record.references[6].title,
            "cDNA cloning and characterization of mitochondrial import stimulation factor (MSF) purified from rat liver cytosol.",
        )
        self.assertEqual(len(record.references[6].references), 1)
        self.assertEqual(record.references[6].references[0], ("MEDLINE", "95122474"))
        self.assertEqual(
            record.references[7].authors, "GAO L., GU X.B., YU D.S., YU R.K., ZENG G."
        )
        self.assertEqual(
            record.references[7].title,
            "Association of a 14-3-3 protein with CMP-NeuAc:GM1 alpha 2,3-sialyltransferase.",
        )
        self.assertEqual(len(record.references[7].references), 1)
        self.assertEqual(record.references[7].references[0], ("MEDLINE", "96280718"))
        self.assertEqual(
            record.references[8].authors, "MCCONNELL J.E., ARMSTRONG J.F., BARD J.B."
        )
        self.assertEqual(
            record.references[8].title,
            "The mouse 14-3-3 epsilon isoform, a kinase regulator whose expression pattern is modulated in mesenchyme and neuronal differentiation.",
        )
        self.assertEqual(len(record.references[8].references), 1)
        self.assertEqual(record.references[8].references[0], ("MEDLINE", "95269876"))
        self.assertEqual(
            record.references[9].authors,
            "TAKIHARA Y., IRIE K., NOMURA M., MOTALEB M., MATSUMOTO K., SHIMADA K.",
        )
        self.assertEqual(record.references[9].title, "")
        self.assertEqual(len(record.references[9].references), 0)
        self.assertEqual(
            record.references[10].authors,
            "JONES J.M., NIIKURA T., PINKE R.M., GUO W., MOLDAY L., LEYKAM J., MCCONNELL D.G.",
        )
        self.assertEqual(
            record.references[10].title,
            "Expression of 14-3-3 proteins in bovine retinal photoreceptors.",
        )
        self.assertEqual(len(record.references[10].references), 0)
        self.assertEqual(
            record.references[11].authors,
            "TOKER A., SELLERS L.A., AMESS B., PATEL Y., HARRIS A., AITKEN A.",
        )
        self.assertEqual(
            record.references[11].title,
            "Multiple isoforms of a protein kinase C inhibitor (KCIP-1/14-3-3) from sheep brain. Amino acid sequence of phosphorylated forms.",
        )
        self.assertEqual(len(record.references[11].references), 1)
        self.assertEqual(record.references[11].references[0], ("MEDLINE", "92283271"))
        self.assertEqual(
            record.references[12].authors,
            "TOKER A., ELLIS C.A., SELLERS L.A., AITKEN A.",
        )
        self.assertEqual(
            record.references[12].title,
            "Protein kinase C inhibitor proteins. Purification from sheep brain and sequence similarity to lipocortins and 14-3-3 protein.",
        )
        self.assertEqual(len(record.references[12].references), 1)
        self.assertEqual(record.references[12].references[0], ("MEDLINE", "90345949"))

        # Check the two parsers agree on the essentials
        self.assertEqual(str(seq_record.seq), record.sequence)
        self.assertEqual(seq_record.description, record.description)
        self.assertEqual(seq_record.name, record.entry_name)
        self.assertIn(seq_record.id, record.accessions)

        # Now try using the iterator - note that all these
        # test cases have only one record.

        # With the SequenceParser
        with open(datafile) as test_handle:
            records = list(SeqIO.parse(test_handle, "swiss"))

        self.assertEqual(len(records), 1)
        self.assertIsInstance(records[0], SeqRecord)

        # Check matches what we got earlier without the iterator:
        self.assertEqual(str(records[0].seq), str(seq_record.seq))
        self.assertEqual(records[0].description, seq_record.description)
        self.assertEqual(records[0].name, seq_record.name)
        self.assertEqual(records[0].id, seq_record.id)

        # With the RecordParser
        with open(datafile) as test_handle:
            records = list(SwissProt.parse(test_handle))

        self.assertEqual(len(records), 1)
        self.assertIsInstance(records[0], SwissProt.Record)

        # Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].sequence, record.sequence)
        self.assertEqual(records[0].description, record.description)
        self.assertEqual(records[0].entry_name, record.entry_name)
        self.assertEqual(records[0].accessions, record.accessions)

    def test_P0A186(self):
        """Parsing SwissProt file P0A186.txt."""
        filename = "P0A186.txt"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)

        with open(datafile) as test_handle:
            seq_record = SeqIO.read(test_handle, "swiss")

        self.assertIsInstance(seq_record, SeqRecord)

        self.assertEqual(seq_record.id, "P0A186")
        self.assertEqual(seq_record.name, "NDOA_PSEU8")
        self.assertEqual(
            seq_record.description,
            "RecName: Full=Naphthalene 1,2-dioxygenase system, ferredoxin component {ECO:0000303|PubMed:8226631};",
        )
        self.assertEqual(
            repr(seq_record.seq),
            "Seq('MTVKWIEAVALSDILEGDVLGVTVEGKELALYEVEGEIYATDNLCTHGSARMSD...DLS')",
        )

        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(record.entry_name, "NDOA_PSEU8")
        self.assertEqual(record.accessions, ["P0A186", "O07829", "P23082", "Q52123"])
        self.assertEqual(record.organism_classification, ["Bacteria", "Proteobacteria"])
        self.assertEqual(record.seqinfo, (104, 11446, "475625DCC3EDCD41"))

        self.assertEqual(len(record.features), 7)
        feature = record.features[0]
        self.assertEqual(feature.type, "INIT_MET")
        self.assertEqual(feature.location.start, 0)
        self.assertEqual(feature.location.end, 1)
        self.assertEqual(feature.qualifiers["note"], "Removed")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000250|UniProtKB:P0A185")
        feature = record.features[1]
        self.assertEqual(feature.id, "PRO_0000201694")
        self.assertEqual(feature.type, "CHAIN")
        self.assertEqual(feature.location.start, 1)
        self.assertEqual(feature.location.end, 104)
        self.assertEqual(
            feature.qualifiers["note"],
            "Naphthalene 1,2-dioxygenase system, ferredoxin component",
        )
        feature = record.features[2]
        self.assertEqual(feature.type, "DOMAIN")
        self.assertEqual(feature.location.start, 5)
        self.assertEqual(feature.location.end, 101)
        self.assertEqual(feature.qualifiers["note"], "Rieske")
        self.assertEqual(
            feature.qualifiers["evidence"], "ECO:0000255|PROSITE-ProRule:PRU00628"
        )
        feature = record.features[3]
        self.assertEqual(feature.type, "METAL")
        self.assertEqual(feature.location.start, 44)
        self.assertEqual(feature.location.end, 45)
        self.assertEqual(feature.qualifiers["note"], "Iron-sulfur (2Fe-2S)")
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000250|UniProtKB:P0A185,ECO:0000255|PROSITE-ProRule:PRU00628",
        )
        self.assertIsNone(feature.id)
        feature = record.features[4]
        self.assertEqual(feature.type, "METAL")
        self.assertEqual(feature.location.start, 46)
        self.assertEqual(feature.location.end, 47)
        self.assertEqual(
            feature.qualifiers["note"], "Iron-sulfur (2Fe-2S); via pros nitrogen"
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000250|UniProtKB:P0A185,ECO:0000255|PROSITE-ProRule:PRU00628",
        )
        self.assertIsNone(feature.id)
        feature = record.features[5]
        self.assertEqual(feature.type, "METAL")
        self.assertEqual(feature.location.start, 63)
        self.assertEqual(feature.location.end, 64)
        self.assertEqual(feature.qualifiers["note"], "Iron-sulfur (2Fe-2S)")
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000250|UniProtKB:P0A185,ECO:0000255|PROSITE-ProRule:PRU00628",
        )
        self.assertIsNone(feature.id)
        feature = record.features[6]
        self.assertEqual(feature.type, "METAL")
        self.assertEqual(feature.location.start, 66)
        self.assertEqual(feature.location.end, 67)
        self.assertEqual(
            feature.qualifiers["note"], "Iron-sulfur (2Fe-2S); via pros nitrogen"
        )
        self.assertIsNone(feature.id)

        self.assertEqual(len(record.references), 1)
        self.assertEqual(
            record.references[0].authors,
            "Denome S.A., Stanley D.C., Olson E.S., Young K.D.",
        )
        self.assertEqual(
            record.references[0].title,
            "Metabolism of dibenzothiophene and naphthalene in Pseudomonas strains: complete DNA sequence of an upper naphthalene catabolic pathway.",
        )
        self.assertEqual(len(record.references[0].references), 2)
        self.assertEqual(record.references[0].references[0], ("PubMed", "8226631"))
        self.assertEqual(
            record.references[0].references[1],
            ("DOI", "10.1128/jb.175.21.6890-6901.1993"),
        )

        # Check the two parsers agree on the essentials
        self.assertEqual(str(seq_record.seq), record.sequence)
        self.assertEqual(seq_record.description, record.description)
        self.assertEqual(seq_record.name, record.entry_name)
        self.assertIn(seq_record.id, record.accessions)

        # Now try using the iterator - note that all these
        # test cases have only one record.

        # With the SequenceParser
        with open(datafile) as test_handle:
            records = list(SeqIO.parse(test_handle, "swiss"))

        self.assertEqual(len(records), 1)
        self.assertIsInstance(records[0], SeqRecord)

        # Check matches what we got earlier without the iterator:
        self.assertEqual(str(records[0].seq), str(seq_record.seq))
        self.assertEqual(records[0].description, seq_record.description)
        self.assertEqual(records[0].name, seq_record.name)
        self.assertEqual(records[0].id, seq_record.id)

        # With the RecordParser
        with open(datafile) as test_handle:
            records = list(SwissProt.parse(test_handle))

        self.assertEqual(len(records), 1)
        self.assertIsInstance(records[0], SwissProt.Record)

        # Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].sequence, record.sequence)
        self.assertEqual(records[0].description, record.description)
        self.assertEqual(records[0].entry_name, record.entry_name)
        self.assertEqual(records[0].accessions, record.accessions)

    def test_sp005(self):
        """Parsing SwissProt file sp005."""
        filename = "sp005"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)

        with open(datafile) as test_handle:
            seq_record = SeqIO.read(test_handle, "swiss")

        self.assertIsInstance(seq_record, SeqRecord)

        self.assertEqual(seq_record.id, "P24973")
        self.assertEqual(seq_record.name, "NU3M_BALPH")
        self.assertEqual(
            seq_record.description,
            "NADH-UBIQUINONE OXIDOREDUCTASE CHAIN 3 (EC 1.6.5.3).",
        )
        self.assertEqual(
            repr(seq_record.seq),
            "Seq('MNLLLTLLTNTTLALLLVFIAFWLPQLNVYAEKTSPYECGFDPMGSARLPFSMK...WAE')",
        )

        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(record.entry_name, "NU3M_BALPH")
        self.assertEqual(record.accessions, ["P24973"])
        self.assertEqual(
            record.organism_classification,
            [
                "Eukaryota",
                "Metazoa",
                "Chordata",
                "Craniata",
                "Vertebrata",
                "Mammalia",
                "Eutheria",
                "Cetartiodactyla",
                "Cetacea",
                "Mysticeti",
                "Balaenopteridae",
                "Balaenoptera",
            ],
        )
        self.assertEqual(record.seqinfo, (115, 13022, "ACF02965"))

        self.assertEqual(len(record.features), 0)

        self.assertEqual(len(record.references), 2)
        self.assertEqual(
            record.references[0].authors, "ARNASON U., GULLBERG A., WIDEGREN B."
        )
        self.assertEqual(
            record.references[0].title,
            "The complete nucleotide sequence of the mitochondrial DNA of the fin whale, Balaenoptera physalus.",
        )
        self.assertEqual(len(record.references[0].references), 1)
        self.assertEqual(record.references[0].references[0], ("MEDLINE", "92139449"))
        self.assertEqual(record.references[1].authors, "ARNASON U., GULLBERG A.")
        self.assertEqual(
            record.references[1].title,
            "Comparison between the complete mtDNA sequences of the blue and the fin whale, two species that can hybridize in nature.",
        )
        self.assertEqual(len(record.references[1].references), 1)
        self.assertEqual(record.references[1].references[0], ("MEDLINE", "94141932"))

        # Check the two parsers agree on the essentials
        self.assertEqual(str(seq_record.seq), record.sequence)
        self.assertEqual(seq_record.description, record.description)
        self.assertEqual(seq_record.name, record.entry_name)
        self.assertIn(seq_record.id, record.accessions)

        # Now try using the iterator - note that all these
        # test cases have only one record.

        # With the SequenceParser
        with open(datafile) as test_handle:
            records = list(SeqIO.parse(test_handle, "swiss"))

        self.assertEqual(len(records), 1)
        self.assertIsInstance(records[0], SeqRecord)

        # Check matches what we got earlier without the iterator:
        self.assertEqual(str(records[0].seq), str(seq_record.seq))
        self.assertEqual(records[0].description, seq_record.description)
        self.assertEqual(records[0].name, seq_record.name)
        self.assertEqual(records[0].id, seq_record.id)

        # With the RecordParser
        with open(datafile) as test_handle:
            records = list(SwissProt.parse(test_handle))

        self.assertEqual(len(records), 1)
        self.assertIsInstance(records[0], SwissProt.Record)

        # Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].sequence, record.sequence)
        self.assertEqual(records[0].description, record.description)
        self.assertEqual(records[0].entry_name, record.entry_name)
        self.assertEqual(records[0].accessions, record.accessions)

    def test_sp006(self):
        """Parsing SwissProt file sp006."""
        filename = "sp006"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)

        with open(datafile) as test_handle:
            seq_record = SeqIO.read(test_handle, "swiss")

        self.assertIsInstance(seq_record, SeqRecord)

        self.assertEqual(seq_record.id, "P39896")
        self.assertEqual(seq_record.name, "TCMO_STRGA")
        self.assertEqual(
            seq_record.description,
            "TETRACENOMYCIN POLYKETIDE SYNTHESIS 8-O-METHYL TRANSFERASE TCMO (EC 2.1.1.-).",
        )
        self.assertEqual(
            repr(seq_record.seq),
            "Seq('MTPHTHVRGPGDILQLTMAFYGSRALISAVELDLFTLLAGKPLPLGELCERAGI...KPR')",
        )

        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(record.entry_name, "TCMO_STRGA")
        self.assertEqual(record.accessions, ["P39896"])
        self.assertEqual(
            record.organism_classification,
            [
                "BACTERIA",
                "FIRMICUTES",
                "ACTINOBACTERIA",
                "ACTINOBACTERIDAE",
                "ACTINOMYCETALES",
                "STREPTOMYCINEAE",
                "STREPTOMYCETACEAE",
                "STREPTOMYCES",
            ],
        )
        self.assertEqual(record.seqinfo, (339, 37035, "848B7337"))

        self.assertEqual(len(record.features), 0)

        self.assertEqual(len(record.references), 1)
        self.assertEqual(
            record.references[0].authors,
            "SUMMERS R.G., WENDT-PIENKOWSKI E., MOTAMEDI H., HUTCHINSON C.R.",
        )
        self.assertEqual(
            record.references[0].title,
            "Nucleotide sequence of the tcmII-tcmIV region of the tetracenomycin C biosynthetic gene cluster of Streptomyces glaucescens and evidence that the tcmN gene encodes a multifunctional cyclase-dehydratase-O-methyl transferase.",
        )
        self.assertEqual(len(record.references[0].references), 1)
        self.assertEqual(record.references[0].references[0], ("MEDLINE", "92193265"))

        # Check the two parsers agree on the essentials
        self.assertEqual(str(seq_record.seq), record.sequence)
        self.assertEqual(seq_record.description, record.description)
        self.assertEqual(seq_record.name, record.entry_name)
        self.assertIn(seq_record.id, record.accessions)

        # Now try using the iterator - note that all these
        # test cases have only one record.

        # With the SequenceParser
        with open(datafile) as test_handle:
            records = list(SeqIO.parse(test_handle, "swiss"))

        self.assertEqual(len(records), 1)
        self.assertIsInstance(records[0], SeqRecord)

        # Check matches what we got earlier without the iterator:
        self.assertEqual(str(records[0].seq), str(seq_record.seq))
        self.assertEqual(records[0].description, seq_record.description)
        self.assertEqual(records[0].name, seq_record.name)
        self.assertEqual(records[0].id, seq_record.id)

        # With the RecordParser
        with open(datafile) as test_handle:
            records = list(SwissProt.parse(test_handle))

        self.assertEqual(len(records), 1)
        self.assertIsInstance(records[0], SwissProt.Record)

        # Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].sequence, record.sequence)
        self.assertEqual(records[0].description, record.description)
        self.assertEqual(records[0].entry_name, record.entry_name)
        self.assertEqual(records[0].accessions, record.accessions)

    def test_sp007(self):
        """Parsing SwissProt file sp007."""
        filename = "sp007"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)

        with open(datafile) as test_handle:
            seq_record = SeqIO.read(test_handle, "swiss")

        self.assertIsInstance(seq_record, SeqRecord)

        self.assertEqual(seq_record.id, "O95832")
        self.assertEqual(seq_record.name, "CLD1_HUMAN")
        self.assertEqual(
            seq_record.description,
            "CLAUDIN-1 (SENESCENCE-ASSOCIATED EPITHELIAL MEMBRANE PROTEIN).",
        )
        self.assertEqual(
            repr(seq_record.seq),
            "Seq('MANAGLQLLGFILAFLGWIGAIVSTALPQWRIYSYAGDNIVTAQAMYEGLWMSC...DYV')",
        )

        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(record.entry_name, "CLD1_HUMAN")
        self.assertEqual(record.accessions, ["O95832"])
        self.assertEqual(
            record.organism_classification,
            [
                "Eukaryota",
                "Metazoa",
                "Chordata",
                "Craniata",
                "Vertebrata",
                "Mammalia",
                "Eutheria",
                "Primates",
                "Catarrhini",
                "Hominidae",
                "Homo",
            ],
        )
        self.assertEqual(record.seqinfo, (211, 22744, "07269000E6C214F0"))

        self.assertEqual(len(record.features), 6)
        feature = record.features[0]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 7)
        self.assertEqual(feature.location.end, 28)
        self.assertEqual(feature.qualifiers["description"], "POTENTIAL.")
        self.assertIsNone(feature.id)
        feature = record.features[1]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 81)
        self.assertEqual(feature.location.end, 102)
        self.assertEqual(feature.qualifiers["description"], "POTENTIAL.")
        self.assertIsNone(feature.id)
        feature = record.features[2]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 115)
        self.assertEqual(feature.location.end, 136)
        self.assertEqual(feature.qualifiers["description"], "POTENTIAL.")
        self.assertIsNone(feature.id)
        feature = record.features[3]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 163)
        self.assertEqual(feature.location.end, 184)
        self.assertEqual(feature.qualifiers["description"], "POTENTIAL.")
        self.assertIsNone(feature.id)
        feature = record.features[4]
        self.assertEqual(feature.type, "CONFLICT")
        self.assertEqual(feature.location.start, 61)
        self.assertEqual(feature.location.end, 62)
        self.assertEqual(feature.qualifiers["description"], "I -> V (IN REF. 2).")
        self.assertIsNone(feature.id)
        feature = record.features[5]
        self.assertEqual(feature.type, "CONFLICT")
        self.assertEqual(feature.location.start, 134)
        self.assertEqual(feature.location.end, 135)
        self.assertEqual(feature.qualifiers["description"], "V -> A (IN REF. 2).")
        self.assertIsNone(feature.id)

        self.assertEqual(len(record.references), 2)
        self.assertEqual(
            record.references[0].authors,
            "Swisshelm K.L., Machl A., Planitzer S., Robertson R., Kubbies M., Hosier S.",
        )
        self.assertEqual(
            record.references[0].title,
            "SEMP1, a senescence-associated cDNA isolated from human mammary epithelial cells, is a member of an epithelial membrane protein superfamily.",
        )
        self.assertEqual(len(record.references[0].references), 1)
        self.assertEqual(record.references[0].references[0], ("MEDLINE", "99132301"))
        self.assertEqual(record.references[1].authors, "Mitic L.M., Anderson J.M.")
        self.assertEqual(
            record.references[1].title, "Human claudin-1 isolated from Caco-2 mRNA."
        )
        self.assertEqual(len(record.references[1].references), 0)

        # Check the two parsers agree on the essentials
        self.assertEqual(str(seq_record.seq), record.sequence)
        self.assertEqual(seq_record.description, record.description)
        self.assertEqual(seq_record.name, record.entry_name)
        self.assertIn(seq_record.id, record.accessions)

        # Now try using the iterator - note that all these
        # test cases have only one record.

        # With the SequenceParser
        with open(datafile) as test_handle:
            records = list(SeqIO.parse(test_handle, "swiss"))

        self.assertEqual(len(records), 1)
        self.assertIsInstance(records[0], SeqRecord)

        # Check matches what we got earlier without the iterator:
        self.assertEqual(str(records[0].seq), str(seq_record.seq))
        self.assertEqual(records[0].description, seq_record.description)
        self.assertEqual(records[0].name, seq_record.name)
        self.assertEqual(records[0].id, seq_record.id)

        # With the RecordParser
        with open(datafile) as test_handle:
            records = list(SwissProt.parse(test_handle))

        self.assertEqual(len(records), 1)
        self.assertIsInstance(records[0], SwissProt.Record)

        # Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].sequence, record.sequence)
        self.assertEqual(records[0].description, record.description)
        self.assertEqual(records[0].entry_name, record.entry_name)
        self.assertEqual(records[0].accessions, record.accessions)

    def test_sp008(self):
        """Parsing SwissProt file sp008."""
        filename = "sp008"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)

        with open(datafile) as test_handle:
            seq_record = SeqIO.read(test_handle, "swiss")

        self.assertIsInstance(seq_record, SeqRecord)

        self.assertEqual(seq_record.id, "P01892")
        self.assertEqual(seq_record.name, "1A02_HUMAN")
        self.assertEqual(
            seq_record.description,
            "HLA CLASS I HISTOCOMPATIBILITY ANTIGEN, A-2 ALPHA CHAIN PRECURSOR.",
        )
        self.assertEqual(
            repr(seq_record.seq),
            "Seq('MAVMAPRTLVLLLSGALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGYVDD...CKV')",
        )

        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(record.entry_name, "1A02_HUMAN")
        self.assertEqual(
            record.accessions,
            [
                "P01892",
                "P06338",
                "P30514",
                "P30444",
                "P30445",
                "P30446",
                "Q29680",
                "Q29899",
                "Q95352",
                "Q29837",
                "Q95380",
            ],
        )
        self.assertEqual(
            record.organism_classification,
            [
                "Eukaryota",
                "Metazoa",
                "Chordata",
                "Craniata",
                "Vertebrata",
                "Mammalia",
                "Eutheria",
                "Primates",
                "Catarrhini",
                "Hominidae",
                "Homo",
            ],
        )
        self.assertEqual(record.seqinfo, (365, 40922, "B54A97B24B337C08"))

        self.assertEqual(len(record.features), 71)
        feature = record.features[0]
        self.assertEqual(feature.type, "SIGNAL")
        self.assertEqual(feature.location.start, 0)
        self.assertEqual(feature.location.end, 24)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[1]
        self.assertEqual(feature.type, "CHAIN")
        self.assertEqual(feature.location.start, 24)
        self.assertEqual(feature.location.end, 365)
        self.assertEqual(
            feature.qualifiers["description"],
            "HLA CLASS I HISTOCOMPATIBILITY ANTIGEN, A-2 ALPHA CHAIN.",
        )
        self.assertIsNone(feature.id)
        feature = record.features[2]
        self.assertEqual(feature.type, "DOMAIN")
        self.assertEqual(feature.location.start, 24)
        self.assertEqual(feature.location.end, 114)
        self.assertEqual(feature.qualifiers["description"], "EXTRACELLULAR ALPHA-1.")
        self.assertIsNone(feature.id)
        feature = record.features[3]
        self.assertEqual(feature.type, "DOMAIN")
        self.assertEqual(feature.location.start, 114)
        self.assertEqual(feature.location.end, 206)
        self.assertEqual(feature.qualifiers["description"], "EXTRACELLULAR ALPHA-2.")
        self.assertIsNone(feature.id)
        feature = record.features[4]
        self.assertEqual(feature.type, "DOMAIN")
        self.assertEqual(feature.location.start, 206)
        self.assertEqual(feature.location.end, 298)
        self.assertEqual(feature.qualifiers["description"], "EXTRACELLULAR ALPHA-3.")
        self.assertIsNone(feature.id)
        feature = record.features[5]
        self.assertEqual(feature.type, "DOMAIN")
        self.assertEqual(feature.location.start, 298)
        self.assertEqual(feature.location.end, 308)
        self.assertEqual(feature.qualifiers["description"], "CONNECTING PEPTIDE.")
        self.assertIsNone(feature.id)
        feature = record.features[6]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 308)
        self.assertEqual(feature.location.end, 332)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[7]
        self.assertEqual(feature.type, "DOMAIN")
        self.assertEqual(feature.location.start, 332)
        self.assertEqual(feature.location.end, 365)
        self.assertEqual(feature.qualifiers["description"], "CYTOPLASMIC TAIL.")
        self.assertIsNone(feature.id)
        feature = record.features[8]
        self.assertEqual(feature.type, "CARBOHYD")
        self.assertEqual(feature.location.start, 109)
        self.assertEqual(feature.location.end, 110)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[9]
        self.assertEqual(feature.type, "DISULFID")
        self.assertEqual(feature.location.start, 124)
        self.assertEqual(feature.location.end, 188)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[10]
        self.assertEqual(feature.type, "DISULFID")
        self.assertEqual(feature.location.start, 226)
        self.assertEqual(feature.location.end, 283)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[11]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 26)
        self.assertEqual(feature.location.end, 36)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[12]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 44)
        self.assertEqual(feature.location.end, 52)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[13]
        self.assertEqual(feature.type, "TURN")
        self.assertEqual(feature.location.start, 52)
        self.assertEqual(feature.location.end, 54)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[14]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 54)
        self.assertEqual(feature.location.end, 61)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[15]
        self.assertEqual(feature.type, "TURN")
        self.assertEqual(feature.location.start, 61)
        self.assertEqual(feature.location.end, 63)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[16]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 69)
        self.assertEqual(feature.location.end, 71)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[17]
        self.assertEqual(feature.type, "HELIX")
        self.assertEqual(feature.location.start, 73)
        self.assertEqual(feature.location.end, 76)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[18]
        self.assertEqual(feature.type, "TURN")
        self.assertEqual(feature.location.start, 76)
        self.assertEqual(feature.location.end, 78)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[19]
        self.assertEqual(feature.type, "HELIX")
        self.assertEqual(feature.location.start, 80)
        self.assertEqual(feature.location.end, 108)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[20]
        self.assertEqual(feature.type, "TURN")
        self.assertEqual(feature.location.start, 108)
        self.assertEqual(feature.location.end, 110)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[21]
        self.assertEqual(feature.type, "TURN")
        self.assertEqual(feature.location.start, 112)
        self.assertEqual(feature.location.end, 114)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[22]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 117)
        self.assertEqual(feature.location.end, 127)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[23]
        self.assertEqual(feature.type, "TURN")
        self.assertEqual(feature.location.start, 128)
        self.assertEqual(feature.location.end, 130)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[24]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 132)
        self.assertEqual(feature.location.end, 142)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[25]
        self.assertEqual(feature.type, "TURN")
        self.assertEqual(feature.location.start, 142)
        self.assertEqual(feature.location.end, 144)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[26]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 144)
        self.assertEqual(feature.location.end, 150)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[27]
        self.assertEqual(feature.type, "TURN")
        self.assertEqual(feature.location.start, 151)
        self.assertEqual(feature.location.end, 153)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[28]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 156)
        self.assertEqual(feature.location.end, 159)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[29]
        self.assertEqual(feature.type, "TURN")
        self.assertEqual(feature.location.start, 162)
        self.assertEqual(feature.location.end, 163)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[30]
        self.assertEqual(feature.type, "HELIX")
        self.assertEqual(feature.location.start, 163)
        self.assertEqual(feature.location.end, 173)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[31]
        self.assertEqual(feature.type, "TURN")
        self.assertEqual(feature.location.start, 173)
        self.assertEqual(feature.location.end, 175)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[32]
        self.assertEqual(feature.type, "HELIX")
        self.assertEqual(feature.location.start, 175)
        self.assertEqual(feature.location.end, 185)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[33]
        self.assertEqual(feature.type, "TURN")
        self.assertEqual(feature.location.start, 185)
        self.assertEqual(feature.location.end, 186)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[34]
        self.assertEqual(feature.type, "HELIX")
        self.assertEqual(feature.location.start, 186)
        self.assertEqual(feature.location.end, 198)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[35]
        self.assertEqual(feature.type, "TURN")
        self.assertEqual(feature.location.start, 198)
        self.assertEqual(feature.location.end, 199)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[36]
        self.assertEqual(feature.type, "HELIX")
        self.assertEqual(feature.location.start, 199)
        self.assertEqual(feature.location.end, 203)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[37]
        self.assertEqual(feature.type, "TURN")
        self.assertEqual(feature.location.start, 203)
        self.assertEqual(feature.location.end, 204)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[38]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 206)
        self.assertEqual(feature.location.end, 207)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[39]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 209)
        self.assertEqual(feature.location.end, 219)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[40]
        self.assertEqual(feature.type, "TURN")
        self.assertEqual(feature.location.start, 219)
        self.assertEqual(feature.location.end, 221)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[41]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 221)
        self.assertEqual(feature.location.end, 233)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[42]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 237)
        self.assertEqual(feature.location.end, 243)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[43]
        self.assertEqual(feature.type, "TURN")
        self.assertEqual(feature.location.start, 243)
        self.assertEqual(feature.location.end, 245)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[44]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 245)
        self.assertEqual(feature.location.end, 247)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[45]
        self.assertEqual(feature.type, "HELIX")
        self.assertEqual(feature.location.start, 248)
        self.assertEqual(feature.location.end, 251)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[46]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 252)
        self.assertEqual(feature.location.end, 254)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[47]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 257)
        self.assertEqual(feature.location.end, 259)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[48]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 264)
        self.assertEqual(feature.location.end, 274)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[49]
        self.assertEqual(feature.type, "TURN")
        self.assertEqual(feature.location.start, 274)
        self.assertEqual(feature.location.end, 276)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[50]
        self.assertEqual(feature.type, "HELIX")
        self.assertEqual(feature.location.start, 277)
        self.assertEqual(feature.location.end, 280)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[51]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 280)
        self.assertEqual(feature.location.end, 286)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[52]
        self.assertEqual(feature.type, "TURN")
        self.assertEqual(feature.location.start, 287)
        self.assertEqual(feature.location.end, 289)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[53]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 293)
        self.assertEqual(feature.location.end, 297)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[54]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 32)
        self.assertEqual(feature.location.end, 33)
        self.assertEqual(
            feature.qualifiers["description"],
            "F -> Y (IN A*0205, A*0206, A*0208, A*0210 AND A*0221).",
        )
        self.assertEqual(feature.id, "VAR_004334")
        feature = record.features[55]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 53)
        self.assertEqual(feature.location.end, 54)
        self.assertEqual(feature.qualifiers["description"], "D -> N (IN A*0221).")
        self.assertEqual(feature.id, "VAR_004335")
        feature = record.features[56]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 66)
        self.assertEqual(feature.location.end, 67)
        self.assertEqual(
            feature.qualifiers["description"], "Q -> R (IN A*0202, A*0205, AND A*0208)."
        )
        self.assertEqual(feature.id, "VAR_004336")
        feature = record.features[57]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 89)
        self.assertEqual(feature.location.end, 90)
        self.assertEqual(
            feature.qualifiers["description"], "K -> N (IN A*0208 AND A*0220)."
        )
        self.assertEqual(feature.id, "VAR_004337")
        feature = record.features[58]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 96)
        self.assertEqual(feature.location.end, 98)
        self.assertEqual(feature.qualifiers["description"], "TH -> ID (IN A*0211).")
        self.assertEqual(feature.id, "VAR_004338")
        feature = record.features[59]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 118)
        self.assertEqual(feature.location.end, 119)
        self.assertEqual(
            feature.qualifiers["description"],
            "V -> L (IN A*0202, A*0205, A*0208 AND A*0217).",
        )
        self.assertEqual(feature.id, "VAR_004339")
        feature = record.features[60]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 120)
        self.assertEqual(feature.location.end, 121)
        self.assertEqual(
            feature.qualifiers["description"], "R -> M (IN A*0204 AND A*0217)."
        )
        self.assertEqual(feature.id, "VAR_004340")
        feature = record.features[61]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 122)
        self.assertEqual(feature.location.end, 123)
        self.assertEqual(
            feature.qualifiers["description"], "Y -> C (IN A*0207 AND A*0218)."
        )
        self.assertEqual(feature.id, "VAR_004341")
        feature = record.features[62]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 122)
        self.assertEqual(feature.location.end, 123)
        self.assertEqual(
            feature.qualifiers["description"], "Y -> F (IN A*0210 AND A*0217)."
        )
        self.assertEqual(feature.id, "VAR_004342")
        feature = record.features[63]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 130)
        self.assertEqual(feature.location.end, 131)
        self.assertEqual(feature.qualifiers["description"], "W -> G (IN A*0210).")
        self.assertEqual(feature.id, "VAR_004343")
        feature = record.features[64]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 161)
        self.assertEqual(feature.location.end, 162)
        self.assertEqual(feature.qualifiers["description"], "M -> K (IN A*0218).")
        self.assertEqual(feature.id, "VAR_004344")
        feature = record.features[65]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 172)
        self.assertEqual(feature.location.end, 173)
        self.assertEqual(feature.qualifiers["description"], "A -> T (IN A*0203).")
        self.assertEqual(feature.id, "VAR_004345")
        feature = record.features[66]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 175)
        self.assertEqual(feature.location.end, 176)
        self.assertEqual(
            feature.qualifiers["description"], "V -> E (IN A*0203 AND A*0213)."
        )
        self.assertEqual(feature.id, "VAR_004346")
        feature = record.features[67]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 179)
        self.assertEqual(feature.location.end, 180)
        self.assertEqual(
            feature.qualifiers["description"],
            "L -> W (IN A*0202, A*0203, A*0205 AND A*0208).",
        )
        self.assertEqual(feature.id, "VAR_004347")
        feature = record.features[68]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 179)
        self.assertEqual(feature.location.end, 180)
        self.assertEqual(
            feature.qualifiers["description"], "L -> Q (IN A*0212 AND A*0213)."
        )
        self.assertEqual(feature.id, "VAR_004348")
        feature = record.features[69]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 186)
        self.assertEqual(feature.location.end, 187)
        self.assertEqual(feature.qualifiers["description"], "T -> E (IN A*0216).")
        self.assertEqual(feature.id, "VAR_004349")
        feature = record.features[70]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 259)
        self.assertEqual(feature.location.end, 260)
        self.assertEqual(feature.qualifiers["description"], "A -> E (IN A*0209).")
        self.assertEqual(feature.id, "VAR_004350")

        self.assertEqual(len(record.references), 27)
        self.assertEqual(record.references[0].authors, "Koller B.H., Orr H.T.")
        self.assertEqual(
            record.references[0].title,
            "Cloning and complete sequence of an HLA-A2 gene: analysis of two HLA-A alleles at the nucleotide level.",
        )
        self.assertEqual(len(record.references[0].references), 1)
        self.assertEqual(record.references[0].references[0], ("MEDLINE", "85132727"))
        self.assertEqual(
            record.references[1].authors,
            "Cianetti L., Testa U., Scotto L., la Valle R., Simeone A., Boccoli G., Giannella G., Peschle C., Boncinelli E.",
        )
        self.assertEqual(
            record.references[1].title,
            "Three new class I HLA alleles: structure of mRNAs and alternative mechanisms of processing.",
        )
        self.assertEqual(len(record.references[1].references), 1)
        self.assertEqual(record.references[1].references[0], ("MEDLINE", "89122144"))
        self.assertEqual(
            record.references[2].authors,
            "Ennis P.D., Zemmour J., Salter R.D., Parham P.",
        )
        self.assertEqual(
            record.references[2].title,
            "Rapid cloning of HLA-A,B cDNA by using the polymerase chain reaction: frequency and nature of errors produced in amplification.",
        )
        self.assertEqual(len(record.references[2].references), 1)
        self.assertEqual(record.references[2].references[0], ("MEDLINE", "90207291"))
        self.assertEqual(
            record.references[3].authors,
            "Belich M.P., Madrigal J.A., Hildebrand W.H., Zemmour J., Williams R.C., Luz R., Petzl-Erler M.L., Parham P.",
        )
        self.assertEqual(
            record.references[3].title,
            "Unusual HLA-B alleles in two tribes of Brazilian Indians.",
        )
        self.assertEqual(len(record.references[3].references), 1)
        self.assertEqual(record.references[3].references[0], ("MEDLINE", "92269955"))
        self.assertEqual(record.references[4].authors, "Krangel M.S.")
        self.assertEqual(
            record.references[4].title,
            "Unusual RNA splicing generates a secreted form of HLA-A2 in a mutagenized B lymphoblastoid cell line.",
        )
        self.assertEqual(len(record.references[4].references), 1)
        self.assertEqual(record.references[4].references[0], ("MEDLINE", "85230571"))
        self.assertEqual(
            record.references[5].authors,
            "Orr H.T., Lopez de Castro J.A., Parham P., Ploegh H.L., Strominger J.L.",
        )
        self.assertEqual(
            record.references[5].title,
            "Comparison of amino acid sequences of two human histocompatibility antigens, HLA-A2 and HLA-B7: location of putative alloantigenic sites.",
        )
        self.assertEqual(len(record.references[5].references), 1)
        self.assertEqual(record.references[5].references[0], ("MEDLINE", "80056745"))
        self.assertEqual(
            record.references[6].authors,
            "Lopez de Castro J.A., Strominger J.L., Strong D.M., Orr H.T.",
        )
        self.assertEqual(
            record.references[6].title,
            "Structure of crossreactive human histocompatibility antigens HLA-A28 and HLA-A2: possible implications for the generation of HLA polymorphism.",
        )
        self.assertEqual(len(record.references[6].references), 1)
        self.assertEqual(record.references[6].references[0], ("MEDLINE", "82247941"))
        self.assertEqual(
            record.references[7].authors,
            "Mattson D.H., Handy D.E., Bradley D.A., Coligan J.E., Cowan E.P., Biddison W.E.",
        )
        self.assertEqual(
            record.references[7].title,
            "DNA sequences of the genes that encode the CTL-defined HLA-A2 variants M7 and DK1.",
        )
        self.assertEqual(len(record.references[7].references), 1)
        self.assertEqual(record.references[7].references[0], ("MEDLINE", "87306734"))
        self.assertEqual(
            record.references[8].authors,
            "Holmes N., Ennis P., Wan A.M., Denney D.W., Parham P.",
        )
        self.assertEqual(
            record.references[8].title,
            "Multiple genetic mechanisms have contributed to the generation of the HLA-A2/A28 family of class I MHC molecules.",
        )
        self.assertEqual(len(record.references[8].references), 1)
        self.assertEqual(record.references[8].references[0], ("MEDLINE", "87252273"))
        self.assertEqual(record.references[9].authors, "Domena J.D.")
        self.assertEqual(record.references[9].title, "")
        self.assertEqual(len(record.references[9].references), 0)
        self.assertEqual(
            record.references[10].authors, "Castano A.R., Lopez de Castro J.A."
        )
        self.assertEqual(
            record.references[10].title,
            "Structure of the HLA-A*0204 antigen, found in South American Indians. Spatial clustering of HLA-A2 subtype polymorphism.",
        )
        self.assertEqual(len(record.references[10].references), 1)
        self.assertEqual(record.references[10].references[0], ("MEDLINE", "92039809"))
        self.assertEqual(
            record.references[11].authors,
            "Watkins D.I., McAdam S.N., Liu X., Stang C.R., Milford E.L., Levine C.G., Garber T.L., Dogon A.L., Lord C.I., Ghim S.H., Troup G.M., Hughes A.L., Letvin N.L.",
        )
        self.assertEqual(
            record.references[11].title,
            "New recombinant HLA-B alleles in a tribe of South American Amerindians indicate rapid evolution of MHC class I loci.",
        )
        self.assertEqual(len(record.references[11].references), 1)
        self.assertEqual(record.references[11].references[0], ("MEDLINE", "92269956"))
        self.assertEqual(
            record.references[12].authors,
            "Parham P., Lawlor D.A., Lomen C.E., Ennis P.D.",
        )
        self.assertEqual(
            record.references[12].title,
            "Diversity and diversification of HLA-A,B,C alleles.",
        )
        self.assertEqual(len(record.references[12].references), 1)
        self.assertEqual(record.references[12].references[0], ("MEDLINE", "89235215"))
        self.assertEqual(
            record.references[13].authors,
            "Ezquerra A., Domenech N., van der Poel J., Strominger J.L., Vega M.A., Lopez de Castro J.A.",
        )
        self.assertEqual(
            record.references[13].title,
            "Molecular analysis of an HLA-A2 functional variant CLA defined by cytolytic T lymphocytes.",
        )
        self.assertEqual(len(record.references[13].references), 1)
        self.assertEqual(record.references[13].references[0], ("MEDLINE", "86305811"))
        self.assertEqual(
            record.references[14].authors,
            "Domenech N., Ezquerra A., Castano R., Lopez de Castro J.A.",
        )
        self.assertEqual(
            record.references[14].title,
            "Structural analysis of HLA-A2.4 functional variant KNE. Implications for the mapping of HLA-A2-specific T-cell epitopes.",
        )
        self.assertEqual(len(record.references[14].references), 1)
        self.assertEqual(record.references[14].references[0], ("MEDLINE", "88113844"))
        self.assertEqual(
            record.references[15].authors,
            "Domenech N., Castano R., Goulmy E., Lopez de Castro J.A.",
        )
        self.assertEqual(
            record.references[15].title,
            "Molecular analysis of HLA-A2.4 functional variant KLO: close structural and evolutionary relatedness to the HLA-A2.2 subtype.",
        )
        self.assertEqual(len(record.references[15].references), 1)
        self.assertEqual(record.references[15].references[0], ("MEDLINE", "88314183"))
        self.assertEqual(
            record.references[16].authors,
            "Castano R., Ezquerra A., Domenech N., Lopez de Castro J.A.",
        )
        self.assertEqual(
            record.references[16].title,
            "An HLA-A2 population variant with structural polymorphism in the alpha 3 region.",
        )
        self.assertEqual(len(record.references[16].references), 1)
        self.assertEqual(record.references[16].references[0], ("MEDLINE", "88186100"))
        self.assertEqual(
            record.references[17].authors, "Epstein H., Kennedy L., Holmes N."
        )
        self.assertEqual(
            record.references[17].title,
            "An Oriental HLA-A2 subtype is closely related to a subset of Caucasoid HLA-A2 alleles.",
        )
        self.assertEqual(len(record.references[17].references), 1)
        self.assertEqual(record.references[17].references[0], ("MEDLINE", "89122133"))
        self.assertEqual(
            record.references[18].authors, "Castano A.R., Lopez de Castro J.A."
        )
        self.assertEqual(
            record.references[18].title,
            "Structure of the HLA-A*0211 (A2.5) subtype: further evidence for selection-driven diversification of HLA-A2 antigens.",
        )
        self.assertEqual(len(record.references[18].references), 1)
        self.assertEqual(record.references[18].references[0], ("MEDLINE", "92218010"))
        self.assertEqual(
            record.references[19].authors,
            "Barber D.F., Fernandez J.M., Lopez de Castro J.A.",
        )
        self.assertEqual(
            record.references[19].title,
            "Primary structure of a new HLA-A2 subtype: HLA-A*0213.",
        )
        self.assertEqual(len(record.references[19].references), 1)
        self.assertEqual(record.references[19].references[0], ("MEDLINE", "94222455"))
        self.assertEqual(
            record.references[20].authors,
            "Barouch D., Krausa P., Bodmer J., Browning M.J., McMichael A.J.",
        )
        self.assertEqual(
            record.references[20].title,
            "Identification of a novel HLA-A2 subtype, HLA-A*0216.",
        )
        self.assertEqual(len(record.references[20].references), 1)
        self.assertEqual(record.references[20].references[0], ("MEDLINE", "95278976"))
        self.assertEqual(
            record.references[21].authors,
            "Selvakumar A., Granja C.B., Salazar M., Alosco S.M., Yunis E.J., Dupont B.",
        )
        self.assertEqual(
            record.references[21].title,
            "A novel subtype of A2 (A*0217) isolated from the South American Indian B-cell line AMALA.",
        )
        self.assertEqual(len(record.references[21].references), 1)
        self.assertEqual(record.references[21].references[0], ("MEDLINE", "95381236"))
        self.assertEqual(
            record.references[22].authors,
            "Kashiwase K., Tokunaga K., Ishikawa Y., Oohashi H., Hashimoto M., Akaza T., Tadokoro K., Juji T.",
        )
        self.assertEqual(
            record.references[22].title, "A new A2 sequence HLA-A2K from Japanese."
        )
        self.assertEqual(len(record.references[22].references), 0)
        self.assertEqual(
            record.references[23].authors,
            "Fleischhauer K., Zino E., Mazzi B., Severini G.M., Benazzi E., Bordignon C.",
        )
        self.assertEqual(
            record.references[23].title,
            "HLA-A*02 subtype distribution in Caucasians from northern Italy: identification of A*0220.",
        )
        self.assertEqual(len(record.references[23].references), 1)
        self.assertEqual(record.references[23].references[0], ("MEDLINE", "97161038"))
        self.assertEqual(record.references[24].authors, "Szmania S., Baxter-Lowe L.A.")
        self.assertEqual(
            record.references[24].title, "Nucleotide sequence of a novel HLA-A2 gene."
        )
        self.assertEqual(len(record.references[24].references), 0)
        self.assertEqual(
            record.references[25].authors,
            "Bjorkman P.J., Saper M.A., Samraoui B., Bennett W.S., Strominger J.L., Wiley D.C.",
        )
        self.assertEqual(
            record.references[25].title,
            "Structure of the human class I histocompatibility antigen, HLA-A2.",
        )
        self.assertEqual(len(record.references[25].references), 1)
        self.assertEqual(record.references[25].references[0], ("MEDLINE", "88014204"))
        self.assertEqual(
            record.references[26].authors, "Saper M.A., Bjorkman P.J., Wiley D.C."
        )
        self.assertEqual(
            record.references[26].title,
            "Refined structure of the human histocompatibility antigen HLA-A2 at 2.6-A resolution.",
        )
        self.assertEqual(len(record.references[26].references), 1)
        self.assertEqual(record.references[26].references[0], ("MEDLINE", "91245570"))

        # Check the two parsers agree on the essentials
        self.assertEqual(str(seq_record.seq), record.sequence)
        self.assertEqual(seq_record.description, record.description)
        self.assertEqual(seq_record.name, record.entry_name)
        self.assertIn(seq_record.id, record.accessions)

        # Now try using the iterator - note that all these
        # test cases have only one record.

        # With the SequenceParser
        with open(datafile) as test_handle:
            records = list(SeqIO.parse(test_handle, "swiss"))

        self.assertEqual(len(records), 1)
        self.assertIsInstance(records[0], SeqRecord)

        # Check matches what we got earlier without the iterator:
        self.assertEqual(str(records[0].seq), str(seq_record.seq))
        self.assertEqual(records[0].description, seq_record.description)
        self.assertEqual(records[0].name, seq_record.name)
        self.assertEqual(records[0].id, seq_record.id)

        # With the RecordParser
        with open(datafile) as test_handle:
            records = list(SwissProt.parse(test_handle))

        self.assertEqual(len(records), 1)
        self.assertIsInstance(records[0], SwissProt.Record)

        # Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].sequence, record.sequence)
        self.assertEqual(records[0].description, record.description)
        self.assertEqual(records[0].entry_name, record.entry_name)
        self.assertEqual(records[0].accessions, record.accessions)

    def test_sp009(self):
        """Parsing SwissProt file sp009."""
        filename = "sp009"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)

        with open(datafile) as test_handle:
            seq_record = SeqIO.read(test_handle, "swiss")

        self.assertIsInstance(seq_record, SeqRecord)

        self.assertEqual(seq_record.id, "O23729")
        self.assertEqual(seq_record.name, "CHS3_BROFI")
        self.assertEqual(
            seq_record.description,
            "CHALCONE SYNTHASE 3 (EC 2.3.1.74) (NARINGENIN-CHALCONE SYNTHASE 3).",
        )
        self.assertEqual(
            repr(seq_record.seq),
            "Seq('MAPAMEEIRQAQRAEGPAAVLAIGTSTPPNALYQADYPDYYFRITKSEHLTELK...GAE')",
        )

        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(record.entry_name, "CHS3_BROFI")
        self.assertEqual(record.accessions, ["O23729"])
        self.assertEqual(
            record.organism_classification,
            [
                "Eukaryota",
                "Viridiplantae",
                "Embryophyta",
                "Tracheophyta",
                "Spermatophyta",
                "Magnoliophyta",
                "Liliopsida",
                "Asparagales",
                "Orchidaceae",
                "Bromheadia",
            ],
        )
        self.assertEqual(record.seqinfo, (394, 42941, "2F8D14AF4870BBB2"))

        self.assertEqual(len(record.features), 1)
        feature = record.features[0]
        self.assertEqual(feature.type, "ACT_SITE")
        self.assertEqual(feature.location.start, 164)
        self.assertEqual(feature.location.end, 165)
        self.assertEqual(feature.qualifiers["description"], "BY SIMILARITY.")
        self.assertIsNone(feature.id)

        self.assertEqual(len(record.references), 1)
        self.assertEqual(
            record.references[0].authors, "Liew C.F., Lim S.H., Loh C.S., Goh C.J."
        )
        self.assertEqual(
            record.references[0].title,
            "Molecular cloning and sequence analysis of chalcone synthase cDNAs of Bromheadia finlaysoniana.",
        )
        self.assertEqual(len(record.references[0].references), 0)

        # Check the two parsers agree on the essentials
        self.assertEqual(str(seq_record.seq), record.sequence)
        self.assertEqual(seq_record.description, record.description)
        self.assertEqual(seq_record.name, record.entry_name)
        self.assertIn(seq_record.id, record.accessions)

        # Now try using the iterator - note that all these
        # test cases have only one record.

        # With the SequenceParser
        with open(datafile) as test_handle:
            records = list(SeqIO.parse(test_handle, "swiss"))

        self.assertEqual(len(records), 1)
        self.assertIsInstance(records[0], SeqRecord)

        # Check matches what we got earlier without the iterator:
        self.assertEqual(str(records[0].seq), str(seq_record.seq))
        self.assertEqual(records[0].description, seq_record.description)
        self.assertEqual(records[0].name, seq_record.name)
        self.assertEqual(records[0].id, seq_record.id)

        # With the RecordParser
        with open(datafile) as test_handle:
            records = list(SwissProt.parse(test_handle))

        self.assertEqual(len(records), 1)
        self.assertIsInstance(records[0], SwissProt.Record)

        # Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].sequence, record.sequence)
        self.assertEqual(records[0].description, record.description)
        self.assertEqual(records[0].entry_name, record.entry_name)
        self.assertEqual(records[0].accessions, record.accessions)

    def test_sp010(self):
        """Parsing SwissProt file sp010."""
        filename = "sp010"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)

        with open(datafile) as test_handle:
            seq_record = SeqIO.read(test_handle, "swiss")

        self.assertIsInstance(seq_record, SeqRecord)

        self.assertEqual(seq_record.id, "Q13639")
        self.assertEqual(seq_record.name, "5H4_HUMAN")
        self.assertEqual(
            seq_record.description,
            "5-HYDROXYTRYPTAMINE 4 RECEPTOR (5-HT-4) (SEROTONIN RECEPTOR) (5-HT4).",
        )
        self.assertEqual(
            repr(seq_record.seq),
            "Seq('MDKLDANVSSEEGFGSVEKVVLLTFLSTVILMAILGNLLVMVAVCWDRQLRKIK...SDT')",
        )

        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(record.entry_name, "5H4_HUMAN")
        self.assertEqual(
            record.accessions,
            ["Q13639", "Q9UBM6", "Q9UQR6", "Q9UE22", "Q9UE23", "Q9UBT4", "Q9NY73"],
        )
        self.assertEqual(
            record.organism_classification,
            [
                "Eukaryota",
                "Metazoa",
                "Chordata",
                "Craniata",
                "Vertebrata",
                "Euteleostomi",
                "Mammalia",
                "Eutheria",
                "Primates",
                "Catarrhini",
                "Hominidae",
                "Homo",
            ],
        )
        self.assertEqual(record.seqinfo, (388, 43761, "7FCFEC60E7BDF560"))

        self.assertEqual(len(record.features), 23)
        feature = record.features[0]
        self.assertEqual(feature.type, "DOMAIN")
        self.assertEqual(feature.location.start, 0)
        self.assertEqual(feature.location.end, 19)
        self.assertEqual(
            feature.qualifiers["description"], "EXTRACELLULAR (POTENTIAL)."
        )
        self.assertIsNone(feature.id)
        feature = record.features[1]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 19)
        self.assertEqual(feature.location.end, 40)
        self.assertEqual(feature.qualifiers["description"], "1 (POTENTIAL).")
        self.assertIsNone(feature.id)
        feature = record.features[2]
        self.assertEqual(feature.type, "DOMAIN")
        self.assertEqual(feature.location.start, 40)
        self.assertEqual(feature.location.end, 58)
        self.assertEqual(feature.qualifiers["description"], "CYTOPLASMIC (POTENTIAL).")
        self.assertIsNone(feature.id)
        feature = record.features[3]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 58)
        self.assertEqual(feature.location.end, 79)
        self.assertEqual(feature.qualifiers["description"], "2 (POTENTIAL).")
        self.assertIsNone(feature.id)
        feature = record.features[4]
        self.assertEqual(feature.type, "DOMAIN")
        self.assertEqual(feature.location.start, 79)
        self.assertEqual(feature.location.end, 93)
        self.assertEqual(
            feature.qualifiers["description"], "EXTRACELLULAR (POTENTIAL)."
        )
        self.assertIsNone(feature.id)
        feature = record.features[5]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 93)
        self.assertEqual(feature.location.end, 116)
        self.assertEqual(feature.qualifiers["description"], "3 (POTENTIAL).")
        self.assertIsNone(feature.id)
        feature = record.features[6]
        self.assertEqual(feature.type, "DOMAIN")
        self.assertEqual(feature.location.start, 116)
        self.assertEqual(feature.location.end, 137)
        self.assertEqual(feature.qualifiers["description"], "CYTOPLASMIC (POTENTIAL).")
        self.assertIsNone(feature.id)
        feature = record.features[7]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 137)
        self.assertEqual(feature.location.end, 158)
        self.assertEqual(feature.qualifiers["description"], "4 (POTENTIAL).")
        self.assertIsNone(feature.id)
        feature = record.features[8]
        self.assertEqual(feature.type, "DOMAIN")
        self.assertEqual(feature.location.start, 158)
        self.assertEqual(feature.location.end, 192)
        self.assertEqual(
            feature.qualifiers["description"], "EXTRACELLULAR (POTENTIAL)."
        )
        self.assertIsNone(feature.id)
        feature = record.features[9]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 192)
        self.assertEqual(feature.location.end, 213)
        self.assertEqual(feature.qualifiers["description"], "5 (POTENTIAL).")
        self.assertIsNone(feature.id)
        feature = record.features[10]
        self.assertEqual(feature.type, "DOMAIN")
        self.assertEqual(feature.location.start, 213)
        self.assertEqual(feature.location.end, 260)
        self.assertEqual(feature.qualifiers["description"], "CYTOPLASMIC (POTENTIAL).")
        self.assertIsNone(feature.id)
        feature = record.features[11]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 260)
        self.assertEqual(feature.location.end, 281)
        self.assertEqual(feature.qualifiers["description"], "6 (POTENTIAL).")
        self.assertIsNone(feature.id)
        feature = record.features[12]
        self.assertEqual(feature.type, "DOMAIN")
        self.assertEqual(feature.location.start, 281)
        self.assertEqual(feature.location.end, 294)
        self.assertEqual(
            feature.qualifiers["description"], "EXTRACELLULAR (POTENTIAL)."
        )
        self.assertIsNone(feature.id)
        feature = record.features[13]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 294)
        self.assertEqual(feature.location.end, 315)
        self.assertEqual(feature.qualifiers["description"], "7 (POTENTIAL).")
        self.assertIsNone(feature.id)
        feature = record.features[14]
        self.assertEqual(feature.type, "DOMAIN")
        self.assertEqual(feature.location.start, 315)
        self.assertEqual(feature.location.end, 388)
        self.assertEqual(feature.qualifiers["description"], "CYTOPLASMIC (POTENTIAL).")
        self.assertIsNone(feature.id)
        feature = record.features[15]
        self.assertEqual(feature.type, "CARBOHYD")
        self.assertEqual(feature.location.start, 6)
        self.assertEqual(feature.location.end, 7)
        self.assertEqual(
            feature.qualifiers["description"], "N-LINKED (GLCNAC...) (POTENTIAL)."
        )
        self.assertIsNone(feature.id)
        feature = record.features[16]
        self.assertEqual(feature.type, "DISULFID")
        self.assertEqual(feature.location.start, 92)
        self.assertEqual(feature.location.end, 184)
        self.assertEqual(feature.qualifiers["description"], "BY SIMILARITY.")
        self.assertIsNone(feature.id)
        feature = record.features[17]
        self.assertEqual(feature.type, "LIPID")
        self.assertEqual(feature.location.start, 328)
        self.assertEqual(feature.location.end, 329)
        self.assertEqual(
            feature.qualifiers["description"], "PALMITATE (BY SIMILARITY)."
        )
        self.assertIsNone(feature.id)
        feature = record.features[18]
        self.assertEqual(feature.type, "VARSPLIC")
        self.assertEqual(feature.location.start, 168)
        self.assertEqual(feature.location.end, 169)
        self.assertEqual(
            feature.qualifiers["description"],
            "L -> LERSLNQGLGQDFHA (IN ISOFORM 5-HT4(F)).",
        )
        self.assertIsNone(feature.id)
        feature = record.features[19]
        self.assertEqual(feature.type, "VARSPLIC")
        self.assertEqual(feature.location.start, 358)
        self.assertEqual(feature.location.end, 388)
        self.assertEqual(
            feature.qualifiers["description"],
            "RDAVECGGQWESQCHPPATSPLVAAQPSDT -> SGCSPVSSFLLLFCNRPVPV (IN ISOFORM 5-HT4(E)).",
        )
        self.assertIsNone(feature.id)
        feature = record.features[20]
        self.assertEqual(feature.type, "VARSPLIC")
        self.assertEqual(feature.location.start, 358)
        self.assertEqual(feature.location.end, 388)
        self.assertEqual(
            feature.qualifiers["description"],
            "RDAVECGGQWESQCHPPATSPLVAAQPSDT -> SSGTETDRRNFGIRKRRLTKPS (IN ISOFORM 5-HT4(D)).",
        )
        self.assertIsNone(feature.id)
        feature = record.features[21]
        self.assertEqual(feature.type, "VARSPLIC")
        self.assertEqual(feature.location.start, 359)
        self.assertEqual(feature.location.end, 388)
        self.assertEqual(
            feature.qualifiers["description"],
            "DAVECGGQWESQCHPPATSPLVAAQPSDT -> F (IN ISOFORM 5-HT4(C)).",
        )
        self.assertIsNone(feature.id)
        feature = record.features[22]
        self.assertEqual(feature.type, "VARSPLIC")
        self.assertEqual(feature.location.start, 359)
        self.assertEqual(feature.location.end, 388)
        self.assertEqual(
            feature.qualifiers["description"],
            "DAVECGGQWESQCHPPATSPLVAAQPSDT -> YTVLHRGHHQELEKLPIHNDPESLESCF (IN ISOFORM 5-HT4(A)).",
        )
        self.assertIsNone(feature.id)
        self.assertEqual(len(record.references), 6)

        self.assertEqual(
            record.references[0].authors,
            "Blondel O., Gastineau M., Dahmoune Y., Langlois M., Fischmeister R.",
        )
        self.assertEqual(
            record.references[0].title,
            "Cloning, expression, and pharmacology of four human 5-hydroxytryptamine receptor isoforms produced by alternative splicing in the carboxyl terminus.",
        )
        self.assertEqual(len(record.references[0].references), 1)
        self.assertEqual(record.references[0].references[0], ("PubMed", "9603189"))
        self.assertEqual(
            record.references[1].authors,
            "Van den Wyngaert I., Gommeren W., Jurzak M., Verhasselt P., Gordon R., Leysen J., Luyten W., Bender E.",
        )
        self.assertEqual(
            record.references[1].title,
            "Cloning and expression of 5-HT4 receptor species and splice variants.",
        )
        self.assertEqual(len(record.references[1].references), 0)
        self.assertEqual(
            record.references[2].authors,
            "Claeysen S., Faye P., Sebben M., Lemaire S., Bockaert J., Dumuis A.",
        )
        self.assertEqual(
            record.references[2].title,
            "Cloning and expression of human 5-HT4S receptors. Effect of receptor density on their coupling to adenylyl cyclase.",
        )
        self.assertEqual(len(record.references[2].references), 1)
        self.assertEqual(record.references[2].references[0], ("PubMed", "9351641"))
        self.assertEqual(
            record.references[3].authors,
            "Claeysen S., Sebben M., Becamel C., Bockaert J., Dumuis A.",
        )
        self.assertEqual(
            record.references[3].title,
            "Novel brain-specific 5-HT4 receptors splice variants show marked constitutive activity: role of the c-terminal intracellular domain.",
        )
        self.assertEqual(len(record.references[3].references), 0)
        self.assertEqual(
            record.references[4].authors,
            "Bender E., Pindon A., van Oers I., Zhang Y.B., Gommeren W., Verhasselt P., Jurzak M., Leysen J., Luyten W.",
        )
        self.assertEqual(
            record.references[4].title,
            "Structure of the human serotonin 5-HT4 receptor gene and cloning of a novel 5-HT4 splice variant.",
        )
        self.assertEqual(len(record.references[4].references), 1)
        self.assertEqual(record.references[4].references[0], ("PubMed", "10646498"))
        self.assertEqual(
            record.references[5].authors,
            "Ullmer C., Schmuck K., Kalkman H.O., Lubbert H.",
        )
        self.assertEqual(
            record.references[5].title,
            "Expression of serotonin receptor mRNAs in blood vessels.",
        )
        self.assertEqual(len(record.references[5].references), 2)
        self.assertEqual(record.references[5].references[0], ("MEDLINE", "95385798"))
        self.assertEqual(record.references[5].references[1], ("PubMed", "7656980"))

        # Check the two parsers agree on the essentials
        self.assertEqual(str(seq_record.seq), record.sequence)
        self.assertEqual(seq_record.description, record.description)
        self.assertEqual(seq_record.name, record.entry_name)
        self.assertIn(seq_record.id, record.accessions)

        # Now try using the iterator - note that all these
        # test cases have only one record.

        # With the SequenceParser
        records = list(SeqIO.parse(datafile, "swiss"))

        self.assertEqual(len(records), 1)
        self.assertIsInstance(records[0], SeqRecord)

        # Check matches what we got earlier without the iterator:
        self.assertEqual(str(records[0].seq), str(seq_record.seq))
        self.assertEqual(records[0].description, seq_record.description)
        self.assertEqual(records[0].name, seq_record.name)
        self.assertEqual(records[0].id, seq_record.id)

        # With the RecordParser
        with open(datafile) as test_handle:
            records = list(SwissProt.parse(test_handle))

        self.assertEqual(len(records), 1)
        self.assertIsInstance(records[0], SwissProt.Record)

        # Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].sequence, record.sequence)
        self.assertEqual(records[0].description, record.description)
        self.assertEqual(records[0].entry_name, record.entry_name)
        self.assertEqual(records[0].accessions, record.accessions)

    def test_Q13639(self):
        """Parsing SwissProt file Q13639."""
        filename = "Q13639.txt"
        # this is a more recent version of the file sp010 above

        datafile = os.path.join("SwissProt", filename)

        with open(datafile) as test_handle:
            seq_record = SeqIO.read(test_handle, "swiss")

        self.assertIsInstance(seq_record, SeqRecord)

        self.assertEqual(seq_record.id, "Q13639")
        self.assertEqual(seq_record.name, "5HT4R_HUMAN")
        self.assertEqual(
            seq_record.description,
            "RecName: Full=5-hydroxytryptamine receptor 4; Short=5-HT-4; Short=5-HT4; AltName: Full=Serotonin receptor 4;",
        )
        self.assertEqual(
            repr(seq_record.seq),
            "Seq('MDKLDANVSSEEGFGSVEKVVLLTFLSTVILMAILGNLLVMVAVCWDRQLRKIK...SDT')",
        )

        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(record.entry_name, "5HT4R_HUMAN")
        self.assertEqual(
            record.accessions,
            [
                "Q13639",
                "Q96KH9",
                "Q96KI0",
                "Q9H199",
                "Q9NY73",
                "Q9UBM6",
                "Q9UBT4",
                "Q9UE22",
                "Q9UE23",
                "Q9UQR6",
            ],
        )
        self.assertEqual(
            record.organism_classification,
            [
                "Eukaryota",
                "Metazoa",
                "Chordata",
                "Craniata",
                "Vertebrata",
                "Euteleostomi",
                "Mammalia",
                "Eutheria",
                "Euarchontoglires",
                "Primates",
                "Haplorrhini",
                "Catarrhini",
                "Hominidae",
                "Homo",
            ],
        )
        self.assertEqual(record.seqinfo, (388, 43761, "7FCFEC60E7BDF560"))

        self.assertEqual(len(record.features), 26)
        feature = record.features[0]
        self.assertEqual(feature.type, "CHAIN")
        self.assertEqual(feature.location.start, 0)
        self.assertEqual(feature.location.end, 388)
        self.assertEqual(
            feature.qualifiers["description"], "5-hydroxytryptamine receptor 4."
        )
        self.assertEqual(feature.id, "PRO_0000068965")
        feature = record.features[1]
        self.assertEqual(feature.type, "TOPO_DOM")
        self.assertEqual(feature.location.start, 0)
        self.assertEqual(feature.location.end, 19)
        self.assertEqual(
            feature.qualifiers["description"], "Extracellular (By similarity)."
        )
        self.assertIsNone(feature.id)
        feature = record.features[2]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 19)
        self.assertEqual(feature.location.end, 40)
        self.assertEqual(feature.qualifiers["description"], "1 (By similarity).")
        self.assertIsNone(feature.id)
        feature = record.features[3]
        self.assertEqual(feature.type, "TOPO_DOM")
        self.assertEqual(feature.location.start, 40)
        self.assertEqual(feature.location.end, 58)
        self.assertEqual(
            feature.qualifiers["description"], "Cytoplasmic (By similarity)."
        )
        self.assertIsNone(feature.id)
        feature = record.features[4]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 58)
        self.assertEqual(feature.location.end, 79)
        self.assertEqual(feature.qualifiers["description"], "2 (By similarity).")
        self.assertIsNone(feature.id)
        feature = record.features[5]
        self.assertEqual(feature.type, "TOPO_DOM")
        self.assertEqual(feature.location.start, 79)
        self.assertEqual(feature.location.end, 93)
        self.assertEqual(
            feature.qualifiers["description"], "Extracellular (By similarity)."
        )
        self.assertIsNone(feature.id)
        feature = record.features[6]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 93)
        self.assertEqual(feature.location.end, 116)
        self.assertEqual(feature.qualifiers["description"], "3 (By similarity).")
        self.assertIsNone(feature.id)
        feature = record.features[7]
        self.assertEqual(feature.type, "TOPO_DOM")
        self.assertEqual(feature.location.start, 116)
        self.assertEqual(feature.location.end, 137)
        self.assertEqual(
            feature.qualifiers["description"], "Cytoplasmic (By similarity)."
        )
        self.assertIsNone(feature.id)
        feature = record.features[8]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 137)
        self.assertEqual(feature.location.end, 158)
        self.assertEqual(feature.qualifiers["description"], "4 (By similarity).")
        self.assertIsNone(feature.id)
        feature = record.features[9]
        self.assertEqual(feature.type, "TOPO_DOM")
        self.assertEqual(feature.location.start, 158)
        self.assertEqual(feature.location.end, 192)
        self.assertEqual(
            feature.qualifiers["description"], "Extracellular (By similarity)."
        )
        self.assertIsNone(feature.id)
        feature = record.features[10]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 192)
        self.assertEqual(feature.location.end, 213)
        self.assertEqual(feature.qualifiers["description"], "5 (By similarity).")
        self.assertIsNone(feature.id)
        feature = record.features[11]
        self.assertEqual(feature.type, "TOPO_DOM")
        self.assertEqual(feature.location.start, 213)
        self.assertEqual(feature.location.end, 260)
        self.assertEqual(
            feature.qualifiers["description"], "Cytoplasmic (By similarity)."
        )
        self.assertIsNone(feature.id)
        feature = record.features[12]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 260)
        self.assertEqual(feature.location.end, 281)
        self.assertEqual(feature.qualifiers["description"], "6 (By similarity).")
        self.assertIsNone(feature.id)
        feature = record.features[13]
        self.assertEqual(feature.type, "TOPO_DOM")
        self.assertEqual(feature.location.start, 281)
        self.assertEqual(feature.location.end, 294)
        self.assertEqual(
            feature.qualifiers["description"], "Extracellular (By similarity)."
        )
        self.assertIsNone(feature.id)
        feature = record.features[14]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 294)
        self.assertEqual(feature.location.end, 315)
        self.assertEqual(feature.qualifiers["description"], "7 (By similarity).")
        self.assertIsNone(feature.id)
        feature = record.features[15]
        self.assertEqual(feature.type, "TOPO_DOM")
        self.assertEqual(feature.location.start, 315)
        self.assertEqual(feature.location.end, 388)
        self.assertEqual(
            feature.qualifiers["description"], "Cytoplasmic (By similarity)."
        )
        self.assertIsNone(feature.id)
        feature = record.features[16]
        self.assertEqual(feature.type, "LIPID")
        self.assertEqual(feature.location.start, 328)
        self.assertEqual(feature.location.end, 329)
        self.assertEqual(
            feature.qualifiers["description"], "S-palmitoyl cysteine (By similarity)."
        )
        self.assertIsNone(feature.id)
        feature = record.features[17]
        self.assertEqual(feature.type, "CARBOHYD")
        self.assertEqual(feature.location.start, 6)
        self.assertEqual(feature.location.end, 7)
        self.assertEqual(
            feature.qualifiers["description"], "N-linked (GlcNAc...) (Potential)."
        )
        self.assertIsNone(feature.id)
        feature = record.features[18]
        self.assertEqual(feature.type, "DISULFID")
        self.assertEqual(feature.location.start, 92)
        self.assertEqual(feature.location.end, 184)
        self.assertEqual(feature.qualifiers["description"], "By similarity.")
        self.assertIsNone(feature.id)
        feature = record.features[19]
        self.assertEqual(feature.type, "VAR_SEQ")
        self.assertEqual(feature.location.start, 168)
        self.assertEqual(feature.location.end, 169)
        self.assertEqual(
            feature.qualifiers["description"],
            "L -> LERSLNQGLGQDFHA (in isoform 5-HT4(F)).",
        )
        self.assertEqual(feature.id, "VSP_001845")
        feature = record.features[20]
        self.assertEqual(feature.type, "VAR_SEQ")
        self.assertEqual(feature.location.start, 358)
        self.assertEqual(feature.location.end, 388)
        self.assertEqual(
            feature.qualifiers["description"],
            "RDAVECGGQWESQCHPPATSPLVAAQPSDT -> SSGTETDRRNFGIRKRRLTKPS (in isoform 5-HT4(D)).",
        )
        self.assertEqual(feature.id, "VSP_001847")
        feature = record.features[21]
        self.assertEqual(feature.type, "VAR_SEQ")
        self.assertEqual(feature.location.start, 358)
        self.assertEqual(feature.location.end, 388)
        self.assertEqual(
            feature.qualifiers["description"],
            "RDAVECGGQWESQCHPPATSPLVAAQPSDT -> SGCSPVSSFLLLFCNRPVPV (in isoform 5-HT4(E)).",
        )
        self.assertEqual(feature.id, "VSP_001846")
        feature = record.features[22]
        self.assertEqual(feature.type, "VAR_SEQ")
        self.assertEqual(feature.location.start, 359)
        self.assertEqual(feature.location.end, 388)
        self.assertEqual(
            feature.qualifiers["description"],
            "DAVECGGQWESQCHPPATSPLVAAQPSDT -> YTVLHRGHHQELEKLPIHNDPESLESCF (in isoform 5-HT4(A)).",
        )
        self.assertEqual(feature.id, "VSP_001849")
        feature = record.features[23]
        self.assertEqual(feature.type, "VAR_SEQ")
        self.assertEqual(feature.location.start, 359)
        self.assertEqual(feature.location.end, 388)
        self.assertEqual(
            feature.qualifiers["description"],
            "DAVECGGQWESQCHPPATSPLVAAQPSDT -> F (in isoform 5-HT4(C)).",
        )
        self.assertEqual(feature.id, "VSP_001848")
        feature = record.features[24]
        self.assertEqual(feature.type, "VAR_SEQ")
        self.assertEqual(feature.location.start, 359)
        self.assertEqual(feature.location.end, 388)
        self.assertEqual(
            feature.qualifiers["description"], "Missing (in isoform 5-HT4(G))."
        )
        self.assertEqual(feature.id, "VSP_001850")
        feature = record.features[25]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 371)
        self.assertEqual(feature.location.end, 372)
        self.assertEqual(
            feature.qualifiers["description"], "C -> Y (in dbSNP:rs34826744)."
        )
        self.assertEqual(feature.id, "VAR_049364")
        self.assertEqual(len(record.references), 8)

        self.assertEqual(
            record.references[0].authors,
            "Blondel O., Gastineau M., Dahmoune Y., Langlois M., Fischmeister R.",
        )
        self.assertEqual(
            record.references[0].title,
            "Cloning, expression, and pharmacology of four human 5-hydroxytryptamine 4 receptor isoforms produced by alternative splicing in the carboxyl terminus.",
        )
        self.assertEqual(len(record.references[0].references), 2)
        self.assertEqual(record.references[0].references[0], ("MEDLINE", "98264328"))
        self.assertEqual(record.references[0].references[1], ("PubMed", "9603189"))
        self.assertEqual(
            record.references[1].authors,
            "Van den Wyngaert I., Gommeren W., Jurzak M., Verhasselt P., Gordon R., Leysen J., Luyten W., Bender E.",
        )
        self.assertEqual(
            record.references[1].title,
            "Cloning and expression of 5-HT4 receptor species and splice variants.",
        )
        self.assertEqual(len(record.references[1].references), 0)
        self.assertEqual(
            record.references[2].authors,
            "Claeysen S., Faye P., Sebben M., Lemaire S., Bockaert J., Dumuis A.",
        )
        self.assertEqual(
            record.references[2].title,
            "Cloning and expression of human 5-HT4S receptors. Effect of receptor density on their coupling to adenylyl cyclase.",
        )
        self.assertEqual(len(record.references[2].references), 2)
        self.assertEqual(record.references[2].references[0], ("MEDLINE", "98012006"))
        self.assertEqual(record.references[2].references[1], ("PubMed", "9351641"))
        self.assertEqual(
            record.references[3].authors,
            "Claeysen S., Sebben M., Becamel C., Bockaert J., Dumuis A.",
        )
        self.assertEqual(
            record.references[3].title,
            "Novel brain-specific 5-HT4 receptor splice variants show marked constitutive activity: role of the C-terminal intracellular domain.",
        )
        self.assertEqual(len(record.references[3].references), 2)
        self.assertEqual(record.references[3].references[0], ("MEDLINE", "99238795"))
        self.assertEqual(record.references[3].references[1], ("PubMed", "10220570"))
        self.assertEqual(
            record.references[4].authors,
            "Vilaro M.T., Domenech T., Palacios J.M., Mengod G.",
        )
        self.assertEqual(
            record.references[4].title,
            "Cloning and characterization of multiple human 5-HT4 receptor variants including a novel variant that lacks the alternatively spliced C-terminal exon.",
        )
        self.assertEqual(
            record.references[4].location,
            "Submitted (SEP-2000) to the EMBL/GenBank/DDBJ databases.",
        )
        self.assertEqual(len(record.references[4].comments), 1)
        self.assertEqual(record.references[4].comments[0], ("TISSUE", "Hippocampus"))
        self.assertEqual(len(record.references[4].positions), 1)
        self.assertEqual(
            record.references[4].positions[0],
            "NUCLEOTIDE SEQUENCE [MRNA] (ISOFORMS 5-HT4(A); 5-HT4(E) AND 5-HT4(G)).",
        )
        self.assertEqual(len(record.references[4].references), 0)
        self.assertEqual(len(record.references[5].positions), 1)
        self.assertEqual(
            record.references[5].positions[0],
            "NUCLEOTIDE SEQUENCE [LARGE SCALE MRNA] (ISOFORM 5-HT4(B)).",
        )
        self.assertEqual(len(record.references[5].references), 2)
        self.assertEqual(record.references[5].references[0], ("PubMed", "15489334"))
        self.assertEqual(
            record.references[5].references[1], ("DOI", "10.1101/gr.2596504")
        )
        self.assertEqual(record.references[5].authors, "The MGC Project Team")
        self.assertEqual(
            record.references[5].title,
            "The status, quality, and expansion of the NIH full-length cDNA project: the Mammalian Gene Collection (MGC).",
        )
        self.assertEqual(
            record.references[5].location, "Genome Res. 14:2121-2127(2004)."
        )
        self.assertEqual(
            record.references[6].authors,
            "Bender E., Pindon A., van Oers I., Zhang Y.B., Gommeren W., Verhasselt P., Jurzak M., Leysen J., Luyten W.",
        )
        self.assertEqual(
            record.references[6].title,
            "Structure of the human serotonin 5-HT4 receptor gene and cloning of a novel 5-HT4 splice variant.",
        )
        self.assertEqual(len(record.references[6].references), 3)
        self.assertEqual(record.references[6].references[0], ("MEDLINE", "20110418"))
        self.assertEqual(record.references[6].references[1], ("PubMed", "10646498"))
        self.assertEqual(
            record.references[6].references[2],
            ("DOI", "10.1046/j.1471-4159.2000.740478.x"),
        )
        self.assertEqual(
            record.references[7].authors,
            "Ullmer C., Schmuck K., Kalkman H.O., Luebbert H.",
        )
        self.assertEqual(
            record.references[7].title,
            "Expression of serotonin receptor mRNAs in blood vessels.",
        )
        self.assertEqual(len(record.references[7].references), 3)
        self.assertEqual(record.references[7].references[0], ("MEDLINE", "95385798"))
        self.assertEqual(record.references[7].references[1], ("PubMed", "7656980"))
        self.assertEqual(
            record.references[7].references[2], ("DOI", "10.1016/0014-5793(95)00828-W")
        )

        # Check the two parsers agree on the essentials
        self.assertEqual(str(seq_record.seq), record.sequence)
        self.assertEqual(seq_record.description, record.description)
        self.assertEqual(seq_record.name, record.entry_name)
        self.assertIn(seq_record.id, record.accessions)

        # Now try using the iterator - note that all these
        # test cases have only one record.

        # With the SequenceParser
        records = list(SeqIO.parse(datafile, "swiss"))

        self.assertEqual(len(records), 1)
        self.assertIsInstance(records[0], SeqRecord)

        # Check matches what we got earlier without the iterator:
        self.assertEqual(str(records[0].seq), str(seq_record.seq))
        self.assertEqual(records[0].description, seq_record.description)
        self.assertEqual(records[0].name, seq_record.name)
        self.assertEqual(records[0].id, seq_record.id)

        # With the RecordParser
        with open(datafile) as test_handle:
            records = list(SwissProt.parse(test_handle))

        self.assertEqual(len(records), 1)
        self.assertIsInstance(records[0], SwissProt.Record)

        # Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].sequence, record.sequence)
        self.assertEqual(records[0].description, record.description)
        self.assertEqual(records[0].entry_name, record.entry_name)
        self.assertEqual(records[0].accessions, record.accessions)

    def test_sp011(self):
        """Parsing SwissProt file sp011."""
        filename = "sp011"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)
        seq_record = SeqIO.read(datafile, "swiss")

        self.assertIsInstance(seq_record, SeqRecord)

        self.assertEqual(seq_record.id, "P16235")
        self.assertEqual(seq_record.name, "LSHR_RAT")
        self.assertEqual(
            seq_record.description,
            "LUTROPIN-CHORIOGONADOTROPIC HORMONE RECEPTOR PRECURSOR (LH/CG-R) (LSH-R) (LUTEINIZING HORMONE RECEPTOR).",
        )
        self.assertEqual(
            repr(seq_record.seq),
            "Seq('MGRRVPALRQLLVLAVLLLKPSQLQSRELSGSRCPEPCDCAPDGALRCPGPRAG...LTH')",
        )

        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(record.entry_name, "LSHR_RAT")
        self.assertEqual(
            record.accessions, ["P16235", "P70646", "Q63807", "Q63808", "Q63809"]
        )
        self.assertEqual(
            record.organism_classification,
            [
                "Eukaryota",
                "Metazoa",
                "Chordata",
                "Craniata",
                "Vertebrata",
                "Euteleostomi",
                "Mammalia",
                "Eutheria",
                "Rodentia",
                "Sciurognathi",
                "Muridae",
                "Murinae",
                "Rattus",
            ],
        )
        self.assertEqual(record.seqinfo, (700, 78035, "31807E73BAC94F1F"))

        self.assertEqual(len(record.features), 52)
        feature = record.features[0]
        self.assertEqual(feature.type, "SIGNAL")
        self.assertEqual(feature.location.start, 0)
        self.assertEqual(feature.location.end, 26)
        self.assertEqual(feature.qualifiers["description"], "")
        self.assertIsNone(feature.id)
        feature = record.features[1]
        self.assertEqual(feature.type, "CHAIN")
        self.assertEqual(feature.location.start, 26)
        self.assertEqual(feature.location.end, 700)
        self.assertEqual(
            feature.qualifiers["description"],
            "LUTROPIN-CHORIOGONADOTROPIC HORMONE RECEPTOR.",
        )
        self.assertIsNone(feature.id)
        feature = record.features[2]
        self.assertEqual(feature.type, "DOMAIN")
        self.assertEqual(feature.location.start, 26)
        self.assertEqual(feature.location.end, 362)
        self.assertEqual(
            feature.qualifiers["description"], "EXTRACELLULAR (POTENTIAL)."
        )
        self.assertIsNone(feature.id)
        feature = record.features[3]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 362)
        self.assertEqual(feature.location.end, 390)
        self.assertEqual(feature.qualifiers["description"], "1 (POTENTIAL).")
        self.assertIsNone(feature.id)
        feature = record.features[4]
        self.assertEqual(feature.type, "DOMAIN")
        self.assertEqual(feature.location.start, 390)
        self.assertEqual(feature.location.end, 399)
        self.assertEqual(feature.qualifiers["description"], "CYTOPLASMIC (POTENTIAL).")
        self.assertIsNone(feature.id)
        feature = record.features[5]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 399)
        self.assertEqual(feature.location.end, 422)
        self.assertEqual(feature.qualifiers["description"], "2 (POTENTIAL).")
        self.assertIsNone(feature.id)
        feature = record.features[6]
        self.assertEqual(feature.type, "DOMAIN")
        self.assertEqual(feature.location.start, 422)
        self.assertEqual(feature.location.end, 443)
        self.assertEqual(
            feature.qualifiers["description"], "EXTRACELLULAR (POTENTIAL)."
        )
        self.assertIsNone(feature.id)
        feature = record.features[7]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 443)
        self.assertEqual(feature.location.end, 466)
        self.assertEqual(feature.qualifiers["description"], "3 (POTENTIAL).")
        self.assertIsNone(feature.id)
        feature = record.features[8]
        self.assertEqual(feature.type, "DOMAIN")
        self.assertEqual(feature.location.start, 466)
        self.assertEqual(feature.location.end, 486)
        self.assertEqual(feature.qualifiers["description"], "CYTOPLASMIC (POTENTIAL).")
        self.assertIsNone(feature.id)
        feature = record.features[9]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 486)
        self.assertEqual(feature.location.end, 509)
        self.assertEqual(feature.qualifiers["description"], "4 (POTENTIAL).")
        self.assertIsNone(feature.id)
        feature = record.features[10]
        self.assertEqual(feature.type, "DOMAIN")
        self.assertEqual(feature.location.start, 509)
        self.assertEqual(feature.location.end, 529)
        self.assertEqual(
            feature.qualifiers["description"], "EXTRACELLULAR (POTENTIAL)."
        )
        self.assertIsNone(feature.id)
        feature = record.features[11]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 529)
        self.assertEqual(feature.location.end, 551)
        self.assertEqual(feature.qualifiers["description"], "5 (POTENTIAL).")
        self.assertIsNone(feature.id)
        feature = record.features[12]
        self.assertEqual(feature.type, "DOMAIN")
        self.assertEqual(feature.location.start, 551)
        self.assertEqual(feature.location.end, 574)
        self.assertEqual(feature.qualifiers["description"], "CYTOPLASMIC (POTENTIAL).")
        self.assertIsNone(feature.id)
        feature = record.features[13]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 574)
        self.assertEqual(feature.location.end, 598)
        self.assertEqual(feature.qualifiers["description"], "6 (POTENTIAL).")
        self.assertIsNone(feature.id)
        feature = record.features[14]
        self.assertEqual(feature.type, "DOMAIN")
        self.assertEqual(feature.location.start, 598)
        self.assertEqual(feature.location.end, 609)
        self.assertEqual(
            feature.qualifiers["description"], "EXTRACELLULAR (POTENTIAL)."
        )
        self.assertIsNone(feature.id)
        feature = record.features[15]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 609)
        self.assertEqual(feature.location.end, 631)
        self.assertEqual(feature.qualifiers["description"], "7 (POTENTIAL).")
        self.assertIsNone(feature.id)
        feature = record.features[16]
        self.assertEqual(feature.type, "DOMAIN")
        self.assertEqual(feature.location.start, 631)
        self.assertEqual(feature.location.end, 700)
        self.assertEqual(feature.qualifiers["description"], "CYTOPLASMIC (POTENTIAL).")
        self.assertIsNone(feature.id)
        feature = record.features[17]
        self.assertEqual(feature.type, "REPEAT")
        self.assertEqual(feature.location.start, 51)
        self.assertEqual(feature.location.end, 75)
        self.assertEqual(feature.qualifiers["description"], "LRR 1.")
        self.assertIsNone(feature.id)
        feature = record.features[18]
        self.assertEqual(feature.type, "REPEAT")
        self.assertEqual(feature.location.start, 125)
        self.assertEqual(feature.location.end, 150)
        self.assertEqual(feature.qualifiers["description"], "LRR 2.")
        self.assertIsNone(feature.id)
        feature = record.features[19]
        self.assertEqual(feature.type, "REPEAT")
        self.assertEqual(feature.location.start, 151)
        self.assertEqual(feature.location.end, 175)
        self.assertEqual(feature.qualifiers["description"], "LRR 3.")
        self.assertIsNone(feature.id)
        feature = record.features[20]
        self.assertEqual(feature.type, "REPEAT")
        self.assertEqual(feature.location.start, 175)
        self.assertEqual(feature.location.end, 200)
        self.assertEqual(feature.qualifiers["description"], "LRR 4.")
        self.assertIsNone(feature.id)
        feature = record.features[21]
        self.assertEqual(feature.type, "REPEAT")
        self.assertEqual(feature.location.start, 201)
        self.assertEqual(feature.location.end, 224)
        self.assertEqual(feature.qualifiers["description"], "LRR 5.")
        self.assertIsNone(feature.id)
        feature = record.features[22]
        self.assertEqual(feature.type, "REPEAT")
        self.assertEqual(feature.location.start, 224)
        self.assertEqual(feature.location.end, 248)
        self.assertEqual(feature.qualifiers["description"], "LRR 6.")
        self.assertIsNone(feature.id)
        feature = record.features[23]
        self.assertEqual(feature.type, "REPEAT")
        self.assertEqual(feature.location.start, 249)
        self.assertEqual(feature.location.end, 271)
        self.assertEqual(feature.qualifiers["description"], "LRR 7.")
        self.assertIsNone(feature.id)
        feature = record.features[24]
        self.assertEqual(feature.type, "DISULFID")
        self.assertEqual(feature.location.start, 442)
        self.assertEqual(feature.location.end, 518)
        self.assertEqual(feature.qualifiers["description"], "BY SIMILARITY.")
        self.assertIsNone(feature.id)
        feature = record.features[25]
        self.assertEqual(feature.type, "CARBOHYD")
        self.assertEqual(feature.location.start, 102)
        self.assertEqual(feature.location.end, 103)
        self.assertEqual(
            feature.qualifiers["description"], "N-LINKED (GLCNAC...) (POTENTIAL)."
        )
        self.assertIsNone(feature.id)
        feature = record.features[26]
        self.assertEqual(feature.type, "CARBOHYD")
        self.assertEqual(feature.location.start, 177)
        self.assertEqual(feature.location.end, 178)
        self.assertEqual(
            feature.qualifiers["description"], "N-LINKED (GLCNAC...) (POTENTIAL)."
        )
        self.assertIsNone(feature.id)
        feature = record.features[27]
        self.assertEqual(feature.type, "CARBOHYD")
        self.assertEqual(feature.location.start, 198)
        self.assertEqual(feature.location.end, 199)
        self.assertEqual(
            feature.qualifiers["description"], "N-LINKED (GLCNAC...) (POTENTIAL)."
        )
        self.assertIsNone(feature.id)
        feature = record.features[28]
        self.assertEqual(feature.type, "CARBOHYD")
        self.assertEqual(feature.location.start, 294)
        self.assertEqual(feature.location.end, 295)
        self.assertEqual(
            feature.qualifiers["description"], "N-LINKED (GLCNAC...) (POTENTIAL)."
        )
        self.assertIsNone(feature.id)
        feature = record.features[29]
        self.assertEqual(feature.type, "CARBOHYD")
        self.assertEqual(feature.location.start, 302)
        self.assertEqual(feature.location.end, 303)
        self.assertEqual(
            feature.qualifiers["description"], "N-LINKED (GLCNAC...) (POTENTIAL)."
        )
        self.assertIsNone(feature.id)
        feature = record.features[30]
        self.assertEqual(feature.type, "CARBOHYD")
        self.assertEqual(feature.location.start, 316)
        self.assertEqual(feature.location.end, 317)
        self.assertEqual(
            feature.qualifiers["description"], "N-LINKED (GLCNAC...) (POTENTIAL)."
        )
        self.assertIsNone(feature.id)
        feature = record.features[31]
        self.assertEqual(feature.type, "VARSPLIC")
        self.assertEqual(feature.location.start, 82)
        self.assertEqual(feature.location.end, 132)
        self.assertEqual(
            feature.qualifiers["description"], "MISSING (IN ISOFORM 1950)."
        )
        self.assertIsNone(feature.id)
        feature = record.features[32]
        self.assertEqual(feature.type, "VARSPLIC")
        self.assertEqual(feature.location.start, 132)
        self.assertEqual(feature.location.end, 157)
        self.assertEqual(
            feature.qualifiers["description"], "MISSING (IN ISOFORM 1759)."
        )
        self.assertIsNone(feature.id)
        feature = record.features[33]
        self.assertEqual(feature.type, "VARSPLIC")
        self.assertEqual(feature.location.start, 183)
        self.assertEqual(feature.location.end, 700)
        self.assertEqual(feature.qualifiers["description"], "MISSING (IN ISOFORM C2).")
        self.assertIsNone(feature.id)
        feature = record.features[34]
        self.assertEqual(feature.type, "VARSPLIC")
        self.assertEqual(feature.location.start, 231)
        self.assertEqual(feature.location.end, 251)
        self.assertEqual(
            feature.qualifiers["description"],
            "DISSTKLQALPSHGLESIQT -> PCRATGWSPFRRSSPCLPTH (IN ISOFORM 2075).",
        )
        self.assertIsNone(feature.id)
        feature = record.features[35]
        self.assertEqual(feature.type, "VARSPLIC")
        self.assertEqual(feature.location.start, 231)
        self.assertEqual(feature.location.end, 293)
        self.assertEqual(
            feature.qualifiers["description"],
            "MISSING (IN ISOFORM E/A2, ISOFORM EB AND ISOFORM B1).",
        )
        self.assertIsNone(feature.id)
        feature = record.features[36]
        self.assertEqual(feature.type, "VARSPLIC")
        self.assertEqual(feature.location.start, 251)
        self.assertEqual(feature.location.end, 700)
        self.assertEqual(
            feature.qualifiers["description"], "MISSING (IN ISOFORM 2075)."
        )
        self.assertIsNone(feature.id)
        feature = record.features[37]
        self.assertEqual(feature.type, "VARSPLIC")
        self.assertEqual(feature.location.start, 293)
        self.assertEqual(feature.location.end, 367)
        self.assertEqual(
            feature.qualifiers["description"],
            "QNFSFSIFENFSKQCESTVRKADNETLYSAIFEENELSGWDYDYGFCSPKTLQCAPEPDAFNPCEDIMGYAFLR -> IFHFPFLKTSPNNAKAQLEKQITRRFIPPSLRRMNSVAGIMIMASVHPRHSNVLQNQMLSTPVKILWAMPSLGS (IN ISOFORM B1 AND ISOFORM B3).",
        )
        self.assertIsNone(feature.id)
        feature = record.features[38]
        self.assertEqual(feature.type, "VARSPLIC")
        self.assertEqual(feature.location.start, 293)
        self.assertEqual(feature.location.end, 294)
        self.assertEqual(feature.qualifiers["description"], "Q -> P (IN ISOFORM C1).")
        self.assertIsNone(feature.id)
        feature = record.features[39]
        self.assertEqual(feature.type, "VARSPLIC")
        self.assertEqual(feature.location.start, 294)
        self.assertEqual(feature.location.end, 700)
        self.assertEqual(feature.qualifiers["description"], "MISSING (IN ISOFORM C1).")
        self.assertIsNone(feature.id)
        feature = record.features[40]
        self.assertEqual(feature.type, "VARSPLIC")
        self.assertEqual(feature.location.start, 320)
        self.assertEqual(feature.location.end, 342)
        self.assertEqual(
            feature.qualifiers["description"],
            "YSAIFEENELSGWDYDYGFCSP -> LHGALPAAHCLRGLPNKRPVL (IN ISOFORM 1834, ISOFORM 1759 AND ISOFORM EB).",
        )
        self.assertIsNone(feature.id)
        feature = record.features[41]
        self.assertEqual(feature.type, "VARSPLIC")
        self.assertEqual(feature.location.start, 342)
        self.assertEqual(feature.location.end, 700)
        self.assertEqual(
            feature.qualifiers["description"],
            "MISSING (IN ISOFORMS 1834, ISOFORM 1759 AND ISOFORM EB).",
        )
        self.assertIsNone(feature.id)
        feature = record.features[42]
        self.assertEqual(feature.type, "VARSPLIC")
        self.assertEqual(feature.location.start, 367)
        self.assertEqual(feature.location.end, 700)
        self.assertEqual(
            feature.qualifiers["description"], "MISSING (IN ISOFORM B1 AND ISOFORM B3)."
        )
        self.assertIsNone(feature.id)
        feature = record.features[43]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 81)
        self.assertEqual(feature.location.end, 82)
        self.assertEqual(feature.qualifiers["description"], "I -> M (IN ISOFORM 1950).")
        self.assertIsNone(feature.id)
        feature = record.features[44]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 178)
        self.assertEqual(feature.location.end, 179)
        self.assertEqual(feature.qualifiers["description"], "E -> G (IN ISOFORM 1759).")
        self.assertIsNone(feature.id)
        feature = record.features[45]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 232)
        self.assertEqual(feature.location.end, 233)
        self.assertEqual(feature.qualifiers["description"], "I -> T (IN ISOFORM 1950).")
        self.assertIsNone(feature.id)
        feature = record.features[46]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 645)
        self.assertEqual(feature.location.end, 646)
        self.assertEqual(feature.qualifiers["description"], "G -> S (IN ISOFORM 1950).")
        self.assertIsNone(feature.id)
        feature = record.features[47]
        self.assertEqual(feature.type, "MUTAGEN")
        self.assertEqual(feature.location.start, 408)
        self.assertEqual(feature.location.end, 409)
        self.assertEqual(
            feature.qualifiers["description"], "D->N: SIGNIFICANT REDUCTION OF BINDING."
        )
        self.assertIsNone(feature.id)
        feature = record.features[48]
        self.assertEqual(feature.type, "MUTAGEN")
        self.assertEqual(feature.location.start, 435)
        self.assertEqual(feature.location.end, 436)
        self.assertEqual(
            feature.qualifiers["description"],
            "D->N: NO CHANGE IN BINDING OR CAMP PROD.",
        )
        self.assertIsNone(feature.id)
        feature = record.features[49]
        self.assertEqual(feature.type, "MUTAGEN")
        self.assertEqual(feature.location.start, 454)
        self.assertEqual(feature.location.end, 455)
        self.assertEqual(
            feature.qualifiers["description"],
            "E->Q: NO CHANGE IN BINDING OR CAMP PROD.",
        )
        self.assertIsNone(feature.id)
        feature = record.features[50]
        self.assertEqual(feature.type, "MUTAGEN")
        self.assertEqual(feature.location.start, 581)
        self.assertEqual(feature.location.end, 582)
        self.assertEqual(
            feature.qualifiers["description"],
            "D->N: NO CHANGE IN BINDING OR CAMP PROD.",
        )
        self.assertIsNone(feature.id)
        feature = record.features[51]
        self.assertEqual(feature.type, "CONFLICT")
        self.assertEqual(feature.location.start, 32)
        self.assertEqual(feature.location.end, 33)
        self.assertEqual(feature.qualifiers["description"], "R -> L (IN REF. 7).")
        self.assertIsNone(feature.id)

        self.assertEqual(len(record.references), 8)
        self.assertEqual(
            record.references[0].authors,
            "McFarland K.C., Sprengel R., Phillips H.S., Koehler M., Rosemblit N., Nikolics K., Segaloff D.L., Seeburg P.H.",
        )
        self.assertEqual(
            record.references[0].title,
            "Lutropin-choriogonadotropin receptor: an unusual member of the G protein-coupled receptor family.",
        )
        self.assertEqual(len(record.references[0].references), 2)
        self.assertEqual(record.references[0].references[0], ("MEDLINE", "89332512"))
        self.assertEqual(record.references[0].references[1], ("PubMed", "2502842"))
        self.assertEqual(
            record.references[1].authors,
            "Aatsinki J.T., Pietila E.M., Lakkakorpi J.T., Rajaniemi H.J.",
        )
        self.assertEqual(
            record.references[1].title,
            "Expression of the LH/CG receptor gene in rat ovarian tissue is regulated by an extensive alternative splicing of the primary transcript.",
        )
        self.assertEqual(len(record.references[1].references), 2)
        self.assertEqual(record.references[1].references[0], ("MEDLINE", "92347604"))
        self.assertEqual(record.references[1].references[1], ("PubMed", "1353463"))
        self.assertEqual(
            record.references[2].authors, "Koo Y.B., Slaughter R.G., Ji T.H."
        )
        self.assertEqual(
            record.references[2].title,
            "Structure of the luteinizing hormone receptor gene and multiple exons of the coding sequence.",
        )
        self.assertEqual(len(record.references[2].references), 2)
        self.assertEqual(record.references[2].references[0], ("MEDLINE", "91209270"))
        self.assertEqual(record.references[2].references[1], ("PubMed", "2019252"))
        self.assertEqual(
            record.references[3].authors, "Bernard M.P., Myers R.V., Moyle W.R."
        )
        self.assertEqual(
            record.references[3].title,
            "Cloning of rat lutropin (LH) receptor analogs lacking the soybean lectin domain.",
        )
        self.assertEqual(len(record.references[3].references), 2)
        self.assertEqual(record.references[3].references[0], ("MEDLINE", "91006819"))
        self.assertEqual(record.references[3].references[1], ("PubMed", "1976554"))
        self.assertEqual(
            record.references[4].authors,
            "Segaloff D.L., Sprengel R., Nikolics K., Ascoli M.",
        )
        self.assertEqual(
            record.references[4].title,
            "Structure of the lutropin/choriogonadotropin receptor.",
        )
        self.assertEqual(len(record.references[4].references), 2)
        self.assertEqual(record.references[4].references[0], ("MEDLINE", "91126285"))
        self.assertEqual(record.references[4].references[1], ("PubMed", "2281186"))
        self.assertEqual(
            record.references[5].authors,
            "Tsai-Morris C.H., Buczko E., Wang W., Dufau M.L.",
        )
        self.assertEqual(
            record.references[5].title,
            "Intronic nature of the rat luteinizing hormone receptor gene defines a soluble receptor subspecies with hormone binding activity.",
        )
        self.assertEqual(len(record.references[5].references), 2)
        self.assertEqual(record.references[5].references[0], ("MEDLINE", "91060531"))
        self.assertEqual(record.references[5].references[1], ("PubMed", "2174034"))
        self.assertEqual(record.references[6].authors, "Roche P.C., Ryan R.J.")
        self.assertEqual(
            record.references[6].title,
            "Purification, characterization, and amino-terminal sequence of rat ovarian receptor for luteinizing hormone/human choriogonadotropin.",
        )
        self.assertEqual(len(record.references[6].references), 2)
        self.assertEqual(record.references[6].references[0], ("MEDLINE", "89174723"))
        self.assertEqual(record.references[6].references[1], ("PubMed", "2925659"))
        self.assertEqual(record.references[7].authors, "Ji I., Ji T.H.")
        self.assertEqual(
            record.references[7].title,
            "Asp383 in the second transmembrane domain of the lutropin receptor is important for high affinity hormone binding and cAMP production.",
        )
        self.assertEqual(len(record.references[7].references), 2)
        self.assertEqual(record.references[7].references[0], ("MEDLINE", "91332007"))
        self.assertEqual(record.references[7].references[1], ("PubMed", "1714448"))

        # Check the two parsers agree on the essentials
        self.assertEqual(str(seq_record.seq), record.sequence)
        self.assertEqual(seq_record.description, record.description)
        self.assertEqual(seq_record.name, record.entry_name)
        self.assertIn(seq_record.id, record.accessions)

        # Now try using the iterator - note that all these
        # test cases have only one record.

        # With the SequenceParser
        records = list(SeqIO.parse(datafile, "swiss"))

        self.assertEqual(len(records), 1)
        self.assertIsInstance(records[0], SeqRecord)

        # Check matches what we got earlier without the iterator:
        self.assertEqual(str(records[0].seq), str(seq_record.seq))
        self.assertEqual(records[0].description, seq_record.description)
        self.assertEqual(records[0].name, seq_record.name)
        self.assertEqual(records[0].id, seq_record.id)

        # With the RecordParser
        with open(datafile) as test_handle:
            records = list(SwissProt.parse(test_handle))

        self.assertEqual(len(records), 1)
        self.assertIsInstance(records[0], SwissProt.Record)

        # Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].sequence, record.sequence)
        self.assertEqual(records[0].description, record.description)
        self.assertEqual(records[0].entry_name, record.entry_name)
        self.assertEqual(records[0].accessions, record.accessions)

    def test_sp012(self):
        """Parsing SwissProt file sp012."""
        filename = "sp012"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)
        seq_record = SeqIO.read(datafile, "swiss")

        self.assertIsInstance(seq_record, SeqRecord)

        self.assertEqual(seq_record.id, "Q9Y736")
        self.assertEqual(seq_record.name, "Q9Y736")
        self.assertEqual(seq_record.description, "UBIQUITIN.")
        self.assertEqual(
            repr(seq_record.seq),
            "Seq('MQIFVKTLTGKTITLEVESSDTIDNVKTKIQDKEGIPPDQQRLIFAGKQLEDGR...GGN')",
        )

        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(record.entry_name, "Q9Y736")
        self.assertEqual(record.accessions, ["Q9Y736"])
        self.assertEqual(
            record.organism_classification,
            [
                "Eukaryota",
                "Fungi",
                "Ascomycota",
                "Pezizomycotina",
                "Eurotiomycetes",
                "Onygenales",
                "Arthrodermataceae",
                "mitosporic Arthrodermataceae",
                "Trichophyton",
            ],
        )
        self.assertEqual(record.seqinfo, (153, 17238, "01153CF30C2DEDFF"))

        self.assertEqual(len(record.features), 0)

        self.assertEqual(len(record.references), 2)
        self.assertEqual(
            record.references[0].authors,
            "Kano R., Nakamura Y., Watanabe S., Hasegawa A.",
        )
        self.assertEqual(
            record.references[0].title,
            "Trichophyton mentagrophytes mRNA for ubiquitin.",
        )
        self.assertEqual(len(record.references[0].references), 0)
        self.assertEqual(record.references[1].authors, "Kano R.")
        self.assertEqual(
            record.references[1].title,
            "Microsporum canis mRNA for ubiquitin, complete cds.",
        )
        self.assertEqual(len(record.references[1].references), 0)

        # Check the two parsers agree on the essentials
        self.assertEqual(str(seq_record.seq), record.sequence)
        self.assertEqual(seq_record.description, record.description)
        self.assertEqual(seq_record.name, record.entry_name)
        self.assertIn(seq_record.id, record.accessions)

        # Now try using the iterator - note that all these
        # test cases have only one record.

        # With the SequenceParser
        records = list(SeqIO.parse(datafile, "swiss"))

        self.assertEqual(len(records), 1)
        self.assertIsInstance(records[0], SeqRecord)

        # Check matches what we got earlier without the iterator:
        self.assertEqual(str(records[0].seq), str(seq_record.seq))
        self.assertEqual(records[0].description, seq_record.description)
        self.assertEqual(records[0].name, seq_record.name)
        self.assertEqual(records[0].id, seq_record.id)

        # With the RecordParser
        with open(datafile) as test_handle:
            records = list(SwissProt.parse(test_handle))

        self.assertEqual(len(records), 1)
        self.assertIsInstance(records[0], SwissProt.Record)

        # Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].sequence, record.sequence)
        self.assertEqual(records[0].description, record.description)
        self.assertEqual(records[0].entry_name, record.entry_name)
        self.assertEqual(records[0].accessions, record.accessions)

    def test_sp013(self):
        """Parsing SwissProt file sp013."""
        filename = "sp013"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)
        seq_record = SeqIO.read(datafile, "swiss")

        self.assertIsInstance(seq_record, SeqRecord)

        self.assertEqual(seq_record.id, "P82909")
        self.assertEqual(seq_record.name, "P82909")
        self.assertEqual(
            seq_record.description, "MITOCHONDRIAL 28S RIBOSOMAL PROTEIN S36 (MRP-S36)."
        )
        self.assertEqual(
            repr(seq_record.seq),
            "Seq('MGSKMASASRVVQVVKPHTPLIRFPDRRDNPKPNVSEALRSAGLPSHSSVISQH...GPE')",
        )

        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(record.entry_name, "P82909")
        self.assertEqual(record.accessions, ["P82909"])
        self.assertEqual(
            record.organism_classification,
            [
                "Eukaryota",
                "Metazoa",
                "Chordata",
                "Craniata",
                "Vertebrata",
                "Euteleostomi",
                "Mammalia",
                "Eutheria",
                "Primates",
                "Catarrhini",
                "Hominidae",
                "Homo",
            ],
        )
        self.assertEqual(record.seqinfo, (102, 11335, "83EF107B42E2FCFD"))

        self.assertEqual(len(record.features), 0)

        self.assertEqual(len(record.references), 2)
        self.assertEqual(record.references[0].authors, "Strausberg R.")
        self.assertEqual(record.references[0].title, "")
        self.assertEqual(len(record.references[0].references), 0)
        self.assertEqual(
            record.references[1].authors,
            "Koc E.C., Burkhart W., Blackburn K., Moseley A., Spremulli L.L.",
        )
        self.assertEqual(
            record.references[1].title,
            "The small subunit of the mammalian mitochondrial ribosome. Identification of the full complement ribosomal proteins present.",
        )
        self.assertEqual(len(record.references[1].references), 0)

        # Check the two parsers agree on the essentials
        self.assertEqual(str(seq_record.seq), record.sequence)
        self.assertEqual(seq_record.description, record.description)
        self.assertEqual(seq_record.name, record.entry_name)
        self.assertIn(seq_record.id, record.accessions)

        # Now try using the iterator - note that all these
        # test cases have only one record.

        # With the SequenceParser
        records = list(SeqIO.parse(datafile, "swiss"))

        self.assertEqual(len(records), 1)
        self.assertIsInstance(records[0], SeqRecord)

        # Check matches what we got earlier without the iterator:
        self.assertEqual(str(records[0].seq), str(seq_record.seq))
        self.assertEqual(records[0].description, seq_record.description)
        self.assertEqual(records[0].name, seq_record.name)
        self.assertEqual(records[0].id, seq_record.id)

        # With the RecordParser
        with open(datafile) as test_handle:
            records = list(SwissProt.parse(test_handle))

        self.assertEqual(len(records), 1)
        self.assertIsInstance(records[0], SwissProt.Record)

        # Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].sequence, record.sequence)
        self.assertEqual(records[0].description, record.description)
        self.assertEqual(records[0].entry_name, record.entry_name)
        self.assertEqual(records[0].accessions, record.accessions)

    def test_P60137(self):
        """Parsing SwissProt file P60137.txt."""
        filename = "P60137.txt"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)
        seq_record = SeqIO.read(datafile, "swiss")

        self.assertIsInstance(seq_record, SeqRecord)

        self.assertEqual(seq_record.id, "P60137")
        self.assertEqual(seq_record.name, "PSBL_ORYSJ")
        self.assertEqual(
            seq_record.description,
            "RecName: Full=Photosystem II reaction center protein L {ECO:0000255|HAMAP-Rule:MF_01317}; Short=PSII-L {ECO:0000255|HAMAP-Rule:MF_01317};",
        )
        self.assertEqual(
            repr(seq_record.seq), "Seq('MTQSNPNEQNVELNRTSLYWGLLLIFVLAVLFSNYFFN')",
        )

        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(record.entry_name, "PSBL_ORYSJ")
        self.assertEqual(
            record.accessions, ["P60137", "O47030", "P12166", "P12167", "Q34007"]
        )
        self.assertEqual(
            record.organism_classification,
            [
                "Eukaryota",
                "Viridiplantae",
                "Streptophyta",
                "Embryophyta",
                "Tracheophyta",
                "Spermatophyta",
                "Magnoliopsida",
                "Liliopsida",
                "Poales",
                "Poaceae",
                "BOP clade",
                "Oryzoideae",
                "Oryzeae",
                "Oryzinae",
                "Oryza",
                "Oryza sativa",
            ],
        )
        self.assertEqual(record.seqinfo, (38, 4497, "55537AEC50D25E8D"))

        self.assertEqual(len(record.features), 2)
        feature = record.features[0]
        self.assertEqual(feature.type, "CHAIN")
        self.assertEqual(feature.location.start, 0)
        self.assertEqual(feature.location.end, 38)
        self.assertEqual(
            feature.qualifiers["note"], "Photosystem II reaction center protein L"
        )
        self.assertEqual(feature.id, "PRO_0000219754")
        feature = record.features[1]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 16)
        self.assertEqual(feature.location.end, 37)
        self.assertEqual(feature.qualifiers["note"], "Helical")
        self.assertEqual(
            feature.qualifiers["evidence"], "ECO:0000255|HAMAP-Rule:MF_01317"
        )

        self.assertEqual(len(record.references), 2)
        self.assertEqual(
            record.references[0].authors,
            "Hiratsuka J., Shimada H., Whittier R., Ishibashi T., Sakamoto M., Mori M., Kondo C., Honji Y., Sun C.-R., Meng B.-Y., Li Y.-Q., Kanno A., Nishizawa Y., Hirai A., Shinozaki K., Sugiura M.",
        )
        self.assertEqual(
            record.references[0].title,
            "The complete sequence of the rice (Oryza sativa) chloroplast genome: intermolecular recombination between distinct tRNA genes accounts for a major plastid DNA inversion during the evolution of the cereals.",
        )
        self.assertEqual(len(record.references[0].references), 2)
        self.assertEqual(record.references[0].references[0], ("PubMed", "2770692"))
        self.assertEqual(
            record.references[0].references[1], ("DOI", "10.1007/bf02464880")
        )
        self.assertEqual(
            record.references[1].authors,
            "Tang J., Xia H., Cao M., Zhang X., Zeng W., Hu S., Tong W., Wang J., Wang J., Yu J., Yang H., Zhu L.",
        )
        self.assertEqual(
            record.references[1].title, "A comparison of rice chloroplast genomes."
        )
        self.assertEqual(len(record.references[1].references), 2)
        self.assertEqual(record.references[1].references[0], ("PubMed", "15122023"))
        self.assertEqual(
            record.references[1].references[1], ("DOI", "10.1104/pp.103.031245")
        )

        # Check the two parsers agree on the essentials
        self.assertEqual(str(seq_record.seq), record.sequence)
        self.assertEqual(seq_record.description, record.description)
        self.assertEqual(seq_record.name, record.entry_name)
        self.assertIn(seq_record.id, record.accessions)

        # Now try using the iterator - note that all these
        # test cases have only one record.

        # With the SequenceParser
        records = list(SeqIO.parse(datafile, "swiss"))

        self.assertEqual(len(records), 1)
        self.assertIsInstance(records[0], SeqRecord)

        # Check matches what we got earlier without the iterator:
        self.assertEqual(str(records[0].seq), str(seq_record.seq))
        self.assertEqual(records[0].description, seq_record.description)
        self.assertEqual(records[0].name, seq_record.name)
        self.assertEqual(records[0].id, seq_record.id)

        # With the RecordParser
        with open(datafile) as test_handle:
            records = list(SwissProt.parse(test_handle))

        self.assertEqual(len(records), 1)
        self.assertIsInstance(records[0], SwissProt.Record)

        # Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].sequence, record.sequence)
        self.assertEqual(records[0].description, record.description)
        self.assertEqual(records[0].entry_name, record.entry_name)
        self.assertEqual(records[0].accessions, record.accessions)

    def test_sp015(self):
        """Parsing SwissProt file sp015."""
        filename = "sp015"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)

        with open(datafile) as test_handle:
            seq_record = SeqIO.read(test_handle, "swiss")

        self.assertIsInstance(seq_record, SeqRecord)

        self.assertEqual(seq_record.id, "IPI00383150")
        self.assertEqual(seq_record.name, "IPI00383150.2")
        self.assertEqual(seq_record.description, "")
        self.assertEqual(
            repr(seq_record.seq),
            "Seq('MSFQAPRRLLELAGQSLLRDQALAISVLDELPRELFPRLFVEAFTSRRCEVLKV...TPC')",
        )

        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(record.entry_name, "IPI00383150.2")
        self.assertEqual(record.accessions, ["IPI00383150"])
        self.assertEqual(
            record.organism_classification,
            [
                "Eukaryota",
                "Metazoa",
                "Chordata",
                "Craniata",
                "Vertebrata",
                "Euteleostomi",
                "Mammalia",
                "Eutheria",
                "Primates",
                "Catarrhini",
                "Hominidae",
                "Homo",
            ],
        )
        self.assertEqual(record.seqinfo, (457, 52856, "5C3151AAADBDE232"))

        self.assertEqual(len(record.features), 0)
        self.assertEqual(len(record.references), 0)

        # Check the two parsers agree on the essentials
        self.assertEqual(str(seq_record.seq), record.sequence)
        self.assertEqual(seq_record.description, record.description)
        self.assertEqual(seq_record.name, record.entry_name)
        self.assertIn(seq_record.id, record.accessions)

        # Now try using the iterator - note that all these
        # test cases have only one record.

        # With the SequenceParser
        with open(datafile) as test_handle:
            records = list(SeqIO.parse(test_handle, "swiss"))

        self.assertEqual(len(records), 1)
        self.assertIsInstance(records[0], SeqRecord)

        # Check matches what we got earlier without the iterator:
        self.assertEqual(str(records[0].seq), str(seq_record.seq))
        self.assertEqual(records[0].description, seq_record.description)
        self.assertEqual(records[0].name, seq_record.name)
        self.assertEqual(records[0].id, seq_record.id)

        # With the RecordParser
        with open(datafile) as test_handle:
            records = list(SwissProt.parse(test_handle))

        self.assertEqual(len(records), 1)
        self.assertIsInstance(records[0], SwissProt.Record)

        # Check matches what we got earlier without the iterator:
        self.assertEqual(records[0].sequence, record.sequence)
        self.assertEqual(records[0].description, record.description)
        self.assertEqual(records[0].entry_name, record.entry_name)
        self.assertEqual(records[0].accessions, record.accessions)

    def test_P0CK95(self):
        """Parsing SwissProt file P0CK95.txt."""
        filename = "P0CK95.txt"
        datafile = os.path.join("SwissProt", filename)
        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)
        # Check the simple variant
        self.assertEqual(
            record.features[5].qualifiers["note"],
            "N -> G (in strain: O15:H- / 83/39 /ETEC)",
        )
        # Check a FT where the 2nd line starts with /
        self.assertEqual(
            record.features[6].qualifiers["note"],
            "DGTPLPEFYSE -> EGELPKFFSD (in strain: O15:H- / 83/39 / ETEC)",
        )

    def test_ft_line(self):
        """Parsing SwissProt file O23729, which has a new-style FT line."""
        filename = "O23729.txt"
        datafile = os.path.join("SwissProt", filename)
        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
