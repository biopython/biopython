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
    def test_Q13454(self):
        """Parsing SwissProt file Q13454.txt."""
        filename = "Q13454.txt"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)

        with open(datafile) as test_handle:
            seq_record = SeqIO.read(test_handle, "swiss")

        self.assertIsInstance(seq_record, SeqRecord)

        self.assertEqual(seq_record.id, "Q13454")
        self.assertEqual(seq_record.name, "TUSC3_HUMAN")
        self.assertEqual(
            seq_record.description,
            "RecName: Full=Tumor suppressor candidate 3; AltName: Full=Dolichyl-diphosphooligosaccharide--protein glycosyltransferase subunit TUSC3; Short=Oligosaccharyl transferase subunit TUSC3; AltName: Full=Magnesium uptake/transporter TUSC3; AltName: Full=Protein N33; Flags: Precursor;",
        )
        self.assertEqual(
            repr(seq_record.seq),
            "Seq('MGARGAPSRRRQAGRRLRYLPTGSFPFLLLLLLLCIQLGGGQKKKENLLAEKVE...DFE')",
        )

        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(record.entry_name, "TUSC3_HUMAN")
        self.assertEqual(
            record.accessions,
            ["Q13454", "A8MSM0", "D3DSP2", "Q14911", "Q14912", "Q96FW0"],
        )
        self.assertEqual(record.gene_name, [{"Name": "TUSC3", "Synonyms": ["N33"]}])
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
        self.assertEqual(record.seqinfo, (348, 39676, "16D97CB1E00C5190"))

        self.assertEqual(len(record.features), 32)
        feature = record.features[0]
        self.assertEqual(feature.type, "SIGNAL")
        self.assertEqual(feature.location.start, 0)
        self.assertEqual(feature.location.end, 41)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        self.assertIsNone(feature.id)
        feature = record.features[1]
        self.assertEqual(feature.type, "CHAIN")
        self.assertEqual(feature.location.start, 41)
        self.assertEqual(feature.location.end, 348)
        self.assertEqual(feature.qualifiers["note"], "Tumor suppressor candidate 3")
        self.assertEqual(feature.id, "PRO_0000215300")
        feature = record.features[2]
        self.assertEqual(feature.type, "TOPO_DOM")
        self.assertEqual(feature.location.start, 41)
        self.assertEqual(feature.location.end, 196)
        self.assertEqual(feature.qualifiers["note"], "Lumenal")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        self.assertIsNone(feature.id)
        feature = record.features[3]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 196)
        self.assertEqual(feature.location.end, 217)
        self.assertEqual(feature.qualifiers["note"], "Helical")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        self.assertIsNone(feature.id)
        feature = record.features[4]
        self.assertEqual(feature.type, "TOPO_DOM")
        self.assertEqual(feature.location.start, 217)
        self.assertEqual(feature.location.end, 221)
        self.assertEqual(feature.qualifiers["note"], "Cytoplasmic")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        self.assertIsNone(feature.id)
        feature = record.features[5]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 221)
        self.assertEqual(feature.location.end, 242)
        self.assertEqual(feature.qualifiers["note"], "Helical")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        self.assertIsNone(feature.id)
        feature = record.features[6]
        self.assertEqual(feature.type, "TOPO_DOM")
        self.assertEqual(feature.location.start, 242)
        self.assertEqual(feature.location.end, 276)
        self.assertEqual(feature.qualifiers["note"], "Lumenal")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        self.assertIsNone(feature.id)
        feature = record.features[7]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 276)
        self.assertEqual(feature.location.end, 297)
        self.assertEqual(feature.qualifiers["note"], "Helical")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        self.assertIsNone(feature.id)
        feature = record.features[8]
        self.assertEqual(feature.type, "TOPO_DOM")
        self.assertEqual(feature.location.start, 297)
        self.assertEqual(feature.location.end, 312)
        self.assertEqual(feature.qualifiers["note"], "Cytoplasmic")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        self.assertIsNone(feature.id)
        feature = record.features[9]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 312)
        self.assertEqual(feature.location.end, 333)
        self.assertEqual(feature.qualifiers["note"], "Helical")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        self.assertIsNone(feature.id)
        feature = record.features[10]
        self.assertEqual(feature.type, "TOPO_DOM")
        self.assertEqual(feature.location.start, 333)
        self.assertEqual(feature.location.end, 348)
        self.assertEqual(feature.qualifiers["note"], "Lumenal")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        self.assertIsNone(feature.id)
        feature = record.features[11]
        self.assertEqual(feature.type, "DOMAIN")
        self.assertEqual(feature.location.start, 58)
        self.assertEqual(feature.location.end, 187)
        self.assertEqual(feature.qualifiers["note"], "Thioredoxin")
        self.assertIsNone(feature.id)
        feature = record.features[12]
        self.assertEqual(feature.type, "DISULFID")
        self.assertEqual(feature.location.start, 98)
        self.assertEqual(feature.location.end, 102)
        self.assertEqual(feature.qualifiers["note"], "Redox-active")
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:24685145,ECO:0007744|PDB:4M8G, ECO:0007744|PDB:4M90",
        )
        self.assertIsNone(feature.id)
        feature = record.features[13]
        self.assertEqual(feature.type, "VAR_SEQ")
        self.assertEqual(feature.location.start, 343)
        self.assertEqual(feature.location.end, 348)
        self.assertEqual(feature.qualifiers["note"], "DLDFE -> FLIK (in isoform 2)")
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000303|PubMed:15489334,ECO:0000303|PubMed:8661104",
        )
        self.assertEqual(feature.id, "VSP_003776")
        feature = record.features[14]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 64)
        self.assertEqual(feature.location.end, 65)
        self.assertEqual(feature.qualifiers["note"], "I -> V (in dbSNP:rs11545035)")
        self.assertEqual(feature.id, "VAR_045836")
        feature = record.features[15]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 246)
        self.assertEqual(feature.location.end, 247)
        self.assertEqual(feature.qualifiers["note"], "M -> V")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000269|PubMed:23033978")
        self.assertEqual(feature.id, "VAR_069369")
        feature = record.features[16]
        self.assertEqual(feature.type, "MUTAGEN")
        self.assertEqual(feature.location.start, 98)
        self.assertEqual(feature.location.end, 99)
        self.assertEqual(
            feature.qualifiers["note"],
            "C->S: Reduces N-glycosylation of cysteine-proximal acceptor sites; when associated with S-102.",
        )
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000269|PubMed:25135935")
        self.assertIsNone(feature.id)
        feature = record.features[17]
        self.assertEqual(feature.type, "MUTAGEN")
        self.assertEqual(feature.location.start, 101)
        self.assertEqual(feature.location.end, 102)
        self.assertEqual(
            feature.qualifiers["note"],
            "C->S: Reduces N-glycosylation of cysteine-proximal acceptor sites; when associated with S-99.",
        )
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000269|PubMed:25135935")
        self.assertIsNone(feature.id)
        feature = record.features[18]
        self.assertEqual(feature.type, "HELIX")
        self.assertEqual(feature.location.start, 43)
        self.assertEqual(feature.location.end, 62)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:4M91")
        self.assertIsNone(feature.id)
        feature = record.features[19]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 63)
        self.assertEqual(feature.location.end, 67)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:4M91")
        self.assertIsNone(feature.id)
        feature = record.features[20]
        self.assertEqual(feature.type, "HELIX")
        self.assertEqual(feature.location.start, 69)
        self.assertEqual(feature.location.end, 76)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:4M91")
        self.assertIsNone(feature.id)
        feature = record.features[21]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 83)
        self.assertEqual(feature.location.end, 91)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:4M91")
        self.assertIsNone(feature.id)
        feature = record.features[22]
        self.assertEqual(feature.type, "HELIX")
        self.assertEqual(feature.location.start, 94)
        self.assertEqual(feature.location.end, 97)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:4M91")
        self.assertIsNone(feature.id)
        feature = record.features[23]
        self.assertEqual(feature.type, "HELIX")
        self.assertEqual(feature.location.start, 99)
        self.assertEqual(feature.location.end, 118)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:4M91")
        self.assertIsNone(feature.id)
        feature = record.features[24]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 125)
        self.assertEqual(feature.location.end, 132)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:4M91")
        self.assertIsNone(feature.id)
        feature = record.features[25]
        self.assertEqual(feature.type, "TURN")
        self.assertEqual(feature.location.start, 132)
        self.assertEqual(feature.location.end, 135)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:4M91")
        self.assertIsNone(feature.id)
        feature = record.features[26]
        self.assertEqual(feature.type, "HELIX")
        self.assertEqual(feature.location.start, 136)
        self.assertEqual(feature.location.end, 142)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:4M91")
        self.assertIsNone(feature.id)
        feature = record.features[27]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 149)
        self.assertEqual(feature.location.end, 154)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:4M91")
        self.assertIsNone(feature.id)
        feature = record.features[28]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 155)
        self.assertEqual(feature.location.end, 158)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:4M91")
        self.assertIsNone(feature.id)
        feature = record.features[29]
        self.assertEqual(feature.type, "HELIX")
        self.assertEqual(feature.location.start, 161)
        self.assertEqual(feature.location.end, 164)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:4M91")
        self.assertIsNone(feature.id)
        feature = record.features[30]
        self.assertEqual(feature.type, "HELIX")
        self.assertEqual(feature.location.start, 167)
        self.assertEqual(feature.location.end, 171)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:4M91")
        self.assertIsNone(feature.id)
        feature = record.features[31]
        self.assertEqual(feature.type, "HELIX")
        self.assertEqual(feature.location.start, 174)
        self.assertEqual(feature.location.end, 186)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:4M91")
        self.assertIsNone(feature.id)
        self.assertEqual(len(record.references), 12)
        self.assertEqual(
            record.references[0].authors,
            "Macgrogan D., Levy A., Bova G.S., Isaacs W.B., Bookstein R.",
        )
        self.assertEqual(
            record.references[0].title,
            "Structure and methylation-associated silencing of a gene within a homozygously deleted region of human chromosome band 8p22.",
        )
        self.assertEqual(
            record.references[0].positions,
            ["NUCLEOTIDE SEQUENCE [GENOMIC DNA / MRNA] (ISOFORMS 1 AND 2)."],
        )
        self.assertEqual(len(record.references[0].references), 2)
        self.assertEqual(record.references[0].references[0], ("PubMed", "8661104"))
        self.assertEqual(
            record.references[0].references[1], ("DOI", "10.1006/geno.1996.0322")
        )
        self.assertEqual(
            record.references[1].authors,
            "Kalnine N., Chen X., Rolfs A., Halleck A., Hines L., Eisenstein S., Koundinya M., Raphael J., Moreira D., Kelley T., LaBaer J., Lin Y., Phelan M., Farmer A.",
        )
        self.assertEqual(
            record.references[1].title,
            "Cloning of human full-length CDSs in BD Creator(TM) system donor vector.",
        )
        self.assertEqual(
            record.references[1].positions,
            ["NUCLEOTIDE SEQUENCE [LARGE SCALE MRNA] (ISOFORM 1)."],
        )
        self.assertEqual(len(record.references[1].references), 0)
        self.assertEqual(
            record.references[2].authors,
            "Nusbaum C., Mikkelsen T.S., Zody M.C., Asakawa S., Taudien S., Garber M., Kodira C.D., Schueler M.G., Shimizu A., Whittaker C.A., Chang J.L., Cuomo C.A., Dewar K., FitzGerald M.G., Yang X., Allen N.R., Anderson S., Asakawa T., Blechschmidt K., Bloom T., Borowsky M.L., Butler J., Cook A., Corum B., DeArellano K., DeCaprio D., Dooley K.T., Dorris L. III, Engels R., Gloeckner G., Hafez N., Hagopian D.S., Hall J.L., Ishikawa S.K., Jaffe D.B., Kamat A., Kudoh J., Lehmann R., Lokitsang T., Macdonald P., Major J.E., Matthews C.D., Mauceli E., Menzel U., Mihalev A.H., Minoshima S., Murayama Y., Naylor J.W., Nicol R., Nguyen C., O'Leary S.B., O'Neill K., Parker S.C.J., Polley A., Raymond C.K., Reichwald K., Rodriguez J., Sasaki T., Schilhabel M., Siddiqui R., Smith C.L., Sneddon T.P., Talamas J.A., Tenzin P., Topham K., Venkataraman V., Wen G., Yamazaki S., Young S.K., Zeng Q., Zimmer A.R., Rosenthal A., Birren B.W., Platzer M., Shimizu N., Lander E.S.",
        )
        self.assertEqual(
            record.references[2].title,
            "DNA sequence and analysis of human chromosome 8.",
        )
        self.assertEqual(
            record.references[2].positions,
            ["NUCLEOTIDE SEQUENCE [LARGE SCALE GENOMIC DNA]."],
        )
        self.assertEqual(len(record.references[2].references), 2)
        self.assertEqual(record.references[2].references[0], ("PubMed", "16421571"))
        self.assertEqual(
            record.references[2].references[1], ("DOI", "10.1038/nature04406")
        )
        self.assertEqual(
            record.references[3].authors,
            "Mural R.J., Istrail S., Sutton G.G., Florea L., Halpern A.L., Mobarry C.M., Lippert R., Walenz B., Shatkay H., Dew I., Miller J.R., Flanigan M.J., Edwards N.J., Bolanos R., Fasulo D., Halldorsson B.V., Hannenhalli S., Turner R., Yooseph S., Lu F., Nusskern D.R., Shue B.C., Zheng X.H., Zhong F., Delcher A.L., Huson D.H., Kravitz S.A., Mouchard L., Reinert K., Remington K.A., Clark A.G., Waterman M.S., Eichler E.E., Adams M.D., Hunkapiller M.W., Myers E.W., Venter J.C.",
        )
        self.assertEqual(
            record.references[3].title,
            "",
        )
        self.assertEqual(
            record.references[3].positions,
            ["NUCLEOTIDE SEQUENCE [LARGE SCALE GENOMIC DNA]."],
        )
        self.assertEqual(len(record.references[3].references), 0)
        self.assertEqual(
            record.references[4].authors,
            "The MGC Project Team",
        )
        self.assertEqual(
            record.references[4].title,
            "The status, quality, and expansion of the NIH full-length cDNA project: the Mammalian Gene Collection (MGC).",
        )
        self.assertEqual(
            record.references[4].positions,
            ["NUCLEOTIDE SEQUENCE [LARGE SCALE MRNA] (ISOFORM 2)."],
        )
        self.assertEqual(len(record.references[4].references), 2)
        self.assertEqual(record.references[4].references[0], ("PubMed", "15489334"))
        self.assertEqual(
            record.references[4].references[1], ("DOI", "10.1101/gr.2596504")
        )
        self.assertEqual(
            record.references[5].authors,
            "Kelleher D.J., Karaoglu D., Mandon E.C., Gilmore R.",
        )
        self.assertEqual(
            record.references[5].title,
            "Oligosaccharyltransferase isoforms that contain different catalytic STT3 subunits have distinct enzymatic properties.",
        )
        self.assertEqual(
            record.references[5].positions,
            [
                "IDENTIFICATION IN THE OLIGOSACCHARYLTRANSFERASE (OST) COMPLEX, AND TISSUE",
                "SPECIFICITY.",
            ],
        )
        self.assertEqual(len(record.references[5].references), 2)
        self.assertEqual(record.references[5].references[0], ("PubMed", "12887896"))
        self.assertEqual(
            record.references[5].references[1], ("DOI", "10.1016/s1097-2765(03)00243-0")
        )
        self.assertEqual(
            record.references[6].authors,
            "Molinari F., Foulquier F., Tarpey P.S., Morelle W., Boissel S., Teague J., Edkins S., Futreal P.A., Stratton M.R., Turner G., Matthijs G., Gecz J., Munnich A., Colleaux L.",
        )
        self.assertEqual(
            record.references[6].title,
            "Oligosaccharyltransferase-subunit mutations in nonsyndromic mental retardation.",
        )
        self.assertEqual(record.references[6].positions, ["INVOLVEMENT IN MRT7."])
        self.assertEqual(len(record.references[6].references), 2)
        self.assertEqual(record.references[6].references[0], ("PubMed", "18455129"))
        self.assertEqual(
            record.references[6].references[1], ("DOI", "10.1016/j.ajhg.2008.03.021")
        )
        self.assertEqual(
            record.references[7].authors,
            "Garshasbi M., Hadavi V., Habibi H., Kahrizi K., Kariminejad R., Behjati F., Tzschach A., Najmabadi H., Ropers H.H., Kuss A.W.",
        )
        self.assertEqual(
            record.references[7].title,
            "A defect in the TUSC3 gene is associated with autosomal recessive mental retardation.",
        )
        self.assertEqual(record.references[7].positions, ["INVOLVEMENT IN MRT7."])
        self.assertEqual(len(record.references[7].references), 2)
        self.assertEqual(record.references[7].references[0], ("PubMed", "18452889"))
        self.assertEqual(
            record.references[7].references[1], ("DOI", "10.1016/j.ajhg.2008.03.018")
        )
        self.assertEqual(
            record.references[8].authors,
            "Zhou H., Clapham D.E.",
        )
        self.assertEqual(
            record.references[8].title,
            "Mammalian MagT1 and TUSC3 are required for cellular magnesium uptake and vertebrate embryonic development.",
        )
        self.assertEqual(
            record.references[8].positions, ["FUNCTION IN MAGNESIUM UPTAKE."]
        )
        self.assertEqual(len(record.references[8].references), 2)
        self.assertEqual(record.references[8].references[0], ("PubMed", "19717468"))
        self.assertEqual(
            record.references[8].references[1], ("DOI", "10.1073/pnas.0908332106")
        )
        self.assertEqual(
            record.references[9].authors,
            "Cherepanova N.A., Shrimal S., Gilmore R.",
        )
        self.assertEqual(
            record.references[9].title,
            "Oxidoreductase activity is necessary for N-glycosylation of cysteine-proximal acceptor sites in glycoproteins.",
        )
        self.assertEqual(
            record.references[9].positions,
            ["FUNCTION, AND MUTAGENESIS OF CYS-99 AND CYS-102."],
        )
        self.assertEqual(len(record.references[9].references), 2)
        self.assertEqual(record.references[9].references[0], ("PubMed", "25135935"))
        self.assertEqual(
            record.references[9].references[1], ("DOI", "10.1083/jcb.201404083")
        )
        self.assertEqual(
            record.references[10].authors,
            "Mohorko E., Owen R.L., Malojcic G., Brozzo M.S., Aebi M., Glockshuber R.",
        )
        self.assertEqual(
            record.references[10].title,
            "Structural basis of substrate specificity of human oligosaccharyl transferase subunit N33/Tusc3 and its role in regulating protein N-glycosylation.",
        )
        self.assertEqual(
            record.references[10].positions,
            [
                "X-RAY CRYSTALLOGRAPHY (1.10 ANGSTROMS) OF 44-194, DISULFIDE BOND, PROPOSED",
                "FUNCTION, AND SUBUNIT.",
            ],
        )
        self.assertEqual(len(record.references[10].references), 2)
        self.assertEqual(record.references[10].references[0], ("PubMed", "24685145"))
        self.assertEqual(
            record.references[10].references[1], ("DOI", "10.1016/j.str.2014.02.013")
        )
        self.assertEqual(
            record.references[11].authors,
            "de Ligt J., Willemsen M.H., van Bon B.W., Kleefstra T., Yntema H.G., Kroes T., Vulto-van Silfhout A.T., Koolen D.A., de Vries P., Gilissen C., del Rosario M., Hoischen A., Scheffer H., de Vries B.B., Brunner H.G., Veltman J.A., Vissers L.E.",
        )
        self.assertEqual(
            record.references[11].title,
            "Diagnostic exome sequencing in persons with severe intellectual disability.",
        )
        self.assertEqual(record.references[11].positions, ["VARIANT VAL-247."])
        self.assertEqual(len(record.references[11].references), 2)
        self.assertEqual(record.references[11].references[0], ("PubMed", "23033978"))
        self.assertEqual(
            record.references[11].references[1], ("DOI", "10.1056/nejmoa1206524")
        )

        # Check that the two parsers agree on the essentials
        self.assertEqual(seq_record.seq, record.sequence)
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
        self.assertEqual(records[0].seq, seq_record.seq)
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

    def test_P60904(self):
        """Parsing SwissProt file P60904.txt."""
        filename = "P60904.txt"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)

        with open(datafile) as test_handle:
            seq_record = SeqIO.read(test_handle, "swiss")

        self.assertIsInstance(seq_record, SeqRecord)

        self.assertEqual(seq_record.id, "P60904")
        self.assertEqual(seq_record.name, "DNJC5_MOUSE")
        self.assertEqual(
            seq_record.description,
            "RecName: Full=DnaJ homolog subfamily C member 5 {ECO:0000305}; AltName: Full=Cysteine string protein {ECO:0000250|UniProtKB:Q9H3Z4}; Short=CSP {ECO:0000250|UniProtKB:Q9H3Z4}; AltName: Full=Cysteine-string protein isoform alpha {ECO:0000303|PubMed:20847230};",
        )
        self.assertEqual(
            repr(seq_record.seq),
            "Seq('MADQRQRSLSTSGESLYHVLGLDKNATSDDIKKSYRKLALKYHPDKNPDNPEAA...GFN')",
        )

        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(record.entry_name, "DNJC5_MOUSE")
        self.assertEqual(record.accessions, ["P60904", "P54101"])
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
                "Glires",
                "Rodentia",
                "Myomorpha",
                "Muroidea",
                "Muridae",
                "Murinae",
                "Mus",
                "Mus",
            ],
        )
        self.assertEqual(record.seqinfo, (198, 22101, "52F98261FBAD978F"))

        self.assertEqual(len(record.features), 23)
        feature = record.features[0]
        self.assertEqual(feature.type, "CHAIN")
        self.assertEqual(feature.location.start, 0)
        self.assertEqual(feature.location.end, 198)
        self.assertEqual(
            feature.qualifiers["note"], "DnaJ homolog subfamily C member 5"
        )
        feature = record.features[1]
        self.assertEqual(feature.type, "DOMAIN")
        self.assertEqual(feature.location.start, 12)
        self.assertEqual(feature.location.end, 82)
        self.assertEqual(feature.qualifiers["note"], "J")
        self.assertEqual(
            feature.qualifiers["evidence"], "ECO:0000255|PROSITE-ProRule:PRU00286"
        )
        feature = record.features[2]
        self.assertEqual(feature.type, "MOD_RES")
        self.assertEqual(feature.location.start, 7)
        self.assertEqual(feature.location.end, 8)
        self.assertEqual(feature.qualifiers["note"], "Phosphoserine")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000250|UniProtKB:Q9H3Z4")
        feature = record.features[3]
        self.assertEqual(feature.type, "MOD_RES")
        self.assertEqual(feature.location.start, 9)
        self.assertEqual(feature.location.end, 10)
        self.assertEqual(feature.qualifiers["note"], "Phosphoserine")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007744|PubMed:19131326")
        feature = record.features[4]
        self.assertEqual(feature.type, "MOD_RES")
        self.assertEqual(feature.location.start, 11)
        self.assertEqual(feature.location.end, 12)
        self.assertEqual(feature.qualifiers["note"], "Phosphoserine")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000250|UniProtKB:Q9H3Z4")
        feature = record.features[5]
        self.assertEqual(feature.type, "MOD_RES")
        self.assertEqual(feature.location.start, 14)
        self.assertEqual(feature.location.end, 15)
        self.assertEqual(feature.qualifiers["note"], "Phosphoserine")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007744|PubMed:21183079")
        feature = record.features[6]
        self.assertEqual(feature.type, "MOD_RES")
        self.assertEqual(feature.location.start, 16)
        self.assertEqual(feature.location.end, 17)
        self.assertEqual(feature.qualifiers["note"], "Phosphotyrosine")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000250|UniProtKB:Q9H3Z4")
        feature = record.features[7]
        self.assertEqual(feature.type, "MOD_RES")
        self.assertEqual(feature.location.start, 55)
        self.assertEqual(feature.location.end, 56)
        self.assertEqual(feature.qualifiers["note"], "N6-acetyllysine")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000250|UniProtKB:Q9H3Z4")
        feature = record.features[8]
        self.assertEqual(feature.type, "MOD_RES")
        self.assertEqual(feature.location.start, 150)
        self.assertEqual(feature.location.end, 151)
        self.assertEqual(feature.qualifiers["note"], "Phosphoserine")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007744|PubMed:21183079")
        feature = record.features[9]
        self.assertEqual(feature.type, "MUTAGEN")
        self.assertEqual(feature.location.start, 9)
        self.assertEqual(feature.location.end, 10)
        self.assertEqual(
            feature.qualifiers["note"],
            "S->D: Reduced interaction with SYT9, but no effect on the interaction with HSC70.",
        )
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000269|PubMed:20847230")
        feature = record.features[10]
        self.assertEqual(feature.type, "MUTAGEN")
        self.assertEqual(feature.location.start, 92)
        self.assertEqual(feature.location.end, 93)
        self.assertEqual(
            feature.qualifiers["note"], "E->V: Reduced interaction with SYT9."
        )
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000269|PubMed:20847230")
        feature = record.features[11]
        self.assertEqual(feature.type, "MUTAGEN")
        self.assertEqual(feature.location.start, 112)
        self.assertEqual(feature.location.end, 113)
        self.assertEqual(
            feature.qualifiers["note"],
            "C->V: No effect on palmitoylation. No change in subcellular location; when associated with G-118 and F-121.",
        )
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000269|PubMed:17034881")
        feature = record.features[12]
        self.assertEqual(feature.type, "MUTAGEN")
        self.assertEqual(feature.location.start, 117)
        self.assertEqual(feature.location.end, 118)
        self.assertEqual(
            feature.qualifiers["note"],
            "C->G: No effect on palmitoylation. No change in subcellular location; when associated with V-113 and F-121.",
        )
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000269|PubMed:17034881")
        feature = record.features[13]
        self.assertEqual(feature.type, "MUTAGEN")
        self.assertEqual(feature.location.start, 120)
        self.assertEqual(feature.location.end, 121)
        self.assertEqual(
            feature.qualifiers["note"],
            "C->F: No effect on palmitoylation. No change in subcellular location; when associated with V-113 and G-118.",
        )
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000269|PubMed:17034881")
        feature = record.features[14]
        self.assertEqual(feature.type, "MUTAGEN")
        self.assertEqual(feature.location.start, 128)
        self.assertEqual(feature.location.end, 129)
        self.assertEqual(
            feature.qualifiers["note"],
            "F->C: No effect on palmitoylation. No change in subcellular location; when associated with H-135.",
        )
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000269|PubMed:17034881")
        feature = record.features[15]
        self.assertEqual(feature.type, "MUTAGEN")
        self.assertEqual(feature.location.start, 134)
        self.assertEqual(feature.location.end, 135)
        self.assertEqual(
            feature.qualifiers["note"],
            "K->H: No effect on palmitoylation. No change in subcellular location; when associated with C-129.",
        )
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000269|PubMed:17034881")
        feature = record.features[16]
        self.assertEqual(feature.type, "HELIX")
        self.assertEqual(feature.location.start, 15)
        self.assertEqual(feature.location.end, 20)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:2CTW")
        feature = record.features[17]
        self.assertEqual(feature.type, "HELIX")
        self.assertEqual(feature.location.start, 27)
        self.assertEqual(feature.location.end, 41)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:2CTW")
        feature = record.features[18]
        self.assertEqual(feature.type, "TURN")
        self.assertEqual(feature.location.start, 43)
        self.assertEqual(feature.location.end, 46)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:2CTW")
        feature = record.features[19]
        self.assertEqual(feature.type, "HELIX")
        self.assertEqual(feature.location.start, 50)
        self.assertEqual(feature.location.end, 67)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:2CTW")
        feature = record.features[20]
        self.assertEqual(feature.type, "HELIX")
        self.assertEqual(feature.location.start, 69)
        self.assertEqual(feature.location.end, 78)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:2CTW")
        feature = record.features[21]
        self.assertEqual(feature.type, "HELIX")
        self.assertEqual(feature.location.start, 80)
        self.assertEqual(feature.location.end, 89)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:2CTW")
        feature = record.features[22]
        self.assertEqual(feature.type, "HELIX")
        self.assertEqual(feature.location.start, 93)
        self.assertEqual(feature.location.end, 100)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:2CTW")
        self.assertIsNone(feature.id)
        self.assertEqual(len(record.references), 16)
        reference = record.references[0]
        self.assertEqual(reference.authors, "Qin N., Lin T., Birnbaumer L.")
        self.assertEqual(reference.title, "")
        self.assertEqual(len(reference.references), 0)
        reference = record.references[1]
        self.assertEqual(
            reference.authors,
            "Carninci P., Kasukawa T., Katayama S., Gough J., Frith M.C., Maeda N., Oyama R., Ravasi T., Lenhard B., Wells C., Kodzius R., Shimokawa K., Bajic V.B., Brenner S.E., Batalov S., Forrest A.R., Zavolan M., Davis M.J., Wilming L.G., Aidinis V., Allen J.E., Ambesi-Impiombato A., Apweiler R., Aturaliya R.N., Bailey T.L., Bansal M., Baxter L., Beisel K.W., Bersano T., Bono H., Chalk A.M., Chiu K.P., Choudhary V., Christoffels A., Clutterbuck D.R., Crowe M.L., Dalla E., Dalrymple B.P., de Bono B., Della Gatta G., di Bernardo D., Down T., Engstrom P., Fagiolini M., Faulkner G., Fletcher C.F., Fukushima T., Furuno M., Futaki S., Gariboldi M., Georgii-Hemming P., Gingeras T.R., Gojobori T., Green R.E., Gustincich S., Harbers M., Hayashi Y., Hensch T.K., Hirokawa N., Hill D., Huminiecki L., Iacono M., Ikeo K., Iwama A., Ishikawa T., Jakt M., Kanapin A., Katoh M., Kawasawa Y., Kelso J., Kitamura H., Kitano H., Kollias G., Krishnan S.P., Kruger A., Kummerfeld S.K., Kurochkin I.V., Lareau L.F., Lazarevic D., Lipovich L., Liu J., Liuni S., McWilliam S., Madan Babu M., Madera M., Marchionni L., Matsuda H., Matsuzawa S., Miki H., Mignone F., Miyake S., Morris K., Mottagui-Tabar S., Mulder N., Nakano N., Nakauchi H., Ng P., Nilsson R., Nishiguchi S., Nishikawa S., Nori F., Ohara O., Okazaki Y., Orlando V., Pang K.C., Pavan W.J., Pavesi G., Pesole G., Petrovsky N., Piazza S., Reed J., Reid J.F., Ring B.Z., Ringwald M., Rost B., Ruan Y., Salzberg S.L., Sandelin A., Schneider C., Schoenbach C., Sekiguchi K., Semple C.A., Seno S., Sessa L., Sheng Y., Shibata Y., Shimada H., Shimada K., Silva D., Sinclair B., Sperling S., Stupka E., Sugiura K., Sultana R., Takenaka Y., Taki K., Tammoja K., Tan S.L., Tang S., Taylor M.S., Tegner J., Teichmann S.A., Ueda H.R., van Nimwegen E., Verardo R., Wei C.L., Yagi K., Yamanishi H., Zabarovsky E., Zhu S., Zimmer A., Hide W., Bult C., Grimmond S.M., Teasdale R.D., Liu E.T., Brusic V., Quackenbush J., Wahlestedt C., Mattick J.S., Hume D.A., Kai C., Sasaki D., Tomaru Y., Fukuda S., Kanamori-Katayama M., Suzuki M., Aoki J., Arakawa T., Iida J., Imamura K., Itoh M., Kato T., Kawaji H., Kawagashira N., Kawashima T., Kojima M., Kondo S., Konno H., Nakano K., Ninomiya N., Nishio T., Okada M., Plessy C., Shibata K., Shiraki T., Suzuki S., Tagami M., Waki K., Watahiki A., Okamura-Oho Y., Suzuki H., Kawai J., Hayashizaki Y.",
        )
        self.assertEqual(
            reference.title, "The transcriptional landscape of the mammalian genome."
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "16141072"))
        self.assertEqual(reference.references[1], ("DOI", "10.1126/science.1112014"))
        reference = record.references[2]
        self.assertEqual(reference.authors, "Lubec G., Kang S.U.")
        self.assertEqual(reference.title, "")
        self.assertEqual(len(reference.references), 0)
        reference = record.references[3]
        self.assertEqual(
            reference.authors,
            "Boal F., Le Pevelen S., Cziepluch C., Scotti P., Lang J.",
        )
        self.assertEqual(
            reference.title,
            "Cysteine-string protein isoform beta (Cspbeta) is targeted to the trans-Golgi network as a non-palmitoylated CSP in clonal beta-cells.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "17034881"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1016/j.bbamcr.2006.08.054")
        )
        reference = record.references[4]
        self.assertEqual(
            reference.authors,
            "Lee J., Xu Y., Chen Y., Sprung R., Kim S.C., Xie S., Zhao Y.",
        )
        self.assertEqual(
            reference.title,
            "Mitochondrial phosphoproteome revealed by an improved IMAC method and MS/MS/MS.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "17208939"))
        self.assertEqual(reference.references[1], ("DOI", "10.1074/mcp.m600218-mcp200"))
        reference = record.references[5]
        self.assertEqual(
            reference.authors, "Villen J., Beausoleil S.A., Gerber S.A., Gygi S.P."
        )
        self.assertEqual(
            reference.title, "Large-scale phosphorylation analysis of mouse liver."
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "17242355"))
        self.assertEqual(reference.references[1], ("DOI", "10.1073/pnas.0609836104"))
        reference = record.references[6]
        self.assertEqual(
            reference.authors,
            "Zhou H., Ye M., Dong J., Han G., Jiang X., Wu R., Zou H.",
        )
        self.assertEqual(
            reference.title,
            "Specific phosphopeptide enrichment with immobilized titanium ion affinity chromatography adsorbent for phosphoproteome analysis.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "18630941"))
        self.assertEqual(reference.references[1], ("DOI", "10.1021/pr800223m"))
        reference = record.references[7]
        self.assertEqual(
            reference.authors,
            "Trost M., English L., Lemieux S., Courcelles M., Desjardins M., Thibault P.",
        )
        self.assertEqual(
            reference.title,
            "The phagosomal proteome in interferon-gamma-activated macrophages.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "19144319"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1016/j.immuni.2008.11.006")
        )
        reference = record.references[8]
        self.assertEqual(
            reference.authors,
            "Sweet S.M., Bailey C.M., Cunningham D.L., Heath J.K., Cooper H.J.",
        )
        self.assertEqual(
            reference.title,
            "Large scale localization of protein phosphorylation by use of electron capture dissociation mass spectrometry.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "19131326"))
        self.assertEqual(reference.references[1], ("DOI", "10.1074/mcp.m800451-mcp200"))
        reference = record.references[9]
        self.assertEqual(
            reference.authors,
            "Huttlin E.L., Jedrychowski M.P., Elias J.E., Goswami T., Rad R., Beausoleil S.A., Villen J., Haas W., Sowa M.E., Gygi S.P.",
        )
        self.assertEqual(
            reference.title,
            "A tissue-specific atlas of mouse protein phosphorylation and expression.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "21183079"))
        self.assertEqual(reference.references[1], ("DOI", "10.1016/j.cell.2010.12.001"))
        reference = record.references[10]
        self.assertEqual(
            reference.authors, "Boal F., Laguerre M., Milochau A., Lang J., Scotti P.A."
        )
        self.assertEqual(
            reference.title,
            "A charged prominence in the linker domain of the cysteine-string protein Cspalpha mediates its regulated interaction with the calcium sensor synaptotagmin 9 during exocytosis.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "20847230"))
        self.assertEqual(reference.references[1], ("DOI", "10.1096/fj.09-152033"))
        reference = record.references[11]
        self.assertEqual(
            reference.authors,
            "Sharma M., Burre J., Bronk P., Zhang Y., Xu W., Suedhof T.C.",
        )
        self.assertEqual(
            reference.title,
            "CSPalpha knockout causes neurodegeneration by impairing SNAP-25 function.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "22187053"))
        self.assertEqual(reference.references[1], ("DOI", "10.1038/emboj.2011.467"))
        reference = record.references[12]
        self.assertEqual(
            reference.authors,
            "Lemonidis K., Gorleku O.A., Sanchez-Perez M.C., Grefen C., Chamberlain L.H.",
        )
        self.assertEqual(
            reference.title,
            "The Golgi S-acylation machinery comprises zDHHC enzymes with major differences in substrate affinity and S-acylation activity.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "25253725"))
        self.assertEqual(reference.references[1], ("DOI", "10.1091/mbc.e14-06-1169"))
        reference = record.references[13]
        self.assertEqual(
            reference.authors, "Lemonidis K., Sanchez-Perez M.C., Chamberlain L.H."
        )
        self.assertEqual(
            reference.title,
            "Identification of a novel sequence motif recognized by the ankyrin repeat domain of zDHHC17/13 S-acyltransferases.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "26198635"))
        self.assertEqual(reference.references[1], ("DOI", "10.1074/jbc.m115.657668"))
        reference = record.references[14]
        self.assertEqual(reference.authors, "Burgoyne R.D., Morgan A.")
        self.assertEqual(
            reference.title,
            "Cysteine string protein (CSP) and its role in preventing neurodegeneration.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "25800794"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1016/j.semcdb.2015.03.008")
        )
        reference = record.references[15]
        self.assertEqual(
            reference.authors, "RIKEN structural genomics initiative (RSGI)"
        )
        self.assertEqual(
            reference.title,
            "Solution structure of J-domain from mouse DnaJ subfamily C member 5.",
        )
        self.assertEqual(len(reference.references), 0)

        # Check the two parsers agree on the essentials
        self.assertEqual(seq_record.seq, record.sequence)
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
        self.assertEqual(records[0].seq, seq_record.seq)
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

    def test_P62258(self):
        """Parsing SwissProt file P62258."""
        filename = "P62258.txt"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)

        with open(datafile) as test_handle:
            seq_record = SeqIO.read(test_handle, "swiss")

        self.assertIsInstance(seq_record, SeqRecord)

        self.assertEqual(seq_record.id, "P62258")
        self.assertEqual(seq_record.name, "1433E_HUMAN")
        self.assertEqual(
            seq_record.description,
            "RecName: Full=14-3-3 protein epsilon; Short=14-3-3E;",
        )
        self.assertEqual(
            repr(seq_record.seq),
            "Seq('MDDREDLVYQAKLAEQAERYDEMVESMKKVAGMDVELTVEERNLLSVAYKNVIG...ENQ')",
        )

        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(record.entry_name, "1433E_HUMAN")
        self.assertEqual(
            record.accessions,
            [
                "P62258",
                "B3KY71",
                "D3DTH5",
                "P29360",
                "P42655",
                "Q4VJB6",
                "Q53XZ5",
                "Q63631",
                "Q7M4R4",
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
        self.assertEqual(record.seqinfo, (255, 29174, "07817CCBD1F75B26"))

        self.assertEqual(len(record.features), 32)
        self.assertEqual(len(record.references), 46)
        reference = record.references[0]
        self.assertEqual(reference.authors, "Conklin D.S., Galaktionov K., Beach D.")
        self.assertEqual(
            reference.title, "14-3-3 proteins associate with cdc25 phosphatases."
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "7644510"))
        self.assertEqual(reference.references[1], ("DOI", "10.1073/pnas.92.17.7892"))
        reference = record.references[1]
        self.assertEqual(
            reference.authors, "Chong S.S., Tanigami A., Roschke A.V., Ledbetter D.H."
        )
        self.assertEqual(
            reference.title,
            "14-3-3 epsilon has no homology to LIS1 and lies telomeric to it on chromosome 17p13.3 outside the Miller-Dieker syndrome chromosome region.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "8858348"))
        self.assertEqual(reference.references[1], ("DOI", "10.1101/gr.6.8.735"))
        reference = record.references[2]
        self.assertEqual(
            reference.authors, "Jin D.-Y., Lyu M.S., Kozak C.A., Jeang K.-T."
        )
        self.assertEqual(reference.title, "Function of 14-3-3 proteins.")
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "8684458"))
        self.assertEqual(reference.references[1], ("DOI", "10.1038/382308a0"))
        reference = record.references[3]
        self.assertEqual(
            reference.authors,
            "Han D., Ye G., Liu T., Chen C., Yang X., Wan B., Pan Y., Yu L.",
        )
        self.assertEqual(
            reference.title,
            "Functional identification of a novel 14-3-3 epsilon splicing variant suggests dimerization is not necessary for 14-3-3 epsilon to inhibit UV-induced apoptosis.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "20417184"))
        self.assertEqual(reference.references[1], ("DOI", "10.1016/j.bbrc.2010.04.104"))
        reference = record.references[4]
        self.assertEqual(reference.authors, "Luk S.C.W., Lee C.Y., Waye M.M.Y.")
        self.assertEqual(
            reference.title, "Sequence determination of human epsilon 14-3-3 protein."
        )
        self.assertEqual(len(reference.references), 0)
        reference = record.references[5]
        self.assertEqual(reference.authors, "Tanigami A., Chong S.S., Ledbetter D.H.")
        self.assertEqual(reference.title, "14-3-3 epsilon genomic sequence.")
        self.assertEqual(len(reference.references), 0)
        reference = record.references[6]
        self.assertEqual(
            reference.authors,
            "Ota T., Suzuki Y., Nishikawa T., Otsuki T., Sugiyama T., Irie R., Wakamatsu A., Hayashi K., Sato H., Nagai K., Kimura K., Makita H., Sekine M., Obayashi M., Nishi T., Shibahara T., Tanaka T., Ishii S., Yamamoto J., Saito K., Kawai Y., Isono Y., Nakamura Y., Nagahari K., Murakami K., Yasuda T., Iwayanagi T., Wagatsuma M., Shiratori A., Sudo H., Hosoiri T., Kaku Y., Kodaira H., Kondo H., Sugawara M., Takahashi M., Kanda K., Yokoi T., Furuya T., Kikkawa E., Omura Y., Abe K., Kamihara K., Katsuta N., Sato K., Tanikawa M., Yamazaki M., Ninomiya K., Ishibashi T., Yamashita H., Murakawa K., Fujimori K., Tanai H., Kimata M., Watanabe M., Hiraoka S., Chiba Y., Ishida S., Ono Y., Takiguchi S., Watanabe S., Yosida M., Hotuta T., Kusano J., Kanehori K., Takahashi-Fujii A., Hara H., Tanase T.-O., Nomura Y., Togiya S., Komai F., Hara R., Takeuchi K., Arita M., Imose N., Musashino K., Yuuki H., Oshima A., Sasaki N., Aotsuka S., Yoshikawa Y., Matsunawa H., Ichihara T., Shiohata N., Sano S., Moriya S., Momiyama H., Satoh N., Takami S., Terashima Y., Suzuki O., Nakagawa S., Senoh A., Mizoguchi H., Goto Y., Shimizu F., Wakebe H., Hishigaki H., Watanabe T., Sugiyama A., Takemoto M., Kawakami B., Yamazaki M., Watanabe K., Kumagai A., Itakura S., Fukuzumi Y., Fujimori Y., Komiyama M., Tashiro H., Tanigami A., Fujiwara T., Ono T., Yamada K., Fujii Y., Ozaki K., Hirao M., Ohmori Y., Kawabata A., Hikiji T., Kobatake N., Inagaki H., Ikema Y., Okamoto S., Okitani R., Kawakami T., Noguchi S., Itoh T., Shigeta K., Senba T., Matsumura K., Nakajima Y., Mizuno T., Morinaga M., Sasaki M., Togashi T., Oyama M., Hata H., Watanabe M., Komatsu T., Mizushima-Sugano J., Satoh T., Shirai Y., Takahashi Y., Nakagawa K., Okumura K., Nagase T., Nomura N., Kikuchi H., Masuho Y., Yamashita R., Nakai K., Yada T., Nakamura Y., Ohara O., Isogai T., Sugano S.",
        )
        self.assertEqual(
            reference.title,
            "Complete sequencing and characterization of 21,243 full-length human cDNAs.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "14702039"))
        self.assertEqual(reference.references[1], ("DOI", "10.1038/ng1285"))
        reference = record.references[7]
        self.assertEqual(
            reference.authors,
            "Kalnine N., Chen X., Rolfs A., Halleck A., Hines L., Eisenstein S., Koundinya M., Raphael J., Moreira D., Kelley T., LaBaer J., Lin Y., Phelan M., Farmer A.",
        )
        self.assertEqual(
            reference.title,
            "Cloning of human full-length CDSs in BD Creator(TM) system donor vector.",
        )
        self.assertEqual(len(reference.references), 0)
        reference = record.references[8]
        self.assertEqual(
            reference.authors,
            "Mural R.J., Istrail S., Sutton G.G., Florea L., Halpern A.L., Mobarry C.M., Lippert R., Walenz B., Shatkay H., Dew I., Miller J.R., Flanigan M.J., Edwards N.J., Bolanos R., Fasulo D., Halldorsson B.V., Hannenhalli S., Turner R., Yooseph S., Lu F., Nusskern D.R., Shue B.C., Zheng X.H., Zhong F., Delcher A.L., Huson D.H., Kravitz S.A., Mouchard L., Reinert K., Remington K.A., Clark A.G., Waterman M.S., Eichler E.E., Adams M.D., Hunkapiller M.W., Myers E.W., Venter J.C.",
        )
        self.assertEqual(reference.title, "")
        self.assertEqual(len(reference.references), 0)
        reference = record.references[9]
        self.assertEqual(reference.authors, "The MGC Project Team")
        self.assertEqual(
            reference.title,
            "The status, quality, and expansion of the NIH full-length cDNA project: the Mammalian Gene Collection (MGC).",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "15489334"))
        self.assertEqual(reference.references[1], ("DOI", "10.1101/gr.2596504"))
        reference = record.references[10]
        self.assertEqual(
            reference.authors,
            "Gevaert K., Goethals M., Martens L., Van Damme J., Staes A., Thomas G.R., Vandekerckhove J.",
        )
        self.assertEqual(
            reference.title,
            "Exploring proteomes and analyzing protein processing by mass spectrometric identification of sorted N-terminal peptides.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "12665801"))
        self.assertEqual(reference.references[1], ("DOI", "10.1038/nbt810"))
        reference = record.references[11]
        self.assertEqual(
            reference.authors,
            "Greninger A.L., Knudsen G.M., Betegon M., Burlingame A.L., DeRisi J.L.",
        )
        self.assertEqual(
            reference.title,
            "ACBD3 interaction with TBC1 domain 22 protein is differentially affected by enteroviral and kobuviral 3A protein binding.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "23572552"))
        self.assertEqual(reference.references[1], ("DOI", "10.1128/mbio.00098-13"))
        reference = record.references[12]
        self.assertEqual(reference.authors, "Bienvenut W.V.")
        self.assertEqual(reference.title, "")
        self.assertEqual(len(reference.references), 0)
        reference = record.references[13]
        self.assertEqual(
            reference.authors,
            "Stewart S., Sundaram M., Zhang Y., Lee J., Han M., Guan K.L.",
        )
        self.assertEqual(
            reference.title,
            "Kinase suppressor of Ras forms a multiprotein signaling complex and modulates MEK localization.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "10409742"))
        self.assertEqual(reference.references[1], ("DOI", "10.1128/mcb.19.8.5523"))
        reference = record.references[14]
        self.assertEqual(
            reference.authors,
            "Demeter J., Medzihradszky D., Kha H., Goetzl E.J., Turck C.W.",
        )
        self.assertEqual(
            reference.title,
            "Isolation and partial characterization of the structures of fibroblast activating factor-related proteins from U937 cells.",
        )
        self.assertEqual(len(reference.references), 1)
        self.assertEqual(reference.references[0], ("PubMed", "2026444"))
        reference = record.references[15]
        self.assertEqual(reference.authors, "Lubec G., Afjehi-Sadat L.")
        self.assertEqual(reference.title, "")
        self.assertEqual(len(reference.references), 0)
        reference = record.references[16]
        self.assertEqual(
            reference.authors, "Aoki H., Hayashi J., Moriyama M., Arakawa Y., Hino O."
        )
        self.assertEqual(
            reference.title,
            "Hepatitis C virus core protein interacts with 14-3-3 protein and activates the kinase Raf-1.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "10644344"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1128/jvi.74.4.1736-1741.2000")
        )
        reference = record.references[17]
        self.assertEqual(
            reference.authors,
            "Ganguly S., Gastel J.A., Weller J.L., Schwartz C., Jaffe H., Namboodiri M.A., Coon S.L., Hickman A.B., Rollag M., Obsil T., Beauverger P., Ferry G., Boutin J.A., Klein D.C.",
        )
        self.assertEqual(
            reference.title,
            "Role of a pineal cAMP-operated arylalkylamine N-acetyltransferase/14-3-3-binding switch in melatonin synthesis.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "11427721"))
        self.assertEqual(reference.references[1], ("DOI", "10.1073/pnas.141118798"))
        reference = record.references[18]
        self.assertEqual(
            reference.authors, "Fujita N., Sato S., Katayama K., Tsuruo T."
        )
        self.assertEqual(
            reference.title,
            "Akt-dependent phosphorylation of p27Kip1 promotes binding to 14-3-3 and cytoplasmic localization.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "12042314"))
        self.assertEqual(reference.references[1], ("DOI", "10.1074/jbc.m203668200"))
        reference = record.references[19]
        self.assertEqual(
            reference.authors, "Wang X., Grammatikakis N., Siganou A., Calderwood S.K."
        )
        self.assertEqual(
            reference.title,
            "Regulation of molecular chaperone gene transcription involves the serine phosphorylation, 14-3-3 epsilon binding, and cytoplasmic sequestration of heat shock factor 1.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "12917326"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1128/mcb.23.17.6013-6026.2003")
        )
        reference = record.references[20]
        self.assertEqual(
            reference.authors,
            "Andersen J.S., Wilkinson C.J., Mayor T., Mortensen P., Nigg E.A., Mann M.",
        )
        self.assertEqual(
            reference.title,
            "Proteomic characterization of the human centrosome by protein correlation profiling.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "14654843"))
        self.assertEqual(reference.references[1], ("DOI", "10.1038/nature02166"))
        reference = record.references[21]
        self.assertEqual(
            reference.authors,
            "Urschel S., Bassermann F., Bai R.Y., Munch S., Peschel C., Duyster J.",
        )
        self.assertEqual(
            reference.title,
            "Phosphorylation of grb10 regulates its interaction with 14-3-3.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "15722337"))
        self.assertEqual(reference.references[1], ("DOI", "10.1074/jbc.m501477200"))
        reference = record.references[22]
        self.assertEqual(
            reference.authors, "Yoshida K., Yamaguchi T., Natsume T., Kufe D., Miki Y."
        )
        self.assertEqual(
            reference.title,
            "JNK phosphorylation of 14-3-3 proteins regulates nuclear targeting of c-Abl in the apoptotic response to DNA damage.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "15696159"))
        self.assertEqual(reference.references[1], ("DOI", "10.1038/ncb1228"))
        reference = record.references[23]
        self.assertEqual(
            reference.authors,
            "Gu Y.-M., Jin Y.-H., Choi J.-K., Baek K.-H., Yeo C.-Y., Lee K.-Y.",
        )
        self.assertEqual(
            reference.title,
            "Protein kinase A phosphorylates and regulates dimerization of 14-3-3 epsilon.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "16376338"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1016/j.febslet.2005.12.024")
        )
        reference = record.references[24]
        self.assertEqual(
            reference.authors,
            "Chi A., Valencia J.C., Hu Z.-Z., Watabe H., Yamaguchi H., Mangini N.J., Huang H., Canfield V.A., Cheng K.C., Yang F., Abe R., Yamagishi S., Shabanowitz J., Hearing V.J., Wu C., Appella E., Hunt D.F.",
        )
        self.assertEqual(
            reference.title,
            "Proteomic and bioinformatic characterization of the biogenesis and function of melanosomes.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "17081065"))
        self.assertEqual(reference.references[1], ("DOI", "10.1021/pr060363j"))
        reference = record.references[25]
        self.assertEqual(
            reference.authors,
            "Linde C.I., Di Leva F., Domi T., Tosatto S.C., Brini M., Carafoli E.",
        )
        self.assertEqual(
            reference.title,
            "Inhibitory interaction of the 14-3-3 proteins with ubiquitous (PMCA1) and tissue-specific (PMCA3) isoforms of the plasma membrane Ca2+ pump.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "18029012"))
        self.assertEqual(reference.references[1], ("DOI", "10.1016/j.ceca.2007.09.003"))
        reference = record.references[26]
        self.assertEqual(
            reference.authors,
            "Brummer T., Larance M., Herrera Abreu M.T., Lyons R.J., Timpson P., Emmerich C.H., Fleuren E.D.G., Lehrbach G.M., Schramek D., Guilhaus M., James D.E., Daly R.J.",
        )
        self.assertEqual(
            reference.title,
            "Phosphorylation-dependent binding of 14-3-3 terminates signalling by the Gab2 docking protein.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "19172738"))
        self.assertEqual(reference.references[1], ("DOI", "10.1038/emboj.2008.159"))
        reference = record.references[27]
        self.assertEqual(
            reference.authors,
            "Han G., Ye M., Zhou H., Jiang X., Feng S., Jiang X., Tian R., Wan D., Zou H., Gu J.",
        )
        self.assertEqual(
            reference.title,
            "Large-scale phosphoproteome analysis of human liver tissue by enrichment and fractionation of phosphopeptides with strong anion exchange chromatography.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "18318008"))
        self.assertEqual(reference.references[1], ("DOI", "10.1002/pmic.200700884"))
        reference = record.references[28]
        self.assertEqual(
            reference.authors,
            "Gauci S., Helbig A.O., Slijper M., Krijgsveld J., Heck A.J., Mohammed S.",
        )
        self.assertEqual(
            reference.title,
            "Lys-N and trypsin cover complementary parts of the phosphoproteome in a refined SCX-based approach.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "19413330"))
        self.assertEqual(reference.references[1], ("DOI", "10.1021/ac9004309"))
        reference = record.references[29]
        self.assertEqual(reference.authors, "Kajiwara Y., Buxbaum J.D., Grice D.E.")
        self.assertEqual(
            reference.title,
            "SLITRK1 binds 14-3-3 and regulates neurite outgrowth in a phosphorylation-dependent manner.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "19640509"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1016/j.biopsych.2009.05.033")
        )
        reference = record.references[30]
        self.assertEqual(
            reference.authors,
            "Jang S.W., Liu X., Fu H., Rees H., Yepes M., Levey A., Ye K.",
        )
        self.assertEqual(
            reference.title,
            "Interaction of Akt-phosphorylated SRPK2 with 14-3-3 mediates cell cycle and cell death in neurons.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "19592491"))
        self.assertEqual(reference.references[1], ("DOI", "10.1074/jbc.m109.026237"))
        reference = record.references[31]
        self.assertEqual(
            reference.authors,
            "Choudhary C., Kumar C., Gnad F., Nielsen M.L., Rehman M., Walther T.C., Olsen J.V., Mann M.",
        )
        self.assertEqual(
            reference.title,
            "Lysine acetylation targets protein complexes and co-regulates major cellular functions.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "19608861"))
        self.assertEqual(reference.references[1], ("DOI", "10.1126/science.1175371"))
        reference = record.references[32]
        self.assertEqual(
            reference.authors,
            "Olsen J.V., Vermeulen M., Santamaria A., Kumar C., Miller M.L., Jensen L.J., Gnad F., Cox J., Jensen T.S., Nigg E.A., Brunak S., Mann M.",
        )
        self.assertEqual(
            reference.title,
            "Quantitative phosphoproteomics reveals widespread full phosphorylation site occupancy during mitosis.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "20068231"))
        self.assertEqual(reference.references[1], ("DOI", "10.1126/scisignal.2000475"))
        reference = record.references[33]
        self.assertEqual(
            reference.authors,
            "Burkard T.R., Planyavsky M., Kaupe I., Breitwieser F.P., Buerckstuemmer T., Bennett K.L., Superti-Furga G., Colinge J.",
        )
        self.assertEqual(
            reference.title, "Initial characterization of the human central proteome."
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "21269460"))
        self.assertEqual(reference.references[1], ("DOI", "10.1186/1752-0509-5-17"))
        reference = record.references[34]
        self.assertEqual(
            reference.authors,
            "Rigbolt K.T., Prokhorova T.A., Akimov V., Henningsen J., Johansen P.T., Kratchmarova I., Kassem M., Mann M., Olsen J.V., Blagoev B.",
        )
        self.assertEqual(
            reference.title,
            "System-wide temporal characterization of the proteome and phosphoproteome of human embryonic stem cell differentiation.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "21406692"))
        self.assertEqual(reference.references[1], ("DOI", "10.1126/scisignal.2001570"))
        reference = record.references[35]
        self.assertEqual(
            reference.authors,
            "Zhou H., Di Palma S., Preisinger C., Peng M., Polat A.N., Heck A.J., Mohammed S.",
        )
        self.assertEqual(
            reference.title,
            "Toward a comprehensive characterization of a human cancer cell phosphoproteome.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "23186163"))
        self.assertEqual(reference.references[1], ("DOI", "10.1021/pr300630k"))
        reference = record.references[36]
        self.assertEqual(
            reference.authors,
            "Bian Y., Song C., Cheng K., Dong M., Wang F., Huang J., Sun D., Wang L., Ye M., Zou H.",
        )
        self.assertEqual(
            reference.title,
            "An enzyme assisted RP-RPLC approach for in-depth analysis of human liver phosphoproteome.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "24275569"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1016/j.jprot.2013.11.014")
        )
        reference = record.references[37]
        self.assertEqual(
            reference.authors,
            "Yuasa K., Ota R., Matsuda S., Isshiki K., Inoue M., Tsuji A.",
        )
        self.assertEqual(
            reference.title,
            "Suppression of death-associated protein kinase 2 by interaction with 14-3-3 proteins.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "26047703"))
        self.assertEqual(reference.references[1], ("DOI", "10.1016/j.bbrc.2015.05.105"))
        reference = record.references[38]
        self.assertEqual(
            reference.authors,
            "Kulasekaran G., Nossova N., Marat A.L., Lund I., Cremer C., Ioannou M.S., McPherson P.S.",
        )
        self.assertEqual(
            reference.title,
            "Phosphorylation-dependent regulation of Connecdenn/DENND1 guanine nucleotide exchange factors.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "26055712"))
        self.assertEqual(reference.references[1], ("DOI", "10.1074/jbc.m115.636712"))
        reference = record.references[39]
        self.assertEqual(
            reference.authors,
            "Gao K., Tang W., Li Y., Zhang P., Wang D., Yu L., Wang C., Wu D.",
        )
        self.assertEqual(
            reference.title,
            "Front-signal-dependent accumulation of the RHOA inhibitor FAM65B at leading edges polarizes neutrophils.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "25588844"))
        self.assertEqual(reference.references[1], ("DOI", "10.1242/jcs.161497"))
        reference = record.references[40]
        self.assertEqual(
            reference.authors,
            "Vaca Jacome A.S., Rabilloud T., Schaeffer-Reiss C., Rompais M., Ayoub D., Lane L., Bairoch A., Van Dorsselaer A., Carapito C.",
        )
        self.assertEqual(
            reference.title, "N-terminome analysis of the human mitochondrial proteome."
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "25944712"))
        self.assertEqual(reference.references[1], ("DOI", "10.1002/pmic.201400617"))
        reference = record.references[41]
        self.assertEqual(
            reference.authors,
            "Masters S.L., Lagou V., Jeru I., Baker P.J., Van Eyck L., Parry D.A., Lawless D., De Nardo D., Garcia-Perez J.E., Dagley L.F., Holley C.L., Dooley J., Moghaddas F., Pasciuto E., Jeandel P.Y., Sciot R., Lyras D., Webb A.I., Nicholson S.E., De Somer L., van Nieuwenhove E., Ruuth-Praz J., Copin B., Cochet E., Medlej-Hashim M., Megarbane A., Schroder K., Savic S., Goris A., Amselem S., Wouters C., Liston A.",
        )
        self.assertEqual(
            reference.title,
            "Familial autoinflammation with neutrophilic dermatosis reveals a regulatory mechanism of pyrin activation.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "27030597"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1126/scitranslmed.aaf1471")
        )
        reference = record.references[42]
        self.assertEqual(
            reference.authors,
            "Hendriks I.A., Lyon D., Young C., Jensen L.J., Vertegaal A.C., Nielsen M.L.",
        )
        self.assertEqual(
            reference.title,
            "Site-specific mapping of the human SUMO proteome reveals co-modification with phosphorylation.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "28112733"))
        self.assertEqual(reference.references[1], ("DOI", "10.1038/nsmb.3366"))
        reference = record.references[43]
        self.assertEqual(
            reference.authors,
            "Sonntag T., Ostojic J., Vaughan J.M., Moresco J.J., Yoon Y.S., Yates J.R. III, Montminy M.",
        )
        self.assertEqual(
            reference.title,
            "Mitogenic Signals Stimulate the CREB Coactivator CRTC3 through PP2A Recruitment.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "30611118"))
        self.assertEqual(reference.references[1], ("DOI", "10.1016/j.isci.2018.12.012"))
        reference = record.references[44]
        self.assertEqual(
            reference.authors, "Chen J., Ou Y., Yang Y., Li W., Xu Y., Xie Y., Liu Y."
        )
        self.assertEqual(
            reference.title,
            "KLHL22 activates amino-acid-dependent mTORC1 signalling to promote tumorigenesis and ageing.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "29769719"))
        self.assertEqual(reference.references[1], ("DOI", "10.1038/s41586-018-0128-9"))
        reference = record.references[45]
        self.assertEqual(
            reference.authors,
            "Yang X., Lee W.H., Sobott F., Papagrigoriou E., Robinson C.V., Grossmann J.G., Sundstroem M., Doyle D.A., Elkins J.M.",
        )
        self.assertEqual(
            reference.title,
            "Structural basis for protein-protein interactions in the 14-3-3 protein family.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "17085597"))
        self.assertEqual(reference.references[1], ("DOI", "10.1073/pnas.0605779103"))

        # Check the two parsers agree on the essentials
        self.assertEqual(seq_record.seq, record.sequence)
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
        self.assertEqual(records[0].seq, seq_record.seq)
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
        self.assertEqual(seq_record.seq, record.sequence)
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
        self.assertEqual(records[0].seq, seq_record.seq)
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

    def test_P68308(self):
        """Parsing SwissProt file P68308.txt."""
        filename = "P68308.txt"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)

        with open(datafile) as test_handle:
            seq_record = SeqIO.read(test_handle, "swiss")

        self.assertIsInstance(seq_record, SeqRecord)

        self.assertEqual(seq_record.id, "P68308")
        self.assertEqual(seq_record.name, "NU3M_BALPH")
        self.assertEqual(
            seq_record.description,
            "RecName: Full=NADH-ubiquinone oxidoreductase chain 3 {ECO:0000250|UniProtKB:P03897}; EC=7.1.1.2 {ECO:0000250|UniProtKB:P03897}; AltName: Full=NADH dehydrogenase subunit 3;",
        )
        self.assertEqual(
            repr(seq_record.seq),
            "Seq('MNLLLTLLTNTTLALLLVFIAFWLPQLNVYAEKTSPYECGFDPMGSARLPFSMK...WAE')",
        )

        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(record.entry_name, "NU3M_BALPH")
        self.assertEqual(record.accessions, ["P68308", "P24973"])
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
                "Laurasiatheria",
                "Artiodactyla",
                "Whippomorpha",
                "Cetacea",
                "Mysticeti",
                "Balaenopteridae",
                "Balaenoptera",
            ],
        )
        self.assertEqual(record.seqinfo, (115, 13022, "405197D2F5D0AC4B"))

        self.assertEqual(len(record.features), 4)
        feature = record.features[0]
        self.assertEqual(feature.type, "CHAIN")
        self.assertEqual(feature.location.start, 0)
        self.assertEqual(feature.location.end, 115)
        self.assertEqual(
            feature.qualifiers["note"], "NADH-ubiquinone oxidoreductase chain 3"
        )
        feature = record.features[1]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 2)
        self.assertEqual(feature.location.end, 23)
        self.assertEqual(feature.qualifiers["note"], "Helical")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[2]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 54)
        self.assertEqual(feature.location.end, 75)
        self.assertEqual(feature.qualifiers["note"], "Helical")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[3]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 83)
        self.assertEqual(feature.location.end, 104)
        self.assertEqual(feature.qualifiers["note"], "Helical")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        self.assertIsNone(feature.id)

        self.assertEqual(len(record.references), 1)
        reference = record.references[0]
        self.assertEqual(reference.authors, "Arnason U., Gullberg A., Widegren B.")
        self.assertEqual(
            reference.title,
            "The complete nucleotide sequence of the mitochondrial DNA of the fin whale, Balaenoptera physalus.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "1779436"))
        self.assertEqual(reference.references[1], ("DOI", "10.1007/bf02102808"))

        # Check the two parsers agree on the essentials
        self.assertEqual(seq_record.seq, record.sequence)
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
        self.assertEqual(records[0].seq, seq_record.seq)
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

    def test_P39896(self):
        """Parsing SwissProt file P39896.txt."""
        filename = "P39896.txt"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)

        with open(datafile) as test_handle:
            seq_record = SeqIO.read(test_handle, "swiss")

        self.assertIsInstance(seq_record, SeqRecord)

        self.assertEqual(seq_record.id, "P39896")
        self.assertEqual(seq_record.name, "TCMO_STRGA")
        self.assertEqual(
            seq_record.description,
            "RecName: Full=Tetracenomycin polyketide synthesis 8-O-methyl transferase TcmO; EC=2.1.1.-;",
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
                "Bacteria",
                "Actinobacteria",
                "Streptomycetales",
                "Streptomycetaceae",
                "Streptomyces",
            ],
        )
        self.assertEqual(record.seqinfo, (339, 37035, "B228B66B24217F80"))

        self.assertEqual(len(record.features), 4)
        feature = record.features[0]
        self.assertEqual(feature.type, "CHAIN")
        self.assertEqual(feature.location.start, 0)
        self.assertEqual(feature.location.end, 339)
        self.assertEqual(
            feature.qualifiers["note"],
            "Tetracenomycin polyketide synthesis 8-O-methyl transferase TcmO",
        )
        feature = record.features[1]
        self.assertEqual(feature.type, "ACT_SITE")
        self.assertEqual(feature.location.start, 245)
        self.assertEqual(feature.location.end, 246)
        self.assertEqual(feature.qualifiers["note"], "Proton acceptor")
        self.assertEqual(
            feature.qualifiers["evidence"], "ECO:0000255|PROSITE-ProRule:PRU01020"
        )
        feature = record.features[2]
        self.assertEqual(feature.type, "BINDING")
        self.assertEqual(feature.location.start, 199)
        self.assertEqual(feature.location.end, 200)
        self.assertEqual(feature.qualifiers["ligand"], "S-adenosyl-L-methionine")
        self.assertEqual(feature.qualifiers["ligand_id"], "ChEBI:CHEBI:59789")
        self.assertEqual(
            feature.qualifiers["evidence"], "ECO:0000255|PROSITE-ProRule:PRU01020"
        )
        feature = record.features[3]
        self.assertEqual(feature.type, "BINDING")
        self.assertEqual(feature.location.start, 225)
        self.assertEqual(feature.location.end, 228)
        self.assertEqual(feature.qualifiers["ligand"], "S-adenosyl-L-methionine")
        self.assertEqual(feature.qualifiers["ligand_id"], "ChEBI:CHEBI:59789")
        self.assertEqual(
            feature.qualifiers["evidence"], "ECO:0000255|PROSITE-ProRule:PRU01020"
        )
        self.assertIsNone(feature.id)

        self.assertEqual(len(record.references), 1)
        reference = record.references[0]
        self.assertEqual(
            reference.authors,
            "Summers R.G., Wendt-Pienkowski E., Motamedi H., Hutchinson C.R.",
        )
        self.assertEqual(
            reference.title,
            "Nucleotide sequence of the tcmII-tcmIV region of the tetracenomycin C biosynthetic gene cluster of Streptomyces glaucescens and evidence that the tcmN gene encodes a multifunctional cyclase-dehydratase-O-methyl transferase.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "1548230"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1128/jb.174.6.1810-1820.1992")
        )

        # Check the two parsers agree on the essentials
        self.assertEqual(seq_record.seq, record.sequence)
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
        self.assertEqual(records[0].seq, seq_record.seq)
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

    def test_O95832(self):
        """Parsing SwissProt file O95832.txt."""
        filename = "O95832.txt"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)

        with open(datafile) as test_handle:
            seq_record = SeqIO.read(test_handle, "swiss")

        self.assertIsInstance(seq_record, SeqRecord)

        self.assertEqual(seq_record.id, "O95832")
        self.assertEqual(seq_record.name, "CLD1_HUMAN")
        self.assertEqual(
            seq_record.description,
            "RecName: Full=Claudin-1; AltName: Full=Senescence-associated epithelial membrane protein;",
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
        self.assertEqual(record.seqinfo, (211, 22744, "07269000E6C214F0"))

        self.assertEqual(len(record.features), 17)
        feature = record.features[0]
        self.assertEqual(feature.type, "CHAIN")
        self.assertEqual(feature.location.start, 0)
        self.assertEqual(feature.location.end, 211)
        self.assertEqual(feature.qualifiers["note"], "Claudin-1")
        feature = record.features[1]
        self.assertEqual(feature.type, "TOPO_DOM")
        self.assertEqual(feature.location.start, 0)
        self.assertEqual(feature.location.end, 7)
        self.assertEqual(feature.qualifiers["note"], "Cytoplasmic")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000250")
        feature = record.features[2]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 7)
        self.assertEqual(feature.location.end, 28)
        self.assertEqual(feature.qualifiers["note"], "Helical")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[3]
        self.assertEqual(feature.type, "TOPO_DOM")
        self.assertEqual(feature.location.start, 28)
        self.assertEqual(feature.location.end, 81)
        self.assertEqual(feature.qualifiers["note"], "Extracellular")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[4]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 81)
        self.assertEqual(feature.location.end, 102)
        self.assertEqual(feature.qualifiers["note"], "Helical")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[5]
        self.assertEqual(feature.type, "TOPO_DOM")
        self.assertEqual(feature.location.start, 102)
        self.assertEqual(feature.location.end, 115)
        self.assertEqual(feature.qualifiers["note"], "Cytoplasmic")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[6]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 115)
        self.assertEqual(feature.location.end, 136)
        self.assertEqual(feature.qualifiers["note"], "Helical")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[7]
        self.assertEqual(feature.type, "TOPO_DOM")
        self.assertEqual(feature.location.start, 136)
        self.assertEqual(feature.location.end, 163)
        self.assertEqual(feature.qualifiers["note"], "Extracellular")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[8]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 163)
        self.assertEqual(feature.location.end, 184)
        self.assertEqual(feature.qualifiers["note"], "Helical")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[9]
        self.assertEqual(feature.type, "TOPO_DOM")
        self.assertEqual(feature.location.start, 184)
        self.assertEqual(feature.location.end, 211)
        self.assertEqual(feature.qualifiers["note"], "Cytoplasmic")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000250")
        feature = record.features[10]
        self.assertEqual(feature.type, "REGION")
        self.assertEqual(feature.location.start, 191)
        self.assertEqual(feature.location.end, 211)
        self.assertEqual(feature.qualifiers["note"], "Disordered")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000256|SAM:MobiDB-lite")
        feature = record.features[11]
        self.assertEqual(feature.type, "REGION")
        self.assertEqual(feature.location.start, 209)
        self.assertEqual(feature.location.end, 211)
        self.assertEqual(
            feature.qualifiers["note"], "Interactions with TJP1, TJP2, TJP3 and PATJ"
        )
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000250")
        feature = record.features[12]
        self.assertEqual(feature.type, "DISULFID")
        self.assertEqual(feature.location.start, 53)
        self.assertEqual(feature.location.end, 64)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000250")
        feature = record.features[13]
        self.assertEqual(feature.type, "MUTAGEN")
        self.assertEqual(feature.location.start, 31)
        self.assertEqual(feature.location.end, 32)
        self.assertEqual(
            feature.qualifiers["note"],
            "I->M: Loss of HCV receptor activity. Significant loss of interaction with CD81. Reduced interaction with OCLN.",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:17325668,ECO:0000269|PubMed:20375010",
        )
        feature = record.features[14]
        self.assertEqual(feature.type, "MUTAGEN")
        self.assertEqual(feature.location.start, 47)
        self.assertEqual(feature.location.end, 48)
        self.assertEqual(
            feature.qualifiers["note"],
            "E->K: Loss of HCV receptor activity. Significant loss of interaction with CD81. Reduced interaction with OCLN. According to PubMed:17325668 no effect observed on HCV infection susceptibility in cell culture.",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:17325668,ECO:0000269|PubMed:20375010",
        )
        feature = record.features[15]
        self.assertEqual(feature.type, "CONFLICT")
        self.assertEqual(feature.location.start, 61)
        self.assertEqual(feature.location.end, 62)
        self.assertEqual(feature.qualifiers["note"], "I -> V (in Ref. 2; AAD22962)")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000305")
        feature = record.features[16]
        self.assertEqual(feature.type, "CONFLICT")
        self.assertEqual(feature.location.start, 134)
        self.assertEqual(feature.location.end, 135)
        self.assertEqual(feature.qualifiers["note"], "V -> A (in Ref. 2; AAD22962)")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000305")
        self.assertIsNone(feature.id)

        self.assertEqual(len(record.references), 17)
        reference = record.references[0]
        self.assertEqual(
            reference.authors,
            "Swisshelm K.L., Machl A., Planitzer S., Robertson R., Kubbies M., Hosier S.",
        )
        self.assertEqual(
            reference.title,
            "SEMP1, a senescence-associated cDNA isolated from human mammary epithelial cells, is a member of an epithelial membrane protein superfamily.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "9931503"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1016/s0378-1119(98)00553-8")
        )
        reference = record.references[1]
        self.assertEqual(reference.authors, "Mitic L.M., Anderson J.M.")
        self.assertEqual(reference.title, "Human claudin-1 isolated from Caco-2 mRNA.")
        self.assertEqual(len(reference.references), 0)
        reference = record.references[2]
        self.assertEqual(
            reference.authors,
            "Halford S., Spencer P., Greenwood J., Winton H., Hunt D.M., Adamson P.",
        )
        self.assertEqual(
            reference.title,
            "Assignment of claudin-1 (CLDN1) to human chromosome 3q28-->q29 with somatic cell hybrids.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "10828592"))
        self.assertEqual(reference.references[1], ("DOI", "10.1159/000015553"))
        reference = record.references[3]
        self.assertEqual(
            reference.authors,
            "Kraemer F., White K., Kubbies M., Swisshelm K.L., Weber B.H.F.",
        )
        self.assertEqual(
            reference.title,
            "Genomic organization of claudin-1 and its assessment in hereditary and sporadic breast cancer.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "11071387"))
        self.assertEqual(reference.references[1], ("DOI", "10.1007/s004390000375"))
        reference = record.references[4]
        self.assertEqual(
            reference.authors,
            "Clark H.F., Gurney A.L., Abaya E., Baker K., Baldwin D.T., Brush J., Chen J., Chow B., Chui C., Crowley C., Currell B., Deuel B., Dowd P., Eaton D., Foster J.S., Grimaldi C., Gu Q., Hass P.E., Heldens S., Huang A., Kim H.S., Klimowski L., Jin Y., Johnson S., Lee J., Lewis L., Liao D., Mark M.R., Robbie E., Sanchez C., Schoenfeld J., Seshagiri S., Simmons L., Singh J., Smith V., Stinson J., Vagts A., Vandlen R.L., Watanabe C., Wieand D., Woods K., Xie M.-H., Yansura D.G., Yi S., Yu G., Yuan J., Zhang M., Zhang Z., Goddard A.D., Wood W.I., Godowski P.J., Gray A.M.",
        )
        self.assertEqual(
            reference.title,
            "The secreted protein discovery initiative (SPDI), a large-scale effort to identify novel human secreted and transmembrane proteins: a bioinformatics assessment.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "12975309"))
        self.assertEqual(reference.references[1], ("DOI", "10.1101/gr.1293003"))
        reference = record.references[5]
        self.assertEqual(
            reference.authors,
            "Ebert L., Schick M., Neubert P., Schatten R., Henze S., Korn B.",
        )
        self.assertEqual(
            reference.title,
            "Cloning of human full open reading frames in Gateway(TM) system entry vector (pDONR201).",
        )
        self.assertEqual(len(reference.references), 0)
        reference = record.references[6]
        self.assertEqual(reference.authors, "The MGC Project Team")
        self.assertEqual(
            reference.title,
            "The status, quality, and expansion of the NIH full-length cDNA project: the Mammalian Gene Collection (MGC).",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "15489334"))
        self.assertEqual(reference.references[1], ("DOI", "10.1101/gr.2596504"))
        reference = record.references[7]
        self.assertEqual(
            reference.authors,
            "Hadj-Rabia S., Baala L., Vabres P., Hamel-Teillac D., Jacquemin E., Fabre M., Lyonnet S., De Prost Y., Munnich A., Hadchouel M., Smahi A.",
        )
        self.assertEqual(
            reference.title,
            "Claudin-1 gene mutations in neonatal sclerosing cholangitis associated with ichthyosis: a tight junction disease.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "15521008"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1053/j.gastro.2004.07.022")
        )
        reference = record.references[8]
        self.assertEqual(
            reference.authors,
            "Feldmeyer L., Huber M., Fellmann F., Beckmann J.S., Frenk E., Hohl D.",
        )
        self.assertEqual(
            reference.title, "Confirmation of the origin of NISCH syndrome."
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "16619213"))
        self.assertEqual(reference.references[1], ("DOI", "10.1002/humu.20333"))
        reference = record.references[9]
        self.assertEqual(
            reference.authors,
            "Evans M.J., von Hahn T., Tscherne D.M., Syder A.J., Panis M., Wolk B., Hatziioannou T., McKeating J.A., Bieniasz P.D., Rice C.M.",
        )
        self.assertEqual(
            reference.title,
            "Claudin-1 is a hepatitis C virus co-receptor required for a late step in entry.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "17325668"))
        self.assertEqual(reference.references[1], ("DOI", "10.1038/nature05654"))
        reference = record.references[10]
        self.assertEqual(
            reference.authors,
            "Harris H.J., Davis C., Mullins J.G., Hu K., Goodall M., Farquhar M.J., Mee C.J., McCaffrey K., Young S., Drummer H., Balfe P., McKeating J.A.",
        )
        self.assertEqual(
            reference.title,
            "Claudin association with CD81 defines hepatitis C virus entry.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "20375010"))
        self.assertEqual(reference.references[1], ("DOI", "10.1074/jbc.m110.104836"))
        reference = record.references[11]
        self.assertEqual(
            reference.authors,
            "Burkard T.R., Planyavsky M., Kaupe I., Breitwieser F.P., Buerckstuemmer T., Bennett K.L., Superti-Furga G., Colinge J.",
        )
        self.assertEqual(
            reference.title, "Initial characterization of the human central proteome."
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "21269460"))
        self.assertEqual(reference.references[1], ("DOI", "10.1186/1752-0509-5-17"))
        reference = record.references[12]
        self.assertEqual(
            reference.authors,
            "Lupberger J., Zeisel M.B., Xiao F., Thumann C., Fofana I., Zona L., Davis C., Mee C.J., Turek M., Gorke S., Royer C., Fischer B., Zahid M.N., Lavillette D., Fresquet J., Cosset F.L., Rothenberg S.M., Pietschmann T., Patel A.H., Pessaux P., Doffoel M., Raffelsberger W., Poch O., McKeating J.A., Brino L., Baumert T.F.",
        )
        self.assertEqual(
            reference.title,
            "EGFR and EphA2 are host factors for hepatitis C virus entry and possible targets for antiviral therapy.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "21516087"))
        self.assertEqual(reference.references[1], ("DOI", "10.1038/nm.2341"))
        reference = record.references[13]
        self.assertEqual(
            reference.authors,
            "Kirschner N., Rosenthal R., Furuse M., Moll I., Fromm M., Brandner J.M.",
        )
        self.assertEqual(
            reference.title,
            "Contribution of tight junction proteins to ion, macromolecule, and water barrier in keratinocytes.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "23407391"))
        self.assertEqual(reference.references[1], ("DOI", "10.1038/jid.2012.507"))
        reference = record.references[14]
        self.assertEqual(reference.authors, "Che P., Tang H., Li Q.")
        self.assertEqual(
            reference.title,
            "The interaction between claudin-1 and dengue viral prM/M protein for its entry.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "24074594"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1016/j.virol.2013.08.009")
        )
        reference = record.references[15]
        self.assertEqual(
            reference.authors,
            "Bonander N., Jamshad M., Oberthuer D., Clare M., Barwell J., Hu K., Farquhar M.J., Stamataki Z., Harris H.J., Dierks K., Dafforn T.R., Betzel C., McKeating J.A., Bill R.M.",
        )
        self.assertEqual(
            reference.title,
            "Production, purification and characterization of recombinant, full-length human claudin-1.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "23704991"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1371/journal.pone.0064517")
        )
        reference = record.references[16]
        self.assertEqual(
            reference.authors,
            "Douam F., Dao Thi V.L., Maurin G., Fresquet J., Mompelat D., Zeisel M.B., Baumert T.F., Cosset F.L., Lavillette D.",
        )
        self.assertEqual(
            reference.title,
            "Critical interaction between E1 and E2 glycoproteins determines binding and fusion properties of hepatitis C virus during cell entry.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "24038151"))
        self.assertEqual(reference.references[1], ("DOI", "10.1002/hep.26733"))

        # Check the two parsers agree on the essentials
        self.assertEqual(seq_record.seq, record.sequence)
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
        self.assertEqual(records[0].seq, seq_record.seq)
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

    def test_P04439(self):
        """Parsing SwissProt file P04439.txt."""
        filename = "P04439.txt"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)

        with open(datafile) as test_handle:
            seq_record = SeqIO.read(test_handle, "swiss")

        self.assertIsInstance(seq_record, SeqRecord)

        self.assertEqual(seq_record.id, "P04439")
        self.assertEqual(seq_record.name, "HLAA_HUMAN")
        self.assertEqual(
            seq_record.description,
            "RecName: Full=HLA class I histocompatibility antigen, A alpha chain; AltName: Full=Human leukocyte antigen A; Short=HLA-A; Flags: Precursor;",
        )
        self.assertEqual(
            repr(seq_record.seq),
            "Seq('MAVMAPRTLLLLLSGALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGYVDD...CKV')",
        )

        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(record.entry_name, "HLAA_HUMAN")
        self.assertEqual(
            record.accessions,
            [
                "P04439",
                "B1PKZ3",
                "O02939",
                "O02954",
                "O02955",
                "O02963",
                "O19509",
                "O19546",
                "O19598",
                "O19605",
                "O19606",
                "O19619",
                "O19647",
                "O19673",
                "O19687",
                "O19695",
                "O19756",
                "O19794",
                "O19795",
                "O43906",
                "O43907",
                "O46874",
                "O62921",
                "O62924",
                "O77937",
                "O77938",
                "O77964",
                "O78073",
                "O78171",
                "O98009",
                "O98010",
                "O98011",
                "O98137",
                "P01891",
                "P01892",
                "P05534",
                "P06338",
                "P10313",
                "P10314",
                "P10315",
                "P10316",
                "P13746",
                "P16188",
                "P16189",
                "P16190",
                "P18462",
                "P30443",
                "P30444",
                "P30445",
                "P30446",
                "P30447",
                "P30448",
                "P30449",
                "P30450",
                "P30451",
                "P30452",
                "P30453",
                "P30454",
                "P30455",
                "P30456",
                "P30457",
                "P30458",
                "P30459",
                "P30512",
                "P30514",
                "P79505",
                "P79562",
                "P79563",
                "Q09160",
                "Q29680",
                "Q29747",
                "Q29835",
                "Q29837",
                "Q29838",
                "Q29899",
                "Q29908",
                "Q29909",
                "Q29910",
                "Q30208",
                "Q31623",
                "Q5S3G1",
                "Q65A82",
                "Q8MHM1",
                "Q8MHN9",
                "Q95352",
                "Q95355",
                "Q95362",
                "Q95377",
                "Q95380",
                "Q95IZ5",
                "Q9BCN0",
                "Q9BD15",
                "Q9BD19",
                "Q9GJE6",
                "Q9GJE7",
                "Q9GJE8",
                "Q9MW42",
                "Q9MY89",
                "Q9MYA3",
                "Q9MYA5",
                "Q9MYC4",
                "Q9MYE6",
                "Q9MYE9",
                "Q9MYG4",
                "Q9MYG5",
                "Q9MYI5",
                "Q9TP25",
                "Q9TPQ3",
                "Q9TPR8",
                "Q9TPX8",
                "Q9TPX9",
                "Q9TPY0",
                "Q9TQ24",
                "Q9TQE8",
                "Q9TQE9",
                "Q9TQF1",
                "Q9TQF5",
                "Q9TQF8",
                "Q9TQF9",
                "Q9TQG0",
                "Q9TQG5",
                "Q9TQG7",
                "Q9TQH5",
                "Q9TQI3",
                "Q9TQK5",
                "Q9TQM6",
                "Q9TQN5",
                "Q9TQP5",
                "Q9TQP6",
                "Q9TQP7",
                "Q9UIN1",
                "Q9UIN2",
                "Q9UIP7",
                "Q9UQU3",
                "Q9UQU6",
                "Q9UQU7",
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
        self.assertEqual(record.seqinfo, (365, 40841, "DEDFCEC4450E0580"))

        self.assertEqual(len(record.features), 161)
        feature = record.features[0]
        self.assertEqual(feature.type, "SIGNAL")
        self.assertEqual(feature.location.start, 0)
        self.assertEqual(feature.location.end, 24)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000269|PubMed:92029")
        feature = record.features[1]
        self.assertEqual(feature.type, "CHAIN")
        self.assertEqual(feature.location.start, 24)
        self.assertEqual(feature.location.end, 365)
        self.assertEqual(
            feature.qualifiers["note"],
            "HLA class I histocompatibility antigen, A alpha chain",
        )
        feature = record.features[2]
        self.assertEqual(feature.type, "TOPO_DOM")
        self.assertEqual(feature.location.start, 24)
        self.assertEqual(feature.location.end, 308)
        self.assertEqual(feature.qualifiers["note"], "Extracellular")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[3]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 308)
        self.assertEqual(feature.location.end, 332)
        self.assertEqual(feature.qualifiers["note"], "Helical")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[4]
        self.assertEqual(feature.type, "TOPO_DOM")
        self.assertEqual(feature.location.start, 332)
        self.assertEqual(feature.location.end, 365)
        self.assertEqual(feature.qualifiers["note"], "Cytoplasmic")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[5]
        self.assertEqual(feature.type, "DOMAIN")
        self.assertEqual(feature.location.start, 208)
        self.assertEqual(feature.location.end, 295)
        self.assertEqual(feature.qualifiers["note"], "Ig-like C1-type")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[6]
        self.assertEqual(feature.type, "REGION")
        self.assertEqual(feature.location.start, 24)
        self.assertEqual(feature.location.end, 114)
        self.assertEqual(feature.qualifiers["note"], "Alpha-1")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[7]
        self.assertEqual(feature.type, "REGION")
        self.assertEqual(feature.location.start, 114)
        self.assertEqual(feature.location.end, 206)
        self.assertEqual(feature.qualifiers["note"], "Alpha-2")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[8]
        self.assertEqual(feature.type, "REGION")
        self.assertEqual(feature.location.start, 206)
        self.assertEqual(feature.location.end, 298)
        self.assertEqual(feature.qualifiers["note"], "Alpha-3")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[9]
        self.assertEqual(feature.type, "REGION")
        self.assertEqual(feature.location.start, 298)
        self.assertEqual(feature.location.end, 308)
        self.assertEqual(feature.qualifiers["note"], "Connecting peptide")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[10]
        self.assertEqual(feature.type, "REGION")
        self.assertEqual(feature.location.start, 338)
        self.assertEqual(feature.location.end, 365)
        self.assertEqual(feature.qualifiers["note"], "Disordered")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000256|SAM:MobiDB-lite")
        feature = record.features[11]
        self.assertEqual(feature.type, "COMPBIAS")
        self.assertEqual(feature.location.start, 340)
        self.assertEqual(feature.location.end, 359)
        self.assertEqual(feature.qualifiers["note"], "Polar residues")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000256|SAM:MobiDB-lite")
        feature = record.features[12]
        self.assertEqual(feature.type, "BINDING")
        self.assertEqual(feature.location.start, 30)
        self.assertEqual(feature.location.end, 31)
        self.assertEqual(feature.qualifiers["ligand"], "a peptide antigen")
        self.assertEqual(feature.qualifiers["ligand_id"], "ChEBI:CHEBI:166823")
        self.assertEqual(feature.qualifiers["ligand_label"], "1")
        self.assertEqual(
            feature.qualifiers["ligand_note"], "pathogen-derived peptide antigen"
        )
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000269|PubMed:21943705")
        feature = record.features[13]
        self.assertEqual(feature.type, "BINDING")
        self.assertEqual(feature.location.start, 96)
        self.assertEqual(feature.location.end, 97)
        self.assertEqual(feature.qualifiers["ligand"], "a peptide antigen")
        self.assertEqual(feature.qualifiers["ligand_id"], "ChEBI:CHEBI:166823")
        self.assertEqual(feature.qualifiers["ligand_label"], "1")
        self.assertEqual(
            feature.qualifiers["ligand_note"], "pathogen-derived peptide antigen"
        )
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000269|PubMed:21943705")
        feature = record.features[14]
        self.assertEqual(feature.type, "BINDING")
        self.assertEqual(feature.location.start, 107)
        self.assertEqual(feature.location.end, 108)
        self.assertEqual(feature.qualifiers["ligand"], "a peptide antigen")
        self.assertEqual(feature.qualifiers["ligand_id"], "ChEBI:CHEBI:166823")
        self.assertEqual(feature.qualifiers["ligand_label"], "1")
        self.assertEqual(
            feature.qualifiers["ligand_note"], "pathogen-derived peptide antigen"
        )
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000269|PubMed:21943705")
        feature = record.features[15]
        self.assertEqual(feature.type, "BINDING")
        self.assertEqual(feature.location.start, 139)
        self.assertEqual(feature.location.end, 140)
        self.assertEqual(feature.qualifiers["ligand"], "a peptide antigen")
        self.assertEqual(feature.qualifiers["ligand_id"], "ChEBI:CHEBI:166823")
        self.assertEqual(feature.qualifiers["ligand_label"], "2")
        self.assertEqual(feature.qualifiers["ligand_note"], "self-peptide antigen")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000269|PubMed:21543847")
        feature = record.features[16]
        self.assertEqual(feature.type, "BINDING")
        self.assertEqual(feature.location.start, 166)
        self.assertEqual(feature.location.end, 167)
        self.assertEqual(feature.qualifiers["ligand"], "a peptide antigen")
        self.assertEqual(feature.qualifiers["ligand_id"], "ChEBI:CHEBI:166823")
        self.assertEqual(feature.qualifiers["ligand_label"], "1")
        self.assertEqual(
            feature.qualifiers["ligand_note"], "pathogen-derived peptide antigen"
        )
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000269|PubMed:21943705")
        feature = record.features[17]
        self.assertEqual(feature.type, "BINDING")
        self.assertEqual(feature.location.start, 169)
        self.assertEqual(feature.location.end, 170)
        self.assertEqual(feature.qualifiers["ligand"], "a peptide antigen")
        self.assertEqual(feature.qualifiers["ligand_id"], "ChEBI:CHEBI:166823")
        self.assertEqual(feature.qualifiers["ligand_label"], "1")
        self.assertEqual(
            feature.qualifiers["ligand_note"], "pathogen-derived peptide antigen"
        )
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000269|PubMed:21943705")
        feature = record.features[18]
        self.assertEqual(feature.type, "BINDING")
        self.assertEqual(feature.location.start, 182)
        self.assertEqual(feature.location.end, 183)
        self.assertEqual(feature.qualifiers["ligand"], "a peptide antigen")
        self.assertEqual(feature.qualifiers["ligand_id"], "ChEBI:CHEBI:166823")
        self.assertEqual(feature.qualifiers["ligand_label"], "1")
        self.assertEqual(
            feature.qualifiers["ligand_note"], "pathogen-derived peptide antigen"
        )
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000269|PubMed:21943705")
        feature = record.features[19]
        self.assertEqual(feature.type, "BINDING")
        self.assertEqual(feature.location.start, 182)
        self.assertEqual(feature.location.end, 183)
        self.assertEqual(feature.qualifiers["ligand"], "a peptide antigen")
        self.assertEqual(feature.qualifiers["ligand_id"], "ChEBI:CHEBI:166823")
        self.assertEqual(feature.qualifiers["ligand_label"], "2")
        self.assertEqual(feature.qualifiers["ligand_note"], "self-peptide antigen")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000269|PubMed:21543847")
        feature = record.features[20]
        self.assertEqual(feature.type, "BINDING")
        self.assertEqual(feature.location.start, 194)
        self.assertEqual(feature.location.end, 195)
        self.assertEqual(feature.qualifiers["ligand"], "a peptide antigen")
        self.assertEqual(feature.qualifiers["ligand_id"], "ChEBI:CHEBI:166823")
        self.assertEqual(feature.qualifiers["ligand_label"], "1")
        self.assertEqual(
            feature.qualifiers["ligand_note"], "pathogen-derived peptide antigen"
        )
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000269|PubMed:21943705")
        feature = record.features[21]
        self.assertEqual(feature.type, "MOD_RES")
        self.assertEqual(feature.location.start, 82)
        self.assertEqual(feature.location.end, 83)
        self.assertEqual(feature.qualifiers["note"], "Sulfotyrosine")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[22]
        self.assertEqual(feature.type, "MOD_RES")
        self.assertEqual(feature.location.start, 342)
        self.assertEqual(feature.location.end, 343)
        self.assertEqual(feature.qualifiers["note"], "Phosphoserine")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007744|PubMed:24275569")
        feature = record.features[23]
        self.assertEqual(feature.type, "MOD_RES")
        self.assertEqual(feature.location.start, 343)
        self.assertEqual(feature.location.end, 344)
        self.assertEqual(feature.qualifiers["note"], "Phosphotyrosine")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007744|PubMed:24275569")
        feature = record.features[24]
        self.assertEqual(feature.type, "MOD_RES")
        self.assertEqual(feature.location.start, 348)
        self.assertEqual(feature.location.end, 349)
        self.assertEqual(feature.qualifiers["note"], "Phosphoserine")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007744|PubMed:24275569")
        feature = record.features[25]
        self.assertEqual(feature.type, "MOD_RES")
        self.assertEqual(feature.location.start, 349)
        self.assertEqual(feature.location.end, 350)
        self.assertEqual(feature.qualifiers["note"], "Phosphoserine")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007744|PubMed:24275569")
        feature = record.features[26]
        self.assertEqual(feature.type, "MOD_RES")
        self.assertEqual(feature.location.start, 351)
        self.assertEqual(feature.location.end, 352)
        self.assertEqual(feature.qualifiers["note"], "Phosphoserine")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007744|PubMed:24275569")
        feature = record.features[27]
        self.assertEqual(feature.type, "MOD_RES")
        self.assertEqual(feature.location.start, 355)
        self.assertEqual(feature.location.end, 356)
        self.assertEqual(feature.qualifiers["note"], "Phosphoserine")
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0007744|PubMed:23186163,ECO:0007744|PubMed:24275569",
        )
        feature = record.features[28]
        self.assertEqual(feature.type, "MOD_RES")
        self.assertEqual(feature.location.start, 358)
        self.assertEqual(feature.location.end, 359)
        self.assertEqual(feature.qualifiers["note"], "Phosphoserine")
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0007744|PubMed:23186163,ECO:0007744|PubMed:24275569",
        )
        feature = record.features[29]
        self.assertEqual(feature.type, "CARBOHYD")
        self.assertEqual(feature.location.start, 109)
        self.assertEqual(feature.location.end, 110)
        self.assertEqual(feature.qualifiers["note"], "N-linked (GlcNAc...) asparagine")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000269|PubMed:19159218")
        feature = record.features[30]
        self.assertEqual(feature.type, "DISULFID")
        self.assertEqual(feature.location.start, 124)
        self.assertEqual(feature.location.end, 188)
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000255|PROSITE-ProRule:PRU00114,ECO:0000269|PubMed:16041067, ECO:0000269|PubMed:19177349,ECO:0000269|PubMed:20844028, ECO:0000269|PubMed:21543847,ECO:0000269|PubMed:21943705, ECO:0000269|PubMed:26758806,ECO:0000269|PubMed:28250417, ECO:0000269|PubMed:7694806",
        )
        feature = record.features[31]
        self.assertEqual(feature.type, "DISULFID")
        self.assertEqual(feature.location.start, 226)
        self.assertEqual(feature.location.end, 283)
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000255|PROSITE-ProRule:PRU00114,ECO:0000269|PubMed:16041067, ECO:0000269|PubMed:19177349,ECO:0000269|PubMed:20844028, ECO:0000269|PubMed:21543847,ECO:0000269|PubMed:21943705, ECO:0000269|PubMed:26758806,ECO:0000269|PubMed:28250417, ECO:0000269|PubMed:7694806",
        )
        feature = record.features[32]
        self.assertEqual(feature.type, "VAR_SEQ")
        self.assertEqual(feature.location.start, 175)
        self.assertEqual(feature.location.end, 187)
        self.assertEqual(
            feature.qualifiers["note"], "EAEQLRAYLDGT -> AAEQQRAYLEGR (in isoform 2)"
        )
        feature = record.features[33]
        self.assertEqual(feature.type, "VAR_SEQ")
        self.assertEqual(feature.location.start, 336)
        self.assertEqual(feature.location.end, 337)
        self.assertEqual(feature.qualifiers["note"], "S -> SGGEGVK (in isoform 2)")
        feature = record.features[34]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 2)
        self.assertEqual(feature.location.end, 3)
        self.assertEqual(feature.qualifiers["note"], "V -> I (in allele A*34:01)")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000269|PubMed:1431115")
        feature = record.features[35]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 4)
        self.assertEqual(feature.location.end, 5)
        self.assertEqual(feature.qualifiers["note"], "A -> P (in allele A*80:01)")
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:8188325,ECO:0000269|PubMed:8284791",
        )
        feature = record.features[36]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 9)
        self.assertEqual(feature.location.end, 10)
        self.assertEqual(
            feature.qualifiers["note"],
            "L -> V (in allele A*02:01, allele A*02:05, allele A*23:01, allele A*24:02, allele A*25:01, allele A*26:01, allele A*34:01, allele A*43:01, allele A*66:01, allele A*68:01 and allele A*69:01)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:1729171, ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2982951, ECO:0000269|PubMed:3496393,ECO:0000269|PubMed:3877632, ECO:0000269|PubMed:7836067,ECO:0000269|PubMed:8026990, ECO:0000269|PubMed:8475492,ECO:0000269|PubMed:9349616",
        )
        feature = record.features[37]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 13)
        self.assertEqual(feature.location.end, 14)
        self.assertEqual(
            feature.qualifiers["note"],
            "S -> L (in allele A*29:02, allele A*31:01, allele A*32:01, allele A*33:01 and allele A*74:01)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:1782566,ECO:0000269|PubMed:2431040, ECO:0000269|PubMed:2478623,ECO:0000269|PubMed:8795145, ECO:0000269|Ref.29",
        )
        feature = record.features[38]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 22)
        self.assertEqual(feature.location.end, 23)
        self.assertEqual(feature.qualifiers["note"], "W -> R (in allele A*74:01)")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000269|PubMed:1431115")
        feature = record.features[39]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 32)
        self.assertEqual(feature.location.end, 33)
        self.assertEqual(
            feature.qualifiers["note"],
            "F -> S (in allele A*23:01, allele A*24:02 and allele A*30:01; dbSNP:rs2075684)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:2478623, ECO:0000269|PubMed:7871528,ECO:0000269|PubMed:9349616",
        )
        feature = record.features[40]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 32)
        self.assertEqual(feature.location.end, 33)
        self.assertEqual(
            feature.qualifiers["note"],
            "F -> T (in allele A*29:02, allele A*31:01 and allele A*33:01; requires 2 nucleotide substitutions)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1782566, ECO:0000269|PubMed:2478623,ECO:0000269|PubMed:8795145",
        )
        feature = record.features[41]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 32)
        self.assertEqual(feature.location.end, 33)
        self.assertEqual(
            feature.qualifiers["note"],
            "F -> Y (in allele A*02:05, allele A*11:01, allele A*25:01, allele A*26:01, allele A*34:01, allele A*43:01, allele A*66:01, allele A*68:01 and allele A*69:01; dbSNP:rs2075684)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:2437024,ECO:0000269|PubMed:2460344, ECO:0000269|PubMed:3496393,ECO:0000269|PubMed:3877632, ECO:0000269|PubMed:8016845,ECO:0000269|PubMed:8026990, ECO:0000269|PubMed:8475492",
        )
        feature = record.features[42]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 40)
        self.assertEqual(feature.location.end, 41)
        self.assertEqual(
            feature.qualifiers["note"], "R -> S (in allele A*30:01; dbSNP:rs1059423)"
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:2478623,ECO:0000269|PubMed:7871528",
        )
        feature = record.features[43]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 54)
        self.assertEqual(feature.location.end, 55)
        self.assertEqual(feature.qualifiers["note"], "T -> S (in allele A*80:01)")
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:8188325,ECO:0000269|PubMed:8284791",
        )
        feature = record.features[44]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 58)
        self.assertEqual(feature.location.end, 59)
        self.assertEqual(feature.qualifiers["note"], "R -> Q (in allele A*80:01)")
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:8188325,ECO:0000269|PubMed:8284791",
        )
        feature = record.features[45]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 66)
        self.assertEqual(feature.location.end, 67)
        self.assertEqual(
            feature.qualifiers["note"], "Q -> R (in allele A*02:05; dbSNP:rs41559117)"
        )
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000269|PubMed:3496393")
        feature = record.features[46]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 67)
        self.assertEqual(feature.location.end, 68)
        self.assertEqual(
            feature.qualifiers["note"], "R -> K (in alleles A*01:01 and allele A*36:01)"
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:2251137, ECO:0000269|PubMed:2715640,ECO:0000269|PubMed:9349617",
        )
        feature = record.features[47]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 79)
        self.assertEqual(feature.location.end, 80)
        self.assertEqual(feature.qualifiers["note"], "G -> E (in allele A*80:01)")
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:8188325,ECO:0000269|PubMed:8284791",
        )
        feature = record.features[48]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 79)
        self.assertEqual(feature.location.end, 80)
        self.assertEqual(
            feature.qualifiers["note"],
            "G -> R (in allele A*30:01 and allele A*31:01; dbSNP:rs1059449)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:2478623, ECO:0000269|PubMed:7871528,ECO:0000269|PubMed:8795145",
        )
        feature = record.features[49]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 85)
        self.assertEqual(feature.location.end, 86)
        self.assertEqual(
            feature.qualifiers["note"],
            "Q -> E (in allele A*23:01, allele 24:02 and allele A*80:01)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:8188325, ECO:0000269|PubMed:8284791,ECO:0000269|PubMed:9349616",
        )
        feature = record.features[50]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 85)
        self.assertEqual(feature.location.end, 86)
        self.assertEqual(
            feature.qualifiers["note"],
            "Q -> G (in allele A*02:01 and allele A*02:05; requires 2 nucleotide substitutions)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2982951, ECO:0000269|PubMed:3496393,ECO:0000269|PubMed:7836067",
        )
        feature = record.features[51]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 85)
        self.assertEqual(feature.location.end, 86)
        self.assertEqual(
            feature.qualifiers["note"], "Q -> L (in alleles A*29:02 and allele A*43:01)"
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:1782566",
        )
        feature = record.features[52]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 85)
        self.assertEqual(feature.location.end, 86)
        self.assertEqual(
            feature.qualifiers["note"],
            "Q -> R (in allele A*25:01, allele A*26:01, allele A*33:01, allele A*34:01, allele A*66:01, allele A*68:01 and allele A*69:01)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:2478623,ECO:0000269|PubMed:3877632, ECO:0000269|PubMed:8026990,ECO:0000269|PubMed:8475492",
        )
        feature = record.features[53]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 86)
        self.assertEqual(feature.location.end, 87)
        self.assertEqual(
            feature.qualifiers["note"],
            "E -> N (in alleles A*25:01, allele A*26:01, allele A*33:01, allele A*34:01, allele A*66:01, allele A*68:01 and allele A*69:01; requires 2 nucleotide substitutions)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:2478623,ECO:0000269|PubMed:3877632, ECO:0000269|PubMed:8026990,ECO:0000269|PubMed:8475492",
        )
        feature = record.features[54]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 86)
        self.assertEqual(feature.location.end, 87)
        self.assertEqual(
            feature.qualifiers["note"], "E -> Q (in allele A*29:02 and allele A*43:01)"
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:1782566",
        )
        feature = record.features[55]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 88)
        self.assertEqual(feature.location.end, 89)
        self.assertEqual(
            feature.qualifiers["note"],
            "R -> G (in allele A*23:01 and allele 24:02; dbSNP:rs199474430)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:9349616",
        )
        feature = record.features[56]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 89)
        self.assertEqual(feature.location.end, 90)
        self.assertEqual(
            feature.qualifiers["note"],
            "N -> K (in allele A*02:01, allele A*02:05, allele A*23:01, allele 24:02 and allele A*34:01; dbSNP:rs199474436)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:1729171, ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2982951, ECO:0000269|PubMed:3496393,ECO:0000269|PubMed:7836067, ECO:0000269|PubMed:9349616",
        )
        feature = record.features[57]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 90)
        self.assertEqual(feature.location.end, 91)
        self.assertEqual(
            feature.qualifiers["note"], "V -> M (in allele A*01:01 and allele A*36:01)"
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:2251137, ECO:0000269|PubMed:2715640,ECO:0000269|PubMed:9349617",
        )
        feature = record.features[58]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 93)
        self.assertEqual(feature.location.end, 94)
        self.assertEqual(
            feature.qualifiers["note"],
            "Q -> H (in allele A*01:01, allele A*02:01, allele A*02:05, allele A*23:01, allele 24:02, allele A*25:01, allele A*26:01, allele A*31:01, allele A*32:01, allele A*33:01, allele A*36:01, allele A*43:01, allele A*74:01 and allele A*80:01; dbSNP:rs78306866)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:2251137, ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2431040, ECO:0000269|PubMed:2478623,ECO:0000269|PubMed:2715640, ECO:0000269|PubMed:2982951,ECO:0000269|PubMed:3496393, ECO:0000269|PubMed:7836067,ECO:0000269|PubMed:8026990, ECO:0000269|PubMed:8188325,ECO:0000269|PubMed:8284791, ECO:0000269|PubMed:8475492,ECO:0000269|PubMed:8795145, ECO:0000269|PubMed:9349616,ECO:0000269|PubMed:9349617, ECO:0000269|Ref.29",
        )
        feature = record.features[59]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 96)
        self.assertEqual(feature.location.end, 97)
        self.assertEqual(
            feature.qualifiers["note"],
            "T -> I (in allele A*31:01 and allele A*33:01; dbSNP:rs199474457)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:2478623, ECO:0000269|PubMed:8795145",
        )
        feature = record.features[60]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 97)
        self.assertEqual(feature.location.end, 98)
        self.assertEqual(
            feature.qualifiers["note"], "D -> H (in allele A*02:01 and allele A*02:05)"
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2982951, ECO:0000269|PubMed:3496393,ECO:0000269|PubMed:7836067",
        )
        feature = record.features[61]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 97)
        self.assertEqual(feature.location.end, 98)
        self.assertEqual(feature.qualifiers["note"], "D -> N (in allele A*80:01)")
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:8188325,ECO:0000269|PubMed:8284791",
        )
        feature = record.features[62]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 99)
        self.assertEqual(feature.location.end, 100)
        self.assertEqual(
            feature.qualifiers["note"],
            "V -> A (in allele A*01:01, allele A*26:01, allele A*29:02, allele A*36:01, allele A*43:01 and allele A*80:01)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:1782566, ECO:0000269|PubMed:2251137,ECO:0000269|PubMed:2715640, ECO:0000269|PubMed:8026990,ECO:0000269|PubMed:8188325, ECO:0000269|PubMed:8284791,ECO:0000269|PubMed:8475492, ECO:0000269|PubMed:9349617",
        )
        feature = record.features[63]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 99)
        self.assertEqual(feature.location.end, 100)
        self.assertEqual(
            feature.qualifiers["note"],
            "V -> E (in allele A*23:01, allele A*24:02, allele A*25:01 and allele A*32:01; dbSNP:rs1071742)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:2431040,ECO:0000269|PubMed:9349616, ECO:0000269|Ref.29",
        )
        feature = record.features[64]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 100)
        self.assertEqual(feature.location.end, 101)
        self.assertEqual(
            feature.qualifiers["note"],
            "D -> N (allele A*01:01, allele A*23:01, allele A*24:02, allele A*26:01, allele A*29:02, allele A*36:01, allele A*43:01 and allele A*80:01; dbSNP:rs1136688)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:1729171, ECO:0000269|PubMed:1782566,ECO:0000269|PubMed:2251137, ECO:0000269|PubMed:2715640,ECO:0000269|PubMed:8026990, ECO:0000269|PubMed:8188325,ECO:0000269|PubMed:8284791, ECO:0000269|PubMed:8475492,ECO:0000269|PubMed:9349616, ECO:0000269|PubMed:9349617",
        )
        feature = record.features[65]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 100)
        self.assertEqual(feature.location.end, 101)
        self.assertEqual(
            feature.qualifiers["note"],
            "D -> S (in allele A*25:01 and allele A*32:01; requires 2 nucleotide substitutions)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2431040, ECO:0000269|Ref.29",
        )
        feature = record.features[66]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 102)
        self.assertEqual(feature.location.end, 107)
        self.assertEqual(
            feature.qualifiers["note"],
            "GTLRG -> RIALR (in allele A*23:01, allele A*24:02, allele A*25:01 and allele A*32:01; Bw4 motif RIALR is involved in the recognition of NK cell inhibitory receptor KIR3DL1)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:17182537,ECO:0000269|PubMed:1729171, ECO:0000269|PubMed:18502829,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:2431040,ECO:0000269|PubMed:9349616, ECO:0000269|Ref.29",
        )
        feature = record.features[67]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 113)
        self.assertEqual(feature.location.end, 114)
        self.assertEqual(
            feature.qualifiers["note"],
            "A -> D (in allele A*01:01, allele A*11:01, allele A*25:01, allele A*26:01, allele A*34:01, allele A*36:01, allele A*43:01, allele A*66:01 and allele A*80:01; dbSNP:rs1136692)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:2251137, ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2437024, ECO:0000269|PubMed:2460344,ECO:0000269|PubMed:2715640, ECO:0000269|PubMed:8016845,ECO:0000269|PubMed:8026990, ECO:0000269|PubMed:8188325,ECO:0000269|PubMed:8284791, ECO:0000269|PubMed:8475492,ECO:0000269|PubMed:9349617",
        )
        feature = record.features[68]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 118)
        self.assertEqual(feature.location.end, 119)
        self.assertEqual(
            feature.qualifiers["note"],
            "I -> L (in allele A*02:05, allele A*23:01 and allele 24:02; dbSNP:rs1071743)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:3496393, ECO:0000269|PubMed:9349616",
        )
        feature = record.features[69]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 118)
        self.assertEqual(feature.location.end, 119)
        self.assertEqual(
            feature.qualifiers["note"], "I -> V (in allele A*02:01 and allele A*69:01)"
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2982951, ECO:0000269|PubMed:3877632,ECO:0000269|PubMed:7836067",
        )
        feature = record.features[70]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 120)
        self.assertEqual(feature.location.end, 121)
        self.assertEqual(
            feature.qualifiers["note"],
            "I -> M (in allele A*23:01, allele 24:02, allele A*29:02, allele A*31:01, allele A*32:01, allele A*33:01, allele A*68:01 and allele A*74:01; dbSNP:rs1136695)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:1782566, ECO:0000269|PubMed:2431040,ECO:0000269|PubMed:2478623, ECO:0000269|PubMed:3877632,ECO:0000269|PubMed:8795145, ECO:0000269|PubMed:9349616,ECO:0000269|Ref.29",
        )
        feature = record.features[71]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 120)
        self.assertEqual(feature.location.end, 121)
        self.assertEqual(
            feature.qualifiers["note"],
            "I -> R (in allele A*02:01, allele A*02:05, allele A*25:01, allele A*26:01, allele A*34:01, allele A*43:01, allele A*66:01 and allele A*69:01)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:2982951,ECO:0000269|PubMed:3496393, ECO:0000269|PubMed:3877632,ECO:0000269|PubMed:7836067, ECO:0000269|PubMed:8026990,ECO:0000269|PubMed:8475492",
        )
        feature = record.features[72]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 122)
        self.assertEqual(feature.location.end, 123)
        self.assertEqual(
            feature.qualifiers["note"],
            "Y -> F (in allele A*23:01, allele 24:02; dbSNP:rs1136697)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:9349616",
        )
        feature = record.features[73]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 128)
        self.assertEqual(feature.location.end, 129)
        self.assertEqual(
            feature.qualifiers["note"],
            "S -> P (in allele A*01:01, allele A*11:01, allele A*25:01, allele A*26:01, allele A*32:01, allele A*34:01, allele A*36:01, allele A*43:01, allele A*66:01 and allele A*74:01; dbSNP:rs1136700)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:2251137, ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2431040, ECO:0000269|PubMed:2437024,ECO:0000269|PubMed:2460344, ECO:0000269|PubMed:2715640,ECO:0000269|PubMed:8016845, ECO:0000269|PubMed:8026990,ECO:0000269|PubMed:8475492, ECO:0000269|PubMed:9349617,ECO:0000269|Ref.29",
        )
        feature = record.features[74]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 130)
        self.assertEqual(feature.location.end, 131)
        self.assertEqual(
            feature.qualifiers["note"],
            "G -> W (in allele A*02:01, allele A*02:05 and allele A*69:01; dbSNP:rs1136702)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2982951, ECO:0000269|PubMed:3496393,ECO:0000269|PubMed:3877632, ECO:0000269|PubMed:7836067",
        )
        feature = record.features[75]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 132)
        self.assertEqual(feature.location.end, 133)
        self.assertEqual(
            feature.qualifiers["note"],
            "F -> L (in allele A*32:01 and allele A*74:01; dbSNP:rs1059488)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:2431040, ECO:0000269|Ref.29",
        )
        feature = record.features[76]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 137)
        self.assertEqual(feature.location.end, 138)
        self.assertEqual(
            feature.qualifiers["note"],
            "R -> E (in allele A*30:01; requires 2 nucleotide substitutions)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:2478623,ECO:0000269|PubMed:7871528",
        )
        feature = record.features[77]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 137)
        self.assertEqual(feature.location.end, 138)
        self.assertEqual(
            feature.qualifiers["note"],
            "R -> H (in allele A*02:01, allele A*02:05, allele A*23:01, allele A*24:02, allele A*69:01)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:2982951,ECO:0000269|PubMed:3496393, ECO:0000269|PubMed:3877632,ECO:0000269|PubMed:7836067, ECO:0000269|PubMed:9349616",
        )
        feature = record.features[78]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 137)
        self.assertEqual(feature.location.end, 138)
        self.assertEqual(
            feature.qualifiers["note"],
            "R -> Q (in allele A*25:01, allele A*26:01, allele A*31:01, allele A*32:01, allele A*33:01, allele A*34:01, allele A*43:01, allele A*66:01, allele A*74:01)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2431040, ECO:0000269|PubMed:2478623,ECO:0000269|PubMed:8026990, ECO:0000269|PubMed:8475492,ECO:0000269|PubMed:8795145, ECO:0000269|Ref.29",
        )
        feature = record.features[79]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 139)
        self.assertEqual(feature.location.end, 140)
        self.assertEqual(feature.qualifiers["note"], "D -> H (in allele A*30:01)")
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:2478623,ECO:0000269|PubMed:7871528",
        )
        feature = record.features[80]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 139)
        self.assertEqual(feature.location.end, 140)
        self.assertEqual(
            feature.qualifiers["note"],
            "D -> Y (in allele A*02:01, allele A*02:05, allele A*23:01, allele A*24:02 and allele A*69:01)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:2982951,ECO:0000269|PubMed:3496393, ECO:0000269|PubMed:3877632,ECO:0000269|PubMed:7836067, ECO:0000269|PubMed:9349616",
        )
        feature = record.features[81]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 150)
        self.assertEqual(feature.location.end, 151)
        self.assertEqual(
            feature.qualifiers["note"],
            "N -> K (in allele A*02:01, allele A*02:05, allele A*23:01, allele A*24:02, allele A*68:01 and allele A*69:01; dbSNP:rs1059509)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:2982951,ECO:0000269|PubMed:3496393, ECO:0000269|PubMed:3877632,ECO:0000269|PubMed:7836067, ECO:0000269|PubMed:9349616",
        )
        feature = record.features[82]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 165)
        self.assertEqual(feature.location.end, 166)
        self.assertEqual(
            feature.qualifiers["note"],
            "I -> T (in allele A*02:01, allele A*02:05, allele A*68:01 and allele A*69:01; dbSNP:rs1059516)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2982951, ECO:0000269|PubMed:3496393,ECO:0000269|PubMed:3877632, ECO:0000269|PubMed:7836067",
        )
        feature = record.features[83]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 167)
        self.assertEqual(feature.location.end, 168)
        self.assertEqual(
            feature.qualifiers["note"],
            "K -> Q (in allele A*23:01, allele A*25:01, allele A*26:01, allele A*29:02, allele A*30:01, allele A*31:01, allele A*32:01, allele A*33:01, allele A*34:01, allele A*43:01, allele A*66:01 and allele A*74:01)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:1782566, ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2431040, ECO:0000269|PubMed:2478623,ECO:0000269|PubMed:7871528, ECO:0000269|PubMed:8026990,ECO:0000269|PubMed:8475492, ECO:0000269|PubMed:8795145,ECO:0000269|Ref.29",
        )
        feature = record.features[84]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 168)
        self.assertEqual(feature.location.end, 169)
        self.assertEqual(
            feature.qualifiers["note"],
            "R -> H (in allele A*02:01, allele A*02:05, allele A*68:01 and allele A*69:01; dbSNP:rs1059520)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2982951, ECO:0000269|PubMed:3496393,ECO:0000269|PubMed:3877632, ECO:0000269|PubMed:7836067",
        )
        feature = record.features[85]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 172)
        self.assertEqual(feature.location.end, 173)
        self.assertEqual(
            feature.qualifiers["note"],
            "A -> T (in allele A*25:01, allele A*26:01, allele A*34:01, allele A*43:01 and allele A*66:01; dbSNP:rs1059526)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:8026990,ECO:0000269|PubMed:8475492",
        )
        feature = record.features[86]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 173)
        self.assertEqual(feature.location.end, 174)
        self.assertEqual(
            feature.qualifiers["note"], "A -> V (in allele A*01:01 and allele A*36:01)"
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:2251137, ECO:0000269|PubMed:2715640,ECO:0000269|PubMed:9349617",
        )
        feature = record.features[87]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 174)
        self.assertEqual(feature.location.end, 175)
        self.assertEqual(
            feature.qualifiers["note"],
            "H -> R (in allele A*23:01, allele A*29:02, allele A*30:01, allele A*31:01, allele A*32:01, allele A*33:01, allele A*74:01 and allele A*80:01)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:1782566, ECO:0000269|PubMed:2431040,ECO:0000269|PubMed:2478623, ECO:0000269|PubMed:7871528,ECO:0000269|PubMed:8188325, ECO:0000269|PubMed:8284791,ECO:0000269|PubMed:8795145, ECO:0000269|Ref.29",
        )
        feature = record.features[88]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 175)
        self.assertEqual(feature.location.end, 176)
        self.assertEqual(
            feature.qualifiers["note"],
            "E -> A (in allele A*01:01, allele A*11:01 and allele A*36:01)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:2251137, ECO:0000269|PubMed:2437024,ECO:0000269|PubMed:2460344, ECO:0000269|PubMed:2715640,ECO:0000269|PubMed:8016845, ECO:0000269|PubMed:9349617",
        )
        feature = record.features[89]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 175)
        self.assertEqual(feature.location.end, 176)
        self.assertEqual(
            feature.qualifiers["note"],
            "E -> R (in allele A*80:01; requires 2 nucleotide substitutions)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:8188325,ECO:0000269|PubMed:8284791",
        )
        feature = record.features[90]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 175)
        self.assertEqual(feature.location.end, 176)
        self.assertEqual(
            feature.qualifiers["note"],
            "E -> V (in allele A*02:01, allele A*02:05, allele A*23:01, allele A*24:02, allele A*29:02, allele A*31:01, allele A*32:01, allele A*33:01, allele A*68:01, allele A*69:01 and allele A*74:01; results in inefficient T cell recognition of epitopes derived from influenza A virus.; dbSNP:rs9256983)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:1782566, ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2431040, ECO:0000269|PubMed:2456340,ECO:0000269|PubMed:2478623, ECO:0000269|PubMed:2982951,ECO:0000269|PubMed:3496393, ECO:0000269|PubMed:3877632,ECO:0000269|PubMed:7836067, ECO:0000269|PubMed:8795145,ECO:0000269|PubMed:9349616, ECO:0000269|Ref.29",
        )
        feature = record.features[91]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 175)
        self.assertEqual(feature.location.end, 176)
        self.assertEqual(
            feature.qualifiers["note"],
            "E -> W (in allele A*30:01; requires 2 nucleotide substitutions)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:2478623,ECO:0000269|PubMed:7871528",
        )
        feature = record.features[92]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 179)
        self.assertEqual(feature.location.end, 180)
        self.assertEqual(
            feature.qualifiers["note"], "L -> Q (in allele A*11:01 and allele A*24:02)"
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:2437024, ECO:0000269|PubMed:2460344,ECO:0000269|PubMed:8016845, ECO:0000269|PubMed:9349616",
        )
        feature = record.features[93]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 179)
        self.assertEqual(feature.location.end, 180)
        self.assertEqual(
            feature.qualifiers["note"], "L -> R (in allele A*01:01 and allele A*36:01)"
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:2251137, ECO:0000269|PubMed:2715640,ECO:0000269|PubMed:9349617",
        )
        feature = record.features[94]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 179)
        self.assertEqual(feature.location.end, 180)
        self.assertEqual(
            feature.qualifiers["note"],
            "L -> W (in allele A*02:05, allele A*25:01, allele A*26:01, allele A*34:01, allele A*43:01, allele A*66:01, allele A*68:01; dbSNP:rs9260156)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:3496393,ECO:0000269|PubMed:3877632, ECO:0000269|PubMed:8026990,ECO:0000269|PubMed:8475492",
        )
        feature = record.features[95]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 181)
        self.assertEqual(feature.location.end, 182)
        self.assertEqual(
            feature.qualifiers["note"], "A -> V (in allele A*01:01 and allele A*36:01)"
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:2251137, ECO:0000269|PubMed:2715640,ECO:0000269|PubMed:9349617",
        )
        feature = record.features[96]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 184)
        self.assertEqual(feature.location.end, 185)
        self.assertEqual(
            feature.qualifiers["note"],
            "D -> E (in allele A*01:01, allele A*02:01, allele A*02:05, allele A*11:01, allele A*23:01, allele A*24:02, allele A*25:01, allele A*26:01, allele A*29:02, allele A*30:01, allele A*31:01, allele A*32:01, allele A*33:01, allele A*34:01, allele A*36:01, allele A*43:01, allele A*66:01, allele A*68:01, allele A*69:01, allele A*74:01 and allele A*80:01; dbSNP:rs1059542)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:1782566, ECO:0000269|PubMed:2251137,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:2431040,ECO:0000269|PubMed:2437024, ECO:0000269|PubMed:2460344,ECO:0000269|PubMed:2478623, ECO:0000269|PubMed:2715640,ECO:0000269|PubMed:2982951, ECO:0000269|PubMed:3496393,ECO:0000269|PubMed:3877632, ECO:0000269|PubMed:7836067,ECO:0000269|PubMed:7871528, ECO:0000269|PubMed:8016845,ECO:0000269|PubMed:8026990, ECO:0000269|PubMed:8188325,ECO:0000269|PubMed:8284791, ECO:0000269|PubMed:8475492,ECO:0000269|PubMed:8795145, ECO:0000269|PubMed:9349616,ECO:0000269|PubMed:9349617, ECO:0000269|Ref.29",
        )
        feature = record.features[97]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 186)
        self.assertEqual(feature.location.end, 187)
        self.assertEqual(
            feature.qualifiers["note"],
            "T -> E (in allele A*80:01; requires 2 nucleotide substitutions)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:8188325,ECO:0000269|PubMed:8284791",
        )
        feature = record.features[98]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 186)
        self.assertEqual(feature.location.end, 187)
        self.assertEqual(
            feature.qualifiers["note"],
            "T -> R (in allele A*01:01, allele A*11:01, allele A*25:01, allele A*26:01, allele A*43:01 and allele A*66:01)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:2251137, ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2437024, ECO:0000269|PubMed:2460344,ECO:0000269|PubMed:2715640, ECO:0000269|PubMed:8016845,ECO:0000269|PubMed:8026990, ECO:0000269|PubMed:8475492,ECO:0000269|PubMed:9349617",
        )
        feature = record.features[99]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 189)
        self.assertEqual(feature.location.end, 190)
        self.assertEqual(
            feature.qualifiers["note"],
            "E -> D (in allele A*01:01, allele A*23:01, allele A*24:02 and allele A*80:01; dbSNP:rs879577815)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:2251137, ECO:0000269|PubMed:2715640,ECO:0000269|PubMed:8188325, ECO:0000269|PubMed:8284791,ECO:0000269|PubMed:9349616, ECO:0000269|PubMed:9349617",
        )
        feature = record.features[100]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 190)
        self.assertEqual(feature.location.end, 191)
        self.assertEqual(
            feature.qualifiers["note"],
            "W -> G (in allele A*01:01, allele A*23:01, allele A*24:02 and allele A*80:01; dbSNP:rs3098019)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:2251137, ECO:0000269|PubMed:2715640,ECO:0000269|PubMed:8188325, ECO:0000269|PubMed:8284791,ECO:0000269|PubMed:9349616, ECO:0000269|PubMed:9349617",
        )
        feature = record.features[101]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 194)
        self.assertEqual(feature.location.end, 195)
        self.assertEqual(feature.qualifiers["note"], "Y -> H (in allele A*33:01)")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000269|PubMed:2478623")
        feature = record.features[102]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 207)
        self.assertEqual(feature.location.end, 208)
        self.assertEqual(
            feature.qualifiers["note"],
            "P -> A (in allele A*02:01, allele A*02:05, allele A*25:01, allele A*26:01, allele A*29:02, allele A*32:01, allele A*34:01, allele A*43:01, allele A*66:01, allele A*68:01, allele A*69:01 and allele A*74:01)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:1782566, ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2431040, ECO:0000269|PubMed:2982951,ECO:0000269|PubMed:3496393, ECO:0000269|PubMed:3877632,ECO:0000269|PubMed:7836067, ECO:0000269|PubMed:8026990,ECO:0000269|PubMed:8475492, ECO:0000269|Ref.29",
        )
        feature = record.features[103]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 209)
        self.assertEqual(feature.location.end, 210)
        self.assertEqual(feature.qualifiers["note"], "K -> R (in allele A*33:01)")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000269|PubMed:2478623")
        feature = record.features[104]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 216)
        self.assertEqual(feature.location.end, 217)
        self.assertEqual(
            feature.qualifiers["note"],
            "P -> A (in allele A*02:01, allele A*02:05, allele A*25:01, allele A*26:01, allele A*29:02, allele A*31:01, allele A*32:01, allele A*33:01, allele A*34:01, allele A*43:01, allele A*66:01, allele A*68:01, allele A*69:01, allele A*74:01)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:1782566,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:2431040,ECO:0000269|PubMed:2478623, ECO:0000269|PubMed:2982951,ECO:0000269|PubMed:3496393, ECO:0000269|PubMed:3877632,ECO:0000269|PubMed:7836067, ECO:0000269|PubMed:8026990,ECO:0000269|PubMed:8475492, ECO:0000269|PubMed:8795145,ECO:0000269|Ref.29",
        )
        feature = record.features[105]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 217)
        self.assertEqual(feature.location.end, 218)
        self.assertEqual(
            feature.qualifiers["note"],
            "I -> V (in allele A*02:01, allele A*02:05, allele A*25:01, allele A*26:01, allele A*29:02, allele A*31:01, allele A*32:01, allele A*33:01, allele A*34:01, allele A*43:01, allele A*66:01, allele A*68:01, allele A*69:01 and allele A*74:01)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:1782566,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:2431040,ECO:0000269|PubMed:2478623, ECO:0000269|PubMed:2982951,ECO:0000269|PubMed:3496393, ECO:0000269|PubMed:3877632,ECO:0000269|PubMed:7836067, ECO:0000269|PubMed:8026990,ECO:0000269|PubMed:8475492, ECO:0000269|PubMed:8795145,ECO:0000269|Ref.29",
        )
        feature = record.features[106]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 230)
        self.assertEqual(feature.location.end, 231)
        self.assertEqual(
            feature.qualifiers["note"],
            "G -> S (in allele A*02:01, allele A*02:05, allele A*25:01, allele A*26:01, allele A*29:02, allele A*31:01, allele A*32:01, allele A*33:01, allele A*34:01, allele A*43:01, allele A*66:01, allele A*68:01, allele A*69:01, allele A*74:01 and allele A*80:01)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:1782566,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:2431040,ECO:0000269|PubMed:2478623, ECO:0000269|PubMed:2982951,ECO:0000269|PubMed:3496393, ECO:0000269|PubMed:3877632,ECO:0000269|PubMed:7836067, ECO:0000269|PubMed:8026990,ECO:0000269|PubMed:8188325, ECO:0000269|PubMed:8284791,ECO:0000269|PubMed:8475492, ECO:0000269|PubMed:8795145,ECO:0000269|Ref.29",
        )
        feature = record.features[107]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 268)
        self.assertEqual(feature.location.end, 269)
        self.assertEqual(
            feature.qualifiers["note"],
            "A -> V (in allele A*68:01; impairs binding to CD8A and reduces recognition by antigen-specific CD8-positive T cells)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:2784196,ECO:0000269|PubMed:3877632",
        )
        feature = record.features[108]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 269)
        self.assertEqual(feature.location.end, 270)
        self.assertEqual(
            feature.qualifiers["note"],
            "A -> S (in allele A*25:01, allele A*26:01, allele A*29:02, allele A*31:01, allele A*32:01, allele A*33:01, allele A*34:01, allele A*43:01, allele A*66:01 and allele A*74:01)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:1782566, ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2431040, ECO:0000269|PubMed:2478623,ECO:0000269|PubMed:8026990, ECO:0000269|PubMed:8475492,ECO:0000269|Ref.29",
        )
        feature = record.features[109]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 276)
        self.assertEqual(feature.location.end, 277)
        self.assertEqual(feature.qualifiers["note"], "E -> K (in allele A*80:01)")
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:8188325,ECO:0000269|PubMed:8284791",
        )
        feature = record.features[110]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 276)
        self.assertEqual(feature.location.end, 277)
        self.assertEqual(
            feature.qualifiers["note"],
            "E -> Q (in allele A*02:01, allele A*02:05, allele A*25:01, allele A*26:01, allele A*29:02, allele A*31:01, allele A*32:01, allele A*33:01, allele A*34:01, allele A*43:01, allele A*66:01, allele A*68:01, allele A*69:01 and allele A*74:01)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:1782566,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:2431040,ECO:0000269|PubMed:2478623, ECO:0000269|PubMed:2982951,ECO:0000269|PubMed:3496393, ECO:0000269|PubMed:3877632,ECO:0000269|PubMed:7836067, ECO:0000269|PubMed:8026990,ECO:0000269|PubMed:8475492, ECO:0000269|PubMed:8795145,ECO:0000269|Ref.29",
        )
        feature = record.features[111]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 278)
        self.assertEqual(feature.location.end, 279)
        self.assertEqual(feature.qualifiers["note"], "Q -> K (in allele A*80:01)")
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:8188325,ECO:0000269|PubMed:8284791",
        )
        feature = record.features[112]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 291)
        self.assertEqual(feature.location.end, 292)
        self.assertEqual(feature.qualifiers["note"], "K -> E (in allele A*80:01)")
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:8188325,ECO:0000269|PubMed:8284791",
        )
        feature = record.features[113]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 299)
        self.assertEqual(feature.location.end, 300)
        self.assertEqual(
            feature.qualifiers["note"],
            "L -> P (in allele A*02:01, allele A*02:05, allele A*23:01, allele A*24:02, allele A*25:01, allele A*26:01, allele A*29:02, allele A*31:01, allele A*32:01, allele A*33:01, allele A*34:01, allele A*43:01, allele A*66:01, allele A*68:01, allele A*69:01, allele A*74:01 and allele A*80:01)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:1782566, ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2431040, ECO:0000269|PubMed:2478623,ECO:0000269|PubMed:2982951, ECO:0000269|PubMed:3496393,ECO:0000269|PubMed:3877632, ECO:0000269|PubMed:7836067,ECO:0000269|PubMed:8026990, ECO:0000269|PubMed:8188325,ECO:0000269|PubMed:8284791, ECO:0000269|PubMed:8475492,ECO:0000269|PubMed:8795145, ECO:0000269|PubMed:9349616,ECO:0000269|Ref.29",
        )
        feature = record.features[114]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 305)
        self.assertEqual(feature.location.end, 306)
        self.assertEqual(
            feature.qualifiers["note"], "I -> V (in allele A*23:01 and allele A*24:02)"
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:9349616",
        )
        feature = record.features[115]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 306)
        self.assertEqual(feature.location.end, 307)
        self.assertEqual(feature.qualifiers["note"], "P -> H (in allele A*23:01)")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000269|PubMed:1729171")
        feature = record.features[116]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 311)
        self.assertEqual(feature.location.end, 312)
        self.assertEqual(feature.qualifiers["note"], "I -> L (in allele A*34:01)")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000269|PubMed:1431115")
        feature = record.features[117]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 317)
        self.assertEqual(feature.location.end, 318)
        self.assertEqual(
            feature.qualifiers["note"],
            "L -> F (in allele A*02:01, allele A*02:05, allele A*25:01, allele A*26:01, allele A*29:02, allele A*31:01, allele A*32:01, allele A*33:01, allele A*34:01, allele A*43:01, allele A*66:01, allele A*68:01, allele A*69:01 and allele A*74:01)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:1782566,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:2431040,ECO:0000269|PubMed:2478623, ECO:0000269|PubMed:2982951,ECO:0000269|PubMed:3496393, ECO:0000269|PubMed:3877632,ECO:0000269|PubMed:7836067, ECO:0000269|PubMed:8026990,ECO:0000269|PubMed:8475492, ECO:0000269|PubMed:8795145,ECO:0000269|Ref.29",
        )
        feature = record.features[118]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 320)
        self.assertEqual(feature.location.end, 321)
        self.assertEqual(
            feature.qualifiers["note"], "V -> M (in allele A*32:01 and allele A*74:01)"
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1431115,ECO:0000269|PubMed:2431040, ECO:0000269|Ref.29",
        )
        feature = record.features[119]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 321)
        self.assertEqual(feature.location.end, 322)
        self.assertEqual(
            feature.qualifiers["note"],
            "I -> F (in allele A*29:02, allele A*31:01, allele A*32:01, allele A*33:01 and allele A*74:01)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:1782566,ECO:0000269|PubMed:2431040, ECO:0000269|PubMed:2478623,ECO:0000269|PubMed:8795145, ECO:0000269|Ref.29",
        )
        feature = record.features[120]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 322)
        self.assertEqual(feature.location.end, 323)
        self.assertEqual(
            feature.qualifiers["note"],
            "T -> A (in allele A*25:01, allele A*26:01, allele A*29:02, allele A*31:01, allele A*32:01, allele A*33:01, allele A*34:01, allele A*43:01, allele A*66:01, allele A*74:01 and allele A*80:01)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:1782566,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:2431040,ECO:0000269|PubMed:2478623, ECO:0000269|PubMed:8026990,ECO:0000269|PubMed:8188325, ECO:0000269|PubMed:8284791,ECO:0000269|PubMed:8475492, ECO:0000269|PubMed:8795145,ECO:0000269|Ref.29",
        )
        feature = record.features[121]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 330)
        self.assertEqual(feature.location.end, 331)
        self.assertEqual(
            feature.qualifiers["note"],
            "M -> R (in allele A*29:02, allele A*31:01, allele A*32:01, allele A*33:01 and allele A*74:01)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:1782566,ECO:0000269|PubMed:2431040, ECO:0000269|PubMed:2478623,ECO:0000269|PubMed:8795145, ECO:0000269|Ref.29",
        )
        feature = record.features[122]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 333)
        self.assertEqual(feature.location.end, 334)
        self.assertEqual(feature.qualifiers["note"], "R -> K (allele A*80:01)")
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:8188325,ECO:0000269|PubMed:8284791",
        )
        feature = record.features[123]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 334)
        self.assertEqual(feature.location.end, 335)
        self.assertEqual(
            feature.qualifiers["note"], "K -> N (in allele A*23:01 and allele A*24:02)"
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:9349616",
        )
        feature = record.features[124]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 337)
        self.assertEqual(feature.location.end, 338)
        self.assertEqual(feature.qualifiers["note"], "D -> V (allele A*80:01)")
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:8188325,ECO:0000269|PubMed:8284791",
        )
        feature = record.features[125]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 344)
        self.assertEqual(feature.location.end, 345)
        self.assertEqual(
            feature.qualifiers["note"],
            "T -> S (in allele A*02:01, allele A*02:05, allele A*23:01, allele A*24:02, allele A*25:01, allele A*26:01, allele A*29:02, allele A*31:01, allele A*32:01, allele A*33:01, allele A*34:01, allele A*43:01, allele A*66:01, allele A*68:01 allele A*69:01, allele A*74:01 and allele A*80:01)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:1729171,ECO:0000269|PubMed:1782566, ECO:0000269|PubMed:2320591,ECO:0000269|PubMed:2431040, ECO:0000269|PubMed:2478623,ECO:0000269|PubMed:2982951, ECO:0000269|PubMed:3496393,ECO:0000269|PubMed:3877632, ECO:0000269|PubMed:7836067,ECO:0000269|PubMed:8026990, ECO:0000269|PubMed:8188325,ECO:0000269|PubMed:8284791, ECO:0000269|PubMed:8475492,ECO:0000269|PubMed:8795145, ECO:0000269|PubMed:9349616,ECO:0000269|Ref.29",
        )
        feature = record.features[126]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 357)
        self.assertEqual(feature.location.end, 358)
        self.assertEqual(
            feature.qualifiers["note"],
            "V -> M (in allele A*25:01, allele A*26:01, allele A*29:02, allele A*31:01, allele A*32:01, allele A*33:01, allele A*34:01, allele A*43:01, allele A*66:01 and allele A*74:01)",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:1317015,ECO:0000269|PubMed:1431115, ECO:0000269|PubMed:1782566,ECO:0000269|PubMed:2320591, ECO:0000269|PubMed:2431040,ECO:0000269|PubMed:2478623, ECO:0000269|PubMed:8026990,ECO:0000269|PubMed:8475492, ECO:0000269|PubMed:8795145,ECO:0000269|Ref.29",
        )
        feature = record.features[127]
        self.assertEqual(feature.type, "MUTAGEN")
        self.assertEqual(feature.location.start, 109)
        self.assertEqual(feature.location.end, 110)
        self.assertEqual(
            feature.qualifiers["note"],
            "N->Q: Impairs the recruitment of HLA-A*02 in the peptide-loading complex.",
        )
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000269|PubMed:21263072")
        feature = record.features[128]
        self.assertEqual(feature.type, "MUTAGEN")
        self.assertEqual(feature.location.start, 155)
        self.assertEqual(feature.location.end, 156)
        self.assertEqual(
            feature.qualifiers["note"],
            "S->C: Impairs the maturation of a peptide-receptive HLA-A*02-B2M complex.",
        )
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000269|PubMed:8805302")
        feature = record.features[129]
        self.assertEqual(feature.type, "MUTAGEN")
        self.assertEqual(feature.location.start, 157)
        self.assertEqual(feature.location.end, 158)
        self.assertEqual(
            feature.qualifiers["note"],
            "T->K: Impairs binding to TAP1-TAP2 transporter, resulting in impaired presentation of intracellular peptides.",
        )
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:8630735,ECO:0000269|PubMed:8805302",
        )
        feature = record.features[130]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 26)
        self.assertEqual(feature.location.end, 36)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:3MRE")
        feature = record.features[131]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 40)
        self.assertEqual(feature.location.end, 43)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:3MRE")
        feature = record.features[132]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 44)
        self.assertEqual(feature.location.end, 52)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:3MRE")
        feature = record.features[133]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 54)
        self.assertEqual(feature.location.end, 61)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:3MRE")
        feature = record.features[134]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 63)
        self.assertEqual(feature.location.end, 66)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:3MRE")
        feature = record.features[135]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 69)
        self.assertEqual(feature.location.end, 73)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:4F7T")
        feature = record.features[136]
        self.assertEqual(feature.type, "HELIX")
        self.assertEqual(feature.location.start, 73)
        self.assertEqual(feature.location.end, 78)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:3MRE")
        feature = record.features[137]
        self.assertEqual(feature.type, "HELIX")
        self.assertEqual(feature.location.start, 80)
        self.assertEqual(feature.location.end, 108)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:3MRE")
        feature = record.features[138]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 112)
        self.assertEqual(feature.location.end, 115)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:3D25")
        feature = record.features[139]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 117)
        self.assertEqual(feature.location.end, 127)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:3MRE")
        feature = record.features[140]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 130)
        self.assertEqual(feature.location.end, 142)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:3MRE")
        feature = record.features[141]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 144)
        self.assertEqual(feature.location.end, 150)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:3MRE")
        feature = record.features[142]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 151)
        self.assertEqual(feature.location.end, 155)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:6EWA")
        feature = record.features[143]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 156)
        self.assertEqual(feature.location.end, 159)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:3MRE")
        feature = record.features[144]
        self.assertEqual(feature.type, "HELIX")
        self.assertEqual(feature.location.start, 161)
        self.assertEqual(feature.location.end, 173)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:3MRE")
        feature = record.features[145]
        self.assertEqual(feature.type, "HELIX")
        self.assertEqual(feature.location.start, 175)
        self.assertEqual(feature.location.end, 184)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:3MRE")
        feature = record.features[146]
        self.assertEqual(feature.type, "HELIX")
        self.assertEqual(feature.location.start, 186)
        self.assertEqual(feature.location.end, 198)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:3MRE")
        feature = record.features[147]
        self.assertEqual(feature.type, "HELIX")
        self.assertEqual(feature.location.start, 199)
        self.assertEqual(feature.location.end, 203)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:3MRE")
        feature = record.features[148]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 209)
        self.assertEqual(feature.location.end, 235)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:3MRE")
        feature = record.features[149]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 237)
        self.assertEqual(feature.location.end, 243)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:3MRE")
        feature = record.features[150]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 245)
        self.assertEqual(feature.location.end, 248)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:3MRE")
        feature = record.features[151]
        self.assertEqual(feature.type, "HELIX")
        self.assertEqual(feature.location.start, 248)
        self.assertEqual(feature.location.end, 251)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:2GTW")
        feature = record.features[152]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 251)
        self.assertEqual(feature.location.end, 254)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:3MRE")
        feature = record.features[153]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 260)
        self.assertEqual(feature.location.end, 263)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:3MRE")
        feature = record.features[154]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 264)
        self.assertEqual(feature.location.end, 274)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:3MRE")
        feature = record.features[155]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 274)
        self.assertEqual(feature.location.end, 277)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:4JFD")
        feature = record.features[156]
        self.assertEqual(feature.type, "HELIX")
        self.assertEqual(feature.location.start, 277)
        self.assertEqual(feature.location.end, 280)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:3MRE")
        feature = record.features[157]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 280)
        self.assertEqual(feature.location.end, 286)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:3MRE")
        feature = record.features[158]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 289)
        self.assertEqual(feature.location.end, 292)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:2V2X")
        feature = record.features[159]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 293)
        self.assertEqual(feature.location.end, 297)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:3MRE")
        feature = record.features[160]
        self.assertEqual(feature.type, "STRAND")
        self.assertEqual(feature.location.start, 347)
        self.assertEqual(feature.location.end, 350)
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0007829|PDB:4EN2")
        self.assertIsNone(feature.id)

        self.assertEqual(len(record.references), 103)
        reference = record.references[0]
        self.assertEqual(reference.authors, "Wan A.M., Ennis P., Parham P., Holmes N.")
        self.assertEqual(
            reference.title,
            "The primary structure of HLA-A32 suggests a region involved in formation of the Bw4/Bw6 epitopes.",
        )
        self.assertEqual(len(reference.references), 1)
        self.assertEqual(reference.references[0], ("PubMed", "2431040"))
        reference = record.references[1]
        self.assertEqual(
            reference.authors, "Holmes N., Ennis P., Wan A.M., Denney D.W., Parham P."
        )
        self.assertEqual(
            reference.title,
            "Multiple genetic mechanisms have contributed to the generation of the HLA-A2/A28 family of class I MHC molecules.",
        )
        self.assertEqual(len(reference.references), 1)
        self.assertEqual(reference.references[0], ("PubMed", "3496393"))
        reference = record.references[2]
        self.assertEqual(
            reference.authors,
            "Mayer W.E., Jonker M., Klein D., Ivanyi P., van Seventer G., Klein J.",
        )
        self.assertEqual(
            reference.title,
            "Nucleotide sequences of chimpanzee MHC class I alleles: evidence for trans-species mode of evolution.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "2460344"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1002/j.1460-2075.1988.tb03131.x")
        )
        reference = record.references[3]
        self.assertEqual(
            reference.authors,
            "Trapani J.A., Mizuno S., Kang S.H., Yang S.Y., Dupont B.",
        )
        self.assertEqual(
            reference.title,
            "Molecular mapping of a new public HLA class I epitope shared by all HLA-B and HLA-C antigens and defined by a monoclonal antibody.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "2461903"))
        self.assertEqual(reference.references[1], ("DOI", "10.1007/bf02341610"))
        reference = record.references[4]
        self.assertEqual(
            reference.authors,
            "Kato K., Trapani J.A., Allopenna J., Dupont B., Yang S.Y.",
        )
        self.assertEqual(
            reference.title,
            "Molecular analysis of the serologically defined HLA-Aw19 antigens. A genetically distinct family of HLA-A antigens comprising A29, A31, A32, and Aw33, but probably not A30.",
        )
        self.assertEqual(len(reference.references), 1)
        self.assertEqual(reference.references[0], ("PubMed", "2478623"))
        reference = record.references[5]
        self.assertEqual(
            reference.authors, "Parham P., Lawlor D.A., Lomen C.E., Ennis P.D."
        )
        self.assertEqual(
            reference.title, "Diversity and diversification of HLA-A,B,C alleles."
        )
        self.assertEqual(len(reference.references), 1)
        self.assertEqual(reference.references[0], ("PubMed", "2715640"))
        reference = record.references[6]
        self.assertEqual(
            reference.authors, "Ennis P.D., Zemmour J., Salter R.D., Parham P."
        )
        self.assertEqual(
            reference.title,
            "Rapid cloning of HLA-A,B cDNA by using the polymerase chain reaction: frequency and nature of errors produced in amplification.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "2320591"))
        self.assertEqual(reference.references[1], ("DOI", "10.1073/pnas.87.7.2833"))
        reference = record.references[7]
        self.assertEqual(
            reference.authors,
            "Tabary T., Prochnicka-Chalufour A., Cornillet P., Lehoang P., Betuel H., Cohen H.M.",
        )
        self.assertEqual(
            reference.title,
            "HLA-A29 sub-types and 'Birdshot' choroido-retinopathy susceptibility: a possible 'resistance motif' in the HLA-A29.1 molecule.",
        )
        self.assertEqual(len(reference.references), 1)
        self.assertEqual(reference.references[0], ("PubMed", "1782566"))
        reference = record.references[8]
        self.assertEqual(reference.authors, "Little A.-M., Madrigal J.A., Parham P.")
        self.assertEqual(
            reference.title,
            "Molecular definition of an elusive third HLA-A9 molecule: HLA-A9.3.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "1729171"))
        self.assertEqual(reference.references[1], ("DOI", "10.1007/bf00216625"))
        reference = record.references[9]
        self.assertEqual(
            reference.authors,
            "Madrigal J.A., Belich M.P., Hildebrand W.H., Benjamin R.J., Little A.-M., Zemmour J., Ennis P.D., Ward F.E., Petzl-Erler M.L., Martell R.W., du Toit E.D., Parham P.",
        )
        self.assertEqual(
            reference.title,
            "Distinctive HLA-A,B antigens of black populations formed by interallelic conversion.",
        )
        self.assertEqual(len(reference.references), 1)
        self.assertEqual(reference.references[0], ("PubMed", "1431115"))
        reference = record.references[10]
        self.assertEqual(
            reference.authors,
            "Belich M.P., Madrigal J.A., Hildebrand W.H., Zemmour J., Williams R.C., Luz R., Petzl-Erler M.L., Parham P.",
        )
        self.assertEqual(
            reference.title, "Unusual HLA-B alleles in two tribes of Brazilian Indians."
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "1317015"))
        self.assertEqual(reference.references[1], ("DOI", "10.1038/357326a0"))
        reference = record.references[11]
        self.assertEqual(
            reference.authors,
            "Madrigal J.A., Hildebrand W.H., Belich M.P., Benjamin R.J., Little A.-M., Zemmour J., Ennis P.D., Ward F.E., Petzl-Erler M.L., du Toit E.D., Parham P.",
        )
        self.assertEqual(
            reference.title,
            "Structural diversity in the HLA-A10 family of alleles: correlations with serology.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "8475492"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1111/j.1399-0039.1993.tb01982.x")
        )
        reference = record.references[12]
        self.assertEqual(
            reference.authors, "Domena J.D., Hildebrand W.H., Bias W.B., Parham P."
        )
        self.assertEqual(
            reference.title, "A sixth family of HLA-A alleles defined by HLA-A*8001."
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "8284791"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1111/j.1399-0039.1993.tb02186.x")
        )
        reference = record.references[13]
        self.assertEqual(
            reference.authors,
            "Ishikawa Y., Tokunaga K., Lin L., Imanishi T., Saitou S., Kashiwase K., Akaza T., Tadokoro K., Juji T.",
        )
        self.assertEqual(
            reference.title,
            "Sequences of four splits of HLA-A10 group. Implications for serologic cross-reactivities and their evolution.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "8026990"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1016/0198-8859(94)90263-1")
        )
        reference = record.references[14]
        self.assertEqual(
            reference.authors,
            "Balas A., Garcia-Sanchez F., Gomez-Reino F., Vicario J.L.",
        )
        self.assertEqual(
            reference.title,
            "Characterization of a new and highly distinguishable HLA-A allele in a Spanish family.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "8188325"))
        self.assertEqual(reference.references[1], ("DOI", "10.1007/bf00176169"))
        reference = record.references[15]
        self.assertEqual(
            reference.authors,
            "Lin L., Tokunaga K., Ishikawa Y., Bannai M., Kashiwase K., Kuwata S., Akaza T., Tadokoro K., Shibata Y., Juji T.",
        )
        self.assertEqual(
            reference.title,
            "Sequence analysis of serological HLA-A11 split antigens, A11.1 and A11.2.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "8016845"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1111/j.1399-0039.1994.tb02304.x")
        )
        reference = record.references[16]
        self.assertEqual(reference.authors, "Olerup O., Daniels T.J., Baxter-Lowe L.")
        self.assertEqual(
            reference.title,
            "Correct sequence of the A*3001 allele obtained by PCR-SSP typing and automated nucleotide sequencing.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "7871528"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1111/j.1399-0039.1994.tb02393.x")
        )
        reference = record.references[17]
        self.assertEqual(reference.authors, "Sun Y., Liu S., Luo Y., Liang F., Xi Y.")
        self.assertEqual(
            reference.title,
            "Identification and frequency of a novel HLA-A allele, A*110104.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "17092262"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1111/j.1399-0039.2006.00687.x")
        )
        reference = record.references[18]
        self.assertEqual(
            reference.authors, "Strachan T., Sodoyer R., Damotte M., Jordan B.R."
        )
        self.assertEqual(
            reference.title,
            "Complete nucleotide sequence of a functional class I HLA gene, HLA-A3: implications for the evolution of HLA genes.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "6609814"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1002/j.1460-2075.1984.tb01901.x")
        )
        reference = record.references[19]
        self.assertEqual(reference.authors, "Holmes N., Parham P.")
        self.assertEqual(
            reference.title,
            "Exon shuffling in vivo can generate novel HLA class I molecules.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "3877632"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1002/j.1460-2075.1985.tb04013.x")
        )
        reference = record.references[20]
        self.assertEqual(reference.authors, "Koller B.H., Orr H.T.")
        self.assertEqual(
            reference.title,
            "Cloning and complete sequence of an HLA-A2 gene: analysis of two HLA-A alleles at the nucleotide level.",
        )
        self.assertEqual(len(reference.references), 1)
        self.assertEqual(reference.references[0], ("PubMed", "2982951"))
        reference = record.references[21]
        self.assertEqual(
            reference.authors, "Cowan E.P., Jelachich M.L., Biddison W.E., Coligan J.E."
        )
        self.assertEqual(
            reference.title,
            "DNA sequence of HLA-A11: remarkable homology with HLA-A3 allows identification of residues involved in epitopes recognized by antibodies and T cells.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "2437024"))
        self.assertEqual(reference.references[1], ("DOI", "10.1007/bf00404694"))
        reference = record.references[22]
        self.assertEqual(reference.authors, "Girdlestone J.")
        self.assertEqual(reference.title, "Nucleotide sequence of an HLA-A1 gene.")
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "2251137"))
        self.assertEqual(reference.references[1], ("DOI", "10.1093/nar/18.22.6701"))
        reference = record.references[23]
        self.assertEqual(
            reference.authors,
            "Balas A., Garcia-Sanchez F., Gomez-Reino F., Vicario J.L.",
        )
        self.assertEqual(
            reference.title,
            "HLA class I allele (HLA-A2) expression defect associated with a mutation in its enhancer B inverted CAT box in two families.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "7836067"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1016/0198-8859(94)90087-6")
        )
        reference = record.references[24]
        self.assertEqual(reference.authors, "Arnett K.L., Adams E.J., Parham P.")
        self.assertEqual(reference.title, "On the sequence of A*3101.")
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "8795145"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1111/j.1399-0039.1996.tb02580.x")
        )
        reference = record.references[25]
        self.assertEqual(
            reference.authors,
            "Laforet M., Froelich N., Parissiadis A., Bausinger H., Pfeiffer B., Tongio M.M.",
        )
        self.assertEqual(
            reference.title,
            "An intronic mutation responsible for a low level of expression of an HLA-A*24 allele.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "9349616"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1111/j.1399-0039.1997.tb02884.x")
        )
        reference = record.references[26]
        self.assertEqual(
            reference.authors,
            "Laforet M., Froelich N., Parissiadis A., Pfeiffer B., Schell A., Faller B., Woehl-Jaegle M.L., Cazenave J.-P., Tongio M.M.",
        )
        self.assertEqual(
            reference.title,
            "A nucleotide insertion in exon 4 is responsible for the absence of expression of an HLA-A*01 allele.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "9349617"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1111/j.1399-0039.1997.tb02885.x")
        )
        reference = record.references[27]
        self.assertEqual(
            reference.authors, "Zhu F., He Y., Zhang W., He J., He J., Xu X., Yan L."
        )
        self.assertEqual(
            reference.title,
            "Analysis of the complete genomic sequence of HLA-A alleles in the Chinese Han population.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "19735485"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1111/j.1744-313x.2009.00874.x")
        )
        reference = record.references[28]
        self.assertEqual(reference.authors, "Domena J.D.")
        self.assertEqual(reference.title, "")
        self.assertEqual(len(reference.references), 0)
        reference = record.references[29]
        self.assertEqual(reference.authors, "Hurley C.K.")
        self.assertEqual(reference.title, "")
        self.assertEqual(len(reference.references), 0)
        reference = record.references[30]
        self.assertEqual(reference.authors, "Ellexson M.E., Hildebrand W.H.")
        self.assertEqual(reference.title, "")
        self.assertEqual(len(reference.references), 0)
        reference = record.references[31]
        self.assertEqual(reference.authors, "Mayor N.P.")
        self.assertEqual(
            reference.title, "Full length sequence of an HLA-A*0301 intron 2 variant."
        )
        self.assertEqual(len(reference.references), 0)
        reference = record.references[32]
        self.assertEqual(reference.authors, "The MGC Project Team")
        self.assertEqual(
            reference.title,
            "The status, quality, and expansion of the NIH full-length cDNA project: the Mammalian Gene Collection (MGC).",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "15489334"))
        self.assertEqual(reference.references[1], ("DOI", "10.1101/gr.2596504"))
        reference = record.references[33]
        self.assertEqual(
            reference.authors,
            "Mungall A.J., Palmer S.A., Sims S.K., Edwards C.A., Ashurst J.L., Wilming L., Jones M.C., Horton R., Hunt S.E., Scott C.E., Gilbert J.G.R., Clamp M.E., Bethel G., Milne S., Ainscough R., Almeida J.P., Ambrose K.D., Andrews T.D., Ashwell R.I.S., Babbage A.K., Bagguley C.L., Bailey J., Banerjee R., Barker D.J., Barlow K.F., Bates K., Beare D.M., Beasley H., Beasley O., Bird C.P., Blakey S.E., Bray-Allen S., Brook J., Brown A.J., Brown J.Y., Burford D.C., Burrill W., Burton J., Carder C., Carter N.P., Chapman J.C., Clark S.Y., Clark G., Clee C.M., Clegg S., Cobley V., Collier R.E., Collins J.E., Colman L.K., Corby N.R., Coville G.J., Culley K.M., Dhami P., Davies J., Dunn M., Earthrowl M.E., Ellington A.E., Evans K.A., Faulkner L., Francis M.D., Frankish A., Frankland J., French L., Garner P., Garnett J., Ghori M.J., Gilby L.M., Gillson C.J., Glithero R.J., Grafham D.V., Grant M., Gribble S., Griffiths C., Griffiths M.N.D., Hall R., Halls K.S., Hammond S., Harley J.L., Hart E.A., Heath P.D., Heathcott R., Holmes S.J., Howden P.J., Howe K.L., Howell G.R., Huckle E., Humphray S.J., Humphries M.D., Hunt A.R., Johnson C.M., Joy A.A., Kay M., Keenan S.J., Kimberley A.M., King A., Laird G.K., Langford C., Lawlor S., Leongamornlert D.A., Leversha M., Lloyd C.R., Lloyd D.M., Loveland J.E., Lovell J., Martin S., Mashreghi-Mohammadi M., Maslen G.L., Matthews L., McCann O.T., McLaren S.J., McLay K., McMurray A., Moore M.J.F., Mullikin J.C., Niblett D., Nickerson T., Novik K.L., Oliver K., Overton-Larty E.K., Parker A., Patel R., Pearce A.V., Peck A.I., Phillimore B.J.C.T., Phillips S., Plumb R.W., Porter K.M., Ramsey Y., Ranby S.A., Rice C.M., Ross M.T., Searle S.M., Sehra H.K., Sheridan E., Skuce C.D., Smith S., Smith M., Spraggon L., Squares S.L., Steward C.A., Sycamore N., Tamlyn-Hall G., Tester J., Theaker A.J., Thomas D.W., Thorpe A., Tracey A., Tromans A., Tubby B., Wall M., Wallis J.M., West A.P., White S.S., Whitehead S.L., Whittaker H., Wild A., Willey D.J., Wilmer T.E., Wood J.M., Wray P.W., Wyatt J.C., Young L., Younger R.M., Bentley D.R., Coulson A., Durbin R.M., Hubbard T., Sulston J.E., Dunham I., Rogers J., Beck S.",
        )
        self.assertEqual(
            reference.title, "The DNA sequence and analysis of human chromosome 6."
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "14574404"))
        self.assertEqual(reference.references[1], ("DOI", "10.1038/nature02055"))
        reference = record.references[34]
        self.assertEqual(
            reference.authors,
            "Orr H.T., Lopez de Castro J.A., Parham P., Ploegh H.L., Strominger J.L.",
        )
        self.assertEqual(
            reference.title,
            "Comparison of amino acid sequences of two human histocompatibility antigens, HLA-A2 and HLA-B7: location of putative alloantigenic sites.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "92029"))
        self.assertEqual(reference.references[1], ("DOI", "10.1073/pnas.76.9.4395"))
        reference = record.references[35]
        self.assertEqual(
            reference.authors,
            "Lopez de Castro J.A., Strominger J.L., Strong D.M., Orr H.T.",
        )
        self.assertEqual(
            reference.title,
            "Structure of crossreactive human histocompatibility antigens HLA-A28 and HLA-A2: possible implications for the generation of HLA polymorphism.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "6179086"))
        self.assertEqual(reference.references[1], ("DOI", "10.1073/pnas.79.12.3813"))
        reference = record.references[36]
        self.assertEqual(
            reference.authors,
            "Jelachich M.L., Cowan E.P., Turner R.V., Coligan J.E., Biddison W.E.",
        )
        self.assertEqual(
            reference.title,
            "Analysis of the molecular basis of HLA-A3 recognition by cytotoxic T cells using defined mutants of the HLA-A3 molecule.",
        )
        self.assertEqual(len(reference.references), 1)
        self.assertEqual(reference.references[0], ("PubMed", "2456340"))
        reference = record.references[37]
        self.assertEqual(
            reference.authors,
            "Salter R.D., Norment A.M., Chen B.P., Clayberger C., Krensky A.M., Littman D.R., Parham P.",
        )
        self.assertEqual(
            reference.title,
            "Polymorphism in the alpha 3 domain of HLA-A molecules affects binding to CD8.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "2784196"))
        self.assertEqual(reference.references[1], ("DOI", "10.1038/338345a0"))
        reference = record.references[38]
        self.assertEqual(
            reference.authors,
            "Traversari C., van der Bruggen P., Luescher I.F., Lurquin C., Chomez P., Van Pel A., De Plaen E., Amar-Costesec A., Boon T.",
        )
        self.assertEqual(
            reference.title,
            "A nonapeptide encoded by human gene MAGE-1 is recognized on HLA-A1 by cytolytic T lymphocytes directed against tumor antigen MZ2-E.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "1402688"))
        self.assertEqual(reference.references[1], ("DOI", "10.1084/jem.176.5.1453"))
        reference = record.references[39]
        self.assertEqual(
            reference.authors,
            "DiBrino M., Tsuchida T., Turner R.V., Parker K.C., Coligan J.E., Biddison W.E.",
        )
        self.assertEqual(
            reference.title,
            "HLA-A1 and HLA-A3 T cell epitopes derived from influenza virus proteins predicted from peptide binding motifs.",
        )
        self.assertEqual(len(reference.references), 1)
        self.assertEqual(reference.references[0], ("PubMed", "7504010"))
        reference = record.references[40]
        self.assertEqual(
            reference.authors,
            "DiBrino M., Parker K.C., Shiloach J., Knierman M., Lukszo J., Turner R.V., Biddison W.E., Coligan J.E.",
        )
        self.assertEqual(
            reference.title,
            "Endogenous peptides bound to HLA-A3 possess a specific combination of anchor residues that permit identification of potential antigenic peptides.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "7679507"))
        self.assertEqual(reference.references[1], ("DOI", "10.1073/pnas.90.4.1508"))
        reference = record.references[41]
        self.assertEqual(
            reference.authors,
            "DiBrino M., Parker K.C., Shiloach J., Turner R.V., Tsuchida T., Garfield M., Biddison W.E., Coligan J.E.",
        )
        self.assertEqual(
            reference.title,
            "Endogenous peptides with distinct amino acid anchor residue motifs bind to HLA-A1 and HLA-B8.",
        )
        self.assertEqual(len(reference.references), 1)
        self.assertEqual(reference.references[0], ("PubMed", "7506728"))
        reference = record.references[42]
        self.assertEqual(
            reference.authors, "Lewis J.W., Neisig A., Neefjes J., Elliott T."
        )
        self.assertEqual(
            reference.title,
            "Point mutations in the alpha 2 domain of HLA-A2.1 define a functionally relevant interaction with TAP.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "8805302"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1016/s0960-9822(02)00611-5")
        )
        reference = record.references[43]
        self.assertEqual(
            reference.authors,
            "Peace-Brewer A.L., Tussey L.G., Matsui M., Li G., Quinn D.G., Frelinger J.A.",
        )
        self.assertEqual(
            reference.title,
            "A point mutation in HLA-A*0201 results in failure to bind the TAP complex and to present virus-derived peptides to CTL.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "8630735"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1016/s1074-7613(00)80416-1")
        )
        reference = record.references[44]
        self.assertEqual(
            reference.authors,
            "Boisgerault F., Khalil I., Tieng V., Connan F., Tabary T., Cohen J.H., Choppin J., Charron D., Toubert A.",
        )
        self.assertEqual(
            reference.title,
            "Definition of the HLA-A29 peptide ligand motif allows prediction of potential T-cell epitopes from the retinal soluble antigen, a candidate autoantigen in birdshot retinopathy.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "8622959"))
        self.assertEqual(reference.references[1], ("DOI", "10.1073/pnas.93.8.3466"))
        reference = record.references[45]
        self.assertEqual(
            reference.authors,
            "Ikeda H., Lethe B.G., Lehmann F., van Baren N., Baurain J.-F., de Smet C., Chambost H., Vitale M., Moretta A., Boon T., Coulie P.G.",
        )
        self.assertEqual(
            reference.title,
            "Characterization of an antigen that is recognized on a melanoma showing partial HLA loss by CTL expressing an NK inhibitory receptor.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "9047241"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1016/s1074-7613(00)80426-4")
        )
        reference = record.references[46]
        self.assertEqual(
            reference.authors,
            "Kawakami Y., Robbins P.F., Wang X., Tupesis J.P., Parkhurst M.R., Kang X., Sakaguchi K., Appella E., Rosenberg S.A.",
        )
        self.assertEqual(
            reference.title,
            "Identification of new melanoma epitopes on melanosomal proteins recognized by tumor infiltrating T lymphocytes restricted by HLA-A1, -A2, and -A3 alleles.",
        )
        self.assertEqual(len(reference.references), 1)
        self.assertEqual(reference.references[0], ("PubMed", "9862734"))
        reference = record.references[47]
        self.assertEqual(
            reference.authors,
            "Fukada K., Chujoh Y., Tomiyama H., Miwa K., Kaneko Y., Oka S., Takiguchi M.",
        )
        self.assertEqual(
            reference.title,
            "HLA-A*1101-restricted cytotoxic T lymphocyte recognition of HIV-1 Pol protein.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "10449296"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1097/00002030-199907300-00021")
        )
        reference = record.references[48]
        self.assertEqual(
            reference.authors,
            "Johnson J.M., Nicot C., Fullen J., Ciminale V., Casareto L., Mulloy J.C., Jacobson S., Franchini G.",
        )
        self.assertEqual(
            reference.title,
            "Free major histocompatibility complex class I heavy chain is preferentially targeted for degradation by human T-cell leukemia/lymphotropic virus type 1 p12(I) protein.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "11390610"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1128/jvi.75.13.6086-6094.2001")
        )
        reference = record.references[49]
        self.assertEqual(
            reference.authors,
            "Hewitt E.W., Duncan L., Mufti D., Baker J., Stevenson P.G., Lehner P.J.",
        )
        self.assertEqual(
            reference.title,
            "Ubiquitylation of MHC class I by the K3 viral protein signals internalization and TSG101-dependent degradation.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "12006494"))
        self.assertEqual(reference.references[1], ("DOI", "10.1093/emboj/21.10.2418"))
        reference = record.references[50]
        self.assertEqual(
            reference.authors,
            "Nagata Y., Ono S., Matsuo M., Gnjatic S., Valmori D., Ritter G., Garrett W., Old L.J., Mellman I.",
        )
        self.assertEqual(
            reference.title,
            "Differential presentation of a soluble exogenous tumor antigen, NY-ESO-1, by distinct human dendritic cell populations.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "12138174"))
        self.assertEqual(reference.references[1], ("DOI", "10.1073/pnas.112331099"))
        reference = record.references[51]
        self.assertEqual(
            reference.authors,
            "Kuzushima K., Hayashi N., Kudoh A., Akatsuka Y., Tsujimura K., Morishima Y., Tsurumi T.",
        )
        self.assertEqual(
            reference.title,
            "Tetramer-assisted identification and characterization of epitopes recognized by HLA A*2402-restricted Epstein-Barr virus-specific CD8+ T cells.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "12393434"))
        self.assertEqual(reference.references[1], ("DOI", "10.1182/blood-2002-04-1240"))
        reference = record.references[52]
        self.assertEqual(
            reference.authors,
            "Satoh M., Takamiya Y., Oka S., Tokunaga K., Takiguchi M.",
        )
        self.assertEqual(
            reference.title,
            "Identification and characterization of HIV-1-specific CD8+ T cell epitopes presented by HLA-A*2601.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "15893615"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1016/j.vaccine.2005.02.022")
        )
        reference = record.references[53]
        self.assertEqual(
            reference.authors,
            "Asemissen A.M., Keilholz U., Tenzer S., Mueller M., Walter S., Stevanovic S., Schild H., Letsch A., Thiel E., Rammensee H.G., Scheibenbogen C.",
        )
        self.assertEqual(
            reference.title,
            "Identification of a highly immunogenic HLA-A*01-binding T cell epitope of WT1.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "17189421"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1158/1078-0432.ccr-06-1337")
        )
        reference = record.references[54]
        self.assertEqual(
            reference.authors,
            "Thananchai H., Gillespie G., Martin M.P., Bashirova A., Yawata N., Yawata M., Easterbrook P., McVicar D.W., Maenaka K., Parham P., Carrington M., Dong T., Rowland-Jones S.",
        )
        self.assertEqual(
            reference.title,
            "Cutting Edge: Allele-specific and peptide-dependent interactions between KIR3DL1 and HLA-A and HLA-B.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "17182537"))
        self.assertEqual(reference.references[1], ("DOI", "10.4049/jimmunol.178.1.33"))
        reference = record.references[55]
        self.assertEqual(
            reference.authors, "Robek M.D., Garcia M.L., Boyd B.S., Chisari F.V."
        )
        self.assertEqual(
            reference.title,
            "Role of immunoproteasome catalytic subunits in the immune response to hepatitis B virus.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "17079320"))
        self.assertEqual(reference.references[1], ("DOI", "10.1128/jvi.01779-06"))
        reference = record.references[56]
        self.assertEqual(
            reference.authors,
            "Stern M., Ruggeri L., Capanni M., Mancusi A., Velardi A.",
        )
        self.assertEqual(
            reference.title,
            "Human leukocyte antigens A23, A24, and A32 but not A25 are ligands for KIR3DL1.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "18502829"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1182/blood-2008-02-137521")
        )
        reference = record.references[57]
        self.assertEqual(reference.authors, "Brennan R.M., Burrows S.R.")
        self.assertEqual(
            reference.title,
            "A mechanism for the HLA-A*01-associated risk for EBV+ Hodgkin lymphoma and infectious mononucleosis.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "18779413"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1182/blood-2008-06-162883")
        )
        reference = record.references[58]
        self.assertEqual(
            reference.authors,
            "Chen R., Jiang X., Sun D., Han G., Wang F., Ye M., Wang L., Zou H.",
        )
        self.assertEqual(
            reference.title,
            "Glycoproteomics analysis of human liver tissue by combination of multiple enzyme digestion and hydrazide chemistry.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "19159218"))
        self.assertEqual(reference.references[1], ("DOI", "10.1021/pr8008012"))
        reference = record.references[59]
        self.assertEqual(
            reference.authors,
            "Hadrup S.R., Bakker A.H., Shu C.J., Andersen R.S., van Veluw J., Hombrink P., Castermans E., Thor Straten P., Blank C., Haanen J.B., Heemskerk M.H., Schumacher T.N.",
        )
        self.assertEqual(
            reference.title,
            "Parallel detection of antigen-specific T-cell responses by multidimensional encoding of MHC multimers.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "19543285"))
        self.assertEqual(reference.references[1], ("DOI", "10.1038/nmeth.1345"))
        reference = record.references[60]
        self.assertEqual(
            reference.authors,
            "Parmentier N., Stroobant V., Colau D., de Diesbach P., Morel S., Chapiro J., van Endert P., Van den Eynde B.J.",
        )
        self.assertEqual(
            reference.title,
            "Production of an antigenic peptide by insulin-degrading enzyme.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "20364150"))
        self.assertEqual(reference.references[1], ("DOI", "10.1038/ni.1862"))
        reference = record.references[61]
        self.assertEqual(
            reference.authors,
            "Marsh S.G., Albert E.D., Bodmer W.F., Bontrop R.E., Dupont B., Erlich H.A., Fernandez-Vina M., Geraghty D.E., Holdsworth R., Hurley C.K., Lau M., Lee K.W., Mach B., Maiers M., Mayr W.R., Mueller C.R., Parham P., Petersdorf E.W., Sasazuki T., Strominger J.L., Svejgaard A., Terasaki P.I., Tiercy J.M., Trowsdale J.",
        )
        self.assertEqual(
            reference.title, "Nomenclature for factors of the HLA system, 2010."
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "20356336"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1111/j.1399-0039.2010.01466.x")
        )
        reference = record.references[62]
        self.assertEqual(
            reference.authors, "Rizvi S.M., Del Cid N., Lybarger L., Raghavan M."
        )
        self.assertEqual(
            reference.title,
            "Distinct functions for the glycans of tapasin and heavy chains in the assembly of MHC class I molecules.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "21263072"))
        self.assertEqual(reference.references[1], ("DOI", "10.4049/jimmunol.1002959"))
        reference = record.references[63]
        self.assertEqual(
            reference.authors,
            "Matthews P.C., Adland E., Listgarten J., Leslie A., Mkhwanazi N., Carlson J.M., Harndahl M., Stryhn A., Payne R.P., Ogwu A., Huang K.H., Frater J., Paioni P., Kloverpris H., Jooste P., Goedhals D., van Vuuren C., Steyn D., Riddell L., Chen F., Luzzi G., Balachandran T., Ndung'u T., Buus S., Carrington M., Shapiro R., Heckerman D., Goulder P.J.",
        )
        self.assertEqual(
            reference.title,
            "HLA-A*7401-mediated control of HIV viremia is independent of its linkage disequilibrium with HLA-B*5703.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "21498667"))
        self.assertEqual(reference.references[1], ("DOI", "10.4049/jimmunol.1003711"))
        reference = record.references[64]
        self.assertEqual(
            reference.authors,
            "Zhou H., Di Palma S., Preisinger C., Peng M., Polat A.N., Heck A.J., Mohammed S.",
        )
        self.assertEqual(
            reference.title,
            "Toward a comprehensive characterization of a human cancer cell phosphoproteome.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "23186163"))
        self.assertEqual(reference.references[1], ("DOI", "10.1021/pr300630k"))
        reference = record.references[65]
        self.assertEqual(
            reference.authors,
            "Shimizu A., Kawana-Tachikawa A., Yamagata A., Han C., Zhu D., Sato Y., Nakamura H., Koibuchi T., Carlson J., Martin E., Brumme C.J., Shi Y., Gao G.F., Brumme Z.L., Fukai S., Iwamoto A.",
        )
        self.assertEqual(
            reference.title,
            "Structure of TCR and antigen complexes at an immunodominant CTL epitope in HIV-1 infection.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "24192765"))
        self.assertEqual(reference.references[1], ("DOI", "10.1038/srep03097"))
        reference = record.references[66]
        self.assertEqual(
            reference.authors,
            "Bian Y., Song C., Cheng K., Dong M., Wang F., Huang J., Sun D., Wang L., Ye M., Zou H.",
        )
        self.assertEqual(
            reference.title,
            "An enzyme assisted RP-RPLC approach for in-depth analysis of human liver phosphoproteome.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "24275569"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1016/j.jprot.2013.11.014")
        )
        reference = record.references[67]
        self.assertEqual(
            reference.authors,
            "Vaca Jacome A.S., Rabilloud T., Schaeffer-Reiss C., Rompais M., Ayoub D., Lane L., Bairoch A., Van Dorsselaer A., Carapito C.",
        )
        self.assertEqual(
            reference.title, "N-terminome analysis of the human mitochondrial proteome."
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "25944712"))
        self.assertEqual(reference.references[1], ("DOI", "10.1002/pmic.201400617"))
        reference = record.references[68]
        self.assertEqual(
            reference.authors,
            "Giam K., Ayala-Perez R., Illing P.T., Schittenhelm R.B., Croft N.P., Purcell A.W., Dudek N.L.",
        )
        self.assertEqual(
            reference.title, "A comprehensive analysis of peptides presented by HLA-A1."
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "25880248"))
        self.assertEqual(reference.references[1], ("DOI", "10.1111/tan.12565"))
        reference = record.references[69]
        self.assertEqual(
            reference.authors,
            "Morozov G.I., Zhao H., Mage M.G., Boyd L.F., Jiang J., Dolan M.A., Venna R., Norcross M.A., McMurtrey C.P., Hildebrand W., Schuck P., Natarajan K., Margulies D.H.",
        )
        self.assertEqual(
            reference.title,
            "Interaction of TAPBPR, a tapasin homolog, with MHC-I molecules promotes peptide editing.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "26869717"))
        self.assertEqual(reference.references[1], ("DOI", "10.1073/pnas.1519894113"))
        reference = record.references[70]
        self.assertEqual(
            reference.authors,
            "Tripathi S.C., Peters H.L., Taguchi A., Katayama H., Wang H., Momin A., Jolly M.K., Celiktas M., Rodriguez-Canales J., Liu H., Behrens C., Wistuba I.I., Ben-Jacob E., Levine H., Molldrem J.J., Hanash S.M., Ostrin E.J.",
        )
        self.assertEqual(
            reference.title,
            "Immunoproteasome deficiency is a feature of non-small cell lung cancer with a mesenchymal phenotype and is associated with a poor outcome.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "26929325"))
        self.assertEqual(reference.references[1], ("DOI", "10.1073/pnas.1521812113"))
        reference = record.references[71]
        self.assertEqual(
            reference.authors,
            "Ebstein F., Textoris-Taube K., Keller C., Golnik R., Vigneron N., Van den Eynde B.J., Schuler-Thurner B., Schadendorf D., Lorenz F.K., Uckert W., Urban S., Lehmann A., Albrecht-Koepke N., Janek K., Henklein P., Niewienda A., Kloetzel P.M., Mishto M.",
        )
        self.assertEqual(
            reference.title,
            "Proteasomes generate spliced epitopes by two different mechanisms and as efficiently as non-spliced epitopes.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "27049119"))
        self.assertEqual(reference.references[1], ("DOI", "10.1038/srep24032"))
        reference = record.references[72]
        self.assertEqual(
            reference.authors,
            "Keib A., Mei Y.F., Cicin-Sain L., Busch D.H., Dennehy K.M.",
        )
        self.assertEqual(
            reference.title,
            "Measuring Antiviral Capacity of T Cell Responses to Adenovirus.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "30530481"))
        self.assertEqual(reference.references[1], ("DOI", "10.4049/jimmunol.1801003"))
        reference = record.references[73]
        self.assertEqual(
            reference.authors,
            "Oxford Immunology Network Covid-19 Response T cell Consortium; ISARIC4C Investigators; Peng Y., Mentzer A.J., Liu G., Yao X., Yin Z., Dong D., Dejnirattisai W., Rostron T., Supasa P., Liu C., Lopez-Camacho C., Slon-Campos J., Zhao Y., Stuart D.I., Paesen G.C., Grimes J.M., Antson A.A., Bayfield O.W., Hawkins D.E.D.P., Ker D.S., Wang B., Turtle L., Subramaniam K., Thomson P., Zhang P., Dold C., Ratcliff J., Simmonds P., de Silva T., Sopp P., Wellington D., Rajapaksa U., Chen Y.L., Salio M., Napolitani G., Paes W., Borrow P., Kessler B.M., Fry J.W., Schwabe N.F., Semple M.G., Baillie J.K., Moore S.C., Openshaw P.J.M., Ansari M.A., Dunachie S., Barnes E., Frater J., Kerr G., Goulder P., Lockett T., Levin R., Zhang Y., Jing R., Ho L.P., Cornall R.J., Conlon C.P., Klenerman P., Screaton G.R., Mongkolsapaya J., McMichael A., Knight J.C., Ogg G., Dong T.",
        )
        self.assertEqual(
            reference.title,
            "Broad and strong memory CD4+ and CD8+ T cells induced by SARS-CoV-2 in UK convalescent individuals following COVID-19.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "32887977"))
        self.assertEqual(reference.references[1], ("DOI", "10.1038/s41590-020-0782-6"))
        reference = record.references[74]
        self.assertEqual(
            reference.authors,
            "Guo H.-C., Jardetzky T.S., Garrett T.P.J., Lane W.S., Strominger J.L., Wiley D.C.",
        )
        self.assertEqual(
            reference.title,
            "Different length peptides bind to HLA-Aw68 similarly at their ends but bulge out in the middle.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "1448153"))
        self.assertEqual(reference.references[1], ("DOI", "10.1038/360364a0"))
        reference = record.references[75]
        self.assertEqual(
            reference.authors, "Silver M.L., Guo H.-C., Strominger J.L., Wiley D.C."
        )
        self.assertEqual(
            reference.title,
            "Atomic structure of a human MHC molecule presenting an influenza virus peptide.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "1448154"))
        self.assertEqual(reference.references[1], ("DOI", "10.1038/360367a0"))
        reference = record.references[76]
        self.assertEqual(reference.authors, "Madden D.R., Garboczi D.N., Wiley D.C.")
        self.assertEqual(
            reference.title,
            "The antigenic identity of peptide-MHC complexes: a comparison of the conformations of five viral peptides presented by HLA-A2.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "7694806"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1016/0092-8674(93)90490-h")
        )
        reference = record.references[77]
        self.assertEqual(reference.authors, "Collins E.J., Garboczi D.N., Wiley D.C.")
        self.assertEqual(
            reference.title,
            "Three-dimensional structure of a peptide extending from one end of a class I MHC binding site.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "7935798"))
        self.assertEqual(reference.references[1], ("DOI", "10.1038/371626a0"))
        reference = record.references[78]
        self.assertEqual(
            reference.authors,
            "Garboczi D.N., Ghosh P., Utz U., Fan Q.R., Biddison W.E., Wiley D.C.",
        )
        self.assertEqual(
            reference.title,
            "Structure of the complex between human T-cell receptor, viral peptide and HLA-A2.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "8906788"))
        self.assertEqual(reference.references[1], ("DOI", "10.1038/384134a0"))
        reference = record.references[79]
        self.assertEqual(
            reference.authors,
            "Gao G.F., Tormo J., Gerth U.C., Wyer J.R., McMichael A.J., Stuart D.I., Bell J.I., Jones E.Y., Jakobsen B.K.",
        )
        self.assertEqual(
            reference.title,
            "Crystal structure of the complex between human CD8alpha(alpha) and HLA-A2.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "9177355"))
        self.assertEqual(reference.references[1], ("DOI", "10.1038/42523"))
        reference = record.references[80]
        self.assertEqual(
            reference.authors,
            "Hillig R.C., Coulie P.G., Stroobant V., Saenger W., Ziegler A., Hulsmeyer M.",
        )
        self.assertEqual(
            reference.title,
            "High-resolution structure of HLA-A*0201 in complex with a tumour-specific antigenic peptide encoded by the MAGE-A4 gene.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "11502003"))
        self.assertEqual(reference.references[1], ("DOI", "10.1006/jmbi.2001.4816"))
        reference = record.references[81]
        self.assertEqual(
            reference.authors,
            "Stewart-Jones G.B.E., McMichael A.J., Bell J.I., Stuart D.I., Jones E.Y.",
        )
        self.assertEqual(
            reference.title,
            "A structural basis for immunodominant human T cell receptor recognition.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "12796775"))
        self.assertEqual(reference.references[1], ("DOI", "10.1038/ni942"))
        reference = record.references[82]
        self.assertEqual(
            reference.authors, "Blicher T., Kastrup J.S., Buus S., Gajhede M."
        )
        self.assertEqual(
            reference.title,
            "High-resolution structure of HLA-A*1101 in complex with SARS nucleocapsid peptide.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "16041067"))
        self.assertEqual(reference.references[1], ("DOI", "10.1107/s0907444905013090"))
        reference = record.references[83]
        self.assertEqual(
            reference.authors,
            "Ishizuka J., Stewart-Jones G.B., van der Merwe A., Bell J.I., McMichael A.J., Jones E.Y.",
        )
        self.assertEqual(
            reference.title,
            "The structural dynamics and energetics of an immunodominant T cell receptor are programmed by its Vbeta domain.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "18275829"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1016/j.immuni.2007.12.018")
        )
        reference = record.references[84]
        self.assertEqual(
            reference.authors,
            "Gras S., Saulquin X., Reiser J.B., Debeaupuis E., Echasserieau K., Kissenpfennig A., Legoux F., Chouquet A., Le Gorrec M., Machillot P., Neveu B., Thielens N., Malissen B., Bonneville M., Housset D.",
        )
        self.assertEqual(
            reference.title,
            "Structural bases for the affinity-driven selection of a public TCR against a dominant human cytomegalovirus epitope.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "19542454"))
        self.assertEqual(reference.references[1], ("DOI", "10.4049/jimmunol.0900556"))
        reference = record.references[85]
        self.assertEqual(
            reference.authors,
            "Kumar P., Vahedi-Faridi A., Saenger W., Ziegler A., Uchanska-Ziegler B.",
        )
        self.assertEqual(
            reference.title,
            "Conformational changes within the HLA-A1:MAGE-A1 complex induced by binding of a recombinant antibody fragment with TCR-like specificity.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "19177349"))
        self.assertEqual(reference.references[1], ("DOI", "10.1002/pro.4"))
        reference = record.references[86]
        self.assertEqual(
            reference.authors,
            "Liu J., Wu P., Gao F., Qi J., Kawana-Tachikawa A., Xie J., Vavricka C.J., Iwamoto A., Li T., Gao G.F.",
        )
        self.assertEqual(
            reference.title,
            "Novel immunodominant peptide presentation strategy: a featured HLA-A*2402-restricted cytotoxic T-lymphocyte epitope stabilized by intrachain hydrogen bonds from severe acute respiratory syndrome coronavirus nucleocapsid protein.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "20844028"))
        self.assertEqual(reference.references[1], ("DOI", "10.1128/jvi.01464-10"))
        reference = record.references[87]
        self.assertEqual(reference.authors, "Borbulevych O.Y., Do P., Baker B.M.")
        self.assertEqual(
            reference.title,
            "Structures of native and affinity-enhanced WT1 epitopes bound to HLA-A*0201: implications for WT1-based cancer therapeutics.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "20619457"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1016/j.molimm.2010.06.005")
        )
        reference = record.references[88]
        self.assertEqual(
            reference.authors,
            "McMahon R.M., Friis L., Siebold C., Friese M.A., Fugger L., Jones E.Y.",
        )
        self.assertEqual(
            reference.title,
            "Structure of HLA-A*0301 in complex with a peptide of proteolipid protein: insights into the role of HLA-A alleles in susceptibility to multiple sclerosis.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "21543847"))
        self.assertEqual(reference.references[1], ("DOI", "10.1107/s0907444911007888"))
        reference = record.references[89]
        self.assertEqual(
            reference.authors,
            "Zhang S., Liu J., Cheng H., Tan S., Qi J., Yan J., Gao G.F.",
        )
        self.assertEqual(
            reference.title,
            "Structural basis of cross-allele presentation by HLA-A*0301 and HLA-A*1101 revealed by two HIV-derived peptide complexes.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "21943705"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1016/j.molimm.2011.08.015")
        )
        reference = record.references[90]
        self.assertEqual(
            reference.authors,
            "Bulek A.M., Cole D.K., Skowera A., Dolton G., Gras S., Madura F., Fuller A., Miles J.J., Gostick E., Price D.A., Drijfhout J.W., Knight R.R., Huang G.C., Lissin N., Molloy P.E., Wooldridge L., Jakobsen B.K., Rossjohn J., Peakman M., Rizkallah P.J., Sewell A.K.",
        )
        self.assertEqual(
            reference.title,
            "Structural basis for the killing of human beta cells by CD8(+) T cells in type 1 diabetes.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "22245737"))
        self.assertEqual(reference.references[1], ("DOI", "10.1038/ni.2206"))
        reference = record.references[91]
        self.assertEqual(
            reference.authors,
            "Quinones-Parra S., Grant E., Loh L., Nguyen T.H., Campbell K.A., Tong S.Y., Miller A., Doherty P.C., Vijaykrishna D., Rossjohn J., Gras S., Kedzierska K.",
        )
        self.assertEqual(
            reference.title,
            "Preexisting CD8+ T-cell immunity to the H7N9 influenza A virus varies across ethnicities.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "24395804"))
        self.assertEqual(reference.references[1], ("DOI", "10.1073/pnas.1322229111"))
        reference = record.references[92]
        self.assertEqual(
            reference.authors,
            "Raman M.C., Rizkallah P.J., Simmons R., Donnellan Z., Dukes J., Bossi G., Le Provost G.S., Todorov P., Baston E., Hickman E., Mahon T., Hassan N., Vuidepot A., Sami M., Cole D.K., Jakobsen B.K.",
        )
        self.assertEqual(
            reference.title,
            "Direct molecular mimicry enables off-target cardiovascular toxicity by an enhanced affinity TCR designed for cancer immunotherapy.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "26758806"))
        self.assertEqual(reference.references[1], ("DOI", "10.1038/srep18851"))
        reference = record.references[93]
        self.assertEqual(
            reference.authors,
            "Song I., Gil A., Mishra R., Ghersi D., Selin L.K., Stern L.J.",
        )
        self.assertEqual(
            reference.title,
            "Broad TCR repertoire and diverse structural solutions for recognition of an immunodominant CD8+ T cell epitope.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "28250417"))
        self.assertEqual(reference.references[1], ("DOI", "10.1038/nsmb.3383"))
        reference = record.references[94]
        self.assertEqual(
            reference.authors,
            "LeHoang P., Ozdemir N., Benhamou A., Tabary T., Edelson C., Betuel H., Semiglia R., Cohen J.H.",
        )
        self.assertEqual(
            reference.title,
            "HLA-A29.2 subtype associated with birdshot retinochoroidopathy.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "1728143"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1016/s0002-9394(14)75749-6")
        )
        reference = record.references[95]
        self.assertEqual(
            reference.authors,
            "Fogdell-Hahn A., Ligers A., Groenning M., Hillert J., Olerup O.",
        )
        self.assertEqual(
            reference.title,
            "Multiple sclerosis: a modifying influence of HLA class I genes in an HLA class II associated autoimmune disease.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "10746785"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1034/j.1399-0039.2000.550205.x")
        )
        reference = record.references[96]
        self.assertEqual(reference.authors, "Nakanishi K., Inoko H.")
        self.assertEqual(
            reference.title,
            "Combination of HLA-A24, -DQA1*03, and -DR9 contributes to acute-onset and early complete beta-cell destruction in type 1 diabetes: longitudinal study of residual beta-cell function.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "16731854"))
        self.assertEqual(reference.references[1], ("DOI", "10.2337/db05-1049"))
        reference = record.references[97]
        self.assertEqual(
            reference.authors,
            "Skowera A., Ellis R.J., Varela-Calvino R., Arif S., Huang G.C., Van-Krinks C., Zaremba A., Rackham C., Allen J.S., Tree T.I., Zhao M., Dayan C.M., Sewell A.K., Unger W.W., Unger W., Drijfhout J.W., Ossendorp F., Roep B.O., Peakman M.",
        )
        self.assertEqual(
            reference.title,
            "CTLs are targeted to kill beta cells in patients with type 1 diabetes through recognition of a glucose-regulated preproinsulin epitope.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "18802479"))
        self.assertEqual(reference.references[1], ("DOI", "10.1172/jci35449"))
        reference = record.references[98]
        self.assertEqual(
            reference.authors,
            "Friese M.A., Jakobsen K.B., Friis L., Etzensperger R., Craner M.J., McMahon R.M., Jensen L.T., Huygelen V., Jones E.Y., Bell J.I., Fugger L.",
        )
        self.assertEqual(
            reference.title,
            "Opposing effects of HLA class I molecules in tuning autoreactive CD8+ T cells in multiple sclerosis.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "18953350"))
        self.assertEqual(reference.references[1], ("DOI", "10.1038/nm.1881"))
        reference = record.references[99]
        self.assertEqual(
            reference.authors,
            "Kronenberg D., Knight R.R., Estorninho M., Ellis R.J., Kester M.G., de Ru A., Eichmann M., Huang G.C., Powrie J., Dayan C.M., Skowera A., van Veelen P.A., Peakman M.",
        )
        self.assertEqual(
            reference.title,
            "Circulating preproinsulin signal peptide-specific CD8 T cells restricted by the susceptibility molecule HLA-A24 are expanded at onset of type 1 diabetes and kill beta-cells.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "22522618"))
        self.assertEqual(reference.references[1], ("DOI", "10.2337/db11-1520"))
        reference = record.references[100]
        self.assertEqual(
            reference.authors,
            "McCormack M., Alfirevic A., Bourgeois S., Farrell J.J., Kasperaviciute D., Carrington M., Sills G.J., Marson T., Jia X., de Bakker P.I., Chinthapalli K., Molokhia M., Johnson M.R., O'Connor G.D., Chaila E., Alhusaini S., Shianna K.V., Radtke R.A., Heinzen E.L., Walley N., Pandolfo M., Pichler W., Park B.K., Depondt C., Sisodiya S.M., Goldstein D.B., Deloukas P., Delanty N., Cavalleri G.L., Pirmohamed M.",
        )
        self.assertEqual(
            reference.title,
            "HLA-A*3101 and carbamazepine-induced hypersensitivity reactions in Europeans.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "21428769"))
        self.assertEqual(reference.references[1], ("DOI", "10.1056/nejmoa1013297"))
        reference = record.references[101]
        self.assertEqual(
            reference.authors,
            "Nakamura J., Meguro A., Ishii G., Mihara T., Takeuchi M., Mizuki Y., Yuda K., Yamane T., Kawagoe T., Ota M., Mizuki N.",
        )
        self.assertEqual(
            reference.title,
            "The association analysis between HLA-A*26 and Behcet's disease.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "30872678"))
        self.assertEqual(reference.references[1], ("DOI", "10.1038/s41598-019-40824-y"))
        reference = record.references[102]
        self.assertEqual(
            reference.authors,
            "Robinson J., Guethlein L.A., Cereb N., Yang S.Y., Norman P.J., Marsh S.G.E., Parham P.",
        )
        self.assertEqual(
            reference.title,
            "Distinguishing functional polymorphism from random variation in the sequences of >10,000 HLA-A, -B and -C alleles.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "28650991"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1371/journal.pgen.1006862")
        )

        # Check the two parsers agree on the essentials
        self.assertEqual(seq_record.seq, record.sequence)
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
        self.assertEqual(records[0].seq, seq_record.seq)
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

    def test_O23729(self):
        """Parsing SwissProt file O23729.txt."""
        filename = "O23729.txt"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)

        with open(datafile) as test_handle:
            seq_record = SeqIO.read(test_handle, "swiss")

        self.assertIsInstance(seq_record, SeqRecord)

        self.assertEqual(seq_record.id, "O23729")
        self.assertEqual(seq_record.name, "CHS3_BROFI")
        self.assertEqual(
            seq_record.description,
            "RecName: Full=Chalcone synthase 3; EC=2.3.1.74; AltName: Full=Naringenin-chalcone synthase 3;",
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
                "Streptophyta",
                "Embryophyta",
                "Tracheophyta",
                "Spermatophyta",
                "Magnoliopsida",
                "Liliopsida",
                "Asparagales",
                "Orchidaceae",
                "Epidendroideae",
                "Vandeae",
                "Adrorhizinae",
                "Bromheadia",
            ],
        )
        self.assertEqual(record.seqinfo, (394, 42942, "2F8D14AF4870BBB2"))

        self.assertEqual(len(record.features), 2)
        feature = record.features[0]
        self.assertEqual(feature.type, "CHAIN")
        self.assertEqual(feature.location.start, 0)
        self.assertEqual(feature.location.end, 394)
        self.assertEqual(feature.qualifiers["note"], "Chalcone synthase 3")
        self.assertEqual(feature.id, "PRO_0000215956")
        feature = record.features[1]
        self.assertEqual(feature.type, "ACT_SITE")
        self.assertEqual(feature.location.start, 164)
        self.assertEqual(feature.location.end, 165)
        self.assertEqual(
            feature.qualifiers["evidence"], "ECO:0000255|PROSITE-ProRule:PRU10023"
        )
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
        self.assertEqual(seq_record.seq, record.sequence)
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
        self.assertEqual(records[0].seq, seq_record.seq)
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
        self.assertEqual(seq_record.seq, record.sequence)
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
        self.assertEqual(records[0].seq, seq_record.seq)
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

    def test_P16235(self):
        """Parsing SwissProt file P16235.txt."""
        filename = "P16235.txt"
        # test the record parser

        datafile = os.path.join("SwissProt", filename)
        seq_record = SeqIO.read(datafile, "swiss")

        self.assertIsInstance(seq_record, SeqRecord)

        self.assertEqual(seq_record.id, "P16235")
        self.assertEqual(seq_record.name, "LSHR_RAT")
        self.assertEqual(
            seq_record.description,
            "RecName: Full=Lutropin-choriogonadotropic hormone receptor; Short=LH/CG-R; AltName: Full=Luteinizing hormone receptor; Short=LSH-R; Flags: Precursor;",
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
            record.accessions,
            ["P16235", "P70646", "Q63807", "Q63808", "Q63809", "Q6LDI7"],
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
                "Glires",
                "Rodentia",
                "Myomorpha",
                "Muroidea",
                "Muridae",
                "Murinae",
                "Rattus",
            ],
        )
        self.assertEqual(record.seqinfo, (700, 78036, "31807E73BAC94F1F"))

        self.assertEqual(len(record.features), 56)
        feature = record.features[0]
        self.assertEqual(feature.type, "SIGNAL")
        self.assertEqual(feature.location.start, 0)
        self.assertEqual(feature.location.end, 26)
        self.assertEqual(
            feature.qualifiers["evidence"],
            "ECO:0000269|PubMed:2601325,ECO:0000269|PubMed:2925659",
        )
        feature = record.features[1]
        self.assertEqual(feature.type, "CHAIN")
        self.assertEqual(feature.location.start, 26)
        self.assertEqual(feature.location.end, 700)
        self.assertEqual(
            feature.qualifiers["note"], "Lutropin-choriogonadotropic hormone receptor"
        )
        feature = record.features[2]
        self.assertEqual(feature.type, "TOPO_DOM")
        self.assertEqual(feature.location.start, 26)
        self.assertEqual(feature.location.end, 362)
        self.assertEqual(feature.qualifiers["note"], "Extracellular")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[3]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 362)
        self.assertEqual(feature.location.end, 390)
        self.assertEqual(feature.qualifiers["note"], "Helical; Name=1")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[4]
        self.assertEqual(feature.type, "TOPO_DOM")
        self.assertEqual(feature.location.start, 390)
        self.assertEqual(feature.location.end, 399)
        self.assertEqual(feature.qualifiers["note"], "Cytoplasmic")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[5]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 399)
        self.assertEqual(feature.location.end, 422)
        self.assertEqual(feature.qualifiers["note"], "Helical; Name=2")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[6]
        self.assertEqual(feature.type, "TOPO_DOM")
        self.assertEqual(feature.location.start, 422)
        self.assertEqual(feature.location.end, 443)
        self.assertEqual(feature.qualifiers["note"], "Extracellular")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[7]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 443)
        self.assertEqual(feature.location.end, 466)
        self.assertEqual(feature.qualifiers["note"], "Helical; Name=3")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[8]
        self.assertEqual(feature.type, "TOPO_DOM")
        self.assertEqual(feature.location.start, 466)
        self.assertEqual(feature.location.end, 486)
        self.assertEqual(feature.qualifiers["note"], "Cytoplasmic")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[9]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 486)
        self.assertEqual(feature.location.end, 509)
        self.assertEqual(feature.qualifiers["note"], "Helical; Name=4")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[10]
        self.assertEqual(feature.type, "TOPO_DOM")
        self.assertEqual(feature.location.start, 509)
        self.assertEqual(feature.location.end, 529)
        self.assertEqual(feature.qualifiers["note"], "Extracellular")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[11]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 529)
        self.assertEqual(feature.location.end, 551)
        self.assertEqual(feature.qualifiers["note"], "Helical; Name=5")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[12]
        self.assertEqual(feature.type, "TOPO_DOM")
        self.assertEqual(feature.location.start, 551)
        self.assertEqual(feature.location.end, 574)
        self.assertEqual(feature.qualifiers["note"], "Cytoplasmic")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[13]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 574)
        self.assertEqual(feature.location.end, 598)
        self.assertEqual(feature.qualifiers["note"], "Helical; Name=6")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[14]
        self.assertEqual(feature.type, "TOPO_DOM")
        self.assertEqual(feature.location.start, 598)
        self.assertEqual(feature.location.end, 609)
        self.assertEqual(feature.qualifiers["note"], "Extracellular")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[15]
        self.assertEqual(feature.type, "TRANSMEM")
        self.assertEqual(feature.location.start, 609)
        self.assertEqual(feature.location.end, 631)
        self.assertEqual(feature.qualifiers["note"], "Helical; Name=7")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[16]
        self.assertEqual(feature.type, "TOPO_DOM")
        self.assertEqual(feature.location.start, 631)
        self.assertEqual(feature.location.end, 700)
        self.assertEqual(feature.qualifiers["note"], "Cytoplasmic")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[17]
        self.assertEqual(feature.type, "REPEAT")
        self.assertEqual(feature.location.start, 125)
        self.assertEqual(feature.location.end, 150)
        self.assertEqual(feature.qualifiers["note"], "LRR 1")
        feature = record.features[18]
        self.assertEqual(feature.type, "REPEAT")
        self.assertEqual(feature.location.start, 151)
        self.assertEqual(feature.location.end, 175)
        self.assertEqual(feature.qualifiers["note"], "LRR 2")
        feature = record.features[19]
        self.assertEqual(feature.type, "REPEAT")
        self.assertEqual(feature.location.start, 175)
        self.assertEqual(feature.location.end, 200)
        self.assertEqual(feature.qualifiers["note"], "LRR 3")
        feature = record.features[20]
        self.assertEqual(feature.type, "REPEAT")
        self.assertEqual(feature.location.start, 201)
        self.assertEqual(feature.location.end, 224)
        self.assertEqual(feature.qualifiers["note"], "LRR 4")
        feature = record.features[21]
        self.assertEqual(feature.type, "REPEAT")
        self.assertEqual(feature.location.start, 224)
        self.assertEqual(feature.location.end, 248)
        self.assertEqual(feature.qualifiers["note"], "LRR 5")
        feature = record.features[22]
        self.assertEqual(feature.type, "REPEAT")
        self.assertEqual(feature.location.start, 249)
        self.assertEqual(feature.location.end, 271)
        self.assertEqual(feature.qualifiers["note"], "LRR 6")
        feature = record.features[23]
        self.assertEqual(feature.type, "MOD_RES")
        self.assertEqual(feature.location.start, 334)
        self.assertEqual(feature.location.end, 335)
        self.assertEqual(feature.qualifiers["note"], "Sulfotyrosine")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000250|UniProtKB:P22888")
        feature = record.features[24]
        self.assertEqual(feature.type, "LIPID")
        self.assertEqual(feature.location.start, 646)
        self.assertEqual(feature.location.end, 647)
        self.assertEqual(feature.qualifiers["note"], "S-palmitoyl cysteine")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000305|PubMed:7776964")
        feature = record.features[25]
        self.assertEqual(feature.type, "LIPID")
        self.assertEqual(feature.location.start, 647)
        self.assertEqual(feature.location.end, 648)
        self.assertEqual(feature.qualifiers["note"], "S-palmitoyl cysteine")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000305|PubMed:7776964")
        feature = record.features[26]
        self.assertEqual(feature.type, "CARBOHYD")
        self.assertEqual(feature.location.start, 102)
        self.assertEqual(feature.location.end, 103)
        self.assertEqual(feature.qualifiers["note"], "N-linked (GlcNAc...) asparagine")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[27]
        self.assertEqual(feature.type, "CARBOHYD")
        self.assertEqual(feature.location.start, 177)
        self.assertEqual(feature.location.end, 178)
        self.assertEqual(feature.qualifiers["note"], "N-linked (GlcNAc...) asparagine")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[28]
        self.assertEqual(feature.type, "CARBOHYD")
        self.assertEqual(feature.location.start, 198)
        self.assertEqual(feature.location.end, 199)
        self.assertEqual(feature.qualifiers["note"], "N-linked (GlcNAc...) asparagine")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[29]
        self.assertEqual(feature.type, "CARBOHYD")
        self.assertEqual(feature.location.start, 294)
        self.assertEqual(feature.location.end, 295)
        self.assertEqual(feature.qualifiers["note"], "N-linked (GlcNAc...) asparagine")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[30]
        self.assertEqual(feature.type, "CARBOHYD")
        self.assertEqual(feature.location.start, 302)
        self.assertEqual(feature.location.end, 303)
        self.assertEqual(feature.qualifiers["note"], "N-linked (GlcNAc...) asparagine")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[31]
        self.assertEqual(feature.type, "CARBOHYD")
        self.assertEqual(feature.location.start, 316)
        self.assertEqual(feature.location.end, 317)
        self.assertEqual(feature.qualifiers["note"], "N-linked (GlcNAc...) asparagine")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000255")
        feature = record.features[32]
        self.assertEqual(feature.type, "DISULFID")
        self.assertEqual(feature.location.start, 442)
        self.assertEqual(feature.location.end, 518)
        self.assertEqual(
            feature.qualifiers["evidence"], "ECO:0000255|PROSITE-ProRule:PRU00521"
        )
        feature = record.features[33]
        self.assertEqual(feature.type, "VAR_SEQ")
        self.assertEqual(feature.location.start, 82)
        self.assertEqual(feature.location.end, 132)
        self.assertEqual(feature.qualifiers["note"], "Missing (in isoform 1950)")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000305")
        feature = record.features[34]
        self.assertEqual(feature.type, "VAR_SEQ")
        self.assertEqual(feature.location.start, 132)
        self.assertEqual(feature.location.end, 157)
        self.assertEqual(feature.qualifiers["note"], "Missing (in isoform 1759)")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000305")
        feature = record.features[35]
        self.assertEqual(feature.type, "VAR_SEQ")
        self.assertEqual(feature.location.start, 183)
        self.assertEqual(feature.location.end, 700)
        self.assertEqual(feature.qualifiers["note"], "Missing (in isoform C2)")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000305")
        feature = record.features[36]
        self.assertEqual(feature.type, "VAR_SEQ")
        self.assertEqual(feature.location.start, 231)
        self.assertEqual(feature.location.end, 293)
        self.assertEqual(
            feature.qualifiers["note"],
            "Missing (in isoform EA2, isoform EB and isoform B1)",
        )
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000305")
        feature = record.features[37]
        self.assertEqual(feature.type, "VAR_SEQ")
        self.assertEqual(feature.location.start, 231)
        self.assertEqual(feature.location.end, 251)
        self.assertEqual(
            feature.qualifiers["note"],
            "DISSTKLQALPSHGLESIQT -> PCRATGWSPFRRSSPCLPTH (in isoform 2075)",
        )
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000305")
        feature = record.features[38]
        self.assertEqual(feature.type, "VAR_SEQ")
        self.assertEqual(feature.location.start, 251)
        self.assertEqual(feature.location.end, 700)
        self.assertEqual(feature.qualifiers["note"], "Missing (in isoform 2075)")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000305")
        feature = record.features[39]
        self.assertEqual(feature.type, "VAR_SEQ")
        self.assertEqual(feature.location.start, 293)
        self.assertEqual(feature.location.end, 367)
        self.assertEqual(
            feature.qualifiers["note"],
            "QNFSFSIFENFSKQCESTVRKADNETLYSAIFEENELSGWDYDYGFCSPKTLQCAPEPDAFNPCEDIMGYAFLR -> IFHFPFLKTSPNNAKAQLEKQITRRFIPPSLRRMNSVAGIMIMASVHPRHSNVLQNQMLSTPVKILWAMPSLGS (in isoform B1 and isoform B3)",
        )
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000305")
        feature = record.features[40]
        self.assertEqual(feature.type, "VAR_SEQ")
        self.assertEqual(feature.location.start, 293)
        self.assertEqual(feature.location.end, 294)
        self.assertEqual(feature.qualifiers["note"], "Q -> P (in isoform C1)")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000305")
        feature = record.features[41]
        self.assertEqual(feature.type, "VAR_SEQ")
        self.assertEqual(feature.location.start, 294)
        self.assertEqual(feature.location.end, 700)
        self.assertEqual(feature.qualifiers["note"], "Missing (in isoform C1)")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000305")
        feature = record.features[42]
        self.assertEqual(feature.type, "VAR_SEQ")
        self.assertEqual(feature.location.start, 320)
        self.assertEqual(feature.location.end, 342)
        self.assertEqual(
            feature.qualifiers["note"],
            "YSAIFEENELSGWDYDYGFCSP -> LHGALPAAHCLRGLPNKRPVL (in isoform 1834, isoform 1759 and isoform EB)",
        )
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000305")
        feature = record.features[43]
        self.assertEqual(feature.type, "VAR_SEQ")
        self.assertEqual(feature.location.start, 342)
        self.assertEqual(feature.location.end, 700)
        self.assertEqual(
            feature.qualifiers["note"],
            "Missing (in isoform 1834, isoform 1759 and isoform EB)",
        )
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000305")
        feature = record.features[44]
        self.assertEqual(feature.type, "VAR_SEQ")
        self.assertEqual(feature.location.start, 367)
        self.assertEqual(feature.location.end, 700)
        self.assertEqual(
            feature.qualifiers["note"], "Missing (in isoform B1 and isoform B3)"
        )
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000305")
        feature = record.features[45]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 81)
        self.assertEqual(feature.location.end, 82)
        self.assertEqual(feature.qualifiers["note"], "I -> M (in isoform 1950)")
        feature = record.features[46]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 178)
        self.assertEqual(feature.location.end, 179)
        self.assertEqual(feature.qualifiers["note"], "E -> G (in isoform 1759)")
        feature = record.features[47]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 232)
        self.assertEqual(feature.location.end, 233)
        self.assertEqual(feature.qualifiers["note"], "I -> T (in isoform 1950)")
        feature = record.features[48]
        self.assertEqual(feature.type, "VARIANT")
        self.assertEqual(feature.location.start, 645)
        self.assertEqual(feature.location.end, 646)
        self.assertEqual(feature.qualifiers["note"], "G -> S (in isoform 1950)")
        feature = record.features[49]
        self.assertEqual(feature.type, "MUTAGEN")
        self.assertEqual(feature.location.start, 408)
        self.assertEqual(feature.location.end, 409)
        self.assertEqual(
            feature.qualifiers["note"], "D->N: Significant reduction of binding."
        )
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000269|PubMed:1714448")
        feature = record.features[50]
        self.assertEqual(feature.type, "MUTAGEN")
        self.assertEqual(feature.location.start, 435)
        self.assertEqual(feature.location.end, 436)
        self.assertEqual(
            feature.qualifiers["note"], "D->N: No change in binding or cAMP production."
        )
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000269|PubMed:1714448")
        feature = record.features[51]
        self.assertEqual(feature.type, "MUTAGEN")
        self.assertEqual(feature.location.start, 454)
        self.assertEqual(feature.location.end, 455)
        self.assertEqual(
            feature.qualifiers["note"], "E->Q: No change in binding or cAMP production."
        )
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000269|PubMed:1714448")
        feature = record.features[52]
        self.assertEqual(feature.type, "MUTAGEN")
        self.assertEqual(feature.location.start, 581)
        self.assertEqual(feature.location.end, 582)
        self.assertEqual(
            feature.qualifiers["note"], "D->N: No change in binding or cAMP production."
        )
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000269|PubMed:1714448")
        feature = record.features[53]
        self.assertEqual(feature.type, "MUTAGEN")
        self.assertEqual(feature.location.start, 646)
        self.assertEqual(feature.location.end, 647)
        self.assertEqual(
            feature.qualifiers["note"],
            "C->A: Trapped intracellularly and does not appear to become mature; when associated with A-648.",
        )
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000269|PubMed:7776964")
        feature = record.features[54]
        self.assertEqual(feature.type, "MUTAGEN")
        self.assertEqual(feature.location.start, 647)
        self.assertEqual(feature.location.end, 648)
        self.assertEqual(
            feature.qualifiers["note"],
            "C->A: Trapped intracellularly and does not appear to become mature; when associated with A-647.",
        )
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000269|PubMed:7776964")
        feature = record.features[55]
        self.assertEqual(feature.type, "CONFLICT")
        self.assertEqual(feature.location.start, 32)
        self.assertEqual(feature.location.end, 33)
        self.assertEqual(feature.qualifiers["note"], "R -> L (in Ref. 9; AA sequence)")
        self.assertEqual(feature.qualifiers["evidence"], "ECO:0000305")
        self.assertIsNone(feature.id)
        self.assertEqual(len(record.references), 11)
        reference = record.references[0]
        self.assertEqual(
            reference.authors,
            "McFarland K.C., Sprengel R., Phillips H.S., Koehler M., Rosemblit N., Nikolics K., Segaloff D.L., Seeburg P.H.",
        )
        self.assertEqual(
            reference.title,
            "Lutropin-choriogonadotropin receptor: an unusual member of the G protein-coupled receptor family.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "2502842"))
        self.assertEqual(reference.references[1], ("DOI", "10.1126/science.2502842"))
        reference = record.references[1]
        self.assertEqual(
            reference.authors,
            "Aatsinki J.T., Pietila E.M., Lakkakorpi J.T., Rajaniemi H.J.",
        )
        self.assertEqual(
            reference.title,
            "Expression of the LH/CG receptor gene in rat ovarian tissue is regulated by an extensive alternative splicing of the primary transcript.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "1353463"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1016/0303-7207(92)90079-l")
        )
        reference = record.references[2]
        self.assertEqual(reference.authors, "Koo Y.B., Slaughter R.G., Ji T.H.")
        self.assertEqual(
            reference.title,
            "Structure of the luteinizing hormone receptor gene and multiple exons of the coding sequence.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "2019252"))
        self.assertEqual(reference.references[1], ("DOI", "10.1210/endo-128-5-2297"))
        reference = record.references[3]
        self.assertEqual(reference.authors, "Bernard M.P., Myers R.V., Moyle W.R.")
        self.assertEqual(
            reference.title,
            "Cloning of rat lutropin (LH) receptor analogs lacking the soybean lectin domain.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "1976554"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1016/0303-7207(90)90034-6")
        )
        reference = record.references[4]
        self.assertEqual(
            reference.authors, "Segaloff D.L., Sprengel R., Nikolics K., Ascoli M."
        )
        self.assertEqual(
            reference.title, "Structure of the lutropin/choriogonadotropin receptor."
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "2281186"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1016/b978-0-12-571146-3.50014-6")
        )
        reference = record.references[5]
        self.assertEqual(
            reference.authors,
            "Tsai-Morris C.H., Buczko E., Wang W., Xie X.-Z., Dufau M.L.",
        )
        self.assertEqual(
            reference.title,
            "Structural organization of the rat luteinizing hormone (LH) receptor gene.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "2040640"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1016/s0021-9258(18)99170-2")
        )
        reference = record.references[6]
        self.assertEqual(
            reference.authors, "Tsai-Morris C.H., Buczko E., Wang W., Dufau M.L."
        )
        self.assertEqual(
            reference.title,
            "Intronic nature of the rat luteinizing hormone receptor gene defines a soluble receptor subspecies with hormone binding activity.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "2174034"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1016/s0021-9258(17)45380-4")
        )
        reference = record.references[7]
        self.assertEqual(
            reference.authors,
            "Dufau M.L., Minegishi T., Buczko E.S., Delgado C.J., Zhang R.",
        )
        self.assertEqual(
            reference.title,
            "Characterization and structure of ovarian and testicular LH/hCG receptors.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "2601325"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1016/0022-4731(89)90482-2")
        )
        reference = record.references[8]
        self.assertEqual(reference.authors, "Roche P.C., Ryan R.J.")
        self.assertEqual(
            reference.title,
            "Purification, characterization, and amino-terminal sequence of rat ovarian receptor for luteinizing hormone/human choriogonadotropin.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "2925659"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1016/s0021-9258(18)83790-5")
        )
        reference = record.references[9]
        self.assertEqual(reference.authors, "Ji I., Ji T.H.")
        self.assertEqual(
            reference.title,
            "Asp383 in the second transmembrane domain of the lutropin receptor is important for high affinity hormone binding and cAMP production.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "1714448"))
        self.assertEqual(
            reference.references[1], ("DOI", "10.1016/s0021-9258(18)98570-4")
        )
        reference = record.references[10]
        self.assertEqual(reference.authors, "Zhu H., Wang H., Ascoli M.")
        self.assertEqual(
            reference.title,
            "The lutropin/choriogonadotropin receptor is palmitoylated at intracellular cysteine residues.",
        )
        self.assertEqual(len(reference.references), 2)
        self.assertEqual(reference.references[0], ("PubMed", "7776964"))
        self.assertEqual(reference.references[1], ("DOI", "10.1210/mend.9.2.7776964"))

        # Check the two parsers agree on the essentials
        self.assertEqual(seq_record.seq, record.sequence)
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
        self.assertEqual(records[0].seq, seq_record.seq)
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
        self.assertEqual(seq_record.seq, record.sequence)
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
        self.assertEqual(records[0].seq, seq_record.seq)
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
        self.assertEqual(seq_record.seq, record.sequence)
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
        self.assertEqual(records[0].seq, seq_record.seq)
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
            repr(seq_record.seq), "Seq('MTQSNPNEQNVELNRTSLYWGLLLIFVLAVLFSNYFFN')"
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
        self.assertEqual(seq_record.seq, record.sequence)
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
        self.assertEqual(records[0].seq, seq_record.seq)
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
        self.assertEqual(seq_record.seq, record.sequence)
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
        self.assertEqual(records[0].seq, seq_record.seq)
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

    def test_Q7Z739(self):
        """Parsing SwissProt file Q7Z739.txt, which has new qualifiers for ligands from Uniprot version 2022_03."""
        filename = "Q7Z739.txt"
        datafile = os.path.join("SwissProt", filename)
        with open(datafile) as test_handle:
            record = SwissProt.read(test_handle)
        self.assertEqual(
            record.gene_name,
            [
                {
                    "Name": "YTHDF3 {ECO:0000303|PubMed:28106072, ECO:0000312|HGNC:HGNC:26465}"
                }
            ],
        )
        # Check the new ligand feature
        self.assertEqual(record.features[10].qualifiers["ligand"], "RNA")
        self.assertEqual(
            record.features[10].qualifiers["ligand_id"], "ChEBI:CHEBI:33697"
        )
        self.assertEqual(
            record.features[10].qualifiers["ligand_part"],
            "N(6)-methyladenosine 5'-phosphate residue",
        )
        self.assertEqual(
            record.features[10].qualifiers["ligand_part_id"], "ChEBI:CHEBI:74449"
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
