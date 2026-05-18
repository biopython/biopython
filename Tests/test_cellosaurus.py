# Copyright 2016 by Stephen Marshall.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for cellosaurus module."""

import unittest

from Bio.ExPASy import cellosaurus


class TestCellosaurus(unittest.TestCase):
    def test_read(self):
        """Test read function."""
        with open("Cellosaurus/cell_lines_1.txt") as handle:
            record = cellosaurus.read(handle)
        self.assertEqual(record["ID"], "#15310-LN")
        self.assertEqual(record["AC"], "CVCL_E548")
        self.assertEqual(
            record["SY"],
            "15310-LN; TER461; TER-461; Ter 461; TER479; TER-479; Ter 479; Extract 519",
        )
        self.assertEqual(record["DR"][0], ("dbMHC", "48439"))
        self.assertEqual(record["DR"][1], ("ECACC", "94050311"))
        self.assertEqual(record["DR"][2], ("IHW", "IHW09326"))
        self.assertEqual(record["DR"][3], ("IPD-IMGT/HLA", "10074"))
        self.assertEqual(record["DR"][4], ("Wikidata", "Q54398957"))
        self.assertEqual(
            record["WW"][0], "http://pathology.ucla.edu/workfiles/360cx.pdf"
        )
        self.assertEqual(
            record["WW"][1], "http://pathology.ucla.edu/workfiles/370cx.pdf"
        )
        self.assertEqual(
            record["CC"][0],
            "Part of: 12th International Histocompatibility Workshop (12IHW) "
            "cell line panel.",
        )
        self.assertEqual(record["CC"][1], "Population: Caucasian; French Canadian.")
        self.assertEqual(
            record["CC"][2],
            "HLA typing: A*03,25; B*37:01:01:01,47:01:01:03; C*06; DPA1*01; DPB1*04:01:01; DQA1*01:01:01:01,"
            "01:03:01; DQB1*05:01:01;06:03:01:02; DRB1*01:01:01,14:17; DRB3*01:01 (IPD-IMGT/HLA=10074).",
        )
        self.assertEqual(
            record["CC"][3],
            "Transformant: NCBI_TaxID; 10376; Epstein-Barr virus (EBV).",
        )
        self.assertEqual(
            record["CC"][4],
            "Derived from site: In situ; Peripheral blood; UBERON=UBERON_0000178.",
        )
        self.assertEqual(record["CC"][5], "Cell type: B-cell; CL=CL_0000236.")
        self.assertEqual(record["OX"][0], "NCBI_TaxID=9606; ! Homo sapiens (Human)")
        self.assertEqual(record["SX"], "Female")
        self.assertEqual(record["AG"], "Age unspecified")
        self.assertEqual(record["CA"], "Transformed cell line")
        self.assertEqual(
            record["DT"], "Created: 22-10-12; Last updated: 30-01-24; Version: 18"
        )

    def test_parse(self):
        """Test parsing function."""
        with open("Cellosaurus/cell_lines_2.txt") as handle:
            records = cellosaurus.parse(handle)
            record = next(records)
            self.assertEqual(record["ID"], "XP3OS")
            self.assertEqual(record["AC"], "CVCL_3245")
            self.assertEqual(record["AS"], "CVCL_F511")
            self.assertEqual(
                record["SY"], "Xeroderma Pigmentosum 3 OSaka; GM04314; GM04314B; GM4314"
            )
            self.assertEqual(record["DR"][0], ("CLO", "CLO_0019557"))
            self.assertEqual(record["DR"][1], ("Coriell", "GM04314"))
            self.assertEqual(record["DR"][2], ("GEO", "GSM1338611"))
            self.assertEqual(record["DR"][3], ("JCRB", "JCRB0303"))
            self.assertEqual(record["DR"][4], ("JCRB", "KURB1002"))
            self.assertEqual(record["DR"][5], ("JCRB", "KURB1003"))
            self.assertEqual(record["DR"][6], ("JCRB", "KURB1004"))
            self.assertEqual(record["DR"][7], ("Wikidata", "Q54994928"))
            self.assertEqual(record["RX"][0], "CelloPub=CLPUB00447;")
            self.assertEqual(record["RX"][1], "PubMed=832273;")
            self.assertEqual(record["RX"][2], "PubMed=1372102;")
            self.assertEqual(record["RX"][3], "PubMed=1702221;")
            self.assertEqual(record["RX"][4], "PubMed=2570806;")
            self.assertEqual(record["RX"][5], "PubMed=7000335;")
            self.assertEqual(record["RX"][6], "PubMed=7830260;")
            self.assertEqual(record["RX"][7], "PubMed=9671271;")
            self.assertEqual(record["RX"][8], "PubMed=27197874;")
            self.assertEqual(len(record["CC"]), 6)
            self.assertEqual(len(record["ST"]), 17)
            self.assertEqual(record["ST"][0], "Source(s): JCRB; PubMed=27197874")
            self.assertEqual(record["ST"][1], "Amelogenin: X")
            self.assertEqual(record["ST"][2], "CSF1PO: 10,11")
            self.assertEqual(record["ST"][3], "D13S317: 9,11")
            self.assertEqual(record["ST"][4], "D16S539: 9,12")
            self.assertEqual(record["ST"][5], "D18S51: 13")
            self.assertEqual(record["ST"][6], "D21S11: 29,30")
            self.assertEqual(record["ST"][7], "D3S1358: 15,16")
            self.assertEqual(record["ST"][8], "D5S818: 10,11")
            self.assertEqual(record["ST"][9], "D7S820: 11,12")
            self.assertEqual(record["ST"][10], "D8S1179: 13,15")
            self.assertEqual(record["ST"][11], "FGA: 22,23")
            self.assertEqual(record["ST"][12], "Penta D: 9")
            self.assertEqual(record["ST"][13], "Penta E: 14,17")
            self.assertEqual(record["ST"][14], "TH01: 7")
            self.assertEqual(record["ST"][15], "TPOX: 8,11")
            self.assertEqual(record["ST"][16], "vWA: 14,16")
            self.assertEqual(
                record["DI"][0],
                "NCIt; C3965; Xeroderma pigmentosum, complementation group A",
            )
            self.assertEqual(
                record["DI"][1],
                "ORDO; Orphanet_910; Xeroderma pigmentosum",
            )
            self.assertEqual(record["OX"][0], "NCBI_TaxID=9606; ! Homo sapiens (Human)")
            self.assertEqual(record["SX"], "Female")
            self.assertEqual(record["AG"], "5Y")
            self.assertEqual(record["CA"], "Finite cell line")
            self.assertEqual(
                record["DT"], "Created: 04-04-12; Last updated: 29-06-23; Version: 20"
            )

            record = next(records)
            self.assertEqual(record["ID"], "1-5c-4")
            self.assertEqual(record["AC"], "CVCL_2260")
            self.assertEqual(
                record["SY"],
                "Clone 1-5c-4; Clone 1-5c-4 WKD of Chang Conjunctiva; "
                "Wong-Kilbourne derivative of Chang conjunctiva; ChWK",
            )
            self.assertEqual(len(record["DR"]), 11)
            self.assertEqual(record["DR"][0], ("CLO", "CLO_0002500"))
            self.assertEqual(record["DR"][1], ("CLO", "CLO_0002501"))
            self.assertEqual(record["DR"][2], ("CLDB", "cl793"))
            self.assertEqual(record["DR"][3], ("CLDB", "cl794"))
            self.assertEqual(record["DR"][4], ("CLDB", "cl795"))
            self.assertEqual(record["DR"][5], ("ATCC", "CCL-20.2"))
            self.assertEqual(record["DR"][6], ("BioSample", "SAMN03151673"))
            self.assertEqual(record["DR"][7], ("ECACC", "88021103"))
            self.assertEqual(record["DR"][8], ("IZSLER", "BS CL 93"))
            self.assertEqual(record["DR"][9], ("KCLB", "10020.2"))
            self.assertEqual(record["DR"][10], ("Wikidata", "Q54399838"))
            self.assertEqual(record["RX"][0], "DOI=10.1007/BF02618370;")
            self.assertEqual(record["RX"][1], "PubMed=566722;")
            self.assertEqual(record["RX"][2], "PubMed=19630270;")
            self.assertEqual(record["RX"][3], "PubMed=20143388;")
            self.assertEqual(
                record["WW"][0],
                "https://iclac.org/wp-content/uploads/Cross-Contaminations_v12_distribution.xlsx",
            )
            self.assertEqual(len(record["CC"]), 6)
            self.assertEqual(
                record["CC"][0],
                "Problematic cell line: Contaminated. Shown to be a HeLa derivative (PubMed=566722; PubMed=20143388). "
                "Originally thought to originate from the conjunctiva of a child.",
            )
            self.assertEqual(
                record["CC"][1],
                "Registration: International Cell Line Authentication Committee, Register of Misidentified Cell Lines; ICLAC-00298.",
            )
            self.assertEqual(record["ST"][0], "Source(s): ATCC; KCLB")
            self.assertEqual(record["ST"][1], "Amelogenin: X")
            self.assertEqual(record["ST"][2], "CSF1PO: 9,10")
            self.assertEqual(record["ST"][3], "D13S317: 12,13.3")
            self.assertEqual(record["ST"][4], "D16S539: 9,10")
            self.assertEqual(record["ST"][5], "D3S1358: 15,18")
            self.assertEqual(record["ST"][6], "D5S818: 11,12")
            self.assertEqual(record["ST"][7], "D7S820: 8,12")
            self.assertEqual(record["ST"][8], "FGA: 18,21")
            self.assertEqual(record["ST"][9], "TH01: 7")
            self.assertEqual(record["ST"][10], "TPOX: 8,12")
            self.assertEqual(record["ST"][11], "vWA: 16,18")
            self.assertEqual(
                record["DI"][0],
                "NCIt; C27677; Human papillomavirus-related cervical adenocarcinoma",
            )
            self.assertEqual(record["OX"][0], "NCBI_TaxID=9606; ! Homo sapiens (Human)")
            self.assertEqual(record["HI"][0], "CVCL_0030 ! HeLa")
            self.assertEqual(record["SX"], "Female")
            self.assertEqual(record["AG"], "30Y6M")
            self.assertEqual(record["CA"], "Cancer cell line")
            self.assertEqual(
                record["DT"], "Created: 04-04-12; Last updated: 05-10-23; Version: 25"
            )
            self.assertRaises(StopIteration, next, records)

    def test__str__(self):
        """Test string function."""
        with open("Cellosaurus/cell_lines_3.txt") as handle:
            record = cellosaurus.read(handle)
        text = (
            "ID: ZZ-R 127 AC: CVCL_5418 AS:  SY: ZZ-R DR: [('CCLV', 'CCLV-RIE 0127'), ('Wikidata', 'Q54996143')] RX: "
            "['PubMed=19656987;', 'PubMed=19941903;', 'PubMed=26082457;', 'PubMed=27795464;', 'PubMed=32851014;', "
            "'PubMed=34737324;'] WW: [] CC: ['Virology: Susceptible to infection by foot-and-mouth disease virus.', "
            "'Derived from site: In situ; Fetal oral cavity, tongue; UBERON=UBERON_0010056.'] ST: [] DI: [] OX: ["
            "'NCBI_TaxID=9925; ! Capra hircus (Goat)'] HI: [] OI: [] SX: Sex unspecified AG: Fetus CA: Spontaneously "
            "immortalized cell line DT: Created: 04-04-12; Last updated: 29-06-23; Version: 15"
        )
        self.assertEqual(str(record), text)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=3)
    unittest.main(testRunner=runner)
