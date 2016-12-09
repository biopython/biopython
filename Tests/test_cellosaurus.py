# Copyright 2016 by Stephen Marshall.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
import unittest

from Bio.ExPASy import cellosaurus


class TestCellosaurus(unittest.TestCase):
    def test_read(self):
        """Test read function"""
        handle = open("Cellosaurus/cell_lines_1.txt")
        record = cellosaurus.read(handle)
        self.assertEquals(record["ID"], '#15310-LN')
        self.assertEqual(record["AC"], 'CVCL_E548')
        self.assertEqual(record["SY"], '15310-LN; TER461')
        self.assertEqual(record["DR"][0], tuple(('dbMHC', '48439')))
        self.assertEqual(record["DR"][1], tuple(('ECACC', '94050311')))
        self.assertEqual(record["DR"][2], tuple(('IHW', 'IHW9326')))
        self.assertEqual(record["DR"][3], tuple(('IMGT/HLA', '10074')))
        self.assertEqual(record["WW"][0], "http://bioinformatics."
                                          "hsanmartino.it"
                                          "/ecbr/cl326.html")
        self.assertEqual(record["CC"][0],
                         "Part of: 12th International "
                         "Histocompatibility Workshop "
                         "(12IHW) cell line panel.")
        self.assertEqual(record["CC"][1],
                         'Transformant: EBV.')
        self.assertEqual(record["OX"][0],
                         'NCBI_TaxID=9606; ! Homo sapiens')
        self.assertEqual(record["SX"], 'Female')
        self.assertEqual(record["CA"], 'Transformed cell line')
        handle.close()

    def test_parse(self):
        """Test parsing function"""
        handle = open("Cellosaurus/cell_lines_2.txt")
        records = cellosaurus.parse(handle)
        record = next(records)
        self.assertEqual(record["ID"], 'XP3OS')
        self.assertEqual(record["AC"], 'CVCL_3245')
        self.assertEqual(record["AS"], 'CVCL_F511')
        self.assertEqual(record["SY"], 'XP30S; GM04314')
        self.assertEqual(record["DR"][0], tuple(('CLO', 'CLO_0019557')))
        self.assertEqual(record["DR"][1], tuple(('Coriell', 'GM04314')))
        self.assertEqual(record["DR"][2], tuple(('JCRB', 'JCRB0303')))
        self.assertEqual(record["DR"][3], tuple(('JCRB', 'KURB1002')))
        self.assertEqual(record["DR"][4], tuple(('JCRB', 'KURB1003')))
        self.assertEqual(record["DR"][5], tuple(('JCRB', 'KURB1004')))
        self.assertEqual(record["RX"][0], 'PubMed=1372102;')
        self.assertEqual(len(record["ST"]), 10)
        self.assertEqual(record["ST"][0], 'Source(s): JCRB')
        self.assertEqual(record["ST"][1], 'Amelogenin: X')
        self.assertEqual(record["ST"][2], 'CSF1PO: 10,11')
        self.assertEqual(record["ST"][3], 'D13S317: 9,11')
        self.assertEqual(record["ST"][4], 'D16S539: 9,12')
        self.assertEqual(record["ST"][5], 'D5S818: 10,11')
        self.assertEqual(record["ST"][6], 'D7S820: 11,12')
        self.assertEqual(record["ST"][7], 'TH01: 7')
        self.assertEqual(record["ST"][8], 'TPOX: 8,11')
        self.assertEqual(record["ST"][9], 'vWA: 14,16')
        self.assertEqual(record["DI"][0], "NCIt; C3965; Xeroderma pigmentosum,"
                                          " complementation group A")
        self.assertEqual(record["OX"][0], 'NCBI_TaxID=9606; ! Homo sapiens')
        self.assertEqual(record["SX"], 'Female')
        self.assertEqual(record["CA"], 'Finite cell line')
        record = next(records)
        self.assertEqual(record["ID"], '1-5c-4')
        self.assertEqual(record["AC"], 'CVCL_2260')
        self.assertEqual(record["SY"],
                         "Clone 1-5c-4; Clone 1-5c-4 WKD of "
                         "Chang Conjunctiva; "
                         "Wong-Kilbourne derivative of "
                         "Chang conjunctiva; ChWK")
        self.assertEqual(len(record["DR"]), 10)
        self.assertEqual(record["DR"][0], tuple(('CLO', 'CLO_0002500')))
        self.assertEqual(record["DR"][1], tuple(('CLO', 'CLO_0002501')))
        self.assertEqual(record["DR"][2], tuple(('CLDB', 'cl793')))
        self.assertEqual(record["DR"][3], tuple(('CLDB', 'cl794')))
        self.assertEqual(record["DR"][4], tuple(('CLDB', 'cl795')))
        self.assertEqual(record["DR"][5], tuple(('ATCC', 'CCL-20.2')))
        self.assertEqual(record["DR"][6], tuple(('BioSample', 'SAMN03151673')))
        self.assertEqual(record["DR"][7], tuple(('ECACC', '88021103')))
        self.assertEqual(record["DR"][8], tuple(('IZSLER', 'BS CL 93')))
        self.assertEqual(record["DR"][9], tuple(('KCLB', '10020.2')))
        self.assertEqual(record["RX"][0], 'PubMed=566722;')
        self.assertEqual(record["RX"][1], 'PubMed=19630270;')
        self.assertEqual(record["RX"][2], 'PubMed=20143388;')
        self.assertEqual(record["WW"][0], "http://iclac.org/"
                                          "wp-content/uploads/"
                                          "Cross-Contaminations-v7_2.pdf")
        self.assertEqual(record["CC"][0],
                         "Problematic cell line: Contaminated. "
                         "Shown to be a HeLa derivative "
                         "(PubMed 566722, PubMed 20143388).")
        self.assertEqual(record["CC"][1], 'Omics: Transcriptome analysis.')
        self.assertEqual(record["ST"][0], 'Source(s): ATCC; KCLB')
        self.assertEqual(record["ST"][1], 'Amelogenin: X')
        self.assertEqual(record["ST"][2], 'CSF1PO: 9,10')
        self.assertEqual(record["ST"][3], 'D13S317: 12,13.3')
        self.assertEqual(record["ST"][4], 'D16S539: 9,10')
        self.assertEqual(record["ST"][5], 'D3S1358: 15,18')
        self.assertEqual(record["ST"][6], 'D5S818: 11,12')
        self.assertEqual(record["ST"][7], 'D7S820: 8,12')
        self.assertEqual(record["ST"][8], 'FGA: 18,21')
        self.assertEqual(record["ST"][9], 'TH01: 7')
        self.assertEqual(record["ST"][10], 'TPOX: 8,12')
        self.assertEqual(record["ST"][11], 'vWA: 16,18')
        self.assertEqual(record["DI"][0], "NCIt; C4029; Cervical "
                                          "adenocarcinoma")
        self.assertEqual(record["OX"][0], 'NCBI_TaxID=9606; ! Homo sapiens')
        self.assertEqual(record["HI"][0], 'CVCL_0030 ! HeLa')
        self.assertEqual(record["SX"], 'Female')
        self.assertEqual(record["CA"], 'Cancer cell line')
        self.assertRaises(StopIteration, next, records)
        handle.close()

    def test__str__(self):
        """Test string function"""
        handle = open("Cellosaurus/cell_lines_3.txt")
        record = cellosaurus.read(handle)
        input = "ID: ZZ-R 127 AC: CVCL_5418 AS:  SY: ZZ-R DR: [('CCLV', 'CCLV-RIE 0127')] " \
                "RX: ['PubMed=19656987;', 'PubMed=19941903;'] WW: [] CC: [] ST: [] DI: [] " \
                "OX: ['NCBI_TaxID=9925; ! Capra hircus'] HI: [] OI: [] SX:  CA: " \
                "Spontaneously immortalized cell line"
        self.assertEquals(record.__str__(), input)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=3)
    unittest.main(testRunner=runner)
