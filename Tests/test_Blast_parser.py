# Copyright 2005 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Test the Blast XML parser."""

import os
import unittest

import numpy as np

from Bio import Blast
from Bio.SeqRecord import SeqRecord


class TestBlastXML(unittest.TestCase):
    """Tests for the Blast XML parser."""

    def test_xml_2212L_blastp_001(self):
        """Parsing BLASTP 2.2.12, gi|49176427|ref|NP_418280.3| (xml_2212L_blastp_001)."""
        filename = "xml_2212L_blastp_001.xml"
        datafile = os.path.join("Blast", filename)
        with open(datafile, "rb") as handle:
            records = Blast.parse(handle)
            self.assertEqual(records.program, "blastp")
            self.assertEqual(records.version, "BLASTP 2.2.12 [Aug-07-2005]")
            self.assertEqual(
                records.reference,
                'Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schäffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.',
            )
            self.assertEqual(records.db, "nr")
            self.assertIsInstance(records.query, SeqRecord)
            self.assertEqual(records.query.id, "gi|49176427|ref|NP_418280.3|")
            self.assertEqual(
                records.query.description,
                "component of Sec-independent translocase [Escherichia coli K12]",
            )
            self.assertEqual(repr(records.query.seq), "Seq(None, length=103)")
            self.assertEqual(len(records.param), 5)
            self.assertEqual(records.param["matrix"], "BLOSUM62")
            self.assertAlmostEqual(records.param["expect"], 10)
            self.assertEqual(records.param["gap-open"], 11)
            self.assertEqual(records.param["gap-extend"], 1)
            self.assertEqual(records.param["filter"], "L;")
            record = next(records)
            self.assertRaises(StopIteration, next, records)
        self.check_xml_2212L_blastp_001(record)

        with open(datafile, "rb") as handle:
            record = Blast.read(handle)
            self.check_xml_2212L_blastp_001(record)

    def check_xml_2212L_blastp_001(self, record):
        self.assertEqual(record.num, 1)
        self.assertIsInstance(record.query, SeqRecord)
        self.assertEqual(record.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            record.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(repr(record.query.seq), "Seq(None, length=103)")
        self.assertEqual(len(record.stat), 7)
        self.assertEqual(record.stat["db-num"], 2934173)
        self.assertEqual(record.stat["db-len"], 1011751523)
        self.assertEqual(record.stat["hsp-len"], 0)
        self.assertAlmostEqual(record.stat["eff-space"], 0)
        self.assertAlmostEqual(record.stat["kappa"], 0.041)
        self.assertAlmostEqual(record.stat["lambda"], 0.267)
        self.assertAlmostEqual(record.stat["entropy"], 0.14)
        hits = record.hits
        self.assertEqual(len(hits), 212)

        hit = hits[0]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(hit.target.name, "NP_418280")
        self.assertEqual(
            hit.target.description,
            "component of Sec-independent translocase [Escherichia coli K12] >gi|26250604|ref|NP_756644.1| Sec-independent protein translocase protein tatA [Escherichia coli CFT073] >gi|30064867|ref|NP_839038.1| hypothetical protein S3840 [Shigella flexneri 2a str. 2457T] >gi|24115132|ref|NP_709642.1| hypothetical protein SF3914 [Shigella flexneri 2a str. 301] >gi|24054404|gb|AAN45349.1| orf, conserved hypothetical protein [Shigella flexneri 2a str. 301] >gi|2367310|gb|AAC76839.1| component of Sec-independent translocase [Escherichia coli K12] >gi|30043127|gb|AAP18849.1| hypothetical protein S3840 [Shigella flexneri 2a str. 2457T] >gi|26111035|gb|AAN83218.1| Sec-independent protein translocase protein tatA [Escherichia coli CFT073] >gi|3193217|gb|AAC19240.1| MttA1 [Escherichia coli] >gi|7444818|pir||E65188 hypothetical 11.3 kD protein in udp-rfaH intergenic region - Escherichia coli (strain K-12)",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=103)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 469.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 185.267)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.20576e-46)
        self.assertEqual(hsp.annotations["identity"], 103)
        self.assertEqual(hsp.annotations["positive"], 103)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0, 103],
                          [  0, 103]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 103))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MRLCLIIIYHRGTCMGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFK...EQV')",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MRLCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFK...EQV')",
        )
        self.assertEqual(hsp.target.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(hsp.target.name, "NP_418280")
        self.assertEqual(
            hsp.target.description,
            "component of Sec-independent translocase [Escherichia coli K12] >gi|26250604|ref|NP_756644.1| Sec-independent protein translocase protein tatA [Escherichia coli CFT073] >gi|30064867|ref|NP_839038.1| hypothetical protein S3840 [Shigella flexneri 2a str. 2457T] >gi|24115132|ref|NP_709642.1| hypothetical protein SF3914 [Shigella flexneri 2a str. 301] >gi|24054404|gb|AAN45349.1| orf, conserved hypothetical protein [Shigella flexneri 2a str. 301] >gi|2367310|gb|AAC76839.1| component of Sec-independent translocase [Escherichia coli K12] >gi|30043127|gb|AAP18849.1| hypothetical protein S3840 [Shigella flexneri 2a str. 2457T] >gi|26111035|gb|AAN83218.1| Sec-independent protein translocase protein tatA [Escherichia coli CFT073] >gi|3193217|gb|AAC19240.1| MttA1 [Escherichia coli] >gi|7444818|pir||E65188 hypothetical 11.3 kD protein in udp-rfaH intergenic region - Escherichia coli (strain K-12)",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MRLCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|491764         0 MRLCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDD
                  0 ||||||||||||||||||||||...........|||||||||||||||||||||||||||
gi|491764         0 MRLCLIIIYHRGTCMGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDD

gi|491764        60 EPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV 103
                 60 ||||||||||||||||||||||||||||||||||||||||||| 103
gi|491764        60 EPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV 103
""",
        )
        hit = hits[1]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|15804428|ref|NP_290468.1|")
        self.assertEqual(hit.target.name, "NP_290468")
        self.assertEqual(
            hit.target.description,
            "twin arginine translocation protein; sec-independent protein export [Escherichia coli O157:H7 EDL933] >gi|12518713|gb|AAG59032.1| twin arginine translocation protein; sec-independent protein export [Escherichia coli O157:H7 EDL933] >gi|25302684|pir||D86071 hypothetical protein tatA [imported] - Escherichia coli  (strain O157:H7, substrain EDL933)",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=103)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 462.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 182.57)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.72609e-45)
        self.assertEqual(hsp.annotations["identity"], 102)
        self.assertEqual(hsp.annotations["positive"], 102)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0, 103],
                          [  0, 103]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 103))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MRLCLIIIYHRGTCMGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFK...EQV')",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MRLCLIIIYHRXTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFK...EQV')",
        )
        self.assertEqual(hsp.target.id, "gi|15804428|ref|NP_290468.1|")
        self.assertEqual(hsp.target.name, "NP_290468")
        self.assertEqual(
            hsp.target.description,
            "twin arginine translocation protein; sec-independent protein export [Escherichia coli O157:H7 EDL933] >gi|12518713|gb|AAG59032.1| twin arginine translocation protein; sec-independent protein export [Escherichia coli O157:H7 EDL933] >gi|25302684|pir||D86071 hypothetical protein tatA [imported] - Escherichia coli  (strain O157:H7, substrain EDL933)",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MRLCLIIIYHR TCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|158044         0 MRLCLIIIYHRXTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDD
                  0 |||||||||||.||||||||||...........|||||||||||||||||||||||||||
gi|491764         0 MRLCLIIIYHRGTCMGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDD

gi|158044        60 EPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV 103
                 60 ||||||||||||||||||||||||||||||||||||||||||| 103
gi|491764        60 EPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV 103
""",
        )
        hit = hits[2]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|74314349|ref|YP_312768.1|")
        self.assertEqual(hit.target.name, "YP_312768")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein SSO_4011 [Shigella sonnei Ss046] >gi|73857826|gb|AAZ90533.1| conserved hypothetical protein [Shigella sonnei Ss046]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=103)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 462.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 182.57)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.72609e-45)
        self.assertEqual(hsp.annotations["identity"], 102)
        self.assertEqual(hsp.annotations["positive"], 102)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0, 103],
                          [  0, 103]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 103))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MRLCLIIIYHRGTCMGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFK...EQV')",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MRPCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFK...EQV')",
        )
        self.assertEqual(hsp.target.id, "gi|74314349|ref|YP_312768.1|")
        self.assertEqual(hsp.target.name, "YP_312768")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein SSO_4011 [Shigella sonnei Ss046] >gi|73857826|gb|AAZ90533.1| conserved hypothetical protein [Shigella sonnei Ss046]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MR CLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|743143         0 MRPCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDD
                  0 ||.|||||||||||||||||||...........|||||||||||||||||||||||||||
gi|491764         0 MRLCLIIIYHRGTCMGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDD

gi|743143        60 EPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV 103
                 60 ||||||||||||||||||||||||||||||||||||||||||| 103
gi|491764        60 EPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV 103
""",
        )
        hit = hits[3]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|75256240|ref|ZP_00727918.1|")
        self.assertEqual(hit.target.name, "ZP_00727918")
        self.assertEqual(
            hit.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Escherichia coli E22] >gi|75238286|ref|ZP_00722285.1| COG1826: Sec-independent protein secretion pathway components [Escherichia coli E110019] >gi|75237629|ref|ZP_00721649.1| COG1826: Sec-independent protein secretion pathway components [Escherichia coli F11] >gi|75230797|ref|ZP_00717260.1| COG1826: Sec-independent protein secretion pathway components [Escherichia coli B7A] >gi|75209537|ref|ZP_00709758.1| COG1826: Sec-independent protein secretion pathway components [Escherichia coli B171] >gi|75197719|ref|ZP_00707789.1| COG1826: Sec-independent protein secretion pathway components [Escherichia coli HS] >gi|75188578|ref|ZP_00701845.1| COG1826: Sec-independent protein secretion pathway components [Escherichia coli E24377A] >gi|75177802|ref|ZP_00697867.1| COG1826: Sec-independent protein secretion pathway components [Shigella boydii BS512] >gi|3123496|emb|CAA06724.1| TatA protein [Escherichia coli] >gi|13364242|dbj|BAB38189.1| Sec-independent protein translocase [Escherichia coli O157:H7] >gi|60415956|sp|P69428|TATA_ECOLI Sec-independent protein translocase protein tatA >gi|60415959|sp|P69431|TATA_SHIFL Sec-independent protein translocase protein tatA >gi|60415958|sp|P69430|TATA_ECO57 Sec-independent protein translocase protein tatA >gi|60415957|sp|P69429|TATA_ECOL6 Sec-independent protein translocase protein tatA >gi|15834020|ref|NP_312793.1| Sec-independent protein translocase [Escherichia coli O157:H7]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=89)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 390.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 154.836)
        self.assertAlmostEqual(hsp.annotations["evalue"], 6.0872e-37)
        self.assertEqual(hsp.annotations["identity"], 89)
        self.assertEqual(hsp.annotations["positive"], 89)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0,  89],
                          [ 14, 103]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 89))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...EQV'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...EQV')",
        )
        self.assertEqual(hsp.target.id, "gi|75256240|ref|ZP_00727918.1|")
        self.assertEqual(hsp.target.name, "ZP_00727918")
        self.assertEqual(
            hsp.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Escherichia coli E22] >gi|75238286|ref|ZP_00722285.1| COG1826: Sec-independent protein secretion pathway components [Escherichia coli E110019] >gi|75237629|ref|ZP_00721649.1| COG1826: Sec-independent protein secretion pathway components [Escherichia coli F11] >gi|75230797|ref|ZP_00717260.1| COG1826: Sec-independent protein secretion pathway components [Escherichia coli B7A] >gi|75209537|ref|ZP_00709758.1| COG1826: Sec-independent protein secretion pathway components [Escherichia coli B171] >gi|75197719|ref|ZP_00707789.1| COG1826: Sec-independent protein secretion pathway components [Escherichia coli HS] >gi|75188578|ref|ZP_00701845.1| COG1826: Sec-independent protein secretion pathway components [Escherichia coli E24377A] >gi|75177802|ref|ZP_00697867.1| COG1826: Sec-independent protein secretion pathway components [Shigella boydii BS512] >gi|3123496|emb|CAA06724.1| TatA protein [Escherichia coli] >gi|13364242|dbj|BAB38189.1| Sec-independent protein translocase [Escherichia coli O157:H7] >gi|60415956|sp|P69428|TATA_ECOLI Sec-independent protein translocase protein tatA >gi|60415959|sp|P69431|TATA_SHIFL Sec-independent protein translocase protein tatA >gi|60415958|sp|P69430|TATA_ECO57 Sec-independent protein translocase protein tatA >gi|60415957|sp|P69429|TATA_ECOL6 Sec-independent protein translocase protein tatA >gi|15834020|ref|NP_312793.1| Sec-independent protein translocase [Escherichia coli O157:H7]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|752562         0 MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT
                  0 ||||||||...........|||||||||||||||||||||||||||||||||||||||||
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|752562        60 AKTIADKQADTNQEQAKTEDAKRHDKEQV  89
                 60 |||||||||||||||||||||||||||||  89
gi|491764        74 AKTIADKQADTNQEQAKTEDAKRHDKEQV 103
""",
        )
        hit = hits[4]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|148236|gb|AAA67633.1|")
        self.assertEqual(hit.target.name, "AAA67633")
        self.assertEqual(hit.target.description, "o261 [Escherichia coli]")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=261)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 312.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 124.79)
        self.assertAlmostEqual(hsp.annotations["evalue"], 6.74582e-28)
        self.assertEqual(hsp.annotations["identity"], 62)
        self.assertEqual(hsp.annotations["positive"], 62)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[18, 80],
                          [33, 95]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 62))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQ...EDA'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({18: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQ...EDA'}, length=261)",
        )
        self.assertEqual(hsp.target.id, "gi|148236|gb|AAA67633.1|")
        self.assertEqual(hsp.target.name, "AAA67633")
        self.assertEqual(hsp.target.description, "o261 [Escherichia coli]")
        self.assertEqual(
            hsp.annotations["midline"],
            "FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDA",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|148236        18 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTE
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTE

gi|148236        78 DA 80
                 60 || 62
gi|491764        93 DA 95
""",
        )
        hit = hits[5]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|29143650|ref|NP_806992.1|")
        self.assertEqual(hit.target.name, "NP_806992")
        self.assertEqual(
            hit.target.description,
            "sec-independent protein translocase protein [Salmonella enterica subsp. enterica serovar Typhi Ty2] >gi|16504464|emb|CAD07919.1| sec-independent protein translocase protein [Salmonella enterica subsp. enterica serovar Typhi] >gi|16422538|gb|AAL22817.1| component of Sec-independent protein secretion pathway [Salmonella typhimurium LT2] >gi|62182441|ref|YP_218858.1| component of Sec-independent protein secretion pathway [Salmonella enterica subsp. enterica serovar Choleraesuis str. SC-B67] >gi|16762161|ref|NP_457778.1| sec-independent protein translocase protein [Salmonella enterica subsp. enterica serovar Typhi str. CT18] >gi|56415827|ref|YP_152902.1| sec-independent protein translocase protein [Salmonella enterica subsp. enterica serovar Paratyphi A str. ATCC 9150] >gi|29139285|gb|AAO70852.1| sec-independent protein translocase protein [Salmonella enterica subsp. enterica serovar Typhi Ty2] >gi|56130084|gb|AAV79590.1| sec-independent protein translocase protein [Salmonella enterica subsp. enterica serovar Paratyphi A str. ATCC 9150] >gi|62130074|gb|AAX67777.1| component of Sec-independent protein secretion pathway [Salmonella enterica subsp. enterica serovar Choleraesuis str. SC-B67] >gi|6960229|gb|AAF33419.1| 83% identity with amino acids 10-79 of E. coli hypothetical protein (YIGT) (SP:P27856) [Salmonella typhimurium LT2] >gi|68064449|sp|P0A2H3|TATA_SALTY Sec-independent protein translocase protein tatA >gi|68064450|sp|P0A2H4|TATA_SALTI Sec-independent protein translocase protein tatA >gi|25302682|pir||AI0915 sec-independent protein translocase protein [imported] - Salmonella enterica subsp. enterica serovar Typhi (strain CT18) >gi|16767243|ref|NP_462858.1| Sec-independent protein secretion pathway component [Salmonella typhimurium LT2]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=84)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 305.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 122.094)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.37251e-27)
        self.assertEqual(hsp.annotations["identity"], 75)
        self.assertEqual(hsp.annotations["positive"], 79)
        self.assertEqual(hsp.annotations["gaps"], 5)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0,  69,  69,  84],
                          [ 14,  83,  88, 103]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 89))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...EQV'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MGGISIWQLLIVAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDDAKQDKTS...EQV')",
        )
        self.assertEqual(hsp.target.id, "gi|29143650|ref|NP_806992.1|")
        self.assertEqual(hsp.target.name, "NP_806992")
        self.assertEqual(
            hsp.target.description,
            "sec-independent protein translocase protein [Salmonella enterica subsp. enterica serovar Typhi Ty2] >gi|16504464|emb|CAD07919.1| sec-independent protein translocase protein [Salmonella enterica subsp. enterica serovar Typhi] >gi|16422538|gb|AAL22817.1| component of Sec-independent protein secretion pathway [Salmonella typhimurium LT2] >gi|62182441|ref|YP_218858.1| component of Sec-independent protein secretion pathway [Salmonella enterica subsp. enterica serovar Choleraesuis str. SC-B67] >gi|16762161|ref|NP_457778.1| sec-independent protein translocase protein [Salmonella enterica subsp. enterica serovar Typhi str. CT18] >gi|56415827|ref|YP_152902.1| sec-independent protein translocase protein [Salmonella enterica subsp. enterica serovar Paratyphi A str. ATCC 9150] >gi|29139285|gb|AAO70852.1| sec-independent protein translocase protein [Salmonella enterica subsp. enterica serovar Typhi Ty2] >gi|56130084|gb|AAV79590.1| sec-independent protein translocase protein [Salmonella enterica subsp. enterica serovar Paratyphi A str. ATCC 9150] >gi|62130074|gb|AAX67777.1| component of Sec-independent protein secretion pathway [Salmonella enterica subsp. enterica serovar Choleraesuis str. SC-B67] >gi|6960229|gb|AAF33419.1| 83% identity with amino acids 10-79 of E. coli hypothetical protein (YIGT) (SP:P27856) [Salmonella typhimurium LT2] >gi|68064449|sp|P0A2H3|TATA_SALTY Sec-independent protein translocase protein tatA >gi|68064450|sp|P0A2H4|TATA_SALTI Sec-independent protein translocase protein tatA >gi|25302682|pir||AI0915 sec-independent protein translocase protein [imported] - Salmonella enterica subsp. enterica serovar Typhi (strain CT18) >gi|16767243|ref|NP_462858.1| Sec-independent protein secretion pathway component [Salmonella typhimurium LT2]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGGISIWQLLI+AVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDD+ KQDKTSQDADFTAK+IADKQ      +AK EDAK  DKEQV",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|291436         0 MGGISIWQLLIVAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDDAKQDKTSQDADFT
                  0 ||||||||...........|||||||||||||||||||||||||||..||||||||||||
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|291436        60 AKSIADKQG-----EAKKEDAKSQDKEQV  84
                 60 ||.|||||.-----.||.||||..|||||  89
gi|491764        74 AKTIADKQADTNQEQAKTEDAKRHDKEQV 103
""",
        )
        hit = hits[6]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|49609685|emb|CAG73118.1|")
        self.assertEqual(hit.target.name, "CAG73118")
        self.assertEqual(
            hit.target.description,
            "sec-independent protein translocase [Erwinia carotovora subsp. atroseptica SCRI1043] >gi|50119159|ref|YP_048326.1| sec-independent protein translocase [Erwinia carotovora subsp. atroseptica SCRI1043]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=86)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 219.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 88.9669)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.10205e-17)
        self.assertEqual(hsp.annotations["identity"], 59)
        self.assertEqual(hsp.annotations["positive"], 67)
        self.assertEqual(hsp.annotations["gaps"], 7)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0,  48,  50,  72,  72,  86],
                          [ 14,  62,  62,  84,  89, 103]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 91))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...EQV'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MGGISLWNLLIIAVIVILLFGTNKLRTLGSDLGASIKGFKKAMGDDQPSTNADK...EQV')",
        )
        self.assertEqual(hsp.target.id, "gi|49609685|emb|CAG73118.1|")
        self.assertEqual(hsp.target.name, "CAG73118")
        self.assertEqual(
            hsp.target.description,
            "sec-independent protein translocase [Erwinia carotovora subsp. atroseptica SCRI1043] >gi|50119159|ref|YP_048326.1| sec-independent protein translocase [Erwinia carotovora subsp. atroseptica SCRI1043]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGGIS+W LLIIAVIV+LLFGT KL ++GSDLGASIKGFKKAM DD+P    DK   DADF+ K+IAD Q+D     AK  D K   KEQV",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|496096         0 MGGISLWNLLIIAVIVILLFGTNKLRTLGSDLGASIKGFKKAMGDDQPSTNADKAQPDAD
                  0 |||||.|............|||.||...|||||||||||||||.||.|--..||...|||
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEP--KQDKTSQDAD

gi|496096        60 FSTKSIADNQSD-----AKQGDTKSQHKEQV  86
                 60 |..|.|||.|.|-----||..|.|...||||  91
gi|491764        72 FTAKTIADKQADTNQEQAKTEDAKRHDKEQV 103
""",
        )
        hit = hits[7]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|37528238|ref|NP_931583.1|")
        self.assertEqual(hit.target.name, "NP_931583")
        self.assertEqual(
            hit.target.description,
            "Sec-independent protein translocase protein [Photorhabdus luminescens subsp. laumondii TTO1] >gi|36787675|emb|CAE16782.1| Sec-independent protein translocase protein [Photorhabdus luminescens subsp. laumondii TTO1]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=86)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 204.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 83.1889)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.25087e-15)
        self.assertEqual(hsp.annotations["identity"], 59)
        self.assertEqual(hsp.annotations["positive"], 71)
        self.assertEqual(hsp.annotations["gaps"], 7)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0,  46,  47,  51,  52,  71,  71,  86],
                          [ 14,  60,  60,  64,  64,  83,  88, 103]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 91))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...EQV'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MGGISIWQLLIIAVIVVLLFGTNKLRTLGSDLGASIKGFKKAIGDDNQPQQAQK...EQV')",
        )
        self.assertEqual(hsp.target.id, "gi|37528238|ref|NP_931583.1|")
        self.assertEqual(hsp.target.name, "NP_931583")
        self.assertEqual(
            hsp.target.description,
            "Sec-independent protein translocase protein [Photorhabdus luminescens subsp. laumondii TTO1] >gi|36787675|emb|CAE16782.1| Sec-independent protein translocase protein [Photorhabdus luminescens subsp. laumondii TTO1]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGGISIWQLLIIAVIVVLLFGT KL ++GSDLGASIKGFKKA+ DD +P+Q  KTS DADF  K I +KQ+      A++E ++  +KEQV",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|375282         0 MGGISIWQLLIIAVIVVLLFGTNKLRTLGSDLGASIKGFKKAIGDDNQPQQAQKTSSDAD
                  0 ||||||||...........|||.||...||||||||||||||..||-.|.|-.|||.|||
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDD-EPKQ-DKTSQDAD

gi|375282        60 FETKNITEKQS-----VAQSETSESKNKEQV  86
                 60 |..|.|..||.-----.|..|......||||  91
gi|491764        72 FTAKTIADKQADTNQEQAKTEDAKRHDKEQV 103
""",
        )
        hit = hits[8]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|59710656|ref|YP_203432.1|")
        self.assertEqual(hit.target.name, "YP_203432")
        self.assertEqual(
            hit.target.description,
            "Sec-independent protein translocase protein TatA [Vibrio fischeri ES114] >gi|71657992|sp|Q5E8V2|TATA_VIBF1 Sec-independent protein translocase protein tatA/E homolog >gi|59478757|gb|AAW84544.1| Sec-independent protein translocase protein TatA [Vibrio fischeri ES114]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=82)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 201.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 82.0333)
        self.assertAlmostEqual(hsp.annotations["evalue"], 5.01441e-15)
        self.assertEqual(hsp.annotations["identity"], 52)
        self.assertEqual(hsp.annotations["positive"], 70)
        self.assertEqual(hsp.annotations["gaps"], 7)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0,  51,  51,  74,  74,  81],
                          [ 14,  65,  66,  89,  95, 102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 88))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...KEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGGISIWQLLIIAVIIVLLFGTKKLRGVGSDLGSAVKGFKKAISEDEPAKEAKK...KEQ'}, length=82)",
        )
        self.assertEqual(hsp.target.id, "gi|59710656|ref|YP_203432.1|")
        self.assertEqual(hsp.target.name, "YP_203432")
        self.assertEqual(
            hsp.target.description,
            "Sec-independent protein translocase protein TatA [Vibrio fischeri ES114] >gi|71657992|sp|Q5E8V2|TATA_VIBF1 Sec-independent protein translocase protein tatA/E homolog >gi|59478757|gb|AAW84544.1| Sec-independent protein translocase protein TatA [Vibrio fischeri ES114]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGGISIWQLLIIAVI+VLLFGTKKL  +GSDLG+++KGFKKA+S+DEP ++   +DADF  + +  K+A+T ++Q      K++DKEQ",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|597106         0 MGGISIWQLLIIAVIIVLLFGTKKLRGVGSDLGSAVKGFKKAISEDEPAKE-AKKDADFV
                  0 ||||||||...........||||||...|||||...||||||.|.|||...-...||||.
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|597106        59 PQNLEKKEAETVEKQ------KQNDKEQ  81
                 60 ......|.|.|...|------|..||||  88
gi|491764        74 AKTIADKQADTNQEQAKTEDAKRHDKEQ 102
""",
        )
        hit = hits[9]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|54307340|ref|YP_128360.1|")
        self.assertEqual(hit.target.name, "YP_128360")
        self.assertEqual(
            hit.target.description,
            "putative TatA protein [Photobacterium profundum SS9] >gi|46911760|emb|CAG18558.1| putative TatA protein [Photobacterium profundum SS9]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=87)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 191.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 78.1814)
        self.assertAlmostEqual(hsp.annotations["evalue"], 7.2408e-14)
        self.assertEqual(hsp.annotations["identity"], 51)
        self.assertEqual(hsp.annotations["positive"], 67)
        self.assertEqual(hsp.annotations["gaps"], 6)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0,  47,  49,  81,  81,  87],
                          [ 14,  61,  61,  93,  97, 103]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 91))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...EQV'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MGGISIWQLLIIALIIVLLFGTKKLRSLGGDLGSAVKGFKKAIGDEELTVKKDN...EQV')",
        )
        self.assertEqual(hsp.target.id, "gi|54307340|ref|YP_128360.1|")
        self.assertEqual(hsp.target.name, "YP_128360")
        self.assertEqual(
            hsp.target.description,
            "putative TatA protein [Photobacterium profundum SS9] >gi|46911760|emb|CAG18558.1| putative TatA protein [Photobacterium profundum SS9]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGGISIWQLLIIA+I+VLLFGTKKL S+G DLG+++KGFKKA+ D+E   K+D T  DADF  KT++ ++  +     K++     DKEQV",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|543073         0 MGGISIWQLLIIALIIVLLFGTKKLRSLGGDLGSAVKGFKKAIGDEELTVKKDNTEADAD
                  0 ||||||||...........||||||.|.|.|||...||||||..|.|--.|.|.|..|||
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDE--PKQDKTSQDAD

gi|543073        60 FEQKTLSKEEQQSEDPVQKSQ----KDKEQV  87
                 60 |..||.............|..----.|||||  91
gi|491764        72 FTAKTIADKQADTNQEQAKTEDAKRHDKEQV 103
""",
        )
        hit = hits[10]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|45437890|gb|AAS63439.1|")
        self.assertEqual(hit.target.name, "AAS63439")
        self.assertEqual(
            hit.target.description,
            "Sec-independent protein translocase protein TatA [Yersinia pestis biovar Medievalis str. 91001] >gi|51594613|ref|YP_068804.1| Sec-independent protein translocase protein TatA [Yersinia pseudotuberculosis IP 32953] >gi|22124367|ref|NP_667790.1| hypothetical protein y0452 [Yersinia pestis KIM] >gi|77633885|ref|ZP_00795999.1| COG1826: Sec-independent protein secretion pathway components [Yersinia pestis Angola] >gi|77630676|ref|ZP_00793262.1| COG1826: Sec-independent protein secretion pathway components [Yersinia pseudotuberculosis IP 31758] >gi|51587895|emb|CAH19498.1| Sec-independent protein translocase protein TatA [Yersinia pseudotuberculosis IP 32953] >gi|15981691|emb|CAC93245.1| Sec-independent protein translocase protein TatA [Yersinia pestis CO92] >gi|21957146|gb|AAM84041.1| hypothetical protein [Yersinia pestis KIM] >gi|45443023|ref|NP_994562.1| Sec-independent protein translocase protein TatA [Yersinia pestis biovar Medievalis str. 91001] >gi|16123912|ref|NP_407225.1| Sec-independent protein translocase protein TatA [Yersinia pestis CO92] >gi|25302681|pir||AI0459 Sec-independent protein translocase protein TatA [imported] - Yersinia pestis (strain CO92) >gi|24212497|sp|Q8ZAM2|TATA_YERPE Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=88)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 188.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 77.0258)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.61308e-13)
        self.assertEqual(hsp.annotations["identity"], 59)
        self.assertEqual(hsp.annotations["positive"], 66)
        self.assertEqual(hsp.annotations["gaps"], 9)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0,  47,  51,  63,  63,  71,  71,  87],
                          [ 14,  61,  61,  73,  74,  82,  86, 102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 92))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...KEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGSIGWAQLLIIAVIVVLLFGTNKLRTLGSDLGASIKGFKKAMGDDSQTPPTNV...KEQ'}, length=88)",
        )
        self.assertEqual(hsp.target.id, "gi|45437890|gb|AAS63439.1|")
        self.assertEqual(hsp.target.name, "AAS63439")
        self.assertEqual(
            hsp.target.description,
            "Sec-independent protein translocase protein TatA [Yersinia pestis biovar Medievalis str. 91001] >gi|51594613|ref|YP_068804.1| Sec-independent protein translocase protein TatA [Yersinia pseudotuberculosis IP 32953] >gi|22124367|ref|NP_667790.1| hypothetical protein y0452 [Yersinia pestis KIM] >gi|77633885|ref|ZP_00795999.1| COG1826: Sec-independent protein secretion pathway components [Yersinia pestis Angola] >gi|77630676|ref|ZP_00793262.1| COG1826: Sec-independent protein secretion pathway components [Yersinia pseudotuberculosis IP 31758] >gi|51587895|emb|CAH19498.1| Sec-independent protein translocase protein TatA [Yersinia pseudotuberculosis IP 32953] >gi|15981691|emb|CAC93245.1| Sec-independent protein translocase protein TatA [Yersinia pestis CO92] >gi|21957146|gb|AAM84041.1| hypothetical protein [Yersinia pestis KIM] >gi|45443023|ref|NP_994562.1| Sec-independent protein translocase protein TatA [Yersinia pestis biovar Medievalis str. 91001] >gi|16123912|ref|NP_407225.1| Sec-independent protein translocase protein TatA [Yersinia pestis CO92] >gi|25302681|pir||AI0459 Sec-independent protein translocase protein TatA [imported] - Yersinia pestis (strain CO92) >gi|24212497|sp|Q8ZAM2|TATA_YERPE Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG I   QLLIIAVIVVLLFGT KL ++GSDLGASIKGFKKAM DD        DKTS DADF AK+I +KQ    Q  AK E++K H+KEQ",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|454378         0 MGSIGWAQLLIIAVIVVLLFGTNKLRTLGSDLGASIKGFKKAMGDDSQTPPTNVDKTSND
                  0 ||.|...|...........|||.||...|||||||||||||||.||.----...||||.|
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDE----PKQDKTSQD

gi|454378        60 ADF-AKSITEKQ----QPVAKAEESKSHEKEQ  87
                 60 |||-||.|..||----|..||.|..|.|.|||  92
gi|491764        70 ADFTAKTIADKQADTNQEQAKTEDAKRHDKEQ 102
""",
        )
        hit = hits[11]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|75856473|ref|ZP_00764101.1|")
        self.assertEqual(hit.target.name, "ZP_00764101")
        self.assertEqual(
            hit.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Vibrio sp. Ex25]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=81)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 187.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 76.6406)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.10675e-13)
        self.assertEqual(hsp.annotations["identity"], 54)
        self.assertEqual(hsp.annotations["positive"], 66)
        self.assertEqual(hsp.annotations["gaps"], 8)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0,  49,  49,  66,  66,  76,  76,  80],
                          [ 14,  63,  65,  82,  86,  96,  98, 102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 88))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...KEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGGISVWQLLIIAVIVVLLFGTKKLRGIGGDLGGAVKGFKKAMSEDEPAKNDKD...KEQ'}, length=81)",
        )
        self.assertEqual(hsp.target.id, "gi|75856473|ref|ZP_00764101.1|")
        self.assertEqual(hsp.target.name, "ZP_00764101")
        self.assertEqual(
            hsp.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Vibrio sp. Ex25]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGGIS+WQLLIIAVIVVLLFGTKKL  IG DLG ++KGFKKAMS+DEP   K  +DADF  K++ ++Q    +++A  E  K  DKEQ",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|758564         0 MGGISVWQLLIIAVIVVLLFGTKKLRGIGGDLGGAVKGFKKAMSEDEPA--KNDKDADFE
                  0 |||||.||...........||||||..||.|||...||||||||.|||.--|...||||.
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|758564        58 PKSLEEQQ----KKEAAPETKK--DKEQ  80
                 60 .|.....|----...|..|..|--||||  88
gi|491764        74 AKTIADKQADTNQEQAKTEDAKRHDKEQ 102
""",
        )
        hit = hits[12]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|75829371|ref|ZP_00758676.1|")
        self.assertEqual(hit.target.name, "ZP_00758676")
        self.assertEqual(
            hit.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Vibrio cholerae MO10] >gi|75828068|ref|ZP_00757503.1| COG1826: Sec-independent protein secretion pathway components [Vibrio cholerae O395] >gi|75823218|ref|ZP_00752737.1| COG1826: Sec-independent protein secretion pathway components [Vibrio cholerae RC385] >gi|75816411|ref|ZP_00746882.1| COG1826: Sec-independent protein secretion pathway components [Vibrio cholerae V52] >gi|9654483|gb|AAF93264.1| tatA protein [Vibrio cholerae O1 biovar eltor str. N16961] >gi|15640118|ref|NP_229745.1| tatA protein [Vibrio cholerae O1 biovar eltor str. N16961] >gi|11356194|pir||G82366 tatA protein VC0086 [imported] - Vibrio cholerae (strain N16961 serogroup O1) >gi|9979012|sp|P57051|TATA_VIBCH Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=82)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 186.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 76.2554)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.7515e-13)
        self.assertEqual(hsp.annotations["identity"], 55)
        self.assertEqual(hsp.annotations["positive"], 63)
        self.assertEqual(hsp.annotations["gaps"], 11)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0,  62,  62,  73,  75,  81],
                          [ 14,  76,  85,  96,  96, 102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 90))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...KEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGGISIWQLLIIAVIVVLLFGTKKLRGIGSDLGSAVKGFKKAMSEEESNSAANQ...KEQ'}, length=82)",
        )
        self.assertEqual(hsp.target.id, "gi|75829371|ref|ZP_00758676.1|")
        self.assertEqual(hsp.target.name, "ZP_00758676")
        self.assertEqual(
            hsp.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Vibrio cholerae MO10] >gi|75828068|ref|ZP_00757503.1| COG1826: Sec-independent protein secretion pathway components [Vibrio cholerae O395] >gi|75823218|ref|ZP_00752737.1| COG1826: Sec-independent protein secretion pathway components [Vibrio cholerae RC385] >gi|75816411|ref|ZP_00746882.1| COG1826: Sec-independent protein secretion pathway components [Vibrio cholerae V52] >gi|9654483|gb|AAF93264.1| tatA protein [Vibrio cholerae O1 biovar eltor str. N16961] >gi|15640118|ref|NP_229745.1| tatA protein [Vibrio cholerae O1 biovar eltor str. N16961] >gi|11356194|pir||G82366 tatA protein VC0086 [imported] - Vibrio cholerae (strain N16961 serogroup O1) >gi|9979012|sp|P57051|TATA_VIBCH Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGGISIWQLLIIAVIVVLLFGTKKL  IGSDLG+++KGFKKAMS++E       +DADF  K         N EQAKT  +   + DKEQ",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|758293         0 MGGISIWQLLIIAVIVVLLFGTKKLRGIGSDLGSAVKGFKKAMSEEESNSAANQKDADFE
                  0 ||||||||...........||||||..||||||...||||||||..|........||||.
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|758293        60 TK---------NLEQAKTNASAEVKKDKEQ  81
                 60 .|---------|.|||||....--..||||  90
gi|491764        74 AKTIADKQADTNQEQAKTEDAK--RHDKEQ 102
""",
        )
        hit = hits[13]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|75820019|ref|ZP_00750077.1|")
        self.assertEqual(hit.target.name, "ZP_00750077")
        self.assertEqual(
            hit.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Vibrio cholerae V51]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=82)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 185.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 75.8702)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.59357e-13)
        self.assertEqual(hsp.annotations["identity"], 55)
        self.assertEqual(hsp.annotations["positive"], 63)
        self.assertEqual(hsp.annotations["gaps"], 11)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0,  62,  62,  72,  74,  81],
                          [ 14,  76,  85,  95,  95, 102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 90))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...KEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGGISIWQLLIIAVIVVLLFGTKKLRGIGSDLGSAVKGFKKAMSEEESNSAANQ...KEQ'}, length=82)",
        )
        self.assertEqual(hsp.target.id, "gi|75820019|ref|ZP_00750077.1|")
        self.assertEqual(hsp.target.name, "ZP_00750077")
        self.assertEqual(
            hsp.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Vibrio cholerae V51]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGGISIWQLLIIAVIVVLLFGTKKL  IGSDLG+++KGFKKAMS++E       +DADF  K         N EQAKT  +   + DKEQ",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|758200         0 MGGISIWQLLIIAVIVVLLFGTKKLRGIGSDLGSAVKGFKKAMSEEESNSAANQKDADFE
                  0 ||||||||...........||||||..||||||...||||||||..|........||||.
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|758200        60 TK---------NLEQAKTNASVEVKKDKEQ  81
                 60 .|---------|.|||||...--...||||  90
gi|491764        74 AKTIADKQADTNQEQAKTEDA--KRHDKEQ 102
""",
        )
        hit = hits[14]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|28896872|ref|NP_796477.1|")
        self.assertEqual(hit.target.name, "NP_796477")
        self.assertEqual(
            hit.target.description,
            "TatA protein [Vibrio parahaemolyticus RIMD 2210633] >gi|28805080|dbj|BAC58361.1| TatA protein [Vibrio parahaemolyticus RIMD 2210633] >gi|32130181|sp|Q87TH1|TATA_VIBPA Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=81)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 178.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 73.1738)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.32928e-12)
        self.assertEqual(hsp.annotations["identity"], 53)
        self.assertEqual(hsp.annotations["positive"], 65)
        self.assertEqual(hsp.annotations["gaps"], 8)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0,  47,  47,  66,  66,  76,  76,  80],
                          [ 14,  61,  63,  82,  86,  96,  98, 102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 88))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...KEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGGISVWQLLIIAVIVVLLFGTKKLRGIGGDLGSAVKGFKKAMSDEDSAKNEKD...KEQ'}, length=81)",
        )
        self.assertEqual(hsp.target.id, "gi|28896872|ref|NP_796477.1|")
        self.assertEqual(hsp.target.name, "NP_796477")
        self.assertEqual(
            hsp.target.description,
            "TatA protein [Vibrio parahaemolyticus RIMD 2210633] >gi|28805080|dbj|BAC58361.1| TatA protein [Vibrio parahaemolyticus RIMD 2210633] >gi|32130181|sp|Q87TH1|TATA_VIBPA Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGGIS+WQLLIIAVIVVLLFGTKKL  IG DLG+++KGFKKAMSD++    K  +DADF  K++  +Q    Q++A  E  K  DKEQ",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|288968         0 MGGISVWQLLIIAVIVVLLFGTKKLRGIGGDLGSAVKGFKKAMSDED--SAKNEKDADFE
                  0 |||||.||...........||||||..||.|||...|||||||||..--..|...||||.
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|288968        58 PKSLEKQQ----QKEAAPETKK--DKEQ  80
                 60 .|.....|----|..|..|..|--||||  88
gi|491764        74 AKTIADKQADTNQEQAKTEDAKRHDKEQ 102
""",
        )
        hit = hits[15]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|27364353|ref|NP_759881.1|")
        self.assertEqual(hit.target.name, "NP_759881")
        self.assertEqual(
            hit.target.description,
            "Sec-independent protein secretion pathway component TatA [Vibrio vulnificus CMCP6] >gi|27360472|gb|AAO09408.1| Sec-independent protein secretion pathway component TatA [Vibrio vulnificus CMCP6] >gi|32130186|sp|Q8DDQ2|TATA_VIBVU Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=78)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 176.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 72.4034)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.97316e-12)
        self.assertEqual(hsp.annotations["identity"], 54)
        self.assertEqual(hsp.annotations["positive"], 65)
        self.assertEqual(hsp.annotations["gaps"], 11)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0,  50,  50,  60,  60,  73,  73,  77],
                          [ 14,  64,  68,  78,  83,  96,  98, 102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 88))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...KEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGGISIWQLLIIAVIVVLLFGTKKLRGIGSDLGGAIKGFKKAMNEEESEKKDAD...KEQ'}, length=78)",
        )
        self.assertEqual(hsp.target.id, "gi|27364353|ref|NP_759881.1|")
        self.assertEqual(hsp.target.name, "NP_759881")
        self.assertEqual(
            hsp.target.description,
            "Sec-independent protein secretion pathway component TatA [Vibrio vulnificus CMCP6] >gi|27360472|gb|AAO09408.1| Sec-independent protein secretion pathway component TatA [Vibrio vulnificus CMCP6] >gi|32130186|sp|Q8DDQ2|TATA_VIBVU Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGGISIWQLLIIAVIVVLLFGTKKL  IGSDLG +IKGFKKAM+++E ++    +DADF  K++     +   +QA TE  K  DKEQ",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|273643         0 MGGISIWQLLIIAVIVVLLFGTKKLRGIGSDLGGAIKGFKKAMNEEESEK----KDADFE
                  0 ||||||||...........||||||..||||||..||||||||...|...----.||||.
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|273643        56 PKSL-----EQQSKQAATESKK--DKEQ  77
                 60 .|..-----.....||.||..|--||||  88
gi|491764        74 AKTIADKQADTNQEQAKTEDAKRHDKEQ 102
""",
        )
        hit = hits[16]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|37678364|ref|NP_932973.1|")
        self.assertEqual(hit.target.name, "NP_932973")
        self.assertEqual(
            hit.target.description,
            "Sec-independent protein secretion pathway component [Vibrio vulnificus YJ016] >gi|37197103|dbj|BAC92944.1| Sec-independent protein secretion pathway component [Vibrio vulnificus YJ016] >gi|71153586|sp|Q7MQ30|TATA_VIBVY Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=78)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 176.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 72.4034)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.97316e-12)
        self.assertEqual(hsp.annotations["identity"], 54)
        self.assertEqual(hsp.annotations["positive"], 65)
        self.assertEqual(hsp.annotations["gaps"], 11)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0,  50,  50,  60,  60,  73,  73,  77],
                          [ 14,  64,  68,  78,  83,  96,  98, 102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 88))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...KEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGGISIWQLLIIAVIVVLLFGTKKLRGIGSDLGGAIKGFKKAMNEEESEKKDAD...KEQ'}, length=78)",
        )
        self.assertEqual(hsp.target.id, "gi|37678364|ref|NP_932973.1|")
        self.assertEqual(hsp.target.name, "NP_932973")
        self.assertEqual(
            hsp.target.description,
            "Sec-independent protein secretion pathway component [Vibrio vulnificus YJ016] >gi|37197103|dbj|BAC92944.1| Sec-independent protein secretion pathway component [Vibrio vulnificus YJ016] >gi|71153586|sp|Q7MQ30|TATA_VIBVY Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGGISIWQLLIIAVIVVLLFGTKKL  IGSDLG +IKGFKKAM+++E ++    +DADF  K++     +   +QA TE  K  DKEQ",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|376783         0 MGGISIWQLLIIAVIVVLLFGTKKLRGIGSDLGGAIKGFKKAMNEEESEK----KDADFE
                  0 ||||||||...........||||||..||||||..||||||||...|...----.||||.
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|376783        56 PKSL-----EQQNKQAATESKK--DKEQ  77
                 60 .|..-----.....||.||..|--||||  88
gi|491764        74 AKTIADKQADTNQEQAKTEDAKRHDKEQ 102
""",
        )
        hit = hits[17]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|71277787|ref|YP_266931.1|")
        self.assertEqual(hit.target.name, "YP_266931")
        self.assertEqual(
            hit.target.description,
            "Sec-independent protein translocase protein TatA [Colwellia psychrerythraea 34H] >gi|71143527|gb|AAZ24000.1| Sec-independent protein translocase protein TatA [Colwellia psychrerythraea 34H]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=85)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 176.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 72.4034)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.97316e-12)
        self.assertEqual(hsp.annotations["identity"], 48)
        self.assertEqual(hsp.annotations["positive"], 64)
        self.assertEqual(hsp.annotations["gaps"], 4)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0,  46,  46,  52,  52,  84],
                          [ 14,  60,  62,  68,  70, 102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 88))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...KEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGGIGIWQLVIVAVIVVLLFGTKKLRNIGGDLGSAIKGFKSAIGEDKEQKNSAE...SDQ'}, length=85)",
        )
        self.assertEqual(hsp.target.id, "gi|71277787|ref|YP_266931.1|")
        self.assertEqual(hsp.target.name, "YP_266931")
        self.assertEqual(
            hsp.target.description,
            "Sec-independent protein translocase protein TatA [Colwellia psychrerythraea 34H] >gi|71143527|gb|AAZ24000.1| Sec-independent protein translocase protein TatA [Colwellia psychrerythraea 34H]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGGI IWQL+I+AVIVVLLFGTKKL +IG DLG++IKGFK A+ +D  K+ K S  A+ T+ T+AD    T +E  KT ++K  + +Q",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|712777         0 MGGIGIWQLVIVAVIVVLLFGTKKLRNIGGDLGSAIKGFKSAIGED--KEQKNS--AEKT
                  0 ||||.|||...........||||||..||.|||..|||||.|...|--|..|.|--|..|
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|712777        56 SDTLADSSKSTTEEVVKTTESKTKESDQ  84
                 60 ..|.||....|..|..||...|.....|  88
gi|491764        74 AKTIADKQADTNQEQAKTEDAKRHDKEQ 102
""",
        )
        hit = hits[18]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|68541995|ref|ZP_00581733.1|")
        self.assertEqual(hit.target.name, "ZP_00581733")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Shewanella baltica OS155] >gi|68520382|gb|EAN43890.1| Twin-arginine translocation protein TatA/E [Shewanella baltica OS155]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=79)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 169.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 69.707)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.57533e-11)
        self.assertEqual(hsp.annotations["identity"], 46)
        self.assertEqual(hsp.annotations["positive"], 60)
        self.assertEqual(hsp.annotations["gaps"], 3)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 63, 66, 79],
                          [14, 77, 77, 90]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 79))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...EQA'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MGGISIWQLLIVALIVVLLFGTKKLRSLGGDLGGAVKGFKNAMSSEEDKKALED...EQA')",
        )
        self.assertEqual(hsp.target.id, "gi|68541995|ref|ZP_00581733.1|")
        self.assertEqual(hsp.target.name, "ZP_00581733")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Shewanella baltica OS155] >gi|68520382|gb|EAN43890.1| Twin-arginine translocation protein TatA/E [Shewanella baltica OS155]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGGISIWQLLI+A+IVVLLFGTKKL S+G DLG ++KGFK AMS +E K+     +A  TA+T     +K+ ++N+EQA",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|685419         0 MGGISIWQLLIVALIVVLLFGTKKLRSLGGDLGGAVKGFKNAMSSEEDKKALEDTEAAKT
                  0 ||||||||...........||||||.|.|.|||...||||.|||..|.|.......|..|
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|685419        60 AQTTQQATEKKPESNKEQA 79
                 60 |.|---...|....|.||| 79
gi|491764        74 AKT---IADKQADTNQEQA 90
""",
        )
        hit = hits[19]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|77813363|ref|ZP_00812641.1|")
        self.assertEqual(hit.target.name, "ZP_00812641")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Shewanella putrefaciens CN-32] >gi|77811811|gb|EAO96181.1| Twin-arginine translocation protein TatA/E [Shewanella putrefaciens CN-32]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=79)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 169.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 69.707)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.57533e-11)
        self.assertEqual(hsp.annotations["identity"], 47)
        self.assertEqual(hsp.annotations["positive"], 60)
        self.assertEqual(hsp.annotations["gaps"], 3)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 63, 66, 79],
                          [14, 77, 77, 90]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 79))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...EQA'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MGGISIWQLLIIALIVVLLFGTKKLRSLGGDLGGAVKGFKNAMSSEEDKKALED...EQA')",
        )
        self.assertEqual(hsp.target.id, "gi|77813363|ref|ZP_00812641.1|")
        self.assertEqual(hsp.target.name, "ZP_00812641")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Shewanella putrefaciens CN-32] >gi|77811811|gb|EAO96181.1| Twin-arginine translocation protein TatA/E [Shewanella putrefaciens CN-32]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGGISIWQLLIIA+IVVLLFGTKKL S+G DLG ++KGFK AMS +E K+     +A  TA+T     +K+ ++N+EQA",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|778133         0 MGGISIWQLLIIALIVVLLFGTKKLRSLGGDLGGAVKGFKNAMSSEEDKKALEDTEAAKT
                  0 ||||||||...........||||||.|.|.|||...||||.|||..|.|.......|..|
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|778133        60 AQTTQQATEKKPESNKEQA 79
                 60 |.|---...|....|.||| 79
gi|491764        74 AKT---IADKQADTNQEQA 90
""",
        )
        hit = hits[20]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|52306607|gb|AAU37107.1|")
        self.assertEqual(hit.target.name, "AAU37107")
        self.assertEqual(
            hit.target.description,
            "TatA protein [Mannheimia succiniciproducens MBEL55E] >gi|52424555|ref|YP_087692.1| TatA protein [Mannheimia succiniciproducens MBEL55E]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=75)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 168.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 69.3218)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.36348e-11)
        self.assertEqual(hsp.annotations["identity"], 49)
        self.assertEqual(hsp.annotations["positive"], 60)
        self.assertEqual(hsp.annotations["gaps"], 14)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0,  49,  49,  57,  57,  70,  70,  74],
                          [ 14,  63,  69,  77,  78,  91,  98, 102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 88))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...KEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGGISIWQLLIIVAIIVLLFGTKKLRTLGTDLGESVKGFKKAMNEDEPKDAEFK...KEQ'}, length=75)",
        )
        self.assertEqual(hsp.target.id, "gi|52306607|gb|AAU37107.1|")
        self.assertEqual(hsp.target.name, "AAU37107")
        self.assertEqual(
            hsp.target.description,
            "TatA protein [Mannheimia succiniciproducens MBEL55E] >gi|52424555|ref|YP_087692.1| TatA protein [Mannheimia succiniciproducens MBEL55E]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGGISIWQLLII  I+VLLFGTKKL ++G+DLG S+KGFKKAM++DEPK      DA+F +    D+ A    E+ K       DKEQ",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|523066         0 MGGISIWQLLIIVAIIVLLFGTKKLRTLGTDLGESVKGFKKAMNEDEPK------DAEFK
                  0 ||||||||...........||||||...|.|||.|.|||||||..||||------||.|.
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|523066        54 SLN-KDESATAGSEKVK-------DKEQ  74
                 60 ...-.|..|....|..|-------||||  88
gi|491764        74 AKTIADKQADTNQEQAKTEDAKRHDKEQ 102
""",
        )
        hit = hits[21]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|75429751|ref|ZP_00732413.1|")
        self.assertEqual(hit.target.name, "ZP_00732413")
        self.assertEqual(
            hit.target.description,
            "sec-independent protein secretion pathway components [Actinobacillus succinogenes 130Z] >gi|74276959|gb|EAO50546.1| sec-independent protein secretion pathway components [Actinobacillus succinogenes 130Z]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=74)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 168.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 69.3218)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.36348e-11)
        self.assertEqual(hsp.annotations["identity"], 48)
        self.assertEqual(hsp.annotations["positive"], 59)
        self.assertEqual(hsp.annotations["gaps"], 7)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 44, 44, 48, 48, 74],
                          [14, 58, 59, 63, 69, 95]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 81))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...EDA'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MGGISIWQLLIIVAIVVLLFGTKKLRTLGSDLGESVKGFKKAMAEEPKDAEFKS...EQA')",
        )
        self.assertEqual(hsp.target.id, "gi|75429751|ref|ZP_00732413.1|")
        self.assertEqual(hsp.target.name, "ZP_00732413")
        self.assertEqual(
            hsp.target.description,
            "sec-independent protein secretion pathway components [Actinobacillus succinogenes 130Z] >gi|74276959|gb|EAO50546.1| sec-independent protein secretion pathway components [Actinobacillus succinogenes 130Z]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGGISIWQLLII  IVVLLFGTKKL ++GSDLG S+KGFKKAM+ +EPK      DA+F +   A+  A T +E+ + E A",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|754297         0 MGGISIWQLLIIVAIVVLLFGTKKLRTLGSDLGESVKGFKKAMA-EEPK------DAEFK
                  0 ||||||||...........||||||...|||||.|.|||||||.-.|||------||.|.
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|754297        53 SLDKAENTAQTKKEEKEKEQA 74
                 60 ....|...|.|..|....|.| 81
gi|491764        74 AKTIADKQADTNQEQAKTEDA 95
""",
        )
        hit = hits[22]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|32033565|ref|ZP_00133892.1|")
        self.assertEqual(hit.target.name, "ZP_00133892")
        self.assertEqual(
            hit.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Actinobacillus pleuropneumoniae serovar 1 str. 4074]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=76)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 165.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 68.1662)
        self.assertAlmostEqual(hsp.annotations["evalue"], 7.49305e-11)
        self.assertEqual(hsp.annotations["identity"], 43)
        self.assertEqual(hsp.annotations["positive"], 60)
        self.assertEqual(hsp.annotations["gaps"], 6)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 47, 47, 74],
                          [14, 61, 67, 94]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 80))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...TED'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGGISIWQLLIIVAIIVLLFGTKKLRTLGTDLGESVKGFKKAMADDKSQPQDAS...EKE'}, length=76)",
        )
        self.assertEqual(hsp.target.id, "gi|32033565|ref|ZP_00133892.1|")
        self.assertEqual(hsp.target.name, "ZP_00133892")
        self.assertEqual(
            hsp.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Actinobacillus pleuropneumoniae serovar 1 str. 4074]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGGISIWQLLII  I+VLLFGTKKL ++G+DLG S+KGFKKAM+DD+      SQ  D + + +  K+A + +++AK ++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|320335         0 MGGISIWQLLIIVAIIVLLFGTKKLRTLGTDLGESVKGFKKAMADDK------SQPQDAS
                  0 ||||||||...........||||||...|.|||.|.|||||||.||.------||..|..
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|320335        54 FEKVEAKEAASTEQKAKEKE 74
                 60 ......|.|......||... 80
gi|491764        74 AKTIADKQADTNQEQAKTED 94
""",
        )
        hit = hits[23]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|12722097|gb|AAK03773.1|")
        self.assertEqual(hit.target.name, "AAK03773")
        self.assertEqual(
            hit.target.description,
            "unknown [Pasteurella multocida subsp. multocida str. Pm70] >gi|15603554|ref|NP_246628.1| hypothetical protein PM1689 [Pasteurella multocida subsp. multocida str. Pm70] >gi|24212511|sp|Q9CKD3|TATA_PASMU Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=76)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 162.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 67.0106)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.66928e-10)
        self.assertEqual(hsp.annotations["identity"], 53)
        self.assertEqual(hsp.annotations["positive"], 63)
        self.assertEqual(hsp.annotations["gaps"], 13)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0,  49,  49,  54,  54,  64,  64,  71,  71,  75],
                          [ 14,  63,  68,  73,  75,  85,  89,  96,  98, 102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 88))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...KEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGGISITQLLIIVAIVVLLFGTKKLRTLGSDLGESVKGFKKAMADDNKEKDAEF...KEQ'}, length=76)",
        )
        self.assertEqual(hsp.target.id, "gi|12722097|gb|AAK03773.1|")
        self.assertEqual(hsp.target.name, "AAK03773")
        self.assertEqual(
            hsp.target.description,
            "unknown [Pasteurella multocida subsp. multocida str. Pm70] >gi|15603554|ref|NP_246628.1| hypothetical protein PM1689 [Pasteurella multocida subsp. multocida str. Pm70] >gi|24212511|sp|Q9CKD3|TATA_PASMU Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGGISI QLLII  IVVLLFGTKKL ++GSDLG S+KGFKKAM+DD  +     +DA+F  K+++D    T    AKTE AK  DKEQ",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|127220         0 MGGISITQLLIIVAIVVLLFGTKKLRTLGSDLGESVKGFKKAMADDNKE-----KDAEF-
                  0 ||||||.|...........||||||...|||||.|.|||||||.||...-----.||.|-
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|127220        54 -KSLSDDSETT----AKTEKAK--DKEQ  75
                 60 -|...|....|----||||.||--||||  88
gi|491764        74 AKTIADKQADTNQEQAKTEDAKRHDKEQ 102
""",
        )
        hit = hits[24]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|68546478|ref|ZP_00586025.1|")
        self.assertEqual(hit.target.name, "ZP_00586025")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Shewanella amazonensis SB2B] >gi|68515892|gb|EAN39606.1| Twin-arginine translocation protein TatA/E [Shewanella amazonensis SB2B]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=77)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 161.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 66.6254)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.18014e-10)
        self.assertEqual(hsp.annotations["identity"], 48)
        self.assertEqual(hsp.annotations["positive"], 56)
        self.assertEqual(hsp.annotations["gaps"], 13)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0,  49,  49,  75],
                          [ 14,  63,  76, 102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 88))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...KEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGGISIWQLLIIALIVVLLFGTKKLRSLGGDLGGAIKGFKNAMSDEEKKALEDK...KKE'}, length=77)",
        )
        self.assertEqual(hsp.target.id, "gi|68546478|ref|ZP_00586025.1|")
        self.assertEqual(hsp.target.name, "ZP_00586025")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Shewanella amazonensis SB2B] >gi|68515892|gb|EAN39606.1| Twin-arginine translocation protein TatA/E [Shewanella amazonensis SB2B]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGGISIWQLLIIA+IVVLLFGTKKL S+G DLG +IKGFK AMSD+E K              + DK+A     Q  TE     DK++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|685464         0 MGGISIWQLLIIALIVVLLFGTKKLRSLGGDLGGAIKGFKNAMSDEEKK-----------
                  0 ||||||||...........||||||.|.|.|||..|||||.||||.|.|-----------
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|685464        49 --ALEDKEAAAQTTQQATEKKPEADKKE  75
                 60 --...||.|.....|..||.....||..  88
gi|491764        74 AKTIADKQADTNQEQAKTEDAKRHDKEQ 102
""",
        )
        hit = hits[25]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|33151888|ref|NP_873241.1|")
        self.assertEqual(hit.target.name, "NP_873241")
        self.assertEqual(
            hit.target.description,
            "sec-independent protein secretion pathway component TatA [Haemophilus ducreyi 35000HP] >gi|33148109|gb|AAP95630.1| sec-independent protein secretion pathway component TatA [Haemophilus ducreyi 35000HP] >gi|71153585|sp|Q7VN63|TATA_HAEDU Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=74)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 159.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 65.855)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.71876e-10)
        self.assertEqual(hsp.annotations["identity"], 45)
        self.assertEqual(hsp.annotations["positive"], 58)
        self.assertEqual(hsp.annotations["gaps"], 7)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 47, 47, 54, 54, 74],
                          [14, 61, 67, 74, 75, 95]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 81))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...EDA'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MGGISIWQLLIIVAIIVLLFGTKKLRTLGTDLGESVKGFKKAMAEDKSQDANFD...EQA')",
        )
        self.assertEqual(hsp.target.id, "gi|33151888|ref|NP_873241.1|")
        self.assertEqual(hsp.target.name, "NP_873241")
        self.assertEqual(
            hsp.target.description,
            "sec-independent protein secretion pathway component TatA [Haemophilus ducreyi 35000HP] >gi|33148109|gb|AAP95630.1| sec-independent protein secretion pathway component TatA [Haemophilus ducreyi 35000HP] >gi|71153585|sp|Q7VN63|TATA_HAEDU Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGGISIWQLLII  I+VLLFGTKKL ++G+DLG S+KGFKKAM++D+      SQDA+F  K  A +   T ++  + E A",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|331518         0 MGGISIWQLLIIVAIIVLLFGTKKLRTLGTDLGESVKGFKKAMAEDK------SQDANFD
                  0 ||||||||...........||||||...|.|||.|.|||||||..|.------||||.|.
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|331518        54 -KVEAKESTSTTEKTKEKEQA 74
                 60 -|..|.....|.......|.| 81
gi|491764        74 AKTIADKQADTNQEQAKTEDA 95
""",
        )
        hit = hits[26]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|24375687|ref|NP_719730.1|")
        self.assertEqual(hit.target.name, "NP_719730")
        self.assertEqual(
            hit.target.description,
            "Sec-independent protein translocase protein TatA [Shewanella oneidensis MR-1] >gi|24350613|gb|AAN57174.1| Sec-independent protein translocase protein TatA [Shewanella oneidensis MR-1]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=88)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 159.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 65.855)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.71876e-10)
        self.assertEqual(hsp.annotations["identity"], 48)
        self.assertEqual(hsp.annotations["positive"], 59)
        self.assertEqual(hsp.annotations["gaps"], 1)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0,  76,  76,  87],
                          [ 14,  90,  91, 102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 88))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...KEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGGISIWQLLIIALIVVLLFGTKKLRSLGGDLGGAVKGFKNAMSSEEDKKALED...KEQ'}, length=88)",
        )
        self.assertEqual(hsp.target.id, "gi|24375687|ref|NP_719730.1|")
        self.assertEqual(hsp.target.name, "NP_719730")
        self.assertEqual(
            hsp.target.description,
            "Sec-independent protein translocase protein TatA [Shewanella oneidensis MR-1] >gi|24350613|gb|AAN57174.1| Sec-independent protein translocase protein TatA [Shewanella oneidensis MR-1]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGGISIWQLLIIA+IVVLLFGTKKL S+G DLG ++KGFK AMS +E K+     +A    +T    Q+    +QA TE     +KEQ",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|243756         0 MGGISIWQLLIIALIVVLLFGTKKLRSLGGDLGGAVKGFKNAMSSEEDKKALEDAEAAKP
                  0 ||||||||...........||||||.|.|.|||...||||.|||..|.|.......|...
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|243756        60 VQTAQTVQSAQPTQQA-TEKKPESNKEQ  87
                 60 ..|....|......||-||......|||  88
gi|491764        74 AKTIADKQADTNQEQAKTEDAKRHDKEQ 102
""",
        )
        hit = hits[27]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|71278553|ref|YP_269744.1|")
        self.assertEqual(hit.target.name, "YP_269744")
        self.assertEqual(
            hit.target.description,
            "Sec-independent protein translocase protein TatA [Colwellia psychrerythraea 34H] >gi|71144293|gb|AAZ24766.1| Sec-independent protein translocase protein TatA [Colwellia psychrerythraea 34H]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=79)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 158.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 65.4698)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.85685e-10)
        self.assertEqual(hsp.annotations["identity"], 47)
        self.assertEqual(hsp.annotations["positive"], 62)
        self.assertEqual(hsp.annotations["gaps"], 10)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0,  56,  56,  70,  71,  78],
                          [ 14,  70,  79,  93,  93, 100]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 87))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...HDK'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGGISIWQLVMIAVIVVLLFGTKKLKNLGADLGSSIKGFKKAITDVDDEDRTES...KDK'}, length=79)",
        )
        self.assertEqual(hsp.target.id, "gi|71278553|ref|YP_269744.1|")
        self.assertEqual(hsp.target.name, "YP_269744")
        self.assertEqual(
            hsp.target.description,
            "Sec-independent protein translocase protein TatA [Colwellia psychrerythraea 34H] >gi|71144293|gb|AAZ24766.1| Sec-independent protein translocase protein TatA [Colwellia psychrerythraea 34H]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGGISIWQL++IAVIVVLLFGTKKL ++G+DLG+SIKGFKKA++D + +    S++         D  ADT+  +A+ E D K  DK",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|712785         0 MGGISIWQLVMIAVIVVLLFGTKKLKNLGADLGSSIKGFKKAITDVDDEDRTESKN----
                  0 ||||||||...........||||||...|.|||.||||||||..|........|..----
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|712785        56 -----DVLADTSNTEARVEVDVKEKDK  78
                 60 -----|..|||....|..|-|.|..||  87
gi|491764        74 AKTIADKQADTNQEQAKTE-DAKRHDK 100
""",
        )
        hit = hits[28]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|69159855|gb|EAN71956.1|")
        self.assertEqual(hit.target.name, "EAN71956")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Shewanella denitrificans OS217] >gi|69941450|ref|ZP_00633384.1| Twin-arginine translocation protein TatA/E [Shewanella denitrificans OS-217]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=79)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 157.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 65.0846)
        self.assertAlmostEqual(hsp.annotations["evalue"], 6.34325e-10)
        self.assertEqual(hsp.annotations["identity"], 47)
        self.assertEqual(hsp.annotations["positive"], 62)
        self.assertEqual(hsp.annotations["gaps"], 12)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0,  52,  52,  59,  60,  74,  74,  78],
                          [ 14,  66,  75,  82,  82,  96,  98, 102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 89))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...KEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGGISIWQLLIIALIVILLFGTKKLRSLGGDLGGAVKGFKNAMTSETSEEEKKA...KEQ'}, length=79)",
        )
        self.assertEqual(hsp.target.id, "gi|69159855|gb|EAN71956.1|")
        self.assertEqual(hsp.target.name, "EAN71956")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Shewanella denitrificans OS217] >gi|69941450|ref|ZP_00633384.1| Twin-arginine translocation protein TatA/E [Shewanella denitrificans OS-217]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGGISIWQLLIIA+IV+LLFGTKKL S+G DLG ++KGFK AM+ +  +++K         K + D Q A T+Q+  K  ++K  DKEQ",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|691598         0 MGGISIWQLLIIALIVILLFGTKKLRSLGGDLGGAVKGFKNAMTSETSEEEK--------
                  0 ||||||||...........||||||.|.|.|||...||||.||........|--------
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|691598        52 -KALEDSQTAQTSQQAEKKPESK--DKEQ  78
                 60 -|...|.|-|.|.|...|....|--||||  89
gi|491764        74 AKTIADKQ-ADTNQEQAKTEDAKRHDKEQ 102
""",
        )
        hit = hits[29]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|69949858|ref|ZP_00637822.1|")
        self.assertEqual(hit.target.name, "ZP_00637822")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Shewanella frigidimarina NCIMB 400] >gi|69166482|gb|EAN75453.1| Twin-arginine translocation protein TatA/E [Shewanella frigidimarina NCIMB 400]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=81)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 155.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 64.3142)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.08199e-09)
        self.assertEqual(hsp.annotations["identity"], 41)
        self.assertEqual(hsp.annotations["positive"], 58)
        self.assertEqual(hsp.annotations["gaps"], 2)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 63, 63, 78],
                          [14, 77, 79, 94]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 80))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...TED'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGGISIWQLLIVALIVILLFGTKKLRSLGGDLGGAVKGFKNAMTPEDENKSLDD...SKD'}, length=81)",
        )
        self.assertEqual(hsp.target.id, "gi|69949858|ref|ZP_00637822.1|")
        self.assertEqual(hsp.target.name, "ZP_00637822")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Shewanella frigidimarina NCIMB 400] >gi|69166482|gb|EAN75453.1| Twin-arginine translocation protein TatA/E [Shewanella frigidimarina NCIMB 400]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGGISIWQLLI+A+IV+LLFGTKKL S+G DLG ++KGFK AM+ ++  +    ++ D TA T   +QA   Q + +++D",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|699498         0 MGGISIWQLLIVALIVILLFGTKKLRSLGGDLGGAVKGFKNAMTPEDENKSLDDKEKDQT
                  0 ||||||||...........||||||.|.|.|||...||||.||..............|.|
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|699498        60 AAT--SQQAAEKQPETESKD 78
                 60 |.|--..||...|......| 80
gi|491764        74 AKTIADKQADTNQEQAKTED 94
""",
        )
        hit = hits[30]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|48863844|ref|ZP_00317737.1|")
        self.assertEqual(hit.target.name, "ZP_00317737")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein Mdeg02000203 [Microbulbifer degradans 2-40]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=83)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 154.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 63.929)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.41313e-09)
        self.assertEqual(hsp.annotations["identity"], 41)
        self.assertEqual(hsp.annotations["positive"], 53)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 2, 81],
                          [14, 93]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 79))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...KTE'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({2: 'LGGISIWQLLIVLVIVLLLFGTKRLKGLGGDLGGAIKGFKKAMSDDEAAKQEAE...KTE'}, length=83)",
        )
        self.assertEqual(hsp.target.id, "gi|48863844|ref|ZP_00317737.1|")
        self.assertEqual(hsp.target.name, "ZP_00317737")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein Mdeg02000203 [Microbulbifer degradans 2-40]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "+GGISIWQLLI+ VIV+LLFGTK+L  +G DLG +IKGFKKAMSDDE  + +  +             A T +++ KTE",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|488638         2 LGGISIWQLLIVLVIVLLLFGTKRLKGLGGDLGGAIKGFKKAMSDDEAAKQEAEEAEQKK
                  0 .|||||||...........||||.|...|.|||..||||||||||||.............
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|488638        62 VAAEEAAAAKTAEQKEKTE 81
                 60 ........|.|.....||| 79
gi|491764        74 AKTIADKQADTNQEQAKTE 93
""",
        )
        hit = hits[31]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|77361831|ref|YP_341406.1|")
        self.assertEqual(hit.target.name, "YP_341406")
        self.assertEqual(
            hit.target.description,
            "twin-arginine translocase subunit, sec-independent protein export [Pseudoalteromonas haloplanktis TAC125] >gi|76876742|emb|CAI87964.1| twin-arginine translocase subunit, sec-independent protein export [Pseudoalteromonas haloplanktis TAC125]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=82)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 145.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 60.4622)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.5624e-08)
        self.assertEqual(hsp.annotations["identity"], 39)
        self.assertEqual(hsp.annotations["positive"], 53)
        self.assertEqual(hsp.annotations["gaps"], 6)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 3, 47, 53, 81],
                          [15, 59, 59, 87]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 78))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({15: 'GGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQ...TNQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({3: 'GGISLWQLLIVLAIIVLLFGTKKLRGIGGDLGGAVKGFKKAMSDEKNTDKEKPE...NNK'}, length=82)",
        )
        self.assertEqual(hsp.target.id, "gi|77361831|ref|YP_341406.1|")
        self.assertEqual(hsp.target.name, "YP_341406")
        self.assertEqual(
            hsp.target.description,
            "twin-arginine translocase subunit, sec-independent protein export [Pseudoalteromonas haloplanktis TAC125] >gi|76876742|emb|CAI87964.1| twin-arginine translocase subunit, sec-independent protein export [Pseudoalteromonas haloplanktis TAC125]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "GGIS+WQLLI+  I+VLLFGTKKL  IG DLG ++KGFKKAMSD      ++P+Q + S+++        +K  D N+",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|773618         3 GGISLWQLLIVLAIIVLLFGTKKLRGIGGDLGGAVKGFKKAMSDEKNTDKEKPEQIQKSE
                  0 ||||.||...........||||||..||.|||...|||||||||------..|.|...|.
gi|491764        15 GGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSD------DEPKQDKTSQ

gi|773618        63 ESAPLDSAHTEKNKDNNK 81
                 60 ...........|..|.|. 78
gi|491764        69 DADFTAKTIADKQADTNQ 87
""",
        )
        hit = hits[32]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|67676224|ref|ZP_00472975.1|")
        self.assertEqual(hit.target.name, "ZP_00472975")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Chromohalobacter salexigens DSM 3043] >gi|67519761|gb|EAM23710.1| Twin-arginine translocation protein TatA/E [Chromohalobacter salexigens DSM 3043]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=90)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 144.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 60.077)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.04055e-08)
        self.assertEqual(hsp.annotations["identity"], 40)
        self.assertEqual(hsp.annotations["positive"], 55)
        self.assertEqual(hsp.annotations["gaps"], 6)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 1, 51, 51, 62, 62, 73],
                          [14, 64, 68, 79, 81, 92]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 78))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...AKT'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({1: 'LGGISIWQLLIVLGIVILIFGTKKLRNLGSDLGGAVKGFKSAVGEEEQKKQDGE...LDT'}, length=90)",
        )
        self.assertEqual(hsp.target.id, "gi|67676224|ref|ZP_00472975.1|")
        self.assertEqual(hsp.target.name, "ZP_00472975")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Chromohalobacter salexigens DSM 3043] >gi|67519761|gb|EAM23710.1| Twin-arginine translocation protein TatA/E [Chromohalobacter salexigens DSM 3043]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "+GGISIWQLLI+  IV+L+FGTKKL ++GSDLG ++KGFK A+ ++E K+    QD + TA+       DTN+E   T",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|676762         1 LGGISIWQLLIVLGIVILIFGTKKLRNLGSDLGGAVKGFKSAVGEEEQKK----QDGEET
                  0 .|||||||...........||||||...|||||...||||.|....|.|.----||...|
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|676762        57 AQHRV--AHDTNEEPLDT 73
                 60 |....--..|||.|...| 78
gi|491764        74 AKTIADKQADTNQEQAKT 92
""",
        )
        hit = hits[33]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|74317722|ref|YP_315462.1|")
        self.assertEqual(hit.target.name, "YP_315462")
        self.assertEqual(
            hit.target.description,
            "twin-arginine translocation protein TatA/E [Thiobacillus denitrificans ATCC 25259] >gi|74057217|gb|AAZ97657.1| twin-arginine translocation protein TatA/E [Thiobacillus denitrificans ATCC 25259] >gi|52007795|ref|ZP_00335172.1| COG1826: Sec-independent protein secretion pathway components [Thiobacillus denitrificans ATCC 25259]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=70)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 142.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 59.3066)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.48066e-08)
        self.assertEqual(hsp.annotations["identity"], 34)
        self.assertEqual(hsp.annotations["positive"], 46)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 68],
                          [14, 82]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 68))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...DKQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGSFSIWHWLIVLVVVLLIFGTKKLRNIGSDLGGAVKGFKEGMKDDAPKISESD...DKQ'}, length=70)",
        )
        self.assertEqual(hsp.target.id, "gi|74317722|ref|YP_315462.1|")
        self.assertEqual(hsp.target.name, "YP_315462")
        self.assertEqual(
            hsp.target.description,
            "twin-arginine translocation protein TatA/E [Thiobacillus denitrificans ATCC 25259] >gi|74057217|gb|AAZ97657.1| twin-arginine translocation protein TatA/E [Thiobacillus denitrificans ATCC 25259] >gi|52007795|ref|ZP_00335172.1| COG1826: Sec-independent protein secretion pathway components [Thiobacillus denitrificans ATCC 25259]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG  SIW  LI+ V+V+L+FGTKKL +IGSDLG ++KGFK+ M DD PK  ++ +        + DKQ",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|743177         0 MGSFSIWHWLIVLVVVLLIFGTKKLRNIGSDLGGAVKGFKEGMKDDAPKISESDKGGHTI
                  0 ||..|||............||||||..||||||...||||..|.||.||...........
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|743177        60 DAEVKDKQ 68
                 60 .....||| 68
gi|491764        74 AKTIADKQ 82
""",
        )
        hit = hits[34]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|77166504|ref|YP_345029.1|")
        self.assertEqual(hit.target.name, "YP_345029")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Nitrosococcus oceani ATCC 19707] >gi|76884818|gb|ABA59499.1| Twin-arginine translocation protein TatA/E [Nitrosococcus oceani ATCC 19707]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=90)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 140.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 58.5362)
        self.assertAlmostEqual(hsp.annotations["evalue"], 5.9371e-08)
        self.assertEqual(hsp.annotations["identity"], 41)
        self.assertEqual(hsp.annotations["positive"], 58)
        self.assertEqual(hsp.annotations["gaps"], 7)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 2, 49, 49, 54, 54, 76],
                          [14, 61, 67, 72, 73, 95]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 81))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...EDA'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({2: 'VGGISIWQLLIVLVIILLLFGTKKLRSIGTDLGSAIKGFRNSLRDEERRDAEEA...QQA'}, length=90)",
        )
        self.assertEqual(hsp.target.id, "gi|77166504|ref|YP_345029.1|")
        self.assertEqual(hsp.target.name, "YP_345029")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Nitrosococcus oceani ATCC 19707] >gi|76884818|gb|ABA59499.1| Twin-arginine translocation protein TatA/E [Nitrosococcus oceani ATCC 19707]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "+GGISIWQLLI+ VI++LLFGTKKL SIG+DLG++IKGF+ ++ D+E       +DA+  A TI  KQA   +  ++ + A",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|771665         2 VGGISIWQLLIVLVIILLLFGTKKLRSIGTDLGSAIKGFRNSLRDEE------RRDAE-E
                  0 .|||||||...........||||||.|||.|||..||||.....|.|------..||.-.
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|771665        55 AATIEHKQAHKAENPSQRQQA 76
                 60 |.||..|||...........| 81
gi|491764        74 AKTIADKQADTNQEQAKTEDA 95
""",
        )
        hit = hits[35]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|16128610|ref|NP_415160.1|")
        self.assertEqual(hit.target.name, "NP_415160")
        self.assertEqual(
            hit.target.description,
            "component of Sec-independent translocase [Escherichia coli K12] >gi|26246607|ref|NP_752647.1| Sec-independent protein translocase protein tatE [Escherichia coli CFT073] >gi|30062184|ref|NP_836355.1| hypothetical protein S0676 [Shigella flexneri 2a str. 2457T] >gi|24112073|ref|NP_706583.1| hypothetical protein SF0654 [Shigella flexneri 2a str. 301] >gi|15800341|ref|NP_286353.1| hypothetical protein Z0772 [Escherichia coli O157:H7 EDL933] >gi|24050900|gb|AAN42290.1| orf, conserved hypothetical protein [Shigella flexneri 2a str. 301] >gi|74311162|ref|YP_309581.1| hypothetical protein SSO_0580 [Shigella sonnei Ss046] >gi|1786845|gb|AAC73728.1| component of Sec-independent translocase [Escherichia coli K12] >gi|75514720|ref|ZP_00736975.1| COG1826: Sec-independent protein secretion pathway components [Escherichia coli 53638] >gi|75258637|ref|ZP_00730044.1| COG1826: Sec-independent protein secretion pathway components [Escherichia coli E22] >gi|75240782|ref|ZP_00724689.1| COG1826: Sec-independent protein secretion pathway components [Escherichia coli F11] >gi|75236608|ref|ZP_00720696.1| COG1826: Sec-independent protein secretion pathway components [Escherichia coli E110019] >gi|75229489|ref|ZP_00716036.1| COG1826: Sec-independent protein secretion pathway components [Escherichia coli B7A] >gi|75210746|ref|ZP_00710878.1| COG1826: Sec-independent protein secretion pathway components [Escherichia coli B171] >gi|75194810|ref|ZP_00704880.1| COG1826: Sec-independent protein secretion pathway components [Escherichia coli HS] >gi|75189665|ref|ZP_00702932.1| COG1826: Sec-independent protein secretion pathway components [Escherichia coli E24377A] >gi|75177456|ref|ZP_00697537.1| COG1826: Sec-independent protein secretion pathway components [Shigella boydii BS512] >gi|73854639|gb|AAZ87346.1| conserved hypothetical protein [Shigella sonnei Ss046] >gi|30040429|gb|AAP16161.1| hypothetical protein S0676 [Shigella flexneri 2a str. 2457T] >gi|13360123|dbj|BAB34088.1| Sec-independent protein translocase [Escherichia coli O157:H7] >gi|1778544|gb|AAB40827.1| HI0187 homolog [Escherichia coli] >gi|26107006|gb|AAN79190.1| Sec-independent protein translocase protein tatE [Escherichia coli CFT073] >gi|4062250|dbj|BAA35270.1| Hypothetical protein (lip 3' region) [Escherichia coli K12] >gi|67473128|sp|P0A843|TATE_ECOLI Sec-independent protein translocase protein tatE >gi|67473131|sp|P0A846|TATE_SHIFL Sec-independent protein translocase protein tatE >gi|67473130|sp|P0A845|TATE_ECO57 Sec-independent protein translocase protein tatE >gi|67473129|sp|P0A844|TATE_ECOL6 Sec-independent protein translocase protein tatE >gi|12513528|gb|AAG54961.1| orf, hypothetical protein [Escherichia coli O157:H7 EDL933] >gi|15829919|ref|NP_308692.1| Sec-independent protein translocase [Escherichia coli O157:H7]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=67)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 137.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 57.3806)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.32265e-07)
        self.assertEqual(hsp.annotations["identity"], 38)
        self.assertEqual(hsp.annotations["positive"], 51)
        self.assertEqual(hsp.annotations["gaps"], 1)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 50, 50, 67],
                          [14, 64, 65, 82]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 68))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...DKQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MGEISITKLLVVAALVVLLFGTKKLRTLGGDLGAAIKGFKKAMNDDDAAAKKGA...HKE')",
        )
        self.assertEqual(hsp.target.id, "gi|16128610|ref|NP_415160.1|")
        self.assertEqual(hsp.target.name, "NP_415160")
        self.assertEqual(
            hsp.target.description,
            "component of Sec-independent translocase [Escherichia coli K12] >gi|26246607|ref|NP_752647.1| Sec-independent protein translocase protein tatE [Escherichia coli CFT073] >gi|30062184|ref|NP_836355.1| hypothetical protein S0676 [Shigella flexneri 2a str. 2457T] >gi|24112073|ref|NP_706583.1| hypothetical protein SF0654 [Shigella flexneri 2a str. 301] >gi|15800341|ref|NP_286353.1| hypothetical protein Z0772 [Escherichia coli O157:H7 EDL933] >gi|24050900|gb|AAN42290.1| orf, conserved hypothetical protein [Shigella flexneri 2a str. 301] >gi|74311162|ref|YP_309581.1| hypothetical protein SSO_0580 [Shigella sonnei Ss046] >gi|1786845|gb|AAC73728.1| component of Sec-independent translocase [Escherichia coli K12] >gi|75514720|ref|ZP_00736975.1| COG1826: Sec-independent protein secretion pathway components [Escherichia coli 53638] >gi|75258637|ref|ZP_00730044.1| COG1826: Sec-independent protein secretion pathway components [Escherichia coli E22] >gi|75240782|ref|ZP_00724689.1| COG1826: Sec-independent protein secretion pathway components [Escherichia coli F11] >gi|75236608|ref|ZP_00720696.1| COG1826: Sec-independent protein secretion pathway components [Escherichia coli E110019] >gi|75229489|ref|ZP_00716036.1| COG1826: Sec-independent protein secretion pathway components [Escherichia coli B7A] >gi|75210746|ref|ZP_00710878.1| COG1826: Sec-independent protein secretion pathway components [Escherichia coli B171] >gi|75194810|ref|ZP_00704880.1| COG1826: Sec-independent protein secretion pathway components [Escherichia coli HS] >gi|75189665|ref|ZP_00702932.1| COG1826: Sec-independent protein secretion pathway components [Escherichia coli E24377A] >gi|75177456|ref|ZP_00697537.1| COG1826: Sec-independent protein secretion pathway components [Shigella boydii BS512] >gi|73854639|gb|AAZ87346.1| conserved hypothetical protein [Shigella sonnei Ss046] >gi|30040429|gb|AAP16161.1| hypothetical protein S0676 [Shigella flexneri 2a str. 2457T] >gi|13360123|dbj|BAB34088.1| Sec-independent protein translocase [Escherichia coli O157:H7] >gi|1778544|gb|AAB40827.1| HI0187 homolog [Escherichia coli] >gi|26107006|gb|AAN79190.1| Sec-independent protein translocase protein tatE [Escherichia coli CFT073] >gi|4062250|dbj|BAA35270.1| Hypothetical protein (lip 3' region) [Escherichia coli K12] >gi|67473128|sp|P0A843|TATE_ECOLI Sec-independent protein translocase protein tatE >gi|67473131|sp|P0A846|TATE_SHIFL Sec-independent protein translocase protein tatE >gi|67473130|sp|P0A845|TATE_ECO57 Sec-independent protein translocase protein tatE >gi|67473129|sp|P0A844|TATE_ECOL6 Sec-independent protein translocase protein tatE >gi|12513528|gb|AAG54961.1| orf, hypothetical protein [Escherichia coli O157:H7 EDL933] >gi|15829919|ref|NP_308692.1| Sec-independent protein translocase [Escherichia coli O157:H7]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG ISI +LL++A +VVLLFGTKKL ++G DLGA+IKGFKKAM+DD+    K   D D  A+ ++ K+",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|161286         0 MGEISITKLLVVAALVVLLFGTKKLRTLGGDLGAAIKGFKKAMNDDDAAA-KKGADVDLQ
                  0 ||.|||.............||||||...|.||||.||||||||.||....-|...|.|..
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|161286        59 AEKLSHKE 67
                 60 |.....|. 68
gi|491764        74 AKTIADKQ 82
""",
        )
        hit = hits[36]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|12831974|emb|CAC29147.1|")
        self.assertEqual(hit.target.name, "CAC29147")
        self.assertEqual(hit.target.description, "TatA protein [Pseudomonas stutzeri]")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=76)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 136.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 56.9954)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.72743e-07)
        self.assertEqual(hsp.annotations["identity"], 36)
        self.assertEqual(hsp.annotations["positive"], 50)
        self.assertEqual(hsp.annotations["gaps"], 3)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 3, 51, 54, 71],
                          [15, 63, 63, 80]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 68))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({15: 'GGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQ...IAD'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({3: 'GGISIWQLLIILLIVIMLFGTKRLKGVGSDLGDAIKGFRKSMGTDEEKPSVEEK...VEE'}, length=76)",
        )
        self.assertEqual(hsp.target.id, "gi|12831974|emb|CAC29147.1|")
        self.assertEqual(hsp.target.name, "CAC29147")
        self.assertEqual(hsp.target.description, "TatA protein [Pseudomonas stutzeri]")
        self.assertEqual(
            hsp.annotations["midline"],
            "GGISIWQLLII +IV++LFGTK+L  +GSDLG +IKGF+K+M  DE K   ++K +   D  A+ + +",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|128319         3 GGISIWQLLIILLIVIMLFGTKRLKGVGSDLGDAIKGFRKSMGTDEEKPSVEEKQNHTID
                  0 |||||||...........||||.|...|||||..||||.|.|..||.|---..|.....|
gi|491764        15 GGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPK---QDKTSQDAD

gi|128319        63 AQARKVEE 71
                 60 ..|..... 68
gi|491764        72 FTAKTIAD 80
""",
        )
        hit = hits[37]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|32029972|ref|ZP_00132908.1|")
        self.assertEqual(hit.target.name, "ZP_00132908")
        self.assertEqual(
            hit.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Haemophilus somnus 2336] >gi|23467861|ref|ZP_00123438.1| COG1826: Sec-independent protein secretion pathway components [Haemophilus somnus 129PT]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=73)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 135.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 56.6102)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.2561e-07)
        self.assertEqual(hsp.annotations["identity"], 38)
        self.assertEqual(hsp.annotations["positive"], 57)
        self.assertEqual(hsp.annotations["gaps"], 8)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 1, 48, 48, 52, 52, 71],
                          [16, 63, 69, 73, 75, 94]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 78))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({16: 'GISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQD...TED'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({1: 'GLSWQQLLILLLVVVVIFGTKKLRNIGSDLGGAVKDFKKAMNDDQPKDAEFKKI...QKE'}, length=73)",
        )
        self.assertEqual(hsp.target.id, "gi|32029972|ref|ZP_00132908.1|")
        self.assertEqual(hsp.target.name, "ZP_00132908")
        self.assertEqual(
            hsp.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Haemophilus somnus 2336] >gi|23467861|ref|ZP_00123438.1| COG1826: Sec-independent protein secretion pathway components [Haemophilus somnus 129PT]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "G+S  QLLI+ ++VV++FGTKKL +IGSDLG ++K FKKAM+DD+PK      DA+F  K I+++   T+ E +K ++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|320299         1 GLSWQQLLILLLVVVVIFGTKKLRNIGSDLGGAVKDFKKAMNDDQPK------DAEF--K
                  0 |.|..|...........||||||..||||||...|.|||||.||.||------||.|--|
gi|491764        16 GISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAK

gi|320299        53 KISEEVEQTSVENSKQKE 71
                 60 .|......|..|..|... 78
gi|491764        76 TIADKQADTNQEQAKTED 94
""",
        )
        hit = hits[38]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|455172|gb|AAA24073.1|")
        self.assertEqual(hit.target.name, "AAA24073")
        self.assertEqual(hit.target.description, "ORF; putative")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=67)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 134.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 56.225)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.94655e-07)
        self.assertEqual(hsp.annotations["identity"], 38)
        self.assertEqual(hsp.annotations["positive"], 50)
        self.assertEqual(hsp.annotations["gaps"], 1)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 50, 50, 67],
                          [14, 64, 65, 82]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 68))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...DKQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MGEISITKLLVVAALVVLLFGTKKLRTLGGDLGAAIKGFKKAMIDDDAAAKKGA...HKE')",
        )
        self.assertEqual(hsp.target.id, "gi|455172|gb|AAA24073.1|")
        self.assertEqual(hsp.target.name, "AAA24073")
        self.assertEqual(hsp.target.description, "ORF; putative")
        self.assertEqual(
            hsp.annotations["midline"],
            "MG ISI +LL++A +VVLLFGTKKL ++G DLGA+IKGFKKAM DD+    K   D D  A+ ++ K+",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|455172         0 MGEISITKLLVVAALVVLLFGTKKLRTLGGDLGAAIKGFKKAMIDDDAAA-KKGADVDLQ
                  0 ||.|||.............||||||...|.||||.||||||||.||....-|...|.|..
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|455172        59 AEKLSHKE 67
                 60 |.....|. 68
gi|491764        74 AKTIADKQ 82
""",
        )
        hit = hits[39]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|1224007|gb|AAA92108.1|")
        self.assertEqual(hit.target.name, "AAA92108")
        self.assertEqual(hit.target.description, "ORF4")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=192)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 133.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 55.8398)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.84832e-07)
        self.assertEqual(hsp.annotations["identity"], 36)
        self.assertEqual(hsp.annotations["positive"], 48)
        self.assertEqual(hsp.annotations["gaps"], 4)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 3, 49, 53, 62],
                          [15, 61, 61, 70]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 59))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({15: 'GGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQD'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({3: 'GGISIWQLLIILLIVVMLFGTKRLKSLGSDLGDAIKGFRKSMDNEENKAPPVEEQKGQD'}, length=192)",
        )
        self.assertEqual(hsp.target.id, "gi|1224007|gb|AAA92108.1|")
        self.assertEqual(hsp.target.name, "AAA92108")
        self.assertEqual(hsp.target.description, "ORF4")
        self.assertEqual(
            hsp.annotations["midline"],
            "GGISIWQLLII +IVV+LFGTK+L S+GSDLG +IKGF+K+M ++E    P +++  QD",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|122400         3 GGISIWQLLIILLIVVMLFGTKRLKSLGSDLGDAIKGFRKSMDNEENKAPPVEEQKGQD
                  0 |||||||...........||||.|.|.|||||..||||.|.|...|----|......||
gi|491764        15 GGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDE----PKQDKTSQD

gi|122400        62
                 59
gi|491764        70
""",
        )
        hit = hits[40]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|68056990|gb|AAX87243.1|")
        self.assertEqual(hit.target.name, "AAX87243")
        self.assertEqual(
            hit.target.description,
            "Sec-independent protein translocase protein TatA/E [Haemophilus influenzae 86-028NP] >gi|68248791|ref|YP_247903.1| Sec-independent protein translocase protein TatA/E [Haemophilus influenzae 86-028NP]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=95)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 131.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 55.0694)
        self.assertAlmostEqual(hsp.annotations["evalue"], 6.56423e-07)
        self.assertEqual(hsp.annotations["identity"], 39)
        self.assertEqual(hsp.annotations["positive"], 60)
        self.assertEqual(hsp.annotations["gaps"], 7)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 7, 63, 63, 74, 74, 93],
                          [ 1, 57, 62, 73, 75, 94]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 93))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({1: 'RLCLIIIYHRGTCMGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKK...TED'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({7: 'RAKFFLFYRTEFIMFGLSPAQLIILLVVILLIFGTKKLRNAGSDLGAAVKGFKK...EKE'}, length=95)",
        )
        self.assertEqual(hsp.target.id, "gi|68056990|gb|AAX87243.1|")
        self.assertEqual(hsp.target.name, "AAX87243")
        self.assertEqual(
            hsp.target.description,
            "Sec-independent protein translocase protein TatA/E [Haemophilus influenzae 86-028NP] >gi|68248791|ref|YP_247903.1| Sec-independent protein translocase protein TatA/E [Haemophilus influenzae 86-028NP]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "R    + Y     M G+S  QL+I+ V+++L+FGTKKL + GSDLGA++KGFKKAM     K+D+  +DA+F  K+I ++ A   +E  K ++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|680569         7 RAKFFLFYRTEFIMFGLSPAQLIILLVVILLIFGTKKLRNAGSDLGAAVKGFKKAM----
                  0 |......|.....|.|.|..|...........||||||...||||||..|||||||----
gi|491764         1 RLCLIIIYHRGTCMGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDE

gi|680569        63 -KEDEKVKDAEF--KSIDNETASAKKENIKEKE 93
                 60 -|.|....||.|--|.|....|....|..|... 93
gi|491764        61 PKQDKTSQDADFTAKTIADKQADTNQEQAKTED 94
""",
        )
        hit = hits[41]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|56461470|ref|YP_156751.1|")
        self.assertEqual(hit.target.name, "YP_156751")
        self.assertEqual(
            hit.target.description,
            "Sec-independent protein secretion pathway component, TatA family [Idiomarina loihiensis L2TR] >gi|56180480|gb|AAV83202.1| Sec-independent protein secretion pathway component, TatA family [Idiomarina loihiensis L2TR]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=73)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 129.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 54.299)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.11969e-06)
        self.assertEqual(hsp.annotations["identity"], 33)
        self.assertEqual(hsp.annotations["positive"], 43)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 1, 51],
                          [16, 66]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 50))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({16: 'GISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDK'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({1: 'GVSPWQLLILLVIVVLLFGTKRLRSLGSDLGNAVKGFKKSMGDEDDSKDK'}, length=73)",
        )
        self.assertEqual(hsp.target.id, "gi|56461470|ref|YP_156751.1|")
        self.assertEqual(hsp.target.name, "YP_156751")
        self.assertEqual(
            hsp.target.description,
            "Sec-independent protein secretion pathway component, TatA family [Idiomarina loihiensis L2TR] >gi|56180480|gb|AAV83202.1| Sec-independent protein secretion pathway component, TatA family [Idiomarina loihiensis L2TR]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "G+S WQLLI+ VIVVLLFGTK+L S+GSDLG ++KGFKK+M D++  +DK",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|564614         1 GVSPWQLLILLVIVVLLFGTKRLRSLGSDLGNAVKGFKKSMGDEDDSKDK 51
                  0 |.|.||...........||||.|.|.|||||...|||||.|.|.....|| 50
gi|491764        16 GISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDK 66
""",
        )
        hit = hits[42]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|76793313|ref|ZP_00775802.1|")
        self.assertEqual(hit.target.name, "ZP_00775802")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Pseudoalteromonas atlantica T6c] >gi|76591415|gb|EAO67616.1| Twin-arginine translocation protein TatA/E [Pseudoalteromonas atlantica T6c]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=84)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 128.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 53.9138)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.46236e-06)
        self.assertEqual(hsp.annotations["identity"], 43)
        self.assertEqual(hsp.annotations["positive"], 57)
        self.assertEqual(hsp.annotations["gaps"], 11)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0,  47,  47,  65,  69,  84],
                          [ 14,  61,  68,  86,  86, 101]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 91))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...DKE'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MFGISPMQLIIILVIVVLLFGTKKLRNMGGDLGSAVKGFKKAVSDDDKDADFKA...KSE')",
        )
        self.assertEqual(hsp.target.id, "gi|76793313|ref|ZP_00775802.1|")
        self.assertEqual(hsp.target.name, "ZP_00775802")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Pseudoalteromonas atlantica T6c] >gi|76591415|gb|EAO67616.1| Twin-arginine translocation protein TatA/E [Pseudoalteromonas atlantica T6c]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "M GIS  QL+II VIVVLLFGTKKL ++G DLG+++KGFKKA+SDD+       +DADF A    +  +  N    Q+  K+E   +   E",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|767933         0 MFGISPMQLIIILVIVVLLFGTKKLRNMGGDLGSAVKGFKKAVSDDD-------KDADFK
                  0 |.|||..|...........||||||...|.|||...||||||.|||.-------.||||.
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|767933        53 ADEKVEDNSKANDQTVQQNVKSESDTKAKSE  84
                 60 |..........|----|...|.|.......|  91
gi|491764        74 AKTIADKQADTN----QEQAKTEDAKRHDKE 101
""",
        )
        hit = hits[43]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|42630489|ref|ZP_00156028.1|")
        self.assertEqual(hit.target.name, "ZP_00156028")
        self.assertEqual(
            hit.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Haemophilus influenzae R2866] >gi|42629145|ref|ZP_00154694.1| COG1826: Sec-independent protein secretion pathway components [Haemophilus influenzae R2846]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=75)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 127.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 53.5286)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.9099e-06)
        self.assertEqual(hsp.annotations["identity"], 37)
        self.assertEqual(hsp.annotations["positive"], 57)
        self.assertEqual(hsp.annotations["gaps"], 7)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 43, 43, 54, 54, 73],
                          [14, 57, 62, 73, 75, 94]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 80))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...TED'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MFGLSPAQLIILLVVILLIFGTKKLRNAGSDLGAAVKGFKKAMKEDEKVKDAEF...EKE'}, length=75)",
        )
        self.assertEqual(hsp.target.id, "gi|42630489|ref|ZP_00156028.1|")
        self.assertEqual(hsp.target.name, "ZP_00156028")
        self.assertEqual(
            hsp.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Haemophilus influenzae R2866] >gi|42629145|ref|ZP_00154694.1| COG1826: Sec-independent protein secretion pathway components [Haemophilus influenzae R2846]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "M G+S  QL+I+ V+++L+FGTKKL + GSDLGA++KGFKKAM     K+D+  +DA+F  K+I ++ A   +E  K ++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|426304         0 MFGLSPAQLIILLVVILLIFGTKKLRNAGSDLGAAVKGFKKAM-----KEDEKVKDAEF-
                  0 |.|.|..|...........||||||...||||||..|||||||-----|.|....||.|-
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|426304        54 -KSIDNETASAKKENIKEKE 73
                 60 -|.|....|....|..|... 80
gi|491764        74 AKTIADKQADTNQEQAKTED 94
""",
        )
        hit = hits[44]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|1074302|pir||B64145")
        self.assertEqual(hit.target.name, "B64145")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein HI0187 - Haemophilus influenzae (strain Rd KW20)",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=109)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 126.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 53.1434)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.49441e-06)
        self.assertEqual(hsp.annotations["identity"], 39)
        self.assertEqual(hsp.annotations["positive"], 59)
        self.assertEqual(hsp.annotations["gaps"], 7)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 7, 63, 63, 74, 74, 92],
                          [ 1, 57, 62, 73, 75, 93]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 92))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({1: 'RLCLIIIYHRGTCMGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKK...KTE'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({7: 'RAKFFLFYRTEFIMFGLSPAQLIILLVVILLIFGTKKLRNAGSDLGAAVKGFKK...KRE'}, length=109)",
        )
        self.assertEqual(hsp.target.id, "gi|1074302|pir||B64145")
        self.assertEqual(hsp.target.name, "B64145")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein HI0187 - Haemophilus influenzae (strain Rd KW20)",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "R    + Y     M G+S  QL+I+ V+++L+FGTKKL + GSDLGA++KGFKKAM     K+D+  +DA+F  K+I ++ A   + + K E",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|107430         7 RAKFFLFYRTEFIMFGLSPAQLIILLVVILLIFGTKKLRNAGSDLGAAVKGFKKAM----
                  0 |......|.....|.|.|..|...........||||||...||||||..|||||||----
gi|491764         1 RLCLIIIYHRGTCMGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDE

gi|107430        63 -KEDEKVKDAEF--KSIDNETASAKKGKYKRE 92
                 60 -|.|....||.|--|.|....|.......|.| 92
gi|491764        61 PKQDKTSQDADFTAKTIADKQADTNQEQAKTE 93
""",
        )
        hit = hits[45]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|67641583|ref|ZP_00440359.1|")
        self.assertEqual(hit.target.name, "ZP_00440359")
        self.assertEqual(
            hit.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Burkholderia mallei GB8 horse 4] >gi|67635791|ref|ZP_00434743.1| COG1826: Sec-independent protein secretion pathway components [Burkholderia mallei 10399] >gi|67760276|ref|ZP_00498999.1| COG1826: Sec-independent protein secretion pathway components [Burkholderia pseudomallei S13] >gi|67756151|ref|ZP_00495038.1| COG1826: Sec-independent protein secretion pathway components [Burkholderia pseudomallei Pasteur] >gi|67737955|ref|ZP_00488671.1| COG1826: Sec-independent protein secretion pathway components [Burkholderia pseudomallei 668] >gi|67716050|ref|ZP_00485411.1| COG1826: Sec-independent protein secretion pathway components [Burkholderia pseudomallei 1710b] >gi|67684327|ref|ZP_00478308.1| COG1826: Sec-independent protein secretion pathway components [Burkholderia pseudomallei 1710a] >gi|67671615|ref|ZP_00468402.1| COG1826: Sec-independent protein secretion pathway components [Burkholderia pseudomallei 1655] >gi|67651364|ref|ZP_00448787.1| COG1826: Sec-independent protein secretion pathway components [Burkholderia mallei SAVP1] >gi|67646939|ref|ZP_00445189.1| COG1826: Sec-independent protein secretion pathway components [Burkholderia mallei NCTC 10247] >gi|67629330|ref|ZP_00429188.1| COG1826: Sec-independent protein secretion pathway components [Burkholderia mallei 10229] >gi|76809713|ref|YP_335046.1| TatA [Burkholderia pseudomallei 1710b] >gi|76579166|gb|ABA48641.1| TatA [Burkholderia pseudomallei 1710b] >gi|52427678|gb|AAU48271.1| twin-arginine translocation protein, TatA/E family [Burkholderia mallei ATCC 23344] >gi|52211149|emb|CAH37138.1| Sec-independent protein translocase protein TatA [Burkholderia pseudomallei K96243] >gi|53720735|ref|YP_109721.1| Sec-independent protein translocase protein TatA [Burkholderia pseudomallei K96243] >gi|69991749|ref|ZP_00643311.1| COG1826: Sec-independent protein secretion pathway components [Burkholderia mallei FMH]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=77)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 125.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 52.7582)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.25779e-06)
        self.assertEqual(hsp.annotations["identity"], 29)
        self.assertEqual(hsp.annotations["positive"], 39)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 51],
                          [14, 65]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 51))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQD'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGGLSIWHWLIVLLIVALVFGTKKLRNIGNDLGSAVKGFKDGMKESEAPAD'}, length=77)",
        )
        self.assertEqual(hsp.target.id, "gi|67641583|ref|ZP_00440359.1|")
        self.assertEqual(hsp.target.name, "ZP_00440359")
        self.assertEqual(
            hsp.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Burkholderia mallei GB8 horse 4] >gi|67635791|ref|ZP_00434743.1| COG1826: Sec-independent protein secretion pathway components [Burkholderia mallei 10399] >gi|67760276|ref|ZP_00498999.1| COG1826: Sec-independent protein secretion pathway components [Burkholderia pseudomallei S13] >gi|67756151|ref|ZP_00495038.1| COG1826: Sec-independent protein secretion pathway components [Burkholderia pseudomallei Pasteur] >gi|67737955|ref|ZP_00488671.1| COG1826: Sec-independent protein secretion pathway components [Burkholderia pseudomallei 668] >gi|67716050|ref|ZP_00485411.1| COG1826: Sec-independent protein secretion pathway components [Burkholderia pseudomallei 1710b] >gi|67684327|ref|ZP_00478308.1| COG1826: Sec-independent protein secretion pathway components [Burkholderia pseudomallei 1710a] >gi|67671615|ref|ZP_00468402.1| COG1826: Sec-independent protein secretion pathway components [Burkholderia pseudomallei 1655] >gi|67651364|ref|ZP_00448787.1| COG1826: Sec-independent protein secretion pathway components [Burkholderia mallei SAVP1] >gi|67646939|ref|ZP_00445189.1| COG1826: Sec-independent protein secretion pathway components [Burkholderia mallei NCTC 10247] >gi|67629330|ref|ZP_00429188.1| COG1826: Sec-independent protein secretion pathway components [Burkholderia mallei 10229] >gi|76809713|ref|YP_335046.1| TatA [Burkholderia pseudomallei 1710b] >gi|76579166|gb|ABA48641.1| TatA [Burkholderia pseudomallei 1710b] >gi|52427678|gb|AAU48271.1| twin-arginine translocation protein, TatA/E family [Burkholderia mallei ATCC 23344] >gi|52211149|emb|CAH37138.1| Sec-independent protein translocase protein TatA [Burkholderia pseudomallei K96243] >gi|53720735|ref|YP_109721.1| Sec-independent protein translocase protein TatA [Burkholderia pseudomallei K96243] >gi|69991749|ref|ZP_00643311.1| COG1826: Sec-independent protein secretion pathway components [Burkholderia mallei FMH]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGG+SIW  LI+ +IV L+FGTKKL +IG+DLG+++KGFK  M + E   D",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|676415         0 MGGLSIWHWLIVLLIVALVFGTKKLRNIGNDLGSAVKGFKDGMKESEAPAD 51
                  0 |||.|||............||||||..||.|||...||||..|...|...| 51
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQD 65
""",
        )
        hit = hits[46]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|67545726|ref|ZP_00423646.1|")
        self.assertEqual(hit.target.name, "ZP_00423646")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Burkholderia vietnamiensis G4] >gi|67533055|gb|EAM29828.1| Twin-arginine translocation protein TatA/E [Burkholderia vietnamiensis G4] >gi|74019444|ref|ZP_00690060.1| Twin-arginine translocation protein TatA/E [Burkholderia ambifaria AMMD] >gi|72607861|gb|EAO43817.1| Twin-arginine translocation protein TatA/E [Burkholderia ambifaria AMMD]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=76)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 125.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 52.7582)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.25779e-06)
        self.assertEqual(hsp.annotations["identity"], 35)
        self.assertEqual(hsp.annotations["positive"], 51)
        self.assertEqual(hsp.annotations["gaps"], 8)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 43, 43, 58, 58, 75],
                          [14, 57, 62, 77, 80, 97]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 83))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...AKR'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGGLSIWHWLIVLLIVALVFGTKKLRNIGNDLGSAVKGFKDGMKEGETPADAQQ...SNK'}, length=76)",
        )
        self.assertEqual(hsp.target.id, "gi|67545726|ref|ZP_00423646.1|")
        self.assertEqual(hsp.target.name, "ZP_00423646")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Burkholderia vietnamiensis G4] >gi|67533055|gb|EAM29828.1| Twin-arginine translocation protein TatA/E [Burkholderia vietnamiensis G4] >gi|74019444|ref|ZP_00690060.1| Twin-arginine translocation protein TatA/E [Burkholderia ambifaria AMMD] >gi|72607861|gb|EAO43817.1| Twin-arginine translocation protein TatA/E [Burkholderia ambifaria AMMD]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGG+SIW  LI+ +IV L+FGTKKL +IG+DLG+++KGFK  M     K+ +T  DA    +T      D N ++    D+ +",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|675457         0 MGGLSIWHWLIVLLIVALVFGTKKLRNIGNDLGSAVKGFKDGM-----KEGETPADAQQL
                  0 |||.|||............||||||..||.|||...||||..|-----|...|..||...
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|675457        55 PRT---GTVDVNAKETTRSDSNK 75
                 60 ..|---...|.|.......|... 83
gi|491764        74 AKTIADKQADTNQEQAKTEDAKR 97
""",
        )
        hit = hits[47]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|45435806|gb|AAS61363.1|")
        self.assertEqual(hit.target.name, "AAS61363")
        self.assertEqual(
            hit.target.description,
            "sec-independent protein translocase protein [Yersinia pestis biovar Medievalis str. 91001] >gi|51595437|ref|YP_069628.1| sec-independent protein translocase protein [Yersinia pseudotuberculosis IP 32953] >gi|22125073|ref|NP_668496.1| hypothetical protein y1170 [Yersinia pestis KIM] >gi|77636535|ref|ZP_00798608.1| COG1826: Sec-independent protein secretion pathway components [Yersinia pestis Angola] >gi|77632039|ref|ZP_00794625.1| COG1826: Sec-independent protein secretion pathway components [Yersinia pseudotuberculosis IP 31758] >gi|51588719|emb|CAH20330.1| sec-independent protein translocase protein [Yersinia pseudotuberculosis IP 32953] >gi|15980584|emb|CAC92840.1| sec-independent protein translocase protein [Yersinia pestis CO92] >gi|21957926|gb|AAM84747.1| hypothetical protein [Yersinia pestis KIM] >gi|45440947|ref|NP_992486.1| sec-independent protein translocase protein [Yersinia pestis biovar Medievalis str. 91001] >gi|16122810|ref|NP_406123.1| sec-independent protein translocase protein [Yersinia pestis CO92] >gi|25302680|pir||AI0316 sec-independent protein translocase protein [imported] - Yersinia pestis (strain CO92) >gi|24212499|sp|Q8ZDH1|TATE_YERPE Sec-independent protein translocase protein tatE",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=85)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 125.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 52.7582)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.25779e-06)
        self.assertEqual(hsp.annotations["identity"], 34)
        self.assertEqual(hsp.annotations["positive"], 58)
        self.assertEqual(hsp.annotations["gaps"], 2)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 61, 63, 84],
                          [14, 75, 75, 96]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 84))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...DAK'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MEGLSITKLLVVGILIVLLFGTSKLRTLGADLGAALKGFKKAMRNDDEVSTSVL...ERK'}, length=85)",
        )
        self.assertEqual(hsp.target.id, "gi|45435806|gb|AAS61363.1|")
        self.assertEqual(hsp.target.name, "AAS61363")
        self.assertEqual(
            hsp.target.description,
            "sec-independent protein translocase protein [Yersinia pestis biovar Medievalis str. 91001] >gi|51595437|ref|YP_069628.1| sec-independent protein translocase protein [Yersinia pseudotuberculosis IP 32953] >gi|22125073|ref|NP_668496.1| hypothetical protein y1170 [Yersinia pestis KIM] >gi|77636535|ref|ZP_00798608.1| COG1826: Sec-independent protein secretion pathway components [Yersinia pestis Angola] >gi|77632039|ref|ZP_00794625.1| COG1826: Sec-independent protein secretion pathway components [Yersinia pseudotuberculosis IP 31758] >gi|51588719|emb|CAH20330.1| sec-independent protein translocase protein [Yersinia pseudotuberculosis IP 32953] >gi|15980584|emb|CAC92840.1| sec-independent protein translocase protein [Yersinia pestis CO92] >gi|21957926|gb|AAM84747.1| hypothetical protein [Yersinia pestis KIM] >gi|45440947|ref|NP_992486.1| sec-independent protein translocase protein [Yersinia pestis biovar Medievalis str. 91001] >gi|16122810|ref|NP_406123.1| sec-independent protein translocase protein [Yersinia pestis CO92] >gi|25302680|pir||AI0316 sec-independent protein translocase protein [imported] - Yersinia pestis (strain CO92) >gi|24212499|sp|Q8ZDH1|TATE_YERPE Sec-independent protein translocase protein tatE",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "M G+SI +LL++ +++VLLFGT KL ++G+DLGA++KGFKKAM +D+        +   +A  KT+A+ +A ++ + A + + K",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|454358         0 MEGLSITKLLVVGILIVLLFGTSKLRTLGADLGAALKGFKKAMRNDDEVSTSVLGETKMS
                  0 |.|.||.............|||.||...|.||||..|||||||..|..............
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|454358        60 AETKTVAETKAASDSQAAASVERK 84
                 60 |--||.|...|......|.....| 84
gi|491764        74 A--KTIADKQADTNQEQAKTEDAK 96
""",
        )
        hit = hits[48]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|49610761|emb|CAG74206.1|")
        self.assertEqual(hit.target.name, "CAG74206")
        self.assertEqual(
            hit.target.description,
            "Sec-independent protein translocase protein [Erwinia carotovora subsp. atroseptica SCRI1043] >gi|50120235|ref|YP_049402.1| Sec-independent protein translocase protein [Erwinia carotovora subsp. atroseptica SCRI1043]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=65)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 125.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 52.7582)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.25779e-06)
        self.assertEqual(hsp.annotations["identity"], 37)
        self.assertEqual(hsp.annotations["positive"], 47)
        self.assertEqual(hsp.annotations["gaps"], 1)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 48, 49, 63],
                          [14, 62, 62, 76]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 63))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...TAK'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MEGISIAKLLVIGALIVLLFGTNKLRSLGGDLGAAIKGFKKAMNDDQSVKTDDT...SRK'}, length=65)",
        )
        self.assertEqual(hsp.target.id, "gi|49610761|emb|CAG74206.1|")
        self.assertEqual(hsp.target.name, "CAG74206")
        self.assertEqual(
            hsp.target.description,
            "Sec-independent protein translocase protein [Erwinia carotovora subsp. atroseptica SCRI1043] >gi|50120235|ref|YP_049402.1| Sec-independent protein translocase protein [Erwinia carotovora subsp. atroseptica SCRI1043]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "M GISI +LL+I  ++VLLFGT KL S+G DLGA+IKGFKKAM+DD+  K D T+   D + K",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|496107         0 MEGISIAKLLVIGALIVLLFGTNKLRSLGGDLGAAIKGFKKAMNDDQSVKTDDTAALNDS
                  0 |.||||.............|||.||.|.|.||||.||||||||.||..-|.|.|....|.
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEP-KQDKTSQDADF

gi|496107        60 SRK 63
                 60 ..| 63
gi|491764        73 TAK 76
""",
        )
        hit = hits[49]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|67663266|ref|ZP_00460549.1|")
        self.assertEqual(hit.target.name, "ZP_00460549")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Burkholderia cenocepacia HI2424] >gi|67659452|ref|ZP_00456814.1| Twin-arginine translocation protein TatA/E [Burkholderia cenocepacia AU 1054] >gi|67103109|gb|EAM20236.1| Twin-arginine translocation protein TatA/E [Burkholderia cenocepacia HI2424] >gi|67093009|gb|EAM10556.1| Twin-arginine translocation protein TatA/E [Burkholderia cenocepacia AU 1054]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=76)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 124.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 52.373)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.25481e-06)
        self.assertEqual(hsp.annotations["identity"], 29)
        self.assertEqual(hsp.annotations["positive"], 39)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 51],
                          [14, 65]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 51))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQD'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGGLSIWHWLIVLLIVALVFGTKKLRNIGNDLGSAVKGFKDGMKEGETPAD'}, length=76)",
        )
        self.assertEqual(hsp.target.id, "gi|67663266|ref|ZP_00460549.1|")
        self.assertEqual(hsp.target.name, "ZP_00460549")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Burkholderia cenocepacia HI2424] >gi|67659452|ref|ZP_00456814.1| Twin-arginine translocation protein TatA/E [Burkholderia cenocepacia AU 1054] >gi|67103109|gb|EAM20236.1| Twin-arginine translocation protein TatA/E [Burkholderia cenocepacia HI2424] >gi|67093009|gb|EAM10556.1| Twin-arginine translocation protein TatA/E [Burkholderia cenocepacia AU 1054]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGG+SIW  LI+ +IV L+FGTKKL +IG+DLG+++KGFK  M + E   D",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|676632         0 MGGLSIWHWLIVLLIVALVFGTKKLRNIGNDLGSAVKGFKDGMKEGETPAD 51
                  0 |||.|||............||||||..||.|||...||||..|...|...| 51
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQD 65
""",
        )
        hit = hits[50]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|33594634|ref|NP_882278.1|")
        self.assertEqual(hit.target.name, "NP_882278")
        self.assertEqual(
            hit.target.description,
            "Sec-independent protein translocase protein [Bordetella pertussis Tohama I] >gi|33598762|ref|NP_886405.1| Sec-independent protein translocase protein [Bordetella parapertussis 12822] >gi|33603836|ref|NP_891396.1| Sec-independent protein translocase protein [Bordetella bronchiseptica RB50] >gi|33577961|emb|CAE35226.1| Sec-independent protein translocase protein [Bordetella bronchiseptica RB50] >gi|33574892|emb|CAE39555.1| Sec-independent protein translocase protein [Bordetella parapertussis] >gi|33564710|emb|CAE44033.1| Sec-independent protein translocase protein [Bordetella pertussis Tohama I]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=75)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 124.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 52.373)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.25481e-06)
        self.assertEqual(hsp.annotations["identity"], 31)
        self.assertEqual(hsp.annotations["positive"], 41)
        self.assertEqual(hsp.annotations["gaps"], 3)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 45, 48, 58],
                          [14, 59, 59, 69]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 58))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGSFSIWHWLIVLVIIALVFGTKKLRNVGSDLGSAVKGFKEGMKDASADKPADQVTQQ'}, length=75)",
        )
        self.assertEqual(hsp.target.id, "gi|33594634|ref|NP_882278.1|")
        self.assertEqual(hsp.target.name, "NP_882278")
        self.assertEqual(
            hsp.target.description,
            "Sec-independent protein translocase protein [Bordetella pertussis Tohama I] >gi|33598762|ref|NP_886405.1| Sec-independent protein translocase protein [Bordetella parapertussis 12822] >gi|33603836|ref|NP_891396.1| Sec-independent protein translocase protein [Bordetella bronchiseptica RB50] >gi|33577961|emb|CAE35226.1| Sec-independent protein translocase protein [Bordetella bronchiseptica RB50] >gi|33574892|emb|CAE39555.1| Sec-independent protein translocase protein [Bordetella parapertussis] >gi|33564710|emb|CAE44033.1| Sec-independent protein translocase protein [Bordetella pertussis Tohama I]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG  SIW  LI+ VI+ L+FGTKKL ++GSDLG+++KGFK+ M D   D+P    T Q",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|335946         0 MGSFSIWHWLIVLVIIALVFGTKKLRNVGSDLGSAVKGFKEGMKDASADKPADQVTQQ 58
                  0 ||..|||............||||||...|||||...||||..|.|---|.|....|.| 58
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSD---DEPKQDKTSQ 69
""",
        )
        hit = hits[51]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|46310681|ref|ZP_00211309.1|")
        self.assertEqual(hit.target.name, "ZP_00211309")
        self.assertEqual(
            hit.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Burkholderia cepacia R18194]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=76)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 124.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 52.373)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.25481e-06)
        self.assertEqual(hsp.annotations["identity"], 29)
        self.assertEqual(hsp.annotations["positive"], 39)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 51],
                          [14, 65]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 51))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQD'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGGLSIWHWLIVLLIVALVFGTKKLRNIGNDLGSAVKGFKDGMKEGEAPAD'}, length=76)",
        )
        self.assertEqual(hsp.target.id, "gi|46310681|ref|ZP_00211309.1|")
        self.assertEqual(hsp.target.name, "ZP_00211309")
        self.assertEqual(
            hsp.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Burkholderia cepacia R18194]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGG+SIW  LI+ +IV L+FGTKKL +IG+DLG+++KGFK  M + E   D",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|463106         0 MGGLSIWHWLIVLLIVALVFGTKKLRNIGNDLGSAVKGFKDGMKEGEAPAD 51
                  0 |||.|||............||||||..||.|||...||||..|...|...| 51
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQD 65
""",
        )
        hit = hits[52]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|58584031|ref|YP_203047.1|")
        self.assertEqual(hit.target.name, "YP_203047")
        self.assertEqual(
            hit.target.description,
            "sec-independent protein translocase [Xanthomonas oryzae pv. oryzae KACC10331] >gi|21244935|ref|NP_644517.1| sec-independent protein translocase [Xanthomonas axonopodis pv. citri str. 306] >gi|21110652|gb|AAM39053.1| sec-independent protein translocase [Xanthomonas axonopodis pv. citri str. 306] >gi|58428625|gb|AAW77662.1| sec-independent protein translocase [Xanthomonas oryzae pv. oryzae KACC10331] >gi|24212477|sp|Q8PEX2|TATA_XANAC Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=75)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 123.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 51.9878)
        self.assertAlmostEqual(hsp.annotations["evalue"], 5.55696e-06)
        self.assertEqual(hsp.annotations["identity"], 34)
        self.assertEqual(hsp.annotations["positive"], 49)
        self.assertEqual(hsp.annotations["gaps"], 1)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 58, 58, 73],
                          [14, 72, 73, 88]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 74))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...NQE'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGGFSIWHWLIVLVIVLLVFGTKRLTSGAKDLGSAVKEFKKGMHDDDKPAGKLG...DRD'}, length=75)",
        )
        self.assertEqual(hsp.target.id, "gi|58584031|ref|YP_203047.1|")
        self.assertEqual(hsp.target.name, "YP_203047")
        self.assertEqual(
            hsp.target.description,
            "sec-independent protein translocase [Xanthomonas oryzae pv. oryzae KACC10331] >gi|21244935|ref|NP_644517.1| sec-independent protein translocase [Xanthomonas axonopodis pv. citri str. 306] >gi|21110652|gb|AAM39053.1| sec-independent protein translocase [Xanthomonas axonopodis pv. citri str. 306] >gi|58428625|gb|AAW77662.1| sec-independent protein translocase [Xanthomonas oryzae pv. oryzae KACC10331] >gi|24212477|sp|Q8PEX2|TATA_XANAC Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGG SIW  LI+ VIV+L+FGTK+L S   DLG+++K FKK M DD+    K   D+  TA+   + QA+ +++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|585840         0 MGGFSIWHWLIVLVIVLLVFGTKRLTSGAKDLGSAVKEFKKGMHDDDKPAGKLGDDSR-T
                  0 |||.|||............||||.|.|...|||...|.|||.|.||.....|...|..-|
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|585840        59 AEQAREAQAERDRD 73
                 60 |......||..... 74
gi|491764        74 AKTIADKQADTNQE 88
""",
        )
        hit = hits[53]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|17429965|emb|CAD16649.1|")
        self.assertEqual(hit.target.name, "CAD16649")
        self.assertEqual(
            hit.target.description,
            "PROBABLE SIGNAL PEPTIDE PROTEIN [Ralstonia solanacearum] >gi|17547661|ref|NP_521063.1| PROBABLE SIGNAL PEPTIDE PROTEIN [Ralstonia solanacearum GMI1000] >gi|24212492|sp|Q8XV89|TATA_RALSO Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=85)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 123.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 51.9878)
        self.assertAlmostEqual(hsp.annotations["evalue"], 5.55696e-06)
        self.assertEqual(hsp.annotations["identity"], 30)
        self.assertEqual(hsp.annotations["positive"], 44)
        self.assertEqual(hsp.annotations["gaps"], 3)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 43, 46, 60],
                          [14, 57, 57, 71]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 60))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDA'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGSFSIWHWLIVLLIIMMVFGTKKLRNIGSDLGSAVKGFKEGMREGSEDKPAGSQQGQQA'}, length=85)",
        )
        self.assertEqual(hsp.target.id, "gi|17429965|emb|CAD16649.1|")
        self.assertEqual(hsp.target.name, "CAD16649")
        self.assertEqual(
            hsp.target.description,
            "PROBABLE SIGNAL PEPTIDE PROTEIN [Ralstonia solanacearum] >gi|17547661|ref|NP_521063.1| PROBABLE SIGNAL PEPTIDE PROTEIN [Ralstonia solanacearum GMI1000] >gi|24212492|sp|Q8XV89|TATA_RALSO Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG  SIW  LI+ +I++++FGTKKL +IGSDLG+++KGFK+ M   S+D+P   +  Q A",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|174299         0 MGSFSIWHWLIVLLIIMMVFGTKKLRNIGSDLGSAVKGFKEGMREGSEDKPAGSQQGQQA
                  0 ||..|||............||||||..||||||...||||..|---|.|.|......|.|
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAM---SDDEPKQDKTSQDA

gi|174299        60 
                 60 
gi|491764        71 
""",
        )
        hit = hits[54]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|47573371|ref|ZP_00243410.1|")
        self.assertEqual(hit.target.name, "ZP_00243410")
        self.assertEqual(
            hit.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Rubrivivax gelatinosus PM1]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=76)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 123.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 51.9878)
        self.assertAlmostEqual(hsp.annotations["evalue"], 5.55696e-06)
        self.assertEqual(hsp.annotations["identity"], 35)
        self.assertEqual(hsp.annotations["positive"], 47)
        self.assertEqual(hsp.annotations["gaps"], 6)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 43, 43, 66],
                          [14, 57, 63, 86]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 72))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...DTN'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGSFSIWHWLIVLLIVVLVFGTKKLKNIGSDLGGAVKGFKDGVRDGSTAPADPA...DAN'}, length=76)",
        )
        self.assertEqual(hsp.target.id, "gi|47573371|ref|ZP_00243410.1|")
        self.assertEqual(hsp.target.name, "ZP_00243410")
        self.assertEqual(
            hsp.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Rubrivivax gelatinosus PM1]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG  SIW  LI+ +IVVL+FGTKKL +IGSDLG ++KGFK  +      +D ++  AD   +  A+K AD N",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|475733         0 MGSFSIWHWLIVLLIVVLVFGTKKLKNIGSDLGGAVKGFKDGV------RDGSTAPADPA
                  0 ||..|||............||||||..||||||...||||...------.|.....||..
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|475733        54 QQVTANKSADAN 66
                 60 ....|.|.||.| 72
gi|491764        74 AKTIADKQADTN 86
""",
        )
        hit = hits[55]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|16273687|ref|NP_438355.1|")
        self.assertEqual(hit.target.name, "NP_438355")
        self.assertEqual(
            hit.target.description,
            "Sec-independent protein secretion pathway component TatA [Haemophilus influenzae Rd KW20] >gi|9988065|sp|P57046|TATA_HAEIN Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=89)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 122.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 51.6026)
        self.assertAlmostEqual(hsp.annotations["evalue"], 7.25761e-06)
        self.assertEqual(hsp.annotations["identity"], 37)
        self.assertEqual(hsp.annotations["positive"], 56)
        self.assertEqual(hsp.annotations["gaps"], 7)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 43, 43, 54, 54, 72],
                          [14, 57, 62, 73, 75, 93]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 79))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...KTE'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MFGLSPAQLIILLVVILLIFGTKKLRNAGSDLGAAVKGFKKAMKEDEKVKDAEF...KRE'}, length=89)",
        )
        self.assertEqual(hsp.target.id, "gi|16273687|ref|NP_438355.1|")
        self.assertEqual(hsp.target.name, "NP_438355")
        self.assertEqual(
            hsp.target.description,
            "Sec-independent protein secretion pathway component TatA [Haemophilus influenzae Rd KW20] >gi|9988065|sp|P57046|TATA_HAEIN Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "M G+S  QL+I+ V+++L+FGTKKL + GSDLGA++KGFKKAM     K+D+  +DA+F  K+I ++ A   + + K E",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|162736         0 MFGLSPAQLIILLVVILLIFGTKKLRNAGSDLGAAVKGFKKAM-----KEDEKVKDAEF-
                  0 |.|.|..|...........||||||...||||||..|||||||-----|.|....||.|-
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|162736        54 -KSIDNETASAKKGKYKRE 72
                 60 -|.|....|.......|.| 79
gi|491764        74 AKTIADKQADTNQEQAKTE 93
""",
        )
        hit = hits[56]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|73542784|ref|YP_297304.1|")
        self.assertEqual(hit.target.name, "YP_297304")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Ralstonia eutropha JMP134] >gi|72120197|gb|AAZ62460.1| Twin-arginine translocation protein TatA/E [Ralstonia eutropha JMP134] >gi|53760639|ref|ZP_00165809.2| COG1826: Sec-independent protein secretion pathway components [Ralstonia eutropha JMP134]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=73)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 121.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 51.2174)
        self.assertAlmostEqual(hsp.annotations["evalue"], 9.47873e-06)
        self.assertEqual(hsp.annotations["identity"], 34)
        self.assertEqual(hsp.annotations["positive"], 44)
        self.assertEqual(hsp.annotations["gaps"], 7)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 49, 56, 70],
                          [14, 63, 63, 77]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 70))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...AKT'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGSFSIWHWLIVLVIVMLVFGTKKLRNIGQDLGGAVKGFKDGMKEGDDKAAPAK...EKT'}, length=73)",
        )
        self.assertEqual(hsp.target.id, "gi|73542784|ref|YP_297304.1|")
        self.assertEqual(hsp.target.name, "YP_297304")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Ralstonia eutropha JMP134] >gi|72120197|gb|AAZ62460.1| Twin-arginine translocation protein TatA/E [Ralstonia eutropha JMP134] >gi|53760639|ref|ZP_00165809.2| COG1826: Sec-independent protein secretion pathway components [Ralstonia eutropha JMP134]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG  SIW  LI+ VIV+L+FGTKKL +IG DLG ++KGFK  M + + K       +D T+ D D   KT",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|735427         0 MGSFSIWHWLIVLVIVMLVFGTKKLRNIGQDLGGAVKGFKDGMKEGDDKAAPAKELRDST
                  0 ||..|||............||||||..||.|||...||||..|.....|-------.|.|
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPK-------QDKT

gi|735427        60 TIDVDAKEKT 70
                 60 ..|.|...|| 70
gi|491764        67 SQDADFTAKT 77
""",
        )
        hit = hits[57]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|26987777|ref|NP_743202.1|")
        self.assertEqual(hit.target.name, "NP_743202")
        self.assertEqual(
            hit.target.description,
            "Sec-independent protein translocase TatA [Pseudomonas putida KT2440] >gi|24982471|gb|AAN66666.1| Sec-independent protein translocase TatA [Pseudomonas putida KT2440]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=77)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 121.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 51.2174)
        self.assertAlmostEqual(hsp.annotations["evalue"], 9.47873e-06)
        self.assertEqual(hsp.annotations["identity"], 35)
        self.assertEqual(hsp.annotations["positive"], 49)
        self.assertEqual(hsp.annotations["gaps"], 8)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 46, 46, 75],
                          [14, 60, 68, 97]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 83))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...AKR'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGGIGIWQLVIVLLIVFLLFGTKRLKGLGSDVGEAIQGFRKSMGGDNDASAADQ...ADR'}, length=77)",
        )
        self.assertEqual(hsp.target.id, "gi|26987777|ref|NP_743202.1|")
        self.assertEqual(hsp.target.name, "NP_743202")
        self.assertEqual(
            hsp.target.description,
            "Sec-independent protein translocase TatA [Pseudomonas putida KT2440] >gi|24982471|gb|AAN66666.1| Sec-independent protein translocase TatA [Pseudomonas putida KT2440]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGGI IWQL+I+ +IV LLFGTK+L  +GSD+G +I+GF+K+M  D         DA    +    +Q   N + A+   A R",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|269877         0 MGGIGIWQLVIVLLIVFLLFGTKRLKGLGSDVGEAIQGFRKSMGGD--------NDASAA
                  0 ||||.|||...........||||.|...|||.|..|.||.|.|..|--------.||...
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|269877        52 DQAHVQQQCPLNGQVAQQSQADR 75
                 60 .......|...|...|....|.| 83
gi|491764        74 AKTIADKQADTNQEQAKTEDAKR 97
""",
        )
        hit = hits[58]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|29142636|ref|NP_805978.1|")
        self.assertEqual(hit.target.name, "NP_805978")
        self.assertEqual(
            hit.target.description,
            "sec-independent protein translocase protein TatE [Salmonella enterica subsp. enterica serovar Typhi Ty2] >gi|16501883|emb|CAD05109.1| sec-independent protein translocase protein TatE [Salmonella enterica subsp. enterica serovar Typhi] >gi|16419144|gb|AAL19583.1| putative Sec-independent protein secretion pathway component [Salmonella typhimurium LT2] >gi|62179231|ref|YP_215648.1| putative Sec-independent protein secretion pathway component [Salmonella enterica subsp. enterica serovar Choleraesuis str. SC-B67] >gi|16759591|ref|NP_455208.1| sec-independent protein translocase protein TatE [Salmonella enterica subsp. enterica serovar Typhi str. CT18] >gi|56414233|ref|YP_151308.1| sec-independent protein translocase protein TatE [Salmonella enterica subsp. enterica serovar Paratyphi A str. ATCC 9150] >gi|29138267|gb|AAO69838.1| sec-independent protein translocase protein TatE [Salmonella enterica subsp. enterica serovar Typhi Ty2] >gi|56128490|gb|AAV77996.1| sec-independent protein translocase protein TatE [Salmonella enterica subsp. enterica serovar Paratyphi A str. ATCC 9150] >gi|62126864|gb|AAX64567.1| putative Sec-independent protein secretion pathway component [Salmonella enterica subsp. enterica serovar Choleraesuis str. SC-B67] >gi|61248019|sp|P0A2H5|TATE_SALTY Sec-independent protein translocase protein tatE >gi|61248020|sp|P0A2H6|TATE_SALTI Sec-independent protein translocase protein tatE >gi|25302676|pir||AC0580 sec-independent protein translocase protein TatE [imported] - Salmonella enterica subsp. enterica serovar Typhi (strain CT18) >gi|16764009|ref|NP_459624.1| Sec-independent protein secretion pathway component [Salmonella typhimurium LT2]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=67)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 119.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 50.447)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.61683e-05)
        self.assertEqual(hsp.annotations["identity"], 35)
        self.assertEqual(hsp.annotations["positive"], 49)
        self.assertEqual(hsp.annotations["gaps"], 1)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 50, 50, 67],
                          [14, 64, 65, 82]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 68))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...DKQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MGEISITKLLVVAALVVLLFGTKKLRTLGGDLGTAIKGFKKAMNDEDAGVKKDV...HKE')",
        )
        self.assertEqual(hsp.target.id, "gi|29142636|ref|NP_805978.1|")
        self.assertEqual(hsp.target.name, "NP_805978")
        self.assertEqual(
            hsp.target.description,
            "sec-independent protein translocase protein TatE [Salmonella enterica subsp. enterica serovar Typhi Ty2] >gi|16501883|emb|CAD05109.1| sec-independent protein translocase protein TatE [Salmonella enterica subsp. enterica serovar Typhi] >gi|16419144|gb|AAL19583.1| putative Sec-independent protein secretion pathway component [Salmonella typhimurium LT2] >gi|62179231|ref|YP_215648.1| putative Sec-independent protein secretion pathway component [Salmonella enterica subsp. enterica serovar Choleraesuis str. SC-B67] >gi|16759591|ref|NP_455208.1| sec-independent protein translocase protein TatE [Salmonella enterica subsp. enterica serovar Typhi str. CT18] >gi|56414233|ref|YP_151308.1| sec-independent protein translocase protein TatE [Salmonella enterica subsp. enterica serovar Paratyphi A str. ATCC 9150] >gi|29138267|gb|AAO69838.1| sec-independent protein translocase protein TatE [Salmonella enterica subsp. enterica serovar Typhi Ty2] >gi|56128490|gb|AAV77996.1| sec-independent protein translocase protein TatE [Salmonella enterica subsp. enterica serovar Paratyphi A str. ATCC 9150] >gi|62126864|gb|AAX64567.1| putative Sec-independent protein secretion pathway component [Salmonella enterica subsp. enterica serovar Choleraesuis str. SC-B67] >gi|61248019|sp|P0A2H5|TATE_SALTY Sec-independent protein translocase protein tatE >gi|61248020|sp|P0A2H6|TATE_SALTI Sec-independent protein translocase protein tatE >gi|25302676|pir||AC0580 sec-independent protein translocase protein TatE [imported] - Salmonella enterica subsp. enterica serovar Typhi (strain CT18) >gi|16764009|ref|NP_459624.1| Sec-independent protein secretion pathway component [Salmonella typhimurium LT2]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG ISI +LL++A +VVLLFGTKKL ++G DLG +IKGFKKAM+D++    K   D    A+ ++ K+",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|291426         0 MGEISITKLLVVAALVVLLFGTKKLRTLGGDLGTAIKGFKKAMNDEDAGV-KKDVDGSVQ
                  0 ||.|||.............||||||...|.|||..||||||||.|.....-|...|....
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|291426        59 AEKLSHKE 67
                 60 |.....|. 68
gi|491764        74 AKTIADKQ 82
""",
        )
        hit = hits[59]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|18389921|gb|AAL68797.1|")
        self.assertEqual(hit.target.name, "AAL68797")
        self.assertEqual(hit.target.description, "TatA [Ralstonia eutropha]")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=77)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 119.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 50.447)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.61683e-05)
        self.assertEqual(hsp.annotations["identity"], 35)
        self.assertEqual(hsp.annotations["positive"], 44)
        self.assertEqual(hsp.annotations["gaps"], 11)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 49, 60, 74],
                          [14, 63, 63, 77]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 74))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...AKT'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGSFSIWHWLIVLVIVMLVFGTKKLRNIGQDLGGAVKGFKDGMKDGEGKAAADP...EKT'}, length=77)",
        )
        self.assertEqual(hsp.target.id, "gi|18389921|gb|AAL68797.1|")
        self.assertEqual(hsp.target.name, "AAL68797")
        self.assertEqual(hsp.target.description, "TatA [Ralstonia eutropha]")
        self.assertEqual(
            hsp.annotations["midline"],
            "MG  SIW  LI+ VIV+L+FGTKKL +IG DLG ++KGFK  M D E K           +D T+ D +   KT",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|183899         0 MGSFSIWHWLIVLVIVMLVFGTKKLRNIGQDLGGAVKGFKDGMKDGEGKAAADPAQSKEL
                  0 ||..|||............||||||..||.|||...||||..|.|.|.|-----------
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPK-----------

gi|183899        60 RDSTTIDVEAKEKT 74
                 60 .|.|..|.....|| 74
gi|491764        63 QDKTSQDADFTAKT 77
""",
        )
        hit = hits[60]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|48781637|ref|ZP_00278228.1|")
        self.assertEqual(hit.target.name, "ZP_00278228")
        self.assertEqual(
            hit.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Burkholderia fungorum LB400]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=79)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 117.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 49.6766)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.75789e-05)
        self.assertEqual(hsp.annotations["identity"], 28)
        self.assertEqual(hsp.annotations["positive"], 39)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 55],
                          [14, 69]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 55))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGSLSIWHWLIVLLIVALVFGTKKLRNIGTDLGGAVKGFKEGMKEAETPAGEAQQ'}, length=79)",
        )
        self.assertEqual(hsp.target.id, "gi|48781637|ref|ZP_00278228.1|")
        self.assertEqual(hsp.target.name, "ZP_00278228")
        self.assertEqual(
            hsp.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Burkholderia fungorum LB400]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG +SIW  LI+ +IV L+FGTKKL +IG+DLG ++KGFK+ M + E    +  Q",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|487816         0 MGSLSIWHWLIVLLIVALVFGTKKLRNIGTDLGGAVKGFKEGMKEAETPAGEAQQ 55
                  0 ||..|||............||||||..||.|||...||||..|...|.......| 55
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQ 69
""",
        )
        hit = hits[61]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|77456610|ref|YP_346115.1|")
        self.assertEqual(hit.target.name, "YP_346115")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Pseudomonas fluorescens PfO-1] >gi|77380613|gb|ABA72126.1| Twin-arginine translocation protein TatA/E [Pseudomonas fluorescens PfO-1]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=92)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 117.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 49.6766)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.75789e-05)
        self.assertEqual(hsp.annotations["identity"], 22)
        self.assertEqual(hsp.annotations["positive"], 28)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[19, 53],
                          [33, 67]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 34))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKT'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({19: 'FGTKKLKNLGTDVGESIKGFRKAMNDDEKPADPT'}, length=92)",
        )
        self.assertEqual(hsp.target.id, "gi|77456610|ref|YP_346115.1|")
        self.assertEqual(hsp.target.name, "YP_346115")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Pseudomonas fluorescens PfO-1] >gi|77380613|gb|ABA72126.1| Twin-arginine translocation protein TatA/E [Pseudomonas fluorescens PfO-1]",
        )
        self.assertEqual(
            hsp.annotations["midline"], "FGTKKL ++G+D+G SIKGF+KAM+DDE   D T"
        )
        self.assertEqual(
            str(hsp),
            """\
gi|774566        19 FGTKKLKNLGTDVGESIKGFRKAMNDDEKPADPT 53
                  0 ||||||...|.|.|.|||||.|||.|||...|.| 34
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKT 67
""",
        )
        hit = hits[62]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|1684735|emb|CAA98158.1|")
        self.assertEqual(hit.target.name, "CAA98158")
        self.assertEqual(
            hit.target.description,
            "ORF57 protein [Pseudomonas stutzeri] >gi|9979022|sp|P95557|TATA_PSEST Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=57)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 117.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 49.6766)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.75789e-05)
        self.assertEqual(hsp.annotations["identity"], 30)
        self.assertEqual(hsp.annotations["positive"], 40)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 1, 46],
                          [16, 61]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 45))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({16: 'GISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDE'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({1: 'GISVWQLLIILLIVVMLFGTKRLRGLGSDLGSAINGFRKSVSDGE'}, length=57)",
        )
        self.assertEqual(hsp.target.id, "gi|1684735|emb|CAA98158.1|")
        self.assertEqual(hsp.target.name, "CAA98158")
        self.assertEqual(
            hsp.target.description,
            "ORF57 protein [Pseudomonas stutzeri] >gi|9979022|sp|P95557|TATA_PSEST Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(
            hsp.annotations["midline"], "GIS+WQLLII +IVV+LFGTK+L  +GSDLG++I GF+K++SD E"
        )
        self.assertEqual(
            str(hsp),
            """\
gi|168473         1 GISVWQLLIILLIVVMLFGTKRLRGLGSDLGSAINGFRKSVSDGE 46
                  0 |||.||...........||||.|...|||||..|.||.|..||.| 45
gi|491764        16 GISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDE 61
""",
        )
        hit = hits[63]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|56476124|ref|YP_157713.1|")
        self.assertEqual(hit.target.name, "YP_157713")
        self.assertEqual(
            hit.target.description,
            "Sec-independent protein translocase subunit A [Azoarcus sp. EbN1] >gi|56312167|emb|CAI06812.1| Sec-independent protein translocase subunit A [Azoarcus sp. EbN1]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=75)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 116.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 49.2914)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.60191e-05)
        self.assertEqual(hsp.annotations["identity"], 30)
        self.assertEqual(hsp.annotations["positive"], 37)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 55],
                          [14, 69]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 55))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGSFSIWHWLIVLVIVVLVFGTKKLRNVGQDLGGAVKGFKDGMRDSEKSGEDVQQ'}, length=75)",
        )
        self.assertEqual(hsp.target.id, "gi|56476124|ref|YP_157713.1|")
        self.assertEqual(hsp.target.name, "YP_157713")
        self.assertEqual(
            hsp.target.description,
            "Sec-independent protein translocase subunit A [Azoarcus sp. EbN1] >gi|56312167|emb|CAI06812.1| Sec-independent protein translocase subunit A [Azoarcus sp. EbN1]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG  SIW  LI+ VIVVL+FGTKKL ++G DLG ++KGFK  M D E   +   Q",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|564761         0 MGSFSIWHWLIVLVIVVLVFGTKKLRNVGQDLGGAVKGFKDGMRDSEKSGEDVQQ 55
                  0 ||..|||............||||||...|.|||...||||..|.|.|.......| 55
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQ 69
""",
        )
        hit = hits[64]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|34496078|ref|NP_900293.1|")
        self.assertEqual(hit.target.name, "NP_900293")
        self.assertEqual(
            hit.target.description,
            "Sec-independent protein translocase protein TatA [Chromobacterium violaceum ATCC 12472] >gi|34101932|gb|AAQ58299.1| Sec-independent protein translocase protein TatA [Chromobacterium violaceum ATCC 12472]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=68)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 116.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 49.2914)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.60191e-05)
        self.assertEqual(hsp.annotations["identity"], 32)
        self.assertEqual(hsp.annotations["positive"], 44)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 68],
                          [14, 82]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 68))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...DKQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MGSLSIWHWLIVLLIVVLVFGTKKLPGIGKDLGNAVKGFKEGMNEGAKDGQPPA...DKK')",
        )
        self.assertEqual(hsp.target.id, "gi|34496078|ref|NP_900293.1|")
        self.assertEqual(hsp.target.name, "NP_900293")
        self.assertEqual(
            hsp.target.description,
            "Sec-independent protein translocase protein TatA [Chromobacterium violaceum ATCC 12472] >gi|34101932|gb|AAQ58299.1| Sec-independent protein translocase protein TatA [Chromobacterium violaceum ATCC 12472]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG +SIW  LI+ +IVVL+FGTKKL  IG DLG ++KGFK+ M++        ++DA       ADK+",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|344960         0 MGSLSIWHWLIVLLIVVLVFGTKKLPGIGKDLGNAVKGFKEGMNEGAKDGQPPAKDAGRI
                  0 ||..|||............||||||..||.|||...||||..|............||...
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|344960        60 IDGEADKK 68
                 60 ....|||. 68
gi|491764        74 AKTIADKQ 82
""",
        )
        hit = hits[65]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|67848115|ref|ZP_00503233.1|")
        self.assertEqual(hit.target.name, "ZP_00503233")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Polaromonas sp. JS666] >gi|67781162|gb|EAM40776.1| Twin-arginine translocation protein TatA/E [Polaromonas sp. JS666]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=83)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 115.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 48.9062)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.70425e-05)
        self.assertEqual(hsp.annotations["identity"], 31)
        self.assertEqual(hsp.annotations["positive"], 42)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 69],
                          [14, 83]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 69))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...KQA'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGSFSIWHWLIVLLIVVMVFGTKKLKNMGSDLGSAVKGFKDGMKDGGQSAAATD...AQA'}, length=83)",
        )
        self.assertEqual(hsp.target.id, "gi|67848115|ref|ZP_00503233.1|")
        self.assertEqual(hsp.target.name, "ZP_00503233")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Polaromonas sp. JS666] >gi|67781162|gb|EAM40776.1| Twin-arginine translocation protein TatA/E [Polaromonas sp. JS666]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG  SIW  LI+ +IVV++FGTKKL ++GSDLG+++KGFK  M D       T       A  + + QA",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|678481         0 MGSFSIWHWLIVLLIVVMVFGTKKLKNMGSDLGSAVKGFKDGMKDGGQSAAATDDKPAAP
                  0 ||..|||............||||||...|||||...||||..|.|.......|.......
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|678481        60 AGQVTNAQA 69
                 60 |......|| 69
gi|491764        74 AKTIADKQA 83
""",
        )
        hit = hits[66]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|26991692|ref|NP_747117.1|")
        self.assertEqual(hit.target.name, "NP_747117")
        self.assertEqual(
            hit.target.description,
            "Sec-independent protein translocase TatA [Pseudomonas putida KT2440] >gi|24986793|gb|AAN70581.1| Sec-independent protein translocase TatA [Pseudomonas putida KT2440]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=90)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 115.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 48.9062)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.70425e-05)
        self.assertEqual(hsp.annotations["identity"], 26)
        self.assertEqual(hsp.annotations["positive"], 40)
        self.assertEqual(hsp.annotations["gaps"], 3)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 19,  47,  47,  85],
                          [ 33,  61,  64, 102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 69))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQ...KEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({19: 'FGTKKLKNFGSDLGESIKGFRKAMNEEETKPAEQTPPPAQPVPPVQNTAQPQQG...VQE'}, length=90)",
        )
        self.assertEqual(hsp.target.id, "gi|26991692|ref|NP_747117.1|")
        self.assertEqual(hsp.target.name, "NP_747117")
        self.assertEqual(
            hsp.target.description,
            "Sec-independent protein translocase TatA [Pseudomonas putida KT2440] >gi|24986793|gb|AAN70581.1| Sec-independent protein translocase TatA [Pseudomonas putida KT2440]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "FGTKKL + GSDLG SIKGF+KAM+++E    K ++     A+ +   Q     +Q  T + + H  ++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|269916        19 FGTKKLKNFGSDLGESIKGFRKAMNEEE---TKPAEQTPPPAQPVPPVQNTAQPQQGHTI
                  0 ||||||...|||||.|||||.|||...|---.|........|......|......|..|.
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTE

gi|269916        76 EGQAHPVQE  85
                 60 ....|....  69
gi|491764        93 DAKRHDKEQ 102
""",
        )
        hit = hits[67]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|15601293|ref|NP_232924.1|")
        self.assertEqual(hit.target.name, "NP_232924")
        self.assertEqual(
            hit.target.description,
            "tatA protein [Vibrio cholerae O1 biovar eltor str. N16961] >gi|9657940|gb|AAF96436.1| tatA protein [Vibrio cholerae O1 biovar eltor str. N16961] >gi|11356195|pir||A82448 tatA protein VCA0533 [imported] - Vibrio cholerae (strain N16961 serogroup O1)",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=78)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 114.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 48.521)
        self.assertAlmostEqual(hsp.annotations["evalue"], 6.14393e-05)
        self.assertEqual(hsp.annotations["identity"], 38)
        self.assertEqual(hsp.annotations["positive"], 54)
        self.assertEqual(hsp.annotations["gaps"], 7)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 47, 47, 53, 53, 63, 63, 77],
                          [14, 61, 65, 71, 73, 83, 84, 98]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 84))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...KRH'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGGISVGKLLILGCIVALVFGTKKLRTIGEDAGYAIRSFQKALRGDEVTTQSST...QRH'}, length=78)",
        )
        self.assertEqual(hsp.target.id, "gi|15601293|ref|NP_232924.1|")
        self.assertEqual(hsp.target.name, "NP_232924")
        self.assertEqual(
            hsp.target.description,
            "tatA protein [Vibrio cholerae O1 biovar eltor str. N16961] >gi|9657940|gb|AAF96436.1| tatA protein [Vibrio cholerae O1 biovar eltor str. N16961] >gi|11356195|pir||A82448 tatA protein VCA0533 [imported] - Vibrio cholerae (strain N16961 serogroup O1)",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGGIS+ +LLI+  IV L+FGTKKL +IG D G +I+ F+KA+  DE     T+Q +  TA+   D  A   QE + + D++RH",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|156012         0 MGGISVGKLLILGCIVALVFGTKKLRTIGEDAGYAIRSFQKALRGDE----VTTQSS--T
                  0 |||||..............||||||..||.|.|..|..|.||...||----.|.|..--|
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|156012        54 AEESVDSFA-IEQESSHSSDSQRH 77
                 60 |....|..|-..||.....|..|| 84
gi|491764        74 AKTIADKQADTNQEQAKTEDAKRH 98
""",
        )
        hit = hits[68]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|66770480|ref|YP_245242.1|")
        self.assertEqual(hit.target.name, "YP_245242")
        self.assertEqual(
            hit.target.description,
            "sec-independent protein translocase [Xanthomonas campestris pv. campestris str. 8004] >gi|21233515|ref|NP_639432.1| sec-independent protein translocase [Xanthomonas campestris pv. campestris str. ATCC 33913] >gi|21115368|gb|AAM43314.1| sec-independent protein translocase [Xanthomonas campestris pv. campestris str. ATCC 33913] >gi|66575812|gb|AAY51222.1| sec-independent protein translocase [Xanthomonas campestris pv. campestris str. 8004] >gi|24212476|sp|Q8P3H8|TATA_XANCP Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=75)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 113.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 48.1358)
        self.assertAlmostEqual(hsp.annotations["evalue"], 8.02423e-05)
        self.assertEqual(hsp.annotations["identity"], 32)
        self.assertEqual(hsp.annotations["positive"], 48)
        self.assertEqual(hsp.annotations["gaps"], 1)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 58, 58, 73],
                          [14, 72, 73, 88]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 74))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...NQE'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGSFSIWHWLIVLVIVLLVFGTKRLTSGAKDLGSAVKEFKKGMHDDDKPAGKLG...DRD'}, length=75)",
        )
        self.assertEqual(hsp.target.id, "gi|66770480|ref|YP_245242.1|")
        self.assertEqual(hsp.target.name, "YP_245242")
        self.assertEqual(
            hsp.target.description,
            "sec-independent protein translocase [Xanthomonas campestris pv. campestris str. 8004] >gi|21233515|ref|NP_639432.1| sec-independent protein translocase [Xanthomonas campestris pv. campestris str. ATCC 33913] >gi|21115368|gb|AAM43314.1| sec-independent protein translocase [Xanthomonas campestris pv. campestris str. ATCC 33913] >gi|66575812|gb|AAY51222.1| sec-independent protein translocase [Xanthomonas campestris pv. campestris str. 8004] >gi|24212476|sp|Q8P3H8|TATA_XANCP Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG  SIW  LI+ VIV+L+FGTK+L S   DLG+++K FKK M DD+    K   D+  +A+   + QA+ +++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|667704         0 MGSFSIWHWLIVLVIVLLVFGTKRLTSGAKDLGSAVKEFKKGMHDDDKPAGKLGDDSR-S
                  0 ||..|||............||||.|.|...|||...|.|||.|.||.....|...|..-.
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|667704        59 AEQAREAQAERDRD 73
                 60 |......||..... 74
gi|491764        74 AKTIADKQADTNQE 88
""",
        )
        hit = hits[69]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|53804435|ref|YP_113945.1|")
        self.assertEqual(hit.target.name, "YP_113945")
        self.assertEqual(
            hit.target.description,
            "Sec-independent protein translocase protein TatA [Methylococcus capsulatus str. Bath] >gi|53758196|gb|AAU92487.1| Sec-independent protein translocase protein TatA [Methylococcus capsulatus str. Bath]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=70)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 113.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 48.1358)
        self.assertAlmostEqual(hsp.annotations["evalue"], 8.02423e-05)
        self.assertEqual(hsp.annotations["identity"], 30)
        self.assertEqual(hsp.annotations["positive"], 51)
        self.assertEqual(hsp.annotations["gaps"], 6)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 1, 46, 46, 53, 53, 67],
                          [16, 61, 64, 71, 74, 88]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 72))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({16: 'GISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQD...NQE'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({1: 'GIGVWELLLLFLIVLVVFGTKRLRNIGGDLGGAIKSFRQAMSENEDKPSEGGAR...KEK'}, length=70)",
        )
        self.assertEqual(hsp.target.id, "gi|53804435|ref|YP_113945.1|")
        self.assertEqual(hsp.target.name, "YP_113945")
        self.assertEqual(
            hsp.target.description,
            "Sec-independent protein translocase protein TatA [Methylococcus capsulatus str. Bath] >gi|53758196|gb|AAU92487.1| Sec-independent protein translocase protein TatA [Methylococcus capsulatus str. Bath]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "GI +W+LL++ +IV+++FGTK+L +IG DLG +IK F++AMS++E   DK S+     A+T+  +  D  ++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|538044         1 GIGVWELLLLFLIVLVVFGTKRLRNIGGDLGGAIKSFRQAMSENE---DKPSEGG---AR
                  0 ||..|............||||.|..||.|||..||.|..|||..|---||.|...---|.
gi|491764        16 GISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAK

gi|538044        55 TLEGEVVDKKEK 67
                 60 |......|.... 72
gi|491764        76 TIADKQADTNQE 88
""",
        )
        hit = hits[70]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|75825357|ref|ZP_00754793.1|")
        self.assertEqual(hit.target.name, "ZP_00754793")
        self.assertEqual(
            hit.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Vibrio cholerae O395]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=80)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 113.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 48.1358)
        self.assertAlmostEqual(hsp.annotations["evalue"], 8.02423e-05)
        self.assertEqual(hsp.annotations["identity"], 33)
        self.assertEqual(hsp.annotations["positive"], 51)
        self.assertEqual(hsp.annotations["gaps"], 5)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 47, 47, 79],
                          [14, 61, 66, 98]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 84))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...KRH'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGGISVGKLLILGCIVALIFGTKKLRTIGEDAGYAIRSFQKALRSDEVTTQSST...QRH'}, length=80)",
        )
        self.assertEqual(hsp.target.id, "gi|75825357|ref|ZP_00754793.1|")
        self.assertEqual(hsp.target.name, "ZP_00754793")
        self.assertEqual(
            hsp.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Vibrio cholerae O395]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGGIS+ +LLI+  IV L+FGTKKL +IG D G +I+ F+KA+  DE      +  +  T +++        QE + + D++RH",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|758253         0 MGGISVGKLLILGCIVALIFGTKKLRTIGEDAGYAIRSFQKALRSDE-----VTTQSSTT
                  0 |||||..............||||||..||.|.|..|..|.||...||-----.......|
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|758253        55 EESVDSFSFAIEQEPSHSSDSQRH 79
                 60 ............||.....|..|| 84
gi|491764        74 AKTIADKQADTNQEQAKTEDAKRH 98
""",
        )
        hit = hits[71]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|71908987|ref|YP_286574.1|")
        self.assertEqual(hit.target.name, "YP_286574")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Dechloromonas aromatica RCB] >gi|71848608|gb|AAZ48104.1| Twin-arginine translocation protein TatA/E [Dechloromonas aromatica RCB] >gi|41723355|ref|ZP_00150282.1| COG1826: Sec-independent protein secretion pathway components [Dechloromonas aromatica RCB]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=76)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 113.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 48.1358)
        self.assertAlmostEqual(hsp.annotations["evalue"], 8.02423e-05)
        self.assertEqual(hsp.annotations["identity"], 31)
        self.assertEqual(hsp.annotations["positive"], 41)
        self.assertEqual(hsp.annotations["gaps"], 2)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 44, 44, 58],
                          [14, 58, 60, 74]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 60))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGSFSIWHWLIVLVIVMLIFGTKKLRNVGQDLGGAVKGFKDGMKEANADKPAEEAQPT'}, length=76)",
        )
        self.assertEqual(hsp.target.id, "gi|71908987|ref|YP_286574.1|")
        self.assertEqual(hsp.target.name, "YP_286574")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Dechloromonas aromatica RCB] >gi|71848608|gb|AAZ48104.1| Twin-arginine translocation protein TatA/E [Dechloromonas aromatica RCB] >gi|41723355|ref|ZP_00150282.1| COG1826: Sec-independent protein secretion pathway components [Dechloromonas aromatica RCB]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG  SIW  LI+ VIV+L+FGTKKL ++G DLG ++KGFK  M   E   DK +++A  T",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|719089         0 MGSFSIWHWLIVLVIVMLIFGTKKLRNVGQDLGGAVKGFKDGMK--EANADKPAEEAQPT
                  0 ||..|||............||||||...|.|||...||||..|.--|...||....|..|
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|719089        58 
                 60 
gi|491764        74 
""",
        )
        hit = hits[72]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|68526571|gb|EAN49542.1|")
        self.assertEqual(hit.target.name, "EAN49542")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Ralstonia metallidurans CH34] >gi|68559139|ref|ZP_00598474.1| Twin-arginine translocation protein TatA/E [Ralstonia metallidurans CH34]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=77)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 112.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 47.7506)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.0001048)
        self.assertEqual(hsp.annotations["identity"], 36)
        self.assertEqual(hsp.annotations["positive"], 47)
        self.assertEqual(hsp.annotations["gaps"], 9)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 45, 48, 51, 57, 77],
                          [14, 59, 59, 62, 62, 82]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 77))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...DKQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MGSFSIWHWLIVLVIVMLVFGTKKLRNIGQDLGGAVKGFKDGMKEGNTDEPATP...RQQ')",
        )
        self.assertEqual(hsp.target.id, "gi|68526571|gb|EAN49542.1|")
        self.assertEqual(hsp.target.name, "EAN49542")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Ralstonia metallidurans CH34] >gi|68559139|ref|ZP_00598474.1| Twin-arginine translocation protein TatA/E [Ralstonia metallidurans CH34]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG  SIW  LI+ VIV+L+FGTKKL +IG DLG ++KGFK  M +   DEP      K+ + S   D  AK  + +Q",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|685265         0 MGSFSIWHWLIVLVIVMLVFGTKKLRNIGQDLGGAVKGFKDGMKEGNTDEPATPTPAKEL
                  0 ||..|||............||||||..||.|||...||||..|..---|||------|..
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSD---DEP------KQD

gi|685265        60 RDSTTIDVEAKEKSRQQ 77
                 60 ..|...|..||.....| 77
gi|491764        65 KTSQDADFTAKTIADKQ 82
""",
        )
        hit = hits[73]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|71736448|ref|YP_272670.1|")
        self.assertEqual(hit.target.name, "YP_272670")
        self.assertEqual(
            hit.target.description,
            "sec-independent protein translocase TatA [Pseudomonas syringae pv. phaseolicola 1448A] >gi|71557001|gb|AAZ36212.1| sec-independent protein translocase TatA [Pseudomonas syringae pv. phaseolicola 1448A]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=91)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 112.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 47.7506)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.0001048)
        self.assertEqual(hsp.annotations["identity"], 30)
        self.assertEqual(hsp.annotations["positive"], 41)
        self.assertEqual(hsp.annotations["gaps"], 6)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 19,  47,  48,  62,  62,  73,  73,  84],
                          [ 33,  61,  61,  75,  78,  89,  91, 102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 70))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQ...KEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({19: 'FGTKKLKGLGSDVGESIKGFRKAMNDDDKPAEQPAPQPQQAQAAPQGSPLNQPH...VDE'}, length=91)",
        )
        self.assertEqual(hsp.target.id, "gi|71736448|ref|YP_272670.1|")
        self.assertEqual(hsp.target.name, "YP_272670")
        self.assertEqual(
            hsp.target.description,
            "sec-independent protein translocase TatA [Pseudomonas syringae pv. phaseolicola 1448A] >gi|71557001|gb|AAZ36212.1| sec-independent protein translocase TatA [Pseudomonas syringae pv. phaseolicola 1448A]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "FGTKKL  +GSD+G SIKGF+KAM+DD+ P +    Q     A   A + +  NQ    T DA+ H  ++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|717364        19 FGTKKLKGLGSDVGESIKGFRKAMNDDDKPAEQPAPQPQQAQA---APQGSPLNQPH--T
                  0 ||||||...|||.|.|||||.|||.||.-|......|.....|---|......||..--|
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDE-PKQDKTSQDADFTAKTIADKQADTNQEQAKT

gi|717364        74 IDAQAHKVDE  84
                 60 .||..|....  70
gi|491764        92 EDAKRHDKEQ 102
""",
        )
        hit = hits[74]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|56460344|ref|YP_155625.1|")
        self.assertEqual(hit.target.name, "YP_155625")
        self.assertEqual(
            hit.target.description,
            "Sec-independent protein secretion pathway component, TatA family [Idiomarina loihiensis L2TR] >gi|56179354|gb|AAV82076.1| Sec-independent protein secretion pathway component, TatA family [Idiomarina loihiensis L2TR]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=72)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 112.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 47.7506)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.0001048)
        self.assertEqual(hsp.annotations["identity"], 32)
        self.assertEqual(hsp.annotations["positive"], 54)
        self.assertEqual(hsp.annotations["gaps"], 17)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  3,  49,  49,  72],
                          [ 15,  61,  78, 101]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 86))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({15: 'GGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQ...DKE'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({3: 'GGMSMVQLGVLLVIVILLFGSKKLRTLGSDLGSAIKGFKSSITDEPSDREKDIQ...HEQ'}, length=72)",
        )
        self.assertEqual(hsp.target.id, "gi|56460344|ref|YP_155625.1|")
        self.assertEqual(hsp.target.name, "YP_155625")
        self.assertEqual(
            hsp.target.description,
            "Sec-independent protein secretion pathway component, TatA family [Idiomarina loihiensis L2TR] >gi|56179354|gb|AAV82076.1| Sec-independent protein secretion pathway component, TatA family [Idiomarina loihiensis L2TR]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "GG+S+ QL ++ VIV+LLFG+KKL ++GSDLG++IKGFK +++D+                  +D++ D  + Q  T + +RH ++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|564603         3 GGMSMVQLGVLLVIVILLFGSKKLRTLGSDLGSAIKGFKSSITDEP--------------
                  0 ||.|..|...........||.|||...|||||..|||||....|..--------------
gi|491764        15 GGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTA

gi|564603        49 ---SDREKDIQRHQNLTSEEERHHEQ  72
                 60 ---.|...|....|..|....||...  86
gi|491764        75 KTIADKQADTNQEQAKTEDAKRHDKE 101
""",
        )
        hit = hits[75]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|68214708|ref|ZP_00566522.1|")
        self.assertEqual(hit.target.name, "ZP_00566522")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Methylobacillus flagellatus KT] >gi|68186844|gb|EAN01544.1| Twin-arginine translocation protein TatA/E [Methylobacillus flagellatus KT]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=72)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 111.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 47.3654)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.000136873)
        self.assertEqual(hsp.annotations["identity"], 29)
        self.assertEqual(hsp.annotations["positive"], 44)
        self.assertEqual(hsp.annotations["gaps"], 2)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 43, 43, 62],
                          [14, 57, 59, 78]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 64))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...KTI'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGSFSIWHWLIVLAIVILIFGTKRLRNLGSDVGGAVKGFKEAVNEGKSAAAALD...QTV'}, length=72)",
        )
        self.assertEqual(hsp.target.id, "gi|68214708|ref|ZP_00566522.1|")
        self.assertEqual(hsp.target.name, "ZP_00566522")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Methylobacillus flagellatus KT] >gi|68186844|gb|EAN01544.1| Twin-arginine translocation protein TatA/E [Methylobacillus flagellatus KT]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG  SIW  LI+  IV+L+FGTK+L ++GSD+G ++KGFK+A+  +E K    + D D   +T+",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|682147         0 MGSFSIWHWLIVLAIVILIFGTKRLRNLGSDVGGAVKGFKEAV--NEGKSAAAALDDDAK
                  0 ||..|||............||||.|...|||.|...||||.|.--.|.|......|.|..
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|682147        58 GQTV 62
                 60 ..|. 64
gi|491764        74 AKTI 78
""",
        )
        hit = hits[76]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|30248650|ref|NP_840720.1|")
        self.assertEqual(hit.target.name, "NP_840720")
        self.assertEqual(
            hit.target.description,
            "mttA/Hcf106 family [Nitrosomonas europaea ATCC 19718] >gi|30180245|emb|CAD84550.1| mttA/Hcf106 family [Nitrosomonas europaea ATCC 19718]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=76)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 111.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 47.3654)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.000136873)
        self.assertEqual(hsp.annotations["identity"], 29)
        self.assertEqual(hsp.annotations["positive"], 43)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 73],
                          [14, 87]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 73))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...TNQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGSFSIWHWLVVLAIVVLVFGTKKLRNLGSDLGGAVRGFKEGMKGAEEESTPPP...KDQ'}, length=76)",
        )
        self.assertEqual(hsp.target.id, "gi|30248650|ref|NP_840720.1|")
        self.assertEqual(hsp.target.name, "NP_840720")
        self.assertEqual(
            hsp.target.description,
            "mttA/Hcf106 family [Nitrosomonas europaea ATCC 19718] >gi|30180245|emb|CAD84550.1| mttA/Hcf106 family [Nitrosomonas europaea ATCC 19718]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG  SIW  L++  IVVL+FGTKKL ++GSDLG +++GFK+ M   E +          T  +I  +  + +Q",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|302486         0 MGSFSIWHWLVVLAIVVLVFGTKKLRNLGSDLGGAVRGFKEGMKGAEEESTPPPPAQQVT
                  0 ||..|||............||||||...|||||....|||..|...|............|
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|302486        60 GHSIKSEIEEKDQ 73
                 60 ...|........| 73
gi|491764        74 AKTIADKQADTNQ 87
""",
        )
        hit = hits[77]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|75822907|ref|ZP_00752458.1|")
        self.assertEqual(hit.target.name, "ZP_00752458")
        self.assertEqual(
            hit.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Vibrio cholerae RC385]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=78)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 111.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 47.3654)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.000136873)
        self.assertEqual(hsp.annotations["identity"], 33)
        self.assertEqual(hsp.annotations["positive"], 50)
        self.assertEqual(hsp.annotations["gaps"], 7)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 64, 64, 77],
                          [14, 78, 85, 98]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 84))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...KRH'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGGISVGKLLILGCIVALVFGTKKLRTIGEDAGYAIRSFQKALRGDEVTTQSST...QRH'}, length=78)",
        )
        self.assertEqual(hsp.target.id, "gi|75822907|ref|ZP_00752458.1|")
        self.assertEqual(hsp.target.name, "ZP_00752458")
        self.assertEqual(
            hsp.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Vibrio cholerae RC385]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGGIS+ +LLI+  IV L+FGTKKL +IG D G +I+ F+KA+  DE     ++ +    +  I        QE + + D++RH",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|758229         0 MGGISVGKLLILGCIVALVFGTKKLRTIGEDAGYAIRSFQKALRGDEVTTQSSTTEESVD
                  0 |||||..............||||||..||.|.|..|..|.||...||.............
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|758229        60 SFAI-------EQESSHSSDSQRH 77
                 60 ...|-------.||.....|..|| 84
gi|491764        74 AKTIADKQADTNQEQAKTEDAKRH 98
""",
        )
        hit = hits[78]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|70733926|ref|YP_257566.1|")
        self.assertEqual(hit.target.name, "YP_257566")
        self.assertEqual(
            hit.target.description,
            "sec-independent protein translocase TatA [Pseudomonas fluorescens Pf-5] >gi|68348225|gb|AAY95831.1| sec-independent protein translocase TatA [Pseudomonas fluorescens Pf-5]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=93)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 111.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 47.3654)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.000136873)
        self.assertEqual(hsp.annotations["identity"], 29)
        self.assertEqual(hsp.annotations["positive"], 40)
        self.assertEqual(hsp.annotations["gaps"], 7)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[19, 68, 69, 74, 80, 92],
                          [33, 82, 82, 87, 87, 99]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 73))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQ...RHD'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({19: 'FGTKKLKNLGTDVGESIKGFRKAMNDDEKPAEPVVPPAAQPVPPVQPQQSAPLN...RKD'}, length=93)",
        )
        self.assertEqual(hsp.target.id, "gi|70733926|ref|YP_257566.1|")
        self.assertEqual(hsp.target.name, "YP_257566")
        self.assertEqual(
            hsp.target.description,
            "sec-independent protein translocase TatA [Pseudomonas fluorescens Pf-5] >gi|68348225|gb|AAY95831.1| sec-independent protein translocase TatA [Pseudomonas fluorescens Pf-5]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "FGTKKL ++G+D+G SIKGF+KAM+DDE   +     A      +  +Q A  NQ      +  K E+  R D",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|707339        19 FGTKKLKNLGTDVGESIKGFRKAMNDDEKPAEPVVPPAAQPVPPVQPQQSAPLNQPHTID
                  0 ||||||...|.|.|.|||||.|||.|||.........|..........|-|..||-----
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQ-ADTNQ-----

gi|707339        79 VQAQKVEEPTRKD 92
                 60 -...|.|...|.| 73
gi|491764        87 -EQAKTEDAKRHD 99
""",
        )
        hit = hits[79]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|63254358|gb|AAY35454.1|")
        self.assertEqual(hit.target.name, "AAY35454")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Pseudomonas syringae pv. syringae B728a] >gi|66043651|ref|YP_233492.1| Twin-arginine translocation protein TatA/E [Pseudomonas syringae pv. syringae B728a]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=91)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 110.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 46.9802)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.000178761)
        self.assertEqual(hsp.annotations["identity"], 27)
        self.assertEqual(hsp.annotations["positive"], 38)
        self.assertEqual(hsp.annotations["gaps"], 4)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 19,  46,  46,  84],
                          [ 33,  60,  64, 102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 69))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQ...KEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({19: 'FGTKKLKGLGSDVGESIKGFRKAMNDDDKPAEQPAPQPQQAQPAPQGSPLNQPH...VDE'}, length=91)",
        )
        self.assertEqual(hsp.target.id, "gi|63254358|gb|AAY35454.1|")
        self.assertEqual(hsp.target.name, "AAY35454")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Pseudomonas syringae pv. syringae B728a] >gi|66043651|ref|YP_233492.1| Twin-arginine translocation protein TatA/E [Pseudomonas syringae pv. syringae B728a]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "FGTKKL  +GSD+G SIKGF+KAM+DD    DK ++      +        +   Q  T DA+ H  ++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|632543        19 FGTKKLKGLGSDVGESIKGFRKAMNDD----DKPAEQPAPQPQQAQPAPQGSPLNQPHTI
                  0 ||||||...|||.|.|||||.|||.||----||......................|..|.
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTE

gi|632543        75 DAQAHKVDE  84
                 60 ||..|....  69
gi|491764        93 DAKRHDKEQ 102
""",
        )
        hit = hits[80]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|73354814|gb|AAZ75668.1|")
        self.assertEqual(hit.target.name, "AAZ75668")
        self.assertEqual(
            hit.target.description, "TatA [Pseudomonas syringae pv. maculicola]"
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=91)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 110.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 46.9802)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.000178761)
        self.assertEqual(hsp.annotations["identity"], 20)
        self.assertEqual(hsp.annotations["positive"], 27)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[19, 52],
                          [33, 66]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 33))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDK'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({19: 'FGTKKLKGLGSDVGESIKGFRKAMHDDDKPEEQ'}, length=91)",
        )
        self.assertEqual(hsp.target.id, "gi|73354814|gb|AAZ75668.1|")
        self.assertEqual(hsp.target.name, "AAZ75668")
        self.assertEqual(
            hsp.target.description, "TatA [Pseudomonas syringae pv. maculicola]"
        )
        self.assertEqual(
            hsp.annotations["midline"], "FGTKKL  +GSD+G SIKGF+KAM DD+  +++"
        )
        self.assertEqual(
            str(hsp),
            """\
gi|733548        19 FGTKKLKGLGSDVGESIKGFRKAMHDDDKPEEQ 52
                  0 ||||||...|||.|.|||||.|||.||...... 33
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDK 66
""",
        )
        hit = hits[81]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|50083761|ref|YP_045271.1|")
        self.assertEqual(hit.target.name, "YP_045271")
        self.assertEqual(
            hit.target.description,
            "Sec-independent protein translocase protein [Acinetobacter sp. ADP1] >gi|49529737|emb|CAG67449.1| Sec-independent protein translocase protein [Acinetobacter sp. ADP1]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=78)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 109.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 46.595)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.000233469)
        self.assertEqual(hsp.annotations["identity"], 30)
        self.assertEqual(hsp.annotations["positive"], 43)
        self.assertEqual(hsp.annotations["gaps"], 2)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 5, 50, 52, 70],
                          [14, 59, 59, 77]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 65))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...AKT'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({5: 'MAGLSIWHVVIFLIVVVLLFGTSKLKNLGKDVGGAVKDFKKSMRDETEENAQLH...VKT'}, length=78)",
        )
        self.assertEqual(hsp.target.id, "gi|50083761|ref|YP_045271.1|")
        self.assertEqual(hsp.target.name, "YP_045271")
        self.assertEqual(
            hsp.target.description,
            "Sec-independent protein translocase protein [Acinetobacter sp. ADP1] >gi|49529737|emb|CAG67449.1| Sec-independent protein translocase protein [Acinetobacter sp. ADP1]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "M G+SIW ++I  ++VVLLFGT KL ++G D+G ++K FKK+M D  +E  Q  T +  D   KT",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|500837         5 MAGLSIWHVVIFLIVVVLLFGTSKLKNLGKDVGGAVKDFKKSMRDETEENAQLHTPRTID
                  0 |.|.|||............|||.||...|.|.|...|.|||.|.|--.|..|..|....|
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSD--DEPKQDKTSQDAD

gi|500837        65 AEVKT 70
                 60 ...|| 65
gi|491764        72 FTAKT 77
""",
        )
        hit = hits[82]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|71548504|ref|ZP_00668728.1|")
        self.assertEqual(hit.target.name, "ZP_00668728")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Nitrosomonas eutropha C71] >gi|71485685|gb|EAO18234.1| Twin-arginine translocation protein TatA/E [Nitrosomonas eutropha C71]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=77)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 108.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 46.2098)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.00030492)
        self.assertEqual(hsp.annotations["identity"], 26)
        self.assertEqual(hsp.annotations["positive"], 35)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 47],
                          [14, 61]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 47))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDE'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGSFSIWHWLVVLAIVVLVFGTKKLRNLGSDLGGAVRGFKEGMKGAE'}, length=77)",
        )
        self.assertEqual(hsp.target.id, "gi|71548504|ref|ZP_00668728.1|")
        self.assertEqual(hsp.target.name, "ZP_00668728")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Nitrosomonas eutropha C71] >gi|71485685|gb|EAO18234.1| Twin-arginine translocation protein TatA/E [Nitrosomonas eutropha C71]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG  SIW  L++  IVVL+FGTKKL ++GSDLG +++GFK+ M   E",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|715485         0 MGSFSIWHWLVVLAIVVLVFGTKKLRNLGSDLGGAVRGFKEGMKGAE 47
                  0 ||..|||............||||||...|||||....|||..|...| 47
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDE 61
""",
        )
        hit = hits[83]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|55247002|gb|EAL42253.1|")
        self.assertEqual(hit.target.name, "EAL42253")
        self.assertEqual(
            hit.target.description,
            "ENSANGP00000028218 [Anopheles gambiae str. PEST] >gi|57965302|ref|XP_561091.1| ENSANGP00000028218 [Anopheles gambiae str. PEST]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=53)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 107.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 45.8246)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.000398238)
        self.assertEqual(hsp.annotations["identity"], 26)
        self.assertEqual(hsp.annotations["positive"], 37)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 2, 49],
                          [14, 61]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 47))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDE'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({2: 'LGGIHIWHLLILLVVVVVVFGSKRLAGAGEDLGTAIRDFRKALRDDE'}, length=53)",
        )
        self.assertEqual(hsp.target.id, "gi|55247002|gb|EAL42253.1|")
        self.assertEqual(hsp.target.name, "EAL42253")
        self.assertEqual(
            hsp.target.description,
            "ENSANGP00000028218 [Anopheles gambiae str. PEST] >gi|57965302|ref|XP_561091.1| ENSANGP00000028218 [Anopheles gambiae str. PEST]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "+GGI IW LLI+ V+VV++FG+K+L   G DLG +I+ F+KA+ DDE",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|552470         2 LGGIHIWHLLILLVVVVVVFGSKRLAGAGEDLGTAIRDFRKALRDDE 49
                  0 .|||.||............||.|.|...|.|||..|..|.||..||| 47
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDE 61
""",
        )
        hit = hits[84]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|50084688|ref|YP_046198.1|")
        self.assertEqual(hit.target.name, "YP_046198")
        self.assertEqual(
            hit.target.description,
            "Sec-independent protein secretion pathway, translocase protein [Acinetobacter sp. ADP1] >gi|49530664|emb|CAG68376.1| Sec-independent protein secretion pathway, translocase protein [Acinetobacter sp. ADP1]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=71)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 106.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 45.4394)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.000520115)
        self.assertEqual(hsp.annotations["identity"], 29)
        self.assertEqual(hsp.annotations["positive"], 41)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 1, 59],
                          [16, 74]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 58))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({16: 'GISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({1: 'GLSIWHILLFAVIVILLFGTSKLKSLGKDLGGAIKDFKDSVNSDDQKLLQNQQYKDIS'}, length=71)",
        )
        self.assertEqual(hsp.target.id, "gi|50084688|ref|YP_046198.1|")
        self.assertEqual(hsp.target.name, "YP_046198")
        self.assertEqual(
            hsp.target.description,
            "Sec-independent protein secretion pathway, translocase protein [Acinetobacter sp. ADP1] >gi|49530664|emb|CAG68376.1| Sec-independent protein secretion pathway, translocase protein [Acinetobacter sp. ADP1]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "G+SIW +L+ AVIV+LLFGT KL S+G DLG +IK FK +++ D+ K  +  Q  D +",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|500846         1 GLSIWHILLFAVIVILLFGTSKLKSLGKDLGGAIKDFKDSVNSDDQKLLQNQQYKDIS 59
                  0 |.|||............|||.||.|.|.|||..||.||.....|..|.....|..|.. 58
gi|491764        16 GISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT 74
""",
        )
        hit = hits[85]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|28872267|ref|NP_794886.1|")
        self.assertEqual(hit.target.name, "NP_794886")
        self.assertEqual(
            hit.target.description,
            "sec-independent protein translocase TatA [Pseudomonas syringae pv. tomato str. DC3000] >gi|28855521|gb|AAO58581.1| sec-independent protein translocase TatA [Pseudomonas syringae pv. tomato str. DC3000]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=91)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 105.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 45.0542)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.000679292)
        self.assertEqual(hsp.annotations["identity"], 26)
        self.assertEqual(hsp.annotations["positive"], 37)
        self.assertEqual(hsp.annotations["gaps"], 4)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 19,  46,  46,  84],
                          [ 33,  60,  64, 102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 69))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQ...KEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({19: 'FGTKKLKGLGTDVGESIKGFRKAMHDDDKPAEQPAPQPQQAQPAPQGSPLNQPH...VDE'}, length=91)",
        )
        self.assertEqual(hsp.target.id, "gi|28872267|ref|NP_794886.1|")
        self.assertEqual(hsp.target.name, "NP_794886")
        self.assertEqual(
            hsp.target.description,
            "sec-independent protein translocase TatA [Pseudomonas syringae pv. tomato str. DC3000] >gi|28855521|gb|AAO58581.1| sec-independent protein translocase TatA [Pseudomonas syringae pv. tomato str. DC3000]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "FGTKKL  +G+D+G SIKGF+KAM DD    DK ++      +        +   Q  T DA+ H  ++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|288722        19 FGTKKLKGLGTDVGESIKGFRKAMHDD----DKPAEQPAPQPQQAQPAPQGSPLNQPHTI
                  0 ||||||...|.|.|.|||||.|||.||----||......................|..|.
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTE

gi|288722        75 DAQAHKVDE  84
                 60 ||..|....  69
gi|491764        93 DAKRHDKEQ 102
""",
        )
        hit = hits[86]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|49082486|gb|AAT50643.1|")
        self.assertEqual(hit.target.name, "AAT50643")
        self.assertEqual(hit.target.description, "PA5068 [synthetic construct]")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=83)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 105.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 45.0542)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.000679292)
        self.assertEqual(hsp.annotations["identity"], 25)
        self.assertEqual(hsp.annotations["positive"], 38)
        self.assertEqual(hsp.annotations["gaps"], 2)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[19, 47, 47, 81],
                          [33, 61, 63, 97]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 64))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQ...AKR'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({19: 'FGTKRLKNLGSDVGEAIKGFRKAVNTEEDDKKDQPAAQPAQPLNQPHTIDAQAQ...ARK'}, length=83)",
        )
        self.assertEqual(hsp.target.id, "gi|49082486|gb|AAT50643.1|")
        self.assertEqual(hsp.target.name, "AAT50643")
        self.assertEqual(hsp.target.description, "PA5068 [synthetic construct]")
        self.assertEqual(
            hsp.annotations["midline"],
            "FGTK+L ++GSD+G +IKGF+KA++ +E   DK  Q A   A+ +        Q Q   E A++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|490824        19 FGTKRLKNLGSDVGEAIKGFRKAVNTEE--DDKKDQPAAQPAQPLNQPHTIDAQAQKVEE
                  0 ||||.|...|||.|..||||.||....|--.||..|.|...|...........|.|...|
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTE

gi|490824        77 PARK 81
                 60 .|.. 64
gi|491764        93 DAKR 97
""",
        )
        hit = hits[87]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|53726598|ref|ZP_00141543.2|")
        self.assertEqual(hit.target.name, "ZP_00141543")
        self.assertEqual(
            hit.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Pseudomonas aeruginosa UCBPP-PA14] >gi|9951361|gb|AAG08453.1| translocation protein TatA [Pseudomonas aeruginosa PAO1] >gi|15600261|ref|NP_253755.1| translocation protein TatA [Pseudomonas aeruginosa PAO1] >gi|11352708|pir||G83011 translocation protein TatA PA5068 [imported] - Pseudomonas aeruginosa (strain PAO1) >gi|24212512|sp|Q9HUB5|TATA_PSEAE Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=82)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 105.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 45.0542)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.000679292)
        self.assertEqual(hsp.annotations["identity"], 25)
        self.assertEqual(hsp.annotations["positive"], 38)
        self.assertEqual(hsp.annotations["gaps"], 2)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[19, 47, 47, 81],
                          [33, 61, 63, 97]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 64))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQ...AKR'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({19: 'FGTKRLKNLGSDVGEAIKGFRKAVNTEEDDKKDQPAAQPAQPLNQPHTIDAQAQ...ARK'}, length=82)",
        )
        self.assertEqual(hsp.target.id, "gi|53726598|ref|ZP_00141543.2|")
        self.assertEqual(hsp.target.name, "ZP_00141543")
        self.assertEqual(
            hsp.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Pseudomonas aeruginosa UCBPP-PA14] >gi|9951361|gb|AAG08453.1| translocation protein TatA [Pseudomonas aeruginosa PAO1] >gi|15600261|ref|NP_253755.1| translocation protein TatA [Pseudomonas aeruginosa PAO1] >gi|11352708|pir||G83011 translocation protein TatA PA5068 [imported] - Pseudomonas aeruginosa (strain PAO1) >gi|24212512|sp|Q9HUB5|TATA_PSEAE Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "FGTK+L ++GSD+G +IKGF+KA++ +E   DK  Q A   A+ +        Q Q   E A++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|537265        19 FGTKRLKNLGSDVGEAIKGFRKAVNTEE--DDKKDQPAAQPAQPLNQPHTIDAQAQKVEE
                  0 ||||.|...|||.|..||||.||....|--.||..|.|...|...........|.|...|
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTE

gi|537265        77 PARK 81
                 60 .|.. 64
gi|491764        93 DAKR 97
""",
        )
        hit = hits[88]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|68213616|ref|ZP_00565447.1|")
        self.assertEqual(hit.target.name, "ZP_00565447")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Methylobacillus flagellatus KT] >gi|68187989|gb|EAN02672.1| Twin-arginine translocation protein TatA/E [Methylobacillus flagellatus KT]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=54)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 104.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 44.669)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.000887182)
        self.assertEqual(hsp.annotations["identity"], 23)
        self.assertEqual(hsp.annotations["positive"], 34)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 47],
                          [14, 61]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 47))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDE'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGSFSIWHWLIVLAVVTMVFGTKRLRNMGSDLGGALKNFKEASKDPD'}, length=54)",
        )
        self.assertEqual(hsp.target.id, "gi|68213616|ref|ZP_00565447.1|")
        self.assertEqual(hsp.target.name, "ZP_00565447")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Methylobacillus flagellatus KT] >gi|68187989|gb|EAN02672.1| Twin-arginine translocation protein TatA/E [Methylobacillus flagellatus KT]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG  SIW  LI+  +V ++FGTK+L ++GSDLG ++K FK+A  D +",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|682136         0 MGSFSIWHWLIVLAVVTMVFGTKRLRNMGSDLGGALKNFKEASKDPD 47
                  0 ||..|||............||||.|...|||||...|.||.|..|.. 47
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDE 61
""",
        )
        hit = hits[89]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|74023810|ref|ZP_00694377.1|")
        self.assertEqual(hit.target.name, "ZP_00694377")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Rhodoferax ferrireducens DSM 15236] >gi|72603427|gb|EAO39432.1| Twin-arginine translocation protein TatA/E [Rhodoferax ferrireducens DSM 15236]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=79)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 103.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 44.2838)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.0011587)
        self.assertEqual(hsp.annotations["identity"], 24)
        self.assertEqual(hsp.annotations["positive"], 32)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[18, 78],
                          [33, 93]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 60))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTE'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({18: 'FGTKKLRNIGSDLGGAVKGFKDGMKEGSDKAADAPAAAPQQVASSATAAKETIDVEAKTK'}, length=79)",
        )
        self.assertEqual(hsp.target.id, "gi|74023810|ref|ZP_00694377.1|")
        self.assertEqual(hsp.target.name, "ZP_00694377")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Rhodoferax ferrireducens DSM 15236] >gi|72603427|gb|EAO39432.1| Twin-arginine translocation protein TatA/E [Rhodoferax ferrireducens DSM 15236]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "FGTKKL +IGSDLG ++KGFK  M +   K       A     + A    +T   +AKT+",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|740238        18 FGTKKLRNIGSDLGGAVKGFKDGMKEGSDKAADAPAAAPQQVASSATAAKETIDVEAKTK
                  0 ||||||..||||||...||||..|.....|.......|.......|.....|....|||.
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTE

gi|740238        78 
                 60 
gi|491764        93 
""",
        )
        hit = hits[90]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|71066554|ref|YP_265281.1|")
        self.assertEqual(hit.target.name, "YP_265281")
        self.assertEqual(
            hit.target.description,
            "twin-arginine translocation protein TatA/E [Psychrobacter arcticus 273-4] >gi|71039539|gb|AAZ19847.1| twin-arginine translocation protein TatA/E [Psychrobacter arcticus 273-4]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=87)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 102.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 43.8986)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.0015133)
        self.assertEqual(hsp.annotations["identity"], 31)
        self.assertEqual(hsp.annotations["positive"], 43)
        self.assertEqual(hsp.annotations["gaps"], 4)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 52, 56, 75],
                          [14, 66, 66, 85]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 75))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...ADT'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGSFSITHWLILLVVVVVVFGTSKLRNAGKDLGGAVKGFKEAVKDENTEHAKKH...SDT'}, length=87)",
        )
        self.assertEqual(hsp.target.id, "gi|71066554|ref|YP_265281.1|")
        self.assertEqual(hsp.target.name, "YP_265281")
        self.assertEqual(
            hsp.target.description,
            "twin-arginine translocation protein TatA/E [Psychrobacter arcticus 273-4] >gi|71039539|gb|AAZ19847.1| twin-arginine translocation protein TatA/E [Psychrobacter arcticus 273-4]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG  SI   LI+ V+VV++FGT KL + G DLG ++KGFK+A+ D+  +  K       DA      I   Q+DT",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|710665         0 MGSFSITHWLILLVVVVVVFGTSKLRNAGKDLGGAVKGFKEAVKDENTEHAKKHVVLDHD
                  0 ||..||.............|||.||...|.|||...||||.|..|......|----...|
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDK----TSQD

gi|710665        60 AGTNPPNITGTQSDT 75
                 60 |......|...|.|| 75
gi|491764        70 ADFTAKTIADKQADT 85
""",
        )
        hit = hits[91]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|15611372|ref|NP_223023.1|")
        self.assertEqual(hit.target.name, "NP_223023")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein jhp0303 [Helicobacter pylori J99] >gi|4154834|gb|AAD05888.1| putative [Helicobacter pylori J99] >gi|7444815|pir||A71948 hypothetical protein jhp0303 - Helicobacter pylori (strain J99) >gi|9979057|sp|Q9ZMB8|TATA_HELPJ Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=79)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 101.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 43.5134)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.00197644)
        self.assertEqual(hsp.annotations["identity"], 31)
        self.assertEqual(hsp.annotations["positive"], 48)
        self.assertEqual(hsp.annotations["gaps"], 3)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0,  4,  5, 59, 59, 78],
                          [14, 18, 18, 72, 74, 93]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 80))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...KTE'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGGFTSIWHWVIVLLVIVLLFGAKKIPELAKGLGSGIKNFKKAVKDDEEEAKNE...KQE'}, length=79)",
        )
        self.assertEqual(hsp.target.id, "gi|15611372|ref|NP_223023.1|")
        self.assertEqual(hsp.target.name, "NP_223023")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein jhp0303 [Helicobacter pylori J99] >gi|4154834|gb|AAD05888.1| putative [Helicobacter pylori J99] >gi|7444815|pir||A71948 hypothetical protein jhp0303 - Helicobacter pylori (strain J99) >gi|9979057|sp|Q9ZMB8|TATA_HELPJ Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGG  SIW  +I+ +++VLLFG KK+  +   LG+ IK FKKA+ DDE +     +  D  A+    K  +T++ ++K E",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|156113         0 MGGFTSIWHWVIVLLVIVLLFGAKKIPELAKGLGSGIKNFKKAVKDDEEEAKNELKTLD-
                  0 |||.-|||............||.||.......||..||.||||..|||..........|-
gi|491764        14 MGGI-SIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADF

gi|156113        59 -AQATQTKVHETSEIKSKQE 78
                 60 -|.....|...|.....|.| 80
gi|491764        73 TAKTIADKQADTNQEQAKTE 93
""",
        )
        hit = hits[92]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|13471183|ref|NP_102752.1|")
        self.assertEqual(hit.target.name, "NP_102752")
        self.assertEqual(
            hit.target.description,
            "sec-independent protein translocase protein tatA/E [Mesorhizobium loti MAFF303099] >gi|14021927|dbj|BAB48538.1| sec-independent protein translocase protein TatA/E [Mesorhizobium loti MAFF303099] >gi|24212508|sp|Q98LC5|TATA_RHILO Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=73)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 99.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 42.743)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.00337129)
        self.assertEqual(hsp.annotations["identity"], 28)
        self.assertEqual(hsp.annotations["positive"], 41)
        self.assertEqual(hsp.annotations["gaps"], 2)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 56, 58, 69],
                          [14, 70, 70, 81]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 69))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...ADK'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGSFSIWHWMIVLVIVLLVFGRGKIPELMGDMAKGIKSFKKGMADDDVADDKRT...KEK'}, length=73)",
        )
        self.assertEqual(hsp.target.id, "gi|13471183|ref|NP_102752.1|")
        self.assertEqual(hsp.target.name, "NP_102752")
        self.assertEqual(
            hsp.target.description,
            "sec-independent protein translocase protein tatA/E [Mesorhizobium loti MAFF303099] >gi|14021927|dbj|BAB48538.1| sec-independent protein translocase protein TatA/E [Mesorhizobium loti MAFF303099] >gi|24212508|sp|Q98LC5|TATA_RHILO Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG  SIW  +I+ VIV+L+FG  K+  +  D+   IK FKK M+DD+   DK + +  AD T   + +K",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|134711         0 MGSFSIWHWMIVLVIVLLVFGRGKIPELMGDMAKGIKSFKKGMADDDVADDKRTVEHRAD
                  0 ||..|||............||..|......|....||.|||.|.||....||....--||
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQD--AD

gi|134711        60 ETVSAVKEK 69
                 60 .|......| 69
gi|491764        72 FTAKTIADK 81
""",
        )
        hit = hits[93]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|42523995|ref|NP_969375.1|")
        self.assertEqual(hit.target.name, "NP_969375")
        self.assertEqual(
            hit.target.description,
            "twin-arginine-dependent translocase protein [Bdellovibrio bacteriovorus HD100] >gi|39576203|emb|CAE80368.1| twin-arginine-dependent translocase protein [Bdellovibrio bacteriovorus HD100]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=81)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 99.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 42.743)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.00337129)
        self.assertEqual(hsp.annotations["identity"], 21)
        self.assertEqual(hsp.annotations["positive"], 36)
        self.assertEqual(hsp.annotations["gaps"], 5)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[18, 55, 55, 76],
                          [33, 70, 75, 96]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 63))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQ...DAK'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({18: 'FGPSKLPNLGRSLGESIRGFKKGLNEDPAADEKEAKQQITQNPQRPVNEQQPQTEEKK'}, length=81)",
        )
        self.assertEqual(hsp.target.id, "gi|42523995|ref|NP_969375.1|")
        self.assertEqual(hsp.target.name, "NP_969375")
        self.assertEqual(
            hsp.target.description,
            "twin-arginine-dependent translocase protein [Bdellovibrio bacteriovorus HD100] >gi|39576203|emb|CAE80368.1| twin-arginine-dependent translocase protein [Bdellovibrio bacteriovorus HD100]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "FG  KL ++G  LG SI+GFKK +++D    +K ++      +   + Q   N++Q +TE+ K",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|425239        18 FGPSKLPNLGRSLGESIRGFKKGLNEDPAADEKEAKQ-----QITQNPQRPVNEQQPQTE
                  0 ||..||...|..||.||.||||....|.....|....-----......|...|..|..||
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTE

gi|425239        73 EKK 76
                 60 ..| 63
gi|491764        93 DAK 96
""",
        )
        hit = hits[94]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|67158086|ref|ZP_00419176.1|")
        self.assertEqual(hit.target.name, "ZP_00419176")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Azotobacter vinelandii AvOP] >gi|67085069|gb|EAM04546.1| Twin-arginine translocation protein TatA/E [Azotobacter vinelandii AvOP]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=85)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 98.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 42.3578)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.00440304)
        self.assertEqual(hsp.annotations["identity"], 21)
        self.assertEqual(hsp.annotations["positive"], 35)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[19, 85],
                          [33, 99]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 66))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQ...RHD'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({19: 'FGTKRLKTLGSDIGEAIKGFRKAVNTEEGENRPAEPQTGTSAGDTLNKTQTIEG...RKD'}, length=85)",
        )
        self.assertEqual(hsp.target.id, "gi|67158086|ref|ZP_00419176.1|")
        self.assertEqual(hsp.target.name, "ZP_00419176")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Azotobacter vinelandii AvOP] >gi|67085069|gb|EAM04546.1| Twin-arginine translocation protein TatA/E [Azotobacter vinelandii AvOP]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "FGTK+L ++GSD+G +IKGF+KA++ +E +          +A    +K      +  K +   R D",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|671580        19 FGTKRLKTLGSDIGEAIKGFRKAVNTEEGENRPAEPQTGTSAGDTLNKTQTIEGQAQKVD
                  0 ||||.|...|||.|..||||.||....|.............|.....|.........|..
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTE

gi|671580        79 TPVRKD 85
                 60 ...|.| 66
gi|491764        93 DAKRHD 99
""",
        )
        hit = hits[95]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|15644948|ref|NP_207118.1|")
        self.assertEqual(hit.target.name, "NP_207118")
        self.assertEqual(
            hit.target.description,
            "conserved hypothetical secreted protein [Helicobacter pylori 26695] >gi|2313428|gb|AAD07397.1| conserved hypothetical secreted protein [Helicobacter pylori 26695] >gi|7429236|pir||H64559 conserved hypothetical secreted protein HP0320 - Helicobacter pylori (strain 26695) >gi|9978986|sp|O25088|TATA_HELPY Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=79)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 98.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 42.3578)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.00440304)
        self.assertEqual(hsp.annotations["identity"], 29)
        self.assertEqual(hsp.annotations["positive"], 47)
        self.assertEqual(hsp.annotations["gaps"], 3)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0,  4,  5, 47, 49, 79],
                          [14, 18, 18, 60, 60, 90]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 79))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...EQA'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MGGFTSIWHWVIVLLVIVLLFGAKKIPELAKGLGSGIKNFKKAVKDDEEEAKNE...QES')",
        )
        self.assertEqual(hsp.target.id, "gi|15644948|ref|NP_207118.1|")
        self.assertEqual(hsp.target.name, "NP_207118")
        self.assertEqual(
            hsp.target.description,
            "conserved hypothetical secreted protein [Helicobacter pylori 26695] >gi|2313428|gb|AAD07397.1| conserved hypothetical secreted protein [Helicobacter pylori 26695] >gi|7429236|pir||H64559 conserved hypothetical secreted protein HP0320 - Helicobacter pylori (strain 26695) >gi|9978986|sp|O25088|TATA_HELPY Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGG  SIW  +I+ +++VLLFG KK+  +   LG+ IK FKKA+ DD  E K +  + DA  T   + +     +++++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|156449         0 MGGFTSIWHWVIVLLVIVLLFGAKKIPELAKGLGSGIKNFKKAVKDDEEEAKNEPKTLDA
                  0 |||.-|||............||.||.......||..||.||||..||--|.|......||
gi|491764        14 MGGI-SIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDD--EPKQDKTSQDA

gi|156449        60 QATQTKVHESSEIKSKQES 79
                 60 ..|................ 79
gi|491764        71 DFTAKTIADKQADTNQEQA 90
""",
        )
        hit = hits[96]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|13277311|emb|CAC34414.1|")
        self.assertEqual(hit.target.name, "CAC34414")
        self.assertEqual(
            hit.target.description, "putative TatA protein [Legionella pneumophila]"
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=61)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 96.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 41.5874)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.00751045)
        self.assertEqual(hsp.annotations["identity"], 26)
        self.assertEqual(hsp.annotations["positive"], 35)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 2, 49],
                          [14, 61]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 47))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDE'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({2: 'LSGISPLSLLLILAIIIALFGTSKLKTIGSDLGEAIKNFRKAMNSEE'}, length=61)",
        )
        self.assertEqual(hsp.target.id, "gi|13277311|emb|CAC34414.1|")
        self.assertEqual(hsp.target.name, "CAC34414")
        self.assertEqual(
            hsp.target.description, "putative TatA protein [Legionella pneumophila]"
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "+ GIS   LL+I  I++ LFGT KL +IGSDLG +IK F+KAM+ +E",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|132773         2 LSGISPLSLLLILAIIIALFGTSKLKTIGSDLGEAIKNFRKAMNSEE 49
                  0 ..|||..............|||.||..||||||..||.|.|||...| 47
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDE 61
""",
        )
        hit = hits[97]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|54298906|ref|YP_125275.1|")
        self.assertEqual(hit.target.name, "YP_125275")
        self.assertEqual(
            hit.target.description,
            "Putative TatA protein(twin arginine translocation) [Legionella pneumophila str. Paris] >gi|54295733|ref|YP_128148.1| Putative TatA protein(twin arginine translocation) [Legionella pneumophila str. Lens] >gi|53755565|emb|CAH17064.1| Putative TatA protein(twin arginine translocation) [Legionella pneumophila str. Lens] >gi|53752691|emb|CAH14126.1| Putative TatA protein(twin arginine translocation) [Legionella pneumophila str. Paris] >gi|47132321|gb|AAT11788.1| TatA protein [Legionella pneumophila]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=61)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 96.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 41.5874)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.00751045)
        self.assertEqual(hsp.annotations["identity"], 27)
        self.assertEqual(hsp.annotations["positive"], 35)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 2, 49],
                          [14, 61]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 47))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDE'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({2: 'LSGISPLSLLLILAIIVALFGTSKLKTIGSDLGEAIKNFRKAMNSEE'}, length=61)",
        )
        self.assertEqual(hsp.target.id, "gi|54298906|ref|YP_125275.1|")
        self.assertEqual(hsp.target.name, "YP_125275")
        self.assertEqual(
            hsp.target.description,
            "Putative TatA protein(twin arginine translocation) [Legionella pneumophila str. Paris] >gi|54295733|ref|YP_128148.1| Putative TatA protein(twin arginine translocation) [Legionella pneumophila str. Lens] >gi|53755565|emb|CAH17064.1| Putative TatA protein(twin arginine translocation) [Legionella pneumophila str. Lens] >gi|53752691|emb|CAH14126.1| Putative TatA protein(twin arginine translocation) [Legionella pneumophila str. Paris] >gi|47132321|gb|AAT11788.1| TatA protein [Legionella pneumophila]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "+ GIS   LL+I  I+V LFGT KL +IGSDLG +IK F+KAM+ +E",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|542989         2 LSGISPLSLLLILAIIVALFGTSKLKTIGSDLGEAIKNFRKAMNSEE 49
                  0 ..|||..............|||.||..||||||..||.|.|||...| 47
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDE 61
""",
        )
        hit = hits[98]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|71363513|ref|ZP_00654157.1|")
        self.assertEqual(hit.target.name, "ZP_00654157")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Psychrobacter cryohalolentis K5] >gi|71161313|gb|EAO11109.1| Twin-arginine translocation protein TatA/E [Psychrobacter cryohalolentis K5]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=94)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 96.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 41.5874)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.00751045)
        self.assertEqual(hsp.annotations["identity"], 25)
        self.assertEqual(hsp.annotations["positive"], 39)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 55],
                          [14, 69]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 55))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGSFSITHWLILLVVVVVVFGTAKLKNAGKDLGGAVKGFKEAVKDEEAEHARKNR'}, length=94)",
        )
        self.assertEqual(hsp.target.id, "gi|71363513|ref|ZP_00654157.1|")
        self.assertEqual(hsp.target.name, "ZP_00654157")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Psychrobacter cryohalolentis K5] >gi|71161313|gb|EAO11109.1| Twin-arginine translocation protein TatA/E [Psychrobacter cryohalolentis K5]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG  SI   LI+ V+VV++FGT KL + G DLG ++KGFK+A+ D+E +  + ++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|713635         0 MGSFSITHWLILLVVVVVVFGTAKLKNAGKDLGGAVKGFKEAVKDEEAEHARKNR 55
                  0 ||..||.............|||.||...|.|||...||||.|..|.|........ 55
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQ 69
""",
        )
        hit = hits[99]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|71362217|ref|ZP_00653377.1|")
        self.assertEqual(hit.target.name, "ZP_00653377")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Psychrobacter cryohalolentis K5] >gi|71162335|gb|EAO12126.1| Twin-arginine translocation protein TatA/E [Psychrobacter cryohalolentis K5]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=80)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 96.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 41.5874)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.00751045)
        self.assertEqual(hsp.annotations["identity"], 27)
        self.assertEqual(hsp.annotations["positive"], 42)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 73],
                          [14, 87]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 73))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...TNQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGSFSITHWLILLVVVVVVFGTSKLRNAGKDLGGAVKGFKEAVKDENTEHAKKQ...TTK'}, length=80)",
        )
        self.assertEqual(hsp.target.id, "gi|71362217|ref|ZP_00653377.1|")
        self.assertEqual(hsp.target.name, "ZP_00653377")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Psychrobacter cryohalolentis K5] >gi|71162335|gb|EAO12126.1| Twin-arginine translocation protein TatA/E [Psychrobacter cryohalolentis K5]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG  SI   LI+ V+VV++FGT KL + G DLG ++KGFK+A+ D+  +  K     D       +++  T +",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|713622         0 MGSFSITHWLILLVVVVVVFGTSKLRNAGKDLGGAVKGFKEAVKDENTEHAKKQVVLDHD
                  0 ||..||.............|||.||...|.|||...||||.|..|......|.....|..
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|713622        60 GNAHTEERPTTTK 73
                 60 ..........|.. 73
gi|491764        74 AKTIADKQADTNQ 87
""",
        )
        hit = hits[100]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|27379862|ref|NP_771391.1|")
        self.assertEqual(hit.target.name, "NP_771391")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein bsl4751 [Bradyrhizobium japonicum USDA 110] >gi|27353015|dbj|BAC50016.1| bsl4751 [Bradyrhizobium japonicum USDA 110]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=77)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 95.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 41.2022)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.00980895)
        self.assertEqual(hsp.annotations["identity"], 27)
        self.assertEqual(hsp.annotations["positive"], 38)
        self.assertEqual(hsp.annotations["gaps"], 5)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 54, 54, 66],
                          [14, 68, 73, 85]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 71))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...ADT'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGSLSIWHWILVIAVVLLLFGRGKISDLMGDVAQGIKAFKKGMQDDDKPADKPE...APT'}, length=77)",
        )
        self.assertEqual(hsp.target.id, "gi|27379862|ref|NP_771391.1|")
        self.assertEqual(hsp.target.name, "NP_771391")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein bsl4751 [Bradyrhizobium japonicum USDA 110] >gi|27353015|dbj|BAC50016.1| bsl4751 [Bradyrhizobium japonicum USDA 110]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG +SIW  +++  +V+LLFG  K+  +  D+   IK FKK M DD+   DK        AK+I    A T",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|273798         0 MGSLSIWHWILVIAVVLLLFGRGKISDLMGDVAQGIKAFKKGMQDDDKPADKPE-----P
                  0 ||..|||............||..|......|....||.|||.|.||....||..-----.
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|273798        55 AKSIEHNAAPT 66
                 60 ||.|....|.| 71
gi|491764        74 AKTIADKQADT 85
""",
        )
        hit = hits[101]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|39935914|ref|NP_948190.1|")
        self.assertEqual(hit.target.name, "NP_948190")
        self.assertEqual(
            hit.target.description,
            "putative sec-independent protein translocase protein tatA/E [Rhodopseudomonas palustris CGA009] >gi|39649768|emb|CAE28290.1| putative sec-independent protein translocase protein tatA/E [Rhodopseudomonas palustris CGA009]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=78)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 94.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 40.817)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.0128109)
        self.assertEqual(hsp.annotations["identity"], 25)
        self.assertEqual(hsp.annotations["positive"], 41)
        self.assertEqual(hsp.annotations["gaps"], 2)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 54, 56, 76],
                          [14, 68, 68, 88]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 76))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...NQE'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGSLSIWHWIVVIAVVLLLFGRGKISDLMGDVAQGIKSFKKGLQDDEKTAEKSE...GSK'}, length=78)",
        )
        self.assertEqual(hsp.target.id, "gi|39935914|ref|NP_948190.1|")
        self.assertEqual(hsp.target.name, "NP_948190")
        self.assertEqual(
            hsp.target.description,
            "putative sec-independent protein translocase protein tatA/E [Rhodopseudomonas palustris CGA009] >gi|39649768|emb|CAE28290.1| putative sec-independent protein translocase protein tatA/E [Rhodopseudomonas palustris CGA009]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG +SIW  +++  +V+LLFG  K+  +  D+   IK FKK + DDE   +K+   +  D T+   A  + D   +",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|399359         0 MGSLSIWHWIVVIAVVLLLFGRGKISDLMGDVAQGIKSFKKGLQDDEKTAEKSEPVKSID
                  0 ||..|||............||..|......|....||.|||...|||....|..--...|
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS--QDAD

gi|399359        60 HTSTPGATNRTDVGSK 76
                 60 .|....|....|.... 76
gi|491764        72 FTAKTIADKQADTNQE 88
""",
        )
        hit = hits[102]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|17935600|ref|NP_532390.1|")
        self.assertEqual(hit.target.name, "NP_532390")
        self.assertEqual(
            hit.target.description,
            "SEC-independent protein translocase protein [Agrobacterium tumefaciens str. C58] >gi|17740143|gb|AAL42706.1| SEC-independent protein translocase protein [Agrobacterium tumefaciens str. C58] >gi|15156802|gb|AAK87479.1| AGR_C_3135p [Agrobacterium tumefaciens str. C58] >gi|15889013|ref|NP_354694.1| hypothetical protein AGR_C_3135 [Agrobacterium tumefaciens str. C58] >gi|25520668|pir||F97565 hypothetical protein AGR_C_3135 [imported] - Agrobacterium tumefaciens (strain C58, Cereon) >gi|25522636|pir||AD2786 SEC-independent protein translocase protein [imported] - Agrobacterium tumefaciens (strain C58, Dupont) >gi|24212489|sp|Q8UEP9|TATA_AGRT5 Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=70)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 94.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 40.817)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.0128109)
        self.assertEqual(hsp.annotations["identity"], 30)
        self.assertEqual(hsp.annotations["positive"], 43)
        self.assertEqual(hsp.annotations["gaps"], 3)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 47, 47, 62, 62, 67],
                          [14, 61, 63, 78, 79, 84]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 70))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...QAD'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGSFSVWHWLIVLVIVLVLFGRGKIPELMGDVAKGIKSFKKGMADEDQTPPPAD...KAD'}, length=70)",
        )
        self.assertEqual(hsp.target.id, "gi|17935600|ref|NP_532390.1|")
        self.assertEqual(hsp.target.name, "NP_532390")
        self.assertEqual(
            hsp.target.description,
            "SEC-independent protein translocase protein [Agrobacterium tumefaciens str. C58] >gi|17740143|gb|AAL42706.1| SEC-independent protein translocase protein [Agrobacterium tumefaciens str. C58] >gi|15156802|gb|AAK87479.1| AGR_C_3135p [Agrobacterium tumefaciens str. C58] >gi|15889013|ref|NP_354694.1| hypothetical protein AGR_C_3135 [Agrobacterium tumefaciens str. C58] >gi|25520668|pir||F97565 hypothetical protein AGR_C_3135 [imported] - Agrobacterium tumefaciens (strain C58, Cereon) >gi|25522636|pir||AD2786 SEC-independent protein translocase protein [imported] - Agrobacterium tumefaciens (strain C58, Dupont) >gi|24212489|sp|Q8UEP9|TATA_AGRT5 Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG  S+W  LI+ VIV++LFG  K+  +  D+   IK FKK M+D++  Q     DA+  AKT+ D +AD",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|179356         0 MGSFSVWHWLIVLVIVLVLFGRGKIPELMGDVAKGIKSFKKGMADED--QTPPPADANAN
                  0 ||..|.|............||..|......|....||.|||.|.|..--|.....||...
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|179356        58 AKTV-DHKAD 67
                 60 |||.-|..|| 70
gi|491764        74 AKTIADKQAD 84
""",
        )
        hit = hits[103]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|62289827|ref|YP_221620.1|")
        self.assertEqual(hit.target.name, "YP_221620")
        self.assertEqual(
            hit.target.description,
            "Sec-independent protein translocase protein TatA, hypothetical [Brucella abortus biovar 1 str. 9-941] >gi|17983053|gb|AAL52265.1| SEC-INDEPENDENT PROTEIN TRANSLOCASE PROTEIN TATA [Brucella melitensis 16M] >gi|62195959|gb|AAX74259.1| Sec-independent protein translocase protein TatA, hypothetical [Brucella abortus biovar 1 str. 9-941] >gi|17987367|ref|NP_540001.1| SEC-INDEPENDENT PROTEIN TRANSLOCASE PROTEIN TATA [Brucella melitensis 16M] >gi|25302679|pir||AF3387 sec-independent protein translocase protein tatA [imported] - Brucella melitensis (strain 16M) >gi|54039707|sp|P66888|TATA_BRUSU Sec-independent protein translocase protein tatA/E homolog >gi|54042236|sp|P66887|TATA_BRUME Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=72)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 93.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 40.4318)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.0167315)
        self.assertEqual(hsp.annotations["identity"], 24)
        self.assertEqual(hsp.annotations["positive"], 36)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 57],
                          [14, 71]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 57))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDA'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGSFSIWHWLIVLAVVLLLFGRGKIPELMGDVAKGIKNFKQGMADEDAKEDPRTIDA'}, length=72)",
        )
        self.assertEqual(hsp.target.id, "gi|62289827|ref|YP_221620.1|")
        self.assertEqual(hsp.target.name, "YP_221620")
        self.assertEqual(
            hsp.target.description,
            "Sec-independent protein translocase protein TatA, hypothetical [Brucella abortus biovar 1 str. 9-941] >gi|17983053|gb|AAL52265.1| SEC-INDEPENDENT PROTEIN TRANSLOCASE PROTEIN TATA [Brucella melitensis 16M] >gi|62195959|gb|AAX74259.1| Sec-independent protein translocase protein TatA, hypothetical [Brucella abortus biovar 1 str. 9-941] >gi|17987367|ref|NP_540001.1| SEC-INDEPENDENT PROTEIN TRANSLOCASE PROTEIN TATA [Brucella melitensis 16M] >gi|25302679|pir||AF3387 sec-independent protein translocase protein tatA [imported] - Brucella melitensis (strain 16M) >gi|54039707|sp|P66888|TATA_BRUSU Sec-independent protein translocase protein tatA/E homolog >gi|54042236|sp|P66887|TATA_BRUME Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG  SIW  LI+  +V+LLFG  K+  +  D+   IK FK+ M+D++ K+D  + DA",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|622898         0 MGSFSIWHWLIVLAVVLLLFGRGKIPELMGDVAKGIKNFKQGMADEDAKEDPRTIDA 57
                  0 ||..|||............||..|......|....||.||..|.|...|.|....|| 57
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDA 71
""",
        )
        hit = hits[104]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|23347697|gb|AAN29810.1|")
        self.assertEqual(hit.target.name, "AAN29810")
        self.assertEqual(
            hit.target.description,
            "Sec-independent protein translocase protein TatA, putative [Brucella suis 1330] >gi|23501768|ref|NP_697895.1| Sec-independent protein translocase protein TatA, putative [Brucella suis 1330]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=80)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 93.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 40.4318)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.0167315)
        self.assertEqual(hsp.annotations["identity"], 24)
        self.assertEqual(hsp.annotations["positive"], 36)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 8, 65],
                          [14, 71]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 57))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDA'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({8: 'MGSFSIWHWLIVLAVVLLLFGRGKIPELMGDVAKGIKNFKQGMADEDAKEDPRTIDA'}, length=80)",
        )
        self.assertEqual(hsp.target.id, "gi|23347697|gb|AAN29810.1|")
        self.assertEqual(hsp.target.name, "AAN29810")
        self.assertEqual(
            hsp.target.description,
            "Sec-independent protein translocase protein TatA, putative [Brucella suis 1330] >gi|23501768|ref|NP_697895.1| Sec-independent protein translocase protein TatA, putative [Brucella suis 1330]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG  SIW  LI+  +V+LLFG  K+  +  D+   IK FK+ M+D++ K+D  + DA",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|233476         8 MGSFSIWHWLIVLAVVLLLFGRGKIPELMGDVAKGIKNFKQGMADEDAKEDPRTIDA 65
                  0 ||..|||............||..|......|....||.||..|.|...|.|....|| 57
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDA 71
""",
        )
        hit = hits[105]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|75675971|ref|YP_318392.1|")
        self.assertEqual(hit.target.name, "YP_318392")
        self.assertEqual(
            hit.target.description,
            "twin-arginine translocation protein TatA/E [Nitrobacter winogradskyi Nb-255] >gi|74420841|gb|ABA05040.1| twin-arginine translocation protein TatA/E [Nitrobacter winogradskyi Nb-255]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=78)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 93.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 40.4318)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.0167315)
        self.assertEqual(hsp.annotations["identity"], 21)
        self.assertEqual(hsp.annotations["positive"], 32)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 52],
                          [14, 66]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 52))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDK'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGSLSIWHWIVVVAVILLLFGRGKISDLMGDVAQGIKAFKKGMKDDEKTAEK'}, length=78)",
        )
        self.assertEqual(hsp.target.id, "gi|75675971|ref|YP_318392.1|")
        self.assertEqual(hsp.target.name, "YP_318392")
        self.assertEqual(
            hsp.target.description,
            "twin-arginine translocation protein TatA/E [Nitrobacter winogradskyi Nb-255] >gi|74420841|gb|ABA05040.1| twin-arginine translocation protein TatA/E [Nitrobacter winogradskyi Nb-255]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG +SIW  +++  +++LLFG  K+  +  D+   IK FKK M DDE   +K",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|756759         0 MGSLSIWHWIVVVAVILLLFGRGKISDLMGDVAQGIKAFKKGMKDDEKTAEK 52
                  0 ||..|||............||..|......|....||.|||.|.|||....| 52
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDK 66
""",
        )
        hit = hits[106]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|69928230|ref|ZP_00625391.1|")
        self.assertEqual(hit.target.name, "ZP_00625391")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Nitrobacter hamburgensis X14] >gi|69143234|gb|EAN61762.1| Twin-arginine translocation protein TatA/E [Nitrobacter hamburgensis X14]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=79)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 93.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 40.4318)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.0167315)
        self.assertEqual(hsp.annotations["identity"], 21)
        self.assertEqual(hsp.annotations["positive"], 32)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 52],
                          [14, 66]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 52))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDK'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGSLSIWHWIVVIAVILLLFGRGKISDLMGDVAQGIKAFKKGMQDDEKTAEK'}, length=79)",
        )
        self.assertEqual(hsp.target.id, "gi|69928230|ref|ZP_00625391.1|")
        self.assertEqual(hsp.target.name, "ZP_00625391")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Nitrobacter hamburgensis X14] >gi|69143234|gb|EAN61762.1| Twin-arginine translocation protein TatA/E [Nitrobacter hamburgensis X14]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG +SIW  +++  +++LLFG  K+  +  D+   IK FKK M DDE   +K",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|699282         0 MGSLSIWHWIVVIAVILLLFGRGKISDLMGDVAQGIKAFKKGMQDDEKTAEK 52
                  0 ||..|||............||..|......|....||.|||.|.|||....| 52
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDK 66
""",
        )
        hit = hits[107]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|77689454|ref|ZP_00804635.1|")
        self.assertEqual(hit.target.name, "ZP_00804635")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Rhodopseudomonas palustris BisB5] >gi|77653821|gb|EAO85595.1| Twin-arginine translocation protein TatA/E [Rhodopseudomonas palustris BisB5]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=79)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 92.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 40.0466)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.0218521)
        self.assertEqual(hsp.annotations["identity"], 25)
        self.assertEqual(hsp.annotations["positive"], 37)
        self.assertEqual(hsp.annotations["gaps"], 2)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 54, 56, 67],
                          [14, 68, 68, 79]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 67))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...TIA'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGSLSIWHWIVVIAVVLLLFGRGKISDLMGDVAQGIKSFKKGLQDDEKTAEKPD...TAA'}, length=79)",
        )
        self.assertEqual(hsp.target.id, "gi|77689454|ref|ZP_00804635.1|")
        self.assertEqual(hsp.target.name, "ZP_00804635")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Rhodopseudomonas palustris BisB5] >gi|77653821|gb|EAO85595.1| Twin-arginine translocation protein TatA/E [Rhodopseudomonas palustris BisB5]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG +SIW  +++  +V+LLFG  K+  +  D+   IK FKK + DDE   +K    +  D  A T A",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|776894         0 MGSLSIWHWIVVIAVVLLLFGRGKISDLMGDVAQGIKSFKKGLQDDEKTAEKPDPVKSID
                  0 ||..|||............||..|......|....||.|||...|||....|..--...|
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS--QDAD

gi|776894        60 HNAPTAA 67
                 60 ..|.|.| 67
gi|491764        72 FTAKTIA 79
""",
        )
        hit = hits[108]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|77743614|ref|ZP_00812071.1|")
        self.assertEqual(hit.target.name, "ZP_00812071")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Rhodopseudomonas palustris BisA53] >gi|77696516|gb|EAO87698.1| Twin-arginine translocation protein TatA/E [Rhodopseudomonas palustris BisA53]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=78)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 91.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 39.6614)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.0285397)
        self.assertEqual(hsp.annotations["identity"], 23)
        self.assertEqual(hsp.annotations["positive"], 34)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 53],
                          [14, 67]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 53))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKT'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGSLSIWHWLVVIVVVLLLFGRGKISDLMGDVAQGIKSFKKGLQDDDKTAEKT'}, length=78)",
        )
        self.assertEqual(hsp.target.id, "gi|77743614|ref|ZP_00812071.1|")
        self.assertEqual(hsp.target.name, "ZP_00812071")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Rhodopseudomonas palustris BisA53] >gi|77696516|gb|EAO87698.1| Twin-arginine translocation protein TatA/E [Rhodopseudomonas palustris BisA53]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG +SIW  L++ V+V+LLFG  K+  +  D+   IK FKK + DD+   +KT",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|777436         0 MGSLSIWHWLVVIVVVLLLFGRGKISDLMGDVAQGIKSFKKGLQDDDKTAEKT 53
                  0 ||..|||............||..|......|....||.|||...||.....|| 53
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKT 67
""",
        )
        hit = hits[109]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|71066141|ref|YP_264868.1|")
        self.assertEqual(hit.target.name, "YP_264868")
        self.assertEqual(
            hit.target.description,
            "twin-arginine translocation protein TatE [Psychrobacter arcticus 273-4] >gi|71039126|gb|AAZ19434.1| twin-arginine translocation protein TatE [Psychrobacter arcticus 273-4]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=89)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 90.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 39.2762)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.037274)
        self.assertEqual(hsp.annotations["identity"], 16)
        self.assertEqual(hsp.annotations["positive"], 22)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[18, 46],
                          [33, 61]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 28))
        self.assertEqual(
            repr(hsp.query.seq), "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDE'}, length=103)"
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq), "Seq({18: 'FGTAKLKNAGKDLGGAVKGFKEAVKDEE'}, length=89)"
        )
        self.assertEqual(hsp.target.id, "gi|71066141|ref|YP_264868.1|")
        self.assertEqual(hsp.target.name, "YP_264868")
        self.assertEqual(
            hsp.target.description,
            "twin-arginine translocation protein TatE [Psychrobacter arcticus 273-4] >gi|71039126|gb|AAZ19434.1| twin-arginine translocation protein TatE [Psychrobacter arcticus 273-4]",
        )
        self.assertEqual(hsp.annotations["midline"], "FGT KL + G DLG ++KGFK+A+ D+E")
        self.assertEqual(
            str(hsp),
            """\
gi|710661        18 FGTAKLKNAGKDLGGAVKGFKEAVKDEE 46
                  0 |||.||...|.|||...||||.|..|.| 28
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDE 61
""",
        )
        hit = hits[110]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|28199457|ref|NP_779771.1|")
        self.assertEqual(hit.target.name, "NP_779771")
        self.assertEqual(
            hit.target.description,
            "SEC-independent protein translocase [Xylella fastidiosa Temecula1] >gi|71900614|ref|ZP_00682740.1| Twin-arginine translocation protein TatA/E [Xylella fastidiosa Ann-1] >gi|71898114|ref|ZP_00680300.1| Twin-arginine translocation protein TatA/E [Xylella fastidiosa Ann-1] >gi|71732088|gb|EAO34144.1| Twin-arginine translocation protein TatA/E [Xylella fastidiosa Ann-1] >gi|71729608|gb|EAO31713.1| Twin-arginine translocation protein TatA/E [Xylella fastidiosa Ann-1] >gi|71274575|ref|ZP_00650863.1| Twin-arginine translocation protein TatA/E [Xylella fastidiosa Dixon] >gi|71164307|gb|EAO14021.1| Twin-arginine translocation protein TatA/E [Xylella fastidiosa Dixon] >gi|28057572|gb|AAO29420.1| SEC-independent protein translocase [Xylella fastidiosa Temecula1] >gi|32130179|sp|Q87B80|TATA_XYLFT Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=71)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 89.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 38.891)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.0486813)
        self.assertEqual(hsp.annotations["identity"], 24)
        self.assertEqual(hsp.annotations["positive"], 40)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 61],
                          [14, 75]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 61))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...FTA'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGSFSLLHWLVVLVIVLLVFGTKRLANGAKDIGSAIKEFKKSLHEDDKPTDQLG...STA'}, length=71)",
        )
        self.assertEqual(hsp.target.id, "gi|28199457|ref|NP_779771.1|")
        self.assertEqual(hsp.target.name, "NP_779771")
        self.assertEqual(
            hsp.target.description,
            "SEC-independent protein translocase [Xylella fastidiosa Temecula1] >gi|71900614|ref|ZP_00682740.1| Twin-arginine translocation protein TatA/E [Xylella fastidiosa Ann-1] >gi|71898114|ref|ZP_00680300.1| Twin-arginine translocation protein TatA/E [Xylella fastidiosa Ann-1] >gi|71732088|gb|EAO34144.1| Twin-arginine translocation protein TatA/E [Xylella fastidiosa Ann-1] >gi|71729608|gb|EAO31713.1| Twin-arginine translocation protein TatA/E [Xylella fastidiosa Ann-1] >gi|71274575|ref|ZP_00650863.1| Twin-arginine translocation protein TatA/E [Xylella fastidiosa Dixon] >gi|71164307|gb|EAO14021.1| Twin-arginine translocation protein TatA/E [Xylella fastidiosa Dixon] >gi|28057572|gb|AAO29420.1| SEC-independent protein translocase [Xylella fastidiosa Temecula1] >gi|32130179|sp|Q87B80|TATA_XYLFT Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG  S+   L++ VIV+L+FGTK+L +   D+G++IK FKK++ +D+   D+    +  TA",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|281994         0 MGSFSLLHWLVVLVIVLLVFGTKRLANGAKDIGSAIKEFKKSLHEDDKPTDQLGSTSQST
                  0 ||..|..............||||.|.....|.|..||.|||....|....|........|
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|281994        60 A 61
                 60 | 61
gi|491764        74 A 75
""",
        )
        hit = hits[111]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|15837166|ref|NP_297854.1|")
        self.assertEqual(hit.target.name, "NP_297854")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein XF0564 [Xylella fastidiosa 9a5c] >gi|9105426|gb|AAF83374.1| conserved hypothetical protein [Xylella fastidiosa 9a5c] >gi|11360714|pir||B82791 conserved hypothetical protein XF0564 [imported] - Xylella fastidiosa (strain 9a5c) >gi|9979030|sp|Q9PFU3|TATA_XYLFA Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=71)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 89.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 38.891)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.0486813)
        self.assertEqual(hsp.annotations["identity"], 24)
        self.assertEqual(hsp.annotations["positive"], 40)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 61],
                          [14, 75]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 61))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...FTA'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGSFSLLHWLVVLVIVLLVFGTKRLANGAKDIGSAIKEFKKSLREDDKPTDQLG...STA'}, length=71)",
        )
        self.assertEqual(hsp.target.id, "gi|15837166|ref|NP_297854.1|")
        self.assertEqual(hsp.target.name, "NP_297854")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein XF0564 [Xylella fastidiosa 9a5c] >gi|9105426|gb|AAF83374.1| conserved hypothetical protein [Xylella fastidiosa 9a5c] >gi|11360714|pir||B82791 conserved hypothetical protein XF0564 [imported] - Xylella fastidiosa (strain 9a5c) >gi|9979030|sp|Q9PFU3|TATA_XYLFA Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG  S+   L++ VIV+L+FGTK+L +   D+G++IK FKK++ +D+   D+    +  TA",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|158371         0 MGSFSLLHWLVVLVIVLLVFGTKRLANGAKDIGSAIKEFKKSLREDDKPTDQLGSTSQST
                  0 ||..|..............||||.|.....|.|..||.|||....|....|........|
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|158371        60 A 61
                 60 | 61
gi|491764        74 A 75
""",
        )
        hit = hits[112]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|15074462|emb|CAC46108.1|")
        self.assertEqual(hit.target.name, "CAC46108")
        self.assertEqual(
            hit.target.description,
            "HYPOTHETICAL TRANSMEMBRANE PROTEIN [Sinorhizobium meliloti] >gi|15965282|ref|NP_385635.1| HYPOTHETICAL TRANSMEMBRANE PROTEIN [Sinorhizobium meliloti 1021] >gi|24212503|sp|Q92Q25|TATA_RHIME Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=68)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 89.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 38.891)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.0486813)
        self.assertEqual(hsp.annotations["identity"], 28)
        self.assertEqual(hsp.annotations["positive"], 35)
        self.assertEqual(hsp.annotations["gaps"], 4)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 47, 47, 63],
                          [14, 61, 65, 81]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 67))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...ADK'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGSFSIWHWLIVLAVVLLLFGRGKIPELMGDVAKGIKNFKKGMGDDEVASADKS...DHK'}, length=68)",
        )
        self.assertEqual(hsp.target.id, "gi|15074462|emb|CAC46108.1|")
        self.assertEqual(hsp.target.name, "CAC46108")
        self.assertEqual(
            hsp.target.description,
            "HYPOTHETICAL TRANSMEMBRANE PROTEIN [Sinorhizobium meliloti] >gi|15965282|ref|NP_385635.1| HYPOTHETICAL TRANSMEMBRANE PROTEIN [Sinorhizobium meliloti 1021] >gi|24212503|sp|Q92Q25|TATA_RHIME Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG  SIW  LI+  +V+LLFG  K+  +  D+   IK FKK M DDE      S D     KT+  K",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|150744         0 MGSFSIWHWLIVLAVVLLLFGRGKIPELMGDVAKGIKNFKKGMGDDE----VASADKSVD
                  0 ||..|||............||..|......|....||.|||.|.|||----..|.|....
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|150744        56 GKTVDHK 63
                 60 .||...| 67
gi|491764        74 AKTIADK 81
""",
        )
        hit = hits[113]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|27462871|gb|AAO15625.1|")
        self.assertEqual(hit.target.name, "AAO15625")
        self.assertEqual(
            hit.target.description,
            "Sec-independent protein translocase TatA [Rhizobium leguminosarum bv. viciae]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=63)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 87.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 38.1206)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.0830377)
        self.assertEqual(hsp.annotations["identity"], 30)
        self.assertEqual(hsp.annotations["positive"], 40)
        self.assertEqual(hsp.annotations["gaps"], 9)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 45, 45, 62],
                          [14, 59, 68, 85]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 71))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...ADT'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGSFSMWHWLIVLVIVLLLFGRGKIPELMGDVAKGIKSFKKGMTDEDAPDTAKT...DET'}, length=63)",
        )
        self.assertEqual(hsp.target.id, "gi|27462871|gb|AAO15625.1|")
        self.assertEqual(hsp.target.name, "AAO15625")
        self.assertEqual(
            hsp.target.description,
            "Sec-independent protein translocase TatA [Rhizobium leguminosarum bv. viciae]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG  S+W  LI+ VIV+LLFG  K+  +  D+   IK FKK M+D         +DA  TAKT+  K  +T",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|274628         0 MGSFSMWHWLIVLVIVLLLFGRGKIPELMGDVAKGIKSFKKGMTD---------EDAPDT
                  0 ||..|.|............||..|......|....||.|||.|.|---------.||..|
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|274628        51 AKTVDHKADET 62
                 60 |||...|...| 71
gi|491764        74 AKTIADKQADT 85
""",
        )
        hit = hits[114]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|35211273|dbj|BAC88652.1|")
        self.assertEqual(hit.target.name, "BAC88652")
        self.assertEqual(
            hit.target.description,
            "gsl0711 [Gloeobacter violaceus PCC 7421] >gi|37520280|ref|NP_923657.1| hypothetical protein gsl0711 [Gloeobacter violaceus PCC 7421]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=72)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 87.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 38.1206)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.0830377)
        self.assertEqual(hsp.annotations["identity"], 17)
        self.assertEqual(hsp.annotations["positive"], 24)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[20, 57],
                          [33, 70]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 37))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQD'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({20: 'FGPKKLPEMGSALGKAIRGFKSGVSDEPAPQQSASKE'}, length=72)",
        )
        self.assertEqual(hsp.target.id, "gi|35211273|dbj|BAC88652.1|")
        self.assertEqual(hsp.target.name, "BAC88652")
        self.assertEqual(
            hsp.target.description,
            "gsl0711 [Gloeobacter violaceus PCC 7421] >gi|37520280|ref|NP_923657.1| hypothetical protein gsl0711 [Gloeobacter violaceus PCC 7421]",
        )
        self.assertEqual(
            hsp.annotations["midline"], "FG KKL  +GS LG +I+GFK  +SD+   Q   S++"
        )
        self.assertEqual(
            str(hsp),
            """\
gi|352112        20 FGPKKLPEMGSALGKAIRGFKSGVSDEPAPQQSASKE 57
                  0 ||.|||...||.||..|.|||...||....|...|.. 37
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQD 70
""",
        )
        hit = hits[115]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|34482347|emb|CAE09348.1|")
        self.assertEqual(hit.target.name, "CAE09348")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein [Wolinella succinogenes] >gi|34556633|ref|NP_906448.1| hypothetical protein WS0187 [Wolinella succinogenes DSM 1740]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=80)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 86.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 37.7354)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.108451)
        self.assertEqual(hsp.annotations["identity"], 18)
        self.assertEqual(hsp.annotations["positive"], 26)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[19, 71],
                          [33, 85]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 52))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADT'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({19: 'FGAKKIPDLAKGLGSGIKNFKKAMKEDEESEVSVTKSEKIEEAPKNEKSAST'}, length=80)",
        )
        self.assertEqual(hsp.target.id, "gi|34482347|emb|CAE09348.1|")
        self.assertEqual(hsp.target.name, "CAE09348")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein [Wolinella succinogenes] >gi|34556633|ref|NP_906448.1| hypothetical protein WS0187 [Wolinella succinogenes DSM 1740]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "FG KK+  +   LG+ IK FKKAM +DE  +   ++          +K A T",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|344823        19 FGAKKIPDLAKGLGSGIKNFKKAMKEDEESEVSVTKSEKIEEAPKNEKSAST 71
                  0 ||.||.......||..||.|||||..||...................|.|.| 52
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADT 85
""",
        )
        hit = hits[116]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|32262257|gb|AAP77305.1|")
        self.assertEqual(hit.target.name, "AAP77305")
        self.assertEqual(
            hit.target.description,
            "component of Sec-independent protein secretion pathway TatA [Helicobacter hepaticus ATCC 51449] >gi|32266207|ref|NP_860239.1| component of Sec-independent protein secretion pathway TatA [Helicobacter hepaticus ATCC 51449]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=82)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 86.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 37.7354)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.108451)
        self.assertEqual(hsp.annotations["identity"], 22)
        self.assertEqual(hsp.annotations["positive"], 37)
        self.assertEqual(hsp.annotations["gaps"], 6)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[19, 47, 47, 54, 54, 74, 75, 81],
                          [33, 61, 63, 70, 73, 93, 93, 99]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 67))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQ...RHD'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({19: 'FGAKKIPELAKGLGSGIKNFKKAVKEDEEDNQSEENTKSQIKQSESKNENVSKT...KQD'}, length=82)",
        )
        self.assertEqual(hsp.target.id, "gi|32262257|gb|AAP77305.1|")
        self.assertEqual(hsp.target.name, "AAP77305")
        self.assertEqual(
            hsp.target.description,
            "component of Sec-independent protein secretion pathway TatA [Helicobacter hepaticus ATCC 51449] >gi|32266207|ref|NP_860239.1| component of Sec-independent protein secretion pathway TatA [Helicobacter hepaticus ATCC 51449]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "FG KK+  +   LG+ IK FKKA+ +DE  +D  S++     K+   +    N+  +KT  D+++ D",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|322622        19 FGAKKIPELAKGLGSGIKNFKKAVKEDE--EDNQSEE---NTKSQIKQSESKNENVSKTH
                  0 ||.||.......||..||.||||...||--.|..|..---..|.........|....||.
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTE

gi|322622        74 TDSQKQD 81
                 60 -|....| 67
gi|491764        93 -DAKRHD 99
""",
        )
        hit = hits[117]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|76261408|ref|ZP_00769019.1|")
        self.assertEqual(hit.target.name, "ZP_00769019")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Chloroflexus aurantiacus J-10-fl] >gi|76163681|gb|EAO57850.1| Twin-arginine translocation protein TatA/E [Chloroflexus aurantiacus J-10-fl]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=62)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 85.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 37.3502)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.141641)
        self.assertEqual(hsp.annotations["identity"], 25)
        self.assertEqual(hsp.annotations["positive"], 39)
        self.assertEqual(hsp.annotations["gaps"], 1)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 1, 44, 45, 57],
                          [14, 57, 57, 69]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 56))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({1: 'IGGLGWGELLIILIIVIAIFGAGKLAGLGGALGSSIREFRKAVKGDDEPRSDAKTE'}, length=62)",
        )
        self.assertEqual(hsp.target.id, "gi|76261408|ref|ZP_00769019.1|")
        self.assertEqual(hsp.target.name, "ZP_00769019")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Chloroflexus aurantiacus J-10-fl] >gi|76163681|gb|EAO57850.1| Twin-arginine translocation protein TatA/E [Chloroflexus aurantiacus J-10-fl]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "+GG+   +LLII +IV+ +FG  KL  +G  LG+SI+ F+KA+  DDEP+ D  ++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|762614         1 IGGLGWGELLIILIIVIAIFGAGKLAGLGGALGSSIREFRKAVKGDDEPRSDAKTE 57
                  0 .||................||..||...|..||.||..|.||.-.||||..|.... 56
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAM-SDDEPKQDKTSQ 69
""",
        )
        hit = hits[118]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|69933726|ref|ZP_00628928.1|")
        self.assertEqual(hit.target.name, "ZP_00628928")
        self.assertEqual(
            hit.target.description,
            "sec-independent translocation protein mttA/Hcf106 [Paracoccus denitrificans PD1222] >gi|69155362|gb|EAN68465.1| sec-independent translocation protein mttA/Hcf106 [Paracoccus denitrificans PD1222]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=159)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 85.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 37.3502)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.141641)
        self.assertEqual(hsp.annotations["identity"], 22)
        self.assertEqual(hsp.annotations["positive"], 34)
        self.assertEqual(hsp.annotations["gaps"], 4)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 87, 125, 125, 134, 135, 149],
                          [ 33,  71,  74,  83,  83,  97]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 65))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQ...AKR'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({87: 'FGRGKVSSLMGEVGKGITAFKKGVKDGAEEYDRAAEDAKKPIGQSPAEDATRAQ...AMR'}, length=159)",
        )
        self.assertEqual(hsp.target.id, "gi|69933726|ref|ZP_00628928.1|")
        self.assertEqual(hsp.target.name, "ZP_00628928")
        self.assertEqual(
            hsp.target.description,
            "sec-independent translocation protein mttA/Hcf106 [Paracoccus denitrificans PD1222] >gi|69155362|gb|EAN68465.1| sec-independent translocation protein mttA/Hcf106 [Paracoccus denitrificans PD1222]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "FG  K+ S+  ++G  I  FKK + D   + D+ ++DA    K I    A D  + QA+T +A R",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|699337        87 FGRGKVSSLMGEVGKGITAFKKGVKDGAEEYDRAAEDA---KKPIGQSPAEDATRAQAQT
                  0 ||..|..|.....|..|..|||...|.....|....||---.|.|....|-|....||.|
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQA-DTNQEQAKT

gi|699337       144 NEAMR 149
                 60 ..|.|  65
gi|491764        92 EDAKR  97
""",
        )
        hit = hits[119]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|15605662|ref|NP_213037.1|")
        self.assertEqual(hit.target.name, "NP_213037")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein aq_064b [Aquifex aeolicus VF5] >gi|2982827|gb|AAC06450.1| hypothetical protein [Aquifex aeolicus VF5] >gi|7444819|pir||B70306 conserved hypothetical protein aq_064b - Aquifex aeolicus >gi|9978994|sp|O66477|TAT2_AQUAE Sec-independent protein translocase tatA/E-like protein 2",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=59)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 84.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 36.965)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.184989)
        self.assertEqual(hsp.annotations["identity"], 19)
        self.assertEqual(hsp.annotations["positive"], 33)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 6, 52],
                          [20, 66]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 46))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({20: 'WQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDK'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({6: 'WQLILILLVILVIFGASKLPEVGKGLGEGIRNFKKALSGEEEEKGK'}, length=59)",
        )
        self.assertEqual(hsp.target.id, "gi|15605662|ref|NP_213037.1|")
        self.assertEqual(hsp.target.name, "NP_213037")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein aq_064b [Aquifex aeolicus VF5] >gi|2982827|gb|AAC06450.1| hypothetical protein [Aquifex aeolicus VF5] >gi|7444819|pir||B70306 conserved hypothetical protein aq_064b - Aquifex aeolicus >gi|9978994|sp|O66477|TAT2_AQUAE Sec-independent protein translocase tatA/E-like protein 2",
        )
        self.assertEqual(
            hsp.annotations["midline"], "WQL++I ++++++FG  KL  +G  LG  I+ FKKA+S +E ++ K"
        )
        self.assertEqual(
            str(hsp),
            """\
gi|156056         6 WQLILILLVILVIFGASKLPEVGKGLGEGIRNFKKALSGEEEEKGK 52
                  0 ||...........||..||...|..||..|..||||.|..|....| 46
gi|491764        20 WQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDK 66
""",
        )
        hit = hits[120]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|68538777|ref|ZP_00578553.1|")
        self.assertEqual(hit.target.name, "ZP_00578553")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Sphingopyxis alaskensis RB2256] >gi|68523891|gb|EAN47017.1| Twin-arginine translocation protein TatA/E [Sphingopyxis alaskensis RB2256]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=77)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 83.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 36.5798)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.241603)
        self.assertEqual(hsp.annotations["identity"], 25)
        self.assertEqual(hsp.annotations["positive"], 37)
        self.assertEqual(hsp.annotations["gaps"], 5)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 47, 52, 66],
                          [14, 61, 61, 75]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 66))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...FTA'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGSFSIWHWLVVGILVLLLFGKGRFSDMMGDVAKGIKSFKKGMSEDDAPTPAPK...LSA'}, length=77)",
        )
        self.assertEqual(hsp.target.id, "gi|68538777|ref|ZP_00578553.1|")
        self.assertEqual(hsp.target.name, "ZP_00578553")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Sphingopyxis alaskensis RB2256] >gi|68523891|gb|EAN47017.1| Twin-arginine translocation protein TatA/E [Sphingopyxis alaskensis RB2256]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG  SIW  L++ ++V+LLFG  +   +  D+   IK FKK MS+D+     PKQ    +  D +A",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|685387         0 MGSFSIWHWLVVGILVLLLFGKGRFSDMMGDVAKGIKSFKKGMSEDDAPTPAPKQIDAQR
                  0 ||..|||............||.........|....||.|||.||.|.-----|||.....
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDE-----PKQDKTSQ

gi|685387        60 APDLSA 66
                 60 ..|..| 66
gi|491764        69 DADFTA 75
""",
        )
        hit = hits[121]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|68136098|ref|ZP_00544086.1|")
        self.assertEqual(hit.target.name, "ZP_00544086")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Methanosarcina barkeri str. fusaro] >gi|72397188|gb|AAZ71461.1| sec-independent protein translocase, protein [Methanosarcina barkeri str. fusaro] >gi|73670026|ref|YP_306041.1| sec-independent protein translocase, protein [Methanosarcina barkeri str. fusaro]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=130)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 82.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 36.1946)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.315543)
        self.assertEqual(hsp.annotations["identity"], 19)
        self.assertEqual(hsp.annotations["positive"], 34)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 17,  84],
                          [ 33, 100]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 67))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQ...HDK'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({17: 'FGPKKLPELARSLGSAMGEFKKAQRASELELNQFDAYTRKAANTVASEEKEKGK...EAK'}, length=130)",
        )
        self.assertEqual(hsp.target.id, "gi|68136098|ref|ZP_00544086.1|")
        self.assertEqual(hsp.target.name, "ZP_00544086")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Methanosarcina barkeri str. fusaro] >gi|72397188|gb|AAZ71461.1| sec-independent protein translocase, protein [Methanosarcina barkeri str. fusaro] >gi|73670026|ref|YP_306041.1| sec-independent protein translocase, protein [Methanosarcina barkeri str. fusaro]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "FG KKL  +   LG+++  FKKA    E + ++        A T+A ++ +  +EQ +T  + +  K",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|681360        17 FGPKKLPELARSLGSAMGEFKKAQRASELELNQFDAYTRKAANTVASEEKEKGKEQKETP
                  0 ||.|||......||.....||||....|.............|.|.|........||..|.
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTE

gi|681360        77 SSSKEAK  84
                 60 ......|  67
gi|491764        93 DAKRHDK 100
""",
        )
        hit = hits[122]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|20259265|gb|AAM14368.1|")
        self.assertEqual(hit.target.name, "AAM14368")
        self.assertEqual(
            hit.target.description,
            "putative Tha4 protein [Arabidopsis thaliana] >gi|15810475|gb|AAL07125.1| putative Tha4 protein [Arabidopsis thaliana] >gi|11762160|gb|AAG40358.1| AT5g28750 [Arabidopsis thaliana] >gi|7682781|gb|AAF67362.1| Hypothetical protein T32B20.e [Arabidopsis thaliana] >gi|15241912|ref|NP_198227.1| thylakoid assembly protein, putative [Arabidopsis thaliana]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=147)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 81.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 35.8094)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.412112)
        self.assertEqual(hsp.annotations["identity"], 25)
        self.assertEqual(hsp.annotations["positive"], 46)
        self.assertEqual(hsp.annotations["gaps"], 4)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 64, 126, 126, 147],
                          [ 16,  78,  82, 103]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 87))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({16: 'GISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQD...EQV'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({64: 'GLGVPELAVIAGVAALLFGPKKLPEIGKSIGKTVKSFQQAAKEFESELKTEPEE...ENV'}, length=147)",
        )
        self.assertEqual(hsp.target.id, "gi|20259265|gb|AAM14368.1|")
        self.assertEqual(hsp.target.name, "AAM14368")
        self.assertEqual(
            hsp.target.description,
            "putative Tha4 protein [Arabidopsis thaliana] >gi|15810475|gb|AAL07125.1| putative Tha4 protein [Arabidopsis thaliana] >gi|11762160|gb|AAG40358.1| AT5g28750 [Arabidopsis thaliana] >gi|7682781|gb|AAF67362.1| Hypothetical protein T32B20.e [Arabidopsis thaliana] >gi|15241912|ref|NP_198227.1| thylakoid assembly protein, putative [Arabidopsis thaliana]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "G+ + +L +IA +  LLFG KKL  IG  +G ++K F++A  + E +     +++   +  +    A +N+E+ K  +     KE V",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|202592        64 GLGVPELAVIAGVAALLFGPKKLPEIGKSIGKTVKSFQQAAKEFESELKTEPEESVAESS
                  0 |................||.|||..||...|...|.|..|....|...............
gi|491764        16 GISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAK

gi|202592       124 QV----ATSNKEEEKKTEVSSSSKENV 147
                 60 ..----|..|.|..|........||.|  87
gi|491764        76 TIADKQADTNQEQAKTEDAKRHDKEQV 103
""",
        )
        hit = hits[123]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|75910646|ref|YP_324942.1|")
        self.assertEqual(hit.target.name, "YP_324942")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Anabaena variabilis ATCC 29413] >gi|75704371|gb|ABA24047.1| Twin-arginine translocation protein TatA/E [Anabaena variabilis ATCC 29413] >gi|46135166|ref|ZP_00162467.2| COG1826: Sec-independent protein secretion pathway components [Anabaena variabilis ATCC 29413]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=90)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 81.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 35.8094)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.412112)
        self.assertEqual(hsp.annotations["identity"], 17)
        self.assertEqual(hsp.annotations["positive"], 34)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[21, 83],
                          [33, 95]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 62))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQ...EDA'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({21: 'FGPKKLPEIGRSLGKAIRGFQEASNEFQSEFKREAEQLEQAVKTTAELEPKQIE...QDS'}, length=90)",
        )
        self.assertEqual(hsp.target.id, "gi|75910646|ref|YP_324942.1|")
        self.assertEqual(hsp.target.name, "YP_324942")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Anabaena variabilis ATCC 29413] >gi|75704371|gb|ABA24047.1| Twin-arginine translocation protein TatA/E [Anabaena variabilis ATCC 29413] >gi|46135166|ref|ZP_00162467.2| COG1826: Sec-independent protein secretion pathway components [Anabaena variabilis ATCC 29413]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "FG KKL  IG  LG +I+GF++A ++ + +  + ++  +   KT A+ +    +     +D+",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|759106        21 FGPKKLPEIGRSLGKAIRGFQEASNEFQSEFKREAEQLEQAVKTTAELEPKQIESVKSEQ
                  0 ||.|||..||..||..|.||..|...................||.|..............
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTE

gi|759106        81 DS 83
                 60 |. 62
gi|491764        93 DA 95
""",
        )
        hit = hits[124]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|39982657|gb|AAR34117.1|")
        self.assertEqual(hit.target.name, "AAR34117")
        self.assertEqual(
            hit.target.description,
            "twin-arginine translocation protein, TatA/E family [Geobacter sulfurreducens PCA] >gi|39995893|ref|NP_951844.1| twin-arginine translocation protein, TatA/E family [Geobacter sulfurreducens PCA]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=57)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 81.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 35.8094)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.412112)
        self.assertEqual(hsp.annotations["identity"], 17)
        self.assertEqual(hsp.annotations["positive"], 23)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[20, 57],
                          [34, 71]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 37))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({34: 'GTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDA'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({20: 'GPSKLPQLGQALGSSIKSFKKGMNEDEVKVINKTNEA'}, length=57)",
        )
        self.assertEqual(hsp.target.id, "gi|39982657|gb|AAR34117.1|")
        self.assertEqual(hsp.target.name, "AAR34117")
        self.assertEqual(
            hsp.target.description,
            "twin-arginine translocation protein, TatA/E family [Geobacter sulfurreducens PCA] >gi|39995893|ref|NP_951844.1| twin-arginine translocation protein, TatA/E family [Geobacter sulfurreducens PCA]",
        )
        self.assertEqual(
            hsp.annotations["midline"], "G  KL  +G  LG+SIK FKK M++DE K    + +A"
        )
        self.assertEqual(
            str(hsp),
            """\
gi|399826        20 GPSKLPQLGQALGSSIKSFKKGMNEDEVKVINKTNEA 57
                  0 |..||...|..||.|||.|||.|..||.|.......| 37
gi|491764        34 GTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDA 71
""",
        )
        hit = hits[125]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|33635687|emb|CAE22011.1|")
        self.assertEqual(hit.target.name, "CAE22011")
        self.assertEqual(
            hit.target.description,
            "mttA/Hcf106 family [Prochlorococcus marinus str. MIT 9313] >gi|33864103|ref|NP_895663.1| mttA/Hcf106 family [Prochlorococcus marinus str. MIT 9313]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=91)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 80.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 35.4242)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.538235)
        self.assertEqual(hsp.annotations["identity"], 27)
        self.assertEqual(hsp.annotations["positive"], 43)
        self.assertEqual(hsp.annotations["gaps"], 2)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 4, 61, 63, 86],
                          [16, 73, 73, 96]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 82))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({16: 'GISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQD...DAK'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({4: 'GIGLPEMAVIGAIALLVFGPKRLPEFGRTLGKTLKGFQAASKEFEREIHKTMAE...EAK'}, length=91)",
        )
        self.assertEqual(hsp.target.id, "gi|33635687|emb|CAE22011.1|")
        self.assertEqual(hsp.target.name, "CAE22011")
        self.assertEqual(
            hsp.target.description,
            "mttA/Hcf106 family [Prochlorococcus marinus str. MIT 9313] >gi|33864103|ref|NP_895663.1| mttA/Hcf106 family [Prochlorococcus marinus str. MIT 9313]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "GI + ++ +I  I +L+FG K+L   G  LG ++KGF+ A  + E +  KT  + +    A T    Q DTN      ++AK",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|336356         4 GIGLPEMAVIGAIALLVFGPKRLPEFGRTLGKTLKGFQAASKEFEREIHKTMAEPESIEQ
                  0 ||...............||.|.|...|..||...|||..|....|....||......--.
gi|491764        16 GISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADF--T

gi|336356        64 AATEESMQQDTNAIGETPQEAK 86
                 60 |.|....|.|||........|| 82
gi|491764        74 AKTIADKQADTNQEQAKTEDAK 96
""",
        )
        hit = hits[126]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|76791934|ref|ZP_00774438.1|")
        self.assertEqual(hit.target.name, "ZP_00774438")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Pseudoalteromonas atlantica T6c] >gi|76592906|gb|EAO69092.1| Twin-arginine translocation protein TatA/E [Pseudoalteromonas atlantica T6c]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=68)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 80.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 35.4242)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.538235)
        self.assertEqual(hsp.annotations["identity"], 22)
        self.assertEqual(hsp.annotations["positive"], 44)
        self.assertEqual(hsp.annotations["gaps"], 3)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 2, 50, 50, 67],
                          [17, 65, 68, 85]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 68))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({17: 'ISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDA...ADT'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({2: 'LSVWQLVLVALLFILLFGRGRIPALMGDLAQGIKSFKRGISEDESQDEQGVEKT...SPT'}, length=68)",
        )
        self.assertEqual(hsp.target.id, "gi|76791934|ref|ZP_00774438.1|")
        self.assertEqual(hsp.target.name, "ZP_00774438")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Pseudoalteromonas atlantica T6c] >gi|76592906|gb|EAO69092.1| Twin-arginine translocation protein TatA/E [Pseudoalteromonas atlantica T6c]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "+S+WQL+++A++ +LLFG  ++ ++  DL   IK FK+ +S+DE + +   Q  + T+    +K++ T",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|767919         2 LSVWQLVLVALLFILLFGRGRIPALMGDLAQGIKSFKRGISEDESQDE---QGVEKTSGE
                  0 .|.||...........||.........||...||.||...|.||....---|....|...
gi|491764        17 ISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKT

gi|767919        59 TLEKKSPT 67
                 60 ...|...| 68
gi|491764        77 IADKQADT 85
""",
        )
        hit = hits[127]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|23129516|ref|ZP_00111343.1|")
        self.assertEqual(hit.target.name, "ZP_00111343")
        self.assertEqual(
            hit.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Nostoc punctiforme PCC 73102]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=91)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 80.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 35.4242)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.538235)
        self.assertEqual(hsp.annotations["identity"], 21)
        self.assertEqual(hsp.annotations["positive"], 40)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 4, 68],
                          [16, 80]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 64))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({16: 'GISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQD...IAD'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({4: 'GIGLPEMAVIMVVALLIFGPKKLPEIGRSVGKTIRSFQEASKEFQSEFQKEAEQ...TAE'}, length=91)",
        )
        self.assertEqual(hsp.target.id, "gi|23129516|ref|ZP_00111343.1|")
        self.assertEqual(hsp.target.name, "ZP_00111343")
        self.assertEqual(
            hsp.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Nostoc punctiforme PCC 73102]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "GI + ++ +I V+ +L+FG KKL  IG  +G +I+ F++A  + + +  K ++  + T KT A+",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|231295         4 GIGLPEMAVIMVVALLIFGPKKLPEIGRSVGKTIRSFQEASKEFQSEFQKEAEQLEETVK
                  0 ||...............||.|||..||...|..|..|..|.........|.......|.|
gi|491764        16 GISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAK

gi|231295        64 TTAE 68
                 60 |.|. 64
gi|491764        76 TIAD 80
""",
        )
        hit = hits[128]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|48764199|ref|ZP_00268751.1|")
        self.assertEqual(hit.target.name, "ZP_00268751")
        self.assertEqual(
            hit.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Rhodospirillum rubrum]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=96)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 79.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 35.039)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.702957)
        self.assertEqual(hsp.annotations["identity"], 22)
        self.assertEqual(hsp.annotations["positive"], 31)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 47],
                          [14, 61]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 47))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDE'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGFSSIWHWIIVLVVVLLLFGAGKIPRLMGDVAKGVKAFKKGMADDE'}, length=96)",
        )
        self.assertEqual(hsp.target.id, "gi|48764199|ref|ZP_00268751.1|")
        self.assertEqual(hsp.target.name, "ZP_00268751")
        self.assertEqual(
            hsp.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Rhodospirillum rubrum]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG  SIW  +I+ V+V+LLFG  K+  +  D+   +K FKK M+DDE",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|487641         0 MGFSSIWHWIIVLVVVLLLFGAGKIPRLMGDVAKGVKAFKKGMADDE 47
                  0 ||..|||............||..|......|.....|.|||.|.||| 47
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDE 61
""",
        )
        hit = hits[129]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|15677995|ref|NP_273645.1|")
        self.assertEqual(hit.target.name, "NP_273645")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein NMB0601 [Neisseria meningitidis MC58] >gi|7379525|emb|CAB84088.1| putative sec-independent protein translocase component [Neisseria meningitidis Z2491] >gi|66731901|gb|AAY52136.1| conserved hypothetical protein [Neisseria meningitidis MC58] >gi|15793779|ref|NP_283601.1| sec-independent protein translocase component [Neisseria meningitidis Z2491] >gi|59800638|ref|YP_207350.1| TatA [Neisseria gonorrhoeae FA 1090] >gi|59717533|gb|AAW88938.1| sec-independent protein translocase [Neisseria gonorrhoeae FA 1090] >gi|11354072|pir||E81925 probable sec-independent protein translocase component NMA0805 [imported] - Neisseria meningitidis (strain Z2491 serogroup A) >gi|54039709|sp|P66892|TATA_NEIMB Sec-independent protein translocase protein tatA/E homolog >gi|54042238|sp|P66891|TATA_NEIMA Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=67)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 79.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 35.039)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.702957)
        self.assertEqual(hsp.annotations["identity"], 21)
        self.assertEqual(hsp.annotations["positive"], 36)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 56],
                          [14, 70]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 56))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQD'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGSFSLTHWIIVLIIVVLIFGTKKLRNVGKDLGGAVHDFKQGLNEGTDGKEAQKDD'}, length=67)",
        )
        self.assertEqual(hsp.target.id, "gi|15677995|ref|NP_273645.1|")
        self.assertEqual(hsp.target.name, "NP_273645")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein NMB0601 [Neisseria meningitidis MC58] >gi|7379525|emb|CAB84088.1| putative sec-independent protein translocase component [Neisseria meningitidis Z2491] >gi|66731901|gb|AAY52136.1| conserved hypothetical protein [Neisseria meningitidis MC58] >gi|15793779|ref|NP_283601.1| sec-independent protein translocase component [Neisseria meningitidis Z2491] >gi|59800638|ref|YP_207350.1| TatA [Neisseria gonorrhoeae FA 1090] >gi|59717533|gb|AAW88938.1| sec-independent protein translocase [Neisseria gonorrhoeae FA 1090] >gi|11354072|pir||E81925 probable sec-independent protein translocase component NMA0805 [imported] - Neisseria meningitidis (strain Z2491 serogroup A) >gi|54039709|sp|P66892|TATA_NEIMB Sec-independent protein translocase protein tatA/E homolog >gi|54042238|sp|P66891|TATA_NEIMA Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG  S+   +I+ +IVVL+FGTKKL ++G DLG ++  FK+ +++    ++    D",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|156779         0 MGSFSLTHWIIVLIIVVLIFGTKKLRNVGKDLGGAVHDFKQGLNEGTDGKEAQKDD 56
                  0 ||..|..............||||||...|.|||.....||...............| 56
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQD 70
""",
        )
        hit = hits[130]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|50917153|ref|XP_468973.1|")
        self.assertEqual(hit.target.name, "XP_468973")
        self.assertEqual(
            hit.target.description,
            "putative sec-independent protein transporter [Oryza sativa (japonica cultivar-group)] >gi|41469490|gb|AAS07275.1| putative sec-independent protein transporter [Oryza sativa (japonica cultivar-group)]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=170)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 79.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 35.039)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.702957)
        self.assertEqual(hsp.annotations["identity"], 20)
        self.assertEqual(hsp.annotations["positive"], 36)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 82, 139],
                          [ 13,  70]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 57))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({13: 'CMGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQD'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({82: 'CLFGLGVPELVVIAGVAALVFGPKQLPEIGRSIGKTVKSFQQAAKEFETELKKESDD'}, length=170)",
        )
        self.assertEqual(hsp.target.id, "gi|50917153|ref|XP_468973.1|")
        self.assertEqual(hsp.target.name, "XP_468973")
        self.assertEqual(
            hsp.target.description,
            "putative sec-independent protein transporter [Oryza sativa (japonica cultivar-group)] >gi|41469490|gb|AAS07275.1| putative sec-independent protein transporter [Oryza sativa (japonica cultivar-group)]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "C+ G+ + +L++IA +  L+FG K+L  IG  +G ++K F++A  + E +  K S D",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|509171        82 CLFGLGVPELVVIAGVAALVFGPKQLPEIGRSIGKTVKSFQQAAKEFETELKKESDD 139
                  0 |..|................||.|.|..||...|...|.|..|....|....|.|.|  57
gi|491764        13 CMGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQD  70
""",
        )
        hit = hits[131]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|16329622|ref|NP_440350.1|")
        self.assertEqual(hit.target.name, "NP_440350")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein slr1046 [Synechocystis sp. PCC 6803] >gi|1652105|dbj|BAA17030.1| slr1046 [Synechocystis sp. PCC 6803] >gi|7470329|pir||S74990 hypothetical protein slr1046 - Synechocystis sp. (strain PCC 6803)",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=126)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 79.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 35.039)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.702957)
        self.assertEqual(hsp.annotations["identity"], 15)
        self.assertEqual(hsp.annotations["positive"], 34)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 62, 123],
                          [ 33,  94]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 61))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQ...TED'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({62: 'FGPKKLPEVGRSLGKALRGFQEASKEFETELKREAQNLEKSVQIKAELEESKTP...SSE'}, length=126)",
        )
        self.assertEqual(hsp.target.id, "gi|16329622|ref|NP_440350.1|")
        self.assertEqual(hsp.target.name, "NP_440350")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein slr1046 [Synechocystis sp. PCC 6803] >gi|1652105|dbj|BAA17030.1| slr1046 [Synechocystis sp. PCC 6803] >gi|7470329|pir||S74990 hypothetical protein slr1046 - Synechocystis sp. (strain PCC 6803)",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "FG KKL  +G  LG +++GF++A  + E +  + +Q+ + + +  A+ +     E + + +",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|163296        62 FGPKKLPEVGRSLGKALRGFQEASKEFETELKREAQNLEKSVQIKAELEESKTPESSSSS
                  0 ||.|||...|..||....||..|....|.......|.........|........|.....
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTE

gi|163296       122 E 123
                 60 .  61
gi|491764        93 D  94
""",
        )
        hit = hits[132]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|71083667|ref|YP_266387.1|")
        self.assertEqual(hit.target.name, "YP_266387")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Candidatus Pelagibacter ubique HTCC1062] >gi|71062780|gb|AAZ21783.1| Twin-arginine translocation protein TatA/E [Candidatus Pelagibacter ubique HTCC1062]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=66)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 79.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 35.039)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.702957)
        self.assertEqual(hsp.annotations["identity"], 24)
        self.assertEqual(hsp.annotations["positive"], 35)
        self.assertEqual(hsp.annotations["gaps"], 3)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 2, 43, 46, 58],
                          [17, 58, 58, 70]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 56))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({17: 'ISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQD'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({2: 'IGIWQIAIVVILVVLLFGRGKISSLMGDVAKGIKSFKKGMATDITDEPEPKNVSEN'}, length=66)",
        )
        self.assertEqual(hsp.target.id, "gi|71083667|ref|YP_266387.1|")
        self.assertEqual(hsp.target.name, "YP_266387")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Candidatus Pelagibacter ubique HTCC1062] >gi|71062780|gb|AAZ21783.1| Twin-arginine translocation protein TatA/E [Candidatus Pelagibacter ubique HTCC1062]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "I IWQ+ I+ ++VVLLFG  K+ S+  D+   IK FKK M+    DEP+    S++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|710836         2 IGIWQIAIVVILVVLLFGRGKISSLMGDVAKGIKSFKKGMATDITDEPEPKNVSEN 58
                  0 |.|||...........||..|..|...|....||.|||.|.---.|||.....|.. 56
gi|491764        17 ISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMS---DDEPKQDKTSQD 70
""",
        )
        hit = hits[133]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|17130190|dbj|BAB72802.1|")
        self.assertEqual(hit.target.name, "BAB72802")
        self.assertEqual(
            hit.target.description,
            "asl0845 [Nostoc sp. PCC 7120] >gi|25530529|pir||AC1912 hypothetical protein asl0845 [imported] - Nostoc sp. (strain PCC 7120) >gi|17228340|ref|NP_484888.1| hypothetical protein asl0845 [Nostoc sp. PCC 7120]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=90)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 79.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 35.039)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.702957)
        self.assertEqual(hsp.annotations["identity"], 16)
        self.assertEqual(hsp.annotations["positive"], 34)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[21, 83],
                          [33, 95]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 62))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQ...EDA'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({21: 'FGPKKLPEIGRSMGKAIRGFQEASNEFQSEFKREAEQLEQAVKTTAELEPKQIE...QDS'}, length=90)",
        )
        self.assertEqual(hsp.target.id, "gi|17130190|dbj|BAB72802.1|")
        self.assertEqual(hsp.target.name, "BAB72802")
        self.assertEqual(
            hsp.target.description,
            "asl0845 [Nostoc sp. PCC 7120] >gi|25530529|pir||AC1912 hypothetical protein asl0845 [imported] - Nostoc sp. (strain PCC 7120) >gi|17228340|ref|NP_484888.1| hypothetical protein asl0845 [Nostoc sp. PCC 7120]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "FG KKL  IG  +G +I+GF++A ++ + +  + ++  +   KT A+ +    +     +D+",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|171301        21 FGPKKLPEIGRSMGKAIRGFQEASNEFQSEFKREAEQLEQAVKTTAELEPKQIESVKSEQ
                  0 ||.|||..||...|..|.||..|...................||.|..............
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTE

gi|171301        81 DS 83
                 60 |. 62
gi|491764        93 DA 95
""",
        )
        hit = hits[134]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|68246031|gb|EAN28138.1|")
        self.assertEqual(hit.target.name, "EAN28138")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Magnetococcus sp. MC-1] >gi|69259242|ref|ZP_00607441.1| Twin-arginine translocation protein TatA/E [Magnetococcus sp. MC-1]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=69)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 78.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 34.6538)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.91809)
        self.assertEqual(hsp.annotations["identity"], 17)
        self.assertEqual(hsp.annotations["positive"], 24)
        self.assertEqual(hsp.annotations["gaps"], 2)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[19, 44, 46, 59],
                          [33, 58, 58, 71]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 40))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDA'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({19: 'FGAGKLPKVMGDMGRGVKSFKKAMNAEDDAPAEPEVSKPA'}, length=69)",
        )
        self.assertEqual(hsp.target.id, "gi|68246031|gb|EAN28138.1|")
        self.assertEqual(hsp.target.name, "EAN28138")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Magnetococcus sp. MC-1] >gi|69259242|ref|ZP_00607441.1| Twin-arginine translocation protein TatA/E [Magnetococcus sp. MC-1]",
        )
        self.assertEqual(
            hsp.annotations["midline"], "FG  KL  +  D+G  +K FKKAM+  DD P + + S+ A"
        )
        self.assertEqual(
            str(hsp),
            """\
gi|682460        19 FGAGKLPKVMGDMGRGVKSFKKAMNAEDDAPAEPEVSKPA 59
                  0 ||..||.....|.|...|.|||||.--||.|.....|..| 40
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMS--DDEPKQDKTSQDA 71
""",
        )
        hit = hits[135]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|15604583|ref|NP_221101.1|")
        self.assertEqual(hit.target.name, "NP_221101")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein RP749 [Rickettsia prowazekii str. Madrid E] >gi|3861278|emb|CAA15177.1| unknown [Rickettsia prowazekii] >gi|7467844|pir||A71635 hypothetical protein RP749 - Rickettsia prowazekii >gi|9979048|sp|Q9ZCJ1|TATA_RICPR Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=54)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 78.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 34.6538)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.91809)
        self.assertEqual(hsp.annotations["identity"], 16)
        self.assertEqual(hsp.annotations["positive"], 20)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[18, 53],
                          [33, 68]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 35))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({18: 'FGAGKLPQVMSDLAKGLKAFKEGMKDDGNDNDKTN'}, length=54)",
        )
        self.assertEqual(hsp.target.id, "gi|15604583|ref|NP_221101.1|")
        self.assertEqual(hsp.target.name, "NP_221101")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein RP749 [Rickettsia prowazekii str. Madrid E] >gi|3861278|emb|CAA15177.1| unknown [Rickettsia prowazekii] >gi|7467844|pir||A71635 hypothetical protein RP749 - Rickettsia prowazekii >gi|9979048|sp|Q9ZCJ1|TATA_RICPR Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(
            hsp.annotations["midline"], "FG  KL  + SDL   +K FK+ M DD    DKT+"
        )
        self.assertEqual(
            str(hsp),
            """\
gi|156045        18 FGAGKLPQVMSDLAKGLKAFKEGMKDDGNDNDKTN 53
                  0 ||..||....|||....|.||..|.||....|||. 35
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS 68
""",
        )
        hit = hits[136]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|77685166|ref|ZP_00800574.1|")
        self.assertEqual(hit.target.name, "ZP_00800574")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Alkaliphilus metalliredigenes QYMF] >gi|77639020|gb|EAO81400.1| Twin-arginine translocation protein TatA/E [Alkaliphilus metalliredigenes QYMF]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=69)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 78.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 34.6538)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.91809)
        self.assertEqual(hsp.annotations["identity"], 19)
        self.assertEqual(hsp.annotations["positive"], 24)
        self.assertEqual(hsp.annotations["gaps"], 3)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[20, 41, 44, 62],
                          [33, 54, 54, 72]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 42))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDAD'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({20: 'FGPSKLPEIGKSLGKSIKEFKKFSKEMKDDLSLEERETKDKD'}, length=69)",
        )
        self.assertEqual(hsp.target.id, "gi|77685166|ref|ZP_00800574.1|")
        self.assertEqual(hsp.target.name, "ZP_00800574")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Alkaliphilus metalliredigenes QYMF] >gi|77639020|gb|EAO81400.1| Twin-arginine translocation protein TatA/E [Alkaliphilus metalliredigenes QYMF]",
        )
        self.assertEqual(
            hsp.annotations["midline"], "FG  KL  IG  LG SIK FK   K M DD   +++ ++D D"
        )
        self.assertEqual(
            str(hsp),
            """\
gi|776851        20 FGPSKLPEIGKSLGKSIKEFKKFSKEMKDDLSLEERETKDKD 62
                  0 ||..||..||..||.|||.||---|.|.||.........|.| 42
gi|491764        33 FGTKKLGSIGSDLGASIKGFK---KAMSDDEPKQDKTSQDAD 72
""",
        )
        hit = hits[137]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|39985226|gb|AAR36581.1|")
        self.assertEqual(hit.target.name, "AAR36581")
        self.assertEqual(
            hit.target.description,
            "twin-arginine translocation protein, TatA/E family [Geobacter sulfurreducens PCA] >gi|39998280|ref|NP_954231.1| twin-arginine translocation protein, TatA/E family [Geobacter sulfurreducens PCA]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=78)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 77.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 34.2686)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.19906)
        self.assertEqual(hsp.annotations["identity"], 26)
        self.assertEqual(hsp.annotations["positive"], 37)
        self.assertEqual(hsp.annotations["gaps"], 11)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[12, 57, 68, 77],
                          [11, 56, 56, 65]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 65))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({11: 'GTCMGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQD'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({12: 'GYCMFGFGMPELIVILIIVLVVFGAGRLPEIGGALGKSIRNFKKASEGKEEIEI...KKD'}, length=78)",
        )
        self.assertEqual(hsp.target.id, "gi|39985226|gb|AAR36581.1|")
        self.assertEqual(hsp.target.name, "AAR36581")
        self.assertEqual(
            hsp.target.description,
            "twin-arginine translocation protein, TatA/E family [Geobacter sulfurreducens PCA] >gi|39998280|ref|NP_954231.1| twin-arginine translocation protein, TatA/E family [Geobacter sulfurreducens PCA]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "G CM G  + +L++I +IV+++FG  +L  IG  LG SI+ FKKA              DEPK+D",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|399852        12 GYCMFGFGMPELIVILIIVLVVFGAGRLPEIGGALGKSIRNFKKASEGKEEIEIKPQKKD
                  0 |.||.|................||...|..||..||.||..||||-----------...|
gi|491764        11 GTCMGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKA-----------MSDD

gi|399852        72 EPKKD 77
                 60 |||.| 65
gi|491764        60 EPKQD 65
""",
        )
        hit = hits[138]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|1825636|gb|AAB42258.1|")
        self.assertEqual(hit.target.name, "AAB42258")
        self.assertEqual(
            hit.target.description,
            "Hypothetical protein ZK354.3 [Caenorhabditis elegans] >gi|17544554|ref|NP_500777.1| ZK354.3 [Caenorhabditis elegans]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=312)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 76.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 33.8834)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.56602)
        self.assertEqual(hsp.annotations["identity"], 15)
        self.assertEqual(hsp.annotations["positive"], 26)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[260, 309],
                          [ 53, 102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 49))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({53: 'KKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({260: 'KKEEKKEEKKDDKKEDDKEKSATKSEDKKSDEKKTEEKKSDEKKNEKSE'}, length=312)",
        )
        self.assertEqual(hsp.target.id, "gi|1825636|gb|AAB42258.1|")
        self.assertEqual(hsp.target.name, "AAB42258")
        self.assertEqual(
            hsp.target.description,
            "Hypothetical protein ZK354.3 [Caenorhabditis elegans] >gi|17544554|ref|NP_500777.1| ZK354.3 [Caenorhabditis elegans]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "KK    +E K DK   D + +A    DK++D  + + K  D K+++K +",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|182563       260 KKEEKKEEKKDDKKEDDKEKSATKSEDKKSDEKKTEEKKSDEKKNEKSE 309
                  0 ||.....|.|.||...|....|....||..|......|..|.|...|..  49
gi|491764        53 KKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQ 102
""",
        )
        hit = hits[139]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|65321915|ref|ZP_00394874.1|")
        self.assertEqual(hit.target.name, "ZP_00394874")
        self.assertEqual(
            hit.target.description,
            "COG5386: Cell surface protein [Bacillus anthracis str. A2012] >gi|47530087|ref|YP_021436.1| lpxtg-motif cell wall anchor domain protein, putative [Bacillus anthracis str. 'Ames Ancestor'] >gi|49187439|ref|YP_030691.1| LPXTG-motif cell wall anchor domain protein, putative [Bacillus anthracis str. Sterne] >gi|47505235|gb|AAT33911.1| LPXTG-motif cell wall anchor domain protein, putative [Bacillus anthracis str. 'Ames Ancestor'] >gi|49181366|gb|AAT56742.1| LPXTG-motif cell wall anchor domain protein, putative [Bacillus anthracis str. Sterne] >gi|30259275|gb|AAP28480.1| LPXTG-motif cell wall anchor domain protein, putative [Bacillus anthracis str. Ames] >gi|30264617|ref|NP_846994.1| LPXTG-motif cell wall anchor domain protein, putative [Bacillus anthracis str. Ames]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=237)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 76.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 33.8834)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.56602)
        self.assertEqual(hsp.annotations["identity"], 18)
        self.assertEqual(hsp.annotations["positive"], 37)
        self.assertEqual(hsp.annotations["gaps"], 2)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[150, 157, 159, 215],
                          [ 40,  47,  47, 103]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 65))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({40: 'SIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTED...EQV'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({150: 'AVGGDNGVAATTKNNDQAKTDTQVKEEKTKVESKETAKEVNKETKNENGKAEKT...ARI'}, length=237)",
        )
        self.assertEqual(hsp.target.id, "gi|65321915|ref|ZP_00394874.1|")
        self.assertEqual(hsp.target.name, "ZP_00394874")
        self.assertEqual(
            hsp.target.description,
            "COG5386: Cell surface protein [Bacillus anthracis str. A2012] >gi|47530087|ref|YP_021436.1| lpxtg-motif cell wall anchor domain protein, putative [Bacillus anthracis str. 'Ames Ancestor'] >gi|49187439|ref|YP_030691.1| LPXTG-motif cell wall anchor domain protein, putative [Bacillus anthracis str. Sterne] >gi|47505235|gb|AAT33911.1| LPXTG-motif cell wall anchor domain protein, putative [Bacillus anthracis str. 'Ames Ancestor'] >gi|49181366|gb|AAT56742.1| LPXTG-motif cell wall anchor domain protein, putative [Bacillus anthracis str. Sterne] >gi|30259275|gb|AAP28480.1| LPXTG-motif cell wall anchor domain protein, putative [Bacillus anthracis str. Ames] >gi|30264617|ref|NP_846994.1| LPXTG-motif cell wall anchor domain protein, putative [Bacillus anthracis str. Ames]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "++G D G  A+ K   +A +D + K++KT  ++  TAK +  +  + N +  KT++ K  D+ ++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|653219       150 AVGGDNGVAATTKNNDQAKTDTQVKEEKTKVESKETAKEVNKETKNENGKAEKTDNPKTG
                  0 ..|.|.|--|..|....|..|...|..||......|||.........|....||...|..
gi|491764        40 SIGSDLG--ASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRH

gi|653219       210 DEARI 215
                 60 |....  65
gi|491764        98 DKEQV 103
""",
        )
        hit = hits[140]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|30022625|ref|NP_834256.1|")
        self.assertEqual(hit.target.name, "NP_834256")
        self.assertEqual(
            hit.target.description,
            "Cell surface protein [Bacillus cereus ATCC 14579] >gi|29898183|gb|AAP11457.1| Cell surface protein [Bacillus cereus ATCC 14579]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=237)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 76.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 33.8834)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.56602)
        self.assertEqual(hsp.annotations["identity"], 18)
        self.assertEqual(hsp.annotations["positive"], 37)
        self.assertEqual(hsp.annotations["gaps"], 2)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[150, 157, 159, 215],
                          [ 40,  47,  47, 103]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 65))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({40: 'SIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTED...EQV'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({150: 'AVGGDNGVAATTKNNDQAKTDTQVKEEKTKVESKETAKEVNKETKNENGKAEKT...ARI'}, length=237)",
        )
        self.assertEqual(hsp.target.id, "gi|30022625|ref|NP_834256.1|")
        self.assertEqual(hsp.target.name, "NP_834256")
        self.assertEqual(
            hsp.target.description,
            "Cell surface protein [Bacillus cereus ATCC 14579] >gi|29898183|gb|AAP11457.1| Cell surface protein [Bacillus cereus ATCC 14579]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "++G D G  A+ K   +A +D + K++KT  ++  TAK +  +  + N +  KT++ K  D+ ++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|300226       150 AVGGDNGVAATTKNNDQAKTDTQVKEEKTKVESKETAKEVNKETKNENGKAEKTDNPKTG
                  0 ..|.|.|--|..|....|..|...|..||......|||.........|....||...|..
gi|491764        40 SIGSDLG--ASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRH

gi|300226       210 DEARI 215
                 60 |....  65
gi|491764        98 DKEQV 103
""",
        )
        hit = hits[141]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|55623442|ref|XP_517520.1|")
        self.assertEqual(hit.target.name, "XP_517520")
        self.assertEqual(
            hit.target.description,
            "PREDICTED: similar to Acidic leucine-rich nuclear phosphoprotein 32 family member C (Tumorigenic protein pp32r1) [Pan troglodytes]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=234)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 76.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 33.8834)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.56602)
        self.assertEqual(hsp.annotations["identity"], 14)
        self.assertEqual(hsp.annotations["positive"], 31)
        self.assertEqual(hsp.annotations["gaps"], 2)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[153, 178, 178, 205],
                          [ 47,  72,  74, 101]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 54))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({47: 'ASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKE'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({153: 'SDIKGYVEGLDDEEEGEDEEEYDEDAQVVEDEEGEDEEEEGEEEDVSGGEEE'}, length=234)",
        )
        self.assertEqual(hsp.target.id, "gi|55623442|ref|XP_517520.1|")
        self.assertEqual(hsp.target.name, "XP_517520")
        self.assertEqual(
            hsp.target.description,
            "PREDICTED: similar to Acidic leucine-rich nuclear phosphoprotein 32 family member C (Tumorigenic protein pp32r1) [Pan troglodytes]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "+ IKG+ + + D+E  +D+   D D  A+ + D++ +  +E+ + ED    ++E",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|556234       153 SDIKGYVEGLDDEEEGEDEEEYDED--AQVVEDEEGEDEEEEGEEEDVSGGEEE 205
                  0 ..|||......|.|...|....|.|--|....|.......|....||......|  54
gi|491764        47 ASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKE 101
""",
        )
        hit = hits[142]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|75762866|ref|ZP_00742681.1|")
        self.assertEqual(hit.target.name, "ZP_00742681")
        self.assertEqual(
            hit.target.description,
            "Cell surface protein [Bacillus thuringiensis serovar israelensis ATCC 35646] >gi|74489647|gb|EAO53048.1| Cell surface protein [Bacillus thuringiensis serovar israelensis ATCC 35646]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=237)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 76.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 33.8834)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.56602)
        self.assertEqual(hsp.annotations["identity"], 18)
        self.assertEqual(hsp.annotations["positive"], 37)
        self.assertEqual(hsp.annotations["gaps"], 2)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[150, 157, 159, 215],
                          [ 40,  47,  47, 103]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 65))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({40: 'SIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTED...EQV'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({150: 'AVGGDNGVAATTKNNDQAKTDTQVKEEKTKVESKETAKEVNKETKNENGKAEKT...ARI'}, length=237)",
        )
        self.assertEqual(hsp.target.id, "gi|75762866|ref|ZP_00742681.1|")
        self.assertEqual(hsp.target.name, "ZP_00742681")
        self.assertEqual(
            hsp.target.description,
            "Cell surface protein [Bacillus thuringiensis serovar israelensis ATCC 35646] >gi|74489647|gb|EAO53048.1| Cell surface protein [Bacillus thuringiensis serovar israelensis ATCC 35646]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "++G D G  A+ K   +A +D + K++KT  ++  TAK +  +  + N +  KT++ K  D+ ++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|757628       150 AVGGDNGVAATTKNNDQAKTDTQVKEEKTKVESKETAKEVNKETKNENGKAEKTDNPKTG
                  0 ..|.|.|--|..|....|..|...|..||......|||.........|....||...|..
gi|491764        40 SIGSDLG--ASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRH

gi|757628       210 DEARI 215
                 60 |....  65
gi|491764        98 DKEQV 103
""",
        )
        hit = hits[143]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|22945598|gb|AAN10511.1|")
        self.assertEqual(hit.target.name, "AAN10511")
        self.assertEqual(
            hit.target.description,
            "CG18497-PC, isoform C [Drosophila melanogaster] >gi|24580583|ref|NP_722616.1| CG18497-PC, isoform C [Drosophila melanogaster] >gi|6715140|gb|AAF26299.1| split ends [Drosophila melanogaster]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=5476)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 76.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 33.8834)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.56602)
        self.assertEqual(hsp.annotations["identity"], 17)
        self.assertEqual(hsp.annotations["positive"], 25)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[1818, 1870],
                          [  50,  102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 52))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({50: 'KGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({1818: 'KSDKNASSSDKRKNSSTSQSSKSATPRIEDDSSEADDTADKAEKNQRHEKEK'}, length=5476)",
        )
        self.assertEqual(hsp.target.id, "gi|22945598|gb|AAN10511.1|")
        self.assertEqual(hsp.target.name, "AAN10511")
        self.assertEqual(
            hsp.target.description,
            "CG18497-PC, isoform C [Drosophila melanogaster] >gi|24580583|ref|NP_722616.1| CG18497-PC, isoform C [Drosophila melanogaster] >gi|6715140|gb|AAF26299.1| split ends [Drosophila melanogaster]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "K  K A S D+ K   TSQ +      I D  ++ +    K E  +RH+KE+",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|229455      1818 KSDKNASSSDKRKNSSTSQSSKSATPRIEDDSSEADDTADKAEKNQRHEKEK 1870
                  0 |..|.|.|.|..|...|||........|.|..........|.|...||.||.   52
gi|491764        50 KGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQ  102
""",
        )
        hit = hits[144]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|10727420|gb|AAF51534.2|")
        self.assertEqual(hit.target.name, "AAF51534")
        self.assertEqual(
            hit.target.description,
            "CG18497-PB, isoform B [Drosophila melanogaster] >gi|24580581|ref|NP_524718.2| CG18497-PB, isoform B [Drosophila melanogaster]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=5533)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 76.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 33.8834)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.56602)
        self.assertEqual(hsp.annotations["identity"], 17)
        self.assertEqual(hsp.annotations["positive"], 25)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[1875, 1927],
                          [  50,  102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 52))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({50: 'KGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({1875: 'KSDKNASSSDKRKNSSTSQSSKSATPRIEDDSSEADDTADKAEKNQRHEKEK'}, length=5533)",
        )
        self.assertEqual(hsp.target.id, "gi|10727420|gb|AAF51534.2|")
        self.assertEqual(hsp.target.name, "AAF51534")
        self.assertEqual(
            hsp.target.description,
            "CG18497-PB, isoform B [Drosophila melanogaster] >gi|24580581|ref|NP_524718.2| CG18497-PB, isoform B [Drosophila melanogaster]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "K  K A S D+ K   TSQ +      I D  ++ +    K E  +RH+KE+",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|107274      1875 KSDKNASSSDKRKNSSTSQSSKSATPRIEDDSSEADDTADKAEKNQRHEKEK 1927
                  0 |..|.|.|.|..|...|||........|.|..........|.|...||.||.   52
gi|491764        50 KGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQ  102
""",
        )
        hit = hits[145]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|10727421|gb|AAF51535.2|")
        self.assertEqual(hit.target.name, "AAF51535")
        self.assertEqual(
            hit.target.description,
            "CG18497-PA, isoform A [Drosophila melanogaster] >gi|24580579|ref|NP_722615.1| CG18497-PA, isoform A [Drosophila melanogaster] >gi|46397733|sp|Q8SX83|SPEN_DROME Split ends protein",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=5560)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 76.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 33.8834)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.56602)
        self.assertEqual(hsp.annotations["identity"], 17)
        self.assertEqual(hsp.annotations["positive"], 25)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[1875, 1927],
                          [  50,  102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 52))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({50: 'KGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({1875: 'KSDKNASSSDKRKNSSTSQSSKSATPRIEDDSSEADDTADKAEKNQRHEKEK'}, length=5560)",
        )
        self.assertEqual(hsp.target.id, "gi|10727421|gb|AAF51535.2|")
        self.assertEqual(hsp.target.name, "AAF51535")
        self.assertEqual(
            hsp.target.description,
            "CG18497-PA, isoform A [Drosophila melanogaster] >gi|24580579|ref|NP_722615.1| CG18497-PA, isoform A [Drosophila melanogaster] >gi|46397733|sp|Q8SX83|SPEN_DROME Split ends protein",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "K  K A S D+ K   TSQ +      I D  ++ +    K E  +RH+KE+",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|107274      1875 KSDKNASSSDKRKNSSTSQSSKSATPRIEDDSSEADDTADKAEKNQRHEKEK 1927
                  0 |..|.|.|.|..|...|||........|.|..........|.|...||.||.   52
gi|491764        50 KGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQ  102
""",
        )
        hit = hits[146]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|71481981|ref|ZP_00661682.1|")
        self.assertEqual(hit.target.name, "ZP_00661682")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Prosthecochloris vibrioformis DSM 265] >gi|71283129|gb|EAO14957.1| Twin-arginine translocation protein TatA/E [Prosthecochloris vibrioformis DSM 265]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=69)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 76.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 33.8834)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.56602)
        self.assertEqual(hsp.annotations["identity"], 19)
        self.assertEqual(hsp.annotations["positive"], 27)
        self.assertEqual(hsp.annotations["gaps"], 1)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[19, 62, 63, 69],
                          [33, 76, 76, 82]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 50))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({19: 'FGAKKLPELAKGLGRGMKEFKKAQTEIEEEFNNAIEDTPAQKKETIKDKE'}, length=69)",
        )
        self.assertEqual(hsp.target.id, "gi|71481981|ref|ZP_00661682.1|")
        self.assertEqual(hsp.target.name, "ZP_00661682")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Prosthecochloris vibrioformis DSM 265] >gi|71283129|gb|EAO14957.1| Twin-arginine translocation protein TatA/E [Prosthecochloris vibrioformis DSM 265]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "FG KKL  +   LG  +K FKKA ++ E + +   +D     K TI DK+",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|714819        19 FGAKKLPELAKGLGRGMKEFKKAQTEIEEEFNNAIEDTPAQKKETIKDKE 69
                  0 ||.|||......||...|.||||....|........|.....|-||.||. 50
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAK-TIADKQ 82
""",
        )
        hit = hits[147]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|71150623|ref|ZP_00649545.1|")
        self.assertEqual(hit.target.name, "ZP_00649545")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Thiomicrospira denitrificans ATCC 33889] >gi|71087818|gb|EAO03597.1| Twin-arginine translocation protein TatA/E [Thiomicrospira denitrificans ATCC 33889]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=81)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 76.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 33.8834)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.56602)
        self.assertEqual(hsp.annotations["identity"], 19)
        self.assertEqual(hsp.annotations["positive"], 31)
        self.assertEqual(hsp.annotations["gaps"], 6)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[19, 54, 54, 76],
                          [33, 68, 74, 96]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 63))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQ...DAK'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({19: 'FGAKKIPELAKGIGKGIKNFKAEMKNEDDDTEVKSASTEAPKKVESAEEVASKESSK'}, length=81)",
        )
        self.assertEqual(hsp.target.id, "gi|71150623|ref|ZP_00649545.1|")
        self.assertEqual(hsp.target.name, "ZP_00649545")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Thiomicrospira denitrificans ATCC 33889] >gi|71087818|gb|EAO03597.1| Twin-arginine translocation protein TatA/E [Thiomicrospira denitrificans ATCC 33889]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "FG KK+  +   +G  IK FK  M +++   +  S      A T A K+ ++ +E A  E +K",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|711506        19 FGAKKIPELAKGIGKGIKNFKAEMKNEDDDTEVKS------ASTEAPKKVESAEEVASKE
                  0 ||.||........|..||.||..|..........|------|.|.|.|......|.|..|
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTE

gi|711506        73 SSK 76
                 60 ..| 63
gi|491764        93 DAK 96
""",
        )
        hit = hits[148]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|20151563|gb|AAM11141.1|")
        self.assertEqual(hit.target.name, "AAM11141")
        self.assertEqual(hit.target.description, "LD15253p [Drosophila melanogaster]")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=1521)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 76.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 33.8834)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.56602)
        self.assertEqual(hsp.annotations["identity"], 17)
        self.assertEqual(hsp.annotations["positive"], 25)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[1391, 1443],
                          [  50,  102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 52))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({50: 'KGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({1391: 'KSDKNASSSDKRKNSSTSQSSKSATPRIEDDSSEADDTADKAEKNQRHEKEK'}, length=1521)",
        )
        self.assertEqual(hsp.target.id, "gi|20151563|gb|AAM11141.1|")
        self.assertEqual(hsp.target.name, "AAM11141")
        self.assertEqual(hsp.target.description, "LD15253p [Drosophila melanogaster]")
        self.assertEqual(
            hsp.annotations["midline"],
            "K  K A S D+ K   TSQ +      I D  ++ +    K E  +RH+KE+",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|201515      1391 KSDKNASSSDKRKNSSTSQSSKSATPRIEDDSSEADDTADKAEKNQRHEKEK 1443
                  0 |..|.|.|.|..|...|||........|.|..........|.|...||.||.   52
gi|491764        50 KGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQ  102
""",
        )
        hit = hits[149]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|6979936|gb|AAF34661.1|")
        self.assertEqual(hit.target.name, "AAF34661")
        self.assertEqual(
            hit.target.description, "split ends long isoform [Drosophila melanogaster]"
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=5554)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 76.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 33.8834)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.56602)
        self.assertEqual(hsp.annotations["identity"], 17)
        self.assertEqual(hsp.annotations["positive"], 25)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[1869, 1921],
                          [  50,  102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 52))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({50: 'KGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({1869: 'KSDKNASSSDKRKNSSTSQSSKSATPRIEDDSSEADDTADKAEKNQRHEKEK'}, length=5554)",
        )
        self.assertEqual(hsp.target.id, "gi|6979936|gb|AAF34661.1|")
        self.assertEqual(hsp.target.name, "AAF34661")
        self.assertEqual(
            hsp.target.description, "split ends long isoform [Drosophila melanogaster]"
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "K  K A S D+ K   TSQ +      I D  ++ +    K E  +RH+KE+",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|697993      1869 KSDKNASSSDKRKNSSTSQSSKSATPRIEDDSSEADDTADKAEKNQRHEKEK 1921
                  0 |..|.|.|.|..|...|||........|.|..........|.|...||.||.   52
gi|491764        50 KGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQ  102
""",
        )
        hit = hits[150]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|6467825|gb|AAF13218.1|")
        self.assertEqual(hit.target.name, "AAF13218")
        self.assertEqual(
            hit.target.description,
            "Spen RNP motif protein long isoform [Drosophila melanogaster]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=5533)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 76.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 33.8834)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.56602)
        self.assertEqual(hsp.annotations["identity"], 17)
        self.assertEqual(hsp.annotations["positive"], 25)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[1875, 1927],
                          [  50,  102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 52))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({50: 'KGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({1875: 'KSDKNASSSDKRKNSSTSQSSKSATPRIEDDSSEADDTADKAEKNQRHEKEK'}, length=5533)",
        )
        self.assertEqual(hsp.target.id, "gi|6467825|gb|AAF13218.1|")
        self.assertEqual(hsp.target.name, "AAF13218")
        self.assertEqual(
            hsp.target.description,
            "Spen RNP motif protein long isoform [Drosophila melanogaster]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "K  K A S D+ K   TSQ +      I D  ++ +    K E  +RH+KE+",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|646782      1875 KSDKNASSSDKRKNSSTSQSSKSATPRIEDDSSEADDTADKAEKNQRHEKEK 1927
                  0 |..|.|.|.|..|...|||........|.|..........|.|...||.||.   52
gi|491764        50 KGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQ  102
""",
        )
        hit = hits[151]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|61102013|ref|ZP_00377467.1|")
        self.assertEqual(hit.target.name, "ZP_00377467")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein ELI2708 [Erythrobacter litoralis HTCC2594] >gi|60736120|gb|EAL74381.1| hypothetical protein ELI2708 [Erythrobacter litoralis HTCC2594]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=80)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 76.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 33.8834)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.56602)
        self.assertEqual(hsp.annotations["identity"], 19)
        self.assertEqual(hsp.annotations["positive"], 33)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 47],
                          [14, 61]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 47))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDE'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGQIGIWQILIIALVILVLFGRGKISDMMGDFGKGVSSFKKGLNEED'}, length=80)",
        )
        self.assertEqual(hsp.target.id, "gi|61102013|ref|ZP_00377467.1|")
        self.assertEqual(hsp.target.name, "ZP_00377467")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein ELI2708 [Erythrobacter litoralis HTCC2594] >gi|60736120|gb|EAL74381.1| hypothetical protein ELI2708 [Erythrobacter litoralis HTCC2594]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG I IWQ+LIIA+++++LFG  K+  +  D G  +  FKK +++++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|611020         0 MGQIGIWQILIIALVILVLFGRGKISDMMGDFGKGVSSFKKGLNEED 47
                  0 ||.|.|||...........||..|......|.|.....|||...... 47
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDE 61
""",
        )
        hit = hits[152]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|68056232|ref|ZP_00540361.1|")
        self.assertEqual(hit.target.name, "ZP_00540361")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Exiguobacterium sp. 255-15] >gi|68007169|gb|EAM86442.1| Twin-arginine translocation protein TatA/E [Exiguobacterium sp. 255-15]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=68)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 75.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 33.4982)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.04529)
        self.assertEqual(hsp.annotations["identity"], 17)
        self.assertEqual(hsp.annotations["positive"], 24)
        self.assertEqual(hsp.annotations["gaps"], 4)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[27, 50, 54, 67],
                          [33, 56, 56, 69]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 40))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({27: 'FGPKKLPELGRAAGQTLKEFKSATHGIMDDDKDKKDETKE'}, length=68)",
        )
        self.assertEqual(hsp.target.id, "gi|68056232|ref|ZP_00540361.1|")
        self.assertEqual(hsp.target.name, "ZP_00540361")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Exiguobacterium sp. 255-15] >gi|68007169|gb|EAM86442.1| Twin-arginine translocation protein TatA/E [Exiguobacterium sp. 255-15]",
        )
        self.assertEqual(
            hsp.annotations["midline"], "FG KKL  +G   G ++K FK A    M DD+ K+D+T +"
        )
        self.assertEqual(
            str(hsp),
            """\
gi|680562        27 FGPKKLPELGRAAGQTLKEFKSATHGIMDDDKDKKDETKE 67
                  0 ||.|||...|...|...|.||.|----|.||..|.|.|.. 40
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKA----MSDDEPKQDKTSQ 69
""",
        )
        hit = hits[153]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|68190120|gb|EAN04781.1|")
        self.assertEqual(hit.target.name, "EAN04781")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Mesorhizobium sp. BNC1] >gi|69279521|ref|ZP_00614376.1| Twin-arginine translocation protein TatA/E [Mesorhizobium sp. BNC1]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=71)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 75.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 33.4982)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.04529)
        self.assertEqual(hsp.annotations["identity"], 24)
        self.assertEqual(hsp.annotations["positive"], 34)
        self.assertEqual(hsp.annotations["gaps"], 1)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 44, 44, 52],
                          [14, 58, 59, 67]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 53))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKT'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGSFSIWHWLIVLVVVLLLFGRGKIPELMGDMAKGIKNFRKGMTDEAEEAKT'}, length=71)",
        )
        self.assertEqual(hsp.target.id, "gi|68190120|gb|EAN04781.1|")
        self.assertEqual(hsp.target.name, "EAN04781")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Mesorhizobium sp. BNC1] >gi|69279521|ref|ZP_00614376.1| Twin-arginine translocation protein TatA/E [Mesorhizobium sp. BNC1]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG  SIW  LI+ V+V+LLFG  K+  +  D+   IK F+K M+ DE ++ KT",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|681901         0 MGSFSIWHWLIVLVVVLLLFGRGKIPELMGDMAKGIKNFRKGMT-DEAEEAKT 52
                  0 ||..|||............||..|......|....||.|.|.|.-||....|| 53
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKT 67
""",
        )
        hit = hits[154]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|15605663|ref|NP_213038.1|")
        self.assertEqual(hit.target.name, "NP_213038")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein aq_064c [Aquifex aeolicus VF5] >gi|2982828|gb|AAC06451.1| hypothetical protein [Aquifex aeolicus VF5] >gi|7444822|pir||C70306 conserved hypothetical protein aq_064c - Aquifex aeolicus >gi|9978997|sp|O66478|TAT1_AQUAE Sec-independent protein translocase tatA/E-like protein 1",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=77)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 75.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 33.4982)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.04529)
        self.assertEqual(hsp.annotations["identity"], 26)
        self.assertEqual(hsp.annotations["positive"], 47)
        self.assertEqual(hsp.annotations["gaps"], 5)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 3, 58, 58, 61, 61, 76],
                          [15, 70, 74, 77, 78, 93]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 78))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({15: 'GGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQ...KTE'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({3: 'GGISMTELIIILAVILLLFGAGRLPEAGRALGEGIRNFRKALSGETEVKEVKAE...KVE'}, length=77)",
        )
        self.assertEqual(hsp.target.id, "gi|15605663|ref|NP_213038.1|")
        self.assertEqual(hsp.target.name, "NP_213038")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein aq_064c [Aquifex aeolicus VF5] >gi|2982828|gb|AAC06451.1| hypothetical protein [Aquifex aeolicus VF5] >gi|7444822|pir||C70306 conserved hypothetical protein aq_064c - Aquifex aeolicus >gi|9978997|sp|O66478|TAT1_AQUAE Sec-independent protein translocase tatA/E-like protein 1",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "GGIS+ +L+II  +++LLFG  +L   G  LG  I+ F+KA+S +   ++  ++D     KT  +++ +  +E+ K E",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|156056         3 GGISMTELIIILAVILLLFGAGRLPEAGRALGEGIRNFRKALSGETEVKEVKAED----V
                  0 ||||..............||...|...|..||..|..|.||.|...........|----.
gi|491764        15 GGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTA

gi|156056        59 KT-EERKEEKKEEKEKVE 76
                 60 ||-.........|..|.| 78
gi|491764        75 KTIADKQADTNQEQAKTE 93
""",
        )
        hit = hits[155]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|60493413|emb|CAH08199.1|")
        self.assertEqual(hit.target.name, "CAH08199")
        self.assertEqual(
            hit.target.description,
            "aerotolerance-related exported protein [Bacteroides fragilis NCTC 9343] >gi|60681979|ref|YP_212123.1| aerotolerance-related exported protein [Bacteroides fragilis NCTC 9343] >gi|4838140|gb|AAD30860.1| BatC [Bacteroides fragilis]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=238)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 75.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 33.4982)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.04529)
        self.assertEqual(hsp.annotations["identity"], 18)
        self.assertEqual(hsp.annotations["positive"], 34)
        self.assertEqual(hsp.annotations["gaps"], 6)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[110, 128, 128, 149, 154, 173],
                          [ 42,  60,  61,  82,  82, 101]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 64))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({42: 'GSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKE'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({110: 'GKDYQKAVEAYKMSLRNNPKDDETRYNLALAQKLLKDQQQNQQNQDQNQDQNKD...DKK'}, length=238)",
        )
        self.assertEqual(hsp.target.id, "gi|60493413|emb|CAH08199.1|")
        self.assertEqual(hsp.target.name, "CAH08199")
        self.assertEqual(
            hsp.target.description,
            "aerotolerance-related exported protein [Bacteroides fragilis NCTC 9343] >gi|60681979|ref|YP_212123.1| aerotolerance-related exported protein [Bacteroides fragilis NCTC 9343] >gi|4838140|gb|AAD30860.1| BatC [Bacteroides fragilis]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "G D   +++ +K ++ ++ PK D+T  +     K + D+Q      D NQ+Q K +  K+ DK+",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|604934       110 GKDYQKAVEAYKMSLRNN-PKDDETRYNLALAQKLLKDQQQNQQNQDQNQDQNKDDQQKQ
                  0 |.|........|......-||.|.|........|...|.|-----.|.||.|.|....|.
gi|491764        42 GSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQ-----ADTNQEQAKTEDAKR

gi|604934       169 QDKK 173
                 60 .||.  64
gi|491764        97 HDKE 101
""",
        )
        hit = hits[156]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|50877510|emb|CAG37350.1|")
        self.assertEqual(hit.target.name, "CAG37350")
        self.assertEqual(
            hit.target.description,
            "related to Sec-independent protein translocase protein TatA [Desulfotalea psychrophila LSv54] >gi|51246473|ref|YP_066357.1| similar to Sec-independent protein translocase protein TatA [Desulfotalea psychrophila LSv54]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=84)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 75.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 33.4982)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.04529)
        self.assertEqual(hsp.annotations["identity"], 15)
        self.assertEqual(hsp.annotations["positive"], 20)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[19, 47],
                          [33, 61]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 28))
        self.assertEqual(
            repr(hsp.query.seq), "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDE'}, length=103)"
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq), "Seq({19: 'FGGKKLPKLGSSLGQAMKSFKQGVSDVE'}, length=84)"
        )
        self.assertEqual(hsp.target.id, "gi|50877510|emb|CAG37350.1|")
        self.assertEqual(hsp.target.name, "CAG37350")
        self.assertEqual(
            hsp.target.description,
            "related to Sec-independent protein translocase protein TatA [Desulfotalea psychrophila LSv54] >gi|51246473|ref|YP_066357.1| similar to Sec-independent protein translocase protein TatA [Desulfotalea psychrophila LSv54]",
        )
        self.assertEqual(hsp.annotations["midline"], "FG KKL  +GS LG ++K FK+ +SD E")
        self.assertEqual(
            str(hsp),
            """\
gi|508775        19 FGGKKLPKLGSSLGQAMKSFKQGVSDVE 47
                  0 ||.|||...||.||...|.||...||.| 28
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDE 61
""",
        )
        hit = hits[157]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|42739647|gb|AAS43573.1|")
        self.assertEqual(hit.target.name, "AAS43573")
        self.assertEqual(
            hit.target.description,
            "conserved domain protein [Bacillus cereus ATCC 10987] >gi|42783718|ref|NP_980965.1| hypothetical protein BCE4672 [Bacillus cereus ATCC 10987]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=236)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 75.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 33.4982)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.04529)
        self.assertEqual(hsp.annotations["identity"], 18)
        self.assertEqual(hsp.annotations["positive"], 37)
        self.assertEqual(hsp.annotations["gaps"], 2)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[150, 157, 159, 215],
                          [ 40,  47,  47, 103]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 65))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({40: 'SIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTED...EQV'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({150: 'AVGGDNGVTAATKNNDQAKTDTQVKEEKTKVESKETAKEASKETKNENGKAEKT...ARI'}, length=236)",
        )
        self.assertEqual(hsp.target.id, "gi|42739647|gb|AAS43573.1|")
        self.assertEqual(hsp.target.name, "AAS43573")
        self.assertEqual(
            hsp.target.description,
            "conserved domain protein [Bacillus cereus ATCC 10987] >gi|42783718|ref|NP_980965.1| hypothetical protein BCE4672 [Bacillus cereus ATCC 10987]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "++G D G  A+ K   +A +D + K++KT  ++  TAK  + +  + N +  KT++ K  D+ ++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|427396       150 AVGGDNGVTAATKNNDQAKTDTQVKEEKTKVESKETAKEASKETKNENGKAEKTDNPKTG
                  0 ..|.|.|--|..|....|..|...|..||......|||.........|....||...|..
gi|491764        40 SIGSDLG--ASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRH

gi|427396       210 DEARI 215
                 60 |....  65
gi|491764        98 DKEQV 103
""",
        )
        hit = hits[158]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|53713708|ref|YP_099700.1|")
        self.assertEqual(hit.target.name, "YP_099700")
        self.assertEqual(
            hit.target.description,
            "conserved hypothetical protein BatC [Bacteroides fragilis YCH46] >gi|52216573|dbj|BAD49166.1| conserved hypothetical protein BatC [Bacteroides fragilis YCH46]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=238)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 75.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 33.4982)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.04529)
        self.assertEqual(hsp.annotations["identity"], 18)
        self.assertEqual(hsp.annotations["positive"], 34)
        self.assertEqual(hsp.annotations["gaps"], 6)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[110, 128, 128, 149, 154, 173],
                          [ 42,  60,  61,  82,  82, 101]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 64))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({42: 'GSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKE'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({110: 'GKDYQKAVEAYKMSLRNNPKDDETRYNLALAQKLLKDQQQNQQNQDQNQDQNKD...DKK'}, length=238)",
        )
        self.assertEqual(hsp.target.id, "gi|53713708|ref|YP_099700.1|")
        self.assertEqual(hsp.target.name, "YP_099700")
        self.assertEqual(
            hsp.target.description,
            "conserved hypothetical protein BatC [Bacteroides fragilis YCH46] >gi|52216573|dbj|BAD49166.1| conserved hypothetical protein BatC [Bacteroides fragilis YCH46]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "G D   +++ +K ++ ++ PK D+T  +     K + D+Q      D NQ+Q K +  K+ DK+",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|537137       110 GKDYQKAVEAYKMSLRNN-PKDDETRYNLALAQKLLKDQQQNQQNQDQNQDQNKDDQQKQ
                  0 |.|........|......-||.|.|........|...|.|-----.|.||.|.|....|.
gi|491764        42 GSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQ-----ADTNQEQAKTEDAKR

gi|537137       169 QDKK 173
                 60 .||.  64
gi|491764        97 HDKE 101
""",
        )
        hit = hits[159]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|33860901|ref|NP_892462.1|")
        self.assertEqual(hit.target.name, "NP_892462")
        self.assertEqual(
            hit.target.description,
            "mttA/Hcf106 family [Prochlorococcus marinus subsp. pastoris str. CCMP1986] >gi|33633843|emb|CAE18802.1| mttA/Hcf106 family [Prochlorococcus marinus subsp. pastoris str. CCMP1986]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=96)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 74.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 33.113)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.67123)
        self.assertEqual(hsp.annotations["identity"], 17)
        self.assertEqual(hsp.annotations["positive"], 34)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[21, 87],
                          [33, 99]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 66))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQ...RHD'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({21: 'FGPKKLPELGKQLGKTLKSLKKASNEFQNEIDQVINEPESENLPKSPQKKFTND...EKD'}, length=96)",
        )
        self.assertEqual(hsp.target.id, "gi|33860901|ref|NP_892462.1|")
        self.assertEqual(hsp.target.name, "NP_892462")
        self.assertEqual(
            hsp.target.description,
            "mttA/Hcf106 family [Prochlorococcus marinus subsp. pastoris str. CCMP1986] >gi|33633843|emb|CAE18802.1| mttA/Hcf106 family [Prochlorococcus marinus subsp. pastoris str. CCMP1986]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "FG KKL  +G  LG ++K  KKA ++ + + D+   + +      + ++  TN     T++++  D",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|338609        21 FGPKKLPELGKQLGKTLKSLKKASNEFQNEIDQVINEPESENLPKSPQKKFTNDLDQVTK
                  0 ||.|||...|..||...|..|||........|...................||.....|.
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTE

gi|338609        81 ESEEKD 87
                 60 .....| 66
gi|491764        93 DAKRHD 99
""",
        )
        hit = hits[160]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|48851224|ref|ZP_00305466.1|")
        self.assertEqual(hit.target.name, "ZP_00305466")
        self.assertEqual(
            hit.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Novosphingobium aromaticivorans DSM 12444]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=83)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 74.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 33.113)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.67123)
        self.assertEqual(hsp.annotations["identity"], 19)
        self.assertEqual(hsp.annotations["positive"], 35)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 49],
                          [14, 63]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 49))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPK'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGGLSLPHLIVLALVVLILFGRGRISEMMGDFGKGIKSFKQGMNDEDSK'}, length=83)",
        )
        self.assertEqual(hsp.target.id, "gi|48851224|ref|ZP_00305466.1|")
        self.assertEqual(hsp.target.name, "ZP_00305466")
        self.assertEqual(
            hsp.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Novosphingobium aromaticivorans DSM 12444]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGG+S+  L+++A++V++LFG  ++  +  D G  IK FK+ M+D++ K",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|488512         0 MGGLSLPHLIVLALVVLILFGRGRISEMMGDFGKGIKSFKQGMNDEDSK 49
                  0 |||.|..............||.........|.|..||.||..|.|...| 49
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPK 63
""",
        )
        hit = hits[161]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|67938449|ref|ZP_00530974.1|")
        self.assertEqual(hit.target.name, "ZP_00530974")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Chlorobium phaeobacteroides BS1] >gi|67915292|gb|EAM64615.1| Twin-arginine translocation protein TatA/E [Chlorobium phaeobacteroides BS1]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=69)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 74.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 33.113)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.67123)
        self.assertEqual(hsp.annotations["identity"], 16)
        self.assertEqual(hsp.annotations["positive"], 21)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[19, 56],
                          [33, 70]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 37))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQD'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({19: 'FGAKKLPELAKGLGKGMKEFKKAQSEIEEELNKAVDD'}, length=69)",
        )
        self.assertEqual(hsp.target.id, "gi|67938449|ref|ZP_00530974.1|")
        self.assertEqual(hsp.target.name, "ZP_00530974")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Chlorobium phaeobacteroides BS1] >gi|67915292|gb|EAM64615.1| Twin-arginine translocation protein TatA/E [Chlorobium phaeobacteroides BS1]",
        )
        self.assertEqual(
            hsp.annotations["midline"], "FG KKL  +   LG  +K FKKA S+ E + +K   D"
        )
        self.assertEqual(
            str(hsp),
            """\
gi|679384        19 FGAKKLPELAKGLGKGMKEFKKAQSEIEEELNKAVDD 56
                  0 ||.|||......||...|.||||.|..|....|...| 37
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQD 70
""",
        )
        hit = hits[162]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|45657833|ref|YP_001919.1|")
        self.assertEqual(hit.target.name, "YP_001919")
        self.assertEqual(
            hit.target.description,
            "sec-independent protein secretion pathway component [Leptospira interrogans serovar Copenhageni str. Fiocruz L1-130] >gi|24214626|ref|NP_712107.1| mttA/Hcf106 family protein [Leptospira interrogans serovar Lai str. 56601] >gi|45601073|gb|AAS70556.1| sec-independent protein secretion pathway component [Leptospira interrogans serovar Copenhageni str. Fiocruz L1-130] >gi|24195601|gb|AAN49125.1| mttA/Hcf106 family protein [Leptospira interrogans serovar lai str. 56601] >gi|59798911|sp|Q8F4W2|TATA_LEPIN Sec-independent protein translocase protein tatA/E homolog >gi|59798776|sp|Q72QX8|TATA_LEPIC Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=90)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 74.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 33.113)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.67123)
        self.assertEqual(hsp.annotations["identity"], 17)
        self.assertEqual(hsp.annotations["positive"], 33)
        self.assertEqual(hsp.annotations["gaps"], 9)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[31, 55, 55, 86],
                          [33, 57, 66, 97]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 64))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQ...AKR'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({31: 'FGGKRLPSLAKDLGDGIRSFRKSLTGESDDSSQQISQEQERSVPKEETKTSKSKK'}, length=90)",
        )
        self.assertEqual(hsp.target.id, "gi|45657833|ref|YP_001919.1|")
        self.assertEqual(hsp.target.name, "YP_001919")
        self.assertEqual(
            hsp.target.description,
            "sec-independent protein secretion pathway component [Leptospira interrogans serovar Copenhageni str. Fiocruz L1-130] >gi|24214626|ref|NP_712107.1| mttA/Hcf106 family protein [Leptospira interrogans serovar Lai str. 56601] >gi|45601073|gb|AAS70556.1| sec-independent protein secretion pathway component [Leptospira interrogans serovar Copenhageni str. Fiocruz L1-130] >gi|24195601|gb|AAN49125.1| mttA/Hcf106 family protein [Leptospira interrogans serovar lai str. 56601] >gi|59798911|sp|Q8F4W2|TATA_LEPIN Sec-independent protein translocase protein tatA/E homolog >gi|59798776|sp|Q72QX8|TATA_LEPIC Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "FG K+L S+  DLG  I+ F+K++         T +  D + +   +++    +E+ KT  +K+",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|456578        31 FGGKRLPSLAKDLGDGIRSFRKSL---------TGESDDSSQQISQEQERSVPKEETKTS
                  0 ||.|.|.|...|||..|..|.|..---------|....|...............|..||.
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTE

gi|456578        82 KSKK 86
                 60 ..|. 64
gi|491764        93 DAKR 97
""",
        )
        hit = hits[163]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|57238048|ref|YP_179297.1|")
        self.assertEqual(hit.target.name, "YP_179297")
        self.assertEqual(
            hit.target.description,
            "twin-arginine translocation protein, TatA/E family [Campylobacter jejuni RM1221] >gi|57166852|gb|AAW35631.1| twin-arginine translocation protein, TatA/E family [Campylobacter jejuni RM1221] >gi|6968609|emb|CAB73430.1| hypothetical protein Cj1176c [Campylobacter jejuni subsp. jejuni NCTC 11168] >gi|11346722|pir||B81323 hypothetical protein Cj1176c [imported] - Campylobacter jejuni (strain NCTC 11168) >gi|9979042|sp|Q9PNB9|TATA_CAMJE Sec-independent protein translocase protein tatA/E homolog >gi|15792500|ref|NP_282323.1| hypothetical protein Cj1176c [Campylobacter jejuni subsp. jejuni NCTC 11168]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=79)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 74.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 33.113)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.67123)
        self.assertEqual(hsp.annotations["identity"], 18)
        self.assertEqual(hsp.annotations["positive"], 30)
        self.assertEqual(hsp.annotations["gaps"], 1)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[20, 48, 48, 77],
                          [33, 61, 62, 91]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 58))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAK'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({20: 'FGAKKIPELAKGLGKGIKTFKDEMNNDDEVAKNTQKIEENKNTTNNTNADASIDETK'}, length=79)",
        )
        self.assertEqual(hsp.target.id, "gi|57238048|ref|YP_179297.1|")
        self.assertEqual(hsp.target.name, "YP_179297")
        self.assertEqual(
            hsp.target.description,
            "twin-arginine translocation protein, TatA/E family [Campylobacter jejuni RM1221] >gi|57166852|gb|AAW35631.1| twin-arginine translocation protein, TatA/E family [Campylobacter jejuni RM1221] >gi|6968609|emb|CAB73430.1| hypothetical protein Cj1176c [Campylobacter jejuni subsp. jejuni NCTC 11168] >gi|11346722|pir||B81323 hypothetical protein Cj1176c [imported] - Campylobacter jejuni (strain NCTC 11168) >gi|9979042|sp|Q9PNB9|TATA_CAMJE Sec-independent protein translocase protein tatA/E homolog >gi|15792500|ref|NP_282323.1| hypothetical protein Cj1176c [Campylobacter jejuni subsp. jejuni NCTC 11168]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "FG KK+  +   LG  IK FK  M++D+ +  K +Q  +    T  +  AD + ++ K",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|572380        20 FGAKKIPELAKGLGKGIKTFKDEMNNDD-EVAKNTQKIEENKNTTNNTNADASIDETK 77
                  0 ||.||.......||..||.||..|..|.-...|..|.......|.....||......| 58
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAK 91
""",
        )
        hit = hits[164]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|56962648|ref|YP_174374.1|")
        self.assertEqual(hit.target.name, "YP_174374")
        self.assertEqual(
            hit.target.description,
            "sec-independent protein translocase protein TatA [Bacillus clausii KSM-K16] >gi|56908886|dbj|BAD63413.1| sec-independent protein translocase protein TatA [Bacillus clausii KSM-K16]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=63)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 74.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 33.113)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.67123)
        self.assertEqual(hsp.annotations["identity"], 20)
        self.assertEqual(hsp.annotations["positive"], 38)
        self.assertEqual(hsp.annotations["gaps"], 3)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 40, 43, 57],
                          [14, 54, 54, 68]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 57))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGGLSVGSVVLIALVALLIFGPKKLPELGKAAGSTLREFKNATKGLADDDDDTKSTN'}, length=63)",
        )
        self.assertEqual(hsp.target.id, "gi|56962648|ref|YP_174374.1|")
        self.assertEqual(hsp.target.name, "YP_174374")
        self.assertEqual(
            hsp.target.description,
            "sec-independent protein translocase protein TatA [Bacillus clausii KSM-K16] >gi|56908886|dbj|BAD63413.1| sec-independent protein translocase protein TatA [Bacillus clausii KSM-K16]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGG+S+  +++IA++ +L+FG KKL  +G   G++++ FK   K ++DD+     T+",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|569626         0 MGGLSVGSVVLIALVALLIFGPKKLPELGKAAGSTLREFKNATKGLADDDDDTKSTN 57
                  0 |||.|..............||.|||...|...|.....||---|...||......|. 57
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFK---KAMSDDEPKQDKTS 68
""",
        )
        hit = hits[165]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|33239734|ref|NP_874676.1|")
        self.assertEqual(hit.target.name, "NP_874676")
        self.assertEqual(
            hit.target.description,
            "Sec-independent protein secretion pathway component TatA [Prochlorococcus marinus subsp. marinus str. CCMP1375] >gi|33237259|gb|AAP99328.1| Sec-independent protein secretion pathway component TatA [Prochlorococcus marinus subsp. marinus str. CCMP1375]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=84)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 74.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 33.113)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.67123)
        self.assertEqual(hsp.annotations["identity"], 24)
        self.assertEqual(hsp.annotations["positive"], 43)
        self.assertEqual(hsp.annotations["gaps"], 8)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 4, 54, 54, 74],
                          [16, 66, 74, 94]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 78))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({16: 'GISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQD...TED'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({4: 'GIGLPEIAVIAAIALFIFGPKKLPALGRGLGKTLRSLQKASSEFENELQKAVSA...HKE'}, length=84)",
        )
        self.assertEqual(hsp.target.id, "gi|33239734|ref|NP_874676.1|")
        self.assertEqual(hsp.target.name, "NP_874676")
        self.assertEqual(
            hsp.target.description,
            "Sec-independent protein secretion pathway component TatA [Prochlorococcus marinus subsp. marinus str. CCMP1375] >gi|33237259|gb|AAP99328.1| Sec-independent protein secretion pathway component TatA [Prochlorococcus marinus subsp. marinus str. CCMP1375]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "GI + ++ +IA I + +FG KKL ++G  LG +++  +KA S+ E +  K        A + +DK    NQ + K ++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|332397         4 GIGLPEIAVIAAIALFIFGPKKLPALGRGLGKTLRSLQKASSEFENELQK--------AV
                  0 ||...............||.|||...|..||.......||.|..|....|--------|.
gi|491764        16 GISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAK

gi|332397        56 SASDKTDQGNQHELKHKE 74
                 60 ...||....||...|... 78
gi|491764        76 TIADKQADTNQEQAKTED 94
""",
        )
        hit = hits[166]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|21674434|ref|NP_662499.1|")
        self.assertEqual(hit.target.name, "NP_662499")
        self.assertEqual(
            hit.target.description,
            "Sec-independent protein translocase protein TatE, putative [Chlorobium tepidum TLS] >gi|21647618|gb|AAM72841.1| Sec-independent protein translocase protein TatE, putative [Chlorobium tepidum TLS]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=67)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 74.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 33.113)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.67123)
        self.assertEqual(hsp.annotations["identity"], 16)
        self.assertEqual(hsp.annotations["positive"], 24)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[19, 62],
                          [33, 76]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 43))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAK'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({19: 'FGAQKLPELAKGLGKGIKEFKKAQNEIEEEFNKATDDSSSKEK'}, length=67)",
        )
        self.assertEqual(hsp.target.id, "gi|21674434|ref|NP_662499.1|")
        self.assertEqual(hsp.target.name, "NP_662499")
        self.assertEqual(
            hsp.target.description,
            "Sec-independent protein translocase protein TatE, putative [Chlorobium tepidum TLS] >gi|21647618|gb|AAM72841.1| Sec-independent protein translocase protein TatE, putative [Chlorobium tepidum TLS]",
        )
        self.assertEqual(
            hsp.annotations["midline"], "FG +KL  +   LG  IK FKKA ++ E + +K + D+    K"
        )
        self.assertEqual(
            str(hsp),
            """\
gi|216744        19 FGAQKLPELAKGLGKGIKEFKKAQNEIEEEFNKATDDSSSKEK 62
                  0 ||..||......||..||.||||....|....|...|.....| 43
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAK 76
""",
        )
        hit = hits[167]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|39968009|ref|XP_365395.1|")
        self.assertEqual(hit.target.name, "XP_365395")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein MG02097.4 [Magnaporthe grisea 70-15]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=823)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 73.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 32.7278)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.48874)
        self.assertEqual(hsp.annotations["identity"], 20)
        self.assertEqual(hsp.annotations["positive"], 32)
        self.assertEqual(hsp.annotations["gaps"], 4)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 67, 124, 128, 139],
                          [ 34,  91,  91, 102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 72))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({34: 'GTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQE...KEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({67: 'GAPLIGKAVDKIVDKVSGGDKKKKEEEEAKKKAEADAAAAAEAEAKKKADAEAA...EEE'}, length=823)",
        )
        self.assertEqual(hsp.target.id, "gi|39968009|ref|XP_365395.1|")
        self.assertEqual(hsp.target.name, "XP_365395")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein MG02097.4 [Magnaporthe grisea 70-15]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "G   +G     +   + G  K   ++E  + K   DA   A+  A K+AD    +AK     ED K+ D+E+",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|399680        67 GAPLIGKAVDKIVDKVSGGDKKKKEEEEAKKKAEADAAAAAEAEAKKKADAEAAEAKKKK
                  0 |....|...........|..|.....|....|...||...|...|.|.||.....||---
gi|491764        34 GTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAK---

gi|399680       127 EAEDKKKKDEEE 139
                 60 -.||.|..|.|.  72
gi|491764        91 -TEDAKRHDKEQ 102
""",
        )
        hit = hits[168]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|4877986|gb|AAD31523.1|")
        self.assertEqual(hit.target.name, "AAD31523")
        self.assertEqual(hit.target.description, "THA9 [Zea mays]")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=169)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 73.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 32.7278)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.48874)
        self.assertEqual(hsp.annotations["identity"], 19)
        self.assertEqual(hsp.annotations["positive"], 33)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 73, 126],
                          [ 13,  66]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 53))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({13: 'CMGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDK'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({73: 'CLFGLGVPELAVIAGVAALVFGPKQLPEIGRSLGKTVKSFQEAAKEFESELKK'}, length=169)",
        )
        self.assertEqual(hsp.target.id, "gi|4877986|gb|AAD31523.1|")
        self.assertEqual(hsp.target.name, "AAD31523")
        self.assertEqual(hsp.target.description, "THA9 [Zea mays]")
        self.assertEqual(
            hsp.annotations["midline"],
            "C+ G+ + +L +IA +  L+FG K+L  IG  LG ++K F++A  + E +  K",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|487798        73 CLFGLGVPELAVIAGVAALVFGPKQLPEIGRSLGKTVKSFQEAAKEFESELKK 126
                  0 |..|................||.|.|..||..||...|.|..|....|....|  53
gi|491764        13 CMGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDK  66
""",
        )
        hit = hits[169]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|67934419|ref|ZP_00527476.1|")
        self.assertEqual(hit.target.name, "ZP_00527476")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Solibacter usitatus Ellin6076] >gi|67858348|gb|EAM53485.1| Twin-arginine translocation protein TatA/E [Solibacter usitatus Ellin6076]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=56)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 73.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 32.7278)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.48874)
        self.assertEqual(hsp.annotations["identity"], 14)
        self.assertEqual(hsp.annotations["positive"], 20)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[19, 55],
                          [33, 69]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 36))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({19: 'FGGKKIPEVAKGLGEGIRNFKSALKGDEEKVEEKKQ'}, length=56)",
        )
        self.assertEqual(hsp.target.id, "gi|67934419|ref|ZP_00527476.1|")
        self.assertEqual(hsp.target.name, "ZP_00527476")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Solibacter usitatus Ellin6076] >gi|67858348|gb|EAM53485.1| Twin-arginine translocation protein TatA/E [Solibacter usitatus Ellin6076]",
        )
        self.assertEqual(
            hsp.annotations["midline"], "FG KK+  +   LG  I+ FK A+  DE K ++  Q"
        )
        self.assertEqual(
            str(hsp),
            """\
gi|679344        19 FGGKKIPEVAKGLGEGIRNFKSALKGDEEKVEEKKQ 55
                  0 ||.||.......||..|..||.|...||.|.....| 36
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQ 69
""",
        )
        hit = hits[170]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|42523658|ref|NP_969038.1|")
        self.assertEqual(hit.target.name, "NP_969038")
        self.assertEqual(
            hit.target.description,
            "twin-argine protein translocase component [Bdellovibrio bacteriovorus HD100] >gi|39575865|emb|CAE80031.1| twin-argine protein translocase component [Bdellovibrio bacteriovorus HD100]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=79)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 73.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 32.7278)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.48874)
        self.assertEqual(hsp.annotations["identity"], 24)
        self.assertEqual(hsp.annotations["positive"], 45)
        self.assertEqual(hsp.annotations["gaps"], 11)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0,  44,  44,  77],
                          [ 14,  58,  69, 102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 88))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...KEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGEFSLTHILLLAVIFLIFFGPSRLPQLGQSMGKAIRGFKQGLNEIDVDAKDIH...ENQ'}, length=79)",
        )
        self.assertEqual(hsp.target.id, "gi|42523658|ref|NP_969038.1|")
        self.assertEqual(hsp.target.name, "NP_969038")
        self.assertEqual(
            hsp.target.description,
            "twin-argine protein translocase component [Bdellovibrio bacteriovorus HD100] >gi|39575865|emb|CAE80031.1| twin-argine protein translocase component [Bdellovibrio bacteriovorus HD100]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG  S+  +L++AVI ++ FG  +L  +G  +G +I+GFK+ ++           + D  AK I D Q  ++Q +      ++  + Q",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|425236         0 MGEFSLTHILLLAVIFLIFFGPSRLPQLGQSMGKAIRGFKQGLN-----------EIDVD
                  0 ||..|..............||...|...|...|..|.|||....-----------..|..
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFT

gi|425236        49 AKDIHDNQQVSHQNKQSMGQTQKQGENQ  77
                 60 ||.|.|.|....|..............|  88
gi|491764        74 AKTIADKQADTNQEQAKTEDAKRHDKEQ 102
""",
        )
        hit = hits[171]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|71546080|ref|ZP_00666945.1|")
        self.assertEqual(hit.target.name, "ZP_00666945")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Syntrophobacter fumaroxidans MPOB] >gi|71488143|gb|EAO20575.1| Twin-arginine translocation protein TatA/E [Syntrophobacter fumaroxidans MPOB]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=73)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 73.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 32.7278)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.48874)
        self.assertEqual(hsp.annotations["identity"], 16)
        self.assertEqual(hsp.annotations["positive"], 20)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[17, 55],
                          [33, 71]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 38))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDA'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({17: 'FGPGKLPDIGKSLGKGIRNFKKATDEESPVPQVPSKTA'}, length=73)",
        )
        self.assertEqual(hsp.target.id, "gi|71546080|ref|ZP_00666945.1|")
        self.assertEqual(hsp.target.name, "ZP_00666945")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Syntrophobacter fumaroxidans MPOB] >gi|71488143|gb|EAO20575.1| Twin-arginine translocation protein TatA/E [Syntrophobacter fumaroxidans MPOB]",
        )
        self.assertEqual(
            hsp.annotations["midline"], "FG  KL  IG  LG  I+ FKKA  ++ P     S+ A"
        )
        self.assertEqual(
            str(hsp),
            """\
gi|715460        17 FGPGKLPDIGKSLGKGIRNFKKATDEESPVPQVPSKTA 55
                  0 ||..||..||..||..|..||||.....|.....|..| 38
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDA 71
""",
        )
        hit = hits[172]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|68002197|ref|ZP_00534828.1|")
        self.assertEqual(hit.target.name, "ZP_00534828")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Geobacter metallireducens GS-15] >gi|67993414|gb|EAM80570.1| Twin-arginine translocation protein TatA/E [Geobacter metallireducens GS-15]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=60)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 72.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 32.3426)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.55643)
        self.assertEqual(hsp.annotations["identity"], 18)
        self.assertEqual(hsp.annotations["positive"], 24)
        self.assertEqual(hsp.annotations["gaps"], 1)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[19, 42, 42, 60],
                          [33, 56, 57, 75]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 42))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTA'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({19: 'FGPAKLPQLGQSLGAGIRNFKKATLEEPEQIATKDKEEGAA'}, length=60)",
        )
        self.assertEqual(hsp.target.id, "gi|68002197|ref|ZP_00534828.1|")
        self.assertEqual(hsp.target.name, "ZP_00534828")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Geobacter metallireducens GS-15] >gi|67993414|gb|EAM80570.1| Twin-arginine translocation protein TatA/E [Geobacter metallireducens GS-15]",
        )
        self.assertEqual(
            hsp.annotations["midline"], "FG  KL  +G  LGA I+ FKKA + +EP+Q  T    +  A"
        )
        self.assertEqual(
            str(hsp),
            """\
gi|680021        19 FGPAKLPQLGQSLGAGIRNFKKA-TLEEPEQIATKDKEEGAA 60
                  0 ||..||...|..|||.|..||||-...||.|..|.......| 42
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTA 75
""",
        )
        hit = hits[173]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|67481641|ref|XP_656170.1|")
        self.assertEqual(hit.target.name, "XP_656170")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein 16.t00066 [Entamoeba histolytica HM-1:IMSS] >gi|56473354|gb|EAL50784.1| hypothetical protein 16.t00066 [Entamoeba histolytica HM-1:IMSS]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=434)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 72.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 32.3426)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.55643)
        self.assertEqual(hsp.annotations["identity"], 21)
        self.assertEqual(hsp.annotations["positive"], 32)
        self.assertEqual(hsp.annotations["gaps"], 3)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[164, 192, 195, 233],
                          [ 35,  63,  63, 101]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 69))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({35: 'TKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQ...DKE'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({164: 'TKGSESSKTDLTKSCKGQKAHTKEDEKKLHSSDSESSDSESSYSDSDSESSDSE...DRE'}, length=434)",
        )
        self.assertEqual(hsp.target.id, "gi|67481641|ref|XP_656170.1|")
        self.assertEqual(hsp.target.name, "XP_656170")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein 16.t00066 [Entamoeba histolytica HM-1:IMSS] >gi|56473354|gb|EAL50784.1| hypothetical protein 16.t00066 [Entamoeba histolytica HM-1:IMSS]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "TK   S  +DL  S KG K    +DE K    D  S D++ +      + +D+  E++ TE     D+E",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|674816       164 TKGSESSKTDLTKSCKGQKAHTKEDEKKLHSSDSESSDSESSYSDSDSESSDSEDEESST
                  0 ||...|...||..|.||.|.....||.|---.|..|.|.............|...|...|
gi|491764        35 TKKLGSIGSDLGASIKGFKKAMSDDEPK---QDKTSQDADFTAKTIADKQADTNQEQAKT

gi|674816       224 ESDSPSDRE 233
                 60 |.....|.|  69
gi|491764        92 EDAKRHDKE 101
""",
        )
        hit = hits[174]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|50935447|ref|XP_477251.1|")
        self.assertEqual(hit.target.name, "XP_477251")
        self.assertEqual(
            hit.target.description,
            "putative Calreticulin precursor [Oryza sativa (japonica cultivar-group)] >gi|51978962|ref|XP_507358.1| PREDICTED OJ1058_C08.28-2 gene product [Oryza sativa (japonica cultivar-group)] >gi|51963340|ref|XP_506239.1| PREDICTED OJ1058_C08.28-2 gene product [Oryza sativa (japonica cultivar-group)] >gi|50509013|dbj|BAD31961.1| putative Calreticulin precursor [Oryza sativa (japonica cultivar-group)] >gi|34393218|dbj|BAC82932.1| putative Calreticulin precursor [Oryza sativa (japonica cultivar-group)]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=424)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 72.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 32.3426)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.55643)
        self.assertEqual(hsp.annotations["identity"], 16)
        self.assertEqual(hsp.annotations["positive"], 27)
        self.assertEqual(hsp.annotations["gaps"], 5)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[382, 401, 401, 423],
                          [ 54,  73,  78, 100]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 46))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({54: 'KAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDK'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({382: 'KAGEDDDDLDDEDAEDEDKADEKADSDAEDGKDSDDEKHDE'}, length=424)",
        )
        self.assertEqual(hsp.target.id, "gi|50935447|ref|XP_477251.1|")
        self.assertEqual(hsp.target.name, "XP_477251")
        self.assertEqual(
            hsp.target.description,
            "putative Calreticulin precursor [Oryza sativa (japonica cultivar-group)] >gi|51978962|ref|XP_507358.1| PREDICTED OJ1058_C08.28-2 gene product [Oryza sativa (japonica cultivar-group)] >gi|51963340|ref|XP_506239.1| PREDICTED OJ1058_C08.28-2 gene product [Oryza sativa (japonica cultivar-group)] >gi|50509013|dbj|BAD31961.1| putative Calreticulin precursor [Oryza sativa (japonica cultivar-group)] >gi|34393218|dbj|BAC82932.1| putative Calreticulin precursor [Oryza sativa (japonica cultivar-group)]",
        )
        self.assertEqual(
            hsp.annotations["midline"], "KA  DD+   D+ ++D D      AD++AD++ E  K  D ++HD+"
        )
        self.assertEqual(
            str(hsp),
            """\
gi|509354       382 KAGEDDDDLDDEDAEDEDK-----ADEKADSDAEDGKDSDDEKHDE 423
                  0 ||..||....|....|.|.-----||..||...|..|..|...||.  46
gi|491764        54 KAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDK 100
""",
        )
        hit = hits[175]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|50978634|ref|NP_001003013.1|")
        self.assertEqual(hit.target.name, "NP_001003013")
        self.assertEqual(
            hit.target.description,
            "acidic (leucine-rich) nuclear phosphoprotein 32 family, member A [Canis familiaris] >gi|26984047|gb|AAN85118.1| inhibitor-1 of protein phosphatase type 2A [Canis familiaris]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=249)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 72.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 32.3426)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.55643)
        self.assertEqual(hsp.annotations["identity"], 13)
        self.assertEqual(hsp.annotations["positive"], 29)
        self.assertEqual(hsp.annotations["gaps"], 2)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[160, 182, 182, 209],
                          [ 50,  72,  74, 101]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 51))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({50: 'KGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKE'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({160: 'EGYVEGLDDDEEDEDEEEYDEDAQVVEDEEDEEEEEEGEEEDVSGEEEE'}, length=249)",
        )
        self.assertEqual(hsp.target.id, "gi|50978634|ref|NP_001003013.1|")
        self.assertEqual(hsp.target.name, "NP_001003013")
        self.assertEqual(
            hsp.target.description,
            "acidic (leucine-rich) nuclear phosphoprotein 32 family, member A [Canis familiaris] >gi|26984047|gb|AAN85118.1| inhibitor-1 of protein phosphatase type 2A [Canis familiaris]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "+G+ + + DDE  +D+   D D  A+ + D++ +  +E+ + ED    ++E",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|509786       160 EGYVEGLDDDEEDEDEEEYDED--AQVVEDEEDEEEEEEGEEEDVSGEEEE 209
                  0 .|......|||...|....|.|--|....|.......|....||......|  51
gi|491764        50 KGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKE 101
""",
        )
        hit = hits[176]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|70936814|ref|XP_739300.1|")
        self.assertEqual(hit.target.name, "XP_739300")
        self.assertEqual(
            hit.target.description,
            "40S ribosomal subunit protein S6, putative [Plasmodium chabaudi] >gi|56516193|emb|CAH81548.1| 40S ribosomal subunit protein S6, putative [Plasmodium chabaudi]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=184)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 72.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 32.3426)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.55643)
        self.assertEqual(hsp.annotations["identity"], 19)
        self.assertEqual(hsp.annotations["positive"], 33)
        self.assertEqual(hsp.annotations["gaps"], 3)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 96, 103, 103, 129, 131, 154],
                          [ 44,  51,  52,  78,  78, 101]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 59))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({44: 'DLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKE'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({96: 'DVGKKLKVKKNISKDAKKNDKEKNSKDDKSKTKKNVDKKSEKSEKSKKVENKKVNEKQ'}, length=184)",
        )
        self.assertEqual(hsp.target.id, "gi|70936814|ref|XP_739300.1|")
        self.assertEqual(hsp.target.name, "XP_739300")
        self.assertEqual(
            hsp.target.description,
            "40S ribosomal subunit protein S6, putative [Plasmodium chabaudi] >gi|56516193|emb|CAH81548.1| 40S ribosomal subunit protein S6, putative [Plasmodium chabaudi]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "D+G  +K  KK +S D  K DK     D  +KT    DK+++ +++  K E+ K ++K+",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|709368        96 DVGKKLK-VKKNISKDAKKNDKEKNSKDDKSKTKKNVDKKSEKSEKSKKVENKKVNEKQ
                  0 |.|...|-.||..|.|..|.||.....|...||.--.||.........|.|..|...|.
gi|491764        44 DLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTI--ADKQADTNQEQAKTEDAKRHDKE

gi|709368       154
                 59
gi|491764       101
""",
        )
        hit = hits[177]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|68075857|ref|XP_679848.1|")
        self.assertEqual(hit.target.name, "XP_679848")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein PB104900.00.0 [Plasmodium berghei] >gi|56500688|emb|CAH97504.1| hypothetical protein PB104900.00.0 [Plasmodium berghei]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=340)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 72.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 32.3426)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.55643)
        self.assertEqual(hsp.annotations["identity"], 23)
        self.assertEqual(hsp.annotations["positive"], 30)
        self.assertEqual(hsp.annotations["gaps"], 3)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 93, 134, 137, 159],
                          [ 37,  78,  78, 100]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 66))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({37: 'KLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAK...HDK'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({93: 'KLVSLGKAHVPKNDKTKEAMKSNNNKEDKKIDDNKEDKKTDDNKEDKKTDDNKE...DKK'}, length=340)",
        )
        self.assertEqual(hsp.target.id, "gi|68075857|ref|XP_679848.1|")
        self.assertEqual(hsp.target.name, "XP_679848")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein PB104900.00.0 [Plasmodium berghei] >gi|56500688|emb|CAH97504.1| hypothetical protein PB104900.00.0 [Plasmodium berghei]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "KL S+G          K+AM  +  K+DK   D     KT     DK+ D N+E  KT+D K   K",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|680758        93 KLVSLGKAHVPKNDKTKEAMKSNNNKEDKKIDDNKEDKKTDDNKEDKKTDDNKEDKKTDD
                  0 ||.|.|..........|.||.....|.||...|.....||.---.||..|.|.|..||.|
gi|491764        37 KLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTI---ADKQADTNQEQAKTED

gi|680758       153 NKEDKK 159
                 60 .|...|  66
gi|491764        94 AKRHDK 100
""",
        )
        hit = hits[178]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|39594005|emb|CAE70115.1|")
        self.assertEqual(hit.target.name, "CAE70115")
        self.assertEqual(
            hit.target.description,
            "Hypothetical protein CBG16567 [Caenorhabditis briggsae]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=192)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 72.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 32.3426)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.55643)
        self.assertEqual(hsp.annotations["identity"], 14)
        self.assertEqual(hsp.annotations["positive"], 24)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[140, 185],
                          [ 57, 102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 45))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({57: 'SDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({140: 'SDDSSSSDDSDDDSDSEAEEDTEKKASTEEDSGKRDDAEERTFER'}, length=192)",
        )
        self.assertEqual(hsp.target.id, "gi|39594005|emb|CAE70115.1|")
        self.assertEqual(hsp.target.name, "CAE70115")
        self.assertEqual(
            hsp.target.description,
            "Hypothetical protein CBG16567 [Caenorhabditis briggsae]",
        )
        self.assertEqual(
            hsp.annotations["midline"], "SDD    D +  D+D  A+   +K+A T ++  K +DA+    E+"
        )
        self.assertEqual(
            str(hsp),
            """\
gi|395940       140 SDDSSSSDDSDDDSDSEAEEDTEKKASTEEDSGKRDDAEERTFER 185
                  0 |||....|....|.|..|.....|.|.|.....|..||.....|.  45
gi|491764        57 SDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQ 102
""",
        )
        hit = hits[179]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|66809957|ref|XP_638702.1|")
        self.assertEqual(hit.target.name, "XP_638702")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein DDB0185950 [Dictyostelium discoideum] >gi|60467300|gb|EAL65333.1| hypothetical protein DDB0185950 [Dictyostelium discoideum]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=721)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 72.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 32.3426)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.55643)
        self.assertEqual(hsp.annotations["identity"], 15)
        self.assertEqual(hsp.annotations["positive"], 27)
        self.assertEqual(hsp.annotations["gaps"], 2)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 23,  41,  41,  70],
                          [ 53,  71,  73, 102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 49))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({53: 'KKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({23: 'EQKMEVDEPKKDATASDSTTADLENKKESSGSDDTATTDVKKDDSKQ'}, length=721)",
        )
        self.assertEqual(hsp.target.id, "gi|66809957|ref|XP_638702.1|")
        self.assertEqual(hsp.target.name, "XP_638702")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein DDB0185950 [Dictyostelium discoideum] >gi|60467300|gb|EAL65333.1| hypothetical protein DDB0185950 [Dictyostelium discoideum]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "++ M  DEPK+D T+ D+  T   + +K+  +  +   T D K+ D +Q",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|668099        23 EQKMEVDEPKKDATASDS--TTADLENKKESSGSDDTATTDVKKDDSKQ  70
                  0 ...|..||||.|.|..|.--|......|..........|.|.|..|..|  49
gi|491764        53 KKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQ 102
""",
        )
        hit = hits[180]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|68550463|ref|ZP_00589911.1|")
        self.assertEqual(hit.target.name, "ZP_00589911")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Pelodictyon phaeoclathratiforme BU-1] >gi|68242628|gb|EAN24841.1| Twin-arginine translocation protein TatA/E [Pelodictyon phaeoclathratiforme BU-1]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=69)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 72.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 32.3426)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.55643)
        self.assertEqual(hsp.annotations["identity"], 17)
        self.assertEqual(hsp.annotations["positive"], 25)
        self.assertEqual(hsp.annotations["gaps"], 10)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[19, 44, 54, 64],
                          [33, 58, 58, 68]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 45))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({19: 'FGAKKLPELARGLGKGMKEFKKAQTEIEEEFNSVVDEKPKKEKTA'}, length=69)",
        )
        self.assertEqual(hsp.target.id, "gi|68550463|ref|ZP_00589911.1|")
        self.assertEqual(hsp.target.name, "ZP_00589911")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Pelodictyon phaeoclathratiforme BU-1] >gi|68242628|gb|EAN24841.1| Twin-arginine translocation protein TatA/E [Pelodictyon phaeoclathratiforme BU-1]",
        )
        self.assertEqual(
            hsp.annotations["midline"], "FG KKL  +   LG  +K FKKA +          D++PK++KT+"
        )
        self.assertEqual(
            str(hsp),
            """\
gi|685504        19 FGAKKLPELARGLGKGMKEFKKAQTEIEEEFNSVVDEKPKKEKTA 64
                  0 ||.|||......||...|.||||..----------|..||..||. 45
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMS----------DDEPKQDKTS 68
""",
        )
        hit = hits[181]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|51473916|ref|YP_067673.1|")
        self.assertEqual(hit.target.name, "YP_067673")
        self.assertEqual(
            hit.target.description,
            "TatA/E-like Sec-independent protein translocase protein [Rickettsia typhi str. Wilmington] >gi|51460228|gb|AAU04191.1| TatA/E-like Sec-independent protein translocase protein [Rickettsia typhi str. Wilmington]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=53)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 72.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 32.3426)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.55643)
        self.assertEqual(hsp.annotations["identity"], 15)
        self.assertEqual(hsp.annotations["positive"], 18)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[18, 51],
                          [33, 66]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 33))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDK'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({18: 'FGAGKLPQVMSDLAKGLKAFKEGMKDDGNDNDK'}, length=53)",
        )
        self.assertEqual(hsp.target.id, "gi|51473916|ref|YP_067673.1|")
        self.assertEqual(hsp.target.name, "YP_067673")
        self.assertEqual(
            hsp.target.description,
            "TatA/E-like Sec-independent protein translocase protein [Rickettsia typhi str. Wilmington] >gi|51460228|gb|AAU04191.1| TatA/E-like Sec-independent protein translocase protein [Rickettsia typhi str. Wilmington]",
        )
        self.assertEqual(
            hsp.annotations["midline"], "FG  KL  + SDL   +K FK+ M DD    DK"
        )
        self.assertEqual(
            str(hsp),
            """\
gi|514739        18 FGAGKLPQVMSDLAKGLKAFKEGMKDDGNDNDK 51
                  0 ||..||....|||....|.||..|.||....|| 33
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDK 66
""",
        )
        hit = hits[182]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|61857708|ref|XP_612559.1|")
        self.assertEqual(hit.target.name, "XP_612559")
        self.assertEqual(
            hit.target.description,
            "PREDICTED: similar to acidic (leucine-rich) nuclear phosphoprotein 32 family, member A [Bos taurus]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=236)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 72.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 32.3426)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.55643)
        self.assertEqual(hsp.annotations["identity"], 13)
        self.assertEqual(hsp.annotations["positive"], 29)
        self.assertEqual(hsp.annotations["gaps"], 2)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[147, 169, 169, 196],
                          [ 50,  72,  74, 101]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 51))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({50: 'KGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKE'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({147: 'EGYVEGLDDDEEDEDEEEYDEDAQVVEDEEDEEEEEEGEEEDVSGEEEE'}, length=236)",
        )
        self.assertEqual(hsp.target.id, "gi|61857708|ref|XP_612559.1|")
        self.assertEqual(hsp.target.name, "XP_612559")
        self.assertEqual(
            hsp.target.description,
            "PREDICTED: similar to acidic (leucine-rich) nuclear phosphoprotein 32 family, member A [Bos taurus]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "+G+ + + DDE  +D+   D D  A+ + D++ +  +E+ + ED    ++E",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|618577       147 EGYVEGLDDDEEDEDEEEYDED--AQVVEDEEDEEEEEEGEEEDVSGEEEE 196
                  0 .|......|||...|....|.|--|....|.......|....||......|  51
gi|491764        50 KGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKE 101
""",
        )
        hit = hits[183]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|39982651|gb|AAR34111.1|")
        self.assertEqual(hit.target.name, "AAR34111")
        self.assertEqual(
            hit.target.description,
            "twin-arginine translocation protein, TatA/E family [Geobacter sulfurreducens PCA] >gi|39995887|ref|NP_951838.1| twin-arginine translocation protein, TatA/E family [Geobacter sulfurreducens PCA]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=59)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 72.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 32.3426)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.55643)
        self.assertEqual(hsp.annotations["identity"], 15)
        self.assertEqual(hsp.annotations["positive"], 19)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[19, 49],
                          [33, 63]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 30))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPK'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({19: 'FGPGKLPQLGQSLGASIRNFKKASLEEPEK'}, length=59)",
        )
        self.assertEqual(hsp.target.id, "gi|39982651|gb|AAR34111.1|")
        self.assertEqual(hsp.target.name, "AAR34111")
        self.assertEqual(
            hsp.target.description,
            "twin-arginine translocation protein, TatA/E family [Geobacter sulfurreducens PCA] >gi|39995887|ref|NP_951838.1| twin-arginine translocation protein, TatA/E family [Geobacter sulfurreducens PCA]",
        )
        self.assertEqual(hsp.annotations["midline"], "FG  KL  +G  LGASI+ FKKA  ++  K")
        self.assertEqual(
            str(hsp),
            """\
gi|399826        19 FGPGKLPQLGQSLGASIRNFKKASLEEPEK 49
                  0 ||..||...|..|||||..||||......| 30
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPK 63
""",
        )
        hit = hits[184]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|50877509|emb|CAG37349.1|")
        self.assertEqual(hit.target.name, "CAG37349")
        self.assertEqual(
            hit.target.description,
            "related to Sec-independent protein translocase protein TatE [Desulfotalea psychrophila LSv54] >gi|51246472|ref|YP_066356.1| similar to Sec-independent protein translocase protein TatE [Desulfotalea psychrophila LSv54]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=66)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 72.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 32.3426)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.55643)
        self.assertEqual(hsp.annotations["identity"], 17)
        self.assertEqual(hsp.annotations["positive"], 24)
        self.assertEqual(hsp.annotations["gaps"], 2)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[19, 47, 49, 60],
                          [33, 61, 61, 72]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 41))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDAD'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({19: 'FGAGKLPEIGSGIGKGIKNFKEATKSEEAAKKIDKSEEKSE'}, length=66)",
        )
        self.assertEqual(hsp.target.id, "gi|50877509|emb|CAG37349.1|")
        self.assertEqual(hsp.target.name, "CAG37349")
        self.assertEqual(
            hsp.target.description,
            "related to Sec-independent protein translocase protein TatE [Desulfotalea psychrophila LSv54] >gi|51246472|ref|YP_066356.1| similar to Sec-independent protein translocase protein TatE [Desulfotalea psychrophila LSv54]",
        )
        self.assertEqual(
            hsp.annotations["midline"], "FG  KL  IGS +G  IK FK+A   +E   K DK+ + ++"
        )
        self.assertEqual(
            str(hsp),
            """\
gi|508775        19 FGAGKLPEIGSGIGKGIKNFKEATKSEEAAKKIDKSEEKSE 60
                  0 ||..||..|||..|..||.||.|....|--.|.||...... 41
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDE--PKQDKTSQDAD 72
""",
        )
        hit = hits[185]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|52699323|ref|ZP_00340731.1|")
        self.assertEqual(hit.target.name, "ZP_00340731")
        self.assertEqual(
            hit.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Rickettsia akari str. Hartford]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=53)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 72.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 32.3426)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.55643)
        self.assertEqual(hsp.annotations["identity"], 15)
        self.assertEqual(hsp.annotations["positive"], 17)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[18, 51],
                          [33, 66]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 33))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDK'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({18: 'FGAGKLPQVMSDLAKGLKAFKNGMQDDGSDNDK'}, length=53)",
        )
        self.assertEqual(hsp.target.id, "gi|52699323|ref|ZP_00340731.1|")
        self.assertEqual(hsp.target.name, "ZP_00340731")
        self.assertEqual(
            hsp.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Rickettsia akari str. Hartford]",
        )
        self.assertEqual(
            hsp.annotations["midline"], "FG  KL  + SDL   +K FK  M DD    DK"
        )
        self.assertEqual(
            str(hsp),
            """\
gi|526993        18 FGAGKLPQVMSDLAKGLKAFKNGMQDDGSDNDK 51
                  0 ||..||....|||....|.||..|.||....|| 33
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDK 66
""",
        )
        hit = hits[186]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|62426215|ref|ZP_00381343.1|")
        self.assertEqual(hit.target.name, "ZP_00381343")
        self.assertEqual(
            hit.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Brevibacterium linens BL2]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=93)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 72.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 32.3426)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.55643)
        self.assertEqual(hsp.annotations["identity"], 21)
        self.assertEqual(hsp.annotations["positive"], 34)
        self.assertEqual(hsp.annotations["gaps"], 8)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 18,  42,  46,  53,  57,  93],
                          [ 33,  57,  57,  64,  64, 100]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 75))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQ...HDK'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({18: 'FGAPKLPALARSLGQSMKIFKSEVKDLRDDDEPKKTEPGELNREATENDTSTTA...HEK'}, length=93)",
        )
        self.assertEqual(hsp.target.id, "gi|62426215|ref|ZP_00381343.1|")
        self.assertEqual(hsp.target.name, "ZP_00381343")
        self.assertEqual(
            hsp.target.description,
            "COG1826: Sec-independent protein secretion pathway components [Brevibacterium linens BL2]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "FG  KL ++   LG S+K FK  +     DDEPK+    +   +  +    T A+      Q   + +D + H+K",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|624262        18 FGAPKLPALARSLGQSMKIFKSEVKDLRDDDEPKKTEPGELNREATENDTSTTAEAARKP
                  0 ||..||......||.|.|.||...----.|||||.----............|.|......
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAM----SDDEPKQ----DKTSQDADFTAKTIADKQADT

gi|624262        78 EQPAKEKQDPQSHEK  93
                 60 .|......|...|.|  75
gi|491764        85 NQEQAKTEDAKRHDK 100
""",
        )
        hit = hits[187]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|11131838|sp|Q9SLY8|CRTC_ORYSA")
        self.assertEqual(hit.target.name, "Q9SLY8")
        self.assertEqual(
            hit.target.description,
            "Calreticulin precursor >gi|6682833|dbj|BAA88900.1| calcium-binding protein [Oryza sativa]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=424)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 72.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 32.3426)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.55643)
        self.assertEqual(hsp.annotations["identity"], 16)
        self.assertEqual(hsp.annotations["positive"], 27)
        self.assertEqual(hsp.annotations["gaps"], 5)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[382, 401, 401, 423],
                          [ 54,  73,  78, 100]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 46))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({54: 'KAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDK'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({382: 'KAGEDDDDLDDEDAEDEDKADEKADSDAEDGKDSDDEKHDE'}, length=424)",
        )
        self.assertEqual(hsp.target.id, "gi|11131838|sp|Q9SLY8|CRTC_ORYSA")
        self.assertEqual(hsp.target.name, "Q9SLY8")
        self.assertEqual(
            hsp.target.description,
            "Calreticulin precursor >gi|6682833|dbj|BAA88900.1| calcium-binding protein [Oryza sativa]",
        )
        self.assertEqual(
            hsp.annotations["midline"], "KA  DD+   D+ ++D D      AD++AD++ E  K  D ++HD+"
        )
        self.assertEqual(
            str(hsp),
            """\
gi|111318       382 KAGEDDDDLDDEDAEDEDK-----ADEKADSDAEDGKDSDDEKHDE 423
                  0 ||..||....|....|.|.-----||..||...|..|..|...||.  46
gi|491764        54 KAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDK 100
""",
        )
        hit = hits[188]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|56543690|gb|AAV89844.1|")
        self.assertEqual(hit.target.name, "AAV89844")
        self.assertEqual(
            hit.target.description,
            "Sec-independent protein secretion pathway component [Zymomonas mobilis subsp. mobilis ZM4] >gi|56552116|ref|YP_162955.1| Sec-independent protein secretion pathway component [Zymomonas mobilis subsp. mobilis ZM4]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=87)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 72.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 32.3426)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.55643)
        self.assertEqual(hsp.annotations["identity"], 19)
        self.assertEqual(hsp.annotations["positive"], 31)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 47],
                          [14, 61]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 47))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDE'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGGMSITHWIVVAVVVMIFFGKGRFSDMMGDVAKGIKSFKKGMSEDD'}, length=87)",
        )
        self.assertEqual(hsp.target.id, "gi|56543690|gb|AAV89844.1|")
        self.assertEqual(hsp.target.name, "AAV89844")
        self.assertEqual(
            hsp.target.description,
            "Sec-independent protein secretion pathway component [Zymomonas mobilis subsp. mobilis ZM4] >gi|56552116|ref|YP_162955.1| Sec-independent protein secretion pathway component [Zymomonas mobilis subsp. mobilis ZM4]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MGG+SI   +++AV+V++ FG  +   +  D+   IK FKK MS+D+",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|565436         0 MGGMSITHWIVVAVVVMIFFGKGRFSDMMGDVAKGIKSFKKGMSEDD 47
                  0 |||.||.............||.........|....||.|||.||.|. 47
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDE 61
""",
        )
        hit = hits[189]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|67923730|ref|ZP_00517196.1|")
        self.assertEqual(hit.target.name, "ZP_00517196")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Crocosphaera watsonii WH 8501] >gi|67854438|gb|EAM49731.1| Twin-arginine translocation protein TatA/E [Crocosphaera watsonii WH 8501]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=95)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 71.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 31.9574)
        self.assertAlmostEqual(hsp.annotations["evalue"], 5.95088)
        self.assertEqual(hsp.annotations["identity"], 19)
        self.assertEqual(hsp.annotations["positive"], 35)
        self.assertEqual(hsp.annotations["gaps"], 9)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[21, 70, 79, 93],
                          [33, 82, 82, 96]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 72))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQ...DAK'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({21: 'FGPKKLPEIGRSLGKAVRGFQDASKEFENEFKREAQQLEESVTMKAELEENKMA...TSE'}, length=95)",
        )
        self.assertEqual(hsp.target.id, "gi|67923730|ref|ZP_00517196.1|")
        self.assertEqual(hsp.target.name, "ZP_00517196")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Crocosphaera watsonii WH 8501] >gi|67854438|gb|EAM49731.1| Twin-arginine translocation protein TatA/E [Crocosphaera watsonii WH 8501]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "FG KKL  IG  LG +++GF+ A  + E +  + +Q  + +    A+ +          D NQ +  T+ ++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|679237        21 FGPKKLPEIGRSLGKAVRGFQDASKEFENEFKREAQQLEESVTMKAELEENKMANPTPMD
                  0 ||.|||..||..||....||..|....|.......|.........|...---------.|
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQ---------AD

gi|679237        81 NNQTKMNTQTSE 93
                 60 .||....|.... 72
gi|491764        84 TNQEQAKTEDAK 96
""",
        )
        hit = hits[190]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|67462585|ref|XP_647954.1|")
        self.assertEqual(hit.target.name, "XP_647954")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein 707.t00002 [Entamoeba histolytica HM-1:IMSS] >gi|56463733|gb|EAL42569.1| hypothetical protein 707.t00002 [Entamoeba histolytica HM-1:IMSS]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=140)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 71.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 31.9574)
        self.assertAlmostEqual(hsp.annotations["evalue"], 5.95088)
        self.assertEqual(hsp.annotations["identity"], 16)
        self.assertEqual(hsp.annotations["positive"], 30)
        self.assertEqual(hsp.annotations["gaps"], 3)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 40,  67,  70,  95],
                          [ 50,  77,  77, 102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 55))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({50: 'KGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({40: 'KNYKKIKTQEENKEETKNEDENYKEKEEIKVLQKQIKRLEEEIKEYQAKRKEEEQ'}, length=140)",
        )
        self.assertEqual(hsp.target.id, "gi|67462585|ref|XP_647954.1|")
        self.assertEqual(hsp.target.name, "XP_647954")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein 707.t00002 [Entamoeba histolytica HM-1:IMSS] >gi|56463733|gb|EAL42569.1| hypothetical protein 707.t00002 [Entamoeba histolytica HM-1:IMSS]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "K +KK  + +E K++  ++D ++  K    +  KQ    +E+ K   AKR ++EQ",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|674625        40 KNYKKIKTQEENKEETKNEDENYKEKEEIKVLQKQIKRLEEEIKEYQAKRKEEEQ  95
                  0 |..||.....|.|......|.....|.---...||.....|..|...|||...||  55
gi|491764        50 KGFKKAMSDDEPKQDKTSQDADFTAKT---IADKQADTNQEQAKTEDAKRHDKEQ 102
""",
        )
        hit = hits[191]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|51970620|dbj|BAD44002.1|")
        self.assertEqual(hit.target.name, "BAD44002")
        self.assertEqual(
            hit.target.description, "unknown protein [Arabidopsis thaliana]"
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=784)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 71.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 31.9574)
        self.assertAlmostEqual(hsp.annotations["evalue"], 5.95088)
        self.assertEqual(hsp.annotations["identity"], 21)
        self.assertEqual(hsp.annotations["positive"], 35)
        self.assertEqual(hsp.annotations["gaps"], 6)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[166, 186, 190, 216, 216, 232],
                          [ 38,  58,  58,  84,  86, 102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 68))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({38: 'LGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKT...KEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({166: 'LDSMGSELGDSRKKRKKKRSEEKQDHEDVDELANEDLEHEESEFSDEESEEEPV...KKK'}, length=784)",
        )
        self.assertEqual(hsp.target.id, "gi|51970620|dbj|BAD44002.1|")
        self.assertEqual(hsp.target.name, "BAD44002")
        self.assertEqual(
            hsp.target.description, "unknown protein [Arabidopsis thaliana]"
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "L S+GS+LG S K  KK  S    D E   +  ++D +      +D++++  +E     D KRH K++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|519706       166 LDSMGSELGDSRKKRKKKRSEEKQDHEDVDELANEDLEHEESEFSDEESE--EEPVGKRD
                  0 |.|.||.||.|.|..||..|----|.|........|.........|....--.|.....|
gi|491764        38 LGSIGSDLGASIKGFKKAMS----DDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTED

gi|519706       224 RKRHKKKK 232
                 60 .|||.|..  68
gi|491764        94 AKRHDKEQ 102
""",
        )
        hit = hits[192]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|34581241|ref|ZP_00142721.1|")
        self.assertEqual(hit.target.name, "ZP_00142721")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein [Rickettsia sibirica 246] >gi|53732344|ref|ZP_00154116.2| COG1826: Sec-independent protein secretion pathway components [Rickettsia rickettsii] >gi|28262626|gb|EAA26130.1| unknown [Rickettsia sibirica 246]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=53)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 71.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 31.9574)
        self.assertAlmostEqual(hsp.annotations["evalue"], 5.95088)
        self.assertEqual(hsp.annotations["identity"], 15)
        self.assertEqual(hsp.annotations["positive"], 17)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[18, 51],
                          [33, 66]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 33))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDK'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({18: 'FGAGKLPQVMSDLAKGLKAFKDGMKDDGSDNDK'}, length=53)",
        )
        self.assertEqual(hsp.target.id, "gi|34581241|ref|ZP_00142721.1|")
        self.assertEqual(hsp.target.name, "ZP_00142721")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein [Rickettsia sibirica 246] >gi|53732344|ref|ZP_00154116.2| COG1826: Sec-independent protein secretion pathway components [Rickettsia rickettsii] >gi|28262626|gb|EAA26130.1| unknown [Rickettsia sibirica 246]",
        )
        self.assertEqual(
            hsp.annotations["midline"], "FG  KL  + SDL   +K FK  M DD    DK"
        )
        self.assertEqual(
            str(hsp),
            """\
gi|345812        18 FGAGKLPQVMSDLAKGLKAFKDGMKDDGSDNDK 51
                  0 ||..||....|||....|.||..|.||....|| 33
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDK 66
""",
        )
        hit = hits[193]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|4877984|gb|AAD31522.1|")
        self.assertEqual(hit.target.name, "AAD31522")
        self.assertEqual(hit.target.description, "THA4 [Zea mays]")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=170)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 71.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 31.9574)
        self.assertAlmostEqual(hsp.annotations["evalue"], 5.95088)
        self.assertEqual(hsp.annotations["identity"], 18)
        self.assertEqual(hsp.annotations["positive"], 33)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 83, 136],
                          [ 13,  66]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 53))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({13: 'CMGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDK'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({83: 'CLFGLGVPELAVIAGVAALVFGPKQLPEIGRSIGKTVKSFQQAAKEFETELKK'}, length=170)",
        )
        self.assertEqual(hsp.target.id, "gi|4877984|gb|AAD31522.1|")
        self.assertEqual(hsp.target.name, "AAD31522")
        self.assertEqual(hsp.target.description, "THA4 [Zea mays]")
        self.assertEqual(
            hsp.annotations["midline"],
            "C+ G+ + +L +IA +  L+FG K+L  IG  +G ++K F++A  + E +  K",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|487798        83 CLFGLGVPELAVIAGVAALVFGPKQLPEIGRSIGKTVKSFQQAAKEFETELKK 136
                  0 |..|................||.|.|..||...|...|.|..|....|....|  53
gi|491764        13 CMGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDK  66
""",
        )
        hit = hits[194]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|9757886|dbj|BAB08393.1|")
        self.assertEqual(hit.target.name, "BAB08393")
        self.assertEqual(
            hit.target.description,
            "unnamed protein product [Arabidopsis thaliana] >gi|15238687|ref|NP_197295.1| MA3 domain-containing protein [Arabidopsis thaliana]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=707)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 71.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 31.9574)
        self.assertAlmostEqual(hsp.annotations["evalue"], 5.95088)
        self.assertEqual(hsp.annotations["identity"], 21)
        self.assertEqual(hsp.annotations["positive"], 35)
        self.assertEqual(hsp.annotations["gaps"], 6)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[146, 166, 170, 196, 196, 212],
                          [ 38,  58,  58,  84,  86, 102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 68))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({38: 'LGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKT...KEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({146: 'LDSMGSELGDSRKKRKKKRSEEKQDHEDVDELANEDLEHEESEFSDEESEEEPV...KKK'}, length=707)",
        )
        self.assertEqual(hsp.target.id, "gi|9757886|dbj|BAB08393.1|")
        self.assertEqual(hsp.target.name, "BAB08393")
        self.assertEqual(
            hsp.target.description,
            "unnamed protein product [Arabidopsis thaliana] >gi|15238687|ref|NP_197295.1| MA3 domain-containing protein [Arabidopsis thaliana]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "L S+GS+LG S K  KK  S    D E   +  ++D +      +D++++  +E     D KRH K++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|975788       146 LDSMGSELGDSRKKRKKKRSEEKQDHEDVDELANEDLEHEESEFSDEESE--EEPVGKRD
                  0 |.|.||.||.|.|..||..|----|.|........|.........|....--.|.....|
gi|491764        38 LGSIGSDLGASIKGFKKAMS----DDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTED

gi|975788       204 RKRHKKKK 212
                 60 .|||.|..  68
gi|491764        94 AKRHDKEQ 102
""",
        )
        hit = hits[195]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|32422107|ref|XP_331497.1|")
        self.assertEqual(hit.target.name, "XP_331497")
        self.assertEqual(
            hit.target.description,
            "predicted protein [Neurospora crassa] >gi|28919664|gb|EAA29078.1| predicted protein [Neurospora crassa]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=216)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 71.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 31.9574)
        self.assertAlmostEqual(hsp.annotations["evalue"], 5.95088)
        self.assertEqual(hsp.annotations["identity"], 13)
        self.assertEqual(hsp.annotations["positive"], 26)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 41,  85],
                          [ 58, 102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 44))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({58: 'DDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({41: 'EDKEEEDKEEENKEEEDKEEEDKEEEDKEEEDKEEENKEEDKEE'}, length=216)",
        )
        self.assertEqual(hsp.target.id, "gi|32422107|ref|XP_331497.1|")
        self.assertEqual(hsp.target.name, "XP_331497")
        self.assertEqual(
            hsp.target.description,
            "predicted protein [Neurospora crassa] >gi|28919664|gb|EAA29078.1| predicted protein [Neurospora crassa]",
        )
        self.assertEqual(
            hsp.annotations["midline"], "+D+ ++DK  ++ +   K   DK+ +  +E+ K E+ K  DKE+"
        )
        self.assertEqual(
            str(hsp),
            """\
gi|324221        41 EDKEEEDKEEENKEEEDKEEEDKEEEDKEEEDKEEENKEEDKEE  85
                  0 .|....||.........|...||......|..|.|..|..|||.  44
gi|491764        58 DDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQ 102
""",
        )
        hit = hits[196]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|68552035|ref|ZP_00591428.1|")
        self.assertEqual(hit.target.name, "ZP_00591428")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Prosthecochloris aestuarii DSM 271] >gi|68241158|gb|EAN23426.1| Twin-arginine translocation protein TatA/E [Prosthecochloris aestuarii DSM 271]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=70)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 71.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 31.9574)
        self.assertAlmostEqual(hsp.annotations["evalue"], 5.95088)
        self.assertEqual(hsp.annotations["identity"], 15)
        self.assertEqual(hsp.annotations["positive"], 21)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[19, 56],
                          [33, 70]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 37))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQD'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({19: 'FGAKKLPELAKGLGKGMKEFKKAQTEIEDEFNKAMSD'}, length=70)",
        )
        self.assertEqual(hsp.target.id, "gi|68552035|ref|ZP_00591428.1|")
        self.assertEqual(hsp.target.name, "ZP_00591428")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Prosthecochloris aestuarii DSM 271] >gi|68241158|gb|EAN23426.1| Twin-arginine translocation protein TatA/E [Prosthecochloris aestuarii DSM 271]",
        )
        self.assertEqual(
            hsp.annotations["midline"], "FG KKL  +   LG  +K FKKA ++ E + +K   D"
        )
        self.assertEqual(
            str(hsp),
            """\
gi|685520        19 FGAKKLPELAKGLGKGMKEFKKAQTEIEDEFNKAMSD 56
                  0 ||.|||......||...|.||||....|....|...| 37
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQD 70
""",
        )
        hit = hits[197]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|68177649|ref|ZP_00550794.1|")
        self.assertEqual(hit.target.name, "ZP_00550794")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Desulfuromonas acetoxidans DSM 684] >gi|67982303|gb|EAM71712.1| Twin-arginine translocation protein TatA/E [Desulfuromonas acetoxidans DSM 684]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=58)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 71.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 31.9574)
        self.assertAlmostEqual(hsp.annotations["evalue"], 5.95088)
        self.assertEqual(hsp.annotations["identity"], 16)
        self.assertEqual(hsp.annotations["positive"], 22)
        self.assertEqual(hsp.annotations["gaps"], 2)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[19, 43, 45, 54],
                          [33, 57, 57, 66]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 35))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDK'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({19: 'FGAKRLPQLGSGLGKGIKNFKQGIKENDTESLEDK'}, length=58)",
        )
        self.assertEqual(hsp.target.id, "gi|68177649|ref|ZP_00550794.1|")
        self.assertEqual(hsp.target.name, "ZP_00550794")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Desulfuromonas acetoxidans DSM 684] >gi|67982303|gb|EAM71712.1| Twin-arginine translocation protein TatA/E [Desulfuromonas acetoxidans DSM 684]",
        )
        self.assertEqual(
            hsp.annotations["midline"], "FG K+L  +GS LG  IK FK+ +  +D E  +DK"
        )
        self.assertEqual(
            str(hsp),
            """\
gi|681776        19 FGAKRLPQLGSGLGKGIKNFKQGIKENDTESLEDK 54
                  0 ||.|.|...||.||..||.||...--.|.|...|| 35
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAM--SDDEPKQDK 66
""",
        )
        hit = hits[198]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|67934756|ref|ZP_00527782.1|")
        self.assertEqual(hit.target.name, "ZP_00527782")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Chlorobium phaeobacteroides DSM 266] >gi|67776249|gb|EAM35909.1| Twin-arginine translocation protein TatA/E [Chlorobium phaeobacteroides DSM 266]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=65)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 71.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 31.9574)
        self.assertAlmostEqual(hsp.annotations["evalue"], 5.95088)
        self.assertEqual(hsp.annotations["identity"], 16)
        self.assertEqual(hsp.annotations["positive"], 24)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[19, 63],
                          [33, 77]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 44))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKT'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({19: 'FGAKKLPELARGLGRGMKEFKKAQTEIEEEFNKVVEEPPAKEKT'}, length=65)",
        )
        self.assertEqual(hsp.target.id, "gi|67934756|ref|ZP_00527782.1|")
        self.assertEqual(hsp.target.name, "ZP_00527782")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Chlorobium phaeobacteroides DSM 266] >gi|67776249|gb|EAM35909.1| Twin-arginine translocation protein TatA/E [Chlorobium phaeobacteroides DSM 266]",
        )
        self.assertEqual(
            hsp.annotations["midline"], "FG KKL  +   LG  +K FKKA ++ E + +K  ++     KT"
        )
        self.assertEqual(
            str(hsp),
            """\
gi|679347        19 FGAKKLPELARGLGRGMKEFKKAQTEIEEEFNKVVEEPPAKEKT 63
                  0 ||.|||......||...|.||||....|....|.........|| 44
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKT 77
""",
        )
        hit = hits[199]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|42550455|gb|EAA73298.1|")
        self.assertEqual(hit.target.name, "EAA73298")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein FG04514.1 [Gibberella zeae PH-1] >gi|46117344|ref|XP_384690.1| hypothetical protein FG04514.1 [Gibberella zeae PH-1]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=297)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 71.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 31.9574)
        self.assertAlmostEqual(hsp.annotations["evalue"], 5.95088)
        self.assertEqual(hsp.annotations["identity"], 22)
        self.assertEqual(hsp.annotations["positive"], 29)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[187, 251],
                          [ 37, 101]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 64))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({37: 'KLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAK...DKE'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({187: 'KLKSHGPDLPASQSFAKKHPRDSEDEDDVSSEEESDTDDDDDDDDDDDKTDEAK...DNE'}, length=297)",
        )
        self.assertEqual(hsp.target.id, "gi|42550455|gb|EAA73298.1|")
        self.assertEqual(hsp.target.name, "EAA73298")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein FG04514.1 [Gibberella zeae PH-1] >gi|46117344|ref|XP_384690.1| hypothetical protein FG04514.1 [Gibberella zeae PH-1]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "KL S G DL AS    KK   D E + D +S++   T     D   D   ++AK+E     D E",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|425504       187 KLKSHGPDLPASQSFAKKHPRDSEDEDDVSSEEESDTDDDDDDDDDDDKTDEAKSEGDTA
                  0 ||.|.|.||.||....||...|.|...|..|.....|.....|...|.....||.|....
gi|491764        37 KLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKR

gi|425504       247 IDNE 251
                 60 .|.|  64
gi|491764        97 HDKE 101
""",
        )
        hit = hits[200]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|15893083|ref|NP_360797.1|")
        self.assertEqual(hit.target.name, "NP_360797")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein RC1160 [Rickettsia conorii str. Malish 7] >gi|15620287|gb|AAL03698.1| unknown [Rickettsia conorii str. Malish 7] >gi|25504411|pir||H97844 hypothetical protein RC1160 [imported] - Rickettsia conorii  (strain Malish 7) >gi|24212502|sp|Q92GG3|TATA_RICCN Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=53)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 71.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 31.9574)
        self.assertAlmostEqual(hsp.annotations["evalue"], 5.95088)
        self.assertEqual(hsp.annotations["identity"], 15)
        self.assertEqual(hsp.annotations["positive"], 17)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[18, 51],
                          [33, 66]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 33))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDK'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({18: 'FGAGKLPQVMSDLAKGLKAFKDGMKDDGSDNDK'}, length=53)",
        )
        self.assertEqual(hsp.target.id, "gi|15893083|ref|NP_360797.1|")
        self.assertEqual(hsp.target.name, "NP_360797")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein RC1160 [Rickettsia conorii str. Malish 7] >gi|15620287|gb|AAL03698.1| unknown [Rickettsia conorii str. Malish 7] >gi|25504411|pir||H97844 hypothetical protein RC1160 [imported] - Rickettsia conorii  (strain Malish 7) >gi|24212502|sp|Q92GG3|TATA_RICCN Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(
            hsp.annotations["midline"], "FG  KL  + SDL   +K FK  M DD    DK"
        )
        self.assertEqual(
            str(hsp),
            """\
gi|158930        18 FGAGKLPQVMSDLAKGLKAFKDGMKDDGSDNDK 51
                  0 ||..||....|||....|.||..|.||....|| 33
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDK 66
""",
        )
        hit = hits[201]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|57233621|ref|YP_182297.1|")
        self.assertEqual(hit.target.name, "YP_182297")
        self.assertEqual(
            hit.target.description,
            "twin-arginine translocation protein, TatA/E family [Dehalococcoides ethenogenes 195] >gi|57224069|gb|AAW39126.1| twin-arginine translocation protein, TatA/E family [Dehalococcoides ethenogenes 195]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=65)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 71.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 31.9574)
        self.assertAlmostEqual(hsp.annotations["evalue"], 5.95088)
        self.assertEqual(hsp.annotations["identity"], 13)
        self.assertEqual(hsp.annotations["positive"], 23)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[19, 56],
                          [33, 70]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 37))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQD'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({19: 'FGIGKLPQVGDAIGKGIRNFRKASSGEEEKEEVETKE'}, length=65)",
        )
        self.assertEqual(hsp.target.id, "gi|57233621|ref|YP_182297.1|")
        self.assertEqual(hsp.target.name, "YP_182297")
        self.assertEqual(
            hsp.target.description,
            "twin-arginine translocation protein, TatA/E family [Dehalococcoides ethenogenes 195] >gi|57224069|gb|AAW39126.1| twin-arginine translocation protein, TatA/E family [Dehalococcoides ethenogenes 195]",
        )
        self.assertEqual(
            hsp.annotations["midline"], "FG  KL  +G  +G  I+ F+KA S +E K++  +++"
        )
        self.assertEqual(
            str(hsp),
            """\
gi|572336        19 FGIGKLPQVGDAIGKGIRNFRKASSGEEEKEEVETKE 56
                  0 ||..||...|...|..|..|.||.|..|.|....... 37
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQD 70
""",
        )
        hit = hits[202]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|75908036|ref|YP_322332.1|")
        self.assertEqual(hit.target.name, "YP_322332")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Anabaena variabilis ATCC 29413] >gi|75701761|gb|ABA21437.1| Twin-arginine translocation protein TatA/E [Anabaena variabilis ATCC 29413] >gi|53764126|ref|ZP_00159794.2| COG1826: Sec-independent protein secretion pathway components [Anabaena variabilis ATCC 29413]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=56)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 71.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 31.9574)
        self.assertAlmostEqual(hsp.annotations["evalue"], 5.95088)
        self.assertEqual(hsp.annotations["identity"], 14)
        self.assertEqual(hsp.annotations["positive"], 24)
        self.assertEqual(hsp.annotations["gaps"], 4)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[19, 43, 47, 56],
                          [33, 57, 57, 66]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 37))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDK'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({19: 'FGPKKIPELGSALGKTLRGFKEELKTPSEDTNPEEEK'}, length=56)",
        )
        self.assertEqual(hsp.target.id, "gi|75908036|ref|YP_322332.1|")
        self.assertEqual(hsp.target.name, "YP_322332")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Anabaena variabilis ATCC 29413] >gi|75701761|gb|ABA21437.1| Twin-arginine translocation protein TatA/E [Anabaena variabilis ATCC 29413] >gi|53764126|ref|ZP_00159794.2| COG1826: Sec-independent protein secretion pathway components [Anabaena variabilis ATCC 29413]",
        )
        self.assertEqual(
            hsp.annotations["midline"], "FG KK+  +GS LG +++GFK+ +     D  P+++K"
        )
        self.assertEqual(
            str(hsp),
            """\
gi|759080        19 FGPKKIPELGSALGKTLRGFKEELKTPSEDTNPEEEK 56
                  0 ||.||....||.||....|||...----.|..|...| 37
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAM----SDDEPKQDK 66
""",
        )
        hit = hits[203]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|72383453|ref|YP_292808.1|")
        self.assertEqual(hit.target.name, "YP_292808")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Prochlorococcus marinus str. NATL2A] >gi|72003303|gb|AAZ59105.1| Twin-arginine translocation protein TatA/E [Prochlorococcus marinus str. NATL2A]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=71)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 71.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 31.9574)
        self.assertAlmostEqual(hsp.annotations["evalue"], 5.95088)
        self.assertEqual(hsp.annotations["identity"], 17)
        self.assertEqual(hsp.annotations["positive"], 39)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 4, 70],
                          [16, 82]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 66))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({16: 'GISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQD...DKQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({4: 'GVGLPEIAVIAGLALVIFGPKRLPELGRTIGKTLKGLQSASTEFEREIKNAMTE...DKE'}, length=71)",
        )
        self.assertEqual(hsp.target.id, "gi|72383453|ref|YP_292808.1|")
        self.assertEqual(hsp.target.name, "YP_292808")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Prochlorococcus marinus str. NATL2A] >gi|72003303|gb|AAZ59105.1| Twin-arginine translocation protein TatA/E [Prochlorococcus marinus str. NATL2A]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "G+ + ++ +IA + +++FG K+L  +G  +G ++KG + A ++ E +      + +   + IADK+",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|723834         4 GVGLPEIAVIAGLALVIFGPKRLPELGRTIGKTLKGLQSASTEFEREIKNAMTEEESANE
                  0 |................||.|.|...|...|...||...|....|...............
gi|491764        16 GISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAK

gi|723834        64 NIADKE 70
                 60 .||||. 66
gi|491764        76 TIADKQ 82
""",
        )
        hit = hits[204]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|1666185|emb|CAB04766.1|")
        self.assertEqual(hit.target.name, "CAB04766")
        self.assertEqual(
            hit.target.description,
            "ORF13(1) [Rhodococcus erythropolis] >gi|9979018|sp|P72267|TATA_RHOER Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=98)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 71.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 31.9574)
        self.assertAlmostEqual(hsp.annotations["evalue"], 5.95088)
        self.assertEqual(hsp.annotations["identity"], 25)
        self.assertEqual(hsp.annotations["positive"], 35)
        self.assertEqual(hsp.annotations["gaps"], 3)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 40, 43, 60],
                          [14, 54, 54, 71]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 60))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDA'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGAMSPWHWAIVALVVVILFGSKKLPDAARGLGRSLRIFKSEVKEMQNDNSTPAPTAQSA'}, length=98)",
        )
        self.assertEqual(hsp.target.id, "gi|1666185|emb|CAB04766.1|")
        self.assertEqual(hsp.target.name, "CAB04766")
        self.assertEqual(
            hsp.target.description,
            "ORF13(1) [Rhodococcus erythropolis] >gi|9979018|sp|P72267|TATA_RHOER Sec-independent protein translocase protein tatA/E homolog",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MG +S W   I+A++VV+LFG+KKL      LG S++ FK   K M +D      T+Q A",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|166618         0 MGAMSPWHWAIVALVVVILFGSKKLPDAARGLGRSLRIFKSEVKEMQNDNSTPAPTAQSA
                  0 ||..|.|............||.|||......||.|...||---|.|..|......|.|.|
gi|491764        14 MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFK---KAMSDDEPKQDKTSQDA

gi|166618        60 
                 60 
gi|491764        71 
""",
        )
        hit = hits[205]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|72138252|ref|XP_800288.1|")
        self.assertEqual(hit.target.name, "XP_800288")
        self.assertEqual(
            hit.target.description,
            "PREDICTED: hypothetical protein XP_795195 [Strongylocentrotus purpuratus]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=946)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 71.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 31.9574)
        self.assertAlmostEqual(hsp.annotations["evalue"], 5.95088)
        self.assertEqual(hsp.annotations["identity"], 13)
        self.assertEqual(hsp.annotations["positive"], 25)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[854, 898],
                          [ 58, 102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 44))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({58: 'DDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({854: 'EDQDKEDKDEEDKDEEDKDKEDKDAEEEEEEDNADEEEDADKDE'}, length=946)",
        )
        self.assertEqual(hsp.target.id, "gi|72138252|ref|XP_800288.1|")
        self.assertEqual(hsp.target.name, "XP_800288")
        self.assertEqual(
            hsp.target.description,
            "PREDICTED: hypothetical protein XP_795195 [Strongylocentrotus purpuratus]",
        )
        self.assertEqual(
            hsp.annotations["midline"], "+D+ K+DK  +D D   K   DK A+  +E+   ++ +  DK++"
        )
        self.assertEqual(
            str(hsp),
            """\
gi|721382       854 EDQDKEDKDEEDKDEEDKDKEDKDAEEEEEEDNADEEEDADKDE 898
                  0 .|..|.||...|.|...|...||.|....|..........||..  44
gi|491764        58 DDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQ 102
""",
        )
        hit = hits[206]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|67923190|ref|ZP_00516678.1|")
        self.assertEqual(hit.target.name, "ZP_00516678")
        self.assertEqual(
            hit.target.description,
            "Twin-arginine translocation protein TatA/E [Crocosphaera watsonii WH 8501] >gi|67854976|gb|EAM50247.1| Twin-arginine translocation protein TatA/E [Crocosphaera watsonii WH 8501]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=50)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 70.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 31.5722)
        self.assertAlmostEqual(hsp.annotations["evalue"], 7.7721)
        self.assertEqual(hsp.annotations["identity"], 11)
        self.assertEqual(hsp.annotations["positive"], 23)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[19, 50],
                          [33, 64]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 31))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({33: 'FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({19: 'FGPKKIPELGSSLGKTLRGFKEGVNNPKDEE'}, length=50)",
        )
        self.assertEqual(hsp.target.id, "gi|67923190|ref|ZP_00516678.1|")
        self.assertEqual(hsp.target.name, "ZP_00516678")
        self.assertEqual(
            hsp.target.description,
            "Twin-arginine translocation protein TatA/E [Crocosphaera watsonii WH 8501] >gi|67854976|gb|EAM50247.1| Twin-arginine translocation protein TatA/E [Crocosphaera watsonii WH 8501]",
        )
        self.assertEqual(hsp.annotations["midline"], "FG KK+  +GS LG +++GFK+ +++ + ++")
        self.assertEqual(
            str(hsp),
            """\
gi|679231        19 FGPKKIPELGSSLGKTLRGFKEGVNNPKDEE 50
                  0 ||.||....||.||....|||.......... 31
gi|491764        33 FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQ 64
""",
        )
        hit = hits[207]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|3329623|gb|AAC26930.1|")
        self.assertEqual(hit.target.name, "AAC26930")
        self.assertEqual(
            hit.target.description,
            "Hypothetical protein F36H12.3 [Caenorhabditis elegans] >gi|17540204|ref|NP_500767.1| F36H12.3 [Caenorhabditis elegans]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=335)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 70.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 31.5722)
        self.assertAlmostEqual(hsp.annotations["evalue"], 7.7721)
        self.assertEqual(hsp.annotations["identity"], 19)
        self.assertEqual(hsp.annotations["positive"], 31)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 74, 141],
                          [ 35, 102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 67))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({35: 'TKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQ...KEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({74: 'SKKSEGKKSDKKEEKKEEKKEEKEDKKEEKKEEKKEDDKKDEKKDEKKDEKEED...KEK'}, length=335)",
        )
        self.assertEqual(hsp.target.id, "gi|3329623|gb|AAC26930.1|")
        self.assertEqual(hsp.target.name, "AAC26930")
        self.assertEqual(
            hsp.target.description,
            "Hypothetical protein F36H12.3 [Caenorhabditis elegans] >gi|17540204|ref|NP_500767.1| F36H12.3 [Caenorhabditis elegans]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "+KK     SD     K  KK   +D+ ++ K  +  D       D++ D  +E  K+ED K  +KE+",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|332962        74 SKKSEGKKSDKKEEKKEEKKEEKEDKKEEKKEEKKEDDKKDEKKDEKKDEKEEDKKSEDK
                  0 .||.....||.....|..||....|.....|.....|.......|...|...|..|.||.
gi|491764        35 TKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDA

gi|332962       134 KEDEKEK 141
                 60 |...||.  67
gi|491764        95 KRHDKEQ 102
""",
        )
        hit = hits[208]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|39597929|emb|CAE68621.1|")
        self.assertEqual(hit.target.name, "CAE68621")
        self.assertEqual(
            hit.target.description,
            "Hypothetical protein CBG14505 [Caenorhabditis briggsae]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=2691)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 70.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 31.5722)
        self.assertAlmostEqual(hsp.annotations["evalue"], 7.7721)
        self.assertEqual(hsp.annotations["identity"], 16)
        self.assertEqual(hsp.annotations["positive"], 27)
        self.assertEqual(hsp.annotations["gaps"], 4)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[1678, 1705, 1709, 1730],
                          [  54,   81,   81,  102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 52))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({54: 'KAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({1678: 'KDESEEKNKEDNMEESKDFTEMTVDNTISLGECSTTQYEEKTIDIDEHEEEQ'}, length=2691)",
        )
        self.assertEqual(hsp.target.id, "gi|39597929|emb|CAE68621.1|")
        self.assertEqual(hsp.target.name, "CAE68621")
        self.assertEqual(
            hsp.target.description,
            "Hypothetical protein CBG14505 [Caenorhabditis briggsae]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "K  S+++ K+D   +  DFT  T+ +     +  T Q + KT D   H++EQ",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|395979      1678 KDESEEKNKEDNMEESKDFTEMTVDNTISLGECSTTQYEEKTIDIDEHEEEQ 1730
                  0 |..|....|.|......|||..|....----...|.|...||.|...|..||   52
gi|491764        54 KAMSDDEPKQDKTSQDADFTAKTIADK----QADTNQEQAKTEDAKRHDKEQ  102
""",
        )
        hit = hits[209]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|68182025|ref|ZP_00555006.1|")
        self.assertEqual(hit.target.name, "ZP_00555006")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein JannDRAFT_2282 [Jannaschia sp. CCS1] >gi|67977679|gb|EAM67298.1| hypothetical protein JannDRAFT_2282 [Jannaschia sp. CCS1]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=438)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 70.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 31.5722)
        self.assertAlmostEqual(hsp.annotations["evalue"], 7.7721)
        self.assertEqual(hsp.annotations["identity"], 14)
        self.assertEqual(hsp.annotations["positive"], 28)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[234, 288],
                          [ 38,  92]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 54))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({38: 'LGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKT'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({234: 'VGAVDQDMDPEARDTDASAADDETISDETTYDDDAEAEAAHDETVEDDTEQTET'}, length=438)",
        )
        self.assertEqual(hsp.target.id, "gi|68182025|ref|ZP_00555006.1|")
        self.assertEqual(hsp.target.name, "ZP_00555006")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein JannDRAFT_2282 [Jannaschia sp. CCS1] >gi|67977679|gb|EAM67298.1| hypothetical protein JannDRAFT_2282 [Jannaschia sp. CCS1]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "+G++  D+    +    + +DDE   D+T+ D D  A+   D+  + + EQ +T",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|681820       234 VGAVDQDMDPEARDTDASAADDETISDETTYDDDAEAEAAHDETVEDDTEQTET 288
                  0 .|....|.............|||...|.|..|.|..|....|.......||..|  54
gi|491764        38 LGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKT  92
""",
        )
        hit = hits[210]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|21204492|dbj|BAB95189.1|")
        self.assertEqual(hit.target.name, "BAB95189")
        self.assertEqual(
            hit.target.description,
            "ebh [Staphylococcus aureus subsp. aureus MW2] >gi|21283053|ref|NP_646141.1| hypothetical protein MW1324 [Staphylococcus aureus subsp. aureus MW2]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=9904)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 70.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 31.5722)
        self.assertAlmostEqual(hsp.annotations["evalue"], 7.7721)
        self.assertEqual(hsp.annotations["identity"], 19)
        self.assertEqual(hsp.annotations["positive"], 36)
        self.assertEqual(hsp.annotations["gaps"], 12)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[5269, 5295, 5305, 5313, 5315, 5342],
                          [  34,   60,   60,   68,   68,   95]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 73))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({34: 'GTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQE...EDA'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({5269: 'GVNQVSTTASELNTAMSNLQRGINDEAATKAAQKYTDADRDKQTAYNDAVTAAK...EQA'}, length=9904)",
        )
        self.assertEqual(hsp.target.id, "gi|21204492|dbj|BAB95189.1|")
        self.assertEqual(hsp.target.name, "BAB95189")
        self.assertEqual(
            hsp.target.description,
            "ebh [Staphylococcus aureus subsp. aureus MW2] >gi|21283053|ref|NP_646141.1| hypothetical protein MW1324 [Staphylococcus aureus subsp. aureus MW2]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "G  ++ +  S+L  ++   ++ ++D+          +  +DK +   DA   AKT+ DK A TN+ +A  E A",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|212044      5269 GVNQVSTTASELNTAMSNLQRGINDEAATKAAQKYTDADRDKQTAYNDAVTAAKTLLDKT
                  0 |........|.|............|.----------....||..--.||...|||..||.
gi|491764        34 GTKKLGSIGSDLGASIKGFKKAMSDD----------EPKQDKTS--QDADFTAKTIADKQ

gi|212044      5329 AGTNENKAAVEQA 5342
                 60 |.||...|..|.|   73
gi|491764        82 ADTNQEQAKTEDA   95
""",
        )
        hit = hits[211]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|39593039|emb|CAE64508.1|")
        self.assertEqual(hit.target.name, "CAE64508")
        self.assertEqual(
            hit.target.description,
            "Hypothetical protein CBG09238 [Caenorhabditis briggsae]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=960)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 70.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 31.5722)
        self.assertAlmostEqual(hsp.annotations["evalue"], 7.7721)
        self.assertEqual(hsp.annotations["identity"], 19)
        self.assertEqual(hsp.annotations["positive"], 33)
        self.assertEqual(hsp.annotations["gaps"], 4)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[409, 420, 423, 441, 441, 459],
                          [ 54,  65,  65,  83,  84, 102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.query.annotations["frame"], 0)
        self.assertEqual(hsp.target.annotations["frame"], 0)
        self.assertEqual(hsp.shape, (2, 51))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({54: 'KAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQ'}, length=103)",
        )
        self.assertEqual(hsp.query.id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(
            hsp.query.description,
            "component of Sec-independent translocase [Escherichia coli K12]",
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({409: 'KKEADDKAKKDLEAKTKKEADEKAKKEADEKAKKEAEAKTKEAEAKTKKE'}, length=960)",
        )
        self.assertEqual(hsp.target.id, "gi|39593039|emb|CAE64508.1|")
        self.assertEqual(hsp.target.name, "CAE64508")
        self.assertEqual(
            hsp.target.description,
            "Hypothetical protein CBG09238 [Caenorhabditis briggsae]",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "K  +DD+ K+D   KT ++AD  AK  AD++A   + +AKT++A+   K++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|395930       409 KKEADDKAKKDLEAKTKKEADEKAKKEADEKA-KKEAEAKTKEAEAKTKKE 459
                  0 |...||..|.|---||...||..||..||..|-.....|||..|....|..  51
gi|491764        54 KAMSDDEPKQD---KTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQ 102
""",
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
