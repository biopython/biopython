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


class TestBlastp(unittest.TestCase):
    """Test the Blast XML parser for blastp output."""

    def test_xml_2218_blastp_002(self):
        """Parsing BLASTP 2.2.18+ (xml_2218_blastp_002.xml)."""
        filename = "xml_2218_blastp_002.xml"
        datafile = os.path.join("Blast", filename)
        with open(datafile, "rb") as handle:
            records = Blast.parse(handle)
            self.assertEqual(records.program, "blastp")
            self.assertEqual(records.version, "BLASTP 2.2.18+")
            self.assertEqual(
                records.reference,
                'Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schäffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.',
            )
            self.assertEqual(records.db, "gpipe/9606/Previous/protein")
            self.assertIsInstance(records.query, SeqRecord)
            self.assertEqual(records.query.id, "gi|585505|sp|Q08386|MOPB_RHOCA")
            self.assertEqual(
                records.query.description,
                "Molybdenum-pterin-binding protein mopB >gi|310278|gb|AAA71913.1| molybdenum-pterin-binding protein",
            )
            self.assertEqual(repr(records.query.seq), "Seq(None, length=270)")
            self.assertEqual(len(records.param), 5)
            self.assertEqual(records.param["matrix"], "BLOSUM62")
            self.assertAlmostEqual(records.param["expect"], 0.01)
            self.assertEqual(records.param["gap-open"], 11)
            self.assertEqual(records.param["gap-extend"], 1)
            self.assertEqual(records.param["filter"], "m L; R -d repeat/repeat_9606;")
            record = next(records)
            self.assertIsInstance(record.query, SeqRecord)
            self.assertEqual(record.query.id, "gi|585505|sp|Q08386|MOPB_RHOCA")
            self.assertEqual(
                record.query.description,
                "Molybdenum-pterin-binding protein mopB >gi|310278|gb|AAA71913.1| molybdenum-pterin-binding protein",
            )
            self.assertEqual(repr(record.query.seq), "Seq(None, length=270)")

            self.assertEqual(len(record.stat), 7)
            self.assertEqual(record.stat["db-num"], 27252)
            self.assertEqual(record.stat["db-len"], 13958303)
            self.assertEqual(record.stat["hsp-len"], 0)
            self.assertAlmostEqual(record.stat["eff-space"], 0.0)
            self.assertAlmostEqual(record.stat["kappa"], 0.041)
            self.assertAlmostEqual(record.stat["lambda"], 0.267)
            self.assertAlmostEqual(record.stat["entropy"], 0.14)
            self.assertEqual(len(record), 0)
            record = next(records)
            self.assertIsInstance(record.query, SeqRecord)
            self.assertEqual(record.query.id, "gi|129628|sp|P07175.1|PARA_AGRTU")
            self.assertEqual(record.query.description, "Protein parA")
            self.assertEqual(repr(record.query.seq), "Seq(None, length=222)")

            self.assertEqual(len(record.stat), 7)
            self.assertEqual(record.stat["db-num"], 27252)
            self.assertEqual(record.stat["db-len"], 13958303)
            self.assertEqual(record.stat["hsp-len"], 0)
            self.assertAlmostEqual(record.stat["eff-space"], 0.0)
            self.assertAlmostEqual(record.stat["kappa"], 0.041)
            self.assertAlmostEqual(record.stat["lambda"], 0.267)
            self.assertAlmostEqual(record.stat["entropy"], 0.14)
            self.assertEqual(len(record), 0)

    def test_xml_2218L_blastp_001(self):
        """Parsing blastp 2.2.18 [Mar-02-2008] (xml_2218L_blastp_001.xml)."""
        filename = "xml_2218L_blastp_001.xml"
        datafile = os.path.join("Blast", filename)
        with open(datafile, "rb") as handle:
            records = Blast.parse(handle)
            self.check_xml_2218L_blastp_001_records(records)

        with Blast.parse(datafile) as records:
            self.check_xml_2218L_blastp_001_records(records)

        with open(datafile, "rb") as handle:
            record = Blast.read(handle)
        self.check_xml_2218L_blastp_001_record(record)

        record = Blast.read(datafile)
        self.check_xml_2218L_blastp_001_record(record)

    def check_xml_2218L_blastp_001_records(self, records):
        self.assertEqual(records.program, "blastp")
        self.assertEqual(records.version, "blastp 2.2.18 [Mar-02-2008]")
        self.assertEqual(
            records.reference,
            '~Reference: Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schaffer, ~Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), ~"Gapped BLAST and PSI-BLAST: a new generation of protein database search~programs",  Nucleic Acids Res. 25:3389-3402.',
        )
        self.assertEqual(
            records.db, "/Users/pjcock/Downloads/Software/blast-2.2.18/data/nr"
        )
        self.assertIsInstance(records.query, SeqRecord)
        self.assertEqual(records.query.id, "lcl|1_0")
        self.assertEqual(records.query.description, "Fake")
        self.assertEqual(repr(records.query.seq), "Seq(None, length=9)")
        self.assertEqual(len(records.param), 5)
        self.assertEqual(records.param["matrix"], "BLOSUM62")
        self.assertAlmostEqual(records.param["expect"], 1e-05)
        self.assertEqual(records.param["gap-open"], 11)
        self.assertEqual(records.param["gap-extend"], 1)
        self.assertEqual(records.param["filter"], "F")
        record = next(records)
        self.assertRaises(StopIteration, next, records)
        self.check_xml_2218L_blastp_001_record(record)

    def check_xml_2218L_blastp_001_record(self, record):
        self.assertEqual(record.num, 1)
        self.assertIsNone(record.query)
        self.assertEqual(len(record.stat), 7)
        self.assertEqual(record.stat["db-num"], 6589360)
        self.assertEqual(record.stat["db-len"], 2253133281)
        self.assertEqual(record.stat["hsp-len"], 0)
        self.assertAlmostEqual(record.stat["eff-space"], 20278200000.0)
        self.assertAlmostEqual(record.stat["kappa"], 0.041)
        self.assertAlmostEqual(record.stat["lambda"], 0.267)
        self.assertAlmostEqual(record.stat["entropy"], 0.14)
        self.assertEqual(len(record), 0)

    def test_xml_2226_blastp_003(self):
        """Parsing BLASTP 2.2.26+ (xml_2226_blastp_003.xml)."""
        filename = "xml_2226_blastp_003.xml"
        datafile = os.path.join("Blast", filename)
        with open(datafile, "rb") as handle:
            records = Blast.parse(handle)
            self.check_xml_2226_blastp_003(records)
        with Blast.parse(datafile) as records:
            self.check_xml_2226_blastp_003(records)

    def check_xml_2226_blastp_003(self, records):
        self.assertEqual(records.program, "blastp")
        self.assertEqual(records.version, "BLASTP 2.2.26+")
        self.assertEqual(
            records.reference,
            'Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.',
        )
        self.assertEqual(records.db, "db/minirefseq_prot")
        self.assertIsInstance(records.query, SeqRecord)
        self.assertEqual(records.query.id, "Query_1")
        self.assertEqual(
            records.query.description,
            "gi|16080617|ref|NP_391444.1| membrane bound lipoprotein [Bacillus subtilis subsp. subtilis str. 168]",
        )
        self.assertEqual(repr(records.query.seq), "Seq(None, length=102)")
        self.assertEqual(len(records.param), 5)
        self.assertEqual(records.param["matrix"], "BLOSUM62")
        self.assertAlmostEqual(records.param["expect"], 10.0)
        self.assertEqual(records.param["gap-open"], 11)
        self.assertEqual(records.param["gap-extend"], 1)
        self.assertEqual(records.param["filter"], "F")
        record = next(records)
        self.assertIsInstance(record.query, SeqRecord)
        self.assertEqual(record.query.id, "Query_1")
        self.assertEqual(
            record.query.description,
            "gi|16080617|ref|NP_391444.1| membrane bound lipoprotein [Bacillus subtilis subsp. subtilis str. 168]",
        )
        self.assertEqual(repr(record.query.seq), "Seq(None, length=102)")

        self.assertEqual(len(record.stat), 7)
        self.assertEqual(record.stat["db-num"], 20)
        self.assertEqual(record.stat["db-len"], 6406)
        self.assertEqual(record.stat["hsp-len"], 38)
        self.assertAlmostEqual(record.stat["eff-space"], 361344.0)
        self.assertAlmostEqual(record.stat["kappa"], 0.041)
        self.assertAlmostEqual(record.stat["lambda"], 0.267)
        self.assertAlmostEqual(record.stat["entropy"], 0.14)
        self.assertEqual(len(record), 5)
        hit = record[0]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gnl|BL_ORD_ID|1")
        self.assertEqual(hit.target.name, "1")
        self.assertEqual(
            hit.target.description,
            "gi|308175296|ref|YP_003922001.1| membrane bound lipoprotein [Bacillus amyloliquefaciens DSM 7]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=100)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 350.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 139.428)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.99275e-46)
        self.assertEqual(hsp.annotations["identity"], 69)
        self.assertEqual(hsp.annotations["positive"], 81)
        self.assertEqual(hsp.annotations["gaps"], 2)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0,  30,  30, 100],
                              [  0,  30,  32, 102]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 102))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MKKFIALLFFILLLSGCGVNSQKSQGEDVSPDSNIETKEGTYVGLADTHTIEVT...RAN')",
        )
        self.assertEqual(hsp.query.id, "Query_1")
        self.assertEqual(
            hsp.query.description,
            "gi|16080617|ref|NP_391444.1| membrane bound lipoprotein [Bacillus subtilis subsp. subtilis str. 168]",
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MKKIFGCLFFILLLAGCGVTNEKSQGEDAGEKLVTKEGTYVGLADTHTIEVTVD...PAN')",
        )
        self.assertEqual(hsp.target.id, "gnl|BL_ORD_ID|1")
        self.assertEqual(hsp.target.name, "1")
        self.assertEqual(
            hsp.target.description,
            "gi|308175296|ref|YP_003922001.1| membrane bound lipoprotein [Bacillus amyloliquefaciens DSM 7]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MKK    LFFILLL+GCGV ++KSQGED      + TKEGTYVGLADTHTIEVTVD+EPVS DITEES  D+   N+G+KVT+ Y+KN +GQL+LKDIE AN",
        )
        self.assertEqual(
            str(hsp),
            """\
gnl|BL_OR         0 MKKIFGCLFFILLLAGCGVTNEKSQGEDAG--EKLVTKEGTYVGLADTHTIEVTVDHEPV
                  0 |||....|||||||.||||...||||||..--....||||||||||||||||||||.|||
Query_1           0 MKKFIALLFFILLLSGCGVNSQKSQGEDVSPDSNIETKEGTYVGLADTHTIEVTVDNEPV

gnl|BL_OR        58 SFDITEESADDVKNLNNGEKVTVKYQKNSKGQLVLKDIEPAN 100
                 60 |.||||||..|....|.|.|||..|.||..|||.|||||.|| 102
Query_1          60 SLDITEESTSDLDKFNSGDKVTITYEKNDEGQLLLKDIERAN 102
""",
        )
        hit = record[1]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gnl|BL_ORD_ID|2")
        self.assertEqual(hit.target.name, "2")
        self.assertEqual(
            hit.target.description,
            "gi|375363999|ref|YP_005132038.1| lytA gene product [Bacillus amyloliquefaciens subsp. plantarum CAU B946]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=105)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 219.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 88.9669)
        self.assertAlmostEqual(hsp.annotations["evalue"], 6.94052e-27)
        self.assertEqual(hsp.annotations["identity"], 48)
        self.assertEqual(hsp.annotations["positive"], 69)
        self.assertEqual(hsp.annotations["gaps"], 5)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0,  13,  17,  32,  32, 104],
                              [  0,  13,  13,  28,  29, 101]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 105))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({0: 'MKKFIALLFFILLLSGCGVNSQKSQGEDVSPDSNIETKEGTYVGLADTHTIEVT...ERA'}, length=102)",
        )
        self.assertEqual(hsp.query.id, "Query_1")
        self.assertEqual(
            hsp.query.description,
            "gi|16080617|ref|NP_391444.1| membrane bound lipoprotein [Bacillus subtilis subsp. subtilis str. 168]",
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MKKTIAASFLILLFSVVLAACGTAEQSKKGSGSSENQAQKETAYYVGMADTHTI...EKA'}, length=105)",
        )
        self.assertEqual(hsp.target.id, "gnl|BL_ORD_ID|2")
        self.assertEqual(hsp.target.name, "2")
        self.assertEqual(
            hsp.target.description,
            "gi|375363999|ref|YP_005132038.1| lytA gene product [Bacillus amyloliquefaciens subsp. plantarum CAU B946]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MKK IA  F ILL    L+ CG   Q  +G   S ++  + +   YVG+ADTHTIEV VD++PVS + +++ +  L+KF+  DKV+ITY  ND+GQ  +K+IE+A",
        )
        self.assertEqual(
            str(hsp),
            """\
gnl|BL_OR         0 MKKTIAASFLILLFSVVLAACGTAEQSKKGSG-SSENQAQKETAYYVGMADTHTIEVKVD
                  0 |||.||..|.|||----|..||...|...|..-|...........|||.||||||||.||
Query_1           0 MKKFIALLFFILL----LSGCGVNSQKSQGEDVSPDSNIETKEGTYVGLADTHTIEVTVD

gnl|BL_OR        59 DQPVSFEFSDDFSDVLNKFSENDKVSITYFTNDKGQKEIKEIEKA 104
                 60 ..|||..........|.||...|||.|||..||.||...|.||.| 105
Query_1          56 NEPVSLDITEESTSDLDKFNSGDKVTITYEKNDEGQLLLKDIERA 101
""",
        )
        hit = record[2]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gnl|BL_ORD_ID|3")
        self.assertEqual(hit.target.name, "3")
        self.assertEqual(
            hit.target.description,
            "gi|154687679|ref|YP_001422840.1| LytA [Bacillus amyloliquefaciens FZB42]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=105)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 219.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 88.9669)
        self.assertAlmostEqual(hsp.annotations["evalue"], 8.41012e-27)
        self.assertEqual(hsp.annotations["identity"], 48)
        self.assertEqual(hsp.annotations["positive"], 69)
        self.assertEqual(hsp.annotations["gaps"], 5)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0,  13,  17,  32,  32, 104],
                              [  0,  13,  13,  28,  29, 101]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 105))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({0: 'MKKFIALLFFILLLSGCGVNSQKSQGEDVSPDSNIETKEGTYVGLADTHTIEVT...ERA'}, length=102)",
        )
        self.assertEqual(hsp.query.id, "Query_1")
        self.assertEqual(
            hsp.query.description,
            "gi|16080617|ref|NP_391444.1| membrane bound lipoprotein [Bacillus subtilis subsp. subtilis str. 168]",
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MKKTIAASFLILLFSVVLAACGTADQSKKGSGSSENQAQKETAYYVGMADTHTI...EKA'}, length=105)",
        )
        self.assertEqual(hsp.target.id, "gnl|BL_ORD_ID|3")
        self.assertEqual(hsp.target.name, "3")
        self.assertEqual(
            hsp.target.description,
            "gi|154687679|ref|YP_001422840.1| LytA [Bacillus amyloliquefaciens FZB42]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MKK IA  F ILL    L+ CG   Q  +G   S ++  + +   YVG+ADTHTIEV VD++PVS + +++ +  L+KF+  DKV+ITY  ND+GQ  +K+IE+A",
        )
        self.assertEqual(
            str(hsp),
            """\
gnl|BL_OR         0 MKKTIAASFLILLFSVVLAACGTADQSKKGSG-SSENQAQKETAYYVGMADTHTIEVKVD
                  0 |||.||..|.|||----|..||...|...|..-|...........|||.||||||||.||
Query_1           0 MKKFIALLFFILL----LSGCGVNSQKSQGEDVSPDSNIETKEGTYVGLADTHTIEVTVD

gnl|BL_OR        59 DQPVSFEFSDDFSDVLNKFSENDKVSITYFTNDKGQKEIKEIEKA 104
                 60 ..|||..........|.||...|||.|||..||.||...|.||.| 105
Query_1          56 NEPVSLDITEESTSDLDKFNSGDKVTITYEKNDEGQLLLKDIERA 101
""",
        )
        hit = record[3]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gnl|BL_ORD_ID|4")
        self.assertEqual(hit.target.name, "4")
        self.assertEqual(
            hit.target.description,
            "gi|311070071|ref|YP_003974994.1| unnamed protein product [Bacillus atrophaeus 1942]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=105)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 204.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 83.1889)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.37847e-24)
        self.assertEqual(hsp.annotations["identity"], 45)
        self.assertEqual(hsp.annotations["positive"], 66)
        self.assertEqual(hsp.annotations["gaps"], 5)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0,  13,  17,  30,  30, 103],
                              [  0,  13,  13,  26,  27, 100]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 104))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({0: 'MKKFIALLFFILLLSGCGVNSQKSQGEDVSPDSNIETKEGTYVGLADTHTIEVT...IER'}, length=102)",
        )
        self.assertEqual(hsp.query.id, "Query_1")
        self.assertEqual(
            hsp.query.description,
            "gi|16080617|ref|NP_391444.1| membrane bound lipoprotein [Bacillus subtilis subsp. subtilis str. 168]",
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MKKNVASSFLILLFSIILAACGTAEQSKEGNGSSSSQVQNETAYYVGMADTHTI...IEK'}, length=105)",
        )
        self.assertEqual(hsp.target.id, "gnl|BL_ORD_ID|4")
        self.assertEqual(hsp.target.name, "4")
        self.assertEqual(
            hsp.target.description,
            "gi|311070071|ref|YP_003974994.1| unnamed protein product [Bacillus atrophaeus 1942]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MKK +A  F ILL    L+ CG   Q  +G + S  S ++ +   YVG+ADTHTIEV +D++PVS + T++ +  L++F   DKV I+Y  ND+GQ  L +IE+",
        )
        self.assertEqual(
            str(hsp),
            """\
gnl|BL_OR         0 MKKNVASSFLILLFSIILAACGTAEQSKEG-NGSSSSQVQNETAYYVGMADTHTIEVKID
                  0 |||..|..|.|||----|..||...|...|-..|..|........|||.||||||||..|
Query_1           0 MKKFIALLFFILL----LSGCGVNSQKSQGEDVSPDSNIETKEGTYVGLADTHTIEVTVD

gnl|BL_OR        59 DQPVSFEFTDDFSEILNEFEENDKVNISYLTNDKGQKELTEIEK 103
                 60 ..|||...|......|..|...|||.|.|..||.||..|..||. 104
Query_1          56 NEPVSLDITEESTSDLDKFNSGDKVTITYEKNDEGQLLLKDIER 100
""",
        )
        hit = record[4]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gnl|BL_ORD_ID|15")
        self.assertEqual(hit.target.name, "15")
        self.assertEqual(
            hit.target.description,
            "gi|332258565|ref|XP_003278367.1| PREDICTED: UPF0764 protein C16orf89-like [Nomascus leucogenys]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=132)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 29.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 15.779)
        self.assertAlmostEqual(hsp.annotations["evalue"], 7.12269)
        self.assertEqual(hsp.annotations["identity"], 7)
        self.assertEqual(hsp.annotations["positive"], 11)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 79, 104],
                              [ 59,  84]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 25))
        self.assertEqual(
            repr(hsp.query.seq), "Seq({59: 'VSLDITEESTSDLDKFNSGDKVTIT'}, length=102)"
        )
        self.assertEqual(hsp.query.id, "Query_1")
        self.assertEqual(
            hsp.query.description,
            "gi|16080617|ref|NP_391444.1| membrane bound lipoprotein [Bacillus subtilis subsp. subtilis str. 168]",
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq), "Seq({79: 'VEMGFLHVGQAGLELVTSGDPPTLT'}, length=132)"
        )
        self.assertEqual(hsp.target.id, "gnl|BL_ORD_ID|15")
        self.assertEqual(hsp.target.name, "15")
        self.assertEqual(
            hsp.target.description,
            "gi|332258565|ref|XP_003278367.1| PREDICTED: UPF0764 protein C16orf89-like [Nomascus leucogenys]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(hsp.annotations["midline"], "V +       + L+   SGD  T+T")
        self.assertEqual(
            str(hsp),
            """\
gnl|BL_OR        79 VEMGFLHVGQAGLELVTSGDPPTLT 104
                  0 |...........|....|||..|.|  25
Query_1          59 VSLDITEESTSDLDKFNSGDKVTIT  84
""",
        )

    def test_xml_2218L_rpsblast_001(self):
        """Parsing blastp 2.2.18 [Mar-02-2008] (xml_2218L_rpsblast_001.xml)."""
        filename = "xml_2218L_rpsblast_001.xml"
        datafile = os.path.join("Blast", filename)
        with open(datafile, "rb") as handle:
            records = Blast.parse(handle)
            self.check_xml_2218L_rpsblast_001(records)
        with Blast.parse(datafile) as records:
            self.check_xml_2218L_rpsblast_001(records)

    def check_xml_2218L_rpsblast_001(self, records):
        self.assertEqual(records.program, "blastp")
        self.assertEqual(records.version, "blastp 2.2.18 [Mar-02-2008]")
        self.assertEqual(
            records.reference,
            '~Reference: Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schaffer, ~Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), ~"Gapped BLAST and PSI-BLAST: a new generation of protein database search~programs",  Nucleic Acids Res. 25:3389-3402.',
        )
        self.assertEqual(records.db, "/opt/BlastDBs/nr")
        self.assertIsInstance(records.query, SeqRecord)
        self.assertEqual(records.query.id, "lcl|QUERY")
        self.assertEqual(records.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(repr(records.query.seq), "Seq(None, length=131)")
        self.assertEqual(len(records.param), 6)
        self.assertEqual(records.param["matrix"], "BLOSUM62")
        self.assertAlmostEqual(records.param["expect"], 10.0)
        self.assertAlmostEqual(records.param["include"], 0.001)
        self.assertEqual(records.param["gap-open"], 11)
        self.assertEqual(records.param["gap-extend"], 1)
        self.assertEqual(records.param["filter"], "S")
        record = next(records)
        self.assertIsNone(record.query)
        self.assertEqual(len(record.stat), 7)
        self.assertEqual(record.stat["db-num"], 2563094)
        self.assertEqual(record.stat["db-len"], 864488805)
        self.assertEqual(record.stat["hsp-len"], 96)
        self.assertAlmostEqual(record.stat["eff-space"], 80278100000.0)
        self.assertAlmostEqual(record.stat["kappa"], 0.041)
        self.assertAlmostEqual(record.stat["lambda"], 0.267)
        self.assertAlmostEqual(record.stat["entropy"], 0.14)
        self.assertEqual(len(record), 11)
        hit = record[0]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|75750454|ref|YP_319893.1|")
        self.assertEqual(hit.target.name, "YP_319893")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein ATV_gp62 [Acidianus two-tailed virus] >gi|74474837|emb|CAI59911.1| hypothetical protein [Acidianus two-tailed virus]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=131)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 680.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 266.544)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.72196e-70)
        self.assertEqual(hsp.annotations["identity"], 131)
        self.assertEqual(hsp.annotations["positive"], 131)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0, 131],
                              [  0, 131]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 131))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEA...MAS')",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEA...MAS')",
        )
        self.assertEqual(hsp.target.id, "gi|75750454|ref|YP_319893.1|")
        self.assertEqual(hsp.target.name, "YP_319893")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein ATV_gp62 [Acidianus two-tailed virus] >gi|74474837|emb|CAI59911.1| hypothetical protein [Acidianus two-tailed virus]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQPIPGEVIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSKGAKKLAQAMAS",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|757504         0 MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVK
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
lcl|QUERY         0 MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVK

gi|757504        60 MISDAMKPYRNKGSGFQSQPIPGEVIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSK
                 60 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
lcl|QUERY        60 MISDAMKPYRNKGSGFQSQPIPGEVIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSK

gi|757504       120 GAKKLAQAMAS 131
                120 ||||||||||| 131
lcl|QUERY       120 GAKKLAQAMAS 131
""",
        )
        hit = record[1]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|51980166|ref|YP_077233.1|")
        self.assertEqual(hit.target.name, "YP_077233")
        self.assertEqual(
            hit.target.description,
            "coat protein [Sulfolobus virus STSV1] >gi|51890299|emb|CAH04223.1| coat protein [Sulfolobus virus STSV1]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=144)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 181.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 74.3294)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.03476e-12)
        self.assertEqual(hsp.annotations["identity"], 36)
        self.assertEqual(hsp.annotations["positive"], 49)
        self.assertEqual(hsp.annotations["gaps"], 1)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 70, 70, 76],
                              [ 0, 70, 71, 77]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 77))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({0: 'MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEA...GFQ'}, length=131)",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MAREEPYKGDYVGGVAKILQGYFANYYGFPNVSLRLAGEEANLSKTGHANAKAI...GFK'}, length=144)",
        )
        self.assertEqual(hsp.target.id, "gi|51980166|ref|YP_077233.1|")
        self.assertEqual(hsp.target.name, "YP_077233")
        self.assertEqual(
            hsp.target.description,
            "coat protein [Sulfolobus virus STSV1] >gi|51890299|emb|CAH04223.1| coat protein [Sulfolobus virus STSV1]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MA+ EP KGDY GG  KIL  +     G+P V+L+LAGEEAN  + G    K  +H ++K+I +A KP R +G GF+",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|519801         0 MAREEPYKGDYVGGVAKILQGYFANYYGFPNVSLRLAGEEANLSKTGHANAKAIVHEMIK
                  0 ||..||.||||.||..|||........|.|.|.|.|||||||....|....|...|...|
lcl|QUERY         0 MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVK

gi|519801        60 VIKEASKPLR-RGKGFK 76
                 60 .|..|.||.|-.|.||. 77
lcl|QUERY        60 MISDAMKPYRNKGSGFQ 77
""",
        )
        hit = record[2]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|70606448|ref|YP_255318.1|")
        self.assertEqual(hit.target.name, "YP_255318")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein Saci_0636 [Sulfolobus acidocaldarius DSM 639] >gi|68567096|gb|AAY80025.1| hypothetical protein Saci_0636 [Sulfolobus acidocaldarius DSM 639]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=90)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 106.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 45.4394)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.00184433)
        self.assertEqual(hsp.annotations["identity"], 22)
        self.assertEqual(hsp.annotations["positive"], 37)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 65],
                              [ 0, 65]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 65))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({0: 'MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEA...SDA'}, length=131)",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MMELKHKKEEYAKELINLLSRFMSGKISYPVLSLRLTGIEAKVHESGFEDLVRL...REA'}, length=90)",
        )
        self.assertEqual(hsp.target.id, "gi|70606448|ref|YP_255318.1|")
        self.assertEqual(hsp.target.name, "YP_255318")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein Saci_0636 [Sulfolobus acidocaldarius DSM 639] >gi|68567096|gb|AAY80025.1| hypothetical protein Saci_0636 [Sulfolobus acidocaldarius DSM 639]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "M + + KK +YA   + +L  F +G++ YP ++L+L G EA    +G E     IH ++K I +A",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|706064         0 MMELKHKKEEYAKELINLLSRFMSGKISYPVLSLRLTGIEAKVHESGFEDLVRLIHELLK
                  0 |.....||..||......|..|..|...||...|.|.|.||.....|.|.....||...|
lcl|QUERY         0 MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVK

gi|706064        60 AIREA 65
                 60 .|..| 65
lcl|QUERY        60 MISDA 65
""",
        )
        hit = record[3]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|75750440|ref|YP_319873.1|")
        self.assertEqual(hit.target.name, "YP_319873")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein ATV_gp42 [Acidianus two-tailed virus] >gi|74474823|emb|CAI59897.1| hypothetical protein [Acidianus two-tailed virus]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=145)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 102.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 43.8986)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.00527737)
        self.assertEqual(hsp.annotations["identity"], 30)
        self.assertEqual(hsp.annotations["positive"], 42)
        self.assertEqual(hsp.annotations["gaps"], 6)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[20, 50, 50, 72, 72, 87, 89, 93],
                              [21, 51, 52, 74, 77, 92, 92, 96]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 77))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({21: 'FENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSG...EYQ'}, length=131)",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({20: 'FYEGVIGYPEIDLRLAGEEAWLKGVNPELAEAVKKIIKTIRRYLEGSPYDGSEK...EVQ'}, length=145)",
        )
        self.assertEqual(hsp.target.id, "gi|75750440|ref|YP_319873.1|")
        self.assertEqual(hsp.target.name, "YP_319873")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein ATV_gp42 [Acidianus two-tailed virus] >gi|74474823|emb|CAI59897.1| hypothetical protein [Acidianus two-tailed virus]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "F  G +GYPE+ L+LAGEEA  +    E   EA+  I+K I   ++     GS    +PIP  +IA++ S   PE Q",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|757504        20 FYEGVIGYPEIDLRLAGEEAWLKGVNPELA-EAVKKIIKTIRRYLEGSPYDGS---EKPI
                  0 |..|..||||..|.||||||.......|..-||...|.|.|..........||---..||
lcl|QUERY        21 FENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQPI

gi|757504        76 PRYIIAEIFSQIAPEVQ 93
                 60 |...||...|.--||.| 77
lcl|QUERY        81 PGEVIAQVTSN--PEYQ 96
""",
        )
        hit = record[4]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|149720683|ref|XP_001495622.1|")
        self.assertEqual(hit.target.name, "XP_001495622")
        self.assertEqual(
            hit.target.description,
            "PREDICTED: similar to P2X7 receptor [Equus caballus]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=595)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 79.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 35.039)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.35127)
        self.assertEqual(hsp.annotations["identity"], 17)
        self.assertEqual(hsp.annotations["positive"], 30)
        self.assertEqual(hsp.annotations["gaps"], 2)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 40,  69,  71,  87],
                              [ 84, 113, 113, 129]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 47))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({84: 'VIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSKGAKKLAQAM'}, length=131)",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({40: 'VIFALVSDKRYQRKEPLISSVHTKVKGIAEVREEIIESGAKKVVQSV'}, length=595)",
        )
        self.assertEqual(hsp.target.id, "gi|149720683|ref|XP_001495622.1|")
        self.assertEqual(hsp.target.name, "XP_001495622")
        self.assertEqual(
            hsp.target.description,
            "PREDICTED: similar to P2X7 receptor [Equus caballus]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "VI  + S+  YQ+ +  ++S  T+V+ I   REE++  GAKK+ Q++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|149720        40 VIFALVSDKRYQRKEPLISSVHTKVKGIAEVREEIIESGAKKVVQSV  87
                  0 ||....|...||.......|..|.|..|.--|||....||||..|..  47
lcl|QUERY        84 VIAQVTSNPEYQQAKAFLASPATQVRNIE--REEVLSKGAKKLAQAM 129
""",
        )
        hit = record[5]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|229220336|ref|ZP_03838482.2|")
        self.assertEqual(hit.target.name, "ZP_03838482")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein CIT292_04594 [Citrobacter youngae ATCC 29220] >gi|228553887|gb|EEK18552.1| hypothetical protein CIT292_04594 [Citrobacter youngae ATCC 29220]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=208)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 78.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 34.6538)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.99494)
        self.assertEqual(hsp.annotations["identity"], 23)
        self.assertEqual(hsp.annotations["positive"], 42)
        self.assertEqual(hsp.annotations["gaps"], 9)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 59,  71,  71,  98, 102, 138],
                              [  8,  20,  25,  52,  52,  88]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 84))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({8: 'GDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMI...IAQ'}, length=131)",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({59: 'GEVVGGAISTALSIGGPQIAKGVSKTERNARRAVQEEAKNLSGESLNSFMKLNE...VSQ'}, length=208)",
        )
        self.assertEqual(hsp.target.id, "gi|229220336|ref|ZP_03838482.2|")
        self.assertEqual(hsp.target.name, "ZP_03838482")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein CIT292_04594 [Citrobacter youngae ATCC 29220] >gi|228553887|gb|EEK18552.1| hypothetical protein CIT292_04594 [Citrobacter youngae ATCC 29220]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "G+  GGA+          +G P++   ++  E NARRA  E  K    E++++ +K+    + P+RN     Q  P  G++++Q",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|229220        59 GEVVGGAISTAL-----SIGGPQIAKGVSKTERNARRAVQEEAKNLSGESLNSFMKLNEK
                  0 |...|||.....-----..|.|.........|.|||||..|..|----|......|....
lcl|QUERY         8 GDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTK----EAIHAIVKMISD

gi|229220       114 KLTPHRNFEILGQCSPAIGQLVSQ 138
                 60 ...|.||.....|..|..|....|  84
lcl|QUERY        64 AMKPYRNKGSGFQSQPIPGEVIAQ  88
""",
        )
        hit = record[6]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|167646253|ref|YP_001683916.1|")
        self.assertEqual(hit.target.name, "YP_001683916")
        self.assertEqual(
            hit.target.description,
            "activator of Hsp90 ATPase 1 family protein [Caulobacter sp. K31] >gi|167348683|gb|ABZ71418.1| Activator of Hsp90 ATPase 1 family protein [Caulobacter sp. K31]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=155)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 78.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 34.6538)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.12253)
        self.assertEqual(hsp.annotations["identity"], 16)
        self.assertEqual(hsp.annotations["positive"], 24)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 44,  92],
                              [ 58, 106]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 48))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({58: 'VKMISDAMKPYRNKGSGFQSQPIPGEVIAQVTSNPEYQQAKAFLASPA'}, length=131)",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({44: 'VKMLVKTTSDRRRKDSGRITQLVPGEAITFTTGNPMLKHTECYTLTPS'}, length=155)",
        )
        self.assertEqual(hsp.target.id, "gi|167646253|ref|YP_001683916.1|")
        self.assertEqual(hsp.target.name, "YP_001683916")
        self.assertEqual(
            hsp.target.description,
            "activator of Hsp90 ATPase 1 family protein [Caulobacter sp. K31] >gi|167348683|gb|ABZ71418.1| Activator of Hsp90 ATPase 1 family protein [Caulobacter sp. K31]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "VKM+       R K SG  +Q +PGE I   T NP  +  + +  +P+",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|167646        44 VKMLVKTTSDRRRKDSGRITQLVPGEAITFTTGNPMLKHTECYTLTPS  92
                  0 |||........|.|.||...|..|||.|...|.||...........|.  48
lcl|QUERY        58 VKMISDAMKPYRNKGSGFQSQPIPGEVIAQVTSNPEYQQAKAFLASPA 106
""",
        )
        hit = record[7]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|46114422|ref|XP_383229.1|")
        self.assertEqual(hit.target.name, "XP_383229")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein FG03053.1 [Gibberella zeae PH-1]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=263)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 77.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 34.2686)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.07814)
        self.assertEqual(hsp.annotations["identity"], 26)
        self.assertEqual(hsp.annotations["positive"], 40)
        self.assertEqual(hsp.annotations["gaps"], 5)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 85, 118, 118, 127, 127, 151, 151, 169],
                              [ 34,  67,  68,  77,  79, 103, 105, 123]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 89))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({34: 'KLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQPIPGEVIAQ...GAK'}, length=131)",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({85: 'KLTTLKRRLRYEAKELRKEALRQIVSKLRSWMKLRDFGCGVKPLIESSLPEIAT...GSR'}, length=263)",
        )
        self.assertEqual(hsp.target.id, "gi|46114422|ref|XP_383229.1|")
        self.assertEqual(hsp.target.name, "XP_383229")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein FG03053.1 [Gibberella zeae PH-1]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "KL   +   R    E  KEA+  IV  +   MK  R+ G G +  P+    + ++ +NP Y Q   +L       R IER   + KG++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|461144        85 KLTTLKRRLRYEAKELRKEALRQIVSKLRSWMK-LRDFGCGVK--PLIESSLPEIATNPH
                  0 ||.......|....|..|||...||......||-.|..|.|..--|...........||.
lcl|QUERY        34 KLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQPIPGEVIAQVTSNPE

gi|461144       142 YNQTSLYLT--MNHPRRIERVSRVIKGSR 169
                 60 |.|....|.--....|.|||.....||..  89
lcl|QUERY        94 YQQAKAFLASPATQVRNIEREEVLSKGAK 123
""",
        )
        hit = record[8]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|16804655|ref|NP_466140.1|")
        self.assertEqual(hit.target.name, "NP_466140")
        self.assertEqual(
            hit.target.description,
            "50S ribosomal protein L6 [Listeria monocytogenes EGD-e] >gi|47097188|ref|ZP_00234753.1| ribosomal protein L6 [Listeria monocytogenes str. 1/2a F6854] >gi|217966184|ref|YP_002351862.1| 50S ribosomal protein L6 (BL10) [Listeria monocytogenes HCC23] >gi|224499290|ref|ZP_03667639.1| 50S ribosomal protein L6 [Listeria monocytogenes Finland 1988] >gi|224503589|ref|ZP_03671896.1| 50S ribosomal protein L6 [Listeria monocytogenes FSL R2-561] >gi|81768735|sp|Q8Y444.1|RL6_LISMO RecName: Full=50S ribosomal protein L6 >gi|16412105|emb|CAD00695.1| ribosomal protein L6 [Listeria monocytogenes] >gi|47014447|gb|EAL05415.1| ribosomal protein L6 [Listeria monocytogenes str. 1/2a F6854] >gi|217335454|gb|ACK41248.1| 50S ribosomal protein L6 (BL10) [Listeria monocytogenes HCC23]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=178)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 77.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 34.2686)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.14677)
        self.assertEqual(hsp.annotations["identity"], 26)
        self.assertEqual(hsp.annotations["positive"], 44)
        self.assertEqual(hsp.annotations["gaps"], 13)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 38,  75,  80,  86,  90, 108, 108, 138],
                              [ 29,  66,  66,  72,  72,  90,  94, 124]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 104))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({29: 'PEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQPIPG...AKK'}, length=131)",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({38: 'PEITIKIEGNEINVSRPTDNKNHRALHGTTRAILNNMVVGVSEGYEKKLELIGV...YNK'}, length=178)",
        )
        self.assertEqual(hsp.target.id, "gi|16804655|ref|NP_466140.1|")
        self.assertEqual(hsp.target.name, "NP_466140")
        self.assertEqual(
            hsp.target.description,
            "50S ribosomal protein L6 [Listeria monocytogenes EGD-e] >gi|47097188|ref|ZP_00234753.1| ribosomal protein L6 [Listeria monocytogenes str. 1/2a F6854] >gi|217966184|ref|YP_002351862.1| 50S ribosomal protein L6 (BL10) [Listeria monocytogenes HCC23] >gi|224499290|ref|ZP_03667639.1| 50S ribosomal protein L6 [Listeria monocytogenes Finland 1988] >gi|224503589|ref|ZP_03671896.1| 50S ribosomal protein L6 [Listeria monocytogenes FSL R2-561] >gi|81768735|sp|Q8Y444.1|RL6_LISMO RecName: Full=50S ribosomal protein L6 >gi|16412105|emb|CAD00695.1| ribosomal protein L6 [Listeria monocytogenes] >gi|47014447|gb|EAL05415.1| ribosomal protein L6 [Listeria monocytogenes str. 1/2a F6854] >gi|217335454|gb|ACK41248.1| 50S ribosomal protein L6 (BL10) [Listeria monocytogenes HCC23]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "PE+T+K+ G E N  R  D +   A+H   + I + M     + Y  K    G G+++Q    +++  V     Y     F+A     +      +V+ KG  K",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|168046        38 PEITIKIEGNEINVSRPTDNKNHRALHGTTRAILNNMVVGVSEGYEKKLELIGVGYRAQK
                  0 ||.|.|..|.|.|..|..|.....|.|.....|...|-----..|..|----|.|...|.
lcl|QUERY        29 PEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAM-----KPYRNK----GSGFQSQP

gi|168046        98 QGDKLVLNVG----YSHPVEFVAPKGVDIEVPANTQVIVKGYNK 138
                 60 ........|.----|.....|.|.............|..||..| 104
lcl|QUERY        80 IPGEVIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSKGAKK 124
""",
        )
        hit = record[9]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|225106043|ref|YP_002674349.1|")
        self.assertEqual(hit.target.name, "YP_002674349")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein OA307_717 [Octadecabacter antarcticus 307] >gi|198251380|gb|EDY75695.1| hypothetical protein OA307_717 [Octadecabacter antarcticus 307]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=202)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 76.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 33.8834)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.81876)
        self.assertEqual(hsp.annotations["identity"], 23)
        self.assertEqual(hsp.annotations["positive"], 38)
        self.assertEqual(hsp.annotations["gaps"], 5)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 8, 16, 21, 92],
                              [ 0,  8,  8, 79]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 84))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({0: 'MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEA...QSQ'}, length=131)",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({8: 'MRKIEPKKMTATERKHWGFGVDAIRMMENELMGYLHEYRALPLEWDEIWKDEDR...QAR'}, length=202)",
        )
        self.assertEqual(hsp.target.id, "gi|225106043|ref|YP_002674349.1|")
        self.assertEqual(hsp.target.name, "YP_002674349")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein OA307_717 [Octadecabacter antarcticus 307] >gi|198251380|gb|EDY75695.1| hypothetical protein OA307_717 [Octadecabacter antarcticus 307]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "M K EPKK       + G  V  + M EN  +GY      L  E     +  D R  +     +++ +D +K +R  G+G+Q++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|225106         8 MRKIEPKKMTATERKHWGFGVDAIRMMENELMGYLHEYRALPLEWDEIWKDEDRRDPKRT
                  0 |.|.||||-----....|..|....|.||...||......|..|........|.|.....
lcl|QUERY         0 MAKYEPKK-----GDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAI

gi|225106        68 RVTIRLDADVVKYFRAMGAGYQAR 92
                 60 ........|..|..|..|.|.|.. 84
lcl|QUERY        55 HAIVKMISDAMKPYRNKGSGFQSQ 79
""",
        )
        hit = record[10]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|16801827|ref|NP_472095.1|")
        self.assertEqual(hit.target.name, "NP_472095")
        self.assertEqual(
            hit.target.description,
            "50S ribosomal protein L6 [Listeria innocua Clip11262] >gi|116873983|ref|YP_850764.1| 50S ribosomal protein L6 [Listeria welshimeri serovar 6b str. SLCC5334] >gi|81774192|sp|Q927M2.1|RL6_LISIN RecName: Full=50S ribosomal protein L6 >gi|123458730|sp|A0ALV3.1|RL6_LISW6 RecName: Full=50S ribosomal protein L6 >gi|16415302|emb|CAC97992.1| ribosomal protein L6 [Listeria innocua] >gi|116742861|emb|CAK21985.1| rplF [Listeria welshimeri serovar 6b str. SLCC5334]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=178)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 74.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 33.113)
        self.assertAlmostEqual(hsp.annotations["evalue"], 9.08516)
        self.assertEqual(hsp.annotations["identity"], 25)
        self.assertEqual(hsp.annotations["positive"], 44)
        self.assertEqual(hsp.annotations["gaps"], 13)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 38,  75,  80,  86,  90, 108, 108, 138],
                              [ 29,  66,  66,  72,  72,  90,  94, 124]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 104))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({29: 'PEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQPIPG...AKK'}, length=131)",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({38: 'PEITINIEGNEINVSRPTDNKNHRALHGTTRAILNNMVVGVSEGYEKKLELIGV...YNK'}, length=178)",
        )
        self.assertEqual(hsp.target.id, "gi|16801827|ref|NP_472095.1|")
        self.assertEqual(hsp.target.name, "NP_472095")
        self.assertEqual(
            hsp.target.description,
            "50S ribosomal protein L6 [Listeria innocua Clip11262] >gi|116873983|ref|YP_850764.1| 50S ribosomal protein L6 [Listeria welshimeri serovar 6b str. SLCC5334] >gi|81774192|sp|Q927M2.1|RL6_LISIN RecName: Full=50S ribosomal protein L6 >gi|123458730|sp|A0ALV3.1|RL6_LISW6 RecName: Full=50S ribosomal protein L6 >gi|16415302|emb|CAC97992.1| ribosomal protein L6 [Listeria innocua] >gi|116742861|emb|CAK21985.1| rplF [Listeria welshimeri serovar 6b str. SLCC5334]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "PE+T+ + G E N  R  D +   A+H   + I + M     + Y  K    G G+++Q    +++  V     Y     F+A    ++      +V+ KG  K",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|168018        38 PEITINIEGNEINVSRPTDNKNHRALHGTTRAILNNMVVGVSEGYEKKLELIGVGYRAQK
                  0 ||.|....|.|.|..|..|.....|.|.....|...|-----..|..|----|.|...|.
lcl|QUERY        29 PEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAM-----KPYRNK----GSGFQSQP

gi|168018        98 QGDKLVLNVG----YSHPVEFVAPKGVEIEVPANTQVIVKGYNK 138
                 60 ........|.----|.....|.|.............|..||..| 104
lcl|QUERY        80 IPGEVIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSKGAKK 124
""",
        )
        record = next(records)
        self.assertIsNone(record.query)
        self.assertEqual(len(record.stat), 7)
        self.assertEqual(record.stat["db-num"], 2563094)
        self.assertEqual(record.stat["db-len"], 864488805)
        self.assertEqual(record.stat["hsp-len"], 96)
        self.assertAlmostEqual(record.stat["eff-space"], 80278100000.0)
        self.assertAlmostEqual(record.stat["kappa"], 0.0413587)
        self.assertAlmostEqual(record.stat["lambda"], 0.267)
        self.assertAlmostEqual(record.stat["entropy"], 0.14)
        self.assertEqual(len(record), 19)
        hit = record[0]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|75750454|ref|YP_319893.1|")
        self.assertEqual(hit.target.name, "YP_319893")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein ATV_gp62 [Acidianus two-tailed virus] >gi|74474837|emb|CAI59911.1| hypothetical protein [Acidianus two-tailed virus]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=131)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 590.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 231.867)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.28615e-59)
        self.assertEqual(hsp.annotations["identity"], 131)
        self.assertEqual(hsp.annotations["positive"], 131)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0, 131],
                              [  0, 131]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 131))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEA...MAS')",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEA...MAS')",
        )
        self.assertEqual(hsp.target.id, "gi|75750454|ref|YP_319893.1|")
        self.assertEqual(hsp.target.name, "YP_319893")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein ATV_gp62 [Acidianus two-tailed virus] >gi|74474837|emb|CAI59911.1| hypothetical protein [Acidianus two-tailed virus]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQPIPGEVIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSKGAKKLAQAMAS",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|757504         0 MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVK
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
lcl|QUERY         0 MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVK

gi|757504        60 MISDAMKPYRNKGSGFQSQPIPGEVIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSK
                 60 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
lcl|QUERY        60 MISDAMKPYRNKGSGFQSQPIPGEVIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSK

gi|757504       120 GAKKLAQAMAS 131
                120 ||||||||||| 131
lcl|QUERY       120 GAKKLAQAMAS 131
""",
        )
        hit = record[1]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|51980166|ref|YP_077233.1|")
        self.assertEqual(hit.target.name, "YP_077233")
        self.assertEqual(
            hit.target.description,
            "coat protein [Sulfolobus virus STSV1] >gi|51890299|emb|CAH04223.1| coat protein [Sulfolobus virus STSV1]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=144)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 346.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 137.878)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.63528e-31)
        self.assertEqual(hsp.annotations["identity"], 36)
        self.assertEqual(hsp.annotations["positive"], 49)
        self.assertEqual(hsp.annotations["gaps"], 1)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 70, 70, 76],
                              [ 0, 70, 71, 77]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 77))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({0: 'MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEA...GFQ'}, length=131)",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MAREEPYKGDYVGGVAKILQGYFANYYGFPNVSLRLAGEEANLSKTGHANAKAI...GFK'}, length=144)",
        )
        self.assertEqual(hsp.target.id, "gi|51980166|ref|YP_077233.1|")
        self.assertEqual(hsp.target.name, "YP_077233")
        self.assertEqual(
            hsp.target.description,
            "coat protein [Sulfolobus virus STSV1] >gi|51890299|emb|CAH04223.1| coat protein [Sulfolobus virus STSV1]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MA+ EP KGDY GG  KIL  +     G+P V+L+LAGEEAN  + G    K  +H ++K+I +A KP R +G GF+",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|519801         0 MAREEPYKGDYVGGVAKILQGYFANYYGFPNVSLRLAGEEANLSKTGHANAKAIVHEMIK
                  0 ||..||.||||.||..|||........|.|.|.|.|||||||....|....|...|...|
lcl|QUERY         0 MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVK

gi|519801        60 VIKEASKPLR-RGKGFK 76
                 60 .|..|.||.|-.|.||. 77
lcl|QUERY        60 MISDAMKPYRNKGSGFQ 77
""",
        )
        hit = record[2]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|75750440|ref|YP_319873.1|")
        self.assertEqual(hit.target.name, "YP_319873")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein ATV_gp42 [Acidianus two-tailed virus] >gi|74474823|emb|CAI59897.1| hypothetical protein [Acidianus two-tailed virus]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=145)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 118.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 50.0528)
        self.assertAlmostEqual(hsp.annotations["evalue"], 7.16661e-05)
        self.assertEqual(hsp.annotations["identity"], 34)
        self.assertEqual(hsp.annotations["positive"], 48)
        self.assertEqual(hsp.annotations["gaps"], 8)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 6,  9,  9, 50, 50, 72, 72, 87, 89, 93],
                              [ 5,  8, 10, 51, 52, 74, 77, 92, 92, 96]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 93))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({5: 'PKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIV...EYQ'}, length=131)",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({6: 'PQKYETHVLDDLMEFYEGVIGYPEIDLRLAGEEAWLKGVNPELAEAVKKIIKTI...EVQ'}, length=145)",
        )
        self.assertEqual(hsp.target.id, "gi|75750440|ref|YP_319873.1|")
        self.assertEqual(hsp.target.name, "YP_319873")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein ATV_gp42 [Acidianus two-tailed virus] >gi|74474823|emb|CAI59897.1| hypothetical protein [Acidianus two-tailed virus]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "P+K  Y    +  L  F  G +GYPE+ L+LAGEEA  +    E   EA+  I+K I   ++     GS    +PIP  +IA++ S   PE Q",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|757504         6 PQK--YETHVLDDLMEFYEGVIGYPEIDLRLAGEEAWLKGVNPELA-EAVKKIIKTIRRY
                  0 |.|--|.......|..|..|..||||..|.||||||.......|..-||...|.|.|...
lcl|QUERY         5 PKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDA

gi|757504        63 LEGSPYDGS---EKPIPRYIIAEIFSQIAPEVQ 93
                 60 .......||---..|||...||...|.--||.| 93
lcl|QUERY        65 MKPYRNKGSGFQSQPIPGEVIAQVTSN--PEYQ 96
""",
        )
        hit = record[3]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|70606448|ref|YP_255318.1|")
        self.assertEqual(hit.target.name, "YP_255318")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein Saci_0636 [Sulfolobus acidocaldarius DSM 639] >gi|68567096|gb|AAY80025.1| hypothetical protein Saci_0636 [Sulfolobus acidocaldarius DSM 639]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=90)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 112.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 47.7416)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.000338308)
        self.assertEqual(hsp.annotations["identity"], 22)
        self.assertEqual(hsp.annotations["positive"], 38)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 67],
                              [ 0, 67]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 67))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({0: 'MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEA...AMK'}, length=131)",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MMELKHKKEEYAKELINLLSRFMSGKISYPVLSLRLTGIEAKVHESGFEDLVRL...AER'}, length=90)",
        )
        self.assertEqual(hsp.target.id, "gi|70606448|ref|YP_255318.1|")
        self.assertEqual(hsp.target.name, "YP_255318")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein Saci_0636 [Sulfolobus acidocaldarius DSM 639] >gi|68567096|gb|AAY80025.1| hypothetical protein Saci_0636 [Sulfolobus acidocaldarius DSM 639]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "M + + KK +YA   + +L  F +G++ YP ++L+L G EA    +G E     IH ++K I +A +",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|706064         0 MMELKHKKEEYAKELINLLSRFMSGKISYPVLSLRLTGIEAKVHESGFEDLVRLIHELLK
                  0 |.....||..||......|..|..|...||...|.|.|.||.....|.|.....||...|
lcl|QUERY         0 MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVK

gi|706064        60 AIREAER 67
                 60 .|..|.. 67
lcl|QUERY        60 MISDAMK 67
""",
        )
        hit = record[4]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|56478923|ref|YP_160512.1|")
        self.assertEqual(hit.target.name, "YP_160512")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein ebA6101 [Aromatoleum aromaticum EbN1] >gi|56314966|emb|CAI09611.1| conserved hypothetical protein [Azoarcus sp. EbN1]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=326)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 85.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 37.3412)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.461067)
        self.assertEqual(hsp.annotations["identity"], 21)
        self.assertEqual(hsp.annotations["positive"], 33)
        self.assertEqual(hsp.annotations["gaps"], 11)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[141, 148, 148, 170, 172, 203, 205, 209, 211, 217],
                              [  8,  15,  20,  42,  42,  73,  73,  77,  77,  83]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 81))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({8: 'GDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMI...IPG'}, length=131)",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({141: 'GDHVGGLGDSAGQMVFPNATIHVAQTDNDFWLSPQMAAQAPAEIQPFFKMAVDA...VPG'}, length=326)",
        )
        self.assertEqual(hsp.target.id, "gi|56478923|ref|YP_160512.1|")
        self.assertEqual(hsp.target.name, "YP_160512")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein ebA6101 [Aromatoleum aromaticum EbN1] >gi|56314966|emb|CAI09611.1| conserved hypothetical protein [Azoarcus sp. EbN1]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "GD+ GG          GQ+ +P  T+ +A  + +         +    I    KM  DA  PY+ +G    F   S+ +PG",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|564789       141 GDHVGGL-----GDSAGQMVFPNATIHVAQTDNDFWLSPQMAAQAPAEIQPFFKMAVDAA
                  0 ||..||.-----....||...|..|...|.....--............|....||..||.
lcl|QUERY         8 GDYAGGAVKILDMFENGQLGYPEVTLKLAGEEAN--ARRAGDERTKEAIHAIVKMISDAM

gi|564789       196 APYQARGQWKTFAEGSEVVPG 217
                 60 .||...|--..|.--|...||  81
lcl|QUERY        66 KPYRNKG--SGFQ--SQPIPG  83
""",
        )
        hit = record[5]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|16082535|ref|NP_394071.1|")
        self.assertEqual(hit.target.name, "NP_394071")
        self.assertEqual(
            hit.target.description,
            "deoxycytidine triphosphate deaminase [Thermoplasma acidophilum DSM 1728] >gi|20141306|sp|Q9HKK0.2|DCD_THEAC RecName: Full=Probable deoxycytidine triphosphate deaminase; Short=dCTP deaminase",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=166)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 82.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 36.1856)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.960826)
        self.assertEqual(hsp.annotations["identity"], 22)
        self.assertEqual(hsp.annotations["positive"], 30)
        self.assertEqual(hsp.annotations["gaps"], 10)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 94, 117, 117, 124, 124, 146, 152, 160],
                              [ 23,  46,  48,  55,  57,  79,  79,  87]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 70))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({23: 'NGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQ...VIA'}, length=131)",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({94: 'AGYHGNLTLSFFNAGSAVNLRRGERIAQIVFVKMIGSAEKPYHIRSGNYQNSRG...QIA'}, length=166)",
        )
        self.assertEqual(hsp.target.id, "gi|16082535|ref|NP_394071.1|")
        self.assertEqual(hsp.target.name, "NP_394071")
        self.assertEqual(
            hsp.target.description,
            "deoxycytidine triphosphate deaminase [Thermoplasma acidophilum DSM 1728] >gi|20141306|sp|Q9HKK0.2|DCD_THEAC RecName: Full=Probable deoxycytidine triphosphate deaminase; Short=dCTP deaminase",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            " G  G   ++   AG   N RR   ER  + +   VKMI  A KPY  +   +Q+       P+ G  IA",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|160825        94 AGYHGNLTLSFFNAGSAVNLRRG--ERIAQIV--FVKMIGSAEKPYHIRSGNYQNSRGIV
                  0 .|..|........||...|.||.--||.....--.||||..|.|||.......|..----
lcl|QUERY        23 NGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQ----

gi|160825       150 KDPVAGSQIA 160
                 60 --|..|..||  70
lcl|QUERY        79 --PIPGEVIA  87
""",
        )
        hit = record[6]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|10639765|emb|CAC11737.1|")
        self.assertEqual(hit.target.name, "CAC11737")
        self.assertEqual(
            hit.target.description,
            "dCTP deaminase related protein [Thermoplasma acidophilum]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=183)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 82.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 36.1856)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.960826)
        self.assertEqual(hsp.annotations["identity"], 22)
        self.assertEqual(hsp.annotations["positive"], 30)
        self.assertEqual(hsp.annotations["gaps"], 10)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[111, 134, 134, 141, 141, 163, 169, 177],
                              [ 23,  46,  48,  55,  57,  79,  79,  87]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 70))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({23: 'NGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQ...VIA'}, length=131)",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({111: 'AGYHGNLTLSFFNAGSAVNLRRGERIAQIVFVKMIGSAEKPYHIRSGNYQNSRG...QIA'}, length=183)",
        )
        self.assertEqual(hsp.target.id, "gi|10639765|emb|CAC11737.1|")
        self.assertEqual(hsp.target.name, "CAC11737")
        self.assertEqual(
            hsp.target.description,
            "dCTP deaminase related protein [Thermoplasma acidophilum]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            " G  G   ++   AG   N RR   ER  + +   VKMI  A KPY  +   +Q+       P+ G  IA",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|106397       111 AGYHGNLTLSFFNAGSAVNLRRG--ERIAQIV--FVKMIGSAEKPYHIRSGNYQNSRGIV
                  0 .|..|........||...|.||.--||.....--.||||..|.|||.......|..----
lcl|QUERY        23 NGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQ----

gi|106397       167 KDPVAGSQIA 177
                 60 --|..|..||  70
lcl|QUERY        79 --PIPGEVIA  87
""",
        )
        hit = record[7]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|38511961|gb|AAH60735.1|")
        self.assertEqual(hit.target.name, "AAH60735")
        self.assertEqual(hit.target.description, "Ttc28 protein [Mus musculus]")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=1370)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 82.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 36.1856)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.01015)
        self.assertEqual(hsp.annotations["identity"], 25)
        self.assertEqual(hsp.annotations["positive"], 45)
        self.assertEqual(hsp.annotations["gaps"], 5)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[738, 782, 784, 811, 811, 838],
                              [ 26,  70,  70,  97, 100, 127]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 103))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({26: 'LGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQP...LAQ'}, length=131)",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({738: 'LGLPNPALQALCKLITASETGEQLISRAVKNMVGMLHQVLVQLQACEKEQDFAS...ASR'}, length=1370)",
        )
        self.assertEqual(hsp.target.id, "gi|38511961|gb|AAH60735.1|")
        self.assertEqual(hsp.target.name, "AAH60735")
        self.assertEqual(hsp.target.description, "Ttc28 protein [Mus musculus]")
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "LG P   L+   +   A   G++    A+  +V M+   +   +   K   F S PIP  +  Q+   P   +   FLA+    +  + +EEV+ K  K+ ++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|385119       738 LGLPNPALQALCKLITASETGEQLISRAVKNMVGMLHQVLVQLQACEKEQDFASAPIPVS
                  0 ||.|...|........|...|......|....|.|.........--.|...|.|.|||..
lcl|QUERY        26 LGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYR--NKGSGFQSQPIPGE

gi|385119       798 LSVQLWRLPGCHE---FLAALGFDLCEVGQEEVILKTGKQASR 838
                 60 ...|....|....---|||...........|||..|..|.... 103
lcl|QUERY        84 VIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSKGAKKLAQ 127
""",
        )
        hit = record[8]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|71153789|sp|Q80XJ3.2|TTC28_MOUSE")
        self.assertEqual(hit.target.name, "Q80XJ3")
        self.assertEqual(
            hit.target.description,
            "RecName: Full=Tetratricopeptide repeat protein 28; Short=TPR repeat protein 28 >gi|55930902|gb|AAH46779.2| Tetratricopeptide repeat domain 28 [Mus musculus]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=1691)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 82.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 36.1856)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.12589)
        self.assertEqual(hsp.annotations["identity"], 25)
        self.assertEqual(hsp.annotations["positive"], 45)
        self.assertEqual(hsp.annotations["gaps"], 5)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[1059, 1103, 1105, 1132, 1132, 1159],
                              [  26,   70,   70,   97,  100,  127]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 103))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({26: 'LGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQP...LAQ'}, length=131)",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({1059: 'LGLPNPALQALCKLITASETGEQLISRAVKNMVGMLHQVLVQLQACEKEQDFAS...ASR'}, length=1691)",
        )
        self.assertEqual(hsp.target.id, "gi|71153789|sp|Q80XJ3.2|TTC28_MOUSE")
        self.assertEqual(hsp.target.name, "Q80XJ3")
        self.assertEqual(
            hsp.target.description,
            "RecName: Full=Tetratricopeptide repeat protein 28; Short=TPR repeat protein 28 >gi|55930902|gb|AAH46779.2| Tetratricopeptide repeat domain 28 [Mus musculus]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "LG P   L+   +   A   G++    A+  +V M+   +   +   K   F S PIP  +  Q+   P   +   FLA+    +  + +EEV+ K  K+ ++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|711537      1059 LGLPNPALQALCKLITASETGEQLISRAVKNMVGMLHQVLVQLQACEKEQDFASAPIPVS
                  0 ||.|...|........|...|......|....|.|.........--.|...|.|.|||..
lcl|QUERY        26 LGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYR--NKGSGFQSQPIPGE

gi|711537      1119 LSVQLWRLPGCHE---FLAALGFDLCEVGQEEVILKTGKQASR 1159
                 60 ...|....|....---|||...........|||..|..|....  103
lcl|QUERY        84 VIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSKGAKKLAQ  127
""",
        )
        hit = record[9]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|149254322|ref|XP_001476594.1|")
        self.assertEqual(hit.target.name, "XP_001476594")
        self.assertEqual(
            hit.target.description,
            "PREDICTED: similar to OTTHUMP00000028696 [Mus musculus]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=2481)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 81.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 35.8004)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.23411)
        self.assertEqual(hsp.annotations["identity"], 25)
        self.assertEqual(hsp.annotations["positive"], 45)
        self.assertEqual(hsp.annotations["gaps"], 5)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[1849, 1893, 1895, 1922, 1922, 1949],
                              [  26,   70,   70,   97,  100,  127]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 103))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({26: 'LGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQP...LAQ'}, length=131)",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({1849: 'LGLPNPALQALCKLITASETGEQLISRAVKNMVGMLHQVLVQLQACEKEQDFAS...ASR'}, length=2481)",
        )
        self.assertEqual(hsp.target.id, "gi|149254322|ref|XP_001476594.1|")
        self.assertEqual(hsp.target.name, "XP_001476594")
        self.assertEqual(
            hsp.target.description,
            "PREDICTED: similar to OTTHUMP00000028696 [Mus musculus]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "LG P   L+   +   A   G++    A+  +V M+   +   +   K   F S PIP  +  Q+   P   +   FLA+    +  + +EEV+ K  K+ ++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|149254      1849 LGLPNPALQALCKLITASETGEQLISRAVKNMVGMLHQVLVQLQACEKEQDFASAPIPVS
                  0 ||.|...|........|...|......|....|.|.........--.|...|.|.|||..
lcl|QUERY        26 LGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYR--NKGSGFQSQPIPGE

gi|149254      1909 LSVQLWRLPGCHE---FLAALGFDLCEVGQEEVILKTGKQASR 1949
                 60 ...|....|....---|||...........|||..|..|....  103
lcl|QUERY        84 VIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSKGAKKLAQ  127
""",
        )
        hit = record[10]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|148688047|gb|EDL19994.1|")
        self.assertEqual(hit.target.name, "EDL19994")
        self.assertEqual(
            hit.target.description, "tetratricopeptide repeat domain 28 [Mus musculus]"
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=2146)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 81.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 35.8004)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.24445)
        self.assertEqual(hsp.annotations["identity"], 25)
        self.assertEqual(hsp.annotations["positive"], 45)
        self.assertEqual(hsp.annotations["gaps"], 5)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[1514, 1558, 1560, 1587, 1587, 1614],
                              [  26,   70,   70,   97,  100,  127]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 103))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({26: 'LGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQP...LAQ'}, length=131)",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({1514: 'LGLPNPALQALCKLITASETGEQLISRAVKNMVGMLHQVLVQLQACEKEQDFAS...ASR'}, length=2146)",
        )
        self.assertEqual(hsp.target.id, "gi|148688047|gb|EDL19994.1|")
        self.assertEqual(hsp.target.name, "EDL19994")
        self.assertEqual(
            hsp.target.description, "tetratricopeptide repeat domain 28 [Mus musculus]"
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "LG P   L+   +   A   G++    A+  +V M+   +   +   K   F S PIP  +  Q+   P   +   FLA+    +  + +EEV+ K  K+ ++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|148688      1514 LGLPNPALQALCKLITASETGEQLISRAVKNMVGMLHQVLVQLQACEKEQDFASAPIPVS
                  0 ||.|...|........|...|......|....|.|.........--.|...|.|.|||..
lcl|QUERY        26 LGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYR--NKGSGFQSQPIPGE

gi|148688      1574 LSVQLWRLPGCHE---FLAALGFDLCEVGQEEVILKTGKQASR 1614
                 60 ...|....|....---|||...........|||..|..|....  103
lcl|QUERY        84 VIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSKGAKKLAQ  127
""",
        )
        hit = record[11]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|149063697|gb|EDM14020.1|")
        self.assertEqual(hit.target.name, "EDM14020")
        self.assertEqual(hit.target.description, "rCG21379 [Rattus norvegicus]")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=2098)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 81.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 35.8004)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.41037)
        self.assertEqual(hsp.annotations["identity"], 25)
        self.assertEqual(hsp.annotations["positive"], 45)
        self.assertEqual(hsp.annotations["gaps"], 5)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[1465, 1509, 1511, 1538, 1538, 1565],
                              [  26,   70,   70,   97,  100,  127]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 103))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({26: 'LGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQP...LAQ'}, length=131)",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({1465: 'LGLPNPALQALCKLITASETGEQLISRAVKNMVGMLHQVLVQLQACEKEQDFAS...ASR'}, length=2098)",
        )
        self.assertEqual(hsp.target.id, "gi|149063697|gb|EDM14020.1|")
        self.assertEqual(hsp.target.name, "EDM14020")
        self.assertEqual(hsp.target.description, "rCG21379 [Rattus norvegicus]")
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "LG P   L+   +   A   G++    A+  +V M+   +   +   K   F S PIP  +  Q+   P   +   FLA+    +  + +EEV+ K  K+ ++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|149063      1465 LGLPNPALQALCKLITASETGEQLISRAVKNMVGMLHQVLVQLQACEKEQDFASAPIPVS
                  0 ||.|...|........|...|......|....|.|.........--.|...|.|.|||..
lcl|QUERY        26 LGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYR--NKGSGFQSQPIPGE

gi|149063      1525 LSVQLWRLPGCHE---FLAALGFDLCEVGQEEVILKTGKQASR 1565
                 60 ...|....|....---|||...........|||..|..|....  103
lcl|QUERY        84 VIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSKGAKKLAQ  127
""",
        )
        hit = record[12]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|163854616|ref|YP_001628914.1|")
        self.assertEqual(hit.target.name, "YP_001628914")
        self.assertEqual(
            hit.target.description,
            "metallo-beta-lactamase superfamily protein [Bordetella petrii DSM 12804] >gi|163258344|emb|CAP40643.1| metallo-beta-lactamase superfamily protein [Bordetella petrii]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=333)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 81.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 35.8004)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.45823)
        self.assertEqual(hsp.annotations["identity"], 15)
        self.assertEqual(hsp.annotations["positive"], 23)
        self.assertEqual(hsp.annotations["gaps"], 6)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[162, 176, 178, 216, 220, 223],
                              [ 28,  42,  42,  80,  80,  83]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 61))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({28: 'YPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQPIPG'}, length=131)",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({162: 'FPNADVYAARAEADFWLSTTVAAQAPADVQPMFKMSRDAAAPYQAAGKFLTYQP...LPG'}, length=333)",
        )
        self.assertEqual(hsp.target.id, "gi|163854616|ref|YP_001628914.1|")
        self.assertEqual(hsp.target.name, "YP_001628914")
        self.assertEqual(
            hsp.target.description,
            "metallo-beta-lactamase superfamily protein [Bordetella petrii DSM 12804] >gi|163258344|emb|CAP40643.1| metallo-beta-lactamase superfamily protein [Bordetella petrii]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "+P   +  A  EA+         +    +  + KM  DA  PY+  G     QP    +PG",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|163854       162 FPNADVYAARAEADFWLSTTVAAQAPADVQPMFKMSRDAAAPYQAAGKFLTYQPDDTLLP
                  0 .|......|..||.--.................||..||..||...|.....||----.|
lcl|QUERY        28 YPEVTLKLAGEEAN--ARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQP----IP

gi|163854       222 G 223
                 60 |  61
lcl|QUERY        82 G  83
""",
        )
        hit = record[13]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|109496198|ref|XP_222260.4|")
        self.assertEqual(hit.target.name, "XP_222260")
        self.assertEqual(
            hit.target.description,
            "PREDICTED: similar to TPR repeat-containing protein KIAA1043 [Rattus norvegicus]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=2570)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 81.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 35.8004)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.47045)
        self.assertEqual(hsp.annotations["identity"], 25)
        self.assertEqual(hsp.annotations["positive"], 45)
        self.assertEqual(hsp.annotations["gaps"], 5)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[1937, 1981, 1983, 2010, 2010, 2037],
                              [  26,   70,   70,   97,  100,  127]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 103))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({26: 'LGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQP...LAQ'}, length=131)",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({1937: 'LGLPNPALQALCKLITASETGEQLISRAVKNMVGMLHQVLVQLQACEKEQDFAS...ASR'}, length=2570)",
        )
        self.assertEqual(hsp.target.id, "gi|109496198|ref|XP_222260.4|")
        self.assertEqual(hsp.target.name, "XP_222260")
        self.assertEqual(
            hsp.target.description,
            "PREDICTED: similar to TPR repeat-containing protein KIAA1043 [Rattus norvegicus]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "LG P   L+   +   A   G++    A+  +V M+   +   +   K   F S PIP  +  Q+   P   +   FLA+    +  + +EEV+ K  K+ ++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|109496      1937 LGLPNPALQALCKLITASETGEQLISRAVKNMVGMLHQVLVQLQACEKEQDFASAPIPVS
                  0 ||.|...|........|...|......|....|.|.........--.|...|.|.|||..
lcl|QUERY        26 LGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYR--NKGSGFQSQPIPGE

gi|109496      1997 LSVQLWRLPGCHE---FLAALGFDLCEVGQEEVILKTGKQASR 2037
                 60 ...|....|....---|||...........|||..|..|....  103
lcl|QUERY        84 VIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSKGAKKLAQ  127
""",
        )
        hit = record[14]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|26327341|dbj|BAC27414.1|")
        self.assertEqual(hit.target.name, "BAC27414")
        self.assertEqual(
            hit.target.description, "unnamed protein product [Mus musculus]"
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=682)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 80.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 35.4152)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.69454)
        self.assertEqual(hsp.annotations["identity"], 25)
        self.assertEqual(hsp.annotations["positive"], 45)
        self.assertEqual(hsp.annotations["gaps"], 5)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 50,  94,  96, 123, 123, 150],
                              [ 26,  70,  70,  97, 100, 127]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 103))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({26: 'LGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQP...LAQ'}, length=131)",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({50: 'LGLPNPALQALCKLITASETGEQLISRAVKNMVGMLHQVLVQLQACEKEQDFAS...ASR'}, length=682)",
        )
        self.assertEqual(hsp.target.id, "gi|26327341|dbj|BAC27414.1|")
        self.assertEqual(hsp.target.name, "BAC27414")
        self.assertEqual(
            hsp.target.description, "unnamed protein product [Mus musculus]"
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "LG P   L+   +   A   G++    A+  +V M+   +   +   K   F S PIP  +  Q+   P   +   FLA+    +  + +EEV+ K  K+ ++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|263273        50 LGLPNPALQALCKLITASETGEQLISRAVKNMVGMLHQVLVQLQACEKEQDFASAPIPVS
                  0 ||.|...|........|...|......|....|.|.........--.|...|.|.|||..
lcl|QUERY        26 LGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYR--NKGSGFQSQPIPGE

gi|263273       110 LSVQLWRLPGCHE---FLAALGFDLCEVGQEEVILKTGKQASR 150
                 60 ...|....|....---|||...........|||..|..|.... 103
lcl|QUERY        84 VIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSKGAKKLAQ 127
""",
        )
        hit = record[15]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|26349137|dbj|BAC38208.1|")
        self.assertEqual(hit.target.name, "BAC38208")
        self.assertEqual(
            hit.target.description, "unnamed protein product [Mus musculus]"
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=641)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 80.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 35.4152)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.842)
        self.assertEqual(hsp.annotations["identity"], 25)
        self.assertEqual(hsp.annotations["positive"], 45)
        self.assertEqual(hsp.annotations["gaps"], 5)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  9,  53,  55,  82,  82, 109],
                              [ 26,  70,  70,  97, 100, 127]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 103))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({26: 'LGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQP...LAQ'}, length=131)",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({9: 'LGLPNPALQALCKLITASETGEQLISRAVKNMVGMLHQVLVQLQACEKEQDFAS...ASR'}, length=641)",
        )
        self.assertEqual(hsp.target.id, "gi|26349137|dbj|BAC38208.1|")
        self.assertEqual(hsp.target.name, "BAC38208")
        self.assertEqual(
            hsp.target.description, "unnamed protein product [Mus musculus]"
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "LG P   L+   +   A   G++    A+  +V M+   +   +   K   F S PIP  +  Q+   P   +   FLA+    +  + +EEV+ K  K+ ++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|263491         9 LGLPNPALQALCKLITASETGEQLISRAVKNMVGMLHQVLVQLQACEKEQDFASAPIPVS
                  0 ||.|...|........|...|......|....|.|.........--.|...|.|.|||..
lcl|QUERY        26 LGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYR--NKGSGFQSQPIPGE

gi|263491        69 LSVQLWRLPGCHE---FLAALGFDLCEVGQEEVILKTGKQASR 109
                 60 ...|....|....---|||...........|||..|..|.... 103
lcl|QUERY        84 VIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSKGAKKLAQ 127
""",
        )
        hit = record[16]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|194291978|ref|YP_002007885.1|")
        self.assertEqual(hit.target.name, "YP_002007885")
        self.assertEqual(
            hit.target.description,
            "non ribosomal peptide synthase, antibiotic synthesis; contains 3 condensation domains, 2 AMP-acid ligases II domains, 2 PP-binding, Phosphopantetheine attachment site [Cupriavidus taiwanensis] >gi|193225882|emb|CAQ71828.1| non ribosomal peptide synthase, antibiotic synthesis; contains 3 condensation domains, 2 AMP-acid ligases II domains, 2 PP-binding, Phosphopantetheine attachment site [Cupriavidus taiwanensis]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=2596)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 80.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 35.4152)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.95279)
        self.assertEqual(hsp.annotations["identity"], 28)
        self.assertEqual(hsp.annotations["positive"], 40)
        self.assertEqual(hsp.annotations["gaps"], 20)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[2334, 2340, 2340, 2354, 2360, 2374, 2376, 2395, 2395, 2402, 2402,
                               2417],
                              [  10,   16,   21,   35,   35,   49,   49,   68,   73,   80,   82,
                                 97]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 95))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({10: 'YAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISD...YQQ'}, length=131)",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({2334: 'YVGGPALARGYLGRAALTADRFVPDPLAGDGARLYRTGDRVRRRADGVFEYLGR...VRQ'}, length=2596)",
        )
        self.assertEqual(hsp.target.id, "gi|194291978|ref|YP_002007885.1|")
        self.assertEqual(hsp.target.name, "YP_002007885")
        self.assertEqual(
            hsp.target.description,
            "non ribosomal peptide synthase, antibiotic synthesis; contains 3 condensation domains, 2 AMP-acid ligases II domains, 2 PP-binding, Phosphopantetheine attachment site [Cupriavidus taiwanensis] >gi|193225882|emb|CAQ71828.1| non ribosomal peptide synthase, antibiotic synthesis; contains 3 condensation domains, 2 AMP-acid ligases II domains, 2 PP-binding, Phosphopantetheine attachment site [Cupriavidus taiwanensis]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "Y GG          G LG   +T        LAG+ A   R GD   R  + +   +  + D +K       GF+ +P  GEV AQV + P+ +Q",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|194291      2334 YVGGPA-----LARGYLGRAALTADRFVPDPLAGDGARLYRTGDRVRRRADGVFEYLGRV
                  0 |.||..-----...|.||....|..------|||..|...|.||.--|............
lcl|QUERY        10 YAGGAVKILDMFENGQLGYPEVTLK------LAGEEANARRAGDE--RTKEAIHAIVKMI

gi|194291      2389 DDQVKI-----RGFRVEP--GEVAAQVAALPQVRQ 2417
                 60 .|..|.-----.||...|--|||.|||...|...|   95
lcl|QUERY        62 SDAMKPYRNKGSGFQSQPIPGEVIAQVTSNPEYQQ   97
""",
        )
        hit = record[17]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|164498981|gb|ABY59058.1|")
        self.assertEqual(hit.target.name, "ABY59058")
        self.assertEqual(hit.target.description, "unknown [Mesorhizobium sp. F28]")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=246)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 77.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 34.2596)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.90317)
        self.assertEqual(hsp.annotations["identity"], 14)
        self.assertEqual(hsp.annotations["positive"], 21)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 7, 55],
                              [ 7, 55]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 48))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({7: 'KGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAI'}, length=131)",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({7: 'RGILVGLVAGILAFAFARVYGEPQVDKAIAFEEQQAQAAGEAPEEEIV'}, length=246)",
        )
        self.assertEqual(hsp.target.id, "gi|164498981|gb|ABY59058.1|")
        self.assertEqual(hsp.target.name, "ABY59058")
        self.assertEqual(hsp.target.description, "unknown [Mesorhizobium sp. F28]")
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "+G   G    IL        G P+V   +A EE  A+ AG+   +E +",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|164498         7 RGILVGLVAGILAFAFARVYGEPQVDKAIAFEEQQAQAAGEAPEEEIV 55
                  0 .|...|....||........|.|.|....|.||..|..||.....|.. 48
lcl|QUERY         7 KGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAI 55
""",
        )
        hit = record[18]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|189198988|ref|XP_001935831.1|")
        self.assertEqual(hit.target.name, "XP_001935831")
        self.assertEqual(
            hit.target.description,
            "predicted protein [Pyrenophora tritici-repentis Pt-1C-BFP] >gi|187982930|gb|EDU48418.1| predicted protein [Pyrenophora tritici-repentis Pt-1C-BFP]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=1026)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 76.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 33.8744)
        self.assertAlmostEqual(hsp.annotations["evalue"], 5.35939)
        self.assertEqual(hsp.annotations["identity"], 14)
        self.assertEqual(hsp.annotations["positive"], 26)
        self.assertEqual(hsp.annotations["gaps"], 4)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[436, 463, 467, 492],
                              [ 28,  55,  55,  80]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 56))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({28: 'YPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQP'}, length=131)",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({436: 'FPEVTMRGGGEESGLNLTQEGNTEDMVPPDYPESFPISSESSDYASEYGASFQSPP'}, length=1026)",
        )
        self.assertEqual(hsp.target.id, "gi|189198988|ref|XP_001935831.1|")
        self.assertEqual(hsp.target.name, "XP_001935831")
        self.assertEqual(
            hsp.target.description,
            "predicted protein [Pyrenophora tritici-repentis Pt-1C-BFP] >gi|187982930|gb|EDU48418.1| predicted protein [Pyrenophora tritici-repentis Pt-1C-BFP]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "+PEVT++  GEE+      +  T++ +         + S++       G+ FQS P",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|189198       436 FPEVTMRGGGEESGLNLTQEGNTEDMVPPDYPESFPISSESSDYASEYGASFQSPP 492
                  0 .||||....|||..........|....----.......|.........|..|||.|  56
lcl|QUERY        28 YPEVTLKLAGEEANARRAGDERTKEAI----HAIVKMISDAMKPYRNKGSGFQSQP  80
""",
        )
        record = next(records)
        self.assertIsNone(record.query)
        self.assertEqual(len(record.stat), 7)
        self.assertEqual(record.stat["db-num"], 2563094)
        self.assertEqual(record.stat["db-len"], 864488805)
        self.assertEqual(record.stat["hsp-len"], 96)
        self.assertAlmostEqual(record.stat["eff-space"], 80278100000.0)
        self.assertAlmostEqual(record.stat["kappa"], 0.0433765)
        self.assertAlmostEqual(record.stat["lambda"], 0.267)
        self.assertAlmostEqual(record.stat["entropy"], 0.14)
        self.assertEqual(len(record), 9)
        hit = record[0]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|75750454|ref|YP_319893.1|")
        self.assertEqual(hit.target.name, "YP_319893")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein ATV_gp62 [Acidianus two-tailed virus] >gi|74474837|emb|CAI59911.1| hypothetical protein [Acidianus two-tailed virus]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=131)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 535.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 210.59)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.43623e-53)
        self.assertEqual(hsp.annotations["identity"], 131)
        self.assertEqual(hsp.annotations["positive"], 131)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0, 131],
                              [  0, 131]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 131))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEA...MAS')",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEA...MAS')",
        )
        self.assertEqual(hsp.target.id, "gi|75750454|ref|YP_319893.1|")
        self.assertEqual(hsp.target.name, "YP_319893")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein ATV_gp62 [Acidianus two-tailed virus] >gi|74474837|emb|CAI59911.1| hypothetical protein [Acidianus two-tailed virus]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQPIPGEVIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSKGAKKLAQAMAS",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|757504         0 MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVK
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
lcl|QUERY         0 MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVK

gi|757504        60 MISDAMKPYRNKGSGFQSQPIPGEVIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSK
                 60 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
lcl|QUERY        60 MISDAMKPYRNKGSGFQSQPIPGEVIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSK

gi|757504       120 GAKKLAQAMAS 131
                120 ||||||||||| 131
lcl|QUERY       120 GAKKLAQAMAS 131
""",
        )
        hit = record[1]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|75750440|ref|YP_319873.1|")
        self.assertEqual(hit.target.name, "YP_319873")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein ATV_gp42 [Acidianus two-tailed virus] >gi|74474823|emb|CAI59897.1| hypothetical protein [Acidianus two-tailed virus]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=145)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 312.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 124.69)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.62891e-27)
        self.assertEqual(hsp.annotations["identity"], 34)
        self.assertEqual(hsp.annotations["positive"], 48)
        self.assertEqual(hsp.annotations["gaps"], 8)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 4,  9,  9, 50, 50, 72, 72, 87, 89, 93],
                              [ 3,  8, 10, 51, 52, 74, 77, 92, 92, 96]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 95))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({3: 'YEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHA...EYQ'}, length=131)",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({4: 'LHPQKYETHVLDDLMEFYEGVIGYPEIDLRLAGEEAWLKGVNPELAEAVKKIIK...EVQ'}, length=145)",
        )
        self.assertEqual(hsp.target.id, "gi|75750440|ref|YP_319873.1|")
        self.assertEqual(hsp.target.name, "YP_319873")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein ATV_gp42 [Acidianus two-tailed virus] >gi|74474823|emb|CAI59897.1| hypothetical protein [Acidianus two-tailed virus]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "  P+K  Y    +  L  F  G +GYPE+ L+LAGEEA  +    E   EA+  I+K I   ++     GS    +PIP  +IA++ S   PE Q",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|757504         4 LHPQK--YETHVLDDLMEFYEGVIGYPEIDLRLAGEEAWLKGVNPELA-EAVKKIIKTIR
                  0 ..|.|--|.......|..|..|..||||..|.||||||.......|..-||...|.|.|.
lcl|QUERY         3 YEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMIS

gi|757504        61 RYLEGSPYDGS---EKPIPRYIIAEIFSQIAPEVQ 93
                 60 .........||---..|||...||...|.--||.| 95
lcl|QUERY        63 DAMKPYRNKGSGFQSQPIPGEVIAQVTSN--PEYQ 96
""",
        )
        hit = record[2]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|51980166|ref|YP_077233.1|")
        self.assertEqual(hit.target.name, "YP_077233")
        self.assertEqual(
            hit.target.description,
            "coat protein [Sulfolobus virus STSV1] >gi|51890299|emb|CAH04223.1| coat protein [Sulfolobus virus STSV1]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=144)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 298.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 119.298)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.0246e-25)
        self.assertEqual(hsp.annotations["identity"], 36)
        self.assertEqual(hsp.annotations["positive"], 49)
        self.assertEqual(hsp.annotations["gaps"], 1)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 70, 70, 77],
                              [ 0, 70, 71, 78]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 78))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({0: 'MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEA...FQS'}, length=131)",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MAREEPYKGDYVGGVAKILQGYFANYYGFPNVSLRLAGEEANLSKTGHANAKAI...FKE'}, length=144)",
        )
        self.assertEqual(hsp.target.id, "gi|51980166|ref|YP_077233.1|")
        self.assertEqual(hsp.target.name, "YP_077233")
        self.assertEqual(
            hsp.target.description,
            "coat protein [Sulfolobus virus STSV1] >gi|51890299|emb|CAH04223.1| coat protein [Sulfolobus virus STSV1]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MA+ EP KGDY GG  KIL  +     G+P V+L+LAGEEAN  + G    K  +H ++K+I +A KP R +G GF+ ",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|519801         0 MAREEPYKGDYVGGVAKILQGYFANYYGFPNVSLRLAGEEANLSKTGHANAKAIVHEMIK
                  0 ||..||.||||.||..|||........|.|.|.|.|||||||....|....|...|...|
lcl|QUERY         0 MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVK

gi|519801        60 VIKEASKPLR-RGKGFKE 77
                 60 .|..|.||.|-.|.||.. 78
lcl|QUERY        60 MISDAMKPYRNKGSGFQS 78
""",
        )
        hit = record[3]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|70606448|ref|YP_255318.1|")
        self.assertEqual(hit.target.name, "YP_255318")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein Saci_0636 [Sulfolobus acidocaldarius DSM 639] >gi|68567096|gb|AAY80025.1| hypothetical protein Saci_0636 [Sulfolobus acidocaldarius DSM 639]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=90)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 242.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 97.7264)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.24619e-19)
        self.assertEqual(hsp.annotations["identity"], 26)
        self.assertEqual(hsp.annotations["positive"], 45)
        self.assertEqual(hsp.annotations["gaps"], 4)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 67, 67, 81],
                              [ 0, 67, 71, 85]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 85))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({0: 'MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEA...GEV'}, length=131)",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MMELKHKKEEYAKELINLLSRFMSGKISYPVLSLRLTGIEAKVHESGFEDLVRL...PYV'}, length=90)",
        )
        self.assertEqual(hsp.target.id, "gi|70606448|ref|YP_255318.1|")
        self.assertEqual(hsp.target.name, "YP_255318")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein Saci_0636 [Sulfolobus acidocaldarius DSM 639] >gi|68567096|gb|AAY80025.1| hypothetical protein Saci_0636 [Sulfolobus acidocaldarius DSM 639]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "M + + KK +YA   + +L  F +G++ YP ++L+L G EA    +G E     IH ++K I +A +    KGS   ++ +   V",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|706064         0 MMELKHKKEEYAKELINLLSRFMSGKISYPVLSLRLTGIEAKVHESGFEDLVRLIHELLK
                  0 |.....||..||......|..|..|...||...|.|.|.||.....|.|.....||...|
lcl|QUERY         0 MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVK

gi|706064        60 AIREAER----KGSNEINKLVRPYV 81
                 60 .|..|..----|||..........| 85
lcl|QUERY        60 MISDAMKPYRNKGSGFQSQPIPGEV 85
""",
        )
        hit = record[4]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|227518463|ref|ZP_03948512.1|")
        self.assertEqual(hit.target.name, "ZP_03948512")
        self.assertEqual(
            hit.target.description,
            "possible histidine kinase [Enterococcus faecalis TX0104] >gi|227074141|gb|EEI12104.1| possible histidine kinase [Enterococcus faecalis TX0104]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=516)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 82.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 36.0945)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.01498)
        self.assertEqual(hsp.annotations["identity"], 25)
        self.assertEqual(hsp.annotations["positive"], 35)
        self.assertEqual(hsp.annotations["gaps"], 5)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[379, 418, 423, 486],
                              [ 29,  68,  68, 131]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 107))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({29: 'PEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQPIPG...MAS'}, length=131)",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({379: 'PILAGFLSGERIKFAEIKTQLAIEIYPEIPPNKRDEDTQNLIAIYRYIHRFLME...EAT'}, length=516)",
        )
        self.assertEqual(hsp.target.id, "gi|227518463|ref|ZP_03948512.1|")
        self.assertEqual(hsp.target.name, "ZP_03948512")
        self.assertEqual(
            hsp.target.description,
            "possible histidine kinase [Enterococcus faecalis TX0104] >gi|227074141|gb|EEI12104.1| possible histidine kinase [Enterococcus faecalis TX0104]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "P +   L+GE         +   E    I     D         YR        QP+P E+I  +   P         A P  Q+   E+E   S  A+ L  A A+",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|227518       379 PILAGFLSGERIKFAEIKTQLAIEIYPEIPPNKRDEDTQNLIAIYRYIHRFLMEQPLPEE
                  0 |.....|.||.............|....|.....|....-----||........||.|.|
lcl|QUERY        29 PEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKP-----YRNKGSGFQSQPIPGE

gi|227518       439 IIETIDYQPGSLTTTYSFAYPKEQLERFEQEFFTSYLARLLENAEAT 486
                 60 .|......|.........|.|..|....|.|...|..|..|..|.|. 107
lcl|QUERY        84 VIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSKGAKKLAQAMAS 131
""",
        )
        hit = record[5]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|78044048|ref|YP_360886.1|")
        self.assertEqual(hit.target.name, "YP_360886")
        self.assertEqual(
            hit.target.description,
            "UDP-N-acetylenolpyruvoylglucosamine reductase [Carboxydothermus hydrogenoformans Z-2901] >gi|90109774|sp|Q3AAE8.1|MURB_CARHZ RecName: Full=UDP-N-acetylenolpyruvoylglucosamine reductase; AltName: Full=UDP-N-acetylmuramate dehydrogenase >gi|77996163|gb|ABB15062.1| UDP-N-acetylenolpyruvoylglucosamine reductase [Carboxydothermus hydrogenoformans Z-2901]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=302)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 78.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 34.5537)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.29147)
        self.assertEqual(hsp.annotations["identity"], 15)
        self.assertEqual(hsp.annotations["positive"], 25)
        self.assertEqual(hsp.annotations["gaps"], 6)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[162, 208, 208, 232],
                              [ 18,  64,  70,  94]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 76))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({18: 'LDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNK...NPE'}, length=131)",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({162: 'LKEFFANECGFKYRSSRFKEEKQWIVKAEFSLNPGDKKEILKKIREFREKRLAS...NPE'}, length=302)",
        )
        self.assertEqual(hsp.target.id, "gi|78044048|ref|YP_360886.1|")
        self.assertEqual(hsp.target.name, "YP_360886")
        self.assertEqual(
            hsp.target.description,
            "UDP-N-acetylenolpyruvoylglucosamine reductase [Carboxydothermus hydrogenoformans Z-2901] >gi|90109774|sp|Q3AAE8.1|MURB_CARHZ RecName: Full=UDP-N-acetylenolpyruvoylglucosamine reductase; AltName: Full=UDP-N-acetylmuramate dehydrogenase >gi|77996163|gb|ABB15062.1| UDP-N-acetylenolpyruvoylglucosamine reductase [Carboxydothermus hydrogenoformans Z-2901]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "L  F   + G+   + +   E+    +A           I+K I +       +     SQP+       V  NPE",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|780440       162 LKEFFANECGFKYRSSRFKEEKQWIVKAEFSLNPGDKKEILKKIRE------FREKRLAS
                  0 |..|.....|..........|......|...........|.|.|..------.......|
lcl|QUERY        18 LDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQS

gi|780440       216 QPLEFPNAGSVFKNPE 232
                 60 ||........|..|||  76
lcl|QUERY        78 QPIPGEVIAQVTSNPE  94
""",
        )
        hit = record[6]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|116750563|ref|YP_847250.1|")
        self.assertEqual(hit.target.name, "YP_847250")
        self.assertEqual(
            hit.target.description,
            "NADH:flavin oxidoreductase/NADH oxidase [Syntrophobacter fumaroxidans MPOB] >gi|116699627|gb|ABK18815.1| NADH:flavin oxidoreductase/NADH oxidase [Syntrophobacter fumaroxidans MPOB]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=645)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 77.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 34.1685)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.67286)
        self.assertEqual(hsp.annotations["identity"], 25)
        self.assertEqual(hsp.annotations["positive"], 35)
        self.assertEqual(hsp.annotations["gaps"], 20)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[417, 430, 430, 449, 449, 461, 461, 497],
                              [  0,  13,  18,  37,  51,  63,  64, 100]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 100))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({0: 'MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEA...AKA'}, length=131)",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({417: 'LAALPPKKGDFGKLVEFYAGELPRLGVDVRLGTAATTELLKSLKADVYVLATGS...GQA'}, length=645)",
        )
        self.assertEqual(hsp.target.id, "gi|116750563|ref|YP_847250.1|")
        self.assertEqual(hsp.target.name, "YP_847250")
        self.assertEqual(
            hsp.target.description,
            "NADH:flavin oxidoreductase/NADH oxidase [Syntrophobacter fumaroxidans MPOB] >gi|116699627|gb|ABK18815.1| NADH:flavin oxidoreductase/NADH oxidase [Syntrophobacter fumaroxidans MPOB]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "+A   PKKGD+       L  F  G+L    V ++L                 A   ++K +  A       GS     PIPG  +  V   PE    +A",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|116750       417 LAALPPKKGDFGK-----LVEFYAGELPRLGVDVRLG--------------TAATTELLK
                  0 .|...|||||...-----|..|..|.|....|...|.--------------..|.....|
lcl|QUERY         0 MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVK

gi|116750       458 SLK-ADVYVLATGSTSSRPPIPGADLPHVFMVPEVLWGQA 497
                 60 ...-|.......||.....||||.....|...||.....| 100
lcl|QUERY        60 MISDAMKPYRNKGSGFQSQPIPGEVIAQVTSNPEYQQAKA 100
""",
        )
        hit = record[7]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|184201116|ref|YP_001855323.1|")
        self.assertEqual(hit.target.name, "YP_001855323")
        self.assertEqual(
            hit.target.description,
            "histidinol dehydrogenase [Kocuria rhizophila DC2201] >gi|183581346|dbj|BAG29817.1| histidinol dehydrogenase [Kocuria rhizophila DC2201]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=507)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 77.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 34.1685)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.67286)
        self.assertEqual(hsp.annotations["identity"], 12)
        self.assertEqual(hsp.annotations["positive"], 24)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 68, 127],
                              [ 30,  89]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 59))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({30: 'EVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQPIPGEVIAQV'}, length=131)",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({68: 'ELAHRFDGVEQQRLRVGPEQIERAVAELAPEVRRALETAIARTRAFAEAQRPRDVEVEV'}, length=507)",
        )
        self.assertEqual(hsp.target.id, "gi|184201116|ref|YP_001855323.1|")
        self.assertEqual(hsp.target.name, "YP_001855323")
        self.assertEqual(
            hsp.target.description,
            "histidinol dehydrogenase [Kocuria rhizophila DC2201] >gi|183581346|dbj|BAG29817.1| histidinol dehydrogenase [Kocuria rhizophila DC2201]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "E+  +  G E    R G E+ + A+  +   +  A++    +   F     P +V  +V",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|184201        68 ELAHRFDGVEQQRLRVGPEQIERAVAELAPEVRRALETAIARTRAFAEAQRPRDVEVEV
                  0 |......|.|....|.|.|....|..........|..........|.....|..|...|
lcl|QUERY        30 EVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQPIPGEVIAQV

gi|184201       127
                 59
lcl|QUERY        89
""",
        )
        hit = record[8]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|46114422|ref|XP_383229.1|")
        self.assertEqual(hit.target.name, "XP_383229")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein FG03053.1 [Gibberella zeae PH-1]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=263)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 74.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 33.0129)
        self.assertAlmostEqual(hsp.annotations["evalue"], 9.81946)
        self.assertEqual(hsp.annotations["identity"], 26)
        self.assertEqual(hsp.annotations["positive"], 40)
        self.assertEqual(hsp.annotations["gaps"], 5)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 85, 118, 118, 126, 126, 151, 151, 169],
                              [ 34,  67,  68,  76,  78, 103, 105, 123]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 89))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({34: 'KLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQPIPGEVIAQ...GAK'}, length=131)",
        )
        self.assertEqual(hsp.query.id, "lcl|QUERY")
        self.assertEqual(hsp.query.description, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({85: 'KLTTLKRRLRYEAKELRKEALRQIVSKLRSWMKLRDFGCGVKPLIESSLPEIAT...GSR'}, length=263)",
        )
        self.assertEqual(hsp.target.id, "gi|46114422|ref|XP_383229.1|")
        self.assertEqual(hsp.target.name, "XP_383229")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein FG03053.1 [Gibberella zeae PH-1]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "KL   +   R    E  KEA+  IV  +   MK  R+ G G   +P+    + ++ +NP Y Q   +L       R IER   + KG++",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|461144        85 KLTTLKRRLRYEAKELRKEALRQIVSKLRSWMK-LRDFGCGV--KPLIESSLPEIATNPH
                  0 ||.......|....|..|||...||......||-.|..|.|.--.|...........||.
lcl|QUERY        34 KLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQPIPGEVIAQVTSNPE

gi|461144       142 YNQTSLYLT--MNHPRRIERVSRVIKGSR 169
                 60 |.|....|.--....|.|||.....||..  89
lcl|QUERY        94 YQQAKAFLASPATQVRNIEREEVLSKGAK 123
""",
        )

    def test_phiblast(self):
        """Parsing BLASTP 2.14.1+ (phiblast.xml)."""
        filename = "phiblast.xml"
        datafile = os.path.join("Blast", filename)
        with open(datafile, "rb") as handle:
            records = Blast.parse(handle)
            self.check_phiblast_records(records)
        with Blast.parse(datafile) as records:
            self.check_phiblast_records(records)
        with open(datafile, "rb") as handle:
            record = Blast.read(handle)
        self.check_phiblast_record(record)
        record = Blast.read(datafile)
        self.check_phiblast_record(record)

    def check_phiblast_records(self, records):
        self.assertEqual(records.program, "blastp")
        self.assertEqual(records.version, "BLASTP 2.14.1+")
        self.assertEqual(
            records.reference,
            'Zheng Zhang, Alejandro A. Sch&auml;ffer, Webb Miller, Thomas L. Madden, David J. Lipman, Eugene V. Koonin, and Stephen F. Altschul (1998), "Protein sequence similarity searches using patterns as seeds", Nucleic Acids Res. 26:3986-3990.',
        )
        self.assertEqual(records.db, "nr")
        self.assertIsInstance(records.query, SeqRecord)
        self.assertEqual(records.query.id, "Query_74414")
        self.assertEqual(records.query.description, "unnamed protein product")
        self.assertEqual(repr(records.query.seq), "Seq(None, length=664)")
        self.assertEqual(len(records.param), 6)
        self.assertEqual(records.param["matrix"], "BLOSUM62")
        self.assertAlmostEqual(records.param["expect"], 0.05)
        self.assertEqual(records.param["gap-open"], 11)
        self.assertEqual(records.param["gap-extend"], 1)
        self.assertEqual(records.param["filter"], "F")
        self.assertEqual(
            records.param["pattern"],
            "[LIVMF]-G-E-x-[GAS]-[LIVM]-x(5,11)-R-[STAQ]-A-x-[LIVMA]-x-[STACV]",
        )
        record = next(records)
        self.assertRaises(StopIteration, next, records)
        self.check_phiblast_record(record)

    def check_phiblast_record(self, record):
        self.assertEqual(record.num, 1)
        self.assertIsInstance(record.query, SeqRecord)
        self.assertEqual(record.query.id, "Query_74414")
        self.assertEqual(record.query.description, "unnamed protein product")
        self.assertEqual(repr(record.query.seq), "Seq(None, length=664)")
        self.assertEqual(len(record.stat), 7)
        self.assertEqual(record.stat["db-num"], 633473216)
        self.assertEqual(record.stat["db-len"], 248084082182)
        self.assertEqual(record.stat["hsp-len"], 0)
        self.assertAlmostEqual(record.stat["eff-space"], 0.0)
        self.assertAlmostEqual(record.stat["kappa"], 0.047)
        self.assertAlmostEqual(record.stat["lambda"], 0.27)
        self.assertAlmostEqual(record.stat["entropy"], 1.0)
        self.assertEqual(len(record), 10)
        hit = record[0]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "ref|NP_001075863.1|")
        self.assertEqual(hit.target.name, "NP_001075863")
        self.assertEqual(
            hit.target.description,
            "cyclic nucleotide-gated olfactory channel [Oryctolagus cuniculus] >emb|CAA42201.1| aorta CNG channel (rACNG) [Oryctolagus cuniculus] >prf||1919268A cyclic nucleotide-gated channel [Oryctolagus cuniculus]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=732)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 3336.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 1290.65)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.0)
        self.assertEqual(hsp.annotations["identity"], 664)
        self.assertEqual(hsp.annotations["positive"], 664)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 68, 732],
                          [  0, 664]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 664))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MTEKSNGVKSSPANNHNNHVPATIKANGKDESRTRSRPQSAADDDTSSELQRLA...EQP')",
        )
        self.assertEqual(hsp.query.id, "Query_74414")
        self.assertEqual(hsp.query.description, "unnamed protein product")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({68: 'MTEKSNGVKSSPANNHNNHVPATIKANGKDESRTRSRPQSAADDDTSSELQRLA...EQP'}, length=732)",
        )
        self.assertEqual(hsp.target.id, "ref|NP_001075863.1|")
        self.assertEqual(hsp.target.name, "NP_001075863")
        self.assertEqual(
            hsp.target.description,
            "cyclic nucleotide-gated olfactory channel [Oryctolagus cuniculus] >emb|CAA42201.1| aorta CNG channel (rACNG) [Oryctolagus cuniculus] >prf||1919268A cyclic nucleotide-gated channel [Oryctolagus cuniculus]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MTEKSNGVKSSPANNHNNHVPATIKANGKDESRTRSRPQSAADDDTSSELQRLAEMDAPQQRRGGFRRIVRLVGVIRQWANRNFREEEARPDSFLERFRGPELQTVTTQQGDGKGDKDGDGKGTKKKFELFVLDPAGDWYYRWLFVIAMPVLYNWCLLVARACFSDLQRGYFLVWLVLDYFSDVVYIADLFIRLRTGFLEQGLLVKDPKKLRDNYIHTLQFKLDVASIIPTDLIYFAVGIHNPELRFNRLLHFARMFEFFDRTETRTSYPNIFRISNLVLYILVIIHWNACIYYAISKSIGFGVDTWVYPNITDPEYGYLAREYIYCLYWSTLTLTTIGETPPPVKDEEYLFVIFDFLIGVLIFATIVGNVGSMISNMNATRAEFQAKIDAVKHYMQFRKVSKEMEAKVIKWFDYLWTNKKTVDEREVLKNLPAKLRAEIAINVHLSTLKKVRIFQDCEAGLLVELVLKLRPQVFSPGDYICRKGDIGKEMYIIKEGKLAVVADDGVTQYALLSAGSCFGEISILNIKGSKMGNRRTANIRSLGYSDLFCLSKDDLMEAVTEYPDAKKVLEERGREILMKEGLLDENEVAASMEVDVQEKLKQLETNMETLYTRFGRLLAEYTGAQQKLKQRITVLEVKMKQNTEDDYLSDGMNSPEPAAAEQP",
        )
        self.assertEqual(
            str(hsp),
            """\
ref|NP_00        68 MTEKSNGVKSSPANNHNNHVPATIKANGKDESRTRSRPQSAADDDTSSELQRLAEMDAPQ
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744         0 MTEKSNGVKSSPANNHNNHVPATIKANGKDESRTRSRPQSAADDDTSSELQRLAEMDAPQ

ref|NP_00       128 QRRGGFRRIVRLVGVIRQWANRNFREEEARPDSFLERFRGPELQTVTTQQGDGKGDKDGD
                 60 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744        60 QRRGGFRRIVRLVGVIRQWANRNFREEEARPDSFLERFRGPELQTVTTQQGDGKGDKDGD

ref|NP_00       188 GKGTKKKFELFVLDPAGDWYYRWLFVIAMPVLYNWCLLVARACFSDLQRGYFLVWLVLDY
                120 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       120 GKGTKKKFELFVLDPAGDWYYRWLFVIAMPVLYNWCLLVARACFSDLQRGYFLVWLVLDY

ref|NP_00       248 FSDVVYIADLFIRLRTGFLEQGLLVKDPKKLRDNYIHTLQFKLDVASIIPTDLIYFAVGI
                180 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       180 FSDVVYIADLFIRLRTGFLEQGLLVKDPKKLRDNYIHTLQFKLDVASIIPTDLIYFAVGI

ref|NP_00       308 HNPELRFNRLLHFARMFEFFDRTETRTSYPNIFRISNLVLYILVIIHWNACIYYAISKSI
                240 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       240 HNPELRFNRLLHFARMFEFFDRTETRTSYPNIFRISNLVLYILVIIHWNACIYYAISKSI

ref|NP_00       368 GFGVDTWVYPNITDPEYGYLAREYIYCLYWSTLTLTTIGETPPPVKDEEYLFVIFDFLIG
                300 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       300 GFGVDTWVYPNITDPEYGYLAREYIYCLYWSTLTLTTIGETPPPVKDEEYLFVIFDFLIG

ref|NP_00       428 VLIFATIVGNVGSMISNMNATRAEFQAKIDAVKHYMQFRKVSKEMEAKVIKWFDYLWTNK
                360 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       360 VLIFATIVGNVGSMISNMNATRAEFQAKIDAVKHYMQFRKVSKEMEAKVIKWFDYLWTNK

ref|NP_00       488 KTVDEREVLKNLPAKLRAEIAINVHLSTLKKVRIFQDCEAGLLVELVLKLRPQVFSPGDY
                420 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       420 KTVDEREVLKNLPAKLRAEIAINVHLSTLKKVRIFQDCEAGLLVELVLKLRPQVFSPGDY

ref|NP_00       548 ICRKGDIGKEMYIIKEGKLAVVADDGVTQYALLSAGSCFGEISILNIKGSKMGNRRTANI
                480 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       480 ICRKGDIGKEMYIIKEGKLAVVADDGVTQYALLSAGSCFGEISILNIKGSKMGNRRTANI

ref|NP_00       608 RSLGYSDLFCLSKDDLMEAVTEYPDAKKVLEERGREILMKEGLLDENEVAASMEVDVQEK
                540 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       540 RSLGYSDLFCLSKDDLMEAVTEYPDAKKVLEERGREILMKEGLLDENEVAASMEVDVQEK

ref|NP_00       668 LKQLETNMETLYTRFGRLLAEYTGAQQKLKQRITVLEVKMKQNTEDDYLSDGMNSPEPAA
                600 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       600 LKQLETNMETLYTRFGRLLAEYTGAQQKLKQRITVLEVKMKQNTEDDYLSDGMNSPEPAA

ref|NP_00       728 AEQP 732
                660 |||| 664
Query_744       660 AEQP 664
""",
        )
        hit = record[1]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "ref|XP_051689802.1|")
        self.assertEqual(hit.target.name, "XP_051689802")
        self.assertEqual(
            hit.target.description,
            "cyclic nucleotide-gated olfactory channel isoform X3 [Oryctolagus cuniculus] >sp|Q28718.1| RecName: Full=Cyclic nucleotide-gated olfactory channel; AltName: Full=Aorta CNG channel; Short=RACNG; AltName: Full=Cyclic nucleotide-gated cation channel 2; AltName: Full=Cyclic nucleotide-gated channel alpha-2; Short=CNG channel alpha-2; Short=CNG-2; Short=CNG2 [Oryctolagus cuniculus]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=664)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 3336.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 1290.65)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.0)
        self.assertEqual(hsp.annotations["identity"], 664)
        self.assertEqual(hsp.annotations["positive"], 664)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0, 664],
                          [  0, 664]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 664))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MTEKSNGVKSSPANNHNNHVPATIKANGKDESRTRSRPQSAADDDTSSELQRLA...EQP')",
        )
        self.assertEqual(hsp.query.id, "Query_74414")
        self.assertEqual(hsp.query.description, "unnamed protein product")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MTEKSNGVKSSPANNHNNHVPATIKANGKDESRTRSRPQSAADDDTSSELQRLA...EQP')",
        )
        self.assertEqual(hsp.target.id, "ref|XP_051689802.1|")
        self.assertEqual(hsp.target.name, "XP_051689802")
        self.assertEqual(
            hsp.target.description,
            "cyclic nucleotide-gated olfactory channel isoform X3 [Oryctolagus cuniculus] >sp|Q28718.1| RecName: Full=Cyclic nucleotide-gated olfactory channel; AltName: Full=Aorta CNG channel; Short=RACNG; AltName: Full=Cyclic nucleotide-gated cation channel 2; AltName: Full=Cyclic nucleotide-gated channel alpha-2; Short=CNG channel alpha-2; Short=CNG-2; Short=CNG2 [Oryctolagus cuniculus]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MTEKSNGVKSSPANNHNNHVPATIKANGKDESRTRSRPQSAADDDTSSELQRLAEMDAPQQRRGGFRRIVRLVGVIRQWANRNFREEEARPDSFLERFRGPELQTVTTQQGDGKGDKDGDGKGTKKKFELFVLDPAGDWYYRWLFVIAMPVLYNWCLLVARACFSDLQRGYFLVWLVLDYFSDVVYIADLFIRLRTGFLEQGLLVKDPKKLRDNYIHTLQFKLDVASIIPTDLIYFAVGIHNPELRFNRLLHFARMFEFFDRTETRTSYPNIFRISNLVLYILVIIHWNACIYYAISKSIGFGVDTWVYPNITDPEYGYLAREYIYCLYWSTLTLTTIGETPPPVKDEEYLFVIFDFLIGVLIFATIVGNVGSMISNMNATRAEFQAKIDAVKHYMQFRKVSKEMEAKVIKWFDYLWTNKKTVDEREVLKNLPAKLRAEIAINVHLSTLKKVRIFQDCEAGLLVELVLKLRPQVFSPGDYICRKGDIGKEMYIIKEGKLAVVADDGVTQYALLSAGSCFGEISILNIKGSKMGNRRTANIRSLGYSDLFCLSKDDLMEAVTEYPDAKKVLEERGREILMKEGLLDENEVAASMEVDVQEKLKQLETNMETLYTRFGRLLAEYTGAQQKLKQRITVLEVKMKQNTEDDYLSDGMNSPEPAAAEQP",
        )
        self.assertEqual(
            str(hsp),
            """\
ref|XP_05         0 MTEKSNGVKSSPANNHNNHVPATIKANGKDESRTRSRPQSAADDDTSSELQRLAEMDAPQ
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744         0 MTEKSNGVKSSPANNHNNHVPATIKANGKDESRTRSRPQSAADDDTSSELQRLAEMDAPQ

ref|XP_05        60 QRRGGFRRIVRLVGVIRQWANRNFREEEARPDSFLERFRGPELQTVTTQQGDGKGDKDGD
                 60 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744        60 QRRGGFRRIVRLVGVIRQWANRNFREEEARPDSFLERFRGPELQTVTTQQGDGKGDKDGD

ref|XP_05       120 GKGTKKKFELFVLDPAGDWYYRWLFVIAMPVLYNWCLLVARACFSDLQRGYFLVWLVLDY
                120 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       120 GKGTKKKFELFVLDPAGDWYYRWLFVIAMPVLYNWCLLVARACFSDLQRGYFLVWLVLDY

ref|XP_05       180 FSDVVYIADLFIRLRTGFLEQGLLVKDPKKLRDNYIHTLQFKLDVASIIPTDLIYFAVGI
                180 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       180 FSDVVYIADLFIRLRTGFLEQGLLVKDPKKLRDNYIHTLQFKLDVASIIPTDLIYFAVGI

ref|XP_05       240 HNPELRFNRLLHFARMFEFFDRTETRTSYPNIFRISNLVLYILVIIHWNACIYYAISKSI
                240 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       240 HNPELRFNRLLHFARMFEFFDRTETRTSYPNIFRISNLVLYILVIIHWNACIYYAISKSI

ref|XP_05       300 GFGVDTWVYPNITDPEYGYLAREYIYCLYWSTLTLTTIGETPPPVKDEEYLFVIFDFLIG
                300 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       300 GFGVDTWVYPNITDPEYGYLAREYIYCLYWSTLTLTTIGETPPPVKDEEYLFVIFDFLIG

ref|XP_05       360 VLIFATIVGNVGSMISNMNATRAEFQAKIDAVKHYMQFRKVSKEMEAKVIKWFDYLWTNK
                360 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       360 VLIFATIVGNVGSMISNMNATRAEFQAKIDAVKHYMQFRKVSKEMEAKVIKWFDYLWTNK

ref|XP_05       420 KTVDEREVLKNLPAKLRAEIAINVHLSTLKKVRIFQDCEAGLLVELVLKLRPQVFSPGDY
                420 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       420 KTVDEREVLKNLPAKLRAEIAINVHLSTLKKVRIFQDCEAGLLVELVLKLRPQVFSPGDY

ref|XP_05       480 ICRKGDIGKEMYIIKEGKLAVVADDGVTQYALLSAGSCFGEISILNIKGSKMGNRRTANI
                480 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       480 ICRKGDIGKEMYIIKEGKLAVVADDGVTQYALLSAGSCFGEISILNIKGSKMGNRRTANI

ref|XP_05       540 RSLGYSDLFCLSKDDLMEAVTEYPDAKKVLEERGREILMKEGLLDENEVAASMEVDVQEK
                540 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       540 RSLGYSDLFCLSKDDLMEAVTEYPDAKKVLEERGREILMKEGLLDENEVAASMEVDVQEK

ref|XP_05       600 LKQLETNMETLYTRFGRLLAEYTGAQQKLKQRITVLEVKMKQNTEDDYLSDGMNSPEPAA
                600 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       600 LKQLETNMETLYTRFGRLLAEYTGAQQKLKQRITVLEVKMKQNTEDDYLSDGMNSPEPAA

ref|XP_05       660 AEQP 664
                660 |||| 664
Query_744       660 AEQP 664
""",
        )
        hit = record[2]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "ref|XP_017206345.1|")
        self.assertEqual(hit.target.name, "XP_017206345")
        self.assertEqual(
            hit.target.description,
            "cyclic nucleotide-gated olfactory channel isoform X2 [Oryctolagus cuniculus]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=677)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 3336.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 1290.65)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.0)
        self.assertEqual(hsp.annotations["identity"], 664)
        self.assertEqual(hsp.annotations["positive"], 664)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 13, 677],
                          [  0, 664]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 664))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MTEKSNGVKSSPANNHNNHVPATIKANGKDESRTRSRPQSAADDDTSSELQRLA...EQP')",
        )
        self.assertEqual(hsp.query.id, "Query_74414")
        self.assertEqual(hsp.query.description, "unnamed protein product")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({13: 'MTEKSNGVKSSPANNHNNHVPATIKANGKDESRTRSRPQSAADDDTSSELQRLA...EQP'}, length=677)",
        )
        self.assertEqual(hsp.target.id, "ref|XP_017206345.1|")
        self.assertEqual(hsp.target.name, "XP_017206345")
        self.assertEqual(
            hsp.target.description,
            "cyclic nucleotide-gated olfactory channel isoform X2 [Oryctolagus cuniculus]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MTEKSNGVKSSPANNHNNHVPATIKANGKDESRTRSRPQSAADDDTSSELQRLAEMDAPQQRRGGFRRIVRLVGVIRQWANRNFREEEARPDSFLERFRGPELQTVTTQQGDGKGDKDGDGKGTKKKFELFVLDPAGDWYYRWLFVIAMPVLYNWCLLVARACFSDLQRGYFLVWLVLDYFSDVVYIADLFIRLRTGFLEQGLLVKDPKKLRDNYIHTLQFKLDVASIIPTDLIYFAVGIHNPELRFNRLLHFARMFEFFDRTETRTSYPNIFRISNLVLYILVIIHWNACIYYAISKSIGFGVDTWVYPNITDPEYGYLAREYIYCLYWSTLTLTTIGETPPPVKDEEYLFVIFDFLIGVLIFATIVGNVGSMISNMNATRAEFQAKIDAVKHYMQFRKVSKEMEAKVIKWFDYLWTNKKTVDEREVLKNLPAKLRAEIAINVHLSTLKKVRIFQDCEAGLLVELVLKLRPQVFSPGDYICRKGDIGKEMYIIKEGKLAVVADDGVTQYALLSAGSCFGEISILNIKGSKMGNRRTANIRSLGYSDLFCLSKDDLMEAVTEYPDAKKVLEERGREILMKEGLLDENEVAASMEVDVQEKLKQLETNMETLYTRFGRLLAEYTGAQQKLKQRITVLEVKMKQNTEDDYLSDGMNSPEPAAAEQP",
        )
        self.assertEqual(
            str(hsp),
            """\
ref|XP_01        13 MTEKSNGVKSSPANNHNNHVPATIKANGKDESRTRSRPQSAADDDTSSELQRLAEMDAPQ
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744         0 MTEKSNGVKSSPANNHNNHVPATIKANGKDESRTRSRPQSAADDDTSSELQRLAEMDAPQ

ref|XP_01        73 QRRGGFRRIVRLVGVIRQWANRNFREEEARPDSFLERFRGPELQTVTTQQGDGKGDKDGD
                 60 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744        60 QRRGGFRRIVRLVGVIRQWANRNFREEEARPDSFLERFRGPELQTVTTQQGDGKGDKDGD

ref|XP_01       133 GKGTKKKFELFVLDPAGDWYYRWLFVIAMPVLYNWCLLVARACFSDLQRGYFLVWLVLDY
                120 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       120 GKGTKKKFELFVLDPAGDWYYRWLFVIAMPVLYNWCLLVARACFSDLQRGYFLVWLVLDY

ref|XP_01       193 FSDVVYIADLFIRLRTGFLEQGLLVKDPKKLRDNYIHTLQFKLDVASIIPTDLIYFAVGI
                180 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       180 FSDVVYIADLFIRLRTGFLEQGLLVKDPKKLRDNYIHTLQFKLDVASIIPTDLIYFAVGI

ref|XP_01       253 HNPELRFNRLLHFARMFEFFDRTETRTSYPNIFRISNLVLYILVIIHWNACIYYAISKSI
                240 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       240 HNPELRFNRLLHFARMFEFFDRTETRTSYPNIFRISNLVLYILVIIHWNACIYYAISKSI

ref|XP_01       313 GFGVDTWVYPNITDPEYGYLAREYIYCLYWSTLTLTTIGETPPPVKDEEYLFVIFDFLIG
                300 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       300 GFGVDTWVYPNITDPEYGYLAREYIYCLYWSTLTLTTIGETPPPVKDEEYLFVIFDFLIG

ref|XP_01       373 VLIFATIVGNVGSMISNMNATRAEFQAKIDAVKHYMQFRKVSKEMEAKVIKWFDYLWTNK
                360 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       360 VLIFATIVGNVGSMISNMNATRAEFQAKIDAVKHYMQFRKVSKEMEAKVIKWFDYLWTNK

ref|XP_01       433 KTVDEREVLKNLPAKLRAEIAINVHLSTLKKVRIFQDCEAGLLVELVLKLRPQVFSPGDY
                420 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       420 KTVDEREVLKNLPAKLRAEIAINVHLSTLKKVRIFQDCEAGLLVELVLKLRPQVFSPGDY

ref|XP_01       493 ICRKGDIGKEMYIIKEGKLAVVADDGVTQYALLSAGSCFGEISILNIKGSKMGNRRTANI
                480 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       480 ICRKGDIGKEMYIIKEGKLAVVADDGVTQYALLSAGSCFGEISILNIKGSKMGNRRTANI

ref|XP_01       553 RSLGYSDLFCLSKDDLMEAVTEYPDAKKVLEERGREILMKEGLLDENEVAASMEVDVQEK
                540 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       540 RSLGYSDLFCLSKDDLMEAVTEYPDAKKVLEERGREILMKEGLLDENEVAASMEVDVQEK

ref|XP_01       613 LKQLETNMETLYTRFGRLLAEYTGAQQKLKQRITVLEVKMKQNTEDDYLSDGMNSPEPAA
                600 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       600 LKQLETNMETLYTRFGRLLAEYTGAQQKLKQRITVLEVKMKQNTEDDYLSDGMNSPEPAA

ref|XP_01       673 AEQP 677
                660 |||| 664
Query_744       660 AEQP 664
""",
        )
        hit = record[3]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "ref|XP_051689801.1|")
        self.assertEqual(hit.target.name, "XP_051689801")
        self.assertEqual(
            hit.target.description,
            "cyclic nucleotide-gated olfactory channel isoform X1 [Oryctolagus cuniculus]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=687)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 3336.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 1290.65)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.0)
        self.assertEqual(hsp.annotations["identity"], 664)
        self.assertEqual(hsp.annotations["positive"], 664)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 23, 687],
                          [  0, 664]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 664))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MTEKSNGVKSSPANNHNNHVPATIKANGKDESRTRSRPQSAADDDTSSELQRLA...EQP')",
        )
        self.assertEqual(hsp.query.id, "Query_74414")
        self.assertEqual(hsp.query.description, "unnamed protein product")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({23: 'MTEKSNGVKSSPANNHNNHVPATIKANGKDESRTRSRPQSAADDDTSSELQRLA...EQP'}, length=687)",
        )
        self.assertEqual(hsp.target.id, "ref|XP_051689801.1|")
        self.assertEqual(hsp.target.name, "XP_051689801")
        self.assertEqual(
            hsp.target.description,
            "cyclic nucleotide-gated olfactory channel isoform X1 [Oryctolagus cuniculus]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MTEKSNGVKSSPANNHNNHVPATIKANGKDESRTRSRPQSAADDDTSSELQRLAEMDAPQQRRGGFRRIVRLVGVIRQWANRNFREEEARPDSFLERFRGPELQTVTTQQGDGKGDKDGDGKGTKKKFELFVLDPAGDWYYRWLFVIAMPVLYNWCLLVARACFSDLQRGYFLVWLVLDYFSDVVYIADLFIRLRTGFLEQGLLVKDPKKLRDNYIHTLQFKLDVASIIPTDLIYFAVGIHNPELRFNRLLHFARMFEFFDRTETRTSYPNIFRISNLVLYILVIIHWNACIYYAISKSIGFGVDTWVYPNITDPEYGYLAREYIYCLYWSTLTLTTIGETPPPVKDEEYLFVIFDFLIGVLIFATIVGNVGSMISNMNATRAEFQAKIDAVKHYMQFRKVSKEMEAKVIKWFDYLWTNKKTVDEREVLKNLPAKLRAEIAINVHLSTLKKVRIFQDCEAGLLVELVLKLRPQVFSPGDYICRKGDIGKEMYIIKEGKLAVVADDGVTQYALLSAGSCFGEISILNIKGSKMGNRRTANIRSLGYSDLFCLSKDDLMEAVTEYPDAKKVLEERGREILMKEGLLDENEVAASMEVDVQEKLKQLETNMETLYTRFGRLLAEYTGAQQKLKQRITVLEVKMKQNTEDDYLSDGMNSPEPAAAEQP",
        )
        self.assertEqual(
            str(hsp),
            """\
ref|XP_05        23 MTEKSNGVKSSPANNHNNHVPATIKANGKDESRTRSRPQSAADDDTSSELQRLAEMDAPQ
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744         0 MTEKSNGVKSSPANNHNNHVPATIKANGKDESRTRSRPQSAADDDTSSELQRLAEMDAPQ

ref|XP_05        83 QRRGGFRRIVRLVGVIRQWANRNFREEEARPDSFLERFRGPELQTVTTQQGDGKGDKDGD
                 60 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744        60 QRRGGFRRIVRLVGVIRQWANRNFREEEARPDSFLERFRGPELQTVTTQQGDGKGDKDGD

ref|XP_05       143 GKGTKKKFELFVLDPAGDWYYRWLFVIAMPVLYNWCLLVARACFSDLQRGYFLVWLVLDY
                120 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       120 GKGTKKKFELFVLDPAGDWYYRWLFVIAMPVLYNWCLLVARACFSDLQRGYFLVWLVLDY

ref|XP_05       203 FSDVVYIADLFIRLRTGFLEQGLLVKDPKKLRDNYIHTLQFKLDVASIIPTDLIYFAVGI
                180 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       180 FSDVVYIADLFIRLRTGFLEQGLLVKDPKKLRDNYIHTLQFKLDVASIIPTDLIYFAVGI

ref|XP_05       263 HNPELRFNRLLHFARMFEFFDRTETRTSYPNIFRISNLVLYILVIIHWNACIYYAISKSI
                240 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       240 HNPELRFNRLLHFARMFEFFDRTETRTSYPNIFRISNLVLYILVIIHWNACIYYAISKSI

ref|XP_05       323 GFGVDTWVYPNITDPEYGYLAREYIYCLYWSTLTLTTIGETPPPVKDEEYLFVIFDFLIG
                300 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       300 GFGVDTWVYPNITDPEYGYLAREYIYCLYWSTLTLTTIGETPPPVKDEEYLFVIFDFLIG

ref|XP_05       383 VLIFATIVGNVGSMISNMNATRAEFQAKIDAVKHYMQFRKVSKEMEAKVIKWFDYLWTNK
                360 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       360 VLIFATIVGNVGSMISNMNATRAEFQAKIDAVKHYMQFRKVSKEMEAKVIKWFDYLWTNK

ref|XP_05       443 KTVDEREVLKNLPAKLRAEIAINVHLSTLKKVRIFQDCEAGLLVELVLKLRPQVFSPGDY
                420 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       420 KTVDEREVLKNLPAKLRAEIAINVHLSTLKKVRIFQDCEAGLLVELVLKLRPQVFSPGDY

ref|XP_05       503 ICRKGDIGKEMYIIKEGKLAVVADDGVTQYALLSAGSCFGEISILNIKGSKMGNRRTANI
                480 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       480 ICRKGDIGKEMYIIKEGKLAVVADDGVTQYALLSAGSCFGEISILNIKGSKMGNRRTANI

ref|XP_05       563 RSLGYSDLFCLSKDDLMEAVTEYPDAKKVLEERGREILMKEGLLDENEVAASMEVDVQEK
                540 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       540 RSLGYSDLFCLSKDDLMEAVTEYPDAKKVLEERGREILMKEGLLDENEVAASMEVDVQEK

ref|XP_05       623 LKQLETNMETLYTRFGRLLAEYTGAQQKLKQRITVLEVKMKQNTEDDYLSDGMNSPEPAA
                600 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       600 LKQLETNMETLYTRFGRLLAEYTGAQQKLKQRITVLEVKMKQNTEDDYLSDGMNSPEPAA

ref|XP_05       683 AEQP 687
                660 |||| 664
Query_744       660 AEQP 664
""",
        )
        hit = record[4]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "ref|XP_004407164.1|")
        self.assertEqual(hit.target.name, "XP_004407164")
        self.assertEqual(
            hit.target.description,
            "PREDICTED: cyclic nucleotide-gated olfactory channel [Odobenus rosmarus divergens]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=664)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 3231.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 1249.79)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.0)
        self.assertEqual(hsp.annotations["identity"], 639)
        self.assertEqual(hsp.annotations["positive"], 652)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0, 664],
                          [  0, 664]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 664))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MTEKSNGVKSSPANNHNNHVPATIKANGKDESRTRSRPQSAADDDTSSELQRLA...EQP')",
        )
        self.assertEqual(hsp.query.id, "Query_74414")
        self.assertEqual(hsp.query.description, "unnamed protein product")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MTEKSNGVKSSPANNHNHHTPPAIKANGKDDHRTNSRPQSAADDDTSSELQRLA...DEP')",
        )
        self.assertEqual(hsp.target.id, "ref|XP_004407164.1|")
        self.assertEqual(hsp.target.name, "XP_004407164")
        self.assertEqual(
            hsp.target.description,
            "PREDICTED: cyclic nucleotide-gated olfactory channel [Odobenus rosmarus divergens]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MTEKSNGVKSSPANNHN+H P  IKANGKD+ RT SRPQSAADDDTSSELQRLAEMDAPQQ RGGFRRIVRLVG+IR+WAN+NFREEE RPDSFLERFRGPELQTVTTQQGDGKGDKDG+GKGTKKKFELFVLDPAGDWYYRWLFVIAMPVLYNWCLLVARACFSDLQRGY+LVWLVLDYFSDVVYI DLFIRLRTGFLEQGLLVKDPKKLRDNYIHTLQFKLDVASIIPTDLIYFAVGIH+PELRFNRLLHFARMFEFFDRTETRTSYPNIFRISNLVLYILVIIHWNACIYYAISKSIGFGVDTWVYPNITDPEYGYLAREYIYCLYWSTLTLTTIGETPPPVKDEEYLFVIFDFLIGVLIFATIVGNVGSMISNMNATRAEFQAKIDAVKHYMQFRKVSKEMEAKVIKWFDYLWTNKK+VDEREVLKNLPAKLRAEIAINVHLSTLKKVRIFQDCEAGLLVELVLKLRPQVFSPGDYICRKGDIGKEMYIIKEGKLAVVADDGVTQYALLSAGSCFGEISILNIKGSKMGNRRTANIRSLGYSDLFCLSKDDLMEAVTEYPDAKKVLEERGREILMKEGLLDENEVAASMEVDVQEKL+QLETNMETLYTRFGRLLAEYTGAQQKLKQRITVLE KMKQN  DDYLSDG+NSPEP AA++P",
        )
        self.assertEqual(
            str(hsp),
            """\
ref|XP_00         0 MTEKSNGVKSSPANNHNHHTPPAIKANGKDDHRTNSRPQSAADDDTSSELQRLAEMDAPQ
                  0 |||||||||||||||||.|.|..|||||||..||.|||||||||||||||||||||||||
Query_744         0 MTEKSNGVKSSPANNHNNHVPATIKANGKDESRTRSRPQSAADDDTSSELQRLAEMDAPQ

ref|XP_00        60 QGRGGFRRIVRLVGIIREWANKNFREEEPRPDSFLERFRGPELQTVTTQQGDGKGDKDGE
                 60 |.||||||||||||.||.|||.||||||.||||||||||||||||||||||||||||||.
Query_744        60 QRRGGFRRIVRLVGVIRQWANRNFREEEARPDSFLERFRGPELQTVTTQQGDGKGDKDGD

ref|XP_00       120 GKGTKKKFELFVLDPAGDWYYRWLFVIAMPVLYNWCLLVARACFSDLQRGYYLVWLVLDY
                120 |||||||||||||||||||||||||||||||||||||||||||||||||||.||||||||
Query_744       120 GKGTKKKFELFVLDPAGDWYYRWLFVIAMPVLYNWCLLVARACFSDLQRGYFLVWLVLDY

ref|XP_00       180 FSDVVYITDLFIRLRTGFLEQGLLVKDPKKLRDNYIHTLQFKLDVASIIPTDLIYFAVGI
                180 |||||||.||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       180 FSDVVYIADLFIRLRTGFLEQGLLVKDPKKLRDNYIHTLQFKLDVASIIPTDLIYFAVGI

ref|XP_00       240 HSPELRFNRLLHFARMFEFFDRTETRTSYPNIFRISNLVLYILVIIHWNACIYYAISKSI
                240 |.||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       240 HNPELRFNRLLHFARMFEFFDRTETRTSYPNIFRISNLVLYILVIIHWNACIYYAISKSI

ref|XP_00       300 GFGVDTWVYPNITDPEYGYLAREYIYCLYWSTLTLTTIGETPPPVKDEEYLFVIFDFLIG
                300 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       300 GFGVDTWVYPNITDPEYGYLAREYIYCLYWSTLTLTTIGETPPPVKDEEYLFVIFDFLIG

ref|XP_00       360 VLIFATIVGNVGSMISNMNATRAEFQAKIDAVKHYMQFRKVSKEMEAKVIKWFDYLWTNK
                360 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       360 VLIFATIVGNVGSMISNMNATRAEFQAKIDAVKHYMQFRKVSKEMEAKVIKWFDYLWTNK

ref|XP_00       420 KSVDEREVLKNLPAKLRAEIAINVHLSTLKKVRIFQDCEAGLLVELVLKLRPQVFSPGDY
                420 |.||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       420 KTVDEREVLKNLPAKLRAEIAINVHLSTLKKVRIFQDCEAGLLVELVLKLRPQVFSPGDY

ref|XP_00       480 ICRKGDIGKEMYIIKEGKLAVVADDGVTQYALLSAGSCFGEISILNIKGSKMGNRRTANI
                480 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       480 ICRKGDIGKEMYIIKEGKLAVVADDGVTQYALLSAGSCFGEISILNIKGSKMGNRRTANI

ref|XP_00       540 RSLGYSDLFCLSKDDLMEAVTEYPDAKKVLEERGREILMKEGLLDENEVAASMEVDVQEK
                540 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       540 RSLGYSDLFCLSKDDLMEAVTEYPDAKKVLEERGREILMKEGLLDENEVAASMEVDVQEK

ref|XP_00       600 LEQLETNMETLYTRFGRLLAEYTGAQQKLKQRITVLETKMKQNNMDDYLSDGVNSPEPTA
                600 |.|||||||||||||||||||||||||||||||||||.|||||..|||||||.|||||.|
Query_744       600 LKQLETNMETLYTRFGRLLAEYTGAQQKLKQRITVLEVKMKQNTEDDYLSDGMNSPEPAA

ref|XP_00       660 ADEP 664
                660 |..| 664
Query_744       660 AEQP 664
""",
        )
        hit = record[5]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "ref|XP_008688471.1|")
        self.assertEqual(hit.target.name, "XP_008688471")
        self.assertEqual(
            hit.target.description,
            "cyclic nucleotide-gated olfactory channel [Ursus maritimus] >ref|XP_026343324.1| cyclic nucleotide-gated olfactory channel [Ursus arctos]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=664)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 3228.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 1248.63)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.0)
        self.assertEqual(hsp.annotations["identity"], 638)
        self.assertEqual(hsp.annotations["positive"], 652)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0, 664],
                          [  0, 664]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 664))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MTEKSNGVKSSPANNHNNHVPATIKANGKDESRTRSRPQSAADDDTSSELQRLA...EQP')",
        )
        self.assertEqual(hsp.query.id, "Query_74414")
        self.assertEqual(hsp.query.description, "unnamed protein product")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MTEKSNGVKSSPANNHNHHAPPAIKANGKDDHRSSSRPQSAVDDDTSSELQRLA...DEP')",
        )
        self.assertEqual(hsp.target.id, "ref|XP_008688471.1|")
        self.assertEqual(hsp.target.name, "XP_008688471")
        self.assertEqual(
            hsp.target.description,
            "cyclic nucleotide-gated olfactory channel [Ursus maritimus] >ref|XP_026343324.1| cyclic nucleotide-gated olfactory channel [Ursus arctos]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MTEKSNGVKSSPANNHN+H P  IKANGKD+ R+ SRPQSA DDDTSSELQRLAEMDAPQ+ RGGFRRIVRLVG+IR WAN+NFREEE RPDSFLERFRGPELQTVTTQQGDGKGDKDG+GKGTKKKFELFVLDPAGDWYYRWLFVIAMPVLYNWCLLVARACFSDLQ+GY+LVWLVLDYFSDVVYI DLFIRLRTGFLEQGLLVKDPKKLRDNYIHTLQFKLDVASIIPTDLIYFAVGIH+PELRFNRLLHFARMFEFFDRTETRTSYPNIFRISNLVLYILVIIHWNACIYYAISKSIGFGVDTWVYPNITDPEYGYLAREYIYCLYWSTLTLTTIGETPPPVKDEEYLFVIFDFLIGVLIFATIVGNVGSMISNMNATRAEFQAKIDAVKHYMQFRKVSKEMEAKVIKWFDYLWTNKK+VDEREVLKNLPAKLRAEIAINVHLSTLKKVRIFQDCEAGLLVELVLKLRPQVFSPGDYICRKGDIGKEMYIIKEGKLAVVADDGVTQYALLSAGSCFGEISILNIKGSKMGNRRTANIRSLGYSDLFCLSKDDLMEAVTEYPDAKKVLEERGREILMKEGLLDENEVAASMEVDVQEKL+QLETNMETLYTRFGRLLAEYTGAQQKLKQRITVLE KMKQN EDDYLSDGMNSPEPAAA++P",
        )
        self.assertEqual(
            str(hsp),
            """\
ref|XP_00         0 MTEKSNGVKSSPANNHNHHAPPAIKANGKDDHRSSSRPQSAVDDDTSSELQRLAEMDAPQ
                  0 |||||||||||||||||.|.|..|||||||..|..||||||.||||||||||||||||||
Query_744         0 MTEKSNGVKSSPANNHNNHVPATIKANGKDESRTRSRPQSAADDDTSSELQRLAEMDAPQ

ref|XP_00        60 RGRGGFRRIVRLVGIIRDWANKNFREEEPRPDSFLERFRGPELQTVTTQQGDGKGDKDGE
                 60 ..||||||||||||.||.|||.||||||.||||||||||||||||||||||||||||||.
Query_744        60 QRRGGFRRIVRLVGVIRQWANRNFREEEARPDSFLERFRGPELQTVTTQQGDGKGDKDGD

ref|XP_00       120 GKGTKKKFELFVLDPAGDWYYRWLFVIAMPVLYNWCLLVARACFSDLQKGYYLVWLVLDY
                120 ||||||||||||||||||||||||||||||||||||||||||||||||.||.||||||||
Query_744       120 GKGTKKKFELFVLDPAGDWYYRWLFVIAMPVLYNWCLLVARACFSDLQRGYFLVWLVLDY

ref|XP_00       180 FSDVVYITDLFIRLRTGFLEQGLLVKDPKKLRDNYIHTLQFKLDVASIIPTDLIYFAVGI
                180 |||||||.||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       180 FSDVVYIADLFIRLRTGFLEQGLLVKDPKKLRDNYIHTLQFKLDVASIIPTDLIYFAVGI

ref|XP_00       240 HSPELRFNRLLHFARMFEFFDRTETRTSYPNIFRISNLVLYILVIIHWNACIYYAISKSI
                240 |.||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       240 HNPELRFNRLLHFARMFEFFDRTETRTSYPNIFRISNLVLYILVIIHWNACIYYAISKSI

ref|XP_00       300 GFGVDTWVYPNITDPEYGYLAREYIYCLYWSTLTLTTIGETPPPVKDEEYLFVIFDFLIG
                300 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       300 GFGVDTWVYPNITDPEYGYLAREYIYCLYWSTLTLTTIGETPPPVKDEEYLFVIFDFLIG

ref|XP_00       360 VLIFATIVGNVGSMISNMNATRAEFQAKIDAVKHYMQFRKVSKEMEAKVIKWFDYLWTNK
                360 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       360 VLIFATIVGNVGSMISNMNATRAEFQAKIDAVKHYMQFRKVSKEMEAKVIKWFDYLWTNK

ref|XP_00       420 KSVDEREVLKNLPAKLRAEIAINVHLSTLKKVRIFQDCEAGLLVELVLKLRPQVFSPGDY
                420 |.||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       420 KTVDEREVLKNLPAKLRAEIAINVHLSTLKKVRIFQDCEAGLLVELVLKLRPQVFSPGDY

ref|XP_00       480 ICRKGDIGKEMYIIKEGKLAVVADDGVTQYALLSAGSCFGEISILNIKGSKMGNRRTANI
                480 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       480 ICRKGDIGKEMYIIKEGKLAVVADDGVTQYALLSAGSCFGEISILNIKGSKMGNRRTANI

ref|XP_00       540 RSLGYSDLFCLSKDDLMEAVTEYPDAKKVLEERGREILMKEGLLDENEVAASMEVDVQEK
                540 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       540 RSLGYSDLFCLSKDDLMEAVTEYPDAKKVLEERGREILMKEGLLDENEVAASMEVDVQEK

ref|XP_00       600 LEQLETNMETLYTRFGRLLAEYTGAQQKLKQRITVLETKMKQNNEDDYLSDGMNSPEPAA
                600 |.|||||||||||||||||||||||||||||||||||.|||||.||||||||||||||||
Query_744       600 LKQLETNMETLYTRFGRLLAEYTGAQQKLKQRITVLEVKMKQNTEDDYLSDGMNSPEPAA

ref|XP_00       660 ADEP 664
                660 |..| 664
Query_744       660 AEQP 664
""",
        )
        hit = record[6]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "ref|XP_011229794.1|")
        self.assertEqual(hit.target.name, "XP_011229794")
        self.assertEqual(
            hit.target.description,
            "cyclic nucleotide-gated olfactory channel [Ailuropoda melanoleuca] >gb|EFB14215.1| hypothetical protein PANDA_013994, partial [Ailuropoda melanoleuca]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=664)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 3227.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 1248.24)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.0)
        self.assertEqual(hsp.annotations["identity"], 638)
        self.assertEqual(hsp.annotations["positive"], 652)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0, 664],
                          [  0, 664]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 664))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MTEKSNGVKSSPANNHNNHVPATIKANGKDESRTRSRPQSAADDDTSSELQRLA...EQP')",
        )
        self.assertEqual(hsp.query.id, "Query_74414")
        self.assertEqual(hsp.query.description, "unnamed protein product")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MTEKSNGVKSSPANNHNHHAPPAIKANGKDDHRSSSRPQSAVDDDTSSELQRLA...DEP')",
        )
        self.assertEqual(hsp.target.id, "ref|XP_011229794.1|")
        self.assertEqual(hsp.target.name, "XP_011229794")
        self.assertEqual(
            hsp.target.description,
            "cyclic nucleotide-gated olfactory channel [Ailuropoda melanoleuca] >gb|EFB14215.1| hypothetical protein PANDA_013994, partial [Ailuropoda melanoleuca]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MTEKSNGVKSSPANNHN+H P  IKANGKD+ R+ SRPQSA DDDTSSELQRLAEMDAPQ+ RGGFRRIVRLVG+IR WAN+NFREEE RPDSFLERFRGPELQTVTTQQGDGKGDKDG+GKGTKKKFELFVLDPAGDWYYRWLFVIAMPVLYNWCLLVARACFSDLQ+GY+LVWLVLDYFSDVVYI DLFIRLRTGFLEQGLLVKDPKKLRDNYIHTLQFKLDVASIIPTDLIYFAVGIH+PELRFNRLLHFARMFEFFDRTETRTSYPNIFRISNLVLYILVIIHWNACIYYAISKSIGFGVDTWVYPNITDPEYGYLAREYIYCLYWSTLTLTTIGETPPPVKDEEYLFVIFDFLIGVLIFATIVGNVGSMISNMNATRAEFQAKIDAVKHYMQFRKVSKEMEAKVIKWFDYLWTNKK+VDEREVLKNLPAKLRAEIAINVHLSTLKKVRIFQDCEAGLLVELVLKLRPQVFSPGDYICRKGDIGKEMYIIKEGKLAVVADDGVTQYALLSAGSCFGEISILNIKGSKMGNRRTANIRSLGYSDLFCLSKDDLMEAVTEYPDAKKVLEERGREILMKEGLLDENEVAASMEVDVQEKL+QLETNMETLYTRFGRLLAEYTGAQQKLKQRITVLE KMKQN EDDYLSDGMNSPEPAAA++P",
        )
        self.assertEqual(
            str(hsp),
            """\
ref|XP_01         0 MTEKSNGVKSSPANNHNHHAPPAIKANGKDDHRSSSRPQSAVDDDTSSELQRLAEMDAPQ
                  0 |||||||||||||||||.|.|..|||||||..|..||||||.||||||||||||||||||
Query_744         0 MTEKSNGVKSSPANNHNNHVPATIKANGKDESRTRSRPQSAADDDTSSELQRLAEMDAPQ

ref|XP_01        60 RGRGGFRRIVRLVGIIRDWANKNFREEEPRPDSFLERFRGPELQTVTTQQGDGKGDKDGE
                 60 ..||||||||||||.||.|||.||||||.||||||||||||||||||||||||||||||.
Query_744        60 QRRGGFRRIVRLVGVIRQWANRNFREEEARPDSFLERFRGPELQTVTTQQGDGKGDKDGD

ref|XP_01       120 GKGTKKKFELFVLDPAGDWYYRWLFVIAMPVLYNWCLLVARACFSDLQKGYYLVWLVLDY
                120 ||||||||||||||||||||||||||||||||||||||||||||||||.||.||||||||
Query_744       120 GKGTKKKFELFVLDPAGDWYYRWLFVIAMPVLYNWCLLVARACFSDLQRGYFLVWLVLDY

ref|XP_01       180 FSDVVYIIDLFIRLRTGFLEQGLLVKDPKKLRDNYIHTLQFKLDVASIIPTDLIYFAVGI
                180 |||||||.||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       180 FSDVVYIADLFIRLRTGFLEQGLLVKDPKKLRDNYIHTLQFKLDVASIIPTDLIYFAVGI

ref|XP_01       240 HSPELRFNRLLHFARMFEFFDRTETRTSYPNIFRISNLVLYILVIIHWNACIYYAISKSI
                240 |.||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       240 HNPELRFNRLLHFARMFEFFDRTETRTSYPNIFRISNLVLYILVIIHWNACIYYAISKSI

ref|XP_01       300 GFGVDTWVYPNITDPEYGYLAREYIYCLYWSTLTLTTIGETPPPVKDEEYLFVIFDFLIG
                300 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       300 GFGVDTWVYPNITDPEYGYLAREYIYCLYWSTLTLTTIGETPPPVKDEEYLFVIFDFLIG

ref|XP_01       360 VLIFATIVGNVGSMISNMNATRAEFQAKIDAVKHYMQFRKVSKEMEAKVIKWFDYLWTNK
                360 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       360 VLIFATIVGNVGSMISNMNATRAEFQAKIDAVKHYMQFRKVSKEMEAKVIKWFDYLWTNK

ref|XP_01       420 KSVDEREVLKNLPAKLRAEIAINVHLSTLKKVRIFQDCEAGLLVELVLKLRPQVFSPGDY
                420 |.||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       420 KTVDEREVLKNLPAKLRAEIAINVHLSTLKKVRIFQDCEAGLLVELVLKLRPQVFSPGDY

ref|XP_01       480 ICRKGDIGKEMYIIKEGKLAVVADDGVTQYALLSAGSCFGEISILNIKGSKMGNRRTANI
                480 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       480 ICRKGDIGKEMYIIKEGKLAVVADDGVTQYALLSAGSCFGEISILNIKGSKMGNRRTANI

ref|XP_01       540 RSLGYSDLFCLSKDDLMEAVTEYPDAKKVLEERGREILMKEGLLDENEVAASMEVDVQEK
                540 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       540 RSLGYSDLFCLSKDDLMEAVTEYPDAKKVLEERGREILMKEGLLDENEVAASMEVDVQEK

ref|XP_01       600 LEQLETNMETLYTRFGRLLAEYTGAQQKLKQRITVLETKMKQNNEDDYLSDGMNSPEPAA
                600 |.|||||||||||||||||||||||||||||||||||.|||||.||||||||||||||||
Query_744       600 LKQLETNMETLYTRFGRLLAEYTGAQQKLKQRITVLEVKMKQNTEDDYLSDGMNSPEPAA

ref|XP_01       660 ADEP 664
                660 |..| 664
Query_744       660 AEQP 664
""",
        )
        hit = record[7]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "ref|XP_045646452.1|")
        self.assertEqual(hit.target.name, "XP_045646452")
        self.assertEqual(
            hit.target.description,
            "cyclic nucleotide-gated olfactory channel [Ursus americanus]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=664)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 3223.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 1246.68)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.0)
        self.assertEqual(hsp.annotations["identity"], 637)
        self.assertEqual(hsp.annotations["positive"], 651)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0, 664],
                          [  0, 664]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 664))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MTEKSNGVKSSPANNHNNHVPATIKANGKDESRTRSRPQSAADDDTSSELQRLA...EQP')",
        )
        self.assertEqual(hsp.query.id, "Query_74414")
        self.assertEqual(hsp.query.description, "unnamed protein product")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MTEKSNGVKCSPANNHNHHAPPAIKANGKDDHRSSSRPQSAVDDDTSSELQRLA...DEP')",
        )
        self.assertEqual(hsp.target.id, "ref|XP_045646452.1|")
        self.assertEqual(hsp.target.name, "XP_045646452")
        self.assertEqual(
            hsp.target.description,
            "cyclic nucleotide-gated olfactory channel [Ursus americanus]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MTEKSNGVK SPANNHN+H P  IKANGKD+ R+ SRPQSA DDDTSSELQRLAEMDAPQ+ RGGFRRIVRLVG+IR WAN+NFREEE RPDSFLERFRGPELQTVTTQQGDGKGDKDG+GKGTKKKFELFVLDPAGDWYYRWLFVIAMPVLYNWCLLVARACFSDLQ+GY+LVWLVLDYFSDVVYI DLFIRLRTGFLEQGLLVKDPKKLRDNYIHTLQFKLDVASIIPTDLIYFAVGIH+PELRFNRLLHFARMFEFFDRTETRTSYPNIFRISNLVLYILVIIHWNACIYYAISKSIGFGVDTWVYPNITDPEYGYLAREYIYCLYWSTLTLTTIGETPPPVKDEEYLFVIFDFLIGVLIFATIVGNVGSMISNMNATRAEFQAKIDAVKHYMQFRKVSKEMEAKVIKWFDYLWTNKK+VDEREVLKNLPAKLRAEIAINVHLSTLKKVRIFQDCEAGLLVELVLKLRPQVFSPGDYICRKGDIGKEMYIIKEGKLAVVADDGVTQYALLSAGSCFGEISILNIKGSKMGNRRTANIRSLGYSDLFCLSKDDLMEAVTEYPDAKKVLEERGREILMKEGLLDENEVAASMEVDVQEKL+QLETNMETLYTRFGRLLAEYTGAQQKLKQRITVLE KMKQN EDDYLSDGMNSPEPAAA++P",
        )
        self.assertEqual(
            str(hsp),
            """\
ref|XP_04         0 MTEKSNGVKCSPANNHNHHAPPAIKANGKDDHRSSSRPQSAVDDDTSSELQRLAEMDAPQ
                  0 |||||||||.|||||||.|.|..|||||||..|..||||||.||||||||||||||||||
Query_744         0 MTEKSNGVKSSPANNHNNHVPATIKANGKDESRTRSRPQSAADDDTSSELQRLAEMDAPQ

ref|XP_04        60 RGRGGFRRIVRLVGIIRDWANKNFREEEPRPDSFLERFRGPELQTVTTQQGDGKGDKDGE
                 60 ..||||||||||||.||.|||.||||||.||||||||||||||||||||||||||||||.
Query_744        60 QRRGGFRRIVRLVGVIRQWANRNFREEEARPDSFLERFRGPELQTVTTQQGDGKGDKDGD

ref|XP_04       120 GKGTKKKFELFVLDPAGDWYYRWLFVIAMPVLYNWCLLVARACFSDLQKGYYLVWLVLDY
                120 ||||||||||||||||||||||||||||||||||||||||||||||||.||.||||||||
Query_744       120 GKGTKKKFELFVLDPAGDWYYRWLFVIAMPVLYNWCLLVARACFSDLQRGYFLVWLVLDY

ref|XP_04       180 FSDVVYITDLFIRLRTGFLEQGLLVKDPKKLRDNYIHTLQFKLDVASIIPTDLIYFAVGI
                180 |||||||.||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       180 FSDVVYIADLFIRLRTGFLEQGLLVKDPKKLRDNYIHTLQFKLDVASIIPTDLIYFAVGI

ref|XP_04       240 HSPELRFNRLLHFARMFEFFDRTETRTSYPNIFRISNLVLYILVIIHWNACIYYAISKSI
                240 |.||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       240 HNPELRFNRLLHFARMFEFFDRTETRTSYPNIFRISNLVLYILVIIHWNACIYYAISKSI

ref|XP_04       300 GFGVDTWVYPNITDPEYGYLAREYIYCLYWSTLTLTTIGETPPPVKDEEYLFVIFDFLIG
                300 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       300 GFGVDTWVYPNITDPEYGYLAREYIYCLYWSTLTLTTIGETPPPVKDEEYLFVIFDFLIG

ref|XP_04       360 VLIFATIVGNVGSMISNMNATRAEFQAKIDAVKHYMQFRKVSKEMEAKVIKWFDYLWTNK
                360 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       360 VLIFATIVGNVGSMISNMNATRAEFQAKIDAVKHYMQFRKVSKEMEAKVIKWFDYLWTNK

ref|XP_04       420 KSVDEREVLKNLPAKLRAEIAINVHLSTLKKVRIFQDCEAGLLVELVLKLRPQVFSPGDY
                420 |.||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       420 KTVDEREVLKNLPAKLRAEIAINVHLSTLKKVRIFQDCEAGLLVELVLKLRPQVFSPGDY

ref|XP_04       480 ICRKGDIGKEMYIIKEGKLAVVADDGVTQYALLSAGSCFGEISILNIKGSKMGNRRTANI
                480 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       480 ICRKGDIGKEMYIIKEGKLAVVADDGVTQYALLSAGSCFGEISILNIKGSKMGNRRTANI

ref|XP_04       540 RSLGYSDLFCLSKDDLMEAVTEYPDAKKVLEERGREILMKEGLLDENEVAASMEVDVQEK
                540 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       540 RSLGYSDLFCLSKDDLMEAVTEYPDAKKVLEERGREILMKEGLLDENEVAASMEVDVQEK

ref|XP_04       600 LEQLETNMETLYTRFGRLLAEYTGAQQKLKQRITVLETKMKQNNEDDYLSDGMNSPEPAA
                600 |.|||||||||||||||||||||||||||||||||||.|||||.||||||||||||||||
Query_744       600 LKQLETNMETLYTRFGRLLAEYTGAQQKLKQRITVLEVKMKQNTEDDYLSDGMNSPEPAA

ref|XP_04       660 ADEP 664
                660 |..| 664
Query_744       660 AEQP 664
""",
        )
        hit = record[8]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "ref|XP_035942617.1|")
        self.assertEqual(hit.target.name, "XP_035942617")
        self.assertEqual(
            hit.target.description,
            "cyclic nucleotide-gated olfactory channel [Halichoerus grypus]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=664)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 3221.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 1245.9)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.0)
        self.assertEqual(hsp.annotations["identity"], 638)
        self.assertEqual(hsp.annotations["positive"], 651)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0, 664],
                          [  0, 664]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 664))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MTEKSNGVKSSPANNHNNHVPATIKANGKDESRTRSRPQSAADDDTSSELQRLA...EQP')",
        )
        self.assertEqual(hsp.query.id, "Query_74414")
        self.assertEqual(hsp.query.description, "unnamed protein product")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MTEKSNGVKSSPANNHNHHAPPVIKANGKDDHRTSSRPQSAADDDTSSELQRLA...DEP')",
        )
        self.assertEqual(hsp.target.id, "ref|XP_035942617.1|")
        self.assertEqual(hsp.target.name, "XP_035942617")
        self.assertEqual(
            hsp.target.description,
            "cyclic nucleotide-gated olfactory channel [Halichoerus grypus]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MTEKSNGVKSSPANNHN+H P  IKANGKD+ RT SRPQSAADDDTSSELQRLAEMD PQQ RGGFRRIVRLVG+IR+WAN+NFREEE RPDSFLERFRGPELQTVTTQQGDGKGDKDG+GKGTKKKFELFVLDPAGDWYYRWLFVIAM VLYNWCLLVARACFSDLQ+GY+LVWLVLDYFSDVVYI DLFIRLRTGFLEQGLLVKDPKKLRDNYIHTLQFKLDVASIIPTDLIYFAVGIH+PELRFNRLLHFARMFEFFDRTETRTSYPNIFRISNLVLYILVIIHWNACIYYAISKSIGFGVDTWVYPNITDPEYGYLAREYIYCLYWSTLTLTTIGETPPPVKDEEYLFVIFDFLIGVLIFATIVGNVGSMISNMNATRAEFQAKIDAVKHYMQFRKVSKEMEAKVIKWFDYLWTNKK+VDEREVLKNLPAKLRAEIAINVHLSTLKKVRIFQDCEAGLLVELVLKLRPQVFSPGDYICRKGDIGKEMYIIKEGKLAVVADDGVTQYALLSAGSCFGEISILNIKGSKMGNRRTANIRSLGYSDLFCLSKDDLMEAVTEYPDAKKVLEERGREILMKEGLLDENEVAASMEVDVQEKL+QLETNMETLYTRFGRLLAEYTGAQQKLKQRITVLE KMKQN  DDYLSDGMNSPEPAAA++P",
        )
        self.assertEqual(
            str(hsp),
            """\
ref|XP_03         0 MTEKSNGVKSSPANNHNHHAPPVIKANGKDDHRTSSRPQSAADDDTSSELQRLAEMDVPQ
                  0 |||||||||||||||||.|.|..|||||||..||.||||||||||||||||||||||.||
Query_744         0 MTEKSNGVKSSPANNHNNHVPATIKANGKDESRTRSRPQSAADDDTSSELQRLAEMDAPQ

ref|XP_03        60 QGRGGFRRIVRLVGIIREWANKNFREEELRPDSFLERFRGPELQTVTTQQGDGKGDKDGE
                 60 |.||||||||||||.||.|||.||||||.||||||||||||||||||||||||||||||.
Query_744        60 QRRGGFRRIVRLVGVIRQWANRNFREEEARPDSFLERFRGPELQTVTTQQGDGKGDKDGD

ref|XP_03       120 GKGTKKKFELFVLDPAGDWYYRWLFVIAMLVLYNWCLLVARACFSDLQKGYYLVWLVLDY
                120 |||||||||||||||||||||||||||||.||||||||||||||||||.||.||||||||
Query_744       120 GKGTKKKFELFVLDPAGDWYYRWLFVIAMPVLYNWCLLVARACFSDLQRGYFLVWLVLDY

ref|XP_03       180 FSDVVYITDLFIRLRTGFLEQGLLVKDPKKLRDNYIHTLQFKLDVASIIPTDLIYFAVGI
                180 |||||||.||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       180 FSDVVYIADLFIRLRTGFLEQGLLVKDPKKLRDNYIHTLQFKLDVASIIPTDLIYFAVGI

ref|XP_03       240 HSPELRFNRLLHFARMFEFFDRTETRTSYPNIFRISNLVLYILVIIHWNACIYYAISKSI
                240 |.||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       240 HNPELRFNRLLHFARMFEFFDRTETRTSYPNIFRISNLVLYILVIIHWNACIYYAISKSI

ref|XP_03       300 GFGVDTWVYPNITDPEYGYLAREYIYCLYWSTLTLTTIGETPPPVKDEEYLFVIFDFLIG
                300 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       300 GFGVDTWVYPNITDPEYGYLAREYIYCLYWSTLTLTTIGETPPPVKDEEYLFVIFDFLIG

ref|XP_03       360 VLIFATIVGNVGSMISNMNATRAEFQAKIDAVKHYMQFRKVSKEMEAKVIKWFDYLWTNK
                360 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       360 VLIFATIVGNVGSMISNMNATRAEFQAKIDAVKHYMQFRKVSKEMEAKVIKWFDYLWTNK

ref|XP_03       420 KSVDEREVLKNLPAKLRAEIAINVHLSTLKKVRIFQDCEAGLLVELVLKLRPQVFSPGDY
                420 |.||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       420 KTVDEREVLKNLPAKLRAEIAINVHLSTLKKVRIFQDCEAGLLVELVLKLRPQVFSPGDY

ref|XP_03       480 ICRKGDIGKEMYIIKEGKLAVVADDGVTQYALLSAGSCFGEISILNIKGSKMGNRRTANI
                480 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       480 ICRKGDIGKEMYIIKEGKLAVVADDGVTQYALLSAGSCFGEISILNIKGSKMGNRRTANI

ref|XP_03       540 RSLGYSDLFCLSKDDLMEAVTEYPDAKKVLEERGREILMKEGLLDENEVAASMEVDVQEK
                540 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       540 RSLGYSDLFCLSKDDLMEAVTEYPDAKKVLEERGREILMKEGLLDENEVAASMEVDVQEK

ref|XP_03       600 LEQLETNMETLYTRFGRLLAEYTGAQQKLKQRITVLETKMKQNNMDDYLSDGMNSPEPAA
                600 |.|||||||||||||||||||||||||||||||||||.|||||..|||||||||||||||
Query_744       600 LKQLETNMETLYTRFGRLLAEYTGAQQKLKQRITVLEVKMKQNTEDDYLSDGMNSPEPAA

ref|XP_03       660 ADEP 664
                660 |..| 664
Query_744       660 AEQP 664
""",
        )
        hit = record[9]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "ref|XP_049729369.1|")
        self.assertEqual(hit.target.name, "XP_049729369")
        self.assertEqual(
            hit.target.description,
            "cyclic nucleotide-gated olfactory channel [Elephas maximus indicus]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=664)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 3219.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 1245.12)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.0)
        self.assertEqual(hsp.annotations["identity"], 635)
        self.assertEqual(hsp.annotations["positive"], 654)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0, 664],
                          [  0, 664]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 664))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MTEKSNGVKSSPANNHNNHVPATIKANGKDESRTRSRPQSAADDDTSSELQRLA...EQP')",
        )
        self.assertEqual(hsp.query.id, "Query_74414")
        self.assertEqual(hsp.query.description, "unnamed protein product")
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MTEKSNGVKSSPANNHNHHVPSTIKANGKDDRRTSSRPQSAADDDTSSELQRLA...EKP')",
        )
        self.assertEqual(hsp.target.id, "ref|XP_049729369.1|")
        self.assertEqual(hsp.target.name, "XP_049729369")
        self.assertEqual(
            hsp.target.description,
            "cyclic nucleotide-gated olfactory channel [Elephas maximus indicus]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MTEKSNGVKSSPANNHN+HVP+TIKANGKD+ RT SRPQSAADDDTSSELQRLAEMDAPQQ RGGFRRI+RLVGVIR+WAN+NFREE+ RPDSFLERFRGPELQTVTTQQGDGK DKDG+GKGTKKKFELFVLDPAGDWYYRWLF IA+PVLYNWCLLVARACFSDLQ+GY+LVWLVLDYFSD+VYIADLFIRLRTGFLEQGLLVKDPKKLRDNYIHT+QFKLDVASIIPTDLIYFAVGIH+PELRFNRLLHFARMFEFFDRTETRTSYPNIFRISNLVLYILVIIHWNACIYYAISKSIGFGVDTWVYPNITDP YGYLAREYIYCLYWSTLTLTTIGETPPPVKDEEYLFVIFDFLIGVLIFATIVGNVGSMISNMNATRAEFQAKIDAVKHYMQFRKVSKEMEAKVIKWFDYLWTNKKTVDEREVLKNLPAKLRAEIAINVHLSTLKKVRIFQDCEAGLLVELVLKLRPQVFSPGDYICRKGDIGKEMYIIKEGKLAVVADDGVTQYALLSAGSCFGEISILNIKGSKMGNRRTANIRSLGYSDLFCLSKDDLMEAVTEYPDAKKVLEERGREILMKEGLLDENEVAA+MEVDVQEKL+QLETNMETLYTRFGRLLAEYTGAQQKLKQRITVLE KMKQN E+DYLSDG+NSPEPAA E+P",
        )
        self.assertEqual(
            str(hsp),
            """\
ref|XP_04         0 MTEKSNGVKSSPANNHNHHVPSTIKANGKDDRRTSSRPQSAADDDTSSELQRLAEMDAPQ
                  0 |||||||||||||||||.|||.||||||||..||.|||||||||||||||||||||||||
Query_744         0 MTEKSNGVKSSPANNHNNHVPATIKANGKDESRTRSRPQSAADDDTSSELQRLAEMDAPQ

ref|XP_04        60 QWRGGFRRIIRLVGVIREWANKNFREEDPRPDSFLERFRGPELQTVTTQQGDGKSDKDGE
                 60 |.|||||||.|||||||.|||.|||||..|||||||||||||||||||||||||.||||.
Query_744        60 QRRGGFRRIVRLVGVIRQWANRNFREEEARPDSFLERFRGPELQTVTTQQGDGKGDKDGD

ref|XP_04       120 GKGTKKKFELFVLDPAGDWYYRWLFFIALPVLYNWCLLVARACFSDLQKGYYLVWLVLDY
                120 |||||||||||||||||||||||||.||.|||||||||||||||||||.||.||||||||
Query_744       120 GKGTKKKFELFVLDPAGDWYYRWLFVIAMPVLYNWCLLVARACFSDLQRGYFLVWLVLDY

ref|XP_04       180 FSDMVYIADLFIRLRTGFLEQGLLVKDPKKLRDNYIHTMQFKLDVASIIPTDLIYFAVGI
                180 |||.||||||||||||||||||||||||||||||||||.|||||||||||||||||||||
Query_744       180 FSDVVYIADLFIRLRTGFLEQGLLVKDPKKLRDNYIHTLQFKLDVASIIPTDLIYFAVGI

ref|XP_04       240 HSPELRFNRLLHFARMFEFFDRTETRTSYPNIFRISNLVLYILVIIHWNACIYYAISKSI
                240 |.||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       240 HNPELRFNRLLHFARMFEFFDRTETRTSYPNIFRISNLVLYILVIIHWNACIYYAISKSI

ref|XP_04       300 GFGVDTWVYPNITDPAYGYLAREYIYCLYWSTLTLTTIGETPPPVKDEEYLFVIFDFLIG
                300 |||||||||||||||.||||||||||||||||||||||||||||||||||||||||||||
Query_744       300 GFGVDTWVYPNITDPEYGYLAREYIYCLYWSTLTLTTIGETPPPVKDEEYLFVIFDFLIG

ref|XP_04       360 VLIFATIVGNVGSMISNMNATRAEFQAKIDAVKHYMQFRKVSKEMEAKVIKWFDYLWTNK
                360 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       360 VLIFATIVGNVGSMISNMNATRAEFQAKIDAVKHYMQFRKVSKEMEAKVIKWFDYLWTNK

ref|XP_04       420 KTVDEREVLKNLPAKLRAEIAINVHLSTLKKVRIFQDCEAGLLVELVLKLRPQVFSPGDY
                420 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       420 KTVDEREVLKNLPAKLRAEIAINVHLSTLKKVRIFQDCEAGLLVELVLKLRPQVFSPGDY

ref|XP_04       480 ICRKGDIGKEMYIIKEGKLAVVADDGVTQYALLSAGSCFGEISILNIKGSKMGNRRTANI
                480 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Query_744       480 ICRKGDIGKEMYIIKEGKLAVVADDGVTQYALLSAGSCFGEISILNIKGSKMGNRRTANI

ref|XP_04       540 RSLGYSDLFCLSKDDLMEAVTEYPDAKKVLEERGREILMKEGLLDENEVAATMEVDVQEK
                540 |||||||||||||||||||||||||||||||||||||||||||||||||||.||||||||
Query_744       540 RSLGYSDLFCLSKDDLMEAVTEYPDAKKVLEERGREILMKEGLLDENEVAASMEVDVQEK

ref|XP_04       600 LEQLETNMETLYTRFGRLLAEYTGAQQKLKQRITVLETKMKQNNEEDYLSDGINSPEPAA
                600 |.|||||||||||||||||||||||||||||||||||.|||||.|.||||||.|||||||
Query_744       600 LKQLETNMETLYTRFGRLLAEYTGAQQKLKQRITVLEVKMKQNTEDDYLSDGMNSPEPAA

ref|XP_04       660 VEKP 664
                660 .|.| 664
Query_744       660 AEQP 664
""",
        )


class TestBlastn(unittest.TestCase):
    """Test the Blast XML parser for blastn output."""

    def test_xml_2900_blastn_001(self):
        """Parsing BLASTN 2.9.0+ (xml_2900_blastn_001.xml)."""
        filename = "xml_2900_blastn_001.xml"
        datafile = os.path.join("Blast", filename)
        with open(datafile, "rb") as handle:
            records = Blast.parse(handle)
            self.check_xml_2900_blastn_001_records(records)
        with Blast.parse(datafile) as records:
            self.check_xml_2900_blastn_001_records(records)
        with open(datafile, "rb") as handle:
            record = Blast.read(handle)
        self.check_xml_2900_blastn_001_record(record)
        record = Blast.read(datafile)
        self.check_xml_2900_blastn_001_record(record)

    def check_xml_2900_blastn_001_records(self, records):
        self.assertEqual(records.program, "blastn")
        self.assertEqual(records.version, "BLASTN 2.9.0+")
        self.assertEqual(
            records.reference,
            'Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.',
        )
        self.assertEqual(
            records.db, "GPIPE/10090/current/all_top_level GPIPE/10090/current/rna"
        )
        self.assertIsInstance(records.query, SeqRecord)
        self.assertEqual(records.query.id, "G26684.1")
        self.assertEqual(
            records.query.description, "human STS STS_D11570, sequence tagged site"
        )
        self.assertEqual(repr(records.query.seq), "Seq(None, length=285)")
        self.assertEqual(len(records.param), 6)
        self.assertAlmostEqual(records.param["expect"], 10.0)
        self.assertEqual(records.param["sc-match"], 2)
        self.assertEqual(records.param["sc-mismatch"], -3)
        self.assertEqual(records.param["gap-open"], 5)
        self.assertEqual(records.param["gap-extend"], 2)
        self.assertEqual(records.param["filter"], "R -d repeatmasker/repeat_9989;m;F;")
        record = next(records)
        self.assertRaises(StopIteration, next, records)
        self.check_xml_2900_blastn_001_record(record)

    def check_xml_2900_blastn_001_record(self, record):
        self.assertEqual(record.num, 1)

        self.assertIsInstance(record.query, SeqRecord)
        self.assertEqual(record.query.id, "G26684.1")
        self.assertEqual(
            record.query.description, "human STS STS_D11570, sequence tagged site"
        )
        self.assertEqual(repr(record.query.seq), "Seq(None, length=285)")

        self.assertEqual(len(record.stat), 7)
        self.assertEqual(record.stat["db-num"], 107382)
        self.assertEqual(record.stat["db-len"], 3164670549)
        self.assertEqual(record.stat["hsp-len"], 0)
        self.assertAlmostEqual(record.stat["eff-space"], 0.0)
        self.assertAlmostEqual(record.stat["kappa"], 0.41)
        self.assertAlmostEqual(record.stat["lambda"], 0.625)
        self.assertAlmostEqual(record.stat["entropy"], 0.78)
        self.assertEqual(len(record), 10)

        hit = record[0]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|372099107|ref|NC_000069.6|")
        self.assertEqual(hit.target.name, "NC_000069")
        self.assertEqual(
            hit.target.description,
            "Mus musculus strain C57BL/6J chromosome 3, GRCm38.p4 C57BL/6J",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=160039680)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 44.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 40.9604)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.375311)
        self.assertEqual(hsp.annotations["identity"], 30)
        self.assertEqual(hsp.annotations["positive"], 30)
        self.assertEqual(hsp.annotations["gaps"], 1)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[101449177, 101449150, 101449149, 101449143],
                          [      133,       160,       160,       166]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 34))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({133: 'GAATCCTAGAGGCTTGATTGGCCCAGGCTGCTG'}, length=285)",
        )
        self.assertEqual(hsp.query.id, "G26684.1")
        self.assertEqual(
            hsp.query.description, "human STS STS_D11570, sequence tagged site"
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({101449143: 'CAGCAGGCCAGGGCCAGTCCAGCCTCTAGGATTC'}, length=160039680)",
        )
        self.assertEqual(hsp.target.id, "gi|372099107|ref|NC_000069.6|")
        self.assertEqual(hsp.target.name, "NC_000069")
        self.assertEqual(
            hsp.target.description,
            "Mus musculus strain C57BL/6J chromosome 3, GRCm38.p4 C57BL/6J",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"], "|||||||||||||| || |||||| || ||||||"
        )
        self.assertEqual(
            str(hsp),
            """\
gi|372099 101449177 GAATCCTAGAGGCTGGACTGGCCCTGGCCTGCTG 101449143
                  0 ||||||||||||||.||.||||||.||-||||||        34
G26684.1        133 GAATCCTAGAGGCTTGATTGGCCCAGG-CTGCTG       166
""",
        )
        hit = record[1]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|372099103|ref|NC_000073.6|")
        self.assertEqual(hit.target.name, "NC_000073")
        self.assertEqual(
            hit.target.description,
            "Mus musculus strain C57BL/6J chromosome 7, GRCm38.p4 C57BL/6J",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=145441459)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 44.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 40.9604)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.375311)
        self.assertEqual(hsp.annotations["identity"], 26)
        self.assertEqual(hsp.annotations["positive"], 26)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[131772185, 131772156],
                          [      204,       233]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 29))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({204: 'GAAAGGAAATNAAAATGGAAAGTTCTTGT'}, length=285)",
        )
        self.assertEqual(hsp.query.id, "G26684.1")
        self.assertEqual(
            hsp.query.description, "human STS STS_D11570, sequence tagged site"
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({131772156: 'ACCAGAACTTTCCATTTTTTTTTCCTTTC'}, length=145441459)",
        )
        self.assertEqual(hsp.target.id, "gi|372099103|ref|NC_000073.6|")
        self.assertEqual(hsp.target.name, "NC_000073")
        self.assertEqual(
            hsp.target.description,
            "Mus musculus strain C57BL/6J chromosome 7, GRCm38.p4 C57BL/6J",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(hsp.annotations["midline"], "|||||||||  ||||||||||||||| ||")
        self.assertEqual(
            str(hsp),
            """\
gi|372099 131772185 GAAAGGAAAAAAAAATGGAAAGTTCTGGT 131772156
                  0 |||||||||..|||||||||||||||.||        29
G26684.1        204 GAAAGGAAATNAAAATGGAAAGTTCTTGT       233
""",
        )
        hit = record[2]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|372099106|ref|NC_000070.6|")
        self.assertEqual(hit.target.name, "NC_000070")
        self.assertEqual(
            hit.target.description,
            "Mus musculus strain C57BL/6J chromosome 4, GRCm38.p4 C57BL/6J",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=156508116)")
        self.assertEqual(len(hit), 2)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 43.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 40.0587)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.30996)
        self.assertEqual(hsp.annotations["identity"], 23)
        self.assertEqual(hsp.annotations["positive"], 23)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[9607562, 9607538],
                          [     61,      85]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 24))
        self.assertEqual(
            repr(hsp.query.seq), "Seq({61: 'CCAACACAGGCCAGCGACTTCTGG'}, length=285)"
        )
        self.assertEqual(hsp.query.id, "G26684.1")
        self.assertEqual(
            hsp.query.description, "human STS STS_D11570, sequence tagged site"
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({9607538: 'CCAGAAGCCGCTGGCCTGTGTTGG'}, length=156508116)",
        )
        self.assertEqual(hsp.target.id, "gi|372099106|ref|NC_000070.6|")
        self.assertEqual(hsp.target.name, "NC_000070")
        self.assertEqual(
            hsp.target.description,
            "Mus musculus strain C57BL/6J chromosome 4, GRCm38.p4 C57BL/6J",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(hsp.annotations["midline"], "|||||||||||||||| |||||||")
        self.assertEqual(
            str(hsp),
            """\
gi|372099   9607562 CCAACACAGGCCAGCGGCTTCTGG 9607538
                  0 ||||||||||||||||.|||||||      24
G26684.1         61 CCAACACAGGCCAGCGACTTCTGG      85
""",
        )
        hsp = hit[1]
        self.assertAlmostEqual(hsp.score, 40.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 37.3537)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.57222)
        self.assertEqual(hsp.annotations["identity"], 28)
        self.assertEqual(hsp.annotations["positive"], 28)
        self.assertEqual(hsp.annotations["gaps"], 1)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[142902531, 142902542, 142902543, 142902563],
                          [      241,       252,       252,       272]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 32))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({241: 'GCCTGACATGGGTAGCTGCTCAATAAATGCT'}, length=285)",
        )
        self.assertEqual(hsp.query.id, "G26684.1")
        self.assertEqual(
            hsp.query.description, "human STS STS_D11570, sequence tagged site"
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({142902531: 'GCCTGGCATGAAGTAACTGCTCAATAAATGCT'}, length=156508116)",
        )
        self.assertEqual(hsp.target.id, "gi|372099106|ref|NC_000070.6|")
        self.assertEqual(hsp.target.name, "NC_000070")
        self.assertEqual(
            hsp.target.description,
            "Mus musculus strain C57BL/6J chromosome 4, GRCm38.p4 C57BL/6J",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(hsp.annotations["midline"], "||||| ||||  ||| ||||||||||||||||")
        self.assertEqual(
            str(hsp),
            """\
gi|372099 142902531 GCCTGGCATGAAGTAACTGCTCAATAAATGCT 142902563
                  0 |||||.||||.-|||.||||||||||||||||        32
G26684.1        241 GCCTGACATGG-GTAGCTGCTCAATAAATGCT       272
""",
        )
        hit = record[3]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|372099108|ref|NC_000068.7|")
        self.assertEqual(hit.target.name, "NC_000068")
        self.assertEqual(
            hit.target.description,
            "Mus musculus strain C57BL/6J chromosome 2, GRCm38.p4 C57BL/6J",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=182113224)")
        self.assertEqual(len(hit), 2)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 42.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 39.157)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.30996)
        self.assertEqual(hsp.annotations["identity"], 27)
        self.assertEqual(hsp.annotations["positive"], 27)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[3799646, 3799677],
                          [    238,     269]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 31))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({238: 'AAGGCCTGACATGGGTAGCTGCTCAATAAAT'}, length=285)",
        )
        self.assertEqual(hsp.query.id, "G26684.1")
        self.assertEqual(
            hsp.query.description, "human STS STS_D11570, sequence tagged site"
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({3799646: 'AAGTCCTGGCATGAGTAGTTGCTCAATAAAT'}, length=182113224)",
        )
        self.assertEqual(hsp.target.id, "gi|372099108|ref|NC_000068.7|")
        self.assertEqual(hsp.target.name, "NC_000068")
        self.assertEqual(
            hsp.target.description,
            "Mus musculus strain C57BL/6J chromosome 2, GRCm38.p4 C57BL/6J",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(hsp.annotations["midline"], "||| |||| |||| |||| ||||||||||||")
        self.assertEqual(
            str(hsp),
            """\
gi|372099   3799646 AAGTCCTGGCATGAGTAGTTGCTCAATAAAT 3799677
                  0 |||.||||.||||.||||.||||||||||||      31
G26684.1        238 AAGGCCTGACATGGGTAGCTGCTCAATAAAT     269
""",
        )
        hsp = hit[1]
        self.assertAlmostEqual(hsp.score, 41.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 38.2554)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.57222)
        self.assertEqual(hsp.annotations["identity"], 23)
        self.assertEqual(hsp.annotations["positive"], 23)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[70278959, 70278984],
                          [     210,      235]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 25))
        self.assertEqual(
            repr(hsp.query.seq), "Seq({210: 'AAATNAAAATGGAAAGTTCTTGTAG'}, length=285)"
        )
        self.assertEqual(hsp.query.id, "G26684.1")
        self.assertEqual(
            hsp.query.description, "human STS STS_D11570, sequence tagged site"
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({70278959: 'AAATGAAAATGGAAAGTTCTTATAG'}, length=182113224)",
        )
        self.assertEqual(hsp.target.id, "gi|372099108|ref|NC_000068.7|")
        self.assertEqual(hsp.target.name, "NC_000068")
        self.assertEqual(
            hsp.target.description,
            "Mus musculus strain C57BL/6J chromosome 2, GRCm38.p4 C57BL/6J",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(hsp.annotations["midline"], "|||| |||||||||||||||| |||")
        self.assertEqual(
            str(hsp),
            """\
gi|372099  70278959 AAATGAAAATGGAAAGTTCTTATAG 70278984
                  0 ||||.||||||||||||||||.|||       25
G26684.1        210 AAATNAAAATGGAAAGTTCTTGTAG      235
""",
        )
        hit = record[4]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|372099097|ref|NC_000079.6|")
        self.assertEqual(hit.target.name, "NC_000079")
        self.assertEqual(
            hit.target.description,
            "Mus musculus strain C57BL/6J chromosome 13, GRCm38.p4 C57BL/6J",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=120421639)")
        self.assertEqual(len(hit), 2)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 42.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 39.157)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.30996)
        self.assertEqual(hsp.annotations["identity"], 25)
        self.assertEqual(hsp.annotations["positive"], 25)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[26806584, 26806556],
                          [     206,      234]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 28))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({206: 'AAGGAAATNAAAATGGAAAGTTCTTGTA'}, length=285)",
        )
        self.assertEqual(hsp.query.id, "G26684.1")
        self.assertEqual(
            hsp.query.description, "human STS STS_D11570, sequence tagged site"
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({26806556: 'TAGAAGAACTTTCCATTTTGATGTCCTT'}, length=120421639)",
        )
        self.assertEqual(hsp.target.id, "gi|372099097|ref|NC_000079.6|")
        self.assertEqual(hsp.target.name, "NC_000079")
        self.assertEqual(
            hsp.target.description,
            "Mus musculus strain C57BL/6J chromosome 13, GRCm38.p4 C57BL/6J",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(hsp.annotations["midline"], "||||| || |||||||||||||||| ||")
        self.assertEqual(
            str(hsp),
            """\
gi|372099  26806584 AAGGACATCAAAATGGAAAGTTCTTCTA 26806556
                  0 |||||.||.||||||||||||||||.||       28
G26684.1        206 AAGGAAATNAAAATGGAAAGTTCTTGTA      234
""",
        )
        hsp = hit[1]
        self.assertAlmostEqual(hsp.score, 40.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 37.3537)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.57222)
        self.assertEqual(hsp.annotations["identity"], 32)
        self.assertEqual(hsp.annotations["positive"], 32)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[56840340, 56840300],
                          [     233,      273]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 40))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({233: 'AGCGCAAGGCCTGACATGGGTAGCTGCTCAATAAATGCTA'}, length=285)",
        )
        self.assertEqual(hsp.query.id, "G26684.1")
        self.assertEqual(
            hsp.query.description, "human STS STS_D11570, sequence tagged site"
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({56840300: 'TAGTATTCACTGAACATTTTCCTATGTCAGGCCTTGCGCT'}, length=120421639)",
        )
        self.assertEqual(hsp.target.id, "gi|372099097|ref|NC_000079.6|")
        self.assertEqual(hsp.target.name, "NC_000079")
        self.assertEqual(
            hsp.target.description,
            "Mus musculus strain C57BL/6J chromosome 13, GRCm38.p4 C57BL/6J",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"], "||||||||||||||||| || |  || ||| | ||| |||"
        )
        self.assertEqual(
            str(hsp),
            """\
gi|372099  56840340 AGCGCAAGGCCTGACATAGGAAAATGTTCAGTGAATACTA 56840300
                  0 |||||||||||||||||.||.|..||.|||.|.|||.|||       40
G26684.1        233 AGCGCAAGGCCTGACATGGGTAGCTGCTCAATAAATGCTA      273
""",
        )
        hit = record[5]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|372099098|ref|NC_000078.6|")
        self.assertEqual(hit.target.name, "NC_000078")
        self.assertEqual(
            hit.target.description,
            "Mus musculus strain C57BL/6J chromosome 12, GRCm38.p4 C57BL/6J",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=120129022)")
        self.assertEqual(len(hit), 2)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 41.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 38.2554)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.57222)
        self.assertEqual(hsp.annotations["identity"], 22)
        self.assertEqual(hsp.annotations["positive"], 22)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[113030662, 113030685],
                          [       48,        71]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 23))
        self.assertEqual(
            repr(hsp.query.seq), "Seq({48: 'CATCCATTCACACCCAACACAGG'}, length=285)"
        )
        self.assertEqual(hsp.query.id, "G26684.1")
        self.assertEqual(
            hsp.query.description, "human STS STS_D11570, sequence tagged site"
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({113030662: 'CATCCATTCACACCCAGCACAGG'}, length=120129022)",
        )
        self.assertEqual(hsp.target.id, "gi|372099098|ref|NC_000078.6|")
        self.assertEqual(hsp.target.name, "NC_000078")
        self.assertEqual(
            hsp.target.description,
            "Mus musculus strain C57BL/6J chromosome 12, GRCm38.p4 C57BL/6J",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(hsp.annotations["midline"], "|||||||||||||||| ||||||")
        self.assertEqual(
            str(hsp),
            """\
gi|372099 113030662 CATCCATTCACACCCAGCACAGG 113030685
                  0 ||||||||||||||||.||||||        23
G26684.1         48 CATCCATTCACACCCAACACAGG        71
""",
        )
        hsp = hit[1]
        self.assertAlmostEqual(hsp.score, 40.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 37.3537)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.57222)
        self.assertEqual(hsp.annotations["identity"], 28)
        self.assertEqual(hsp.annotations["positive"], 28)
        self.assertEqual(hsp.annotations["gaps"], 1)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[108990272, 108990248, 108990248, 108990241],
                          [      230,       254,       255,       262]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 32))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({230: 'TGTAGCGCAAGGCCTGACATGGGTAGCTGCTC'}, length=285)",
        )
        self.assertEqual(hsp.query.id, "G26684.1")
        self.assertEqual(
            hsp.query.description, "human STS STS_D11570, sequence tagged site"
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({108990241: 'GACCAGCACCCATGTCAGGCCTAGAGCTACA'}, length=120129022)",
        )
        self.assertEqual(hsp.target.id, "gi|372099098|ref|NC_000078.6|")
        self.assertEqual(hsp.target.name, "NC_000078")
        self.assertEqual(
            hsp.target.description,
            "Mus musculus strain C57BL/6J chromosome 12, GRCm38.p4 C57BL/6J",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(hsp.annotations["midline"], "|||||| | ||||||||||||||| |||| ||")
        self.assertEqual(
            str(hsp),
            """\
gi|372099 108990272 TGTAGCTCTAGGCCTGACATGGGT-GCTGGTC 108990241
                  0 ||||||.|.|||||||||||||||-||||.||        32
G26684.1        230 TGTAGCGCAAGGCCTGACATGGGTAGCTGCTC       262
""",
        )
        hit = record[6]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|372099109|ref|NC_000067.6|")
        self.assertEqual(hit.target.name, "NC_000067")
        self.assertEqual(
            hit.target.description,
            "Mus musculus strain C57BL/6J chromosome 1, GRCm38.p4 C57BL/6J",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=195471971)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 40.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 37.3537)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.57222)
        self.assertEqual(hsp.annotations["identity"], 35)
        self.assertEqual(hsp.annotations["positive"], 35)
        self.assertEqual(hsp.annotations["gaps"], 2)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[65190107, 65190128, 65190128, 65190144, 65190144, 65190148],
                          [      86,      107,      108,      124,      125,      129]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 43))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({86: 'GCTCAGCCACAGACATGGTTTGTNACTNTTGAGCTTCTGTTCC'}, length=285)",
        )
        self.assertEqual(hsp.query.id, "G26684.1")
        self.assertEqual(
            hsp.query.description, "human STS STS_D11570, sequence tagged site"
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({65190107: 'GCTCAGCCACATACATGGTTTTAAGTGTTGAGGCTCTTTCC'}, length=195471971)",
        )
        self.assertEqual(hsp.target.id, "gi|372099109|ref|NC_000067.6|")
        self.assertEqual(hsp.target.name, "NC_000067")
        self.assertEqual(
            hsp.target.description,
            "Mus musculus strain C57BL/6J chromosome 1, GRCm38.p4 C57BL/6J",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"], "||||||||||| ||||||||| | | | |||||  ||| ||||"
        )
        self.assertEqual(
            str(hsp),
            """\
gi|372099  65190107 GCTCAGCCACATACATGGTTT-TAAGTGTTGAGGCTCT-TTCC 65190148
                  0 |||||||||||.|||||||||-|.|.|.|||||..|||-||||       43
G26684.1         86 GCTCAGCCACAGACATGGTTTGTNACTNTTGAGCTTCTGTTCC      129
""",
        )
        hit = record[7]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|372099101|ref|NC_000075.6|")
        self.assertEqual(hit.target.name, "NC_000075")
        self.assertEqual(
            hit.target.description,
            "Mus musculus strain C57BL/6J chromosome 9, GRCm38.p4 C57BL/6J",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=124595110)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 40.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 37.3537)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.57222)
        self.assertEqual(hsp.annotations["identity"], 36)
        self.assertEqual(hsp.annotations["positive"], 36)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[58227241, 58227194],
                          [     237,      284]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 47))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({237: 'CAAGGCCTGACATGGGTAGCTGCTCAATAAATGCTAGTNTGTTATTT'}, length=285)",
        )
        self.assertEqual(hsp.query.id, "G26684.1")
        self.assertEqual(
            hsp.query.description, "human STS STS_D11570, sequence tagged site"
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({58227194: 'AAAAAAAAAAATAGTATTTATTGAGCAGTCATACCTGTCAGGCTTTG'}, length=124595110)",
        )
        self.assertEqual(hsp.target.id, "gi|372099101|ref|NC_000075.6|")
        self.assertEqual(hsp.target.name, "NC_000075")
        self.assertEqual(
            hsp.target.description,
            "Mus musculus strain C57BL/6J chromosome 9, GRCm38.p4 C57BL/6J",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "||| |||||||| |  |  ||||||||||||| ||| | | || |||",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|372099  58227241 CAAAGCCTGACAGGTATGACTGCTCAATAAATACTATTTTTTTTTTT 58227194
                  0 |||.||||||||.|..|..|||||||||||||.|||.|.|.||.|||       47
G26684.1        237 CAAGGCCTGACATGGGTAGCTGCTCAATAAATGCTAGTNTGTTATTT      284
""",
        )
        hit = record[8]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|372099100|ref|NC_000076.6|")
        self.assertEqual(hit.target.name, "NC_000076")
        self.assertEqual(
            hit.target.description,
            "Mus musculus strain C57BL/6J chromosome 10, GRCm38.p4 C57BL/6J",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=130694993)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 40.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 37.3537)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.57222)
        self.assertEqual(hsp.annotations["identity"], 20)
        self.assertEqual(hsp.annotations["positive"], 20)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[119337185, 119337205],
                          [      254,       274]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 20))
        self.assertEqual(
            repr(hsp.query.seq), "Seq({254: 'AGCTGCTCAATAAATGCTAG'}, length=285)"
        )
        self.assertEqual(hsp.query.id, "G26684.1")
        self.assertEqual(
            hsp.query.description, "human STS STS_D11570, sequence tagged site"
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({119337185: 'AGCTGCTCAATAAATGCTAG'}, length=130694993)",
        )
        self.assertEqual(hsp.target.id, "gi|372099100|ref|NC_000076.6|")
        self.assertEqual(hsp.target.name, "NC_000076")
        self.assertEqual(
            hsp.target.description,
            "Mus musculus strain C57BL/6J chromosome 10, GRCm38.p4 C57BL/6J",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(hsp.annotations["midline"], "||||||||||||||||||||")
        self.assertEqual(
            str(hsp),
            """\
gi|372099 119337185 AGCTGCTCAATAAATGCTAG 119337205
                  0 ||||||||||||||||||||        20
G26684.1        254 AGCTGCTCAATAAATGCTAG       274
""",
        )
        hit = record[9]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|372099094|ref|NC_000082.6|")
        self.assertEqual(hit.target.name, "NC_000082")
        self.assertEqual(
            hit.target.description,
            "Mus musculus strain C57BL/6J chromosome 16, GRCm38.p4 C57BL/6J",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=98207768)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 40.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 37.3537)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.57222)
        self.assertEqual(hsp.annotations["identity"], 43)
        self.assertEqual(hsp.annotations["positive"], 43)
        self.assertEqual(hsp.annotations["gaps"], 2)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[18854779, 18854803, 18854804, 18854812, 18854813, 18854835],
                          [     174,      198,      198,      206,      206,      228]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 56))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({174: 'GGAGGCAAAGAATCCCTACCTCCTAGGGGTGAAAGGAAATNAAAATGGAAAGTT'}, length=285)",
        )
        self.assertEqual(hsp.query.id, "G26684.1")
        self.assertEqual(
            hsp.query.description, "human STS STS_D11570, sequence tagged site"
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({18854779: 'GGAGGCAAAGAATCCCTACATTGTGACAGCTGATAAAGAAGGTAAAATGGAAAATT'}, length=98207768)",
        )
        self.assertEqual(hsp.target.id, "gi|372099094|ref|NC_000082.6|")
        self.assertEqual(hsp.target.name, "NC_000082")
        self.assertEqual(
            hsp.target.description,
            "Mus musculus strain C57BL/6J chromosome 16, GRCm38.p4 C57BL/6J",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "||||||||||||||||||| |  | |  | ||| || |||   |||||||||| ||",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|372099  18854779 GGAGGCAAAGAATCCCTACATTGTGACAGCTGATAAAGAAGGTAAAATGGAAAATT
                  0 |||||||||||||||||||.|..|-|..|.|||-||.|||...||||||||||.||
G26684.1        174 GGAGGCAAAGAATCCCTACCTCCT-AGGGGTGA-AAGGAAATNAAAATGGAAAGTT

gi|372099  18854835
                 56
G26684.1        228
""",
        )

    def test_megablast_legacy(self):
        """Parsing megablast 2.2.26 [Sep-21-2011] (megablast_legacy.xml)."""
        filename = "megablast_legacy.xml"
        datafile = os.path.join("Blast", filename)
        with open(datafile, "rb") as handle:
            records = Blast.parse(handle)
            self.check_megablast_legacy_records(records)
        with Blast.parse(datafile) as records:
            self.check_megablast_legacy_records(records)
        with open(datafile, "rb") as handle:
            record = Blast.read(handle)
        self.check_megablast_legacy_record(record)
        record = Blast.read(datafile)
        self.check_megablast_legacy_record(record)

    def check_megablast_legacy_records(self, records):
        self.assertEqual(records.program, "megablast")
        self.assertEqual(records.version, "megablast 2.2.26 [Sep-21-2011]")
        self.assertEqual(
            records.reference,
            '~Reference: Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schaffer, ~Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), ~"Gapped BLAST and PSI-BLAST: a new generation of protein database search~programs",  Nucleic Acids Res. 25:3389-3402.',
        )
        self.assertEqual(records.db, "m_cold.fasta")
        self.assertIsInstance(records.query, SeqRecord)
        self.assertEqual(records.query.id, "lcl|1_")
        self.assertEqual(
            records.query.description,
            "gi|8332116|gb|BE037100.1|BE037100 MP14H09 MP Mesembryanthemum crystallinum cDNA 5' similar to cold acclimation protein, mRNA sequence",
        )
        self.assertEqual(repr(records.query.seq), "Seq(None, length=1111)")
        self.assertEqual(len(records.param), 6)
        self.assertAlmostEqual(records.param["expect"], 10.0)
        self.assertEqual(records.param["sc-match"], 1)
        self.assertEqual(records.param["sc-mismatch"], -3)
        self.assertEqual(records.param["gap-open"], 0)
        self.assertEqual(records.param["gap-extend"], 0)
        self.assertEqual(records.param["filter"], "D")
        record = next(records)
        self.assertRaises(StopIteration, next, records)
        self.check_megablast_legacy_record(record)
        self.assertEqual(len(records.mbstat), 7)
        self.assertEqual(records.mbstat["db-num"], 1)
        self.assertEqual(records.mbstat["db-len"], 1111)
        self.assertEqual(records.mbstat["hsp-len"], 10)
        self.assertAlmostEqual(records.mbstat["eff-space"], 1212200.0)
        self.assertAlmostEqual(records.mbstat["kappa"], 0.710603)
        self.assertAlmostEqual(records.mbstat["lambda"], 1.37406)
        self.assertAlmostEqual(records.mbstat["entropy"], 1.30725)

    def check_megablast_legacy_record(self, record):
        self.assertEqual(record.num, 0)
        self.assertIsInstance(record.query, SeqRecord)
        self.assertEqual(record.query.id, "lcl|1_")
        self.assertEqual(
            record.query.description,
            "gi|8332116|gb|BE037100.1|BE037100 MP14H09 MP Mesembryanthemum crystallinum cDNA 5' similar to cold acclimation protein, mRNA sequence",
        )
        self.assertEqual(repr(record.query.seq), "Seq(None, length=1111)")
        self.assertEqual(len(record), 1)
        hit = record[0]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gnl|BL_ORD_ID|0")
        self.assertEqual(hit.target.name, "0")
        self.assertEqual(
            hit.target.description,
            "gi|8332116|gb|BE037100.1|BE037100 MP14H09 MP Mesembryanthemum crystallinum cDNA 5' similar to cold acclimation protein, mRNA sequence",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=1111)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 788.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 1562.59)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.0)
        self.assertEqual(hsp.annotations["identity"], 797)
        self.assertEqual(hsp.annotations["positive"], 797)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0, 797],
                          [  0, 797]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 797))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({0: 'CACTAGTACTCGAGCGTNCTGCACCAATTCGGCACGAGCAAGTGACTACGTTNT...GTG'}, length=1111)",
        )
        self.assertEqual(hsp.query.id, "lcl|1_")
        self.assertEqual(
            hsp.query.description,
            "gi|8332116|gb|BE037100.1|BE037100 MP14H09 MP Mesembryanthemum crystallinum cDNA 5' similar to cold acclimation protein, mRNA sequence",
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'CACTAGTACTCGAGCGTNCTGCACCAATTCGGCACGAGCAAGTGACTACGTTNT...GTG'}, length=1111)",
        )
        self.assertEqual(hsp.target.id, "gnl|BL_ORD_ID|0")
        self.assertEqual(hsp.target.name, "0")
        self.assertEqual(
            hsp.target.description,
            "gi|8332116|gb|BE037100.1|BE037100 MP14H09 MP Mesembryanthemum crystallinum cDNA 5' similar to cold acclimation protein, mRNA sequence",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||",
        )
        self.assertEqual(
            str(hsp),
            """\
gnl|BL_OR         0 CACTAGTACTCGAGCGTNCTGCACCAATTCGGCACGAGCAAGTGACTACGTTNTGTGAAC
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
lcl|1_            0 CACTAGTACTCGAGCGTNCTGCACCAATTCGGCACGAGCAAGTGACTACGTTNTGTGAAC

gnl|BL_OR        60 AGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTA
                 60 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
lcl|1_           60 AGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTA

gnl|BL_OR       120 ATATGATCGATTCCGATATCAATGAGCTTAAAATGGCAACAATGAGGCTCATCAATGATG
                120 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
lcl|1_          120 ATATGATCGATTCCGATATCAATGAGCTTAAAATGGCAACAATGAGGCTCATCAATGATG

gnl|BL_OR       180 CTAGTATGCTCGGTCATTACGGGTTTGGCACTCATTTCCTCAAATGGCTCGCCTGCCTTG
                180 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
lcl|1_          180 CTAGTATGCTCGGTCATTACGGGTTTGGCACTCATTTCCTCAAATGGCTCGCCTGCCTTG

gnl|BL_OR       240 CGGCTATTTACTTGTTGATATTGGATCGAACAAACTGGAGAACCAACATGCTCACGTCAC
                240 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
lcl|1_          240 CGGCTATTTACTTGTTGATATTGGATCGAACAAACTGGAGAACCAACATGCTCACGTCAC

gnl|BL_OR       300 TTTTAGTCCCTTACATATTCCTCAGTCTTCCATCCGGGCCATTTCATCTGTTCAGAGGCG
                300 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
lcl|1_          300 TTTTAGTCCCTTACATATTCCTCAGTCTTCCATCCGGGCCATTTCATCTGTTCAGAGGCG

gnl|BL_OR       360 AGGTCGGGAAATGGATTGCCATCATTGCAGTCGTGTTAAGGCTGTTCTTCAACCGGCATT
                360 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
lcl|1_          360 AGGTCGGGAAATGGATTGCCATCATTGCAGTCGTGTTAAGGCTGTTCTTCAACCGGCATT

gnl|BL_OR       420 TCCCAGTTTGGCTGGAAATGCCTGGATCGTTGATACTCCTCCTGGTGGTGGCACCAGACT
                420 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
lcl|1_          420 TCCCAGTTTGGCTGGAAATGCCTGGATCGTTGATACTCCTCCTGGTGGTGGCACCAGACT

gnl|BL_OR       480 TCTTTACACACAAAGTGAAGGAGAGCTGGATCGGAATTGCAATTATGATAGCGATAGGGT
                480 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
lcl|1_          480 TCTTTACACACAAAGTGAAGGAGAGCTGGATCGGAATTGCAATTATGATAGCGATAGGGT

gnl|BL_OR       540 GTCACCTGATGCAAGAACATATCAGAGCCACTGGTGGCTTTTGGAATTCCTTCACACAGA
                540 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
lcl|1_          540 GTCACCTGATGCAAGAACATATCAGAGCCACTGGTGGCTTTTGGAATTCCTTCACACAGA

gnl|BL_OR       600 GCCACGGAACTTTTAACACAATTGGGCTTATCCTTCTACTGGCTTACCCTGTCTGTTTAT
                600 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
lcl|1_          600 GCCACGGAACTTTTAACACAATTGGGCTTATCCTTCTACTGGCTTACCCTGTCTGTTTAT

gnl|BL_OR       660 GGTCATCTTCATGATGTAGTAGCTTAGTCTTGATCCTAATCCTCAAATNTACTTTTCCAG
                660 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
lcl|1_          660 GGTCATCTTCATGATGTAGTAGCTTAGTCTTGATCCTAATCCTCAAATNTACTTTTCCAG

gnl|BL_OR       720 CTCTTTCGACGCTCTTGCTAAAGCCCATTCAATTCGCCCCATATTTCGCACACATTCATT
                720 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
lcl|1_          720 CTCTTTCGACGCTCTTGCTAAAGCCCATTCAATTCGCCCCATATTTCGCACACATTCATT

gnl|BL_OR       780 TCACCACCCAATACGTG 797
                780 ||||||||||||||||| 797
lcl|1_          780 TCACCACCCAATACGTG 797
""",
        )


class TestBlastx(unittest.TestCase):
    """Test the Blast XML parser for blastx output."""

    def test_xml_2222_blastx_001(self):
        """Parsing BLASTX 2.2.22+ (xml_2222_blastx_001.xml)."""
        filename = "xml_2222_blastx_001.xml"
        datafile = os.path.join("Blast", filename)
        with open(datafile, "rb") as handle:
            records = Blast.parse(handle)
            self.assertEqual(records.program, "blastx")
            self.assertEqual(records.version, "BLASTX 2.2.22+")
            self.assertEqual(
                records.reference,
                'Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.',
            )
            self.assertEqual(records.db, "nr")
            self.assertIsInstance(records.query, SeqRecord)
            self.assertEqual(records.query.id, "1")
            self.assertEqual(
                records.query.description,
                "gi|4104054|gb|AH007193.1|SEG_CVIGS Centaurea vallesiaca 18S ribosomal RNA gene, partial sequence",
            )
            self.assertEqual(repr(records.query.seq), "Seq(None, length=1002)")
            self.assertEqual(len(records.param), 5)
            self.assertEqual(records.param["matrix"], "BLOSUM62")
            self.assertAlmostEqual(records.param["expect"], 0.0001)
            self.assertEqual(records.param["gap-open"], 11)
            self.assertEqual(records.param["gap-extend"], 1)
            self.assertEqual(records.param["filter"], "L;")
            record = next(records)
            self.assertIsInstance(record.query, SeqRecord)
            self.assertEqual(record.query.id, "1")
            self.assertEqual(
                record.query.description,
                "gi|4104054|gb|AH007193.1|SEG_CVIGS Centaurea vallesiaca 18S ribosomal RNA gene, partial sequence",
            )
            self.assertEqual(repr(record.query.seq), "Seq(None, length=1002)")

            self.assertEqual(len(record.stat), 7)
            self.assertEqual(record.stat["db-num"], 8994603)
            self.assertEqual(record.stat["db-len"], -1216159329)
            self.assertEqual(record.stat["hsp-len"], 0)
            self.assertAlmostEqual(record.stat["eff-space"], 367397307882.0)
            self.assertAlmostEqual(record.stat["kappa"], 0.041)
            self.assertAlmostEqual(record.stat["lambda"], 0.267)
            self.assertAlmostEqual(record.stat["entropy"], 0.14)
            self.assertEqual(len(record), 1)
            hit = record[0]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|149390769|gb|ABR25402.1|")
            self.assertEqual(hit.target.name, "ABR25402")
            self.assertEqual(
                hit.target.description, "unknown [Oryza sativa (indica cultivar-group)]"
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=26)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 129.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 54.2989775733826)
            self.assertAlmostEqual(hsp.annotations["evalue"], 1.83262460293058e-05)
            self.assertEqual(hsp.annotations["identity"], 24)
            self.assertEqual(hsp.annotations["positive"], 25)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[ 0, 26],
                              [ 0, 26]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 26))
            self.assertEqual(repr(hsp.query.seq), "Seq('HMLVSKIKPCMCKYEQIQTVKLRMAH')")
            self.assertEqual(hsp.query.id, "1")
            self.assertEqual(
                hsp.query.description,
                "gi|4104054|gb|AH007193.1|SEG_CVIGS Centaurea vallesiaca 18S ribosomal RNA gene, partial sequence",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(26))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "1:911..988")
            self.assertEqual(repr(hsp.target.seq), "Seq('HMLVSKIKPCMCKYELIRTVKLRMAH')")
            self.assertEqual(hsp.target.id, "gi|149390769|gb|ABR25402.1|")
            self.assertEqual(hsp.target.name, "ABR25402")
            self.assertEqual(
                hsp.target.description, "unknown [Oryza sativa (indica cultivar-group)]"
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(hsp.annotations["midline"], "HMLVSKIKPCMCKYE I+TVKLRMAH")
            self.assertEqual(
                str(hsp),
                """\
gi|149390         0 HMLVSKIKPCMCKYELIRTVKLRMAH 26
                  0 |||||||||||||||.|.|||||||| 26
1                 0 HMLVSKIKPCMCKYEQIQTVKLRMAH 26
""",
            )
            record = next(records)
            self.assertIsInstance(record.query, SeqRecord)
            self.assertEqual(record.query.id, "2")
            self.assertEqual(
                record.query.description,
                "gi|4218935|gb|AF074388.1|AF074388 Sambucus nigra hevein-like protein HLPf gene, partial cds",
            )
            self.assertEqual(repr(record.query.seq), "Seq(None, length=2050)")

            self.assertEqual(len(record.stat), 7)
            self.assertEqual(record.stat["db-num"], 8994603)
            self.assertEqual(record.stat["db-len"], -1216159329)
            self.assertEqual(record.stat["hsp-len"], 0)
            self.assertAlmostEqual(record.stat["eff-space"], 967993058520.0)
            self.assertAlmostEqual(record.stat["kappa"], 0.041)
            self.assertAlmostEqual(record.stat["lambda"], 0.267)
            self.assertAlmostEqual(record.stat["entropy"], 0.14)
            self.assertEqual(len(record), 10)
            hit = record[0]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|4218936|gb|AAD12237.1|")
            self.assertEqual(hit.target.name, "AAD12237")
            self.assertEqual(
                hit.target.description, "hevein-like protein HLPf [Sambucus nigra]"
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=333)")
            self.assertEqual(len(hit), 2)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 1053.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 410.223385721017)
            self.assertAlmostEqual(hsp.annotations["evalue"], 3.48406066731465e-112)
            self.assertEqual(hsp.annotations["identity"], 199)
            self.assertEqual(hsp.annotations["positive"], 200)
            self.assertEqual(hsp.annotations["gaps"], 33)
            hsp = hit[1]
            self.assertAlmostEqual(hsp.score, 683.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 267.699542631596)
            self.assertAlmostEqual(hsp.annotations["evalue"], 2.79278546744412e-69)
            self.assertEqual(hsp.annotations["identity"], 127)
            self.assertEqual(hsp.annotations["positive"], 127)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[206, 333],
                              [  0, 127]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 127))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('NYNYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINA...SVV')",
            )
            self.assertEqual(hsp.query.id, "2")
            self.assertEqual(
                hsp.query.description,
                "gi|4218935|gb|AF074388.1|AF074388 Sambucus nigra hevein-like protein HLPf gene, partial cds",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(127))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "2:1669..2049")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq({206: 'NYNYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINA...SVV'}, length=333)",
            )
            self.assertEqual(hsp.target.id, "gi|4218936|gb|AAD12237.1|")
            self.assertEqual(hsp.target.name, "AAD12237")
            self.assertEqual(
                hsp.target.description, "hevein-like protein HLPf [Sambucus nigra]"
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "NYNYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINANSEASDQVPSYGVISKIINSNFGHQSGLDTITTSIGYYKRYCDMLEVSYGDNLENWFDETPFTKVAHIKMSVV",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|421893       206 NYNYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINANSEASD
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
2                 0 NYNYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINANSEASD

gi|421893       266 QVPSYGVISKIINSNFGHQSGLDTITTSIGYYKRYCDMLEVSYGDNLENWFDETPFTKVA
                 60 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
2                60 QVPSYGVISKIINSNFGHQSGLDTITTSIGYYKRYCDMLEVSYGDNLENWFDETPFTKVA

gi|421893       326 HIKMSVV 333
                120 ||||||| 127
2               120 HIKMSVV 127
""",
            )
            hit = record[1]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|4206074|gb|AAD11408.1|")
            self.assertEqual(hit.target.name, "AAD11408")
            self.assertEqual(
                hit.target.description, "hevein-like protein [Sambucus nigra]"
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=333)")
            self.assertEqual(len(hit), 2)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 1043.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 406.371389961843)
            self.assertAlmostEqual(hsp.annotations["evalue"], 5.03097287018806e-111)
            self.assertEqual(hsp.annotations["identity"], 198)
            self.assertEqual(hsp.annotations["positive"], 199)
            self.assertEqual(hsp.annotations["gaps"], 33)
            hsp = hit[1]
            self.assertAlmostEqual(hsp.score, 672.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 263.462347296505)
            self.assertAlmostEqual(hsp.annotations["evalue"], 5.26696544712228e-68)
            self.assertEqual(hsp.annotations["identity"], 125)
            self.assertEqual(hsp.annotations["positive"], 126)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[206, 333],
                              [  0, 127]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 127))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('NYNYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINA...SVV')",
            )
            self.assertEqual(hsp.query.id, "2")
            self.assertEqual(
                hsp.query.description,
                "gi|4218935|gb|AF074388.1|AF074388 Sambucus nigra hevein-like protein HLPf gene, partial cds",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(127))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "2:1669..2049")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq({206: 'NYYYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINA...SLV'}, length=333)",
            )
            self.assertEqual(hsp.target.id, "gi|4206074|gb|AAD11408.1|")
            self.assertEqual(hsp.target.name, "AAD11408")
            self.assertEqual(
                hsp.target.description, "hevein-like protein [Sambucus nigra]"
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "NY YGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINANSEASDQVPSYGVISKIINSNFGHQSGLDTITTSIGYYKRYCDMLEVSYGDNLENWFDETPFTKVAHIKMS+V",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|420607       206 NYYYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINANSEASD
                  0 ||.|||||||||||||||||||||||||||||||||||||||||||||||||||||||||
2                 0 NYNYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINANSEASD

gi|420607       266 QVPSYGVISKIINSNFGHQSGLDTITTSIGYYKRYCDMLEVSYGDNLENWFDETPFTKVA
                 60 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
2                60 QVPSYGVISKIINSNFGHQSGLDTITTSIGYYKRYCDMLEVSYGDNLENWFDETPFTKVA

gi|420607       326 HIKMSLV 333
                120 |||||.| 127
2               120 HIKMSVV 127
""",
            )
            hit = record[2]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|4206070|gb|AAD11406.1|")
            self.assertEqual(hit.target.name, "AAD11406")
            self.assertEqual(
                hit.target.description, "hevein-like protein [Sambucus nigra]"
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=333)")
            self.assertEqual(len(hit), 2)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 1043.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 406.371389961843)
            self.assertAlmostEqual(hsp.annotations["evalue"], 5.03097287018806e-111)
            self.assertEqual(hsp.annotations["identity"], 198)
            self.assertEqual(hsp.annotations["positive"], 199)
            self.assertEqual(hsp.annotations["gaps"], 33)
            hsp = hit[1]
            self.assertAlmostEqual(hsp.score, 680.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 266.543943903844)
            self.assertAlmostEqual(hsp.annotations["evalue"], 6.22167692942359e-69)
            self.assertEqual(hsp.annotations["identity"], 126)
            self.assertEqual(hsp.annotations["positive"], 127)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[206, 333],
                              [  0, 127]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 127))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('NYNYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINA...SVV')",
            )
            self.assertEqual(hsp.query.id, "2")
            self.assertEqual(
                hsp.query.description,
                "gi|4218935|gb|AF074388.1|AF074388 Sambucus nigra hevein-like protein HLPf gene, partial cds",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(127))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "2:1669..2049")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq({206: 'NYNYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINA...SLV'}, length=333)",
            )
            self.assertEqual(hsp.target.id, "gi|4206070|gb|AAD11406.1|")
            self.assertEqual(hsp.target.name, "AAD11406")
            self.assertEqual(
                hsp.target.description, "hevein-like protein [Sambucus nigra]"
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "NYNYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINANSEASDQVPSYGVISKIINSNFGHQSGLDTITTSIGYYKRYCDMLEVSYGDNLENWFDETPFTKVAHIKMS+V",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|420607       206 NYNYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINANSEASD
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
2                 0 NYNYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINANSEASD

gi|420607       266 QVPSYGVISKIINSNFGHQSGLDTITTSIGYYKRYCDMLEVSYGDNLENWFDETPFTKVA
                 60 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
2                60 QVPSYGVISKIINSNFGHQSGLDTITTSIGYYKRYCDMLEVSYGDNLENWFDETPFTKVA

gi|420607       326 HIKMSLV 333
                120 |||||.| 127
2               120 HIKMSVV 127
""",
            )
            hit = record[3]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|4206072|gb|AAD11407.1|")
            self.assertEqual(hit.target.name, "AAD11407")
            self.assertEqual(
                hit.target.description, "hevein-like protein [Sambucus nigra]"
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=333)")
            self.assertEqual(len(hit), 2)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 1016.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 395.971001412075)
            self.assertAlmostEqual(hsp.annotations["evalue"], 6.7995613312017e-108)
            self.assertEqual(hsp.annotations["identity"], 193)
            self.assertEqual(hsp.annotations["positive"], 195)
            self.assertEqual(hsp.annotations["gaps"], 33)
            hsp = hit[1]
            self.assertAlmostEqual(hsp.score, 646.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 253.447158322654)
            self.assertAlmostEqual(hsp.annotations["evalue"], 5.45045505347399e-65)
            self.assertEqual(hsp.annotations["identity"], 120)
            self.assertEqual(hsp.annotations["positive"], 124)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[206, 333],
                              [  0, 127]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 127))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('NYNYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINA...SVV')",
            )
            self.assertEqual(hsp.query.id, "2")
            self.assertEqual(
                hsp.query.description,
                "gi|4218935|gb|AF074388.1|AF074388 Sambucus nigra hevein-like protein HLPf gene, partial cds",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(127))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "2:1669..2049")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq({206: 'NYYYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINA...SVV'}, length=333)",
            )
            self.assertEqual(hsp.target.id, "gi|4206072|gb|AAD11407.1|")
            self.assertEqual(hsp.target.name, "AAD11407")
            self.assertEqual(
                hsp.target.description, "hevein-like protein [Sambucus nigra]"
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "NY YGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINANSEASDQVPSYGVIS+II+SN GHQS LDTITTSIGYYKRYCDMLEVSYGDNL+NWFDETPF+KVAHIKMSVV",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|420607       206 NYYYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINANSEASD
                  0 ||.|||||||||||||||||||||||||||||||||||||||||||||||||||||||||
2                 0 NYNYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINANSEASD

gi|420607       266 QVPSYGVISEIIDSNIGHQSSLDTITTSIGYYKRYCDMLEVSYGDNLKNWFDETPFSKVA
                 60 |||||||||.||.||.||||.||||||||||||||||||||||||||.||||||||.|||
2                60 QVPSYGVISKIINSNFGHQSGLDTITTSIGYYKRYCDMLEVSYGDNLENWFDETPFTKVA

gi|420607       326 HIKMSVV 333
                120 ||||||| 127
2               120 HIKMSVV 127
""",
            )
            hit = record[4]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|16903131|gb|AAL30421.1|AF434174_1")
            self.assertEqual(hit.target.name, "AAL30421")
            self.assertEqual(
                hit.target.description, "hevein-like protein [Sambucus nigra]"
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=330)")
            self.assertEqual(len(hit), 2)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 986.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 384.415014134554)
            self.assertAlmostEqual(hsp.annotations["evalue"], 2.04729155722083e-104)
            self.assertEqual(hsp.annotations["identity"], 190)
            self.assertEqual(hsp.annotations["positive"], 191)
            self.assertEqual(hsp.annotations["gaps"], 36)
            hsp = hit[1]
            self.assertAlmostEqual(hsp.score, 679.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 266.158744327927)
            self.assertAlmostEqual(hsp.annotations["evalue"], 8.12576171382949e-69)
            self.assertEqual(hsp.annotations["identity"], 126)
            self.assertEqual(hsp.annotations["positive"], 126)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[203, 330],
                              [  0, 127]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 127))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('NYNYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINA...SVV')",
            )
            self.assertEqual(hsp.query.id, "2")
            self.assertEqual(
                hsp.query.description,
                "gi|4218935|gb|AF074388.1|AF074388 Sambucus nigra hevein-like protein HLPf gene, partial cds",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(127))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "2:1669..2049")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq({203: 'NYNYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINA...SVV'}, length=330)",
            )
            self.assertEqual(hsp.target.id, "gi|16903131|gb|AAL30421.1|AF434174_1")
            self.assertEqual(hsp.target.name, "AAL30421")
            self.assertEqual(
                hsp.target.description, "hevein-like protein [Sambucus nigra]"
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "NYNYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINANSEASDQVPSYGVISKIINSNFGHQSGLDTITTSIGYYKRYCDMLEVSYGDNLENWFDETPFTKV HIKMSVV",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|169031       203 NYNYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINANSEASD
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
2                 0 NYNYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINANSEASD

gi|169031       263 QVPSYGVISKIINSNFGHQSGLDTITTSIGYYKRYCDMLEVSYGDNLENWFDETPFTKVV
                 60 |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||.
2                60 QVPSYGVISKIINSNFGHQSGLDTITTSIGYYKRYCDMLEVSYGDNLENWFDETPFTKVA

gi|169031       323 HIKMSVV 330
                120 ||||||| 127
2               120 HIKMSVV 127
""",
            )
            hit = record[5]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|16903133|gb|AAL30422.1|AF434175_1")
            self.assertEqual(hit.target.name, "AAL30422")
            self.assertEqual(
                hit.target.description, "hevein-like protein [Sambucus nigra]"
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=336)")
            self.assertEqual(len(hit), 2)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 713.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 279.255529909117)
            self.assertAlmostEqual(hsp.annotations["evalue"], 9.27553088557319e-73)
            self.assertEqual(hsp.annotations["identity"], 148)
            self.assertEqual(hsp.annotations["positive"], 162)
            self.assertEqual(hsp.annotations["gaps"], 40)
            hsp = hit[1]
            self.assertAlmostEqual(hsp.score, 620.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 243.431969348803)
            self.assertAlmostEqual(hsp.annotations["evalue"], 5.64033703812707e-62)
            self.assertEqual(hsp.annotations["identity"], 115)
            self.assertEqual(hsp.annotations["positive"], 120)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[209, 336],
                              [  0, 127]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 127))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('NYNYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINA...SVV')",
            )
            self.assertEqual(hsp.query.id, "2")
            self.assertEqual(
                hsp.query.description,
                "gi|4218935|gb|AF074388.1|AF074388 Sambucus nigra hevein-like protein HLPf gene, partial cds",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(127))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "2:1669..2049")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq({209: 'NYNYGLAGEAIGIDLVNDPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINA...HVV'}, length=336)",
            )
            self.assertEqual(hsp.target.id, "gi|16903133|gb|AAL30422.1|AF434175_1")
            self.assertEqual(hsp.target.name, "AAL30422")
            self.assertEqual(
                hsp.target.description, "hevein-like protein [Sambucus nigra]"
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "NYNYGLAGEA+GIDLVN PDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINAN EASDQVPSYGV+S IINSN GH+SGLD ITTSIGYYKRYCDMLEVSYGDNL+NWFDETPF+KVA IKM VV",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|169031       209 NYNYGLAGEAIGIDLVNDPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINANYEASD
                  0 ||||||||||.||||||.|||||||||||||||||||||||||||||||||||||.||||
2                 0 NYNYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINANSEASD

gi|169031       269 QVPSYGVLSNIINSNSGHKSGLDIITTSIGYYKRYCDMLEVSYGDNLKNWFDETPFSKVA
                 60 |||||||.|.|||||.||.||||.|||||||||||||||||||||||.||||||||.|||
2                60 QVPSYGVISKIINSNFGHQSGLDTITTSIGYYKRYCDMLEVSYGDNLENWFDETPFTKVA

gi|169031       329 RIKMHVV 336
                120 .|||.|| 127
2               120 HIKMSVV 127
""",
            )
            hit = record[6]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|30691147|gb|AAO17294.1|")
            self.assertEqual(hit.target.name, "AAO17294")
            self.assertEqual(hit.target.description, "chitinase [Ficus carica]")
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=321)")
            self.assertEqual(len(hit), 2)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 481.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 189.889228296291)
            self.assertAlmostEqual(hsp.annotations["evalue"], 7.40075731140555e-46)
            self.assertEqual(hsp.annotations["identity"], 113)
            self.assertEqual(hsp.annotations["positive"], 138)
            self.assertEqual(hsp.annotations["gaps"], 49)
            hsp = hit[1]
            self.assertAlmostEqual(hsp.score, 426.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 168.703251620836)
            self.assertAlmostEqual(hsp.annotations["evalue"], 1.76559312729927e-39)
            self.assertEqual(hsp.annotations["identity"], 81)
            self.assertEqual(hsp.annotations["positive"], 99)
            self.assertEqual(hsp.annotations["gaps"], 10)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[202, 261, 261, 266, 266, 308, 309, 320],
                              [  0,  59,  60,  65,  73, 115, 115, 126]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 127))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('NYNYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINA...MSV')",
            )
            self.assertEqual(hsp.query.id, "2")
            self.assertEqual(
                hsp.query.description,
                "gi|4218935|gb|AF074388.1|AF074388 Sambucus nigra hevein-like protein HLPf gene, partial cds",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(126))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "2:1669..2046")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq({202: 'NHNYGLVGEALGIDLINNPDLVATDPVVSFKTAIWFWMTRHQNKPSFHGVIINA...MPV'}, length=321)",
            )
            self.assertEqual(hsp.target.id, "gi|30691147|gb|AAO17294.1|")
            self.assertEqual(hsp.target.name, "AAO17294")
            self.assertEqual(hsp.target.description, "chitinase [Ficus carica]")
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "N+NYGL GEALGIDL+N+PDLVATDP+VSFKTAIWFWMT+H N PS H ++INANSE S  +P++        SNFG +S LD +  SIGYYKRYCDML+VS+GDNL+ W+D TP F+ V+ I M V",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|306911       202 NHNYGLVGEALGIDLINNPDLVATDPVVSFKTAIWFWMTRHQNKPSFHGVIINANSEPS-
                  0 |.||||.||||||||.|.||||||||.||||||||||||.|.|.||.|...||||||.|-
2                 0 NYNYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINANSEASD

gi|306911       261 HIPNH--------SNFGQESVLDVVNRSIGYYKRYCDMLKVSFGDNLKYWYDGTPNFSDV
                 60 ..|..--------||||..|.||....||||||||||||.||.||||..|.|.||-|..|
2                60 QVPSYGVISKIINSNFGHQSGLDTITTSIGYYKRYCDMLEVSYGDNLENWFDETP-FTKV

gi|306911       313 SRIGMPV 320
                120 ..|.|.| 127
2               119 AHIKMSV 126
""",
            )
            hit = record[7]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|222139388|gb|ACM45713.1|")
            self.assertEqual(hit.target.name, "ACM45713")
            self.assertEqual(
                hit.target.description, "class I chitinase [Pyrus pyrifolia]"
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=317)")
            self.assertEqual(len(hit), 2)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 469.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 185.266833385283)
            self.assertAlmostEqual(hsp.annotations["evalue"], 1.82286993845418e-44)
            self.assertEqual(hsp.annotations["identity"], 111)
            self.assertEqual(hsp.annotations["positive"], 137)
            self.assertEqual(hsp.annotations["gaps"], 50)
            hsp = hit[1]
            self.assertAlmostEqual(hsp.score, 318.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 127.101697421762)
            self.assertAlmostEqual(hsp.annotations["evalue"], 5.89123449548921e-27)
            self.assertEqual(hsp.annotations["identity"], 62)
            self.assertEqual(hsp.annotations["positive"], 84)
            self.assertEqual(hsp.annotations["gaps"], 7)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[195, 247, 252, 284, 285, 309, 309, 316],
                              [  0,  52,  52,  84,  84, 108, 109, 116]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 122))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('NYNYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINA...TPF')",
            )
            self.assertEqual(hsp.query.id, "2")
            self.assertEqual(
                hsp.query.description,
                "gi|4218935|gb|AF074388.1|AF074388 Sambucus nigra hevein-like protein HLPf gene, partial cds",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(116))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "2:1669..2016")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq({195: 'NYNYGQAGKAIGKDLINNPDLVATDPVVSFKTAIWFWMTPQGNKPSSHDVITGR...RPF'}, length=317)",
            )
            self.assertEqual(hsp.target.id, "gi|222139388|gb|ACM45713.1|")
            self.assertEqual(hsp.target.name, "ACM45713")
            self.assertEqual(
                hsp.target.description, "class I chitinase [Pyrus pyrifolia]"
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "NYNYG AG+A+G DL+N+PDLVATDP+VSFKTAIWFWMT   N PS HD++      +    ++ +VP YGVI+ IIN       G D  + + IG+Y+RYC +L V+ GDNL+  +++ PF",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|222139       195 NYNYGQAGKAIGKDLINNPDLVATDPVVSFKTAIWFWMTPQGNKPSSHDVITGRWSPSTA
                  0 |||||.||.|.|.||.|.||||||||.||||||||||||...|.||.||...-----...
2                 0 NYNYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILI-----NAN

gi|222139       255 DRSAGRVPGYGVITNIINGGVECGKGQDARVASRIGFYRRYCQILGVNPGDNLD-CYNQR
                 60 ......||.||||..|||.......|.|.-....||.|.|||..|.|..||||.-.....
2                55 SEASDQVPSYGVISKIINSNFGHQSGLDT-ITTSIGYYKRYCDMLEVSYGDNLENWFDET

gi|222139       314 PF 316
                120 || 122
2               114 PF 116
""",
            )
            hit = record[8]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|23496435|dbj|BAB40817.2|")
            self.assertEqual(hit.target.name, "BAB40817")
            self.assertEqual(
                hit.target.description, "endochitinase MCHT-2 [Cucumis melo]"
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=311)")
            self.assertEqual(len(hit), 2)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 460.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 181.800037202027)
            self.assertAlmostEqual(hsp.annotations["evalue"], 2.01541888137674e-43)
            self.assertEqual(hsp.annotations["identity"], 109)
            self.assertEqual(hsp.annotations["positive"], 132)
            self.assertEqual(hsp.annotations["gaps"], 54)
            hsp = hit[1]
            self.assertAlmostEqual(hsp.score, 285.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 114.39011141649)
            self.assertAlmostEqual(hsp.annotations["evalue"], 3.95161831690075e-23)
            self.assertEqual(hsp.annotations["identity"], 56)
            self.assertEqual(hsp.annotations["positive"], 75)
            self.assertEqual(hsp.annotations["gaps"], 7)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[191, 211, 211, 242, 247, 277, 278, 304],
                              [  0,  20,  21,  52,  52,  82,  82, 108]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 114))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('NYNYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINA...NLE')",
            )
            self.assertEqual(hsp.query.id, "2")
            self.assertEqual(
                hsp.query.description,
                "gi|4218935|gb|AF074388.1|AF074388 Sambucus nigra hevein-like protein HLPf gene, partial cds",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(108))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "2:1669..1992")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq({191: 'NYNYGPAGKAIGAPLLTNPDTATDPVTSFKTALWFWMTAQGNKPSCHNVITGNW...NLD'}, length=311)",
            )
            self.assertEqual(hsp.target.id, "gi|23496435|dbj|BAB40817.2|")
            self.assertEqual(hsp.target.name, "BAB40817")
            self.assertEqual(
                hsp.target.description, "endochitinase MCHT-2 [Cucumis melo]"
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "NYNYG AG+A+G  L+ +PD  ATDP+ SFKTA+WFWMT   N PS H+++      ++   A+ +VP YGVI+ IIN       G  D +   IG+YKRYCDML + YG+NL+",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|234964       191 NYNYGPAGKAIGAPLLTNPD-TATDPVTSFKTALWFWMTAQGNKPSCHNVITGNWQPSSA
                  0 |||||.||.|.|..|...||-.||||..|||||.|||||...|.||.|....-----...
2                 0 NYNYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILI-----NAN

gi|234964       250 DNAAGRVPGYGVITNIINGGLECGRGPDDRVKDRIGFYKRYCDMLGIGYGNNLD 304
                 60 ..|...||.||||..|||.......|.-|.....||.||||||||...||.||. 114
2                55 SEASDQVPSYGVISKIINSNFGHQSGL-DTITTSIGYYKRYCDMLEVSYGDNLE 108
""",
            )
            hit = record[9]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|82621253|gb|ABB86300.1|")
            self.assertEqual(hit.target.name, "ABB86300")
            self.assertEqual(hit.target.description, "chitinase [Ficus awkeotsang]")
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=301)")
            self.assertEqual(len(hit), 2)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 459.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 181.414837626109)
            self.assertAlmostEqual(hsp.annotations["evalue"], 2.6322185753765e-43)
            self.assertEqual(hsp.annotations["identity"], 114)
            self.assertEqual(hsp.annotations["positive"], 134)
            self.assertEqual(hsp.annotations["gaps"], 50)
            hsp = hit[1]
            self.assertAlmostEqual(hsp.score, 359.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 142.894880034374)
            self.assertAlmostEqual(hsp.annotations["evalue"], 1.03749166509001e-31)
            self.assertEqual(hsp.annotations["identity"], 67)
            self.assertEqual(hsp.annotations["positive"], 83)
            self.assertEqual(hsp.annotations["gaps"], 9)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[203, 263, 263, 268, 268, 301],
                              [  0,  60,  61,  66,  74, 107]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 107))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('NYNYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINA...DNL')",
            )
            self.assertEqual(hsp.query.id, "2")
            self.assertEqual(
                hsp.query.description,
                "gi|4218935|gb|AF074388.1|AF074388 Sambucus nigra hevein-like protein HLPf gene, partial cds",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(107))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "2:1669..1989")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq({203: 'NHNYGLVGEALGIDLINNPELVATDPVISFKTAIWFWMARYEDKPSFHDVIINA...DNL'}, length=301)",
            )
            self.assertEqual(hsp.target.id, "gi|82621253|gb|ABB86300.1|")
            self.assertEqual(hsp.target.name, "ABB86300")
            self.assertEqual(hsp.target.description, "chitinase [Ficus awkeotsang]")
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "N+NYGL GEALGIDL+N+P+LVATDP++SFKTAIWFWM ++++ PS HD++INAN EASD +P +G        N G +S LD +  SIGYYKRYCDML VS  DNL",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|826212       203 NHNYGLVGEALGIDLINNPELVATDPVISFKTAIWFWMARYEDKPSFHDVIINANFEASD
                  0 |.||||.||||||||.|.|.||||||..||||||||||......||.||..||||.||||
2                 0 NYNYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINANSEASD

gi|826212       263 -IPYHG--------NSGQESSLDVVNRSIGYYKRYCDMLGVSCEDNL 301
                 60 -.|..|--------|.|..|.||....||||||||||||.||..||| 107
2                60 QVPSYGVISKIINSNFGHQSGLDTITTSIGYYKRYCDMLEVSYGDNL 107
""",
            )
            record = next(records)
            self.assertIsInstance(record.query, SeqRecord)
            self.assertEqual(record.query.id, "3")
            self.assertEqual(
                record.query.description,
                "gi|5690369|gb|AF158246.1|AF158246 Cricetulus griseus glucose phosphate isomerase (GPI) gene, partial intron sequence",
            )
            self.assertEqual(repr(record.query.seq), "Seq(None, length=550)")

            self.assertEqual(len(record.stat), 7)
            self.assertEqual(record.stat["db-num"], 8994603)
            self.assertEqual(record.stat["db-len"], -1216159329)
            self.assertEqual(record.stat["hsp-len"], 0)
            self.assertAlmostEqual(record.stat["eff-space"], 108443629616.0)
            self.assertAlmostEqual(record.stat["kappa"], 0.041)
            self.assertAlmostEqual(record.stat["lambda"], 0.267)
            self.assertAlmostEqual(record.stat["entropy"], 0.14)
            self.assertEqual(len(record), 0)
            record = next(records)
            self.assertIsInstance(record.query, SeqRecord)
            self.assertEqual(record.query.id, "4")
            self.assertEqual(
                record.query.description,
                "gi|5049839|gb|AI730987.1|AI730987 BNLGHi8354 Six-day Cotton fiber Gossypium hirsutum cDNA 5' similar to TUBULIN BETA-1 CHAIN gi|486734|pir|S35142 tubulin beta chain - white lupine gi|402636 (X70184) Beta tubulin 1 [Lupinus albus], mRNA sequence",
            )
            self.assertEqual(repr(record.query.seq), "Seq(None, length=655)")

            self.assertEqual(len(record.stat), 7)
            self.assertEqual(record.stat["db-num"], 8994603)
            self.assertEqual(record.stat["db-len"], -1216159329)
            self.assertEqual(record.stat["hsp-len"], 0)
            self.assertAlmostEqual(record.stat["eff-space"], 165344802738.0)
            self.assertAlmostEqual(record.stat["kappa"], 0.041)
            self.assertAlmostEqual(record.stat["lambda"], 0.267)
            self.assertAlmostEqual(record.stat["entropy"], 0.14)
            self.assertEqual(len(record), 10)
            hit = record[0]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|166343825|gb|ABY86655.1|")
            self.assertEqual(hit.target.name, "ABY86655")
            self.assertEqual(
                hit.target.description, "beta-tubulin 4 [Gossypium hirsutum]"
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=448)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 1048.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 408.29738784143)
            self.assertAlmostEqual(hsp.annotations["evalue"], 2.26145081918239e-112)
            self.assertEqual(hsp.annotations["identity"], 196)
            self.assertEqual(hsp.annotations["positive"], 197)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[  0, 201],
                              [  0, 201]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 201))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('MREILHIQAGQCGNQIGANFWEVVCAEHGINSTGRYQGDNDLQLERVNVYYNEA...CMV')",
            )
            self.assertEqual(hsp.query.id, "4")
            self.assertEqual(
                hsp.query.description,
                "gi|5049839|gb|AI730987.1|AI730987 BNLGHi8354 Six-day Cotton fiber Gossypium hirsutum cDNA 5' similar to TUBULIN BETA-1 CHAIN gi|486734|pir|S35142 tubulin beta chain - white lupine gi|402636 (X70184) Beta tubulin 1 [Lupinus albus], mRNA sequence",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(201))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "4:50..652")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq({0: 'MREILHIQAGQCGNQIGAKFWEVVCAEHGIDSTGRYQGDNDLQLERVNVYYNEA...CMV'}, length=448)",
            )
            self.assertEqual(hsp.target.id, "gi|166343825|gb|ABY86655.1|")
            self.assertEqual(hsp.target.name, "ABY86655")
            self.assertEqual(
                hsp.target.description, "beta-tubulin 4 [Gossypium hirsutum]"
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "MREILHIQAGQCGNQIGA FWEVVCAEHGI+STGRYQGDNDLQLERVNVYYNEASCGRFVPRAVLMDLEPGTMDSVRSGPYGQIFRPDNFVFGQSGAGNNWAKGHYTEGAELIDS LDVVRKEAENCDCLQGFQVCHSLG GTGSGMGTLLISKIREEYPDRMMLTFSVFPSPKVSDTVVEPYNATLSVH LVENADECMV",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|166343         0 MREILHIQAGQCGNQIGAKFWEVVCAEHGIDSTGRYQGDNDLQLERVNVYYNEASCGRFV
                  0 ||||||||||||||||||.|||||||||||.|||||||||||||||||||||||||||||
4                 0 MREILHIQAGQCGNQIGANFWEVVCAEHGINSTGRYQGDNDLQLERVNVYYNEASCGRFV

gi|166343        60 PRAVLMDLEPGTMDSVRSGPYGQIFRPDNFVFGQSGAGNNWAKGHYTEGAELIDSVLDVV
                 60 |||||||||||||||||||||||||||||||||||||||||||||||||||||||.||||
4                60 PRAVLMDLEPGTMDSVRSGPYGQIFRPDNFVFGQSGAGNNWAKGHYTEGAELIDSXLDVV

gi|166343       120 RKEAENCDCLQGFQVCHSLGGGTGSGMGTLLISKIREEYPDRMMLTFSVFPSPKVSDTVV
                120 ||||||||||||||||||||.|||||||||||||||||||||||||||||||||||||||
4               120 RKEAENCDCLQGFQVCHSLGRGTGSGMGTLLISKIREEYPDRMMLTFSVFPSPKVSDTVV

gi|166343       180 EPYNATLSVHQLVENADECMV 201
                180 ||||||||||.|||||||||| 201
4               180 EPYNATLSVHXLVENADECMV 201
""",
            )
            hit = record[1]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|223549899|gb|EEF51386.1|")
            self.assertEqual(hit.target.name, "EEF51386")
            self.assertEqual(
                hit.target.description,
                "tubulin beta chain, putative [Ricinus communis]",
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=448)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 1044.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 406.756589537761)
            self.assertAlmostEqual(hsp.annotations["evalue"], 6.57981456092236e-112)
            self.assertEqual(hsp.annotations["identity"], 195)
            self.assertEqual(hsp.annotations["positive"], 196)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[  0, 201],
                              [  0, 201]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 201))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('MREILHIQAGQCGNQIGANFWEVVCAEHGINSTGRYQGDNDLQLERVNVYYNEA...CMV')",
            )
            self.assertEqual(hsp.query.id, "4")
            self.assertEqual(
                hsp.query.description,
                "gi|5049839|gb|AI730987.1|AI730987 BNLGHi8354 Six-day Cotton fiber Gossypium hirsutum cDNA 5' similar to TUBULIN BETA-1 CHAIN gi|486734|pir|S35142 tubulin beta chain - white lupine gi|402636 (X70184) Beta tubulin 1 [Lupinus albus], mRNA sequence",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(201))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "4:50..652")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq({0: 'MREILHIQGGQCGNQIGAKFWEVVCAEHGIDSTGRYQGDNDLQLERVNVYYNEA...CMV'}, length=448)",
            )
            self.assertEqual(hsp.target.id, "gi|223549899|gb|EEF51386.1|")
            self.assertEqual(hsp.target.name, "EEF51386")
            self.assertEqual(
                hsp.target.description,
                "tubulin beta chain, putative [Ricinus communis]",
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "MREILHIQ GQCGNQIGA FWEVVCAEHGI+STGRYQGDNDLQLERVNVYYNEASCGRFVPRAVLMDLEPGTMDSVRSGPYGQIFRPDNFVFGQSGAGNNWAKGHYTEGAELIDS LDVVRKEAENCDCLQGFQVCHSLG GTGSGMGTLLISKIREEYPDRMMLTFSVFPSPKVSDTVVEPYNATLSVH LVENADECMV",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|223549         0 MREILHIQGGQCGNQIGAKFWEVVCAEHGIDSTGRYQGDNDLQLERVNVYYNEASCGRFV
                  0 ||||||||.|||||||||.|||||||||||.|||||||||||||||||||||||||||||
4                 0 MREILHIQAGQCGNQIGANFWEVVCAEHGINSTGRYQGDNDLQLERVNVYYNEASCGRFV

gi|223549        60 PRAVLMDLEPGTMDSVRSGPYGQIFRPDNFVFGQSGAGNNWAKGHYTEGAELIDSVLDVV
                 60 |||||||||||||||||||||||||||||||||||||||||||||||||||||||.||||
4                60 PRAVLMDLEPGTMDSVRSGPYGQIFRPDNFVFGQSGAGNNWAKGHYTEGAELIDSXLDVV

gi|223549       120 RKEAENCDCLQGFQVCHSLGGGTGSGMGTLLISKIREEYPDRMMLTFSVFPSPKVSDTVV
                120 ||||||||||||||||||||.|||||||||||||||||||||||||||||||||||||||
4               120 RKEAENCDCLQGFQVCHSLGRGTGSGMGTLLISKIREEYPDRMMLTFSVFPSPKVSDTVV

gi|223549       180 EPYNATLSVHQLVENADECMV 201
                180 ||||||||||.|||||||||| 201
4               180 EPYNATLSVHXLVENADECMV 201
""",
            )
            hit = record[2]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|18420724|ref|NP_568437.1|")
            self.assertEqual(hit.target.name, "NP_568437")
            self.assertEqual(
                hit.target.description,
                "TUB8 (tubulin beta-8) [Arabidopsis thaliana] >gi|27735261|sp|P29516.2|TBB8_ARATH RecName: Full=Tubulin beta-8 chain; AltName: Full=Beta-8-tubulin >gi|10176853|dbj|BAB10059.1| beta tubulin [Arabidopsis thaliana]",
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=449)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 1040.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 405.215791234091)
            self.assertAlmostEqual(hsp.annotations["evalue"], 1.91443295113426e-111)
            self.assertEqual(hsp.annotations["identity"], 194)
            self.assertEqual(hsp.annotations["positive"], 196)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[  0, 201],
                              [  0, 201]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 201))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('MREILHIQAGQCGNQIGANFWEVVCAEHGINSTGRYQGDNDLQLERVNVYYNEA...CMV')",
            )
            self.assertEqual(hsp.query.id, "4")
            self.assertEqual(
                hsp.query.description,
                "gi|5049839|gb|AI730987.1|AI730987 BNLGHi8354 Six-day Cotton fiber Gossypium hirsutum cDNA 5' similar to TUBULIN BETA-1 CHAIN gi|486734|pir|S35142 tubulin beta chain - white lupine gi|402636 (X70184) Beta tubulin 1 [Lupinus albus], mRNA sequence",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(201))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "4:50..652")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq({0: 'MREILHIQGGQCGNQIGAKFWEVVCAEHGIDSTGRYQGENDLQLERVNVYYNEA...CMV'}, length=449)",
            )
            self.assertEqual(hsp.target.id, "gi|18420724|ref|NP_568437.1|")
            self.assertEqual(hsp.target.name, "NP_568437")
            self.assertEqual(
                hsp.target.description,
                "TUB8 (tubulin beta-8) [Arabidopsis thaliana] >gi|27735261|sp|P29516.2|TBB8_ARATH RecName: Full=Tubulin beta-8 chain; AltName: Full=Beta-8-tubulin >gi|10176853|dbj|BAB10059.1| beta tubulin [Arabidopsis thaliana]",
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "MREILHIQ GQCGNQIGA FWEVVCAEHGI+STGRYQG+NDLQLERVNVYYNEASCGRFVPRAVLMDLEPGTMDSVRSGPYGQIFRPDNFVFGQSGAGNNWAKGHYTEGAELIDS LDVVRKEAENCDCLQGFQVCHSLG GTGSGMGTLLISKIREEYPDRMMLTFSVFPSPKVSDTVVEPYNATLSVH LVENADECMV",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|184207         0 MREILHIQGGQCGNQIGAKFWEVVCAEHGIDSTGRYQGENDLQLERVNVYYNEASCGRFV
                  0 ||||||||.|||||||||.|||||||||||.|||||||.|||||||||||||||||||||
4                 0 MREILHIQAGQCGNQIGANFWEVVCAEHGINSTGRYQGDNDLQLERVNVYYNEASCGRFV

gi|184207        60 PRAVLMDLEPGTMDSVRSGPYGQIFRPDNFVFGQSGAGNNWAKGHYTEGAELIDSVLDVV
                 60 |||||||||||||||||||||||||||||||||||||||||||||||||||||||.||||
4                60 PRAVLMDLEPGTMDSVRSGPYGQIFRPDNFVFGQSGAGNNWAKGHYTEGAELIDSXLDVV

gi|184207       120 RKEAENCDCLQGFQVCHSLGGGTGSGMGTLLISKIREEYPDRMMLTFSVFPSPKVSDTVV
                120 ||||||||||||||||||||.|||||||||||||||||||||||||||||||||||||||
4               120 RKEAENCDCLQGFQVCHSLGRGTGSGMGTLLISKIREEYPDRMMLTFSVFPSPKVSDTVV

gi|184207       180 EPYNATLSVHQLVENADECMV 201
                180 ||||||||||.|||||||||| 201
4               180 EPYNATLSVHXLVENADECMV 201
""",
            )
            hit = record[3]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|225426385|ref|XP_002271992.1|")
            self.assertEqual(hit.target.name, "XP_002271992")
            self.assertEqual(
                hit.target.description,
                "PREDICTED: hypothetical protein [Vitis vinifera] >gi|157356601|emb|CAO62796.1| unnamed protein product [Vitis vinifera]",
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=447)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 1034.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 402.904593778587)
            self.assertAlmostEqual(hsp.annotations["evalue"], 9.50123195540709e-111)
            self.assertEqual(hsp.annotations["identity"], 193)
            self.assertEqual(hsp.annotations["positive"], 195)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[  0, 201],
                              [  0, 201]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 201))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('MREILHIQAGQCGNQIGANFWEVVCAEHGINSTGRYQGDNDLQLERVNVYYNEA...CMV')",
            )
            self.assertEqual(hsp.query.id, "4")
            self.assertEqual(
                hsp.query.description,
                "gi|5049839|gb|AI730987.1|AI730987 BNLGHi8354 Six-day Cotton fiber Gossypium hirsutum cDNA 5' similar to TUBULIN BETA-1 CHAIN gi|486734|pir|S35142 tubulin beta chain - white lupine gi|402636 (X70184) Beta tubulin 1 [Lupinus albus], mRNA sequence",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(201))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "4:50..652")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq({0: 'MREILHIQGGQCGNQIGAKFWEVVCAEHGIDSTGRYHGDSDLQLERVNVYYNEA...CMV'}, length=447)",
            )
            self.assertEqual(hsp.target.id, "gi|225426385|ref|XP_002271992.1|")
            self.assertEqual(hsp.target.name, "XP_002271992")
            self.assertEqual(
                hsp.target.description,
                "PREDICTED: hypothetical protein [Vitis vinifera] >gi|157356601|emb|CAO62796.1| unnamed protein product [Vitis vinifera]",
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "MREILHIQ GQCGNQIGA FWEVVCAEHGI+STGRY GD+DLQLERVNVYYNEASCGRFVPRAVLMDLEPGTMDSVRSGPYGQIFRPDNFVFGQSGAGNNWAKGHYTEGAELIDS LDVVRKEAENCDCLQGFQVCHSLG GTGSGMGTLLISKIREEYPDRMMLTFSVFPSPKVSDTVVEPYNATLSVH LVENADECMV",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|225426         0 MREILHIQGGQCGNQIGAKFWEVVCAEHGIDSTGRYHGDSDLQLERVNVYYNEASCGRFV
                  0 ||||||||.|||||||||.|||||||||||.|||||.||.||||||||||||||||||||
4                 0 MREILHIQAGQCGNQIGANFWEVVCAEHGINSTGRYQGDNDLQLERVNVYYNEASCGRFV

gi|225426        60 PRAVLMDLEPGTMDSVRSGPYGQIFRPDNFVFGQSGAGNNWAKGHYTEGAELIDSVLDVV
                 60 |||||||||||||||||||||||||||||||||||||||||||||||||||||||.||||
4                60 PRAVLMDLEPGTMDSVRSGPYGQIFRPDNFVFGQSGAGNNWAKGHYTEGAELIDSXLDVV

gi|225426       120 RKEAENCDCLQGFQVCHSLGGGTGSGMGTLLISKIREEYPDRMMLTFSVFPSPKVSDTVV
                120 ||||||||||||||||||||.|||||||||||||||||||||||||||||||||||||||
4               120 RKEAENCDCLQGFQVCHSLGRGTGSGMGTLLISKIREEYPDRMMLTFSVFPSPKVSDTVV

gi|225426       180 EPYNATLSVHQLVENADECMV 201
                180 ||||||||||.|||||||||| 201
4               180 EPYNATLSVHXLVENADECMV 201
""",
            )
            hit = record[4]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|15451226|gb|AAK96884.1|")
            self.assertEqual(hit.target.name, "AAK96884")
            self.assertEqual(
                hit.target.description,
                "beta tubulin [Arabidopsis thaliana] >gi|20148289|gb|AAM10035.1| beta tubulin [Arabidopsis thaliana]",
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=449)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 1034.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 402.904593778587)
            self.assertAlmostEqual(hsp.annotations["evalue"], 9.50123195540709e-111)
            self.assertEqual(hsp.annotations["identity"], 193)
            self.assertEqual(hsp.annotations["positive"], 195)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[  0, 201],
                              [  0, 201]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 201))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('MREILHIQAGQCGNQIGANFWEVVCAEHGINSTGRYQGDNDLQLERVNVYYNEA...CMV')",
            )
            self.assertEqual(hsp.query.id, "4")
            self.assertEqual(
                hsp.query.description,
                "gi|5049839|gb|AI730987.1|AI730987 BNLGHi8354 Six-day Cotton fiber Gossypium hirsutum cDNA 5' similar to TUBULIN BETA-1 CHAIN gi|486734|pir|S35142 tubulin beta chain - white lupine gi|402636 (X70184) Beta tubulin 1 [Lupinus albus], mRNA sequence",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(201))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "4:50..652")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq({0: 'MREILHIQGGQCGNQIGAKFWEVVCAEHGIDSTGRYQGEKDLQLERVNVYYNEA...CMV'}, length=449)",
            )
            self.assertEqual(hsp.target.id, "gi|15451226|gb|AAK96884.1|")
            self.assertEqual(hsp.target.name, "AAK96884")
            self.assertEqual(
                hsp.target.description,
                "beta tubulin [Arabidopsis thaliana] >gi|20148289|gb|AAM10035.1| beta tubulin [Arabidopsis thaliana]",
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "MREILHIQ GQCGNQIGA FWEVVCAEHGI+STGRYQG+ DLQLERVNVYYNEASCGRFVPRAVLMDLEPGTMDSVRSGPYGQIFRPDNFVFGQSGAGNNWAKGHYTEGAELIDS LDVVRKEAENCDCLQGFQVCHSLG GTGSGMGTLLISKIREEYPDRMMLTFSVFPSPKVSDTVVEPYNATLSVH LVENADECMV",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|154512         0 MREILHIQGGQCGNQIGAKFWEVVCAEHGIDSTGRYQGEKDLQLERVNVYYNEASCGRFV
                  0 ||||||||.|||||||||.|||||||||||.|||||||..||||||||||||||||||||
4                 0 MREILHIQAGQCGNQIGANFWEVVCAEHGINSTGRYQGDNDLQLERVNVYYNEASCGRFV

gi|154512        60 PRAVLMDLEPGTMDSVRSGPYGQIFRPDNFVFGQSGAGNNWAKGHYTEGAELIDSVLDVV
                 60 |||||||||||||||||||||||||||||||||||||||||||||||||||||||.||||
4                60 PRAVLMDLEPGTMDSVRSGPYGQIFRPDNFVFGQSGAGNNWAKGHYTEGAELIDSXLDVV

gi|154512       120 RKEAENCDCLQGFQVCHSLGGGTGSGMGTLLISKIREEYPDRMMLTFSVFPSPKVSDTVV
                120 ||||||||||||||||||||.|||||||||||||||||||||||||||||||||||||||
4               120 RKEAENCDCLQGFQVCHSLGRGTGSGMGTLLISKIREEYPDRMMLTFSVFPSPKVSDTVV

gi|154512       180 EPYNATLSVHQLVENADECMV 201
                180 ||||||||||.|||||||||| 201
4               180 EPYNATLSVHXLVENADECMV 201
""",
            )
            hit = record[5]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|225470745|ref|XP_002267380.1|")
            self.assertEqual(hit.target.name, "XP_002267380")
            self.assertEqual(
                hit.target.description,
                "PREDICTED: hypothetical protein [Vitis vinifera] >gi|157327486|emb|CAO15467.1| unnamed protein product [Vitis vinifera]",
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=449)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 1033.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 402.51939420267)
            self.assertAlmostEqual(hsp.annotations["evalue"], 1.24089932237309e-110)
            self.assertEqual(hsp.annotations["identity"], 192)
            self.assertEqual(hsp.annotations["positive"], 195)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[  0, 201],
                              [  0, 201]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 201))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('MREILHIQAGQCGNQIGANFWEVVCAEHGINSTGRYQGDNDLQLERVNVYYNEA...CMV')",
            )
            self.assertEqual(hsp.query.id, "4")
            self.assertEqual(
                hsp.query.description,
                "gi|5049839|gb|AI730987.1|AI730987 BNLGHi8354 Six-day Cotton fiber Gossypium hirsutum cDNA 5' similar to TUBULIN BETA-1 CHAIN gi|486734|pir|S35142 tubulin beta chain - white lupine gi|402636 (X70184) Beta tubulin 1 [Lupinus albus], mRNA sequence",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(201))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "4:50..652")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq({0: 'MREILHVQGGQCGNQIGAKFWEVVCAEHGIDSTGRYQGDTELQLERVNVYYNEA...CMV'}, length=449)",
            )
            self.assertEqual(hsp.target.id, "gi|225470745|ref|XP_002267380.1|")
            self.assertEqual(hsp.target.name, "XP_002267380")
            self.assertEqual(
                hsp.target.description,
                "PREDICTED: hypothetical protein [Vitis vinifera] >gi|157327486|emb|CAO15467.1| unnamed protein product [Vitis vinifera]",
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "MREILH+Q GQCGNQIGA FWEVVCAEHGI+STGRYQGD +LQLERVNVYYNEASCGRFVPRAVLMDLEPGTMDSVRSGPYGQIFRPDNFVFGQSGAGNNWAKGHYTEGAELIDS LDVVRKEAENCDCLQGFQVCHSLG GTGSGMGTLLISKIREEYPDRMMLTFSVFPSPKVSDTVVEPYNATLSVH LVENADECMV",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|225470         0 MREILHVQGGQCGNQIGAKFWEVVCAEHGIDSTGRYQGDTELQLERVNVYYNEASCGRFV
                  0 ||||||.|.|||||||||.|||||||||||.||||||||..|||||||||||||||||||
4                 0 MREILHIQAGQCGNQIGANFWEVVCAEHGINSTGRYQGDNDLQLERVNVYYNEASCGRFV

gi|225470        60 PRAVLMDLEPGTMDSVRSGPYGQIFRPDNFVFGQSGAGNNWAKGHYTEGAELIDSVLDVV
                 60 |||||||||||||||||||||||||||||||||||||||||||||||||||||||.||||
4                60 PRAVLMDLEPGTMDSVRSGPYGQIFRPDNFVFGQSGAGNNWAKGHYTEGAELIDSXLDVV

gi|225470       120 RKEAENCDCLQGFQVCHSLGGGTGSGMGTLLISKIREEYPDRMMLTFSVFPSPKVSDTVV
                120 ||||||||||||||||||||.|||||||||||||||||||||||||||||||||||||||
4               120 RKEAENCDCLQGFQVCHSLGRGTGSGMGTLLISKIREEYPDRMMLTFSVFPSPKVSDTVV

gi|225470       180 EPYNATLSVHQLVENADECMV 201
                180 ||||||||||.|||||||||| 201
4               180 EPYNATLSVHXLVENADECMV 201
""",
            )
            hit = record[6]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|586076|sp|P37392.1|TBB1_LUPAL")
            self.assertEqual(hit.target.name, "P37392")
            self.assertEqual(
                hit.target.description,
                "RecName: Full=Tubulin beta-1 chain; AltName: Full=Beta-1-tubulin >gi|402636|emb|CAA49736.1| Beta tubulin 1 [Lupinus albus]",
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=447)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 1033.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 402.51939420267)
            self.assertAlmostEqual(hsp.annotations["evalue"], 1.24089932237309e-110)
            self.assertEqual(hsp.annotations["identity"], 193)
            self.assertEqual(hsp.annotations["positive"], 195)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[  0, 201],
                              [  0, 201]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 201))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('MREILHIQAGQCGNQIGANFWEVVCAEHGINSTGRYQGDNDLQLERVNVYYNEA...CMV')",
            )
            self.assertEqual(hsp.query.id, "4")
            self.assertEqual(
                hsp.query.description,
                "gi|5049839|gb|AI730987.1|AI730987 BNLGHi8354 Six-day Cotton fiber Gossypium hirsutum cDNA 5' similar to TUBULIN BETA-1 CHAIN gi|486734|pir|S35142 tubulin beta chain - white lupine gi|402636 (X70184) Beta tubulin 1 [Lupinus albus], mRNA sequence",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(201))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "4:50..652")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq({0: 'MREILHIQGGQCGNQIGAKFWEVVCAEHGIDSTGRYGGDNELQLERVNVYYNEA...CMV'}, length=447)",
            )
            self.assertEqual(hsp.target.id, "gi|586076|sp|P37392.1|TBB1_LUPAL")
            self.assertEqual(hsp.target.name, "P37392")
            self.assertEqual(
                hsp.target.description,
                "RecName: Full=Tubulin beta-1 chain; AltName: Full=Beta-1-tubulin >gi|402636|emb|CAA49736.1| Beta tubulin 1 [Lupinus albus]",
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "MREILHIQ GQCGNQIGA FWEVVCAEHGI+STGRY GDN+LQLERVNVYYNEASCGRFVPRAVLMDLEPGTMDSVRSGPYGQIFRPDNFVFGQSGAGNNWAKGHYTEGAELIDS LDVVRKEAENCDCLQGFQVCHSLG GTGSGMGTLLISKIREEYPDRMMLTFSVFPSPKVSDTVVEPYNATLSVH LVENADECMV",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|586076         0 MREILHIQGGQCGNQIGAKFWEVVCAEHGIDSTGRYGGDNELQLERVNVYYNEASCGRFV
                  0 ||||||||.|||||||||.|||||||||||.|||||.|||.|||||||||||||||||||
4                 0 MREILHIQAGQCGNQIGANFWEVVCAEHGINSTGRYQGDNDLQLERVNVYYNEASCGRFV

gi|586076        60 PRAVLMDLEPGTMDSVRSGPYGQIFRPDNFVFGQSGAGNNWAKGHYTEGAELIDSVLDVV
                 60 |||||||||||||||||||||||||||||||||||||||||||||||||||||||.||||
4                60 PRAVLMDLEPGTMDSVRSGPYGQIFRPDNFVFGQSGAGNNWAKGHYTEGAELIDSXLDVV

gi|586076       120 RKEAENCDCLQGFQVCHSLGGGTGSGMGTLLISKIREEYPDRMMLTFSVFPSPKVSDTVV
                120 ||||||||||||||||||||.|||||||||||||||||||||||||||||||||||||||
4               120 RKEAENCDCLQGFQVCHSLGRGTGSGMGTLLISKIREEYPDRMMLTFSVFPSPKVSDTVV

gi|586076       180 EPYNATLSVHQLVENADECMV 201
                180 ||||||||||.|||||||||| 201
4               180 EPYNATLSVHXLVENADECMV 201
""",
            )
            hit = record[7]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|224104341|ref|XP_002313404.1|")
            self.assertEqual(hit.target.name, "XP_002313404")
            self.assertEqual(
                hit.target.description,
                "tubulin, beta chain [Populus trichocarpa] >gi|222849812|gb|EEE87359.1| tubulin, beta chain [Populus trichocarpa]",
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=451)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 1031.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 401.748995050835)
            self.assertAlmostEqual(hsp.annotations["evalue"], 2.1166536544662e-110)
            self.assertEqual(hsp.annotations["identity"], 192)
            self.assertEqual(hsp.annotations["positive"], 195)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[  0, 201],
                              [  0, 201]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 201))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('MREILHIQAGQCGNQIGANFWEVVCAEHGINSTGRYQGDNDLQLERVNVYYNEA...CMV')",
            )
            self.assertEqual(hsp.query.id, "4")
            self.assertEqual(
                hsp.query.description,
                "gi|5049839|gb|AI730987.1|AI730987 BNLGHi8354 Six-day Cotton fiber Gossypium hirsutum cDNA 5' similar to TUBULIN BETA-1 CHAIN gi|486734|pir|S35142 tubulin beta chain - white lupine gi|402636 (X70184) Beta tubulin 1 [Lupinus albus], mRNA sequence",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(201))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "4:50..652")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq({0: 'MREILHIQGGQCGNQIGAKFWEVVCAEHGIDSTGRYQGDSPLQLERINVYYNEA...CMV'}, length=451)",
            )
            self.assertEqual(hsp.target.id, "gi|224104341|ref|XP_002313404.1|")
            self.assertEqual(hsp.target.name, "XP_002313404")
            self.assertEqual(
                hsp.target.description,
                "tubulin, beta chain [Populus trichocarpa] >gi|222849812|gb|EEE87359.1| tubulin, beta chain [Populus trichocarpa]",
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "MREILHIQ GQCGNQIGA FWEVVCAEHGI+STGRYQGD+ LQLER+NVYYNEASCGRFVPRAVLMDLEPGTMDSVRSGPYGQIFRPDNFVFGQSGAGNNWAKGHYTEGAELIDS LDVVRKEAENCDCLQGFQVCHSLG GTGSGMGTLLISKIREEYPDRMMLTFSVFPSPKVSDTVVEPYNATLSVH LVENADECMV",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|224104         0 MREILHIQGGQCGNQIGAKFWEVVCAEHGIDSTGRYQGDSPLQLERINVYYNEASCGRFV
                  0 ||||||||.|||||||||.|||||||||||.||||||||..|||||.|||||||||||||
4                 0 MREILHIQAGQCGNQIGANFWEVVCAEHGINSTGRYQGDNDLQLERVNVYYNEASCGRFV

gi|224104        60 PRAVLMDLEPGTMDSVRSGPYGQIFRPDNFVFGQSGAGNNWAKGHYTEGAELIDSVLDVV
                 60 |||||||||||||||||||||||||||||||||||||||||||||||||||||||.||||
4                60 PRAVLMDLEPGTMDSVRSGPYGQIFRPDNFVFGQSGAGNNWAKGHYTEGAELIDSXLDVV

gi|224104       120 RKEAENCDCLQGFQVCHSLGGGTGSGMGTLLISKIREEYPDRMMLTFSVFPSPKVSDTVV
                120 ||||||||||||||||||||.|||||||||||||||||||||||||||||||||||||||
4               120 RKEAENCDCLQGFQVCHSLGRGTGSGMGTLLISKIREEYPDRMMLTFSVFPSPKVSDTVV

gi|224104       180 EPYNATLSVHQLVENADECMV 201
                180 ||||||||||.|||||||||| 201
4               180 EPYNATLSVHXLVENADECMV 201
""",
            )
            hit = record[8]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|223549679|gb|EEF51167.1|")
            self.assertEqual(hit.target.name, "EEF51167")
            self.assertEqual(
                hit.target.description,
                "tubulin beta chain, putative [Ricinus communis]",
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=446)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 1029.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 400.978595899)
            self.assertAlmostEqual(hsp.annotations["evalue"], 3.61046429165375e-110)
            self.assertEqual(hsp.annotations["identity"], 191)
            self.assertEqual(hsp.annotations["positive"], 194)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[  0, 201],
                              [  0, 201]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 201))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('MREILHIQAGQCGNQIGANFWEVVCAEHGINSTGRYQGDNDLQLERVNVYYNEA...CMV')",
            )
            self.assertEqual(hsp.query.id, "4")
            self.assertEqual(
                hsp.query.description,
                "gi|5049839|gb|AI730987.1|AI730987 BNLGHi8354 Six-day Cotton fiber Gossypium hirsutum cDNA 5' similar to TUBULIN BETA-1 CHAIN gi|486734|pir|S35142 tubulin beta chain - white lupine gi|402636 (X70184) Beta tubulin 1 [Lupinus albus], mRNA sequence",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(201))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "4:50..652")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq({0: 'MREILHVQGGQCGNQIGAKFWEVVCAEHGIDSTGRYHGDTDLQLERVNVYYNEA...CMV'}, length=446)",
            )
            self.assertEqual(hsp.target.id, "gi|223549679|gb|EEF51167.1|")
            self.assertEqual(hsp.target.name, "EEF51167")
            self.assertEqual(
                hsp.target.description,
                "tubulin beta chain, putative [Ricinus communis]",
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "MREILH+Q GQCGNQIGA FWEVVCAEHGI+STGRY GD DLQLERVNVYYNEASCGRFVPRAVLMDLEPGTMDSVRSGPYGQIFRPDNFVFGQSGAGNNWAKGHYTEGAELIDS LDVVRKEAENCDCLQGFQVCHSLG GTGSGMGTLLISK+REEYPDRMMLTFSVFPSPKVSDTVVEPYNATLSVH LVENADECMV",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|223549         0 MREILHVQGGQCGNQIGAKFWEVVCAEHGIDSTGRYHGDTDLQLERVNVYYNEASCGRFV
                  0 ||||||.|.|||||||||.|||||||||||.|||||.||.||||||||||||||||||||
4                 0 MREILHIQAGQCGNQIGANFWEVVCAEHGINSTGRYQGDNDLQLERVNVYYNEASCGRFV

gi|223549        60 PRAVLMDLEPGTMDSVRSGPYGQIFRPDNFVFGQSGAGNNWAKGHYTEGAELIDSVLDVV
                 60 |||||||||||||||||||||||||||||||||||||||||||||||||||||||.||||
4                60 PRAVLMDLEPGTMDSVRSGPYGQIFRPDNFVFGQSGAGNNWAKGHYTEGAELIDSXLDVV

gi|223549       120 RKEAENCDCLQGFQVCHSLGGGTGSGMGTLLISKMREEYPDRMMLTFSVFPSPKVSDTVV
                120 ||||||||||||||||||||.|||||||||||||.|||||||||||||||||||||||||
4               120 RKEAENCDCLQGFQVCHSLGRGTGSGMGTLLISKIREEYPDRMMLTFSVFPSPKVSDTVV

gi|223549       180 EPYNATLSVHQLVENADECMV 201
                180 ||||||||||.|||||||||| 201
4               180 EPYNATLSVHXLVENADECMV 201
""",
            )
            hit = record[9]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|224058553|ref|XP_002299541.1|")
            self.assertEqual(hit.target.name, "XP_002299541")
            self.assertEqual(
                hit.target.description,
                "tubulin, beta chain [Populus trichocarpa] >gi|222846799|gb|EEE84346.1| tubulin, beta chain [Populus trichocarpa]",
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=447)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 1029.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 400.978595899)
            self.assertAlmostEqual(hsp.annotations["evalue"], 3.61046429165375e-110)
            self.assertEqual(hsp.annotations["identity"], 192)
            self.assertEqual(hsp.annotations["positive"], 195)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[  0, 201],
                              [  0, 201]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 201))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('MREILHIQAGQCGNQIGANFWEVVCAEHGINSTGRYQGDNDLQLERVNVYYNEA...CMV')",
            )
            self.assertEqual(hsp.query.id, "4")
            self.assertEqual(
                hsp.query.description,
                "gi|5049839|gb|AI730987.1|AI730987 BNLGHi8354 Six-day Cotton fiber Gossypium hirsutum cDNA 5' similar to TUBULIN BETA-1 CHAIN gi|486734|pir|S35142 tubulin beta chain - white lupine gi|402636 (X70184) Beta tubulin 1 [Lupinus albus], mRNA sequence",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(201))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "4:50..652")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq({0: 'MREILHIQGGQCGNQIGAKFWEVVCAEHGIDSTGRYQGDSALQIERVNVYYNEA...CMV'}, length=447)",
            )
            self.assertEqual(hsp.target.id, "gi|224058553|ref|XP_002299541.1|")
            self.assertEqual(hsp.target.name, "XP_002299541")
            self.assertEqual(
                hsp.target.description,
                "tubulin, beta chain [Populus trichocarpa] >gi|222846799|gb|EEE84346.1| tubulin, beta chain [Populus trichocarpa]",
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "MREILHIQ GQCGNQIGA FWEVVCAEHGI+STGRYQGD+ LQ+ERVNVYYNEASCGRFVPRAVLMDLEPGTMDSVRSGPYGQIFRPDNFVFGQSGAGNNWAKGHYTEGAELIDS LDVVRKEAENCDCLQGFQVCHSLG GTGSGMGTLLISKIREEYPDRMMLTFSVFPSPKVSDTVVEPYNATLSVH LVENADECMV",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|224058         0 MREILHIQGGQCGNQIGAKFWEVVCAEHGIDSTGRYQGDSALQIERVNVYYNEASCGRFV
                  0 ||||||||.|||||||||.|||||||||||.||||||||..||.||||||||||||||||
4                 0 MREILHIQAGQCGNQIGANFWEVVCAEHGINSTGRYQGDNDLQLERVNVYYNEASCGRFV

gi|224058        60 PRAVLMDLEPGTMDSVRSGPYGQIFRPDNFVFGQSGAGNNWAKGHYTEGAELIDSVLDVV
                 60 |||||||||||||||||||||||||||||||||||||||||||||||||||||||.||||
4                60 PRAVLMDLEPGTMDSVRSGPYGQIFRPDNFVFGQSGAGNNWAKGHYTEGAELIDSXLDVV

gi|224058       120 RKEAENCDCLQGFQVCHSLGGGTGSGMGTLLISKIREEYPDRMMLTFSVFPSPKVSDTVV
                120 ||||||||||||||||||||.|||||||||||||||||||||||||||||||||||||||
4               120 RKEAENCDCLQGFQVCHSLGRGTGSGMGTLLISKIREEYPDRMMLTFSVFPSPKVSDTVV

gi|224058       180 EPYNATLSVHQLVENADECMV 201
                180 ||||||||||.|||||||||| 201
4               180 EPYNATLSVHXLVENADECMV 201
""",
            )
            record = next(records)
            self.assertIsInstance(record.query, SeqRecord)
            self.assertEqual(record.query.id, "5")
            self.assertEqual(
                record.query.description,
                "gi|5052071|gb|AF067555.1|AF067555 Phlox stansburyi internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2, complete sequence",
            )
            self.assertEqual(repr(record.query.seq), "Seq(None, length=623)")

            self.assertEqual(len(record.stat), 7)
            self.assertEqual(record.stat["db-num"], 8994603)
            self.assertEqual(record.stat["db-len"], -1216159329)
            self.assertEqual(record.stat["hsp-len"], 0)
            self.assertAlmostEqual(record.stat["eff-space"], 147032237429.0)
            self.assertAlmostEqual(record.stat["kappa"], 0.041)
            self.assertAlmostEqual(record.stat["lambda"], 0.267)
            self.assertAlmostEqual(record.stat["entropy"], 0.14)
            self.assertEqual(len(record), 10)
            hit = record[0]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|110740644|dbj|BAE98425.1|")
            self.assertEqual(hit.target.name, "BAE98425")
            self.assertEqual(
                hit.target.description, "hypothetical protein [Arabidopsis thaliana]"
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=80)")
            self.assertEqual(len(hit), 2)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 231.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 93.5893343169526)
            self.assertAlmostEqual(hsp.annotations["evalue"], 5.57283114448317e-19)
            self.assertEqual(hsp.annotations["identity"], 42)
            self.assertEqual(hsp.annotations["positive"], 45)
            self.assertEqual(hsp.annotations["gaps"], 1)
            hsp = hit[1]
            self.assertAlmostEqual(hsp.score, 53.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 25.0238098036637)
            self.assertAlmostEqual(hsp.annotations["evalue"], 5.57283114448317e-19)
            self.assertEqual(hsp.annotations["identity"], 13)
            self.assertEqual(hsp.annotations["positive"], 13)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[53, 70],
                              [ 0, 17]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 17))
            self.assertEqual(repr(hsp.query.seq), "Seq('RKLVSRALRCAVGLNKS')")
            self.assertEqual(hsp.query.id, "5")
            self.assertEqual(
                hsp.query.description,
                "gi|5052071|gb|AF067555.1|AF067555 Phlox stansburyi internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2, complete sequence",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(17))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "5:455..505")
            self.assertEqual(
                repr(hsp.target.seq), "Seq({53: 'RKLVSRVLPHAVGLNPS'}, length=80)"
            )
            self.assertEqual(hsp.target.id, "gi|110740644|dbj|BAE98425.1|")
            self.assertEqual(hsp.target.name, "BAE98425")
            self.assertEqual(
                hsp.target.description, "hypothetical protein [Arabidopsis thaliana]"
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(hsp.annotations["midline"], "RKLVSR L  AVGLN S")
            self.assertEqual(
                str(hsp),
                """\
gi|110740        53 RKLVSRVLPHAVGLNPS 70
                  0 ||||||.|..|||||.| 17
5                 0 RKLVSRALRCAVGLNKS 17
""",
            )
            hit = record[1]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|226453533|gb|EEH50844.1|")
            self.assertEqual(hit.target.name, "EEH50844")
            self.assertEqual(
                hit.target.description,
                "predicted protein [Micromonas pusilla CCMP1545]",
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=81)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 238.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 96.2857313483741)
            self.assertAlmostEqual(hsp.annotations["evalue"], 1.69151855577931e-18)
            self.assertEqual(hsp.annotations["identity"], 42)
            self.assertEqual(hsp.annotations["positive"], 45)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[ 0, 49],
                              [ 0, 49]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 49))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('MKNVAKCDTWCELQNPVNHRVFERKLRPKPLGRGHVCLGVSHRVAPNPF')",
            )
            self.assertEqual(hsp.query.id, "5")
            self.assertEqual(
                hsp.query.description,
                "gi|5052071|gb|AF067555.1|AF067555 Phlox stansburyi internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2, complete sequence",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(49))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "5:283..429")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq({0: 'MKNVAKCDTWCELQNPVNHRVFERKLRPKPSGRGHVCLGVTNRRPPSSF'}, length=81)",
            )
            self.assertEqual(hsp.target.id, "gi|226453533|gb|EEH50844.1|")
            self.assertEqual(hsp.target.name, "EEH50844")
            self.assertEqual(
                hsp.target.description,
                "predicted protein [Micromonas pusilla CCMP1545]",
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "MKNVAKCDTWCELQNPVNHRVFERKLRPKP GRGHVCLGV++R  P+ F",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|226453         0 MKNVAKCDTWCELQNPVNHRVFERKLRPKPSGRGHVCLGVTNRRPPSSF 49
                  0 ||||||||||||||||||||||||||||||.|||||||||..|..|..| 49
5                 0 MKNVAKCDTWCELQNPVNHRVFERKLRPKPLGRGHVCLGVSHRVAPNPF 49
""",
            )
            hit = record[2]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|168069582|ref|XP_001786502.1|")
            self.assertEqual(hit.target.name, "XP_001786502")
            self.assertEqual(
                hit.target.description,
                "predicted protein [Physcomitrella patens subsp. patens] >gi|162661153|gb|EDQ48685.1| predicted protein [Physcomitrella patens subsp. patens]",
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=88)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 183.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 75.0997546729196)
            self.assertAlmostEqual(hsp.annotations["evalue"], 4.03544314604194e-12)
            self.assertEqual(hsp.annotations["identity"], 37)
            self.assertEqual(hsp.annotations["positive"], 39)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[ 2, 44],
                              [ 0, 42]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 42))
            self.assertEqual(
                repr(hsp.query.seq), "Seq('ASGATCVQKLDGSRDSAIHTKYRISLRSSSMREPRYPLPRVV')"
            )
            self.assertEqual(hsp.query.id, "5")
            self.assertEqual(
                hsp.query.description,
                "gi|5052071|gb|AF067555.1|AF067555 Phlox stansburyi internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2, complete sequence",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(42))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "complement(5:245..370)")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq({2: 'ASGATCVQKLDDSRDSAIHTTYRISLRSSSLQEPRYPLLRVV'}, length=88)",
            )
            self.assertEqual(hsp.target.id, "gi|168069582|ref|XP_001786502.1|")
            self.assertEqual(hsp.target.name, "XP_001786502")
            self.assertEqual(
                hsp.target.description,
                "predicted protein [Physcomitrella patens subsp. patens] >gi|162661153|gb|EDQ48685.1| predicted protein [Physcomitrella patens subsp. patens]",
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"], "ASGATCVQKLD SRDSAIHT YRISLRSSS++EPRYPL RVV"
            )
            self.assertEqual(
                str(hsp),
                """\
gi|168069         2 ASGATCVQKLDDSRDSAIHTTYRISLRSSSLQEPRYPLLRVV 44
                  0 |||||||||||.||||||||.|||||||||..||||||.||| 42
5                 0 ASGATCVQKLDGSRDSAIHTKYRISLRSSSMREPRYPLPRVV 42
""",
            )
            hit = record[3]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|168068558|ref|XP_001786120.1|")
            self.assertEqual(hit.target.name, "XP_001786120")
            self.assertEqual(
                hit.target.description,
                "predicted protein [Physcomitrella patens subsp. patens] >gi|162662102|gb|EDQ49068.1| predicted protein [Physcomitrella patens subsp. patens]",
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=130)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 178.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 73.1737567933329)
            self.assertAlmostEqual(hsp.annotations["evalue"], 1.53346675969648e-11)
            self.assertEqual(hsp.annotations["identity"], 36)
            self.assertEqual(hsp.annotations["positive"], 39)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[ 2, 44],
                              [ 0, 42]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 42))
            self.assertEqual(
                repr(hsp.query.seq), "Seq('ASGATCVQKLDGSRDSAIHTKYRISLRSSSMREPRYPLPRVV')"
            )
            self.assertEqual(hsp.query.id, "5")
            self.assertEqual(
                hsp.query.description,
                "gi|5052071|gb|AF067555.1|AF067555 Phlox stansburyi internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2, complete sequence",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(42))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "complement(5:245..370)")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq({2: 'ASGATCVQKLDDSRNSAIHTTYRISLRSSSLQEPRYPLLRVV'}, length=130)",
            )
            self.assertEqual(hsp.target.id, "gi|168068558|ref|XP_001786120.1|")
            self.assertEqual(hsp.target.name, "XP_001786120")
            self.assertEqual(
                hsp.target.description,
                "predicted protein [Physcomitrella patens subsp. patens] >gi|162662102|gb|EDQ49068.1| predicted protein [Physcomitrella patens subsp. patens]",
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"], "ASGATCVQKLD SR+SAIHT YRISLRSSS++EPRYPL RVV"
            )
            self.assertEqual(
                str(hsp),
                """\
gi|168068         2 ASGATCVQKLDDSRNSAIHTTYRISLRSSSLQEPRYPLLRVV 44
                  0 |||||||||||.||.|||||.|||||||||..||||||.||| 42
5                 0 ASGATCVQKLDGSRDSAIHTKYRISLRSSSMREPRYPLPRVV 42
""",
            )
            hit = record[4]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|168068926|ref|XP_001786259.1|")
            self.assertEqual(hit.target.name, "XP_001786259")
            self.assertEqual(
                hit.target.description,
                "predicted protein [Physcomitrella patens subsp. patens] >gi|168069965|ref|XP_001786641.1| predicted protein [Physcomitrella patens subsp. patens] >gi|162660807|gb|EDQ48545.1| predicted protein [Physcomitrella patens subsp. patens] >gi|162661808|gb|EDQ48929.1| predicted protein [Physcomitrella patens subsp. patens]",
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=148)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 178.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 73.1737567933329)
            self.assertAlmostEqual(hsp.annotations["evalue"], 1.53346675969648e-11)
            self.assertEqual(hsp.annotations["identity"], 36)
            self.assertEqual(hsp.annotations["positive"], 39)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[ 2, 44],
                              [ 0, 42]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 42))
            self.assertEqual(
                repr(hsp.query.seq), "Seq('ASGATCVQKLDGSRDSAIHTKYRISLRSSSMREPRYPLPRVV')"
            )
            self.assertEqual(hsp.query.id, "5")
            self.assertEqual(
                hsp.query.description,
                "gi|5052071|gb|AF067555.1|AF067555 Phlox stansburyi internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2, complete sequence",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(42))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "complement(5:245..370)")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq({2: 'ASGATCVQKLDDSRNSAIHTTYRISLRSSSLQEPRYPLLRVV'}, length=148)",
            )
            self.assertEqual(hsp.target.id, "gi|168068926|ref|XP_001786259.1|")
            self.assertEqual(hsp.target.name, "XP_001786259")
            self.assertEqual(
                hsp.target.description,
                "predicted protein [Physcomitrella patens subsp. patens] >gi|168069965|ref|XP_001786641.1| predicted protein [Physcomitrella patens subsp. patens] >gi|162660807|gb|EDQ48545.1| predicted protein [Physcomitrella patens subsp. patens] >gi|162661808|gb|EDQ48929.1| predicted protein [Physcomitrella patens subsp. patens]",
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"], "ASGATCVQKLD SR+SAIHT YRISLRSSS++EPRYPL RVV"
            )
            self.assertEqual(
                str(hsp),
                """\
gi|168068         2 ASGATCVQKLDDSRNSAIHTTYRISLRSSSLQEPRYPLLRVV 44
                  0 |||||||||||.||.|||||.|||||||||..||||||.||| 42
5                 0 ASGATCVQKLDGSRDSAIHTKYRISLRSSSMREPRYPLPRVV 42
""",
            )
            hit = record[5]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|168070288|ref|XP_001786759.1|")
            self.assertEqual(hit.target.name, "XP_001786759")
            self.assertEqual(
                hit.target.description,
                "predicted protein [Physcomitrella patens subsp. patens] >gi|162660550|gb|EDQ48427.1| predicted protein [Physcomitrella patens subsp. patens]",
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=148)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 178.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 73.1737567933329)
            self.assertAlmostEqual(hsp.annotations["evalue"], 1.53346675969648e-11)
            self.assertEqual(hsp.annotations["identity"], 36)
            self.assertEqual(hsp.annotations["positive"], 39)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[ 2, 44],
                              [ 0, 42]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 42))
            self.assertEqual(
                repr(hsp.query.seq), "Seq('ASGATCVQKLDGSRDSAIHTKYRISLRSSSMREPRYPLPRVV')"
            )
            self.assertEqual(hsp.query.id, "5")
            self.assertEqual(
                hsp.query.description,
                "gi|5052071|gb|AF067555.1|AF067555 Phlox stansburyi internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2, complete sequence",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(42))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "complement(5:245..370)")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq({2: 'ASGATCVQKLDDSRNSAIHTTYRISLRSSSLQEPRYPLLRVV'}, length=148)",
            )
            self.assertEqual(hsp.target.id, "gi|168070288|ref|XP_001786759.1|")
            self.assertEqual(hsp.target.name, "XP_001786759")
            self.assertEqual(
                hsp.target.description,
                "predicted protein [Physcomitrella patens subsp. patens] >gi|162660550|gb|EDQ48427.1| predicted protein [Physcomitrella patens subsp. patens]",
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"], "ASGATCVQKLD SR+SAIHT YRISLRSSS++EPRYPL RVV"
            )
            self.assertEqual(
                str(hsp),
                """\
gi|168070         2 ASGATCVQKLDDSRNSAIHTTYRISLRSSSLQEPRYPLLRVV 44
                  0 |||||||||||.||.|||||.|||||||||..||||||.||| 42
5                 0 ASGATCVQKLDGSRDSAIHTKYRISLRSSSMREPRYPLPRVV 42
""",
            )
            hit = record[6]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|168068591|ref|XP_001786133.1|")
            self.assertEqual(hit.target.name, "XP_001786133")
            self.assertEqual(
                hit.target.description,
                "predicted protein [Physcomitrella patens subsp. patens] >gi|162662081|gb|EDQ49057.1| predicted protein [Physcomitrella patens subsp. patens]",
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=220)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 172.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 70.8625593378288)
            self.assertAlmostEqual(hsp.annotations["evalue"], 7.61051640442713e-11)
            self.assertEqual(hsp.annotations["identity"], 42)
            self.assertEqual(hsp.annotations["positive"], 50)
            self.assertEqual(hsp.annotations["gaps"], 8)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[142, 169, 172, 183, 183, 220],
                              [  0,  27,  27,  38,  43,  80]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 83))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('RPTAHRSARETNFRSQTVESRRKWVGGDAM*DAQADVPSA*WLRAQLAFKNSMV...IRC')",
            )
            self.assertEqual(hsp.query.id, "5")
            self.assertEqual(
                hsp.query.description,
                "gi|5052071|gb|AF067555.1|AF067555 Phlox stansburyi internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2, complete sequence",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(80))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "complement(5:256..495)")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq({142: 'RGLCHHADSDGQFHSTLPIKDIKRIGGCRDDALAGMPSDEPRAQLAFKNSMIHG...IRC'}, length=220)",
            )
            self.assertEqual(hsp.target.id, "gi|168068591|ref|XP_001786133.1|")
            self.assertEqual(hsp.target.name, "XP_001786133")
            self.assertEqual(
                hsp.target.description,
                "predicted protein [Physcomitrella patens subsp. patens] >gi|162662081|gb|EDQ49057.1| predicted protein [Physcomitrella patens subsp. patens]",
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "R   H +  +  F S       K +GG   DA+    +D P     RAQLAFKNSM+HGILQFT  IAFR VLHRC+S+DIRC",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|168068       142 RGLCHHADSDGQFHSTLPIKDIKRIGGCRDDALAGMPSDEP-----RAQLAFKNSMIHGI
                  0 |...|.......|.|.......|..||---||......|.|-----||||||||||.|||
5                 0 RPTAHRSARETNFRSQTVESRRKWVGG---DAM*DAQADVPSA*WLRAQLAFKNSMVHGI

gi|168068       197 LQFTLRIAFRCVLHRCKSQDIRC 220
                 60 ||||..||||.|||||.|.||||  83
5                57 LQFTPSIAFRYVLHRCESRDIRC  80
""",
            )
            hit = record[7]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|74622391|sp|Q8TGM5|ART3_YEAST")
            self.assertEqual(hit.target.name, "Q8TGM5")
            self.assertEqual(
                hit.target.description,
                "Uncharacterized protein ART3 (Antisense to ribosomal RNA transcript protein 3) >gi|18767126|gb|AAL79278.1| unknown [Saccharomyces cerevisiae]",
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=67)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 141.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 58.9213724843908)
            self.assertAlmostEqual(hsp.annotations["evalue"], 2.99274389212967e-07)
            self.assertEqual(hsp.annotations["identity"], 29)
            self.assertEqual(hsp.annotations["positive"], 32)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[ 7, 46],
                              [ 0, 39]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 39))
            self.assertEqual(
                repr(hsp.query.seq), "Seq('GATCVQKLDGSRDSAIHTKYRISLRSSSMREPRYPLPRV')"
            )
            self.assertEqual(hsp.query.id, "5")
            self.assertEqual(
                hsp.query.description,
                "gi|5052071|gb|AF067555.1|AF067555 Phlox stansburyi internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2, complete sequence",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(39))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "complement(5:248..364)")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq({7: 'GAMCVQRFDDSRNSAIHITYRISLRSSSMREPRDPLLKV'}, length=67)",
            )
            self.assertEqual(hsp.target.id, "gi|74622391|sp|Q8TGM5|ART3_YEAST")
            self.assertEqual(hsp.target.name, "Q8TGM5")
            self.assertEqual(
                hsp.target.description,
                "Uncharacterized protein ART3 (Antisense to ribosomal RNA transcript protein 3) >gi|18767126|gb|AAL79278.1| unknown [Saccharomyces cerevisiae]",
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"], "GA CVQ+ D SR+SAIH  YRISLRSSSMREPR PL +V"
            )
            self.assertEqual(
                str(hsp),
                """\
gi|746223         7 GAMCVQRFDDSRNSAIHITYRISLRSSSMREPRDPLLKV 46
                  0 ||.|||..|.||.||||..||||||||||||||.||..| 39
5                 0 GATCVQKLDGSRDSAIHTKYRISLRSSSMREPRYPLPRV 39
""",
            )
            hit = record[8]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|168069944|ref|XP_001786634.1|")
            self.assertEqual(hit.target.name, "XP_001786634")
            self.assertEqual(
                hit.target.description,
                "predicted protein [Physcomitrella patens subsp. patens] >gi|162660825|gb|EDQ48552.1| predicted protein [Physcomitrella patens subsp. patens]",
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=138)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 137.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 57.3805741807214)
            self.assertAlmostEqual(hsp.annotations["evalue"], 8.70755166175354e-07)
            self.assertEqual(hsp.annotations["identity"], 28)
            self.assertEqual(hsp.annotations["positive"], 31)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[ 0, 34],
                              [ 0, 34]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 34))
            self.assertEqual(
                repr(hsp.query.seq), "Seq('KLDGSRDSAIHTKYRISLRSSSMREPRYPLPRVV')"
            )
            self.assertEqual(hsp.query.id, "5")
            self.assertEqual(
                hsp.query.description,
                "gi|5052071|gb|AF067555.1|AF067555 Phlox stansburyi internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2, complete sequence",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(34))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "complement(5:245..346)")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq({0: 'KLDDSRNSAIHTTYRISLRSSSLQEPRYPLLRVV'}, length=138)",
            )
            self.assertEqual(hsp.target.id, "gi|168069944|ref|XP_001786634.1|")
            self.assertEqual(hsp.target.name, "XP_001786634")
            self.assertEqual(
                hsp.target.description,
                "predicted protein [Physcomitrella patens subsp. patens] >gi|162660825|gb|EDQ48552.1| predicted protein [Physcomitrella patens subsp. patens]",
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"], "KLD SR+SAIHT YRISLRSSS++EPRYPL RVV"
            )
            self.assertEqual(
                str(hsp),
                """\
gi|168069         0 KLDDSRNSAIHTTYRISLRSSSLQEPRYPLLRVV 34
                  0 |||.||.|||||.|||||||||..||||||.||| 34
5                 0 KLDGSRDSAIHTKYRISLRSSSMREPRYPLPRVV 34
""",
            )
            hit = record[9]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|50307717|ref|XP_453851.1|")
            self.assertEqual(hit.target.name, "XP_453851")
            self.assertEqual(
                hit.target.description, "unnamed protein product [Kluyveromyces lactis]"
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=54)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 134.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 56.2249754529693)
            self.assertAlmostEqual(hsp.annotations["evalue"], 1.93984013155423e-06)
            self.assertEqual(hsp.annotations["identity"], 28)
            self.assertEqual(hsp.annotations["positive"], 31)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[ 7, 47],
                              [ 0, 40]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 40))
            self.assertEqual(
                repr(hsp.query.seq), "Seq('GATCVQKLDGSRDSAIHTKYRISLRSSSMREPRYPLPRVV')"
            )
            self.assertEqual(hsp.query.id, "5")
            self.assertEqual(
                hsp.query.description,
                "gi|5052071|gb|AF067555.1|AF067555 Phlox stansburyi internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2, complete sequence",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(40))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "complement(5:245..364)")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq({7: 'GAMCVQRFDDSRKSAIHNTYRNSLRSSSMREPRDPLLKVL'}, length=54)",
            )
            self.assertEqual(hsp.target.id, "gi|50307717|ref|XP_453851.1|")
            self.assertEqual(hsp.target.name, "XP_453851")
            self.assertEqual(
                hsp.target.description, "unnamed protein product [Kluyveromyces lactis]"
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"], "GA CVQ+ D SR SAIH  YR SLRSSSMREPR PL +V+"
            )
            self.assertEqual(
                str(hsp),
                """\
gi|503077         7 GAMCVQRFDDSRKSAIHNTYRNSLRSSSMREPRDPLLKVL 47
                  0 ||.|||..|.||.||||..||.|||||||||||.||..|. 40
5                 0 GATCVQKLDGSRDSAIHTKYRISLRSSSMREPRYPLPRVV 40
""",
            )
            record = next(records)
            self.assertIsInstance(record.query, SeqRecord)
            self.assertEqual(record.query.id, "6")
            self.assertEqual(
                record.query.description,
                "gi|3176602|gb|U78617.1|LOU78617 Lathyrus odoratus phytochrome A (PHYA) gene, partial cds",
            )
            self.assertEqual(repr(record.query.seq), "Seq(None, length=309)")

            self.assertEqual(len(record.stat), 7)
            self.assertEqual(record.stat["db-num"], 8994603)
            self.assertEqual(record.stat["db-len"], -1216159329)
            self.assertEqual(record.stat["hsp-len"], 0)
            self.assertAlmostEqual(record.stat["eff-space"], 75367093081.0)
            self.assertAlmostEqual(record.stat["kappa"], 0.041)
            self.assertAlmostEqual(record.stat["lambda"], 0.267)
            self.assertAlmostEqual(record.stat["entropy"], 0.14)
            self.assertEqual(len(record), 10)
            hit = record[0]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|3176603|gb|AAC18749.1|")
            self.assertEqual(hit.target.name, "AAC18749")
            self.assertEqual(
                hit.target.description, "phytochrome A [Lathyrus odoratus]"
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=103)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 543.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 213.771602003167)
            self.assertAlmostEqual(hsp.annotations["evalue"], 3.7262743863676e-54)
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
            self.assertEqual(hsp.shape, (2, 103))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMA...RFV')",
            )
            self.assertEqual(hsp.query.id, "6")
            self.assertEqual(
                hsp.query.description,
                "gi|3176602|gb|U78617.1|LOU78617 Lathyrus odoratus phytochrome A (PHYA) gene, partial cds",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(103))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "6:1..309")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq('QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMA...RFV')",
            )
            self.assertEqual(hsp.target.id, "gi|3176603|gb|AAC18749.1|")
            self.assertEqual(hsp.target.name, "AAC18749")
            self.assertEqual(
                hsp.target.description, "phytochrome A [Lathyrus odoratus]"
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIASLVMAVVVNDSDEDGDSRDAVLPQKKKRLWGLVVCHNTTPRFV",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|317660         0 QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
6                 0 QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA

gi|317660        60 SLVMAVVVNDSDEDGDSRDAVLPQKKKRLWGLVVCHNTTPRFV 103
                 60 ||||||||||||||||||||||||||||||||||||||||||| 103
6                60 SLVMAVVVNDSDEDGDSRDAVLPQKKKRLWGLVVCHNTTPRFV 103
""",
            )
            hit = record[1]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|130188|sp|P15001.1|PHYA_PEA")
            self.assertEqual(hit.target.name, "P15001")
            self.assertEqual(
                hit.target.description,
                "RecName: Full=Phytochrome A >gi|169132|gb|AAA33682.1| phytochrome [Pisum sativum] >gi|295830|emb|CAA32242.1| phytochrome apoprotein [Pisum sativum] >gi|51173514|gb|AAT97643.1| phytochrome A apoprotein [Pisum sativum] >gi|226757|prf||1604466A phytochrome",
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=1124)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 530.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 208.764007516241)
            self.assertAlmostEqual(hsp.annotations["evalue"], 1.1987013044853e-52)
            self.assertEqual(hsp.annotations["identity"], 101)
            self.assertEqual(hsp.annotations["positive"], 102)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[275, 378],
                              [  0, 103]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 103))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMA...RFV')",
            )
            self.assertEqual(hsp.query.id, "6")
            self.assertEqual(
                hsp.query.description,
                "gi|3176602|gb|U78617.1|LOU78617 Lathyrus odoratus phytochrome A (PHYA) gene, partial cds",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(103))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "6:1..309")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq({275: 'QAARFLFMKNKVRMIVDCNAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMA...RFV'}, length=1124)",
            )
            self.assertEqual(hsp.target.id, "gi|130188|sp|P15001.1|PHYA_PEA")
            self.assertEqual(hsp.target.name, "P15001")
            self.assertEqual(
                hsp.target.description,
                "RecName: Full=Phytochrome A >gi|169132|gb|AAA33682.1| phytochrome [Pisum sativum] >gi|295830|emb|CAA32242.1| phytochrome apoprotein [Pisum sativum] >gi|51173514|gb|AAT97643.1| phytochrome A apoprotein [Pisum sativum] >gi|226757|prf||1604466A phytochrome",
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "QAARFLFMKNKVRMIVDC+AKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIASLVMAVVVNDSDEDGDS DAVLPQKKKRLWGLVVCHNTTPRFV",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|130188       275 QAARFLFMKNKVRMIVDCNAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA
                  0 ||||||||||||||||||.|||||||||||||||||||||||||||||||||||||||||
6                 0 QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA

gi|130188       335 SLVMAVVVNDSDEDGDSADAVLPQKKKRLWGLVVCHNTTPRFV 378
                 60 |||||||||||||||||.||||||||||||||||||||||||| 103
6                60 SLVMAVVVNDSDEDGDSRDAVLPQKKKRLWGLVVCHNTTPRFV 103
""",
            )
            hit = record[2]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|2499555|sp|P93673.1|PHYA_LATSA")
            self.assertEqual(hit.target.name, "P93673")
            self.assertEqual(
                hit.target.description,
                "RecName: Full=Phytochrome type A >gi|1848273|gb|AAB47994.1| phytochrome type A [Lathyrus sativus]",
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=1124)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 530.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 208.764007516241)
            self.assertAlmostEqual(hsp.annotations["evalue"], 1.1987013044853e-52)
            self.assertEqual(hsp.annotations["identity"], 101)
            self.assertEqual(hsp.annotations["positive"], 102)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[275, 378],
                              [  0, 103]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 103))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMA...RFV')",
            )
            self.assertEqual(hsp.query.id, "6")
            self.assertEqual(
                hsp.query.description,
                "gi|3176602|gb|U78617.1|LOU78617 Lathyrus odoratus phytochrome A (PHYA) gene, partial cds",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(103))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "6:1..309")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq({275: 'QAARFLFMKNKVRMIVDCNAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMA...RFV'}, length=1124)",
            )
            self.assertEqual(hsp.target.id, "gi|2499555|sp|P93673.1|PHYA_LATSA")
            self.assertEqual(hsp.target.name, "P93673")
            self.assertEqual(
                hsp.target.description,
                "RecName: Full=Phytochrome type A >gi|1848273|gb|AAB47994.1| phytochrome type A [Lathyrus sativus]",
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "QAARFLFMKNKVRMIVDC+AKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIASLVMAVVVNDSDEDGDS DAVLPQKKKRLWGLVVCHNTTPRFV",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|249955       275 QAARFLFMKNKVRMIVDCNAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA
                  0 ||||||||||||||||||.|||||||||||||||||||||||||||||||||||||||||
6                 0 QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA

gi|249955       335 SLVMAVVVNDSDEDGDSADAVLPQKKKRLWGLVVCHNTTPRFV 378
                 60 |||||||||||||||||.||||||||||||||||||||||||| 103
6                60 SLVMAVVVNDSDEDGDSRDAVLPQKKKRLWGLVVCHNTTPRFV 103
""",
            )
            hit = record[3]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|3176595|gb|AAC18745.1|")
            self.assertEqual(hit.target.name, "AAC18745")
            self.assertEqual(
                hit.target.description,
                "phytochrome A [Lennea melanocarpa] >gi|3176597|gb|AAC18746.1| phytochrome A [Hebestigma cubense] >gi|3176609|gb|AAC18752.1| phytochrome A [Sesbania cochichinensis] >gi|3176611|gb|AAC18753.1| phytochrome A [Sesbania emerus]",
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=103)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 528.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 207.993608364407)
            self.assertAlmostEqual(hsp.annotations["evalue"], 2.04467473791515e-52)
            self.assertEqual(hsp.annotations["identity"], 100)
            self.assertEqual(hsp.annotations["positive"], 101)
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
            self.assertEqual(hsp.shape, (2, 103))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMA...RFV')",
            )
            self.assertEqual(hsp.query.id, "6")
            self.assertEqual(
                hsp.query.description,
                "gi|3176602|gb|U78617.1|LOU78617 Lathyrus odoratus phytochrome A (PHYA) gene, partial cds",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(103))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "6:1..309")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq('QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMA...RFV')",
            )
            self.assertEqual(hsp.target.id, "gi|3176595|gb|AAC18745.1|")
            self.assertEqual(hsp.target.name, "AAC18745")
            self.assertEqual(
                hsp.target.description,
                "phytochrome A [Lennea melanocarpa] >gi|3176597|gb|AAC18746.1| phytochrome A [Hebestigma cubense] >gi|3176609|gb|AAC18752.1| phytochrome A [Sesbania cochichinensis] >gi|3176611|gb|AAC18753.1| phytochrome A [Sesbania emerus]",
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIASLVMAVVVNDSDEDGDS DAV PQK+KRLWGLVVCHNTTPRFV",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|317659         0 QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
6                 0 QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA

gi|317659        60 SLVMAVVVNDSDEDGDSSDAVQPQKRKRLWGLVVCHNTTPRFV 103
                 60 |||||||||||||||||.|||.|||.||||||||||||||||| 103
6                60 SLVMAVVVNDSDEDGDSRDAVLPQKKKRLWGLVVCHNTTPRFV 103
""",
            )
            hit = record[4]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|1711106|gb|AAC18675.1|")
            self.assertEqual(hit.target.name, "AAC18675")
            self.assertEqual(hit.target.description, "phytochrome A [Sophora affinis]")
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=210)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 528.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 207.993608364407)
            self.assertAlmostEqual(hsp.annotations["evalue"], 2.04467473791515e-52)
            self.assertEqual(hsp.annotations["identity"], 100)
            self.assertEqual(hsp.annotations["positive"], 101)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[ 40, 143],
                              [  0, 103]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 103))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMA...RFV')",
            )
            self.assertEqual(hsp.query.id, "6")
            self.assertEqual(
                hsp.query.description,
                "gi|3176602|gb|U78617.1|LOU78617 Lathyrus odoratus phytochrome A (PHYA) gene, partial cds",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(103))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "6:1..309")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq({40: 'QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMA...RFV'}, length=210)",
            )
            self.assertEqual(hsp.target.id, "gi|1711106|gb|AAC18675.1|")
            self.assertEqual(hsp.target.name, "AAC18675")
            self.assertEqual(hsp.target.description, "phytochrome A [Sophora affinis]")
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIASLVMAVVVNDSDEDGDS DAV PQK+KRLWGLVVCHNTTPRFV",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|171110        40 QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
6                 0 QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA

gi|171110       100 SLVMAVVVNDSDEDGDSSDAVQPQKRKRLWGLVVCHNTTPRFV 143
                 60 |||||||||||||||||.|||.|||.||||||||||||||||| 103
6                60 SLVMAVVVNDSDEDGDSRDAVLPQKKKRLWGLVVCHNTTPRFV 103
""",
            )
            hit = record[5]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|1711090|gb|AAC18670.1|")
            self.assertEqual(hit.target.name, "AAC18670")
            self.assertEqual(
                hit.target.description, "phytochrome A [Myrospermum sousanum]"
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=210)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 525.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 206.838009636654)
            self.assertAlmostEqual(hsp.annotations["evalue"], 4.55506009801166e-52)
            self.assertEqual(hsp.annotations["identity"], 99)
            self.assertEqual(hsp.annotations["positive"], 101)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[ 40, 143],
                              [  0, 103]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 103))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMA...RFV')",
            )
            self.assertEqual(hsp.query.id, "6")
            self.assertEqual(
                hsp.query.description,
                "gi|3176602|gb|U78617.1|LOU78617 Lathyrus odoratus phytochrome A (PHYA) gene, partial cds",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(103))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "6:1..309")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq({40: 'QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMA...RFV'}, length=210)",
            )
            self.assertEqual(hsp.target.id, "gi|1711090|gb|AAC18670.1|")
            self.assertEqual(hsp.target.name, "AAC18670")
            self.assertEqual(
                hsp.target.description, "phytochrome A [Myrospermum sousanum]"
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIASLV+AVVVNDSDEDGDS DAV PQK+KRLWGLVVCHNTTPRFV",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|171109        40 QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
6                 0 QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA

gi|171109       100 SLVLAVVVNDSDEDGDSSDAVQPQKRKRLWGLVVCHNTTPRFV 143
                 60 |||.|||||||||||||.|||.|||.||||||||||||||||| 103
6                60 SLVMAVVVNDSDEDGDSRDAVLPQKKKRLWGLVVCHNTTPRFV 103
""",
            )
            hit = record[6]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|3176605|gb|AAC18750.1|")
            self.assertEqual(hit.target.name, "AAC18750")
            self.assertEqual(
                hit.target.description, "phytochrome A [Hybosema robustum]"
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=103)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 524.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 206.452810060737)
            self.assertAlmostEqual(hsp.annotations["evalue"], 5.94909272347008e-52)
            self.assertEqual(hsp.annotations["identity"], 99)
            self.assertEqual(hsp.annotations["positive"], 100)
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
            self.assertEqual(hsp.shape, (2, 103))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMA...RFV')",
            )
            self.assertEqual(hsp.query.id, "6")
            self.assertEqual(
                hsp.query.description,
                "gi|3176602|gb|U78617.1|LOU78617 Lathyrus odoratus phytochrome A (PHYA) gene, partial cds",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(103))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "6:1..309")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq('QATRFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMA...RFV')",
            )
            self.assertEqual(hsp.target.id, "gi|3176605|gb|AAC18750.1|")
            self.assertEqual(hsp.target.name, "AAC18750")
            self.assertEqual(
                hsp.target.description, "phytochrome A [Hybosema robustum]"
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "QA RFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIASLVMAVVVNDSDEDGDS DAV PQK+KRLWGLVVCHNTTPRFV",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|317660         0 QATRFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA
                  0 ||.|||||||||||||||||||||||||||||||||||||||||||||||||||||||||
6                 0 QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA

gi|317660        60 SLVMAVVVNDSDEDGDSSDAVQPQKRKRLWGLVVCHNTTPRFV 103
                 60 |||||||||||||||||.|||.|||.||||||||||||||||| 103
6                60 SLVMAVVVNDSDEDGDSRDAVLPQKKKRLWGLVVCHNTTPRFV 103
""",
            )
            hit = record[7]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|3176454|gb|AAC18668.1|")
            self.assertEqual(hit.target.name, "AAC18668")
            self.assertEqual(
                hit.target.description, "phytochrome A [Cyclolobium nutans]"
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=207)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 523.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 206.06761048482)
            self.assertAlmostEqual(hsp.annotations["evalue"], 7.76975571582328e-52)
            self.assertEqual(hsp.annotations["identity"], 99)
            self.assertEqual(hsp.annotations["positive"], 101)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[ 37, 140],
                              [  0, 103]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 103))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMA...RFV')",
            )
            self.assertEqual(hsp.query.id, "6")
            self.assertEqual(
                hsp.query.description,
                "gi|3176602|gb|U78617.1|LOU78617 Lathyrus odoratus phytochrome A (PHYA) gene, partial cds",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(103))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "6:1..309")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq({37: 'QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMA...RFV'}, length=207)",
            )
            self.assertEqual(hsp.target.id, "gi|3176454|gb|AAC18668.1|")
            self.assertEqual(hsp.target.name, "AAC18668")
            self.assertEqual(
                hsp.target.description, "phytochrome A [Cyclolobium nutans]"
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIASLVMAVVVNDSDEDG+S DAV PQK+KRLWGLVVCHNTTPRFV",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|317645        37 QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
6                 0 QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA

gi|317645        97 SLVMAVVVNDSDEDGNSSDAVQPQKRKRLWGLVVCHNTTPRFV 140
                 60 |||||||||||||||.|.|||.|||.||||||||||||||||| 103
6                60 SLVMAVVVNDSDEDGDSRDAVLPQKKKRLWGLVVCHNTTPRFV 103
""",
            )
            hit = record[8]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|3176523|gb|AAC18709.1|")
            self.assertEqual(hit.target.name, "AAC18709")
            self.assertEqual(
                hit.target.description, "phytochrome A [Millettia richardiana]"
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=139)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 521.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 205.297211332985)
            self.assertAlmostEqual(hsp.annotations["evalue"], 1.3253195915005e-51)
            self.assertEqual(hsp.annotations["identity"], 98)
            self.assertEqual(hsp.annotations["positive"], 101)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[ 36, 139],
                              [  0, 103]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 103))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMA...RFV')",
            )
            self.assertEqual(hsp.query.id, "6")
            self.assertEqual(
                hsp.query.description,
                "gi|3176602|gb|U78617.1|LOU78617 Lathyrus odoratus phytochrome A (PHYA) gene, partial cds",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(103))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "6:1..309")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq({36: 'QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMA...RFV'}, length=139)",
            )
            self.assertEqual(hsp.target.id, "gi|3176523|gb|AAC18709.1|")
            self.assertEqual(hsp.target.name, "AAC18709")
            self.assertEqual(
                hsp.target.description, "phytochrome A [Millettia richardiana]"
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIASLVMAVVVND++EDGDS DAV PQK+KRLWGLVVCHNTTPRFV",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|317652        36 QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
6                 0 QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA

gi|317652        96 SLVMAVVVNDNEEDGDSSDAVQPQKRKRLWGLVVCHNTTPRFV 139
                 60 ||||||||||..|||||.|||.|||.||||||||||||||||| 103
6                60 SLVMAVVVNDSDEDGDSRDAVLPQKKKRLWGLVVCHNTTPRFV 103
""",
            )
            hit = record[9]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|3176494|gb|AAC18693.1|")
            self.assertEqual(hit.target.name, "AAC18693")
            self.assertEqual(
                hit.target.description, "phytochrome A [Callerya atropurpurea]"
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=177)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 520.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 204.912011757068)
            self.assertAlmostEqual(hsp.annotations["evalue"], 1.73092099081406e-51)
            self.assertEqual(hsp.annotations["identity"], 98)
            self.assertEqual(hsp.annotations["positive"], 101)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[  7, 110],
                              [  0, 103]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 103))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMA...RFV')",
            )
            self.assertEqual(hsp.query.id, "6")
            self.assertEqual(
                hsp.query.description,
                "gi|3176602|gb|U78617.1|LOU78617 Lathyrus odoratus phytochrome A (PHYA) gene, partial cds",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(103))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "6:1..309")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq({7: 'QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMA...RFV'}, length=177)",
            )
            self.assertEqual(hsp.target.id, "gi|3176494|gb|AAC18693.1|")
            self.assertEqual(hsp.target.name, "AAC18693")
            self.assertEqual(
                hsp.target.description, "phytochrome A [Callerya atropurpurea]"
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIASLVMAVVVNDS+EDGDS +AV PQK+KRLWGLVVCHNTTPRFV",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|317649         7 QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
6                 0 QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA

gi|317649        67 SLVMAVVVNDSEEDGDSSEAVQPQKRKRLWGLVVCHNTTPRFV 110
                 60 |||||||||||.|||||..||.|||.||||||||||||||||| 103
6                60 SLVMAVVVNDSDEDGDSRDAVLPQKKKRLWGLVVCHNTTPRFV 103
""",
            )
            record = next(records)
            self.assertIsInstance(record.query, SeqRecord)
            self.assertEqual(record.query.id, "7")
            self.assertEqual(
                record.query.description,
                "gi|5817701|gb|AF142731.1|AF142731 Wisteria frutescens maturase-like protein (matK) gene, complete cds; chloroplast gene for chloroplast product",
            )
            self.assertEqual(repr(record.query.seq), "Seq(None, length=2551)")

            self.assertEqual(len(record.stat), 7)
            self.assertEqual(record.stat["db-num"], 8994603)
            self.assertEqual(record.stat["db-len"], -1216159329)
            self.assertEqual(record.stat["hsp-len"], 0)
            self.assertAlmostEqual(record.stat["eff-space"], 1251086325060.0)
            self.assertAlmostEqual(record.stat["kappa"], 0.041)
            self.assertAlmostEqual(record.stat["lambda"], 0.267)
            self.assertAlmostEqual(record.stat["entropy"], 0.14)
            self.assertEqual(len(record), 10)
            hit = record[0]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|27805603|sp|Q9TKP6.1|MATK_WISFR")
            self.assertEqual(hit.target.name, "Q9TKP6")
            self.assertEqual(
                hit.target.description,
                "RecName: Full=Maturase K; AltName: Full=Intron maturase >gi|5817759|gb|AAD52902.1|AF142731_1 maturase-like protein [Wisteria frutescens]",
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=506)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 2451.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 948.732392853477)
            self.assertAlmostEqual(hsp.annotations["evalue"], 0.0)
            self.assertEqual(hsp.annotations["identity"], 506)
            self.assertEqual(hsp.annotations["positive"], 506)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[  0, 506],
                              [  0, 506]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 506))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDNKSSL...NHE')",
            )
            self.assertEqual(hsp.query.id, "7")
            self.assertEqual(
                hsp.query.description,
                "gi|5817701|gb|AF142731.1|AF142731 Wisteria frutescens maturase-like protein (matK) gene, complete cds; chloroplast gene for chloroplast product",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(506))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "7:727..2244")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq('MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDNKSSL...NHE')",
            )
            self.assertEqual(hsp.target.id, "gi|27805603|sp|Q9TKP6.1|MATK_WISFR")
            self.assertEqual(hsp.target.name, "Q9TKP6")
            self.assertEqual(
                hsp.target.description,
                "RecName: Full=Maturase K; AltName: Full=Intron maturase >gi|5817759|gb|AAD52902.1|AF142731_1 maturase-like protein [Wisteria frutescens]",
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDNKSSLLIVKRLITRMYQQNHLIISANDSNKNPFLGYNKNFYSQIISDGFAVVVEIPFFLQLSSSLEEAEIVKSYHNLRSIHSIFPFLEDKFTYLNYVSDIRIPYPIHLEILVQILRYWVKDASFFHLLRFFLYHFSNRNSLITPKKSISTFSKSNPRLFLFLYNFYVCEYESIFRFLRNQSSHLRLKSFSVFFERIFFYAKREHLVKVFPKDFSSTLTFFKDPFIHYVRYQGKSILASKNAPLLMNKWKHYFIHLWQCFFDVWSQPGTIHINQLSEHSFHFLGYFSNVRLNRSVVRSQMLQNTFLIEIVIKKLDIIVPIIPLIRSLAKAKFCNVLGHPLSKSVWADSSDFDIIDRFLRICRNLSHYYNGSSKKKNLYRIKYILRLSCIKTLACKHKSTVRAFLKKSGSEELLEEFFTEEEEILSLIFPRTSSTLQRLHRNRIWYLDILFSNDLVNHE",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|278056         0 MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDNKSSLLIVKRL
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
7                 0 MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDNKSSLLIVKRL

gi|278056        60 ITRMYQQNHLIISANDSNKNPFLGYNKNFYSQIISDGFAVVVEIPFFLQLSSSLEEAEIV
                 60 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
7                60 ITRMYQQNHLIISANDSNKNPFLGYNKNFYSQIISDGFAVVVEIPFFLQLSSSLEEAEIV

gi|278056       120 KSYHNLRSIHSIFPFLEDKFTYLNYVSDIRIPYPIHLEILVQILRYWVKDASFFHLLRFF
                120 ||||||||||||||||||||||||||||||||||||||||||||||||||||........
7               120 KSYHNLRSIHSIFPFLEDKFTYLNYVSDIRIPYPIHLEILVQILRYWVKDASXXXXXXXX

gi|278056       180 LYHFSNRNSLITPKKSISTFSKSNPRLFLFLYNFYVCEYESIFRFLRNQSSHLRLKSFSV
                180 ....||||||||||||||||||||||||||||||||||||||||||||||||||||||||
7               180 XXXXSNRNSLITPKKSISTFSKSNPRLFLFLYNFYVCEYESIFRFLRNQSSHLRLKSFSV

gi|278056       240 FFERIFFYAKREHLVKVFPKDFSSTLTFFKDPFIHYVRYQGKSILASKNAPLLMNKWKHY
                240 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
7               240 FFERIFFYAKREHLVKVFPKDFSSTLTFFKDPFIHYVRYQGKSILASKNAPLLMNKWKHY

gi|278056       300 FIHLWQCFFDVWSQPGTIHINQLSEHSFHFLGYFSNVRLNRSVVRSQMLQNTFLIEIVIK
                300 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
7               300 FIHLWQCFFDVWSQPGTIHINQLSEHSFHFLGYFSNVRLNRSVVRSQMLQNTFLIEIVIK

gi|278056       360 KLDIIVPIIPLIRSLAKAKFCNVLGHPLSKSVWADSSDFDIIDRFLRICRNLSHYYNGSS
                360 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
7               360 KLDIIVPIIPLIRSLAKAKFCNVLGHPLSKSVWADSSDFDIIDRFLRICRNLSHYYNGSS

gi|278056       420 KKKNLYRIKYILRLSCIKTLACKHKSTVRAFLKKSGSEELLEEFFTEEEEILSLIFPRTS
                420 ||||||||||||||||||||||||||||||||||||....................||||
7               420 KKKNLYRIKYILRLSCIKTLACKHKSTVRAFLKKSGXXXXXXXXXXXXXXXXXXXXPRTS

gi|278056       480 STLQRLHRNRIWYLDILFSNDLVNHE 506
                480 |||||||||||||||||||||||||| 506
7               480 STLQRLHRNRIWYLDILFSNDLVNHE 506
""",
            )
            hit = record[1]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|171909144|gb|ACB58148.1|")
            self.assertEqual(hit.target.name, "ACB58148")
            self.assertEqual(hit.target.description, "maturase K [Wisteria frutescens]")
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=506)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 2445.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 946.421195397973)
            self.assertAlmostEqual(hsp.annotations["evalue"], 0.0)
            self.assertEqual(hsp.annotations["identity"], 505)
            self.assertEqual(hsp.annotations["positive"], 505)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[  0, 506],
                              [  0, 506]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 506))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDNKSSL...NHE')",
            )
            self.assertEqual(hsp.query.id, "7")
            self.assertEqual(
                hsp.query.description,
                "gi|5817701|gb|AF142731.1|AF142731 Wisteria frutescens maturase-like protein (matK) gene, complete cds; chloroplast gene for chloroplast product",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(506))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "7:727..2244")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq('MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDNKYSL...NHE')",
            )
            self.assertEqual(hsp.target.id, "gi|171909144|gb|ACB58148.1|")
            self.assertEqual(hsp.target.name, "ACB58148")
            self.assertEqual(hsp.target.description, "maturase K [Wisteria frutescens]")
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDNK SLLIVKRLITRMYQQNHLIISANDSNKNPFLGYNKNFYSQIISDGFAVVVEIPFFLQLSSSLEEAEIVKSYHNLRSIHSIFPFLEDKFTYLNYVSDIRIPYPIHLEILVQILRYWVKDASFFHLLRFFLYHFSNRNSLITPKKSISTFSKSNPRLFLFLYNFYVCEYESIFRFLRNQSSHLRLKSFSVFFERIFFYAKREHLVKVFPKDFSSTLTFFKDPFIHYVRYQGKSILASKNAPLLMNKWKHYFIHLWQCFFDVWSQPGTIHINQLSEHSFHFLGYFSNVRLNRSVVRSQMLQNTFLIEIVIKKLDIIVPIIPLIRSLAKAKFCNVLGHPLSKSVWADSSDFDIIDRFLRICRNLSHYYNGSSKKKNLYRIKYILRLSCIKTLACKHKSTVRAFLKKSGSEELLEEFFTEEEEILSLIFPRTSSTLQRLHRNRIWYLDILFSNDLVNHE",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|171909         0 MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDNKYSLLIVKRL
                  0 |||||||||||||||||||||||||||||||||||||||||||||||||||.||||||||
7                 0 MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDNKSSLLIVKRL

gi|171909        60 ITRMYQQNHLIISANDSNKNPFLGYNKNFYSQIISDGFAVVVEIPFFLQLSSSLEEAEIV
                 60 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
7                60 ITRMYQQNHLIISANDSNKNPFLGYNKNFYSQIISDGFAVVVEIPFFLQLSSSLEEAEIV

gi|171909       120 KSYHNLRSIHSIFPFLEDKFTYLNYVSDIRIPYPIHLEILVQILRYWVKDASFFHLLRFF
                120 ||||||||||||||||||||||||||||||||||||||||||||||||||||........
7               120 KSYHNLRSIHSIFPFLEDKFTYLNYVSDIRIPYPIHLEILVQILRYWVKDASXXXXXXXX

gi|171909       180 LYHFSNRNSLITPKKSISTFSKSNPRLFLFLYNFYVCEYESIFRFLRNQSSHLRLKSFSV
                180 ....||||||||||||||||||||||||||||||||||||||||||||||||||||||||
7               180 XXXXSNRNSLITPKKSISTFSKSNPRLFLFLYNFYVCEYESIFRFLRNQSSHLRLKSFSV

gi|171909       240 FFERIFFYAKREHLVKVFPKDFSSTLTFFKDPFIHYVRYQGKSILASKNAPLLMNKWKHY
                240 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
7               240 FFERIFFYAKREHLVKVFPKDFSSTLTFFKDPFIHYVRYQGKSILASKNAPLLMNKWKHY

gi|171909       300 FIHLWQCFFDVWSQPGTIHINQLSEHSFHFLGYFSNVRLNRSVVRSQMLQNTFLIEIVIK
                300 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
7               300 FIHLWQCFFDVWSQPGTIHINQLSEHSFHFLGYFSNVRLNRSVVRSQMLQNTFLIEIVIK

gi|171909       360 KLDIIVPIIPLIRSLAKAKFCNVLGHPLSKSVWADSSDFDIIDRFLRICRNLSHYYNGSS
                360 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
7               360 KLDIIVPIIPLIRSLAKAKFCNVLGHPLSKSVWADSSDFDIIDRFLRICRNLSHYYNGSS

gi|171909       420 KKKNLYRIKYILRLSCIKTLACKHKSTVRAFLKKSGSEELLEEFFTEEEEILSLIFPRTS
                420 ||||||||||||||||||||||||||||||||||||....................||||
7               420 KKKNLYRIKYILRLSCIKTLACKHKSTVRAFLKKSGXXXXXXXXXXXXXXXXXXXXPRTS

gi|171909       480 STLQRLHRNRIWYLDILFSNDLVNHE 506
                480 |||||||||||||||||||||||||| 506
7               480 STLQRLHRNRIWYLDILFSNDLVNHE 506
""",
            )
            hit = record[2]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|171909146|gb|ACB58149.1|")
            self.assertEqual(hit.target.name, "ACB58149")
            self.assertEqual(
                hit.target.description,
                "maturase K [Wisteria frutescens] >gi|171909148|gb|ACB58150.1| maturase K [Wisteria frutescens var. macrostachya]",
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=506)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 2443.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 945.650796246138)
            self.assertAlmostEqual(hsp.annotations["evalue"], 0.0)
            self.assertEqual(hsp.annotations["identity"], 505)
            self.assertEqual(hsp.annotations["positive"], 505)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[  0, 506],
                              [  0, 506]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 506))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDNKSSL...NHE')",
            )
            self.assertEqual(hsp.query.id, "7")
            self.assertEqual(
                hsp.query.description,
                "gi|5817701|gb|AF142731.1|AF142731 Wisteria frutescens maturase-like protein (matK) gene, complete cds; chloroplast gene for chloroplast product",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(506))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "7:727..2244")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq('MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDNKSSL...NHE')",
            )
            self.assertEqual(hsp.target.id, "gi|171909146|gb|ACB58149.1|")
            self.assertEqual(hsp.target.name, "ACB58149")
            self.assertEqual(
                hsp.target.description,
                "maturase K [Wisteria frutescens] >gi|171909148|gb|ACB58150.1| maturase K [Wisteria frutescens var. macrostachya]",
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDNKSSLLIVKRLITRMYQQNHLIISANDSNKNPFLGYNKNFYSQIISDGFAVVVEIPFFLQLSSSLEEAEIVKSYHNLRSIHSIFPFLEDKFTYLNYVSDIRIPYPIHLEILVQILRYWVKDASFFHLLRFFLYHFSNRNSLITP KSISTFSKSNPRLFLFLYNFYVCEYESIFRFLRNQSSHLRLKSFSVFFERIFFYAKREHLVKVFPKDFSSTLTFFKDPFIHYVRYQGKSILASKNAPLLMNKWKHYFIHLWQCFFDVWSQPGTIHINQLSEHSFHFLGYFSNVRLNRSVVRSQMLQNTFLIEIVIKKLDIIVPIIPLIRSLAKAKFCNVLGHPLSKSVWADSSDFDIIDRFLRICRNLSHYYNGSSKKKNLYRIKYILRLSCIKTLACKHKSTVRAFLKKSGSEELLEEFFTEEEEILSLIFPRTSSTLQRLHRNRIWYLDILFSNDLVNHE",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|171909         0 MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDNKSSLLIVKRL
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
7                 0 MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDNKSSLLIVKRL

gi|171909        60 ITRMYQQNHLIISANDSNKNPFLGYNKNFYSQIISDGFAVVVEIPFFLQLSSSLEEAEIV
                 60 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
7                60 ITRMYQQNHLIISANDSNKNPFLGYNKNFYSQIISDGFAVVVEIPFFLQLSSSLEEAEIV

gi|171909       120 KSYHNLRSIHSIFPFLEDKFTYLNYVSDIRIPYPIHLEILVQILRYWVKDASFFHLLRFF
                120 ||||||||||||||||||||||||||||||||||||||||||||||||||||........
7               120 KSYHNLRSIHSIFPFLEDKFTYLNYVSDIRIPYPIHLEILVQILRYWVKDASXXXXXXXX

gi|171909       180 LYHFSNRNSLITPIKSISTFSKSNPRLFLFLYNFYVCEYESIFRFLRNQSSHLRLKSFSV
                180 ....|||||||||.||||||||||||||||||||||||||||||||||||||||||||||
7               180 XXXXSNRNSLITPKKSISTFSKSNPRLFLFLYNFYVCEYESIFRFLRNQSSHLRLKSFSV

gi|171909       240 FFERIFFYAKREHLVKVFPKDFSSTLTFFKDPFIHYVRYQGKSILASKNAPLLMNKWKHY
                240 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
7               240 FFERIFFYAKREHLVKVFPKDFSSTLTFFKDPFIHYVRYQGKSILASKNAPLLMNKWKHY

gi|171909       300 FIHLWQCFFDVWSQPGTIHINQLSEHSFHFLGYFSNVRLNRSVVRSQMLQNTFLIEIVIK
                300 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
7               300 FIHLWQCFFDVWSQPGTIHINQLSEHSFHFLGYFSNVRLNRSVVRSQMLQNTFLIEIVIK

gi|171909       360 KLDIIVPIIPLIRSLAKAKFCNVLGHPLSKSVWADSSDFDIIDRFLRICRNLSHYYNGSS
                360 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
7               360 KLDIIVPIIPLIRSLAKAKFCNVLGHPLSKSVWADSSDFDIIDRFLRICRNLSHYYNGSS

gi|171909       420 KKKNLYRIKYILRLSCIKTLACKHKSTVRAFLKKSGSEELLEEFFTEEEEILSLIFPRTS
                420 ||||||||||||||||||||||||||||||||||||....................||||
7               420 KKKNLYRIKYILRLSCIKTLACKHKSTVRAFLKKSGXXXXXXXXXXXXXXXXXXXXPRTS

gi|171909       480 STLQRLHRNRIWYLDILFSNDLVNHE 506
                480 |||||||||||||||||||||||||| 506
7               480 STLQRLHRNRIWYLDILFSNDLVNHE 506
""",
            )
            hit = record[3]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|171909132|gb|ACB58142.1|")
            self.assertEqual(hit.target.name, "ACB58142")
            self.assertEqual(hit.target.description, "maturase K [Callerya megasperma]")
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=506)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 2439.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 944.109997942469)
            self.assertAlmostEqual(hsp.annotations["evalue"], 0.0)
            self.assertEqual(hsp.annotations["identity"], 501)
            self.assertEqual(hsp.annotations["positive"], 504)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[  0, 506],
                              [  0, 506]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 506))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDNKSSL...NHE')",
            )
            self.assertEqual(hsp.query.id, "7")
            self.assertEqual(
                hsp.query.description,
                "gi|5817701|gb|AF142731.1|AF142731 Wisteria frutescens maturase-like protein (matK) gene, complete cds; chloroplast gene for chloroplast product",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(506))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "7:727..2244")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq('MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDNKSSL...NHE')",
            )
            self.assertEqual(hsp.target.id, "gi|171909132|gb|ACB58142.1|")
            self.assertEqual(hsp.target.name, "ACB58142")
            self.assertEqual(hsp.target.description, "maturase K [Callerya megasperma]")
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDNKSSLLIVKRLITRMYQQNHLIISANDSNKNPFLGYNKNFYSQIISDGFAVVVEIPFFLQLSSSLEEAEIVKSYHNLRSIHSIFPFLEDKFTYLNYVSDIRIPYPIHLEILVQILRYWVKDASFFHLLRFFLY++ NRNSLITPKKSISTFSKSNPRLFLFLYNFYVCEYESIFRFLRNQSSHLRLKSFSVFFERIFFYAKREHLVKVFPKDFSSTLTFFKDPFIHYVRYQGKSILASKNAPLLMNKWKHYFIHLWQCFFDVWSQPGTIHINQLSEHSFHFLGYFSNVRLNRSVVRSQMLQNTFLIEIVIKKLDIIVPIIPLIRSLAKAKFCNVLGHP+SKSVWADSSDFDIIDRFLRICRNLSHYYNGSSKKKNLYRIKYILRLSCIKTLACKHKSTVRAFLKKSGSEELLEEFFTEEEEILSLIFPR SSTLQRLHRNRIWYLDILFSNDLVNHE",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|171909         0 MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDNKSSLLIVKRL
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
7                 0 MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDNKSSLLIVKRL

gi|171909        60 ITRMYQQNHLIISANDSNKNPFLGYNKNFYSQIISDGFAVVVEIPFFLQLSSSLEEAEIV
                 60 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
7                60 ITRMYQQNHLIISANDSNKNPFLGYNKNFYSQIISDGFAVVVEIPFFLQLSSSLEEAEIV

gi|171909       120 KSYHNLRSIHSIFPFLEDKFTYLNYVSDIRIPYPIHLEILVQILRYWVKDASFFHLLRFF
                120 ||||||||||||||||||||||||||||||||||||||||||||||||||||........
7               120 KSYHNLRSIHSIFPFLEDKFTYLNYVSDIRIPYPIHLEILVQILRYWVKDASXXXXXXXX

gi|171909       180 LYNYCNRNSLITPKKSISTFSKSNPRLFLFLYNFYVCEYESIFRFLRNQSSHLRLKSFSV
                180 .....|||||||||||||||||||||||||||||||||||||||||||||||||||||||
7               180 XXXXSNRNSLITPKKSISTFSKSNPRLFLFLYNFYVCEYESIFRFLRNQSSHLRLKSFSV

gi|171909       240 FFERIFFYAKREHLVKVFPKDFSSTLTFFKDPFIHYVRYQGKSILASKNAPLLMNKWKHY
                240 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
7               240 FFERIFFYAKREHLVKVFPKDFSSTLTFFKDPFIHYVRYQGKSILASKNAPLLMNKWKHY

gi|171909       300 FIHLWQCFFDVWSQPGTIHINQLSEHSFHFLGYFSNVRLNRSVVRSQMLQNTFLIEIVIK
                300 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
7               300 FIHLWQCFFDVWSQPGTIHINQLSEHSFHFLGYFSNVRLNRSVVRSQMLQNTFLIEIVIK

gi|171909       360 KLDIIVPIIPLIRSLAKAKFCNVLGHPISKSVWADSSDFDIIDRFLRICRNLSHYYNGSS
                360 |||||||||||||||||||||||||||.||||||||||||||||||||||||||||||||
7               360 KLDIIVPIIPLIRSLAKAKFCNVLGHPLSKSVWADSSDFDIIDRFLRICRNLSHYYNGSS

gi|171909       420 KKKNLYRIKYILRLSCIKTLACKHKSTVRAFLKKSGSEELLEEFFTEEEEILSLIFPRAS
                420 ||||||||||||||||||||||||||||||||||||....................||.|
7               420 KKKNLYRIKYILRLSCIKTLACKHKSTVRAFLKKSGXXXXXXXXXXXXXXXXXXXXPRTS

gi|171909       480 STLQRLHRNRIWYLDILFSNDLVNHE 506
                480 |||||||||||||||||||||||||| 506
7               480 STLQRLHRNRIWYLDILFSNDLVNHE 506
""",
            )
            hit = record[4]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|5817760|gb|AAD52903.1|AF142732_1")
            self.assertEqual(hit.target.name, "AAD52903")
            self.assertEqual(
                hit.target.description,
                "maturase-like protein [Wisteria sinensis] >gi|171909136|gb|ACB58144.1| maturase K [Wisteria brachybotrys] >gi|171909138|gb|ACB58145.1| maturase K [Wisteria floribunda] >gi|171909140|gb|ACB58146.1| maturase K [Wisteria floribunda] >gi|171909142|gb|ACB58147.1| maturase K [Wisteria floribunda] >gi|171909150|gb|ACB58151.1| maturase K [Wisteria sinensis] >gi|171909152|gb|ACB58152.1| maturase K [Wisteria sinensis] >gi|171909154|gb|ACB58153.1| maturase K [Wisteria sinensis] >gi|171909156|gb|ACB58154.1| maturase K [Wisteria villosa] >gi|171909158|gb|ACB58155.1| maturase K [Wisteria villosa] >gi|171909160|gb|ACB58156.1| maturase K [Wisteria villosa]",
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=506)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 2418.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 936.020806848204)
            self.assertAlmostEqual(hsp.annotations["evalue"], 0.0)
            self.assertEqual(hsp.annotations["identity"], 498)
            self.assertEqual(hsp.annotations["positive"], 500)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[  0, 506],
                              [  0, 506]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 506))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDNKSSL...NHE')",
            )
            self.assertEqual(hsp.query.id, "7")
            self.assertEqual(
                hsp.query.description,
                "gi|5817701|gb|AF142731.1|AF142731 Wisteria frutescens maturase-like protein (matK) gene, complete cds; chloroplast gene for chloroplast product",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(506))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "7:727..2244")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq('MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDKKSSL...NHE')",
            )
            self.assertEqual(hsp.target.id, "gi|5817760|gb|AAD52903.1|AF142732_1")
            self.assertEqual(hsp.target.name, "AAD52903")
            self.assertEqual(
                hsp.target.description,
                "maturase-like protein [Wisteria sinensis] >gi|171909136|gb|ACB58144.1| maturase K [Wisteria brachybotrys] >gi|171909138|gb|ACB58145.1| maturase K [Wisteria floribunda] >gi|171909140|gb|ACB58146.1| maturase K [Wisteria floribunda] >gi|171909142|gb|ACB58147.1| maturase K [Wisteria floribunda] >gi|171909150|gb|ACB58151.1| maturase K [Wisteria sinensis] >gi|171909152|gb|ACB58152.1| maturase K [Wisteria sinensis] >gi|171909154|gb|ACB58153.1| maturase K [Wisteria sinensis] >gi|171909156|gb|ACB58154.1| maturase K [Wisteria villosa] >gi|171909158|gb|ACB58155.1| maturase K [Wisteria villosa] >gi|171909160|gb|ACB58156.1| maturase K [Wisteria villosa]",
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYD KSSLLIVKRLITRMYQQNHLIISANDSNKNPFLGYN NFYSQIISDGFAVVVEIPFFLQLSSSLEEAEIVKSYHNLRSIHSIFPFLEDK TY NYVSDIRIPYPIHLEILVQILRYWVKDASFFHLLRFFLY+F NRNSLITPKKSISTFSKSNPRLFLFLYNFYVCEYESIFRFLRNQSSHLRLKSFSVFFERIFFYAKREHLVKVFPKDFSSTLTFFKDPFIHYVRYQGKSILASKNAPLLMNKWKHYFIHLWQCFFDVWSQPGTIHINQLSEHSFHFLGYFSNVRLNRSVVRSQMLQNTFLIEIVIKKLDIIVPIIPLIRSLAKAKFCNVLGHP+SKSVWADSSDFDIIDRFLRICRNLSHYYNGSSKKKNLYRIKYILRLSCIKTLACKHKSTVRAFLKKSGSEELLEEFFTEEEEILSLIFPR SSTLQRLHRNRIWYLDILFSNDLVNHE",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|581776         0 MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDKKSSLLIVKRL
                  0 |||||||||||||||||||||||||||||||||||||||||||||||||.||||||||||
7                 0 MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDNKSSLLIVKRL

gi|581776        60 ITRMYQQNHLIISANDSNKNPFLGYNNNFYSQIISDGFAVVVEIPFFLQLSSSLEEAEIV
                 60 ||||||||||||||||||||||||||.|||||||||||||||||||||||||||||||||
7                60 ITRMYQQNHLIISANDSNKNPFLGYNKNFYSQIISDGFAVVVEIPFFLQLSSSLEEAEIV

gi|581776       120 KSYHNLRSIHSIFPFLEDKLTYFNYVSDIRIPYPIHLEILVQILRYWVKDASFFHLLRFF
                120 |||||||||||||||||||.||.|||||||||||||||||||||||||||||........
7               120 KSYHNLRSIHSIFPFLEDKFTYLNYVSDIRIPYPIHLEILVQILRYWVKDASXXXXXXXX

gi|581776       180 LYNFCNRNSLITPKKSISTFSKSNPRLFLFLYNFYVCEYESIFRFLRNQSSHLRLKSFSV
                180 .....|||||||||||||||||||||||||||||||||||||||||||||||||||||||
7               180 XXXXSNRNSLITPKKSISTFSKSNPRLFLFLYNFYVCEYESIFRFLRNQSSHLRLKSFSV

gi|581776       240 FFERIFFYAKREHLVKVFPKDFSSTLTFFKDPFIHYVRYQGKSILASKNAPLLMNKWKHY
                240 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
7               240 FFERIFFYAKREHLVKVFPKDFSSTLTFFKDPFIHYVRYQGKSILASKNAPLLMNKWKHY

gi|581776       300 FIHLWQCFFDVWSQPGTIHINQLSEHSFHFLGYFSNVRLNRSVVRSQMLQNTFLIEIVIK
                300 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
7               300 FIHLWQCFFDVWSQPGTIHINQLSEHSFHFLGYFSNVRLNRSVVRSQMLQNTFLIEIVIK

gi|581776       360 KLDIIVPIIPLIRSLAKAKFCNVLGHPISKSVWADSSDFDIIDRFLRICRNLSHYYNGSS
                360 |||||||||||||||||||||||||||.||||||||||||||||||||||||||||||||
7               360 KLDIIVPIIPLIRSLAKAKFCNVLGHPLSKSVWADSSDFDIIDRFLRICRNLSHYYNGSS

gi|581776       420 KKKNLYRIKYILRLSCIKTLACKHKSTVRAFLKKSGSEELLEEFFTEEEEILSLIFPRAS
                420 ||||||||||||||||||||||||||||||||||||....................||.|
7               420 KKKNLYRIKYILRLSCIKTLACKHKSTVRAFLKKSGXXXXXXXXXXXXXXXXXXXXPRTS

gi|581776       480 STLQRLHRNRIWYLDILFSNDLVNHE 506
                480 |||||||||||||||||||||||||| 506
7               480 STLQRLHRNRIWYLDILFSNDLVNHE 506
""",
            )
            hit = record[5]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|171909134|gb|ACB58143.1|")
            self.assertEqual(hit.target.name, "ACB58143")
            self.assertEqual(
                hit.target.description, "maturase K [Wisteria brachybotrys]"
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=506)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 2398.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 928.316815329857)
            self.assertAlmostEqual(hsp.annotations["evalue"], 0.0)
            self.assertEqual(hsp.annotations["identity"], 496)
            self.assertEqual(hsp.annotations["positive"], 498)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[  0, 506],
                              [  0, 506]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 506))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDNKSSL...NHE')",
            )
            self.assertEqual(hsp.query.id, "7")
            self.assertEqual(
                hsp.query.description,
                "gi|5817701|gb|AF142731.1|AF142731 Wisteria frutescens maturase-like protein (matK) gene, complete cds; chloroplast gene for chloroplast product",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(506))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "7:727..2244")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq('MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDKKSSL...NHE')",
            )
            self.assertEqual(hsp.target.id, "gi|171909134|gb|ACB58143.1|")
            self.assertEqual(hsp.target.name, "ACB58143")
            self.assertEqual(
                hsp.target.description, "maturase K [Wisteria brachybotrys]"
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYD KSSLLIVKRLITRMYQQNHLIISANDSNKNPFLGYN  FYSQIISDGFAVVVEIPFFLQLSSSLEEAEIVKSYHNLRSIHSIFPFLEDK TY NYVSDIRIPYPIHLEILVQILRY VKDASFFHLLRFFLY+F NRNSLITPKKSISTFSKSNPRLFLFLYNFYVCEYESIFRFLRNQSSHLRLKSFSVFFERIFFYAKREHLVKVFPKDFSSTLTFFKDPFIHYVRYQGKSILASKNAPLLMNKWKHYFIHLWQCFFDVWSQPGTIHINQLSEHSFHFLGYFSNVRLNRSVVRSQMLQNTFLIEIVIKKLDIIVPIIPLIRSLAKAKFCNVLGHP+SKSVWADSSDFDIIDRFLRICRNLSHYYNGSSKKKNLYRIKYILRLSCIKTLACKHKSTVRAFLKKSGSEELLEEFFTEEEEILSLIFPR SSTLQRLHRNRIWYLDILFSNDLVNHE",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|171909         0 MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDKKSSLLIVKRL
                  0 |||||||||||||||||||||||||||||||||||||||||||||||||.||||||||||
7                 0 MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDNKSSLLIVKRL

gi|171909        60 ITRMYQQNHLIISANDSNKNPFLGYNNKFYSQIISDGFAVVVEIPFFLQLSSSLEEAEIV
                 60 ||||||||||||||||||||||||||..||||||||||||||||||||||||||||||||
7                60 ITRMYQQNHLIISANDSNKNPFLGYNKNFYSQIISDGFAVVVEIPFFLQLSSSLEEAEIV

gi|171909       120 KSYHNLRSIHSIFPFLEDKLTYFNYVSDIRIPYPIHLEILVQILRYRVKDASFFHLLRFF
                120 |||||||||||||||||||.||.|||||||||||||||||||||||.|||||........
7               120 KSYHNLRSIHSIFPFLEDKFTYLNYVSDIRIPYPIHLEILVQILRYWVKDASXXXXXXXX

gi|171909       180 LYNFCNRNSLITPKKSISTFSKSNPRLFLFLYNFYVCEYESIFRFLRNQSSHLRLKSFSV
                180 .....|||||||||||||||||||||||||||||||||||||||||||||||||||||||
7               180 XXXXSNRNSLITPKKSISTFSKSNPRLFLFLYNFYVCEYESIFRFLRNQSSHLRLKSFSV

gi|171909       240 FFERIFFYAKREHLVKVFPKDFSSTLTFFKDPFIHYVRYQGKSILASKNAPLLMNKWKHY
                240 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
7               240 FFERIFFYAKREHLVKVFPKDFSSTLTFFKDPFIHYVRYQGKSILASKNAPLLMNKWKHY

gi|171909       300 FIHLWQCFFDVWSQPGTIHINQLSEHSFHFLGYFSNVRLNRSVVRSQMLQNTFLIEIVIK
                300 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
7               300 FIHLWQCFFDVWSQPGTIHINQLSEHSFHFLGYFSNVRLNRSVVRSQMLQNTFLIEIVIK

gi|171909       360 KLDIIVPIIPLIRSLAKAKFCNVLGHPISKSVWADSSDFDIIDRFLRICRNLSHYYNGSS
                360 |||||||||||||||||||||||||||.||||||||||||||||||||||||||||||||
7               360 KLDIIVPIIPLIRSLAKAKFCNVLGHPLSKSVWADSSDFDIIDRFLRICRNLSHYYNGSS

gi|171909       420 KKKNLYRIKYILRLSCIKTLACKHKSTVRAFLKKSGSEELLEEFFTEEEEILSLIFPRAS
                420 ||||||||||||||||||||||||||||||||||||....................||.|
7               420 KKKNLYRIKYILRLSCIKTLACKHKSTVRAFLKKSGXXXXXXXXXXXXXXXXXXXXPRTS

gi|171909       480 STLQRLHRNRIWYLDILFSNDLVNHE 506
                480 |||||||||||||||||||||||||| 506
7               480 STLQRLHRNRIWYLDILFSNDLVNHE 506
""",
            )
            hit = record[6]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|5817761|gb|AAD52904.1|AF142733_1")
            self.assertEqual(hit.target.name, "AAD52904")
            self.assertEqual(
                hit.target.description, "maturase-like protein [Callerya reticulata]"
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=506)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 2390.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 925.235218722518)
            self.assertAlmostEqual(hsp.annotations["evalue"], 0.0)
            self.assertEqual(hsp.annotations["identity"], 493)
            self.assertEqual(hsp.annotations["positive"], 498)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[  0, 506],
                              [  0, 506]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 506))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDNKSSL...NHE')",
            )
            self.assertEqual(hsp.query.id, "7")
            self.assertEqual(
                hsp.query.description,
                "gi|5817701|gb|AF142731.1|AF142731 Wisteria frutescens maturase-like protein (matK) gene, complete cds; chloroplast gene for chloroplast product",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(506))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "7:727..2244")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq('MKEYQAYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDNKSSL...NHE')",
            )
            self.assertEqual(hsp.target.id, "gi|5817761|gb|AAD52904.1|AF142733_1")
            self.assertEqual(hsp.target.name, "AAD52904")
            self.assertEqual(
                hsp.target.description, "maturase-like protein [Callerya reticulata]"
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "MKEYQ YLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDNKSSLLIVKRLITRMYQQNHLIISANDSNKNPFLGYNKNFYSQIISDGFAVVVEIPFFLQL+SSLEEAEIVKSYHNLRSIHSIFPFLEDKFTYLNYVSDIRIPYPIHLEILVQILRYWVKDASFFHLLRFFLY+F NRNSLITPKKSISTFSK NPRLFLFLYNFYV EYESIF FLRNQSSHLR KSFSVFFERIFFYAKREHL+KVFPKDFSSTLTFFKDPFIHYVRYQ KSILASKNAPLLMNKWKHYFIHLWQCFFDVWSQPGTIHINQLSEHSFHFLGYFSNVRLNRSVVRSQMLQNTFLIEIVIKKLDIIVPIIPLIRSLAKAKFCNVLGHP+SKSVWADSSDFDIIDRFLRICRNLSHYYNGSSKKKNLYRIKYILRLSCIKTLACKHKSTVRAFLKKSGSEELLEEFFTEEEEILSLIFPR SSTL+RLHRNRIWYLDILFSNDLVNHE",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|581776         0 MKEYQAYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDNKSSLLIVKRL
                  0 |||||.||||||||||||||||||||||||||||||||||||||||||||||||||||||
7                 0 MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDNKSSLLIVKRL

gi|581776        60 ITRMYQQNHLIISANDSNKNPFLGYNKNFYSQIISDGFAVVVEIPFFLQLNSSLEEAEIV
                 60 ||||||||||||||||||||||||||||||||||||||||||||||||||.|||||||||
7                60 ITRMYQQNHLIISANDSNKNPFLGYNKNFYSQIISDGFAVVVEIPFFLQLSSSLEEAEIV

gi|581776       120 KSYHNLRSIHSIFPFLEDKFTYLNYVSDIRIPYPIHLEILVQILRYWVKDASFFHLLRFF
                120 ||||||||||||||||||||||||||||||||||||||||||||||||||||........
7               120 KSYHNLRSIHSIFPFLEDKFTYLNYVSDIRIPYPIHLEILVQILRYWVKDASXXXXXXXX

gi|581776       180 LYNFCNRNSLITPKKSISTFSKCNPRLFLFLYNFYVWEYESIFLFLRNQSSHLRFKSFSV
                180 .....|||||||||||||||||.|||||||||||||.||||||.||||||||||.|||||
7               180 XXXXSNRNSLITPKKSISTFSKSNPRLFLFLYNFYVCEYESIFRFLRNQSSHLRLKSFSV

gi|581776       240 FFERIFFYAKREHLLKVFPKDFSSTLTFFKDPFIHYVRYQEKSILASKNAPLLMNKWKHY
                240 ||||||||||||||.|||||||||||||||||||||||||.|||||||||||||||||||
7               240 FFERIFFYAKREHLVKVFPKDFSSTLTFFKDPFIHYVRYQGKSILASKNAPLLMNKWKHY

gi|581776       300 FIHLWQCFFDVWSQPGTIHINQLSEHSFHFLGYFSNVRLNRSVVRSQMLQNTFLIEIVIK
                300 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
7               300 FIHLWQCFFDVWSQPGTIHINQLSEHSFHFLGYFSNVRLNRSVVRSQMLQNTFLIEIVIK

gi|581776       360 KLDIIVPIIPLIRSLAKAKFCNVLGHPISKSVWADSSDFDIIDRFLRICRNLSHYYNGSS
                360 |||||||||||||||||||||||||||.||||||||||||||||||||||||||||||||
7               360 KLDIIVPIIPLIRSLAKAKFCNVLGHPLSKSVWADSSDFDIIDRFLRICRNLSHYYNGSS

gi|581776       420 KKKNLYRIKYILRLSCIKTLACKHKSTVRAFLKKSGSEELLEEFFTEEEEILSLIFPRAS
                420 ||||||||||||||||||||||||||||||||||||....................||.|
7               420 KKKNLYRIKYILRLSCIKTLACKHKSTVRAFLKKSGXXXXXXXXXXXXXXXXXXXXPRTS

gi|581776       480 STLKRLHRNRIWYLDILFSNDLVNHE 506
                480 |||.|||||||||||||||||||||| 506
7               480 STLQRLHRNRIWYLDILFSNDLVNHE 506
""",
            )
            hit = record[7]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|5817762|gb|AAD52905.1|AF142734_1")
            self.assertEqual(hit.target.name, "AAD52905")
            self.assertEqual(
                hit.target.description, "maturase-like protein [Callerya atropurpurea]"
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=506)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 2301.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 890.952456465874)
            self.assertAlmostEqual(hsp.annotations["evalue"], 0.0)
            self.assertEqual(hsp.annotations["identity"], 472)
            self.assertEqual(hsp.annotations["positive"], 488)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[  0, 506],
                              [  0, 506]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 506))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDNKSSL...NHE')",
            )
            self.assertEqual(hsp.query.id, "7")
            self.assertEqual(
                hsp.query.description,
                "gi|5817701|gb|AF142731.1|AF142731 Wisteria frutescens maturase-like protein (matK) gene, complete cds; chloroplast gene for chloroplast product",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(506))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "7:727..2244")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq('MKEYQVYLERDRSRQQDFLYPLIFREYTYGLAYSHDFNRSIFVENVGYDNKSSL...NHE')",
            )
            self.assertEqual(hsp.target.id, "gi|5817762|gb|AAD52905.1|AF142734_1")
            self.assertEqual(hsp.target.name, "AAD52905")
            self.assertEqual(
                hsp.target.description, "maturase-like protein [Callerya atropurpurea]"
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "MKEYQVYLERDRSRQQDFLYPLIFREY YGLAYSHDFNRSIFVENVGYDNKSSLLIVKRLITRMYQQNHLIIS NDSNKNPFLGYNKNFYSQIIS+ FA+V EIPFF QLSSSLE+AEIVKSYHNLRSIHSIFPFLEDKFTYLNYVSDIRIPYPIHLEILVQILRYWVKDA FFHLLR FLY+F N N++ TPKKSISTFS+SNPR FLFLYNFYVCEYESIF FLRN+SSHLRLKSFSVFFERIFFYAKREHLV+VF KDFSSTLTFFKDPFIHYVRYQGKSILASKNAPLLMNKWKHYFIHLWQ FFDVWSQPGTIHINQLSEHSFHFLGYFSNVRLNRSVVRSQMLQNTF+IEI IKKLDIIVPIIPLIRSLAKAKFCNVLGHP+SK VWADSSDFDII+RFLRICRNLSHYYNGSSKKK+LYRIKYILRLSCIKTLACKHKSTVRAFLK+ GSEELLEEFFTEEEEILSLIFPR SSTLQ+LHRNRIWYLDILF+NDLVNHE",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|581776         0 MKEYQVYLERDRSRQQDFLYPLIFREYTYGLAYSHDFNRSIFVENVGYDNKSSLLIVKRL
                  0 |||||||||||||||||||||||||||.||||||||||||||||||||||||||||||||
7                 0 MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDNKSSLLIVKRL

gi|581776        60 ITRMYQQNHLIISPNDSNKNPFLGYNKNFYSQIISEVFAIVAEIPFFRQLSSSLEQAEIV
                 60 |||||||||||||.|||||||||||||||||||||..||.|.|||||.|||||||.||||
7                60 ITRMYQQNHLIISANDSNKNPFLGYNKNFYSQIISDGFAVVVEIPFFLQLSSSLEEAEIV

gi|581776       120 KSYHNLRSIHSIFPFLEDKFTYLNYVSDIRIPYPIHLEILVQILRYWVKDAPFFHLLRLF
                120 |||||||||||||||||||||||||||||||||||||||||||||||||||.........
7               120 KSYHNLRSIHSIFPFLEDKFTYLNYVSDIRIPYPIHLEILVQILRYWVKDASXXXXXXXX

gi|581776       180 LYNFCNWNTVFTPKKSISTFSRSNPRFFLFLYNFYVCEYESIFLFLRNKSSHLRLKSFSV
                180 .....|.|...||||||||||.||||.||||||||||||||||.||||.|||||||||||
7               180 XXXXSNRNSLITPKKSISTFSKSNPRLFLFLYNFYVCEYESIFRFLRNQSSHLRLKSFSV

gi|581776       240 FFERIFFYAKREHLVEVFAKDFSSTLTFFKDPFIHYVRYQGKSILASKNAPLLMNKWKHY
                240 |||||||||||||||.||.|||||||||||||||||||||||||||||||||||||||||
7               240 FFERIFFYAKREHLVKVFPKDFSSTLTFFKDPFIHYVRYQGKSILASKNAPLLMNKWKHY

gi|581776       300 FIHLWQSFFDVWSQPGTIHINQLSEHSFHFLGYFSNVRLNRSVVRSQMLQNTFIIEIGIK
                300 ||||||.||||||||||||||||||||||||||||||||||||||||||||||.|||.||
7               300 FIHLWQCFFDVWSQPGTIHINQLSEHSFHFLGYFSNVRLNRSVVRSQMLQNTFLIEIVIK

gi|581776       360 KLDIIVPIIPLIRSLAKAKFCNVLGHPISKPVWADSSDFDIIERFLRICRNLSHYYNGSS
                360 |||||||||||||||||||||||||||.||.|||||||||||.|||||||||||||||||
7               360 KLDIIVPIIPLIRSLAKAKFCNVLGHPLSKSVWADSSDFDIIDRFLRICRNLSHYYNGSS

gi|581776       420 KKKSLYRIKYILRLSCIKTLACKHKSTVRAFLKRLGSEELLEEFFTEEEEILSLIFPRAS
                420 |||.|||||||||||||||||||||||||||||..|....................||.|
7               420 KKKNLYRIKYILRLSCIKTLACKHKSTVRAFLKKSGXXXXXXXXXXXXXXXXXXXXPRTS

gi|581776       480 STLQKLHRNRIWYLDILFTNDLVNHE 506
                480 ||||.|||||||||||||.||||||| 506
7               480 STLQRLHRNRIWYLDILFSNDLVNHE 506
""",
            )
            hit = record[8]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|152014012|gb|ABS20107.1|")
            self.assertEqual(hit.target.name, "ABS20107")
            self.assertEqual(
                hit.target.description, "maturase-like protein [Astragalus uliginosus]"
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=506)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 2293.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 887.870859858535)
            self.assertAlmostEqual(hsp.annotations["evalue"], 0.0)
            self.assertEqual(hsp.annotations["identity"], 470)
            self.assertEqual(hsp.annotations["positive"], 487)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[  0, 506],
                              [  0, 506]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 506))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDNKSSL...NHE')",
            )
            self.assertEqual(hsp.query.id, "7")
            self.assertEqual(
                hsp.query.description,
                "gi|5817701|gb|AF142731.1|AF142731 Wisteria frutescens maturase-like protein (matK) gene, complete cds; chloroplast gene for chloroplast product",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(506))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "7:727..2244")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq('MKEYQVFLERDRSRQQDFLYPLIFREYVYGLAYSHDFNRSTFVENVGYDNKYSL...NHE')",
            )
            self.assertEqual(hsp.target.id, "gi|152014012|gb|ABS20107.1|")
            self.assertEqual(hsp.target.name, "ABS20107")
            self.assertEqual(
                hsp.target.description, "maturase-like protein [Astragalus uliginosus]"
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "MKEYQV+LERDRSRQQDFLYPLIFREY+YGLAYSHDFNRS FVENVGYDNK SLLIVKRLITRMYQQNHLIISANDS KNPFLGYNKNFYSQIIS+GFA+VVEIPFFLQ SSSL+EAEIVKSY NLRSIHSIFPFLEDKF YLNYVSDIRIPYPIHLEILVQILRYWVKDA FFHLLR FLY+F NRNS +TPKKSISTFSKSNPRLFLFLYNFYVCEYESIF FLR +SSHLRLKSFSVFFERIFFYAKREHLV+VF KDFSSTLTFFKDP IHYVRYQGKSILASKNAPLLMNKWKHYFIHLW+CFFDVWSQPGTIHI QLSEHSF+ LGYFSNVRLNRSVVRSQMLQNTFLIEIV KKLDIIVPIIP+IRSLAKAKFCNVLGHP+SK+VWADSSDFDIIDRFLRICRNLSHYYNGSSKKK+LYRIKYILRLSCIKTLACKHKSTVRAFLK+SGSEELLEEFFTEEEEILSLIFPR SSTLQ+LH NRIWYLDILFSNDLVNHE",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|152014         0 MKEYQVFLERDRSRQQDFLYPLIFREYVYGLAYSHDFNRSTFVENVGYDNKYSLLIVKRL
                  0 ||||||.||||||||||||||||||||.||||||||||||.||||||||||.||||||||
7                 0 MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDNKSSLLIVKRL

gi|152014        60 ITRMYQQNHLIISANDSKKNPFLGYNKNFYSQIISEGFAIVVEIPFFLQFSSSLKEAEIV
                 60 |||||||||||||||||.|||||||||||||||||.|||.|||||||||.||||.|||||
7                60 ITRMYQQNHLIISANDSNKNPFLGYNKNFYSQIISDGFAVVVEIPFFLQLSSSLEEAEIV

gi|152014       120 KSYKNLRSIHSIFPFLEDKFPYLNYVSDIRIPYPIHLEILVQILRYWVKDAPFFHLLRLF
                120 |||.||||||||||||||||.||||||||||||||||||||||||||||||.........
7               120 KSYHNLRSIHSIFPFLEDKFTYLNYVSDIRIPYPIHLEILVQILRYWVKDASXXXXXXXX

gi|152014       180 LYNFCNRNSFLTPKKSISTFSKSNPRLFLFLYNFYVCEYESIFLFLRKKSSHLRLKSFSV
                180 .....||||..||||||||||||||||||||||||||||||||.|||..|||||||||||
7               180 XXXXSNRNSLITPKKSISTFSKSNPRLFLFLYNFYVCEYESIFRFLRNQSSHLRLKSFSV

gi|152014       240 FFERIFFYAKREHLVEVFAKDFSSTLTFFKDPLIHYVRYQGKSILASKNAPLLMNKWKHY
                240 |||||||||||||||.||.|||||||||||||.|||||||||||||||||||||||||||
7               240 FFERIFFYAKREHLVKVFPKDFSSTLTFFKDPFIHYVRYQGKSILASKNAPLLMNKWKHY

gi|152014       300 FIHLWECFFDVWSQPGTIHIKQLSEHSFYLLGYFSNVRLNRSVVRSQMLQNTFLIEIVSK
                300 |||||.||||||||||||||.|||||||..||||||||||||||||||||||||||||.|
7               300 FIHLWQCFFDVWSQPGTIHINQLSEHSFHFLGYFSNVRLNRSVVRSQMLQNTFLIEIVIK

gi|152014       360 KLDIIVPIIPIIRSLAKAKFCNVLGHPISKAVWADSSDFDIIDRFLRICRNLSHYYNGSS
                360 ||||||||||.||||||||||||||||.||.|||||||||||||||||||||||||||||
7               360 KLDIIVPIIPLIRSLAKAKFCNVLGHPLSKSVWADSSDFDIIDRFLRICRNLSHYYNGSS

gi|152014       420 KKKSLYRIKYILRLSCIKTLACKHKSTVRAFLKRSGSEELLEEFFTEEEEILSLIFPRAS
                420 |||.|||||||||||||||||||||||||||||.||....................||.|
7               420 KKKNLYRIKYILRLSCIKTLACKHKSTVRAFLKKSGXXXXXXXXXXXXXXXXXXXXPRTS

gi|152014       480 STLQKLHGNRIWYLDILFSNDLVNHE 506
                480 ||||.||.|||||||||||||||||| 506
7               480 STLQRLHRNRIWYLDILFSNDLVNHE 506
""",
            )
            hit = record[9]
            self.assertIsInstance(hit.target, SeqRecord)
            self.assertEqual(hit.target.id, "gi|146197442|dbj|BAF57483.1|")
            self.assertEqual(hit.target.name, "BAF57483")
            self.assertEqual(
                hit.target.description,
                "maturase [Glycyrrhiza uralensis] >gi|146197444|dbj|BAF57484.1| maturase [Glycyrrhiza glabra] >gi|152014018|gb|ABS20110.1| maturase-like protein [Glycyrrhiza pallidiflora]",
            )
            self.assertEqual(repr(hit.target.seq), "Seq(None, length=506)")
            self.assertEqual(len(hit), 1)
            hsp = hit[0]
            self.assertAlmostEqual(hsp.score, 2292.0)
            self.assertAlmostEqual(hsp.annotations["bit score"], 887.485660282618)
            self.assertAlmostEqual(hsp.annotations["evalue"], 0.0)
            self.assertEqual(hsp.annotations["identity"], 471)
            self.assertEqual(hsp.annotations["positive"], 489)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[  0, 506],
                              [  0, 506]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 506))
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq('MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDNKSSL...NHE')",
            )
            self.assertEqual(hsp.query.id, "7")
            self.assertEqual(
                hsp.query.description,
                "gi|5817701|gb|AF142731.1|AF142731 Wisteria frutescens maturase-like protein (matK) gene, complete cds; chloroplast gene for chloroplast product",
            )
            self.assertEqual(len(hsp.query.features), 1)
            feature = hsp.query.features[0]
            self.assertEqual(
                repr(feature.location),
                "SimpleLocation(ExactPosition(0), ExactPosition(506))",
            )
            self.assertEqual(feature.type, "CDS")
            self.assertEqual(len(feature.qualifiers), 1)
            self.assertEqual(feature.qualifiers["coded_by"], "7:727..2244")
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq('MKEYQVYLERDRSRQQDFLYPLIFREYIYGIAYSHNLNRSIFVENVGYDNKFSL...NNE')",
            )
            self.assertEqual(hsp.target.id, "gi|146197442|dbj|BAF57483.1|")
            self.assertEqual(hsp.target.name, "BAF57483")
            self.assertEqual(
                hsp.target.description,
                "maturase [Glycyrrhiza uralensis] >gi|146197444|dbj|BAF57484.1| maturase [Glycyrrhiza glabra] >gi|152014018|gb|ABS20110.1| maturase-like protein [Glycyrrhiza pallidiflora]",
            )
            self.assertEqual(len(hsp.target.features), 0)
            self.assertEqual(
                hsp.annotations["midline"],
                "MKEYQVYLERDRSRQQDFLYPLIFREYIYG+AYSH+ NRSIFVENVGYDNK SLLIVKRLITRMYQQNHLIISANDSNKNPF GYNKN YSQ+ISDGFAVVVEIPFFLQ SSSLEEAEIVKSY+NLRSIHSIFPFLEDKFTYLNYVSDIRIPYPIHLEILVQILRYWVKDA FFHLLR FLY+F N NSLITPKKSISTFSKSNPRLFLFLYNFYVCEYESIF FLRN+SSHLRLKSFSVFFERIFFYAKREHLV VF KD+S TLT FKDPFIHYVRYQGK+ILAS+NAPLLMNKWKHYFIHLWQCFFDVWSQPGTIHINQLSEHSFHFLGYFSNVRLNRSVVRSQMLQNTFLIEIVIKKLDIIVPIIPLIRSLAKAKFCNVLGHP+SK VWADSSDF+II+RFLRICRNLSHYY+GSSKKK+LYRIKYILRLSCIKTLACKHKSTVRAFLK+ GSEELLEEFFTEEEEILSLIFP+ SSTLQ+LHRNRIWYLDILFSNDLVN+E",
            )
            self.assertEqual(
                str(hsp),
                """\
gi|146197         0 MKEYQVYLERDRSRQQDFLYPLIFREYIYGIAYSHNLNRSIFVENVGYDNKFSLLIVKRL
                  0 ||||||||||||||||||||||||||||||.||||..||||||||||||||.||||||||
7                 0 MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYDNKSSLLIVKRL

gi|146197        60 ITRMYQQNHLIISANDSNKNPFSGYNKNIYSQLISDGFAVVVEIPFFLQFSSSLEEAEIV
                 60 ||||||||||||||||||||||.|||||.|||.||||||||||||||||.||||||||||
7                60 ITRMYQQNHLIISANDSNKNPFLGYNKNFYSQIISDGFAVVVEIPFFLQLSSSLEEAEIV

gi|146197       120 KSYNNLRSIHSIFPFLEDKFTYLNYVSDIRIPYPIHLEILVQILRYWVKDAPFFHLLRLF
                120 |||.|||||||||||||||||||||||||||||||||||||||||||||||.........
7               120 KSYHNLRSIHSIFPFLEDKFTYLNYVSDIRIPYPIHLEILVQILRYWVKDASXXXXXXXX

gi|146197       180 LYNFCNWNSLITPKKSISTFSKSNPRLFLFLYNFYVCEYESIFLFLRNKSSHLRLKSFSV
                180 .....|.||||||||||||||||||||||||||||||||||||.||||.|||||||||||
7               180 XXXXSNRNSLITPKKSISTFSKSNPRLFLFLYNFYVCEYESIFRFLRNQSSHLRLKSFSV

gi|146197       240 FFERIFFYAKREHLVDVFAKDYSPTLTLFKDPFIHYVRYQGKAILASRNAPLLMNKWKHY
                240 |||||||||||||||.||.||.|.|||.||||||||||||||.||||.||||||||||||
7               240 FFERIFFYAKREHLVKVFPKDFSSTLTFFKDPFIHYVRYQGKSILASKNAPLLMNKWKHY

gi|146197       300 FIHLWQCFFDVWSQPGTIHINQLSEHSFHFLGYFSNVRLNRSVVRSQMLQNTFLIEIVIK
                300 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
7               300 FIHLWQCFFDVWSQPGTIHINQLSEHSFHFLGYFSNVRLNRSVVRSQMLQNTFLIEIVIK

gi|146197       360 KLDIIVPIIPLIRSLAKAKFCNVLGHPISKPVWADSSDFEIIERFLRICRNLSHYYSGSS
                360 |||||||||||||||||||||||||||.||.||||||||.||.|||||||||||||.|||
7               360 KLDIIVPIIPLIRSLAKAKFCNVLGHPLSKSVWADSSDFDIIDRFLRICRNLSHYYNGSS

gi|146197       420 KKKSLYRIKYILRLSCIKTLACKHKSTVRAFLKRLGSEELLEEFFTEEEEILSLIFPKAS
                420 |||.|||||||||||||||||||||||||||||..|....................|..|
7               420 KKKNLYRIKYILRLSCIKTLACKHKSTVRAFLKKSGXXXXXXXXXXXXXXXXXXXXPRTS

gi|146197       480 STLQKLHRNRIWYLDILFSNDLVNNE 506
                480 ||||.|||||||||||||||||||.| 506
7               480 STLQRLHRNRIWYLDILFSNDLVNHE 506
""",
            )

    def test_xml_2900_blastx_001(self):
        """Parsing BLASTX 2.9.0+ (xml_2900_blastx_001.xml)."""
        filename = "xml_2900_blastx_001.xml"
        datafile = os.path.join("Blast", filename)
        with open(datafile, "rb") as handle:
            records = Blast.parse(handle)
            self.check_xml_2900_blastx_001_records(records)
        with Blast.parse(datafile) as records:
            self.check_xml_2900_blastx_001_records(records)
        with open(datafile, "rb") as handle:
            record = Blast.read(handle)
        self.check_xml_2900_blastx_001_record(record)
        record = Blast.read(datafile)
        self.check_xml_2900_blastx_001_record(record)

    def check_xml_2900_blastx_001_records(self, records):
        self.assertEqual(records.program, "blastx")
        self.assertEqual(records.version, "BLASTX 2.9.0+")
        self.assertEqual(
            records.reference,
            'Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.',
        )
        self.assertEqual(records.db, "nr")
        self.assertIsInstance(records.query, SeqRecord)
        self.assertEqual(records.query.id, "AI021773.1")
        self.assertEqual(
            records.query.description,
            "MAAD0534.RAR Schistosoma mansoni, adult worm (J.C.Parra) Schistosoma mansoni cDNA clone MAAD0534.RAR 5' end similar to S. mansoni actin mRNA, complete cds, mRNA sequence",
        )
        self.assertEqual(repr(records.query.seq), "Seq(None, length=365)")
        self.assertEqual(len(records.param), 5)
        self.assertEqual(records.param["matrix"], "BLOSUM62")
        self.assertAlmostEqual(records.param["expect"], 10.0)
        self.assertEqual(records.param["gap-open"], 11)
        self.assertEqual(records.param["gap-extend"], 1)
        self.assertEqual(records.param["filter"], "F")
        record = next(records)
        self.assertRaises(StopIteration, next, records)
        self.check_xml_2900_blastx_001_record(record)

    def check_xml_2900_blastx_001_record(self, record):
        self.assertEqual(record.num, 1)

        self.assertIsInstance(record.query, SeqRecord)
        self.assertEqual(record.query.id, "AI021773.1")
        self.assertEqual(
            record.query.description,
            "MAAD0534.RAR Schistosoma mansoni, adult worm (J.C.Parra) Schistosoma mansoni cDNA clone MAAD0534.RAR 5' end similar to S. mansoni actin mRNA, complete cds, mRNA sequence",
        )
        self.assertEqual(repr(record.query.seq), "Seq(None, length=365)")

        self.assertEqual(len(record.stat), 7)
        self.assertEqual(record.stat["db-num"], 197469652)
        self.assertEqual(record.stat["db-len"], 71133367251)
        self.assertEqual(record.stat["hsp-len"], 0)
        self.assertAlmostEqual(record.stat["eff-space"], 0.0)
        self.assertAlmostEqual(record.stat["kappa"], 0.041)
        self.assertAlmostEqual(record.stat["lambda"], 0.267)
        self.assertAlmostEqual(record.stat["entropy"], 0.14)
        self.assertEqual(len(record), 10)

        hit = record[0]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|1530504495|emb|VDM03167.1|")
        self.assertEqual(hit.target.name, "VDM03167")
        self.assertEqual(
            hit.target.description,
            "unnamed protein product, partial [Schistocephalus solidus]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=132)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 408.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 161.77)
        self.assertAlmostEqual(hsp.annotations["evalue"], 8.11609e-49)
        self.assertEqual(hsp.annotations["identity"], 81)
        self.assertEqual(hsp.annotations["positive"], 83)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0, 108],
                          [  0, 108]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 108))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSY...LTE')",
        )
        self.assertEqual(hsp.query.id, "AI021773.1")
        self.assertEqual(
            hsp.query.description,
            "MAAD0534.RAR Schistosoma mansoni, adult worm (J.C.Parra) Schistosoma mansoni cDNA clone MAAD0534.RAR 5' end similar to S. mansoni actin mRNA, complete cds, mRNA sequence",
        )
        self.assertEqual(len(hsp.query.features), 1)
        feature = hsp.query.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(108))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(feature.qualifiers["coded_by"], "AI021773.1:20..343")
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MADEEVQALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSY...LTE'}, length=132)",
        )
        self.assertEqual(hsp.target.id, "gi|1530504495|emb|VDM03167.1|")
        self.assertEqual(hsp.target.name, "VDM03167")
        self.assertEqual(
            hsp.target.description,
            "unnamed protein product, partial [Schistocephalus solidus]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MADEEVQALVVDNGSGMCKAG       ++  P               G KDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|153050         0 MADEEVQALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQ
                  0 |||||||||||||||||||||...........|...............|.||||||||||
AI021773.         0 MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQ

gi|153050        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE 108
                 60 |||||||||||||||||||||||||||||||||||||||||||||||| 108
AI021773.        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE 108
""",
        )
        hit = record[1]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|510859078|gb|EPB74633.1|")
        self.assertEqual(hit.target.name, "EPB74633")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein ANCCEY_06263 [Ancylostoma ceylanicum]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=119)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 405.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 160.614)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.40046e-48)
        self.assertEqual(hsp.annotations["identity"], 81)
        self.assertEqual(hsp.annotations["positive"], 85)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0, 115],
                          [  0, 115]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 115))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSY...RKP')",
        )
        self.assertEqual(hsp.query.id, "AI021773.1")
        self.assertEqual(
            hsp.query.description,
            "MAAD0534.RAR Schistosoma mansoni, adult worm (J.C.Parra) Schistosoma mansoni cDNA clone MAAD0534.RAR 5' end similar to S. mansoni actin mRNA, complete cds, mRNA sequence",
        )
        self.assertEqual(len(hsp.query.features), 1)
        feature = hsp.query.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(115))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(feature.qualifiers["coded_by"], "AI021773.1:20..364")
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MCDDDVAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSY...LKP'}, length=119)",
        )
        self.assertEqual(hsp.target.id, "gi|510859078|gb|EPB74633.1|")
        self.assertEqual(hsp.target.name, "EPB74633")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein ANCCEY_06263 [Ancylostoma ceylanicum]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "M D++V ALVVDNGSGMCKAG       ++  P               G KDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE H I KP",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|510859         0 MCDDDVAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQ
                  0 |.|..|.||||||||||||||...........|...............|.||||||||||
AI021773.         0 MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQ

gi|510859        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTEAHSILKP 115
                 60 ||||||||||||||||||||||||||||||||||||||||||||||||.|.|.|| 115
AI021773.        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTELHCIRKP 115
""",
        )
        hit = record[2]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|684409690|ref|XP_009175831.1|")
        self.assertEqual(hit.target.name, "XP_009175831")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein T265_11027 [Opisthorchis viverrini] >gi|663044098|gb|KER20427.1| hypothetical protein T265_11027 [Opisthorchis viverrini]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=246)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 413.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 163.696)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.40404e-48)
        self.assertEqual(hsp.annotations["identity"], 81)
        self.assertEqual(hsp.annotations["positive"], 83)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0, 108],
                          [  0, 108]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 108))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSY...LTE')",
        )
        self.assertEqual(hsp.query.id, "AI021773.1")
        self.assertEqual(
            hsp.query.description,
            "MAAD0534.RAR Schistosoma mansoni, adult worm (J.C.Parra) Schistosoma mansoni cDNA clone MAAD0534.RAR 5' end similar to S. mansoni actin mRNA, complete cds, mRNA sequence",
        )
        self.assertEqual(len(hsp.query.features), 1)
        feature = hsp.query.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(108))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(feature.qualifiers["coded_by"], "AI021773.1:20..343")
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MADEEVQALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSY...LTE'}, length=246)",
        )
        self.assertEqual(hsp.target.id, "gi|684409690|ref|XP_009175831.1|")
        self.assertEqual(hsp.target.name, "XP_009175831")
        self.assertEqual(
            hsp.target.description,
            "hypothetical protein T265_11027 [Opisthorchis viverrini] >gi|663044098|gb|KER20427.1| hypothetical protein T265_11027 [Opisthorchis viverrini]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MADEEVQALVVDNGSGMCKAG       ++  P               G KDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|684409         0 MADEEVQALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQ
                  0 |||||||||||||||||||||...........|...............|.||||||||||
AI021773.         0 MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQ

gi|684409        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE 108
                 60 |||||||||||||||||||||||||||||||||||||||||||||||| 108
AI021773.        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE 108
""",
        )
        hit = record[3]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|449710331|gb|EMD49430.1|")
        self.assertEqual(hit.target.name, "EMD49430")
        self.assertEqual(
            hit.target.description,
            "actin, putative, partial [Entamoeba histolytica KU27]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=124)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 401.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 159.073)
        self.assertAlmostEqual(hsp.annotations["evalue"], 9.0486e-48)
        self.assertEqual(hsp.annotations["identity"], 78)
        self.assertEqual(hsp.annotations["positive"], 81)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0, 108],
                          [  0, 108]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 108))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSY...LTE')",
        )
        self.assertEqual(hsp.query.id, "AI021773.1")
        self.assertEqual(
            hsp.query.description,
            "MAAD0534.RAR Schistosoma mansoni, adult worm (J.C.Parra) Schistosoma mansoni cDNA clone MAAD0534.RAR 5' end similar to S. mansoni actin mRNA, complete cds, mRNA sequence",
        )
        self.assertEqual(len(hsp.query.features), 1)
        feature = hsp.query.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(108))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(feature.qualifiers["coded_by"], "AI021773.1:20..343")
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGDEEVQALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHVSVMAGMGQKDAY...LTE'}, length=124)",
        )
        self.assertEqual(hsp.target.id, "gi|449710331|gb|EMD49430.1|")
        self.assertEqual(hsp.target.name, "EMD49430")
        self.assertEqual(
            hsp.target.description,
            "actin, putative, partial [Entamoeba histolytica KU27]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "M DEEVQALVVDNGSGMCKAG       ++  P               G KD+YVGDEAQSKRGILTLKYPIEHGIV NWDDMEKIWHHTFYNELRVAPEEHPVLLTE",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|449710         0 MGDEEVQALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHVSVMAGMGQKDAYVGDEAQ
                  0 |.|||||||||||||||||||...........|...............|.||.|||||||
AI021773.         0 MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQ

gi|449710        60 SKRGILTLKYPIEHGIVNNWDDMEKIWHHTFYNELRVAPEEHPVLLTE 108
                 60 |||||||||||||||||.|||||||||||||||||||||||||||||| 108
AI021773.        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE 108
""",
        )
        hit = record[4]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|257215766|emb|CAX83035.1|")
        self.assertEqual(hit.target.name, "CAX83035")
        self.assertEqual(
            hit.target.description, "Actin-2, partial [Schistosoma japonicum]"
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=252)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 411.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 162.925)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.00219e-47)
        self.assertEqual(hsp.annotations["identity"], 81)
        self.assertEqual(hsp.annotations["positive"], 83)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0, 108],
                          [  0, 108]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 108))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSY...LTE')",
        )
        self.assertEqual(hsp.query.id, "AI021773.1")
        self.assertEqual(
            hsp.query.description,
            "MAAD0534.RAR Schistosoma mansoni, adult worm (J.C.Parra) Schistosoma mansoni cDNA clone MAAD0534.RAR 5' end similar to S. mansoni actin mRNA, complete cds, mRNA sequence",
        )
        self.assertEqual(len(hsp.query.features), 1)
        feature = hsp.query.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(108))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(feature.qualifiers["coded_by"], "AI021773.1:20..343")
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MADEEVQALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSY...LTE'}, length=252)",
        )
        self.assertEqual(hsp.target.id, "gi|257215766|emb|CAX83035.1|")
        self.assertEqual(hsp.target.name, "CAX83035")
        self.assertEqual(
            hsp.target.description, "Actin-2, partial [Schistosoma japonicum]"
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MADEEVQALVVDNGSGMCKAG       ++  P               G KDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|257215         0 MADEEVQALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQ
                  0 |||||||||||||||||||||...........|...............|.||||||||||
AI021773.         0 MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQ

gi|257215        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE 108
                 60 |||||||||||||||||||||||||||||||||||||||||||||||| 108
AI021773.        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE 108
""",
        )
        hit = record[5]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|1535393712|emb|VDP83060.1|")
        self.assertEqual(hit.target.name, "VDP83060")
        self.assertEqual(
            hit.target.description,
            "unnamed protein product, partial [Echinostoma caproni]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=209)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 407.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 161.384)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.16397e-47)
        self.assertEqual(hsp.annotations["identity"], 80)
        self.assertEqual(hsp.annotations["positive"], 83)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0, 108],
                          [  0, 108]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 108))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSY...LTE')",
        )
        self.assertEqual(hsp.query.id, "AI021773.1")
        self.assertEqual(
            hsp.query.description,
            "MAAD0534.RAR Schistosoma mansoni, adult worm (J.C.Parra) Schistosoma mansoni cDNA clone MAAD0534.RAR 5' end similar to S. mansoni actin mRNA, complete cds, mRNA sequence",
        )
        self.assertEqual(len(hsp.query.features), 1)
        feature = hsp.query.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(108))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(feature.qualifiers["coded_by"], "AI021773.1:20..343")
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MADDEVQALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSY...LTE'}, length=209)",
        )
        self.assertEqual(hsp.target.id, "gi|1535393712|emb|VDP83060.1|")
        self.assertEqual(hsp.target.name, "VDP83060")
        self.assertEqual(
            hsp.target.description,
            "unnamed protein product, partial [Echinostoma caproni]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MAD+EVQALVVDNGSGMCKAG       ++  P               G KDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|153539         0 MADDEVQALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQ
                  0 |||.|||||||||||||||||...........|...............|.||||||||||
AI021773.         0 MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQ

gi|153539        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE 108
                 60 |||||||||||||||||||||||||||||||||||||||||||||||| 108
AI021773.        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE 108
""",
        )
        hit = record[6]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|312773|emb|CAA50205.1|")
        self.assertEqual(hit.target.name, "CAA50205")
        self.assertEqual(
            hit.target.description, "actin, partial [Entamoeba histolytica]"
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=137)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 401.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 159.073)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.25869e-47)
        self.assertEqual(hsp.annotations["identity"], 78)
        self.assertEqual(hsp.annotations["positive"], 81)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0, 108],
                          [  0, 108]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 108))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSY...LTE')",
        )
        self.assertEqual(hsp.query.id, "AI021773.1")
        self.assertEqual(
            hsp.query.description,
            "MAAD0534.RAR Schistosoma mansoni, adult worm (J.C.Parra) Schistosoma mansoni cDNA clone MAAD0534.RAR 5' end similar to S. mansoni actin mRNA, complete cds, mRNA sequence",
        )
        self.assertEqual(len(hsp.query.features), 1)
        feature = hsp.query.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(108))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(feature.qualifiers["coded_by"], "AI021773.1:20..343")
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGDEEVQALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHVSVMAGMGQKDAY...LTE'}, length=137)",
        )
        self.assertEqual(hsp.target.id, "gi|312773|emb|CAA50205.1|")
        self.assertEqual(hsp.target.name, "CAA50205")
        self.assertEqual(
            hsp.target.description, "actin, partial [Entamoeba histolytica]"
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "M DEEVQALVVDNGSGMCKAG       ++  P               G KD+YVGDEAQSKRGILTLKYPIEHGIV NWDDMEKIWHHTFYNELRVAPEEHPVLLTE",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|312773         0 MGDEEVQALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHVSVMAGMGQKDAYVGDEAQ
                  0 |.|||||||||||||||||||...........|...............|.||.|||||||
AI021773.         0 MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQ

gi|312773        60 SKRGILTLKYPIEHGIVNNWDDMEKIWHHTFYNELRVAPEEHPVLLTE 108
                 60 |||||||||||||||||.|||||||||||||||||||||||||||||| 108
AI021773.        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE 108
""",
        )
        hit = record[7]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|1530341495|emb|VDN44756.1|")
        self.assertEqual(hit.target.name, "VDN44756")
        self.assertEqual(
            hit.target.description,
            "unnamed protein product, partial [Dibothriocephalus latus]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=145)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 400.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 158.688)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.78336e-47)
        self.assertEqual(hsp.annotations["identity"], 78)
        self.assertEqual(hsp.annotations["positive"], 82)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0, 108],
                          [  0, 108]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 108))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSY...LTE')",
        )
        self.assertEqual(hsp.query.id, "AI021773.1")
        self.assertEqual(
            hsp.query.description,
            "MAAD0534.RAR Schistosoma mansoni, adult worm (J.C.Parra) Schistosoma mansoni cDNA clone MAAD0534.RAR 5' end similar to S. mansoni actin mRNA, complete cds, mRNA sequence",
        )
        self.assertEqual(len(hsp.query.features), 1)
        feature = hsp.query.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(108))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(feature.qualifiers["coded_by"], "AI021773.1:20..343")
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MGDEDVQALVIDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSY...LTE'}, length=145)",
        )
        self.assertEqual(hsp.target.id, "gi|1530341495|emb|VDN44756.1|")
        self.assertEqual(hsp.target.name, "VDN44756")
        self.assertEqual(
            hsp.target.description,
            "unnamed protein product, partial [Dibothriocephalus latus]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "M DE+VQALV+DNGSGMCKAG       ++  P               G KDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|153034         0 MGDEDVQALVIDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQ
                  0 |.||.|||||.||||||||||...........|...............|.||||||||||
AI021773.         0 MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQ

gi|153034        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE 108
                 60 |||||||||||||||||||||||||||||||||||||||||||||||| 108
AI021773.        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE 108
""",
        )
        hit = record[8]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|1524877828|ref|XP_027046469.1|")
        self.assertEqual(hit.target.name, "XP_027046469")
        self.assertEqual(
            hit.target.description, "actin-1, partial [Pocillopora damicornis]"
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=122)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 398.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 157.918)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.93331e-47)
        self.assertEqual(hsp.annotations["identity"], 78)
        self.assertEqual(hsp.annotations["positive"], 82)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0, 108],
                          [  0, 108]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 108))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSY...LTE')",
        )
        self.assertEqual(hsp.query.id, "AI021773.1")
        self.assertEqual(
            hsp.query.description,
            "MAAD0534.RAR Schistosoma mansoni, adult worm (J.C.Parra) Schistosoma mansoni cDNA clone MAAD0534.RAR 5' end similar to S. mansoni actin mRNA, complete cds, mRNA sequence",
        )
        self.assertEqual(len(hsp.query.features), 1)
        feature = hsp.query.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(108))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(feature.qualifiers["coded_by"], "AI021773.1:20..343")
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MADEEVAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSY...LTE'}, length=122)",
        )
        self.assertEqual(hsp.target.id, "gi|1524877828|ref|XP_027046469.1|")
        self.assertEqual(hsp.target.name, "XP_027046469")
        self.assertEqual(
            hsp.target.description, "actin-1, partial [Pocillopora damicornis]"
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MADEEV ALVVDNGSGMCKAG       ++  P               G KDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELR+APEEHP+LLTE",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|152487         0 MADEEVAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQ
                  0 ||||||.||||||||||||||...........|...............|.||||||||||
AI021773.         0 MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQ

gi|152487        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRIAPEEHPILLTE 108
                 60 ||||||||||||||||||||||||||||||||||||.||||||.|||| 108
AI021773.        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE 108
""",
        )
        hit = record[9]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|1524877860|ref|XP_027046487.1|")
        self.assertEqual(hit.target.name, "XP_027046487")
        self.assertEqual(
            hit.target.description, "actin-1-like [Pocillopora damicornis]"
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=134)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 399.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 158.303)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.36088e-47)
        self.assertEqual(hsp.annotations["identity"], 79)
        self.assertEqual(hsp.annotations["positive"], 82)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0, 108],
                          [  0, 108]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 108))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSY...LTE')",
        )
        self.assertEqual(hsp.query.id, "AI021773.1")
        self.assertEqual(
            hsp.query.description,
            "MAAD0534.RAR Schistosoma mansoni, adult worm (J.C.Parra) Schistosoma mansoni cDNA clone MAAD0534.RAR 5' end similar to S. mansoni actin mRNA, complete cds, mRNA sequence",
        )
        self.assertEqual(len(hsp.query.features), 1)
        feature = hsp.query.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(108))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(feature.qualifiers["coded_by"], "AI021773.1:20..343")
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({0: 'MADEDVAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSY...LTE'}, length=134)",
        )
        self.assertEqual(hsp.target.id, "gi|1524877860|ref|XP_027046487.1|")
        self.assertEqual(hsp.target.name, "XP_027046487")
        self.assertEqual(
            hsp.target.description, "actin-1-like [Pocillopora damicornis]"
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MADE+V ALVVDNGSGMCKAG       ++  P               G KDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|152487         0 MADEDVAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQ
                  0 ||||.|.||||||||||||||...........|...............|.||||||||||
AI021773.         0 MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQ

gi|152487        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE 108
                 60 |||||||||||||||||||||||||||||||||||||||||||||||| 108
AI021773.        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE 108
""",
        )


class TestTBlastn(unittest.TestCase):
    """Test the Blast XML parser for tblastn output."""

    def test_xml_2900_tblastn_001(self):
        """Parsing TBLASTN 2.9.0+ (xml_2900_tblastn_001.xml)."""
        filename = "xml_2900_tblastn_001.xml"
        datafile = os.path.join("Blast", filename)
        with open(datafile, "rb") as handle:
            records = Blast.parse(handle)
            self.check_xml_2900_tblastn_001_records(records)
        with Blast.parse(datafile) as records:
            self.check_xml_2900_tblastn_001_records(records)
        with open(datafile, "rb") as handle:
            record = Blast.read(handle)
        self.check_xml_2900_tblastn_001_record(record)
        record = Blast.read(datafile)
        self.check_xml_2900_tblastn_001_record(record)

    def check_xml_2900_tblastn_001_records(self, records):
        self.assertEqual(records.program, "tblastn")
        self.assertEqual(records.version, "TBLASTN 2.9.0+")
        self.assertEqual(
            records.reference,
            'Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.',
        )
        self.assertEqual(records.db, "nr")
        self.assertIsInstance(records.query, SeqRecord)
        self.assertEqual(records.query.id, "CAJ99216.1")
        self.assertEqual(
            records.query.description, "tim [Helicobacter acinonychis str. Sheeba]"
        )
        self.assertEqual(repr(records.query.seq), "Seq(None, length=234)")
        self.assertEqual(len(records.param), 5)
        self.assertEqual(records.param["matrix"], "BLOSUM62")
        self.assertAlmostEqual(records.param["expect"], 10.0)
        self.assertEqual(records.param["gap-open"], 11)
        self.assertEqual(records.param["gap-extend"], 1)
        self.assertEqual(records.param["filter"], "F")
        record = next(records)
        self.assertRaises(StopIteration, next, records)
        self.check_xml_2900_tblastn_001_record(record)

    def check_xml_2900_tblastn_001_record(self, record):
        self.assertEqual(record.num, 1)

        self.assertIsInstance(record.query, SeqRecord)
        self.assertEqual(record.query.id, "CAJ99216.1")
        self.assertEqual(
            record.query.description, "tim [Helicobacter acinonychis str. Sheeba]"
        )
        self.assertEqual(repr(record.query.seq), "Seq(None, length=234)")

        self.assertEqual(len(record.stat), 7)
        self.assertEqual(record.stat["db-num"], 51302476)
        self.assertEqual(record.stat["db-len"], 204305995260)
        self.assertEqual(record.stat["hsp-len"], 0)
        self.assertAlmostEqual(record.stat["eff-space"], 0.0)
        self.assertAlmostEqual(record.stat["kappa"], 0.041)
        self.assertAlmostEqual(record.stat["lambda"], 0.267)
        self.assertAlmostEqual(record.stat["entropy"], 0.14)
        self.assertEqual(len(record), 10)

        hit = record[0]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|109713861|emb|AM260522.1|")
        self.assertEqual(hit.target.name, "AM260522")
        self.assertEqual(
            hit.target.description,
            "Helicobacter acinonychis str. Sheeba complete genome, strain Sheeba",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=1553927)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 1228.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 477.633)
        self.assertAlmostEqual(hsp.annotations["evalue"], 8.53508e-152)
        self.assertEqual(hsp.annotations["identity"], 234)
        self.assertEqual(hsp.annotations["positive"], 234)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0, 234],
                          [  0, 234]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 234))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHF...SFL')",
        )
        self.assertEqual(hsp.query.id, "CAJ99216.1")
        self.assertEqual(
            hsp.query.description, "tim [Helicobacter acinonychis str. Sheeba]"
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHF...SFL')",
        )
        self.assertEqual(hsp.target.id, "gi|109713861|emb|AM260522.1|")
        self.assertEqual(hsp.target.name, "AM260522")
        self.assertEqual(
            hsp.target.description,
            "Helicobacter acinonychis str. Sheeba complete genome, strain Sheeba",
        )
        self.assertEqual(len(hsp.target.features), 1)
        feature = hsp.target.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(1553927))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(
            feature.qualifiers["coded_by"],
            "gi|109713861|emb|AM260522.1|:325802..326503",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQNAYPKDCGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPNFLKEKFDFFKDKKFKIVYCIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHGFLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|109713         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQN
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
CAJ99216.         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQN

gi|109713        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPNFLKEKFDFFKDKKFKIVY
                 60 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
CAJ99216.        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPNFLKEKFDFFKDKKFKIVY

gi|109713       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHG
                120 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
CAJ99216.       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHG

gi|109713       180 FLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL 234
                180 |||||||||||||||||||||||||||||||||||||||||||||||||||||| 234
CAJ99216.       180 FLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL 234
""",
        )
        hit = record[1]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|1336928286|emb|LT900055.1|")
        self.assertEqual(hit.target.name, "LT900055")
        self.assertEqual(
            hit.target.description,
            "Helicobacter acinonychis isolate 212_9 genome assembly, chromosome: I",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=1550239)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 1219.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 474.167)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.37046e-150)
        self.assertEqual(hsp.annotations["identity"], 232)
        self.assertEqual(hsp.annotations["positive"], 232)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0, 234],
                          [  0, 234]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 234))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHF...SFL')",
        )
        self.assertEqual(hsp.query.id, "CAJ99216.1")
        self.assertEqual(
            hsp.query.description, "tim [Helicobacter acinonychis str. Sheeba]"
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHF...SFL')",
        )
        self.assertEqual(hsp.target.id, "gi|1336928286|emb|LT900055.1|")
        self.assertEqual(hsp.target.name, "LT900055")
        self.assertEqual(
            hsp.target.description,
            "Helicobacter acinonychis isolate 212_9 genome assembly, chromosome: I",
        )
        self.assertEqual(len(hsp.target.features), 1)
        feature = hsp.target.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(1550239))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(
            feature.qualifiers["coded_by"],
            "gi|1336928286|emb|LT900055.1|:325704..326405",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQNAYPKDCGAFTGEITSKHLEELKINTLLIGHSERR LLKESPNFLKEKFDFFKDKKFKIVYCIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHGFLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGS SLELENFKTIISFL",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|133692         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQN
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
CAJ99216.         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQN

gi|133692        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRALLKESPNFLKEKFDFFKDKKFKIVY
                 60 ||||||||||||||||||||||||||||||||||.|||||||||||||||||||||||||
CAJ99216.        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPNFLKEKFDFFKDKKFKIVY

gi|133692       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHG
                120 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
CAJ99216.       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHG

gi|133692       180 FLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSASLELENFKTIISFL 234
                180 |||||||||||||||||||||||||||||||||||||||.|||||||||||||| 234
CAJ99216.       180 FLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL 234
""",
        )
        hit = record[2]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|1033012332|gb|CP011486.1|")
        self.assertEqual(hit.target.name, "CP011486")
        self.assertEqual(
            hit.target.description, "Helicobacter pylori strain K26A1, complete genome"
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=1570310)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 1164.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 452.981)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.45778e-143)
        self.assertEqual(hsp.annotations["identity"], 221)
        self.assertEqual(hsp.annotations["positive"], 224)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0, 234],
                          [  0, 234]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 234))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHF...SFL')",
        )
        self.assertEqual(hsp.query.id, "CAJ99216.1")
        self.assertEqual(
            hsp.query.description, "tim [Helicobacter acinonychis str. Sheeba]"
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHF...SFL')",
        )
        self.assertEqual(hsp.target.id, "gi|1033012332|gb|CP011486.1|")
        self.assertEqual(hsp.target.name, "CP011486")
        self.assertEqual(
            hsp.target.description, "Helicobacter pylori strain K26A1, complete genome"
        )
        self.assertEqual(len(hsp.target.features), 1)
        feature = hsp.target.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(1570310))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(
            feature.qualifiers["coded_by"],
            "gi|1033012332|gb|CP011486.1|:196043..196744",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLG QNAYPKDCGAFTGEITSKHLEELKINTLLIGHSERR LLKESP+FLKEKFDFFKDK FKI+YCIGEDLKTREKGL AVKEFLNEQLENIDL Y NLIVAYEPIWAIGT KSASLEDIYLTHGFLKQ LNQK PLLYGGSVNTQNAKEILGIDSVDGLL+GS SLELENFKTIISFL",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|103301         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGAQN
                  0 |||||||||||||||||||||||||||||||||||||||||||||||||||||||||.||
CAJ99216.         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQN

gi|103301        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRTLLKESPSFLKEKFDFFKDKNFKIIY
                 60 ||||||||||||||||||||||||||||||||||.||||||.||||||||||||.|||.|
CAJ99216.        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPNFLKEKFDFFKDKKFKIVY

gi|103301       120 CIGEDLKTREKGLAAVKEFLNEQLENIDLSYHNLIVAYEPIWAIGTKKSASLEDIYLTHG
                120 |||||||||||||.|||||||||||||||.|.||||||||||||||.|||||||||||||
CAJ99216.       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHG

gi|103301       180 FLKQILNQKTPLLYGGSVNTQNAKEILGIDSVDGLLVGSASLELENFKTIISFL 234
                180 ||||.||||.||||||||||||||||||||||||||.||.|||||||||||||| 234
CAJ99216.       180 FLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL 234
""",
        )
        hit = record[3]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|1559948212|dbj|AP017633.1|")
        self.assertEqual(hit.target.name, "AP017633")
        self.assertEqual(
            hit.target.description,
            "Helicobacter pylori DNA, complete genome, strain: PMSS1",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=1603093)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 1136.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 442.195)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.08706e-139)
        self.assertEqual(hsp.annotations["identity"], 218)
        self.assertEqual(hsp.annotations["positive"], 223)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0, 234],
                          [  0, 234]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 234))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHF...SFL')",
        )
        self.assertEqual(hsp.query.id, "CAJ99216.1")
        self.assertEqual(
            hsp.query.description, "tim [Helicobacter acinonychis str. Sheeba]"
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHFDRVFVFPDFLGLLPNSFLHF...SFL')",
        )
        self.assertEqual(hsp.target.id, "gi|1559948212|dbj|AP017633.1|")
        self.assertEqual(hsp.target.name, "AP017633")
        self.assertEqual(
            hsp.target.description,
            "Helicobacter pylori DNA, complete genome, strain: PMSS1",
        )
        self.assertEqual(len(hsp.target.features), 1)
        feature = hsp.target.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(1603093))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(
            feature.qualifiers["coded_by"],
            "gi|1559948212|dbj|AP017633.1|:1330168..1330869",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQH DRVFVFPDFLGLLPN+FLHFTLGVQNAYP+DCGAFTGEITSKHLEELKI+TLLIGHSERRVLLKESP+FLKEKFDFFKDK FKIVYCIGEDL TREKG  AVKEFLNEQLENIDL+Y NLIVAYEPIWAIGT KSASLEDIYLTHGFLKQ LNQK PLLYGGSVNTQNAKEILGIDSVDGLLIGS S ELENFKTIISFL",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|155994         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHFDRVFVFPDFLGLLPNSFLHFTLGVQN
                  0 |||||||||||||||||||||||||||||||||.|||||||||||||||.||||||||||
CAJ99216.         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQN

gi|155994        60 AYPRDCGAFTGEITSKHLEELKIHTLLIGHSERRVLLKESPSFLKEKFDFFKDKNFKIVY
                 60 |||.|||||||||||||||||||.|||||||||||||||||.||||||||||||.|||||
CAJ99216.        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPNFLKEKFDFFKDKKFKIVY

gi|155994       120 CIGEDLTTREKGFKAVKEFLNEQLENIDLNYSNLIVAYEPIWAIGTKKSASLEDIYLTHG
                120 ||||||.|||||..|||||||||||||||.|.||||||||||||||.|||||||||||||
CAJ99216.       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHG

gi|155994       180 FLKQILNQKTPLLYGGSVNTQNAKEILGIDSVDGLLIGSASWELENFKTIISFL 234
                180 ||||.||||.|||||||||||||||||||||||||||||.|.|||||||||||| 234
CAJ99216.       180 FLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL 234
""",
        )
        hit = record[4]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|1143706535|gb|CP018823.1|")
        self.assertEqual(hit.target.name, "CP018823")
        self.assertEqual(
            hit.target.description, "Helicobacter pylori strain PMSS1 complete genome"
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=1618480)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 1136.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 442.195)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.08707e-139)
        self.assertEqual(hsp.annotations["identity"], 218)
        self.assertEqual(hsp.annotations["positive"], 223)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0, 234],
                          [  0, 234]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 234))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHF...SFL')",
        )
        self.assertEqual(hsp.query.id, "CAJ99216.1")
        self.assertEqual(
            hsp.query.description, "tim [Helicobacter acinonychis str. Sheeba]"
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHFDRVFVFPDFLGLLPNSFLHF...SFL')",
        )
        self.assertEqual(hsp.target.id, "gi|1143706535|gb|CP018823.1|")
        self.assertEqual(hsp.target.name, "CP018823")
        self.assertEqual(
            hsp.target.description, "Helicobacter pylori strain PMSS1 complete genome"
        )
        self.assertEqual(len(hsp.target.features), 1)
        feature = hsp.target.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(1618480))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(
            feature.qualifiers["coded_by"],
            "gi|1143706535|gb|CP018823.1|:190464..191165",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQH DRVFVFPDFLGLLPN+FLHFTLGVQNAYP+DCGAFTGEITSKHLEELKI+TLLIGHSERRVLLKESP+FLKEKFDFFKDK FKIVYCIGEDL TREKG  AVKEFLNEQLENIDL+Y NLIVAYEPIWAIGT KSASLEDIYLTHGFLKQ LNQK PLLYGGSVNTQNAKEILGIDSVDGLLIGS S ELENFKTIISFL",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|114370         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHFDRVFVFPDFLGLLPNSFLHFTLGVQN
                  0 |||||||||||||||||||||||||||||||||.|||||||||||||||.||||||||||
CAJ99216.         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQN

gi|114370        60 AYPRDCGAFTGEITSKHLEELKIHTLLIGHSERRVLLKESPSFLKEKFDFFKDKNFKIVY
                 60 |||.|||||||||||||||||||.|||||||||||||||||.||||||||||||.|||||
CAJ99216.        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPNFLKEKFDFFKDKKFKIVY

gi|114370       120 CIGEDLTTREKGFKAVKEFLNEQLENIDLNYSNLIVAYEPIWAIGTKKSASLEDIYLTHG
                120 ||||||.|||||..|||||||||||||||.|.||||||||||||||.|||||||||||||
CAJ99216.       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHG

gi|114370       180 FLKQILNQKTPLLYGGSVNTQNAKEILGIDSVDGLLIGSASWELENFKTIISFL 234
                180 ||||.||||.|||||||||||||||||||||||||||||.|.|||||||||||| 234
CAJ99216.       180 FLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL 234
""",
        )
        hit = record[5]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|1149039824|gb|CP009259.1|")
        self.assertEqual(hit.target.name, "CP009259")
        self.assertEqual(
            hit.target.description, "Helicobacter pylori SS1, complete genome"
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=1619098)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 1136.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 442.195)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.08707e-139)
        self.assertEqual(hsp.annotations["identity"], 218)
        self.assertEqual(hsp.annotations["positive"], 223)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0, 234],
                          [  0, 234]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 234))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHF...SFL')",
        )
        self.assertEqual(hsp.query.id, "CAJ99216.1")
        self.assertEqual(
            hsp.query.description, "tim [Helicobacter acinonychis str. Sheeba]"
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHFDRVFVFPDFLGLLPNSFLHF...SFL')",
        )
        self.assertEqual(hsp.target.id, "gi|1149039824|gb|CP009259.1|")
        self.assertEqual(hsp.target.name, "CP009259")
        self.assertEqual(
            hsp.target.description, "Helicobacter pylori SS1, complete genome"
        )
        self.assertEqual(len(hsp.target.features), 1)
        feature = hsp.target.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(1619098))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(
            feature.qualifiers["coded_by"],
            "gi|1149039824|gb|CP009259.1|:190465..191166",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQH DRVFVFPDFLGLLPN+FLHFTLGVQNAYP+DCGAFTGEITSKHLEELKI+TLLIGHSERRVLLKESP+FLKEKFDFFKDK FKIVYCIGEDL TREKG  AVKEFLNEQLENIDL+Y NLIVAYEPIWAIGT KSASLEDIYLTHGFLKQ LNQK PLLYGGSVNTQNAKEILGIDSVDGLLIGS S ELENFKTIISFL",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|114903         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHFDRVFVFPDFLGLLPNSFLHFTLGVQN
                  0 |||||||||||||||||||||||||||||||||.|||||||||||||||.||||||||||
CAJ99216.         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQN

gi|114903        60 AYPRDCGAFTGEITSKHLEELKIHTLLIGHSERRVLLKESPSFLKEKFDFFKDKNFKIVY
                 60 |||.|||||||||||||||||||.|||||||||||||||||.||||||||||||.|||||
CAJ99216.        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPNFLKEKFDFFKDKKFKIVY

gi|114903       120 CIGEDLTTREKGFKAVKEFLNEQLENIDLNYSNLIVAYEPIWAIGTKKSASLEDIYLTHG
                120 ||||||.|||||..|||||||||||||||.|.||||||||||||||.|||||||||||||
CAJ99216.       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHG

gi|114903       180 FLKQILNQKTPLLYGGSVNTQNAKEILGIDSVDGLLIGSASWELENFKTIISFL 234
                180 ||||.||||.|||||||||||||||||||||||||||||.|.|||||||||||| 234
CAJ99216.       180 FLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL 234
""",
        )
        hit = record[6]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|1560269298|gb|CP035106.1|")
        self.assertEqual(hit.target.name, "CP035106")
        self.assertEqual(
            hit.target.description,
            "Helicobacter pylori strain Hpbs3 chromosome, complete genome",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=1514674)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 1133.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 441.039)
        self.assertAlmostEqual(hsp.annotations["evalue"], 5.81854e-139)
        self.assertEqual(hsp.annotations["identity"], 215)
        self.assertEqual(hsp.annotations["positive"], 224)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0, 234],
                          [  0, 234]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 234))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHF...SFL')",
        )
        self.assertEqual(hsp.query.id, "CAJ99216.1")
        self.assertEqual(
            hsp.query.description, "tim [Helicobacter acinonychis str. Sheeba]"
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MTKIAMANFKSAMPIFKSHAYLEELEKTLKPQHSDRVFVFPDFLGLLPNSFLHF...SFL')",
        )
        self.assertEqual(hsp.target.id, "gi|1560269298|gb|CP035106.1|")
        self.assertEqual(hsp.target.name, "CP035106")
        self.assertEqual(
            hsp.target.description,
            "Helicobacter pylori strain Hpbs3 chromosome, complete genome",
        )
        self.assertEqual(len(hsp.target.features), 1)
        feature = hsp.target.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(1514674))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(
            feature.qualifiers["coded_by"],
            "gi|1560269298|gb|CP035106.1|:1231590..1232291",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MTKIAMANFKSAMPIFKSHAYL+ELEKTLKPQH DRVFVFPDFLGLLPN+FLHFTLGVQNAYP+DCGAFTGEITS+HLEELKINTLLIGHSERR+LLKESP+FLKEKFDFFK K FKIVYCIGE+L TREKG  AVKEFLNEQLENIDL+Y NL+VAYEPIWAIGT KSASLEDIYLTHGFLKQ LNQK PLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|156026         0 MTKIAMANFKSAMPIFKSHAYLEELEKTLKPQHSDRVFVFPDFLGLLPNSFLHFTLGVQN
                  0 ||||||||||||||||||||||.||||||||||.|||||||||||||||.||||||||||
CAJ99216.         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQN

gi|156026        60 AYPRDCGAFTGEITSQHLEELKINTLLIGHSERRLLLKESPSFLKEKFDFFKSKNFKIVY
                 60 |||.|||||||||||.||||||||||||||||||.||||||.||||||||||.|.|||||
CAJ99216.        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPNFLKEKFDFFKDKKFKIVY

gi|156026       120 CIGEELTTREKGFKAVKEFLNEQLENIDLNYPNLVVAYEPIWAIGTKKSASLEDIYLTHG
                120 ||||.|.|||||..|||||||||||||||.|.||.|||||||||||.|||||||||||||
CAJ99216.       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHG

gi|156026       180 FLKQVLNQKTPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL 234
                180 ||||.||||.|||||||||||||||||||||||||||||||||||||||||||| 234
CAJ99216.       180 FLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL 234
""",
        )
        hit = record[7]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|1540258434|emb|LR134517.1|")
        self.assertEqual(hit.target.name, "LR134517")
        self.assertEqual(
            hit.target.description,
            "Helicobacter pylori strain NCTC13345 genome assembly, chromosome: 1",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=1644315)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 1131.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 440.269)
        self.assertAlmostEqual(hsp.annotations["evalue"], 9.62198e-139)
        self.assertEqual(hsp.annotations["identity"], 216)
        self.assertEqual(hsp.annotations["positive"], 222)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0, 234],
                          [  0, 234]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 234))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHF...SFL')",
        )
        self.assertEqual(hsp.query.id, "CAJ99216.1")
        self.assertEqual(
            hsp.query.description, "tim [Helicobacter acinonychis str. Sheeba]"
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHFDRVFVFPDFLGLLPNAFLHF...SFL')",
        )
        self.assertEqual(hsp.target.id, "gi|1540258434|emb|LR134517.1|")
        self.assertEqual(hsp.target.name, "LR134517")
        self.assertEqual(
            hsp.target.description,
            "Helicobacter pylori strain NCTC13345 genome assembly, chromosome: 1",
        )
        self.assertEqual(len(hsp.target.features), 1)
        feature = hsp.target.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(1644315))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(
            feature.qualifiers["coded_by"],
            "gi|1540258434|emb|LR134517.1|:363062..363763",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQH DRVFVFPDFLGLLPNAFLHFTLG QNAYP+DCGAFTGEITSKHLEELKI+TLLIGHSERR LLKESP+FLKEKFDFFKDK FKIVYC+GEDL TREKG  AVKEFL+EQLENIDL+Y NLIVAYEPIWAIGT KSASLEDIYLTHGFLKQ LNQK PLLYGGSVNTQNAKEILGIDSVDGLLIGS SLELENFKTIISFL",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|154025         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHFDRVFVFPDFLGLLPNAFLHFTLGAQN
                  0 |||||||||||||||||||||||||||||||||.|||||||||||||||||||||||.||
CAJ99216.         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQN

gi|154025        60 AYPRDCGAFTGEITSKHLEELKIHTLLIGHSERRALLKESPSFLKEKFDFFKDKNFKIVY
                 60 |||.|||||||||||||||||||.||||||||||.||||||.||||||||||||.|||||
CAJ99216.        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPNFLKEKFDFFKDKKFKIVY

gi|154025       120 CVGEDLTTREKGFRAVKEFLSEQLENIDLNYSNLIVAYEPIWAIGTKKSASLEDIYLTHG
                120 |.||||.|||||..||||||.||||||||.|.||||||||||||||.|||||||||||||
CAJ99216.       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHG

gi|154025       180 FLKQILNQKTPLLYGGSVNTQNAKEILGIDSVDGLLIGSASLELENFKTIISFL 234
                180 ||||.||||.|||||||||||||||||||||||||||||.|||||||||||||| 234
CAJ99216.       180 FLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL 234
""",
        )
        hit = record[8]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|1033005233|gb|CP011487.1|")
        self.assertEqual(hit.target.name, "CP011487")
        self.assertEqual(
            hit.target.description, "Helicobacter pylori strain PNG84A, complete genome"
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=1531450)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 1129.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 439.499)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.21067e-138)
        self.assertEqual(hsp.annotations["identity"], 215)
        self.assertEqual(hsp.annotations["positive"], 221)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0, 234],
                          [  0, 234]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 234))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHF...SFL')",
        )
        self.assertEqual(hsp.query.id, "CAJ99216.1")
        self.assertEqual(
            hsp.query.description, "tim [Helicobacter acinonychis str. Sheeba]"
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHFDRVFVFPDFLGLLPNAFLHF...SFL')",
        )
        self.assertEqual(hsp.target.id, "gi|1033005233|gb|CP011487.1|")
        self.assertEqual(hsp.target.name, "CP011487")
        self.assertEqual(
            hsp.target.description, "Helicobacter pylori strain PNG84A, complete genome"
        )
        self.assertEqual(len(hsp.target.features), 1)
        feature = hsp.target.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(1531450))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(
            feature.qualifiers["coded_by"],
            "gi|1033005233|gb|CP011487.1|:161806..162507",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQH DRVFVFPDFLGLLPNAFLHFTLGVQNAYP+DCGAFTGEITSKHLEELKINTLLIGHSERR+LLKESP+FLKEKFDFFK K FKI+YCIGE+L TREK   AVKEFLNEQLENIDLDY NL+VAYEPIWAIG  KSASLEDIYLTHGFLKQ LNQK PLLYGGSVNTQNAKEILGIDSVDGLLIGS SLELENFKTIISFL",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|103300         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHFDRVFVFPDFLGLLPNAFLHFTLGVQN
                  0 |||||||||||||||||||||||||||||||||.||||||||||||||||||||||||||
CAJ99216.         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQN

gi|103300        60 AYPRDCGAFTGEITSKHLEELKINTLLIGHSERRMLLKESPSFLKEKFDFFKSKNFKIIY
                 60 |||.||||||||||||||||||||||||||||||.||||||.||||||||||.|.|||.|
CAJ99216.        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPNFLKEKFDFFKDKKFKIVY

gi|103300       120 CIGEELTTREKSFKAVKEFLNEQLENIDLDYPNLVVAYEPIWAIGAKKSASLEDIYLTHG
                120 ||||.|.||||...|||||||||||||||||.||.||||||||||..|||||||||||||
CAJ99216.       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHG

gi|103300       180 FLKQILNQKTPLLYGGSVNTQNAKEILGIDSVDGLLIGSASLELENFKTIISFL 234
                180 ||||.||||.|||||||||||||||||||||||||||||.|||||||||||||| 234
CAJ99216.       180 FLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL 234
""",
        )
        hit = record[9]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|1509197163|gb|CP025748.1|")
        self.assertEqual(hit.target.name, "CP025748")
        self.assertEqual(
            hit.target.description, "Helicobacter pylori strain Hp_TH2099 chromosome"
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=1667396)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 1128.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 439.113)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.2539e-138)
        self.assertEqual(hsp.annotations["identity"], 216)
        self.assertEqual(hsp.annotations["positive"], 220)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0, 234],
                          [  0, 234]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 234))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHF...SFL')",
        )
        self.assertEqual(hsp.query.id, "CAJ99216.1")
        self.assertEqual(
            hsp.query.description, "tim [Helicobacter acinonychis str. Sheeba]"
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHFDRVFVFPDFLGLLPNSFLHF...SFL')",
        )
        self.assertEqual(hsp.target.id, "gi|1509197163|gb|CP025748.1|")
        self.assertEqual(hsp.target.name, "CP025748")
        self.assertEqual(
            hsp.target.description, "Helicobacter pylori strain Hp_TH2099 chromosome"
        )
        self.assertEqual(len(hsp.target.features), 1)
        feature = hsp.target.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(1667396))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(
            feature.qualifiers["coded_by"],
            "gi|1509197163|gb|CP025748.1|:200718..201419",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQH DRVFVFPDFLGLLPN+FLHFTLG QNAYPKDCGAFTGEITSKHLEELKINTLLIGHSERR LLKESPNFLKEKFDFFK K FKIVYCIGE+L TREKG  AVKEFLNEQLENIDL+Y NL+VAYEPIWAIGT KSASLEDIYLTHGFLKQ LNQK PLLYGGSVN QNAKEILGIDSVDGLLIGS SLELENFKTIISFL",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|150919         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHFDRVFVFPDFLGLLPNSFLHFTLGAQN
                  0 |||||||||||||||||||||||||||||||||.|||||||||||||||.|||||||.||
CAJ99216.         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQN

gi|150919        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRALLKESPNFLKEKFDFFKSKNFKIVY
                 60 ||||||||||||||||||||||||||||||||||.|||||||||||||||||.|.|||||
CAJ99216.        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPNFLKEKFDFFKDKKFKIVY

gi|150919       120 CIGEELTTREKGFKAVKEFLNEQLENIDLNYPNLVVAYEPIWAIGTKKSASLEDIYLTHG
                120 ||||.|.|||||..|||||||||||||||.|.||.|||||||||||.|||||||||||||
CAJ99216.       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHG

gi|150919       180 FLKQILNQKTPLLYGGSVNVQNAKEILGIDSVDGLLIGSASLELENFKTIISFL 234
                180 ||||.||||.|||||||||.|||||||||||||||||||.|||||||||||||| 234
CAJ99216.       180 FLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL 234
""",
        )


class TestTBlastx(unittest.TestCase):
    """Test the Blast XML parser for tblastx output."""

    def test_xml_2226_tblastx_004(self):
        """Parsing TBLASTX 2.2.26+ (xml_2226_tblastx_004.xml)."""
        filename = "xml_2226_tblastx_004.xml"
        datafile = os.path.join("Blast", filename)
        with open(datafile, "rb") as handle:
            records = Blast.parse(handle)
            self.check_xml_2226_tblastx_004(records)
        with Blast.parse(datafile) as records:
            self.check_xml_2226_tblastx_004(records)

    def check_xml_2226_tblastx_004(self, records):
        self.assertEqual(records.program, "tblastx")
        self.assertEqual(records.version, "TBLASTX 2.2.26+")
        self.assertEqual(
            records.reference,
            'Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.',
        )
        self.assertEqual(records.db, "refseq_rna")
        self.assertIsInstance(records.query, SeqRecord)
        self.assertEqual(records.query.id, "Query_1")
        self.assertEqual(records.query.description, "random_s00")
        self.assertEqual(repr(records.query.seq), "Seq(None, length=128)")
        self.assertEqual(len(records.param), 5)
        self.assertEqual(records.param["matrix"], "BLOSUM62")
        self.assertAlmostEqual(records.param["expect"], 10.0)
        self.assertEqual(records.param["gap-open"], 11)
        self.assertEqual(records.param["gap-extend"], 1)
        self.assertEqual(records.param["filter"], "L;")
        record = next(records)
        self.assertIsInstance(record.query, SeqRecord)
        self.assertEqual(record.query.id, "Query_1")
        self.assertEqual(record.query.description, "random_s00")
        self.assertEqual(repr(record.query.seq), "Seq(None, length=128)")

        self.assertEqual(len(record.stat), 7)
        self.assertEqual(record.stat["db-num"], 2933984)
        self.assertEqual(record.stat["db-len"], 4726730735)
        self.assertEqual(record.stat["hsp-len"], 0)
        self.assertAlmostEqual(record.stat["eff-space"], 0.0)
        self.assertAlmostEqual(record.stat["kappa"], 0.0)
        self.assertAlmostEqual(record.stat["lambda"], 0.0)
        self.assertAlmostEqual(record.stat["entropy"], 0.0)
        self.assertEqual(len(record), 0)
        record = next(records)
        self.assertIsInstance(record.query, SeqRecord)
        self.assertEqual(record.query.id, "Query_2")
        self.assertEqual(
            record.query.description,
            "gi|296147483:1-350 Saccharomyces cerevisiae S288c Mon2p (MON2) mRNA, complete cds",
        )
        self.assertEqual(repr(record.query.seq), "Seq(None, length=350)")

        self.assertEqual(len(record.stat), 7)
        self.assertEqual(record.stat["db-num"], 2933984)
        self.assertEqual(record.stat["db-len"], 4726730735)
        self.assertEqual(record.stat["hsp-len"], 0)
        self.assertAlmostEqual(record.stat["eff-space"], 0.0)
        self.assertAlmostEqual(record.stat["kappa"], 0.0)
        self.assertAlmostEqual(record.stat["lambda"], 0.0)
        self.assertAlmostEqual(record.stat["entropy"], 0.0)
        self.assertEqual(len(record), 5)
        hit = record[0]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|296147483|ref|NM_001183135.1|")
        self.assertEqual(hit.target.name, "NM_001183135")
        self.assertEqual(
            hit.target.description,
            "Saccharomyces cerevisiae S288c Mon2p (MON2) mRNA, complete cds >gi|116616412|gb|EF059095.1| Synthetic construct Saccharomyces cerevisiae clone FLH203015.01X MON2, complete sequence",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=4911)")
        self.assertEqual(len(hit), 7)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 626.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 289.739)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.05874e-76)
        self.assertEqual(hsp.annotations["identity"], 116)
        self.assertEqual(hsp.annotations["positive"], 116)
        self.assertEqual(hsp.annotations["gaps"], 0)
        hsp = hit[1]
        self.assertAlmostEqual(hsp.score, 602.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 278.742)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.16381e-73)
        self.assertEqual(hsp.annotations["identity"], 116)
        self.assertEqual(hsp.annotations["positive"], 116)
        self.assertEqual(hsp.annotations["gaps"], 0)
        hsp = hit[2]
        self.assertAlmostEqual(hsp.score, 593.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 274.618)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.77251e-72)
        self.assertEqual(hsp.annotations["identity"], 116)
        self.assertEqual(hsp.annotations["positive"], 116)
        self.assertEqual(hsp.annotations["gaps"], 0)
        hsp = hit[3]
        self.assertAlmostEqual(hsp.score, 583.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 270.036)
        self.assertAlmostEqual(hsp.annotations["evalue"], 9.03598e-71)
        self.assertEqual(hsp.annotations["identity"], 116)
        self.assertEqual(hsp.annotations["positive"], 116)
        self.assertEqual(hsp.annotations["gaps"], 0)
        hsp = hit[4]
        self.assertAlmostEqual(hsp.score, 495.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 229.713)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.24226e-58)
        self.assertEqual(hsp.annotations["identity"], 116)
        self.assertEqual(hsp.annotations["positive"], 116)
        self.assertEqual(hsp.annotations["gaps"], 0)
        hsp = hit[5]
        self.assertAlmostEqual(hsp.score, 425.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 197.639)
        self.assertAlmostEqual(hsp.annotations["evalue"], 9.12288e-54)
        self.assertEqual(hsp.annotations["identity"], 85)
        self.assertEqual(hsp.annotations["positive"], 85)
        self.assertEqual(hsp.annotations["gaps"], 0)
        hsp = hit[6]
        self.assertAlmostEqual(hsp.score, 73.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 36.3494)
        self.assertAlmostEqual(hsp.annotations["evalue"], 9.12288e-54)
        self.assertEqual(hsp.annotations["identity"], 14)
        self.assertEqual(hsp.annotations["positive"], 14)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 14],
                              [ 0, 14]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 14))
        self.assertEqual(repr(hsp.query.seq), "Seq('MAMNTGGFDSMQRQ')")
        self.assertEqual(hsp.query.id, "Query_2")
        self.assertEqual(
            hsp.query.description,
            "gi|296147483:1-350 Saccharomyces cerevisiae S288c Mon2p (MON2) mRNA, complete cds",
        )
        self.assertEqual(len(hsp.query.features), 1)
        feature = hsp.query.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(14))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(feature.qualifiers["coded_by"], "Query_2:1..42")
        self.assertEqual(repr(hsp.target.seq), "Seq('MAMNTGGFDSMQRQ')")
        self.assertEqual(hsp.target.id, "gi|296147483|ref|NM_001183135.1|")
        self.assertEqual(hsp.target.name, "NM_001183135")
        self.assertEqual(
            hsp.target.description,
            "Saccharomyces cerevisiae S288c Mon2p (MON2) mRNA, complete cds >gi|116616412|gb|EF059095.1| Synthetic construct Saccharomyces cerevisiae clone FLH203015.01X MON2, complete sequence",
        )
        self.assertEqual(len(hsp.target.features), 1)
        feature = hsp.target.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(4911))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(
            feature.qualifiers["coded_by"], "gi|296147483|ref|NM_001183135.1|:1..42"
        )
        self.assertEqual(hsp.annotations["midline"], "MAMNTGGFDSMQRQ")
        self.assertEqual(
            str(hsp),
            """\
gi|296147         0 MAMNTGGFDSMQRQ 14
                  0 |||||||||||||| 14
Query_2           0 MAMNTGGFDSMQRQ 14
""",
        )
        hit = record[1]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|365982352|ref|XM_003667962.1|")
        self.assertEqual(hit.target.name, "XM_003667962")
        self.assertEqual(
            hit.target.description,
            "Naumovozyma dairenensis CBS 421 hypothetical protein (NDAI0A06120), mRNA",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=4932)")
        self.assertEqual(len(hit), 6)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 327.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 152.734)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.38069e-37)
        self.assertEqual(hsp.annotations["identity"], 62)
        self.assertEqual(hsp.annotations["positive"], 73)
        self.assertEqual(hsp.annotations["gaps"], 0)
        hsp = hit[1]
        self.assertAlmostEqual(hsp.score, 51.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 26.2688)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.38069e-37)
        self.assertEqual(hsp.annotations["identity"], 11)
        self.assertEqual(hsp.annotations["positive"], 11)
        self.assertEqual(hsp.annotations["gaps"], 0)
        hsp = hit[2]
        self.assertAlmostEqual(hsp.score, 142.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 67.9658)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.80116e-20)
        self.assertEqual(hsp.annotations["identity"], 34)
        self.assertEqual(hsp.annotations["positive"], 38)
        self.assertEqual(hsp.annotations["gaps"], 0)
        hsp = hit[3]
        self.assertAlmostEqual(hsp.score, 109.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 52.8449)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.80116e-20)
        self.assertEqual(hsp.annotations["identity"], 24)
        self.assertEqual(hsp.annotations["positive"], 29)
        self.assertEqual(hsp.annotations["gaps"], 0)
        hsp = hit[4]
        self.assertAlmostEqual(hsp.score, 127.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 61.0927)
        self.assertAlmostEqual(hsp.annotations["evalue"], 7.14684e-08)
        self.assertEqual(hsp.annotations["identity"], 36)
        self.assertEqual(hsp.annotations["positive"], 52)
        self.assertEqual(hsp.annotations["gaps"], 0)
        hsp = hit[5]
        self.assertAlmostEqual(hsp.score, 87.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 42.7643)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.0235231)
        self.assertEqual(hsp.annotations["identity"], 28)
        self.assertEqual(hsp.annotations["positive"], 36)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 65],
                              [ 0, 65]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 65))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('G*QSL*ALHCQGRHFSIP*LASQHERECEIRMSF*LLKTMYSFQYLNGFITSMA...FGR')",
        )
        self.assertEqual(hsp.query.id, "Query_2")
        self.assertEqual(
            hsp.query.description,
            "gi|296147483:1-350 Saccharomyces cerevisiae S288c Mon2p (MON2) mRNA, complete cds",
        )
        self.assertEqual(len(hsp.query.features), 1)
        feature = hsp.query.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(65))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(feature.qualifiers["coded_by"], "complement(Query_2:67..261)")
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('GRQSL*ALHGN*S*FGISR*TGQDQRCNKIRVPYKFFYILNSL*NID*FVASMF...F*R')",
        )
        self.assertEqual(hsp.target.id, "gi|365982352|ref|XM_003667962.1|")
        self.assertEqual(hsp.target.name, "XM_003667962")
        self.assertEqual(
            hsp.target.description,
            "Naumovozyma dairenensis CBS 421 hypothetical protein (NDAI0A06120), mRNA",
        )
        self.assertEqual(len(hsp.target.features), 1)
        feature = hsp.target.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(4932))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(
            feature.qualifiers["coded_by"],
            "complement(gi|365982352|ref|XM_003667962.1|:61..255)",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "G QSL*ALH     F I     Q +R  +IR+ +     + S   ++ F+ SM NG*ISSFRF R",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|365982         0 GRQSL*ALHGN*S*FGISR*TGQDQRCNKIRVPYKFFYILNSL*NID*FVASMFNG*ISS
                  0 |.|||||||.....|.|.....|..|...||..........|......|..||.||||||
Query_2           0 G*QSL*ALHCQGRHFSIP*LASQHERECEIRMSF*LLKTMYSFQYLNGFITSMANG*ISS

gi|365982        60 FRF*R 65
                 60 |||.| 65
Query_2          60 FRFGR 65
""",
        )
        hit = record[2]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|366988334|ref|XM_003673886.1|")
        self.assertEqual(hit.target.name, "XM_003673886")
        self.assertEqual(
            hit.target.description,
            "Naumovozyma castellii CBS 4309 hypothetical protein (NCAS0A09950) mRNA, complete cds",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=4938)")
        self.assertEqual(len(hit), 4)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 306.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 143.112)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.45826e-32)
        self.assertEqual(hsp.annotations["identity"], 58)
        self.assertEqual(hsp.annotations["positive"], 71)
        self.assertEqual(hsp.annotations["gaps"], 0)
        hsp = hit[1]
        self.assertAlmostEqual(hsp.score, 130.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 62.4673)
        self.assertAlmostEqual(hsp.annotations["evalue"], 5.61057e-16)
        self.assertEqual(hsp.annotations["identity"], 30)
        self.assertEqual(hsp.annotations["positive"], 36)
        self.assertEqual(hsp.annotations["gaps"], 0)
        hsp = hit[2]
        self.assertAlmostEqual(hsp.score, 91.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 44.5971)
        self.assertAlmostEqual(hsp.annotations["evalue"], 5.61057e-16)
        self.assertEqual(hsp.annotations["identity"], 20)
        self.assertEqual(hsp.annotations["positive"], 24)
        self.assertEqual(hsp.annotations["gaps"], 0)
        hsp = hit[3]
        self.assertAlmostEqual(hsp.score, 112.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 54.2195)
        self.assertAlmostEqual(hsp.annotations["evalue"], 8.37784e-06)
        self.assertEqual(hsp.annotations["identity"], 38)
        self.assertEqual(hsp.annotations["positive"], 58)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0, 115],
                              [  0, 115]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 115))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('*LNLHREMSSLNEGIQNFRQPASRNRWNG*QSL*ALHCQGRHFSIP*LASQHER...HGH')",
        )
        self.assertEqual(hsp.query.id, "Query_2")
        self.assertEqual(
            hsp.query.description,
            "gi|296147483:1-350 Saccharomyces cerevisiae S288c Mon2p (MON2) mRNA, complete cds",
        )
        self.assertEqual(len(hsp.query.features), 1)
        feature = hsp.query.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(115))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(feature.qualifiers["coded_by"], "complement(Query_2:1..345)")
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('*LYFNCKMCRFYKSIENSR*SAFRYGWCRC*AL*TLHSNRGQLSVP**T*FYQW...NSH')",
        )
        self.assertEqual(hsp.target.id, "gi|366988334|ref|XM_003673886.1|")
        self.assertEqual(hsp.target.name, "XM_003673886")
        self.assertEqual(
            hsp.target.description,
            "Naumovozyma castellii CBS 4309 hypothetical protein (NCAS0A09950) mRNA, complete cds",
        )
        self.assertEqual(len(hsp.target.features), 1)
        feature = hsp.target.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(4938))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(
            feature.qualifiers["coded_by"],
            "complement(gi|366988334|ref|XM_003673886.1|:1..345)",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "*L  + +M    + I+N R  A R  W    +L* LH      S+P*    ++   +IRM+  +   + SF  L+  +T M +  + SFRFGR   QF  +L L G+K  S++ H",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|366988         0 *LYFNCKMCRFYKSIENSR*SAFRYGWCRC*AL*TLHSNRGQLSVP**T*FYQWSDKIRM
                  0 ||.....|......|.|.|..|.|..|.....||.||......|.||..........|||
Query_2           0 *LNLHREMSSLNEGIQNFRQPASRNRWNG*QSL*ALHCQGRHFSIP*LASQHERECEIRM

gi|366988        60 ACQIFNVLDSF*DLDRLVTCMLDRSVPSFRFGRQ*MQF*VQLFLKGLKTCSLNSH 115
                 60 .........||..|....|.|......||||||...||...|.|.|.|..|...| 115
Query_2          60 SF*LLKTMYSFQYLNGFITSMANG*ISSFRFGR*RTQFCFKLPLHGVKPSSVHGH 115
""",
        )
        hit = record[3]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|255710474|ref|XM_002551475.1|")
        self.assertEqual(hit.target.name, "XM_002551475")
        self.assertEqual(
            hit.target.description,
            "Lachancea thermotolerans CBS 6340 KLTH0A01342p (KLTH0A01342g) mRNA, complete cds",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=4845)")
        self.assertEqual(len(hit), 2)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 303.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 141.737)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.78129e-32)
        self.assertEqual(hsp.annotations["identity"], 55)
        self.assertEqual(hsp.annotations["positive"], 71)
        self.assertEqual(hsp.annotations["gaps"], 0)
        hsp = hit[1]
        self.assertAlmostEqual(hsp.score, 92.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 45.0554)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.00480643)
        self.assertEqual(hsp.annotations["identity"], 25)
        self.assertEqual(hsp.annotations["positive"], 29)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 55],
                              [ 0, 55]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 55))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('TFN*ISIAR*VASMKASKISDSRLRGIDGTVDSPCRHCIARVVILAFLDWQANTK')",
        )
        self.assertEqual(hsp.query.id, "Query_2")
        self.assertEqual(
            hsp.query.description,
            "gi|296147483:1-350 Saccharomyces cerevisiae S288c Mon2p (MON2) mRNA, complete cds",
        )
        self.assertEqual(len(hsp.query.features), 1)
        feature = hsp.query.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(55))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(feature.qualifiers["coded_by"], "complement(Query_2:186..350)")
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('TLSCISAARWVESMNASSTSSILSSGMQLTDDIFCRHCTETLVSLAFLEAHERTK')",
        )
        self.assertEqual(hsp.target.id, "gi|255710474|ref|XM_002551475.1|")
        self.assertEqual(hsp.target.name, "XM_002551475")
        self.assertEqual(
            hsp.target.description,
            "Lachancea thermotolerans CBS 6340 KLTH0A01342p (KLTH0A01342g) mRNA, complete cds",
        )
        self.assertEqual(len(hsp.target.features), 1)
        feature = hsp.target.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(4845))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(
            feature.qualifiers["coded_by"],
            "complement(gi|255710474|ref|XM_002551475.1|:183..347)",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "T + IS AR V SM AS  S     G+  T D  CRHC   +V LAFL+    TK",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|255710         0 TLSCISAARWVESMNASSTSSILSSGMQLTDDIFCRHCTETLVSLAFLEAHERTK 55
                  0 |...||.||.|.||.||..|.....|...|.|..||||....|.||||.....|| 55
Query_2           0 TFN*ISIAR*VASMKASKISDSRLRGIDGTVDSPCRHCIARVVILAFLDWQANTK 55
""",
        )
        hit = record[4]
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|254579534|ref|XM_002495708.1|")
        self.assertEqual(hit.target.name, "XM_002495708")
        self.assertEqual(
            hit.target.description,
            "Zygosaccharomyces rouxii hypothetical protein (ZYRO0C02266g) mRNA, complete cds",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=4866)")
        self.assertEqual(len(hit), 4)
        hsp = hit[0]
        self.assertAlmostEqual(hsp.score, 302.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 141.279)
        self.assertAlmostEqual(hsp.annotations["evalue"], 5.19486e-32)
        self.assertEqual(hsp.annotations["identity"], 57)
        self.assertEqual(hsp.annotations["positive"], 72)
        self.assertEqual(hsp.annotations["gaps"], 0)
        hsp = hit[1]
        self.assertAlmostEqual(hsp.score, 105.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 51.0121)
        self.assertAlmostEqual(hsp.annotations["evalue"], 8.66978e-12)
        self.assertEqual(hsp.annotations["identity"], 27)
        self.assertEqual(hsp.annotations["positive"], 33)
        self.assertEqual(hsp.annotations["gaps"], 0)
        hsp = hit[2]
        self.assertAlmostEqual(hsp.score, 85.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 41.8479)
        self.assertAlmostEqual(hsp.annotations["evalue"], 8.66978e-12)
        self.assertEqual(hsp.annotations["identity"], 20)
        self.assertEqual(hsp.annotations["positive"], 25)
        self.assertEqual(hsp.annotations["gaps"], 0)
        hsp = hit[3]
        self.assertAlmostEqual(hsp.score, 92.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 45.0554)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.00480643)
        self.assertEqual(hsp.annotations["identity"], 31)
        self.assertEqual(hsp.annotations["positive"], 53)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 99],
                              [ 0, 99]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 99))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('RIAFFIFRIEKKKFNHSPC***IH*DIEKST*F*GARKTSGFRTPFRVGLPIKE...SIK')",
        )
        self.assertEqual(hsp.query.id, "Query_2")
        self.assertEqual(
            hsp.query.description,
            "gi|296147483:1-350 Saccharomyces cerevisiae S288c Mon2p (MON2) mRNA, complete cds",
        )
        self.assertEqual(len(hsp.query.features), 1)
        feature = hsp.query.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(99))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(feature.qualifiers["coded_by"], "Query_2:51..347")
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('RVAFTFIRVEEKEYCY*EC*R*IN*NFESSSQL*GIIKTSRFYSTASDVMCIQE...KTK')",
        )
        self.assertEqual(hsp.target.id, "gi|254579534|ref|XM_002495708.1|")
        self.assertEqual(hsp.target.name, "XM_002495708")
        self.assertEqual(
            hsp.target.description,
            "Zygosaccharomyces rouxii hypothetical protein (ZYRO0C02266g) mRNA, complete cds",
        )
        self.assertEqual(len(hsp.target.features), 1)
        feature = hsp.target.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(4866))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(
            feature.qualifiers["coded_by"], "gi|254579534|ref|XM_002495708.1|:51..347"
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "R+AF   R+E+K++ +  C* *I+*+ E S+  *G  KTS F +     + I+EC  D   NAM     + ++Y+ +  +    C++ G  + +G   K",
        )
        self.assertEqual(
            str(hsp),
            """\
gi|254579         0 RVAFTFIRVEEKEYCY*EC*R*IN*NFESSSQL*GIIKTSRFYSTASDVMCIQECQIDYY
                  0 |.||...|.|.|......||.||.|..|.|...||..|||.|.........|.||..|..
Query_2           0 RIAFFIFRIEKKKFNHSPC***IH*DIEKST*F*GARKTSGFRTPFRVGLPIKEC*NDDP

gi|254579        60 INAMFPKIGHSAMYTGR*TLRRT*CVYRGQPAGNGYKTK 99
                 60 .|||.........|..........|...|.....|...| 99
Query_2          60 GNAMPTGTVNRSIYSSKPAV*NFGCLH*GYSSRDGDSIK 99
""",
        )


class TestBlastErrors(unittest.TestCase):
    """Tests if the Blast XML parser raises the appropriate Exception."""

    def test_not_xml(self):
        """Try to parse a FASTA file."""
        message = "Failed to parse the XML data (syntax error: line 1, column 0). Please make sure that the input data are in XML format."
        filename = "wisteria.nu"
        path = os.path.join("Fasta", filename)
        with open(path, "rb") as stream:
            with self.assertRaises(Blast.NotXMLError) as cm:
                records = Blast.parse(stream)
            self.assertEqual(str(cm.exception), message)
        with self.assertRaises(Blast.NotXMLError) as cm:
            records = Blast.parse(path)
        self.assertEqual(str(cm.exception), message)

    def test_premature_end_header(self):
        """Try to parse an XML file terminating in the header."""
        message = r"^premature end of XML file \(after reading [1-9]\d* bytes\)$"
        filename = "broken1.xml"
        path = os.path.join("Blast", filename)
        with open(path, "rb") as stream:
            with self.assertRaises(ValueError) as cm:
                records = Blast.parse(stream)
            self.assertRegex(str(cm.exception), message)
        with self.assertRaises(ValueError) as cm:
            records = Blast.parse(path)
        self.assertRegex(str(cm.exception), message)
        with open(path, "rb") as stream:
            with self.assertRaises(ValueError) as cm:
                records = Blast.read(stream)
            self.assertRegex(str(cm.exception), message)
        with self.assertRaises(ValueError) as cm:
            records = Blast.read(path)
        self.assertRegex(str(cm.exception), message)

    def test_premature_end_first_block(self):
        """Try to parse an XML file terminating within the first block."""
        message = r"^premature end of XML file \(after reading [1-9]\d* bytes\)$"
        filename = "broken2.xml"
        path = os.path.join("Blast", filename)
        with open(path, "rb") as stream:
            records = Blast.parse(stream)
            with self.assertRaises(ValueError) as cm:
                record = next(records)
            self.assertRegex(str(cm.exception), message)
        with Blast.parse(path) as records:
            with self.assertRaises(ValueError) as cm:
                record = next(records)
            self.assertRegex(str(cm.exception), message)
        with open(path, "rb") as stream:
            with self.assertRaises(ValueError) as cm:
                record = Blast.read(stream)
            self.assertRegex(str(cm.exception), message)
        with self.assertRaises(ValueError) as cm:
            record = Blast.read(path)
        self.assertRegex(str(cm.exception), message)

    def test_premature_end_second_block(self):
        """Try to parse an XML file terminating in the second block."""
        message = r"^premature end of XML file \(after reading [1-9]\d* bytes\)$"
        filename = "broken3.xml"
        path = os.path.join("Blast", filename)
        with open(path, "rb") as stream:
            records = Blast.parse(stream)
            with self.assertRaises(ValueError) as cm:
                record = next(records)
            self.assertRegex(str(cm.exception), message)
        with Blast.parse(path) as records:
            with self.assertRaises(ValueError) as cm:
                record = next(records)
            self.assertRegex(str(cm.exception), message)
        with open(path, "rb") as stream:
            with self.assertRaises(ValueError) as cm:
                record = Blast.read(stream)
            self.assertRegex(str(cm.exception), message)
        with self.assertRaises(ValueError) as cm:
            record = Blast.read(path)
        self.assertRegex(str(cm.exception), message)

    def test_premature_end_after_one_record(self):
        """Try to parse an XML file terminating after the first record."""
        message = r"^premature end of XML file \(after reading [1-9]\d* bytes\)$"
        filename = "broken4.xml"
        path = os.path.join("Blast", filename)
        with open(path, "rb") as stream:
            records = Blast.parse(stream)
            record = next(records)
            with self.assertRaises(ValueError) as cm:
                record = next(records)
            self.assertRegex(str(cm.exception), message)
        with Blast.parse(path) as records:
            record = next(records)
            with self.assertRaises(ValueError) as cm:
                record = next(records)
            self.assertRegex(str(cm.exception), message)
        with open(path, "rb") as stream:
            with self.assertRaises(ValueError) as cm:
                record = Blast.read(stream)
            self.assertRegex(str(cm.exception), message)
        with self.assertRaises(ValueError) as cm:
            record = Blast.read(path)
        self.assertRegex(str(cm.exception), message)

    def test_corrupt_xml(self):
        """Try to parse a broken XML file."""
        message = "Failed to parse the XML data (not well-formed (invalid token): line 10, column 2). Please make sure that the input data are not corrupted."

        filename = "broken5.xml"
        path = os.path.join("Blast", filename)
        with open(path, "rb") as stream:
            with self.assertRaises(Blast.CorruptedXMLError) as cm:
                records = Blast.parse(stream)
            self.assertEqual(str(cm.exception), message)
        with self.assertRaises(Blast.CorruptedXMLError) as cm:
            with Blast.parse(path) as records:
                pass
            self.assertEqual(str(cm.exception), message)
        with open(path, "rb") as stream:
            with self.assertRaises(Blast.CorruptedXMLError) as cm:
                record = Blast.read(stream)
            self.assertEqual(str(cm.exception), message)
        with self.assertRaises(Blast.CorruptedXMLError) as cm:
            record = Blast.read(path)
        self.assertEqual(str(cm.exception), message)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
