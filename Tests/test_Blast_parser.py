# Copyright 2005 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Test the Blast XML parser."""

import io
import os
import unittest

import numpy as np

from Bio import Blast
from Bio import StreamModeError
from Bio.SeqRecord import SeqRecord


class TestBlastp(unittest.TestCase):
    """Test the Blast XML parser for blastp output."""

    def check_xml_2218_blastp_002_header(self, records):
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

    def check_xml_2218_blastp_002_record_0(self, record):
        self.assertEqual(record.num, 1)
        self.assertEqual(
            repr(record),
            "<Bio.Blast.Record query.id='gi|585505|sp|Q08386|MOPB_RHOCA'; no hits>",
        )
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
        self.assertEqual(record.stat["eff-space"], 0)
        self.assertAlmostEqual(record.stat["kappa"], 0.041)
        self.assertAlmostEqual(record.stat["lambda"], 0.267)
        self.assertAlmostEqual(record.stat["entropy"], 0.14)
        self.assertEqual(len(record), 0)

    def check_xml_2218_blastp_002_record_1(self, record):
        self.assertEqual(record.num, 2)
        self.assertIsInstance(record.query, SeqRecord)
        self.assertEqual(record.query.id, "gi|129628|sp|P07175.1|PARA_AGRTU")
        self.assertEqual(record.query.description, "Protein parA")
        self.assertEqual(repr(record.query.seq), "Seq(None, length=222)")
        self.assertEqual(len(record.stat), 7)
        self.assertEqual(record.stat["db-num"], 27252)
        self.assertEqual(record.stat["db-len"], 13958303)
        self.assertEqual(record.stat["hsp-len"], 0)
        self.assertEqual(record.stat["eff-space"], 0)
        self.assertAlmostEqual(record.stat["kappa"], 0.041)
        self.assertAlmostEqual(record.stat["lambda"], 0.267)
        self.assertAlmostEqual(record.stat["entropy"], 0.14)
        self.assertEqual(len(record), 0)

    def test_xml_2218_blastp_002_iterator(self):
        """Parsing BLASTP 2.2.18+ (xml_2218_blastp_002.xml) by iteration."""
        filename = "xml_2218_blastp_002.xml"
        path = os.path.join("Blast", filename)
        with open(path, "rb") as stream:
            records = Blast.parse(stream)
            self.check_xml_2218_blastp_002_header(records)
            record = next(records)
            self.check_xml_2218_blastp_002_record_0(record)
            record = next(records)
            self.check_xml_2218_blastp_002_record_1(record)
        with open(path, "rb") as stream:
            records = Blast.parse(stream)
            records = records[:]
            self.check_xml_2218_blastp_002_header(records)
            record = next(records)
            self.check_xml_2218_blastp_002_record_0(record)
            record = next(records)
            self.check_xml_2218_blastp_002_record_1(record)
            self.assertRaises(StopIteration, next, records)
        with open(path) as stream:
            with self.assertRaises(StreamModeError) as cm:
                Blast.parse(stream)
            self.assertEqual(
                str(cm.exception), "BLAST output files must be opened in binary mode."
            )

    def test_xml_2218_blastp_002_list(self):
        """Parsing BLASTP 2.2.18+ (xml_2218_blastp_002.xml) as a list."""
        filename = "xml_2218_blastp_002.xml"
        path = os.path.join("Blast", filename)
        with open(path, "rb") as stream:
            records = Blast.parse(stream)
            self.check_xml_2218_blastp_002_header(records)
            record = records[0]
        # everything should have been read in by now
        self.check_xml_2218_blastp_002_record_0(record)
        record = records[1]
        self.check_xml_2218_blastp_002_record_1(record)
        # header should still be OK
        self.check_xml_2218_blastp_002_header(records)
        self.assertEqual(len(records), 2)
        self.assertEqual(
            str(records),
            """\
Program: BLASTP 2.2.18+
     db: gpipe/9606/Previous/protein

  Query: gi|585505|sp|Q08386|MOPB_RHOCA (length=270)
         Molybdenum-pterin-binding protein mopB >gi|310278|gb|AAA71913.1|
         molybdenum-pterin-binding protein
   Hits: No hits found

  Query: gi|129628|sp|P07175.1|PARA_AGRTU (length=222)
         Protein parA
   Hits: No hits found""",
        )
        # check if converting the records to a list does not lose the header:
        with open(path, "rb") as stream:
            records = Blast.parse(stream)
            records = records[:]
        self.check_xml_2218_blastp_002_header(records)

    def test_xml_2218_blastp_002_writer(self):
        """Writing BLASTP 2.2.18+ (xml_2218_blastp_002.xml)."""
        filename = "xml_2218_blastp_002.xml"
        path = os.path.join("Blast", filename)
        with open(path, "rb") as stream:
            records = Blast.parse(stream)
            records = records[:]
        stream = io.BytesIO()
        n = Blast.write(records, stream)
        self.assertEqual(n, 2)
        stream.seek(0)
        written_records = Blast.parse(stream)
        self.check_xml_2218_blastp_002_header(written_records)
        record = next(written_records)
        self.check_xml_2218_blastp_002_record_0(record)
        record = next(written_records)
        self.check_xml_2218_blastp_002_record_1(record)
        self.assertRaises(StopIteration, next, written_records)

    def test_xml_2218L_blastp_001_parser(self):
        """Parsing blastp 2.2.18 [Mar-02-2008] (xml_2218L_blastp_001.xml)."""
        filename = "xml_2218L_blastp_001.xml"
        path = os.path.join("Blast", filename)
        with open(path, "rb") as stream:
            records = Blast.parse(stream)
            self.check_xml_2218L_blastp_001_records(records)
            self.check_xml_2218L_blastp_001_str(records)

        with Blast.parse(path) as records:
            self.check_xml_2218L_blastp_001_records(records)
            self.check_xml_2218L_blastp_001_str(records)

        with open(path, "rb") as stream:
            record = Blast.read(stream)
        self.check_xml_2218L_blastp_001_record(record)

        record = Blast.read(path)
        self.check_xml_2218L_blastp_001_record(record)

    def check_xml_2218L_blastp_001_str(self, records):
        self.assertEqual(
            str(records),
            """\
Program: blastp 2.2.18 [Mar-02-2008]
     db: /Users/pjcock/Downloads/Software/blast-2.2.18/data/nr

   Hits: No hits found""",
        )

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
        self.assertAlmostEqual(record.stat["eff-space"], 2.02782e10)
        self.assertAlmostEqual(record.stat["kappa"], 0.041)
        self.assertAlmostEqual(record.stat["lambda"], 0.267)
        self.assertAlmostEqual(record.stat["entropy"], 0.14)
        self.assertEqual(len(record), 0)

    def test_xml_2218L_blastp_001_writer(self):
        """Writing blastp 2.2.18 [Mar-02-2008] (xml_2218L_blastp_001.xml)."""
        filename = "xml_2218L_blastp_001.xml"
        path = os.path.join("Blast", filename)
        with Blast.parse(path) as records:
            stream = io.BytesIO()
            n = Blast.write(records, stream)
            self.assertEqual(n, 1)
            stream.seek(0)
            written_records = Blast.parse(stream)
            self.check_xml_2218L_blastp_001_records(written_records)

    def test_xml_2226_blastp_003(self):
        """Parsing BLASTP 2.2.26+ (xml_2226_blastp_003.xml)."""
        filename = "xml_2226_blastp_003.xml"
        path = os.path.join("Blast", filename)
        with open(path, "rb") as stream:
            records = Blast.parse(stream)
            self.check_xml_2226_blastp_003(records)
        with Blast.parse(path) as records:
            self.check_xml_2226_blastp_003(records)
        with Blast.parse(path) as records:
            self.assertEqual(
                str(records),
                """\
Program: BLASTP 2.2.26+
     db: db/minirefseq_prot

  Query: Query_1 (length=102)
         gi|16080617|ref|NP_391444.1| membrane bound lipoprotein [Bacillus
         subtilis subsp. subtilis str. 168]
   Hits: ----  -----  ----------------------------------------------------------
            #  # HSP  ID + description
         ----  -----  ----------------------------------------------------------
            0      1  gnl|BL_ORD_ID|1  gi|308175296|ref|YP_003922001.1| membr...
            1      1  gnl|BL_ORD_ID|2  gi|375363999|ref|YP_005132038.1| lytA ...
            2      1  gnl|BL_ORD_ID|3  gi|154687679|ref|YP_001422840.1| LytA ...
            3      1  gnl|BL_ORD_ID|4  gi|311070071|ref|YP_003974994.1| unnam...
            4      1  gnl|BL_ORD_ID|15  gi|332258565|ref|XP_003278367.1| PRED...""",
            )
        record = Blast.read(path)
        self.assertEqual(
            str(record),
            """\
Program: BLASTP 2.2.26+
     db: db/minirefseq_prot
  Query: Query_1 (length=102)
         gi|16080617|ref|NP_391444.1| membrane bound lipoprotein [Bacillus
         subtilis subsp. subtilis str. 168]
   Hits: ----  -----  ----------------------------------------------------------
            #  # HSP  ID + description
         ----  -----  ----------------------------------------------------------
            0      1  gnl|BL_ORD_ID|1  gi|308175296|ref|YP_003922001.1| membr...
            1      1  gnl|BL_ORD_ID|2  gi|375363999|ref|YP_005132038.1| lytA ...
            2      1  gnl|BL_ORD_ID|3  gi|154687679|ref|YP_001422840.1| LytA ...
            3      1  gnl|BL_ORD_ID|4  gi|311070071|ref|YP_003974994.1| unnam...
            4      1  gnl|BL_ORD_ID|15  gi|332258565|ref|XP_003278367.1| PRED...""",
        )

    def check_xml_2226_blastp_003(self, records):
        self.assertEqual(records.program, "blastp")
        self.assertEqual(records.version, "BLASTP 2.2.26+")
        self.assertEqual(
            records.reference,
            'Stephen F. Altschul, Thomas L. Madden, Alejandro A. Schäffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.',
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
        self.assertEqual(record.num, 1)
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
        self.assertEqual(record.stat["eff-space"], 361344)
        self.assertAlmostEqual(record.stat["kappa"], 0.041)
        self.assertAlmostEqual(record.stat["lambda"], 0.267)
        self.assertAlmostEqual(record.stat["entropy"], 0.14)
        self.assertEqual(len(record), 5)
        hit = record[0]
        self.assertEqual(hit.num, 1)
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
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 350.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 139.428)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.99275e-46, places=51)
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
            repr(hsp),
            "<Bio.Blast.HSP target.id='gnl|BL_ORD_ID|1' query.id='Query_1'; 2 rows x 102 columns>",
        )
        self.assertEqual(
            str(hsp),
            """\
Query : Query_1 Length: 102 Strand: Plus
        gi|16080617|ref|NP_391444.1| membrane bound lipoprotein [Bacillus
        subtilis subsp. subtilis str. 168]
Target: gnl|BL_ORD_ID|1 Length: 100 Strand: Plus
        gi|308175296|ref|YP_003922001.1| membrane bound lipoprotein [Bacillus
        amyloliquefaciens DSM 7]

Score:139 bits(350), Expect:2e-46,
Identities:69/102(68%),  Positives:81/102(79%),  Gaps:2.102(2%)

gnl|BL_OR         0 MKKIFGCLFFILLLAGCGVTNEKSQGEDAG--EKLVTKEGTYVGLADTHTIEVTVDHEPV
                  0 |||....|||||||.||||...||||||..--....||||||||||||||||||||.|||
Query_1           0 MKKFIALLFFILLLSGCGVNSQKSQGEDVSPDSNIETKEGTYVGLADTHTIEVTVDNEPV

gnl|BL_OR        58 SFDITEESADDVKNLNNGEKVTVKYQKNSKGQLVLKDIEPAN 100
                 60 |.||||||..|....|.|.|||..|.||..|||.|||||.|| 102
Query_1          60 SLDITEESTSDLDKFNSGDKVTITYEKNDEGQLLLKDIERAN 102

""",
        )
        hit = record[1]
        self.assertEqual(hit.num, 2)
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
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 219.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 88.9669)
        self.assertAlmostEqual(hsp.annotations["evalue"], 6.94052e-27, places=32)
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
            repr(hsp),
            "<Bio.Blast.HSP target.id='gnl|BL_ORD_ID|2' query.id='Query_1'; 2 rows x 105 columns>",
        )
        self.assertEqual(
            str(hsp),
            """\
Query : Query_1 Length: 102 Strand: Plus
        gi|16080617|ref|NP_391444.1| membrane bound lipoprotein [Bacillus
        subtilis subsp. subtilis str. 168]
Target: gnl|BL_ORD_ID|2 Length: 105 Strand: Plus
        gi|375363999|ref|YP_005132038.1| lytA gene product [Bacillus
        amyloliquefaciens subsp. plantarum CAU B946]

Score:88 bits(219), Expect:7e-27,
Identities:48/105(46%),  Positives:69/105(66%),  Gaps:5.105(5%)

gnl|BL_OR         0 MKKTIAASFLILLFSVVLAACGTAEQSKKGSG-SSENQAQKETAYYVGMADTHTIEVKVD
                  0 |||.||..|.|||----|..||...|...|..-|...........|||.||||||||.||
Query_1           0 MKKFIALLFFILL----LSGCGVNSQKSQGEDVSPDSNIETKEGTYVGLADTHTIEVTVD

gnl|BL_OR        59 DQPVSFEFSDDFSDVLNKFSENDKVSITYFTNDKGQKEIKEIEKA 104
                 60 ..|||..........|.||...|||.|||..||.||...|.||.| 105
Query_1          56 NEPVSLDITEESTSDLDKFNSGDKVTITYEKNDEGQLLLKDIERA 101

""",
        )
        hit = record[2]
        self.assertEqual(hit.num, 3)
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
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 219.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 88.9669)
        self.assertAlmostEqual(hsp.annotations["evalue"], 8.41012e-27, places=32)
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
            repr(hsp),
            "<Bio.Blast.HSP target.id='gnl|BL_ORD_ID|3' query.id='Query_1'; 2 rows x 105 columns>",
        )
        self.assertEqual(
            str(hsp),
            """\
Query : Query_1 Length: 102 Strand: Plus
        gi|16080617|ref|NP_391444.1| membrane bound lipoprotein [Bacillus
        subtilis subsp. subtilis str. 168]
Target: gnl|BL_ORD_ID|3 Length: 105 Strand: Plus
        gi|154687679|ref|YP_001422840.1| LytA [Bacillus amyloliquefaciens FZB42]

Score:88 bits(219), Expect:8e-27,
Identities:48/105(46%),  Positives:69/105(66%),  Gaps:5.105(5%)

gnl|BL_OR         0 MKKTIAASFLILLFSVVLAACGTADQSKKGSG-SSENQAQKETAYYVGMADTHTIEVKVD
                  0 |||.||..|.|||----|..||...|...|..-|...........|||.||||||||.||
Query_1           0 MKKFIALLFFILL----LSGCGVNSQKSQGEDVSPDSNIETKEGTYVGLADTHTIEVTVD

gnl|BL_OR        59 DQPVSFEFSDDFSDVLNKFSENDKVSITYFTNDKGQKEIKEIEKA 104
                 60 ..|||..........|.||...|||.|||..||.||...|.||.| 105
Query_1          56 NEPVSLDITEESTSDLDKFNSGDKVTITYEKNDEGQLLLKDIERA 101

""",
        )
        hit = record[3]
        self.assertEqual(hit.num, 4)
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
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 204.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 83.1889)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.37847e-24, places=29)
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
            repr(hsp),
            "<Bio.Blast.HSP target.id='gnl|BL_ORD_ID|4' query.id='Query_1'; 2 rows x 104 columns>",
        )
        self.assertEqual(
            str(hsp),
            """\
Query : Query_1 Length: 102 Strand: Plus
        gi|16080617|ref|NP_391444.1| membrane bound lipoprotein [Bacillus
        subtilis subsp. subtilis str. 168]
Target: gnl|BL_ORD_ID|4 Length: 105 Strand: Plus
        gi|311070071|ref|YP_003974994.1| unnamed protein product [Bacillus
        atrophaeus 1942]

Score:83 bits(204), Expect:1e-24,
Identities:45/104(43%),  Positives:66/104(63%),  Gaps:5.104(5%)

gnl|BL_OR         0 MKKNVASSFLILLFSIILAACGTAEQSKEG-NGSSSSQVQNETAYYVGMADTHTIEVKID
                  0 |||..|..|.|||----|..||...|...|-..|..|........|||.||||||||..|
Query_1           0 MKKFIALLFFILL----LSGCGVNSQKSQGEDVSPDSNIETKEGTYVGLADTHTIEVTVD

gnl|BL_OR        59 DQPVSFEFTDDFSEILNEFEENDKVNISYLTNDKGQKELTEIEK 103
                 60 ..|||...|......|..|...|||.|.|..||.||..|..||. 104
Query_1          56 NEPVSLDITEESTSDLDKFNSGDKVTITYEKNDEGQLLLKDIER 100

""",
        )
        hit = record[4]
        self.assertEqual(hit.num, 5)
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
        self.assertEqual(hsp.num, 1)
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
            repr(hsp),
            "<Bio.Blast.HSP target.id='gnl|BL_ORD_ID|15' query.id='Query_1'; 2 rows x 25 columns>",
        )
        self.assertEqual(
            str(hsp),
            """\
Query : Query_1 Length: 102 Strand: Plus
        gi|16080617|ref|NP_391444.1| membrane bound lipoprotein [Bacillus
        subtilis subsp. subtilis str. 168]
Target: gnl|BL_ORD_ID|15 Length: 132 Strand: Plus
        gi|332258565|ref|XP_003278367.1| PREDICTED: UPF0764 protein
        C16orf89-like [Nomascus leucogenys]

Score:15 bits(29), Expect:7,
Identities:7/25(28%),  Positives:11/25(44%),  Gaps:0.25(0%)

gnl|BL_OR        79 VEMGFLHVGQAGLELVTSGDPPTLT 104
                  0 |...........|....|||..|.|  25
Query_1          59 VSLDITEESTSDLDKFNSGDKVTIT  84

""",
        )
        with self.assertRaises(IndexError) as cm:
            record[5]
        self.assertEqual(str(cm.exception), "index out of range")
        with self.assertRaises(TypeError) as cm:
            record[None]
        self.assertEqual(str(cm.exception), "key must be an integer, slice, or str")
        with self.assertRaises(KeyError) as cm:
            record["weird_key"]
        self.assertEqual(str(cm.exception), "'weird_key'")
        target_id = "gnl|BL_ORD_ID|4"
        self.assertIn(target_id, record)
        self.assertNotIn("weird_id", record)
        self.assertEqual(record[target_id].target.id, target_id)
        self.assertEqual(record.index(target_id), 3)
        with self.assertRaises(ValueError) as cm:
            record.index("weird_id")
        self.assertEqual(str(cm.exception), "'weird_id' not found")
        self.assertEqual(
            repr(hit),
            "<Bio.Blast.Hit target.id='gnl|BL_ORD_ID|15' query.id='Query_1'; 1 HSP>",
        )
        self.assertEqual(
            repr(hit[:0]), "<Bio.Blast.Hit target.id='gnl|BL_ORD_ID|15'; no hits>"
        )
        self.assertEqual(
            record.keys(),
            [
                "gnl|BL_ORD_ID|1",
                "gnl|BL_ORD_ID|2",
                "gnl|BL_ORD_ID|3",
                "gnl|BL_ORD_ID|4",
                "gnl|BL_ORD_ID|15",
            ],
        )

    def test_xml_2226_blastp_003_writer(self):
        """Writing BLASTP 2.2.26+ (xml_2226_blastp_003.xml)."""
        filename = "xml_2226_blastp_003.xml"
        path = os.path.join("Blast", filename)
        with Blast.parse(path) as records:
            stream = io.BytesIO()
            n = Blast.write(records, stream)
            self.assertEqual(n, 1)
            stream.seek(0)
            written_records = Blast.parse(stream)
            self.check_xml_2226_blastp_003(written_records)

    def test_phiblast_parser(self):
        """Parsing BLASTP 2.14.1+ (phiblast.xml)."""
        filename = "phiblast.xml"
        path = os.path.join("Blast", filename)
        with open(path, "rb") as stream:
            records = Blast.parse(stream)
            self.check_phiblast_records(records)
        with Blast.parse(path) as records:
            self.check_phiblast_records(records)
        with open(path, "rb") as stream:
            record = Blast.read(stream)
        self.check_phiblast_record(record)
        record = Blast.read(path)
        self.check_phiblast_record(record)

    def check_phiblast_records(self, records):
        self.assertEqual(records.program, "blastp")
        self.assertEqual(records.version, "BLASTP 2.14.1+")
        self.assertEqual(
            records.reference,
            'Zheng Zhang, Alejandro A. Schäffer, Webb Miller, Thomas L. Madden, David J. Lipman, Eugene V. Koonin, and Stephen F. Altschul (1998), "Protein sequence similarity searches using patterns as seeds", Nucleic Acids Res. 26:3986-3990.',
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
        self.assertEqual(record.num, 1)
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
        self.assertEqual(record.stat["eff-space"], 0)
        self.assertAlmostEqual(record.stat["kappa"], 0.047)
        self.assertAlmostEqual(record.stat["lambda"], 0.27)
        self.assertAlmostEqual(record.stat["entropy"], 1.0)
        self.assertEqual(len(record), 10)
        hit = record[0]
        self.assertEqual(hit.num, 1)
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
        self.assertEqual(hsp.num, 1)
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
Query : Query_74414 Length: 664 Strand: Plus
        unnamed protein product
Target: ref|NP_001075863.1| Length: 732 Strand: Plus
        cyclic nucleotide-gated olfactory channel [Oryctolagus cuniculus]
        >emb|CAA42201.1| aorta CNG channel (rACNG) [Oryctolagus cuniculus]
        >prf||1919268A cyclic nucleotide-gated channel [Oryctolagus cuniculus]

Score:1290 bits(3336), Expect:0,
Identities:664/664(100%),  Positives:664/664(100%),  Gaps:0.664(0%)

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
        self.assertEqual(hit.num, 2)
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
        self.assertEqual(hsp.num, 1)
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
Query : Query_74414 Length: 664 Strand: Plus
        unnamed protein product
Target: ref|XP_051689802.1| Length: 664 Strand: Plus
        cyclic nucleotide-gated olfactory channel isoform X3 [Oryctolagus
        cuniculus] >sp|Q28718.1| RecName: Full=Cyclic nucleotide-gated olfactory
        channel; AltName: Full=Aorta CNG channel; Short=RACNG; AltName:
        Full=Cyclic nucleotide-gated cation channel 2; AltName: Full=Cyclic
        nucleotide-gated channel alpha-2; Short=CNG channel alpha-2;
        Short=CNG-2; Short=CNG2 [Oryctolagus cuniculus]

Score:1290 bits(3336), Expect:0,
Identities:664/664(100%),  Positives:664/664(100%),  Gaps:0.664(0%)

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
        self.assertEqual(hit.num, 3)
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
        self.assertEqual(hsp.num, 1)
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
Query : Query_74414 Length: 664 Strand: Plus
        unnamed protein product
Target: ref|XP_017206345.1| Length: 677 Strand: Plus
        cyclic nucleotide-gated olfactory channel isoform X2 [Oryctolagus
        cuniculus]

Score:1290 bits(3336), Expect:0,
Identities:664/664(100%),  Positives:664/664(100%),  Gaps:0.664(0%)

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
        self.assertEqual(hit.num, 4)
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
        self.assertEqual(hsp.num, 1)
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
Query : Query_74414 Length: 664 Strand: Plus
        unnamed protein product
Target: ref|XP_051689801.1| Length: 687 Strand: Plus
        cyclic nucleotide-gated olfactory channel isoform X1 [Oryctolagus
        cuniculus]

Score:1290 bits(3336), Expect:0,
Identities:664/664(100%),  Positives:664/664(100%),  Gaps:0.664(0%)

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
        self.assertEqual(hit.num, 5)
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
        self.assertEqual(hsp.num, 1)
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
Query : Query_74414 Length: 664 Strand: Plus
        unnamed protein product
Target: ref|XP_004407164.1| Length: 664 Strand: Plus
        PREDICTED: cyclic nucleotide-gated olfactory channel [Odobenus rosmarus
        divergens]

Score:1249 bits(3231), Expect:0,
Identities:639/664(96%),  Positives:652/664(98%),  Gaps:0.664(0%)

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
        self.assertEqual(hit.num, 6)
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
        self.assertEqual(hsp.num, 1)
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
Query : Query_74414 Length: 664 Strand: Plus
        unnamed protein product
Target: ref|XP_008688471.1| Length: 664 Strand: Plus
        cyclic nucleotide-gated olfactory channel [Ursus maritimus]
        >ref|XP_026343324.1| cyclic nucleotide-gated olfactory channel [Ursus
        arctos]

Score:1248 bits(3228), Expect:0,
Identities:638/664(96%),  Positives:652/664(98%),  Gaps:0.664(0%)

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
        self.assertEqual(hit.num, 7)
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
        self.assertEqual(hsp.num, 1)
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
Query : Query_74414 Length: 664 Strand: Plus
        unnamed protein product
Target: ref|XP_011229794.1| Length: 664 Strand: Plus
        cyclic nucleotide-gated olfactory channel [Ailuropoda melanoleuca]
        >gb|EFB14215.1| hypothetical protein PANDA_013994, partial [Ailuropoda
        melanoleuca]

Score:1248 bits(3227), Expect:0,
Identities:638/664(96%),  Positives:652/664(98%),  Gaps:0.664(0%)

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
        self.assertEqual(hit.num, 8)
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
        self.assertEqual(hsp.num, 1)
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
Query : Query_74414 Length: 664 Strand: Plus
        unnamed protein product
Target: ref|XP_045646452.1| Length: 664 Strand: Plus
        cyclic nucleotide-gated olfactory channel [Ursus americanus]

Score:1246 bits(3223), Expect:0,
Identities:637/664(96%),  Positives:651/664(98%),  Gaps:0.664(0%)

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
        self.assertEqual(hit.num, 9)
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
        self.assertEqual(hsp.num, 1)
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
Query : Query_74414 Length: 664 Strand: Plus
        unnamed protein product
Target: ref|XP_035942617.1| Length: 664 Strand: Plus
        cyclic nucleotide-gated olfactory channel [Halichoerus grypus]

Score:1245 bits(3221), Expect:0,
Identities:638/664(96%),  Positives:651/664(98%),  Gaps:0.664(0%)

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
        self.assertEqual(hit.num, 10)
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
        self.assertEqual(hsp.num, 1)
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
Query : Query_74414 Length: 664 Strand: Plus
        unnamed protein product
Target: ref|XP_049729369.1| Length: 664 Strand: Plus
        cyclic nucleotide-gated olfactory channel [Elephas maximus indicus]

Score:1245 bits(3219), Expect:0,
Identities:635/664(96%),  Positives:654/664(98%),  Gaps:0.664(0%)

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

    def test_phiblast_writer(self):
        """Writing BLASTP 2.14.1+ (phiblast.xml)."""
        filename = "phiblast.xml"
        path = os.path.join("Blast", filename)
        with Blast.parse(path) as records:
            stream = io.BytesIO()
            n = Blast.write(records, stream)
            self.assertEqual(n, 1)
            stream.seek(0)
            written_records = Blast.parse(stream)
            self.check_phiblast_records(written_records)

    def test_xml_21500_blastp_001_parser(self):
        """Parsing BLASTP 2.15.0+ (xml_21500_blastp_001.xml)."""
        filename = "xml_21500_blastp_001.xml"
        path = os.path.join("Blast", filename)
        with open(path, "rb") as stream:
            records = Blast.parse(stream)
            self.check_xml_21500_blastp_001_records(records)
        with Blast.parse(path) as records:
            self.check_xml_21500_blastp_001_records(records)
        with open(path, "rb") as stream:
            record = Blast.read(stream)
        self.check_xml_21500_blastp_001_record(record)
        record = Blast.read(path)
        self.check_xml_21500_blastp_001_record(record)
        with Blast.parse(path) as records:
            self.assertEqual(
                str(records),
                """\
Program: BLASTP 2.15.0+
     db: nr

  Query: WXX52402.1 (length=239)
         RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]
   Hits: ----  -----  ----------------------------------------------------------
            #  # HSP  ID + description
         ----  -----  ----------------------------------------------------------
            0      1  ref|WP_003221446.1|  MULTISPECIES: RNA polymerase sporu...
            1      1  dbj|BAI85158.2|  sporulation sigma factor SigE [Bacillu...
            2      1  ref|WP_120028072.1|  RNA polymerase sporulation sigma f...
            3      1  ref|WP_326121348.1|  RNA polymerase sporulation sigma f...
            4      1  ref|WP_128473893.1|  RNA polymerase sporulation sigma f...
            5      1  ref|WP_174228079.1|  RNA polymerase sporulation sigma f...
            6      1  ref|WP_315947263.1|  RNA polymerase sporulation sigma f...
            7      1  ref|WP_219912761.1|  RNA polymerase sporulation sigma f...
            8      1  ref|WP_038828182.1|  RNA polymerase sporulation sigma f...
            9      1  ref|WP_326211018.1|  RNA polymerase sporulation sigma f...""",
            )
        record = Blast.read(path)
        self.assertEqual(
            str(record),
            """\
Program: BLASTP 2.15.0+
     db: nr
  Query: WXX52402.1 (length=239)
         RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]
   Hits: ----  -----  ----------------------------------------------------------
            #  # HSP  ID + description
         ----  -----  ----------------------------------------------------------
            0      1  ref|WP_003221446.1|  MULTISPECIES: RNA polymerase sporu...
            1      1  dbj|BAI85158.2|  sporulation sigma factor SigE [Bacillu...
            2      1  ref|WP_120028072.1|  RNA polymerase sporulation sigma f...
            3      1  ref|WP_326121348.1|  RNA polymerase sporulation sigma f...
            4      1  ref|WP_128473893.1|  RNA polymerase sporulation sigma f...
            5      1  ref|WP_174228079.1|  RNA polymerase sporulation sigma f...
            6      1  ref|WP_315947263.1|  RNA polymerase sporulation sigma f...
            7      1  ref|WP_219912761.1|  RNA polymerase sporulation sigma f...
            8      1  ref|WP_038828182.1|  RNA polymerase sporulation sigma f...
            9      1  ref|WP_326211018.1|  RNA polymerase sporulation sigma f...""",
        )

    def test_xml2_21500_blastp_001_parser(self):
        """Parsing BLASTP 2.15.0+ (xml2_21500_blastp_001.xml)."""
        filename = "xml2_21500_blastp_001.xml"
        path = os.path.join("Blast", filename)
        with open(path, "rb") as stream:
            records = Blast.parse(stream)
            self.check_xml_21500_blastp_001_records(records, xml2=True)
        with Blast.parse(path) as records:
            self.check_xml_21500_blastp_001_records(records, xml2=True)
        with open(path, "rb") as stream:
            record = Blast.read(stream)
        self.check_xml_21500_blastp_001_record(record, xml2=True)
        record = Blast.read(path)
        self.check_xml_21500_blastp_001_record(record, xml2=True)

    def check_xml_21500_blastp_001_records(self, records, xml2=False):
        self.assertEqual(records.program, "blastp")
        self.assertEqual(records.version, "BLASTP 2.15.0+")
        self.assertEqual(
            records.reference,
            'Stephen F. Altschul, Thomas L. Madden, Alejandro A. Schäffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.',
        )
        self.assertEqual(records.db, "nr")
        if not xml2:
            self.assertIsInstance(records.query, SeqRecord)
            self.assertEqual(records.query.id, "WXX52402.1")
            self.assertEqual(
                records.query.description,
                "RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]",
            )
            self.assertEqual(repr(records.query.seq), "Seq(None, length=239)")
        self.assertEqual(records.param["matrix"], "BLOSUM62")
        self.assertAlmostEqual(records.param["expect"], 0.05)
        self.assertEqual(records.param["gap-open"], 11)
        self.assertEqual(records.param["gap-extend"], 1)
        self.assertEqual(records.param["filter"], "F")
        if xml2:
            self.assertEqual(records.param["cbs"], 2)
            self.assertEqual(len(records.param), 6)
        else:
            self.assertEqual(len(records.param), 5)
        record = next(records)
        self.assertRaises(StopIteration, next, records)
        self.check_xml_21500_blastp_001_record(record, xml2=xml2)

    def check_xml_21500_blastp_001_record(self, record, xml2=False):
        hit = record[0]
        self.assertEqual(hit.num, 1)
        target = hit.target
        self.assertIsInstance(target, SeqRecord)
        self.assertEqual(target.id, "ref|WP_003221446.1|")
        self.assertEqual(target.name, "WP_003221446")
        seq = target.seq
        self.assertEqual(repr(seq), "Seq(None, length=239)")
        if xml2:
            self.assertEqual(
                target.description,
                "MULTISPECIES: RNA polymerase sporulation sigma factor SigE [Bacillales]",
            )
            self.assertEqual(target.annotations["taxid"], 1385)
            self.assertEqual(target.annotations["sciname"], "Bacillales")
            self.assertIs(target, hit.targets[0])
            self.assertEqual(len(hit.targets), 8)
            target = hit.targets[1]
            self.assertIsInstance(target, SeqRecord)
            self.assertEqual(target.id, "ref|NP_389415.2|")
            self.assertEqual(target.name, "NP_389415")
            self.assertIs(target.seq, seq)
            self.assertEqual(
                target.description,
                "RNA polymerase sporulation-specific sigma-29 factor (sigma-E) [Bacillus subtilis subsp. subtilis str. 168]",
            )
            self.assertEqual(target.annotations["taxid"], 224308)
            self.assertEqual(
                target.annotations["sciname"],
                "Bacillus subtilis subsp. subtilis str. 168",
            )
            target = hit.targets[2]
            self.assertIsInstance(target, SeqRecord)
            self.assertEqual(target.id, "sp|P06222.1|")
            self.assertEqual(target.name, "P06222")
            self.assertIs(target.seq, seq)
            self.assertEqual(
                target.description,
                "RecName: Full=RNA polymerase sigma-E factor; AltName: Full=P31; AltName: Full=Sigma-29; AltName: Full=Stage II sporulation protein GB; Flags: Precursor [Bacillus subtilis subsp. subtilis str. 168]",
            )
            self.assertEqual(target.annotations["taxid"], 224308)
            self.assertEqual(
                target.annotations["sciname"],
                "Bacillus subtilis subsp. subtilis str. 168",
            )
            target = hit.targets[3]
            self.assertIsInstance(target, SeqRecord)
            self.assertEqual(target.id, "gb|KFI04694.1|")
            self.assertEqual(target.name, "KFI04694")
            self.assertIs(target.seq, seq)
            self.assertEqual(
                target.description,
                "sporulation sigma factor SigE [Bacillus sp. BSC154]",
            )
            self.assertEqual(target.annotations["taxid"], 1549811)
            self.assertEqual(target.annotations["sciname"], "Bacillus sp. BSC154")
            target = hit.targets[4]
            self.assertIsInstance(target, SeqRecord)
            self.assertEqual(target.id, "gb|MDZ5720185.1|")
            self.assertEqual(target.name, "MDZ5720185")
            self.assertIs(target.seq, seq)
            self.assertEqual(
                target.description,
                "RNA polymerase sporulation sigma factor SigE [Bacillus sp. X(2023)]",
            )
            self.assertEqual(target.annotations["taxid"], 3106047)
            self.assertEqual(target.annotations["sciname"], "Bacillus sp. X(2023)")
            target = hit.targets[5]
            self.assertIsInstance(target, SeqRecord)
            self.assertEqual(target.id, "gb|POO83984.1|")
            self.assertEqual(target.name, "POO83984")
            self.assertIs(target.seq, seq)
            self.assertEqual(
                target.description,
                "RNA polymerase sporulation sigma factor SigE [Bacillus sp. MBGLi97]",
            )
            self.assertEqual(target.annotations["taxid"], 2070760)
            self.assertEqual(target.annotations["sciname"], "Bacillus sp. MBGLi97")
            target = hit.targets[6]
            self.assertIsInstance(target, SeqRecord)
            self.assertEqual(target.id, "dbj|BAM52179.1|")
            self.assertEqual(target.name, "BAM52179")
            self.assertIs(target.seq, seq)
            self.assertEqual(
                target.description,
                "sporulation sigma factor SigE [Bacillus subtilis BEST7613]",
            )
            self.assertEqual(target.annotations["taxid"], 1204343)
            self.assertEqual(
                target.annotations["sciname"], "Bacillus subtilis BEST7613"
            )
            target = hit.targets[7]
            self.assertIsInstance(target, SeqRecord)
            self.assertEqual(target.id, "gb|ADM37626.1|")
            self.assertEqual(target.name, "ADM37626")
            self.assertIs(target.seq, seq)
            self.assertEqual(
                target.description,
                "sporulation-specific sigma factor sigma-E [Bacillus spizizenii str. W23]",
            )
            self.assertEqual(target.annotations["taxid"], 655816)
            self.assertEqual(
                target.annotations["sciname"], "Bacillus spizizenii str. W23"
            )
        else:
            self.assertEqual(
                target.description,
                "MULTISPECIES: RNA polymerase sporulation sigma factor SigE [Bacillales] >ref|NP_389415.2| RNA polymerase sporulation-specific sigma-29 factor (sigma-E) [Bacillus subtilis subsp. subtilis str. 168] >sp|P06222.1| RecName: Full=RNA polymerase sigma-E factor; AltName: Full=P31; AltName: Full=Sigma-29; AltName: Full=Stage II sporulation protein GB; Flags: Precursor [Bacillus subtilis subsp. subtilis str. 168] >gb|KFI04694.1| sporulation sigma factor SigE [Bacillus sp. BSC154] >gb|MDZ5720185.1| RNA polymerase sporulation sigma factor SigE [Bacillus sp. X(2023)] >gb|POO83984.1| RNA polymerase sporulation sigma factor SigE [Bacillus sp. MBGLi97] >dbj|BAM52179.1| sporulation sigma factor SigE [Bacillus subtilis BEST7613] >gb|ADM37626.1| sporulation-specific sigma factor sigma-E [Bacillus spizizenii str. W23]",
            )
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 1227.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 477.248)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.44722e-169, places=174)
        self.assertEqual(hsp.annotations["identity"], 239)
        self.assertEqual(hsp.annotations["positive"], 239)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[0, 239],
                          [0, 239]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 239))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNG...KMV')",
        )
        self.assertEqual(hsp.query.id, "WXX52402.1")
        self.assertEqual(
            hsp.query.description,
            "RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]",
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNG...KMV')",
        )
        self.assertEqual(hsp.target.id, hit.target.id)
        self.assertEqual(hsp.target.name, hit.target.name)
        self.assertEqual(hsp.target.description, hit.target.description)
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARAILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIENEILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKALEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV",
        )
        self.assertEqual(
            repr(hsp),
            "<Bio.Blast.HSP target.id='ref|WP_003221446.1|' query.id='WXX52402.1'; 2 rows x 239 columns>",
        )
        if xml2:
            self.assertEqual(
                str(hsp),
                """\
Query : WXX52402.1 Length: 239 Strand: Plus
        RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]
Target: ref|WP_003221446.1| Length: 239 Strand: Plus
        MULTISPECIES: RNA polymerase sporulation sigma factor SigE [Bacillales]

Score:477 bits(1227), Expect:2e-169,
Identities:239/239(100%),  Positives:239/239(100%),  Gaps:0.239(0%)

ref|WP_00         0 MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARA
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.         0 MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARA

ref|WP_00        60 ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIEN
                 60 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.        60 ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIEN

ref|WP_00       120 EILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKA
                120 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.       120 EILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKA

ref|WP_00       180 LEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV
                180 |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.       180 LEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV

ref|WP_00       239
                239
WXX52402.       239

""",
            )
        else:
            self.assertEqual(
                str(hsp),
                """\
Query : WXX52402.1 Length: 239 Strand: Plus
        RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]
Target: ref|WP_003221446.1| Length: 239 Strand: Plus
        MULTISPECIES: RNA polymerase sporulation sigma factor SigE [Bacillales]
        >ref|NP_389415.2| RNA polymerase sporulation-specific sigma-29 factor
        (sigma-E) [Bacillus subtilis subsp. subtilis str. 168] >sp|P06222.1|
        RecName: Full=RNA polymerase sigma-E factor; AltName: Full=P31; AltName:
        Full=Sigma-29; AltName: Full=Stage II sporulation protein GB; Flags:
        Precursor [Bacillus subtilis subsp. subtilis str. 168] >gb|KFI04694.1|
        sporulation sigma factor SigE [Bacillus sp. BSC154] >gb|MDZ5720185.1|
        RNA polymerase sporulation sigma factor SigE [Bacillus sp. X(2023)]
        >gb|POO83984.1| RNA polymerase sporulation sigma factor SigE [Bacillus
        sp. MBGLi97] >dbj|BAM52179.1| sporulation sigma factor SigE [Bacillus
        subtilis BEST7613] >gb|ADM37626.1| sporulation-specific sigma factor
        sigma-E [Bacillus spizizenii str. W23]

Score:477 bits(1227), Expect:2e-169,
Identities:239/239(100%),  Positives:239/239(100%),  Gaps:0.239(0%)

ref|WP_00         0 MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARA
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.         0 MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARA

ref|WP_00        60 ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIEN
                 60 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.        60 ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIEN

ref|WP_00       120 EILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKA
                120 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.       120 EILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKA

ref|WP_00       180 LEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV
                180 |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.       180 LEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV

ref|WP_00       239
                239
WXX52402.       239

""",
            )
        hit = record[1]
        self.assertEqual(hit.num, 2)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "dbj|BAI85158.2|")
        self.assertEqual(hit.target.name, "BAI85158")
        self.assertEqual(
            hit.target.description,
            "sporulation sigma factor SigE [Bacillus subtilis subsp. natto BEST195]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=260)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 1227.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 477.248)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.91488e-169, places=174)
        self.assertEqual(hsp.annotations["identity"], 239)
        self.assertEqual(hsp.annotations["positive"], 239)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[21, 260],
                          [ 0, 239]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 239))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNG...KMV')",
        )
        self.assertEqual(hsp.query.id, "WXX52402.1")
        self.assertEqual(
            hsp.query.description,
            "RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]",
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({21: 'MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNG...KMV'}, length=260)",
        )
        self.assertEqual(hsp.target.id, "dbj|BAI85158.2|")
        self.assertEqual(hsp.target.name, "BAI85158")
        self.assertEqual(
            hsp.target.description,
            "sporulation sigma factor SigE [Bacillus subtilis subsp. natto BEST195]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARAILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIENEILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKALEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV",
        )
        self.assertEqual(
            repr(hsp),
            "<Bio.Blast.HSP target.id='dbj|BAI85158.2|' query.id='WXX52402.1'; 2 rows x 239 columns>",
        )
        self.assertEqual(
            str(hsp),
            """\
Query : WXX52402.1 Length: 239 Strand: Plus
        RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]
Target: dbj|BAI85158.2| Length: 260 Strand: Plus
        sporulation sigma factor SigE [Bacillus subtilis subsp. natto BEST195]

Score:477 bits(1227), Expect:4e-169,
Identities:239/239(100%),  Positives:239/239(100%),  Gaps:0.239(0%)

dbj|BAI85        21 MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARA
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.         0 MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARA

dbj|BAI85        81 ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIEN
                 60 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.        60 ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIEN

dbj|BAI85       141 EILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKA
                120 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.       120 EILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKA

dbj|BAI85       201 LEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV
                180 |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.       180 LEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV

dbj|BAI85       260
                239
WXX52402.       239

""",
        )
        hit = record[2]
        self.assertEqual(hit.num, 3)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "ref|WP_120028072.1|")
        self.assertEqual(hit.target.name, "WP_120028072")
        if xml2:
            self.assertEqual(
                hit.target.description,
                "RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]",
            )
        else:
            self.assertEqual(
                hit.target.description,
                "RNA polymerase sporulation sigma factor SigE [Bacillus subtilis] >gb|RJS52520.1| RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]",
            )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=239)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 1225.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 476.478)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.00949e-169, places=174)
        self.assertEqual(hsp.annotations["identity"], 238)
        self.assertEqual(hsp.annotations["positive"], 239)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[0, 239],
                          [0, 239]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 239))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNG...KMV')",
        )
        self.assertEqual(hsp.query.id, "WXX52402.1")
        self.assertEqual(
            hsp.query.description,
            "RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]",
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNG...KMV')",
        )
        self.assertEqual(hsp.target.id, hit.target.id)
        self.assertEqual(hsp.target.name, hit.target.name)
        self.assertEqual(hsp.target.description, hit.target.description)
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARAILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIENEILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKALEQLNEREKQIMELRFGL+GEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV",
        )
        self.assertEqual(
            repr(hsp),
            "<Bio.Blast.HSP target.id='ref|WP_120028072.1|' query.id='WXX52402.1'; 2 rows x 239 columns>",
        )
        if xml2:
            self.assertEqual(
                str(hsp),
                """\
Query : WXX52402.1 Length: 239 Strand: Plus
        RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]
Target: ref|WP_120028072.1| Length: 239 Strand: Plus
        RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]

Score:476 bits(1225), Expect:4e-169,
Identities:238/239(100%),  Positives:239/239(100%),  Gaps:0.239(0%)

ref|WP_12         0 MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARA
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.         0 MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARA

ref|WP_12        60 ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIEN
                 60 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.        60 ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIEN

ref|WP_12       120 EILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKA
                120 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.       120 EILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKA

ref|WP_12       180 LEQLNEREKQIMELRFGLIGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV
                180 ||||||||||||||||||.||||||||||||||||||||||||||||||||||||||||
WXX52402.       180 LEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV

ref|WP_12       239
                239
WXX52402.       239

""",
            )
        else:
            self.assertEqual(
                str(hsp),
                """\
Query : WXX52402.1 Length: 239 Strand: Plus
        RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]
Target: ref|WP_120028072.1| Length: 239 Strand: Plus
        RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]
        >gb|RJS52520.1| RNA polymerase sporulation sigma factor SigE [Bacillus
        subtilis]

Score:476 bits(1225), Expect:4e-169,
Identities:238/239(100%),  Positives:239/239(100%),  Gaps:0.239(0%)

ref|WP_12         0 MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARA
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.         0 MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARA

ref|WP_12        60 ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIEN
                 60 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.        60 ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIEN

ref|WP_12       120 EILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKA
                120 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.       120 EILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKA

ref|WP_12       180 LEQLNEREKQIMELRFGLIGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV
                180 ||||||||||||||||||.||||||||||||||||||||||||||||||||||||||||
WXX52402.       180 LEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV

ref|WP_12       239
                239
WXX52402.       239

""",
            )
        hit = record[3]
        self.assertEqual(hit.num, 4)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "ref|WP_326121348.1|")
        self.assertEqual(hit.target.name, "WP_326121348")
        if xml2:
            self.assertEqual(
                hit.target.description,
                "RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]",
            )
        else:
            self.assertEqual(
                hit.target.description,
                "RNA polymerase sporulation sigma factor SigE [Bacillus subtilis] >gb|MEC0320584.1| RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]",
            )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=239)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 1224.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 476.093)
        self.assertAlmostEqual(hsp.annotations["evalue"], 5.63375e-169, places=174)
        self.assertEqual(hsp.annotations["identity"], 238)
        self.assertEqual(hsp.annotations["positive"], 239)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[0, 239],
                          [0, 239]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 239))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNG...KMV')",
        )
        self.assertEqual(hsp.query.id, "WXX52402.1")
        self.assertEqual(
            hsp.query.description,
            "RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]",
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLIMKLPNG...KMV')",
        )
        self.assertEqual(hsp.target.id, hit.target.id)
        self.assertEqual(hsp.target.name, hit.target.name)
        self.assertEqual(hsp.target.description, hit.target.description)
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVL+MKLPNGDQAARAILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIENEILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKALEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV",
        )
        self.assertEqual(
            repr(hsp),
            "<Bio.Blast.HSP target.id='ref|WP_326121348.1|' query.id='WXX52402.1'; 2 rows x 239 columns>",
        )
        if xml2:
            self.assertEqual(
                str(hsp),
                """\
Query : WXX52402.1 Length: 239 Strand: Plus
        RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]
Target: ref|WP_326121348.1| Length: 239 Strand: Plus
        RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]

Score:476 bits(1224), Expect:6e-169,
Identities:238/239(100%),  Positives:239/239(100%),  Gaps:0.239(0%)

ref|WP_32         0 MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLIMKLPNGDQAARA
                  0 |||||||||||||||||||||||||||||||||||||||||||||||.||||||||||||
WXX52402.         0 MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARA

ref|WP_32        60 ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIEN
                 60 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.        60 ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIEN

ref|WP_32       120 EILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKA
                120 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.       120 EILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKA

ref|WP_32       180 LEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV
                180 |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.       180 LEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV

ref|WP_32       239
                239
WXX52402.       239

""",
            )
        else:
            self.assertEqual(
                str(hsp),
                """\
Query : WXX52402.1 Length: 239 Strand: Plus
        RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]
Target: ref|WP_326121348.1| Length: 239 Strand: Plus
        RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]
        >gb|MEC0320584.1| RNA polymerase sporulation sigma factor SigE [Bacillus
        subtilis]

Score:476 bits(1224), Expect:6e-169,
Identities:238/239(100%),  Positives:239/239(100%),  Gaps:0.239(0%)

ref|WP_32         0 MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLIMKLPNGDQAARA
                  0 |||||||||||||||||||||||||||||||||||||||||||||||.||||||||||||
WXX52402.         0 MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARA

ref|WP_32        60 ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIEN
                 60 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.        60 ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIEN

ref|WP_32       120 EILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKA
                120 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.       120 EILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKA

ref|WP_32       180 LEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV
                180 |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.       180 LEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV

ref|WP_32       239
                239
WXX52402.       239

""",
            )
        hit = record[4]
        self.assertEqual(hit.num, 5)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "ref|WP_128473893.1|")
        self.assertEqual(hit.target.name, "WP_128473893")
        if xml2:
            self.assertEqual(
                hit.target.description,
                "RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]",
            )
        else:
            self.assertEqual(
                hit.target.description,
                "RNA polymerase sporulation sigma factor SigE [Bacillus subtilis] >gb|QAR61557.1| RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]",
            )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=239)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 1224.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 476.093)
        self.assertAlmostEqual(hsp.annotations["evalue"], 5.6959e-169, places=173)
        self.assertEqual(hsp.annotations["identity"], 238)
        self.assertEqual(hsp.annotations["positive"], 239)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[0, 239],
                          [0, 239]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 239))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNG...KMV')",
        )
        self.assertEqual(hsp.query.id, "WXX52402.1")
        self.assertEqual(
            hsp.query.description,
            "RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]",
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MKKLKLRLTHLWYKLLMKLGMKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNG...KMV')",
        )
        self.assertEqual(hsp.target.id, hit.target.id)
        self.assertEqual(hsp.target.name, hit.target.name)
        self.assertEqual(hsp.target.description, hsp.target.description)
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MKKLKLRLTHLWYKLLMKLG+KSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARAILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIENEILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKALEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV",
        )
        self.assertEqual(
            repr(hsp),
            "<Bio.Blast.HSP target.id='ref|WP_128473893.1|' query.id='WXX52402.1'; 2 rows x 239 columns>",
        )
        if xml2:
            self.assertEqual(
                str(hsp),
                """\
Query : WXX52402.1 Length: 239 Strand: Plus
        RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]
Target: ref|WP_128473893.1| Length: 239 Strand: Plus
        RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]

Score:476 bits(1224), Expect:6e-169,
Identities:238/239(100%),  Positives:239/239(100%),  Gaps:0.239(0%)

ref|WP_12         0 MKKLKLRLTHLWYKLLMKLGMKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARA
                  0 ||||||||||||||||||||.|||||||||||||||||||||||||||||||||||||||
WXX52402.         0 MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARA

ref|WP_12        60 ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIEN
                 60 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.        60 ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIEN

ref|WP_12       120 EILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKA
                120 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.       120 EILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKA

ref|WP_12       180 LEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV
                180 |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.       180 LEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV

ref|WP_12       239
                239
WXX52402.       239

""",
            )
        else:
            self.assertEqual(
                str(hsp),
                """\
Query : WXX52402.1 Length: 239 Strand: Plus
        RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]
Target: ref|WP_128473893.1| Length: 239 Strand: Plus
        RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]
        >gb|QAR61557.1| RNA polymerase sporulation sigma factor SigE [Bacillus
        subtilis]

Score:476 bits(1224), Expect:6e-169,
Identities:238/239(100%),  Positives:239/239(100%),  Gaps:0.239(0%)

ref|WP_12         0 MKKLKLRLTHLWYKLLMKLGMKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARA
                  0 ||||||||||||||||||||.|||||||||||||||||||||||||||||||||||||||
WXX52402.         0 MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARA

ref|WP_12        60 ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIEN
                 60 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.        60 ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIEN

ref|WP_12       120 EILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKA
                120 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.       120 EILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKA

ref|WP_12       180 LEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV
                180 |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.       180 LEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV

ref|WP_12       239
                239
WXX52402.       239

""",
            )
        hit = record[5]
        self.assertEqual(hit.num, 6)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "ref|WP_174228079.1|")
        self.assertEqual(hit.target.name, "WP_174228079")
        if xml2:
            self.assertEqual(
                hit.target.description,
                "RNA polymerase sporulation sigma factor SigE [Bacillus tequilensis]",
            )
        else:
            self.assertEqual(
                hit.target.description,
                "RNA polymerase sporulation sigma factor SigE [Bacillus tequilensis] >gb|NTU26438.1| RNA polymerase sporulation sigma factor SigE [Bacillus tequilensis]",
            )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=239)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 1224.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 476.093)
        self.assertAlmostEqual(hsp.annotations["evalue"], 7.17171e-169, places=174)
        self.assertEqual(hsp.annotations["identity"], 238)
        self.assertEqual(hsp.annotations["positive"], 239)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[0, 239],
                          [0, 239]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 239))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNG...KMV')",
        )
        self.assertEqual(hsp.query.id, "WXX52402.1")
        self.assertEqual(
            hsp.query.description,
            "RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]",
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MKKLKLRLTHLWYRLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNG...KMV')",
        )
        self.assertEqual(hsp.target.id, hit.target.id)
        self.assertEqual(hsp.target.name, hit.target.name)
        self.assertEqual(hsp.target.description, hit.target.description)
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MKKLKLRLTHLWY+LLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARAILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIENEILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKALEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV",
        )
        self.assertEqual(
            repr(hsp),
            "<Bio.Blast.HSP target.id='ref|WP_174228079.1|' query.id='WXX52402.1'; 2 rows x 239 columns>",
        )
        if xml2:
            self.assertEqual(
                str(hsp),
                """\
Query : WXX52402.1 Length: 239 Strand: Plus
        RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]
Target: ref|WP_174228079.1| Length: 239 Strand: Plus
        RNA polymerase sporulation sigma factor SigE [Bacillus tequilensis]

Score:476 bits(1224), Expect:7e-169,
Identities:238/239(100%),  Positives:239/239(100%),  Gaps:0.239(0%)

ref|WP_17         0 MKKLKLRLTHLWYRLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARA
                  0 |||||||||||||.||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.         0 MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARA

ref|WP_17        60 ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIEN
                 60 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.        60 ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIEN

ref|WP_17       120 EILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKA
                120 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.       120 EILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKA

ref|WP_17       180 LEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV
                180 |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.       180 LEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV

ref|WP_17       239
                239
WXX52402.       239

""",
            )
        else:
            self.assertEqual(
                str(hsp),
                """\
Query : WXX52402.1 Length: 239 Strand: Plus
        RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]
Target: ref|WP_174228079.1| Length: 239 Strand: Plus
        RNA polymerase sporulation sigma factor SigE [Bacillus tequilensis]
        >gb|NTU26438.1| RNA polymerase sporulation sigma factor SigE [Bacillus
        tequilensis]

Score:476 bits(1224), Expect:7e-169,
Identities:238/239(100%),  Positives:239/239(100%),  Gaps:0.239(0%)

ref|WP_17         0 MKKLKLRLTHLWYRLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARA
                  0 |||||||||||||.||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.         0 MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARA

ref|WP_17        60 ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIEN
                 60 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.        60 ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIEN

ref|WP_17       120 EILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKA
                120 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.       120 EILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKA

ref|WP_17       180 LEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV
                180 |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.       180 LEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV

ref|WP_17       239
                239
WXX52402.       239

""",
            )
        hit = record[6]
        self.assertEqual(hit.num, 7)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "ref|WP_315947263.1|")
        self.assertEqual(hit.target.name, "WP_315947263")
        if xml2:
            self.assertEqual(
                hit.target.description,
                "RNA polymerase sporulation sigma factor SigE [Bacillus cabrialesii]",
            )
        else:
            self.assertEqual(
                hit.target.description,
                "RNA polymerase sporulation sigma factor SigE [Bacillus cabrialesii] >gb|MDU0153406.1| RNA polymerase sporulation sigma factor SigE [Bacillus cabrialesii]",
            )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=239)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 1224.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 476.093)
        self.assertAlmostEqual(hsp.annotations["evalue"], 7.17171e-169, places=174)
        self.assertEqual(hsp.annotations["identity"], 238)
        self.assertEqual(hsp.annotations["positive"], 239)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[0, 239],
                          [0, 239]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 239))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNG...KMV')",
        )
        self.assertEqual(hsp.query.id, "WXX52402.1")
        self.assertEqual(
            hsp.query.description,
            "RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]",
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MKKLKLRLTHLWYKLLMRLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNG...KMV')",
        )
        self.assertEqual(hsp.target.id, hit.target.id)
        self.assertEqual(hsp.target.name, hit.target.name)
        self.assertEqual(hsp.target.description, hit.target.description)
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MKKLKLRLTHLWYKLLM+LGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARAILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIENEILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKALEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV",
        )
        self.assertEqual(
            repr(hsp),
            "<Bio.Blast.HSP target.id='ref|WP_315947263.1|' query.id='WXX52402.1'; 2 rows x 239 columns>",
        )
        if xml2:
            self.assertEqual(
                str(hsp),
                """\
Query : WXX52402.1 Length: 239 Strand: Plus
        RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]
Target: ref|WP_315947263.1| Length: 239 Strand: Plus
        RNA polymerase sporulation sigma factor SigE [Bacillus cabrialesii]

Score:476 bits(1224), Expect:7e-169,
Identities:238/239(100%),  Positives:239/239(100%),  Gaps:0.239(0%)

ref|WP_31         0 MKKLKLRLTHLWYKLLMRLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARA
                  0 |||||||||||||||||.||||||||||||||||||||||||||||||||||||||||||
WXX52402.         0 MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARA

ref|WP_31        60 ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIEN
                 60 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.        60 ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIEN

ref|WP_31       120 EILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKA
                120 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.       120 EILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKA

ref|WP_31       180 LEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV
                180 |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.       180 LEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV

ref|WP_31       239
                239
WXX52402.       239

""",
            )
        else:
            self.assertEqual(
                str(hsp),
                """\
Query : WXX52402.1 Length: 239 Strand: Plus
        RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]
Target: ref|WP_315947263.1| Length: 239 Strand: Plus
        RNA polymerase sporulation sigma factor SigE [Bacillus cabrialesii]
        >gb|MDU0153406.1| RNA polymerase sporulation sigma factor SigE [Bacillus
        cabrialesii]

Score:476 bits(1224), Expect:7e-169,
Identities:238/239(100%),  Positives:239/239(100%),  Gaps:0.239(0%)

ref|WP_31         0 MKKLKLRLTHLWYKLLMRLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARA
                  0 |||||||||||||||||.||||||||||||||||||||||||||||||||||||||||||
WXX52402.         0 MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARA

ref|WP_31        60 ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIEN
                 60 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.        60 ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIEN

ref|WP_31       120 EILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKA
                120 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.       120 EILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKA

ref|WP_31       180 LEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV
                180 |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.       180 LEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV

ref|WP_31       239
                239
WXX52402.       239

""",
            )
        hit = record[7]
        self.assertEqual(hit.num, 8)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "ref|WP_219912761.1|")
        self.assertEqual(hit.target.name, "WP_219912761")
        if xml2:
            self.assertEqual(
                hit.target.description,
                "RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]",
            )
        else:
            self.assertEqual(
                hit.target.description,
                "RNA polymerase sporulation sigma factor SigE [Bacillus subtilis] >emb|COM82207.1| RNA polymerase sigma factor RpoD [Bacillus subtilis]",
            )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=239)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 1223.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 475.707)
        self.assertAlmostEqual(hsp.annotations["evalue"], 7.49345e-169, places=174)
        self.assertEqual(hsp.annotations["identity"], 238)
        self.assertEqual(hsp.annotations["positive"], 239)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[0, 239],
                          [0, 239]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 239))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNG...KMV')",
        )
        self.assertEqual(hsp.query.id, "WXX52402.1")
        self.assertEqual(
            hsp.query.description,
            "RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]",
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNG...KMV')",
        )
        self.assertEqual(hsp.target.id, hit.target.id)
        self.assertEqual(hsp.target.name, hit.target.name)
        self.assertEqual(hsp.target.description, hit.target.description)
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARAILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATY+SRCIENEILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKALEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV",
        )
        self.assertEqual(
            repr(hsp),
            "<Bio.Blast.HSP target.id='ref|WP_219912761.1|' query.id='WXX52402.1'; 2 rows x 239 columns>",
        )
        if xml2:
            self.assertEqual(
                str(hsp),
                """\
Query : WXX52402.1 Length: 239 Strand: Plus
        RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]
Target: ref|WP_219912761.1| Length: 239 Strand: Plus
        RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]

Score:475 bits(1223), Expect:7e-169,
Identities:238/239(100%),  Positives:239/239(100%),  Gaps:0.239(0%)

ref|WP_21         0 MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARA
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.         0 MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARA

ref|WP_21        60 ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYSSRCIEN
                 60 |||||||||||||||||||||||||||||||||||||||||||||||||||||.||||||
WXX52402.        60 ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIEN

ref|WP_21       120 EILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKA
                120 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.       120 EILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKA

ref|WP_21       180 LEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV
                180 |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.       180 LEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV

ref|WP_21       239
                239
WXX52402.       239

""",
            )
        else:
            self.assertEqual(
                str(hsp),
                """\
Query : WXX52402.1 Length: 239 Strand: Plus
        RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]
Target: ref|WP_219912761.1| Length: 239 Strand: Plus
        RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]
        >emb|COM82207.1| RNA polymerase sigma factor RpoD [Bacillus subtilis]

Score:475 bits(1223), Expect:7e-169,
Identities:238/239(100%),  Positives:239/239(100%),  Gaps:0.239(0%)

ref|WP_21         0 MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARA
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.         0 MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARA

ref|WP_21        60 ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYSSRCIEN
                 60 |||||||||||||||||||||||||||||||||||||||||||||||||||||.||||||
WXX52402.        60 ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIEN

ref|WP_21       120 EILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKA
                120 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.       120 EILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKA

ref|WP_21       180 LEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV
                180 |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.       180 LEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV

ref|WP_21       239
                239
WXX52402.       239

""",
            )
        hit = record[8]
        self.assertEqual(hit.num, 9)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "ref|WP_038828182.1|")
        self.assertEqual(hit.target.name, "WP_038828182")
        self.assertEqual(
            hit.target.description,
            "RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=239)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 1223.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 475.707)
        self.assertAlmostEqual(hsp.annotations["evalue"], 7.49345e-169, places=174)
        self.assertEqual(hsp.annotations["identity"], 238)
        self.assertEqual(hsp.annotations["positive"], 239)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[0, 239],
                          [0, 239]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 239))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNG...KMV')",
        )
        self.assertEqual(hsp.query.id, "WXX52402.1")
        self.assertEqual(
            hsp.query.description,
            "RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]",
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNG...KMV')",
        )
        self.assertEqual(hsp.target.id, "ref|WP_038828182.1|")
        self.assertEqual(hsp.target.name, "WP_038828182")
        self.assertEqual(
            hsp.target.description,
            "RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAAR+ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIENEILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKALEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV",
        )
        self.assertEqual(
            repr(hsp),
            "<Bio.Blast.HSP target.id='ref|WP_038828182.1|' query.id='WXX52402.1'; 2 rows x 239 columns>",
        )
        self.assertEqual(
            str(hsp),
            """\
Query : WXX52402.1 Length: 239 Strand: Plus
        RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]
Target: ref|WP_038828182.1| Length: 239 Strand: Plus
        RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]

Score:475 bits(1223), Expect:7e-169,
Identities:238/239(100%),  Positives:239/239(100%),  Gaps:0.239(0%)

ref|WP_03         0 MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARS
                  0 |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||.
WXX52402.         0 MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARA

ref|WP_03        60 ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIEN
                 60 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.        60 ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIEN

ref|WP_03       120 EILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKA
                120 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.       120 EILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKA

ref|WP_03       180 LEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV
                180 |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.       180 LEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV

ref|WP_03       239
                239
WXX52402.       239

""",
        )
        hit = record[9]
        self.assertEqual(hit.num, 10)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "ref|WP_326211018.1|")
        self.assertEqual(hit.target.name, "WP_326211018")
        if xml2:
            self.assertEqual(
                hit.target.description,
                "RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]",
            )
        else:
            self.assertEqual(
                hit.target.description,
                "RNA polymerase sporulation sigma factor SigE [Bacillus subtilis] >gb|MEC1541690.1| RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]",
            )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=239)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 1223.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 475.707)
        self.assertAlmostEqual(hsp.annotations["evalue"], 8.36237e-169, places=174)
        self.assertEqual(hsp.annotations["identity"], 238)
        self.assertEqual(hsp.annotations["positive"], 239)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[0, 239],
                          [0, 239]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 239))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNG...KMV')",
        )
        self.assertEqual(hsp.query.id, "WXX52402.1")
        self.assertEqual(
            hsp.query.description,
            "RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]",
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNG...KMM')",
        )
        self.assertEqual(hsp.target.id, hit.target.id)
        self.assertEqual(hsp.target.name, hit.target.name)
        self.assertEqual(hsp.target.description, hit.target.description)
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARAILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIENEILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKALEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKM+",
        )
        self.assertEqual(
            repr(hsp),
            "<Bio.Blast.HSP target.id='ref|WP_326211018.1|' query.id='WXX52402.1'; 2 rows x 239 columns>",
        )
        if xml2:
            self.assertEqual(
                str(hsp),
                """\
Query : WXX52402.1 Length: 239 Strand: Plus
        RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]
Target: ref|WP_326211018.1| Length: 239 Strand: Plus
        RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]

Score:475 bits(1223), Expect:8e-169,
Identities:238/239(100%),  Positives:239/239(100%),  Gaps:0.239(0%)

ref|WP_32         0 MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARA
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.         0 MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARA

ref|WP_32        60 ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIEN
                 60 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.        60 ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIEN

ref|WP_32       120 EILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKA
                120 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.       120 EILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKA

ref|WP_32       180 LEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMM
                180 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||.
WXX52402.       180 LEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV

ref|WP_32       239
                239
WXX52402.       239

""",
            )
        else:
            self.assertEqual(
                str(hsp),
                """\
Query : WXX52402.1 Length: 239 Strand: Plus
        RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]
Target: ref|WP_326211018.1| Length: 239 Strand: Plus
        RNA polymerase sporulation sigma factor SigE [Bacillus subtilis]
        >gb|MEC1541690.1| RNA polymerase sporulation sigma factor SigE [Bacillus
        subtilis]

Score:475 bits(1223), Expect:8e-169,
Identities:238/239(100%),  Positives:239/239(100%),  Gaps:0.239(0%)

ref|WP_32         0 MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARA
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.         0 MKKLKLRLTHLWYKLLMKLGLKSDEVYYIGGSEALPPPLSKDEEQVLLMKLPNGDQAARA

ref|WP_32        60 ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIEN
                 60 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.        60 ILIERNLRLVVYIARKFENTGINIEDLISIGTIGLIKAVNTFNPEKKIKLATYASRCIEN

ref|WP_32       120 EILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKA
                120 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
WXX52402.       120 EILMYLRRNNKIRSEVSFDEPLNIDWDGNELLLSDVLGTDDDIITKDIEANVDKKLLKKA

ref|WP_32       180 LEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMM
                180 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||.
WXX52402.       180 LEQLNEREKQIMELRFGLVGEEEKTQKDVADMMGISQSYISRLEKRIIKRLRKEFNKMV

ref|WP_32       239
                239
WXX52402.       239

""",
            )
        with self.assertRaises(IndexError) as cm:
            record[12]
        self.assertEqual(str(cm.exception), "index out of range")
        with self.assertRaises(TypeError) as cm:
            record[None]
        self.assertEqual(str(cm.exception), "key must be an integer, slice, or str")
        with self.assertRaises(KeyError) as cm:
            record["weird_key"]
        self.assertEqual(str(cm.exception), "'weird_key'")
        target_id = "ref|WP_326121348.1|"
        self.assertIn(target_id, record)
        self.assertNotIn("weird_id", record)
        self.assertEqual(record[target_id].target.id, target_id)
        self.assertEqual(record.index(target_id), 3)
        with self.assertRaises(ValueError) as cm:
            record.index("weird_id")
        self.assertEqual(str(cm.exception), "'weird_id' not found")
        self.assertEqual(
            repr(hit),
            "<Bio.Blast.Hit target.id='ref|WP_326211018.1|' query.id='WXX52402.1'; 1 HSP>",
        )
        self.assertEqual(
            repr(hit[:0]), "<Bio.Blast.Hit target.id='ref|WP_326211018.1|'; no hits>"
        )
        self.assertEqual(
            record.keys(),
            [
                "ref|WP_003221446.1|",
                "dbj|BAI85158.2|",
                "ref|WP_120028072.1|",
                "ref|WP_326121348.1|",
                "ref|WP_128473893.1|",
                "ref|WP_174228079.1|",
                "ref|WP_315947263.1|",
                "ref|WP_219912761.1|",
                "ref|WP_038828182.1|",
                "ref|WP_326211018.1|",
            ],
        )
        self.assertEqual(len(record.stat), 7)
        self.assertEqual(record.stat["db-num"], 718658123)
        self.assertEqual(record.stat["db-len"], 277418989154)
        if xml2:
            self.assertEqual(record.stat["hsp-len"], 161)
            self.assertEqual(record.stat["eff-space"], 12613772445378)
        else:
            self.assertEqual(record.stat["hsp-len"], 0)
            self.assertEqual(record.stat["eff-space"], 0)
        self.assertAlmostEqual(record.stat["kappa"], 0.041)
        self.assertAlmostEqual(record.stat["lambda"], 0.267)
        self.assertAlmostEqual(record.stat["entropy"], 0.14)

    def test_xml_21500_blastp_001_writer(self):
        """Writing BLASTP 2.15.0+ (xml_21500_blastp_001.xml)."""
        filename = "xml_21500_blastp_001.xml"
        path = os.path.join("Blast", filename)
        with Blast.parse(path) as records:
            stream = io.BytesIO()
            n = Blast.write(records, stream)
            self.assertEqual(n, 1)
            stream.seek(0)
            written_records = Blast.parse(stream)
            self.check_xml_21500_blastp_001_records(written_records)

    def test_xml2_21500_blastp_001_writer(self):
        """Writing BLASTP 2.15.0+ (xml2_21500_blastp_001.xml)."""
        filename = "xml2_21500_blastp_001.xml"
        path = os.path.join("Blast", filename)
        with Blast.parse(path) as records:
            stream = io.BytesIO()
            n = Blast.write(records, stream, fmt="XML2")
            self.assertEqual(n, 1)
            stream.seek(0)
            written_records = Blast.parse(stream)
            self.check_xml_21500_blastp_001_records(written_records, xml2=True)


class TestBlastn(unittest.TestCase):
    """Test the Blast XML parser for blastn output."""

    def test_xml_21500_blastn_001_parser(self):
        """Parsing BLASTN 2.15.0+ (xml_21500_blastn_001.xml)."""
        filename = "xml_21500_blastn_001.xml"
        path = os.path.join("Blast", filename)
        with open(path, "rb") as stream:
            records = Blast.parse(stream)
            self.check_xml_21500_blastn_001_records(records)
        with Blast.parse(path) as records:
            self.check_xml_21500_blastn_001_records(records)
        with open(path, "rb") as stream:
            record = Blast.read(stream)
        self.check_xml_21500_blastn_001_record(record)
        record = Blast.read(path)
        self.check_xml_21500_blastn_001_record(record)
        with Blast.parse(path) as records:
            self.assertEqual(
                str(records),
                """\
Program: BLASTN 2.15.0+
     db: genomic/10090/GCF_000001635.26

  Query: Query_78041 (length=285)
         G26684.1 human STS STS_D11570, sequence tagged site
   Hits: ----  -----  ----------------------------------------------------------
            #  # HSP  ID + description
         ----  -----  ----------------------------------------------------------
            0      1  gi|372099107|ref|NC_000069.6|  Mus musculus strain C57B...
            1      1  gi|372099103|ref|NC_000073.6|  Mus musculus strain C57B...
            2      2  gi|372099106|ref|NC_000070.6|  Mus musculus strain C57B...
            3      2  gi|372099108|ref|NC_000068.7|  Mus musculus strain C57B...
            4      2  gi|372099097|ref|NC_000079.6|  Mus musculus strain C57B...
            5      2  gi|372099098|ref|NC_000078.6|  Mus musculus strain C57B...
            6      1  gi|372099049|ref|NT_187008.1|  Mus musculus strain 129S...
            7      1  gi|372099109|ref|NC_000067.6|  Mus musculus strain C57B...
            8      1  gi|372099101|ref|NC_000075.6|  Mus musculus strain C57B...
            9      1  gi|372099100|ref|NC_000076.6|  Mus musculus strain C57B...
           10      1  gi|372099094|ref|NC_000082.6|  Mus musculus strain C57B...""",
            )

    def test_xml2_21500_blastn_001_parser(self):
        """Parsing BLASTN 2.15.0+ (xml2_21500_blastn_001.xml)."""
        filename = "xml2_21500_blastn_001.xml"
        path = os.path.join("Blast", filename)
        with open(path, "rb") as stream:
            records = Blast.parse(stream)
            self.check_xml_21500_blastn_001_records(records, xml2=True)
        with Blast.parse(path) as records:
            self.check_xml_21500_blastn_001_records(records, xml2=True)
        with open(path, "rb") as stream:
            record = Blast.read(stream)
        self.check_xml_21500_blastn_001_record(record, xml2=True)
        record = Blast.read(path)
        self.check_xml_21500_blastn_001_record(record, xml2=True)

    def check_xml_21500_blastn_001_records(self, records, xml2=False):
        self.assertEqual(records.program, "blastn")
        self.assertEqual(records.version, "BLASTN 2.15.0+")
        self.assertEqual(
            records.reference,
            'Stephen F. Altschul, Thomas L. Madden, Alejandro A. Schäffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.',
        )
        self.assertEqual(records.db, "genomic/10090/GCF_000001635.26")
        if not xml2:
            self.assertIsInstance(records.query, SeqRecord)
            self.assertEqual(records.query.id, "Query_78041")
            self.assertEqual(
                records.query.description,
                "G26684.1 human STS STS_D11570, sequence tagged site",
            )
            self.assertEqual(repr(records.query.seq), "Seq(None, length=285)")
        self.assertEqual(len(records.param), 6)
        self.assertAlmostEqual(records.param["expect"], 10.0)
        self.assertEqual(records.param["sc-match"], 2)
        self.assertEqual(records.param["sc-mismatch"], -3)
        self.assertEqual(records.param["gap-open"], 5)
        self.assertEqual(records.param["gap-extend"], 2)
        self.assertEqual(records.param["filter"], "L;m;")
        record = next(records)
        self.assertRaises(StopIteration, next, records)
        self.check_xml_21500_blastn_001_record(record, xml2=xml2)

    def check_xml_21500_blastn_001_record(self, record, xml2=False):
        if not xml2:
            self.assertEqual(record.num, 1)
        self.assertIsInstance(record.query, SeqRecord)
        self.assertEqual(record.query.id, "Query_78041")
        self.assertEqual(
            record.query.description,
            "G26684.1 human STS STS_D11570, sequence tagged site",
        )
        self.assertEqual(repr(record.query.seq), "Seq(None, length=285)")

        self.assertEqual(len(record.stat), 7)
        self.assertEqual(record.stat["db-num"], 239)
        self.assertEqual(record.stat["db-len"], 2818974565)
        if xml2:
            self.assertEqual(record.stat["hsp-len"], 31)
            self.assertEqual(record.stat["eff-space"], 716017657624)
        else:
            self.assertEqual(record.stat["hsp-len"], 0)
            self.assertEqual(record.stat["eff-space"], 0)
        self.assertAlmostEqual(record.stat["kappa"], 0.41)
        self.assertAlmostEqual(record.stat["lambda"], 0.625)
        self.assertAlmostEqual(record.stat["entropy"], 0.78)
        self.assertEqual(len(record), 11)
        hit = record[0]
        self.assertEqual(hit.num, 1)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|372099107|ref|NC_000069.6|")
        self.assertEqual(hit.target.name, "NC_000069")
        self.assertEqual(
            hit.target.description,
            "Mus musculus strain C57BL/6J chromosome 3, GRCm38.p6 C57BL/6J",
        )
        if xml2:
            self.assertEqual(hit.target.annotations["taxid"], 10090)
            self.assertEqual(hit.target.annotations["sciname"], "Mus musculus")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=160039680)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 44.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 40.9604)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.334664)
        self.assertEqual(hsp.annotations["identity"], 30)
        if not xml2:
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
        self.assertEqual(hsp.query.id, "Query_78041")
        self.assertEqual(
            hsp.query.description, "G26684.1 human STS STS_D11570, sequence tagged site"
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
            "Mus musculus strain C57BL/6J chromosome 3, GRCm38.p6 C57BL/6J",
        )
        if xml2:
            self.assertEqual(hit.target.annotations["taxid"], 10090)
            self.assertEqual(hit.target.annotations["sciname"], "Mus musculus")
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"], "|||||||||||||| || |||||| || ||||||"
        )
        if xml2:
            self.assertEqual(
                str(hsp),
                """\
Query : Query_78041 Length: 285 Strand: Plus
        G26684.1 human STS STS_D11570, sequence tagged site
Target: gi|372099107|ref|NC_000069.6| Length: 160039680 Strand: Minus
        Mus musculus strain C57BL/6J chromosome 3, GRCm38.p6 C57BL/6J

Score:40 bits(44), Expect:0.3,
Identities:30/34(88%),  Gaps:1.34(3%)

gi|372099 101449177 GAATCCTAGAGGCTGGACTGGCCCTGGCCTGCTG 101449143
                  0 ||||||||||||||.||.||||||.||-||||||        34
Query_780       133 GAATCCTAGAGGCTTGATTGGCCCAGG-CTGCTG       166

""",
            )
        else:
            self.assertEqual(
                str(hsp),
                """\
Query : Query_78041 Length: 285 Strand: Plus
        G26684.1 human STS STS_D11570, sequence tagged site
Target: gi|372099107|ref|NC_000069.6| Length: 160039680 Strand: Minus
        Mus musculus strain C57BL/6J chromosome 3, GRCm38.p6 C57BL/6J

Score:40 bits(44), Expect:0.3,
Identities:30/34(88%),  Positives:30/34(88%),  Gaps:1.34(3%)

gi|372099 101449177 GAATCCTAGAGGCTGGACTGGCCCTGGCCTGCTG 101449143
                  0 ||||||||||||||.||.||||||.||-||||||        34
Query_780       133 GAATCCTAGAGGCTTGATTGGCCCAGG-CTGCTG       166

""",
            )
        hit = record[1]
        self.assertEqual(hit.num, 2)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|372099103|ref|NC_000073.6|")
        self.assertEqual(hit.target.name, "NC_000073")
        self.assertEqual(
            hit.target.description,
            "Mus musculus strain C57BL/6J chromosome 7, GRCm38.p6 C57BL/6J",
        )
        if xml2:
            self.assertEqual(hit.target.annotations["taxid"], 10090)
            self.assertEqual(hit.target.annotations["sciname"], "Mus musculus")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=145441459)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 44.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 40.9604)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.334664)
        self.assertEqual(hsp.annotations["identity"], 26)
        if not xml2:
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
        self.assertEqual(hsp.query.id, "Query_78041")
        self.assertEqual(
            hsp.query.description, "G26684.1 human STS STS_D11570, sequence tagged site"
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
            "Mus musculus strain C57BL/6J chromosome 7, GRCm38.p6 C57BL/6J",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(hsp.annotations["midline"], "|||||||||  ||||||||||||||| ||")
        if xml2:
            self.assertEqual(
                str(hsp),
                """\
Query : Query_78041 Length: 285 Strand: Plus
        G26684.1 human STS STS_D11570, sequence tagged site
Target: gi|372099103|ref|NC_000073.6| Length: 145441459 Strand: Minus
        Mus musculus strain C57BL/6J chromosome 7, GRCm38.p6 C57BL/6J

Score:40 bits(44), Expect:0.3,
Identities:26/29(90%),  Gaps:0.29(0%)

gi|372099 131772185 GAAAGGAAAAAAAAATGGAAAGTTCTGGT 131772156
                  0 |||||||||..|||||||||||||||.||        29
Query_780       204 GAAAGGAAATNAAAATGGAAAGTTCTTGT       233

""",
            )
        else:
            self.assertEqual(
                str(hsp),
                """\
Query : Query_78041 Length: 285 Strand: Plus
        G26684.1 human STS STS_D11570, sequence tagged site
Target: gi|372099103|ref|NC_000073.6| Length: 145441459 Strand: Minus
        Mus musculus strain C57BL/6J chromosome 7, GRCm38.p6 C57BL/6J

Score:40 bits(44), Expect:0.3,
Identities:26/29(90%),  Positives:26/29(90%),  Gaps:0.29(0%)

gi|372099 131772185 GAAAGGAAAAAAAAATGGAAAGTTCTGGT 131772156
                  0 |||||||||..|||||||||||||||.||        29
Query_780       204 GAAAGGAAATNAAAATGGAAAGTTCTTGT       233

""",
            )
        hit = record[2]
        self.assertEqual(hit.num, 3)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|372099106|ref|NC_000070.6|")
        self.assertEqual(hit.target.name, "NC_000070")
        self.assertEqual(
            hit.target.description,
            "Mus musculus strain C57BL/6J chromosome 4, GRCm38.p6 C57BL/6J",
        )
        if xml2:
            self.assertEqual(hit.target.annotations["taxid"], 10090)
            self.assertEqual(hit.target.annotations["sciname"], "Mus musculus")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=156508116)")
        self.assertEqual(len(hit), 2)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 43.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 40.0587)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.16809)
        self.assertEqual(hsp.annotations["identity"], 23)
        if not xml2:
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
        self.assertEqual(hsp.query.id, "Query_78041")
        self.assertEqual(
            hsp.query.description, "G26684.1 human STS STS_D11570, sequence tagged site"
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
            "Mus musculus strain C57BL/6J chromosome 4, GRCm38.p6 C57BL/6J",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(hsp.annotations["midline"], "|||||||||||||||| |||||||")
        if xml2:
            self.assertEqual(
                str(hsp),
                """\
Query : Query_78041 Length: 285 Strand: Plus
        G26684.1 human STS STS_D11570, sequence tagged site
Target: gi|372099106|ref|NC_000070.6| Length: 156508116 Strand: Minus
        Mus musculus strain C57BL/6J chromosome 4, GRCm38.p6 C57BL/6J

Score:40 bits(43), Expect:1,
Identities:23/24(96%),  Gaps:0.24(0%)

gi|372099   9607562 CCAACACAGGCCAGCGGCTTCTGG 9607538
                  0 ||||||||||||||||.|||||||      24
Query_780        61 CCAACACAGGCCAGCGACTTCTGG      85

""",
            )
        else:
            self.assertEqual(
                str(hsp),
                """\
Query : Query_78041 Length: 285 Strand: Plus
        G26684.1 human STS STS_D11570, sequence tagged site
Target: gi|372099106|ref|NC_000070.6| Length: 156508116 Strand: Minus
        Mus musculus strain C57BL/6J chromosome 4, GRCm38.p6 C57BL/6J

Score:40 bits(43), Expect:1,
Identities:23/24(96%),  Positives:23/24(96%),  Gaps:0.24(0%)

gi|372099   9607562 CCAACACAGGCCAGCGGCTTCTGG 9607538
                  0 ||||||||||||||||.|||||||      24
Query_780        61 CCAACACAGGCCAGCGACTTCTGG      85

""",
            )
        hsp = hit[1]
        self.assertEqual(hsp.num, 2)
        self.assertAlmostEqual(hsp.score, 40.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 37.3537)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.07705)
        self.assertEqual(hsp.annotations["identity"], 28)
        if not xml2:
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
        self.assertEqual(hsp.query.id, "Query_78041")
        self.assertEqual(
            hsp.query.description, "G26684.1 human STS STS_D11570, sequence tagged site"
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
            "Mus musculus strain C57BL/6J chromosome 4, GRCm38.p6 C57BL/6J",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(hsp.annotations["midline"], "||||| ||||  ||| ||||||||||||||||")
        if xml2:
            self.assertEqual(
                str(hsp),
                """\
Query : Query_78041 Length: 285 Strand: Plus
        G26684.1 human STS STS_D11570, sequence tagged site
Target: gi|372099106|ref|NC_000070.6| Length: 156508116 Strand: Plus
        Mus musculus strain C57BL/6J chromosome 4, GRCm38.p6 C57BL/6J

Score:37 bits(40), Expect:4,
Identities:28/32(88%),  Gaps:1.32(3%)

gi|372099 142902531 GCCTGGCATGAAGTAACTGCTCAATAAATGCT 142902563
                  0 |||||.||||.-|||.||||||||||||||||        32
Query_780       241 GCCTGACATGG-GTAGCTGCTCAATAAATGCT       272

""",
            )
        else:
            self.assertEqual(
                str(hsp),
                """\
Query : Query_78041 Length: 285 Strand: Plus
        G26684.1 human STS STS_D11570, sequence tagged site
Target: gi|372099106|ref|NC_000070.6| Length: 156508116 Strand: Plus
        Mus musculus strain C57BL/6J chromosome 4, GRCm38.p6 C57BL/6J

Score:37 bits(40), Expect:4,
Identities:28/32(88%),  Positives:28/32(88%),  Gaps:1.32(3%)

gi|372099 142902531 GCCTGGCATGAAGTAACTGCTCAATAAATGCT 142902563
                  0 |||||.||||.-|||.||||||||||||||||        32
Query_780       241 GCCTGACATGG-GTAGCTGCTCAATAAATGCT       272

""",
            )
        hit = record[3]
        self.assertEqual(hit.num, 4)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|372099108|ref|NC_000068.7|")
        self.assertEqual(hit.target.name, "NC_000068")
        self.assertEqual(
            hit.target.description,
            "Mus musculus strain C57BL/6J chromosome 2, GRCm38.p6 C57BL/6J",
        )
        if xml2:
            self.assertEqual(hit.target.annotations["taxid"], 10090)
            self.assertEqual(hit.target.annotations["sciname"], "Mus musculus")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=182113224)")
        self.assertEqual(len(hit), 2)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 42.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 39.157)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.16809)
        self.assertEqual(hsp.annotations["identity"], 27)
        if not xml2:
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
        self.assertEqual(hsp.query.id, "Query_78041")
        self.assertEqual(
            hsp.query.description, "G26684.1 human STS STS_D11570, sequence tagged site"
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
            "Mus musculus strain C57BL/6J chromosome 2, GRCm38.p6 C57BL/6J",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(hsp.annotations["midline"], "||| |||| |||| |||| ||||||||||||")
        if xml2:
            self.assertEqual(
                str(hsp),
                """\
Query : Query_78041 Length: 285 Strand: Plus
        G26684.1 human STS STS_D11570, sequence tagged site
Target: gi|372099108|ref|NC_000068.7| Length: 182113224 Strand: Plus
        Mus musculus strain C57BL/6J chromosome 2, GRCm38.p6 C57BL/6J

Score:39 bits(42), Expect:1,
Identities:27/31(87%),  Gaps:0.31(0%)

gi|372099   3799646 AAGTCCTGGCATGAGTAGTTGCTCAATAAAT 3799677
                  0 |||.||||.||||.||||.||||||||||||      31
Query_780       238 AAGGCCTGACATGGGTAGCTGCTCAATAAAT     269

""",
            )
        else:
            self.assertEqual(
                str(hsp),
                """\
Query : Query_78041 Length: 285 Strand: Plus
        G26684.1 human STS STS_D11570, sequence tagged site
Target: gi|372099108|ref|NC_000068.7| Length: 182113224 Strand: Plus
        Mus musculus strain C57BL/6J chromosome 2, GRCm38.p6 C57BL/6J

Score:39 bits(42), Expect:1,
Identities:27/31(87%),  Positives:27/31(87%),  Gaps:0.31(0%)

gi|372099   3799646 AAGTCCTGGCATGAGTAGTTGCTCAATAAAT 3799677
                  0 |||.||||.||||.||||.||||||||||||      31
Query_780       238 AAGGCCTGACATGGGTAGCTGCTCAATAAAT     269

""",
            )
        hsp = hit[1]
        self.assertEqual(hsp.num, 2)
        self.assertAlmostEqual(hsp.score, 41.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 38.2554)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.07705)
        self.assertEqual(hsp.annotations["identity"], 23)
        if not xml2:
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
        self.assertEqual(hsp.query.id, "Query_78041")
        self.assertEqual(
            hsp.query.description, "G26684.1 human STS STS_D11570, sequence tagged site"
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
            "Mus musculus strain C57BL/6J chromosome 2, GRCm38.p6 C57BL/6J",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(hsp.annotations["midline"], "|||| |||||||||||||||| |||")
        if xml2:
            self.assertEqual(
                str(hsp),
                """\
Query : Query_78041 Length: 285 Strand: Plus
        G26684.1 human STS STS_D11570, sequence tagged site
Target: gi|372099108|ref|NC_000068.7| Length: 182113224 Strand: Plus
        Mus musculus strain C57BL/6J chromosome 2, GRCm38.p6 C57BL/6J

Score:38 bits(41), Expect:4,
Identities:23/25(92%),  Gaps:0.25(0%)

gi|372099  70278959 AAATGAAAATGGAAAGTTCTTATAG 70278984
                  0 ||||.||||||||||||||||.|||       25
Query_780       210 AAATNAAAATGGAAAGTTCTTGTAG      235

""",
            )
        else:
            self.assertEqual(
                str(hsp),
                """\
Query : Query_78041 Length: 285 Strand: Plus
        G26684.1 human STS STS_D11570, sequence tagged site
Target: gi|372099108|ref|NC_000068.7| Length: 182113224 Strand: Plus
        Mus musculus strain C57BL/6J chromosome 2, GRCm38.p6 C57BL/6J

Score:38 bits(41), Expect:4,
Identities:23/25(92%),  Positives:23/25(92%),  Gaps:0.25(0%)

gi|372099  70278959 AAATGAAAATGGAAAGTTCTTATAG 70278984
                  0 ||||.||||||||||||||||.|||       25
Query_780       210 AAATNAAAATGGAAAGTTCTTGTAG      235

""",
            )
        hit = record[4]
        self.assertEqual(hit.num, 5)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|372099097|ref|NC_000079.6|")
        self.assertEqual(hit.target.name, "NC_000079")
        self.assertEqual(
            hit.target.description,
            "Mus musculus strain C57BL/6J chromosome 13, GRCm38.p6 C57BL/6J",
        )
        if xml2:
            self.assertEqual(hit.target.annotations["taxid"], 10090)
            self.assertEqual(hit.target.annotations["sciname"], "Mus musculus")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=120421639)")
        self.assertEqual(len(hit), 2)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 42.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 39.157)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.16809)
        self.assertEqual(hsp.annotations["identity"], 25)
        if not xml2:
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
        self.assertEqual(hsp.query.id, "Query_78041")
        self.assertEqual(
            hsp.query.description, "G26684.1 human STS STS_D11570, sequence tagged site"
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
            "Mus musculus strain C57BL/6J chromosome 13, GRCm38.p6 C57BL/6J",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(hsp.annotations["midline"], "||||| || |||||||||||||||| ||")
        if xml2:
            self.assertEqual(
                str(hsp),
                """\
Query : Query_78041 Length: 285 Strand: Plus
        G26684.1 human STS STS_D11570, sequence tagged site
Target: gi|372099097|ref|NC_000079.6| Length: 120421639 Strand: Minus
        Mus musculus strain C57BL/6J chromosome 13, GRCm38.p6 C57BL/6J

Score:39 bits(42), Expect:1,
Identities:25/28(89%),  Gaps:0.28(0%)

gi|372099  26806584 AAGGACATCAAAATGGAAAGTTCTTCTA 26806556
                  0 |||||.||.||||||||||||||||.||       28
Query_780       206 AAGGAAATNAAAATGGAAAGTTCTTGTA      234

""",
            )
        else:
            self.assertEqual(
                str(hsp),
                """\
Query : Query_78041 Length: 285 Strand: Plus
        G26684.1 human STS STS_D11570, sequence tagged site
Target: gi|372099097|ref|NC_000079.6| Length: 120421639 Strand: Minus
        Mus musculus strain C57BL/6J chromosome 13, GRCm38.p6 C57BL/6J

Score:39 bits(42), Expect:1,
Identities:25/28(89%),  Positives:25/28(89%),  Gaps:0.28(0%)

gi|372099  26806584 AAGGACATCAAAATGGAAAGTTCTTCTA 26806556
                  0 |||||.||.||||||||||||||||.||       28
Query_780       206 AAGGAAATNAAAATGGAAAGTTCTTGTA      234

""",
            )
        hsp = hit[1]
        self.assertEqual(hsp.num, 2)
        self.assertAlmostEqual(hsp.score, 40.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 37.3537)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.07705)
        self.assertEqual(hsp.annotations["identity"], 32)
        if not xml2:
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
        self.assertEqual(hsp.query.id, "Query_78041")
        self.assertEqual(
            hsp.query.description, "G26684.1 human STS STS_D11570, sequence tagged site"
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
            "Mus musculus strain C57BL/6J chromosome 13, GRCm38.p6 C57BL/6J",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"], "||||||||||||||||| || |  || ||| | ||| |||"
        )
        if xml2:
            self.assertEqual(
                str(hsp),
                """\
Query : Query_78041 Length: 285 Strand: Plus
        G26684.1 human STS STS_D11570, sequence tagged site
Target: gi|372099097|ref|NC_000079.6| Length: 120421639 Strand: Minus
        Mus musculus strain C57BL/6J chromosome 13, GRCm38.p6 C57BL/6J

Score:37 bits(40), Expect:4,
Identities:32/40(80%),  Gaps:0.40(0%)

gi|372099  56840340 AGCGCAAGGCCTGACATAGGAAAATGTTCAGTGAATACTA 56840300
                  0 |||||||||||||||||.||.|..||.|||.|.|||.|||       40
Query_780       233 AGCGCAAGGCCTGACATGGGTAGCTGCTCAATAAATGCTA      273

""",
            )
        else:
            self.assertEqual(
                str(hsp),
                """\
Query : Query_78041 Length: 285 Strand: Plus
        G26684.1 human STS STS_D11570, sequence tagged site
Target: gi|372099097|ref|NC_000079.6| Length: 120421639 Strand: Minus
        Mus musculus strain C57BL/6J chromosome 13, GRCm38.p6 C57BL/6J

Score:37 bits(40), Expect:4,
Identities:32/40(80%),  Positives:32/40(80%),  Gaps:0.40(0%)

gi|372099  56840340 AGCGCAAGGCCTGACATAGGAAAATGTTCAGTGAATACTA 56840300
                  0 |||||||||||||||||.||.|..||.|||.|.|||.|||       40
Query_780       233 AGCGCAAGGCCTGACATGGGTAGCTGCTCAATAAATGCTA      273

""",
            )
        hit = record[5]
        self.assertEqual(hit.num, 6)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|372099098|ref|NC_000078.6|")
        self.assertEqual(hit.target.name, "NC_000078")
        self.assertEqual(
            hit.target.description,
            "Mus musculus strain C57BL/6J chromosome 12, GRCm38.p6 C57BL/6J",
        )
        if xml2:
            self.assertEqual(hit.target.annotations["taxid"], 10090)
            self.assertEqual(hit.target.annotations["sciname"], "Mus musculus")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=120129022)")
        self.assertEqual(len(hit), 2)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 41.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 38.2554)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.07705)
        self.assertEqual(hsp.annotations["identity"], 22)
        if not xml2:
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
        self.assertEqual(hsp.query.id, "Query_78041")
        self.assertEqual(
            hsp.query.description, "G26684.1 human STS STS_D11570, sequence tagged site"
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
            "Mus musculus strain C57BL/6J chromosome 12, GRCm38.p6 C57BL/6J",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(hsp.annotations["midline"], "|||||||||||||||| ||||||")
        if xml2:
            self.assertEqual(
                str(hsp),
                """\
Query : Query_78041 Length: 285 Strand: Plus
        G26684.1 human STS STS_D11570, sequence tagged site
Target: gi|372099098|ref|NC_000078.6| Length: 120129022 Strand: Plus
        Mus musculus strain C57BL/6J chromosome 12, GRCm38.p6 C57BL/6J

Score:38 bits(41), Expect:4,
Identities:22/23(96%),  Gaps:0.23(0%)

gi|372099 113030662 CATCCATTCACACCCAGCACAGG 113030685
                  0 ||||||||||||||||.||||||        23
Query_780        48 CATCCATTCACACCCAACACAGG        71

""",
            )
        else:
            self.assertEqual(
                str(hsp),
                """\
Query : Query_78041 Length: 285 Strand: Plus
        G26684.1 human STS STS_D11570, sequence tagged site
Target: gi|372099098|ref|NC_000078.6| Length: 120129022 Strand: Plus
        Mus musculus strain C57BL/6J chromosome 12, GRCm38.p6 C57BL/6J

Score:38 bits(41), Expect:4,
Identities:22/23(96%),  Positives:22/23(96%),  Gaps:0.23(0%)

gi|372099 113030662 CATCCATTCACACCCAGCACAGG 113030685
                  0 ||||||||||||||||.||||||        23
Query_780        48 CATCCATTCACACCCAACACAGG        71

""",
            )
        hsp = hit[1]
        self.assertEqual(hsp.num, 2)
        self.assertAlmostEqual(hsp.score, 40.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 37.3537)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.07705)
        self.assertEqual(hsp.annotations["identity"], 28)
        if not xml2:
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
        self.assertEqual(hsp.query.id, "Query_78041")
        self.assertEqual(
            hsp.query.description, "G26684.1 human STS STS_D11570, sequence tagged site"
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
            "Mus musculus strain C57BL/6J chromosome 12, GRCm38.p6 C57BL/6J",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(hsp.annotations["midline"], "|||||| | ||||||||||||||| |||| ||")
        if xml2:
            self.assertEqual(
                str(hsp),
                """\
Query : Query_78041 Length: 285 Strand: Plus
        G26684.1 human STS STS_D11570, sequence tagged site
Target: gi|372099098|ref|NC_000078.6| Length: 120129022 Strand: Minus
        Mus musculus strain C57BL/6J chromosome 12, GRCm38.p6 C57BL/6J

Score:37 bits(40), Expect:4,
Identities:28/32(88%),  Gaps:1.32(3%)

gi|372099 108990272 TGTAGCTCTAGGCCTGACATGGGT-GCTGGTC 108990241
                  0 ||||||.|.|||||||||||||||-||||.||        32
Query_780       230 TGTAGCGCAAGGCCTGACATGGGTAGCTGCTC       262

""",
            )
        else:
            self.assertEqual(
                str(hsp),
                """\
Query : Query_78041 Length: 285 Strand: Plus
        G26684.1 human STS STS_D11570, sequence tagged site
Target: gi|372099098|ref|NC_000078.6| Length: 120129022 Strand: Minus
        Mus musculus strain C57BL/6J chromosome 12, GRCm38.p6 C57BL/6J

Score:37 bits(40), Expect:4,
Identities:28/32(88%),  Positives:28/32(88%),  Gaps:1.32(3%)

gi|372099 108990272 TGTAGCTCTAGGCCTGACATGGGT-GCTGGTC 108990241
                  0 ||||||.|.|||||||||||||||-||||.||        32
Query_780       230 TGTAGCGCAAGGCCTGACATGGGTAGCTGCTC       262

""",
            )
        hit = record[6]
        self.assertEqual(hit.num, 7)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|372099049|ref|NT_187008.1|")
        self.assertEqual(hit.target.name, "NT_187008")
        self.assertEqual(
            hit.target.description,
            "Mus musculus strain 129S1/SvImJ chromosome 16 genomic scaffold, GRCm38.p6 alternate locus group 129S1/SvImJ 129S1/SVIMJ_MMCHR16_CTG2",
        )
        if xml2:
            self.assertEqual(hit.target.annotations["taxid"], 10090)
            self.assertEqual(hit.target.annotations["sciname"], "Mus musculus")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=250595)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 40.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 37.3537)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.07705)
        self.assertEqual(hsp.annotations["identity"], 43)
        if not xml2:
            self.assertEqual(hsp.annotations["positive"], 43)
        self.assertEqual(hsp.annotations["gaps"], 2)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[158458, 158482, 158483, 158491, 158492, 158514],
                          [   174,    198,    198,    206,    206,    228]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 56))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({174: 'GGAGGCAAAGAATCCCTACCTCCTAGGGGTGAAAGGAAATNAAAATGGAAAGTT'}, length=285)",
        )
        self.assertEqual(hsp.query.id, "Query_78041")
        self.assertEqual(
            hsp.query.description, "G26684.1 human STS STS_D11570, sequence tagged site"
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({158458: 'GGAGGCAAAGAATCCCTACATTGTGACAGCTGATAAAGAAGGTAAAATGGAAAATT'}, length=250595)",
        )
        self.assertEqual(hsp.target.id, "gi|372099049|ref|NT_187008.1|")
        self.assertEqual(hsp.target.name, "NT_187008")
        self.assertEqual(
            hsp.target.description,
            "Mus musculus strain 129S1/SvImJ chromosome 16 genomic scaffold, GRCm38.p6 alternate locus group 129S1/SvImJ 129S1/SVIMJ_MMCHR16_CTG2",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "||||||||||||||||||| |  | |  | ||| || |||   |||||||||| ||",
        )
        if xml2:
            self.assertEqual(
                str(hsp),
                """\
Query : Query_78041 Length: 285 Strand: Plus
        G26684.1 human STS STS_D11570, sequence tagged site
Target: gi|372099049|ref|NT_187008.1| Length: 250595 Strand: Plus
        Mus musculus strain 129S1/SvImJ chromosome 16 genomic scaffold,
        GRCm38.p6 alternate locus group 129S1/SvImJ 129S1/SVIMJ_MMCHR16_CTG2

Score:37 bits(40), Expect:4,
Identities:43/56(77%),  Gaps:2.56(4%)

gi|372099    158458 GGAGGCAAAGAATCCCTACATTGTGACAGCTGATAAAGAAGGTAAAATGGAAAATT
                  0 |||||||||||||||||||.|..|-|..|.|||-||.|||...||||||||||.||
Query_780       174 GGAGGCAAAGAATCCCTACCTCCT-AGGGGTGA-AAGGAAATNAAAATGGAAAGTT

gi|372099    158514
                 56
Query_780       228

""",
            )
        else:
            self.assertEqual(
                str(hsp),
                """\
Query : Query_78041 Length: 285 Strand: Plus
        G26684.1 human STS STS_D11570, sequence tagged site
Target: gi|372099049|ref|NT_187008.1| Length: 250595 Strand: Plus
        Mus musculus strain 129S1/SvImJ chromosome 16 genomic scaffold,
        GRCm38.p6 alternate locus group 129S1/SvImJ 129S1/SVIMJ_MMCHR16_CTG2

Score:37 bits(40), Expect:4,
Identities:43/56(77%),  Positives:43/56(77%),  Gaps:2.56(4%)

gi|372099    158458 GGAGGCAAAGAATCCCTACATTGTGACAGCTGATAAAGAAGGTAAAATGGAAAATT
                  0 |||||||||||||||||||.|..|-|..|.|||-||.|||...||||||||||.||
Query_780       174 GGAGGCAAAGAATCCCTACCTCCT-AGGGGTGA-AAGGAAATNAAAATGGAAAGTT

gi|372099    158514
                 56
Query_780       228

""",
            )
        hit = record[7]
        self.assertEqual(hit.num, 8)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|372099109|ref|NC_000067.6|")
        self.assertEqual(hit.target.name, "NC_000067")
        self.assertEqual(
            hit.target.description,
            "Mus musculus strain C57BL/6J chromosome 1, GRCm38.p6 C57BL/6J",
        )
        if xml2:
            self.assertEqual(hit.target.annotations["taxid"], 10090)
            self.assertEqual(hit.target.annotations["sciname"], "Mus musculus")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=195471971)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 40.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 37.3537)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.07705)
        self.assertEqual(hsp.annotations["identity"], 35)
        if not xml2:
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
        self.assertEqual(hsp.query.id, "Query_78041")
        self.assertEqual(
            hsp.query.description, "G26684.1 human STS STS_D11570, sequence tagged site"
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
            "Mus musculus strain C57BL/6J chromosome 1, GRCm38.p6 C57BL/6J",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"], "||||||||||| ||||||||| | | | |||||  ||| ||||"
        )
        if xml2:
            self.assertEqual(
                str(hsp),
                """\
Query : Query_78041 Length: 285 Strand: Plus
        G26684.1 human STS STS_D11570, sequence tagged site
Target: gi|372099109|ref|NC_000067.6| Length: 195471971 Strand: Plus
        Mus musculus strain C57BL/6J chromosome 1, GRCm38.p6 C57BL/6J

Score:37 bits(40), Expect:4,
Identities:35/43(81%),  Gaps:2.43(5%)

gi|372099  65190107 GCTCAGCCACATACATGGTTT-TAAGTGTTGAGGCTCT-TTCC 65190148
                  0 |||||||||||.|||||||||-|.|.|.|||||..|||-||||       43
Query_780        86 GCTCAGCCACAGACATGGTTTGTNACTNTTGAGCTTCTGTTCC      129

""",
            )
        else:
            self.assertEqual(
                str(hsp),
                """\
Query : Query_78041 Length: 285 Strand: Plus
        G26684.1 human STS STS_D11570, sequence tagged site
Target: gi|372099109|ref|NC_000067.6| Length: 195471971 Strand: Plus
        Mus musculus strain C57BL/6J chromosome 1, GRCm38.p6 C57BL/6J

Score:37 bits(40), Expect:4,
Identities:35/43(81%),  Positives:35/43(81%),  Gaps:2.43(5%)

gi|372099  65190107 GCTCAGCCACATACATGGTTT-TAAGTGTTGAGGCTCT-TTCC 65190148
                  0 |||||||||||.|||||||||-|.|.|.|||||..|||-||||       43
Query_780        86 GCTCAGCCACAGACATGGTTTGTNACTNTTGAGCTTCTGTTCC      129

""",
            )
        hit = record[8]
        self.assertEqual(hit.num, 9)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|372099101|ref|NC_000075.6|")
        self.assertEqual(hit.target.name, "NC_000075")
        self.assertEqual(
            hit.target.description,
            "Mus musculus strain C57BL/6J chromosome 9, GRCm38.p6 C57BL/6J",
        )
        if xml2:
            self.assertEqual(hit.target.annotations["taxid"], 10090)
            self.assertEqual(hit.target.annotations["sciname"], "Mus musculus")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=124595110)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 40.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 37.3537)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.07705)
        self.assertEqual(hsp.annotations["identity"], 36)
        if not xml2:
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
        self.assertEqual(hsp.query.id, "Query_78041")
        self.assertEqual(
            hsp.query.description, "G26684.1 human STS STS_D11570, sequence tagged site"
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
            "Mus musculus strain C57BL/6J chromosome 9, GRCm38.p6 C57BL/6J",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "||| |||||||| |  |  ||||||||||||| ||| | | || |||",
        )
        if xml2:
            self.assertEqual(
                str(hsp),
                """\
Query : Query_78041 Length: 285 Strand: Plus
        G26684.1 human STS STS_D11570, sequence tagged site
Target: gi|372099101|ref|NC_000075.6| Length: 124595110 Strand: Minus
        Mus musculus strain C57BL/6J chromosome 9, GRCm38.p6 C57BL/6J

Score:37 bits(40), Expect:4,
Identities:36/47(77%),  Gaps:0.47(0%)

gi|372099  58227241 CAAAGCCTGACAGGTATGACTGCTCAATAAATACTATTTTTTTTTTT 58227194
                  0 |||.||||||||.|..|..|||||||||||||.|||.|.|.||.|||       47
Query_780       237 CAAGGCCTGACATGGGTAGCTGCTCAATAAATGCTAGTNTGTTATTT      284

""",
            )
        else:
            self.assertEqual(
                str(hsp),
                """\
Query : Query_78041 Length: 285 Strand: Plus
        G26684.1 human STS STS_D11570, sequence tagged site
Target: gi|372099101|ref|NC_000075.6| Length: 124595110 Strand: Minus
        Mus musculus strain C57BL/6J chromosome 9, GRCm38.p6 C57BL/6J

Score:37 bits(40), Expect:4,
Identities:36/47(77%),  Positives:36/47(77%),  Gaps:0.47(0%)

gi|372099  58227241 CAAAGCCTGACAGGTATGACTGCTCAATAAATACTATTTTTTTTTTT 58227194
                  0 |||.||||||||.|..|..|||||||||||||.|||.|.|.||.|||       47
Query_780       237 CAAGGCCTGACATGGGTAGCTGCTCAATAAATGCTAGTNTGTTATTT      284

""",
            )
        hit = record[9]
        self.assertEqual(hit.num, 10)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|372099100|ref|NC_000076.6|")
        self.assertEqual(hit.target.name, "NC_000076")
        self.assertEqual(
            hit.target.description,
            "Mus musculus strain C57BL/6J chromosome 10, GRCm38.p6 C57BL/6J",
        )
        if xml2:
            self.assertEqual(hit.target.annotations["taxid"], 10090)
            self.assertEqual(hit.target.annotations["sciname"], "Mus musculus")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=130694993)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 40.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 37.3537)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.07705)
        self.assertEqual(hsp.annotations["identity"], 20)
        if not xml2:
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
        self.assertEqual(hsp.query.id, "Query_78041")
        self.assertEqual(
            hsp.query.description, "G26684.1 human STS STS_D11570, sequence tagged site"
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
            "Mus musculus strain C57BL/6J chromosome 10, GRCm38.p6 C57BL/6J",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(hsp.annotations["midline"], "||||||||||||||||||||")
        if xml2:
            self.assertEqual(
                str(hsp),
                """\
Query : Query_78041 Length: 285 Strand: Plus
        G26684.1 human STS STS_D11570, sequence tagged site
Target: gi|372099100|ref|NC_000076.6| Length: 130694993 Strand: Plus
        Mus musculus strain C57BL/6J chromosome 10, GRCm38.p6 C57BL/6J

Score:37 bits(40), Expect:4,
Identities:20/20(100%),  Gaps:0.20(0%)

gi|372099 119337185 AGCTGCTCAATAAATGCTAG 119337205
                  0 ||||||||||||||||||||        20
Query_780       254 AGCTGCTCAATAAATGCTAG       274

""",
            )
        else:
            self.assertEqual(
                str(hsp),
                """\
Query : Query_78041 Length: 285 Strand: Plus
        G26684.1 human STS STS_D11570, sequence tagged site
Target: gi|372099100|ref|NC_000076.6| Length: 130694993 Strand: Plus
        Mus musculus strain C57BL/6J chromosome 10, GRCm38.p6 C57BL/6J

Score:37 bits(40), Expect:4,
Identities:20/20(100%),  Positives:20/20(100%),  Gaps:0.20(0%)

gi|372099 119337185 AGCTGCTCAATAAATGCTAG 119337205
                  0 ||||||||||||||||||||        20
Query_780       254 AGCTGCTCAATAAATGCTAG       274

""",
            )
        hit = record[10]
        self.assertEqual(hit.num, 11)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|372099094|ref|NC_000082.6|")
        self.assertEqual(hit.target.name, "NC_000082")
        self.assertEqual(
            hit.target.description,
            "Mus musculus strain C57BL/6J chromosome 16, GRCm38.p6 C57BL/6J",
        )
        if xml2:
            self.assertEqual(hit.target.annotations["taxid"], 10090)
            self.assertEqual(hit.target.annotations["sciname"], "Mus musculus")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=98207768)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 40.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 37.3537)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.07705)
        self.assertEqual(hsp.annotations["identity"], 43)
        if not xml2:
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
        self.assertEqual(hsp.query.id, "Query_78041")
        self.assertEqual(
            hsp.query.description, "G26684.1 human STS STS_D11570, sequence tagged site"
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
            "Mus musculus strain C57BL/6J chromosome 16, GRCm38.p6 C57BL/6J",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "||||||||||||||||||| |  | |  | ||| || |||   |||||||||| ||",
        )
        if xml2:
            self.assertEqual(
                str(hsp),
                """\
Query : Query_78041 Length: 285 Strand: Plus
        G26684.1 human STS STS_D11570, sequence tagged site
Target: gi|372099094|ref|NC_000082.6| Length: 98207768 Strand: Plus
        Mus musculus strain C57BL/6J chromosome 16, GRCm38.p6 C57BL/6J

Score:37 bits(40), Expect:4,
Identities:43/56(77%),  Gaps:2.56(4%)

gi|372099  18854779 GGAGGCAAAGAATCCCTACATTGTGACAGCTGATAAAGAAGGTAAAATGGAAAATT
                  0 |||||||||||||||||||.|..|-|..|.|||-||.|||...||||||||||.||
Query_780       174 GGAGGCAAAGAATCCCTACCTCCT-AGGGGTGA-AAGGAAATNAAAATGGAAAGTT

gi|372099  18854835
                 56
Query_780       228

""",
            )
        else:
            self.assertEqual(
                str(hsp),
                """\
Query : Query_78041 Length: 285 Strand: Plus
        G26684.1 human STS STS_D11570, sequence tagged site
Target: gi|372099094|ref|NC_000082.6| Length: 98207768 Strand: Plus
        Mus musculus strain C57BL/6J chromosome 16, GRCm38.p6 C57BL/6J

Score:37 bits(40), Expect:4,
Identities:43/56(77%),  Positives:43/56(77%),  Gaps:2.56(4%)

gi|372099  18854779 GGAGGCAAAGAATCCCTACATTGTGACAGCTGATAAAGAAGGTAAAATGGAAAATT
                  0 |||||||||||||||||||.|..|-|..|.|||-||.|||...||||||||||.||
Query_780       174 GGAGGCAAAGAATCCCTACCTCCT-AGGGGTGA-AAGGAAATNAAAATGGAAAGTT

gi|372099  18854835
                 56
Query_780       228

""",
            )

    def test_xml_21500_blastn_001_writer(self):
        """Writing BLASTN 2.15.0+ (xml_21500_blastn_001.xml)."""
        filename = "xml_21500_blastn_001.xml"
        path = os.path.join("Blast", filename)
        with Blast.parse(path) as records:
            stream = io.BytesIO()
            n = Blast.write(records, stream)
            self.assertEqual(n, 1)
            stream.seek(0)
            written_records = Blast.parse(stream)
            self.check_xml_21500_blastn_001_records(written_records)

    def test_xml2_21500_blastn_001_writer(self):
        """Writing BLASTN 2.15.0+ XML2 (xml2_21500_blastn_001.xml)."""
        filename = "xml2_21500_blastn_001.xml"
        path = os.path.join("Blast", filename)
        with Blast.parse(path) as records:
            stream = io.BytesIO()
            n = Blast.write(records, stream, fmt="XML2")
            self.assertEqual(n, 1)
            stream.seek(0)
            written_records = Blast.parse(stream)
            self.check_xml_21500_blastn_001_records(written_records, xml2=True)

    def test_megablast_legacy(self):
        """Parsing megablast 2.2.26 [Sep-21-2011] (megablast_legacy.xml)."""
        filename = "megablast_legacy.xml"
        path = os.path.join("Blast", filename)
        with open(path, "rb") as stream:
            records = Blast.parse(stream)
            self.check_megablast_legacy_records(records)
        with Blast.parse(path) as records:
            self.check_megablast_legacy_records(records)
        with open(path, "rb") as stream:
            record = Blast.read(stream)
        self.check_megablast_legacy_record(record)
        record = Blast.read(path)
        self.check_megablast_legacy_record(record)
        self.assertEqual(
            str(record[1::2]),
            """\
Program: megablast 2.2.26 [Sep-21-2011]
     db: m_cold.fasta
  Query: lcl|1_ (length=1111)
         gi|8332116|gb|BE037100.1|BE037100 MP14H09 MP Mesembryanthemum
         crystallinum cDNA 5' similar to cold acclimation protein, mRNA sequence
   Hits: No hits found""",
        )

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
        self.assertEqual(repr(record), "<Bio.Blast.Record query.id='lcl|1_'; 1 hit>")
        self.assertEqual(len(record), 1)
        hit = record[0]
        self.assertEqual(hit.num, 1)
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
        self.assertEqual(hsp.num, 1)
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
Query : lcl|1_ Length: 1111 Strand: Plus
        gi|8332116|gb|BE037100.1|BE037100 MP14H09 MP Mesembryanthemum
        crystallinum cDNA 5' similar to cold acclimation protein, mRNA sequence
Target: gnl|BL_ORD_ID|0 Length: 1111 Strand: Plus
        gi|8332116|gb|BE037100.1|BE037100 MP14H09 MP Mesembryanthemum
        crystallinum cDNA 5' similar to cold acclimation protein, mRNA sequence

Score:1562 bits(788), Expect:0,
Identities:797/797(100%),  Positives:797/797(100%)

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

    def test_megablast_legacy_writer(self):
        """Writing megablast 2.2.26 [Sep-21-2011] (megablast_legacy.xml)."""
        filename = "megablast_legacy.xml"
        path = os.path.join("Blast", filename)
        with Blast.parse(path) as records:
            stream = io.BytesIO()
            n = Blast.write(records, stream)
            self.assertEqual(n, 1)
            stream.seek(0)
            written_records = Blast.parse(stream)
            self.check_megablast_legacy_records(written_records)


class TestBlastx(unittest.TestCase):
    """Test the Blast XML parser for blastx output."""

    def test_xml_2222_blastx_001_parser(self):
        """Parsing BLASTX 2.2.22+ (xml_2222_blastx_001.xml)."""
        filename = "xml_2222_blastx_001.xml"
        path = os.path.join("Blast", filename)
        with open(path, "rb") as stream:
            records = Blast.parse(stream)
            self.check_xml_2222_blastx_001(records)
        with open(path, "rb") as stream:
            records = Blast.parse(stream)
            self.assertEqual(
                str(records),
                """\
Program: BLASTX 2.2.22+
     db: nr

  Query: 1 (length=1002)
         gi|4104054|gb|AH007193.1|SEG_CVIGS Centaurea vallesiaca 18S ribosomal
         RNA gene, partial sequence
   Hits: ----  -----  ----------------------------------------------------------
            #  # HSP  ID + description
         ----  -----  ----------------------------------------------------------
            0      1  gi|149390769|gb|ABR25402.1|  unknown [Oryza sativa (ind...

  Query: 2 (length=2050)
         gi|4218935|gb|AF074388.1|AF074388 Sambucus nigra hevein-like protein
         HLPf gene, partial cds
   Hits: ----  -----  ----------------------------------------------------------
            #  # HSP  ID + description
         ----  -----  ----------------------------------------------------------
            0      2  gi|4218936|gb|AAD12237.1|  hevein-like protein HLPf [Sa...
            1      2  gi|4206074|gb|AAD11408.1|  hevein-like protein [Sambucu...
            2      2  gi|4206070|gb|AAD11406.1|  hevein-like protein [Sambucu...
            3      2  gi|4206072|gb|AAD11407.1|  hevein-like protein [Sambucu...
            4      2  gi|16903131|gb|AAL30421.1|AF434174_1  hevein-like prote...
            5      2  gi|16903133|gb|AAL30422.1|AF434175_1  hevein-like prote...
            6      2  gi|30691147|gb|AAO17294.1|  chitinase [Ficus carica]
            7      2  gi|222139388|gb|ACM45713.1|  class I chitinase [Pyrus p...
            8      2  gi|23496435|dbj|BAB40817.2|  endochitinase MCHT-2 [Cucu...
            9      2  gi|82621253|gb|ABB86300.1|  chitinase [Ficus awkeotsang]

  Query: 3 (length=550)
         gi|5690369|gb|AF158246.1|AF158246 Cricetulus griseus glucose phosphate
         isomerase (GPI) gene, partial intron sequence
   Hits: No hits found

  Query: 4 (length=655)
         gi|5049839|gb|AI730987.1|AI730987 BNLGHi8354 Six-day Cotton fiber
         Gossypium hirsutum cDNA 5' similar to TUBULIN BETA-1 CHAIN
         gi|486734|pir|S35142 tubulin beta chain - white lupine gi|402636
         (X70184) Beta tubulin 1 [Lupinus albus], mRNA sequence
   Hits: ----  -----  ----------------------------------------------------------
            #  # HSP  ID + description
         ----  -----  ----------------------------------------------------------
            0      1  gi|166343825|gb|ABY86655.1|  beta-tubulin 4 [Gossypium ...
            1      1  gi|223549899|gb|EEF51386.1|  tubulin beta chain, putati...
            2      1  gi|18420724|ref|NP_568437.1|  TUB8 (tubulin beta-8) [Ar...
            3      1  gi|225426385|ref|XP_002271992.1|  PREDICTED: hypothetic...
            4      1  gi|15451226|gb|AAK96884.1|  beta tubulin [Arabidopsis t...
            5      1  gi|225470745|ref|XP_002267380.1|  PREDICTED: hypothetic...
            6      1  gi|586076|sp|P37392.1|TBB1_LUPAL  RecName: Full=Tubulin...
            7      1  gi|224104341|ref|XP_002313404.1|  tubulin, beta chain [...
            8      1  gi|223549679|gb|EEF51167.1|  tubulin beta chain, putati...
            9      1  gi|224058553|ref|XP_002299541.1|  tubulin, beta chain [...

  Query: 5 (length=623)
         gi|5052071|gb|AF067555.1|AF067555 Phlox stansburyi internal transcribed
         spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2,
         complete sequence
   Hits: ----  -----  ----------------------------------------------------------
            #  # HSP  ID + description
         ----  -----  ----------------------------------------------------------
            0      2  gi|110740644|dbj|BAE98425.1|  hypothetical protein [Ara...
            1      1  gi|226453533|gb|EEH50844.1|  predicted protein [Micromo...
            2      1  gi|168069582|ref|XP_001786502.1|  predicted protein [Ph...
            3      1  gi|168068558|ref|XP_001786120.1|  predicted protein [Ph...
            4      1  gi|168068926|ref|XP_001786259.1|  predicted protein [Ph...
            5      1  gi|168070288|ref|XP_001786759.1|  predicted protein [Ph...
            6      1  gi|168068591|ref|XP_001786133.1|  predicted protein [Ph...
            7      1  gi|74622391|sp|Q8TGM5|ART3_YEAST  Uncharacterized prote...
            8      1  gi|168069944|ref|XP_001786634.1|  predicted protein [Ph...
            9      1  gi|50307717|ref|XP_453851.1|  unnamed protein product [...

  Query: 6 (length=309)
         gi|3176602|gb|U78617.1|LOU78617 Lathyrus odoratus phytochrome A (PHYA)
         gene, partial cds
   Hits: ----  -----  ----------------------------------------------------------
            #  # HSP  ID + description
         ----  -----  ----------------------------------------------------------
            0      1  gi|3176603|gb|AAC18749.1|  phytochrome A [Lathyrus odor...
            1      1  gi|130188|sp|P15001.1|PHYA_PEA  RecName: Full=Phytochro...
            2      1  gi|2499555|sp|P93673.1|PHYA_LATSA  RecName: Full=Phytoc...
            3      1  gi|3176595|gb|AAC18745.1|  phytochrome A [Lennea melano...
            4      1  gi|1711106|gb|AAC18675.1|  phytochrome A [Sophora affinis]
            5      1  gi|1711090|gb|AAC18670.1|  phytochrome A [Myrospermum s...
            6      1  gi|3176605|gb|AAC18750.1|  phytochrome A [Hybosema robu...
            7      1  gi|3176454|gb|AAC18668.1|  phytochrome A [Cyclolobium n...
            8      1  gi|3176523|gb|AAC18709.1|  phytochrome A [Millettia ric...
            9      1  gi|3176494|gb|AAC18693.1|  phytochrome A [Callerya atro...

  Query: 7 (length=2551)
         gi|5817701|gb|AF142731.1|AF142731 Wisteria frutescens maturase-like
         protein (matK) gene, complete cds; chloroplast gene for chloroplast
         product
   Hits: ----  -----  ----------------------------------------------------------
            #  # HSP  ID + description
         ----  -----  ----------------------------------------------------------
            0      1  gi|27805603|sp|Q9TKP6.1|MATK_WISFR  RecName: Full=Matur...
            1      1  gi|171909144|gb|ACB58148.1|  maturase K [Wisteria frute...
            2      1  gi|171909146|gb|ACB58149.1|  maturase K [Wisteria frute...
            3      1  gi|171909132|gb|ACB58142.1|  maturase K [Callerya megas...
            4      1  gi|5817760|gb|AAD52903.1|AF142732_1  maturase-like prot...
            5      1  gi|171909134|gb|ACB58143.1|  maturase K [Wisteria brach...
            6      1  gi|5817761|gb|AAD52904.1|AF142733_1  maturase-like prot...
            7      1  gi|5817762|gb|AAD52905.1|AF142734_1  maturase-like prot...
            8      1  gi|152014012|gb|ABS20107.1|  maturase-like protein [Ast...
            9      1  gi|146197442|dbj|BAF57483.1|  maturase [Glycyrrhiza ura...""",
            )

    def check_xml_2222_blastx_001(self, records):
        self.assertEqual(records.program, "blastx")
        self.assertEqual(records.version, "BLASTX 2.2.22+")
        self.assertEqual(
            records.reference,
            'Stephen F. Altschul, Thomas L. Madden, Alejandro A. Schäffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.',
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
        self.assertEqual(record.num, 1)
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
        self.assertEqual(record.stat["eff-space"], 367397307882)
        self.assertAlmostEqual(record.stat["kappa"], 0.041)
        self.assertAlmostEqual(record.stat["lambda"], 0.267)
        self.assertAlmostEqual(record.stat["entropy"], 0.14)
        self.assertEqual(len(record), 1)
        hit = record[0]
        self.assertEqual(hit.num, 1)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|149390769|gb|ABR25402.1|")
        self.assertEqual(hit.target.name, "ABR25402")
        self.assertEqual(
            hit.target.description, "unknown [Oryza sativa (indica cultivar-group)]"
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=26)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 129.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 54.2989775733826)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 1.83262460293058e-05, places=19
        )
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
Query : 1 Length: 26 Strand: Plus
        gi|4104054|gb|AH007193.1|SEG_CVIGS Centaurea vallesiaca 18S ribosomal
        RNA gene, partial sequence
Target: gi|149390769|gb|ABR25402.1| Length: 26 Strand: Plus
        unknown [Oryza sativa (indica cultivar-group)]

Score:54 bits(129), Expect:2e-05,
Identities:24/26(92%),  Positives:25/26(96%),  Gaps:0.26(0%)

gi|149390         0 HMLVSKIKPCMCKYELIRTVKLRMAH 26
                  0 |||||||||||||||.|.|||||||| 26
1                 0 HMLVSKIKPCMCKYEQIQTVKLRMAH 26

""",
        )
        record = next(records)
        self.assertEqual(record.num, 2)
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
        self.assertEqual(record.stat["eff-space"], 967993058520)
        self.assertAlmostEqual(record.stat["kappa"], 0.041)
        self.assertAlmostEqual(record.stat["lambda"], 0.267)
        self.assertAlmostEqual(record.stat["entropy"], 0.14)
        self.assertEqual(len(record), 10)
        hit = record[0]
        self.assertEqual(hit.num, 1)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|4218936|gb|AAD12237.1|")
        self.assertEqual(hit.target.name, "AAD12237")
        self.assertEqual(
            hit.target.description, "hevein-like protein HLPf [Sambucus nigra]"
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=333)")
        self.assertEqual(len(hit), 2)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 1053.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 410.223385721017)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 3.48406066731465e-112, places=126
        )
        self.assertEqual(hsp.annotations["identity"], 199)
        self.assertEqual(hsp.annotations["positive"], 200)
        self.assertEqual(hsp.annotations["gaps"], 33)
        hsp = hit[1]
        self.assertEqual(hsp.num, 2)
        self.assertAlmostEqual(hsp.score, 683.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 267.699542631596)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 2.79278546744412e-69, places=83
        )
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
Query : 2 Length: 127 Strand: Plus
        gi|4218935|gb|AF074388.1|AF074388 Sambucus nigra hevein-like protein
        HLPf gene, partial cds
Target: gi|4218936|gb|AAD12237.1| Length: 333 Strand: Plus
        hevein-like protein HLPf [Sambucus nigra]

Score:267 bits(683), Expect:3e-69,
Identities:127/127(100%),  Positives:127/127(100%),  Gaps:0.127(0%)

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
        self.assertEqual(hit.num, 2)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|4206074|gb|AAD11408.1|")
        self.assertEqual(hit.target.name, "AAD11408")
        self.assertEqual(hit.target.description, "hevein-like protein [Sambucus nigra]")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=333)")
        self.assertEqual(len(hit), 2)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 1043.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 406.371389961843)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 5.03097287018806e-111, places=125
        )
        self.assertEqual(hsp.annotations["identity"], 198)
        self.assertEqual(hsp.annotations["positive"], 199)
        self.assertEqual(hsp.annotations["gaps"], 33)
        hsp = hit[1]
        self.assertEqual(hsp.num, 2)
        self.assertAlmostEqual(hsp.score, 672.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 263.462347296505)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 5.26696544712228e-68, places=82
        )
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
        self.assertEqual(hsp.target.description, "hevein-like protein [Sambucus nigra]")
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "NY YGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINANSEASDQVPSYGVISKIINSNFGHQSGLDTITTSIGYYKRYCDMLEVSYGDNLENWFDETPFTKVAHIKMS+V",
        )
        self.assertEqual(
            str(hsp),
            """\
Query : 2 Length: 127 Strand: Plus
        gi|4218935|gb|AF074388.1|AF074388 Sambucus nigra hevein-like protein
        HLPf gene, partial cds
Target: gi|4206074|gb|AAD11408.1| Length: 333 Strand: Plus
        hevein-like protein [Sambucus nigra]

Score:263 bits(672), Expect:5e-68,
Identities:125/127(98%),  Positives:126/127(99%),  Gaps:0.127(0%)

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
        self.assertEqual(hit.num, 3)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|4206070|gb|AAD11406.1|")
        self.assertEqual(hit.target.name, "AAD11406")
        self.assertEqual(hit.target.description, "hevein-like protein [Sambucus nigra]")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=333)")
        self.assertEqual(len(hit), 2)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 1043.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 406.371389961843)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 5.03097287018806e-111, places=125
        )
        self.assertEqual(hsp.annotations["identity"], 198)
        self.assertEqual(hsp.annotations["positive"], 199)
        self.assertEqual(hsp.annotations["gaps"], 33)
        hsp = hit[1]
        self.assertEqual(hsp.num, 2)
        self.assertAlmostEqual(hsp.score, 680.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 266.543943903844)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 6.22167692942359e-69, places=83
        )
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
        self.assertEqual(hsp.target.description, "hevein-like protein [Sambucus nigra]")
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "NYNYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINANSEASDQVPSYGVISKIINSNFGHQSGLDTITTSIGYYKRYCDMLEVSYGDNLENWFDETPFTKVAHIKMS+V",
        )
        self.assertEqual(
            str(hsp),
            """\
Query : 2 Length: 127 Strand: Plus
        gi|4218935|gb|AF074388.1|AF074388 Sambucus nigra hevein-like protein
        HLPf gene, partial cds
Target: gi|4206070|gb|AAD11406.1| Length: 333 Strand: Plus
        hevein-like protein [Sambucus nigra]

Score:266 bits(680), Expect:6e-69,
Identities:126/127(99%),  Positives:127/127(100%),  Gaps:0.127(0%)

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
        self.assertEqual(hit.num, 4)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|4206072|gb|AAD11407.1|")
        self.assertEqual(hit.target.name, "AAD11407")
        self.assertEqual(hit.target.description, "hevein-like protein [Sambucus nigra]")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=333)")
        self.assertEqual(len(hit), 2)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 1016.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 395.971001412075)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 6.7995613312017e-108, places=121
        )
        self.assertEqual(hsp.annotations["identity"], 193)
        self.assertEqual(hsp.annotations["positive"], 195)
        self.assertEqual(hsp.annotations["gaps"], 33)
        hsp = hit[1]
        self.assertEqual(hsp.num, 2)
        self.assertAlmostEqual(hsp.score, 646.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 253.447158322654)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 5.45045505347399e-65, places=79
        )
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
        self.assertEqual(hsp.target.description, "hevein-like protein [Sambucus nigra]")
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "NY YGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINANSEASDQVPSYGVIS+II+SN GHQS LDTITTSIGYYKRYCDMLEVSYGDNL+NWFDETPF+KVAHIKMSVV",
        )
        self.assertEqual(
            str(hsp),
            """\
Query : 2 Length: 127 Strand: Plus
        gi|4218935|gb|AF074388.1|AF074388 Sambucus nigra hevein-like protein
        HLPf gene, partial cds
Target: gi|4206072|gb|AAD11407.1| Length: 333 Strand: Plus
        hevein-like protein [Sambucus nigra]

Score:253 bits(646), Expect:5e-65,
Identities:120/127(94%),  Positives:124/127(98%),  Gaps:0.127(0%)

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
        self.assertEqual(hit.num, 5)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|16903131|gb|AAL30421.1|AF434174_1")
        self.assertEqual(hit.target.name, "AAL30421")
        self.assertEqual(hit.target.description, "hevein-like protein [Sambucus nigra]")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=330)")
        self.assertEqual(len(hit), 2)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 986.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 384.415014134554)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 2.04729155722083e-104, places=118
        )
        self.assertEqual(hsp.annotations["identity"], 190)
        self.assertEqual(hsp.annotations["positive"], 191)
        self.assertEqual(hsp.annotations["gaps"], 36)
        hsp = hit[1]
        self.assertEqual(hsp.num, 2)
        self.assertAlmostEqual(hsp.score, 679.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 266.158744327927)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 8.12576171382949e-69, places=83
        )
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
        self.assertEqual(hsp.target.description, "hevein-like protein [Sambucus nigra]")
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "NYNYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINANSEASDQVPSYGVISKIINSNFGHQSGLDTITTSIGYYKRYCDMLEVSYGDNLENWFDETPFTKV HIKMSVV",
        )
        self.assertEqual(
            str(hsp),
            """\
Query : 2 Length: 127 Strand: Plus
        gi|4218935|gb|AF074388.1|AF074388 Sambucus nigra hevein-like protein
        HLPf gene, partial cds
Target: gi|16903131|gb|AAL30421.1|AF434174_1 Length: 330 Strand: Plus
        hevein-like protein [Sambucus nigra]

Score:266 bits(679), Expect:8e-69,
Identities:126/127(99%),  Positives:126/127(99%),  Gaps:0.127(0%)

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
        self.assertEqual(hit.num, 6)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|16903133|gb|AAL30422.1|AF434175_1")
        self.assertEqual(hit.target.name, "AAL30422")
        self.assertEqual(hit.target.description, "hevein-like protein [Sambucus nigra]")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=336)")
        self.assertEqual(len(hit), 2)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 713.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 279.255529909117)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 9.27553088557319e-73, places=87
        )
        self.assertEqual(hsp.annotations["identity"], 148)
        self.assertEqual(hsp.annotations["positive"], 162)
        self.assertEqual(hsp.annotations["gaps"], 40)
        hsp = hit[1]
        self.assertEqual(hsp.num, 2)
        self.assertAlmostEqual(hsp.score, 620.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 243.431969348803)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 5.64033703812707e-62, places=76
        )
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
        self.assertEqual(hsp.target.description, "hevein-like protein [Sambucus nigra]")
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "NYNYGLAGEA+GIDLVN PDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINAN EASDQVPSYGV+S IINSN GH+SGLD ITTSIGYYKRYCDMLEVSYGDNL+NWFDETPF+KVA IKM VV",
        )
        self.assertEqual(
            str(hsp),
            """\
Query : 2 Length: 127 Strand: Plus
        gi|4218935|gb|AF074388.1|AF074388 Sambucus nigra hevein-like protein
        HLPf gene, partial cds
Target: gi|16903133|gb|AAL30422.1|AF434175_1 Length: 336 Strand: Plus
        hevein-like protein [Sambucus nigra]

Score:243 bits(620), Expect:6e-62,
Identities:115/127(91%),  Positives:120/127(94%),  Gaps:0.127(0%)

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
        self.assertEqual(hit.num, 7)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|30691147|gb|AAO17294.1|")
        self.assertEqual(hit.target.name, "AAO17294")
        self.assertEqual(hit.target.description, "chitinase [Ficus carica]")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=321)")
        self.assertEqual(len(hit), 2)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 481.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 189.889228296291)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 7.40075731140555e-46, places=60
        )
        self.assertEqual(hsp.annotations["identity"], 113)
        self.assertEqual(hsp.annotations["positive"], 138)
        self.assertEqual(hsp.annotations["gaps"], 49)
        hsp = hit[1]
        self.assertEqual(hsp.num, 2)
        self.assertAlmostEqual(hsp.score, 426.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 168.703251620836)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 1.76559312729927e-39, places=53
        )
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
Query : 2 Length: 126 Strand: Plus
        gi|4218935|gb|AF074388.1|AF074388 Sambucus nigra hevein-like protein
        HLPf gene, partial cds
Target: gi|30691147|gb|AAO17294.1| Length: 321 Strand: Plus
        chitinase [Ficus carica]

Score:168 bits(426), Expect:2e-39,
Identities:81/127(64%),  Positives:99/127(78%),  Gaps:10.127(8%)

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
        self.assertEqual(hit.num, 8)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|222139388|gb|ACM45713.1|")
        self.assertEqual(hit.target.name, "ACM45713")
        self.assertEqual(hit.target.description, "class I chitinase [Pyrus pyrifolia]")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=317)")
        self.assertEqual(len(hit), 2)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 469.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 185.266833385283)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 1.82286993845418e-44, places=58
        )
        self.assertEqual(hsp.annotations["identity"], 111)
        self.assertEqual(hsp.annotations["positive"], 137)
        self.assertEqual(hsp.annotations["gaps"], 50)
        hsp = hit[1]
        self.assertEqual(hsp.num, 2)
        self.assertAlmostEqual(hsp.score, 318.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 127.101697421762)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 5.89123449548921e-27, places=40
        )
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
        self.assertEqual(hsp.target.description, "class I chitinase [Pyrus pyrifolia]")
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "NYNYG AG+A+G DL+N+PDLVATDP+VSFKTAIWFWMT   N PS HD++      +    ++ +VP YGVI+ IIN       G D  + + IG+Y+RYC +L V+ GDNL+  +++ PF",
        )
        self.assertEqual(
            str(hsp),
            """\
Query : 2 Length: 116 Strand: Plus
        gi|4218935|gb|AF074388.1|AF074388 Sambucus nigra hevein-like protein
        HLPf gene, partial cds
Target: gi|222139388|gb|ACM45713.1| Length: 317 Strand: Plus
        class I chitinase [Pyrus pyrifolia]

Score:127 bits(318), Expect:6e-27,
Identities:62/122(51%),  Positives:84/122(69%),  Gaps:7.122(6%)

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
        self.assertEqual(hit.num, 9)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|23496435|dbj|BAB40817.2|")
        self.assertEqual(hit.target.name, "BAB40817")
        self.assertEqual(hit.target.description, "endochitinase MCHT-2 [Cucumis melo]")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=311)")
        self.assertEqual(len(hit), 2)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 460.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 181.800037202027)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 2.01541888137674e-43, places=57
        )
        self.assertEqual(hsp.annotations["identity"], 109)
        self.assertEqual(hsp.annotations["positive"], 132)
        self.assertEqual(hsp.annotations["gaps"], 54)
        hsp = hit[1]
        self.assertEqual(hsp.num, 2)
        self.assertAlmostEqual(hsp.score, 285.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 114.39011141649)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 3.95161831690075e-23, places=37
        )
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
        self.assertEqual(hsp.target.description, "endochitinase MCHT-2 [Cucumis melo]")
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "NYNYG AG+A+G  L+ +PD  ATDP+ SFKTA+WFWMT   N PS H+++      ++   A+ +VP YGVI+ IIN       G  D +   IG+YKRYCDML + YG+NL+",
        )
        self.assertEqual(
            str(hsp),
            """\
Query : 2 Length: 108 Strand: Plus
        gi|4218935|gb|AF074388.1|AF074388 Sambucus nigra hevein-like protein
        HLPf gene, partial cds
Target: gi|23496435|dbj|BAB40817.2| Length: 311 Strand: Plus
        endochitinase MCHT-2 [Cucumis melo]

Score:114 bits(285), Expect:4e-23,
Identities:56/114(49%),  Positives:75/114(66%),  Gaps:7.114(6%)

gi|234964       191 NYNYGPAGKAIGAPLLTNPD-TATDPVTSFKTALWFWMTAQGNKPSCHNVITGNWQPSSA
                  0 |||||.||.|.|..|...||-.||||..|||||.|||||...|.||.|....-----...
2                 0 NYNYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILI-----NAN

gi|234964       250 DNAAGRVPGYGVITNIINGGLECGRGPDDRVKDRIGFYKRYCDMLGIGYGNNLD 304
                 60 ..|...||.||||..|||.......|.-|.....||.||||||||...||.||. 114
2                55 SEASDQVPSYGVISKIINSNFGHQSGL-DTITTSIGYYKRYCDMLEVSYGDNLE 108

""",
        )
        hit = record[9]
        self.assertEqual(hit.num, 10)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|82621253|gb|ABB86300.1|")
        self.assertEqual(hit.target.name, "ABB86300")
        self.assertEqual(hit.target.description, "chitinase [Ficus awkeotsang]")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=301)")
        self.assertEqual(len(hit), 2)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 459.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 181.414837626109)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 2.6322185753765e-43, places=56
        )
        self.assertEqual(hsp.annotations["identity"], 114)
        self.assertEqual(hsp.annotations["positive"], 134)
        self.assertEqual(hsp.annotations["gaps"], 50)
        hsp = hit[1]
        self.assertEqual(hsp.num, 2)
        self.assertAlmostEqual(hsp.score, 359.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 142.894880034374)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 1.03749166509001e-31, places=45
        )
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
Query : 2 Length: 107 Strand: Plus
        gi|4218935|gb|AF074388.1|AF074388 Sambucus nigra hevein-like protein
        HLPf gene, partial cds
Target: gi|82621253|gb|ABB86300.1| Length: 301 Strand: Plus
        chitinase [Ficus awkeotsang]

Score:142 bits(359), Expect:1e-31,
Identities:67/107(63%),  Positives:83/107(78%),  Gaps:9.107(8%)

gi|826212       203 NHNYGLVGEALGIDLINNPELVATDPVISFKTAIWFWMARYEDKPSFHDVIINANFEASD
                  0 |.||||.||||||||.|.|.||||||..||||||||||......||.||..||||.||||
2                 0 NYNYGLAGEALGIDLVNHPDLVATDPIVSFKTAIWFWMTQHDNNPSLHDILINANSEASD

gi|826212       263 -IPYHG--------NSGQESSLDVVNRSIGYYKRYCDMLGVSCEDNL 301
                 60 -.|..|--------|.|..|.||....||||||||||||.||..||| 107
2                60 QVPSYGVISKIINSNFGHQSGLDTITTSIGYYKRYCDMLEVSYGDNL 107

""",
        )
        record = next(records)
        self.assertEqual(record.num, 3)
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
        self.assertEqual(record.stat["eff-space"], 108443629616)
        self.assertAlmostEqual(record.stat["kappa"], 0.041)
        self.assertAlmostEqual(record.stat["lambda"], 0.267)
        self.assertAlmostEqual(record.stat["entropy"], 0.14)
        self.assertEqual(len(record), 0)
        record = next(records)
        self.assertEqual(record.num, 4)
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
        self.assertAlmostEqual(record.stat["eff-space"], 165344802738)
        self.assertAlmostEqual(record.stat["kappa"], 0.041)
        self.assertAlmostEqual(record.stat["lambda"], 0.267)
        self.assertAlmostEqual(record.stat["entropy"], 0.14)
        self.assertEqual(len(record), 10)
        hit = record[0]
        self.assertEqual(hit.num, 1)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|166343825|gb|ABY86655.1|")
        self.assertEqual(hit.target.name, "ABY86655")
        self.assertEqual(hit.target.description, "beta-tubulin 4 [Gossypium hirsutum]")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=448)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 1048.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 408.29738784143)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 2.26145081918239e-112, places=126
        )
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
        self.assertEqual(hsp.target.description, "beta-tubulin 4 [Gossypium hirsutum]")
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MREILHIQAGQCGNQIGA FWEVVCAEHGI+STGRYQGDNDLQLERVNVYYNEASCGRFVPRAVLMDLEPGTMDSVRSGPYGQIFRPDNFVFGQSGAGNNWAKGHYTEGAELIDS LDVVRKEAENCDCLQGFQVCHSLG GTGSGMGTLLISKIREEYPDRMMLTFSVFPSPKVSDTVVEPYNATLSVH LVENADECMV",
        )
        self.assertEqual(
            str(hsp),
            """\
Query : 4 Length: 201 Strand: Plus
        gi|5049839|gb|AI730987.1|AI730987 BNLGHi8354 Six-day Cotton fiber
        Gossypium hirsutum cDNA 5' similar to TUBULIN BETA-1 CHAIN
        gi|486734|pir|S35142 tubulin beta chain - white lupine gi|402636
        (X70184) Beta tubulin 1 [Lupinus albus], mRNA sequence
Target: gi|166343825|gb|ABY86655.1| Length: 448 Strand: Plus
        beta-tubulin 4 [Gossypium hirsutum]

Score:408 bits(1048), Expect:2e-112,
Identities:196/201(98%),  Positives:197/201(98%),  Gaps:0.201(0%)

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
        self.assertEqual(hit.num, 2)
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
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 1044.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 406.756589537761)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 6.57981456092236e-112, places=126
        )
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
Query : 4 Length: 201 Strand: Plus
        gi|5049839|gb|AI730987.1|AI730987 BNLGHi8354 Six-day Cotton fiber
        Gossypium hirsutum cDNA 5' similar to TUBULIN BETA-1 CHAIN
        gi|486734|pir|S35142 tubulin beta chain - white lupine gi|402636
        (X70184) Beta tubulin 1 [Lupinus albus], mRNA sequence
Target: gi|223549899|gb|EEF51386.1| Length: 448 Strand: Plus
        tubulin beta chain, putative [Ricinus communis]

Score:406 bits(1044), Expect:7e-112,
Identities:195/201(97%),  Positives:196/201(98%),  Gaps:0.201(0%)

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
        self.assertEqual(hit.num, 3)
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
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 1040.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 405.215791234091)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 1.91443295113426e-111, places=125
        )
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
Query : 4 Length: 201 Strand: Plus
        gi|5049839|gb|AI730987.1|AI730987 BNLGHi8354 Six-day Cotton fiber
        Gossypium hirsutum cDNA 5' similar to TUBULIN BETA-1 CHAIN
        gi|486734|pir|S35142 tubulin beta chain - white lupine gi|402636
        (X70184) Beta tubulin 1 [Lupinus albus], mRNA sequence
Target: gi|18420724|ref|NP_568437.1| Length: 449 Strand: Plus
        TUB8 (tubulin beta-8) [Arabidopsis thaliana]
        >gi|27735261|sp|P29516.2|TBB8_ARATH RecName: Full=Tubulin beta-8 chain;
        AltName: Full=Beta-8-tubulin >gi|10176853|dbj|BAB10059.1| beta tubulin
        [Arabidopsis thaliana]

Score:405 bits(1040), Expect:2e-111,
Identities:194/201(97%),  Positives:196/201(98%),  Gaps:0.201(0%)

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
        self.assertEqual(hit.num, 4)
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
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 1034.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 402.904593778587)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 9.50123195540709e-111, places=125
        )
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
Query : 4 Length: 201 Strand: Plus
        gi|5049839|gb|AI730987.1|AI730987 BNLGHi8354 Six-day Cotton fiber
        Gossypium hirsutum cDNA 5' similar to TUBULIN BETA-1 CHAIN
        gi|486734|pir|S35142 tubulin beta chain - white lupine gi|402636
        (X70184) Beta tubulin 1 [Lupinus albus], mRNA sequence
Target: gi|225426385|ref|XP_002271992.1| Length: 447 Strand: Plus
        PREDICTED: hypothetical protein [Vitis vinifera]
        >gi|157356601|emb|CAO62796.1| unnamed protein product [Vitis vinifera]

Score:402 bits(1034), Expect:1e-110,
Identities:193/201(96%),  Positives:195/201(97%),  Gaps:0.201(0%)

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
        self.assertEqual(hit.num, 5)
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
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 1034.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 402.904593778587)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 9.50123195540709e-111, places=125
        )
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
Query : 4 Length: 201 Strand: Plus
        gi|5049839|gb|AI730987.1|AI730987 BNLGHi8354 Six-day Cotton fiber
        Gossypium hirsutum cDNA 5' similar to TUBULIN BETA-1 CHAIN
        gi|486734|pir|S35142 tubulin beta chain - white lupine gi|402636
        (X70184) Beta tubulin 1 [Lupinus albus], mRNA sequence
Target: gi|15451226|gb|AAK96884.1| Length: 449 Strand: Plus
        beta tubulin [Arabidopsis thaliana] >gi|20148289|gb|AAM10035.1| beta
        tubulin [Arabidopsis thaliana]

Score:402 bits(1034), Expect:1e-110,
Identities:193/201(96%),  Positives:195/201(97%),  Gaps:0.201(0%)

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
        self.assertEqual(hit.num, 6)
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
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 1033.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 402.51939420267)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 1.24089932237309e-110, places=124
        )
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
Query : 4 Length: 201 Strand: Plus
        gi|5049839|gb|AI730987.1|AI730987 BNLGHi8354 Six-day Cotton fiber
        Gossypium hirsutum cDNA 5' similar to TUBULIN BETA-1 CHAIN
        gi|486734|pir|S35142 tubulin beta chain - white lupine gi|402636
        (X70184) Beta tubulin 1 [Lupinus albus], mRNA sequence
Target: gi|225470745|ref|XP_002267380.1| Length: 449 Strand: Plus
        PREDICTED: hypothetical protein [Vitis vinifera]
        >gi|157327486|emb|CAO15467.1| unnamed protein product [Vitis vinifera]

Score:402 bits(1033), Expect:1e-110,
Identities:192/201(96%),  Positives:195/201(97%),  Gaps:0.201(0%)

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
        self.assertEqual(hit.num, 7)
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
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 1033.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 402.51939420267)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 1.24089932237309e-110, places=124
        )
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
Query : 4 Length: 201 Strand: Plus
        gi|5049839|gb|AI730987.1|AI730987 BNLGHi8354 Six-day Cotton fiber
        Gossypium hirsutum cDNA 5' similar to TUBULIN BETA-1 CHAIN
        gi|486734|pir|S35142 tubulin beta chain - white lupine gi|402636
        (X70184) Beta tubulin 1 [Lupinus albus], mRNA sequence
Target: gi|586076|sp|P37392.1|TBB1_LUPAL Length: 447 Strand: Plus
        RecName: Full=Tubulin beta-1 chain; AltName: Full=Beta-1-tubulin
        >gi|402636|emb|CAA49736.1| Beta tubulin 1 [Lupinus albus]

Score:402 bits(1033), Expect:1e-110,
Identities:193/201(96%),  Positives:195/201(97%),  Gaps:0.201(0%)

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
        self.assertEqual(hit.num, 8)
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
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 1031.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 401.748995050835)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 2.1166536544662e-110, places=123
        )
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
Query : 4 Length: 201 Strand: Plus
        gi|5049839|gb|AI730987.1|AI730987 BNLGHi8354 Six-day Cotton fiber
        Gossypium hirsutum cDNA 5' similar to TUBULIN BETA-1 CHAIN
        gi|486734|pir|S35142 tubulin beta chain - white lupine gi|402636
        (X70184) Beta tubulin 1 [Lupinus albus], mRNA sequence
Target: gi|224104341|ref|XP_002313404.1| Length: 451 Strand: Plus
        tubulin, beta chain [Populus trichocarpa] >gi|222849812|gb|EEE87359.1|
        tubulin, beta chain [Populus trichocarpa]

Score:401 bits(1031), Expect:2e-110,
Identities:192/201(96%),  Positives:195/201(97%),  Gaps:0.201(0%)

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
        self.assertEqual(hit.num, 9)
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
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 1029.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 400.978595899)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 3.61046429165375e-110, places=124
        )
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
Query : 4 Length: 201 Strand: Plus
        gi|5049839|gb|AI730987.1|AI730987 BNLGHi8354 Six-day Cotton fiber
        Gossypium hirsutum cDNA 5' similar to TUBULIN BETA-1 CHAIN
        gi|486734|pir|S35142 tubulin beta chain - white lupine gi|402636
        (X70184) Beta tubulin 1 [Lupinus albus], mRNA sequence
Target: gi|223549679|gb|EEF51167.1| Length: 446 Strand: Plus
        tubulin beta chain, putative [Ricinus communis]

Score:400 bits(1029), Expect:4e-110,
Identities:191/201(95%),  Positives:194/201(97%),  Gaps:0.201(0%)

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
        self.assertEqual(hit.num, 10)
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
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 1029.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 400.978595899)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 3.61046429165375e-110, places=124
        )
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
Query : 4 Length: 201 Strand: Plus
        gi|5049839|gb|AI730987.1|AI730987 BNLGHi8354 Six-day Cotton fiber
        Gossypium hirsutum cDNA 5' similar to TUBULIN BETA-1 CHAIN
        gi|486734|pir|S35142 tubulin beta chain - white lupine gi|402636
        (X70184) Beta tubulin 1 [Lupinus albus], mRNA sequence
Target: gi|224058553|ref|XP_002299541.1| Length: 447 Strand: Plus
        tubulin, beta chain [Populus trichocarpa] >gi|222846799|gb|EEE84346.1|
        tubulin, beta chain [Populus trichocarpa]

Score:400 bits(1029), Expect:4e-110,
Identities:192/201(96%),  Positives:195/201(97%),  Gaps:0.201(0%)

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
        self.assertEqual(record.num, 5)
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
        self.assertEqual(record.stat["eff-space"], 147032237429)
        self.assertAlmostEqual(record.stat["kappa"], 0.041)
        self.assertAlmostEqual(record.stat["lambda"], 0.267)
        self.assertAlmostEqual(record.stat["entropy"], 0.14)
        self.assertEqual(len(record), 10)
        hit = record[0]
        self.assertEqual(hit.num, 1)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|110740644|dbj|BAE98425.1|")
        self.assertEqual(hit.target.name, "BAE98425")
        self.assertEqual(
            hit.target.description, "hypothetical protein [Arabidopsis thaliana]"
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=80)")
        self.assertEqual(len(hit), 2)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 231.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 93.5893343169526)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 5.57283114448317e-19, places=33
        )
        self.assertEqual(hsp.annotations["identity"], 42)
        self.assertEqual(hsp.annotations["positive"], 45)
        self.assertEqual(hsp.annotations["gaps"], 1)
        hsp = hit[1]
        self.assertEqual(hsp.num, 2)
        self.assertAlmostEqual(hsp.score, 53.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 25.0238098036637)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 5.57283114448317e-19, places=33
        )
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
Query : 5 Length: 17 Strand: Plus
        gi|5052071|gb|AF067555.1|AF067555 Phlox stansburyi internal transcribed
        spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2,
        complete sequence
Target: gi|110740644|dbj|BAE98425.1| Length: 80 Strand: Plus
        hypothetical protein [Arabidopsis thaliana]

Score:25 bits(53), Expect:6e-19,
Identities:13/17(76%),  Positives:13/17(76%),  Gaps:0.17(0%)

gi|110740        53 RKLVSRVLPHAVGLNPS 70
                  0 ||||||.|..|||||.| 17
5                 0 RKLVSRALRCAVGLNKS 17

""",
        )
        hit = record[1]
        self.assertEqual(hit.num, 2)
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
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 238.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 96.2857313483741)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 1.69151855577931e-18, places=32
        )
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
Query : 5 Length: 49 Strand: Plus
        gi|5052071|gb|AF067555.1|AF067555 Phlox stansburyi internal transcribed
        spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2,
        complete sequence
Target: gi|226453533|gb|EEH50844.1| Length: 81 Strand: Plus
        predicted protein [Micromonas pusilla CCMP1545]

Score:96 bits(238), Expect:2e-18,
Identities:42/49(86%),  Positives:45/49(92%),  Gaps:0.49(0%)

gi|226453         0 MKNVAKCDTWCELQNPVNHRVFERKLRPKPSGRGHVCLGVTNRRPPSSF 49
                  0 ||||||||||||||||||||||||||||||.|||||||||..|..|..| 49
5                 0 MKNVAKCDTWCELQNPVNHRVFERKLRPKPLGRGHVCLGVSHRVAPNPF 49

""",
        )
        hit = record[2]
        self.assertEqual(hit.num, 3)
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
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 183.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 75.0997546729196)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 4.03544314604194e-12, places=26
        )
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
Query : 5 Length: 42 Strand: Plus
        gi|5052071|gb|AF067555.1|AF067555 Phlox stansburyi internal transcribed
        spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2,
        complete sequence
Target: gi|168069582|ref|XP_001786502.1| Length: 88 Strand: Plus
        predicted protein [Physcomitrella patens subsp. patens]
        >gi|162661153|gb|EDQ48685.1| predicted protein [Physcomitrella patens
        subsp. patens]

Score:75 bits(183), Expect:4e-12,
Identities:37/42(88%),  Positives:39/42(93%),  Gaps:0.42(0%)

gi|168069         2 ASGATCVQKLDDSRDSAIHTTYRISLRSSSLQEPRYPLLRVV 44
                  0 |||||||||||.||||||||.|||||||||..||||||.||| 42
5                 0 ASGATCVQKLDGSRDSAIHTKYRISLRSSSMREPRYPLPRVV 42

""",
        )
        hit = record[3]
        self.assertEqual(hit.num, 4)
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
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 178.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 73.1737567933329)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 1.53346675969648e-11, places=25
        )
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
Query : 5 Length: 42 Strand: Plus
        gi|5052071|gb|AF067555.1|AF067555 Phlox stansburyi internal transcribed
        spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2,
        complete sequence
Target: gi|168068558|ref|XP_001786120.1| Length: 130 Strand: Plus
        predicted protein [Physcomitrella patens subsp. patens]
        >gi|162662102|gb|EDQ49068.1| predicted protein [Physcomitrella patens
        subsp. patens]

Score:73 bits(178), Expect:2e-11,
Identities:36/42(86%),  Positives:39/42(93%),  Gaps:0.42(0%)

gi|168068         2 ASGATCVQKLDDSRNSAIHTTYRISLRSSSLQEPRYPLLRVV 44
                  0 |||||||||||.||.|||||.|||||||||..||||||.||| 42
5                 0 ASGATCVQKLDGSRDSAIHTKYRISLRSSSMREPRYPLPRVV 42

""",
        )
        hit = record[4]
        self.assertEqual(hit.num, 5)
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
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 178.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 73.1737567933329)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 1.53346675969648e-11, places=25
        )
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
Query : 5 Length: 42 Strand: Plus
        gi|5052071|gb|AF067555.1|AF067555 Phlox stansburyi internal transcribed
        spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2,
        complete sequence
Target: gi|168068926|ref|XP_001786259.1| Length: 148 Strand: Plus
        predicted protein [Physcomitrella patens subsp. patens]
        >gi|168069965|ref|XP_001786641.1| predicted protein [Physcomitrella
        patens subsp. patens] >gi|162660807|gb|EDQ48545.1| predicted protein
        [Physcomitrella patens subsp. patens] >gi|162661808|gb|EDQ48929.1|
        predicted protein [Physcomitrella patens subsp. patens]

Score:73 bits(178), Expect:2e-11,
Identities:36/42(86%),  Positives:39/42(93%),  Gaps:0.42(0%)

gi|168068         2 ASGATCVQKLDDSRNSAIHTTYRISLRSSSLQEPRYPLLRVV 44
                  0 |||||||||||.||.|||||.|||||||||..||||||.||| 42
5                 0 ASGATCVQKLDGSRDSAIHTKYRISLRSSSMREPRYPLPRVV 42

""",
        )
        hit = record[5]
        self.assertEqual(hit.num, 6)
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
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 178.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 73.1737567933329)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 1.53346675969648e-11, places=25
        )
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
Query : 5 Length: 42 Strand: Plus
        gi|5052071|gb|AF067555.1|AF067555 Phlox stansburyi internal transcribed
        spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2,
        complete sequence
Target: gi|168070288|ref|XP_001786759.1| Length: 148 Strand: Plus
        predicted protein [Physcomitrella patens subsp. patens]
        >gi|162660550|gb|EDQ48427.1| predicted protein [Physcomitrella patens
        subsp. patens]

Score:73 bits(178), Expect:2e-11,
Identities:36/42(86%),  Positives:39/42(93%),  Gaps:0.42(0%)

gi|168070         2 ASGATCVQKLDDSRNSAIHTTYRISLRSSSLQEPRYPLLRVV 44
                  0 |||||||||||.||.|||||.|||||||||..||||||.||| 42
5                 0 ASGATCVQKLDGSRDSAIHTKYRISLRSSSMREPRYPLPRVV 42

""",
        )
        hit = record[6]
        self.assertEqual(hit.num, 7)
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
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 172.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 70.8625593378288)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 7.61051640442713e-11, places=25
        )
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
Query : 5 Length: 80 Strand: Plus
        gi|5052071|gb|AF067555.1|AF067555 Phlox stansburyi internal transcribed
        spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2,
        complete sequence
Target: gi|168068591|ref|XP_001786133.1| Length: 220 Strand: Plus
        predicted protein [Physcomitrella patens subsp. patens]
        >gi|162662081|gb|EDQ49057.1| predicted protein [Physcomitrella patens
        subsp. patens]

Score:70 bits(172), Expect:8e-11,
Identities:42/83(51%),  Positives:50/83(60%),  Gaps:8.83(10%)

gi|168068       142 RGLCHHADSDGQFHSTLPIKDIKRIGGCRDDALAGMPSDEP-----RAQLAFKNSMIHGI
                  0 |...|.......|.|.......|..||---||......|.|-----||||||||||.|||
5                 0 RPTAHRSARETNFRSQTVESRRKWVGG---DAM*DAQADVPSA*WLRAQLAFKNSMVHGI

gi|168068       197 LQFTLRIAFRCVLHRCKSQDIRC 220
                 60 ||||..||||.|||||.|.||||  83
5                57 LQFTPSIAFRYVLHRCESRDIRC  80

""",
        )
        hit = record[7]
        self.assertEqual(hit.num, 8)
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
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 141.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 58.9213724843908)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 2.99274389212967e-07, places=21
        )
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
Query : 5 Length: 39 Strand: Plus
        gi|5052071|gb|AF067555.1|AF067555 Phlox stansburyi internal transcribed
        spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2,
        complete sequence
Target: gi|74622391|sp|Q8TGM5|ART3_YEAST Length: 67 Strand: Plus
        Uncharacterized protein ART3 (Antisense to ribosomal RNA transcript
        protein 3) >gi|18767126|gb|AAL79278.1| unknown [Saccharomyces
        cerevisiae]

Score:58 bits(141), Expect:3e-07,
Identities:29/39(74%),  Positives:32/39(82%),  Gaps:0.39(0%)

gi|746223         7 GAMCVQRFDDSRNSAIHITYRISLRSSSMREPRDPLLKV 46
                  0 ||.|||..|.||.||||..||||||||||||||.||..| 39
5                 0 GATCVQKLDGSRDSAIHTKYRISLRSSSMREPRYPLPRV 39

""",
        )
        hit = record[8]
        self.assertEqual(hit.num, 9)
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
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 137.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 57.3805741807214)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 8.70755166175354e-07, places=21
        )
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
Query : 5 Length: 34 Strand: Plus
        gi|5052071|gb|AF067555.1|AF067555 Phlox stansburyi internal transcribed
        spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2,
        complete sequence
Target: gi|168069944|ref|XP_001786634.1| Length: 138 Strand: Plus
        predicted protein [Physcomitrella patens subsp. patens]
        >gi|162660825|gb|EDQ48552.1| predicted protein [Physcomitrella patens
        subsp. patens]

Score:57 bits(137), Expect:9e-07,
Identities:28/34(82%),  Positives:31/34(91%),  Gaps:0.34(0%)

gi|168069         0 KLDDSRNSAIHTTYRISLRSSSLQEPRYPLLRVV 34
                  0 |||.||.|||||.|||||||||..||||||.||| 34
5                 0 KLDGSRDSAIHTKYRISLRSSSMREPRYPLPRVV 34

""",
        )
        hit = record[9]
        self.assertEqual(hit.num, 10)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|50307717|ref|XP_453851.1|")
        self.assertEqual(hit.target.name, "XP_453851")
        self.assertEqual(
            hit.target.description, "unnamed protein product [Kluyveromyces lactis]"
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=54)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 134.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 56.2249754529693)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 1.93984013155423e-06, places=20
        )
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
Query : 5 Length: 40 Strand: Plus
        gi|5052071|gb|AF067555.1|AF067555 Phlox stansburyi internal transcribed
        spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2,
        complete sequence
Target: gi|50307717|ref|XP_453851.1| Length: 54 Strand: Plus
        unnamed protein product [Kluyveromyces lactis]

Score:56 bits(134), Expect:2e-06,
Identities:28/40(70%),  Positives:31/40(78%),  Gaps:0.40(0%)

gi|503077         7 GAMCVQRFDDSRKSAIHNTYRNSLRSSSMREPRDPLLKVL 47
                  0 ||.|||..|.||.||||..||.|||||||||||.||..|. 40
5                 0 GATCVQKLDGSRDSAIHTKYRISLRSSSMREPRYPLPRVV 40

""",
        )
        record = next(records)
        self.assertEqual(record.num, 6)
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
        self.assertEqual(record.stat["eff-space"], 75367093081)
        self.assertAlmostEqual(record.stat["kappa"], 0.041)
        self.assertAlmostEqual(record.stat["lambda"], 0.267)
        self.assertAlmostEqual(record.stat["entropy"], 0.14)
        self.assertEqual(len(record), 10)
        hit = record[0]
        self.assertEqual(hit.num, 1)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|3176603|gb|AAC18749.1|")
        self.assertEqual(hit.target.name, "AAC18749")
        self.assertEqual(hit.target.description, "phytochrome A [Lathyrus odoratus]")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=103)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 543.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 213.771602003167)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 3.7262743863676e-54, places=67
        )
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
        self.assertEqual(hsp.target.description, "phytochrome A [Lathyrus odoratus]")
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIASLVMAVVVNDSDEDGDSRDAVLPQKKKRLWGLVVCHNTTPRFV",
        )
        self.assertEqual(
            str(hsp),
            """\
Query : 6 Length: 103 Strand: Plus
        gi|3176602|gb|U78617.1|LOU78617 Lathyrus odoratus phytochrome A (PHYA)
        gene, partial cds
Target: gi|3176603|gb|AAC18749.1| Length: 103 Strand: Plus
        phytochrome A [Lathyrus odoratus]

Score:213 bits(543), Expect:4e-54,
Identities:103/103(100%),  Positives:103/103(100%),  Gaps:0.103(0%)

gi|317660         0 QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
6                 0 QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA

gi|317660        60 SLVMAVVVNDSDEDGDSRDAVLPQKKKRLWGLVVCHNTTPRFV 103
                 60 ||||||||||||||||||||||||||||||||||||||||||| 103
6                60 SLVMAVVVNDSDEDGDSRDAVLPQKKKRLWGLVVCHNTTPRFV 103

""",
        )
        hit = record[1]
        self.assertEqual(hit.num, 2)
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
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 530.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 208.764007516241)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 1.1987013044853e-52, places=65
        )
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
Query : 6 Length: 103 Strand: Plus
        gi|3176602|gb|U78617.1|LOU78617 Lathyrus odoratus phytochrome A (PHYA)
        gene, partial cds
Target: gi|130188|sp|P15001.1|PHYA_PEA Length: 1124 Strand: Plus
        RecName: Full=Phytochrome A >gi|169132|gb|AAA33682.1| phytochrome [Pisum
        sativum] >gi|295830|emb|CAA32242.1| phytochrome apoprotein [Pisum
        sativum] >gi|51173514|gb|AAT97643.1| phytochrome A apoprotein [Pisum
        sativum] >gi|226757|prf||1604466A phytochrome

Score:208 bits(530), Expect:1e-52,
Identities:101/103(98%),  Positives:102/103(99%),  Gaps:0.103(0%)

gi|130188       275 QAARFLFMKNKVRMIVDCNAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA
                  0 ||||||||||||||||||.|||||||||||||||||||||||||||||||||||||||||
6                 0 QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA

gi|130188       335 SLVMAVVVNDSDEDGDSADAVLPQKKKRLWGLVVCHNTTPRFV 378
                 60 |||||||||||||||||.||||||||||||||||||||||||| 103
6                60 SLVMAVVVNDSDEDGDSRDAVLPQKKKRLWGLVVCHNTTPRFV 103

""",
        )
        hit = record[2]
        self.assertEqual(hit.num, 3)
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
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 530.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 208.764007516241)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 1.1987013044853e-52, places=65
        )
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
Query : 6 Length: 103 Strand: Plus
        gi|3176602|gb|U78617.1|LOU78617 Lathyrus odoratus phytochrome A (PHYA)
        gene, partial cds
Target: gi|2499555|sp|P93673.1|PHYA_LATSA Length: 1124 Strand: Plus
        RecName: Full=Phytochrome type A >gi|1848273|gb|AAB47994.1| phytochrome
        type A [Lathyrus sativus]

Score:208 bits(530), Expect:1e-52,
Identities:101/103(98%),  Positives:102/103(99%),  Gaps:0.103(0%)

gi|249955       275 QAARFLFMKNKVRMIVDCNAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA
                  0 ||||||||||||||||||.|||||||||||||||||||||||||||||||||||||||||
6                 0 QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA

gi|249955       335 SLVMAVVVNDSDEDGDSADAVLPQKKKRLWGLVVCHNTTPRFV 378
                 60 |||||||||||||||||.||||||||||||||||||||||||| 103
6                60 SLVMAVVVNDSDEDGDSRDAVLPQKKKRLWGLVVCHNTTPRFV 103

""",
        )
        hit = record[3]
        self.assertEqual(hit.num, 4)
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
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 528.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 207.993608364407)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 2.04467473791515e-52, places=66
        )
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
Query : 6 Length: 103 Strand: Plus
        gi|3176602|gb|U78617.1|LOU78617 Lathyrus odoratus phytochrome A (PHYA)
        gene, partial cds
Target: gi|3176595|gb|AAC18745.1| Length: 103 Strand: Plus
        phytochrome A [Lennea melanocarpa] >gi|3176597|gb|AAC18746.1|
        phytochrome A [Hebestigma cubense] >gi|3176609|gb|AAC18752.1|
        phytochrome A [Sesbania cochichinensis] >gi|3176611|gb|AAC18753.1|
        phytochrome A [Sesbania emerus]

Score:207 bits(528), Expect:2e-52,
Identities:100/103(97%),  Positives:101/103(98%),  Gaps:0.103(0%)

gi|317659         0 QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
6                 0 QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA

gi|317659        60 SLVMAVVVNDSDEDGDSSDAVQPQKRKRLWGLVVCHNTTPRFV 103
                 60 |||||||||||||||||.|||.|||.||||||||||||||||| 103
6                60 SLVMAVVVNDSDEDGDSRDAVLPQKKKRLWGLVVCHNTTPRFV 103

""",
        )
        hit = record[4]
        self.assertEqual(hit.num, 5)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|1711106|gb|AAC18675.1|")
        self.assertEqual(hit.target.name, "AAC18675")
        self.assertEqual(hit.target.description, "phytochrome A [Sophora affinis]")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=210)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 528.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 207.993608364407)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 2.04467473791515e-52, places=66
        )
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
Query : 6 Length: 103 Strand: Plus
        gi|3176602|gb|U78617.1|LOU78617 Lathyrus odoratus phytochrome A (PHYA)
        gene, partial cds
Target: gi|1711106|gb|AAC18675.1| Length: 210 Strand: Plus
        phytochrome A [Sophora affinis]

Score:207 bits(528), Expect:2e-52,
Identities:100/103(97%),  Positives:101/103(98%),  Gaps:0.103(0%)

gi|171110        40 QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
6                 0 QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA

gi|171110       100 SLVMAVVVNDSDEDGDSSDAVQPQKRKRLWGLVVCHNTTPRFV 143
                 60 |||||||||||||||||.|||.|||.||||||||||||||||| 103
6                60 SLVMAVVVNDSDEDGDSRDAVLPQKKKRLWGLVVCHNTTPRFV 103

""",
        )
        hit = record[5]
        self.assertEqual(hit.num, 6)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|1711090|gb|AAC18670.1|")
        self.assertEqual(hit.target.name, "AAC18670")
        self.assertEqual(hit.target.description, "phytochrome A [Myrospermum sousanum]")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=210)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 525.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 206.838009636654)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 4.55506009801166e-52, places=66
        )
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
        self.assertEqual(hsp.target.description, "phytochrome A [Myrospermum sousanum]")
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIASLV+AVVVNDSDEDGDS DAV PQK+KRLWGLVVCHNTTPRFV",
        )
        self.assertEqual(
            str(hsp),
            """\
Query : 6 Length: 103 Strand: Plus
        gi|3176602|gb|U78617.1|LOU78617 Lathyrus odoratus phytochrome A (PHYA)
        gene, partial cds
Target: gi|1711090|gb|AAC18670.1| Length: 210 Strand: Plus
        phytochrome A [Myrospermum sousanum]

Score:206 bits(525), Expect:5e-52,
Identities:99/103(96%),  Positives:101/103(98%),  Gaps:0.103(0%)

gi|171109        40 QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
6                 0 QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA

gi|171109       100 SLVLAVVVNDSDEDGDSSDAVQPQKRKRLWGLVVCHNTTPRFV 143
                 60 |||.|||||||||||||.|||.|||.||||||||||||||||| 103
6                60 SLVMAVVVNDSDEDGDSRDAVLPQKKKRLWGLVVCHNTTPRFV 103

""",
        )
        hit = record[6]
        self.assertEqual(hit.num, 7)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|3176605|gb|AAC18750.1|")
        self.assertEqual(hit.target.name, "AAC18750")
        self.assertEqual(hit.target.description, "phytochrome A [Hybosema robustum]")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=103)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 524.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 206.452810060737)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 5.94909272347008e-52, places=66
        )
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
        self.assertEqual(hsp.target.description, "phytochrome A [Hybosema robustum]")
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "QA RFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIASLVMAVVVNDSDEDGDS DAV PQK+KRLWGLVVCHNTTPRFV",
        )
        self.assertEqual(
            str(hsp),
            """\
Query : 6 Length: 103 Strand: Plus
        gi|3176602|gb|U78617.1|LOU78617 Lathyrus odoratus phytochrome A (PHYA)
        gene, partial cds
Target: gi|3176605|gb|AAC18750.1| Length: 103 Strand: Plus
        phytochrome A [Hybosema robustum]

Score:206 bits(524), Expect:6e-52,
Identities:99/103(96%),  Positives:100/103(97%),  Gaps:0.103(0%)

gi|317660         0 QATRFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA
                  0 ||.|||||||||||||||||||||||||||||||||||||||||||||||||||||||||
6                 0 QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA

gi|317660        60 SLVMAVVVNDSDEDGDSSDAVQPQKRKRLWGLVVCHNTTPRFV 103
                 60 |||||||||||||||||.|||.|||.||||||||||||||||| 103
6                60 SLVMAVVVNDSDEDGDSRDAVLPQKKKRLWGLVVCHNTTPRFV 103

""",
        )
        hit = record[7]
        self.assertEqual(hit.num, 8)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|3176454|gb|AAC18668.1|")
        self.assertEqual(hit.target.name, "AAC18668")
        self.assertEqual(hit.target.description, "phytochrome A [Cyclolobium nutans]")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=207)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 523.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 206.06761048482)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 7.76975571582328e-52, places=66
        )
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
        self.assertEqual(hsp.target.description, "phytochrome A [Cyclolobium nutans]")
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIASLVMAVVVNDSDEDG+S DAV PQK+KRLWGLVVCHNTTPRFV",
        )
        self.assertEqual(
            str(hsp),
            """\
Query : 6 Length: 103 Strand: Plus
        gi|3176602|gb|U78617.1|LOU78617 Lathyrus odoratus phytochrome A (PHYA)
        gene, partial cds
Target: gi|3176454|gb|AAC18668.1| Length: 207 Strand: Plus
        phytochrome A [Cyclolobium nutans]

Score:206 bits(523), Expect:8e-52,
Identities:99/103(96%),  Positives:101/103(98%),  Gaps:0.103(0%)

gi|317645        37 QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
6                 0 QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA

gi|317645        97 SLVMAVVVNDSDEDGNSSDAVQPQKRKRLWGLVVCHNTTPRFV 140
                 60 |||||||||||||||.|.|||.|||.||||||||||||||||| 103
6                60 SLVMAVVVNDSDEDGDSRDAVLPQKKKRLWGLVVCHNTTPRFV 103

""",
        )
        hit = record[8]
        self.assertEqual(hit.num, 9)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|3176523|gb|AAC18709.1|")
        self.assertEqual(hit.target.name, "AAC18709")
        self.assertEqual(
            hit.target.description, "phytochrome A [Millettia richardiana]"
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=139)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 521.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 205.297211332985)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 1.3253195915005e-51, places=64
        )
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
Query : 6 Length: 103 Strand: Plus
        gi|3176602|gb|U78617.1|LOU78617 Lathyrus odoratus phytochrome A (PHYA)
        gene, partial cds
Target: gi|3176523|gb|AAC18709.1| Length: 139 Strand: Plus
        phytochrome A [Millettia richardiana]

Score:205 bits(521), Expect:1e-51,
Identities:98/103(95%),  Positives:101/103(98%),  Gaps:0.103(0%)

gi|317652        36 QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
6                 0 QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA

gi|317652        96 SLVMAVVVNDNEEDGDSSDAVQPQKRKRLWGLVVCHNTTPRFV 139
                 60 ||||||||||..|||||.|||.|||.||||||||||||||||| 103
6                60 SLVMAVVVNDSDEDGDSRDAVLPQKKKRLWGLVVCHNTTPRFV 103

""",
        )
        hit = record[9]
        self.assertEqual(hit.num, 10)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|3176494|gb|AAC18693.1|")
        self.assertEqual(hit.target.name, "AAC18693")
        self.assertEqual(
            hit.target.description, "phytochrome A [Callerya atropurpurea]"
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=177)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 520.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 204.912011757068)
        self.assertAlmostEqual(
            hsp.annotations["evalue"], 1.73092099081406e-51, places=65
        )
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
Query : 6 Length: 103 Strand: Plus
        gi|3176602|gb|U78617.1|LOU78617 Lathyrus odoratus phytochrome A (PHYA)
        gene, partial cds
Target: gi|3176494|gb|AAC18693.1| Length: 177 Strand: Plus
        phytochrome A [Callerya atropurpurea]

Score:204 bits(520), Expect:2e-51,
Identities:98/103(95%),  Positives:101/103(98%),  Gaps:0.103(0%)

gi|317649         7 QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
6                 0 QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA

gi|317649        67 SLVMAVVVNDSEEDGDSSEAVQPQKRKRLWGLVVCHNTTPRFV 110
                 60 |||||||||||.|||||..||.|||.||||||||||||||||| 103
6                60 SLVMAVVVNDSDEDGDSRDAVLPQKKKRLWGLVVCHNTTPRFV 103

""",
        )
        record = next(records)
        self.assertEqual(record.num, 7)
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
        self.assertEqual(record.stat["eff-space"], 1251086325060)
        self.assertAlmostEqual(record.stat["kappa"], 0.041)
        self.assertAlmostEqual(record.stat["lambda"], 0.267)
        self.assertAlmostEqual(record.stat["entropy"], 0.14)
        self.assertEqual(len(record), 10)
        hit = record[0]
        self.assertEqual(hit.num, 1)
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
        self.assertEqual(hsp.num, 1)
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
Query : 7 Length: 506 Strand: Plus
        gi|5817701|gb|AF142731.1|AF142731 Wisteria frutescens maturase-like
        protein (matK) gene, complete cds; chloroplast gene for chloroplast
        product
Target: gi|27805603|sp|Q9TKP6.1|MATK_WISFR Length: 506 Strand: Plus
        RecName: Full=Maturase K; AltName: Full=Intron maturase
        >gi|5817759|gb|AAD52902.1|AF142731_1 maturase-like protein [Wisteria
        frutescens]

Score:948 bits(2451), Expect:0,
Identities:506/506(100%),  Positives:506/506(100%),  Gaps:0.506(0%)

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
        self.assertEqual(hit.num, 2)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|171909144|gb|ACB58148.1|")
        self.assertEqual(hit.target.name, "ACB58148")
        self.assertEqual(hit.target.description, "maturase K [Wisteria frutescens]")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=506)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
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
Query : 7 Length: 506 Strand: Plus
        gi|5817701|gb|AF142731.1|AF142731 Wisteria frutescens maturase-like
        protein (matK) gene, complete cds; chloroplast gene for chloroplast
        product
Target: gi|171909144|gb|ACB58148.1| Length: 506 Strand: Plus
        maturase K [Wisteria frutescens]

Score:946 bits(2445), Expect:0,
Identities:505/506(100%),  Positives:505/506(100%),  Gaps:0.506(0%)

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
        self.assertEqual(hit.num, 3)
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
        self.assertEqual(hsp.num, 1)
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
Query : 7 Length: 506 Strand: Plus
        gi|5817701|gb|AF142731.1|AF142731 Wisteria frutescens maturase-like
        protein (matK) gene, complete cds; chloroplast gene for chloroplast
        product
Target: gi|171909146|gb|ACB58149.1| Length: 506 Strand: Plus
        maturase K [Wisteria frutescens] >gi|171909148|gb|ACB58150.1| maturase K
        [Wisteria frutescens var. macrostachya]

Score:945 bits(2443), Expect:0,
Identities:505/506(100%),  Positives:505/506(100%),  Gaps:0.506(0%)

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
        self.assertEqual(hit.num, 4)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|171909132|gb|ACB58142.1|")
        self.assertEqual(hit.target.name, "ACB58142")
        self.assertEqual(hit.target.description, "maturase K [Callerya megasperma]")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=506)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
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
Query : 7 Length: 506 Strand: Plus
        gi|5817701|gb|AF142731.1|AF142731 Wisteria frutescens maturase-like
        protein (matK) gene, complete cds; chloroplast gene for chloroplast
        product
Target: gi|171909132|gb|ACB58142.1| Length: 506 Strand: Plus
        maturase K [Callerya megasperma]

Score:944 bits(2439), Expect:0,
Identities:501/506(99%),  Positives:504/506(100%),  Gaps:0.506(0%)

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
        self.assertEqual(hit.num, 5)
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
        self.assertEqual(hsp.num, 1)
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
Query : 7 Length: 506 Strand: Plus
        gi|5817701|gb|AF142731.1|AF142731 Wisteria frutescens maturase-like
        protein (matK) gene, complete cds; chloroplast gene for chloroplast
        product
Target: gi|5817760|gb|AAD52903.1|AF142732_1 Length: 506 Strand: Plus
        maturase-like protein [Wisteria sinensis] >gi|171909136|gb|ACB58144.1|
        maturase K [Wisteria brachybotrys] >gi|171909138|gb|ACB58145.1| maturase
        K [Wisteria floribunda] >gi|171909140|gb|ACB58146.1| maturase K
        [Wisteria floribunda] >gi|171909142|gb|ACB58147.1| maturase K [Wisteria
        floribunda] >gi|171909150|gb|ACB58151.1| maturase K [Wisteria sinensis]
        >gi|171909152|gb|ACB58152.1| maturase K [Wisteria sinensis]
        >gi|171909154|gb|ACB58153.1| maturase K [Wisteria sinensis]
        >gi|171909156|gb|ACB58154.1| maturase K [Wisteria villosa]
        >gi|171909158|gb|ACB58155.1| maturase K [Wisteria villosa]
        >gi|171909160|gb|ACB58156.1| maturase K [Wisteria villosa]

Score:936 bits(2418), Expect:0,
Identities:498/506(98%),  Positives:500/506(99%),  Gaps:0.506(0%)

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
        self.assertEqual(hit.num, 6)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|171909134|gb|ACB58143.1|")
        self.assertEqual(hit.target.name, "ACB58143")
        self.assertEqual(hit.target.description, "maturase K [Wisteria brachybotrys]")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=506)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
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
        self.assertEqual(hsp.target.description, "maturase K [Wisteria brachybotrys]")
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MKEYQVYLERDRSRQQDFLYPLIFREYIYGLAYSHDFNRSIFVENVGYD KSSLLIVKRLITRMYQQNHLIISANDSNKNPFLGYN  FYSQIISDGFAVVVEIPFFLQLSSSLEEAEIVKSYHNLRSIHSIFPFLEDK TY NYVSDIRIPYPIHLEILVQILRY VKDASFFHLLRFFLY+F NRNSLITPKKSISTFSKSNPRLFLFLYNFYVCEYESIFRFLRNQSSHLRLKSFSVFFERIFFYAKREHLVKVFPKDFSSTLTFFKDPFIHYVRYQGKSILASKNAPLLMNKWKHYFIHLWQCFFDVWSQPGTIHINQLSEHSFHFLGYFSNVRLNRSVVRSQMLQNTFLIEIVIKKLDIIVPIIPLIRSLAKAKFCNVLGHP+SKSVWADSSDFDIIDRFLRICRNLSHYYNGSSKKKNLYRIKYILRLSCIKTLACKHKSTVRAFLKKSGSEELLEEFFTEEEEILSLIFPR SSTLQRLHRNRIWYLDILFSNDLVNHE",
        )
        self.assertEqual(
            str(hsp),
            """\
Query : 7 Length: 506 Strand: Plus
        gi|5817701|gb|AF142731.1|AF142731 Wisteria frutescens maturase-like
        protein (matK) gene, complete cds; chloroplast gene for chloroplast
        product
Target: gi|171909134|gb|ACB58143.1| Length: 506 Strand: Plus
        maturase K [Wisteria brachybotrys]

Score:928 bits(2398), Expect:0,
Identities:496/506(98%),  Positives:498/506(98%),  Gaps:0.506(0%)

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
        self.assertEqual(hit.num, 7)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|5817761|gb|AAD52904.1|AF142733_1")
        self.assertEqual(hit.target.name, "AAD52904")
        self.assertEqual(
            hit.target.description, "maturase-like protein [Callerya reticulata]"
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=506)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
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
Query : 7 Length: 506 Strand: Plus
        gi|5817701|gb|AF142731.1|AF142731 Wisteria frutescens maturase-like
        protein (matK) gene, complete cds; chloroplast gene for chloroplast
        product
Target: gi|5817761|gb|AAD52904.1|AF142733_1 Length: 506 Strand: Plus
        maturase-like protein [Callerya reticulata]

Score:925 bits(2390), Expect:0,
Identities:493/506(97%),  Positives:498/506(98%),  Gaps:0.506(0%)

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
        self.assertEqual(hit.num, 8)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|5817762|gb|AAD52905.1|AF142734_1")
        self.assertEqual(hit.target.name, "AAD52905")
        self.assertEqual(
            hit.target.description, "maturase-like protein [Callerya atropurpurea]"
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=506)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
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
Query : 7 Length: 506 Strand: Plus
        gi|5817701|gb|AF142731.1|AF142731 Wisteria frutescens maturase-like
        protein (matK) gene, complete cds; chloroplast gene for chloroplast
        product
Target: gi|5817762|gb|AAD52905.1|AF142734_1 Length: 506 Strand: Plus
        maturase-like protein [Callerya atropurpurea]

Score:890 bits(2301), Expect:0,
Identities:472/506(93%),  Positives:488/506(96%),  Gaps:0.506(0%)

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
        self.assertEqual(hit.num, 9)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|152014012|gb|ABS20107.1|")
        self.assertEqual(hit.target.name, "ABS20107")
        self.assertEqual(
            hit.target.description, "maturase-like protein [Astragalus uliginosus]"
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=506)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
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
Query : 7 Length: 506 Strand: Plus
        gi|5817701|gb|AF142731.1|AF142731 Wisteria frutescens maturase-like
        protein (matK) gene, complete cds; chloroplast gene for chloroplast
        product
Target: gi|152014012|gb|ABS20107.1| Length: 506 Strand: Plus
        maturase-like protein [Astragalus uliginosus]

Score:887 bits(2293), Expect:0,
Identities:470/506(93%),  Positives:487/506(96%),  Gaps:0.506(0%)

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
        self.assertEqual(hit.num, 10)
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
        self.assertEqual(hsp.num, 1)
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
Query : 7 Length: 506 Strand: Plus
        gi|5817701|gb|AF142731.1|AF142731 Wisteria frutescens maturase-like
        protein (matK) gene, complete cds; chloroplast gene for chloroplast
        product
Target: gi|146197442|dbj|BAF57483.1| Length: 506 Strand: Plus
        maturase [Glycyrrhiza uralensis] >gi|146197444|dbj|BAF57484.1| maturase
        [Glycyrrhiza glabra] >gi|152014018|gb|ABS20110.1| maturase-like protein
        [Glycyrrhiza pallidiflora]

Score:887 bits(2292), Expect:0,
Identities:471/506(93%),  Positives:489/506(97%),  Gaps:0.506(0%)

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

    def test_xml_2222_blastx_001_writer(self):
        """Writing BLASTX 2.2.22+ (xml_2222_blastx_001.xml)."""
        filename = "xml_2222_blastx_001.xml"
        path = os.path.join("Blast", filename)
        with Blast.parse(path) as records:
            stream = io.BytesIO()
            n = Blast.write(records, stream)
            self.assertEqual(n, 7)
            stream.seek(0)
            written_records = Blast.parse(stream)
            self.check_xml_2222_blastx_001(written_records)

    def test_xml_21500_blastx_001_parser(self):
        """Parsing BLASTX 2.15.0+ (xml_21500_blastx_001.xml)."""
        filename = "xml_21500_blastx_001.xml"
        path = os.path.join("Blast", filename)
        with open(path, "rb") as stream:
            records = Blast.parse(stream)
            self.check_xml_21500_blastx_001_records(records)
        with Blast.parse(path) as records:
            self.check_xml_21500_blastx_001_records(records)
        with open(path, "rb") as stream:
            record = Blast.read(stream)
        self.check_xml_21500_blastx_001_record(record, xml2=False)
        record = Blast.read(path)
        self.check_xml_21500_blastx_001_record(record, xml2=False)
        with Blast.parse(path) as records:
            self.assertEqual(
                str(records),
                """\
Program: BLASTX 2.15.0+
     db: nr

  Query: AI021773.1 (length=365)
         MAAD0534.RAR Schistosoma mansoni, adult worm (J.C.Parra) Schistosoma
         mansoni cDNA clone MAAD0534.RAR 5' end similar to S. mansoni actin
         mRNA, complete cds, mRNA sequence
   Hits: ----  -----  ----------------------------------------------------------
            #  # HSP  ID + description
         ----  -----  ----------------------------------------------------------
            0      1  emb|VDM03167.1|  unnamed protein product, partial [Schi...
            1      1  gb|EPB74633.1|  hypothetical protein ANCCEY_06263 [Ancy...
            2      1  ref|XP_009175831.1|  hypothetical protein T265_11027 [O...
            3      1  gb|EMD49430.1|  actin, putative, partial [Entamoeba his...
            4      1  emb|CAX83035.1|  Actin-2, partial [Schistosoma japonicum]
            5      1  emb|VDP83060.1|  unnamed protein product, partial [Echi...
            6      1  emb|CAA50205.1|  actin, partial [Entamoeba histolytica]
            7      1  emb|VDN44756.1|  unnamed protein product, partial [Dibo...
            8      1  ref|XP_027046469.1|  actin-1, partial [Pocillopora dami...
            9      1  ref|XP_027046487.1|  actin-1-like [Pocillopora damicornis]""",
            )

    def check_xml_21500_blastx_001_records(self, records, xml2=False):
        self.assertEqual(records.program, "blastx")
        self.assertEqual(records.version, "BLASTX 2.15.0+")
        self.assertEqual(
            records.reference,
            'Stephen F. Altschul, Thomas L. Madden, Alejandro A. Schäffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.',
        )
        self.assertEqual(records.db, "nr")
        if xml2:
            self.assertEqual(len(records.param), 7)
        else:
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
        if xml2:
            self.assertEqual(records.param["cbs"], 2)
            self.assertEqual(records.param["query-gencode"], 1)
        record = next(records)
        self.assertRaises(StopIteration, next, records)
        self.check_xml_21500_blastx_001_record(record, xml2)

    def check_xml_21500_blastx_001_record(self, record, xml2):
        if not xml2:
            self.assertEqual(record.num, 1)
        self.assertIsInstance(record.query, SeqRecord)
        self.assertEqual(record.query.id, "AI021773.1")
        self.assertEqual(
            record.query.description,
            "MAAD0534.RAR Schistosoma mansoni, adult worm (J.C.Parra) Schistosoma mansoni cDNA clone MAAD0534.RAR 5' end similar to S. mansoni actin mRNA, complete cds, mRNA sequence",
        )
        self.assertEqual(repr(record.query.seq), "Seq(None, length=365)")

        self.assertEqual(len(record.stat), 7)
        self.assertEqual(record.stat["db-num"], 718367499)
        self.assertEqual(record.stat["db-len"], 277248733561)
        if xml2:
            self.assertEqual(record.stat["hsp-len"], 89)
            self.assertEqual(record.stat["eff-space"], 6826048836800)
        else:
            self.assertEqual(record.stat["hsp-len"], 0)
            self.assertEqual(record.stat["eff-space"], 0)
        self.assertAlmostEqual(record.stat["kappa"], 0.041)
        self.assertAlmostEqual(record.stat["lambda"], 0.267)
        self.assertAlmostEqual(record.stat["entropy"], 0.14)
        self.assertEqual(len(record), 10)
        hit = record[0]
        self.assertEqual(hit.num, 1)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "emb|VDM03167.1|")
        self.assertEqual(hit.target.name, "VDM03167")
        self.assertEqual(
            hit.target.description,
            "unnamed protein product, partial [Schistocephalus solidus]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=132)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 408.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 161.77)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.13203e-48, places=53)
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
        self.assertEqual(hsp.target.id, "emb|VDM03167.1|")
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
        self.maxDiff = None
        self.assertEqual(
            str(hsp),
            """\
Query : AI021773.1 Length: 108 Strand: Plus
        MAAD0534.RAR Schistosoma mansoni, adult worm (J.C.Parra) Schistosoma
        mansoni cDNA clone MAAD0534.RAR 5' end similar to S. mansoni actin mRNA,
        complete cds, mRNA sequence
Target: emb|VDM03167.1| Length: 132 Strand: Plus
        unnamed protein product, partial [Schistocephalus solidus]

Score:161 bits(408), Expect:3e-48,
Identities:81/108(75%),  Positives:83/108(77%),  Gaps:0.108(0%)

emb|VDM03         0 MADEEVQALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQ
                  0 |||||||||||||||||||||...........|...............|.||||||||||
AI021773.         0 MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQ

emb|VDM03        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE 108
                 60 |||||||||||||||||||||||||||||||||||||||||||||||| 108
AI021773.        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE 108

""",
        )
        hit = record[1]
        self.assertEqual(hit.num, 2)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gb|EPB74633.1|")
        self.assertEqual(hit.target.name, "EPB74633")
        self.assertEqual(
            hit.target.description,
            "hypothetical protein ANCCEY_06263 [Ancylostoma ceylanicum]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=119)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 405.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 160.614)
        self.assertAlmostEqual(hsp.annotations["evalue"], 5.40441e-48, places=53)
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
        self.assertEqual(hsp.target.id, "gb|EPB74633.1|")
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
Query : AI021773.1 Length: 115 Strand: Plus
        MAAD0534.RAR Schistosoma mansoni, adult worm (J.C.Parra) Schistosoma
        mansoni cDNA clone MAAD0534.RAR 5' end similar to S. mansoni actin mRNA,
        complete cds, mRNA sequence
Target: gb|EPB74633.1| Length: 119 Strand: Plus
        hypothetical protein ANCCEY_06263 [Ancylostoma ceylanicum]

Score:160 bits(405), Expect:5e-48,
Identities:81/115(70%),  Positives:85/115(74%),  Gaps:0.115(0%)

gb|EPB746         0 MCDDDVAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQ
                  0 |.|..|.||||||||||||||...........|...............|.||||||||||
AI021773.         0 MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQ

gb|EPB746        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTEAHSILKP 115
                 60 ||||||||||||||||||||||||||||||||||||||||||||||||.|.|.|| 115
AI021773.        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTELHCIRKP 115

""",
        )
        hit = record[2]
        self.assertEqual(hit.num, 3)
        target = hit.target
        self.assertIsInstance(target, SeqRecord)
        self.assertEqual(target.id, "ref|XP_009175831.1|")
        self.assertEqual(target.name, "XP_009175831")
        seq = target.seq
        self.assertEqual(repr(seq), "Seq(None, length=246)")
        if xml2:
            self.assertEqual(
                target.description,
                "hypothetical protein T265_11027 [Opisthorchis viverrini]",
            )
            self.assertIs(target, hit.targets[0])
            self.assertEqual(len(hit.targets), 2)
            target = hit.targets[1]
            self.assertIsInstance(target, SeqRecord)
            self.assertEqual(target.id, "gb|KER20427.1|")
            self.assertEqual(target.name, "KER20427")
            self.assertIs(target.seq, seq)
            self.assertEqual(
                target.description,
                "hypothetical protein T265_11027 [Opisthorchis viverrini]",
            )
        else:
            self.assertEqual(
                target.description,
                "hypothetical protein T265_11027 [Opisthorchis viverrini] >gb|KER20427.1| hypothetical protein T265_11027 [Opisthorchis viverrini]",
            )
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 413.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 163.696)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.69953e-47, places=52)
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
        self.assertEqual(hsp.target.id, hit.target.id)
        self.assertEqual(hsp.target.name, hit.target.name)
        self.assertEqual(hsp.target.description, hit.target.description)
        self.assertEqual(len(hit), 1)
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MADEEVQALVVDNGSGMCKAG       ++  P               G KDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE",
        )
        if xml2:
            self.assertEqual(
                str(hsp),
                """\
Query : AI021773.1 Length: 108 Strand: Plus
        MAAD0534.RAR Schistosoma mansoni, adult worm (J.C.Parra) Schistosoma
        mansoni cDNA clone MAAD0534.RAR 5' end similar to S. mansoni actin mRNA,
        complete cds, mRNA sequence
Target: ref|XP_009175831.1| Length: 246 Strand: Plus
        hypothetical protein T265_11027 [Opisthorchis viverrini]

Score:163 bits(413), Expect:2e-47,
Identities:81/108(75%),  Positives:83/108(77%),  Gaps:0.108(0%)

ref|XP_00         0 MADEEVQALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQ
                  0 |||||||||||||||||||||...........|...............|.||||||||||
AI021773.         0 MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQ

ref|XP_00        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE 108
                 60 |||||||||||||||||||||||||||||||||||||||||||||||| 108
AI021773.        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE 108

""",
            )
        else:
            self.assertEqual(
                str(hsp),
                """\
Query : AI021773.1 Length: 108 Strand: Plus
        MAAD0534.RAR Schistosoma mansoni, adult worm (J.C.Parra) Schistosoma
        mansoni cDNA clone MAAD0534.RAR 5' end similar to S. mansoni actin mRNA,
        complete cds, mRNA sequence
Target: ref|XP_009175831.1| Length: 246 Strand: Plus
        hypothetical protein T265_11027 [Opisthorchis viverrini] >gb|KER20427.1|
        hypothetical protein T265_11027 [Opisthorchis viverrini]

Score:163 bits(413), Expect:2e-47,
Identities:81/108(75%),  Positives:83/108(77%),  Gaps:0.108(0%)

ref|XP_00         0 MADEEVQALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQ
                  0 |||||||||||||||||||||...........|...............|.||||||||||
AI021773.         0 MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQ

ref|XP_00        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE 108
                 60 |||||||||||||||||||||||||||||||||||||||||||||||| 108
AI021773.        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE 108

""",
            )
        hit = record[3]
        self.assertEqual(hit.num, 4)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gb|EMD49430.1|")
        self.assertEqual(hit.target.name, "EMD49430")
        self.assertEqual(
            hit.target.description,
            "actin, putative, partial [Entamoeba histolytica KU27]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=124)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 401.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 159.073)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.49189e-47, places=52)
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
        self.assertEqual(hsp.target.id, "gb|EMD49430.1|")
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
Query : AI021773.1 Length: 108 Strand: Plus
        MAAD0534.RAR Schistosoma mansoni, adult worm (J.C.Parra) Schistosoma
        mansoni cDNA clone MAAD0534.RAR 5' end similar to S. mansoni actin mRNA,
        complete cds, mRNA sequence
Target: gb|EMD49430.1| Length: 124 Strand: Plus
        actin, putative, partial [Entamoeba histolytica KU27]

Score:159 bits(401), Expect:3e-47,
Identities:78/108(72%),  Positives:81/108(75%),  Gaps:0.108(0%)

gb|EMD494         0 MGDEEVQALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHVSVMAGMGQKDAYVGDEAQ
                  0 |.|||||||||||||||||||...........|...............|.||.|||||||
AI021773.         0 MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQ

gb|EMD494        60 SKRGILTLKYPIEHGIVNNWDDMEKIWHHTFYNELRVAPEEHPVLLTE 108
                 60 |||||||||||||||||.|||||||||||||||||||||||||||||| 108
AI021773.        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE 108

""",
        )
        hit = record[4]
        self.assertEqual(hit.num, 5)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "emb|CAX83035.1|")
        self.assertEqual(hit.target.name, "CAX83035")
        self.assertEqual(
            hit.target.description, "Actin-2, partial [Schistosoma japonicum]"
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=252)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 411.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 162.925)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.86747e-47, places=52)
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
        self.assertEqual(hsp.target.id, "emb|CAX83035.1|")
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
Query : AI021773.1 Length: 108 Strand: Plus
        MAAD0534.RAR Schistosoma mansoni, adult worm (J.C.Parra) Schistosoma
        mansoni cDNA clone MAAD0534.RAR 5' end similar to S. mansoni actin mRNA,
        complete cds, mRNA sequence
Target: emb|CAX83035.1| Length: 252 Strand: Plus
        Actin-2, partial [Schistosoma japonicum]

Score:162 bits(411), Expect:4e-47,
Identities:81/108(75%),  Positives:83/108(77%),  Gaps:0.108(0%)

emb|CAX83         0 MADEEVQALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQ
                  0 |||||||||||||||||||||...........|...............|.||||||||||
AI021773.         0 MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQ

emb|CAX83        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE 108
                 60 |||||||||||||||||||||||||||||||||||||||||||||||| 108
AI021773.        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE 108

""",
        )
        hit = record[5]
        self.assertEqual(hit.num, 6)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "emb|VDP83060.1|")
        self.assertEqual(hit.target.name, "VDP83060")
        self.assertEqual(
            hit.target.description,
            "unnamed protein product, partial [Echinostoma caproni]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=209)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 407.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 161.384)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.49182e-47, places=52)
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
        self.assertEqual(hsp.target.id, "emb|VDP83060.1|")
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
Query : AI021773.1 Length: 108 Strand: Plus
        MAAD0534.RAR Schistosoma mansoni, adult worm (J.C.Parra) Schistosoma
        mansoni cDNA clone MAAD0534.RAR 5' end similar to S. mansoni actin mRNA,
        complete cds, mRNA sequence
Target: emb|VDP83060.1| Length: 209 Strand: Plus
        unnamed protein product, partial [Echinostoma caproni]

Score:161 bits(407), Expect:4e-47,
Identities:80/108(74%),  Positives:83/108(77%),  Gaps:0.108(0%)

emb|VDP83         0 MADDEVQALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQ
                  0 |||.|||||||||||||||||...........|...............|.||||||||||
AI021773.         0 MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQ

emb|VDP83        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE 108
                 60 |||||||||||||||||||||||||||||||||||||||||||||||| 108
AI021773.        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE 108

""",
        )
        hit = record[6]
        self.assertEqual(hit.num, 7)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "emb|CAA50205.1|")
        self.assertEqual(hit.target.name, "CAA50205")
        self.assertEqual(
            hit.target.description, "actin, partial [Entamoeba histolytica]"
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=137)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 401.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 159.073)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.85734e-47, places=52)
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
        self.assertEqual(hsp.target.id, "emb|CAA50205.1|")
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
Query : AI021773.1 Length: 108 Strand: Plus
        MAAD0534.RAR Schistosoma mansoni, adult worm (J.C.Parra) Schistosoma
        mansoni cDNA clone MAAD0534.RAR 5' end similar to S. mansoni actin mRNA,
        complete cds, mRNA sequence
Target: emb|CAA50205.1| Length: 137 Strand: Plus
        actin, partial [Entamoeba histolytica]

Score:159 bits(401), Expect:5e-47,
Identities:78/108(72%),  Positives:81/108(75%),  Gaps:0.108(0%)

emb|CAA50         0 MGDEEVQALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHVSVMAGMGQKDAYVGDEAQ
                  0 |.|||||||||||||||||||...........|...............|.||.|||||||
AI021773.         0 MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQ

emb|CAA50        60 SKRGILTLKYPIEHGIVNNWDDMEKIWHHTFYNELRVAPEEHPVLLTE 108
                 60 |||||||||||||||||.|||||||||||||||||||||||||||||| 108
AI021773.        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE 108

""",
        )
        hit = record[7]
        self.assertEqual(hit.num, 8)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "emb|VDN44756.1|")
        self.assertEqual(hit.target.name, "VDN44756")
        self.assertEqual(
            hit.target.description,
            "unnamed protein product, partial [Dibothriocephalus latus]",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=145)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 400.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 158.688)
        self.assertAlmostEqual(hsp.annotations["evalue"], 6.88203e-47, places=52)
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
        self.assertEqual(hsp.target.id, "emb|VDN44756.1|")
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
Query : AI021773.1 Length: 108 Strand: Plus
        MAAD0534.RAR Schistosoma mansoni, adult worm (J.C.Parra) Schistosoma
        mansoni cDNA clone MAAD0534.RAR 5' end similar to S. mansoni actin mRNA,
        complete cds, mRNA sequence
Target: emb|VDN44756.1| Length: 145 Strand: Plus
        unnamed protein product, partial [Dibothriocephalus latus]

Score:158 bits(400), Expect:7e-47,
Identities:78/108(72%),  Positives:82/108(76%),  Gaps:0.108(0%)

emb|VDN44         0 MGDEDVQALVIDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQ
                  0 |.||.|||||.||||||||||...........|...............|.||||||||||
AI021773.         0 MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQ

emb|VDN44        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE 108
                 60 |||||||||||||||||||||||||||||||||||||||||||||||| 108
AI021773.        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE 108

""",
        )
        hit = record[8]
        self.assertEqual(hit.num, 9)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "ref|XP_027046469.1|")
        self.assertEqual(hit.target.name, "XP_027046469")
        self.assertEqual(
            hit.target.description, "actin-1, partial [Pocillopora damicornis]"
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=122)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 398.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 157.918)
        self.assertAlmostEqual(hsp.annotations["evalue"], 7.4607e-47, places=52)
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
        self.assertEqual(hsp.target.id, "ref|XP_027046469.1|")
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
Query : AI021773.1 Length: 108 Strand: Plus
        MAAD0534.RAR Schistosoma mansoni, adult worm (J.C.Parra) Schistosoma
        mansoni cDNA clone MAAD0534.RAR 5' end similar to S. mansoni actin mRNA,
        complete cds, mRNA sequence
Target: ref|XP_027046469.1| Length: 122 Strand: Plus
        actin-1, partial [Pocillopora damicornis]

Score:157 bits(398), Expect:7e-47,
Identities:78/108(72%),  Positives:82/108(76%),  Gaps:0.108(0%)

ref|XP_02         0 MADEEVAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQ
                  0 ||||||.||||||||||||||...........|...............|.||||||||||
AI021773.         0 MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQ

ref|XP_02        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRIAPEEHPILLTE 108
                 60 ||||||||||||||||||||||||||||||||||||.||||||.|||| 108
AI021773.        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE 108

""",
        )
        hit = record[9]
        self.assertEqual(hit.num, 10)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "ref|XP_027046487.1|")
        self.assertEqual(hit.target.name, "XP_027046487")
        self.assertEqual(
            hit.target.description, "actin-1-like [Pocillopora damicornis]"
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=134)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 399.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 158.303)
        self.assertAlmostEqual(hsp.annotations["evalue"], 9.11071e-47, places=52)
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
        self.assertEqual(hsp.target.id, "ref|XP_027046487.1|")
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
Query : AI021773.1 Length: 108 Strand: Plus
        MAAD0534.RAR Schistosoma mansoni, adult worm (J.C.Parra) Schistosoma
        mansoni cDNA clone MAAD0534.RAR 5' end similar to S. mansoni actin mRNA,
        complete cds, mRNA sequence
Target: ref|XP_027046487.1| Length: 134 Strand: Plus
        actin-1-like [Pocillopora damicornis]

Score:158 bits(399), Expect:9e-47,
Identities:79/108(73%),  Positives:82/108(76%),  Gaps:0.108(0%)

ref|XP_02         0 MADEDVAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQ
                  0 ||||.|.||||||||||||||...........|...............|.||||||||||
AI021773.         0 MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQ

ref|XP_02        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE 108
                 60 |||||||||||||||||||||||||||||||||||||||||||||||| 108
AI021773.        60 SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE 108

""",
        )

    def test_xml2_21500_blastx_001_parser(self):
        """Parsing BLASTX 2.15.0+ (xml2_21500_blastx_001.xml)."""
        filename = "xml2_21500_blastx_001.xml"
        path = os.path.join("Blast", filename)
        with open(path, "rb") as stream:
            records = Blast.parse(stream)
            self.check_xml_21500_blastx_001_records(records, xml2=True)
        with Blast.parse(path) as records:
            self.check_xml_21500_blastx_001_records(records, xml2=True)
        with open(path, "rb") as stream:
            record = Blast.read(stream)
        self.check_xml_21500_blastx_001_record(record, xml2=True)
        record = Blast.read(path)
        self.check_xml_21500_blastx_001_record(record, xml2=True)

    def test_xml_21500_blastx_001_writer(self):
        """Writing BLASTX 2.15.0+ (xml_21500_blastx_001.xml)."""
        filename = "xml_21500_blastx_001.xml"
        path = os.path.join("Blast", filename)
        with Blast.parse(path) as records:
            stream = io.BytesIO()
            n = Blast.write(records, stream)
            self.assertEqual(n, 1)
            stream.seek(0)
            written_records = Blast.parse(stream)
            self.check_xml_21500_blastx_001_records(written_records, xml2=False)

    def test_xml2_21500_blastx_001_writer(self):
        """Writing BLASTX 2.15.0+ XML2 (xml2_21500_blastx_001.xml)."""
        filename = "xml2_21500_blastx_001.xml"
        path = os.path.join("Blast", filename)
        with Blast.parse(path) as records:
            stream = io.BytesIO()
            n = Blast.write(records, stream, fmt="XML2")
            self.assertEqual(n, 1)
            stream.seek(0)
            written_records = Blast.parse(stream)
            self.check_xml_21500_blastx_001_records(written_records, xml2=True)


class TestTBlastn(unittest.TestCase):
    """Test the Blast XML parser for tblastn output."""

    def test_xml_21500_tblastn_001_parser(self):
        """Parsing TBLASTN 2.15.0+ (xml_21500_tblastn_001.xml)."""
        filename = "xml_21500_tblastn_001.xml"
        path = os.path.join("Blast", filename)
        with open(path, "rb") as stream:
            records = Blast.parse(stream)
            self.check_xml_21500_tblastn_001_records(records)
        with Blast.parse(path) as records:
            self.check_xml_21500_tblastn_001_records(records)
        with open(path, "rb") as stream:
            record = Blast.read(stream)
        self.check_xml_21500_tblastn_001_record(record, xml2=False)
        record = Blast.read(path)
        self.check_xml_21500_tblastn_001_record(record, xml2=False)
        with Blast.parse(path) as records:
            self.assertEqual(
                str(records),
                """\
Program: TBLASTN 2.15.0+
     db: nt

  Query: CAJ99216.1 (length=234)
         tim [Helicobacter acinonychis str. Sheeba]
   Hits: ----  -----  ----------------------------------------------------------
            #  # HSP  ID + description
         ----  -----  ----------------------------------------------------------
            0      1  gi|109713861|emb|AM260522.1|  Helicobacter acinonychis ...
            1      1  gi|1336928286|emb|LT900055.1|  Helicobacter acinonychis...
            2      1  gi|1033012332|gb|CP011486.1|  Helicobacter pylori strai...
            3      1  gi|2641533851|gb|CP078169.1|  Helicobacter pylori strai...
            4      1  gi|2641529532|gb|CP078166.1|  Helicobacter pylori strai...
            5      1  gi|1033009499|gb|CP011484.1|  Helicobacter pylori strai...
            6      1  gi|317010283|gb|CP002336.1|  Helicobacter pylori SouthA...
            7      1  gi|2641538240|gb|CP078172.1|  Helicobacter pylori strai...
            8      1  gi|2640367186|gb|CP079077.1|  Helicobacter pylori strai...
            9      1  gi|532105813|gb|CP006691.1|  Helicobacter pylori SouthA...""",
            )

    def test_xml2_21500_tblastn_001_parser(self):
        """Parsing TBLASTN 2.15.0+ (xml2_21500_tblastn_001.xml)."""
        filename = "xml2_21500_tblastn_001.xml"
        path = os.path.join("Blast", filename)
        with open(path, "rb") as stream:
            records = Blast.parse(stream)
            self.check_xml_21500_tblastn_001_records(records, xml2=True)
        with Blast.parse(path) as records:
            self.check_xml_21500_tblastn_001_records(records, xml2=True)
        with open(path, "rb") as stream:
            record = Blast.read(stream)
        self.check_xml_21500_tblastn_001_record(record, xml2=True)
        record = Blast.read(path)
        self.check_xml_21500_tblastn_001_record(record, xml2=True)

    def check_xml_21500_tblastn_001_records(self, records, xml2=False):
        self.assertEqual(records.program, "tblastn")
        self.assertEqual(records.version, "TBLASTN 2.15.0+")
        self.assertEqual(
            records.reference,
            'Stephen F. Altschul, Thomas L. Madden, Alejandro A. Schäffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.',
        )
        self.assertEqual(records.db, "nt")
        if xml2:
            self.assertEqual(len(records.param), 7)
        else:
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
        self.assertEqual(records.param["filter"], "L;")
        if xml2:
            self.assertEqual(records.param["cbs"], 2)
            self.assertEqual(records.param["db-gencode"], 1)
        record = next(records)
        self.assertRaises(StopIteration, next, records)
        self.check_xml_21500_tblastn_001_record(record, xml2)

    def check_xml_21500_tblastn_001_record(self, record, xml2):
        if not xml2:
            self.assertEqual(record.num, 1)
        self.assertIsInstance(record.query, SeqRecord)
        self.assertEqual(record.query.id, "CAJ99216.1")
        self.assertEqual(
            record.query.description, "tim [Helicobacter acinonychis str. Sheeba]"
        )
        self.assertEqual(repr(record.query.seq), "Seq(None, length=234)")
        if xml2:
            self.assertEqual(len(record.query.features), 1)
            feature = record.query.features[0]
            self.assertEqual(feature.type, "masking")
            location = feature.location
            self.assertEqual(
                repr(location), "SimpleLocation(ExactPosition(101), ExactPosition(116))"
            )
        self.assertEqual(len(record.stat), 7)
        self.assertEqual(record.stat["db-num"], 104338276)
        self.assertEqual(record.stat["db-len"], 348058846)
        if xml2:
            self.assertEqual(record.stat["hsp-len"], 168)
            self.assertEqual(record.stat["eff-space"], 34567702523838)
        else:
            self.assertEqual(record.stat["hsp-len"], 0)
            self.assertEqual(record.stat["eff-space"], 0)
        self.assertAlmostEqual(record.stat["kappa"], 0.041)
        self.assertAlmostEqual(record.stat["lambda"], 0.267)
        self.assertAlmostEqual(record.stat["entropy"], 0.14)
        self.assertEqual(len(record), 10)

        hit = record[0]
        self.assertEqual(hit.num, 1)
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
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 1137.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 442.58)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.28996e-138, places=143)
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
        if xml2:
            # xml2 is True
            self.assertEqual(
                str(hsp),
                """\
Query : CAJ99216.1 Length: 234 Strand: Plus
        tim [Helicobacter acinonychis str. Sheeba]
Target: gi|109713861|emb|AM260522.1| Length: 234 Strand: Plus
        Helicobacter acinonychis str. Sheeba complete genome, strain Sheeba

Score:442 bits(1137), Expect:1e-138,
Identities:234/234(100%),  Positives:234/234(100%),  Gaps:0.234(0%)

gi|109713         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQN
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
CAJ99216.         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQN

gi|109713        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPNFLKEKFDFFKDKKFKIVY
                 60 ||||||||||||||||||||||||||||||||||||||||||...............|||
CAJ99216.        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPNflkekfdffkdkkfkIVY

gi|109713       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHG
                120 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
CAJ99216.       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHG

gi|109713       180 FLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL 234
                180 |||||||||||||||||||||||||||||||||||||||||||||||||||||| 234
CAJ99216.       180 FLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL 234

""",
            )
        else:
            # xml2 is False
            self.assertEqual(
                str(hsp),
                """\
Query : CAJ99216.1 Length: 234 Strand: Plus
        tim [Helicobacter acinonychis str. Sheeba]
Target: gi|109713861|emb|AM260522.1| Length: 234 Strand: Plus
        Helicobacter acinonychis str. Sheeba complete genome, strain Sheeba

Score:442 bits(1137), Expect:1e-138,
Identities:234/234(100%),  Positives:234/234(100%),  Gaps:0.234(0%)

gi|109713         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQN
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
CAJ99216.         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQN

gi|109713        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPNFLKEKFDFFKDKKFKIVY
                 60 ||||||||||||||||||||||||||||||||||||||||||...............|||
CAJ99216.        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPNXXXXXXXXXXXXXXXIVY

gi|109713       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHG
                120 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
CAJ99216.       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHG

gi|109713       180 FLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL 234
                180 |||||||||||||||||||||||||||||||||||||||||||||||||||||| 234
CAJ99216.       180 FLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL 234

""",
            )
        hit = record[1]
        self.assertEqual(hit.num, 2)
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
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 1128.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 439.113)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.05117e-137, places=142)
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
        if xml2:
            # xml2 is True
            self.assertEqual(
                str(hsp),
                """\
Query : CAJ99216.1 Length: 234 Strand: Plus
        tim [Helicobacter acinonychis str. Sheeba]
Target: gi|1336928286|emb|LT900055.1| Length: 234 Strand: Plus
        Helicobacter acinonychis isolate 212_9 genome assembly, chromosome: I

Score:439 bits(1128), Expect:2e-137,
Identities:232/234(99%),  Positives:232/234(99%),  Gaps:0.234(0%)

gi|133692         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQN
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
CAJ99216.         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQN

gi|133692        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRALLKESPNFLKEKFDFFKDKKFKIVY
                 60 ||||||||||||||||||||||||||||||||||.|||||||...............|||
CAJ99216.        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPNflkekfdffkdkkfkIVY

gi|133692       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHG
                120 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
CAJ99216.       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHG

gi|133692       180 FLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSASLELENFKTIISFL 234
                180 |||||||||||||||||||||||||||||||||||||||.|||||||||||||| 234
CAJ99216.       180 FLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL 234

""",
            )
        else:
            # xml2 is False
            self.assertEqual(
                str(hsp),
                """\
Query : CAJ99216.1 Length: 234 Strand: Plus
        tim [Helicobacter acinonychis str. Sheeba]
Target: gi|1336928286|emb|LT900055.1| Length: 234 Strand: Plus
        Helicobacter acinonychis isolate 212_9 genome assembly, chromosome: I

Score:439 bits(1128), Expect:2e-137,
Identities:232/234(99%),  Positives:232/234(99%),  Gaps:0.234(0%)

gi|133692         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQN
                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
CAJ99216.         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQN

gi|133692        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRALLKESPNFLKEKFDFFKDKKFKIVY
                 60 ||||||||||||||||||||||||||||||||||.|||||||...............|||
CAJ99216.        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPNXXXXXXXXXXXXXXXIVY

gi|133692       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHG
                120 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
CAJ99216.       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHG

gi|133692       180 FLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSASLELENFKTIISFL 234
                180 |||||||||||||||||||||||||||||||||||||||.|||||||||||||| 234
CAJ99216.       180 FLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL 234

""",
            )
        hit = record[2]
        self.assertEqual(hit.num, 3)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|1033012332|gb|CP011486.1|")
        self.assertEqual(hit.target.name, "CP011486")
        self.assertEqual(
            hit.target.description, "Helicobacter pylori strain K26A1, complete genome"
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=1570310)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 1076.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 419.083)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.02413e-130, places=135)
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
        if xml2:
            # xml2 is True
            self.assertEqual(
                str(hsp),
                """\
Query : CAJ99216.1 Length: 234 Strand: Plus
        tim [Helicobacter acinonychis str. Sheeba]
Target: gi|1033012332|gb|CP011486.1| Length: 234 Strand: Plus
        Helicobacter pylori strain K26A1, complete genome

Score:419 bits(1076), Expect:2e-130,
Identities:221/234(94%),  Positives:224/234(96%),  Gaps:0.234(0%)

gi|103301         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGAQN
                  0 |||||||||||||||||||||||||||||||||||||||||||||||||||||||||.||
CAJ99216.         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQN

gi|103301        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRTLLKESPSFLKEKFDFFKDKNFKIIY
                 60 ||||||||||||||||||||||||||||||||||.||||||................|.|
CAJ99216.        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPNflkekfdffkdkkfkIVY

gi|103301       120 CIGEDLKTREKGLAAVKEFLNEQLENIDLSYHNLIVAYEPIWAIGTKKSASLEDIYLTHG
                120 |||||||||||||.|||||||||||||||.|.||||||||||||||.|||||||||||||
CAJ99216.       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHG

gi|103301       180 FLKQILNQKTPLLYGGSVNTQNAKEILGIDSVDGLLVGSASLELENFKTIISFL 234
                180 ||||.||||.||||||||||||||||||||||||||.||.|||||||||||||| 234
CAJ99216.       180 FLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL 234

""",
            )
        else:
            # xml2 is False
            self.assertEqual(
                str(hsp),
                """\
Query : CAJ99216.1 Length: 234 Strand: Plus
        tim [Helicobacter acinonychis str. Sheeba]
Target: gi|1033012332|gb|CP011486.1| Length: 234 Strand: Plus
        Helicobacter pylori strain K26A1, complete genome

Score:419 bits(1076), Expect:2e-130,
Identities:221/234(94%),  Positives:224/234(96%),  Gaps:0.234(0%)

gi|103301         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGAQN
                  0 |||||||||||||||||||||||||||||||||||||||||||||||||||||||||.||
CAJ99216.         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQN

gi|103301        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRTLLKESPSFLKEKFDFFKDKNFKIIY
                 60 ||||||||||||||||||||||||||||||||||.||||||................|.|
CAJ99216.        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPNXXXXXXXXXXXXXXXIVY

gi|103301       120 CIGEDLKTREKGLAAVKEFLNEQLENIDLSYHNLIVAYEPIWAIGTKKSASLEDIYLTHG
                120 |||||||||||||.|||||||||||||||.|.||||||||||||||.|||||||||||||
CAJ99216.       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHG

gi|103301       180 FLKQILNQKTPLLYGGSVNTQNAKEILGIDSVDGLLVGSASLELENFKTIISFL 234
                180 ||||.||||.||||||||||||||||||||||||||.||.|||||||||||||| 234
CAJ99216.       180 FLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL 234

""",
            )
        hit = record[3]
        self.assertEqual(hit.num, 4)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|2641533851|gb|CP078169.1|")
        self.assertEqual(hit.target.name, "CP078169")
        self.assertEqual(
            hit.target.description,
            "Helicobacter pylori strain HpGP-ZAF-006 chromosome, complete genome",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=1641225)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 1069.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 416.387)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.58851e-129, places=134)
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
            "Seq('MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDKVFVFPDFLGLLPNAFLHF...SFL')",
        )
        self.assertEqual(hsp.target.id, "gi|2641533851|gb|CP078169.1|")
        self.assertEqual(hsp.target.name, "CP078169")
        self.assertEqual(
            hsp.target.description,
            "Helicobacter pylori strain HpGP-ZAF-006 chromosome, complete genome",
        )
        self.assertEqual(len(hsp.target.features), 1)
        feature = hsp.target.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(1641225))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(
            feature.qualifiers["coded_by"],
            "gi|2641533851|gb|CP078169.1|:260420..261121",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCD+VFVFPDFLGLLPNAFLHFTLG QNAYPKDCGAFTGEITSKHLEELKI+TLLIGHSERR LLKESP+FLKEKFDFFKDK FKIVYCIGEDLKTREKGLGAVKEFLNEQLENIDL Y NLIVAYEPIWAIGT KSASLEDIYLTHGFLKQ LNQK PLLYGGSVN QNAKEILGIDSVDGLLIGS SLELENFKTIISFL",
        )
        if xml2:
            # xml2 is True
            self.assertEqual(
                str(hsp),
                """\
Query : CAJ99216.1 Length: 234 Strand: Plus
        tim [Helicobacter acinonychis str. Sheeba]
Target: gi|2641533851|gb|CP078169.1| Length: 234 Strand: Plus
        Helicobacter pylori strain HpGP-ZAF-006 chromosome, complete genome

Score:416 bits(1069), Expect:2e-129,
Identities:221/234(94%),  Positives:224/234(96%),  Gaps:0.234(0%)

gi|264153         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDKVFVFPDFLGLLPNAFLHFTLGAQN
                  0 |||||||||||||||||||||||||||||||||||.|||||||||||||||||||||.||
CAJ99216.         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQN

gi|264153        60 AYPKDCGAFTGEITSKHLEELKIHTLLIGHSERRALLKESPSFLKEKFDFFKDKNFKIVY
                 60 |||||||||||||||||||||||.||||||||||.||||||................|||
CAJ99216.        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPNflkekfdffkdkkfkIVY

gi|264153       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLSYHNLIVAYEPIWAIGTKKSASLEDIYLTHG
                120 |||||||||||||||||||||||||||||.|.||||||||||||||.|||||||||||||
CAJ99216.       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHG

gi|264153       180 FLKQILNQKTPLLYGGSVNAQNAKEILGIDSVDGLLIGSASLELENFKTIISFL 234
                180 ||||.||||.|||||||||.|||||||||||||||||||.|||||||||||||| 234
CAJ99216.       180 FLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL 234

""",
            )
        else:
            # xml2 is False
            self.assertEqual(
                str(hsp),
                """\
Query : CAJ99216.1 Length: 234 Strand: Plus
        tim [Helicobacter acinonychis str. Sheeba]
Target: gi|2641533851|gb|CP078169.1| Length: 234 Strand: Plus
        Helicobacter pylori strain HpGP-ZAF-006 chromosome, complete genome

Score:416 bits(1069), Expect:2e-129,
Identities:221/234(94%),  Positives:224/234(96%),  Gaps:0.234(0%)

gi|264153         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDKVFVFPDFLGLLPNAFLHFTLGAQN
                  0 |||||||||||||||||||||||||||||||||||.|||||||||||||||||||||.||
CAJ99216.         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQN

gi|264153        60 AYPKDCGAFTGEITSKHLEELKIHTLLIGHSERRALLKESPSFLKEKFDFFKDKNFKIVY
                 60 |||||||||||||||||||||||.||||||||||.||||||................|||
CAJ99216.        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPNXXXXXXXXXXXXXXXIVY

gi|264153       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLSYHNLIVAYEPIWAIGTKKSASLEDIYLTHG
                120 |||||||||||||||||||||||||||||.|.||||||||||||||.|||||||||||||
CAJ99216.       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHG

gi|264153       180 FLKQILNQKTPLLYGGSVNAQNAKEILGIDSVDGLLIGSASLELENFKTIISFL 234
                180 ||||.||||.|||||||||.|||||||||||||||||||.|||||||||||||| 234
CAJ99216.       180 FLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL 234

""",
            )
        hit = record[4]
        self.assertEqual(hit.num, 5)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|2641529532|gb|CP078166.1|")
        self.assertEqual(hit.target.name, "CP078166")
        self.assertEqual(
            hit.target.description,
            "Helicobacter pylori strain HpGP-ZAF-009 chromosome, complete genome",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=1681930)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 1067.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 415.616)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.21842e-129, places=134)
        self.assertEqual(hsp.annotations["identity"], 220)
        self.assertEqual(hsp.annotations["positive"], 226)
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
            "Seq('MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDKVFVFPDFLGLLPNAFLHF...SFL')",
        )
        self.assertEqual(hsp.target.id, "gi|2641529532|gb|CP078166.1|")
        self.assertEqual(hsp.target.name, "CP078166")
        self.assertEqual(
            hsp.target.description,
            "Helicobacter pylori strain HpGP-ZAF-009 chromosome, complete genome",
        )
        self.assertEqual(len(hsp.target.features), 1)
        feature = hsp.target.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(1681930))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(
            feature.qualifiers["coded_by"],
            "gi|2641529532|gb|CP078166.1|:198327..199028",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCD+VFVFPDFLGLLPNAFLHFTLG QNAYPKDCGAFTGEITSKHLEELKI+TLLIGHSERR LLKESP+FLKEKFDFFKDKKFKI+YCIGEDLKTREKGLGAVKEFLNEQLENIDL Y +LIVAYEPIWAIGT KSASLEDIYLTHGFLKQ LNQK PLLYGGSVNTQNAKEILGIDSVDGLL+GS SLELENFKTIISFL",
        )
        if xml2:
            # xml2 is True
            self.assertEqual(
                str(hsp),
                """\
Query : CAJ99216.1 Length: 234 Strand: Plus
        tim [Helicobacter acinonychis str. Sheeba]
Target: gi|2641529532|gb|CP078166.1| Length: 234 Strand: Plus
        Helicobacter pylori strain HpGP-ZAF-009 chromosome, complete genome

Score:415 bits(1067), Expect:3e-129,
Identities:220/234(94%),  Positives:226/234(97%),  Gaps:0.234(0%)

gi|264152         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDKVFVFPDFLGLLPNAFLHFTLGAQN
                  0 |||||||||||||||||||||||||||||||||||.|||||||||||||||||||||.||
CAJ99216.         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQN

gi|264152        60 AYPKDCGAFTGEITSKHLEELKIHTLLIGHSERRTLLKESPSFLKEKFDFFKDKKFKIIY
                 60 |||||||||||||||||||||||.||||||||||.||||||................|.|
CAJ99216.        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPNflkekfdffkdkkfkIVY

gi|264152       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLSYHHLIVAYEPIWAIGTKKSASLEDIYLTHG
                120 |||||||||||||||||||||||||||||.|..|||||||||||||.|||||||||||||
CAJ99216.       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHG

gi|264152       180 FLKQILNQKTPLLYGGSVNTQNAKEILGIDSVDGLLVGSASLELENFKTIISFL 234
                180 ||||.||||.||||||||||||||||||||||||||.||.|||||||||||||| 234
CAJ99216.       180 FLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL 234

""",
            )
        else:
            # xml2 is False
            self.assertEqual(
                str(hsp),
                """\
Query : CAJ99216.1 Length: 234 Strand: Plus
        tim [Helicobacter acinonychis str. Sheeba]
Target: gi|2641529532|gb|CP078166.1| Length: 234 Strand: Plus
        Helicobacter pylori strain HpGP-ZAF-009 chromosome, complete genome

Score:415 bits(1067), Expect:3e-129,
Identities:220/234(94%),  Positives:226/234(97%),  Gaps:0.234(0%)

gi|264152         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDKVFVFPDFLGLLPNAFLHFTLGAQN
                  0 |||||||||||||||||||||||||||||||||||.|||||||||||||||||||||.||
CAJ99216.         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQN

gi|264152        60 AYPKDCGAFTGEITSKHLEELKIHTLLIGHSERRTLLKESPSFLKEKFDFFKDKKFKIIY
                 60 |||||||||||||||||||||||.||||||||||.||||||................|.|
CAJ99216.        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPNXXXXXXXXXXXXXXXIVY

gi|264152       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLSYHHLIVAYEPIWAIGTKKSASLEDIYLTHG
                120 |||||||||||||||||||||||||||||.|..|||||||||||||.|||||||||||||
CAJ99216.       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHG

gi|264152       180 FLKQILNQKTPLLYGGSVNTQNAKEILGIDSVDGLLVGSASLELENFKTIISFL 234
                180 ||||.||||.||||||||||||||||||||||||||.||.|||||||||||||| 234
CAJ99216.       180 FLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL 234

""",
            )
        hit = record[5]
        self.assertEqual(hit.num, 6)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|1033009499|gb|CP011484.1|")
        self.assertEqual(hit.target.name, "CP011484")
        self.assertEqual(
            hit.target.description, "Helicobacter pylori strain CC33C, complete genome"
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=1659899)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 1066.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 415.231)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.05935e-129, places=134)
        self.assertEqual(hsp.annotations["identity"], 220)
        self.assertEqual(hsp.annotations["positive"], 227)
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
            "Seq('MTKIAMANFKSSMPIFKSHAYLKELEKTLKPQHCDKVFVFPDFLGLLPNAFLHF...SFL')",
        )
        self.assertEqual(hsp.target.id, "gi|1033009499|gb|CP011484.1|")
        self.assertEqual(hsp.target.name, "CP011484")
        self.assertEqual(
            hsp.target.description, "Helicobacter pylori strain CC33C, complete genome"
        )
        self.assertEqual(len(hsp.target.features), 1)
        feature = hsp.target.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(1659899))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(
            feature.qualifiers["coded_by"],
            "gi|1033009499|gb|CP011484.1|:248047..248748",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MTKIAMANFKS+MPIFKSHAYLKELEKTLKPQHCD+VFVFPDFLGLLPNAFLHFTLG QNAYPKDCGAFTGEITSKHLEELKI+TLLIGHSERRVLLKESP+FLKEKFDFFKDKKFKIVYCIGEDLKTREKGLGAVKEFLNEQLENIDL+Y +LIVAYEPIWAIGT KSASLEDIYLTHGFLKQ LNQK PLLYGGSVN QNAKEILGIDSVDGLL+GS SLELENFKTIISFL",
        )
        if xml2:
            # xml2 is True
            self.assertEqual(
                str(hsp),
                """\
Query : CAJ99216.1 Length: 234 Strand: Plus
        tim [Helicobacter acinonychis str. Sheeba]
Target: gi|1033009499|gb|CP011484.1| Length: 234 Strand: Plus
        Helicobacter pylori strain CC33C, complete genome

Score:415 bits(1066), Expect:4e-129,
Identities:220/234(94%),  Positives:227/234(97%),  Gaps:0.234(0%)

gi|103300         0 MTKIAMANFKSSMPIFKSHAYLKELEKTLKPQHCDKVFVFPDFLGLLPNAFLHFTLGAQN
                  0 |||||||||||.|||||||||||||||||||||||.|||||||||||||||||||||.||
CAJ99216.         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQN

gi|103300        60 AYPKDCGAFTGEITSKHLEELKIHTLLIGHSERRVLLKESPSFLKEKFDFFKDKKFKIVY
                 60 |||||||||||||||||||||||.|||||||||||||||||................|||
CAJ99216.        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPNflkekfdffkdkkfkIVY

gi|103300       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLNYHHLIVAYEPIWAIGTKKSASLEDIYLTHG
                120 |||||||||||||||||||||||||||||.|..|||||||||||||.|||||||||||||
CAJ99216.       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHG

gi|103300       180 FLKQILNQKTPLLYGGSVNAQNAKEILGIDSVDGLLVGSASLELENFKTIISFL 234
                180 ||||.||||.|||||||||.||||||||||||||||.||.|||||||||||||| 234
CAJ99216.       180 FLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL 234

""",
            )
        else:
            # xml2 is False
            self.assertEqual(
                str(hsp),
                """\
Query : CAJ99216.1 Length: 234 Strand: Plus
        tim [Helicobacter acinonychis str. Sheeba]
Target: gi|1033009499|gb|CP011484.1| Length: 234 Strand: Plus
        Helicobacter pylori strain CC33C, complete genome

Score:415 bits(1066), Expect:4e-129,
Identities:220/234(94%),  Positives:227/234(97%),  Gaps:0.234(0%)

gi|103300         0 MTKIAMANFKSSMPIFKSHAYLKELEKTLKPQHCDKVFVFPDFLGLLPNAFLHFTLGAQN
                  0 |||||||||||.|||||||||||||||||||||||.|||||||||||||||||||||.||
CAJ99216.         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQN

gi|103300        60 AYPKDCGAFTGEITSKHLEELKIHTLLIGHSERRVLLKESPSFLKEKFDFFKDKKFKIVY
                 60 |||||||||||||||||||||||.|||||||||||||||||................|||
CAJ99216.        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPNXXXXXXXXXXXXXXXIVY

gi|103300       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLNYHHLIVAYEPIWAIGTKKSASLEDIYLTHG
                120 |||||||||||||||||||||||||||||.|..|||||||||||||.|||||||||||||
CAJ99216.       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHG

gi|103300       180 FLKQILNQKTPLLYGGSVNAQNAKEILGIDSVDGLLVGSASLELENFKTIISFL 234
                180 ||||.||||.|||||||||.||||||||||||||||.||.|||||||||||||| 234
CAJ99216.       180 FLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL 234

""",
            )
        hit = record[6]
        self.assertEqual(hit.num, 7)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|317010283|gb|CP002336.1|")
        self.assertEqual(hit.target.name, "CP002336")
        self.assertEqual(
            hit.target.description,
            "Helicobacter pylori SouthAfrica7, complete genome",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=1653913)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 1066.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 415.231)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.73879e-129, places=134)
        self.assertEqual(hsp.annotations["identity"], 221)
        self.assertEqual(hsp.annotations["positive"], 226)
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
            "Seq('MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDKVFVFPDFLGLLPNAFLHF...SFL')",
        )
        self.assertEqual(hsp.target.id, "gi|317010283|gb|CP002336.1|")
        self.assertEqual(hsp.target.name, "CP002336")
        self.assertEqual(
            hsp.target.description,
            "Helicobacter pylori SouthAfrica7, complete genome",
        )
        self.assertEqual(len(hsp.target.features), 1)
        feature = hsp.target.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(1653913))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(
            feature.qualifiers["coded_by"],
            "gi|317010283|gb|CP002336.1|:194539..195240",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCD+VFVFPDFLGLLPNAFLHFTLG QNAYPKD GAFTGEITSKHLEELKINTLLIGHSERRVLLKESP+FLKEKFDFFKDKKFKI+YCIGEDLKTREKGLGAVKEFLNEQLENIDL Y +LIVAYEPIWAIGT KSASLEDIYLTHGFLKQ LNQK PLLYGGSVNTQNAKEILGIDSVDGLL+GS SLELENFKTIISFL",
        )
        if xml2:
            # xml2 is True
            self.assertEqual(
                str(hsp),
                """\
Query : CAJ99216.1 Length: 234 Strand: Plus
        tim [Helicobacter acinonychis str. Sheeba]
Target: gi|317010283|gb|CP002336.1| Length: 234 Strand: Plus
        Helicobacter pylori SouthAfrica7, complete genome

Score:415 bits(1066), Expect:5e-129,
Identities:221/234(94%),  Positives:226/234(97%),  Gaps:0.234(0%)

gi|317010         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDKVFVFPDFLGLLPNAFLHFTLGAQN
                  0 |||||||||||||||||||||||||||||||||||.|||||||||||||||||||||.||
CAJ99216.         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQN

gi|317010        60 AYPKDYGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPSFLKEKFDFFKDKKFKIIY
                 60 |||||.|||||||||||||||||||||||||||||||||||................|.|
CAJ99216.        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPNflkekfdffkdkkfkIVY

gi|317010       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLSYHHLIVAYEPIWAIGTKKSASLEDIYLTHG
                120 |||||||||||||||||||||||||||||.|..|||||||||||||.|||||||||||||
CAJ99216.       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHG

gi|317010       180 FLKQILNQKTPLLYGGSVNTQNAKEILGIDSVDGLLVGSASLELENFKTIISFL 234
                180 ||||.||||.||||||||||||||||||||||||||.||.|||||||||||||| 234
CAJ99216.       180 FLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL 234

""",
            )
        else:
            # xml2 is False
            self.assertEqual(
                str(hsp),
                """\
Query : CAJ99216.1 Length: 234 Strand: Plus
        tim [Helicobacter acinonychis str. Sheeba]
Target: gi|317010283|gb|CP002336.1| Length: 234 Strand: Plus
        Helicobacter pylori SouthAfrica7, complete genome

Score:415 bits(1066), Expect:5e-129,
Identities:221/234(94%),  Positives:226/234(97%),  Gaps:0.234(0%)

gi|317010         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDKVFVFPDFLGLLPNAFLHFTLGAQN
                  0 |||||||||||||||||||||||||||||||||||.|||||||||||||||||||||.||
CAJ99216.         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQN

gi|317010        60 AYPKDYGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPSFLKEKFDFFKDKKFKIIY
                 60 |||||.|||||||||||||||||||||||||||||||||||................|.|
CAJ99216.        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPNXXXXXXXXXXXXXXXIVY

gi|317010       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLSYHHLIVAYEPIWAIGTKKSASLEDIYLTHG
                120 |||||||||||||||||||||||||||||.|..|||||||||||||.|||||||||||||
CAJ99216.       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHG

gi|317010       180 FLKQILNQKTPLLYGGSVNTQNAKEILGIDSVDGLLVGSASLELENFKTIISFL 234
                180 ||||.||||.||||||||||||||||||||||||||.||.|||||||||||||| 234
CAJ99216.       180 FLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL 234

""",
            )
        hit = record[7]
        self.assertEqual(hit.num, 8)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|2641538240|gb|CP078172.1|")
        self.assertEqual(hit.target.name, "CP078172")
        self.assertEqual(
            hit.target.description,
            "Helicobacter pylori strain HpGP-ZAF-001 chromosome, complete genome",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=1714499)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 1062.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 413.69)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.37323e-128, places=133)
        self.assertEqual(hsp.annotations["identity"], 219)
        self.assertEqual(hsp.annotations["positive"], 225)
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
            "Seq('MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDKVFVFPDFLGLLPNAFLHF...SFL')",
        )
        self.assertEqual(hsp.target.id, "gi|2641538240|gb|CP078172.1|")
        self.assertEqual(hsp.target.name, "CP078172")
        self.assertEqual(
            hsp.target.description,
            "Helicobacter pylori strain HpGP-ZAF-001 chromosome, complete genome",
        )
        self.assertEqual(len(hsp.target.features), 1)
        feature = hsp.target.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(1714499))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(
            feature.qualifiers["coded_by"],
            "complement(gi|2641538240|gb|CP078172.1|:1021670..1022371)",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCD+VFVFPDFLGLLPNAFLHFTLG QNAYPKDCGAFTGEITSKHLEELKI+TLLIGHSERR LLKESP+FLKEKFDFFKDKKFKI+YCIGEDLKTREKGLGAVKEFLNEQLENIDL Y +LIVAYEPIWAIGT KSASLEDIYLTHGFLKQ LNQK PLLYGGSVN QNAKEILGIDSVDGLL+GS SLELENFKTIISFL",
        )
        if xml2:
            # xml2 is True
            self.assertEqual(
                str(hsp),
                """\
Query : CAJ99216.1 Length: 234 Strand: Plus
        tim [Helicobacter acinonychis str. Sheeba]
Target: gi|2641538240|gb|CP078172.1| Length: 234 Strand: Plus
        Helicobacter pylori strain HpGP-ZAF-001 chromosome, complete genome

Score:413 bits(1062), Expect:1e-128,
Identities:219/234(94%),  Positives:225/234(96%),  Gaps:0.234(0%)

gi|264153         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDKVFVFPDFLGLLPNAFLHFTLGAQN
                  0 |||||||||||||||||||||||||||||||||||.|||||||||||||||||||||.||
CAJ99216.         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQN

gi|264153        60 AYPKDCGAFTGEITSKHLEELKIHTLLIGHSERRTLLKESPSFLKEKFDFFKDKKFKIIY
                 60 |||||||||||||||||||||||.||||||||||.||||||................|.|
CAJ99216.        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPNflkekfdffkdkkfkIVY

gi|264153       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLSYHHLIVAYEPIWAIGTKKSASLEDIYLTHG
                120 |||||||||||||||||||||||||||||.|..|||||||||||||.|||||||||||||
CAJ99216.       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHG

gi|264153       180 FLKQILNQKTPLLYGGSVNAQNAKEILGIDSVDGLLVGSASLELENFKTIISFL 234
                180 ||||.||||.|||||||||.||||||||||||||||.||.|||||||||||||| 234
CAJ99216.       180 FLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL 234

""",
            )
        else:
            # xml2 is False
            self.assertEqual(
                str(hsp),
                """\
Query : CAJ99216.1 Length: 234 Strand: Plus
        tim [Helicobacter acinonychis str. Sheeba]
Target: gi|2641538240|gb|CP078172.1| Length: 234 Strand: Plus
        Helicobacter pylori strain HpGP-ZAF-001 chromosome, complete genome

Score:413 bits(1062), Expect:1e-128,
Identities:219/234(94%),  Positives:225/234(96%),  Gaps:0.234(0%)

gi|264153         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDKVFVFPDFLGLLPNAFLHFTLGAQN
                  0 |||||||||||||||||||||||||||||||||||.|||||||||||||||||||||.||
CAJ99216.         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQN

gi|264153        60 AYPKDCGAFTGEITSKHLEELKIHTLLIGHSERRTLLKESPSFLKEKFDFFKDKKFKIIY
                 60 |||||||||||||||||||||||.||||||||||.||||||................|.|
CAJ99216.        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPNXXXXXXXXXXXXXXXIVY

gi|264153       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLSYHHLIVAYEPIWAIGTKKSASLEDIYLTHG
                120 |||||||||||||||||||||||||||||.|..|||||||||||||.|||||||||||||
CAJ99216.       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHG

gi|264153       180 FLKQILNQKTPLLYGGSVNAQNAKEILGIDSVDGLLVGSASLELENFKTIISFL 234
                180 ||||.||||.|||||||||.||||||||||||||||.||.|||||||||||||| 234
CAJ99216.       180 FLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL 234

""",
            )
        hit = record[8]
        self.assertEqual(hit.num, 9)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|2640367186|gb|CP079077.1|")
        self.assertEqual(hit.target.name, "CP079077")
        self.assertEqual(
            hit.target.description,
            "Helicobacter pylori strain HpGP-ARG-001 chromosome, complete genome",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=1624657)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 1059.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 412.535)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.5777e-128, places=132)
        self.assertEqual(hsp.annotations["identity"], 220)
        self.assertEqual(hsp.annotations["positive"], 225)
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
        self.assertEqual(hsp.target.id, "gi|2640367186|gb|CP079077.1|")
        self.assertEqual(hsp.target.name, "CP079077")
        self.assertEqual(
            hsp.target.description,
            "Helicobacter pylori strain HpGP-ARG-001 chromosome, complete genome",
        )
        self.assertEqual(len(hsp.target.features), 1)
        feature = hsp.target.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(1624657))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(
            feature.qualifiers["coded_by"],
            "gi|2640367186|gb|CP079077.1|:197481..198182",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQH DRVFVFPDFLGLLPNAFLHFTLGVQNAYP+DCGAFTGEITSKHLEELKI+TLLIGHSERRVLLKESP+FLKEKFDFFKDK FKIVYCIGEDLKTREKG  AVKEFL+EQLENIDL+Y NLIVAYEPIWAIGT KSASLEDIYLTHGFLKQ LNQK PLLYGGSVNTQNAKEILGIDSVDGLLIGS SLELENFKTIISFL",
        )
        if xml2:
            # xml2 is True
            self.assertEqual(
                str(hsp),
                """\
Query : CAJ99216.1 Length: 234 Strand: Plus
        tim [Helicobacter acinonychis str. Sheeba]
Target: gi|2640367186|gb|CP079077.1| Length: 234 Strand: Plus
        Helicobacter pylori strain HpGP-ARG-001 chromosome, complete genome

Score:412 bits(1059), Expect:4e-128,
Identities:220/234(94%),  Positives:225/234(96%),  Gaps:0.234(0%)

gi|264036         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHFDRVFVFPDFLGLLPNAFLHFTLGVQN
                  0 |||||||||||||||||||||||||||||||||.||||||||||||||||||||||||||
CAJ99216.         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQN

gi|264036        60 AYPRDCGAFTGEITSKHLEELKIHTLLIGHSERRVLLKESPSFLKEKFDFFKDKNFKIVY
                 60 |||.|||||||||||||||||||.|||||||||||||||||................|||
CAJ99216.        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPNflkekfdffkdkkfkIVY

gi|264036       120 CIGEDLKTREKGFKAVKEFLSEQLENIDLNYSNLIVAYEPIWAIGTKKSASLEDIYLTHG
                120 ||||||||||||..||||||.||||||||.|.||||||||||||||.|||||||||||||
CAJ99216.       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHG

gi|264036       180 FLKQILNQKTPLLYGGSVNTQNAKEILGIDSVDGLLIGSASLELENFKTIISFL 234
                180 ||||.||||.|||||||||||||||||||||||||||||.|||||||||||||| 234
CAJ99216.       180 FLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL 234

""",
            )
        else:
            self.assertEqual(
                str(hsp),
                """\
Query : CAJ99216.1 Length: 234 Strand: Plus
        tim [Helicobacter acinonychis str. Sheeba]
Target: gi|2640367186|gb|CP079077.1| Length: 234 Strand: Plus
        Helicobacter pylori strain HpGP-ARG-001 chromosome, complete genome

Score:412 bits(1059), Expect:4e-128,
Identities:220/234(94%),  Positives:225/234(96%),  Gaps:0.234(0%)

gi|264036         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHFDRVFVFPDFLGLLPNAFLHFTLGVQN
                  0 |||||||||||||||||||||||||||||||||.||||||||||||||||||||||||||
CAJ99216.         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQN

gi|264036        60 AYPRDCGAFTGEITSKHLEELKIHTLLIGHSERRVLLKESPSFLKEKFDFFKDKNFKIVY
                 60 |||.|||||||||||||||||||.|||||||||||||||||................|||
CAJ99216.        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPNXXXXXXXXXXXXXXXIVY

gi|264036       120 CIGEDLKTREKGFKAVKEFLSEQLENIDLNYSNLIVAYEPIWAIGTKKSASLEDIYLTHG
                120 ||||||||||||..||||||.||||||||.|.||||||||||||||.|||||||||||||
CAJ99216.       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHG

gi|264036       180 FLKQILNQKTPLLYGGSVNTQNAKEILGIDSVDGLLIGSASLELENFKTIISFL 234
                180 ||||.||||.|||||||||||||||||||||||||||||.|||||||||||||| 234
CAJ99216.       180 FLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL 234

""",
            )
        hit = record[9]
        self.assertEqual(hit.num, 10)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|532105813|gb|CP006691.1|")
        self.assertEqual(hit.target.name, "CP006691")
        self.assertEqual(
            hit.target.description, "Helicobacter pylori SouthAfrica20, complete genome"
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=1622903)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 1057.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 411.764)
        self.assertAlmostEqual(hsp.annotations["evalue"], 7.68167e-128, places=133)
        self.assertEqual(hsp.annotations["identity"], 218)
        self.assertEqual(hsp.annotations["positive"], 225)
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
            "Seq('MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDKVFVFPDFLGLLPNAFLHF...SFL')",
        )
        self.assertEqual(hsp.target.id, "gi|532105813|gb|CP006691.1|")
        self.assertEqual(hsp.target.name, "CP006691")
        self.assertEqual(
            hsp.target.description, "Helicobacter pylori SouthAfrica20, complete genome"
        )
        self.assertEqual(len(hsp.target.features), 1)
        feature = hsp.target.features[0]
        self.assertEqual(
            repr(feature.location),
            "SimpleLocation(ExactPosition(0), ExactPosition(1622903))",
        )
        self.assertEqual(feature.type, "CDS")
        self.assertEqual(len(feature.qualifiers), 1)
        self.assertEqual(
            feature.qualifiers["coded_by"],
            "gi|532105813|gb|CP006691.1|:197836..198537",
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCD+VFVFPDFLGLLPNAFLHFTLG QNAYPKDCGAFTGEITSKHLEELKI+TLLIGHSERRVLLKESP+FLKEKFDFFKDKKFKI+YCIGEDLKTREKG  AVKEFLNEQLENIDL+Y +LIVAYEPIWAIGT KSASLEDIYLTHGFLKQ LNQK PLLYGGSVN QNAKEILGIDSVDGLL+GS SLELENFKTIISFL",
        )
        if xml2:
            # xml2 is True
            self.assertEqual(
                str(hsp),
                """\
Query : CAJ99216.1 Length: 234 Strand: Plus
        tim [Helicobacter acinonychis str. Sheeba]
Target: gi|532105813|gb|CP006691.1| Length: 234 Strand: Plus
        Helicobacter pylori SouthAfrica20, complete genome

Score:411 bits(1057), Expect:8e-128,
Identities:218/234(93%),  Positives:225/234(96%),  Gaps:0.234(0%)

gi|532105         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDKVFVFPDFLGLLPNAFLHFTLGAQN
                  0 |||||||||||||||||||||||||||||||||||.|||||||||||||||||||||.||
CAJ99216.         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQN

gi|532105        60 AYPKDCGAFTGEITSKHLEELKIHTLLIGHSERRVLLKESPSFLKEKFDFFKDKKFKIIY
                 60 |||||||||||||||||||||||.|||||||||||||||||................|.|
CAJ99216.        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPNflkekfdffkdkkfkIVY

gi|532105       120 CIGEDLKTREKGFNAVKEFLNEQLENIDLNYHHLIVAYEPIWAIGTKKSASLEDIYLTHG
                120 ||||||||||||..|||||||||||||||.|..|||||||||||||.|||||||||||||
CAJ99216.       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHG

gi|532105       180 FLKQILNQKTPLLYGGSVNAQNAKEILGIDSVDGLLVGSASLELENFKTIISFL 234
                180 ||||.||||.|||||||||.||||||||||||||||.||.|||||||||||||| 234
CAJ99216.       180 FLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL 234

""",
            )
        else:
            self.assertEqual(
                str(hsp),
                """\
Query : CAJ99216.1 Length: 234 Strand: Plus
        tim [Helicobacter acinonychis str. Sheeba]
Target: gi|532105813|gb|CP006691.1| Length: 234 Strand: Plus
        Helicobacter pylori SouthAfrica20, complete genome

Score:411 bits(1057), Expect:8e-128,
Identities:218/234(93%),  Positives:225/234(96%),  Gaps:0.234(0%)

gi|532105         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDKVFVFPDFLGLLPNAFLHFTLGAQN
                  0 |||||||||||||||||||||||||||||||||||.|||||||||||||||||||||.||
CAJ99216.         0 MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQN

gi|532105        60 AYPKDCGAFTGEITSKHLEELKIHTLLIGHSERRVLLKESPSFLKEKFDFFKDKKFKIIY
                 60 |||||||||||||||||||||||.|||||||||||||||||................|.|
CAJ99216.        60 AYPKDCGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPNXXXXXXXXXXXXXXXIVY

gi|532105       120 CIGEDLKTREKGFNAVKEFLNEQLENIDLNYHHLIVAYEPIWAIGTKKSASLEDIYLTHG
                120 ||||||||||||..|||||||||||||||.|..|||||||||||||.|||||||||||||
CAJ99216.       120 CIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHG

gi|532105       180 FLKQILNQKTPLLYGGSVNAQNAKEILGIDSVDGLLVGSASLELENFKTIISFL 234
                180 ||||.||||.|||||||||.||||||||||||||||.||.|||||||||||||| 234
CAJ99216.       180 FLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL 234

""",
            )

    def test_xml_21500_tblastn_001_writer(self):
        """Writing TBLASTN 2.15.0+ (xml_21500_tblastn_001.xml)."""
        filename = "xml_21500_tblastn_001.xml"
        path = os.path.join("Blast", filename)
        with Blast.parse(path) as records:
            stream = io.BytesIO()
            n = Blast.write(records, stream)
            self.assertEqual(n, 1)
            stream.seek(0)
            written_records = Blast.parse(stream)
            self.check_xml_21500_tblastn_001_records(written_records)

    def test_xml2_21500_tblastn_001_writer(self):
        """Writing TBLASTN 2.15.0+ XML2 (xml2_21500_tblastn_001.xml)."""
        filename = "xml2_21500_tblastn_001.xml"
        path = os.path.join("Blast", filename)
        with Blast.parse(path) as records:
            stream = io.BytesIO()
            n = Blast.write(records, stream, fmt="XML2")
            self.assertEqual(n, 1)
            stream.seek(0)
            written_records = Blast.parse(stream)
            self.check_xml_21500_tblastn_001_records(written_records, xml2=True)


class TestTBlastx(unittest.TestCase):
    """Test the Blast XML parser for tblastx output."""

    def test_xml_2226_tblastx_004(self):
        """Parsing TBLASTX 2.2.26+ (xml_2226_tblastx_004.xml)."""
        filename = "xml_2226_tblastx_004.xml"
        path = os.path.join("Blast", filename)
        with open(path, "rb") as stream:
            records = Blast.parse(stream)
            self.check_xml_2226_tblastx_004(records)
        with Blast.parse(path) as records:
            self.check_xml_2226_tblastx_004(records)
        with Blast.parse(path) as records:
            self.assertEqual(
                str(records),
                """\
Program: TBLASTX 2.2.26+
     db: refseq_rna

  Query: Query_1 (length=128)
         random_s00
   Hits: No hits found

  Query: Query_2 (length=350)
         gi|296147483:1-350 Saccharomyces cerevisiae S288c Mon2p (MON2) mRNA,
         complete cds
   Hits: ----  -----  ----------------------------------------------------------
            #  # HSP  ID + description
         ----  -----  ----------------------------------------------------------
            0      7  gi|296147483|ref|NM_001183135.1|  Saccharomyces cerevis...
            1      6  gi|365982352|ref|XM_003667962.1|  Naumovozyma dairenens...
            2      4  gi|366988334|ref|XM_003673886.1|  Naumovozyma castellii...
            3      2  gi|255710474|ref|XM_002551475.1|  Lachancea thermotoler...
            4      4  gi|254579534|ref|XM_002495708.1|  Zygosaccharomyces rou...""",
            )

    def check_xml_2226_tblastx_004(self, records):
        self.assertEqual(records.program, "tblastx")
        self.assertEqual(records.version, "TBLASTX 2.2.26+")
        self.assertEqual(
            records.reference,
            'Stephen F. Altschul, Thomas L. Madden, Alejandro A. Schäffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.',
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
        self.assertEqual(record.num, 1)
        self.assertIsInstance(record.query, SeqRecord)
        self.assertEqual(record.query.id, "Query_1")
        self.assertEqual(record.query.description, "random_s00")
        self.assertEqual(repr(record.query.seq), "Seq(None, length=128)")

        self.assertEqual(len(record.stat), 7)
        self.assertEqual(record.stat["db-num"], 2933984)
        self.assertEqual(record.stat["db-len"], 4726730735)
        self.assertEqual(record.stat["hsp-len"], 0)
        self.assertEqual(record.stat["eff-space"], 0)
        self.assertAlmostEqual(record.stat["kappa"], 0.0)
        self.assertAlmostEqual(record.stat["lambda"], 0.0)
        self.assertAlmostEqual(record.stat["entropy"], 0.0)
        self.assertEqual(len(record), 0)
        record = next(records)
        self.assertEqual(record.num, 2)
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
        self.assertEqual(record.stat["eff-space"], 0)
        self.assertAlmostEqual(record.stat["kappa"], 0.0)
        self.assertAlmostEqual(record.stat["lambda"], 0.0)
        self.assertAlmostEqual(record.stat["entropy"], 0.0)
        self.assertEqual(len(record), 5)
        hit = record[0]
        self.assertEqual(hit.num, 1)
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
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 626.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 289.739)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.05874e-76, places=81)
        self.assertEqual(hsp.annotations["identity"], 116)
        self.assertEqual(hsp.annotations["positive"], 116)
        self.assertEqual(hsp.annotations["gaps"], 0)
        hsp = hit[1]
        self.assertEqual(hsp.num, 2)
        self.assertAlmostEqual(hsp.score, 602.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 278.742)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.16381e-73, places=78)
        self.assertEqual(hsp.annotations["identity"], 116)
        self.assertEqual(hsp.annotations["positive"], 116)
        self.assertEqual(hsp.annotations["gaps"], 0)
        hsp = hit[2]
        self.assertEqual(hsp.num, 3)
        self.assertAlmostEqual(hsp.score, 593.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 274.618)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.77251e-72, places=77)
        self.assertEqual(hsp.annotations["identity"], 116)
        self.assertEqual(hsp.annotations["positive"], 116)
        self.assertEqual(hsp.annotations["gaps"], 0)
        hsp = hit[3]
        self.assertEqual(hsp.num, 4)
        self.assertAlmostEqual(hsp.score, 583.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 270.036)
        self.assertAlmostEqual(hsp.annotations["evalue"], 9.03598e-71, places=76)
        self.assertEqual(hsp.annotations["identity"], 116)
        self.assertEqual(hsp.annotations["positive"], 116)
        self.assertEqual(hsp.annotations["gaps"], 0)
        hsp = hit[4]
        self.assertEqual(hsp.num, 5)
        self.assertAlmostEqual(hsp.score, 495.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 229.713)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.24226e-58, places=63)
        self.assertEqual(hsp.annotations["identity"], 116)
        self.assertEqual(hsp.annotations["positive"], 116)
        self.assertEqual(hsp.annotations["gaps"], 0)
        hsp = hit[5]
        self.assertEqual(hsp.num, 6)
        self.assertAlmostEqual(hsp.score, 425.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 197.639)
        self.assertAlmostEqual(hsp.annotations["evalue"], 9.12288e-54, places=59)
        self.assertEqual(hsp.annotations["identity"], 85)
        self.assertEqual(hsp.annotations["positive"], 85)
        self.assertEqual(hsp.annotations["gaps"], 0)
        hsp = hit[6]
        self.assertEqual(hsp.num, 7)
        self.assertAlmostEqual(hsp.score, 73.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 36.3494)
        self.assertAlmostEqual(hsp.annotations["evalue"], 9.12288e-54, places=59)
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
Query : Query_2 Length: 14 Strand: Plus
        gi|296147483:1-350 Saccharomyces cerevisiae S288c Mon2p (MON2) mRNA,
        complete cds
Target: gi|296147483|ref|NM_001183135.1| Length: 14 Strand: Plus
        Saccharomyces cerevisiae S288c Mon2p (MON2) mRNA, complete cds
        >gi|116616412|gb|EF059095.1| Synthetic construct Saccharomyces
        cerevisiae clone FLH203015.01X MON2, complete sequence

Score:36 bits(73), Expect:9e-54,
Identities:14/14(100%),  Positives:14/14(100%),  Gaps:0.14(0%)

gi|296147         0 MAMNTGGFDSMQRQ 14
                  0 |||||||||||||| 14
Query_2           0 MAMNTGGFDSMQRQ 14

""",
        )
        hit = record[1]
        self.assertEqual(hit.num, 2)
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
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 327.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 152.734)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.38069e-37, places=42)
        self.assertEqual(hsp.annotations["identity"], 62)
        self.assertEqual(hsp.annotations["positive"], 73)
        self.assertEqual(hsp.annotations["gaps"], 0)
        hsp = hit[1]
        self.assertEqual(hsp.num, 2)
        self.assertAlmostEqual(hsp.score, 51.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 26.2688)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.38069e-37, places=42)
        self.assertEqual(hsp.annotations["identity"], 11)
        self.assertEqual(hsp.annotations["positive"], 11)
        self.assertEqual(hsp.annotations["gaps"], 0)
        hsp = hit[2]
        self.assertEqual(hsp.num, 3)
        self.assertAlmostEqual(hsp.score, 142.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 67.9658)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.80116e-20, places=25)
        self.assertEqual(hsp.annotations["identity"], 34)
        self.assertEqual(hsp.annotations["positive"], 38)
        self.assertEqual(hsp.annotations["gaps"], 0)
        hsp = hit[3]
        self.assertEqual(hsp.num, 4)
        self.assertAlmostEqual(hsp.score, 109.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 52.8449)
        self.assertAlmostEqual(hsp.annotations["evalue"], 4.80116e-20, places=25)
        self.assertEqual(hsp.annotations["identity"], 24)
        self.assertEqual(hsp.annotations["positive"], 29)
        self.assertEqual(hsp.annotations["gaps"], 0)
        hsp = hit[4]
        self.assertEqual(hsp.num, 5)
        self.assertAlmostEqual(hsp.score, 127.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 61.0927)
        self.assertAlmostEqual(hsp.annotations["evalue"], 7.14684e-08, places=13)
        self.assertEqual(hsp.annotations["identity"], 36)
        self.assertEqual(hsp.annotations["positive"], 52)
        self.assertEqual(hsp.annotations["gaps"], 0)
        hsp = hit[5]
        self.assertEqual(hsp.num, 6)
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
Query : Query_2 Length: 65 Strand: Plus
        gi|296147483:1-350 Saccharomyces cerevisiae S288c Mon2p (MON2) mRNA,
        complete cds
Target: gi|365982352|ref|XM_003667962.1| Length: 65 Strand: Plus
        Naumovozyma dairenensis CBS 421 hypothetical protein (NDAI0A06120), mRNA

Score:42 bits(87), Expect:0.02,
Identities:28/65(43%),  Positives:36/65(55%),  Gaps:0.65(0%)

gi|365982         0 GRQSL*ALHGN*S*FGISR*TGQDQRCNKIRVPYKFFYILNSL*NID*FVASMFNG*ISS
                  0 |.|||||||.....|.|.....|..|...||..........|......|..||.||||||
Query_2           0 G*QSL*ALHCQGRHFSIP*LASQHERECEIRMSF*LLKTMYSFQYLNGFITSMANG*ISS

gi|365982        60 FRF*R 65
                 60 |||.| 65
Query_2          60 FRFGR 65

""",
        )
        hit = record[2]
        self.assertEqual(hit.num, 3)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|366988334|ref|XM_003673886.1|")
        self.assertEqual(hit.target.name, "XM_003673886")
        self.assertEqual(
            hit.target.description,
            "Naumovozyma castellii CBS 4309 hypothetical protein (NCAS0A09950) mRNA, complete cds",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=4938)")
        self.assertEqual(len(hit), 4)
        self.assertEqual(
            repr(hit),
            "<Bio.Blast.Hit target.id='gi|366988334|ref|XM_003673886.1|' query.id='Query_2'; 4 HSPs>",
        )
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 306.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 143.112)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.45826e-32, places=37)
        self.assertEqual(hsp.annotations["identity"], 58)
        self.assertEqual(hsp.annotations["positive"], 71)
        self.assertEqual(hsp.annotations["gaps"], 0)
        hsp = hit[1]
        self.assertEqual(hsp.num, 2)
        self.assertAlmostEqual(hsp.score, 130.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 62.4673)
        self.assertAlmostEqual(hsp.annotations["evalue"], 5.61057e-16, places=21)
        self.assertEqual(hsp.annotations["identity"], 30)
        self.assertEqual(hsp.annotations["positive"], 36)
        self.assertEqual(hsp.annotations["gaps"], 0)
        hsp = hit[2]
        self.assertEqual(hsp.num, 3)
        self.assertAlmostEqual(hsp.score, 91.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 44.5971)
        self.assertAlmostEqual(hsp.annotations["evalue"], 5.61057e-16, places=21)
        self.assertEqual(hsp.annotations["identity"], 20)
        self.assertEqual(hsp.annotations["positive"], 24)
        self.assertEqual(hsp.annotations["gaps"], 0)
        hsp = hit[3]
        self.assertEqual(hsp.num, 4)
        self.assertAlmostEqual(hsp.score, 112.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 54.2195)
        self.assertAlmostEqual(hsp.annotations["evalue"], 8.37784e-06, places=11)
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
Query : Query_2 Length: 115 Strand: Plus
        gi|296147483:1-350 Saccharomyces cerevisiae S288c Mon2p (MON2) mRNA,
        complete cds
Target: gi|366988334|ref|XM_003673886.1| Length: 115 Strand: Plus
        Naumovozyma castellii CBS 4309 hypothetical protein (NCAS0A09950) mRNA,
        complete cds

Score:54 bits(112), Expect:8e-06,
Identities:38/115(33%),  Positives:58/115(50%),  Gaps:0.115(0%)

gi|366988         0 *LYFNCKMCRFYKSIENSR*SAFRYGWCRC*AL*TLHSNRGQLSVP**T*FYQWSDKIRM
                  0 ||.....|......|.|.|..|.|..|.....||.||......|.||..........|||
Query_2           0 *LNLHREMSSLNEGIQNFRQPASRNRWNG*QSL*ALHCQGRHFSIP*LASQHERECEIRM

gi|366988        60 ACQIFNVLDSF*DLDRLVTCMLDRSVPSFRFGRQ*MQF*VQLFLKGLKTCSLNSH 115
                 60 .........||..|....|.|......||||||...||...|.|.|.|..|...| 115
Query_2          60 SF*LLKTMYSFQYLNGFITSMANG*ISSFRFGR*RTQFCFKLPLHGVKPSSVHGH 115

""",
        )
        hit = record[3]
        self.assertEqual(hit.num, 4)
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
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 303.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 141.737)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.78129e-32, places=37)
        self.assertEqual(hsp.annotations["identity"], 55)
        self.assertEqual(hsp.annotations["positive"], 71)
        self.assertEqual(hsp.annotations["gaps"], 0)
        hsp = hit[1]
        self.assertEqual(hsp.num, 2)
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
Query : Query_2 Length: 55 Strand: Plus
        gi|296147483:1-350 Saccharomyces cerevisiae S288c Mon2p (MON2) mRNA,
        complete cds
Target: gi|255710474|ref|XM_002551475.1| Length: 55 Strand: Plus
        Lachancea thermotolerans CBS 6340 KLTH0A01342p (KLTH0A01342g) mRNA,
        complete cds

Score:45 bits(92), Expect:0.005,
Identities:25/55(45%),  Positives:29/55(53%),  Gaps:0.55(0%)

gi|255710         0 TLSCISAARWVESMNASSTSSILSSGMQLTDDIFCRHCTETLVSLAFLEAHERTK 55
                  0 |...||.||.|.||.||..|.....|...|.|..||||....|.||||.....|| 55
Query_2           0 TFN*ISIAR*VASMKASKISDSRLRGIDGTVDSPCRHCIARVVILAFLDWQANTK 55

""",
        )
        hit = record[4]
        self.assertEqual(hit.num, 5)
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
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 302.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 141.279)
        self.assertAlmostEqual(hsp.annotations["evalue"], 5.19486e-32, places=37)
        self.assertEqual(hsp.annotations["identity"], 57)
        self.assertEqual(hsp.annotations["positive"], 72)
        self.assertEqual(hsp.annotations["gaps"], 0)
        hsp = hit[1]
        self.assertEqual(hsp.num, 2)
        self.assertAlmostEqual(hsp.score, 105.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 51.0121)
        self.assertAlmostEqual(hsp.annotations["evalue"], 8.66978e-12, places=17)
        self.assertEqual(hsp.annotations["identity"], 27)
        self.assertEqual(hsp.annotations["positive"], 33)
        self.assertEqual(hsp.annotations["gaps"], 0)
        hsp = hit[2]
        self.assertEqual(hsp.num, 3)
        self.assertAlmostEqual(hsp.score, 85.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 41.8479)
        self.assertAlmostEqual(hsp.annotations["evalue"], 8.66978e-12, places=17)
        self.assertEqual(hsp.annotations["identity"], 20)
        self.assertEqual(hsp.annotations["positive"], 25)
        self.assertEqual(hsp.annotations["gaps"], 0)
        hsp = hit[3]
        self.assertEqual(hsp.num, 4)
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
Query : Query_2 Length: 99 Strand: Plus
        gi|296147483:1-350 Saccharomyces cerevisiae S288c Mon2p (MON2) mRNA,
        complete cds
Target: gi|254579534|ref|XM_002495708.1| Length: 99 Strand: Plus
        Zygosaccharomyces rouxii hypothetical protein (ZYRO0C02266g) mRNA,
        complete cds

Score:45 bits(92), Expect:0.005,
Identities:31/99(31%),  Positives:53/99(54%),  Gaps:0.99(0%)

gi|254579         0 RVAFTFIRVEEKEYCY*EC*R*IN*NFESSSQL*GIIKTSRFYSTASDVMCIQECQIDYY
                  0 |.||...|.|.|......||.||.|..|.|...||..|||.|.........|.||..|..
Query_2           0 RIAFFIFRIEKKKFNHSPC***IH*DIEKST*F*GARKTSGFRTPFRVGLPIKEC*NDDP

gi|254579        60 INAMFPKIGHSAMYTGR*TLRRT*CVYRGQPAGNGYKTK 99
                 60 .|||.........|..........|...|.....|...| 99
Query_2          60 GNAMPTGTVNRSIYSSKPAV*NFGCLH*GYSSRDGDSIK 99

""",
        )
        hit = record[0]
        self.assertEqual(hit.num, 1)
        hsps = hit[1:5:2]
        self.assertEqual(
            str(hsps),
            """\
Query: Query_2
       gi|296147483:1-350 Saccharomyces cerevisiae S288c Mon2p (MON2) mRNA,
       complete cds
  Hit: gi|296147483|ref|NM_001183135.1| (length=4911)
       Saccharomyces cerevisiae S288c Mon2p (MON2) mRNA, complete cds
       >gi|116616412|gb|EF059095.1| Synthetic construct Saccharomyces cerevisiae
       clone FLH203015.01X MON2, complete sequence
 HSPs: ----  --------  ---------  ------  ---------------  ---------------------
          #   E-value  Bit score    Span      Query range              Hit range
       ----  --------  ---------  ------  ---------------  ---------------------
          0   2.2e-73     278.74     116          [0:116]                [0:116]
          1     9e-71     270.04     116          [0:116]                [0:116]""",
        )

    def test_xml_2226_tblastx_004_writer(self):
        """Writing TBLASTX 2.2.26+ (xml_2226_tblastx_004.xml)."""
        filename = "xml_2226_tblastx_004.xml"
        path = os.path.join("Blast", filename)
        with Blast.parse(path) as records:
            stream = io.BytesIO()
            n = Blast.write(records, stream)
            self.assertEqual(n, 2)
            stream.seek(0)
            written_records = Blast.parse(stream)
            self.check_xml_2226_tblastx_004(written_records)

    def test_xml_21500_tblastx_001_parser(self):
        """Parsing TBLASTX 2.15.0+ (xml_21500_tblastx_001.xml)."""
        filename = "xml_21500_tblastx_001.xml"
        path = os.path.join("Blast", filename)
        with open(path, "rb") as stream:
            records = Blast.parse(stream)
            self.check_xml_21500_tblastx_001_records(records)
        with Blast.parse(path) as records:
            self.check_xml_21500_tblastx_001_records(records)
        with open(path, "rb") as stream:
            record = Blast.read(stream)
        self.check_xml_21500_tblastx_001_record(record, xml2=False)
        record = Blast.read(path)
        self.check_xml_21500_tblastx_001_record(record, xml2=False)
        with Blast.parse(path) as records:
            self.assertEqual(
                str(records),
                """\
Program: TBLASTX 2.15.0+
     db: refseq_rna

  Query: Query_949527 (length=804)
         NC_000964.3:2518023-2518826 Bacillus subtilis subsp. subtilis str. 168
         complete genome
   Hits: ----  -----  ----------------------------------------------------------
            #  # HSP  ID + description
         ----  -----  ----------------------------------------------------------
            0      1  gi|1048027041|ref|XM_017618728.1|  PREDICTED: Rhagoleti...
            1      2  gi|2670929493|ref|XM_062937224.1|  Kwoniella shivajii u...
            2      1  gi|2432161024|ref|XM_053088494.1|  Dioszegia hungarica ...
            3      1  gi|799335188|ref|XM_012195534.1|  Cryptococcus neoforma...
            4      1  gi|1799711371|ref|XM_032001855.1|  Kwoniella shandongen...
            5      1  gi|1102541390|ref|XM_019193349.1|  Kwoniella bestiolae ...
            6      1  gi|1799711369|ref|XM_032001854.1|  Kwoniella shandongen...
            7      1  gi|2592096353|ref|XM_060229193.1|  PREDICTED: Ylistrum ...
            8      1  gi|2044197324|ref|XM_041762518.1|  PREDICTED: Vulpes la...
            9      1  gi|1101784196|ref|XM_019135081.1|  Cryptococcus amylole...""",
            )

    def test_xml2_21500_tblastx_001_parser(self):
        """Parsing TBLASTX 2.15.0+ (xml2_21500_tblastx_001.xml)."""
        filename = "xml2_21500_tblastx_001.xml"
        path = os.path.join("Blast", filename)
        with open(path, "rb") as stream:
            records = Blast.parse(stream)
            self.check_xml_21500_tblastx_001_records(records, xml2=True)
        with Blast.parse(path) as records:
            self.check_xml_21500_tblastx_001_records(records, xml2=True)
        with open(path, "rb") as stream:
            record = Blast.read(stream)
        self.check_xml_21500_tblastx_001_record(record, xml2=True)
        record = Blast.read(path)
        self.check_xml_21500_tblastx_001_record(record, xml2=True)

    def check_xml_21500_tblastx_001_records(self, records, xml2=False):
        self.assertEqual(records.program, "tblastx")
        self.assertEqual(records.version, "TBLASTX 2.15.0+")
        self.assertEqual(
            records.reference,
            'Stephen F. Altschul, Thomas L. Madden, Alejandro A. Schäffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.',
        )
        self.assertEqual(records.db, "refseq_rna")
        if not xml2:
            self.assertIsInstance(records.query, SeqRecord)
            self.assertEqual(records.query.id, "Query_949527")
            self.assertEqual(
                records.query.description,
                "NC_000964.3:2518023-2518826 Bacillus subtilis subsp. subtilis str. 168 complete genome",
            )
            self.assertEqual(repr(records.query.seq), "Seq(None, length=804)")
        self.assertEqual(len(records.param), 5)
        self.assertEqual(records.param["matrix"], "BLOSUM62")
        self.assertAlmostEqual(records.param["expect"], 0.05)
        self.assertEqual(records.param["filter"], "L;")
        if xml2:
            self.assertEqual(records.param["query-gencode"], 1)
            self.assertEqual(records.param["db-gencode"], 1)
        else:
            self.assertEqual(records.param["gap-open"], 11)
            self.assertEqual(records.param["gap-extend"], 1)
        record = next(records)
        self.assertRaises(StopIteration, next, records)
        self.check_xml_21500_tblastx_001_record(record, xml2)

    def check_xml_21500_tblastx_001_record(self, record, xml2):
        if not xml2:
            self.assertEqual(record.num, 1)
        self.assertIsInstance(record.query, SeqRecord)
        self.assertEqual(record.query.id, "Query_949527")
        self.assertEqual(
            record.query.description,
            "NC_000964.3:2518023-2518826 Bacillus subtilis subsp. subtilis str. 168 complete genome",
        )
        self.assertEqual(repr(record.query.seq), "Seq(None, length=804)")
        if xml2:
            self.assertEqual(len(record.query.features), 1)
            feature = record.query.features[0]
            self.assertEqual(feature.type, "masking")
            location = feature.location
            self.assertEqual(
                repr(location), "SimpleLocation(ExactPosition(259), ExactPosition(293))"
            )
        self.assertEqual(len(record.stat), 7)
        self.assertEqual(record.stat["db-num"], 60837734)
        self.assertEqual(record.stat["db-len"], 161330429008)
        if xml2:
            self.assertEqual(record.stat["hsp-len"], 66)
            self.assertEqual(record.stat["eff-space"], 10051826883450)
            self.assertEqual(record.stat["kappa"], -1)
            self.assertEqual(record.stat["lambda"], -1)
            self.assertEqual(record.stat["entropy"], -1)
        else:
            self.assertEqual(record.stat["hsp-len"], 0)
            self.assertEqual(record.stat["eff-space"], 0)
            self.assertAlmostEqual(record.stat["kappa"], 0.133956144488482)
            self.assertAlmostEqual(record.stat["lambda"], 0.317605957635731)
            self.assertAlmostEqual(record.stat["entropy"], 0.401214524497119)
        self.assertEqual(len(record), 10)
        hit = record[0]
        self.assertEqual(hit.num, 1)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|1048027041|ref|XM_017618728.1|")
        self.assertEqual(hit.target.name, "XM_017618728")
        self.assertEqual(
            hit.target.description,
            "PREDICTED: Rhagoletis zephyria response regulator PleD-like (LOC108364862), partial mRNA",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=917)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 150.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 71.6314)
        self.assertAlmostEqual(hsp.annotations["evalue"], 5.49617e-09, places=14)
        self.assertEqual(hsp.annotations["identity"], 29)
        self.assertEqual(hsp.annotations["positive"], 52)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 89],
                          [ 0, 89]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 89))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('AYNGQECLSLFKEKDPDVLVLDIIMPHLDGLAVLERLRESDLKKQPNVIMLTAF...QVS')",
        )
        self.assertEqual(hsp.query.id, "Query_949527")
        self.assertEqual(
            hsp.query.description,
            "NC_000964.3:2518023-2518826 Bacillus subtilis subsp. subtilis str. 168 complete genome",
        )
        self.assertEqual(len(hsp.query.features), 1)
        feature = hsp.query.features[0]
        self.assertEqual(feature.type, "CDS")
        location = feature.location
        self.assertEqual(
            repr(location), "SimpleLocation(ExactPosition(0), ExactPosition(89))"
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('ACGGQEALNLVVEQQPDIVLLDVMMPGIDGFMVCKEMKDNPMTTHIPVVMVTAL...SLT')",
        )
        self.assertEqual(hsp.target.id, "gi|1048027041|ref|XM_017618728.1|")
        self.assertEqual(hsp.target.name, "XM_017618728")
        self.assertEqual(
            hsp.target.description,
            "PREDICTED: Rhagoletis zephyria response regulator PleD-like (LOC108364862), partial mRNA",
        )
        if xml2:
            self.assertEqual(len(record.query.features), 1)
            feature = record.query.features[0]
            self.assertEqual(feature.type, "masking")
            location = feature.location
            self.assertEqual(
                repr(location), "SimpleLocation(ExactPosition(259), ExactPosition(293))"
            )
        self.assertEqual(
            hsp.annotations["midline"],
            "A  GQE L+L  E+ PD+++LD++MP +DG  V + ++++ +     V+M+TA    +   K ++ GA  F+ KP D   L   I+ ++",
        )
        self.assertEqual(
            repr(hsp),
            "<Bio.Blast.HSP target.id='gi|1048027041|ref|XM_017618728.1|' query.id='Query_949527'; 2 rows x 89 columns>",
        )
        self.assertEqual(
            str(hsp),
            """\
Query : Query_949527 Length: 89 Strand: Plus
        NC_000964.3:2518023-2518826 Bacillus subtilis subsp. subtilis str. 168
        complete genome
Target: gi|1048027041|ref|XM_017618728.1| Length: 89 Strand: Plus
        PREDICTED: Rhagoletis zephyria response regulator PleD-like
        (LOC108364862), partial mRNA

Score:71 bits(150), Expect:5e-09,
Identities:29/89(33%),  Positives:52/89(58%),  Gaps:0.89(0%)

gi|104802         0 ACGGQEALNLVVEQQPDIVLLDVMMPGIDGFMVCKEMKDNPMTTHIPVVMVTALHDTEDR
                  0 |..|||.|.|..|..||...||..||..||..|..............|.|.||.......
Query_949         0 AYNGQECLSLFKEKDPDVLVLDIIMPHLDGLAVLERLRESDLKKQPNVIMLTAFGQEDVT

gi|104802        60 VKGINAGADDFLTKPIDETALSARIKSLT 89
                 60 .|....||..|..||.|...|...|.... 89
Query_949        60 KKAVDLGASYFILKPFDMENLVGHIRQVS 89

""",
        )
        hit = record[1]
        self.assertEqual(hit.num, 2)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|2670929493|ref|XM_062937224.1|")
        self.assertEqual(hit.target.name, "XM_062937224")
        self.assertEqual(
            hit.target.description,
            "Kwoniella shivajii uncharacterized protein (IL334_005512), partial mRNA",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=9882)")
        self.assertEqual(len(hit), 2)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 146.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 69.7986)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.95794e-08, places=13)
        self.assertEqual(hsp.annotations["identity"], 32)
        self.assertEqual(hsp.annotations["positive"], 51)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 94],
                          [ 0, 94]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 94))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('IEGQEDMEVIGVAYNGQECLSLFKEKDPDVLVLDIIMPHLDGLAVLERLRESDL...NLV')",
        )
        self.assertEqual(hsp.query.id, "Query_949527")
        self.assertEqual(
            hsp.query.description,
            "NC_000964.3:2518023-2518826 Bacillus subtilis subsp. subtilis str. 168 complete genome",
        )
        self.assertEqual(len(hsp.query.features), 1)
        feature = hsp.query.features[0]
        self.assertEqual(feature.type, "CDS")
        location = feature.location
        self.assertEqual(
            repr(location), "SimpleLocation(ExactPosition(0), ExactPosition(94))"
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('VDDSYDMRKVIEAHDGQEALELCSRALPDLIISDIMMPRLDGFGLLQALKSSSN...ELL')",
        )
        self.assertEqual(hsp.target.id, "gi|2670929493|ref|XM_062937224.1|")
        self.assertEqual(hsp.target.name, "XM_062937224")
        self.assertEqual(
            hsp.target.description,
            "Kwoniella shivajii uncharacterized protein (IL334_005512), partial mRNA",
        )
        self.assertEqual(len(hsp.target.features), 1)
        feature = hsp.target.features[0]
        self.assertEqual(feature.type, "CDS")
        location = feature.location
        self.assertEqual(
            repr(location), "SimpleLocation(ExactPosition(0), ExactPosition(9882))"
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "++   DM  +  A++GQE L L     PD+++ DI+MP LDG  +L+ L+ S       +I+LTA G +D     +  GA  ++ KPF    L+",
        )
        self.assertEqual(
            repr(hsp),
            "<Bio.Blast.HSP target.id='gi|2670929493|ref|XM_062937224.1|' query.id='Query_949527'; 2 rows x 94 columns>",
        )
        self.assertEqual(
            str(hsp),
            """\
Query : Query_949527 Length: 94 Strand: Plus
        NC_000964.3:2518023-2518826 Bacillus subtilis subsp. subtilis str. 168
        complete genome
Target: gi|2670929493|ref|XM_062937224.1| Length: 94 Strand: Plus
        Kwoniella shivajii uncharacterized protein (IL334_005512), partial mRNA

Score:69 bits(146), Expect:2e-08,
Identities:32/94(34%),  Positives:51/94(54%),  Gaps:0.94(0%)

gi|267092         0 VDDSYDMRKVIEAHDGQEALELCSRALPDLIISDIMMPRLDGFGLLQALKSSSNLISVPI
                  0 .....||.....|..|||.|.|.....||....||.||.|||...|..|..|........
Query_949         0 IEGQEDMEVIGVAYNGQECLSLFKEKDPDVLVLDIIMPHLDGLAVLERLRESDLKKQPNV

gi|267092        60 ILLTARGDQDFKVGGLMSGAEDYLSKPFSTPELL 94
                 60 |.|||.|..|........||.....|||....|. 94
Query_949        60 IMLTAFGQEDVTKKAVDLGASYFILKPFDMENLV 94

""",
        )
        hsp = hit[1]
        self.assertEqual(hsp.num, 2)
        self.assertAlmostEqual(hsp.score, 118.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 56.9688)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.000142549)
        self.assertEqual(hsp.annotations["identity"], 25)
        self.assertEqual(hsp.annotations["positive"], 42)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 68],
                          [ 0, 68]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 68))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('AYNGQECLSLFKEKDPDVLVLDIIMPHLDGLAVLERLRESDLKKQPNVIMLTAF...LGA')",
        )
        self.assertEqual(hsp.query.id, "Query_949527")
        self.assertEqual(
            hsp.query.description,
            "NC_000964.3:2518023-2518826 Bacillus subtilis subsp. subtilis str. 168 complete genome",
        )
        self.assertEqual(len(hsp.query.features), 1)
        feature = hsp.query.features[0]
        self.assertEqual(feature.type, "CDS")
        location = feature.location
        self.assertEqual(
            repr(location), "SimpleLocation(ExactPosition(0), ExactPosition(68))"
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('ARDGQEALALCAQQAPDLIISDVMMPNLDGFGLLRALKQSKKLAIIPIIMLTAR...AGA')",
        )
        self.assertEqual(hsp.target.id, "gi|2670929493|ref|XM_062937224.1|")
        self.assertEqual(hsp.target.name, "XM_062937224")
        self.assertEqual(
            hsp.target.description,
            "Kwoniella shivajii uncharacterized protein (IL334_005512), partial mRNA",
        )
        self.assertEqual(len(hsp.target.features), 1)
        feature = hsp.target.features[0]
        self.assertEqual(feature.type, "CDS")
        location = feature.location
        self.assertEqual(
            repr(location), "SimpleLocation(ExactPosition(0), ExactPosition(9882))"
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "A +GQE L+L  ++ PD+++ D++MP+LDG  +L  L++S       +IMLTA G ++     +  GA",
        )
        self.assertEqual(
            repr(hsp),
            "<Bio.Blast.HSP target.id='gi|2670929493|ref|XM_062937224.1|' query.id='Query_949527'; 2 rows x 68 columns>",
        )
        self.assertEqual(
            str(hsp),
            """\
Query : Query_949527 Length: 68 Strand: Plus
        NC_000964.3:2518023-2518826 Bacillus subtilis subsp. subtilis str. 168
        complete genome
Target: gi|2670929493|ref|XM_062937224.1| Length: 68 Strand: Plus
        Kwoniella shivajii uncharacterized protein (IL334_005512), partial mRNA

Score:56 bits(118), Expect:0.0001,
Identities:25/68(37%),  Positives:42/68(62%),  Gaps:0.68(0%)

gi|267092         0 ARDGQEALALCAQQAPDLIISDVMMPNLDGFGLLRALKQSKKLAIIPIIMLTARGGDEAR
                  0 |..|||.|.|.....||....|..||.|||...|..|..|........|||||.|.....
Query_949         0 AYNGQECLSLFKEKDPDVLVLDIIMPHLDGLAVLERLRESDLKKQPNVIMLTAFGQEDVT

gi|267092        60 VDGILAGA 68
                 60 ......|| 68
Query_949        60 KKAVDLGA 68

""",
        )
        hit = record[2]
        self.assertEqual(hit.num, 3)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|2432161024|ref|XM_053088494.1|")
        self.assertEqual(hit.target.name, "XM_053088494")
        self.assertEqual(
            hit.target.description,
            "Dioszegia hungarica CnHHK4 (MKK02DRAFT_32386), partial mRNA",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=5601)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 146.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 69.7986)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.95794e-08, places=13)
        self.assertEqual(hsp.annotations["identity"], 32)
        self.assertEqual(hsp.annotations["positive"], 47)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 82],
                          [ 0, 82]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 82))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('AYNGQECLSLFKEKDPDVLVLDIIMPHLDGLAVLERLRESDLKKQPNVIMLTAF...NLV')",
        )
        self.assertEqual(hsp.query.id, "Query_949527")
        self.assertEqual(
            hsp.query.description,
            "NC_000964.3:2518023-2518826 Bacillus subtilis subsp. subtilis str. 168 complete genome",
        )
        self.assertEqual(len(hsp.query.features), 1)
        feature = hsp.query.features[0]
        self.assertEqual(feature.type, "CDS")
        location = feature.location
        self.assertEqual(
            repr(location), "SimpleLocation(ExactPosition(0), ExactPosition(82))"
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('ARDGQEALEMCKKSLPHVVITDVMMPNLDGFGLLAALKEDPKLSMVPVIMLTAR...ELI')",
        )
        self.assertEqual(hsp.target.id, "gi|2432161024|ref|XM_053088494.1|")
        self.assertEqual(hsp.target.name, "XM_053088494")
        self.assertEqual(
            hsp.target.description,
            "Dioszegia hungarica CnHHK4 (MKK02DRAFT_32386), partial mRNA",
        )
        self.assertEqual(len(hsp.target.features), 1)
        feature = hsp.target.features[0]
        self.assertEqual(feature.type, "CDS")
        location = feature.location
        self.assertEqual(
            repr(location), "SimpleLocation(ExactPosition(0), ExactPosition(5601))"
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "A +GQE L + K+  P V++ D++MP+LDG  +L  L+E        VIMLTA G E+     +  GA  +I KPF+   L+",
        )
        self.assertEqual(
            repr(hsp),
            "<Bio.Blast.HSP target.id='gi|2432161024|ref|XM_053088494.1|' query.id='Query_949527'; 2 rows x 82 columns>",
        )
        self.assertEqual(
            str(hsp),
            """\
Query : Query_949527 Length: 82 Strand: Plus
        NC_000964.3:2518023-2518826 Bacillus subtilis subsp. subtilis str. 168
        complete genome
Target: gi|2432161024|ref|XM_053088494.1| Length: 82 Strand: Plus
        Dioszegia hungarica CnHHK4 (MKK02DRAFT_32386), partial mRNA

Score:69 bits(146), Expect:2e-08,
Identities:32/82(39%),  Positives:47/82(57%),  Gaps:0.82(0%)

gi|243216         0 ARDGQEALEMCKKSLPHVVITDVMMPNLDGFGLLAALKEDPKLSMVPVIMLTARGGEEAK
                  0 |..|||.|...|...|.|...|..||.|||...|..|.|........||||||.|.|...
Query_949         0 AYNGQECLSLFKEKDPDVLVLDIIMPHLDGLAVLERLRESDLKKQPNVIMLTAFGQEDVT

gi|243216        60 VDGLLAGADDYIAKPFNARELI 82
                 60 ......||...|.|||....|. 82
Query_949        60 KKAVDLGASYFILKPFDMENLV 82

""",
        )
        hit = record[3]
        self.assertEqual(hit.num, 4)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|799335188|ref|XM_012195534.1|")
        self.assertEqual(hit.target.name, "XM_012195534")
        self.assertEqual(
            hit.target.description,
            "Cryptococcus neoformans var. grubii H99 hypothetical protein (CNAG_03355), mRNA",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=5764)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 146.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 69.7986)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.95794e-08, places=13)
        self.assertEqual(hsp.annotations["identity"], 31)
        self.assertEqual(hsp.annotations["positive"], 48)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 82],
                          [ 0, 82]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 82))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('AYNGQECLSLFKEKDPDVLVLDIIMPHLDGLAVLERLRESDLKKQPNVIMLTAF...NLV')",
        )
        self.assertEqual(hsp.query.id, "Query_949527")
        self.assertEqual(
            hsp.query.description,
            "NC_000964.3:2518023-2518826 Bacillus subtilis subsp. subtilis str. 168 complete genome",
        )
        self.assertEqual(len(hsp.query.features), 1)
        feature = hsp.query.features[0]
        self.assertEqual(feature.type, "CDS")
        location = feature.location
        self.assertEqual(
            repr(location), "SimpleLocation(ExactPosition(0), ExactPosition(82))"
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('ARDGQEALELCQKSVPDLIISDVMMPHLNGFELLVALKRSKDLKMVPVIMLTAR...EIV')",
        )
        self.assertEqual(hsp.target.id, "gi|799335188|ref|XM_012195534.1|")
        self.assertEqual(hsp.target.name, "XM_012195534")
        self.assertEqual(
            hsp.target.description,
            "Cryptococcus neoformans var. grubii H99 hypothetical protein (CNAG_03355), mRNA",
        )
        self.assertEqual(len(hsp.target.features), 1)
        feature = hsp.target.features[0]
        self.assertEqual(feature.type, "CDS")
        location = feature.location
        self.assertEqual(
            repr(location), "SimpleLocation(ExactPosition(0), ExactPosition(5764))"
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "A +GQE L L ++  PD+++ D++MPHL+G  +L  L+ S   K   VIMLTA G ++     +  GA  ++ KPF    +V",
        )
        self.assertEqual(
            repr(hsp),
            "<Bio.Blast.HSP target.id='gi|799335188|ref|XM_012195534.1|' query.id='Query_949527'; 2 rows x 82 columns>",
        )
        self.assertEqual(
            str(hsp),
            """\
Query : Query_949527 Length: 82 Strand: Plus
        NC_000964.3:2518023-2518826 Bacillus subtilis subsp. subtilis str. 168
        complete genome
Target: gi|799335188|ref|XM_012195534.1| Length: 82 Strand: Plus
        Cryptococcus neoformans var. grubii H99 hypothetical protein
        (CNAG_03355), mRNA

Score:69 bits(146), Expect:2e-08,
Identities:31/82(38%),  Positives:48/82(59%),  Gaps:0.82(0%)

gi|799335         0 ARDGQEALELCQKSVPDLIISDVMMPHLNGFELLVALKRSKDLKMVPVIMLTARGADESK
                  0 |..|||.|.|.....||....|..||||.|...|..|..|...|...||||||.|.....
Query_949         0 AYNGQECLSLFKEKDPDVLVLDIIMPHLDGLAVLERLRESDLKKQPNVIMLTAFGQEDVT

gi|799335        60 VDGIMAGAEDYLAKPFSAREIV 82
                 60 ......||.....|||.....| 82
Query_949        60 KKAVDLGASYFILKPFDMENLV 82

""",
        )
        hit = record[4]
        self.assertEqual(hit.num, 5)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|1799711371|ref|XM_032001855.1|")
        self.assertEqual(hit.target.name, "XM_032001855")
        self.assertEqual(
            hit.target.description,
            "Kwoniella shandongensis uncharacterized protein (CI109_000662), partial mRNA",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=5538)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 145.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 69.3404)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.68988e-08, places=13)
        self.assertEqual(hsp.annotations["identity"], 31)
        self.assertEqual(hsp.annotations["positive"], 51)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 82],
                          [ 0, 82]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 82))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('AYNGQECLSLFKEKDPDVLVLDIIMPHLDGLAVLERLRESDLKKQPNVIMLTAF...NLV')",
        )
        self.assertEqual(hsp.query.id, "Query_949527")
        self.assertEqual(
            hsp.query.description,
            "NC_000964.3:2518023-2518826 Bacillus subtilis subsp. subtilis str. 168 complete genome",
        )
        self.assertEqual(len(hsp.query.features), 1)
        feature = hsp.query.features[0]
        self.assertEqual(feature.type, "CDS")
        location = feature.location
        self.assertEqual(
            repr(location), "SimpleLocation(ExactPosition(0), ExactPosition(82))"
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('ARDGVEALQLCKKQLPNLIITDVMMPNLDGFGLLAALKESRAMKVIPVIMLTAR...EIV')",
        )
        self.assertEqual(hsp.target.id, "gi|1799711371|ref|XM_032001855.1|")
        self.assertEqual(hsp.target.name, "XM_032001855")
        self.assertEqual(
            hsp.target.description,
            "Kwoniella shandongensis uncharacterized protein (CI109_000662), partial mRNA",
        )
        self.assertEqual(len(hsp.target.features), 1)
        feature = hsp.target.features[0]
        self.assertEqual(feature.type, "CDS")
        location = feature.location
        self.assertEqual(
            repr(location), "SimpleLocation(ExactPosition(0), ExactPosition(5538))"
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "A +G E L L K++ P++++ D++MP+LDG  +L  L+ES   K   VIMLTA G ++   + +  GA  ++ KPF+   +V",
        )
        self.assertEqual(
            repr(hsp),
            "<Bio.Blast.HSP target.id='gi|1799711371|ref|XM_032001855.1|' query.id='Query_949527'; 2 rows x 82 columns>",
        )
        self.assertEqual(
            str(hsp),
            """\
Query : Query_949527 Length: 82 Strand: Plus
        NC_000964.3:2518023-2518826 Bacillus subtilis subsp. subtilis str. 168
        complete genome
Target: gi|1799711371|ref|XM_032001855.1| Length: 82 Strand: Plus
        Kwoniella shandongensis uncharacterized protein (CI109_000662), partial
        mRNA

Score:69 bits(145), Expect:3e-08,
Identities:31/82(38%),  Positives:51/82(62%),  Gaps:0.82(0%)

gi|179971         0 ARDGVEALQLCKKQLPNLIITDVMMPNLDGFGLLAALKESRAMKVIPVIMLTARGGDESK
                  0 |..|.|.|.|.|...|.....|..||.|||...|..|.||...|...||||||.|.....
Query_949         0 AYNGQECLSLFKEKDPDVLVLDIIMPHLDGLAVLERLRESDLKKQPNVIMLTAFGQEDVT

gi|179971        60 VEGILAGADDYLAKPFNAREIV 82
                 60 ......||.....|||.....| 82
Query_949        60 KKAVDLGASYFILKPFDMENLV 82

""",
        )
        hit = record[5]
        self.assertEqual(hit.num, 6)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|1102541390|ref|XM_019193349.1|")
        self.assertEqual(hit.target.name, "XM_019193349")
        self.assertEqual(
            hit.target.description,
            "Kwoniella bestiolae CBS 10118 hypothetical protein partial mRNA",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=5319)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 144.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 68.8822)
        self.assertAlmostEqual(hsp.annotations["evalue"], 3.69545e-08, places=13)
        self.assertEqual(hsp.annotations["identity"], 32)
        self.assertEqual(hsp.annotations["positive"], 49)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 82],
                          [ 0, 82]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 82))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('AYNGQECLSLFKEKDPDVLVLDIIMPHLDGLAVLERLRESDLKKQPNVIMLTAF...NLV')",
        )
        self.assertEqual(hsp.query.id, "Query_949527")
        self.assertEqual(
            hsp.query.description,
            "NC_000964.3:2518023-2518826 Bacillus subtilis subsp. subtilis str. 168 complete genome",
        )
        self.assertEqual(len(hsp.query.features), 1)
        feature = hsp.query.features[0]
        self.assertEqual(feature.type, "CDS")
        location = feature.location
        self.assertEqual(
            repr(location), "SimpleLocation(ExactPosition(0), ExactPosition(82))"
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('ARDGQEALEMCGKKMPDLIISDVMMPNLDGFGLLEALKASKELSIIPVIMLTAR...ELV')",
        )
        self.assertEqual(hsp.target.id, "gi|1102541390|ref|XM_019193349.1|")
        self.assertEqual(hsp.target.name, "XM_019193349")
        self.assertEqual(
            hsp.target.description,
            "Kwoniella bestiolae CBS 10118 hypothetical protein partial mRNA",
        )
        self.assertEqual(len(hsp.target.features), 1)
        feature = hsp.target.features[0]
        self.assertEqual(feature.type, "CDS")
        location = feature.location
        self.assertEqual(
            repr(location), "SimpleLocation(ExactPosition(0), ExactPosition(5319))"
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "A +GQE L +  +K PD+++ D++MP+LDG  +LE L+ S       VIMLTA G ++     +  GA  ++ KPF+   LV",
        )
        self.assertEqual(
            repr(hsp),
            "<Bio.Blast.HSP target.id='gi|1102541390|ref|XM_019193349.1|' query.id='Query_949527'; 2 rows x 82 columns>",
        )
        self.assertEqual(
            str(hsp),
            """\
Query : Query_949527 Length: 82 Strand: Plus
        NC_000964.3:2518023-2518826 Bacillus subtilis subsp. subtilis str. 168
        complete genome
Target: gi|1102541390|ref|XM_019193349.1| Length: 82 Strand: Plus
        Kwoniella bestiolae CBS 10118 hypothetical protein partial mRNA

Score:68 bits(144), Expect:4e-08,
Identities:32/82(39%),  Positives:49/82(60%),  Gaps:0.82(0%)

gi|110254         0 ARDGQEALEMCGKKMPDLIISDVMMPNLDGFGLLEALKASKELSIIPVIMLTARGGDEAK
                  0 |..|||.|.....|.||....|..||.|||...||.|..|.......||||||.|.....
Query_949         0 AYNGQECLSLFKEKDPDVLVLDIIMPHLDGLAVLERLRESDLKKQPNVIMLTAFGQEDVT

gi|110254        60 VDGLLAGADDYLAKPFNSRELV 82
                 60 ......||.....|||....|| 82
Query_949        60 KKAVDLGASYFILKPFDMENLV 82

""",
        )
        hit = record[6]
        self.assertEqual(hit.num, 7)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|1799711369|ref|XM_032001854.1|")
        self.assertEqual(hit.target.name, "XM_032001854")
        self.assertEqual(
            hit.target.description,
            "Kwoniella shandongensis uncharacterized protein (CI109_000661), partial mRNA",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=3591)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 143.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 68.424)
        self.assertAlmostEqual(hsp.annotations["evalue"], 5.07694e-08, places=13)
        self.assertEqual(hsp.annotations["identity"], 31)
        self.assertEqual(hsp.annotations["positive"], 53)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 86],
                          [ 0, 86]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 86))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('VIGVAYNGQECLSLFKEKDPDVLVLDIIMPHLDGLAVLERLRESDLKKQPNVIM...NLV')",
        )
        self.assertEqual(hsp.query.id, "Query_949527")
        self.assertEqual(
            hsp.query.description,
            "NC_000964.3:2518023-2518826 Bacillus subtilis subsp. subtilis str. 168 complete genome",
        )
        self.assertEqual(len(hsp.query.features), 1)
        feature = hsp.query.features[0]
        self.assertEqual(feature.type, "CDS")
        location = feature.location
        self.assertEqual(
            repr(location), "SimpleLocation(ExactPosition(0), ExactPosition(86))"
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('VVVEARDGQEALDKCMQIRPDLIITDVMMPVLDGFGLLRALKQSDELKAIPVIM...ELI')",
        )
        self.assertEqual(hsp.target.id, "gi|1799711369|ref|XM_032001854.1|")
        self.assertEqual(hsp.target.name, "XM_032001854")
        self.assertEqual(
            hsp.target.description,
            "Kwoniella shandongensis uncharacterized protein (CI109_000661), partial mRNA",
        )
        self.assertEqual(len(hsp.target.features), 1)
        feature = hsp.target.features[0]
        self.assertEqual(feature.type, "CDS")
        location = feature.location
        self.assertEqual(
            repr(location), "SimpleLocation(ExactPosition(0), ExactPosition(3591))"
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "V+  A +GQE L    +  PD+++ D++MP LDG  +L  L++SD  K   VIM+TA   ++   +A+  GA  +++KPF++  L+",
        )
        self.assertEqual(
            repr(hsp),
            "<Bio.Blast.HSP target.id='gi|1799711369|ref|XM_032001854.1|' query.id='Query_949527'; 2 rows x 86 columns>",
        )
        self.assertEqual(
            str(hsp),
            """\
Query : Query_949527 Length: 86 Strand: Plus
        NC_000964.3:2518023-2518826 Bacillus subtilis subsp. subtilis str. 168
        complete genome
Target: gi|1799711369|ref|XM_032001854.1| Length: 86 Strand: Plus
        Kwoniella shandongensis uncharacterized protein (CI109_000661), partial
        mRNA

Score:68 bits(143), Expect:5e-08,
Identities:31/86(36%),  Positives:53/86(62%),  Gaps:0.86(0%)

gi|179971         0 VVVEARDGQEALDKCMQIRPDLIITDVMMPVLDGFGLLRALKQSDELKAIPVIMVTAHDG
                  0 |...|..|||.|.......||....|..||.|||...|..|..||..|...|||.||...
Query_949         0 VIGVAYNGQECLSLFKEKDPDVLVLDIIMPHLDGLAVLERLRESDLKKQPNVIMLTAFGQ

gi|179971        60 DEAKVEALLGGADDYMVKPFNVRELI 86
                 60 ......|...||.....|||....|. 86
Query_949        60 EDVTKKAVDLGASYFILKPFDMENLV 86

""",
        )
        hit = record[7]
        self.assertEqual(hit.num, 8)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|2592096353|ref|XM_060229193.1|")
        self.assertEqual(hit.target.name, "XM_060229193")
        self.assertEqual(
            hit.target.description,
            "PREDICTED: Ylistrum balloti signal transduction histidine-protein kinase BarA-like (LOC132564539), mRNA",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=1374)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 142.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 67.9658)
        self.assertAlmostEqual(hsp.annotations["evalue"], 6.97488e-08, places=13)
        self.assertEqual(hsp.annotations["identity"], 31)
        self.assertEqual(hsp.annotations["positive"], 51)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 84],
                          [ 0, 84]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 84))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('IGVAYNGQECLSLFKEKDPDVLVLDIIMPHLDGLAVLERLRESDLKKQPNVIML...ENL')",
        )
        self.assertEqual(hsp.query.id, "Query_949527")
        self.assertEqual(
            hsp.query.description,
            "NC_000964.3:2518023-2518826 Bacillus subtilis subsp. subtilis str. 168 complete genome",
        )
        self.assertEqual(len(hsp.query.features), 1)
        feature = hsp.query.features[0]
        self.assertEqual(feature.type, "CDS")
        location = feature.location
        self.assertEqual(
            repr(location), "SimpleLocation(ExactPosition(0), ExactPosition(84))"
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('IHLASNGQEAIQKTKVLDPDLIFLDLHMPDMSGLEVIEILRSITAYRDTPIIIL...DKL')",
        )
        self.assertEqual(hsp.target.id, "gi|2592096353|ref|XM_060229193.1|")
        self.assertEqual(hsp.target.name, "XM_060229193")
        self.assertEqual(
            hsp.target.description,
            "PREDICTED: Ylistrum balloti signal transduction histidine-protein kinase BarA-like (LOC132564539), mRNA",
        )
        self.assertEqual(len(hsp.target.features), 1)
        feature = hsp.target.features[0]
        self.assertEqual(feature.type, "CDS")
        location = feature.location
        self.assertEqual(
            repr(location), "SimpleLocation(ExactPosition(0), ExactPosition(1374))"
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "I +A NGQE +   K  DPD++ LD+ MP + GL V+E LR     +   +I+L+A    +  +KA+ +GAS ++ KP +++ L",
        )
        self.assertEqual(
            repr(hsp),
            "<Bio.Blast.HSP target.id='gi|2592096353|ref|XM_060229193.1|' query.id='Query_949527'; 2 rows x 84 columns>",
        )
        self.assertEqual(
            str(hsp),
            """\
Query : Query_949527 Length: 84 Strand: Plus
        NC_000964.3:2518023-2518826 Bacillus subtilis subsp. subtilis str. 168
        complete genome
Target: gi|2592096353|ref|XM_060229193.1| Length: 84 Strand: Plus
        PREDICTED: Ylistrum balloti signal transduction histidine-protein kinase
        BarA-like (LOC132564539), mRNA

Score:67 bits(142), Expect:7e-08,
Identities:31/84(37%),  Positives:51/84(61%),  Gaps:0.84(0%)

gi|259209         0 IHLASNGQEAIQKTKVLDPDLIFLDLHMPDMSGLEVIEILRSITAYRDTPIIILSADAII
                  0 |..|.||||.....|..|||...||..||...||.|.|.||..........|.|.|....
Query_949         0 IGVAYNGQECLSLFKEKDPDVLVLDIIMPHLDGLAVLERLRESDLKKQPNVIMLTAFGQE

gi|259209        60 EQQEKALAVGASAYLTKPIEIDKL 84
                 60 ....||...|||....||.....| 84
Query_949        60 DVTKKAVDLGASYFILKPFDMENL 84

""",
        )
        hit = record[8]
        self.assertEqual(hit.num, 9)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|2044197324|ref|XM_041762518.1|")
        self.assertEqual(hit.target.name, "XM_041762518")
        self.assertEqual(
            hit.target.description,
            "PREDICTED: Vulpes lagopus sensory transduction protein RegX3-like (LOC121495284), mRNA",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=681)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 142.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 67.9658)
        self.assertAlmostEqual(hsp.annotations["evalue"], 6.97488e-08, places=13)
        self.assertEqual(hsp.annotations["identity"], 30)
        self.assertEqual(hsp.annotations["positive"], 45)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 73],
                          [ 0, 73]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 73))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('PDVLVLDIIMPHLDGLAVLERLRESDLKKQPNVIMLTAFGQEDVTKKAVDLGAS...RQV')",
        )
        self.assertEqual(hsp.query.id, "Query_949527")
        self.assertEqual(
            hsp.query.description,
            "NC_000964.3:2518023-2518826 Bacillus subtilis subsp. subtilis str. 168 complete genome",
        )
        self.assertEqual(len(hsp.query.features), 1)
        feature = hsp.query.features[0]
        self.assertEqual(feature.type, "CDS")
        location = feature.location
        self.assertEqual(
            repr(location), "SimpleLocation(ExactPosition(0), ExactPosition(73))"
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('PDLILLDIMLPGSDGLSSLEELRAKKSSESIPVIMATAKGTEFDKVKGLDMGAD...KAV')",
        )
        self.assertEqual(hsp.target.id, "gi|2044197324|ref|XM_041762518.1|")
        self.assertEqual(hsp.target.name, "XM_041762518")
        self.assertEqual(
            hsp.target.description,
            "PREDICTED: Vulpes lagopus sensory transduction protein RegX3-like (LOC121495284), mRNA",
        )
        self.assertEqual(len(hsp.target.features), 1)
        feature = hsp.target.features[0]
        self.assertEqual(feature.type, "CDS")
        location = feature.location
        self.assertEqual(
            repr(location), "SimpleLocation(ExactPosition(0), ExactPosition(681))"
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "PD+++LDI++P  DGL+ LE LR     +   VIM TA G E    K +D+GA  +++KPF M  ++  I+ V",
        )
        self.assertEqual(
            repr(hsp),
            "<Bio.Blast.HSP target.id='gi|2044197324|ref|XM_041762518.1|' query.id='Query_949527'; 2 rows x 73 columns>",
        )
        self.assertEqual(
            str(hsp),
            """\
Query : Query_949527 Length: 73 Strand: Plus
        NC_000964.3:2518023-2518826 Bacillus subtilis subsp. subtilis str. 168
        complete genome
Target: gi|2044197324|ref|XM_041762518.1| Length: 73 Strand: Plus
        PREDICTED: Vulpes lagopus sensory transduction protein RegX3-like
        (LOC121495284), mRNA

Score:67 bits(142), Expect:7e-08,
Identities:30/73(41%),  Positives:45/73(62%),  Gaps:0.73(0%)

gi|204419         0 PDLILLDIMLPGSDGLSSLEELRAKKSSESIPVIMATAKGTEFDKVKGLDMGADDYLVKP
                  0 ||...|||..|..|||..||.||.........|||.||.|.|....|..|.||.....||
Query_949         0 PDVLVLDIIMPHLDGLAVLERLRESDLKKQPNVIMLTAFGQEDVTKKAVDLGASYFILKP

gi|204419        60 FGMMEMISRIKAV 73
                 60 |.|......|..| 73
Query_949        60 FDMENLVGHIRQV 73

""",
        )
        hit = record[9]
        self.assertEqual(hit.num, 10)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gi|1101784196|ref|XM_019135081.1|")
        self.assertEqual(hit.target.name, "XM_019135081")
        self.assertEqual(
            hit.target.description,
            "Cryptococcus amylolentus CBS 6039 hypothetical protein (L202_01642), partial mRNA",
        )
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=5487)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 142.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 67.9658)
        self.assertAlmostEqual(hsp.annotations["evalue"], 6.97488e-08, places=13)
        self.assertEqual(hsp.annotations["identity"], 28)
        self.assertEqual(hsp.annotations["positive"], 50)
        self.assertEqual(hsp.annotations["gaps"], 0)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[ 0, 82],
                          [ 0, 82]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 82))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq('AYNGQECLSLFKEKDPDVLVLDIIMPHLDGLAVLERLRESDLKKQPNVIMLTAF...NLV')",
        )
        self.assertEqual(hsp.query.id, "Query_949527")
        self.assertEqual(
            hsp.query.description,
            "NC_000964.3:2518023-2518826 Bacillus subtilis subsp. subtilis str. 168 complete genome",
        )
        self.assertEqual(len(hsp.query.features), 1)
        feature = hsp.query.features[0]
        self.assertEqual(feature.type, "CDS")
        location = feature.location
        self.assertEqual(
            repr(location), "SimpleLocation(ExactPosition(0), ExactPosition(82))"
        )
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('ARDGREALELCAKQKPNLIISDVMMPHVDGFELLTTLKDSSEFRMIPVIMLTAR...EIV')",
        )
        self.assertEqual(hsp.target.id, "gi|1101784196|ref|XM_019135081.1|")
        self.assertEqual(hsp.target.name, "XM_019135081")
        self.assertEqual(
            hsp.target.description,
            "Cryptococcus amylolentus CBS 6039 hypothetical protein (L202_01642), partial mRNA",
        )
        self.assertEqual(len(hsp.target.features), 1)
        feature = hsp.target.features[0]
        self.assertEqual(feature.type, "CDS")
        location = feature.location
        self.assertEqual(
            repr(location), "SimpleLocation(ExactPosition(0), ExactPosition(5487))"
        )
        self.assertEqual(
            hsp.annotations["midline"],
            "A +G+E L L  ++ P++++ D++MPH+DG  +L  L++S   +   VIMLTA G ++     +  GA  ++ KPF+   +V",
        )
        self.assertEqual(
            repr(hsp),
            "<Bio.Blast.HSP target.id='gi|1101784196|ref|XM_019135081.1|' query.id='Query_949527'; 2 rows x 82 columns>",
        )
        self.assertEqual(
            str(hsp),
            """\
Query : Query_949527 Length: 82 Strand: Plus
        NC_000964.3:2518023-2518826 Bacillus subtilis subsp. subtilis str. 168
        complete genome
Target: gi|1101784196|ref|XM_019135081.1| Length: 82 Strand: Plus
        Cryptococcus amylolentus CBS 6039 hypothetical protein (L202_01642),
        partial mRNA

Score:67 bits(142), Expect:7e-08,
Identities:28/82(34%),  Positives:50/82(61%),  Gaps:0.82(0%)

gi|110178         0 ARDGREALELCAKQKPNLIISDVMMPHVDGFELLTTLKDSSEFRMIPVIMLTARGADESK
                  0 |..|.|.|.|.....|.....|..|||.||...|..|..|.......||||||.|.....
Query_949         0 AYNGQECLSLFKEKDPDVLVLDIIMPHLDGLAVLERLRESDLKKQPNVIMLTAFGQEDVT

gi|110178        60 VSGIMAGAEDYLAKPFNAREIV 82
                 60 ......||.....|||.....| 82
Query_949        60 KKAVDLGASYFILKPFDMENLV 82

""",
        )

    def test_xml_21500_tblastx_001_writer(self):
        """Writing TBLASTX 2.15.0+ (xml_21500_tblastx_001.xml)."""
        filename = "xml_21500_tblastx_001.xml"
        path = os.path.join("Blast", filename)
        with Blast.parse(path) as records:
            stream = io.BytesIO()
            n = Blast.write(records, stream)
            self.assertEqual(n, 1)
            stream.seek(0)
            written_records = Blast.parse(stream)
            self.check_xml_21500_tblastx_001_records(written_records)

    def test_xml2_21500_tblastx_001_writer(self):
        """Writing TBLASTX 2.15.0+ XML2 (xml2_21500_tblastx_001.xml)."""
        filename = "xml2_21500_tblastx_001.xml"
        path = os.path.join("Blast", filename)
        with Blast.parse(path) as records:
            stream = io.BytesIO()
            n = Blast.write(records, stream, fmt="XML2")
            self.assertEqual(n, 1)
            stream.seek(0)
            written_records = Blast.parse(stream)
            self.check_xml_21500_tblastx_001_records(written_records, xml2=True)


class TestRPSBlast(unittest.TestCase):
    """Test the Blast XML parser for rpsblast output."""

    def test_xml_21500_rpsblast_001_parser(self):
        """Parsing RPSBLAST 2.15.0+ (xml_21500_rpsblast_001.xml)."""
        filename = "xml_21500_rpsblast_001.xml"
        path = os.path.join("Blast", filename)
        with open(path, "rb") as stream:
            records = Blast.parse(stream)
            self.check_xml_21500_rpsblast_001_records(records)
        with Blast.parse(path) as records:
            self.check_xml_21500_rpsblast_001_records(records)
        with open(path, "rb") as stream:
            record = Blast.read(stream)
        self.check_xml_21500_rpsblast_001_record(record)
        record = Blast.read(path)
        self.check_xml_21500_rpsblast_001_record(record)

    def test_xml2_21500_rpsblast_001_parser(self):
        """Parsing RPSBLAST 2.15.0+ (xml2_21500_rpsblast_001.xml)."""
        filename = "xml2_21500_rpsblast_001.xml"
        path = os.path.join("Blast", filename)
        with open(path, "rb") as stream:
            records = Blast.parse(stream)
            self.check_xml_21500_rpsblast_001_records(records, xml2=True)
        with Blast.parse(path) as records:
            self.check_xml_21500_rpsblast_001_records(records, xml2=True)
        with open(path, "rb") as stream:
            record = Blast.read(stream)
        self.check_xml_21500_rpsblast_001_record(record, xml2=True)
        record = Blast.read(path)
        self.check_xml_21500_rpsblast_001_record(record, xml2=True)

    def test_xml_21500_rpsblast_001_writer(self):
        """Writing rpsblast 2.15.0+ (xml_21500_rpsblast_001.xml)."""
        filename = "xml_21500_rpsblast_001.xml"
        path = os.path.join("Blast", filename)
        with Blast.parse(path) as records:
            stream = io.BytesIO()
            n = Blast.write(records, stream)
            self.assertEqual(n, 1)
            stream.seek(0)
            written_records = Blast.parse(stream)
            self.check_xml_21500_rpsblast_001_records(written_records)

    def test_xml2_21500_rpsblast_001_writer(self):
        """Writing rpsblast 2.9.0+ XML2 (xml2_21500_rpsblast_001_v2.xml)."""
        filename = "xml2_21500_rpsblast_001.xml"
        path = os.path.join("Blast", filename)
        with Blast.parse(path) as records:
            stream = io.BytesIO()
            n = Blast.write(records, stream, fmt="XML2")
            self.assertEqual(n, 1)
            stream.seek(0)
            written_records = Blast.parse(stream)
            self.check_xml_21500_rpsblast_001_records(written_records, xml2=True)

    def check_xml_21500_rpsblast_001_records(self, records, xml2=False):
        self.assertEqual(records.program, "rpsblast")
        self.assertEqual(records.version, "RPSBLAST 2.15.0+")
        self.assertEqual(
            records.reference,
            'Stephen F. Altschul, Thomas L. Madden, Alejandro A. Schäffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.',
        )
        self.assertEqual(records.db, "Cdd")
        if not xml2:
            self.assertIsInstance(records.query, SeqRecord)
            self.assertEqual(records.query.id, "Query_1")
            self.assertEqual(
                records.query.description,
                "pdb|3FAJ|A Chain A, Structure of the structural protein P131 of the archaeal virus Acidianus Two-tailed virus (ATV)",
            )
            self.assertEqual(repr(records.query.seq), "Seq(None, length=151)")
        self.assertEqual(records.param["matrix"], "BLOSUM62")
        self.assertAlmostEqual(records.param["expect"], 1)
        self.assertEqual(records.param["gap-open"], 11)
        self.assertEqual(records.param["gap-extend"], 1)
        self.assertEqual(records.param["filter"], "F")
        if xml2:
            self.assertEqual(records.param["cbs"], 1)
            self.assertEqual(len(records.param), 6)
        else:
            self.assertEqual(len(records.param), 5)
        record = next(records)
        self.assertRaises(StopIteration, next, records)
        self.check_xml_21500_rpsblast_001_record(record, xml2)

    def check_xml_21500_rpsblast_001_record(self, record, xml2=False):
        if not xml2:
            self.assertEqual(record.num, 1)
        self.assertEqual(
            repr(record),
            "<Bio.Blast.Record query.id='Query_1'; 2 hits>",
        )
        self.assertIsInstance(record.query, SeqRecord)
        self.assertEqual(record.query.id, "Query_1")
        self.assertEqual(
            record.query.description,
            "pdb|3FAJ|A Chain A, Structure of the structural protein P131 of the archaeal virus Acidianus Two-tailed virus (ATV)",
        )
        self.assertEqual(repr(record.query.seq), "Seq(None, length=151)")
        self.assertEqual(len(record.stat), 7)
        self.assertEqual(record.stat["db-num"], 59693)
        self.assertEqual(record.stat["db-len"], 13521240)
        self.assertEqual(record.stat["hsp-len"], 89)
        self.assertEqual(record.stat["eff-space"], 508930906)
        self.assertAlmostEqual(record.stat["kappa"], 0.041)
        self.assertAlmostEqual(record.stat["lambda"], 0.267)
        self.assertAlmostEqual(record.stat["entropy"], 0.14)
        self.assertEqual(len(record), 2)
        hit = record[0]
        self.assertEqual(hit.num, 1)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gnl|CDD|165101")
        self.assertEqual(
            hit.target.description,
            "PHA02734, PHA02734, coat protein; Provisional.",
        )
        self.assertEqual(hit.target.name, "165101")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=149)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 520.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 204.685)
        self.assertAlmostEqual(hsp.annotations["evalue"], 9.29691e-69, places=74)
        self.assertEqual(hsp.annotations["identity"], 94)
        self.assertEqual(hsp.annotations["positive"], 103)
        self.assertEqual(hsp.annotations["gaps"], 18)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[  0, 77, 93, 110, 112, 149],
                          [ 20, 97, 97, 114, 114, 151]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 149))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({20: 'MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEA...MAS'}, length=151)",
        )
        self.assertEqual(hsp.query.id, "Query_1")
        self.assertEqual(
            hsp.query.description,
            "pdb|3FAJ|A Chain A, Structure of the structural protein P131 of the archaeal virus Acidianus Two-tailed virus (ATV)",
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq('MAKEEPKKGDYAGGAAKILDGFEAGQLGFPEVSLKLAGEEANARKAGDANAKAA...AAM')",
        )
        self.assertEqual(hsp.target.id, "gnl|CDD|165101")
        self.assertEqual(
            hsp.target.description,
            "PHA02734, PHA02734, coat protein; Provisional.",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(
            hsp.annotations["midline"],
            "MAK EPKKGDYAGGA KILD FE GQLG+PEV+LKLAGEEANAR+AGD   K AIHAI+KMI DAMKP RNKG GF+                SQ IPGE+ AQV +  E  YQQAKAFLA+PA   R   + E LSKGAK LA A A ",
        )
        hit = record[1]
        self.assertEqual(hit.num, 2)
        self.assertIsInstance(hit.target, SeqRecord)
        self.assertEqual(hit.target.id, "gnl|CDD|410801")
        self.assertEqual(
            hit.target.description,
            "cd20027, FH_FOXL1, Forkhead (FH) domain found in Forkhead box protein L1 (FOXL1) and similar proteins.  FOXL1, also called Forkhead-related protein FKHL11 or Forkhead-related transcription factor 7 (FREAC-7), acts as a transcription factor required for proper proliferation and differentiation in the gastrointestinal epithelium. It may play a critical role in suppressing tumorigenesis. The FH domain is a winged helix DNA-binding domain. FOX transcription factors recognize the core sequence 5'-(A/C)AA(C/T)A-3'.",
        )
        self.assertEqual(hit.target.name, "410801")
        self.assertEqual(repr(hit.target.seq), "Seq(None, length=98)")
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 65.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 29.0095)
        self.assertAlmostEqual(hsp.annotations["evalue"], 0.517343)
        self.assertEqual(hsp.annotations["identity"], 13)
        self.assertEqual(hsp.annotations["positive"], 15)
        self.assertEqual(hsp.annotations["gaps"], 4)
        self.assertTrue(
            np.array_equal(
                hsp.coordinates,
                # fmt: off
                np.array([[63, 66, 66, 90],
                          [14, 17, 21, 45]])
                # fmt: on
            )
        )
        self.assertEqual(hsp.shape, (2, 31))
        self.assertEqual(
            repr(hsp.query.seq),
            "Seq({14: 'VPRGSHMAKYEPKKGDYAGGAVKILDMFENG'}, length=151)",
        )
        self.assertEqual(hsp.query.id, "Query_1")
        self.assertEqual(
            hsp.query.description,
            "pdb|3FAJ|A Chain A, Structure of the structural protein P131 of the archaeal virus Acidianus Two-tailed virus (ATV)",
        )
        self.assertEqual(len(hsp.query.features), 0)
        self.assertEqual(
            repr(hsp.target.seq),
            "Seq({63: 'VPREKGRPGKGNYWTLDPDCEEMFENG'}, length=98)",
        )
        self.assertEqual(hsp.target.id, "gnl|CDD|410801")
        self.assertEqual(
            hsp.target.description,
            "cd20027, FH_FOXL1, Forkhead (FH) domain found in Forkhead box protein L1 (FOXL1) and similar proteins.  FOXL1, also called Forkhead-related protein FKHL11 or Forkhead-related transcription factor 7 (FREAC-7), acts as a transcription factor required for proper proliferation and differentiation in the gastrointestinal epithelium. It may play a critical role in suppressing tumorigenesis. The FH domain is a winged helix DNA-binding domain. FOX transcription factors recognize the core sequence 5'-(A/C)AA(C/T)A-3'.",
        )
        self.assertEqual(len(hsp.target.features), 0)
        self.assertEqual(hsp.annotations["midline"], "VPR     K  P KG+Y        +MFENG")


class TestPSIBlast(unittest.TestCase):
    """Test the Blast XML parser for psiblast output."""

    def test_xml_21500_psiblast_001_parser(self):
        """Parsing PSIBLAST 2.15.0+ (xml_21500_psiblast_001.xml)."""
        filename = "xml_21500_psiblast_001.xml"
        path = os.path.join("Blast", filename)
        with open(path, "rb") as stream:
            records = Blast.parse(stream)
            self.check_xml_21500_psiblast_001_records(records)
        with Blast.parse(path) as records:
            self.check_xml_21500_psiblast_001_records(records)
        with open(path, "rb") as stream:
            record = Blast.read(stream)
        self.check_xml_21500_psiblast_001_record(record)
        record = Blast.read(path)
        self.check_xml_21500_psiblast_001_record(record)

    def test_xml2_21500_psiblast_001_parser(self):
        """Parsing PSIBLAST 2.15.0+ (xml2_21500_psiblast_001.xml)."""
        filename = "xml2_21500_psiblast_001.xml"
        path = os.path.join("Blast", filename)
        with open(path, "rb") as stream:
            records = Blast.parse(stream)
            self.check_xml_21500_psiblast_001_records(records, xml2=True)
        with Blast.parse(path) as records:
            self.check_xml_21500_psiblast_001_records(records, xml2=True)
        with open(path, "rb") as stream:
            record = Blast.read(stream)
        self.check_xml_21500_psiblast_001_record(record, xml2=True)
        record = Blast.read(path)
        self.check_xml_21500_psiblast_001_record(record, xml2=True)

    def test_xml_21500_psiblast_001_writer(self):
        """Writing psiblast 2.15.0+ (xml_21500_psiblast_001.xml)."""
        filename = "xml_21500_psiblast_001.xml"
        path = os.path.join("Blast", filename)
        with Blast.parse(path) as records:
            stream = io.BytesIO()
            n = Blast.write(records, stream)
            self.assertEqual(n, 1)
            stream.seek(0)
            written_records = Blast.parse(stream)
            self.check_xml_21500_psiblast_001_records(written_records)

    def test_xml2_21500_psiblast_001_writer(self):
        """Writing psiblast 2.9.0+ XML2 (xml2_21500_psiblast_001_v2.xml)."""
        filename = "xml2_21500_psiblast_001.xml"
        path = os.path.join("Blast", filename)
        with Blast.parse(path) as records:
            stream = io.BytesIO()
            n = Blast.write(records, stream, fmt="XML2")
            self.assertEqual(n, 1)
            stream.seek(0)
            written_records = Blast.parse(stream)
            self.check_xml_21500_psiblast_001_records(written_records, xml2=True)

    def check_xml_21500_psiblast_001_records(self, records, xml2=False):
        self.assertEqual(records.program, "psiblast")
        self.assertEqual(records.version, "PSIBLAST 2.15.0+")
        self.assertEqual(
            records.reference,
            'Alejandro A. Schäffer, L. Aravind, Thomas L. Madden, Sergei Shavirin, John L. Spouge, Yuri I. Wolf, Eugene V. Koonin, and Stephen F. Altschul (2001), "Improving the accuracy of PSI-BLAST protein database searches with composition-based statistics and other refinements", Nucleic Acids Res. 29:2994-3005.',
        )
        self.assertEqual(records.db, "swissprot")
        if not xml2:
            self.assertIsInstance(records.query, SeqRecord)
            self.assertEqual(records.query.id, "Query_1")
            self.assertEqual(
                records.query.description,
                "WP_001234791.1 Sec-independent protein translocase subunit TatA [Shigella flexneri]",
            )
            self.assertEqual(repr(records.query.seq), "Seq(None, length=103)")
        self.assertEqual(records.param["matrix"], "BLOSUM62")
        self.assertAlmostEqual(records.param["expect"], 1e-30, places=36)
        self.assertEqual(records.param["gap-open"], 11)
        self.assertEqual(records.param["gap-extend"], 1)
        self.assertEqual(records.param["filter"], "F")
        if xml2:
            self.assertEqual(records.param["cbs"], 2)
            self.assertEqual(len(records.param), 6)
        else:
            self.assertEqual(len(records.param), 5)
        record = next(records)
        self.assertRaises(StopIteration, next, records)
        self.check_xml_21500_psiblast_001_record(record, xml2)

    def check_xml_21500_psiblast_001_record(self, record, xml2=False):
        self.assertEqual(record.num, 1)
        self.assertEqual(
            repr(record),
            "<Bio.Blast.Record query.id='%s'; 2 hits>" % record.query.id,
        )
        self.assertIsInstance(record.query, SeqRecord)
        if xml2:
            self.assertEqual(record.query.id, "lcl|Query_1")
        else:
            self.assertEqual(record.query.id, "Query_1")
            self.assertEqual(
                record.query.description,
                "WP_001234791.1 Sec-independent protein translocase subunit TatA [Shigella flexneri]",
            )
        self.assertEqual(repr(record.query.seq), "Seq(None, length=103)")
        self.assertEqual(len(record.stat), 7)
        self.assertEqual(record.stat["db-num"], 482816)
        self.assertEqual(record.stat["db-len"], 183558113)
        self.assertEqual(record.stat["hsp-len"], 72)
        self.assertEqual(record.stat["eff-space"], 4627826878)
        self.assertAlmostEqual(record.stat["kappa"], 0.041)
        self.assertAlmostEqual(record.stat["lambda"], 0.267)
        self.assertAlmostEqual(record.stat["entropy"], 0.14)
        self.assertEqual(len(record), 2)
        hit = record[0]
        self.assertEqual(hit.num, 1)
        target = hit.target
        self.assertIsInstance(target, SeqRecord)
        self.assertEqual(target.id, "sp|P69428.1|")
        self.assertEqual(target.name, "P69428")
        seq = target.seq
        self.assertEqual(repr(seq), "Seq(None, length=89)")
        if xml2:
            self.assertEqual(
                target.description,
                "RecName: Full=Sec-independent protein translocase protein TatA [Escherichia coli K-12]",
            )
            self.assertEqual(target.annotations["taxid"], 83333)
            self.assertIs(target, hit.targets[0])
            self.assertEqual(len(hit.targets), 4)
            target = hit.targets[1]
            self.assertIsInstance(target, SeqRecord)
            self.assertEqual(target.id, "sp|P69429.1|")
            self.assertEqual(target.name, "P69429")
            self.assertIs(target.seq, seq)
            self.assertEqual(
                target.description,
                "RecName: Full=Sec-independent protein translocase protein TatA [Escherichia coli CFT073]",
            )
            self.assertEqual(target.annotations["taxid"], 199310)
            target = hit.targets[2]
            self.assertIsInstance(target, SeqRecord)
            self.assertEqual(target.id, "sp|P69430.1|")
            self.assertEqual(target.name, "P69430")
            self.assertIs(target.seq, seq)
            self.assertEqual(
                target.description,
                "RecName: Full=Sec-independent protein translocase protein TatA [Escherichia coli O157:H7]",
            )
            self.assertEqual(target.annotations["taxid"], 83334)
            target = hit.targets[3]
            self.assertIsInstance(target, SeqRecord)
            self.assertEqual(target.id, "sp|P69431.1|")
            self.assertEqual(target.name, "P69431")
            self.assertIs(target.seq, seq)
            self.assertEqual(
                target.description,
                "RecName: Full=Sec-independent protein translocase protein TatA [Shigella flexneri]",
            )
            self.assertEqual(target.annotations["taxid"], 623)
        else:
            self.assertEqual(
                target.description,
                "RecName: Full=Sec-independent protein translocase protein TatA [Escherichia coli K-12] >sp|P69429.1| RecName: Full=Sec-independent protein translocase protein TatA [Escherichia coli CFT073] >sp|P69430.1| RecName: Full=Sec-independent protein translocase protein TatA [Escherichia coli O157:H7] >sp|P69431.1| RecName: Full=Sec-independent protein translocase protein TatA [Shigella flexneri]",
            )
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 448.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 177.178)
        self.assertAlmostEqual(hsp.annotations["evalue"], 2.3039e-58, places=64)
        self.assertEqual(hsp.annotations["identity"], 89)
        if xml2:
            self.assertEqual(hsp.shape, (2, 0))
        else:
            self.assertEqual(hsp.annotations["positive"], 89)
            self.assertEqual(hsp.annotations["gaps"], 0)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[ 0,  89],
                              [14, 103]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 89))
        if xml2:
            self.assertEqual(repr(hsp.query.seq), "Seq(None, length=103)")
            self.assertEqual(hsp.query.id, "lcl|Query_1")
        else:
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq({14: 'MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...EQV'}, length=103)",
            )
            self.assertEqual(hsp.query.id, "Query_1")
            self.assertEqual(
                hsp.query.description,
                "WP_001234791.1 Sec-independent protein translocase subunit TatA [Shigella flexneri]",
            )
        self.assertEqual(len(hsp.query.features), 0)
        if xml2:
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq(None, length=89)",
            )
        else:
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq('MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...EQV')",
            )
        self.assertEqual(hsp.target.id, "sp|P69428.1|")
        if xml2:
            self.assertEqual(
                hsp.target.description,
                "RecName: Full=Sec-independent protein translocase protein TatA [Escherichia coli K-12]",
            )
        else:
            self.assertEqual(
                hsp.target.description,
                "RecName: Full=Sec-independent protein translocase protein TatA [Escherichia coli K-12] >sp|P69429.1| RecName: Full=Sec-independent protein translocase protein TatA [Escherichia coli CFT073] >sp|P69430.1| RecName: Full=Sec-independent protein translocase protein TatA [Escherichia coli O157:H7] >sp|P69431.1| RecName: Full=Sec-independent protein translocase protein TatA [Shigella flexneri]",
            )
        self.assertEqual(len(hsp.target.features), 0)
        if not xml2:
            self.assertEqual(
                hsp.annotations["midline"],
                "MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV",
            )
        hit = record[1]
        self.assertEqual(hit.num, 2)
        target = hit.target
        self.assertIsInstance(target, SeqRecord)
        self.assertEqual(target.id, "sp|P0A2H3.1|")
        self.assertEqual(target.name, "P0A2H3")
        seq = target.seq
        self.assertEqual(repr(seq), "Seq(None, length=84)")
        if xml2:
            self.assertEqual(
                target.description,
                "RecName: Full=Sec-independent protein translocase protein TatA [Salmonella enterica subsp. enterica serovar Typhimurium str. LT2]",
            )
            self.assertEqual(target.annotations["taxid"], 99287)
            self.assertIs(target, hit.targets[0])
            self.assertEqual(len(hit.targets), 2)
            target = hit.targets[1]
            self.assertIsInstance(target, SeqRecord)
            self.assertEqual(target.id, "sp|P0A2H4.1|")
            self.assertEqual(target.name, "P0A2H4")
            self.assertIs(target.seq, seq)
            self.assertEqual(
                target.description,
                "RecName: Full=Sec-independent protein translocase protein TatA [Salmonella enterica subsp. enterica serovar Typhi]",
            )
            self.assertEqual(target.annotations["taxid"], 90370)
        else:
            self.assertEqual(
                hit.target.description,
                "RecName: Full=Sec-independent protein translocase protein TatA [Salmonella enterica subsp. enterica serovar Typhimurium str. LT2] >sp|P0A2H4.1| RecName: Full=Sec-independent protein translocase protein TatA [Salmonella enterica subsp. enterica serovar Typhi]",
            )
        self.assertEqual(len(hit), 1)
        hsp = hit[0]
        self.assertEqual(hsp.num, 1)
        self.assertAlmostEqual(hsp.score, 358.0)
        self.assertAlmostEqual(hsp.annotations["bit score"], 142.51)
        self.assertAlmostEqual(hsp.annotations["evalue"], 1.0691e-44, places=50)
        self.assertEqual(hsp.annotations["identity"], 75)
        if xml2:
            self.assertEqual(hsp.shape, (2, 0))
        else:
            self.assertEqual(hsp.annotations["positive"], 79)
            self.assertEqual(hsp.annotations["gaps"], 5)
            self.assertTrue(
                np.array_equal(
                    hsp.coordinates,
                    # fmt: off
                    np.array([[  0, 69, 69,  84],
                              [ 14, 83, 88, 103]])
                    # fmt: on
                )
            )
            self.assertEqual(hsp.shape, (2, 89))
        if xml2:
            self.assertEqual(repr(hsp.query.seq), "Seq(None, length=103)")
        else:
            self.assertEqual(
                repr(hsp.query.seq),
                "Seq({14: 'MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTS...EQV'}, length=103)",
            )
        if xml2:
            self.assertEqual(hsp.query.id, "lcl|Query_1")
        else:
            self.assertEqual(hsp.query.id, "Query_1")
            self.assertEqual(
                hsp.query.description,
                "WP_001234791.1 Sec-independent protein translocase subunit TatA [Shigella flexneri]",
            )
        self.assertEqual(len(hsp.query.features), 0)
        if xml2:
            self.assertEqual(repr(hsp.target.seq), "Seq(None, length=84)")
        else:
            self.assertEqual(
                repr(hsp.target.seq),
                "Seq('MGGISIWQLLIVAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDDAKQDKTS...EQV')",
            )
        self.assertEqual(hsp.target.id, "sp|P0A2H3.1|")
        if xml2:
            self.assertEqual(
                hsp.target.description,
                "RecName: Full=Sec-independent protein translocase protein TatA [Salmonella enterica subsp. enterica serovar Typhimurium str. LT2]",
            )
        else:
            self.assertEqual(
                hsp.target.description,
                "RecName: Full=Sec-independent protein translocase protein TatA [Salmonella enterica subsp. enterica serovar Typhimurium str. LT2] >sp|P0A2H4.1| RecName: Full=Sec-independent protein translocase protein TatA [Salmonella enterica subsp. enterica serovar Typhi]",
            )
        self.assertEqual(len(hsp.target.features), 0)
        if not xml2:
            self.assertEqual(
                hsp.annotations["midline"],
                "MGGISIWQLLI+AVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDD+ KQDKTSQDADFTAK+IADKQ      +AK EDAK  DKEQV",
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
        message = r"^premature end of XML file: line [0-9]\d*, column [0-9]\d*$"
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
        message = r"^premature end of XML file: line [0-9]\d*, column [0-9]\d*$"
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
        message = r"^premature end of XML file: line [0-9]\d*, column [0-9]\d*$"
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
        message = r"^premature end of XML file: line [0-9]\d*, column [0-9]\d*$"
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
            with Blast.parse(path):
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
