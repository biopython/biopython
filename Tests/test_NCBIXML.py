# Copyright 2005 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Test the NCBI XML parser."""

import os
import unittest
from Bio.Blast import NCBIXML

E_VALUE_THRESH = 1e-10


class TestNCBIXML(unittest.TestCase):
    """Tests for the NCBI XML parser."""

    def test_xml_2212L_blastp_001(self):
        """Parsing BLASTP 2.2.12, gi|49176427|ref|NP_418280.3| (xml_2212L_blastp_001)."""
        filename = "xml_2212L_blastp_001.xml"
        datafile = os.path.join("Blast", filename)
        with open(datafile, "rb") as handle:
            records = NCBIXML.parse(handle)
            record = next(records)
            self.assertRaises(StopIteration, next, records)

        alignments = record.alignments
        self.assertEqual(len(alignments), 212)
        self.assertEqual(record.query_id, "gi|49176427|ref|NP_418280.3|")
        self.assertEqual(sum(len(a.hsps) for a in alignments), 212)

        alignment = alignments[0]
        self.assertEqual(
            alignment.title[:50], "gi|49176427|ref|NP_418280.3| component of Sec-inde"
        )
        self.assertEqual(alignment.length, 103)
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertTrue(hsp.expect, 4.20576e-46)
        self.assertEqual(
            hsp.query[:75],
            "MRLCLIIIYHRGTCMGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTA",
        )
        self.assertEqual(
            hsp.match[:75],
            "MRLCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTA",
        )
        self.assertEqual(
            hsp.sbjct[:75],
            "MRLCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTA",
        )

        alignment = alignments[1]
        self.assertEqual(
            alignment.title[:50], "gi|15804428|ref|NP_290468.1| twin arginine translo"
        )
        self.assertEqual(alignment.length, 103)
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertTrue(hsp.expect, 2.72609e-45)
        self.assertEqual(
            hsp.query[:75],
            "MRLCLIIIYHRGTCMGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTA",
        )
        self.assertEqual(
            hsp.match[:75],
            "MRLCLIIIYHR TCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTA",
        )
        self.assertEqual(
            hsp.sbjct[:75],
            "MRLCLIIIYHRXTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTA",
        )

        alignment = alignments[2]
        self.assertEqual(
            alignment.title[:50], "gi|74314349|ref|YP_312768.1| hypothetical protein "
        )
        self.assertEqual(alignment.length, 103)
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertTrue(hsp.expect, 2.72609e-45)
        self.assertEqual(
            hsp.query[:75],
            "MRLCLIIIYHRGTCMGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTA",
        )
        self.assertEqual(
            hsp.match[:75],
            "MR CLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTA",
        )
        self.assertEqual(
            hsp.sbjct[:75],
            "MRPCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTA",
        )

        alignment = alignments[3]
        self.assertEqual(
            alignment.title[:50], "gi|75256240|ref|ZP_00727918.1| COG1826: Sec-indepe"
        )
        self.assertEqual(alignment.length, 89)
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertTrue(hsp.expect, 6.0872e-37)
        self.assertEqual(
            hsp.query[:75],
            "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQ",
        )
        self.assertEqual(
            hsp.match[:75],
            "MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQ",
        )
        self.assertEqual(
            hsp.sbjct[:75],
            "MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQ",
        )

        alignment = alignments[4]
        self.assertEqual(
            alignment.title[:50], "gi|148236|gb|AAA67633.1| o261 [Escherichia coli]"
        )
        self.assertEqual(alignment.length, 261)
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertTrue(hsp.expect, 6.74582e-28)
        self.assertEqual(
            hsp.query[:75],
            "FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDA",
        )
        self.assertEqual(
            hsp.match[:75],
            "FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDA",
        )
        self.assertEqual(
            hsp.sbjct[:75],
            "FGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDA",
        )

        alignment = alignments[5]
        self.assertEqual(
            alignment.title[:50], "gi|29143650|ref|NP_806992.1| sec-independent prote"
        )
        self.assertEqual(alignment.length, 84)
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertTrue(hsp.expect, 4.37251e-27)
        self.assertEqual(
            hsp.query[:75],
            "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQ",
        )
        self.assertEqual(
            hsp.match[:75],
            "MGGISIWQLLI+AVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDD+ KQDKTSQDADFTAK+IADKQ      +",
        )
        self.assertEqual(
            hsp.sbjct[:75],
            "MGGISIWQLLIVAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDDAKQDKTSQDADFTAKSIADKQG-----E",
        )

        alignment = alignments[6]
        self.assertEqual(
            alignment.title[:50], "gi|49609685|emb|CAG73118.1| sec-independent protei"
        )
        self.assertEqual(alignment.length, 86)
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertTrue(hsp.expect, 4.10205e-17)
        self.assertEqual(
            hsp.query[:75],
            "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEP--KQDKTSQDADFTAKTIADKQADTNQ",
        )
        self.assertEqual(
            hsp.match[:75],
            "MGGIS+W LLIIAVIV+LLFGT KL ++GSDLGASIKGFKKAM DD+P    DK   DADF+ K+IAD Q+D   ",
        )
        self.assertEqual(
            hsp.sbjct[:75],
            "MGGISLWNLLIIAVIVILLFGTNKLRTLGSDLGASIKGFKKAMGDDQPSTNADKAQPDADFSTKSIADNQSD---",
        )

        alignment = alignments[7]
        self.assertEqual(
            alignment.title[:50], "gi|37528238|ref|NP_931583.1| Sec-independent prote"
        )
        self.assertEqual(alignment.length, 86)
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertTrue(hsp.expect, 2.25087e-15)
        self.assertEqual(
            hsp.query[:75],
            "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDD-EPKQ-DKTSQDADFTAKTIADKQADTNQ",
        )
        self.assertEqual(
            hsp.match[:75],
            "MGGISIWQLLIIAVIVVLLFGT KL ++GSDLGASIKGFKKA+ DD +P+Q  KTS DADF  K I +KQ+    ",
        )
        self.assertEqual(
            hsp.sbjct[:75],
            "MGGISIWQLLIIAVIVVLLFGTNKLRTLGSDLGASIKGFKKAIGDDNQPQQAQKTSSDADFETKNITEKQS----",
        )

        alignment = alignments[8]
        self.assertEqual(
            alignment.title[:50], "gi|59710656|ref|YP_203432.1| Sec-independent prote"
        )
        self.assertEqual(alignment.length, 82)
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertTrue(hsp.expect, 5.01441e-15)
        self.assertEqual(
            hsp.query[:75],
            "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQ",
        )
        self.assertEqual(
            hsp.match[:75],
            "MGGISIWQLLIIAVI+VLLFGTKKL  +GSDLG+++KGFKKA+S+DEP ++   +DADF  + +  K+A+T ++Q",
        )
        self.assertEqual(
            hsp.sbjct[:75],
            "MGGISIWQLLIIAVIIVLLFGTKKLRGVGSDLGSAVKGFKKAISEDEPAKE-AKKDADFVPQNLEKKEAETVEKQ",
        )

        alignment = alignments[9]
        self.assertEqual(
            alignment.title[:50], "gi|54307340|ref|YP_128360.1| putative TatA protein"
        )
        self.assertEqual(alignment.length, 87)
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertTrue(hsp.expect, 7.2408e-14)
        self.assertEqual(
            hsp.query[:75],
            "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDE--PKQDKTSQDADFTAKTIADKQADTNQ",
        )
        self.assertEqual(
            hsp.match[:75],
            "MGGISIWQLLIIA+I+VLLFGTKKL S+G DLG+++KGFKKA+ D+E   K+D T  DADF  KT++ ++  +  ",
        )
        self.assertEqual(
            hsp.sbjct[:75],
            "MGGISIWQLLIIALIIVLLFGTKKLRSLGGDLGSAVKGFKKAIGDEELTVKKDNTEADADFEQKTLSKEEQQSED",
        )

        alignment = alignments[10]
        self.assertEqual(
            alignment.title[:50], "gi|45437890|gb|AAS63439.1| Sec-independent protein"
        )
        self.assertEqual(alignment.length, 88)
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertTrue(hsp.expect, 1.61308e-13)
        self.assertEqual(
            hsp.query[:75],
            "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDE----PKQDKTSQDADFTAKTIADKQADT",
        )
        self.assertEqual(
            hsp.match[:75],
            "MG I   QLLIIAVIVVLLFGT KL ++GSDLGASIKGFKKAM DD        DKTS DADF AK+I +KQ   ",
        )
        self.assertEqual(
            hsp.sbjct[:75],
            "MGSIGWAQLLIIAVIVVLLFGTNKLRTLGSDLGASIKGFKKAMGDDSQTPPTNVDKTSNDADF-AKSITEKQ---",
        )

        alignment = alignments[11]
        self.assertEqual(
            alignment.title[:50], "gi|75856473|ref|ZP_00764101.1| COG1826: Sec-indepe"
        )
        self.assertEqual(alignment.length, 81)
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertTrue(hsp.expect, 2.10675e-13)
        self.assertEqual(
            hsp.query[:75],
            "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQ",
        )
        self.assertEqual(
            hsp.match[:75],
            "MGGIS+WQLLIIAVIVVLLFGTKKL  IG DLG ++KGFKKAMS+DEP   K  +DADF  K++ ++Q    +++",
        )
        self.assertEqual(
            hsp.sbjct[:75],
            "MGGISVWQLLIIAVIVVLLFGTKKLRGIGGDLGGAVKGFKKAMSEDEPA--KNDKDADFEPKSLEEQQ----KKE",
        )

        alignment = alignments[12]
        self.assertEqual(
            alignment.title[:50], "gi|75829371|ref|ZP_00758676.1| COG1826: Sec-indepe"
        )
        self.assertEqual(alignment.length, 82)
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertTrue(hsp.expect, 2.7515e-13)
        self.assertEqual(
            hsp.query[:75],
            "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQ",
        )
        self.assertEqual(
            hsp.match[:75],
            "MGGISIWQLLIIAVIVVLLFGTKKL  IGSDLG+++KGFKKAMS++E       +DADF  K         N EQ",
        )
        self.assertEqual(
            hsp.sbjct[:75],
            "MGGISIWQLLIIAVIVVLLFGTKKLRGIGSDLGSAVKGFKKAMSEEESNSAANQKDADFETK---------NLEQ",
        )

        alignment = alignments[13]
        self.assertEqual(
            alignment.title[:50], "gi|75820019|ref|ZP_00750077.1| COG1826: Sec-indepe"
        )
        self.assertEqual(alignment.length, 82)
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertTrue(hsp.expect, 3.59357e-13)
        self.assertEqual(
            hsp.query[:75],
            "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQ",
        )
        self.assertEqual(
            hsp.match[:75],
            "MGGISIWQLLIIAVIVVLLFGTKKL  IGSDLG+++KGFKKAMS++E       +DADF  K         N EQ",
        )
        self.assertEqual(
            hsp.sbjct[:75],
            "MGGISIWQLLIIAVIVVLLFGTKKLRGIGSDLGSAVKGFKKAMSEEESNSAANQKDADFETK---------NLEQ",
        )

        alignment = alignments[14]
        self.assertEqual(
            alignment.title[:50], "gi|28896872|ref|NP_796477.1| TatA protein [Vibrio "
        )
        self.assertEqual(alignment.length, 81)
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertTrue(hsp.expect, 2.32928e-12)
        self.assertEqual(
            hsp.query[:75],
            "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQ",
        )
        self.assertEqual(
            hsp.match[:75],
            "MGGIS+WQLLIIAVIVVLLFGTKKL  IG DLG+++KGFKKAMSD++    K  +DADF  K++  +Q    Q++",
        )
        self.assertEqual(
            hsp.sbjct[:75],
            "MGGISVWQLLIIAVIVVLLFGTKKLRGIGGDLGSAVKGFKKAMSDED--SAKNEKDADFEPKSLEKQQ----QKE",
        )

        alignment = alignments[15]
        self.assertEqual(
            alignment.title[:50], "gi|27364353|ref|NP_759881.1| Sec-independent prote"
        )
        self.assertEqual(alignment.length, 78)
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertTrue(hsp.expect, 3.97316e-12)
        self.assertEqual(
            hsp.query[:75],
            "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQ",
        )
        self.assertEqual(
            hsp.match[:75],
            "MGGISIWQLLIIAVIVVLLFGTKKL  IGSDLG +IKGFKKAM+++E ++    +DADF  K++     +   +Q",
        )
        self.assertEqual(
            hsp.sbjct[:75],
            "MGGISIWQLLIIAVIVVLLFGTKKLRGIGSDLGGAIKGFKKAMNEEESEK----KDADFEPKSL-----EQQSKQ",
        )

        alignment = alignments[16]
        self.assertEqual(
            alignment.title[:50], "gi|37678364|ref|NP_932973.1| Sec-independent prote"
        )
        self.assertEqual(alignment.length, 78)
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertTrue(hsp.expect, 3.97316e-12)
        self.assertEqual(
            hsp.query[:75],
            "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQ",
        )
        self.assertEqual(
            hsp.match[:75],
            "MGGISIWQLLIIAVIVVLLFGTKKL  IGSDLG +IKGFKKAM+++E ++    +DADF  K++     +   +Q",
        )
        self.assertEqual(
            hsp.sbjct[:75],
            "MGGISIWQLLIIAVIVVLLFGTKKLRGIGSDLGGAIKGFKKAMNEEESEK----KDADFEPKSL-----EQQNKQ",
        )

        alignment = alignments[17]
        self.assertEqual(
            alignment.title[:50], "gi|71277787|ref|YP_266931.1| Sec-independent prote"
        )
        self.assertEqual(alignment.length, 85)
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertTrue(hsp.expect, 3.97316e-12)
        self.assertEqual(
            hsp.query[:75],
            "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQ",
        )
        self.assertEqual(
            hsp.match[:75],
            "MGGI IWQL+I+AVIVVLLFGTKKL +IG DLG++IKGFK A+ +D  K+ K S  A+ T+ T+AD    T +E ",
        )
        self.assertEqual(
            hsp.sbjct[:75],
            "MGGIGIWQLVIVAVIVVLLFGTKKLRNIGGDLGSAIKGFKSAIGED--KEQKNS--AEKTSDTLADSSKSTTEEV",
        )

        alignment = alignments[18]
        self.assertEqual(
            alignment.title[:50], "gi|68541995|ref|ZP_00581733.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 79)
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertTrue(hsp.expect, 2.57533e-11)
        self.assertEqual(
            hsp.query[:75],
            "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKT---IADKQADTN",
        )
        self.assertEqual(
            hsp.match[:75],
            "MGGISIWQLLI+A+IVVLLFGTKKL S+G DLG ++KGFK AMS +E K+     +A  TA+T     +K+ ++N",
        )
        self.assertEqual(
            hsp.sbjct[:75],
            "MGGISIWQLLIVALIVVLLFGTKKLRSLGGDLGGAVKGFKNAMSSEEDKKALEDTEAAKTAQTTQQATEKKPESN",
        )

        alignment = alignments[19]
        self.assertEqual(
            alignment.title[:50], "gi|77813363|ref|ZP_00812641.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 79)
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertTrue(hsp.expect, 2.57533e-11)
        self.assertEqual(
            hsp.query[:75],
            "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKT---IADKQADTN",
        )
        self.assertEqual(
            hsp.match[:75],
            "MGGISIWQLLIIA+IVVLLFGTKKL S+G DLG ++KGFK AMS +E K+     +A  TA+T     +K+ ++N",
        )
        self.assertEqual(
            hsp.sbjct[:75],
            "MGGISIWQLLIIALIVVLLFGTKKLRSLGGDLGGAVKGFKNAMSSEEDKKALEDTEAAKTAQTTQQATEKKPESN",
        )

        alignment = alignments[20]
        self.assertEqual(
            alignment.title[:50], "gi|52306607|gb|AAU37107.1| TatA protein [Mannheimi"
        )
        self.assertEqual(alignment.length, 75)
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertTrue(hsp.expect, 3.36348e-11)
        self.assertEqual(
            hsp.query[:75],
            "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQ",
        )
        self.assertEqual(
            hsp.match[:75],
            "MGGISIWQLLII  I+VLLFGTKKL ++G+DLG S+KGFKKAM++DEPK      DA+F +    D+ A    E+",
        )
        self.assertEqual(
            hsp.sbjct[:75],
            "MGGISIWQLLIIVAIIVLLFGTKKLRTLGTDLGESVKGFKKAMNEDEPK------DAEFKSLN-KDESATAGSEK",
        )

        alignment = alignments[21]
        self.assertEqual(
            alignment.title[:50], "gi|75429751|ref|ZP_00732413.1| sec-independent pro"
        )
        self.assertEqual(alignment.length, 74)
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertTrue(hsp.expect, 3.36348e-11)
        self.assertEqual(
            hsp.query[:75],
            "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQ",
        )
        self.assertEqual(
            hsp.match[:75],
            "MGGISIWQLLII  IVVLLFGTKKL ++GSDLG S+KGFKKAM+ +EPK      DA+F +   A+  A T +E+",
        )
        self.assertEqual(
            hsp.sbjct[:75],
            "MGGISIWQLLIIVAIVVLLFGTKKLRTLGSDLGESVKGFKKAMA-EEPK------DAEFKSLDKAENTAQTKKEE",
        )

        alignment = alignments[22]
        self.assertEqual(
            alignment.title[:50], "gi|32033565|ref|ZP_00133892.1| COG1826: Sec-indepe"
        )
        self.assertEqual(alignment.length, 76)
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertTrue(hsp.expect, 7.49305e-11)
        self.assertEqual(
            hsp.query[:75],
            "MGGISIWQXXXXXXXXXXXFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQ",
        )
        self.assertEqual(
            hsp.match[:75],
            "MGGISIWQLLII  I+VLLFGTKKL ++G+DLG S+KGFKKAM+DD+      SQ  D + + +  K+A + +++",
        )
        self.assertEqual(
            hsp.sbjct[:75],
            "MGGISIWQLLIIVAIIVLLFGTKKLRTLGTDLGESVKGFKKAMADDK------SQPQDASFEKVEAKEAASTEQK",
        )

        alignment = alignments[23]
        self.assertEqual(
            alignment.title[:50], "gi|12722097|gb|AAK03773.1| unknown [Pasteurella mu"
        )
        self.assertEqual(alignment.length, 76)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[24]
        self.assertEqual(
            alignment.title[:50], "gi|68546478|ref|ZP_00586025.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 77)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[25]
        self.assertEqual(
            alignment.title[:50], "gi|33151888|ref|NP_873241.1| sec-independent prote"
        )
        self.assertEqual(alignment.length, 74)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[26]
        self.assertEqual(
            alignment.title[:50], "gi|24375687|ref|NP_719730.1| Sec-independent prote"
        )
        self.assertEqual(alignment.length, 88)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[27]
        self.assertEqual(
            alignment.title[:50], "gi|71278553|ref|YP_269744.1| Sec-independent prote"
        )
        self.assertEqual(alignment.length, 79)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[28]
        self.assertEqual(
            alignment.title[:50], "gi|69159855|gb|EAN71956.1| Twin-arginine transloca"
        )
        self.assertEqual(alignment.length, 79)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[29]
        self.assertEqual(
            alignment.title[:50], "gi|69949858|ref|ZP_00637822.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 81)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[30]
        self.assertEqual(
            alignment.title[:50], "gi|48863844|ref|ZP_00317737.1| hypothetical protei"
        )
        self.assertEqual(alignment.length, 83)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[31]
        self.assertEqual(
            alignment.title[:50], "gi|77361831|ref|YP_341406.1| twin-arginine translo"
        )
        self.assertEqual(alignment.length, 82)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[32]
        self.assertEqual(
            alignment.title[:50], "gi|67676224|ref|ZP_00472975.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 90)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[33]
        self.assertEqual(
            alignment.title[:50], "gi|74317722|ref|YP_315462.1| twin-arginine translo"
        )
        self.assertEqual(alignment.length, 70)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[34]
        self.assertEqual(
            alignment.title[:50], "gi|77166504|ref|YP_345029.1| Twin-arginine translo"
        )
        self.assertEqual(alignment.length, 90)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[35]
        self.assertEqual(
            alignment.title[:50], "gi|16128610|ref|NP_415160.1| component of Sec-inde"
        )
        self.assertEqual(alignment.length, 67)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[36]
        self.assertEqual(
            alignment.title[:50], "gi|12831974|emb|CAC29147.1| TatA protein [Pseudomo"
        )
        self.assertEqual(alignment.length, 76)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[37]
        self.assertEqual(
            alignment.title[:50], "gi|32029972|ref|ZP_00132908.1| COG1826: Sec-indepe"
        )
        self.assertEqual(alignment.length, 73)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[38]
        self.assertEqual(alignment.title[:50], "gi|455172|gb|AAA24073.1| ORF; putative")
        self.assertEqual(alignment.length, 67)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[39]
        self.assertEqual(alignment.title[:50], "gi|1224007|gb|AAA92108.1| ORF4")
        self.assertEqual(alignment.length, 192)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[40]
        self.assertEqual(
            alignment.title[:50], "gi|68056990|gb|AAX87243.1| Sec-independent protein"
        )
        self.assertEqual(alignment.length, 95)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[41]
        self.assertEqual(
            alignment.title[:50], "gi|56461470|ref|YP_156751.1| Sec-independent prote"
        )
        self.assertEqual(alignment.length, 73)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[42]
        self.assertEqual(
            alignment.title[:50], "gi|76793313|ref|ZP_00775802.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 84)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[43]
        self.assertEqual(
            alignment.title[:50], "gi|42630489|ref|ZP_00156028.1| COG1826: Sec-indepe"
        )
        self.assertEqual(alignment.length, 75)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[44]
        self.assertEqual(
            alignment.title[:50], "gi|1074302|pir||B64145 hypothetical protein HI0187"
        )
        self.assertEqual(alignment.length, 109)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[45]
        self.assertEqual(
            alignment.title[:50], "gi|67641583|ref|ZP_00440359.1| COG1826: Sec-indepe"
        )
        self.assertEqual(alignment.length, 77)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[46]
        self.assertEqual(
            alignment.title[:50], "gi|67545726|ref|ZP_00423646.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 76)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[47]
        self.assertEqual(
            alignment.title[:50], "gi|45435806|gb|AAS61363.1| sec-independent protein"
        )
        self.assertEqual(alignment.length, 85)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[48]
        self.assertEqual(
            alignment.title[:50], "gi|49610761|emb|CAG74206.1| Sec-independent protei"
        )
        self.assertEqual(alignment.length, 65)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[49]
        self.assertEqual(
            alignment.title[:50], "gi|67663266|ref|ZP_00460549.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 76)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[50]
        self.assertEqual(
            alignment.title[:50], "gi|33594634|ref|NP_882278.1| Sec-independent prote"
        )
        self.assertEqual(alignment.length, 75)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[51]
        self.assertEqual(
            alignment.title[:50], "gi|46310681|ref|ZP_00211309.1| COG1826: Sec-indepe"
        )
        self.assertEqual(alignment.length, 76)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[52]
        self.assertEqual(
            alignment.title[:50], "gi|58584031|ref|YP_203047.1| sec-independent prote"
        )
        self.assertEqual(alignment.length, 75)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[53]
        self.assertEqual(
            alignment.title[:50], "gi|17429965|emb|CAD16649.1| PROBABLE SIGNAL PEPTID"
        )
        self.assertEqual(alignment.length, 85)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[54]
        self.assertEqual(
            alignment.title[:50], "gi|47573371|ref|ZP_00243410.1| COG1826: Sec-indepe"
        )
        self.assertEqual(alignment.length, 76)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[55]
        self.assertEqual(
            alignment.title[:50], "gi|16273687|ref|NP_438355.1| Sec-independent prote"
        )
        self.assertEqual(alignment.length, 89)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[56]
        self.assertEqual(
            alignment.title[:50], "gi|73542784|ref|YP_297304.1| Twin-arginine translo"
        )
        self.assertEqual(alignment.length, 73)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[57]
        self.assertEqual(
            alignment.title[:50], "gi|26987777|ref|NP_743202.1| Sec-independent prote"
        )
        self.assertEqual(alignment.length, 77)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[58]
        self.assertEqual(
            alignment.title[:50], "gi|29142636|ref|NP_805978.1| sec-independent prote"
        )
        self.assertEqual(alignment.length, 67)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[59]
        self.assertEqual(
            alignment.title[:50], "gi|18389921|gb|AAL68797.1| TatA [Ralstonia eutroph"
        )
        self.assertEqual(alignment.length, 77)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[60]
        self.assertEqual(
            alignment.title[:50], "gi|48781637|ref|ZP_00278228.1| COG1826: Sec-indepe"
        )
        self.assertEqual(alignment.length, 79)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[61]
        self.assertEqual(
            alignment.title[:50], "gi|77456610|ref|YP_346115.1| Twin-arginine translo"
        )
        self.assertEqual(alignment.length, 92)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[62]
        self.assertEqual(
            alignment.title[:50], "gi|1684735|emb|CAA98158.1| ORF57 protein [Pseudomo"
        )
        self.assertEqual(alignment.length, 57)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[63]
        self.assertEqual(
            alignment.title[:50], "gi|56476124|ref|YP_157713.1| Sec-independent prote"
        )
        self.assertEqual(alignment.length, 75)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[64]
        self.assertEqual(
            alignment.title[:50], "gi|34496078|ref|NP_900293.1| Sec-independent prote"
        )
        self.assertEqual(alignment.length, 68)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[65]
        self.assertEqual(
            alignment.title[:50], "gi|67848115|ref|ZP_00503233.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 83)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[66]
        self.assertEqual(
            alignment.title[:50], "gi|26991692|ref|NP_747117.1| Sec-independent prote"
        )
        self.assertEqual(alignment.length, 90)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[67]
        self.assertEqual(
            alignment.title[:50], "gi|15601293|ref|NP_232924.1| tatA protein [Vibrio "
        )
        self.assertEqual(alignment.length, 78)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[68]
        self.assertEqual(
            alignment.title[:50], "gi|66770480|ref|YP_245242.1| sec-independent prote"
        )
        self.assertEqual(alignment.length, 75)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[69]
        self.assertEqual(
            alignment.title[:50], "gi|53804435|ref|YP_113945.1| Sec-independent prote"
        )
        self.assertEqual(alignment.length, 70)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[70]
        self.assertEqual(
            alignment.title[:50], "gi|75825357|ref|ZP_00754793.1| COG1826: Sec-indepe"
        )
        self.assertEqual(alignment.length, 80)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[71]
        self.assertEqual(
            alignment.title[:50], "gi|71908987|ref|YP_286574.1| Twin-arginine translo"
        )
        self.assertEqual(alignment.length, 76)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[72]
        self.assertEqual(
            alignment.title[:50], "gi|68526571|gb|EAN49542.1| Twin-arginine transloca"
        )
        self.assertEqual(alignment.length, 77)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[73]
        self.assertEqual(
            alignment.title[:50], "gi|71736448|ref|YP_272670.1| sec-independent prote"
        )
        self.assertEqual(alignment.length, 91)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[74]
        self.assertEqual(
            alignment.title[:50], "gi|56460344|ref|YP_155625.1| Sec-independent prote"
        )
        self.assertEqual(alignment.length, 72)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[75]
        self.assertEqual(
            alignment.title[:50], "gi|68214708|ref|ZP_00566522.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 72)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[76]
        self.assertEqual(
            alignment.title[:50], "gi|30248650|ref|NP_840720.1| mttA/Hcf106 family [N"
        )
        self.assertEqual(alignment.length, 76)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[77]
        self.assertEqual(
            alignment.title[:50], "gi|75822907|ref|ZP_00752458.1| COG1826: Sec-indepe"
        )
        self.assertEqual(alignment.length, 78)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[78]
        self.assertEqual(
            alignment.title[:50], "gi|70733926|ref|YP_257566.1| sec-independent prote"
        )
        self.assertEqual(alignment.length, 93)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[79]
        self.assertEqual(
            alignment.title[:50], "gi|63254358|gb|AAY35454.1| Twin-arginine transloca"
        )
        self.assertEqual(alignment.length, 91)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[80]
        self.assertEqual(
            alignment.title[:50], "gi|73354814|gb|AAZ75668.1| TatA [Pseudomonas syrin"
        )
        self.assertEqual(alignment.length, 91)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[81]
        self.assertEqual(
            alignment.title[:50], "gi|50083761|ref|YP_045271.1| Sec-independent prote"
        )
        self.assertEqual(alignment.length, 78)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[82]
        self.assertEqual(
            alignment.title[:50], "gi|71548504|ref|ZP_00668728.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 77)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[83]
        self.assertEqual(
            alignment.title[:50], "gi|55247002|gb|EAL42253.1| ENSANGP00000028218 [Ano"
        )
        self.assertEqual(alignment.length, 53)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[84]
        self.assertEqual(
            alignment.title[:50], "gi|50084688|ref|YP_046198.1| Sec-independent prote"
        )
        self.assertEqual(alignment.length, 71)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[85]
        self.assertEqual(
            alignment.title[:50], "gi|28872267|ref|NP_794886.1| sec-independent prote"
        )
        self.assertEqual(alignment.length, 91)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[86]
        self.assertEqual(
            alignment.title[:50], "gi|49082486|gb|AAT50643.1| PA5068 [synthetic const"
        )
        self.assertEqual(alignment.length, 83)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[87]
        self.assertEqual(
            alignment.title[:50], "gi|53726598|ref|ZP_00141543.2| COG1826: Sec-indepe"
        )
        self.assertEqual(alignment.length, 82)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[88]
        self.assertEqual(
            alignment.title[:50], "gi|68213616|ref|ZP_00565447.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 54)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[89]
        self.assertEqual(
            alignment.title[:50], "gi|74023810|ref|ZP_00694377.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 79)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[90]
        self.assertEqual(
            alignment.title[:50], "gi|71066554|ref|YP_265281.1| twin-arginine translo"
        )
        self.assertEqual(alignment.length, 87)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[91]
        self.assertEqual(
            alignment.title[:50], "gi|15611372|ref|NP_223023.1| hypothetical protein "
        )
        self.assertEqual(alignment.length, 79)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[92]
        self.assertEqual(
            alignment.title[:50], "gi|13471183|ref|NP_102752.1| sec-independent prote"
        )
        self.assertEqual(alignment.length, 73)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[93]
        self.assertEqual(
            alignment.title[:50], "gi|42523995|ref|NP_969375.1| twin-arginine-depende"
        )
        self.assertEqual(alignment.length, 81)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[94]
        self.assertEqual(
            alignment.title[:50], "gi|67158086|ref|ZP_00419176.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 85)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[95]
        self.assertEqual(
            alignment.title[:50], "gi|15644948|ref|NP_207118.1| conserved hypothetica"
        )
        self.assertEqual(alignment.length, 79)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[96]
        self.assertEqual(
            alignment.title[:50], "gi|13277311|emb|CAC34414.1| putative TatA protein "
        )
        self.assertEqual(alignment.length, 61)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[97]
        self.assertEqual(
            alignment.title[:50], "gi|54298906|ref|YP_125275.1| Putative TatA protein"
        )
        self.assertEqual(alignment.length, 61)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[98]
        self.assertEqual(
            alignment.title[:50], "gi|71363513|ref|ZP_00654157.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 94)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[99]
        self.assertEqual(
            alignment.title[:50], "gi|71362217|ref|ZP_00653377.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 80)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[100]
        self.assertEqual(
            alignment.title[:50], "gi|27379862|ref|NP_771391.1| hypothetical protein "
        )
        self.assertEqual(alignment.length, 77)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[101]
        self.assertEqual(
            alignment.title[:50], "gi|39935914|ref|NP_948190.1| putative sec-independ"
        )
        self.assertEqual(alignment.length, 78)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[102]
        self.assertEqual(
            alignment.title[:50], "gi|17935600|ref|NP_532390.1| SEC-independent prote"
        )
        self.assertEqual(alignment.length, 70)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[103]
        self.assertEqual(
            alignment.title[:50], "gi|62289827|ref|YP_221620.1| Sec-independent prote"
        )
        self.assertEqual(alignment.length, 72)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[104]
        self.assertEqual(
            alignment.title[:50], "gi|23347697|gb|AAN29810.1| Sec-independent protein"
        )
        self.assertEqual(alignment.length, 80)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[105]
        self.assertEqual(
            alignment.title[:50], "gi|75675971|ref|YP_318392.1| twin-arginine translo"
        )
        self.assertEqual(alignment.length, 78)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[106]
        self.assertEqual(
            alignment.title[:50], "gi|69928230|ref|ZP_00625391.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 79)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[107]
        self.assertEqual(
            alignment.title[:50], "gi|77689454|ref|ZP_00804635.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 79)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[108]
        self.assertEqual(
            alignment.title[:50], "gi|77743614|ref|ZP_00812071.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 78)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[109]
        self.assertEqual(
            alignment.title[:50], "gi|71066141|ref|YP_264868.1| twin-arginine translo"
        )
        self.assertEqual(alignment.length, 89)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[110]
        self.assertEqual(
            alignment.title[:50], "gi|28199457|ref|NP_779771.1| SEC-independent prote"
        )
        self.assertEqual(alignment.length, 71)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[111]
        self.assertEqual(
            alignment.title[:50], "gi|15837166|ref|NP_297854.1| hypothetical protein "
        )
        self.assertEqual(alignment.length, 71)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[112]
        self.assertEqual(
            alignment.title[:50], "gi|15074462|emb|CAC46108.1| HYPOTHETICAL TRANSMEMB"
        )
        self.assertEqual(alignment.length, 68)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[113]
        self.assertEqual(
            alignment.title[:50], "gi|27462871|gb|AAO15625.1| Sec-independent protein"
        )
        self.assertEqual(alignment.length, 63)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[114]
        self.assertEqual(
            alignment.title[:50], "gi|35211273|dbj|BAC88652.1| gsl0711 [Gloeobacter v"
        )
        self.assertEqual(alignment.length, 72)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[115]
        self.assertEqual(
            alignment.title[:50], "gi|34482347|emb|CAE09348.1| hypothetical protein ["
        )
        self.assertEqual(alignment.length, 80)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[116]
        self.assertEqual(
            alignment.title[:50], "gi|32262257|gb|AAP77305.1| component of Sec-indepe"
        )
        self.assertEqual(alignment.length, 82)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[117]
        self.assertEqual(
            alignment.title[:50], "gi|76261408|ref|ZP_00769019.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 62)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[118]
        self.assertEqual(
            alignment.title[:50], "gi|69933726|ref|ZP_00628928.1| sec-independent tra"
        )
        self.assertEqual(alignment.length, 159)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[119]
        self.assertEqual(
            alignment.title[:50], "gi|15605662|ref|NP_213037.1| hypothetical protein "
        )
        self.assertEqual(alignment.length, 59)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[120]
        self.assertEqual(
            alignment.title[:50], "gi|68538777|ref|ZP_00578553.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 77)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[121]
        self.assertEqual(
            alignment.title[:50], "gi|68136098|ref|ZP_00544086.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 130)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[122]
        self.assertEqual(
            alignment.title[:50], "gi|20259265|gb|AAM14368.1| putative Tha4 protein ["
        )
        self.assertEqual(alignment.length, 147)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[123]
        self.assertEqual(
            alignment.title[:50], "gi|75910646|ref|YP_324942.1| Twin-arginine translo"
        )
        self.assertEqual(alignment.length, 90)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[124]
        self.assertEqual(
            alignment.title[:50], "gi|39982657|gb|AAR34117.1| twin-arginine transloca"
        )
        self.assertEqual(alignment.length, 57)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[125]
        self.assertEqual(
            alignment.title[:50], "gi|33635687|emb|CAE22011.1| mttA/Hcf106 family [Pr"
        )
        self.assertEqual(alignment.length, 91)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[126]
        self.assertEqual(
            alignment.title[:50], "gi|76791934|ref|ZP_00774438.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 68)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[127]
        self.assertEqual(
            alignment.title[:50], "gi|23129516|ref|ZP_00111343.1| COG1826: Sec-indepe"
        )
        self.assertEqual(alignment.length, 91)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[128]
        self.assertEqual(
            alignment.title[:50], "gi|48764199|ref|ZP_00268751.1| COG1826: Sec-indepe"
        )
        self.assertEqual(alignment.length, 96)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[129]
        self.assertEqual(
            alignment.title[:50], "gi|15677995|ref|NP_273645.1| hypothetical protein "
        )
        self.assertEqual(alignment.length, 67)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[130]
        self.assertEqual(
            alignment.title[:50], "gi|50917153|ref|XP_468973.1| putative sec-independ"
        )
        self.assertEqual(alignment.length, 170)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[131]
        self.assertEqual(
            alignment.title[:50], "gi|16329622|ref|NP_440350.1| hypothetical protein "
        )
        self.assertEqual(alignment.length, 126)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[132]
        self.assertEqual(
            alignment.title[:50], "gi|71083667|ref|YP_266387.1| Twin-arginine translo"
        )
        self.assertEqual(alignment.length, 66)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[133]
        self.assertEqual(
            alignment.title[:50], "gi|17130190|dbj|BAB72802.1| asl0845 [Nostoc sp. PC"
        )
        self.assertEqual(alignment.length, 90)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[134]
        self.assertEqual(
            alignment.title[:50], "gi|68246031|gb|EAN28138.1| Twin-arginine transloca"
        )
        self.assertEqual(alignment.length, 69)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[135]
        self.assertEqual(
            alignment.title[:50], "gi|15604583|ref|NP_221101.1| hypothetical protein "
        )
        self.assertEqual(alignment.length, 54)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[136]
        self.assertEqual(
            alignment.title[:50], "gi|77685166|ref|ZP_00800574.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 69)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[137]
        self.assertEqual(
            alignment.title[:50], "gi|39985226|gb|AAR36581.1| twin-arginine transloca"
        )
        self.assertEqual(alignment.length, 78)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[138]
        self.assertEqual(
            alignment.title[:50], "gi|1825636|gb|AAB42258.1| Hypothetical protein ZK3"
        )
        self.assertEqual(alignment.length, 312)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[139]
        self.assertEqual(
            alignment.title[:50], "gi|65321915|ref|ZP_00394874.1| COG5386: Cell surfa"
        )
        self.assertEqual(alignment.length, 237)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[140]
        self.assertEqual(
            alignment.title[:50], "gi|30022625|ref|NP_834256.1| Cell surface protein "
        )
        self.assertEqual(alignment.length, 237)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[141]
        self.assertEqual(
            alignment.title[:50], "gi|55623442|ref|XP_517520.1| PREDICTED: similar to"
        )
        self.assertEqual(alignment.length, 234)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[142]
        self.assertEqual(
            alignment.title[:50], "gi|75762866|ref|ZP_00742681.1| Cell surface protei"
        )
        self.assertEqual(alignment.length, 237)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[143]
        self.assertEqual(
            alignment.title[:50], "gi|22945598|gb|AAN10511.1| CG18497-PC, isoform C ["
        )
        self.assertEqual(alignment.length, 5476)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[144]
        self.assertEqual(
            alignment.title[:50], "gi|10727420|gb|AAF51534.2| CG18497-PB, isoform B ["
        )
        self.assertEqual(alignment.length, 5533)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[145]
        self.assertEqual(
            alignment.title[:50], "gi|10727421|gb|AAF51535.2| CG18497-PA, isoform A ["
        )
        self.assertEqual(alignment.length, 5560)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[146]
        self.assertEqual(
            alignment.title[:50], "gi|71481981|ref|ZP_00661682.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 69)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[147]
        self.assertEqual(
            alignment.title[:50], "gi|71150623|ref|ZP_00649545.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 81)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[148]
        self.assertEqual(
            alignment.title[:50], "gi|20151563|gb|AAM11141.1| LD15253p [Drosophila me"
        )
        self.assertEqual(alignment.length, 1521)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[149]
        self.assertEqual(
            alignment.title[:50], "gi|6979936|gb|AAF34661.1| split ends long isoform "
        )
        self.assertEqual(alignment.length, 5554)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[150]
        self.assertEqual(
            alignment.title[:50], "gi|6467825|gb|AAF13218.1| Spen RNP motif protein l"
        )
        self.assertEqual(alignment.length, 5533)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[151]
        self.assertEqual(
            alignment.title[:50], "gi|61102013|ref|ZP_00377467.1| hypothetical protei"
        )
        self.assertEqual(alignment.length, 80)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[152]
        self.assertEqual(
            alignment.title[:50], "gi|68056232|ref|ZP_00540361.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 68)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[153]
        self.assertEqual(
            alignment.title[:50], "gi|68190120|gb|EAN04781.1| Twin-arginine transloca"
        )
        self.assertEqual(alignment.length, 71)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[154]
        self.assertEqual(
            alignment.title[:50], "gi|15605663|ref|NP_213038.1| hypothetical protein "
        )
        self.assertEqual(alignment.length, 77)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[155]
        self.assertEqual(
            alignment.title[:50], "gi|60493413|emb|CAH08199.1| aerotolerance-related "
        )
        self.assertEqual(alignment.length, 238)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[156]
        self.assertEqual(
            alignment.title[:50], "gi|50877510|emb|CAG37350.1| related to Sec-indepen"
        )
        self.assertEqual(alignment.length, 84)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[157]
        self.assertEqual(
            alignment.title[:50], "gi|42739647|gb|AAS43573.1| conserved domain protei"
        )
        self.assertEqual(alignment.length, 236)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[158]
        self.assertEqual(
            alignment.title[:50], "gi|53713708|ref|YP_099700.1| conserved hypothetica"
        )
        self.assertEqual(alignment.length, 238)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[159]
        self.assertEqual(
            alignment.title[:50], "gi|33860901|ref|NP_892462.1| mttA/Hcf106 family [P"
        )
        self.assertEqual(alignment.length, 96)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[160]
        self.assertEqual(
            alignment.title[:50], "gi|48851224|ref|ZP_00305466.1| COG1826: Sec-indepe"
        )
        self.assertEqual(alignment.length, 83)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[161]
        self.assertEqual(
            alignment.title[:50], "gi|67938449|ref|ZP_00530974.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 69)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[162]
        self.assertEqual(
            alignment.title[:50], "gi|45657833|ref|YP_001919.1| sec-independent prote"
        )
        self.assertEqual(alignment.length, 90)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[163]
        self.assertEqual(
            alignment.title[:50], "gi|57238048|ref|YP_179297.1| twin-arginine translo"
        )
        self.assertEqual(alignment.length, 79)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[164]
        self.assertEqual(
            alignment.title[:50], "gi|56962648|ref|YP_174374.1| sec-independent prote"
        )
        self.assertEqual(alignment.length, 63)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[165]
        self.assertEqual(
            alignment.title[:50], "gi|33239734|ref|NP_874676.1| Sec-independent prote"
        )
        self.assertEqual(alignment.length, 84)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[166]
        self.assertEqual(
            alignment.title[:50], "gi|21674434|ref|NP_662499.1| Sec-independent prote"
        )
        self.assertEqual(alignment.length, 67)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[167]
        self.assertEqual(
            alignment.title[:50], "gi|39968009|ref|XP_365395.1| hypothetical protein "
        )
        self.assertEqual(alignment.length, 823)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[168]
        self.assertEqual(
            alignment.title[:50], "gi|4877986|gb|AAD31523.1| THA9 [Zea mays]"
        )
        self.assertEqual(alignment.length, 169)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[169]
        self.assertEqual(
            alignment.title[:50], "gi|67934419|ref|ZP_00527476.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 56)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[170]
        self.assertEqual(
            alignment.title[:50], "gi|42523658|ref|NP_969038.1| twin-argine protein t"
        )
        self.assertEqual(alignment.length, 79)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[171]
        self.assertEqual(
            alignment.title[:50], "gi|71546080|ref|ZP_00666945.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 73)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[172]
        self.assertEqual(
            alignment.title[:50], "gi|68002197|ref|ZP_00534828.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 60)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[173]
        self.assertEqual(
            alignment.title[:50], "gi|67481641|ref|XP_656170.1| hypothetical protein "
        )
        self.assertEqual(alignment.length, 434)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[174]
        self.assertEqual(
            alignment.title[:50], "gi|50935447|ref|XP_477251.1| putative Calreticulin"
        )
        self.assertEqual(alignment.length, 424)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[175]
        self.assertEqual(
            alignment.title[:50], "gi|50978634|ref|NP_001003013.1| acidic (leucine-ri"
        )
        self.assertEqual(alignment.length, 249)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[176]
        self.assertEqual(
            alignment.title[:50], "gi|70936814|ref|XP_739300.1| 40S ribosomal subunit"
        )
        self.assertEqual(alignment.length, 184)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[177]
        self.assertEqual(
            alignment.title[:50], "gi|68075857|ref|XP_679848.1| hypothetical protein "
        )
        self.assertEqual(alignment.length, 340)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[178]
        self.assertEqual(
            alignment.title[:50], "gi|39594005|emb|CAE70115.1| Hypothetical protein C"
        )
        self.assertEqual(alignment.length, 192)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[179]
        self.assertEqual(
            alignment.title[:50], "gi|66809957|ref|XP_638702.1| hypothetical protein "
        )
        self.assertEqual(alignment.length, 721)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[180]
        self.assertEqual(
            alignment.title[:50], "gi|68550463|ref|ZP_00589911.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 69)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[181]
        self.assertEqual(
            alignment.title[:50], "gi|51473916|ref|YP_067673.1| TatA/E-like Sec-indep"
        )
        self.assertEqual(alignment.length, 53)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[182]
        self.assertEqual(
            alignment.title[:50], "gi|61857708|ref|XP_612559.1| PREDICTED: similar to"
        )
        self.assertEqual(alignment.length, 236)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[183]
        self.assertEqual(
            alignment.title[:50], "gi|39982651|gb|AAR34111.1| twin-arginine transloca"
        )
        self.assertEqual(alignment.length, 59)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[184]
        self.assertEqual(
            alignment.title[:50], "gi|50877509|emb|CAG37349.1| related to Sec-indepen"
        )
        self.assertEqual(alignment.length, 66)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[185]
        self.assertEqual(
            alignment.title[:50], "gi|52699323|ref|ZP_00340731.1| COG1826: Sec-indepe"
        )
        self.assertEqual(alignment.length, 53)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[186]
        self.assertEqual(
            alignment.title[:50], "gi|62426215|ref|ZP_00381343.1| COG1826: Sec-indepe"
        )
        self.assertEqual(alignment.length, 93)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[187]
        self.assertEqual(
            alignment.title[:50], "gi|11131838|sp|Q9SLY8|CRTC_ORYSA Calreticulin prec"
        )
        self.assertEqual(alignment.length, 424)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[188]
        self.assertEqual(
            alignment.title[:50], "gi|56543690|gb|AAV89844.1| Sec-independent protein"
        )
        self.assertEqual(alignment.length, 87)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[189]
        self.assertEqual(
            alignment.title[:50], "gi|67923730|ref|ZP_00517196.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 95)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[190]
        self.assertEqual(
            alignment.title[:50], "gi|67462585|ref|XP_647954.1| hypothetical protein "
        )
        self.assertEqual(alignment.length, 140)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[191]
        self.assertEqual(
            alignment.title[:50], "gi|51970620|dbj|BAD44002.1| unknown protein [Arabi"
        )
        self.assertEqual(alignment.length, 784)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[192]
        self.assertEqual(
            alignment.title[:50], "gi|34581241|ref|ZP_00142721.1| hypothetical protei"
        )
        self.assertEqual(alignment.length, 53)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[193]
        self.assertEqual(
            alignment.title[:50], "gi|4877984|gb|AAD31522.1| THA4 [Zea mays]"
        )
        self.assertEqual(alignment.length, 170)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[194]
        self.assertEqual(
            alignment.title[:50], "gi|9757886|dbj|BAB08393.1| unnamed protein product"
        )
        self.assertEqual(alignment.length, 707)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[195]
        self.assertEqual(
            alignment.title[:50], "gi|32422107|ref|XP_331497.1| predicted protein [Ne"
        )
        self.assertEqual(alignment.length, 216)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[196]
        self.assertEqual(
            alignment.title[:50], "gi|68552035|ref|ZP_00591428.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 70)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[197]
        self.assertEqual(
            alignment.title[:50], "gi|68177649|ref|ZP_00550794.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 58)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[198]
        self.assertEqual(
            alignment.title[:50], "gi|67934756|ref|ZP_00527782.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 65)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[199]
        self.assertEqual(
            alignment.title[:50], "gi|42550455|gb|EAA73298.1| hypothetical protein FG"
        )
        self.assertEqual(alignment.length, 297)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[200]
        self.assertEqual(
            alignment.title[:50], "gi|15893083|ref|NP_360797.1| hypothetical protein "
        )
        self.assertEqual(alignment.length, 53)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[201]
        self.assertEqual(
            alignment.title[:50], "gi|57233621|ref|YP_182297.1| twin-arginine translo"
        )
        self.assertEqual(alignment.length, 65)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[202]
        self.assertEqual(
            alignment.title[:50], "gi|75908036|ref|YP_322332.1| Twin-arginine translo"
        )
        self.assertEqual(alignment.length, 56)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[203]
        self.assertEqual(
            alignment.title[:50], "gi|72383453|ref|YP_292808.1| Twin-arginine translo"
        )
        self.assertEqual(alignment.length, 71)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[204]
        self.assertEqual(
            alignment.title[:50], "gi|1666185|emb|CAB04766.1| ORF13(1) [Rhodococcus e"
        )
        self.assertEqual(alignment.length, 98)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[205]
        self.assertEqual(
            alignment.title[:50], "gi|72138252|ref|XP_800288.1| PREDICTED: hypothetic"
        )
        self.assertEqual(alignment.length, 946)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[206]
        self.assertEqual(
            alignment.title[:50], "gi|67923190|ref|ZP_00516678.1| Twin-arginine trans"
        )
        self.assertEqual(alignment.length, 50)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[207]
        self.assertEqual(
            alignment.title[:50], "gi|3329623|gb|AAC26930.1| Hypothetical protein F36"
        )
        self.assertEqual(alignment.length, 335)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[208]
        self.assertEqual(
            alignment.title[:50], "gi|39597929|emb|CAE68621.1| Hypothetical protein C"
        )
        self.assertEqual(alignment.length, 2691)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[209]
        self.assertEqual(
            alignment.title[:50], "gi|68182025|ref|ZP_00555006.1| hypothetical protei"
        )
        self.assertEqual(alignment.length, 438)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[210]
        self.assertEqual(
            alignment.title[:50], "gi|21204492|dbj|BAB95189.1| ebh [Staphylococcus au"
        )
        self.assertEqual(alignment.length, 9904)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        alignment = alignments[211]
        self.assertEqual(
            alignment.title[:50], "gi|39593039|emb|CAE64508.1| Hypothetical protein C"
        )
        self.assertEqual(alignment.length, 960)
        self.assertEqual(len(alignment.hsps), 1)
        self.assertGreater(alignment.hsps[0].expect, E_VALUE_THRESH)

        with open(datafile, "rb") as handle:
            record = NCBIXML.read(handle)

    def test_xml_2212L_blastn_001(self):
        """Parsing BLASTN 2.2.12, gi|1348916|gb|G26684.1|G26684 (xml_2212L_blastn_001)."""
        filename = "xml_2212L_blastn_001.xml"
        datafile = os.path.join("Blast", filename)
        with open(datafile, "rb") as handle:
            records = NCBIXML.parse(handle)
            record = next(records)
            self.assertRaises(StopIteration, next, records)
        alignments = record.alignments

        self.assertEqual(record.query_id, "gi|1348916|gb|G26684.1|G26684")
        self.assertEqual(len(alignments), 2)
        self.assertEqual(sum(len(a.hsps) for a in alignments), 2)
        self.assertEqual(
            alignments[0].title[:50],
            "gi|9950606|gb|AE004854.1| Pseudomonas aeruginosa P",
        )
        self.assertEqual(alignments[0].length, 11884)
        self.assertEqual(len(alignments[0].hsps), 1)
        self.assertGreater(alignments[0].hsps[0].expect, E_VALUE_THRESH)
        self.assertEqual(
            alignments[1].title[:50],
            "gi|15073988|emb|AL591786.1|SME591786 Sinorhizobium",
        )
        self.assertEqual(alignments[1].length, 299350)
        self.assertEqual(len(alignments[1].hsps), 1)
        self.assertGreater(alignments[1].hsps[0].expect, E_VALUE_THRESH)

        with open(datafile, "rb") as handle:
            record = NCBIXML.read(handle)

    def test_xml_2212L_blastx_001(self):
        """Parsing BLASTX 2.2.12, gi|1347369|gb|G25137.1|G25137 (xml_2212L_blastx_001)."""
        filename = "xml_2212L_blastx_001.xml"
        datafile = os.path.join("Blast", filename)

        with open(datafile, "rb") as handle:
            records = NCBIXML.parse(handle)
            record = next(records)
            self.assertRaises(StopIteration, next, records)
        alignments = record.alignments
        self.assertEqual(record.query_id, "gi|1347369|gb|G25137.1|G25137")
        self.assertEqual(len(alignments), 78)
        self.assertEqual(sum(len(a.hsps) for a in alignments), 84)

        with open(datafile, "rb") as handle:
            record = NCBIXML.read(handle)

        hsp = record.alignments[0].hsps[0]
        self.assertEqual(hsp.score, 630.0)
        self.assertEqual(hsp.bits, 247.284)

    def test_xml_2212L_tblastn_001(self):
        """Parsing TBLASTN 2.2.12, gi|729325|sp|P39483|DHG2_BACME (xml_2212L_tblastn_001)."""
        filename = "xml_2212L_tblastn_001.xml"
        datafile = os.path.join("Blast", filename)

        with open(datafile, "rb") as handle:
            records = NCBIXML.parse(handle)
            record = next(records)
            self.assertRaises(StopIteration, next, records)

        alignments = record.alignments
        self.assertEqual(record.query_id, "gi|729325|sp|P39483|DHG2_BACME")
        self.assertEqual(len(alignments), 100)
        self.assertEqual(sum(len(a.hsps) for a in alignments), 127)

        with open(datafile, "rb") as handle:
            record = NCBIXML.read(handle)

    def test_xml_2212L_tblastx_001(self):
        """Parsing TBLASTX 2.2.12, gi|1348853|gb|G26621.1|G26621, BLOSUM80 (xml_2212L_tblastx_001)."""
        filename = "xml_2212L_tblastx_001.xml"
        datafile = os.path.join("Blast", filename)

        with open(datafile, "rb") as handle:
            records = NCBIXML.parse(handle)
            record = next(records)
            self.assertRaises(StopIteration, next, records)

        alignments = record.alignments
        self.assertEqual(record.query_id, "gi|1348853|gb|G26621.1|G26621")
        self.assertEqual(len(alignments), 10)
        self.assertEqual(sum(len(a.hsps) for a in alignments), 102)

        with open(datafile, "rb") as handle:
            record = NCBIXML.read(handle)

    def test_xml_2218_blastp_001(self):
        """Parsing BLASTP 2.2.18+, gi|160837788|ref|NP_075631.2| (xml_2218_blastp_001)."""
        # NOTE - no date in version field, downloaded 2008/05/08

        filename = "xml_2218_blastp_001.xml"
        datafile = os.path.join("Blast", filename)

        with open(datafile, "rb") as handle:
            records = NCBIXML.parse(handle)
            record = next(records)
            self.assertRaises(StopIteration, next, records)

        alignments = record.alignments
        self.assertEqual(record.query_id, "31493")
        self.assertEqual(len(alignments), 10)
        self.assertEqual(sum(len(a.hsps) for a in alignments), 14)
        self.assertEqual(
            alignments[0].title[:50],
            "gi|151942244|gb|EDN60600.1| cytosolic iron-sulfur ",
        )
        self.assertEqual(alignments[0].length, 330)
        self.assertEqual(len(alignments[0].hsps), 1)
        self.assertGreater(alignments[0].hsps[0].expect, E_VALUE_THRESH)
        self.assertEqual(
            alignments[1].title[:50],
            "gi|476059|emb|CAA55606.1| YBR0832 [Saccharomyces c",
        )
        self.assertEqual(alignments[1].length, 535)
        self.assertEqual(len(alignments[1].hsps), 1)
        self.assertGreater(alignments[1].hsps[0].expect, E_VALUE_THRESH)
        self.assertEqual(
            alignments[2].title[:50],
            "gi|6320473|ref|NP_010553.1| Essential protein invo",
        )
        self.assertEqual(alignments[2].length, 330)
        self.assertEqual(len(alignments[2].hsps), 1)
        self.assertGreater(alignments[2].hsps[0].expect, E_VALUE_THRESH)
        self.assertEqual(
            alignments[3].title[:50],
            "gi|61679798|pdb|1R5M|A Chain A, Crystal Structure ",
        )
        self.assertEqual(alignments[3].length, 425)
        self.assertEqual(len(alignments[3].hsps), 1)
        self.assertGreater(alignments[3].hsps[0].expect, E_VALUE_THRESH)
        self.assertEqual(
            alignments[4].title[:50],
            "gi|6319579|ref|NP_009661.1| WD40 repeat-containing",
        )
        self.assertEqual(alignments[4].length, 535)
        self.assertEqual(len(alignments[4].hsps), 1)
        self.assertGreater(alignments[4].hsps[0].expect, E_VALUE_THRESH)
        self.assertEqual(
            alignments[5].title[:50],
            "gi|151946495|gb|EDN64717.1| Sir4p-interacting fact",
        )
        self.assertEqual(alignments[5].length, 535)
        self.assertEqual(len(alignments[5].hsps), 1)
        self.assertGreater(alignments[5].hsps[0].expect, E_VALUE_THRESH)
        self.assertEqual(
            alignments[6].title[:50],
            "gi|151943708|gb|EDN62018.1| nuclear pore complex s",
        )
        self.assertEqual(alignments[6].length, 349)
        self.assertEqual(len(alignments[6].hsps), 2)
        self.assertGreater(alignments[6].hsps[0].expect, E_VALUE_THRESH)
        self.assertGreater(alignments[6].hsps[1].expect, E_VALUE_THRESH)
        self.assertEqual(
            alignments[7].title[:50],
            "gi|151567866|pdb|2PM7|B Chain B, Crystal Structure",
        )
        self.assertEqual(alignments[7].length, 297)
        self.assertEqual(len(alignments[7].hsps), 2)
        self.assertGreater(alignments[7].hsps[0].expect, E_VALUE_THRESH)
        self.assertGreater(alignments[7].hsps[1].expect, E_VALUE_THRESH)
        self.assertEqual(
            alignments[8].title[:50],
            "gi|6321338|ref|NP_011415.1| Nuclear pore protein t",
        )
        self.assertEqual(alignments[8].length, 349)
        self.assertEqual(len(alignments[8].hsps), 2)
        self.assertGreater(alignments[8].hsps[0].expect, E_VALUE_THRESH)
        self.assertGreater(alignments[8].hsps[1].expect, E_VALUE_THRESH)
        self.assertEqual(
            alignments[9].title[:50],
            "gi|151567870|pdb|2PM9|B Chain B, Crystal Structure",
        )
        self.assertEqual(alignments[9].length, 297)
        self.assertEqual(len(alignments[9].hsps), 2)
        self.assertGreater(alignments[9].hsps[0].expect, E_VALUE_THRESH)
        self.assertGreater(alignments[9].hsps[1].expect, E_VALUE_THRESH)

        with open(datafile, "rb") as handle:
            record = NCBIXML.read(handle)
            handle.close()

    def test_xml_2218_blastp_002(self):
        """Parsing BLASTP 2.2.18+, SwissProt Q08386 and P07175, no hits (xml_2218_blastp_002)."""
        filename = "xml_2218_blastp_002.xml"
        datafile = os.path.join("Blast", filename)
        with open(datafile, "rb") as handle:
            records = NCBIXML.parse(handle)
            record = next(records)
            self.assertEqual(record.query_id, "gi|585505|sp|Q08386|MOPB_RHOCA")
            self.assertEqual(len(record.alignments), 0)
            record = next(records)
            self.assertEqual(record.query_id, "gi|129628|sp|P07175.1|PARA_AGRTU")
            self.assertEqual(len(record.alignments), 0)
            self.assertRaises(StopIteration, next, records)

    def test_xml_2218L_blastp_001(self):
        """Parsing BLASTP 2.2.18, Fake query (xml_2218L_blastp_001)."""
        filename = "xml_2218L_blastp_001.xml"
        datafile = os.path.join("Blast", filename)

        with open(datafile, "rb") as handle:
            records = NCBIXML.parse(handle)
            record = next(records)
            self.assertRaises(StopIteration, next, records)
        self.assertEqual(record.query_id, "lcl|1_0")
        alignments = record.alignments
        self.assertEqual(len(alignments), 0)

        with open(datafile, "rb") as handle:
            record = NCBIXML.read(handle)

    def test_xml_2222_blastx_001(self):
        """Parsing BLASTX 2.2.22+, multiple queries against NR (xml_2222_blastx_001)."""
        # See also plain text file bt081.txt (matching output from blastx tool)

        filename = "xml_2222_blastx_001.xml"
        datafile = os.path.join("Blast", filename)

        with open(datafile, "rb") as handle:
            records = NCBIXML.parse(handle)

            record = next(records)
            self.assertEqual(record.application, "BLASTX")
            self.assertEqual(record.version, "2.2.22+")
            self.assertEqual(record.date, "")
            self.assertEqual(
                record.query,
                "gi|4104054|gb|AH007193.1|SEG_CVIGS Centaurea vallesiaca 18S ribosomal RNA gene, partial sequence",
            )
            self.assertEqual(record.query_letters, 1002)
            self.assertEqual(record.database, "nr")
            self.assertEqual(record.num_sequences_in_database, 8994603)
            self.assertEqual(record.database_sequences, 8994603)
            # self.assertEqual(record.database_length, 3078807967)
            self.assertEqual(record.database_length, -1216159329)  # NCBI bug!
            self.assertEqual(len(record.descriptions), 1)
            self.assertEqual(len(record.alignments), 1)
            self.assertEqual(len(record.alignments[0].hsps), 1)

            record = next(records)
            self.assertEqual(record.application, "BLASTX")
            self.assertEqual(record.version, "2.2.22+")
            self.assertEqual(record.date, "")
            self.assertEqual(
                record.query,
                "gi|4218935|gb|AF074388.1|AF074388 Sambucus nigra hevein-like protein HLPf gene, partial cds",
            )
            self.assertEqual(record.query_letters, 2050)
            self.assertEqual(record.database, "nr")
            self.assertEqual(record.num_sequences_in_database, 8994603)
            self.assertEqual(record.database_sequences, 8994603)
            # self.assertEqual(record.database_length, 3078807967)
            self.assertEqual(record.database_length, -1216159329)  # NCBI bug!
            # I used -num_descriptions 10 and -num_alignments 1
            self.assertEqual(len(record.descriptions), 10)
            self.assertEqual(len(record.alignments), 10)
            self.assertEqual(len(record.alignments[0].hsps), 2)
            self.assertEqual(len(record.alignments[1].hsps), 2)
            self.assertEqual(len(record.alignments[9].hsps), 2)

            record = next(records)
            self.assertEqual(record.application, "BLASTX")
            self.assertEqual(record.version, "2.2.22+")
            self.assertEqual(record.date, "")
            self.assertEqual(
                record.query,
                "gi|5690369|gb|AF158246.1|AF158246 Cricetulus griseus glucose phosphate isomerase (GPI) gene, partial intron sequence",
            )
            self.assertEqual(record.query_letters, 550)
            self.assertEqual(record.database, "nr")
            self.assertEqual(record.num_sequences_in_database, 8994603)
            self.assertEqual(record.database_sequences, 8994603)
            # self.assertEqual(record.database_length, 3078807967)
            self.assertEqual(record.database_length, -1216159329)  # NCBI bug!
            self.assertEqual(len(record.descriptions), 0)
            self.assertEqual(len(record.alignments), 0)

            record = next(records)
            self.assertEqual(record.application, "BLASTX")
            self.assertEqual(record.version, "2.2.22+")
            self.assertEqual(record.date, "")
            self.assertEqual(
                record.query,
                "gi|5049839|gb|AI730987.1|AI730987 BNLGHi8354 Six-day Cotton fiber Gossypium hirsutum cDNA 5' similar to TUBULIN BETA-1 CHAIN gi|486734|pir|S35142 tubulin beta chain - white lupine gi|402636 (X70184) Beta tubulin 1 [Lupinus albus], mRNA sequence",
            )
            self.assertEqual(record.query_letters, 655)
            self.assertEqual(record.database, "nr")
            self.assertEqual(record.num_sequences_in_database, 8994603)
            self.assertEqual(record.database_sequences, 8994603)
            # self.assertEqual(record.database_length, 3078807967)
            self.assertEqual(record.database_length, -1216159329)  # NCBI bug!
            self.assertEqual(len(record.descriptions), 10)
            self.assertEqual(len(record.alignments), 10)
            self.assertEqual(len(record.alignments[0].hsps), 1)
            self.assertEqual(len(record.alignments[9].hsps), 1)

            record = next(records)
            self.assertEqual(record.application, "BLASTX")
            self.assertEqual(record.version, "2.2.22+")
            self.assertEqual(record.date, "")
            self.assertEqual(
                record.query,
                "gi|5052071|gb|AF067555.1|AF067555 Phlox stansburyi internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2, complete sequence",
            )
            self.assertEqual(record.query_letters, 623)
            self.assertEqual(record.database, "nr")
            self.assertEqual(record.num_sequences_in_database, 8994603)
            self.assertEqual(record.database_sequences, 8994603)
            # self.assertEqual(record.database_length, 3078807967)
            self.assertEqual(record.database_length, -1216159329)  # NCBI bug!
            self.assertEqual(len(record.descriptions), 10)
            self.assertEqual(len(record.alignments), 10)
            self.assertEqual(len(record.alignments[0].hsps), 2)
            self.assertEqual(len(record.alignments[9].hsps), 1)

            record = next(records)
            self.assertEqual(record.application, "BLASTX")
            self.assertEqual(record.version, "2.2.22+")
            self.assertEqual(record.date, "")
            self.assertEqual(
                record.query,
                "gi|3176602|gb|U78617.1|LOU78617 Lathyrus odoratus phytochrome A (PHYA) gene, partial cds",
            )
            self.assertEqual(record.query_letters, 309)
            self.assertEqual(record.database, "nr")
            self.assertEqual(record.num_sequences_in_database, 8994603)
            self.assertEqual(record.database_sequences, 8994603)
            # self.assertEqual(record.database_length, 3078807967)
            self.assertEqual(record.database_length, -1216159329)  # NCBI bug!
            self.assertEqual(len(record.descriptions), 10)
            self.assertEqual(len(record.alignments), 10)
            self.assertEqual(len(record.alignments[0].hsps), 1)
            self.assertEqual(len(record.alignments[9].hsps), 1)

            record = next(records)
            self.assertEqual(record.application, "BLASTX")
            self.assertEqual(record.version, "2.2.22+")
            self.assertEqual(record.date, "")
            self.assertEqual(
                record.query,
                "gi|5817701|gb|AF142731.1|AF142731 Wisteria frutescens maturase-like protein (matK) gene, complete cds; chloroplast gene for chloroplast product",
            )
            self.assertEqual(record.query_letters, 2551)
            self.assertEqual(record.database, "nr")
            self.assertEqual(record.num_sequences_in_database, 8994603)
            self.assertEqual(record.database_sequences, 8994603)
            # self.assertEqual(record.database_length, 3078807967)
            self.assertEqual(record.database_length, -1216159329)  # NCBI bug!
            self.assertEqual(len(record.descriptions), 10)
            self.assertEqual(len(record.alignments), 10)
            self.assertEqual(len(record.alignments[0].hsps), 1)
            self.assertEqual(len(record.alignments[9].hsps), 1)

            self.assertRaises(StopIteration, next, records)

    def test_xml_2222_blastp_001(self):
        """Parsing BLASTP 2.2.22+, multiple queries against NR (xml_2222_blastp_001)."""
        # This is from blastp NOT blastall

        filename = "xml_2222_blastp_001.xml"
        datafile = os.path.join("Blast", filename)

        with open(datafile, "rb") as handle:
            records = NCBIXML.parse(handle)

            record = next(records)
            self.assertEqual(record.application, "BLASTP")
            self.assertEqual(record.version, "2.2.22+")
            self.assertEqual(record.date, "")
            self.assertEqual(record.query, "gi|3298468|dbj|BAA31520.1| SAMIPF")
            self.assertEqual(record.query_letters, 107)
            self.assertEqual(record.database, "nr")
            self.assertEqual(record.num_sequences_in_database, 8994603)
            self.assertEqual(record.database_sequences, 8994603)
            # self.assertEqual(record.database_length, 3078807967)
            self.assertEqual(record.database_length, -1216159329)  # NCBI bug!
            self.assertEqual(len(record.descriptions), 10)
            self.assertEqual(len(record.alignments), 10)
            self.assertEqual(len(record.alignments[0].hsps), 1)

            record = next(records)
            self.assertEqual(
                record.query,
                "gi|2781234|pdb|1JLY|B Chain B, Crystal Structure Of Amaranthus Caudatus Agglutinin",
            )
            self.assertEqual(record.query_letters, 304)

            record = next(records)
            self.assertEqual(
                record.query,
                "gi|4959044|gb|AAD34209.1|AF069992_1 LIM domain interacting RING finger protein",
            )
            self.assertEqual(record.query_letters, 600)

            record = next(records)
            self.assertEqual(
                record.query, "gi|671626|emb|CAA85685.1| rubisco large subunit"
            )
            self.assertEqual(record.query_letters, 473)

            self.assertRaises(StopIteration, next, records)

    def test_xml_2218L_rpsblast_001(self):
        """Parsing PSI-BLASTP 2.2.18, single query which converges in 3 iterations (xml_2218L_rpsblast_001)."""
        # This is from old pgpblast command line tool, NOT new psiblast
        # NOTE - The parser currently returns three BLAST record objects.
        # The old text parser would return a single PSI BLAST record object with three rounds.
        # This may change... although it may require a PSI BLAST specific XML parser.

        filename = "xml_2218L_rpsblast_001.xml"
        datafile = os.path.join("Blast", filename)

        with open(datafile, "rb") as handle:
            records = NCBIXML.parse(handle)

            record = next(records)
            self.assertEqual(record.application, "BLASTP")
            self.assertEqual(record.version, "2.2.18")
            self.assertEqual(record.date, "Mar-02-2008")
            self.assertEqual(record.query, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
            self.assertEqual(record.query_letters, 131)
            self.assertEqual(record.database, "/opt/BlastDBs/nr")
            self.assertEqual(record.num_sequences_in_database, 2563094)
            self.assertEqual(record.database_sequences, 2563094)
            self.assertEqual(record.database_length, 864488805)
            self.assertEqual(len(record.descriptions), 11)
            self.assertEqual(len(record.alignments), 11)
            self.assertEqual(len(record.alignments[0].hsps), 1)
            hsp = record.alignments[0].hsps[0]
            self.assertEqual(hsp.align_length, 131)
            self.assertEqual(hsp.identities, 131)
            self.assertEqual(hsp.positives, 131)
            self.assertEqual(
                hsp.query,
                "MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQPIPGEVIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSKGAKKLAQAMAS",
            )
            self.assertEqual(
                hsp.sbjct,
                "MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQPIPGEVIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSKGAKKLAQAMAS",
            )
            self.assertEqual(
                hsp.match,
                "MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQPIPGEVIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSKGAKKLAQAMAS",
            )
            self.assertEqual(hsp.score, 680)
            self.assertEqual(hsp.expect, 4.72196e-70)
            self.assertEqual(hsp.query_start, 1)
            self.assertEqual(hsp.query_end, 131)
            self.assertEqual(hsp.sbjct_start, 1)
            self.assertEqual(hsp.sbjct_end, 131)
            self.assertEqual(len(record.alignments[1].hsps), 1)
            hsp = record.alignments[1].hsps[0]
            self.assertEqual(hsp.align_length, 77)
            self.assertEqual(hsp.identities, 36)
            self.assertEqual(hsp.positives, 49)
            self.assertEqual(
                hsp.query,
                "MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQ",
            )
            self.assertEqual(
                hsp.sbjct,
                "MAREEPYKGDYVGGVAKILQGYFANYYGFPNVSLRLAGEEANLSKTGHANAKAIVHEMIKVIKEASKPLR-RGKGFK",
            )
            self.assertEqual(
                hsp.match,
                "MA+ EP KGDY GG  KIL  +     G+P V+L+LAGEEAN  + G    K  +H ++K+I +A KP R +G GF+",
            )
            self.assertEqual(hsp.score, 181)
            self.assertEqual(hsp.expect, 3.03476e-12)
            self.assertEqual(hsp.query_start, 1)
            self.assertEqual(hsp.query_end, 77)
            self.assertEqual(hsp.sbjct_start, 1)
            self.assertEqual(hsp.sbjct_end, 76)

            record = next(records)
            self.assertEqual(record.application, "BLASTP")
            self.assertEqual(record.version, "2.2.18")
            self.assertEqual(record.date, "Mar-02-2008")
            self.assertEqual(record.query, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
            self.assertEqual(record.query_letters, 131)
            self.assertEqual(record.database, "/opt/BlastDBs/nr")
            self.assertEqual(record.num_sequences_in_database, 2563094)
            self.assertEqual(record.database_sequences, 2563094)
            self.assertEqual(record.database_length, 864488805)
            self.assertEqual(len(record.descriptions), 19)
            self.assertEqual(len(record.alignments), 19)
            self.assertEqual(len(record.alignments[0].hsps), 1)
            hsp = record.alignments[0].hsps[0]
            self.assertEqual(hsp.align_length, 131)
            self.assertEqual(hsp.identities, 131)
            self.assertEqual(hsp.positives, 131)
            self.assertEqual(
                hsp.query,
                "MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQPIPGEVIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSKGAKKLAQAMAS",
            )
            self.assertEqual(
                hsp.sbjct,
                "MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQPIPGEVIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSKGAKKLAQAMAS",
            )
            self.assertEqual(
                hsp.match,
                "MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQPIPGEVIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSKGAKKLAQAMAS",
            )
            self.assertEqual(hsp.score, 590)
            self.assertEqual(hsp.expect, 1.28615e-59)
            self.assertEqual(hsp.query_start, 1)
            self.assertEqual(hsp.query_end, 131)
            self.assertEqual(hsp.sbjct_start, 1)
            self.assertEqual(hsp.sbjct_end, 131)

            record = next(records)
            self.assertEqual(record.application, "BLASTP")
            self.assertEqual(record.version, "2.2.18")
            self.assertEqual(record.date, "Mar-02-2008")
            self.assertEqual(record.query, "tr|Q3V4Q3|Q3V4Q3_9VIRU")
            self.assertEqual(record.query_letters, 131)
            self.assertEqual(record.database, "/opt/BlastDBs/nr")
            self.assertEqual(record.num_sequences_in_database, 2563094)
            self.assertEqual(record.database_sequences, 2563094)
            self.assertEqual(record.database_length, 864488805)
            self.assertEqual(len(record.descriptions), 9)
            self.assertEqual(len(record.alignments), 9)
            self.assertEqual(len(record.alignments[0].hsps), 1)
            hsp = record.alignments[0].hsps[0]
            self.assertEqual(hsp.align_length, 131)
            self.assertEqual(hsp.identities, 131)
            self.assertEqual(hsp.positives, 131)
            self.assertEqual(
                hsp.query,
                "MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQPIPGEVIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSKGAKKLAQAMAS",
            )
            self.assertEqual(
                hsp.sbjct,
                "MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQPIPGEVIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSKGAKKLAQAMAS",
            )
            self.assertEqual(
                hsp.match,
                "MAKYEPKKGDYAGGAVKILDMFENGQLGYPEVTLKLAGEEANARRAGDERTKEAIHAIVKMISDAMKPYRNKGSGFQSQPIPGEVIAQVTSNPEYQQAKAFLASPATQVRNIEREEVLSKGAKKLAQAMAS",
            )
            self.assertEqual(hsp.score, 535)
            self.assertEqual(hsp.expect, 3.43623e-53)
            self.assertEqual(hsp.query_start, 1)
            self.assertEqual(hsp.query_end, 131)
            self.assertEqual(hsp.sbjct_start, 1)
            self.assertEqual(hsp.sbjct_end, 131)

            # TODO - Can we detect the convergence status:
            # <Iteration_message>CONVERGED</Iteration_message>
            self.assertRaises(StopIteration, next, records)

    def test_xml_2900_blastp_001_v1(self):
        record = self._test_xml_2900_blastp_001("xml_2900_blastp_001.xml")

        description = record.descriptions[0]
        self.assertEqual(len(description.title), 4706)
        self.assertEqual(
            description.title[:300],
            "gi|447157535|ref|WP_001234791.1| MULTISPECIES: Sec-independent protein translocase subunit TatA [Shigella] >gi|24115132|ref|NP_709642.1| twin-arginine translocation protein TatA [Shigella flexneri 2a str. 301] >gi|82778983|ref|YP_405332.1| twin-arginine translocation protein TatA [Shigella dysenteri",
        )
        description = record.descriptions[1]
        self.assertEqual(len(description.title), 106)
        self.assertEqual(
            description.title,
            "gi|91074959|gb|ABE09840.1| sec-independent twin-arginine translocase subunit tatA [Escherichia coli UTI89]",
        )
        description = record.descriptions[2]
        self.assertEqual(len(description.title), 979)
        self.assertEqual(
            description.title[:300],
            "gi|73857826|gb|AAZ90533.1| conserved hypothetical protein [Shigella sonnei Ss046] >gi|331067587|gb|EGI38991.1| Sec-independent protein translocase protein TatA [Escherichia coli TA280] >gi|412965214|emb|CCK49144.1| sec-independent protein translocase protein tata/e homolog 2 [Escherichia coli chi712",
        )
        description = record.descriptions[3]
        self.assertEqual(len(description.title), 260)
        self.assertEqual(
            description.title,
            "gi|25302684|pir||D86071 hypothetical protein tatA [imported] - Escherichia coli (strain O157:H7, substrain EDL933) >gi|12518713|gb|AAG59032.1|AE005614_12 twin arginine translocation protein; sec-independent protein export [Escherichia coli O157:H7 str. EDL933]",
        )
        description = record.descriptions[4]
        self.assertEqual(len(description.title), 100)
        self.assertEqual(
            description.title,
            "gi|331072495|gb|EGI43827.1| Sec-independent protein translocase protein TatA [Escherichia coli H591]",
        )
        description = record.descriptions[5]
        self.assertEqual(len(description.title), 120)
        self.assertEqual(
            description.title,
            "gi|808042844|pdb|2MN7|A Chain A, Solution structure of monomeric TatA of twin-arginine translocation system from E. coli",
        )
        description = record.descriptions[6]
        self.assertEqual(len(description.title), 774516)
        self.assertEqual(
            description.title[:300],
            "gi|481023661|ref|WP_001295260.1| MULTISPECIES: Sec-independent protein translocase subunit TatA [Proteobacteria] >gi|15834020|ref|NP_312793.1| TatABCE protein translocation system subunit TatA [Escherichia coli O157:H7 str. Sakai] >gi|90111653|ref|NP_418280.4| twin arginine protein translocation sys",
        )
        description = record.descriptions[7]
        self.assertEqual(len(description.title), 238)
        self.assertEqual(
            description.title,
            "gi|808042842|pdb|2MN6|B Chain B, Solution structure of dimeric TatA of twin-arginine translocation system from E. coli >gi|808042843|pdb|2MN6|A Chain A, Solution structure of dimeric TatA of twin-arginine translocation system from E. coli",
        )
        description = record.descriptions[8]
        self.assertEqual(len(description.title), 2111)
        self.assertEqual(
            description.title[:300],
            "gi|491167042|ref|WP_005025412.1| MULTISPECIES: twin-arginine translocase subunit TatA [Enterobacteriaceae] >gi|320176770|gb|EFW51804.1| Twin-arginine translocation protein TatA [Shigella dysenteriae CDC 74-1112] >gi|391297993|gb|EIQ56018.1| twin arginine-targeting translocase, TatA/E family protein ",
        )
        description = record.descriptions[9]
        self.assertEqual(len(description.title), 384)
        self.assertEqual(
            description.title[:300],
            "gi|332996954|gb|EGK16572.1| sec-independent translocase protein tatA [Shigella flexneri VA-6] >gi|391245379|gb|EIQ04650.1| twin arginine-targeting translocase, TatA/E family protein [Shigella flexneri K-1770] >gi|1411457050|emb|SRN34259.1| twin arginine translocase protein A [Shigella flexneri] >gi|",
        )

    def test_xml_2900_blastp_001_v2(self):
        record = self._test_xml_2900_blastp_001("xml_2900_blastp_001_v2.xml")

        alignment = record.alignments[0]
        self.assertEqual(
            alignment.title,
            "gi|447157535|ref|WP_001234791.1| MULTISPECIES: Sec-independent protein translocase subunit TatA [Shigella]",
        )

        description = record.descriptions[0]
        self.assertEqual(
            description.title,
            "gi|447157535|ref|WP_001234791.1| MULTISPECIES: Sec-independent protein translocase subunit TatA [Shigella]",
        )
        self.assertEqual(len(description.items), 48)
        description_item = description.items[0]
        self.assertEqual(description_item.id, "gi|447157535|ref|WP_001234791.1|")
        self.assertEqual(description_item.accession, "WP_001234791")
        self.assertEqual(
            description_item.title,
            "MULTISPECIES: Sec-independent protein translocase subunit TatA [Shigella]",
        )
        self.assertEqual(description_item.taxid, 620)
        self.assertEqual(description_item.sciname, "Shigella")

    def _test_xml_2900_blastp_001(self, filename):
        datafile = os.path.join("Blast", filename)
        with open(datafile, "rb") as handle:
            records = NCBIXML.parse(handle)
            record = next(records)
            self.assertRaises(StopIteration, next, records)

        self.assertEqual(record.application, "BLASTP")
        self.assertEqual(record.version, "2.9.0+")
        self.assertEqual(
            record.query, "twin argininte translocase protein A [Escherichia coli K12]"
        )
        self.assertEqual(record.query_letters, 103)
        self.assertEqual(record.database, "nr")
        self.assertEqual(record.num_sequences_in_database, 194611632)
        self.assertEqual(record.database_sequences, 194611632)
        self.assertEqual(record.database_length, 2104817704)

        alignment = record.alignments[0]
        self.assertEqual(alignment.hit_id, "gi|447157535|ref|WP_001234791.1|")
        self.assertEqual(alignment.accession, "WP_001234791")
        self.assertEqual(alignment.length, 103)
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertEqual(hsp.align_length, 103)
        self.assertEqual(hsp.identities, 103)
        self.assertEqual(hsp.positives, 103)
        self.assertEqual(
            hsp.query,
            "MRLCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV",
        )
        self.assertEqual(
            hsp.sbjct,
            "MRLCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV",
        )
        self.assertEqual(
            hsp.match,
            "MRLCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV",
        )
        self.assertEqual(hsp.score, 532.0)
        self.assertEqual(hsp.expect, 1.71849e-68)
        self.assertEqual(hsp.query_start, 1)
        self.assertEqual(hsp.query_end, 103)
        self.assertEqual(hsp.sbjct_start, 1)
        self.assertEqual(hsp.sbjct_end, 103)

        alignment = record.alignments[1]
        self.assertEqual(alignment.hit_id, "gi|91074959|gb|ABE09840.1|")
        self.assertEqual(alignment.accession, "ABE09840")
        self.assertEqual(
            alignment.title,
            "gi|91074959|gb|ABE09840.1| sec-independent twin-arginine translocase subunit tatA [Escherichia coli UTI89]",
        )
        self.assertEqual(alignment.length, 103)
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertEqual(hsp.align_length, 103)
        self.assertEqual(hsp.identities, 102)
        self.assertEqual(hsp.positives, 102)
        self.assertEqual(
            hsp.query,
            "MRLCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV",
        )
        self.assertEqual(
            hsp.sbjct,
            "MRLCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKIEDAKRHDKEQV",
        )
        self.assertEqual(
            hsp.match,
            "MRLCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAK EDAKRHDKEQV",
        )
        self.assertEqual(hsp.score, 526.0)
        self.assertEqual(hsp.expect, 1.44614e-67)
        self.assertEqual(hsp.query_start, 1)
        self.assertEqual(hsp.query_end, 103)
        self.assertEqual(hsp.sbjct_start, 1)
        self.assertEqual(hsp.sbjct_end, 103)

        alignment = record.alignments[2]
        self.assertEqual(alignment.hit_id, "gi|73857826|gb|AAZ90533.1|")
        self.assertEqual(alignment.accession, "AAZ90533")
        self.assertEqual(alignment.length, 103)
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertEqual(hsp.align_length, 103)
        self.assertEqual(hsp.identities, 102)
        self.assertEqual(hsp.positives, 102)
        self.assertEqual(
            hsp.query,
            "MRLCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV",
        )
        self.assertEqual(
            hsp.sbjct,
            "MRPCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV",
        )
        self.assertEqual(
            hsp.match,
            "MR CLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV",
        )
        self.assertEqual(hsp.score, 525.0)
        self.assertEqual(hsp.expect, 1.80126e-67)
        self.assertEqual(hsp.query_start, 1)
        self.assertEqual(hsp.query_end, 103)
        self.assertEqual(hsp.sbjct_start, 1)
        self.assertEqual(hsp.sbjct_end, 103)

        alignment = record.alignments[3]
        self.assertEqual(alignment.hit_id, "gi|25302684|pir||D86071")
        self.assertEqual(alignment.accession, "D86071")
        self.assertEqual(alignment.length, 103)
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertEqual(hsp.align_length, 103)
        self.assertEqual(hsp.identities, 102)
        self.assertEqual(hsp.positives, 102)
        self.assertEqual(
            hsp.query,
            "MRLCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV",
        )
        self.assertEqual(
            hsp.sbjct,
            "MRLCLIIIYHRXTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV",
        )
        self.assertEqual(
            hsp.match,
            "MRLCLIIIYHR TCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV",
        )
        self.assertEqual(hsp.score, 525.0)
        self.assertEqual(hsp.expect, 2.10054e-67)
        self.assertEqual(hsp.query_start, 1)
        self.assertEqual(hsp.query_end, 103)
        self.assertEqual(hsp.sbjct_start, 1)
        self.assertEqual(hsp.sbjct_end, 103)

        alignment = record.alignments[4]
        self.assertEqual(alignment.hit_id, "gi|331072495|gb|EGI43827.1|")
        self.assertEqual(alignment.accession, "EGI43827")
        self.assertEqual(alignment.length, 103)
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertEqual(hsp.align_length, 103)
        self.assertEqual(hsp.identities, 101)
        self.assertEqual(hsp.positives, 102)
        self.assertEqual(
            hsp.query,
            "MRLCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV",
        )
        self.assertEqual(
            hsp.sbjct,
            "MRPCLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQANTNQEQAKTEDAKRHDKEQV",
        )
        self.assertEqual(
            hsp.match,
            "MR CLIIIYHRGTCMGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQA+TNQEQAKTEDAKRHDKEQV",
        )
        self.assertEqual(hsp.score, 522.0)
        self.assertEqual(hsp.expect, 5.11164e-67)
        self.assertEqual(hsp.query_start, 1)
        self.assertEqual(hsp.query_end, 103)
        self.assertEqual(hsp.sbjct_start, 1)
        self.assertEqual(hsp.sbjct_end, 103)

        alignment = record.alignments[5]
        self.assertEqual(alignment.hit_id, "gi|808042844|pdb|2MN7|A")
        self.assertEqual(alignment.accession, "2MN7_A")
        self.assertEqual(alignment.length, 97)
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertEqual(hsp.align_length, 89)
        self.assertEqual(hsp.identities, 89)
        self.assertEqual(hsp.positives, 89)
        self.assertEqual(
            hsp.query,
            "MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV",
        )
        self.assertEqual(
            hsp.sbjct,
            "MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV",
        )
        self.assertEqual(
            hsp.match,
            "MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV",
        )
        self.assertEqual(hsp.score, 449.0)
        self.assertEqual(hsp.expect, 5.99855e-56)
        self.assertEqual(hsp.query_start, 15)
        self.assertEqual(hsp.query_end, 103)
        self.assertEqual(hsp.sbjct_start, 1)
        self.assertEqual(hsp.sbjct_end, 89)

        alignment = record.alignments[6]
        self.assertEqual(alignment.hit_id, "gi|481023661|ref|WP_001295260.1|")
        self.assertEqual(alignment.accession, "WP_001295260")
        self.assertEqual(alignment.length, 89)
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertEqual(hsp.align_length, 89)
        self.assertEqual(hsp.identities, 89)
        self.assertEqual(hsp.positives, 89)
        self.assertEqual(
            hsp.query,
            "MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV",
        )
        self.assertEqual(
            hsp.sbjct,
            "MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV",
        )
        self.assertEqual(
            hsp.match,
            "MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV",
        )
        self.assertEqual(hsp.score, 448.0)
        self.assertEqual(hsp.expect, 8.86198e-56)
        self.assertEqual(hsp.query_start, 15)
        self.assertEqual(hsp.query_end, 103)
        self.assertEqual(hsp.sbjct_start, 1)
        self.assertEqual(hsp.sbjct_end, 89)

        alignment = record.alignments[7]
        self.assertEqual(alignment.hit_id, "gi|808042842|pdb|2MN6|B")
        self.assertEqual(alignment.accession, "2MN6_B")
        self.assertEqual(alignment.length, 100)
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertEqual(hsp.align_length, 89)
        self.assertEqual(hsp.identities, 89)
        self.assertEqual(hsp.positives, 89)
        self.assertEqual(
            hsp.query,
            "MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV",
        )
        self.assertEqual(
            hsp.sbjct,
            "MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV",
        )
        self.assertEqual(
            hsp.match,
            "MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV",
        )
        self.assertEqual(hsp.score, 447.0)
        self.assertEqual(hsp.expect, 1.55019e-55)
        self.assertEqual(hsp.query_start, 15)
        self.assertEqual(hsp.query_end, 103)
        self.assertEqual(hsp.sbjct_start, 4)
        self.assertEqual(hsp.sbjct_end, 92)

        alignment = record.alignments[8]
        self.assertEqual(alignment.hit_id, "gi|491167042|ref|WP_005025412.1|")
        self.assertEqual(alignment.accession, "WP_005025412")
        self.assertEqual(alignment.length, 89)
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertEqual(hsp.align_length, 89)
        self.assertEqual(hsp.identities, 88)
        self.assertEqual(hsp.positives, 89)
        self.assertEqual(
            hsp.query,
            "MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV",
        )
        self.assertEqual(
            hsp.sbjct,
            "MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDSDFTAKTIADKQADTNQEQAKTEDAKRHDKEQV",
        )
        self.assertEqual(
            hsp.match,
            "MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQD+DFTAKTIADKQADTNQEQAKTEDAKRHDKEQV",
        )
        self.assertEqual(hsp.score, 445.0)
        self.assertEqual(hsp.expect, 2.25306e-55)
        self.assertEqual(hsp.query_start, 15)
        self.assertEqual(hsp.query_end, 103)
        self.assertEqual(hsp.sbjct_start, 1)
        self.assertEqual(hsp.sbjct_end, 89)

        alignment = record.alignments[9]
        self.assertEqual(alignment.hit_id, "gi|332996954|gb|EGK16572.1|")
        self.assertEqual(alignment.accession, "EGK16572")
        self.assertEqual(alignment.length, 89)
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertEqual(hsp.align_length, 89)
        self.assertEqual(hsp.identities, 88)
        self.assertEqual(hsp.positives, 89)
        self.assertEqual(
            hsp.query,
            "MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTIADKQADTNQEQAKTEDAKRHDKEQV",
        )
        self.assertEqual(
            hsp.sbjct,
            "MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTISDKQADTNQEQAKTEDAKRHDKEQV",
        )
        self.assertEqual(
            hsp.match,
            "MGGISIWQLLIIAVIVVLLFGTKKLGSIGSDLGASIKGFKKAMSDDEPKQDKTSQDADFTAKTI+DKQADTNQEQAKTEDAKRHDKEQV",
        )
        self.assertEqual(hsp.score, 445.0)
        self.assertEqual(hsp.expect, 2.25306e-55)
        self.assertEqual(hsp.query_start, 15)
        self.assertEqual(hsp.query_end, 103)
        self.assertEqual(hsp.sbjct_start, 1)
        self.assertEqual(hsp.sbjct_end, 89)

        self.assertEqual(len(record.descriptions), 10)

        description = record.descriptions[0]
        self.assertEqual(description.score, 532.0)
        self.assertEqual(description.e, 1.71849e-68)
        self.assertEqual(description.num_alignments, 1)
        description = record.descriptions[1]
        self.assertEqual(description.score, 526.0)
        self.assertEqual(description.e, 1.44614e-67)
        self.assertEqual(description.num_alignments, 1)
        description = record.descriptions[2]
        self.assertEqual(description.score, 525.0)
        self.assertEqual(description.e, 1.80126e-67)
        self.assertEqual(description.num_alignments, 1)
        description = record.descriptions[3]
        self.assertEqual(description.score, 525.0)
        self.assertEqual(description.e, 2.10054e-67)
        self.assertEqual(description.num_alignments, 1)
        description = record.descriptions[4]
        self.assertEqual(description.score, 522.0)
        self.assertEqual(description.e, 5.11164e-67)
        self.assertEqual(description.num_alignments, 1)
        description = record.descriptions[5]
        self.assertEqual(description.score, 449.0)
        self.assertEqual(description.e, 5.99855e-56)
        self.assertEqual(description.num_alignments, 1)
        description = record.descriptions[6]
        self.assertEqual(description.score, 448.0)
        self.assertEqual(description.e, 8.86198e-56)
        self.assertEqual(description.num_alignments, 1)
        description = record.descriptions[7]
        self.assertEqual(description.score, 447.0)
        self.assertEqual(description.e, 1.55019e-55)
        self.assertEqual(description.num_alignments, 1)
        description = record.descriptions[8]
        self.assertEqual(description.score, 445.0)
        self.assertEqual(description.e, 2.25306e-55)
        self.assertEqual(description.num_alignments, 1)
        description = record.descriptions[9]
        self.assertEqual(description.score, 445.0)
        self.assertEqual(description.e, 2.25306e-55)
        self.assertEqual(description.num_alignments, 1)

        return record

    def test_xml_2900_blastn_001_v1(self):
        self._test_xml_2900_blastn_001("xml_2900_blastn_001.xml")

    def test_xml_2900_blastn_001_v2(self):
        self._test_xml_2900_blastn_001("xml_2900_blastn_001_v2.xml")

    def _test_xml_2900_blastn_001(self, filename):
        datafile = os.path.join("Blast", filename)
        with open(datafile, "rb") as handle:
            records = NCBIXML.parse(handle)
            record = next(records)
            self.assertRaises(StopIteration, next, records)

        self.assertEqual(record.application, "BLASTN")
        self.assertEqual(record.version, "2.9.0+")
        self.assertEqual(record.query, "human STS STS_D11570, sequence tagged site")
        self.assertEqual(record.query_letters, 285)
        self.assertEqual(
            record.database, "GPIPE/10090/current/all_top_level GPIPE/10090/current/rna"
        )
        self.assertEqual(record.num_sequences_in_database, 107382)
        self.assertEqual(record.database_sequences, 107382)
        self.assertEqual(record.database_length, 3164670549)

        alignment = record.alignments[0]
        self.assertEqual(alignment.hit_id, "gi|372099107|ref|NC_000069.6|")
        self.assertEqual(alignment.accession, "NC_000069")
        self.assertEqual(alignment.length, 160039680)
        self.assertEqual(
            alignment.title,
            "gi|372099107|ref|NC_000069.6| Mus musculus strain C57BL/6J chromosome 3, GRCm38.p4 C57BL/6J",
        )
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertEqual(hsp.align_length, 34)
        self.assertEqual(hsp.identities, 30)
        self.assertEqual(hsp.positives, 30)
        self.assertEqual(hsp.query, "GAATCCTAGAGGCTTGATTGGCCCAGG-CTGCTG")
        self.assertEqual(hsp.sbjct, "GAATCCTAGAGGCTGGACTGGCCCTGGCCTGCTG")
        self.assertEqual(hsp.match, "|||||||||||||| || |||||| || ||||||")
        self.assertEqual(hsp.score, 44.0)
        self.assertEqual(hsp.expect, 0.375311)
        self.assertEqual(hsp.query_start, 134)
        self.assertEqual(hsp.query_end, 166)
        self.assertEqual(hsp.sbjct_start, 101449177)
        self.assertEqual(hsp.sbjct_end, 101449144)
        self.assertEqual(hsp.frame, (1, -1))
        self.assertEqual(hsp.strand, ("Plus", "Minus"))

        alignment = record.alignments[1]
        self.assertEqual(alignment.hit_id, "gi|372099103|ref|NC_000073.6|")
        self.assertEqual(alignment.accession, "NC_000073")
        self.assertEqual(alignment.length, 145441459)
        self.assertEqual(
            alignment.title,
            "gi|372099103|ref|NC_000073.6| Mus musculus strain C57BL/6J chromosome 7, GRCm38.p4 C57BL/6J",
        )
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertEqual(hsp.align_length, 29)
        self.assertEqual(hsp.identities, 26)
        self.assertEqual(hsp.positives, 26)
        self.assertEqual(hsp.query, "GAAAGGAAATNAAAATGGAAAGTTCTTGT")
        self.assertEqual(hsp.sbjct, "GAAAGGAAAAAAAAATGGAAAGTTCTGGT")
        self.assertEqual(hsp.match, "|||||||||  ||||||||||||||| ||")
        self.assertEqual(hsp.score, 44.0)
        self.assertEqual(hsp.expect, 0.375311)
        self.assertEqual(hsp.query_start, 205)
        self.assertEqual(hsp.query_end, 233)
        self.assertEqual(hsp.sbjct_start, 131772185)
        self.assertEqual(hsp.sbjct_end, 131772157)
        self.assertEqual(hsp.frame, (1, -1))
        self.assertEqual(hsp.strand, ("Plus", "Minus"))

        alignment = record.alignments[2]
        self.assertEqual(alignment.hit_id, "gi|372099106|ref|NC_000070.6|")
        self.assertEqual(alignment.accession, "NC_000070")
        self.assertEqual(alignment.length, 156508116)
        self.assertEqual(
            alignment.title,
            "gi|372099106|ref|NC_000070.6| Mus musculus strain C57BL/6J chromosome 4, GRCm38.p4 C57BL/6J",
        )
        self.assertEqual(len(alignment.hsps), 2)
        hsp = alignment.hsps[0]
        self.assertEqual(hsp.align_length, 24)
        self.assertEqual(hsp.identities, 23)
        self.assertEqual(hsp.positives, 23)
        self.assertEqual(hsp.query, "CCAACACAGGCCAGCGACTTCTGG")
        self.assertEqual(hsp.sbjct, "CCAACACAGGCCAGCGGCTTCTGG")
        self.assertEqual(hsp.match, "|||||||||||||||| |||||||")
        self.assertEqual(hsp.score, 43.0)
        self.assertEqual(hsp.expect, 1.30996)
        self.assertEqual(hsp.query_start, 62)
        self.assertEqual(hsp.query_end, 85)
        self.assertEqual(hsp.sbjct_start, 9607562)
        self.assertEqual(hsp.sbjct_end, 9607539)
        self.assertEqual(hsp.frame, (1, -1))
        self.assertEqual(hsp.strand, ("Plus", "Minus"))
        hsp = alignment.hsps[1]
        self.assertEqual(hsp.align_length, 32)
        self.assertEqual(hsp.identities, 28)
        self.assertEqual(hsp.positives, 28)
        self.assertEqual(hsp.query, "GCCTGACATGG-GTAGCTGCTCAATAAATGCT")
        self.assertEqual(hsp.sbjct, "GCCTGGCATGAAGTAACTGCTCAATAAATGCT")
        self.assertEqual(hsp.match, "||||| ||||  ||| ||||||||||||||||")
        self.assertEqual(hsp.score, 40.0)
        self.assertEqual(hsp.expect, 4.57222)
        self.assertEqual(hsp.query_start, 242)
        self.assertEqual(hsp.query_end, 272)
        self.assertEqual(hsp.sbjct_start, 142902532)
        self.assertEqual(hsp.sbjct_end, 142902563)
        self.assertEqual(hsp.frame, (1, 1))
        self.assertEqual(hsp.strand, ("Plus", "Plus"))

        alignment = record.alignments[3]
        self.assertEqual(alignment.hit_id, "gi|372099108|ref|NC_000068.7|")
        self.assertEqual(alignment.accession, "NC_000068")
        self.assertEqual(alignment.length, 182113224)
        self.assertEqual(
            alignment.title,
            "gi|372099108|ref|NC_000068.7| Mus musculus strain C57BL/6J chromosome 2, GRCm38.p4 C57BL/6J",
        )
        self.assertEqual(len(alignment.hsps), 2)
        hsp = alignment.hsps[0]
        self.assertEqual(hsp.align_length, 31)
        self.assertEqual(hsp.identities, 27)
        self.assertEqual(hsp.positives, 27)
        self.assertEqual(hsp.query, "AAGGCCTGACATGGGTAGCTGCTCAATAAAT")
        self.assertEqual(hsp.sbjct, "AAGTCCTGGCATGAGTAGTTGCTCAATAAAT")
        self.assertEqual(hsp.match, "||| |||| |||| |||| ||||||||||||")
        self.assertEqual(hsp.score, 42.0)
        self.assertEqual(hsp.expect, 1.30996)
        self.assertEqual(hsp.query_start, 239)
        self.assertEqual(hsp.query_end, 269)
        self.assertEqual(hsp.sbjct_start, 3799647)
        self.assertEqual(hsp.sbjct_end, 3799677)
        self.assertEqual(hsp.frame, (1, 1))
        self.assertEqual(hsp.strand, ("Plus", "Plus"))
        hsp = alignment.hsps[1]
        self.assertEqual(hsp.align_length, 25)
        self.assertEqual(hsp.identities, 23)
        self.assertEqual(hsp.positives, 23)
        self.assertEqual(hsp.query, "AAATNAAAATGGAAAGTTCTTGTAG")
        self.assertEqual(hsp.sbjct, "AAATGAAAATGGAAAGTTCTTATAG")
        self.assertEqual(hsp.match, "|||| |||||||||||||||| |||")
        self.assertEqual(hsp.score, 41.0)
        self.assertEqual(hsp.expect, 4.57222)
        self.assertEqual(hsp.query_start, 211)
        self.assertEqual(hsp.query_end, 235)
        self.assertEqual(hsp.sbjct_start, 70278960)
        self.assertEqual(hsp.sbjct_end, 70278984)
        self.assertEqual(hsp.frame, (1, 1))
        self.assertEqual(hsp.strand, ("Plus", "Plus"))

        alignment = record.alignments[4]
        self.assertEqual(alignment.hit_id, "gi|372099097|ref|NC_000079.6|")
        self.assertEqual(alignment.accession, "NC_000079")
        self.assertEqual(alignment.length, 120421639)
        self.assertEqual(
            alignment.title,
            "gi|372099097|ref|NC_000079.6| Mus musculus strain C57BL/6J chromosome 13, GRCm38.p4 C57BL/6J",
        )
        self.assertEqual(len(alignment.hsps), 2)
        hsp = alignment.hsps[0]
        self.assertEqual(hsp.align_length, 28)
        self.assertEqual(hsp.identities, 25)
        self.assertEqual(hsp.positives, 25)
        self.assertEqual(hsp.query, "AAGGAAATNAAAATGGAAAGTTCTTGTA")
        self.assertEqual(hsp.sbjct, "AAGGACATCAAAATGGAAAGTTCTTCTA")
        self.assertEqual(hsp.match, "||||| || |||||||||||||||| ||")
        self.assertEqual(hsp.score, 42.0)
        self.assertEqual(hsp.expect, 1.30996)
        self.assertEqual(hsp.query_start, 207)
        self.assertEqual(hsp.query_end, 234)
        self.assertEqual(hsp.sbjct_start, 26806584)
        self.assertEqual(hsp.sbjct_end, 26806557)
        self.assertEqual(hsp.frame, (1, -1))
        self.assertEqual(hsp.strand, ("Plus", "Minus"))
        hsp = alignment.hsps[1]
        self.assertEqual(hsp.align_length, 40)
        self.assertEqual(hsp.identities, 32)
        self.assertEqual(hsp.positives, 32)
        self.assertEqual(hsp.query, "AGCGCAAGGCCTGACATGGGTAGCTGCTCAATAAATGCTA")
        self.assertEqual(hsp.sbjct, "AGCGCAAGGCCTGACATAGGAAAATGTTCAGTGAATACTA")
        self.assertEqual(hsp.match, "||||||||||||||||| || |  || ||| | ||| |||")
        self.assertEqual(hsp.score, 40.0)
        self.assertEqual(hsp.expect, 4.57222)
        self.assertEqual(hsp.query_start, 234)
        self.assertEqual(hsp.query_end, 273)
        self.assertEqual(hsp.sbjct_start, 56840340)
        self.assertEqual(hsp.sbjct_end, 56840301)
        self.assertEqual(hsp.frame, (1, -1))
        self.assertEqual(hsp.strand, ("Plus", "Minus"))

        alignment = record.alignments[5]
        self.assertEqual(alignment.hit_id, "gi|372099098|ref|NC_000078.6|")
        self.assertEqual(alignment.accession, "NC_000078")
        self.assertEqual(alignment.length, 120129022)
        self.assertEqual(
            alignment.title,
            "gi|372099098|ref|NC_000078.6| Mus musculus strain C57BL/6J chromosome 12, GRCm38.p4 C57BL/6J",
        )
        self.assertEqual(len(alignment.hsps), 2)
        hsp = alignment.hsps[0]
        self.assertEqual(hsp.align_length, 23)
        self.assertEqual(hsp.identities, 22)
        self.assertEqual(hsp.positives, 22)
        self.assertEqual(hsp.query, "CATCCATTCACACCCAACACAGG")
        self.assertEqual(hsp.sbjct, "CATCCATTCACACCCAGCACAGG")
        self.assertEqual(hsp.match, "|||||||||||||||| ||||||")
        self.assertEqual(hsp.score, 41.0)
        self.assertEqual(hsp.expect, 4.57222)
        self.assertEqual(hsp.query_start, 49)
        self.assertEqual(hsp.query_end, 71)
        self.assertEqual(hsp.sbjct_start, 113030663)
        self.assertEqual(hsp.sbjct_end, 113030685)
        self.assertEqual(hsp.frame, (1, 1))
        self.assertEqual(hsp.strand, ("Plus", "Plus"))
        hsp = alignment.hsps[1]
        self.assertEqual(hsp.align_length, 32)
        self.assertEqual(hsp.identities, 28)
        self.assertEqual(hsp.positives, 28)
        self.assertEqual(hsp.query, "TGTAGCGCAAGGCCTGACATGGGTAGCTGCTC")
        self.assertEqual(hsp.sbjct, "TGTAGCTCTAGGCCTGACATGGGT-GCTGGTC")
        self.assertEqual(hsp.match, "|||||| | ||||||||||||||| |||| ||")
        self.assertEqual(hsp.score, 40.0)
        self.assertEqual(hsp.expect, 4.57222)
        self.assertEqual(hsp.query_start, 231)
        self.assertEqual(hsp.query_end, 262)
        self.assertEqual(hsp.sbjct_start, 108990272)
        self.assertEqual(hsp.sbjct_end, 108990242)
        self.assertEqual(hsp.frame, (1, -1))
        self.assertEqual(hsp.strand, ("Plus", "Minus"))

        alignment = record.alignments[6]
        self.assertEqual(alignment.hit_id, "gi|372099109|ref|NC_000067.6|")
        self.assertEqual(alignment.accession, "NC_000067")
        self.assertEqual(alignment.length, 195471971)
        self.assertEqual(
            alignment.title,
            "gi|372099109|ref|NC_000067.6| Mus musculus strain C57BL/6J chromosome 1, GRCm38.p4 C57BL/6J",
        )
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertEqual(hsp.align_length, 43)
        self.assertEqual(hsp.identities, 35)
        self.assertEqual(hsp.positives, 35)
        self.assertEqual(hsp.query, "GCTCAGCCACAGACATGGTTTGTNACTNTTGAGCTTCTGTTCC")
        self.assertEqual(hsp.sbjct, "GCTCAGCCACATACATGGTTT-TAAGTGTTGAGGCTCT-TTCC")
        self.assertEqual(hsp.match, "||||||||||| ||||||||| | | | |||||  ||| ||||")
        self.assertEqual(hsp.score, 40.0)
        self.assertEqual(hsp.expect, 4.57222)
        self.assertEqual(hsp.query_start, 87)
        self.assertEqual(hsp.query_end, 129)
        self.assertEqual(hsp.sbjct_start, 65190108)
        self.assertEqual(hsp.sbjct_end, 65190148)
        self.assertEqual(hsp.frame, (1, 1))
        self.assertEqual(hsp.strand, ("Plus", "Plus"))

        alignment = record.alignments[7]
        self.assertEqual(alignment.hit_id, "gi|372099101|ref|NC_000075.6|")
        self.assertEqual(alignment.accession, "NC_000075")
        self.assertEqual(alignment.length, 124595110)
        self.assertEqual(
            alignment.title,
            "gi|372099101|ref|NC_000075.6| Mus musculus strain C57BL/6J chromosome 9, GRCm38.p4 C57BL/6J",
        )
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertEqual(hsp.align_length, 47)
        self.assertEqual(hsp.identities, 36)
        self.assertEqual(hsp.positives, 36)
        self.assertEqual(hsp.query, "CAAGGCCTGACATGGGTAGCTGCTCAATAAATGCTAGTNTGTTATTT")
        self.assertEqual(hsp.sbjct, "CAAAGCCTGACAGGTATGACTGCTCAATAAATACTATTTTTTTTTTT")
        self.assertEqual(hsp.match, "||| |||||||| |  |  ||||||||||||| ||| | | || |||")
        self.assertEqual(hsp.score, 40.0)
        self.assertEqual(hsp.expect, 4.57222)
        self.assertEqual(hsp.query_start, 238)
        self.assertEqual(hsp.query_end, 284)
        self.assertEqual(hsp.sbjct_start, 58227241)
        self.assertEqual(hsp.sbjct_end, 58227195)
        self.assertEqual(hsp.frame, (1, -1))
        self.assertEqual(hsp.strand, ("Plus", "Minus"))

        alignment = record.alignments[8]
        self.assertEqual(alignment.hit_id, "gi|372099100|ref|NC_000076.6|")
        self.assertEqual(alignment.accession, "NC_000076")
        self.assertEqual(alignment.length, 130694993)
        self.assertEqual(
            alignment.title,
            "gi|372099100|ref|NC_000076.6| Mus musculus strain C57BL/6J chromosome 10, GRCm38.p4 C57BL/6J",
        )
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertEqual(hsp.align_length, 20)
        self.assertEqual(hsp.identities, 20)
        self.assertEqual(hsp.positives, 20)
        self.assertEqual(hsp.query, "AGCTGCTCAATAAATGCTAG")
        self.assertEqual(hsp.sbjct, "AGCTGCTCAATAAATGCTAG")
        self.assertEqual(hsp.match, "||||||||||||||||||||")
        self.assertEqual(hsp.score, 40.0)
        self.assertEqual(hsp.expect, 4.57222)
        self.assertEqual(hsp.query_start, 255)
        self.assertEqual(hsp.query_end, 274)
        self.assertEqual(hsp.sbjct_start, 119337186)
        self.assertEqual(hsp.sbjct_end, 119337205)
        self.assertEqual(hsp.frame, (1, 1))
        self.assertEqual(hsp.strand, ("Plus", "Plus"))

        alignment = record.alignments[9]
        self.assertEqual(alignment.hit_id, "gi|372099094|ref|NC_000082.6|")
        self.assertEqual(alignment.accession, "NC_000082")
        self.assertEqual(alignment.length, 98207768)
        self.assertEqual(
            alignment.title,
            "gi|372099094|ref|NC_000082.6| Mus musculus strain C57BL/6J chromosome 16, GRCm38.p4 C57BL/6J",
        )
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertEqual(hsp.align_length, 56)
        self.assertEqual(hsp.identities, 43)
        self.assertEqual(hsp.positives, 43)
        self.assertEqual(
            hsp.query, "GGAGGCAAAGAATCCCTACCTCCT-AGGGGTGA-AAGGAAATNAAAATGGAAAGTT"
        )
        self.assertEqual(
            hsp.sbjct, "GGAGGCAAAGAATCCCTACATTGTGACAGCTGATAAAGAAGGTAAAATGGAAAATT"
        )
        self.assertEqual(
            hsp.match, "||||||||||||||||||| |  | |  | ||| || |||   |||||||||| ||"
        )
        self.assertEqual(hsp.score, 40.0)
        self.assertEqual(hsp.expect, 4.57222)
        self.assertEqual(hsp.query_start, 175)
        self.assertEqual(hsp.query_end, 228)
        self.assertEqual(hsp.sbjct_start, 18854780)
        self.assertEqual(hsp.sbjct_end, 18854835)
        self.assertEqual(hsp.frame, (1, 1))
        self.assertEqual(hsp.strand, ("Plus", "Plus"))
        self.assertEqual(len(record.descriptions), 10)
        description = record.descriptions[0]
        self.assertEqual(description.score, 44.0)
        self.assertEqual(description.e, 0.375311)
        self.assertEqual(description.num_alignments, 1)
        description = record.descriptions[1]
        self.assertEqual(description.score, 44.0)
        self.assertEqual(description.e, 0.375311)
        self.assertEqual(description.num_alignments, 1)
        description = record.descriptions[2]
        self.assertEqual(description.score, 43.0)
        self.assertEqual(description.e, 1.30996)
        self.assertEqual(description.num_alignments, 2)
        description = record.descriptions[3]
        self.assertEqual(description.score, 42.0)
        self.assertEqual(description.e, 1.30996)
        self.assertEqual(description.num_alignments, 2)
        description = record.descriptions[4]
        self.assertEqual(description.score, 42.0)
        self.assertEqual(description.e, 1.30996)
        self.assertEqual(description.num_alignments, 2)
        description = record.descriptions[5]
        self.assertEqual(description.score, 41.0)
        self.assertEqual(description.e, 4.57222)
        self.assertEqual(description.num_alignments, 2)
        description = record.descriptions[6]
        self.assertEqual(description.score, 40.0)
        self.assertEqual(description.e, 4.57222)
        self.assertEqual(description.num_alignments, 1)
        description = record.descriptions[7]
        self.assertEqual(description.score, 40.0)
        self.assertEqual(description.e, 4.57222)
        self.assertEqual(description.num_alignments, 1)
        description = record.descriptions[8]
        self.assertEqual(description.score, 40.0)
        self.assertEqual(description.e, 4.57222)
        self.assertEqual(description.num_alignments, 1)
        description = record.descriptions[9]
        self.assertEqual(description.score, 40.0)
        self.assertEqual(description.e, 4.57222)
        self.assertEqual(description.num_alignments, 1)

    def test_xml_2900_blastx_001_v1(self):
        self._test_xml_2900_blastx_001("xml_2900_blastx_001.xml")

    def test_xml_2900_blastx_001_v2(self):
        self._test_xml_2900_blastx_001("xml_2900_blastx_001_v2.xml")

    def _test_xml_2900_blastx_001(self, filename):
        datafile = os.path.join("Blast", filename)
        with open(datafile, "rb") as handle:
            records = NCBIXML.parse(handle)
            record = next(records)

        self.assertEqual(record.application, "BLASTX")
        self.assertEqual(record.version, "2.9.0+")
        self.assertEqual(
            record.query,
            "MAAD0534.RAR Schistosoma mansoni, adult worm (J.C.Parra) Schistosoma mansoni cDNA clone MAAD0534.RAR 5' end similar to S. mansoni actin mRNA, complete cds, mRNA sequence",
        )
        self.assertEqual(record.query_letters, 365)
        self.assertEqual(record.database, "nr")

        alignment = record.alignments[0]
        self.assertEqual(alignment.hit_id, "gi|1530504495|emb|VDM03167.1|")
        self.assertEqual(alignment.accession, "VDM03167")
        self.assertEqual(alignment.length, 132)
        self.assertEqual(
            alignment.title,
            "gi|1530504495|emb|VDM03167.1| unnamed protein product, partial [Schistocephalus solidus]",
        )
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertEqual(hsp.align_length, 108)
        self.assertEqual(hsp.identities, 81)
        self.assertEqual(hsp.positives, 83)
        self.assertEqual(
            hsp.query,
            "MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE",
        )
        self.assertEqual(
            hsp.sbjct,
            "MADEEVQALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE",
        )
        self.assertEqual(
            hsp.match,
            "MADEEVQALVVDNGSGMCKAG       ++  P               G KDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE",
        )
        self.assertEqual(hsp.score, 408.0)
        self.assertEqual(hsp.bits, 161.77)
        self.assertEqual(hsp.expect, 8.11609e-49)
        self.assertEqual(hsp.query_start, 20)
        self.assertEqual(hsp.query_end, 343)
        self.assertEqual(hsp.sbjct_start, 1)
        self.assertEqual(hsp.sbjct_end, 108)
        self.assertEqual(hsp.frame, (2, 0))

        alignment = record.alignments[1]
        self.assertEqual(alignment.hit_id, "gi|510859078|gb|EPB74633.1|")
        self.assertEqual(alignment.accession, "EPB74633")
        self.assertEqual(alignment.length, 119)
        self.assertEqual(
            alignment.title,
            "gi|510859078|gb|EPB74633.1| hypothetical protein ANCCEY_06263 [Ancylostoma ceylanicum]",
        )
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertEqual(hsp.align_length, 115)
        self.assertEqual(hsp.identities, 81)
        self.assertEqual(hsp.positives, 85)
        self.assertEqual(
            hsp.query,
            "MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTELHCIRKP",
        )
        self.assertEqual(
            hsp.sbjct,
            "MCDDDVAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTEAHSILKP",
        )
        self.assertEqual(
            hsp.match,
            "M D++V ALVVDNGSGMCKAG       ++  P               G KDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE H I KP",
        )
        self.assertEqual(hsp.score, 405.0)
        self.assertEqual(hsp.expect, 1.40046e-48)
        self.assertEqual(hsp.query_start, 20)
        self.assertEqual(hsp.query_end, 364)
        self.assertEqual(hsp.sbjct_start, 1)
        self.assertEqual(hsp.sbjct_end, 115)
        self.assertEqual(hsp.frame, (2, 0))

        alignment = record.alignments[2]
        self.assertEqual(alignment.hit_id, "gi|684409690|ref|XP_009175831.1|")
        self.assertEqual(alignment.accession, "XP_009175831")
        self.assertEqual(alignment.length, 246)
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertEqual(hsp.align_length, 108)
        self.assertEqual(hsp.identities, 81)
        self.assertEqual(hsp.positives, 83)
        self.assertEqual(
            hsp.query,
            "MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE",
        )
        self.assertEqual(
            hsp.sbjct,
            "MADEEVQALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE",
        )
        self.assertEqual(
            hsp.match,
            "MADEEVQALVVDNGSGMCKAG       ++  P               G KDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE",
        )
        self.assertEqual(hsp.score, 413.0)
        self.assertEqual(hsp.expect, 4.40404e-48)
        self.assertEqual(hsp.query_start, 20)
        self.assertEqual(hsp.query_end, 343)
        self.assertEqual(hsp.sbjct_start, 1)
        self.assertEqual(hsp.sbjct_end, 108)
        self.assertEqual(hsp.frame, (2, 0))

        alignment = record.alignments[3]
        self.assertEqual(alignment.hit_id, "gi|449710331|gb|EMD49430.1|")
        self.assertEqual(alignment.accession, "EMD49430")
        self.assertEqual(alignment.length, 124)
        self.assertEqual(
            alignment.title,
            "gi|449710331|gb|EMD49430.1| actin, putative, partial [Entamoeba histolytica KU27]",
        )
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertEqual(hsp.align_length, 108)
        self.assertEqual(hsp.identities, 78)
        self.assertEqual(hsp.positives, 81)
        self.assertEqual(
            hsp.query,
            "MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE",
        )
        self.assertEqual(
            hsp.sbjct,
            "MGDEEVQALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHVSVMAGMGQKDAYVGDEAQSKRGILTLKYPIEHGIVNNWDDMEKIWHHTFYNELRVAPEEHPVLLTE",
        )
        self.assertEqual(
            hsp.match,
            "M DEEVQALVVDNGSGMCKAG       ++  P               G KD+YVGDEAQSKRGILTLKYPIEHGIV NWDDMEKIWHHTFYNELRVAPEEHPVLLTE",
        )
        self.assertEqual(hsp.score, 401.0)
        self.assertEqual(hsp.expect, 9.0486e-48)
        self.assertEqual(hsp.query_start, 20)
        self.assertEqual(hsp.query_end, 343)
        self.assertEqual(hsp.sbjct_start, 1)
        self.assertEqual(hsp.sbjct_end, 108)
        self.assertEqual(hsp.frame, (2, 0))

        alignment = record.alignments[4]
        self.assertEqual(alignment.hit_id, "gi|257215766|emb|CAX83035.1|")
        self.assertEqual(alignment.accession, "CAX83035")
        self.assertEqual(alignment.length, 252)
        self.assertEqual(
            alignment.title,
            "gi|257215766|emb|CAX83035.1| Actin-2, partial [Schistosoma japonicum]",
        )
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertEqual(hsp.align_length, 108)
        self.assertEqual(hsp.identities, 81)
        self.assertEqual(hsp.positives, 83)
        self.assertEqual(
            hsp.query,
            "MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE",
        )
        self.assertEqual(
            hsp.sbjct,
            "MADEEVQALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE",
        )
        self.assertEqual(
            hsp.match,
            "MADEEVQALVVDNGSGMCKAG       ++  P               G KDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE",
        )
        self.assertEqual(hsp.score, 411.0)
        self.assertEqual(hsp.expect, 1.00219e-47)
        self.assertEqual(hsp.query_start, 20)
        self.assertEqual(hsp.query_end, 343)
        self.assertEqual(hsp.sbjct_start, 1)
        self.assertEqual(hsp.sbjct_end, 108)
        self.assertEqual(hsp.frame, (2, 0))

        alignment = record.alignments[5]
        self.assertEqual(alignment.hit_id, "gi|1535393712|emb|VDP83060.1|")
        self.assertEqual(alignment.accession, "VDP83060")
        self.assertEqual(alignment.length, 209)
        self.assertEqual(
            alignment.title,
            "gi|1535393712|emb|VDP83060.1| unnamed protein product, partial [Echinostoma caproni]",
        )
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertEqual(hsp.align_length, 108)
        self.assertEqual(hsp.identities, 80)
        self.assertEqual(hsp.positives, 83)
        self.assertEqual(
            hsp.query,
            "MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE",
        )
        self.assertEqual(
            hsp.sbjct,
            "MADDEVQALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE",
        )
        self.assertEqual(
            hsp.match,
            "MAD+EVQALVVDNGSGMCKAG       ++  P               G KDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE",
        )
        self.assertEqual(hsp.score, 407.0)
        self.assertEqual(hsp.expect, 1.16397e-47)
        self.assertEqual(hsp.query_start, 20)
        self.assertEqual(hsp.query_end, 343)
        self.assertEqual(hsp.sbjct_start, 1)
        self.assertEqual(hsp.sbjct_end, 108)
        self.assertEqual(hsp.frame, (2, 0))

        alignment = record.alignments[6]
        self.assertEqual(alignment.hit_id, "gi|312773|emb|CAA50205.1|")
        self.assertEqual(alignment.accession, "CAA50205")
        self.assertEqual(alignment.length, 137)
        self.assertEqual(
            alignment.title,
            "gi|312773|emb|CAA50205.1| actin, partial [Entamoeba histolytica]",
        )
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertEqual(hsp.align_length, 108)
        self.assertEqual(hsp.identities, 78)
        self.assertEqual(hsp.positives, 81)
        self.assertEqual(
            hsp.query,
            "MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE",
        )
        self.assertEqual(
            hsp.sbjct,
            "MGDEEVQALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHVSVMAGMGQKDAYVGDEAQSKRGILTLKYPIEHGIVNNWDDMEKIWHHTFYNELRVAPEEHPVLLTE",
        )
        self.assertEqual(
            hsp.match,
            "M DEEVQALVVDNGSGMCKAG       ++  P               G KD+YVGDEAQSKRGILTLKYPIEHGIV NWDDMEKIWHHTFYNELRVAPEEHPVLLTE",
        )
        self.assertEqual(hsp.score, 401.0)
        self.assertEqual(hsp.expect, 1.25869e-47)
        self.assertEqual(hsp.query_start, 20)
        self.assertEqual(hsp.query_end, 343)
        self.assertEqual(hsp.sbjct_start, 1)
        self.assertEqual(hsp.sbjct_end, 108)
        self.assertEqual(hsp.frame, (2, 0))

        alignment = record.alignments[7]
        self.assertEqual(alignment.hit_id, "gi|1530341495|emb|VDN44756.1|")
        self.assertEqual(alignment.accession, "VDN44756")
        self.assertEqual(alignment.length, 145)
        self.assertEqual(
            alignment.title,
            "gi|1530341495|emb|VDN44756.1| unnamed protein product, partial [Dibothriocephalus latus]",
        )
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertEqual(hsp.align_length, 108)
        self.assertEqual(hsp.identities, 78)
        self.assertEqual(hsp.positives, 82)
        self.assertEqual(
            hsp.query,
            "MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE",
        )
        self.assertEqual(
            hsp.sbjct,
            "MGDEDVQALVIDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE",
        )
        self.assertEqual(
            hsp.match,
            "M DE+VQALV+DNGSGMCKAG       ++  P               G KDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE",
        )
        self.assertEqual(hsp.score, 400.0)
        self.assertEqual(hsp.expect, 1.78336e-47)
        self.assertEqual(hsp.query_start, 20)
        self.assertEqual(hsp.query_end, 343)
        self.assertEqual(hsp.sbjct_start, 1)
        self.assertEqual(hsp.sbjct_end, 108)
        self.assertEqual(hsp.frame, (2, 0))

        alignment = record.alignments[8]
        self.assertEqual(alignment.hit_id, "gi|1524877828|ref|XP_027046469.1|")
        self.assertEqual(alignment.accession, "XP_027046469")
        self.assertEqual(alignment.length, 122)
        self.assertEqual(
            alignment.title,
            "gi|1524877828|ref|XP_027046469.1| actin-1, partial [Pocillopora damicornis]",
        )
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertEqual(hsp.align_length, 108)
        self.assertEqual(hsp.identities, 78)
        self.assertEqual(hsp.positives, 82)
        self.assertEqual(
            hsp.query,
            "MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE",
        )
        self.assertEqual(
            hsp.sbjct,
            "MADEEVAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRIAPEEHPILLTE",
        )
        self.assertEqual(
            hsp.match,
            "MADEEV ALVVDNGSGMCKAG       ++  P               G KDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELR+APEEHP+LLTE",
        )
        self.assertEqual(hsp.score, 398.0)
        self.assertEqual(hsp.expect, 1.93331e-47)
        self.assertEqual(hsp.query_start, 20)
        self.assertEqual(hsp.query_end, 343)
        self.assertEqual(hsp.sbjct_start, 1)
        self.assertEqual(hsp.sbjct_end, 108)
        self.assertEqual(hsp.frame, (2, 0))

        alignment = record.alignments[9]
        self.assertEqual(alignment.hit_id, "gi|1524877860|ref|XP_027046487.1|")
        self.assertEqual(alignment.accession, "XP_027046487")
        self.assertEqual(alignment.length, 134)
        self.assertEqual(
            alignment.title,
            "gi|1524877860|ref|XP_027046487.1| actin-1-like [Pocillopora damicornis]",
        )
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertEqual(hsp.align_length, 108)
        self.assertEqual(hsp.identities, 79)
        self.assertEqual(hsp.positives, 82)
        self.assertEqual(
            hsp.query,
            "MADEEVQALVVDNGSGMCKAGIRW**CTKSSIPFHRWTTSTSRCDGWYGSKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE",
        )
        self.assertEqual(
            hsp.sbjct,
            "MADEDVAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE",
        )
        self.assertEqual(
            hsp.match,
            "MADE+V ALVVDNGSGMCKAG       ++  P               G KDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE",
        )
        self.assertEqual(hsp.score, 399.0)
        self.assertEqual(hsp.expect, 2.36088e-47)
        self.assertEqual(hsp.query_start, 20)
        self.assertEqual(hsp.query_end, 343)
        self.assertEqual(hsp.sbjct_start, 1)
        self.assertEqual(hsp.sbjct_end, 108)
        self.assertEqual(hsp.frame, (2, 0))

        self.assertEqual(len(record.descriptions), 10)
        description = record.descriptions[0]
        self.assertEqual(
            description.title,
            "gi|1530504495|emb|VDM03167.1| unnamed protein product, partial [Schistocephalus solidus]",
        )
        self.assertEqual(description.score, 408.0)
        self.assertEqual(description.e, 8.11609e-49)
        self.assertEqual(description.num_alignments, 1)
        description = record.descriptions[1]
        self.assertEqual(
            description.title,
            "gi|510859078|gb|EPB74633.1| hypothetical protein ANCCEY_06263 [Ancylostoma ceylanicum]",
        )
        self.assertEqual(description.score, 405.0)
        self.assertEqual(description.e, 1.40046e-48)
        self.assertEqual(description.num_alignments, 1)
        description = record.descriptions[2]
        self.assertEqual(description.score, 413.0)
        self.assertEqual(description.e, 4.40404e-48)
        self.assertEqual(description.num_alignments, 1)
        description = record.descriptions[3]
        self.assertEqual(
            description.title,
            "gi|449710331|gb|EMD49430.1| actin, putative, partial [Entamoeba histolytica KU27]",
        )
        self.assertEqual(description.score, 401.0)
        self.assertEqual(description.e, 9.0486e-48)
        self.assertEqual(description.num_alignments, 1)
        description = record.descriptions[4]
        self.assertEqual(
            description.title,
            "gi|257215766|emb|CAX83035.1| Actin-2, partial [Schistosoma japonicum]",
        )
        self.assertEqual(description.score, 411.0)
        self.assertEqual(description.e, 1.00219e-47)
        self.assertEqual(description.num_alignments, 1)
        description = record.descriptions[5]
        self.assertEqual(
            description.title,
            "gi|1535393712|emb|VDP83060.1| unnamed protein product, partial [Echinostoma caproni]",
        )
        self.assertEqual(description.score, 407.0)
        self.assertEqual(description.e, 1.16397e-47)
        self.assertEqual(description.num_alignments, 1)
        description = record.descriptions[6]
        self.assertEqual(
            description.title,
            "gi|312773|emb|CAA50205.1| actin, partial [Entamoeba histolytica]",
        )
        self.assertEqual(description.score, 401.0)
        self.assertEqual(description.e, 1.25869e-47)
        self.assertEqual(description.num_alignments, 1)
        description = record.descriptions[7]
        self.assertEqual(
            description.title,
            "gi|1530341495|emb|VDN44756.1| unnamed protein product, partial [Dibothriocephalus latus]",
        )
        self.assertEqual(description.score, 400.0)
        self.assertEqual(description.e, 1.78336e-47)
        self.assertEqual(description.num_alignments, 1)
        description = record.descriptions[8]
        self.assertEqual(
            description.title,
            "gi|1524877828|ref|XP_027046469.1| actin-1, partial [Pocillopora damicornis]",
        )
        self.assertEqual(description.score, 398.0)
        self.assertEqual(description.e, 1.93331e-47)
        self.assertEqual(description.num_alignments, 1)
        description = record.descriptions[9]
        self.assertEqual(
            description.title,
            "gi|1524877860|ref|XP_027046487.1| actin-1-like [Pocillopora damicornis]",
        )
        self.assertEqual(description.score, 399.0)
        self.assertEqual(description.e, 2.36088e-47)
        self.assertEqual(description.num_alignments, 1)

    def test_xml_2900_tblastn_001_v1(self):
        self._test_xml_2900_tblastn_001("xml_2900_tblastn_001.xml")

    def test_xml_2900_tblastn_001_v2(self):
        self._test_xml_2900_tblastn_001("xml_2900_tblastn_001_v2.xml")

    def _test_xml_2900_tblastn_001(self, filename):
        datafile = os.path.join("Blast", filename)
        with open(datafile, "rb") as handle:
            records = NCBIXML.parse(handle)
            record = next(records)
            self.assertRaises(StopIteration, next, records)

        self.assertEqual(record.application, "TBLASTN")
        self.assertEqual(record.version, "2.9.0+")
        self.assertEqual(record.query, "tim [Helicobacter acinonychis str. Sheeba]")
        self.assertEqual(record.query_letters, 234)
        self.assertEqual(record.database, "nr")

        alignment = record.alignments[4]
        self.assertEqual(alignment.hit_id, "gi|1143706535|gb|CP018823.1|")
        self.assertEqual(alignment.accession, "CP018823")
        self.assertEqual(alignment.length, 1618480)
        self.assertEqual(
            alignment.title,
            "gi|1143706535|gb|CP018823.1| Helicobacter pylori strain PMSS1 complete genome",
        )
        self.assertEqual(len(alignment.hsps), 1)
        hsp = alignment.hsps[0]
        self.assertEqual(hsp.align_length, 234)
        self.assertEqual(hsp.identities, 218)
        self.assertEqual(hsp.positives, 223)
        self.assertEqual(
            hsp.query,
            "MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHCDRVFVFPDFLGLLPNAFLHFTLGVQNAYPKDCGAFTGEITSKHLEELKINTLLIGHSERRVLLKESPNFLKEKFDFFKDKKFKIVYCIGEDLKTREKGLGAVKEFLNEQLENIDLDYQNLIVAYEPIWAIGTGKSASLEDIYLTHGFLKQHLNQKMPLLYGGSVNTQNAKEILGIDSVDGLLIGSTSLELENFKTIISFL",
        )
        self.assertEqual(
            hsp.sbjct,
            "MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQHFDRVFVFPDFLGLLPNSFLHFTLGVQNAYPRDCGAFTGEITSKHLEELKIHTLLIGHSERRVLLKESPSFLKEKFDFFKDKNFKIVYCIGEDLTTREKGFKAVKEFLNEQLENIDLNYSNLIVAYEPIWAIGTKKSASLEDIYLTHGFLKQILNQKTPLLYGGSVNTQNAKEILGIDSVDGLLIGSASWELENFKTIISFL",
        )
        self.assertEqual(
            hsp.match,
            "MTKIAMANFKSAMPIFKSHAYLKELEKTLKPQH DRVFVFPDFLGLLPN+FLHFTLGVQNAYP+DCGAFTGEITSKHLEELKI+TLLIGHSERRVLLKESP+FLKEKFDFFKDK FKIVYCIGEDL TREKG  AVKEFLNEQLENIDL+Y NLIVAYEPIWAIGT KSASLEDIYLTHGFLKQ LNQK PLLYGGSVNTQNAKEILGIDSVDGLLIGS S ELENFKTIISFL",
        )
        self.assertEqual(hsp.score, 1136.0)
        self.assertEqual(hsp.bits, 442.195)
        self.assertEqual(hsp.expect, 2.08707e-139)
        self.assertEqual(hsp.query_start, 1)
        self.assertEqual(hsp.query_end, 234)
        self.assertEqual(hsp.sbjct_start, 190464)
        self.assertEqual(hsp.sbjct_end, 191165)
        self.assertEqual(hsp.frame, (0, 3))

        self.assertEqual(len(record.descriptions), 10)
        description = record.descriptions[4]
        self.assertEqual(
            description.title,
            "gi|1143706535|gb|CP018823.1| Helicobacter pylori strain PMSS1 complete genome",
        )
        self.assertEqual(description.score, 1136.0)
        self.assertEqual(description.e, 2.08707e-139)
        self.assertEqual(description.num_alignments, 1)

    def test_xml_2900_tblastx_001_v1(self):
        self._test_xml_2900_tblastx_001("xml_2900_tblastx_001.xml")

    def test_xml_2900_tblastx_001_v2(self):
        self._test_xml_2900_tblastx_001("xml_2900_tblastx_001_v2.xml")

    def _test_xml_2900_tblastx_001(self, filename):
        datafile = os.path.join("Blast", filename)
        with open(datafile, "rb") as handle:
            records = NCBIXML.parse(handle)
            record = next(records)

        self.assertEqual(record.application, "TBLASTX")
        self.assertEqual(record.version, "2.9.0+")
        self.assertEqual(
            record.query,
            "MAAD0534.RAR Schistosoma mansoni, adult worm (J.C.Parra) Schistosoma mansoni cDNA clone MAAD0534.RAR 5' end similar to S. mansoni actin mRNA, complete cds, mRNA sequence",
        )
        self.assertEqual(record.query_letters, 365)
        self.assertEqual(record.database, "nr")

        alignment = record.alignments[4]
        self.assertEqual(alignment.hit_id, "gi|1590279025|ref|XM_028307351.1|")
        self.assertEqual(alignment.accession, "XM_028307351")
        self.assertEqual(alignment.length, 1599)
        self.assertEqual(
            alignment.title,
            "gi|1590279025|ref|XM_028307351.1| PREDICTED: Ostrinia furnacalis actin, cytoplasmic A3a (LOC114354791), mRNA",
        )
        self.assertEqual(len(alignment.hsps), 12)
        hsp = alignment.hsps[0]
        self.assertEqual(hsp.align_length, 60)
        self.assertEqual(hsp.identities, 59)
        self.assertEqual(hsp.positives, 59)
        self.assertEqual(
            hsp.query, "GSKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE"
        )
        self.assertEqual(
            hsp.sbjct, "GQKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE"
        )
        self.assertEqual(
            hsp.match, "G KDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTE"
        )
        self.assertEqual(hsp.score, 325.0)
        self.assertEqual(hsp.expect, 7.83523e-61)
        self.assertEqual(hsp.query_start, 164)
        self.assertEqual(hsp.query_end, 343)
        self.assertEqual(hsp.sbjct_start, 389)
        self.assertEqual(hsp.sbjct_end, 568)
        self.assertEqual(hsp.frame, (2, 2))
        hsp = alignment.hsps[1]
        self.assertEqual(hsp.align_length, 36)
        self.assertEqual(hsp.identities, 33)
        self.assertEqual(hsp.positives, 33)
        self.assertEqual(hsp.query, "GCAKLGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQK")
        self.assertEqual(hsp.sbjct, "GMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQK")
        self.assertEqual(hsp.match, "G  K GFAGDDAPRAVFPSIVGRPRHQGVMVGMGQK")
        self.assertEqual(hsp.score, 174.0)
        self.assertEqual(hsp.expect, 7.83523e-61)
        self.assertEqual(hsp.query_start, 66)
        self.assertEqual(hsp.query_end, 173)
        self.assertEqual(hsp.sbjct_start, 290)
        self.assertEqual(hsp.sbjct_end, 397)
        self.assertEqual(hsp.frame, (3, 2))
        hsp = alignment.hsps[2]
        self.assertEqual(hsp.align_length, 21)
        self.assertEqual(hsp.identities, 19)
        self.assertEqual(hsp.positives, 19)
        self.assertEqual(hsp.query, "MADEEVQALVVDNGSGMCKAG")
        self.assertEqual(hsp.sbjct, "MCDEEVAALVVDNGSGMCKAG")
        self.assertEqual(hsp.match, "M DEEV ALVVDNGSGMCKAG")
        self.assertEqual(hsp.score, 97.0)
        self.assertEqual(hsp.expect, 7.83523e-61)
        self.assertEqual(hsp.query_start, 20)
        self.assertEqual(hsp.query_end, 82)
        self.assertEqual(hsp.sbjct_start, 245)
        self.assertEqual(hsp.sbjct_end, 307)
        self.assertEqual(hsp.frame, (2, 2))
        hsp = alignment.hsps[3]
        self.assertEqual(hsp.align_length, 61)
        self.assertEqual(hsp.identities, 48)
        self.assertEqual(hsp.positives, 52)
        self.assertEqual(
            hsp.query, "SSVNRTGCSSGATRNSL*NV*CQIFSMSSQFVTIPCSIGYFSVRIPRFDCASSPT*LSFDP"
        )
        self.assertEqual(
            hsp.sbjct, "ASVRRTGCSSGATRSSL*KV*CQIFSMSSQFVTIPCSMGYLSVRMPLLLWASSPT*ESFCP"
        )
        self.assertEqual(
            hsp.match, "+SV RTGCSSGATR+SL* V*CQIFSMSSQFVTIPCS+GY SVR+P    ASSPT* SF P"
        )
        self.assertEqual(hsp.score, 225.0)
        self.assertEqual(hsp.expect, 3.26469e-39)
        self.assertEqual(hsp.query_start, 163)
        self.assertEqual(hsp.query_end, 345)
        self.assertEqual(hsp.sbjct_start, 388)
        self.assertEqual(hsp.sbjct_end, 570)
        self.assertEqual(hsp.frame, (-3, -1))
        hsp = alignment.hsps[4]
        self.assertEqual(hsp.align_length, 36)
        self.assertEqual(hsp.identities, 24)
        self.assertEqual(hsp.positives, 28)
        self.assertEqual(hsp.query, "F*PIPTITP*CRGRPTMEGNTALGASSPANPSFAHP")
        self.assertEqual(hsp.sbjct, "FCPMPTITPWWRGRPTIDGNTARGASSPAKPALHIP")
        self.assertEqual(hsp.match, "F P+PTITP  RGRPT++GNTA GASSPA P+   P")
        self.assertEqual(hsp.score, 121.0)
        self.assertEqual(hsp.expect, 3.26469e-39)
        self.assertEqual(hsp.query_start, 65)
        self.assertEqual(hsp.query_end, 172)
        self.assertEqual(hsp.sbjct_start, 289)
        self.assertEqual(hsp.sbjct_end, 396)
        self.assertEqual(hsp.frame, (-2, -1))
        hsp = alignment.hsps[5]
        self.assertEqual(hsp.align_length, 26)
        self.assertEqual(hsp.identities, 18)
        self.assertEqual(hsp.positives, 22)
        self.assertEqual(hsp.query, "PALHIPDPLSTTRA*TSSSAMIIFQL")
        self.assertEqual(hsp.sbjct, "PALHIPDPLSTTNAATSSSHILVYLL")
        self.assertEqual(hsp.match, "PALHIPDPLSTT A TSSS ++++ L")
        self.assertEqual(hsp.score, 91.0)
        self.assertEqual(hsp.expect, 3.26469e-39)
        self.assertEqual(hsp.query_start, 4)
        self.assertEqual(hsp.query_end, 81)
        self.assertEqual(hsp.sbjct_start, 229)
        self.assertEqual(hsp.sbjct_end, 306)
        self.assertEqual(hsp.frame, (-3, -1))
        hsp = alignment.hsps[6]
        self.assertEqual(hsp.align_length, 57)
        self.assertEqual(hsp.identities, 34)
        self.assertEqual(hsp.positives, 43)
        self.assertEqual(
            hsp.query, "GQQDRVFFWSHTQFIVECVMPDLLHVIPVRHNTVFDWVFQCENTTFRLCFITDVAVF"
        )
        self.assertEqual(
            hsp.sbjct, "GKKDWVFLGSDTEFIVEGVMPDLLHVIPVCDDSVFDGVFECEDASLALGLISDVGVF"
        )
        self.assertEqual(
            hsp.match, "G++D VF  S T+FIVE VMPDLLHVIPV  ++VFD VF+CE+ +  L  I+DV VF"
        )
        self.assertEqual(hsp.score, 167.0)
        self.assertEqual(hsp.expect, 2.27039e-16)
        self.assertEqual(hsp.query_start, 170)
        self.assertEqual(hsp.query_end, 340)
        self.assertEqual(hsp.sbjct_start, 395)
        self.assertEqual(hsp.sbjct_end, 565)
        self.assertEqual(hsp.frame, (-2, -3))
        hsp = alignment.hsps[7]
        self.assertEqual(hsp.align_length, 32)
        self.assertEqual(hsp.identities, 13)
        self.assertEqual(hsp.positives, 23)
        self.assertEqual(hsp.query, "LTHTNHHTLMSRSSNDGREYCSWCIITSESQL")
        self.assertEqual(hsp.sbjct, "LSHANHHTLVARAAHDRWEHGAGRVIARETGL")
        self.assertEqual(hsp.match, "L+H NHHTL++R+++D  E+ +  +I  E+ L")
        self.assertEqual(hsp.score, 69.0)
        self.assertEqual(hsp.expect, 2.27039e-16)
        self.assertEqual(hsp.query_start, 75)
        self.assertEqual(hsp.query_end, 170)
        self.assertEqual(hsp.sbjct_start, 299)
        self.assertEqual(hsp.sbjct_end, 394)
        self.assertEqual(hsp.frame, (-1, -3))
        hsp = alignment.hsps[8]
        self.assertEqual(hsp.align_length, 64)
        self.assertEqual(hsp.identities, 36)
        self.assertEqual(hsp.positives, 44)
        self.assertEqual(
            hsp.query,
            "KRQLRR**STIETWYSHTEIPNRTRYCDELG*HGEDLASHILQ*IACGSRRTPCPVDRAPLYPK",
        )
        self.assertEqual(
            hsp.sbjct,
            "KRLLRRR*GPEQERHPHTQIPHRTRNRHKLG*HGEDLASHLLQ*TPCRSRGTPSPSYRSPPEPQ",
        )
        self.assertEqual(
            hsp.match,
            "KR LRR *   +  + HT+IP+RTR   +LG*HGEDLASH+LQ*  C SR TP P  R+P  P+",
        )
        self.assertEqual(hsp.score, 160.0)
        self.assertEqual(hsp.expect, 2.96635e-13)
        self.assertEqual(hsp.query_start, 169)
        self.assertEqual(hsp.query_end, 360)
        self.assertEqual(hsp.sbjct_start, 394)
        self.assertEqual(hsp.sbjct_end, 585)
        self.assertEqual(hsp.frame, (1, 1))
        hsp = alignment.hsps[9]
        self.assertEqual(hsp.align_length, 20)
        self.assertEqual(hsp.identities, 8)
        self.assertEqual(hsp.positives, 12)
        self.assertEqual(hsp.query, "IPFHRWTTSTSRCDGWYGSK")
        self.assertEqual(hsp.sbjct, "VPIDRGPPAPPGCDGWHGTK")
        self.assertEqual(hsp.match, "+P  R   +   CDGW+G+K")
        self.assertEqual(hsp.score, 53.0)
        self.assertEqual(hsp.expect, 2.96635e-13)
        self.assertEqual(hsp.query_start, 113)
        self.assertEqual(hsp.query_end, 172)
        self.assertEqual(hsp.sbjct_start, 337)
        self.assertEqual(hsp.sbjct_end, 396)
        self.assertEqual(hsp.frame, (2, 1))
        hsp = alignment.hsps[10]
        self.assertEqual(hsp.align_length, 48)
        self.assertEqual(hsp.identities, 29)
        self.assertEqual(hsp.positives, 33)
        self.assertEqual(hsp.query, "FGYNGARSTGQGVLLEPHAIHCRMCDARSSPCHPSSSQYRVRLGISV*")
        self.assertEqual(hsp.sbjct, "WGSGGLR*EGLGVPRERHGVHCRRCDARSSPCHPSL*RFRVRWGI*V*")
        self.assertEqual(hsp.match, "+G  G R  G GV  E H +HCR CDARSSPCHPS  ++RVR GI V*")
        self.assertEqual(hsp.score, 142.0)
        self.assertEqual(hsp.expect, 2.61074e-08)
        self.assertEqual(hsp.query_start, 216)
        self.assertEqual(hsp.query_end, 359)
        self.assertEqual(hsp.sbjct_start, 441)
        self.assertEqual(hsp.sbjct_end, 584)
        self.assertEqual(hsp.frame, (-1, -2))
        hsp = alignment.hsps[11]
        self.assertEqual(hsp.align_length, 58)
        self.assertEqual(hsp.identities, 33)
        self.assertEqual(hsp.positives, 37)
        self.assertEqual(
            hsp.query, "QKTATSVMKHNRNVVFSH*NTQSNTVL*RTGMTWRRSGITHSTMNCVWLQKNTLSC*P"
        )
        self.assertEqual(
            hsp.sbjct, "KKTPTSEMRPRAREASSHSNTPSNTESSQTGMTWRRSGITPSTMNSVSLPRNTQSFLP"
        )
        self.assertEqual(
            hsp.match, "+KT TS M+       SH NT SNT   +TGMTWRRSGIT STMN V L +NT S  P"
        )
        self.assertEqual(hsp.score, 141.0)
        self.assertEqual(hsp.expect, 3.58672e-08)
        self.assertEqual(hsp.query_start, 168)
        self.assertEqual(hsp.query_end, 341)
        self.assertEqual(hsp.sbjct_start, 393)
        self.assertEqual(hsp.sbjct_end, 566)
        self.assertEqual(hsp.frame, (3, 3))
        self.assertEqual(len(record.descriptions), 10)
        description = record.descriptions[4]
        self.assertEqual(
            description.title,
            "gi|1590279025|ref|XM_028307351.1| PREDICTED: Ostrinia furnacalis actin, cytoplasmic A3a "
            "(LOC114354791), mRNA",
        )
        self.assertEqual(description.score, 325.0)
        self.assertEqual(description.e, 7.83523e-61)
        self.assertEqual(description.num_alignments, 12)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
