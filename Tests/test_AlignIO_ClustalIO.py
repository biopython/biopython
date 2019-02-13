# Copyright 2006-2014 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.AlignIO.ClustalIO"""

import unittest

from Bio._py3k import StringIO

from Bio.AlignIO.ClustalIO import ClustalIterator, ClustalWriter

# This is a truncated version of the example in Tests/cw02.aln
# Notice the inclusion of sequence numbers (right hand side)
aln_example1 = \
"""CLUSTAL W (1.81) multiple sequence alignment


gi|4959044|gb|AAD34209.1|AF069      MENSDSNDKGSDQSAAQRRSQMDRLDREEAFYQFVNNLSEEDYRLMRDNN 50
gi|671626|emb|CAA85685.1|           ---------MSPQTETKASVGFKAGVKEYKLTYYTPEYETKDTDILAAFR 41
                                              * *: ::    :.   :*  :  :. : . :*  ::   .

gi|4959044|gb|AAD34209.1|AF069      LLGTPGESTEEELLRRLQQIKEGPPPQSPDENRAGESSDDVTNSDSIIDW 100
gi|671626|emb|CAA85685.1|           VTPQPG-----------------VPPEEAGAAVAAESSTGT--------- 65
                                    :   **                  **:...   *.*** ..         

gi|4959044|gb|AAD34209.1|AF069      LNSVRQTGNTTRSRQRGNQSWRAVSRTNPNSGDFRFSLEINVNRNNGSQT 150
gi|671626|emb|CAA85685.1|           WTTVWTDGLTSLDRYKG-----RCYHIEPVPG------------------ 92
                                     .:*   * *: .* :*        : :* .*                  

gi|4959044|gb|AAD34209.1|AF069      SENESEPSTRRLSVENMESSSQRQMENSASESASARPSRAERNSTEAVTE 200
gi|671626|emb|CAA85685.1|           -EKDQCICYVAYPLDLFEEGSVTNMFTSIVGNVFGFKALRALRLEDLRIP 141
                                     *::.  .    .:: :*..*  :* .*   .. .  :    .  :    

gi|4959044|gb|AAD34209.1|AF069      VPTTRAQRRA 210
gi|671626|emb|CAA85685.1|           VAYVKTFQGP 151
                                    *. .:: : .
"""  # noqa : W291

# This example is a truncated version of the dataset used here:
# http://virgil.ruc.dk/kurser/Sekvens/Treedraw.htm
# with the last record repeated twice (deliberate toture test)
aln_example2 = \
"""CLUSTAL X (1.83) multiple sequence alignment


V_Harveyi_PATH                 --MKNWIKVAVAAIA--LSAA------------------TVQAATEVKVG
B_subtilis_YXEM                MKMKKWTVLVVAALLAVLSACG------------NGNSSSKEDDNVLHVG
B_subtilis_GlnH_homo_YCKK      MKKALLALFMVVSIAALAACGAGNDNQSKDNAKDGDLWASIKKKGVLTVG
YA80_HAEIN                     MKKLLFTTALLTGAIAFSTF-----------SHAGEIADRVEKTKTLLVG
FLIY_ECOLI                     MKLAHLGRQALMGVMAVALVAG---MSVKSFADEG-LLNKVKERGTLLVG
E_coli_GlnH                    --MKSVLKVSLAALTLAFAVS------------------SHAADKKLVVA
Deinococcus_radiodurans        -MKKSLLSLKLSGLLVPSVLALS--------LSACSSPSSTLNQGTLKIA
HISJ_E_COLI                    MKKLVLSLSLVLAFSSATAAF-------------------AAIPQNIRIG
HISJ_E_COLI                    MKKLVLSLSLVLAFSSATAAF-------------------AAIPQNIRIG
                                         : .                                 : :.

V_Harveyi_PATH                 MSGRYFPFTFVKQ--DKLQGFEVDMWDEIGKRNDYKIEYVTANFSGLFGL
B_subtilis_YXEM                ATGQSYPFAYKEN--GKLTGFDVEVMEAVAKKIDMKLDWKLLEFSGLMGE
B_subtilis_GlnH_homo_YCKK      TEGTYEPFTYHDKDTDKLTGYDVEVITEVAKRLGLKVDFKETQWGSMFAG
YA80_HAEIN                     TEGTYAPFTFHDK-SGKLTGFDVEVIRKVAEKLGLKVEFKETQWDAMYAG
FLIY_ECOLI                     LEGTYPPFSFQGD-DGKLTGFEVEFAQQLAKHLGVEASLKPTKWDGMLAS
E_coli_GlnH                    TDTAFVPFEFKQG--DKYVGFDVDLWAAIAKELKLDYELKPMDFSGIIPA
Deinococcus_radiodurans        MEGTYPPFTSKNE-QGELVGFDVDIAKAVAQKLNLKPEFVLTEWSGILAG
HISJ_E_COLI                    TDPTYAPFESKNS-QGELVGFDIDLAKELCKRINTQCTFVENPLDALIPS
HISJ_E_COLI                    TDPTYAPFESKNS-QGELVGFDIDLAKELCKRINTQCTFVENPLDALIPS
                                     **       .:  *::::.   : :.   .        ..:   

V_Harveyi_PATH                 LETGRIDTISNQITMTDARKAKYLFADPYVVDG-AQI
B_subtilis_YXEM                LQTGKLDTISNQVAVTDERKETYNFTKPYAYAG-TQI
B_subtilis_GlnH_homo_YCKK      LNSKRFDVVANQVG-KTDREDKYDFSDKYTTSR-AVV
YA80_HAEIN                     LNAKRFDVIANQTNPSPERLKKYSFTTPYNYSG-GVI
FLIY_ECOLI                     LDSKRIDVVINQVTISDERKKKYDFSTPYTISGIQAL
E_coli_GlnH                    LQTKNVDLALAGITITDERKKAIDFSDGYYKSG-LLV
Deinococcus_radiodurans        LQANKYDVIVNQVGITPERQNSIGFSQPYAYSRPEII
HISJ_E_COLI                    LKAKKIDAIMSSLSITEKRQQEIAFTDKLYAADSRLV
HISJ_E_COLI                    LKAKKIDAIMSSLSITEKRQQEIAFTDKLYAADSRLV
                               *.: . *        .  *     *:          :

"""  # noqa : W291


aln_example3 = \
"""CLUSTAL 2.0.9 multiple sequence alignment


Test1seq             ------------------------------------------------------------
AT3G20900.1-SEQ      ATGAACAAAGTAGCGAGGAAGAACAAAACATCAGGTGAACAAAAAAAAAACTCAATCCAC
AT3G20900.1-CDS      ------------------------------------------------------------
                                                                                 

Test1seq             -----AGTTACAATAACTGACGAAGCTAAGTAGGCTACTAATTAACGTCATCAACCTAAT
AT3G20900.1-SEQ      ATCAAAGTTACAATAACTGACGAAGCTAAGTAGGCTAGAAATTAAAGTCATCAACCTAAT
AT3G20900.1-CDS      ------------------------------------------------------------
                                                                                 

Test1seq             ACATAGCACTTAGAAAAAAGTGAAGTAAGAAAATATAAAATAATAAAAGGGTGGGTTATC
AT3G20900.1-SEQ      ACATAGCACTTAGAAAAAAGTGAAGCAAGAAAATATAAAATAATAAAAGGGTGGGTTATC
AT3G20900.1-CDS      ------------------------------------------------------------
                                                                                 

Test1seq             AATTGATAGTGTAAATCATCGTATTCCGGTGATATACCCTACCACAAAAACTCAAACCGA
AT3G20900.1-SEQ      AATTGATAGTGTAAATCATAGTTGATTTTTGATATACCCTACCACAAAAACTCAAACCGA
AT3G20900.1-CDS      ------------------------------------------------------------
                                                                                 

Test1seq             CTTGATTCAAATCATCTCAATAAATTAGCGCCAAAATAATGAAAAAAATAATAACAAACA
AT3G20900.1-SEQ      CTTGATTCAAATCATCTCAAAAAACAAGCGCCAAAATAATGAAAAAAATAATAACAAAAA
AT3G20900.1-CDS      ------------------------------------------------------------
                                                                                 

Test1seq             AAAACAAACCAAAATAAGAAAAAACATTACGCAAAACATAATAATTTACTCTTCGTTATT
AT3G20900.1-SEQ      CAAACAAACCAAAATAAGAAAAAACATTACGCAAAACATAATAATTTACTCTTCGTTATT
AT3G20900.1-CDS      ------------------------------------------------------------
                                                                                 

Test1seq             GTATTAACAAATCAAAGAGCTGAATTTTGATCACCTGCTAATACTACTTTCTGTATTGAT
AT3G20900.1-SEQ      GTATTAACAAATCAAAGAGATGAATTTTGATCACCTGCTAATACTACTTTCTGTATTGAT
AT3G20900.1-CDS      ------------------------------------------------------------
                                                                                 

Test1seq             CCTATATCAACGTAAACAAAGATACTAATAATTAACTAAAAGTACGTTCATCGATCGTGT
AT3G20900.1-SEQ      CCTATATCAAAAAAAAAAAAGATACTAATAATTAACTAAAAGTACGTTCATCGATCGTGT
AT3G20900.1-CDS      ------------------------------------------------------ATGAAC
                                                                             *   

Test1seq             TCGTTGACGAAGAAGAGCTCTATCTCCGGCGGAGCAAAGAAAACGATCTGTCTCCGTCGT
AT3G20900.1-SEQ      GCGTTGACGAAGAAGAGCTCTATCTCCGGCGGAGCAAAGAAAACGATCTGTCTCCGTCGT
AT3G20900.1-CDS      AAAGTAGCGAGGAAGAACAAAACATC------AGCAAAGAAAACGATCTGTCTCCGTCGT
                         *  *** ***** *   *  **      ****************************

Test1seq             AACACACGGTCGCTAGAGAAACTTTGCTTCTTCGGCGCCGGTGGACACGTCAGCATCTCC
AT3G20900.1-SEQ      AACACACAGTTTTTCGAGACCCTTTGCTTCTTCGGCGCCGGTGGACACGTCAGCATCTCC
AT3G20900.1-CDS      AACACACAGTTTTTCGAGACCCTTTGCTTCTTCGGCGCCGGTGGACACGTCAGCATCTCC
                     ******* **   * ****  ***************************************

Test1seq             GGTATCCTAGACTTCTTGGCTTTCGGGGTACAACAACCGCGTGGTGACGTCAGCACCGCT
AT3G20900.1-SEQ      GGTATCCTAGACTTCTTGGCTTTCGGGGTACAACAACCGCCTGGTGACGTCAGCACCGCT
AT3G20900.1-CDS      GGTATCCTAGACTTCTTGGCTTTCGGGGTACAACAACCGCCTGGTGACGTCAGCACCGCT
                     **************************************** *******************

Test1seq             GCTGGGGATGGAGAGGGAACAGAGTT-
AT3G20900.1-SEQ      GCTGGGGATGGAGAGGGAACAGAGTAG
AT3G20900.1-CDS      GCTGGGGATGGAGAGGGAACAGAGTAG
                     *************************  
"""  # noqa : W291

aln_example4 = \
"""Kalign (2.0) alignment in ClustalW format

Test1seq             GCTGGGGATGGAGAGGGAACAGAGTT-
AT3G20900.1-SEQ      GCTGGGGATGGAGAGGGAACAGAGTAG

"""


class TestClustalIO(unittest.TestCase):

    def test_one(self):
        alignments = list(ClustalIterator(StringIO(aln_example1)))
        self.assertEqual(1, len(alignments))
        self.assertEqual(alignments[0]._version, "1.81")
        alignment = alignments[0]
        self.assertEqual(2, len(alignment))
        self.assertEqual(alignment[0].id, "gi|4959044|gb|AAD34209.1|AF069")
        self.assertEqual(alignment[1].id, "gi|671626|emb|CAA85685.1|")
        self.assertEqual(str(alignment[0].seq),
                         "MENSDSNDKGSDQSAAQRRSQMDRLDREEAFYQFVNNLSEEDYRLMRDNN"
                         "LLGTPGESTEEELLRRLQQIKEGPPPQSPDENRAGESSDDVTNSDSIIDW"
                         "LNSVRQTGNTTRSRQRGNQSWRAVSRTNPNSGDFRFSLEINVNRNNGSQT"
                         "SENESEPSTRRLSVENMESSSQRQMENSASESASARPSRAERNSTEAVTE"
                         "VPTTRAQRRA")

    def test_two(self):
        alignments = list(ClustalIterator(StringIO(aln_example2)))
        self.assertEqual(1, len(alignments))
        self.assertEqual(alignments[0]._version, "1.83")
        alignment = alignments[0]
        self.assertEqual(9, len(alignment))
        self.assertEqual(alignment[-1].id, "HISJ_E_COLI")
        self.assertEqual(str(alignment[-1].seq),
                         "MKKLVLSLSLVLAFSSATAAF-------------------AAIPQNIRIG"
                         "TDPTYAPFESKNS-QGELVGFDIDLAKELCKRINTQCTFVENPLDALIPS"
                         "LKAKKIDAIMSSLSITEKRQQEIAFTDKLYAADSRLV")

    def test_cat_one_two(self):
        alignments = list(ClustalIterator(StringIO(aln_example2 + aln_example1)))
        self.assertEqual(2, len(alignments))
        self.assertEqual(9, len(alignments[0]))
        self.assertEqual(137, alignments[0].get_alignment_length())
        self.assertEqual(2, len(alignments[1]))
        self.assertEqual(210, alignments[1].get_alignment_length())

    def test_empy(self):
        """Checking empty file."""
        self.assertEqual(0, len(list(ClustalIterator(StringIO("")))))

    def test_write_read(self):
        """Checking write/read."""
        alignments = (list(ClustalIterator(StringIO(aln_example1)))
                      + list(ClustalIterator(StringIO(aln_example2))) * 2)
        handle = StringIO()
        self.assertEqual(3, ClustalWriter(handle).write_file(alignments))
        handle.seek(0)
        for i, a in enumerate(ClustalIterator(handle)):
            self.assertEqual(a.get_alignment_length(), alignments[i].get_alignment_length())

    def test_write_read_single(self):
        """Testing write/read when there is only one sequence."""
        alignment = next(ClustalIterator(StringIO(aln_example1)))
        # Now thae just the first row as a new alignment:
        alignment = alignment[0:1]
        handle = StringIO()
        ClustalWriter(handle).write_file([alignment])
        handle.seek(0)
        for i, a in enumerate(ClustalIterator(handle)):
            self.assertEqual(a.get_alignment_length(), alignment.get_alignment_length())
            self.assertEqual(len(a), 1)

    def test_three(self):
        alignments = list(ClustalIterator(StringIO(aln_example3)))
        self.assertEqual(1, len(alignments))
        self.assertEqual(alignments[0]._version, "2.0.9")

    def test_kalign_header(self):
        """Make sure we can parse the Kalign header."""
        alignments = next(ClustalIterator(StringIO(aln_example4)))
        self.assertEqual(2, len(alignments))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
