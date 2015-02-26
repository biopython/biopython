# Copyright 2006-2014 by Peter Cock.  All rights reserved.
# Revisions copyright 2011 Brandon Invergo. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for Bio.AlignIO.PhylipIO"""

import unittest

from Bio._py3k import StringIO

from Bio.AlignIO.PhylipIO import PhylipIterator, PhylipWriter

phylip_text = """     8    286
V_Harveyi_ --MKNWIKVA VAAIA--LSA A--------- ---------T VQAATEVKVG
B_subtilis MKMKKWTVLV VAALLAVLSA CG-------- ----NGNSSS KEDDNVLHVG
B_subtilis MKKALLALFM VVSIAALAAC GAGNDNQSKD NAKDGDLWAS IKKKGVLTVG
YA80_HAEIN MKKLLFTTAL LTGAIAFSTF ---------- -SHAGEIADR VEKTKTLLVG
FLIY_ECOLI MKLAHLGRQA LMGVMAVALV AG---MSVKS FADEG-LLNK VKERGTLLVG
E_coli_Gln --MKSVLKVS LAALTLAFAV S--------- ---------S HAADKKLVVA
Deinococcu -MKKSLLSLK LSGLLVPSVL ALS------- -LSACSSPSS TLNQGTLKIA
HISJ_E_COL MKKLVLSLSL VLAFSSATAA F--------- ---------- AAIPQNIRIG

           MSGRYFPFTF VKQ--DKLQG FEVDMWDEIG KRNDYKIEYV TANFSGLFGL
           ATGQSYPFAY KEN--GKLTG FDVEVMEAVA KKIDMKLDWK LLEFSGLMGE
           TEGTYEPFTY HDKDTDKLTG YDVEVITEVA KRLGLKVDFK ETQWGSMFAG
           TEGTYAPFTF HDK-SGKLTG FDVEVIRKVA EKLGLKVEFK ETQWDAMYAG
           LEGTYPPFSF QGD-DGKLTG FEVEFAQQLA KHLGVEASLK PTKWDGMLAS
           TDTAFVPFEF KQG--DKYVG FDVDLWAAIA KELKLDYELK PMDFSGIIPA
           MEGTYPPFTS KNE-QGELVG FDVDIAKAVA QKLNLKPEFV LTEWSGILAG
           TDPTYAPFES KNS-QGELVG FDIDLAKELC KRINTQCTFV ENPLDALIPS

           LETGRIDTIS NQITMTDARK AKYLFADPYV VDG-AQITVR KGNDSIQGVE
           LQTGKLDTIS NQVAVTDERK ETYNFTKPYA YAG-TQIVVK KDNTDIKSVD
           LNSKRFDVVA NQVG-KTDRE DKYDFSDKYT TSR-AVVVTK KDNNDIKSEA
           LNAKRFDVIA NQTNPSPERL KKYSFTTPYN YSG-GVIVTK SSDNSIKSFE
           LDSKRIDVVI NQVTISDERK KKYDFSTPYT ISGIQALVKK GNEGTIKTAD
           LQTKNVDLAL AGITITDERK KAIDFSDGYY KSG-LLVMVK ANNNDVKSVK
           LQANKYDVIV NQVGITPERQ NSIGFSQPYA YSRPEIIVAK NNTFNPQSLA
           LKAKKIDAIM SSLSITEKRQ QEIAFTDKLY AADSRLVVAK NSDIQP-TVE

           DLAGKTVAVN LGSNFEQLLR DYDKDGKINI KTYDT--GIE HDVALGRADA
           DLKGKTVAAV LGSNHAKNLE SKDPDKKINI KTYETQEGTL KDVAYGRVDA
           DVKGKTSAQS LTSNYNKLAT N----AGAKV EGVEGMAQAL QMIQQARVDM
           DLKGRKSAQS ATSNWGKDAK A----AGAQI LVVDGLAQSL ELIKQGRAEA
           DLKGKKVGVG LGTNYEEWLR QNV--QGVDV RTYDDDPTKY QDLRVGRIDA
           DLDGKVVAVK SGTGSVDYAK AN--IKTKDL RQFPNIDNAY MELGTNRADA
           DLKGKRVGST LGSNYEKQLI DTG---DIKI VTYPGAPEIL ADLVAGRIDA
           SLKGKRVGVL QGTTQETFGN EHWAPKGIEI VSYQGQDNIY SDLTAGRIDA

           FIMDRLSALE -LIKKT-GLP LQLAGEPFET I-----QNAW PFVDNEKGRK
           YVNSRTVLIA -QIKKT-GLP LKLAGDPIVY E-----QVAF PFAKDDAHDK
           TYNDKLAVLN -YLKTSGNKN VKIAFETGEP Q-----STYF TFRKGS--GE
           TINDKLAVLD -YFKQHPNSG LKIAYDRGDK T-----PTAF AFLQGE--DA
           ILVDRLAALD -LVKKT-NDT LAVTGEAFSR Q-----ESGV ALRKGN--ED
           VLHDTPNILY -FIKTAGNGQ FKAVGDSLEA Q-----QYGI AFPKGS--DE
           AYNDRLVVNY -IINDQ-KLP VRGAGQIGDA A-----PVGI ALKKGN--SA
           AFQDEVAASE GFLKQPVGKD YKFGGPSVKD EKLFGVGTGM GLRKED--NE

           LQAEVNKALA EMRADGTVEK ISVKWFGADI TK----
           LRKKVNKALD ELRKDGTLKK LSEKYFNEDI TVEQKH
           VVDQVNKALK EMKEDGTLSK ISKKWFGEDV SK----
           LITKFNQVLE ALRQDGTLKQ ISIEWFGYDI TQ----
           LLKAVNDAIA EMQKDGTLQA LSEKWFGADV TK----
           LRDKVNGALK TLRENGTYNE IYKKWFGTEP K-----
           LKDQIDKALT EMRSDGTFEK ISQKWFGQDV GQP---
           LREALNKAFA EMRADGTYEK LAKKYFDFDV YGG---
"""

# From here:
# http://atgc.lirmm.fr/phyml/usersguide.html
phylip_text2 = """5 60
Tax1        CCATCTCACGGTCGGTACGATACACCTGCTTTTGGCAG
Tax2        CCATCTCACGGTCAGTAAGATACACCTGCTTTTGGCGG
Tax3        CCATCTCCCGCTCAGTAAGATACCCCTGCTGTTGGCGG
Tax4        TCATCTCATGGTCAATAAGATACTCCTGCTTTTGGCGG
Tax5        CCATCTCACGGTCGGTAAGATACACCTGCTTTTGGCGG

GAAATGGTCAATATTACAAGGT
GAAATGGTCAACATTAAAAGAT
GAAATCGTCAATATTAAAAGGT
GAAATGGTCAATCTTAAAAGGT
GAAATGGTCAATATTAAAAGGT"""

phylip_text3 = """5 60
Tax1        CCATCTCACGGTCGGTACGATACACCTGCTTTTGGCAGGAAATGGTCAATATTACAAGGT
Tax2        CCATCTCACGGTCAGTAAGATACACCTGCTTTTGGCGGGAAATGGTCAACATTAAAAGAT
Tax3        CCATCTCCCGCTCAGTAAGATACCCCTGCTGTTGGCGGGAAATCGTCAATATTAAAAGGT
Tax4        TCATCTCATGGTCAATAAGATACTCCTGCTTTTGGCGGGAAATGGTCAATCTTAAAAGGT
Tax5        CCATCTCACGGTCGGTAAGATACACCTGCTTTTGGCGGGAAATGGTCAATATTAAAAGGT"""

# From here:
# http://evolution.genetics.washington.edu/phylip/doc/sequence.html
# Note the lack of any white space between names 2 and 3 and their seqs.
phylip_text4 = """  5    42
Turkey    AAGCTNGGGC ATTTCAGGGT
Salmo gairAAGCCTTGGC AGTGCAGGGT
H. SapiensACCGGTTGGC CGTTCAGGGT
Chimp     AAACCCTTGC CGTTACGCTT
Gorilla   AAACCCTTGC CGGTACGCTT

GAGCCCGGGC AATACAGGGT AT
GAGCCGTGGC CGGGCACGGT AT
ACAGGTTGGC CGTTCAGGGT AA
AAACCGAGGC CGGGACACTC AT
AAACCATTGC CGGTACGCTT AA"""

# From here:
# http://evolution.genetics.washington.edu/phylip/doc/sequence.html
phylip_text5 = """  5    42
Turkey    AAGCTNGGGC ATTTCAGGGT
GAGCCCGGGC AATACAGGGT AT
Salmo gairAAGCCTTGGC AGTGCAGGGT
GAGCCGTGGC CGGGCACGGT AT
H. SapiensACCGGTTGGC CGTTCAGGGT
ACAGGTTGGC CGTTCAGGGT AA
Chimp     AAACCCTTGC CGTTACGCTT
AAACCGAGGC CGGGACACTC AT
Gorilla   AAACCCTTGC CGGTACGCTT
AAACCATTGC CGGTACGCTT AA"""

phylip_text5a = """  5    42
Turkey    AAGCTNGGGC ATTTCAGGGT GAGCCCGGGC AATACAGGGT AT
Salmo gairAAGCCTTGGC AGTGCAGGGT GAGCCGTGGC CGGGCACGGT AT
H. SapiensACCGGTTGGC CGTTCAGGGT ACAGGTTGGC CGTTCAGGGT AA
Chimp     AAACCCTTGC CGTTACGCTT AAACCGAGGC CGGGACACTC AT
Gorilla   AAACCCTTGC CGGTACGCTT AAACCATTGC CGGTACGCTT AA"""


class TestPhylipIO(unittest.TestCase):

    def test_one(self):
        handle = StringIO(phylip_text)
        ids = []
        for alignment in PhylipIterator(handle):
            for record in alignment:
                ids.append(record.id)
        self.assertEqual(ids, ['V_Harveyi_', 'B_subtilis', 'B_subtilis',
                               'YA80_HAEIN', 'FLIY_ECOLI', 'E_coli_Gln',
                               'Deinococcu', 'HISJ_E_COL'])

        expected = """mkklvlslsl vlafssataa faaipqniri gtdptyapfe sknsqgelvg
        fdidlakelc krintqctfv enpldalips lkakkidaim sslsitekrq qeiaftdkly
        aadsrlvvak nsdiqptves lkgkrvgvlq gttqetfgne hwapkgieiv syqgqdniys
        dltagridaafqdevaaseg flkqpvgkdy kfggpsvkde klfgvgtgmg lrkednelre
        alnkafaemradgtyeklak kyfdfdvygg""".replace(" ", "").replace("\n", "").upper()
        self.assertEqual(str(record.seq).replace("-", ""), expected)

    def test_two_and_three(self):
        handle = StringIO(phylip_text2)
        list2 = list(PhylipIterator(handle))
        handle.close()
        self.assertEqual(len(list2), 1)
        self.assertEqual(len(list2[0]), 5)

        handle = StringIO(phylip_text3)
        list3 = list(PhylipIterator(handle))
        handle.close()
        self.assertEqual(len(list3), 1)
        self.assertEqual(len(list3[0]), 5)

        for i in range(0, 5):
            self.assertEqual(list2[0][i].id, list3[0][i].id)
            self.assertEqual(str(list2[0][i].seq), str(list3[0][i].seq))

    def test_four(self):
        handle = StringIO(phylip_text4)
        list4 = list(PhylipIterator(handle))
        handle.close()
        self.assertEqual(len(list4), 1)
        self.assertEqual(len(list4[0]), 5)

    def test_five(self):
        handle = StringIO(phylip_text5)
        self.assertRaises(ValueError, list, PhylipIterator(handle))
        handle.close()

    def test_five_a(self):
        handle = StringIO(phylip_text5a)
        list5 = list(PhylipIterator(handle))
        handle.close()
        self.assertEqual(len(list5), 1)

    def test_concatenation(self):
        handle = StringIO(phylip_text4 + "\n" + phylip_text4)
        self.assertEqual(len(list(PhylipIterator(handle))), 2)

        handle = StringIO(phylip_text3 + "\n" + phylip_text4 + "\n\n\n" + phylip_text)
        self.assertEqual(len(list(PhylipIterator(handle))), 3)

    def test_write_read(self):
        handle = StringIO(phylip_text5a)
        list5 = list(PhylipIterator(handle))
        handle.close()

        handle = StringIO()
        PhylipWriter(handle).write_file(list5)
        handle.seek(0)
        list6 = list(PhylipIterator(handle))

        self.assertEqual(len(list5), len(list6))
        for a1, a2 in zip(list5, list6):
            self.assertEqual(len(a1), len(a2))
            for r1, r2 in zip(a1, a2):
                self.assertEqual(r1.id, r2.id)
                self.assertEqual(str(r1.seq), str(r2.seq))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
