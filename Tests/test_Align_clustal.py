# Copyright 2006-2014 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.Align.clustal module."""
import unittest

from io import StringIO

from Bio import Align

import numpy


class TestClustalReadingWriting(unittest.TestCase):
    def check_reading_writing(self, path):
        alignments = Align.parse(path, "clustal")
        stream = StringIO()
        n = Align.write(alignments, stream, "clustal")
        self.assertEqual(n, 1)
        alignments = Align.parse(path, "clustal")
        alignment = next(alignments)
        stream.seek(0)
        saved_alignments = Align.parse(stream, "clustal")
        self.assertEqual(saved_alignments.metadata, alignments.metadata)
        saved_alignment = next(saved_alignments)
        with self.assertRaises(StopIteration):
            next(saved_alignments)
        self.assertEqual(len(alignment), len(saved_alignment))
        for i, (sequence, saved_sequence) in enumerate(
            zip(alignment.sequences, saved_alignment.sequences)
        ):
            self.assertEqual(sequence.id, saved_sequence.id)
            self.assertEqual(sequence.seq, saved_sequence.seq)
            self.assertEqual(alignment[i], saved_alignment[i])

    def test_clustalw(self):
        path = "Clustalw/clustalw.aln"
        # includes the sequence length on the right hand side of each line
        with open(path) as stream:
            alignments = Align.parse(stream, "clustal")
            self.assertEqual(alignments.metadata["Program"], "CLUSTAL")
            self.assertEqual(alignments.metadata["Version"], "1.81")
            alignment = next(alignments)
            with self.assertRaises(StopIteration):
                next(alignments)
        self.assertEqual(
            repr(alignment),
            "<Alignment object (2 rows x 601 columns) at 0x%x>" % id(alignment),
        )
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.sequences[0].id, "gi|4959044|gb|AAD34209.1|AF069")
        self.assertEqual(alignment.sequences[1].id, "gi|671626|emb|CAA85685.1|")
        self.assertEqual(
            alignment.sequences[0].seq,
            "MENSDSNDKGSDQSAAQRRSQMDRLDREEAFYQFVNNLSEEDYRLMRDNNLLGTPGESTEEELLRRLQQIKEGPPPQSPDENRAGESSDDVTNSDSIIDWLNSVRQTGNTTRSRQRGNQSWRAVSRTNPNSGDFRFSLEINVNRNNGSQTSENESEPSTRRLSVENMESSSQRQMENSASESASARPSRAERNSTEAVTEVPTTRAQRRARSRSPEHRRTRARAERSMSPLQPTSEIPRRAPTLEQSSENEPEGSSRTRHHVTLRQQISGPELLGRGLFAASGSRNPSQGTSSSDTGSNSESSGSGQRPPTIVLDLQVRRVRPGEYRQRDSIASRTRSRSQAPNNTVTYESERGGFRRTFSRSERAGVRTYVSTIRIPIRRILNTGLSETTSVAIQTMLRQIMTGFGELSYFMYSDSDSEPSASVSSRNVERVESRNGRGSSGGGNSSGSSSSSSPSPSSSGESSESSSKMFEGSSEGGSSGPSRKDGRHRAPVTFDESGSLPFFSLAQFFLLNEDDEDQPRGLTKEQIDNLAMRSFGENDALKTCSVCITEYTEGDKLRKLPCSHEFHVHCIDRWLSENSTCPICRRAVLSSGNRESVV",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "MSPQTETKASVGFKAGVKEYKLTYYTPEYETKDTDILAAFRVTPQPGVPPEEAGAAVAAESSTGTWTTVWTDGLTSLDRYKGRCYHIEPVPGEKDQCICYVAYPLDLFEEGSVTNMFTSIVGNVFGFKALRALRLEDLRIPVAYVKTFQGPPHGIQVERDKLNKYGRPLLGCTIKPKLGLSAKNYGRAVYECLRGGLDFTKDDENVNSQPFMRWRDRFLFCAEAIYKAQAETGEIKGHYLNATAGTCEEMIKRAIFARELGVPIVMHDYLTGGFTANTSLAHYCRDNGLLLHIHRAMHAVIDRQKNHGMHFRVLAKALRLSGGDHIHSGTVVGKLEGERDITLGFVDLLRDDFIEKDRSRGIYFTQDWVSLPGVIPVASGGIHVWHMPALTEIFGDDSVLQFGGGTLGHPWGNAPGAVANRVAVEACVKARNEGRDLAAEGNAIIREACKWSPELAAACEVWKEIKFEFPAMD",
        )
        self.assertEqual(
            alignment[0],
            "MENSDSNDKGSDQSAAQRRSQMDRLDREEAFYQFVNNLSEEDYRLMRDNNLLGTPGESTEEELLRRLQQIKEGPPPQSPDENRAGESSDDVTNSDSIIDWLNSVRQTGNTTRSRQRGNQSWRAVSRTNPNSGDFRFSLEINVNRNNGSQTSENESEPSTRRLSVENMESSSQRQMENSASESASARPSRAERNSTEAVTEVPTTRAQRRARSRSPEHRRTRARAERSMSPLQPTSEIPRRAPTLEQSSENEPEGSSRTRHHVTLRQQISGPELLGRGLFAASGSRNPSQGTSSSDTGSNSESSGSGQRPPTIVLDLQVRRVRPGEYRQRDSIASRTRSRSQAPNNTVTYESERGGFRRTFSRSERAGVRTYVSTIRIPIRRILNTGLSETTSVAIQTMLRQIMTGFGELSYFMYSDSDSEPSASVSSRNVERVESRNGRGSSGGGNSSGSSSSSSPSPSSSGESSESSSKMFEGSSEGGSSGPSRKDGRHRAPVTFDESGSLPFFSLAQFFLLNEDDEDQPRGLTKEQIDNLAMRSFGENDALKTCSVCITEYTEGDKLRKLPCSHEFHVHCIDRWLSE-NSTCPICRRAVLSSGNRESVV",
        )
        self.assertEqual(
            alignment[1],
            "---------MSPQTETKASVGFKAGVKEYKLTYYTPEYETKDTDILAAFRVTPQPG-----------------VPPEEAGAAVAAESSTGT---------WTTVWTDGLTSLDRYKG-----RCYHIEPVPG-------------------EKDQCICYVAYPLDLFEEGSVTNMFTSIVGNVFGFKALRALRLEDLRIPVAYVKTFQGPPHGIQVERDKLNKYGRPLLGCTIKPKLGLSAKNYGRAVYECLRGGLDFTKDDENVNSQPFMRWRDRFLFCAEAIYKAQAETGEIKGHYLNATAG-----------------------TCEEMIKRAIFARELGVPIVMHDYLTGGFTANTSLAHYCRDNGLLLHIHRAMHAVIDRQKNHGMHFRVLAKALRLSGGDHIHSGTVVGKLEGERDITLGFVDLLRDDFIEKDRSRGIYFTQDWVSLPGVIPVASG-----------------------------GIHVWHMPALTEIFGDDSVLQFGGGTLGHPWGNAPGAVANRVA-----------VEACVKARNEG---RDLAAEGNAIIREACKWSPELAAACEVWKEIKFEFPAMD---",
        )
        self.assertEqual(
            format(alignment, "clustal"),
            """\
gi|4959044|gb|AAD34209.1|AF069      MENSDSNDKGSDQSAAQRRSQMDRLDREEAFYQFVNNLSEEDYRLMRDNN
gi|671626|emb|CAA85685.1|           ---------MSPQTETKASVGFKAGVKEYKLTYYTPEYETKDTDILAAFR
                                              * *: ::    :.   :*  :  :. : . :*  ::   .

gi|4959044|gb|AAD34209.1|AF069      LLGTPGESTEEELLRRLQQIKEGPPPQSPDENRAGESSDDVTNSDSIIDW
gi|671626|emb|CAA85685.1|           VTPQPG-----------------VPPEEAGAAVAAESSTGT---------
                                    :   **                  **:...   *.*** ..         

gi|4959044|gb|AAD34209.1|AF069      LNSVRQTGNTTRSRQRGNQSWRAVSRTNPNSGDFRFSLEINVNRNNGSQT
gi|671626|emb|CAA85685.1|           WTTVWTDGLTSLDRYKG-----RCYHIEPVPG------------------
                                     .:*   * *: .* :*        : :* .*                  

gi|4959044|gb|AAD34209.1|AF069      SENESEPSTRRLSVENMESSSQRQMENSASESASARPSRAERNSTEAVTE
gi|671626|emb|CAA85685.1|           -EKDQCICYVAYPLDLFEEGSVTNMFTSIVGNVFGFKALRALRLEDLRIP
                                     *::.  .    .:: :*..*  :* .*   .. .  :    .  :    

gi|4959044|gb|AAD34209.1|AF069      VPTTRAQRRARSRSPEHRRTRARAERSMSPLQPTSEIPRRAPTLEQSSEN
gi|671626|emb|CAA85685.1|           VAYVKTFQGPPHGIQVERDKLNKYGRPLLGCTIKPKLGLSAKNYGRAVYE
                                    *. .:: : .      .* .  :  *.:     ..::   * .  ::  :

gi|4959044|gb|AAD34209.1|AF069      EPEGSSRTRHHVTLRQQISGPELLGRGLFAASGSRNPSQGTSSSDTGSNS
gi|671626|emb|CAA85685.1|           CLRGGLDFTKDDENVNSQPFMRWRDRFLFCAEAIYKAQAETGEIKGHYLN
                                      .*.    :.    :. .  .  .* **.*..  :..  *.. .    .

gi|4959044|gb|AAD34209.1|AF069      ESSGSGQRPPTIVLDLQVRRVRPGEYRQRDSIASRTRSRSQAPNNTVTYE
gi|671626|emb|CAA85685.1|           ATAG-----------------------TCEEMIKRAIFARELGVPIVMHD
                                     ::*                         :.: .*:    :     * ::

gi|4959044|gb|AAD34209.1|AF069      SERGGFRRTFSRSERAGVRTYVSTIRIPIRRILNTGLSETTSVAIQTMLR
gi|671626|emb|CAA85685.1|           YLTGGFTANTSLAHYCRDNGLLLHIHRAMHAVIDRQKNHGMHFRVLAKAL
                                       ***  . * :. .  .  :  *: .:: :::   ..   . : :   

gi|4959044|gb|AAD34209.1|AF069      QIMTGFGELSYFMYSDSDSEPSASVSSRNVERVESRNGRGSSGGGNSSGS
gi|671626|emb|CAA85685.1|           RLSGGDHIHSGTVVGKLEGERDITLGFVDLLRDDFIEKDRSRGIYFTQDW
                                    ::  *    *  : .. :.* . ::.  :: * :  :   * *   :.. 

gi|4959044|gb|AAD34209.1|AF069      SSSSSPSPSSSGESSESSSKMFEGSSEGGSSGPSRKDGRHRAPVTFDESG
gi|671626|emb|CAA85685.1|           VSLPGVIPVASG-----------------------------GIHVWHMPA
                                     * ..  * :**                             .  .:. ..

gi|4959044|gb|AAD34209.1|AF069      SLPFFSLAQFFLLNEDDEDQPRGLTKEQIDNLAMRSFGENDALKTCSVCI
gi|671626|emb|CAA85685.1|           LTEIFGDDSVLQFGGGTLGHPWGNAPGAVANRVA-----------VEACV
                                       :*.  ..: :. .  .:* * :   : * .             ..*:

gi|4959044|gb|AAD34209.1|AF069      TEYTEGDKLRKLPCSHEFHVHCIDRWLSE-NSTCPICRRAVLSSGNRESV
gi|671626|emb|CAA85685.1|           KARNEG---RDLAAEGNAIIREACKWSPELAAACEVWKEIKFEFPAMD--
                                    .  .**   *.*... :  ::   :* .*  ::* : :.  :.    :  

gi|4959044|gb|AAD34209.1|AF069      V
gi|671626|emb|CAA85685.1|           -
                                     


""",
        )
        self.check_reading_writing(path)

    def test_msaprobs(self):
        path = "Clustalw/msaprobs.aln"
        # This example was obtained from
        # http://virgil.ruc.dk/kurser/Sekvens/Treedraw.htm
        with open(path) as stream:
            alignments = Align.parse(stream, "clustal")
            self.assertEqual(alignments.metadata["Program"], "MSAPROBS")
            self.assertEqual(alignments.metadata["Version"], "0.9.7")
            alignment = next(alignments)
            with self.assertRaises(StopIteration):
                next(alignments)
        self.assertEqual(
            repr(alignment),
            "<Alignment object (8 rows x 298 columns) at 0x%x>" % id(alignment),
        )
        self.assertEqual(len(alignment), 8)
        self.assertEqual(alignment.shape, (8, 298))
        self.assertEqual(alignment.sequences[0].id, "V_Harveyi_PATH")
        self.assertEqual(alignment.sequences[1].id, "B_subtilis_YXEM")
        self.assertEqual(alignment.sequences[2].id, "FLIY_ECOLI")
        self.assertEqual(alignment.sequences[3].id, "Deinococcus_radiodurans")
        self.assertEqual(alignment.sequences[4].id, "B_subtilis_GlnH_homo_YCKK")
        self.assertEqual(alignment.sequences[5].id, "YA80_HAEIN")
        self.assertEqual(alignment.sequences[6].id, "E_coli_GlnH")
        self.assertEqual(alignment.sequences[7].id, "HISJ_E_COLI")
        self.assertEqual(
            alignment.sequences[0].seq,
            "MKNWIKVAVAAIALSAATVQAATEVKVGMSGRYFPFTFVKQDKLQGFEVDMWDEIGKRNDYKIEYVTANFSGLFGLLETGRIDTISNQITMTDARKAKYLFADPYVVDGAQITVRKGNDSIQGVEDLAGKTVAVNLGSNFEQLLRDYDKDGKINIKTYDTGIEHDVALGRADAFIMDRLSALELIKKTGLPLQLAGEPFETIQNAWPFVDNEKGRKLQAEVNKALAEMRADGTVEKISVKWFGADITK",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "MKMKKWTVLVVAALLAVLSACGNGNSSSKEDDNVLHVGATGQSYPFAYKENGKLTGFDVEVMEAVAKKIDMKLDWKLLEFSGLMGELQTGKLDTISNQVAVTDERKETYNFTKPYAYAGTQIVVKKDNTDIKSVDDLKGKTVAAVLGSNHAKNLESKDPDKKINIKTYETQEGTLKDVAYGRVDAYVNSRTVLIAQIKKTGLPLKLAGDPIVYEQVAFPFAKDDAHDKLRKKVNKALDELRKDGTLKKLSEKYFNEDITVEQKH",
        )
        self.assertEqual(
            alignment.sequences[2].seq,
            "MKLAHLGRQALMGVMAVALVAGMSVKSFADEGLLNKVKERGTLLVGLEGTYPPFSFQGDDGKLTGFEVEFAQQLAKHLGVEASLKPTKWDGMLASLDSKRIDVVINQVTISDERKKKYDFSTPYTISGIQALVKKGNEGTIKTADDLKGKKVGVGLGTNYEEWLRQNVQGVDVRTYDDDPTKYQDLRVGRIDAILVDRLAALDLVKKTNDTLAVTGEAFSRQESGVALRKGNEDLLKAVNDAIAEMQKDGTLQALSEKWFGADVTK",
        )
        self.assertEqual(
            alignment.sequences[3].seq,
            "MKKSLLSLKLSGLLVPSVLALSLSACSSPSSTLNQGTLKIAMEGTYPPFTSKNEQGELVGFDVDIAKAVAQKLNLKPEFVLTEWSGILAGLQANKYDVIVNQVGITPERQNSIGFSQPYAYSRPEIIVAKNNTFNPQSLADLKGKRVGSTLGSNYEKQLIDTGDIKIVTYPGAPEILADLVAGRIDAAYNDRLVVNYIINDQKLPVRGAGQIGDAAPVGIALKKGNSALKDQIDKALTEMRSDGTFEKISQKWFGQDVGQP",
        )
        self.assertEqual(
            alignment.sequences[4].seq,
            "MKKALLALFMVVSIAALAACGAGNDNQSKDNAKDGDLWASIKKKGVLTVGTEGTYEPFTYHDKDTDKLTGYDVEVITEVAKRLGLKVDFKETQWGSMFAGLNSKRFDVVANQVGKTDREDKYDFSDKYTTSRAVVVTKKDNNDIKSEADVKGKTSAQSLTSNYNKLATNAGAKVEGVEGMAQALQMIQQARVDMTYNDKLAVLNYLKTSGNKNVKIAFETGEPQSTYFTFRKGSGEVVDQVNKALKEMKEDGTLSKISKKWFGEDVSK",
        )
        self.assertEqual(
            alignment.sequences[5].seq,
            "MKKLLFTTALLTGAIAFSTFSHAGEIADRVEKTKTLLVGTEGTYAPFTFHDKSGKLTGFDVEVIRKVAEKLGLKVEFKETQWDAMYAGLNAKRFDVIANQTNPSPERLKKYSFTTPYNYSGGVIVTKSSDNSIKSFEDLKGRKSAQSATSNWGKDAKAAGAQILVVDGLAQSLELIKQGRAEATINDKLAVLDYFKQHPNSGLKIAYDRGDKTPTAFAFLQGEDALITKFNQVLEALRQDGTLKQISIEWFGYDITQ",
        )
        self.assertEqual(
            alignment.sequences[6].seq,
            "MKSVLKVSLAALTLAFAVSSHAADKKLVVATDTAFVPFEFKQGDKYVGFDVDLWAAIAKELKLDYELKPMDFSGIIPALQTKNVDLALAGITITDERKKAIDFSDGYYKSGLLVMVKANNNDVKSVKDLDGKVVAVKSGTGSVDYAKANIKTKDLRQFPNIDNAYMELGTNRADAVLHDTPNILYFIKTAGNGQFKAVGDSLEAQQYGIAFPKGSDELRDKVNGALKTLRENGTYNEIYKKWFGTEPK",
        )
        self.assertEqual(
            alignment.sequences[7].seq,
            "MKKLVLSLSLVLAFSSATAAFAAIPQNIRIGTDPTYAPFESKNSQGELVGFDIDLAKELCKRINTQCTFVENPLDALIPSLKAKKIDAIMSSLSITEKRQQEIAFTDKLYAADSRLVVAKNSDIQPTVESLKGKRVGVLQGTTQETFGNEHWAPKGIEIVSYQGQDNIYSDLTAGRIDAAFQDEVAASEGFLKQPVGKDYKFGGPSVKDEKLFGVGTGMGLRKEDNELREALNKAFAEMRADGTYEKLAKKYFDFDVYGG",
        )
        self.assertEqual(
            alignment[0],
            "MKNW--------IKV----AVAAI-A--LSAA-------------------TVQAATEVKVGMSGRYFPFTFVK--QDKLQGFEVDMWDEIGKRNDYKIEYVTANFSGLFGLLETGRIDTISNQITMTDARKAKYLFADPYVVDGAQITVRK-GNDSIQGVEDLAGKTVAVNLGSNFEQLLRDYDKDGKINIKTYDT--GIEHDVALGRADAFIMDRLSALE-LIKKTG-LPLQLAGEPFE-----TIQNAWPFVDNEKGRKLQAEVNKALAEMRADGTVEKISVKWFGADITK----",
        )
        self.assertEqual(
            alignment[1],
            "MKMKKW------TVL----VVAALLA-VLSACGN------------G-NSSSKEDDNVLHVGATGQSYPFAYKE--NGKLTGFDVEVMEAVAKKIDMKLDWKLLEFSGLMGELQTGKLDTISNQVAVTDERKETYNFTKPYAYAGTQIVVKK-DNTDIKSVDDLKGKTVAAVLGSNHAKNLESKDPDKKINIKTYETQEGTLKDVAYGRVDAYVNSRTVLIA-QIKKTG-LPLKLAGDPIV-----YEQVAFPFAKDDAHDKLRKKVNKALDELRKDGTLKKLSEKYFNEDITVEQKH",
        )
        self.assertEqual(
            alignment[2],
            "MKLAHLGRQALMGVM----AVALVAG--MSVKSF---------ADEG-LLNKVKERGTLLVGLEGTYPPFSFQGD-DGKLTGFEVEFAQQLAKHLGVEASLKPTKWDGMLASLDSKRIDVVINQVTISDERKKKYDFSTPYTISGIQALVKKGNEGTIKTADDLKGKKVGVGLGTNYEEWLRQN--VQGVDVRTYDDDPTKYQDLRVGRIDAILVDRLAALD-LVKKTN-DTLAVTGEAFS-----RQESGVALRK--GNEDLLKAVNDAIAEMQKDGTLQALSEKWFGADVTK----",
        )
        self.assertEqual(
            alignment[3],
            "MKKSLL------SLKLSGLLVPSVLALSLSACSS---------------PSSTLNQGTLKIAMEGTYPPFTSKNE-QGELVGFDVDIAKAVAQKLNLKPEFVLTEWSGILAGLQANKYDVIVNQVGITPERQNSIGFSQPYAYSRPEIIVAKNNTFNPQSLADLKGKRVGSTLGSNYEKQLI-D--TGDIKIVTYPGAPEILADLVAGRIDAAYNDRLVVNY-IIND-QKLPVRGAGQIGD-----AAPVGIALKK--GNSALKDQIDKALTEMRSDGTFEKISQKWFGQDVGQ---P",
        )
        self.assertEqual(
            alignment[4],
            "MKKALL------ALF----MVVSIAA--LAACGAGNDNQSKDNAKDGDLWASIKKKGVLTVGTEGTYEPFTYHDKDTDKLTGYDVEVITEVAKRLGLKVDFKETQWGSMFAGLNSKRFDVVANQVGKTD-REDKYDFSDKYTTSRAVVVTKK-DNNDIKSEADVKGKTSAQSLTSNYNKLAT-N--A-GAKVEGVEGMAQALQMIQQARVDMTYNDKLAVLN-YLKTSGNKNVKIAFETGE-----PQSTYFTFRK--GSGEVVDQVNKALKEMKEDGTLSKISKKWFGEDVSK----",
        )
        self.assertEqual(
            alignment[5],
            "MKKLLF------TTA----LLTGAIA--FSTFS-----------HAGEIADRVEKTKTLLVGTEGTYAPFTFHDK-SGKLTGFDVEVIRKVAEKLGLKVEFKETQWDAMYAGLNAKRFDVIANQTNPSPERLKKYSFTTPYNYSGGVIVTKS-SDNSIKSFEDLKGRKSAQSATSNWGKDAK-A--A-GAQILVVDGLAQSLELIKQGRAEATINDKLAVLD-YFKQHPNSGLKIAYDRGD-----KTPTAFAFLQ--GEDALITKFNQVLEALRQDGTLKQISIEWFGYDITQ----",
        )
        self.assertEqual(
            alignment[6],
            "MKSVL-------KVS----LAALTLA--FAVSSH---------A----------ADKKLVVATDTAFVPFEFKQ--GDKYVGFDVDLWAAIAKELKLDYELKPMDFSGIIPALQTKNVDLALAGITITDERKKAIDFSDGYYKSGLLVMVKAN-NNDVKSVKDLDGKVVAVKSGTGSVDYAKAN--IKTKDLRQFPNIDNAYMELGTNRADAVLHDTPNILY-FIKTAGNGQFKAVGDSLE-----AQQYGIAFPK--GSDELRDKVNGALKTLRENGTYNEIYKKWFGTEP-K----",
        )
        self.assertEqual(
            alignment[7],
            "MKKLVL------SLS----LV---LA--FSSATA---------------A-FAAIPQNIRIGTDPTYAPFESKNS-QGELVGFDIDLAKELCKRINTQCTFVENPLDALIPSLKAKKIDAIMSSLSITEKRQQEIAFTDKLYAADSRLVVAK-NSDIQPTVESLKGKRVGVLQGTTQETFGNEHWAPKGIEIVSYQGQDNIYSDLTAGRIDAAFQDEVAASEGFLKQPVGKDYKFGGPSVKDEKLFGVGTGMGLRK--EDNELREALNKAFAEMRADGTYEKLAKKYFDFDVYG---G",
        )
        self.assertEqual(
            alignment.column_annotations["clustal_consensus"],
            "**                       .  ::             *. *:          : :.      **    .  .:  *::::.   : :.   .        ..:   *.: . *        :  *     *:           .  ..        .: *:  .    :               .:            :   * :    .        .:                           : .::    :   .: .:  :: :** . :  ::*. :       ",
        )
        self.assertEqual(
            format(alignment, "clustal"),
            """\
V_Harveyi_PATH                      MKNW--------IKV----AVAAI-A--LSAA------------------
B_subtilis_YXEM                     MKMKKW------TVL----VVAALLA-VLSACGN------------G-NS
FLIY_ECOLI                          MKLAHLGRQALMGVM----AVALVAG--MSVKSF---------ADEG-LL
Deinococcus_radiodurans             MKKSLL------SLKLSGLLVPSVLALSLSACSS---------------P
B_subtilis_GlnH_homo_YCKK           MKKALL------ALF----MVVSIAA--LAACGAGNDNQSKDNAKDGDLW
YA80_HAEIN                          MKKLLF------TTA----LLTGAIA--FSTFS-----------HAGEIA
E_coli_GlnH                         MKSVL-------KVS----LAALTLA--FAVSSH---------A------
HISJ_E_COLI                         MKKLVL------SLS----LV---LA--FSSATA---------------A
                                    **                       .  ::             *. *:  

V_Harveyi_PATH                      -TVQAATEVKVGMSGRYFPFTFVK--QDKLQGFEVDMWDEIGKRNDYKIE
B_subtilis_YXEM                     SSKEDDNVLHVGATGQSYPFAYKE--NGKLTGFDVEVMEAVAKKIDMKLD
FLIY_ECOLI                          NKVKERGTLLVGLEGTYPPFSFQGD-DGKLTGFEVEFAQQLAKHLGVEAS
Deinococcus_radiodurans             SSTLNQGTLKIAMEGTYPPFTSKNE-QGELVGFDVDIAKAVAQKLNLKPE
B_subtilis_GlnH_homo_YCKK           ASIKKKGVLTVGTEGTYEPFTYHDKDTDKLTGYDVEVITEVAKRLGLKVD
YA80_HAEIN                          DRVEKTKTLLVGTEGTYAPFTFHDK-SGKLTGFDVEVIRKVAEKLGLKVE
E_coli_GlnH                         ----ADKKLVVATDTAFVPFEFKQ--GDKYVGFDVDLWAAIAKELKLDYE
HISJ_E_COLI                         -FAAIPQNIRIGTDPTYAPFESKNS-QGELVGFDIDLAKELCKRINTQCT
                                            : :.      **    .  .:  *::::.   : :.   .  

V_Harveyi_PATH                      YVTANFSGLFGLLETGRIDTISNQITMTDARKAKYLFADPYVVDGAQITV
B_subtilis_YXEM                     WKLLEFSGLMGELQTGKLDTISNQVAVTDERKETYNFTKPYAYAGTQIVV
FLIY_ECOLI                          LKPTKWDGMLASLDSKRIDVVINQVTISDERKKKYDFSTPYTISGIQALV
Deinococcus_radiodurans             FVLTEWSGILAGLQANKYDVIVNQVGITPERQNSIGFSQPYAYSRPEIIV
B_subtilis_GlnH_homo_YCKK           FKETQWGSMFAGLNSKRFDVVANQVGKTD-REDKYDFSDKYTTSRAVVVT
YA80_HAEIN                          FKETQWDAMYAGLNAKRFDVIANQTNPSPERLKKYSFTTPYNYSGGVIVT
E_coli_GlnH                         LKPMDFSGIIPALQTKNVDLALAGITITDERKKAIDFSDGYYKSGLLVMV
HISJ_E_COLI                         FVENPLDALIPSLKAKKIDAIMSSLSITEKRQQEIAFTDKLYAADSRLVV
                                          ..:   *.: . *        :  *     *:           .

V_Harveyi_PATH                      RK-GNDSIQGVEDLAGKTVAVNLGSNFEQLLRDYDKDGKINIKTYDT--G
B_subtilis_YXEM                     KK-DNTDIKSVDDLKGKTVAAVLGSNHAKNLESKDPDKKINIKTYETQEG
FLIY_ECOLI                          KKGNEGTIKTADDLKGKKVGVGLGTNYEEWLRQN--VQGVDVRTYDDDPT
Deinococcus_radiodurans             AKNNTFNPQSLADLKGKRVGSTLGSNYEKQLI-D--TGDIKIVTYPGAPE
B_subtilis_GlnH_homo_YCKK           KK-DNNDIKSEADVKGKTSAQSLTSNYNKLAT-N--A-GAKVEGVEGMAQ
YA80_HAEIN                          KS-SDNSIKSFEDLKGRKSAQSATSNWGKDAK-A--A-GAQILVVDGLAQ
E_coli_GlnH                         KAN-NNDVKSVKDLDGKVVAVKSGTGSVDYAKAN--IKTKDLRQFPNIDN
HISJ_E_COLI                         AK-NSDIQPTVESLKGKRVGVLQGTTQETFGNEHWAPKGIEIVSYQGQDN
                                      ..        .: *:  .    :               .:        

V_Harveyi_PATH                      IEHDVALGRADAFIMDRLSALE-LIKKTG-LPLQLAGEPFE-----TIQN
B_subtilis_YXEM                     TLKDVAYGRVDAYVNSRTVLIA-QIKKTG-LPLKLAGDPIV-----YEQV
FLIY_ECOLI                          KYQDLRVGRIDAILVDRLAALD-LVKKTN-DTLAVTGEAFS-----RQES
Deinococcus_radiodurans             ILADLVAGRIDAAYNDRLVVNY-IIND-QKLPVRGAGQIGD-----AAPV
B_subtilis_GlnH_homo_YCKK           ALQMIQQARVDMTYNDKLAVLN-YLKTSGNKNVKIAFETGE-----PQST
YA80_HAEIN                          SLELIKQGRAEATINDKLAVLD-YFKQHPNSGLKIAYDRGD-----KTPT
E_coli_GlnH                         AYMELGTNRADAVLHDTPNILY-FIKTAGNGQFKAVGDSLE-----AQQY
HISJ_E_COLI                         IYSDLTAGRIDAAFQDEVAASEGFLKQPVGKDYKFGGPSVKDEKLFGVGT
                                        :   * :    .        .:                        

V_Harveyi_PATH                      AWPFVDNEKGRKLQAEVNKALAEMRADGTVEKISVKWFGADITK----
B_subtilis_YXEM                     AFPFAKDDAHDKLRKKVNKALDELRKDGTLKKLSEKYFNEDITVEQKH
FLIY_ECOLI                          GVALRK--GNEDLLKAVNDAIAEMQKDGTLQALSEKWFGADVTK----
Deinococcus_radiodurans             GIALKK--GNSALKDQIDKALTEMRSDGTFEKISQKWFGQDVGQ---P
B_subtilis_GlnH_homo_YCKK           YFTFRK--GSGEVVDQVNKALKEMKEDGTLSKISKKWFGEDVSK----
YA80_HAEIN                          AFAFLQ--GEDALITKFNQVLEALRQDGTLKQISIEWFGYDITQ----
E_coli_GlnH                         GIAFPK--GSDELRDKVNGALKTLRENGTYNEIYKKWFGTEP-K----
HISJ_E_COLI                         GMGLRK--EDNELREALNKAFAEMRADGTYEKLAKKYFDFDVYG---G
                                       : .::    :   .: .:  :: :** . :  ::*. :       


""",
        )
        self.check_reading_writing(path)

    def test_muscle(self):
        path = "Clustalw/muscle.aln"
        # includes the sequence length on the right hand side of each line
        with open(path) as stream:
            alignments = Align.parse(stream, "clustal")
            self.assertEqual(alignments.metadata["Program"], "MUSCLE")
            self.assertEqual(alignments.metadata["Version"], "3.8")
            alignment = next(alignments)
            with self.assertRaises(StopIteration):
                next(alignments)
        self.assertEqual(
            repr(alignment),
            "<Alignment object (3 rows x 687 columns) at 0x%x>" % id(alignment),
        )
        self.assertEqual(len(alignment), 3)
        self.assertEqual(alignment.sequences[0].id, "Test1seq")
        self.assertEqual(alignment.sequences[1].id, "AT3G20900.1-SEQ")
        self.assertEqual(alignment.sequences[2].id, "AT3G20900.1-CDS")
        self.assertEqual(
            alignment.sequences[0].seq,
            "AGTTACAATAACTGACGAAGCTAAGTAGGCTACTAATTAACGTCATCAACCTAATACATAGCACTTAGAAAAAAGTGAAGTAAGAAAATATAAAATAATAAAAGGGTGGGTTATCAATTGATAGTGTAAATCATCGTATTCCGGTGATATACCCTACCACAAAAACTCAAACCGACTTGATTCAAATCATCTCAATAAATTAGCGCCAAAATAATGAAAAAAATAATAACAAACAAAAACAAACCAAAATAAGAAAAAACATTACGCAAAACATAATAATTTACTCTTCGTTATTGTATTAACAAATCAAAGAGCTGAATTTTGATCACCTGCTAATACTACTTTCTGTATTGATCCTATATCAACGTAAACAAAGATACTAATAATTAACTAAAAGTACGTTCATCGATCGTGTTCGTTGACGAAGAAGAGCTCTATCTCCGGCGGAGCAAAGAAAACGATCTGTCTCCGTCGTAACACACGGTCGCTAGAGAAACTTTGCTTCTTCGGCGCCGGTGGACACGTCAGCATCTCCGGTATCCTAGACTTCTTGGCTTTCGGGGTACAACAACCGCGTGGTGACGTCAGCACCGCTGCTGGGGATGGAGAGGGAACAGAGTT",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "ATGAACAAAGTAGCGAGGAAGAACAAAACATCAGGTGAACAAAAAAAAAACTCAATCCACATCAAAGTTACAATAACTGACGAAGCTAAGTAGGCTAGAAATTAAAGTCATCAACCTAATACATAGCACTTAGAAAAAAGTGAAGCAAGAAAATATAAAATAATAAAAGGGTGGGTTATCAATTGATAGTGTAAATCATAGTTGATTTTTGATATACCCTACCACAAAAACTCAAACCGACTTGATTCAAATCATCTCAAAAAACAAGCGCCAAAATAATGAAAAAAATAATAACAAAAACAAACAAACCAAAATAAGAAAAAACATTACGCAAAACATAATAATTTACTCTTCGTTATTGTATTAACAAATCAAAGAGATGAATTTTGATCACCTGCTAATACTACTTTCTGTATTGATCCTATATCAAAAAAAAAAAAGATACTAATAATTAACTAAAAGTACGTTCATCGATCGTGTGCGTTGACGAAGAAGAGCTCTATCTCCGGCGGAGCAAAGAAAACGATCTGTCTCCGTCGTAACACACAGTTTTTCGAGACCCTTTGCTTCTTCGGCGCCGGTGGACACGTCAGCATCTCCGGTATCCTAGACTTCTTGGCTTTCGGGGTACAACAACCGCCTGGTGACGTCAGCACCGCTGCTGGGGATGGAGAGGGAACAGAGTAG",
        )
        self.assertEqual(
            alignment.sequences[2].seq,
            "ATGAACAAAGTAGCGAGGAAGAACAAAACATCAGCAAAGAAAACGATCTGTCTCCGTCGTAACACACAGTTTTTCGAGACCCTTTGCTTCTTCGGCGCCGGTGGACACGTCAGCATCTCCGGTATCCTAGACTTCTTGGCTTTCGGGGTACAACAACCGCCTGGTGACGTCAGCACCGCTGCTGGGGATGGAGAGGGAACAGAGTAG",
        )
        self.assertEqual(
            alignment[0],
            "-----------------------------------------------------------------AGTTACAATAACTGACGAAGCTAAGTAGGCTACTAATTAACGTCATCAACCTAATACATAGCACTTAGAAAAAAGTGAAGTAAGAAAATATAAAATAATAAAAGGGTGGGTTATCAATTGATAGTGTAAATCATCGTATTCCGGTGATATACCCTACCACAAAAACTCAAACCGACTTGATTCAAATCATCTCAATAAATTAGCGCCAAAATAATGAAAAAAATAATAACAAACAAAAACAAACCAAAATAAGAAAAAACATTACGCAAAACATAATAATTTACTCTTCGTTATTGTATTAACAAATCAAAGAGCTGAATTTTGATCACCTGCTAATACTACTTTCTGTATTGATCCTATATCAACGTAAACAAAGATACTAATAATTAACTAAAAGTACGTTCATCGATCGTGTTCGTTGACGAAGAAGAGCTCTATCTCCGGCGGAGCAAAGAAAACGATCTGTCTCCGTCGTAACACACGGTCGCTAGAGAAACTTTGCTTCTTCGGCGCCGGTGGACACGTCAGCATCTCCGGTATCCTAGACTTCTTGGCTTTCGGGGTACAACAACCGCGTGGTGACGTCAGCACCGCTGCTGGGGATGGAGAGGGAACAGAGTT-",
        )
        self.assertEqual(
            alignment[1],
            "ATGAACAAAGTAGCGAGGAAGAACAAAACATCAGGTGAACAAAAAAAAAACTCAATCCACATCAAAGTTACAATAACTGACGAAGCTAAGTAGGCTAGAAATTAAAGTCATCAACCTAATACATAGCACTTAGAAAAAAGTGAAGCAAGAAAATATAAAATAATAAAAGGGTGGGTTATCAATTGATAGTGTAAATCATAGTTGATTTTTGATATACCCTACCACAAAAACTCAAACCGACTTGATTCAAATCATCTCAAAAAACAAGCGCCAAAATAATGAAAAAAATAATAACAAAAACAAACAAACCAAAATAAGAAAAAACATTACGCAAAACATAATAATTTACTCTTCGTTATTGTATTAACAAATCAAAGAGATGAATTTTGATCACCTGCTAATACTACTTTCTGTATTGATCCTATATCAAAAAAAAAAAAGATACTAATAATTAACTAAAAGTACGTTCATCGATCGTGTGCGTTGACGAAGAAGAGCTCTATCTCCGGCGGAGCAAAGAAAACGATCTGTCTCCGTCGTAACACACAGTTTTTCGAGACCCTTTGCTTCTTCGGCGCCGGTGGACACGTCAGCATCTCCGGTATCCTAGACTTCTTGGCTTTCGGGGTACAACAACCGCCTGGTGACGTCAGCACCGCTGCTGGGGATGGAGAGGGAACAGAGTAG",
        )
        self.assertEqual(
            alignment[2],
            "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ATGAACAAAGTAGCGAGGAAGAA------------------------------CAAAACATC----------------------------------------------------------------------------------------------------------------------------------------------------------------------------AGCAAAGAAAACGATCTGTCTCCGTCGTAACACACAGTTTTTCGAGACCCTTTGCTTCTTCGGCGCCGGTGGACACGTCAGCATCTCCGGTATCCTAGACTTCTTGGCTTTCGGGGTACAACAACCGCCTGGTGACGTCAGCACCGCTGCTGGGGATGGAGAGGGAACAGAGTAG",
        )
        self.assertEqual(
            alignment.column_annotations["clustal_consensus"],
            "                                                                                                                                                                                                                                                                                      ***** *** **   *  ** *                               ********                                                                                                                                                                             *********************************** **   * ****  ******************************************************************************* ********************************************  ",
        )
        self.assertEqual(
            format(alignment, "clustal"),
            """\
Test1seq                            --------------------------------------------------
AT3G20900.1-SEQ                     ATGAACAAAGTAGCGAGGAAGAACAAAACATCAGGTGAACAAAAAAAAAA
AT3G20900.1-CDS                     --------------------------------------------------
                                                                                      

Test1seq                            ---------------AGTTACAATAACTGACGAAGCTAAGTAGGCTACTA
AT3G20900.1-SEQ                     CTCAATCCACATCAAAGTTACAATAACTGACGAAGCTAAGTAGGCTAGAA
AT3G20900.1-CDS                     --------------------------------------------------
                                                                                      

Test1seq                            ATTAACGTCATCAACCTAATACATAGCACTTAGAAAAAAGTGAAGTAAGA
AT3G20900.1-SEQ                     ATTAAAGTCATCAACCTAATACATAGCACTTAGAAAAAAGTGAAGCAAGA
AT3G20900.1-CDS                     --------------------------------------------------
                                                                                      

Test1seq                            AAATATAAAATAATAAAAGGGTGGGTTATCAATTGATAGTGTAAATCATC
AT3G20900.1-SEQ                     AAATATAAAATAATAAAAGGGTGGGTTATCAATTGATAGTGTAAATCATA
AT3G20900.1-CDS                     --------------------------------------------------
                                                                                      

Test1seq                            GTATTCCGGTGATATACCCTACCACAAAAACTCAAACCGACTTGATTCAA
AT3G20900.1-SEQ                     GTTGATTTTTGATATACCCTACCACAAAAACTCAAACCGACTTGATTCAA
AT3G20900.1-CDS                     --------------------------------------------------
                                                                                      

Test1seq                            ATCATCTCAATAAATTAGCGCCAAAATAATGAAAAAAATAATAACAAACA
AT3G20900.1-SEQ                     ATCATCTCAAAAAACAAGCGCCAAAATAATGAAAAAAATAATAACAAAAA
AT3G20900.1-CDS                     ----------------------------ATGAACAAAGTAGCGAGGAAGA
                                                                ***** *** **   *  ** *

Test1seq                            AAAACAAACCAAAATAAGAAAAAACATTACGCAAAACATAATAATTTACT
AT3G20900.1-SEQ                     CAAACAAACCAAAATAAGAAAAAACATTACGCAAAACATAATAATTTACT
AT3G20900.1-CDS                     A------------------------------CAAAACATC----------
                                                                   ********           

Test1seq                            CTTCGTTATTGTATTAACAAATCAAAGAGCTGAATTTTGATCACCTGCTA
AT3G20900.1-SEQ                     CTTCGTTATTGTATTAACAAATCAAAGAGATGAATTTTGATCACCTGCTA
AT3G20900.1-CDS                     --------------------------------------------------
                                                                                      

Test1seq                            ATACTACTTTCTGTATTGATCCTATATCAACGTAAACAAAGATACTAATA
AT3G20900.1-SEQ                     ATACTACTTTCTGTATTGATCCTATATCAAAAAAAAAAAAGATACTAATA
AT3G20900.1-CDS                     --------------------------------------------------
                                                                                      

Test1seq                            ATTAACTAAAAGTACGTTCATCGATCGTGTTCGTTGACGAAGAAGAGCTC
AT3G20900.1-SEQ                     ATTAACTAAAAGTACGTTCATCGATCGTGTGCGTTGACGAAGAAGAGCTC
AT3G20900.1-CDS                     --------------------------------------------------
                                                                                      

Test1seq                            TATCTCCGGCGGAGCAAAGAAAACGATCTGTCTCCGTCGTAACACACGGT
AT3G20900.1-SEQ                     TATCTCCGGCGGAGCAAAGAAAACGATCTGTCTCCGTCGTAACACACAGT
AT3G20900.1-CDS                     ------------AGCAAAGAAAACGATCTGTCTCCGTCGTAACACACAGT
                                                *********************************** **

Test1seq                            CGCTAGAGAAACTTTGCTTCTTCGGCGCCGGTGGACACGTCAGCATCTCC
AT3G20900.1-SEQ                     TTTTCGAGACCCTTTGCTTCTTCGGCGCCGGTGGACACGTCAGCATCTCC
AT3G20900.1-CDS                     TTTTCGAGACCCTTTGCTTCTTCGGCGCCGGTGGACACGTCAGCATCTCC
                                       * ****  ***************************************

Test1seq                            GGTATCCTAGACTTCTTGGCTTTCGGGGTACAACAACCGCGTGGTGACGT
AT3G20900.1-SEQ                     GGTATCCTAGACTTCTTGGCTTTCGGGGTACAACAACCGCCTGGTGACGT
AT3G20900.1-CDS                     GGTATCCTAGACTTCTTGGCTTTCGGGGTACAACAACCGCCTGGTGACGT
                                    **************************************** *********

Test1seq                            CAGCACCGCTGCTGGGGATGGAGAGGGAACAGAGTT-
AT3G20900.1-SEQ                     CAGCACCGCTGCTGGGGATGGAGAGGGAACAGAGTAG
AT3G20900.1-CDS                     CAGCACCGCTGCTGGGGATGGAGAGGGAACAGAGTAG
                                    ***********************************  


""",
        )
        self.check_reading_writing(path)

    def test_kalign(self):
        """Make sure we can parse the Kalign header."""
        path = "Clustalw/kalign.aln"
        with open(path) as stream:
            alignments = Align.parse(stream, "clustal")
            self.assertEqual(alignments.metadata["Program"], "Kalign")
            self.assertEqual(alignments.metadata["Version"], "2.0")
            alignment = next(alignments)
            with self.assertRaises(StopIteration):
                next(alignments)
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['G', 'C', 'T', 'G', 'G', 'G', 'G', 'A', 'T', 'G', 'G', 'A', 'G',
              'A', 'G', 'G', 'G', 'A', 'A', 'C', 'A', 'G', 'A', 'G', 'T', '-',
              'T'],
             ['G', 'C', 'T', 'G', 'G', 'G', 'G', 'A', 'T', 'G', 'G', 'A', 'G',
              'A', 'G', 'G', 'G', 'A', 'A', 'C', 'A', 'G', 'A', 'G', 'T', 'A',
              'G']], dtype='U')
                # fmt: on
            )
        )
        self.assertEqual(
            repr(alignment),
            "<Alignment object (2 rows x 27 columns) at 0x%x>" % id(alignment),
        )
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.sequences[0].id, "Test1seq")
        self.assertEqual(alignment.sequences[1].id, "AT3G20900")
        self.assertEqual(alignment.sequences[0].seq, "GCTGGGGATGGAGAGGGAACAGAGTT")
        self.assertEqual(alignment.sequences[1].seq, "GCTGGGGATGGAGAGGGAACAGAGTAG")
        self.assertEqual(alignment[0], "GCTGGGGATGGAGAGGGAACAGAGT-T")
        self.assertEqual(alignment[1], "GCTGGGGATGGAGAGGGAACAGAGTAG")
        self.assertEqual(
            format(alignment, "clustal"),
            """\
Test1seq                            GCTGGGGATGGAGAGGGAACAGAGT-T
AT3G20900                           GCTGGGGATGGAGAGGGAACAGAGTAG


""",
        )
        self.check_reading_writing(path)

    def test_probcons(self):
        path = "Clustalw/probcons.aln"
        # example taken from the PROBCONS documentation
        with open(path) as stream:
            alignments = Align.parse(stream, "clustal")
            self.assertEqual(alignments.metadata["Program"], "PROBCONS")
            self.assertEqual(alignments.metadata["Version"], "1.12")
            alignment = next(alignments)
            with self.assertRaises(StopIteration):
                next(alignments)
        self.assertTrue(
            numpy.array_equal(
                numpy.array(alignment, "U"),
                # fmt: off
# flake8: noqa
numpy.array([['D', '-', 'V', 'L', 'L', 'G', 'A', 'N', 'G', 'G', 'V', 'L', 'V',
              'F', 'E', 'P', 'N', 'D', 'F', 'S', 'V', 'K', 'A', 'G', 'E', 'T',
              'I', 'T', 'F', 'K', 'N', 'N', 'A', 'G', 'Y', 'P', 'H', 'N', 'V',
              'V', 'F', 'D', 'E', 'D', 'A', 'V', 'P', 'S', 'G', '-', 'V', 'D',
              '-', 'V', 'S', 'K', 'I', 'S', 'Q', 'E', 'E', 'Y', 'L', 'T', 'A',
              'P', 'G', 'E', 'T', 'F', 'S', 'V', 'T', 'L', 'T', 'V', '-', '-',
              '-', 'P', 'G', 'T', 'Y', 'G', 'F', 'Y', 'C', 'E', 'P', 'H', 'A',
              'G', 'A', 'G', 'M', 'V', 'G', 'K', 'V', 'T', 'V'],
             ['-', '-', 'V', 'K', 'L', 'G', 'A', 'D', 'S', 'G', 'A', 'L', 'E',
              'F', 'V', 'P', 'K', 'T', 'L', 'T', 'I', 'K', 'S', 'G', 'E', 'T',
              'V', 'N', 'F', 'V', 'N', 'N', 'A', 'G', 'F', 'P', 'H', 'N', 'I',
              'V', 'F', 'D', 'E', 'D', 'A', 'I', 'P', 'S', 'G', '-', 'V', 'N',
              '-', 'A', 'D', 'A', 'I', 'S', 'R', 'D', 'D', 'Y', 'L', 'N', 'A',
              'P', 'G', 'E', 'T', 'Y', 'S', 'V', 'K', 'L', 'T', 'A', '-', '-',
              '-', 'A', 'G', 'E', 'Y', 'G', 'Y', 'Y', 'C', 'E', 'P', 'H', 'Q',
              'G', 'A', 'G', 'M', 'V', 'G', 'K', 'I', 'I', 'V'],
             ['-', '-', 'V', 'K', 'L', 'G', 'S', 'D', 'K', 'G', 'L', 'L', 'V',
              'F', 'E', 'P', 'A', 'K', 'L', 'T', 'I', 'K', 'P', 'G', 'D', 'T',
              'V', 'E', 'F', 'L', 'N', 'N', 'K', 'V', 'P', 'P', 'H', 'N', 'V',
              'V', 'F', 'D', 'A', 'A', 'L', 'N', 'P', 'A', 'K', 'S', 'A', 'D',
              'L', 'A', 'K', 'S', 'L', 'S', 'H', 'K', 'Q', 'L', 'L', 'M', 'S',
              'P', 'G', 'Q', 'S', 'T', 'S', 'T', 'T', 'F', 'P', 'A', 'D', 'A',
              'P', 'A', 'G', 'E', 'Y', 'T', 'F', 'Y', 'C', 'E', 'P', 'H', 'R',
              'G', 'A', 'G', 'M', 'V', 'G', 'K', 'I', 'T', 'V'],
             ['V', 'Q', 'I', 'K', 'M', 'G', 'T', 'D', 'K', 'Y', 'A', 'P', 'L',
              'Y', 'E', 'P', 'K', 'A', 'L', 'S', 'I', 'S', 'A', 'G', 'D', 'T',
              'V', 'E', 'F', 'V', 'M', 'N', 'K', 'V', 'G', 'P', 'H', 'N', 'V',
              'I', 'F', 'D', 'K', '-', '-', 'V', 'P', 'A', 'G', '-', 'E', 'S',
              '-', 'A', 'P', 'A', 'L', 'S', 'N', 'T', 'K', 'L', 'R', 'I', 'A',
              'P', 'G', 'S', 'F', 'Y', 'S', 'V', 'T', 'L', 'G', 'T', '-', '-',
              '-', 'P', 'G', 'T', 'Y', 'S', 'F', 'Y', 'C', 'T', 'P', 'H', 'R',
              'G', 'A', 'G', 'M', 'V', 'G', 'T', 'I', 'T', 'V'],
             ['V', 'H', 'M', 'L', 'N', 'K', 'G', 'K', 'D', 'G', 'A', 'M', 'V',
              'F', 'E', 'P', 'A', 'S', 'L', 'K', 'V', 'A', 'P', 'G', 'D', 'T',
              'V', 'T', 'F', 'I', 'P', 'T', 'D', 'K', '-', 'G', 'H', 'N', 'V',
              'E', 'T', 'I', 'K', 'G', 'M', 'I', 'P', 'D', 'G', '-', 'A', 'E',
              '-', 'A', '-', '-', '-', '-', '-', '-', '-', 'F', 'K', 'S', 'K',
              'I', 'N', 'E', 'N', 'Y', 'K', 'V', 'T', 'F', 'T', 'A', '-', '-',
              '-', 'P', 'G', 'V', 'Y', 'G', 'V', 'K', 'C', 'T', 'P', 'H', 'Y',
              'G', 'M', 'G', 'M', 'V', 'G', 'V', 'V', 'E', 'V']], dtype='U')
                # fmt: on
            )
        )
        self.assertEqual(
            repr(alignment),
            "<Alignment object (5 rows x 101 columns) at 0x%x>" % id(alignment),
        )
        self.assertEqual(len(alignment), 5)
        self.assertEqual(alignment.sequences[0].id, "plas_horvu")
        self.assertEqual(alignment.sequences[1].id, "plas_chlre")
        self.assertEqual(alignment.sequences[2].id, "plas_anava")
        self.assertEqual(alignment.sequences[3].id, "plas_proho")
        self.assertEqual(alignment.sequences[4].id, "azup_achcy")
        self.assertEqual(
            alignment.sequences[0].seq,
            "DVLLGANGGVLVFEPNDFSVKAGETITFKNNAGYPHNVVFDEDAVPSGVDVSKISQEEYLTAPGETFSVTLTVPGTYGFYCEPHAGAGMVGKVTV",
        )
        self.assertEqual(
            alignment.sequences[1].seq,
            "VKLGADSGALEFVPKTLTIKSGETVNFVNNAGFPHNIVFDEDAIPSGVNADAISRDDYLNAPGETYSVKLTAAGEYGYYCEPHQGAGMVGKIIV",
        )
        self.assertEqual(
            alignment.sequences[2].seq,
            "VKLGSDKGLLVFEPAKLTIKPGDTVEFLNNKVPPHNVVFDAALNPAKSADLAKSLSHKQLLMSPGQSTSTTFPADAPAGEYTFYCEPHRGAGMVGKITV",
        )
        self.assertEqual(
            alignment.sequences[3].seq,
            "VQIKMGTDKYAPLYEPKALSISAGDTVEFVMNKVGPHNVIFDKVPAGESAPALSNTKLRIAPGSFYSVTLGTPGTYSFYCTPHRGAGMVGTITV",
        )
        self.assertEqual(
            alignment.sequences[4].seq,
            "VHMLNKGKDGAMVFEPASLKVAPGDTVTFIPTDKGHNVETIKGMIPDGAEAFKSKINENYKVTFTAPGVYGVKCTPHYGMGMVGVVEV",
        )
        self.assertEqual(
            alignment.column_annotations["clustal_consensus"],
            " ::    .     : *  :.: .*:*: *  .    **:       *    . .  :*. .     ..  ...: .   .* *   * ** * **** : *",
        )
        self.assertEqual(
            alignment[0],
            "D-VLLGANGGVLVFEPNDFSVKAGETITFKNNAGYPHNVVFDEDAVPSG-VD-VSKISQEEYLTAPGETFSVTLTV---PGTYGFYCEPHAGAGMVGKVTV",
        )
        self.assertEqual(
            alignment[1],
            "--VKLGADSGALEFVPKTLTIKSGETVNFVNNAGFPHNIVFDEDAIPSG-VN-ADAISRDDYLNAPGETYSVKLTA---AGEYGYYCEPHQGAGMVGKIIV",
        )
        self.assertEqual(
            alignment[2],
            "--VKLGSDKGLLVFEPAKLTIKPGDTVEFLNNKVPPHNVVFDAALNPAKSADLAKSLSHKQLLMSPGQSTSTTFPADAPAGEYTFYCEPHRGAGMVGKITV",
        )
        self.assertEqual(
            alignment[3],
            "VQIKMGTDKYAPLYEPKALSISAGDTVEFVMNKVGPHNVIFDK--VPAG-ES-APALSNTKLRIAPGSFYSVTLGT---PGTYSFYCTPHRGAGMVGTITV",
        )
        self.assertEqual(
            alignment[4],
            "VHMLNKGKDGAMVFEPASLKVAPGDTVTFIPTDK-GHNVETIKGMIPDG-AE-A-------FKSKINENYKVTFTA---PGVYGVKCTPHYGMGMVGVVEV",
        )
        self.assertEqual(
            format(alignment, "clustal"),
            """\
plas_horvu                          D-VLLGANGGVLVFEPNDFSVKAGETITFKNNAGYPHNVVFDEDAVPSG-
plas_chlre                          --VKLGADSGALEFVPKTLTIKSGETVNFVNNAGFPHNIVFDEDAIPSG-
plas_anava                          --VKLGSDKGLLVFEPAKLTIKPGDTVEFLNNKVPPHNVVFDAALNPAKS
plas_proho                          VQIKMGTDKYAPLYEPKALSISAGDTVEFVMNKVGPHNVIFDK--VPAG-
azup_achcy                          VHMLNKGKDGAMVFEPASLKVAPGDTVTFIPTDK-GHNVETIKGMIPDG-
                                     ::    .     : *  :.: .*:*: *  .    **:       *   

plas_horvu                          VD-VSKISQEEYLTAPGETFSVTLTV---PGTYGFYCEPHAGAGMVGKVT
plas_chlre                          VN-ADAISRDDYLNAPGETYSVKLTA---AGEYGYYCEPHQGAGMVGKII
plas_anava                          ADLAKSLSHKQLLMSPGQSTSTTFPADAPAGEYTFYCEPHRGAGMVGKIT
plas_proho                          ES-APALSNTKLRIAPGSFYSVTLGT---PGTYSFYCTPHRGAGMVGTIT
azup_achcy                          AE-A-------FKSKINENYKVTFTA---PGVYGVKCTPHYGMGMVGVVE
                                     . .  :*. .     ..  ...: .   .* *   * ** * **** : 

plas_horvu                          V
plas_chlre                          V
plas_anava                          V
plas_proho                          V
azup_achcy                          V
                                    *


""",
        )
        self.check_reading_writing(path)

    def test_empty(self):
        """Checking empty file."""
        stream = StringIO()
        with self.assertRaises(ValueError):
            Align.parse(stream, "clustal")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
