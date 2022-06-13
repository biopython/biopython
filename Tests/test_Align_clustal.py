# Copyright 2006-2014 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.Align.clustal module."""
import unittest
import warnings

from io import StringIO

from Bio import BiopythonExperimentalWarning

with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonExperimentalWarning)
    from Bio.Align.clustal import AlignmentIterator
    from Bio.Align.clustal import AlignmentWriter


class TestClustalReadingWriting(unittest.TestCase):
    def check_reading_writing(self, path):
        alignments = AlignmentIterator(path)
        stream = StringIO()
        writer = AlignmentWriter(stream)
        n = writer.write_file(alignments, mincount=1, maxcount=1)
        self.assertEqual(n, 1)
        alignments = AlignmentIterator(path)
        alignment = next(alignments)
        stream.seek(0)
        saved_alignments = AlignmentIterator(stream)
        self.assertEqual(saved_alignments.program, alignments.program)
        self.assertEqual(saved_alignments.version, alignments.version)
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
            alignments = AlignmentIterator(stream)
            self.assertEqual(alignments.program, "CLUSTAL")
            self.assertEqual(alignments.version, "1.81")
            alignment = next(alignments)
            with self.assertRaises(StopIteration):
                next(alignments)
        self.assertEqual(
            repr(alignment),
            "<Bio.Align.Alignment object (2 rows x 601 columns) at 0x%x>"
            % id(alignment),
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
        self.check_reading_writing(path)

    def test_msaprobs(self):
        path = "Clustalw/msaprobs.aln"
        # This example was obtained from
        # http://virgil.ruc.dk/kurser/Sekvens/Treedraw.htm
        with open(path) as stream:
            alignments = AlignmentIterator(stream)
            self.assertEqual(alignments.program, "MSAPROBS")
            self.assertEqual(alignments.version, "0.9.7")
            alignment = next(alignments)
            with self.assertRaises(StopIteration):
                next(alignments)
        self.assertEqual(
            repr(alignment),
            "<Bio.Align.Alignment object (8 rows x 298 columns) at 0x%x>"
            % id(alignment),
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
        self.check_reading_writing(path)

    def test_muscle(self):
        path = "Clustalw/muscle.aln"
        # includes the sequence length on the right hand side of each line
        with open(path) as stream:
            alignments = AlignmentIterator(stream)
            self.assertEqual(alignments.program, "MUSCLE")
            self.assertEqual(alignments.version, "3.8")
            alignment = next(alignments)
            with self.assertRaises(StopIteration):
                next(alignments)
        self.assertEqual(
            repr(alignment),
            "<Bio.Align.Alignment object (3 rows x 687 columns) at 0x%x>"
            % id(alignment),
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
        self.check_reading_writing(path)

    def test_kalign(self):
        """Make sure we can parse the Kalign header."""
        path = "Clustalw/kalign.aln"
        with open(path) as stream:
            alignments = AlignmentIterator(stream)
            self.assertEqual(alignments.program, "Kalign")
            self.assertEqual(alignments.version, "2.0")
            alignment = next(alignments)
            with self.assertRaises(StopIteration):
                next(alignments)
        self.assertEqual(
            repr(alignment),
            "<Bio.Align.Alignment object (2 rows x 27 columns) at 0x%x>"
            % id(alignment),
        )
        self.assertEqual(len(alignment), 2)
        self.assertEqual(alignment.sequences[0].id, "Test1seq")
        self.assertEqual(alignment.sequences[1].id, "AT3G20900")
        self.assertEqual(alignment.sequences[0].seq, "GCTGGGGATGGAGAGGGAACAGAGTT")
        self.assertEqual(alignment.sequences[1].seq, "GCTGGGGATGGAGAGGGAACAGAGTAG")
        self.assertEqual(alignment[0], "GCTGGGGATGGAGAGGGAACAGAGT-T")
        self.assertEqual(alignment[1], "GCTGGGGATGGAGAGGGAACAGAGTAG")
        self.check_reading_writing(path)

    def test_probcons(self):
        path = "Clustalw/probcons.aln"
        # example taken from the PROBCONS documentation
        with open(path) as stream:
            alignments = AlignmentIterator(stream)
            self.assertEqual(alignments.program, "PROBCONS")
            self.assertEqual(alignments.version, "1.12")
            alignment = next(alignments)
            with self.assertRaises(StopIteration):
                next(alignments)
        self.assertEqual(
            repr(alignment),
            "<Bio.Align.Alignment object (5 rows x 101 columns) at 0x%x>"
            % id(alignment),
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
        self.check_reading_writing(path)

    def test_empty(self):
        """Checking empty file."""
        stream = StringIO()
        with self.assertRaises(ValueError):
            AlignmentIterator(stream)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
