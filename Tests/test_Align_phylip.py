# Copyright 2006-2014 by Peter Cock.  All rights reserved.
# Revisions copyright 2011 Brandon Invergo. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.Align.phylip module."""
import unittest

from io import StringIO

from Bio.Align.phylip import AlignmentIterator
from Bio.Align.phylip import AlignmentWriter


class TestPhylipReading(unittest.TestCase):

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

    def test_one(self):
        path = "Phylip/one.dat"
        with open(path) as stream:
            alignments = AlignmentIterator(stream)
            alignment = next(alignments)
            with self.assertRaises(StopIteration):
                next(alignments)
        self.assertEqual(
            repr(alignment),
            "<Bio.Align.Alignment object (8 rows x 286 columns) at 0x%x>"
            % id(alignment),
        )
        self.assertEqual(len(alignment), 8)
        self.assertEqual(alignment.sequences[0].id, "V_Harveyi_")
        self.assertEqual(alignment.sequences[1].id, "B_subtilis")
        self.assertEqual(alignment.sequences[2].id, "B_subtilis")
        self.assertEqual(alignment.sequences[3].id, "YA80_HAEIN")
        self.assertEqual(alignment.sequences[4].id, "FLIY_ECOLI")
        self.assertEqual(alignment.sequences[5].id, "E_coli_Gln")
        self.assertEqual(alignment.sequences[6].id, "Deinococcu")
        self.assertEqual(alignment.sequences[7].id, "HISJ_E_COL")
        self.assertEqual(alignment.sequences[0].seq, "MKNWIKVAVAAIALSAATVQAATEVKVGMSGRYFPFTFVKQDKLQGFEVDMWDEIGKRNDYKIEYVTANFSGLFGLLETGRIDTISNQITMTDARKAKYLFADPYVVDGAQITVRKGNDSIQGVEDLAGKTVAVNLGSNFEQLLRDYDKDGKINIKTYDTGIEHDVALGRADAFIMDRLSALELIKKTGLPLQLAGEPFETIQNAWPFVDNEKGRKLQAEVNKALAEMRADGTVEKISVKWFGADITK")
        self.assertEqual(alignment.sequences[1].seq, "MKMKKWTVLVVAALLAVLSACGNGNSSSKEDDNVLHVGATGQSYPFAYKENGKLTGFDVEVMEAVAKKIDMKLDWKLLEFSGLMGELQTGKLDTISNQVAVTDERKETYNFTKPYAYAGTQIVVKKDNTDIKSVDDLKGKTVAAVLGSNHAKNLESKDPDKKINIKTYETQEGTLKDVAYGRVDAYVNSRTVLIAQIKKTGLPLKLAGDPIVYEQVAFPFAKDDAHDKLRKKVNKALDELRKDGTLKKLSEKYFNEDITVEQKH")
        self.assertEqual(alignment.sequences[2].seq, "MKKALLALFMVVSIAALAACGAGNDNQSKDNAKDGDLWASIKKKGVLTVGTEGTYEPFTYHDKDTDKLTGYDVEVITEVAKRLGLKVDFKETQWGSMFAGLNSKRFDVVANQVGKTDREDKYDFSDKYTTSRAVVVTKKDNNDIKSEADVKGKTSAQSLTSNYNKLATNAGAKVEGVEGMAQALQMIQQARVDMTYNDKLAVLNYLKTSGNKNVKIAFETGEPQSTYFTFRKGSGEVVDQVNKALKEMKEDGTLSKISKKWFGEDVSK")
        self.assertEqual(alignment.sequences[3].seq, "MKKLLFTTALLTGAIAFSTFSHAGEIADRVEKTKTLLVGTEGTYAPFTFHDKSGKLTGFDVEVIRKVAEKLGLKVEFKETQWDAMYAGLNAKRFDVIANQTNPSPERLKKYSFTTPYNYSGGVIVTKSSDNSIKSFEDLKGRKSAQSATSNWGKDAKAAGAQILVVDGLAQSLELIKQGRAEATINDKLAVLDYFKQHPNSGLKIAYDRGDKTPTAFAFLQGEDALITKFNQVLEALRQDGTLKQISIEWFGYDITQ")
        self.assertEqual(alignment.sequences[4].seq, "MKLAHLGRQALMGVMAVALVAGMSVKSFADEGLLNKVKERGTLLVGLEGTYPPFSFQGDDGKLTGFEVEFAQQLAKHLGVEASLKPTKWDGMLASLDSKRIDVVINQVTISDERKKKYDFSTPYTISGIQALVKKGNEGTIKTADDLKGKKVGVGLGTNYEEWLRQNVQGVDVRTYDDDPTKYQDLRVGRIDAILVDRLAALDLVKKTNDTLAVTGEAFSRQESGVALRKGNEDLLKAVNDAIAEMQKDGTLQALSEKWFGADVTK")
        self.assertEqual(alignment.sequences[5].seq, "MKSVLKVSLAALTLAFAVSSHAADKKLVVATDTAFVPFEFKQGDKYVGFDVDLWAAIAKELKLDYELKPMDFSGIIPALQTKNVDLALAGITITDERKKAIDFSDGYYKSGLLVMVKANNNDVKSVKDLDGKVVAVKSGTGSVDYAKANIKTKDLRQFPNIDNAYMELGTNRADAVLHDTPNILYFIKTAGNGQFKAVGDSLEAQQYGIAFPKGSDELRDKVNGALKTLRENGTYNEIYKKWFGTEPK")
        self.assertEqual(alignment.sequences[6].seq, "MKKSLLSLKLSGLLVPSVLALSLSACSSPSSTLNQGTLKIAMEGTYPPFTSKNEQGELVGFDVDIAKAVAQKLNLKPEFVLTEWSGILAGLQANKYDVIVNQVGITPERQNSIGFSQPYAYSRPEIIVAKNNTFNPQSLADLKGKRVGSTLGSNYEKQLIDTGDIKIVTYPGAPEILADLVAGRIDAAYNDRLVVNYIINDQKLPVRGAGQIGDAAPVGIALKKGNSALKDQIDKALTEMRSDGTFEKISQKWFGQDVGQP")
        self.assertEqual(alignment.sequences[7].seq, "MKKLVLSLSLVLAFSSATAAFAAIPQNIRIGTDPTYAPFESKNSQGELVGFDIDLAKELCKRINTQCTFVENPLDALIPSLKAKKIDAIMSSLSITEKRQQEIAFTDKLYAADSRLVVAKNSDIQPTVESLKGKRVGVLQGTTQETFGNEHWAPKGIEIVSYQGQDNIYSDLTAGRIDAAFQDEVAASEGFLKQPVGKDYKFGGPSVKDEKLFGVGTGMGLRKEDNELREALNKAFAEMRADGTYEKLAKKYFDFDVYGG")
        self.assertEqual(alignment[0], "--MKNWIKVAVAAIA--LSAA------------------TVQAATEVKVGMSGRYFPFTFVKQ--DKLQGFEVDMWDEIGKRNDYKIEYVTANFSGLFGLLETGRIDTISNQITMTDARKAKYLFADPYVVDG-AQITVRKGNDSIQGVEDLAGKTVAVNLGSNFEQLLRDYDKDGKINIKTYDT--GIEHDVALGRADAFIMDRLSALE-LIKKT-GLPLQLAGEPFETI-----QNAWPFVDNEKGRKLQAEVNKALAEMRADGTVEKISVKWFGADITK----")
        self.assertEqual(alignment[1], "MKMKKWTVLVVAALLAVLSACG------------NGNSSSKEDDNVLHVGATGQSYPFAYKEN--GKLTGFDVEVMEAVAKKIDMKLDWKLLEFSGLMGELQTGKLDTISNQVAVTDERKETYNFTKPYAYAG-TQIVVKKDNTDIKSVDDLKGKTVAAVLGSNHAKNLESKDPDKKINIKTYETQEGTLKDVAYGRVDAYVNSRTVLIA-QIKKT-GLPLKLAGDPIVYE-----QVAFPFAKDDAHDKLRKKVNKALDELRKDGTLKKLSEKYFNEDITVEQKH")
        self.assertEqual(alignment[2], "MKKALLALFMVVSIAALAACGAGNDNQSKDNAKDGDLWASIKKKGVLTVGTEGTYEPFTYHDKDTDKLTGYDVEVITEVAKRLGLKVDFKETQWGSMFAGLNSKRFDVVANQVG-KTDREDKYDFSDKYTTSR-AVVVTKKDNNDIKSEADVKGKTSAQSLTSNYNKLATN----AGAKVEGVEGMAQALQMIQQARVDMTYNDKLAVLN-YLKTSGNKNVKIAFETGEPQ-----STYFTFRKGS--GEVVDQVNKALKEMKEDGTLSKISKKWFGEDVSK----")
        self.assertEqual(alignment[3], "MKKLLFTTALLTGAIAFSTF-----------SHAGEIADRVEKTKTLLVGTEGTYAPFTFHDK-SGKLTGFDVEVIRKVAEKLGLKVEFKETQWDAMYAGLNAKRFDVIANQTNPSPERLKKYSFTTPYNYSG-GVIVTKSSDNSIKSFEDLKGRKSAQSATSNWGKDAKA----AGAQILVVDGLAQSLELIKQGRAEATINDKLAVLD-YFKQHPNSGLKIAYDRGDKT-----PTAFAFLQGE--DALITKFNQVLEALRQDGTLKQISIEWFGYDITQ----")
        self.assertEqual(alignment[4], "MKLAHLGRQALMGVMAVALVAG---MSVKSFADEG-LLNKVKERGTLLVGLEGTYPPFSFQGD-DGKLTGFEVEFAQQLAKHLGVEASLKPTKWDGMLASLDSKRIDVVINQVTISDERKKKYDFSTPYTISGIQALVKKGNEGTIKTADDLKGKKVGVGLGTNYEEWLRQNV--QGVDVRTYDDDPTKYQDLRVGRIDAILVDRLAALD-LVKKT-NDTLAVTGEAFSRQ-----ESGVALRKGN--EDLLKAVNDAIAEMQKDGTLQALSEKWFGADVTK----")
        self.assertEqual(alignment[5], "--MKSVLKVSLAALTLAFAVS------------------SHAADKKLVVATDTAFVPFEFKQG--DKYVGFDVDLWAAIAKELKLDYELKPMDFSGIIPALQTKNVDLALAGITITDERKKAIDFSDGYYKSG-LLVMVKANNNDVKSVKDLDGKVVAVKSGTGSVDYAKAN--IKTKDLRQFPNIDNAYMELGTNRADAVLHDTPNILY-FIKTAGNGQFKAVGDSLEAQ-----QYGIAFPKGS--DELRDKVNGALKTLRENGTYNEIYKKWFGTEPK-----")
        self.assertEqual(alignment[6], "-MKKSLLSLKLSGLLVPSVLALS--------LSACSSPSSTLNQGTLKIAMEGTYPPFTSKNE-QGELVGFDVDIAKAVAQKLNLKPEFVLTEWSGILAGLQANKYDVIVNQVGITPERQNSIGFSQPYAYSRPEIIVAKNNTFNPQSLADLKGKRVGSTLGSNYEKQLIDTG---DIKIVTYPGAPEILADLVAGRIDAAYNDRLVVNY-IINDQ-KLPVRGAGQIGDAA-----PVGIALKKGN--SALKDQIDKALTEMRSDGTFEKISQKWFGQDVGQP---")
        self.assertEqual(alignment[7], "MKKLVLSLSLVLAFSSATAAF-------------------AAIPQNIRIGTDPTYAPFESKNS-QGELVGFDIDLAKELCKRINTQCTFVENPLDALIPSLKAKKIDAIMSSLSITEKRQQEIAFTDKLYAADSRLVVAKNSDIQP-TVESLKGKRVGVLQGTTQETFGNEHWAPKGIEIVSYQGQDNIYSDLTAGRIDAAFQDEVAASEGFLKQPVGKDYKFGGPSVKDEKLFGVGTGMGLRKED--NELREALNKAFAEMRADGTYEKLAKKYFDFDVYGG---")
        self.check_reading_writing(path)

    def test_two_and_three(self):
        paths = ("Phylip/two.dat", "Phylip/three.dat")
        # derived from http://atgc.lirmm.fr/phyml/usersguide.html
        for path in paths:
            with open(path) as stream:
                alignments = AlignmentIterator(stream)
                alignment = next(alignments)
                with self.assertRaises(StopIteration):
                    next(alignments)
            self.assertEqual(
                repr(alignment),
                "<Bio.Align.Alignment object (5 rows x 60 columns) at 0x%x>"
                % id(alignment),
            )
            self.assertEqual(len(alignment), 5)
            self.assertEqual(alignment.sequences[0].id, "Tax1")
            self.assertEqual(alignment.sequences[1].id, "Tax2")
            self.assertEqual(alignment.sequences[2].id, "Tax3")
            self.assertEqual(alignment.sequences[3].id, "Tax4")
            self.assertEqual(alignment.sequences[4].id, "Tax5")
            self.assertEqual(
                alignment.sequences[0].seq, "CCATCTCACGGTCGGTACGATACACCTGCTTTTGGCAGGAAATGGTCAATATTACAAGGT")
            self.assertEqual(alignment.sequences[1].seq, "CCATCTCACGGTCAGTAAGATACACCTGCTTTTGGCGGGAAATGGTCAACATTAAAAGAT")
            self.assertEqual(alignment.sequences[2].seq, "CCATCTCCCGCTCAGTAAGATACCCCTGCTGTTGGCGGGAAATCGTCAATATTAAAAGGT")
            self.assertEqual(alignment.sequences[3].seq, "TCATCTCATGGTCAATAAGATACTCCTGCTTTTGGCGGGAAATGGTCAATCTTAAAAGGT")
            self.assertEqual(alignment.sequences[4].seq, "CCATCTCACGGTCGGTAAGATACACCTGCTTTTGGCGGGAAATGGTCAATATTAAAAGGT")
            self.assertEqual(
                alignment[0], "CCATCTCACGGTCGGTACGATACACCTGCTTTTGGCAGGAAATGGTCAATATTACAAGGT")
            self.assertEqual(alignment[1], "CCATCTCACGGTCAGTAAGATACACCTGCTTTTGGCGGGAAATGGTCAACATTAAAAGAT")
            self.assertEqual(alignment[2], "CCATCTCCCGCTCAGTAAGATACCCCTGCTGTTGGCGGGAAATCGTCAATATTAAAAGGT")
            self.assertEqual(alignment[3], "TCATCTCATGGTCAATAAGATACTCCTGCTTTTGGCGGGAAATGGTCAATCTTAAAAGGT")
            self.assertEqual(alignment[4], "CCATCTCACGGTCGGTAAGATACACCTGCTTTTGGCGGGAAATGGTCAATATTAAAAGGT")
            self.check_reading_writing(path)

    def test_four(self):
        path = "Phylip/four.dat"
        # File derived from here:
        # http://evolution.genetics.washington.edu/phylip/doc/sequence.html
        # Note the lack of any white space between names 2 and 3 and their seqs.
        with open(path) as stream:
            alignments = AlignmentIterator(stream)
            alignment = next(alignments)
            with self.assertRaises(StopIteration):
                next(alignments)
        self.assertEqual(
            repr(alignment),
            "<Bio.Align.Alignment object (5 rows x 42 columns) at 0x%x>"
            % id(alignment),
        )
        self.assertEqual(len(alignment), 5)
        self.assertEqual(alignment.sequences[0].id, "Turkey")
        self.assertEqual(alignment.sequences[1].id, "Salmo gair")
        self.assertEqual(alignment.sequences[2].id, "H. Sapiens")
        self.assertEqual(alignment.sequences[3].id, "Chimp")
        self.assertEqual(alignment.sequences[4].id, "Gorilla")
        self.assertEqual(alignment.sequences[0].seq, "AAGCTNGGGCATTTCAGGGTGAGCCCGGGCAATACAGGGTAT")
        self.assertEqual(alignment.sequences[1].seq, "AAGCCTTGGCAGTGCAGGGTGAGCCGTGGCCGGGCACGGTAT")
        self.assertEqual(alignment.sequences[2].seq, "ACCGGTTGGCCGTTCAGGGTACAGGTTGGCCGTTCAGGGTAA")
        self.assertEqual(alignment.sequences[3].seq, "AAACCCTTGCCGTTACGCTTAAACCGAGGCCGGGACACTCAT")
        self.assertEqual(alignment.sequences[4].seq, "AAACCCTTGCCGGTACGCTTAAACCATTGCCGGTACGCTTAA")
        self.assertEqual(alignment[0], "AAGCTNGGGCATTTCAGGGTGAGCCCGGGCAATACAGGGTAT")
        self.assertEqual(alignment[1], "AAGCCTTGGCAGTGCAGGGTGAGCCGTGGCCGGGCACGGTAT")
        self.assertEqual(alignment[2], "ACCGGTTGGCCGTTCAGGGTACAGGTTGGCCGTTCAGGGTAA")
        self.assertEqual(alignment[3], "AAACCCTTGCCGTTACGCTTAAACCGAGGCCGGGACACTCAT")
        self.assertEqual(alignment[4], "AAACCCTTGCCGGTACGCTTAAACCATTGCCGGTACGCTTAA")
        self.check_reading_writing(path)

    def test_five(self):
        # File derived from here:
        # http://evolution.genetics.washington.edu/phylip/doc/sequence.html
        path = "Phylip/five.dat"
        with open(path) as stream:
            self.assertRaises(ValueError, next, AlignmentIterator(stream))

    def test_five_a(self):
        # File derived from here:
        # http://evolution.genetics.washington.edu/phylip/doc/sequence.html
        path = "Phylip/five_a.dat"
        with open(path) as stream:
            alignments = AlignmentIterator(stream)
            alignment = next(alignments)
            with self.assertRaises(StopIteration):
                next(alignments)
        self.assertEqual(
            repr(alignment),
            "<Bio.Align.Alignment object (5 rows x 42 columns) at 0x%x>"
            % id(alignment),
        )
        self.assertEqual(len(alignment), 5)
        self.assertEqual(alignment.sequences[0].id, "Turkey")
        self.assertEqual(alignment.sequences[1].id, "Salmo gair")
        self.assertEqual(alignment.sequences[2].id, "H. Sapiens")
        self.assertEqual(alignment.sequences[3].id, "Chimp")
        self.assertEqual(alignment.sequences[4].id, "Gorilla")
        self.assertEqual(alignment.sequences[0].seq, "AAGCTNGGGCATTTCAGGGTGAGCCCGGGCAATACAGGGTAT")
        self.assertEqual(alignment.sequences[1].seq, "AAGCCTTGGCAGTGCAGGGTGAGCCGTGGCCGGGCACGGTAT")
        self.assertEqual(alignment.sequences[2].seq, "ACCGGTTGGCCGTTCAGGGTACAGGTTGGCCGTTCAGGGTAA")
        self.assertEqual(alignment.sequences[3].seq, "AAACCCTTGCCGTTACGCTTAAACCGAGGCCGGGACACTCAT")
        self.assertEqual(alignment.sequences[4].seq, "AAACCCTTGCCGGTACGCTTAAACCATTGCCGGTACGCTTAA")
        self.assertEqual(alignment[0], "AAGCTNGGGCATTTCAGGGTGAGCCCGGGCAATACAGGGTAT")
        self.assertEqual(alignment[1], "AAGCCTTGGCAGTGCAGGGTGAGCCGTGGCCGGGCACGGTAT")
        self.assertEqual(alignment[2], "ACCGGTTGGCCGTTCAGGGTACAGGTTGGCCGTTCAGGGTAA")
        self.assertEqual(alignment[3], "AAACCCTTGCCGTTACGCTTAAACCGAGGCCGGGACACTCAT")
        self.assertEqual(alignment[4], "AAACCCTTGCCGGTACGCTTAAACCATTGCCGGTACGCTTAA")
        self.check_reading_writing(path)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
