# Copyright 2006-2014 by Peter Cock.  All rights reserved.
# Revisions copyright 2011 Brandon Invergo. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.Align.phylip module."""
import unittest
import warnings

from io import StringIO

from Bio import BiopythonExperimentalWarning

with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonExperimentalWarning)
    from Bio.Align.phylip import AlignmentIterator
    from Bio.Align.phylip import AlignmentWriter


try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align.phylip."
    ) from None


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
        self.assertTrue(
            numpy.array_equal(alignment.coordinates, saved_alignment.coordinates)
        )

    def test_one(self):
        path = "Phylip/one.dat"
        with open(path) as stream:
            alignments = AlignmentIterator(stream)
            alignment = next(alignments)
            with self.assertRaises(StopIteration):
                next(alignments)
        self.assertEqual(
            repr(alignment),
            "<Alignment object (8 rows x 286 columns) at 0x%x>"
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
            "MKKALLALFMVVSIAALAACGAGNDNQSKDNAKDGDLWASIKKKGVLTVGTEGTYEPFTYHDKDTDKLTGYDVEVITEVAKRLGLKVDFKETQWGSMFAGLNSKRFDVVANQVGKTDREDKYDFSDKYTTSRAVVVTKKDNNDIKSEADVKGKTSAQSLTSNYNKLATNAGAKVEGVEGMAQALQMIQQARVDMTYNDKLAVLNYLKTSGNKNVKIAFETGEPQSTYFTFRKGSGEVVDQVNKALKEMKEDGTLSKISKKWFGEDVSK",
        )
        self.assertEqual(
            alignment.sequences[3].seq,
            "MKKLLFTTALLTGAIAFSTFSHAGEIADRVEKTKTLLVGTEGTYAPFTFHDKSGKLTGFDVEVIRKVAEKLGLKVEFKETQWDAMYAGLNAKRFDVIANQTNPSPERLKKYSFTTPYNYSGGVIVTKSSDNSIKSFEDLKGRKSAQSATSNWGKDAKAAGAQILVVDGLAQSLELIKQGRAEATINDKLAVLDYFKQHPNSGLKIAYDRGDKTPTAFAFLQGEDALITKFNQVLEALRQDGTLKQISIEWFGYDITQ",
        )
        self.assertEqual(
            alignment.sequences[4].seq,
            "MKLAHLGRQALMGVMAVALVAGMSVKSFADEGLLNKVKERGTLLVGLEGTYPPFSFQGDDGKLTGFEVEFAQQLAKHLGVEASLKPTKWDGMLASLDSKRIDVVINQVTISDERKKKYDFSTPYTISGIQALVKKGNEGTIKTADDLKGKKVGVGLGTNYEEWLRQNVQGVDVRTYDDDPTKYQDLRVGRIDAILVDRLAALDLVKKTNDTLAVTGEAFSRQESGVALRKGNEDLLKAVNDAIAEMQKDGTLQALSEKWFGADVTK",
        )
        self.assertEqual(
            alignment.sequences[5].seq,
            "MKSVLKVSLAALTLAFAVSSHAADKKLVVATDTAFVPFEFKQGDKYVGFDVDLWAAIAKELKLDYELKPMDFSGIIPALQTKNVDLALAGITITDERKKAIDFSDGYYKSGLLVMVKANNNDVKSVKDLDGKVVAVKSGTGSVDYAKANIKTKDLRQFPNIDNAYMELGTNRADAVLHDTPNILYFIKTAGNGQFKAVGDSLEAQQYGIAFPKGSDELRDKVNGALKTLRENGTYNEIYKKWFGTEPK",
        )
        self.assertEqual(
            alignment.sequences[6].seq,
            "MKKSLLSLKLSGLLVPSVLALSLSACSSPSSTLNQGTLKIAMEGTYPPFTSKNEQGELVGFDVDIAKAVAQKLNLKPEFVLTEWSGILAGLQANKYDVIVNQVGITPERQNSIGFSQPYAYSRPEIIVAKNNTFNPQSLADLKGKRVGSTLGSNYEKQLIDTGDIKIVTYPGAPEILADLVAGRIDAAYNDRLVVNYIINDQKLPVRGAGQIGDAAPVGIALKKGNSALKDQIDKALTEMRSDGTFEKISQKWFGQDVGQP",
        )
        self.assertEqual(
            alignment.sequences[7].seq,
            "MKKLVLSLSLVLAFSSATAAFAAIPQNIRIGTDPTYAPFESKNSQGELVGFDIDLAKELCKRINTQCTFVENPLDALIPSLKAKKIDAIMSSLSITEKRQQEIAFTDKLYAADSRLVVAKNSDIQPTVESLKGKRVGVLQGTTQETFGNEHWAPKGIEIVSYQGQDNIYSDLTAGRIDAAFQDEVAASEGFLKQPVGKDYKFGGPSVKDEKLFGVGTGMGLRKEDNELREALNKAFAEMRADGTYEKLAKKYFDFDVYGG",
        )
        self.assertEqual(
            alignment[0],
            "--MKNWIKVAVAAIA--LSAA------------------TVQAATEVKVGMSGRYFPFTFVKQ--DKLQGFEVDMWDEIGKRNDYKIEYVTANFSGLFGLLETGRIDTISNQITMTDARKAKYLFADPYVVDG-AQITVRKGNDSIQGVEDLAGKTVAVNLGSNFEQLLRDYDKDGKINIKTYDT--GIEHDVALGRADAFIMDRLSALE-LIKKT-GLPLQLAGEPFETI-----QNAWPFVDNEKGRKLQAEVNKALAEMRADGTVEKISVKWFGADITK----",
        )
        self.assertEqual(
            alignment[1],
            "MKMKKWTVLVVAALLAVLSACG------------NGNSSSKEDDNVLHVGATGQSYPFAYKEN--GKLTGFDVEVMEAVAKKIDMKLDWKLLEFSGLMGELQTGKLDTISNQVAVTDERKETYNFTKPYAYAG-TQIVVKKDNTDIKSVDDLKGKTVAAVLGSNHAKNLESKDPDKKINIKTYETQEGTLKDVAYGRVDAYVNSRTVLIA-QIKKT-GLPLKLAGDPIVYE-----QVAFPFAKDDAHDKLRKKVNKALDELRKDGTLKKLSEKYFNEDITVEQKH",
        )
        self.assertEqual(
            alignment[2],
            "MKKALLALFMVVSIAALAACGAGNDNQSKDNAKDGDLWASIKKKGVLTVGTEGTYEPFTYHDKDTDKLTGYDVEVITEVAKRLGLKVDFKETQWGSMFAGLNSKRFDVVANQVG-KTDREDKYDFSDKYTTSR-AVVVTKKDNNDIKSEADVKGKTSAQSLTSNYNKLATN----AGAKVEGVEGMAQALQMIQQARVDMTYNDKLAVLN-YLKTSGNKNVKIAFETGEPQ-----STYFTFRKGS--GEVVDQVNKALKEMKEDGTLSKISKKWFGEDVSK----",
        )
        self.assertEqual(
            alignment[3],
            "MKKLLFTTALLTGAIAFSTF-----------SHAGEIADRVEKTKTLLVGTEGTYAPFTFHDK-SGKLTGFDVEVIRKVAEKLGLKVEFKETQWDAMYAGLNAKRFDVIANQTNPSPERLKKYSFTTPYNYSG-GVIVTKSSDNSIKSFEDLKGRKSAQSATSNWGKDAKA----AGAQILVVDGLAQSLELIKQGRAEATINDKLAVLD-YFKQHPNSGLKIAYDRGDKT-----PTAFAFLQGE--DALITKFNQVLEALRQDGTLKQISIEWFGYDITQ----",
        )
        self.assertEqual(
            alignment[4],
            "MKLAHLGRQALMGVMAVALVAG---MSVKSFADEG-LLNKVKERGTLLVGLEGTYPPFSFQGD-DGKLTGFEVEFAQQLAKHLGVEASLKPTKWDGMLASLDSKRIDVVINQVTISDERKKKYDFSTPYTISGIQALVKKGNEGTIKTADDLKGKKVGVGLGTNYEEWLRQNV--QGVDVRTYDDDPTKYQDLRVGRIDAILVDRLAALD-LVKKT-NDTLAVTGEAFSRQ-----ESGVALRKGN--EDLLKAVNDAIAEMQKDGTLQALSEKWFGADVTK----",
        )
        self.assertEqual(
            alignment[5],
            "--MKSVLKVSLAALTLAFAVS------------------SHAADKKLVVATDTAFVPFEFKQG--DKYVGFDVDLWAAIAKELKLDYELKPMDFSGIIPALQTKNVDLALAGITITDERKKAIDFSDGYYKSG-LLVMVKANNNDVKSVKDLDGKVVAVKSGTGSVDYAKAN--IKTKDLRQFPNIDNAYMELGTNRADAVLHDTPNILY-FIKTAGNGQFKAVGDSLEAQ-----QYGIAFPKGS--DELRDKVNGALKTLRENGTYNEIYKKWFGTEPK-----",
        )
        self.assertEqual(
            alignment[6],
            "-MKKSLLSLKLSGLLVPSVLALS--------LSACSSPSSTLNQGTLKIAMEGTYPPFTSKNE-QGELVGFDVDIAKAVAQKLNLKPEFVLTEWSGILAGLQANKYDVIVNQVGITPERQNSIGFSQPYAYSRPEIIVAKNNTFNPQSLADLKGKRVGSTLGSNYEKQLIDTG---DIKIVTYPGAPEILADLVAGRIDAAYNDRLVVNY-IINDQ-KLPVRGAGQIGDAA-----PVGIALKKGN--SALKDQIDKALTEMRSDGTFEKISQKWFGQDVGQP---",
        )
        self.assertEqual(
            alignment[7],
            "MKKLVLSLSLVLAFSSATAAF-------------------AAIPQNIRIGTDPTYAPFESKNS-QGELVGFDIDLAKELCKRINTQCTFVENPLDALIPSLKAKKIDAIMSSLSITEKRQQEIAFTDKLYAADSRLVVAKNSDIQP-TVESLKGKRVGVLQGTTQETFGNEHWAPKGIEIVSYQGQDNIYSDLTAGRIDAAFQDEVAASEGFLKQPVGKDYKFGGPSVKDEKLFGVGTGMGLRKED--NELREALNKAFAEMRADGTYEKLAKKYFDFDVYGG---",
        )
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [
                        [
                            0,
                            0,
                            0,
                            13,
                            13,
                            16,
                            17,
                            17,
                            17,
                            17,
                            17,
                            17,
                            17,
                            17,
                            17,
                            18,
                            41,
                            41,
                            41,
                            90,
                            91,
                            109,
                            109,
                            121,
                            122,
                            146,
                            147,
                            148,
                            149,
                            150,
                            151,
                            160,
                            160,
                            183,
                            183,
                            188,
                            188,
                            202,
                            202,
                            212,
                            214,
                            247,
                            248,
                            248,
                            248,
                        ],
                        [
                            0,
                            1,
                            2,
                            15,
                            17,
                            20,
                            21,
                            22,
                            22,
                            22,
                            22,
                            22,
                            23,
                            24,
                            27,
                            28,
                            51,
                            51,
                            51,
                            100,
                            101,
                            119,
                            119,
                            131,
                            132,
                            156,
                            157,
                            158,
                            159,
                            160,
                            161,
                            170,
                            172,
                            195,
                            195,
                            200,
                            200,
                            214,
                            214,
                            224,
                            226,
                            259,
                            260,
                            261,
                            264,
                        ],
                        [
                            0,
                            1,
                            2,
                            15,
                            17,
                            20,
                            21,
                            22,
                            23,
                            25,
                            31,
                            34,
                            35,
                            36,
                            39,
                            40,
                            63,
                            64,
                            65,
                            114,
                            114,
                            132,
                            132,
                            144,
                            145,
                            169,
                            169,
                            169,
                            169,
                            169,
                            170,
                            179,
                            181,
                            204,
                            204,
                            209,
                            210,
                            224,
                            224,
                            234,
                            234,
                            267,
                            268,
                            268,
                            268,
                        ],
                        [
                            0,
                            1,
                            2,
                            15,
                            17,
                            20,
                            20,
                            20,
                            20,
                            20,
                            20,
                            23,
                            24,
                            25,
                            28,
                            29,
                            52,
                            52,
                            53,
                            102,
                            103,
                            121,
                            121,
                            133,
                            134,
                            158,
                            158,
                            158,
                            158,
                            158,
                            159,
                            168,
                            170,
                            193,
                            193,
                            198,
                            199,
                            213,
                            213,
                            223,
                            223,
                            256,
                            257,
                            257,
                            257,
                        ],
                        [
                            0,
                            1,
                            2,
                            15,
                            17,
                            20,
                            21,
                            22,
                            22,
                            22,
                            28,
                            31,
                            32,
                            32,
                            35,
                            36,
                            59,
                            59,
                            60,
                            109,
                            110,
                            128,
                            129,
                            141,
                            142,
                            166,
                            167,
                            168,
                            168,
                            168,
                            169,
                            178,
                            180,
                            203,
                            203,
                            208,
                            208,
                            222,
                            222,
                            232,
                            232,
                            265,
                            266,
                            266,
                            266,
                        ],
                        [
                            0,
                            0,
                            0,
                            13,
                            15,
                            18,
                            19,
                            19,
                            19,
                            19,
                            19,
                            19,
                            19,
                            19,
                            19,
                            20,
                            43,
                            43,
                            43,
                            92,
                            93,
                            111,
                            111,
                            123,
                            124,
                            148,
                            149,
                            149,
                            149,
                            150,
                            151,
                            160,
                            162,
                            185,
                            185,
                            190,
                            191,
                            205,
                            205,
                            215,
                            215,
                            248,
                            248,
                            248,
                            248,
                        ],
                        [
                            0,
                            0,
                            1,
                            14,
                            16,
                            19,
                            20,
                            21,
                            22,
                            22,
                            22,
                            25,
                            26,
                            27,
                            30,
                            31,
                            54,
                            54,
                            55,
                            104,
                            105,
                            123,
                            124,
                            136,
                            137,
                            161,
                            162,
                            163,
                            163,
                            163,
                            163,
                            172,
                            174,
                            197,
                            197,
                            202,
                            202,
                            216,
                            216,
                            226,
                            226,
                            259,
                            260,
                            261,
                            261,
                        ],
                        [
                            0,
                            1,
                            2,
                            15,
                            17,
                            20,
                            21,
                            21,
                            21,
                            21,
                            21,
                            21,
                            21,
                            21,
                            21,
                            21,
                            44,
                            44,
                            45,
                            94,
                            95,
                            113,
                            114,
                            126,
                            126,
                            150,
                            151,
                            152,
                            153,
                            154,
                            155,
                            164,
                            166,
                            189,
                            190,
                            195,
                            196,
                            210,
                            215,
                            225,
                            225,
                            258,
                            259,
                            260,
                            260,
                        ],
                    ]
                ),
            )
        )
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
                "<Alignment object (5 rows x 60 columns) at 0x%x>"
                % id(alignment),
            )
            self.assertEqual(len(alignment), 5)
            self.assertEqual(alignment.sequences[0].id, "Tax1")
            self.assertEqual(alignment.sequences[1].id, "Tax2")
            self.assertEqual(alignment.sequences[2].id, "Tax3")
            self.assertEqual(alignment.sequences[3].id, "Tax4")
            self.assertEqual(alignment.sequences[4].id, "Tax5")
            self.assertEqual(
                alignment.sequences[0].seq,
                "CCATCTCACGGTCGGTACGATACACCTGCTTTTGGCAGGAAATGGTCAATATTACAAGGT",
            )
            self.assertEqual(
                alignment.sequences[1].seq,
                "CCATCTCACGGTCAGTAAGATACACCTGCTTTTGGCGGGAAATGGTCAACATTAAAAGAT",
            )
            self.assertEqual(
                alignment.sequences[2].seq,
                "CCATCTCCCGCTCAGTAAGATACCCCTGCTGTTGGCGGGAAATCGTCAATATTAAAAGGT",
            )
            self.assertEqual(
                alignment.sequences[3].seq,
                "TCATCTCATGGTCAATAAGATACTCCTGCTTTTGGCGGGAAATGGTCAATCTTAAAAGGT",
            )
            self.assertEqual(
                alignment.sequences[4].seq,
                "CCATCTCACGGTCGGTAAGATACACCTGCTTTTGGCGGGAAATGGTCAATATTAAAAGGT",
            )
            self.assertEqual(
                alignment[0],
                "CCATCTCACGGTCGGTACGATACACCTGCTTTTGGCAGGAAATGGTCAATATTACAAGGT",
            )
            self.assertEqual(
                alignment[1],
                "CCATCTCACGGTCAGTAAGATACACCTGCTTTTGGCGGGAAATGGTCAACATTAAAAGAT",
            )
            self.assertEqual(
                alignment[2],
                "CCATCTCCCGCTCAGTAAGATACCCCTGCTGTTGGCGGGAAATCGTCAATATTAAAAGGT",
            )
            self.assertEqual(
                alignment[3],
                "TCATCTCATGGTCAATAAGATACTCCTGCTTTTGGCGGGAAATGGTCAATCTTAAAAGGT",
            )
            self.assertEqual(
                alignment[4],
                "CCATCTCACGGTCGGTAAGATACACCTGCTTTTGGCGGGAAATGGTCAATATTAAAAGGT",
            )
            self.check_reading_writing(path)
            self.assertTrue(
                numpy.array_equal(
                    alignment.coordinates,
                    numpy.array([[0, 60], [0, 60], [0, 60], [0, 60], [0, 60]]),
                )
            )

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
            "<Alignment object (5 rows x 42 columns) at 0x%x>"
            % id(alignment),
        )
        self.assertEqual(len(alignment), 5)
        self.assertEqual(alignment.sequences[0].id, "Turkey")
        self.assertEqual(alignment.sequences[1].id, "Salmo gair")
        self.assertEqual(alignment.sequences[2].id, "H. Sapiens")
        self.assertEqual(alignment.sequences[3].id, "Chimp")
        self.assertEqual(alignment.sequences[4].id, "Gorilla")
        self.assertEqual(
            alignment.sequences[0].seq, "AAGCTNGGGCATTTCAGGGTGAGCCCGGGCAATACAGGGTAT"
        )
        self.assertEqual(
            alignment.sequences[1].seq, "AAGCCTTGGCAGTGCAGGGTGAGCCGTGGCCGGGCACGGTAT"
        )
        self.assertEqual(
            alignment.sequences[2].seq, "ACCGGTTGGCCGTTCAGGGTACAGGTTGGCCGTTCAGGGTAA"
        )
        self.assertEqual(
            alignment.sequences[3].seq, "AAACCCTTGCCGTTACGCTTAAACCGAGGCCGGGACACTCAT"
        )
        self.assertEqual(
            alignment.sequences[4].seq, "AAACCCTTGCCGGTACGCTTAAACCATTGCCGGTACGCTTAA"
        )
        self.assertEqual(alignment[0], "AAGCTNGGGCATTTCAGGGTGAGCCCGGGCAATACAGGGTAT")
        self.assertEqual(alignment[1], "AAGCCTTGGCAGTGCAGGGTGAGCCGTGGCCGGGCACGGTAT")
        self.assertEqual(alignment[2], "ACCGGTTGGCCGTTCAGGGTACAGGTTGGCCGTTCAGGGTAA")
        self.assertEqual(alignment[3], "AAACCCTTGCCGTTACGCTTAAACCGAGGCCGGGACACTCAT")
        self.assertEqual(alignment[4], "AAACCCTTGCCGGTACGCTTAAACCATTGCCGGTACGCTTAA")
        self.check_reading_writing(path)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[0, 42], [0, 42], [0, 42], [0, 42], [0, 42]]),
            )
        )

    def test_five_and_six(self):
        paths = ("Phylip/five.dat", "Phylip/six.dat")
        # http://evolution.genetics.washington.edu/phylip/doc/sequence.html
        for path in paths:
            with open(path) as stream:
                alignments = AlignmentIterator(stream)
                alignment = next(alignments)
                with self.assertRaises(StopIteration):
                    next(alignments)
            self.assertEqual(
                repr(alignment),
                "<Alignment object (5 rows x 42 columns) at 0x%x>"
                % id(alignment),
            )
            self.assertEqual(len(alignment), 5)
            self.assertEqual(alignment.sequences[0].id, "Turkey")
            self.assertEqual(alignment.sequences[1].id, "Salmo gair")
            self.assertEqual(alignment.sequences[2].id, "H. Sapiens")
            self.assertEqual(alignment.sequences[3].id, "Chimp")
            self.assertEqual(alignment.sequences[4].id, "Gorilla")
            self.assertEqual(
                alignment.sequences[0].seq, "AAGCTNGGGCATTTCAGGGTGAGCCCGGGCAATACAGGGTAT"
            )
            self.assertEqual(
                alignment.sequences[1].seq, "AAGCCTTGGCAGTGCAGGGTGAGCCGTGGCCGGGCACGGTAT"
            )
            self.assertEqual(
                alignment.sequences[2].seq, "ACCGGTTGGCCGTTCAGGGTACAGGTTGGCCGTTCAGGGTAA"
            )
            self.assertEqual(
                alignment.sequences[3].seq, "AAACCCTTGCCGTTACGCTTAAACCGAGGCCGGGACACTCAT"
            )
            self.assertEqual(
                alignment.sequences[4].seq, "AAACCCTTGCCGGTACGCTTAAACCATTGCCGGTACGCTTAA"
            )
            self.assertEqual(alignment[0], "AAGCTNGGGCATTTCAGGGTGAGCCCGGGCAATACAGGGTAT")
            self.assertEqual(alignment[1], "AAGCCTTGGCAGTGCAGGGTGAGCCGTGGCCGGGCACGGTAT")
            self.assertEqual(alignment[2], "ACCGGTTGGCCGTTCAGGGTACAGGTTGGCCGTTCAGGGTAA")
            self.assertEqual(alignment[3], "AAACCCTTGCCGTTACGCTTAAACCGAGGCCGGGACACTCAT")
            self.assertEqual(alignment[4], "AAACCCTTGCCGGTACGCTTAAACCATTGCCGGTACGCTTAA")
            self.assertTrue(
                numpy.array_equal(
                    alignment.coordinates,
                    numpy.array([[0, 42], [0, 42], [0, 42], [0, 42], [0, 42]]),
                )
            )
            self.check_reading_writing(path)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
